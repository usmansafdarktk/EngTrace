"""
annotate_utils.py — Shared utilities for EngTrace LLM annotation scripts.

Provides:
  - SYSTEM_PROMPT  : Full annotation manual baked into a prompt string
  - USER_PROMPT    : Per-trace prompt builder
  - OUTPUT_SCHEMA  : Expected JSON fields
  - process_file() : File-safe, concurrent, resumable annotation runner
"""

import json
import os
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

from tqdm import tqdm
from dotenv import load_dotenv

load_dotenv()

# ---------------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------------

INPUT_DIR = Path(__file__).parent.parent / "samples"

MODEL_FILES = [
    "gemini-3.1-pro-preview_annotation_sample.jsonl",
    "deepseek-r1_annotation_sample.jsonl",
    "meta-llama-3.1-70b_annotation_sample.jsonl",
    "gpt-5-mini_annotation_sample.jsonl",
    "gpt-5_annotation_sample.jsonl",
    "deepseek-v4-pro_annotation_sample.jsonl",
    "claude-opus-4-7_annotation_sample.jsonl",
    "qwen2.5-72b_annotation_sample.jsonl",
    "gemma-3-27b_annotation_sample.jsonl",
    "mathstral-7b_annotation_sample.jsonl",
    "wizardmath-7b_annotation_sample.jsonl",
]

# ---------------------------------------------------------------------------
# Full annotation system prompt — every rule from the EngTrace manual
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """
You are an expert engineering annotator working on the EngTrace benchmark.
Your task is to identify the root-cause error in an AI model's reasoning for
an engineering problem that the model answered INCORRECTLY.

══════════════════════════════════════════════════════════════════
ANNOTATION MANUAL — READ CAREFULLY BEFORE ANNOTATING
══════════════════════════════════════════════════════════════════

## STEP 1 — What you receive
You will be given:
  • problem_statement  : The full engineering problem text.
  • gold_answer        : The correct reference solution. Treat this as ground truth.
  • model_reasoning    : The model's chain-of-thought. THIS is what you annotate.

The model's final answer is confirmed WRONG (final_answer_acc = 0.0).

## STEP 2 — Six error categories (ordered by severity)

| # | Category                  | Definition                                                                                      | Severity |
|---|---------------------------|-------------------------------------------------------------------------------------------------|----------|
| 1 | Hallucination             | Model uses a value, constant, or equation that does NOT exist in any real engineering reference — it was fabricated. | Critical |
| 2 | Setup / Assumption Error  | Problem is framed incorrectly BEFORE any equation is applied — wrong boundary condition, ignored stated constraint, or misidentified system configuration. | High |
| 3 | Formula / Principle Error | Setup is correct but the WRONG real governing equation, law, or principle is selected.          | High |
| 4 | Unit / Dimensional Error  | Correct formula and setup, but dimensional consistency fails — missing conversion, mixed unit systems, or value substituted in wrong unit. | Medium |
| 5 | Sign / Direction Error    | Correct formula, units, and setup — but WRONG sign or physical direction assigned to a quantity, violating the reference frame or process direction. | Medium |
| 6 | Calculation Error         | Everything is correct — setup, equation, units, sign — but the arithmetic or algebraic execution produces the wrong number. | Low |

## STEP 3 — Decision tree (work through in order, stop at first YES)

Q1. Does the model use a value, constant, or equation that does NOT exist in
    any real engineering reference (i.e., it was fabricated)?
    → YES → Category 1: Hallucination

Q2. Is the problem framed incorrectly BEFORE any equation is applied?
    (wrong boundary condition, ignored stated constraint, misidentified
    system type, wrong process assumption not supported by the problem)
    → YES → Category 2: Setup / Assumption Error

Q3. Is the selected governing equation or principle WRONG, given that the
    setup is correct?
    (e.g., using ideal gas law when van der Waals is required; using CSTR
    design equation for a PFR; applying Bernoulli to a viscous-flow problem;
    using Carnot efficiency for a non-Carnot cycle)
    → YES → Category 3: Formula / Principle Error

Q4. Does dimensional analysis fail?
    (missing unit conversion, mixing SI and imperial, value substituted in
    wrong unit — e.g., kPa used where Pa is required, °C used directly in
    an absolute-temperature formula, RPM not converted to rad/s, cm not
    converted to m in a formula)
    → YES → Category 4: Unit / Dimensional Error

Q5. Is the SIGN or physical direction of any key quantity wrong, given that
    formula, units, and setup are all correct?
    (e.g., heat addition given a negative sign, work done BY system given
    positive sign, wrong sign on enthalpy of exothermic reaction, wrong
    current direction, wrong moment direction)
    → YES → Category 5: Sign / Direction Error

Q6. Does the arithmetic or algebraic execution produce the wrong number,
    while setup, formula, units, and signs are all correct?
    → YES → Category 6: Calculation Error

IMPORTANT PRIORITY RULE: A higher-level error invalidates everything below
it. Always assign the FIRST matching category. Never label something a
Calculation Error if a more fundamental error exists higher up.

## STEP 4 — triggering_text_span

Copy ONE short verbatim excerpt from model_reasoning — one sentence or one
equation — that directly shows the error. This is the evidence for your label.
Keep it under 300 characters. Do NOT paraphrase; copy exactly from the text.

## STEP 5 — annotator_notes (optional)

Use this field to flag anything ambiguous, borderline, or worth explaining
about your labelling decision. Leave empty string "" if nothing to add.

══════════════════════════════════════════════════════════════════
OUTPUT FORMAT
══════════════════════════════════════════════════════════════════

Return ONLY a valid JSON object — no markdown, no explanation outside JSON.
The object must contain exactly these six fields:

{
  "decision_step_triggered": "Q1" | "Q2" | "Q3" | "Q4" | "Q5" | "Q6",
  "error_category": 1 | 2 | 3 | 4 | 5 | 6,
  "error_category_name": "Hallucination" | "Setup / Assumption Error" | "Formula / Principle Error" | "Unit / Dimensional Error" | "Sign / Direction Error" | "Calculation Error",
  "triggering_text_span": "<verbatim excerpt from model_reasoning, ≤300 chars>",
  "annotator_id": "<your annotator name>",
  "annotator_notes": "<optional explanation, or empty string>"
}

The three category fields (decision_step_triggered, error_category,
error_category_name) must be CONSISTENT with each other.
""".strip()


def build_user_prompt(trace: dict, annotator_id: str) -> str:
    """Build the per-trace user message."""
    return f"""Please annotate the following engineering reasoning trace.

--- PROBLEM STATEMENT ---
{trace.get('problem_statement', '').strip()}

--- GOLD ANSWER (correct solution) ---
{trace.get('gold_answer', '').strip()}

--- MODEL REASONING (what you are annotating) ---
{trace.get('model_reasoning', '').strip()}

---
Metadata: difficulty={trace.get('difficulty','')}, branch={trace.get('branch','')}, domain={trace.get('domain','')}
Your annotator_id should be: "{annotator_id}"

Apply the Q1→Q6 decision tree and return ONLY the JSON object.
""".strip()


# ---------------------------------------------------------------------------
# Expected output fields (for validation / fallback)
# ---------------------------------------------------------------------------

VALID_CATEGORIES = {
    1: ("Q1", "Hallucination"),
    2: ("Q2", "Setup / Assumption Error"),
    3: ("Q3", "Formula / Principle Error"),
    4: ("Q4", "Unit / Dimensional Error"),
    5: ("Q5", "Sign / Direction Error"),
    6: ("Q6", "Calculation Error"),
}


def validate_and_fix(ann: dict, annotator_id: str) -> dict:
    """Validate annotation JSON; coerce minor issues; return cleaned dict."""
    # Ensure error_category is int
    cat = ann.get("error_category")
    try:
        cat = int(cat)
    except (TypeError, ValueError):
        cat = 6
    cat = max(1, min(6, cat))

    q_step, cat_name = VALID_CATEGORIES[cat]

    return {
        "decision_step_triggered": ann.get("decision_step_triggered", q_step) or q_step,
        "error_category": cat,
        "error_category_name": ann.get("error_category_name", cat_name) or cat_name,
        "triggering_text_span": str(ann.get("triggering_text_span", ""))[:350],
        "annotator_id": annotator_id,
        "annotator_notes": str(ann.get("annotator_notes", "") or ""),
        "adjudicated": "",
        "final_label": "",
    }


# ---------------------------------------------------------------------------
# Resume helper — returns set of already-done annotation_ids in output file
# ---------------------------------------------------------------------------

def load_done_ids(output_path: Path) -> set:
    done = set()
    if not output_path.exists():
        return done
    with open(output_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                aid = obj.get("annotation_id")
                # Only count as done if annotation fields are actually filled
                if aid and obj.get("error_category") and obj.get("triggering_text_span"):
                    done.add(aid)
            except json.JSONDecodeError:
                pass
    return done


# ---------------------------------------------------------------------------
# Core file processor
# ---------------------------------------------------------------------------

def process_file(
    input_path: Path,
    output_path: Path,
    call_model_fn: Callable[[dict], dict],
    annotator_id: str,
    n_workers: int = 15,
    max_retries: int = 3,
    retry_delay: float = 5.0,
):
    """
    Annotate one JSONL file with concurrent workers.

    Parameters
    ----------
    input_path     : path to unannotated input file
    output_path    : path to write annotated output (appended to, not overwritten)
    call_model_fn  : function(trace_dict) → annotation_dict; raises on failure
    annotator_id   : string to put in annotator_id field
    n_workers      : number of concurrent threads
    max_retries    : API retry attempts per trace
    retry_delay    : base delay (seconds) between retries (doubles each attempt)
    """
    # Load all traces
    traces = []
    with open(input_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                traces.append(json.loads(line))
            except json.JSONDecodeError:
                pass

    if not traces:
        print(f"  [WARN] No traces found in {input_path.name}")
        return

    # Resume: find already-completed annotation_ids
    output_path.parent.mkdir(parents=True, exist_ok=True)
    done_ids = load_done_ids(output_path)
    pending = [t for t in traces if t.get("annotation_id") not in done_ids]

    if not pending:
        print(f"  [SKIP] {input_path.name} — all {len(traces)} traces already annotated.")
        return

    print(f"\n  File : {input_path.name}")
    print(f"  Total: {len(traces)} | Done: {len(done_ids)} | Remaining: {len(pending)}")

    # Thread-safe file writer
    write_lock = threading.Lock()

    def write_result(trace: dict, annotation: dict):
        merged = {**trace, **annotation}
        with write_lock:
            with open(output_path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(merged, ensure_ascii=False) + "\n")

    # Worker function with retries
    def annotate_with_retry(trace: dict) -> tuple:
        """Returns (trace, annotation_dict) or raises after max_retries."""
        delay = retry_delay
        last_exc = None
        for attempt in range(1, max_retries + 1):
            try:
                raw = call_model_fn(trace)
                ann = validate_and_fix(raw, annotator_id)
                return trace, ann
            except Exception as exc:
                last_exc = exc
                if attempt < max_retries:
                    time.sleep(delay)
                    delay *= 2
        raise last_exc

    # Run with progress bar
    failed = []
    with tqdm(
        total=len(pending),
        desc=f"  {input_path.stem[:40]}",
        unit="trace",
        dynamic_ncols=True,
        leave=True,
    ) as pbar:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(annotate_with_retry, t): t for t in pending}
            for future in as_completed(futures):
                trace = futures[future]
                try:
                    _, ann = future.result()
                    write_result(trace, ann)
                except Exception as exc:
                    failed.append(trace.get("annotation_id", "?"))
                    # Write a fallback record so the file stays consistent
                    fallback = {
                        "decision_step_triggered": "Q6",
                        "error_category": 6,
                        "error_category_name": "Calculation Error",
                        "triggering_text_span": "[annotation failed — API error]",
                        "annotator_id": annotator_id,
                        "annotator_notes": f"API error after {max_retries} retries: {exc}",
                        "adjudicated": "",
                        "final_label": "",
                    }
                    write_result(trace, fallback)
                finally:
                    pbar.update(1)

    if failed:
        print(f"  [WARN] {len(failed)} traces failed after retries: {failed[:10]}")
    else:
        print(f"  [OK]   {input_path.name} — {len(pending)} traces annotated.")


# ---------------------------------------------------------------------------
# Top-level runner — iterates all model files sequentially
# ---------------------------------------------------------------------------

def run_all_files(
    output_dir: Path,
    call_model_fn: Callable[[dict], dict],
    annotator_id: str,
    n_workers: int = 15,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"  EngTrace Annotation — {annotator_id}")
    print(f"  Output : {output_dir}")
    print(f"  Workers: {n_workers}")
    print(f"{'='*60}")

    for fname in MODEL_FILES:
        inp = INPUT_DIR / fname
        out = output_dir / fname
        if not inp.exists():
            print(f"\n  [MISSING] {fname} — skipping.")
            continue
        process_file(
            input_path=inp,
            output_path=out,
            call_model_fn=call_model_fn,
            annotator_id=annotator_id,
            n_workers=n_workers,
        )

    print(f"\n{'='*60}")
    print(f"  All files processed for {annotator_id}.")
    print(f"{'='*60}\n")
