"""
rerun_failed_gemini.py — Targeted rerun for failed entries in Gemini annotator folders.

Scans both Gemini annotator folders for entries where the annotation failed,
re-annotates only those entries using the correct Gemini model, and writes
corrected files back in-place.

Usage:
    # Rerun Gemini 2.0 Flash failures:
    python rerun_failed_gemini.py --annotator flash

    # Rerun Gemini 3.1 Pro failures:
    python rerun_failed_gemini.py --annotator pro

    # Dry-run (just count failures, no API calls):
    python rerun_failed_gemini.py --annotator flash --dry-run
    python rerun_failed_gemini.py --annotator pro   --dry-run
"""

import argparse
import json
import os
import re
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm
from dotenv import load_dotenv
from google import genai
from google.genai import types as genai_types

from annotate_utils import SYSTEM_PROMPT, build_user_prompt, validate_and_fix, INPUT_DIR

load_dotenv()

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE     = Path(__file__).parent.parent          # Error Analysis Annotations/
SCRIPTS  = Path(__file__).parent                 # Error Analysis Annotations/scripts/

ANNOTATORS = {
    "flash": {
        "folder":       BASE / "Annotator C",
        "annotator_id": "Annotator C",
        "model_name":   "gemini-2.0-flash",
    },
    "pro": {
        "folder":       BASE / "Annotator B",
        "annotator_id": "Annotator B",
        "model_name":   "gemini-2.5-pro",
    },
    "pro25": {
        "folder":       SCRIPTS / "Annotator Gemini 2.5 Pro",
        "annotator_id": "Annotator Gemini",
        "model_name":   "gemini-2.5-pro",
    },
}

N_WORKERS   = 10
MAX_RETRIES = 3
RETRY_DELAY = 5.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_failed(trace: dict) -> bool:
    span  = str(trace.get("triggering_text_span", ""))
    notes = str(trace.get("annotator_notes", ""))
    return "[annotation failed" in span or "API error" in notes


_BS = chr(92)  # literal backslash
_BAD_ESCAPE_RE = re.compile(r'(\\\\)|\\([^"\\\/bfnrtu])')


def _fix_bad_escapes(s: str) -> str:
    """Escape bare backslashes that form invalid JSON escape sequences.

    Gemini sometimes outputs LaTeX like \\, or \\vec inside JSON string values.
    These are invalid JSON escapes; we double the lone backslash so json.loads works.
    Already-valid \\\\ pairs are left untouched.
    """
    def _replace(m: re.Match) -> str:
        if m.group(1) is not None:          # matched valid \\\\ pair — keep as-is
            return m.group(1)
        return _BS + _BS + m.group(2)       # bare \\X → \\\\X

    return _BAD_ESCAPE_RE.sub(_replace, s)


def extract_json(text: str) -> dict:
    s = text.strip()
    # Strip markdown fences if present
    s = re.sub(r"^```(?:json)?\s*", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s*```$", "", s)
    # Try direct parse first
    try:
        return json.loads(s)
    except json.JSONDecodeError:
        pass
    # Fix invalid backslash escapes (e.g. LaTeX \, \vec \hat inside JSON strings)
    try:
        return json.loads(_fix_bad_escapes(s))
    except json.JSONDecodeError:
        pass
    # Last resort: extract first {...} block and fix escapes
    m = re.search(r"\{[\s\S]+\}", s)
    if m:
        block = m.group(0)
        try:
            return json.loads(_fix_bad_escapes(block))
        except json.JSONDecodeError:
            return json.loads(block)
    raise ValueError(f"Could not parse JSON from response:\n{text[:300]}")


# ---------------------------------------------------------------------------
# Gemini caller
# ---------------------------------------------------------------------------

def make_gemini_caller(model_name: str):
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        raise EnvironmentError("GEMINI_API_KEY not found in .env")
    client = genai.Client(api_key=api_key)

    def call(trace: dict, annotator_id: str) -> dict:
        user_msg = build_user_prompt(trace, annotator_id)
        response = client.models.generate_content(
            model=model_name,
            contents=user_msg,
            config=genai_types.GenerateContentConfig(
                system_instruction=SYSTEM_PROMPT,
                temperature=0.1,
                response_mime_type="application/json",
            ),
        )
        return extract_json(response.text)

    return call


# ---------------------------------------------------------------------------
# Rerun one file in-place
# ---------------------------------------------------------------------------

def rerun_file(
    filepath: Path,
    call_fn,
    annotator_id: str,
    dry_run: bool = False,
) -> tuple:
    """Returns (n_failed_found, n_fixed)."""

    # Load current output file
    traces = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                traces.append(json.loads(line))
            except json.JSONDecodeError:
                pass

    failed_indices = [i for i, t in enumerate(traces) if is_failed(t)]
    n_failed = len(failed_indices)

    if n_failed == 0:
        return 0, 0

    if dry_run:
        print(f"    [DRY-RUN] {filepath.name}: {n_failed} failed entries.")
        return n_failed, 0

    # Load original input traces for full problem/gold/reasoning fields
    input_file = INPUT_DIR / filepath.name
    orig_by_id = {}
    if input_file.exists():
        with open(input_file, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    t = json.loads(line)
                    orig_by_id[t.get("annotation_id")] = t
                except json.JSONDecodeError:
                    pass

    write_lock = threading.Lock()
    fixed_map  = {}

    def fix_one(idx_trace: tuple):
        idx, trace = idx_trace
        aid    = trace.get("annotation_id")
        source = orig_by_id.get(aid, trace)
        delay  = RETRY_DELAY
        last_exc = None
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                raw = call_fn(source, annotator_id)
                ann = validate_and_fix(raw, annotator_id)
                return idx, ann, None
            except Exception as exc:
                last_exc = exc
                if attempt < MAX_RETRIES:
                    time.sleep(delay)
                    delay *= 2
        return idx, None, last_exc

    pending = [(i, traces[i]) for i in failed_indices]

    with tqdm(
        total=n_failed,
        desc=f"  {filepath.name[:50]}",
        unit="trace",
        dynamic_ncols=True,
        leave=True,
    ) as pbar:
        with ThreadPoolExecutor(max_workers=N_WORKERS) as executor:
            futures = {executor.submit(fix_one, item): item for item in pending}
            for future in as_completed(futures):
                idx, ann, exc = future.result()
                if ann is not None:
                    with write_lock:
                        fixed_map[idx] = ann
                pbar.update(1)

    # Merge fixes into traces
    for idx, ann in fixed_map.items():
        traces[idx].update(ann)

    # Write back in-place via a temp file (safe against partial writes)
    tmp = filepath.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        for t in traces:
            fh.write(json.dumps(t, ensure_ascii=False) + "\n")
    tmp.replace(filepath)

    n_fixed       = len(fixed_map)
    still_failing = n_failed - n_fixed
    if still_failing:
        print(f"  [WARN] {filepath.name}: {still_failing} entries still failed after retries.")
    return n_failed, n_fixed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Targeted rerun of failed Gemini annotation entries."
    )
    parser.add_argument(
        "--annotator",
        required=True,
        choices=["flash", "pro", "pro25"],
        help="flash = Gemini 2.0 Flash | pro = Gemini 2.5 Pro (Annotator B) | pro25 = Gemini 2.5 Pro (scripts/Annotator Gemini 2.5 Pro)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Count failures only; do not call the API.",
    )
    args = parser.parse_args()

    cfg          = ANNOTATORS[args.annotator]
    folder       = cfg["folder"]
    annotator_id = cfg["annotator_id"]
    model_name   = cfg["model_name"]

    if not folder.exists():
        print(f"[ERROR] Folder not found: {folder}")
        sys.exit(1)

    call_fn = None if args.dry_run else make_gemini_caller(model_name)

    files = sorted(folder.glob("*.jsonl"))
    if not files:
        print(f"[ERROR] No .jsonl files in {folder}")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  Rerun  : {annotator_id}")
    print(f"  Model  : {model_name}")
    print(f"  Mode   : {'DRY RUN' if args.dry_run else 'LIVE'}")
    print(f"  Files  : {len(files)}")
    print(f"{'='*60}")

    total_failed = 0
    total_fixed  = 0

    for fpath in files:
        n_fail, n_fix = rerun_file(
            filepath=fpath,
            call_fn=call_fn,
            annotator_id=annotator_id,
            dry_run=args.dry_run,
        )
        if n_fail > 0 and not args.dry_run:
            print(f"  {fpath.name}: {n_fail} failed → {n_fix} fixed")
        total_failed += n_fail
        total_fixed  += n_fix

    print(f"\n{'='*60}")
    if args.dry_run:
        print(f"  Total failed entries: {total_failed}")
    else:
        print(f"  Total failed : {total_failed}")
        print(f"  Total fixed  : {total_fixed}")
        print(f"  Still failing: {total_failed - total_fixed}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
