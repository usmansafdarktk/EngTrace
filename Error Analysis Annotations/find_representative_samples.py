"""
Find representative error samples where all 3 annotators agree,
covering diverse error categories, with clear reasoning traces.
"""

import json
import os
import numpy as np
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ANNOTATORS = ["Annotator A", "Annotator B", "Annotator C"]
MODELS = [
    "claude-opus-4-7", "deepseek-r1", "deepseek-v4-pro",
    "gemini-3.1-pro-preview", "gemma-3-27b", "gpt-5", "gpt-5-mini",
    "mathstral-7b", "meta-llama-3.1-70b", "qwen2.5-72b", "wizardmath-7b",
]
CATEGORY_NAMES = {
    0: "No Error", 1: "Hallucination", 2: "Setup/Assumption Error",
    3: "Formula/Principle Error", 4: "Unit/Dimensional Error",
    5: "Sign/Direction Error", 6: "Calculation Error",
}
NAME_TO_CODE = {v.lower(): k for k, v in CATEGORY_NAMES.items()}
NAME_TO_CODE.update({"no error": 0, "hallucination": 1, "setup": 2,
                     "formula": 3, "unit": 4, "sign": 5, "calculation": 6})


def resolve_cat(cat):
    if cat is None or cat == "" or (isinstance(cat, float) and np.isnan(cat)):
        return 0
    if isinstance(cat, str):
        key = cat.lower().strip()
        for name, code in NAME_TO_CODE.items():
            if name in key or key in name:
                return code
        return 0
    return int(cat)


def load_all():
    # aid -> {annotator: (record, category)}
    by_id = defaultdict(dict)
    for annotator in ANNOTATORS:
        for model in MODELS:
            fpath = os.path.join(BASE_DIR, annotator, f"{model}_annotation_sample.jsonl")
            if not os.path.exists(fpath):
                continue
            with open(fpath, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    aid = rec["annotation_id"]
                    cat = resolve_cat(rec.get("error_category"))
                    by_id[aid][annotator] = (rec, cat)
    return by_id


def find_candidates(by_id):
    # Group unanimous (all 3 agree) samples by category
    unanimous = defaultdict(list)
    for aid, ann_map in by_id.items():
        if len(ann_map) < 3:
            continue
        cats = [ann_map[a][1] for a in ANNOTATORS]
        if len(set(cats)) == 1 and cats[0] != 0:  # skip "no error"
            category = cats[0]
            rec = ann_map[ANNOTATORS[0]][0]
            notes = [ann_map[a][0].get("annotator_notes", "") for a in ANNOTATORS]
            # Quality score: longer model_reasoning + clear notes + non-trivial problem
            reasoning_len = len(rec.get("model_reasoning", ""))
            notes_len = sum(len(n) for n in notes if n)
            score = reasoning_len + notes_len * 2
            unanimous[category].append((score, aid, rec, notes))

    # Sort each category by score descending
    for cat in unanimous:
        unanimous[cat].sort(key=lambda x: -x[0])

    return unanimous


def print_sample(idx, rec, notes, category):
    sep = "=" * 72
    print(f"\n{sep}")
    print(f"  SAMPLE {idx}  |  Error Category: {category}: {CATEGORY_NAMES[category]}")
    print(f"  Model: {rec['model']}  |  ID: {rec['annotation_id']}")
    print(f"  Difficulty: {rec.get('difficulty','')}  |  Domain: {rec.get('domain','')}")
    print(sep)
    print("\n[PROBLEM STATEMENT]")
    print(rec.get("problem_statement", "").strip())
    print("\n[GOLD ANSWER (excerpt, first 600 chars)]")
    print(rec.get("gold_answer", "").strip()[:600])
    print("\n[MODEL REASONING (excerpt, first 1200 chars)]")
    print(rec.get("model_reasoning", "").strip()[:1200])
    print(f"\n[TRIGGERING TEXT SPAN]")
    print(rec.get("triggering_text_span", ""))
    print(f"\n[ANNOTATOR NOTES]")
    for ann, note in zip(ANNOTATORS, notes):
        print(f"  {ann}: {note}")
    print()


OUT_FILE = os.path.join(os.path.dirname(BASE_DIR), "error_analysis_results_metrics",
                        "representative_samples.txt")


def main():
    by_id = load_all()
    candidates = find_candidates(by_id)

    preferred_cats = [3, 4, 5, 2, 6, 1]
    chosen = []
    used_cats = set()
    for cat in preferred_cats:
        if cat in candidates and candidates[cat] and cat not in used_cats:
            chosen.append((cat, candidates[cat][0]))
            used_cats.add(cat)
        if len(chosen) == 3:
            break

    import io
    buf = io.StringIO()

    def w(s=""):
        buf.write(s + "\n")

    sep = "=" * 72
    w(sep)
    w("  THREE REPRESENTATIVE SAMPLES FOR THE PAPER")
    w("  (unanimous agreement across all 3 annotators)")
    w(sep)

    for i, (cat, (score, aid, rec, notes)) in enumerate(chosen, 1):
        w()
        w(sep)
        w(f"  SAMPLE {i}  |  Error Category {cat}: {CATEGORY_NAMES[cat]}")
        w(f"  Model     : {rec['model']}")
        w(f"  Item ID   : {rec['annotation_id']}")
        w(f"  Difficulty: {rec.get('difficulty','')}  |  Domain: {rec.get('domain','')}")
        w(sep)
        w()
        w("[PROBLEM STATEMENT]")
        w(rec.get("problem_statement", "").strip())
        w()
        w("[GOLD ANSWER]")
        w(rec.get("gold_answer", "").strip())
        w()
        w("[MODEL REASONING]")
        w(rec.get("model_reasoning", "").strip())
        w()
        w("[TRIGGERING TEXT SPAN]")
        w(f"  \"{rec.get('triggering_text_span', '')}\"")
        w()
        w("[ANNOTATOR NOTES — all three independent annotators]")
        for ann, note in zip(ANNOTATORS, notes):
            w(f"  {ann}: {note}")
        w()

    content = buf.getvalue()
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    with open(OUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Saved: {OUT_FILE}")

    # Also print category distribution
    print("\nUnanimous sample counts by category:")
    for cat in sorted(candidates):
        print(f"  {cat}: {CATEGORY_NAMES[cat]:30s} -> {len(candidates[cat])}")


if __name__ == "__main__":
    main()
