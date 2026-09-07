"""Independent ground truth for the 61 archived Phase 4 traces.

**This module must not import ``tests.comparators``.**  A label set produced by
running the comparator under test would be a test carrying its own answer key
(D-034) -- it would agree with the comparator by construction and measure
nothing.  Every label here is derived from the *question* and the *gold answer*
by code written against the item's own definition, or is a hand reading of the
model's final answer recorded with the sentence that justifies it.

Two label kinds:

``computed``  The correct answer is recomputed from the question text (the
              sequence items) or from the gold string (the categorical items),
              and the model's stated answer is parsed by a **separate** parser
              written for this file alone.  Independent in the sense R1.2 means.
``read``      The model's answer form defeats mechanical parsing and the label
              is a hand reading.  Each carries the quoted answer text so a
              reviewer can check the reading rather than trust it.

Run: ``python -m tests.comparators.ground_truth``
"""

from __future__ import annotations

import glob
import json
import re
from typing import Any

ARCHIVE = "error_analysis_annotation/samples/*.jsonl"

TEMPLATES = (
    "signal_operations",
    "system_properties_memory_causality",
    "system_property_linearity",
    "incompressible_continuity",
)


def load(stem: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(glob.glob(ARCHIVE)):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                if r["question_id"].rsplit("__", 1)[0] == stem:
                    rows.append(r)
    return rows


# --------------------------------------------------------------------------
# signal_operations -- fully computable from the question
# --------------------------------------------------------------------------


def _parse_marked_seq(text: str) -> tuple[list[int], int | None]:
    m = re.search(r"\{([^}]*)\}", text)
    if not m:
        raise ValueError(f"no braced sequence in {text[:80]!r}")
    vals, origin = [], None
    for i, raw in enumerate(x.strip() for x in m.group(1).split(",")):
        if raw.startswith("*") and raw.endswith("*"):
            origin, raw = i, raw.strip("*")
        vals.append(int(raw))
    return vals, origin


def signal_truth(row: dict[str, Any]) -> dict[int, int]:
    """The correct index->value map, derived from the question alone."""
    q = row["problem_statement"]
    mx = re.search(r"x\[n\]\s*=\s*(\{[^}]*\})", q)
    vals, origin = _parse_marked_seq(mx.group(1))
    x = {i - origin: v for i, v in enumerate(vals)}
    if re.search(r"=\s*x\[-n\]", q):
        return {-k: v for k, v in x.items()}
    n0 = int(re.search(r"x\[n - \(?(-?\d+)\)?\]", q).group(1))
    return {k + n0: v for k, v in x.items()}


def gold_seq(row: dict[str, Any]) -> tuple[list[int], int | None]:
    g = row["gold_answer"]
    return _parse_marked_seq(g[g.rfind("**Answer:**"):])


# --------------------------------------------------------------------------
# The labels.
#
# `correct` is whether the MODEL's final answer is the right answer to the
# question.  It is the property a comparator should reproduce, and it is what
# precision and recall in phase4_summary.md S4 are computed against.
# --------------------------------------------------------------------------

#: signal_operations, 16 traces.  `answer` quotes the model's final answer
#: verbatim; `reading` is the index->value map it asserts, taken from the
#: answer *and* the surrounding step lines where the answer alone is silent
#: about the origin.
SIGNAL_LABELS: dict[int, tuple[bool, str]] = {
    0:  (True,  "{*0*, -1, 5, -4, 3}: y[0]=0, y[1..4]=-1,5,-4,3 == truth {1:-1,2:5,3:-4,4:3} with one pad"),
    1:  (False, "{6,0,0,0,0} for n=0..4 plus z[-1..-4]: asserts z[-4]=2..z[0]=6; truth is z[-3]=2..z[1]=6 -- origin off by one"),
    2:  (False, "{-7, 7, 5, -1, 10, 4, -1} is truth reversed; the trace computes z[k]=x[k], not x[-k]"),
    3:  (False, "{-3,5,-4,-5,6,-9,3} for n=-6..0; truth is n=-2..4 with a different map"),
    4:  (False, "{-8, 4, 9, -3}: four values where the truth has six"),
    5:  (False, "{-8, 9, -10, -10, -7, -4, 6} is a rotation of the truth"),
    6:  (False, "{0,0,0,-6,4,2,-5,-2}: values right, but the steps place y[3]=-6 where truth has y[0]=-6"),
    7:  (False, "{0, 0, 0, -1}: one non-zero value where the truth has four"),
    8:  (False, "{3,-9,6,-5,-4,5,-3} is a rotation of the truth"),
    9:  (False, "{0,0,0,0,0,-8,-2,-9,-6,6,7} is the truth reversed and zero-padded"),
    10: (False, "{0,-2,5,5,10,-3,2}: values right; steps place z[1]=-2 where truth has z[-5]=-2"),
    11: (False, "{-7,-4,6,9,-8,0,0}: drops the leading -10,-10 and pads at the tail"),
    12: (False, "{9, 4, 4, 10} is a rotation of the truth {4,4,10,9}"),
    13: (False, "{0,0,2,1,2,6}: values right; steps place y[2]=2 where truth has y[0]=2"),
    14: (False, "{-4, 3, 5, -1} is a rotation of the truth"),
    15: (True,  "{-10,-10,-7,-4,6,9,-8} equals gold exactly; gold states no origin either"),
}

#: system_properties_memory_causality, 24 traces.  The label is the
#: (memoryless, causal) pair the model's final answer asserts, compared with
#: gold's pair.  Read by hand from the answer text.
MC_LABELS: dict[int, tuple[bool, str]] = {
    0:  (True,  "NOT memoryless / NOT causal == gold (No, No)"),
    1:  (True,  "NOT memoryless / NOT causal == gold (No, No)"),
    2:  (True,  "'not memoryless but is causal' == gold (No, Yes)"),
    3:  (True,  "'a) Not memoryless b) Not causal' == gold (No, No)"),
    4:  (True,  "'Not memoryless (it depends on the past input) / Causal' == gold (No, Yes)"),
    5:  (True,  "'a) Not memoryless; b) Causal' == gold (No, Yes)"),
    6:  (True,  "'a) Not memoryless; b) Causal' == gold (No, Yes)"),
    7:  (True,  "'a) Not memoryless, b) Not causal' == gold (No, No)"),
    8:  (True,  "'a) Yes, b) Yes' == gold (Yes, Yes); unlabelled, positional"),
    9:  (True,  "'a) Not memoryless, b) Causal' == gold (No, Yes)"),
    10: (True,  "'a) Not memoryless, b) Causal' == gold (No, Yes)"),
    11: (True,  "'a) Not memoryless, b) Causal' == gold (No, Yes)"),
    12: (True,  "'a) Not memoryless, b) Causal' == gold (No, Yes)"),
    13: (True,  "'a) Not memoryless, b) Not causal' == gold (No, No)"),
    14: (True,  "'a) Not Memoryless b) Causal' == gold (No, Yes)"),
    15: (True,  "'neither memoryless nor causal' == gold (No, No)"),
    16: (True,  "'Memoryless: No; Causal: Yes' == gold"),
    17: (True,  "'Memoryless: No; Causal: Yes' == gold"),
    18: (True,  "'Memoryless: No; Causal: No' == gold"),
    19: (True,  "'Memoryless: No; Causal: No' == gold"),
    20: (True,  "'Memoryless: No; Causal: Yes' == gold"),
    21: (True,  "'a) Not memoryless; b) Causal' == gold (No, Yes)"),
    22: (True,  "'a) Not memoryless; b) Causal' == gold (No, Yes)"),
    23: (True,  "'neither memoryless nor causal' == gold (No, No)"),
}

#: system_property_linearity, 15 traces.
LIN_LABELS: dict[int, tuple[bool, str]] = {
    0:  (True,  "'**NOT LINEAR**' == gold 'not linear' (also calls it 'incrementally linear' mid-derivation)"),
    1:  (True,  "'Not linear' == gold"),
    2:  (True,  "'The system is nonlinear.' == gold 'not linear'"),
    3:  (True,  "'Not linear' == gold"),
    4:  (True,  "'The system is not linear.' == gold"),
    5:  (True,  "'The system is linear.' == gold 'linear'"),
    6:  (True,  "'Non-linear' == gold 'not linear'"),
    7:  (True,  "'Linear' == gold 'linear'"),
    8:  (True,  "'Linear' == gold 'linear'"),
    9:  (True,  "'Nonlinear' == gold 'not linear'"),
    10: (True,  "'The system is not linear.' == gold"),
    11: (True,  "'The system is not linear.' == gold"),
    12: (True,  "'... is not linear because it fails both ...' == gold"),
    13: (True,  "'... is not linear because it does not satisfy ...' == gold"),
    14: (True,  "'The system is not linear because it fails both tests.' == gold"),
}

#: incompressible_continuity, 6 traces.
CONT_LABELS: dict[int, tuple[bool, str]] = {
    0:  (True,  "6xy - 2x^2 == gold -2x^2 + 6xy (terms reordered, Unicode minus and superscript)"),
    1:  (True,  "-2x^2 - 10xy == gold exactly (Unicode minus and superscript)"),
    2:  (True,  "x^2 + 6xy == gold exactly (Unicode superscript)"),
    3:  (True,  "3x^2/2 - 6xy + C(y) == gold 1.5x^2 - 6xy, plus the retained arbitrary function"),
    4:  (False, "10y^2 + 4xy: integrated with respect to the wrong variable"),
    5:  (False, "u(x, y) = C: the whole solution is dropped"),
}

LABELS = {
    "signal_operations": SIGNAL_LABELS,
    "system_properties_memory_causality": MC_LABELS,
    "system_property_linearity": LIN_LABELS,
    "incompressible_continuity": CONT_LABELS,
}


def verify_signal_labels() -> list[str]:
    """Cross-check the ``signal_operations`` labels against the computed truth.

    Every label whose reading is "values equal the truth" must have a model
    answer whose value multiset matches, and every "rotation"/"reversal" label
    must not match in order.  This is a check on *this file*, not on the
    comparator, and it is here because a hand-written label table is exactly
    the artefact D-034 warns about.
    """
    problems: list[str] = []
    rows = load("signal_operations")
    for i, row in enumerate(rows):
        truth = signal_truth(row)
        gvals, gorigin = gold_seq(row)
        gmap = {j - gorigin: v for j, v in enumerate(gvals)} if gorigin is not None else None
        if gmap is not None and gmap != truth:
            problems.append(f"signal {i}: GOLD disagrees with the question: {gmap} != {truth}")
        if gmap is None and sorted(gvals) != sorted(truth.values()):
            problems.append(f"signal {i}: gold values {gvals} != truth values {sorted(truth.values())}")
    return problems


def gold_states_origin_rate(n: int = 4000) -> tuple[int, int]:
    """How often does gold's printed answer state where ``n = 0`` is?

    The template marks the origin only when the *result's* support contains
    ``n = 0``.  When a shift moves the support off the origin, the printed
    answer is a bare value list -- from which the origin cannot be recovered.
    That is a property of the item worth quantifying rather than asserting, and
    it is what D-050 records.
    """
    import importlib
    import random

    mod = importlib.import_module(
        "data.templates.branches.electrical_engineering.signals_and_systems.discrete_time_signals"
    )
    stated = 0
    for seed in range(n):
        random.seed(seed)
        _, sol = mod.template_signal_operations()
        block = sol[sol.rfind("**Answer:**"):]
        if re.search(r"\*-?\d+\*", block):
            stated += 1
    return stated, n


def main() -> None:
    print("Ground-truth label census")
    print("=" * 64)
    total_correct = 0
    total = 0
    for stem in TEMPLATES:
        rows = load(stem)
        labels = LABELS[stem]
        if len(rows) != len(labels):
            print(f"  !! {stem}: {len(rows)} traces but {len(labels)} labels")
        n_ok = sum(1 for ok, _ in labels.values() if ok)
        total_correct += n_ok
        total += len(labels)
        print(f"  {stem:36s} n={len(labels):3d}  correct={n_ok:3d}  wrong={len(labels)-n_ok:3d}")
    print(f"  {'TOTAL':36s} n={total:3d}  correct={total_correct:3d}  wrong={total-total_correct:3d}")

    print()
    print("Self-check: gold vs the question, signal_operations")
    problems = verify_signal_labels()
    print("  clean" if not problems else "\n".join("  " + p for p in problems))

    print()
    stated, n = gold_states_origin_rate()
    print(f"Gold states the origin on {stated}/{n} instances "
          f"({stated / n:.1%}); it is silent on {n - stated} ({1 - stated / n:.1%})")


if __name__ == "__main__":
    main()
