"""D4.6 -- the comparator for D3.4 |S|7, and the first implementation of it.

Phase 3 shipped |S|7 unexercised and said so: every verifier it built implements
|S|6, gold checking, and **the rules that actually decide a model's score had no
corpus at all** (``phase3_summary.md`` |S|13.3, both schema reviewers
independently).  This module implements |S|7.1-|S|7.4 and
``build_conformance.py`` supplies the corpus.

The shape of the answer, and why it is two verdicts rather than one:

    ``answer``   did the candidate reach gold's answer?
    ``process``  did it execute the stated scheme?

|S|7.2 and |S|7.3 exist precisely because these come apart, and in *opposite*
directions for the two node types:

* ``iteration`` (``cardinality: incidental``) -- a correct answer reached in
  four updates instead of three is correct.  The count is an observation.
* ``decision`` (``cardinality: answer_bearing``) -- a different element count
  **is** a different answer, because ``n`` is the answer.  But two different
  valid packings can share an ``n``, so the answer can be right while the
  process is wrong.

A verifier built from one merged node type gets one of these wrong, which is
D-038 and the finding that shaped the whole schema.

Usage:
    python -m tests.trace_schema.candidate_7 <candidates.json>
"""

from __future__ import annotations

import json
import sys
from decimal import Decimal
from typing import Any

from tests.trace_schema.reviewer_d2_verifier_15 import verify_node

MATCH = "MATCH"
MISMATCH = "MISMATCH"
UNRESOLVED = "UNRESOLVED"


class NodeVerdict:
    """Answer credit and process credit, kept separate on purpose."""

    def __init__(self, label: str = "") -> None:
        self.label = label
        self.answer = MATCH
        self.process = MATCH
        self.answer_reason = ""
        self.process_reason = ""
        self.observations: list[str] = []
        self.clauses: list[str] = []

    def fail_answer(self, clause: str, reason: str) -> None:
        if self.answer == MATCH:
            self.answer, self.answer_reason = MISMATCH, reason
        self.clauses.append(clause)

    def fail_process(self, clause: str, reason: str) -> None:
        if self.process == MATCH:
            self.process, self.process_reason = MISMATCH, reason
        self.clauses.append(clause)

    def unresolve_answer(self, clause: str, reason: str) -> None:
        if self.answer == MATCH:
            self.answer, self.answer_reason = UNRESOLVED, reason
        self.clauses.append(clause)

    def observe(self, note: str) -> None:
        self.observations.append(note)

    def __str__(self) -> str:  # pragma: no cover - display
        return (f"answer={self.answer:10s} process={self.process:10s} "
                f"clauses={sorted(set(self.clauses)) or '-'}")


# --------------------------------------------------------------------------
# |S|7.1 -- numeric tolerance, always from the GOLD node
# --------------------------------------------------------------------------


def tolerance(gold: dict, symbol: str, dp: int | None = None) -> Decimal:
    """Half a unit in the last place **gold** displays ``symbol``.

    ``symbol_precision`` is read from ``gold`` and never from the candidate.
    D3.4 |S|7.1 states this and Reviewer D2 asked for it to be said out loud,
    because a careless implementer reads "the last place it is displayed" as
    the candidate's rendering -- and then a candidate that declares coarser
    precisions widens the tolerance it is judged against.  That is a candidate
    grading itself, and ``C10`` in the corpus is the trace that does it.
    """
    if dp is None:
        prec = gold.get("symbol_precision", {})
        if symbol not in prec:
            raise KeyError(f"8A.9: {symbol!r} absent from the gold symbol_precision")
        dp = int(prec[symbol])
    return Decimal("0.5") * (Decimal(10) ** -dp)


def within(gold_val: Any, cand_val: Any, tol: Decimal) -> bool:
    """Boundary **inclusive**, per |S|7.1."""
    return abs(Decimal(str(cand_val)) - Decimal(str(gold_val))) <= tol


# --------------------------------------------------------------------------
# The comparator
# --------------------------------------------------------------------------


def compare(gold: dict, cand: dict, label: str = "") -> NodeVerdict:
    """Score ``cand`` against ``gold`` under D3.4 |S|7."""
    v = NodeVerdict(label)

    if gold.get("node_type") != cand.get("node_type"):
        v.fail_answer("7.0", f"node_type {cand.get('node_type')!r} != gold "
                             f"{gold.get('node_type')!r}")
        v.fail_process("7.0", "node types differ")
        return v

    # A candidate must still be a well-formed node.  |S|6.0 says the mode is the
    # caller's, so it is checked in `candidate` mode -- under which a budget
    # overrun and `satisfied: false` are observations rather than failures,
    # and everything structural still binds.
    rep = verify_node(cand, label=label or "candidate", mode="candidate")
    structural = [f for f in getattr(rep, "failures", [])]
    if structural:
        for f in structural:
            v.fail_process(str(getattr(f, "clause", "6")), str(f))

    _result_credit(gold, cand, v)

    if gold["cardinality"] == "incidental":
        _incidental(gold, cand, v)
    else:
        _answer_bearing(gold, cand, v)
    return v


def _result_credit(gold: dict, cand: dict, v: NodeVerdict) -> None:
    """|S|7.1 + |S|7.4 -- does the candidate's result equal gold's?"""
    g_res, c_res = gold.get("result"), cand.get("result")
    if g_res is None:
        return
    if c_res is None:
        v.unresolve_answer("7.1", "candidate carries no result")
        return
    dp = g_res.get("dp")
    tol = tolerance(gold, "result", dp=dp)

    # |S|7.4 -- when the question prescribes a method, the method's own stated
    # tolerance governs if it is the looser of the two.  Marking a direct
    # solver wrong for one display unit is the same error as marking a
    # four-update solver wrong for taking four updates, and |S|7.2 already
    # refuses the second.
    method_tol = gold.get("termination", {}).get("tolerance")
    if method_tol is not None and gold["termination"].get("kind") == "convergence":
        mt = Decimal(str(method_tol))
        if mt > tol:
            v.observe(f"7.4: tolerance widened {tol} -> {mt} by the prescribed method")
            tol = mt

    if within(g_res["value"], c_res["value"], tol):
        v.observe(f"result {c_res['value']} within {tol} of gold {g_res['value']}")
    else:
        v.fail_answer("7.1", f"result {c_res['value']} != gold {g_res['value']} "
                             f"(tolerance {tol})")


# --------------------------------------------------------------------------
# |S|7.2 -- cardinality: incidental
# --------------------------------------------------------------------------


def _incidental(gold: dict, cand: dict, v: NodeVerdict) -> None:
    g_els, c_els = gold["elements"], cand["elements"]

    # Row 1: a different count is NOT a defect.  Reported, never scored.
    if len(g_els) != len(c_els):
        v.observe(f"7.2: element count gold {len(g_els)}, candidate {len(c_els)} "
                  f"(incidental; not scored)")

    # Row 5: a budget overrun is reported, not failed.  The budget is the
    # generator's guard; a solver taking six updates to the same depth has not
    # made an error.
    max_el = gold.get("termination", {}).get("max_elements")
    if max_el is not None and len(c_els) > int(max_el):
        v.observe(f"7.2: candidate used {len(c_els)} elements against a "
                  f"max_elements of {max_el} (reported, not failed)")

    # Row 4: the termination predicate must hold on the candidate's LAST
    # element.  A candidate that stopped early with a change above tolerance
    # has not executed the stated scheme, even if its answer happens to match.
    term = cand.get("termination", {})
    if term.get("satisfied") is not True:
        v.fail_process("7.2", "candidate's termination predicate does not hold on "
                              "its last element (early stop)")
    else:
        roles = cand.get("roles", {})
        qty = roles.get(term.get("quantity_role", ""), "")
        if c_els and qty in c_els[-1]:
            tol = Decimal(str(gold["termination"]["tolerance"]))
            last = abs(Decimal(str(c_els[-1][qty])))
            if not last < tol:
                v.fail_process("7.2", f"termination quantity {last} is not below "
                                      f"tolerance {tol} on the last element")

    # Row 6: the preamble is part of the QUESTION.  A candidate that starts
    # from a different bracket has solved a different problem, however well.
    g_pre, c_pre = gold.get("preamble"), cand.get("preamble")
    if g_pre is not None and c_pre is not None:
        if len(g_pre) != len(c_pre):
            v.fail_process("7.2", f"preamble has {len(c_pre)} frames, gold has {len(g_pre)}")
        else:
            for i, (gf, cf) in enumerate(zip(g_pre, c_pre)):
                for sym in gold.get("evaluation_symbols", []):
                    if sym not in gf or sym not in cf:
                        continue
                    if not within(gf[sym], cf[sym], tolerance(gold, sym)):
                        v.fail_process(
                            "7.2", f"preamble frame {i} symbol {sym}: {cf[sym]} != "
                                   f"gold {gf[sym]}; the starting trials are part of the question")

    # Row 3: process credit requires every candidate element to satisfy |S|4.3,
    # **including |S|4.3.9's frame check**.  Through 1.3 that check did not
    # exist, so a candidate emitting arithmetically consistent nonsense took
    # full process credit -- the mechanism rewarding what it was built to
    # detect.  It is carried by the structural pass in `compare`, and this is
    # the observation that says so.
    v.observe("7.2: process credit rests on the S4.3 local invariants, S4.3.9 included")


# --------------------------------------------------------------------------
# |S|7.3 -- cardinality: answer_bearing
# --------------------------------------------------------------------------


def _answer_bearing(gold: dict, cand: dict, v: NodeVerdict) -> None:
    g_els, c_els = gold["elements"], cand["elements"]
    sym = gold.get("cardinality_symbol", "n")

    # Row 1: a different element count IS a different answer, and it is
    # reported as `cardinality_mismatch` on `cardinality_symbol` -- not as a
    # formatting or style difference.
    if len(g_els) != len(c_els):
        v.fail_answer("7.3", f"cardinality_mismatch on {sym}: gold {len(g_els)}, "
                             f"candidate {len(c_els)}")
        v.fail_process("7.3", "element counts differ; the rule replay cannot align")
        return

    g_roles, c_roles = gold["roles"], cand["roles"]
    g_key, c_key = g_roles["committed"], c_roles["committed"]

    # Row 2: same count, different `committed` sets -- process failure, and the
    # answer may still be correct, because two different valid packings can
    # share an `n`.  Score the answer on `n`; score the process on the replay.
    # Row 3: compared as SETS for milestone credit.
    for i, (ge, ce) in enumerate(zip(g_els, c_els)):
        if set(ge[g_key]) != set(ce[c_key]):
            v.fail_process(
                "7.3", f"element {i}: committed set {sorted(ce[c_key])} != gold "
                       f"{sorted(ge[g_key])}; a different packing with the same n")
            break

    # Row 3 continued: as ORDERED lists when replaying `selection`, whose
    # output is ordered.  Reported separately -- a set match with an order
    # mismatch is a weaker process failure than a set mismatch, and collapsing
    # them loses that.
    for i, (ge, ce) in enumerate(zip(g_els, c_els)):
        if set(ge[g_key]) == set(ce[c_key]) and list(ge[g_key]) != list(ce[c_key]):
            v.observe(f"7.3: element {i} commits gold's set in a different order "
                      f"({list(ce[c_key])} vs {list(ge[g_key])})")

    # Row 4: a trial violating |S|5.4.4 is a process failure located at
    # (station_id, trial index).  Carried by the structural pass; located here.
    for f in getattr(getattr(v, "_rep", None), "failures", []) or []:
        if "5.4.4" in str(f):
            v.fail_process("7.3", f"selection violated: {f}")


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    with open(argv[1], encoding="utf-8") as fh:
        corpus = json.load(fh)

    width = max(len(c["id"]) for c in corpus)
    ok = bad = 0
    print(f"{'case':{width}s}  {'answer':>10s} {'process':>10s}   expected     result")
    print("-" * (width + 56))
    for case in corpus:
        v = compare(case["gold"], case["candidate"], label=case["id"])
        exp_a, exp_p = case["expect_answer"], case["expect_process"]
        good = (v.answer == exp_a) and (v.process == exp_p)
        ok, bad = ok + good, bad + (not good)
        print(f"{case['id']:{width}s}  {v.answer:>10s} {v.process:>10s}   "
              f"{exp_a[:5]}/{exp_p[:5]:<6s} {'ok' if good else 'MISMATCH'}")
        if not good:
            print(f"{'':{width}s}    expected {exp_a}/{exp_p}; "
                  f"answer_reason={v.answer_reason!r} process_reason={v.process_reason!r}")
    print("-" * (width + 56))
    print(f"{ok} as expected, {bad} not")

    # |S|7.1/|S|7.4 tolerance semantics, exercised directly.  F7-3 records why
    # one of them cannot be reached through a node at all.
    from tests.trace_schema.build_conformance_7 import FINDINGS, TOLERANCE_CASES

    print("\nS7.1 / S7.4 tolerance semantics, at the function level")
    print("-" * (width + 56))
    t_bad = 0
    for tid, gv, dp, cv, mt, expected in TOLERANCE_CASES:
        tol = Decimal("0.5") * (Decimal(10) ** -dp)
        if mt is not None and Decimal(str(mt)) > tol:
            tol = Decimal(str(mt))
        got = within(gv, cv, tol)
        t_bad += got != expected
        print(f"  {tid:38s} gold {gv} dp{dp} cand {cv} tol {tol}  "
              f"{'accept' if got else 'reject':6s} {'ok' if got == expected else 'MISMATCH'}")
    print(f"  {len(TOLERANCE_CASES) - t_bad}/{len(TOLERANCE_CASES)} as expected")

    covered = set()
    findings = set()
    for case in corpus:
        covered.update(case.get("exercises", []))
        if case.get("finding"):
            findings.add(case["finding"])
    missing = sorted(set(SEVEN_DISPOSITIONS) - covered)
    print(f"\nS7 dispositions exercised: {len(covered)}/{len(SEVEN_DISPOSITIONS)}")
    print("NOT exercised: " + (", ".join(missing) if missing else "none -- S7 is fully covered"))

    if findings:
        print("\n" + "=" * (width + 56))
        print("FINDINGS AGAINST D3.4 S7 -- the specification, not the verifier")
        print("=" * (width + 56))
        for fid in sorted(findings):
            print(f"\n{fid}:")
            body = FINDINGS[fid]
            line = ""
            for word in body.split():
                if len(line) + len(word) > 74:
                    print("  " + line)
                    line = word
                else:
                    line = (line + " " + word).strip()
            if line:
                print("  " + line)
    return 0 if bad == 0 and t_bad == 0 and not missing else 1


#: Every disposition D3.4 |S|7 states, enumerated so the corpus can be checked
#: for coverage rather than assumed to have it.  A conformance corpus that does
#: not know what it is meant to cover is a list of remembered cases (D4.11).
SEVEN_DISPOSITIONS = (
    "7.1-tolerance-inclusive",
    "7.1-tolerance-exceeded",
    "7.1-gold-precision-governs",
    "7.2-count-differs-not-scored",
    "7.2-answer-credit",
    "7.2-process-credit",
    "7.2-early-stop",
    "7.2-budget-overrun-reported",
    "7.2-preamble-differs",
    "7.2-fabricated-frame",
    "7.3-cardinality-mismatch",
    "7.3-same-n-different-packing",
    "7.3-order-within-set",
    "7.3-trial-violates-selection",
    "7.4-method-tolerance-governs",
)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
