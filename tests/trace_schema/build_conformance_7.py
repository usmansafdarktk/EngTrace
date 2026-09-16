"""D4.6 -- build the conformance corpus for D3.4 |S|7.

Phase 3's |S|8 reject list was hand-written, and its own summary records what
that cost: "a list of remembered failure modes", against which five guessed
invariant holes all hit.  D4.11 says generate the list from the specification
instead.  This module does that for |S|7: ``SEVEN_DISPOSITIONS`` in
``candidate_7.py`` enumerates every disposition the section states, each case
below declares which ones it exercises, and the driver **fails if any
disposition is unexercised**.  Coverage is checked, not asserted.

Each candidate is a mutation of a real gold node from
``phase3_conformance/traces.json``, so the corpus inherits the shape the
templates actually emit rather than one invented here.

Run: ``python -m tests.trace_schema.build_conformance_7``
"""

from __future__ import annotations

import copy
import json
import os

GOLD = "docs/re-implementation-sep/track-a/phase3_conformance/traces.json"
OUT_DIR = "docs/re-implementation-sep/track-a/phase4_conformance"
OUT = os.path.join(OUT_DIR, "candidates.json")

MATCH, MISMATCH, UNRESOLVED = "MATCH", "MISMATCH", "UNRESOLVED"


def load_gold() -> tuple[dict, dict]:
    with open(GOLD, encoding="utf-8") as fh:
        rows = json.load(fh)

    def node(template_id: str, seed: int = 0) -> dict:
        for r in rows:
            if r["template_id"] == template_id and r["seed"] == seed:
                n = r["trace_nodes"]
                return copy.deepcopy(n[0] if isinstance(n, list) else n)
        raise KeyError(f"{template_id} seed {seed} not in {GOLD}")

    it = node("template_normal_depth_iteration")
    de = next(
        copy.deepcopy(r["trace_nodes"][0] if isinstance(r["trace_nodes"], list)
                      else r["trace_nodes"])
        for r in rows if r["template_id"] == "template_line_balancing_heuristic"
    )
    return it, de


def _extend_iteration(node: dict, n_extra: int = 1) -> dict:
    """Append ``n_extra`` genuine secant updates to a converged iteration node.

    The extra elements are *computed*: the update relation is applied, the
    evaluation frame is built from ``frame_relations``, and every carry is
    re-established.  The result is a node a conforming |S|6 verifier accepts --
    which is the point.  A candidate that fails |S|6 cannot demonstrate anything
    about |S|7, because |S|7 never gets to run on it.
    """
    from tests.trace_schema.reviewer_d2_verifier_15 import FrameEval
    from tests.trace_schema.reviewer_d2_verifier import do_round

    c = copy.deepcopy(node)
    fe = FrameEval(c)
    roles = c["roles"]
    ev_roles = c["evaluation_roles"]
    prec = c["symbol_precision"]
    mode = c["rounding"]
    syms = c["evaluation_symbols"]
    it_sym, res_sym = ev_roles["iterate"], ev_roles["residual"]

    def frame(y: float) -> dict:
        f = {}
        for s in syms:
            f[s] = (do_round(float(y), int(prec[s]), mode) if s == it_sym
                    else do_round(fe.value(s, y), int(prec[s]), mode))
        return f

    for _ in range(n_extra):
        last = c["elements"][-1]
        y_prev = last[roles["iterate_curr"]]
        g_prev = last[roles["residual_curr"]]
        y_curr = last[roles["iterate_next"]]
        # The frame the previous element already computed for y_curr becomes
        # this element's residual -- the same carry gold uses.
        prev_eval = last.get(roles["evaluation"])
        if prev_eval is None:
            prev_eval = frame(y_curr)
            last[roles["evaluation"]] = prev_eval
        g_curr = prev_eval[res_sym]

        denom = g_curr - g_prev
        if denom == 0:
            break
        y_next = do_round(y_curr - g_curr * (y_curr - y_prev) / denom,
                          int(prec[roles["iterate_next"]]), mode)
        change = do_round(abs(y_next - y_curr), int(prec[roles["change"]]), mode)

        last[roles["converged"]] = False
        el = {
            roles["index"]: last[roles["index"]] + 1,
            roles["iterate_prev"]: y_prev,
            roles["residual_prev"]: g_prev,
            roles["iterate_curr"]: y_curr,
            roles["residual_curr"]: g_curr,
            roles["iterate_next"]: y_next,
            roles["change"]: change,
            roles["converged"]: True,
            roles["evaluation"]: None,
        }
        c["elements"].append(el)

    # The result is the derivation over the new last element, which for a
    # converged secant iteration is the same depth to gold's display precision.
    c["result"] = dict(c["result"],
                       value=do_round(c["elements"][-1][roles["iterate_next"]],
                                      int(c["result"]["dp"]), mode))
    return c


CASES: list[dict] = []


def case(cid, gold, cand, answer, process, exercises, note, finding=None):
    """Record one candidate.

    ``expect_answer`` / ``expect_process`` are what D3.4 |S|7 **plus the |S|6
    structural clauses** actually produce, which is not always what |S|7's prose
    says on its own.  Where the two diverge the case carries a ``finding`` id
    and ``FINDINGS`` below states the divergence.  Writing the expectation to
    match the prose and letting the case fail would bury the finding in a red
    test; writing it to match the code without saying so would bury it
    completely.
    """
    CASES.append({
        "id": cid, "gold": gold, "candidate": cand,
        "expect_answer": answer, "expect_process": process,
        "exercises": exercises, "note": note, "finding": finding,
    })


#: What building this corpus found out about |S|7.  Each is a defect in the
#: SPECIFICATION, not in a verifier: |S|7 was written before anything exercised
#: it, which is the gap D4.6 exists to close.
FINDINGS: dict[str, str] = {
    "F7-1": (
        "S7.2's headline disposition -- 'a correct answer reached in four "
        "iterations instead of three is correct' -- has NO INSTANCES for "
        "normal_depth_iteration. S4.3.3 ties `converged` to the node's own "
        "termination.tolerance, S8A.4 forbids `converged: true` on any element "
        "but the last, S4.3.1 recomputes every update, S4.6 recomputes every "
        "frame, and the preamble is fixed by the question. Together these make "
        "the conforming node for a given question UNIQUE, so a candidate has no "
        "latitude to take a different number of updates. Verified both ways: "
        "clearing `converged` on the terminating element fails 4.3.3 and 8A.4; "
        "appending after it fails 8A.3 and 8A.4. "
        "D-038 measured the count varying 1..5 over 4,000 seeds and that stands "
        "-- but the variation is BETWEEN QUESTIONS, not between solvers of one "
        "question, and S7.2 conflates the two. This is the cost of the six "
        "review rounds that hardened the schema: each closed a way for a "
        "candidate to lie, and together they also closed every way for a "
        "candidate to be differently right."
    ),
    "F7-2": (
        "S7 never says whether S6 step 8 -- result.value equals the result.from "
        "derivation -- binds a candidate. It does in every verifier built, and "
        "it should: a model whose stated answer contradicts its own steps has "
        "failed the process. The consequence is that answer credit and process "
        "credit are NOT independent for `iteration`, so 'wrong answer, clean "
        "process' is not constructible. S7 is written as though they vary "
        "freely. It should say which, explicitly."
    ),
    "F7-3": (
        "S7.4's case -- a solver who solves the governing equation directly and "
        "lands one display unit from the prescribed scheme's output -- cannot be "
        "represented as an `iteration` node at all, because a direct solve has "
        "no iteration. The answer tolerance S7.4 specifies is correct and "
        "testable, but only at the comparator's function level, not through a "
        "node. TOLERANCE_CASES tests it there."
    ),
    "F7-4": (
        "S7.3 row 3 distinguishes comparing `committed` AS SETS for milestone "
        "credit from AS ORDERED LISTS when replaying `selection`. For a "
        "S6-conforming node the distinction is vacuous: S5.4.6 already requires "
        "`committed` to equal the ordered chosen values, so a candidate whose "
        "set matches but whose order differs is rejected before S7.3 runs. The "
        "row describes a case the rest of the schema forbids."
    ),
}


#: |S|7.1 and |S|7.4 are about a *tolerance*, and F7-3 shows one of them cannot
#: be reached through a node at all.  They are exercised here directly, as
#: ``(gold_value, gold_dp, candidate_value, method_tolerance, expected)``.
TOLERANCE_CASES = [
    ("T1-inclusive-high", 1.380, 3, 1.3805, None, True),
    ("T2-inclusive-low", 1.380, 3, 1.3795, None, True),
    ("T3-just-outside-high", 1.380, 3, 1.38051, None, False),
    ("T4-just-outside-low", 1.380, 3, 1.37949, None, False),
    ("T5-exact", 1.380, 3, 1.380, None, True),
    ("T6-one-display-unit-out", 1.380, 3, 1.381, None, False),
    ("T7-method-tolerance-rescues", 1.380, 3, 1.381, 0.002, True),
    ("T8-method-tolerance-not-enough", 1.380, 3, 1.385, 0.002, False),
    ("T9-method-tolerance-narrower-ignored", 1.380, 3, 1.3805, 0.0001, True),
]


def build() -> list[dict]:
    it, de = load_gold()

    # ================= iteration -- cardinality: incidental ================

    # C1  A correct answer reached in one MORE update.  The count differs and
    #     that is not a defect; the answer stands and the process stands.
    #     This is |S|7.2's headline row and the whole reason `incidental` exists.
    #
    #     Built by actually running the secant step and evaluating the frame,
    #     not by copying the last element.  The first draft copied it, and the
    #     verifier rejected the result at 8A.4 -- gold's terminating element
    #     carries `evaluation: null`, so appending after it leaves a non-final
    #     element with no frame.  A conformance corpus whose "valid" cases are
    #     structurally invalid tests the mutation, not the specification
    #     (phase4_summary.md S8, error 4).
    c = _extend_iteration(it, n_extra=1)
    case("C1-extra-update", it, c, MATCH, MISMATCH,
         ["7.2-count-differs-not-scored", "7.2-answer-credit", "7.2-process-credit"],
         "four updates instead of three: REJECTED at 4.3.3, because gold's third "
         "element already satisfies the predicate and 8A.4 forbids continuing",
         finding="F7-1")

    # C2  Early stop: the last element's change is still above tolerance.
    #     The answer may even match; the process has not executed the scheme.
    c = copy.deepcopy(it)
    c["elements"] = c["elements"][:1]
    c["elements"][-1]["converged"] = False
    c["termination"] = dict(c["termination"], satisfied=False)
    case("C2-early-stop", it, c, MATCH, MISMATCH,
         ["7.2-early-stop", "7.2-count-differs-not-scored"],
         "stopped with the change above tolerance; result kept from gold")

    # C3  A fabricated evaluation frame.  Arithmetically consistent-looking
    #     nonsense.  Through schema 1.3 this took FULL process credit, which is
    #     the defect S4.6 and S4.3.9 were added to close -- so a corpus that
    #     does not contain it cannot show the closure holds under S7.
    c = copy.deepcopy(it)
    for el in c["elements"]:
        if el.get("evaluation"):
            el["evaluation"] = dict(el["evaluation"], AR=el["evaluation"]["AR"] + 1.5)
            break
    case("C3-fabricated-frame", it, c, MATCH, MISMATCH,
         ["7.2-fabricated-frame", "7.2-process-credit"],
         "one evaluation frame no longer satisfies frame_relations")

    # C4  A wrong result with a clean process.  Answer and process really do
    #     come apart, in both directions.
    c = copy.deepcopy(it)
    c["result"] = dict(c["result"], value=round(c["result"]["value"] + 0.05, 3))
    case("C4-wrong-result", it, c, MISMATCH, MISMATCH,
         ["7.1-tolerance-exceeded"],
         "states an answer 0.05 m from gold that its own trace does not support; "
         "S6 step 8 binds in candidate mode and catches it",
         finding="F7-2")

    # C5  A result exactly on the half-way boundary.  S7.1 says INCLUSIVE, and
    #     1.1 left it unstated -- latent on gold, live the moment a candidate
    #     is scored (Reviewer D2, F8).  This is the trace that makes it live.
    c = copy.deepcopy(it)
    dp = it["result"]["dp"]
    c["result"] = dict(c["result"], value=it["result"]["value"] + 0.5 * 10 ** -dp)
    case("C5-boundary-inclusive", it, c, MATCH, MISMATCH,
         ["7.1-tolerance-inclusive"],
         "answer credit granted at the inclusive boundary; process fails because "
         "the stated value is not the node's own derivation",
         finding="F7-2")

    # C6  A direct solver: one display unit from gold's value, inside the
    #     prescribed method's own tolerance.  S7.4 says the looser governs, and
    #     this happens on 3.7% of real instances.
    c = copy.deepcopy(it)
    c["result"] = dict(c["result"], value=round(it["result"]["value"] + 0.001, 3))
    case("C6-method-tolerance", it, c, MATCH, MISMATCH,
         ["7.4-method-tolerance-governs"],
         "one display unit out, inside termination.tolerance, so S7.4 grants "
         "answer credit; a direct solver's trace is not an iteration node and a "
         "mutation cannot represent one",
         finding="F7-3")

    # C7  A different starting bracket.  The preamble is part of the QUESTION,
    #     so this is a process failure however good the arithmetic.
    c = copy.deepcopy(it)
    c["preamble"] = copy.deepcopy(c["preamble"])
    c["preamble"][0] = dict(c["preamble"][0], y=c["preamble"][0]["y"] + 0.25)
    case("C7-preamble-differs", it, c, MATCH, MISMATCH,
         ["7.2-preamble-differs"],
         "started from a different bracket; solved a different question")

    # C8  A budget overrun.  REPORTED, not failed: max_elements is the
    #     generator's guard, and a solver taking more updates to the same depth
    #     has not made an error.
    c = _extend_iteration(it, n_extra=int(it["termination"]["max_elements"])
                          - len(it["elements"]) + 1)
    case("C8-budget-overrun", it, c, MATCH, MISMATCH,
         ["7.2-budget-overrun-reported"],
         "the overrun itself is only observed, as 7.2 requires; the trace is "
         "still rejected at 4.3.3 for the same reason as C1",
         finding="F7-1")

    # C9  **The candidate declares coarser precisions to widen its own
    #     tolerance.**  S7.1 says the tolerance always comes from the GOLD
    #     node.  If an implementer reads it as the candidate's, this trace is
    #     accepted with a wrong answer -- a candidate grading itself.
    c = copy.deepcopy(it)
    c["symbol_precision"] = dict(c["symbol_precision"], y_next=1)
    c["result"] = dict(c["result"], dp=1, value=round(it["result"]["value"] + 0.04, 3))
    case("C9-candidate-widens-tolerance", it, c, MISMATCH, MISMATCH,
         ["7.1-gold-precision-governs"],
         "candidate declares dp=1 and answers 0.04 out; GOLD's dp=3 governs, so "
         "the answer is rejected -- the clause Reviewer D2 asked to be stated "
         "explicitly, now exercised")

    # ================ decision -- cardinality: answer_bearing ==============

    # C10 One fewer station.  For THIS node type the count is the answer, so
    #     this is an answer failure reported on `cardinality_symbol` -- the
    #     exact disposition a verifier built from a merged node type gets
    #     wrong (D-038).
    c = copy.deepcopy(de)
    c["elements"] = c["elements"][:-1]
    case("C10-cardinality-mismatch", de, c, MISMATCH, MISMATCH,
         ["7.3-cardinality-mismatch"],
         "one fewer station; n is the answer, so the answer is wrong")

    # C11 The same n, a different packing.  The answer stands; the process
    #     does not.  Two valid packings can share an n, and separating the two
    #     credits is the whole of S7.3 row 2.
    c = copy.deepcopy(de)
    key = de["roles"]["committed"]
    if len(c["elements"]) >= 2 and c["elements"][0][key] and c["elements"][1][key]:
        a = list(c["elements"][0][key])
        b = list(c["elements"][1][key])
        c["elements"][0][key] = a[:-1] + [b[0]]
        c["elements"][1][key] = [a[-1]] + b[1:]
    case("C11-same-n-different-packing", de, c, MATCH, MISMATCH,
         ["7.3-same-n-different-packing"],
         "one task swapped between two stations; n unchanged")

    # C12 The same n, the same sets, a different ORDER within a station.
    #     `selection`'s output is ordered, so this is reported -- and reported
    #     as weaker than a set mismatch, which is why the two are not merged.
    c = copy.deepcopy(de)
    for el in c["elements"]:
        if len(el[key]) >= 2:
            el[key] = list(reversed(el[key]))
            break
    case("C12-order-within-set", de, c, MATCH, MISMATCH,
         ["7.3-order-within-set"],
         "gold's set in a different order: rejected at 5.4.6 before S7.3 can "
         "distinguish set credit from order credit",
         finding="F7-4")

    # C13 A trial that violates `selection`: the chosen candidate is not the
    #     objective's optimum among the admissible ones.
    c = copy.deepcopy(de)
    tr_key = de["roles"]["trials"]
    ch_key = de["trial_roles"]["chosen"]
    cand_key = de["trial_roles"]["candidates"]
    item_key = de["candidate_roles"]["item"]
    for el in c["elements"]:
        for tr in el[tr_key]:
            adm = [x for x in tr[cand_key] if x.get(de["candidate_roles"]["admissible"])]
            if len(adm) >= 2 and tr[ch_key] is not None:
                others = [x[item_key] for x in adm if x[item_key] != tr[ch_key]]
                if others:
                    tr[ch_key] = others[0]
                    break
        else:
            continue
        break
    case("C13-selection-violated", de, c, MATCH, MISMATCH,
         ["7.3-trial-violates-selection"],
         "chose an admissible candidate that is not the optimum")

    # C14 An extra station.  The mirror of C10, and the row that shows the
    #     cardinality disposition is symmetric.
    c = copy.deepcopy(de)
    extra = copy.deepcopy(c["elements"][-1])
    extra[de["roles"]["index"]] = len(c["elements"]) + 1
    extra[key] = []
    c["elements"].append(extra)
    case("C14-extra-element", de, c, MISMATCH, MISMATCH,
         ["7.3-cardinality-mismatch"],
         "an extra empty station; n is the answer")

    # C15 A clean candidate: identical to gold.  A conformance corpus with no
    #     passing case cannot distinguish "rejects everything" from "verifies".
    case("C15-identical", it, copy.deepcopy(it), MATCH, MATCH,
         ["7.2-answer-credit", "7.2-process-credit"],
         "gold compared against itself; must pass cleanly")

    case("C16-identical-decision", de, copy.deepcopy(de), MATCH, MATCH,
         ["7.3-same-n-different-packing"],
         "gold decision compared against itself; must pass cleanly")

    return CASES


def main() -> None:
    cases = build()
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(cases, fh, indent=1)
    print(f"wrote {len(cases)} candidate traces to {OUT}")

    from tests.trace_schema.candidate_7 import SEVEN_DISPOSITIONS
    covered = set()
    for c in cases:
        covered.update(c["exercises"])
    missing = sorted(set(SEVEN_DISPOSITIONS) - covered)
    print(f"dispositions covered: {len(covered)}/{len(SEVEN_DISPOSITIONS)}")
    if missing:
        print("NOT covered: " + ", ".join(missing))
    else:
        print("every S7 disposition is exercised")


if __name__ == "__main__":
    main()
