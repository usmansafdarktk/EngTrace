"""Reviewer D2 -- negative cases for section 8, plus the reusability probe.

Every clause of section 8 ("what a conforming verifier must REJECT") is
exercised by mutating a gold node and asserting that
`reviewer_d2_verifier.verify_node` reports the expected clause.

The reusability probe builds, from the spec alone, an `iteration` node from a
DIFFERENT (synthetic) template whose iterate is called `S`, whose residual is
called `r`, whose index is called `n`, and whose `update_relation` is a
*damped* secant, and checks that the same unmodified verifier accepts it.

Usage:
    python -m tests.trace_schema.reviewer_d2_negatives \\
        docs/re-implementation-sep/phase3_conformance/traces.json
"""

from __future__ import annotations

import copy
import json
import sys
from decimal import Decimal, ROUND_HALF_UP

from tests.trace_schema.reviewer_d2_verifier import verify_node

RESULTS = []


def r4(x):
    return float(Decimal(repr(x)).quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP))


def r3(x):
    return float(Decimal(repr(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def check(name, clause, node, expect_clause, lenient=False):
    rep = verify_node(node, name, strict_roles=not lenient)
    clauses = [c for c, _ in rep.failures]
    hit = any(expect_clause in c for c in clauses)
    RESULTS.append((clause, name, "REJECTED" if hit else
                    ("MISSED (accepted)" if rep.passed else
                     "rejected for the WRONG reason: %s" % clauses),
                    hit))
    return hit


# --------------------------------------------------------------------------
# The synthetic alien-symbol iteration node -- built from the spec, not from
# any template source.  f(S) = S^2 - 2; damped secant with a 0.5 literal.
# --------------------------------------------------------------------------

def build_alien_node():
    def f(s):
        return r3(s * s - 2.0)

    s_prev, s_curr = 1.0, 2.0
    r_prev, r_curr = f(s_prev), f(s_curr)
    tol, elements, k = 0.002, [], 1
    while True:
        num = r_curr * (s_curr - s_prev)
        den = r_curr - r_prev
        s_next = r4(s_curr - 0.5 * num / den)
        delta = r4(abs(s_next - s_curr))
        done = delta < tol
        frame = None if done else {"S": s_next, "r": f(s_next)}
        elements.append({"n": k, "S_prev": s_prev, "r_prev": r_prev,
                         "S_curr": s_curr, "r_curr": r_curr, "S_next": s_next,
                         "delta": delta, "done": done, "frame": frame})
        if done:
            break
        s_prev, r_prev = s_curr, r_curr
        s_curr, r_curr = s_next, frame["r"]
        k += 1
        if k > 60:
            raise SystemExit("alien node failed to converge")
    return {
        "schema_version": "1.1",
        "node_type": "iteration",
        "node_id": "tXX_damped_secant_sqrt2",
        "cardinality": "incidental",
        "roles": {"index": "n", "iterate_prev": "S_prev", "iterate_curr": "S_curr",
                  "iterate_next": "S_next", "residual_prev": "r_prev",
                  "residual_curr": "r_curr", "change": "delta",
                  "converged": "done", "evaluation": "frame"},
        "update_relation": "iterate_curr - 0.5 * residual_curr * "
                           "(iterate_curr - iterate_prev) / "
                           "(residual_curr - residual_prev)",
        "termination": {"kind": "convergence", "quantity_role": "change",
                        "comparison": "lt", "tolerance": tol,
                        "predicate_prose": "abs(S_next - S_curr) < 0.002",
                        "max_elements": 60, "satisfied": True},
        "rounding": "decimal-half-up",
        "symbol_precision": {"S_prev": 4, "S_curr": 4, "S_next": 4, "delta": 4,
                             "r_prev": 3, "r_curr": 3, "S": 4, "r": 3},
        "preamble_symbols": ["S", "r"],
        "evaluation_symbols": ["S", "r"],
        "evaluation_roles": {"iterate": "S", "residual": "r"},
        "element_symbols": ["n", "S_prev", "r_prev", "S_curr", "r_curr",
                            "S_next", "delta", "done", "frame"],
        "carry": {"S_prev": "S_curr", "r_prev": "r_curr",
                  "S_curr": "S_next", "r_curr": "frame.r"},
        "carry_notes": {},
        "preamble": [{"S": 1.0, "r": -1.0}, {"S": 2.0, "r": 2.0}],
        "elements": elements,
        "result": {"symbol": "Sroot", "value": r3(elements[-1]["S_next"]),
                   "dp": 3, "from": "last.iterate_next"},
    }


def main(argv):
    rows = json.load(open(argv[1], encoding="utf-8"))
    it = next(r["trace_nodes"] for r in rows
              if r["trace_nodes"]["node_type"] == "iteration")
    de = next(r["trace_nodes"] for r in rows
              if r["trace_nodes"]["node_type"] == "decision")

    def I():
        return copy.deepcopy(it)

    def D():
        return copy.deepcopy(de)

    # --- 8.1 -------------------------------------------------------------
    n = I(); n["termination"]["satisfied"] = False
    check("8.1 iteration satisfied:false", "8.1", n, "8.1")
    n = D(); n["termination"]["satisfied"] = False
    check("8.1 decision satisfied:false", "8.1", n, "8.1", lenient=True)

    # --- 8.2 (all five declared key sets) --------------------------------
    n = I(); n["elements"][0]["bogus"] = 1
    check("8.2 element extra key", "8.2", n, "8.2")
    n = I(); del n["elements"][0]["change"]
    check("8.2 element missing key", "8.2", n, "8.2")
    n = I(); del n["elements"][0]["evaluation"]["A"]
    check("8.2 evaluation frame missing key", "8.2", n, "8.2")
    n = I(); n["preamble"][0]["bogus"] = 1
    check("8.2 preamble frame extra key", "8.2", n, "8.2")
    n = D(); n["elements"][0]["trials"][0]["bogus"] = 1
    check("8.2 trial extra key", "8.2", n, "8.2", lenient=True)
    n = D(); del n["elements"][0]["trials"][0]["eligible"][0]["fits"]
    check("8.2 eligible entry missing key", "8.2", n, "8.2", lenient=True)

    # --- 8.3 -------------------------------------------------------------
    n = I(); n["elements"][1]["y_prev"] = 9.9999
    check("8.3 carry relation broken", "8.3", n, "8.3")
    n = I(); n["carry"]["y_prev"] = "y_curr[0]"
    check("8.3 carry path not well-formed (index)", "8.3", n, "8.3")
    n = I(); n["carry"]["y_prev"] = "evaluation .g"
    check("8.3 carry path not well-formed (space)", "8.3", n, "8.3")
    n = I(); n["carry"]["y_prev"] = "nosuch.g"
    check("8.3 carry path root not an element symbol", "8.3", n, "8.3")

    # --- 8.4 (a verifier behaviour, not a trace mutation) ----------------
    n = I(); n["carry_notes"] = {"y_curr": "refined each pass"}
    rep = verify_node(n, "8.4 carry_notes present")
    ok = rep.passed and rep.status == "PASS+UNCHECKED" and len(rep.unchecked) == 1
    RESULTS.append(("8.4", "carry_notes surfaced on the unchecked channel",
                    "REJECTED-as-unchecked (status %s)" % rep.status
                    if ok else "MISSED (status %s)" % rep.status, ok))

    # --- 8.5 -------------------------------------------------------------
    n = I(); n["elements"][0]["converged"] = True
    check("8.5 converged true before the last element", "8.5", n, "8.5")
    n = I()
    n["elements"][-1]["evaluation"] = dict(n["elements"][0]["evaluation"])
    check("8.5 evaluation non-null on the terminating element", "8.5", n, "8.5")
    n = I(); n["elements"][0]["evaluation"]["y"] = 42.0
    check("8.5 evaluation frame bound to the wrong iterate", "8.5", n, "8.5")
    n = I(); n["elements"][0]["g_prev"] = n["elements"][0]["g_curr"]
    check("8.5 zero update denominator", "8.5", n, "4.3.7")
    n = I(); n["elements"][-1]["k"] = 99
    check("8.5 non-contiguous index", "8.5", n, "4.3.8")

    # --- 8.6 -------------------------------------------------------------
    n = D(); n["elements"][1]["assigned"] = list(n["elements"][0]["assigned"])
    check("8.6 overlapping assigned sets", "8.6", n, "5.3", lenient=True)
    n = D(); n["termination"]["universe"] = n["termination"]["universe"] + ["z"]
    n["termination"]["precedences"]["z"] = []
    check("8.6 union != universe", "8.6", n, "5.3", lenient=True)
    n = D()
    for el in n["elements"]:
        for tr in el["trials"]:
            if len(tr["eligible"]) > 1 and tr["chosen"] is not None:
                cands = sorted((e for e in tr["eligible"] if e["fits"]),
                               key=lambda e: e["duration"])
                if len(cands) > 1 and cands[0]["task"] != tr["chosen"]:
                    tr["chosen"] = cands[0]["task"]
                    break
    check("8.6 chosen does not optimise selection", "8.6", n, "5.4.4", lenient=True)
    n = D()
    for el in n["elements"]:
        for tr in el["trials"]:
            if len(tr["eligible"]) > 1:
                del tr["eligible"][1]
                break
        else:
            continue
        break
    check("8.6 eligible omits a precedence-satisfied item", "8.6", n,
          "5.4.3b", lenient=True)
    n = D()
    last = n["termination"]["universe"][-1]
    tr0 = n["elements"][0]["trials"][0]
    dur = {e["task"]: e["duration"] for el in n["elements"]
           for t in el["trials"] for e in t["eligible"]}
    tr0["eligible"].append({"task": last, "duration": dur[last], "fits": False})
    check("8.6 eligible includes an item whose precedences are unmet", "8.6", n,
          "5.4.3b", lenient=True)
    n = D(); n["elements"][-1]["station_id"] = 99
    check("8.6 non-contiguous station_id", "8.6", n, "5.3", lenient=True)
    n = D(); n["elements"][-1]["capacity"] = n["elements"][-1]["capacity"] + 1
    check("8.6 varying capacity under capacity_constant", "8.6", n,
          "5.4.8", lenient=True)

    # --- 8.7 -------------------------------------------------------------
    n = D(); del n["cardinality_symbol"]
    check("8.7 answer_bearing with no cardinality_symbol", "8.7", n, "8.7",
          lenient=True)
    n = I(); n["cardinality_symbol"] = "k"
    check("8.7 incidental carrying a cardinality_symbol", "8.7", n, "8.7")

    # --- 8.8 -------------------------------------------------------------
    n = I(); del n["roles"]["change"]
    check("8.8 mandatory role missing from roles", "8.8", n, "8.8")
    n = I(); n["roles"]["change"] = "not_a_symbol"
    check("8.8 roles value not in element_symbols", "8.8", n, "8.8")
    n = I(); n["evaluation_roles"]["residual"] = "not_a_symbol"
    check("8.8 evaluation_roles value not in evaluation_symbols", "8.8", n, "8.8")

    # --- 8.9 -------------------------------------------------------------
    n = I(); n["schema_version"] = "2.0"
    check("8.9 unknown major schema version", "8.9", n, "8.9")

    # --- extra: section 6.4 lower bound ----------------------------------
    n = I(); n["elements"] = []
    check("6.4 zero-element sequence", "6.4", n, "6.4")

    # --- report ----------------------------------------------------------
    print("=" * 78)
    print("SECTION 8 NEGATIVE CASES")
    print("=" * 78)
    for clause, name, outcome, ok in RESULTS:
        print("%-5s %-4s %-58s %s" % ("ok" if ok else "MISS", clause, name, outcome))
    missed = [r for r in RESULTS if not r[3]]
    print("\n%d negative cases, %d rejected as expected, %d missed"
          % (len(RESULTS), len(RESULTS) - len(missed), len(missed)))

    # --- the reusability probe -------------------------------------------
    print()
    print("=" * 78)
    print("REUSABILITY PROBE -- an `iteration` node from a different template")
    print("=" * 78)
    alien = build_alien_node()
    print("iterate symbol: %s   residual symbol: %s   index symbol: %s"
          % (alien["roles"]["iterate_curr"], alien["roles"]["residual_curr"],
             alien["roles"]["index"]))
    print("update_relation: %s" % alien["update_relation"])
    print("elements: %d   result: %r" % (len(alien["elements"]), alien["result"]))
    rep = verify_node(alien, "alien iteration node")
    print("verdict: %s" % rep.status)
    for c, m in rep.failures:
        print("   [%s] %s" % (c, m))
    # and the same node must still be rejected when corrupted
    bad = copy.deepcopy(alien)
    bad["elements"][0]["S_next"] = bad["elements"][0]["S_next"] + 0.01
    rep2 = verify_node(bad, "alien node, corrupted update")
    print("corrupted alien node verdict: %s  %s"
          % (rep2.status, [c for c, _ in rep2.failures]))

    return 1 if missed or not rep.passed or rep2.passed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
