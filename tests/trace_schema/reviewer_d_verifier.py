"""Reviewer D's independent verifier for the D3.3 `iteration`/`decision` node types.

Written from `docs/re-implementation-sep/track-a/phase3_node_types.md` alone, as the
Phase 3 exit gate requires ("Reviewer D implemented a working verifier from D3.3
alone").  The template sources and `tests/trace_schema/extract.py` were NOT read.

Every place the spec did not determine the implementation is marked with an
`ASSUMPTION (D-An)` comment and cross-referenced from
`docs/re-implementation-sep/reviews/phase3_reviewer_d_schema.md`.

Usage:
    python -m tests.trace_schema.reviewer_d_verifier <traces.json>
    python -m tests.trace_schema.reviewer_d_verifier <traces.json> --self-test
"""

from __future__ import annotations

import copy
import json
import re
import sys

# --------------------------------------------------------------------------
# Assumptions the spec forced.  Each is a finding in the review report.
# --------------------------------------------------------------------------
#
# D-A1  Container shape.  Section 9 says traces.json "carries 40 seeds of each
#       template as {seed, question, solution, trace_nodes}".  It never says
#       whether `trace_nodes` (plural) is one node, a list of nodes, or a
#       node_id->node map.  Observed: a single node object.  Accept all three.
#
# D-A2  `iteration` element schema.  Section 5 gives `decision` a normative
#       element and trial table.  Section 4 gives `iteration` no such table: its
#       invariants name y_prev/y_curr/y_next/g_prev/g_curr/change/converged/
#       evaluation as if they were type-level names, but section 3 says
#       `element_symbols` is the template's own domain.  Assumed: those eight
#       names plus `k` are type-level for `iteration`.  A template that called
#       its iterate `h` would not verify.
#
# D-A3  Update relation.  Section 4.1 states the relation only "for the
#       reference template" (the secant step).  Nothing in the node encodes it,
#       so a verifier cannot be generic.  Hardcoded below.
#
# D-A4  Per-symbol display precision.  Sections 4.2 / 7.1 speak of "the display
#       precision of y" and "the element symbol's display precision", but the
#       node declares `dp` only on `result`.  Inferred per symbol from the data.
#
# D-A5  Path vs prose in `carry`.  Section 3.2 says a RHS "may also be a prose
#       description"; no syntax distinguishes them.  Heuristic below.
#
# D-A6  Predicate evaluation.  `termination.predicate` is "in the template's own
#       notation" with free names (`tol`, `tasks`).  No grammar or binding
#       environment is specified.  Bound by hand.
#
# D-A7  Precedence.  Section 5 global invariants and section 8.6 require
#       rejecting "a precedence violated at the moment of choice", but the node
#       carries no precedence relation at all.  Reported UNCHECKABLE, never as a
#       pass.

ITERATION_ELEMENT_SYMBOLS = [
    "k", "y_prev", "g_prev", "y_curr", "g_curr",
    "y_next", "change", "converged", "evaluation",
]                                                              # D-A2
DECISION_ELEMENT_SYMBOLS = [
    "station_id", "capacity", "remaining_initial",
    "trials", "assigned", "remaining_final",
]                                                              # section 5 table
DECISION_TRIAL_SYMBOLS = ["eligible", "chosen", "remaining_after", "closes"]
ELIGIBLE_KEYS = {"task", "duration", "fits"}                   # section 5, prose only

# D-A4: display precision per element symbol, inferred not declared.
DISPLAY_DP = {
    "y_prev": 4, "y_curr": 4, "y_next": 4, "change": 4,
    "g_prev": 3, "g_curr": 3,
    "y": 4, "A": 3, "P": 3, "AR": 3, "g": 3,
}
DEFAULT_DP = 4

_PATH_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)*$")


class Report:
    """Collects failures (reject) and unchecked relations (section 8.4)."""

    def __init__(self, label):
        self.label = label
        self.failures = []
        self.unchecked = []

    def fail(self, code, msg):
        self.failures.append((code, msg))

    def uncheckable(self, code, msg):
        self.unchecked.append((code, msg))

    @property
    def ok(self):
        return not self.failures


def _tol_for(dp):
    """Section 7.1: half a unit in the last place gold displays it."""
    return 0.5 * 10 ** (-dp) + 1e-12


def _close(a, b, dp):
    return abs(a - b) <= _tol_for(dp)


def _is_path(rhs):
    """D-A5: treat a RHS as a checkable path iff it is a dotted identifier."""
    return bool(_PATH_RE.match(rhs.strip()))


def _resolve(obj, path):
    cur = obj
    for seg in path.split("."):
        if not isinstance(cur, dict) or seg not in cur:
            return (False, None)
        cur = cur[seg]
    return (True, cur)


# --------------------------------------------------------------------------
# Section 6.1 Shape / 6.2 Homogeneity / 6.3 Budget  (common)
# --------------------------------------------------------------------------

REQUIRED_BASE = ["node_type", "node_id", "cardinality", "termination",
                 "element_symbols", "elements", "carry", "result"]
REQUIRED_TERM = ["kind", "quantity", "predicate", "satisfied", "max_elements"]

FIXED_BY_TYPE = {
    "iteration": ("incidental", "convergence"),      # section 4
    "decision": ("answer_bearing", "exhaustion"),    # section 5
}


def check_shape(node, rep):
    for f in REQUIRED_BASE:
        if f not in node:
            rep.fail("SHAPE", "missing required field %r (s3)" % f)
    nt = node.get("node_type")
    if nt not in FIXED_BY_TYPE:
        rep.fail("SHAPE", "unrecognised node_type %r (s3)" % (nt,))
        return False

    want_card, want_kind = FIXED_BY_TYPE[nt]
    if node.get("cardinality") != want_card:
        rep.fail("SHAPE", "node_type %r fixes cardinality %r, got %r (s4/s5)"
                 % (nt, want_card, node.get("cardinality")))

    term = node.get("termination")
    if not isinstance(term, dict):
        rep.fail("SHAPE", "termination missing or not an object (s3.1)")
        return False
    for f in REQUIRED_TERM:
        if f not in term:
            rep.fail("SHAPE", "termination missing required field %r (s3.1)" % f)
    if term.get("kind") != want_kind:
        rep.fail("SHAPE", "node_type %r fixes termination.kind %r, got %r (s4/s5)"
                 % (nt, want_kind, term.get("kind")))
    if term.get("kind") == "convergence" and "tolerance" not in term:
        rep.fail("SHAPE", "convergence requires termination.tolerance (s3.1)")
    if term.get("kind") == "exhaustion" and "universe" not in term:
        rep.fail("SHAPE", "exhaustion requires termination.universe (s3.1)")

    # section 8.7 / section 3: cardinality_symbol present iff answer_bearing
    has_cs = "cardinality_symbol" in node and node["cardinality_symbol"] is not None
    if node.get("cardinality") == "answer_bearing" and not has_cs:
        rep.fail("CARD_SYMBOL", "answer_bearing without cardinality_symbol (s8.7)")
    if node.get("cardinality") == "incidental" and has_cs:
        rep.fail("CARD_SYMBOL", "incidental with a cardinality_symbol (s8.7)")

    # section 3: node_id "must not encode the element count"
    nid = node.get("node_id")
    if isinstance(nid, str) and isinstance(node.get("elements"), list):
        if re.search(r"(?<!\d)%d(?!\d)" % len(node["elements"]), nid):
            rep.fail("NODE_ID",
                     "node_id %r appears to encode the element count (s3)" % nid)

    # section 3: conditionally-required companions
    if nt == "decision":
        for f in ("rule", "trial_symbols"):
            if f not in node:
                rep.fail("SHAPE", "decision requires %r (s3)" % f)
    if "preamble" in node and "preamble_symbols" not in node:
        rep.fail("SHAPE", "preamble present without preamble_symbols (s3)")

    # section 8.1: a gold node must not carry satisfied: false
    if term.get("satisfied") is not True:
        rep.fail("SATISFIED", "gold node with termination.satisfied=%r (s3.1, s8.1)"
                 % (term.get("satisfied"),))
    return True


def check_homogeneity(node, rep):
    """Section 6.2 / 8.2 -- key set equality in BOTH directions."""
    syms = node.get("element_symbols") or []
    for i, el in enumerate(node.get("elements") or []):
        if not isinstance(el, dict):
            rep.fail("HOMOGENEITY", "element %d is not an object" % i)
            continue
        if set(el.keys()) != set(syms):
            missing = sorted(set(syms) - set(el.keys()))
            extra = sorted(set(el.keys()) - set(syms))
            rep.fail("HOMOGENEITY",
                     "element %d key set != element_symbols (missing=%s extra=%s) (s6.2, s8.2)"
                     % (i + 1, missing, extra))
    if "preamble" in node:
        psyms = set(node.get("preamble_symbols") or [])
        for i, fr in enumerate(node["preamble"]):
            if set(fr.keys()) != psyms:
                rep.fail("HOMOGENEITY",
                         "preamble frame %d key set != preamble_symbols (s6.2)" % (i + 1))
    if node.get("node_type") == "decision":
        tsyms = set(node.get("trial_symbols") or [])
        for i, el in enumerate(node.get("elements") or []):
            for j, tr in enumerate(el.get("trials") or []):
                if set(tr.keys()) != tsyms:
                    rep.fail("HOMOGENEITY",
                             "element %d trial %d key set != trial_symbols (s6.2)"
                             % (i + 1, j + 1))
                for x in tr.get("eligible") or []:
                    if set(x.keys()) != ELIGIBLE_KEYS:
                        rep.fail("HOMOGENEITY",
                                 "element %d trial %d eligible entry keys %s != %s (s5)"
                                 % (i + 1, j + 1, sorted(x.keys()), sorted(ELIGIBLE_KEYS)))


def check_budget(node, rep):
    """Section 6.3 -- 1 <= len(elements) <= max_elements."""
    els = node.get("elements")
    if not isinstance(els, list):
        rep.fail("BUDGET", "elements is not a list")
        return
    term = node.get("termination") or {}
    if len(els) < 1:
        # section 3 allows empty only if satisfied at entry; gold never is.
        rep.fail("BUDGET", "elements is empty (s6.3)")
    mx = term.get("max_elements")
    if isinstance(mx, int) and len(els) > mx:
        rep.fail("BUDGET",
                 "len(elements)=%d > max_elements=%d (s6.3, invariant not hope)"
                 % (len(els), mx))


def check_carry(node, rep):
    """Section 6.4 / 8.3 / 8.4."""
    carry = node.get("carry")
    if not isinstance(carry, dict):
        rep.fail("CARRY", "carry missing or not an object (s3.2)")
        return
    els = node.get("elements") or []
    for lhs, rhs in carry.items():
        if not isinstance(rhs, str) or not _is_path(rhs):
            # section 8.4: must be REPORTED as unchecked, never silently passed.
            rep.uncheckable("CARRY_UNCHECKED",
                            "carry[%r] = %r is not a path; advisory and UNCHECKED (s3.2, s8.4)"
                            % (lhs, rhs))
            continue
        for k in range(len(els) - 1):
            if lhs not in els[k + 1]:
                rep.fail("CARRY", "carry LHS %r absent from element %d (s6.4)"
                         % (lhs, k + 2))
                continue
            okr, want = _resolve(els[k], rhs)
            if not okr:
                rep.fail("CARRY", "carry path %r does not resolve in element %d (s6.4)"
                         % (rhs, k + 1))
                continue
            got = els[k + 1][lhs]
            same = (got == want)
            if (not same and isinstance(got, (int, float))
                    and isinstance(want, (int, float))
                    and not isinstance(got, bool) and not isinstance(want, bool)):
                same = _close(float(got), float(want), DISPLAY_DP.get(lhs, DEFAULT_DP))
            if not same:
                rep.fail("CARRY",
                         "carry %r<-%r broken between elements %d/%d: %r != %r (s6.4, s8.3)"
                         % (lhs, rhs, k + 1, k + 2, got, want))


# --------------------------------------------------------------------------
# Section 4 -- iteration local invariants
# --------------------------------------------------------------------------

def check_iteration(node, rep):
    els = node["elements"]
    term = node["termination"]
    tol = term.get("tolerance")
    last = len(els) - 1
    for i, e in enumerate(els):
        tag = "element %d" % (i + 1)
        bad = False
        for s in ("y_prev", "g_prev", "y_curr", "g_curr", "y_next", "change"):
            if not isinstance(e.get(s), (int, float)) or isinstance(e.get(s), bool):
                rep.fail("ITER_SHAPE", "%s: %r is not numeric (s4, D-A2)" % (tag, s))
                bad = True
        if bad:
            continue

        # section 4.6 / 8.5 -- non-zero update denominator.
        den = e["g_curr"] - e["g_prev"]
        if den == 0:
            rep.fail("ITER_DENOM",
                     "%s: zero update denominator g_curr-g_prev (s4.6, s8.5)" % tag)
        else:
            # section 4.1 -- the update relation holds, to section 7.1 tolerance.  D-A3.
            pred = e["y_curr"] - e["g_curr"] * (e["y_curr"] - e["y_prev"]) / den
            if not _close(pred, e["y_next"], DISPLAY_DP["y_next"]):
                rep.fail("ITER_UPDATE",
                         "%s: secant step gives %.6f, trace says %r (s4.1)"
                         % (tag, pred, e["y_next"]))

        # section 4.2 -- change == |y_next - y_curr| at display precision of y.
        want = abs(e["y_next"] - e["y_curr"])
        if not _close(want, e["change"], DISPLAY_DP["change"]):
            rep.fail("ITER_CHANGE", "%s: change=%r but |y_next-y_curr|=%.6f (s4.2)"
                     % (tag, e["change"], want))

        # section 4.3 -- converged == (change < tolerance)
        conv = e.get("converged")
        if not isinstance(conv, bool):
            rep.fail("ITER_CONV", "%s: converged is not boolean (s4.3)" % tag)
        else:
            if conv != (e["change"] < tol):
                rep.fail("ITER_CONV",
                         "%s: converged=%r but change=%r vs tolerance=%r (s4.3)"
                         % (tag, conv, e["change"], tol))
            # section 4.4 / 8.5 -- true on the last element only.
            if conv and i != last:
                rep.fail("ITER_CONV",
                         "%s: converged=true on a non-final element (s4.4, s8.5)" % tag)
            if not conv and i == last:
                rep.fail("ITER_CONV",
                         "%s: converged=false on the final element (s4.4)" % tag)

        # section 4.5 / 8.5 -- evaluation is null iff converged.
        ev_null = e.get("evaluation") is None
        if isinstance(conv, bool) and ev_null != conv:
            rep.fail("ITER_EVAL",
                     "%s: evaluation %s but converged=%r (s4.5, s8.5)"
                     % (tag, "is null" if ev_null else "is non-null", conv))

    # section 4 preamble -- prescribed by the question, not free.
    if "preamble" in node:
        pre = node["preamble"]
        if len(pre) >= 2 and els and isinstance(els[0].get("y_prev"), (int, float)):
            for idx, (pk, ek) in enumerate([("y", "y_prev"), ("y", "y_curr")]):
                pass
            checks = [(0, "y", "y_prev"), (1, "y", "y_curr"),
                      (0, "g", "g_prev"), (1, "g", "g_curr")]
            for j, pk, ek in checks:
                if pk not in pre[j] or ek not in els[0]:
                    continue
                if not _close(float(pre[j][pk]), float(els[0][ek]),
                              DISPLAY_DP.get(ek, DEFAULT_DP)):
                    rep.fail("ITER_PREAMBLE",
                             "preamble[%d].%s=%r does not seed element 1 %s=%r "
                             "(s4 Preamble)" % (j, pk, pre[j][pk], ek, els[0][ek]))


# --------------------------------------------------------------------------
# Section 5 -- decision local + global invariants
# --------------------------------------------------------------------------

def check_decision(node, rep):
    els = node["elements"]
    term = node["termination"]
    universe = list(term.get("universe") or [])
    last = len(els) - 1

    for i, e in enumerate(els):
        tag = "station %s" % e.get("station_id")

        # section 5.1
        if e.get("remaining_initial") != e.get("capacity"):
            rep.fail("DEC_OPEN", "%s: remaining_initial=%r != capacity=%r (s5.1)"
                     % (tag, e.get("remaining_initial"), e.get("capacity")))

        trials = e.get("trials")
        if not isinstance(trials, list) or not trials:
            rep.fail("DEC_TRIALS", "%s: trials missing or empty (s5 'Never empty')" % tag)
            continue

        budget = e.get("remaining_initial")
        chosen_seq = []
        for j, tr in enumerate(trials):
            elig = tr.get("eligible")
            if not isinstance(elig, list):
                rep.fail("DEC_TRIALS", "%s trial %d: eligible is not a list (s5)"
                         % (tag, j + 1))
                continue
            by_task = {}
            for x in elig:
                if not isinstance(x, dict) or "task" not in x:
                    continue
                by_task[x["task"]] = x
                # `fits` is reported separately so the verifier checks it (s5).
                want_fits = x["duration"] <= budget
                if bool(x["fits"]) != want_fits:
                    rep.fail("DEC_FITS",
                             "%s trial %d: task %r fits=%r but duration %r vs budget %r (s5)"
                             % (tag, j + 1, x["task"], x["fits"], x["duration"], budget))

            ch = tr.get("chosen")
            fitting = [x for x in elig if x.get("fits")]
            if ch is None:
                # section 5.3 contrapositive: null only when nothing fits.
                if fitting:
                    rep.fail("DEC_RULE",
                             "%s trial %d: chosen=null but %d candidate(s) fit (s5.3, s8.6)"
                             % (tag, j + 1, len(fitting)))
                # section 5.4
                if tr.get("remaining_after") != budget:
                    rep.fail("DEC_BUDGET",
                             "%s trial %d: chosen=null must leave budget %r, got %r (s5.4)"
                             % (tag, j + 1, budget, tr.get("remaining_after")))
                # section 5 Trial: closes true exactly where chosen is null
                if tr.get("closes") is not True:
                    rep.fail("DEC_CLOSES",
                             "%s trial %d: chosen=null must have closes=true (s5 Trial)"
                             % (tag, j + 1))
            else:
                # section 5.3 -- c in eligible, fits, and is the longest fitting.
                if ch not in by_task:
                    rep.fail("DEC_RULE", "%s trial %d: chosen %r not in eligible (s5.3, s8.6)"
                             % (tag, j + 1, ch))
                    budget = tr.get("remaining_after")
                    continue
                if not by_task[ch].get("fits"):
                    rep.fail("DEC_RULE", "%s trial %d: chosen %r does not fit (s5.3, s8.6)"
                             % (tag, j + 1, ch))
                if fitting:
                    best = max(x["duration"] for x in fitting)
                    if by_task[ch]["duration"] != best:
                        rep.fail("DEC_RULE",
                                 "%s trial %d: chose %r (%r) but longest fitting is %r "
                                 "(s5.3, s8.6)"
                                 % (tag, j + 1, ch, by_task[ch]["duration"], best))
                    else:
                        tied = [x["task"] for x in fitting if x["duration"] == best]
                        if len(tied) > 1 and tied[0] != ch:
                            # tie-break by scan order (`rule`); eligible is the
                            # deterministic scan order per section 5 Determinism.
                            rep.fail("DEC_TIE",
                                     "%s trial %d: tie %s broken to %r, scan order gives %r "
                                     "(s5.3)" % (tag, j + 1, tied, ch, tied[0]))
                if tr.get("closes") is not False:
                    rep.fail("DEC_CLOSES",
                             "%s trial %d: a committing trial must have closes=false (s5.5)"
                             % (tag, j + 1))
                # section 5.4
                want = budget - by_task[ch]["duration"]
                if tr.get("remaining_after") != want:
                    rep.fail("DEC_BUDGET",
                             "%s trial %d: remaining_after=%r, expected %r (s5.4)"
                             % (tag, j + 1, tr.get("remaining_after"), want))
                chosen_seq.append(ch)

            # section 5.2 -- trials are consecutive.
            budget = tr.get("remaining_after")

        # section 5.5 -- closes true on the last trial only.
        for j, tr in enumerate(trials[:-1]):
            if tr.get("closes"):
                rep.fail("DEC_CLOSES", "%s trial %d: closes=true before the last trial (s5.5)"
                         % (tag, j + 1))
        if trials[-1].get("closes") is not True:
            rep.fail("DEC_CLOSES", "%s: last trial does not close (s5.5)" % tag)
        if not trials[-1].get("eligible") and i != last:
            rep.fail("DEC_CLOSES",
                     "%s: empty eligible on a closing trial is permitted only on the "
                     "final element (s5.5)" % tag)

        # section 5.6
        if list(e.get("assigned") or []) != chosen_seq:
            rep.fail("DEC_ASSIGNED",
                     "%s: assigned=%r != ordered chosen values %r (s5.6)"
                     % (tag, e.get("assigned"), chosen_seq))

        # section 5.7
        durs = {}
        for tr in trials:
            for x in tr.get("eligible") or []:
                if isinstance(x, dict) and "task" in x:
                    durs[x["task"]] = x["duration"]
        tot = sum(durs.get(t, 0) for t in (e.get("assigned") or []))
        if e.get("remaining_final") != e.get("capacity") - tot:
            rep.fail("DEC_FINAL",
                     "%s: remaining_final=%r != capacity-sum(durations)=%r (s5.7)"
                     % (tag, e.get("remaining_final"), e.get("capacity") - tot))

    # Global invariants (section 5).
    ids = [e.get("station_id") for e in els]
    if ids != list(range(1, len(els) + 1)):
        rep.fail("DEC_IDS", "station_id %r is not 1..%d contiguous ascending (s5, s8.6)"
                 % (ids, len(els)))
    seen, dupes = set(), []
    for e in els:
        for t in e.get("assigned") or []:
            if t in seen:
                dupes.append(t)
            seen.add(t)
    if dupes:
        rep.fail("DEC_DISJOINT", "assigned sets overlap on %r (s5, s8.6)"
                 % sorted(set(dupes)))
    if seen != set(universe):
        rep.fail("DEC_UNIVERSE",
                 "union of assigned %r != termination.universe %r (s5, s8.6)"
                 % (sorted(seen), sorted(universe)))

    # section 8.6 requires rejecting a precedence violation, but the node carries
    # no precedence relation.  Never silently pass it (s8.4's principle).  D-A7.
    rep.uncheckable("PRECEDENCE_UNCHECKABLE",
                    "s5 global invariant / s8.6 'a precedence violated at the moment of "
                    "choice' is UNCHECKABLE: the node exposes no precedence relation, "
                    "only prose in `rule`")


# --------------------------------------------------------------------------
# Section 6.6 Termination / 6.7 Result
# --------------------------------------------------------------------------

def check_termination(node, rep):
    term = node["termination"]
    els = node["elements"]
    if not els:
        return
    if term.get("kind") == "convergence":
        # D-A6: predicate string is not machine-evaluable; bind by hand via
        # `change` and `tolerance`, which section 4.3 defines the predicate to be.
        tol = term.get("tolerance")
        for i, e in enumerate(els):
            v = e.get("change")
            holds = isinstance(v, (int, float)) and not isinstance(v, bool) and v < tol
            if holds and i != len(els) - 1:
                rep.fail("TERM", "predicate holds on non-final element %d (s6.6)" % (i + 1))
            if not holds and i == len(els) - 1:
                rep.fail("TERM", "predicate does not hold on the last element (s6.6)")
        lastv = els[-1].get("change")
        agrees = isinstance(lastv, (int, float)) and lastv < tol
        if term.get("satisfied") is not agrees:
            rep.fail("TERM", "satisfied=%r disagrees with the predicate (s6.6)"
                     % term.get("satisfied"))
    else:
        # "len(assigned) == len(tasks)" with `tasks` bound to universe (D-A6).
        got = sum(len(e.get("assigned") or []) for e in els)
        want = len(term.get("universe") or [])
        if term.get("satisfied") is not (got == want):
            rep.fail("TERM", "satisfied=%r but %d/%d of the universe assigned (s6.6)"
                     % (term.get("satisfied"), got, want))
        for i in range(len(els) - 1):
            run = sum(len(x.get("assigned") or []) for x in els[:i + 1])
            if run == want:
                rep.fail("TERM", "universe exhausted at element %d, before the last (s6.6)"
                         % (i + 1))


def check_result(node, rep):
    res = node.get("result")
    if not isinstance(res, dict) or "symbol" not in res or "value" not in res:
        rep.fail("RESULT", "result must carry at least {symbol, value} (s3)")
        return
    els = node["elements"]
    if not els:
        return
    if node["node_type"] == "iteration":
        dp = res.get("dp")
        if dp is None:
            rep.uncheckable("RESULT_DP", "result.dp absent; s6.7 needs it to round (s3)")
            return
        yn = els[-1].get("y_next")
        if not isinstance(yn, (int, float)) or isinstance(yn, bool):
            rep.fail("RESULT", "last element has no numeric y_next to round (s6.7)")
            return
        want = round(yn, dp)
        if not _close(float(res["value"]), want, dp + 1):
            rep.fail("RESULT", "result.value=%r != round(last.y_next, %d)=%r (s6.7)"
                     % (res["value"], dp, want))
    else:
        if res["value"] != len(els):
            rep.fail("RESULT", "result.value=%r != len(elements)=%d (s6.7)"
                     % (res["value"], len(els)))
        cs = node.get("cardinality_symbol")
        if cs is not None and res["symbol"] != cs:
            rep.fail("RESULT", "result.symbol=%r != cardinality_symbol=%r (s3, s7.3)"
                     % (res["symbol"], cs))


# --------------------------------------------------------------------------
# Section 6.8 Prose agreement (partial -- see report finding)
# --------------------------------------------------------------------------

_NUM_RE = re.compile(r"-?\d+(?:\.\d+)?")


def check_prose(node, solution, rep):
    """Section 6.8: every numeric token the node carries appears in the prose at
    its stated precision.  Only the unambiguous half is implemented -- the
    result value.  See report finding F10 for why the general form is not."""
    if not solution:
        return
    toks = set(_NUM_RE.findall(solution))
    res = node.get("result") or {}
    v = res.get("value")
    if v is None:
        return
    dp = res.get("dp")
    cands = {str(v)}
    if isinstance(v, float) and dp is not None:
        cands.add(("%%.%df" % dp) % v)
    if isinstance(v, (int, float)) and float(v) == int(v):
        cands.add(str(int(v)))
    if not (cands & toks):
        rep.fail("PROSE",
                 "result.value %r does not appear in the printed solution (s6.8)" % v)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def verify_node(node, label, solution=None):
    rep = Report(label)
    if not isinstance(node, dict):
        rep.fail("SHAPE", "node is not an object")
        return rep
    if not check_shape(node, rep):
        return rep
    check_homogeneity(node, rep)
    check_budget(node, rep)
    check_carry(node, rep)
    if isinstance(node.get("elements"), list) and node["elements"]:
        try:
            if node["node_type"] == "iteration":
                check_iteration(node, rep)
            else:
                check_decision(node, rep)
            check_termination(node, rep)
            check_result(node, rep)
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            # A verifier must be total over traces a model might produce (s5).
            rep.fail("CRASH", "verifier could not complete: %r" % (exc,))
    if solution is not None:
        check_prose(node, solution, rep)
    return rep


def nodes_of(row):
    """D-A1: accept a single node, a list, or a node_id->node map."""
    tn = row.get("trace_nodes")
    if isinstance(tn, dict) and "node_type" in tn:
        return [tn]
    if isinstance(tn, list):
        return list(tn)
    if isinstance(tn, dict):
        return list(tn.values())
    return []


def run(path, verbose=False):
    rows = json.load(open(path, encoding="utf-8"))
    npass = nfail = 0
    unchecked = {}
    for row in rows:
        label = "%s seed=%s" % (row.get("template_id"), row.get("seed"))
        for node in nodes_of(row):
            rep = verify_node(node, label, row.get("solution"))
            for code, msg in rep.unchecked:
                unchecked.setdefault((code, msg), 0)
                unchecked[(code, msg)] += 1
            if rep.ok:
                npass += 1
                if verbose:
                    print("PASS  %s" % label)
            else:
                nfail += 1
                print("FAIL  %s" % label)
                for code, msg in rep.failures:
                    print("        [%s] %s" % (code, msg))
    print("")
    print("--- s8.4 relations reported UNCHECKED (not passed) ---")
    for (code, msg), n in sorted(unchecked.items()):
        print("  [%s] x%d %s" % (code, n, msg))
    print("")
    print("%d passed, %d failed, %d nodes total" % (npass, nfail, npass + nfail))
    return 0 if nfail == 0 else 1


# --------------------------------------------------------------------------
# Negative cases -- sections 8.1 .. 8.7, plus mutations of 4/5/6
# --------------------------------------------------------------------------

def _find(rows, kind):
    key = ("template_normal_depth_iteration" if kind == "iteration"
           else "template_line_balancing_heuristic")
    for r in rows:
        if r["template_id"] == key:
            return copy.deepcopy(r["trace_nodes"])
    raise SystemExit("no %s node in corpus" % kind)


def _find_where(rows, kind, pred):
    key = ("template_normal_depth_iteration" if kind == "iteration"
           else "template_line_balancing_heuristic")
    for r in rows:
        if r["template_id"] == key and pred(r["trace_nodes"]):
            return copy.deepcopy(r["trace_nodes"])
    raise SystemExit("no matching %s node" % kind)


def mutations(rows):
    """(spec clause, description, mutated node, expected failure code)"""
    out = []

    # --- 8.1 gold node with satisfied: false
    n = _find(rows, "iteration"); n["termination"]["satisfied"] = False
    out.append(("s8.1", "gold iteration with termination.satisfied=false", n, "SATISFIED"))
    n = _find(rows, "decision"); n["termination"]["satisfied"] = False
    out.append(("s8.1", "gold decision with termination.satisfied=false", n, "SATISFIED"))

    # --- 8.2 key set differs, in either direction
    n = _find(rows, "iteration"); del n["elements"][0]["change"]
    out.append(("s8.2", "iteration element missing a declared symbol", n, "HOMOGENEITY"))
    n = _find(rows, "iteration"); n["elements"][0]["y_star"] = 1.0
    out.append(("s8.2", "iteration element with an undeclared extra symbol", n, "HOMOGENEITY"))
    n = _find(rows, "decision"); del n["elements"][0]["trials"][0]["closes"]
    out.append(("s8.2", "decision trial missing a declared trial symbol", n, "HOMOGENEITY"))

    # --- 8.3 carry path that does not hold
    n = _find_where(rows, "iteration", lambda x: len(x["elements"]) >= 2)
    n["elements"][1]["y_prev"] = n["elements"][1]["y_prev"] + 0.5
    out.append(("s8.3", "carry y_prev<-y_curr broken between elements 1 and 2", n, "CARRY"))
    n = _find_where(rows, "iteration", lambda x: len(x["elements"]) >= 2)
    n["elements"][1]["g_curr"] = n["elements"][1]["g_curr"] + 1.0
    out.append(("s8.3", "carry g_curr<-evaluation.g broken (dotted path)", n, "CARRY"))

    # --- 8.5 iteration
    n = _find_where(rows, "iteration", lambda x: len(x["elements"]) >= 3)
    n["elements"][0]["converged"] = True
    out.append(("s8.5", "converged=true on a non-final element", n, "ITER_CONV"))
    n = _find(rows, "iteration")
    n["elements"][-1]["evaluation"] = {"y": 1.0, "A": 1.0, "P": 1.0, "AR": 1.0, "g": 0.0}
    out.append(("s8.5", "evaluation non-null on the terminating element", n, "ITER_EVAL"))
    n = _find(rows, "iteration"); n["elements"][0]["g_prev"] = n["elements"][0]["g_curr"]
    out.append(("s8.5", "zero update denominator (g_curr == g_prev)", n, "ITER_DENOM"))
    n = _find(rows, "iteration"); n["elements"][0]["y_next"] = n["elements"][0]["y_next"] + 0.01
    out.append(("s4.1", "secant update relation violated", n, "ITER_UPDATE"))
    n = _find(rows, "iteration"); n["elements"][0]["change"] = 9.9999
    out.append(("s4.2", "change != |y_next - y_curr|", n, "ITER_CHANGE"))
    n = _find(rows, "iteration"); n["preamble"][0]["y"] = 0.25
    out.append(("s4", "preamble start value not the one the question prescribes", n,
                "ITER_PREAMBLE"))

    # --- 8.6 decision
    n = _find(rows, "decision")
    n["elements"][1]["assigned"] = list(n["elements"][0]["assigned"])
    out.append(("s8.6", "overlapping assigned sets", n, "DEC_DISJOINT"))
    n = _find(rows, "decision")
    n["termination"]["universe"] = list(n["termination"]["universe"]) + ["z"]
    out.append(("s8.6", "union of assigned != termination.universe", n, "DEC_UNIVERSE"))
    n = _find_where(rows, "decision",
                    lambda x: any(len([y for y in t["eligible"] if y["fits"]]) >= 2
                                  for e in x["elements"] for t in e["trials"]))
    done = False
    for e in n["elements"]:
        for t in e["trials"]:
            f = [y for y in t["eligible"] if y["fits"]]
            if len(f) >= 2:
                t["chosen"] = min(f, key=lambda y: y["duration"])["task"]
                done = True
                break
        if done:
            break
    out.append(("s8.6", "chosen is not the longest fitting eligible candidate", n, "DEC_RULE"))
    n = _find(rows, "decision"); n["elements"][-1]["station_id"] = 99
    out.append(("s8.6", "non-contiguous station_id", n, "DEC_IDS"))
    n = _find(rows, "decision"); n["elements"][0]["trials"][0]["remaining_after"] += 5
    out.append(("s5.4", "remaining_after != remaining_before - duration(chosen)", n,
                "DEC_BUDGET"))
    n = _find(rows, "decision"); n["elements"][0]["remaining_final"] += 7
    out.append(("s5.7", "remaining_final != capacity - sum(durations)", n, "DEC_FINAL"))
    n = _find(rows, "decision"); n["elements"][0]["remaining_initial"] += 1
    out.append(("s5.1", "remaining_initial != capacity", n, "DEC_OPEN"))
    n = _find(rows, "decision")
    done = False
    for e in n["elements"]:
        for t in e["trials"]:
            if t["chosen"] is None and t["eligible"]:
                t["eligible"][0]["fits"] = True
                done = True
                break
        if done:
            break
    out.append(("s5.3", "closed a station while a candidate still fits", n, "DEC_RULE"))
    n = _find_where(rows, "decision", lambda x: len(x["elements"][0]["assigned"]) >= 2)
    n["elements"][0]["assigned"] = list(reversed(n["elements"][0]["assigned"]))
    out.append(("s5.6", "assigned is not the ordered chosen sequence", n, "DEC_ASSIGNED"))

    # --- 8.7 cardinality_symbol
    n = _find(rows, "decision"); del n["cardinality_symbol"]
    out.append(("s8.7", "answer_bearing with no cardinality_symbol", n, "CARD_SYMBOL"))
    n = _find(rows, "iteration"); n["cardinality_symbol"] = "k"
    out.append(("s8.7", "incidental with a cardinality_symbol", n, "CARD_SYMBOL"))

    # --- 6.1 type/cardinality consistency
    n = _find(rows, "iteration"); n["cardinality"] = "answer_bearing"
    out.append(("s6.1", "iteration claiming cardinality answer_bearing", n, "SHAPE"))
    n = _find(rows, "decision"); n["termination"]["kind"] = "convergence"
    out.append(("s6.1", "decision claiming termination.kind convergence", n, "SHAPE"))
    n = _find(rows, "iteration"); n["node_type"] = "recursion"
    out.append(("s6.1", "unrecognised node_type", n, "SHAPE"))

    # --- 6.3 budget
    n = _find(rows, "iteration"); n["termination"]["max_elements"] = 1
    out.append(("s6.3", "len(elements) > max_elements", n, "BUDGET"))
    n = _find(rows, "iteration"); n["elements"] = []
    out.append(("s6.3", "empty elements on a satisfied node", n, "BUDGET"))

    # --- 6.7 result
    n = _find(rows, "iteration"); n["result"]["value"] = 9.99
    out.append(("s6.7", "iteration result.value != round(last.y_next, dp)", n, "RESULT"))
    n = _find(rows, "decision"); n["result"]["value"] = 99
    out.append(("s6.7", "decision result.value != len(elements)", n, "RESULT"))

    # --- 6.6 termination predicate position
    n = _find_where(rows, "iteration", lambda x: len(x["elements"]) >= 3)
    n["elements"] = n["elements"][:-1]
    out.append(("s6.6", "truncated: predicate does not hold on the last element", n, "TERM"))

    # --- 3 node_id must not encode the element count
    n = _find_where(rows, "decision", lambda x: len(x["elements"]) == 4)
    n["node_id"] = "t19_greedy_4_stations"
    out.append(("s3", "node_id encodes the element count", n, "NODE_ID"))

    return out


def self_test(path):
    rows = json.load(open(path, encoding="utf-8"))
    muts = mutations(rows)
    rejected = missed = 0
    print("=== s8 negative cases ===")
    for clause, desc, node, want in muts:
        rep = verify_node(node, desc)
        codes = {c for c, _ in rep.failures}
        if rep.ok:
            missed += 1
            print("  MISSED  %-6s %s" % (clause, desc))
        elif want not in codes:
            rejected += 1
            print("  REJECT* %-6s %s  (rejected as %s, expected %s)"
                  % (clause, desc, sorted(codes), want))
        else:
            rejected += 1
            print("  REJECT  %-6s %s  [%s]" % (clause, desc, want))
    print("")
    print("%d/%d negative cases rejected, %d missed" % (rejected, len(muts), missed))

    # Control: the unmutated gold nodes must still pass.
    print("")
    print("=== control: unmutated gold ===")
    ctrl_fail = 0
    total = 0
    for row in rows:
        for node in nodes_of(row):
            total += 1
            if not verify_node(node, "", row.get("solution")).ok:
                ctrl_fail += 1
    print("  %d/%d gold nodes pass" % (total - ctrl_fail, total))
    return 0 if missed == 0 and ctrl_fail == 0 else 1


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = {a for a in argv[1:] if a.startswith("--")}
    path = args[0] if args else "docs/re-implementation-sep/track-a/phase3_conformance/traces.json"
    if "--self-test" in flags:
        return self_test(path)
    return run(path, verbose="--verbose" in flags)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
