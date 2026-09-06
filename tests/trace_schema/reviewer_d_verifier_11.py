"""Reviewer D round-2 verifier for schema 1.1.

Strictly roles-driven: no literal template symbol name appears anywhere in the
checking logic for `iteration`.  Used to test the coordinator's round-2 claims,
in particular F3 (is the type reusable?) and F1 (are precedences sufficient?).

    python -m tests.trace_schema.reviewer_d_verifier_11 <traces.json>
    python -m tests.trace_schema.reviewer_d_verifier_11 <traces.json> --rename
    python -m tests.trace_schema.reviewer_d_verifier_11 <traces.json> --prec
"""

from __future__ import annotations

import ast
import copy
import decimal
import json
import sys

ITER_ROLES = ["index", "iterate_prev", "iterate_curr", "residual_prev",
              "residual_curr", "iterate_next", "change", "converged", "evaluation"]
DEC_ROLES = ["station_id", "capacity", "remaining_initial", "trials",
             "assigned", "remaining_final"]


class Rep:
    def __init__(self):
        self.failures = []
        self.unchecked = []

    def fail(self, c, m):
        self.failures.append((c, m))

    @property
    def ok(self):
        return not self.failures


def half_up(x, dp):
    """s6 `rounding: decimal-half-up` -- explicitly not binary round()."""
    q = decimal.Decimal(1).scaleb(-dp)
    return float(decimal.Decimal(repr(float(x))).quantize(q, decimal.ROUND_HALF_UP))


def tol_of(node, symbol):
    """s7.1 -- from symbol_precision, never inferred from data."""
    sp = node.get("symbol_precision") or {}
    if symbol not in sp:
        return None
    return 0.5 * 10 ** (-sp[symbol]) + 1e-12


# --- s4.2 update_relation: closed grammar over role names -------------------

_OPS = (ast.Add, ast.Sub, ast.Mult, ast.Div)


def eval_relation(expr, env, rep, tag):
    """Role identifiers, decimal literals, + - * /, unary -, parens.  Nothing else."""
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError:
        rep.fail("REL_GRAMMAR", "%s: update_relation does not parse" % tag)
        return None

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.BinOp):
            if not isinstance(n.op, _OPS):
                raise ValueError("operator not in grammar")
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if b == 0:
                raise ZeroDivisionError("zero denominator")
            return a / b
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return -ev(n.operand)
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return float(n.value)
        if isinstance(n, ast.Name):
            if n.id not in env:
                raise ValueError("name %r is not a role of this node" % n.id)
            return float(env[n.id])
        raise ValueError("construct not in grammar: %s" % type(n).__name__)

    try:
        return ev(tree)
    except ZeroDivisionError:
        rep.fail("ITER_DENOM", "%s: zero update denominator (s4.3.7, s8.5)" % tag)
        return None
    except (ValueError, TypeError) as e:
        rep.fail("REL_GRAMMAR", "%s: %s (s4.2)" % (tag, e))
        return None


def keys_eq(got, want):
    return set(got) == set(want)


def verify(node):
    rep = Rep()
    if not isinstance(node, dict):
        rep.fail("SHAPE", "not an object")
        return rep

    # s6.1 version + shape
    sv = str(node.get("schema_version", ""))
    if not sv:
        rep.fail("SHAPE", "schema_version missing (s3.1)")
        return rep
    if sv.split(".")[0] != "1":
        rep.fail("VERSION", "unknown major version %r (s8.9)" % sv)
        return rep
    nt = node.get("node_type")
    if nt not in ("iteration", "decision"):
        rep.fail("SHAPE", "unrecognised node_type %r" % (nt,))
        return rep

    want_card = "incidental" if nt == "iteration" else "answer_bearing"
    want_kind = "convergence" if nt == "iteration" else "exhaustion"
    term = node.get("termination") or {}
    if node.get("cardinality") != want_card:
        rep.fail("SHAPE", "cardinality != %r" % want_card)
    if term.get("kind") != want_kind:
        rep.fail("SHAPE", "termination.kind != %r" % want_kind)
    has_cs = node.get("cardinality_symbol") is not None
    if (node.get("cardinality") == "answer_bearing") != has_cs:
        rep.fail("CARD_SYMBOL", "cardinality_symbol must be present iff answer_bearing (s8.7)")
    for f in (("roles",) if nt == "iteration" else ()) + ("rounding", "symbol_precision", "element_symbols",
              "elements", "carry", "carry_notes", "result"):
        if f not in node:
            rep.fail("SHAPE", "missing required field %r (s3.1)" % f)
    if term.get("satisfied") is not True:
        rep.fail("SATISFIED", "gold node with satisfied=%r (s8.1)" % term.get("satisfied"))

    els = node.get("elements")
    if not isinstance(els, list) or not els:
        rep.fail("BUDGET", "elements must be a list of length >= 1 (s6.4)")
        return rep
    if len(els) > term.get("max_elements", 10 ** 9):
        rep.fail("BUDGET", "len(elements) > max_elements (s6.4)")

    esyms = node.get("element_symbols") or []
    roles = node.get("roles") or {}

    # s6.2 roles
    mand = ITER_ROLES if nt == "iteration" else DEC_ROLES
    if nt == "iteration":
        for r in mand:
            if r not in roles:
                rep.fail("ROLES", "mandatory role %r missing from roles (s8.8)" % r)
        for r, s in roles.items():
            if s not in esyms:
                rep.fail("ROLES", "roles[%r]=%r not in element_symbols (s8.8)" % (r, s))

    def RV(el, role):
        return el[roles[role]]

    # s6.3 homogeneity
    for i, e in enumerate(els):
        if not keys_eq(e, esyms):
            rep.fail("HOMOGENEITY", "element %d key set != element_symbols (s8.2)" % (i + 1))
    if "preamble" in node:
        for i, f in enumerate(node["preamble"]):
            if not keys_eq(f, node.get("preamble_symbols") or []):
                rep.fail("HOMOGENEITY", "preamble frame %d != preamble_symbols" % (i + 1))

    # s6.5 carry (paths only) + carry_notes on the unchecked channel
    for lhs, rhs in (node.get("carry") or {}).items():
        segs = rhs.split(".")
        if segs[0] not in esyms:
            rep.fail("CARRY", "carry path %r does not start at an element_symbol (s3.5)" % rhs)
            continue
        for k in range(len(els) - 1):
            cur = els[k]
            ok = True
            for s in segs:
                if not isinstance(cur, dict) or s not in cur:
                    ok = False
                    break
                cur = cur[s]
            if not ok:
                rep.fail("CARRY", "carry path %r unresolvable at element %d" % (rhs, k + 1))
                continue
            got = els[k + 1].get(lhs)
            t = tol_of(node, lhs)
            same = got == cur
            if not same and t is not None and isinstance(got, (int, float)) \
                    and isinstance(cur, (int, float)):
                same = abs(got - cur) <= t
            if not same:
                rep.fail("CARRY", "carry %r<-%r broken at %d/%d (s8.3)" % (lhs, rhs, k + 1, k + 2))
    for lhs, txt in (node.get("carry_notes") or {}).items():
        rep.unchecked.append(("CARRY_NOTE", "carry_notes[%r] UNCHECKED (s8.4)" % lhs))

    if nt == "iteration":
        check_iter(node, rep, roles, RV)
    else:
        check_dec(node, rep)
    check_result(node, rep, roles)
    return rep


def check_iter(node, rep, roles, RV):
    els = node["elements"]
    term = node["termination"]
    tol = term.get("tolerance")
    cmpr = term.get("comparison", "lt")
    esyms_eval = node.get("evaluation_symbols") or []
    eroles = node.get("evaluation_roles") or {}
    for i, e in enumerate(els):
        tag = "element %d" % (i + 1)
        env = {r: e[s] for r, s in roles.items()
               if isinstance(e.get(s), (int, float)) and not isinstance(e.get(s), bool)}
        # s4.3.1
        pred = eval_relation(node.get("update_relation", ""), env, rep, tag)
        t_next = tol_of(node, roles["iterate_next"])
        if pred is not None and t_next is not None:
            if abs(pred - RV(e, "iterate_next")) > t_next:
                rep.fail("ITER_UPDATE", "%s: update_relation != iterate_next (s4.3.1)" % tag)
        # s4.3.2
        want = abs(RV(e, "iterate_next") - RV(e, "iterate_curr"))
        t_ch = tol_of(node, roles["change"])
        if t_ch is not None and abs(want - RV(e, "change")) > t_ch:
            rep.fail("ITER_CHANGE", "%s: change != |iterate_next-iterate_curr| (s4.3.2)" % tag)
        # s4.3.3/4
        conv = RV(e, "converged")
        want_conv = (RV(e, "change") < tol) if cmpr == "lt" else (RV(e, "change") <= tol)
        if conv != want_conv:
            rep.fail("ITER_CONV", "%s: converged != predicate (s4.3.3)" % tag)
        if conv != (i == len(els) - 1):
            rep.fail("ITER_CONV", "%s: converged true only on the last element (s4.3.4)" % tag)
        # s4.3.5/6
        ev = RV(e, "evaluation")
        if (ev is None) != bool(conv):
            rep.fail("ITER_EVAL", "%s: evaluation null iff converged (s4.3.5)" % tag)
        if ev is not None:
            if not keys_eq(ev, esyms_eval):
                rep.fail("HOMOGENEITY", "%s: evaluation frame != evaluation_symbols (s6.3)" % tag)
            if "iterate" in eroles and eroles["iterate"] in ev:
                if ev[eroles["iterate"]] != RV(e, "iterate_next"):
                    rep.fail("ITER_FRAME",
                             "%s: evaluation.iterate != iterate_next (s4.3.6, s8.5)" % tag)
        # s4.3.8
        if RV(e, "index") != i + 1:
            rep.fail("ITER_INDEX", "%s: index not contiguous 1..n (s4.3.8, s8.5)" % tag)
    # s6.7 termination position
    for i, e in enumerate(els):
        holds = RV(e, "change") < tol if cmpr == "lt" else RV(e, "change") <= tol
        if holds != (i == len(els) - 1):
            rep.fail("TERM", "predicate holds on element %d only if last (s6.7)" % (i + 1))
    # s4.4 preamble seeds the first element
    if "preamble" in node and len(node["preamble"]) >= 2 and "iterate" in (
            node.get("evaluation_roles") or {}):
        ir = node["evaluation_roles"]["iterate"]
        rr = node["evaluation_roles"].get("residual")
        p = node["preamble"]
        for j, role in ((0, "iterate_prev"), (1, "iterate_curr")):
            if p[j].get(ir) != RV(els[0], role):
                rep.fail("ITER_PREAMBLE", "preamble[%d].%s does not seed %s (s4.4)"
                         % (j, ir, role))
        if rr:
            for j, role in ((0, "residual_prev"), (1, "residual_curr")):
                if p[j].get(rr) != RV(els[0], role):
                    rep.fail("ITER_PREAMBLE", "preamble[%d].%s does not seed %s (s4.4)"
                             % (j, rr, role))


def check_dec(node, rep):
    els = node["elements"]
    term = node["termination"]
    uni = list(term.get("universe") or [])
    prec = term.get("precedences") or {}
    sel = node.get("selection") or {}
    filt, crit = sel.get("filter"), sel.get("criterion")
    obj, tb = sel.get("objective", "max"), sel.get("tie_break", "universe_order")
    tsyms = node.get("trial_symbols") or []
    elsyms = node.get("eligible_symbols") or []

    committed = []
    for i, e in enumerate(els):
        tag = "station %s" % e.get("station_id")
        if e.get("remaining_initial") != e.get("capacity"):
            rep.fail("DEC_OPEN", "%s: remaining_initial != capacity (s5.4.1)" % tag)
        trials = e.get("trials") or []
        if not trials:
            rep.fail("DEC_TRIALS", "%s: trials empty (s5.1)" % tag)
            continue
        budget = e.get("remaining_initial")
        seq = []
        for j, tr in enumerate(trials):
            ttag = "%s trial %d" % (tag, j + 1)
            if not keys_eq(tr, tsyms):
                rep.fail("HOMOGENEITY", "%s: keys != trial_symbols (s6.3)" % ttag)
            elig = tr.get("eligible") or []
            for x in elig:
                if not keys_eq(x, elsyms):
                    rep.fail("HOMOGENEITY", "%s: eligible entry != eligible_symbols (s6.3)" % ttag)
            by = {x["task"]: x for x in elig if "task" in x}
            # filter symbol is honest
            for x in elig:
                if bool(x.get(filt)) != (x.get(crit) <= budget):
                    rep.fail("DEC_FITS", "%s: %r filter symbol disagrees with budget (s5.1)"
                             % (ttag, x.get("task")))
            # --- s5.4.3b: eligible checked in BOTH directions
            done = set(committed) | set(seq)
            for x in elig:
                miss = [p for p in prec.get(x["task"], []) if p not in done]
                if miss:
                    rep.fail("DEC_PREC_IN",
                             "%s: %r is eligible but its precedences %r are uncommitted (s5.4.3b)"
                             % (ttag, x["task"], miss))
            for item in uni:
                if item in done:
                    continue
                if all(p in done for p in prec.get(item, [])) and item not in by:
                    rep.fail("DEC_PREC_OUT",
                             "%s: %r has all precedences satisfied but is absent from eligible "
                             "(s5.4.3b)" % (ttag, item))
            ch = tr.get("chosen")
            fitting = [x for x in elig if x.get(filt)]
            if ch is None:
                if fitting:
                    rep.fail("DEC_RULE", "%s: chose nothing while candidates fit (s5.4.4)" % ttag)
                if tr.get("remaining_after") != budget:
                    rep.fail("DEC_BUDGET", "%s: budget moved on a null choice (s5.4.5)" % ttag)
            else:
                if ch not in by or not by[ch].get(filt):
                    rep.fail("DEC_RULE", "%s: chosen absent or unfiltered (s5.4.3)" % ttag)
                elif fitting:
                    best = (max if obj == "max" else min)(x[crit] for x in fitting)
                    if by[ch][crit] != best:
                        rep.fail("DEC_RULE", "%s: chosen does not optimise selection (s5.4.4)"
                                 % ttag)
                    else:
                        tied = [x["task"] for x in fitting if x[crit] == best]
                        if len(tied) > 1:
                            key = (uni.index if tb == "universe_order"
                                   else [x["task"] for x in elig].index)
                            if ch != min(tied, key=key):
                                rep.fail("DEC_TIE", "%s: tie_break not applied (s5.4.4)" % ttag)
                    if tr.get("remaining_after") != budget - by[ch][crit]:
                        rep.fail("DEC_BUDGET", "%s: remaining_after wrong (s5.4.5)" % ttag)
                seq.append(ch)
            # s5.1: closes == (chosen is null), one rule
            if tr.get("closes") != (ch is None):
                rep.fail("DEC_CLOSES", "%s: closes != (chosen is null) (s5.1)" % ttag)
            budget = tr.get("remaining_after")
        if list(e.get("assigned") or []) != seq:
            rep.fail("DEC_ASSIGNED", "%s: assigned != ordered chosen (s5.4.6)" % tag)
        durs = {x["task"]: x[crit] for tr in trials for x in (tr.get("eligible") or [])
                if "task" in x}
        if e.get("remaining_final") != e.get("capacity") - sum(
                durs.get(t, 0) for t in (e.get("assigned") or [])):
            rep.fail("DEC_FINAL", "%s: remaining_final wrong (s5.4.7)" % tag)
        committed.extend(e.get("assigned") or [])

    # s5.3 global
    if [e.get("station_id") for e in els] != list(range(1, len(els) + 1)):
        rep.fail("DEC_IDS", "station_id not 1..n (s5.3)")
    if len(committed) != len(set(committed)):
        rep.fail("DEC_DISJOINT", "assigned sets overlap (s5.3)")
    if set(committed) != set(uni):
        rep.fail("DEC_UNIVERSE", "union != universe (s5.3)")
    if term.get("max_elements", 0) < len(uni):
        rep.fail("DEC_BUDGETREL", "max_elements < len(universe) (s5.3)")
    if node.get("capacity_constant") is True:
        if len({e.get("capacity") for e in els}) > 1:
            rep.fail("DEC_CAPCONST", "capacity varies under capacity_constant (s5.4.8, s8.6)")
    # s6.7 exhaustion, cumulative scope
    if term.get("scope") == "cumulative":
        acc = term.get("accumulator_symbol")
        for i in range(len(els) - 1):
            run = set()
            for x in els[:i + 1]:
                run |= set(x.get(acc) or [])
            if run == set(uni):
                rep.fail("TERM", "universe exhausted before the last element (s6.7)")


def check_result(node, rep, roles):
    res = node.get("result") or {}
    frm = res.get("from")
    els = node["elements"]
    dp = res.get("dp")
    if frm == "len(elements)":
        if node.get("cardinality") != "answer_bearing":
            rep.fail("RESULT", "from=len(elements) requires answer_bearing (s3.6)")
        if res.get("value") != len(els):
            rep.fail("RESULT", "result.value != len(elements) (s6.8)")
    elif isinstance(frm, str) and frm.startswith("last."):
        role = frm[5:]
        if role not in roles:
            rep.fail("RESULT", "result.from names role %r absent from roles (s3.6)" % role)
            return
        want = half_up(els[-1][roles[role]], dp)
        if abs(float(res.get("value")) - want) > 0.5 * 10 ** (-(dp + 1)):
            rep.fail("RESULT", "result.value != %s rounded to dp (s6.8)" % frm)
    else:
        rep.fail("RESULT", "result.from %r not in the closed set (s3.6)" % (frm,))


def run(path):
    rows = json.load(open(path, encoding="utf-8"))
    p = f = 0
    unch = {}
    for r in rows:
        rep = verify(r["trace_nodes"])
        for c, m in rep.unchecked:
            unch[(c, m)] = unch.get((c, m), 0) + 1
        if rep.ok:
            p += 1
        else:
            f += 1
            print("FAIL %s seed=%s" % (r["template_id"], r["seed"]))
            for c, m in rep.failures[:6]:
                print("       [%s] %s" % (c, m))
    print("")
    for (c, m), n in sorted(unch.items()):
        print("  [%s] x%d %s" % (c, n, m))
    print("")
    print("%d passed, %d failed, %d total" % (p, f, p + f))
    return 0 if f == 0 else 1


# --------------------------------------------------------------------------
# Break attempt 1 (F3): rename every symbol of an iteration node.
# --------------------------------------------------------------------------

def rename_test(path):
    rows = json.load(open(path, encoding="utf-8"))
    src = next(r for r in rows if r["template_id"] == "template_normal_depth_iteration")
    n = copy.deepcopy(src["trace_nodes"])
    # A hypothetical routing template: iterate `S`, residual `phi`, index `step`.
    m = {"k": "step", "y_prev": "S_old", "y_curr": "S_cur", "y_next": "S_new",
         "g_prev": "phi_old", "g_curr": "phi_cur", "change": "delta",
         "converged": "done", "evaluation": "probe"}
    em = {"y": "S", "A": "area", "P": "wetted", "AR": "AR", "g": "phi"}
    n["node_id"] = "tXX_reservoir_routing_step"
    n["element_symbols"] = [m[s] for s in n["element_symbols"]]
    n["roles"] = {r: m[s] for r, s in n["roles"].items()}
    n["evaluation_symbols"] = [em[s] for s in n["evaluation_symbols"]]
    n["preamble_symbols"] = [em[s] for s in n["preamble_symbols"]]
    n["evaluation_roles"] = {r: em[s] for r, s in n["evaluation_roles"].items()}
    n["symbol_precision"] = {m.get(k, em.get(k, k)): v
                             for k, v in n["symbol_precision"].items()}
    n["carry"] = {m[k]: ".".join([m[v.split(".")[0]]] +
                                 [em[x] for x in v.split(".")[1:]])
                  for k, v in n["carry"].items()}
    n["elements"] = [{m[k]: ({em[a]: b for a, b in v.items()} if k == "evaluation" and v
                             else v) for k, v in e.items()} for e in n["elements"]]
    n["preamble"] = [{em[k]: v for k, v in f.items()} for f in n["preamble"]]
    n["result"] = dict(n["result"], symbol="S_out")
    # update_relation is untouched: it is written over ROLES, not symbols.
    rep = verify(n)
    print("=== F3 break attempt: iteration node with every symbol renamed ===")
    print("    element_symbols :", n["element_symbols"])
    print("    update_relation :", n["update_relation"])
    print("    result          :", n["result"])
    print("    verdict         :", "PASSES (type is reusable)" if rep.ok
          else "FAILS -- F3 not addressed")
    for c, m_ in rep.failures[:8]:
        print("      [%s] %s" % (c, m_))

    # Control: the same rename must still catch a corrupted update.
    n2 = copy.deepcopy(n)
    n2["elements"][0]["S_new"] += 0.01
    r2 = verify(n2)
    print("    control (broken update on the renamed node):",
          "rejected [%s]" % sorted({c for c, _ in r2.failures}) if not r2.ok
          else "ACCEPTED -- verifier is vacuous")
    return 0


# --------------------------------------------------------------------------
# Break attempt 2 (F1): is `precedences` sufficient, and does it agree?
# --------------------------------------------------------------------------

def prec_test(path):
    rows = json.load(open(path, encoding="utf-8"))
    lb = [r for r in rows if r["template_id"] == "template_line_balancing_heuristic"]
    print("=== F1 break attempt: precedences sufficiency and agreement ===")
    precs = {json.dumps(r["trace_nodes"]["termination"].get("precedences"), sort_keys=True)
             for r in lb}
    print("    distinct precedence relations across 40 traces:", len(precs))
    for p in precs:
        print("      ", p)
    # Does the DAG cover the universe, and is it acyclic?
    bad = 0
    for r in lb:
        t = r["trace_nodes"]["termination"]
        pr, uni = t.get("precedences") or {}, t["universe"]
        if set(pr) != set(uni):
            bad += 1
        # topological order check against the recorded assignment order
        order, seen = [], set()
        for e in r["trace_nodes"]["elements"]:
            order.extend(e["assigned"])
        for it in order:
            if any(p not in seen for p in pr.get(it, [])):
                bad += 1
            seen.add(it)
    print("    traces whose precedences do not key exactly the universe, or whose")
    print("    recorded assignment order violates them:", bad, "of", len(lb))

    # Both directions, over every trial in the corpus.
    fwd = rev = trials = 0
    for r in lb:
        n = r["trace_nodes"]
        pr, uni = n["termination"]["precedences"], n["termination"]["universe"]
        done = []
        for e in n["elements"]:
            for tr in e["trials"]:
                trials += 1
                d = set(done)
                present = {x["task"] for x in tr["eligible"]}
                for x in tr["eligible"]:
                    if any(p not in d for p in pr.get(x["task"], [])):
                        fwd += 1
                for it in uni:
                    if it in d:
                        continue
                    if all(p in d for p in pr.get(it, [])) and it not in present:
                        rev += 1
                if tr["chosen"] is not None:
                    done.append(tr["chosen"])
    print("    trials examined:", trials)
    print("    direction A (eligible entry with unmet precedences):", fwd)
    print("    direction B (ready item absent from eligible)      :", rev)

    # Negative controls: both directions must be detectable.
    src = copy.deepcopy(lb[0]["trace_nodes"])
    n1 = copy.deepcopy(src)
    tgt = [i for i in n1["termination"]["universe"]
           if n1["termination"]["precedences"].get(i)][-1]
    n1["elements"][0]["trials"][0]["eligible"].append(
        {"task": tgt, "duration": 1, "fits": True})
    a = verify(n1)
    n2 = copy.deepcopy(src)
    n2["elements"][0]["trials"][0]["eligible"] = [
        x for x in n2["elements"][0]["trials"][0]["eligible"]
        if x["task"] != n2["elements"][0]["trials"][0]["eligible"][0]["task"]]
    b = verify(n2)
    print("    control A (inject a not-yet-ready item into eligible):",
          "rejected %s" % sorted({c for c, _ in a.failures}) if not a.ok else "ACCEPTED")
    print("    control B (delete a ready item from eligible)       :",
          "rejected %s" % sorted({c for c, _ in b.failures}) if not b.ok else "ACCEPTED")
    return 0


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = {a for a in argv[1:] if a.startswith("--")}
    path = args[0]
    if "--rename" in flags:
        return rename_test(path)
    if "--prec" in flags:
        return prec_test(path)
    return run(path)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
