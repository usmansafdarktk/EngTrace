"""Reviewer D2 round 3 -- verifier for D3.3 schema version 1.5.

Adapted from `reviewer_d2_verifier.py` (round 2, schema 1.1), which is left
untouched as the round-2 record.  Grammar and rounding helpers are reused from
it; every check below is rewritten against 1.2's role maps.

Conforms to 8B.11: no element, trial or candidate value is resolved by a
literal symbol name.  Every access goes through `roles`, `trial_roles`,
`candidate_roles` or `evaluation_roles`.

Usage:
    python -m tests.trace_schema.reviewer_d2_verifier_15 <traces.json> [--mode gold|candidate]
"""

from __future__ import annotations

import json
import sys

from tests.trace_schema.reviewer_d2_verifier import (
    Failure, Report, Parser, tokenize, eval_expr, parse_path, select_path,
    do_round, close, keyset_check, tolerance_for,
)

KNOWN_MAJOR = "1"

ITERATION_ROLES = (
    "index", "iterate_prev", "iterate_curr", "iterate_next",
    "residual_prev", "residual_curr", "change", "converged", "evaluation",
)
ITERATION_EVAL_ROLES = ("iterate", "residual")
# Section 5.1, 1.2 -- ROLES now, not symbols.  The third column of those tables
# ("Reference symbol") is deliberately not reproduced here: 8B.11 forbids using it.
DECISION_ROLES = ("index", "budget_total", "budget_initial",
                  "trials", "committed", "budget_final")
DECISION_TRIAL_ROLES = ("candidates", "chosen", "budget_after", "closes")
DECISION_CANDIDATE_ROLES = ("item", "measure", "admissible")

COMPARATORS = ("<=", ">=", "==", "<", ">")


def parse_filter_relation(src):
    """Section 5.2 -- a §4.2 expression, a comparator, a second §4.2 expression."""
    for op in COMPARATORS:
        idx = src.find(op)
        if idx == -1:
            continue
        # avoid splitting "<=" at "<"
        if op in ("<", ">") and idx + 1 < len(src) and src[idx + 1] == "=":
            continue
        lhs, rhs = src[:idx], src[idx + len(op):]
        return (Parser(tokenize(lhs)).parse(), op, Parser(tokenize(rhs)).parse())
    raise Failure("5.2: no comparator in filter_relation %r" % (src,))


def eval_filter_relation(tree, env):
    lhs, op, rhs = tree
    a, b = eval_expr(lhs, env), eval_expr(rhs, env)
    return {"<=": a <= b, "<": a < b, ">=": a >= b, ">": a > b,
            "==": a == b}[op]


def verify_node(node, label="node", mode="gold"):
    rep = Report(label)
    try:
        _verify(node, rep, mode)
    except Failure as e:
        if not rep.failures:
            rep.fail("6", str(e))
    except (KeyError, TypeError, IndexError, ValueError) as e:
        rep.fail("6", "%s while verifying: %s" % (type(e).__name__, e))
    return rep


def _verify(node, rep, mode):
    if not isinstance(node, dict):
        rep.fail("3.0", "node is not a JSON object")
        return

    # --- step 1: version and shape ---------------------------------------
    sv = node.get("schema_version")
    if not isinstance(sv, str):
        rep.fail("8A.8", "schema_version missing or not a string")
        return
    if sv.split(".")[0] != KNOWN_MAJOR:
        rep.fail("8A.8", "unknown major schema version %r" % sv)
        return

    ntype = node.get("node_type")
    if ntype not in ("iteration", "decision"):
        rep.fail("6.1", "unrecognised node_type %r" % (ntype,))
        return

    required = ["schema_version", "node_type", "node_id", "cardinality", "roles",
                "termination", "rounding", "symbol_precision", "element_symbols",
                "elements", "carry", "carry_notes", "result"] + {
        "iteration": ["evaluation_symbols", "evaluation_roles", "update_relation",
                      "constants", "frame_relations"],
        "decision": ["trial_symbols", "trial_roles", "candidate_roles",
                     "eligible_symbols", "selection"],
    }[ntype]
    for f in required:
        if f not in node:
            rep.fail("6.1", "required field %r absent" % f)
    if rep.failures:
        return

    fixed = {"iteration": ("incidental", "convergence"),
             "decision": ("answer_bearing", "exhaustion")}[ntype]
    if node["cardinality"] != fixed[0]:
        rep.fail("6.1", "cardinality %r; %s fixes %r"
                 % (node["cardinality"], ntype, fixed[0]))
    term = node["termination"]
    if term.get("kind") != fixed[1]:
        rep.fail("6.1", "termination.kind %r; %s fixes %r"
                 % (term.get("kind"), ntype, fixed[1]))

    ab = node["cardinality"] == "answer_bearing"
    if ab != ("cardinality_symbol" in node):
        rep.fail("8A.6", "cardinality_symbol presence disagrees with cardinality")
    if ("preamble" in node) != ("preamble_symbols" in node):
        rep.fail("3.1", "preamble and preamble_symbols must appear together")
    for f in ("satisfied", "max_elements"):
        if f not in term:
            rep.fail("3.4", "termination.%s absent" % f)
    for f in {"convergence": ("quantity_role", "comparison", "tolerance"),
              "exhaustion": ("scope", "accumulator_role", "universe",
                             "precedences", "item_measures", "budget")}.get(term.get("kind"), ()):
        if f not in term:
            rep.fail("3.4", "termination.%s absent for kind %r" % (f, term["kind"]))
    if rep.failures:
        return

    esyms, elements = node["element_symbols"], node["elements"]

    # --- step 2: roles ----------------------------------------------------
    roles = node["roles"]
    maps = [("roles", roles,
             ITERATION_ROLES if ntype == "iteration" else DECISION_ROLES, esyms)]
    if ntype == "iteration":
        maps.append(("evaluation_roles", node["evaluation_roles"],
                     ITERATION_EVAL_ROLES, node["evaluation_symbols"]))
    else:
        maps.append(("trial_roles", node["trial_roles"],
                     DECISION_TRIAL_ROLES, node["trial_symbols"]))
        maps.append(("candidate_roles", node["candidate_roles"],
                     DECISION_CANDIDATE_ROLES, node["eligible_symbols"]))
    for mname, m, mandatory, syms in maps:
        for r in mandatory:
            if r not in m:
                rep.fail("8A.7", "mandatory role %r missing from %s" % (r, mname))
        for r, sym in m.items():
            if sym not in syms:
                rep.fail("8A.7", "%s[%r] = %r not in the declared symbol set"
                         % (mname, r, sym))
    if rep.failures:
        return

    def E(el, role):
        return el[roles[role]]

    # --- step 3: homogeneity ---------------------------------------------
    for i, el in enumerate(elements):
        keyset_check(rep, "8A.2", "element[%d]" % i, el, esyms)
    if "preamble" in node:
        for i, fr in enumerate(node["preamble"]):
            keyset_check(rep, "8A.2", "preamble[%d]" % i, fr,
                         node["preamble_symbols"])
    if rep.failures:
        return
    if ntype == "iteration":
        for i, el in enumerate(elements):
            ev = E(el, "evaluation")
            if ev is not None:
                keyset_check(rep, "8A.2", "element[%d].evaluation" % i, ev,
                             node["evaluation_symbols"])
    else:
        tr_r, ca_r = node["trial_roles"], node["candidate_roles"]
        for i, el in enumerate(elements):
            trials = E(el, "trials")
            if not isinstance(trials, list) or not trials:
                rep.fail("5.1", "element[%d] trials empty or not a list" % i)
                continue
            for j, tr in enumerate(trials):
                keyset_check(rep, "8A.2", "element[%d].trial[%d]" % (i, j), tr,
                             node["trial_symbols"])
                cands = tr.get(tr_r["candidates"])
                if isinstance(cands, list):
                    for m, en in enumerate(cands):
                        keyset_check(rep, "8A.2",
                                     "element[%d].trial[%d].candidate[%d]" % (i, j, m),
                                     en, node["eligible_symbols"])
    if rep.failures:
        return

    # --- step 4: budget ---------------------------------------------------
    if len(elements) < 1:
        rep.fail("6.4", "elements is empty; step 4 requires >= 1")
    elif len(elements) > term["max_elements"]:
        # Section 6.0: an invariant in gold mode, an observation in candidate mode.
        if mode == "gold":
            rep.fail("6.4", "len(elements)=%d exceeds max_elements=%d"
                     % (len(elements), term["max_elements"]))
        else:
            rep.note_unchecked("7.2: len(elements)=%d over budget %d, reported "
                               "not failed" % (len(elements), term["max_elements"]))
    if rep.failures:
        return

    # --- step 5: carry ----------------------------------------------------
    for target, path in node["carry"].items():
        if target not in esyms:
            rep.fail("8A.3", "carry target %r is not an element symbol" % target)
            continue
        try:
            parts = parse_path(path, esyms)
        except Failure as e:
            rep.fail("8A.3", str(e))
            continue
        for k in range(len(elements) - 1):
            try:
                got = select_path(parts, elements[k])
            except Failure as e:
                rep.fail("8A.3", "element[%d]: %s" % (k, e))
                continue
            want = elements[k + 1][target]
            ok = (got == want) if (isinstance(got, bool)
                                   or not isinstance(got, (int, float))) \
                else (isinstance(want, (int, float)) and float(got) == float(want))
            if not ok:
                rep.fail("8A.3", "carry %r <- %r broken between element %d and %d"
                         % (target, path, k, k + 1))
    for k, v in node["carry_notes"].items():
        rep.note_unchecked("carry_notes[%r] = %r (8B.10: unchecked)" % (k, v))

    # --- step 6: local invariants ----------------------------------------
    (_iteration_locals if ntype == "iteration" else _decision_locals)(node, rep, roles)
    if rep.failures:
        return

    # --- step 7: termination ---------------------------------------------
    _termination(node, rep, roles, ntype, mode)
    if rep.failures:
        return

    # --- step 8: result ---------------------------------------------------
    _result(node, rep, roles)



# --------------------------------------------------------------------------
# Section 4.6 -- `constants` and `frame_relations`.
# Grammar: 4.2 plus `**`, sqrt(x), round(x, n).  Names resolve to constants,
# to the frame's `iterate` role, or to an EARLIER frame_relations symbol.
# --------------------------------------------------------------------------

def tok46(s):
    toks, i = [], 0
    while i < len(s):
        c = s[i]
        if c in " \t\n\r":
            i += 1
        elif s.startswith("**", i):
            toks.append(("**", "**")); i += 2
        elif c in "+-*/(),":
            toks.append((c, c)); i += 1
        elif c.isdigit() or (c == "." and i + 1 < len(s) and s[i + 1].isdigit()):
            j = i
            while j < len(s) and (s[j].isdigit() or s[j] == "."):
                j += 1
            toks.append(("num", s[i:j])); i = j
        elif c.isalpha() or c == "_":
            j = i
            while j < len(s) and (s[j].isalnum() or s[j] == "_"):
                j += 1
            toks.append(("name", s[i:j])); i = j
        else:
            raise Failure("4.6: character %r is not in the grammar" % c)
    return toks


class P46:
    def __init__(self, toks):
        self.t, self.i = toks, 0

    def peek(self):
        return self.t[self.i][0] if self.i < len(self.t) else None

    def take(self):
        x = self.t[self.i]; self.i += 1; return x

    def parse(self):
        n = self.expr()
        if self.i != len(self.t):
            raise Failure("4.6: trailing tokens")
        return n

    def expr(self):
        n = self.term()
        while self.peek() in ("+", "-"):
            n = (self.take()[0], n, self.term())
        return n

    def term(self):
        n = self.unary()
        while self.peek() in ("*", "/"):
            n = (self.take()[0], n, self.unary())
        return n

    def unary(self):
        if self.peek() == "-":
            self.take(); return ("neg", self.unary())
        return self.power()

    def power(self):
        b = self.atom()
        if self.peek() == "**":
            self.take()
            return ("**", b, self.unary())
        return b

    def atom(self):
        if self.peek() is None:
            raise Failure("4.6: unexpected end of expression")
        k, x = self.take()
        if k == "(":
            n = self.expr()
            if self.peek() != ")":
                raise Failure("4.6: unbalanced parenthesis")
            self.take(); return n
        if k == "num":
            return ("lit", float(x))
        if k == "name":
            if self.peek() == "(":
                self.take()
                a = self.expr()
                if x == "sqrt":
                    if self.peek() != ")":
                        raise Failure("4.6: sqrt takes one argument")
                    self.take(); return ("sqrt", a)
                if x == "round":
                    if self.peek() != ")":
                        raise Failure("4.6: round takes exactly one argument "
                                      "(1.5: the digit comes from symbol_precision)")
                    self.take()
                    if a[0] != "name":
                        raise Failure("4.6: round(sym) takes a symbol, not an "
                                      "expression -- its digit is "
                                      "symbol_precision[sym]")
                    return ("round", a)
                raise Failure("4.6: %r is not a permitted function" % x)
            return ("name", x)
        raise Failure("4.6: token %r is not in the grammar" % (x,))


class FrameEval:
    """Evaluates 4.6 relations.  A bare name is the UNROUNDED value (full
    substitution down to the iterate and the constants); round(x, n) is the
    displayed one."""

    def __init__(self, node):
        self.consts = node["constants"]
        self.iterate_sym = node["evaluation_roles"]["iterate"]
        self.mode = node["rounding"]
        self.prec = node["symbol_precision"]
        self.order, self.trees, self.seen = [], {}, set()
        for pair in node["frame_relations"]:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                raise Failure("4.6: frame_relations entry is not [symbol, expression]")
            sym, expr = pair
            tree = P46(tok46(expr)).parse()
            self._names(tree, sym)
            self.order.append(sym); self.trees[sym] = tree; self.seen.add(sym)

    def _names(self, t, sym):
        if t[0] == "name":
            n = t[1]
            if n not in self.consts and n != self.iterate_sym and n not in self.seen:
                raise Failure("4.6: %r resolves to nothing (not a constant, not the "
                              "iterate, not an EARLIER frame_relations symbol) in the "
                              "relation for %r" % (n, sym))
            return
        for sub in t[1:]:
            if isinstance(sub, tuple):
                self._names(sub, sym)

    def value(self, sym, y):
        return self._ev(self.trees[sym], y)

    def _ev(self, t, y):
        k = t[0]
        if k == "lit":
            return t[1]
        if k == "name":
            n = t[1]
            if n == self.iterate_sym:
                return float(y)
            if n in self.consts:
                return float(self.consts[n])
            return self._ev(self.trees[n], y)
        if k == "neg":
            return -self._ev(t[1], y)
        if k == "sqrt":
            v = self._ev(t[1], y)
            if v < 0:
                raise Failure("4.6: sqrt of a negative value")
            return v ** 0.5
        if k == "round":
            sym = t[1][1]
            if sym not in self.prec:
                raise Failure("4.6/8A.9: round(%s) has no symbol_precision entry"
                              % sym)
            return do_round(self._ev(t[1], y), int(self.prec[sym]), self.mode)
        a = self._ev(t[1], y); b = self._ev(t[2], y)
        if k == "+":
            return a + b
        if k == "-":
            return a - b
        if k == "*":
            return a * b
        if k == "/":
            if b == 0:
                raise Failure("4.6: division by zero")
            return a / b
        if k == "**":
            try:
                return a ** b
            except (ValueError, OverflowError, ZeroDivisionError) as e:
                raise Failure("4.6: ** failed (%s)" % e)
        raise Failure("4.6: bad operator")


def _coverage_ok(node, rep):
    """4.6 coverage: frame_relations names every evaluation_symbols member
    except the iterate role, exactly once."""
    isym = node["evaluation_roles"]["iterate"]
    named = []
    for pair in node["frame_relations"]:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            rep.fail("6.1/8A.4", "frame_relations entry is not [symbol, expression]")
            return False
        named.append(pair[0])
    want = [s for s in node["evaluation_symbols"] if s != isym]
    missing = [s for s in want if s not in named]
    dupes = sorted({s for s in named if named.count(s) > 1})
    extra = [s for s in named if s not in want]
    ok = True
    if missing:
        rep.fail("6.1/8A.4", "frame_relations does not cover %s (4.6 coverage: "
                 "every evaluation symbol except the iterate, exactly once)"
                 % sorted(missing))
        ok = False
    if dupes:
        rep.fail("6.1/8A.4", "frame_relations names %s more than once" % dupes)
        ok = False
    if extra:
        rep.fail("6.1/8A.4", "frame_relations names %s, which is not a non-iterate "
                 "evaluation symbol" % sorted(set(extra)))
        ok = False
    return ok


def _frames_ok(node, rep, roles):
    """4.3.9 -- every frame, preamble and evaluation alike, satisfies
    `frame_relations`."""
    if not _coverage_ok(node, rep):
        return
    try:
        fe = FrameEval(node)
    except Failure as e:
        rep.fail("4.6", str(e)); return
    isym = node["evaluation_roles"]["iterate"]
    frames = [("preamble[%d]" % i, f) for i, f in enumerate(node.get("preamble", []))]
    frames += [("element[%d].evaluation" % i, el[roles["evaluation"]])
               for i, el in enumerate(node["elements"])
               if el[roles["evaluation"]] is not None]
    for label, fr in frames:
        y = fr.get(isym)
        for sym in [s for s in node["evaluation_symbols"] if s != isym]:
            if sym not in fr:
                rep.fail("4.3.9", "%s has no stored %r to check" % (label, sym))
                continue
            try:
                raw = fe.value(sym, y)
            except Failure as e:
                rep.fail("4.3.9/8A.4", "%s: %s" % (label, e))
                continue
            if sym not in node["symbol_precision"]:
                rep.fail("8A.9", "%s: no symbol_precision for %r" % (label, sym))
                continue
            got = do_round(raw, int(node["symbol_precision"][sym]), node["rounding"])
            if abs(got - float(fr[sym])) > 1e-9:
                rep.fail("4.3.9/8A.4", "%s: %s recomputes to %r (raw %.10g), stored %r"
                         % (label, sym, got, raw, fr[sym]))


def _iteration_locals(node, rep, roles):
    elements, term = node["elements"], node["termination"]
    eroles = node["evaluation_roles"]
    try:
        tree = Parser(tokenize(node["update_relation"])).parse()
    except Failure as e:
        rep.fail("4.2", str(e))
        return

    # 4.4 preamble_binding (1.3: {frame, frame_role}, required iff preamble)
    if "preamble" in node and "preamble_binding" not in node:
        rep.fail("3.1", "preamble present but preamble_binding absent "
                        "(1.3 requires it iff preamble is present)")
    for role, spec in node.get("preamble_binding", {}).items():
        if role not in roles:
            rep.fail("4.4", "preamble_binding names unknown role %r" % role)
            continue
        if not isinstance(spec, dict) or set(spec) != {"frame", "frame_role"}:
            rep.fail("4.4", "preamble_binding[%r] is not {frame, frame_role}" % role)
            continue
        idx, frole_name = spec["frame"], spec["frame_role"]
        if not isinstance(idx, int) or not (0 <= idx < len(node.get("preamble", []))):
            rep.fail("8A.4", "preamble_binding[%r] frame %r out of range" % (role, idx))
            continue
        if frole_name not in eroles:
            rep.fail("4.4", "preamble_binding[%r] frame_role %r is not an "
                     "evaluation role" % (role, frole_name))
            continue
        got = elements[0][roles[role]]
        want = node["preamble"][idx][eroles[frole_name]]
        if got != want:
            rep.fail("8A.4", "preamble_binding %r -> frame %d.%s: element[0] has "
                     "%r, frame has %r" % (role, idx, frole_name, got, want))

    _frames_ok(node, rep, roles)

    for i, el in enumerate(elements):
        env = {r: el[s] for r, s in roles.items()}
        if el[roles["index"]] != i + 1:
            rep.fail("8A.4", "element[%d] index %r, expected %d"
                     % (i, el[roles["index"]], i + 1))
        try:
            got = eval_expr(tree, env)
            tol = tolerance_for(node, roles["iterate_next"], rep, "8A.9")
            if not close(got, el[roles["iterate_next"]], tol):
                rep.fail("4.3.1", "element[%d] update_relation gives %.10g vs "
                         "iterate_next %r" % (i, got, el[roles["iterate_next"]]))
        except Failure as e:
            rep.fail("4.3.1/4.3.7", "element[%d]: %s" % (i, e))
        try:
            tol = tolerance_for(node, roles["change"], rep, "8A.9")
            want = abs(float(el[roles["iterate_next"]])
                       - float(el[roles["iterate_curr"]]))
            if not close(want, el[roles["change"]], tol):
                rep.fail("4.3.2", "element[%d] change %r vs |next-curr| %.10g"
                         % (i, el[roles["change"]], want))
        except Failure:
            pass
        q = el[roles[term["quantity_role"]]]
        pred = (q < term["tolerance"]) if term["comparison"] == "lt" \
            else (q <= term["tolerance"])
        if bool(el[roles["converged"]]) != pred:
            rep.fail("4.3.3", "element[%d] converged disagrees with the predicate" % i)
        if bool(el[roles["converged"]]) != (i == len(elements) - 1):
            rep.fail("8A.4", "element[%d] converged=%r; only the last may be true"
                     % (i, el[roles["converged"]]))
        ev = el[roles["evaluation"]]
        if (ev is None) != bool(el[roles["converged"]]):
            rep.fail("8A.4", "element[%d] evaluation nullity disagrees with converged" % i)
        if ev is not None and ev[eroles["iterate"]] != el[roles["iterate_next"]]:
            rep.fail("8A.4", "element[%d] evaluation bound to the wrong iterate" % i)


def _decision_locals(node, rep, roles):
    elements, term, sel = node["elements"], node["termination"], node["selection"]
    tr_r, ca_r = node["trial_roles"], node["candidate_roles"]
    universe, prec = list(term["universe"]), term["precedences"]
    im = term["item_measures"]
    filt_role, crit_role = sel["filter"], sel["criterion"]
    for r in (filt_role, crit_role):
        if r not in ca_r:
            rep.fail("5.2", "selection names %r, which is not a candidate role" % r)
            return
    try:
        frel = parse_filter_relation(sel["filter_relation"])
    except Failure as e:
        rep.fail("5.2", str(e))
        return

    def C(en, role):
        return en[ca_r[role]]

    def T(tr, role):
        return tr[tr_r[role]]

    # 5.3 global
    for i, el in enumerate(elements):
        if el[roles["index"]] != i + 1:
            rep.fail("8A.5", "element[%d] index %r, expected %d"
                     % (i, el[roles["index"]], i + 1))
    seen = set()
    for el in elements:
        for it in el[roles["committed"]]:
            if it in seen:
                rep.fail("8A.5", "%r committed in more than one element" % it)
            seen.add(it)
    if seen != set(universe):
        rep.fail("8A.5", "union of committed %s != universe %s"
                 % (sorted(seen), sorted(set(universe))))
    if term["max_elements"] < len(universe):
        rep.fail("5.3", "max_elements < len(universe)")
    for i, el in enumerate(elements):
        if float(el[roles["budget_total"]]) != float(term["budget"]):
            rep.fail("8A.5", "element[%d] budget_total %r != termination.budget %r"
                     % (i, el[roles["budget_total"]], term["budget"]))

    committed = []
    for i, el in enumerate(elements):
        total = el[roles["budget_total"]]
        if el[roles["budget_initial"]] != total:
            rep.fail("5.4.1", "element[%d] budget_initial != budget_total" % i)
        budget = el[roles["budget_initial"]]
        chosen_seq, trials = [], el[roles["trials"]]
        for j, tr in enumerate(trials):
            cands, chosen = T(tr, "candidates"), T(tr, "chosen")
            by_item = {C(en, "item"): en for en in cands}
            if bool(T(tr, "closes")) != (chosen is None):
                rep.fail("8A.5", "element[%d] trial[%d] closes disagrees with chosen"
                         % (i, j))
            # 5.4.9 (new in 1.2)
            if bool(T(tr, "closes")) and j != len(trials) - 1:
                rep.fail("8A.5", "element[%d] trial[%d] closes but is not last" % (i, j))
            # 5.4.3b precedences, both directions
            comm = set(committed) | set(chosen_seq)
            for en in cands:
                for p in prec.get(C(en, "item"), []):
                    if p not in comm:
                        rep.fail("8A.5", "element[%d] trial[%d] candidate %r has "
                                 "uncommitted prerequisite %r"
                                 % (i, j, C(en, "item"), p))
            for u in universe:
                if u not in comm and all(p in comm for p in prec.get(u, [])) \
                        and u not in by_item:
                    rep.fail("8A.5", "element[%d] trial[%d] omits candidate %r whose "
                             "precedences are satisfied" % (i, j, u))
            # 5.4.3c -- recompute `admissible` for EVERY candidate (new in 1.2)
            for en in cands:
                env = {r: en[s] for r, s in ca_r.items()}
                env["budget_before"] = budget
                try:
                    want = eval_filter_relation(frel, env)
                except Failure as e:
                    rep.fail("5.4.3c", "element[%d] trial[%d] candidate %r: %s"
                             % (i, j, C(en, "item"), e))
                    continue
                if bool(C(en, "admissible")) != want:
                    rep.fail("5.4.3c/8A.5", "element[%d] trial[%d] candidate %r has "
                             "admissible=%r but filter_relation %r gives %r"
                             % (i, j, C(en, "item"), C(en, "admissible"),
                                sel["filter_relation"], want))
            # 5.4.3d -- measure declared once on the node (new in 1.3)
            for en in cands:
                it_, mv = C(en, "item"), C(en, "measure")
                if it_ not in im:
                    rep.fail("5.4.3d/8A.5", "element[%d] trial[%d] candidate %r not "
                             "in item_measures" % (i, j, it_))
                elif float(mv) != float(im[it_]):
                    rep.fail("5.4.3d/8A.5", "element[%d] trial[%d] candidate %r has "
                             "measure %r, item_measures declares %r"
                             % (i, j, it_, mv, im[it_]))
            # 5.4.4 selection replay
            adm = [en for en in cands if bool(en[ca_r[filt_role]])]
            if not adm:
                if chosen is not None:
                    rep.fail("8A.5", "element[%d] trial[%d] chose %r with no "
                             "admissible candidate" % (i, j, chosen))
            else:
                best = (max if sel["objective"] == "max" else min)(
                    float(en[ca_r[crit_role]]) for en in adm)
                tied = [en for en in adm if float(en[ca_r[crit_role]]) == best]
                want = min(tied, key=lambda en: universe.index(C(en, "item"))) \
                    if sel["tie_break"] == "universe_order" else tied[0]
                if chosen != C(want, "item"):
                    rep.fail("8A.5", "element[%d] trial[%d] chose %r; selection "
                             "requires %r" % (i, j, chosen, C(want, "item")))
            if chosen is not None:
                if chosen not in by_item:
                    rep.fail("5.4.3", "element[%d] trial[%d] chosen %r not a candidate"
                             % (i, j, chosen))
                elif not bool(by_item[chosen][ca_r[filt_role]]):
                    rep.fail("5.4.3", "element[%d] trial[%d] chosen %r is inadmissible"
                             % (i, j, chosen))
            spent = float(by_item[chosen][ca_r[crit_role]]) if chosen in by_item else 0.0
            if not (chosen is not None and chosen not in by_item):
                if not close(float(budget) - spent, T(tr, "budget_after"), 1e-9):
                    rep.fail("5.4.5", "element[%d] trial[%d] budget_after wrong" % (i, j))
            budget = T(tr, "budget_after")
            if chosen is not None:
                chosen_seq.append(chosen)
        if list(el[roles["committed"]]) != chosen_seq:
            rep.fail("5.4.6", "element[%d] committed != ordered chosen values" % i)
        measures = {}
        for tr in trials:
            for en in T(tr, "candidates"):
                measures[C(en, "item")] = float(C(en, "measure"))
        missing = [t for t in el[roles["committed"]] if t not in measures]
        if missing:
            rep.fail("5.4.7", "element[%d] cannot price %r" % (i, missing))
        elif not close(float(total) - sum(measures[t] for t in el[roles["committed"]]),
                       el[roles["budget_final"]], 1e-9):
            rep.fail("5.4.7", "element[%d] budget_final wrong" % i)
        committed.extend(el[roles["committed"]])


def _termination(node, rep, roles, ntype, mode):
    term, elements = node["termination"], node["elements"]
    if not term["satisfied"]:
        if mode == "gold":
            rep.fail("8A.1", "termination.satisfied is false in gold mode")
        else:
            rep.note_unchecked("3.4/6.0: satisfied is false; legitimate in "
                               "candidate mode, reported")
    if ntype == "iteration":
        q = term["quantity_role"]
        if q not in roles:
            rep.fail("3.4", "quantity_role %r is not a role" % q)
            return
        hold = [(el[roles[q]] < term["tolerance"]) if term["comparison"] == "lt"
                else (el[roles[q]] <= term["tolerance"]) for el in elements]
    else:
        if term["scope"] != "cumulative":
            rep.fail("3.7", "unsupported scope %r" % term["scope"])
            return
        # ASSUMPTION[5] SURVIVES INTO 1.2: §3.7 still reads "the union over ALL
        # elements", under which step 7's "and on no earlier one" is undecidable.
        # Still read as the union over the PREFIX ending at element k.  See F7.
        if term["accumulator_role"] not in roles:
            rep.fail("3.4", "accumulator_role %r is not an element role"
                     % term["accumulator_role"])
            return
        acc = roles[term["accumulator_role"]]
        uni, run, hold = set(term["universe"]), set(), []
        for el in elements:
            run |= set(el[acc])
            hold.append(run == uni)
    if not hold[-1]:
        rep.fail("6.7", "termination predicate does not hold on the last element")
    if any(hold[:-1]):
        rep.fail("6.7", "termination predicate holds before the last element")


def _result(node, rep, roles):
    res = node["result"]
    if res["from"] == "len(elements)":
        got = float(len(node["elements"]))
    elif res["from"].startswith("last."):
        role = res["from"][5:]
        if role not in roles:
            rep.fail("3.6", "result.from names unknown role %r" % role)
            return
        got = do_round(node["elements"][-1][roles[role]], res["dp"], node["rounding"])
    else:
        rep.fail("3.6", "result.from %r not in the closed set" % res["from"])
        return
    if abs(got - float(res["value"])) > 1e-9:
        rep.fail("6.8", "result.value %r != derivation %r" % (res["value"], got))


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    mode = "candidate" if "--mode" in argv and "candidate" in argv else "gold"
    rows = json.load(open(args[0], encoding="utf-8"))
    npass = nfail = nunchecked = 0
    per = {}
    for r in rows:
        rep = verify_node(r["trace_nodes"],
                          "%s seed %s" % (r.get("template_id"), r.get("seed")), mode)
        d = per.setdefault(r.get("template_id"), [0, 0, set()])
        if rep.passed:
            npass += 1
            d[0] += 1
        else:
            nfail += 1
            d[1] += 1
            d[2].update(c for c, _ in rep.failures)
            if d[1] <= 2:
                print("FAIL %s" % rep.label)
                for c, m in rep.failures[:4]:
                    print("      [%s] %s" % (c, m))
        if rep.unchecked:
            nunchecked += 1
    print()
    for t, (p, f, cl) in sorted(per.items()):
        print("%-40s pass %2d  fail %2d  %s"
              % (t, p, f, ("clauses " + ",".join(sorted(cl))) if cl else ""))
    print("TOTAL %d rows: %d pass, %d fail; %d carried unchecked notes"
          % (len(rows), npass, nfail, nunchecked))
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
