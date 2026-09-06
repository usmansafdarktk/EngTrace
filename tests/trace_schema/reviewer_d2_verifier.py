"""Reviewer D2 -- independent verifier for D3.3 schema version 1.1.

Written from `docs/re-implementation-sep/phase3_node_types.md` and the
conformance corpus ALONE. No template source, no extractor, and no round-1
verifier was read. Every place the spec forced an assumption carries an
`ASSUMPTION[n]` comment; the numbering matches the review report
`docs/re-implementation-sep/reviews/phase3_reviewer_d2_schema.md`.

Usage:
    python -m tests.trace_schema.reviewer_d2_verifier <traces.json>
    python -m tests.trace_schema.reviewer_d2_verifier <traces.json> --lenient-roles
    python -m tests.trace_schema.reviewer_d2_verifier <traces.json> --check-filter
"""

from __future__ import annotations

import json
import sys
from decimal import Decimal, ROUND_HALF_UP

KNOWN_MAJOR = "1"

# Section 4.1 -- the mandatory type-level roles of an `iteration` node.
ITERATION_ROLES = (
    "index", "iterate_prev", "iterate_curr", "iterate_next",
    "residual_prev", "residual_curr", "change", "converged", "evaluation",
)
# Section 4.1 prose -- the two roles `evaluation_roles` must bind.
ITERATION_EVAL_ROLES = ("iterate", "residual")

# Section 5.1 -- the `decision` tables are written as literal SYMBOL names,
# never as roles.  ASSUMPTION[6]: no role table for `decision` exists anywhere
# in the spec, so these names are hardcoded.  This is the ONLY place this
# verifier reaches for a literal symbol name, and it does so because the spec
# left no alternative.  See finding F2.
DECISION_ELEMENT_SYMBOLS = (
    "station_id", "capacity", "remaining_initial",
    "trials", "assigned", "remaining_final",
)
DECISION_TRIAL_SYMBOLS = ("eligible", "chosen", "remaining_after", "closes")


class Failure(Exception):
    pass


class Report:
    def __init__(self, label):
        self.label = label
        self.failures = []      # (clause, message)
        self.unchecked = []     # section 6.5 / 8.4 channel

    def fail(self, clause, message):
        self.failures.append((clause, message))

    def note_unchecked(self, what):
        self.unchecked.append(what)

    @property
    def passed(self):
        return not self.failures

    @property
    def status(self):
        if self.failures:
            return "FAIL"
        # Section 8.4: an unchecked relation must never silently count as a
        # satisfied one, so it is surfaced in the verdict itself.
        return "PASS+UNCHECKED" if self.unchecked else "PASS"


# --------------------------------------------------------------------------
# Section 3.1 `rounding`; section 6 "every rounding uses `rounding`".
# --------------------------------------------------------------------------

def do_round(value, dp, mode):
    if mode != "decimal-half-up":
        raise Failure("6: unknown rounding mode %r" % (mode,))
    q = Decimal(1).scaleb(-int(dp))
    return float(Decimal(repr(float(value))).quantize(q, rounding=ROUND_HALF_UP))


def tolerance_for(node, symbol, rep, clause):
    """Section 7.1 -- half a unit in the last place gold displays the symbol.

    Taken from `symbol_precision` only; never inferred from the data.
    """
    sp = node.get("symbol_precision")
    if not isinstance(sp, dict) or symbol not in sp:
        rep.fail(clause, "7.1: no symbol_precision entry for %r; tolerance cannot "
                         "be derived without inferring it" % (symbol,))
        raise Failure(symbol)
    return 0.5 * (10.0 ** -int(sp[symbol]))


def close(a, b, tol):
    # ASSUMPTION[1]: section 7.1 says "agrees to half a unit in the last place"
    # but does not say whether the boundary itself matches.  Read as inclusive,
    # because section 4.5 reports a worst gold residual of 4.985e-05 "against a
    # 5.0e-05 tolerance" and calls that conforming.
    return abs(float(a) - float(b)) <= tol + 1e-12


# --------------------------------------------------------------------------
# Section 4.2 -- `update_relation` grammar: role identifiers, decimal
# literals, binary + - * /, unary -, parentheses.  Nothing else.
# Recursive descent; never `eval`.
# --------------------------------------------------------------------------

def tokenize(src):
    toks, i = [], 0
    while i < len(src):
        c = src[i]
        if c in " \t\n\r":
            i += 1
        elif c in "+-*/()":
            toks.append((c, c))
            i += 1
        elif c.isdigit() or (c == "." and i + 1 < len(src) and src[i + 1].isdigit()):
            j = i
            while j < len(src) and (src[j].isdigit() or src[j] == "."):
                j += 1
            toks.append(("num", src[i:j]))
            i = j
        elif c.isalpha() or c == "_":
            j = i
            while j < len(src) and (src[j].isalnum() or src[j] == "_"):
                j += 1
            toks.append(("name", src[i:j]))
            i = j
        else:
            raise Failure("4.2: character %r is not in the update_relation grammar" % c)
    return toks


class Parser:
    def __init__(self, toks):
        self.toks, self.i = toks, 0

    def peek(self):
        return self.toks[self.i][0] if self.i < len(self.toks) else None

    def take(self):
        t = self.toks[self.i]
        self.i += 1
        return t

    def parse(self):
        node = self.expr()
        if self.i != len(self.toks):
            raise Failure("4.2: trailing tokens in update_relation")
        return node

    def expr(self):
        node = self.term()
        while self.peek() in ("+", "-"):
            op = self.take()[0]
            node = (op, node, self.term())
        return node

    def term(self):
        node = self.unary()
        while self.peek() in ("*", "/"):
            op = self.take()[0]
            node = (op, node, self.unary())
        return node

    def unary(self):
        if self.peek() == "-":
            self.take()
            return ("neg", self.unary())
        return self.atom()

    def atom(self):
        if self.peek() is None:
            raise Failure("4.2: unexpected end of update_relation")
        kind, text = self.take()
        if kind == "(":
            node = self.expr()
            if self.peek() != ")":
                raise Failure("4.2: unbalanced parenthesis in update_relation")
            self.take()
            return node
        if kind == "num":
            return ("lit", float(text))
        if kind == "name":
            return ("role", text)
        raise Failure("4.2: token %r is not in the update_relation grammar" % (text,))


def eval_expr(node, env):
    kind = node[0]
    if kind == "lit":
        return node[1]
    if kind == "role":
        if node[1] not in env:
            # Section 4.2: "no names that are not roles of this node".
            raise Failure("4.2: %r is not a role of this node" % (node[1],))
        v = env[node[1]]
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise Failure("4.2: role %r is not numeric" % (node[1],))
        return float(v)
    if kind == "neg":
        return -eval_expr(node[1], env)
    a = eval_expr(node[1], env)
    b = eval_expr(node[2], env)
    if kind == "+":
        return a + b
    if kind == "-":
        return a - b
    if kind == "*":
        return a * b
    if kind == "/":
        if b == 0:
            # Sections 4.3.7 / 8.5 -- a zero update denominator.
            raise Failure("4.3.7: zero denominator in update_relation")
        return a / b
    raise Failure("4.2: bad operator")


# --------------------------------------------------------------------------
# Section 3.5 -- `carry` path grammar and evaluation.
# --------------------------------------------------------------------------

def parse_path(path, element_symbols):
    if not isinstance(path, str) or not path:
        raise Failure("3.5: carry path must be a non-empty string")
    parts = path.split(".")
    for p in parts:
        if not p or not (p[0].isalpha() or p[0] == "_"):
            raise Failure("3.5: %r is not a well-formed path" % (path,))
        if not all(ch.isalnum() or ch == "_" for ch in p):
            raise Failure("3.5: %r is not a well-formed path" % (path,))
    if parts[0] not in element_symbols:
        raise Failure("3.5: path root %r is not in element_symbols" % (parts[0],))
    return parts


def select_path(parts, element):
    cur = element
    for p in parts:
        if not isinstance(cur, dict):
            # ASSUMPTION[2]: the path grammar does not say what happens when a
            # path traverses a null (`evaluation` is null on the terminating
            # element).  Treated as a carry failure.
            raise Failure("3.5: path segment %r indexes into a non-object" % p)
        if p not in cur:
            raise Failure("3.5: path segment %r absent" % p)
        cur = cur[p]
    return cur


# --------------------------------------------------------------------------
# Section 6 -- the normative algorithm.
# --------------------------------------------------------------------------

def keyset_check(rep, clause, what, obj, symbols):
    if not isinstance(obj, dict):
        rep.fail(clause, "%s is not an object" % what)
        return
    have, want = set(obj), set(symbols)
    if have != want:
        rep.fail(clause, "%s key set %s != declared %s (missing %s, extra %s)"
                 % (what, sorted(have), sorted(want),
                    sorted(want - have), sorted(have - want)))


def verify_node(node, label="node", strict_roles=True, check_filter=False):
    rep = Report(label)
    try:
        _verify(node, rep, strict_roles, check_filter)
    except Failure as e:
        if not rep.failures:
            rep.fail("6", str(e))
    except (KeyError, TypeError, IndexError, ValueError) as e:
        rep.fail("6", "%s while verifying: %s" % (type(e).__name__, e))
    return rep


def _verify(node, rep, strict_roles, check_filter):
    if not isinstance(node, dict):
        rep.fail("3.0", "node is not a JSON object")
        return

    # --- 6.1 version and shape -------------------------------------------
    sv = node.get("schema_version")
    if not isinstance(sv, str):
        rep.fail("8.9", "schema_version missing or not a string")
        return
    if sv.split(".")[0] != KNOWN_MAJOR:
        rep.fail("8.9", "unknown major schema version %r" % sv)
        return

    ntype = node.get("node_type")
    if ntype not in ("iteration", "decision"):
        rep.fail("6.1", "unrecognised node_type %r" % (ntype,))
        return

    base_required = ["schema_version", "node_type", "node_id", "cardinality",
                     "roles", "termination", "rounding", "symbol_precision",
                     "element_symbols", "elements", "carry", "carry_notes",
                     "result"]
    cond = {"iteration": ["evaluation_symbols", "evaluation_roles", "update_relation"],
            "decision": ["trial_symbols", "eligible_symbols", "selection",
                         "capacity_constant"]}[ntype]
    required = base_required + cond
    if not strict_roles:
        # See finding F1: all 40 gold `decision` nodes carry no `roles` key,
        # although section 3.1 marks it required for every node.  This switch
        # exists only so the remaining checks can be exercised at all.
        required = [f for f in required if f != "roles"]
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
    if not isinstance(term, dict):
        rep.fail("3.4", "termination is not an object")
        return
    if term.get("kind") != fixed[1]:
        rep.fail("6.1", "termination.kind %r; %s fixes %r"
                 % (term.get("kind"), ntype, fixed[1]))

    # Section 8.7
    ab = node["cardinality"] == "answer_bearing"
    if ab and "cardinality_symbol" not in node:
        rep.fail("8.7", "answer_bearing with no cardinality_symbol")
    if not ab and "cardinality_symbol" in node:
        rep.fail("8.7", "incidental with a cardinality_symbol")

    if ("preamble" in node) != ("preamble_symbols" in node):
        rep.fail("3.1", "preamble and preamble_symbols must appear together")

    for f in ("satisfied", "max_elements"):
        if f not in term:
            rep.fail("3.4", "termination.%s absent" % f)
    tcond = {"convergence": ("quantity_role", "comparison", "tolerance"),
             "exhaustion": ("scope", "accumulator_symbol", "universe", "precedences")}
    for f in tcond.get(term.get("kind"), ()):
        if f not in term:
            rep.fail("3.4", "termination.%s absent (required for kind %r)"
                     % (f, term.get("kind")))
    if rep.failures:
        return

    esyms = node["element_symbols"]
    elements = node["elements"]
    if not isinstance(elements, list):
        rep.fail("3.1", "elements is not a list")
        return

    # --- 6.2 roles --------------------------------------------------------
    roles = node.get("roles", {})
    if strict_roles or "roles" in node:
        if not isinstance(roles, dict):
            rep.fail("6.2", "roles is not an object")
            return
        mandatory = ITERATION_ROLES if ntype == "iteration" else DECISION_ELEMENT_SYMBOLS
        for r in mandatory:
            if r not in roles:
                rep.fail("8.8", "mandatory role %r missing from roles" % r)
        for r, sym in roles.items():
            if sym not in esyms:
                rep.fail("8.8", "roles[%r] = %r is not in element_symbols" % (r, sym))
    if ntype == "decision" and not roles:
        # Only reachable under --lenient-roles.
        roles = {s: s for s in DECISION_ELEMENT_SYMBOLS}

    if ntype == "iteration":
        eroles = node["evaluation_roles"]
        esym_eval = node["evaluation_symbols"]
        for r in ITERATION_EVAL_ROLES:
            if r not in eroles:
                rep.fail("8.8", "mandatory evaluation role %r missing" % r)
        for r, sym in eroles.items():
            if sym not in esym_eval:
                rep.fail("8.8", "evaluation_roles[%r] = %r not in evaluation_symbols"
                         % (r, sym))
    if rep.failures:
        return

    def R(el, role):
        return el[roles[role]]

    # --- 6.3 homogeneity --------------------------------------------------
    for i, el in enumerate(elements):
        keyset_check(rep, "8.2", "element[%d]" % i, el, esyms)
    if "preamble" in node:
        for i, fr in enumerate(node["preamble"]):
            keyset_check(rep, "8.2", "preamble[%d]" % i, fr, node["preamble_symbols"])
        if ntype == "iteration" and \
                set(node["preamble_symbols"]) != set(node["evaluation_symbols"]):
            # Section 4.4 states this as fact, not as a MUST.  ASSUMPTION[3]:
            # read as normative.
            rep.fail("4.4", "preamble_symbols != evaluation_symbols")
    if rep.failures:
        return
    if ntype == "iteration":
        for i, el in enumerate(elements):
            ev = R(el, "evaluation")
            if ev is not None:
                keyset_check(rep, "8.2", "element[%d].evaluation" % i, ev,
                             node["evaluation_symbols"])
    else:
        for i, el in enumerate(elements):
            trials = el[roles["trials"]]
            if not isinstance(trials, list) or not trials:
                rep.fail("5.1", "element[%d].trials empty or not a list" % i)
                continue
            for j, tr in enumerate(trials):
                keyset_check(rep, "8.2", "element[%d].trial[%d]" % (i, j),
                             tr, node["trial_symbols"])
                if isinstance(tr, dict) and isinstance(tr.get("eligible"), list):
                    for m, en in enumerate(tr["eligible"]):
                        keyset_check(rep, "8.2",
                                     "element[%d].trial[%d].eligible[%d]" % (i, j, m),
                                     en, node["eligible_symbols"])
    if rep.failures:
        return

    # --- 6.4 budget -------------------------------------------------------
    if len(elements) < 1:
        rep.fail("6.4", "elements is empty; section 6.4 requires >= 1")
    if len(elements) > term["max_elements"]:
        rep.fail("6.4", "len(elements)=%d exceeds max_elements=%d"
                 % (len(elements), term["max_elements"]))
    if rep.failures:
        return

    # --- 6.5 carry --------------------------------------------------------
    carry = node["carry"]
    for target, path in carry.items():
        if target not in esyms:
            rep.fail("8.3", "carry target %r is not an element symbol" % target)
            continue
        try:
            parts = parse_path(path, esyms)
        except Failure as e:
            rep.fail("8.3", str(e))
            continue
        for k in range(len(elements) - 1):
            try:
                got = select_path(parts, elements[k])
            except Failure as e:
                rep.fail("8.3", "element[%d]: %s" % (k, e))
                continue
            want = elements[k + 1][target]
            if isinstance(got, bool) or not isinstance(got, (int, float)):
                ok = got == want
            else:
                ok = isinstance(want, (int, float)) and float(got) == float(want)
            if not ok:
                rep.fail("8.3", "carry %r <- %r broken between element %d and %d "
                                "(%r vs %r)" % (target, path, k, k + 1, got, want))
    for k, v in node["carry_notes"].items():
        rep.note_unchecked("carry_notes[%r] = %r (section 8.4: reported unchecked, "
                           "never counted as satisfied)" % (k, v))

    # --- 6.6 local invariants --------------------------------------------
    if ntype == "iteration":
        _iteration_locals(node, rep, roles)
    else:
        _decision_locals(node, rep, roles, check_filter)
    if rep.failures:
        return

    # --- 6.7 termination --------------------------------------------------
    _termination(node, rep, roles, ntype)
    if rep.failures:
        return

    # --- 6.8 result -------------------------------------------------------
    _result(node, rep, roles)


def _iteration_locals(node, rep, roles):
    elements = node["elements"]
    term = node["termination"]
    try:
        tree = Parser(tokenize(node["update_relation"])).parse()
    except Failure as e:
        rep.fail("4.2", str(e))
        return

    for i, el in enumerate(elements):
        env = {r: el[sym] for r, sym in roles.items()}
        # 4.3.8 -- index is 1..len(elements) in order.
        if el[roles["index"]] != i + 1:
            rep.fail("4.3.8/8.5", "element[%d] index is %r, expected %d"
                     % (i, el[roles["index"]], i + 1))
        # 4.3.1 -- update_relation reproduces iterate_next.
        try:
            got = eval_expr(tree, env)
            tol = tolerance_for(node, roles["iterate_next"], rep, "4.3.1")
            if not close(got, el[roles["iterate_next"]], tol):
                rep.fail("4.3.1", "element[%d] update_relation gives %.10g, "
                                  "iterate_next is %r (tol %g)"
                         % (i, got, el[roles["iterate_next"]], tol))
        except Failure as e:
            rep.fail("4.3.1/4.3.7", "element[%d]: %s" % (i, e))
        # 4.3.2 -- change is |next - curr|.
        try:
            tol = tolerance_for(node, roles["change"], rep, "4.3.2")
            want = abs(float(el[roles["iterate_next"]]) - float(el[roles["iterate_curr"]]))
            if not close(want, el[roles["change"]], tol):
                rep.fail("4.3.2", "element[%d] change is %r, |next-curr| is %.10g"
                         % (i, el[roles["change"]], want))
        except Failure:
            pass
        # 4.3.3 -- converged is the termination test's value on this element.
        q = el[roles[term["quantity_role"]]]
        pred = (q < term["tolerance"]) if term["comparison"] == "lt" \
            else (q <= term["tolerance"])
        if bool(el[roles["converged"]]) != pred:
            rep.fail("4.3.3", "element[%d] converged=%r but predicate is %r"
                     % (i, el[roles["converged"]], pred))
        # 4.3.4 -- converged true on the last element only.
        if bool(el[roles["converged"]]) != (i == len(elements) - 1):
            rep.fail("4.3.4/8.5", "element[%d] converged=%r; only the last element "
                     "may be true" % (i, el[roles["converged"]]))
        # 4.3.5 -- evaluation is null iff converged.
        ev = el[roles["evaluation"]]
        if (ev is None) != bool(el[roles["converged"]]):
            rep.fail("4.3.5/8.5", "element[%d] evaluation %s but converged=%r"
                     % (i, "null" if ev is None else "non-null",
                        el[roles["converged"]]))
        # 4.3.6 -- the frame is bound to the iterate that produced it.
        if ev is not None:
            ir = node["evaluation_roles"]["iterate"]
            if ev[ir] != el[roles["iterate_next"]]:
                rep.fail("4.3.6/8.5", "element[%d] evaluation iterate %r != "
                         "iterate_next %r" % (i, ev[ir], el[roles["iterate_next"]]))


def _identity_symbol(node, rep):
    """The `eligible` entry's identity key.

    ASSUMPTION[7]: section 5.2 declares `filter` and `criterion` as data but
    NEVER declares which `eligible_symbols` member carries the candidate's
    name -- yet `chosen` is compared against it and against `universe`.
    Derived here as the unique remaining member of `eligible_symbols`, which
    works only because that list has exactly three entries.  See finding F5.
    """
    sel = node["selection"]
    rest = [s for s in node["eligible_symbols"]
            if s not in (sel["filter"], sel["criterion"])]
    if len(rest) != 1:
        rep.fail("5.2", "cannot determine the eligible identity symbol from "
                        "eligible_symbols %r minus filter/criterion (candidates %r)"
                 % (node["eligible_symbols"], rest))
        raise Failure("identity symbol undeterminable")
    return rest[0]


def _decision_locals(node, rep, roles, check_filter):
    elements = node["elements"]
    term = node["termination"]
    sel = node["selection"]
    universe = list(term["universe"])
    prec = term["precedences"]
    filt, crit = sel["filter"], sel["criterion"]
    ident = _identity_symbol(node, rep)

    # --- 5.3 global -------------------------------------------------------
    for i, el in enumerate(elements):
        if el[roles["station_id"]] != i + 1:
            rep.fail("5.3/8.6", "element[%d] station_id %r, expected %d"
                     % (i, el[roles["station_id"]], i + 1))
    seen = set()
    for el in elements:
        for item in el[roles["assigned"]]:
            if item in seen:
                rep.fail("5.3/8.6", "%r assigned in more than one element" % item)
            seen.add(item)
    if seen != set(universe):
        rep.fail("5.3/8.6", "union of assigned %s != universe %s"
                 % (sorted(seen), sorted(set(universe))))
    if term["max_elements"] < len(universe):
        rep.fail("5.3", "max_elements=%d < len(universe)=%d"
                 % (term["max_elements"], len(universe)))
    if node["capacity_constant"]:
        caps = {el[roles["capacity"]] for el in elements}
        if len(caps) > 1:
            rep.fail("5.4.8/8.6", "capacity_constant but capacities are %s"
                     % sorted(caps))

    # --- 5.4 local --------------------------------------------------------
    committed = []
    for i, el in enumerate(elements):
        cap = el[roles["capacity"]]
        if el[roles["remaining_initial"]] != cap:
            rep.fail("5.4.1", "element[%d] remaining_initial %r != capacity %r"
                     % (i, el[roles["remaining_initial"]], cap))
        budget = el[roles["remaining_initial"]]
        chosen_seq = []
        trials = el[roles["trials"]]
        for j, tr in enumerate(trials):
            elig = tr["eligible"]
            chosen = tr["chosen"]
            by_id = {en[ident]: en for en in elig}
            # 5.1 -- closes exactly when nothing was chosen.
            if bool(tr["closes"]) != (chosen is None):
                rep.fail("5.1/8.6", "element[%d] trial[%d] closes=%r but chosen=%r"
                         % (i, j, tr["closes"], chosen))
            # 5.4.3b -- precedences, both directions.
            comm = set(committed) | set(chosen_seq)
            for en in elig:
                for p in prec.get(en[ident], []):
                    if p not in comm:
                        rep.fail("5.4.3b/8.6", "element[%d] trial[%d] eligible %r has "
                                 "uncommitted prerequisite %r" % (i, j, en[ident], p))
            for u in universe:
                if u in comm:
                    continue
                if all(p in comm for p in prec.get(u, [])) and u not in by_id:
                    rep.fail("5.4.3b/8.6", "element[%d] trial[%d] omits eligible item "
                             "%r whose precedences are satisfied" % (i, j, u))
            # 5.4.4 -- replay `selection` over the FILTERED candidates.
            cands = [en for en in elig if bool(en[filt])]
            if check_filter:
                # NOT NORMATIVE -- the spec never says how `filter` is computed.
                # See finding F4.
                for en in elig:
                    want = float(en[crit]) <= float(budget) + 1e-12
                    if bool(en[filt]) != want:
                        rep.fail("F4(non-normative)", "element[%d] trial[%d] %r has "
                                 "%s=%r but criterion %r against budget %r"
                                 % (i, j, en[ident], filt, en[filt], en[crit], budget))
            if not cands:
                if chosen is not None:
                    rep.fail("5.4.4/8.6", "element[%d] trial[%d] chose %r with no "
                             "selectable candidate" % (i, j, chosen))
            else:
                best = (max if sel["objective"] == "max" else min)(
                    float(en[crit]) for en in cands)
                tied = [en for en in cands if float(en[crit]) == best]
                if sel["tie_break"] == "universe_order":
                    want = min(tied, key=lambda en: universe.index(en[ident]))
                elif sel["tie_break"] == "first":
                    want = tied[0]
                else:
                    rep.fail("5.2", "unknown tie_break %r" % sel["tie_break"])
                    want = tied[0]
                if chosen != want[ident]:
                    rep.fail("5.4.4/8.6", "element[%d] trial[%d] chose %r; selection "
                             "requires %r" % (i, j, chosen, want[ident]))
            # 5.4.3 -- the chosen candidate is eligible and passes `filter`.
            if chosen is not None:
                if chosen not in by_id:
                    rep.fail("5.4.3", "element[%d] trial[%d] chosen %r is not in "
                             "eligible" % (i, j, chosen))
                elif not bool(by_id[chosen][filt]):
                    rep.fail("5.4.3", "element[%d] trial[%d] chosen %r has %s false"
                             % (i, j, chosen, filt))
            # 5.4.5 / 5.4.2 -- budget arithmetic, threaded across trials.
            spent = float(by_id[chosen][crit]) if chosen in by_id else 0.0
            if chosen is not None and chosen not in by_id:
                pass  # already reported by 5.4.3
            elif not close(float(budget) - spent, tr["remaining_after"], 1e-9):
                rep.fail("5.4.5", "element[%d] trial[%d] remaining_after %r != %r - %r"
                         % (i, j, tr["remaining_after"], budget, spent))
            budget = tr["remaining_after"]
            if chosen is not None:
                chosen_seq.append(chosen)
            if tr["closes"] and j != len(trials) - 1:
                # ASSUMPTION[4]: the spec never says a closing trial must be
                # the last one, but section 5's "opened, filled, closed" and
                # the budget chain leave no other reading.
                rep.fail("5.1", "element[%d] trial[%d] closes but is not the last"
                         % (i, j))
        # 5.4.6
        if list(el[roles["assigned"]]) != chosen_seq:
            rep.fail("5.4.6", "element[%d] assigned %r != ordered chosen values %r"
                     % (i, el[roles["assigned"]], chosen_seq))
        # 5.4.7
        durs = {}
        for tr in trials:
            for en in tr["eligible"]:
                durs[en[ident]] = float(en[crit])
        missing = [t for t in el[roles["assigned"]] if t not in durs]
        if missing:
            rep.fail("5.4.7", "element[%d] cannot price assigned items %r" % (i, missing))
        else:
            tot = sum(durs[t] for t in el[roles["assigned"]])
            if not close(float(cap) - tot, el[roles["remaining_final"]], 1e-9):
                rep.fail("5.4.7", "element[%d] remaining_final %r != capacity %r - %r"
                         % (i, el[roles["remaining_final"]], cap, tot))
        committed.extend(el[roles["assigned"]])


def _termination(node, rep, roles, ntype):
    term = node["termination"]
    elements = node["elements"]
    # 8.1 -- a gold node with satisfied:false is malformed.
    if not term["satisfied"]:
        rep.fail("8.1", "termination.satisfied is false on a gold node")
    if ntype == "iteration":
        q = term["quantity_role"]
        if q not in roles:
            rep.fail("3.4", "termination.quantity_role %r is not a role" % q)
            return
        vals = [el[roles[q]] for el in elements]
        hold = [(v < term["tolerance"]) if term["comparison"] == "lt"
                else (v <= term["tolerance"]) for v in vals]
    else:
        # ASSUMPTION[5]: section 3.7 defines `cumulative` as "the union over
        # ALL elements", which can never be false at an earlier element, so
        # section 6.7's "and on no earlier one" would be undecidable.  Read as
        # the union over the PREFIX ending at element k.  See finding F3.
        if term["scope"] != "cumulative":
            rep.fail("3.7", "unsupported termination.scope %r" % term["scope"])
            return
        acc = roles.get(term["accumulator_symbol"], term["accumulator_symbol"])
        uni, run, hold = set(term["universe"]), set(), []
        for el in elements:
            run |= set(el[acc])
            hold.append(run == uni)
    if not hold[-1]:
        rep.fail("6.7", "termination predicate does not hold on the last element")
    if any(hold[:-1]):
        rep.fail("6.7", "termination predicate holds on element %d before the last"
                 % hold.index(True))


def _result(node, rep, roles):
    res = node["result"]
    frm, mode = res["from"], node["rounding"]
    if frm == "len(elements)":
        if node["cardinality"] != "answer_bearing":
            rep.fail("3.6", "result.from len(elements) on an incidental node")
        got = float(len(node["elements"]))
    elif frm.startswith("last."):
        role = frm[5:]
        if role not in roles:
            rep.fail("3.6", "result.from names unknown role %r" % role)
            return
        got = do_round(node["elements"][-1][roles[role]], res["dp"], mode)
    else:
        rep.fail("3.6", "result.from %r is not one of the closed set of section 3.6"
                 % frm)
        return
    if abs(got - float(res["value"])) > 1e-9:
        rep.fail("6.8", "result.value %r != the result.from derivation %r"
                 % (res["value"], got))


# --------------------------------------------------------------------------

def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    lenient = "--lenient-roles" in argv
    check_filter = "--check-filter" in argv
    rows = json.load(open(args[0], encoding="utf-8"))
    npass = nfail = nunchecked = 0
    per_template = {}
    for r in rows:
        label = "%s seed %s" % (r.get("template_id"), r.get("seed"))
        rep = verify_node(r["trace_nodes"], label,
                          strict_roles=not lenient, check_filter=check_filter)
        d = per_template.setdefault(r.get("template_id"), [0, 0, set()])
        if rep.passed:
            npass += 1
            d[0] += 1
        else:
            nfail += 1
            d[1] += 1
            d[2].update(c for c, _ in rep.failures)
            if d[1] <= 2:
                print("FAIL %s" % label)
                for c, m in rep.failures[:4]:
                    print("      [%s] %s" % (c, m))
        if rep.unchecked:
            nunchecked += 1
    print()
    for t, (p, f, cl) in sorted(per_template.items()):
        print("%-40s pass %2d  fail %2d  %s"
              % (t, p, f, ("clauses " + ",".join(sorted(cl))) if cl else ""))
    print("TOTAL %d rows: %d pass, %d fail; %d rows carried unchecked carry_notes"
          % (len(rows), npass, nfail, nunchecked))
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
