"""Run §3.8's audit against the nodes actually shipping (D3.3 v1.5).

    python -m tests.trace_schema.audit_3_8

§3.8 says: every value a check consumes is either declared once on the node and
checked against something, or recomputed from something that is. Its three audit
questions have caught three shipped defects — `budget_total` restated per element
(1.2, 1.3), evaluation frames unchecked (1.0–1.3), and `frame_relations` with no
coverage requirement (1.4). Every one was found *after* shipping, twice by a
reviewer.

This runs the audit mechanically instead, so it is cheap enough to run before.
It is not a conformance verifier — the reviewers' verifiers do that — it answers
only §3.8's questions:

  Q0  can every value be resolved through a role map, or only by literal name?
  Q1  is any value a check reads restated per element rather than declared once?
  Q2  is any DERIVED value declared once but recomputed by nothing?
  Q3  does a recomputation cover every value it should, or only those present?
"""
from __future__ import annotations

import json
import sys

# Values the spec classifies as problem data: declaring once is the whole
# requirement, because nothing on the node binds the node to its question
# (§3.8, D-039). Anything else a check reads must be recomputed.
PROBLEM_DATA = {
    'constants', 'universe', 'precedences', 'item_measures', 'budget',
    'tolerance', 'max_elements', 'symbol_precision', 'rounding',
    'roles', 'trial_roles', 'candidate_roles', 'evaluation_roles',
    'element_symbols', 'trial_symbols', 'eligible_symbols',
    'evaluation_symbols', 'preamble_symbols',
}


def _r(node: dict, kind: str, role: str) -> str:
    """Resolve a type-level role to this node's own symbol name (§3.3).

    Reaching for a literal symbol name is non-conforming under §8B.11, and the
    first version of this file did exactly that for `decision` -- it raised
    KeyError on a renamed node, including a renamed node that was corrupt, so it
    produced no finding for precisely the class of node the role machinery
    exists to serve (Reviewer D2, round 6). It failed safe rather than green,
    which is why it was not blocking, but silently inapplicable is not applicable.
    """
    return node[kind][role]


def audit(node: dict) -> list[str]:
    out: list[str] = []
    t = node['node_type']
    term = node['termination']

    for kind in ('roles',) + (('trial_roles', 'candidate_roles') if t == 'decision'
                              else ('evaluation_roles',)):
        if kind not in node:
            out.append(f'Q0 FAIL: no {kind}; values can only be resolved by '
                       f'literal symbol name, which is non-conforming (§8B.11)')
    if out:
        return out

    # Q1 - restated per element rather than declared once
    if t == 'decision':
        budget_total = _r(node, 'roles', 'budget_total')
        item = _r(node, 'candidate_roles', 'item')
        measure = _r(node, 'candidate_roles', 'measure')
        candidates = _r(node, 'trial_roles', 'candidates')
        trials = _r(node, 'roles', 'trials')

        budgets = {e[budget_total] for e in node['elements']}
        if 'budget' not in term:
            out.append('Q1 FAIL: no termination.budget; the budget is restated '
                       'per element and checked against nothing')
        elif budgets != {term['budget']}:
            out.append(f'Q1 FAIL: element budgets {budgets} disagree with '
                       f'termination.budget {term["budget"]}')

        measures: dict = {}
        for e in node['elements']:
            for tr in e[trials]:
                for c in tr[candidates]:
                    measures.setdefault(c[item], set()).add(c[measure])
        multi = {k: v for k, v in measures.items() if len(v) > 1}
        if multi:
            out.append(f'Q1 FAIL: measures restated inconsistently: {multi}')
        if 'item_measures' not in term:
            out.append('Q1 FAIL: no termination.item_measures')
        else:
            off = {k: v for k, v in measures.items()
                   if term['item_measures'].get(k) not in v}
            if off:
                out.append(f'Q1 FAIL: candidate measures disagree with '
                           f'item_measures: {sorted(off)}')

    # Q2 - derived values that nothing recomputes
    if t == 'iteration':
        if 'update_relation' not in node:
            out.append('Q2 FAIL: the update is derived and undeclared')
        if 'frame_relations' not in node:
            out.append('Q2 FAIL: frame symbols are derived and unchecked')
        if 'preamble' in node and 'preamble_binding' not in node:
            out.append('Q2 FAIL: a preamble is present but unbound to elements[0]')
    else:
        if 'filter_relation' not in node.get('selection', {}):
            out.append('Q2 FAIL: admissibility is derived and only asserted')
    if 'from' not in node.get('result', {}):
        out.append('Q2 FAIL: result.value is derived with no stated derivation')

    # Q3 - coverage, not merely presence
    if t == 'iteration' and 'frame_relations' in node:
        iterate = _r(node, 'evaluation_roles', 'iterate')
        need = {s for s in node['evaluation_symbols'] if s != iterate}
        have = [k for k, _ in node['frame_relations']]
        missing = need - set(have)
        dupes = {k for k in have if have.count(k) > 1}
        extra = set(have) - need
        if missing:
            out.append(f'Q3 FAIL: frame symbols with no relation: {sorted(missing)} '
                       f'-- they are exempt from §4.3.9, which is variant 5')
        if dupes:
            out.append(f'Q3 FAIL: symbols with more than one relation: {sorted(dupes)}')
        if extra:
            out.append(f'Q3 FAIL: relations for non-frame symbols: {sorted(extra)}')
        for r in ('iterate', 'residual'):
            if _r(node, 'evaluation_roles', r) not in node['evaluation_symbols']:
                out.append(f'Q3 FAIL: evaluation role {r!r} names a symbol absent '
                           f'from evaluation_symbols, so coverage cannot demand it')
    elif t == 'decision':
        universe = set(term['universe'])
        if set(term.get('item_measures', {})) != universe:
            out.append('Q3 FAIL: item_measures does not cover the universe exactly')
        if set(term.get('precedences', {})) != universe:
            out.append('Q3 FAIL: precedences does not cover the universe exactly')
    return out


def main(path: str) -> int:
    rows = json.load(open(path, encoding='utf-8'))
    failures = 0
    seen = {}
    for r in rows:
        problems = audit(r['trace_nodes'])
        if problems:
            failures += 1
            key = r['template_id']
            if key not in seen:
                seen[key] = problems
                print(f'{key} seed {r["seed"]}:')
                for p in problems:
                    print(f'   {p}')
    print(f'\n3.8 audit over {len(rows)} nodes: '
          f'{len(rows) - failures} clean, {failures} with findings')
    return 1 if failures else 0


if __name__ == '__main__':
    p = (sys.argv[1] if len(sys.argv) > 1
         else 'docs/re-implementation-sep/phase3_conformance/traces.json')
    raise SystemExit(main(p))
