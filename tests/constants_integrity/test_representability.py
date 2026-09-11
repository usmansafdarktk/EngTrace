"""C3.6 - D-025's representability assertion, a gate item on S2.

    python -m tests.constants_integrity.test_representability [--seeds N]
    python -m tests.constants_integrity.test_representability --selftest

D-025 (Phase 1 Reviewer A, F-5): a Phase 1 oracle re-solves an item from the
numbers its QUESTION prints, so a constant the template consumes at a precision
other than the one the question states makes gold and question disagree - and
two latent cases were hidden only by an accident of the tables. Its action:
"every constant consumed by a Phase 1 template is exactly representable at the
precision its question states it to".

The two cases are two different shapes, and the check is written for both:

  * STATED ROUNDED - the question itself prints the constant lossily. AISC SI Ix
    at 2 dp is printed `.1f` ("34.4 x 10^6 mm^4") while the template consumes
    3.441e-5 m^4.
  * CONSUMED ROUNDED - the question prints the constant exactly but the trace
    consumes a rounding of it. torsion.py prints "G = 77.15 GPa" and computes
    with `_as_printed(G * 1e9, ".2e")` = 7.72e10 Pa.

SCOPE. The Phase 1 templates are the ones T2 has an oracle for; every run prints
that set against instance_dump.TEMPLATES. A table is consumed by one when the
census's static consumer analysis says so.

METHOD, per numeric LEAF of such a table, per seed, per consumer. Nudge that one
leaf (census._nudge) and regenerate. The number tokens that DISAPPEAR from the
question and from the solution are the printed forms of the leaf and of whatever
was computed from it. Against the leaf's value v (exact decimal of repr(v)):

  exact token   - equals v times a power of ten (a unit prefix moves the point)
  rounded token - a STRICT rounding of v times a power of ten, to the token's own
                  last place, with at least MIN_SF significant figures (shorter
                  tokens coincide with rounded derived numbers far too easily)

  question has a rounded token and no exact one ...... STATED-ROUNDED     fail
  question has both .................................. STATED-TWO-WAYS    fail
  question exact, solution has a rounded token ....... CONSUMED-ROUNDED   fail
  question exact ..................................... STATED             pass
  question changes, no token represents v ............ DERIVED-ONLY       pass
      (v shapes printed numbers without being printed - a span chosen from a
      depth, or a rejection loop accepting a different row - so no oracle reads v)
  question unchanged, solution changed:
      an exact token of v is in the question ......... COARSE             pass
      (the nudge is below the display step of an exactly-stated value)
      otherwise ...................................... HIDDEN             fail
  the nudge makes generation raise ................... ERROR              fail

A leaf no seed uses is UNCOVERED: counted and printed, never passed as checked.

LIMITS, stated. A rounding shorter than MIN_SF significant figures is not
recognised, so a constant printed lossily at 1-2 s.f. reads DERIVED-ONLY. COARSE
and STATED accept any exact token among those that disappeared, so a derived
number that coincides exactly with v could mask a rounding elsewhere; STATED-TWO-
WAYS catches that only when the rounding is also in the question.
"""
from __future__ import annotations

import argparse
import ast
import copy
import importlib
import os
import sys
from collections import Counter
from decimal import (ROUND_HALF_DOWN, ROUND_HALF_EVEN, ROUND_HALF_UP, Decimal,
                     InvalidOperation, localcontext)

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, NUM_RE, Patch, _gen, _leaves, _nudge, _reload_all, census)

MIN_SF = 3
FAILING = ('STATED-ROUNDED', 'STATED-TWO-WAYS', 'CONSUMED-ROUNDED', 'HIDDEN', 'ERROR')


def _dec(x):
    if isinstance(x, str):
        return Decimal(x.replace(',', '').rstrip('.')).copy_abs()
    return Decimal(repr(x) if isinstance(x, float) else str(x)).copy_abs()


def exact(token, value):
    """The token is the value times an integer power of ten."""
    try:
        t, d = _dec(token), _dec(value)
    except InvalidOperation:
        return False
    if d == 0 or t == 0:
        return d == t
    with localcontext() as ctx:
        ctx.prec = 60
        r = t / d
        return r == Decimal(1).scaleb(r.adjusted())


def rounded(token, value, min_sf=MIN_SF):
    """The token is a strict rounding of the value times a power of ten, to the
    token's own last written place, carrying at least min_sf significant figures."""
    try:
        t, d = _dec(token), _dec(value)
    except InvalidOperation:
        return False
    if d == 0 or t == 0:
        return False
    digits = t.as_tuple().digits
    if len(digits) < min_sf:
        return False
    place = Decimal(1).scaleb(t.as_tuple().exponent)
    with localcontext() as ctx:
        ctx.prec = 60
        for k in (t.adjusted() - d.adjusted(), t.adjusted() - d.adjusted() - 1):
            scaled = d.scaleb(k)
            if scaled == t:
                return False
            for mode in (ROUND_HALF_UP, ROUND_HALF_EVEN, ROUND_HALF_DOWN):
                if scaled.quantize(place, rounding=mode) == t:
                    return True
    return False


def verdict(v, q0, s0, q1, s1):
    """The table above, for one leaf value and one seed's before/after texts."""
    if q1 is None:
        return 'ERROR', 'the nudge makes generation raise'
    if q0 == q1:
        if s0 == s1:
            return None, ''
        if any(exact(t, v) for t in NUM_RE.findall(q0)):
            return 'COARSE', ''
        return 'HIDDEN', 'the nudge changes the solution and never the question'
    gq = list((Counter(NUM_RE.findall(q0)) - Counter(NUM_RE.findall(q1))).elements())
    gs = list((Counter(NUM_RE.findall(s0)) - Counter(NUM_RE.findall(s1))).elements())
    eq = [t for t in gq if exact(t, v)]
    rq = [t for t in gq if rounded(t, v)]
    rs = [t for t in gs if rounded(t, v)]
    if rq and not eq:
        return 'STATED-ROUNDED', f'the question prints it as {rq[:3]}'
    if rq and eq:
        return 'STATED-TWO-WAYS', f'the question prints it as {eq[:2]} and as {rq[:2]}'
    if eq and rs:
        return 'CONSUMED-ROUNDED', f'the question states {eq[0]}, the solution consumes {rs[:3]}'
    if eq:
        return 'STATED', ''
    return 'DERIVED-ONLY', ''


def perturb_leaf(obj, target, path=()):
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        return _nudge(obj) if path == target else obj
    if isinstance(obj, dict):
        return {k: perturb_leaf(v, target, path + (('d', k),)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [perturb_leaf(v, target, path + (('l', i),)) for i, v in enumerate(obj)]
    if isinstance(obj, tuple):
        return tuple(perturb_leaf(v, target, path + (('t', i),)) for i, v in enumerate(obj))
    return obj


def leaf_name(path):
    return ''.join(f'[{k!r}]' for _kind, k in path)


def check_table(cmod, name, table, consumers, seeds):
    """Classify every (leaf, seed) for one table object - the real one or a plant -
    against its Phase 1 consumers {template_id: module}. Returns (results, uncovered)."""
    orig = getattr(cmod, name)
    modnames = sorted(set(consumers.values()))
    res = {tid: dict(counts=Counter(), failures=[], base_errors=0) for tid in consumers}
    seen = set()
    leaves = list(_leaves(table))
    try:
        with Patch(cmod, orig, table):
            mods, err = _reload_all(modnames)
            if err:
                raise RuntimeError(err)
            base = {tid: [_gen(getattr(mods[m], tid), s) for s in range(seeds)]
                    for tid, m in consumers.items()}
        for tid in consumers:
            res[tid]['base_errors'] = sum(q is None for q, _s in base[tid])
        for path, v in leaves:
            with Patch(cmod, orig, perturb_leaf(table, path)):
                mods, err = _reload_all(modnames)
                for tid, m in consumers.items():
                    r = res[tid]
                    if err:
                        r['counts']['ERROR'] += 1
                        r['failures'].append((None, path, v, 'ERROR', f'reload failed: {err}'))
                        continue
                    fn = getattr(mods[m], tid)
                    for s in range(seeds):
                        q0, s0 = base[tid][s]
                        if q0 is None:
                            continue
                        q1, s1 = _gen(fn, s)
                        cls, why = verdict(v, q0, s0, q1, s1)
                        if cls is None:
                            continue
                        seen.add(path)
                        r['counts'][cls] += 1
                        if cls in FAILING:
                            r['failures'].append((s, path, v, cls, why))
    finally:
        _reload_all(modnames)
    uncovered = [leaf_name(p) for p, _v in leaves if p not in seen]
    return res, uncovered


def phase1_consumers(branches=BRANCHES):
    """(oracles, [(branch, table, {template_id: module})]) for Phase 1 consumption."""
    from tests.template_integrity.checks.t2_roundtrip import available_oracles
    oracles = available_oracles()
    out = []
    for t in census(seeds=0, probe_on=False, branches=branches):
        cons = {tid: u['module'] for tid, u in t['consumers'].items()
                if tid in oracles and not any(str(x).startswith('copied literal')
                                              for x in u.get('via', []))}
        if cons:
            out.append((t['branch'], t['name'], cons))
    return oracles, out


def _dump_list():
    """instance_dump.TEMPLATES, read without importing: that module dumps at import."""
    path = os.path.join(REPO, 'tests', 'template_integrity', 'instance_dump.py')
    for node in ast.parse(open(path, encoding='utf-8').read()).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == 'TEMPLATES':
            return set(ast.literal_eval(node.value))
    return set()


def _print(branch, name, res, uncovered, limit=4):
    nfail = 0
    for tid, r in res.items():
        nfail += len(r['failures'])
        c = r['counts']
        print(f"  [{'ok' if not r['failures'] else 'FAIL'}] {branch[:5]}.{name:22s} {tid:40s} "
              + '  '.join(f'{k} {c[k]}' for k in ('STATED', 'DERIVED-ONLY', 'COARSE') if c[k])
              + (f"  | {'  '.join(f'{k} {c[k]}' for k in FAILING if c[k])}" if r['failures'] else '')
              + (f"  baseline errors {r['base_errors']}" if r['base_errors'] else ''))
        for s, path, v, cls, why in r['failures'][:limit]:
            print(f'        seed {s} {leaf_name(path)} = {v!r}: {cls}: {why}')
    if uncovered:
        print(f"        {len(uncovered)} leaf/leaves UNCOVERED at this seed count: "
              f"{', '.join(uncovered[:6])}{' ...' if len(uncovered) > 6 else ''}")
    return nfail


def run(seeds):
    oracles, pairs = phase1_consumers()
    listed = _dump_list()
    print(f'Phase 1 templates: {len(oracles)} with an oracle, {len(listed)} in '
          f'instance_dump.TEMPLATES; oracle-only {sorted(set(oracles) - listed)}, '
          f'dump-only {sorted(listed - set(oracles))}')
    failures = 0
    for branch, name, cons in pairs:
        cmod = importlib.import_module(f'data.templates.branches.{branch}.constants')
        res, uncovered = check_table(cmod, name, getattr(cmod, name), cons, seeds)
        failures += _print(branch, name, res, uncovered)
    print(f'{len(pairs)} tables consumed by a Phase 1 template, {seeds} seeds')
    print('all pass' if not failures else f'{failures} FAILURES')
    return 1 if failures else 0


def selftest(seeds=60):
    """D-025's own two demonstrations, applied in memory, through check_table.

    Written from D-025's record rather than from this detector, and of the two
    shapes it names: a shear-modulus family gaining a significant figure behind a
    question that prints it exactly and a trace that consumes it `.2e`
    (CONSUMED-ROUNDED), and AISC SI Ix gaining a decimal place behind a question
    that prints it `.1f` (STATED-ROUNDED). A plant counts only if the unplanted
    table has no failure and the plant has the failure class it was written for.
    """
    _oracles, pairs = phase1_consumers(branches=('mechanical_engineering', 'civil_engineering'))
    by = {(b, n): c for b, n, c in pairs}
    bad = []

    def plant_shear(tbl):
        return {k: round(v - 0.05, 2) for k, v in tbl.items()}

    def plant_ix(tbl):
        out = copy.deepcopy(tbl)
        for row in out.values():
            row['si']['Ix'] = round(row['si']['Ix'] + 0.01, 2)
        return out

    for label, branch, name, make, want in (
            ('SHEAR_MODULUS_VALUES one s.f. more', 'mechanical_engineering',
             'SHEAR_MODULUS_VALUES', plant_shear, 'CONSUMED-ROUNDED'),
            ('AISC SI Ix one d.p. more', 'civil_engineering', 'AISC_W_SHAPES', plant_ix,
             'STATED-ROUNDED')):
        cmod = importlib.import_module(f'data.templates.branches.{branch}.constants')
        cons = by.get((branch, name))
        if not cons:
            bad.append(f'{label}: no Phase 1 consumer found - nothing to plant against')
            continue
        clean, _u = check_table(cmod, name, getattr(cmod, name), cons, seeds)
        planted, _u = check_table(cmod, name, make(getattr(cmod, name)), cons, seeds)
        n_clean = sum(len(r['failures']) for r in clean.values())
        hits = [f for r in planted.values() for f in r['failures'] if f[3] == want]
        ok = n_clean == 0 and bool(hits)
        print(f"  [{'ok' if ok else 'FAIL'}] {label}: clean {n_clean} failure(s), planted "
              f"{len(hits)} {want}"
              + (f' - e.g. seed {hits[0][0]} {leaf_name(hits[0][1])} = {hits[0][2]!r}: {hits[0][4][:70]}'
                 if hits else ''))
        if not ok:
            bad.append(label)
    cases = (('exact', '7.72e+10', 77.2, True), ('exact', '7.72e+10', 77.15, False),
             ('exact', '29,000', 29000, True), ('exact', '0.3', 0.30, True),
             ('exact', '0', 0, True), ('exact', '1.0', 0, False),
             ('rounded', '7.72e+10', 77.15, True), ('rounded', '7.72e+10', 77.2, False),
             ('rounded', '34.4', 34.41, True), ('rounded', '10.0', 9.996, True),
             ('rounded', '1.0', 0.9982, False), ('rounded', '0.022', 0.0216, False))
    wrong = [c for c in cases if (exact if c[0] == 'exact' else rounded)(c[1], c[2]) != c[3]]
    print(f"  [{'ok' if not wrong else 'FAIL'}] exact()/rounded(): {len(cases)} decimal cases"
          + (f' - wrong: {wrong}' if wrong else ''))
    bad += [f'{c[0]}({c[1]!r}, {c[2]!r}) is not {c[3]}' for c in wrong]
    for b in bad:
        print('  - ' + b)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=100)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    sys.exit(selftest() if a.selftest else run(a.seeds))
