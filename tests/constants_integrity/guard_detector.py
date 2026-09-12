"""C3.10 - a STATIC guard detector: is every read of a table, in a consumer, a guard?

    python -m tests.constants_integrity.guard_detector
    python -m tests.constants_integrity.guard_detector --selftest

The census's perturbation probe cannot tell a guard from a crash (C1 Reviewer G, G-2):
a nudged table that makes a template raise reads as ERROR, and a table read only in a
screen reads as NO-EFFECT or ERROR. C1 closed that gap with human declarations,
`# @given: guard (evidence)`. Spec C3.10 asks for a static detector to replace those
declarations where it can.

A READ is every load - inside the template function and every module-level function
it calls (census._reach) - of the table, a stand-in alias of it (census module and
constants aliases, through census.template_reaches), or a LOCAL name whose value
derives from the table (census._Flow). Each read is classified by where it sits:

  GUARD     anywhere inside an `assert` statement - its test or its message, since a
            failure message is not an item; inside the test of an `if`, a `while` or
            a conditional expression; inside a comprehension's `if` clause
  TRANSFER  the whole value of an assignment whose every target is itself a derived
            local (`lo, hi = TABLE`): neither - the targets' own reads decide
  VALUE     anything else: an expression that computes, prints or draws, and the
            iterable of a `for` loop (`for m in range(T[0], T[1] + 1)` generates values)

A (table, consumer) pair is STATIC-GUARD when it has at least one GUARD read and no
VALUE read; VALUE when any read is VALUE; UNUSED when every read is TRANSFER.

EVIDENCE, NOT PROOF. Derived locals come from the census flow walk's `ever` map,
which over-approximates (a name ever bound from the table counts at every read), so
within what the walk can SEE the detector errs toward VALUE. What it cannot see can
still hide a VALUE read and produce a false STATIC-GUARD. The first run found one:
attributes_control_charts.py builds `_T27_PLANS = _t27_admissible_plans()` at module
level from SPC_NUM_SUBGROUPS and P_CHART_SUBGROUP_N, and the p-chart template draws
`random.choice(_T27_PLANS)`. The census walk took a module-level name's tables from
the names its value loads, `{_t27_admissible_plans}`, never from what that function
reads, so the draw was invisible and both tables read as guard-only. (This file's
first docstring said the detector "cannot call a value GUARD". It could.) Access the
walk still cannot follow - getattr, exec, a function stored in a container - would
hide a read the same way. It proves a syntactic property only; whether a guard's
outcome can still move an item (a screen that redraws changes the draw sequence) is
the probe's question, not this one.

WHAT IT REPORTS. Each `@given: guard` declaration is VERIFIED when the table has
census consumers and every one is STATIC-GUARD or UNUSED, and is otherwise left to the
committed evidence it cites. Every other (table, consumer) pair that is STATIC-GUARD
is listed too - for a table the probe could not settle, a candidate for a declaration
the detector already backs. Gates nothing.
"""
from __future__ import annotations

import ast
import os
import shutil
import sys
import tempfile
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, _Flow, _free_loads, constants_aliases, header_fields,
    numeric_tables, template_reaches)

GUARD, TRANSFER, VALUE = 'GUARD', 'TRANSFER', 'VALUE'


def _parents(fn):
    par = {}
    for node in ast.walk(fn):
        for child in ast.iter_child_nodes(node):
            par[child] = node
    return par


def _all_derived(target, derived):
    if isinstance(target, ast.Name):
        return target.id in derived
    if isinstance(target, (ast.Tuple, ast.List)):
        return all(_all_derived(e, derived) for e in target.elts)
    if isinstance(target, ast.Starred):
        return _all_derived(target.value, derived)
    return False


def context(name_node, par, derived):
    """GUARD, TRANSFER or VALUE for one load, from its ancestors up to its statement."""
    child, node = name_node, par.get(name_node)
    while node is not None:
        if isinstance(node, ast.Assert):
            return GUARD
        if isinstance(node, (ast.If, ast.While, ast.IfExp)) and child is node.test:
            return GUARD
        if isinstance(node, ast.comprehension) and child in node.ifs:
            return GUARD
        if isinstance(node, (ast.For, ast.AsyncFor)) and child is node.iter:
            return VALUE
        if isinstance(node, ast.Assign) and child is node.value:
            return TRANSFER if all(_all_derived(t, derived) for t in node.targets) else VALUE
        if isinstance(node, ast.AnnAssign) and child is node.value:
            return TRANSFER if _all_derived(node.target, derived) else VALUE
        if isinstance(node, ast.stmt):
            return VALUE
        child, node = node, par.get(node)
    return VALUE


def analyse(branch_dir, module_root, tables, calias):
    """{(table, template): {verdict, module, reads: [(line, name, context)]}}."""
    out = {}
    for modname, name, reach, sources, here, calls in template_reaches(branch_dir, module_root, tables, calias):
        free = set().union(*(_free_loads(r) for r in reach))
        flow = _Flow(sources, calls)
        for r in reach:
            flow.walk(r.body, {})
        parents = [(r, _parents(r)) for r in reach]
        for t in tables:
            stand_ins = {t} | {a for a, s in here.items() if t in s}
            if not stand_ins & free:
                continue
            derived = {n for n, ts in flow.ever.items() if t in ts} - stand_ins
            watch = stand_ins | derived
            reads = []
            for r, par in parents:
                for node in ast.walk(r):
                    if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load) and node.id in watch:
                        reads.append((node.lineno, node.id, context(node, par, derived)))
            kinds = {c for _l, _n, c in reads}
            verdict = VALUE if VALUE in kinds else ('STATIC-GUARD' if GUARD in kinds else 'UNUSED')
            out[(t, name)] = dict(verdict=verdict, module=modname, reads=sorted(reads))
    return out


def run():
    verified = unverified = 0
    for branch in BRANCHES:
        src = open(os.path.join(BRANCHES_DIR, branch, 'constants.py'), encoding='utf-8').read()
        tables = numeric_tables(src)
        names = sorted(t['name'] for t in tables)
        res = analyse(os.path.join(BRANCHES_DIR, branch), REPO, names, constants_aliases(src, names))
        by_table = defaultdict(dict)
        for (t, tid), r in res.items():
            by_table[t][tid] = r
        for t in tables:
            given = header_fields(t['header']).get('given', '')
            cons = by_table.get(t['name'], {})
            if given.split(' ')[0] == 'guard':
                values = {tid for tid, r in cons.items() if r['verdict'] == VALUE}
                ok = bool(cons) and not values
                verified += ok
                unverified += not ok
                print(f"  @given: guard  {branch[:5]}.{t['name']}: "
                      + ('VERIFIED statically' if ok else
                         f'NOT verified (VALUE in {sorted(values) or "no consumer"}) - keeps its cited evidence'))
        for t, cons in sorted(by_table.items()):
            guards = sorted(tid for tid, r in cons.items() if r['verdict'] == 'STATIC-GUARD')
            if guards:
                others = sorted(f"{tid} {r['verdict']}" for tid, r in cons.items() if r['verdict'] != 'STATIC-GUARD')
                print(f'  {branch[:5]}.{t}: STATIC-GUARD in {guards}' + (f'; also {others}' if others else ''))
    print(f'@given: guard declarations: {verified} verified statically, {unverified} left to their evidence')
    return 0


_SELFTEST_CONSTANTS = '''
FLOOR = 1.5
BAND = (0.01, 1.0)
LIMIT = 3000.0
WIN = (20, 30)
SCALE = 2.0
MIXED = 5.0
LAUNDERED = (0.2, 0.8)
PLANWIN = (20, 30)
'''

_SELFTEST_TEMPLATES = '''
import random
from pkg.constants import FLOOR, BAND, LIMIT, WIN, SCALE, MIXED, LAUNDERED, PLANWIN


def template_direct_if():
    while True:
        s = random.uniform(0, 3)
        if s < FLOOR:
            continue
        return f"s = {s}", "a"


def template_unpacked_assert():
    k = random.uniform(0.02, 0.09)
    lo, hi = BAND
    assert lo <= k <= hi, f"k outside {lo}-{hi}"
    return f"k = {k}", "a"


def template_assert_message():
    t = random.uniform(500, 2500)
    assert t <= LIMIT, f"{LIMIT} K exceeded"
    return f"T = {t}", "a"


def template_range_loop():
    plans = []
    for m in range(WIN[0], WIN[1] + 1):
        plans.append(m)
    return f"m = {random.choice(plans)}", "a"


def template_value():
    x = random.uniform(1, 2)
    return f"x = {x}", f"{x * SCALE}"


def template_mixed():
    x = random.uniform(1, 9)
    if x > MIXED:
        x = x - 1
    return f"x = {x}, cap {MIXED}", "a"


def template_unpacked_then_printed():
    p = random.uniform(0.3, 0.7)
    lo, hi = LAUNDERED
    assert lo <= p <= hi
    return f"p = {p}, at most {hi}", "a"


def _admissible():
    out = []
    for m in range(PLANWIN[0], PLANWIN[1] + 1):
        out.append(m)
    return out


_PLANS = _admissible()


def template_module_plans():
    m = random.choice(_PLANS)
    m_lo, m_hi = PLANWIN
    assert m_lo <= m <= m_hi
    return f"m = {m}", "a"
'''


def selftest():
    """Plants written from the class definition: three guard shapes of different
    surface form (a direct `if` screen, an unpacked alias read by an `assert`, an
    `assert` whose message prints the table) and four value shapes (arithmetic, a
    range loop that generates values, a guard plus a print of the same table, and an
    unpacked alias that one read guards and another prints - TRANSFER must not
    launder the print)."""
    tmp = tempfile.mkdtemp(prefix='guard_selftest_')
    try:
        pkg = os.path.join(tmp, 'pkg')
        os.makedirs(pkg)
        open(os.path.join(pkg, '__init__.py'), 'w').close()
        open(os.path.join(pkg, 'constants.py'), 'w', encoding='utf-8').write(_SELFTEST_CONSTANTS)
        open(os.path.join(pkg, 'tmpl.py'), 'w', encoding='utf-8').write(_SELFTEST_TEMPLATES)
        names = sorted(t['name'] for t in numeric_tables(_SELFTEST_CONSTANTS))
        res = analyse(pkg, tmp, names, constants_aliases(_SELFTEST_CONSTANTS, names))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    want = {
        ('FLOOR', 'template_direct_if'): 'STATIC-GUARD',
        ('BAND', 'template_unpacked_assert'): 'STATIC-GUARD',
        ('LIMIT', 'template_assert_message'): 'STATIC-GUARD',
        ('WIN', 'template_range_loop'): VALUE,
        ('SCALE', 'template_value'): VALUE,
        ('MIXED', 'template_mixed'): VALUE,
        ('LAUNDERED', 'template_unpacked_then_printed'): VALUE,
        # the _T27_PLANS shape: an assert on the table AND a draw from a module-level
        # list built by calling a function that reads it - VALUE, not STATIC-GUARD
        ('PLANWIN', 'template_module_plans'): VALUE,
    }
    bad = []
    for key, verdict in want.items():
        got = res.get(key, {}).get('verdict')
        ok = got == verdict
        print(f"  [{'ok' if ok else 'FAIL'}] {key[0]} in {key[1]}: {got} (planted {verdict})")
        if not ok:
            bad.append(key)
    extra = sorted(k for k in res if k not in want)
    if extra:
        print(f'  [FAIL] pairs the fixture does not plant: {extra}')
        bad.append('extra')
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
