"""C1.4 / C1.5 - every numeric table declares what it is and what units it is in.

    python -m tests.constants_integrity.test_table_metadata
    python -m tests.constants_integrity.test_table_metadata --selftest

Spec C1.4 asked for "a machine-readable unit declaration per table" and C1
skipped it without anyone noticing, because nothing tested it. Reviewer G's
words: *"an unenforced deliverable is a suggestion."* This is the enforcement.

  M1  every numeric table (census P-TABLE) declares `# @kind:` from census.KINDS
  M2  every numeric table declares `# @units:`; each declared key names at least
      one field of the live table; every field is covered by a declaration; and
      every unit parses under the grammar below
  M3  every numeric table declares `# @domain:` - `key=value unit`,
      `key=lo..hi unit`, `key=word`, or `none (reason)` - OR is listed in
      domain_worklist.txt, the C3 worklist. The list RATCHETS: a listed table
      that now declares a domain fails until it is removed, and a listed name
      that is not a numeric table fails.

The fields are written in the table's contiguous comment header (see
docs/references/README.md). A FIELD is a census.fields() path - every leaf
reached by the same path, row keys wildcarded. A declared key maps onto fields
as: `E_GPa` -> `['E_GPa']`, `us.Ix` -> `['us']['Ix']`, `[1]` -> `[1]`, matched
after the row wildcard and as a prefix, so `e` covers both `['e']` and a
range-valued `['e'][0]`.
"""
from __future__ import annotations

import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, KINDS, fields, header_fields, static_tables)

WORKLIST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'domain_worklist.txt')

# --------------------------------------------------------------------------
# The unit grammar. ONE grammar for constants and - Phase 6, D6.11 - answers.
#
#   unit    := '1' | 'per-row(key)' | expr
#   expr    := term (('*' | '/') term)*
#   term    := factor ('^' exp)?
#   factor  := SYMBOL | SCALE | '1' | '(' expr ')'
#   exp     := INT | '-' INT | '(' INT '/' INT ')' | 'n'
#   SCALE   := 1e<k>, only as the first factor of a product
# --------------------------------------------------------------------------

SYMBOLS = frozenset("""
    m cm mm um micron km in ft g kg lb lbf mol s min h yr month week shift sample
    K degC Pa kPa MPa GPa bar atm psi ksi N kN J kJ W V A F Hz MHz rad deg ohm
    USD angstrom mL L %
""".split())

_TOK = re.compile(r'\s*(1e[-+]?\d+|\d+|[A-Za-z%]+|[*/^()\-])')


def _tokens(s):
    pos, out = 0, []
    while pos < len(s):
        m = _TOK.match(s, pos)
        if not m or m.end() == pos:
            raise ValueError(f'cannot tokenise {s[pos:]!r}')
        out.append(m.group(1))
        pos = m.end()
    return out


def parse_unit(s):
    """Raise ValueError unless `s` is a well-formed unit."""
    s = s.strip()
    if s == '1' or s == 'per-row(key)':
        return True
    toks = _tokens(s)
    i = 0

    def peek():
        return toks[i] if i < len(toks) else None

    def take(expected=None):
        nonlocal i
        t = peek()
        if t is None or (expected is not None and t != expected):
            raise ValueError(f'expected {expected!r} at {toks[i:]!r} in {s!r}')
        i += 1
        return t

    def factor(first):
        t = peek()
        if t == '(':
            take('(')
            expr()
            take(')')
        elif t is not None and t.startswith('1e'):
            if not first:
                raise ValueError(f'scale factor {t} not leading in {s!r}')
            take()
            if peek() != '*':
                raise ValueError(f'scale factor {t} must multiply a unit in {s!r}')
        elif t == '1':
            take()
        elif t in SYMBOLS:
            take()
        else:
            raise ValueError(f'unknown unit symbol {t!r} in {s!r}')

    def exponent():
        t = peek()
        if t == '(':
            take('(')
            if not take().isdigit():
                raise ValueError(f'bad fractional exponent in {s!r}')
            take('/')
            if not take().isdigit():
                raise ValueError(f'bad fractional exponent in {s!r}')
            take(')')
        elif t == '-':
            take('-')
            if not take().isdigit():
                raise ValueError(f'bad negative exponent in {s!r}')
        elif t == 'n' or (t is not None and t.isdigit()):
            take()
        else:
            raise ValueError(f'bad exponent {t!r} in {s!r}')

    def term(first):
        factor(first)
        if peek() == '^':
            take('^')
            exponent()

    def expr():
        term(True)
        while peek() in ('*', '/'):
            take()
            term(False)

    expr()
    if i != len(toks):
        raise ValueError(f'trailing {toks[i:]!r} in {s!r}')
    return True


def _split_top(s):
    parts, depth, cur = [], 0, ''
    for ch in s:
        depth += ch == '('
        depth -= ch == ')'
        if ch == ',' and depth == 0:
            parts.append(cur.strip())
            cur = ''
        else:
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def _keypath(key):
    if key.startswith('['):
        return key
    return ''.join(f"[{k!r}]" for k in key.split('.'))


def _row_relative(field, obj):
    if isinstance(obj, (dict, list)) and field.startswith('[*]'):
        return field[3:]
    return field


def check_units(name, decl, obj):
    """M2 for one table. Returns a list of failure strings."""
    out = []
    flds = fields(obj)
    if not flds:
        return [f'{name}: no numeric field in the live table']
    parts = _split_top(decl)
    keyed = [p for p in parts if '=' in p.split('(')[0]]
    if keyed and len(keyed) != len(parts):
        return [f'{name}: @units mixes a table-wide unit with keyed units: {decl!r}']
    if not keyed:
        if len(parts) != 1:
            return [f'{name}: @units has {len(parts)} table-wide units: {decl!r}']
        try:
            parse_unit(parts[0])
        except ValueError as exc:
            out.append(f'{name}: {exc}')
        if parts[0] == 'per-row(key)':
            out.extend(_per_row(name, obj))
        return out
    covered = set()
    for p in keyed:
        key, unit = (x.strip() for x in p.split('=', 1))
        kp = _keypath(key)
        hit = {f for f in flds if _row_relative(f, obj) == kp
               or _row_relative(f, obj).startswith(kp + '[')}
        if not hit:
            out.append(f'{name}: @units key {key!r} names no field '
                       f'(fields: {sorted(flds)})')
        covered |= hit
        try:
            parse_unit(unit)
        except ValueError as exc:
            out.append(f'{name}: {key}: {exc}')
        if unit == 'per-row(key)':
            out.extend(_per_row(name, obj))
    for f in sorted(set(flds) - covered):
        out.append(f'{name}: field {f} has no declared unit')
    return out


def _per_row(name, obj):
    out = []
    for key in (obj if isinstance(obj, dict) else []):
        m = re.search(r'\(([^()]+)\)\s*$', str(key))
        if not m:
            out.append(f'{name}: per-row(key) but row {key!r} names no unit in parentheses')
            continue
        try:
            parse_unit(m.group(1))
        except ValueError as exc:
            out.append(f'{name}: row {key!r}: {exc}')
    return out


_DOMAIN_PAIR = re.compile(
    r'^(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<val>-?\d+(?:\.\d+)?(?:e-?\d+)?'
    r'(?:\.\.-?\d+(?:\.\d+)?(?:e-?\d+)?)?|[A-Za-z][\w-]*)(?:\s+(?P<unit>.+))?$')


def check_domain(name, decl):
    decl = decl.strip()
    if re.match(r'^none \(.+\)$', decl):
        return []
    out = []
    for p in _split_top(decl):
        m = _DOMAIN_PAIR.match(p)
        if not m:
            out.append(f'{name}: @domain entry {p!r} is not key=value [unit]')
            continue
        if m.group('unit'):
            try:
                parse_unit(m.group('unit'))
            except ValueError as exc:
                out.append(f'{name}: @domain {m.group("key")}: {exc}')
    return out


def read_worklist(path=WORKLIST):
    if not os.path.exists(path):
        return set()
    return {ln.strip() for ln in open(path, encoding='utf-8')
            if ln.strip() and not ln.lstrip().startswith('#')}


def check_branch(branch, src, worklist):
    """M1-M3 for one constants.py source. Returns (checks, failures)."""
    ns = {}
    exec(compile(src, f'{branch}/constants.py', 'exec'), ns)     # noqa: S102
    failures, checks = [], 0
    names = set()
    for t in static_tables(src):
        name = t['name']
        names.add(name)
        label = f'{branch}.{name}'
        d = header_fields(t['header'])
        checks += 3
        kind = d.get('kind', '').split(' ')[0]
        if not kind:
            failures.append(f'M1 {label}: no @kind')
        elif kind not in KINDS:
            failures.append(f'M1 {label}: @kind {kind!r} is not one of {KINDS}')
        if 'units' not in d:
            failures.append(f'M2 {label}: no @units')
        else:
            failures.extend(f'M2 {branch}.{f}' for f in check_units(name, d['units'], ns[name]))
        if 'domain' in d:
            failures.extend(f'M3 {branch}.{f}' for f in check_domain(name, d['domain']))
            if label in worklist:
                failures.append(f'M3 {label}: declares @domain but is still on the '
                                f'worklist - remove it (ratchet)')
        elif label not in worklist:
            failures.append(f'M3 {label}: no @domain and not on the C3 worklist')
    return checks, failures, {f'{branch}.{n}' for n in names}


def run():
    worklist = read_worklist()
    failures, checks, seen = [], 0, set()
    for branch in BRANCHES:
        src = open(os.path.join(BRANCHES_DIR, branch, 'constants.py'),
                   encoding='utf-8').read().replace('\r\n', '\n')
        c, f, names = check_branch(branch, src, worklist)
        checks, seen = checks + c, seen | names
        failures += f
    for stale in sorted(worklist - seen):
        failures.append(f'M3 worklist names {stale}, which is not a numeric table')
    checks += 1
    for f in failures:
        print('  - ' + f)
    print(f'{checks} checks over {len(seen)} tables; '
          f'{len(worklist & seen)} tables on the C3 @domain worklist')
    print('all pass' if not failures else f'{len(failures)} FAILURES')
    return 1 if failures else 0


# --------------------------------------------------------------------------
# Self-test. Each plant is written from M1-M3's definitions, not from the
# code above, and each property gets two plants of different surface form.
# --------------------------------------------------------------------------

_GOOD = '''
# @kind: property
# @units: E_GPa=GPa, nu=1
# @domain: T=293.15 K, form=wrought
MAT = {"Steel": {"E_GPa": 200.0, "nu": 0.3}}

# @kind: range
# @units: Hz
# @domain: none (a sampling window)
FREQ = (50, 2000)

# @kind: property
# @units: target=per-row(key), frac=1
SPC = {"diameter (mm)": {"target": (10, 80), "frac": (0.001, 0.01)}}
'''

_PLANTS = {
    # M1 - absent, and present-but-not-a-kind
    'M1 no @kind': ('# @kind: property\n# @units: E_GPa', '# @units: E_GPa'),
    'M1 unknown @kind': ('# @kind: range\n', '# @kind: sampling-window\n'),
    # M2 - absent, a key naming no field, an uncovered field, a unit that does
    #      not parse, and a per-row unit whose row names none
    'M2 no @units': ('# @units: Hz\n', ''),
    'M2 key names no field': ('E_GPa=GPa, nu=1', 'E_Gpa=GPa, nu=1'),
    'M2 field without unit': ('E_GPa=GPa, nu=1', 'E_GPa=GPa'),
    'M2 unit does not parse': ('# @units: Hz', '# @units: kg/m3'),
    'M2 unknown symbol': ('E_GPa=GPa', 'E_GPa=furlong'),
    'M2 per-row key names no unit': ('"diameter (mm)"', '"diameter"'),
    # M3 - absent and not listed; malformed
    'M3 no @domain, not on worklist': ('# @domain: none (a sampling window)\n', ''),
    'M3 malformed @domain': ('T=293.15 K, form=wrought', 'at room temperature'),
}


def selftest():
    bad = []
    _c, clean, _ = check_branch('plant', _GOOD, {'plant.SPC'})
    if clean:
        # A plant is judged by the failure it ADDS. Over a failing fixture every
        # plant of the same code "passes" whether or not it was detected - which
        # is exactly how this self-test's first version passed all three M3
        # plants while its own fixture was being rejected for an uppercase key.
        print('  the clean fixture itself fails - no plant can be judged:')
        for x in clean:
            print('    ' + x)
        return 1
    for label, (old, new) in _PLANTS.items():
        assert _GOOD.count(old) == 1, f'plant {label!r} anchor not unique'
        _c, f, _ = check_branch('plant', _GOOD.replace(old, new), {'plant.SPC'})
        code = label.split()[0]
        if not any(x.startswith(code) and x not in clean for x in f):
            bad.append(f'{label}: NOT detected (got {f})')
    # the ratchet, both directions
    _c, f, _ = check_branch('plant', _GOOD, {'plant.SPC', 'plant.MAT'})
    if not any('still on the worklist' in x for x in f):
        bad.append('M3 ratchet: a table with @domain left on the worklist was not caught')
    names = check_branch('plant', _GOOD, set())[2]
    if 'plant.NOT_A_TABLE' in names:
        bad.append('check_branch reported a name that is not a table')
    for label in list(_PLANTS) + ['M3 ratchet']:
        print(f'  [{"FAIL" if any(b.startswith(label) for b in bad) else "ok"}] {label}')
    for b in bad:
        print('  - ' + b)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
