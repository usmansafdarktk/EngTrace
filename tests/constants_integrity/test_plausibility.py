"""C3.2 - plausibility suites: cross-property consistency, checked against artefacts.

    python -m tests.constants_integrity.test_plausibility
    python -m tests.constants_integrity.test_plausibility --selftest

C2's plausibility suite was green while never looking where the data was used (D-032);
this one asks a different question that suite could not: **do a table's own columns agree
with each other, and is every number a tag invents printed somewhere on disk?** A table
can resolve leaf by leaf and still be internally inconsistent - `MATERIAL_PROPERTIES`
carries an SI and a US column for the same modulus, and nothing until now compared them.

  P1  MATERIAL_PROPERTIES' E_GPa and E_ksi are companion values of ONE modulus. They
      must agree through the SP 811 ksi->MPa factor, to within the rounding slack their
      OWN precision implies: half a unit in the last significant digit of each. That
      bound is derived per row, not chosen - a threshold picked to make rows pass is
      not an oracle.
  P2  every `scale=` in every citation tag is either a pure power of ten (a decimal
      re-scaling, which needs no source) or a factor an artefact PRINTS. SP 811 is
      searched numerically, so the check does not depend on a page number being right.
      A scale nobody prints is a number invented inside a tag.
  P3  SPECIFIC_GRAVITY_RANGES lies inside the only bound Das states for Gs - "most of
      the values fall within a range of 2.6 to 2.9" (PGE §2.6). LOCAL-ONLY: reported
      UNRESOLVABLE-FROM-CLONE, never failed, when the book is absent (C1.2).

The flame-temperature NOTE that C3.2 also owns lives in its own suite
(test_chemical_thermochemistry.py F-2), where the template it measures is.
"""
from __future__ import annotations

import json
import math
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.test_citations_resolve import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, REFS, extract_tags, parse_citation)

SP811 = os.path.join(REFS, 'nist_sp811', 'nistspecialpublication811e2008.pdf')
KSI_TO_MPA = 6.894757          # SP 811: kip per square inch -> pascal, 6.894 757 E+06
DAS_GS_LO, DAS_GS_HI = 2.6, 2.9


# ---------------------------------------------------------------- helpers
def slack(x):
    """Half a unit in the last significant digit of `x` as written.

    29000 is written to 2 s.f. (the trailing zeros are placeholders), so its slack is
    500; 68.9 is written to 3 s.f. and its slack is 0.05. repr() recovers the digits the
    author typed, which is what the rounding slack is a property of.
    """
    s = repr(float(x))
    if 'e' in s or 'E' in s:
        mant, _, exp = s.partition('e')
        digits = len(mant.replace('-', '').replace('.', '').rstrip('0')) or 1
        return 0.5 * 10 ** (int(exp) - digits + 1)
    if '.' in s:
        frac = len(s.split('.')[1].rstrip('0'))
        if frac:
            return 0.5 * 10 ** -frac
        s = s.split('.')[0]
    s = s.rstrip('0')
    trailing = len(repr(int(float(x)))) - len(s)
    return 0.5 * 10 ** trailing


_SP811_NUMBERS = None


def sp811_numbers():
    """Every `d.ddd ddd E+dd` factor SP 811 prints, as {value: page}."""
    global _SP811_NUMBERS
    if _SP811_NUMBERS is None:
        from pypdf import PdfReader
        out = {}
        reader = PdfReader(SP811)
        for i, page in enumerate(reader.pages, start=1):
            text = ' '.join((page.extract_text() or '').split())
            for m in re.finditer(r'(\d[\d ]*\.?[\d ]*)E\s*([+−-])\s*(\d+)', text):
                try:
                    v = float(m.group(1).replace(' ', ''))
                    e = int(m.group(3)) * (-1 if m.group(2) in '-−' else 1)
                except ValueError:
                    continue
                out.setdefault(v * 10.0 ** e, i)
        _SP811_NUMBERS = out
    return _SP811_NUMBERS


def printed_on(value, table=None):
    """The SP 811 page printing `value` (to 1e-6 relative), or None."""
    for v, page in (table if table is not None else sp811_numbers()).items():
        if v and math.isclose(v, value, rel_tol=1e-6):
            return page
    return None


def is_power_of_ten(x):
    try:
        lg = math.log10(abs(float(x)))
    except (ValueError, TypeError):
        return False
    return abs(lg - round(lg)) < 1e-12


# ---------------------------------------------------------------- P1
def p1_material_companions(mat, factor=KSI_TO_MPA):
    """[(row, E_GPa, E_ksi, implied_GPa, bound, ok)] for every dual-unit row."""
    out = []
    for row, d in mat.items():
        if not isinstance(d, dict) or 'E_GPa' not in d or 'E_ksi' not in d:
            continue
        gpa, ksi = d['E_GPa'], d['E_ksi']
        implied = ksi * factor / 1000.0
        bound = slack(gpa) + slack(ksi) * factor / 1000.0
        out.append((row, gpa, ksi, implied, bound, abs(gpa - implied) <= bound))
    return out


# ---------------------------------------------------------------- P2
def p2_scales(branches=BRANCHES, branches_dir=BRANCHES_DIR, table=None):
    """[(where, scale, verdict, detail)] for every `scale=` in every tag."""
    out = []
    for branch in branches:
        path = os.path.join(branches_dir, branch, 'constants.py')
        src = open(path, encoding='utf-8').read().replace('\r\n', '\n')
        for tag in extract_tags(src):
            _artefact, kv, _flags, _note = parse_citation(tag['body'])
            if 'scale' not in kv:
                continue
            where = f"{branch}.{tag['table']} (line {tag['line']})"
            try:
                val = float(kv['scale'])
            except ValueError:
                out.append((where, kv['scale'], 'MALFORMED', 'not a number'))
                continue
            if is_power_of_ten(val):
                out.append((where, val, 'POWER-OF-TEN', 'a decimal re-scaling; no source needed'))
                continue
            page = printed_on(val, table)
            out.append((where, val, 'PRINTED' if page else 'UNPRINTED',
                        f'SP 811 p.{page}' if page else 'no artefact on disk prints this factor'))
    return out


# ---------------------------------------------------------------- P3
def p3_specific_gravity(ranges, lo=DAS_GS_LO, hi=DAS_GS_HI):
    """[(row, value, ok)] for every endpoint of SPECIFIC_GRAVITY_RANGES."""
    out = []
    for row, span in ranges.items():
        for v in (span if isinstance(span, (tuple, list)) else (span,)):
            out.append((row, v, lo <= v <= hi))
    return out


def das_bound():
    """(lo, hi, where) read from the book, or None when it is not on this machine."""
    man = json.load(open(os.path.join(REFS, 'MANIFEST.json'), encoding='utf-8'))
    path = next((f['path'] for f in man.get('local_only_copyrighted', {}).get('files', [])
                 if 'Das, Sobhan' in f['path']), None)
    if not path or not os.path.exists(os.path.join(REPO, path)):
        return None
    from pypdf import PdfReader
    reader = PdfReader(os.path.join(REPO, path))
    for page in (63, 64):
        text = ' '.join((reader.pages[page - 1].extract_text() or '').split())
        m = re.search(r'range of (\d\.\d) to (\d\.\d)', text)
        if m:
            return float(m.group(1)), float(m.group(2)), f'PGE PDF p.{page}'
    return None


# ---------------------------------------------------------------- run
def run():
    import importlib
    failures, checks = [], 0

    print('P1  MATERIAL_PROPERTIES: E_GPa against E_ksi through SP 811 6.894757')
    mech = importlib.import_module('data.templates.branches.mechanical_engineering.constants')
    rows = p1_material_companions(mech.MATERIAL_PROPERTIES)
    for row, gpa, ksi, implied, bound, ok in rows:
        checks += 1
        if not ok:
            failures.append(f'P1 {row}: E_GPa {gpa} vs E_ksi {ksi} -> {implied:.4f} GPa, '
                            f'off by {abs(gpa - implied):.4f} > the {bound:.4f} its own '
                            f'precision allows ({100 * (gpa - implied) / implied:+.2f}%)')
    print(f'    {len(rows)} dual-unit rows, {sum(1 for r in rows if not r[5])} disagree')

    print('P2  every scale= is a printed factor or a power of ten')
    scales = p2_scales()
    for where, val, verdict, detail in scales:
        checks += 1
        if verdict in ('UNPRINTED', 'MALFORMED'):
            failures.append(f'P2 {where}: scale={val} {detail}')
    kinds = {}
    for _w, _v, verdict, _d in scales:
        kinds[verdict] = kinds.get(verdict, 0) + 1
    print(f'    {len(scales)} scale= tags: {kinds}')

    print('P3  SPECIFIC_GRAVITY_RANGES inside the bound Das states')
    civil = importlib.import_module('data.templates.branches.civil_engineering.constants')
    bound = das_bound()
    if bound is None:
        print('    UNRESOLVABLE-FROM-CLONE: Das PGE is not on this machine - not a failure')
    else:
        lo, hi, where = bound
        for row, v, ok in p3_specific_gravity(civil.SPECIFIC_GRAVITY_RANGES, lo, hi):
            checks += 1
            if not ok:
                failures.append(f'P3 {row}: Gs {v} outside Das {lo}-{hi} ({where})')
        print(f'    Das states {lo}-{hi} ({where}); '
              f'{len(civil.SPECIFIC_GRAVITY_RANGES)} rows checked')

    print(f'{checks} checks')
    for f in failures:
        print('  - ' + f)
    print('all pass' if not failures else f'{len(failures)} FAILURES')
    return 1 if failures else 0


# ---------------------------------------------------------------- self-test
_MAT = {
    'Steel': {'E_GPa': 200, 'E_ksi': 29000, 'nu': 0.30},
    'Aluminum 6061-T6': {'E_GPa': 68.9, 'E_ksi': 10000, 'nu': 0.33},
    'Nylon': {'E_GPa': 2.1, 'E_ksi': 300, 'nu': 0.39},
}
_FAKE_SP811 = {6.894757: 63, 27679.9: 66, 47.88026: 65}


def selftest():
    bad = []
    clean = [r for r in p1_material_companions(_MAT) if not r[5]]
    if clean:
        print('  the clean fixture itself fails - no plant can be judged:')
        for r in clean:
            print(f'    {r}')
        return 1

    cases = []
    # P1: two ways for a companion pair to be wrong - a decimal slip and a transposition
    slipped = dict(_MAT, Steel={'E_GPa': 20.0, 'E_ksi': 29000, 'nu': 0.30})
    cases.append(('P1 an SI value off by a factor of ten',
                  any(not r[5] and r[0] == 'Steel' for r in p1_material_companions(slipped))))
    transposed = dict(_MAT, Nylon={'E_GPa': 2.1, 'E_ksi': 30, 'nu': 0.39})
    cases.append(('P1 a transposed US value (300 -> 30)',
                  any(not r[5] and r[0] == 'Nylon' for r in p1_material_companions(transposed))))
    # and a control: a pair that disagrees by less than its own rounding slack
    rounded = dict(_MAT, Steel={'E_GPa': 199.9, 'E_ksi': 29000, 'nu': 0.30})
    cases.append(('P1 control: a pair inside its own rounding slack',
                  all(r[5] for r in p1_material_companions(rounded))))

    # P2: a factor nobody prints, and one that is a power of ten (must pass)
    import tempfile
    plants = {
        'unprinted': '# [ON-DISK] a.pdf @ page=1 scale=6.8948 precision=3sf\nX = 1\n',
        'powerten': '# [ON-DISK] a.pdf @ page=1 scale=1e3 precision=3sf\nY = 1\n',
        'printed': '# [ON-DISK] a.pdf @ page=1 scale=6.894757 precision=3sf\nZ = 1\n',
        'malformed': '# [ON-DISK] a.pdf @ page=1 scale=six precision=3sf\nW = 1\n',
    }
    with tempfile.TemporaryDirectory() as td:
        for name, src in plants.items():
            os.makedirs(os.path.join(td, name), exist_ok=True)
            with open(os.path.join(td, name, 'constants.py'), 'w', encoding='utf-8') as fh:
                fh.write('# @kind: standard\n# @units: 1\n' + src)
        got = {n: p2_scales([n], td, _FAKE_SP811)[0][2] for n in plants}
    cases += [
        ('P2 a scale no artefact prints', got['unprinted'] == 'UNPRINTED'),
        ('P2 a scale that is not a number', got['malformed'] == 'MALFORMED'),
        ('P2 control: a power of ten needs no source', got['powerten'] == 'POWER-OF-TEN'),
        ('P2 control: a printed factor is found', got['printed'] == 'PRINTED'),
    ]

    # P3: above the bound and below it
    cases += [
        ('P3 a Gs window above Das 2.9',
         any(not ok for _r, _v, ok in p3_specific_gravity({'x': (2.95, 3.1)}))),
        ('P3 a Gs window below Das 2.6',
         any(not ok for _r, _v, ok in p3_specific_gravity({'x': (2.1, 2.4)}))),
        ('P3 control: a window inside the bound',
         all(ok for _r, _v, ok in p3_specific_gravity({'x': (2.65, 2.80)}))),
    ]

    for label, ok in cases:
        print(f"  [{'ok' if ok else 'FAIL'}] {label}")
        if not ok:
            bad.append(label)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
