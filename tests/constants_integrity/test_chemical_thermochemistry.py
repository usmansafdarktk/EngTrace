"""C2.2 — physical-plausibility suite for the chemical thermochemistry tables.

Guards the three tables Phase C2 re-grounded — `CP_PARAMS`,
`HEATS_OF_FORMATION` and `COMBUSTION_REACTIONS` — against the four defect
classes that were actually found in them, so none can return silently:

  Class 1  a 10x exponent slip, from transcribing the source's column heading
           (10^3 B, 10^6 C, 10^-5 D) as part of the value
  Class 2  a row whose coefficients are simply wrong (O2, C2H2, C2H5OH, C3H6O)
  Class 3  a column shift (NO2: D's mantissa written into C)
  Class 4  gas-phase coefficients under a liquid key (benzene, toluene, hexane,
           acetone)

Every check compares against `docs/references/nist_webbook/`, which holds NIST
Chemistry WebBook data retrieved 2026-09-06. NIST publishes the **Shomate**
form; the constants use the **Smith-Van Ness** form. They are different fits to
the same underlying thermochemistry, so agreement between them is evidence
rather than a tautology - and the tolerances below are set by how well two good
fits of the same data can be expected to agree, not by what the current numbers
happen to give.

    python -m tests.constants_integrity.test_chemical_thermochemistry
"""
from __future__ import annotations

import json
import math
import os
import re
import sys
from collections import Counter

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from data.templates.branches.chemical_engineering.constants import (  # noqa: E402
    COMBUSTION_REACTIONS, CP_PARAMS, HEATS_OF_FORMATION, REACTIONS,
)

R = 8.314462618          # CODATA 2018 molar gas constant, J/(mol K)
REF_PATH = os.path.join(REPO, 'docs', 'references', 'nist_webbook',
                        'shomate_coefficients.json')

# Two independent fits of the same thermochemistry agree to a few percent. The
# gases sit inside 3.4%; the one row above that is NaCl(s) at 5.6%, a genuine
# disagreement between the source fits rather than a transcription defect
# (see the phase summary). 8% therefore passes everything correct and still
# catches every defect found: the smallest was NO2 at 14.8%.
CP_TOL_PCT = 8.0

# Rows that cannot be checked against NIST, each for a stated reason.
CP_NO_REFERENCE = {
    'He(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Ar(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Ne(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Air(g)': 'a mixture; checked against its own components instead',
    'H2SO4(l)': '[UNVERIFIED] NIST condensed-phase data is paywalled',
    'CaCO3(s)': '[UNVERIFIED] NIST free tier carries no Cp for calcite',
}

# NIST / CODATA standard enthalpies of formation, kJ/mol, with the stated
# uncertainty. Where sources differ by more than their uncertainty the tolerance
# is floored at 2 kJ/mol - immaterial here, since a heat of reaction is
# thousands of kJ/mol.
DHF_REF = {
    'CH4(g)': (-74.6, 0.3), 'C2H2(g)': (226.73, 0.8), 'C2H6(g)': (-84.0, 0.7),
    'C3H8(g)': (-104.7, 0.5), 'C4H10(g)': (-125.6, 0.67),
    'C8H18(g)': (-208.4, 0.67), 'C6H6(l)': (49.0, 0.9),
    'CH3OH(l)': (-239.5, 0.2), 'CH3OH(g)': (-205.0, 10.0),
    'C2H5OH(l)': (-276.0, 2.0), 'C2H5OH(g)': (-234.0, 2.0),
    'O2(g)': (0.0, 0.0), 'H2(g)': (0.0, 0.0), 'N2(g)': (0.0, 0.0),
    'CO(g)': (-110.53, 0.17), 'CO2(g)': (-393.51, 0.13),
    'H2O(g)': (-241.826, 0.04), 'H2O(l)': (-285.83, 0.04),
    'NH3(g)': (-45.94, 0.35), 'NO(g)': (90.29, 0.2), 'NO2(g)': (33.10, 0.2),
}
DHF_TOL_FLOOR = 2.0

# Dry air, mole fractions.
AIR = {'N2(g)': 0.78084, 'O2(g)': 0.20946, 'Ar(g)': 0.00934}
N2_PER_O2 = 3.76          # theoretical air


def cp_svn(row, T):
    """Cp from the Smith-Van Ness form the constants use, J/(mol K)."""
    return R * (row['A'] + row['B'] * T + row['C'] * T * T
                + (row['D'] / (T * T) if row['D'] else 0.0))


def cp_nist(entry, T):
    """Cp from NIST's Shomate form, or a tabulated value at 298.15 K."""
    for r in entry.get('ranges', []):
        if r['T_min'] <= T <= r['T_max']:
            t = T / 1000.0
            return (r['A'] + r['B'] * t + r['C'] * t * t + r['D'] * t ** 3
                    + r['E'] / (t * t))
    if abs(T - 298.15) < 0.5:
        for k in ('cp_298_tabulated', 'cp_298_liquid_tabulated'):
            if k in entry:
                return entry[k]
    for k in ('cp_table_TRC1997', 'cp_table_Chao1986'):
        if k in entry:
            for tk, v in entry[k].items():
                if abs(float(tk) - T) < 0.5:
                    return v
    return None


def _atoms(species):
    out = Counter()
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', species.split('(')[0]):
        if el:
            out[el] += int(n or 1)
    return out


def run():
    ref = json.load(open(REF_PATH, encoding='utf-8'))
    failures = []
    checks = 0

    # -- CP_PARAMS against NIST ------------------------------------------
    for sp, row in CP_PARAMS.items():
        if sp in CP_NO_REFERENCE:
            continue
        entry = ref.get(sp)
        if entry is None:
            failures.append(f'CP_PARAMS[{sp}] has no NIST reference and is not '
                            f'listed in CP_NO_REFERENCE with a reason')
            continue
        compared = False
        for T in (298.15, 500.0, 800.0, 1200.0):
            n = cp_nist(entry, T)
            if n is None:
                continue
            compared = True
            checks += 1
            err = 100.0 * (cp_svn(row, T) - n) / n
            if abs(err) > CP_TOL_PCT:
                failures.append(
                    f'CP_PARAMS[{sp}] Cp({T:.0f}K) = {cp_svn(row, T):.2f} vs '
                    f'NIST {n:.2f} ({err:+.1f}%, tol {CP_TOL_PCT}%)')
        if not compared:
            failures.append(f'CP_PARAMS[{sp}] reference carries no usable Cp')

    # -- monatomic gases are exact ----------------------------------------
    for sp in ('He(g)', 'Ar(g)', 'Ne(g)'):
        checks += 1
        cp = cp_svn(CP_PARAMS[sp], 298.15)
        if abs(cp - 2.5 * R) > 1e-9:
            failures.append(f'CP_PARAMS[{sp}] = {cp} but a monatomic ideal gas '
                            f'is exactly 5R/2 = {2.5 * R}')

    # -- air against its own components ------------------------------------
    for T in (298.15, 600.0, 1000.0):
        checks += 1
        mix = sum(x * cp_svn(CP_PARAMS[s], T) for s, x in AIR.items())
        err = 100.0 * (cp_svn(CP_PARAMS['Air(g)'], T) - mix) / mix
        if abs(err) > 2.0:
            failures.append(f'CP_PARAMS[Air(g)] at {T:.0f}K is {err:+.1f}% from '
                            f'the mole-weighted N2/O2/Ar average')

    # -- Class 4 guard: a liquid row must not hold gas-phase values --------
    for sp in [s for s in CP_PARAMS if s.endswith('(l)')]:
        entry = ref.get(sp)
        if entry and 'cp_298_liquid_tabulated' in entry:
            checks += 1
            n = entry['cp_298_liquid_tabulated']
            err = 100.0 * (cp_svn(CP_PARAMS[sp], 298.15) - n) / n
            if abs(err) > CP_TOL_PCT:
                failures.append(
                    f'CP_PARAMS[{sp}] is keyed as a liquid but its Cp(298) is '
                    f'{err:+.1f}% from the NIST LIQUID value - the row may hold '
                    f'gas-phase coefficients')

    # -- HEATS_OF_FORMATION -------------------------------------------------
    for sp, v in HEATS_OF_FORMATION.items():
        if sp not in DHF_REF:
            failures.append(f'HEATS_OF_FORMATION[{sp}] has no reference value')
            continue
        checks += 1
        r, u = DHF_REF[sp]
        tol = max(u, DHF_TOL_FLOOR)
        if abs(v - r) > tol:
            failures.append(f'HEATS_OF_FORMATION[{sp}] = {v} vs NIST {r} '
                            f'(diff {v - r:+.2f}, tol {tol})')

    # -- reaction stoichiometry --------------------------------------------
    for rxns, label in ((COMBUSTION_REACTIONS, 'COMBUSTION_REACTIONS'),
                        (REACTIONS, 'REACTIONS')):
        for rx in rxns:
            checks += 1
            left, right = Counter(), Counter()
            for sp, c in rx['reactants'].items():
                for e, n in _atoms(sp).items():
                    left[e] += n * c
            for sp, c in rx['products'].items():
                for e, n in _atoms(sp).items():
                    right[e] += n * c
            for e in set(left) | set(right):
                if abs(left[e] - right[e]) > 1e-9:
                    failures.append(
                        f'{label} "{rx["name"]}" does not balance for {e}: '
                        f'{left[e]} -> {right[e]}')
            n2, o2 = rx['reactants'].get('N2(g)'), rx['reactants'].get('O2(g)')
            if n2 and o2 and abs(n2 / o2 - N2_PER_O2) > 0.02:
                failures.append(
                    f'{label} "{rx["name"]}" has N2/O2 = {n2 / o2:.3f}, not the '
                    f'theoretical-air {N2_PER_O2}')

    # -- every species a reaction names must be priced ---------------------
    for rx in COMBUSTION_REACTIONS + REACTIONS:
        for sp in list(rx['reactants']) + list(rx['products']):
            checks += 1
            if sp not in HEATS_OF_FORMATION:
                failures.append(f'"{rx["name"]}" names {sp}, which has no heat '
                                f'of formation')

    print(f'{checks} checks')
    if failures:
        print(f'{len(failures)} FAILED:')
        for f in failures:
            print(f'  - {f}')
        return 1
    print('all pass')
    return 0


if __name__ == '__main__':
    sys.exit(run())
