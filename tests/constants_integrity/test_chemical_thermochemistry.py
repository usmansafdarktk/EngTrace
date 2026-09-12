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
    AIR_COMPOSITION, COMBUSTION_REACTIONS, CP_PARAMS, CP_VALID_T_MAX,
    HEATS_OF_FORMATION, REACTIONS,
)

R = 8.314462618          # CODATA 2018 molar gas constant, J/(mol K)
REF_PATH = os.path.join(REPO, 'docs', 'references', 'nist_webbook',
                        'shomate_coefficients.json')

# Two independent fits of the same thermochemistry agree to a few percent; every
# row here now sits inside 3.4% at 298 K. 8% therefore passes everything correct
# and still catches every defect found - the smallest was NO2 at 14.8%.
#
# Checked at 298-1200 K. That is NOT the full range the templates use; see the
# validity note at the end of run(), which reports where a template integrates
# past CP_VALID_T_MAX.
CP_TOL_PCT = 8.0

# Rows that cannot be checked against NIST, each for a stated reason.
CP_NO_REFERENCE = {
    'He(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Ar(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Ne(g)': 'monatomic ideal gas, Cp = 5R/2 exactly',
    'Air(g)': 'a mixture; checked against its own components instead',
}

# Heats of formation are NOT tabulated here. They are read from the on-disk
# NIST artefact at check time.
#
# Reviewer G traced three findings to the fact that this file used to carry its
# own answer key: a hardcoded DHF_REF whose values, uncertainties and even the
# set of species could drift from the artefact the citations point at. It did
# drift - C8H18, CH3OH(g) and C2H5OH(l) had reference values here and NOWHERE
# on disk, so the suite certified them against itself, and C2H6's uncertainty
# had been widened from NIST's 0.4 to 0.7, exactly enough to cover its
# deviation. A test that carries its own answers cannot fail the way it needs
# to. See DECISIONS D-034.
#
# The old check also floored every tolerance at 2 kJ/mol, which let three rows
# pass while sitting outside the uncertainty they themselves cited. There is no
# floor now: a row must land inside the stated uncertainty of a specific
# measurement NIST publishes.
#
# ROUNDING. The table is presented to 1 decimal place, so a row may sit up to
# 0.05 from the value it transcribes for that reason alone.
DHF_ROUNDING = 0.05

AIR = AIR_COMPOSITION     # declared in constants.py, not duplicated here
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
    for k, tbl in entry.items():
        if k.startswith('cp_table_') and isinstance(tbl, dict):
            for tk, v in tbl.items():
                if abs(float(tk) - T) < 0.5:
                    return v
    return None


def dhf_measurements(entry):
    """Every heat of formation the artefact records for this species.

    Returns dicts of {value, uncertainty, source}. NIST lists several
    independent determinations for most hydrocarbons and they differ by more
    than any one of them claims, so "the NIST value" is not well defined - a
    constants row is checked against the whole set and must match one.
    """
    for field, ms in entry.items():
        if field.startswith('dHf_') and field.endswith('_measurements'):
            return [m for m in ms if m.get('value') is not None]
    out = []
    for k in ('dHf_gas_kJ_per_mol', 'dHf_liquid_kJ_per_mol',
              'dHf_solid_kJ_per_mol'):
        if k in entry:
            out.append({'value': entry[k],
                        'uncertainty': entry.get('dHf_uncertainty'),
                        'source': entry.get('dHf_source', '(sole listed)')})
            break
    if 'dHf_alternative' in entry:
        out.append({'value': entry['dHf_alternative'], 'uncertainty': None,
                    'source': entry.get('dHf_alternative_source', '(alt)')})
    return out


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

    # -- how thin is the evidence behind each validity ceiling? -------------
    # G-8: C3H8 and C4H10 are NIST-checked at one temperature (298.15 K, the
    # only thing NIST publishes for them) yet CP_VALID_T_MAX licenses 1500 K.
    # That is one point licensing a 1200 K extrapolation. Reported, not failed:
    # the coefficients come from Smith-Van Ness, which states the range - but
    # the range is the SOURCE's claim, not something this repo verified.
    single_point = []
    for sp, row in CP_PARAMS.items():
        if sp in CP_NO_REFERENCE:
            continue
        entry = ref.get(sp) or {}
        pts = set()
        for r in entry.get('ranges', []):
            pts.add((r['T_min'], r['T_max']))
        for k, tbl in entry.items():
            if k.startswith('cp_table_') and isinstance(tbl, dict):
                pts.update(tbl)
        if not pts and any(k.startswith('cp_298') for k in entry):
            single_point.append((sp, CP_VALID_T_MAX.get(sp)))
    if single_point:
        print('  NOTE: verified at 298.15 K only, but CP_VALID_T_MAX licenses '
              'far beyond it (source-stated range, not verified here): '
              + ', '.join(f'{sp} to {t:.0f} K' for sp, t in single_point))

    # -- HEATS_OF_FORMATION, against the on-disk artefact -------------------
    for sp, v in HEATS_OF_FORMATION.items():
        entry = ref.get(sp)
        if entry is None:
            failures.append(f'HEATS_OF_FORMATION[{sp}] has no entry in '
                            f'{os.path.basename(REF_PATH)} - its citation '
                            f'resolves to nothing')
            continue
        cands = dhf_measurements(entry)
        if not cands:
            failures.append(f'HEATS_OF_FORMATION[{sp}] has a reference entry '
                            f'but no heat of formation in it')
            continue
        checks += 1
        best = None
        for m in cands:
            tol = max(m.get('uncertainty') or 0.0, DHF_ROUNDING)
            dev = abs(v - m['value'])
            if dev <= tol and (best is None or dev < best[0]):
                best = (dev, m)
        if best is None:
            listed = ', '.join(
                f"{m['value']}+/-{m.get('uncertainty')} ({m.get('source', '?')})"
                for m in cands)
            failures.append(
                f'HEATS_OF_FORMATION[{sp}] = {v} is outside the stated '
                f'uncertainty of every measurement NIST lists: {listed}')

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

    # -- F-2: does any consuming template integrate Cp past its validity? --
    #
    # This check exists because the suite was previously green while never
    # looking above 1200 K, and template_adiabatic_flame_temperature integrates
    # to ~2900 K. A check that stops short of where the data is actually used
    # reports on a range nobody consumes.
    import random
    import re as _re
    import data.templates.branches.chemical_engineering.thermodynamics.heat_effects as _he

    reached = 0.0
    for s_ in range(200):
        random.seed(s_)
        try:
            q, sol = _he.template_adiabatic_flame_temperature()
        except Exception:
            continue
        i = sol.find('**Answer')
        for tok in _re.findall(r'([\d,]+\.?\d*)\s*K', sol[i:] if i >= 0 else sol):
            reached = max(reached, float(tok.replace(',', '')))
    # C3.2 REPAIR. This NOTE measured the wrong table. It compared `reached` with
    # CP_VALID_T_MAX (1500 K) and computed the error from CP_PARAMS - but since D-036
    # split the table by consumer, this template reads CP_PARAMS_COMBUSTION, whose
    # ceiling is CP_COMBUSTION_VALID_T_MAX (3000 K). So it reported a template as
    # ~90% past a validity limit that has not applied to it since D-036, and quoted a
    # Cp error from polynomials the template no longer uses. The check now measures the
    # table the template actually reads, and fails rather than notes: being outside the
    # range you were refitted for is a defect, not a footnote.
    from data.templates.branches.chemical_engineering.constants import (
        CP_COMBUSTION_VALID_T_MAX, CP_PARAMS_COMBUSTION)
    checks += 1
    prod_max = CP_COMBUSTION_VALID_T_MAX
    worst = 0.0
    for sp in ('CO2(g)', 'H2O(g)', 'N2(g)', 'O2(g)'):
        e = ref.get(sp)
        n = cp_nist(e, reached) if e else None
        if n:
            worst = max(worst, abs(100.0 * (cp_svn(CP_PARAMS_COMBUSTION[sp], reached) - n) / n))
    if reached > prod_max:
        over = 100.0 * (reached - prod_max) / prod_max
        failures.append(f'adiabatic_flame_temperature reaches {reached:.0f} K, '
                        f'{over:.1f}% past the {prod_max:.0f} K validity of '
                        f'CP_PARAMS_COMBUSTION, the table it reads (D-036)')
    else:
        print(f'  F-2: adiabatic_flame_temperature reaches {reached:.0f} K, inside the '
              f'{prod_max:.0f} K validity of CP_PARAMS_COMBUSTION (D-036); worst Cp '
              f'error against NIST there is {worst:.1f}%')

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
