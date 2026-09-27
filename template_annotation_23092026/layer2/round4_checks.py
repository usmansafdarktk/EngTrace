"""Round 4 evidence: the two widened chemical templates, measured before the experts see them.

    python -m template_annotation_23092026.layer2.round4_checks     # writes round4_checks.md beside this file

heat_of_reaction_formation (D-115) draws from HESS_REACTIONS, 25 reactions: the 4 shared
REACTIONS rows plus 21, in a table only this template reads.
Checked: every reaction balances atom by atom, every species is priced in HEATS_OF_FORMATION,
the printed equation says exactly what the dicts say, and no two reactions are the same.

adiabatic_flame_temperature (D-115) now samples 0-100 percent excess air, 11 fuels x 101
levels = 1,111 cases. Its docstring and its Step 6 make claims measured at theoretical air
only, so each is re-measured over every case:
  - the six-pass iteration with its displayed intermediates against the exact fixed point,
    and whether both round to the same whole kelvin;
  - the distance of the exact fixed point from a half-kelvin rounding boundary;
  - the contraction of the iteration map at the fixed point;
  - the flame temperature against a direct NIST Shomate solve, the basis of Step 6's
    "about +/-2 K". NIST's water range starts at 500 K, so below it the lowest range is
    extrapolated, as the CP_PARAMS_COMBUSTION refit did;
  - the template itself, over 3,000 seeds, against this script's replica of its iteration.
"""
from __future__ import annotations

import json
import math
import random
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from data.templates.branches.chemical_engineering.constants import (  # noqa: E402
    COMBUSTION_REACTIONS, CP_COMBUSTION_VALID_T_MAX, CP_PARAMS_COMBUSTION,
    HEATS_OF_FORMATION, HESS_REACTIONS as REACTIONS)
from data.templates.branches.chemical_engineering.thermodynamics import heat_effects as he  # noqa: E402
from tests.constants_integrity.test_chemical_thermochemistry import _atoms  # noqa: E402

T0, R, GUESS, PASSES = 298.15, 8.314, 2000.0, 6
NIST = json.loads((REPO / 'docs' / 'references' / 'nist_webbook' / 'shomate_coefficients.json')
                  .read_text(encoding='utf-8'))
TERM = re.compile(r'^(\d+(?:\.\d+)?)?(.+)$')


def parse_side(side: str) -> dict[str, float]:
    out = {}
    for term in side.split(' + '):
        m = TERM.match(term.strip())
        out[m.group(2)] = float(m.group(1)) if m.group(1) else 1.0
    return out


def check_reactions() -> list[str]:
    lines, bad = [], []
    for rx in REACTIONS:
        left, right = Counter(), Counter()
        for sp, c in rx['reactants'].items():
            for e, n in _atoms(sp).items():
                left[e] += n * c
        for sp, c in rx['products'].items():
            for e, n in _atoms(sp).items():
                right[e] += n * c
        if any(abs(left[e] - right[e]) > 1e-9 for e in set(left) | set(right)):
            bad.append(f"{rx['name']}: does not balance")
        for sp in list(rx['reactants']) + list(rx['products']):
            if sp not in HEATS_OF_FORMATION:
                bad.append(f"{rx['name']}: {sp} has no heat of formation")
        lhs, rhs = rx['equation'].split(' → ')
        if (parse_side(lhs) != {k: float(v) for k, v in rx['reactants'].items()}
                or parse_side(rhs) != {k: float(v) for k, v in rx['products'].items()}):
            bad.append(f"{rx['name']}: printed equation differs from its dicts")
    n_eq = len({rx['equation'] for rx in REACTIONS})
    lines.append(f'| reactions | {len(REACTIONS)} |')
    lines.append(f'| distinct equations | {n_eq} |')
    lines.append(f'| unbalanced, unpriced or misprinted | {len(bad)} |')
    return lines + [f'| defect | {b} |' for b in bad]


def mean_cp_over_r(p, T):
    return (p['A'] + (p['B'] / 2) * (T + T0) + (p['C'] / 3) * (T * T + T * T0 + T0 ** 2)
            + (p['D'] / (T * T0) if p['D'] else 0.0))


def mix_cp(products, T):
    return sum(nu * R * mean_cp_over_r(CP_PARAMS_COMBUSTION[s], T) for s, nu in products.items())


def replica(rx, pct):
    """The template's own arithmetic: display-bound intermediates, six passes."""
    reactants, products, _eq, _air = he._excess_air_mixture(rx, pct)
    pe = sum(nu * HEATS_OF_FORMATION[s] for s, nu in products.items())
    re_ = sum(nu * HEATS_OF_FORMATION[s] for s, nu in reactants.items())
    q = -he._as_printed(pe - re_, '.2f') * 1000.0
    T = GUESS
    for _ in range(PASSES):
        cp = he._as_printed(mix_cp(products, T), '.2f')
        T_next = he._as_printed(T0 + q / cp, '.2f')
        step = abs(T_next - T)
        T = T_next
    return products, q, T, step


def exact_fixed_point(products, q):
    T = GUESS
    for _ in range(300):
        T = T0 + q / mix_cp(products, T)
    return T


def cp_nist(sp, T):
    ranges = NIST[sp]['ranges']
    r = next((r for r in ranges if r['T_min'] <= T <= r['T_max']),
             min(ranges, key=lambda r: r['T_min']) if T < ranges[0]['T_min'] else max(
                 ranges, key=lambda r: r['T_max']))
    t = T / 1000.0
    return r['A'] + r['B'] * t + r['C'] * t * t + r['D'] * t ** 3 + r['E'] / (t * t)


def enthalpy_tables(t_max=3200.0, h=0.5):
    """H(T) - H(T0) per species on a grid, trapezoid on Cp, J/mol."""
    grid = [T0 + i * h for i in range(int((t_max - T0) / h) + 1)]
    tables = {}
    for sp in CP_PARAMS_COMBUSTION:
        acc, prev, vals = 0.0, cp_nist(sp, grid[0]), [0.0]
        for T in grid[1:]:
            cur = cp_nist(sp, T)
            acc += 0.5 * (prev + cur) * h
            vals.append(acc)
            prev = cur
        tables[sp] = vals
    return grid, tables


def nist_flame(products, q, grid, tables):
    def total(i):
        return sum(nu * tables[s][i] for s, nu in products.items())
    lo, hi = 0, len(grid) - 1
    if total(hi) < q:
        return float('nan')
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if total(mid) < q:
            lo = mid
        else:
            hi = mid
    a, b = total(lo), total(hi)
    return grid[lo] + (q - a) / (b - a) * (grid[hi] - grid[lo])


def check_flame() -> list[str]:
    grid, tables = enthalpy_tables()
    cases, table = {}, []
    for rx in COMBUSTION_REACTIONS:
        for pct in range(101):
            products, q, T6, step = replica(rx, pct)
            Tx = exact_fixed_point(products, q)
            h = 0.01
            g = lambda T: T0 + q / mix_cp(products, T)  # noqa: E731
            slope = abs(g(Tx + h) - g(Tx - h)) / (2 * h)
            Tn = nist_flame(products, q, grid, tables)
            cases[(rx['fuel'], pct)] = round(T6)
            table.append({'fuel': rx['fuel'], 'pct': pct, 'T6': T6, 'Tx': Tx, 'step': step,
                          'slope': slope, 'Tn': Tn, 'ans': round(T6)})
    worst_fp = max(abs(r['T6'] - r['Tx']) for r in table)
    mism = [r for r in table if round(r['T6']) != round(r['Tx'])]
    ties = [r for r in table if abs(r['T6'] - math.floor(r['T6']) - 0.5) < 1e-9]
    excluded = {(r['fuel'], r['pct']) for r in mism + ties}
    kept = [r for r in table if (r['fuel'], r['pct']) not in excluded]
    margin = min(table, key=lambda r: abs(r['Tx'] - (int(r['Tx']) + 0.5)))
    dev = max(kept, key=lambda r: abs(r['ans'] - r['Tn']))
    dev0 = max((r for r in kept if r['pct'] == 0), key=lambda r: abs(r['ans'] - r['Tn']))
    lines = [
        f'| cases | {len(table)} ({len(COMBUSTION_REACTIONS)} fuels x 101 levels) |',
        f"| flame temperature range | {min(r['ans'] for r in table)} to {max(r['ans'] for r in table)} K; "
        f"validity ceiling {CP_COMBUSTION_VALID_T_MAX:.0f} K |",
        f'| six passes vs exact fixed point, worst | {worst_fp:.4f} K |',
        f'| cases where the two round differently | {len(mism)} |',
        f'| cases whose last pass displays an exact half kelvin | {len(ties)} |',
        f'| cases the template redraws, either reason | {len(excluded)} |',
        f"| last-pass change, worst (the assert allows < 0.5 K) | {max(r['step'] for r in table):.2f} K |",
        f"| closest exact fixed point to a half-kelvin boundary | {abs(margin['Tx'] - (int(margin['Tx']) + 0.5)):.4f} K "
        f"({margin['fuel']}, {margin['pct']}%) |",
        f"| contraction at the fixed point | {min(r['slope'] for r in table):.3f} to {max(r['slope'] for r in table):.3f} |",
        f"| answer vs NIST Shomate solve, worst at theoretical air | {abs(dev0['ans'] - dev0['Tn']):.2f} K "
        f"({dev0['fuel']}: {dev0['ans']} against {dev0['Tn']:.2f}) |",
        f"| answer vs NIST Shomate solve, worst over the cases the template emits | {abs(dev['ans'] - dev['Tn']):.2f} K "
        f"({dev['fuel']}, {dev['pct']}%: {dev['ans']} against {dev['Tn']:.2f}) |",
    ]
    bands = [(0, 0), (1, 25), (26, 50), (51, 75), (76, 100)]
    for lo, hi in bands:
        sub = [abs(r['ans'] - r['Tn']) for r in kept if lo <= r['pct'] <= hi]
        lines.append(f'| worst NIST deviation, {lo}-{hi}% excess | {max(sub):.2f} K |')
    agree, n, seen, hit = 0, 3000, set(), 0
    for seed in range(n):
        random.seed(seed)
        qtext, sol = he.template_adiabatic_flame_temperature()
        fuel = re.match(r'(.+?) gas enters', qtext).group(1)
        m = re.search(r'with (\d+)% excess', qtext)
        pct = int(m.group(1)) if m else 0
        ans = int(re.search(r'\*\*(\d+) K\*\*', sol.split('**Answer:**')[1]).group(1))
        agree += ans == cases[(fuel, pct)]
        hit += (fuel, pct) in excluded
        seen.add((fuel, pct))
    lines.append(f'| template answers equal to this replica, {n} seeds | {agree} of {n} '
                 f'({len(seen)} distinct cases drawn) |')
    lines.append(f'| instances drawn from a case the guard excludes, {n} seeds | {hit} |')
    return lines


def main() -> int:
    out = ['# Layer 2 round 4: the two widened chemical templates, checked', '',
           'Generated by `round4_checks.py`; the definitions are in its docstring (D-115).', '',
           '## heat_of_reaction_formation: HESS_REACTIONS', '', '| | |', '|---|---|']
    out += check_reactions()
    out += ['', '## adiabatic_flame_temperature: 0-100% excess air', '', '| | |', '|---|---|']
    out += check_flame()
    text = '\n'.join(out) + '\n'
    (HERE / 'round4_checks.md').write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
