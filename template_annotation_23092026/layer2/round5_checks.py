"""Round 5: the two repaired chemical templates, measured.

    python -m template_annotation_23092026.layer2.round5_checks              # seeds 0-499; writes round5_checks.md/.json
    python -m template_annotation_23092026.layer2.round5_checks --seeds 50
    python -m template_annotation_23092026.layer2.round5_checks --pool       # also the frozen pool's items (local, counts only)
    python -m template_annotation_23092026.layer2.round5_checks --kit 6      # also a round's kit instances (local)
    python -m template_annotation_23092026.layer2.round5_checks --check      # self-test of the parsers and solvers

Round 5 (fixes_round5.md) changes the wording of two templates: work_isothermal_virial states its
reading (a closed system, W = -∫P dV, the volume form of the two-term virial equation, B from the
Pitzer correlation with Abbott's B0 and B1), and adiabatic_flame_temperature prints the data its
energy balance consumes. One behaviour changes besides: virial redraws a draw whose printed work lies
more than half the tolerance from the exact work for its stated data (helium, at a few kelvin).
Three measurements, at the layer-0 gate's seeds 0..N-1:

  unchanged   the current code against the pre-repair code (BASE, the last commit that touched
              either file, which is also the freeze commit), seed by seed. Virial: a seed either
              states the same substance, P1, P2, T, Tc, Pc and omega, and then its solution from
              Step 2 on must be identical, or it was redrawn by the guard, and then the pre-repair
              instance must lie more than half the tolerance from its exact work. Flame: the
              solution is identical and the question is the old question with the data block
              appended; the data block prints exactly the floats the code consumes.
  re-solved   each question solved again from its own text, at full precision, by code that shares
              nothing with the template: virial's W in SI units from the stated P1, P2, T, Tc, Pc and
              omega; flame's T_ad from the stated fuel, air and data, with the stoichiometry balanced
              from the fuel's formula and the energy balance integrated exactly. Reported: the
              largest relative distance from the gold answer, against the answer check's tolerance
              (REL = 0.002), and for flame whether the exact solution rounds to the gold kelvin.
  readings    for virial, where the readings the experts named would land: the closed system with
              the pressure-explicit form, and steady-flow shaft work, ∫V dP, under either form.
              Reported: the share of instances each puts outside the tolerance of the gold.

Writes round5_checks.md and round5_checks.json beside this file. With --pool, the frozen pool's
items of the two templates are re-solved too and only counts are written (the pool is private).
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import re
import subprocess
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.core import discover, generate, seed_all  # noqa: E402
from data.templates.branches.chemical_engineering.constants import (  # noqa: E402
    CP_PARAMS_COMBUSTION, HEATS_OF_FORMATION)

BASE = '03b1dd7f0bab084b1ac2b6640b0e9a6e9b1e31bd'   # round 4 (D-115) and the freeze commit
REL = 0.002                                           # evaluator_pilot_17092026/evaluators/answer.py
VIRIAL = ('template_work_isothermal_virial',
          'data/templates/branches/chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.py')
FLAME = ('template_adiabatic_flame_temperature',
         'data/templates/branches/chemical_engineering/thermodynamics/heat_effects.py')
POOL = REPO / 'full_run_28092026' / 'pool'
NUM = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'

# ------------------------------------------------------------------ parsing the texts

V_NEW = re.compile(rf'One mole of (?P<sub>.+?) gas in a closed system .*? from (?P<P1>{NUM}) bar to '
                   rf'(?P<P2>{NUM}) bar at a constant temperature of (?P<T>{NUM}) K\.', re.S)
V_OLD = re.compile(rf'compress 1 mole of (?P<sub>.+?) from (?P<P1>{NUM}) bar to (?P<P2>{NUM}) bar at a '
                   rf'constant temperature of (?P<T>{NUM}) K\.', re.S)
V_PROPS = re.compile(rf'Tc = (?P<Tc>{NUM}) K, Pc = (?P<Pc>{NUM}) bar, ω = (?P<w>{NUM})$')
V_GOLD = re.compile(rf'\*\*Answer:\*\*.*?\*\*(?P<W>{NUM}) J/mol\*\*')

F_HEAD = re.compile(rf'^(?P<fuel>.+?) gas enters a furnace at (?P<T0>{NUM}) K and is burned completely '
                    rf'with (?:the theoretical amount of dry air|(?P<ex>\d+)% excess dry air)')
F_AIR = re.compile(rf'Dry air is (?P<ratio>{NUM}) mol N2 per mol O2')
F_HF = re.compile(rf'Standard heats of formation at (?P<Tref>{NUM}) K, in kJ/mol \(water as vapour\): '
                  rf'(?P<list>.+)$', re.M)
F_R = re.compile(rf'R = (?P<R>{NUM}) J/\(mol·K\)')
F_CP = re.compile(rf'^- (?P<sp>\S+): A = (?P<A>{NUM}), B = (?P<B>{NUM}), C = (?P<C>{NUM}), '
                  rf'D = (?P<D>{NUM})$', re.M)
F_GOLD = re.compile(r'\*\*Answer:\*\*.*?\*\*(?P<T>\d+) K\*\*')


def glyphs(text: str) -> str:
    """The solution with its operator glyphs mapped back (round 6 prints · and ^ where it printed * and **)."""
    return text.replace('·', '*').replace('^', '**')


def parse_virial(q: str, s: str, new: bool = True) -> dict:
    m = (V_NEW if new else V_OLD).search(q)
    p = V_PROPS.search(q.strip().splitlines()[-1])
    g = V_GOLD.search(s)
    if not (m and p and g):
        raise ValueError('virial text does not parse')
    return {'sub': m['sub'], 'P1': float(m['P1']), 'P2': float(m['P2']), 'T': float(m['T']),
            'Tc': float(p['Tc']), 'Pc': float(p['Pc']), 'w': float(p['w']), 'gold': float(g['W'])}


def parse_flame(q: str, s: str) -> dict:
    h, a, hf, r, g = F_HEAD.search(q), F_AIR.search(q), F_HF.search(q), F_R.search(q), F_GOLD.search(s)
    if not (h and a and hf and r and g):
        raise ValueError('flame text does not parse')
    hf_vals = {}
    for part in hf['list'].split(', '):
        sp, val = part.split(': ')
        hf_vals[sp] = float(val)
    cp = {m['sp']: {k: float(m[k]) for k in 'ABCD'} for m in F_CP.finditer(q)}
    if not cp:
        raise ValueError('flame question carries no heat-capacity rows')
    return {'fuel': h['fuel'], 'T0': float(h['T0']), 'excess': int(h['ex'] or 0),
            'ratio': float(a['ratio']), 'Tref': float(hf['Tref']), 'hf': hf_vals, 'R': float(r['R']),
            'cp': cp, 'gold': int(g['T'])}


# ------------------------------------------------------------------ independent solvers

def virial_readings(d: dict) -> dict:
    """W on the gas, J/mol, in SI units, under the stated reading and under the other three."""
    R = 8.314                                    # J/(mol K)
    T, P1, P2 = d['T'], d['P1'] * 1e5, d['P2'] * 1e5
    Tr = T / d['Tc']
    B0 = 0.083 - 0.422 / Tr ** 1.6
    B1 = 0.139 - 0.172 / Tr ** 4.2
    B = R * d['Tc'] / (d['Pc'] * 1e5) * (B0 + d['w'] * B1)   # m3/mol
    RT = R * T

    def v_volume_form(P):                        # P V^2 - RT V - RT B = 0, the gas root
        return (RT + math.sqrt(RT * RT + 4 * P * RT * B)) / (2 * P)
    V1, V2 = v_volume_form(P1), v_volume_form(P2)
    closed_volume = -(RT * math.log(V2 / V1) - B * RT * (1 / V2 - 1 / V1))
    closed_pressure = RT * math.log(P2 / P1)     # V = RT/P + B, so -∫P dV = ∫RT/P dP
    flow_pressure = RT * math.log(P2 / P1) + B * (P2 - P1)
    flow_volume = closed_volume + B * RT * (1 / V2 - 1 / V1)   # ∫V dP = Δ(PV) - ∫P dV
    return {'stated': closed_volume, 'closed, pressure-explicit form': closed_pressure,
            'steady flow, volume form': flow_volume, 'steady flow, pressure-explicit form': flow_pressure,
            '_B': B, '_V': (V1, V2)}


def atoms(species: str) -> dict:
    formula = species.split('(')[0]
    out: dict[str, int] = {}
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def flame_solve(d: dict) -> float:
    """T_ad from the question's own data, per mole of fuel, energy balance integrated exactly."""
    products_listed = set(d['cp'])
    fuel = [sp for sp in d['hf'] if sp not in products_listed and sp not in ('O2(g)', 'N2(g)')]
    if len(fuel) != 1:
        raise ValueError(f'cannot tell the fuel among {list(d["hf"])}')
    fuel = fuel[0]
    a = atoms(fuel)
    C, H, O, N = (a.get(k, 0) for k in 'CHON')
    o2_th = C + H / 4 - O / 2
    o2_sup = o2_th * (1 + d['excess'] / 100)
    n2_air = d['ratio'] * o2_sup
    prod = {'CO2(g)': C, 'H2O(g)': H / 2, 'N2(g)': N / 2 + n2_air, 'O2(g)': o2_sup - o2_th}
    prod = {k: v for k, v in prod.items() if v > 1e-12}
    if set(prod) != products_listed:
        raise ValueError(f'balanced products {sorted(prod)} differ from the listed {sorted(products_listed)}')
    dH = (sum(n * d['hf'][sp] for sp, n in prod.items())
          - d['hf'][fuel] - o2_sup * d['hf']['O2(g)'] - n2_air * d['hf']['N2(g)']) * 1000.0   # J
    T0, R = d['T0'], d['R']

    def sensible(T):
        return sum(n * R * (c['A'] * (T - T0) + c['B'] / 2 * (T * T - T0 * T0)
                            + c['C'] / 3 * (T ** 3 - T0 ** 3) - c['D'] * (1 / T - 1 / T0))
                   for sp, n in prod.items() for c in [d['cp'][sp]])
    lo, hi = T0, 6000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if sensible(mid) + dH < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ------------------------------------------------------------------ the pre-repair code

def base_function(tid: str, path: str):
    src = subprocess.run(['git', 'show', f'{BASE}:{path}'], cwd=REPO, capture_output=True,
                         text=True, encoding='utf-8', check=True).stdout
    mod = types.ModuleType(f'_round5_base_{tid}')
    exec(compile(src, f'{BASE[:7]}:{path}', 'exec'), mod.__dict__)   # noqa: S102 - the repo's own code
    return getattr(mod, tid)


def run(seeds: int) -> dict:
    refs = {r.template_id: r for r in discover(None)}
    out = {'base': BASE, 'seeds': seeds, 'rel_tol': REL}

    # virial. A seed gives either the same problem as the pre-repair code (every stated number
    # identical), or a redraw, which is right only if the pre-repair instance lay more than half the
    # tolerance from its exact work (the helium guard) or was an organic above the cap (round 6).
    from data.templates.branches.chemical_engineering.thermodynamics.volumetric_properties_pure_fluids import (
        _VIRIAL_ORGANICS, _VIRIAL_ORGANIC_T_MAX)
    old_fn = base_function(*VIRIAL)
    same = same_tail = parsed = redrawn = redrawn_explained = old_outside = within_half = 0
    hot_organics, by_guard, by_cap = 0, 0, 0
    substances = collections.Counter()
    worst, worst_seed, outside = 0.0, None, []
    alt = {k: 0 for k in ('closed, pressure-explicit form', 'steady flow, volume form',
                          'steady flow, pressure-explicit form')}
    redrawn_subs, old_outside_subs, old_worst_by_sub = {}, {}, {}
    step1 = set()
    for s in range(seeds):
        new = generate(refs[VIRIAL[0]], s, capture=False)
        seed_all(s)
        oq, osol = old_fn()
        step1.add(new.solution.split('**Step 2:**', 1)[0])
        dn, do = parse_virial(new.question, new.solution), parse_virial(oq, osol, new=False)
        parsed += 1
        wo = virial_readings(do)['stated']
        old_rel = abs(wo - do['gold']) / abs(do['gold'])
        old_worst_by_sub[do['sub']] = max(old_worst_by_sub.get(do['sub'], 0.0), old_rel)
        if old_rel > REL:
            old_outside += 1
            old_outside_subs[do['sub']] = old_outside_subs.get(do['sub'], 0) + 1
        substances[dn['sub']] += 1
        hot_organics += dn['sub'] in _VIRIAL_ORGANICS and dn['T'] > _VIRIAL_ORGANIC_T_MAX
        if dn == do:
            same += 1
            same_tail += glyphs(new.solution.split('**Step 2:**', 1)[1]) == glyphs(osol.split('**Step 2:**', 1)[1])
        else:
            redrawn += 1
            redrawn_subs[do['sub']] = redrawn_subs.get(do['sub'], 0) + 1
            guard = abs(wo - do['gold']) > (REL / 2) * abs(wo)
            cap = do['sub'] in _VIRIAL_ORGANICS and do['T'] > _VIRIAL_ORGANIC_T_MAX
            by_guard += guard
            by_cap += cap
            redrawn_explained += guard or cap
        w = virial_readings(dn)
        rel = abs(w['stated'] - dn['gold']) / abs(dn['gold'])
        within_half += rel <= REL / 2
        if rel > worst:
            worst, worst_seed = rel, s
        if rel > REL:
            outside.append({'seed': s, 'substance': dn['sub'], 'T': dn['T'], 'gold': dn['gold'],
                            'exact': round(w['stated'], 3), 'rel': round(rel, 5)})
        for k in alt:
            alt[k] += abs(w[k] - dn['gold']) > REL * abs(dn['gold'])
    out['virial'] = {
        'generated': parsed, 'same_problem': same, 'same_problem_solution_from_step2_identical': same_tail,
        'redrawn': redrawn, 'redrawn_explained': redrawn_explained, 'redrawn_old_beyond_half_tol': by_guard,
        'redrawn_old_organic_above_cap': by_cap, 'redrawn_by_substance': redrawn_subs,
        'organic_cap_K': _VIRIAL_ORGANIC_T_MAX, 'organics_above_cap': hot_organics,
        'substances_drawn': len(substances), 'organics_drawn': sorted(x for x in substances if x in _VIRIAL_ORGANICS),
        'distinct_step1_blocks': len(step1),
        'pre_repair_outside_tol': old_outside, 'pre_repair_outside_tol_by_substance': old_outside_subs,
        'pre_repair_worst_rel_other_substances': round(max(
            (r for sub, r in old_worst_by_sub.items() if sub not in old_outside_subs), default=0.0), 6),
        'resolved_within_tol': parsed - len(outside), 'resolved_within_half_tol': within_half,
        'resolved_worst_rel': round(worst, 6), 'resolved_worst_seed': worst_seed,
        'outside_tol': outside[:20], 'other_readings_outside_tol': alt}

    # flame
    old_fn = base_function(*FLAME)
    same_sol = prefix = data_exact = parsed = rounds_same = 0
    worst_k, worst_rel, worst_seed, zero_excess = 0.0, 0.0, None, 0
    for s in range(seeds):
        new = generate(refs[FLAME[0]], s, capture=False)
        seed_all(s)
        oq, osol = old_fn()
        same_sol += new.solution == osol
        prefix += oq.endswith('in Kelvin.') and new.question.startswith(oq[:-1] + ', reported to')
        d = parse_flame(new.question, new.solution)
        parsed += 1
        zero_excess += d['excess'] == 0
        data_exact += (all(d['cp'][sp] == {k: CP_PARAMS_COMBUSTION[sp][k] for k in 'ABCD'} for sp in d['cp'])
                       and all(v == HEATS_OF_FORMATION[sp] for sp, v in d['hf'].items()))
        T = flame_solve(d)
        rounds_same += round(T) == d['gold']
        if abs(T - d['gold']) > worst_k:
            worst_k, worst_rel, worst_seed = abs(T - d['gold']), abs(T - d['gold']) / d['gold'], s
    out['flame'] = {
        'generated': parsed, 'at_zero_excess_air': zero_excess, 'solution_identical': same_sol,
        'question_is_old_plus_data_block': prefix, 'data_block_exact': data_exact,
        'resolved_rounds_to_gold': rounds_same, 'resolved_worst_K': round(worst_k, 4),
        'resolved_worst_rel': round(worst_rel, 6), 'resolved_worst_seed': worst_seed}
    return out


def run_pool() -> dict:
    """The frozen pool's items of the two templates, re-solved from their text. Counts only."""
    res = {}
    for tid, _ in (VIRIAL, FLAME):
        short = tid[len('template_'):]
        items = [json.loads(x) for p in sorted(POOL.rglob('*.jsonl'))
                 for x in p.read_text(encoding='utf-8').splitlines() if json.loads(x)['id'] == short]
        ok, worst = 0, 0.0
        for it in items:
            if tid == VIRIAL[0]:
                d = parse_virial(it['question'], it['solution'])
                rel = abs(virial_readings(d)['stated'] - d['gold']) / abs(d['gold'])
                ok += rel <= REL
            else:
                d = parse_flame(it['question'], it['solution'])
                T = flame_solve(d)
                rel = abs(T - d['gold']) / d['gold']
                ok += round(T) == d['gold']
            worst = max(worst, rel)
        res[short] = {'items': len(items), 'within': ok, 'worst_rel': round(worst, 6)}
    return res


def run_kit(round_no: int) -> dict:
    """A round's kit instances (tasks_round<N>/pool.json, local), re-solved from their text."""
    pool = json.loads((HERE / f'tasks_round{round_no}' / 'pool.json').read_text(encoding='utf8'))
    res = {'round': round_no, 'instances': 0, 'within': 0,
           'seeds': sorted({i['seed'] for x in pool.values() for i in x['instances']})}
    for item in pool.values():
        for inst in item['instances']:
            q, s = inst['question'], inst['solution']
            res['instances'] += 1
            if V_NEW.search(q):
                d = parse_virial(q, s)
                res['within'] += abs(virial_readings(d)['stated'] - d['gold']) <= REL / 2 * abs(d['gold'])
            else:
                d = parse_flame(q, s)
                res['within'] += round(flame_solve(d)) == d['gold']
    return res


def report(r: dict) -> str:
    v, f = r['virial'], r['flame']
    n = r['seeds']
    lines = [
        '# Round 5 checks: the two repaired chemical templates', '',
        f'Generated by `round5_checks.py` at seeds 0-{n - 1} (the layer-0 gate\'s seeds), against the',
        f'pre-repair code at `{r["base"][:7]}`. Answer tolerance {r["rel_tol"]} relative (answer.py `REL`).', '',
        '## work_isothermal_virial', '',
        '| check | result |', '|---|---|',
        f'| generated | {v["generated"]} of {n} |',
        f'| same problem as the pre-repair code (substance, P1, P2, T, Tc, Pc, omega identical) | '
        f'{v["same_problem"]} of {n} |',
        f'| of these, solution from Step 2 on identical but for the operator glyphs (· for *, ^ for **) | '
        f'{v["same_problem_solution_from_step2_identical"]} of {v["same_problem"]} |',
        f'| redrawn | {v["redrawn"]} of {n}, each explained by a round-5 change in {v["redrawn_explained"]}: the '
        f'pre-repair instance lay more than half the tolerance from its exact work in {v["redrawn_old_beyond_half_tol"]}, '
        f'was an organic above {v["organic_cap_K"]:.0f} K in {v["redrawn_old_organic_above_cap"]} '
        f'({", ".join(f"{k} {c}" for k, c in sorted(v["redrawn_by_substance"].items())) or "none"}) |',
        f'| organics above {v["organic_cap_K"]:.0f} K | {v["organics_above_cap"]} of {n} |',
        f'| substances drawn | {v["substances_drawn"]}; organics among them: {", ".join(v["organics_drawn"]) or "none"} |',
        f'| distinct Step 1 texts (the new reading and form) | {v["distinct_step1_blocks"]} |',
        f'| pre-repair code: exact work under the stated reading outside tolerance of the gold | '
        f'{v["pre_repair_outside_tol"]} of {n} ({", ".join(f"{k} {c}" for k, c in v["pre_repair_outside_tol_by_substance"].items()) or "none"}); '
        f'every other substance within {v["pre_repair_worst_rel_other_substances"]:.6f} |',
        f'| re-solved from the question (SI units, full precision) within tolerance of the gold | '
        f'{v["resolved_within_tol"]} of {n}; within half of it {v["resolved_within_half_tol"]} of {n}; '
        f'worst {v["resolved_worst_rel"]:.6f} (seed {v["resolved_worst_seed"]}) |',
    ]
    for k, c in v['other_readings_outside_tol'].items():
        lines.append(f'| reading "{k}" outside tolerance of the gold | {c} of {n} |')
    if v['outside_tol']:
        lines += ['', 'Instances whose exact answer under the stated reading is outside the tolerance of the gold:',
                  '', '| seed | substance | T (K) | gold (J/mol) | exact (J/mol) | relative |', '|---:|---|---:|---:|---:|---:|']
        lines += [f'| {x["seed"]} | {x["substance"]} | {x["T"]} | {x["gold"]} | {x["exact"]} | {x["rel"]} |'
                  for x in v['outside_tol']]
    lines += [
        '', '## adiabatic_flame_temperature', '',
        '| check | result |', '|---|---|',
        f'| generated | {f["generated"]} of {n} ({f["at_zero_excess_air"]} at 0% excess air) |',
        f'| solution identical to the pre-repair code | {f["solution_identical"]} of {n} |',
        f'| question is the pre-repair question with the data block appended | {f["question_is_old_plus_data_block"]} of {n} |',
        f'| data block prints exactly the Cp/R coefficients and heats of formation the code consumes | {f["data_block_exact"]} of {n} |',
        f'| re-solved from the question (stoichiometry from the formula, exact integral) rounds to the gold kelvin | '
        f'{f["resolved_rounds_to_gold"]} of {n}; worst {f["resolved_worst_K"]} K, {f["resolved_worst_rel"]:.6f} relative '
        f'(seed {f["resolved_worst_seed"]}) |',
    ]
    if 'kit' in r:
        k = r['kit']
        lines += ['', f"## The round-{k['round']} kit (tasks_round{k['round']}/, local)", '',
                  f"Instances at seeds {k['seeds'][0]}-{k['seeds'][-1]}, re-solved from their text: {k['within']} of "
                  f"{k['instances']} reach the gold (virial within half the tolerance, flame to the kelvin)."]
    if 'pool' in r:
        lines += ['', '## The frozen pool (local; counts only)', '', '| template | items | within | worst relative |',
                  '|---|---:|---:|---:|']
        lines += [f'| {k} | {x["items"]} | {x["within"]} | {x["worst_rel"]} |' for k, x in r['pool'].items()]
    return '\n'.join(lines) + '\n'


def self_test() -> int:
    fails = []
    if atoms('CH3OH(g)') != {'C': 1, 'H': 4, 'O': 1}:
        fails.append('atoms CH3OH')
    if atoms('C2H5OH(g)') != {'C': 2, 'H': 6, 'O': 1} or atoms('NH3(g)') != {'N': 1, 'H': 3}:
        fails.append('atoms C2H5OH / NH3')
    # methane, theoretical air, the table's data: the docstring's documented 2327 K
    d = {'T0': 298.15, 'excess': 0, 'ratio': 3.76, 'R': 8.314,
         'hf': {sp: HEATS_OF_FORMATION[sp] for sp in ('CH4(g)', 'CO2(g)', 'H2O(g)', 'N2(g)', 'O2(g)')},
         'cp': {sp: dict(CP_PARAMS_COMBUSTION[sp]) for sp in ('CO2(g)', 'H2O(g)', 'N2(g)')}}
    if round(flame_solve(d)) != 2327:
        fails.append(f'methane at theoretical air gives {flame_solve(d):.2f} K, not 2327')
    # virial: B = 0 is the ideal gas, where every reading gives RT ln(P2/P1)
    ideal = {'T': 400.0, 'P1': 10.0, 'P2': 30.0, 'Tc': 300.0, 'Pc': 40.0, 'w': 0.0}
    w = virial_readings(ideal)
    V1, V2 = w['_V']
    for P, V in ((10e5, V1), (30e5, V2)):
        if abs(8.314 * 400 / V * (1 + w['_B'] / V) - P) > 1e-6 * P:
            fails.append('volume-form root does not satisfy P = RT/V (1 + B/V)')
    if abs(w['steady flow, volume form'] - w['stated'] - w['_B'] * 8.314 * 400 * (1 / V2 - 1 / V1)) > 1e-9:
        fails.append('flow work is not closed work plus the change in PV')
    # parsers, on the current code at seed 0
    refs = {r.template_id: r for r in discover(None)}
    iv = generate(refs[VIRIAL[0]], 0, capture=False)
    pv = parse_virial(iv.question, iv.solution)
    if not (pv['P1'] < pv['P2'] and pv['Tc'] > 0 and pv['gold'] > 0):
        fails.append('virial parse at seed 0')
    fl = generate(refs[FLAME[0]], 0, capture=False)
    pf = parse_flame(fl.question, fl.solution)
    if round(flame_solve(pf)) != pf['gold']:
        fails.append('flame re-solve at seed 0')
    print('SELF-TEST ' + ('FAILED: ' + '; '.join(fails) if fails else 'OK'))
    return 1 if fails else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--seeds', type=int, default=500)
    ap.add_argument('--pool', action='store_true', help='also re-solve the frozen pool items (local)')
    ap.add_argument('--kit', type=int, nargs='?', const=5, default=None, metavar='ROUND',
                    help="also re-solve a round's kit instances (local); default round 5")
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        return self_test()
    r = run(a.seeds)
    if a.kit:
        r['kit'] = run_kit(a.kit)
    if a.pool:
        r['pool'] = run_pool()
    (HERE / 'round5_checks.json').write_text(json.dumps(r, indent=1, ensure_ascii=False) + '\n',
                                             encoding='utf-8', newline='\n')
    md = report(r)
    (HERE / 'round5_checks.md').write_text(md, encoding='utf-8', newline='\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print(md)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
