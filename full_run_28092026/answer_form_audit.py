"""Three readings of the answer check that the full run's "incorrect" verdicts point at, measured
before any change, on the gold, the pilot's 300 expert-labelled traces and every full-run trace.

    python -m full_run_28092026.answer_form_audit                   # FREE: writes ANSWER_FORM_AUDIT.md beside this file
    python -m full_run_28092026.answer_form_audit --models kimi-k3  # some models only; writes nothing
    python -m full_run_28092026.answer_form_audit --fix             # FREE: answer.py at PRE_FIX against POST_FIX, both from
                                                                    # git; writes ANSWER_FORM_FIX.md (the record of the change)

The readings below were the design; they are applied in memory to answer.py as it stands at the time of
the run. Once the change is committed, `--fix` measures the shipped code between two pinned commits the
way parser_fix.py and digit_fix.py do, and ANSWER_FORM_FIX.md is the record the paper cites. Readings P,
A (at unit scale) and N were adopted (D-169); C was not.

WHY. Twelve "incorrect" verdicts drawn at random from the top five models' answered traces were read by
hand on 2026-10-02. Six were right answers in another form: two wrote the numeric value of a symbolic
gold (`0.583 * Q(0.968)` as 0.0971), three wrote a pi-fraction the check misreads, one answered in
radians where the gold's line ends in degrees. Across the top six models 269 answered traces are
incorrect, and 52 of them sit on `continuous_to_discrete_conversion`, 36 on `ber_estimation_mary`, 26 on
`decimation_aliasing_analysis` and 13 on `composite_shafts_series`. None of these forms occurs in the
pilot's 15 templates, so the experts never arbitrated them. Nothing here changes answer.py.

READING P   pi-fractions. `answer.PI_EXPR` reads `841*pi` and stops at the parenthesis of `(841*pi)/2447`,
            so the gold's target becomes 841*pi, a value the template also computes, so the milestone
            filter keeps it, instead of 841*pi/2447; and LaTeX's `\\pi` is never joined to its
            coefficient, so `\\dfrac{841\\pi}{2447}` reads as pi. Under P the text is normalised before
            `values` reads it: `\\pi` and the symbol become `pi`, `\\frac{A}{B}` becomes `A/B`, and
            `(a*pi)/b` and `(a/b)*pi` become `a*pi/b`, which the existing rule reads as one value.
            Applied to the gold's targets and to the trace alike.
READING A   alternative renderings, on top of P. A scalar gold line that states one quantity twice
            ("0.088960 radians, or 5.0970 degrees") yields one target, the last; a trace that answers in
            the other unit is marked wrong. Under A a scalar item is also correct when the trace matches
            another computed value the gold's answer line states that is the SAME quantity in another
            unit: its ratio to the target is 180/pi, 2*pi, a power of ten, or their inverses. A first
            version accepted any computed value on the line, the `any` mode of the phase-4 comparator
            bindings; read on the full run it credited traces right on one of two different quantities
            (`vdw_solve_for_volume`'s liquid and vapour volumes, `levenspiel_plot_interpretation`'s two
            reactor volumes, both typed scalar) and an array stated in part (`finite_convolution`), so it
            was narrowed to renderings. Those templates appear below instead, in the list of scalar and
            array items whose gold line states more than one computed value: the check verifies only the
            last of them, so a trace right on the other alone is incorrect and one right on the last alone
            is correct. The rendering is matched at unit scale, without the check's unit factors: with
            them, a cut-off Gemma trace holding `pi d^4/32` and no answer matched 0.5315 rad through
            32 +- 1 at the 1/60 factor.
READING C   not adopted, measured for the record: a number that is pi's coefficient (`0.7995 pi`) is not
            a standalone value. It removes the one wrong credit P gives (gpt-5.4-mini writes
            "omega = 0.7995 pi = 2.511 rad/sample" for 0.7987 rad/sample and is credited through the bare
            0.7995) and two accidental credits of wrong answers, but it also strips the accidental credit
            of four right answers on `cd_dc_system_analysis` (a phase written -0.75 pi for a 0.0075 s
            shift) and of two on `continuous_to_discrete_conversion#10`.
READING N   on top of P+A: a pi-fraction a*pi/b on the gold's answer line whose numerator and denominator
            are both targets is one target, the fraction's value. Eight of `continuous_to_discrete_conversion`'s
            fifteen items carry only `num` and `den` as milestones, so their targets are two integers and
            a trace stating the right omega as a decimal or an unreduced fraction gets no credit.
READING S   symbolic answers evaluated numerically. Counted, not ruled (D-138: symbolic answers are
            scored by the numbers they state): on `ber_estimation_mary` the gold is `a * Q(b)`; how many
            of its incorrect and partial verdicts state a*Q(b) itself, within 1% or one last digit.

Each reading is measured as D-137's was: the gold must stay 2,250 of 2,250 correct; on the pilot no
verdict should move away from the experts; on the full run the changed verdicts are counted per model
and per template, and written with their answer segments to scores/answer_form_audit_changes.jsonl,
gitignored with the traces, to be read. The current verdict is recomputed beside the store's label as a
check that the audit reads the same traces the store scored.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import itertools
import json
import math
import re
import subprocess
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402
import milestones  # noqa: E402

from full_run_28092026 import score  # noqa: E402

CHANGES = score.SCORES / 'answer_form_audit_changes.jsonl'
BER_TEMPLATE = 'template_ber_estimation_mary'
BER = re.compile(r'(\d+\.\d+)\s*\*\s*Q\(\s*(\d+\.\d+)\s*\)')
FRAC = re.compile(r'\\[dt]?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}')
PAREN_PI = re.compile(r'\(\s*(\d+(?:\.\d+)?)\s*\*?\s*pi\s*\)\s*/\s*\(?\s*(\d+(?:\.\d+)?)\s*\)?')
FRAC_PI = re.compile(r'\(\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)\s*\)\s*\*?\s*pi\b')


def normalise(t: str) -> str:
    """Reading P's text normalisation: pi-fractions in one shape the existing rule reads."""
    t = t.replace('\\pi', ' pi ').replace('\u03c0', ' pi ')
    t = FRAC.sub(r' \1/\2 ', t)
    t = PAREN_PI.sub(r' \1*pi/\2 ', t)
    t = FRAC_PI.sub(r' \1*pi/\2 ', t)
    return t


COEF = re.compile(r'(?<![\w.])(\d+(?:\.\d+)?)\s*\*?\s*pi\b')


def drop_coefficients(seg: str, vals: list) -> list:
    """Reading C: a number that is pi's coefficient (`0.7995 pi`) states the value 0.7995*pi, which
    PI_EXPR reads, and not 0.7995 as well. One (value, ulp) entry is removed per coefficient found."""
    out = list(vals)
    for m in COEF.finditer(seg):
        v, u = float(m.group(1)), answer._ulp(m.group(1))
        for k, (x, xu) in enumerate(out):
            if x == v and xu == u:
                del out[k]
                break
    return out


PIFRAC = re.compile(r'(?<![\w.])(\d+(?:\.\d+)?)\s*\*?\s*pi\s*/\s*(\d+(?:\.\d+)?)')


def merge_pifrac_targets(P, want: dict, item: dict) -> dict:
    """Reading N: on the gold's answer line a pi-fraction a*pi/b whose numerator and denominator are both
    targets (both computed milestones, as `continuous_to_discrete_conversion` exposes the reduced
    fraction's `num` and `den`) is one quantity, so the pair of integers is replaced by the fraction's
    value. A trace stating that value as a decimal or as an unreduced fraction then matches."""
    nums = list(want['numbers'])
    for m in PIFRAC.finditer(normalise(P.segment(item['solution']))):
        a, b = float(m.group(1)), float(m.group(2))
        ia = [k for k, (v, _) in enumerate(nums) if milestones.close(v, a, 1e-9)]
        ib = [k for k, (v, _) in enumerate(nums) if milestones.close(v, b, 1e-9)]
        if ia and ib and ia[0] != ib[0]:
            nums = [x for k, x in enumerate(nums) if k not in (ia[0], ib[0])] + [(a * math.pi / b, 0.0)]
    return {**want, 'numbers': nums}


def variant_module(coef: bool = False, pifrac: bool = False):
    """answer.py as it stands, loaded as a second module whose `values` normalises first (reading P);
    with `coef`, pi's coefficients are not standalone values either (reading C); with `pifrac`, a
    pi-fraction of two targets becomes one target (reading N)."""
    name = 'answer_P' + ('C' if coef else '') + ('N' if pifrac else '')
    src = (EVAL / 'answer.py').read_text(encoding='utf-8')
    mod = types.ModuleType(name)
    mod.__file__ = str(EVAL / 'answer.py')
    exec(compile(src, name + '.py', 'exec'), mod.__dict__)
    orig = mod.values
    if coef:
        mod.values = lambda seg: drop_coefficients(normalise(seg), orig(normalise(seg)))
    else:
        mod.values = lambda seg: orig(normalise(seg))
    if pifrac:
        orig_targets = mod.targets
        mod.targets = lambda item, ms=None: merge_pifrac_targets(mod, orig_targets(item, ms), item)
    return mod


RENDERINGS = (180.0 / math.pi, 2.0 * math.pi) + tuple(10.0 ** k for k in range(1, 13))


def same_quantity(a: float, b: float) -> bool:
    """Are a and b one quantity in two units: related by 180/pi, 2*pi or a power of ten, either way?"""
    if not a or not b:
        return False
    r = abs(a / b)
    return any(abs(r - f) <= 1e-3 * f or abs(r - 1 / f) <= 1e-3 / f for f in RENDERINGS)


def computed_on_line(P, item: dict, vals: tuple) -> list[tuple[float, float]]:
    """The distinct computed values the gold's answer line states, in the order `values` reads them."""
    comps = []
    for v, u in P.values(P.segment(item['solution'])):
        if milestones.scaled_match(v, list(vals), milestones.DISPLAY_TOL) is None:
            continue
        if any(milestones.close(v, x, 1e-9) for x, _ in comps):
            continue
        comps.append((v, u))
    return comps


def match_at_unit_scale(P, gold: float, have: list, exact: bool, gold_ulp: float) -> bool:
    """`answer.match` with no unit factor: the rendering already is the other unit, so a trace that
    answers in it states the value as it is. With the factors, a cut-off trace holding `pi d^4/32` and
    no answer matched 0.5315 rad through 32 +- 1 at the 1/60 factor (`composite_shafts_series#16`)."""
    for v, u in have:
        if exact:
            if abs(v - gold) <= 1e-9 * max(1.0, abs(gold)):
                return True
        elif abs(v - gold) <= max(P.REL * abs(gold), u if abs(v) > u else 0.0, gold_ulp) * (1 + P.SLACK):
            return True
    return False


def verdict_pa(P, text: str, item: dict, vals: tuple) -> str:
    """Reading A on top of P: a scalar item is correct when the trace matches another computed value the
    gold's answer line states that is the target in another unit, matched at unit scale."""
    lab, det = P.verdict(text, item, None, vals)
    if lab == 'correct' or item.get('answer_type') != 'scalar' or not det['targets']['numbers']:
        return lab
    target = det['targets']['numbers'][-1][0]
    comps = [(v, u) for v, u in computed_on_line(P, item, vals)
             if not milestones.close(v, target, 1e-9) and same_quantity(v, target)]
    if not comps:
        return lab
    have = P.values(P.segment(text))
    absolute = [(abs(v), u) for v, u in have]
    exact = det['targets']['exact']
    for v, gu in comps:
        if match_at_unit_scale(P, v, have, exact, gu) or match_at_unit_scale(P, abs(v), absolute, exact, gu):
            return 'correct'
    return lab


def multi_value(P, items: dict, ms_all: dict) -> dict:
    """Scalar and array templates whose gold answer line states more than one computed value, and on how
    many of their items every pair of values is one quantity in two units."""
    out = {}
    for iid, it in items.items():
        if it['answer_type'] not in ('scalar', 'array'):
            continue
        comps = computed_on_line(P, it, tuple(m['value'] for m in ms_all[iid]))
        if len(comps) < 2:
            continue
        r = out.setdefault(it['template_id'], {'type': it['answer_type'], 'items': 0, 'max': 0, 'renderings': 0})
        r['items'] += 1
        r['max'] = max(r['max'], len(comps))
        r['renderings'] += all(same_quantity(a, b) for (a, _), (b, _) in itertools.combinations(comps, 2))
    return out


def ber_numeric(item: dict, text: str) -> bool | None:
    """Reading S: does the trace state the numeric value of the gold's `a * Q(b)`? None when the gold
    line is not of that form."""
    m = BER.search(answer.segment(item['solution']))
    if not m:
        return None
    a, b = float(m.group(1)), float(m.group(2))
    val = a * 0.5 * math.erfc(b / math.sqrt(2))
    return any(abs(v - val) <= max(0.01 * val, u) for v, u in answer.values(answer.segment(text)))


# ------------------------------------------------------------------ gold and pilot

def gold(P, items: dict, ms_all: dict) -> collections.Counter:
    bad = collections.Counter()
    for iid, it in items.items():
        vals = tuple(m['value'] for m in ms_all[iid])
        if P.verdict(it['solution'], it, None, vals)[0] != 'correct':
            bad[it['template_id']] += 1
    return bad


def pilot(P, PC, PN) -> dict:
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, items, texts = V.pilot_inputs()
    labs = {tag: {} for tag in ('current',) + READINGS}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        text, item = texts[key], items[key[1]]
        vals = tuple(m['value'] for m in t['milestones'])
        labs['current'][code] = answer.verdict(text, item, None, vals)[0]
        p = P.verdict(text, item, None, vals)[0]
        labs['P'][code] = p
        labs['P+A'][code] = verdict_pa(P, text, item, vals) if p != 'correct' else 'correct'
        labs['P+A+C'][code] = verdict_pa(PC, text, item, vals)
        labs['P+A+N'][code] = verdict_pa(PN, text, item, vals)
    codes = sorted(labs['current'])
    exp = {c: truth[c]['final_answer'] for c in codes}
    npc = [c for c in codes if exp[c] != 'partial']
    res = {'traces': len(codes)}
    for tag, L in labs.items():
        moved = [c for c in codes if L[c] != labs['current'][c]]
        res[tag] = {'non_partial': sum((L[c] == 'correct') == (exp[c] == 'correct') for c in npc) / len(npc),
                    'three_way': sum(L[c] == exp[c] for c in codes) / len(codes),
                    'moved': len(moved), 'to_experts': sum(L[c] == exp[c] for c in moved),
                    'from_experts': sum(labs['current'][c] == exp[c] for c in moved),
                    'transitions': collections.Counter(f"{labs['current'][c]} -> {L[c]} (experts: {exp[c]})"
                                                       for c in moved)}
    return res


# ------------------------------------------------------------------ the full run

READINGS = ('P', 'P+A', 'P+A+C', 'P+A+N')


def _model(key: str):
    P = variant_module()
    PC = variant_module(coef=True)
    PN = variant_module(pifrac=True)
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    store = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                           .read_text(encoding='utf-8').splitlines())}
    tr = {tag: collections.Counter() for tag in READINGS}
    by_t = {tag: collections.Counter() for tag in READINGS}
    to_correct_t = {tag: collections.Counter() for tag in READINGS}
    c_vs_pa = collections.Counter()
    sums = collections.Counter()
    ber = collections.Counter()
    changes = []
    n = mismatch = 0
    for row in score.trace_rows('main', key):
        s = store[row['item_id']]
        n += 1
        cur = 'unusable' if s['unusable'] else s['answer']['label']
        sums['current'] += score.SCORE_OF.get(cur, 0.0)
        if cur == 'incorrect':
            sums['incorrect_now'] += 1
        if s['unusable']:
            continue
        item = items[row['item_id']]
        vals = tuple(m['value'] for m in ms_all[row['item_id']])
        text = row.get('text') or ''
        if answer.verdict(text, item, None, vals)[0] != cur:
            mismatch += 1
        p = P.verdict(text, item, None, vals)[0]
        pa = verdict_pa(P, text, item, vals) if p != 'correct' else 'correct'
        pac = verdict_pa(PC, text, item, vals)
        pan = verdict_pa(PN, text, item, vals)
        for tag, lab in (('P', p), ('P+A', pa), ('P+A+C', pac), ('P+A+N', pan)):
            if lab != cur:
                tr[tag][(cur, lab)] += 1
                by_t[tag][item['template_id']] += 1
                if lab == 'correct':
                    to_correct_t[tag][item['template_id']] += 1
        if pac != pa:
            c_vs_pa[('C', item['template_id'], pa, pac)] += 1
        if pan != pa:
            c_vs_pa[('N', item['template_id'], pa, pan)] += 1
        if pa != cur or pac != cur or pan != cur:
            changes.append({'model_key': key, 'item_id': row['item_id'], 'current': cur, 'P': p, 'P+A': pa,
                            'P+A+C': pac, 'P+A+N': pan, 'segment': answer.segment(text)})
        if item['template_id'] == BER_TEMPLATE and cur in ('incorrect', 'partial'):
            hit = ber_numeric(item, text)
            if hit is not None:
                ber['verdicts'] += 1
                ber['numeric'] += hit
    res = {'traces': n, 'mismatch': mismatch, 'incorrect_now': sums['incorrect_now'],
           'score_current': sums['current'] / n, 'ber': dict(ber),
           'c_vs_pa': {f'{r} {t}: {a} -> {b}': v for (r, t, a, b), v in sorted(c_vs_pa.items())}}
    for tag in READINGS:
        delta = sum((score.SCORE_OF[w] - score.SCORE_OF.get(o, 0.0)) * v for (o, w), v in tr[tag].items())
        res[tag] = {'changed': sum(tr[tag].values()),
                    'to_correct': sum(v for (o, w), v in tr[tag].items() if w == 'correct'),
                    'from_correct': sum(v for (o, w), v in tr[tag].items() if o == 'correct'),
                    'score': (sums['current'] + delta) / n,
                    'transitions': {f'{o} -> {w}': v for (o, w), v in sorted(tr[tag].items())},
                    'by_template': dict(by_t[tag]), 'to_correct_by_template': dict(to_correct_t[tag])}
    return key, res, changes


def full_run(keys: list[str], write: bool) -> dict:
    out, changes = {}, []
    with cf.ProcessPoolExecutor(max_workers=min(8, len(keys))) as pool:
        for key, res, ch in pool.map(_model, keys):
            out[key] = res
            changes.extend(ch)
    if write:
        with open(CHANGES, 'w', encoding='utf-8', newline='\n') as fh:
            for c in changes:
                fh.write(json.dumps(c, ensure_ascii=False) + '\n')
    return out


# ------------------------------------------------------------------ report

def report(g: collections.Counter, p: dict, f: dict, keys: list[str], mv: dict) -> list[str]:
    L = ['# Three readings of the answer check, measured before any change', '',
         'Generated by `answer_form_audit.py`; the readings and the comparison are defined in its docstring. '
         'Counts only; `answer.py` is unchanged. P reads pi-fractions as one value; A, on top of P, accepts the '
         'target of a scalar gold line in another unit the line also states; S only counts symbolic answers written '
         'as their numeric value.', '',
         '## Gold, all 2,250 items, under P, P+A+C and P+A+N', '',
         '| check | result |', '|---|---:|',
         f"| gold items scored correct under every reading | {2250 - sum(g.values())} of 2250 |"]
    if g:
        L += ['', 'Templates whose gold fails: ' + ', '.join(f'`{t}` {n}' for t, n in sorted(g.items()))]
    L += ['', f"## The pilot's {p['traces']} traces against the experts", '',
          '| reading | non-partial agreement | three-way agreement | verdicts moved | to the experts | away from the experts |',
          '|---|---:|---:|---:|---:|---:|']
    for tag in ('current',) + READINGS:
        r = p[tag]
        L.append(f"| {tag} | {r['non_partial']:.3f} | {r['three_way']:.3f} | {r['moved']} | {r['to_experts']} | "
                 f"{r['from_experts']} |")
    for tag in READINGS:
        if p[tag]['transitions']:
            L += ['', f'Transitions under {tag}: ' + '; '.join(f'{k} {v}' for k, v in sorted(p[tag]['transitions'].items()))]
    L += ['', '## The full run', '',
          'Per model, over its 2,250 final traces (unusable rows unchanged): the answered incorrect verdicts now, '
          'the verdicts each reading changes, and the answer score now and under each reading. "Mismatch" is the '
          "number of rows whose recomputed current verdict differs from the store's label; it should be 0. C and N "
          'are each read on top of P+A: under C a number that is pi\'s coefficient is not a standalone value; under '
          'N a pi-fraction of two targets is one target, its value.', '',
          '| model | incorrect now | mismatch | P: changed | to correct | from correct | P+A: changed | to correct | '
          'from correct | P+A+C: changed | to correct | from correct | P+A+N: changed | to correct | from correct | '
          'answer score now | under P | under P+A | under P+A+C | under P+A+N |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for k in keys:
        r = f[k]
        L.append(f"| `{k}` | {r['incorrect_now']} | {r['mismatch']} | " + ' | '.join(
            f"{r[tag]['changed']} | {r[tag]['to_correct']} | {r[tag]['from_correct']}" for tag in READINGS)
            + f" | {r['score_current']:.3f} | " + ' | '.join(f"{r[tag]['score']:.3f}" for tag in READINGS) + ' |')
    L += ['', 'Transitions per model:', '']
    for k in keys:
        L.append(f"- `{k}`: " + ' | '.join(
            f"{tag}: " + ('; '.join(f'{t} {v}' for t, v in f[k][tag]['transitions'].items()) or 'none')
            for tag in READINGS))
    L += ['', 'What C and N each change against P+A, per model (reading, template: P+A verdict -> verdict, count):', '']
    for k in keys:
        L.append(f"- `{k}`: " + ('; '.join(f'{t} {v}' for t, v in f[k]['c_vs_pa'].items()) or 'nothing'))
    by_t = {tag: collections.Counter() for tag in READINGS}
    to_c = {tag: collections.Counter() for tag in READINGS}
    for k in keys:
        for tag in READINGS:
            by_t[tag].update(f[k][tag]['by_template'])
            to_c[tag].update(f[k][tag]['to_correct_by_template'])
    L += ['', '### By template, over the models above', '',
          '| template | ' + ' | '.join(f'{tag}: changed | to correct' for tag in READINGS) + ' |',
          '|---|' + '---:|---:|' * len(READINGS)]
    for t, n in sum((by_t[tag] for tag in READINGS), collections.Counter()).most_common():
        L.append(f"| `{t}` | " + ' | '.join(f"{by_t[tag][t]} | {to_c[tag][t]}" for tag in READINGS) + ' |')
    L += ['', '### Scalar and array templates whose gold answer line states more than one computed value', '',
          'The check verifies the last of them (`answer.targets`). "Renderings" counts the items on which every '
          'pair of values is one quantity in two units, which reading A accepts; on the other items the line '
          'states different quantities, which a scalar or array type cannot score part by part: a trace right on '
          'the last alone is correct and one right on another alone is incorrect.', '',
          '| template | answer type | items with 2+ values | most values on a line | items whose values are '
          'renderings of one quantity |', '|---|---|---:|---:|---:|']
    for t, r in sorted(mv.items()):
        L.append(f"| `{t}` | {r['type']} | {r['items']} | {r['max']} | {r['renderings']} |")
    L += ['', '## Reading S: `ber_estimation_mary`, symbolic gold written as its numeric value', '',
          'Incorrect and partial verdicts on the template whose gold line is `a * Q(b)`, and how many of them state '
          'a*Q(b) itself (within 1% or one unit of the last digit shown). Counted only; D-138 scores a symbolic answer '
          'by the numbers it states.', '',
          '| model | incorrect or partial | state the numeric value |', '|---|---:|---:|']
    for k in keys:
        b = f[k]['ber']
        L.append(f"| `{k}` | {b.get('verdicts', 0)} | {b.get('numeric', 0)} |")
    L.append(f"| all | {sum(f[k]['ber'].get('verdicts', 0) for k in keys)} | "
             f"{sum(f[k]['ber'].get('numeric', 0) for k in keys)} |")
    return L


# ------------------------------------------------------------------ the fix, between pinned commits

PRE_FIX = 'e116b4f'       # the last commit before D-169's change to answer.py
POST_FIX = None           # the commit that holds the change; set when it exists (see --fix)
FIX_CHANGES = score.SCORES / 'answer_form_fix_changes.jsonl'


def module_at(commit: str):
    """answer.py as it was at a commit, with milestones.py as it stands (D-169 does not touch it)."""
    src = subprocess.run(['git', 'show', f'{commit}:evaluator_pilot_17092026/evaluators/answer.py'],
                         cwd=REPO, capture_output=True, check=True).stdout.decode('utf-8')
    mod = types.ModuleType(f'answer_{commit[:7]}')
    mod.__file__ = str(EVAL / 'answer.py')
    exec(compile(src, f'answer_{commit[:7]}.py', 'exec'), mod.__dict__)
    return mod


def _model_fix(key: str):
    """One model's traces under answer.py at PRE_FIX and at POST_FIX."""
    OLD, NEW = module_at(PRE_FIX), module_at(POST_FIX)
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    store = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                           .read_text(encoding='utf-8').splitlines())}
    tr, by_t, to_c = collections.Counter(), collections.Counter(), collections.Counter()
    sums = collections.Counter()
    changes = []
    n = mismatch = 0
    for row in score.trace_rows('main', key):
        s = store[row['item_id']]
        n += 1
        cur = 'unusable' if s['unusable'] else s['answer']['label']
        sums['old'] += score.SCORE_OF.get(cur, 0.0)
        if s['unusable']:
            continue
        item = items[row['item_id']]
        vals = tuple(m['value'] for m in ms_all[row['item_id']])
        text = row.get('text') or ''
        old = OLD.verdict(text, item, None, vals)[0]
        new = NEW.verdict(text, item, None, vals)[0]
        if old != cur:
            mismatch += 1
        if new != old:
            tr[(old, new)] += 1
            by_t[item['template_id']] += 1
            to_c[item['template_id']] += new == 'correct'
            changes.append({'model_key': key, 'item_id': row['item_id'], 'old': old, 'new': new,
                            'segment': NEW.segment(text)})
    delta = sum((score.SCORE_OF[w] - score.SCORE_OF.get(o, 0.0)) * v for (o, w), v in tr.items())
    return key, {'traces': n, 'mismatch': mismatch, 'changed': sum(tr.values()),
                 'to_correct': sum(v for (o, w), v in tr.items() if w == 'correct'),
                 'from_correct': sum(v for (o, w), v in tr.items() if o == 'correct'),
                 'score_old': sums['old'] / n, 'score_new': (sums['old'] + delta) / n,
                 'transitions': {f'{o} -> {w}': v for (o, w), v in sorted(tr.items())},
                 'by_template': dict(by_t), 'to_correct_by_template': dict(to_c)}, changes


def fix_report(keys: list[str]) -> list[str]:
    if not POST_FIX:
        raise SystemExit('set POST_FIX to the commit that holds the change first')
    OLD, NEW = module_at(PRE_FIX), module_at(POST_FIX)
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    g = collections.Counter()
    for iid, it in items.items():
        vals = tuple(m['value'] for m in ms_all[iid])
        for tag, M in (('before', OLD), ('after', NEW)):
            if M.verdict(it['solution'], it, None, vals)[0] != 'correct':
                g[(tag, it['template_id'])] += 1
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, pitems, texts = V.pilot_inputs()
    labs = {'before': {}, 'after': {}}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in pitems or key not in texts:
            continue
        vals = tuple(m['value'] for m in t['milestones'])
        labs['before'][code] = OLD.verdict(texts[key], pitems[key[1]], None, vals)[0]
        labs['after'][code] = NEW.verdict(texts[key], pitems[key[1]], None, vals)[0]
    codes = sorted(labs['before'])
    exp = {c: truth[c]['final_answer'] for c in codes}
    npc = [c for c in codes if exp[c] != 'partial']
    moved = [c for c in codes if labs['before'][c] != labs['after'][c]]
    out, changes = {}, []
    with cf.ProcessPoolExecutor(max_workers=min(8, len(keys))) as pool:
        for key, res, ch in pool.map(_model_fix, keys):
            out[key] = res
            changes.extend(ch)
    with open(FIX_CHANGES, 'w', encoding='utf-8', newline='\n') as fh:
        for c in changes:
            fh.write(json.dumps(c, ensure_ascii=False) + '\n')
    L = ['# The answer check before and after D-169, between pinned commits', '',
         f'Generated by `answer_form_audit.py --fix`: `answer.py` at `{PRE_FIX}` against `{POST_FIX}`, both loaded '
         'from git, on the gold, the pilot\'s 300 expert-labelled traces and every full-run trace, with the same '
         'milestones. Counts only; the changed verdicts and their answer segments are in '
         '`scores/answer_form_fix_changes.jsonl`, local. The readings the change was designed from are measured in '
         '`ANSWER_FORM_AUDIT.md`.', '',
         '## Gold, all 2,250 items', '', '| | before | after |', '|---|---:|---:|',
         f"| gold items scored correct | {2250 - sum(v for (t, _), v in g.items() if t == 'before')} of 2250 | "
         f"{2250 - sum(v for (t, _), v in g.items() if t == 'after')} of 2250 |", '',
         f"## The pilot's {len(codes)} traces against the experts", '',
         '| | before | after |', '|---|---:|---:|']
    for name, f_ in (('non-partial agreement', lambda L_: sum((L_[c] == 'correct') == (exp[c] == 'correct') for c in npc) / len(npc)),
                     ('three-way agreement', lambda L_: sum(L_[c] == exp[c] for c in codes) / len(codes))):
        L.append(f"| {name} | {f_(labs['before']):.3f} | {f_(labs['after']):.3f} |")
    L += [f"| verdicts moved | | {len(moved)} ({sum(labs['after'][c] == exp[c] for c in moved)} to the experts, "
          f"{sum(labs['before'][c] == exp[c] for c in moved)} away) |", '',
          '## The full run', '',
          'Per model, over its 2,250 final traces (unusable rows unchanged). "Mismatch" counts rows whose verdict '
          "under the code at PRE_FIX differs from the store's label; it should be 0.", '',
          '| model | mismatch | verdicts changed | to correct | from correct | answer score before | after |',
          '|---|---:|---:|---:|---:|---:|---:|']
    for k in keys:
        r = out[k]
        L.append(f"| `{k}` | {r['mismatch']} | {r['changed']} | {r['to_correct']} | {r['from_correct']} | "
                 f"{r['score_old']:.3f} | {r['score_new']:.3f} |")
    L += ['', 'Transitions per model:', '']
    L += [f"- `{k}`: " + ('; '.join(f'{t} {v}' for t, v in out[k]['transitions'].items()) or 'none') for k in keys]
    by_t, to_c = collections.Counter(), collections.Counter()
    for k in keys:
        by_t.update(out[k]['by_template'])
        to_c.update(out[k]['to_correct_by_template'])
    L += ['', '### By template', '', '| template | verdicts changed | to correct |', '|---|---:|---:|']
    L += [f"| `{t}` | {n} | {to_c[t]} |" for t, n in by_t.most_common()]
    return L


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', nargs='*', help='model keys; with it nothing is written')
    ap.add_argument('--fix', action='store_true', help='answer.py at PRE_FIX against POST_FIX; writes ANSWER_FORM_FIX.md')
    a = ap.parse_args()
    from full_run_28092026.analyze import ROSTER
    keys = a.models or ROSTER
    write = not a.models
    if a.fix:
        L = fix_report(keys)
        if write:
            (HERE / 'ANSWER_FORM_FIX.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
        print('\n'.join(L))
        return 0
    P = variant_module()
    PC = variant_module(coef=True)
    PN = variant_module(pifrac=True)
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    g = gold(P, items, ms_all) + gold(PC, items, ms_all) + gold(PN, items, ms_all)
    mv = multi_value(P, items, ms_all)
    p = pilot(P, PC, PN)
    f = full_run(keys, write)
    L = report(g, p, f, keys, mv)
    if write:
        (HERE / 'ANSWER_FORM_AUDIT.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
