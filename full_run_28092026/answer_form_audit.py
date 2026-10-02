"""Three readings of the answer check that the full run's "incorrect" verdicts point at, measured
before any change, on the gold, the pilot's 300 expert-labelled traces and every full-run trace.

    python -m full_run_28092026.answer_form_audit                   # FREE: writes ANSWER_FORM_AUDIT.md beside this file
    python -m full_run_28092026.answer_form_audit --models kimi-k3  # some models only; writes nothing

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
            is correct.
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


def variant_module():
    """answer.py as it stands, loaded as a second module whose `values` normalises first (reading P)."""
    src = (EVAL / 'answer.py').read_text(encoding='utf-8')
    mod = types.ModuleType('answer_P')
    mod.__file__ = str(EVAL / 'answer.py')
    exec(compile(src, 'answer_P.py', 'exec'), mod.__dict__)
    orig = mod.values
    mod.values = lambda seg: orig(normalise(seg))
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


def verdict_pa(P, text: str, item: dict, vals: tuple) -> str:
    """Reading A on top of P: a scalar item is correct when the trace matches another computed value the
    gold's answer line states that is the target in another unit."""
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
        if P.match(v, have, P.REL, exact, gu) or P.match(abs(v), absolute, P.REL, exact, gu):
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


def pilot(P) -> dict:
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, items, texts = V.pilot_inputs()
    labs = {'current': {}, 'P': {}, 'P+A': {}}
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

def _model(key: str):
    P = variant_module()
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    store = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                           .read_text(encoding='utf-8').splitlines())}
    tr = {'P': collections.Counter(), 'P+A': collections.Counter()}
    by_t = {'P': collections.Counter(), 'P+A': collections.Counter()}
    to_correct_t = {'P': collections.Counter(), 'P+A': collections.Counter()}
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
        for tag, lab in (('P', p), ('P+A', pa)):
            if lab != cur:
                tr[tag][(cur, lab)] += 1
                by_t[tag][item['template_id']] += 1
                if lab == 'correct':
                    to_correct_t[tag][item['template_id']] += 1
        if pa != cur:
            changes.append({'model_key': key, 'item_id': row['item_id'], 'current': cur, 'P': p, 'P+A': pa,
                            'segment': answer.segment(text)})
        if item['template_id'] == BER_TEMPLATE and cur in ('incorrect', 'partial'):
            hit = ber_numeric(item, text)
            if hit is not None:
                ber['verdicts'] += 1
                ber['numeric'] += hit
    res = {'traces': n, 'mismatch': mismatch, 'incorrect_now': sums['incorrect_now'],
           'score_current': sums['current'] / n, 'ber': dict(ber)}
    for tag in ('P', 'P+A'):
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
         '## Gold, all 2,250 items, under P', '',
         '| check | result |', '|---|---:|',
         f"| gold items scored correct under P | {2250 - sum(g.values())} of 2250 |"]
    if g:
        L += ['', 'Templates whose gold fails under P: ' + ', '.join(f'`{t}` {n}' for t, n in sorted(g.items()))]
    L += ['', f"## The pilot's {p['traces']} traces against the experts", '',
          '| reading | non-partial agreement | three-way agreement | verdicts moved | to the experts | away from the experts |',
          '|---|---:|---:|---:|---:|---:|']
    for tag in ('current', 'P', 'P+A'):
        r = p[tag]
        L.append(f"| {tag} | {r['non_partial']:.3f} | {r['three_way']:.3f} | {r['moved']} | {r['to_experts']} | "
                 f"{r['from_experts']} |")
    for tag in ('P', 'P+A'):
        if p[tag]['transitions']:
            L += ['', f'Transitions under {tag}: ' + '; '.join(f'{k} {v}' for k, v in sorted(p[tag]['transitions'].items()))]
    L += ['', '## The full run', '',
          'Per model, over its 2,250 final traces (unusable rows unchanged): the answered incorrect verdicts now, '
          'the verdicts each reading changes, and the answer score now and under each reading. "Mismatch" is the '
          "number of rows whose recomputed current verdict differs from the store's label; it should be 0.", '',
          '| model | incorrect now | mismatch | P: changed | to correct | from correct | P+A: changed | to correct | '
          'from correct | answer score now | under P | under P+A |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for k in keys:
        r = f[k]
        L.append(f"| `{k}` | {r['incorrect_now']} | {r['mismatch']} | {r['P']['changed']} | {r['P']['to_correct']} | "
                 f"{r['P']['from_correct']} | {r['P+A']['changed']} | {r['P+A']['to_correct']} | "
                 f"{r['P+A']['from_correct']} | {r['score_current']:.3f} | {r['P']['score']:.3f} | "
                 f"{r['P+A']['score']:.3f} |")
    L += ['', 'Transitions per model:', '']
    for k in keys:
        L.append(f"- `{k}`: P: " + ('; '.join(f'{t} {v}' for t, v in f[k]['P']['transitions'].items()) or 'none')
                 + ' | P+A: ' + ('; '.join(f'{t} {v}' for t, v in f[k]['P+A']['transitions'].items()) or 'none'))
    by_t = {tag: collections.Counter() for tag in ('P', 'P+A')}
    to_c = {tag: collections.Counter() for tag in ('P', 'P+A')}
    for k in keys:
        for tag in ('P', 'P+A'):
            by_t[tag].update(f[k][tag]['by_template'])
            to_c[tag].update(f[k][tag]['to_correct_by_template'])
    L += ['', '### By template, over the models above', '',
          '| template | P: changed | to correct | P+A: changed | to correct |', '|---|---:|---:|---:|---:|']
    for t, n in by_t['P+A'].most_common():
        L.append(f"| `{t}` | {by_t['P'][t]} | {to_c['P'][t]} | {n} | {to_c['P+A'][t]} |")
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', nargs='*', help='model keys; with it nothing is written')
    a = ap.parse_args()
    from full_run_28092026.analyze import ROSTER
    keys = a.models or ROSTER
    write = not a.models
    P = variant_module()
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    g = gold(P, items, ms_all)
    mv = multi_value(P, items, ms_all)
    p = pilot(P)
    f = full_run(keys, write)
    L = report(g, p, f, keys, mv)
    if write:
        (HERE / 'ANSWER_FORM_AUDIT.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
