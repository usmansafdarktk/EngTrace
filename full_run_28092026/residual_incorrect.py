"""How far the remaining "incorrect" answers sit from their gold target, after D-169.

    python -m full_run_28092026.residual_incorrect                    # FREE: writes RESIDUAL_INCORRECT.md beside this file
    python -m full_run_28092026.residual_incorrect --sample 6 --templates template_work_isothermal_virial ...
                                                                      # prints answer segments to stdout for reading; writes nothing

WHY. D-168 read twelve "incorrect" verdicts of the top five models and found six to be right answers in
another form; D-169 fixed those forms and moved 238 verdicts. Nobody has yet measured what the verdicts
that remain look like. This script does not judge them; it measures, for every answered, usable trace whose
label is `incorrect`, the closest approach of any number on the trace's Answer segment to any numeric
target, under the check's own unit scales (`answer.values`, `answer.segment`, `answer.SCALES`), and buckets
that relative distance. A trace whose closest approach is under 1% is a near miss, which on an iterative
template may be a convergence difference rather than a wrong method; one over 20% is a wrong answer by any
reading; a segment with no number is a trace that stated no numeric answer where one was asked. Items whose
targets are words (classification) are counted apart, as are symbolic items, whose targets are the numbers
inside the gold's expression (D-138) and whose incorrect verdicts D-168's reading S characterised.

Counts only go into the report; the segments printed by `--sample` stay on the terminal.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import score  # noqa: E402

MODELS = ['deepseek-v4.1-flash', 'kimi-k3', 'claude-sonnet-5', 'glm-5.3-flash', 'muse-glimmer-30b',
          'glm-5.3', 'qwen3-235b-a22b-2507', 'gemini-3.1-flash-lite', 'gemma-4-26b-a4b', 'gpt-5.4-mini',
          'gpt-oss-20b']
TOP5 = MODELS[:5]
BUCKETS = ['<=0.2%', '0.2-1%', '1-5%', '5-20%', '>20%', 'no number', 'word target', 'symbolic']
OUT = HERE / 'RESIDUAL_INCORRECT.md'


def closest(numbers: list[float], have: list[tuple[float, float]]) -> float | None:
    best = None
    for g in numbers:
        for v, _u in have:
            for sc in answer.SCALES:
                d = abs(v * sc - g) / max(abs(g), 1e-12)
                best = d if best is None or d < best else best
    return best


def bucket(d: float | None) -> str:
    if d is None:
        return 'no number'
    if d <= 0.002:
        return '<=0.2%'
    if d <= 0.01:
        return '0.2-1%'
    if d <= 0.05:
        return '1-5%'
    if d <= 0.20:
        return '5-20%'
    return '>20%'


def rows(key: str):
    store = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                           .read_text(encoding='utf-8').splitlines())}
    for row in score.trace_rows('main', key):
        s = store[row['item_id']]
        if s['status'] != 'answered' or s['unusable'] or s['answer']['label'] != 'incorrect':
            continue
        yield s, row.get('text') or ''


def measure(keys: list[str]):
    by_model = {k: collections.Counter() for k in keys}
    by_template = collections.defaultdict(collections.Counter)
    by_template_top5 = collections.defaultdict(collections.Counter)
    within = collections.defaultdict(collections.Counter)   # the <=0.2% rows, by why the check still said no
    for key in keys:
        for s, text in rows(key):
            t = s['answer']['targets']
            if s['answer_type'] == 'symbolic':
                b = 'symbolic'
            elif not t['numbers']:
                b = 'word target'
            else:
                b = bucket(closest(t['numbers'], answer.values(answer.segment(text))))
            by_model[key][b] += 1
            by_template[s['template_id']][b] += 1
            if key in TOP5:
                by_template_top5[s['template_id']][b] += 1
            if b == '<=0.2%':
                why = 'exact digits required' if t['exact'] else (
                    'a word or labelled part' if t['words'] or t['labeled'] else 'several numeric targets' if len(t['numbers']) > 1 else 'other')
                within[s['template_id']][why] += 1
    return by_model, by_template, by_template_top5, within


VIRIAL = 'template_work_isothermal_virial'
FLAME = 'template_adiabatic_flame_temperature'
import re  # noqa: E402

P_RANGE = re.compile(r'from\s+([\d.]+)\s*bar\s+to\s+([\d.]+)\s*bar')
B_LINE = re.compile(r'\bB\s*=\s*\(R.*?=\s*(-?[\d.]+)\s*L/mol')
W_IDEAL = re.compile(r'W_ideal\s*=.*?=\s*(-?[\d.]+)\s*J/mol')


def virial_readings(keys: list[str], items: dict) -> dict:
    """On the virial template's incorrect traces: does the stated answer equal the FLOW work under the
    pressure-explicit two-term form, W = RT ln(P2/P1) + B (P2 - P1), or the ideal-gas work alone? Both are
    computed from the gold's own stated B and W_ideal, and matched with the check's own rule (0.2% or one
    unit of the last digit), so the count says how many "wrong" answers are a different, stated reading of
    an under-specified question rather than a slip."""
    out = collections.Counter()
    for key in keys:
        for s, text in rows(key):
            if s['template_id'] != VIRIAL:
                continue
            it = items[s['item_id']]
            p = P_RANGE.search(it['question'])
            b = B_LINE.search(it['solution'])
            w = W_IDEAL.search(it['solution'])
            out['incorrect'] += 1
            if not (p and b and w):
                out['not parsed'] += 1
                continue
            p1, p2 = float(p.group(1)), float(p.group(2))
            w_flow = float(w.group(1)) + float(b.group(1)) * (p2 - p1) * 100.0
            have = answer.values(answer.segment(text))
            if answer.match(w_flow, have):
                out['flow work, pressure-explicit form'] += 1
            elif answer.match(float(w.group(1)), have):
                out['ideal-gas work'] += 1
            else:
                out['neither'] += 1
    return out


def cliff_without(keys: list[str]) -> list[str]:
    """Each model's Easy and Advanced template means with and without the two templates, descriptive, no interval."""
    out = ['### The Easy-to-Advanced gap without the two chemical templates (descriptive, no interval)', '',
           'Both are Advanced templates; the gap is the mean of the Easy template means minus the mean of the Advanced '
           'ones, as Q2 defines it, with and without `work_isothermal_virial` and `adiabatic_flame_temperature` (32 '
           'Advanced templates remain). Points lost on the two: 30 minus the sum of the item scores on their 30 items.', '',
           '| model | Easy | Advanced | gap | Advanced without the two | gap without | points lost on the two, of 30 |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for key in keys:
        tm = collections.defaultdict(list)
        lv = {}
        lost = 0.0
        for line in (score.SCORES / 'main' / f'{key}.jsonl').read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            tm[r['template_id']].append(r['score'])
            lv[r['template_id']] = r['level']
            if r['template_id'] in (VIRIAL, FLAME):
                lost += 1 - r['score']
        def mean(ts):
            ms = [sum(v) / len(v) for t, v in tm.items() if t in ts]
            return sum(ms) / len(ms)
        easy = [t for t in tm if lv[t] == 'Easy']
        adv = [t for t in tm if lv[t] == 'Advanced']
        adv2 = [t for t in adv if t not in (VIRIAL, FLAME)]
        out.append(f'| `{key}` | {mean(easy):.3f} | {mean(adv):.3f} | {mean(easy) - mean(adv):+.3f} | {mean(adv2):.3f} | '
                   f'{mean(easy) - mean(adv2):+.3f} | {lost:.1f} |')
    out.append('')
    return out


def exact_templates(keys: list[str]) -> tuple[int, int]:
    """How many templates carry an exact-digits target, and how many of the roster's incorrect verdicts fall on them."""
    ts = set()
    inc = 0
    for key in keys:
        for line in (score.SCORES / 'main' / f'{key}.jsonl').read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            if r['status'] == 'answered' and not r['unusable'] and r['answer']['targets'].get('exact'):
                ts.add(r['template_id'])
                if r['answer']['label'] == 'incorrect':
                    inc += 1
    return len(ts), inc


def table(title: str, counters: dict, order: list[str]) -> list[str]:
    out = [f'### {title}', '', '| | n | ' + ' | '.join(BUCKETS) + ' |', '|---|---:|' + '---:|' * len(BUCKETS)]
    tot = collections.Counter()
    for k in order:
        c = counters[k]
        tot.update(c)
        out.append(f'| `{k}` | {sum(c.values())} | ' + ' | '.join(str(c[b]) for b in BUCKETS) + ' |')
    out.append(f'| all | {sum(tot.values())} | ' + ' | '.join(str(tot[b]) for b in BUCKETS) + ' |')
    out.append('')
    return out


def report(keys: list[str]) -> str:
    by_model, by_template, by_template_top5, within = measure(keys)
    lines = ['# The remaining "incorrect" verdicts: closest approach to the gold target', '',
             'Generated by `residual_incorrect.py`; the measure is defined in its docstring. Over every answered, usable '
             'trace the store labels `incorrect` (the store as re-scored under D-169), the smallest relative distance '
             'between any number on the Answer segment and any numeric target, under the check\'s unit scales. '
             'Counts only. This measures how far the stated answer is from the gold; it does not say whether a near miss '
             'is the model\'s error or the check\'s tolerance.', '']
    lines += table('By model', by_model, keys)
    top = sorted(by_template_top5, key=lambda t: -sum(by_template_top5[t].values()))[:15]
    lines += table('Top five models, by template (the fifteen templates with most)', by_template_top5, top)
    top_all = sorted(by_template, key=lambda t: -sum(by_template[t].values()))[:15]
    lines += table('All eleven models, by template (the fifteen templates with most)', by_template, top_all)
    lines += ['### Within 0.2% and still incorrect: why, by template (all eleven models)', '',
              'A number within the relative tolerance fails when the item requires exact digits (the question prescribes '
              'the rounding; `exact`), when a word or labelled part of the answer did not match, or when the item has '
              'several numeric targets and only some were within reach.', '',
              '| template | n | exact digits required | a word or labelled part | several numeric targets | other |',
              '|---|---:|---:|---:|---:|---:|']
    for t in sorted(within, key=lambda t: -sum(within[t].values())):
        c = within[t]
        lines.append(f"| `{t}` | {sum(c.values())} | {c['exact digits required']} | {c['a word or labelled part']} | "
                     f"{c['several numeric targets']} | {c['other']} |")
    lines.append('')
    n_t, n_inc = exact_templates(keys)
    lines += [f'Templates with an exact-digits target (the question prescribes the rounding): {n_t} of 150; incorrect '
              f'verdicts on them across the roster: {n_inc}.', '']
    items = score.pool_items()
    v = virial_readings(keys, items)
    lines += ['### `work_isothermal_virial`: what the incorrect answers state', '',
              'The question asks for the work to compress one mole isothermally and reversibly "based on the virial '
              'equation truncated to two terms", and names neither the closed-system form (W = -∫P dV with Z = 1 + B/V, '
              'which the gold uses) nor the flow form (W = ∫V dP with Z = 1 + BP/RT, which equals RT ln(P2/P1) + B(P2 - P1)). '
              'For each incorrect trace the flow-work value is computed from the gold\'s own stated B and ideal-gas work and '
              'matched with the check\'s rule (`answer.match`, 0.2% or one unit of the last digit).', '',
              '| incorrect traces | state the flow work under the pressure-explicit form | state the ideal-gas work | neither | not parsed |',
              '|---:|---:|---:|---:|---:|',
              f"| {v['incorrect']} | {v['flow work, pressure-explicit form']} | {v['ideal-gas work']} | {v['neither']} | {v['not parsed']} |", '']
    lines += cliff_without(keys)
    return '\n'.join(lines)


def sample(keys: list[str], templates: list[str], n: int) -> None:
    import re
    for key in keys:
        shown = collections.Counter()
        for s, text in rows(key):
            t = s['template_id']
            if templates and t not in templates:
                continue
            if shown[t] >= n:
                continue
            shown[t] += 1
            seg = re.sub(r'\s+', ' ', answer.segment(text))[:260]
            tg = s['answer']['targets']
            d = closest(tg['numbers'], answer.values(answer.segment(text))) if tg['numbers'] else None
            print(f"{key} | {s['item_id']} | {s['answer_type']} | targets {tg['numbers']} {tg['words']} | "
                  f"closest {d if d is None else f'{d:.4f}'} | {s['finish_reason']}\n    {seg}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', nargs='*', default=MODELS)
    ap.add_argument('--sample', type=int, default=0, help='print N answer segments per template per model; writes nothing')
    ap.add_argument('--templates', nargs='*', default=[])
    a = ap.parse_args()
    if a.sample:
        sample(a.models, a.templates, a.sample)
        return 0
    text = report(a.models)
    OUT.write_text(text + '\n', encoding='utf-8')
    print(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
