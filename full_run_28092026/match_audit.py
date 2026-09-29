"""D-138: stray digits and degenerate precision in the answer check, measured before any change.

    python -m full_run_28092026.match_audit     # FREE: writes MATCH_AUDIT.md beside this file

THE WEAKNESS. `answer.match` accepts a stated number within one unit of its own last digit, after
scaling by a unit factor up to 1e12. A number whose one-unit window reaches zero - a bare `0` or
`1`, `0.1` - therefore matches every gold at some scale: `1` at 1e6 is 1,000,000 +- 1,000,000.
And `answer.values` reads subscripts as numbers, so `p_1`, `x_{1}` put such a `1` in the answer.

TWO CANDIDATE RULES, each measured alone and together, against the D-137 check:
  S   a digit written as a subscript (`p_1`, `x_{1}`) is not a value, as E3's reader already has it
  R   a stated number vouches for a gold at its own precision only when its one-unit window excludes
      zero (|v| > ulp); the relative tolerance and the gold's own precision are unchanged
Measured on the gold (every one must stay correct), on the pilot's 300 traces against the experts,
and on every full-run trace. The variants are applied in memory, each defined here in full, so the
audit reads the same after answer.py adopted S+R (D-138) as before.
"""
from __future__ import annotations

import collections
import concurrent.futures as cf
import re
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

NUM = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')                  # answer.NUM at D-137
NUM_S = re.compile(r'(?<![\d.])(?<!_)(?<!_\{)[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')


def match_d137(gold, have, rel=answer.REL, exact=False, gold_ulp=0.0, unit=1.0):
    """answer.match at D-137 (`unit` is accepted and ignored: answer.verdict passes it since D-147)."""
    for v, u in have:
        for sc in answer.SCALES:
            t, tu = v * sc, u * sc
            if exact:
                if abs(t - gold) <= 1e-9 * max(1.0, abs(gold)):
                    return True
            elif abs(t - gold) <= max(rel * abs(gold), tu, gold_ulp):
                return True
    return False


MATCH = match_d137


def match_r(gold, have, rel=answer.REL, exact=False, gold_ulp=0.0, unit=1.0):
    for v, u in have:
        for sc in answer.SCALES:
            t, tu = v * sc, u * sc
            if exact:
                if abs(t - gold) <= 1e-9 * max(1.0, abs(gold)):
                    return True
                continue
            tol = max(rel * abs(gold), gold_ulp)
            if abs(v) > u:
                tol = max(tol, tu)
            if abs(t - gold) <= tol:
                return True
    return False


VARIANTS = {'D-137': (NUM, MATCH), 'S': (NUM_S, MATCH), 'R': (NUM, match_r), 'S+R': (NUM_S, match_r)}


def use(name):
    from full_run_28092026.word_audit import word_hit
    answer.NUM, answer.match = VARIANTS[name]
    answer._word_hit = lambda low, w: word_hit(low, w, False, False)     # the word rule before D-139


def label(text, item, vals):
    if not score.readable(text):
        return 'unusable'
    return answer.verdict(text, item, None, vals)[0]


def _gold(name):
    use(name)
    items = score.pool_items()
    ms = score.milestone_sets(items)
    return name, sum(answer.verdict(it['solution'], it, None, tuple(m['value'] for m in ms[i]))[0] == 'correct'
                     for i, it in items.items())


def _pilot(name):
    use(name)
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, items, texts = V.pilot_inputs()
    lab, exp = {}, {}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] in items and key in texts:
            lab[code] = answer.verdict(texts[key], items[key[1]], None, tuple(m['value'] for m in t['milestones']))[0]
            exp[code] = t['final_answer']
    np_codes = [c for c in lab if exp[c] != 'partial']
    return name, {'non_partial': sum((lab[c] == 'correct') == (exp[c] == 'correct') for c in np_codes) / len(np_codes),
                  'three_way': sum(lab[c] == exp[c] for c in lab) / len(lab), 'labels': lab, 'experts': exp}


def _model(key):
    items = score.pool_items()
    ms = score.milestone_sets(items)
    rows = score.trace_rows('main', key)
    out = {}
    for name in VARIANTS:
        use(name)
        out[name] = {r['item_id']: (label(r.get('text') or '', items[r['item_id']],
                                          tuple(m['value'] for m in ms[r['item_id']]))
                                    if r.get('status') == 'answered' else 'unusable') for r in rows}
    return key, out


def main() -> int:
    from full_run_28092026.run_traces import config as roster
    keys = [m['key'] for m in roster()['models'] if (score.TRACES / f"{m['key']}.jsonl").exists()]
    with cf.ProcessPoolExecutor(max_workers=8) as pool:
        gold = dict(pool.map(_gold, VARIANTS))
        pilot = dict(pool.map(_pilot, VARIANTS))
        runs = dict(pool.map(_model, keys))
    base = pilot['D-137']['labels']
    exp = pilot['D-137']['experts']
    L = ['# Stray digits and degenerate precision in the answer check (D-138)', '',
         'Generated by `match_audit.py`; the weakness and the candidate rules are defined in its docstring. '
         'Counts only.', '', '## The gold and the pilot', '',
         '| check | gold scored correct | pilot non-partial agreement | pilot three-way agreement | '
         'pilot verdicts changed | to the experts\' verdict | away from it |', '|---|---:|---:|---:|---:|---:|---:|']
    for name in VARIANTS:
        lab = pilot[name]['labels']
        moved = [c for c in lab if lab[c] != base[c]]
        L.append(f"| {name} | {gold[name]} of 2250 | {pilot[name]['non_partial']:.3f} | {pilot[name]['three_way']:.3f} | "
                 f"{len(moved)} | {sum(lab[c] == exp[c] for c in moved)} | {sum(base[c] == exp[c] for c in moved)} |")
    L += ['', '## The full run', '',
          'Per model, the traces whose verdict differs from the D-137 check, and the answer score under each rule '
          '(correct 1, partial 0.5, otherwise 0).', '',
          '| model | D-137 score | S: changed | S score | R: changed | R score | S+R: changed | S+R score |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    val = {'correct': 1.0, 'partial': 0.5}
    for k, v in runs.items():
        cells = [f"{sum(val.get(x, 0) for x in v['D-137'].values()) / len(v['D-137']):.3f}"]
        for name in ('S', 'R', 'S+R'):
            ch = sum(v[name][i] != v['D-137'][i] for i in v['D-137'])
            cells += [str(ch), f"{sum(val.get(x, 0) for x in v[name].values()) / len(v[name]):.3f}"]
        L.append(f'| `{k}` | ' + ' | '.join(cells) + ' |')
    tr, by_type, by_template = collections.Counter(), collections.Counter(), collections.Counter()
    items = score.pool_items()
    for v in runs.values():
        for i in v['D-137']:
            if v['S+R'][i] != v['D-137'][i]:
                tr[f"{v['D-137'][i]} -> {v['S+R'][i]}"] += 1
                by_type[items[i]['answer_type']] += 1
                by_template[items[i]['template_id'].removeprefix('template_')] += 1
    L += ['', f'S+R over all {len(runs)} models: {sum(tr.values())} verdicts change: '
          + '; '.join(f'{t} {n}' for t, n in tr.most_common()) + '.', '',
          'By answer type: ' + ', '.join(f'{t} {n}' for t, n in by_type.most_common()) + '.', '',
          f'By template, {len(by_template)} templates: '
          + ', '.join(f'`{t}` {n}' for t, n in by_template.most_common()) + '.']
    (HERE / 'MATCH_AUDIT.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
