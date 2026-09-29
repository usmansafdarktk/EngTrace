"""D-147: the answer check's last-digit windows, measured before the boundary rule changed.

    python -m full_run_28092026.boundary_audit     # FREE: writes BOUNDARY_AUDIT.md beside this file

THE WEAKNESS. `answer.match` accepts a stated number within one unit of its own last digit, or of the
gold's. A value exactly one unit off sits on the window's edge, and in binary floating point the edge
is not exact: 0.063 - 0.062 is 0.0010000000000000009, so whether such a value passed depended on how
the difference happened to round (the second review, 2026-09-29). The windows are documented as
inclusive ("within one unit"), so the rule adopted is the inclusive one with a relative slack of 1e-9
(answer.SLACK): the verdict follows the rule, not the rounding.

THREE READINGS, each measured against the check as committed before this change ("before", loaded
from git at the commit given):
  inclusive    the rule adopted: the same windows, the boundary inclusive whatever the rounding
  half-unit    the stricter reading, a correct rounding at the precision shown (half a unit either
               way): a sensitivity, not the score, because the experts validated the one-unit rule
  whole trace  a numeric part the Answer line leaves out is credited when the trace states it
               anywhere, for a trace whose Answer line already matches at least one part: a
               quantity the question asks for, computed in the body and left off the Answer line,
               then counts. A sensitivity, not the score, because the prompt asks for the final
               result on the Answer line and the experts validated the check on that line; the
               pilot holds only two such traces, which they split one each way. A first version of
               this reading re-read every part from the whole trace, which credited wrong Answer
               lines whose working held the right number: 1,089 incorrect verdicts became correct
               and 27 pilot verdicts moved away from the experts, so it was narrowed to omissions
Measured on the gold (every one must stay correct), on the pilot's 300 traces against the experts,
and on every full-run trace. Counts only.
"""
from __future__ import annotations

import collections
import concurrent.futures as cf
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import parser_fix, score  # noqa: E402

BEFORE = '5b87a15'          # the last commit of answer.py before D-147
READINGS = ('before', 'inclusive', 'half-unit', 'whole trace')


def _before():
    return parser_fix.modules_at(BEFORE, 'd147_before')[0]


def label(name, A, text, item, vals, check_readable=True):
    if check_readable and not score.readable(text):
        return 'unusable'
    if name == 'before':
        return A.verdict(text, item, None, vals)[0]
    if name == 'inclusive':
        return answer.verdict(text, item, None, vals)[0]
    if name == 'half-unit':
        return answer.verdict(text, item, None, vals, unit=0.5)[0]
    return answer.verdict(text, item, None, vals, whole=True)[0]


def _gold(name):
    A = _before() if name == 'before' else None
    items = score.pool_items()
    ms = score.milestone_sets(items)
    return name, sum(label(name, A, it['solution'], it, tuple(m['value'] for m in ms[i]), False) == 'correct'
                     for i, it in items.items())


def _pilot(name):
    A = _before() if name == 'before' else None
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, items, texts = V.pilot_inputs()
    lab, exp = {}, {}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] in items and key in texts:
            lab[code] = label(name, A, texts[key], items[key[1]], tuple(m['value'] for m in t['milestones']), False)
            exp[code] = t['final_answer']
    np_codes = [c for c in lab if exp[c] != 'partial']
    return name, {'non_partial': sum((lab[c] == 'correct') == (exp[c] == 'correct') for c in np_codes) / len(np_codes),
                  'three_way': sum(lab[c] == exp[c] for c in lab) / len(lab), 'labels': lab, 'experts': exp}


def _model(key):
    A = _before()
    items = score.pool_items()
    ms = score.milestone_sets(items)
    rows = score.trace_rows('main', key)
    out = {}
    for name in READINGS:
        out[name] = {r['item_id']: (label(name, A, r.get('text') or '', items[r['item_id']],
                                          tuple(m['value'] for m in ms[r['item_id']]))
                                    if r.get('status') == 'answered' else 'unusable') for r in rows}
    return key, out


def main() -> int:
    from full_run_28092026.run_traces import config as roster
    keys = [m['key'] for m in roster()['models'] if (score.TRACES / f"{m['key']}.jsonl").exists()]
    with cf.ProcessPoolExecutor(max_workers=8) as pool:
        gold = dict(pool.map(_gold, READINGS))
        pilot = dict(pool.map(_pilot, READINGS))
        runs = dict(pool.map(_model, keys))
    base, exp = pilot['before']['labels'], pilot['before']['experts']
    val = {'correct': 1.0, 'partial': 0.5}
    L = ['# The answer check\'s last-digit windows: boundary, half-unit and whole-trace readings (D-147)', '',
         'Generated by `boundary_audit.py`; the readings are defined in its docstring. Counts only; "before" is '
         f'`answer.py` at `{BEFORE}`.', '', '## The gold and the pilot', '',
         '| reading | gold scored correct | pilot non-partial agreement | pilot three-way agreement | '
         'pilot verdicts changed | to the experts\' verdict | away from it |', '|---|---:|---:|---:|---:|---:|---:|']
    for name in READINGS:
        lab = pilot[name]['labels']
        moved = [c for c in lab if lab[c] != base[c]]
        L.append(f"| {name} | {gold[name]} of 2250 | {pilot[name]['non_partial']:.3f} | {pilot[name]['three_way']:.3f} | "
                 f"{len(moved)} | {sum(lab[c] == exp[c] for c in moved)} | {sum(base[c] == exp[c] for c in moved)} |")
    L += ['', '## The full run', '',
          'Per model, the traces whose verdict differs from the check before, and the answer score under each reading '
          '(correct 1, partial 0.5, otherwise 0).', '',
          '| model | before | inclusive: changed | score | half-unit: changed | score | whole trace: changed | score |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for k, v in runs.items():
        cells = [f"{sum(val.get(x, 0) for x in v['before'].values()) / len(v['before']):.3f}"]
        for name in READINGS[1:]:
            ch = sum(v[name][i] != v['before'][i] for i in v['before'])
            cells += [str(ch), f"{sum(val.get(x, 0) for x in v[name].values()) / len(v[name]):.3f}"]
        L.append(f'| `{k}` | ' + ' | '.join(cells) + ' |')
    items = score.pool_items()
    for name in READINGS[1:]:
        tr, by_type, by_template = collections.Counter(), collections.Counter(), collections.Counter()
        for v in runs.values():
            for i in v['before']:
                if v[name][i] != v['before'][i]:
                    tr[f"{v['before'][i]} -> {v[name][i]}"] += 1
                    by_type[items[i]['answer_type']] += 1
                    by_template[items[i]['template_id'].removeprefix('template_')] += 1
        L += ['', f'**{name}** over all {len(runs)} models: {sum(tr.values())} verdicts change'
              + (': ' + '; '.join(f'{t} {n}' for t, n in tr.most_common()) if tr else '') + '.',
              'By answer type: ' + (', '.join(f'{t} {n}' for t, n in by_type.most_common()) or 'none') + '. '
              f'By template, {len(by_template)}: '
              + (', '.join(f'`{t}` {n}' for t, n in by_template.most_common(12)) or 'none')
              + (' ...' if len(by_template) > 12 else '') + '.']
    (HERE / 'BOUNDARY_AUDIT.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
