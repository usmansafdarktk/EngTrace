"""Does the pilot's trace population span what a small-model roster would produce?

Every candidate evaluator reads structure, not prose: extract_steps() needs step
markers, the final-answer check needs a parseable number, E3 will need (value, unit)
pairs per step. If all five pilot models write the same tidy structure, an
evaluator validated here says nothing about a 27B model that does not.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import glob
import json
import os
import re
import sys

ROOT = _REPO
sys.path.insert(0, os.path.join(ROOT, 'evaluation'))
from engineering_parser import extract_steps            # noqa: E402

P = os.path.join(ROOT, 'evaluator_pilot_17092026')
items = {json.loads(l)['item_id']: json.loads(l)
         for l in open(os.path.join(P, 'slice', 'manifest.jsonl'), encoding='utf-8')}

print('%-17s %6s %7s %8s %9s %9s %10s %9s' %
      ('model', 'traces', 'steps', 'no mark', 'has Answer', 'final=None', 'med chars', 'num/step'))
gold_steps = []
for f in sorted(glob.glob(os.path.join(P, 'traces', '*.jsonl'))):
    rows = [json.loads(l) for l in open(f, encoding='utf-8')]
    ok = {r['item_id']: r for r in rows if r['ok'] and
          r.get('finish_reason') in ('stop', 'end_turn', 'eos', None)}
    key = os.path.basename(f)[:-6]
    n_steps, unmarked, has_ans, no_final, lens, numeric = [], 0, 0, 0, [], []
    for iid, r in ok.items():
        t = r['text']
        steps, vals, final = extract_steps(t)
        marked = bool(re.search(r'(?:^|\n)\s*(?:\*\*|#*)\s*(?:Step\s*\d+|\d+\.)', t, re.I))
        n_steps.append(len(steps))
        unmarked += (not marked)
        has_ans += ('**Answer:**' in t)
        no_final += (final is None)
        lens.append(len(t))
        numeric.append(sum(v is not None for v in vals) / max(len(steps), 1))
    med = lambda a: sorted(a)[len(a) // 2]
    print('%-17s %6d %7.1f %8s %9s %9s %10d %9.2f'
          % (key, len(ok), sum(n_steps) / len(n_steps),
             '%d (%.0f%%)' % (unmarked, 100 * unmarked / len(ok)),
             '%d (%.0f%%)' % (has_ans, 100 * has_ans / len(ok)),
             '%d (%.0f%%)' % (no_final, 100 * no_final / len(ok)),
             med(lens), sum(numeric) / len(numeric)))

for iid, it in items.items():
    gold_steps.append(len(extract_steps(it['solution'])[0]))
print('\ngold solutions: %.1f steps on average' % (sum(gold_steps) / len(gold_steps)))
