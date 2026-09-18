"""Decompose E3's 0.258 null baseline.

Two very different things could produce it:
  (a) SHARED VALUES - a sibling item's milestone equals this item's own milestone
      or given (same substance, same constant). Then "reaching" it is not
      coincidence at all; the null is simply contaminated.
  (b) COINCIDENCE - a trace number lands within 2% of an unrelated value, helped
      by 13 unit factors. That is E3 matching noise, and it inflates real scores.
Only (b) is a defect. This measures each, and how much unit scaling contributes.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import glob
import json
import os
import statistics as st
import sys
from collections import defaultdict

ROOT = _REPO
P = os.path.join(ROOT, 'evaluator_pilot_17092026')
sys.path[:0] = [os.path.join(P, 'evaluators'), ROOT]
import milestones as ms            # noqa: E402

TOL = 0.02
items = {json.loads(l)['item_id']: json.loads(l) for l in open(os.path.join(P, 'slice', 'manifest.jsonl'), encoding='utf-8')}
M = ms.build_all(os.path.join(P, 'slice', 'manifest.jsonl'))
traces = {}
for f in glob.glob(os.path.join(P, 'traces', '*.jsonl')):
    k = os.path.basename(f)[:-6]
    if k in ('qwen3.8-27b', 'gemma-4-31b'):
        continue
    for l in open(f, encoding='utf-8'):
        r = json.loads(l)
        if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
            traces[(k, r['item_id'])] = r

by_t = defaultdict(list)
for iid, it in items.items():
    by_t[it['template_id']].append(iid)

res = defaultdict(list)
per_t = defaultdict(list)
for (mk, iid), tr in traces.items():
    nums = ms.numbers(tr['text'])
    own = [m['value'] for m in M[iid]['milestones']] + ms.numbers(items[iid]['question'])
    for other in by_t[items[iid]['template_id']]:
        if other == iid:
            continue
        sib = M[other]['milestones']
        if not sib:
            continue
        # (a) drop sibling milestones that ARE one of this item's own values/givens
        distinct = [m for m in sib if ms.scaled_match(m['value'], own, ms.DISPLAY_TOL) is None]
        shared = len(sib) - len(distinct)
        res['shared_fraction'].append(shared / len(sib))
        if not distinct:
            continue
        any_scale = sum(ms.scaled_match(m['value'], nums, TOL) is not None for m in distinct) / len(distinct)
        unit_only = sum(any(ms.close(n, m['value'], TOL) for n in nums) for m in distinct) / len(distinct)
        res['coincidence_any_scale'].append(any_scale)
        res['coincidence_scale_1_only'].append(unit_only)
        per_t[items[iid]['template_id']].append(any_scale)

print('sibling milestones that are genuinely SHARED with this item: %.1f%%'
      % (100 * st.mean(res['shared_fraction'])))
print()
print('pure coincidence rate, shared values removed:')
print('  with all 13 unit factors  %.3f' % st.mean(res['coincidence_any_scale']))
print('  exact units only (x1)     %.3f' % st.mean(res['coincidence_scale_1_only']))
print()
print('coincidence by template (all factors), worst first:')
for t, v in sorted(per_t.items(), key=lambda kv: -st.mean(kv[1])):
    print('  %-40s %.3f' % (t.replace('template_', ''), st.mean(v)))
