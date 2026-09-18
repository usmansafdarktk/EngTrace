"""E3 separation vs tolerance and unit scaling: real coverage, null coverage, and
the correlation with final-answer correctness under each setting.

The choice is made on PRINCIPLE, then measured, not picked to maximise a number
over 300 traces: the question is only whether the null drops while real coverage
survives. Separation = real - null.
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

items = {json.loads(l)['item_id']: json.loads(l) for l in open(os.path.join(P, 'slice', 'manifest.jsonl'), encoding='utf-8')}
M = ms.build_all(os.path.join(P, 'slice', 'manifest.jsonl'))
e0 = {}
for f in glob.glob(os.path.join(P, 'scores', 'e0', '*.jsonl')):
    for l in open(f, encoding='utf-8'):
        r = json.loads(l)
        if not r.get('error'):
            e0[(r['model_key'], r['item_id'])] = r
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


def hit(v, nums, tol, scaled):
    if scaled:
        return ms.scaled_match(v, nums, tol) is not None
    return any(ms.close(n, v, tol) for n in nums)


def corr(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** .5
    return num / den if den else float('nan')


print('%-9s %-8s %7s %7s %10s %12s' % ('tol', 'scaling', 'real', 'null', 'separation', 'corr w/ FAC'))
for tol in (0.02, 0.01, 0.005, 0.002):
    for scaled in (True, False):
        real, null, fac = [], [], []
        for (mk, iid), tr in traces.items():
            nums = ms.numbers(tr['text'])
            mset = M[iid]['milestones']
            if mset:
                real.append(sum(hit(m['value'], nums, tol, scaled) for m in mset) / len(mset))
                fac.append(e0[(mk, iid)]['scores']['final_answer_acc'] if (mk, iid) in e0 else 0)
            own = [m['value'] for m in mset] + ms.numbers(items[iid]['question'])
            for other in by_t[items[iid]['template_id']]:
                if other == iid:
                    continue
                d = [m for m in M[other]['milestones'] if ms.scaled_match(m['value'], own, ms.DISPLAY_TOL) is None]
                if d:
                    null.append(sum(hit(m['value'], nums, tol, scaled) for m in d) / len(d))
        print('%-9s %-8s %7.3f %7.3f %10.3f %12.3f' % (
            '%.1f%%' % (tol * 100), 'all' if scaled else 'none', st.mean(real), st.mean(null),
            st.mean(real) - st.mean(null), corr(real, fac)))
