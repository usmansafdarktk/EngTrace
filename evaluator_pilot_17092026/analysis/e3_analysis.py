"""E3 against E0, and the null baseline that says whether E3 measures anything.

NULL BASELINE. Score each trace against the milestones of a DIFFERENT item of the
SAME template - same question shape, same kind of numbers, different values. A
trace cannot legitimately reach another instance's intermediates, so any coverage
there is coincidence from tolerance and unit scaling. If that number is not near
zero, E3's real scores are inflated by the same amount and cannot be trusted.
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
import e3_milestones as e3         # noqa: E402


def load(d):
    """Rows of the evaluator's CURRENT config only.

    An earlier version kept the last row per trace whatever its config, and after
    the E0 re-seed that silently mixed rows from two configs into one table.
    """
    rows = [json.loads(l) for f in glob.glob(os.path.join(P, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


e0r = load('e0')
# E3 also covers the robustness cohort; compare on the gold five E0 scores.
e3r = {k: v for k, v in load('e3').items() if k in e0r}
items = {json.loads(l)['item_id']: json.loads(l) for l in open(os.path.join(P, 'slice', 'manifest.jsonl'), encoding='utf-8')}
M = ms.build_all(os.path.join(P, 'slice', 'manifest.jsonl'))
traces = {}
for f in glob.glob(os.path.join(P, 'traces', '*.jsonl')):
    k = os.path.basename(f)[:-6]
    for l in open(f, encoding='utf-8'):
        r = json.loads(l)
        if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
            traces[(k, r['item_id'])] = r

cov = lambda r: r['scores']['milestone_coverage']

print('E3 milestone coverage vs E0, per model')
print('%-16s %8s %8s %8s' % ('model', 'E3 cov', 'E0 F1', 'E0 FAC'))
by = defaultdict(list)
for k, r in e3r.items():
    by[k[0]].append(k)
for m in sorted(by, key=lambda m: -st.mean(cov(e3r[k]) for k in by[m])):
    ks = by[m]
    print('%-16s %8.3f %8.3f %8.3f' % (m, st.mean(cov(e3r[k]) for k in ks),
                                     st.mean(e0r[k]['scores']['recovered_f1'] for k in ks if k in e0r),
                                     st.mean(e0r[k]['scores']['final_answer_acc'] for k in ks if k in e0r)))

# ---- null baseline -------------------------------------------------------
by_t = defaultdict(list)
for iid, it in items.items():
    by_t[it['template_id']].append(iid)
null, real = [], []
for (mk, iid), tr in traces.items():
    if (mk, iid) not in e3r:
        continue
    real.append(cov(e3r[(mk, iid)]))
    for other in by_t[items[iid]['template_id']]:
        if other == iid:
            continue
        mset = M[other]['milestones']
        if not mset:
            continue
        hits = e3.reach(mset, tr['text'])
        null.append(sum(h['reached'] for h in hits) / len(hits))
print()
print('NULL BASELINE (trace vs a sibling item\'s milestones): mean coverage %.3f over %d pairings'
      % (st.mean(null), len(null)))
print('REAL                                                : mean coverage %.3f over %d traces'
      % (st.mean(real), len(real)))

# ---- does E3 agree with correctness? ---------------------------------------
def corr(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** .5
    return num / den if den else float('nan')

common = [k for k in e3r if k in e0r]
c3 = [cov(e3r[k]) for k in common]
f0 = [e0r[k]['scores']['recovered_f1'] for k in common]
fa = [e0r[k]['scores']['final_answer_acc'] for k in common]
print()
print('corr(E3 coverage, E0 reasoning F1)      = %.3f' % corr(c3, f0))
print('corr(E3 coverage, E0 final-answer acc)  = %.3f' % corr(c3, fa))

scaled = sum(e3r[k]['meta']['reached_by_unit_scaling'] for k in e3r)
reached = sum(e3r[k]['scores']['milestones_reached'] for k in e3r)
print()
print('milestones reached: %d, of which via unit scaling: %d (%.1f%%)' % (reached, scaled, 100 * scaled / max(reached, 1)))

print()
print('by answer type:  E3 cov | E0 F1')
bt = defaultdict(list)
for k in common:
    bt[e3r[k]['answer_type']].append(k)
for t in sorted(bt, key=lambda t: -st.mean(cov(e3r[k]) for k in bt[t])):
    ks = bt[t]
    print('  %-15s %.3f | %.3f  (n=%d)' % (t, st.mean(cov(e3r[k]) for k in ks),
                                          st.mean(e0r[k]['scores']['recovered_f1'] for k in ks), len(ks)))

# E0-F1 check: does E3 get rackett right where E0 scored 0/20?
rk = [k for k in common if k[1].startswith('rackett')]
print()
print('rackett_equation_volume (E0 final-answer acc 0 of 20 by parser bug):')
print('  E3 mean coverage %.3f; traces reaching V_sat: %d of %d' % (
    st.mean(cov(e3r[k]) for k in rk),
    sum(any(h['id'] == 'V_sat' and h['reached'] for h in e3r[k]['meta']['milestones']) for k in rk), len(rk)))
