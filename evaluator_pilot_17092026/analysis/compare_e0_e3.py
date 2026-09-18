"""E0 (shared-seed) vs E3 v2, on CURRENT rows only.

The earlier table mixed E0 rows from two configs because the loader kept the last
row per trace regardless of config. Here each evaluator's rows are filtered to the
config of its most recent run, and the script refuses to report unless every
(model, item) has exactly one current row.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import glob
import json
import os
import statistics as st
from collections import Counter, defaultdict

P = _PILOT_DIR


def current(d):
    rows = [json.loads(l) for f in glob.glob(os.path.join(P, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    cur = {}
    for r in rows:
        if r['config_sha256'] == latest:
            cur[(r['model_key'], r['item_id'])] = r
    return latest, cur


s0, e0 = current('e0')
s3, e3 = current('e3')
# E3 also scores the robustness cohort (gemma). The comparison is over the gold
# five, the only models E0 scores, so restrict E3 to E0's traces.
e3 = {k: v for k, v in e3.items() if k in e0}
assert len(e0) == 300 and len(e3) == 300, (len(e0), len(e3))
print('E0 config %s (300 rows)   E3 config %s (300 rows)' % (s0[:12], s3[:12]))
cov = lambda r: r['scores']['milestone_coverage']
f1 = lambda r: r['scores']['recovered_f1']
fac = lambda r: r['scores']['final_answer_acc']


def corr(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** .5
    return num / den if den else float('nan')


keys = sorted(e0)
print()
print('%-16s %8s %8s %8s %12s' % ('model', 'E0 FAC', 'E0 F1', 'E3 cov', 'E0 judged'))
by = defaultdict(list)
for k in keys:
    by[k[0]].append(k)
for m in sorted(by, key=lambda m: -st.mean(fac(e0[k]) for k in by[m])):
    ks = by[m]
    print('%-16s %8.3f %8.3f %8.3f %12d' % (m, st.mean(fac(e0[k]) for k in ks), st.mean(f1(e0[k]) for k in ks),
                                          st.mean(cov(e3[k]) for k in ks), sum(bool(e0[k]['calls']) for k in ks)))
print('%-16s %8.3f %8.3f %8.3f %12d' % ('ALL', st.mean(fac(e0[k]) for k in keys), st.mean(f1(e0[k]) for k in keys),
                                      st.mean(cov(e3[k]) for k in keys), sum(bool(e0[k]['calls']) for k in keys)))
print()
a, b, c = [fac(e0[k]) for k in keys], [f1(e0[k]) for k in keys], [cov(e3[k]) for k in keys]
print('corr(E0 F1, E0 FAC) = %.3f' % corr(b, a))
print('corr(E3 cov, E0 FAC) = %.3f' % corr(c, a))
print('corr(E3 cov, E0 F1)  = %.3f' % corr(c, b))

# Where they disagree most: E3 says the reasoning is there, E0 says it is not, and vice versa
print()
dis = sorted(keys, key=lambda k: cov(e3[k]) - f1(e0[k]))
print('E0 high, E3 low (E0 F1 - E3 cov largest):')
for k in dis[:4]:
    print('  %-16s %-32s E0 F1 %.2f  E3 %.2f  FAC %d' % (k[0], k[1], f1(e0[k]), cov(e3[k]), fac(e0[k])))
print('E3 high, E0 low:')
for k in dis[-4:]:
    print('  %-16s %-32s E0 F1 %.2f  E3 %.2f  FAC %d' % (k[0], k[1], f1(e0[k]), cov(e3[k]), fac(e0[k])))

# Does the E0 re-run's own sample now match what e0_3j's would be? (seed check)
old = [json.loads(l) for f in glob.glob(os.path.join(P, 'scores', 'e0', '*.jsonl')) for l in open(f, encoding='utf-8')]
seeds = Counter(r['config_sha256'][:12] for r in old)
print()
print('e0 rows on disk by config: %s' % dict(seeds))
