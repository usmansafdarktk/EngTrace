"""E0 test-retest: how much does E0 move between two runs of the same code?

The two E0 runs differ in exactly one configured thing, the seed of the 20%
wrong-answer sample. But the judges are also stochastic (GPT-5 takes no
temperature; the framework sets none for Claude). So a changed score has two
possible causes, and they can be separated:

  SAMPLING   the trace was judged in one run and not the other
  JUDGES     the trace was judged in BOTH runs and still scored differently

Any change on a trace Tier 1 decided alone (judged in neither run) would mean
something deterministic moved, which should be impossible - it is checked too.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import glob
import json
import os
import statistics as st
from collections import defaultdict

P = _PILOT_DIR
rows = [json.loads(l) for f in glob.glob(os.path.join(P, 'scores', 'e0', '*.jsonl'))
        for l in open(f, encoding='utf-8')]
by_cfg = defaultdict(dict)
for r in rows:
    if not r.get('error'):
        by_cfg[r['config_sha256']][(r['model_key'], r['item_id'])] = r
(ca, A), (cb, B) = sorted(by_cfg.items(), key=lambda kv: min(r['ts'] for r in kv[1].values()))
print('run 1 %s   run 2 %s' % (ca[:12], cb[:12]))
f1 = lambda r: r['scores']['recovered_f1']
jud = lambda r: bool(r['calls'])

cats = defaultdict(list)
for k in A:
    a, b = A[k], B[k]
    d = f1(b) - f1(a)
    if jud(a) and jud(b):
        cats['judged in both'].append(d)
    elif jud(a) or jud(b):
        cats['judged in one only (sampling)'].append(d)
    else:
        cats['judged in neither (Tier 1 only)'].append(d)

print()
print('%-34s %5s %9s %14s' % ('', 'n', 'changed', 'mean |delta|'))
for c, ds in cats.items():
    ch = [d for d in ds if abs(d) > 1e-9]
    print('%-34s %5d %9d %14.3f' % (c, len(ds), len(ch), st.mean(abs(d) for d in ch) if ch else 0))

both = [k for k in A if jud(A[k]) and jud(B[k])]
ch = [k for k in both if abs(f1(A[k]) - f1(B[k])) > 1e-9]
print()
print('JUDGE NON-DETERMINISM: %d of %d traces judged in both runs scored differently (%.0f%%)'
      % (len(ch), len(both), 100 * len(ch) / len(both)))

# per model: how much of each model's mean moved, and from which cause
print()
print('%-16s %8s %8s %8s   %s' % ('model', 'run 1', 'run 2', 'delta', 'from judges / from sampling'))
bym = defaultdict(list)
for k in A:
    bym[k[0]].append(k)
for m in sorted(bym):
    ks = bym[m]
    r1 = st.mean(f1(A[k]) for k in ks)
    r2 = st.mean(f1(B[k]) for k in ks)
    dj = sum(f1(B[k]) - f1(A[k]) for k in ks if jud(A[k]) and jud(B[k])) / len(ks)
    ds = sum(f1(B[k]) - f1(A[k]) for k in ks if jud(A[k]) != jud(B[k])) / len(ks)
    print('%-16s %8.3f %8.3f %+8.3f   %+.3f / %+.3f' % (m, r1, r2, r2 - r1, dj, ds))
