"""E4: what the arithmetic check adds over E3, and whether its flags are real.

Reports, on current-config rows only:
  - per model: E3 coverage, E4 coverage, verified coverage, arithmetic consistency
  - how many milestones E4 classes as contradicted, and on which models
  - the parse coverage behind every consistency figure (a checker that parses 5%
    of a trace and reports it consistent is not evidence of anything)
  - the contradicted cases themselves, so each can be read and judged: a genuine
    right-number-wrong-reason trace, or a checker false positive
  - correlation of each score with E0's final-answer correctness
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)

import glob
import json
import statistics as st
from collections import defaultdict

P = _PILOT_DIR


def load(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(P, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def corr(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** .5
    return num / den if den else float('nan')


e4 = load('e4')
e0 = load('e0')
S = lambda r, k: r['scores'][k]

by = defaultdict(list)
for k in e4:
    by[k[0]].append(k)

print('%-16s %6s %6s %9s %11s %9s %8s' % ('model', 'E3', 'E4', 'verified', 'arith cons', 'claims/tr', 'contra'))
for m in sorted(by, key=lambda m: -st.mean(S(e4[k], 'e4_coverage') for k in by[m])):
    ks = by[m]
    ac = [S(e4[k], 'arith_consistency') for k in ks if S(e4[k], 'arith_consistency') is not None]
    print('%-16s %6.3f %6.3f %9.3f %11s %9.1f %8d' % (
        m, st.mean(S(e4[k], 'e3_coverage') for k in ks), st.mean(S(e4[k], 'e4_coverage') for k in ks),
        st.mean(S(e4[k], 'verified_coverage') for k in ks),
        '%.3f' % st.mean(ac) if ac else '-',
        st.mean(e4[k]['meta']['claims_checked'] for k in ks),
        sum(S(e4[k], 'milestones_contradicted') for k in ks)))

seg = defaultdict(int)
for r in e4.values():
    for kk, v in r['meta']['segments'].items():
        seg[kk] += v
tot = sum(seg.values())
print()
print('segments across all traces: evaluated %d (%.0f%%), symbolic %d (%.0f%%), unparseable %d (%.1f%%)'
      % (seg['num'], 100 * seg['num'] / tot, seg['symbolic'], 100 * seg['symbolic'] / tot,
         seg['unparseable'], 100 * seg['unparseable'] / tot))
status = defaultdict(int)
for r in e4.values():
    for h in r['meta']['milestones']:
        status[h['status']] += 1
print('milestone status across all traces:', dict(status))

gold = [k for k in e4 if k in e0]
print()
fac = [S(e0[k], 'final_answer_acc') for k in gold]
for name in ('e3_coverage', 'e4_coverage', 'verified_coverage'):
    print('corr(%-18s, E0 final-answer correctness) = %.3f' % (name, corr([S(e4[k], name) for k in gold], fac)))

print()
print('CONTRADICTED milestones - read each one:')
n = 0
for k in sorted(e4):
    for h in e4[k]['meta']['milestones']:
        if h['status'] == 'contradicted' and n < 12:
            ev = h['evidence'][0]
            print('  %-15s %-32s %-14s milestone %.6g' % (k[0], k[1], h['id'], h['value']))
            print('        shown: %s = %s   (evaluates %s vs %s)' % (
                ev['left'][:55], ev['right'][:25],
                ['%.5g' % x for x in ev['left_value']], ['%.5g' % x for x in ev['right_value']]))
            n += 1
