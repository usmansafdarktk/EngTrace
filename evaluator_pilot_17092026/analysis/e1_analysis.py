"""E1 against E0, and inter-judge agreement for both panels.

Reports, on each evaluator's CURRENT config only:
  1. E1 run health: rows, traces judged, judges called and parsed per trace,
     failures, truncations, cost, whether every reply came from the fetched store.
  2. Reasoning F1 per model under E1, E0 (shared seed) and E0-3J.
     E1 vs E0 share the same wrong-answer sample (D-082); E0-3J was run before the
     seed fix, so E1 vs E0-3J is compared on the traces both judged.
  3. INTER-JUDGE AGREEMENT - what E1 exists partly to report, and what reviewer
     9W1B asked for and E0 never reported. Each judge's per-step verdict is parsed
     from its raw reply with the framework's own logic; agreement is measured per
     (trace, step) both judges ruled on, on the 4 categories and on the binary
     "Alternative Correct" vs error. Pairwise % agreement, Cohen's kappa, and Fleiss'
     kappa across all three. Computed for E1's panel AND for E0-3J's original panel.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import glob
import itertools
import json
import statistics as st
import sys
from collections import Counter, defaultdict

sys.path.insert(0, _ANALYSIS)
from judge_probe import parse  # noqa: E402  - the framework's own reply parsing, copied

CATS = ('alternative correct', 'calculation error', 'conceptual error', 'other')


def current(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    if not rows:
        return None, {}
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return latest, {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def norm(c):
    c = str(c).lower()
    for k in CATS:
        if k.split()[0] in c:
            return k
    return 'other'


def verdicts(rows):
    """{(trace, step): {judge: category}} from each call's raw reply."""
    out = defaultdict(dict)
    for k, r in rows.items():
        for c in r['calls']:
            res = parse(c.get('text') or '') if c.get('ok') else None
            for x in res or []:
                if isinstance(x, dict) and isinstance(x.get('step_index'), int):
                    out[(k, x['step_index'])][c['model']] = norm(x.get('category'))
    return out


def kappa(pairs):
    n = len(pairs)
    if not n:
        return float('nan'), 0.0
    po = sum(a == b for a, b in pairs) / n
    ca, cb = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)
    pe = sum(ca[k] * cb[k] for k in set(ca) | set(cb)) / (n * n)
    return ((po - pe) / (1 - pe) if pe < 1 else float('nan')), po


def fleiss(items, cats):
    items = [v for v in items if len(v) >= 2]
    if not items:
        return float('nan')
    m = len(items[0])
    items = [v for v in items if len(v) == m]
    N = len(items)
    tot = Counter(c for v in items for c in v)
    pj = {c: tot[c] / (N * m) for c in cats}
    Pi = [(sum(Counter(v)[c] ** 2 for c in cats) - m) / (m * (m - 1)) for v in items]
    Pbar, Pe = st.mean(Pi), sum(p * p for p in pj.values())
    return (Pbar - Pe) / (1 - Pe) if Pe < 1 else float('nan')


def agreement(name, rows):
    v = verdicts(rows)
    judges = sorted({j for d in v.values() for j in d})
    print('\n%s panel: %s' % (name, ', '.join(judges)))
    dist = {j: Counter(d[j] for d in v.values() if j in d) for j in judges}
    for j in judges:
        n = sum(dist[j].values())
        print('  %-30s %4d step verdicts, Alternative Correct %.1f%%' % (
            j, n, 100 * dist[j]['alternative correct'] / max(n, 1)))
    print('  %-52s %9s %9s %10s %10s' % ('pair', 'n', 'agree-4', 'kappa-4', 'kappa-bin'))
    for a, b in itertools.combinations(judges, 2):
        pairs = [(d[a], d[b]) for d in v.values() if a in d and b in d]
        k4, p4 = kappa(pairs)
        kb, _ = kappa([(x == 'alternative correct', y == 'alternative correct') for x, y in pairs])
        print('  %-52s %9d %8.1f%% %10.3f %10.3f' % ('%s / %s' % (a, b), len(pairs), 100 * p4, k4, kb))
    full = [[d[j] for j in judges] for d in v.values() if all(j in d for j in judges)]
    print('  Fleiss kappa, all %d judges, %d steps all ruled on: 4-way %.3f, binary %.3f' % (
        len(judges), len(full), fleiss(full, CATS),
        fleiss([['alt' if c == 'alternative correct' else 'err' for c in x] for x in full], ('alt', 'err'))))


s1, e1 = current('e1')
s0, e0 = current('e0')
s3, e03 = current('e0_3j')
if not e1:
    raise SystemExit('no E1 rows yet')

# ---- 1. run health ------------------------------------------------------------
judged = [r for r in e1.values() if r['calls']]
calls = [c for r in judged for c in r['calls']]
print('E1 config %s: %d rows, %d judged' % (s1[:12], len(e1), len(judged)))
print('  judges called per judged trace: %s' % dict(Counter(len(r['meta']['judges_called']) for r in judged)))
print('  judges parsed per judged trace: %s' % dict(Counter(len(r['meta']['judges_parsed']) for r in judged)))
print('  judge failures: %d   truncated: %d   calls not from the fetched store: %d' % (
    sum(len(r['meta']['judge_failures']) for r in judged),
    sum(len(r['meta']['judges_truncated']) for r in judged),
    sum(1 for c in calls if not c.get('prefetched'))))
print('  cost: $%.2f   (%s)' % (sum(c.get('cost_usd', 0) for c in calls),
                               ', '.join('%s $%.2f' % (m, sum(c.get('cost_usd', 0) for c in calls if c['model'] == m))
                                         for m in sorted({c['model'] for c in calls}))))

# ---- 2. scores ------------------------------------------------------------------
f1 = lambda r: r['scores']['recovered_f1']
by = defaultdict(list)
for k in e1:
    by[k[0]].append(k)
both3 = {k for k in e1 if k in e03 and e1[k]['calls'] and e03[k]['calls']}
print('\nReasoning F1 per model')
print('%-16s %8s %8s   %s' % ('model', 'E0', 'E1', 'on traces E1 and E0-3J both judged: E0-3J -> E1'))
for m in sorted(by):
    ks = by[m]
    kk = [k for k in ks if k in both3]
    print('%-16s %8.3f %8.3f   n=%-3d %.3f -> %.3f' % (
        m, st.mean(f1(e0[k]) for k in ks if k in e0), st.mean(f1(e1[k]) for k in ks),
        len(kk), st.mean(f1(e03[k]) for k in kk) if kk else 0, st.mean(f1(e1[k]) for k in kk) if kk else 0))
same = [k for k in e1 if k in e0]
print('%-16s %8.3f %8.3f' % ('ALL', st.mean(f1(e0[k]) for k in same), st.mean(f1(e1[k]) for k in same)))
jd = [k for k in e1 if k in e0 and bool(e1[k]['calls']) == bool(e0[k]['calls'])]
print('traces with the same judged/not-judged status under E0 and E1: %d of %d (shared seed)' % (len(jd), len(same)))
moved = [k for k in both3 if abs(f1(e1[k]) - f1(e03[k])) > 1e-9]
print('E1 vs E0-3J on %d traces both judged: F1 changed on %d (%d up, %d down)' % (
    len(both3), len(moved), sum(f1(e1[k]) > f1(e03[k]) for k in moved), sum(f1(e1[k]) < f1(e03[k]) for k in moved)))

# ---- 3. agreement ---------------------------------------------------------------
agreement('E1', e1)
if e03:
    agreement('E0-3J (the original judges, all three connected)', e03)
