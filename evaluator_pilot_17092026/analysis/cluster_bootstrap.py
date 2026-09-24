"""X5 - what survives clustering, and what this slice could ever have detected.

    python evaluator_pilot_17092026/analysis/cluster_bootstrap.py [--labels DIR] [-B 2000]

Every interval in RESULTS_X1 comes from resampling TRACES. That treats 300 traces as 300
independent observations, and they are not: they are 15 templates x 4 instances x 5 models.
Two traces of `rackett_equation_volume#0` and `#3` are the same problem with different
numbers, sharing a gold solution, a milestone structure and whatever quirk of phrasing an
evaluator reacts to. A reviewer who notices that will say the intervals are too narrow, and
they will be right. This script answers it with the data already in hand.

  DESIGN EFFECT   the same statistic, bootstrapped by resampling traces and by resampling
                  TEMPLATES (all 20 traces of a drawn template come with it). deff is the
                  ratio of the variances, and n_eff = 300 / deff is the number of
                  independent traces the slice is worth for that statistic.
  CLUSTER CIs     every headline comparison re-reported with template-level intervals.
  MDE             the smallest difference this design could detect at 80% power, given the
                  clustered standard error. A difference below it is not "no difference" -
                  it is a difference this slice was never able to see, and saying which is
                  the point of reporting power rather than leaving it implied.

Clusters are templates, not items: the item is the instance, the template is the thing the
next benchmark run would draw a new instance of. Item-level clustering is reported too,
because it is the weaker assumption and the difference between them is informative.

15 clusters is few. A cluster bootstrap on 15 draws is itself coarse, and its intervals
are lumpy rather than smooth - that is a property of the design, not of the estimator, and
it is exactly the thing being reported.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import json
import random
import statistics as st
import sys
from collections import defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _ANALYSIS]
import score_against_labels as S  # noqa: E402
import x1_analysis as X  # noqa: E402

Z = 1.959964 + 0.8416212          # 95% two-sided, 80% power


def load(labels_dir):
    truth, _ = S.build_truth(S.read_labels(labels_dir), S.read_consensus(labels_dir))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    tpl = {}
    for line in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8'):
        it = json.loads(line)
        tpl[it['item_id']] = it['template_id']
    rows = {d: S.current(d) for d in {d for _, d, _ in X.HEADS} | {'e0'}}
    return truth, keyfile, tpl, rows


def draws(codes, groups, B, seed=1):
    """B resamples. `groups` maps code -> cluster id; resampling draws clusters."""
    rnd = random.Random(seed)
    by = defaultdict(list)
    for c in codes:
        by[groups[c]].append(c)
    keys = sorted(by)
    for _ in range(B):
        out = []
        for _ in range(len(keys)):
            out.extend(by[keys[rnd.randrange(len(keys))]])
        yield out


def spread(fn, codes, groups, B):
    """(point estimate, se, 2.5th, 97.5th) of fn over resamples of `groups`."""
    vals = sorted(v for v in (fn(s) for s in draws(codes, groups, B)) if v == v)
    if len(vals) < 20:
        return float('nan'), float('nan'), float('nan'), float('nan')
    return (fn(codes), st.pstdev(vals), vals[int(.025 * len(vals))], vals[int(.975 * len(vals)) - 1])


def scores_for(rows, keyfile, codes, name, d, f):
    return {c: rows[d].get(keyfile[c], {}).get('scores', {}).get(f) for c in codes}


def report(title, codes, y, score, groupings, B):
    print('\n%s  (%d traces, %d positive)' % (title, len(codes), sum(y.values())))
    print('  %-24s %7s | %-22s | %-22s | %6s %7s'
          % ('evaluator', 'AUROC', 'trace bootstrap 95%', 'template bootstrap 95%', 'deff', 'n_eff'))
    base = None
    for name in score:
        ok = [c for c in codes if score[name][c] is not None]
        fn = lambda ks, nm=name: X.auroc([(score[nm][c], y[c]) for c in ks])
        pt, se_t, lo_t, hi_t = spread(fn, ok, groupings['trace'], B)
        _, se_c, lo_c, hi_c = spread(fn, ok, groupings['template'], B)
        deff = (se_c / se_t) ** 2 if se_t else float('nan')
        print('  %-24s %7.3f | (%.3f, %.3f)   %5.3f | (%.3f, %.3f)   %5.3f | %6.2f %7.0f'
              % (name, pt, lo_t, hi_t, se_t, lo_c, hi_c, se_c, deff, len(ok) / deff if deff else float('nan')))
        if name == 'E0':
            base = name
    if base:
        print('\n  difference from E0, with template-level intervals and the detectable difference')
        print('  %-24s %8s | %-22s | %10s %s'
              % ('evaluator', 'minus E0', 'template bootstrap 95%', 'MDE(80%)', 'verdict'))
        for name in score:
            if name == base:
                continue
            ok = [c for c in codes if score[name][c] is not None and score[base][c] is not None]
            fn = lambda ks, nm=name: (X.auroc([(score[nm][c], y[c]) for c in ks])
                                      - X.auroc([(score[base][c], y[c]) for c in ks]))
            pt, se_c, lo_c, hi_c = spread(fn, ok, groupings['template'], B)
            mde = Z * se_c
            if lo_c > 0 or hi_c < 0:
                verdict = 'separates'
            elif abs(pt) < mde:
                verdict = 'below what this design can detect'
            else:
                verdict = 'no separation'
            print('  %-24s %+8.3f | (%+.3f, %+.3f)        | %10.3f %s'
                  % (name, pt, lo_c, hi_c, mde, verdict))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    ap.add_argument('-B', type=int, default=2000)
    a = ap.parse_args()
    truth, keyfile, tpl, rows = load(a.labels)
    codes = sorted(truth)
    groupings = {'trace': {c: c for c in codes},
                 'item': {c: keyfile[c][1] for c in codes},
                 'template': {c: tpl[keyfile[c][1]] for c in codes}}
    print('%d traces, %d items, %d templates' % (
        len(codes), len({groupings['item'][c] for c in codes}), len({groupings['template'][c] for c in codes})))
    print('Clusters are templates: an instance is a re-draw of the same problem, so the '
          'independent unit is\nthe template, of which there are 15.')

    score = {n: scores_for(rows, keyfile, codes, n, d, f) for n, d, f in X.HEADS}
    score['baseline: expert answer verdict'] = {
        c: 1.0 if truth[c]['final_answer'] == 'correct' else (0.5 if truth[c]['final_answer'] == 'partial' else 0.0)
        for c in codes}
    y = {c: truth[c]['reasoning_sound'] == 'yes' for c in codes}
    report('1. TRACE LEVEL - separating sound from unsound reasoning', codes, y, score, groupings, a.B)

    right = [c for c in codes if truth[c]['final_answer'] == 'correct']
    flawed = {c: not any(s['label'] == 'incorrect' for s in truth[c]['steps']) for c in right}
    score2 = {n: {c: score[n][c] for c in right} for n in score if not n.startswith('baseline')}
    report('2. THE HARD CASE - correct answer, does the trace carry an incorrect step',
           right, flawed, score2, groupings, a.B)

    print('\n3. ITEM-LEVEL CLUSTERING, for comparison (the weaker assumption: 60 clusters)')
    print('  %-24s %7s | %-22s | %6s' % ('evaluator', 'AUROC', 'item bootstrap 95%', 'deff'))
    for name in ('E0', 'E2 (72B min)', 'E5'):
        ok = [c for c in codes if score[name][c] is not None]
        fn = lambda ks, nm=name: X.auroc([(score[nm][c], y[c]) for c in ks])
        pt, se_t, _, _ = spread(fn, ok, groupings['trace'], a.B)
        _, se_i, lo_i, hi_i = spread(fn, ok, groupings['item'], a.B)
        print('  %-24s %7.3f | (%.3f, %.3f)   %5.3f | %6.2f'
              % (name, pt, lo_i, hi_i, se_i, (se_i / se_t) ** 2 if se_t else float('nan')))

    print('\nMDE is the smallest difference detectable at 80% power and 95% confidence given the')
    print('clustered standard error: 2.80 x SE. It is the honest reading of every null result')
    print('here - the pilot can rule out differences larger than this, and not smaller ones.')


if __name__ == '__main__':
    main()
