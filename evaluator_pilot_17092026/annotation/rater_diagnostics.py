"""Which rater, if any, is holding a branch's agreement down.

    python evaluator_pilot_17092026/annotation/rater_diagnostics.py [--labels DIR]
    ... --drop che-3.1 --add ../version_1/labels/superseded/che-3.jsonl   # what a bad set looks like

A low branch kappa has two very different causes, and they need different fixes:

  ONE RATER IS OUT OF STEP   their label distribution differs from the other two, they
                             are the odd one out in most 2-1 splits, the adjudication
                             goes against them, and dropping them lifts the branch's
                             kappa a long way. That is a rater to re-run (che-3 was).

  THE BRANCH IS HARDER       all three disagree with each other about equally, no
                             single rater is the outlier, and leave-one-out barely
                             moves the kappa. Re-running one expert fixes nothing.

Per rater it reports: the share of steps given each label, mean pairwise Cohen kappa
with the other two, how often they are the odd one out in a 2-1 split, how often the
adjudication (the experts' own blind second look) ended up against their original
label, and their intra-rater agreement from the verification round.
"""
import argparse
import glob
import json
import os
import sys
from collections import Counter, defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
sys.path[:0] = [_HERE]
import score_against_labels as S  # noqa: E402
import verification_report as V  # noqa: E402

DEFAULT = os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels')


def by_branch(labels):
    out = defaultdict(lambda: defaultdict(dict))
    for r in labels:
        if r.get('for_truth'):
            out[r['branch']][r['code']][r['annotator']] = r
    return out


def fleiss_steps(codes_map, raters=None):
    rows = []
    for code, per in codes_map.items():
        rs = [r for a, r in per.items() if raters is None or a in raters]
        if len(rs) < 2:
            continue
        for i in range(rs[0]['n_steps']):
            rows.append(Counter(r['steps'][i]['label'] for r in rs))
    return S.fleiss(rows, S.STEP_CATS), len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=DEFAULT)
    ap.add_argument('--drop', action='append', default=[], help='annotator id to exclude')
    ap.add_argument('--add', action='append', default=[], help='extra label file to include')
    a = ap.parse_args()

    labels = [r for r in S.read_labels(a.labels) if r['annotator'] not in a.drop]
    for p in a.add:
        p = p if os.path.isabs(p) else os.path.join(a.labels, p)
        labels += [json.loads(l) for l in open(p, encoding='utf-8') if l.strip()]
    branches = by_branch(labels)
    relabels = V.read_relabels(a.labels)
    first = {(r['annotator'], r['code']): r for r in labels}
    cons = {}
    for f in sorted(glob.glob(os.path.join(a.labels, 'adjudication', '*-consensus.jsonl'))):
        for line in open(f, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                cons[(r['code'], r['index'])] = r

    print('labels: %s%s\n' % (os.path.relpath(a.labels, _PILOT),
                              (' (dropped %s, added %s)' % (','.join(a.drop), ','.join(a.add))
                               if a.drop or a.add else '')))
    for branch in sorted(branches):
        codes = branches[branch]
        raters = sorted({a_ for per in codes.values() for a_ in per})
        k_all, nsteps = fleiss_steps(codes)
        print('%s - %d traces, %d rater-steps, Fleiss kappa %.3f'
              % (branch.replace('_engineering', ''), len(codes), nsteps, k_all))
        print('  %-9s %6s %6s %6s %6s %8s %9s %10s %9s'
              % ('rater', 'corr', 'alt', 'inc', 'n/a', 'pair k', 'odd/splits', 'adj against', 'self k'))
        for r_ in raters:
            dist = Counter()
            pair, odd, splits, against = [], 0, 0, 0
            for code, per in codes.items():
                if r_ not in per:
                    continue
                mine = per[r_]
                for i in range(mine['n_steps']):
                    dist[mine['steps'][i]['label']] += 1
                    labs = {a2: p['steps'][i]['label'] for a2, p in per.items()}
                    if len(labs) == 3:
                        c = Counter(labs.values())
                        if len(c) == 2 and min(c.values()) == 1:
                            splits += 1
                            if c[labs[r_]] == 1:
                                odd += 1
                            adj = cons.get((code, i))
                            if adj and adj['consensus']['label'] != labs[r_]:
                                against += 1
                for a2, other in per.items():
                    if a2 == r_:
                        continue
                    pair.append(V.cohen([(mine['steps'][i]['label'], other['steps'][i]['label'])
                                         for i in range(mine['n_steps'])], S.STEP_CATS))
            tot = sum(dist.values())
            rl = [(first[(r_, c)], v) for (a2, c), v in relabels.items()
                  if a2 == r_ and (r_, c) in first]
            selfsteps = [(x['steps'][i]['label'], y['steps'][i]['label'])
                         for x, y in rl for i in range(x['n_steps'])]
            pk = [p for p in pair if p == p]
            print('  %-9s %6.3f %6.3f %6.3f %6.3f %8.3f %9s %10s %9s'
                  % (r_, dist['correct'] / tot, dist['alternative_correct'] / tot,
                     dist['incorrect'] / tot, dist['not_a_claim'] / tot,
                     sum(pk) / len(pk) if pk else float('nan'),
                     '%d/%d' % (odd, splits),
                     '%d/%d' % (against, splits) if splits else '-',
                     '%.3f' % V.cohen(selfsteps, S.STEP_CATS) if selfsteps else '-'))
        for r_ in raters:
            k_wo, _ = fleiss_steps(codes, raters=[x for x in raters if x != r_])
            print('  without %-9s the other two agree at Cohen kappa %.3f  (%+.3f)'
                  % (r_, k_wo, k_wo - k_all))
        print()


if __name__ == '__main__':
    main()
