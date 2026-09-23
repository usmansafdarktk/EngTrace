"""X1: the evaluators against the expert labels - with uncertainty, a baseline, and the hard case.

    python evaluator_pilot_17092026/analysis/x1_analysis.py [--labels DIR]

annotation/score_against_labels.py gives the headline comparison. This script asks
what those numbers can and cannot support:

1. UNCERTAINTY. Only a minority of the 300 traces are unsound, so AUROC differences of
   a few points may be noise. 95% bootstrap intervals (2,000 resamples of traces), and
   each evaluator's difference from E0 with its own interval, on the same resamples.
2. A BASELINE. A reasoning evaluator must beat "is the final answer right?". Two
   versions: E0's own final-answer check (what a benchmark would have), and the
   experts' final-answer verdict (the best any answer check could do).
3. THE HARD CASE. Among traces whose final answer the experts call correct, which
   evaluator finds the unsound reasoning? That is what a reasoning evaluator is for;
   on wrong-answer traces an answer check already does the job.
4. PER MODEL. Do the evaluators order the five models as the experts do?
5. SENSITIVITY. Everything again without chemical engineering, whose experts agree
   least (step kappa 0.31).
6. STEP LEVEL. E2's rewards against the experts' step labels, over all traces and
   inside correct-answer traces - the one place reasoning is judged apart from the
   answer. Bootstrap resamples traces, since steps within a trace are not independent.
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

sys.path[:0] = [_os.path.join(_PILOT, 'annotation')]
import score_against_labels as S  # noqa: E402

HEADS = [('E0', 'e0', 'recovered_f1'), ('E0-3J', 'e0_3j', 'recovered_f1'),
         ('E1', 'e1', 'recovered_f1'), ('E2 (72B frac)', 'e2', 'e2'),
         ('E2 (72B min)', 'e2', 'qwen72_min'), ('E3', 'e3', 'milestone_coverage'),
         ('E4', 'e4', 'e4_coverage'), ('E5', 'e5', 'e5_strict')]
MODELS = ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b')
B = 2000


def auroc(pairs):
    """Mann-Whitney AUROC from average ranks - equal to the pairwise count with ties
    as a half, in O(n log n), which the 2,000-resample bootstrap needs."""
    n1 = sum(1 for _, y in pairs if y)
    n0 = len(pairs) - n1
    if not n1 or not n0:
        return float('nan')
    order = sorted(range(len(pairs)), key=lambda i: pairs[i][0])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pairs[order[j + 1]][0] == pairs[order[i]][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    r1 = sum(r for r, (_, y) in zip(ranks, pairs) if y)
    return (r1 - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot(keys, fn, seed=1):
    rnd = random.Random(seed)
    out = []
    for _ in range(B):
        out.append(fn([keys[rnd.randrange(len(keys))] for _ in keys]))
    out = sorted(x for x in out if x == x)
    return out[int(.025 * len(out))], out[int(.975 * len(out)) - 1]


def table(rows_by_code, truth, keyfile, title, codes=None):
    codes = sorted(codes if codes is not None else truth)
    y = {c: truth[c]['reasoning_sound'] == 'yes' for c in codes}
    npos, nneg = sum(y.values()), len(codes) - sum(y.values())
    print('\n%s  (%d traces: %d sound, %d unsound)' % (title, len(codes), npos, nneg))
    if nneg < 5 or npos < 5:
        print('  too few of one class to estimate AUROC')
        return
    score = {}
    for name, d, f in HEADS:
        score[name] = {c: rows_by_code[d].get(keyfile[c], {}).get('scores', {}).get(f) for c in codes}
    e0 = rows_by_code['e0']
    score['baseline: E0 answer check'] = {c: e0.get(keyfile[c], {}).get('scores', {}).get('final_answer_acc') for c in codes}
    score['baseline: expert answer verdict'] = {c: 1.0 if truth[c]['final_answer'] == 'correct' else
                                                (0.5 if truth[c]['final_answer'] == 'partial' else 0.0) for c in codes}
    print('  %-32s %7s  %-15s  %-22s' % ('evaluator', 'AUROC', '95% CI', 'minus E0 (95% CI)'))
    for name in list(score):
        ok = [c for c in codes if score[name][c] is not None]
        a = auroc([(score[name][c], y[c]) for c in ok])
        lo, hi = boot(ok, lambda ks: auroc([(score[name][c], y[c]) for c in ks]))
        both = [c for c in ok if score['E0'][c] is not None]
        if name == 'E0':
            d = '-'
        else:
            dd = auroc([(score[name][c], y[c]) for c in both]) - auroc([(score['E0'][c], y[c]) for c in both])
            dlo, dhi = boot(both, lambda ks: auroc([(score[name][c], y[c]) for c in ks])
                            - auroc([(score['E0'][c], y[c]) for c in ks]))
            d = '%+.3f (%+.3f, %+.3f)%s' % (dd, dlo, dhi, '' if dlo <= 0 <= dhi else ' *')
        print('  %-32s %7.3f  (%.3f, %.3f)   %s' % (name, a, lo, hi, d))


def per_model(rows_by_code, truth, keyfile, codes):
    print('\nPER MODEL - share of traces the experts call sound, and each evaluator\'s mean')
    by = defaultdict(list)
    for c in codes:
        by[keyfile[c][0]].append(c)
    cols = [('experts: sound', lambda c: 1.0 if truth[c]['reasoning_sound'] == 'yes' else 0.0),
            ('experts: steps ok', lambda c: st.mean([s['label'] != 'incorrect' for s in truth[c]['steps']]))]
    for name, d, f in HEADS:
        cols.append((name, lambda c, d=d, f=f: rows_by_code[d].get(keyfile[c], {}).get('scores', {}).get(f)))
    print('  %-16s' % 'model' + ''.join('%15s' % n[:14] for n, _ in cols))
    means = {}
    for m in MODELS:
        vals = []
        for n, fn in cols:
            xs = [fn(c) for c in by[m] if fn(c) is not None]
            vals.append(st.mean(xs) if xs else float('nan'))
        means[m] = vals
        print('  %-16s' % m + ''.join('%15.3f' % v for v in vals))
    order = lambda j: [m for m in sorted(MODELS, key=lambda m: -means[m][j])]
    print('  expert order (sound): %s' % ' > '.join(order(0)))
    for j, (n, _) in enumerate(cols[2:], 2):
        print('  %-18s %s' % (n + ':', ' > '.join(order(j))))


def steps(truth, keyfile, codes, title):
    """E2 at the level it works at: each step's reward against the experts' label.
    Bootstrap resamples TRACES (steps within a trace are not independent)."""
    e2 = S.current('e2')
    per = {}
    for c in codes:
        row = e2.get(keyfile[c])
        if not row:
            continue
        per[c] = {p: row['meta']['prm'][p]['step_probs'] for p in S.PRMS}
    bad = sum(1 for c in per for s in truth[c]['steps'] if s['label'] == 'incorrect')
    good = sum(1 for c in per for s in truth[c]['steps'] if s['label'] in ('correct', 'alternative_correct'))
    print('\n%s  (%d traces, %d steps the experts call incorrect, %d correct)' % (title, len(per), bad, good))
    if bad < 5:
        print('  too few incorrect steps')
        return
    print('  %-8s %7s  %-15s %6s %6s %6s' % ('PRM', 'AUROC', '95% CI', 'prec', 'rec', 'F1'))
    for p in S.PRMS:
        def pairs(ks):
            out = []
            for c in ks:
                pr = per[c][p]
                if pr is None or len(pr) != len(truth[c]['steps']):
                    continue
                for s, v in zip(truth[c]['steps'], pr):
                    if s['label'] in ('correct', 'alternative_correct', 'incorrect'):
                        out.append((1 - v, s['label'] == 'incorrect'))   # low reward = flagged
            return out
        ks = sorted(per)
        a = auroc(pairs(ks))
        lo, hi = boot(ks, lambda kk: auroc(pairs(kk)))
        pp = pairs(ks)
        tp = sum(1 for v, y in pp if y and v > 0.5)
        fp = sum(1 for v, y in pp if not y and v > 0.5)
        fn = sum(1 for v, y in pp if y and v <= 0.5)
        pr_, rc, f1 = S.prf(tp, fp, fn)
        print('  %-8s %7.3f  (%.3f, %.3f)  %6.3f %6.3f %6.3f' % (p, a, lo, hi, pr_, rc, f1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    a = ap.parse_args()
    labels = S.read_labels(a.labels)
    truth, disputed = S.build_truth(labels, S.read_consensus(a.labels))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    rows_by_code = {d: S.current(d) for d in {d for _, d, _ in HEADS}}
    print('expert ground truth: %d traces, disputed %s' % (len(truth), dict(disputed)))

    table(rows_by_code, truth, keyfile, '1-2. ALL TRACES - separating sound from unsound reasoning')
    right = [c for c in truth if truth[c]['final_answer'] == 'correct']
    table(rows_by_code, truth, keyfile, '3. THE HARD CASE - final answer correct, per the experts', right)
    steps(truth, keyfile, list(truth), '6. STEP LEVEL - E2 against the experts, all traces')
    steps(truth, keyfile, right, '6b. STEP LEVEL - inside correct-answer traces (the hard case, per step)')
    per_model(rows_by_code, truth, keyfile, list(truth))
    nochem = [c for c in truth if truth[c]['branch'] != 'chemical_engineering']
    table(rows_by_code, truth, keyfile, '5. SENSITIVITY - without chemical engineering', nochem)
    table(rows_by_code, truth, keyfile, '5b. SENSITIVITY - hard case without chemical engineering',
          [c for c in right if c in nochem])
    print('\n* = the difference from E0 excludes zero at 95%')


if __name__ == '__main__':
    main()
