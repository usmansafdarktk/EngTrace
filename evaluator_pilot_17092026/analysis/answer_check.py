"""X4 - the corrected final-answer check, measured against the expert verdicts.

    python evaluator_pilot_17092026/analysis/answer_check.py [--labels DIR]

RESULTS_X1 Finding 1: the final answer decides the trace-level verdict almost entirely
(the experts' own answer verdict predicts their soundness verdict at AUROC 0.974), and
E0's check agrees with the experts on only 76% of traces - erring almost always in one
direction, calling a correct answer wrong. E0_RERUN.md is the task; `evaluators/answer.py`
is the replacement; this is its evidence.

Reported here:
  AGREEMENT   E0 against the new check, on the same 300 traces, binary (correct or not)
              and three-way (correct / partial / incorrect - 19 traces are partial and a
              binary check cannot represent them at all).
  SPLIT-HALF  the relative tolerance is one fitted number, so it is chosen on one half of
              the traces and reported on the other, both ways round.
  PER MODEL   what the check does to the benchmark's headline accuracy and its ranking.
  RESIDUAL    every remaining disagreement, by template, so what is left is visible
              rather than summarised.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators')]
import answer as A  # noqa: E402
import score_against_labels as S  # noqa: E402

TOLS = (0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.02)


def load(labels_dir):
    truth, _ = S.build_truth(S.read_labels(labels_dir), S.read_consensus(labels_dir))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    items = {}
    for line in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8'):
        it = json.loads(line)
        items[it['item_id']] = it
    texts = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for line in open(_os.path.join(_PILOT, 'traces', f), encoding='utf-8'):
                r = json.loads(line)
                if r.get('ok'):
                    texts[(r['model_key'], r['item_id'])] = r['text']
    return truth, keyfile, items, texts


def verdicts(truth, keyfile, items, texts, tol=None):
    """The new check's verdict per trace, with the milestones the labels already carry."""
    out = {}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        ms = tuple(m['value'] for m in t['milestones'])
        out[code] = A.verdict(texts[key], items[key[1]], tol, ms)[0]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    a = ap.parse_args()
    truth, keyfile, items, texts = load(a.labels)
    e0 = S.current('e0')
    got = verdicts(truth, keyfile, items, texts)
    codes = sorted(got)
    exp = {c: truth[c]['final_answer'] for c in codes}

    three = sum(got[c] == exp[c] for c in codes)
    binary = sum((got[c] == 'correct') == (exp[c] == 'correct') for c in codes)
    e0ok = sum(1 for c in codes if e0.get(keyfile[c])
               and (e0[keyfile[c]]['scores']['final_answer_acc'] == 1.0) == (exp[c] == 'correct'))
    nonp = [c for c in codes if exp[c] != 'partial']
    print('AGREEMENT with the experts, %d traces\n' % len(codes))
    print('  %-34s %8s' % ('check', 'agrees'))
    print('  %-34s %8.3f' % ('E0, the published framework', e0ok / len(codes)))
    print('  %-34s %8.3f' % ('new check, correct-or-not', binary / len(codes)))
    print('  %-34s %8.3f' % ('new check, three-way', three / len(codes)))
    print('  %-34s %8.3f   (E0: %.3f)  <- E0_RERUN target 0.90'
          % ('new check, non-partial traces only',
             sum((got[c] == 'correct') == (exp[c] == 'correct') for c in nonp) / len(nonp),
             sum(1 for c in nonp if e0.get(keyfile[c])
                 and (e0[keyfile[c]]['scores']['final_answer_acc'] == 1.0) == (exp[c] == 'correct')) / len(nonp)))
    print('\n  confusion, experts -> new check')
    conf = Counter((exp[c], got[c]) for c in codes)
    for (e, g), n in conf.most_common():
        print('    %-10s -> %-10s %4d%s' % (e, g, n, '' if e == g else '   <- disagreement'))

    print('\nSPLIT-HALF - the tolerance is fitted, so it is reported off the half it was not fitted on')
    half = lambda c: int(hashlib.blake2b(c.encode(), digest_size=4).hexdigest(), 16) % 2
    acc = {}
    for tol in TOLS:
        g = verdicts(truth, keyfile, items, texts, tol)
        for h in (0, 1):
            ks = [c for c in codes if half(c) == h]
            acc[(tol, h)] = sum(g[c] == exp[c] for c in ks) / len(ks)
    for h in (0, 1):
        best = max(TOLS, key=lambda t: acc[(t, h)])
        print('  fit on half %d -> %.4f; held out on half %d: %.3f'
              % (h, best, 1 - h, acc[(best, 1 - h)]))
    print('  (the module ships REL = %.4f)' % A.REL)

    print('\nPER MODEL - answer accuracy, the benchmark\'s headline number')
    per = defaultdict(Counter)
    for c in codes:
        k = keyfile[c]
        per[k[0]]['n'] += 1
        per[k[0]]['exp'] += exp[c] == 'correct'
        per[k[0]]['new'] += got[c] == 'correct'
        if e0.get(k):
            per[k[0]]['e0'] += e0[k]['scores']['final_answer_acc'] == 1.0
    print('  %-16s %8s %8s %8s' % ('model', 'experts', 'E0', 'new'))
    tot = Counter()
    for m, c in sorted(per.items(), key=lambda kv: -kv[1]['exp'] / kv[1]['n']):
        tot.update(c)
        print('  %-16s %8.3f %8.3f %8.3f' % (m, c['exp'] / c['n'], c['e0'] / c['n'], c['new'] / c['n']))
    print('  %-16s %8.3f %8.3f %8.3f' % ('ALL', tot['exp'] / tot['n'], tot['e0'] / tot['n'],
                                         tot['new'] / tot['n']))
    for name, key in (('experts', 'exp'), ('E0', 'e0'), ('new', 'new')):
        print('  %-8s ranking: %s' % (name, ' > '.join(sorted(per, key=lambda m: -per[m][key] / per[m]['n']))))

    print('\nRESIDUAL disagreements, by template')
    left = defaultdict(list)
    for c in codes:
        if got[c] != exp[c]:
            left[items[keyfile[c][1]]['template_id']].append((keyfile[c][0], exp[c], got[c]))
    for tid, rows in sorted(left.items(), key=lambda kv: -len(kv[1])):
        kinds = Counter('%s->%s' % (e, g) for _, e, g in rows)
        print('  %-40s %3d  %s' % (tid.replace('template_', '')[:40], len(rows), dict(kinds)))
    print('  %d disagreements in total, of %d traces' % (sum(len(v) for v in left.values()), len(codes)))


if __name__ == '__main__':
    main()
