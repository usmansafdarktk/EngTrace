"""The verification and adjudication rounds: how consistent the experts are, and what
the adjudication changed.

    python evaluator_pilot_17092026/annotation/verification_report.py [--labels DIR]

VERIFICATION (intra-rater). Each expert re-labelled about 10% of their own traces,
blind to their first pass. Agreement between the two passes is that expert's own
consistency - the ceiling on what any agreement number between different experts can
mean. Reported as raw agreement and Cohen's kappa on step labels, plus the verdict and
final-answer agreement, per expert and per branch.

ADJUDICATION. Each step where a branch's three experts split 2-1 went back to all
three, blind and shuffled, and the majority of the second pass is the consensus label.
This reports how many steps that settled, how often the reviewers came back unanimous,
and how many ground-truth step labels the consensus actually changes - the second is
the number that matters, since a consensus that only confirms the majority changes
nothing downstream.
"""
import argparse
import glob
import json
import os
import statistics as st
import sys
from collections import Counter, defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
sys.path[:0] = [_HERE]
import score_against_labels as S  # noqa: E402

DEFAULT = os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels')


def cohen(pairs, cats):
    """Cohen's kappa for two passes over the same items."""
    n = len(pairs)
    if not n:
        return float('nan')
    po = sum(a == b for a, b in pairs) / n
    pe = sum((sum(a == k for a, _ in pairs) / n) * (sum(b == k for _, b in pairs) / n) for k in cats)
    return (po - pe) / (1 - pe) if pe < 1 else float('nan')


def read_relabels(labels_dir):
    out = {}
    for f in sorted(glob.glob(os.path.join(labels_dir, 'verification', '*', '*-relabels.jsonl'))):
        for line in open(f, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                out[(r['annotator'], r['code'])] = r
    return out


def verification(labels_dir):
    first = {(r['annotator'], r['code']): r for r in S.read_labels(labels_dir)}
    again = read_relabels(labels_dir)
    if not again:
        print('no verification round in %s' % labels_dir)
        return
    print('VERIFICATION - each expert against their own first pass')
    print('  %-8s %7s %7s %8s %7s %8s %8s %9s' % ('expert', 'traces', 'steps', 'step agr',
                                                  'kappa', 'verdict', 'answer', 'min/trace'))
    per_branch = defaultdict(list)
    rows = []
    current = {a for a, _ in first}
    for (ann, code), r2 in sorted(again.items()):
        if ann not in current:
            continue          # a re-label belonging to a set that has since been superseded
        r1 = first.get((ann, code))
        if not r1:
            print('  %s %s: no first pass to compare' % (ann, code))
            continue
        rows.append((ann, r1, r2))
    for ann in sorted({a for a, _, _ in rows}):
        mine = [(r1, r2) for a, r1, r2 in rows if a == ann]
        steps = [(s1['label'], s2['label']) for r1, r2 in mine
                 for s1, s2 in zip(r1['steps'], r2['steps'])]
        verdict = [(r1['reasoning_sound'], r2['reasoning_sound']) for r1, r2 in mine]
        answer = [(r1['final_answer'], r2['final_answer']) for r1, r2 in mine]
        agr = sum(a == b for a, b in steps) / len(steps)
        k = cohen(steps, S.STEP_CATS)
        secs = [r2['seconds'] for r1, r2 in mine if r2.get('seconds')]
        print('  %-8s %7d %7d %8.3f %7.3f %8.3f %8.3f %9.1f'
              % (ann, len(mine), len(steps), agr, k,
                 sum(a == b for a, b in verdict) / len(verdict),
                 sum(a == b for a, b in answer) / len(answer),
                 st.median(secs) / 60 if secs else float('nan')))
        per_branch[mine[0][0]['branch']].append((steps, verdict))
    print('  %-8s %7s %7s %8s %7s' % ('branch', 'experts', 'steps', 'step agr', 'kappa'))
    allsteps = []
    for br, lots in sorted(per_branch.items()):
        steps = [p for s, _ in lots for p in s]
        allsteps += steps
        print('  %-8s %7d %7d %8.3f %7.3f' % (br.replace('_engineering', ''), len(lots), len(steps),
                                              sum(a == b for a, b in steps) / len(steps),
                                              cohen(steps, S.STEP_CATS)))
    print('  %-8s %7d %7d %8.3f %7.3f' % ('POOLED', len(rows), len(allsteps),
                                          sum(a == b for a, b in allsteps) / len(allsteps),
                                          cohen(allsteps, S.STEP_CATS)))
    disagree = Counter((a, b) for a, b in allsteps if a != b)
    if disagree:
        print('  where the two passes differ: %s'
              % ', '.join('%s->%s %d' % (a, b, n) for (a, b), n in disagree.most_common(6)))


def adjudication(labels_dir):
    cons = S.read_consensus(labels_dir)
    if not cons:
        print('\nno adjudication round in %s' % labels_dir)
        return
    labels = S.read_labels(labels_dir)
    majority_only, _ = S.build_truth(labels)
    settled, _ = S.build_truth(labels, cons)
    rows = [json.loads(l) for f in sorted(glob.glob(os.path.join(labels_dir, 'adjudication', '*-consensus.jsonl')))
            for l in open(f, encoding='utf-8') if l.strip()]
    print('\nADJUDICATION - the steps the branch experts split 2-1 on')
    print('  %d split steps reviewed, in %d traces' % (len(rows), len({r['code'] for r in rows})))
    print('  second pass came back: %s' % dict(Counter(r['agreement'] for r in rows)))
    by_branch = Counter(r['code'][:0] or 'all' for r in rows)
    changed = same = was_tied = 0
    for (code, i), c in cons.items():
        if code not in majority_only:
            continue
        before = majority_only[code]['steps'][i]['label']
        after = settled[code]['steps'][i]['label']
        if before is None:
            was_tied += 1
        elif before != after:
            changed += 1
        else:
            same += 1
    print('  effect on ground truth: %d step labels changed, %d confirmed the majority, '
          '%d had no majority before' % (changed, same, was_tied))
    print('  consensus labels: %s' % dict(Counter(c['label'] for c in cons.values())))
    withnote = sum(1 for c in cons.values() if c.get('note'))
    print('  consensus rows carrying a reason: %d of %d' % (withnote, len(cons)))
    bad_before = sum(1 for t in majority_only.values() for s in t['steps'] if s['label'] == 'incorrect')
    bad_after = sum(1 for t in settled.values() for s in t['steps'] if s['label'] == 'incorrect')
    tie_before = sum(1 for t in majority_only.values() for s in t['steps'] if s['label'] is None)
    tie_after = sum(1 for t in settled.values() for s in t['steps'] if s['label'] is None)
    print('  steps called incorrect: %d before, %d after; steps left with no label: %d before, %d after'
          % (bad_before, bad_after, tie_before, tie_after))
    print('  (by_branch check: %s)' % dict(by_branch) if False else '', end='')


def reasons(labels_dir):
    """The other thing this round added: a reason on every step called incorrect."""
    labels = S.read_labels(labels_dir)
    bad = [s for r in labels for s in r['steps'] if s['label'] == 'incorrect']
    withnote = sum(1 for s in bad if (s.get('note') or '').strip())
    print('\nREASONS - steps an expert called incorrect')
    print('  %d such steps in the submitted labels, %d carry a written reason (%.1f%%)'
          % (len(bad), withnote, 100 * withnote / len(bad) if bad else float('nan')))
    per = defaultdict(lambda: [0, 0])
    for r in labels:
        for s in r['steps']:
            if s['label'] == 'incorrect':
                per[r['annotator']][0] += 1
                per[r['annotator']][1] += bool((s.get('note') or '').strip())
    short = [(a, n, w) for a, (n, w) in sorted(per.items()) if w < n]
    if short:
        print('  still missing a reason: %s'
              % ', '.join('%s %d/%d' % (a, n - w, n) for a, n, w in short))
    else:
        print('  every expert gave a reason on every step they called incorrect')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=DEFAULT)
    a = ap.parse_args()
    print('labels: %s\n' % os.path.relpath(a.labels, _PILOT))
    verification(a.labels)
    adjudication(a.labels)
    reasons(a.labels)


if __name__ == '__main__':
    main()
