"""X4 - what a corrected final-answer check is worth.

    python evaluator_pilot_17092026/analysis/answer_check.py [--labels DIR]

RESULTS_X1 Finding 1: the final answer decides the trace-level verdict almost entirely
(the experts' own answer verdict predicts their soundness verdict at AUROC 0.974), and
E0's answer check agrees with the experts on only 76% of traces. So the check is the
largest single accuracy gain in the pilot, and it needs no model. This measures how much
of that gap simple, deterministic fixes close, against the experts' verdict on 300
traces.

The two defects it has to fix (FINDINGS E0-F1, E0-F2):

  E0-F1  the gold's answer line restates the inputs - "...at 283.81 K is 113.55 cm3/mol"
         - and the first number on it is read as the gold answer.
  E0-F2  8 of 15 templates put no number on the answer line at all, because the answer
         is a word: turbulent, laminar, safe.

Variants, each a strictly bigger fix than the last:

  E0        what the published framework scores today
  last      the trace's last numbers against the gold's last number
  segment   read the trace's ANSWER SEGMENT (after its final Answer marker), and take
            the gold value from the gold's last computed number, not its answer line
  +words    segment, plus: when the gold's answer is qualitative, require the word
  +parts    +words, plus: a multi-part answer (a)/(b) is correct only if every part is
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import json
import re
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]
import milestones  # noqa: E402
import score_against_labels as S  # noqa: E402

TOL = 0.02
ANSWER = re.compile(r'(?i)(?:#+\s*final\s+answer|\*{0,2}answer\b)')
PART = re.compile(r'(?i)\(?\b([ab])\)\s*[:.\-]?')
WORDS = ('turbulent', 'laminar', 'transitional', 'safe', 'unsafe', 'yes', 'no', 'feasible',
         'infeasible', 'stable', 'unstable', 'acceptable', 'adequate', 'inadequate',
         'subsonic', 'supersonic', 'saturated', 'superheated', 'critical', 'subcritical')


def segment(text):
    """Everything after the trace's last Answer marker - what a reader would call the answer."""
    hits = list(ANSWER.finditer(text))
    return text[hits[-1].end():] if hits else text[-200:]


def gold_answer(item):
    """The gold's final value, and any qualitative word its answer states."""
    nums = milestones.numbers(item['solution'])
    seg = segment(item['solution']).lower()
    return (nums[-1] if nums else None), [w for w in WORDS if re.search(r'\b%s\b' % w, seg)]


def check(text, gold_value, gold_words, mode):
    seg = segment(text)
    nums = milestones.numbers(seg if mode != 'last' else text)
    nums = nums[-6:] if mode != 'last' else nums[-3:]
    ok = gold_value is not None and bool(nums) and milestones.scaled_match(gold_value, nums, TOL) is not None
    if mode in ('+words', '+parts') and gold_words:
        ok = ok and all(re.search(r'\b%s\b' % w, seg.lower()) for w in gold_words)
        if not [n for n in nums] and gold_words:      # purely qualitative answer
            ok = all(re.search(r'\b%s\b' % w, seg.lower()) for w in gold_words)
    if mode == '+parts' and len(set(m.group(1).lower() for m in PART.finditer(seg))) > 1:
        ok = ok and (not gold_words or all(re.search(r'\b%s\b' % w, seg.lower()) for w in gold_words))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    a = ap.parse_args()
    truth, _ = S.build_truth(S.read_labels(a.labels), S.read_consensus(a.labels))
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
    e0 = S.current('e0')

    modes = ('E0', 'last', 'segment', '+words', '+parts')
    res = {m: Counter() for m in modes}
    per_template = defaultdict(Counter)
    for code, t in truth.items():
        key = keyfile[code]
        item, text = items.get(key[1]), texts.get(key)
        if not item or not text or t['final_answer'] is None:
            continue
        want = t['final_answer'] == 'correct'
        gv, gw = gold_answer(item)
        for m in modes:
            if m == 'E0':
                row = e0.get(key)
                if not row:
                    continue
                got = row['scores']['final_answer_acc'] == 1.0
            else:
                got = check(text, gv, gw, m)
            res[m]['ok' if got == want else ('false_correct' if got else 'false_wrong')] += 1
            if m == '+parts' and got != want:
                per_template[item['template_id']][key[0]] += 1
    print('Against the experts\' final-answer verdict, %d traces\n' % len(truth))
    print('  %-9s %8s %14s %14s   %s' % ('check', 'accuracy', 'says correct/', 'says wrong/', 'what it adds'))
    print('  %-9s %8s %14s %14s' % ('', '', 'experts wrong', 'experts right'))
    notes = {'E0': 'the published framework today', 'last': 'no answer-line parsing at all',
             'segment': 'fixes E0-F1 on both sides', '+words': 'fixes E0-F2 (qualitative answers)',
             '+parts': 'multi-part answers must be right in every part'}
    for m in modes:
        c = res[m]
        n = sum(c.values())
        print('  %-9s %8.3f %14d %14d   %s' % (m, c['ok'] / n, c['false_correct'], c['false_wrong'], notes[m]))
    print('\n  ceiling: the experts\' own verdict predicts their soundness verdict at AUROC 0.974;')
    print('  E0\'s check, used the same way, scores 0.812.')
    print('\n  WHAT THE CHECK DECIDES: answer accuracy per model, the benchmark\'s headline number')
    print('    %-16s %7s %10s %10s %10s' % ('model', 'traces', 'experts', 'E0', 'corrected'))
    per_model = defaultdict(Counter)
    for code, t in truth.items():
        key = keyfile[code]
        item, text = items.get(key[1]), texts.get(key)
        if not item or not text or t['final_answer'] is None:
            continue
        gv, gw = gold_answer(item)
        row = e0.get(key)
        c = per_model[key[0]]
        c['n'] += 1
        c['expert'] += t['final_answer'] == 'correct'
        c['fixed'] += check(text, gv, gw, '+words')
        if row:
            c['e0n'] += 1
            c['e0'] += row['scores']['final_answer_acc'] == 1.0
    order = {}
    for m, c in sorted(per_model.items(), key=lambda kv: -kv[1]['expert'] / kv[1]['n']):
        print('    %-16s %7d %10.3f %10.3f %10.3f'
              % (m, c['n'], c['expert'] / c['n'], c['e0'] / c['e0n'] if c['e0n'] else float('nan'),
                 c['fixed'] / c['n']))
        order[m] = (c['expert'] / c['n'], c['e0'] / c['e0n'] if c['e0n'] else -1, c['fixed'] / c['n'])
    tot = Counter()
    for c in per_model.values():
        tot.update(c)
    print('    %-16s %7d %10.3f %10.3f %10.3f' % ('ALL', tot['n'], tot['expert'] / tot['n'],
                                                  tot['e0'] / tot['e0n'], tot['fixed'] / tot['n']))
    for j, name in ((0, 'experts'), (1, 'E0'), (2, 'corrected')):
        print('    %-10s ranking: %s' % (name, ' > '.join(sorted(order, key=lambda m: -order[m][j]))))

    if per_template:
        print('\n  where the best variant still disagrees, by template:')
        for tid, c in sorted(per_template.items(), key=lambda kv: -sum(kv[1].values()))[:6]:
            print('    %-38s %d  (%s)' % (tid[:38], sum(c.values()),
                                          ', '.join('%s %d' % kv for kv in c.most_common(3))))


if __name__ == '__main__':
    main()
