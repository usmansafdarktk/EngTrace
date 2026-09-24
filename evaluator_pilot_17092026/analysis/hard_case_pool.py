"""X2 - the hard case: flawed reasoning behind a correct final answer.

    python evaluator_pilot_17092026/analysis/hard_case_pool.py [--labels DIR] [--top 60]

X1 could not run the comparison a reasoning evaluator exists for. Of the 228 traces
whose final answer the experts call correct, their holistic verdict calls 3 unsound -
too few to separate any two evaluators. This script does two things about that.

A. THE POOL ALREADY IN THE LABELS. The experts' STEP labels disagree with their own
   verdict: 56 of those correct-answer traces carry at least one step the experts mark
   incorrect. Scoring the evaluators on "does this trace contain an incorrect step,
   given the answer is right" is the hard-case comparison with 56 positives instead of
   3. It needs no new traces and no new labelling - only the decision that this, not
   the holistic verdict, is the trace-level target when the answer is correct.

B. CANDIDATES FOR A SECOND ROUND. To widen the pool, traces nobody has labelled
   (gemma-4-31b, qwen3.8-27b: the robustness cohort, same 60 items) are ranked by
   signals that cost nothing to compute - the arithmetic checker (E4's), the 72B PRM's
   minimum step reward (E2) and milestone coverage (E3). Every signal's precision is
   measured on the labelled traces first, so the yield of a labelling round is
   estimated rather than assumed. Candidates are written to
   analysis/out/hard_case_candidates.jsonl for build_tasks.py.

The deterministic answer check used for unlabelled traces is the one below, not E0's:
its agreement with the experts' verdict is printed so it is not taken on faith.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import glob
import json
import statistics as st
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]
import arith  # noqa: E402
import milestones  # noqa: E402
import score_against_labels as S  # noqa: E402
import x1_analysis as X  # noqa: E402

UNLABELLED = ('gemma-4-31b', 'qwen3.8-27b')
ANSWER_TOL = 0.02          # E0's final-answer tolerance, not the 0.5% display tolerance
OUT = _os.path.join(_ANALYSIS, 'out')


def trace_text():
    """Every generated trace that came back ok, keyed as the score rows are."""
    out = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        for line in open(f, encoding='utf-8'):
            r = json.loads(line)
            if r.get('ok') and r.get('text'):
                out[(r['model_key'], r['item_id'])] = r
    return out


def gold_final():
    """The gold's final value per item: the last number its solution writes."""
    out = {}
    for line in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8'):
        item = json.loads(line)
        nums = milestones.numbers(item['solution'])
        if nums:
            out[item['item_id']] = nums[-1]
    return out


def answer_ok(text, gold):
    """Does the trace end on the gold's final value? Last 3 numbers, unit factors allowed."""
    nums = milestones.numbers(text)[-3:]
    return bool(nums) and milestones.scaled_match(gold, nums, ANSWER_TOL) is not None


def arith_flags(text):
    """Failed arithmetic claims, as E4's checker counts them, under both of its
    rules: the digit rule E4 now scores on (D-097) and the 1% tolerance it shipped
    with. The two are reported side by side because they select different traces."""
    rep = arith.check(text)
    return rep.checked - rep.consistent_digit, rep.checked - rep.consistent, rep.checked


def describe(truth, keyfile, codes):
    print('\nA. THE POOL IN THE LABELS - correct-answer traces carrying an incorrect step')
    bad = [c for c in codes if any(s['label'] == 'incorrect' for s in truth[c]['steps'])]
    nsteps = sum(sum(s['label'] == 'incorrect' for s in truth[c]['steps']) for c in bad)
    verdict = Counter(truth[c]['reasoning_sound'] for c in bad)
    print('  %d of %d correct-answer traces, %d incorrect steps in them' % (len(bad), len(codes), nsteps))
    print('  the experts\' own verdict on those %d traces: %s' % (len(bad), dict(verdict)))
    print('  (the verdict calls %d unsound; the step labels call %d flawed - the gap this set closes)'
          % (sum(1 for c in codes if truth[c]['reasoning_sound'] == 'no'), len(bad)))
    for field, label in (('branch', 'by branch'), (None, 'by model')):
        if field:
            c = Counter(truth[x][field] for x in bad)
            tot = Counter(truth[x][field] for x in codes)
        else:
            c = Counter(keyfile[x][0] for x in bad)
            tot = Counter(keyfile[x][0] for x in codes)
        print('  %-10s %s' % (label, ', '.join('%s %d/%d' % (k.replace('_engineering', ''), v, tot[k])
                                               for k, v in sorted(c.items()))))
    err = Counter(s['error_type'] for c in bad for s in truth[c]['steps']
                  if s['label'] == 'incorrect' and s['error_type'])
    print('  error types  %s' % ', '.join('%s %d' % kv for kv in err.most_common()))
    return bad


def table(rows_by_code, truth, keyfile, codes, y, title):
    """Every evaluator against a yes/no target, with bootstrap intervals and E0 deltas."""
    npos = sum(y.values())
    print('\n%s  (%d traces: %d flawed, %d clean)' % (title, len(codes), npos, len(codes) - npos))
    if npos < 5 or len(codes) - npos < 5:
        print('  too few of one class'); return
    score = {}
    for name, d, f in X.HEADS:
        score[name] = {c: rows_by_code[d].get(keyfile[c], {}).get('scores', {}).get(f) for c in codes}
    score['baseline: E0 answer check'] = {
        c: rows_by_code['e0'].get(keyfile[c], {}).get('scores', {}).get('final_answer_acc') for c in codes}
    # `n`: how many of these traces the score exists for (E4's arithmetic score is
    # None where the checker found nothing to check). The E0 delta uses the overlap.
    print('  %-32s %7s  %-15s %5s  %-22s' % ('evaluator', 'AUROC', '95% CI', 'n', 'minus E0 (95% CI)'))
    for name in list(score):
        ok = [c for c in codes if score[name][c] is not None]
        # a flawed trace should score LOW, so AUROC is computed on "clean" as the positive
        a = X.auroc([(score[name][c], not y[c]) for c in ok])
        lo, hi = X.boot(ok, lambda ks: X.auroc([(score[name][c], not y[c]) for c in ks]))
        both = [c for c in ok if score['E0'][c] is not None]
        if name == 'E0':
            d = '-'
        else:
            dd = (X.auroc([(score[name][c], not y[c]) for c in both])
                  - X.auroc([(score['E0'][c], not y[c]) for c in both]))
            dlo, dhi = X.boot(both, lambda ks: X.auroc([(score[name][c], not y[c]) for c in ks])
                              - X.auroc([(score['E0'][c], not y[c]) for c in ks]))
            d = '%+.3f (%+.3f, %+.3f)%s' % (dd, dlo, dhi, '' if dlo <= 0 <= dhi else ' *')
        print('  %-32s %7.3f  (%.3f, %.3f) %5d  %s' % (name, a, lo, hi, len(ok), d))
    print('  * = the difference from E0 excludes zero at 95%')


def signals(rows_by_code, truth, keyfile, texts, codes, y):
    """What each cheap signal is worth as a filter: precision against the base rate."""
    print('\nB1. CHEAP SIGNALS, measured on the labelled correct-answer traces')
    base = sum(y.values()) / len(codes)
    print('  base rate (label at random): %.1f%% of traces carry an incorrect step' % (100 * base))
    print('  %-34s %6s %6s %6s %7s' % ('filter', 'n', 'prec', 'rec', 'lift'))
    rows = []
    for c in codes:
        t = texts.get(keyfile[c])
        bad, bad1, checked = arith_flags(t['text']) if t else (0, 0, 0)
        e2 = rows_by_code['e2'].get(keyfile[c], {}).get('scores', {})
        e3 = rows_by_code['e3'].get(keyfile[c], {}).get('scores', {})
        rows.append((c, {'arith_failed': bad, 'arith_failed_tol1': bad1, 'arith_checked': checked,
                         'prm_min': e2.get('qwen72_min'), 'prm_frac': e2.get('qwen72_frac_ok'),
                         'e3_cov': e3.get('milestone_coverage')}))
    tests = [('arithmetic 1%: >=1 failed claim', lambda s: s['arith_failed_tol1'] >= 1),
             ('arithmetic 1%: >=2 failed claims', lambda s: s['arith_failed_tol1'] >= 2),
             ('arithmetic digit: >=1 failed claim', lambda s: s['arith_failed'] >= 1),
             ('arithmetic digit: >=2 failed claims', lambda s: s['arith_failed'] >= 2),
             ('PRM 72B: min reward < 0.20', lambda s: s['prm_min'] is not None and s['prm_min'] < 0.20),
             ('PRM 72B: min reward < 0.05', lambda s: s['prm_min'] is not None and s['prm_min'] < 0.05),
             ('PRM 72B: frac_ok < 0.90', lambda s: s['prm_frac'] is not None and s['prm_frac'] < 0.90),
             ('E3: milestone coverage < 1.0', lambda s: s['e3_cov'] is not None and s['e3_cov'] < 1.0),
             ('arith digit OR PRM min < 0.20', lambda s: s['arith_failed'] >= 1 or
              (s['prm_min'] is not None and s['prm_min'] < 0.20)),
             ('arith digit AND PRM min < 0.20', lambda s: s['arith_failed'] >= 1 and
              (s['prm_min'] is not None and s['prm_min'] < 0.20))]
    for name, fn in tests:
        sel = [c for c, s in rows if fn(s)]
        if not sel:
            print('  %-34s %6d %6s %6s %7s' % (name, 0, '-', '-', '-')); continue
        tp = sum(y[c] for c in sel)
        prec, rec = tp / len(sel), tp / sum(y.values())
        print('  %-34s %6d %6.3f %6.3f %6.2fx' % (name, len(sel), prec, rec, prec / base))
    return dict(rows)


def mine(texts, golds, rows_by_code, labelled, top):
    """Rank the unlabelled traces. Answer right by the deterministic check, reasoning suspect."""
    print('\nB2. CANDIDATES from the traces nobody has labelled')
    cand = []
    for (model, item), t in texts.items():
        if model not in UNLABELLED or (model, item) in labelled or item not in golds:
            continue
        if not answer_ok(t['text'], golds[item]):
            continue
        bad, bad1, checked = arith_flags(t['text'])
        e2 = rows_by_code['e2'].get((model, item), {}).get('scores', {})
        e3 = rows_by_code['e3'].get((model, item), {}).get('scores', {})
        prm, cov = e2.get('qwen72_min'), e3.get('milestone_coverage')
        rank = (bad > 0) * 2 + (prm is not None and prm < 0.20) + (cov is not None and cov < 1.0)
        cand.append({'model_key': model, 'item_id': item, 'branch': t['branch'], 'level': t['level'],
                     'arith_failed': bad, 'arith_failed_tol1': bad1, 'arith_checked': checked,
                     'prm_min': prm, 'milestone_coverage': cov, 'rank': rank})
    cand.sort(key=lambda r: (-r['rank'], r['prm_min'] if r['prm_min'] is not None else 1.0))
    per = Counter(r['model_key'] for r in cand)
    print('  answer correct by the deterministic check: %s' % dict(per))
    print('  %-14s %6s %5s %5s %6s' % ('signals fired', 'n', 'chem', 'adv', 'PRM?'))
    for r_ in sorted({c['rank'] for c in cand}, reverse=True):
        sel = [c for c in cand if c['rank'] == r_]
        print('  %-14s %6d %5d %5d %6d' % (r_, len(sel),
                                           sum(c['branch'] == 'chemical_engineering' for c in sel),
                                           sum(c['level'] == 'Advanced' for c in sel),
                                           sum(c['prm_min'] is not None for c in sel)))
    _os.makedirs(OUT, exist_ok=True)
    path = _os.path.join(OUT, 'hard_case_candidates.jsonl')
    with open(path, 'w', encoding='utf-8') as fh:
        for c in cand[:top]:
            fh.write(json.dumps(c) + '\n')
    print('  wrote the top %d to %s' % (min(top, len(cand)), _os.path.relpath(path, _PILOT)))
    return cand


def check_answer_check(texts, golds, truth, keyfile):
    """How far the deterministic answer check can be trusted: against the experts' verdict."""
    agree = n = 0
    for c, t in truth.items():
        tr = texts.get(keyfile[c])
        if not tr or keyfile[c][1] not in golds or t['final_answer'] is None:
            continue
        n += 1
        agree += answer_ok(tr['text'], golds[keyfile[c][1]]) == (t['final_answer'] == 'correct')
    e0 = S.current('e0')
    e0_agree = sum(1 for c, t in truth.items()
                   if e0.get(keyfile[c]) and
                   (e0[keyfile[c]]['scores']['final_answer_acc'] == 1.0) == (t['final_answer'] == 'correct'))
    print('\n  answer check against the experts\' verdict: this script %.3f (n=%d), E0 %.3f'
          % (agree / n, n, e0_agree / len(truth)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    ap.add_argument('--top', type=int, default=60)
    a = ap.parse_args()
    truth, _ = S.build_truth(S.read_labels(a.labels), S.read_consensus(a.labels))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    rows_by_code = {d: S.current(d) for d in {d for _, d, _ in X.HEADS} | {'e0', 'e2', 'e3'}}
    texts, golds = trace_text(), gold_final()
    right = [c for c in truth if truth[c]['final_answer'] == 'correct']
    print('%d traces labelled, %d with a correct final answer' % (len(truth), len(right)))

    describe(truth, keyfile, right)
    y = {c: any(s['label'] == 'incorrect' for s in truth[c]['steps']) for c in right}
    table(rows_by_code, truth, keyfile, right, y,
          'A2. THE HARD CASE, step-label target - which evaluator finds the flawed trace')
    yall = {c: any(s['label'] == 'incorrect' for s in truth[c]['steps']) for c in truth}
    table(rows_by_code, truth, keyfile, list(truth), yall,
          'A3. The same target over all 300 traces, for comparison')
    signals(rows_by_code, truth, keyfile, texts, right, y)
    check_answer_check(texts, golds, truth, keyfile)
    mine(texts, golds, rows_by_code, {keyfile[c] for c in truth}, a.top)


if __name__ == '__main__':
    main()
