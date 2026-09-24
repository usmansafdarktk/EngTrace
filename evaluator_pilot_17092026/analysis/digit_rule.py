"""X3 - the hard case, attacked with the experts' own rule instead of a model.

    python evaluator_pilot_17092026/analysis/digit_rule.py [--labels DIR]

Finding 5 of RESULTS_X1 says no evaluator detects a flawed step behind a correct final
answer. It also says what those flaws are: 164 of the 167 are calculation slips, not
conceptual errors. The annotation guide's rule for them is mechanical, and the
adjudication notes quote it constantly - *rounding is not an error, a wrong digit is*.

So apply that rule by machine. For every arithmetic claim a trace writes, recompute the
left side from the numbers the trace itself shows and ask whether the displayed right
side is a correct rounding of it at the precision shown. `0.92979^(2/3) = 0.95300` is
flagged because the value rounds to 0.95263 at five decimals; `= 0.9526` is not flagged,
and neither is any legitimate rounding.

This is E4's checker with its tolerance replaced by the displayed precision. E4 asked
whether a claim is within 1% - enough to catch a fabricated number, blind to every slip
the experts actually marked - and now scores on the rule this script measured (D-097).
Four rules are compared so every difference between them stays visible:

    tol1    relative tolerance 1%, what E4 shipped with
    tol01   relative tolerance 0.1%
    digit   the experts' rule: the last digit shown must be a correct rounding
    e4      the same rule AS THE EVALUATOR NOW SHIPS IT (arith.Claim.ok_digit)

`digit` is the rule as measured for D-097: the displayed precision alone, units
ignored. `e4` is what `evaluators/arith.py` applies now, and differs in two ways,
both found by re-running the gold validation (analysis/arith_gold_validation.py):
a unit conversion is compared after its unit factor, and a result is not held to a
precision its own displayed operands cannot pin down. Both loosen the rule, so `e4`
flags fewer claims than `digit`; the two lines below say what that costs.

Nothing here calls a model, so the cost is zero and the flag is auditable: it names the
claim, the value shown and the value recomputed.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import json
import re
import sys
from collections import Counter

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]
import arith  # noqa: E402
import e2_prm  # noqa: E402
import score_against_labels as S  # noqa: E402
import x1_analysis as X  # noqa: E402

NUM = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')
RULES = ('tol1', 'tol01', 'digit', 'e4')


def ulp(lit):
    """Place value of the last digit shown: '0.95300' -> 1e-5, '1.2e3' -> 100."""
    m = re.match(r'[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?$', lit.strip())
    if not m:
        return None
    return 10.0 ** ((int(m.group(3)) if m.group(3) else 0) - len(m.group(2)))


def flagged(text, rule):
    """Claims in this text the rule rejects, with what they said and what they compute to."""
    out = []
    for c in arith.check(text).claims:
        if rule == 'e4':                      # the evaluator's own verdict, unaltered
            if not c.ok_digit:
                out.append((c.left, c.right, c.left_value[0]))
            continue
        lit = NUM.search(c.right.replace(',', ''))
        if not lit or not c.left_value or not c.right_value:
            continue
        u = ulp(lit.group(0))
        if u is None:
            continue
        pairs = [(lv, rv) for lv in c.left_value for rv in c.right_value]
        if rule == 'tol1':
            ok = any(arith.agree(lv, rv, 0.01) for lv, rv in pairs)
        elif rule == 'tol01':
            ok = any(arith.agree(lv, rv, 0.001) for lv, rv in pairs)
        else:
            ok = any(abs(lv - rv) <= 0.5 * u * (1 + 1e-7) for lv, rv in pairs)
        if not ok:
            out.append((c.left, c.right, c.left_value[0]))
    return out


def load(labels_dir):
    labels = S.read_labels(labels_dir)
    truth, _ = S.build_truth(labels, S.read_consensus(labels_dir))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    texts = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for line in open(_os.path.join(_PILOT, 'traces', f), encoding='utf-8'):
                r = json.loads(line)
                if r.get('ok'):
                    texts[(r['model_key'], r['item_id'])] = r['text']
    return truth, keyfile, texts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    ap.add_argument('--examples', type=int, default=5)
    a = ap.parse_args()
    truth, keyfile, texts = load(a.labels)

    per_step = {r: {} for r in RULES}          # code -> [flag count per step]
    skipped = 0
    for code, t in truth.items():
        text = texts.get(keyfile[code])
        steps = e2_prm.steps_of(text) if text else []
        if len(steps) != len(t['steps']):
            skipped += 1
            continue
        for rule in RULES:
            per_step[rule][code] = [len(flagged(s, rule)) for s in steps]
    print('%d traces scored, %d skipped (no text or step-count mismatch)' % (len(per_step['digit']), skipped))

    right = [c for c in truth if truth[c]['final_answer'] == 'correct']
    for title, codes in (('ALL 300 TRACES', list(per_step['digit'])),
                         ('INSIDE CORRECT-ANSWER TRACES - the hard case', right)):
        print('\nSTEP LEVEL, %s' % title)
        print('  %-7s %5s %5s %5s %7s %7s %7s' % ('rule', 'tp', 'fp', 'fn', 'prec', 'rec', 'F1'))
        for rule in RULES:
            c = Counter()
            for code in codes:
                if code not in per_step[rule]:
                    continue
                for n, st in zip(per_step[rule][code], truth[code]['steps']):
                    if st['label'] == 'incorrect':
                        c['tp' if n else 'fn'] += 1
                    elif st['label'] in ('correct', 'alternative_correct'):
                        c['fp' if n else 'tn'] += 1
            p, r, f = S.prf(c['tp'], c['fp'], c['fn'])
            print('  %-7s %5d %5d %5d %7.3f %7.3f %7.3f' % (rule, c['tp'], c['fp'], c['fn'], p, r, f))

    print('\nTRACE LEVEL, correct-answer traces: does the rule find the flawed trace?')
    print('  %-7s %8s %7s %7s %7s %8s  %-15s' % ('rule', 'flagged', 'prec', 'rec', 'F1', 'AUROC', '95% CI'))
    codes = [c for c in right if c in per_step['digit']]
    y = {c: any(s['label'] == 'incorrect' for s in truth[c]['steps']) for c in codes}
    for rule in RULES:
        hit = {c: sum(per_step[rule][c]) for c in codes}
        tp = sum(1 for c in codes if hit[c] and y[c])
        fp = sum(1 for c in codes if hit[c] and not y[c])
        fn = sum(1 for c in codes if not hit[c] and y[c])
        p, r, f = S.prf(tp, fp, fn)
        auc = X.auroc([(-hit[c], not y[c]) for c in codes])       # more flags = flawed
        lo, hi = X.boot(codes, lambda ks: X.auroc([(-hit[c], not y[c]) for c in ks]))
        print('  %-7s %8d %7.3f %7.3f %7.3f %8.3f  (%.3f, %.3f)'
              % (rule, tp + fp, p, r, f, auc, lo, hi))
    print('  for comparison, the best evaluator on this target is E2 (72B min) at AUROC 0.584')

    print('\nWHAT THE DIGIT RULE FLAGS - first %d in correct-answer traces the experts also marked'
          % a.examples)
    shown = 0
    for code in right:
        if shown >= a.examples or code not in per_step['digit']:
            continue
        text = texts.get(keyfile[code])
        for i, (n, st) in enumerate(zip(per_step['digit'][code], truth[code]['steps'])):
            if n and st['label'] == 'incorrect' and shown < a.examples:
                left, shown_val, computed = flagged(e2_prm.steps_of(text)[i], 'digit')[0]
                print('  %s step %d: %s = %s, computes to %.6g' % (code, i, left[:44], shown_val[:22], computed))
                shown += 1


if __name__ == '__main__':
    main()
