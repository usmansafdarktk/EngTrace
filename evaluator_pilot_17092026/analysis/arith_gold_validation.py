"""The checker's gold validation, run under BOTH of its rules.

    python evaluator_pilot_17092026/analysis/arith_gold_validation.py [--rule digit|tol1|both]

RESULTS_E4 step 2: gold arithmetic is correct by construction, so every claim the
checker marks inconsistent on a gold solution is a checker bug, not an error in the
gold. That validation was run at the 1% tolerance and reached 100%. D-097 replaces
the tolerance with the experts' rule - a displayed value must be a correct rounding
at the precision it shows - which is far stricter, so the validation has to be run
again: a rule that flags gold is flagging its own parsing.

This prints, over the 60 frozen gold solutions:
  - claims checked, and the consistency rate under each rule
  - the digit rule in three forms, so the two fixes it needed to pass gold are each
    worth a number rather than an assertion:
      bare      the displayed precision alone, as D-097 measured it
      +units    a conversion is compared after its unit factor
      +shown    and a result is not held to a precision its own displayed operands
                cannot pin down (`arith.shown_uncertainty`) - what E4 now ships
  - every claim still rejected, with the value shown, the value recomputed and the
    place value that judged it, so each one can be read rather than counted

Nothing here calls a model or the network; it reads slice/manifest.jsonl only.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import json
import sys
from collections import Counter

sys.path[:0] = [_os.path.join(_PILOT, 'evaluators')]
import arith  # noqa: E402


def gold_items():
    path = _os.path.join(_PILOT, 'slice', 'manifest.jsonl')
    return [json.loads(l) for l in open(path, encoding='utf-8')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rule', choices=('tol1', 'digit', 'both'), default='both')
    ap.add_argument('--show', type=int, default=40, help='how many flagged claims to print')
    a = ap.parse_args()

    items = gold_items()
    claims, flags = [], {'tol1': [], 'digit': [], 'digit_bare': [], 'digit_units': []}
    seg = Counter()
    for it in items:
        rep = arith.check(it['solution'])
        seg.update(rep.segments)
        for c in rep.claims:
            claims.append(c)
            if not c.ok:
                flags['tol1'].append((it['item_id'], c))
            if not c.ok_digit:
                flags['digit'].append((it['item_id'], c))
            if c.ulp is None:
                continue
            # The same rule without each fix, to price the fix.
            bare = any(abs(a - b) <= 0.5 * c.ulp * (1 + 1e-7)
                       for a in c.left_value for b in c.right_value)
            if not bare:
                flags['digit_bare'].append((it['item_id'], c))
            if not arith.agree_any_digit(c.left_value, c.right_value,
                                         c.left_unit, c.right_unit, c.ulp):
                flags['digit_units'].append((it['item_id'], c))

    n = len(claims)
    judged = sum(1 for c in claims if c.ulp is not None)
    print('%d gold solutions, %d arithmetic claims checked '
          '(%d of them with a readable displayed precision)' % (len(items), n, judged))
    print('segments: evaluated %d, symbolic %d, unparseable %d'
          % (seg['num'], seg['symbolic'], seg['unparseable']))
    print()
    print('  %-28s %8s %9s %12s' % ('rule', 'flagged', 'rate', 'consistent'))
    order = (('tol1', '1% relative (ARITH_TOL)'),
             ('digit_bare', 'digit rule, bare'),
             ('digit_units', 'digit rule, +units'),
             ('digit', 'digit rule, +shown (E4)'))
    for rule, label in order:
        if a.rule == 'tol1' and rule != 'tol1':
            continue
        if a.rule == 'digit' and rule == 'tol1':
            continue
        f = len(flags[rule])
        print('  %-28s %8d %8.2f%% %11.3f' % (label, f, 100.0 * f / n if n else 0.0,
                                              (n - f) / n if n else float('nan')))
    print()
    print('Anything above zero on gold is a checker bug: the gold arithmetic is correct')
    print('by construction. Read every flag before reporting any rate the rule produces.')

    for rule in ('tol1', 'digit_bare', 'digit'):
        if (a.rule == 'tol1') != (rule == 'tol1') or not flags[rule]:
            if a.rule != 'both':
                continue
        if not flags[rule]:
            continue
        print('\nFLAGGED ON GOLD BY %s - each one is a parsing bug to fix or exclude:' % rule)
        for item_id, c in flags[rule][:a.show]:
            print('  %-34s line %-4d %s = %s' % (item_id, c.line, c.left[:46], c.right[:26]))
            print('       recomputed %s   shown %s   last digit %s'
                  % (['%.8g' % x for x in c.left_value], ['%.8g' % x for x in c.right_value],
                     ('%.1e' % c.ulp) if c.ulp is not None else '-'))
        if len(flags[rule]) > a.show:
            print('  ... %d more' % (len(flags[rule]) - a.show))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
