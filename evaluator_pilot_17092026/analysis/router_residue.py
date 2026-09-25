"""How much would a step router actually send to a judge, and what would that cost?

    python evaluator_pilot_17092026/analysis/router_residue.py

The router RESULTS_X1 recommendation 5 argues for verifies what it can deterministically and
sends only the residue to a judge. Its cost is therefore not a judgement call but a
measurement: how many steps carry nothing the checker can verify.

A step is VERIFIABLE here when `arith.check` finds at least one claim in it whose displayed
value can be recomputed from the numbers the step itself shows - exactly the condition under
which the digit rule scores 1.000 on planted defects and outside which it scores 0.000
(Finding 7). Everything else is residue: prose, a value stated with no working, a setup step.

Two rates matter, because they price differently:

  per trace   a router that batches a trace's residue into one prompt, as the framework's
              own Tribunal batches its mismatched steps, pays once per trace that has any
              residue at all.
  per step    the size of that prompt, and what a per-step router would pay.

Costs are quoted at MiMo-V2.5-Pro's measured rate, since that is the judge E5 uses and the
one that survives the judge/judged objection: $0.435 and $0.870 per million tokens, about
$0.0034 a call as JUDGE_SELECTION records it.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import json
import statistics as st
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]
import arith  # noqa: E402
import e2_prm  # noqa: E402
import score_against_labels as S  # noqa: E402

PER_CALL = 0.0034          # MiMo, measured over the pilot's E5 calls
FULL_RUN = 12 * 2250       # the roster x the benchmark


def verifiable(step):
    """Does this step contain a claim whose displayed value can be recomputed?"""
    for c in arith.check(step).claims:
        if c.left_value and c.right_value:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    a = ap.parse_args()
    truth, _ = S.build_truth(S.read_labels(a.labels), S.read_consensus(a.labels))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    texts = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for line in open(_os.path.join(_PILOT, 'traces', f), encoding='utf-8'):
                r = json.loads(line)
                if r.get('ok'):
                    texts[(r['model_key'], r['item_id'])] = r['text']

    per_trace, by_model = [], defaultdict(list)
    label_of_residue = Counter()
    steps_total = residue_total = 0
    for code, t in truth.items():
        text = texts.get(keyfile[code])
        if not text:
            continue
        steps = e2_prm.steps_of(text)
        if len(steps) != len(t['steps']):
            continue
        res = [i for i, s in enumerate(steps) if not verifiable(s)]
        steps_total += len(steps)
        residue_total += len(res)
        per_trace.append(len(res))
        by_model[keyfile[code][0]].append(len(res))
        for i in res:
            label_of_residue[t['steps'][i]['label']] += 1

    n = len(per_trace)
    with_res = sum(1 for x in per_trace if x)
    print('RESIDUE - steps a deterministic checker cannot verify (%d traces, %d steps)\n'
          % (n, steps_total))
    print('  steps that are residue          %d of %d  (%.1f%%)'
          % (residue_total, steps_total, 100 * residue_total / steps_total))
    print('  traces with at least one        %d of %d  (%.1f%%)' % (with_res, n, 100 * with_res / n))
    print('  residue steps per trace         median %.1f, mean %.1f, max %d'
          % (st.median(per_trace), st.mean(per_trace), max(per_trace)))
    print('\n  by model')
    for m, v in sorted(by_model.items()):
        print('    %-16s %5.1f residue steps per trace, %3d%% of its traces have one'
              % (m, st.mean(v), round(100 * sum(1 for x in v if x) / len(v))))

    print('\n  WHAT THE EXPERTS SAID ABOUT THOSE STEPS')
    tot = sum(label_of_residue.values())
    for k, v in label_of_residue.most_common():
        print('    %-20s %5d  (%.1f%%)' % (k, v, 100 * v / tot))
    bad = label_of_residue['incorrect']
    print('    so %.1f%% of the residue is a step the experts call incorrect - the fraction of'
          % (100 * bad / tot))
    print('    a judge\'s attention that would land on a real error')

    print('\n  COST AT FULL SCALE (%d traces, MiMo at $%.4f a call)' % (FULL_RUN, PER_CALL))
    rate = with_res / n
    print('    one call per trace with residue:   %.0f calls  =  $%.0f'
          % (FULL_RUN * rate, FULL_RUN * rate * PER_CALL))
    print('    one call per residue STEP:         %.0f calls  =  $%.0f'
          % (FULL_RUN * st.mean(per_trace), FULL_RUN * st.mean(per_trace) * PER_CALL))
    print('\n  The per-trace figure is the one to budget: the framework already batches a')
    print('  trace\'s mismatched steps into a single prompt, and a router would do the same.')


if __name__ == '__main__':
    main()
