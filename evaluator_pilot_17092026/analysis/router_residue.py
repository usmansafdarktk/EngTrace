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

The two rates are two different routers, and only one of them has been measured for what
it catches. The planted-defect probe (planted_judges.py, router_planted.py) asked a judge
about ONE step per call, which is the per-step design. The per-trace figure assumes a
trace's residue steps are batched into one prompt, which is cheaper and has not been
tested: more context per call, and more steps competing for the judge's attention.

The full run is --models 11, the roster D-110 fixed.
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
ITEMS = 2250               # the benchmark: 150 templates x 15 instances
MODELS = 11                # the roster, D-110


def verifiable(step):
    """Does this step contain a claim whose displayed value can be recomputed?"""
    for c in arith.check(step).claims:
        if c.left_value and c.right_value:
            return True
    return False


def measure(labels_dir, models=MODELS):
    """The residue over the labelled 300 and what routing it would cost, as one dict."""
    truth, _ = S.build_truth(S.read_labels(labels_dir), S.read_consensus(labels_dir))
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

    full = models * ITEMS
    with_res = sum(1 for x in per_trace if x)
    return {'traces': len(per_trace), 'steps': steps_total, 'residue_steps': residue_total,
            'traces_with_residue': with_res, 'per_trace': per_trace, 'by_model': dict(by_model),
            'labels': dict(label_of_residue), 'models': models, 'full_run_traces': full,
            'per_call': PER_CALL,
            'calls_batched': full * with_res / len(per_trace),
            'calls_per_step': full * st.mean(per_trace),
            'cost_batched': full * with_res / len(per_trace) * PER_CALL,
            'cost_per_step': full * st.mean(per_trace) * PER_CALL}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    ap.add_argument('--models', type=int, default=MODELS, help='roster size (D-110: 11)')
    a = ap.parse_args()
    m = measure(a.labels, a.models)
    per_trace, n, with_res = m['per_trace'], m['traces'], m['traces_with_residue']
    print('RESIDUE - steps a deterministic checker cannot verify (%d traces, %d steps)\n'
          % (n, m['steps']))
    print('  steps that are residue          %d of %d  (%.1f%%)'
          % (m['residue_steps'], m['steps'], 100 * m['residue_steps'] / m['steps']))
    print('  traces with at least one        %d of %d  (%.1f%%)' % (with_res, n, 100 * with_res / n))
    print('  residue steps per trace         median %.1f, mean %.1f, max %d'
          % (st.median(per_trace), st.mean(per_trace), max(per_trace)))
    print('\n  by model')
    for model, v in sorted(m['by_model'].items()):
        print('    %-16s %5.1f residue steps per trace, %3d%% of its traces have one'
              % (model, st.mean(v), round(100 * sum(1 for x in v if x) / len(v))))

    print('\n  WHAT THE EXPERTS SAID ABOUT THOSE STEPS')
    labels = Counter(m['labels'])
    tot = sum(labels.values())
    for k, v in labels.most_common():
        print('    %-20s %5d  (%.1f%%)' % (k, v, 100 * v / tot))
    print('    so %.1f%% of the residue is a step the experts call incorrect - the fraction of'
          % (100 * labels['incorrect'] / tot))
    print('    a judge\'s attention that would land on a real error')

    print('\n  COST AT FULL SCALE (%d models x %d items = %d traces, MiMo at $%.4f a call)'
          % (m['models'], ITEMS, m['full_run_traces'], PER_CALL))
    print('    one call per trace with residue:   %.0f calls  =  $%.0f   (batched; untested)'
          % (m['calls_batched'], m['cost_batched']))
    print('    one call per residue STEP:         %.0f calls  =  $%.0f   (the design the probe measured)'
          % (m['calls_per_step'], m['cost_per_step']))
    print('\n  The framework already batches a trace\'s mismatched steps into a single prompt, and a')
    print('  router would do the same - but the detection rates in router_planted.py come from')
    print('  asking about one step at a time, so they describe the per-step design until a')
    print('  batched prompt is measured.')


if __name__ == '__main__':
    main()
