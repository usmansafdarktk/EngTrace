"""Every value the pilot summary's figures draw, computed from the analyses, never typed in.

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/summary_numbers.py

The five figures in PILOT_SUMMARY.md used to carry their values as literals in
make_summary_figures.py, transcribed by hand from analysis output. That is how a router bar
nobody had measured (shown 1.000, caught 0.350) and a stale interval reached the summary.
This computes each value with the script that owns it and writes them to
figures/summary_numbers.json; make_summary_figures.py draws from that file and nothing else.

  trace     cluster_bootstrap.trace_scores and spread, templates resampled (B = 2,000), and
            cluster_bootstrap.margins_over: the experts' answer verdict minus every evaluator
  steps     digit_rule.flagged against the step labels inside correct-answer traces, and the
            72B PRM's step rewards at its 0.5 cut-off, read as x1_analysis.steps reads them
  planted   planted_judges.verdicts, the matched judge probe; router_planted.measure
  routing   router_planted.measure: forwarded, caught when asked, and the two jointly
  cost      judge_cost.rates on the D-110 roster (8 open-weight, 3 closed), and
            router_residue.measure at 11 models
  text      numbers the summary's prose states that no other script prints: the
            annotation app's recorded time, the slice's answer types against the pool's,
            steps per branch, the planted slips the shipped digit rule gives up, Wilson
            intervals on the judge rates, and the PRM's GPU time at full scale

It needs what every analysis here needs - traces/, scores/, experts_filled_labels/ and
analysis/out/planted/, all local and gitignored - and the .venv, for planted_judges. The
JSON it writes holds aggregate numbers only: no label, no trace, nothing per annotator.
"""
import glob
import io
import json
import math
import os as _os
import sys
from collections import Counter
from contextlib import redirect_stdout

_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS, _PILOT]

import cluster_bootstrap as CB  # noqa: E402
import digit_rule as DR  # noqa: E402
import e2_prm  # noqa: E402
import judge_cost as JC  # noqa: E402
import router_residue as RR  # noqa: E402
import score_against_labels as S  # noqa: E402
import x1_analysis as X  # noqa: E402

with redirect_stdout(io.StringIO()):
    import planted_judges as PJ  # noqa: E402
    import router_planted as RP  # noqa: E402

LABELS = _os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels')
OUT = _os.path.join(_PILOT, 'figures', 'summary_numbers.json')
B = 2000
OPEN, CLOSED = 8, 3                      # the roster, D-110

TRACE_ROWS = [('E0 (published)', 'E0'), ('E0 + 3rd judge', 'E0-3J'), ('E1 (judges swapped)', 'E1'),
              ('E2 (PRM, fraction)', 'E2 (72B frac)'), ('E2 (PRM, minimum)', 'E2 (72B min)'),
              ('E3 (milestones)', 'E3'), ('E4 (+ arithmetic)', 'E4'),
              ('E5 (milestones + judge)', 'E5'), ('the experts\' answer verdict', CB.BASELINE)]


def trace():
    truth, keyfile, tpl, rows = CB.load(LABELS)
    codes = sorted(truth)
    groups = {'template': {c: tpl[keyfile[c][1]] for c in codes}}
    score, y = CB.trace_scores(codes, truth, rows, keyfile)
    out = []
    for label, name in TRACE_ROWS:
        ok = [c for c in codes if score[name][c] is not None]
        fn = lambda ks, nm=name: X.auroc([(score[nm][c], y[c]) for c in ks])
        pt, _se, lo, hi = CB.spread(fn, ok, groups['template'], B)
        out.append({'label': label, 'head': name, 'auroc': pt, 'lo': lo, 'hi': hi})
    return {'rows': out, 'margins': CB.margins_over(codes, y, score, groups, B)}


def _prf(c):
    p, r, f = S.prf(c['tp'], c['fp'], c['fn'])
    return {'precision': p, 'recall': r, 'f1': f, 'tp': c['tp'], 'fp': c['fp'], 'fn': c['fn']}


def steps():
    """Inside correct-answer traces, per step, against the experts' labels."""
    truth, keyfile, texts = DR.load(LABELS)
    right = [c for c in truth if truth[c]['final_answer'] == 'correct']
    out = {}
    for rule in ('tol1', 'e4'):
        c = Counter()
        for code in right:
            text = texts.get(keyfile[code])
            split = e2_prm.steps_of(text) if text else []
            if len(split) != len(truth[code]['steps']):
                continue
            for s, lab in zip(split, truth[code]['steps']):
                n = len(DR.flagged(s, rule))
                if lab['label'] == 'incorrect':
                    c['tp' if n else 'fn'] += 1
                elif lab['label'] in ('correct', 'alternative_correct'):
                    c['fp' if n else 'tn'] += 1
        out['digit_' + rule] = _prf(c)
    e2 = S.current('e2')
    c = Counter()
    for code in right:
        row = e2.get(keyfile[code])
        probs = row['meta']['prm']['qwen72']['step_probs'] if row else None
        if probs is None or len(probs) != len(truth[code]['steps']):
            continue
        for lab, v in zip(truth[code]['steps'], probs):
            flagged = v < e2_prm.VERDICT_AT
            if lab['label'] == 'incorrect':
                c['tp' if flagged else 'fn'] += 1
            elif lab['label'] in ('correct', 'alternative_correct'):
                c['fp' if flagged else 'tn'] += 1
    out['prm_qwen72'] = _prf(c)
    return out


def planted(rp):
    probe, got = PJ.verdicts()
    fam = {p['set_id']: p['family'] for p in probe}
    judges = {}
    for j, _model, _price in PJ.JUDGES:
        judges[j] = {}
        for f in ('conceptual', 'arithmetic'):
            sets = [s for (jj, s) in got if jj == j and fam.get(s) == f and len(got[(jj, s)]) == 2]
            judges[j][f] = {
                'n': len(sets),
                'caught': sum(PJ._wrong(got[(j, s)]['planted']['verdict'])
                              and not PJ._wrong(got[(j, s)]['original']['verdict']) for s in sets),
                'flags_original': sum(PJ._wrong(got[(j, s)]['original']['verdict']) for s in sets)}
    return {'judges': judges,
            'digit_rule': {f: {'n': rp[f]['n'], 'caught': rp[f]['digit_rule']}
                           for f in ('conceptual', 'arithmetic')},
            'digit_false_alarm_on_judged_steps': rp['digit_false_alarm_on_judged_steps']}


def cost():
    rates = JC.rates()
    residue = RR.measure(LABELS, OPEN + CLOSED)
    return {'models': OPEN + CLOSED, 'open': OPEN, 'closed': CLOSED,
            'traces': (OPEN + CLOSED) * JC.ITEMS,
            'e0': JC.roster_cost(rates['E0'], OPEN, CLOSED),
            'e0_range': [JC.ITEMS * (OPEN + CLOSED) * rates['E0']['W'],
                         JC.ITEMS * (OPEN + CLOSED) * rates['E0']['F']],
            'e5': JC.roster_cost(rates['E5'], OPEN, CLOSED),
            'router_batched': residue['cost_batched'],
            'router_per_step': residue['cost_per_step']}


def wilson(k, n, z=1.959964):
    """95% Wilson interval for k of n."""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [centre - half, centre + half]


# RESULTS_E2: the 72B scored all 360 traces in 6 minutes on 2 GPUs.
PRM_GPU_MINUTES_PER_TRACE = 6 * 2 / 360


def text(planted_numbers, cost_numbers):
    import freeze as F                       # the slice's own answer-type inventory
    # the annotation app's timer: seconds a trace was open before it was submitted
    latest = {}
    for f in glob.glob(_os.path.join(LABELS, '*.jsonl')):
        for line in open(f, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                k = (r['annotator'], r['code'])
                if k not in latest or (r.get('ts') or '') >= (latest[k].get('ts') or ''):
                    latest[k] = r
    seconds = [r['seconds'] for r in latest.values() if isinstance(r.get('seconds'), (int, float))]
    truth, _ = S.build_truth(S.read_labels(LABELS), S.read_consensus(LABELS))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    branch = {}
    for line in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8'):
        it = json.loads(line)
        branch[it['item_id']] = it['branch']
    per_branch = Counter()
    for c, t in truth.items():
        per_branch[branch[keyfile[c][1]]] += len(t['steps'])
    freeze = json.load(open(_os.path.join(_PILOT, 'slice', 'FREEZE.json'), encoding='utf-8'))
    J = planted_numbers['judges']
    return {
        'annotation': {'submissions': len(latest), 'with_seconds': len(seconds),
                       'hours': sum(seconds) / 3600, 'median_seconds': sorted(seconds)[len(seconds) // 2]},
        'answer_types': {'slice_items': freeze['by_answer_type'],
                         'pool_templates': dict(Counter(F.answer_types().values()))},
        'steps_by_branch': dict(per_branch),
        'digit_rule_lost_relative_errors': RP.lost_to_shipped_rule(),
        'wilson': {j: {f: wilson(J[j][f]['caught'], J[j][f]['n']) for f in J[j]} for j in J},
        'prm_gpu_hours_full_run': cost_numbers['traces'] * PRM_GPU_MINUTES_PER_TRACE / 60,
    }


def main():
    with redirect_stdout(io.StringIO()):
        rp = RP.measure()
    numbers = {'trace': trace(), 'steps': steps(), 'planted': planted(rp),
               'routing': {f: {'n': rp[f]['n'], 'forwarded': rp[f]['forwarded'],
                               'e0': rp[f]['e0'], 'router': rp[f]['router']}
                           for f in ('conceptual', 'arithmetic')},
               'stack_arithmetic': rp['arithmetic']['stack'],
               'cost': cost()}
    numbers['text'] = text(numbers['planted'], numbers['cost'])
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(numbers, fh, indent=1, sort_keys=True)
        fh.write('\n')
    print('wrote', _os.path.relpath(OUT, _PILOT))
    for r in numbers['trace']['rows']:
        print('  %-28s %.3f  (%.3f, %.3f)' % (r['label'], r['auroc'], r['lo'], r['hi']))
    m = numbers['trace']['margins']['E5']
    print('  experts minus E5  %+.3f  (%+.3f, %+.3f)' % (m['margin'], m['lo'], m['hi']))
    for k, v in numbers['steps'].items():
        print('  %-12s precision %.3f  recall %.3f' % (k, v['precision'], v['recall']))
    c = numbers['cost']
    print('  cost  E0 %.2f  E5 %.2f  router batched %.2f  per step %.2f'
          % (c['e0'], c['e5'], c['router_batched'], c['router_per_step']))
    t = numbers['text']
    print('  annotation  %d submissions, %.1f hours recorded, median %.0f s'
          % (t['annotation']['submissions'], t['annotation']['hours'], t['annotation']['median_seconds']))
    print('  answer types  slice items %s  pool templates %s' % (t['answer_types']['slice_items'],
                                                               t['answer_types']['pool_templates']))
    print('  steps by branch  %s' % t['steps_by_branch'])
    print('  digit rule, planted slips the shipped rule gives up: %s'
          % ', '.join('%.1e' % x for x in t['digit_rule_lost_relative_errors']))
    print('  Wilson 95%% MiMo conceptual (%.3f, %.3f)' % tuple(t['wilson']['mimo-v2.5-pro']['conceptual']))
    print('  PRM GPU-hours at full scale %.1f' % t['prm_gpu_hours_full_run'])


if __name__ == '__main__':
    main()
