"""Would a verify-first step router show a judge the planted defect, and what would it catch?

    python evaluator_pilot_17092026/analysis/router_planted.py

RESULTS_X1 recommendation 5 argues for a router: verify what a checker can verify, and send
only the residue to a judge. `router_residue.py` measures what that costs. What it buys was
first stated from two rates multiplied together, and one of them - that the router always
shows the judge the corrupted step - was assumed, not measured. This measures both, per
planted defect, from what the pilot already holds:

  forwarded   `router_residue.verifiable`, the router's own rule, applied to the step the
              defect sits in. A step holding a claim the checker can recompute is settled
              by the checker and never reaches a judge; everything else is residue.
  caught      the matched verdicts `planted_judges.py` stored: a judge catches a defect
              when it flags the planted step and does not flag the same step untouched.
  end to end  forwarded AND caught, counted defect by defect - not a product of rates.

The published framework's routing is put on the same basis, from the rows
`planted_routing.py` recorded: the planted step shown to the Tribunal AND caught by either
of E0's two judges (with two votes a split falls to `min()`, so either flag is decisive).

Two more comparisons the pilot summary makes, measured here so neither is a transcription:
  - how often the digit rule, as E4 ships it, flags the very steps the judges were shown
    untouched - a per-step false-alarm rate, the unit the judges' 0.000 is in;
  - what the full stack catches of the planted ARITHMETIC defects: the digit rule on the
    steps the router settles, and the router's judge on the residue.

No model is called: every verdict is read from `scores/planted_judges/replies.jsonl`.
The planted set is `analysis/out/planted/planted.jsonl` (rebuilt by `planted.py`).
"""
import io
import os as _os
import sys
from contextlib import redirect_stdout

_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]

import json  # noqa: E402

import e2_prm  # noqa: E402
from router_residue import verifiable  # noqa: E402

with redirect_stdout(io.StringIO()):
    import planted as P  # noqa: E402
    import planted_judges as PJ  # noqa: E402

ROUTING = _os.path.join(_PILOT, 'scores', 'planted_routing', 'routing.jsonl')
E0_JUDGES = ('gpt-5', 'opus-4.5')
ROUTER_JUDGES = ('mimo-v2.5-pro', 'gpt-5')     # the stack's judge first


def _rows(path):
    return [json.loads(l) for l in open(path, encoding='utf-8') if l.strip()]


def originals():
    out = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for r in _rows(_os.path.join(_PILOT, 'traces', f)):
                if r.get('ok'):
                    out[(r['model_key'], r['item_id'])] = r['text']
    return out


def measure():
    """Everything this script reports, as one dict, so the summary figures can read it."""
    plants = [r for r in _rows(PJ.PLANTED) if r['family'] in ('arithmetic', 'conceptual')]
    orig = originals()
    _probe, got = PJ.verdicts()

    def caught(judge, sid):
        a = got.get((judge, sid))
        if not a or len(a) != 2:
            return None                  # the judge never returned on one of the two arms
        return PJ._wrong(a['planted']['verdict']) and not PJ._wrong(a['original']['verdict'])

    shown = {}
    for r in _rows(ROUTING):
        if r['arm'] == 'planted':
            shown[r['set_id']] = bool(r['shown'])

    rec = {}
    for r in plants:
        steps = e2_prm.steps_of(r['text'])
        i = r['step_index']
        if i is None or i >= len(steps) or r['after'] not in steps[i]:
            raise SystemExit('%s: the planted text is not in step %s - step spaces disagree'
                             % (r['set_id'], i))
        before = P.digit_flags(orig[(r['model_key'], r['item_id'])], 'e4')
        after = P.digit_flags(r['text'], 'e4')
        rec[r['set_id']] = {
            'family': r['family'], 'kind': r.get('defect_kind') or r.get('basis'),
            'forwarded': not verifiable(steps[i]),
            'e0_shown': shown.get(r['set_id']),
            'digit_new_flag': i < len(after) and i < len(before) and len(after[i]) > len(before[i]),
            'digit_flags_original_step': i < len(before) and bool(before[i]),
            'caught': {j: caught(j, r['set_id']) for j in set(E0_JUDGES) | set(ROUTER_JUDGES)},
        }

    out = {}
    for fam in ('conceptual', 'arithmetic'):
        ks = [k for k, v in rec.items() if v['family'] == fam]
        f = {'n': len(ks),
             'forwarded': sum(rec[k]['forwarded'] for k in ks),
             'forwarded_by_kind': {}}
        for kind in sorted({rec[k]['kind'] for k in ks}):
            kk = [k for k in ks if rec[k]['kind'] == kind]
            f['forwarded_by_kind'][kind] = [sum(rec[k]['forwarded'] for k in kk), len(kk)]

        # E0's routing, measured jointly, on the defects both E0 judges returned on
        e0 = [k for k in ks if rec[k]['e0_shown'] is not None
              and all(rec[k]['caught'][j] is not None for j in E0_JUDGES)]
        either = lambda k: any(rec[k]['caught'][j] for j in E0_JUDGES)
        f['e0'] = {'n': len(e0),
                   'shown': sum(rec[k]['e0_shown'] for k in e0),
                   'caught_when_asked': sum(either(k) for k in e0),
                   'end_to_end': sum(rec[k]['e0_shown'] and either(k) for k in e0)}

        f['router'] = {}
        for j in ROUTER_JUDGES:
            ok = [k for k in ks if rec[k]['caught'][j] is not None]
            f['router'][j] = {'n': len(ok),
                              'forwarded': sum(rec[k]['forwarded'] for k in ok),
                              'caught_when_asked': sum(rec[k]['caught'][j] for k in ok),
                              'end_to_end': sum(rec[k]['forwarded'] and rec[k]['caught'][j] for k in ok)}
        f['digit_rule'] = sum(rec[k]['digit_new_flag'] for k in ks)
        if fam == 'arithmetic':
            j = ROUTER_JUDGES[0]
            ok = [k for k in ks if rec[k]['caught'][j] is not None]
            f['stack'] = {'n': len(ok), 'judge': j,
                          'caught': sum(rec[k]['digit_new_flag']
                                        or (rec[k]['forwarded'] and rec[k]['caught'][j]) for k in ok),
                          'missed_by_both': sum(not rec[k]['digit_new_flag'] and
                                                not (rec[k]['forwarded'] and rec[k]['caught'][j])
                                                for k in ok)}
        out[fam] = f
    out['digit_false_alarm_on_judged_steps'] = [sum(v['digit_flags_original_step'] for v in rec.values()),
                                                len(rec)]
    return out


def lost_to_shipped_rule():
    """Relative size of the planted arithmetic slips, inside a parseable claim, that the bare
    digit rule catches and the rule E4 ships does not - what its two gold-validation
    corrections cost on defects that are certainly present."""
    orig = originals()
    out = []
    for r in _rows(PJ.PLANTED):
        if r['family'] != 'arithmetic' or r.get('basis') != 'claim':
            continue
        base, i = orig[(r['model_key'], r['item_id'])], r['step_index']
        new_flag = {}
        for rule in ('digit', 'e4'):
            a, b = P.digit_flags(r['text'], rule), P.digit_flags(base, rule)
            new_flag[rule] = len(a[i]) > len(b[i])
        if new_flag['digit'] and not new_flag['e4']:
            out.append(r['relative_error'])
    return sorted(out)


def main():
    m = measure()
    rate = lambda a, b: '%3d of %-3d = %.3f' % (a, b, a / b if b else float('nan'))
    print('A VERIFY-FIRST ROUTER ON THE PLANTED DEFECTS  (no call made; stored verdicts)\n')
    for fam in ('conceptual', 'arithmetic'):
        f = m[fam]
        print('%s, %d defects' % (fam.upper(), f['n']))
        print('  forwarded to a judge by the router      %s' % rate(f['forwarded'], f['n']))
        for kind, (a, b) in f['forwarded_by_kind'].items():
            print('      %-20s                 %s' % (kind, rate(a, b)))
        e = f['e0']
        print('  published routing, E0\'s two judges (n = %d defects both judges returned on)' % e['n'])
        print('    planted step shown                    %s' % rate(e['shown'], e['n']))
        print('    caught by either judge when asked     %s' % rate(e['caught_when_asked'], e['n']))
        print('    end to end, shown AND caught          %s' % rate(e['end_to_end'], e['n']))
        for j, r in f['router'].items():
            print('  router, judge %s (n = %d defects it returned on)' % (j, r['n']))
            print('    forwarded                             %s' % rate(r['forwarded'], r['n']))
            print('    caught when asked                     %s' % rate(r['caught_when_asked'], r['n']))
            print('    end to end, forwarded AND caught      %s' % rate(r['end_to_end'], r['n']))
        print('  digit rule as E4 ships it, new flag on the planted step  %s' % rate(f['digit_rule'], f['n']))
        if fam == 'arithmetic':
            s = f['stack']
            print('  the stack: digit rule, or forwarded and caught by %s   %s'
                  % (s['judge'], rate(s['caught'], s['n'])))
            print('    missed by both                        %s' % rate(s['missed_by_both'], s['n']))
        print()
    a, b = m['digit_false_alarm_on_judged_steps']
    print('FALSE ALARMS on the steps the judges were shown untouched')
    print('  digit rule as E4 ships it flags the untouched step   %s' % rate(a, b))
    print('  (the judges flagged none of them; see planted_judges.py report)')
    lost = lost_to_shipped_rule()
    print('\nWHAT THE SHIPPED DIGIT RULE GIVES UP, inside a parseable claim: %d planted slips' % len(lost))
    print('  relative size: %s' % ', '.join('%.1e' % x for x in lost))


if __name__ == '__main__':
    main()
