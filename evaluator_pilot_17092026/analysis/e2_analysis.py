"""E2: do the PRMs discriminate, do they agree, and what do they say about the models?

    python evaluator_pilot_17092026/analysis/e2_analysis.py

Reads the harness rows (scores/e2, current config) - no GPU, no calls.

1. VALIDATION on steps whose status is known without any model (D-087): the judge
   probe's labelled steps (analysis/judge_probe.py - SLIP = shown arithmetic wrong,
   confirmed by reading; CLEAN = a step of a right-answer trace whose checkable claims
   all hold; RELABEL applied). A PRM that discriminates gives SLIPs lower rewards:
   reported per PRM as AUROC (P(reward_clean > reward_slip)), slips flagged (< 0.5)
   and clean steps passed (>= 0.5). The probe indexes E0's raw extract_steps list;
   E2 drops empty steps, so each index is mapped, and the step text is checked.
2. PER MODEL: mean frac_ok, mean min reward, share of traces with no flagged step,
   per PRM, for the gold five and gemma.
3. AGREEMENT between PRMs: Cohen's kappa on step verdicts, Spearman on trace frac_ok.
4. AGAINST OTHER EVALUATORS: per-trace Pearson with E0's final-answer correctness,
   E3 coverage and E5-strict. As in RESULTS_E5, a correlation with E0's own
   correctness check does not rank evaluators; only expert labels can.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import glob
import json
import statistics as st
import sys
from collections import defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'evaluators'), _os.path.join(_REPO, 'evaluation'), _ANALYSIS]
import e2_prm  # noqa: E402

PRMS = e2_prm.PRMS
GOLD = ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b')
ORDER = GOLD + ('gemma-4-31b',)


def current(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    if not rows:
        return {}
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def mean(xs):
    xs = [x for x in xs if x is not None]
    return st.mean(xs) if xs else float('nan')


def auroc(pos, neg):
    """P(a random positive scores above a random negative); ties count half."""
    if not pos or not neg:
        return float('nan')
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def pearson(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return float('nan')
    a, b = zip(*pairs)
    try:
        return st.correlation(a, b)
    except st.StatisticsError:
        return float('nan')


def spearman(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return float('nan')
    a, b = zip(*pairs)
    return st.correlation(a, b, method='ranked')


def kappa(a, b):
    n = len(a)
    po = sum(x == y for x, y in zip(a, b)) / n
    pa, pb = sum(a) / n, sum(b) / n
    pe = pa * pb + (1 - pa) * (1 - pb)
    return (po - pe) / (1 - pe) if pe < 1 else float('nan')


def validate(e2):
    from engineering_parser import extract_steps
    import judge_probe as jp
    probe = json.load(open(jp.PROBE, encoding='utf-8'))
    traces = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        k = _os.path.basename(f)[:-6]
        for l in open(f, encoding='utf-8'):
            r = json.loads(l)
            if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
                traces[(k, r['item_id'])] = r          # the trace E2 and the probe used
    rewards = {p: {'SLIP': [], 'CLEAN': []} for p in PRMS}
    rows_out, skipped = [], []
    for e in probe:
        label = jp.RELABEL.get((e['model'], e['item_id']), (e['label'],))[0]
        key = (e['model'], e['item_id'])
        if key not in e2:
            skipped.append('%s/%s: no E2 row' % key)
            continue
        raw = extract_steps(traces[key]['text'])[0]
        if raw[e['step']].strip() == '':
            skipped.append('%s/%s: labelled step is empty' % key)
            continue
        idx = sum(1 for s in raw[:e['step']] if s.strip())          # index after dropping empties
        steps = e2_prm.steps_of(traces[key]['text'])
        if steps[idx] != raw[e['step']].strip():
            skipped.append('%s/%s: step text mismatch after mapping' % key)
            continue
        got = {}
        for p in PRMS:
            probs = e2[key]['meta']['prm'][p]['step_probs']
            if probs is None:
                continue
            got[p] = probs[idx]
            rewards[p][label].append(probs[idx])
        rows_out.append((label, key, idx, got))
    return rewards, rows_out, skipped


def main():
    e2 = current('e2')
    if not e2:
        raise SystemExit('no E2 rows - run: python -m evaluator_pilot_17092026.run_evaluator e2')
    print('E2 rows (current config): %d' % len(e2))
    over = {p: sum(1 for r in e2.values() if r['meta']['prm'][p]['over_length']) for p in PRMS}
    print('over the context limit (not scored): %s' % over)

    print('\n1. VALIDATION on the judge probe\'s known-label steps')
    rewards, rows, skipped = validate(e2)
    for s in skipped:
        print('  skipped ' + s)
    print('  %-8s %5s %5s  %-10s %-10s  %-15s %-15s' % ('PRM', 'SLIP', 'CLEAN', 'mean SLIP', 'mean CLEAN',
                                                        'slips flagged', 'clean passed'))
    for p in PRMS:
        s, c = rewards[p]['SLIP'], rewards[p]['CLEAN']
        print('  %-8s %5d %5d  %-10.3f %-10.3f  %2d/%-2d (%3.0f%%)    %2d/%-2d (%3.0f%%)    AUROC %.3f'
              % (p, len(s), len(c), mean(s), mean(c),
                 sum(x < 0.5 for x in s), len(s), 100 * sum(x < 0.5 for x in s) / max(len(s), 1),
                 sum(x >= 0.5 for x in c), len(c), 100 * sum(x >= 0.5 for x in c) / max(len(c), 1),
                 auroc(c, s)))
    print('  per step:')
    for label, key, idx, got in sorted(rows):
        print('    %-5s %-14s %-32s step %2d  %s' % (label, key[0], key[1], idx,
                                                  '  '.join('%s %.3f' % (p, got[p]) for p in PRMS if p in got)))

    print('\n2. PER MODEL')
    by = defaultdict(list)
    for (mk, iid), r in e2.items():
        by[mk].append(r)
    for p in PRMS:
        print('  %s' % p)
        print('    %-16s %4s  %-8s %-8s %-8s' % ('model', 'n', 'frac_ok', 'min', 'clean'))
        for mk in ORDER:
            rs = by.get(mk, [])
            if rs:
                print('    %-16s %4d  %-8.3f %-8.3f %-8.3f' % (
                    mk, len(rs), mean(r['scores']['%s_frac_ok' % p] for r in rs),
                    mean(r['scores']['%s_min' % p] for r in rs),
                    mean(r['scores']['%s_clean' % p] for r in rs)))

    print('\n3. AGREEMENT between PRMs')
    keys = sorted(e2)
    for i, a in enumerate(PRMS):
        for b in PRMS[i + 1:]:
            va, vb = [], []
            for k in keys:
                pa, pb = e2[k]['meta']['prm'][a]['step_probs'], e2[k]['meta']['prm'][b]['step_probs']
                if pa is None or pb is None:
                    continue
                va += [x >= 0.5 for x in pa]
                vb += [x >= 0.5 for x in pb]
            rho = spearman([e2[k]['scores']['%s_frac_ok' % a] for k in keys],
                           [e2[k]['scores']['%s_frac_ok' % b] for k in keys])
            print('  %-7s vs %-7s  step verdicts: n=%d, agree %.3f, kappa %.3f   trace frac_ok Spearman %.3f'
                  % (a, b, len(va), sum(x == y for x, y in zip(va, vb)) / len(va), kappa(va, vb), rho))

    print('\n4. AGAINST OTHER EVALUATORS (gold five; per-trace Pearson)')
    e0, e3, e5 = current('e0'), current('e3'), current('e5')
    gk = [k for k in keys if k[0] in GOLD]
    refs = {'E0 final-answer correct': [e0[k]['scores']['final_answer_acc'] if k in e0 else None for k in gk],
            'E3 coverage': [e3[k]['scores']['milestone_coverage'] if k in e3 else None for k in gk],
            'E5 strict': [e5[k]['scores']['e5_strict'] if k in e5 else None for k in gk]}
    for p in PRMS:
        for f in ('frac_ok', 'min'):
            xs = [e2[k]['scores']['%s_%s' % (p, f)] for k in gk]
            print('  %-7s %-8s ' % (p, f) + '   '.join('%s %.3f' % (n, pearson(xs, ys)) for n, ys in refs.items()))


if __name__ == '__main__':
    main()
