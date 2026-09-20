"""Stage 3: turn expert labels into ground truth, then score every evaluator on it.

    python evaluator_pilot_17092026/annotation/score_against_labels.py            # real labels
    python evaluator_pilot_17092026/annotation/score_against_labels.py --simulate # pipeline test

This is the referee the pilot has been missing. Until it runs on real labels, every
comparison between E0-E5 is evaluators against each other (RESULTS_E5, RESULTS_E2).

GROUND TRUTH. Three experts of the trace's own branch label each trace (`for_truth`).
Per step, the majority label wins; with no majority the step is DISPUTED and is
excluded from scoring and reported. Labels given outside an annotator's branch are
calibration only and never enter the truth.

AGREEMENT. Fleiss' kappa over the three raters on step labels (4 categories),
milestone status (3) and the final-answer verdict (4), per branch and pooled. The
shared calibration set is scored separately, across branches: it measures whether the
branches apply the same standard, not the ground truth.

SCORING, at the level each evaluator actually works at:
  step        E2's three PRMs: a step is flagged when its reward < 0.5. Precision,
              recall, F1 and AUROC against "incorrect".
  milestone   E3, E4, E5: milestone reached vs the experts' "obtained".
  trace       every evaluator's headline score vs the experts' overall verdict
              (AUROC for separating sound from unsound traces), and each evaluator's
              final-answer check against the experts' final-answer verdict.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import random
import statistics as st
import sys
from collections import Counter, defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
sys.path[:0] = [os.path.join(_PILOT, 'evaluators')]

LABELS = os.path.join(_HERE, 'labels')
TASKS = os.path.join(_HERE, 'tasks')
STEP_CATS = ('correct', 'alternative_correct', 'incorrect', 'not_a_claim')
MS_CATS = ('reached', 'not_reached', 'not_needed')
FINAL_CATS = ('correct', 'incorrect', 'partial', 'not_stated')
PRMS = ('qwen72', 'versa', 'qwen7')


def current(d):
    rows = [json.loads(l) for f in glob.glob(os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    if not rows:
        return {}
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def read_labels(d=LABELS):
    """Last row per (annotator, code) wins."""
    out = {}
    for f in sorted(glob.glob(os.path.join(d, '*.jsonl'))):
        for line in open(f, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                out[(r['annotator'], r['code'])] = r
    return list(out.values())


def fleiss(rows_of_counts, cats):
    """Fleiss' kappa. rows_of_counts: one Counter per item, over `cats`."""
    rows = [c for c in rows_of_counts if sum(c.values()) > 1]
    if not rows:
        return float('nan')
    n = sum(rows[0].values())
    if any(sum(c.values()) != n for c in rows):          # unequal raters: drop to common n
        n = min(sum(c.values()) for c in rows)
        rows = [c for c in rows if sum(c.values()) == n]
    N = len(rows)
    p_j = [sum(c[k] for c in rows) / (N * n) for k in cats]
    P_i = [(sum(c[k] ** 2 for k in cats) - n) / (n * (n - 1)) for c in rows]
    P_bar, P_e = sum(P_i) / N, sum(p * p for p in p_j)
    return (P_bar - P_e) / (1 - P_e) if P_e < 1 else float('nan')


def auroc(pos, neg):
    if not pos or not neg:
        return float('nan')
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else float('nan')
    r = tp / (tp + fn) if tp + fn else float('nan')
    f = 2 * p * r / (p + r) if p == p and r == r and p + r else float('nan')
    return p, r, f


def majority(values):
    """The value held by more than half, else None (disputed)."""
    c = Counter(values)
    top, n = c.most_common(1)[0]
    return top if n * 2 > len(values) else None


def build_truth(labels):
    """Majority ground truth per trace, from that trace's own-branch experts."""
    by_code = defaultdict(list)
    for r in labels:
        if r.get('for_truth'):
            by_code[r['code']].append(r)
    truth, disputed = {}, Counter()
    for code, rs in by_code.items():
        steps = []
        for i in range(rs[0]['n_steps']):
            lab = majority([r['steps'][i]['label'] for r in rs])
            disputed['step'] += lab is None
            errs = [r['steps'][i].get('error_type') for r in rs if r['steps'][i]['label'] == 'incorrect']
            steps.append({'label': lab, 'error_type': majority(errs) if errs else None,
                          'n_raters': len(rs)})
        ms = []
        for j, m in enumerate(rs[0]['milestones']):
            statm = majority([r['milestones'][j]['status'] for r in rs])
            disputed['milestone'] += statm is None
            ms.append({'id': m['id'], 'value': m['value'], 'status': statm})
        fin = majority([r['final_answer'] for r in rs])
        sound = majority([r['reasoning_sound'] for r in rs])
        disputed['final'] += fin is None
        disputed['sound'] += sound is None
        truth[code] = {'code': code, 'branch': rs[0]['branch'], 'n_raters': len(rs),
                       'steps': steps, 'milestones': ms, 'final_answer': fin,
                       'reasoning_sound': sound,
                       'seconds': [r['seconds'] for r in rs]}
    return truth, disputed


def report_agreement(labels):
    print('\nAGREEMENT (Fleiss kappa over the experts of each trace\'s own branch)')
    by_code = defaultdict(list)
    for r in labels:
        if r.get('for_truth'):
            by_code[r['code']].append(r)
    def kappas(codes):
        s = [Counter(r['steps'][i]['label'] for r in by_code[c])
             for c in codes for i in range(by_code[c][0]['n_steps'])]
        m = [Counter(r['milestones'][j]['status'] for r in by_code[c])
             for c in codes for j in range(len(by_code[c][0]['milestones']))]
        f = [Counter(r['final_answer'] for r in by_code[c]) for c in codes]
        return fleiss(s, STEP_CATS), fleiss(m, MS_CATS), fleiss(f, FINAL_CATS), len(s)
    print('  %-24s %8s %10s %8s %8s' % ('scope', 'steps', 'milestones', 'final', 'n steps'))
    for br in sorted({r['branch'] for r in labels if r.get('for_truth')}):
        codes = [c for c in by_code if by_code[c][0]['branch'] == br]
        k = kappas(codes)
        print('  %-24s %8.3f %10.3f %8.3f %8d' % (br.replace('_', ' '), *k))
    k = kappas(list(by_code))
    print('  %-24s %8.3f %10.3f %8.3f %8d' % ('POOLED', *k))

    calib = defaultdict(list)
    for r in labels:
        if r.get('calibration'):
            calib[r['code']].append(r)
    if calib:
        s = [Counter(r['steps'][i]['label'] for r in calib[c])
             for c in calib for i in range(calib[c][0]['n_steps'])]
        f = [Counter(r['final_answer'] for r in calib[c]) for c in calib]
        print('  %-24s %8.3f %10s %8.3f %8d   (all annotators, across branches)'
              % ('CALIBRATION SET', fleiss(s, STEP_CATS), '-', fleiss(f, FINAL_CATS), len(s)))


def report_steps(truth, keyfile):
    print('\nSTEP LEVEL - E2\'s process reward models against the experts')
    e2 = current('e2')
    if not e2:
        print('  no E2 rows found'); return
    print('  %-8s %6s %6s %9s %8s %8s %8s' % ('PRM', 'TP', 'FP', 'FN', 'prec', 'rec', 'F1'))
    for p in PRMS:
        tp = fp = fn = 0
        pos, neg = [], []
        for code, t in truth.items():
            key = keyfile.get(code)
            row = e2.get(key)
            if not row:
                continue
            probs = row['meta']['prm'][p]['step_probs']
            if probs is None or len(probs) != len(t['steps']):
                continue
            for s, pr in zip(t['steps'], probs):
                if s['label'] in (None, 'not_a_claim'):
                    continue
                bad = s['label'] == 'incorrect'
                flagged = pr < 0.5
                tp += bad and flagged
                fp += (not bad) and flagged
                fn += bad and not flagged
                (pos if bad else neg).append(pr)
        pr_, rc, f1 = prf(tp, fp, fn)
        print('  %-8s %6d %6d %9d %8.3f %8.3f %8.3f   AUROC %.3f'
              % (p, tp, fp, fn, pr_, rc, f1, 1 - auroc(pos, neg)))


def report_milestones(truth, keyfile):
    print('\nMILESTONE LEVEL - E3, E4 and E5 against the experts')
    tables = {'e3': current('e3'), 'e4': current('e4'), 'e5': current('e5')}
    print('  %-6s %6s %6s %6s %8s %8s %8s' % ('eval', 'TP', 'FP', 'FN', 'prec', 'rec', 'F1'))
    for name, rows in tables.items():
        if not rows:
            continue
        tp = fp = fn = 0
        for code, t in truth.items():
            row = rows.get(keyfile.get(code))
            if not row:
                continue
            got = {m['id']: m for m in row['meta']['milestones']}
            for m in t['milestones']:
                if m['status'] is None or m['id'] not in got:
                    continue
                said = bool(got[m['id']].get('reached'))
                if name == 'e5':
                    # E5 rows carry source 'e3' (found deterministically) or the judge's
                    # verdict string; e5_strict credits E3's matches plus REACHED.
                    said = got[m['id']].get('source') in ('e3', 'REACHED')
                truly = m['status'] == 'reached'
                tp += truly and said
                fp += (not truly) and said
                fn += truly and not said
        p, r, f = prf(tp, fp, fn)
        print('  %-6s %6d %6d %6d %8.3f %8.3f %8.3f' % (name, tp, fp, fn, p, r, f))


def report_traces(truth, keyfile):
    print('\nTRACE LEVEL - does the evaluator\'s score separate sound from unsound traces?')
    heads = [('e0', 'recovered_f1'), ('e0_3j', 'recovered_f1'), ('e1', 'recovered_f1'),
             ('e3', 'milestone_coverage'), ('e4', 'e4_coverage'), ('e5', 'e5_strict'),
             ('e2', 'e2'), ('e2', 'qwen72_min')]
    print('  %-6s %-18s %8s %8s %8s' % ('eval', 'score', 'AUROC', 'n sound', 'n unsound'))
    for name, field in heads:
        rows = current(name)
        if not rows:
            continue
        pos, neg = [], []
        for code, t in truth.items():
            row = rows.get(keyfile.get(code))
            if not row or t['reasoning_sound'] is None:
                continue
            v = row['scores'].get(field)
            if v is None:
                continue
            (pos if t['reasoning_sound'] == 'yes' else neg).append(v)
        print('  %-6s %-18s %8.3f %8d %8d' % (name, field, auroc(pos, neg), len(pos), len(neg)))

    print('\n  Each evaluator\'s final-answer check against the experts\' verdict')
    print('  %-6s %8s %8s %8s' % ('eval', 'agree', 'n', 'note'))
    for name in ('e0', 'e1'):
        rows = current(name)
        if not rows:
            continue
        ok = n = 0
        for code, t in truth.items():
            row = rows.get(keyfile.get(code))
            if not row or t['final_answer'] is None:
                continue
            n += 1
            ok += int(row['scores']['final_answer_acc'] == 1) == (t['final_answer'] == 'correct')
        print('  %-6s %8.3f %8d   %s' % (name, ok / n if n else float('nan'), n,
                                         'E0-F1/F2 predict disagreement here'))


def simulate():
    """Write synthetic labels so the pipeline can be tested before experts start.

    Labels are drawn from E2's primary PRM plus noise - NOT ground truth, and never
    written to annotation/labels/. Purpose: prove this script runs end to end.
    """
    pool = json.load(open(os.path.join(TASKS, 'pool.json'), encoding='utf-8'))
    assignment = json.load(open(os.path.join(TASKS, 'assignment.json'), encoding='utf-8'))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(os.path.join(TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    e2, e3 = current('e2'), current('e3')
    out = os.path.join(_HERE, 'labels_simulated')
    os.makedirs(out, exist_ok=True)
    for f in glob.glob(os.path.join(out, '*.jsonl')):
        os.remove(f)
    rnd = random.Random(7)
    for who, me in assignment.items():
        rows = []
        for code in me['order']:
            t = pool[code]
            key = keyfile[code]
            probs = (e2.get(key) or {}).get('meta', {}).get('prm', {}).get('qwen72', {}).get('step_probs')
            steps = []
            for i in range(len(t['steps'])):
                pr = probs[i] if probs and i < len(probs) else 0.9
                lab = 'incorrect' if (pr < 0.5) != (rnd.random() < 0.15) else 'correct'
                steps.append({'index': i, 'label': lab,
                              'error_type': 'calculation' if lab == 'incorrect' else None, 'note': ''})
            ms3 = {m['id']: m for m in (e3.get(key) or {}).get('meta', {}).get('milestones', [])}
            ms = [{'id': m['id'], 'value': m['value'],
                   'status': 'reached' if ms3.get(m['id'], {}).get('reached', rnd.random() < 0.7)
                             else 'not_reached'} for m in t['milestones']]
            rows.append({
                'annotator': who, 'branch_of_annotator': me['branch'], 'code': code,
                'branch': t['branch'], 'level': t['level'],
                'for_truth': code in me['for_truth'], 'calibration': code in me['calibration'],
                'item_sha256': t['item_sha256'], 'trace_sha256': t['trace_sha256'],
                'n_steps': len(t['steps']), 'steps': steps, 'milestones': ms,
                'final_answer': 'correct' if all(s['label'] != 'incorrect' for s in steps) else 'incorrect',
                'reasoning_sound': 'yes' if sum(s['label'] == 'incorrect' for s in steps) == 0 else 'no',
                'confidence': 'medium', 'comment': '', 'flagged': False,
                'seconds': rnd.uniform(300, 800), 'ts': '2026-01-01T00:00:00Z'})
        with open(os.path.join(out, who + '.jsonl'), 'w', encoding='utf-8', newline='\n') as fh:
            for r in rows:
                fh.write(json.dumps(r) + '\n')
    print('SIMULATED labels in %s - synthetic, not ground truth' % out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--simulate', action='store_true',
                    help='generate synthetic labels and run on them, to test the pipeline')
    a = ap.parse_args()

    src = simulate() if a.simulate else LABELS
    labels = read_labels(src)
    if not labels:
        raise SystemExit('no labels in %s yet' % src)
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(os.path.join(TASKS, 'keyfile.jsonl'), encoding='utf-8'))}

    done = {(r['annotator'], r['code']) for r in labels}
    print('labels: %d submissions from %d annotators, %d distinct traces'
          % (len(done), len({a_ for a_, _ in done}), len({c for _, c in done})))
    secs = [r['seconds'] for r in labels if r.get('seconds')]
    if secs:
        print('time per trace: median %.1f min, total %.1f hours'
              % (st.median(secs) / 60, sum(secs) / 3600))

    truth, disputed = build_truth(labels)
    full = sum(1 for t in truth.values() if t['n_raters'] >= 3)
    print('ground truth: %d traces (%d with all three experts); disputed: %s'
          % (len(truth), full, dict(disputed)))
    if a.simulate:
        print('*** SYNTHETIC LABELS - every number below is a pipeline test, not a result ***')

    report_agreement(labels)
    report_steps(truth, keyfile)
    report_milestones(truth, keyfile)
    report_traces(truth, keyfile)


if __name__ == '__main__':
    main()
