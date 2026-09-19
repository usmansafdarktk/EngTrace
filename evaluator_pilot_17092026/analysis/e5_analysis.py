"""E5: is its judge trustworthy, and what does it add over E3?

    python evaluator_pilot_17092026/analysis/e5_analysis.py validate   # paid, ~$0.30
    python evaluator_pilot_17092026/analysis/e5_analysis.py report

VALIDATE - before any E5 number is trusted (D-087). E5's judge is asked about the
milestones E3 could not find, and a judge can say "yes, reached" to anything. So it
is tested, in E5's own prompt shape, on two sets whose answer is known without it:

  SENSITIVITY  milestones E3 DID find in the trace, with their true values.
               The right verdict is REACHED.
  SPECIFICITY  the same milestones with each value multiplied by 1.37 - a quantity
               the trace never states. The right verdict is MISSING.

A judge that accepts the perturbed values is rubber-stamping, and E5's residual
credit would be worthless. Fixed sample: 40 traces (seed 11), up to 3 milestones each.

REPORT - on current configs only: E5 strict / lenient against E3, E4 and E0, the
judged fraction, verdict counts, correlation with final-answer correctness, cost.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import concurrent.futures as cf
import glob
import json
import random
import statistics as st
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'evaluators'), _REPO]
import e1_panel as e1      # noqa: E402
import e3_milestones as e3  # noqa: E402
import e5_hybrid as e5     # noqa: E402

VALID = _os.path.join(_PILOT, 'scores', 'e5_validation.jsonl')
PERTURB = 1.37


def current(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    if not rows:
        return {}
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def validate():
    from dotenv import load_dotenv
    load_dotenv(_os.path.join(_REPO, '.env'))
    items = {json.loads(l)['item_id']: json.loads(l)
             for l in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8')}
    state = e3.setup()
    traces = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        k = _os.path.basename(f)[:-6]
        if k in ('qwen3.8-27b', 'gemma-4-31b'):
            continue
        for l in open(f, encoding='utf-8'):
            r = json.loads(l)
            if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
                traces[(k, r['item_id'])] = r
    pool = []
    for (mk, iid), tr in sorted(traces.items()):
        hits = e3.reach(state['milestones'][iid]['milestones'], tr['text'])
        got = [h for h in hits if h['reached']]
        if got:
            pool.append((mk, iid, got[:3]))
    random.seed(11)
    sample = random.sample(pool, 40)
    jobs = []
    for mk, iid, got in sample:
        jobs.append(('SENS', mk, iid, got, e5.prompt_for(items[iid], traces[(mk, iid)]['text'], got)))
        fake = [dict(h, value=h['value'] * PERTURB) for h in got]
        jobs.append(('SPEC', mk, iid, fake, e5.prompt_for(items[iid], traces[(mk, iid)]['text'], fake)))
    print('%d validation prompts (%d sensitivity, %d specificity)' % (
        len(jobs), sum(j[0] == 'SENS' for j in jobs), sum(j[0] == 'SPEC' for j in jobs)))
    with cf.ThreadPoolExecutor(max_workers=24) as ex, open(VALID, 'w', encoding='utf-8') as fh:
        for job, rep in zip(jobs, ex.map(lambda j: e1.fetch(e5.JUDGE, j[4]), jobs)):
            kind, mk, iid, ms_, _ = job
            fh.write(json.dumps({'kind': kind, 'model': mk, 'item_id': iid,
                                 'milestones': [m['id'] for m in ms_], 'reply': rep}, ensure_ascii=False) + '\n')
    print('written', VALID)


def report_validation():
    if not _os.path.exists(VALID):
        print('(no validation run yet)')
        return
    rows = [json.loads(l) for l in open(VALID, encoding='utf-8')]
    tally = {'SENS': Counter(), 'SPEC': Counter()}
    cost_in = cost_out = 0
    for r in rows:
        v = e5.parse(r['reply'].get('text')) if r['reply'].get('ok') else {}
        for m in r['milestones']:
            tally[r['kind']][v.get(m, 'UNPARSED')] += 1
        cost_in += r['reply'].get('billed_prompt_tokens') or 0
        cost_out += r['reply'].get('billed_completion_tokens') or 0
    print('VALIDITY OF THE E5 JUDGE (%s)' % e5.JUDGE)
    for kind, want, label in (('SENS', 'REACHED', 'true values - should be REACHED'),
                              ('SPEC', 'MISSING', 'values x1.37 - should be MISSING')):
        t = tally[kind]
        n = sum(t.values())
        print('  %-38s %3d/%-3d right (%.0f%%)   %s' % (label, t[want], n, 100 * t[want] / max(n, 1), dict(t)))
    print('  validation tokens: %d in, %d out' % (cost_in, cost_out))


def report():
    report_validation()
    e5r, e3r, e4r, e0r = current('e5'), current('e3'), current('e4'), current('e0')
    if not e5r:
        print('\n(no E5 rows yet)')
        return
    gold = [k for k in e5r if k in e0r]
    S = lambda r, k: r['scores'][k]

    def corr(a, b):
        ma, mb = st.mean(a), st.mean(b)
        num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
        den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** .5
        return num / den if den else float('nan')

    judged = [r for r in e5r.values() if r['calls']]
    cost = sum(c.get('cost_usd', 0) for r in judged for c in r['calls'])
    print('\nE5 over %d traces: %d sent to the judge, %d judge failures, $%.2f' % (
        len(e5r), len(judged), sum(1 for r in judged if r['meta'].get('judge_failed')), cost))
    req = sum(S(r, 'milestones_required') for r in e5r.values())
    jf = sum(S(r, 'judged_fraction') * S(r, 'milestones_required') for r in e5r.values())
    print('judged fraction: %.0f of %d milestones (%.1f%%) needed the judge; the rest were settled by E3'
          % (jf, req, 100 * jf / req))
    v = Counter()
    for r in e5r.values():
        v.update({'REACHED': S(r, 'judge_reached'), 'NOT_NEEDED': S(r, 'judge_not_needed'),
                  'MISSING': S(r, 'judge_missing'), 'UNJUDGED': S(r, 'unjudged')})
    print('judge verdicts on E3\'s misses: %s' % dict(v))

    print('\n%-16s %6s %6s %9s %10s %7s %7s' % ('model', 'E3', 'E4', 'E5 strict', 'E5 lenient', 'E0 F1', 'E0 FAC'))
    by = defaultdict(list)
    for k in gold:
        by[k[0]].append(k)
    for m in sorted(by, key=lambda m: -st.mean(S(e5r[k], 'e5_strict') for k in by[m])):
        ks = by[m]
        print('%-16s %6.3f %6.3f %9.3f %10.3f %7.3f %7.3f' % (
            m, st.mean(S(e3r[k], 'milestone_coverage') for k in ks), st.mean(S(e4r[k], 'e4_coverage') for k in ks),
            st.mean(S(e5r[k], 'e5_strict') for k in ks), st.mean(S(e5r[k], 'e5_lenient') for k in ks),
            st.mean(S(e0r[k], 'recovered_f1') for k in ks), st.mean(S(e0r[k], 'final_answer_acc') for k in ks)))
    fac = [S(e0r[k], 'final_answer_acc') for k in gold]
    print()
    for name, rows, key in (('E3', e3r, 'milestone_coverage'), ('E5 strict', e5r, 'e5_strict'),
                            ('E5 lenient', e5r, 'e5_lenient'), ('E0 reasoning F1', e0r, 'recovered_f1')):
        print('corr(%-16s, final-answer correctness) = %.3f' % (name, corr([S(rows[k], key) for k in gold], fac)))


if __name__ == '__main__':
    {'validate': validate, 'report': report}[sys.argv[1] if len(sys.argv) > 1 else 'report']()
