"""Does the full-run scorer reproduce what was validated? Gold, and the pilot's expert agreement.

    python -m full_run_28092026.validate_scorer      # writes SCORER_VALIDATION.md beside this file

1. GOLD. score.gold(): every gold solution scored correct at all three tolerances, none unusable,
   every milestone found, no digit-rule flag.
2. THE PILOT'S 300 TRACES, scored by score.score_trace exactly as the full run is, against the
   experts' labels (experts_filled_labels/version_2, local), with the pilot's own conventions:
     answer check   agreement with the experts' final-answer verdict: on the traces they judged
                    fully right or wrong, and three-way (RESULTS_X1 Finding 1b); the milestones
                    are the ones the labels carry, as answer_check.py passes them
     E3             reached against the experts' "obtained", per milestone
                    (score_against_labels.report_milestones; RESULTS_X1 Finding 2)
     digit rule     per step against the experts' step labels, a trace skipped when its step count
                    differs from the labels' (digit_rule.py, its e4 row; RESULTS_X1 Finding 5)
   Each figure is also recomputed with the pilot's own functions on the same inputs, so a
   difference would be the scorer's plumbing, not a change in the evaluators.
Nothing here calls a model or the network. It writes figures only.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PILOT = REPO / 'evaluator_pilot_17092026'
for p in (str(REPO), str(PILOT / 'evaluators'), str(PILOT / 'annotation'), str(PILOT / 'analysis'),
          str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer as A  # noqa: E402
import arith  # noqa: E402
import e2_prm  # noqa: E402
import score_against_labels as S  # noqa: E402

from full_run_28092026 import score  # noqa: E402

LABELS = PILOT / 'experts_filled_labels' / 'version_2' / 'labels'
PUBLISHED = {   # RESULTS_X1 Findings 1b, 2 and 5; digit_rule.py's e4 rows
    'answer, non-partial agreement': 0.947, 'answer, three-way agreement': 0.893,
    'E3 precision': 0.926, 'E3 recall': 0.915, 'E3 F1': 0.921,
    'digit rule, all traces, precision': 0.825, 'digit rule, all traces, recall': 0.255,
    'digit rule, hard case, precision': 0.750, 'digit rule, hard case, recall': 0.320,
    'digit rule, hard case, F1': 0.449,
}


def pilot_inputs():
    truth, _ = S.build_truth(S.read_labels(str(LABELS)), S.read_consensus(str(LABELS)))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(Path(S.TASKS) / 'keyfile.jsonl', encoding='utf-8'))}
    items = {json.loads(l)['item_id']: json.loads(l)
             for l in (PILOT / 'slice' / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()}
    texts = {}
    for f in (PILOT / 'traces').glob('*.jsonl'):
        for line in f.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            if r.get('ok'):
                texts[(r['model_key'], r['item_id'])] = r['text']
    return truth, keyfile, items, texts


def prf(c):
    p, r, f = S.prf(c['tp'], c['fp'], c['fn'])
    return round(p, 3), round(r, 3), round(f, 3)


def main() -> int:
    g = score.gold(8)
    truth, keyfile, items, texts = pilot_inputs()
    rows, pilot_labels = {}, {}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        it = dict(items[key[1]])
        it.setdefault('instance_index', it.get('instance_index', 0))
        for k in ('branch', 'domain', 'level', 'answer_type', 'repeat_of', 'single_path'):
            it.setdefault(k, None)
        ms = [{'id': m['id'], 'value': m['value']} for m in t['milestones']]
        rows[code] = score.score_trace(it, {'status': 'answered', 'text': texts[key]}, ms)
        pilot_labels[code] = A.verdict(texts[key], items[key[1]], None, tuple(m['value'] for m in ms))[0]

    # answer check
    codes = sorted(rows)
    same_as_pilot = sum(rows[c]['answer']['label'] == pilot_labels[c] for c in codes)
    # answer_check.py's definition: correct-or-not, over the traces the experts did not call partial
    np_codes = [c for c in codes if truth[c]['final_answer'] != 'partial']
    agree_np = sum((rows[c]['answer']['label'] == 'correct') == (truth[c]['final_answer'] == 'correct')
                   for c in np_codes) / len(np_codes)
    agree_3 = sum(rows[c]['answer']['label'] == truth[c]['final_answer'] for c in codes) / len(codes)

    # E3 against "obtained"
    m = Counter()
    for c in codes:
        for reached, tm in zip(rows[c]['e3']['reached'], truth[c]['milestones']):
            if tm['status'] is None:
                continue
            truly = tm['status'] == 'reached'
            m['tp'] += truly and reached
            m['fp'] += (not truly) and reached
            m['fn'] += truly and not reached

    # digit rule, per step, with digit_rule.py's own flag for comparison
    step_c = {'all': Counter(), 'hard': Counter()}
    own = {'all': Counter(), 'hard': Counter()}
    scored = skipped = 0
    for c in codes:
        steps = rows[c]['steps']
        if len(steps) != len(truth[c]['steps']):
            skipped += 1
            continue
        scored += 1
        text_steps = e2_prm.steps_of(texts[keyfile[c]])
        for scope in ('all', 'hard'):
            if scope == 'hard' and truth[c]['final_answer'] != 'correct':
                continue
            for s, ts, st in zip(steps, text_steps, truth[c]['steps']):
                flag = s['digit_flags'] > 0
                pilot_flag = any(not x.ok_digit for x in arith.check(ts).claims)
                for cnt, f in ((step_c[scope], flag), (own[scope], pilot_flag)):
                    if st['label'] == 'incorrect':
                        cnt['tp' if f else 'fn'] += 1
                    elif st['label'] in ('correct', 'alternative_correct'):
                        cnt['fp' if f else 'tn'] += 1

    ep, er, ef = prf(m)
    ap_, ar, af = prf(step_c['all'])
    hp, hr, hf = prf(step_c['hard'])
    got = {'answer, non-partial agreement': round(agree_np, 3), 'answer, three-way agreement': round(agree_3, 3),
           'E3 precision': ep, 'E3 recall': er, 'E3 F1': ef,
           'digit rule, all traces, precision': ap_, 'digit rule, all traces, recall': ar,
           'digit rule, hard case, precision': hp, 'digit rule, hard case, recall': hr,
           'digit rule, hard case, F1': hf}
    L = ['# Scorer validation', '', 'Generated by `validate_scorer.py`; the checks are defined in its docstring.', '',
         '## Gold, all 2,250 items', '', '| check | result |', '|---|---|',
         f"| scored correct at all three tolerances | {g['answer_correct_all_three_tols']} of {g['items']} |",
         f"| unusable | {g['unusable']} |",
         f"| every milestone found | {g['e3_complete']}; {g['e3_no_milestones']} items have no milestones |",
         f"| digit-rule flags | {g['digit_flags']} of {g['claims_checked']} claims in {g['steps']} steps |", '',
         f'## The pilot\'s {len(codes)} traces against the experts', '',
         f'Answer labels identical to the pilot\'s own `answer.verdict` call: {same_as_pilot} of {len(codes)}. '
         f'Step-level scoring covers {scored} traces; {skipped} are skipped for a step count that differs '
         f'from the labels\', as in `digit_rule.py`, whose own flags give tp/fp/fn '
         f"{own['all']['tp']}/{own['all']['fp']}/{own['all']['fn']} over all traces and "
         f"{own['hard']['tp']}/{own['hard']['fp']}/{own['hard']['fn']} in the hard case, against the scorer's "
         f"{step_c['all']['tp']}/{step_c['all']['fp']}/{step_c['all']['fn']} and "
         f"{step_c['hard']['tp']}/{step_c['hard']['fp']}/{step_c['hard']['fn']}.", '',
         '| figure | published | scorer | same |', '|---|---:|---:|---|']
    ok = True
    for k, v in PUBLISHED.items():
        same = abs(got[k] - v) < 5e-4
        ok &= same
        L.append(f'| {k} | {v:.3f} | {got[k]:.3f} | {"yes" if same else "NO"} |')
    L += ['', 'All figures reproduced.' if ok else 'A figure differs: see the rows marked NO.']
    (HERE / 'SCORER_VALIDATION.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
