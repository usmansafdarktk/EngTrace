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
3. THE CODE CHANGED SINCE. The answer check and E3 read numbers differently since D-137, so the
   published figures are reproduced with the code they were published with (answer.py and
   milestones.py at parser_fix.PRE_FIX, loaded from git), and the figures the current code gives
   are reported beside them. The digit rule reads more of the traces since D-156, so its published
   figures are reproduced the same way, with arith.py at digit_fix.PRE_FIX. The check passes when the
   old code reproduces every published figure and the gold is clean.
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

from full_run_28092026 import parser_fix, score  # noqa: E402

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


def answer_figures(labels, truth, codes):
    """answer_check.py's definitions: correct-or-not over the traces the experts did not call
    partial, and three-way over all."""
    np_codes = [c for c in codes if truth[c]['final_answer'] != 'partial']
    return (round(sum((labels[c] == 'correct') == (truth[c]['final_answer'] == 'correct')
                      for c in np_codes) / len(np_codes), 3),
            round(sum(labels[c] == truth[c]['final_answer'] for c in codes) / len(codes), 3))


def e3_counts(hits, truth, codes):
    """E3 against the experts' "obtained", per milestone (score_against_labels.report_milestones)."""
    m = Counter()
    for c in codes:
        for reached, tm in zip(hits[c], truth[c]['milestones']):
            if tm['status'] is None:
                continue
            truly = tm['status'] == 'reached'
            m['tp'] += truly and reached
            m['fp'] += (not truly) and reached
            m['fn'] += truly and not reached
    return m


def main() -> int:
    g = score.gold(8)
    truth, keyfile, items, texts = pilot_inputs()
    old_a, old_m = parser_fix.pre_fix_modules()
    rows, pilot_labels, old_labels, old_hits = {}, {}, {}, {}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        it = dict(items[key[1]])
        it.setdefault('instance_index', it.get('instance_index', 0))
        for k in ('branch', 'domain', 'level', 'answer_type', 'repeat_of', 'single_path'):
            it.setdefault(k, None)
        ms = [{'id': m['id'], 'value': m['value']} for m in t['milestones']]
        vals = tuple(m['value'] for m in ms)
        rows[code] = score.score_trace(it, {'status': 'answered', 'text': texts[key]}, ms)
        pilot_labels[code] = A.verdict(texts[key], items[key[1]], None, vals)[0]
        old_labels[code] = old_a.verdict(texts[key], items[key[1]], None, vals)[0]
        old_hits[code] = parser_fix.reached(old_m, ms, texts[key])

    # answer check, now and with the code the figures were published with
    codes = sorted(rows)
    same_as_pilot = sum(rows[c]['answer']['label'] == pilot_labels[c] for c in codes)
    agree_np, agree_3 = answer_figures({c: rows[c]['answer']['label'] for c in codes}, truth, codes)
    old_np, old_3 = answer_figures(old_labels, truth, codes)

    # E3 against "obtained"
    m = e3_counts({c: rows[c]['e3']['reached'] for c in codes}, truth, codes)
    m_old = e3_counts(old_hits, truth, codes)

    # digit rule, per step, with digit_rule.py's own flag for comparison, and with arith.py as the
    # figures were published (before D-156) to reproduce them
    from full_run_28092026 import digit_fix
    pre_arith = digit_fix.pre_fix_arith()
    step_c = {'all': Counter(), 'hard': Counter()}
    own = {'all': Counter(), 'hard': Counter()}
    then_c = {'all': Counter(), 'hard': Counter()}
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
                pre_flag = any(not x.ok_digit for x in pre_arith.check(ts).claims)
                for cnt, f in ((step_c[scope], flag), (own[scope], pilot_flag), (then_c[scope], pre_flag)):
                    if st['label'] == 'incorrect':
                        cnt['tp' if f else 'fn'] += 1
                    elif st['label'] in ('correct', 'alternative_correct'):
                        cnt['fp' if f else 'tn'] += 1

    ep, er, ef = prf(m)
    op, orr, of = prf(m_old)
    ap_, ar, af = prf(step_c['all'])
    hp, hr, hf = prf(step_c['hard'])
    digit = {'digit rule, all traces, precision': ap_, 'digit rule, all traces, recall': ar,
             'digit rule, hard case, precision': hp, 'digit rule, hard case, recall': hr,
             'digit rule, hard case, F1': hf}
    tap, tar, _ = prf(then_c['all'])
    thp, thr, thf = prf(then_c['hard'])
    digit_then = {'digit rule, all traces, precision': tap, 'digit rule, all traces, recall': tar,
                  'digit rule, hard case, precision': thp, 'digit rule, hard case, recall': thr,
                  'digit rule, hard case, F1': thf}
    got = {'answer, non-partial agreement': agree_np, 'answer, three-way agreement': agree_3,
           'E3 precision': ep, 'E3 recall': er, 'E3 F1': ef, **digit}
    then = {'answer, non-partial agreement': old_np, 'answer, three-way agreement': old_3,
            'E3 precision': op, 'E3 recall': orr, 'E3 F1': of, **digit_then}
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
         f'"Published code" is the scorer with `answer.py` and `milestones.py` as they were published, at '
         f'`{parser_fix.PRE_FIX[:7]}`; "now" is the current code: LaTeX numbers read (D-137, `PARSER_FIX.md`), '
         'no credit from a subscript or a bare 0 or 1 (D-138, `MATCH_AUDIT.md`), a verdict word decided within '
         'its family (D-139, `WORD_AUDIT.md`) and the last-digit windows decided by the rule rather than by binary '
         'rounding (D-147, `BOUNDARY_AUDIT.md`), and pi-fractions, a scalar gold line\'s second unit and a '
         'pi-fraction of two computed milestones read (D-169, `ANSWER_FORM_FIX.md`); the last four change none of '
         'these verdicts. The digit rule\'s '
         f'published figures are reproduced with `arith.py` at `{digit_fix.PRE_FIX[:7]}`; "now" is the rule as '
         'fixed twice after a domain expert read the full run\'s flags (D-156, D-159; `DIGIT_FIX.md`, '
         '`DIGIT_FIX_2.md`) and corrected after two review agents checked both fixes (D-160; `DIGIT_FIX_3.md`).', '',
         '| figure | published | published code | reproduced | now |', '|---|---:|---:|---|---:|']
    ok = g['answer_correct_all_three_tols'] == g['items'] and not g['unusable'] and not g['digit_flags']
    for k, v in PUBLISHED.items():
        same = abs(then[k] - v) < 5e-4
        ok &= same
        L.append(f'| {k} | {v:.3f} | {then[k]:.3f} | {"yes" if same else "NO"} | {got[k]:.3f} |')
    L += ['', 'The gold is clean and the published code reproduces every figure.' if ok
          else 'A check failed: see the gold table and the rows marked NO.']
    (HERE / 'SCORER_VALIDATION.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
