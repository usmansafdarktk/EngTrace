"""Error attribution validated by type: each component's flags against the experts' step error types, on the pilot's
300 labelled traces (docs/EVALUATION_NEXT_STEPS.md A6; D-176).

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/attribution.py       # FREE, offline: RESULTS_ATTRIBUTION.md, analysis/out/attribution.json

WHY. The full run's Q3 reports, on the answered wrong answers, the share with a digit-rule flag, with a milestone E5
rules MISSING and with a step the router's judge flags, and the notes say to read those as shares of traces flagged,
never as counts of error types. This measures what each flag is, on the only traces with a ground truth: the experts'
step labels, each incorrect step typed calculation, conceptual or unsupported (annotation/guide.md), with the
adjudication applied.

WHAT IT CROSS-TABULATES.
  steps   the digit rule as shipped (`digit_rule.flagged(step, 'e4')`) and the router's judge by its category
          (Calculation Error, Conceptual Error, Other; from the batched run's stored replies, `router_batched.verdicts`)
          against the expert label of the same step (correct or alternative correct; incorrect, by type). Precision by
          type: of a component's flags, the share that land on an expert-incorrect step of each type. Recall by type:
          of the expert-incorrect steps of each type, the share the component flags. Over all traces and inside
          correct-answer traces.
  traces  a MISSING milestone under E5 (the pilot's E5 rows), any digit flag and any router flag, against whether the
          trace holds an expert-incorrect step of each type, split by the experts' final-answer verdict.
The step spaces are the labels' (the framework's non-empty steps, `e2_prm.steps_of`); a trace whose split does not
match the labels' step count is left out and counted.
"""
from __future__ import annotations

import collections
import glob
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_REPO = os.path.dirname(_PILOT)
for p in (os.path.join(_PILOT, 'annotation'), os.path.join(_PILOT, 'evaluators'), _HERE, os.path.join(_REPO, 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import digit_rule as DR  # noqa: E402
import e2_prm  # noqa: E402
import router_batched as RB  # noqa: E402
import score_against_labels as S  # noqa: E402

OUT_MD = os.path.join(_PILOT, 'RESULTS_ATTRIBUTION.md')
OUT_JSON = os.path.join(_HERE, 'out', 'attribution.json')
LABELS = os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels')
TYPES = ('calculation', 'conceptual', 'unsupported')
CATS = ('Calculation Error', 'Conceptual Error', 'Other')


def rows(path):
    return [json.loads(l) for l in open(path, encoding='utf-8') if l.strip()]


def current(d):
    rs = [r for f in glob.glob(os.path.join(_PILOT, 'scores', d, '*.jsonl')) for r in rows(f) if not r.get('error')]
    latest = max(rs, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rs if r['config_sha256'] == latest}


def expert_type(step: dict):
    """'correct' for a correct or alternative-correct step, the error type for an incorrect one, None otherwise."""
    if step['label'] in ('correct', 'alternative_correct'):
        return 'correct'
    if step['label'] == 'incorrect':
        return step.get('error_type') if step.get('error_type') in TYPES else 'untyped'
    return None


def main() -> int:
    truth, _ = S.build_truth(S.read_labels(LABELS), S.read_consensus(LABELS))
    keyfile = {r['code']: (r['model_key'], r['item_id']) for r in rows(os.path.join(S.TASKS, 'keyfile.jsonl'))}
    texts = RB.texts()
    _probe, judged, _ = RB.verdicts()
    e5 = current('e5')
    per_step = []          # (code, j, expert type, digit flag, judge category or None, judge sent)
    skipped = 0
    for code, t in truth.items():
        text = texts.get(keyfile[code])
        if not text:
            skipped += 1
            continue
        steps = e2_prm.steps_of(text)
        if len(steps) != len(t['steps']):
            skipped += 1
            continue
        rv = judged.get(('labelled', code, 'trace'))
        for j, st in enumerate(t['steps']):
            ty = expert_type(st)
            if ty is None:
                continue
            cat = rv.get(j) if rv else None
            sent = rv is not None and j in rv
            per_step.append((code, j, ty, bool(DR.flagged(steps[j], 'e4')), cat, sent))
    right = {c for c in truth if truth[c]['final_answer'] == 'correct'}

    def judge_cat(cat):
        if cat is None:
            return None
        c = cat.lower()
        return 'Calculation Error' if 'calculation' in c else 'Conceptual Error' if 'conceptual' in c else 'Other' if 'error' in c or 'other' in c else None

    res = {'steps_labelled': len(per_step), 'traces_skipped': skipped, 'scopes': {}}
    for scope, pick in (('all', lambda c: True), ('correct_answer', lambda c: c in right)):
        rows_ = [r for r in per_step if pick(r[0])]
        base = collections.Counter(r[2] for r in rows_)
        comps = {'digit rule': [r for r in rows_ if r[3]]}
        for cat in CATS:
            comps[f'router judge: {cat}'] = [r for r in rows_ if judge_cat(r[4]) == cat]
        comps['router judge: any flag'] = [r for r in rows_ if judge_cat(r[4]) is not None]
        comps['router: rule or judge'] = [r for r in rows_ if r[3] or judge_cat(r[4]) is not None]
        out = {'steps_by_expert_type': dict(base), 'components': {}}
        for name, flagged in comps.items():
            by = collections.Counter(r[2] for r in flagged)
            n = len(flagged)
            out['components'][name] = {
                'flags': n,
                'landing_on': {ty: by[ty] for ty in ('correct', 'untyped') + TYPES},
                'precision_any_error': (n - by['correct']) / n if n else None,
                'precision_by_type': {ty: by[ty] / n if n else None for ty in TYPES},
                'recall_by_type': {ty: by[ty] / base[ty] if base[ty] else None for ty in TYPES},
                'recall_correct_flagged': by['correct'] / base['correct'] if base['correct'] else None}
        res['scopes'][scope] = out
    # trace level
    by_code = collections.defaultdict(list)
    for r in per_step:
        by_code[r[0]].append(r)
    tr = {}
    for scope, pick in (('all', lambda c: True), ('correct_answer', lambda c: c in right), ('wrong_answer', lambda c: c not in right)):
        codes = [c for c in by_code if pick(c)]
        has = {c: {ty: any(r[2] == ty for r in by_code[c]) for ty in TYPES} for c in codes}
        flags = {}
        for c in codes:
            e = e5.get(keyfile[c])
            missing = bool(e and any(m.get('source') == 'MISSING' for m in e['meta']['milestones']))
            flags[c] = {'E5 MISSING': missing, 'digit rule': any(r[3] for r in by_code[c]),
                        'router judge': any(judge_cat(r[4]) is not None for r in by_code[c])}
        out = {'traces': len(codes), 'with_expert_error_by_type': {ty: sum(has[c][ty] for c in codes) for ty in TYPES},
               'components': {}}
        for comp in ('E5 MISSING', 'digit rule', 'router judge'):
            fl = [c for c in codes if flags[c][comp]]
            out['components'][comp] = {
                'traces_flagged': len(fl),
                'of_flagged_with_error_type': {ty: sum(has[c][ty] for c in fl) / len(fl) if fl else None for ty in TYPES},
                'of_flagged_with_no_incorrect_step': sum(not any(has[c].values()) for c in fl) / len(fl) if fl else None,
                'recall_by_type': {ty: (sum(1 for c in fl if has[c][ty]) / sum(has[c][ty] for c in codes))
                                   if sum(has[c][ty] for c in codes) else None for ty in TYPES}}
        tr[scope] = out
    res['traces'] = tr
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    json.dump(res, open(OUT_JSON, 'w', encoding='utf-8'), indent=1)
    text = render(res)
    open(OUT_MD, 'w', encoding='utf-8', newline='\n').write(text)
    print(text)
    return 0


def f3(v):
    return '-' if v is None else f'{v:.3f}'


def render(res) -> str:
    L = ['# Error attribution validated by type: the components\' flags against the experts\' step error types', '',
         'Generated by `analysis/attribution.py`; the cross-tabulations are defined in its docstring. The pilot\'s 300 '
         f"expert-labelled traces (version 2 with the adjudication); {res['steps_labelled']} labelled steps in the "
         f"{300 - res['traces_skipped']} traces whose step split matches the labels' ({res['traces_skipped']} left out). "
         'A flag\'s precision by type is the share of the component\'s flags that land on an expert-incorrect step of that '
         'type; its recall by type is the share of such steps it flags (D-176; next steps A6).', '']
    for scope, title in (('all', 'All traces'), ('correct_answer', 'Inside correct-answer traces (the hard case)')):
        s = res['scopes'][scope]
        b = s['steps_by_expert_type']
        L += [f'## Steps, {title.lower()}', '',
              f"Expert steps: correct or alternative correct {b.get('correct', 0)}; incorrect: calculation {b.get('calculation', 0)}, "
              f"conceptual {b.get('conceptual', 0)}, unsupported {b.get('unsupported', 0)}, untyped {b.get('untyped', 0)}.", '',
              '| component | flags | land on: correct | calculation | conceptual | unsupported | precision, any error | '
              'precision: calculation / conceptual | recall: calculation / conceptual / unsupported | correct steps flagged |',
              '|---|---:|---:|---:|---:|---:|---:|---|---|---:|']
        for name, c in s['components'].items():
            lo, pb, rb = c['landing_on'], c['precision_by_type'], c['recall_by_type']
            L.append(f"| {name} | {c['flags']} | {lo['correct']} | {lo['calculation']} | {lo['conceptual']} | {lo['unsupported']} | "
                     f"{f3(c['precision_any_error'])} | {f3(pb['calculation'])} / {f3(pb['conceptual'])} | "
                     f"{f3(rb['calculation'])} / {f3(rb['conceptual'])} / {f3(rb['unsupported'])} | {f3(c['recall_correct_flagged'])} |")
        L.append('')
    L += ['## Traces', '', 'Per trace: whether the component flags it (E5: a milestone the judge rules MISSING; the digit rule: any '
          'flagged step; the router\'s judge: any flagged step) against whether the experts mark an incorrect step of each type '
          'in it. "Of flagged" columns are shares of the flagged traces; recall is the share of the traces with an error of that '
          'type the component flags.', '',
          '| traces | n | with a calculation / conceptual / unsupported error | component | flagged | of flagged: calculation / '
          'conceptual / unsupported | of flagged: no incorrect step | recall: calculation / conceptual |',
          '|---|---:|---|---|---:|---|---:|---|']
    for scope, title in (('all', 'all'), ('correct_answer', 'correct answer'), ('wrong_answer', 'wrong answer')):
        t = res['traces'][scope]
        w = t['with_expert_error_by_type']
        for comp, c in t['components'].items():
            o, r = c['of_flagged_with_error_type'], c['recall_by_type']
            L.append(f"| {title} | {t['traces']} | {w['calculation']} / {w['conceptual']} / {w['unsupported']} | {comp} | "
                     f"{c['traces_flagged']} | {f3(o['calculation'])} / {f3(o['conceptual'])} / {f3(o['unsupported'])} | "
                     f"{f3(c['of_flagged_with_no_incorrect_step'])} | {f3(r['calculation'])} / {f3(r['conceptual'])} |")
    L += ['', '## How to read the full run\'s attribution table through this', '',
          'Q3\'s "what points at a wrong answer" gives, per model, the share of answered wrong answers with a digit-rule flag, '
          'an E5 MISSING milestone and a router-judge flag. The tables above say what such a flag is on labelled traces: the '
          'precision column is the share of flags that sit on a step the experts call incorrect, and the type columns split '
          'those by cause. A share of flagged traces multiplied by the component\'s precision is the share with a confirmed '
          'error of some type; it is never a count of calculation or conceptual errors, because recall differs by type and the '
          'conceptual steps are few (58 of 388 incorrect steps).', '']
    return '\n'.join(L)


if __name__ == '__main__':
    sys.exit(main())
