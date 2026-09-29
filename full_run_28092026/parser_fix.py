"""D-137: the number-reading fix, measured before against after, on the pilot's 300 expert-labelled
traces and on every full-run trace.

    python -m full_run_28092026.parser_fix      # FREE: writes PARSER_FIX.md beside this file

THE FIX. `answer.values` (the answer check) and `milestones.numbers` (E3) split LaTeX's
`3.04 \\times 10^{5}` and `\\cdot 10^{5}` into 3.04, 10 and 5: the exponent rule ran before
`\\times` was rewritten, and E3's reader never rewrote it. Both now read one number, and both join
thousands grouped as `28\\,570`, `11{,}003` or with a narrow space. Nothing else changed.

THE COMPARISON. The code before the fix and after it are both loaded from git, at PRE_FIX and at
POST_FIX, so the two sides differ in that code alone, run on the same traces with the same
milestones, and a later change to answer.py (D-138) leaves this comparison as it was:
  the pilot     agreement with the experts, answer check and E3, exactly as validate_scorer.py
                computes it (experts_filled_labels/version_2, local)
  the full run  per model, the answer verdicts that change, the answer score, the fully-solved rate
                and E3 coverage, before and after
The changed full-run verdicts are listed, with their answer text, in scores/parser_fix_changes.jsonl,
gitignored with the traces, to be read.
"""
from __future__ import annotations

import collections
import concurrent.futures as cf
import json
import subprocess
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402
import e3_milestones as e3  # noqa: E402
import milestones  # noqa: E402

from full_run_28092026 import score  # noqa: E402

PRE_FIX = '3a7f247954f487c360dd7b3a716bae794d9649f4'   # the last commit of answer.py before D-137
POST_FIX = 'b73c711'                                     # the D-137 commit
CHANGES = score.SCORES / 'parser_fix_changes.jsonl'


def modules_at(commit: str, tag: str):
    """answer.py and milestones.py as they were at a commit, that answer reading that milestones."""
    def src(name):
        return subprocess.run(['git', 'show', f'{commit}:evaluator_pilot_17092026/evaluators/{name}'],
                              cwd=REPO, capture_output=True, check=True).stdout.decode('utf-8')
    ms = types.ModuleType(f'milestones_{tag}')
    ms.__file__ = str(EVAL / 'milestones.py')
    exec(compile(src('milestones.py'), f'milestones_{tag}.py', 'exec'), ms.__dict__)
    saved = sys.modules['milestones']
    sys.modules['milestones'] = ms
    try:
        ans = types.ModuleType(f'answer_{tag}')
        ans.__file__ = str(EVAL / 'answer.py')
        exec(compile(src('answer.py'), f'answer_{tag}.py', 'exec'), ans.__dict__)
    finally:
        sys.modules['milestones'] = saved
    return ans, ms


def pre_fix_modules():
    """The code the pilot's figures were published with."""
    return modules_at(PRE_FIX, 'd120')


def readable(A, text: str) -> bool:
    """score.readable, under either version of the answer check."""
    if A.HEADING.search(text) or A.ANSWER.search(text):
        return True
    tail = text[-A.WINDOW:]
    words = set(A.VERDICT_WORDS) | set(A.LINEAR)
    return bool(A.values(tail)) or any(w in words for w in A.WORD.findall(tail.lower()))


def reached(M, ms: list[dict], text: str) -> list[bool]:
    """e3_milestones.reach, under either version of the number reader."""
    nums = M.numbers(text)
    return [M.scaled_match(m['value'], nums, e3.STEP_TOL) is not None for m in ms]


def verdict3(A, text: str, item: dict, vals: tuple) -> str:
    return A.verdict(text, item, None, vals)[0] if readable(A, text) else 'unusable'


# ------------------------------------------------------------------ the pilot

def pilot(old_a, old_m) -> dict:
    from full_run_28092026 import validate_scorer as V
    import score_against_labels as S
    new_a, new_m = modules_at(POST_FIX, 'd137')
    truth, keyfile, items, texts = V.pilot_inputs()
    lab = {'old': {}, 'new': {}}
    e3c = {'old': collections.Counter(), 'new': collections.Counter()}
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        text, item = texts[key], items[key[1]]
        vals = tuple(m['value'] for m in t['milestones'])
        ms = [{'id': m['id'], 'value': m['value']} for m in t['milestones']]
        for side, A, M in (('old', old_a, old_m), ('new', new_a, new_m)):
            lab[side][code] = A.verdict(text, item, None, vals)[0]
            for hit, tm in zip(reached(M, ms, text), t['milestones']):
                if tm['status'] is None:
                    continue
                truly = tm['status'] == 'reached'
                e3c[side]['tp'] += truly and hit
                e3c[side]['fp'] += (not truly) and hit
                e3c[side]['fn'] += truly and not hit
    codes = sorted(lab['new'])
    exp = {c: truth[c]['final_answer'] for c in codes}
    np_codes = [c for c in codes if exp[c] != 'partial']
    out = {'traces': len(codes)}
    for side in ('old', 'new'):
        L = lab[side]
        p, r, f = S.prf(e3c[side]['tp'], e3c[side]['fp'], e3c[side]['fn'])
        out[side] = {'non_partial': sum((L[c] == 'correct') == (exp[c] == 'correct') for c in np_codes) / len(np_codes),
                     'three_way': sum(L[c] == exp[c] for c in codes) / len(codes),
                     'e3_precision': p, 'e3_recall': r, 'e3_f1': f}
    moved = [c for c in codes if lab['old'][c] != lab['new'][c]]
    out['changed'] = len(moved)
    out['changed_to_experts'] = sum(lab['new'][c] == exp[c] for c in moved)
    out['changed_from_experts'] = sum(lab['old'][c] == exp[c] for c in moved)
    out['transitions'] = collections.Counter(f"{lab['old'][c]} -> {lab['new'][c]} (experts: {exp[c]})"
                                             for c in moved)
    return out


# ------------------------------------------------------------------ the full run

def _model(key: str):
    old_a, old_m = pre_fix_modules()
    new_a, new_m = modules_at(POST_FIX, 'd137')
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    tr = collections.Counter()
    sums = collections.Counter()
    e3n = 0
    changes = []
    for row in score.trace_rows('main', key):
        item, ms = items[row['item_id']], ms_all[row['item_id']]
        vals = tuple(m['value'] for m in ms)
        text = row.get('text') or ''
        if row.get('status') != 'answered':
            old = new = 'unusable'
            cov_old = cov_new = 0.0 if ms else None
        else:
            old, new = verdict3(old_a, text, item, vals), verdict3(new_a, text, item, vals)
            ro, rn = reached(old_m, ms, text), reached(new_m, ms, text)
            cov_old = sum(ro) / len(ro) if ms else None
            cov_new = sum(rn) / len(rn) if ms else None
            sums['e3_sets_changed'] += ro != rn
        tr[(old, new)] += 1
        for side, v in (('old', old), ('new', new)):
            sums[f'score_{side}'] += score.SCORE_OF.get(v, 0.0)
            sums[f'fully_{side}'] += v == 'correct'
        if cov_old is not None:
            e3n += 1
            sums['cov_old'] += cov_old
            sums['cov_new'] += cov_new
        if old != new:
            changes.append({'model_key': key, 'item_id': row['item_id'], 'old': old, 'new': new,
                            'segment': answer.segment(text)})
    n = sum(tr.values())
    return key, {'traces': n, 'changed': sum(v for (o, w), v in tr.items() if o != w),
                 'to_correct': sum(v for (o, w), v in tr.items() if o != w and w == 'correct'),
                 'from_correct': sum(v for (o, w), v in tr.items() if o != w and o == 'correct'),
                 'transitions': {f'{o} -> {w}': v for (o, w), v in sorted(tr.items()) if o != w},
                 'score_old': sums['score_old'] / n, 'score_new': sums['score_new'] / n,
                 'fully_old': sums['fully_old'] / n, 'fully_new': sums['fully_new'] / n,
                 'e3_old': sums['cov_old'] / e3n, 'e3_new': sums['cov_new'] / e3n,
                 'e3_sets_changed': sums['e3_sets_changed']}, changes


def full_run(keys: list[str]) -> dict:
    out, changes = {}, []
    with cf.ProcessPoolExecutor(max_workers=min(8, len(keys))) as pool:
        for key, res, ch in pool.map(_model, keys):
            out[key] = res
            changes.extend(ch)
    with open(CHANGES, 'w', encoding='utf-8', newline='\n') as fh:
        for c in changes:
            fh.write(json.dumps(c, ensure_ascii=False) + '\n')
    return out


def main() -> int:
    from full_run_28092026.run_traces import config as roster
    keys = [m['key'] for m in roster()['models'] if (score.TRACES / f"{m['key']}.jsonl").exists()]
    old_a, old_m = pre_fix_modules()
    p = pilot(old_a, old_m)
    f = full_run(keys)
    L = ['# The number-reading fix, before and after (D-137)', '',
         'Generated by `parser_fix.py`; the fix and the comparison are defined in its docstring. Counts only; '
         f'the code before the fix is `answer.py` and `milestones.py` at `{PRE_FIX[:7]}`, after it at `{POST_FIX}`.', '',
         f"## The pilot's {p['traces']} traces against the experts", '',
         '| figure | before | after |', '|---|---:|---:|',
         f"| answer check, non-partial agreement | {p['old']['non_partial']:.3f} | {p['new']['non_partial']:.3f} |",
         f"| answer check, three-way agreement | {p['old']['three_way']:.3f} | {p['new']['three_way']:.3f} |",
         f"| E3 precision | {p['old']['e3_precision']:.3f} | {p['new']['e3_precision']:.3f} |",
         f"| E3 recall | {p['old']['e3_recall']:.3f} | {p['new']['e3_recall']:.3f} |",
         f"| E3 F1 | {p['old']['e3_f1']:.3f} | {p['new']['e3_f1']:.3f} |", '',
         f"{p['changed']} answer verdicts change. {p['changed_to_experts']} of them now equal the experts' verdict "
         f"and {p['changed_from_experts']} equalled it before:", '']
    L += [f'- {k}: {v}' for k, v in sorted(p['transitions'].items())]
    L += ['', '## The full run', '',
          'Per model, over its 2,250 final traces: the verdicts that change (unusable counts as a verdict), and '
          'the answer score, fully-solved rate and E3 coverage before and after. E3 coverage is the mean over '
          'traces whose item has milestones.', '',
          '| model | verdicts changed | to correct | from correct | answer score | fully solved | E3 coverage | '
          'traces whose E3 hits change |', '|---|---:|---:|---:|---:|---:|---:|---:|']
    for k, r in f.items():
        L.append(f"| `{k}` | {r['changed']} | {r['to_correct']} | {r['from_correct']} | "
                 f"{r['score_old']:.3f} -> {r['score_new']:.3f} | {r['fully_old']:.3f} -> {r['fully_new']:.3f} | "
                 f"{r['e3_old']:.3f} -> {r['e3_new']:.3f} | {r['e3_sets_changed']} |")
    L += ['', 'Transitions per model:', '']
    L += [f"- `{k}`: " + '; '.join(f'{t} {v}' for t, v in r['transitions'].items()) for k, r in f.items()]
    (HERE / 'PARSER_FIX.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
