"""What the ANSWER FINAL answer check changes in the main store, before WS-G re-scores it (WS-E).

    python -m full_run_28092026.symbolic.rescore_preview      # FREE: RESCORE_PREVIEW.md beside this file; no store written

The answer check changed at ANSWER FINAL in three ways: the symbolic equivalence step on five templates (D5, raises
only), an answer's own last digit bounded at 1% of the target (D11d, lowers only) and an exponent after a bracket no
longer read as a value (D11b). This scores every answered, usable row of the main store's eleven models with
answer.verdict, the call score.score_answer makes for the headline label, with the cached milestones, and counts the
verdicts that move against the stored ones, by model, template and cause (the symbolic step, or the number rule).

The main store was scored with the answer.py these three changes edit (its CONFIG.json records that file's
sha256), so each change's own count comes from adding them in turn: D11b alone (the symbolic step off, the cap
lifted), then D11d, then the symbolic step, which is the new check. A verdict can move at two steps and back, so the
three counts need not sum to the net.

The milestone reader's change (D11a, milestones.numbers) moves E3, not the verdict: the responses whose reached
milestones change (e3_milestones.reach, as score.score_trace calls it, against the stored `e3.reached`), and of
them the ones that keep a missed milestone, whose judge prompt (judge.job) changes and is sent again at the re-score.

It also counts the answers in the error analysis's current B2 sample (the readings expert_kits.reading_status keeps,
top-up rounds included; the kits are local) that the new check no longer scores incorrect: `expert_kits.py
--score`, re-run after the re-score, leaves them out. Counts only.
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402
import e3_milestones as e3  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

OUT = HERE / 'RESCORE_PREVIEW.md'
SCORE_OF = {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}
CHANGES = ('an exponent after a bracket is not a value (D11b)', "the answer's own last digit capped at 1% (D11d)",
           'the symbolic step (D5)')


def in_turn(text: str, item: dict, vals: tuple) -> tuple[list[str], dict]:
    """The verdict with D11b alone, then D11d added, then the symbolic step added (the new check), and its detail."""
    enabled, cap = answer.SYMBOLIC_EQUIVALENCE_TEMPLATES, answer.OWN_DIGIT_CAP
    try:
        answer.SYMBOLIC_EQUIVALENCE_TEMPLATES, answer.OWN_DIGIT_CAP = (), float('inf')
        first = answer.verdict(text, item, score.TOLS['fitted'], vals)[0]
        answer.OWN_DIGIT_CAP = cap
        second = answer.verdict(text, item, score.TOLS['fitted'], vals)[0]
    finally:
        answer.SYMBOLIC_EQUIVALENCE_TEMPLATES, answer.OWN_DIGIT_CAP = enabled, cap
    last, det = answer.verdict(text, item, score.TOLS['fitted'], vals)
    return [first, second, last], det


def main() -> int:
    items = score.pool_items()
    ms = score.milestone_sets(items)
    moved = collections.Counter()
    by_change = collections.defaultdict(collections.Counter)
    by_model = collections.defaultdict(collections.Counter)
    by_template = collections.defaultdict(collections.Counter)
    mean_old, mean_new, n_rows = collections.Counter(), collections.Counter(), collections.Counter()
    e3_moved = collections.defaultdict(collections.Counter)
    new_label = {}
    for key in ROSTER:
        rows = [json.loads(l) for l in (score.SCORES / 'main' / f'{key}.jsonl').read_text(encoding='utf-8').splitlines()]
        texts = score.texts_matching('main', key, rows)
        for r in rows:
            n_rows[key] += 1
            mean_old[key] += r['score']
            if r['status'] != 'answered':
                continue
            it, m = items[r['item_id']], ms[r['item_id']]
            reached, was = [h['reached'] for h in e3.reach(m, texts[r['item_id']])], r['e3']['reached']
            if reached != was:
                c = e3_moved[key]
                c['responses'] += 1
                c['gained'] += sum(a and not b for a, b in zip(reached, was))
                c['lost'] += sum(b and not a for a, b in zip(reached, was))
                c['judge jobs'] += not all(reached)
            if r['unusable']:
                continue
            old = r['answer']['label']
            chain, det = in_turn(texts[r['item_id']], it, tuple(x['value'] for x in m))
            for change, a, b in zip(CHANGES, [old] + chain[:-1], chain):
                if a != b:
                    by_change[change]['raised' if SCORE_OF[b] > SCORE_OF[a] else 'lowered'] += 1
            lab = chain[-1]
            new_label[(key, r['item_id'])] = lab
            mean_new[key] += SCORE_OF[lab]
            if lab != old:
                cause = 'symbolic step' if (det.get('symbolic') or {}).get('equivalent') else 'number rule'
                moved[(old, lab, cause)] += 1
                by_model[key][(old, lab)] += 1
                by_template[r['template_id']][cause] += 1
    L = ['# Re-score preview: what the ANSWER FINAL answer check changes in the main store', '',
         'Generated by `rescore_preview.py` (its docstring defines the counts). No store is written; WS-G\'s re-score '
         'writes them. Answered, usable rows of the eleven models; the score is correct 1, partial 0.5, incorrect 0 and '
         'unusable 0, as stored.', '', '| from | to | cause | verdicts |', '|---|---|---|---:|']
    for (old, new, cause), n in sorted(moved.items(), key=lambda kv: (-kv[1], kv[0])):
        L.append(f'| {old} | {new} | {cause} | {n} |')
    L += [f"| all | | | {sum(moved.values())} |", '', 'Each change on its own, added in this order (a verdict can move at '
          'two of them, so the rows need not sum to the net above):', '',
          '| change | verdicts raised | verdicts lowered |', '|---|---:|---:|']
    for change in CHANGES:
        L.append(f"| {change} | {by_change[change]['raised']} | {by_change[change]['lowered']} |")
    L += ['', '| model | verdicts moved | mean score stored | mean score now | change |', '|---|---:|---:|---:|---:|']
    for key in ROSTER:
        old, new = mean_old[key] / n_rows[key], mean_new[key] / n_rows[key]
        L.append(f"| `{key}` | {sum(by_model[key].values())} | {old:.4f} | {new:.4f} | {new - old:+.4f} |")
    L += ['', '| template | moved by the symbolic step | moved by the number rule |', '|---|---:|---:|']
    for tid, c in sorted(by_template.items(), key=lambda kv: -sum(kv[1].values())):
        L.append(f"| `{tid.replace('template_', '')}` | {c['symbolic step']} | {c['number rule']} |")
    L += ['', 'The milestone reader (D11a) on answered rows: the responses whose reached milestones change, the '
          'milestones newly reached and no longer reached, and the judge jobs whose prompt changes (a milestone is '
          'still missed), which the re-score sends again:', '',
          '| model | responses | milestones newly reached | no longer reached | judge jobs sent again |',
          '|---|---:|---:|---:|---:|']
    for key in ROSTER:
        c = e3_moved[key]
        L.append(f"| `{key}` | {c['responses']} | {c['gained']} | {c['lost']} | {c['judge jobs']} |")
    tot = sum(e3_moved.values(), collections.Counter())
    L.append(f"| all | {tot['responses']} | {tot['gained']} | {tot['lost']} | {tot['judge jobs']} |")
    from full_run_28092026 import expert_kits as K
    if (K.OUT / 'tasks' / 'keyfile.json').exists():
        # the current sample: the readings the store still holds (expert_kits.reading_status), top-up rounds included
        keyfile, shown = json.loads((K.OUT / 'tasks' / 'keyfile.json').read_text(encoding='utf-8')), K.shown_pool(K.OUT)
        returns = K.read_returns(K.RETURNED, K.OUT)[2] if K.RETURNED.exists() else {}
        for top in K.topup_folders(K.OUT):
            keyfile.update(json.loads((top / 'tasks' / 'keyfile.json').read_text(encoding='utf-8')))
            shown.update(K.shown_pool(top))
            found = [d for d in [top / 'returned'] + sorted(top.glob('experts_filled*')) if d.exists()]
            if found:                                  # where the owner files a round's returns, as expert_kits.current
                returns.update(K.read_returns(found[0], top)[2])
        majority = collections.defaultdict(collections.Counter)
        for (_reader, c), r in returns.items():
            if r['answers'].get('category'):           # B2 readings; B1 and B3 rows answer other questions
                majority[c][r['answers']['category'].split(':')[0]] += 1
        store = K.store_rows(sorted({k['model'] for k in keyfile.values() if k.get('kind') == 'error'}))
        sample = {c: (k['model'], k['item_id']) for c, k in keyfile.items() if k.get('kind') == 'error'
                  and K.reading_status(k, shown.get(c), items, store) == K.KEPT}
        leaving = [c for c, mi in sample.items() if new_label.get(mi) not in (None, 'incorrect')]
        per_model = collections.Counter(sample[c][0] for c in leaving)
        said = collections.Counter('no majority' if not majority[c] or majority[c].most_common(1)[0][1] < 2 else
                                   ('"No error"' if majority[c].most_common(1)[0][0] == 'No error' else 'an error')
                                   for c in leaving)
        L += ['', f"Error analysis (B2): of the {len(sample)} answers in its current sample, {len(leaving)} are "
              'no longer scored incorrect by the new check and leave the sample when `expert_kits.py --score` is re-run '
              'after the re-score: ' + (', '.join(f'`{m}` {n}' for m, n in sorted(per_model.items())) or 'none') + '. '
              'The three readers\' majority on them: ' + (', '.join(f'{w} {n}' for w, n in said.most_common()) or 'none')
              + '.']
    L.append('')
    OUT.write_text('\n'.join(L), encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    sys.exit(main())
