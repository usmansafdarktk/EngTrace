"""Gold validation of the full run's deterministic evaluators, on all 2,250 frozen items.

    python -m full_run_28092026.gold_validation     # writes GOLD_VALIDATION.md and gold_validation.json

The gold solutions are correct by construction and certified by the experts, so on gold each
evaluator has a known right answer, and anything it gets wrong is a gap in the evaluator:

  1. reproduction   every item regenerates byte-identically from its template at its seed;
                    `milestones.build` refuses otherwise. The pool is the templates' output.
  2. answer check   every gold solution is `correct` against itself, in all six answer types.
  3. E3             every milestone is found in its own gold (coverage 1.0). Against a sibling
                    item of the same template, the coverage is the rate at which E3 credits a
                    number the trace never derived (raw: shared values count as hits).
  4. digit rule     as E4 ships it (`arith.check`, `ok_digit`), no arithmetic claim in any gold
                    solution is flagged.
  5. milestones     every template has some: with none or one, E3 is an answer check.

The evaluators are the pilot's, imported unmodified from evaluator_pilot_17092026/evaluators/.
Nothing here calls a model or the network. It reads pool/ (local) and manifest.jsonl, and writes
counts, template ids and item ids only; the text of anything that fails goes to
pool/_gold_validation_details.txt, which is gitignored with the pool.
"""
from __future__ import annotations

import collections
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL)):
    if p not in sys.path:
        sys.path.insert(0, p)

import tests.template_integrity.core as core  # noqa: E402

_REFS = core.discover()
core.discover = lambda *a, **k: _REFS          # milestones.build re-discovers per item

import answer  # noqa: E402
import arith  # noqa: E402
import e3_milestones as e3  # noqa: E402
import milestones  # noqa: E402

POOL = HERE / 'pool'
DETAILS = POOL / '_gold_validation_details.txt'


def load_items() -> list[dict]:
    rows = {json.loads(l)['item_id']: json.loads(l)
            for l in (HERE / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()}
    items = []
    for path in sorted(POOL.rglob('*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            r['template_id'] = 'template_' + r['id']
            r['answer_type'] = rows[r['item_id']]['answer_type']
            items.append(r)
    return items


def main() -> int:
    items = load_items()
    detail = []
    by_t = collections.defaultdict(list)
    for it in items:
        by_t[it['template_id']].append(it)

    built, not_reproduced = {}, []
    for it in items:
        try:
            built[it['item_id']] = milestones.build(it)['milestones']
        except ValueError:
            not_reproduced.append(it['item_id'])

    ans = collections.Counter()
    ans_fail = collections.defaultdict(list)
    e3_short = collections.defaultdict(list)
    null = collections.defaultdict(list)
    digit_flags = collections.defaultdict(list)
    claims_checked = 0
    for tid, its in by_t.items():
        for k, it in enumerate(its):
            ms = built.get(it['item_id'])
            if ms is None:
                continue
            label, info = answer.verdict(it['solution'], it, ms=[m['value'] for m in ms])
            ans[(it['answer_type'], label)] += 1
            if label != 'correct':
                ans_fail[tid].append(it['item_id'])
                detail.append(f"ANSWER {it['item_id']} {label}: matched {info.get('matched')} "
                              f"of {info.get('of')}; targets {info.get('targets')}")
            hits = e3.reach(ms, it['solution'])
            if ms and not all(h['reached'] for h in hits):
                e3_short[tid].append(it['item_id'])
                detail.append(f"E3 {it['item_id']}: missed "
                              f"{[h['id'] for h in hits if not h['reached']]}")
            sib = its[(k + 1) % len(its)]
            if ms and sib['question'] != it['question']:
                sh = e3.reach(ms, sib['solution'])
                null[tid].append(sum(h['reached'] for h in sh) / len(sh))
            rep = arith.check(it['solution'])
            claims_checked += len(rep.claims)
            for c in rep.claims:
                if not c.ok_digit:
                    digit_flags[tid].append(it['item_id'])
                    detail.append(f"DIGIT {it['item_id']} line {c.line}: {c.left} = {c.right}")

    counts = {tid: [len(built[i['item_id']]) for i in its if i['item_id'] in built]
              for tid, its in by_t.items()}
    thin = sorted(t for t, c in counts.items() if c and min(c) <= 1)
    none_ = sorted(t for t, c in counts.items() if c and max(c) == 0)
    null_all = [x for v in null.values() for x in v]
    by_type = collections.defaultdict(lambda: collections.Counter())
    for (typ, label), n in ans.items():
        by_type[typ][label] += n

    out = {
        'items': len(items), 'templates': len(by_t),
        'not_reproduced': not_reproduced,
        'answer_check': {typ: dict(c) for typ, c in sorted(by_type.items())},
        'answer_check_failures': {t: v for t, v in sorted(ans_fail.items())},
        'e3_gold_short': {t: v for t, v in sorted(e3_short.items())},
        'e3_null_mean': round(statistics.mean(null_all), 4) if null_all else None,
        'e3_null_by_template_max': (max((round(statistics.mean(v), 4), t) for t, v in null.items())
                                    if null else None),
        'digit_rule': {'claims_checked': claims_checked,
                       'flagged_claims': sum(len(v) for v in digit_flags.values()),
                       'templates_flagged': {t: len(v) for t, v in sorted(digit_flags.items())}},
        'milestones_per_item': {'min': min(min(c) for c in counts.values() if c),
                                'median': statistics.median([x for c in counts.values() for x in c]),
                                'max': max(max(c) for c in counts.values() if c)},
        'templates_with_an_item_of_at_most_one_milestone': thin,
        'templates_without_milestones': none_,
    }
    (HERE / 'gold_validation.json').write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8',
                                               newline='\n')
    DETAILS.write_text('\n'.join(detail) + '\n', encoding='utf-8')

    L = ['# Gold validation of the deterministic evaluators', '',
         'Generated by `gold_validation.py` over the frozen pool; the checks are defined in its '
         'docstring. Counts and ids only.', '', '| check | result |', '|---|---|',
         f"| items regenerating byte-identically from their templates | {len(items) - len(not_reproduced)} of {len(items)} |"]
    for typ, c in sorted(by_type.items()):
        tot = sum(c.values())
        L.append(f"| answer check, {typ}: gold scored correct | {c.get('correct', 0)} of {tot} |")
    L += [f"| templates where the answer check fails some gold | {len(ans_fail)} |",
          f"| templates where E3 misses a milestone in its own gold | {len(e3_short)} |",
          f"| E3 coverage of a sibling item's gold, mean (the null) | {out['e3_null_mean']} |",
          f"| digit rule: claims checked on gold, claims flagged | {claims_checked}, {out['digit_rule']['flagged_claims']} |",
          f"| templates with a flagged gold claim | {len(digit_flags)} |",
          f"| milestones per item: min, median, max | {out['milestones_per_item']['min']}, "
          f"{out['milestones_per_item']['median']}, {out['milestones_per_item']['max']} |",
          f"| templates with an item of at most one milestone | {len(thin)} |", '']
    for title, d in (('Answer check fails on gold', ans_fail), ('E3 misses a milestone in its own gold', e3_short),
                     ('Digit rule flags a gold claim', digit_flags)):
        if d:
            L.append(f'## {title}')
            L.append('')
            for t, v in sorted(d.items()):
                L.append(f"- `{t}`: {len(v)} item(s)")
            L.append('')
    if thin:
        L += ['## Templates with an item of at most one milestone', '']
        L += [f"- `{t}`: counts {sorted(set(counts[t]))}" for t in thin] + ['']
    (HERE / 'GOLD_VALIDATION.md').write_text('\n'.join(L), encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
