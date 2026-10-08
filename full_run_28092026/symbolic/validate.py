"""E3 (WS-E): the symbolic equivalence check against the experts, and its reach over the main store.

    python -m full_run_28092026.symbolic.validate              # FREE: validation.json and SYMBOLIC_CHECK.md
    python -m full_run_28092026.symbolic.validate --selftest   # FREE: the measures on synthetic grades

WHAT IT MEASURES, for each template equivalence.SPECS covers:
  golds      every gold of the template's 15 items is equivalent to itself, and each other item's gold (another
             question) is not; a pair that matches is listed (two values that coincide within the precision shown)
  grades     against the E2 grades (grades.json, written by score_grades.py): on every graded item (the grade agreed,
             not split), positive = "equivalent", and "not equivalent" or "unreadable" negative; precision and recall
             per template, on all graded items and on the ones the check can move (the answer check's incorrect and
             partial verdicts; it never lowers a correct one), under the check's window (half a unit of the gold's
             last digits) and, as a sensitivity, a whole unit. A template is ENABLED when, on the verdicts it can move,
             precision and recall are both at least 0.9, it moves at least one graded verdict, and every disagreement
             on it is explained (EXPLAINED below; an unexplained one keeps it off). answer.SYMBOLIC_EQUIVALENCE_TEMPLATES
             names the enabled templates.
  readings   the earlier readings that touch these templates (expert_request/, local), made before the check existed:
             B1, the experts' verdict on the final answer (two readers; correct is positive, any other verdict
             negative, a split skipped) and B2, the error category of a wrong answer (three readers; a "No error"
             majority is positive, an error or "Incomplete" majority negative). They judged correctness, not
             equivalence to the reference.
  reach      over the main store's answered, usable rows of the eleven models: the answer check's verdict by its number
             rule (answer.verdict with the symbolic step off, not the stored label, so the counts are the same before
             and after the stores are re-scored) against the equivalence check's: the verdicts the step moves. For the
             symbolic templates without a spec, the usable rows the number rule scores incorrect or partial.
Counts only; no response text and no expert's words are written. Nothing here calls a model.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULL = HERE.parent
REPO = FULL.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402
from full_run_28092026.symbolic import equivalence as EQ  # noqa: E402

GRADES = HERE / 'grades.json'
OUT_JSON = HERE / 'validation.json'
OUT_MD = HERE / 'SYMBOLIC_CHECK.md'
BAR = 0.9
ACTING = ('incorrect', 'partial')
# Each graded item on which the check and the experts disagree, and why: the facts of the case, in our words, with
# no number of the item (the pool stays private, and so do the grades per item). Keyed by template, model and
# verdict group; a key that two disagreements share explains neither.
EXPLAINED = {
    ('template_ber_estimation_mary', 'claude-sonnet-5', 'incorrect'):
        "States a value that does not round to the reference's: it lies about 5% above it, and about 3% above the "
        'value the exact argument gives. The expert accepted it as a rounding of the steep Q-function, which the guide '
        'left to their judgement; the check holds a stated value to rounding.',
    ('template_ber_estimation_mary', 'claude-sonnet-5', 'partial'):
        'Writes the linear Eb/N0 inside the expression two units off in its last digit (a conversion of the decibel '
        "value), so its value lies 0.28% from the reference's, outside the reference's printed precision. The expert "
        'accepted the rounding.',
    ('template_phasor_addition', 'claude-sonnet-5', 'control'):
        'A magnitude one unit off in its last digit and a phase given to one decimal fewer than the reference: a '
        'rounding the expert accepted here. A control (the answer check already scores it correct) on a template the '
        'check is not applied to.',
    ('template_phasor_addition', 'gemma-4-26b-a4b', 'control'):
        'The magnitude to three decimals where the reference prints two and truncates the third, so the correct '
        "rounding lies outside the reference's printed digits. A control on a template the check is not applied to.",
}


def check_safe(item: dict, text: str, gold_unit: float | None = None) -> tuple[bool, dict]:
    try:
        return EQ.check(item, text, gold_unit=gold_unit)
    except Exception as e:  # noqa: BLE001  (a parse that SymPy cannot finish counts as not equivalent)
        return False, {'how': f'error {type(e).__name__}'}


def gold_checks(items: dict) -> dict:
    out = {}
    for tid in EQ.SPECS:
        its = sorted((i for i in items.values() if i['template_id'] == tid), key=lambda i: i['instance_index'])
        self_ok = sum(check_safe(it, it['solution'])[0] for it in its)
        pairs, matched = 0, []
        for it in its:
            for other in its:
                if other is it or other['question'] == it['question']:
                    continue
                pairs += 1
                if check_safe(it, other['solution'])[0]:
                    matched.append([it['item_id'], other['item_id']])
        out[tid] = {'items': len(its), 'gold_equivalent_to_itself': self_ok, 'other_golds_tried': pairs,
                    'other_golds_accepted': len(matched), 'accepted_pairs': matched}
    return out


def prf(tp: int, fp: int, fn: int) -> tuple[float | None, float | None]:
    return (tp / (tp + fp) if tp + fp else None), (tp / (tp + fn) if tp + fn else None)


def against(labels: list[tuple[str, bool, bool, str]]) -> dict:
    """labels: (template, check_says_equivalent, reference_says_equivalent, group) -> per template and 'all'."""
    out = {}
    by_t = collections.defaultdict(list)
    for t, c, r, g in labels:
        by_t[t].append((c, r, g))
        by_t['all'].append((c, r, g))
    for t, rows in by_t.items():
        tp = sum(c and r for c, r, _ in rows)
        fp = sum(c and not r for c, r, _ in rows)
        fn = sum(r and not c for c, r, _ in rows)
        tn = sum(not c and not r for c, r, _ in rows)
        p, rcl = prf(tp, fp, fn)
        out[t] = {'n': len(rows), 'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
                  'precision': None if p is None else round(p, 3), 'recall': None if rcl is None else round(rcl, 3)}
    return out


def store_rows_and_texts(models: list[str]) -> tuple[dict, dict]:
    rows, texts = {}, {}
    for key in models:
        store = [json.loads(l) for l in (score.SCORES / 'main' / f'{key}.jsonl').read_text(encoding='utf-8').splitlines()]
        texts[key] = score.texts_matching('main', key, store)
        rows[key] = {r['item_id']: r for r in store}
    return rows, texts


def with_grades(items: dict, texts: dict) -> dict | None:
    if not GRADES.exists():
        return None
    g = json.loads(GRADES.read_text(encoding='utf-8'))
    graded, skipped = [], collections.Counter()
    for it in g['items'].values():
        if it['template_id'] not in EQ.SPECS:
            continue
        if it['grade'] not in ('equivalent', 'not equivalent', 'unreadable'):
            skipped[it['grade'] or 'ungraded'] += 1
            continue
        text = texts[it['model']][it['item_id']]
        ok, d = check_safe(items[it['item_id']], text)
        ok1, _d1 = check_safe(items[it['item_id']], text, gold_unit=1.0)
        graded.append((it, ok, ok1, d))

    def table(rule: int, acting: bool) -> dict:
        return against([(it['template_id'], (ok, ok1)[rule], it['grade'] == 'equivalent', it['group'])
                        for it, ok, ok1, _d in graded if not acting or it['group'] in ACTING])

    disagreements = []
    differ = [(it, ok, d) for it, ok, _ok1, d in graded if ok != (it['grade'] == 'equivalent')]
    shared = collections.Counter((it['template_id'], it['model'], it['group']) for it, _ok, _d in differ)
    for it, ok, d in differ:
        key = (it['template_id'], it['model'], it['group'])
        disagreements.append({'template_id': it['template_id'], 'item_id': it['item_id'], 'model': it['model'],
                              'group': it['group'], 'grade': it['grade'],
                              'check': 'equivalent' if ok else 'not equivalent', 'detail': d.get('how'),
                              'margin': d.get('margin'),
                              'explanation': EXPLAINED.get(key) if shared[key] == 1 else None})
    acting, every = table(0, True), table(0, False)
    enabled = []
    for t, v in acting.items():
        if t == 'all' or v['precision'] is None or v['recall'] is None:
            continue
        unexplained = [x for x in disagreements if x['template_id'] == t and not x['explanation']]
        if v['precision'] >= BAR and v['recall'] >= BAR and v['tp'] + v['fp'] >= 1 and not unexplained:
            enabled.append(t)
    return {'source': 'grades.json', 'scored_at_utc': g.get('scored_at_utc'), 'skipped': dict(skipped),
            'window': {'gold_unit': EQ.GOLD_UNIT}, 'acting': acting, 'all_graded': every,
            'acting_one_unit': table(1, True), 'all_graded_one_unit': table(1, False),
            'enabled': sorted(enabled), 'disagreements': disagreements,
            'unexplained': [x for x in disagreements if not x['explanation']]}


def with_readings(items: dict) -> dict | None:
    """The B1 and B2 readings of the 2026-10-02 request that fall on these templates (local files)."""
    from full_run_28092026 import expert_kits as K
    if not (K.OUT / 'tasks' / 'keyfile.json').exists() or not K.RETURNED.exists():
        return None
    keyfile, _queues, rows = K.read_returns(K.RETURNED, K.OUT)
    shown = K.shown_pool(K.OUT)
    by_code = collections.defaultdict(list)
    for (aid, c), r in rows.items():
        by_code[c].append(r)
    labels, disagreements, skipped = [], [], collections.Counter()
    for c, k in keyfile.items():
        if k['kind'] not in ('answer', 'error') or k['template_id'] not in EQ.SPECS or not by_code.get(c):
            continue
        rs = by_code[c]
        if k['kind'] == 'answer':
            votes = collections.Counter(r['answers']['verdict'] == 'correct' for r in rs)
            group = f"B1, check {k['check']}"
        else:
            votes = collections.Counter(r['answers']['category'].startswith('No error') for r in rs)
            group = 'B2, check incorrect'
        top, n = votes.most_common(1)[0]
        if n * 2 <= sum(votes.values()):
            skipped[k['kind']] += 1
            continue
        ok, d = check_safe(items[k['item_id']], shown[c]['trace'])
        labels.append((k['template_id'], ok, top, group))
        if ok != top:
            disagreements.append({'item_id': k['item_id'], 'model': k['model'], 'reading': group,
                                  'experts': 'right' if top else 'wrong', 'check': 'equivalent' if ok else 'not equivalent',
                                  'detail': d.get('how'), 'margin': d.get('margin')})
    return {'source': 'expert_request/ (B1, B2)', 'skipped_split': dict(skipped), 'by_template': against(labels),
            'disagreements': disagreements}


def numbers_label(item: dict, text: str, ms: list[dict]) -> str:
    """The answer check's verdict by its number rule: answer.verdict as score.score_answer calls it, with the symbolic
    step off, so the reach reads the same before and after the stores are re-scored with the step."""
    enabled = answer.SYMBOLIC_EQUIVALENCE_TEMPLATES
    try:
        answer.SYMBOLIC_EQUIVALENCE_TEMPLATES = ()
        return answer.verdict(text, item, score.TOLS['fitted'], tuple(m['value'] for m in ms))[0]
    finally:
        answer.SYMBOLIC_EQUIVALENCE_TEMPLATES = enabled


def usable(rows: dict, tid: str):
    for key in rows:
        for i, r in rows[key].items():
            if r['template_id'] == tid and r['status'] == 'answered' and not r['unusable']:
                yield key, i


def reach(items: dict, rows: dict, texts: dict, ms: dict) -> dict:
    out = {}
    for tid in EQ.SPECS:
        cross = collections.Counter()
        moved_by_model = collections.Counter()
        for key, i in usable(rows, tid):
            label = numbers_label(items[i], texts[key][i], ms[i])
            ok, _d = check_safe(items[i], texts[key][i])
            cross[(label, ok)] += 1
            if ok and label != 'correct':
                moved_by_model[key] += 1
        out[tid] = {'correct': {'equivalent': cross[('correct', True)], 'not': cross[('correct', False)]},
                    'partial': {'equivalent': cross[('partial', True)], 'not': cross[('partial', False)]},
                    'incorrect': {'equivalent': cross[('incorrect', True)], 'not': cross[('incorrect', False)]},
                    'would_move_by_model': dict(sorted(moved_by_model.items()))}
    return out


def no_spec(items: dict, rows: dict, texts: dict, ms: dict, templates: list[str]) -> dict:
    """The symbolic templates without a spec: their usable responses, and how many the number rule scores incorrect or
    partial (the verdicts a spec could move)."""
    out = {}
    for tid in templates:
        labels = collections.Counter(numbers_label(items[i], texts[key][i], ms[i]) for key, i in usable(rows, tid))
        out[tid] = {'usable': sum(labels.values()), 'incorrect_or_partial': labels['incorrect'] + labels['partial']}
    return out


def fmt(x) -> str:
    return '-' if x is None else f'{x:.3f}'


def short(t: str) -> str:
    return t.replace('template_', '')


PREAMBLE = """# The symbolic equivalence check

Generated by `validate.py` from the check (`evaluator_pilot_17092026/evaluators/symbolic_equivalence.py`, which
`answer.py` calls and pins; `equivalence.py` here re-exports it), `grades.json` (the experts' E2 grades) and the main
store. Counts only.

## The rule

Nine templates answer with an expression, and the answer check scores them by the numbers they state (the numbers
the gold computed, `answer.targets`). An answer written in another but equal form then scores wrong. The check
decides equivalence instead, and the definition is the one the experts graded against: the same quantity as the
reference, to the precision the reference shows.

- **Reading.** The gold's answer line and every expression on the response's final-answer segment (`answer.segment`,
  the reader the score uses) are parsed into SymPy, from LaTeX, plain text and unicode; units and words inside an
  expression are dropped. Names the question gives a value to (Eb/N0, M, k) take it, or the answer's own Eb/N0 where
  it states a different one. A number the question gives is exact in the answer, and so is every integer, but for an
  angle in whole degrees and the mantissa of a power of ten. A sequence the answer says is zero before n = 0 is read
  with u[n].
- **Equivalence.** At every sample point of the free variable (twenty over two periods of a signal; n = -3 to 11 for
  a sequence; a grid for a velocity field; the support for the autocorrelation; one point for the bit error rate),
  the candidate lies within half a unit of the last digit of every decimal the gold writes, carried through the
  gold's expression. An answer coarser than the reference, or off in a digit the reference shows, is not equivalent;
  an exact form is (1 + sqrt 2 for 2.41, y^2/2 for 0.5y^2). For the bit error rate, whose reference c Q(a) prints no
  digits of its value, a single written number is also equivalent when the reference's value rounds to it.
- **The answer as a whole.** It is equivalent when an expression on its segment is, unless it also states the asked
  quantity as something the answer check's own number rule calls different: an answer with two values is not the
  reference. A bare 0 contradicts nothing. Where a template states a support, the answer's stated support must be the
  gold's.
- **Scope.** On an enabled template, an answer the answer check scores incorrect or partial is scored correct when it
  is equivalent. No verdict is lowered. Only the templates in `answer.SYMBOLIC_EQUIVALENCE_TEMPLATES` are checked; the
  others keep the by-numbers rule.

How the rule was settled. A first version accepted each gold number within 0.2% or a unit of its last digit (the answer
check's own number rule, carried through the expression) and any expression on the segment. Against the grades it had
recall 1.000 and precision 0.909: every one of its 22 disagreements was an answer it accepted and the experts did not.
Reading them gave the definition above (the experts held answers to the reference's printed precision, as the guide
asked) and five reader defects: an answer stating two different values; a decimal the question gives (10^2.3 for
23 dB) treated as rounded; an answer's own Eb/N0 ignored; a sequence not checked before n = 0; and parsing artefacts
(a decibel unit dropped, a percentage read as a factor, a display read as fragments). The figures below are of the
settled rule on the same grades, so they are not an estimate on unseen answers; the B1 and B2 readings further down
were made before the check existed.
"""


def report(res: dict) -> str:
    L = [PREAMBLE.rstrip(), '', '## The golds', '',
         '| template | items | gold equivalent to itself | other golds tried | accepted |', '|---|---:|---:|---:|---:|']
    for t, v in res['golds'].items():
        L.append(f"| `{short(t)}` | {v['items']} | {v['gold_equivalent_to_itself']} | {v['other_golds_tried']} | "
                 f"{v['other_golds_accepted']} |")
    acc = collections.Counter(t for t, v in res['golds'].items() for _p in v['accepted_pairs'])
    if acc:
        L += ['', 'Accepted pairs (two questions whose values coincide within the precision shown; their item ids are in '
              'validation.json, local): ' + '; '.join(f'`{short(t)}` {n}' for t, n in acc.items()) + '.']
    L += ['', '## Against the experts', '']
    g = res.get('grades')
    if g is None:
        L += ['The E2 grades have not come back; `grades.json` does not exist yet. The enabled templates are decided '
              'on them.', '']
    else:
        L += [f"From `grades.json` (scored {g['scored_at_utc']}): positive is graded equivalent; items split between "
              f"their two readers are left out ({sum(g['skipped'].values())}). \"Can move\" is the answer check's "
              'incorrect and partial verdicts, the only ones the check acts on; "all" adds the controls, verdicts the '
              'answer check scores correct.', '',
              '| template | can move: graded | TP | FP | FN | TN | precision | recall | all: graded | precision | '
              'recall | enabled |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
        for t, v in sorted(g['acting'].items(), key=lambda kv: (kv[0] == 'all', kv[0])):
            a = g['all_graded'][t]
            en = '-' if t == 'all' else ('yes' if t in g['enabled'] else 'no')
            L.append(f"| {'all' if t == 'all' else '`' + short(t) + '`'} | {v['n']} | {v['tp']} | {v['fp']} | {v['fn']} | "
                     f"{v['tn']} | {fmt(v['precision'])} | {fmt(v['recall'])} | {a['n']} | {fmt(a['precision'])} | "
                     f"{fmt(a['recall'])} | {en} |")
        L += ['', 'A template is enabled when, on the verdicts it can move, precision and recall are both at least 0.9, '
              'it moves at least one graded verdict, and every disagreement on it is explained. Phasor addition and the '
              'undamped oscillator move no verdict the experts grade equivalent, so they stay by the numbers.', '',
              '**Sensitivity: a whole unit of the gold\'s last digits** (the answer check\'s own last-digit window), on '
              'the verdicts the check can move: ' + '; '.join(
                  f"{'all' if t == 'all' else '`' + short(t) + '`'} {fmt(v['precision'])} / {fmt(v['recall'])}"
                  for t, v in sorted(g['acting_one_unit'].items(), key=lambda kv: (kv[0] == 'all', kv[0])))
              + ' (precision / recall). The wider window accepts values that do not round to the reference, which '
              'the experts did not.', '',
              f"**Disagreements under the rule: {len(g['disagreements'])}**, {len(g['unexplained'])} unexplained.", '']
        for x in g['disagreements']:                  # the template, not the item: the per-item grades stay local
            L.append(f"- `{short(x['template_id'])}`, `{x['model']}` ({x['group']}; graded {x['grade']}, the check "
                     f"{x['check']}): {x['explanation'] or 'UNEXPLAINED'}")
        L.append('')
    r = res.get('readings')
    L += ['## The B1 and B2 readings on these templates', '']
    if r is None:
        L += ['The readings are not on this machine.', '']
    else:
        L += ['Made before the check existed, they judged correctness rather than equivalence to the reference: B1 (two '
              'readers, the final answer correct or not) and B2 (three readers, "No error" or an error). A split is left '
              f"out ({sum(r['skipped_split'].values())}).", '',
              '| template | readings | precision | recall | TP | FP | FN | TN |', '|---|---:|---:|---:|---:|---:|---:|---:|']
        for t, v in sorted(r['by_template'].items(), key=lambda kv: (kv[0] == 'all', kv[0])):
            L.append(f"| {'all' if t == 'all' else '`' + short(t) + '`'} | {v['n']} | {fmt(v['precision'])} | "
                     f"{fmt(v['recall'])} | {v['tp']} | {v['fp']} | {v['fn']} | {v['tn']} |")
        L += ['', f"Disagreements: {len(r['disagreements'])}.", '']
    L += ['## Reach over the main store', '',
          'The eleven models\' answered, usable responses: the answer check\'s verdict by its number rule (the symbolic '
          'step off, as `answer.verdict` gives it today, so the table reads the same before and after the stores are '
          're-scored) against the equivalence check. An enabled template scores the incorrect and partial verdicts the '
          'check reads as equivalent correct; a correct verdict the check does not recognise stays correct.', '',
          '| template | enabled | correct: equivalent / not | partial: equivalent / not | incorrect: equivalent / not |',
          '|---|---|---|---|---|']
    enabled = set((g or {}).get('enabled', []))
    moved = collections.Counter()
    for t, v in res['reach'].items():
        L.append(f"| `{short(t)}` | {'yes' if t in enabled else 'no'} | {v['correct']['equivalent']} / {v['correct']['not']} | "
                 f"{v['partial']['equivalent']} / {v['partial']['not']} | {v['incorrect']['equivalent']} / "
                 f"{v['incorrect']['not']} |")
        if t in enabled:
            moved.update(v['would_move_by_model'])
    L += ['', f"Verdicts the enabled templates move to correct: {sum(moved.values())}; by model: "
          + (', '.join(f'`{m}` {n}' for m, n in sorted(moved.items(), key=lambda kv: (-kv[1], kv[0]))) or 'none') + '.']
    L += ['', '## Templates that stay by the numbers', '']
    L += [f"- `{short(t)}`: not enabled (see the table above)." for t in EQ.SPECS if t not in enabled]
    L += [f"- `{short(t)}`: no spec; the number rule scores {v['incorrect_or_partial'] or 'none'} of its {v['usable']} "
          'usable responses incorrect or partial.' for t, v in res['no_spec'].items()]
    L.append('')
    return '\n'.join(L)


def run(write: bool = True) -> dict:
    items = score.pool_items()
    ms = score.milestone_sets(items)
    rows, texts = store_rows_and_texts(list(ROSTER))
    sym = sorted({i['template_id'] for i in items.values() if i['answer_type'] == 'symbolic'})
    res = {'generated_by': 'full_run_28092026/symbolic/validate.py',
           'run_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
           'golds': gold_checks(items), 'grades': with_grades(items, texts), 'readings': with_readings(items),
           'reach': reach(items, rows, texts, ms),
           'no_spec': no_spec(items, rows, texts, ms, [t for t in sym if t not in EQ.SPECS])}
    if write:
        OUT_JSON.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8', newline='\n')
        OUT_MD.write_text(report(res), encoding='utf-8', newline='\n')
    return res


def selftest() -> int:
    rng = random.Random(1)
    labels = [('a', True, True, 'incorrect')] * 18 + [('a', True, False, 'partial')] + [('a', False, True, 'partial')] \
        + [('a', False, False, 'incorrect')] * 5 + [('b', True, True, 'partial')] * 3 + [('b', False, True, 'partial')] * 2
    rng.shuffle(labels)
    r = against(labels)
    assert r['a']['precision'] == round(18 / 19, 3) and r['a']['recall'] == round(18 / 19, 3), r['a']
    assert r['b']['precision'] == 1.0 and r['b']['recall'] == 0.6, r['b']
    assert r['all']['n'] == len(labels)
    assert prf(0, 0, 0) == (None, None)
    print('SELFTEST OK: precision and recall per template on synthetic labels')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    res = run()
    print(report(res))
    return 0


if __name__ == '__main__':
    sys.exit(main())
