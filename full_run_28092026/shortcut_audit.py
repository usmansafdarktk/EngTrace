"""The corpus-wide surface-shortcut audit: D-057's measurement on every template (docs/EVALUATION_NEXT_STEPS.md A7; D-177).

    python -m full_run_28092026.shortcut_audit [--seeds 500]      # FREE: SHORTCUT_AUDIT.md, results/shortcut_audit.json

TWO MEASUREMENTS, both of the question's surface and neither of what any model did (the paraphrase arm, Q5, is the
behavioural half).

  (a) Classification templates, on public draws. D-057 found two templates whose label a depth-2 decision tree reads
      off the question's wording with 100% held-out accuracy against blind-guess floors of 50% and 34%. The same
      measurement here on every classification template in the pool: `--seeds` draws of the public template code
      (seeds 0 to N-1; the pool's private seed is never used), the label read from the gold's Answer segment by the
      answer check's own word and labelled-part readers, the question's numbers masked and its remaining word and
      symbol tokens a binary bag, a depth-2 tree (sklearn) fitted on the even seeds and scored on the odd ones; the
      blind-guess floor is the share of the odd seeds that carry the even seeds' majority label, and the lift is the
      held-out accuracy minus that floor. A lift near D-057's two (+50 and +66 points) marks a lookup; a small lift
      means the label needs the numbers.
  (b) Every template, on the frozen pool's items as the store scores them: how often the check's numeric targets are
      numbers the question already states. A target counts as restated when a question number lies within the
      check's relative tolerance of it, and as rescaled when the ratio to a question number is a power of ten
      (within the same tolerance). Per template, the share of items with every target restated or rescaled and the
      share with at least one. An answer that restates an input is not a reasoning item; one that rescales an input
      is a unit conversion. Items whose targets are words only (classification) have no numeric target and are
      counted apart.

WHAT IS FLAGGED. A classification template with a lift of 0.25 or more; a template with every target restated or
rescaled on half or more of its items. The four templates the pilot excluded for shortcuts (D-057, D-046, D-066) are
marked wherever they appear; the Sensitivity table of RESULTS.md already gives the headline without them. The audit
ends with the answer score per model without any newly flagged template, from the store, beside the headline.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER, SHORTCUT  # noqa: E402
from full_run_28092026.diversity import NUM  # noqa: E402
from full_run_28092026.freeze import answer_types  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402

OUT_MD = HERE / 'SHORTCUT_AUDIT.md'
OUT_JSON = HERE / 'results' / 'shortcut_audit.json'
LIFT_FLAG = 0.25
RESTATED_FLAG = 0.5
TOKEN = r"(?u)\b\w[\w\-]*\b|[\^\*\+\-/=\(\)\[\]<>]"


def tid_of(ref) -> str:
    t = ref.template_id
    return t if t.startswith('template_') else 'template_' + t


FROUDE = 'template_critical_depth_froude_classification'


def label_of(solution: str, tid: str) -> tuple:
    seg = answer.segment(solution)
    if tid == FROUDE:                                   # D-118: the answer is Fr, the regime follows from it
        vals = [v for v, _u in answer.values(seg)]
        if not vals:
            return ()
        fr = vals[0]
        return ('supercritical' if fr > 1 else 'subcritical' if fr < 1 else 'critical',)
    words = tuple(answer._words(seg, 'classification'))
    pairs = tuple(f'{a}:{b}' for a, b in answer._labeled(seg, 'classification'))
    return words + pairs


def surface(question: str) -> str:
    return NUM.sub(' <num> ', question)


def classification_audit(seeds: int, types: dict) -> list[dict]:
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.tree import DecisionTreeClassifier
    out = []
    for ref in discover(None):
        tid = tid_of(ref)
        if types.get(tid) != 'classification':
            continue
        draws = [generate(ref, s, capture=False) for s in range(seeds)]
        draws = [d for d in draws if d.ok]
        labelled = [(d.seed, surface(d.question), label_of(d.solution, tid)) for d in draws]
        unlabelled = sum(1 for _, _, lab in labelled if not lab)
        labelled = [x for x in labelled if x[2]]
        train = [x for x in labelled if x[0] % 2 == 0]
        test = [x for x in labelled if x[0] % 2 == 1]
        if len(train) < 10 or len(test) < 10:
            out.append({'template': tid, 'draws': len(draws), 'unlabelled': unlabelled, 'error': 'too few labelled draws'})
            continue
        vec = CountVectorizer(binary=True, lowercase=True, token_pattern=TOKEN, min_df=2)
        Xtr = vec.fit_transform([q for _, q, _ in train])
        Xte = vec.transform([q for _, q, _ in test])
        ytr = ['|'.join(lab) for _, _, lab in train]
        yte = ['|'.join(lab) for _, _, lab in test]
        majority = collections.Counter(ytr).most_common(1)[0][0]
        floor = sum(y == majority for y in yte) / len(yte)
        tree = DecisionTreeClassifier(max_depth=2, random_state=7).fit(Xtr, ytr)
        acc = float(tree.score(Xte, yte))
        names = vec.get_feature_names_out()
        used = sorted({str(names[f]) for f in tree.tree_.feature if f >= 0})
        out.append({'template': tid, 'draws': len(draws), 'unlabelled': unlabelled, 'train': len(train), 'test': len(test),
                    'labels': len(set(ytr) | set(yte)), 'label_shares_train': dict(collections.Counter(ytr)),
                    'floor': float(floor), 'held_out_accuracy': acc, 'lift': acc - float(floor), 'splitting_tokens': used,
                    'known_shortcut': tid in SHORTCUT, 'flagged': acc - float(floor) >= LIFT_FLAG})
    return out


def power_of_ten_ratio(g: float, v: float, tol: float) -> bool:
    if g == 0 or v == 0:
        return False
    r = abs(g / v)
    k = round(math.log10(r))
    return k != 0 and abs(r - 10 ** k) <= tol * 10 ** k


def restated_audit(types: dict) -> dict:
    items = score.pool_items()
    rows = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{ROSTER[0]}.jsonl')
                                            .read_text(encoding='utf-8').splitlines())}
    per_t = collections.defaultdict(lambda: collections.Counter())
    for iid, r in rows.items():
        t = r['template_id']
        targets = r['answer']['targets']['numbers'] if r['answer'] and r['answer'].get('targets') else []
        per_t[t]['items'] += 1
        if not targets:
            per_t[t]['no_numeric_target'] += 1
            continue
        qvals = [v for v, _u in answer.values(items[iid]['question'])]
        hits = []
        for g in targets:
            restated = any(abs(v - g) <= answer.REL * abs(g) for v in qvals)
            rescaled = (not restated) and any(power_of_ten_ratio(g, v, answer.REL) for v in qvals)
            hits.append('restated' if restated else 'rescaled' if rescaled else None)
        per_t[t]['targets'] += len(hits)
        per_t[t]['targets_restated'] += sum(h == 'restated' for h in hits)
        per_t[t]['targets_rescaled'] += sum(h == 'rescaled' for h in hits)
        per_t[t]['items_all'] += all(h is not None for h in hits)
        per_t[t]['items_any'] += any(h is not None for h in hits)
    out = {}
    for t, c in per_t.items():
        n = c['items'] - c['no_numeric_target']
        out[t] = {'answer_type': types.get(t), 'items': c['items'], 'items_with_numeric_targets': n,
                  'targets': c['targets'], 'targets_restated': c['targets_restated'], 'targets_rescaled': c['targets_rescaled'],
                  'share_items_all_restated_or_rescaled': c['items_all'] / n if n else None,
                  'share_items_any_restated_or_rescaled': c['items_any'] / n if n else None,
                  'known_shortcut': t in SHORTCUT,
                  'flagged': bool(n and c['items_all'] / n >= RESTATED_FLAG)}
    return out


def without(flagged: set[str]) -> list[dict]:
    """Each model's answer score as the mean of template means, with and without the flagged templates."""
    out = []
    for k in ROSTER:
        rows = map(json.loads, (score.SCORES / 'main' / f'{k}.jsonl').read_text(encoding='utf-8').splitlines())
        by_t = collections.defaultdict(list)
        for r in rows:
            by_t[r['template_id']].append(r['score'])
        means = {t: float(np.mean(v)) for t, v in by_t.items()}
        out.append({'model': k, 'headline': float(np.mean(list(means.values()))),
                    'without_flagged': float(np.mean([m for t, m in means.items() if t not in flagged])) if flagged else None})
    return out


def render(res: dict) -> str:
    L = ['# The surface-shortcut audit, corpus-wide', '',
         'Generated by `shortcut_audit.py`; the two measurements and the flags are defined in its docstring. Part (a) is on '
         f"{res['seeds']} public draws per classification template (the pool's seed is never used); part (b) is on the frozen "
         'pool\'s items as the main store scores them. Counts and shares only (D-177; next steps A7).', '',
         '## (a) Classification templates: can a depth-2 tree read the label off the question\'s wording?', '',
         '| template | draws (unlabelled) | labels | blind-guess floor | depth-2 held-out accuracy | lift | splitting tokens | known shortcut | flagged |',
         '|---|---:|---:|---:|---:|---:|---|---|---|']
    for r in res['classification']:
        if r.get('error'):
            L.append(f"| `{r['template'].removeprefix('template_')}` | {r['draws']} ({r['unlabelled']}) | | | | | {r['error']} | | |")
            continue
        L.append(f"| `{r['template'].removeprefix('template_')}` | {r['draws']} ({r['unlabelled']}) | {r['labels']} | {r['floor']:.3f} | "
                 f"{r['held_out_accuracy']:.3f} | {r['lift']:+.3f} | {', '.join(f'`{t}`' for t in r['splitting_tokens']) or '-'} | "
                 f"{'yes' if r['known_shortcut'] else ''} | {'yes' if r['flagged'] else ''} |")
    rs = res['restated']
    flagged_b = sorted(t for t, v in rs.items() if v['flagged'])
    any_b = sorted((t for t, v in rs.items() if v['items_with_numeric_targets'] and v['share_items_any_restated_or_rescaled']),
                   key=lambda t: -rs[t]['share_items_all_restated_or_rescaled'])
    n_num = sum(1 for v in rs.values() if v['items_with_numeric_targets'])
    L += ['', '## (b) Every template: are the answer\'s numbers already in the question?', '',
          f"{n_num} of {len(rs)} templates have numeric targets. Templates with any item whose target is a question number "
          f"restated or rescaled by a power of ten: {len(any_b)}. Flagged (every target restated or rescaled on at least "
          f"{RESTATED_FLAG:.0%} of items): {len(flagged_b)}.", '',
          '| template | answer type | items with numeric targets | targets | restated | rescaled | items, all targets | items, any target | known shortcut | flagged |',
          '|---|---|---:|---:|---:|---:|---:|---:|---|---|']
    for t in any_b:
        v = rs[t]
        L.append(f"| `{t.removeprefix('template_')}` | {v['answer_type']} | {v['items_with_numeric_targets']} | {v['targets']} | "
                 f"{v['targets_restated']} | {v['targets_rescaled']} | {v['share_items_all_restated_or_rescaled']:.2f} | "
                 f"{v['share_items_any_restated_or_rescaled']:.2f} | {'yes' if v['known_shortcut'] else ''} | {'yes' if v['flagged'] else ''} |")
    newly = sorted(set(flagged_b) | {r['template'] for r in res['classification'] if r.get('flagged')})
    newly = [t for t in newly if t not in SHORTCUT]
    L += ['', '## Flagged templates, and the headline without them', '',
          'Known shortcut templates (already a Sensitivity row): ' + ', '.join(f'`{t.removeprefix("template_")}`' for t in SHORTCUT) + '. '
          + ('Newly flagged: ' + ', '.join(f'`{t.removeprefix("template_")}`' for t in newly) + '.' if newly else 'No template is newly flagged.'), '']
    if res['without']:
        L += ['| model | headline | without the newly flagged |', '|---|---:|---:|']
        for r in res['without']:
            L.append(f"| `{r['model']}` | {r['headline']:.3f} | {r['without_flagged']:.3f} |")
    L += ['', 'What this measures: predictability of the label, or the presence of the answer, in the question\'s surface. It does '
          'not say what the models did; the paraphrase arm (Q5) is that measurement. A rescaled target is a unit conversion the '
          'item asks for, not necessarily a shortcut; the flag is a prompt to read the template.', '']
    return '\n'.join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=500)
    a = ap.parse_args()
    types = answer_types()
    types = {(t if t.startswith('template_') else 'template_' + t): v for t, v in types.items()}
    res = {'seeds': a.seeds, 'classification': classification_audit(a.seeds, types), 'restated': restated_audit(types)}
    newly = {t for t, v in res['restated'].items() if v['flagged']} | {r['template'] for r in res['classification'] if r.get('flagged')}
    newly -= set(SHORTCUT)
    res['newly_flagged'] = sorted(newly)
    res['without'] = without(newly) if newly else []
    OUT_JSON.parent.mkdir(exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1), encoding='utf-8')
    text = render(res)
    OUT_MD.write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
