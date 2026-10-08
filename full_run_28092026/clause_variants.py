"""Scoring-rule variants of the final-answer check, recomputed offline (WS-C3 step 1; review 1 M3 and Q4,
review 2 W8, W10, Q5 and Q13).

    python -m full_run_28092026.clause_variants [--quick]   # writes results/sensitivity_variants.json and
                                                           # results/sections/sensitivity_variants.md
    python -m full_run_28092026.clause_variants --selftest  # the variants on hand-made cases; writes nothing

WHAT IT ASKS. How many verdicts rest on each clause of the final-answer rule, and on the digits the
prescribing templates ask for. Each variant changes one clause and leaves the rest of the rule as it is:
  abs_clause_off       a stated value matches only with its own sign: the |v| reading (a deflection stated
                       negative against a gold stated positive) is dropped, for the target and its renderings
  last_digit_bounded   the window of one unit in the response's last displayed digit, u(y^), is capped at 1%
                       of the gold, min(u(y^), 0.01|y|), so a coarse rounding (`3.4 x 10^3` for 3,364) no
                       longer vouches for itself; the relative tolerance and the gold's own last-digit window stay
  prescribed_relaxed   the questions that prescribe their rounding (answer.ROUNDING, 17 templates) are scored at
                       the tolerance like every other question, instead of on the exact digits
  per_part_credit      a partial answer scores matched / of instead of 0.5; correct 1 and incorrect 0 unchanged
An unusable row (no readable answer) scores 0 under every variant, as in the headline.

HOW. answer.py has no parameter for the first three clauses (its `tol`, `unit` and `whole` are other
readings), so this file holds a local copy of answer.verdict and answer.match with three switches; every
other piece (targets, segment, values, the word and label rules, the tolerance, the unit factors) is
answer.py's own, imported. The copy with every switch off must reproduce the stored label, `matched` and
`of` of every readable row: the run stops if one row differs, so the copy cannot drift from the rule the
store was scored with. per_part_credit needs no re-read: a partial answer's stored `matched` and `of`.

THE SYMBOLIC CHECK (WS-E, at ANSWER FINAL). Once answer.py carries it, verdict() raises an incorrect or partial
label to correct on the enabled templates when the final answer is equivalent to the gold
(symbolic_equivalence.apply, after the rule's own label; no verdict lowered). This file reads the module and the
enabled list from answer.py (`answer.symbolic_equivalence`, `answer.SYMBOLIC_EQUIVALENCE_TEMPLATES`, the names
WS-E's report gives) and applies the same step after each variant's label, so every variant is the full rule with
one clause changed; today neither name exists and nothing is applied. A row the step raised keeps the by-numbers
`matched` and `of` out of the reproduction check, scores 1 under per-part credit, and is counted apart among the
accepted answers (`symbolic`: no stated number was read). If answer.py adopts other names, the reproduction check
stops the run on the first raised row; adapt symbolic_hook().

WHAT IT REPORTS, per store (main: the eleven roster models; each reasoning store named in
results/matched_config.json, a sibling with `run: false` skipped; and `matched`, each model at its reasoning
store where it has one, else at main): FAC under each variant (the item mean: correct 1, partial 0.5 or
matched / of, incorrect and unusable 0) with its template bootstrap interval (percentile, templates
resampled), the same interval for the variant minus the headline, the number of verdicts each variant moves up
and down, and Kendall's tau between the variant's ordering of the models and the headline's. Among accepted
answers (readable, headline label correct): the relative error of the stated value the rule accepted for each
numeric target (the closest accepted one, under the unit factors, the sign and the renderings the rule allows),
binned by the answer's worst numeric part: <= 0.2% (the tolerance), 0.2 to 1%, 1 to 5%, over 5%. Answers whose
targets are words or labels only are counted apart (`no_numeric`).

Reads the score stores, the traces (every answered row's recorded hash is checked, score.texts_matching) and
the milestone cache (used only when its sidecar names the current manifest and milestones.py; never rebuilt
here). Calls no model. `--quick` (or ENGTRACE_QUICK=1) draws 1,000 resamples instead of 10,000.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import datetime as dt
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

RESULTS = HERE / 'results'
SECTIONS = RESULTS / 'sections'
OUT_JSON = RESULTS / 'sensitivity_variants.json'
OUT_MD = SECTIONS / 'sensitivity_variants.md'
VARIANTS = ('abs_clause_off', 'last_digit_bounded', 'prescribed_relaxed', 'per_part_credit')
LABELS = {'headline': 'headline', 'abs_clause_off': 'absolute-value clause off',
          'last_digit_bounded': 'last digit bounded', 'prescribed_relaxed': 'prescribed digits relaxed',
          'per_part_credit': 'per-part credit'}
SCORE_OF = {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}
BOUND = 0.01                     # the cap on the response's last-digit window, as a fraction of the gold
BINS = (('le_0.2', answer.REL), ('0.2_1', 0.01), ('1_5', 0.05), ('gt_5', math.inf))


def quick_mode(flag: bool) -> bool:
    return bool(flag or os.environ.get('ENGTRACE_QUICK') == '1')


# ------------------------------------------------------------------ the stores

def matched_config() -> dict:
    """results/matched_config.json, read fresh: model -> its reasoning store (None where it has none). An
    entry with `run: false` or no default store (the sibling that was not run) is skipped."""
    cfg = json.loads((RESULTS / 'matched_config.json').read_text(encoding='utf-8'))
    return {m['model']: m.get('reasoning_store') for m in cfg['models']
            if m.get('run', True) and m.get('default_store')}


def stores() -> dict[str, list[str]]:
    """Each store with the models it holds for this analysis: main's roster, then each reasoning store."""
    out = {'main': [k for k in ROSTER if (score.SCORES / 'main' / f'{k}.jsonl').exists()]}
    for k, s in matched_config().items():
        if s and (score.SCORES / s / f'{k}.jsonl').exists():
            out.setdefault(s, []).append(k)
    return out


def load_rows(store: str, key: str) -> list[dict]:
    p = score.SCORES / store / f'{key}.jsonl'
    return [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines()]


def milestone_cache(items: dict) -> dict[str, list[dict]]:
    """The milestone cache score.py wrote, trusted only when its sidecar names the current manifest and
    milestones.py, as score.milestone_sets trusts it; unlike that function this one never rebuilds it."""
    meta = json.loads(score.MILESTONES_META.read_text(encoding='utf-8'))
    want = {'manifest_sha256': score.sha((HERE / 'manifest.jsonl').read_bytes()),
            'milestones_py_sha256_lf': score.norm_sha(score.EVAL / 'milestones.py'), 'items': len(items)}
    if {k: meta.get(k) for k in want} != want:
        raise SystemExit('scores/milestones.json was built from another manifest or milestones.py; '
                         're-score the stores before running this analysis')
    cache = json.loads(score.MILESTONES.read_text(encoding='utf-8'))
    if set(cache) != set(items):
        raise SystemExit('scores/milestones.json does not hold the pool items')
    return cache


# ------------------------------------------------------------------ the rule, with its clauses switchable

def _match(gold, have, rel, exact, gold_ulp, scales=answer.SCALES, bounded=False, errs=None):
    """answer.match at unit 1.0 (the validated rule), with the response's last-digit window optionally capped
    at BOUND x |gold|. With `errs`, every accepted (value, unit factor) adds its relative error, and the
    search does not stop at the first."""
    hit = False
    for v, u in have:
        for sc in scales:
            t, tu = v * sc, u * sc
            if exact:
                ok = abs(t - gold) <= 1e-9 * max(1.0, abs(gold))
            else:
                own = tu if abs(v) > u else 0.0
                if bounded:
                    own = min(own, BOUND * abs(gold))
                ok = abs(t - gold) <= max(rel * abs(gold), own, gold_ulp) * (1 + answer.SLACK)
            if ok:
                if errs is None:
                    return True
                hit = True
                errs.append(abs(t - gold) / abs(gold) if gold else (0.0 if t == gold else math.inf))
    return hit


def _hits(want, have, low, abs_clause=True, bounded=False, errs=None):
    """answer.verdict's per-target hits at the fitted tolerance (unit 1.0, whole False), the absolute-value
    clause and the last-digit cap switchable. With `errs` (a list per numeric target), every accepted reading
    of each target adds its relative error."""
    absolute = [(abs(v), u) for v, u in have]
    rel, exact = answer.REL, want['exact']
    hits = []
    for k, (v, gu) in enumerate(want['numbers']):
        e = None if errs is None else errs[k]
        parts = [_match(v, have, rel, exact, gu, bounded=bounded, errs=e)]
        if abs_clause:
            parts.append(_match(abs(v), absolute, rel, exact, gu, bounded=bounded, errs=e))
        for w, wu in want.get('renderings', []):
            parts.append(_match(w, have, rel, exact, wu, scales=(1.0,), bounded=bounded, errs=e))
            if abs_clause:
                parts.append(_match(abs(w), absolute, rel, exact, wu, scales=(1.0,), bounded=bounded, errs=e))
            if errs is None and any(parts):
                break
        hits.append(any(parts))
    for w in want['words']:
        hits.append(answer._word_hit(low, w))
    for label, value in want.get('labeled', []):
        hits.append(answer._stance(low, label) == value)
    return hits


def _label(hits) -> tuple[str, int, int]:
    if not hits:
        return 'incorrect', 0, 0
    n = sum(hits)
    return ('correct' if n == len(hits) else ('partial' if n else 'incorrect')), n, len(hits)


def symbolic_hook():
    """answer.py's symbolic equivalence step, once WS-E enables it: (label, item, text) -> label. None while
    answer.py does not carry it (as on 2026-10-07)."""
    mod = getattr(answer, 'symbolic_equivalence', None)
    enabled = getattr(answer, 'SYMBOLIC_EQUIVALENCE_TEMPLATES', None)
    if mod is None or not enabled:
        return None
    return lambda label, item, text: mod.apply(label, item, text, enabled, answer.REL, 1.0, answer.segment)[0]


def variant_verdicts(text: str, item: dict, ms, hook=None) -> dict:
    """The headline verdict and the three re-read variants of one response, with the relative errors of an
    accepted answer's numeric targets. With `hook` (symbolic_hook), each label then goes through the symbolic step,
    which depends on the response alone, so it is asked once."""
    want = answer.targets(item, ms)
    seg = answer.segment(text)
    have = answer.values(seg)
    low = answer._norm_linear(seg).lower()
    out = {'headline': _label(_hits(want, have, low)),
           'abs_clause_off': _label(_hits(want, have, low, abs_clause=False)),
           'last_digit_bounded': _label(_hits(want, have, low, bounded=True)),
           'prescribed_relaxed': _label(_hits(dict(want, exact=False), have, low)),
           'exact': want['exact'], 'rel_errors': None, 'raised': False}
    if out['headline'][0] == 'correct' and want['numbers']:
        errs = [[] for _ in want['numbers']]
        _hits(want, have, low, errs=errs)
        out['rel_errors'] = [min(e) if e else math.inf for e in errs]
    out['by_numbers'] = out['headline']
    if hook is not None and any(out[v][0] != 'correct' for v in ('headline',) + VARIANTS[:3]):
        if hook('incorrect', item, text) == 'correct':
            out['raised'] = out['headline'][0] != 'correct'
            for v in ('headline',) + VARIANTS[:3]:
                out[v] = ('correct',) + tuple(out[v][1:])
    return out


def bin_of(err: float) -> str:
    return next(name for name, hi in BINS if err <= hi * (1 + answer.SLACK))


# ------------------------------------------------------------------ one model in one store

def score_model(args) -> dict:
    """Every row of one model in one store under the headline and each variant (a worker's unit)."""
    store, key, items, ms = args
    rows = load_rows(store, key)
    texts = score.texts_matching(store, key, rows)
    hook = symbolic_hook()
    out, bad, raised = [], [], 0
    for r in rows:
        rec = {'template_id': r['template_id'], 'answer_type': r['answer_type'], 'headline': r['score']}
        if r['unusable']:
            rec.update({v: 0.0 for v in VARIANTS})
            rec.update({'bin': None})
            out.append(rec)
            continue
        a = r['answer']
        vals = tuple(m['value'] for m in ms[r['item_id']])          # what score.score_answer passes
        vv = variant_verdicts(texts[r['item_id']], items[r['item_id']], vals, hook)
        raised += vv['raised']
        # the label always; matched and of where the symbolic step did not raise the label
        same = vv['headline'][0] == a['label'] and (vv['raised'] or vv['headline'][1:] == (a['matched'], a['of']))
        if not same:
            bad.append((r['item_id'], vv['headline'], (a['label'], a['matched'], a['of'])))
        for v in ('abs_clause_off', 'last_digit_bounded', 'prescribed_relaxed'):
            rec[v] = SCORE_OF[vv[v][0]]
        rec['per_part_credit'] = {'correct': 1.0, 'incorrect': 0.0}.get(
            a['label'], (a['matched'] / a['of']) if a['of'] else 0.0)
        if vv['raised']:
            rec['bin'] = 'symbolic'
        elif vv['rel_errors'] is None:
            rec['bin'] = 'no_numeric' if a['label'] == 'correct' else None
        else:
            rec['bin'] = bin_of(max(vv['rel_errors']))
        out.append(rec)
    return {'store': store, 'key': key, 'rows': out, 'mismatch': bad, 'raised': raised,
            'readable': sum(not r['unusable'] for r in rows)}


# ------------------------------------------------------------------ statistics

def template_boot(by_t: dict[str, list[float]], draws: int, seed: int) -> list[float]:
    """Percentile interval of the item mean (equal-size templates: the mean of template means), templates
    resampled."""
    means = np.array([np.mean(v) for v in by_t.values()])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(means), size=(draws, len(means)))
    d = means[idx].mean(axis=1)
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def summarise(res: dict, draws: int, seed: int) -> dict:
    rows = res['rows']
    by_t = {c: collections.defaultdict(list) for c in ('headline',) + VARIANTS}
    diff = {v: collections.defaultdict(list) for v in VARIANTS}
    changed = {v: {'up': 0, 'down': 0} for v in VARIANTS}
    where = {v: collections.Counter() for v in VARIANTS}
    for r in rows:
        for c in ('headline',) + VARIANTS:
            by_t[c][r['template_id']].append(r[c])
        for v in VARIANTS:
            diff[v][r['template_id']].append(r[v] - r['headline'])
            if r[v] > r['headline'] + 1e-12:
                changed[v]['up'] += 1
                where[v][r['template_id']] += 1
            elif r[v] < r['headline'] - 1e-12:
                changed[v]['down'] += 1
                where[v][r['template_id']] += 1
    m = {c: float(np.mean([r[c] for r in rows])) for c in ('headline',) + VARIANTS}
    m['ci'] = {c: template_boot(by_t[c], draws, seed + k) for k, c in enumerate(('headline',) + VARIANTS)}
    m['delta_ci'] = {v: template_boot(diff[v], draws, seed + 10 + k) for k, v in enumerate(VARIANTS)}
    m['changed'] = changed
    m['changed_templates'] = {v: dict(where[v]) for v in VARIANTS}
    m['rows'], m['readable'], m['raised_by_symbolic'] = len(rows), res['readable'], res.get('raised', 0)
    return m


def where_changed(models: dict, top: int = 5) -> dict:
    """Per variant, the templates holding most of the changed verdicts over a store's models, and how many
    templates hold any."""
    out = {}
    for v in VARIANTS:
        c = collections.Counter()
        for m in models.values():
            c.update(m['changed_templates'][v])
        out[v] = {'templates': len(c), 'verdicts': sum(c.values()), 'top': dict(c.most_common(top))}
    return out


def bins_of(rows) -> dict:
    c = collections.Counter(r['bin'] for r in rows if r['bin'] is not None)
    numeric = sum(c[n] for n, _ in BINS)
    out = {n: c[n] for n, _ in BINS}
    out.update({'accepted': numeric + c['no_numeric'] + c['symbolic'], 'numeric': numeric,
                'no_numeric': c['no_numeric'], 'symbolic': c['symbolic'],
                'share': {n: (c[n] / numeric if numeric else None) for n, _ in BINS}})
    # where an accepted error above 1% comes from: its templates and answer kinds, and how many the
    # last-digit cap still accepts (through the gold's own last digit or another reading)
    over = [r for r in rows if r['bin'] in ('1_5', 'gt_5')]
    out['over_1pct'] = {'n': len(over), 'kept_when_bounded': sum(r['last_digit_bounded'] == 1.0 for r in over),
                        'by_answer_type': dict(collections.Counter(r['answer_type'] for r in over).most_common()),
                        'top_templates': dict(collections.Counter(r['template_id'] for r in over).most_common(5))}
    return out


def tau_block(models: dict) -> dict:
    keys = list(models)
    if len(keys) < 2:
        return {v: None for v in VARIANTS}
    head = [models[k]['headline'] for k in keys]
    return {v: float(stats.kendalltau(head, [models[k][v] for k in keys]).statistic) for v in VARIANTS}


# ------------------------------------------------------------------ the run

def run(quick: bool, workers: int) -> dict:
    draws = 1_000 if quick else 10_000
    if quick:
        print('QUICK: 1,000 resamples')
    items = score.pool_items()
    ms = milestone_cache(items)
    plan = stores()
    jobs = [(s, k, items, ms) for s, ks in plan.items() for k in ks]
    done = {}
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for res in pool.map(score_model, jobs):
            if res['mismatch']:
                raise SystemExit(f"{res['store']}/{res['key']}: the local copy of the rule differs from the stored "
                                 f"verdict on {len(res['mismatch'])} rows (first: {res['mismatch'][:3]}); the store "
                                 f"was scored with another answer.py, or the copy has drifted")
            done[(res['store'], res['key'])] = res
            print(f"{res['store']}/{res['key']}: {len(res['rows'])} rows, every readable verdict reproduced", flush=True)
    out = {'quick': quick, 'generated_by': 'python -m full_run_28092026.clause_variants' + (' --quick' if quick else ''),
           'written_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
           'resamples': draws, 'prescribing_templates': sorted({it['template_id'] for it in items.values()
                                                                if answer.ROUNDING.search(it['question'] or '')}),
           'symbolic_check_templates': sorted(getattr(answer, 'SYMBOLIC_EQUIVALENCE_TEMPLATES', None) or [])
           if symbolic_hook() else [],
           'stores': {}}
    seed = 3100
    for s, ks in plan.items():
        models = {}
        for k in ks:
            seed += 100
            models[k] = summarise(done[(s, k)], draws, seed)
        out['stores'][s] = {'models': models, 'tau_with_headline': tau_block(models),
                            'where_changed': where_changed(models),
                            'relative_error_bins': {k: bins_of(done[(s, k)]['rows']) for k in ks}}
        out['stores'][s]['relative_error_bins']['all'] = bins_of([r for k in ks for r in done[(s, k)]['rows']])
    # the matched configuration: each model at its reasoning store where it has one, else at main
    cfg = matched_config()
    mm, bins = {}, {}
    for k in plan['main']:
        s = cfg.get(k) if cfg.get(k) in plan and k in plan[cfg.get(k)] else 'main'
        mm[k] = dict(out['stores'][s]['models'][k], store=s)
        bins[k] = out['stores'][s]['relative_error_bins'][k]
    out['stores']['matched'] = {'models': mm, 'tau_with_headline': tau_block(mm), 'where_changed': where_changed(mm),
                                'relative_error_bins': bins}
    out['relative_error_bins'] = out['stores']['main']['relative_error_bins']
    out['n_prescribing_templates'] = len(out['prescribing_templates'])
    return out


# ------------------------------------------------------------------ writing

def f3(x) -> str:
    return '-' if x is None else f'{x:.3f}'


def render(res: dict) -> str:
    q = ' QUICK (1,000 resamples).' if res['quick'] else ''
    L = ['## Scoring-rule variants', '',
         f"Generated by `{res['generated_by']}`.{q} FAC (the item mean: correct 1, partial 0.5, incorrect and "
         f"unusable 0) under the headline rule and under four variants, each changing one clause: the absolute-value "
         f"clause off (a value matches only with its own sign); the response's last-digit window capped at "
         f"min(u(y^), 0.01|y|); the exact digits of the {res['n_prescribing_templates']} templates whose question "
         f"prescribes its rounding relaxed to the tolerance; a partial answer credited matched / of instead of 0.5. "
         f"Beside each variant, the verdicts it moves up / down. Intervals: 95% template bootstrap, "
         f"{res['resamples']:,} resamples; Δ is the variant minus the headline. The last row of each table is "
         f"Kendall's tau between the variant's ordering of the models and the headline's.", '']
    sym = res.get('symbolic_check_templates') or []
    L += [(f"The symbolic equivalence check is part of the rule on {len(sym)} templates and is applied after every "
           f"variant's label: " + ', '.join(f'`{t}`' for t in sym) + '.') if sym else
          'The rule has no symbolic equivalence step (answer.py carries none): symbolic answers are scored by the '
          'numbers they state.', '']
    titles = {'main': 'Main run (provider defaults)', 'matched': 'Matched settings (reasoning store where one exists)'}
    for s, blk in res['stores'].items():
        L += [f"### {titles.get(s, s)}", '',
              '| model | headline | ' + ' | '.join(f'{LABELS[v]} | Δ 95% CI | up / down' for v in VARIANTS) + ' |',
              '|---|---:|' + '---:|---:|---:|' * len(VARIANTS)]
        for k, m in blk['models'].items():
            cells = [f"{m[v]:.3f} | {m['delta_ci'][v][0]:+.3f} to {m['delta_ci'][v][1]:+.3f} | "
                     f"{m['changed'][v]['up']} / {m['changed'][v]['down']}" for v in VARIANTS]
            tag = f" ({m['store']})" if s == 'matched' and m.get('store') != 'main' else ''
            L.append(f"| `{k}`{tag} | {m['headline']:.3f} | " + ' | '.join(cells) + ' |')
        t = blk['tau_with_headline']
        L += ['| tau with the headline | 1 | ' + ' | '.join(f'{f3(t[v])} | | ' for v in VARIANTS) + '|', '']
        w = blk['where_changed']
        L += ['Where the changed verdicts sit: ' + '; '.join(
            f"{LABELS[v]}, {w[v]['verdicts']} verdicts on {w[v]['templates']} templates"
            + (' (most: ' + ', '.join(f'`{k}` {n}' for k, n in w[v]['top'].items()) + ')' if w[v]['top'] else '')
            for v in VARIANTS) + '.', '']
    L += ['### Relative error among accepted answers (main run)', '',
          'Readable responses the headline rule accepts, binned by the worst numeric part: for each numeric target, '
          'the relative error of the closest stated value the rule accepted (under the unit factors, the sign and '
          'the renderings it allows). `no numeric` counts accepted answers whose targets are words or labels only'
          + (f"; {res['relative_error_bins']['all']['symbolic']} more were accepted by the symbolic check alone and "
             'are left out of the bins' if sym else '') + '.', '',
          '| model | accepted | no numeric | ≤ 0.2% | 0.2 to 1% | 1 to 5% | > 5% |', '|---|---:|---:|---:|---:|---:|---:|']
    for k, b in res['relative_error_bins'].items():
        L.append(f"| {'all' if k == 'all' else '`' + k + '`'} | {b['accepted']} | {b['no_numeric']} | "
                 + ' | '.join(f"{b[n]} ({100 * b['share'][n]:.1f}%)" if b['share'][n] is not None else '-'
                              for n, _ in BINS) + ' |')
    o = res['relative_error_bins']['all']['over_1pct']
    L += ['', f"Of the {o['n']} accepted answers above 1%, {o['kept_when_bounded']} stay accepted with the last-digit "
          f"window capped (through the gold's own last digit or another reading). By answer kind: "
          + ', '.join(f'{k} {v}' for k, v in o['by_answer_type'].items()) + '. Most frequent templates: '
          + ', '.join(f'`{k}` {v}' for k, v in o['top_templates'].items()) + '.', '']
    return '\n'.join(L)


def write(res: dict) -> None:
    SECTIONS.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    OUT_MD.write_text(render(res), encoding='utf-8')
    print(f'wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}')


# ------------------------------------------------------------------ self-test

def selftest() -> int:
    """Each clause flipped on a hand-made case, and the copy with every switch off equal to answer.verdict."""
    item = lambda sol, q='', t='scalar': {'item_id': 'x#0', 'template_id': 'template_x', 'solution': sol,
                                          'question': q, 'answer_type': t}
    cases = [
        # (gold line, question, type, milestones, response, expected label per check)
        ('**Answer:** The deflection at the free end is 7.3 mm', '', 'scalar', (7.3,), '**Answer:** -7.326 mm',
         {'headline': 'correct', 'abs_clause_off': 'incorrect', 'last_digit_bounded': 'correct',
          'prescribed_relaxed': 'correct'}),
        ('**Answer:** The volume is 3,364 m3.', '', 'scalar', (3364.0,), '**Answer:** about 3.4 x 10^3 m3',
         {'headline': 'correct', 'abs_clause_off': 'correct', 'last_digit_bounded': 'incorrect',
          'prescribed_relaxed': 'correct'}),
        ('**Answer:** The probability is 0.4987.', '', 'scalar', (0.4987,), '**Answer:** P = 0.5',
         {'headline': 'correct', 'abs_clause_off': 'correct', 'last_digit_bounded': 'correct',
          'prescribed_relaxed': 'correct'}),
        ('**Answer:** The fraction is 0.1235.', 'Give it to 4 decimals, round half up.', 'scalar', (0.1235,),
         '**Answer:** 0.1234',
         {'headline': 'incorrect', 'abs_clause_off': 'incorrect', 'last_digit_bounded': 'incorrect',
          'prescribed_relaxed': 'correct'}),
        ('**Answer:** k = 74254 N/m, omega_n = 30.54 rad/s and zeta = 0.120', '', 'multipart', (74254.0, 30.54, 0.12),
         '**Answer:** k = 74254 N/m, omega_n = 30.54 rad/s, zeta = 0.31',
         {'headline': 'partial', 'abs_clause_off': 'partial', 'last_digit_bounded': 'partial',
          'prescribed_relaxed': 'partial'}),
    ]
    bad = 0
    for sol, q, typ, ms, text, want in cases:
        it = item(sol, q, typ)
        got = variant_verdicts(text, it, ms)
        ref = answer.verdict(text, it, answer.REL, ms)
        if got['headline'][0] != ref[0] or got['headline'][1:] != (ref[1]['matched'], ref[1]['of']):
            bad += 1
            print(f'FAIL copy != answer.verdict on {text!r}: {got["headline"]} vs {ref[0]}')
        for k, lab in want.items():
            if got[k][0] != lab:
                bad += 1
                print(f'FAIL {k}: want {lab} got {got[k][0]} | {text}')
    # per-part credit: 2 of 3 parts is 2/3, not 0.5
    got = variant_verdicts(cases[4][4], item(*cases[4][:3]), cases[4][3])['headline']
    if abs(got[1] / got[2] - 2 / 3) > 1e-12:
        bad += 1
        print(f'FAIL per-part credit: {got}')
    # the symbolic step (WS-E, at ANSWER FINAL): a response it accepts is correct under every variant; one it does
    # not accept keeps each variant's own label
    labels = ('headline', 'abs_clause_off', 'last_digit_bounded', 'prescribed_relaxed')
    vv = variant_verdicts(cases[4][4], item(*cases[4][:3]), cases[4][3], lambda lab, it, text: 'correct')
    if not (vv['raised'] and all(vv[k][0] == 'correct' for k in labels)):
        bad += 1
        print(f'FAIL symbolic step raising: {[vv[k][0] for k in labels]}')
    vv = variant_verdicts(cases[0][4], item(*cases[0][:3]), cases[0][3], lambda lab, it, text: lab)
    if vv['raised'] or [vv[k][0] for k in labels] != ['correct', 'incorrect', 'correct', 'correct']:
        bad += 1
        print(f'FAIL symbolic step declining: {[vv[k][0] for k in labels]}')
    # the relative error of an accepted answer: 3.4e3 for 3,364 is 1.07%, -7.326 for 7.3 is 0.36%
    for (sol, q, typ, ms, text, _), b in zip(cases[:3], ('0.2_1', '1_5', '0.2_1')):
        e = variant_verdicts(text, item(sol, q, typ), ms)['rel_errors']
        if not e or bin_of(max(e)) != b:
            bad += 1
            print(f'FAIL relative-error bin: want {b} got {e} | {text}')
    # the summary: one flip up and one down are counted, and the template interval brackets the mean
    rows = [{'template_id': f't{i % 3}', 'headline': 1.0, 'abs_clause_off': 0.0 if i == 0 else 1.0,
             'last_digit_bounded': 1.0, 'prescribed_relaxed': 1.0, 'per_part_credit': 1.0, 'bin': 'le_0.2'}
            for i in range(9)]
    rows.append({'template_id': 't0', 'headline': 0.5, 'abs_clause_off': 0.5, 'last_digit_bounded': 0.5,
                 'prescribed_relaxed': 0.5, 'per_part_credit': 2 / 3, 'bin': None})
    s = summarise({'rows': rows, 'readable': 10}, 200, 1)
    if s['changed']['abs_clause_off'] != {'up': 0, 'down': 1} or s['changed']['per_part_credit'] != {'up': 1, 'down': 0}:
        bad += 1
        print(f'FAIL change counts: {s["changed"]}')
    total = len(cases) * 5 + 1 + 2 + 3 + 1
    print(f'{total - bad}/{total} checks pass')
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--quick', action='store_true', help='1,000 resamples (also ENGTRACE_QUICK=1)')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    write(run(quick_mode(a.quick), a.workers))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
