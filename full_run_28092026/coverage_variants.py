"""Milestone Coverage variants and their tests, computed from the score store (WS-C2; review 1 M1 a to f,
review 2 W3, Q3, Q4). What MC measures, and whether its one top-tier separation survives its own controls.

    python -m full_run_28092026.coverage_variants              # writes results/coverage_variants.json and
                                                               # results/sections/coverage_variants.md
    python -m full_run_28092026.coverage_variants --quick      # 1,000 resamples, 10,000 permutations; QUICK in both files
    python -m full_run_28092026.coverage_variants --selftest   # e5_strict reproduced for every stored response, the
                                                               # variants, the regression and the matched family on
                                                               # hand-made data; writes nothing
    python -m full_run_28092026.coverage_variants --render     # the section again from the JSON written; computes nothing

Reads scores/<store>/<model>.jsonl (score.py: E3's `reached` per milestone, the answer row's `targets`, the token
counts), scores/<store>/e5/<model>.jsonl (judge.py: each milestone's `sources`), scores/milestones.json, the traces'
response text (through score.texts_matching, which refuses a trace that no longer matches its store row) and
results/matched_config.json (the models, their default stores and any reasoning store; an entry without a default
store, the Qwen Thinking sibling, is skipped). Nothing here calls a model or writes a store.

FOUR READINGS OF ONE RESPONSE'S COVERAGE, over the n milestones of its item:
  - as scored: (e3 + REACHED) / n, judge.py's e5_strict. Recomputed from the sources and asserted equal to the
    stored e5_strict and to analyze.coverage_value, the value Table 1 prints, for every response.
  - matching alone: E3's coverage, the milestones found by deterministic matching, with no judge.
  - route-adjusted: (e3 + REACHED) / (n - NOT_NEEDED). A milestone the judge rules not needed on the response's own
    route leaves the denominator; MISSING and UNJUDGED stay in it as not reached. A response whose milestones the
    judge rules all not needed (E3 found none of them) has no route-adjusted value: it is left out and counted.
  - intermediate-only: as scored, over the milestones that are not answer targets. THE TARGET RULE.
    evaluation's milestones.py records no answer line, so a target milestone is identified the way answer.targets
    makes a gold answer value a target in the first place: a milestone is a target when some number in the answer
    row's targets.numbers equals its value under one of the answer check's unit factors (answer.SCALES) within the
    milestone display tolerance (milestones.DISPLAY_TOL, 0.5%). The targets are the item's, the same in every row
    that has them (asserted). Where the gold's answer value is not among the milestones (answer.targets' fallback),
    no milestone is a target and all stay. An instance whose milestones are all targets is left without milestones.
A response has no value under any reading when its item has no milestones, or when the judge's call got no reply
(D-148, as analyze.coverage_value has it). An empty or unreadable response scores what it reached, an empty one
nothing, as Table 1 has it (D-170).

PER MODEL (each reading): the mean of the model's template means over the templates on which every model has a
value under that reading, with the template bootstrap (analyze.boot_mean), which for "as scored" is Table 1's MC and
its interval (analyze.q3_coverage, the same seeds); and on wrong answers (readable responses scoring 0, as
q3_coverage's coverage_wrong) the item mean with the template bootstrap of that ratio (analyze.cluster_mean). The
same for each model's reasoning store where matched_config.json names one, over the same templates.

THE PAIRS (each reading, default stores): q3_coverage's machinery imported from analyze.py, the per-template paired
differences over the common templates, the paired template bootstrap, the sign-flip test, Holm over the 55 pairs of
that reading, the detectable difference at 80% power (plain and at the strictest Holm step). Every reading uses the
same seeds as q3_coverage, so "as scored" reproduces its pairs and the readings differ by their values alone.

THE SEPARATION RULE (the results sentence): Claude Sonnet 5 and DeepSeek V4.1 Flash separate when the pair's Holm-
adjusted p is below 0.05 under matching alone and under route-adjusted MC, with the difference in the direction it
has as scored. The p under each reading is printed beside the outcome.

THE VERBOSITY MODEL: response-level least squares of MC as scored on the number of numeric values the response
displays and on its visible completion tokens, with template fixed effects (within-template demeaning), pooled over
the default stores' models, over the readable responses with a value. The count is arith's reader: the literals
arith.ROUNDED_LITERAL finds in the text after arith.normalise (thousands groupings and scientific notation folded
into one number, a subscript digit not a number), exponents after `**` left out. Visible completion tokens are
completion_tokens minus reasoning_tokens, floored at zero (a few rows of one endpoint report more reasoning than
completion tokens; they are counted). The slopes' 95% intervals resample templates, every response of a drawn
template coming with it (B draws; the fixed-effects fit is re-solved from per-template cross-products). Beside them:
the same model with template-by-model fixed effects (the slope within a model and template), the within-model
Spearman of MC against the count, and the Spearman across the models' means.

MC AGAINST REASONING TOKENS (descriptive): per model, its responses with a value cut into quartiles of
reasoning_tokens within the model (ranks, ties in row order), MC as scored per quartile, the median tokens and the
share of empty responses per quartile, at the default store where the model reasons there (reasoning tokens on at
least 90% of its rows) and otherwise at its reasoning store.

MATCHED SETTINGS (the headline configuration, key `matched`): each model at its reasoning store where
matched_config.json names one, else at its default store, whose responses are the same in both configurations. The
family again on that selection, in the default family's shapes: the templates on which every model has a value under
each reading, the per-model summaries (seeds by the model's place in the roster, so an unchanged model's values equal
its default ones, asserted), the 55 pairs with Holm per reading, the separation rule, the verbosity model over the
matched responses and the reasoning-token quartiles, each model at its matched store. The default family stays as it is.

QUICK MODE (--quick or ENGTRACE_QUICK=1): 1,000 bootstrap draws and 10,000 sign flips instead of analyze's 10,000
and 100,000; both outputs say QUICK. The integration pass runs without it.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import analyze, score  # noqa: E402  score puts the evaluators on sys.path
import answer  # noqa: E402
import arith  # noqa: E402
import milestones  # noqa: E402

OUT = HERE / 'results'
JSON_OUT = OUT / 'coverage_variants.json'
MD_OUT = OUT / 'sections' / 'coverage_variants.md'
MATCHED = OUT / 'matched_config.json'
VARIANTS = ('as_scored', 'matching_only', 'route_adjusted', 'intermediate_only')
LABEL = {'as_scored': 'as scored (E3 + judge)', 'matching_only': 'matching alone (E3)',
         'route_adjusted': 'route-adjusted', 'intermediate_only': 'intermediate-only'}
PAIR = ('claude-sonnet-5', 'deepseek-v4.1-flash')
REASONS = 0.9               # a store "reasons" when reasoning_tokens > 0 on at least this share of its rows
EXPONENT = re.compile(r'\*\*\s*\(?\s*$')
QUICK_B, QUICK_TEST = 1_000, 10_000


# ------------------------------------------------------------------ the roster

def roster() -> list[dict]:
    """matched_config.json's models with a default store present, in its order; no hard-coded eleven."""
    cfg = json.loads(MATCHED.read_text(encoding='utf-8'))
    out = []
    for m in cfg['models']:
        if not m.get('default_store') or m.get('run') is False:
            continue
        if not (score.SCORES / m['default_store'] / f"{m['model']}.jsonl").exists():
            continue
        rs = m.get('reasoning_store')
        out.append({'key': m['model'], 'default': m['default_store'],
                    'reasoning': rs if rs and (score.SCORES / rs / f"{m['model']}.jsonl").exists() else None})
    return out


def matched_store(m: dict) -> str:
    """A roster entry's store at matched settings: its reasoning store where there is one, else its default."""
    return m['reasoning'] or m['default']


# ------------------------------------------------------------------ one response

def target_flags(ms: list[dict], numbers) -> list[bool]:
    """Which of an item's milestones are answer targets: equal to a target number under one of the answer
    check's unit factors, within the milestone display tolerance (the docstring's target rule)."""
    return [any(milestones.close(m['value'], v * sc, milestones.DISPLAY_TOL) for v in numbers for sc in answer.SCALES)
            for m in ms]


def variants(r: dict, e5row: dict | None, is_target: list[bool]) -> dict:
    """The four readings of one response's coverage; None where a reading has no value."""
    none = dict.fromkeys(VARIANTS)
    base = analyze.coverage_value(r, e5row)
    if base is None:
        return none
    n = r['milestones_required']
    reached = r['e3']['reached']
    src = e5row['sources'] if e5row else ['e3' if x else 'UNJUDGED' for x in reached]
    assert len(src) == n == len(reached) == len(is_target), r['item_id']
    assert [s == 'e3' for s in src] == list(reached), f"{r['item_id']}: E5's e3 sources differ from E3's reached"
    got = [s in ('e3', 'REACHED') for s in src]
    strict = sum(got) / n
    assert abs(strict - base) < 1e-12, f"{r['item_id']}: as scored {strict} != coverage_value {base}"
    if e5row:
        assert abs(strict - e5row['e5_strict']) < 1e-12, f"{r['item_id']}: as scored {strict} != e5_strict"
    needed = n - src.count('NOT_NEEDED')
    keep = [not t for t in is_target]
    return {'as_scored': strict, 'matching_only': r['e3']['coverage'],
            'route_adjusted': sum(got) / needed if needed else None,
            'intermediate_only': (sum(g for g, k in zip(got, keep) if k) / sum(keep)) if any(keep) else None}


def numbers_shown(text: str) -> int:
    """The numeric values a response displays, read by arith's reader (the docstring)."""
    t = arith.normalise(text or '')
    return sum(1 for m in arith.ROUNDED_LITERAL.finditer(t) if not EXPONENT.search(t[max(0, m.start() - 4):m.start()]))


def visible_tokens(r: dict) -> float:
    return max(0.0, float((r.get('completion_tokens') or 0) - (r.get('reasoning_tokens') or 0)))


# ------------------------------------------------------------------ loading

def item_targets(stores: list[tuple[str, str]]) -> dict[str, tuple]:
    """The answer targets per item from every answered row of the stores read; the same in each (asserted)."""
    seen = {}
    for store, key in stores:
        for i, r in analyze.load(store, key).items():
            if r.get('answer'):
                t = tuple(r['answer']['targets']['numbers'])
                assert seen.setdefault(i, t) == t, f'{i}: targets differ between rows ({seen[i]} vs {t})'
    return seen


def store_values(store: str, key: str, ms_all: dict, targets: dict) -> tuple[dict, dict]:
    rows = analyze.load(store, key)
    e5 = analyze.load_stage('e5', key, store)
    vals = {i: variants(r, e5.get(i) if e5 else None, target_flags(ms_all[i], targets.get(i, ()))) for i, r in rows.items()}
    return rows, vals


def source_counts(store: str, key: str, vals: dict) -> dict:
    """The judge's sources over the responses with a value: what route adjustment takes out of the denominator."""
    e5 = analyze.load_stage('e5', key, store) or {}
    c = collections.Counter(s for i, x in vals.items() if x['as_scored'] is not None and i in e5 for s in e5[i]['sources'])
    unmatched = sum(n for s, n in c.items() if s != 'e3')
    return {**{s: c[s] for s in ('e3', 'REACHED', 'NOT_NEEDED', 'MISSING', 'UNJUDGED')}, 'unmatched': unmatched,
            'not_needed_share_of_unmatched': c['NOT_NEEDED'] / unmatched if unmatched else None}


def target_rule(ms_all: dict, targets: dict, rows: dict) -> dict:
    items = set(rows)
    no_ms = [i for i in items if not ms_all[i]]
    flags = {i: target_flags(ms_all[i], targets.get(i, ())) for i in items if ms_all[i]}
    all_t = [i for i, f in flags.items() if all(f)]
    by_t = collections.defaultdict(list)
    for i, r in rows.items():
        by_t[r['template_id']].append(bool(ms_all[i]) and not all(flags.get(i, [True])))
    return {'targets': ('a milestone is an answer target when a number in the answer row\'s targets.numbers equals its '
                        'value under one of the answer check\'s unit factors (answer.SCALES) within the milestone '
                        'display tolerance (0.5%), the test answer.targets uses to make a gold answer value a target; '
                        'milestones.py records no answer line'),
            'instances': len(items), 'instances_without_milestones': len(no_ms) + len(all_t),
            'no_milestones_at_all': len(no_ms), 'all_milestones_are_targets': len(all_t),
            'target_milestones': sum(sum(f) for f in flags.values()),
            'milestones': sum(len(f) for f in flags.values()),
            'answer_value_not_a_milestone': sum(1 for i, f in flags.items() if not any(f)),
            'items_without_targets': sum(1 for i in items if i not in targets),
            'templates_without_intermediate_milestones': sum(1 for v in by_t.values() if not any(v))}


# ------------------------------------------------------------------ per model and pairs

def matrix(rows: dict, vals: dict, v: str, templates: list[str]) -> np.ndarray:
    by_t = analyze.per_template(rows.values(), lambda r: vals[r['item_id']][v],
                                lambda r: vals[r['item_id']][v] is not None)
    return np.array([np.mean(by_t[t]) if by_t.get(t) else np.nan for t in templates])


def summary(rows: dict, vals: dict, c: np.ndarray, use: np.ndarray, v: str, i: int) -> dict:
    """One model under one reading: q3_coverage's estimands and seeds (5000 + i, 5150 + i)."""
    f = lambda r: vals[r['item_id']][v]
    read = [r for r in rows.values() if not r['unusable'] and f(r) is not None]
    w, w_ci, w_n = analyze.cluster_mean(analyze.per_template(read, f, lambda r: r['score'] == 0), 5150 + i)
    return {'all': float(np.nanmean(c[use])), 'ci': analyze.boot_mean(c[use], 5000 + i),
            'wrong': w, 'wrong_ci': w_ci, 'wrong_n': w_n,
            'responses': sum(1 for r in rows.values() if f(r) is not None)}


def left_out(vals: dict) -> int:
    """Responses with a value as scored and none route-adjusted (the judge rules all their milestones not needed)."""
    return sum(1 for x in vals.values() if x['as_scored'] is not None and x['route_adjusted'] is None)


def family(keys: list[str], got: dict, templates: list[str], draws: int) -> tuple[dict, dict, dict, dict]:
    """One configuration (got: model -> (rows, vals)): per reading, the models-by-templates matrix, the templates on
    which every model has a value, each model's summaries (seeds by its place in keys) and the pairs."""
    C = {v: np.vstack([matrix(*got[k], v, templates) for k in keys]) for v in VARIANTS}
    use = {v: ~np.isnan(C[v]).any(axis=0) for v in VARIANTS}
    summ = {k: {v: summary(*got[k], C[v][i], use[v], v, i) for v in VARIANTS} for i, k in enumerate(keys)}
    return C, use, summ, {v: pairs(C[v], use[v], keys, draws) for v in VARIANTS}


def pairs(C: np.ndarray, use: np.ndarray, keys: list[str], draws: int) -> list[dict]:
    """q3_coverage's pairs: seeds 6000 + n (interval) and 7000 + n (sign flips), Holm over the family."""
    out = []
    n_pairs = len(keys) * (len(keys) - 1) // 2
    for n, (i, j) in enumerate(itertools.combinations(range(len(keys)), 2)):
        d = (C[i] - C[j])[use]
        dt = analyze.detectable_paired(d, n_pairs)
        out.append({'a': keys[i], 'b': keys[j], 'diff': float(np.mean(d)), 'ci': analyze.boot_mean(d, 6000 + n),
                    'p': analyze.sign_flip_p(d, 7000 + n, draws=draws),
                    'detectable': dt['detectable'], 'detectable_holm': dt['detectable_holm']})
    for o, ph in zip(out, analyze.holm([p['p'] for p in out])):
        o['p_holm'] = ph
    return out


def find_pair(ps: list[dict], a: str, b: str) -> dict | None:
    """The pair as a minus b, whichever order the family holds it in."""
    for p in ps:
        if (p['a'], p['b']) == (a, b):
            return p
        if (p['a'], p['b']) == (b, a):
            return {**p, 'a': a, 'b': b, 'diff': -p['diff'], 'ci': [-p['ci'][1], -p['ci'][0]]}
    return None


def separation(fam: dict) -> dict:
    got = {v: find_pair(fam[v], *PAIR) for v in VARIANTS}
    if any(g is None for g in got.values()):
        return {'claude_vs_deepseek': None}
    sign = np.sign(got['as_scored']['diff'])
    holds = all(got[v]['p_holm'] < 0.05 and np.sign(got[v]['diff']) == sign
                for v in ('as_scored', 'matching_only', 'route_adjusted'))
    return {'claude_vs_deepseek': {**{v: got[v]['p_holm'] for v in VARIANTS},
                                   'diff': {v: got[v]['diff'] for v in VARIANTS},
                                   'holds': bool(holds)}}


# ------------------------------------------------------------------ the verbosity model

def fe_crossproducts(y: np.ndarray, X: np.ndarray, cell: np.ndarray, cluster: np.ndarray, n_clusters: int):
    """Per cluster, the cross-products of X and y demeaned within each fixed-effect cell (cells nest in clusters),
    so a fit on any resample of clusters is one solve of their weighted sums."""
    k = X.shape[1]
    Sxx, Sxy = np.zeros((n_clusters, k, k)), np.zeros((n_clusters, k))
    order = np.argsort(cell, kind='stable')
    bounds = np.flatnonzero(np.diff(cell[order])) + 1
    for g in np.split(order, bounds):
        Xd, yd = X[g] - X[g].mean(axis=0), y[g] - y[g].mean()
        c = cluster[g[0]]
        Sxx[c] += Xd.T @ Xd
        Sxy[c] += Xd.T @ yd
    return Sxx, Sxy


def fe_fit(y, X, cell, cluster, n_clusters, draws: int, seed: int) -> dict:
    Sxx, Sxy = fe_crossproducts(y, X, cell, cluster, n_clusters)
    beta = np.linalg.solve(Sxx.sum(axis=0), Sxy.sum(axis=0))
    rng = np.random.default_rng(seed)
    W = np.stack([np.bincount(rng.integers(0, n_clusters, n_clusters), minlength=n_clusters) for _ in range(draws)])
    bb = np.linalg.solve(np.einsum('bt,tij->bij', W, Sxx), np.einsum('bt,ti->bi', W, Sxy)[..., None])[..., 0]
    return {'beta': beta.tolist(), 'ci': [analyze.pct(bb[:, j]) for j in range(X.shape[1])]}


def verbosity_rows(key: str, rows: dict, vals: dict, texts: dict) -> list[dict]:
    """The verbosity model's rows for one model: its readable responses with a value, in store order."""
    return [{'model': key, 'template': r['template_id'], 'mc': vals[i]['as_scored'], 'numbers': numbers_shown(texts[i]),
             'tokens': visible_tokens(r)} for i, r in rows.items() if not r['unusable'] and vals[i]['as_scored'] is not None]


def floored(got: dict, keys: list[str]) -> dict:
    return {k: sum(1 for r in got[k][0].values() if (r.get('completion_tokens') or 0) < (r.get('reasoning_tokens') or 0))
            for k in keys}


def verbosity(data: list[dict], draws: int) -> dict:
    """data: one dict per readable response with a value (model, template, mc, numbers, tokens)."""
    tids = {t: j for j, t in enumerate(sorted({d['template'] for d in data}))}
    keys = sorted({d['model'] for d in data})
    y = np.array([d['mc'] for d in data])
    X = np.column_stack([[d['numbers'] / 10 for d in data], [d['tokens'] / 1000 for d in data]])
    tpl = np.array([tids[d['template']] for d in data])
    mdl = np.array([keys.index(d['model']) for d in data])
    pooled = fe_fit(y, X, tpl, tpl, len(tids), draws, 9100)
    within = fe_fit(y, X, tpl * len(keys) + mdl, tpl, len(tids), draws, 9200)
    sp, sp_t, means = {}, {}, {}
    for j, k in enumerate(keys):
        s = mdl == j
        sp[k] = float(stats.spearmanr(X[s, 0], y[s]).statistic)
        sp_t[k] = float(stats.spearmanr(X[s, 1], y[s]).statistic)
        means[k] = (float(y[s].mean()), float(X[s, 0].mean() * 10), float(np.median(X[s, 0]) * 10),
                    float(np.median(X[s, 1]) * 1000))
    between = stats.spearmanr([means[k][0] for k in keys], [means[k][1] for k in keys]).statistic
    return {'slope_numbers': pooled['beta'][0], 'ci': pooled['ci'][0],
            'slope_tokens': pooled['beta'][1], 'ci_tokens': pooled['ci'][1],
            'units': {'numbers': 'MC per 10 numeric values shown', 'tokens': 'MC per 1,000 visible completion tokens'},
            'within_model': {'slope_numbers': within['beta'][0], 'ci': within['ci'][0],
                             'slope_tokens': within['beta'][1], 'ci_tokens': within['ci'][1]},
            'spearman_within': sp, 'spearman_within_tokens': sp_t,
            'spearman_between_models': float(between),
            'per_model': {k: {'mc': m[0], 'numbers_mean': m[1], 'numbers_median': m[2], 'tokens_median': m[3]}
                          for k, m in means.items()},
            'responses': len(data), 'templates': len(tids),
            'iqr_numbers': [float(np.percentile(X[:, 0], 25) * 10), float(np.percentile(X[:, 0], 75) * 10)],
            'iqr_tokens': [float(np.percentile(X[:, 1], 25) * 1000), float(np.percentile(X[:, 1], 75) * 1000)]}


def reasoning_quartiles(rows: dict, vals: dict) -> dict | None:
    sel = [r for r in rows.values() if vals[r['item_id']]['as_scored'] is not None]
    tok = np.array([float(r.get('reasoning_tokens') or 0) for r in sel])
    if not sel or np.mean(tok > 0) < REASONS:
        return None
    q = np.empty(len(sel), dtype=int)
    q[np.argsort(tok, kind='stable')] = np.arange(len(sel)) * 4 // len(sel)
    mc = np.array([vals[r['item_id']]['as_scored'] for r in sel])
    empty = np.array([r['status'] != 'answered' for r in sel])
    out = {f'q{j + 1}': float(mc[q == j].mean()) for j in range(4)}
    out.update({'tokens_median': [float(np.median(tok[q == j])) for j in range(4)],
                'empty_share': [float(empty[q == j].mean()) for j in range(4)],
                'responses': len(sel), 'spearman': float(stats.spearmanr(tok, mc).statistic)})
    return out


# ------------------------------------------------------------------ the run

def sha_file(p: Path) -> str | None:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def provenance(stores: set[str]) -> dict:
    git = lambda *a: subprocess.run(['git', *a], cwd=REPO, capture_output=True, text=True).stdout.strip()
    return {'git': git('rev-parse', 'HEAD'), 'script_sha256_lf': score.norm_sha(Path(__file__)),
            'matched_config_sha256': sha_file(MATCHED),
            'stores': {s: {'config_sha256': sha_file(score.SCORES / s / 'CONFIG.json'),
                           'e5_config_sha256': sha_file(score.SCORES / s / 'e5' / 'CONFIG.json')} for s in sorted(stores)},
            'milestones_sha256': sha_file(score.MILESTONES),
            'written_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}


def run(quick: bool) -> dict:
    if quick:
        analyze.B = QUICK_B                 # boot_mean and cluster_mean read the module's B at call time
    draws_test = QUICK_TEST if quick else analyze.B_TEST
    ros = roster()
    keys = [m['key'] for m in ros]
    ms_all = json.loads(score.MILESTONES.read_text(encoding='utf-8'))
    stores = [(m['default'], m['key']) for m in ros] + [(m['reasoning'], m['key']) for m in ros if m['reasoning']]
    targets = item_targets(stores)
    default = {m['key']: store_values(m['default'], m['key'], ms_all, targets) for m in ros}
    ref = next(iter(default.values()))[0]
    for k, (rows, _v) in default.items():
        if set(rows) != set(ref):
            raise SystemExit(f'{k}: its default store holds a different item set')
    templates = sorted({r['template_id'] for r in ref.values()})
    C, use, summ, fam = family(keys, default, templates, draws_test)

    reasoning = {m['key']: store_values(m['reasoning'], m['key'], ms_all, targets) for m in ros if m['reasoning']}
    for k, (rows, _v) in reasoning.items():
        if set(rows) != set(ref):
            raise SystemExit(f'{k}: its reasoning store holds a different item set')
    models = {}
    for i, m in enumerate(ros):
        k = m['key']
        rows, vals = default[k]
        block = {**summ[k], 'store': m['default'], 'route_adjusted_left_out': left_out(vals),
                 'sources': source_counts(m['default'], k, vals)}
        if m['reasoning']:
            rr, rv = reasoning[k]
            block['matched_store'] = {'store': m['reasoning'],
                                      **{v: summary(rr, rv, matrix(rr, rv, v, templates), use[v], v, i) for v in VARIANTS}}
        else:
            block['matched_store'] = None
        models[k] = block

    vrows = {m['key']: verbosity_rows(m['key'], *default[m['key']],
                                      score.texts_matching(m['default'], m['key'], default[m['key']][0].values()))
             for m in ros}
    verb = verbosity([d for k in keys for d in vrows[k]], analyze.B)
    verb['tokens_floored'] = floored(default, keys)

    rq = {}
    for m in ros:
        k = m['key']
        got = reasoning_quartiles(*default[k])
        if got is not None:
            rq[k] = {'store': m['default'], **got}
        elif m['reasoning']:
            got = reasoning_quartiles(*reasoning[k])
            rq[k] = {'store': m['reasoning'], **got} if got else None
        else:
            rq[k] = None

    mstores = {m['key']: matched_store(m) for m in ros}
    matched = {k: reasoning[k] if k in reasoning else default[k] for k in keys}
    _Cm, use_m, summ_m, fam_m = family(keys, matched, templates, draws_test)
    m_models = {}
    for k in keys:
        rows, vals = matched[k]
        if k in reasoning:
            same = [v for v in VARIANTS if (use_m[v] == use[v]).all()]
            assert all(jsonable(summ_m[k][v]) == jsonable(models[k]['matched_store'][v]) for v in same), k
        else:
            assert jsonable(summ_m[k]) == jsonable(summ[k]), f'{k}: unchanged at matched settings, yet its values differ'
        m_models[k] = {**summ_m[k], 'store': mstores[k], 'route_adjusted_left_out': left_out(vals),
                       'sources': source_counts(mstores[k], k, vals)}
    m_vrows = {k: verbosity_rows(k, *matched[k], score.texts_matching(mstores[k], k, matched[k][0].values()))
               if k in reasoning else vrows[k] for k in keys}
    m_verb = verbosity([d for k in keys for d in m_vrows[k]], analyze.B)
    m_verb['tokens_floored'] = floored(matched, keys)
    m_rq = {}
    for k in keys:
        got = reasoning_quartiles(*matched[k])
        m_rq[k] = {'store': mstores[k], **got} if got else None

    return {'quick': quick, 'draws': {'bootstrap': analyze.B, 'sign_flips': draws_test},
            'rule': {**target_rule(ms_all, targets, ref),
                     'route_adjusted': ('(e3 + REACHED) / (n - NOT_NEEDED); a response whose milestones the judge '
                                        'rules all not needed has no value and is left out (count per model)')},
            'keys': keys,
            'templates': {v: int(use[v].sum()) for v in VARIANTS},
            'models': models,
            'pairs': fam,
            'separation_rule': separation(fam),
            'verbosity': verb,
            'reasoning_tokens': rq,
            'matched': {'templates': {v: int(use_m[v].sum()) for v in VARIANTS},
                        'models': m_models,
                        'pairs': fam_m,
                        'separation_rule': separation(fam_m),
                        'verbosity': m_verb,
                        'reasoning_tokens': m_rq,
                        'stores': mstores},
            'provenance': provenance({s for s, _k in stores})}


# ------------------------------------------------------------------ the section

def f3(x) -> str:
    return '' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:.3f}'


def ci(c) -> str:
    return f'[{c[0]:.3f}, {c[1]:.3f}]'


def pfmt(p) -> str:
    return f'{p:.4f}' if p >= 0.0001 else f'{p:.1e}'


def coverage_rows(models: dict, keys: list[str], store: bool = False) -> list[str]:
    """The per-model table's rows: each reading with its interval, then the wrong-answer means."""
    return [f'| {k} | ' + (f'{models[k]["store"]} | ' if store else '')
            + ' | '.join(f'{f3(models[k][v]["all"])} {ci(models[k][v]["ci"])}' for v in VARIANTS) + ' | '
            + ' | '.join(f3(models[k][v]['wrong']) for v in VARIANTS) + f' | {models[k]["as_scored"]["wrong_n"]} |'
            for k in keys]


def separation_lines(sep: dict | None, where: str) -> list[str]:
    if sep is None:
        return [f'Claude Sonnet 5 or DeepSeek V4.1 Flash is missing from {where}; the rule is not evaluated.']
    L = [f'Claude Sonnet 5 minus DeepSeek V4.1 Flash, Holm-adjusted p over each reading\'s 55 pairs: '
         + '; '.join(f'{LABEL[v]} {sep["diff"][v]:+.3f}, p = {pfmt(sep[v])}' for v in VARIANTS) + '. ']
    sign = np.sign(sep['diff']['as_scored'])
    ok = {v: sep[v] < 0.05 and np.sign(sep['diff'][v]) == sign for v in ('as_scored', 'matching_only', 'route_adjusted')}
    said = lambda vs: ' and '.join(f'{LABEL[v]} (p = {pfmt(sep[v])})' for v in vs)
    yes, no = [v for v in ok if ok[v]], [v for v in ok if not ok[v]]
    if sep['holds']:
        L.append(f'**The separation holds**: the pair separates after Holm, in one direction, under {said(yes)}.')
    else:
        got = f'separates after Holm under {said(yes)} but not under ' if yes else 'does not separate after Holm under '
        L.append(f'**The separation does not hold**: the pair {got}{said(no)}, so the results sentence does not '
                 'claim it.')
    return L


def pair_lines(P: dict) -> list[str]:
    """The 55 pairs under each reading: the count separated, the table, the detectable difference."""
    L = ['Difference a minus b in the mean of template means, Holm-adjusted sign-flip p within each reading\'s family. '
         'Pairs separated at 0.05: ' + ', '.join(f'{LABEL[v]} {sum(p["p_holm"] < 0.05 for p in P[v])}'
                                                 for v in VARIANTS) + '.', '',
         '| a | b | ' + ' | '.join(f'{LABEL[v]}: diff | p Holm' for v in VARIANTS) + ' |',
         '|---|---|' + '---|---|' * len(VARIANTS)]
    for n in range(len(P['as_scored'])):
        ps = [P[v][n] for v in VARIANTS]
        L.append(f'| {ps[0]["a"]} | {ps[0]["b"]} | ' + ' | '.join(
            f'{p["diff"]:+.3f} | {pfmt(p["p_holm"])}{"*" if p["p_holm"] < 0.05 else ""}' for p in ps) + ' |')
    L += ['', 'Smallest detectable difference (80% power, at the strictest Holm step), median over the 55 pairs: '
          + ', '.join(f'{LABEL[v]} {np.median([p["detectable_holm"] for p in P[v]]):.3f}' for v in VARIANTS) + '.']
    return L


def verbosity_text(vb: dict) -> str:
    w = vb['within_model']
    return (f'Response-level least squares of MC as scored on the numeric values shown and the visible completion tokens, '
            f'template fixed effects, pooled over {len(vb["per_model"])} models ({vb["responses"]:,} readable responses, '
            f'{vb["templates"]} templates); 95% intervals resample templates. Per 10 numeric values: '
            f'{vb["slope_numbers"]:+.4f} {ci4(vb["ci"])}; per 1,000 visible tokens: {vb["slope_tokens"]:+.4f} '
            f'{ci4(vb["ci_tokens"])}. Within model and template (template-by-model fixed effects): '
            f'{w["slope_numbers"]:+.4f} {ci4(w["ci"])} and {w["slope_tokens"]:+.4f} {ci4(w["ci_tokens"])}. '
            f'Interquartile range of the count: {vb["iqr_numbers"][0]:.0f} to {vb["iqr_numbers"][1]:.0f}; of visible '
            f'tokens: {vb["iqr_tokens"][0]:,.0f} to {vb["iqr_tokens"][1]:,.0f}. Spearman across the models\' means '
            f'(MC against the mean count): {vb["spearman_between_models"]:+.3f}. Rows with more reasoning than completion '
            f'tokens, floored at zero visible tokens: '
            + ', '.join(f'{k} {n}' for k, n in vb['tokens_floored'].items() if n) + '.')


def matched_lines(res: dict) -> list[str]:
    """The "Matched settings" subsection: the family again with each model at its matched store."""
    M, keys = res['matched'], res['keys']
    mm, D = M['models'], res['pairs']
    rerun = [k for k in keys if M['stores'][k] != res['models'][k]['store']]
    add = [mm[k]['as_scored']['all'] - mm[k]['matching_only']['all'] for k in keys]
    L = ['', '### Matched settings', '',
         'Each model at its reasoning setting where its endpoint offers one, else at its default: '
         + ', '.join(f'{k} at `{M["stores"][k]}`' for k in rerun) + '; the other '
         f'{len(keys) - len(rerun)} at their default store, the same responses as above (their values are asserted '
         'equal to the default family\'s). The same readings, seeds and tests as above. Templates in each reading\'s '
         'comparison: ' + ', '.join(f'{LABEL[v]} {M["templates"][v]}' for v in VARIANTS) + '. The judge adds '
         f'{min(add):.3f} to {max(add):.3f} to a model\'s MC over matching alone. Route-adjusted, responses left out: '
         + ', '.join(f'{k} {mm[k]["route_adjusted_left_out"]}' for k in keys) + '.', '',
         '| model | store | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong: matching | wrong: route-adj. | wrong: interm. | wrong n |',
         '|---|---|---|---|---|---|---|---|---|---|---|'] + coverage_rows(mm, keys, store=True)
    L += ['', '#### The separation rule at matched settings', '']
    L += separation_lines(M['separation_rule']['claude_vs_deepseek'], 'the matched stores')
    L += ['', '#### The 55 pairs at matched settings', '']
    for v in VARIANTS:
        sep = {(p['a'], p['b']) for p in M['pairs'][v] if p['p_holm'] < 0.05}
        was = {(p['a'], p['b']) for p in D[v] if p['p_holm'] < 0.05}
        name = lambda s: ', '.join(f'{a} vs {b}' for a, b in sorted(s)) or 'none'
        L.append(f'- {LABEL[v]}: {len(sep)} of {len(M["pairs"][v])} separated after Holm; separated here but not in '
                 f'the default family: {name(sep - was)}; there but not here: {name(was - sep)}.')
    L += [''] + pair_lines(M['pairs'])
    L += ['', '#### Verbosity at matched settings', '', verbosity_text(M['verbosity'])]
    if M['reasoning_tokens'] == res['reasoning_tokens']:
        L += ['', 'MC against reasoning tokens: the quartiles above, each already at its model\'s matched store.']
    else:
        L += ['', '#### MC against reasoning tokens at matched settings', ''] + rq_lines(M['reasoning_tokens'], keys)
    return L


def rq_lines(rq: dict, keys: list[str]) -> list[str]:
    L = ['| model | store | Q1 | Q2 | Q3 | Q4 | median tokens Q1 to Q4 | empty Q1 to Q4 | Spearman |',
         '|---|---|---|---|---|---|---|---|---|']
    for k in keys:
        x = rq[k]
        if x is None:
            L.append(f'| {k} | | no reasoning tokens in any store | | | | | | |')
            continue
        L.append(f'| {k} | {x["store"]} | ' + ' | '.join(f3(x[f"q{j}"]) for j in range(1, 5)) + ' | '
                 + ', '.join(f'{t:,.0f}' for t in x['tokens_median']) + ' | '
                 + ', '.join(f'{e:.3f}' for e in x['empty_share']) + f' | {x["spearman"]:+.3f} |')
    return L


def render(res: dict) -> str:
    q = ' QUICK' if res['quick'] else ''
    R, keys = res['rule'], res['keys']
    L = [f'## Coverage variants{q}', '',
         f'Printed by `coverage_variants.py`{" in quick mode (" + format(res["draws"]["bootstrap"], ",") + " bootstrap draws, " + format(res["draws"]["sign_flips"], ",") + " sign flips): QUICK, not for the paper" if res["quick"] else ""}; '
         'the method is in its docstring. MC as scored is Table 1\'s Milestone Coverage, recomputed from the judge\'s '
         'per-milestone sources and checked equal to the stored `e5_strict` for every response. Each model\'s value '
         'is the mean of its template means over the templates on which every model has a value under that reading, '
         'with a 95% template-bootstrap interval; "wrong" is the item mean over readable responses scoring 0.', '',
         f'**Target rule.** {R["targets"][0].upper() + R["targets"][1:]}. Of {R["instances"]:,} instances, '
         f'{R["no_milestones_at_all"]} have no milestones and {R["all_milestones_are_targets"]} have only target '
         f'milestones, so {R["instances_without_milestones"]} are left without milestones under intermediate-only MC '
         f'({R["templates_without_intermediate_milestones"]} templates have none left). {R["target_milestones"]:,} of '
         f'{R["milestones"]:,} milestones are targets; on {R["answer_value_not_a_milestone"]} instances the gold\'s '
         'answer value is not a milestone, so all their milestones stay. '
         f'Templates in each reading\'s comparison: ' + ', '.join(f'{LABEL[v]} {res["templates"][v]}' for v in VARIANTS) + '.', '',
         f'**Route-adjusted.** {R["route_adjusted"][0].upper() + R["route_adjusted"][1:]}: '
         + ', '.join(f'{k} {res["models"][k]["route_adjusted_left_out"]}' for k in keys) + '. '
         'Of the milestones E3 leaves unmatched, the share the judge rules not needed (and the count unmatched): '
         + ', '.join(f'{k} {res["models"][k]["sources"]["not_needed_share_of_unmatched"]:.3f} '
                     f'({res["models"][k]["sources"]["unmatched"]:,})' for k in keys) + '.', '',
         '| model | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong: matching | wrong: route-adj. | wrong: interm. | wrong n |',
         '|---|---|---|---|---|---|---|---|---|---|'] + coverage_rows(res['models'], keys)
    rs = [k for k in keys if res['models'][k]['matched_store']]
    if rs:
        L += ['', 'At the reasoning store `matched_config.json` names (same templates):', '',
              '| model | store | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong n |',
              '|---|---|---|---|---|---|---|---|']
        for k in rs:
            b = res['models'][k]['matched_store']
            L.append(f'| {k} | {b["store"]} | ' + ' | '.join(f'{f3(b[v]["all"])} {ci(b[v]["ci"])}' for v in VARIANTS)
                     + f' | {f3(b["as_scored"]["wrong"])} | {b["as_scored"]["wrong_n"]} |')
    L += ['', '### The separation rule', ''] + separation_lines(res['separation_rule']['claude_vs_deepseek'],
                                                                 'the default stores')
    L += ['', '### The 55 pairs under each reading', ''] + pair_lines(res['pairs'])
    vb = res['verbosity']
    L += ['', '### Verbosity', '', verbosity_text(vb), '',
          'The within-model Spearman is over a model\'s readable responses with no template control, so it mixes in '
          'the templates\' difficulty (a harder template draws a longer response and a lower MC); the regression\'s '
          'fixed effects remove that.', '',
          '| model | MC (readable) | numeric values, mean | median | visible tokens, median | Spearman MC vs count | vs tokens |',
          '|---|---|---|---|---|---|---|']
    for k in keys:
        pm = vb['per_model'][k]
        L.append(f'| {k} | {pm["mc"]:.3f} | {pm["numbers_mean"]:.1f} | {pm["numbers_median"]:.0f} | '
                 f'{pm["tokens_median"]:,.0f} | {vb["spearman_within"][k]:+.3f} | {vb["spearman_within_tokens"][k]:+.3f} |')
    L += ['', '### MC against reasoning tokens (descriptive)', '',
          'Quartiles of reasoning tokens within each model; MC as scored per quartile (item mean), the median reasoning '
          'tokens and the share of empty responses per quartile.', ''] + rq_lines(res['reasoning_tokens'], keys)
    if res.get('matched'):
        L += matched_lines(res)
    return '\n'.join(L) + '\n'


def ci4(c) -> str:
    return f'[{c[0]:+.4f}, {c[1]:+.4f}]'


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ------------------------------------------------------------------ self-test

def selftest() -> int:
    # 1. Every stored response: as scored equals e5_strict and coverage_value (variants() asserts both).
    ros = roster()
    ms_all = json.loads(score.MILESTONES.read_text(encoding='utf-8'))
    stores = [(m['default'], m['key']) for m in ros] + [(m['reasoning'], m['key']) for m in ros if m['reasoning']]
    targets = item_targets(stores)
    n = 0
    for store, key in stores:
        _rows, vals = store_values(store, key, ms_all, targets)
        n += sum(v['as_scored'] is not None for v in vals.values())
    print(f'as scored = e5_strict = coverage_value on {n:,} responses in {len(stores)} model stores')

    # 2. The four readings on a hand-made response: five milestones, the fifth the answer (3.0 m = 3000 mm).
    ms = [{'id': 'a', 'value': 1.5}, {'id': 'b', 'value': 2.5}, {'id': 'c', 'value': 7.0},
          {'id': 'd', 'value': 9.0}, {'id': 'e', 'value': 3000.0}]
    flags = target_flags(ms, (3.0,))
    assert flags == [False, False, False, False, True], flags
    r = {'item_id': 'x#0', 'milestones_required': 5, 'e3': {'coverage': 0.4, 'reached': [True, False, False, False, True]}}
    e5 = {'sources': ['e3', 'REACHED', 'NOT_NEEDED', 'MISSING', 'e3'], 'e5_strict': 0.6, 'sent': True, 'reply_ok': True}
    got = variants(r, e5, flags)
    want = {'as_scored': 3 / 5, 'matching_only': 2 / 5, 'route_adjusted': 3 / 4, 'intermediate_only': 2 / 4}
    assert all(abs(got[k] - want[k]) < 1e-12 for k in want), got
    e5u = {**e5, 'sources': ['e3', 'REACHED', 'NOT_NEEDED', 'UNJUDGED', 'e3']}
    assert variants(r, e5u, flags)['route_adjusted'] == 3 / 4          # UNJUDGED stays in the denominator
    r0 = {'item_id': 'x#1', 'milestones_required': 2, 'e3': {'coverage': 0.0, 'reached': [False, False]}}
    e50 = {'sources': ['NOT_NEEDED', 'NOT_NEEDED'], 'e5_strict': 0.0, 'sent': True, 'reply_ok': True}
    got0 = variants(r0, e50, [True, True])
    assert got0['as_scored'] == 0.0 and got0['route_adjusted'] is None and got0['intermediate_only'] is None, got0
    assert variants(r, {**e5, 'reply_ok': False}, flags) == dict.fromkeys(VARIANTS)      # no reply: no value (D-148)
    assert variants({**r, 'milestones_required': 0, 'e3': {'coverage': None, 'reached': []}}, None, []) == dict.fromkeys(VARIANTS)
    try:
        variants(r, {**e5, 'e5_strict': 0.4}, flags)
    except AssertionError:
        pass
    else:
        raise AssertionError('a stored e5_strict that disagrees with the sources passed')
    print('the four readings on the hand-made response: as scored 0.6, matching 0.4, route-adjusted 0.75, '
          'intermediate-only 0.5; all not needed and all targets give no value; no reply gives none')

    # 3. The reader: a thousands grouping and scientific notation are one number each; a subscript and an
    #    exponent are none: 13855, 2.1506e5 and 3.2.
    text = 'V = 13 855 m^3, p = 2.1506 \\times 10^{5} Pa, x_1 = 3.2^2'
    assert numbers_shown(text) == 3, (numbers_shown(text), arith.normalise(text))

    # 4. The fixed-effects fit recovers known slopes, and a resample's solve equals a direct fit on it.
    rng = np.random.default_rng(1)
    T, per = 40, 30
    tpl = np.repeat(np.arange(T), per)
    X = rng.normal(size=(T * per, 2)) + tpl[:, None] * 0.1
    y = rng.normal(size=T)[tpl] + X @ np.array([0.03, -0.02]) + rng.normal(scale=0.01, size=T * per)
    fit = fe_fit(y, X, tpl, tpl, T, 200, 3)
    assert np.allclose(fit['beta'], [0.03, -0.02], atol=2e-3), fit['beta']
    assert fit['ci'][0][0] < 0.03 < fit['ci'][0][1], fit['ci']
    pick = np.array([0, 0, 5, 7])
    rows_ = np.concatenate([np.flatnonzero(tpl == t) for t in pick])
    cell = np.concatenate([np.full(per, j) for j in range(len(pick))])
    direct = fe_fit(y[rows_], X[rows_], cell, cell, len(pick), 1, 0)['beta']
    Sxx, Sxy = fe_crossproducts(y, X, tpl, tpl, T)
    w = np.bincount(pick, minlength=T)
    assert np.allclose(np.linalg.solve(np.einsum('t,tij->ij', w, Sxx), np.einsum('t,ti->i', w, Sxy)), direct)
    print('the fixed-effects fit recovers (0.03, -0.02); a resample solved from cross-products equals a direct fit')

    # 5. The quartiles and the pair lookup.
    rows_q = {f'i{j}': {'item_id': f'i{j}', 'reasoning_tokens': j + 1, 'status': 'answered'} for j in range(8)}
    vals_q = {f'i{j}': {'as_scored': j / 7} for j in range(8)}
    rq = reasoning_quartiles(rows_q, vals_q)
    assert np.isclose(rq['q1'], 0.5 / 7) and np.isclose(rq['q4'], 6.5 / 7) and rq['tokens_median'] == [1.5, 3.5, 5.5, 7.5], rq
    assert reasoning_quartiles({k: {**v, 'reasoning_tokens': 0} for k, v in rows_q.items()}, vals_q) is None
    fam = [{'a': PAIR[1], 'b': PAIR[0], 'diff': -0.02, 'ci': [-0.03, -0.01], 'p_holm': 0.01}]
    assert find_pair(fam, *PAIR)['diff'] == 0.02 and find_pair(fam, *PAIR)['ci'] == [0.01, 0.03]
    print('quartiles and the pair lookup: ok')

    # 6. The matched family on hand-made data: three models over twelve templates, the third re-run with reasoning
    #    (higher coverage, reasoning tokens). Template t00 has no route-adjusted value anywhere, so that reading
    #    compares eleven. The two unchanged models keep their summaries and their pair; the third's change.
    def store6(shift: float, seed: int, think: bool) -> tuple[dict, dict]:
        g = np.random.default_rng(seed)
        rows, vals = {}, {}
        for t in range(12):
            for j in range(3):
                i, x = f't{t:02d}#{j}', float(np.clip(g.uniform(0.2, 0.8) + shift, 0, 1))
                rows[i] = {'item_id': i, 'template_id': f't{t:02d}', 'unusable': False, 'score': int(j == 0),
                           'status': 'answered', 'completion_tokens': 400 + 50 * j * j + (900 if think else 0),
                           'reasoning_tokens': 900 + 10 * t if think else 0}
                vals[i] = {'as_scored': x, 'matching_only': 0.8 * x, 'route_adjusted': x if t else None,
                           'intermediate_only': 0.9 * x}
        return rows, vals
    ros6 = [{'key': PAIR[0], 'default': 'main', 'reasoning': None},
            {'key': PAIR[1], 'default': 'main', 'reasoning': None},
            {'key': 'gpt-5.4-mini', 'default': 'main', 'reasoning': 'reasoning-medium-full'}]
    keys6 = [m['key'] for m in ros6]
    assert [matched_store(m) for m in ros6] == ['main', 'main', 'reasoning-medium-full']
    default6 = {PAIR[0]: store6(0.1, 1, True), PAIR[1]: store6(0.0, 2, True), 'gpt-5.4-mini': store6(-0.2, 3, False)}
    matched6 = {**default6, 'gpt-5.4-mini': store6(0.0, 4, True)}
    tpl6 = sorted({r['template_id'] for r in default6[PAIR[0]][0].values()})
    _C, use_d, sd, fd = family(keys6, default6, tpl6, 2_000)
    Cm, use_m, sm, fm = family(keys6, matched6, tpl6, 2_000)
    assert all(use_d[v].sum() == use_m[v].sum() == (11 if v == 'route_adjusted' else 12) for v in VARIANTS)
    assert all(jsonable(sm[k]) == jsonable(sd[k]) for k in PAIR), 'an unchanged model moved at matched settings'
    rr, rv = matched6['gpt-5.4-mini']
    direct = np.mean([np.mean([rv[i]['as_scored'] for i in rr if rr[i]['template_id'] == t]) for t in tpl6])
    assert np.isclose(sm['gpt-5.4-mini']['as_scored']['all'], direct)
    assert sm['gpt-5.4-mini']['as_scored']['all'] > sd['gpt-5.4-mini']['as_scored']['all'] + 0.1
    for v in VARIANTS:
        a, b = fd[v][0], fm[v][0]               # pair 0 is the two unchanged models
        assert (a['a'], a['b']) == PAIR and all(a[x] == b[x] for x in ('diff', 'ci', 'p')), (v, a, b)
        assert fd[v][1]['diff'] != fm[v][1]['diff']
    sep6 = separation(fm)['claude_vs_deepseek']
    assert all(sep6[v] == find_pair(fm[v], *PAIR)['p_holm'] for v in VARIANTS)
    texts6 = {k: {i: ', '.join(f'{"abcde"[n]} = {n + 1}.5 m' for n in range(int(i[-1]) + 1 + h)) for i in matched6[k][0]}
              for h, k in enumerate(keys6)}
    vd = {k: verbosity_rows(k, *default6[k], texts6[k]) for k in keys6}
    vm = {**vd, 'gpt-5.4-mini': verbosity_rows('gpt-5.4-mini', *matched6['gpt-5.4-mini'], texts6['gpt-5.4-mini'])}
    assert [d['numbers'] for d in vm['gpt-5.4-mini']] == [3, 4, 5] * 12 and vm['gpt-5.4-mini'][2]['tokens'] == 600
    vb_d = verbosity([d for k in keys6 for d in vd[k]], 200)
    vb_m = verbosity([d for k in keys6 for d in vm[k]], 200)
    assert all(vb_d['per_model'][k] == vb_m['per_model'][k] for k in PAIR) and vb_m['responses'] == 108
    assert vb_m['per_model']['gpt-5.4-mini']['mc'] > vb_d['per_model']['gpt-5.4-mini']['mc']
    assert floored(matched6, keys6) == dict.fromkeys(keys6, 0)
    assert reasoning_quartiles(*default6['gpt-5.4-mini']) is None and reasoning_quartiles(*matched6['gpt-5.4-mini'])
    print('the matched family on hand-made data: the unchanged models keep their summaries, their pair (diff, '
          'interval, p) and their verbosity rows; the re-run model\'s matrix row, pairs, fit and quartiles change')
    print('selftest passed')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--quick', action='store_true', help='1,000 resamples and 10,000 sign flips; the outputs say QUICK')
    ap.add_argument('--selftest', '--check', action='store_true', dest='selftest',
                    help='reproduce e5_strict for every stored response and the variants on a hand-made example')
    ap.add_argument('--render', action='store_true', help='rewrite the section from the JSON already written; no computing')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.render:
        MD_OUT.write_text(render(json.loads(JSON_OUT.read_text(encoding='utf-8'))), encoding='utf-8')
        print(f'wrote {MD_OUT.relative_to(REPO)} from {JSON_OUT.relative_to(REPO)}')
        return 0
    quick = a.quick or os.environ.get('ENGTRACE_QUICK') == '1'
    res = jsonable(run(quick))
    OUT.mkdir(exist_ok=True)
    MD_OUT.parent.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    MD_OUT.write_text(render(res), encoding='utf-8')
    sep = res['separation_rule']['claude_vs_deepseek']
    print(f'wrote {JSON_OUT.relative_to(REPO)} and {MD_OUT.relative_to(REPO)}{" (QUICK)" if quick else ""}')
    if sep:
        print('Claude Sonnet 5 vs DeepSeek V4.1 Flash, p Holm: '
              + ', '.join(f'{v} {sep[v]:.4f}' for v in VARIANTS) + f'; separation holds: {sep["holds"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
