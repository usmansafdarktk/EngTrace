"""The analysis plan (ANALYSIS_PLAN.md, D-117), computed from the score store.

    python -m full_run_28092026.analyze              # writes results/RESULTS.md, results.json and per_template.csv
    python -m full_run_28092026.analyze --selftest   # the statistics on synthetic data; writes nothing

Reads scores/main/<model>.jsonl (score.py) for the eleven roster models, and scores/paraphrase/ and
scores/repeat1..3/ once those runs exist. The set-aside Qwen3.8-27B (D-132) is reported on its own
line, outside every comparison. Nothing here reads a trace or calls a model.

WHAT THE PLAN FIXES. The unit (the template), the intervals (template bootstrap, B = 10,000,
percentile), the scoring (correct 1, partial 0.5, incorrect and unusable 0), the tests and Holm's
correction within each family. The script stops unless the store holds what the plan describes:
150 templates x 15 items per model, 58 Easy and 34 Advanced templates, 58 single-path ones.

WHAT THIS SCRIPT SETS WHERE THE PLAN IS SILENT, written before any result was read (D-136):
  - Q2's p-value: a permutation test of the tier labels among the 92 Easy and Advanced templates,
    beside the plan's within-tier bootstrap interval.
  - Q3's "wrong-answer traces" are those scoring 0 (incorrect or unusable), its "correct-answer
    traces" those fully solved; a partial answer is in neither. A trace counts as flagged by the digit
    rule when any of its steps is. Over a subset of traces the estimate is the item mean, its interval
    the template bootstrap of that ratio.
  - Q1's fully-solved check "gives the same verdict" on a pair when both tests or neither hold at
    Holm-adjusted 0.05, and when both hold, in the same direction.
  - Q5's paired difference is the item mean; the sign flips act on each template's summed
    difference, so the test statistic is that same mean. Kendall's tau's interval resamples templates.
  - Every p-value from resampling is (count + 1) / (draws + 1), so none is reported as 0.
  - Each test draws from its own seeded generator, so a re-run prints the same numbers.
Per-label accuracy on classification templates takes the label from the gold's answer targets; the
Froude template's answer is the Froude number, its label the regime it implies (D-118).

CORRECTED AND ADDED AFTER THE FIRST RESULTS WERE READ, each labelled where it appears (D-146, D-148,
D-149; the two independent reviews of 2026-09-29):
  - Q2. The tier-label permutation of the raw difference is liberal when the smaller tier has the
    larger spread, and Advanced template means spread two to four times as widely as Easy ones; its
    count of models also moved with the seed (7, 7, 5 of 11). A Welch t-test, valid under unequal
    variances, is the test the count now rests on; the planned permutation is printed beside it. The
    resampling tests draw B_TEST = 100,000 permutations, so the Holm floor over 55 pairs is 0.0006
    rather than 0.0055, and the Monte Carlo standard error of a p-value near 0.05 is 0.0007. The
    detectable gap is also given at the strictest Holm step from the Welch standard error.
  - Q2 is also shown with unusable rows left out of the template means (glm-5.3's gap is an output-
    ceiling effect: 100 of its 123 wrong answers are empty rows) and without the symbolic templates.
  - Q1 adds a template-level sign-flip test on the fully-solved rate beside McNemar's item-level one,
    the per-template SD the July rebuttal promised (within a template, over its 15 items, and between
    template means), and the split of "unusable" into empty and unreadable.
  - Q3's digit-rule rate on wrong answers is over the answered wrong answers: an empty trace has no
    step to flag. It adds the 1% rate beside the digit rule, E3's chance floor per model (the trace
    against a sibling item's milestones), the position of the first flagged step, and the wrong-answer
    rate against the item's milestone count.
  - A stage with calls left unanswered is incomplete: its rates are computed over the traces that
    were answered, the count without a reply is printed beside them, and the header says so.
  - Sensitivity adds three readings of the answer check that the experts could not arbitrate: the
    half-unit window, the whole trace instead of the Answer line (D-147), and the pool without the
    nine symbolic templates (D-138); and a noise floor for Kendall's tau (the ordering on one half of
    each template's items against the other).
  - Provider: the score by serving endpoint, matched on template, because dispatch order confounds a
    raw per-provider mean with the templates each endpoint happened to serve (D-133).
  - Q5 adds the paired E3 coverage difference (and E5-strict, once E5 has run on both arms) and the
    same-provider pairs alone.
  - results/per_template.csv: one row per template and model, aggregates only, so figures and a
    reader's own bootstrap can be redone without the private store.
  - results.json records what produced it: this script's commit and LF hash, the store's CONFIG and
    every stage CONFIG, checked against the store they were built from.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

SCORES = HERE / 'scores'
OUT = HERE / 'results'
B = 10_000              # bootstrap draws for every interval (the plan)
B_TEST = 100_000        # permutations / sign flips for every test (D-146; D-136 had 10,000)
EQUIV_MARGIN = 0.05     # Q5's equivalence margin in answer-score points: the plan's largest detectable paired
                        # difference (5.1 points at 15% discordance), fixed 2026-10-01 after the point estimates
                        # were known and before the 90% intervals were computed (D-163)
ROSTER = ['gpt-oss-20b', 'gemma-4-26b-a4b', 'deepseek-v4.1-flash', 'qwen3-235b-a22b-2507', 'glm-5.3-flash',
          'glm-5.3', 'muse-glimmer-30b', 'kimi-k3', 'gpt-5.4-mini', 'gemini-3.1-flash-lite', 'claude-sonnet-5']
SET_ASIDE = ['qwen3.8-27b']
# the four templates the pilot excluded for shortcuts: D-057 (two), D-046, D-066
SHORTCUT = ['template_system_property_linearity', 'template_system_properties_memory_causality',
            'template_line_balancing_heuristic', 'template_levenspiel_plot_interpretation']
FROUDE = 'template_critical_depth_froude_classification'
REPEATS = ['repeat1', 'repeat2', 'repeat3']
PLAN_COUNTS = {'templates': 150, 'per_template': 15, 'Easy': 58, 'Advanced': 34, 'single_path': 58}
MS_BUCKETS = (('0', 0, 0), ('1', 1, 1), ('2', 2, 2), ('3', 3, 3), ('4-5', 4, 5), ('6+', 6, 10 ** 6))
Z80 = float(stats.norm.ppf(0.8))


# ------------------------------------------------------------------ loading

def load(variant: str, key: str) -> dict[str, dict] | None:
    p = SCORES / variant / f'{key}.jsonl'
    if not p.exists():
        return None
    return {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())}


def load_stage(stage: str, key: str, store: str = 'main') -> dict[str, dict] | None:
    """judge.py's E5 rows or router.py's rows for one model, when that stage has run."""
    p = SCORES / store / stage / f'{key}.jsonl'
    if not p.exists():
        return None
    return {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())}


def read_json(p: Path):
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


def fully(r: dict) -> float:
    return 1.0 if (not r['unusable'] and r['answer']['label'] == 'correct') else 0.0


def score_with(r: dict, field: str) -> float | None:
    """The score under another answer label; None when the store predates that label."""
    if r['unusable']:
        return 0.0
    lab = r['answer'].get(field)
    return None if lab is None else {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}[lab]


def flagged(r: dict) -> float:
    return 1.0 if any(s['digit_flags'] for s in r['steps']) else 0.0


def flagged_tol1(r: dict) -> float:
    return 1.0 if any(s.get('tol1_flags') for s in r['steps']) else 0.0


def verdict3(r: dict) -> str:
    return 'unusable' if r['unusable'] else r['answer']['label']


def per_template(rows, fn, keep=lambda r: True) -> dict[str, list[float]]:
    out = collections.defaultdict(list)
    for r in rows:
        if keep(r):
            out[r['template_id']].append(fn(r))
    return out


def template_matrix(runs: dict, fn, templates: list[str], keep=lambda r: True) -> np.ndarray:
    """models x templates: the mean of fn(row) over each template's kept items; NaN when none kept."""
    M = np.full((len(runs), len(templates)), np.nan)
    for i, rows in enumerate(runs.values()):
        by_t = per_template(rows.values(), fn, keep)
        M[i] = [np.mean(by_t[t]) if by_t.get(t) else np.nan for t in templates]
    return M


def check_store(runs: dict) -> tuple[list[str], dict, dict]:
    """The store must hold what the plan describes, the same for every model."""
    ref = next(iter(runs.values()))
    templates = sorted({r['template_id'] for r in ref.values()})
    levels = {r['template_id']: r['level'] for r in ref.values()}
    single = {r['template_id']: r['single_path'] for r in ref.values()}
    tier = collections.Counter(levels.values())
    got = {'templates': len(templates), 'Easy': tier['Easy'], 'Advanced': tier['Advanced'],
           'single_path': sum(single.values()),
           'per_template': min(collections.Counter(r['template_id'] for r in ref.values()).values())}
    bad = {k: (v, got[k]) for k, v in PLAN_COUNTS.items() if got[k] != v}
    for k, rows in runs.items():
        if set(rows) != set(ref):
            bad[k] = 'holds a different item set'
    missing = [t for t in SHORTCUT + [FROUDE] if t not in levels]
    if bad or missing:
        raise SystemExit(f'the store is not the pool the plan describes: {bad} {missing}')
    return templates, levels, single


# ------------------------------------------------------------------ statistics

def pct(a) -> list[float]:
    return [float(np.nanpercentile(a, 2.5)), float(np.nanpercentile(a, 97.5))]


def boot_mean(v: np.ndarray, seed: int) -> list[float]:
    """Percentile interval of a mean over templates, templates resampled."""
    rng = np.random.default_rng(seed)
    v = v[~np.isnan(v)]
    return pct(v[rng.integers(0, len(v), size=(B, len(v)))].mean(axis=1))


def cluster_mean(groups: dict[str, list[float]], seed: int) -> tuple[float, list[float], int]:
    """The item mean over a subset of traces, with the template bootstrap of that ratio; the
    resampled units are the templates with a qualifying trace."""
    if not groups:
        return float('nan'), [float('nan'), float('nan')], 0
    tot = np.array([sum(v) for v in groups.values()])
    cnt = np.array([len(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(tot), size=(B, len(tot)))
    return float(tot.sum() / cnt.sum()), pct(tot[idx].sum(axis=1) / cnt[idx].sum(axis=1)), int(cnt.sum())


def cluster_boot(groups: dict[str, list[float]], seed: int) -> tuple[float, np.ndarray, int]:
    """cluster_mean's point estimate with its bootstrap draws (the same draws, so the same 95% interval),
    for an interval at another level as well (D-163)."""
    if not groups:
        return float('nan'), np.array([]), 0
    tot = np.array([sum(v) for v in groups.values()])
    cnt = np.array([len(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(tot), size=(B, len(tot)))
    return float(tot.sum() / cnt.sum()), tot[idx].sum(axis=1) / cnt[idx].sum(axis=1), int(cnt.sum())


def sign_flip_p(d: np.ndarray, seed: int, draws: int = B_TEST) -> float:
    rng = np.random.default_rng(seed)
    d = d[~np.isnan(d)]
    obs = abs(d.mean())
    null = np.abs((rng.choice([-1.0, 1.0], size=(draws, len(d))) * d).mean(axis=1))
    return float((np.sum(null >= obs - 1e-12) + 1) / (draws + 1))


def holm(ps: list[float]) -> list[float]:
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, (len(ps) - rank) * ps[i])
        adj[i] = min(1.0, run)
    return adj


def mcnemar_p(a: np.ndarray, b: np.ndarray) -> float:
    """Exact McNemar on paired 0/1 verdicts."""
    n01 = int(np.sum((a == 1) & (b == 0)))
    n10 = int(np.sum((a == 0) & (b == 1)))
    return 1.0 if n01 + n10 == 0 else float(stats.binomtest(n01, n01 + n10, 0.5).pvalue)


def tier_perm_p(x: np.ndarray, is_easy: np.ndarray, seed: int, draws: int = B_TEST) -> float:
    """Easy minus Advanced, the tier labels shuffled among the templates of both tiers (the planned
    test, D-136; liberal under unequal variances, D-146)."""
    rng = np.random.default_rng(seed)
    obs = abs(x[is_easy].mean() - x[~is_easy].mean())
    n_e = int(is_easy.sum())
    null = np.empty(draws)
    for start in range(0, draws, 10_000):                       # in slices, to bound memory
        n = min(10_000, draws - start)
        px = x[np.argsort(rng.random((n, len(x))), axis=1)]
        null[start:start + n] = np.abs(px[:, :n_e].mean(axis=1) - px[:, n_e:].mean(axis=1))
    return float((np.sum(null >= obs - 1e-12) + 1) / (draws + 1))


def welch(xe: np.ndarray, xa: np.ndarray) -> dict:
    """Welch's t-test, two-sided, and the gap detectable at 80% power: at two-sided 0.05 (the plan's
    2.80 sigma) and at the strictest Holm step, 0.05 / 11, both from the Welch standard error."""
    xe, xa = xe[~np.isnan(xe)], xa[~np.isnan(xa)]
    se = float(np.sqrt(xe.var(ddof=1) / len(xe) + xa.var(ddof=1) / len(xa)))
    p = float(stats.ttest_ind(xe, xa, equal_var=False).pvalue)
    z_plain, z_holm = float(stats.norm.ppf(1 - 0.05 / 2)), float(stats.norm.ppf(1 - 0.05 / (2 * len(ROSTER))))
    return {'p_welch': p, 'se_welch': se, 'sd_easy': float(xe.std(ddof=1)), 'sd_advanced': float(xa.std(ddof=1)),
            'detectable_welch': (z_plain + Z80) * se, 'detectable_holm': (z_holm + Z80) * se}


def boot_gap(xe: np.ndarray, xa: np.ndarray, seed: int) -> list[float]:
    """Easy minus Advanced, templates resampled within each tier."""
    rng = np.random.default_rng(seed)
    xe, xa = xe[~np.isnan(xe)], xa[~np.isnan(xa)]
    me = xe[rng.integers(0, len(xe), size=(B, len(xe)))].mean(axis=1)
    ma = xa[rng.integers(0, len(xa), size=(B, len(xa)))].mean(axis=1)
    return pct(me - ma)


def kendall_boot(S: np.ndarray, N: np.ndarray, seed: int) -> dict:
    """Kendall's tau between the models' scores in two arms, templates resampled.
    S: arms x models x templates of summed scores; N: the item counts, the same shape."""
    def tau(j):
        a, p = S[0][:, j].sum(axis=1) / N[0][:, j].sum(axis=1), S[1][:, j].sum(axis=1) / N[1][:, j].sum(axis=1)
        return stats.kendalltau(a, p).statistic
    rng = np.random.default_rng(seed)
    T = S.shape[2]
    taus = np.array([tau(rng.integers(0, T, size=T)) for _ in range(B)])
    return {'tau': float(tau(np.arange(T))), 'ci': [float(np.nanpercentile(taus, 2.5)),
                                                    float(np.nanpercentile(taus, 97.5))]}


# ------------------------------------------------------------------ the questions

def q1(runs, templates, keys):
    S = template_matrix(runs, lambda r: r['score'], templates)
    F = template_matrix(runs, fully, templates)
    models = []
    for i, k in enumerate(keys):
        rows = list(runs[k].values())
        split = collections.Counter(verdict3(r) for r in rows)
        within = [np.std(v, ddof=1) for v in per_template(rows, lambda r: r['score']).values() if len(v) > 1]
        capped = [r for r in rows if r['status'] == 'answered' and r.get('finish_reason') == 'length']
        odd = [r for r in rows if r['status'] == 'answered' and r.get('finish_reason') not in ('stop', 'length')]
        models.append({'model': k, 'score': float(S[i].mean()), 'ci': boot_mean(S[i], 100 + i),
                       'within_template_sd': float(np.mean(within)), 'between_template_sd': float(S[i].std(ddof=1)),
                       'fully_solved': float(F[i].mean()), 'fully_ci': boot_mean(F[i], 200 + i),
                       **{v: split[v] for v in ('correct', 'partial', 'incorrect', 'unusable')},
                       'empty': sum(r['status'] != 'answered' for r in rows),
                       'unreadable': sum(r['status'] == 'answered' and r['unusable'] for r in rows),
                       'capped_scored': len(capped),
                       'capped_verdicts': dict(collections.Counter(verdict3(r) for r in capped)),
                       'odd_finish_scored': len(odd)})
    ids = sorted(runs[keys[0]])
    V = {k: np.array([fully(runs[k][x]) for x in ids]) for k in keys}
    pairs = []
    for n, (i, j) in enumerate(itertools.combinations(range(len(keys)), 2)):
        d, df = S[i] - S[j], F[i] - F[j]
        pairs.append({'a': keys[i], 'b': keys[j], 'diff': float(d.mean()), 'ci': boot_mean(d, 1000 + n),
                      'p': sign_flip_p(d, 2000 + n), 'fully_diff': float(df.mean()),
                      'fully_ci': boot_mean(df, 3000 + n), 'fully_p': sign_flip_p(df, 4000 + n),
                      'mcnemar_p': mcnemar_p(V[keys[i]], V[keys[j]])})
    for o, ph, fh, qh in zip(pairs, holm([p['p'] for p in pairs]), holm([p['fully_p'] for p in pairs]),
                             holm([p['mcnemar_p'] for p in pairs])):
        o['p_holm'], o['fully_p_holm'], o['mcnemar_p_holm'] = ph, fh, qh
        s1, s2, s3 = ph < 0.05, fh < 0.05, qh < 0.05
        agree = lambda t: bool((s1 == t) and (not s1 or np.sign(o['diff']) == np.sign(o['fully_diff'])))
        o['same_verdict'] = agree(s3)                      # answer score vs McNemar, as D-136 defined it
        o['same_verdict_template'] = agree(s2)             # answer score vs the template-level fully-solved test
    return {'models': models, 'pairs': pairs}


def q2(runs, templates, levels, keys, symbolic: set[str]):
    tiers = np.array([levels[t] for t in templates])
    both = (tiers == 'Easy') | (tiers == 'Advanced')
    variants = {'as_scored': (lambda r: True, set(templates)),
                'unusable_excluded': (lambda r: not r['unusable'], set(templates)),
                'without_symbolic': (lambda r: True, set(templates) - symbolic)}
    out = []
    for i, k in enumerate(keys):
        row = {'model': k}
        for name, (keep, tset) in variants.items():
            S = template_matrix({k: runs[k]}, lambda r: r['score'], templates, keep)[0]
            use = np.array([t in tset for t in templates])
            xe, xa = S[(tiers == 'Easy') & use], S[(tiers == 'Advanced') & use]
            w = welch(xe, xa)
            res = {'easy': float(np.nanmean(xe)), 'advanced': float(np.nanmean(xa)),
                   'gap': float(np.nanmean(xe) - np.nanmean(xa)), 'ci': boot_gap(xe, xa, 400 + i), **w,
                   'templates': int(use.sum())}
            if name == 'as_scored':
                xe_, xa_ = xe[~np.isnan(xe)], xa[~np.isnan(xa)]
                sigma = float(np.sqrt(((len(xe_) - 1) * xe_.var(ddof=1) + (len(xa_) - 1) * xa_.var(ddof=1))
                                      / (len(xe_) + len(xa_) - 2)))
                res.update({'p_perm': tier_perm_p(S[both], tiers[both] == 'Easy', 500 + i), 'sigma': sigma,
                            'detectable_planned': 2.80 * sigma * float(np.sqrt(1 / len(xe_) + 1 / len(xa_)))})
                row.update(res)
            else:
                row[name] = res
        out.append(row)
    for o, ph in zip(out, holm([o['p_perm'] for o in out])):
        o['p_perm_holm'] = ph
    for o, ph in zip(out, holm([o['p_welch'] for o in out])):
        o['p_welch_holm'] = ph
    for name in ('unusable_excluded', 'without_symbolic'):
        for o, ph in zip(out, holm([o[name]['p_welch'] for o in out])):
            o[name]['p_welch_holm'] = ph
    return out


def stage_q3(rows: list[dict], e5: dict | None, router: dict | None, seed: int) -> dict:
    """Q3's judged columns for one model: E5 on the wrong-answer traces, its judged and unjudged
    fractions, the router's flags, and which component points at a wrong answer (also reported).
    A trace whose call got no reply is left out of every rate and counted (D-148)."""
    wrong = lambda r: r['score'] == 0.0
    out = {}
    if e5:
        ok = lambda r: not (e5[r['item_id']]['sent'] and e5[r['item_id']]['reply_ok'] is False)
        answered = [e5[r['item_id']] for r in rows if r['status'] == 'answered' and e5[r['item_id']]['milestones_required']]
        sent = [x for x in answered if x['sent']]
        replied = [x for x in sent if x['reply_ok']]
        judged = sum(x['milestones_required'] - x['by_e3'] for x in replied)
        cov, cov_ci, n = cluster_mean(per_template(rows, lambda r: e5[r['item_id']]['e5_strict'],
                                                   lambda r: wrong(r) and ok(r) and e5[r['item_id']]['e5_strict'] is not None),
                                      seed)
        rcov, rcov_ci, _n = cluster_mean(per_template(rows, lambda r: e5[r['item_id']]['e5_strict'],
                                                      lambda r: wrong(r) and not r['unusable'] and ok(r)
                                                      and e5[r['item_id']]['e5_strict'] is not None), seed + 1)
        out.update({'e5_coverage_on_wrong': cov, 'e5_ci': cov_ci, 'e5_coverage_on_readable_wrong': rcov,
                    'e5_readable_ci': rcov_ci, 'e5_wrong_used': n,
                    'e5_judged_fraction': judged / max(1, sum(x['milestones_required'] for x in answered
                                                             if x['reply_ok'] is not False)),
                    'e5_unjudged_rate': sum(x['unjudged'] for x in replied) / max(1, judged),
                    'e5_reached_share': sum(x['reached'] for x in replied) / max(1, judged),
                    'e5_calls': len(sent), 'e5_without_reply': len(sent) - len(replied),
                    'e5_incomplete': len(sent) - len(replied) > 0})
    if router:
        ok_r = lambda r: router[r['item_id']]['reply_ok'] is not False
        flag_any = lambda r: 1.0 if router[r['item_id']]['router_flagged'] else 0.0
        judge_any = lambda r: 1.0 if router[r['item_id']]['judge_flagged'] else 0.0
        answered = [r for r in rows if r['status'] == 'answered']
        replied = [r for r in answered if ok_r(r)]
        fs, fs_ci, _ = cluster_mean(per_template(rows, flag_any, lambda r: fully(r) == 1 and ok_r(r)), seed + 2)
        js, js_ci, _ = cluster_mean(per_template(rows, judge_any, lambda r: fully(r) == 1 and ok_r(r)), seed + 3)
        ws, ws_ci, _ = cluster_mean(per_template(rows, flag_any, lambda r: wrong(r) and r['status'] == 'answered' and ok_r(r)),
                                    seed + 4)
        sent = sum(len(router[r['item_id']]['sent']) for r in replied)
        out.update({'router_rate_on_fully_solved': fs, 'router_ci': fs_ci,
                    'router_judge_rate_on_fully_solved': js, 'router_judge_ci': js_ci,
                    'router_rate_on_wrong': ws, 'router_wrong_ci': ws_ci,
                    'router_steps_flagged_per_trace': float(np.mean([len(router[r['item_id']]['router_flagged'])
                                                                     for r in replied])) if replied else float('nan'),
                    'router_unjudged_rate': sum(router[r['item_id']]['unjudged'] for r in replied) / max(1, sent),
                    'router_calls': sum(bool(router[r['item_id']]['sent']) for r in answered),
                    'router_without_reply': sum(router[r['item_id']]['reply_ok'] is False for r in answered),
                    'router_incomplete': any(router[r['item_id']]['reply_ok'] is False for r in answered)})
    if e5 or router:
        w = [r for r in rows if wrong(r) and r['status'] == 'answered'
             and (not e5 or e5[r['item_id']]['reply_ok'] is not False) and (not router or router[r['item_id']]['reply_ok'] is not False)]
        out['attribution_on_wrong'] = {
            'traces': len(w), 'digit_rule': float(np.mean([flagged(r) for r in w])) if w else None,
            'e5_missing': (float(np.mean([e5[r['item_id']]['missing'] > 0 for r in w])) if w else None) if e5 else None,
            'router_judge': (float(np.mean([bool(router[r['item_id']]['judge_flagged']) for r in w])) if w else None)
            if router else None}
    return out


def first_flag_position(rows) -> dict:
    """Among fully solved traces the digit rule flags: the index of the first flagged step as a
    fraction of the trace's steps (0 = first step), its median and quartiles."""
    pos = []
    for r in rows:
        if fully(r) != 1 or not r['steps']:
            continue
        idx = next((i for i, s in enumerate(r['steps']) if s['digit_flags']), None)
        if idx is not None:
            pos.append(idx / max(1, len(r['steps']) - 1) if len(r['steps']) > 1 else 0.0)
    return {'traces': len(pos), 'median': float(np.median(pos)) if pos else None,
            'q1': float(np.percentile(pos, 25)) if pos else None, 'q3': float(np.percentile(pos, 75)) if pos else None}


def by_milestone_count(rows) -> dict:
    """The wrong-answer rate of a model against the item's milestone count, in buckets."""
    out = {}
    for name, lo, hi in MS_BUCKETS:
        sel = [r for r in rows if lo <= r['milestones_required'] <= hi]
        out[name] = {'items': len(sel), 'wrong_rate': float(np.mean([r['score'] == 0.0 for r in sel])) if sel else None}
    return out


def q3(runs, keys, store='main'):
    wrong = lambda r: r['score'] == 0.0
    wrong_answered = lambda r: wrong(r) and r['status'] == 'answered'
    out = []
    for i, k in enumerate(keys):
        rows = list(runs[k].values())
        cov, cov_ci, n_cov = cluster_mean(per_template(rows, lambda r: r['e3']['coverage'],
                                                       lambda r: wrong(r) and r['e3']['coverage'] is not None),
                                          600 + i)
        rcov, rcov_ci, n_rcov = cluster_mean(per_template(rows, lambda r: r['e3']['coverage'],
                                                          lambda r: wrong(r) and not r['unusable']
                                                          and r['e3']['coverage'] is not None), 650 + i)
        fc, fc_ci, n_right = cluster_mean(per_template(rows, flagged, lambda r: fully(r) == 1), 700 + i)
        fw, fw_ci, n_wrong_a = cluster_mean(per_template(rows, flagged, wrong_answered), 750 + i)
        t1, t1_ci, _ = cluster_mean(per_template(rows, flagged_tol1, lambda r: fully(r) == 1), 760 + i)
        null_w = [r['e3'].get('null_coverage') for r in rows if wrong_answered(r) and r['e3'].get('null_coverage') is not None]
        null_c = [r['e3'].get('null_coverage') for r in rows if fully(r) == 1 and r['e3'].get('null_coverage') is not None]
        answered = [r for r in rows if r['status'] == 'answered']
        claims = [sum(s['claims'] for s in r['steps']) for r in answered]
        out.append({'model': k, 'claims_per_trace': float(np.mean(claims)),
                    'traces_with_a_claim': float(np.mean([c > 0 for c in claims])),
                    'wrong': sum(wrong(r) for r in rows), 'wrong_unusable': sum(wrong(r) and r['unusable'] for r in rows),
                    'wrong_with_milestones': n_cov, 'e3_coverage_on_wrong': cov, 'e3_ci': cov_ci,
                    'readable_wrong_with_milestones': n_rcov, 'e3_coverage_on_readable_wrong': rcov,
                    'e3_readable_ci': rcov_ci,
                    'e3_null_on_readable_wrong': float(np.mean(null_w)) if null_w else None,
                    'e3_null_on_fully_solved': float(np.mean(null_c)) if null_c else None,
                    'fully_solved': n_right, 'digit_flag_rate_on_fully_solved': fc, 'digit_ci': fc_ci,
                    'tol1_flag_rate_on_fully_solved': t1, 'tol1_ci': t1_ci,
                    'wrong_answered': n_wrong_a, 'digit_flag_rate_on_wrong': fw, 'digit_wrong_ci': fw_ci,
                    'first_flag_position': first_flag_position(rows), 'by_milestone_count': by_milestone_count(rows),
                    **stage_q3(rows, load_stage('e5', k, store), load_stage('router', k, store), 1100 + 10 * i)})
    return out


def q4(runs, templates, single, keys):
    out = []
    for i, k in enumerate(keys):
        by_t = per_template(runs[k].values(), fully)
        row = {'model': k}
        for g, (name, group) in enumerate((('single_path', [t for t in templates if single[t]]),
                                           ('multi_path', [t for t in templates if not single[t]]))):
            n = np.array([sum(by_t[t]) for t in group])
            row[name] = {'templates': len(group)}
            for c, (cat, mask) in enumerate((('all', n == 15), ('some', (n > 0) & (n < 15)), ('none', n == 0))):
                v = mask.astype(float)
                row[name][cat] = float(v.mean())
                row[name][cat + '_ci'] = boot_mean(v, 800 + 10 * i + 3 * g + c)
        out.append(row)
    return out


def accepted_pairs() -> dict | None:
    """The experts' check (paraphrase_kit.py --score): the pairs to keep, which are every assigned pair
    the experts did not reject, with the counts returned, rejected and outstanding; None before any
    return has been scored."""
    p = HERE / 'paraphrase' / 'accepted.json'
    if not p.exists():
        return None
    acc = json.loads(p.read_text(encoding='utf-8'))
    return {'keep': {i for i, r in acc.items() if r.get('kept') is not False},
            'returned': sum(1 for r in acc.values() if r.get('returned', True)),
            'rejected': sum(1 for r in acc.values() if r.get('kept') is False),
            'outstanding': sum(1 for r in acc.values() if not r.get('returned', True))}


def paired(main_rows, para_rows, ids, fn, seed):
    """Second arm minus main, paired by item: the item mean with its template bootstrap at 95% and, from the
    same draws, at 90%, which is the two-one-sided-tests reading against EQUIV_MARGIN (D-163)."""
    diff = per_template([{'template_id': main_rows[x]['template_id'], 'd': fn(para_rows[x]) - fn(main_rows[x])}
                         for x in ids], lambda r: r['d'])
    point, draws, n = cluster_boot(diff, seed)
    nan = [float('nan'), float('nan')]
    c = pct(draws) if n else nan
    c90 = [float(np.nanpercentile(draws, 5)), float(np.nanpercentile(draws, 95))] if n else nan
    return {'items': n, 'templates': len(diff), 'diff': point, 'ci': c, 'ci90': c90,
            'within_margin': bool(n and -EQUIV_MARGIN <= c90[0] and c90[1] <= EQUIV_MARGIN),
            'p': sign_flip_p(np.array([sum(v) for v in diff.values()]), seed + 50)}


def q5(main, para, keys, keep=None, e5_main=None, e5_para=None):
    """Paired paraphrase minus original, per model, on the items both arms hold; a pair the expert
    rejected leaves both arms (`keep`, when the check has returned). Besides the answer score: the
    E3 coverage difference on items with milestones, E5-strict once both arms carry it, and the
    answer score on the pairs both arms served from the same provider."""
    out, tested = [], []
    for i, k in enumerate(keys):
        if not para.get(k):
            continue
        ids = sorted(set(para[k]) & set(main[k]) & (keep if keep is not None else set(para[k])))
        row = {'model': k, **paired(main[k], para[k], ids, lambda r: r['score'], 900 + i),
               'mcnemar_p': mcnemar_p(np.array([fully(para[k][x]) for x in ids]),
                                      np.array([fully(main[k][x]) for x in ids]))}
        ms_ids = [x for x in ids if main[k][x]['e3']['coverage'] is not None and para[k][x]['e3']['coverage'] is not None]
        row['e3'] = paired(main[k], para[k], ms_ids, lambda r: r['e3']['coverage'], 1900 + i) if ms_ids else None
        if e5_main and e5_para and e5_main.get(k) and e5_para.get(k):
            e5_ids = [x for x in ms_ids if e5_main[k].get(x, {}).get('e5_strict') is not None
                      and e5_para[k].get(x, {}).get('e5_strict') is not None
                      and e5_main[k][x]['reply_ok'] is not False and e5_para[k][x]['reply_ok'] is not False]
            m5 = {x: {'template_id': main[k][x]['template_id'], 'v': e5_main[k][x]['e5_strict']} for x in e5_ids}
            p5 = {x: {'template_id': main[k][x]['template_id'], 'v': e5_para[k][x]['e5_strict']} for x in e5_ids}
            row['e5'] = paired(m5, p5, e5_ids, lambda r: r['v'], 2900 + i) if e5_ids else None
        else:
            row['e5'] = None
        same = [x for x in ids if main[k][x].get('provider') and main[k][x].get('provider') == para[k][x].get('provider')]
        row['same_provider'] = paired(main[k], para[k], same, lambda r: r['score'], 3900 + i) if same else None
        out.append(row)
        tested.append(k)
    for o, ph, qh in zip(out, holm([o['p'] for o in out]), holm([o['mcnemar_p'] for o in out])):
        o['p_holm'], o['mcnemar_p_holm'] = ph, qh
    tau = None
    if len(tested) >= 3:
        common = set.intersection(*(set(para[k]) & set(main[k]) for k in tested))
        if keep is not None:                 # a rejected pair leaves both arms here too (D-162)
            common &= keep
        common = sorted(common)
        ts = sorted({main[tested[0]][x]['template_id'] for x in common})
        col = {t: j for j, t in enumerate(ts)}
        S = np.zeros((2, len(tested), len(ts)))
        N = np.zeros((2, len(tested), len(ts)))
        for m, k in enumerate(tested):
            for x in common:
                j = col[main[k][x]['template_id']]
                S[0, m, j] += main[k][x]['score']
                S[1, m, j] += para[k][x]['score']
                N[:, m, j] += 1
        tau = {**kendall_boot(S, N, 990), 'items': len(common), 'noise_arm': tau_noise_arm(main, tested, common)}
    return {'models': out, 'tau': tau, 'expert_check': None if keep is None else len(keep),
            'margin': EQUIV_MARGIN, 'vs_repeats': vs_repeats(main, para, tested, keep, out)}


def tau_noise_arm(main, keys, ids, splits: int = 200, seed: int = 4243) -> dict:
    """The tau noise floor at the paraphrase arm's own size and template mix: for each template with k
    pairs in the arm, 2k of the main run's items drawn without replacement and split into two halves of
    k (k halved where a template has fewer than 2k items); the models ordered on each half by the
    item-weighted mean, as kendall_boot orders the two arms; tau between the two orderings (D-163)."""
    rng = np.random.default_rng(seed)
    ref = main[keys[0]]
    k_by_t = collections.Counter(ref[x]['template_id'] for x in ids)
    items_by_t = collections.defaultdict(list)
    for x in sorted(ref):
        items_by_t[ref[x]['template_id']].append(x)
    taus, n_half = [], 0
    for _ in range(splits):
        a_ids, b_ids = [], []
        for t, k in sorted(k_by_t.items()):
            pool = items_by_t[t]
            k = min(k, len(pool) // 2)
            if k == 0:
                continue
            pick = rng.choice(len(pool), size=2 * k, replace=False)
            a_ids += [pool[i] for i in pick[:k]]
            b_ids += [pool[i] for i in pick[k:]]
        n_half = len(a_ids)
        a = [np.mean([main[m][x]['score'] for x in a_ids]) for m in keys]
        b = [np.mean([main[m][x]['score'] for x in b_ids]) for m in keys]
        taus.append(stats.kendalltau(a, b).statistic)
    taus = np.array(taus)
    return {'splits': splits, 'items_per_half': n_half, 'templates': len(k_by_t), 'median': float(np.median(taus)),
            'q1': float(np.percentile(taus, 25)), 'q3': float(np.percentile(taus, 75)),
            'p5': float(np.percentile(taus, 5)), 'p95': float(np.percentile(taus, 95))}


def vs_repeats(main, para, keys, keep, q5_rows) -> dict:
    """The paraphrase difference beside the run-to-run differences of the models with decoding repeats:
    repeat minus main on the repeat items, paired as Q5 pairs; and both on the kept pairs the two
    subsamples share (D-163)."""
    out = {}
    by_model = {r['model']: r for r in q5_rows}
    for i, k in enumerate(keys):
        arms = {v: rows for v in REPEATS if (rows := load(v, k))}
        if not arms or not para.get(k) or k not in by_model:
            continue
        rep_ids = sorted(set.intersection(*(set(r) for r in arms.values())) & set(main[k]))
        para_ids = set(para[k]) & set(main[k]) & (keep if keep is not None else set(para[k]))
        common = sorted(set(rep_ids) & para_ids)
        reps = {v: paired(main[k], rows, rep_ids, lambda r: r['score'], 5900 + i) for v, rows in arms.items()}
        pr = by_model[k]
        out[k] = {'paraphrase': {x: pr[x] for x in ('items', 'diff', 'ci')}, 'repeats': reps,
                  'common_items': len(common),
                  'paraphrase_on_common': paired(main[k], para[k], common, lambda r: r['score'], 7900 + i) if common else None,
                  'repeats_on_common': {v: paired(main[k], rows, common, lambda r: r['score'], 6900 + i)
                                        for v, rows in arms.items()} if common else {},
                  'abs_within_repeat_spread': bool(abs(pr['diff']) <= max(abs(r['diff']) for r in reps.values()))}
    return out


def tau_noise(runs, keys, splits: int = 200, seed: int = 4242) -> dict:
    """Kendall's tau between the models' orderings on two halves of each template's items, over
    random halves: what tau looks like when only sampling noise separates the two arms."""
    rng = np.random.default_rng(seed)
    by_t = {k: collections.defaultdict(list) for k in keys}
    for k in keys:
        for r in runs[k].values():
            by_t[k][r['template_id']].append(r['score'])
    templates = sorted(by_t[keys[0]])
    taus = []
    for _ in range(splits):
        masks = {t: rng.permutation(len(by_t[keys[0]][t])) < len(by_t[keys[0]][t]) // 2 for t in templates}
        a = [np.mean([np.mean(np.array(by_t[k][t])[masks[t]]) for t in templates]) for k in keys]
        b = [np.mean([np.mean(np.array(by_t[k][t])[~masks[t]]) for t in templates]) for k in keys]
        taus.append(stats.kendalltau(a, b).statistic)
    return {'splits': splits, 'median': float(np.median(taus)), 'q1': float(np.percentile(taus, 25)),
            'q3': float(np.percentile(taus, 75))}


def sensitivity(runs, templates, keys, symbolic: set[str]):
    keep = set(templates) - set(SHORTCUT)
    out = []
    for k in keys:
        rows = list(runs[k].values())
        usable = [r for r in rows if not r['unusable']]
        def under(field):
            v = [score_with(r, field) for r in rows]
            return None if any(x is None for x in v) else float(np.mean(v))
        out.append({'model': k,
                    'half_tol': under('label_half_tol'),
                    'fitted': float(np.mean([r['score'] for r in rows])),
                    'double_tol': under('label_double_tol'),
                    'fully_solved': float(np.mean([fully(r) for r in rows])),
                    'unusable_excluded': float(np.mean([r['score'] for r in usable])),
                    'usable': len(usable),
                    'without_shortcut_templates': float(np.mean([r['score'] for r in rows if r['template_id'] in keep])),
                    'without_symbolic_templates': float(np.mean([r['score'] for r in rows if r['template_id'] not in symbolic])),
                    'half_unit': under('label_half_unit'),
                    'whole_trace': under('label_whole_trace')})
    head = [o['fitted'] for o in out]
    cols = ('half_tol', 'double_tol', 'fully_solved', 'unusable_excluded', 'without_shortcut_templates',
            'without_symbolic_templates', 'half_unit', 'whole_trace')
    tau = {c: (float(stats.kendalltau(head, [o[c] for o in out]).statistic) if all(o[c] is not None for o in out) else None)
           for c in cols}
    return {'models': out, 'tau_with_headline': tau, 'tau_noise': tau_noise(runs, keys)}


def class_label(r: dict) -> str:
    t = r['answer']['targets'] if r['answer'] else None
    if t is None:
        return '?'
    if r['template_id'] == FROUDE and t['numbers']:     # D-118: the answer is Fr, the regime follows
        fr = t['numbers'][0]
        return 'supercritical' if fr > 1 else 'subcritical' if fr < 1 else 'critical'
    if t['labeled']:
        return ', '.join(f'{a} {b}' for a, b in t['labeled'])
    return ', '.join(t['words']) or '(number)'


def provider_table(runs, keys) -> dict:
    """Per model with more than one serving endpoint: rows, raw mean and unusable share per provider,
    and the template-matched difference: over the templates that provider served, the mean of (its
    mean on the template minus the other providers' mean on the same template)."""
    out = {}
    for k in keys:
        rows = list(runs[k].values())
        provs = collections.Counter(r.get('provider') or '?' for r in rows)
        if len(provs) < 2:
            continue
        by_tp = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            by_tp[r['template_id']][r.get('provider') or '?'].append(r['score'])
        out[k] = {}
        for p, n in provs.most_common():
            diffs = []
            for t, d in by_tp.items():
                others = [s for q, v in d.items() if q != p for s in v]
                if p in d and others:
                    diffs.append(np.mean(d[p]) - np.mean(others))
            mine = [r for r in rows if (r.get('provider') or '?') == p]
            out[k][p] = {'rows': n, 'score': float(np.mean([r['score'] for r in mine])),
                         'unusable': float(np.mean([r['unusable'] for r in mine])),
                         'matched_diff': float(np.mean(diffs)) if diffs else None, 'templates_matched': len(diffs)}
    return out


def reported(runs, keys):
    """Also reported, not tested: branch, domain, answer type, tokens against score, per-label accuracy,
    the provider table."""
    gold_label = {}
    for rows in runs.values():                      # the label is the gold's; any model's row carries it
        for x, r in rows.items():
            if r['answer_type'] == 'classification' and r['answer'] and x not in gold_label:
                gold_label[x] = class_label(r)
    cls = collections.defaultdict(list)
    for x, lab in gold_label.items():
        cls[(runs[keys[0]][x]['template_id'], lab)].append(x)
    out = {'models': {}, 'classification': {}, 'providers': provider_table(runs, keys)}
    for k in keys:
        rows = list(runs[k].values())
        by = lambda f: {g: float(np.mean([r['score'] for r in rows if r[f] == g])) for g in sorted({r[f] for r in rows})}
        tok = lambda keep: float(np.median([r['completion_tokens'] or 0 for r in rows if keep(r)] or [np.nan]))
        out['models'][k] = {'branch': by('branch'), 'domain': by('domain'), 'level': by('level'),
                            'answer_type': by('answer_type'),
                            'median_tokens': tok(lambda r: True),
                            'median_tokens_fully_solved': tok(lambda r: fully(r) == 1),
                            'median_tokens_not_fully_solved': tok(lambda r: fully(r) == 0),
                            'score': float(np.mean([r['score'] for r in rows]))}
    for (t, lab), xs in sorted(cls.items()):
        out['classification'][f'{t}|{lab}'] = {'items': len(xs),
                                                **{k: float(np.mean([fully(runs[k][x]) for x in xs])) for k in keys}}
    return out


def repeats(main, keys):
    out = {}
    for k in keys:
        arms = {v: rows for v in REPEATS if (rows := load(v, k))}
        if not arms:
            continue
        ids = sorted(set.intersection(*(set(r) for r in arms.values())))
        scores = {v: float(np.mean([rows[x]['score'] for x in ids])) for v, rows in arms.items()}
        out[k] = {'items': len(ids), 'scores': scores,
                  'main_on_same_items': float(np.mean([main[k][x]['score'] for x in ids])),
                  'sd': float(np.std(list(scores.values()), ddof=1)) if len(scores) > 1 else None,
                  'range': float(max(scores.values()) - min(scores.values())),
                  'same_verdict_every_repeat': float(np.mean([len({verdict3(rows[x]) for rows in arms.values()}) == 1
                                                              for x in ids]))}
    return out


def set_aside(store='main'):
    out = {}
    for k in SET_ASIDE:
        rows = load(store, k)
        if rows:
            s = np.array([np.mean(v) for v in per_template(rows.values(), lambda r: r['score']).values()])
            out[k] = {'score': float(s.mean()), 'ci': boot_mean(s, 1999),
                      'fully_solved': float(np.mean([fully(r) for r in rows.values()])),
                      'unusable': sum(r['unusable'] for r in rows.values())}
    return out


def milestone_facts(rows) -> dict:
    n = per_template(rows, lambda r: r['milestones_required'])
    return {'items_without': sum(v.count(0) for v in n.values()),
            'templates_all_without': sum(all(c == 0 for c in v) for v in n.values()),
            'templates_some_without': sum(0 in v and any(v) for v in n.values()),
            'templates_with_item_at_most_one': sum(min(v) <= 1 for v in n.values()),
            'templates_with_item_exactly_one': sum(1 in v for v in n.values())}


def per_template_rows(runs, keys, levels, single) -> list[dict]:
    """One row per template and model, aggregates only (D-149)."""
    out = []
    for k in keys:
        by_t = per_template(runs[k].values(), lambda r: r)
        for t, rs in sorted(by_t.items()):
            cov = [r['e3']['coverage'] for r in rs if r['e3']['coverage'] is not None]
            out.append({'template_id': t, 'branch': rs[0]['branch'], 'domain': rs[0]['domain'], 'level': levels[t],
                        'answer_type': rs[0]['answer_type'], 'single_path': int(single[t]), 'model': k, 'items': len(rs),
                        'answer_score': round(float(np.mean([r['score'] for r in rs])), 6),
                        'fully_solved': round(float(np.mean([fully(r) for r in rs])), 6),
                        'correct': sum(verdict3(r) == 'correct' for r in rs), 'partial': sum(verdict3(r) == 'partial' for r in rs),
                        'incorrect': sum(verdict3(r) == 'incorrect' for r in rs), 'unusable': sum(r['unusable'] for r in rs),
                        'e3_coverage': round(float(np.mean(cov)), 6) if cov else '',
                        'digit_flag_rate': round(float(np.mean([flagged(r) for r in rs if r['status'] == 'answered'] or [0])), 6),
                        'median_completion_tokens': int(np.median([r['completion_tokens'] or 0 for r in rs]))})
    return out


# ------------------------------------------------------------------ provenance

def git(*args) -> str:
    return subprocess.run(['git', *args], cwd=REPO, capture_output=True, text=True).stdout.strip()


def provenance(store: str) -> dict:
    """This script's commit and LF hash, the store's CONFIG, and each stage's CONFIG, each checked
    against the digest of the store CONFIG it was built from (D-144)."""
    me = Path(__file__)
    cfg_path = SCORES / store / 'CONFIG.json'
    cfg = read_json(cfg_path)
    out = {'analyze': {'git': git('rev-parse', 'HEAD'), 'tag': git('describe', '--tags', '--exact-match') or None,
                       'dirty': bool(git('status', '--porcelain', '--', me.relative_to(REPO).as_posix())),
                       'sha256_lf': hashlib.sha256(me.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                       'B': B, 'B_TEST': B_TEST},
           'store': {'name': store, **(cfg or {})}, 'stages': {}}
    digest = hashlib.sha256(cfg_path.read_bytes()).hexdigest() if cfg_path.exists() else None
    for stage in ('e5', 'router'):
        sc = read_json(SCORES / store / stage / 'CONFIG.json')
        if sc is None:
            continue
        if sc.get('store_config_sha256') and sc['store_config_sha256'] != digest:
            raise SystemExit(f'scores/{store}/{stage} was built from another score store than the one on disk: '
                             f'run `{"judge" if stage == "e5" else "router"} --score` to rebuild its rows first')
        out['stages'][stage] = sc
    return out


# ------------------------------------------------------------------ report

def ci(c) -> str:
    return f'{c[0]:.3f} to {c[1]:.3f}'


def f3(v) -> str:
    return '' if v is None else f'{v:.3f}'


def render(res) -> str:
    m1, pairs = res['q1']['models'], res['q1']['pairs']
    sig = sum(p['p_holm'] < 0.05 for p in pairs)
    same = sum(p['same_verdict'] for p in pairs)
    same_t = sum(p['same_verdict_template'] for p in pairs)
    differ = [p for p in pairs if not p['same_verdict']]
    only_mcnemar = sum(p['mcnemar_p_holm'] < 0.05 <= p['p_holm'] for p in differ)
    mf = res['milestones']
    q3r = res['q3']
    has_e5 = any('e5_coverage_on_wrong' in r for r in q3r)
    has_router = any('router_rate_on_fully_solved' in r for r in q3r)
    incomplete = [r['model'] for r in q3r if r.get('e5_incomplete') or r.get('router_incomplete')]
    stack = ("The answer check, E3 and E4's digit rule" + (', E5' if has_e5 else '') + (' and the step router' if has_router else '')
             + ('' if has_e5 else "; E5's judge has not run, so its columns are empty"))
    L = ['# Results of the full run' + (': the deterministic stack' if not has_e5 else ''), '',
         'Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the '
         'method choices the plan leaves open are fixed in the script\'s docstring, and the corrections and '
         f'additions made after the first results were read are labelled where they appear (D-146 to D-149). {stack}. '
         'Every interval is 95% and resamples templates (B = 10,000): all 150 for a model\'s score, and the templates '
         'with a qualifying trace for a rate over a subset of traces. Every resampling test draws 100,000 permutations, '
         'so its Holm-adjusted floor over 55 pairs is 0.0006.'
         + (f' **Incomplete judged stages:** {", ".join(f"`{m}`" for m in incomplete)} have calls without a reply; '
            'their judged rates are over the answered calls only.' if incomplete else ''), '',
         '## Q1. Answer score per model', '',
         'Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0. SD within: the mean over '
         'templates of the SD of the score across a template\'s 15 items; SD between: the SD of the 150 template means '
         '(the per-template variability the July rebuttal promised). Unusable is empty plus unreadable; capped: answered '
         'rows that stopped at the output cap and are scored on what they state (D-117).', '',
         '| model | answer score | 95% CI | SD within | SD between | fully solved | 95% CI | correct | partial | incorrect | '
         'unusable (empty + unreadable) | capped, scored |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in sorted(m1, key=lambda r: -r['score']):
        L.append(f"| `{r['model']}` | {r['score']:.3f} | {ci(r['ci'])} | {r['within_template_sd']:.3f} | "
                 f"{r['between_template_sd']:.3f} | {r['fully_solved']:.3f} | {ci(r['fully_ci'])} | "
                 f"{r['correct']} | {r['partial']} | {r['incorrect']} | {r['unusable']} ({r['empty']} + {r['unreadable']}) | "
                 f"{r['capped_scored']} |")
    odd = sum(r['odd_finish_scored'] for r in m1)
    L += ['', f'Of the 55 pairs, {sig} differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation '
          f'over the 150 per-template differences). The same template-level test on the fully-solved rate gives the '
          f'same verdict on {same_t} of 55. McNemar\'s exact test on the paired item verdicts (Holm) gives the same '
          f'verdict on {same} of 55; on {only_mcnemar} of the {len(differ)} others McNemar holds and the template-level '
          'test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), '
          f'so a claim rests on the template-level tests. {odd} answered rows across the roster ended with a finish '
          'reason other than stop or length (a provider fault inside a 200) and are scored on what they state; the '
          'harness now retries such a reply (D-148).', '',
          '| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | template p (Holm) | McNemar p (Holm) | same verdict |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for p in sorted(pairs, key=lambda p: (p['p_holm'], -abs(p['diff']))):
        L.append(f"| `{p['a']}` | `{p['b']}` | {p['diff']:+.3f} | {ci(p['ci'])} | {p['p_holm']:.4f} | "
                 f"{p['fully_diff']:+.3f} | {ci(p['fully_ci'])} | {p['fully_p_holm']:.4f} | {p['mcnemar_p_holm']:.4f} | "
                 f"{'yes' if p['same_verdict_template'] else 'no'} |")
    n_welch = sum(r['p_welch_holm'] < 0.05 for r in res['q2'])
    n_perm = sum(r['p_perm_holm'] < 0.05 for r in res['q2'])
    L += ['', '## Q2. The complexity cliff: Easy minus Advanced', '',
          'Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier. '
          '**Corrected after the first results (D-146):** the planned test permuted the tier labels on the raw '
          'difference, which is liberal when the smaller tier has the larger spread, and Advanced template means '
          'spread two to four times as widely as Easy ones (the two SD columns); its count of models also moved with '
          f'the seed. The count now rests on Welch\'s t-test, Holm across the eleven models: **{n_welch} of 11** hold, '
          f'against {n_perm} under the planned permutation, printed beside it. The detectable gap at 80% power is given '
          'as the plan defined it (2.80 x pooled sigma x sqrt(1/58 + 1/34)) and at the strictest Holm step from the '
          'Welch standard error.', '',
          '| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | planned p (Holm) | SD Easy | SD Adv | '
          'detectable, planned | detectable, Holm |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in res['q2']:
        L.append(f"| `{r['model']}` | {r['easy']:.3f} | {r['advanced']:.3f} | {r['gap']:+.3f} | {ci(r['ci'])} | "
                 f"{r['p_welch_holm']:.4f} | {r['p_perm_holm']:.4f} | {r['sd_easy']:.3f} | {r['sd_advanced']:.3f} | "
                 f"{r['detectable_planned']:.3f} | {r['detectable_holm']:.3f} |")
    L += ['', '**The cliff under two variations.** Unusable rows left out of the template means, because an empty '
          'row at the output cap measures finishing within the ceiling as well as solving, and most such rows fall on '
          'Advanced templates; and without the nine symbolic templates (D-138). Welch p, Holm across models.', '',
          '| model | gap, as scored | gap, unusable left out | 95% CI | Welch p (Holm) | gap, no symbolic | 95% CI | Welch p (Holm) |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in res['q2']:
        u, s = r['unusable_excluded'], r['without_symbolic']
        L.append(f"| `{r['model']}` | {r['gap']:+.3f} | {u['gap']:+.3f} | {ci(u['ci'])} | {u['p_welch_holm']:.4f} | "
                 f"{s['gap']:+.3f} | {ci(s['ci'])} | {s['p_welch_holm']:.4f} |")
    L += ['', '## Q3. What the process scores add beyond the answer', '',
          'Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in '
          f"neither. E3 coverage leaves out the {mf['items_without']} items with no milestones (all 15 of "
          f"{mf['templates_all_without']} templates and some items of {mf['templates_some_without']} more); "
          f"{mf['templates_with_item_at_most_one']} templates have an item with at most one milestone "
          f"({mf['templates_with_item_exactly_one']} with exactly one), where coverage is close to an answer check. "
          'E3 is the deterministic part of E5, which adds the judge\'s verdict on the milestones E3 does not find. '
          'The floor is the same trace scored against a sibling item\'s milestones (D-149): what coverage a trace '
          'reaches by chance, per model, on the readable wrong answers. '
          "The digit rule's flag counts a trace when any step is flagged. Against the experts' step labels on the "
          "pilot's fully solved traces it has precision 0.817 and recall 0.427 (0.750 and 0.320 before D-156; "
          "0.800 and 0.427 between D-156 and D-159, unchanged by D-160; SCORER_VALIDATION.md). *Added 2026-09-30 "
          'and 2026-10-01 (D-154 to D-163), after the first results were read:* on this roster a domain expert found 111 of 220 '
          'sampled flags real before the first fix, 0.505 (FLAG_REVIEW.md), and 155 of 206 real after it, 0.752 '
          "(FLAG_REVIEW_2.md). The second fix, built from those notes, took the flag off 42 of their 51 misread "
          "steps and left it on 154 of the 155 slips' steps (DIGIT_FIX_2.md). Two review agents then checked both "
          'fixes; the corrections they led to, with LaTeX control spaces now read, took the flag off 82 full-run '
          "steps and put it on 370, and kept every slip's step the second fix had kept (D-160, DIGIT_FIX_3.md). "
          "The same expert then read 190 flags of the rule as it now stands: 171 of the 189 decided are real, "
          "0.905 (95% Wilson 0.854 to 0.939; per model 0.632 to 1.000; FLAG_REVIEW_3.md), and by the owner's "
          "rule no fourth reading follows. So its rate is still not a count of slips; and it "
          "reads a different amount of arithmetic in each model's traces (claims checked per answered trace, "
          'shown), so a low rate can mean little was read.', '',
          '**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it '
          'stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.', '',
          '| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | floor, readable | '
          'E5 coverage |', '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in q3r:
        e5_cell = (f"{r['e5_coverage_on_wrong']:.3f} ({ci(r['e5_ci'])}); readable {r['e5_coverage_on_readable_wrong']:.3f}"
                   + (f" (incomplete: {r['e5_without_reply']} without a reply)" if r.get('e5_incomplete') else '')
                   if 'e5_coverage_on_wrong' in r else 'pending')
        L.append(f"| `{r['model']}` | {r['wrong']} | {r['wrong_unusable']} | {r['e3_coverage_on_wrong']:.3f} | "
                 f"{ci(r['e3_ci'])} | {r['e3_coverage_on_readable_wrong']:.3f} | {ci(r['e3_readable_ci'])} | "
                 f"{f3(r['e3_null_on_readable_wrong'])} | {e5_cell} |")
    if has_e5:
        L += ['', "**E5's judge**, MiMo-V2.5-Pro on the milestones E3 did not find (`judge.py`). E5 coverage is "
              "E5-strict: E3's milestones plus those the judge rules REACHED. On the pilot the judge never called a "
              'fabricated value REACHED and found 76% of true ones, so the score is conservative (RESULTS_E5). '
              'Judged fraction: milestones sent to the judge over those required, on answered traces; unjudged: '
              'milestones the reply did not name, over those sent. A trace whose call got no reply is left out of '
              'every rate and counted (D-148).', '',
              '| model | calls | without a reply | judged fraction | of the judged, REACHED | unjudged |',
              '|---|---:|---:|---:|---:|---:|']
        for r in q3r:
            if 'e5_coverage_on_wrong' in r:
                L.append(f"| `{r['model']}` | {r['e5_calls']} | {r['e5_without_reply']} | {r['e5_judged_fraction']:.3f} | "
                         f"{r['e5_reached_share']:.3f} | {r['e5_unjudged_rate']:.3f} |")
    L += ['', '**The digit rule.** The wrong-answer rate is over the answered wrong answers (D-149: an empty trace '
          'has no step to flag). Beside it, the 1% tolerance E4 shipped with, blind to most slips (RESULTS_X1), and '
          'where the first flag falls in a fully solved trace (0 = first step, 1 = last).', '',
          '| model | wrong answers, answered | flag rate | 95% CI | fully solved traces | flag rate | 95% CI | '
          'at 1% | claims checked per trace | traces with a claim | first flag: traces, median position |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in q3r:
        ff = r['first_flag_position']
        L.append(f"| `{r['model']}` | {r['wrong_answered']} | {r['digit_flag_rate_on_wrong']:.3f} | {ci(r['digit_wrong_ci'])} | "
                 f"{r['fully_solved']} | {r['digit_flag_rate_on_fully_solved']:.3f} | {ci(r['digit_ci'])} | "
                 f"{r['tol1_flag_rate_on_fully_solved']:.3f} | {r['claims_per_trace']:.2f} | {r['traces_with_a_claim']:.3f} | "
                 f"{ff['traces']}, {f3(ff['median'])} |")
    L += ['', '**Wrong-answer rate against the item\'s milestone count** (D-149): how failure grows with the depth '
          'of the gold derivation.', '',
          '| model | ' + ' | '.join(f'{b[0]} ({q3r[0]["by_milestone_count"][b[0]]["items"]} items)' for b in MS_BUCKETS) + ' |',
          '|---|' + '---:|' * len(MS_BUCKETS)]
    for r in q3r:
        L.append(f"| `{r['model']}` | " + ' | '.join(f3(r['by_milestone_count'][b[0]]['wrong_rate']) for b in MS_BUCKETS) + ' |')
    if has_router:
        L += ['', "**The step router** (`router.py`): the digit rule's flags, and MiMo-V2.5-Pro on every other step in "
              'one batched call per trace. On the pilot\'s labelled traces it had precision 0.707 and recall 0.603 over '
              'all steps, and 0.703 and 0.360 inside correct-answer traces (ROUTER_VALIDATION.md), so its rates are '
              'flags, not counts of errors. A trace counts as flagged when any step is; a trace whose call got no reply '
              'is left out and counted.', '',
              '| model | calls | without a reply | flagged, fully solved | 95% CI | by the judge | 95% CI | '
              'flagged, wrong answers | 95% CI | steps flagged per trace | unjudged steps |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in q3r:
            if 'router_rate_on_fully_solved' in r:
                L.append(f"| `{r['model']}` | {r['router_calls']} | {r['router_without_reply']} | "
                         f"{r['router_rate_on_fully_solved']:.3f} | {ci(r['router_ci'])} | "
                         f"{r['router_judge_rate_on_fully_solved']:.3f} | {ci(r['router_judge_ci'])} | "
                         f"{r['router_rate_on_wrong']:.3f} | {ci(r['router_wrong_ci'])} | "
                         f"{r['router_steps_flagged_per_trace']:.2f} | {r['router_unjudged_rate']:.3f} |")
    if any('attribution_on_wrong' in r for r in q3r):
        L += ['', '**Also reported, not tested: what points at a wrong answer.** On the answered traces that score 0, '
              'the share with a digit-rule flag, with a milestone E5 rules MISSING, and with a step the router\'s '
              'judge flags. A trace can be in several columns, or in none.', '',
              '| model | answered wrong-answer traces | digit rule | E5 MISSING | router judge |', '|---|---:|---:|---:|---:|']
        for r in q3r:
            a = r.get('attribution_on_wrong')
            if a:
                L.append(f"| `{r['model']}` | {a['traces']} | {f3(a['digit_rule'])} | {f3(a['e5_missing'])} | "
                         f"{f3(a['router_judge'])} |")
    s0 = res['q4'][0]
    L += ['', '## Q4. Consistency within a template', '',
          'The share of templates fully solved on all 15 instances, on some, and on none: the '
          f"{s0['single_path']['templates']} single-path templates (one reasoning path across their instances, "
          f"`diversity.py`'s lower reading) and the {s0['multi_path']['templates']} others. Intervals are in "
          '`results.json`.', '',
          '| model | single: all | some | none | others: all | some | none |', '|---|---:|---:|---:|---:|---:|---:|']
    for r in res['q4']:
        s, m = r['single_path'], r['multi_path']
        L.append(f"| `{r['model']}` | {s['all']:.3f} | {s['some']:.3f} | {s['none']:.3f} | {m['all']:.3f} | "
                 f"{m['some']:.3f} | {m['none']:.3f} |")
    L += ['', '## Q5. Paraphrase robustness', '']
    if res['q5']['models']:
        es = res['q5'].get('expert_stats')
        L += ['**Provisional: the experts\' check of the paraphrases has not returned**, so no pair has been '
              'dropped yet.' if es is None else
              f"The experts' check: {es['returned']} pairs returned, {es['rejected']} rejected and dropped from both "
              f"arms, {es['outstanding']} not yet returned and kept provisionally"
              + (' (**provisional until they return**).' if es['outstanding'] else '.'), '',
              'Paraphrase minus original, paired by item; besides the answer score, E3 coverage on the items with '
              'milestones, E5-strict where both arms carry it, and the answer score on the pairs both arms served '
              'from the same endpoint (D-149).', '',
              '| model | items | answer score diff | 95% CI | p (Holm) | McNemar p (Holm) | E3 coverage diff | 95% CI | '
              'E5 diff | 95% CI | same provider: pairs, diff |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in res['q5']['models']:
            e3c, e5c, sp = r.get('e3'), r.get('e5'), r.get('same_provider')
            e3_cell = (f"{e3c['diff']:+.3f} | {ci(e3c['ci'])}") if e3c else ' | '
            e5_cell = (f"{e5c['diff']:+.3f} | {ci(e5c['ci'])}") if e5c else 'pending | '
            sp_cell = f"{sp['items']}, {sp['diff']:+.3f}" if sp else ''
            L.append(f"| `{r['model']}` | {r['items']} | {r['diff']:+.3f} | {ci(r['ci'])} | {r['p_holm']:.4f} | "
                     f"{r['mcnemar_p_holm']:.4f} | {e3_cell} | {e5_cell} | {sp_cell} |")
        if res['q5']['tau']:
            t = res['q5']['tau']
            L += ['', f"Kendall's tau between the models' answer scores on the originals and on the paraphrases, over "
                  f"the {t['items']} items every tested model holds: {t['tau']:.3f}, 95% CI {ci(t['ci'])}; the noise "
                  f"floor for tau on this roster is {res['sensitivity']['tau_noise']['median']:.3f} (below)."]
            na = t.get('noise_arm')
            if na:
                L += ['', f"The noise floor at the arm's own size and template mix (D-163): two disjoint draws of "
                      f"{na['items_per_half']} main-run items with the arm's per-template counts over its {na['templates']} "
                      f"templates, the models ordered on each, {na['splits']} draws: median tau {na['median']:.3f}, quartiles "
                      f"{na['q1']:.3f} to {na['q3']:.3f}, 5th to 95th percentile {na['p5']:.3f} to {na['p95']:.3f}. This, not the "
                      f"roster-wide floor (halves of 7 or 8 items per template), is the comparison for the arm's tau."]
        margin = res['q5'].get('margin')
        if margin is not None:
            L += ['', f"**Bounds (D-163).** The same bootstrap's 90% interval per model, read against a margin of "
                  f"±{margin:.2f} in answer score: a model is within the margin when the whole interval is (the "
                  f"two-one-sided-tests rule at 5%). The margin is the plan's largest detectable paired difference "
                  f"(5.1 points at 15% discordance); it was fixed after the point estimates were known and before "
                  f"these intervals were computed.", '',
                  '| model | answer score diff | 90% CI | within the margin |', '|---|---:|---:|---|']
            for r in res['q5']['models']:
                L.append(f"| `{r['model']}` | {r['diff']:+.3f} | {ci(r['ci90'])} | {'yes' if r['within_margin'] else 'no'} |")
        vr = res['q5'].get('vs_repeats')
        if vr:
            L += ['', "**Against run-to-run noise (D-163).** For the models with decoding repeats: the paraphrase "
                  "difference beside each repeat minus the main run on the 300 repeat items, each paired by item with "
                  "the same template bootstrap; and both on the kept pairs the two subsamples share.", '',
                  '| model | paraphrase − original (pairs) | ' + ' | '.join(f'{v} − main (300)' for v in REPEATS)
                  + ' | shared items | paraphrase on them | repeats on them | |paraphrase| within the repeats\' spread |',
                  '|---|---:|' + '---:|' * len(REPEATS) + '---:|---:|---:|---|']
            for k, v in vr.items():
                p = v['paraphrase']
                reps = ' | '.join((f"{v['repeats'][x]['diff']:+.3f} {ci(v['repeats'][x]['ci'])}" if x in v['repeats'] else '')
                                  for x in REPEATS)
                pc = v['paraphrase_on_common']
                pc_cell = f"{pc['diff']:+.3f} {ci(pc['ci'])}" if pc else ''
                rc_cell = ', '.join(f"{v['repeats_on_common'][x]['diff']:+.3f}" for x in REPEATS if x in v['repeats_on_common'])
                L.append(f"| `{k}` | {p['diff']:+.3f} {ci(p['ci'])} ({p['items']}) | {reps} | {v['common_items']} | {pc_cell} | "
                         f"{rc_cell} | {'yes' if v['abs_within_repeat_spread'] else 'no'} |")
    else:
        L.append('Not run yet.')
    sens = res['sensitivity']
    tn = sens['tau_noise']
    L += ['', '## Sensitivity', '',
          "Answer score under each variation; the last row is Kendall's tau between that ordering of the models "
          'and the headline one. The three added readings (D-147, D-149) the experts could not arbitrate: the '
          'half-unit window requires a correct rounding at the precision shown where the rule accepts one unit '
          'either way; the whole-trace reading credits a quantity the question asks for when it is stated in the '
          'body and left off the Answer line (the prompt asks for it there); the pool without the nine symbolic '
          'templates, whose answers the check scores by the numbers they state (D-138). "Unusable excluded" is an item '
          'mean over the usable rows, not a mean of template means. '
          f'Tau\'s noise floor: the ordering on one random half of each template\'s items against the other, median '
          f"{tn['median']:.3f} over {tn['splits']} splits (quartiles {tn['q1']:.3f} to {tn['q3']:.3f}); the top five models "
          'lie within 0.012 of each other, so tau falls below 1 from sampling noise alone.', '',
          '| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | '
          'without the 4 shortcut templates | without the 9 symbolic templates | half-unit window | whole trace |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in sens['models']:
        L.append(f"| `{r['model']}` | {f3(r['half_tol'])} | {r['fitted']:.3f} | {f3(r['double_tol'])} | "
                 f"{r['fully_solved']:.3f} | {r['unusable_excluded']:.3f} | {r['without_shortcut_templates']:.3f} | "
                 f"{r['without_symbolic_templates']:.3f} | {f3(r['half_unit'])} | {f3(r['whole_trace'])} |")
    t = sens['tau_with_headline']
    L += [f"| tau with the headline | {f3(t['half_tol'])} | 1 | {f3(t['double_tol'])} | {f3(t['fully_solved'])} | "
          f"{f3(t['unusable_excluded'])} | {f3(t['without_shortcut_templates'])} | {f3(t['without_symbolic_templates'])} | "
          f"{f3(t['half_unit'])} | {f3(t['whole_trace'])} |", '',
          "The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had "
          'not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).']
    rep = res['reported']
    first = next(iter(rep['models'].values()))
    L += ['', '## Also reported, not tested', '', '### By branch and level', '',
          '| model | ' + ' | '.join(b.replace('_engineering', '') for b in first['branch']) +
          ' | Easy | Intermediate | Advanced |', '|---|' + '---:|' * (len(first['branch']) + 3)]
    for k, r in rep['models'].items():
        L.append(f'| `{k}` | ' + ' | '.join(f'{v:.3f}' for v in r['branch'].values()) +
                 f" | {r['level']['Easy']:.3f} | {r['level']['Intermediate']:.3f} | {r['level']['Advanced']:.3f} |")
    keys = list(rep['models'])
    for title, field in (('By domain', 'domain'), ('By answer type', 'answer_type')):
        L += ['', f'### {title}', '', '| ' + field.replace('_', ' ') + ' | ' + ' | '.join(f'`{k}`' for k in keys) + ' |',
              '|---|' + '---:|' * len(keys)]
        for g in first[field]:
            L.append(f'| {g} | ' + ' | '.join(f"{rep['models'][k][field][g]:.3f}" for k in keys) + ' |')
    L += ['', '### Tokens against score', '', 'Median completion tokens as billed.', '',
          '| model | answer score | median tokens | on fully solved | on the rest |', '|---|---:|---:|---:|---:|']
    for k, r in rep['models'].items():
        L.append(f"| `{k}` | {r['score']:.3f} | {r['median_tokens']:.0f} | {r['median_tokens_fully_solved']:.0f} | "
                 f"{r['median_tokens_not_fully_solved']:.0f} |")
    L += ['', '### By serving endpoint', '',
          'Open-weight models were served by several endpoints under the fp8-or-better rule (D-133). The harness '
          'dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is '
          'confounded with the templates each endpoint happened to serve; the matched difference compares an endpoint '
          'with the other endpoints on the same templates (D-149).', '',
          '| model | endpoint | rows | raw score | unusable | matched difference | templates matched |',
          '|---|---|---:|---:|---:|---:|---:|']
    for k, provs in rep['providers'].items():
        for p, v in provs.items():
            md = '' if v['matched_diff'] is None else f"{v['matched_diff']:+.3f}"
            L.append(f"| `{k}` | {p} | {v['rows']} | {v['score']:.3f} | {v['unusable']:.3f} | {md} | {v['templates_matched']} |")
    L += ['', '### Classification templates: fully solved by the gold label', '',
          'The pool over-represents the rare labels by design (D-116), so these rates are per label, not a pooled '
          'accuracy.', '', '| template | gold label | items | ' + ' | '.join(f'`{k}`' for k in keys) + ' |',
          '|---|---|---:|' + '---:|' * len(keys)]
    for key, r in rep['classification'].items():
        t, lab = key.split('|')
        L.append(f"| `{t.removeprefix('template_')}` | {lab} | {r['items']} | " +
                 ' | '.join(f'{r[k]:.3f}' for k in keys) + ' |')
    if res['repeats']:
        L += ['', '## Decoding repeats', '', '| model | items | ' + ' | '.join(REPEATS) +
              ' | main run, same items | SD | range | same verdict in every repeat |',
              '|---|---:|' + '---:|' * (len(REPEATS) + 4)]
        for k, r in res['repeats'].items():
            sd = '' if r['sd'] is None else f"{r['sd']:.3f}"
            L.append(f"| `{k}` | {r['items']} | " +
                     ' | '.join(f"{r['scores'][v]:.3f}" if v in r['scores'] else '' for v in REPEATS) +
                     f" | {r['main_on_same_items']:.3f} | {sd} | {r['range']:.3f} | {r['same_verdict_every_repeat']:.3f} |")
    if res['set_aside']:
        L += ['', '## Set aside (D-132), outside every comparison', '']
        L += [f"- `{k}`: answer score {v['score']:.3f} (95% CI {ci(v['ci'])}), fully solved {v['fully_solved']:.3f}, "
              f"unusable {v['unusable']}" for k, v in res['set_aside'].items()]
    pv, st = res['provenance'], res['provenance']['store']
    stages = ', '.join(f"{s} at `{str(c.get('git', '?'))[:7]}`" for s, c in pv['stages'].items())
    L += ['', '## Provenance', '',
          f"`analyze.py` at commit `{pv['analyze']['git'][:7]}`" + (f" (tag `{pv['analyze']['tag']}`)" if pv['analyze']['tag'] else '')
          + (', dirty' if pv['analyze']['dirty'] else '') + f"; the score store `{st['name']}` scored at commit "
          f"`{str(st.get('git', '?'))[:7]}`" + (f" (tag `{st['tag']}`)" if st.get('tag') else '')
          + (', dirty' if st.get('dirty') else '') + f" on {st.get('written_at_utc', st.get('scored_at_utc', '?'))}"
          + (f'; stages: {stages}' if stages else '')
          + '. The evaluator hashes and the per-model trace hashes are in `results.json`.']
    return '\n'.join(L) + '\n'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--store', default='main',
                    help='the scores/ directory to read; another one is written as a labelled record')
    ap.add_argument('--note', default='', help='a line saying what that other store is')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    runs = {k: load(a.store, k) for k in ROSTER}
    missing = [k for k, v in runs.items() if v is None]
    if missing:
        raise SystemExit(f'no {a.store} scores for {missing}: run score.py first')
    templates, levels, single = check_store(runs)
    symbolic = {r['template_id'] for r in runs[ROSTER[0]].values() if r['answer_type'] == 'symbolic'}
    pv = provenance(a.store)
    para = {k: load('paraphrase', k) for k in ROSTER} if a.store == 'main' else {}
    e5_main = {k: load_stage('e5', k, 'main') for k in ROSTER} if a.store == 'main' else None
    e5_para = {k: load_stage('e5', k, 'paraphrase') for k in ROSTER} if a.store == 'main' else None
    check = accepted_pairs()
    res = {'q1': q1(runs, templates, ROSTER), 'q2': q2(runs, templates, levels, ROSTER, symbolic),
           'q3': q3(runs, ROSTER, a.store), 'q4': q4(runs, templates, single, ROSTER),
           'q5': {**q5(runs, para, ROSTER, check['keep'] if check else None, e5_main, e5_para),
                  'expert_stats': {k: v for k, v in check.items() if k != 'keep'} if check else None},
           'sensitivity': sensitivity(runs, templates, ROSTER, symbolic), 'reported': reported(runs, ROSTER),
           'repeats': repeats(runs, ROSTER) if a.store == 'main' else {},
           'set_aside': set_aside(a.store), 'milestones': milestone_facts(runs[ROSTER[0]].values()),
           'symbolic_templates': sorted(symbolic), 'provenance': pv}
    out = OUT if a.store == 'main' else OUT / a.store
    out.mkdir(parents=True, exist_ok=True)
    (out / 'results.json').write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8', newline='\n')
    rows = per_template_rows(runs, ROSTER, levels, single)
    with open(out / 'per_template.csv', 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    text = render(res)
    if a.store != 'main':
        text = text.replace('\n', f'\n\n**A record, not the results:** computed from `scores/{a.store}`. {a.note}\n', 1)
    (out / 'RESULTS.md').write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0


# ------------------------------------------------------------------ self-test

def fake(tid, n, score, level='Easy', single=True, provider='p'):
    lab = {1.0: 'correct', 0.5: 'partial', 0.0: 'incorrect'}[score]
    return {'item_id': f'{tid}-{n}', 'template_id': tid, 'score': score, 'unusable': False, 'status': 'answered',
            'finish_reason': 'stop', 'provider': provider, 'milestones_required': 2, 'completion_tokens': 10,
            'answer': {'label': lab, 'label_half_tol': lab, 'label_double_tol': lab, 'label_half_unit': lab,
                       'label_whole_trace': lab, 'targets': None},
            'e3': {'coverage': 0.5, 'null_coverage': 0.1}, 'steps': [], 'level': level, 'single_path': single,
            'branch': 'b', 'domain': 'd', 'answer_type': 'scalar'}


def selftest() -> int:
    """Known answers on synthetic data. Writes nothing."""
    bad = []
    rng = np.random.default_rng(0)
    h = holm([0.01, 0.04, 0.03])
    if not (abs(h[0] - 0.03) < 1e-12 and abs(h[1] - 0.06) < 1e-12 and abs(h[2] - 0.06) < 1e-12):
        bad.append(f'holm {h}')
    null = rng.normal(0, 1, 150)
    null -= null.mean()
    if sign_flip_p(null, 1) < 0.5:
        bad.append('sign flip rejected a zero-mean difference')
    if sign_flip_p(null + 1.0, 2) > 0.001:
        bad.append('sign flip missed a shift of 1 SD')
    if mcnemar_p(np.array([1] * 30 + [0] * 70), np.array([1] * 30 + [0] * 70)) != 1.0:
        bad.append('mcnemar with no discordant pair')
    if abs(mcnemar_p(np.ones(10), np.zeros(10)) - 2 * 0.5 ** 10) > 1e-12:
        bad.append('mcnemar exact tail')
    if boot_mean(np.ones(150), 3) != [1.0, 1.0]:
        bad.append('bootstrap of a constant')
    p, _c, n = cluster_mean({'a': [1, 1, 0], 'b': [0]}, 5)
    if (p, n) != (0.5, 4):
        bad.append(f'cluster mean {p} {n}')
    x = np.r_[np.full(58, 0.9), np.full(34, 0.5)] + rng.normal(0, 0.01, 92)
    easy = np.r_[np.ones(58, bool), np.zeros(34, bool)]
    if tier_perm_p(x, easy, 4) > 0.001:
        bad.append('tier permutation missed a 0.4 gap')
    if tier_perm_p(rng.permutation(np.r_[np.full(46, 0.9), np.full(46, 0.5)]), easy, 6) < 0.01:
        bad.append('tier permutation rejected shuffled tiers')
    # Welch: a real gap is found; equal means with unequal spreads are not called different at 0.05 more
    # often than chance (checked on the average over 200 draws), where the raw-difference permutation is liberal
    w = welch(np.full(58, 0.9) + rng.normal(0, 0.05, 58), np.full(34, 0.5) + rng.normal(0, 0.2, 34))
    if w['p_welch'] > 1e-6 or not (w['detectable_holm'] > w['detectable_welch'] > 0):
        bad.append(f'welch {w}')
    rej_w = rej_p = 0
    for d in range(200):
        xe, xa = rng.normal(0.8, 0.05, 58), rng.normal(0.8, 0.30, 34)
        rej_w += welch(xe, xa)['p_welch'] < 0.05
        rej_p += tier_perm_p(np.r_[xe, xa], easy, 10 + d, draws=2000) < 0.05
    if rej_w > 25 or rej_p < rej_w:
        bad.append(f'welch rejects {rej_w} of 200 nulls with unequal spread; the permutation {rej_p}')
    # Q1: a model against itself differs by nothing
    rows = {f't{t:03d}-{n}': fake(f't{t:03d}', n, float(rng.integers(0, 2))) for t in range(150) for n in range(15)}
    r1 = q1({'a': rows, 'b': dict(rows)}, sorted({r['template_id'] for r in rows.values()}), ['a', 'b'])
    p1 = r1['pairs'][0]
    if p1['diff'] != 0 or p1['p'] != 1.0 or p1['fully_p'] != 1.0 or p1['mcnemar_p'] != 1.0:
        bad.append(f'q1 self-comparison {p1}')
    if abs(r1['models'][0]['between_template_sd'] - np.std([np.mean([rows[f"t{t:03d}-{n}"]["score"] for n in range(15)])
                                                            for t in range(150)], ddof=1)) > 1e-12:
        bad.append('between-template SD')
    # Q5: eleven models, three items per template; models 0-2 lose some solved items on the paraphrase
    main_, para_ = {}, {}
    for m in range(11):
        k = f'm{m}'
        main_[k], para_[k] = {}, {}
        for t in range(150):
            for n in range(3):
                base = float(rng.random() < 0.3 + 0.05 * m)
                drop = m < 3 and base == 1.0 and rng.random() < 0.3
                main_[k][f't{t}-{n}'] = fake(f't{t}', n, base, provider='p' if n else 'q')
                para_[k][f't{t}-{n}'] = fake(f't{t}', n, 0.0 if drop else base)
    r5 = q5(main_, para_, list(main_))
    got = {o['model']: o for o in r5['models']}
    if not all(got[f'm{m}']['p_holm'] < 0.05 and got[f'm{m}']['diff'] < 0 for m in range(3)):
        bad.append('q5 missed a paraphrase loss')
    if not all(got[f'm{m}']['diff'] == 0 and got[f'm{m}']['p'] == 1.0 for m in range(3, 11)):
        bad.append('q5 found a loss that is not there')
    if not (r5['tau']['tau'] > 0.6 and r5['tau']['items'] == 450):
        bad.append(f"q5 tau {r5['tau']}")
    if not (all(got[f'm{m}']['within_margin'] for m in range(3, 11)) and not any(got[f'm{m}']['within_margin'] for m in range(3))):
        bad.append('q5 bounds: a zero difference must sit within the margin and a 9-point loss outside it')
    if not np.isfinite(r5['tau']['noise_arm']['median']):
        bad.append(f"q5 arm-size noise floor {r5['tau']['noise_arm']}")
    if got['m0']['same_provider']['items'] != 300 or got['m0']['e3']['diff'] != 0.0:
        bad.append(f"q5 same-provider or e3 columns {got['m0']['same_provider']} {got['m0']['e3']}")
    # Q3's judged columns on two templates of a solved and a wrong trace each, answers worked by hand;
    # then the same with one call unanswered, which must leave the rates and be counted (D-148)
    rows3 = []
    for t in (0, 1):
        for n, sc in ((0, 1.0), (1, 0.0)):
            r = fake(f'u{t}', n, sc)
            r.update(status='answered', steps=[{'digit_flags': int(t == 0 and n == 1), 'tol1_flags': 0}, {'digit_flags': 0, 'tol1_flags': 0}])
            rows3.append(r)
    e5r = {'u0-0': dict(milestones_required=2, by_e3=2, sent=False, reached=0, missing=0, unjudged=0, e5_strict=1.0, reply_ok=None),
           'u1-0': dict(milestones_required=2, by_e3=2, sent=False, reached=0, missing=0, unjudged=0, e5_strict=1.0, reply_ok=None),
           'u0-1': dict(milestones_required=4, by_e3=1, sent=True, reached=1, missing=2, unjudged=0, e5_strict=0.5, reply_ok=True),
           'u1-1': dict(milestones_required=4, by_e3=1, sent=True, reached=0, missing=2, unjudged=1, e5_strict=0.25, reply_ok=True)}
    rtr = {'u0-0': dict(router_flagged=[1], judge_flagged=[1], sent=[0, 1], unjudged=0, reply_ok=True),
           'u1-0': dict(router_flagged=[], judge_flagged=[], sent=[0], unjudged=1, reply_ok=True),
           'u0-1': dict(router_flagged=[0], judge_flagged=[], sent=[1], unjudged=0, reply_ok=True),
           'u1-1': dict(router_flagged=[], judge_flagged=[], sent=[0, 1], unjudged=0, reply_ok=True)}
    s3 = stage_q3(rows3, e5r, rtr, 7)
    want3 = {'e5_coverage_on_wrong': 0.375, 'e5_judged_fraction': 0.5, 'e5_unjudged_rate': 1 / 6,
             'e5_reached_share': 1 / 6, 'router_rate_on_fully_solved': 0.5, 'router_judge_rate_on_fully_solved': 0.5,
             'router_rate_on_wrong': 0.5, 'router_unjudged_rate': 1 / 6, 'router_steps_flagged_per_trace': 0.5}
    if any(abs(s3[k] - v) > 1e-12 for k, v in want3.items()) or \
            s3['attribution_on_wrong'] != {'traces': 2, 'digit_rule': 0.5, 'e5_missing': 1.0, 'router_judge': 0.0} or \
            s3['e5_incomplete'] or s3['router_incomplete']:
        bad.append(f'stage_q3 {s3}')
    e5x = dict(e5r, **{'u1-1': dict(milestones_required=4, by_e3=1, sent=True, reached=0, missing=0, unjudged=3, e5_strict=0.25, reply_ok=False)})
    rtx = dict(rtr, **{'u1-1': dict(router_flagged=[], judge_flagged=[], sent=[0, 1], unjudged=2, reply_ok=False)})
    s3x = stage_q3(rows3, e5x, rtx, 7)
    if not (s3x['e5_incomplete'] and s3x['e5_without_reply'] == 1 and abs(s3x['e5_coverage_on_wrong'] - 0.5) < 1e-12
            and s3x['e5_wrong_used'] == 1 and s3x['router_incomplete'] and s3x['router_without_reply'] == 1
            and abs(s3x['router_rate_on_wrong'] - 1.0) < 1e-12 and s3x['attribution_on_wrong']['traces'] == 1):
        bad.append(f'stage_q3 did not leave the unanswered call out: {s3x}')
    keep = {i for i in main_['m0'] if not i.startswith('t0-')}          # the experts reject template t0
    r5k = q5(main_, para_, list(main_), keep)
    if not all(o['items'] == 447 for o in r5k['models']) or r5k['expert_check'] != 447:
        bad.append('q5 did not drop the rejected pairs from both arms')
    # first-flag position and the milestone buckets on a hand case
    rf = [fake('v', 0, 1.0), fake('v', 1, 1.0), fake('v', 2, 0.0)]
    rf[0]['steps'] = [{'digit_flags': 0, 'tol1_flags': 0}, {'digit_flags': 1, 'tol1_flags': 0}, {'digit_flags': 0, 'tol1_flags': 0}]
    rf[1]['steps'] = [{'digit_flags': 1, 'tol1_flags': 1}]
    rf[2]['milestones_required'] = 7
    ff = first_flag_position(rf)
    if ff['traces'] != 2 or ff['median'] != 0.25:
        bad.append(f'first flag position {ff}')
    bm = by_milestone_count(rf)
    if bm['2']['wrong_rate'] != 0.0 or bm['6+']['wrong_rate'] != 1.0 or bm['0']['items'] != 0:
        bad.append(f'milestone buckets {bm}')
    # the template-matched provider difference: provider q is 0.2 better on every template it shares
    runs_p = {'m': {}}
    for t in range(20):
        for n in range(4):
            runs_p['m'][f'w{t}-{n}'] = fake(f'w{t}', n, 1.0 if n == 0 else 0.0, provider='q' if n == 0 else 'p')
    pt = provider_table(runs_p, ['m'])['m']
    if abs(pt['q']['matched_diff'] - 1.0) > 1e-12 or abs(pt['p']['matched_diff'] + 1.0) > 1e-12 or pt['q']['templates_matched'] != 20:
        bad.append(f'provider table {pt}')
    print('selftest:', 'all pass' if not bad else bad)
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())
