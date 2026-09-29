"""The analysis plan (ANALYSIS_PLAN.md, D-117), computed from the score store.

    python -m full_run_28092026.analyze              # writes results/RESULTS.md and results/results.json
    python -m full_run_28092026.analyze --selftest   # the statistics on synthetic data; writes nothing

Reads scores/main/<model>.jsonl (score.py) for the eleven roster models, and scores/paraphrase/ and
scores/repeat1..3/ once those runs exist. The set-aside Qwen3.8-27B (D-132) is reported on its own
line, outside every comparison. Nothing here reads a trace or calls a model.

WHAT THE PLAN FIXES. The unit (the template), the intervals (template bootstrap, B = 10,000,
percentile), the scoring (correct 1, partial 0.5, incorrect and unusable 0), the tests and Holm's
correction within each family. The script stops unless the store holds what the plan describes:
150 templates x 15 items per model, 58 Easy and 34 Advanced templates, 58 single-path ones.

WHAT THIS SCRIPT SETS WHERE THE PLAN IS SILENT, written before any result was read:
  - Q2's p-value: a permutation test of the tier labels among the 92 Easy and Advanced templates
    (10,000 shuffles), beside the plan's within-tier bootstrap interval.
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
"""
from __future__ import annotations

import argparse
import collections
import itertools
import json
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
B = 10_000
ROSTER = ['gpt-oss-20b', 'gemma-4-26b-a4b', 'deepseek-v4.1-flash', 'qwen3-235b-a22b-2507', 'glm-5.3-flash',
          'glm-5.3', 'muse-glimmer-30b', 'kimi-k3', 'gpt-5.4-mini', 'gemini-3.1-flash-lite', 'claude-sonnet-5']
SET_ASIDE = ['qwen3.8-27b']
# the four templates the pilot excluded for shortcuts: D-057 (two), D-046, D-066
SHORTCUT = ['template_system_property_linearity', 'template_system_properties_memory_causality',
            'template_line_balancing_heuristic', 'template_levenspiel_plot_interpretation']
FROUDE = 'template_critical_depth_froude_classification'
REPEATS = ['repeat1', 'repeat2', 'repeat3']
PLAN_COUNTS = {'templates': 150, 'per_template': 15, 'Easy': 58, 'Advanced': 34, 'single_path': 58}


# ------------------------------------------------------------------ loading

def load(variant: str, key: str) -> dict[str, dict] | None:
    p = SCORES / variant / f'{key}.jsonl'
    if not p.exists():
        return None
    return {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())}


def fully(r: dict) -> float:
    return 1.0 if (not r['unusable'] and r['answer']['label'] == 'correct') else 0.0


def score_with(r: dict, field: str) -> float:
    if r['unusable']:
        return 0.0
    return {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}[r['answer'][field]]


def flagged(r: dict) -> float:
    return 1.0 if any(s['digit_flags'] for s in r['steps']) else 0.0


def verdict3(r: dict) -> str:
    return 'unusable' if r['unusable'] else r['answer']['label']


def per_template(rows, fn, keep=lambda r: True) -> dict[str, list[float]]:
    out = collections.defaultdict(list)
    for r in rows:
        if keep(r):
            out[r['template_id']].append(fn(r))
    return out


def template_matrix(runs: dict, fn, templates: list[str]) -> np.ndarray:
    """models x templates: the mean of fn(row) over each template's items."""
    M = np.zeros((len(runs), len(templates)))
    for i, rows in enumerate(runs.values()):
        by_t = per_template(rows.values(), fn)
        M[i] = [np.mean(by_t[t]) for t in templates]
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
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def boot_mean(v: np.ndarray, seed: int) -> list[float]:
    """Percentile interval of a mean over templates, templates resampled."""
    rng = np.random.default_rng(seed)
    return pct(v[rng.integers(0, len(v), size=(B, len(v)))].mean(axis=1))


def cluster_mean(groups: dict[str, list[float]], seed: int) -> tuple[float, list[float], int]:
    """The item mean over a subset of traces, with the template bootstrap of that ratio."""
    if not groups:
        return float('nan'), [float('nan'), float('nan')], 0
    tot = np.array([sum(v) for v in groups.values()])
    cnt = np.array([len(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(tot), size=(B, len(tot)))
    return float(tot.sum() / cnt.sum()), pct(tot[idx].sum(axis=1) / cnt[idx].sum(axis=1)), int(cnt.sum())


def sign_flip_p(d: np.ndarray, seed: int) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(d.mean())
    null = np.abs((rng.choice([-1.0, 1.0], size=(B, len(d))) * d).mean(axis=1))
    return float((np.sum(null >= obs - 1e-12) + 1) / (B + 1))


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


def tier_perm_p(x: np.ndarray, is_easy: np.ndarray, seed: int) -> float:
    """Easy minus Advanced, the tier labels shuffled among the templates of both tiers."""
    rng = np.random.default_rng(seed)
    obs = abs(x[is_easy].mean() - x[~is_easy].mean())
    n_e = int(is_easy.sum())
    px = x[np.argsort(rng.random((B, len(x))), axis=1)]
    null = np.abs(px[:, :n_e].mean(axis=1) - px[:, n_e:].mean(axis=1))
    return float((np.sum(null >= obs - 1e-12) + 1) / (B + 1))


def boot_gap(xe: np.ndarray, xa: np.ndarray, seed: int) -> list[float]:
    """Easy minus Advanced, templates resampled within each tier."""
    rng = np.random.default_rng(seed)
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
        split = collections.Counter(verdict3(r) for r in runs[k].values())
        models.append({'model': k, 'score': float(S[i].mean()), 'ci': boot_mean(S[i], 100 + i),
                       'fully_solved': float(F[i].mean()), 'fully_ci': boot_mean(F[i], 200 + i),
                       **{v: split[v] for v in ('correct', 'partial', 'incorrect', 'unusable')}})
    ids = sorted(runs[keys[0]])
    V = {k: np.array([fully(runs[k][x]) for x in ids]) for k in keys}
    pairs = []
    for n, (i, j) in enumerate(itertools.combinations(range(len(keys)), 2)):
        d, df = S[i] - S[j], F[i] - F[j]
        pairs.append({'a': keys[i], 'b': keys[j], 'diff': float(d.mean()), 'ci': boot_mean(d, 1000 + n),
                      'p': sign_flip_p(d, 2000 + n), 'fully_diff': float(df.mean()),
                      'fully_ci': boot_mean(df, 3000 + n), 'mcnemar_p': mcnemar_p(V[keys[i]], V[keys[j]])})
    for o, ph, qh in zip(pairs, holm([p['p'] for p in pairs]), holm([p['mcnemar_p'] for p in pairs])):
        o['p_holm'], o['mcnemar_p_holm'] = ph, qh
        s1, s2 = ph < 0.05, qh < 0.05
        o['same_verdict'] = bool((s1 == s2) and (not s1 or np.sign(o['diff']) == np.sign(o['fully_diff'])))
    return {'models': models, 'pairs': pairs}


def q2(runs, templates, levels, keys):
    S = template_matrix(runs, lambda r: r['score'], templates)
    tiers = np.array([levels[t] for t in templates])
    both = (tiers == 'Easy') | (tiers == 'Advanced')
    out = []
    for i, k in enumerate(keys):
        xe, xa = S[i][tiers == 'Easy'], S[i][tiers == 'Advanced']
        sigma = float(np.sqrt(((len(xe) - 1) * xe.var(ddof=1) + (len(xa) - 1) * xa.var(ddof=1))
                              / (len(xe) + len(xa) - 2)))
        out.append({'model': k, 'easy': float(xe.mean()), 'advanced': float(xa.mean()),
                    'gap': float(xe.mean() - xa.mean()), 'ci': boot_gap(xe, xa, 400 + i),
                    'p': tier_perm_p(S[i][both], tiers[both] == 'Easy', 500 + i), 'sigma': sigma,
                    'detectable': 2.80 * sigma * float(np.sqrt(1 / len(xe) + 1 / len(xa)))})
    for o, ph in zip(out, holm([o['p'] for o in out])):
        o['p_holm'] = ph
    return out


def q3(runs, keys):
    wrong = lambda r: r['score'] == 0.0
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
        fw, fw_ci, n_wrong = cluster_mean(per_template(rows, flagged, wrong), 750 + i)
        answered = [r for r in rows if r['status'] == 'answered']
        claims = [sum(s['claims'] for s in r['steps']) for r in answered]
        out.append({'model': k, 'claims_per_trace': float(np.mean(claims)),
                    'traces_with_a_claim': float(np.mean([c > 0 for c in claims])),
                    'wrong': n_wrong, 'wrong_unusable': sum(wrong(r) and r['unusable'] for r in rows),
                    'wrong_with_milestones': n_cov, 'e3_coverage_on_wrong': cov, 'e3_ci': cov_ci,
                    'readable_wrong_with_milestones': n_rcov, 'e3_coverage_on_readable_wrong': rcov,
                    'e3_readable_ci': rcov_ci,
                    'fully_solved': n_right, 'digit_flag_rate_on_fully_solved': fc, 'digit_ci': fc_ci,
                    'digit_flag_rate_on_wrong': fw, 'digit_wrong_ci': fw_ci,
                    'e5_coverage_on_wrong': None, 'e5_judged_fraction': None, 'e5_unjudged_rate': None})
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


def q5(main, para, keys):
    """Paired paraphrase minus original, per model, on the items both arms hold."""
    out, tested = [], []
    for i, k in enumerate(keys):
        if not para.get(k):
            continue
        ids = sorted(set(para[k]) & set(main[k]))
        diff = per_template([{'template_id': main[k][x]['template_id'],
                              'd': para[k][x]['score'] - main[k][x]['score']} for x in ids], lambda r: r['d'])
        point, c, n = cluster_mean(diff, 900 + i)
        out.append({'model': k, 'items': n, 'templates': len(diff), 'diff': point, 'ci': c,
                    'p': sign_flip_p(np.array([sum(v) for v in diff.values()]), 950 + i),
                    'mcnemar_p': mcnemar_p(np.array([fully(para[k][x]) for x in ids]),
                                           np.array([fully(main[k][x]) for x in ids]))})
        tested.append(k)
    for o, ph, qh in zip(out, holm([o['p'] for o in out]), holm([o['mcnemar_p'] for o in out])):
        o['p_holm'], o['mcnemar_p_holm'] = ph, qh
    tau = None
    if len(tested) >= 3:
        common = sorted(set.intersection(*(set(para[k]) & set(main[k]) for k in tested)))
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
        tau = {**kendall_boot(S, N, 990), 'items': len(common)}
    return {'models': out, 'tau': tau}


def sensitivity(runs, templates, keys):
    keep = set(templates) - set(SHORTCUT)
    out = []
    for k in keys:
        rows = list(runs[k].values())
        usable = [r for r in rows if not r['unusable']]
        out.append({'model': k,
                    'half_tol': float(np.mean([score_with(r, 'label_half_tol') for r in rows])),
                    'fitted': float(np.mean([r['score'] for r in rows])),
                    'double_tol': float(np.mean([score_with(r, 'label_double_tol') for r in rows])),
                    'fully_solved': float(np.mean([fully(r) for r in rows])),
                    'unusable_excluded': float(np.mean([r['score'] for r in usable])),
                    'usable': len(usable),
                    'without_shortcut_templates': float(np.mean([r['score'] for r in rows
                                                                 if r['template_id'] in keep]))})
    head = [o['fitted'] for o in out]
    tau = {c: float(stats.kendalltau(head, [o[c] for o in out]).statistic)
           for c in ('half_tol', 'double_tol', 'fully_solved', 'unusable_excluded', 'without_shortcut_templates')}
    return {'models': out, 'tau_with_headline': tau}


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


def reported(runs, keys):
    """Also reported, not tested: branch, domain, answer type, tokens against score, per-label accuracy."""
    gold_label = {}
    for rows in runs.values():                      # the label is the gold's; any model's row carries it
        for x, r in rows.items():
            if r['answer_type'] == 'classification' and r['answer'] and x not in gold_label:
                gold_label[x] = class_label(r)
    cls = collections.defaultdict(list)
    for x, lab in gold_label.items():
        cls[(runs[keys[0]][x]['template_id'], lab)].append(x)
    out = {'models': {}, 'classification': {}}
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


# ------------------------------------------------------------------ report

def ci(c) -> str:
    return f'{c[0]:.3f} to {c[1]:.3f}'


def render(res) -> str:
    m1, pairs = res['q1']['models'], res['q1']['pairs']
    sig = sum(p['p_holm'] < 0.05 for p in pairs)
    same = sum(p['same_verdict'] for p in pairs)
    differ = [p for p in pairs if not p['same_verdict']]
    only_mcnemar = sum(p['mcnemar_p_holm'] < 0.05 <= p['p_holm'] for p in differ)
    mf = res['milestones']
    L = ['# Results of the full run: the deterministic stack', '',
         'Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the '
         'method choices the plan leaves open are fixed in the script\'s docstring. The answer check, E3 and '
         "E4's digit rule only: E5's judge has not run, so its columns are empty. Every interval is 95% and "
         'resamples the 150 templates (B = 10,000).', '',
         '## Q1. Answer score per model', '',
         'Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0.', '',
         '| model | answer score | 95% CI | fully solved | 95% CI | correct | partial | incorrect | unusable |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in sorted(m1, key=lambda r: -r['score']):
        L.append(f"| `{r['model']}` | {r['score']:.3f} | {ci(r['ci'])} | {r['fully_solved']:.3f} | {ci(r['fully_ci'])} | "
                 f"{r['correct']} | {r['partial']} | {r['incorrect']} | {r['unusable']} |")
    L += ['', f'Of the 55 pairs, {sig} differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip '
          'permutation over the 150 per-template differences). The fully-solved check (McNemar exact, Holm) '
          f'gives the same verdict on {same} of 55: both tests or neither hold, and when both hold, in the same '
          f'direction. On {only_mcnemar} of the {len(differ)} others McNemar holds and the template-level test does '
          'not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so '
          'a claim rests on the template-level test.', '',
          '| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | McNemar p (Holm) | same verdict |',
          '|---|---|---:|---:|---:|---:|---:|---:|---|']
    for p in sorted(pairs, key=lambda p: (p['p_holm'], -abs(p['diff']))):
        L.append(f"| `{p['a']}` | `{p['b']}` | {p['diff']:+.3f} | {ci(p['ci'])} | {p['p_holm']:.4f} | "
                 f"{p['fully_diff']:+.3f} | {ci(p['fully_ci'])} | {p['mcnemar_p_holm']:.4f} | "
                 f"{'yes' if p['same_verdict'] else 'no'} |")
    L += ['', '## Q2. The complexity cliff: Easy minus Advanced', '',
          'Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier; '
          'p from permuting the tier labels, Holm across the eleven models. Sigma is the within-tier SD of the '
          'template-mean score; the detectable gap, 2.80 x sigma x sqrt(1/58 + 1/34), is the smallest this '
          'test finds at 80% power.', '',
          '| model | Easy | Advanced | gap | 95% CI | p (Holm) | sigma | detectable gap |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in res['q2']:
        L.append(f"| `{r['model']}` | {r['easy']:.3f} | {r['advanced']:.3f} | {r['gap']:+.3f} | {ci(r['ci'])} | "
                 f"{r['p_holm']:.4f} | {r['sigma']:.3f} | {r['detectable']:.3f} |")
    L += ['', '## Q3. What the process scores add beyond the answer', '',
          'Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in '
          f"neither. E3 coverage leaves out the {mf['items_without']} items with no milestones (all 15 of "
          f"{mf['templates_all_without']} templates and some items of {mf['templates_some_without']} more); "
          f"{mf['templates_with_item_at_most_one']} templates have an item with at most one milestone "
          f"({mf['templates_with_item_exactly_one']} with exactly one), where coverage is close to an answer check. "
          'E3 is the deterministic part of E5, which adds the judge\'s verdict on the milestones E3 does not find. '
          "The digit rule's flag counts a trace when any step is flagged. Against the experts, on fully solved traces, "
          'it had precision 0.750 and recall 0.320 (SCORER_VALIDATION.md), so its rate is not a count of '
          'slips; and it reads a different amount of arithmetic in each model\'s traces (claims checked per '
          'answered trace, shown), so a low rate can mean little was read.', '',
          '**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it '
          'stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.', '',
          '| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | '
          'E5 coverage |', '|---|---:|---:|---:|---:|---:|---:|---|']
    for r in res['q3']:
        L.append(f"| `{r['model']}` | {r['wrong']} | {r['wrong_unusable']} | {r['e3_coverage_on_wrong']:.3f} | "
                 f"{ci(r['e3_ci'])} | {r['e3_coverage_on_readable_wrong']:.3f} | {ci(r['e3_readable_ci'])} | pending |")
    L += ['', '**The digit rule.**', '',
          '| model | flag rate, wrong answers | 95% CI | fully solved traces | flag rate, fully solved | 95% CI | '
          'claims checked per trace | traces with a claim |', '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in res['q3']:
        L.append(f"| `{r['model']}` | {r['digit_flag_rate_on_wrong']:.3f} | {ci(r['digit_wrong_ci'])} | "
                 f"{r['fully_solved']} | {r['digit_flag_rate_on_fully_solved']:.3f} | {ci(r['digit_ci'])} | "
                 f"{r['claims_per_trace']:.2f} | {r['traces_with_a_claim']:.3f} |")
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
        L += ['| model | items | paraphrase - original | 95% CI | p (Holm) | McNemar p (Holm) |', '|---|---:|---:|---:|---:|---:|']
        L += [f"| `{r['model']}` | {r['items']} | {r['diff']:+.3f} | {ci(r['ci'])} | {r['p_holm']:.4f} | "
              f"{r['mcnemar_p_holm']:.4f} |" for r in res['q5']['models']]
        if res['q5']['tau']:
            t = res['q5']['tau']
            L += ['', f"Kendall's tau between the models' answer scores on the originals and on the paraphrases, over "
                  f"the {t['items']} items every tested model holds: {t['tau']:.3f}, 95% CI {ci(t['ci'])}."]
    else:
        L.append('Not run yet.')
    sens = res['sensitivity']
    L += ['', '## Sensitivity', '',
          "Answer score under each variation; the last row is Kendall's tau between that ordering of the models "
          'and the headline one.', '',
          '| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | '
          'without the 4 shortcut templates |', '|---|---:|---:|---:|---:|---:|---:|']
    for r in sens['models']:
        L.append(f"| `{r['model']}` | {r['half_tol']:.3f} | {r['fitted']:.3f} | {r['double_tol']:.3f} | "
                 f"{r['fully_solved']:.3f} | {r['unusable_excluded']:.3f} | {r['without_shortcut_templates']:.3f} |")
    t = sens['tau_with_headline']
    L += [f"| tau with the headline | {t['half_tol']:.3f} | 1 | {t['double_tol']:.3f} | {t['fully_solved']:.3f} | "
          f"{t['unusable_excluded']:.3f} | {t['without_shortcut_templates']:.3f} |", '',
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
    para = {k: load('paraphrase', k) for k in ROSTER} if a.store == 'main' else {}
    res = {'q1': q1(runs, templates, ROSTER), 'q2': q2(runs, templates, levels, ROSTER), 'q3': q3(runs, ROSTER),
           'q4': q4(runs, templates, single, ROSTER), 'q5': q5(runs, para, ROSTER),
           'sensitivity': sensitivity(runs, templates, ROSTER), 'reported': reported(runs, ROSTER),
           'repeats': repeats(runs, ROSTER) if a.store == 'main' else {},
           'set_aside': set_aside(a.store), 'milestones': milestone_facts(runs[ROSTER[0]].values()),
           'store': {'name': a.store, **json.loads((SCORES / a.store / 'CONFIG.json').read_text(encoding='utf-8'))}}
    out = OUT if a.store == 'main' else OUT / a.store
    out.mkdir(parents=True, exist_ok=True)
    (out / 'results.json').write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8', newline='\n')
    text = render(res)
    if a.store != 'main':
        text = text.replace('\n', f'\n\n**A record, not the results:** computed from `scores/{a.store}`. {a.note}\n', 1)
    (out / 'RESULTS.md').write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0


# ------------------------------------------------------------------ self-test

def fake(tid, n, score, level='Easy', single=True):
    lab = {1.0: 'correct', 0.5: 'partial', 0.0: 'incorrect'}[score]
    return {'item_id': f'{tid}-{n}', 'template_id': tid, 'score': score, 'unusable': False,
            'answer': {'label': lab, 'label_half_tol': lab, 'label_double_tol': lab, 'targets': None},
            'e3': {'coverage': 0.5}, 'steps': [], 'level': level, 'single_path': single}


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
    # Q1: a model against itself differs by nothing
    rows = {f't{t:03d}-{n}': fake(f't{t:03d}', n, float(rng.integers(0, 2))) for t in range(150) for n in range(15)}
    r1 = q1({'a': rows, 'b': dict(rows)}, sorted({r['template_id'] for r in rows.values()}), ['a', 'b'])['pairs'][0]
    if r1['diff'] != 0 or r1['p'] != 1.0 or r1['mcnemar_p'] != 1.0:
        bad.append(f'q1 self-comparison {r1}')
    # Q5: eleven models, three items per template; models 0-2 lose some solved items on the paraphrase
    main_, para_ = {}, {}
    for m in range(11):
        k = f'm{m}'
        main_[k], para_[k] = {}, {}
        for t in range(150):
            for n in range(3):
                base = float(rng.random() < 0.3 + 0.05 * m)
                drop = m < 3 and base == 1.0 and rng.random() < 0.3
                main_[k][f't{t}-{n}'] = fake(f't{t}', n, base)
                para_[k][f't{t}-{n}'] = fake(f't{t}', n, 0.0 if drop else base)
    r5 = q5(main_, para_, list(main_))
    got = {o['model']: o for o in r5['models']}
    if not all(got[f'm{m}']['p_holm'] < 0.05 and got[f'm{m}']['diff'] < 0 for m in range(3)):
        bad.append('q5 missed a paraphrase loss')
    if not all(got[f'm{m}']['diff'] == 0 and got[f'm{m}']['p'] == 1.0 for m in range(3, 11)):
        bad.append('q5 found a loss that is not there')
    if not (r5['tau']['tau'] > 0.6 and r5['tau']['items'] == 450):
        bad.append(f"q5 tau {r5['tau']}")
    print('selftest:', 'all pass' if not bad else bad)
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())
