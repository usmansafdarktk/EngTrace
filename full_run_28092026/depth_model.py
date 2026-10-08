"""Wrong answers against derivation depth, controlled for answer kind (WS-C3 step 3; review 1 M6).

    python -m full_run_28092026.depth_model [--quick]   # writes results/depth_model.json and results/sections/depth_model.md
    python -m full_run_28092026.depth_model --selftest  # the model on synthetic depth tables; writes nothing

WHAT IT ASKS. Does a model answer wrongly more often the more milestones the gold derivation has, once the kind
of answer (scalar, multipart, vector, array, symbolic, classification) is held fixed? analyze.py's Q3 prints the
raw wrong-answer rate by milestone count over every row; that rate mixes in empty and unreadable responses, whose
share also moves with depth (longer derivations hit output ceilings), and the answer kinds, which are not spread
evenly over depth.

THE CHOICES, fixed before the results were read:
  - Rows: readable responses (an empty or unreadable one is left out) whose item has at least one milestone; the
    milestone count is the gold trace's (`milestones_required`, the items' milestones as score.py derived them).
  - Wrong: the answer check's label is incorrect, as analyze.py's Q3 counts a wrong answer (a score of 0); a
    partial answer is not wrong.
  - Per bin of milestone count (1, 2, 3, 4-5, 6+): the item mean of wrong with a template bootstrap interval
    (percentile; the templates with a row in that bin resampled, the ratio of sums), and the rows.
  - Per model: a logistic regression of wrong on the milestone count (linear in the log odds) with answer kind as
    a categorical covariate (scalar the reference), fitted by maximum likelihood (statsmodels GLM, binomial), its
    covariance cluster-robust over templates. The slope is the change in log odds per milestone; its interval is
    the Wald 95% interval from the clustered standard error, and p is the two-sided Wald test. Holm's correction
    runs over the models of one configuration. A slope "holds" when its Holm-adjusted p is below 0.05.
  - Separation: where a model has no wrong answer (or only wrong answers) on one answer kind, that kind's
    coefficient has no finite estimate. Its rows are then left out of that model's fit, which gives the slope the
    full likelihood tends to (the separated rows contribute nothing to it in the limit); the kinds left out are
    recorded per model.
  - Pooled: the same model over every model's rows with a model fixed effect, clustered by template, as one
    summary beside the per-model slopes (outside the Holm family).
Configurations: main (the eleven roster models at provider defaults) and matched (each model at its reasoning
store from results/matched_config.json where it has one; the sibling with `run: false` skipped).

Reads only the score stores. Calls no model. `--quick` (or ENGTRACE_QUICK=1) draws 1,000 resamples for the bin
intervals instead of 10,000.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026.analyze import holm  # noqa: E402
from full_run_28092026.clause_variants import load_rows, matched_config, quick_mode, stores  # noqa: E402

RESULTS = HERE / 'results'
OUT_JSON = RESULTS / 'depth_model.json'
OUT_MD = RESULTS / 'sections' / 'depth_model.md'
BINS = (('1', 1, 1), ('2', 2, 2), ('3', 3, 3), ('4-5', 4, 5), ('6+', 6, 10 ** 6))
REFERENCE = 'scalar'
Z = 1.959963984540054


def depth_rows(rows: list[dict]) -> list[dict]:
    """The rows the model reads: readable, with at least one milestone; wrong = incorrect."""
    return [{'template_id': r['template_id'], 'answer_type': r['answer_type'], 'm': int(r['milestones_required']),
             'wrong': 1.0 if r['answer']['label'] == 'incorrect' else 0.0}
            for r in rows if not r['unusable'] and r['milestones_required'] >= 1]


def ratio_boot(groups: dict[str, list[float]], draws: int, seed: int) -> list[float]:
    tot = np.array([sum(v) for v in groups.values()])
    cnt = np.array([len(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(tot), size=(draws, len(tot)))
    d = tot[idx].sum(axis=1) / cnt[idx].sum(axis=1)
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def bins_of(rows: list[dict], draws: int, seed: int) -> dict:
    out = {}
    for k, (name, lo, hi) in enumerate(BINS):
        g = collections.defaultdict(list)
        for r in rows:
            if lo <= r['m'] <= hi:
                g[r['template_id']].append(r['wrong'])
        n = sum(len(v) for v in g.values())
        out[name] = {'rate': (sum(sum(v) for v in g.values()) / n) if n else None,
                     'ci': ratio_boot(g, draws, seed + k) if n else None, 'n': n, 'templates': len(g)}
    return out


def fit(rows: list[dict], extra: list[str] | None = None) -> dict:
    """Logistic regression of wrong on milestone count and answer kind (and `extra`, a categorical column such as
    the model), covariance clustered by template. Kinds with no variation in wrong are left out first."""
    import statsmodels.api as sm
    by_kind = collections.defaultdict(list)
    for r in rows:
        by_kind[r['answer_type']].append(r['wrong'])
    dropped = sorted(k for k, v in by_kind.items() if min(v) == max(v))
    use = [r for r in rows if r['answer_type'] not in dropped]
    wrong = sum(r['wrong'] for r in use)
    base = {'n': len(use), 'wrong': int(wrong), 'templates': len({r['template_id'] for r in use}),
            'kinds_dropped': dropped, 'm_range': [min(r['m'] for r in use), max(r['m'] for r in use)] if use else None}
    if not use or wrong in (0, len(use)):
        return {**base, 'slope': None, 'se': None, 'ci': None, 'p': None, 'odds_ratio': None, 'or_ci': None}
    kinds = sorted({r['answer_type'] for r in use} - {REFERENCE})
    cols = [np.ones(len(use)), np.array([r['m'] for r in use], float)]
    cols += [np.array([r['answer_type'] == k for r in use], float) for k in kinds]
    names = ['const', 'milestones'] + [f'kind[{k}]' for k in kinds]
    if REFERENCE not in {r['answer_type'] for r in use} and kinds:   # no reference rows: the first kind is it
        cols.pop(2)
        names.pop(2)
    for col in extra or []:
        levels = sorted({r[col] for r in use})
        for lv in levels[1:]:
            cols.append(np.array([r[col] == lv for r in use], float))
            names.append(f'{col}[{lv}]')
    X = np.column_stack(cols)
    y = np.array([r['wrong'] for r in use])
    groups = np.unique([r['template_id'] for r in use], return_inverse=True)[1]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit(cov_type='cluster', cov_kwds={'groups': groups})
    b, se = float(res.params[1]), float(res.bse[1])
    p = float(res.pvalues[1])
    ci = [b - Z * se, b + Z * se]
    return {**base, 'slope': b, 'se': se, 'ci': ci, 'p': p, 'odds_ratio': float(np.exp(b)),
            'or_ci': [float(np.exp(ci[0])), float(np.exp(ci[1]))], 'converged': bool(res.converged),
            'terms': names}


def configuration(name: str, model_stores: dict[str, str], draws: int) -> dict:
    models, pooled = {}, []
    for i, (k, s) in enumerate(model_stores.items()):
        rows = depth_rows(load_rows(s, k))
        m = fit(rows)
        m['bins'] = bins_of(rows, draws, 5100 + 10 * i)
        m['store'], m['rows'], m['wrong_rows'] = s, len(rows), int(sum(r['wrong'] for r in rows))
        models[k] = m
        pooled += [dict(r, model=k) for r in rows]
    keys = [k for k in models if models[k]['p'] is not None]
    for k, ph in zip(keys, holm([models[k]['p'] for k in keys])):
        models[k]['p_holm'] = ph
        models[k]['holds'] = bool(ph < 0.05)
    for k in models:
        models[k].setdefault('p_holm', None)
        models[k].setdefault('holds', None)
    return {'models': models, 'pooled': fit(pooled, extra=['model']), 'holm_over': len(keys),
            'holds': sum(bool(models[k]['holds']) for k in keys)}


def run(quick: bool) -> dict:
    draws = 1_000 if quick else 10_000
    if quick:
        print('QUICK: 1,000 resamples')
    plan = stores()
    main = {k: 'main' for k in plan['main']}
    cfg = matched_config()
    matched = {k: (cfg[k] if cfg.get(k) in plan and k in plan[cfg[k]] else 'main') for k in plan['main']}
    out = {'quick': quick, 'generated_by': 'python -m full_run_28092026.depth_model' + (' --quick' if quick else ''),
           'written_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
           'resamples': draws, 'configurations': {}}
    for name, ms in (('main', main), ('matched', matched)):
        out['configurations'][name] = configuration(name, ms, draws)
        c = out['configurations'][name]
        print(f"{name}: {c['holds']} of {c['holm_over']} slopes hold after Holm; pooled slope "
              f"{c['pooled']['slope']:.3f} ({c['pooled']['ci'][0]:.3f} to {c['pooled']['ci'][1]:.3f})")
    out['models'] = out['configurations']['main']['models']
    return out


# ------------------------------------------------------------------ writing

def render(res: dict) -> str:
    q = ' QUICK (1,000 resamples).' if res['quick'] else ''
    L = ['## Depth, controlled', '',
         f"Generated by `{res['generated_by']}`.{q} Readable responses whose item has at least one milestone. Wrong is "
         f"the incorrect label (a partial answer is not wrong). Per bin of the gold trace's milestone count: the "
         f"wrong-answer rate, its 95% template bootstrap interval ({res['resamples']:,} resamples) and the rows. The "
         f"slope is the change in log odds of a wrong answer per milestone in a logistic model with answer kind as a "
         f"covariate, its 95% interval and p from standard errors clustered by template; Holm over the models of the "
         f"configuration. Answer kinds on which a model has no wrong answer are left out of its fit (they have no "
         f"finite coefficient), listed in the last column.", '']
    for name, c in res['configurations'].items():
        title = 'Main run (provider defaults)' if name == 'main' else 'Matched settings (reasoning store where one exists)'
        L += [f'### {title}', '',
              '| model | rows | wrong | ' + ' | '.join(f'milestones {b}' for b, _, _ in BINS)
              + ' | slope (95% CI) | odds ratio | p | p Holm | rows fitted | kinds left out |',
              '|---|---:|---:|' + '---:|' * len(BINS) + '---:|---:|---:|---:|---:|---|']
        for k, m in c['models'].items():
            cells = []
            for b, _, _ in BINS:
                x = m['bins'][b]
                cells.append('-' if x['rate'] is None else
                             f"{x['rate']:.3f} [{x['ci'][0]:.3f}, {x['ci'][1]:.3f}] n={x['n']}")
            tag = f" ({m['store']})" if m['store'] != 'main' else ''
            if m['slope'] is None:
                fitted = '- | - | - | -'
            else:
                fitted = (f"{m['slope']:+.3f} ({m['ci'][0]:+.3f} to {m['ci'][1]:+.3f}) | {m['odds_ratio']:.2f} | "
                          f"{m['p']:.2g} | {m['p_holm']:.2g}{' (holds)' if m['holds'] else ''}")
            L.append(f"| `{k}`{tag} | {m['rows']} | {m['wrong_rows']} | " + ' | '.join(cells) + f" | {fitted} | "
                     f"{m['n']} | {', '.join(m['kinds_dropped']) or '-'} |")
        p = c['pooled']
        L += ['', f"{c['holds']} of {c['holm_over']} slopes hold after Holm. Pooled over the models (a model fixed "
              f"effect, clustered by template): slope {p['slope']:+.3f} (95% CI {p['ci'][0]:+.3f} to {p['ci'][1]:+.3f}), "
              f"odds ratio {p['odds_ratio']:.2f} per milestone, p = {p['p']:.2g}, over {p['n']} rows and {p['templates']} "
              f"templates.", '']
    return '\n'.join(L)


def write(res: dict) -> None:
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    OUT_MD.write_text(render(res), encoding='utf-8')
    print(f'wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}')


# ------------------------------------------------------------------ self-test

def synthetic(slope: float, seed: int, kinds=('scalar', 'multipart', 'symbolic')) -> list[dict]:
    """150 templates x 15 items: a template's milestone count from 1 to 10, its kind, a template effect, and a
    wrong answer drawn from the logistic model with the given slope."""
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(150):
        m0, kind, u = int(rng.integers(1, 11)), kinds[t % len(kinds)], rng.normal(0, 0.5)
        shift = {'scalar': 0.0, 'multipart': 0.5, 'symbolic': -0.5}.get(kind, 0.0)
        for i in range(15):
            m = max(1, m0 + int(rng.integers(-1, 2)))
            eta = -3.0 + slope * m + shift + u
            rows.append({'template_id': f't{t}', 'answer_type': kind, 'm': m,
                         'wrong': float(rng.random() < 1 / (1 + np.exp(-eta)))})
    return rows


def selftest() -> int:
    bad = 0
    # a hand-made depth table: the bins count and average exactly
    table = [('a', 1, 0), ('a', 1, 1), ('b', 2, 0), ('b', 3, 1), ('c', 4, 1), ('c', 5, 0), ('c', 7, 1), ('d', 6, 1)]
    rows = [{'template_id': t, 'answer_type': 'scalar', 'm': m, 'wrong': float(w)} for t, m, w in table]
    b = bins_of(rows, 200, 1)
    want = {'1': (0.5, 2), '2': (0.0, 1), '3': (1.0, 1), '4-5': (0.5, 2), '6+': (1.0, 2)}
    for k, (rate, n) in want.items():
        if abs(b[k]['rate'] - rate) > 1e-12 or b[k]['n'] != n:
            bad += 1
            print(f'FAIL bin {k}: {b[k]} want rate {rate} n {n}')
    # the rows the model reads: readable, at least one milestone, wrong = incorrect only
    store = [{'template_id': 't', 'answer_type': 'scalar', 'milestones_required': 3, 'unusable': False,
              'answer': {'label': 'incorrect'}},
             {'template_id': 't', 'answer_type': 'scalar', 'milestones_required': 3, 'unusable': False,
              'answer': {'label': 'partial'}},
             {'template_id': 't', 'answer_type': 'scalar', 'milestones_required': 0, 'unusable': False,
              'answer': {'label': 'incorrect'}},
             {'template_id': 't', 'answer_type': 'scalar', 'milestones_required': 3, 'unusable': True, 'answer': None}]
    if [r['wrong'] for r in depth_rows(store)] != [1.0, 0.0]:
        bad += 1
        print('FAIL depth_rows')
    # a real slope is recovered and detected; a null one is not detected
    pos = fit(synthetic(0.3, 7))
    if not (pos['ci'][0] < 0.3 < pos['ci'][1] and pos['p'] < 0.001):
        bad += 1
        print(f'FAIL slope 0.3: {pos["slope"]:.3f} {pos["ci"]} p {pos["p"]}')
    null = fit(synthetic(0.0, 8))
    if not (null['ci'][0] < 0.0 < null['ci'][1]):
        bad += 1
        print(f'FAIL null slope: {null["slope"]:.3f} {null["ci"]}')
    # separation: a kind with no wrong answer is left out, and the slope equals the fit without its rows
    sep = synthetic(0.3, 9) + [{'template_id': f's{t}', 'answer_type': 'classification', 'm': 1 + t % 9, 'wrong': 0.0}
                               for t in range(60)]
    a, b2 = fit(sep), fit([r for r in sep if r['answer_type'] != 'classification'])
    if a['kinds_dropped'] != ['classification'] or abs(a['slope'] - b2['slope']) > 1e-9:
        bad += 1
        print(f'FAIL separation: {a["kinds_dropped"]} {a["slope"]} vs {b2["slope"]}')
    # Holm marks the slopes that hold
    if holm([0.001, 0.04, 0.5]) != [0.003, 0.08, 0.5]:
        bad += 1
        print('FAIL holm')
    total = len(want) + 5
    print(f'{total - bad}/{total} checks pass')
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--quick', action='store_true', help='1,000 resamples (also ENGTRACE_QUICK=1)')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    write(run(quick_mode(a.quick)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
