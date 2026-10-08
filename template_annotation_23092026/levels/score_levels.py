"""E1 (WS-E): score the experts' difficulty ratings.

    python -m template_annotation_23092026.levels.score_levels                # FREE: experts_filled_levels/ ->
                                                                             #       agreement.json, RESULTS.md
    python -m template_annotation_23092026.levels.score_levels --returns DIR  # FREE: the returns from another folder
    python -m template_annotation_23092026.levels.score_levels --selftest     # FREE: simulated returns, temp folder

WHAT IT READS. The keyfile and the queues build_kit.py wrote to tasks/, and every <id>.jsonl under the returns
folder (any layout). A row whose code is not in the keyfile, or was not assigned to its expert, stops the run. The
last submission per expert and template counts (the app appends; a resubmission replaces).

WHAT IT COMPUTES, per branch and over all 150 templates:
  agreement      Fleiss' kappa on the three levels, and the weighted variant for ordered levels, with linear weights
                 (adjacent levels 0.5) and quadratic weights (adjacent 0.75): Gwet's generalisation (Handbook of
                 Inter-Rater Reliability, 4th ed., ch. 2), which reduces to Fleiss' kappa under identity weights, and
                 for two raters to Scott's pi. 95% intervals: percentile bootstrap over templates (2,000 resamples,
                 seeded; the overall interval resamples within branches).
  majority       per template, the level more than half of its raters chose (two of three); none when the three
                 differ. The templates whose majority differs from the current label, by one step or by two. The
                 majority against the current label: the share equal, and Cohen's kappa (unweighted and linear).
  formula        per template, the majority answer to "the problem states the governing formula or names the method";
                 the yes count per branch and current level (and per majority level); Fleiss' kappa on the answer.
WHAT IT WRITES. agreement.json (every figure above, and per template its rating counts and majority: anonymous
counts, no expert ids, no notes) and RESULTS.md (counts only). Nothing calls a model or the network.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import random
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from template_annotation_23092026.levels import build_kit  # noqa: E402

RETURNS = HERE / 'experts_filled_levels'
AGREEMENT = HERE / 'agreement.json'
LEVELS = build_kit.LEVELS
FORMULA = ['yes', 'no']
BOOT = 2000
BOOT_SEED = 20261007


# ------------------------------------------------------------------ statistics

def weights(kind: str, q: int = 3) -> list[list[float]]:
    if kind == 'identity':
        return [[1.0 if k == l else 0.0 for l in range(q)] for k in range(q)]
    if kind == 'linear':
        return [[1.0 - abs(k - l) / (q - 1) for l in range(q)] for k in range(q)]
    if kind == 'quadratic':
        return [[1.0 - (k - l) ** 2 / (q - 1) ** 2 for l in range(q)] for k in range(q)]
    raise ValueError(kind)


def kappa(rows: list[collections.Counter], cats: list[str], w: list[list[float]] | None = None) -> float | None:
    """Fleiss' kappa over items rated by two or more raters; with `w`, Gwet's weighted form of it.

    p_a = mean over items of sum_k r_ik (r*_ik - 1) / (r_i (r_i - 1)), r*_ik = sum_l w_kl r_il;
    pi_k = mean over items of r_ik / r_i;  p_e = sum_kl w_kl pi_k pi_l;  kappa = (p_a - p_e) / (1 - p_e).
    None when fewer than two items qualify or p_e is 1 (every rating in one category: kappa is 0/0)."""
    w = w or weights('identity', len(cats))
    rows = [r for r in rows if sum(r[c] for c in cats) >= 2]
    if len(rows) < 2:
        return None
    n = len(rows)
    pi = [sum(r[c] / sum(r[x] for x in cats) for r in rows) / n for c in cats]
    pa = 0.0
    for r in rows:
        ri = sum(r[c] for c in cats)
        s = sum(r[ck] * (sum(w[k][l] * r[cl] for l, cl in enumerate(cats)) - 1) for k, ck in enumerate(cats))
        pa += s / (ri * (ri - 1))
    pa /= n
    pe = sum(w[k][l] * pi[k] * pi[l] for k in range(len(cats)) for l in range(len(cats)))
    if pe >= 1 - 1e-12:
        return None
    return (pa - pe) / (1 - pe)


def cohen(pairs: list[tuple[str, str]], cats: list[str], w: list[list[float]] | None = None) -> float | None:
    """Cohen's kappa between two labellings (each with its own marginals); weighted with `w`."""
    w = w or weights('identity', len(cats))
    n = len(pairs)
    if n < 2:
        return None
    idx = {c: i for i, c in enumerate(cats)}
    po = sum(w[idx[a]][idx[b]] for a, b in pairs) / n
    pa = [sum(a == c for a, _ in pairs) / n for c in cats]
    pb = [sum(b == c for _, b in pairs) / n for c in cats]
    pe = sum(w[k][l] * pa[k] * pb[l] for k in range(len(cats)) for l in range(len(cats)))
    return None if pe >= 1 - 1e-12 else (po - pe) / (1 - pe)


def boot_ci(groups: list[list[collections.Counter]], cats: list[str], w, reps: int = BOOT,
            seed: int = BOOT_SEED) -> list[float] | None:
    """Percentile 95% interval of kappa, resampling items with replacement within each group."""
    rng = random.Random(seed)
    vals = []
    for _ in range(reps):
        sample = [g[rng.randrange(len(g))] for g in groups for _ in range(len(g))]
        k = kappa(sample, cats, w)
        if k is not None:
            vals.append(k)
    if len(vals) < reps * 0.9:
        return None
    vals.sort()
    return [round(vals[int(0.025 * len(vals))], 3), round(vals[min(len(vals) - 1, int(0.975 * len(vals)))], 3)]


def majority(c: collections.Counter, cats: list[str]) -> str | None:
    n = sum(c[x] for x in cats)
    top = max(cats, key=lambda x: c[x])
    return top if n and c[top] * 2 > n else None


def median_level(c: collections.Counter) -> str | None:
    seq = [lv for lv in LEVELS for _ in range(c[lv])]
    return seq[(len(seq) - 1) // 2] if seq else None


# ------------------------------------------------------------------ the returns

def read_returns(returns: Path, tasks: Path) -> tuple[dict, dict, dict]:
    keyfile = json.loads((tasks / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((tasks / 'assignment.json').read_text(encoding='utf-8'))
    owners = collections.defaultdict(set)
    for aid, q in queues.items():
        for c in q['codes']:
            owners[c].add(aid)
    rows = {}
    for f in sorted(returns.rglob('*.jsonl')):
        for line in f.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r['code'] not in keyfile:
                raise SystemExit(f"{f.name}: code {r['code']} is not in this build's keyfile")
            if r['annotator_id'] not in owners[r['code']]:
                raise SystemExit(f"{f.name}: {r['code']} was not assigned to {r['annotator_id']}")
            if r['level'] not in LEVELS or r['formula_stated'] not in FORMULA:
                raise SystemExit(f"{f.name}: {r['code']} has level {r['level']!r}, formula {r['formula_stated']!r}")
            rows[(r['annotator_id'], r['code'])] = r          # a later submission replaces an earlier one
    return keyfile, queues, rows


def summarise(keyfile: dict, queues: dict, rows: dict, build: dict) -> dict:
    lv_rows = collections.defaultdict(collections.Counter)
    fm_rows = collections.defaultdict(collections.Counter)
    for (aid, c), r in rows.items():
        lv_rows[c][r['level']] += 1
        fm_rows[c][r['formula_stated']] += 1
    branches = sorted({k['branch'] for k in keyfile.values()})
    W = {'fleiss': weights('identity'), 'weighted_linear': weights('linear'), 'weighted_quadratic': weights('quadratic')}

    def block(codes: list[str], groups: list[list[str]]) -> dict:
        items = [lv_rows[c] for c in codes]
        out = {'templates': sum(sum(lv_rows[c].values()) >= 2 for c in codes),
               'ratings': sum(sum(lv_rows[c].values()) for c in codes)}
        for name, w in W.items():
            k = kappa(items, LEVELS, w)
            out[name] = None if k is None else round(k, 3)
            if name != 'weighted_quadratic':
                out[name + '_ci95'] = boot_ci([[lv_rows[c] for c in g] for g in groups], LEVELS, w)
        f = kappa([fm_rows[c] for c in codes], FORMULA)
        out['formula_fleiss'] = None if f is None else round(f, 3)
        return out

    by_branch_codes = {b: sorted(c for c, k in keyfile.items() if k['branch'] == b) for b in branches}
    agreement = {'overall': block(sorted(keyfile), list(by_branch_codes.values())),
                 'by_branch': {b: block(cs, [cs]) for b, cs in by_branch_codes.items()}}

    templates = {}
    for c, k in sorted(keyfile.items(), key=lambda kv: kv[1]['template_id']):
        cnt = lv_rows[c]
        maj = majority(cnt, LEVELS)
        fmaj = majority(fm_rows[c], FORMULA)
        templates[k['template_id']] = {
            'branch': k['branch'], 'current': k['level'], 'ratings': {lv: cnt[lv] for lv in LEVELS},
            'majority': maj, 'median': median_level(cnt),
            'step': None if maj is None else abs(LEVELS.index(maj) - LEVELS.index(k['level'])),
            'formula': {x: fm_rows[c][x] for x in FORMULA}, 'formula_majority': fmaj}

    def maj_block(ts: list[dict]) -> dict:
        rated = [t for t in ts if sum(t['ratings'].values()) >= 2]
        with_maj = [t for t in rated if t['majority']]
        steps = collections.Counter(t['step'] for t in with_maj)
        pairs = [(t['current'], t['majority']) for t in with_maj]
        ck = cohen(pairs, LEVELS)
        ckl = cohen(pairs, LEVELS, weights('linear'))
        return {'templates_rated': len(rated), 'with_majority': len(with_maj), 'no_majority': len(rated) - len(with_maj),
                'same_as_current': steps[0], 'differs_by_one_step': steps[1], 'differs_by_two_steps': steps[2],
                'majority_levels': {lv: sum(t['majority'] == lv for t in with_maj) for lv in LEVELS},
                'agreement_with_current': round(steps[0] / len(with_maj), 3) if with_maj else None,
                'cohen_kappa_with_current': None if ck is None else round(ck, 3),
                'cohen_kappa_linear_with_current': None if ckl is None else round(ckl, 3)}

    tl = list(templates.values())
    majority_summary = {'overall': maj_block(tl), 'by_branch': {b: maj_block([t for t in tl if t['branch'] == b])
                                                               for b in branches}}

    def formula_block(ts: list[dict], level_key: str) -> dict:
        out = {}
        for lv in LEVELS:
            sel = [t for t in ts if t[level_key] == lv]
            out[lv] = {'templates': len(sel), 'yes': sum(t['formula_majority'] == 'yes' for t in sel),
                       'no': sum(t['formula_majority'] == 'no' for t in sel),
                       'no_majority': sum(t['formula_majority'] is None for t in sel)}
        return out

    formula = {'by_current_level': formula_block(tl, 'current'),
               'by_majority_level': formula_block([t for t in tl if t['majority']], 'majority'),
               'by_branch_current_level': {b: formula_block([t for t in tl if t['branch'] == b], 'current')
                                           for b in branches},
               'yes_templates': sum(t['formula_majority'] == 'yes' for t in tl)}
    dist = {'overall': {lv: sum(r['level'] == lv for r in rows.values()) for lv in LEVELS},
            'by_branch': {b: {lv: sum(r['level'] == lv for (a, c), r in rows.items() if keyfile[c]['branch'] == b)
                              for lv in LEVELS} for b in branches}}
    returned = sorted({aid for aid, _ in rows})
    complete = sum(sum(lv_rows[c].values()) == 3 for c in keyfile)
    return {'generated_by': 'template_annotation_23092026/levels/score_levels.py',
            'scored_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
            'build': {k: build.get(k) for k in ('built', 'seed', 'manifest_sha256')},
            'returns': {'ratings': len(rows), 'assigned': sum(len(q['codes']) for q in queues.values()),
                        'experts_returned': len(returned), 'experts': len(queues),
                        'experts_complete': sum(sum(1 for (a, _c) in rows if a == aid) == len(queues[aid]['codes'])
                                                for aid in returned),
                        'templates_with_three_ratings': complete},
            'method': {'kappa': "Fleiss' kappa; weighted: Gwet's generalisation, linear and quadratic weights",
                       'ci': f'percentile bootstrap over templates, {BOOT} resamples, seed {BOOT_SEED}; overall within branches',
                       'majority': 'the level chosen by more than half of the raters; none when the three differ'},
            'agreement': agreement, 'majority': majority_summary, 'formula_stated': formula,
            'rating_distribution': dist, 'templates': templates}


def fmt(x) -> str:
    return '-' if x is None else f'{x:.3f}'


def ci(x) -> str:
    return '-' if not x else f'[{x[0]:.3f}, {x[1]:.3f}]'


def result_lines(res: dict) -> list[str]:
    r = res['returns']
    L = ['## What came back', '',
         f"Ratings returned {r['ratings']} of {r['assigned']} assigned, from {r['experts_returned']} of {r['experts']} "
         f"experts ({r['experts_complete']} complete); templates with three ratings: {r['templates_with_three_ratings']} "
         'of 150.', '',
         '### Agreement on the level', '',
         "| | templates | Fleiss' kappa | 95% interval | weighted, linear | 95% interval | weighted, quadratic |",
         '|---|---:|---:|---|---:|---|---:|']
    a = res['agreement']
    for name, v in list(a['by_branch'].items()) + [('all branches', a['overall'])]:
        L.append(f"| {name.replace('_engineering', '')} | {v['templates']} | {fmt(v['fleiss'])} | {ci(v['fleiss_ci95'])} | "
                 f"{fmt(v['weighted_linear'])} | {ci(v['weighted_linear_ci95'])} | {fmt(v['weighted_quadratic'])} |")
    L += ['', '### The majority label against the current label', '',
          '| | rated | with a majority | same | one step apart | two steps apart | share same | Cohen\'s kappa (linear) |',
          '|---|---:|---:|---:|---:|---:|---:|---|']
    m = res['majority']
    for name, v in list(m['by_branch'].items()) + [('all branches', m['overall'])]:
        L.append(f"| {name.replace('_engineering', '')} | {v['templates_rated']} | {v['with_majority']} | "
                 f"{v['same_as_current']} | {v['differs_by_one_step']} | {v['differs_by_two_steps']} | "
                 f"{fmt(v['agreement_with_current'])} | {fmt(v['cohen_kappa_with_current'])} "
                 f"({fmt(v['cohen_kappa_linear_with_current'])}) |")
    ov = m['overall']['majority_levels']
    L += ['', f"Majority levels over the templates with a majority: Easy {ov['Easy']}, Intermediate {ov['Intermediate']}, "
              f"Advanced {ov['Advanced']}.", '']
    d = res['rating_distribution']['overall']
    L += [f"All ratings: Easy {d['Easy']}, Intermediate {d['Intermediate']}, Advanced {d['Advanced']}.", '',
          '### The problem states the governing formula or names the method', '',
          'Templates whose raters\' majority answered yes, by the current level (and by the majority level).', '',
          '| | Easy | Intermediate | Advanced |', '|---|---|---|---|']
    f = res['formula_stated']
    for b, v in f['by_branch_current_level'].items():
        L.append(f"| {b.replace('_engineering', '')} | " + ' | '.join(f"{v[lv]['yes']} of {v[lv]['templates']}"
                                                                     for lv in LEVELS) + ' |')
    L.append('| all branches, current level | ' + ' | '.join(f"{f['by_current_level'][lv]['yes']} of "
                                                         f"{f['by_current_level'][lv]['templates']}" for lv in LEVELS) + ' |')
    L.append('| all branches, majority level | ' + ' | '.join(f"{f['by_majority_level'][lv]['yes']} of "
                                                          f"{f['by_majority_level'][lv]['templates']}" for lv in LEVELS) + ' |')
    L += ['', f"Templates with a yes majority: {f['yes_templates']} of 150. Fleiss' kappa on the answer, all branches: "
              f"{fmt(a['overall']['formula_fleiss'])}.", '']
    return L


def score(returns: Path, tasks: Path, out_json: Path | None, out_md: Path | None) -> dict:
    keyfile, queues, rows = read_returns(returns, tasks)
    build = json.loads((tasks / 'build.json').read_text(encoding='utf-8'))
    res = summarise(keyfile, queues, rows, build)
    if out_json:
        out_json.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8', newline='\n')
    lines = build_kit.sent_lines(build) + result_lines(res)
    if out_md:
        out_md.write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(result_lines(res)))
    return res


# ------------------------------------------------------------------ selftest

def simulate(tasks: Path, out: Path, p_same: float, seed: int = 1, drop: str | None = None) -> int:
    """Synthetic returns (layer2/simulate.py's pattern): each rater keeps the current label with probability p_same,
    else moves one step; answers the formula question yes with probability 0.4 / 0.25 / 0.1 by level."""
    keyfile = json.loads((tasks / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((tasks / 'assignment.json').read_text(encoding='utf-8'))
    rng = random.Random(seed)
    t0 = dt.datetime(2026, 10, 8, 9, 0, tzinfo=dt.timezone.utc)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for aid, q in queues.items():
        if aid == drop:
            continue
        with (out / f'{aid}.jsonl').open('w', encoding='utf-8') as fh:
            for pos, c in enumerate(q['codes'], 1):
                cur = LEVELS.index(keyfile[c]['level'])
                lv = cur if rng.random() < p_same else min(2, max(0, cur + rng.choice([-1, 1])))
                yes = rng.random() < (0.4, 0.25, 0.1)[cur]
                opened = t0 + dt.timedelta(minutes=pos)
                fh.write(json.dumps({'annotator_id': aid, 'branch': q['branch'], 'code': c, 'position': pos,
                                     'opened_at': opened.isoformat(), 'level': LEVELS[lv],
                                     'formula_stated': 'yes' if yes else 'no', 'note': '',
                                     'submitted_at': (opened + dt.timedelta(seconds=50)).isoformat(),
                                     'source': 'simulated'}) + '\n')
                n += 1
    return n


def selftest() -> int:
    # 1. the statistics against published or independent values
    wiki = [[0, 0, 0, 0, 14], [0, 2, 6, 4, 2], [0, 0, 3, 5, 6], [0, 3, 9, 2, 0], [2, 2, 8, 1, 1],
            [7, 7, 0, 0, 0], [3, 2, 6, 3, 0], [2, 5, 3, 2, 2], [6, 5, 2, 1, 0], [0, 2, 2, 3, 7]]
    cats5 = list('abcde')
    k = kappa([collections.Counter(dict(zip(cats5, r))) for r in wiki], cats5)
    assert abs(k - 0.2099) < 5e-4, k                          # Fleiss (1971) worked example: kappa 0.210
    rng = random.Random(7)
    pairs = [(rng.choice(LEVELS), rng.choice(LEVELS)) for _ in range(40)]
    for kind in ('identity', 'linear', 'quadratic'):          # two raters: the Gwet form is Scott's pi
        w = weights(kind)
        idx = {c: i for i, c in enumerate(LEVELS)}
        po = sum(w[idx[a]][idx[b]] for a, b in pairs) / len(pairs)
        pi = [sum((a == c) + (b == c) for a, b in pairs) / (2 * len(pairs)) for c in LEVELS]
        pe = sum(w[i][j] * pi[i] * pi[j] for i in range(3) for j in range(3))
        scott = (po - pe) / (1 - pe)
        got = kappa([collections.Counter([a, b]) for a, b in pairs], LEVELS, w)
        assert abs(got - scott) < 1e-12, (kind, got, scott)
    sys.path.insert(0, str(REPO / 'full_run_28092026'))
    from full_run_28092026.expert_kits import fleiss as pooled_fleiss
    three = [collections.Counter(rng.choice(LEVELS) for _ in range(3)) for _ in range(60)]
    assert abs(kappa(three, LEVELS) - pooled_fleiss(three)) < 1e-12   # equal rater counts: Fleiss' pooled form
    assert kappa([collections.Counter({'Easy': 3})] * 5, LEVELS) is None   # one category throughout: 0/0
    print("statistics: Fleiss' worked example 0.210 reproduced; two raters give Scott's pi under all three weightings; "
          "equal to the pooled form at three raters; undefined where every rating agrees on one level")
    # 2. the pipeline on simulated returns
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-levels-score-'))
    try:
        tasks = tmp / 'tasks'
        build_kit.build(tasks, tmp / 'dist', build_kit.SEED, None)
        n = simulate(tasks, tmp / 'perfect', 1.0)
        res = score(tmp / 'perfect', tasks, tmp / 'a.json', tmp / 'a.md')
        assert n == 450 and res['returns']['ratings'] == 450, n
        assert res['agreement']['overall']['fleiss'] == 1.0, res['agreement']['overall']
        assert res['majority']['overall']['same_as_current'] == 150, res['majority']['overall']
        n = simulate(tasks, tmp / 'noisy', 0.7, seed=3)
        res = score(tmp / 'noisy', tasks, tmp / 'b.json', tmp / 'b.md')
        o = res['agreement']['overall']
        assert 0 < o['fleiss'] < o['weighted_linear'] < o['weighted_quadratic'] < 1, o   # one-step noise: weights help
        assert o['fleiss_ci95'][0] < o['fleiss'] < o['fleiss_ci95'][1], o
        m = res['majority']['overall']
        assert m['same_as_current'] + m['differs_by_one_step'] + m['differs_by_two_steps'] == m['with_majority'], m
        assert sum(t['ratings'][lv] for t in res['templates'].values() for lv in LEVELS) == 450
        assert 'ele-' not in (tmp / 'b.json').read_text(encoding='utf-8'), 'an expert id reached agreement.json'
        res = score_partial(tasks, tmp)
        bad = tmp / 'bad'
        bad.mkdir()
        (bad / 'x.jsonl').write_text(json.dumps({'annotator_id': 'ele-1', 'code': 'LV-nope', 'level': 'Easy',
                                                 'formula_stated': 'no'}) + '\n', encoding='utf-8')
        try:
            score(bad, tasks, None, None)
            raise AssertionError('an unknown code was accepted')
        except SystemExit:
            pass
        print('pipeline: perfect returns give kappa 1 and no change; noisy returns give 0 < kappa < linear < quadratic < 1 '
              'inside its interval; one expert missing is scored on the rest; an unknown code stops the run')
        print('SELFTEST OK')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def score_partial(tasks: Path, tmp: Path) -> dict:
    simulate(tasks, tmp / 'partial', 0.8, seed=5, drop='che-2')
    res = score(tmp / 'partial', tasks, None, None)
    assert res['returns']['experts_returned'] == 14 and res['returns']['ratings'] == 420, res['returns']
    assert res['agreement']['by_branch']['chemical_engineering']['templates'] == 30   # two raters still count
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--returns', type=Path, default=RETURNS)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.returns.exists():
        raise SystemExit(f'no returns folder: {a.returns}')
    score(a.returns, build_kit.TASKS, AGREEMENT, build_kit.RESULTS)
    return 0


if __name__ == '__main__':
    sys.exit(main())
