"""Consensus and agreement for one screening pass.

    python -m template_annotation_23092026.screen.analyze_screen --pass 1

Reads screen/pass<N>/replies.jsonl and writes, beside it:

  summary.csv   one row per template: category, medians, sigma_max, votes, each judge's scores
  flagged.md    the templates the panel flags, with every judge's explanation
  stats.md      the numbers the paper needs, all computed here

The consensus rule is the paper's (section 3.3, equations 1-2), applied unchanged:
a template passes when the median of the three judges reaches 4 on every scored
dimension AND fewer than two judges raise the review flag; among the survivors,
sigma_max (the largest population standard deviation across the three dimensions)
above 0.5 marks the template controversial and routes it to manual review. The
categories are the December run's: Critical Failure, Controversial, Pass.

Agreement. The July 2026 rebuttal promised the aggregate sigma_max and pairwise
inter-judge agreement (Gwet's AC1) on the binary certification decision. Both are
here, with Fleiss' kappa and raw percent agreement beside them so a reader can see
what the chance correction does. AC1 follows Gwet (2008): for q categories and r
raters, p_a = mean over items of sum_k r_ik (r_ik - 1) / (r (r - 1)), p_e =
sum_k pi_k (1 - pi_k) / (q - 1), AC1 = (p_a - p_e) / (1 - p_e). Self-tested below.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCORE_DIMS = ['physical_plausibility_score', 'mathematical_correctness_score', 'pedagogical_clarity_score']
SHORT = {'physical_plausibility_score': 'phys', 'mathematical_correctness_score': 'math',
         'pedagogical_clarity_score': 'ped'}
THETA = 4          # median threshold, every dimension (eq. 1)
SIGMA_MAX = 0.5    # disagreement threshold (eq. 2)


# ------------------------------------------------------------------ agreement

def gwet_ac1(ratings: list[list[int]], q: int = 2) -> float:
    """ratings[i] = the r raters' category labels for item i (all items fully rated)."""
    r = len(ratings[0])
    n = len(ratings)
    pa = 0.0
    counts = collections.Counter()
    for row in ratings:
        c = collections.Counter(row)
        pa += sum(v * (v - 1) for v in c.values()) / (r * (r - 1))
        counts.update(c)
    pa /= n
    pis = [counts[k] / (n * r) for k in range(q)]
    pe = sum(p * (1 - p) for p in pis) / (q - 1)
    return (pa - pe) / (1 - pe) if pe < 1 else 1.0


def fleiss_kappa(ratings: list[list[int]], q: int = 2) -> float:
    r = len(ratings[0])
    n = len(ratings)
    pa, counts = 0.0, collections.Counter()
    for row in ratings:
        c = collections.Counter(row)
        pa += sum(v * (v - 1) for v in c.values()) / (r * (r - 1))
        counts.update(c)
    pa /= n
    pe = sum((counts[k] / (n * r)) ** 2 for k in range(q))
    return (pa - pe) / (1 - pe) if pe < 1 else 1.0


def percent_agreement(ratings: list[list[int]]) -> float:
    return sum(1 for row in ratings if len(set(row)) == 1) / len(ratings)


def _selftest() -> None:
    same = [[0, 0]] * 10
    assert abs(gwet_ac1(same) - 1.0) < 1e-12
    opposite = [[0, 1]] * 10
    assert abs(gwet_ac1(opposite) - (-1.0)) < 1e-12
    # Gwet (2008) worked example shape: high agreement, skewed prevalence.
    skew = [[0, 0]] * 9 + [[0, 1]]
    ac1, kap = gwet_ac1(skew), fleiss_kappa(skew)
    assert ac1 > kap, (ac1, kap)      # AC1 does not collapse under prevalence
    assert abs(percent_agreement(skew) - 0.9) < 1e-12


_selftest()


# ------------------------------------------------------------------ consensus

def pstdev(xs: list[int]) -> float:
    return statistics.pstdev(xs) if len(xs) > 1 else 0.0


def consensus(rows: dict[str, dict]) -> dict:
    """rows: judge key -> parsed scores. Returns the paper's certification decision."""
    judges = sorted(rows)
    med = {d: statistics.median([rows[j][d] for j in judges]) for d in SCORE_DIMS}
    sig = {d: pstdev([rows[j][d] for j in judges]) for d in SCORE_DIMS}
    flags = sum(1 for j in judges if rows[j]['human_review_flag'])
    sigma_max = max(sig.values())
    if flags >= 2 or any(m < THETA for m in med.values()):
        cat = 'Critical Failure'
    elif sigma_max > SIGMA_MAX:
        cat = 'Controversial'
    else:
        cat = 'Pass'
    return {'category': cat, 'medians': med, 'sigmas': sig, 'sigma_max': sigma_max,
            'flags': flags, 'n_judges': len(judges)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--pass', dest='n', type=int, default=1)
    a = ap.parse_args()
    d = HERE / f'pass{a.n}'
    path = d / 'replies.jsonl'
    if not path.exists():
        raise SystemExit(f'no replies at {path}')
    cfg = json.loads((d / 'config.json').read_text(encoding='utf8'))
    judges = [j['key'] for j in cfg['panel']]
    rows = [json.loads(ln) for ln in path.open(encoding='utf8') if ln.strip()]
    by = collections.defaultdict(dict)
    meta = {}
    for r in rows:
        meta[r['template_id']] = (r['branch'], r['file'])
        if r['ok']:
            by[r['template_id']][r['judge']] = r
    complete = {t: v for t, v in by.items() if len(v) == len(judges)}
    incomplete = sorted(t for t in meta if t not in complete)

    # ---- per template
    out = []
    for t in sorted(complete):
        scores = {j: complete[t][j]['scores'] for j in judges}
        c = consensus(scores)
        out.append({'template_id': t, 'branch': meta[t][0], **c, 'scores': scores,
                    'explanations': {j: scores[j]['explanation'] for j in judges}})
    with (d / 'summary.csv').open('w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(['template_id', 'branch', 'category', 'med_phys', 'med_math', 'med_ped', 'sigma_max',
                    'flags'] + [f'{j}_{SHORT[dim]}' for j in judges for dim in SCORE_DIMS]
                   + [f'{j}_flag' for j in judges] + [f'{j}_conf' for j in judges])
        for o in out:
            w.writerow([o['template_id'], o['branch'], o['category']]
                       + [o['medians'][dim] for dim in SCORE_DIMS] + [round(o['sigma_max'], 3), o['flags']]
                       + [o['scores'][j][dim] for j in judges for dim in SCORE_DIMS]
                       + [int(o['scores'][j]['human_review_flag']) for j in judges]
                       + [o['scores'][j]['confidence_score'] for j in judges])

    flagged = [o for o in out if o['category'] != 'Pass']
    with (d / 'flagged.md').open('w', encoding='utf8') as fh:
        fh.write(f'# Pass {a.n} - templates the panel flags\n\n')
        fh.write(f'{len(flagged)} of {len(out)} templates; consensus rule as in the paper, section 3.3. '
                 f'Each judge\'s one-sentence explanation is quoted verbatim.\n\n')
        for o in sorted(flagged, key=lambda o: (o['category'], o['branch'], o['template_id'])):
            fh.write(f"## {o['template_id']}  [{o['branch']}]  {o['category']}\n\n")
            fh.write(f"medians phys/math/ped {o['medians'][SCORE_DIMS[0]]:g}/{o['medians'][SCORE_DIMS[1]]:g}/"
                     f"{o['medians'][SCORE_DIMS[2]]:g}, sigma_max {o['sigma_max']:.2f}, "
                     f"review flags {o['flags']}/{o['n_judges']}\n\n")
            for j in judges:
                s = o['scores'][j]
                fh.write(f"- **{j}** ({s[SCORE_DIMS[0]]}/{s[SCORE_DIMS[1]]}/{s[SCORE_DIMS[2]]}, "
                         f"flag={'yes' if s['human_review_flag'] else 'no'}, conf {s['confidence_score']}): "
                         f"{s['explanation']}\n")
            fh.write('\n')

    # ---- statistics
    cats = collections.Counter(o['category'] for o in out)
    by_branch = collections.defaultdict(collections.Counter)
    for o in out:
        by_branch[o['branch']][o['category']] += 1
    sig = [o['sigma_max'] for o in out]
    flag_ratings = [[int(o['scores'][j]['human_review_flag']) for j in judges] for o in out]
    pair_ac1 = {}
    for i in range(len(judges)):
        for k in range(i + 1, len(judges)):
            pr = [[row[i], row[k]] for row in flag_ratings]
            pair_ac1[(judges[i], judges[k])] = (gwet_ac1(pr), fleiss_kappa(pr), percent_agreement(pr))
    dim_exact = {dim: sum(1 for o in out if len({o['scores'][j][dim] for j in judges}) == 1) / len(out)
                 for dim in SCORE_DIMS}
    dim_ac1 = {dim: gwet_ac1([[o['scores'][j][dim] - 1 for j in judges] for o in out], q=5)
               for dim in SCORE_DIMS}
    cost = sum(r.get('cost_usd') or 0 for r in rows)
    served = collections.Counter((r['judge'], r.get('served_model')) for r in rows if r['ok'])
    attempts = collections.Counter((r['judge'], len(r.get('attempts') or [])) for r in rows)
    ok_rows = [r for r in rows if r['ok']]

    with (d / 'stats.md').open('w', encoding='utf8') as fh:
        w = fh.write
        w(f'# Pass {a.n} - the numbers\n\n')
        w(f'Computed by `analyze_screen.py` from `replies.jsonl` ({len(rows)} rows, {len(ok_rows)} ok). '
          f'Judges: {", ".join(judges)}. Corpus at git `{cfg["git_head"][:10]}`, prompt template '
          f'`{cfg["prompt_template_sha256"][:10]}`, instance seeds {cfg["seeds"]}.\n\n')
        w('## Certification outcome (paper rule, section 3.3)\n\n')
        w('| Category | Templates |\n|---|---:|\n')
        for c in ('Pass', 'Controversial', 'Critical Failure'):
            w(f'| {c} | {cats.get(c, 0)} |\n')
        w(f'| judged by all {len(judges)} | {len(out)} |\n')
        w(f'| incomplete (a judge failed) | {len(incomplete)} |\n\n')
        w('| Branch | Pass | Controversial | Critical Failure |\n|---|---:|---:|---:|\n')
        for b in sorted(by_branch):
            c = by_branch[b]
            w(f"| {b} | {c.get('Pass', 0)} | {c.get('Controversial', 0)} | {c.get('Critical Failure', 0)} |\n")
        w('\n## Disagreement among the judges (sigma_max, equation 2)\n\n')
        w('| Statistic | Value |\n|---|---:|\n')
        w(f'| mean sigma_max | {statistics.mean(sig):.3f} |\n')
        w(f'| median sigma_max | {statistics.median(sig):.3f} |\n')
        w(f'| max sigma_max | {max(sig):.3f} |\n')
        w(f'| templates with sigma_max <= {SIGMA_MAX} | {sum(1 for s in sig if s <= SIGMA_MAX)} of {len(sig)} |\n\n')
        w('## Inter-judge agreement on the review flag (binary)\n\n')
        allr = flag_ratings
        w(f'All {len(judges)} judges: Gwet AC1 {gwet_ac1(allr):.3f}, Fleiss kappa {fleiss_kappa(allr):.3f}, '
          f'percent agreement {percent_agreement(allr):.1%}.\n\n')
        w('| Pair | Gwet AC1 | Cohen/Fleiss kappa | Percent agreement |\n|---|---:|---:|---:|\n')
        for (x, y), (ac1, kap, pa) in pair_ac1.items():
            w(f'| {x} vs {y} | {ac1:.3f} | {kap:.3f} | {pa:.1%} |\n')
        w('\nPer judge: flag rate and mean scores.\n\n')
        w('| Judge | Flag rate | mean phys | mean math | mean ped | mean confidence |\n|---|---:|---:|---:|---:|---:|\n')
        for j in judges:
            js = [o['scores'][j] for o in out]
            w(f"| {j} | {sum(s['human_review_flag'] for s in js) / len(js):.1%} | "
              + ' | '.join(f"{statistics.mean(s[dim] for s in js):.2f}" for dim in SCORE_DIMS)
              + f" | {statistics.mean(s['confidence_score'] for s in js):.2f} |\n")
        w('\n## Agreement on the 1-5 scores\n\n')
        w('| Dimension | All three identical | Gwet AC1 (5 categories, unweighted) |\n|---|---:|---:|\n')
        for dim in SCORE_DIMS:
            w(f'| {SHORT[dim]} | {dim_exact[dim]:.1%} | {dim_ac1[dim]:.3f} |\n')
        w('\n## The run itself\n\n')
        carried = [r for r in rows if r.get('carried_from')]
        w(f'Cost: ${cost:.3f} for {len(rows) - len(carried)} judged rows '
          f'({sum(1 for r in rows if r.get("cost_source") == "openrouter" and not r.get("carried_from"))} priced by '
          f'OpenRouter\'s own usage.cost, the rest from the catalogue).'
          + (f' {len(carried)} rows were carried forward from pass {carried[0]["carried_from"]} because the '
             f'template\'s prompt text is byte-identical there; their original cost was '
             f'${sum(r.get("carried_cost_usd") or 0 for r in carried):.3f}.' if carried else '')
          + '\n\n')
        w('| Judge | Served model id | Rows |\n|---|---|---:|\n')
        for (j, m), c in sorted(served.items()):
            w(f'| {j} | {m} | {c} |\n')
        w('\n| Judge | Attempts needed | Rows |\n|---|---:|---:|\n')
        for (j, k), c in sorted(attempts.items()):
            w(f'| {j} | {k} | {c} |\n')
        if incomplete:
            w('\nIncomplete templates: ' + ', '.join(incomplete) + '\n')
    print(f'{len(out)} templates: {dict(cats)}; {len(incomplete)} incomplete; ${cost:.3f}')
    print(f'wrote {d / "summary.csv"}, {d / "flagged.md"}, {d / "stats.md"}')


if __name__ == '__main__':
    main()
