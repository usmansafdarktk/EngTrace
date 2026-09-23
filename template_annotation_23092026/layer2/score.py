"""Layer 2 - read the experts' labels back: agreement, plant detection, hand checks, and the fix list.

    python -m template_annotation_23092026.layer2.score

Reads labels/*.jsonl (app and workbook rows alike), tasks/keyfile.jsonl (which codes
are plants) and screen/pass2/summary.csv (the panel's verdict on the same corpus).
Writes layer2/RESULTS.md. Every number the paper will quote about human certification
comes from here.

What it reports, and why each exists:
  * plant detection per expert and overall - the sensitivity of the review, which a
    set of approvals alone cannot show (the December 2025 record had none);
  * hand-check agreement - how often an expert's own answer matched the template's,
    and where it did not, whether they rejected;
  * agreement among the three experts of a branch on Approve/Reject (Fleiss kappa,
    Gwet AC1, percent) and on the 1-5 scores (Gwet AC2 with quadratic weights,
    Krippendorff-style alpha is deliberately not repeated: see Appendix K's discussion);
  * the panel-vs-experts false-positive rate: templates the screen passed that the
    experts reject by majority, and MAD between expert median scores and the panel's;
  * dwell time from the app's timestamps (workbook rows carry none);
  * the fix list: every template rejected by at least one expert, with the notes.
"""
from __future__ import annotations

import collections
import csv
import datetime as dt
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
LABELS = HERE / 'labels'
SCREEN = HERE.parent / 'screen' / 'pass2' / 'summary.csv'
DIMS = ['physical_plausibility', 'mathematical_correctness', 'pedagogical_clarity']


# ------------------------------------------------------------- agreement

def _counts(ratings: list[list[int]], q: int):
    r, n = len(ratings[0]), len(ratings)
    pa, marg = 0.0, collections.Counter()
    for row in ratings:
        c = collections.Counter(row)
        pa += sum(v * (v - 1) for v in c.values()) / (r * (r - 1))
        marg.update(c)
    return pa / n, [marg[k] / (n * r) for k in range(q)]


def fleiss_kappa(ratings, q=2):
    pa, pis = _counts(ratings, q)
    pe = sum(p * p for p in pis)
    return (pa - pe) / (1 - pe) if pe < 1 else 1.0


def gwet_ac1(ratings, q=2):
    pa, pis = _counts(ratings, q)
    pe = sum(p * (1 - p) for p in pis) / (q - 1)
    return (pa - pe) / (1 - pe) if pe < 1 else 1.0


def gwet_ac2(ratings: list[list[int]], q: int = 5) -> float:
    """Gwet's AC2 with quadratic weights, categories 0..q-1 (Gwet 2014, ch. 4)."""
    w = [[1 - ((k - l) ** 2) / ((q - 1) ** 2) for l in range(q)] for k in range(q)]
    r, n = len(ratings[0]), len(ratings)
    pa = 0.0
    marg = collections.Counter()
    for row in ratings:
        c = collections.Counter(row)
        marg.update(c)
        for k in range(q):
            rk = c.get(k, 0)
            rstar = sum(w[k][l] * c.get(l, 0) for l in range(q))
            pa += rk * (rstar - 1) / (r * (r - 1))
    pa /= n
    pis = [marg[k] / (n * r) for k in range(q)]
    tw = sum(w[k][l] for k in range(q) for l in range(q) if k != l)
    pe = tw / (q * (q - 1)) * sum(p * (1 - p) for p in pis)
    return (pa - pe) / (1 - pe) if pe < 1 else 1.0


def pct(ratings):
    return sum(1 for row in ratings if len(set(row)) == 1) / len(ratings)


# ------------------------------------------------------------------ main

def main() -> None:
    key = {}
    for ln in (TASKS / 'keyfile.jsonl').open(encoding='utf8'):
        if ln.strip():
            k = json.loads(ln)
            key[k['code']] = k
    rows = []
    for p in sorted(LABELS.glob('*.jsonl')):
        for ln in p.open(encoding='utf8'):
            if ln.strip():
                rows.append(json.loads(ln))
    if not rows:
        raise SystemExit('no labels yet')
    latest = {}
    for r in rows:                                   # a resubmission replaces the earlier row
        latest[(r['annotator_id'], r['code'])] = r
    rows = list(latest.values())
    by_code = collections.defaultdict(dict)
    for r in rows:
        by_code[r['code']][r['annotator_id']] = r
    experts = sorted({r['annotator_id'] for r in rows})
    branch_of = {r['annotator_id']: r['branch'] for r in rows}

    out = []
    w = out.append
    w('# Layer 2 - human certification results\n')
    w(f'Generated {dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")} by `score.py` from '
      f'{len(rows)} label rows by {len(experts)} experts.\n')

    # ---- plants
    w('## Plant detection (sensitivity)\n')
    w('| Expert | Branch | Plants seen | Rejected | Detection | Real templates approved |\n|---|---|---:|---:|---:|---:|')
    tot_seen = tot_hit = 0
    for e in experts:
        mine = [r for r in rows if r['annotator_id'] == e]
        plants = [r for r in mine if key[r['code']]['kind'] == 'plant']
        hit = sum(1 for r in plants if r['decision'] == 'Reject')
        real = [r for r in mine if key[r['code']]['kind'] == 'template']
        appr = sum(1 for r in real if r['decision'] == 'Approve')
        tot_seen += len(plants)
        tot_hit += hit
        w(f'| {e} | {branch_of[e].split("_")[0]} | {len(plants)} | {hit} | '
          f'{(hit / len(plants)) if plants else float("nan"):.0%} | {appr} of {len(real)} |')
    w(f'\nOverall: {tot_hit} of {tot_seen} planted defects rejected ({(tot_hit / tot_seen) if tot_seen else 0:.0%}).\n')
    w('| Plant | Defect | Experts who saw it | Rejected by |\n|---|---|---:|---|')
    for code, k in key.items():
        if k['kind'] != 'plant' or code not in by_code:
            continue
        seen = by_code[code]
        rej = [e for e, r in seen.items() if r['decision'] == 'Reject']
        w(f"| {k['plant_id']} ({k['base']}) | {k['defect_class']}: {k['description']} | {len(seen)} | {', '.join(rej) or 'nobody'} |")

    # ---- hand checks
    w('\n## Hand checks\n')
    hc = [r for r in rows if r.get('hand_match') is not None]
    if hc:
        m = sum(1 for r in hc if r['hand_match'])
        w(f'{len(hc)} hand checks with a comparable number: {m} matched the template within 1% ({m / len(hc):.0%}).')
        mism = [r for r in hc if not r['hand_match']]
        rej = sum(1 for r in mism if r['decision'] == 'Reject')
        w(f'Of the {len(mism)} mismatches, {rej} ended in a rejection and {len(mism) - rej} in an approval '
          f'(the expert found their own slip, or judged the difference immaterial).')
        real_m = [r for r in hc if key[r['code']]['kind'] == 'template' and not r['hand_match'] and r['decision'] == 'Approve']
        if real_m:
            w('\nReal templates approved despite a mismatch (check these):')
            for r in real_m:
                w(f"- {key[r['code']]['template_id']} by {r['annotator_id']}: entered `{r['hand_answer'][:60]}`; note: {r['feedback'][:120] or '-'}")
    else:
        w('No app rows with a comparable hand check yet.')

    # ---- agreement, real templates only, three experts
    w('\n## Agreement among the three experts of a branch (real templates)\n')
    per_branch = collections.defaultdict(list)
    for code, seen in by_code.items():
        if key[code]['kind'] == 'template' and len(seen) == 3:
            per_branch[key[code]['branch']].append(seen)
    w('| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | '
      'AC2 phys | AC2 math | AC2 ped |\n|---|---:|---:|---:|---:|---:|---:|---:|')
    all_bin, all_dim = [], {d: [] for d in DIMS}
    for b in sorted(per_branch):
        seen_list = per_branch[b]
        bin_ = [[int(r['decision'] == 'Reject') for r in seen.values()] for seen in seen_list]
        dims = {d: [[r['scores'][d] - 1 for r in seen.values()] for seen in seen_list] for d in DIMS}
        all_bin += bin_
        for d in DIMS:
            all_dim[d] += dims[d]
        w(f'| {b.split("_")[0]} | {len(seen_list)} | {fleiss_kappa(bin_):.3f} | {gwet_ac1(bin_):.3f} | {pct(bin_):.0%} | '
          + ' | '.join(f'{gwet_ac2(dims[d]):.3f}' for d in DIMS) + ' |')
    if all_bin:
        w(f'| all | {len(all_bin)} | {fleiss_kappa(all_bin):.3f} | {gwet_ac1(all_bin):.3f} | {pct(all_bin):.0%} | '
          + ' | '.join(f'{gwet_ac2(all_dim[d]):.3f}' for d in DIMS) + ' |')
    w('\nKappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); '
      'AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.')

    # ---- against the screen
    w('\n## Against the screening panel (pass 2)\n')
    if SCREEN.exists():
        screen = {r['template_id']: r for r in csv.DictReader(SCREEN.open(encoding='utf8'))}
        fp, n_pass, mad = 0, 0, {d: [] for d in DIMS}
        smap = {'physical_plausibility': 'med_phys', 'mathematical_correctness': 'med_math', 'pedagogical_clarity': 'med_ped'}
        rejected_by_experts = []
        for code, seen in by_code.items():
            k = key[code]
            if k['kind'] != 'template' or len(seen) < 3:
                continue
            s = screen.get(k['template_id'])
            if not s:
                continue
            maj_reject = sum(1 for r in seen.values() if r['decision'] == 'Reject') >= 2
            if s['category'] == 'Pass':
                n_pass += 1
                if maj_reject:
                    fp += 1
                    rejected_by_experts.append(k['template_id'])
            for d in DIMS:
                med = statistics.median(r['scores'][d] for r in seen.values())
                mad[d].append(abs(med - float(s[smap[d]])))
        w(f'Templates the panel passed and the experts reject by majority: {fp} of {n_pass} '
          f'(false-positive rate {(fp / n_pass) if n_pass else 0:.1%}).'
          + (f' They are: {", ".join(rejected_by_experts)}.' if rejected_by_experts else ''))
        w('\n| Dimension | MAD, expert median vs panel median |\n|---|---:|')
        for d in DIMS:
            if mad[d]:
                w(f'| {d} | {statistics.mean(mad[d]):.3f} |')
    else:
        w('No screen summary found.')

    # ---- dwell
    w('\n## Time spent (app rows only)\n')
    dw = []
    for r in rows:
        if r.get('opened_at') and r.get('submitted_at'):
            a = dt.datetime.fromisoformat(r['opened_at'])
            b = dt.datetime.fromisoformat(r['submitted_at'])
            dw.append((r['annotator_id'], (b - a).total_seconds() / 60))
    if dw:
        w('| Expert | Items | Median minutes | Under 2 min |\n|---|---:|---:|---:|')
        for e in experts:
            mine = [m for a, m in dw if a == e]
            if mine:
                w(f'| {e} | {len(mine)} | {statistics.median(mine):.1f} | {sum(1 for m in mine if m < 2)} |')
    else:
        w('No timestamps (workbook route only).')

    # ---- fix list
    w('\n## Fix list: real templates rejected by at least one expert\n')
    any_rej = False
    for code, seen in sorted(by_code.items(), key=lambda kv: key[kv[0]].get('template_id', '')):
        k = key[code]
        if k['kind'] != 'template':
            continue
        rej = {e: r for e, r in seen.items() if r['decision'] == 'Reject'}
        if rej:
            any_rej = True
            w(f"- **{k['template_id']}** rejected by {len(rej)} of {len(seen)}:")
            for e, r in rej.items():
                w(f"  - {e} [{', '.join(r['defects'])}]: {r['feedback']}")
    if not any_rej:
        w('none')

    (HERE / 'RESULTS.md').write_text('\n'.join(out) + '\n', encoding='utf8')
    print(f'wrote {HERE / "RESULTS.md"}: {len(rows)} rows, {len(experts)} experts, '
          f'{tot_hit}/{tot_seen} plants detected')


if __name__ == '__main__':
    main()
