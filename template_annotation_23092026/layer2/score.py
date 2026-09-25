"""Layer 2 - read the experts' labels back: agreement, plant detection, hand checks, and the fix list.

    python -m template_annotation_23092026.layer2.score --labels <round-1 folder>
    python -m template_annotation_23092026.layer2.score --round 2 --labels <round-2 folder> \
        --prev-labels <round-1 folder>

Reads the experts' <id>.jsonl files (app and workbook rows alike; default layer2/labels/),
tasks/keyfile.jsonl (which codes are plants) and screen/pass2/summary.csv (the panel's
verdict on the same corpus). Writes layer2/RESULTS.md. Every number the paper will quote
about human certification comes from here.

A later round (--round N) reads tasks_round<N>/keyfile.jsonl and writes RESULTS_round<N>.md,
so round 1's record is never overwritten. It has no plants (sensitivity was measured in
round 1) and no panel comparison (the panel judged the templates before the fixes), and
with --prev-labels it adds the round-over-round table: each template's verdicts in the
previous round beside this one, expert by expert.

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


# ------------------------------------------------------------------ input

def tasks_dir(n: int) -> Path:
    return TASKS if n == 1 else HERE / f'tasks_round{n}'


def read_key(n: int) -> dict:
    key = {}
    for ln in (tasks_dir(n) / 'keyfile.jsonl').open(encoding='utf8'):
        if ln.strip():
            k = json.loads(ln)
            key[k['code']] = k
    return key


def read_rows(folder: Path) -> list[dict]:
    rows = []
    for p in sorted(folder.glob('*.jsonl')):
        for ln in p.open(encoding='utf8'):
            if ln.strip():
                rows.append(json.loads(ln))
    latest = {}
    for r in rows:                                   # a resubmission replaces the earlier row
        latest[(r['annotator_id'], r['code'])] = r
    return list(latest.values())


# ------------------------------------------------------------------ main

def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=None, help='folder of <id>.jsonl files (default: layer2/labels/)')
    ap.add_argument('--round', type=int, default=1, help='read tasks_round<N>/, write RESULTS_round<N>.md')
    ap.add_argument('--prev-labels', default=None,
                    help="a later round: the previous round's label folder, for the round-over-round table")
    a, _ = ap.parse_known_args()
    global LABELS
    if a.labels:
        LABELS = Path(a.labels)
    key = read_key(a.round)
    rows = read_rows(LABELS)
    if not rows:
        raise SystemExit('no labels yet')
    unknown = sorted({r['code'] for r in rows if r['code'] not in key})
    if unknown:
        raise SystemExit(f'{len(unknown)} label codes are not in {tasks_dir(a.round).name}/keyfile.jsonl '
                         f'(wrong --round, or the tasks were rebuilt after the kits went out): {unknown[:5]}')
    by_code = collections.defaultdict(dict)
    for r in rows:
        by_code[r['code']][r['annotator_id']] = r
    experts = sorted({r['annotator_id'] for r in rows})
    branch_of = {r['annotator_id']: r['branch'] for r in rows}

    out = []
    w = out.append
    w('# Layer 2 - human certification results' + (f', round {a.round}' if a.round > 1 else '') + '\n')
    w(f'Generated {dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")} by `score.py` from '
      f'{len(rows)} label rows by {len(experts)} experts.\n')
    if a.round > 1:
        build = json.loads((tasks_dir(a.round) / 'BUILD.json').read_text(encoding='utf8'))
        w(f"Round {a.round} re-certifies the {build['templates']} templates changed after round {a.round - 1}, "
          f"built at git `{build['git_head'][:7]}`: no planted defects, and a fresh hand-check instance "
          f"(seed {build['instance_seeds'][0]}) that no expert saw before.\n")

    # ---- plants
    tot_seen = tot_hit = 0
    if any(k['kind'] == 'plant' for k in key.values()):
        w('## Plant detection (sensitivity)\n')
        w('| Expert | Branch | Plants seen | Rejected | Detection | Real templates approved |\n|---|---|---:|---:|---:|---:|')
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
    else:
        w('## Verdicts per expert\n')
        w("No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).\n")
        w('| Expert | Branch | Items | Approved | Rejected |\n|---|---|---:|---:|---:|')
        for e in experts:
            mine = [r for r in rows if r['annotator_id'] == e]
            appr = sum(1 for r in mine if r['decision'] == 'Approve')
            w(f'| {e} | {branch_of[e].split("_")[0]} | {len(mine)} | {appr} | {len(mine) - appr} |')

    # ---- round over round
    prev = collections.defaultdict(dict)                 # template_id -> {expert: previous-round row}
    if a.round > 1 and a.prev_labels:
        prev_key = read_key(a.round - 1)
        for r in read_rows(Path(a.prev_labels)):
            k = prev_key.get(r['code'])
            if k and k['kind'] == 'template':
                prev[k['template_id']][r['annotator_id']] = r
        w(f'\n## Round {a.round - 1} against round {a.round}\n')
        w("Each expert's verdict on the same template in both rounds: A approve, R reject.\n")
        w(f'| Template | Round {a.round - 1} | Round {a.round} | Outcome |\n|---|---|---|---|')
        trans, outcomes = collections.Counter(), collections.Counter()
        for code, seen in sorted(by_code.items(), key=lambda kv: key[kv[0]].get('template_id', '')):
            k = key[code]
            if k['kind'] != 'template':
                continue
            before = prev.get(k['template_id'], {})
            ids = sorted(seen)
            fmt = lambda d: ' '.join(f"{e} {d[e]['decision'][0] if e in d else '-'}" for e in ids)  # noqa: E731
            n_rej = sum(1 for r in seen.values() if r['decision'] == 'Reject')
            outcome = ('approved by all' if n_rej == 0 else
                       f'approved by majority, rejected by {n_rej}' if 2 * n_rej < len(seen) else
                       f'rejected by majority, {n_rej} of {len(seen)}')
            outcomes[outcome.split(',')[0]] += 1
            for e, r in seen.items():
                if e in before:
                    trans[(before[e]['decision'], r['decision'])] += 1
            w(f"| {k['template_id']} | {fmt(before)} | {fmt(seen)} | {outcome} |")
        n_prev_rej = trans[('Reject', 'Approve')] + trans[('Reject', 'Reject')]
        w(f"\nTemplates: {outcomes['approved by all']} approved by all three, "
          f"{outcomes['approved by majority']} approved by majority, {outcomes['rejected by majority']} rejected by majority.")
        w(f"Of the {n_prev_rej} round-{a.round - 1} rejections of these templates, {trans[('Reject', 'Approve')]} became "
          f"approvals and {trans[('Reject', 'Reject')]} stayed rejections; {trans[('Approve', 'Reject')]} verdicts went "
          f"the other way, from approve to reject.")

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
    if a.round > 1:
        w('Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.')
    elif SCREEN.exists():
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
            t_open = dt.datetime.fromisoformat(r['opened_at'])
            t_sub = dt.datetime.fromisoformat(r['submitted_at'])
            dw.append((r['annotator_id'], (t_sub - t_open).total_seconds() / 60))
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
                was = prev.get(k['template_id'], {}).get(e)
                tag = f" (round {a.round - 1}: {was['decision'].lower()})" if was else ''
                w(f"  - {e}{tag} [{', '.join(r['defects'])}]: {r['feedback']}")
    if not any_rej:
        w('none')

    out_path = HERE / ('RESULTS.md' if a.round == 1 else f'RESULTS_round{a.round}.md')
    out_path.write_text('\n'.join(out) + '\n', encoding='utf8')
    print(f'wrote {out_path}: {len(rows)} rows, {len(experts)} experts, '
          f'{tot_hit}/{tot_seen} plants detected')


if __name__ == '__main__':
    main()
