"""C2, the judge-swap robustness check (docs/EVALUATION_NEXT_STEPS.md C2; D-181): a second judge from another family on
a per-model sample of the traces E5 sends, against MiMo-V2.5-Pro's verdicts on the same milestones.

    python -m full_run_28092026.judge_swap [--judge x-ai/grok-4.6] [--variant main]   # FREE: JUDGE_SWAP.md, results/judge_swap.json

Reads scores/<variant>/e5/ (MiMo, the stage the results use) and scores/<variant>/e5_<slug>/ (the second judge, written
by `judge.py --judge ... --sample ...` on the same prompts), and for the traces both hold: per model, the milestones the
two judges both ruled on, their agreement on the three verdicts and on REACHED against not, E5-strict coverage over the
sampled traces under each judge with the paired difference and a template bootstrap, and the share of judged milestones
each judge rules REACHED. The question it answers is whether a model's E5-strict coverage figure (Q3) would move by more
than its interval under a judge from another family, on a sample of what E5 sends (20 traces drawn per model); it does
not check the 55 pairwise coverage differences, and it does not say which judge is right, which the pilot's expert
labels do for MiMo (E5_VALIDATION.md, RESULTS_E5).

A trace is compared only when both judges were sent the same prompt: the prompt lists the milestones E3 did not
reach, so the two rows must leave the same milestones to the judge. A trace whose E3 matches changed after the second
judge ruled (the milestone reader's fix at ANSWER FINAL) carries a new MiMo prompt and is left out and counted, as are
the drawn traces whose rows are no longer in the second judge's folder (the items re-run in round 5).
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER, boot_mean  # noqa: E402

B = 10_000


def rows_of(path: Path) -> dict:
    return {r['item_id']: r for r in map(json.loads, path.read_text(encoding='utf-8').splitlines())} if path.exists() else {}


def kappa(pairs, cls=lambda v: v):
    """Cohen's kappa over a Counter of (base verdict, other verdict) -> count, with `cls` mapping a verdict to the class
    compared (identity for three-way; REACHED against not for the binary reading). None when chance agreement is 1."""
    n = sum(pairs.values())
    if not n:
        return None
    po = sum(c for (a, b), c in pairs.items() if cls(a) == cls(b)) / n
    ma, mb = collections.Counter(), collections.Counter()
    for (a, b), c in pairs.items():
        ma[cls(a)] += c
        mb[cls(b)] += c
    pe = sum(ma[k] * mb[k] for k in set(ma) | set(mb)) / n ** 2
    return (po - pe) / (1 - pe) if pe < 1 else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--judge', default='x-ai/grok-4.6')
    ap.add_argument('--variant', default='main')
    a = ap.parse_args()
    slug = 'e5_' + a.judge.split('/')[-1].replace('.', '-')
    base_dir, other_dir = score.SCORES / a.variant / 'e5', score.SCORES / a.variant / slug
    if not other_dir.exists():
        raise SystemExit(f'{other_dir} does not exist: run judge.py --judge {a.judge} --sample N first')
    cfg_other = json.loads((other_dir / 'CONFIG.json').read_text(encoding='utf-8'))
    out = {'judge': a.judge, 'base_judge': json.loads((base_dir / 'CONFIG.json').read_text(encoding='utf-8')).get('judge'),
           'variant': a.variant, 'sample': len(cfg_other.get('sample') or []), 'models': {}}
    same_prompt = lambda b, o: [s == 'e3' for s in b['sources']] == [s == 'e3' for s in o['sources']]  # noqa: E731
    for key in ROSTER:
        base, other = rows_of(base_dir / f'{key}.jsonl'), rows_of(other_dir / f'{key}.jsonl')
        drawn = [i for i in other if i in base and other[i].get('sent')]
        common = [i for i in drawn if base[i].get('sent') and same_prompt(base[i], other[i])]
        prompt_changed = [i for i in drawn if i not in set(common)]       # a new prompt, or none (E3 now reaches all)
        if not common:
            continue
        agree3 = agree2 = n_ms = 0
        pair = collections.Counter()
        reached = collections.Counter()
        by_t = collections.defaultdict(list)
        no_reply = 0
        for i in common:
            b, o = base[i], other[i]
            if b.get('reply_ok') is False or o.get('reply_ok') is False:
                no_reply += 1
                continue
            for sb, so in zip(b['sources'], o['sources']):
                if sb == 'e3' or so == 'e3':
                    continue
                if 'UNJUDGED' in (sb, so):
                    continue
                n_ms += 1
                agree3 += sb == so
                agree2 += (sb == 'REACHED') == (so == 'REACHED')
                pair[(sb, so)] += 1
                reached['base'] += sb == 'REACHED'
                reached['other'] += so == 'REACHED'
            if b.get('e5_strict') is not None and o.get('e5_strict') is not None:
                by_t[b['template_id']].append(o['e5_strict'] - b['e5_strict'])
        d_all = np.array([x for v in by_t.values() for x in v])
        t_means = np.array([np.mean(v) for v in by_t.values()])
        out['models'][key] = {
            'traces': len(common), 'prompt_changed': len(prompt_changed), 'without_reply': no_reply, 'milestones_both_judged': n_ms,
            'agreement_three_way': agree3 / n_ms if n_ms else None, 'agreement_reached_vs_not': agree2 / n_ms if n_ms else None,
            'kappa_three_way': kappa(pair), 'kappa_reached_vs_not': kappa(pair, lambda v: v == 'REACHED'),
            'reached_share_base': reached['base'] / n_ms if n_ms else None, 'reached_share_other': reached['other'] / n_ms if n_ms else None,
            'pairs': {f'{sb}->{so}': n for (sb, so), n in sorted(pair.items())},
            'e5_strict_base': float(np.mean([base[i]['e5_strict'] for i in common if base[i].get('e5_strict') is not None])),
            'e5_strict_other': float(np.mean([other[i]['e5_strict'] for i in common if other[i].get('e5_strict') is not None])),
            'diff': float(d_all.mean()) if len(d_all) else None,
            'diff_ci_templates': boot_mean(t_means, 777) if len(t_means) > 1 else None,
            'templates': len(by_t)}
    OUT_JSON = HERE / 'results' / f'judge_swap_{a.variant}.json'
    ms_ = out['models'].values()
    out['traces'] = sum(m['traces'] for m in ms_)
    out['prompt_changed'] = sum(m['prompt_changed'] for m in ms_)
    out['not_in_rows'] = out['sample'] - out['traces'] - out['prompt_changed']
    out['traces_per_model'] = [min(m['traces'] for m in ms_), max(m['traces'] for m in ms_)]
    out['templates_per_model'] = [min(m['templates'] for m in ms_), max(m['templates'] for m in ms_)]
    OUT_JSON.write_text(json.dumps(out, indent=1), encoding='utf-8')
    span = lambda lo_hi: f'{lo_hi[0]}' if lo_hi[0] == lo_hi[1] else f'{lo_hi[0]} to {lo_hi[1]}'  # noqa: E731
    L = ['# C2: the judge swap, a second judge on a sample of what E5 sends', '',
         f"Generated by `judge_swap.py`; what it compares is in its docstring. `{out['judge']}` against `{out['base_judge']}` on "
         f"{out['traces']} sampled traces of the `{a.variant}` store, each sent to both judges with the same prompt (D-181); of the "
         f"{out['sample']} traces drawn, {out['not_in_rows']} are no longer in the second judge's rows (items re-run since) and "
         f"{out['prompt_changed']} carry a prompt that changed after the second judge ruled, and both are left out. Agreement is "
         'over the milestones both judges ruled on; coverage is E5-strict over the sampled traces under each judge, and the '
         'difference is the other judge minus MiMo with a template bootstrap. Which judge is right is not measured here; MiMo\'s '
         'validation against the experts is.', '',
         '| model | traces | milestones both judged | agreement, three-way | agreement, REACHED vs not | kappa, three-way | '
         'kappa, REACHED vs not | REACHED share: MiMo / other | '
         'E5-strict: MiMo / other | difference | 95% CI, templates |', '|---|---:|---:|---:|---:|---:|---:|---|---|---:|---:|']
    fk = lambda v: '-' if v is None else f'{v:.3f}'
    for k, m in out['models'].items():
        ci = m['diff_ci_templates']
        L.append(f"| `{k}` | {m['traces']} | {m['milestones_both_judged']} | {m['agreement_three_way']:.3f} | {m['agreement_reached_vs_not']:.3f} | "
                 f"{fk(m['kappa_three_way'])} | {fk(m['kappa_reached_vs_not'])} | "
                 f"{m['reached_share_base']:.3f} / {m['reached_share_other']:.3f} | {m['e5_strict_base']:.3f} / {m['e5_strict_other']:.3f} | "
                 f"{m['diff']:+.3f} | " + (f"{ci[0]:+.3f} to {ci[1]:+.3f}" if ci else '-') + ' |')
    pairs = collections.Counter()
    for m in out['models'].values():
        for k, n in m['pairs'].items():
            pairs[k] += n
    pooled = collections.Counter({tuple(k.split('->')): n for k, n in pairs.items()})
    out['pooled'] = {'milestones': sum(pooled.values()),
                     'agreement_three_way': sum(c for (a, b), c in pooled.items() if a == b) / max(1, sum(pooled.values())),
                     'agreement_reached_vs_not': sum(c for (a, b), c in pooled.items() if (a == 'REACHED') == (b == 'REACHED')) / max(1, sum(pooled.values())),
                     'kappa_three_way': kappa(pooled), 'kappa_reached_vs_not': kappa(pooled, lambda v: v == 'REACHED')}
    OUT_JSON.write_text(json.dumps(out, indent=1), encoding='utf-8')
    pk = out['pooled']
    L += ['', 'Verdict pairs over all sampled milestones (MiMo -> other): ' + ', '.join(f'{k} {n}' for k, n in pairs.most_common()) + '.', '',
          f"Pooled over the {pk['milestones']} milestones: agreement {pk['agreement_three_way']:.3f} three-way and "
          f"{pk['agreement_reached_vs_not']:.3f} on REACHED against not; Cohen's kappa {fk(pk['kappa_three_way'])} and "
          f"{fk(pk['kappa_reached_vs_not'])}. Kappa is chance-corrected and reads lower than the raw agreement where one verdict "
          'dominates: on what E5 sends (milestones E3 did not find), MISSING and NOT_NEEDED are the common verdicts and REACHED the '
          "minority (MiMo's REACHED share 0.09 to 0.36 per model above). The per-model coverage difference with its interval is "
          'the figure that answers the question the swap asks; kappa is given so that a reader can see the prevalence effect '
          f"rather than suspect it. The sample is what E5 sends, {span(out['traces_per_model'])} traces per model over "
          f"{span(out['templates_per_model'])} templates: it checks the per-model coverage figure, not the "
          f"{len(out['models']) * (len(out['models']) - 1) // 2} pairwise coverage differences of Q3.", '']
    (HERE / 'JUDGE_SWAP.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    sys.exit(main())
