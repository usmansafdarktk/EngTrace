"""The numbers sheet for the two writing sessions (WS-G, step G5): every number the writers need, each with the file
and key that hold it, grouped by section of the paper. The headline is matched settings (results/matched.json and the
`matched` families of the other result files); the three re-run models' default rows follow each table as the secondary
comparison; the error analysis, the further conditions and the serving endpoints stay at the providers' defaults.

    python -m full_run_28092026.numbers_sheet      # FREE: writes docs/mock_review_workstreams/reports/NUMBERS_SHEET.md

Reads only files that committed scripts write (named beside each value): results/*.json (analyze.py, coverage_variants.py,
clause_variants.py, flag_precision.py, depth_model.py, judge_swap.py, decoding_table.py, rescore_diff.py,
run_traces.py --matched-config), EXPERT_REQUEST.md (expert_kits.py), PARAPHRASE_REVIEW.md (paraphrase_kit.py),
SCORER_VALIDATION.md (validate_scorer.py), symbolic/SYMBOLIC_CHECK.md (symbolic/validate.py), the certification
records under template_annotation_23092026/layer2/ (score.py, certification.py), manifest.jsonl (freeze.py), the
pilot's RESULTS_X1.md and PILOT_SUMMARY.md, and answer.py (the enabled templates, the own-digit cap). It computes
nothing beyond picking a value, a range over models or a count, and says so where it does; the query dates are the
first and last timestamps of the final trace rows (score.trace_rows). Values are printed as the files hold them,
rounded to three decimals unless the label says otherwise. Counts only: no expert's file, no per-item reading.
"""
from __future__ import annotations

import collections
import datetime as dt
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUT = REPO / 'docs' / 'mock_review_workstreams' / 'reports' / 'NUMBERS_SHEET.md'
SCORE_OF = {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}
RES = HERE / 'results'
L2 = REPO / 'template_annotation_23092026' / 'layer2'
PILOT = REPO / 'evaluator_pilot_17092026'

NAME = {'deepseek-v4.1-flash': 'DeepSeek V4.1 Flash', 'kimi-k3': 'Kimi K3', 'claude-sonnet-5': 'Claude Sonnet 5',
        'glm-5.3-flash': 'GLM-5.3-Flash', 'muse-glimmer-30b': 'Muse Glimmer 30B', 'glm-5.3': 'GLM-5.3',
        'qwen3-235b-a22b-2507': 'Qwen3-235B-2507', 'gemini-3.1-flash-lite': 'Gemini 3.1 Flash-Lite',
        'gemma-4-26b-a4b': 'Gemma 4 26B', 'gpt-5.4-mini': 'GPT-5.4 mini', 'gpt-oss-20b': 'gpt-oss-20b',
        'qwen3-235b-a22b-thinking-2507': 'Qwen3-235B-A22B-Thinking-2507'}


def load(p: Path):
    return json.loads(p.read_text(encoding='utf-8'))


def rel(p: Path) -> str:
    return p.resolve().relative_to(REPO).as_posix()


def f3(x) -> str:
    return '-' if x is None else f'{x:.3f}'


def ci(c) -> str:
    return '-' if not c else f'[{c[0]:.3f}, {c[1]:.3f}]'


def pv(p) -> str:
    return '-' if p is None else (f'{p:.1e}' if p < 0.001 else f'{p:.4f}')


def nm(k: str) -> str:
    return NAME.get(k, k)


def md_rows(text: str, header: str) -> list[list[str]]:
    lines = text.split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith('|') and re.search(header, l)) + 2
    rows = []
    while i < len(lines) and lines[i].startswith('|'):
        rows.append([c.strip().strip('`*') for c in lines[i].strip().strip('|').split('|')])
        i += 1
    return rows


def section_text(text: str, heading: str) -> str:
    m = re.search(rf'^#+ {re.escape(heading)}\n(.*?)(?=^#|\Z)', text, re.S | re.M)
    assert m, heading
    return m.group(1)


class Sheet:
    def __init__(self):
        self.lines = []

    def h(self, title: str):
        self.lines += ['', f'## {title}', '']

    def h3(self, title: str):
        self.lines += ['', f'### {title}', '']

    def put(self, label: str, value, source: str):
        self.lines.append(f'- {label}: **{value}** (`{source}`)')

    def note(self, text: str):
        self.lines.append(text)

    def table(self, header: list[str], rows: list[list], source: str):
        self.lines += ['', '| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
        self.lines += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
        self.lines += ['', f'Source: `{source}`.', '']


def final_ts(path: Path) -> tuple[str, str] | None:
    """First and last timestamp of the final rows (the last answered or empty row per item, as score.trace_rows)."""
    if not path.exists():
        return None
    final = {}
    for ln in path.read_text(encoding='utf-8').splitlines():
        try:
            r = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if r.get('status') in ('answered', 'empty'):
            final[r['item_id']] = r.get('ts')
    ts = sorted(t for t in final.values() if t)
    return (ts[0], ts[-1]) if ts else None


def main() -> int:
    S = Sheet()
    res = load(RES / 'results.json')
    mat = load(RES / 'matched.json')
    cfg = load(RES / 'matched_config.json')
    dec = {d['model_key']: d for d in load(RES / 'decoding_table.json')}
    dec_r = {d['model_key']: d for d in load(RES / 'decoding_table_reasoning-medium-full.json')}
    cov = load(RES / 'coverage_variants.json')
    sv = load(RES / 'sensitivity_variants.json')
    dep = load(RES / 'depth_model.json')
    fp = load(RES / 'flag_precision.json')
    sp = load(RES / 'single_path.json')
    prov = load(RES / 'providers.json')
    js = load(RES / 'judge_swap_main.json')
    rd = load(RES / 'rescore_diff.json')
    for name, d in (('matched.json', mat), ('coverage_variants.json', cov), ('sensitivity_variants.json', sv),
                    ('depth_model.json', dep), ('flag_precision.json', fp), ('single_path.json', sp)):
        assert not d.get('quick') and not d.get('stand_in'), f'{name} is a quick-mode or stand-in file'
    assert res['provenance']['analyze']['B'] == 10000 and res['provenance']['analyze']['B_TEST'] == 100000
    q1 = {m['model']: m for m in res['q1']['models']}
    order = sorted(q1, key=lambda k: -q1[k]['score'])
    q2 = {m['model']: m for m in res['q2']}
    q3 = {m['model']: m for m in res['q3']}
    qc = {m['model']: m for m in res['q3_coverage']['models']}
    R, M = 'full_run_28092026/results/results.json', 'full_run_28092026/results/matched.json'

    ev = load(HERE / 'scores' / 'main' / 'CONFIG.json')['evaluators']
    ev_sha = {k.rsplit('/', 1)[-1]: v['sha256_lf'][:12] for k, v in ev.items()}
    S.lines += ['# Numbers sheet', '',
                'Generated by `python -m full_run_28092026.numbers_sheet` (WS-G, step G5) from the final stores: every '
                f"store re-scored under the final evaluator (`answer.py` {ev_sha['answer.py']}, `milestones.py` "
                f"{ev_sha['milestones.py']}, as `scores/main/CONFIG.json` records them), the judge's new rulings "
                'in, every analysis at full resolution (10,000 bootstrap draws, 100,000 sign flips). Each value names the '
                'file that holds it; regenerate the file, then this sheet, rather than copying a number by hand. Scores '
                'are Final Answer Accuracy (FAC) and Milestone Coverage (MC), as the paper defines them.', '',
                'The paper reports matched settings everywhere (each model at its reasoning setting where the endpoint '
                "offers one, else its default); the three re-run models' default rows are the secondary comparison. The "
                "error analysis, the further conditions and the serving endpoints stay at the providers' defaults by design."]

    # ------------------------------------------------------------------ setup
    S.h('Setup (Section 5.1 and the models appendix)')
    open_ = [k for k in order if dec[k]['weights'] == 'open']
    S.put('Models evaluated', f'{len(q1)} ({len(open_)} open-weights, {len(q1) - len(open_)} closed)',
          'full_run_28092026/results/decoding_table.json: weights')
    rows_main = sum(dec[k]['rows'] for k in q1)
    rerun = [m['model'] for m in cfg['models'] if m.get('run', True) and m.get('reasoning_store')]
    rows_m = sum((dec_r if k in rerun else dec)[k]['rows'] for k in q1)
    S.put('Responses at matched settings (one per instance and model; the three re-run models at reasoning on)',
          f'{rows_m:,}', 'full_run_28092026/results/decoding_table.json, decoding_table_reasoning-medium-full.json: rows '
          '(summed here)')
    S.put('Re-run with reasoning on (matched settings)', ', '.join(nm(k) for k in rerun) + f" at {cfg['effort']}",
          'full_run_28092026/results/matched_config.json: models[].reasoning_store, effort')
    S.put('Responses in the main run at the providers\' defaults (one per instance and model)', f'{rows_main:,}',
          'full_run_28092026/results/decoding_table.json: rows')
    S.put('Reasoning-on re-run responses (the three models, every instance)',
          f"{sum(dec_r[k]['rows'] for k in rerun):,}", 'full_run_28092026/results/decoding_table_reasoning-medium-full.json: rows')
    none_offered = [m['model'] for m in cfg['models'] if m.get('reasoning_setting') == 'none offered']
    S.put('Endpoint offering no reasoning setting', ', '.join(nm(k) for k in none_offered),
          'full_run_28092026/results/matched_config.json: reasoning_setting "none offered"')
    by_default = [m['model'] for m in cfg['models'] if m.get('run', True) and m.get('reasoning_setting') == 'reasons by default']
    S.put('Models reasoning at their providers\' defaults', f'{len(by_default)} (' + ', '.join(nm(k) for k in by_default) + ')',
          'full_run_28092026/results/matched_config.json: reasoning_setting "reasons by default"')
    share = {k: dec[k]['reasoning_tokens']['share_above_zero'] for k in by_default}
    S.put('Their share of main-run responses with reasoning tokens', f'{min(share.values()):.3f} to {max(share.values()):.3f} '
          '(range over the models)', 'full_run_28092026/results/decoding_table.json: reasoning_tokens.share_above_zero')
    twelfth = [m for m in cfg['models'] if m.get('run') is False]
    S.put('Twelfth model', 'none run' + (f" ({', '.join(nm(m['model']) for m in twelfth)} is listed with `run: false`, "
          'not evaluated)' if twelfth else ''), 'full_run_28092026/results/matched_config.json: run')
    ceil = collections.Counter(dec[k]['max_tokens'] for k in q1)
    S.put('Output ceiling', ', '.join(f'{c:,} tokens ({n} models)' for c, n in ceil.most_common()) + '; reasoning-on '
          're-runs ' + ', '.join(f'{c:,}' for c in sorted({d["max_tokens"] for d in dec_r.values()})),
          'full_run_28092026/results/decoding_table*.json: max_tokens')
    mq1 = {m['model']: m for m in mat['q1']['models']}
    unusable_m, unusable = sum(mq1[k]['unusable'] for k in mq1), sum(q1[k]['unusable'] for k in q1)
    S.put('Responses with no readable final answer (score 0), matched settings; default (main run)',
          f'{unusable_m} ({100 * unusable_m / rows_m:.1f}%); {unusable} ({100 * unusable / rows_main:.1f}%)',
          f'{M}: q1.models[].unusable; {R}: q1.models[].unusable')
    S.table(['model (matched order)', 'setting', 'empty, matched', 'no readable answer, matched', 'empty, default',
             'no readable answer, default'],
            [[nm(k), 'reasoning on' if k in rerun else 'default', mq1[k]['empty'], mq1[k]['unusable'], q1[k]['empty'],
              q1[k]['unusable']] for k in mat['order']], f'{M}: q1.models[].empty, unusable; {R}: q1.models[].empty, unusable')
    dates = {}
    for k in q1:
        t = final_ts(HERE / 'traces' / f'{k}.jsonl')
        if t:
            dates[k] = t
    first, last = min(t[0] for t in dates.values()), max(t[1] for t in dates.values())
    S.put('Query dates, main run (first and last final row, UTC)', f'{first[:10]} to {last[:10]}',
          'full_run_28092026/traces/<model>.jsonl: ts (local, gitignored)')
    dr = [final_ts(HERE / 'traces' / 'reasoning-medium-full' / f'{k}.jsonl') for k in rerun]
    dr = [t for t in dr if t]
    if dr:
        S.put('Query dates, reasoning-on re-runs', f'{min(t[0] for t in dr)[:10]} to {max(t[1] for t in dr)[:10]}',
              'full_run_28092026/traces/reasoning-medium-full/<model>.jsonl: ts (local, gitignored)')

    # ------------------------------------------------------------------ results
    S.h('Results (Section 6 and the results appendix)')
    mm, morder, dflt = mat['models'], mat['order'], mat['default']
    rr = [k for k in morder if mm[k].get('configuration') == 'reasoning']
    assert sorted(rr) == sorted(rerun), (rr, rerun)
    cf = lambda k: 'reasoning on' if k in rr else 'default'   # noqa: E731
    ci4 = lambda c: '-' if not c else f'[{c[0]:+.4f}, {c[1]:+.4f}]'   # noqa: E731
    mq2 = {m['model']: m for m in mat['q2']}
    mbl = {m['model']: m for m in mat['branches_levels']}
    mq3 = {k: mm[k]['q3'] for k in morder}
    bl = {m['model']: m for m in res['branches_levels']}
    sig = lambda ps: sum(1 for p in ps if p['p_holm'] < 0.05)   # noqa: E731
    DR = 'Default rows of the three re-run models (the secondary comparison):'
    S.note(f"Headline: matched settings, every table in the matched order (`{M}: order`). Each model runs at its reasoning "
           f"setting where the endpoint offers one, else at its default: {', '.join(nm(k) for k in rr)} at reasoning on "
           f"({cfg['effort']} effort), the other eight at their defaults (the same responses in both configurations). Under "
           'each matched table, the default rows of the three re-run models follow as the secondary comparison.')
    S.h3('Table 1 (matched settings)')
    T1 = ['FAC [95% CI]', 'tier', 'MC [95% CI]', 'MC tier', 'MC by matching alone [95% CI]', 'tier', 'judge-decided share',
          'arithmetic flags on correct answers [95% CI]', 'calculations per response', 'no readable answer (share)']
    S.table(['model (matched order)', 'setting'] + T1,
            [[nm(k), cf(k), f"{f3(m['fac'])} {ci(m['fac_ci'])}", m['letter'], f"{f3(m['mc_strict'])} {ci(m['mc_strict_ci'])}",
              m['mc_letter'], f"{f3(m['mc_e3'])} {ci(m['mc_e3_ci'])}", m['mc_e3_letter'], f3(m['judge_decided_share']),
              f"{f3(m['digit_flag_rate'])} {ci(m['digit_flag_ci'])}", f"{m['claims_per_trace']:.1f}", f3(m['unreadable'])]
             for k, m in ((k, mm[k]) for k in morder)],
            f'{M}: models[] (fac, fac_ci, letter, mc_strict, mc_strict_ci, mc_letter, mc_e3, mc_e3_ci, mc_e3_letter, '
            'judge_decided_share, digit_flag_rate, digit_flag_ci, claims_per_trace, unreadable)')
    S.note(DR)
    e3d = mat['e3_default']['models']
    S.table(['model', 'setting'] + T1,
            [[nm(k), 'default', f"{f3(q1[k]['score'])} {ci(q1[k]['ci'])}", dflt['letters'][k],
              f"{f3(qc[k]['coverage'])} {ci(qc[k]['ci'])}", dflt['mc_letters'][k], f"{f3(e3d[k]['coverage'])} {ci(e3d[k]['ci'])}",
              dflt['e3_letters'][k], f3(q3[k]['e5_judged_fraction']),
              f"{f3(q3[k]['digit_flag_rate_on_fully_solved'])} {ci(q3[k]['digit_ci'])}", f"{q3[k]['claims_per_trace']:.1f}",
              f3(q1[k]['unusable'] / dec[k]['rows'])] for k in rr],
            f'{R}: q1.models[] (score, ci, unusable), q3_coverage.models[] (coverage, ci), q3[] (e5_judged_fraction, '
            f'digit_flag_rate_on_fully_solved, digit_ci, claims_per_trace); {M}: default.letters, mc_letters, e3_letters (tiers '
            'within the default family), e3_default.models; full_run_28092026/results/decoding_table.json: rows (share computed here)')
    S.put('Tiers within the default family, all eleven in the default order (FAC / MC)',
          ', '.join(f"{nm(k)} {dflt['letters'][k]} / {dflt['mc_letters'][k]}" for k in dflt['order']),
          f'{M}: default.order, letters, mc_letters')
    ps = mat['pairs_separated']
    assert (ps['pairs'], ps['fac'], ps['mc'], ps['mc_e3']) == (len(mat['pairs']), sig(mat['pairs']), sig(mat['mc_pairs']),
                                                               sig(mat['mc_e3_pairs'])), ps
    S.put(f"Pairs separated after Holm, matched settings, of {ps['pairs']}: FAC by sign flips; strict FAC (fully solved "
          'instances); McNemar', f"{ps['fac']}; {ps['fac_strict']} (the same call as FAC on {ps['fac_strict_agree']} pairs); "
          f"{ps['fac_mcnemar']}", f'{M}: pairs_separated (fac, fac_strict, fac_strict_agree, fac_mcnemar)')
    S.put('Pairs separated after Holm, matched settings: MC by sign flips; Wilcoxon; MC by matching alone',
          f"{ps['mc']}; {ps['mc_wilcoxon']} (the same call as sign flips on {ps['mc_wilcoxon_agree']} pairs); {ps['mc_e3']}",
          f'{M}: pairs_separated (mc, mc_wilcoxon, mc_wilcoxon_agree, mc_e3)')
    S.put('Pairs separated after Holm, default configuration (secondary; FAC / MC / MC by matching alone, of '
          f"{len(res['q1']['pairs'])})", f"{dflt['pairs_significant']} / {dflt['mc_pairs_significant']} / "
          f"{dflt['e3_pairs_significant']}", f'{M}: default.pairs_significant, mc_pairs_significant, e3_pairs_significant')
    S.put('First five on FAC, matched settings; FAC of the fifth and the sixth to four decimals',
          ', '.join(nm(k) for k in morder[:5]) + f"; {nm(morder[4])} {mm[morder[4]]['fac']:.4f}, {nm(morder[5])} "
          f"{mm[morder[5]]['fac']:.4f}", f'{M}: order, models[].fac')
    S.put('First five on FAC, default (secondary)', ', '.join(nm(k) for k in order[:5]), f'{R}: q1.models[].score (ordered here)')
    t = mat['tau_default_vs_matched']
    S.put("Kendall's tau between the default and matched FAC orderings", f"{t['tau']:.3f} [{t['ci'][0]:.3f}, {t['ci'][1]:.3f}]",
          f'{M}: tau_default_vs_matched')
    S.h3('Paired change per re-run model (reasoning on minus default, all instances; the secondary comparison)')
    S.table(['model', 'FAC default -> on', 'change [95% CI]', '90% CI', 'p', 'detectable (Holm)', 'within the margin',
             'MC change'],
            [[nm(k), f"{f3(c['fac_default'])} -> {f3(c['fac_reasoning'])}", f"{c['change']:+.3f} {ci(c['ci'])}", ci(c['ci90']),
              pv(c.get('p_holm', c.get('p'))), f3(c['detectable_holm']), 'yes' if c['within_margin'] else 'no',
              f"{c['mc_change']:+.3f}" if c.get('mc_change') is not None else '-']
             for k, c in mat['paired_change'].items()], f'{M}: paired_change')
    S.h3('Stability and sensitivity (matched settings; default beside)')
    stm, st = mat['sensitivity']['tau_with_headline'], res['sensitivity']['tau_with_headline']
    tw = lambda *ks: ', '.join(f'{f3(stm[x])} (default {f3(st[x])})' for x in ks)   # noqa: E731
    S.put('Readable-only reordering: tau with FAC as scored when unusable responses are left out', tw('unusable_excluded'),
          f'{M}: sensitivity.tau_with_headline.unusable_excluded; {R}: the same')
    S.put('Tau under half and double the tolerance; counting fully solved instances only',
          tw('half_tol', 'double_tol', 'fully_solved'), f'{M}: sensitivity.tau_with_headline; {R}: the same')
    S.put('Tau without the symbolic templates, without the read-off templates, half-unit window, whole trace',
          tw('without_symbolic_templates', 'without_shortcut_templates', 'half_unit', 'whole_trace'),
          f'{M}: sensitivity.tau_with_headline; {R}: the same')
    tn, tnd = mat['sensitivity']['tau_noise'], res['sensitivity']['tau_noise']
    S.put(f"Tau noise floor (two random halves of each template's instances, {tn['splits']} splits): median, quartiles",
          f"{f3(tn['median'])}, {f3(tn['q1'])} to {f3(tn['q3'])} (default {f3(tnd['median'])}, {f3(tnd['q1'])} to {f3(tnd['q3'])})",
          f'{M}: sensitivity.tau_noise; {R}: the same')
    for lab, spx, src in (('matched settings', {k: mm[k]['single_path'] for k in morder}, f'{M}: models[].single_path'),
                          ('default', sp, 'full_run_28092026/results/single_path.json: <model>')):
        nones = [(k, g, spx[k][g]['none']) for k in morder for g in ('single', 'others')]
        worst = max(nones, key=lambda x: x[2])
        S.put(f"Single-path table, {lab}: share of templates solved on no instance (range over models and both groups; "
              f"{sp['n_single']} single-path, {sp['n_others']} other templates)", f"{min(x[2] for x in nones):.3f} to {worst[2]:.3f} "
              f"(highest: {nm(worst[0])}, {worst[1]})", f'{src}.single/others.none')
    S.h3('Branches, domains, levels (matched settings)')
    br = lambda b: b.replace('_engineering', '')   # noqa: E731

    def bd_row(b, dom):
        lo, hi = (f(b['branch'], key=lambda x: b['branch'][x]['mean']) for f in (min, max))
        dlo = min(dom, key=dom.get)
        return [f"{b['branch'][lo]['mean']:.3f} ({br(lo)})", f"{b['branch'][hi]['mean']:.3f} ({br(hi)})",
                f"{b['branch'][hi]['mean'] - b['branch'][lo]['mean']:.3f}", f'{dom[dlo]:.3f} ({dlo})']
    BH = ['model', 'setting', 'branch means: lowest', 'highest', 'spread', 'lowest domain mean']
    S.table(BH, [[nm(k), cf(k)] + bd_row(mbl[k], mat['reported']['models'][k]['domain']) for k in morder],
            f'{M}: branches_levels[].branch.<branch>.mean (spread computed here); reported.models.<model>.domain (minimum)')
    S.note(DR)
    S.table(BH, [[nm(k), 'default'] + bd_row(bl[k], res['reported']['models'][k]['domain']) for k in rr],
            f'{R}: branches_levels[].branch.<branch>.mean (spread computed here); reported.models.<model>.domain (minimum)')
    yn = lambda p: 'yes' if p < 0.05 else 'no'   # noqa: E731
    lv_row = lambda q, b: [f3(q['easy']), f3(b['level']['Intermediate']['mean']), f3(q['advanced']),   # noqa: E731
                           f"{q['gap']:+.3f} {ci(q['ci'])}", pv(q['p_welch_holm']), pv(q['p_perm_holm']),
                           f"{yn(q['p_welch_holm'])} / {yn(q['p_perm_holm'])}",
                           f"{q['without_two_chemical']['gap']:+.3f}, p {pv(q['without_two_chemical']['p_welch_holm'])}"]
    LH = ['model', 'setting', 'Easy', 'Intermediate', 'Advanced', 'gap Easy - Advanced [95% CI]', 'Welch p (Holm)',
          'permutation p (Holm)', 'gap holds after Holm (Welch / permutation)',
          'without the two chemical templates: gap, Welch p (Holm)']
    S.table(LH, [[nm(k), cf(k)] + lv_row(mq2[k], mbl[k]) for k in morder],
            f'{M}: q2[] (easy, advanced, gap, ci, p_welch_holm, p_perm_holm, without_two_chemical); branches_levels[].level')
    S.note(DR)
    S.table(LH, [[nm(k), 'default'] + lv_row(q2[k], bl[k]) for k in rr],
            f'{R}: q2[] (easy, advanced, gap, ci, p_welch_holm, p_perm_holm, without_two_chemical); branches_levels[].level')
    for lab, qq, B, src in (('matched settings', mq2, mbl, M), ('default (secondary)', q2, bl, R)):
        S.put(f'Models whose Easy-minus-Advanced gap holds after Holm, {lab}: Welch; permutation; Welch without the two '
              'chemical templates', '; '.join(', '.join(nm(k) for k in morder if f(qq[k]) < 0.05) or 'none' for f in (
                  lambda q: q['p_welch_holm'], lambda q: q['p_perm_holm'], lambda q: q['without_two_chemical']['p_welch_holm'])),
              f'{src}: q2[].p_welch_holm, p_perm_holm, without_two_chemical.p_welch_holm')
        hold = [f"{nm(k)}: {br(p['a'])} minus {br(p['b'])} {p['diff']:+.3f}" for k in morder for p in B[k]['pairs']
                if p['p_holm'] < 0.05]
        S.put(f'Branch pairs whose difference holds after Holm, {lab}', '; '.join(hold) or 'none',
              f'{src}: branches_levels[].pairs[] with p_holm < 0.05')
    S.h3('Depth (wrong-answer odds per milestone; matched settings, default beside)')
    DP = 'full_run_28092026/results/depth_model.json'
    for conf, lab in (('matched', 'matched settings'), ('main', 'default (secondary)')):
        c = dep['configurations'][conf]
        holds = [k for k in morder if c['models'].get(k, {}).get('holds')]
        assert c.get('holds', len(holds)) == len(holds), conf
        S.put(f'{lab}: models whose slope holds after Holm', f"{len(holds)} of {len(c['models'])} (" +
              ', '.join(nm(k) for k in holds) + ')', f'{DP}: configurations.{conf}.models[].holds')
        po = c['pooled']
        S.put(f'{lab}: pooled slope per milestone [95% CI], odds ratio, p', f"{po['slope']:+.3f} {ci(po['ci'])}, "
              f"{po.get('odds_ratio', float('nan')):.2f}, {pv(po.get('p'))}", f'{DP}: configurations.{conf}.pooled')
    drow = lambda d: [f"{d['wrong']} / {d['rows']}", f"{d['slope']:+.3f} {ci(d['ci'])}", f"{d['odds_ratio']:.2f}",   # noqa: E731
                      pv(d['p_holm']), 'yes' if d['holds'] else 'no']
    DH = ['model', 'setting', 'wrong / rows', 'slope per milestone [95% CI]', 'odds ratio', 'p (Holm)', 'holds']
    dmm, dmd = dep['configurations']['matched']['models'], dep['configurations']['main']['models']
    S.table(DH, [[nm(k), cf(k)] + drow(dmm[k]) for k in morder if k in dmm], f'{DP}: configurations.matched.models')
    S.note(DR)
    S.table(DH, [[nm(k), 'default'] + drow(dmd[k]) for k in rr if k in dmd], f'{DP}: configurations.main.models')
    S.h3('Milestone Coverage under four readings (matched settings, all responses)')
    CV = 'full_run_28092026/results/coverage_variants.json'
    rd_ = ('as_scored', 'matching_only', 'route_adjusted', 'intermediate_only')
    cvm = cov['matched']
    assert cvm['stores'] == {k: mm[k]['store'] for k in cvm['stores']}, cvm['stores']
    cvr = lambda d: [f"{f3(d[r]['all'])} {ci(d[r].get('ci'))}" for r in rd_]   # noqa: E731
    S.table(['model', 'setting'] + [f"{r} ({cvm['templates'][r]} templates)" for r in rd_],
            [[nm(k), cf(k)] + cvr(cvm['models'][k]) for k in morder], f'{CV}: matched.models.<model>.<reading>.all, ci; matched.templates')
    S.note(DR)
    S.table(['model', 'setting'] + [f"{r} ({cov['templates'][r]} templates)" for r in rd_],
            [[nm(k), 'default'] + cvr(cov['models'][k]) for k in rr], f'{CV}: models.<model>.<reading>.all, ci; templates')
    for lab, fam, pre in (('matched settings', cvm, 'matched.'), ('default (secondary)', cov, '')):
        sepr = fam['separation_rule']['claude_vs_deepseek']
        S.put(f'Claude Sonnet 5 minus DeepSeek V4.1 Flash on MC, Holm p per reading, {lab}',
              '; '.join(f"{r}: {sepr['diff'][r]:+.3f}, p {pv(sepr[r])}" for r in rd_) + '; the separation holds under every '
              f"reading: {'yes' if sepr['holds'] else 'no'}", f'{CV}: {pre}separation_rule.claude_vs_deepseek')
        S.put(f'Pairs separated after Holm per reading, {lab}', ', '.join(
            f"{r} {sum(1 for p in fam['pairs'][r] if p['p_holm'] < 0.05)}" for r in rd_ if r in fam['pairs']),
            f'{CV}: {pre}pairs.<reading>[] (counted here)')
        vb = fam['verbosity']
        S.put(f"Verbosity, {lab}: MC slope ({vb['units']['numbers']}), pooled [95% CI]; within model and template [95% CI]; "
              'responses', f"{vb['slope_numbers']:+.4f} {ci4(vb['ci'])}; {vb['within_model']['slope_numbers']:+.4f} "
              f"{ci4(vb['within_model']['ci'])}; {vb['responses']:,}", f'{CV}: {pre}verbosity.slope_numbers, ci, within_model, responses')
    S.put('Judge-decided share of milestones (range over models), matched settings',
          f"{min(mm[k]['judge_decided_share'] for k in morder):.3f} to {max(mm[k]['judge_decided_share'] for k in morder):.3f} "
          f"(default {min(q3[k]['e5_judged_fraction'] for k in order):.3f} to {max(q3[k]['e5_judged_fraction'] for k in order):.3f})",
          f'{M}: models[].judge_decided_share; {R}: q3[].e5_judged_fraction')
    S.h3('Judge swap (a second judge on a sample of what the judge is sent)')
    S.note("- At the providers' defaults by design: the sample is drawn from the main run's judge jobs.")
    S.put('Responses compared', f"{js['traces']} ({js['traces_per_model'][0]} to {js['traces_per_model'][1]} per model, "
          f"{js['templates_per_model'][0]} to {js['templates_per_model'][1]} templates per model); of {js['sample']} drawn, "
          f"{js['not_in_rows']} re-run since and {js['prompt_changed']} with a changed prompt are left out",
          'full_run_28092026/results/judge_swap_main.json: traces, traces_per_model, templates_per_model, sample, not_in_rows, prompt_changed')
    po = js['pooled']
    S.put('Pooled agreement, three-way and reached against not; Cohen\'s kappa, both', f"{po['agreement_three_way']:.3f}, "
          f"{po['agreement_reached_vs_not']:.3f}; {po['kappa_three_way']:.3f}, {po['kappa_reached_vs_not']:.3f} over "
          f"{po['milestones']} milestones", 'full_run_28092026/results/judge_swap_main.json: pooled')
    diffs = [m['diff'] for m in js['models'].values()]
    S.put('MC shift under the second judge (range over models)', f'{min(diffs):+.3f} to {max(diffs):+.3f}',
          'full_run_28092026/results/judge_swap_main.json: models.<model>.diff')
    S.h3('Arithmetic flags and the judged step check (matched settings; default beside)')
    FP = 'full_run_28092026/results/flag_precision.json'
    FK = ('flags', 'carried_precision', 'other_truncated', 'other_no_upstream', 'other_not_reproduced')
    fps = fp['stores']
    assert fps['main']['all'] == fp['all'] and all(sum(fps['main']['models'][k][x] for k in morder) == fp['all'][x] for x in FK)
    fm = collections.Counter()
    for k in morder:
        for x in FK:
            fm[x] += fps['reasoning-medium-full' if k in rr else 'main']['models'][k][x]
    fr = fps['reasoning-medium-full']['all']
    cp = lambda a: (f"{a['carried_precision']} of {a['flags']} ({100 * a['carried_precision'] / a['flags']:.1f}%); other: "   # noqa: E731
                    f"truncated {a.get('other_truncated')}, no upstream value {a.get('other_no_upstream')}, not reproduced "
                    f"{a.get('other_not_reproduced')}")
    S.put('Arithmetic flags on correct answers that are carried precision, matched settings', cp(fm),
          f'{FP}: stores.main.models (the eight), stores.reasoning-medium-full.models (the re-run three); summed here')
    S.put('The same, the three re-run models at reasoning on', cp(fr), f'{FP}: stores.reasoning-medium-full.all')
    S.put('The same, default (main run, all eleven)', cp(fp['all']), f'{FP}: all')
    ex = fp['expert_confirmed']
    S.put('Slips the domain expert confirmed that are carried precision (expert readings of default-setting responses)',
          f"{ex['carried_precision']} of {ex['slips']} ({ex['step_changed']} in steps whose text changed with the repaired items)",
          f'{FP}: expert_confirmed')
    q3h = ['model', 'setting', 'arithmetic flags on correct answers [95% CI]',
           'any step flagged on correct answers, arithmetic or judge [95% CI] (the judge alone: tab:judged_steps)',
           'any step flagged on wrong answers [95% CI]', 'steps flagged per response']
    q3r = lambda q: [f"{f3(q['digit_flag_rate_on_fully_solved'])} {ci(q['digit_ci'])}",   # noqa: E731
                     f"{f3(q['router_rate_on_fully_solved'])} {ci(q['router_ci'])}",
                     f"{f3(q['router_rate_on_wrong'])} {ci(q['router_wrong_ci'])}", f"{q['router_steps_flagged_per_trace']:.2f}"]
    S.table(q3h, [[nm(k), cf(k)] + q3r(mq3[k]) for k in morder],
            f'{M}: models[].q3 (digit_flag_rate_on_fully_solved, router_rate_on_fully_solved, router_rate_on_wrong, ...)')
    S.note(DR)
    S.table(q3h, [[nm(k), 'default'] + q3r(q3[k]) for k in rr],
            f'{R}: q3[] (digit_flag_rate_on_fully_solved, router_rate_on_fully_solved, router_rate_on_wrong, ...)')
    rv = {r[0]: r for r in md_rows((HERE / 'ROUTER_VALIDATION.md').read_text(encoding='utf-8'), r'^\| steps the experts call incorrect')}
    for row in ('all traces, the router', 'inside correct-answer traces, the router'):
        S.put(f'Judged step and arithmetic checks against the experts\' step labels on the pilot\'s traces, {row.split(",")[0]}: '
              'tp / fp / fn, precision, recall, F1', f'{rv[row][1]}, {rv[row][2]}, {rv[row][3]}, {rv[row][4]}',
              'full_run_28092026/ROUTER_VALIDATION.md')
    S.h3('Scoring-rule variants (matched settings; verdicts up / down)')
    SV = 'full_run_28092026/results/sensitivity_variants.json'
    svs, svm = sv['stores']['matched'], sv['stores']['main']
    V = ('abs_clause_off', 'last_digit_unbounded', 'prescribed_relaxed', 'per_part_credit')
    vrow = lambda d, k: [f3(d['models'][k]['headline'])] + [   # noqa: E731
        f"{f3(d['models'][k][v])} ({d['models'][k]['changed'][v]['up']}/{d['models'][k]['changed'][v]['down']})" for v in V]
    S.table(['model', 'setting', 'FAC'] + list(V), [[nm(k), cf(k)] + vrow(svs, k) for k in morder], f'{SV}: stores.matched.models')
    S.note(DR)
    S.table(['model', 'setting', 'FAC'] + list(V), [[nm(k), 'default'] + vrow(svm, k) for k in rr], f'{SV}: stores.main.models')
    for v in V:
        w, wd = svs['where_changed'][v], svm['where_changed'][v]
        S.put(f'{v}: verdicts moved, templates; tau with FAC as scored (matched settings; default beside)',
              f"{w['verdicts']} on {w['templates']} templates; {f3(svs['tau_with_headline'][v])} (default: {wd['verdicts']} on "
              f"{wd['templates']} templates; {f3(svm['tau_with_headline'][v])})",
              f'{SV}: stores.matched.where_changed.{v}, tau_with_headline; stores.main, the same')
    BN = ('le_0.2', '0.2_1', '1_5', 'gt_5')
    bm = collections.Counter()
    for k in morder:
        for x in BN + ('numeric', 'symbolic'):
            bm[x] += svs['relative_error_bins'][k][x]
    S.put('Accepted numeric answers by relative error, matched settings: within 0.2%, 0.2 to 1%, 1 to 5%, above 5%',
          ', '.join(f"{bm[n]} ({100 * bm[n] / bm['numeric']:.1f}%)" for n in BN) + f" of {bm['numeric']}; "
          f"{bm['symbolic']} more accepted by the symbolic step",
          f'{SV}: stores.matched.relative_error_bins.<model> (summed here, shares computed here; the file holds no matched total)')
    b = sv['relative_error_bins']['all']
    S.put('The same, default (main run)', ', '.join(f"{b[n]} ({100 * b['share'][n]:.1f}%)" for n in BN) + f" of {b['numeric']}; "
          f"{b['symbolic']} more accepted by the symbolic step", f'{SV}: relative_error_bins.all')
    S.h3('Repeats, paraphrases and the further conditions')
    rpm, rp = mat['repeats'], res['repeats']
    no_rep = [k for k in rr if k not in rpm and k not in rp]
    S.put('Decoding repeats at matched settings: models and instances', ', '.join(f'{nm(k)} ({cf(k)})' for k in rpm) +
          f"; {', '.join(sorted({str(v['items']) for v in rpm.values()}))} instances each, "
          f"{', '.join(sorted({str(len(v['scores'])) for v in rpm.values()}))} repeats" +
          (f"; no repeat runs for {', '.join(nm(k) for k in no_rep)}" if no_rep else ''), f'{M}: repeats; {R}: repeats')
    RH = ['model', 'setting', 'repeat scores', 'SD', 'same verdict in every repeat']
    rprow = lambda v: [', '.join(f3(s) for s in v['scores'].values()), f"{v['sd']:.4f}", f3(v['same_verdict_every_repeat'])]   # noqa: E731
    S.table(RH, [[nm(k), cf(k)] + rprow(v) for k, v in rpm.items()], f'{M}: repeats')
    S.note(DR)
    S.table(RH, [[nm(k), 'default'] + rprow(rp[k]) for k in rr if k in rp], f'{R}: repeats')
    q5m, q5 = mat['q5'], res['q5']
    S.put('Paraphrase pairs kept by the experts', q5m['expert_check'], f'{M}: q5.expert_check; full_run_28092026/PARAPHRASE_REVIEW.md')
    Q5H = ['change (paraphrase minus original) [95% CI]', '90% CI', 'p (Holm)', f"within +/-{q5m['margin']}"]
    q5r = lambda m: [f"{m['diff']:+.3f} {ci(m['ci'])}", ci(m['ci90']), pv(m['p_holm']), 'yes' if m['within_margin'] else 'no']   # noqa: E731
    q5mm, q5d = {m['model']: m for m in q5m['models']}, {m['model']: m for m in q5['models']}
    S.table(['model', 'setting (paraphrase store)'] + Q5H, [[nm(k), f"{cf(k)} ({q5m['stores'][k]})"] + q5r(q5mm[k]) for k in morder],
            f'{M}: q5.models, q5.stores')
    S.note(DR)
    S.table(['model', 'setting'] + Q5H, [[nm(k), 'default'] + q5r(q5d[k]) for k in rr], f'{R}: q5.models')
    for lab, qq, src in (('matched settings', q5m, M), ('default (secondary)', q5, R)):
        inm = [m for m in qq['models'] if m['within_margin']]
        S.put(f"Paraphrase, {lab}: within +/-{qq['margin']}; lowest p (Holm)",
              f"{len(inm)} of {len(qq['models'])}; {pv(min(m['p_holm'] for m in qq['models']))}", f'{src}: q5.models[].within_margin, p_holm')
        t, na = qq['tau'], qq['tau']['noise_arm']
        S.put(f'Paraphrase-against-original tau [95% CI], {lab}; noise floor at the size and template mix of the paraphrase arm: '
              'median (quartiles)', f"{t['tau']:.3f} {ci(t['ci'])}; {na['median']:.3f} ({na['q1']:.3f} to {na['q3']:.3f})",
              f'{src}: q5.tau, q5.tau.noise_arm')
    S.note("- The further conditions stay at the providers' defaults by design: each against the main run on the same instances.")
    S.table(['condition', 'model', 'instances', 'main on the same instances', 'under the condition', 'change [95% CI]'],
            [[a['arm'], nm(a['model']), a['items'], f3(a['main_score_on_items']), f3(a['arm_score']), f"{a['diff']:+.3f} {ci(a['ci'])}"]
             for a in res['reasoning_arms']], f'{R}: reasoning_arms (the subset includes the repaired items, re-run in every condition)')

    # ------------------------------------------------------------------ error analysis
    S.h('Error analysis (B2, the expert reading of wrong answers)')
    S.note("- At the providers' defaults by design: the experts read main-run responses.")
    er =(HERE / 'EXPERT_REQUEST.md').read_text(encoding='utf-8')
    b2 = section_text(er, 'B2 now')
    sample = re.search(r'The sample: (.*?); (\d+) items on (\d+) templates', b2)
    S.put('Sample (default-setting responses; no re-reading for the re-run models)', f'{sample.group(2)} items on {sample.group(3)} '
          f'templates: {sample.group(1)}', 'full_run_28092026/EXPERT_REQUEST.md: B2 now')
    kap = re.search(r'Readings (\d+) on (\d+) items; Fleiss\' kappa over the three readers, all models: (0\.\d+)', b2)
    S.put('Readings, items, Fleiss\' kappa', f'{kap.group(1)}, {kap.group(2)}, {kap.group(3)}', 'full_run_28092026/EXPERT_REQUEST.md: B2 now')
    S.table(['model', '"No error" share', 'Fleiss kappa', 'item majorities'],
            [[r[0], r[3], r[4], r[2]] for r in md_rows(b2, r'^\| model \| readings by category')],
            'full_run_28092026/EXPERT_REQUEST.md: B2 now')
    b4 = [r for r in md_rows(er, r'^\| template \| ') if len(r) > 1 and r[1] == 'near'] if '| near |' in er else []
    near = sorted({r[0] for r in b4})
    S.put('Templates whose wording does not pin the answer, still in the set (B4 "near" templates)', f'{len(near)}: ' + ', '.join(near),
          'full_run_28092026/EXPERT_REQUEST.md: B4, the templates')

    # ------------------------------------------------------------------ certification
    S.h('Certification (the certification appendix)')
    cert = (L2 / 'CERTIFICATION.md').read_text(encoding='utf-8')
    rounds = md_rows(cert, r'^\| Round \| Labels')
    S.put('Rounds; templates and verdicts per round', f'{len(rounds)}; ' + '; '.join(f'round {r[0]}: {r[3]} templates, {r[4]} verdicts'
          for r in rounds), 'template_annotation_23092026/layer2/CERTIFICATION.md')
    status = dict((r[0], r[1]) for r in md_rows(cert, r'^\| \| Templates'))
    S.put('Certified', status.get('certified'), 'template_annotation_23092026/layer2/CERTIFICATION.md')
    for n in (5, 6):
        t5 = (L2 / f'RESULTS_round{n}.md').read_text(encoding='utf-8')
        verdicts = '; '.join(f'{r[0]} {r[3]} approved, {r[4]} rejected' for r in md_rows(t5, r'^\| Expert \| Branch \| Items'))
        hc = re.search(r'(\d+) hand checks with a comparable number: (\d+) matched', t5)
        S.put(f'Round {n}: verdicts per expert; hand checks matched', f'{verdicts}; {hc.group(2)} of {hc.group(1)}',
              f'template_annotation_23092026/layer2/RESULTS_round{n}.md')
    r1 = (L2 / 'RESULTS.md').read_text(encoding='utf-8')
    hc1 = re.search(r'(\d+) hand checks with a comparable number: (\d+) matched', r1)
    lab = re.search(r'from (\d+) label rows by (\d+) experts', r1)
    S.put('Round 1 hand checks: with a comparable number, of all; matched', f'{hc1.group(1)} of {lab.group(1)} (the other '
          f'{int(lab.group(1)) - int(hc1.group(1))} give no number to compare); {hc1.group(2)} matched',
          'template_annotation_23092026/layer2/RESULTS.md')
    S.note('- What rounds 5 and 6 re-certified: the two chemical templates whose questions now state their reading and their '
           'data (round 5), then the virial template with organic compounds held below their decomposition temperature (round 6) '
           '(`docs/mock_review_workstreams/reports/WS-A_report.md`).')
    ag = md_rows(r1, r'^\| Branch \| Templates with 3 verdicts')
    undefined = [r[0] for r in ag if r[0] != 'all' and r[4] == '100%']
    S.put('Fleiss\' kappa undefined (every verdict approves, 0/0; the report prints 1.000)', ', '.join(undefined),
          'template_annotation_23092026/layer2/RESULTS.md; docs/appendix_certification.py prints "--"')

    # ------------------------------------------------------------------ difficulty
    S.h('Difficulty labels')
    man = [json.loads(l) for l in (HERE / 'manifest.jsonl').read_text(encoding='utf-8').splitlines() if l.strip()]
    lv = {}
    for r in man:
        lv[r['template_id']] = (r['branch'], r['level'])
    tot = collections.Counter(v[1] for v in lv.values())
    S.put('Templates per level (the domain experts\' labels, unchanged)', f"Easy {tot['Easy']}, Intermediate {tot['Intermediate']}, "
          f"Advanced {tot['Advanced']} of {len(lv)}", 'full_run_28092026/manifest.jsonl (also `python -m full_run_28092026.export_levels`)')
    S.table(['branch', 'Easy', 'Intermediate', 'Advanced'],
            [[b.replace('_engineering', ''), *(sum(1 for v in lv.values() if v == (b, L)) for L in ('Easy', 'Intermediate', 'Advanced'))]
             for b in sorted({v[0] for v in lv.values()})], 'full_run_28092026/manifest.jsonl')
    S.note('- Labelling procedure: `template_annotation_23092026/levels/difficulty_labelling_protocol.md`. No agreement figure '
           '(D7).')

    # ------------------------------------------------------------------ symbolic
    S.h('Symbolic answers')
    ans = (REPO / 'evaluator_pilot_17092026' / 'evaluators' / 'answer.py').read_text(encoding='utf-8')
    enabled = re.findall(r"'template_(\w+)'", re.search(r'SYMBOLIC_EQUIVALENCE_TEMPLATES = \((.*?)\)', ans, re.S).group(1))
    symbolic = sorted({r['template_id'].replace('template_', '') for r in man if r['answer_type'] == 'symbolic'})
    S.put('Templates with symbolic answers; with the equivalence step enabled', f"{len(symbolic)}; {len(enabled)} ("
          + ', '.join(enabled) + ')', 'full_run_28092026/manifest.jsonl: answer_type; answer.py: SYMBOLIC_EQUIVALENCE_TEMPLATES')
    S.put('Scored by the numbers they state', ', '.join(t for t in symbolic if t not in enabled), 'the same')
    S.put('Own-digit cap (an answer\'s last digit vouches for at most this share of the target)',
          re.search(r'^OWN_DIGIT_CAP = ([\d.]+)', ans, re.M).group(1), 'answer.py: OWN_DIGIT_CAP')
    sc = (HERE / 'symbolic' / 'SYMBOLIC_CHECK.md').read_text(encoding='utf-8')
    g = re.search(r'^\| all \| (\d+) \| \d+ \| \d+ \| \d+ \| \d+ \| (\d\.\d+) \| (\d\.\d+) \|', sc, re.M)
    e = re.search(r'^\| all \| (\d+) \| (\d\.\d+) \| (\d\.\d+) \| \d+ \| \d+ \| \d+ \| \d+ \|', sc, re.M)
    S.put('Validation on the graded verdicts it can move (in-sample: the rule was settled on them): precision, recall, verdicts',
          f'{g.group(2)}, {g.group(3)}, {g.group(1)}', 'full_run_28092026/symbolic/SYMBOLIC_CHECK.md')
    S.put('On the earlier B1 and B2 readings (made before the rule): precision, recall, readings', f'{e.group(2)}, {e.group(3)}, '
          f'{e.group(1)}', 'full_run_28092026/symbolic/SYMBOLIC_CHECK.md')
    main_store = next(s for s in rd['stores'] if s['variant'] == 'main')
    eleven = [k for k in order]
    mv = collections.Counter()
    for k in eleven:
        for m in main_store['models'][k]['moves']:
            mv[(m['from'], m['to'], m['cause'])] += m['n']
    up = lambda a, b: SCORE_OF[b] > SCORE_OF[a]   # noqa: E731
    sym = sum(n for (_a, _b, c), n in mv.items() if c == 'symbolic step')
    rule_down = sum(n for (a, b, c), n in mv.items() if c == 'number rule' and not up(a, b))
    rule_up = sum(n for (a, b, c), n in mv.items() if c == 'number rule' and up(a, b))
    S.put('Verdicts the final evaluator moved over the eleven models (main run, against the store before ANSWER FINAL)',
          f'{sym + rule_down + rule_up}: {sym} raised by the symbolic step; {rule_down} lowered and {rule_up} raised by the '
          'number rule (the last-digit bounds, and the answer targets the milestone matcher feeds)',
          f"full_run_28092026/results/rescore_diff.json (since {rd.get('since')}): stores[main].models[].moves")
    S.put('Change in mean score per model (main run)', ', '.join(f"{nm(k)} {main_store['models'][k]['mean_new'] - main_store['models'][k]['mean_old']:+.4f}"
          for k in eleven), 'full_run_28092026/results/rescore_diff.json: mean_old, mean_new')

    # ------------------------------------------------------------------ validation
    S.h('Evaluator validation (the validation appendix)')
    svd = (HERE / 'SCORER_VALIDATION.md').read_text(encoding='utf-8')
    now = {r[0]: r[4] for r in md_rows(svd, r'^\| figure \| published \| published code')}
    S.table(['figure', 'the final evaluator'], [[k, v] for k, v in now.items()], 'full_run_28092026/SCORER_VALIDATION.md (the "now" column)')
    x1 = (PILOT / 'RESULTS_X1.md').read_text(encoding='utf-8')
    m17 = re.search(r'(\d+) of its (\d+) calls never returned', x1)
    if m17:
        S.put('The step-matching judge\'s calls that never returned after retries (MiMo)', f'{m17.group(1)} of {m17.group(2)}',
              'evaluator_pilot_17092026/RESULTS_X1.md')
    psum = (PILOT / 'PILOT_SUMMARY.md').read_text(encoding='utf-8').replace('**', '')
    f1j = re.search(r'\| deterministic, then a judge on the residue \(E5\) \| (0\.\d+) \| (0\.\d+) \| (0\.\d+) \|', psum)
    S.put('Milestone F1: matching alone (final evaluator) and with the judge (pilot)', f"{now.get('E3 F1')} and {f1j.group(3)}",
          'full_run_28092026/SCORER_VALIDATION.md; evaluator_pilot_17092026/PILOT_SUMMARY.md')
    b1 = section_text(er, 'B1 now')
    b3 = section_text(er, 'B3 now')
    S.table(['reading', 'the experts said'],
            [[r[0], r[1]] for r in md_rows(b1, r'^\| \| \|') if r[0]] + [[r[0], r[1]] for r in md_rows(b3, r'^\| \| \|') if r[0]],
            'full_run_28092026/EXPERT_REQUEST.md: B1 now, B3 now (the readings the store still holds; B1 against the current verdicts)')

    # ------------------------------------------------------------------ the answer check by answer kind
    aud = load(HERE / 'results' / 'answer_audit_paper.json')
    src = 'full_run_28092026/results/answer_audit_paper.json (answer_audit.py --write)'
    S.h3('The final-answer check by answer kind (tab:answer_kinds) and its discrepancies')
    S.note("- On the main run (the providers' defaults) and the experts' readings of it; the file holds no matched-settings "
           'version of the bound below.')
    S.table(['answer kind', 'study responses', 'three-way', 'not partial', 'readings agreeing'],
            [[k, v['study'][1], f"{v['study'][0] / v['study'][1]:.3f}", f"{v['study'][2] / v['study'][3]:.3f}",
              f"{v['readings'][0]} of {v['readings'][1]}"] for k, v in aud['by_kind'].items()], src)
    dr = aud['readings_direction']
    S.put('Readings that differ from the check: the check credits less / more / experts split or judge partial credit',
          f"{sum(dr.values())}: {dr['strict']} / {dr['lenient']} / {dr['split']}", src)
    S.table(['check > experts', 'causes (readings)'],
            [[c, ', '.join(f'{k} {n}' for k, n in sorted(v.items(), key=lambda kv: -kv[1]))]
             for c, v in aud['readings_by_category'].items()], src)
    S.put('Templates of other kinds whose question also asks for a label no target checks', aud['label_part_templates'], src)
    b = aud['bound']
    S.put('Every stated part and label part a target, unit factor only where the written unit allows it, absolute value '
          'only where the sign is a convention, targets near zero without their window, notation and label wording read: '
          'largest FAC change; bound; tau; verdicts moved; pairs reversed (gap as scored)',
          f"{b['combined']['largest_change']:+.4f}; no more than {b['combined']['at_most']:.3f}; {b['combined']['tau']:.3f}; "
          f"{b['combined']['moved']}; " + (', '.join(f'{x} and {y} ({g:.4f})' for x, y, g in b['combined']['swaps']) or 'none'),
          src + ' from results/answer_audit.json (matched settings)')
    S.put('Requiring a correct rounding (half a unit): largest FAC change; rounded; tau; verdicts moved; pairs reversed',
          f"{b['half_unit']['largest_change']:+.4f}; {b['half_unit']['largest_rounded']:.3f}; {b['half_unit']['tau']:.3f}; "
          f"{b['half_unit']['moved']}; " + (', '.join(f'{x} and {y} ({g:.4f})' for x, y, g in b['half_unit']['swaps']) or 'none'),
          src + ' from results/answer_audit.json (matched settings)')
    a, c = aud['cases']['lenient'], aud['cases']['strict']
    S.put('tab:answer_cases, credited more than the experts (their readings are in the generated table, not here)',
          f"{a['item']}: {a['stated']} for {a['gold']} (unrounded {a['unrounded']}); check {a['check']}", src)
    S.put('tab:answer_cases, credited less than the experts (their readings are in the generated table, not here)',
          f"{c['item']}: cos(pi n/{c['den']}) for cos({c['gold_num']} pi n/{c['den']}); check {c['check']} "
          f"({c['matched']} of {c['of']})", src)

    # ------------------------------------------------------------------ providers
    S.h('Serving endpoints')
    S.note("- At the providers' defaults by design: the main run's endpoints.")
    well = [r for r in prov if r.get('matched_diff') is not None and not r['few_matched']]
    lo, hi = min(well, key=lambda r: r['matched_diff']), max(well, key=lambda r: r['matched_diff'])
    S.put('Template-matched difference over endpoints matched on at least 20 templates (range)', f"{lo['matched_diff']:+.3f} ({nm(lo['model'])}, "
          f"{lo['endpoint']}) to {hi['matched_diff']:+.3f} ({nm(hi['model'])}, {hi['endpoint']}), over {len(well)} endpoints",
          'full_run_28092026/results/providers.json: matched_diff, few_matched')
    S.put('Endpoints serving fewer than 20 templates; with a matched difference on fewer than 20', f"{sum(r['few_templates'] for r in prov)}; "
          f"{sum(r['few_matched'] for r in prov)} (of {len(prov)} endpoint rows)", 'full_run_28092026/results/providers.json')

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('\n'.join(S.lines) + '\n', encoding='utf-8', newline='\n')
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    print(f'wrote {rel(OUT)} ({len(S.lines)} lines)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
