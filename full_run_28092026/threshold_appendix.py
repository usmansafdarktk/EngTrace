"""The threshold appendix: every tolerance and cut-off in the stack, with the sensitivity each was measured under,
collated from the scripts that measured them (docs/EVALUATION_NEXT_STEPS.md A5; D-175).

    python -m full_run_28092026.threshold_appendix            # FREE: runs the pilot's analyses, reads results.json;
                                                              #       writes THRESHOLD_APPENDIX.md, results/threshold_appendix.json
    python -m full_run_28092026.threshold_appendix --cached   # FREE: renders again from the outputs captured last time

No number here is typed in. The pilot's scripts are run as they stand, under the pilot's pinned environment and, where
milestones are derived from templates, under `pinned_templates` (the templates as they were at the freeze), and the
tables they print are parsed; each parse asserts the rows it expects. The full run's figures are read from
results/results.json and gold_validation.json. The captured outputs are kept under scores/_threshold_logs/, local,
because one of them prints claim text.

WHAT IT COLLATES.
  answer check   the relative tolerance (answer.REL), fitted once and reported split-half on the pilot
                 (analysis/answer_check.py); the full run under half and double that tolerance, the half-unit window
                 and the whole-trace reading (results.json, Sensitivity), with Kendall's tau against the headline.
  E3             the milestone tolerance (milestones.DISPLAY_TOL) and unit scaling: real coverage, null coverage,
                 their separation and the correlation with the framework's answer check at 2%, 1%, 0.5% and 0.2%,
                 with and without unit scaling (analysis/e3_grid.py under the pin); the full run's chance floor on
                 the gold (gold_validation.json) and per model (results.json, Q3).
  digit rule     the four readings of the arithmetic checker on the pilot's labels, 1%, 0.1%, the bare digit rule
                 and the rule as shipped, step level inside correct-answer traces and trace level
                 (analysis/digit_rule.py); the full run's flag rates under the shipped rule and at 1% (results.json).
  PRM            E2's 0.5 cut-off fitted on one half of the traces and reported on the other
                 (analysis/prm_threshold.py).
  E0             the published framework's own constants, read from its source, for the record.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PILOT = REPO / 'evaluator_pilot_17092026'
VENV = PILOT / '.venv' / 'Scripts' / 'python.exe'
LOGS = HERE / 'scores' / '_threshold_logs'
OUT_MD = HERE / 'THRESHOLD_APPENDIX.md'
OUT_JSON = HERE / 'results' / 'threshold_appendix.json'

RUNS = {
    'e3_grid': ['-m', 'evaluator_pilot_17092026.pinned_templates', 'analysis.e3_grid'],
    'digit_rule': [str(PILOT / 'analysis' / 'digit_rule.py')],
    'prm_threshold': [str(PILOT / 'analysis' / 'prm_threshold.py')],
    'answer_check': [str(PILOT / 'analysis' / 'answer_check.py')],
}


def run(name: str, cached: bool) -> str:
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS / f'{name}.txt'
    if cached and log.exists():
        return log.read_text(encoding='utf-8')
    py = str(VENV) if VENV.exists() else sys.executable
    p = subprocess.run([py] + RUNS[name], cwd=str(REPO), capture_output=True, text=True, encoding='utf-8',
                       errors='replace', env={**__import__('os').environ, 'PYTHONIOENCODING': 'utf-8'})
    text = p.stdout + ('\n[stderr]\n' + p.stderr if p.returncode else '')
    log.write_text(text, encoding='utf-8')
    if p.returncode:
        raise SystemExit(f'{name} failed ({p.returncode}); see {log}\n{p.stderr[-2000:]}')
    return p.stdout


def must(pattern: str, text: str, flags=re.M) -> list:
    found = re.findall(pattern, text, flags)
    if not found:
        raise SystemExit(f'expected a line matching {pattern!r}; not found')
    return found


def parse_e3(text: str) -> list[dict]:
    rows = must(r'^\s*([\d.]+)%\s+(all|none)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+|nan)\s*$', text)
    return [{'tolerance': float(t) / 100, 'unit_scaling': sc == 'all', 'real': float(r), 'null': float(n),
             'separation': float(s), 'corr_with_fac': float(c)} for t, sc, r, n, s, c in rows]


def parse_digit(text: str) -> dict:
    out = {}
    for title, key in (('ALL 300 TRACES', 'all'), ('INSIDE CORRECT-ANSWER TRACES - the hard case', 'correct_answer')):
        block = text.split('STEP LEVEL, ' + title, 1)[1].split('\n\n', 1)[0]
        rows = must(r'^\s+(tol1|tol01|digit|e4)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', block)
        out[key] = {r[0]: {'tp': int(r[1]), 'fp': int(r[2]), 'fn': int(r[3]), 'precision': float(r[4]),
                           'recall': float(r[5]), 'f1': float(r[6])} for r in rows}
    block = text.split('TRACE LEVEL, correct-answer traces', 1)[1].split('\n\n', 1)[0]
    rows = must(r'^\s+(tol1|tol01|digit|e4)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+\(([\d.]+), ([\d.]+)\)', block)
    out['trace_level'] = {r[0]: {'flagged': int(r[1]), 'precision': float(r[2]), 'recall': float(r[3]), 'f1': float(r[4]),
                                 'auroc': float(r[5]), 'ci': [float(r[6]), float(r[7])]} for r in rows}
    m = must(r'^(\d+) traces scored, (\d+) skipped', text)[0]
    out['traces_scored'], out['traces_skipped'] = int(m[0]), int(m[1])
    return out


def parse_prm(text: str) -> list[dict]:
    block = text.split('VERDICT - does calibrating the threshold buy anything, held out?', 1)[1]
    rows = must(r'^\s+(.+?) / (qwen72|versa|qwen7)\s+([\d.]+)\s+([\d.]+)\s+([-+][\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', block)
    return [{'subset': s.strip(), 'prm': p, 'f1_at_fitted': float(a), 'f1_at_0.5': float(b), 'gain': float(g),
             'fitted_threshold_halves': [float(t0), float(t1)]} for s, p, a, b, g, t0, t1 in rows]


def parse_answer(text: str) -> dict:
    halves = must(r'fit on half (\d) -> ([\d.]+); held out on half (\d): ([\d.]+)', text)
    rel = must(r'\(the module ships REL = ([\d.]+)\)', text)[0]
    agree = must(r'^\s+(E0, the published framework|new check, correct-or-not|new check, three-way)\s+([\d.]+)\s*$', text)
    return {'shipped_rel': float(rel),
            'split_half': [{'fit_half': int(a), 'fitted_rel': float(b), 'held_out_half': int(c), 'held_out_agreement': float(d)}
                           for a, b, c, d in halves],
            'agreement': {k: float(v) for k, v in agree}}


def framework_constants() -> dict:
    src = (REPO / 'engtrace_evaluation_framework.py').read_text(encoding='utf-8', errors='replace')
    out = {}
    for name in ('STEP_TOLERANCE', 'SEMANTIC_TAU', 'ERROR_SAMPLE_RATE', 'TRIBUNAL_TRIGGER_THRESHOLD', 'FINAL_TOLERANCE'):
        m = re.search(rf'^{name}\s*=\s*([\d.]+)', src, re.M)
        out[name] = float(m.group(1)) if m else None
    return out


def stack_constants() -> dict:
    ev = PILOT / 'evaluators'
    a = re.search(r'^REL\s*=\s*([\d.]+)', (ev / 'answer.py').read_text(encoding='utf-8'), re.M)
    d = re.search(r'^DISPLAY_TOL\s*=\s*([\d.]+)', (ev / 'milestones.py').read_text(encoding='utf-8'), re.M)
    return {'answer_rel': float(a.group(1)) if a else None, 'milestone_tol': float(d.group(1)) if d else None}


def full_run() -> dict:
    res = json.loads((HERE / 'results' / 'results.json').read_text(encoding='utf-8'))
    gold = json.loads((HERE / 'gold_validation.json').read_text(encoding='utf-8'))
    sens = res['sensitivity']
    q3 = {r['model']: r for r in res['q3']}
    return {'sensitivity': sens['models'], 'tau_with_headline': sens['tau_with_headline'],
            'digit_vs_tol1_on_fully_solved': {m: {'digit_rule': r['digit_flag_rate_on_fully_solved'],
                                                  'tol1': r['tol1_flag_rate_on_fully_solved']} for m, r in q3.items()},
            'e3_floor_readable_wrong': {m: r['e3_null_on_readable_wrong'] for m, r in q3.items()},
            'gold_null': next((v for k, v in gold.items() if 'null' in k.lower()), None),
            'gold_null_key': next((k for k in gold if 'null' in k.lower()), None),
            'store_commit': str(res['provenance']['store'].get('git', '?'))[:7]}


def render(d: dict) -> str:
    def f3(v):
        return '-' if v is None else f'{v:.3f}'
    sc, fw, fr = d['stack_constants'], d['framework_constants'], d['full_run']
    L = ['# The thresholds in the stack, and what each was measured under', '',
         'Generated by `threshold_appendix.py`; what it collates and how is in its docstring. Every row is parsed from a '
         'committed script\'s printed output or read from `results/results.json`; nothing is typed in. The pilot figures are '
         'on the 300 expert-labelled traces (15 templates, five models), from the pilot\'s scripts run again on 2026-10-03 with '
         'the evaluators as they now stand (the answer check after D-137 to D-169, the digit rule after D-156 to D-160); where '
         'a figure differs from the pilot\'s RESULTS files, those were computed with the code of their date, and '
         '`SCORER_VALIDATION.md` reproduces the published ones with the published code. The full-run figures are on the main '
         f"store at `{fr['store_commit']}`.", '',
         '## 1. The answer check: one relative tolerance', '',
         f"The check accepts a value within a relative tolerance of {sc['answer_rel']} of the gold's, or within one unit of the "
         'last digit shown, or of the gold\'s (`answer.match`); the inclusive boundary is D-147. The tolerance is the check\'s '
         'one fitted number, chosen on one half of the pilot\'s traces and reported on the other:', '',
         '| fitted on half | fitted tolerance | held out on half | three-way agreement with the experts |', '|---:|---:|---:|---:|']
    for h in d['answer_check']['split_half']:
        L.append(f"| {h['fit_half']} | {h['fitted_rel']:.4f} | {h['held_out_half']} | {h['held_out_agreement']:.3f} |")
    ag = d['answer_check']['agreement']
    L += ['', f"Agreement with the experts on the 300 traces: the published check {f3(ag.get('E0, the published framework'))}, "
          f"the corrected check {f3(ag.get('new check, correct-or-not'))} correct-or-not and {f3(ag.get('new check, three-way'))} "
          'three-way (RESULTS_X1 Finding 1b).', '',
          'On the full run, the answer score under the tolerance halved and doubled, the fully-solved rate, the half-unit window '
          'and the whole-trace reading (the Sensitivity table of RESULTS.md), with Kendall\'s tau of each ordering against the '
          'headline:', '',
          '| model | half | fitted | double | half-unit window | whole trace |', '|---|---:|---:|---:|---:|---:|']
    for r in fr['sensitivity']:
        L.append(f"| `{r['model']}` | {f3(r['half_tol'])} | {r['fitted']:.3f} | {f3(r['double_tol'])} | {f3(r['half_unit'])} | {f3(r['whole_trace'])} |")
    t = fr['tau_with_headline']
    L += [f"| tau with the headline | {f3(t['half_tol'])} | 1 | {f3(t['double_tol'])} | {f3(t['half_unit'])} | {f3(t['whole_trace'])} |", '',
          '## 2. E3: the milestone tolerance and unit scaling', '',
          f"A trace reaches a milestone when it states the value within {sc['milestone_tol']} relative (the template's display tolerance), "
          'under exact units or one of the unit factors. The grid below is the pilot\'s: real coverage over the 300 traces, the '
          'null (the same traces against a sibling item\'s milestones, shared values removed), their separation, and the '
          'correlation with the framework\'s answer check, at four tolerances with and without unit scaling (`analysis/e3_grid.py`, '
          'templates pinned at the freeze). The 0.5% setting sits on a plateau.', '',
          '| tolerance | unit scaling | real | null | separation | corr. with FAC |', '|---:|---|---:|---:|---:|---:|']
    for r in d['e3_grid']:
        L.append(f"| {r['tolerance'] * 100:g}% | {'yes' if r['unit_scaling'] else 'no'} | {r['real']:.3f} | {r['null']:.3f} | {r['separation']:.3f} | {r['corr_with_fac']:.3f} |")
    L += ['', f"On the full run the null on the gold solutions themselves is {f3(fr['gold_null'])} (`gold_validation.json`, `{fr['gold_null_key']}`), "
          'and per model on the readable wrong answers: ' + ', '.join(f"`{m}` {f3(v)}" for m, v in fr['e3_floor_readable_wrong'].items()) + '.', '',
          '## 3. The digit rule: four readings of the arithmetic checker', '',
          'On the pilot\'s labels, inside correct-answer traces (the hard case) and over all 300: the 1% relative tolerance E4 '
          'shipped with, 0.1%, the bare digit rule (the displayed precision alone) and the rule as the evaluator now ships it '
          f"(`analysis/digit_rule.py`; {d['digit_rule']['traces_scored']} traces scored, {d['digit_rule']['traces_skipped']} skipped for a step-count mismatch).", '',
          '| reading | hard case: tp / fp / fn | precision | recall | F1 | all traces: precision | recall | F1 | trace AUROC (95% CI) |',
          '|---|---|---:|---:|---:|---:|---:|---:|---|']
    names = {'tol1': '1% tolerance', 'tol01': '0.1% tolerance', 'digit': 'digit rule, bare', 'e4': 'digit rule, as shipped'}
    for k in ('tol1', 'tol01', 'digit', 'e4'):
        h, a, tl = d['digit_rule']['correct_answer'][k], d['digit_rule']['all'][k], d['digit_rule']['trace_level'][k]
        L.append(f"| {names[k]} | {h['tp']} / {h['fp']} / {h['fn']} | {h['precision']:.3f} | {h['recall']:.3f} | {h['f1']:.3f} | "
                 f"{a['precision']:.3f} | {a['recall']:.3f} | {a['f1']:.3f} | {tl['auroc']:.3f} ({tl['ci'][0]:.3f} to {tl['ci'][1]:.3f}) |")
    L += ['', 'On the full run, the share of fully solved traces with a flag under the shipped rule and at the 1% reading:', '',
          '| model | digit rule | at 1% |', '|---|---:|---:|']
    for m, v in fr['digit_vs_tol1_on_fully_solved'].items():
        L.append(f"| `{m}` | {v['digit_rule']:.3f} | {v['tol1']:.3f} |")
    L += ['', '## 4. E2: the process reward models\' 0.5 cut-off', '',
          'The threshold that maximises step F1 against the experts\' labels, fitted on one half of the traces and reported on '
          'the other, against the stock 0.5 (`analysis/prm_threshold.py`; the gain is the mean over both directions of the split).', '',
          '| subset | PRM | F1 at the fitted threshold | F1 at 0.5 | gain | fitted threshold (half 0, half 1) |', '|---|---|---:|---:|---:|---|']
    for r in d['prm_threshold']:
        L.append(f"| {r['subset']} | {r['prm']} | {r['f1_at_fitted']:.3f} | {r['f1_at_0.5']:.3f} | {r['gain']:+.3f} | "
                 f"{r['fitted_threshold_halves'][0]:.3f}, {r['fitted_threshold_halves'][1]:.3f} |")
    L += ['', '## 5. The judged stages, and the published framework\'s constants', '',
          'E5\'s judge and the step router carry no numeric threshold: a milestone is REACHED, NOT_NEEDED or MISSING, and a step '
          'is flagged or not; the strict reading (REACHED only) was chosen after the judge gave NOT_NEEDED to a quarter of fabricated '
          'values (RESULTS_E5). For the record, the published framework\'s constants as its source states them: step tolerance '
          f"{fw['STEP_TOLERANCE']}, cross-encoder threshold {fw['SEMANTIC_TAU']}, tribunal trigger {fw['TRIBUNAL_TRIGGER_THRESHOLD']}, "
          f"wrong-answer sample rate {fw['ERROR_SAMPLE_RATE']}, final-answer tolerance {fw['FINAL_TOLERANCE']}; its cross-encoder and "
          'alignment-ratio thresholds no longer exist in the stack, so the sensitivity the July rebuttal promised for them is moot.', '',
          'What cannot be redone: the tolerance\'s split-half fit rests on the pilot\'s labels; the full run has none, so no threshold '
          'was re-tuned there (next steps, "Not to do").', '']
    return '\n'.join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--cached', action='store_true')
    a = ap.parse_args()
    d = {'stack_constants': stack_constants(), 'framework_constants': framework_constants(), 'full_run': full_run(),
         'answer_check': parse_answer(run('answer_check', a.cached)),
         'e3_grid': parse_e3(run('e3_grid', a.cached)),
         'digit_rule': parse_digit(run('digit_rule', a.cached)),
         'prm_threshold': parse_prm(run('prm_threshold', a.cached))}
    OUT_JSON.parent.mkdir(exist_ok=True)
    OUT_JSON.write_text(json.dumps(d, indent=1), encoding='utf-8')
    text = render(d)
    OUT_MD.write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
