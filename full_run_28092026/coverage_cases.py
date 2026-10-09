"""Two coverage cases from saved outputs, written as appendices/coverage_cases.tex (the supervisor's comment of
7 October on coverage and alternative derivations).

    python -m full_run_28092026.coverage_cases                  # FREE: print the block
    python -m full_run_28092026.coverage_cases --write          # FREE: write into overleaf_source_04102026/
    python -m full_run_28092026.coverage_cases --selftest       # FREE: recompute every printed value

CASE A, a valid route that bypasses milestones: a correct answer from the main store whose derivation never states
two of the three gold milestones, because it differentiates symbolically where the gold trace evaluates the two
partial derivatives at the point. Matching finds the one milestone the response states, the judge rules the other
two not needed, and MC credits only reached, so m(r) = 1/3 for a correct, valid derivation. Chosen from the
responses of the first five models with a correct answer, no milestone ruled missing and at least one ruled not
needed (1,279 on the final store), the shortest question plus response among those with three milestones.

CASE B, matched values behind a misstated rule: a planted conceptual defect of the validation study
(evaluator_pilot_17092026/scores/planted_judges/probe_set.json, record conc:T-4d854e). The step's rule text is
changed and no digit is: "R = A / P" becomes "R = P / A" while the numbers still compute A / P. Every milestone
matches, the arithmetic check passes, the final answer is correct, and the four one-step judge probes of the study
(replies.jsonl) show which judges call the step a conceptual error.

Every value printed is read from the store rows, the trace, the probe record and the judges' replies;
`--selftest` recomputes the scores and checks the quoted text against the trace.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import worked_example as we  # noqa: E402  (gather, tex, num, factor, wrap, NAME)
from full_run_28092026 import score  # noqa: E402

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

MODEL_A, ITEM_A = 'kimi-k3', 'vorticity_check#6'
PROBE_B = 'conc:T-4d854e'
PILOT = REPO / 'evaluator_pilot_17092026'
PROBES = PILOT / 'scores' / 'planted_judges' / 'probe_set.json'
REPLIES = PILOT / 'scores' / 'planted_judges' / 'replies.jsonl'
SLICE = PILOT / 'slice' / 'manifest.jsonl'
TEX_OUT = Path('appendices') / 'coverage_cases.tex'
LABEL = 'box:coverage_cases'
MARK = ('% BEGIN GENERATED {name} (full_run_28092026/coverage_cases.py --write)\n{body}\n'
        '% END GENERATED {name}')
JUDGE_NAME = {'gpt-5': 'GPT-5', 'opus-4.5': 'Claude Opus 4.5', 'mimo-v2.5-pro': 'MiMo-V2.5-Pro', 'grok-4.6': 'Grok 4.6'}
JUDGE_ORDER = ['gpt-5', 'opus-4.5', 'mimo-v2.5-pro', 'grok-4.6']
PAPER_JUDGE = 'mimo-v2.5-pro'
STUDY_NAME = {'claude-opus-4.7': 'Claude Opus 4.7', 'gpt-5': 'GPT-5', 'gemini-3.1-pro': 'Gemini 3.1 Pro',
              'deepseek-r1': 'DeepSeek-R1', 'llama-3.1-70b': 'Llama 3.1 70B'}


# ------------------------------------------------------------------ case A

def case_a() -> dict:
    d = we.gather(MODEL_A, ITEM_A)
    text = d['text']
    # the route: the response works symbolically and evaluates the vorticity at the point in one product
    symbolic = re.search(r'=\s*2x\s*-\s*\(-x\)\s*=\s*3x', text)
    evaluated = re.search(r'3\(1\.3\)\s*=\s*3\.9', text)
    assert symbolic and evaluated, 'the route the box describes is not in the trace'
    not_needed = [s for s in d['settled'] if s['source'] == 'NOT_NEEDED']
    matched = [s for s in d['settled'] if s['source'] == 'e3']
    assert d['label'] == 'correct' and d['v'] == 1.0 and len(not_needed) == 2 and len(matched) == 1 and d['n'] == 3, d['settled']
    assert not any(s['source'] in ('MISSING', 'UNJUDGED') for s in d['settled'])
    route_adjusted = d['reached'] / (d['n'] - len(not_needed))
    return {**d, 'not_needed': not_needed, 'matched': matched, 'route_adjusted': route_adjusted}


# ------------------------------------------------------------------ case B

def case_b() -> dict:
    probes = json.loads(PROBES.read_text(encoding='utf-8'))
    p = next(x for x in probes if x['set_id'] == PROBE_B and x['arm'] == 'planted')
    assert p['family'] == 'conceptual' and p['step'].count(p['after']) == 1, p
    untouched = p['step'].replace(p['after'], p['before'], 1)
    digits = lambda s: re.findall(r'\d', s)
    assert digits(untouched) == digits(p['step']), 'a conceptual plant changed a digit'
    # the arithmetic check passes the planted step: its numbers still compute the right quantity
    claims = score.arith.check(p['step']).claims
    assert claims and all(c.ok_digit for c in claims), [c.__dict__ for c in claims]
    verdicts = {}
    for line in REPLIES.read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get('set_id') != PROBE_B or r.get('arm') != 'planted':
            continue
        try:
            cats = [x['category'] for x in json.loads(r['text'])['results']]
        except (KeyError, ValueError, TypeError):
            cats = None
        verdicts[r['judge']] = cats
    question = next((json.loads(l)['question'] for l in SLICE.read_text(encoding='utf-8').splitlines()
                     if l.strip() and json.loads(l).get('item_id') == p['item_id']), None)
    assert question, f"{p['item_id']} is not in the validation study's manifest"
    clean = lambda s: re.sub(r'\s+', ' ', s.replace('**', '')).strip()
    return {'probe': p, 'untouched': clean(untouched), 'planted': clean(p['step']), 'verdicts': verdicts,
            'question': question, 'claim': claims[0]}


# ------------------------------------------------------------------ LaTeX

def render(a: dict, b: dict) -> str:
    tex, num, wrap = we.tex, we.num, we.wrap
    L = [f"% Case A: {a['model']} on {a['item_id']} (scores/main); case B: probe {PROBE_B} of the",
         '% validation study\'s planted defects; generated by full_run_28092026/coverage_cases.py from the',
         '% store rows, the trace, the probe record and the judges\' replies.',
         r'\begin{figure}[t]', r'\centering',
         r'\begin{tcolorbox}[colback=gray!5, colframe=gray!40, boxrule=0.4pt, arc=2pt, left=5pt, right=5pt,',
         r'    top=4pt, bottom=4pt]', r'\footnotesize', r'\setlength{\tabcolsep}{3pt}',
         r'\renewcommand{\arraystretch}{1.05}']
    # --- A
    item = a['item']
    L.append(wrap(r'\textbf{A. A valid route that bypasses milestones.} ' + f"A {we.BRANCH[item['branch']]}-engineering "
                  f"question ({item['level']}) gives a velocity field and asks for the vorticity at a point:\\\\"))
    qlines = [l.strip() for l in item['question'].splitlines() if l.strip()]
    for i, q in enumerate(qlines):
        L.append(wrap(tex(q) + (r'\\' if i < len(qlines) - 1 else '')))
    L += ['', r'\nopagebreak\smallskip', r'\noindent\begin{tabular}{@{}l r l@{}}', r'\toprule',
          r'\textbf{Milestone} & \textbf{Gold value} & \textbf{Settled by} \\', r'\midrule']
    for s in a['settled']:
        c = '1' if abs(s.get('c', 1.0) - 1.0) < 1e-12 else we.factor(s['c'])
        how = (f"matching: states ${num(s['stated'])}$, $c = {c}$" if s['source'] == 'e3'
               else we.RULING.get(s['source'], s['source'].lower()))
        L.append(f"\\texttt{{{tex(s['id'])}}} & ${num(s['value'])}$ & {how} \\\\")
    L += [r'\bottomrule', r'\end{tabular}', '', r'\smallskip']
    nn = ' and '.join(f"\\texttt{{{tex(s['id'])}}}" for s in a['not_needed'])
    L.append(wrap(r'\noindent ' + f"The gold trace evaluates the two partial derivatives at the point; the response "
                  f"differentiates symbolically, $\\omega_z = 2x - (-x) = 3x$, and evaluates $3(1.3) = 3.9$ in one "
                  f"step, so it never states {nn}. The answer is correct, $v(r) = {num(a['v'])}$; MC credits only "
                  f"reached, so $m(r) = {a['reached']}/{a['n']}$, and route-adjusted coverage, with the not-needed "
                  f"milestones out of the denominator, is {a['reached']}/{a['n'] - len(a['not_needed'])}."))
    # --- B
    p, cl = b['probe'], b['claim']
    L += ['', r'\medskip', r'\hrule', r'\smallskip']
    L.append(wrap(r'\noindent\textbf{B. Matched values behind a misstated rule.} ' + f"A planted defect of the validation "
                  f"study: a clean response by \\texttt{{{STUDY_NAME.get(p['model_key'], p['model_key'])}}} to a "
                  f"rectangular-channel discharge question (Manning's equation) with one step's rule text changed and "
                  f"no digit changed.\\\\"))
    L.append(wrap(r'\textit{Untouched step:} ' + tex(b['untouched']) + r'\\'))
    L.append(wrap(r'\textit{Planted step:} ' + tex(b['planted'])))
    L += ['', r'\smallskip']
    judged = [k for k in JUDGE_ORDER if b['verdicts'].get(k)]
    caught = [k for k in judged if 'Conceptual Error' in b['verdicts'][k]]
    missed = [k for k in judged if k not in caught]
    names = lambda ks: ' and '.join(f"\\texttt{{{JUDGE_NAME[k]}}}" for k in ks)
    paper_judge_note = (' (the judge of this paper)' if PAPER_JUDGE in missed else '')
    unit = we.unit_tex(cl.right_unit)
    right = cl.right.strip()
    if cl.right_unit and right.endswith(cl.right_unit):
        right = right[:-len(cl.right_unit)].strip()
    L.append(wrap(r'\noindent ' + f"The hydraulic radius is $A/P$; the planted step states $P/A$ while its numbers "
                  f"still compute $A/P$ (${tex(cl.left.strip())} = {tex(right)}$" + (f"\\,{unit}" if unit else '')
                  + " passes the arithmetic check), every milestone matches, and the final answer is correct, so no "
                  f"deterministic check flags it. Shown the step alone, {names(caught)} "
                  f"{'labels' if len(caught) == 1 else 'label'} it a conceptual error and {names(missed)}"
                  f"{paper_judge_note} {'does' if len(missed) == 1 else 'do'} not."))
    L += [r'\end{tcolorbox}']
    caption = (r"\textbf{What coverage does and does not measure.} A: a correct, valid derivation that bypasses two of "
               f"three gold milestones scores $m(r) = {a['reached']}/{a['n']}$, since Milestone Coverage credits the "
               "gold route. B: a misstated rule whose numbers are those of the right rule passes matching, the "
               "arithmetic check and the answer check; only a judge can see it, and not every judge does.")
    L.append(wrap(caption, indent=' ' * 9).replace(' ' * 9, '\\caption{', 1) + '}')
    L += [f'\\label{{{LABEL}}}', r'\end{figure}']
    return MARK.format(name=LABEL, body='\n'.join(L))


# ------------------------------------------------------------------ main

def selftest() -> int:
    a, b = case_a(), case_b()
    body = render(a, b)
    assert all(len(l) <= 100 for l in body.splitlines() if not l.startswith('%')), [l for l in body.splitlines() if len(l) > 100]
    assert ITEM_A not in body.split('\\begin{figure}')[1] and MODEL_A not in body.split('\\begin{figure}')[1]
    print(f"SELFTEST OK: case A {MODEL_A} on {ITEM_A}: label {a['label']}, sources "
          f"{[s['source'] for s in a['settled']]}, m(r) = {a['reached']}/{a['n']}, route-adjusted {a['route_adjusted']:.2f}; "
          f"case B {PROBE_B}: digits unchanged, arithmetic check passes, judges {b['verdicts']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--write', action='store_true', help='write into overleaf_source_04102026/appendices/')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    body = render(case_a(), case_b()) + '\n'
    if args.write:
        target = we.SRC / TEX_OUT
        target.write_text(body, encoding='utf-8', newline='\n')
        print(f'wrote {target} ({len(body.splitlines())} lines)')
    else:
        print(body)
    return 0


if __name__ == '__main__':
    sys.exit(main())
