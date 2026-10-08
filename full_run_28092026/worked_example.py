"""A worked scoring example from the main store, written as appendices/worked_example.tex (WS-C5, step C9).

    python -m full_run_28092026.worked_example --out DIR        # FREE: write DIR/appendices/worked_example.tex
    python -m full_run_28092026.worked_example                  # FREE: print the block
    python -m full_run_28092026.worked_example --write          # FREE: write into overleaf_source_04102026/ (Phase 2)
    python -m full_run_28092026.worked_example --candidates     # FREE: the responses that show every mechanism
    python -m full_run_28092026.worked_example --reader-gap     # FREE: size of the reader gap noted below, per model
    python -m full_run_28092026.worked_example --selftest       # FREE: recompute every printed value from the trace

THE RESPONSE. One real response from scores/main, chosen so that one box shows every mechanism of the scoring
process: a wrong final answer by a lower-tier model, at least one milestone matched deterministically under a
unit factor, one judge ruling of each kind (reached, not needed, missing), and an arithmetic flag. `--candidates`
lists every main-store response that meets those conditions (14 on the store of 2026-10-07), shortest first.
The default, MODEL on ITEM below, is the one among them whose rulings all hold up on reading:
  - matching settles `MT` under the factor 10^3 (the response works in kg/mol where the correlation takes g/mol,
    so it states M*T a thousand times smaller than the gold);
  - the judge rules `T_star` and `epsilon` not needed (the question gives the collision integral), `sigma2` and
    `denominator` reached (the response states them in m^2, a factor of 10^-20 the fixed unit factors do not
    hold), and `numerator` and `viscosity` missing (both wrong, by the same unit slip);
  - the arithmetic check flags the final division, whose printed result is 10^18 times smaller than its printed
    operands give.
The sibling response on `gas_viscosity_kinetic_theory#7` was passed over: its second deterministic match is an
artefact of the milestone reader, which drops a unicode superscript exponent (`1.152 × 10⁻¹⁹` reads as 1.152)
and so matches the numerator under the 3,600 factor; on `#13` the judge credits a wrong numerator as reached.
`--reader-gap` sizes that artefact over the main run (answer.values and arith.normalise fold such exponents;
milestones.numbers does not): a fix belongs to the evaluator code and a re-score, not to this script.

WHAT IS PRINTED, every value from the store rows and the trace: the question as the model received it
(line breaks kept, nothing cut), the gold milestones with their values, how each was settled (matched: the
stated value and the unit factor c of eq. match; else the judge's ruling), the flagged calculation with its
printed operands, its printed result, the value recomputed from the operands and the place value that judged
it, the stated final answer against its target, then v(r) and m(r). The item id and model key are in a LaTeX
comment, not in the box. The block is a `figure` float holding a `tcolorbox`, labelled `box:worked_example`.

SOURCES. scores/main/<model>.jsonl (answer, e3, steps, score), scores/main/e5/<model>.jsonl (sources,
e5_strict), scores/main/router/<model>.jsonl (digit_flagged), scores/milestones.json, traces/<model>.jsonl
(the text), pool/ (the question). `--selftest` re-runs the answer check, E3 and the arithmetic check on the
trace through score.py's own functions and asserts that they reproduce the stored label, matches, flags,
score and e5_strict, and that the printed v(r) and m(r) equal `score` and `e5_strict`.
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

from full_run_28092026 import score  # noqa: E402  (puts the evaluators on sys.path)

import milestones as ms  # noqa: E402

answer, arith, e2_prm, e3 = score.answer, score.arith, score.e2_prm, score.e3

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

MODEL = 'qwen3-235b-a22b-2507'
ITEM = 'gas_viscosity_kinetic_theory#0'
STORE = score.SCORES / 'main'
TEX_OUT = Path('appendices') / 'worked_example.tex'
SRC = REPO / 'overleaf_source_04102026'
LABEL = 'box:worked_example'
MARK = ('% BEGIN GENERATED {name} (full_run_28092026/worked_example.py --write)\n{body}\n'
        '% END GENERATED {name}')
NAME = {  # display names as paper_results.NAME spells them (copied; that module computes every table on import)
    'deepseek-v4.1-flash': 'DeepSeek V4.1 Flash', 'kimi-k3': 'Kimi K3', 'claude-sonnet-5': 'Claude Sonnet 5',
    'glm-5.3-flash': 'GLM-5.3-Flash', 'muse-glimmer-30b': 'Muse Glimmer 30B', 'glm-5.3': 'GLM-5.3',
    'qwen3-235b-a22b-2507': 'Qwen3-235B-2507', 'gemini-3.1-flash-lite': 'Gemini 3.1 Flash-Lite',
    'gemma-4-26b-a4b': 'Gemma 4 26B', 'gpt-5.4-mini': 'GPT-5.4 mini', 'gpt-oss-20b': 'gpt-oss-20b',
}
BRANCH = {'chemical_engineering': 'chemical', 'civil_engineering': 'civil', 'electrical_engineering': 'electrical',
          'industrial_engineering': 'industrial', 'mechanical_engineering': 'mechanical'}
RULING = {'REACHED': 'judge: reached', 'NOT_NEEDED': 'judge: not needed', 'MISSING': 'judge: missing',
          'UNJUDGED': 'judge: no ruling'}


# ------------------------------------------------------------------ the rows

def find(path: Path, item_id: str) -> dict | None:
    with path.open(encoding='utf-8') as fh:
        for line in fh:
            if line.strip() and json.loads(line)['item_id'] == item_id:
                return json.loads(line)
    return None


def rows(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {r['item_id']: r for r in (json.loads(l) for l in path.read_text(encoding='utf-8').splitlines() if l.strip())}


def gather(model: str, item_id: str) -> dict:
    sc = find(STORE / f'{model}.jsonl', item_id)
    e5 = find(STORE / 'e5' / f'{model}.jsonl', item_id)
    rt = find(STORE / 'router' / f'{model}.jsonl', item_id)
    tr = find(score.TRACES / f'{model}.jsonl', item_id)
    if not (sc and e5 and rt and tr):
        raise SystemExit(f'{model} / {item_id}: a store row is missing (score {bool(sc)}, e5 {bool(e5)}, '
                         f'router {bool(rt)}, trace {bool(tr)})')
    milestones = json.loads(score.MILESTONES.read_text(encoding='utf-8'))[item_id]
    item = score.pool_items()[item_id]
    text = tr['text']
    nums = ms.numbers(text)
    settled = []
    for m, src, scale in zip(milestones, e5['sources'], sc['e3']['scale']):
        entry = {'id': m['id'], 'value': m['value'], 'source': src}
        if src == 'e3':
            stated = next((n for n in nums if ms.close(n, m['value'] * scale, e3.STEP_TOL)), None)
            entry.update({'stated': stated, 'scale': scale, 'c': 1.0 / scale})   # eq. match: c * stated = gold
        settled.append(entry)
    steps = e2_prm.steps_of(text)
    flags = []
    for k in rt['digit_flagged']:
        for c in arith.check(steps[k]).claims:
            if not c.ok_digit:
                flags.append({'step': k, 'left': c.left, 'right': c.right, 'left_value': c.left_value[0],
                              'right_value': c.right_value[0], 'ulp': c.ulp, 'unit': c.right_unit})
    seg = answer.segment(text)
    stated = [v for v, _u in answer.values(seg)]
    targets = sc['answer']['targets']['numbers']
    n = len(milestones)
    reached = sum(1 for s in e5['sources'] if s in ('e3', 'REACHED'))
    return {'model': model, 'item_id': item_id, 'item': item, 'text': text, 'score_row': sc, 'e5_row': e5,
            'router_row': rt, 'milestones': milestones, 'settled': settled, 'steps': steps, 'flags': flags,
            'stated': stated, 'targets': targets, 'label': sc['answer']['label'], 'v': sc['score'],
            'reached': reached, 'n': n, 'm': reached / n if n else None,
            'claims': sum(s['claims'] for s in sc['steps'])}


# ------------------------------------------------------------------ LaTeX

UNICODE = {'μ': r'$\mu$', 'σ': r'$\sigma$', 'Ω': r'$\Omega$', 'ε': r'$\varepsilon$', 'ρ': r'$\rho$', 'τ': r'$\tau$',
           'θ': r'$\theta$', 'ω': r'$\omega$', 'π': r'$\pi$', 'λ': r'$\lambda$', 'γ': r'$\gamma$', 'Δ': r'$\Delta$',
           'α': r'$\alpha$', 'β': r'$\beta$', 'φ': r'$\phi$', 'η': r'$\eta$', 'ν': r'$\nu$', 'Å': r'\AA{}',
           '·': r'$\cdot$', '×': r'$\times$', '≈': r'$\approx$', '≤': r'$\le$', '≥': r'$\ge$', '−': '--', '–': '--',
           '—': '---', '°': r'$^\circ$', '√': r'$\surd$', '∞': r'$\infty$', '…': r'\ldots{}', '’': "'",
           '“': '``', '”': "''", ' ': '~'}
SUB = {c: str(i) for i, c in enumerate('₀₁₂₃₄₅₆₇₈₉')}
SUP = {**{c: str(i) for i, c in enumerate('⁰¹²³⁴⁵⁶⁷⁸⁹')}, '⁻': '-', '⁺': '+'}


TOKEN = re.compile('[₀-₉]+|[⁰¹²³⁴-⁹⁻⁺]+|.', re.S)


def tex(s: str) -> str:
    """Plain text for LaTeX: the special characters escaped, the unicode the traces use rendered (a run of
    subscript or superscript digits becomes one math group)."""
    out = []
    for tok in TOKEN.findall(s):
        if tok[0] in SUB:
            out.append('$_{' + ''.join(SUB[c] for c in tok) + '}$')
        elif tok[0] in SUP:
            out.append('$^{' + ''.join(SUP[c] for c in tok) + '}$')
        elif tok in UNICODE:
            out.append(UNICODE[tok])
        elif tok in '&%$#_{}':
            out.append('\\' + tok)
        elif tok == '\\':
            out.append(r'\textbackslash{}')
        elif tok == '^':
            out.append(r'\^{}')
        elif tok == '~':
            out.append(r'\textasciitilde{}')
        elif ord(tok) > 127:
            print(f'WARNING: no LaTeX rendering for U+{ord(tok):04X} {tok!r}; printed as ?', file=sys.stderr)
            out.append('?')
        else:
            out.append(tok)
    return ''.join(out)


def num(v: float, sig: int = 6) -> str:
    """A value in math mode as the judge prints it (%.6g): scientific notation as a power of ten, thousands
    grouped with a comma (braced, so math mode sets no space after it)."""
    s = f'{v:.{sig}g}'
    if 'e' in s:
        mant, exp = s.split('e')
        return f'{mant}\\times10^{{{int(exp)}}}'
    if abs(v) >= 1e4:
        whole, _, frac = s.partition('.')
        s = f'{int(whole):,}'.replace(',', '{,}') + ('.' + frac if frac else '')
    return s


def factor(c: float) -> str:
    e = round(_log10(c))
    return f'10^{{{e}}}' if abs(c - 10 ** e) < 1e-9 * c else num(c, 4)


def _log10(x: float) -> float:
    import math
    return math.log10(x)


def claim_tex(expr: str, unit: str = '') -> str:
    """arith's normalised text of one side of a claim, in math mode, its unit tail removed:
    `((2.1155e-5)) / ((3.128e-19))` or `(6.762e-5) Pa*s`."""
    s = expr.strip()
    if unit and s.endswith(unit):
        s = s[:-len(unit)].strip()
    s = re.sub(r'\(\(([^()]*)\)\)', r'(\1)', s)
    s = re.sub(r'(\d)e([-+]?\d+)', lambda m: f'{m.group(1)}\\times10^{{{int(m.group(2))}}}', s)
    s = s.replace('**', '^').replace('*', r'\cdot ')
    return s.strip()


def wrap(text: str, width: int = 100, indent: str = '') -> str:
    import textwrap
    return '\n'.join(textwrap.wrap(text, width=width, initial_indent=indent, subsequent_indent=indent,
                                   break_long_words=False, break_on_hyphens=False))


def unit_tex(u: str) -> str:
    return tex(u.replace('*', '·')) if u else ''


def render(d: dict) -> str:
    item, sc = d['item'], d['score_row']
    name = NAME.get(d['model'], d['model'])
    L = [f"% Worked example: {d['model']} on {d['item_id']} (scores/main), generated by",
         '% full_run_28092026/worked_example.py; every value is read from the store rows and the trace.',
         r'\begin{figure}[t]', r'\centering',
         r'\begin{tcolorbox}[colback=gray!5, colframe=gray!40, boxrule=0.4pt, arc=2pt, left=5pt, right=5pt,',
         r'    top=4pt, bottom=4pt]', r'\footnotesize', r'\setlength{\tabcolsep}{3pt}',
         r'\renewcommand{\arraystretch}{1.05}']
    # the question
    L.append(r'\textbf{Question} (as posed; ' + f"{BRANCH[item['branch']]} engineering, {item['level']}).\\\\")
    qlines = [l.strip() for l in item['question'].splitlines() if l.strip()]
    for i, q in enumerate(qlines):
        body = tex(q[2:]) if q.startswith('- ') else tex(q)
        lead = '--~' if q.startswith('- ') else ''
        L.append(wrap(lead + body + (r'\\' if i < len(qlines) - 1 else '')))
    # the milestones
    L += ['', r'\medskip', r'\noindent\textbf{Gold milestones} ($\mathcal{M}$, ' + f"{d['n']}" + ') and how each was settled.',
          '', r'\nopagebreak\smallskip', r'\noindent\begin{tabular}{@{}l r l@{}}', r'\toprule',
          r'\textbf{Milestone} & \textbf{Gold value} & \textbf{Settled by} \\', r'\midrule']
    for s in d['settled']:
        if s['source'] == 'e3':
            how = f"matching: states ${num(s['stated'])}$, $c = {factor(s['c'])}$"
        else:
            how = RULING.get(s['source'], s['source'].lower())
        L.append(f"\\texttt{{{tex(s['id'])}}} & ${num(s['value'])}$ & {how} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    # the arithmetic flag
    if d['flags']:
        f = d['flags'][0]
        unit = unit_tex(f['unit'])
        unit_txt = f'\\,{unit}' if unit else ''
        L += ['', r'\medskip']
        L.append(wrap(r'\noindent\textbf{Arithmetic check.} ' + f"The flagged step prints "
                      f"${claim_tex(f['left'], f['unit'])} = {claim_tex(f['right'], f['unit'])}${unit_txt}; "
                      f"the printed operands give ${num(f['left_value'], 4)}$, which does not round to the printed "
                      f"result at its last digit (${num(f['ulp'], 1)}$), so the step is flagged"
                      + (f" ({len(d['flags'])} flags in {len(d['steps'])} steps, {d['claims']} calculations parsed)."
                         if len(d['flags']) > 1 else
                         f" (the one flag in {len(d['steps'])} steps, {d['claims']} calculations parsed).")))
    # the final answer
    y = d['targets'][0]
    yhat = d['stated'][0] if d['stated'] else None
    gap = abs(yhat - y) / abs(y) if yhat is not None else None
    sentence = (r'\noindent\textbf{Final answer.} ' + f"The response states $\\hat{{y}} = {num(yhat)}$ against the "
                f"target $y = {num(y)}$: at $c = 1$ the gap is {100 * gap:.0f}\\% of $y$, and no other factor in "
                r'$\mathcal{C}$ brings it within $\epsilon$ or one unit of either last digit (\autoref{eq:match}), '
                f"so the answer is {tex(d['label'])} and $v(r) = {num(d['v'])}$.")
    L += ['', r'\medskip', wrap(sentence)]
    n_e3 = sum(1 for s in d['settled'] if s['source'] == 'e3')
    n_judge = sum(1 for s in d['settled'] if s['source'] == 'REACHED')
    L += ['', r'\smallskip']
    L.append(wrap(r'\noindent\textbf{Coverage.} ' + f"$|\\mathcal{{R}}(r)| = {n_e3}$ matched $+\\,{n_judge}$ ruled "
                  f"reached $= {d['reached']}$ of $|\\mathcal{{M}}| = {d['n']}$, so "
                  f"$m(r) = {d['reached']}/{d['n']} = {d['m']:.2f}$ (\\autoref{{eq:coverage}})."))
    L += [r'\end{tcolorbox}']
    caption = (f"A worked scoring example: a wrong answer by \\texttt{{{name}}} on an {item['level']} "
               f"{BRANCH[item['branch']]}-engineering instance. Matching settles {n_e3} of the {d['n']} milestones "
               f"under a unit factor, the judge rules the other {d['n'] - n_e3}, the arithmetic check flags "
               f"{len(d['flags'])} of the {len(d['steps'])} steps, and the response scores "
               f"$v(r) = {num(d['v'])}$ and $m(r) = {d['reached']}/{d['n']}$.")
    L.append(wrap(caption, indent=' ' * 9).replace(' ' * 9, '\\caption{', 1) + '}')
    L += [f'\\label{{{LABEL}}}', r'\end{figure}']
    return MARK.format(name=LABEL, body='\n'.join(L))


# ------------------------------------------------------------------ candidates

def candidates() -> list[dict]:
    """Every main-store response with a readable wrong answer, a deterministic match, all three judge rulings
    and an arithmetic flag, shortest (question plus response) first."""
    qlen = {}
    for p in sorted(score.POOL.rglob('*.jsonl')):
        for l in p.read_text(encoding='utf-8').splitlines():
            r = json.loads(l)
            qlen[r['item_id']] = len(r['question'])
    out = []
    for f in sorted(STORE.glob('*.jsonl')):
        model = f.stem
        e5, rt = rows(STORE / 'e5' / f.name), rows(STORE / 'router' / f.name)
        for iid, r in rows(f).items():
            a = r.get('answer') or {}
            e, t = e5.get(iid), rt.get(iid)
            if r.get('status') != 'answered' or not r.get('readable') or a.get('label') != 'incorrect' or not e or not t:
                continue
            kinds = set(e['sources'])
            flags = sum(s.get('digit_flags', 0) for s in r['steps'])
            if 'e3' not in kinds or flags < 1 or not {'REACHED', 'NOT_NEEDED', 'MISSING'} <= kinds or 'UNJUDGED' in kinds:
                continue
            out.append({'model': model, 'item_id': iid, 'milestones': len(e['sources']),
                        'matched': e['sources'].count('e3'),
                        'unit_factor': sum(1 for s in r['e3']['scale'] if s not in (None, 1.0)),
                        'flags': flags, 'steps': len(r['steps']), 'tokens': r.get('completion_tokens') or 0,
                        'question_chars': qlen.get(iid, 0), 'e5_strict': e['e5_strict'], 'level': r['level']})
    out.sort(key=lambda c: c['question_chars'] + 2 * c['tokens'])
    return out


# ------------------------------------------------------------------ the reader gap behind the passed-over response

SUP_DIGITS = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
SUP_POWER = re.compile('10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)')


def fold_superscripts(text: str) -> str:
    """`3.467 × 10⁻¹⁰` as `3.467 × 10^{-10}`, the form milestones.numbers folds into one number."""
    return SUP_POWER.sub(lambda m: '10^{' + m.group(1).translate(SUP_DIGITS) + '}', text)


def reader_gap() -> list[dict]:
    """Per model on the main traces: the answered responses that write a power of ten with unicode superscripts
    (which milestones.numbers reads as the mantissa and 10), and what E3 would do with the exponent folded in:
    milestones newly reached, matches lost (the mantissa no longer matches under a unit factor), responses whose
    E3 row changes, and judge jobs whose missed-milestone list, hence prompt, would change. A diagnostic only:
    nothing here re-scores a store."""
    MS = json.loads(score.MILESTONES.read_text(encoding='utf-8'))
    out = []
    for f in sorted(score.TRACES.glob('*.jsonl')):
        if not (STORE / f.name).exists():
            continue
        n = gain = loss = changed = jobs = 0
        with f.open(encoding='utf-8') as fh:
            for line in fh:
                r = json.loads(line)
                text = r.get('text') or ''
                if r.get('status') != 'answered' or not SUP_POWER.search(text):
                    continue
                n += 1
                mil = MS.get(r['item_id']) or []
                if not mil:
                    continue
                a, b = e3.reach(mil, text), e3.reach(mil, fold_superscripts(text))
                g = sum(1 for x, y in zip(a, b) if not x['reached'] and y['reached'])
                lo = sum(1 for x, y in zip(a, b) if x['reached'] and not y['reached'])
                gain, loss = gain + g, loss + lo
                if g or lo:
                    changed += 1
                    jobs += any(not x['reached'] for x in a)
        out.append({'model': f.stem, 'with_form': n, 'gain': gain, 'loss': loss, 'changed': changed, 'jobs': jobs})
    return out


# ------------------------------------------------------------------ self-test

def selftest(model: str = MODEL, item_id: str = ITEM) -> int:
    d = gather(model, item_id)
    item, sc, text = d['item'], d['score_row'], d['text']
    ans = score.score_answer(item, text, d['milestones'])
    assert ans['label'] == sc['answer']['label'] == d['label'], (ans['label'], sc['answer']['label'])
    assert ans['matched'] == sc['answer']['matched'] and ans['of'] == sc['answer']['of'], ans
    assert ans['targets']['numbers'] == sc['answer']['targets']['numbers'], ans['targets']
    v = 0.0 if not score.readable(text) else score.SCORE_OF[ans['label']]
    assert v == sc['score'] == d['v'], (v, sc['score'])
    hits = e3.reach(d['milestones'], text)
    assert [h['reached'] for h in hits] == sc['e3']['reached'], hits
    assert [h['scale'] for h in hits] == sc['e3']['scale'], hits
    for s in d['settled']:
        if s['source'] == 'e3':
            assert s['stated'] is not None and ms.close(s['stated'] * s['c'], s['value'], e3.STEP_TOL), s
    steps = score.score_steps(text)
    assert [s['digit_flags'] for s in steps] == [s['digit_flags'] for s in sc['steps']], steps
    assert sorted(d['router_row']['digit_flagged']) == [k for k, s in enumerate(steps) if s['digit_flags']], steps
    assert d['flags'], 'no flagged claim found on the flagged steps'
    e5 = d['e5_row']
    assert len(e5['sources']) == len(d['milestones']) == e5['milestones_required']
    assert abs(d['m'] - e5['e5_strict']) < 1e-12, (d['m'], e5['e5_strict'])
    assert d['stated'], 'no value read from the answer segment'
    assert not answer.match(d['targets'][0], answer.values(answer.segment(text))), 'the stated answer matches'
    body = render(d)
    for needle in (f"$v(r) = {num(d['v'])}$", f"$m(r) = {d['reached']}/{d['n']} = {d['m']:.2f}$",
                   f"\\label{{{LABEL}}}", f"% BEGIN GENERATED {LABEL}"):
        assert needle in body, needle
    assert item_id not in body.split('\\begin{figure}')[1] and model not in body.split('\\begin{figure}')[1], \
        'the item id or model key is inside the box'
    assert all(len(l) <= 100 for l in body.splitlines() if not l.startswith('%')), \
        [l for l in body.splitlines() if len(l) > 100]
    print(f"SELFTEST OK: {model} on {item_id}: label {ans['label']} (score {v}), E3 {sc['e3']['reached']} with scales "
          f"{sc['e3']['scale']}, digit flags {[s['digit_flags'] for s in steps]} and e5_strict {e5['e5_strict']:.4f} "
          f"reproduced from the trace; the box prints v(r) = {num(d['v'])} and m(r) = {d['reached']}/{d['n']}")
    return 0


# ------------------------------------------------------------------ main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--model', default=MODEL)
    ap.add_argument('--item', default=ITEM)
    ap.add_argument('--out', type=Path, help='write DIR/appendices/worked_example.tex (a copy of the tex tree)')
    ap.add_argument('--write', action='store_true', help='write into overleaf_source_04102026/appendices/ (Phase 2)')
    ap.add_argument('--candidates', action='store_true')
    ap.add_argument('--reader-gap', action='store_true',
                    help='size the milestone reader\'s unicode-superscript gap on E3 (diagnostic; no store is changed)')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest(a.model, a.item)
    if a.reader_gap:
        print('Answered main-run responses that write a power of ten with unicode superscripts, and what E3 would do '
              'with the exponent folded in (milestones.numbers reads `1.152 × 10⁻¹⁹` as 1.152 and 10):')
        print('| model | responses with the form | milestones newly reached | matches lost | responses whose E3 changes '
              '| judge jobs whose prompt would change |')
        print('|---|---:|---:|---:|---:|---:|')
        for g in reader_gap():
            print(f"| {g['model']} | {g['with_form']} | {g['gain']} | {g['loss']} | {g['changed']} | {g['jobs']} |")
        return 0
    if a.candidates:
        cs = candidates()
        print(f'{len(cs)} responses show every mechanism (readable wrong answer, a deterministic match, '
              f'reached, not needed and missing rulings, an arithmetic flag), shortest first:')
        print('| model | item | milestones | matched | with a unit factor | flags | steps | tokens | question chars '
              '| e5_strict | level |')
        print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|')
        for c in cs:
            print(f"| {c['model']} | {c['item_id']} | {c['milestones']} | {c['matched']} | {c['unit_factor']} | "
                  f"{c['flags']} | {c['steps']} | {c['tokens']} | {c['question_chars']} | {c['e5_strict']:.3f} | "
                  f"{c['level']} |")
        return 0
    body = render(gather(a.model, a.item)) + '\n'
    if a.out or a.write:
        target = (a.out if a.out else SRC) / TEX_OUT
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(body, encoding='utf-8', newline='\n')
        print(f'wrote {target} ({len(body.splitlines())} lines)')
    else:
        print(body)
    return 0


if __name__ == '__main__':
    sys.exit(main())
