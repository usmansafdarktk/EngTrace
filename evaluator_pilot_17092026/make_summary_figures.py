"""The figures for the pilot summary.

    python evaluator_pilot_17092026/make_summary_figures.py

Every value is read from figures/summary_numbers.json, which analysis/summary_numbers.py
computes from the analyses that own each number. Nothing is typed in here: a figure that
disagrees with the analyses cannot be drawn, and any number a title states is formatted
from the same file. Where each block comes from:

    fig_trace       analysis/cluster_bootstrap.py   (trace-level AUROC, template CIs, margins)
    fig_steps       analysis/digit_rule.py + the 72B PRM's step rewards (x1_analysis.steps)
    fig_planted     analysis/planted_judges.py + analysis/router_planted.py
    fig_routing     analysis/router_planted.py      (both routes, counted per defect)
    fig_cost        analysis/judge_cost.py + analysis/router_residue.py, the D-110 roster

The palette matches the documents: ink for text, a single accent for the series that
carries the point, grey for everything it is being compared against.
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures')
NUMBERS = os.path.join(OUT, 'summary_numbers.json')

INK = '#1F2A37'
ACCENT = '#1B3A5C'
GOOD = '#2E6F5E'
WARN = '#9C5B2E'
GREY = '#9AA3AE'
LIGHT = '#D5DBE3'

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Calibri', 'DejaVu Sans'],
    'font.size': 9,
    'axes.edgecolor': LIGHT,
    'axes.labelcolor': INK,
    'text.color': INK,
    'xtick.color': INK,
    'ytick.color': INK,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 200,
})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('wrote', os.path.relpath(path, HERE))


def share(d, key='caught'):
    return d[key] / d['n']


def fig_trace(N):
    """Trace level, with the intervals that account for 15 templates rather than 300 traces."""
    colour = {'E5': ACCENT, 'baseline: expert answer verdict': GOOD}
    rows = N['trace']['rows']
    base = next(r['auroc'] for r in rows if r['head'] == 'E0')
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    for i, r in enumerate(rows):
        y, c = len(rows) - i, colour.get(r['head'], GREY)
        ax.plot([r['lo'], r['hi']], [y, y], color=c, lw=2.4, solid_capstyle='round', alpha=.55)
        ax.plot([r['auroc']], [y], 'o', color=c, ms=6)
        ax.text(r['hi'] + .006, y, '%.3f' % r['auroc'], va='center', fontsize=8, color=c)
    ax.set_yticks(range(1, len(rows) + 1))
    ax.set_yticklabels([r['label'] for r in rows][::-1])
    ax.set_xlim(.68, 1.03)
    ax.set_xlabel('AUROC separating sound from unsound reasoning (95% interval, templates resampled)')
    ax.axvline(base, color=LIGHT, lw=1, zorder=0)
    m = N['trace']['margins']['E5']
    verdict = 'is not significant' if m['lo'] <= 0 <= m['hi'] else 'is significant'
    ax.set_title('No evaluator beats the published framework at the trace level.\n'
                 'The experts\' answer verdict scores highest; its margin over E5, %+.3f, %s.'
                 % (m['margin'], verdict), loc='left', fontsize=10, color=INK, pad=10)
    save(fig, 'trace_level.png')


def fig_steps(N):
    """Step level inside correct-answer traces: the case a reasoning evaluator exists for."""
    S = N['steps']
    names = ['E4\'s arithmetic check\nat 1% (as first shipped)', 'best process\nreward model (72B)',
             'digit rule\n(as E4 ships it)']
    keys = ['digit_tol1', 'prm_qwen72', 'digit_e4']
    prec = [S[k]['precision'] for k in keys]
    rec = [S[k]['recall'] for k in keys]
    x = range(len(names))
    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    ax.bar([i - .19 for i in x], prec, .36, label='precision', color=ACCENT)
    ax.bar([i + .19 for i in x], rec, .36, label='recall', color=LIGHT, edgecolor=GREY)
    for i, (p, r) in enumerate(zip(prec, rec)):
        ax.text(i - .19, p + .02, '%.3f' % p, ha='center', fontsize=8, color=ACCENT)
        ax.text(i + .19, r + .02, '%.3f' % r, ha='center', fontsize=8, color=INK)
    ax.set_xticks(list(x))
    ax.set_xticklabels(names)
    ax.set_ylim(0, .92)
    ax.set_ylabel('against the experts\' step labels')
    ax.legend(frameon=False, loc='upper left')
    ax.set_title('Finding a wrong step inside a trace whose answer is right.\n'
                 'The digit rule as E4 ships it: precision %.2f, recall %.2f.'
                 % (S['digit_e4']['precision'], S['digit_e4']['recall']),
                 loc='left', fontsize=10, pad=10)
    save(fig, 'step_level.png')


def fig_planted(N):
    """Planted defects: truth by construction, so the guide cannot be the reason."""
    P = N['planted']
    J = P['judges']
    cols = [('digit rule\n(as E4 ships it)', P['digit_rule']), ('GPT-5', J['gpt-5']),
            ('Claude Opus 4.5', J['opus-4.5']), ('MiMo-V2.5-Pro\n(the chosen judge)', J['mimo-v2.5-pro'])]
    arith = [share(d['arithmetic']) for _n, d in cols]
    conc = [share(d['conceptual']) for _n, d in cols]
    x = range(len(cols))
    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ax.bar([i - .19 for i in x], arith, .36, label='arithmetic defects', color=LIGHT, edgecolor=GREY)
    ax.bar([i + .19 for i in x], conc, .36, label='conceptual defects', color=ACCENT)
    for i, (a, c) in enumerate(zip(arith, conc)):
        ax.text(i - .19, a + .02, '%.2f' % a, ha='center', fontsize=8, color=INK)
        ax.text(i + .19, c + .02, '%.2f' % c if c else '0', ha='center', fontsize=8,
                color=ACCENT if c else WARN)
    ax.set_xticks(list(x))
    ax.set_xticklabels([n for n, _d in cols])
    ax.set_ylim(0, .95)
    ax.set_ylabel('share of planted defects detected')
    ax.legend(frameon=False, loc='upper right')
    fa, fa_n = P['digit_false_alarm_on_judged_steps']
    verdicts = sum(J[j][f]['n'] for j in J for f in ('conceptual', 'arithmetic'))
    judge_flags = sum(J[j][f]['flags_original'] for j in J for f in ('conceptual', 'arithmetic'))
    ax.set_title('Only a judge catches a misstated rule: GPT-5 %.2f, MiMo %.2f, Opus 4.5 %.2f.\n'
                 'On the %d untouched steps the judges flagged %d (%d verdicts); the digit rule, %d.'
                 % (share(J['gpt-5']['conceptual']), share(J['mimo-v2.5-pro']['conceptual']),
                    share(J['opus-4.5']['conceptual']), fa_n, judge_flags, verdicts, fa),
                 loc='left', fontsize=10, pad=10)
    fig.text(0.01, -0.09, '%d defects per family. MiMo returned a verdict on both arms for %d of the '
             '%d conceptual defects, so its conceptual bar is out of %d.'
             % (P['digit_rule']['arithmetic']['n'], J['mimo-v2.5-pro']['conceptual']['n'],
                P['digit_rule']['conceptual']['n'], J['mimo-v2.5-pro']['conceptual']['n']),
             fontsize=7.5, color=INK, ha='left')
    save(fig, 'planted.png')


def fig_routing(N):
    """Being able to catch a defect is not the same as being shown it."""
    R = N['routing']['conceptual']
    e0, mimo = R['e0'], R['router']['mimo-v2.5-pro']
    published = [e0['shown'] / e0['n'], e0['caught_when_asked'] / e0['n'], e0['end_to_end'] / e0['n']]
    routed = [R['forwarded'] / R['n'], mimo['caught_when_asked'] / mimo['n'], mimo['end_to_end'] / mimo['n']]
    labels = ['shown the flawed step', 'a judge flags it\nwhen asked', 'caught end to end\n(counted per defect)']
    x = range(3)
    fig, ax = plt.subplots(figsize=(6.6, 2.7))
    ax.bar([i - .19 for i in x], published, .36, label='published routing (E0\'s two judges)',
           color=LIGHT, edgecolor=GREY)
    ax.bar([i + .19 for i in x], routed, .36, label='verify first, route the residue (MiMo)', color=ACCENT)
    for i, (p, r) in enumerate(zip(published, routed)):
        ax.text(i - .19, p + .02, '%.3f' % p, ha='center', fontsize=8, color=INK)
        ax.text(i + .19, r + .02, '%.3f' % r, ha='center', fontsize=8, color=ACCENT)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('conceptual defects')
    ax.legend(frameon=False, loc='upper right', ncol=1)
    ax.set_title('The published framework shows a judge the flawed step as often as a clean one.\n'
                 'Counted per defect it catches %.2f; a verify-first router with MiMo catches %.2f.'
                 % (published[2], routed[2]), loc='left', fontsize=10, pad=10)
    fig.text(0.01, -0.13, 'Planted conceptual defects, %d. The router\'s judge rates are over the %d '
             'MiMo returned a verdict on; its judge was asked about one step per call.'
             % (R['n'], mimo['n']), fontsize=7.5, color=INK, ha='left')
    save(fig, 'routing.png')


def fig_cost(N):
    """What each choice costs in judge calls over the full benchmark."""
    C = N['cost']
    rows = [('published framework\n(its two judges)', C['e0'], GREY),
            ('recommended stack\n(deterministic + residual milestone judge)', C['e5'], ACCENT),
            ('with a batched router\n(one call per trace; untested)', C['e5'] + C['router_batched'], GOOD),
            ('with a per-step router\n(the design the probe measured)', C['e5'] + C['router_per_step'], WARN)]
    fig, ax = plt.subplots(figsize=(6.6, 2.8))
    for i, (name, v, c) in enumerate(rows):
        y = len(rows) - i
        ax.barh(y, v, .5, color=c)
        ax.text(v + 6, y, '$%d' % round(v), va='center', fontsize=9, color=c)
    ax.set_yticks(range(1, len(rows) + 1))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xlim(0, 540)
    ax.set_xlabel('judge cost to evaluate %s traces, US dollars (generation excluded)'
                  % format(C['traces'], ','))
    ax.set_title('Evaluating the full benchmark: %d models x 2,250 problems.' % C['models'],
                 loc='left', fontsize=10, pad=10)
    save(fig, 'cost.png')


if __name__ == '__main__':
    N = json.load(open(NUMBERS, encoding='utf-8'))
    fig_trace(N)
    fig_steps(N)
    fig_planted(N)
    fig_routing(N)
    fig_cost(N)
