"""The figures for the pilot summary.

    python evaluator_pilot_17092026/make_summary_figures.py

Every value below is a measured result, and the comment on each block names the analysis
that produced it, so a figure can always be traced back to the script that computed it:

    fig_trace       analysis/cluster_bootstrap.py   (trace-level AUROC, template CIs)
    fig_steps       analysis/digit_rule.py + annotation/score_against_labels.py
    fig_planted     analysis/planted.py + analysis/planted_judges.py
    fig_routing     analysis/planted_routing.py
    fig_cost        analysis/judge_cost.py at the 11-model roster

The palette matches the documents: ink for text, a single accent for the series that
carries the point, grey for everything it is being compared against.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures')

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


def fig_trace():
    """Trace level, with the intervals that account for 15 templates rather than 300 traces."""
    rows = [('E0 (published)', 0.850, 0.760, 0.927, GREY),
            ('E0 + 3rd judge', 0.810, 0.731, 0.885, GREY),
            ('E1 (judges swapped)', 0.849, 0.758, 0.927, GREY),
            ('E2 (PRM, fraction)', 0.862, 0.785, 0.943, GREY),
            ('E2 (PRM, minimum)', 0.878, 0.789, 0.956, GREY),
            ('E3 (milestones)', 0.834, 0.726, 0.941, GREY),
            ('E4 (+ arithmetic)', 0.835, 0.725, 0.941, GREY),
            ('E5 (milestones + judge)', 0.886, 0.784, 0.978, ACCENT),
            ('the experts\' answer verdict', 0.974, 0.945, 0.993, GOOD)]
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    for i, (name, v, lo, hi, c) in enumerate(rows):
        y = len(rows) - i
        ax.plot([lo, hi], [y, y], color=c, lw=2.4, solid_capstyle='round', alpha=.55)
        ax.plot([v], [y], 'o', color=c, ms=6)
        ax.text(hi + .006, y, '%.3f' % v, va='center', fontsize=8, color=c)
    ax.set_yticks(range(1, len(rows) + 1))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xlim(.68, 1.03)
    ax.set_xlabel('AUROC separating sound from unsound reasoning (95% interval, templates resampled)')
    ax.axvline(.850, color=LIGHT, lw=1, zorder=0)
    ax.set_title('No evaluator separates from the published framework.\n'
                 'The experts\' own answer verdict beats all of them.',
                 loc='left', fontsize=10, color=INK, pad=10)
    save(fig, 'trace_level.png')


def fig_steps():
    """Step level inside correct-answer traces: the case a reasoning evaluator exists for."""
    names = ['arithmetic check\nas shipped (1%)', 'best process\nreward model (72B)',
             'digit rule\n(as E4 now ships it)']
    prec = [0.154, 0.246, 0.750]
    rec = [0.034, 0.255, 0.320]
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
                 'Three flags in four are real once the tolerance is the displayed precision.',
                 loc='left', fontsize=10, pad=10)
    save(fig, 'step_level.png')


def fig_planted():
    """Planted defects: truth by construction, so the guide cannot be the reason."""
    groups = ['deterministic\nchecks', 'GPT-5', 'Claude Opus 4.5', 'MiMo-V2.5-Pro\n(the chosen judge)']
    arith = [0.683, 0.717, 0.717, 0.717]
    conc = [0.000, 0.333, 0.133, 0.314]
    x = range(len(groups))
    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ax.bar([i - .19 for i in x], arith, .36, label='arithmetic defects', color=LIGHT, edgecolor=GREY)
    ax.bar([i + .19 for i in x], conc, .36, label='conceptual defects', color=ACCENT)
    for i, (a, c) in enumerate(zip(arith, conc)):
        ax.text(i - .19, a + .02, '%.2f' % a, ha='center', fontsize=8, color=INK)
        ax.text(i + .19, c + .02, '%.2f' % c if c else '0', ha='center', fontsize=8,
                color=ACCENT if c else WARN)
    ax.set_xticks(list(x))
    ax.set_xticklabels(groups)
    ax.set_ylim(0, .95)
    ax.set_ylabel('share of 60 planted defects detected')
    ax.legend(frameon=False, loc='upper right')
    ax.set_title('Only a judge catches a misstated rule, and it catches about a third.\n'
                 'No false alarm on any untouched step, from any judge.',
                 loc='left', fontsize=10, pad=10)
    save(fig, 'planted.png')


def fig_routing():
    """Being able to catch a defect is not the same as being shown it."""
    fig, ax = plt.subplots(figsize=(6.6, 2.5))
    labels = ['shown the flawed step', 'judge catches it\nwhen shown', 'caught end to end']
    published = [0.500, 0.350, 0.175]
    routed = [1.000, 0.350, 0.350]
    x = range(3)
    ax.bar([i - .19 for i in x], published, .36, label='published routing', color=LIGHT, edgecolor=GREY)
    ax.bar([i + .19 for i in x], routed, .36, label='verify first, then route the residue', color=ACCENT)
    for i, (p, r) in enumerate(zip(published, routed)):
        ax.text(i - .19, p + .02, '%.3f' % p, ha='center', fontsize=8, color=INK)
        ax.text(i + .19, r + .02, '%.3f' % r, ha='center', fontsize=8, color=ACCENT)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel('conceptual defects')
    ax.legend(frameon=False, loc='upper left', ncol=1)
    ax.set_title('The published framework shows a judge the flawed step as often as a clean one.\n'
                 'Routing what cannot be verified doubles what is caught.',
                 loc='left', fontsize=10, pad=10)
    save(fig, 'routing.png')


def fig_cost():
    """What each choice costs over the full benchmark, 11 models x 2,250 problems."""
    rows = [('published framework\n(judges on every trace)', 334, GREY),
            ('recommended stack\n(deterministic + residual judge)', 77, ACCENT),
            ('with the router\n(judge on the unverifiable residue)', 156, GOOD)]
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    for i, (name, v, c) in enumerate(rows):
        y = len(rows) - i
        ax.barh(y, v, .5, color=c)
        ax.text(v + 6, y, '$%d' % v, va='center', fontsize=9, color=c)
    ax.set_yticks(range(1, len(rows) + 1))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xlim(0, 400)
    ax.set_xlabel('cost to evaluate 24,750 traces (generation excluded)')
    ax.set_title('Evaluating the full benchmark, three ways.', loc='left', fontsize=10, pad=10)
    save(fig, 'cost.png')


if __name__ == '__main__':
    fig_trace()
    fig_steps()
    fig_planted()
    fig_routing()
    fig_cost()
