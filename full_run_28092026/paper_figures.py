"""The paper's figures, drawn from the result data that paper_results.py loads. Extracted from paper_results.py on
2026-10-07 (WS-F step F1, done by WS-C4): the same drawing code, each function taking the data it needs as arguments
and the folder it writes into; nothing here reads the tex tree or a result file.

    draw_all(results, figs_dir, out_dir=None)   # the placed figures (FIGURES) into out_dir, or figs_dir when out_dir is None

`results` is the dict paper_results.figure_data() assembles: order (model keys by FAC), q1, q2, q3, q5 (per-model dicts
of results.json), bl (branches_levels per model), rep (reported.models), by_model (the B2 readings per model),
categories (the (full label, short label) list), incomplete (the full label of the incomplete option), b2_models,
representative (the four models the figures show), branch (key -> name), levels, n_templates, margin (the paraphrase
margin), fig_name and name (model display names for figures and text), branch_of (domain -> branch).

Style: one hue, text in ink, recessive axes; marker fill is the second encoding, so the figures read in greyscale.
"""
from __future__ import annotations

from pathlib import Path

BLUE, DARK, LIGHT, INK, MUTED, BAND = "#2a78d6", "#0d366b", "#cde2fb", "#0b0b0b", "#52514e", "#eceae6"
HATCH = "#d3d3d3"  # the hatch lines of the error-category figure
COLUMN = 3.03  # the ACL column width in inches
SERIF = {"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 7, "legend.fontsize": 7}  # the figures' type

FIGURES = (  # the four placed figures; fig_domain_radar also draws the unlabeled radar beside the labeled one
    "branch-bars.pdf", "error-categories.pdf", "level-bars.pdf", "domain_radar/domain-radar-labeled.pdf")


def wrap(name: str, width: int = 15) -> str:
    """A figure label over two lines, broken at its last space, when it is longer than width characters."""
    return name[::-1].replace(" ", "\n", 1)[::-1] if len(name) > width and " " in name else name


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif", "font.size": 6.5, "axes.linewidth": 0.4,
                         "xtick.major.width": 0.4, "ytick.major.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 0,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.labelsize": 6.5, "legend.fontsize": 6,
                         "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": INK, "axes.labelcolor": INK,
                         "hatch.linewidth": 0.4, "savefig.dpi": 300})
    return plt


def vector_hatches(fig) -> None:
    """Redraw each hatched patch's hatch as clipped vector lines. Matplotlib writes hatches as PDF tiling patterns, which PDF
    viewers rasterize at low resolution; the lines keep the pattern's geometry, colour and width (72 pt cells anchored at the
    page's top left), so the figure looks the same, only sharp."""
    import numpy as np
    import matplotlib as mpl
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path as MPath
    fig.canvas.draw()  # fixes every extent, the legends' included
    to_inches, dy = fig.dpi_scale_trans.inverted(), fig.get_figheight() % 1.0
    legends = list(fig.legends) + [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
    owner = {id(h): lg for lg in legends for h in lg.get_patches()}
    for patch in fig.findobj(lambda a: isinstance(a, Patch) and bool(a.get_hatch())):
        (x0, y0), (x1, y1) = to_inches.transform(patch.get_window_extent().get_points())
        cell = MPath.hatch(patch.get_hatch())
        tiles = [cell.vertices + (i, j + dy) for i in range(int(np.floor(x0)) - 1, int(np.ceil(x1)) + 1)
                 for j in range(int(np.floor(y0 - dy)) - 1, int(np.ceil(y1 - dy)) + 1)]
        color = patch.get_hatchcolor()
        width = patch.get_hatch_linewidth() if hasattr(patch, "get_hatch_linewidth") else mpl.rcParams["hatch.linewidth"]
        legend = owner.get(id(patch))
        lines = PathPatch(MPath(np.concatenate(tiles), np.concatenate([cell.codes] * len(tiles))), transform=fig.dpi_scale_trans,
                          facecolor=color, edgecolor=color, linewidth=width,
                          zorder=(legend.get_zorder() if legend is not None else patch.get_zorder()) + 0.01)
        patch.set_hatch(None)
        parent = legend.axes if legend is not None and legend.axes is not None else (fig if legend is not None else patch.axes or fig)
        parent.add_artist(lines)
        lines.set_clip_path(patch)  # after adding: an axes clips a new artist to itself when it has no clip path of its own


def save(fig, figs_dir: Path, name: str) -> None:
    figs_dir.mkdir(parents=True, exist_ok=True)
    vector_hatches(fig)
    fig.savefig(figs_dir / name, metadata={"CreationDate": None})
    import matplotlib.pyplot as plt
    plt.close(fig)


def recess(ax, left: bool = True) -> None:
    for s in ("top", "right") + (() if left else ("left",)):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e7e6e2", linewidth=0.4)
    ax.set_axisbelow(True)


def rows_axes(plt, n: int, height: float, width: float = COLUMN, left: float = 0.36, right: float = 0.97, bottom: float = 0.13,
              top: float = 0.90, ncols: int = 1, wspace: float = 0.08):
    fig, axes = plt.subplots(1, ncols, figsize=(width, height), sharey=True,
                             gridspec_kw={"left": left, "right": right, "bottom": bottom, "top": top, "wspace": wspace})
    axes = list(axes) if ncols > 1 else [axes]
    for ax in axes:
        ax.set_ylim(n - 0.5, -0.5)
        recess(ax, left=ax is axes[0])
    return fig, axes


def interval_rows(ax, items, lw: float = 1.0) -> None:
    """items: (y, lo, hi, x, filled, marker) rows: an interval line with a marker, filled or hollow."""
    for y, lo, hi, x, filled, marker in items:
        ax.plot([lo, hi], [y, y], color=BLUE, linewidth=lw, solid_capstyle="butt", zorder=2)
        ax.plot([x], [y], marker=marker, markersize=4.2, markeredgewidth=0.8, markeredgecolor=BLUE,
                markerfacecolor=BLUE if filled else "white", linestyle="none", zorder=3)


def gradient_rows(ax, items, height: float = 0.38, steps: int = 80) -> None:
    """items: (y, lo, hi, x, color) rows: an interval band in color, deepest at the estimate x and paler toward its ends. Each
    slice runs on under the next, which is drawn over it, so no seam shows between slices."""
    from matplotlib.colors import to_rgb
    for y, lo, hi, x, color in items:
        rgb, span, step = to_rgb(color), max(x - lo, hi - x), (hi - lo) / steps
        for j in range(steps):
            t = abs(lo + (j + 0.5) * step - x) / span  # 0 at the estimate, 1 at the farther end
            ax.barh(y, step * (2 if j < steps - 1 else 1), left=lo + j * step, height=height, color=[c + (1 - c) * 0.72 * t for c in rgb],
                    linewidth=0, zorder=2)


def top_legend(fig, handles, ncol: int, y: float = 0.995) -> None:
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, frameon=False, handlelength=1.4,
               columnspacing=1.0, handletextpad=0.5)


def fig_level_gap(order: list[str], q2: dict, fig_name: dict, figs_dir: Path) -> None:
    """The Easy-minus-Advanced gap per model with its interval (level-gap.pdf; not placed in the paper)."""
    plt = _plt()
    from matplotlib.lines import Line2D
    # The interval as a band, deepest at the gap; the gaps that hold in strong blue with a dark filled marker, the rest paler and hollow.
    hold, rest, rest_edge = BLUE, "#93b6e2", "#5b8fd3"
    holds = {k: q2[k]["p_welch_holm"] < 0.05 for k in order}
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(order), 2.6, left=0.318, right=0.97, bottom=0.12, top=0.985)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([fig_name[k] for k in order], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        ax.axvline(0, color=MUTED, linewidth=0.6, zorder=1)
        ax.text(-0.004, -0.42, "No gap", rotation=90, ha="right", va="top", fontsize=5.5, style="italic", color=MUTED)
        gradient_rows(ax, [(i, q2[k]["ci"][0], q2[k]["ci"][1], q2[k]["gap"], hold if holds[k] else rest) for i, k in enumerate(order)])
        for i, k in enumerate(order):
            ax.plot([q2[k]["gap"]], [i], marker="o", markersize=4.6, linestyle="none", zorder=3,
                    **({"markerfacecolor": DARK, "markeredgecolor": "white", "markeredgewidth": 0.6} if holds[k] else
                       {"markerfacecolor": "white", "markeredgecolor": rest_edge, "markeredgewidth": 0.9}))
        ax.set_xlim(-0.04, 0.36)
        ax.set_xlabel("Final Answer Accuracy drop from Easy to Advanced", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="Holds after Holm"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=rest_edge, markeredgewidth=0.9, markersize=4.6,
                          linestyle="none", label="Does not hold"),
                   Line2D([], [], color="#7fb0ea", linewidth=4.5, solid_capstyle="butt", label="95% interval")]
        legend = ax.legend(handles=handles, loc="upper right", frameon=True, fancybox=False, framealpha=1, edgecolor="black", fontsize=6,
                           borderaxespad=0.3, borderpad=0.35, handlelength=1.3, handletextpad=0.4, labelspacing=0.3)
        legend.get_frame().set_linewidth(0.5)
        save(fig, figs_dir, "level-gap.pdf")


def fig_coverage_wrong(order: list[str], q3: dict, fig_name: dict, figs_dir: Path) -> None:
    """Milestone Coverage of the wrong answers per model (coverage-wrong.pdf; not placed in the paper)."""
    plt = _plt()
    from matplotlib.lines import Line2D
    # A band from the chance floor, palest, to the coverage with the judge, deepest; matching alone sits on it.
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(order), 2.8, left=0.318, right=0.91, bottom=0.111, top=0.877)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([fig_name[k] for k in order], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        rows = [(i, q3[k]["e3_null_on_readable_wrong"], q3[k]["e3_coverage_on_readable_wrong"], q3[k]["e5_coverage_on_readable_wrong"])
                for i, k in enumerate(order)]
        gradient_rows(ax, [(i, f, e5, e5, BLUE) for i, f, _, e5 in rows])
        for i, f, e3, e5 in rows:
            ax.plot([f], [i], marker="s", markersize=3.8, color=MUTED, markeredgecolor="white", markeredgewidth=0.5, linestyle="none", zorder=3)
            ax.plot([e3], [i], marker="o", markersize=4.6, markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=0.9, linestyle="none",
                    zorder=3)
            ax.plot([e5], [i], marker="o", markersize=4.6, markerfacecolor=DARK, markeredgecolor="white", markeredgewidth=0.6, linestyle="none",
                    zorder=4)
        for i, k in enumerate(order):
            ax.text(1.025, i, str(q3[k]["readable_wrong_with_milestones"]), va="center", ha="left", fontsize=6.5, color="black", clip_on=False)
        ax.text(1.025, -0.8, "n", ha="left", va="center", fontsize=6.5, style="italic", color="black", clip_on=False)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("Milestone Coverage of wrong answers", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="s", color=MUTED, markersize=3.8, linestyle="none", label="Chance floor"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=0.9, markersize=4.6, linestyle="none",
                          label="Matching alone"),
                   Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="With the judge")]
        legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=((0.318 + 0.91) / 2, 0.995), ncol=3, frameon=True,
                            fancybox=False, framealpha=1, edgecolor="black", fontsize=6, borderpad=0.35, handlelength=1.0,
                            handletextpad=0.4, columnspacing=1.2)
        legend.get_frame().set_linewidth(0.5)
        save(fig, figs_dir, "coverage-wrong.pdf")


def band_axes(plt, n: int, height: float, xlabel: str, margin: float, left: float = 0.36, top: float = 0.88):
    fig, (ax,) = rows_axes(plt, n, height, left=left, top=top)
    ax.axvspan(-margin, margin, color=BAND, zorder=0)
    ax.axvline(0, color=MUTED, linewidth=0.5, linestyle=(0, (2, 2)), zorder=1)
    ax.set_xlabel(xlabel)
    return fig, ax


def fig_paraphrase(order: list[str], q5: dict, margin: float, fig_name: dict, figs_dir: Path) -> None:
    """The FAC change under paraphrase per model with its 90% interval (paraphrase.pdf; not placed in the paper)."""
    plt = _plt()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    # The 90% interval as a band, deepest at the change; within the margin in strong blue with a dark filled marker, the rest paler and hollow.
    rest, rest_edge = "#93b6e2", "#5b8fd3"
    within = {k: q5[k]["within_margin"] for k in order}
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(order), 2.6, left=0.318, right=0.72, bottom=0.12, top=0.985)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([fig_name[k] for k in order], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        ax.axvspan(-margin, margin, color=BAND, zorder=0)
        ax.axvline(0, color=MUTED, linewidth=0.6, zorder=1)
        ax.text(0.0015, 0.5, "No change", rotation=90, ha="left", va="center", fontsize=5.5, style="italic", color=MUTED)
        gradient_rows(ax, [(i, q5[k]["ci90"][0], q5[k]["ci90"][1], q5[k]["diff"], BLUE if within[k] else rest) for i, k in enumerate(order)])
        for i, k in enumerate(order):
            ax.plot([q5[k]["diff"]], [i], marker="o", markersize=4.6, linestyle="none", zorder=3,
                    **({"markerfacecolor": DARK, "markeredgecolor": "white", "markeredgewidth": 0.6} if within[k] else
                       {"markerfacecolor": "white", "markeredgecolor": rest_edge, "markeredgewidth": 0.9}))
        ax.set_xlim(-0.075, 0.061)
        ax.set_xticks([-0.05, 0, 0.05])
        ax.set_xticklabels(["−0.05", "0", "0.05"])
        ax.set_xlabel("Final Answer Accuracy change under paraphrase", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="Within the margin"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=rest_edge, markeredgewidth=0.9, markersize=4.6,
                          linestyle="none", label="Not within"),
                   Line2D([], [], color="#7fb0ea", linewidth=4.5, solid_capstyle="butt", label="90% interval"),
                   Patch(color=BAND, label=f"±{margin:.2f} margin")]
        legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=True, fancybox=False, framealpha=1,
                           edgecolor="black", fontsize=6, borderaxespad=0, borderpad=0.35, handlelength=1.1, handletextpad=0.4,
                           labelspacing=0.3)
        legend.get_frame().set_linewidth(0.5)
        save(fig, figs_dir, "paraphrase.pdf")


def fig_error_categories(by_model: dict, categories: list, incomplete: str, b2_models: list[str], figs_dir: Path) -> None:
    """The experts' error categories of the wrong answers read, one stacked bar per model (error-categories.pdf)."""
    plt = _plt()
    from matplotlib.patches import Patch
    cats = [c for c in categories if c[0] != incomplete]  # no "Incomplete" reading was given
    # The May figure's palette: the errors before calculation in browns, darker the more fundamental, calculation in hatched
    # mustard, no error in green; the lightness order and the hatching carry the categories in greyscale.
    fills = ["#5c3b24", "#7f5638", "#a87c5b", "#b78f71", "#c3a084", "#d9b86c", "#b3dda0"]
    hatches = ["...", "...", "...", "...", "...", "//", ""]
    labels = {"claude-sonnet-5": "Claude\nSonnet 5", "gpt-5.4-mini": "GPT-5.4\nmini", "gemma-4-26b-a4b": "Gemma 4\n26B", "gpt-oss-20b": "GPT OSS\n20B"}
    b2_models = sorted(b2_models, key=list(labels).index)  # the committed bar order, whatever order the data arrives in
    with plt.rc_context({"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 7, "legend.fontsize": 7}):
        fig, ax = plt.subplots(figsize=(COLUMN, 2.3), gridspec_kw={"left": 0.11, "right": 0.75, "bottom": 0.15, "top": 0.97})
        for x, m in enumerate(b2_models):
            n, bottom = sum(by_model[m].values()), 0.0
            for (full, _), fill, hatch in reversed(list(zip(cats, fills, hatches))):  # no error at the base, the most fundamental on top
                v = by_model[m].get(full, 0) / n
                ax.bar(x, v, bottom=bottom, width=0.6, color=fill, linewidth=0, zorder=2)
                if hatch:
                    ax.bar(x, v, bottom=bottom, width=0.6, fill=False, hatch=hatch, edgecolor=HATCH, linewidth=0, zorder=3)
                ax.bar(x, v, bottom=bottom, width=0.6, fill=False, edgecolor="white", linewidth=0.6, zorder=4)
                if v >= 0.1:
                    ax.text(x, bottom + v / 2, f"{v * 100:.0f}", ha="center", va="center", fontsize=7, fontweight="bold", zorder=5,
                            color="white" if fill in fills[:2] else "black", bbox={"facecolor": fill, "edgecolor": "none", "pad": 0.6})
                bottom += v
        ax.set_xticks(range(len(b2_models)))
        ax.set_xticklabels([labels[m] for m in b2_models], fontweight="bold")
        ax.tick_params(axis="x", length=0, colors="black", labelsize=7)
        ax.set_xlim(-0.42, len(b2_models) - 0.58)
        ax.tick_params(axis="y", colors="black", labelsize=7)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0", "25", "50", "75", "100"])
        ax.set_ylabel("Share of the readings (%)", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        handles = [Patch(facecolor=fill, hatch=hatch, edgecolor=HATCH, linewidth=0, label=short.split(" or ")[0])
                   for (_, short), fill, hatch in zip(cats, fills, hatches)]
        legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.025, 1.0), frameon=True, fancybox=False, framealpha=1,
                           edgecolor="black", fontsize=6, borderaxespad=0, borderpad=0.35, handlelength=1.1, handleheight=0.9,
                           handletextpad=0.4, labelspacing=0.25)
        legend.get_frame().set_linewidth(0.5)
        vector_hatches(fig)  # the hatching as vector lines, so it stays sharp in PDF viewers; nothing else changes
        save(fig, figs_dir, "error-categories.pdf")


MODEL_COLORS = {"deepseek-v4.1-flash": "#2a78d6", "claude-sonnet-5": "#1baf7a", "gpt-5.4-mini": "#4a3aa7", "gpt-oss-20b": "#eb6834"}
HATCHES = ["", "//", "..", "xx", "--"]


FIG_NAME_2 = {"deepseek-v4.1-flash": "DeepSeek\nV4.1 Flash", "claude-sonnet-5": "Claude\nSonnet 5", "gpt-5.4-mini": "GPT-5.4\nmini",
              "gpt-oss-20b": "GPT OSS\n20B"}  # the figure labels on two lines, for the column width
BAR_LEFT, BAR_RIGHT = 0.10, 0.995


def column_bar_axes(plt):
    return plt.subplots(figsize=(COLUMN, 2.45), gridspec_kw={"left": BAR_LEFT, "right": BAR_RIGHT, "bottom": 0.155, "top": 0.80})


def finish_column_bars(fig, ax, representative: list[str], labels: list[str], hatches: list[str], legend_title: str) -> None:
    """The shared frame of the two bar figures: model names under the bars, the FAC axis, the boxed legend of the hatches."""
    from matplotlib.patches import Patch
    ax.set_xticks(range(len(representative)))
    ax.set_xticklabels([FIG_NAME_2[m] for m in representative], fontsize=7, fontweight="bold", linespacing=1.0)
    ax.set_xlim(-0.55, len(representative) - 0.45)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("Final Answer Accuracy", fontsize=6.5, fontweight="bold", labelpad=2)
    ax.grid(axis="y", linestyle="--", color="#d9d8d4", linewidth=0.4)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=6, pad=1.5)
    handles = [Patch(facecolor="#8c8c8c" if not h else "#e6e6e6", edgecolor=INK, hatch=h, linewidth=0.5, label=l) for l, h in zip(labels, hatches)]
    legend = fig.legend(handles=handles, title=legend_title, loc="upper center", bbox_to_anchor=((BAR_LEFT + BAR_RIGHT) / 2, 1.0),
                        ncol=len(labels), frameon=True, edgecolor=INK, fancybox=False, fontsize=6, title_fontsize=6.5, handlelength=1.6,
                        handleheight=0.9, handletextpad=0.4, columnspacing=0.9, borderpad=0.35)
    legend.get_title().set_fontweight("bold")
    legend.get_frame().set_linewidth(0.5)
    vector_hatches(fig)


def grouped_bars(name: str, representative: list[str], groups: list[tuple[str, str]], legend_title: str, values, figs_dir: Path) -> None:
    """The May layout at the ACL column width: one color per representative model, one hatch per group, the value above every
    bar, a boxed legend; values(model, group label) returns (mean, (lo, hi))."""
    plt = _plt()
    import numpy as np
    from matplotlib.colors import to_rgb
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "hatch.linewidth": 0.4})
    fig, ax = column_bar_axes(plt)
    x = np.arange(len(representative))
    width = 0.8 / len(groups)
    for i, model in enumerate(representative):
        color = MODEL_COLORS[model]
        tint = tuple(1 - 0.3 * (1 - c) for c in to_rgb(color))
        for j, (label, hatch) in enumerate(groups):
            mean, (lo, hi) = values(model, label)
            pos = x[i] - 0.4 + width * (j + 0.5)
            ax.bar(pos, mean, width=width * 0.92, facecolor=color if not hatch else tint, edgecolor=color, hatch=hatch, linewidth=0.5, zorder=2)
            ax.errorbar(pos, mean, yerr=[[mean - lo], [hi - mean]], fmt="none", ecolor=INK, elinewidth=0.45, capsize=1.0, capthick=0.45, zorder=3)
            ax.text(pos, hi + 0.012, f"{mean:.2f}", ha="center", va="bottom", fontsize=5.5, fontweight="bold", color=INK)
    finish_column_bars(fig, ax, representative, [g[0] for g in groups], [g[1] for g in groups], legend_title)
    save(fig, figs_dir, name)


def fig_level_bars(representative: list[str], bl: dict, levels: list[str], figs_dir: Path) -> None:
    """FAC by difficulty level for the four representative models (level-bars.pdf)."""
    grouped_bars("level-bars.pdf", representative, list(zip(levels, HATCHES)), "Difficulty Level",
                 lambda m, lv: (bl[m]["level"][lv]["mean"], bl[m]["level"][lv]["ci"]), figs_dir)


def text_on(fill) -> str:
    """White or black, whichever contrasts more with the fill (WCAG relative luminance)."""
    from matplotlib.colors import to_rgb
    r, g, b = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in to_rgb(fill)]
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if 1.05 / (y + 0.05) > (y + 0.05) / 0.05 else "black"


def fig_branch_bars(representative: list[str], bl: dict, q1: dict, branch: dict, n_templates: int, figs_dir: Path) -> None:
    """FAC by branch for the four representative models at the ACL column width, as stacked bars. The five branches hold 30
    templates each, so a segment is one branch's share of the model's FAC (the branch mean divided by five) and a bar's height
    is the model's FAC, printed above it; the label in a segment gives the branch mean. The colors and hatches are those of the
    level figure."""
    plt = _plt()
    from matplotlib.colors import to_rgb
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "hatch.linewidth": 0.4})
    fig, ax = column_bar_axes(plt)
    width, n = 0.56, len(branch)
    for i, model in enumerate(representative):
        color = MODEL_COLORS[model]
        tint = tuple(1 - 0.3 * (1 - c) for c in to_rgb(color))
        bottom = 0.0
        for (key, _), hatch in zip(branch.items(), HATCHES):
            assert bl[model]["branch"][key]["templates"] * n == n_templates
            mean = bl[model]["branch"][key]["mean"]
            fill = color if not hatch else tint
            ax.bar(i, mean / n, bottom=bottom, width=width, facecolor=fill, edgecolor=color, hatch=hatch, linewidth=0.5, zorder=2)
            if bottom > 0:  # a white seam between segments
                ax.plot([i - width / 2, i + width / 2], [bottom, bottom], color="white", linewidth=0.8, solid_capstyle="butt", zorder=3)
            ax.text(i, bottom + mean / n / 2, f"{mean:.2f}", ha="center", va="center", fontsize=6, fontweight="bold", zorder=4,
                    color=text_on(fill), bbox={"facecolor": fill, "edgecolor": "none", "pad": 0.5})
            bottom += mean / n
        assert abs(bottom - q1[model]["score"]) < 1e-9  # the stack is the model's FAC
        ax.text(i, bottom + 0.012, f"{bottom:.2f}", ha="center", va="bottom", fontsize=6, fontweight="bold", color=INK)
    finish_column_bars(fig, ax, representative, list(branch.values()), HATCHES[:n], "Engineering Branch")
    save(fig, figs_dir, "branch-bars.pdf")


DOMAIN_LABEL = {  # the paper's domain names, broken for the radar's rim
    "reaction_kinetics": "Reaction\nKinetics", "thermodynamics": "Thermodynamics", "transport_phenomena": "Transport\nPhenomena",
    "digital_communications": "Digital\nCommunications", "electromagnetics_and_waves": "Electromagnetics\nand Waves",
    "signals_and_systems": "Signals and\nSystems", "fluid_mechanics": "Fluid\nMechanics", "mechanics_of_materials": "Mechanics of\nMaterials",
    "vibrations_and_acoustics": "Vibrations and\nAcoustics", "geotechnical_engineering": "Geotechnical\nEngineering",
    "structural_analysis": "Structural\nAnalysis", "water_resources": "Water\nResources", "production_and_inventory": "Production and\nInventory",
    "quality_and_reliability_control": "Quality and\nReliability Control", "stochastic_operations": "Stochastic\nOperations",
}
RADAR_BRANCH_ORDER = ["chemical_engineering", "electrical_engineering", "mechanical_engineering", "civil_engineering", "industrial_engineering"]


def fig_domain_radar(representative: list[str], rep: dict, branch_of: dict, fig_name: dict, name: dict, figs_dir: Path) -> None:
    """FAC by domain for the four representative models, domains grouped by branch around the rim, in two versions:
    with the domain names (domain-radar-labeled.pdf) and without (domain-radar-unlabeled.pdf), under figs_dir/domain_radar/."""
    import numpy as np
    from matplotlib.lines import Line2D
    plt = _plt()
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"]})
    domains = [d for b in RADAR_BRANCH_ORDER for d in sorted((d for d, br in branch_of.items() if br == b), key=lambda d: DOMAIN_LABEL[d])]
    assert len(domains) == len(DOMAIN_LABEL) == 15
    angles = np.linspace(0, 2 * np.pi, len(domains), endpoint=False)
    closed = np.append(angles, angles[0])
    styles = [((0, (1, 1.2)), "o"), ((0, (4, 1, 1, 1, 1, 1)), "^"), ((0, (4, 1.5, 1, 1.5)), "s"), ("-", "D")]
    (figs_dir / "domain_radar").mkdir(parents=True, exist_ok=True)
    for labeled in (True, False):
        fig = plt.figure(figsize=(4.8, 4.2))
        ax = fig.add_axes([0.20, 0.13, 0.60, 0.65], polar=True)
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        handles = []
        for m, (ls, mk) in zip(representative, styles):
            vals = [rep[m]["domain"][d] for d in domains]
            vals = np.append(vals, vals[0])
            color = MODEL_COLORS[m]
            ax.fill(closed, vals, color=color, alpha=0.05, zorder=1)
            ax.plot(closed, vals, linestyle=ls, linewidth=1.1, color=color, marker=mk, markersize=3.6, markerfacecolor="white",
                    markeredgewidth=0.9, zorder=3)
            handles.append(Line2D([], [], linestyle=ls, linewidth=1.1, color=color, marker=mk, markersize=3.6, markerfacecolor="white",
                                  markeredgewidth=0.9, label=fig_name.get(m, name[m])))
        ax.set_xticks(angles)
        ax.set_xticklabels([DOMAIN_LABEL[d] for d in domains] if labeled else [], fontsize=6.5, fontweight="bold", color=INK)
        ax.tick_params(axis="x", pad=-1)
        for label, theta in zip(ax.get_xticklabels(), angles):  # align each rim label away from the circle
            x, y = np.cos(np.pi / 2 - theta), np.sin(np.pi / 2 - theta)
            label.set_ha("left" if x > 0.15 else "right" if x < -0.15 else "center")
            label.set_va("bottom" if y > 0.15 else "top" if y < -0.15 else "center")
        ax.set_ylim(0.4, 1.07)  # the rim sits outside the 1.0 ring, so the polygons never touch it
        ax.set_yticks([0.6, 0.8, 1.0])
        ax.set_yticklabels(["0.6", "0.8", "1.0"], fontsize=5.5, color=MUTED)
        ax.set_rlabel_position(90)
        for label in ax.get_yticklabels():
            label.set_bbox({"facecolor": "white", "edgecolor": "none", "pad": 0.4})
        ax.grid(color="#cfcecb", linewidth=0.5)
        ax.spines["polar"].set_color(INK)
        ax.spines["polar"].set_linewidth(0.8)
        legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=4, frameon=True, fancybox=False,
                            edgecolor=INK, fontsize=7, handlelength=2.2, columnspacing=1.0)
        legend.get_frame().set_linewidth(0.6)
        save(fig, figs_dir, f"domain_radar/domain-radar-{'labeled' if labeled else 'unlabeled'}.pdf")


def draw_all(results: dict, figs_dir: Path, out_dir: Path | None = None) -> tuple[str, ...]:
    """Draw the placed figures (FIGURES) into out_dir, or into figs_dir when out_dir is None, in the order they were drawn
    before the extraction; returns the file names drawn. `results` is described in the module docstring."""
    dest = Path(out_dir) if out_dir is not None else Path(figs_dir)
    r = results
    fig_branch_bars(r["representative"], r["bl"], r["q1"], r["branch"], r["n_templates"], dest)
    fig_error_categories(r["by_model"], r["categories"], r["incomplete"], r["b2_models"], dest)
    fig_level_bars(r["representative"], r["bl"], r["levels"], dest)
    fig_domain_radar(r["representative"], r["rep"], r["branch_of"], r["fig_name"], r["name"], dest)
    return FIGURES
