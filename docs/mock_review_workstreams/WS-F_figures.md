# WS-F: Figures (Phase 1)

## Amendment, 2026-10-08 (owner): refresh the data only

This section replaces the mission and steps F2 to F8 below. F1 is done (FIGURES EXTRACTED, by WS-C4).

Every figure in `overleaf_source_04102026/figs/` keeps its current design exactly; only its data is brought to the
final results. No redesign and no new figure: the heatmap and the all-models level figure are dropped, and the owner
commissions any further figure. `figs/app_candidates/` stays as it is: do not touch or remove it.

The figures, all drawn by `full_run_28092026/paper_figures.py`:
- `error-categories.pdf` (`fig:error_categories`): the PDF still shows the earlier error-analysis readings, while the
  caption block WS-G wrote gives the current ones (B2: 137 items). It must change.
- `branch-bars.pdf` (`fig:branch_bars`), `level-bars.pdf` (`fig:level_bars`) and
  `domain_radar/domain-radar-labeled.pdf` (`fig:domain_radar`).
- `level-gap.pdf`, `coverage-wrong.pdf` and `paraphrase.pdf` (in `figs/`, not placed in the paper).
- `domain_radar/domain-radar-unlabeled.pdf` (not placed; `fig_domain_radar` draws it beside the labelled one).

Steps:
1. Draw into a scratch copy of the tree: `python full_run_28092026/paper_results.py --out <scratch> --headline default
   --repaired` (without `--text-only`). The real tree is not touched.
2. Keep `paper_figures.py`'s design as it is: layout, colours, hatching, fonts, labels ("GPT OSS 20B"). If a figure
   cannot be drawn from the final data without a code change, make only the data-side change and report it.
3. Render each current and refreshed PDF to PNG and compare. The design must be identical. List what moved in the
   data (which bars or segments, by how much); a figure whose data did not change stays byte-identical.
4. Show the owner the before-and-after renders and wait for approval. Before approval: no tex edit, no commit.
5. On approval, copy the refreshed PDFs to the same paths in `figs/`.
6. Read each placed figure's caption block (WS-G regenerated them) and list any sentence that no longer fits the
   refreshed figure, for WS-D1 (it edits every caption in `paper_results.py`). Write no caption text yourself.

Signal: `SIGNAL: FIGURES REFRESHED <date time>` at the top of `reports/WS-F_report.md`, with, per figure: design
identical (yes or no), what changed in the data, the approval date, and the caption sentences that no longer fit.

## Mission

First, move the figure-drawing code out of `full_run_28092026/paper_results.py` into a module of its own, so
that WS-C can rework the tables in parallel. Then redraw the figures the reviewers found hard to read, in a form
that reads in greyscale, show them to the owner, and hand the approved PDFs and proposed captions over. No tex is
edited; captions are generated blocks that Phase 2 places.

Closes: review 1 §5 "Figures" (Figures 4, 5, 6; the model label), m7 (the figure side); review 2 W16 (naming),
§5 (Fig. 3, Fig. 4, the radar).

## Context (read `PLAN_CONTEXT.md` sections 1 to 3 first)

- The figures the paper uses: `overleaf_source_04102026/figs/level-bars.pdf` (FAC by level, four models),
  `error-categories.pdf` (stacked error categories, four models), `branch-bars.pdf` (FAC by branch as stacked
  fifths), `domain_radar/domain-radar-labeled.pdf` (FAC by domain, radial axis from 0.4), and the appendix
  figures `coverage-wrong.pdf`, `level-gap.pdf`, `paraphrase.pdf`. `paper_results.py --write` draws them (the
  drawing section starts near its line 1043, `import matplotlib`); `--write --text-only` leaves them.
- What the reviewers said: the error-category figure uses five shades of brown with dot hatching that cannot be
  told apart, and its "Unit" category is empty for every model; the branch figure's stacked fifths make the
  branches look equal and cannot be compared by eye (use grouped bars or a dot plot); the radar's axis from 0.4
  exaggerates differences and the table already carries the numbers (drop it or use a heatmap); the level figure
  shows four models and suggests a level effect for all, while the gap is not significant for the top tier after
  correction (show all eleven compactly in the appendix, or mark the non-significant gaps); figures say
  "GPT OSS 20B" while the text says `gpt-oss-20b`.
- The authors' earlier choices: figures label the model "GPT OSS 20B" on purpose (the code name reads oddly in a
  figure); the same label must appear in every figure. Decision D9 in `00_ORCHESTRATION.md` says whether a caption
  note or a switch resolves the inconsistency.
- Figures must read in greyscale (ACL); the column is 3.03 inches wide.

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, branch `master` only. Never create, switch to, list or inspect
  another branch.
- Other sessions work in this same working tree at the same time. Edit only the files under "Files you own".
  `paper_results.py` you edit exactly once (step F1), then never again; WS-C owns it after your signal.
- No paid API calls in this stream.
- No commits or pushes unless the owner asks in this session. When asked: one short line, no body, push at once.
- Figures: draw, render, show, wait. No tex edits, no other figure touched than the ones listed, no caption
  written into a tex file. Quality is not redesign: keep the figures' content and the paper's visual language;
  change what the reviewers named.
- Finish by writing `docs/mock_review_workstreams/reports/WS-F_report.md` from the template at the end. Post the
  signals at its top as soon as they are true.

## Before you start

- Read `paper_results.py` from the drawing section to the end, `figures_oct_12/make_overview5_large.py` (the
  taxonomy figure's style, for consistency), the figure captions in `overleaf_source_04102026/6_results.tex` and
  `appendices/branch_domain.tex`, and the dataviz guidance available in this environment (the `dataviz` skill),
  including its greyscale and colour-blind checks.
- The data every figure draws from: `full_run_28092026/results/results.json` (Q1, Q2, Q3, the branch and level
  tables) and `EXPERT_REQUEST.md` (B2 counts). WS-C may add keys; the final redraw happens at integration from
  the final `results.json`, so your drawing functions must read the data fresh, never embed numbers.

## Files you own

- `full_run_28092026/paper_figures.py` (new).
- `full_run_28092026/paper_results.py`: one edit, in step F1 only.
- `overleaf_source_04102026/figs/*.pdf` and the PNG renders you show (put renders in your scratchpad or
  `figs/_renders/`, gitignored if you create it).

Not yours: every `.tex`; `figures_oct_12/` (the taxonomy figure is out of scope); `results.json` (WS-C).

## Steps

### F1. Extract the drawing code (first 30 minutes)

1. Move every drawing function and the matplotlib setup from `paper_results.py` into `paper_figures.py`, with one
   entry point `draw_all(results: dict, figs_dir: Path, out_dir: Path | None = None)` and one function per
   figure, each taking the data it needs as arguments (no module-level reads of the tex tree).
2. In `paper_results.py`, replace the moved section by `from full_run_28092026 import paper_figures` and the call
   in `--write` (respecting `--text-only`). Nothing else in `paper_results.py` changes.
3. Regenerate the current figures and confirm they are unchanged (byte-identical, or identical renders if the
   PDF metadata differs). `paper_results.py --check` still passes.
4. Post **FIGURES EXTRACTED** in your report. From now on do not touch `paper_results.py`.

### F2. Error categories (`error-categories.pdf`, `fig:error_categories`)

Distinct hues for the categories, chosen so that they also separate in greyscale (vary lightness) and for
colour-blind readers; keep a light hatching only where two neighbours would otherwise merge; merge the empty
"Unit or dimension" category into "Other" (or drop it and say so in the caption); print the share in every
segment of at least 10% as now; label the hallucination segment for `gpt-oss-20b` even if under 10%.

### F3. Branch figure (`branch-bars.pdf`, `fig:branch_bars`)

Replace the stacked fifths with grouped bars (five branches × four models) or a dot plot with the model's
overall FAC as a reference line; draw the 95% template intervals per branch (they exist per model in the
branch table) and mark the one branch pair that differs after correction (`gpt-oss-20b`, electrical above
civil). The reader must be able to compare branches within a model by eye.

### F4. Domain figure (replace `domain_radar`, new label `fig:domain_heatmap`)

A heatmap of FAC by domain (rows: 15 domains grouped by branch; columns: the four representative models, or all
eleven if it stays legible at column width), greyscale-safe sequential palette, the numbers printed in the
cells, no axis that starts at 0.4. If the owner prefers to drop the figure (the table carries the numbers),
report that choice; default: the heatmap.

### F5. Level figure (`level-bars.pdf`, `fig:level_bars`)

Keep the four models, and mark each model's Easy-minus-Advanced gap as holding or not after Holm's correction
(a small marker or a hatched Advanced bar for "not significant"), with the detectable gap drawn as a bracket;
add a compact appendix version with all eleven models (`level-bars-all.pdf`, `fig:level_bars_all`), one row of
small multiples, same scale.

### F6. Labels

Per D9: keep "GPT OSS 20B" and provide the caption note text ("GPT OSS 20B is `gpt-oss-20b`") for every figure
that shows it, or switch the label in every figure to `gpt-oss-20b` in typewriter style. Consistent across all
figures, including the appendix ones you did not redraw (check and report).

### F7. Render, check, show, wait

Render every figure to PNG at the column width, convert to greyscale, and view both. Show the owner the
before-and-after renders and wait for approval before replacing the PDFs in `figs/`. If a figure is sent back,
revise only what was asked.

### F8. Captions

For each approved figure, write the proposed caption (sentence case, full stop, what the figure shows, the
sample it is drawn from, what is and is not tested) into your report under the figure's label. Phase 2 places
them in the generated figure blocks.

## Acceptance

- `paper_figures.py` draws every figure from `results.json` and the counts files; `paper_results.py --write`
  calls it; `--check` passes; the extraction regenerated the current figures identically.
- Every redrawn figure approved by the owner, legible at column width, legible in greyscale, consistent labels.
- The appendix figures not redrawn still render from the module.

## Signals

`SIGNAL: FIGURES EXTRACTED <date time>` after F1; `SIGNAL: CAPTIONS <date time>` after F8.

## Report template (`reports/WS-F_report.md`)

```
SIGNAL: FIGURES EXTRACTED <date time>
SIGNAL: CAPTIONS <date time>

# WS-F report

## Extraction
- Functions moved; the one-line change in paper_results.py; regeneration check (identical yes/no).

## Figures
| label | file | what changed | approved (date) |
- Dropped or added figures (fig:domain_heatmap replaces fig:domain_radar; fig:level_bars_all added).
- Label decision applied (D9).

## Captions (for Phase 2)
### fig:error_categories
<caption text>
### fig:branch_bars
...
## Open items
```
