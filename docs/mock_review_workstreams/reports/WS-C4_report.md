SIGNAL: FIGURES EXTRACTED 2026-10-07 14:35 UTC (WS-F's step F1 done by WS-C4; `paper_figures.py` is WS-F's from here on)

# WS-C4 report: Table 1 and the generated blocks

Done, 2026-10-07. Nothing is committed (none asked). No paid call. No `.tex` under `overleaf_source_04102026/` or
`current_overleaf_project/` was written: every block went into scratch copies under `--out DIR`. No other stream's
file was edited, with two exceptions that the brief allows or implies: `paper_figures.py` was created once (step F1,
on WS-F's behalf) and `.gitignore` gained one line for the stand-in folder.

## Extraction (WS-F step F1)

- Moved from `paper_results.py` into `full_run_28092026/paper_figures.py`: the style constants, `_plt`,
  `vector_hatches` (the file had two identical copies; one kept), `save`, `recess`, `rows_axes`, `interval_rows`,
  `gradient_rows`, `top_legend`, `fig_level_gap`, `fig_coverage_wrong`, `band_axes`, `fig_paraphrase`,
  `fig_error_categories`, `column_bar_axes`, `finish_column_bars`, `grouped_bars`, `fig_level_bars`, `text_on`,
  `fig_branch_bars`, `fig_domain_radar`, `MODEL_COLORS`, `HATCHES`, `FIG_NAME_2`, `DOMAIN_LABEL`,
  `RADAR_BRANCH_ORDER`, and the dead `wrap()`. Every function now takes the data it needs and the folder it writes
  into; nothing in the module reads the tex tree or a result file.
- Entry point: `paper_figures.draw_all(results, figs_dir, out_dir=None)`, drawing `paper_figures.FIGURES` (the four
  placed figures; the radar function also draws the unlabeled radar). `results` is the dict
  `paper_results.figure_data()` assembles (keys documented in the module docstring).
- In `paper_results.py`: `from full_run_28092026 import paper_figures`, `figure_data()`, and the call in `--write`
  (`--text-only` respected). `--check` looks for `paper_figures.FIGURES`.
- Regeneration check: the five PDFs (`branch-bars`, `error-categories`, `level-bars`, the labeled and unlabeled radar)
  drawn through the extracted module are byte-identical to those drawn by the unmodified code, which are themselves
  byte-identical to the PDFs in the tree (matplotlib 3.11.1). The tree's `figs/` was not touched.

## Scripts and flags

| script | what it does | reads | writes | check |
|---|---|---|---|---|
| `full_run_28092026/paper_results.py` | every table and figure block of §5.3, §5.4 and their appendices, the phrases the prose must contain, the stand-ins | `results/results.json`, `matched_config.json`, the seven WS-C files (or their stand-ins), `expert_request/scored_current.json`, the Markdown records, `subsamples.py`, the answer module's symbolic list, `appendices/validation.tex` | the blocks (real tree with `--write`; a copy with `--out DIR`), `DIR/generated/*.tex` for blocks without markers, `DIR/generated/phrases.txt`, the figures through `paper_figures` | `--check [--out DIR]` |
| `full_run_28092026/paper_figures.py` (new, WS-F's) | the drawing code | the dict `figure_data()` passes | the four placed figures | regeneration byte-identical (above) |

Flags added to `paper_results.py` (the registry's names): `--out DIR`, `--headline default|matched`, `--repaired`,
`--text-only` (existed), `--judged-in-table`, `--stand-in`. `--check` takes `--out DIR` to check the copy.

- `--out DIR` copies `overleaf_source_04102026/` to DIR (`copytree`, existing files overwritten), writes every block
  wherever in the copy its markers are (any `.tex` file, not only the default one, so the writers can move a marker
  pair and the block follows), writes a block that has no markers to `DIR/generated/<label>.tex` (full marker-wrapped
  block, ready to paste), writes every phrase to `DIR/generated/phrases.txt`, and draws the figures into `DIR/figs`
  unless `--text-only`. `--out` refuses a DIR inside the real tree.
- `--headline default`: Table 1 holds the eleven at the providers' defaults (from `results.json`) and a marked block
  "With reasoning at medium effort" with the three re-run models from `matched.json`, tier letter `--` (their letters
  belong to the matched family in `tab:matched`). `--headline matched`: every model at its matched store from
  `matched.json`, re-run models marked `‡`, a twelfth model (any key in `matched.json` not in `results.json`) appended
  under "Added model"; `tab:matched` holds both configurations under either headline.
- `--repaired`: drops the "with and without the two chemical templates" phrase, the B4-reading phrases on the two
  templates, the "whose wording does not pin the answer" clauses, the `tab:level_gap` column "Without two chemical"
  and its caption clause, and the claim checks on the B4 answers.
- `--judged-in-table`: the judged step flags stay as Table 1's last column (they are in `tab:judged_steps` otherwise).
- `--stand-in`: writes `results/stand_in/<name>.json` for the seven files, each with `"stand_in": true` and a note;
  `results/stand_in/` is now in `.gitignore`. A real file marked `stand_in` under `results/` stops the script.
- Every block that reads a stand-in or a quick-mode file starts its caption with "STAND-IN values" or "QUICK-mode
  values (from <file>)"; every phrase whose numbers come from such a file starts with "STAND-IN ". Nothing of the
  kind can reach the paper unnoticed.
- `claim(cond, text)` and `drift(cond, text)` replace the hard-coded assertions: a prose claim the data no longer
  supports prints `CLAIM FAILS`, a disagreement between sources prints `SOURCE DRIFT`; both fail `--check` and never
  stop generation, so the writers see the numbers and the failing claims together. Structural assertions (file
  shapes) still stop the script. No count is typed: the eleven, 55 pairs, 2,250 items, 150 templates, 480 readings,
  277 pairs and 316 returns are all derived.

## Blocks written (19; the registry's ten in bold)

| label | default file | source file(s) | state at 2026-10-07 20:20 local |
|---|---|---|---|
| **`tab:main_results`** | `6_results.tex` | `results.json`; `matched.json` for the reasoning rows or the matched headline | QUICK (matched.json) |
| **`tab:matched`** | `appendices/results.tex` | `matched.json`, `results.json` | QUICK |
| **`tab:judged_steps`** | `appendices/results.tex` | `results.json`; `matched.json` for the re-run rows | QUICK |
| **`tab:single_path`** | `appendices/results.tex` | `single_path.json` | QUICK |
| **`tab:coverage_variants`** | `appendices/results.tex` | `coverage_variants.json` | full resolution |
| **`tab:scoring_variants`** | `appendices/results.tex` | `sensitivity_variants.json` (store `main`, or `matched` under `--headline matched`) | QUICK |
| **`tab:providers`** | `appendices/results.tex` | `providers.json` | QUICK |
| **`tab:depth_model`** | `appendices/results.tex` | `depth_model.json` | QUICK |
| **`tab:flag_precision`** | `appendices/error_analysis.tex` | `flag_precision.json` | **STAND-IN** (WS-C3 has not written it) |
| **`tab:errors`** (extended) | `appendices/error_analysis.tex` | `scored_current.json` | current figures |
| `tab:level_gap`, `tab:coverage`, `tab:branch_domain`, `tab:paraphrase`, `tab:experiments` | as before | `results.json` and the Markdown records | as before (`tab:level_gap` loses a column under `--repaired`) |
| `fig:level_bars`, `fig:error_categories`, `fig:branch_bars`, `fig:domain_radar` | as before | `figure_data()` | the error-categories caption and figure now read the current B2 figures |

Table 1 as redesigned: Model; FAC with interval; Tier (compact letter display from `letters()`, no bold "best"); No
readable answer (share of the 2,250 responses, empty or without a stated answer); MC with the judge (interval); MC
by matching alone (interval); Judge-decided share; Arithmetic flags on correct answers (interval); Calculations
parsed per response; a `Formula flags` column appears automatically when `results.json` carries
`formula_flag_rate_on_fully_solved` (and `matched.json` `formula_flag_rate`); Judged step flags only with
`--judged-in-table`. Captions of the redesigned and new blocks are in sentence case, end with a full stop, and say
what the numbers are and are not; the untouched appendix tables keep their bold title-case lead-ins for Phase 2.

Phrases: every phrase of the previous version is kept with its numbers from the data, except these changes: the B2
sample sentence and the error-appendix sample phrases are generated from `scored_current.json` (140 wrong answers
read of a sample of 148; 23, 38, 40 and 39 per model; 420 readings), the "40 from each of four models" assertion is
gone; the repeat sentence names the four models (`results.json`'s `repeats`, the count checked against
`subsamples.repeat_ids()`); the branch spread prints three decimals; "Four experiments" (the number of distinct arms);
new phrases for the matched-settings sentence (changes per re-run model, top-tier sizes under both configurations,
Kendall's tau), the readable-only reordering (main text), the single-path pointer (`\autoref{tab:single_path}`), the
coverage-variant separation rule and verbosity slope, the depth model, the scoring-rule variants, the carried
precision, the full-set reasoning changes (experiments appendix), and "N of the nine templates with symbolic answers
are scored by the numbers they state" read from `answer.SYMBOLIC_EQUIVALENCE_TEMPLATES` (empty until ANSWER FINAL,
so today it says all nine).

## Phase 2 command lines (after WS-D1 and D2 place the markers; from the repository root)

```
python full_run_28092026/paper_results.py --stand-in          # only if a WS-C file is still missing; never commit results/stand_in/
python full_run_28092026/paper_results.py --write --headline <D2: default|matched> [--repaired] [--judged-in-table]
python full_run_28092026/paper_results.py --check
# while editing prose (WS-D1):
python full_run_28092026/paper_results.py --write --text-only --headline <D2> [--repaired]
python full_run_28092026/paper_results.py --check
# a dry run into a copy at any time:
python full_run_28092026/paper_results.py --out <scratch> --headline <D2> [--repaired]
python full_run_28092026/paper_results.py --check --out <scratch> --headline <D2> [--repaired]
```

`--write` without `--out` stops if any block has no markers in the tree and names the default file; the labels and
default files are the table above (the writers may put the markers in any `.tex` of the tree).

## Check outcomes

- `--out <scratch> --headline matched --repaired` and `--out <scratch> --headline default` both write all 19 blocks
  (11 into the copied files, 8 to `generated/`) and the phrases; the four figures are drawn into the copy. A
  column-count check over every generated tabular (spec against each row) finds no mismatch under either headline
  or with `--judged-in-table`. No LaTeX compiler is on this machine, so the blocks were not compiled.
- `--check --out <scratch>` on each copy: 4 failures, the rest pending Phase 2 (8 blocks in `generated/`, 28 phrases
  not yet in the prose, 18 prose numbers that no phrase generates any more, e.g. 480, 160, 40). The 4 failures are
  the two source drifts and the two failing claims below, which are real.
- `--check` on the untouched real tree: 59 failures: 4 STALE blocks (`tab:main_results`, `tab:errors`,
  `tab:level_gap` only under `--repaired`, `fig:error_categories`), 8 UNPLACED blocks (the new ones), 26 MISSING
  phrases (the changed and new ones), 17 NOT GENERATED numbers (the old sample sentence's), plus the same 2 drifts and
  2 claims. The brief's acceptance "passes on the untouched tree" cannot hold once Table 1 is redesigned and the B2
  sentence is regenerated from the current figures: the tree still holds the old Table 1 and the old sentence. Every
  item the check lists on the real tree is a Phase 2 placement or wording edit, and the check distinguishes them
  from genuine failures (citations, labels, long lines, figures, drift, claims), of which the real tree has only the
  four below.
- The `--check` of the previous version crashed on import (a hard-coded paraphrase count, 316 returns, which
  `PARAPHRASE.md` no longer holds); it runs now.

**Source drift found (fails every check until fixed, by design):**
1. `PARAPHRASE.md` says 314 paraphrases pass the checks; `PARAPHRASE_REVIEW.md` still says 316 were returned to the
   experts and 277 kept. `results.json`'s Q5 also still runs on 277 pairs. WS-B's and WS-A's figure is 275 kept
   (314 less 39 rejected). `PARAPHRASE_REVIEW.md` needs regenerating (`paraphrase_kit.py --score`, WS-E or WS-G) and
   `analyze.py`'s main pass re-running (WS-G).
2. `results.json` (4 October, before the round-5 re-run) and `matched.json` (today) disagree on the default-store
   FAC of 7 models. `--headline default` therefore prints the old eleven beside new reasoning rows until WS-G runs
   `analyze.py` without `--quick` (the drift line names the models).

**Claims the current B2 figures no longer support (printed as `CLAIM FAILS`; WS-D1 rewords):**
1. "no error is Claude Sonnet 5's largest majority label": with the readings on the two repaired templates excluded,
   Claude's majorities are calculation 10 and no error 9 (`scored_current.json`).
2. "the Advanced column's no-error share comes largely from Claude Sonnet 5": the lower bound the old check used
   (Claude's no-error readings minus its non-Advanced readings) is now negative, so the data cannot support it as
   stated. The generated phrase still carries the number (21%) for the writers to qualify.

## Schema against C1 to C3's files (all read as written; nothing blocking)

- `matched.json` (C1): the brief's keys plus the extras C1 lists. Used beyond the brief: `judge_decided_share`,
  `claims_per_trace`, `digit_flag_ci`, `router_flag_ci`, `setting`, `empty_n` (the dagger note under
  `--headline matched`), `tau_default_vs_matched.ci`. `unreadable` is a share; the loader also accepts a count
  (`share()` divides a value above 1 by the items). A `formula_flag_rate` is picked up when present.
- `single_path.json` (C1): as specified (model keys at the top level with `quick`, `n_single`, `n_others`).
- `providers.json` (C1): a list, `quick` on every row; `few_templates` and `few_matched` both mark the endpoint
  (`§`: served or matched fewer than 20 templates), as C1 advises; a `matched_diff` of null prints `--`. Only the
  eight models with more than one endpoint appear, and the caption says so.
- `coverage_variants.json` (C2): the brief's keys plus C2's extras; the verbosity slope's unit is read from
  `verbosity.units.numbers` ("MC per 10 numeric values shown"), the rule text and the counts of instances left
  without milestones from `rule`. Letters per variant are drawn by `letters()` from `pairs[<variant>]`.
- `sensitivity_variants.json` (C3): `stores` has `main`, `reasoning-medium-full` and `matched`; the table prints
  `main`, or `matched` under `--headline matched`. `relative_error_bins` (counts) printed as counts.
- `depth_model.json` (C3): `models[k].bins[b]` with `rate`, `ci`, `n`; `slope`, `ci`, `p_holm`; as specified.
- `flag_precision.json` (C3): **missing**; the stand-in is read and every block and phrase that uses it says
  STAND-IN. The expected schema is the brief's (`models[k]` with `flags`, `carried_precision`, `other`;
  `expert_confirmed` with `slips`, `carried_precision`, `other`); `drift()` checks `slips` against
  `FLAG_REVIEW_3.md`'s 171.

## Other decisions and notes

- `NO_REASONING` (the `∗` mark) is derived from `matched_config.json` (`reasoning_setting` other than "reasons by
  default"), `RERUN` from its `reasoning_store`; a `run: false` entry is skipped. Model order and the pair count come
  from the stores present (`results.json`'s Q1), never a hard-coded eleven; the twelfth model needs a display name in
  `NAME` (`qwen3-235b-a22b-thinking-2507` is there; add it to `paper_setup.py` if it ever runs).
- `tab:errors` gained, per level, the share with the no-error readings removed (`readings (share / share without
  no-error readings)`), the sample composition by level from `scored_current.json`'s `composition`, and the readings
  per column; its caption states sample 148, read 140, per-model counts.
- `fig:error_categories` is now drawn from the current B2 figures (different readings per model, so the caption
  gives the range), which is why the real tree's copy of its block reads as stale. WS-F redraws at integration.
- `6_results.tex` in the real tree was modified at 20:01 local today by another session (it was already modified in
  the session-start `git status`); not by this script, which never wrote to the real tree.
- For WS-G (amendment 3): `docs/appendix_evaluation.py` still hard-codes the B2 sample (150 and 100 around its lines
  177 to 182); not edited here.

## Open items

- WS-C3: write `results/flag_precision.json`; until then `tab:flag_precision` and its phrase are STAND-IN.
- WS-E or WS-G: regenerate `PARAPHRASE_REVIEW.md` (275 kept) so the drift clears; WS-G: `analyze.py` full pass so
  `results.json` matches the re-scored store (the second drift) and the Q5 pair count moves to 275.
- WS-D1: reword the two failing claims; place the markers for the new blocks and the phrases; the `%% FILL` numbers
  can be taken from `generated/phrases.txt` of a `--out` run.
- WS-F: `fig:error_categories` caption text changed (readings per model now vary); the figure's data is the current
  B2 figures.
- Optional: `tab:judged_steps` prints only the correct-answer flag rate for the re-run models (the other columns are
  not in `matched.json`); if the writers want the full row, C1 would add `router_steps_flagged_per_trace`,
  `router_rate_on_wrong` and its interval to `matched.json`'s model entries.
