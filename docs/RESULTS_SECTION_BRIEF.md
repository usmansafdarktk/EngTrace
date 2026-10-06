# Results and Error Analysis: structure brief (2026-10-05)

The agreed layout of Sections 5.3 and 5.4 and of the appendix sections that carry their detail. Files are in
`overleaf_source_04102026/`; every number, table and figure block comes from `full_run_28092026/paper_results.py`
(`--write` rewrites the generated blocks, `--check` must pass). Four figures are approved and no other figure is
placed.

**Update, 2026-10-05 (owner, D-193):** the level bars figure goes in the main text (5.3.3) and the branch bars figure
in the appendix, beside the domain radar; both are single-column figures drawn at the column width, the branch bars
as stacked bars. Where the text below says otherwise, this update holds.

**Update, 2026-10-06 (D-194):** the results text and its appendices are rewritten to the owner's concision rules; the
appendix holds six tables, listed in D-194. Where the text below lists more, D-194 holds.

**Update, 2026-10-06 (D-195):** the results and its appendices are rewritten around the process measures: Table 1
carries the answer measures beside the derivation measures (the per-level columns move to the appendix level table);
5.3.1 runs four paragraphs; 5.3.4 Reasoning, Equations, and Tools is added for the three experiments;
`appendices/conditions.tex` is renamed `further_experiments.tex` (label `appendix:experiments`). Where the text below
differs, D-195 holds.

## Main text

**5.1 Evaluated Model Suite, 5.2 Experimental Setup** (6_experiments.tex, other session): unchanged.

**5.3 Results** (6_results.tex), three numbered subsubsections in the May layout.

- **5.3.1 Overall Model Performance.** Table 1: FAC and MC per model with 95% intervals and tier letters. Paragraph 1:
  the top tier of five within 0.009 and not separable; GLM-5.3 next, with its 100 empty responses; the lower five
  separated from all six; 32 of 55 pairs differ; the four models without reasoning tokens and their reasoning-on
  scores on the subset; MC orders the models differently and separates one top pair; no trade-off read into the two
  orderings. Paragraph 2: one sentence each, with its appendix pointer, on coverage on wrong answers against the
  chance floor, the two diagnostics as flags with stated precision that do not rank models, the paraphrase bound on
  277 expert-kept pairs, and the flagship anchors, governing equations and code tool on the 450-instance subset.
- **5.3.2 Performance by Engineering Branch and Domain.** Figure (full width): `figs/branch-bars.pdf`, FAC by branch
  for the four representative models. One paragraph: branch means support no order of difficulty; what 30 templates
  per branch detect; the one branch pair that separates; thermodynamics the lowest domain mean for eight of eleven
  models; one sentence reading the figure; pointer to the appendix.
- **5.3.3 Performance by Difficulty Level.** No figure. One paragraph: every model lower on Advanced by 0.055 to
  0.199; significant for four of eleven as scored and for none without the two Advanced chemical templates; the
  chemical experts' reading of those templates; GLM-5.3's gap as an output-ceiling effect; failure growing with the
  depth of the gold derivation; pointer to the appendix level figure and tables.

**5.4 Error Analysis** (6_results.tex). Figure (column width): `figs/error-categories.pdf`. Paragraph 1: the expert
reading (160 wrong answers of four models, three readers, kappa 0.930; Claude Sonnet 5 mostly "no error", the other
three mostly calculation slips; formula or principle errors 18% on Intermediate and Advanced against 3% on Easy).
Paragraph 2: what the benchmark discriminates (122 of the top five's 171 incorrect verdicts are form or question; the
headroom; the limit of verification). Targets: about 950 words for 5.3 and 340 for 5.4.

## Appendix, in order after Experimental Setup

- **Full Results** (`appendices/results.tex`), no figures. Scores and Pairwise Comparisons; Branch, Level, Domain, and
  Answer Kind; Sensitivity; Consistency, Instance Variance, Depth, Repeats, and Tokens; Coverage and Diagnostics.
- **Difficulty-Level and Domain Performance** (`appendices/level_domain.tex`, new), the May shape: one paragraph that
  reads both figures, then `figs/level-bars.pdf` (full width) and `figs/domain_radar/domain-radar-labeled.pdf` (0.74 of
  the text width, its drawn size). The paragraph names the four models, says the level gap holds for the four
  lowest-scoring models and for none without the two chemical templates, and that domain means carry no interval and
  no test.
- **Paraphrase Test** (`appendices/paraphrase.tex`), no figure: design, scripted and expert checks, funnel table,
  results table, tau against its noise floor, paraphrase beside decoding noise.
- **Additional Conditions** (`appendices/conditions.tex`), no figure: the four designs, the paired table, the anchors.
- **Error Analysis** (`appendices/error_analysis.tex`): taxonomy and protocol, sample, agreement and tables, the
  templates whose wording does not pin the answer.

## Figures and representative models

| Figure | Where | File |
|---|---|---|
| FAC by branch, four models (stacked) | appendix, column width | `figs/branch-bars.pdf` |
| Error categories by model and level | 5.4, column | `figs/error-categories.pdf` |
| FAC by level, four models | 5.3.3, column width | `figs/level-bars.pdf` |
| FAC by domain, four models | appendix, 0.74 text width | `figs/domain_radar/domain-radar-labeled.pdf` |

The four representative models are DeepSeek V4.1 Flash (FAC leader, open-weights), Claude Sonnet 5 (MC leader,
closed), GPT-5.4 mini (closed, no reasoning tokens, in the error analysis and all conditions) and gpt-oss-20b
(lowest FAC, open-weights; labeled "GPT OSS 20B" in the figures). The same four appear in all three bar and radar
figures. `figs/level-gap.pdf`, `figs/coverage-wrong.pdf`, `figs/paraphrase.pdf` and the unlabeled radar are not
placed; their blocks leave the text, and the files stay until the authors decide.

## Rules that bind every sentence

The plan's do-not-state list (`docs/PAPER_PLAN_OCT2026.md`, 3.2): tiers not an order, no branch or domain order, no
trade-off between the measures, the level gap always with both counts, no "robust to paraphrasing" (a bound), the
two closed models beside the tool result. "Paraphrase" is the one term for that experiment. Captions open with a
bold lead-in and end with a full stop; every figure and table is read in the text; one appendix pointer per main-text
paragraph; `WRITING_RULES.md` for form.

## Figure placement in detail

| Label | File | Section and position | Environment and width | Caption lead-in |
|---|---|---|---|---|
| `fig:branch_bars` | `figs/branch-bars.pdf` | appendix, beside the domain radar | `figure`, `\columnwidth` | Final Answer Accuracy by Engineering Branch for Four Representative Models. |
| `fig:error_categories` | `figs/error-categories.pdf` | 5.4, block right after `\label{subsec:error_analysis}` (already there) | `figure`, `\columnwidth` | Error Categories of the Wrong Answers Read. |
| `fig:level_bars` | `figs/level-bars.pdf` | 5.3.3, block right after `\label{subsec:difficulty_level}`, before the paragraph | `figure`, `\columnwidth` | Final Answer Accuracy by Difficulty Level for Four Representative Models. |
| `fig:domain_radar` | `figs/domain_radar/domain-radar-labeled.pdf` | new appendix section, after `fig:level_bars` | `figure*`, `0.74\textwidth` | Final Answer Accuracy by Domain for Four Representative Models. |

Each block is generated: `% BEGIN GENERATED <label> (full_run_28092026/paper_results.py --write)` to
`% END GENERATED <label>`, with the figure environment, caption and label inside, so captions carry data-driven
counts (58, 58 and 34 templates; 30 templates per branch; the four lowest-scoring models) without being typed. Every
figure is read in the text in one sentence; each main-text paragraph keeps one appendix pointer. Blocks to remove:
`fig:level_gap` from 5.3.3, `fig:coverage_wrong` from Full Results, `fig:paraphrase` from the Paraphrase Test; the
sentences that read them lose their figure reference and keep their table reference.

## Steps on approval

These steps overwrite the committed files in place: the current text is the base, the generated blocks are
replaced by the script's output, and only `appendices/level_domain.tex` is new.

1. 6_results.tex: add the branch-bars block and reading sentence to 5.3.2; remove the gap-plot block from 5.3.3 and
   point its sentence to the appendix figure.
2. appendices: new `level_domain.tex` with the paragraph and two blocks; add it to `7_appendix.tex` after
   `results`; remove the coverage-floor block from `results.tex` and the paraphrase block from `paraphrase.tex`.
3. `paper_results.py`: move the three blocks to their files; drop the three unplaced blocks; `FIGURES` keeps only
   the four placed files.
4. Rewrite the generated text blocks without redrawing any PDF; wrap; `--check` passes.
5. Commit only on the owner's word, since `paper_results.py` also carries the other session's uncommitted edits.

## Prompt for the session that carries this out

```
Place the four approved figures and finish Sections 5.3 and 5.4 of our ARR paper as docs/RESULTS_SECTION_BRIEF.md
specifies. Repo: C:\Users\ayesha.gull01\EngTrace, branch master. Another session works in the same repository on
full_run_28092026/paper_results.py and on figs/level-gap.pdf, figs/coverage-wrong.pdf and figs/paraphrase.pdf:
never redraw, delete or restore those PDFs, change paper_results.py only by targeted replacements (never rewrite the
file), and never touch 6_experiments.tex, 6_experiments_results.tex or appendices/models.tex. Commit with one short
line and push only when I say so; fetch before each push. Be token-efficient: reuse what exists, check only what you
need.

RULES, in order:
1. The four figures are approved as drawn. Do not redraw them; do not place any other figure.
2. The main text keeps the May layout: 5.3.1 Overall Model Performance (Table 1), 5.3.2 Performance by Engineering
   Branch and Domain (figs/branch-bars.pdf), 5.3.3 Performance by Difficulty Level (no figure), 5.4 Error Analysis
   (figs/error-categories.pdf). Each subsubsection is one or two paragraphs; every paragraph leads with its finding
   and ends with one appendix pointer; the main text carries headline numbers only.
3. The appendix gains one section, "Difficulty-Level and Domain Performance", in the shape of the May appendix's
   "Branch and Domain-Level Performance": one paragraph that reads both figures, then figs/level-bars.pdf
   (figure*, \textwidth) and figs/domain_radar/domain-radar-labeled.pdf (figure*, 0.74\textwidth). It goes right
   after Full Results in appendices/7_appendix.tex.
4. Every number and every figure block comes from paper_results.py; nothing is typed. Rewrite the generated text
   blocks without drawing: write_blocks() also redraws every FIGURES entry, so use a blocks-only rewrite (or add a
   --text-only flag) and leave the PDFs alone.
5. Keep every interpretation of the plan (docs/PAPER_PLAN_OCT2026.md, 3.2 and 8): tiers not an order; the four
   models without reasoning tokens beside any closed-against-open or level sentence; the level gap with both counts
   (4 of 11; 0 of 11 without the two Advanced chemical templates); no branch or domain order; "paraphrase", never
   "rewording"; a bound of ±5 points, never "robust"; the two closed models beside the tool result.

WHAT TO OVERWRITE: the committed 6_results.tex, appendices/results.tex, paraphrase.tex and 7_appendix.tex are
the base and are changed in place; the generated blocks are rewritten by the script, replacing their contents;
the only new file is appendices/level_domain.tex. Make no copies, no alternative versions and no new sections
beyond that one; the three removed figure blocks are deleted from the text, while their PDF files stay in figs/.

READ FIRST: docs/RESULTS_SECTION_BRIEF.md; overleaf_source_04102026/WRITING_RULES.md; 6_results.tex and
appendices/results.tex, paraphrase.tex, conditions.tex, error_analysis.tex as committed; the generated-block
mechanism and the phrase lists in paper_results.py (--check verifies blocks, phrases, every number in the prose,
citations, labels, figure files).

DO:
- 6_results.tex: add the fig:branch_bars block after \label{subsec:branch_domain} and one sentence reading it in the
  paragraph; remove the fig:level_gap block from 5.3.3 and point that paragraph's figure reference to the appendix
  level figure.
- appendices/level_domain.tex (new): the section, its paragraph, the two blocks; add \input{sections/appendices/level_domain}
  to 7_appendix.tex after results. Remove the fig:coverage_wrong block from results.tex and the fig:paraphrase block
  from paraphrase.tex, keeping their table references.
- paper_results.py: move blocks["main"]["fig:level_bars"] to the new file's key, drop the blocks for fig:level_gap,
  fig:coverage_wrong and fig:paraphrase, keep FIGURES to the four placed files (the drawing functions of the other
  figures may stay until the other session is done).
- Wrap prose at 100 characters, one sentence per line; run paper_results.py --check (exit 0) and
  docs/check_plan_claims.py (123 of 123).

REPORT: word counts for 5.3 and 5.4; each figure's label, file and section; what left the text; the check results;
anything you could not verify. Then wait for my word before committing.
```
