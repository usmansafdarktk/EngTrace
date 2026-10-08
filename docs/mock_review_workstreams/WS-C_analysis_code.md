> **Split on 2026-10-07 into five parallel sessions, each three to four hours.** This file stays the shared
> context (mission, rules, methods, amendments); each session opens its own brief and owns only the files it
> names: `WS-C1_matched_family_and_exports.md` (analyze.py), `WS-C2_coverage_variants.md`
> (coverage_variants.py), `WS-C3_scoring_variants_precision_depth.md` (clause_variants.py, flag_precision.py,
> depth_model.py), `WS-C4_paper_machinery.md` (paper_results.py), `WS-C5_statistics_tables_and_worked_example.md`
> (docs/appendix_statistics.py, worked_example.py). Every script writes its own result file under
> `full_run_28092026/results/` and its own section under `results/sections/`, so no two sessions write one file.
> Every script has a `--quick` mode for development; the orchestrator runs the full-resolution pass afterwards.

# WS-C: Analysis code, tables and the matched-settings family (Phase 1)

## Mission

Turn every analysis the reviews ask for into committed, checked code, so that at integration one command
regenerates every table and figure block from the final stores: the single-path consistency table, Milestone
Coverage variants with their pairwise tests, the cross-model verbosity model, scoring-clause and
prescribed-digit sensitivities, proportional partial credit, carried-precision classification of arithmetic
flags, a controlled depth model, the matched-settings family, Table 1 redesigned, per-provider and
no-readable-answer tables, templates per area, the no-error-removed error shares, and a worked scoring example.
No tex is edited in Phase 1: scripts write into a scratch copy of the tex tree (`--out DIR`).

Closes: review 1 M1 (a to f), M3 (counts), M4 (readable-only, providers), M6 (depth), M9a, M9c, M9d, M9g, m6,
m7, m13; review 2 W3, W7 (table), W8, W10, W12, W13 (the single-path table), Q3, Q4, Q5, Q8, Q9, Q13, Q14.

## Context (read `PLAN_CONTEXT.md` sections 1 to 3 first)

What already exists and only needs printing (`full_run_28092026/results/RESULTS.md`, `results.json`):
- Q4 "Consistency within a template": per model, the share of the 58 single-path templates solved on all, some
  and no instances, and the same for the 92 others (`results.json["q4"]`). The main text cites it; no table
  prints it.
- "Milestone coverage over all traces" (D-170): E3 (matching alone) and E5-strict per model, all and readable,
  with template intervals.
- "By serving endpoint" (D-133): raw score, unusable share and a template-matched difference per endpoint.
- Q1: unusable (empty plus unreadable) per model. Sensitivity: the readable-only ordering (tau 0.636).
- The score store: `scores/<variant>/<model>.jsonl` rows carry `answer` (label, `matched`, `of`, `targets`
  with the numeric targets), `e3` (`coverage`, `reached` per milestone, `scale`, `null_coverage`), `steps`
  (claims, digit flags per step), `provider`, `reasoning_tokens`, `completion_tokens`, `billed_usd`.
  `scores/<variant>/e5/<model>.jsonl` rows carry per-milestone `sources` (`e3`, `REACHED`, `NOT_NEEDED`,
  `MISSING`, `UNJUDGED`), `e5_strict`, `e5_lenient`. `scores/<variant>/router/<model>.jsonl` rows carry the
  judged step flags. `scores/milestones.json` is the milestone cache (values per item).
- Matching and the arithmetic check read the response text only (`score.py` scores `row['text']`; the
  endpoint's reasoning text is stored separately by `run_traces.py`).
- `paper_results.py` already has `letters()` (compact letter display) and a `figure()` block writer; its drawing
  code is being moved by WS-F into `paper_figures.py` (wait for FIGURES EXTRACTED before editing
  `paper_results.py`).

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, branch `master` only. Never create, switch to, list or inspect
  another branch.
- Other sessions work in this same working tree at the same time. Edit only the files under "Files you own". If a
  step needs another stream's file, write the need into your report and continue.
- Paid API calls: none, except the optional step-check stability sample, cap **$3**, after a dry run.
- No commits or pushes unless the owner asks in this session. When asked: one short line, no body, push at once.
- Every script lives in the repository, has a `--check` or self-test, and prints into `results/RESULTS.md` under
  its own heading and into `results/results.json` under its own key; no typed numbers anywhere.
- Write new code freely within your files; refactor where it makes a step cleaner; keep every existing `--check`
  passing on the current store.
- No `.tex` file is edited in Phase 1. `paper_results.py --write` and `docs/appendix_statistics.py --write` go to
  a scratch copy via `--out DIR` until Phase 2.
- Finish by writing `docs/mock_review_workstreams/reports/WS-C_report.md` from the template at the end.

## Before you start

- Read: `full_run_28092026/analyze.py` (its docstring: the method choices; `stage_q3`, `q3_coverage`, the
  pairwise machinery, the sensitivity block, the endpoint block), `paper_results.py` (blocks, phrases, `letters`,
  `table`, `figure`, `check`), `ANALYSIS_PLAN.md`, `results/RESULTS.md` end to end, `score.py` (row fields),
  `judge.py` (row fields), `router.py` (row fields), `expert_kits.py` (what `--score` prints into
  `EXPERT_REQUEST.md`), `RESIDUAL_INCORRECT.md`, `FLAG_REVIEW_3.md`, `JUDGE_SWAP.md`, `docs/appendix_statistics.py`,
  `evaluation/milestones.py` (how milestones are derived; whether the final-answer line is marked).
- The registry of names in `00_ORCHESTRATION.md` section 7 (labels, flags, variant names). Use exactly those.
- Signals you depend on: FIGURES EXTRACTED (WS-F) before any edit of `paper_results.py`; FREEZE DONE (WS-A)
  changes the manifest (your scripts must read it fresh, never cache it); RUNS SCORED and
  `results/matched_config.json` (WS-B) for the matched family on real stores (until then, test it on the
  450-item `reasoning-medium` stores as a stand-in).

## Files you own

- `full_run_28092026/analyze.py`, `paper_results.py` (after FIGURES EXTRACTED; everything except the drawing
  module), `results/RESULTS.md`, `results/results.json`, `results/per_template.csv`.
- New: `full_run_28092026/coverage_variants.py`, `clause_variants.py`, `flag_precision.py`, `depth_model.py`,
  `worked_example.py`, `single_path_table` (inside `paper_results.py`), and their outputs under `results/`.
- `docs/appendix_statistics.py` (add `--write` for marked blocks and the new tables).
- The optional stability run's rows: `scores/main/router_stability/` (new).

Not yours: `answer.py`, `expert_kits.py`, `EXPERT_REQUEST.md` (WS-E; read only); `score.py`, `judge.py`,
`router.py` (run only); `freeze.py`, `manifest.jsonl` (WS-A; read only); `run_traces.py`, `models.json` (WS-B);
`paper_figures.py` and the PDFs (WS-F); every `.tex`.

## Steps

### C1. `analyze.py`: the matched-settings family and the exports (day 1, no `paper_results.py` yet)

1. A `matched` family: given `results/matched_config.json` (model → default store, reasoning store or null),
   build the configuration "each model at its reasoning-on store where one exists, else its default store", and
   compute on it everything Q1 and Q3 compute on the default configuration: per-model FAC with template
   intervals, the 55 (or 66 with a twelfth model) pairwise sign-flip tests with Holm, the letters, the share with
   no readable answer, the level means and gaps, MC with and without the judge and its pairwise family, the flag
   rates. Also the paired change per re-run model against its default rows (the full-set version of the
   reasoning-effort experiment: paired by item, template bootstrap, sign-flip, detectable difference, the
   90% interval against ±0.05). Write `RESULTS.md` section "Matched settings" and `results.json["matched"]`.
   Until RUNS SCORED, test with a stand-in config that points the two closed models at the 450-item
   `reasoning-medium` store (restricted to its items) and label the output as a stand-in.
2. Export the Q4 single-path numbers with intervals as a table-ready structure (they are in `results.json`
   already; make sure the "none" column and the 92-others columns are there).
3. Export the endpoint table (`reported.providers`) with the matched difference and the number of templates
   matched, and a flag for endpoints that served fewer than 20 templates.
4. Pairwise family on E3 coverage (matching alone), as Q3 does for E5-strict.
5. Drop from the judge-swap table the sampled responses that WS-A reports as sitting on the two repaired
   templates (read the count from WS-A's report when it arrives; until then, parameterize).
6. A twelfth model: `ORDER`-style lists must come from the stores present plus `matched_config.json`, not from a
   hard-coded eleven.

### C2. `coverage_variants.py` (new)

From the E3 rows and the judge rows of a store, per response:
- MC as scored (reproduce `e5_strict`; assert equality as the self-test);
- MC by matching alone (E3);
- route-adjusted MC: milestones ruled `NOT_NEEDED` leave the denominator (`(e3 + REACHED) / (n - NOT_NEEDED)`);
- intermediate-only MC: the answer-target milestones leave the set. Identify a target milestone by matching its
  value to the answer row's `targets.numbers` under the check's unit factors, or by the final-answer line if
  `evaluation/milestones.py` records which milestones the gold answer line states (read it; prefer the recorded
  rule). Report the rule used and how many instances are then left without milestones.
Per model: means with template bootstrap, for all responses and for wrong answers; then feed the 55-pair tests
(sign-flip, Holm) for each variant through `analyze.py`. Also:
- the cross-model verbosity model: response-level regression of MC on the number of numeric values the response
  displays (count them with the same reader the arithmetic check uses) and on visible completion tokens, with
  template fixed effects, pooled over models; report the slope with a template-cluster bootstrap interval; keep
  the within-model Spearman beside it;
- MC against reasoning tokens per model (descriptive: means by reasoning-token quartile within model).
Write `RESULTS.md` "Coverage variants" and `results.json["coverage_variants"]`; `--check` on a stand-in store.

### C3. `clause_variants.py` (new)

Re-apply the final-answer rule offline over the main store (and the matched stores) with:
(a) the absolute-value clause off; (b) the last-digit term bounded, min(u(ŷ), 0.01|y|); (c) the exact-digit
requirement of the 17 prescribing templates relaxed to the tolerance; (d) proportional partial credit
(`matched / of` instead of 0.5). Per model: FAC under each, the number of verdicts that change and in which
direction, and the distribution of relative error among accepted answers (bins ≤0.2%, 0.2 to 1%, 1 to 5%,
larger). `answer.py` is WS-E's: do not edit it; call `answer.verdict` with its parameters if they exist, else
implement the variant on the stored `targets` and the response text with a local copy of the rule, documented in
the docstring. Write `RESULTS.md` "Scoring-rule variants" (extend the Sensitivity table) and
`results.json["sensitivity"]`.

### C4. `flag_precision.py` (new)

Classify every arithmetic flag on a correct-answer response into "carried precision" (the printed result is a
correct rounding of the value recomputed from the unrounded upstream values that appear earlier in the response)
and "other", per model, using the arithmetic module's parser (`arith.py`) and the stored step claims. On the
flags the domain expert read (`FLAG_REVIEW_3.md`; the item list is local), report how many of the 171 confirmed
slips are carried-precision cases. Write `RESULTS.md` "Arithmetic flags: carried precision" and
`results.json["flag_precision"]`.

### C5. `depth_model.py` (new)

Wrong-answer rate against the gold trace's milestone count over readable responses only (empty responses
excluded), per model: per-bin means with template intervals, and a template-clustered logistic model of wrong
answer on milestone count with answer kind as a covariate (cluster-robust or bootstrap over templates). Report
the slope per model and whether it holds after Holm over models. Write `RESULTS.md` "Depth, controlled" and
`results.json["depth_model"]`.

### C6. Error shares with "no error" removed

From the counts `expert_kits.py --score` prints into `EXPERT_REQUEST.md` (read only; WS-E regenerates it after
the top-up), compute the by-level shares with the no-error readings removed, for the `tab:errors` block. Keep the
code reading the file, so the top-up flows through.

### C7. `paper_results.py`: Table 1 and the new blocks (after FIGURES EXTRACTED)

1. `--out DIR`: write every generated block and figure into a copy of `overleaf_source_04102026/` under DIR;
   `--check` unchanged. Never write to the real tree in Phase 1.
2. `--headline default|matched` and `--repaired` flags (registry, section 7 of `00_ORCHESTRATION.md`).
3. Table 1 (`tab:main_results`): FAC with interval; a tier letter from `letters()` (no bold "best");
   no-readable-answer share; MC with the judge; MC by matching alone; the judge-decided share; arithmetic flags
   on correct answers with "calculations parsed per response"; the judged step flags move to `tab:judged_steps`
   in the appendix (keep a `--judged-in-table` switch in case the writers want them back with a caveat). Rows per
   `--headline`: `default` prints the eleven at provider defaults plus a marked block "with reasoning at medium
   effort" for the re-run models; `matched` prints each model at its matched setting with a footnote mark on the
   re-run models and a twelfth row if present, and the default rows go to `tab:matched`.
4. New blocks: `tab:single_path`, `tab:matched` (the matched family: scores, letters, the paired changes),
   `tab:coverage_variants`, `tab:scoring_variants`, `tab:providers`, `tab:depth_model`, `tab:flag_precision`,
   `tab:judged_steps`, the extended `tab:errors`. Captions in sentence case, ending with a full stop, saying what
   the numbers are and are not.
5. The phrases section: keep every phrase generated from `results.json`; add phrases for the matched-settings
   sentence, the readable-only reordering, the single-path pointer, the four repeat models' names, the branch
   spread at three decimals, "four experiments"; remove the "with and without the two chemical templates"
   phrases under `--repaired`. Phase 2 (WS-D1) will edit wording; your job is that every number in a phrase comes
   from the data.
6. `--check` passes on the real tree (unchanged) and `--out` produces the full new set.

### C8. `docs/appendix_statistics.py`

Add `--write` that fills marked blocks (the same BEGIN/END convention) and these tables: `tab:area` (templates
per area and level from `manifest.jsonl`'s `area` field, grouped by domain), `tab:levels_agreement` (from
`template_annotation_23092026/levels/agreement.json` when present: kappa per branch and overall, templates where
the majority differs from the label, the formula-stated count), and the single-path counts if the writers want
them there. `--check` keeps passing; `--out DIR` as above.

### C9. `worked_example.py` (new)

Pick one real response from the main store that shows every mechanism (a lower-tier wrong answer with a judge
ruling and an arithmetic flag; or two short ones) and generate `appendices/worked_example.tex` (under `--out` in
Phase 1): the question in brief, the gold milestones with values, which matched deterministically (with the unit
factor), which the judge ruled reached, not needed or missing, the arithmetic flag with the operands and the
recomputed value, then v(r) and m(r). Every value from the store rows. Keep it to one column of a page.

### C10. Optional: step-check stability sample ($3 cap)

Re-run `router.py` twice on the 220 judge-swap responses at the provider's default settings into
`scores/main/router_stability/` and report the flag agreement across the reruns (per response and per step).
Dry run first; skip if the stage script cannot be pointed at a sample without editing it.

### C11. Run everything on the current store and record the numbers

Run every script end to end on the stores as they are today. Put the headline findings into your report so the
orchestrator can read them before integration: in particular whether Claude Sonnet 5 and DeepSeek V4.1 Flash
separate after Holm under matching alone and under route-adjusted MC (the rule for the results sentence), the
clause counts, the carried-precision split, the depth slopes, and the stand-in matched family.

## Acceptance

- Each new script: in the repository, `--check` or self-test passes, `RESULTS.md` section and `results.json` key
  written, docstring states the method and the choices.
- `paper_results.py --check` passes on the real (unchanged) tree; `--out DIR` writes every block listed in the
  registry; `--headline` and `--repaired` both work.
- `docs/appendix_statistics.py --check` passes; `--write --out DIR` fills `tab:area` (and
  `tab:levels_agreement` when the file exists).
- The twelfth model, if present in `matched_config.json`, flows through every table.
- No `.tex` under `overleaf_source_04102026/` or `current_overleaf_project/` modified (check with `git status`).

## Report template (`reports/WS-C_report.md`)

```
# WS-C report

## Scripts and flags
| script | what it computes | RESULTS.md section | results.json key | check |
- paper_results.py: flags added; blocks added (labels); the exact Phase 2 command lines (--write with flags).
- docs/appendix_statistics.py: --write blocks.

## Findings on the current store (numbers with their source)
- Claude vs DeepSeek on MC: as scored / matching alone / route-adjusted (p Holm each) -> rule outcome.
- Intermediate-only MC: instances left without milestones; per-model means.
- Verbosity model slope and interval; Spearman range.
- Clause counts per model (absolute-value, last-digit bound, prescribed digits, per-part credit).
- Carried precision: share of flags; of the 171 confirmed slips.
- Depth model: slopes per model, Holm outcome.
- Matched family (stand-in or real): tiers under both configurations.
- Single-path table: the "none" column range.

## Worked example
- Item chosen and why; file written under --out.

## Open items
- Anything that needs another stream (e.g., a parameter in answer.py), with the exact need.
```

## Amendments, 2026-10-07 (from the WS-A, WS-B and WS-E reports; read before C1)

1. **State of the stores.** FREEZE DONE, STORES UPDATED and RUNS SCORED are posted. The manifest is final
   (`f2c1dd2f...`); every store's `CONFIG.json` is current against it (the `dirty` flag only says the manifest is
   uncommitted). `results/matched_config.json` exists with 12 entries: three models carry a reasoning store
   (`reasoning-medium-full`: GPT-5.4 mini, Gemini 3.1 Flash-Lite, Gemma 4 26B), one entry is the Qwen Thinking
   sibling with `run: false` and no store, and it must be skipped everywhere; Qwen3-235B-2507 is "none offered".
   The cascade stores are `paraphrase-reasoning-medium` (275 pairs, judged) and `repeat1..3-reasoning-medium`
   (300 items, scored only, Gemini and Gemma).
2. **Numbers that changed.** Paraphrase kept pairs: 275 (not 277). Judge-swap sample: 9 of the 220 responses sat on
   the two repaired templates and are archived; 211 remain. Three virial item ids changed
   (`work_isothermal_virial#12`, `#19`, `#24` out; `#14`, `#15`, `#16` in). The set-aside `qwen3.8-27b` store covers
   2,220 items. The main run's unusable counts moved slightly (GLM-5.3 now 99 empty).
3. **Error-analysis sample.** `paper_results.py` asserts the as-returned B2 sample (480 readings, 40 per model;
   around its line 433) and words "40 from each of four models". WS-E now writes the current figures to
   `full_run_28092026/expert_request/scored_current.json` and to a new section of `EXPERT_REQUEST.md`, "The
   current figures" (B2: 148 items, Claude Sonnet 5 has 28 wrong answers in all; B1 and B3 with reduced
   denominators). Read the current figures, drop the assertion, and generate the sample sentence from the data.
   `docs/appendix_evaluation.py` hard-codes 150 and 100 (its lines 177 to 182): note it for WS-G, do not edit it.
4. **The symbolic check.** WS-E's equivalence check will move up to 218 verdicts on up to seven symbolic templates
   when enabled at ANSWER FINAL; the enabled list is in `full_run_28092026/symbolic/SYMBOLIC_CHECK.md`. Every
   "symbolic answers scored by the numbers they state" phrase and the symbolic sensitivity row must read that list
   and `symbolic/grades.json` (keys in WS-E's report) rather than assume nine by-numbers templates. Run nothing
   that re-scores a store; WS-G does the one re-score after ANSWER FINAL.
5. **A fourth deterministic diagnostic may arrive.** A formula-consistency check (the formula a response states,
   evaluated with the values it gives, against the number it prints) is under discussion as WS-H. Write the
   flag-rate code generically over the step-row fields (`digit_flags` today; a `formula_flags` field if it lands),
   so a column and a validation row can be added without reworking Table 1 or `analyze.py`.
6. **Figures.** If FIGURES EXTRACTED has not been posted when you reach C7, do the extraction yourself exactly as
   WS-F's brief describes in its step F1 (move the drawing section of `paper_results.py` into
   `paper_figures.py`, one import and call left behind, identical regeneration), post FIGURES EXTRACTED in your
   own report, and leave `paper_figures.py` to WS-F from then on.
7. **Trace reviews.** `trace_review.py` now recognises the three new variant families (patched by the
   orchestrator); `decoding_table.py` heads any `reasoning-` variant with the old "(C1)" note, which is cosmetic.
