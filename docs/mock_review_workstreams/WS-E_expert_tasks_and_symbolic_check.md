# WS-E: Expert tasks and the symbolic check (Phase 1)

## Mission

Build, score and integrate the expert tasks this cycle needs: independent difficulty ratings (E1), grading of
the incorrect symbolic verdicts (E2), the error-analysis top-up after the two chemical templates are repaired
(E3b), and, if the matched configuration becomes the headline, the error-analysis re-read for the re-run models
(E3c). Implement a symbolic equivalence check in the final-answer rule, validated against E2, and the script
that adopts majority difficulty labels if the owner decides so. The owner dispatches kits to the experts; you
build and score them. No tex is edited.

Closes: review 1 M2 (symbolic answers), M6 (labels, method hints), m10, Q2, Q9 (part), Q15; review 2 W1
(labels), W2 (symbolic), Q6, Q14 (the counts feed WS-C).

## Context (read `PLAN_CONTEXT.md` sections 1 to 3 first)

- Difficulty labels: §3.1 says domain experts assigned Easy, Intermediate or Advanced along conceptual
  complexity, mathematical sophistication and procedural depth; 58 / 58 / 34 templates. No protocol or agreement
  is reported, the gap survives correction for none of the top models, the showcased Advanced listing is a
  two-formula problem, and one reviewer suspects Easy questions state the formula (the showcased Easy question
  does).
- Symbolic answers: nine templates have symbolic answers and the check scores them "by the numbers they state";
  `ber_estimation_mary` alone yields 34 of the top five's 171 incorrect verdicts, 68 across the eleven models;
  `RESIDUAL_INCORRECT.md`'s "symbolic" column counts 90 incorrect symbolic verdicts in all. The experts' B1
  readings called 12 of 98 incorrect verdicts correct.
- Error analysis: three readers read 160 wrong answers (40 per model: Claude Sonnet 5, GPT-5.4 mini, Gemma 4 26B,
  gpt-oss-20b) from the default-setting run (`EXPERT_REQUEST.md` B2; kits by `expert_kits.py --build`, scored by
  `--score`). For Claude, 25 of its 26 "no error" majorities fall on three contested templates (the two chemical
  ones and the symbolic one). After WS-A's repair, the wrong answers on the two chemical templates no longer
  exist; after a symbolic fix, many of the symbolic "wrong answers" become correct verdicts.
- The scoring provenance: `score.py` records the hash of `answer.py` in every store's `CONFIG.json`; any change
  to `answer.py` means every store is re-scored once at integration (WS-G). So finish `answer.py` early and post
  ANSWER FINAL.
- Where levels live: `freeze.py` calls `difficulty_map()`; find its source (`grep -rn "def difficulty_map"`).

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, branch `master` only. Never create, switch to, list or inspect
  another branch.
- Other sessions work in this same working tree at the same time. Edit only the files under "Files you own". If a
  step needs another stream's file, write the need into your report and continue.
- Paid API calls: none in this stream.
- No commits or pushes unless the owner asks in this session. When asked: one short line, no body, push at once.
- The experts' filled files stay local and gitignored (follow the existing pattern: `experts_filled_*` folders);
  reports and `EXPERT_REQUEST.md` carry counts only; no expert names anywhere.
- Scratch files go to the session's scratchpad; scripts and kit builders live in the repository.
- Write new code freely within your files. No `.tex` file is edited in Phase 1.
- Finish by writing `docs/mock_review_workstreams/reports/WS-E_report.md` from the template at the end. Post the
  signals at its top as soon as they are true.

## Before you start

- Read: `full_run_28092026/expert_kits.py` (`draw_b2`, `build`, `composition`, `--score`), `reading_app.py`,
  `EXPERT_READING_GUIDE.md`, `EXPERT_REQUEST.md`, `RESIDUAL_INCORRECT.md`, `residual_incorrect.py`, `answer.py`
  (the reader, the six answer kinds, how symbolic answers are scored today), `validate_scorer.py`,
  `SCORER_VALIDATION.md`, `score.py` (`EVALUATOR_FILES`), `template_annotation_23092026/layer2/app.py`,
  `workbooks.py`, `make_kits.py` (kit patterns), `docs/appendix_statistics.py` (how levels are counted),
  `freeze.py` (`difficulty_map`).
- The decision D5 (symbolic) and D7 (labels) in `00_ORCHESTRATION.md`; D2 decides whether E3c runs.
- Signals you depend on: STORES UPDATED (WS-A) before E3b; RUNS SCORED (WS-B) before E3c.

## Files you own

- `template_annotation_23092026/levels/` (new): `build_kit.py`, `score_levels.py`, `apply_labels.py`,
  `agreement.json`, `RESULTS.md`, `dist/`, `experts_filled_levels/` (local).
- `full_run_28092026/symbolic/` (new): `build_kit.py`, `score_grades.py`, `validate.py`, `grades.json`,
  `GRADES.md`, `SYMBOLIC_CHECK.md`, `dist/`, `experts_filled_symbolic/` (local).
- `full_run_28092026/answer.py`, `validate_scorer.py`.
- `full_run_28092026/expert_kits.py`, `EXPERT_REQUEST.md`, `expert_request/` (kits and returns, local).
- The level source that `difficulty_map()` reads (only through `apply_labels.py`, and only when D7 says so).

Not yours: `score.py`, `judge.py`, `router.py` (run only); the stores (you read them; WS-A and WS-B write them);
`freeze.py`, `manifest.jsonl` (WS-A); `analyze.py`, `paper_results.py` (WS-C); every `.tex`.

## Steps

### E1. Difficulty ratings kit

1. `levels/build_kit.py`: per branch, the 30 templates in a random order under opaque codes, each with one
   evaluation-set instance (question and gold solution) and the rubric of §3.1 (the three dimensions, one line
   each, with the level definitions). The expert picks Easy, Intermediate or Advanced, and ticks one box:
   "the question states the governing formula or names the method". The current label is not shown. Reuse the
   layer-2 kit format (a workbook or the small app) so the experts see a familiar interface. Write `dist/` per
   branch for the owner to send.
2. `levels/score_levels.py`: Fleiss' kappa per branch and overall on the three levels (and the weighted variant,
   since the levels are ordered), the majority label per template, the templates where the majority differs from
   the current label (and whether by one step or two), the formula-stated count per branch and level. Write
   `agreement.json` and `levels/RESULTS.md` (counts only). Post **KAPPA IN**.
3. `levels/apply_labels.py` (run only if D7 adopts the majority): update the level source that `difficulty_map()`
   reads, print the diff (template, old, new), and state in its output that `freeze.py` must be re-run by WS-G to
   refresh the manifest's `level` field and `FREEZE.json`'s `by_level` (instances do not change; only the label).
   Self-test on a copy.

### E2. Symbolic grading kit

1. `symbolic/build_kit.py`: every incorrect verdict on a symbolic template across the eleven models (the 90 rows
   of `RESIDUAL_INCORRECT.md`'s symbolic column; derive them from the store, not from the Markdown), grouped by
   template: the question, the gold answer line (from `pool/`), the response's final-answer segment as the check
   reads it (use `answer.py`'s reader so the expert sees what the check saw), under opaque codes, model hidden.
   The expert marks: equivalent to the gold, not equivalent, or unreadable, with an optional note. One kit per
   branch that has symbolic templates. Write `dist/`.
2. `symbolic/score_grades.py`: per template and per model, the counts; `grades.json` keyed by item id and model;
   `GRADES.md`. Post **GRADES IN**.

### E3. The symbolic equivalence check in `answer.py`

1. Read `ber_estimation_mary` (and the other eight symbolic templates) to learn the forms the gold answers take
   (functions such as Q, erfc, sqrt, log2; parameters such as M, Eb/N0). Design a parser for the response's
   final expression: LaTeX (`\frac`, `\sqrt`, `\log_2`, `Q(\cdot)`, `\operatorname{erfc}`) and plain text, into
   SymPy, with a small table of function aliases (Q in terms of erfc, and the like). Equivalence: `simplify`
   of the difference, then numeric sampling over the free symbols in their physical ranges; equivalent if both
   agree, or if sampling agrees within 1e-6 at twenty points when `simplify` is inconclusive.
2. Apply it only to templates where `symbolic/validate.py` shows it agrees with the expert grades (report
   precision and recall per template on the graded items; a template is enabled when both are at least 0.9 and
   every disagreement is explained). Keep the by-numbers rule for the rest, through an explicit list in
   `answer.py` (`SYMBOLIC_EQUIVALENCE_TEMPLATES`), documented.
3. Extend `validate_scorer.py` so the symbolic branch is covered by its self-test (gold answers of the enabled
   templates score correct; the graded items score as graded).
4. Write `symbolic/SYMBOLIC_CHECK.md`: the rule, the enabled templates, the validation numbers, what remains
   by-numbers and why. Post **ANSWER FINAL**. If the parser cannot reach the bar by the end of day 3, post
   ANSWER FINAL with no change to `answer.py` and say so: the expert-adjudicated score then goes in as a
   sensitivity (WS-C reads `grades.json`).

### E3b. Error-analysis top-up (after STORES UPDATED)

1. In `expert_kits.py` add `--b2-topup`: for each of the four error-analysis models, drop the read items that sit
   on the two repaired templates (and, if E3 enabled the symbolic check for a template, the read items whose
   verdict is now correct); draw replacements from the re-scored main store among wrong answers not yet read,
   keeping the level proportion rule where possible; if a model now has fewer than 40 wrong answers in all, take
   all of them. Three readers per item, the same six-category hierarchy and app. Write the kits; the owner
   dispatches.
2. `--score` merges the returned readings with the existing ones (the dropped items excluded), recomputes Fleiss'
   kappa per model and overall, and rewrites the B2 tables in `EXPERT_REQUEST.md`. Also recompute B1 and B3
   excluding the items on the two repaired templates (they were read against the old questions) and state the
   new denominators.
3. If D2 makes the matched configuration the headline, **E3c**: build B2 kits for the re-run models that are in
   the error analysis (GPT-5.4 mini; Gemma 4 if it ran) from WS-B's reasoning-on store (40 wrong answers each
   or all if fewer), three readers; score into a separate B2 section "matched configuration" so the writers can
   choose. The judge-ruling readings (B3) stay as they are and the paper says they are of the default-setting
   responses.

### E4. Hand the counts over

Everything the paper will cite from this stream must be in `agreement.json`, `grades.json`,
`SYMBOLIC_CHECK.md` and `EXPERT_REQUEST.md`, so that WS-C's and WS-G's scripts read them; nothing lives only in
your report.

## Acceptance

- E1 kits built for five branches; `score_levels.py` self-test passes on simulated returns (reuse the pattern of
  `layer2/simulate.py`); after the returns: `agreement.json` and `levels/RESULTS.md` written.
- E2 kits built; `score_grades.py` self-test passes; after the returns: `grades.json`, `GRADES.md`.
- `answer.py`: the equivalence check enabled only for validated templates; `validate_scorer.py` passes; every
  gold answer of the 2,250 still scores correct (run `score.py --check` or its equivalent self-test; do not
  re-score the stores yourself, WS-G does).
- `expert_kits.py --b2-topup` builds; `--score` merges; `EXPERT_REQUEST.md` regenerated with the new counts.
- `apply_labels.py` written and self-tested, run only on D7.

## Signals

`SIGNAL: KAPPA IN <date time>`; `SIGNAL: GRADES IN <date time>`; `SIGNAL: ANSWER FINAL <date time>`.

## Report template (`reports/WS-E_report.md`)

```
SIGNAL: KAPPA IN <date time>
SIGNAL: GRADES IN <date time>
SIGNAL: ANSWER FINAL <date time>

# WS-E report

## E1 difficulty ratings
- Kits sent (date), returned (date, how many of 15).
- Fleiss' kappa per branch and overall (unweighted / weighted); templates where the majority differs (count; by one step / two); formula-stated count per level.
- D7 recommendation from the data; apply_labels.py run or not.

## E2 symbolic grading
- Items graded (per template, per model); equivalent / not / unreadable counts.

## The symbolic check
- Enabled templates; validation precision and recall per template; what stays by-numbers and why; answer.py changed (yes/no) -> WS-G must re-score every store.

## E3b top-up (and E3c if run)
- Items dropped per model and why; items added; readers returned (date); kappa per model; new B1/B3 denominators.

## Files other streams read
- Paths and keys.

## Open items
```

## Amendments, 2026-10-07 (evening): evaluator defects to settle before ANSWER FINAL

The five WS-C sessions found three defects in the evaluator's readers and one rule the reviews attacked with
reason. All four belong to `answer.py` or the milestone reader, so they are yours, and every one of them is
covered by the single re-score WS-G runs after ANSWER FINAL. Do them before posting that signal; re-run the
expert-study validation (`validate_scorer.py`, `SCORER_VALIDATION.md`) afterwards, since the agreement figures were
measured with the readers as they were.

1. **Unicode superscript exponents in the milestone reader** (WS-C5, open item 1). `milestones.numbers` folds
   `x10^{-19}`, `x 10^-19` and the LaTeX `\times 10^{-19}` into one number but reads the unicode form
   `1.152 x 10^-19` written with superscript digits as 1.152 and 10; `answer.values` and `arith.normalise` fold
   that form. So matching misses milestones stated this way (19 to 24% of Claude Sonnet 5, GLM-5.3 and
   GLM-5.3-Flash responses use it; 0 to 1.5% of DeepSeek, Gemma, Gemini and GPT-5.4 mini), and can match a bare
   mantissa under a unit factor. Fold the form the way the other two readers do. Size: at most 70 responses of
   GLM-5.3-Flash and 17 of Claude change E3, and up to 183 judge jobs over the roster get a shorter
   missed-milestone list (new calls at the judge's rate, about $1; WS-G runs them). Decision D11a in
   `00_ORCHESTRATION.md`; the recommendation is to fix it, because it biases "matching alone" against the two
   models that write exponents this way, which is the comparison the results sentence rests on.
2. **A unit exponent read as a number** (WS-C3, "what drives the accepted answers above 5%"). `answer.values`
   drops a unit exponent only when a letter precedes it, so the `2` of `(signal units)^2` or `\text{(...)}^2` is
   read as a value, and `match` then credits a target through the percent factor with a one-unit window of
   plus or minus 100. Fix the reader. Decision D11b; the recommendation is to fix it (a defect, not a rule).
3. **A number on a non-final segment matched at a unit factor** when a response has no final-answer heading
   (C3's aliasing example: 0.5 x 10^3 within 100 against 540). Decide with C3's sample in hand whether the
   reader's fallback window is too wide; fix or document. Decision D11c.
4. **The last-digit window** (review 1 M3; C3's `last_digit_bounded` variant). The rule accepts a difference up
   to one unit of the response's last digit, so "2" matches 1.01. Bounding the term by min(u(y_hat), 0.01|y|)
   removes 121 main-run verdicts (tau with the headline 0.964) and every acceptance above 5% relative error that
   rests on a coarse last digit. Decision D11d: adopt the bound as the rule (recommended; it is the reviewers'
   own proposal, and the re-score is free) or keep the rule and report the variant. If adopted, state it in the
   scoring appendix's match rule and re-validate.
5. **Names C3's script expects at ANSWER FINAL**: `answer.symbolic_equivalence.apply(label, item, text, enabled,
   tol, unit, segment)` and `answer.SYMBOLIC_EQUIVALENCE_TEMPLATES`. If you use other names, tell the orchestrator
   so `clause_variants.py`'s `symbolic_hook()` is adapted before the full pass.
6. **E3c** is needed if D2 is the matched headline: GPT-5.4 mini has 72 wrong answers with reasoning on
   (`reasoning-medium-full`), Gemma 4 about 150; build the kits (`--b2-matched reasoning-medium-full`) as soon as
   D2 is recorded.
7. **Optional, E6**: the absolute-value clause carries 434 main-run verdicts on 22 templates (C3). If the paper
   is to call the clause a sign-convention allowance, one reader per item on a sample of about 40 (stratified by
   template) must say convention or wrong sign. Build only if the owner asks.
