# WS-D2: Appendices, bibliography, listings and the small corrections (Phase 2, after INTEGRATED)

## Mission

Write the prose around every regenerated appendix block, the new release appendix, the extended related work,
the certification appendix's round-5 account, the validation and scoring appendices' clarifications, the
taxonomy appendix's selection rule and citations, the template-listing swap, the bibliography, and every small
correction both reviews listed. Every number comes from `reports/NUMBERS_SHEET.md`; a sentence that waits on a
number carries a FILL marker.

Closes: review 1 C1 (release, dates), M2 (the certification sentence; symbolic and prescribed-digit rules
stated), M5 (text), M7 (appendix wording, citations), M8 (release policy), M9b, M9c, m1, m2, m3, m4, m5, m9, m12,
m13 (caption), §5 tables, citations and typography; review 2 W2 (text), W7 (caption), W9 (appendix), W11 (release,
dates), W13 (text), W14 (text), Q11, Q12, Q13, §5 comments on the appendices.

## Context (read `PLAN_CONTEXT.md` sections 1 to 3, then `reports/NUMBERS_SHEET.md`)

What the appendices must now carry, and why:
- Certification: five rounds; round 5 re-certified the two chemical templates after their questions were made
  to state the system, the equation form and the property data (one sentence, no history words); the planted
  defects test numeric correctness and not under-specification; 504 of the 510 hand checks compared a number;
  kappa is undefined where every verdict approves; the round-3 approvals despite two mismatched hand checks and
  the round-4 "extension" explained in one sentence each from `layer2/fixes_round4.md` and `RESULTS_round3.md`.
- Scoring: the symbolic rule as it now stands (equivalence on the enabled templates, by the numbers they state
  on the rest); the prescribed-digit policy and its relaxed sensitivity; the absolute-value and last-digit clauses
  with the counts of verdicts that depend on them; proportional partial credit as a sensitivity; the worked
  example placed; the judged step check at the provider's default settings.
- Validation: how the 15 study templates were chosen (E4); which settings were chosen on the labels and which
  were cross-fitted; the judge's unreturned calls (51 of 60, 111 of 120); the two F1 values are two measurements
  (recorded matches against live matching), in the caption in plain words; Claude Opus 4.5 in the panel against
  Opus 4.7 in the study, one sentence; the per-category readings with their new denominators; the judged step
  flags' precision beside the appendix table that now holds them; the PRM threshold sentence.
- Models and setup: both configurations and their settings; query dates; the provider table and its sentence;
  the twelfth model; the statistical protocol gains the matched family.
- Full results: the new tables (single path, matched, coverage variants, scoring variants, depth model, flag
  precision, judged steps) each with one paragraph saying what it shows and does not; the judge-swap shifts in the
  same paragraph as the top-tier MC gaps; the Holm note for the level-gap table's columns.
- Taxonomy: the selection rule in the same words as §3.1 (D10); the LLM panel's role stated once; the per-area
  table and the agreement table with their sentences; "30 templates per branch" instead of "perfectly balanced";
  "ART-DEIT" cited or named plainly; IEEE 100's role stated or dropped.
- Template examples: an Advanced listing that meets §3.1's criteria and an Easy listing whose question does not
  state the formula; the precision policy sentence.
- Error analysis: the readings after the top-up; the paragraph on templates whose wording does not pin the
  answer reduced to the six near-miss templates; the four read-off templates named; the error shares with
  no-error readings removed.
- Paraphrase: the four repeat models and their 300 instances named; `Mistral Large 3` cited; the new kept count
  if pairs were dropped; the margin's basis stated (it equals the detectable paired difference).
- Further experiments: "four experiments"; the repaired instances in or out of each experiment, as the sheet says.
- Release (new): what is released and when, what is withheld and why, query dates, the script behind each
  table, the contamination policy (D8).
- Extended related work (new, optional): the fifteen closest benchmarks from `RELATED_WORK_v3.md`'s appendix part.

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, `master` only.
- WS-D1 works at the same time on the main text. You own the files listed below and nothing else. `custom.bib`:
  add entries only under a `% --- added by WS-D2` marker at the end and remove or correct existing entries only
  for the citations the appendices use (the society entries, `piepho2004letter`, `Mistral Large 3`); never edit
  WS-D1's block. `paper_results.py` is WS-D1's in Phase 2: send caption changes for appendix blocks through your
  report (or the owner) rather than editing it.
- Writing rules: `overleaf_source_04102026/WRITING_RULES.md` and `docs/PAPER_PLAN_OCT2026.md` section 8;
  appendices give implementation detail in lists where lists read better; final state only; every number from a
  script whose `--check` passes.
- A sentence that waits on a number: `%% FILL: <what>` on the line above.
- No commits or pushes unless the owner asks. When asked: one short line, no body, push at once.
- Finish by writing `docs/mock_review_workstreams/reports/WS-D2_report.md` from the template at the end.

## Before you start

- `reports/NUMBERS_SHEET.md`, the decisions D5, D6, D8, D9, D10 and the E4 facts (`00_ORCHESTRATION.md`).
- Read every appendix file end to end with its regenerated blocks, `7_appendix.tex`, `WRITING_RULES.md` (the
  label list at its end, which you keep current), `docs/appendix_certification.py`, `docs/appendix_evaluation.py`,
  `docs/appendix_listings.py`, `docs/appendix_statistics.py` (what each checks), WS-A's, WS-E's and WS-F's
  reports, `docs/related_work_oct2026/RELATED_WORK_v3.md` (the appendix part).

## Files you own

`overleaf_source_04102026/appendices/*.tex` (prose; the generated blocks stay as WS-G wrote them),
`appendices/release.tex` and `appendices/related_work.tex` (new), `7_appendix.tex`, `custom.bib` under your marker
and the entries named above, `docs/appendix_listings.py` (the template choice), `docs/appendix_evaluation.py` and
`docs/appendix_certification.py` (the sentences they check, after WS-G updated their data), `WRITING_RULES.md`
(the label list).

Not yours: `main.tex` and every main-text file, `paper_results.py`, `paper_setup.py` (WS-D1; send the
`appendices/models.tex` prose it checks to WS-D1 through your report if it must change).

## Steps

### D2-1. Certification (`appendices/certification.tex`)

Five rounds in the "Results" paragraph; the round-5 sentence; the planted-defect scope sentence; "504 of the 510
hand checks compared a number"; the mismatch split into real and planted items if `layer2/RESULTS.md` gives it;
the round-3 and round-4 sentences; the agreement table prints "undefined (all approve)" for civil and industrial
(edit the table and the caption; the numbers are in `appendix_certification.py`'s data); the electrical kappa
0.365 sentence. `python docs/appendix_certification.py --check` passes.

### D2-2. Scoring (`appendices/scoring.tex`)

The symbolic rule paragraph; the prescribed-digit paragraph with a pointer to `tab:scoring_variants`; the clause
counts; the proportional-credit sentence; the judged step check's settings; place the worked example
(`\input{sections/appendices/worked_example}`) with one introducing sentence; `appendix_evaluation.py --check`
passes.

### D2-3. Validation (`appendices/validation.tex`)

The study-template selection sentence (E4; FILL if missing); the chosen-versus-cross-fitted sentence; the
footnote on the judge's unreturned calls; the caption's two-measurements sentence; the Opus 4.5 sentence; the PRM
threshold sentence; the readings paragraph with the new denominators (B1, B3 after the repair; B2 after the
top-up); the judged step flags' table (`tab:judged_steps`) introduced with its precision.

### D2-4. Models and setup (`appendices/models.tex`, `appendices/results.tex`)

The configurations paragraph (both runs; the settings; the endpoints with no setting; the twelfth model; query
dates); the provider table's paragraph (`tab:providers`): matched differences, the endpoints with few templates
named as such; the statistical protocol gains one bullet for the matched family. In `results.tex`: one paragraph
per new table (what it shows, what it does not); the judge-swap shifts beside the top-tier MC gaps; the Holm
note on the level-gap table; the single-path table's sentence; the matched family's paragraph (tiers under both
configurations; the paired changes).

### D2-5. Taxonomy and statistics (`appendices/taxonomy_content.tex`)

The domain-selection paragraph rewritten to the actual process (D10), in the same words §3.1 uses (WS-D1 sends
its sentence; align); the per-area table's sentence; the agreement table's sentence (kappa, majority
differences, formula-stated count; D7's outcome); the statistics caption; the data-source table's unclear
entries. `appendix_statistics.py --check` passes.

### D2-6. Template examples (`appendices/template_examples.tex`, `docs/appendix_listings.py`)

Choose an Advanced template that meets §3.1's criteria (several coupled quantities, an iterative or implicit
solve, six or more gold milestones; take the milestone counts from `scores/milestones.json` read-only) and an
Easy template whose question does not state the formula, from different branches; regenerate the listings;
rewrite the introductory paragraph; one sentence on the precision policy (which templates prescribe rounding and
why). `appendix_listings.py --check` passes.

### D2-7. Error analysis, paraphrase, further experiments, branch and domain

`appendices/error_analysis.tex`: the readings paragraph after the top-up (and the matched-configuration
readings if E3c ran, in their own paragraph with the configuration named); the "wording does not pin the answer"
paragraph reduced to the six near-miss templates; the four read-off templates named with the sensitivity;
the no-error-removed shares introduced. `appendices/paraphrase.tex`: the repeat models named; `Mistral Large 3`
cited; the kept count; the margin's basis. `appendices/further_experiments.tex`: four; the instance counts per
experiment after the repair. `appendices/branch_domain.tex`: the heatmap's paragraph replacing the radar's; the
grouped-bar figure's paragraph.

### D2-8. Release appendix (`appendices/release.tex`, label `appendix:release`)

Per D8: what is released (templates and generator, evaluator code, the evaluation set with its seed at
publication, the responses with scores and judge rulings, the paraphrase pairs' hashes), what is withheld and
why (names; the experts' per-rater files), query dates, the script behind each table (one list), the
contamination policy (fresh sets from a new seed; a held-out seed kept). Add its `\input` to `7_appendix.tex`
and tell WS-D1 the label for the Limitations and Ethics pointers.

### D2-9. Extended related work (optional, `appendices/related_work.tex`)

From `RELATED_WORK_v3.md`'s appendix part: the fifteen closest benchmarks and the further-related-work
paragraphs, in the paper's register; add the `\input`; tell WS-D1 the label for the §2 pointer.

### D2-10. Bibliography (`custom.bib`)

Replace `aiche2025constitution`, `asme2025vision` and `ieee2025standards` with curricular documents (the NCEES FE
exam specifications for the five disciplines, or the ABET program criteria per discipline), each entry checked by
hand per the rules (published version, `series`, `address`, no `month`; protect capitals only where needed);
keep ABET, ASCE BOK3 and IISE BOK; add `Mistral Large 3`; remove `piepho2004letter` if nothing cites it; add
every key the appendices use; report the new keys to WS-D1 for §3.1.

### D2-11. Small corrections

Across your files: "physically-validated" to "physically validated" if it occurs in a caption you own; matching
quotation marks inside `\ttfamily` boxes; captions in sentence case ending with a full stop; "verdict" for scorer
outcomes only; one spelling of every model name; the figure label note per D9 in the captions you send to WS-D1;
Table 13's Holm note; `WRITING_RULES.md`'s label list updated with every new label.

### D2-12. Checks

`docs/appendix_certification.py --check`, `docs/appendix_evaluation.py --check`, `docs/appendix_statistics.py
--check`, `docs/appendix_listings.py --check` all pass; the vocabulary grep returns nothing in your files; every
`\autoref` in your files points at an existing label.

## Acceptance

- Every appendix file rewritten per the steps; the new appendices input from `7_appendix.tex`; the four check
  scripts pass; FILL markers left are listed with what they need.

## Report template (`reports/WS-D2_report.md`)

```
# WS-D2 report
- Files rewritten; one line each on what changed.
- New appendices and their labels (for WS-D1's pointers).
- Bib entries added (keys) and removed; the keys WS-D1 needs for §3.1.
- Caption changes for WS-D1 to apply in paper_results.py (label, new caption text).
- Listings chosen (templates, why).
- FILL markers left (file, line, need).
- Check outcomes.
- Open items.
```
