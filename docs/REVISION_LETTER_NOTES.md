# Notes for the revision letter (October 2026 submission)

Split out of `docs/PAPER_PLAN_OCT2026.md` on 2026-10-03. Low priority until the paper is written; the letter is
drafted after it. The letter is the only place where "what changed and why" is said; the paper describes EngTrace
as it is.

## Spine: one row per reviewer point

The letter's body is a table with one row per numbered point from both cycles (January: gFWV, nWW3, cqGs and the
AC; July: yAYU, 9W1B, ynoK and the meta-review's three suggested revisions): the point in one line, where the paper
now answers it (section or appendix), and the answer in one or two sentences with its number. The record for the
table is `docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` section 4 (concern, who, where answered, state) and
`docs/NEXT_CYCLE_REVIEW.md` section 1 (read with its 2026-10-03 correction on authorship). The rebuttal texts are
`docs/EngTrace Rebuttals Jan 2026.docx` and `docs/EngTrace_Rebuttal_Jul2026 (1).docx`; the meta-reviews
`docs/meta-reviews.md`.

## 1. Why the numbers moved and cannot be compared

The evaluator the May paper used was measured against 15 experts' labels on 300 traces and found wrong on 72 of
them (68 correct answers called wrong), ranking the strongest model fourth where the experts rank it first; its
judge panel was two judges where the paper said three; its 20% wrong-answer sample reordered the models between
identical runs; its routing showed a judge the flawed step as often as a clean one (pilot `RESULTS_X1.md`,
`FINDINGS.md`). It was replaced, not patched. The May evaluation set is not reproducible: eight published templates'
gold traces did not reproduce their own answers, 77 templates moved during the rework, the generations are gone, and
the roster changed (`docs/re-implementation-sep/`, D-105, D-111). So no October number is placed beside a May number.
`full_run_28092026/RESULTS_PAPER_NOTES.md`, "Metric continuity", says what became of each May column (Final Answer
Accuracy, Reasoning F1, BERTScore and ROUGE).

## 2. Promises kept, in their new form

- Variance and significance (yAYU 2, 9W1B 4): per-template SD, template-level intervals, 55 Holm-corrected pairs,
  Welch's test on the level gap, detectable differences beside every null (`results/RESULTS.md` Q1, Q2).
- Wilcoxon on the continuous reasoning score beside McNemar on the answer (the July rebuttal's own promise): the
  coverage comparison over the 55 pairs (Q3, "Coverage compared across models").
- The judge-exclusion ablation with placebo and untouched controls (meta-review 1; yAYU 1; 9W1B 2; gFWV 4): replayed
  offline on the panel design, no family effect at a detectable 0.003 to 0.025 (pilot `RESULTS_LOJO.md`, D-174); plus
  the per-judge bias (uniform leniency), the judge swap on the new evaluator (`JUDGE_SWAP.md`, D-181) and a judge from
  no roster family.
- Threshold sensitivity (yAYU 4, 9W1B 3): `THRESHOLD_APPENDIX.md`; the cross-encoder and alignment-ratio thresholds
  the rebuttal promised to vary no longer exist in the evaluator.
- Inter-judge agreement (9W1B): the May panel 0.725 / 0.776, the replacement panel 0.563 / 0.622 (pilot `RESULTS_E1.md`).
- A fourth branch (ynoK 1; January AC): civil and industrial, 30 templates each, certified.
- Linguistic diversity and template exploitation (meta-review 3; 9W1B 1; cqGs 2; nWW3 2): the paraphrase arm as a
  bound of ±5 points for ten of eleven models (`PARAPHRASE_PAPER_NOTES.md`); the surface-shortcut audit
  (`SHORTCUT_AUDIT.md`).
- Tool or retrieval baselines (cqGs 1 and 4; ynoK 3; January AC): the open-book condition and the tool condition, in
  whatever state they are at submission (D-183, D-184).
- Error analysis too coarse, causes of the cliff (nWW3 3 and 4; cqGs 3): attribution validated by error type on the
  experts' labels (pilot `RESULTS_ATTRIBUTION.md`), the human reading of 160 wrong answers (`EXPERT_REQUEST.md` B2),
  failure against derivation depth, the two-template sensitivity of the level gap.
- Validation in the main text (9W1B 6; meta-review 2): the evaluation framework's validation subsection.
- Proofreading and the broken reference: moot, the text is new.

## 3. Promises replaced with a reason

- Branch-level reporting for the math-specialised models and the math-pretraining claim (meta-review 2; yAYU 5;
  9W1B 5): the claim is dropped; no math-specialised model is on the roster (the roster rule and the budget set it,
  D-110, D-117). A claim carried without the models would be worse than no claim.
- An interval on rho = 0.632 and DeepSeek V4 Pro's category: those numbers no longer exist; the 100-response study is
  replaced by the 300-trace study with 15 experts.

## 4. Objections answered by construction rather than by argument

- gFWV 1 (perfect kappa unconvincing): the certification rejected 53 real templates over its rounds, caught 60 of 60
  planted defects, and records hand checks and timestamps.
- gFWV 3 (instance-level verification shallow): the integrity checks at 500 seeds, the validation of every
  evaluation-set item's gold solution, the display-tie census.
- gFWV 5 (authorship): stated exactly for each branch. The 90 chemical, electrical and mechanical templates were
  written by the team's domain experts; the 60 civil and industrial templates by a colleague of the authors; all 150
  passed the same certification. (The September plan's note that the 60 were "AI-drafted" was a misreading, corrected
  2026-10-03; nothing of the kind goes in the letter.)
- gFWV 4 and 9W1B 2 (validated by the same kind of system): expert labels, deterministic components, a judge from
  outside the roster.
- nWW3 1 and 2 (comprehensiveness; 90 templates are 90 problems): five branches, 150 templates, structural variation
  measured (72 templates change the governing equation, step count or computed quantities with a sampled parameter),
  the evaluation set's coverage of reasoning paths.
- cqGs 3 (causes of the cliff): depth, attribution, the two-template sensitivity, the human reading by level.

## 5. What the letter must disclose before a reader of the repository finds it

The September audit's findings on the published templates; the evaluator's measured defects; the Markdown rendering
during certification (65 templates judged through a renderer that dropped `*` and `$`); the answer-check corrections
after the run (238 verdicts moved, 234 to correct, every one read; one known wrong credit); and that the new numbers
are high because the models of September 2026 solve the evaluation set at the top on final answers.
