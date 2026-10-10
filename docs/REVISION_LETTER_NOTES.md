# Notes for the revision letter (October 2026 submission)

Split out of `docs/PAPER_PLAN_OCT2026.md` on 2026-10-03. Low priority until the paper is written; the letter is
drafted after it. The letter is the only place where "what changed and why" is said; the paper describes EngTrace
as it is.

Updated 2026-10-10 for the October revision (section 7). The paper's headline is the matched configuration: every
model at its reasoning setting where the endpoint offers one (GPT-5.4 mini, Gemini 3.1 Flash-Lite and Gemma 4 26B at
medium effort; Qwen3-235B-2507's endpoint offers none), the three models' default scores beside as the comparison.
The expert readings (error analysis, flag precision, judge swap) and the further experiments are of the default
responses and the paper says so. Every figure below is from `docs/mock_review_workstreams/reports/NUMBERS_SHEET.md`
(each line there names its source file); cost figures stay out of the letter.

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
  Welch's test on the level gap, detectable differences beside every null (`results/RESULTS.md` Q1, Q2). Now: 36 of
  55 FAC pairs separated after Holm at matched settings (38 at the providers' defaults); decoding repeats for four
  models (300 instances, three repeats each; SD at most 0.008).
- Wilcoxon on the continuous reasoning score beside McNemar on the answer (the July rebuttal's own promise): the
  coverage comparison over the 55 pairs (Q3, "Coverage compared across models").
- The judge-exclusion ablation with placebo and untouched controls (meta-review 1; yAYU 1; 9W1B 2; gFWV 4): replayed
  offline on the panel design, no family effect at a detectable 0.003 to 0.025 (pilot `RESULTS_LOJO.md`, D-174); plus
  the per-judge bias (uniform leniency), the judge swap on the new evaluator (`JUDGE_SWAP.md`, D-181; now 201 responses on unchanged prompts, Cohen's
  kappa 0.720 three-way, MC shifts of -0.018 to +0.019) and a judge from
  no roster family.
- Threshold sensitivity (yAYU 4, 9W1B 3): `THRESHOLD_APPENDIX.md`; the cross-encoder and alignment-ratio thresholds
  the rebuttal promised to vary no longer exist in the evaluator. Now also the scoring-rule variants
  (`tab:scoring_variants`, matched settings): tolerance halved or doubled, tau 0.964 / 0.927; the absolute-value
  clause off, 310 verdicts (tau 0.782); prescribed digits relaxed, 136 (tau 0.917); the last digit uncapped, 63 (tau
  0.964); proportional part credit, 90 (tau 1.000).
- Inter-judge agreement (9W1B): the May panel 0.725 / 0.776, the replacement panel 0.563 / 0.622 (pilot `RESULTS_E1.md`).
- A fourth branch (ynoK 1; January AC): civil and industrial, 30 templates each, certified.
- Linguistic diversity and template exploitation (meta-review 3; 9W1B 1; cqGs 2; nWW3 2): the paraphrase arm on 275
  expert-kept pairs over 114 templates, nine of eleven models within ±0.05 (gpt-oss-20b and Qwen3-235B-2507 are not); the surface-shortcut audit
  (`SHORTCUT_AUDIT.md`).
- Tool or retrieval baselines (cqGs 1 and 4; ynoK 3; January AC): the open-book condition (405 instances; it lifts
  gpt-oss-20b, +0.044 on the 388 instances answered in both runs, and leaves the two closed models within the margin)
  and the tool condition for the two closed models (Claude Sonnet 5 +0.010, GPT-5.4 mini -0.003) (D-183, D-184).
- Error analysis too coarse, causes of the cliff (nWW3 3 and 4; cqGs 3): attribution validated by error type on the
  experts' labels (pilot `RESULTS_ATTRIBUTION.md`), the human reading of 137 wrong answers of four models at provider defaults, three experts each (411 readings, Fleiss'
  kappa 0.912; `EXPERT_REQUEST.md` B2 now), failure against derivation depth (wrong-answer odds 1.21 per milestone
  pooled, 95% interval 1.09 to 1.33; the slope holds for 3 of 11 models after Holm), the level gap after Holm (1 of 11
  models by Welch's test, gpt-oss-20b; 4 by permutation).
- Validation in the main text (9W1B 6; meta-review 2): the evaluation framework's validation subsection; the final
  evaluator agrees with the experts on 0.930 of verdicts three-way and 0.986 of those not marked partial; milestone
  F1 0.926 by matching alone, 0.958 with the judge.
- Proofreading and the broken reference: moot, the text is new.

## 3. Promises replaced with a reason

- Branch-level reporting for the math-specialised models and the math-pretraining claim (meta-review 2; yAYU 5;
  9W1B 5): the claim is dropped; no math-specialised model is on the roster (the roster rule and the budget set it,
  D-110, D-117). A claim carried without the models would be worse than no claim.
- An interval on rho = 0.632 and DeepSeek V4 Pro's category: those numbers no longer exist; the 100-response study is
  replaced by the 300-trace study with 15 experts.

## 4. Objections answered by construction rather than by argument

- gFWV 1 (perfect kappa unconvincing): the certification ran six rounds (150 of 150 certified), caught 60 of 60
  planted defects, and records hand checks (459 of the 504 comparable ones match within 1%) and timestamps; kappa is
  printed as undefined where every verdict approves. The count of 53 rejected real templates is not in
  `layer2/CERTIFICATION.md`: confirm its source before the letter uses it.
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
- cqGs 3 (causes of the cliff): depth (the depth model), attribution, the human reading by level; the two templates
  behind the earlier sensitivity are repaired (section 7).

## 5. What the letter must disclose before a reader of the repository finds it

The September audit's findings on the published templates; the evaluator's measured defects; the answer-check corrections
after the run (238 verdicts moved, 234 to correct, every one read; one known wrong credit) and the two refinements of
the arithmetic check on an expert's readings before its final held-out reading (0.505, 0.752, then 0.905 on a fresh
sample); and that the new numbers are high because the models of September 2026 solve the evaluation set at the top on
final answers. These are letter material: the paper describes the final check and its validation (decided 2026-10-04).

Added for October (section 7): the two chemical templates repaired after the run, because their wording did not pin
the answer, re-certified in rounds 5 and 6, and their 30 items re-run for every model and condition; and the one
re-score with the final evaluator, which moved 276 verdicts over the eleven models (198 raised by the symbolic step,
77 lowered and 1 raised by the number rule) and the mean score per model by -0.0064 to +0.0073. The symbolic rule was
settled on the experts' grades of these same verdicts, so its 1.000 / 0.990 is in-sample; the earlier readings it did
not see give 1.000 / 0.840.

## 6. The two retrieval conditions

Both ran on 2026-10-03 (D-183, D-184). The paper reports the corrected open-book condition (`openbook2`, 405
instances) and the tool condition for the two closed models. The letter may say that a first open-book build was
discarded after a review found its equation blocks defective, and that the open-weight model could not be served with
a tool through the endpoints the provider rule admits (three attempts recorded under `traces/tool/_*`, local).

## 7. The October revision (6 to 9 October), point by point

An internal mock review of the draft set the work; each change is listed with where the paper
holds it, its number, and the numbered points of the earlier cycles it also answers. "Internal" marks a change that
answers no numbered point; the letter lists those under further changes.

| Change | Where | Number | Points |
|---|---|---|---|
| Two chemical templates repaired after the experts' B4 reading (wording did not pin the answer), re-certified, their 30 items re-run everywhere | section 3.3, `appendix:certification` | round 5: 2 templates, one rejection; round 6: unanimous; 150 of 150 certified; six near-miss templates named in the error analysis | gFWV 1, gFWV 3 |
| Matched settings as the headline: the three models that return no reasoning tokens at their providers' defaults evaluated with reasoning at medium effort on all 2,250 items, and every table, figure and sentence reporting that configuration; their default rows beside | sections 5.1 and 6, Table 1, `tab:matched`, `appendix:models` | reasoning adds +0.118 (GPT-5.4 mini), +0.070 (Gemma 4 26B), +0.038 (Gemini 3.1 Flash-Lite); GPT-5.4 mini enters the first five (0.969); tau 0.709 between the orderings; 36 of 55 pairs differ; Qwen3-235B-2507 has no setting | R1 M4, R2 W6, Q7 (unequal inference conditions) |
| Symbolic answers graded by experts; an equivalence step on five templates, the other four scored by the numbers they state | `appendix:scoring` | 302 verdicts, 362 readings (two-reader agreement 0.950, kappa 0.882); 1.000 / 0.990 in-sample, 1.000 / 0.840 on earlier readings | 9W1B 6, meta-review 2 |
| Final evaluator: 1% cap on the answer's own last digit, two reader fixes, the milestone floor; one re-score; re-validation | section 4, `appendix:validation` | 0.930 three-way, 0.986 non-partial; milestone F1 0.926 (0.958 with the judge); experts' readings: 144 items, 0.847 three-way with the current verdict | 9W1B 6, meta-review 2 |
| Scoring-rule sensitivities | `tab:scoring_variants` | as in section 2 (threshold sensitivity) | yAYU 4, 9W1B 3 |
| Final-answer check audited by answer kind on both expert populations, every disagreement given a cause, the stricter readings bounded on the saved responses; the rule unchanged, no re-score | section 4, `appendix:validation` (`tab:answer_kinds`, `tab:answer_cases`), `appendix:scoring` (What Must Match, Sensitivity of the Rule), Limitations | study: every scalar, array and symbolic response agrees, multipart least (0.800); 44 discrepant readings: the check stricter in 28, more lenient in 12, experts split in 4; stricter readings move no model's FAC by more than 0.008, only GPT-5.4 mini and Muse Glimmer 30B (0.0002 apart) trade places (tau 0.964); a correct rounding up to 0.026, same exchange; matched settings (`answer_audit.py`, `ANSWER_AUDIT.md`) | internal (supervisor's comment) |
| Milestone Coverage under four readings; verbosity; single-path counts | section 6, `appendix:results` | Claude Sonnet 5 above DeepSeek V4.1 Flash: by matching alone (Holm p 0.0005) and intermediate-only (0.032), not as scored (0.069) nor route-adjusted (1.000): route conformity; MC separates 15 of 55 pairs, matching alone 20; verbosity +0.0020 MC per 10 numbers shown (0.0005 to 0.0041); no-instance share 0.000 to 0.069 | internal |
| Judge swap on unchanged prompts | `appendix:results` | 201 responses; kappa 0.720; MC shifts -0.018 to +0.019 | meta-review 1, yAYU 1, 9W1B 2, gFWV 4 |
| Error-analysis sample brought current | `appendix:error_analysis`, `tab:errors` | 137 items on 77 templates, 411 readings, kappa 0.912; shares with the no-error readings removed; the readings are of the default-setting responses, stated in the paper | nWW3 3, nWW3 4, cqGs 3 |
| Depth model; level gap after Holm | section 6, `tab:level_gap` (with the depth slopes) | odds ratio 1.21 per milestone pooled (1.09 to 1.33), 3 of 11 models; gap: Welch 1 of 11 (gpt-oss-20b), permutation 4 | cqGs 3, nWW3 4, yAYU 2, 9W1B 4 |
| Paraphrase bound | `appendix:paraphrase` | 275 pairs over 114 templates; nine of eleven within ±0.05 | meta-review 3, 9W1B 1, cqGs 2, nWW3 2 |
| Taxonomy and coverage: domains from the NCEES FE exam specifications through a three-LLM panel; per-area table; uncovered core areas named | section 3.1, `tab:area`, `appendix:taxonomy` | 15 domains, three per branch; 42 areas | nWW3 1, meta-review 3 |
| Framing: title, "process evaluation", synthetic scope; Limitations, Ethics, Conclusion | front matter, main text (`conclusion.tex`, `limitations.tex`, `ethics.tex`) | eleven scope sentences in Limitations | meta-review 3 |
| Related work: dynamic and functional benchmarks, reference-free evaluators, PRMBench, five science benchmarks; ThermoQA and FinChain deltas | section 2 | the extended comparison table was dropped on 2026-10-09 (external figures, not asked for) | internal |
| Contamination: the generator draws fresh instances from new seeds | Limitations | | internal |
| Difficulty: the domain experts' original labels and their protocol | section 3.1 | 58 / 58 / 34 | internal |
| Worked example; listings for one template per level | `appendix:scoring`, `appendix:template_examples` | | meta-review 3 (presentation) |
| Figures: the same designs, data brought to the final results | figures | | internal |

Claims the letter must not carry over from the earlier cycles:
- The "complexity cliff" is now model-specific: at matched settings the Easy-minus-Advanced gap holds after Holm for
  1 of 11 models by Welch's test (gpt-oss-20b) and 4 by permutation.
- The two tiers of the earlier draft are a result at the providers' defaults; at matched settings the models form a
  graded order (36 of 55 pairs differ) and only gpt-oss-20b differs from every other model.
- The contrast between frontier and open-weight models is removed from the paper; the July meta-review's summary
  mentions it, so the letter says it was removed (the authors give the reason).
- The math-pretraining claim stays dropped (section 3).
