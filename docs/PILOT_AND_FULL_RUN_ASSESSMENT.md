# The evaluator pilot and the full run, assessed against the two review cycles

Written 2026-10-02, for the authors, ten days before the ARR October deadline (12 October 2026, AoE). It reads
the evaluator pilot (`evaluator_pilot_17092026/`, D-078 to D-113) and the full run (`full_run_28092026/`,
D-114 to D-170) against the January and July reviews (`docs/EngTrace Rebuttals Jan 2026.docx`,
`docs/EngTrace_Rebuttal_Jul2026 (1).docx`) and the two meta-reviews (`docs/meta-reviews.md`), and asks three
things: what the two pieces of work establish, what they still leave a reviewer to object to, and what to do in
the time left. It adds one measurement of its own, `full_run_28092026/residual_incorrect.py`
(`RESIDUAL_INCORRECT.md`, D-171), and otherwise cites the committed outputs. Every number below names its source;
none is new unless it names that script.

## 1. In brief

1. **The evaluation is in far better shape than the May version's.** The May paper's reasoning metric rested
   on an answer check that was wrong on 72 of 300 expert-labelled traces, a Tribunal that called two judges
   while describing three, and a 20% sample that reordered models between identical runs (FINDINGS E0-F1 to
   E0-F7). The replacement stack is deterministic where a check exists, validated against 15 experts' labels on
   300 traces, and spends a judge from no roster family only on what no check settles. The full run's analysis
   plan was written before the data (D-117), the template is the unit, every interval resamples templates, every
   family of tests is Holm-corrected, and every figure in `results/RESULTS.md` is printed by a committed script
   from a store whose commit and evaluator hashes are recorded. The July meta-review's first two asks, variance
   and significance for headline claims and a direct answer to judge-and-judged overlap, are answered in
   substance.
2. **Four things a careful reviewer will still find.** (a) The benchmark is saturated at the top: five models
   score 0.967 to 0.976 and are not separable, and this assessment's measurement says the top tier's remaining
   "wrong" answers are mostly not reasoning failures (section 5). (b) The literal promise of the July rebuttal,
   leave-one-judge-out on the Tribunal with placebo and untouched controls, has not been run; the structural
   answer (a judge from no roster family; swapping the whole panel moved scores by at most 0.004) is better but
   is not the experiment that was promised. (c) The stack was validated on five models, none on the roster, and
   on 15 templates; the answer check then had to be corrected five times on forms the pilot templates never
   contained, the last time after the run, and the only roster-specific validation is a domain expert's reading
   of the digit rule's flags. (d) Four of the eleven models ran with no reasoning tokens, so the roster is not
   compared at equal effort and the paper's closed-versus-open sentences have to be written around it.
3. **One finding of this review changes the order of the remaining work.** Of the top five models' 171
   remaining incorrect verdicts, 38 are symbolic items the check scores by their numbers (34 on one template), 57 sit on two Advanced
   chemical templates whose questions do not pin the answer to the check's 0.2% tolerance (53 of the 97
   roster-wide wrong answers on the virial template state, exactly, the flow-work reading of an ambiguous
   question), and 27 are last-digit deviations on items that prescribe their rounding. Only 33 of the 171 are
   more than 5% from the gold. This does not move the headline (the sensitivity table absorbs it), but it
   roughly halves the top tier's Easy-to-Advanced gap, it empties the top models' wrong-answer subsets that Q3's
   process-score analyses rest on, and it makes the expert reading of the answer check on this roster (next
   steps B1) the most urgent item rather than a batch item. Section 5 and D-171.

## 2. What the pilot established, and how far it carries

**What it did** (`PILOT_SUMMARY.md`, `RESULTS_X1.md`). 60 problems, 15 templates over five branches and three
levels, five models, 300 traces; 15 experts, three of the trace's own branch per trace, with a reason on every
incorrect step, a blind re-labelling round and a blind adjudication of every 2-1 split. Seven evaluator
candidates scored against those labels, plus 120 planted defects whose truth is known by construction.

**What it supports, with the number that supports it:**

| claim | evidence | source |
|---|---|---|
| the labels are reliable | Fleiss kappa 0.781 between experts on steps, 0.828 within an expert on a blind re-label; 0.880 on milestones, 0.966 on the final answer | RESULTS_X1, ground truth |
| the published check was the dominant defect | wrong on 72 of 300 traces, 68 of them correct answers called wrong; agreement 0.747 against 0.982 for the corrected check (non-partial), 0.927 three-way | RESULTS_X1 Finding 1b, SCORER_VALIDATION |
| the final answer nearly decides the trace verdict | experts' answer verdict predicts their soundness verdict at AUROC 0.974; 3 of 228 correct-answer traces called unsound | Finding 1 |
| E5 is the most accurate milestone evaluator | precision 0.930, recall 0.989, F1 0.958 against the experts' milestone labels; its judge gave REACHED to 0 of 88 fabricated values | Finding 2, RESULTS_E5 |
| milestone coverage is blind to a flawed step behind a right answer | 93 of 228 correct-answer traces carry an incorrect step; E3/E4/E5 score 0.39 to 0.43 AUROC on finding them, below E0's 0.54 | Finding 5 |
| arithmetic slips can be flagged without a model | the digit rule: precision 0.817, recall 0.427 inside correct-answer traces (as shipped); 45 of 45 planted slips inside a parseable claim | Finding 5, Finding 7 |
| no deterministic check sees a misstated rule; a judge sees about a third | 0 of 60 planted conceptual defects for every checker; MiMo 16 of 51, GPT-5 20 of 60, no false alarm on untouched steps | Findings 7 and 8 |
| swapping the whole judge panel does not move scores | E1 (Grok, MiniMax, MiMo) against E0: pooled 0.378 against 0.379, no model by more than 0.004 | RESULTS_E1 |
| the published Tribunal's inter-judge agreement, never reported before | Fleiss kappa 0.725 (four-way) and 0.776 (binary) with all three original judges connected | RESULTS_E1 |
| the batched router buys step-error recall | precision 0.707, recall 0.603 over all steps against the digit rule's 0.825 / 0.255; inside correct-answer traces 0.703 / 0.360; batching loses no planted detection | RESULTS_X1, the router batched |

**Where it stops, and the paper has to say so:**

- **Power.** Clustered by template the 300 traces are worth 70 to 115 independent ones (design effect 2.6 to
  4.2); the slice detects AUROC differences between evaluators of 0.12 to 0.19 and no smaller. "No evaluator
  beats E0 at the trace level" is a statement about power. The one positive hard-case result, the digit rule's
  +0.149 over E0, has a template-level interval of +0.000 to +0.300 (Finding 6).
- **Population.** Four of the five pilot models are frontier models of September 2026 and none is on the
  roster; the roster is mostly smaller open-weight models (D-110). The stack's agreement with experts on the
  roster's traces is measured for one component only, the digit rule (section 3). The pilot's limits section
  says this; the paper's evaluator section must too.
- **Templates.** 15 of 150. The forms the answer check has since been corrected for (LaTeX numbers, stray
  digits, verdict words, last-digit boundaries, pi-fractions, two-unit lines) occur in none of the 15, so the
  experts' 0.982 does not vouch for them (D-137 to D-139, D-147, D-169), and section 5 adds two more forms the
  slice never saw.
- **The hard case is arithmetic.** 175 of the 178 flawed steps behind a correct answer are calculation slips;
  conceptual error behind a right answer is measured only on planted defects, which say what is caught when a
  defect is present, not how often models produce one.
- **Circularity, mild and stated.** The digit rule implements the annotation guide's own rule, and the answer
  check's tolerance was fitted on the same 300 traces (split-half 0.876 / 0.905); the planted set is the
  independent check for the first, nothing is for the second.
- **Reproducibility.** The labels are the experts' own and are not in the repository (D-099, the owner's
  rule); traces and scores are not committed. RESULTS_X1 cannot be reproduced from a clone. An anonymised truth
  file (codes and majority labels, no per-rater rows) is still a pending decision (PILOT_SUMMARY open item 5),
  and reviewer 9W1B already listed the artefacts as unavailable once.

## 3. What the full run established

**The run** (`README.md`, `TRACE_REVIEW.md`, D-114 to D-133). 150 certified templates × 15 instances from a
private seed, frozen with per-item hashes and a seed commitment; eleven models, 24,750 traces, one per item, at
each provider's default decoding with a 32,768-token ceiling (16,384 for Muse); served model ids recorded on
every row; a twelfth model run and set aside under the roster rule (D-132). 246 unusable traces (1.0%), 243 of
them empty at the ceiling, 153 on six templates, four iterative by construction.

**The scoring** (`SCORER_VALIDATION.md`, `GOLD_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`).
The scorer reproduces all ten pilot agreement figures with the pilot's code and the gold scores clean (2,250 of
2,250 correct, 0 digit-rule flags on 6,875 gold claims). E5 and the router replay the pilot's stored replies and
reproduce it exactly. Judged stages: E5 8,032 calls, every one answered, $36.30; the router 24,503 of 24,506,
$104.16.

**The results** (`results/RESULTS.md` at store commit `2950875`; `RESULTS_PAPER_NOTES.md` says how to report
each):

| | finding | the test behind it |
|---|---|---|
| Q1 | answer score 0.814 (gpt-oss-20b) to 0.976 (DeepSeek V4.1 Flash); the top five within 0.009; 32 of 55 pairs differ; none among the top five | template-level paired bootstrap and sign-flip permutation, Holm over 55; fully-solved test agrees on 54 of 55 |
| Q2 | every model lower on Advanced than Easy, by +0.055 to +0.199; significant for 4 of 11 after Holm (the four lowest-scoring models) | Welch's t with Holm, adopted after the planned permutation proved liberal (D-146); 9 of 11 under the planned test, printed beside it |
| Q3 | E5-strict coverage 0.816 to 0.923, top eight overlap; on wrong answers models still reach 46% to 85% of milestones against a chance floor of 10% to 22%; digit rule flags 0.5% to 22.5% of fully solved traces (precision 0.905 on this roster, 171 of 189 read); router judge flags 1.1% to 22.0% | descriptive with template intervals, as planned; the all-traces coverage table added after the data and labelled (D-170) |
| Q4 | the strongest models solve all 15 instances of 85% to 91% of single-path templates; the weakest 36% to 48% | descriptive |
| Q5 | 277 expert-kept paraphrase pairs over 115 templates: paired change −0.031 to +0.020, none significant after Holm; 10 of 11 within ±5 points at 90%; tau 0.673 between arms against a noise-floor median of 0.782 at the arm's size | sign-flip over templates with Holm; TOST at the plan's detectable difference, margin fixed after the point estimates and said so (D-165) |
| repeats | four cheap models × 300 items × 3: SD 0.004 to 0.013; the same verdict every time on 86% to 94% of items | descriptive; no figure for the other seven |

**What was done right that the May version lacked**, and should be said in the paper's method section: the plan
before the data; the template as the unit; the Holm families; an unusable trace scored 0 by a stated rule with the
rate reported; every correction to the check measured on the gold, the labelled traces and every full-run trace
before adoption, with the changed verdicts read and the between-commit record kept (`PARSER_FIX.md`,
`MATCH_AUDIT.md`, `WORD_AUDIT.md`, `BOUNDARY_AUDIT.md`, `ANSWER_FORM_FIX.md`); the digit rule's precision read by
a domain expert on this roster across three rounds (0.505, 0.752, 0.905) with each fix recorded; the paraphrase
arm with an expert check that rejected 12% of script-passed pairs; and a provenance line on every results file.

**Two housekeeping items.** `results/RESULTS.md`'s provenance line reads "`analyze.py` at commit `2f388c3`,
dirty": the last regeneration ran with uncommitted edits to the script. The committed script is now what produced
it, but the line should be regenerated from a clean tree before any table is copied into the paper. And the
spend: by the figures the records state, the round is at about $535 (pilot about $49; inference $304.62 in the
rows including the set-aside model's $56.71; E5 $36.30; router $104.16; repeats $1.86; paraphrase $39.15), against
the supervisor's roughly $500 per round, so every item in section C of the next steps is new money.

## 4. Against the reviews: what each concern got

The July meta-review's three suggested revisions, the July reviewers' numbered points, and the January points
that bear on evaluation. "Where" names the committed evidence; "state" is how it stands today.

| concern | who | where it is answered | state |
|---|---|---|---|
| variance and significance for headline claims; per-template SD column; template bootstrap on the cliff; paired tests | meta-review 1; yAYU 2; 9W1B 4 | Q1 intervals, SD within / between, 55 Holm-corrected pairs; Q2 with Welch and the detectable gap; McNemar printed as a check | **done**, stronger than promised |
| Wilcoxon signed-rank on the continuous reasoning score, promised beside McNemar | yAYU 2, 9W1B 4 (rebuttal) | not run; next steps A1 (on E5-strict coverage, 55 pairs, template bootstrap beside it) | **done 2026-10-03** (D-173: the coverage comparison in Q3, sign-flip and Wilcoxon over the 55 pairs, Holm) |
| judge-exclusion ablation with placebo and untouched controls, "promised for camera-ready" | meta-review 1; yAYU 1; 9W1B 2; gFWV 4 (Jan) | the structural answer: MiMo-V2.5-Pro shares no family with the roster (D-110); E1 swapped the whole panel, ≤0.004; MiMo matches GPT-5 on planted defects. The literal LOJO is next steps A3 (on E0-3J's stored votes), X2 bias is A4 | **done 2026-10-03** (D-174: replayed offline, reproduces all 179 judged scores; no family effect at a detectable 0.003 to 0.025 F1; X2 finds every judge lenient the same way) |
| inter-judge agreement within the Tribunal | 9W1B (comments) | Fleiss kappa 0.725 / 0.776 for the original panel, 0.563 / 0.622 for E1's (RESULTS_E1) | **done** for the Tribunal; the full run's judge is single, so the equivalent is its validation (REACHED never given to a fake value; planted-defect rates) |
| threshold sensitivity: numeric tolerance, cross-encoder, alignment ratio | yAYU 4; 9W1B 3 | tolerance at half and double, half-unit window, whole-trace reading (Sensitivity); E3's grid 2% / 1% / 0.5% / 0.2% (RESULTS_E3); the digit rule at 1% and 0.1% (RESULTS_X1); the cross-encoder and alignment ratio no longer exist | **done 2026-10-03** (D-175, `THRESHOLD_APPENDIX.md`, collated from the scripts that measured each threshold) |
| the ρ = 0.632 validation on 100 responses, no interval | yAYU 3; gFWV 2 and 4 (Jan) | superseded: 300 traces, 15 experts, three per trace, kappa reported, cluster-robust intervals, power stated | **done**; the paper must say the old study is replaced and why, not add an interval to it |
| branch-level reporting for the math-specialised models; the causal claim | meta-review 2; yAYU 5; 9W1B 5 | the claim is dropped: no math-specialised model on the roster (D-110, D-117) | **dropped with a reason**; say so in one sentence; C5 only if the claim returns |
| branch-level reporting in general | meta-review 2 | by-branch and by-domain tables, descriptive, no intervals | **done 2026-10-03** (D-173): intervals and Welch pairs; one branch pair holds in the whole roster, so no "hardest branch" sentence |
| promote validation and thresholds to the main text | meta-review 2; 9W1B 6 | the material exists | **writing** |
| framing: synthetic scope, physical versus linguistic diversity | meta-review 3; 9W1B 1; cqGs 2 (Jan) | the paraphrase arm (Q5) is a measured answer to the linguistic half | **done beyond what was promised**; state it as a bound (±5 points for 10 of 11), not as "robust" |
| template exploitation / surface shortcuts | cqGs 2 (Jan), 9W1B 1 | the four known shortcut templates are a sensitivity row (tau 1.000); the corpus-wide audit is A7 | **done 2026-10-03** (D-177): two templates newly flagged, the headline unchanged by more than 0.002 |
| a fourth branch | ynoK 1 | civil and industrial added, 30 templates each, certified | **done**, twice over |
| tool-augmented or retrieval baselines | ynoK 3; cqGs 1 and 4 (Jan) | nothing run; C4 | **open**, paid; a third "future work" answer is the risk |
| error analysis too coarse; causes of the cliff | nWW3 3 and 4, cqGs 3 (Jan) | Q3's attribution table and failure-against-depth (descriptive); A6 validates attribution by error type on the pilot's labels; B2 is the human sample | **partly**: A6 done (D-176, each flag's precision by error type); B2, the human sample on the new run, is in the experts' request (D-172) |
| the evaluation framework validated by the same kind of system it uses | gFWV 4 (Jan) | human labels now; deterministic components; the judge from outside the roster | **done** |
| proofreading, broken reference | all | — | writing |

Two promises were made about numbers that no longer exist (an interval on ρ = 0.632; DeepSeek-V4-Pro's
category). The paper should not try to keep them; it should say the evaluation was rebuilt and why (the
metric-continuity table in `RESULTS_PAPER_NOTES.md`).

## 5. The remaining incorrect verdicts, measured (new; `residual_incorrect.py`, D-171)

D-168 read twelve of the top five models' "incorrect" verdicts and found six to be right answers in another
form; D-169 fixed those forms and moved 238 verdicts. Nobody had asked what the verdicts that remain look like.
`residual_incorrect.py` measures, for every answered, usable trace the store labels incorrect, the closest
approach of any number on its Answer segment to any numeric target under the check's own unit scales, and
buckets the distance. It judges nothing; it says how far the stated answer is from the gold.

**Across the roster** (1,558 incorrect verdicts): 120 within 0.2%, 406 at 0.2% to 1%, 336 at 1% to 5%, 228 at 5% to
20%, 365 beyond 20%, 90 on symbolic items, 11 on word-only targets, 2 with no number stated. So 862 of 1,558 (55%)
lie within 5% of the gold, and every one of the 120 within the relative tolerance sits on one of the 17 templates
whose question prescribes the rounding (`exact`), where the experts' rule requires the digits.

**For the top five models** (171 verdicts): 38 symbolic (34 on `ber_estimation_mary`, whose gold is `a · Q(b)` and
whose incorrect verdicts D-168's reading S found mostly state the right numeric value); 27 within 0.2%, all
exact-digit items; 31 at 0.2% to 1%; 40 at 1% to 5%; 12 at 5% to 20%; 21 beyond 20%; 2 word-only. By template:
`work_isothermal_virial` 38, `ber_estimation_mary` 34, `adiabatic_flame_temperature` 19, `qr_policy_one_iteration`
13, then 8 or fewer each.

Three of those templates were read:

- **`work_isothermal_virial`** (Advanced, chemical; 97 incorrect roster-wide, 38 top-five, 31 of them at 1% to 5%).
  The question asks for "the work in J/mol required to isothermally and reversibly compress 1 mole ... based on the
  virial equation of state truncated to two terms". It names neither the system (closed, W = −∫P dV, which the gold
  uses with Z = 1 + B/V) nor the flow form (W = ∫V dP with Z = 1 + BP/RT, which equals RT ln(P2/P1) + B(P2 − P1)).
  Computing the flow-work value from the gold's own stated B and ideal-gas work and matching it with the check's
  own rule: **53 of the 97 incorrect traces state the flow-work value**, 0 the ideal-gas value, 44 neither. Three
  top models give 11,867 J/mol to five figures on one item whose gold says 11,590; both are textbook readings of
  that sentence. Three chemical experts certified the template (round 2, A A A) after an objection about the
  pressure range; the ambiguity of the wording was not raised.
- **`adiabatic_flame_temperature`** (Advanced, chemical; 89 incorrect roster-wide, 19 top-five). The question
  gives no heat-capacity data and asks for an estimate; the gold iterates a polynomial Cp with coefficients it does
  not cite. 48 of the 89 wrong answers lie at 0.2% to 5% from the gold (top five: 15 of 19), which is where a
  different property table lands. Not read trace by trace; an expert should.
- **The exact-digit templates** (17 of 150; 311 incorrect verdicts roster-wide, 120 of them within 0.2%). The
  (Q, R) question prescribes every intermediate rounding ("carry each rounded value into all subsequent steps");
  the top models compute at full precision and land 1 to 4 units off an integer lot size. `chart_pair_selection`,
  `single_sampling_oc_point` and `newsvendor_normal_demand` are the same: a last digit set by a prescribed rounding
  chain the trace did not follow. These are legitimate zeros under the experts' rule and the question's own terms,
  and they are a finding about the models (prescribed-precision compliance), not about the check; the paper should
  name the 17 templates and the rule.

**What follows.**

1. **Q1.** The headline does not move: the sensitivity table already drops the nine symbolic templates (+0.002 to
   +0.011 per model), and the two chemical templates are worth under a point to each of the top five (6 to 22 of
   30 items lost, on a mean over 150 templates). But the honest
   reading of the top tier is that its genuine wrong answers on this pool number a few dozen out of 2,250, not 26
   to 50: the saturation in section 6 is stronger than Q1 shows.
2. **Q2.** Both templates are Advanced, and the roster loses 237 of 330 points on them. Without the two, the
   top five's Easy-to-Advanced gap falls from +0.055 to +0.073 to +0.029 to +0.047, GLM-5.3's from +0.134 to
   +0.087, and the four significant models' from +0.151 to +0.199 to +0.105 to +0.157 (descriptive, no interval;
   `RESIDUAL_INCORRECT.md`). The cliff for the weak models stands; for the top tier about half of it is these two
   templates, and the paper's cliff paragraph should say so or carry the sensitivity.
3. **Q3.** The top models' wrong-answer subsets (26 to 50 traces) are what "coverage on wrong answers", the
   attribution table and failure-against-depth rest on for those models. For the top five most of those traces are
   a complete, defensible derivation to a value the check cannot credit, which is what their coverage on wrong
   answers (0.54 to 0.75, readable) also shows. State the subset sizes and this reading beside those rows, or
   restrict the Q3 wrong-answer analyses to the models with more than 100 wrong answers.
4. **The experts.** B1 (an expert reading of the answer check's verdicts on this roster) stops being a batch item:
   the two chemical templates need one chemical expert to say which reading the question supports, and whether
   the flame question can pin an answer to 0.2% without the data. The pool is frozen and inference is done, so the
   outcome is a stated limitation and a sensitivity row, not a re-score, unless the owner decides otherwise. The
   next-steps file now carries this as B4, in the same request as B1 to B3.
5. **The method.** Near misses concentrating on a template, with several models agreeing against the gold, is a
   check the certification did not make and a free script can: `RESIDUAL_INCORRECT.md`'s per-template table lists
   the next candidates (`vdw_solve_for_volume`, `pfr_volume_changing_rate`, `best_hydraulic_rectangular_section`,
   `manning_rectangular_discharge`, `annulus_flowrate`, `pitzer_correlation_z`). Worth a sentence in the paper's
   validation section and a line in the limitations.

## 6. What a reviewer will raise, and what the paper can say

- **Saturation.** Five models at 0.967 to 0.976, not separable at 150 templates, and section 5 says the residue is
  mostly form. The paper cannot present EngTrace as stress-testing the frontier on final answers; its headroom is
  the Advanced tier (top five 0.92 to 0.93, and about half of that gap is two templates), depth (wrong-answer rate
  0.037 to 0.074 on items with six or more milestones for the top five, against 0.249 to 0.291 for the four weakest),
  prescribed-precision compliance, and the process scores. Say it as a finding: the pool is solved at the top by
  models of September 2026, and the benchmark's discrimination now lies in where and how models fail rather than
  whether. No flagship is on the roster (D-110's rule and the budget); C3's 450-item anchor is the only way to
  put one there, about $40 to $55, the supervisor's call.
- **The roster is not compared at equal effort.** Four models wrote no reasoning tokens and are four of the five
  lowest scores (D-168). Appendix P must say this per model (A8, two hours), and no sentence may read the closed
  tier's scores as a ceiling. C1 (a reasoning-on variant for GPT-5.4 mini and Gemini 3.1 Flash-Lite) is the one
  paid item that changes how the roster reads; it needs a dry run and approval.
- **The check was corrected after the run, twice.** Every correction is measured, read and recorded, none moved
  a pilot verdict against the experts, and the forms were absent from the pilot's templates. That is defensible
  only if the paper says it plainly and gives the counts (D-169: 238 of 24,750, 234 to correct). Section 5 adds
  that two more forms remain unarbitrated and names them.
- **Validation population.** The agreement figures are on five other models and 15 templates; the digit rule's
  0.905 is the one roster-specific figure. B1 and B3 close most of this with an hour of expert time each; the
  sampling scripts are the long pole and should be built first.
- **A single judge.** MiMo-V2.5-Pro decides 11% to 25% of required milestones (10% to 29% of those REACHED) and
  every router-judge flag. Its defence is the pilot's validation and its independence from the roster; C2 (a
  second judge from another family on a stratified sample, $3 to $5) is cheap insurance if the objection is
  expected.
- **Exact-digit items and symbolic items.** 17 templates require the digits; nine symbolic templates are scored
  by the numbers they state. Name both sets, give the sensitivity without them, and say that a last-digit miss on a
  prescribed rounding scores 0 by the experts' rule.
- **The hard case.** The paper's reasoning claim must be the modest one: coverage measures progress, the digit
  rule flags slips at 0.905 precision and under half recall, and only the judged router has any signal on a
  misstated rule, at about a third of planted defects. `RESULTS_PAPER_NOTES.md` has the wording; keep to it.
- **Reproducibility.** The private pool and seed are released at publication; the labels are not. Decide the
  anonymised truth file before submission, or state in the paper why it is withheld.

## 7. The order for the ten days left

*Update 2026-10-03: items 0 and A1 to A11 of the next steps are done (D-173 to D-178), the experts' kits are built and
awaiting the send (D-172), and C1 to C4 remain, each needing a dry run and the owner's approval. The order below is as written on
2026-10-02.*

The next-steps file's order is right in substance; this review moves two things forward.

1. **Today, free:** regenerate `results/RESULTS.md` from a clean tree (the "dirty" provenance line); A8 (the
   decoding table) and A9 (power statements), hours each.
2. **Tomorrow, free:** the sampling scripts for B1, B3 and the new B4 (section 5), so the experts get one request
   as early as possible; their time is the long pole. B4's reading list is two templates and a question each.
3. **Then, free, in this order:** A1 (the promised Wilcoxon on coverage), A2 (branch and level intervals, with
   the two-template sensitivity for Q2), A3 and A4 together (the promised LOJO and the X2 bias on the pilot's
   stored votes), A5 (the threshold appendix), A10 and A11 (coverage against verbosity and against the verdict).
   A6 and A7 as time allows.
4. **Decisions for the owner and the supervisor, each with a dry run:** C1 first (it changes the roster's
   reading), then C2 (cheap), then C3 and C4 (the standing objections). The round is already at about $535.
5. **In the paper, whatever else happens:** the saturation finding stated as such; the four no-reasoning models;
   the two corrected check rounds with their counts; the 17 exact-digit and nine symbolic templates; the two
   chemical templates as a stated limitation with the cliff sensitivity; the pilot's validation population; tiers,
   not an order, at the top.
