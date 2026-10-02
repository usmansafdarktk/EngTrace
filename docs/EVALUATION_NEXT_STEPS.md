# The evaluation: what is still to do, and how

Written 2026-10-02 (D-170), after the full run closed: inference, E5, the router, the repeats and the
paraphrase arm are done, the evaluator corrections of D-137 to D-169 are in, and `full_run_28092026/results/RESULTS.md`
is regenerated at `2950875`. This file replaces section 10 of `NEXT_CYCLE_REVIEW.md`, which is no longer
maintained. Each item says what it answers, how it is done, at what unit, what it costs, and where its
output lands. The rules that bind every item: the template is the unit of analysis; every interval
resamples templates; any test added now is labelled exploratory in RESULTS.md, with the date and the
D-entry, because the plan (D-117) was fixed before the data and these were not in it; every number the
paper quotes is printed by a committed script; nothing paid runs without the owner's approval of its
estimate.

The deadline is 12 October 2026 (ARR, anywhere on Earth). Section A is free and fits in two to three
days of work; section B needs the experts, who should be asked once, together; section C needs spend
and a decision each.

## A. Analyses on the data in hand, all free

**A1. The reasoning score compared across models, as the July rebuttal promised.** *Answers:* yAYU 2 and
9W1B 4 (no variance or significance), and the rebuttal's own promise of "Wilcoxon on milestone coverage"
beside McNemar on the answer. *Method:* on the 150 per-template mean E5-strict coverages, for each of the
55 model pairs, (a) the paired template bootstrap interval of the difference and (b) a Wilcoxon
signed-rank test, with the sign-flip permutation printed beside it (the two should agree; the plan's test
for Q1 is the sign flip), Holm over the 55 pairs; the per-template SD columns (within and between) as Q1
has them. Beside it, the stability of the two orderings: Kendall's tau between the models' answer-score
order and their coverage order, with a template bootstrap, and the per-model paired difference between the
two scores is *not* a comparison to make (different scales). *Unit:* the template. *Caveat to print:*
coverage credits stated intermediates; a terser model scores lower without reasoning worse, so the
comparison ranks what traces state, not how soundly they reason. *Where:* a new Q3 block, "Coverage
compared across models (exploratory)"; `analyze.py`. *Effort:* half a day.

**A2. Branch, level and domain scores with intervals and tests.** *Answers:* the meta-review's
"branch-level reporting", and any sentence of the form "branch X is hardest". *Method:* template-level
bootstrap intervals for each model's branch and level means (30 templates per branch, 58 / 58 / 34 per
level); within a model, the ten pairwise branch differences with a sign-flip test and Holm; the smallest
branch difference 30 templates detect at 80% power, from the within-branch template SD. Report the
intervals in the branch table; make no ordering claim that does not survive the test. Domain means stay
descriptive (10 templates each). Also the Q2 gap without `work_isothermal_virial` and
`adiabatic_flame_temperature`, with its interval (D-171; the descriptive values are in `RESIDUAL_INCORRECT.md`). *Where:* "By branch and level" gains interval columns; an appendix table
of branch pairs. *Effort:* half a day.

**A3. Leave-one-judge-out on the published Tribunal, with controls.** *Answers:* the judge-exclusion
ablation the July meta-review names as promised, and yAYU 1 / 9W1B 2 / gFWV 4 on judge-and-judged overlap.
*Method:* on the pilot's 300 traces, from the stored per-judge votes of E0-3J (`evaluator_pilot_17092026/scores/e0_3j/`,
parsed by `e1_analysis.py`): for GPT-5 traces drop the GPT-5 judge, for Gemini traces the Gemini judge,
re-vote with the remaining two under the framework's conservative tie-break, and report the change in the
per-model reasoning score and in agreement with the experts' step labels; a placebo (drop a non-family
judge) and the untouched models (DeepSeek R1, Claude Opus 4.7, Llama) as controls. *Unit:* the trace, with
template-resampled intervals and the pilot's detectable difference stated. *Where:* a new
`analysis/lojo.py` and `RESULTS_LOJO.md` under the pilot. *Effort:* one day. *Caveat:* a historical
analysis of E0, which the full run does not use; it discharges a literal promise and should be one
paragraph in the paper beside the E1 result (swapping the whole panel moved scores by at most 0.004) and
the structural fact that MiMo shares no family with the roster.

**A4. X2, the per-judge bias estimate.** *Answers:* the same objection, as the pilot's `JUDGE_SELECTION.md`
item 4 proposed and `RESULTS_E1.md` left open ("not yet X2"). *Method:* per judge and per generator family,
the mean of the judge's step verdict minus the experts' label on the pilot's 300 traces (E0's two judges
and E1's three), with intervals; a null result settles the objection as well as a positive one. *Where:*
`analysis/x2_bias.py`, reported in RESULTS_E1. *Effort:* half a day, after A3 (same parsed votes).

**A5. The threshold appendix, collated.** *Answers:* yAYU 4 and 9W1B 3 (thresholds never varied). *Method:*
one table from existing outputs: E3's milestone tolerance grid (RESULTS_E3: real, null and separation at
2%, 1%, 0.5%, 0.2%), the digit rule against the 1% and 0.1% readings (RESULTS_X1, the digit-rule table;
the 1% column is already in Q3), E2's 0.5 threshold held out (D-100), the answer check's relative
tolerance at half and double, the half-unit window and the whole-trace reading (Sensitivity). State that
the tolerance's split-half fit (0.876 / 0.905 on the pilot) cannot be redone on the full run, which has no
expert labels. *Where:* an appendix table assembled by a small script that reads the existing JSON
outputs, so no number is typed. *Effort:* half a day.

**A6. Error attribution validated by type.** *Answers:* nWW3 3, cqGs 3 (error analysis coarse; causes).
*Method:* on the pilot's 300 labelled traces, cross-tabulate each component's flags (the digit rule, an E5
MISSING milestone, a router-judge flag by its category) against the experts' step error types
(calculation / conceptual / unsupported; `annotation/guide.md`): precision per type and per component.
Then read the full run's attribution table (Q3, "what points at a wrong answer") through those precisions
as shares of traces flagged, never as counts of error types. *Where:* `analysis/attribution.py` under the
pilot; a sentence in Q3. *Effort:* one day.

**A7. The corpus-wide surface-shortcut audit.** *Answers:* nWW3 2, cqGs 2, 9W1B 1 (template exploitation).
*Method:* D-057's held-out depth-2 rule on all 150 templates' public draws: the blind-guess floor per
classification template and the surface-model lift per template; flag any template above a stated lift
and give the headline answer score without it (the Sensitivity table already drops the four known ones:
tau 1.000 with the headline). *Where:* a `tools/` or `full_run_28092026/shortcut_audit.py` and an appendix
table. *Effort:* one day. *Caveat:* it measures predictability from the question surface, not what the
models did; the paraphrase test is the behavioural half.

**A8. The decoding-settings table (Appendix P), printed from the traces.** *Answers:* the review's finding
that four models wrote no reasoning tokens (D-168), and 9W1B's reproducibility ask. *Method:* per model,
from the trace rows: the request parameters, the token ceiling, the median and the 90th percentile of
completion and reasoning tokens, the share of traces with any reasoning tokens, the served model id and
the providers (`TRACE_REVIEW.md` has most of it). *Where:* extend `trace_review.py` or a small
`decoding_table.py`; one appendix table. *Effort:* two hours. *Then the paper says*: the closed tier's
scores are not read as the models' ceiling (RESULTS_PAPER_NOTES item 7).

**A9. Power statements where the paper makes a null claim.** *Answers:* yAYU 2 / 9W1B 4, and the
"pre-registration" question (9.3 item 1). *Method:* beside every null (the top five, the paraphrase
arm, the coverage comparisons of A1), the smallest difference the design detects at 80% power from the
paired per-template SD, as Q2 and Q5 already do. *Where:* `analyze.py`, one column. *Effort:* two hours.

**A10. Does coverage track verbosity?** *Answers:* the caveat A1 prints, with a number. *Method:*
exploratory: per model, the correlation across templates between mean E5-strict coverage and the mean
number of steps and of claims per trace; and the coverage of a model's fully solved traces against its
wrong ones. If coverage rises with the number of stated steps, say so where coverage is reported. *Where:*
`analyze.py`, descriptive paragraph. *Effort:* two hours.

**A11. Verdict against coverage at the trace level.** *Answers:* "what verification adds" (NEXT_CYCLE_REVIEW
4.2) on this roster. *Method:* per model, the share of fully solved traces with coverage below 0.5 and of
wrong-answer traces with coverage 1.0, with template intervals; the paper's claim that process scores earn
their place on the wrong answers rests on the second. *Where:* Q3. *Effort:* two hours.

## B. Needs the experts, no API spend; one request, batched

**B1. An expert reading of the answer check on this roster.** *Answers:* 9.3 item 4 (evaluator validity on
the deployed roster): the pilot validated the check on other models and 15 templates, and D-168 found
roster-specific misreadings by reading verdicts. *Method:* a stratified sample of about 150 verdicts from
the top models, half incorrect or partial and half correct, in a workbook like the digit rule's
(`flag_sample.py --reader-copy` is the pattern: the trace's answer segment, the gold's answer line, no
model name), read by one or two experts for "correct / partial / incorrect"; report agreement with the
check's verdict and the forms behind any disagreement; fix the D-137 way only if a pattern appears.
*Effort:* a day to build, an hour per expert.

**B2. The human error-analysis sample on the new run.** *Answers:* nWW3 3 and 4, cqGs 3; continuity with the
May version's six-category taxonomy (Appendices R to T). *Method:* 50 to 100 wrong-answer traces for each of
three or four contrastive models (a top model, GLM-5.3 for its empties, a non-reasoning model, gpt-oss-20b),
three annotators, the same decision hierarchy and agreement reporting as before; smaller than May's 2,200
because A6's automated layer covers the population. *Effort:* expert time; the sampling script is an hour.

**B3. A spot-check of the judge on this roster.** *Answers:* the judge's REACHED verdict was validated on the
pilot (it credited none of 88 fabricated values) but not on these models' traces. *Method:* about 100 judged
milestones, half REACHED and half MISSING, shown with the trace to an expert of the branch; report the
judge's precision on each. *Effort:* an hour per expert; the sampling script two hours.

**B4. The two under-specified chemical templates, and the near-miss list** (*added 2026-10-02, D-171*). *Answers:*
what the top models' remaining wrong answers are. `residual_incorrect.py` measured every remaining incorrect verdict's
distance from the gold: 53 of the 97 wrong answers on `work_isothermal_virial` state the flow-work reading of a
question that names neither closed-system nor flow work, and 48 of the 89 on `adiabatic_flame_temperature` lie within
5% of a gold whose heat-capacity data the question does not give. *Method:* one chemical expert reads each question and
gold with three traces, and says which reading the wording supports and whether the estimate can be held to 0.2%
without the data; a glance at the six next templates on `RESIDUAL_INCORRECT.md`'s per-template table. *Outcome:* a
stated limitation and a Q2 sensitivity row (the pool is frozen; no re-score unless the owner decides). *Effort:* an
hour of one expert; the reading sheet half an hour. Goes in the same request as B1 to B3.

## C. Needs spend, and a decision each

**C1. A reasoning-on run for GPT-5.4 mini and Gemini 3.1 Flash-Lite.** *Answers:* the four-models-without-
reasoning finding, which otherwise colours every closed-versus-open sentence. *Method:* the same prompt,
pool and scoring with reasoning effort set explicitly, as a labelled variant (`run_traces.py --variant`),
reported beside the main run, not in place of it; a dry run prices it (the two models' main runs cost $7.74
and $2.48; reasoning multiplies output tokens). *Decision:* run it, or state the setting and leave the
comparison unmade (RESULTS_PAPER_NOTES item 7).

**C2. A judge-swap robustness check.** *Answers:* yAYU 1 / 9W1B 2 on the single judge. *Method:* E5 with a
second judge from another family on a stratified sample of the residual milestones (about $3 to $5), and
Grok 4.6 on the planted set (about $3.40, `JUDGE_SELECTION.md`); report whether any Q3 figure moves by more
than its interval. *Decision:* cheap; worth it if a reviewer's objection to one judge is expected.

**C3. A flagship anchor on the 450-item subsample.** *Answers:* 9.3 item 9 (no frontier model on the
roster). *Method:* two or three flagship models on `subsamples.py`'s 450 items, scored by the same stack,
reported as an anchor with 450-item intervals, outside the pairwise family. *Cost:* about $40 to $55 at the
pilot's per-trace rates. *Decision:* the supervisor's; the paper otherwise needs a sentence on why the
flagship tier is absent.

**C4. A tool-use condition, or the open-book fallback.** *Answers:* cqGs 1 and 4, ynoK 3, the January AC.
*Method:* a Python tool on a stratified subsample for two or three models, same prompt otherwise, same
stack; the digit rule's flag rate under the tool is itself the measurement. Or, cheaper, an open-book
condition that supplies the governing equations in the prompt, bounding what retrieval could add. *Cost:*
$65 to $150 for the tool condition on three models; the open-book condition about the inference cost of
the subsample. *Decision:* the supervisor's; a third "future work" answer should be avoided.

**C5. The Qwen2.5-Math pair.** Only if the math-pretraining claim returns; the plan drops it (D-117).

## Not to do

- No E0 at scale: its answer check is wrong on a quarter of the pilot traces and unlabelled full-run
  numbers would add cost without evidence (D-105).
- No RAG baseline: a corpus the paper cannot release and an attack surface of its own; the open-book
  condition is the honest cheap answer (NEXT_CYCLE_REVIEW 5, 9.3 item 11).
- No Wilcoxon on coverage as a test of a "decoupling" claim: the claim is not made. A1 uses the test to
  compare models on coverage, which is what the rebuttal promised.
- No re-tuning of any threshold on the full run: there are no expert labels there; A5 reports the grids.

## Order

A8 and A9 first (hours; they change sentences in the paper); then A1 and A2 (the reasoning comparisons and
the branch intervals, which decide what the results section may say); then A3 and A4 together (the pilot's
votes); A6, A7, A10, A11 as time allows. Send B1 to B4 to the experts in one request as soon as the
sampling scripts exist, since their time is the long pole; D-171 makes B1 and B4 the most urgent items here, so
their scripts come before A1. C1 to C4 each need a dry-run estimate and the
owner's approval; C1 is the one that changes how the roster is read, C3 and C4 answer standing objections.
