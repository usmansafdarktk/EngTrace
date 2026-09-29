# X1 — the evaluators against the expert labels

Scored 2026-09-23 on **version 2 of the expert labels**. 300 traces, each labelled by
three experts from the trace's own branch (15 experts, 1,020 submissions for the ground
truth plus the shared calibration set). The labels live outside git, in
`experts_filled_labels/version_2/labels/`; version 1 is kept beside it.

Version 2 adds three things version 1 did not have:

- **A reason on every step called incorrect.** 1,042 of 1,042 (100%).
- **A verification round.** Each expert re-labelled 7 of their own traces (about 10%),
  blind to their first pass, so intra-rater consistency is measured rather than assumed.
- **An adjudication round.** Every step where a branch's three experts split 2–1 went
  back to all three, blind and shuffled; the majority of that second pass is the
  consensus label and replaces the majority for that step alone.

Two sets were re-annotated after the verification round showed them least consistent,
and both originals are kept under `superseded/`: chemical's `che-3.1`, and electrical's
`ele-3.1` (2026-09-24). Electrical's between-rater step κ rose from 0.584 to **0.762** and
the pooled figure from 0.750 to 0.781. The re-annotation moved the two patterns the
diagnostics identified: "not a claim" fell from 6.1% of that expert's steps to 0.4%, and
steps marked incorrect rose from 10.3% to 22.2%, in line with the branch's other two
experts.

**Replacing the set moved nothing; re-adjudicating it moved 12 steps.** Rebuilding the
truth with `ele-3.1` in place of `ele-3`, holding the adjudication fixed, leaves all 2,091
step labels and all 300 verdicts identical - every step where the replaced set could have
swung a majority had already been settled by the branch's three experts reviewing it
blind. What the re-annotation did change is which steps are *disputed*: electrical's 2-1
splits fell from 77 to 47, 12 of them new. Those 12 went back through the same blind
review, which called all 12 incorrect, so the truth now holds **388 incorrect steps rather
than 376**. Verdicts and final answers are unchanged, so every trace-level result in this
document stands; the step-level numbers are re-computed on the 388.

That two-stage result is worth reporting as a robustness property: an entire expert's set -
a fifteenth of the annotation effort, and the least self-consistent one - was replaced, and
the ground truth moved by 12 of 2,091 labels, none of them a verdict.

Scripts: `annotation/score_against_labels.py` (headline comparison),
`annotation/verification_report.py` (the two new rounds), `analysis/x1_analysis.py`
(intervals, baselines, step level), `analysis/hard_case_pool.py` (Finding 5).

## The ground truth

Between raters, Fleiss κ over each trace's own three experts:

| | steps | milestones | final answer |
|---|---|---|---|
| pooled over branches | **0.781** | 0.880 | 0.966 |
| chemical / civil / electrical / industrial / mechanical | 0.772 / 0.645 / 0.762 / 0.867 / 0.810 | 0.931 / 0.989 / 0.899 / 0.874 / 0.659 | 0.974 / 0.968 / 0.964 / 0.980 / 0.904 |
| calibration set, all 15 experts across branches | 0.775 | — | 0.869 |

Within raters, from the verification round (792 re-labelled steps):

| | chemical | civil | electrical | industrial | mechanical | pooled |
|---|---|---|---|---|---|---|
| step agreement with their own first pass | 0.972 | 0.935 | 0.945 | 0.915 | 0.991 | **0.947** |
| Cohen's κ | 0.913 | 0.519 | 0.845 | 0.820 | 0.947 | **0.828** |

An expert agrees with themselves on 95% of steps (κ 0.83), and three experts agree with
each other at κ 0.78. The between-rater figure is therefore close to the ceiling their own
consistency sets, which is the point of measuring both. All fifteen sets are covered:
`ele-3.1` was verified on the same seven traces as the branch's other two experts, and the
round that had been run on the set it replaced is kept under `verification/ele/superseded/`.

The adjudication settled **272 split steps across 134 traces** (255 distinct steps; the
shared calibration traces are reviewed by more than one branch). The three reviewers came
back unanimous on 240 of them and 2–1 on 32. It **changes 85 ground-truth step labels** and
confirms the majority on 170, and every consensus row carries a written reason. After
adjudication no step, milestone or verdict is left disputed, and the count of steps the
experts call incorrect rises from 322 to **388** of 2,091 — adjudication mostly resolved
splits *towards* an error being real. Electrical was adjudicated twice: the round run
against `ele-3` is kept under `adjudication/superseded/`, and the current round covers all
47 of the splits that `ele-3.1` leaves.

The experts judge **229 traces sound and 71 unsound**.

## Finding 1 — at the trace level, the final answer decides almost everything

**The experts' own final-answer verdict predicts their "reasoning sound" verdict with
AUROC 0.974**, better than every evaluator. Of the 228 traces with a correct final
answer, the experts' verdict calls only **3** unsound.

| Separating sound from unsound traces | AUROC | 95% CI | minus E0 |
|---|---|---|---|
| E0 — the published framework | 0.850 | 0.804–0.892 | — |
| E0-3J | 0.810 | 0.759–0.857 | −0.040 |
| E1 | 0.849 | 0.803–0.890 | −0.002 (significant, negligible) |
| E2, 72B, fraction of steps ok | 0.862 | 0.810–0.907 | +0.012 |
| E2, 72B, lowest step reward | 0.878 | 0.833–0.916 | +0.029 |
| E3 | 0.834 | 0.773–0.888 | −0.016 |
| E4 | 0.834 | 0.773–0.888 | −0.016 |
| E5 | **0.886** | 0.835–0.934 | +0.036 |
| *baseline: E0's own answer check* | 0.812 | 0.768–0.852 | **−0.038** (significant) |
| *baseline: the experts' answer verdict* | **0.974** | 0.945–0.996 | **+0.124** (significant) |

Three consequences:

- **No evaluator beats E0 at the trace level.** E5 ranks highest but its interval
  overlaps E0's. E1's −0.002 is statistically distinguishable and practically nothing.
  Without chemical engineering, E0-3J, E3 and E4 are significantly *worse* (−0.05,
  −0.09, −0.09) and the ordering is otherwise unchanged.
- **The holistic verdict cannot rank reasoning evaluators on this slice.** "Right
  answer, flawed reasoning" is the case a reasoning evaluator exists for, and the
  experts' verdict marks 3 such traces. Their *step* labels mark 93 — see Finding 5.
- **The single largest available improvement is a correct final-answer check.** E0's
  check agrees with the experts on **76%** of traces, which puts a number on E0-F1 (a
  restated input read as the answer) and E0-F2 (no number on the Answer line). As a
  predictor of soundness it scores 0.812, where an accurate answer check scores 0.974.

  **Most of that gap is closed** (`evaluators/answer.py`, measured by
  `analysis/answer_check.py`). E0 errs almost entirely in one direction - of its 72
  disagreements with the experts, 68 are traces the experts call **correct** and E0 calls
  wrong, and the corrected check calls 56 of those 68 correct. A deterministic check that reads the trace's answer segment, takes each target
  from a quantity the gold **computed**, scores every part of a multi-part answer and
  returns correct / partial / incorrect agrees with the experts on **0.947** of the 281
  non-partial traces where E0 manages 0.747, and on 0.893 of all 300 three ways. See
  Finding 1b. E0 itself has not been re-scored with it (D-098): the number it fixes is a
  property of the traces, not of E0, and re-running E0 alone would put it on a different
  check from E0-3J and E1.

### Finding 1b — what a corrected answer check is worth

| Answer check, against the experts' verdict on 300 traces | agrees |
|---|---|
| E0, the published framework | 0.747 |
| **corrected check, non-partial traces** | **0.947** |
| corrected check, three-way (correct / partial / incorrect) | 0.893 |

| Answer accuracy per model | experts | E0 reports | corrected |
|---|---|---|---|
| gpt-5 | 0.950 | 0.617 | **0.950** |
| claude-opus-4.7 | 0.917 | 0.700 | 0.950 |
| deepseek-r1 | 0.917 | 0.633 | 0.817 |
| gemini-3.1-pro | 0.867 | 0.633 | 0.867 |
| llama-3.1-70b | 0.150 | 0.150 | 0.183 |
| **all 300** | **0.760** | **0.547** | **0.753** |

*Addendum, 2026-09-29 (D-137).* The corrected check read LaTeX's `3.04 \times 10^{5}` as three
numbers, because its exponent rule ran before `\times` was rewritten; that is fixed, with thousands
written `28\,570` or `11{,}003`. Re-run on these 300 traces (`analysis/answer_check.py`), 10 verdicts
change, all 10 to the experts' verdict. Agreement is now **0.982** on the non-partial traces and
**0.927** three ways. Per model the check gives gpt-5 0.950, claude-opus-4.7 0.950, deepseek-r1 0.950
(0.817 above), gemini-3.1-pro 0.900 (0.867), llama-3.1-70b 0.183; all 300, 0.787. The split-half
tolerance fit still gives 0.0015 and 0.0020. The figures above are kept as published;
`full_run_28092026/validate_scorer.py` reproduces them with the code they were published with.

E0 understates accuracy by 21 points overall - 33 for GPT-5, 28 for DeepSeek R1, 23 for
Gemini, 22 for Claude, none for Llama - and ranks **GPT-5 fourth**; the experts and the
corrected check both put it first. These are sizes on this slice, which was built to
over-represent the answer types the comparator reads worst (12 of its 60 items are scalar,
against 89 of the benchmark's 150 templates), so they do not carry over to the full
benchmark: the shift there has to be measured when its table is regenerated. Five defects
account for it, each measured:

1. **The gold value was the last number in the solution.** For `manning_rectangular_discharge`
   that is the **3 in `m^3/s`**, so a trace scored correct if it wrote its unit in ASCII and
   wrong if it wrote `m3/s` in unicode. This one defect explains every manning and
   fluid-acceleration disagreement. Targets are now quantities the gold *computed*.
2. **Restated inputs** were read as the answer (E0-F1).
3. **Multi-part answers** write one marker per part, so reading after the *last* marker kept
   only part (b) and discarded the value.
4. **Notation**: `5.52 x 10^-5` in unicode superscripts, LaTeX `\times`, `\frac{5}{2}`,
   subscripts and thousands separators all read as wrong numbers or none.
5. **One relative tolerance cannot work**: `114` for 113.55 is a correct rounding at the
   precision shown (0.40% out) while `50.20` for 50.05 is wrong (0.30% out). A value now
   passes within a relative tolerance, or within one unit of its own last digit, or of the
   gold's.

The relative tolerance is the one fitted number, so it is chosen on one half of the traces
and reported on the other: 0.876 and 0.905. A 10-case self-test in `evaluators/answer.py`
pins each defect. 32 disagreements remain, 15 of them `aoq_ati_rectifying`, whose question
asks six quantities while its gold states one on the answer line, and 7 `lorentz_force`.

## Finding 2 — milestone coverage is accurate, and E5's judge earns its place

| Milestones, against the experts' "obtained" | precision | recall | F1 |
|---|---|---|---|
| E3 — deterministic matching | 0.926 | 0.915 | 0.921 |
| E4 — E3 plus stated-arithmetic checks | 0.926 | 0.915 | 0.921 |
| E4 — not crediting a milestone whose arithmetic fails | 0.929 | 0.909 | 0.919 |
| **E5 — E3, then a judge on what E3 misses** | **0.930** | **0.989** | **0.958** |

E5's judge recovers nearly every milestone E3's number-matching misses (recall 0.915 →
0.989) at no cost in precision. Refusing to credit a milestone whose shown arithmetic does
not check out buys nothing here either (0.919): milestones are mostly right in traces that
state them, whatever the arithmetic around them does. That confirms, against experts, what the synthetic
validation in RESULTS_E5 suggested: its REACHED credits are sound. E4 adds nothing over
E3 at this level, as RESULTS_E4 found.

## Finding 3 — the 72B process reward model ranks steps well, and over-flags where it matters most

| Steps, against the experts' "incorrect" | AUROC | 95% CI | precision | recall | F1 |
|---|---|---|---|---|---|
| Qwen2.5-Math-PRM-72B | **0.825** | 0.798–0.850 | 0.539 | 0.515 | 0.527 |
| Qwen2.5-Math-PRM-7B | 0.767 | 0.735–0.799 | 0.441 | 0.416 | 0.428 |
| VersaPRM | 0.629 | 0.593–0.666 | 0.385 | 0.064 | 0.110 |
| *inside correct-answer traces only (178 incorrect steps):* | | | | | |
| Qwen2.5-Math-PRM-72B | 0.752 | 0.718–0.785 | **0.246** | 0.255 | 0.250 |

- **The 72B ranks steps well overall** (AUROC 0.83) and now flags slightly fewer steps than
  the experts mark (358 flags against 388 incorrect steps), with just over half its flags
  real. The richer labels moved its precision from 0.47 to 0.54 without moving its AUROC
  much: the model was partly being scored against missing labels, not failing.
- **Inside correct-answer traces, where a step-level evaluator would add most, 75% of
  its flags are false alarms** and it finds 26% of the real errors.
- **VersaPRM's failure is confirmed:** it finds 6% of the incorrect steps. Raising its
  threshold to ~0.9 restores its recall and still leaves it at F1 0.345 against the 72B's
  0.520, so the failure is discrimination, not scale.
- **The 0.5 cut-off is not tuned.** Fitting it on half the traces and reporting on the
  other, both ways round, moves the 72B's held-out F1 by −0.018 and +0.006, intervals
  straddling zero: the stock threshold already sits on the flat top of the curve
  (`analysis/prm_threshold.py`). Inside correct-answer traces no PRM reaches precision 0.80
  at any threshold, so the over-flagging in Finding 5 cannot be tuned away.

## Finding 4 — the experts contradict E2's reordering of the frontier

| | experts: traces sound | E0 | E1 | E5 | E2 (72B) |
|---|---|---|---|---|---|
| gpt-5 | **0.983** | 0.474 | 0.470 | 0.962 | 0.821 |
| claude-opus-4.7 | 0.917 | 0.449 | 0.448 | 0.946 | 0.864 |
| gemini-3.1-pro | 0.867 | 0.454 | 0.452 | 0.929 | 0.876 |
| deepseek-r1 | 0.933 | 0.402 | 0.401 | 0.923 | 0.921 |
| llama-3.1-70b | 0.117 | 0.115 | 0.117 | 0.321 | 0.560 |

The experts put **GPT-5 first**, and so do E0, E1 and E5. E2 put it **fourth** —
RESULTS_E2 flagged that as possibly a penalty on terse steps, and the experts say it
was. Every evaluator separates Llama from the frontier. Among the four frontier models
the experts' figures sit within about 12 points of each other, too close for this slice
to rank the evaluators on them.

## Finding 5 — the hard case, scored on the step labels

The experts' step labels mark at least one incorrect step in **93 of the 228
correct-answer traces** (178 steps), while their holistic verdict calls 90 of those same
traces sound. Scoring the trace-level question as "does this trace contain an incorrect
step, given the answer is correct" gives the hard-case comparison 93 positives instead
of 3 (`analysis/hard_case_pool.py`).

| Hard case: 228 correct-answer traces, 93 flawed | AUROC | 95% CI | minus E0 |
|---|---|---|---|
| E0 | 0.542 | 0.466–0.615 | — |
| E0-3J | 0.561 | 0.482–0.635 | +0.019 |
| E1 | 0.546 | 0.470–0.619 | +0.004 |
| E2, 72B, fraction of steps ok | 0.549 | 0.475–0.620 | +0.010 |
| E2, 72B, lowest step reward | **0.583** | 0.505–0.656 | +0.044 |
| E3 | 0.392 | 0.335–0.453 | **−0.150** (significant) |
| E4 | 0.397 | 0.339–0.457 | **−0.145** (significant) |
| E5 | 0.426 | 0.382–0.472 | **−0.116** (significant) |
| *baseline: E0's own answer check* | 0.521 | 0.462–0.580 | −0.021 |

- **No *judge-based or reward-model* evaluator detects flawed reasoning behind a correct
  answer.** Every AUROC sits near chance; the best, E2's lowest step reward at 0.583, does
  not separate from E0.
- **E4's arithmetic score does, once its tolerance is fixed** (D-101): 0.661 against E0's
  0.542, +0.149 with a trace-level interval excluding zero. It is the only evaluator score
  in the pilot that separates on this target. Under template clustering the interval is
  (+0.000, +0.300) — it touches zero and sits below the 0.211 this design can detect
  (Finding 6), so the direction is clear and the certification is not. At the 1% tolerance
  E4 shipped with, the same score is 0.494: chance.
- **E3, E4 and E5 score well below chance-level E0 here** (0.39–0.43 against 0.54),
  because they score a trace by milestones a correct answer already implies. Their strength
  at the milestone level (Finding 2) is not a strength at this question. The deficit is
  significant when traces are resampled, but **not** when templates are — see Finding 6.
  The direction is consistent and the magnitude is large; the design cannot certify it.
- **175 of the 178 incorrect steps are calculation slips and 3 are conceptual**, so what
  this set mostly holds is arithmetic that does not change the answer. That is worth
  stating plainly in the paper rather than presenting the set as deep reasoning failure.

A second labelling round to enlarge this set was costed and rejected (D-096): the
signals that would select candidates barely enrich (the 72B's minimum reward under 0.20
yields 28% against a 25% base rate), so ~170 new traces would need labelling to add ~50
hard cases, and the intervals would not tighten enough to change any conclusion.

### The rule that does work (`analysis/digit_rule.py`)

Since 175 of the 178 flaws are calculation slips, the guide's own rule for them can be
applied by machine: *rounding is not an error, a wrong digit is*. For every arithmetic
claim a trace writes, recompute the left side from the numbers the trace itself shows
and ask whether the displayed right side is a correct rounding at the precision shown.
That is E4's checker with its 1% tolerance replaced by the displayed precision.

| Hard case, 228 correct-answer traces | step precision | step recall | step F1 | trace AUROC |
|---|---|---|---|---|
| E4 at the 1% tolerance it shipped with | 0.154 | 0.034 | 0.055 | 0.479 |
| tolerance 0.1% | 0.346 | 0.101 | 0.157 | 0.526 |
| the digit rule, as D-097 measured it | 0.506 | 0.472 | 0.488 | 0.655 (0.594–0.717) |
| **the digit rule as E4 now ships it** | **0.750** | 0.320 | 0.449 | 0.639 (0.587–0.692) |
| *E2, 72B, for comparison* | 0.246 | 0.255 | 0.250 | 0.583 |

*Re-run 2026-09-28 after the full pool's gold validation (D-120) fixed how the checker splits clauses: the 1% and 0.1% readings' trace AUROC moved by 0.001; the digit rule as shipped did not move.*

The rule E4 ships (D-101) is the measured one with two corrections the gold validation
forced: a unit conversion is compared after its unit factor (`0.09024 hours = 5.41 minutes`
is not an arithmetic error), and a result is not held to a precision its own displayed
operands cannot pin down. Both loosen it, trading recall for precision — **three flags in
four are real errors, against one in four for the best PRM** — and both were required to
get the false-flag rate on gold solutions, whose arithmetic is correct by construction,
from 15.9% to zero.

The digit rule doubles the 72B PRM at step level and is the only thing in the pilot
whose hard-case interval clears chance, at no API cost and with an auditable flag: it
names the claim, the value shown and the value recomputed (`ln(0.17911) = -1.71918`,
computes to -1.71976). E4's tolerance, not its design, was the problem.

Two limits. Recall is 0.320 as E4 ships the rule (0.472 before the gold-validation
corrections), bounded by what the checker can parse, so extending parse coverage is the
next gain. And 19 flagged steps in correct-answer traces are not marked incorrect by the
experts; a sample of those should go to one expert in the adjudication
format already used, since each is either a rounding chain the rule should tolerate or a
slip the experts missed - and both answers are worth having.

## Finding 6 — what this slice can detect, once clustering is accounted for

Every interval above comes from resampling traces. The 300 traces are not 300 independent
observations: they are **15 templates × 4 instances × 5 models**, and two instances of the
same template are the same problem with different numbers, sharing a gold solution, a
milestone structure and whatever phrasing an evaluator reacts to. Re-running every
comparison while resampling **templates** gives the honest precision
(`analysis/cluster_bootstrap.py`).

| Trace level | AUROC | resampling traces | resampling templates | design effect | effective n |
|---|---|---|---|---|---|
| E0 | 0.850 | 0.804–0.892 | 0.760–0.927 | 3.7 | 82 |
| E0-3J | 0.810 | 0.759–0.857 | 0.731–0.885 | 2.6 | 115 |
| E1 | 0.849 | 0.803–0.890 | 0.758–0.927 | 3.7 | 81 |
| E2, 72B, fraction of steps ok | 0.862 | 0.810–0.907 | 0.785–0.943 | 2.8 | 107 |
| E2, 72B, lowest step reward | 0.878 | 0.833–0.916 | 0.789–0.956 | 4.2 | 72 |
| E3 | 0.834 | 0.773–0.888 | 0.726–0.941 | 3.9 | 78 |
| E4 | 0.835 | 0.773–0.890 | 0.725–0.946 | 3.9 | 76 |
| E5 | 0.886 | 0.835–0.934 | 0.784–0.978 | 4.2 | 71 |
| *baseline: the experts' answer verdict* | 0.974 | 0.945–0.996 | **0.945–0.993** | **0.9** | 326 |

**The slice is worth about 70 to 115 independent traces, not 300**, for separating
evaluators, depending on the evaluator (82 for E0). Clustering at the item rather than the
template — the weaker assumption — gives a design effect of about 1.5–1.8 instead, so the
true penalty sits between the two.

The minimum difference this design could detect, at 80% power and 95% confidence, follows
from the clustered standard error:

| Comparison | difference from E0 | template-level 95% | smallest detectable |
|---|---|---|---|
| E5 | +0.036 | −0.093 to +0.173 | 0.192 |
| E2, 72B, lowest step reward | +0.029 | −0.056 to +0.109 | 0.118 |
| E2, 72B, fraction of steps ok | +0.012 | −0.075 to +0.106 | 0.133 |
| E1 | −0.002 | −0.005 to −0.000 | 0.004 |
| E0-3J | −0.040 | −0.084 to +0.003 | 0.062 |
| E3 | −0.016 | −0.143 to +0.115 | 0.189 |
| E4 | −0.016 | −0.144 to +0.119 | 0.193 |
| *the experts' answer verdict* | **+0.124** | **+0.043 to +0.218** | 0.127 |

The last row is the answer verdict's margin over E0. Its margin over each evaluator
(`cluster_bootstrap.py`, section 1b), which is what "the final answer beats every
evaluator" has to rest on:

| the experts' answer verdict minus | margin | template-level 95% |
|---|---|---|
| E0 | +0.124 | +0.043 to +0.218 |
| E0-3J | +0.164 | +0.085 to +0.247 |
| E1 | +0.126 | +0.044 to +0.218 |
| E2, 72B, fraction of steps ok | +0.112 | +0.026 to +0.190 |
| E2, 72B, lowest step reward | +0.096 | +0.011 to +0.185 |
| E3 | +0.140 | +0.033 to +0.242 |
| E4 | +0.140 | +0.030 to +0.243 |
| **E5** | **+0.088** | **−0.002 to +0.181** |

Three things follow, and they should be stated in the paper rather than left implied:

- **The headline survives clustering against every evaluator but the best.** The experts'
  answer verdict beats E0 by +0.124 and every other evaluator with a template-level interval
  that excludes zero, except E5, whose interval includes it. That the final answer
  dominates the trace-level verdict is not an artifact of treating instances as
  independent; that it beats the best evaluator is not established at this power.
- **"No evaluator beats E0" is a statement about power, not about equality.** Every
  evaluator difference is smaller than what this design can detect (0.12–0.19 AUROC for
  the evaluators built differently from E0; 0.06 for E0-3J and 0.004 for E1, which track
  it trace by trace). The pilot rules out large differences between evaluators; it never
  could have resolved small ones. E1's −0.002 is the exception that proves the point: E1
  tracks E0 trace by trace, so its paired standard error is tiny and a difference of 0.002
  is "significant" and meaningless.
- **The one positive result on the hard case does not clear the bar either.** E4's
  arithmetic score is +0.149 over E0 there, and clustered by template that interval runs
  (+0.000, +0.300) against a detectable difference of 0.211. It is the largest evaluator
  effect the pilot found and the design still cannot certify it.
- **On the hard case, the deterministic evaluators' deficit is no longer significant.**
  E3, E4 and E5 score 0.39–0.43 against E0's 0.54, but with 15 clusters the difference
  (−0.150 for E3) carries an interval of −0.320 to +0.024. The point estimate is large and
  in the direction Finding 5 describes; the design cannot certify it. The hard case is the
  thinnest part of the slice: effective n falls to 34–93 there.

This is the limit that more analysis cannot fix. It is a property of 15 templates, and the
only remedy is more templates — which the full benchmark run provides for model-level
claims, and which a future evaluator study would need for evaluator-level ones.

## Finding 7 — planted defects: what is caught when the answer is known in advance

Findings 5 and 6 rest on the experts' labels, and two objections follow them. The digit
rule implements the same rule the annotation guide gave the experts, so agreement between
them is partly built in. And conceptual error behind a correct answer cannot be measured at
all on this corpus: the experts found **3** such steps against 175 calculation slips.

`analysis/planted.py` answers both by building a set whose ground truth is true **by
construction**. It starts from the 129 traces the experts labelled clean under both their
own verdict and the deterministic answer check, and plants exactly one defect in each:

| | n | what is changed |
|---|---|---|
| **arithmetic** | 60 | one digit of one displayed intermediate, one character; later steps and the final answer keep their original values |
| **conceptual** | 60 | the reasoning a step states — a criterion flipped, a rule misstated, a stated formula contradicting the one used — with **no digit anywhere in the trace changed** |
| **control** | 60 | nothing |

Every plant is verified: the answer is still correct, an arithmetic plant is proved wrong by
recomputing the claim at full precision or against the gold's own milestone value — never by
the checker under test — and a conceptual plant leaves the trace's digit string
byte-identical. The build is seeded and reproduces byte-identically.

| evaluator | arithmetic | conceptual | false alarm on controls |
|---|---|---|---|
| digit rule, as Finding 5 measures it | **0.750** | **0.000** | 0.267 |
| digit rule, as E4 now ships it | 0.683 | **0.000** | **0.117** |
| E4 arithmetic, the 1% reading | 0.533 | **0.000** | 0.083 |
| E4 milestone "contradicted" | 0.033 | **0.000** | 0.017 |
| E3 milestone coverage falls | 0.017 | **0.000** | — |

**No deterministic evaluator detects a conceptual defect. 0 of 60.** A trace that computes
correctly, states the wrong rule for what it is doing, and lands the right answer passes
every checker the pilot has - which is unsurprising, since none of them reads prose. The
evaluators that could catch it are the judges, and they were asked directly (Finding 8).

**The digit rule, as Finding 5 measures it, is perfect where it can parse and blind where
it cannot.** Split by where the defect sits: **45 of 45** inside a claim `arith.py`
parses, **0 of 15** on a stated value with no parseable working behind it. As E4 ships it,
with the two gold-validation corrections, it catches **40 of 45** and 1 of 15. So the bare
rule's recall of 0.472 in Finding 5 is a parse ceiling, not a rule failure, and the shipped
rule's 0.320 is that ceiling less what the corrections give up - extending parse coverage
is the remaining gain.

**E4's 1% tolerance is confirmed as the defect, by construction:** 1.000 on planted errors
of 1% or more, **0.000** on the 14 below it, where the digit rule scores 1.000.

**The gold-validation corrections cost 5 detections and remove 20 of 30 false flags.** The
five are last-digit slips at about 1e-7 to 3e-5 relative (`router_planted.py` prints them);
the false-alarm rate on expert-clean traces falls from 26.7% to 11.7%. That trade now has a
number on it. The 11.7% is per trace - any flag anywhere in an untouched clean trace. On
the 120 single steps the judges were shown untouched (Finding 8), the shipped rule flags 3.

What this set cannot support: it says what an evaluator catches when a defect of a given
shape is present, not how often models produce such defects or in what mix. The conceptual
defects come from a 43-rule catalogue, varied in kind but not in phrasing, so a detector
could in principle learn the catalogue — it is a diagnostic, not a held-out test. And the
arithmetic row for the digit rule is close to an upper bound rather than an estimate, since
a plantable site must itself be a parseable claim; the unbiased half of that family is the
0.000 on values with no working.

## What this means

1. **Fix the final-answer check first.** It is the dominant trace-level signal, and E0's
   parser is wrong on a quarter of traces.
2. **For reasoning progress, milestone coverage with a residual judge (E5) is the
   strongest candidate:** F1 0.958 against experts, and cheap ($0.47 for 300 traces) —
   but Finding 5 bounds the claim: it tracks progress, it does not catch a flawed step
   behind a right answer.
3. **Step-level error detection is not solved by an off-the-shelf process reward model,
   and does not need one.** The 72B is the best PRM available and flags one real error in
   four inside correct-answer traces; its threshold is not the problem (D-100: calibrating
   it on held-out halves gains −0.006). A deterministic digit check reaches three in four
   on the same steps and costs nothing, and now ships inside E4 (D-101).
4. **Report the hard case with its two halves separated.** Arithmetic flaws behind a
   correct answer can be flagged deterministically — three flags in four are real, at no
   cost, though the shipped rule finds a third of them (recall 0.320) and 40 of 45 planted
   slips it can parse. Conceptual flaws are detected by no checker at all (0 of 60) and by
   the best judges about a third of the time, with no false alarms (Findings 7 and 8).
5. **The evaluator this points to is a router, not a new judge.** The digit rule, as
   Finding 5 measures it, catches every arithmetic defect it can parse and none it cannot;
   E0's two judges catch 87% of exactly those it cannot (MiMo 67%), and the judges are the
   only thing with any conceptual signal. Sending what the checker cannot verify to a judge
   covers much of both, and costs a judge call only on the residue — which is what E5
   already does for milestones, applied to steps. The measured routing makes the case
   concrete: E0 shows a judge the corrupted step half the time for conceptual defects, no
   more often than it shows the clean one, and catches 9 of 60 end to end; a verify-first
   router with MiMo catches 15 of 51 (Finding 8).
6. **Report the power, not just the intervals.** Clustered by template, this slice can
   detect AUROC differences from E0 of about 0.12–0.19 for the evaluators built differently
   from it, and no smaller; its own variants, which track it trace by trace, are resolved
   more finely (0.06 for E0-3J, 0.004 for E1; Finding 6). Every null result here should be read as "this design rules out a large difference",
   which is a claim the data supports, rather than "the evaluators are equivalent", which
   it does not.

## Finding 8 — the judges, asked directly about a planted step

Finding 7's checkers read numbers, not prose, so their failure on conceptual defects cannot
carry a general claim. `analysis/planted_judges.py` asks the evaluators that do read the
step: E0's own two judges, on the framework's own Tribunal prompt, built by its own code,
with exactly one step under review.

**The probe is matched.** Each of the 120 planted defects is judged twice — once in the
planted trace, once in the *same step of the untouched original*. Same question, same gold,
one character different. A judge that flags both is flagging the step, not detecting the
defect, so only the difference counts. 480 calls, $4.67, no failures.

| judge | conceptual | arithmetic | flags the original |
|---|---|---|---|
| GPT-5 | 0.333 (60) | 0.717 (60) | **0.000** |
| Claude Opus 4.5 | 0.133 (60), 0.483 counting *Other* | 0.717 (60) | **0.000** |
| **MiMo-V2.5-Pro** — the judge E5 uses | **0.314** (51) | **0.717** (60) | **0.000** |
| *either of E0's two* | 0.350 | 0.800 | — |

MiMo matters more than the other two. GPT-5 and Claude share families with models on the
evaluated roster, which is the judge/judged objection E1 exists to answer; MiMo does not, and
it performs like them — 0.314 against GPT-5's 0.333 on conceptual defects, the same 0.717 on
arithmetic, and no false alarm on the 111 untouched steps it returned a verdict on. Its
verdicts are also the most decisive: of 51
planted steps it named 11 a Conceptual Error where Opus retreated to *Other* 22 times out of
60. (17 of its 240 calls never returned, after retries, so its conceptual row rests on 51 of
the 60 matched sets.)

**The judges do detect conceptual defects, and the deterministic evaluators never will.**
A third of them, against zero for every checker. That is the result the pilot was missing,
and it changes the recommendation from "this is unsolved" to "only a judge can do this, and
it does it a third of the time".

**Neither judge produced a single false alarm.** 120 untouched steps each, 120 verdicts of
*Alternative Correct*, both judges. Whatever E0's judges are criticised for — RESULTS_E0
found they answer *Alternative Correct* to 93% of what they are sent — over-flagging a
clean step is not it. The problem is sensitivity, not precision.

**Judges and the digit rule are complementary, and the split is exact:**

| arithmetic defect sits… | digit rule, bare | digit rule, as E4 ships it | GPT-5 | Opus 4.5 | MiMo |
|---|---|---|---|---|---|
| inside a claim the checker parses (45) | **1.000** | 0.889 | 0.667 | 0.667 | 0.733 |
| in a stated value with no parseable working (15) | **0.000** | 0.067 | **0.867** | **0.867** | 0.667 |

The bare digit rule is perfect where it can parse and blind where it cannot; the judges are
strongest exactly where it is blind. A checker that routes what it cannot parse to a judge
would cover much of both, and neither component is the one the framework currently spends
on. The split is not exact at the step level: a router settles a whole step when any claim in
it verifies, so a stated value sitting beside a verified claim never reaches a judge - 3 of
the 15 here, of which the digit rule still catches 1 (`router_planted.py`).

### Would E0 show the judge that step? (`analysis/planted_routing.py`)

Finding 8 asks what a judge does when shown a step. This asks whether E0 shows it — Tier 1
run for real on all 120 planted defects and their unmodified originals, with the Tribunal
replaced by a recorder, so the routing is measured and nothing is spent (163 minutes, $0).

| | triggered | **planted step shown** | *same step, original* | end to end, counted per defect |
|---|---|---|---|---|
| conceptual | 0.650 | **0.500** | *0.500* | 9 of 60 = **0.150** |
| arithmetic | 0.833 | 0.767 | *0.683* | 36 of 60 = **0.600** |

End to end counts a defect as caught only when it was shown AND flagged by either of E0's
two judges (`analysis/router_planted.py`). This section first reported the product of the
two rates, 0.500 × 0.350 = 0.175 and 0.767 × 0.800 = 0.613; the product assumes the judges
catch the shown defects as often as the rest, and they do not - 9 of the 30 conceptual
defects shown, against 21 of all 60.

**For a conceptual defect the routing carries no signal whatsoever.** The corrupted step is
sent to a judge exactly as often as the untouched one — 0.500 against 0.500. Tier 1 forwards
it because it cannot match the step to the gold at all, not because anything about it is
wrong. Whether E0 catches a misstated rule is therefore decided by a coin it was already
tossing, and it catches **9 of 60** end to end.

For an arithmetic defect the routing is mildly informative (0.767 against 0.683), since a
corrupted number is a little harder for Tier 1 to match, and it catches 36 of 60 end to end.

One more thing this exposed: E0's answer check calls the final answer **wrong on 32 of the
120 planted traces**, though the plants never touch it and the originals were all correct.
Those traces enter the wrong-answer sample (D3, probability 0.20), so part of E0's routing
today is driven by its own broken answer check — the defect D-098 fixed for accuracy
reporting but deliberately did not re-run E0 against.

### What a router would cost and buy (`analysis/router_residue.py`)

The routing above is E0's. A router would instead send a judge exactly what the checker
cannot verify. That residue is measurable, and it is large:

| | |
|---|---|
| steps with no claim the checker can recompute | **1,324 of 2,091 (63.3%)** |
| traces with at least one | 282 of 300 (**94%**) |
| residue steps per trace | median 4, mean 4.4, max 16 |
| of that residue, steps the experts call incorrect | **14.0%** |
| full run, 11 models: one batched call per trace with residue | 23,265 calls, **$79** at MiMo's rate (untested design) |
| full run, 11 models: one call per residue step | 109,230 calls, **$371** (the design the probe measured) |

So the router is not a small targeted sample: at trace level it forwards almost everything.
Its selectivity is *within* the trace — 63% of steps instead of Tier 1's coin toss, and none
of the 37% already settled deterministically.

What it buys, measured on the planted set (`analysis/router_planted.py`): it forwards 51 of
the 60 conceptual defects, not all of them - the other 9 sit in steps holding a claim that
verifies, which the router settles without a judge. With MiMo as its judge it catches **15 of
51** end to end (0.294; 16 of 60 with GPT-5), against E0's 9 of 60. This section first said
the router shows a judge the corrupted step every time and catches about 0.31; the first was
assumed, and both are replaced by the counts above. On arithmetic the stack - the digit rule
where a step verifies, the router's judge where it does not - catches 49 of 60, against E0's
36 of 60.

**The routing rule has a gap, measured on the experts' labels (D-112).** A step counts as
verifiable when any one of its claims can be recomputed, so a step can be settled by a claim
that passes while its error sits in another. Placing every step the experts call incorrect
(`router_residue.py`, last section):

| incorrect steps | inside correct-answer traces (178) | over all 300 traces (388) |
|---|---|---|
| the digit rule, as E4 ships it, flags the step | 57 | 99 |
| residue: the router sends it to a judge | 63 | 186 |
| settled as verified; the bare digit rule would flag it | 27 | 52 |
| settled as verified; no rule flags it | 31 | 51 |
| **never flagged and never shown to a judge** | **58 (32.6%)** | **103 (26.5%)** |

The 27 are steps the shipped rule declines to flag because its operands' displayed precision
cannot pin the result down; the 31 hold their error somewhere the checker does not reach,
such as a link inside a chained equality (`aoq_ati_rectifying#3` is one). As specified, the
router would never show a judge a third of the slips behind a correct answer. A rule that
forwards every step the checker has not positively cleared would put those steps in front of
a judge, which catches a share of them, at the cost of a larger residue; which rule, and at
what cost, belongs to the router's design (open).

### The router, batched: two smoke checks (`analysis/router_batched.py`, D-113)

The full run can afford the router only batched, one call per trace; one call per step is out of
the budget (D-113). Both checks run rule C - every step the digit rule, as E4 ships it, does not
flag goes to MiMo, a trace's steps in one call - on the framework's own Tribunal prompt, called as
the single-step probe called it. $2.30 for both; every call returned after retries.

**Planted defects, matched as before** (240 calls, $0.92, 6.4 steps per call):

| | sent to the judge | caught, batched | the same defects, one step at a time |
|---|---|---|---|
| conceptual | 58 of 60 | **19 of 58** (0.328); 19 of 60 end to end | 16 of 49 batched, **16 of 49** single |
| arithmetic | 19 of 60 (41 flagged by the digit rule) | **14 of 19** | 14 batched, **12** single, of 19 |
| false alarms on the other steps under review | | **2 of 842** (0.002) | |

Batching costs no detection: on the defects judged both ways the batched judge catches as many
conceptual defects and more arithmetic ones, and it flags almost no clean step. End to end on this
set the router catches 19 of 60 conceptual defects against E0's 9 of 60, and the stack - the digit
rule where it flags, the batched judge on the rest - catches 55 of 60 arithmetic ones against E0's
36.

**The 300 labelled traces** (299 calls, $1.38, 6.6 steps per call; 8 of the 1,971 steps sent were
left unjudged):

| steps the experts call incorrect | precision | recall | F1 |
|---|---|---|---|
| all traces, the digit rule alone | 0.825 | 0.255 | 0.390 |
| all traces, the digit rule plus the batched judge | 0.707 | **0.603** | **0.651** |
| inside correct-answer traces, the digit rule alone | 0.750 | 0.320 | 0.449 |
| inside correct-answer traces, the digit rule plus the batched judge | 0.703 | 0.360 | 0.476 |

The judge more than doubles step-error recall overall, and almost all of that gain is in traces
whose answer is not correct: inside correct-answer traces it adds 7 true flags (64 against 57) for
8 false ones. Ranking the 228 correct-answer traces by flagged steps gives AUROC 0.675 (0.622 to
0.730, traces resampled) against the digit rule's 0.639 and E0's 0.542; template-level intervals
were not computed for it. For comparison, the 72B PRM reaches F1 0.527 over all traces and 0.250
inside correct-answer ones (Finding 3).

**What a batched call costs.** $0.0038 on the planted traces and $0.0046 on the labelled ones,
against the $0.0034 the $79 estimate assumed. At the labelled rate, one call for each of the full
run's traces with a step to send comes to about $114 ($95 at the planted rate).

**What remains untested.** Both checks ran on the pilot's five models; the roster's steps per call,
and so its cost per call, are estimates. E1's three-judge panel was not probed as a panel. And the
planted set is a diagnostic: these rates describe what an evaluator does with defects of this
shape, not how often models produce them.

## Limitations, and what this pilot does not show

Stated here so a reader does not have to derive them, and so the claims above can be read
against them.

**The design is 15 templates.** 300 traces is 15 templates x 4 instances x 5 models, and
instances of a template are the same problem with different numbers. Clustered by template
the slice is worth about **70 to 115 independent traces**, depending on the evaluator, and it
can detect AUROC differences from E0 of roughly **0.12 to 0.19** for the evaluators built
differently from it, and no smaller (Finding 6; E0's own variants track it closely enough to
be measured more finely). Every null result here means "no difference this large", not "no
difference". The finding that matters most clears the bar against all but the best
evaluator: the final answer dominates the trace-level verdict, +0.124 over E0 with a
clustered interval of +0.043 to +0.218, but only +0.088 over E5, with an interval of −0.002
to +0.181 that includes zero.

**The hard-case result is about arithmetic.** 175 of the 178 flawed steps behind a correct
answer are calculation slips. Three are conceptual, which is too few to study. The planted
set (Finding 7) supplies conceptual defects by construction and **no evaluator detects any
of them**, but planted defects say what is caught when a defect of a given shape is present,
not how often models produce one. Nothing here measures how common conceptual error behind
a correct answer is in the wild.

**The digit rule shares its rule with the annotation guide.** The experts were told that
rounding is not an error and a wrong digit is; the rule implements that. Their agreement is
therefore partly by construction, which is why Finding 7 exists - on planted arithmetic
defects, where truth does not come from the guide, it scores 45 of 45 inside a claim the
checker can parse and 0 of 15 outside one.

**The evaluators are not all scored under one configuration.** E0, E0-3J and E1 are scored
with the published framework's own final-answer check, which the expert labels show is wrong
on 72 of 300 traces. The corrected accuracy is reported beside them from
`analysis/answer_check.py` rather than by re-running them, deliberately (D-098): re-running
E0 alone would break the controlled comparison between the three, which requires one check
and not a correct one. A reader comparing the accuracy column with the reasoning column is
comparing two configurations, and the documents say so at each point.

**Model coverage is five models plus a two-model robustness cohort.** Four are frontier
models within 12 points of each other on expert soundness, and Llama 3.1 70B is the only
weak model, at 0.117. The slice can separate frontier from weak; it cannot rank frontier
models, and no evaluator claim here should be read as one.

**The labels are not in the repository.** They are the experts' own annotations and are kept
local by decision, so RESULTS_X1 cannot be reproduced from a clone alone: the scripts are
committed and the numbers stated, but the label files must come from the authors. An
anonymised truth file - codes and majority labels, no per-rater rows - would close that, and
is not yet a decision.

**Two label sets were re-annotated.** Chemical's and electrical's third sets were replaced
after the verification round identified them as least self-consistent, and electrical was
then re-adjudicated. Replacing an entire expert's set changed no ground-truth label;
re-adjudicating the disputes it created changed 12 of electrical's 406 step labels, all to
incorrect, and no verdict (D-099). Both originals are kept under `superseded/`, and the
reliability figures above are the post-replacement ones.

**Two branches invert the reliability ceiling.** Pooled, an expert agrees with their own
first pass (κ 0.83) more than with the other experts (κ 0.78). In civil (0.519 within
against 0.645 between) and industrial (0.820 against 0.867) it is the other way round. Both
are the branches added in September, and each within-rater figure rests on seven
re-labelled traces per expert.

**What the judges were and were not asked.** E0's two judges were asked directly about
every planted step, matched against the same step unmodified (Finding 8, $4.67): they catch
a third of the conceptual defects and four fifths of the arithmetic ones, with no false
alarms. E0's own routing was then measured separately (`analysis/planted_routing.py`, $0): it
shows the judge the corrupted step half the time for conceptual defects — exactly as often
as the untouched step, so the routing carries no signal — and, counted defect by defect,
catches 9 of 60 conceptual defects end to end and 36 of 60 arithmetic ones. E5's judge,
MiMo, was probed the same way (Finding 8); E1's panel was not probed as a panel, which would
cost about $3 on this set.

### What the pilot supports, and what it does not

| Supported | Not supported |
|---|---|
| The final answer nearly determines the trace-level verdict, and beats E0 by a margin that survives clustering | That it beats the best evaluator, E5: +0.088, interval −0.002 to +0.181 |
| E5 is the most accurate milestone evaluator, at $0.47 in judge calls against $6.28–6.53 for an E0 run | That E5 beats E0 at the trace level; E0 has no milestone score to compare |
| An off-the-shelf PRM over-flags inside correct-answer traces, and its threshold is not the cause | That a better-calibrated PRM could not do better |
| A digit-level arithmetic check flags slips at three real errors in four, deterministically and free | That it finds most of them (recall 0.320), or any conceptual error — it finds none |
| No deterministic evaluator detects a conceptual defect behind a correct answer | That the judges cannot: asked directly, E0's two judges catch 35% between them, MiMo 31% |
| E0's judges do not over-flag: 240 clean steps, 240 clean verdicts | That E0's routing is selective: for conceptual defects it shows the corrupted step no more often than the clean one |
| End to end, E0 catches 15% of conceptual and 60% of arithmetic defects of this shape; a verify-first router with MiMo catches 29% of the conceptual ones, asked one step at a time | That these rates transfer to defects models actually make, or to a batched router prompt |
| The expert labels are reliable at kappa 0.78 between raters and 0.83 within, pooled | That every branch meets that ceiling: civil and industrial invert it; or that 300 traces from 15 templates can separate evaluators finely |

## Follow-ups this turned up

- ~~The judge probe labelled Claude's `aoq_ati_rectifying#3` step as CLEAN~~ **Done
  2026-09-27 (D-112).** The experts unanimously found a real slip in it: 0.93^49 =
  0.0285538, not 0.02857, so P(X=1) is 0.0999, not 0.1000. It is relabelled SLIP in
  `judge_probe.RELABEL`, and the probe and E2's validation are re-scored (JUDGE_SELECTION,
  RESULTS_E2). Kimi K3 and Grok 4.6 now lead the probe at 0.96, with MiMo at 0.92; the 72B's
  probe AUROC is 0.923. The same step also shows a gap in the recommended stack: neither
  digit rule flags it, because the slip sits inside a chained equality the checker does not
  compare, and the router would settle the step without a judge, because another claim in
  it verifies (Finding 8, D-112).
- ~~The guide does not say how to treat a notational slip whose computed result is
  right~~ **Withdrawn 2026-09-27 (D-112): the case was misread.** GPT-5's
  `normal_depth_iteration#1` step writes a wrong sign in its displayed formula, but its
  value is wrong too: 1.86921 where the formula gives 1.86942, because 0.046272/2.8172 is
  0.016425, not 0.01621. All three experts mark it a calculation error after adjudication,
  with that arithmetic as the reason; only one expert's first pass called it correct. The
  guide's existing rule covers it ("a wrong digit is a calculation error, even if the final
  answer survives it"), the digit rule flags it too, and the pilot holds no step with a
  notational slip and a right value, so nothing in the guide changes.
- Electrical was the weakest branch when this was first written (between-rater step κ 0.584)
  and its third set was re-annotated (D-099). It is civil now: the lowest between-rater
  step κ (0.645) and the lowest within-rater κ (0.519), which is below its between-rater
  figure; industrial inverts too (0.820 within against 0.867 between). If any branch needs
  a second look, it is civil, and industrial's verification round is worth widening past
  seven traces per expert.
