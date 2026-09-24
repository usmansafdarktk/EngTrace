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
  experts' verdict marks 3 such traces. Their *step* labels mark 87 — see Finding 5.
- **The single largest available improvement is a correct final-answer check.** E0's
  check agrees with the experts on **76%** of traces, which puts a number on E0-F1 (a
  restated input read as the answer) and E0-F2 (no number on the Answer line). As a
  predictor of soundness it scores 0.812, where an accurate answer check scores 0.974.

  **That gap is now closed** (`evaluators/answer.py`, measured by
  `analysis/answer_check.py`). E0 errs almost entirely in one direction - of its 72
  disagreements with the experts, 68 are traces the experts call **correct** and E0 calls
  wrong. A deterministic check that reads the trace's answer segment, takes each target
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

E0 understates every model by about 21 points and ranks **GPT-5 fourth**; the experts and
the corrected check both put it first. Five defects account for it, each measured:

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
| E3 — deterministic matching | 0.921 | 0.916 | 0.919 |
| E4 — E3 plus stated-arithmetic checks | 0.921 | 0.916 | 0.919 |
| **E5 — E3, then a judge on what E3 misses** | **0.925** | **0.991** | **0.957** |

E5's judge recovers nearly every milestone E3's number-matching misses (recall 0.916 →
0.991) at no cost in precision. That confirms, against experts, what the synthetic
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
- **VersaPRM's failure is confirmed:** it finds 7% of the incorrect steps.

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

- **No *evaluator* detects flawed reasoning behind a correct answer.** Every AUROC sits
  near chance; the best, E2's lowest step reward at 0.583, does not separate from E0. A
  deterministic digit check does - see below.
- **E3, E4 and E5 are significantly worse than chance-level E0 here**, because they
  score a trace by milestones a correct answer already implies. Their strength at the
  milestone level (Finding 2) is not a strength at this question.
- **175 of the 178 incorrect steps are calculation slips and 3 are conceptual**, so what
  this set mostly holds is arithmetic that does not change the answer. That is worth
  stating plainly in the paper rather than presenting the set as deep reasoning failure.

A second labelling round to enlarge this set was costed and rejected (D-096): the
signals that would select candidates barely enrich (the 72B's minimum reward under 0.20
yields 28% against a 25% base rate), so ~170 new traces would need labelling to add ~50
hard cases, and the intervals would not tighten enough to change any conclusion.

### The rule that does work (`analysis/digit_rule.py`)

Since 164 of the 167 flaws are calculation slips, the guide's own rule for them can be
applied by machine: *rounding is not an error, a wrong digit is*. For every arithmetic
claim a trace writes, recompute the left side from the numbers the trace itself shows
and ask whether the displayed right side is a correct rounding at the precision shown.
That is E4's checker with its 1% tolerance replaced by the displayed precision.

| Hard case, 228 correct-answer traces | step precision | step recall | step F1 | trace AUROC |
|---|---|---|---|---|
| E4 as it ships (1% tolerance) | 0.154 | 0.034 | 0.055 | 0.480 |
| tolerance 0.1% | 0.346 | 0.101 | 0.157 | 0.527 |
| **the digit rule** | **0.506** | **0.472** | **0.488** | **0.655** (0.594–0.717) |
| *E2, 72B, for comparison* | 0.246 | 0.255 | 0.250 | 0.583 |

The digit rule doubles the 72B PRM at step level and is the only thing in the pilot
whose hard-case interval clears chance, at no API cost and with an auditable flag: it
names the claim, the value shown and the value recomputed (`ln(0.17911) = -1.71918`,
computes to -1.71976). E4's tolerance, not its design, was the problem.

Two limits. Recall is 0.472, bounded by what the checker can parse, so extending parse
coverage is the next gain. And 82 flagged steps in correct-answer traces are not marked
incorrect by the experts; a sample of those should go to one expert in the adjudication
format already used, since each is either a rounding chain the rule should tolerate or a
slip the experts missed - and both answers are worth having.

## What this means

1. **Fix the final-answer check first.** It is the dominant trace-level signal, and E0's
   parser is wrong on a quarter of traces.
2. **For reasoning progress, milestone coverage with a residual judge (E5) is the
   strongest candidate:** F1 0.957 against experts, and cheap ($0.47 for 300 traces) —
   but Finding 5 bounds the claim: it tracks progress, it does not catch a flawed step
   behind a right answer.
3. **Step-level error detection is not solved by an off-the-shelf process reward
   model.** The 72B is the best available. Recalibrating its threshold against these
   labels, on a split so the threshold is not fitted and reported on the same data, is
   the obvious next experiment.
4. **Report the hard case as a negative result.** It is honest, it is measured on 87
   traces, and it is the clearest open problem the pilot identifies.

## Follow-ups this turned up

- The judge probe labelled Claude's `aoq_ati_rectifying#3` step as CLEAN. The experts
  unanimously found a real slip in it: 0.93^49 = 0.0285538, not 0.02857, so P(X=1) is
  0.0999, not 0.1000. It should be relabelled SLIP in `judge_probe.RELABEL`, and E2's
  probe validation re-run.
- The guide does not say how to treat a notational slip whose computed result is right
  (GPT-5's `normal_depth_iteration#1`: a dropped minus sign in a displayed expression,
  correct value). The experts called it correct; the guide should say so explicitly.
- Electrical engineering has the lowest between-rater step κ (0.584) and the lowest
  intra-rater agreement (0.891). If any branch needs a second look, it is that one.
