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
- **VersaPRM's failure is confirmed:** it finds 7% of the incorrect steps. Raising its
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

Since 164 of the 167 flaws are calculation slips, the guide's own rule for them can be
applied by machine: *rounding is not an error, a wrong digit is*. For every arithmetic
claim a trace writes, recompute the left side from the numbers the trace itself shows
and ask whether the displayed right side is a correct rounding at the precision shown.
That is E4's checker with its 1% tolerance replaced by the displayed precision.

| Hard case, 228 correct-answer traces | step precision | step recall | step F1 | trace AUROC |
|---|---|---|---|---|
| E4 at the 1% tolerance it shipped with | 0.154 | 0.034 | 0.055 | 0.480 |
| tolerance 0.1% | 0.346 | 0.101 | 0.157 | 0.527 |
| the digit rule, as D-097 measured it | 0.506 | 0.472 | 0.488 | 0.655 (0.594–0.717) |
| **the digit rule as E4 now ships it** | **0.750** | 0.320 | 0.449 | 0.639 (0.587–0.692) |
| *E2, 72B, for comparison* | 0.246 | 0.255 | 0.250 | 0.583 |

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
| E2, 72B, lowest step reward | 0.878 | 0.833–0.916 | 0.789–0.956 | 4.2 | 72 |
| E5 | 0.886 | 0.835–0.934 | 0.784–0.978 | 4.2 | 71 |
| E3 / E4 | 0.834 | 0.773–0.888 | 0.726–0.941 | 3.9 | 78 |
| *baseline: the experts' answer verdict* | 0.974 | 0.945–0.996 | **0.945–0.993** | **0.9** | 326 |

**The slice is worth about 80 independent traces, not 300**, for separating evaluators.
Clustering at the item rather than the template — the weaker assumption — gives a design
effect of about 1.5–1.8 instead, so the true penalty sits between the two.

The minimum difference this design could detect, at 80% power and 95% confidence, follows
from the clustered standard error:

| Comparison | difference from E0 | template-level 95% | smallest detectable |
|---|---|---|---|
| E5 | +0.036 | −0.093 to +0.173 | 0.192 |
| E2, 72B, lowest step reward | +0.029 | −0.056 to +0.109 | 0.118 |
| E0-3J | −0.040 | −0.084 to +0.003 | 0.062 |
| E3 / E4 | −0.016 | −0.144 to +0.115 | 0.190 |
| *the experts' answer verdict* | **+0.124** | **+0.043 to +0.218** | 0.127 |

Three things follow, and they should be stated in the paper rather than left implied:

- **The headline finding survives clustering.** The experts' answer verdict beats every
  evaluator by +0.124 with a template-level interval that excludes zero. That the final
  answer dominates the trace-level verdict is not an artifact of treating instances as
  independent.
- **"No evaluator beats E0" is a statement about power, not about equality.** Every
  evaluator difference is smaller than what this design can detect (0.06–0.19 AUROC). The
  pilot rules out large differences between evaluators; it never could have resolved small
  ones. E1's −0.002 is the exception that proves the point: E1 tracks E0 trace by trace, so
  its paired standard error is tiny and a difference of 0.002 is "significant" and
  meaningless.
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

**Nothing detects a conceptual defect. 0 of 60, on every free evaluator.** A trace that
computes correctly, states the wrong rule for what it is doing, and lands the right answer
passes every deterministic evaluator the pilot has. That is the cleanest statement of the
open problem, and it is now independent of the annotation guide.

**The digit rule is perfect where it can parse and blind where it cannot.** Split by where
the defect sits: **45 of 45** inside a claim `arith.py` parses, **0 of 15** on a stated
value with no parseable working behind it. So the recall of 0.472 in Finding 5 is a parse
ceiling, not a rule failure — extending parse coverage is the remaining gain, and changing
the rule is not.

**E4's 1% tolerance is confirmed as the defect, by construction:** 1.000 on planted errors
of 1% or more, **0.000** on the 14 below it, where the digit rule scores 1.000.

**The gold-validation corrections cost 5 detections and remove 20 of 30 false flags.** The
five are last-digit slips at 1e-7 to 1e-5 relative; the false-alarm rate on expert-clean
traces falls from 26.7% to 11.7%. That trade now has a number on it.

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
   strongest candidate:** F1 0.957 against experts, and cheap ($0.47 for 300 traces) —
   but Finding 5 bounds the claim: it tracks progress, it does not catch a flawed step
   behind a right answer.
3. **Step-level error detection is not solved by an off-the-shelf process reward model,
   and does not need one.** The 72B is the best PRM available and flags one real error in
   four inside correct-answer traces; its threshold is not the problem (D-100: calibrating
   it on held-out halves gains −0.006). A deterministic digit check reaches three in four
   on the same steps and costs nothing, and now ships inside E4 (D-101).
4. **Report the hard case as a negative result, and scope it.** Arithmetic flaws behind a
   correct answer are deterministically detectable — three flags in four are real, at no
   cost. Conceptual flaws behind a correct answer are detected by nothing the pilot has:
   0 of 60 planted defects, and only 3 such steps exist in the labelled corpus to study
   (Finding 7). That is the clearest open problem the pilot identifies, and it is now
   stated independently of the annotation guide.
5. **Report the power, not just the intervals.** Clustered by template, this slice can
   detect AUROC differences of about 0.12–0.19 between evaluators and no smaller (Finding
   6). Every null result here should be read as "this design rules out a large difference",
   which is a claim the data supports, rather than "the evaluators are equivalent", which
   it does not.

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
