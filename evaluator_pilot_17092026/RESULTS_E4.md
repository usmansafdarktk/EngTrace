# E4 — E3 plus a check of the arithmetic each trace shows

Run 2026-09-19, **re-run 2026-09-24 under the digit rule (D-097)**.
**300 gold traces plus 60 Gemma-4-31B, $0.00, deterministic.**
Evaluator `evaluators/e4_arith.py`; checker `evaluators/arith.py`; every figure below
is printed by `analysis/e4_analysis.py`, `analysis/e4_claim_audit.py`,
`analysis/arith_gold_validation.py` and `analysis/digit_rule.py`.

> **The 2026-09-19 result below was a tolerance artifact.** E4 shipped at a 1%
> relative tolerance, which asks whether a number was fabricated and cannot see a
> slip in the last digit — which is what 175 of the 178 flawed steps behind a correct
> answer are (`hard_case_pool.py`, current labels). The re-run scores the arithmetic on the experts' own rule instead, and
> the last section records what changed. The headline is that the *milestone* score
> barely moves and the *arithmetic* score stops being decorative: on the hard case
> its AUROC goes from 0.494 (0.459, 0.532) — chance — to **0.661 (0.599, 0.726)**,
> the best number any evaluator in this pilot reaches there.

## What E4 checks, and why not symbolic equivalence

E3 asks whether a trace **states** each milestone. E4 adds: does the arithmetic the
trace **shows** actually produce the number it states? For every milestone E3
found, E4 looks for a computation in the trace that ends in that value and
classifies the milestone:

| Status | Meaning |
|---|---|
| verified | a shown computation ends in this value, and it checks out |
| **contradicted** | computations end in this value and none checks out — **right number, shown work does not produce it** |
| stated | the value appears with no checkable computation behind it |

The Suggested Actions phrase E4 as algebraic equivalence to the gold formula. That
needs a trace's symbols aligned with the gold's (`V_sat`, `V_{sat}`, `V_L`), which
free text does not support reliably (D-085). Checking the shown arithmetic is what
catches the failure E4 exists for: a right number reached by the wrong route.

## Results

| Model | E3 | **E4** | verified | contradicted | arithmetic consistency |
|---|---|---|---|---|---|
| claude-opus-4.7 | 0.877 | **0.877** | 0.520 | 0 | 1.000 |
| gpt-5 | 0.867 | **0.863** | 0.507 | 1 | 0.990 |
| gemini-3.1-pro | 0.860 | **0.860** | 0.445 | 0 | 0.984 |
| deepseek-r1 | 0.855 | **0.855** | 0.360 | 0 | 0.985 |
| gemma-4-31b *(robustness)* | 0.836 | **0.836** | 0.281 | 0 | 0.992 |
| llama-3.1-70b | 0.285 | **0.281** | 0.056 | 1 | 0.812 |

Correlation with E0's final-answer correctness: E3 0.471, **E4 0.470**, verified
coverage alone 0.190.

## The finding: on this corpus, E4 adds almost nothing to E3

Of **1,494** milestones checked across 360 traces, **2** are contradicted — 0.13%.
Both were read and both are genuine:

- **gpt-5, `normal_depth_iteration#1`**: writes
  `1.853 − (0.2892 × (−0.160))/(−2.8172) = 1.86921`. As written it evaluates to
  1.8366; 1.86921 is what a **plus** gives. The sign in the shown formula is wrong,
  the reported value is right (milestone 1.87).
- **llama-3.1-70b, `normal_depth_iteration#1`**: writes
  `(3.8 + 2×1.500×2) × 1.500 = 14.100 m²`. That evaluates to 14.7; the stated 14.100
  matches the milestone (14.0998). The shown coefficient is wrong, the value is not.

So the failure mode E4 targets — right number, wrong work — **exists but is rare**:
two cases in 1,141 reached milestones. That is a result, not a disappointment. It
says E3 alone is not being fooled by copied or guessed values here, and that E4's
extra machinery would change a score in 2 of 360 traces. Whether that stays true on
weaker models or harder items is exactly what the robustness cohort and the expert
labels can test.

## The side metric, and how far to trust it

**Arithmetic consistency** is every numeric claim in a trace, milestone or not. It
is informative — Llama's 0.812 against 0.98-1.00 for everything else is a real
difference — but it is **not** a validated per-model measure, and should be quoted
with this caveat.

A reproducible sample of 30 of the 71 inconsistent claims (`e4_claim_audit.py 30`,
seed 7), classified by hand:

| | n |
|---|---|
| **Real slips in the shown work** | **19** |
| Checker false positives | 10 |
| Unclear | 1 |

**Precision ≈ 2/3.** The real slips are not only Llama's: Gemma writes
`0.4525 × 0.089 × 80 / 100 = 0.32218` (10× too large), DeepSeek twice drops the
factor 2 from a wetted perimeter, and GPT-5 has the sign slip above. These are
exactly what a step-level annotator would mark "Calculation Error", so the positions
of these flags are worth comparing with the expert first-error labels in X1.

The false positives that remain, and why they were not "fixed":

| Cause | n | Why not fixed |
|---|---|---|
| cancellation: `7.71731 − 7.65047`, a 0.007% rounding amplified to ~1% by subtracting near-equal numbers | 4 | an absolute tolerance near zero would also pass genuine small errors in probabilities |
| prose word read as a unit (`0.091 for x = 0`) | 2 | widening the prose list is whack-a-mole |
| a sum cut at a line break (`a + b + c = a + b`) | 2 | needs cross-line parsing |
| `ppm` conversion, `8 E I` read as 8 | 2 | one-offs |

Each fix available here trades false positives for false negatives. The score does
not depend on them — it rests only on the two milestone contradictions, both read —
so the side metric is reported with its measured precision instead.

## How the checker was validated, in order

This is the part worth reusing, because the first validation was not enough.

1. **Plants.** 38 cases in `arith.py`'s self-test (`python evaluators/arith.py`),
   each expecting true, false or no-claim.
2. **Gold.** Gold arithmetic is correct by construction, so every inconsistency on a
   gold solution is a checker bug. First cut: **61.5%** consistent (111 bugs).
   After ten fixes: **100%** (227 of 227 claims).
3. **Every milestone flag on real traces, read.** At 100% on gold the checker still
   flagged **15** milestones as contradicted, and **13 were checker bugs** — gold is
   formatted uniformly by templates, traces are not. Five more fixes left **2**, both
   genuine, and both are now plants so a regression is caught.
4. **A sample of the non-milestone flags, read.** That found `\mathrm{m}^2` losing
   its unit, Greek `μ`, `%` and `deg/rad` conversions, clause fragments and a lost
   binary minus. Claude's consistency went 0.962 → **1.000**: every "inconsistency"
   in its traces had been the checker.

**Validating on gold is necessary and not sufficient.** Read what a checker flags on
real output before reporting any number it produces.

## Parse coverage

Across all traces: 32% of equation segments evaluated, 44% symbolic (skipped, as
intended — they contain variables), 24% unparseable. The unparseable share rose as
the rules tightened, which is the safe direction: a segment the checker cannot read
is skipped, not guessed at. (Unchanged by the re-run: the rule below changes how a
parsed claim is judged, not what parses.)

---

# The re-run, 2026-09-24: the tolerance was the defect (D-097)

`arith.ARITH_TOL` is 1%. Everything above is measured at that tolerance, and it
answers one question well — *is this number fabricated?* The experts were answering
another: the annotation guide's rule for a calculation slip is **rounding is not an
error, a wrong digit is**, and a wrong digit moves a number by far less than 1%. At
1% the checker reaches **recall 0.034** on the steps inside correct-answer traces.

So `arith.py` now judges every claim twice and carries both verdicts:

| field | rule | question |
|---|---|---|
| `ok` | within 1% relative | is this number fabricated? |
| `ok_digit` | a correct rounding at the precision displayed, of the numbers the trace itself shows | is this digit wrong? |

E4 scores on `ok_digit`; `arith_consistency_tol1`, `milestones_contradicted_tol1`
and each claim's `ok` keep the 1% reading reportable. **E3 is untouched**, by
construction and by check: `e3_coverage` equals E3's own `milestone_coverage` on all
360 rows, and every `reached` flag in E4's meta is E3's.

## The gold validation, re-run — and the two bugs it found

Gold arithmetic is correct by construction, so anything flagged on a gold solution is
a checker bug. `analysis/arith_gold_validation.py`, over 227 claims in the 60 gold
solutions:

| rule | flagged on gold | consistent |
|---|---|---|
| 1% relative | 0 | 1.000 |
| digit rule, bare (as D-097 measured it) | **36** (15.9%) | 0.841 |
| digit rule + units | 12 (5.3%) | 0.947 |
| digit rule + units + shown precision — **what E4 ships** | **0** | **1.000** |

The 36 fall into two classes, both parsing rather than gold arithmetic, and both were
fixed rather than excluded — no gold claim is skipped:

1. **24 unit conversions** (`2.47 mm = 2.47e-03 m`, `121.0×10⁶ mm⁴ = 1.2100e-04 m⁴`,
   `200 GPa = 200000000 kN/m²`). The bare rule compares the displayed precision
   before the unit factor, so a conversion looks like a wrong digit. The digit rule
   now reuses the checker's existing unit handling and compares **after** the factor.
2. **12 cross-product claims in `lorentz_force`** —
   `(-21)·(-5.06e-02) − (-27)·(-6.42e-02) = -0.6699`. The gold states B to three
   figures and computes with the unrounded value, so the operands it *shows*
   reproduce -0.6708. Holding a result to a precision its own displayed operands
   cannot pin down is a checker bug, so `arith.shown_uncertainty` moves each rounded
   literal by half its last digit and widens the tolerance by how far the result
   moves (first order, worst case; correct under cancellation, which a relative
   tolerance is not).
   A bare integer is exact and does not move, or `2/3 = 0.66` would pass.

Both fixes only widen the rule, so they cost recall and buy precision — measured, in
the table below, not assumed. The bare rule is still available and still reported, as
`digit` in `analysis/digit_rule.py`.

## What changed in E4

| model | E3 | E4 before | **E4 now** | verified | arith cons. before | **now** | contradicted |
|---|---|---|---|---|---|---|---|
| claude-opus-4.7 | 0.877 | 0.877 | **0.861** | 0.504 | 1.000 | **0.949** | 0 → 4 |
| gpt-5 | 0.867 | 0.863 | **0.861** | 0.506 | 0.990 | **0.971** | 1 → 2 |
| gemini-3.1-pro | 0.860 | 0.860 | **0.858** | 0.442 | 0.984 | **0.966** | 0 → 1 |
| deepseek-r1 | 0.855 | 0.855 | **0.855** | 0.360 | 0.985 | **0.951** | 0 → 0 |
| gemma-4-31b *(robustness)* | 0.836 | 0.836 | **0.834** | 0.279 | 0.992 | **0.957** | 0 → 1 |
| llama-3.1-70b | 0.285 | 0.281 | **0.279** | 0.055 | 0.812 | **0.719** | 1 → 2 |

Of 2,466 claims the 1% rule flags 71; the digit rule flags 186. Of 1,494 milestones,
contradicted goes **2 → 10** and verified 584 → 576. All ten were read
(`e4_analysis.py` prints them): the two originals, three more clear slips
(claude's `(0.77647)^(2/3) = 0.8407` for 0.84479, llama's squares in
`sqrt(517.448² + 395.568²)`, claude's `2×577.72×19.135 = 22110.5` for 22109.3), and
five marginal ones a human might not mark — a truncation where the rule wants a round
(`= 0.9614` for 0.96149, `= 0.84478` for 0.844794), or a value out by a few units in
the last place shown. That ratio is the measured precision of the rule, not a
surprise: on the expert-labelled steps it is 0.750.

## Against the expert labels

Milestone level (`annotation/score_against_labels.py`), 300 labelled traces:

| | TP | FP | FN | prec | rec | F1 |
|---|---|---|---|---|---|---|
| E3 | 881 | 70 | 82 | 0.926 | 0.915 | 0.921 |
| E4 (a milestone the trace states) | 881 | 70 | 82 | 0.926 | 0.915 | 0.921 |
| E4-ok (contradicted is not credited) | 875 | 67 | 88 | 0.929 | 0.909 | 0.919 |

**Unchanged, and it has to be**: the experts' milestone label asks whether the trace
*obtained* the quantity, which is E3's question. Refusing credit for the ten
contradictions trades 6 true reaches for 3 false ones — noise.

Trace level (`analysis/x1_analysis.py`, 2,000-resample intervals):

| target | score | 1% | digit rule |
|---|---|---|---|
| sound vs unsound, 300 traces | `e4_coverage` | 0.834 (0.773, 0.888) | 0.835 (0.773, 0.890) |
| sound vs unsound, 300 traces | `arith_consistency` | 0.657 (0.586, 0.732) | 0.686 (0.601, 0.772) |
| **hard case**: has an incorrect step, 228 correct-answer traces | `e4_coverage` | 0.397 (0.339, 0.457) | 0.411 (0.352, 0.474) |
| **hard case** | `arith_consistency` | 0.494 (0.459, 0.532) | **0.661 (0.599, 0.726)** |

Both `arith_consistency` rows are scored only on the traces where the checker found
something to check: it is None on 84 of 360 rows (23%), so the hard-case rows cover
**178 of the 228** correct-answer traces. That is the abstention the parse coverage
above predicts, the tables now print the `n`, and every difference from E0 is taken
on the traces both scores cover.

The hard-case row is the result. It is +0.149 over E0 with the difference excluding
zero, against E2's 0.583 and E5's 0.426, and it costs nothing to compute. As a filter
(`analysis/hard_case_pool.py`, base rate 40.8%):

| filter | n | precision | recall | lift |
|---|---|---|---|---|
| arithmetic 1%: ≥1 failed claim | 12 | 0.417 | 0.054 | 1.02× |
| **arithmetic digit: ≥1 failed claim** | **43** | **0.767** | **0.355** | **1.88×** |
| arithmetic digit: ≥2 failed claims | 21 | 0.857 | 0.194 | 2.10× |
| PRM 72B: min reward < 0.05 | 33 | 0.515 | 0.183 | 1.26× |

Per step, inside correct-answer traces (`analysis/digit_rule.py`, 178 steps the
experts mark incorrect):

| rule | prec | rec | F1 | trace AUROC (flag count) |
|---|---|---|---|---|
| 1% (what E4 shipped) | 0.154 | 0.034 | 0.055 | 0.480 (0.437, 0.528) |
| digit, bare | 0.506 | 0.472 | 0.488 | 0.655 (0.594, 0.717) |
| **digit as E4 ships it** | **0.750** | **0.320** | **0.449** | **0.639 (0.587, 0.692)** |

(That last column ranks traces by *how many* claims the rule flags, over all 228
correct-answer traces; the 0.661 above ranks them by E4's `arith_consistency` rate,
which abstains where nothing parsed. Two statistics, both reported, neither mixed.)

The two gold fixes cost 0.15 of recall and buy 0.24 of precision; the trace-level
intervals overlap, so the choice between them is not settled by this corpus. The
evaluator takes the one that passes its own gold validation.

## Does this change the finding?

**Partly, and the part that changes is the part that mattered.**

What stands: E4's *milestone* score still adds almost nothing to E3. Contradicted
milestones go from 2 to 10 of 1,494 — 0.7% — `e4_coverage` moves by at most 0.016 on
any model, its AUROC on sound-vs-unsound moves 0.834 → 0.835, and at milestone level
against the experts E4 is identical to E3. Right-number-wrong-work remains rare, and
nothing here rescues milestone coverage as a hard-case signal: at 0.411 it is still
below chance, because a trace with a slip usually still reaches every milestone.

What does not stand: "E4 adds nothing" was a statement about **E4's arithmetic
check**, and at 1% that check was answering the wrong question. The same machinery,
judged by the experts' rule, is the only signal in this pilot that separates a flawed
correct-answer trace from a clean one — AUROC 0.661, precision 0.767 as a filter,
$0.00, no model in the loop. The original section "on this corpus, E4 adds almost
nothing to E3" should be read as: *E4's milestone arithmetic adds almost nothing; the
per-claim arithmetic was never the milestone score, and it is worth more than either*.

## What the hard-case result can and cannot carry

`arith_consistency` at 0.661 against E0's 0.542 is the only evaluator score in the pilot
that separates on the hard case, and the difference (+0.149) excludes zero when traces are
resampled. It does not when **templates** are: the interval is (+0.000, +0.300) against a
minimum detectable difference of 0.211 (RESULTS_X1 Finding 6). The 300 traces are 15
templates, and this comparison runs on the 178 of 228 correct-answer traces where the
checker finds something to check, so its effective sample is smaller still.

The claim this supports is "a deterministic digit check is the most promising signal the
pilot found for flawed reasoning behind a correct answer, at three real flags in four, and
the design cannot certify the margin". It does not support "E4 solves the hard case".
