# E4 — E3 plus a check of the arithmetic each trace shows

Run 2026-09-19. **300 gold traces plus 60 Gemma-4-31B, $0.00, deterministic.**
Evaluator `evaluators/e4_arith.py`; checker `evaluators/arith.py`; every figure below
is printed by `analysis/e4_analysis.py` and `analysis/e4_claim_audit.py`.

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
weaker models or harder items is exactly what the robustness cohort and the human
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
of these flags are worth comparing with the human first-error labels in X1.

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
is skipped, not guessed at.
