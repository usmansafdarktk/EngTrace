# E3 — deterministic milestone verification

Run 2026-09-18. **300 of 300 gold traces plus 60 Gemma-4-31B traces, $0.00, no
model in the loop.** Scores in `scores/e3/`; the evaluator is
`evaluators/e3_milestones.py`, the milestone rule `evaluators/milestones.py`.

This is not yet a verdict on E3. That needs the human labels (X1). What it
establishes is that E3 is **measuring something, deterministically, at zero
marginal cost**, and where it differs from E0.

## How E3 decides

A trace reaches a milestone if it states that quantity within 0.5%, under exact
units or one of 13 unit factors. The score is coverage: milestones reached over
milestones required, order-free.

A milestone is **a value the template computed, the gold states, and the question
does not give** — derived by rule from each frozen item's own internals. Those are
recovered by regenerating the item at its recorded seed with the repo's frame-local
capture, and all 60 items reproduce **byte-identically** before a single value is
used. No template was edited; the freeze cannot move.

| Template | Milestones per item |
|---|---|
| lorentz_force | 7 |
| fluid_particle_acceleration | 6-12 |
| phasor_addition | 7-12 |
| manning_rectangular_discharge | 5-6 |
| server_configuration_selection | 4-5 |
| aoq_ati_rectifying | 4 |
| normal_depth_iteration | 2-5 (converged quantities only) |
| rackett_equation_volume, gas_phase_concentration, two_state_steady_state | 3 |
| incompressible_continuity | 2-3 |
| damping_classification | 2 |
| cantilever_double_integration | 1-2 |
| coaxial_capacitance, reynolds_number_flow_regime | 1 |

## Results

| Model | E3 coverage | E0 final-answer acc | E0 reasoning F1 |
|---|---|---|---|
| claude-opus-4.7 | **0.877** | 0.700 | 0.449 |
| gpt-5 | 0.867 | 0.617 | 0.474 |
| gemini-3.1-pro | 0.860 | 0.633 | 0.454 |
| deepseek-r1 | 0.855 | 0.633 | 0.402 |
| gemma-4-31b *(robustness)* | 0.836 | — | — |
| llama-3.1-70b | 0.285 | 0.150 | 0.115 |

E0 figures are the shared-seed run (`04a27d893f49`). Correlation with E0's
final-answer correctness: E3 **0.471**. Correlation between E3 and E0's reasoning
F1: 0.511.

**On the small-model question:** Gemma-4-31B scores 0.836, just below the four
frontier models and far above Llama. E3 functions on a small open-weight model's
traces and gives it a plausible score. Qwen3.8-27B is not scored: its column fails
T1 and T7 (FINDINGS R-F2), and the harness refuses unverified columns.

## The null baseline, and the fix it forced

Every score above would be worthless if a trace could reach milestones by
coincidence. So each trace was also scored against the milestones of a **sibling
item** — same template, different values — where no match is legitimate.

The first version failed that test: traces "reached" **23%** of a sibling's
milestones. Only 3.8% of that was values genuinely shared between siblings; the rest
was coincidence, from a 2% tolerance (borrowed from E0) and unit scaling, on
numbers that crowd into the same ranges.

| Tolerance | Real | Null | Separation | Corr w/ correctness |
|---|---|---|---|---|
| 2.0% (v1) | 0.799 | 0.230 | 0.569 | 0.414 |
| 1.0% | 0.776 | 0.120 | 0.656 | 0.454 |
| **0.5% (v2)** | **0.749** | **0.040** | **0.709** | **0.471** |
| 0.2% | 0.720 | 0.023 | 0.697 | 0.446 |

v2 uses 0.5%, which is the display tolerance the milestone rule already uses to
decide the gold states a value. It sits on a plateau (0.2% gives nearly the same
result), so it is not a peak picked from the table. The tightening hit Llama
hardest, 0.421 to 0.285: coincidental matches had been inflating the weakest model
most.

Unit scaling accounts for 7% of the milestones reached on the gold 300 (67 of 951), and at 0.5% it
adds 0.073 to real coverage for 0.021 of null.

## Where E3 and E0 disagree

**E3 high, E0 low.** Every one of the four largest gaps is a known E0 defect:
- `rackett_equation_volume`: GPT-5 traces reach every milestone (E3 1.00), and E0
  scores reasoning F1 0.00, because of its final-answer parser bug (E0-F1).
- `reynolds_number_flow_regime`: E3 1.00, E0 0.00. It is a classification item with
  no number on the Answer line (E0-F2), so E0's final-answer check compares against
  an incidental value.

**E0 high, E3 low.** These show E3's own weakness:
- `coaxial_capacitance` and `reynolds_number_flow_regime` have **one** milestone
  each, so E3 is all-or-nothing there. One missed value scores 0.
- `phasor_addition#0`, `cantilever_double_integration#1`: E3 0.50 against E0 F1
  above 0.8, on traces with correct final answers. They are candidates for the
  "alternative path" check — a correct route that states different intermediates.

## Known limits

- **Values, not reasoning.** A right number reached for the wrong reason counts as
  reached. E4 (symbolic equation checking) exists to close that gap.
- **Few-milestone items are brittle.** Where the gold states only one intermediate,
  E3 reduces to a final-answer check.
- **The alternative-path question is open.** E3 is order-free by design. Whether it
  is *path*-free is exactly what the human labels' "correct-but-alternative-path"
  category will test.

## Why E3 is cheap to trust

It is deterministic: re-running it gives identical scores. E0 is not, and the next
section of FINDINGS.md (E0-F7) measures by how much.
