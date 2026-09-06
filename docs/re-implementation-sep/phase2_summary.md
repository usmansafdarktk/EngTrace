# Phase 2 — determinism and solver-in-the-loop

**Date:** 2026-09-06 · **Branch:** `redesign/phase2-determinism`
**Gated on:** Phase C2 (sync point S1), merged as `9151d8c`
**Scope:** four templates · **Decisions:** D-036, D-037

---

## 1. What this phase was for

A benchmark whose items cannot be regenerated from their recorded seeds cannot
be audited, replicated, or corrected after publication. One template in the
corpus was in exactly that state, and it is in published results.

`levenspiel_plot_interpretation` sampled its measurement noise from `np.random`,
which has a generator of its own that `random.seed()` does not touch. The
production path seeds only `random`. **The recorded seed never regenerated the
item.** It was the corpus's only T3 failure, and T3 now passes 150/150.

The second half of the phase is subtler and matters more for where this project
is going. Two templates outsourced their reasoning to a library — `fsolve` and
`quad` — and then printed *the name of the library* as a solution step. A trace
that says "we use a numerical solver" and states an answer is not a reasoning
trace. There is nothing in it for a reader to follow, nothing for a
deterministic verifier to check, and its reproducibility is a solver tolerance
rather than arithmetic. **Both are now stated algorithms with every intermediate
printed.**

---

## 2. The four templates

| Template | What was wrong | What it is now |
|---|---|---|
| `levenspiel_plot_interpretation` | unseeded `np.random`; six blanket `except Exception`; steps numbered `1,2,3,1,2,3,4,5,6` | stdlib RNG only, numpy removed from the file; zero fallbacks; steps 1–8 contiguous |
| `pfr_volume_changing_rate` | printed `quad`'s own error estimate — a machine-dependent number in a gold trace; question promised a non-integer order while `n = 2.0` on 9.6% of draws | closed form `[(1−X)^(1−n) − 1]/(n−1)`, every step hand-checkable; wording corrected |
| `adiabatic_flame_temperature` | `fsolve`; blanket `except`; a `while/else` returning a tuple with no steps and no answer; Cp extrapolated ~1300 K past validity | 6-pass mean-heat-capacity iteration, all intermediates printed; preconditions raise; `CP_PARAMS_COMBUSTION` (D-036) |
| `vibration_isolator_design` | root selection by `max()`; Step 4 hardcoded its verdict; two answer-less returns; `{:,.0f}` broke the parser | selection is a stated predicate; the isolation condition is evaluated; both returns are asserts; plain formatting |

---

## 3. The defect the spec did not list

Spec §2.1 names four defects for `vibration_isolator_design`. It does not name
this one: **the trace failed to reproduce its own answer on 272 of 300
instances.** `omega_n = omega / r` was computed at full precision and printed
from rounded operands, so a reader following the printed numbers did not reach
the printed answer.

That is the P2 violation Phase 1 existed to remove, in a template Phase 1 did
not take. Fixed here, because P1 governs every template a phase edits — but it
is worth recording that a template can sit in a phase's scope for a *different*
reason and carry the previous phase's defect unnoticed. **The spec's per-template
defect lists are a starting point, not an inventory.** Phase 1 and Phase C2 each
found the same thing.

I also introduced this defect myself, once, and T1 caught it within a minute: my
first PFR rewrite printed `(1−X)^(1−n)` rounded and then computed the integral
from the unrounded value. The harness earning its keep is worth noting.

---

## 4. The stated algorithm (D2.4)

`fsolve` is replaced by the mean-heat-capacity iteration from
Smith–Van Ness–Abbott — the same source `CP_PARAMS` follows:

```
T_(k+1) = T0 + (−ΔH_rxn) / Σ_i n_i·<Cp>_i(T0, T_k)

<Cp>/R = A + (B/2)(T + T0) + (C/3)(T² + T·T0 + T0²) + D/(T·T0)
```

- **Initial guess:** 2000 K, stated in the trace.
- **Iterations:** fixed at 6, stated in the trace. Six passes bring all eleven
  reactions within **0.06 K** of their fixed point, so the answer is quoted to
  the nearest kelvin and no run-time convergence test is needed — the count is
  part of the stated method, not a search.
- **Every pass is printed**, with its divisor, as one division.
- The divisor is bound through its display, so each printed line closes exactly.

Worked example, carbon monoxide (the trace prints all six passes):

```
Pass 1: T = 2000.00 K → Σ n_i·<Cp>_i = 115.63 J/K
  T = 298.15 + 283000.0 / 115.63 = 2745.61 K
Pass 2: T = 2745.61 K → Σ n_i·<Cp>_i = 120.09 J/K
  T = 298.15 + 283000.0 / 120.09 = 2654.72 K
...
Pass 6: T = 2662.20 K → Σ n_i·<Cp>_i = 119.71 J/K
  T = 298.15 + 283000.0 / 119.71 = 2662.20 K
```

Every number on those lines is one a reader can check with a calculator. That
was the point.

The PFR integral needed no algorithm at all: `∫₀ˣ(1−X)⁻ⁿdX` is elementary for
`n ≠ 1`, and the sampled range is `n ∈ [1.5, 2.5]`. The closed form agrees with
`quad` to **2.7e-15** over 3000 seeds — so the black box was never buying
accuracy, only opacity.

---

## 5. D-032 decided (D2.7)

Phase C2 handed on the Cp-extrapolation trade rather than deciding it. Decided
as **D-036: split the table by consumer.**

`CP_PARAMS` is fitted to 1500 K; this template integrates to 2844 K. Restricting
the range is not available — the flame temperature is an **output**, not a
sampled input. Refitting `CP_PARAMS` would fix this template and damage
`sensible_heat_temp_dependent_cp`, which lives at 298–1000 K where a wide-range
fit is worse. So `CP_PARAMS` is untouched and a second table,
`CP_PARAMS_COMBUSTION`, is fitted to NIST over 298–3000 K for the four product
species and read only here.

**Result:** all eleven flame temperatures now agree with a direct NIST Shomate
solve to within **1.9 K** (was 5–65 K low). Methane computes **2327 K** against
the 2326.35 K complete-combustion reference — **0.03%**, from 0.65%.

**A near miss.** My first evidence script integrated a single Shomate range
across a span covering several and reported methane at 2191 K — which would have
said the extrapolation made the answer *better*. It surfaced only because
Reviewer H had independently derived 2327.2 K in Phase C2 and the numbers
disagreed. Full account in D-036.

---

## 6. Display ties, and a limit on D-016 (D-037)

D-016 says a half-way display tie is removed, not resolved, and Phase 1
implemented that as resampling. **Resampling cannot work in this template.**
Levenspiel's trapezoid terms come from a 2-dp table, and `(y_i + y_i+1)/2` lands
on a half-cent whenever the sum is odd in its last digit — 50% of intervals by
construction, so **99.53% of instances carry at least one tie.**

The tie was an artefact of the *display*, not the data: a half-sum of two 2-dp
values is exact in 3 dp, and that times a 2-dp width is exact in 5 dp. Displayed
so, nothing rounds, so no tie can arise. The template now reports **zero T1
marginals**, from 38%.

D-037 records the general rule: before resampling a display tie, check whether
the quantity is exactly representable at a slightly longer display. Resampling
is for quantities that are irrational at any precision.

---

## 7. Verification

Corpus, measured against a `master` worktree — **not** against a recalled number:

| Check | master | Phase 2 | Δ |
|---|---:|---:|---:|
| T1 printed-arithmetic closure | 32 | 30 | −2 |
| T2 round-trip oracle | 0 | 0 | — |
| **T3 determinism** | **1** | **0** | **−1** |
| T4 output contract | 4 | 3 | −1 |
| T5 binding / rounding | 71 | 68 | −3 |
| T6 distribution | 142 | 142 | — (see §9) |
| T7 invariant asserts | 87 | 83 | −4 |

All four templates pass T1, T2, T3, T4, T5 and T7 individually. **Zero
newly-failing templates anywhere in the corpus.**

**Determinism (D2.2).** All four byte-identical across **Python 3.13.13 and
3.12.13**, 200 seeds, each run as its own process. Two versions on one machine —
the second *machine* is Reviewer A's, and is not something I can supply.

**Fallbacks (D2.3).** Ten removed, each measured over 5000 seeds before removal
(resolving to ~0.06%). Every one fired on 0.00%; two are provably unreachable
with the argument recorded. Full table in
[`phase2_item_pool_impact.md`](phase2_item_pool_impact.md) §4.

**Constants (D2.6).** C2 landed first (`9151d8c`). The C2.2 suite still passes
222 checks and C2.5 55 checks with `CP_PARAMS_COMBUSTION` added.

---

## 8. Item-pool impact (D2.5)

Full account in [`phase2_item_pool_impact.md`](phase2_item_pool_impact.md).

| Template | Questions changed | Answers changed | Median shift |
|---|---:|---:|---:|
| `levenspiel_plot_interpretation` | 300/300 | 300/300 | whole pool |
| `pfr_volume_changing_rate` | 300/300 (wording) | 79/300 | 0.041% |
| `adiabatic_flame_temperature` | 0/300 | 300/300 | 0.845% |
| `vibration_isolator_design` | 0/300 | 300/300 | 0.012% |

Levenspiel's pool changes entirely and unavoidably: it was never a function of
the seed, so there is no "before" to preserve. The flame temperatures move
because they were **wrong by 5–65 K**. The isolator's 300 answers were
additionally unparseable by the scorer, because `{:,.0f}` emits the thousands
separators that D-003 mis-reads.

---

## 9. Residual risk register

| Item | Status |
|---|---|
| **T6 baseline is stale corpus-wide** | 142/150 fail against it **on `master`**. D2.5 is delivered as a direct before/after instance diff instead. Regenerating the baseline is deliberately **not** done here — a baseline refreshed by the phase it gates is not a gate. `SPEC-CHANGE`, raised not fixed |
| two tables for four species | `CP_PARAMS` and `CP_PARAMS_COMBUSTION` can drift apart. Mitigated by a header stating who reads which and why, and a suite check that the flame template stays inside `CP_COMBUSTION_VALID_T_MAX` (D-036) |
| wide-range fit at the low end | `CP_PARAMS_COMBUSTION` is −5.4% on CO₂ at 298 K. Acceptable because the flame integral has almost no mass there (≤1.3% above 1000 K) and the eleven flame temperatures land within 1.9 K of NIST — but it is a real cost, stated |
| 5-dp trapezoid terms | heavier to read than 3 dp. The trade is exactness for legibility (D-037) |
| `levenspiel` in Phase 5 | the spec notes this template is also in Phase 5's scope for step renumbering. **Done here** — steps are now contiguous 1–8. Phase 5 should confirm rather than redo |
| second machine for T3 | not supplied by me; Reviewer A's gate |

---

## 10. What Phase 3 inherits

Phase 3 takes the two templates whose traces are *shaped* by iteration —
`normal_depth_iteration` and `line_balancing_heuristic` — and is the phase where
the schema, not the template, may be what changes.

This phase is directly relevant to it. `adiabatic_flame_temperature` **was** an
iterative solve, and it was resolved *without* a schema change by fixing the
iteration count and printing every pass. That is exactly the "alternative if the
schema route is rejected" the spec offers for §3.1, and it worked here —
but only because the iteration count could be fixed without giving anything
away. Phase 3's objection to fixing the count in `normal_depth_iteration` is
that it pins the pedagogy to "do exactly three secant updates", and in
`line_balancing_heuristic` that the station count **is part of the answer**.
Neither objection applies to a flame temperature, where nobody is being tested
on how many passes convergence takes.

**So this phase supplies a worked precedent for the cheap route and a clean
statement of its limit:** fix the iteration count when the count is not part of
what the item tests. That is a sharper criterion than "cheaper vs better", and
Phase 3 should apply it to each of its two templates separately rather than
choosing one route for both.

---

## 11. R4 triage — Reviewer A (determinism)

*Filed in [`reviews/phase2_reviewer_a_determinism.md`](reviews/phase2_reviewer_a_determinism.md).*

## 12. R4 triage — Reviewer C (numerical methods)

*Filed in [`reviews/phase2_reviewer_c_numerical.md`](reviews/phase2_reviewer_c_numerical.md).*
