# Reviewer F — triage (R4)

**Date:** 2026-09-16 · Report filed unmodified at `8a180dd` **before** any of this.

Every finding and every §4 suggestion is dispositioned. An untriaged suggestion blocks the gate
exactly as a CONFIRMED finding does.

F's verdict in one line: *"Right conclusions, unsupported evidence — which for a gate that rests on
its own honesty is the finding."* **Accepted in full.** Three of my claims are corrected below.

---

## F1 — calibration reproduces · ACCEPTED, no defect

Confirmed independently: 6166/6166, all five branches +0 at `67d41f4`, `GATE: PASS`.

**Action taken:** record the reader-trap F found. The audit carries **two** per-branch inline
tables. The gate matches `template_audit_report.md:85` (chemical 65 / civil 37 / electrical 116 /
industrial 125 / mechanical 89). `:69` gives different figures, and a reviewer checking against it
would wrongly conclude the calibration failed on all five branches. **`:85` is the authoritative
table.** The strict total 449 vs the audit's 432 is the +3.9% more-inclusive predicate disclosed at
`phase0_baseline.md:43`.

## F2 — 107 cells · ACCEPTED, no defect

## F3 — the attribution was a category error, and D6.7 did not move that cell · **ACCEPTED, MY CLAIM RETRACTED**

Commit `9c8622a` said `template_system_property_linearity`'s `78 → 50` was a D6.7 predicate
artefact. **That is wrong twice over, and F measured both.**

- **Category error.** `pct_step_values_recoverable` is computed *dynamically* — by generating
  instances and classifying `step_tokens()` through `match_value()`. The Call-node-vs-Name-node
  argument is a *static* AST fact feeding only `n_inline_computed`. It cannot move the dynamic
  column at all.
- **Refuted by measurement.** Running the pre-D6.7 source through the identical harness gives a
  byte-identical `tokens=157 recovery=50% {'MISSING': 79, 'EXACT': 33, 'SCALED': 45}`.
  **D6.7 moved that cell by zero.**

**The real mechanism** is F's: `NUM_RE` harvests *index digits out of symbolic identifiers*. The
unrecoverable lines are `T{x1[n] + x2[n]} = T{x1[n]} + T{x2[n]}`, `y1[n] = (x1[n])^3` — the 1/2/3
are subscripts and exponents, not quantities. The whole `78 → 50` gap is harness-vs-audit, and it
is **not a trace-quality statement at all**.

**What survives:** the *static* half. `{2*C}`×2 plus `{C}`×2 becoming four `signed_term(...)` Calls
gives `n_inline_computed 2 → 4` exactly as described. The helper-call → inline-count correspondence
holds; the recoverability claim does not.

## F4 — the attribution table was incomplete · **ACCEPTED, "one-to-one" RETRACTED**

D6.7 edited **ten** template modules; my table listed **six**. Omitted: `discrete_time_signals.py`
(the largest change, 11 lines — and the very file `system_property_linearity` lives in),
`continuous_time_signals.py` (7), `volumetric_properties_pure_fluids.py` (3),
`waves_and_phasors.py` (1). **"The correspondence is one-to-one" is withdrawn**: it was asserted
from a subset and F is right that a subset cannot establish it.

## F5 — the class-D claim overstated the judgement limbs · **ACCEPTED, AND NOW SETTLED BY MEASUREMENT**

F is correct that limb 2 — *"the trace's own chain does not reproduce its own answer"* — is exactly
what `t1_closure.py` checks by its own docstring, and that I filed it under judgement without
applying it.

**Applied it. All 16 class-D rows PASS T1, 0 failures each:**

`levenspiel_plot_interpretation` · `pfr_volume_changing_rate` · `adiabatic_flame_temperature` ·
`annulus_flowrate` · `cantilever_double_integration` · `normal_depth_iteration` · `mean_variance` ·
`signal_operations` · `system_properties_memory_causality` · `system_property_linearity` ·
`line_balancing_heuristic` · `incompressible_continuity` · `rotating_unbalance` ·
`vibration_isolator_design` · `vibration_transmissibility` · `damping_classification`

**So limb 2 holds none of them.** Combined with F's finding that not one row is held by the
mechanical `D:no-numeric-content` limb, the corrected statement is:

> **All 16 rows sit in class D on inherited judgement alone.** The two limbs that *can* be decided
> mechanically — limb 1 and limb 2 — clear every one of them.

That is stronger and more useful than "the gate cannot be evaluated". It also concedes F's sharper
point: **`CLASS D: 16 → 16` was a tautology of `classify_restated()`**, which opens
`if inherited_class == 'D': return 'D'` and is structurally incapable of removing a row. Presenting
it beside measured figures implied a measurement it was not.

## F6 — "34 of 34, zero hard signal" does not reproduce · **ACCEPTED, RETRACTED**

Neither artefact yields it. `--diagnose` gives **37** newly-"yes" rows, 37 prose-only. The CSV's 36
moves give **28** B-limb moves, of which **26 prose-only and 2 on a hard step-count signal**. So
*"zero hard signal"* is **false**. The corrected figure is **26 of 28 (93%)**; the rhetorical point
survives, the number does not. Corrected in `phase6_residual_register.md` and in this triage.

## F7 — nothing regressed · ACCEPTED, no action

Independently reproduced: all eight ceilings, `ci_ratchet` exit 0.

## F8 — T1 marginals unratcheted · **ACCEPTED as a residual**

97 templates sit within tolerance but on the rounding boundary (worst:
`critical_depth_froude_classification` 76%, `influence_line_max_reaction` 74%,
`upward_seepage_quick_condition` 69%). Drift there trips nothing. **Registered as the largest
unguarded surface in the suite**, with F's suggestion 5 adopted below.

---

## Suggestions (§4), triaged separately

| # | disposition |
|---|---|
| 1 | **DONE** — this document lists all ten files and replaces the `78 → 50` explanation with the symbolic-token mechanism. |
| 2 | **ADOPT-PHASE-7** — an `--at-rev` mode for the *dynamic* columns. F is right that the static/dynamic asymmetry is what produced F3: `measure_static_at_rev` exists, its dynamic counterpart does not, so "did this commit move this cell" was an argument rather than a command. Highest-leverage tooling gap found this phase. |
| 3 | **DONE** — T1 applied to all 16 D rows; per-row result published above. |
| 4 | **DONE** — 26-of-28 published, "34 of 34" retracted. |
| 5 | **ADOPTED, registered** — ratchet the T1 marginal count (97). Deferred to the same tranche as 2, since both touch the instrument rather than the corpus. |
| 6 | **REGISTERED, needs a decision** — excluding symbolic-identifier index digits from `step_tokens()` would change what the column *means* corpus-wide, not just where it looks worst. It is the right diagnosis; it is not a change to make inside a triage. |

---

## What this costs the phase

Nothing measured changes. The corpus, the ratchet, the P6 measurements and the pool all stand —
F reproduced the ceilings independently. What changes is **three claims in two commit messages**,
and the class-D story gets *better*: from "not decidable" to "both decidable limbs clear all 16;
D rests on inherited judgement".

The pattern worth carrying forward: every one of F3, F4 and F6 is the same error — **a mechanism
asserted from a subset and never A/B'd.** The measurement that would have caught each was cheap,
and in F3's case ~30 lines. That is now suggestion 2's job.
