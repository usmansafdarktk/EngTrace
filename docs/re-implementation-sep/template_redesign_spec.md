# Template Redesign — Implementation Specification

**Companion to:** [`template_audit_report.md`](template_audit_report.md), [`template_inventory.csv`](template_inventory.csv)
**Date:** 2026-09-05
**Scope:** the 16 class-D templates, plus 3 class-B/C templates that share the class-D round-trip defect, plus 11 templates with output-contract defects (Phase 5, adjacent scope). 29 distinct templates in all.
**Two tracks.** **Track A** (Phases 0–6) is template integrity. **Track B** (Phases C1–C3) is constants re-grounding, running in parallel with two hard sync points.
**Document map.** [`template_audit_report.md`](template_audit_report.md) findings · [`template_inventory.csv`](template_inventory.csv) per-template classification · **this file** the plan · [`DECISIONS.md`](DECISIONS.md) decisions and pivots, append-only · [`phase0_summary.md`](phase0_summary.md) phase close-out and R4 triage · [`phase0_baseline.md`](phase0_baseline.md) measurements · [`reviews/`](reviews/) independent reviews.
**Out of scope, referenced where it gates:** the parser fix and the value-extractor build. Those remain separate tracks in the action-item list; only their dependencies appear here.

---

## 0. Governing principles

Six rules. Every phase is judged against them.

**P1 — The trace must reproduce its own answer.** A solver who reads only the question, follows only the printed steps, and uses only the printed operands must arrive at the printed answer within the stated tolerance. This is the property most class-D templates violate and it is the reason the class exists.

**P2 — Round then recompute, and break the tie in decimal.** Any quantity that is displayed and then consumed downstream must be rounded to its display precision *before* it is consumed, so the stored value and the printed value are the same value. **That is necessary but not sufficient** (D-012): rounding to display precision removes the *double* rounding but can leave an exact half-way tie, which `round()` then resolves on the **binary** value. The tie must also be broken deterministically in decimal:

```python
delta = round(delta, 5)                                    # necessary
delta_mm = Decimal(f"{delta:.5f}") * 1000                  # and sufficient:
delta_mm = float(delta_mm.quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
#   NOT round(delta * 1000, 1) — `0.01185 * 1000` is 11.849999… in binary,
#   so round() gives 11.8 where a decimal reader gets 11.9.
```

Phase 0 found that the template originally cited here as the worked exemplar, `template_beam_deflection_formula`, applies the first half correctly and **still fails closure on 5.89% of its SI instances** for exactly this reason. Applied as originally written, Phase 1 would have reproduced that residue eight more times. Fix that template first; it is the pattern the others copy.

**P3 — Given values are stated first, then used.** If the question states a rounded value, the solution chain consumes the rounded value, never the unrounded pre-image it was derived from. Back-solving to construct a problem is permitted; back-solving and then *using the pre-image* is not.

**P4 — Determinism.** `random.seed(s)` must fully determine the output. No unseeded generators, no solver tolerances or library versions in the critical path, no machine-dependent values printed.

**P5 — No silent fallbacks.** A template either produces a conforming `(question, solution)` or raises. Blanket `except Exception` and answer-less return tuples are prohibited.

**P6 — Pedagogy is preserved unless explicitly traded.** These are benchmark items. Any change that alters what an item tests, or its difficulty, is a scoping decision that must be recorded and approved — not a side effect of a correctness fix. Every phase carries a distribution non-regression gate for this reason.

### Working conventions

- All work on branch `redesign/template-integrity`, off `master`. One PR per phase. **No phase merges before its exit gate passes.**
- No template is edited before Phase 0 is complete and its baseline is committed.
- Each phase's independent review is performed by agents that **have not seen the implementer's reasoning** — they receive the spec, the diff, and read access to the repo, and must re-derive the acceptance numbers themselves.
- **Every phase ends with a mandatory independent review** that files both findings *and* forward-looking suggestions, each of which is triaged before the phase closes. See **R0–R5**. This applies to phases with no code changes and to phases whose automated checks are all green.
- Review reports live in `docs/re-implementation-sep/reviews/`. This spec is itself under review: reviewers may return `SPEC-CHANGE` items against it.
- **Scope every review to R6 before dispatching it.** One mandatory gate task, a stated time box, no measurement commissioned twice, tooling supplied, and everything already settled named as out of scope. Over-scoping a review does not buy assurance — it stalls the review and leaves the gate unchecked. Run the R6 pre-dispatch checklist every time.

> **Measurement status — RESOLVED by Phase 0.** The defect rates quoted below came from the per-branch audit sweeps. Phase 0 re-measured all 22 testable claims independently; **20 reproduced**, several exactly. Two corrections stand: audit claim 4c is false as stated (the `damping_classification` flips are ~50/50 Overdamped/Underdamped, not all Overdamped), and the headline worst-case figures are extreme values of a 200-seed sample rather than template properties. Full reconciliation in [`phase0_baseline.md`](phase0_baseline.md); decisions in [`DECISIONS.md`](DECISIONS.md).

---

## Phase 0 — Verification infrastructure and baseline

**Nothing is edited in this phase.** Its entire purpose is to make later phases falsifiable. If Phase 0 is skipped or rushed, no later gate means anything.

### 0.1 Deliverables

| # | Deliverable | Path |
|---|---|---|
| D0.1 | Test harness, seven checks (§0.2) | `tests/template_integrity/` |
| D0.2 | Baseline snapshot: 1,000 seeded instances per template (150 × 1,000), hashed | `tests/template_integrity/baseline/` |
| D0.3 | Baseline measurement report — every audit claim re-measured | `docs/re-implementation-sep/phase0_baseline.md` |
| D0.4 | CI entry point running T1–T7 on changed templates | `tests/template_integrity/run.py` |
| D0.5 | Defect-rate reconciliation table: audit-claimed vs. re-measured, with discrepancies flagged | in D0.3 |

The audit's probe scripts (AST interpolation classifier, `sys.settrace` locals-capture harness, instance/skeleton differ) are working prototypes and should be promoted into `tests/template_integrity/` rather than rewritten.

### 0.2 The seven checks

**T1 — Printed-arithmetic closure (generic).**
For every emitted line matching `<numeric expression> = <value>`, evaluate the expression from the *printed operands* and assert it equals the printed value within ±0.5 of the last displayed digit.
*This single check catches the Phase-1 defects directly:* `rotating_unbalance` Step 1, `cantilever_double_integration` Step 6, `annulus_flowrate`'s final line, `mean_variance`'s per-term products.
Must report **coverage** (% of `=` lines parsed) and list every skipped line. Coverage below 60% on any template is itself a finding — it means the trace is mostly prose.

**T2 — Round-trip oracle (per-template).**
A hand-written `recompute(question_givens) -> answer` for each template in Phases 1–3. Asserts the answer derived *only* from stated givens matches the gold within the template's declared tolerance. T1 catches broken lines; **T2 catches broken chains** — notably `damping_classification`, where every individual line is fine and the classification is still unrecoverable.

**T3 — Determinism.**
For seeds 0–199: generate twice **in separate processes**, assert byte-identical. Separate processes is not optional — it is what catches unseeded `np.random` and hash-ordering effects. The whole corpus goes through **two child processes in total**, not two per template (SPEC-CHANGE 2). The harness must seed *only* `random` in T3's children: seeding numpy there makes an unseeded `np.random` draw reproduce across both processes and silently disarms the check.

**T4 — Output-contract conformance.**
Strict `\*\*Step (\d+):\*\*` markers, contiguous from 1, no duplicates/restarts; exactly one terminal answer marker from an approved set; no return path yielding a tuple without steps and an answer; no non-ASCII in numeric fields.

**T5 — Value binding and rounding discipline** (SPEC-CHANGE 1).
**T5a**: no quantity appearing as a **step result** may be computed inline in an f-string. Operand restatements (`{2*L}`, `{L/2}`) remain permitted and are reported separately.
**T5b**: a variable printed at N decimal places, assigned unrounded, and then consumed by a later computation is a **P2 violation** — confirmed against runtime values, so a constant whose display rounding is lossless is not flagged. T5b exists because T1 structurally cannot see this class: in `template_cantilever_double_integration` every printed line closes within its display tolerance while the trace is still wrong.

**T6 — Distribution non-regression.**
Over 1,000 seeds, before vs. after: distinct-answer count, answer magnitude quantiles (p5/p50/p95), step-count distribution, branch-coverage proportions, difficulty label. **This is the P6 guard.** Default tolerance: distinct-answer count must not fall by >10%; branch proportions must stay within ±5 points; median answer magnitude within ±2%. Any breach requires explicit sign-off, not a tolerance adjustment.

**T7 — Invariant asserts present.**
Every touched template must carry ≥2 physical/mathematical asserts in the civil/industrial style (bounds, ordering, or a conservation identity). Templates without them fail the phase.

### 0.3 Independent review (Phase 0)

One reviewer, **Harness Adversary**. Given the spec and the harness, **not** the implementer's notes. Tasks:

1. Write ≥5 deliberately broken template variants (a wrong rounding, an unseeded RNG, a non-closing printed line, an answer-less return, an inline step result) and confirm the harness **fails** each. A harness that cannot detect a planted defect cannot certify a fix.
2. Independently re-derive three audit claims from source (e.g. the 4.95% `cantilever` rate) and compare to D0.5.
3. Report any check whose implementation is weaker than its stated intent.

### 0.4 Exit gate

- [ ] T1–T7 run green on the **unmodified** repo except for the known-defective templates, and the failures are exactly the expected set
- [ ] All ≥5 planted defects detected
- [ ] D0.5 reconciliation complete; every discrepancy >20% relative either explained or the claim withdrawn
- [ ] Baseline committed and hash-stable across two machines
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 25–40 h.** *Do not compress this phase.*

---

## Phase 1 — Round-trip integrity (12 templates)

The highest-value phase. These are **defects producing incorrect gold traces in published results**, not design limitations.

> **Re-scoped by Phase 0.** The table in §1.1 is the *original* nine as commissioned. Phase 0's measurements changed both the size and the shape of this phase, and §1.1a below is authoritative. In short: three templates leave (not defective), three join (same defect shape, found by D0.5), and the remainder split into **chain breaks** and **display defects**, which need *different* fixes. See [`DECISIONS.md`](DECISIONS.md) D-011 and D-013, and [`phase0_summary.md`](phase0_summary.md) §8.

### 1.1a Revised scope — authoritative

| Category | Templates | Fix | Measured |
|---|---|---|---|
| **Chain break** — answer wrong from the stated givens | `mean_variance`, `rotating_unbalance`, `vibration_transmissibility` | **P3** (state, then recompute forward) | 56.8% / 33.7% / 39% of instances |
| **Display defect** — answer round-trips, an intermediate line does not | `beam_deflection_formula` ← **fix first**, `cantilever_double_integration`, `annulus_flowrate`, `statically_indeterminate_shaft`, `shaft_design_power`, `composite_shafts_series` | **P2 as amended** (decimal tie-break) | 5.89% / 3.7% / 94% / 47% / 30% / 89.9% |
| **Ill-posed, not wrong** | `damping_classification` | state `c` to more digits; drop `elif zeta == 1` (D-010) | 33.7% flips, all within stated precision |
| **Leaving Phase 1 — not defective** | `poissons_ratio`, `logarithmic_decrement`, `system_properties` | none (record in D6.6) | 0% round-trip failure |

**Verify fixes against T5b, not T1 alone** — T1's ±0.5-ulp rule structurally cannot see the half-way-tie class (adversary F3).

Three defects found during Phase 0 that must be fixed alongside: `mean_variance` probabilities are never renormalised (sum 0.999–1.002); `rotating_unbalance` never states whether the unbalance mass is included in the stated total; and some `rotating_unbalance` answers print to one significant figure (`0.002 mm`), ungradeable regardless of the chain fix.

### 1.1 Original scope (superseded by §1.1a — retained for provenance)

| Template | Branch | Class | Defect | Claimed rate |
|---|---|---|---|---|
| `template_rotating_unbalance` | mech | D | RPM stated rounded, ω used unrounded | 7.77% answer error |
| `template_vibration_transmissibility` | mech | D | same construction | 0.87% TR error |
| `template_damping_classification` | mech | D | ζ unrecoverable from stated `c`; `zeta == 1` float equality | 175/500 class flips |
| `template_cantilever_double_integration` | civil | D | `delta` full precision → Step 6 ×1000 | 4.95% |
| `template_annulus_flowrate` | chem | D | ΔP back-solved from hidden `Re_target`, re-rounded | final line does not close |
| `template_mean_variance` | elec | D | per-term values built in a comprehension inside an f-string | printed mean does not recompute |
| `template_poissons_ratio` | mech | C | load back-solved from hidden `target_strain`; sign leak | 76/300 sign mismatch |
| `template_logarithmic_decrement` | mech | B | back-solved from hidden `zeta_actual` | by design — see 1.3 |
| `template_system_properties` | mech | C | `c` stated rounded, ζ printed from pre-image | 4.5e-6 |

### 1.2 Required change

Uniform transformation, applied per template:

```python
# BEFORE — pre-image consumed downstream
omega = freq_ratio_r * omega_n
operating_speed_rpm = omega * 60 / (2 * math.pi)
# question states round(operating_speed_rpm, 0); solution uses unrounded omega

# AFTER — P3: state, then recompute forward from the stated value
operating_speed_rpm = round(freq_ratio_r * omega_n * 60 / (2 * math.pi), 0)
omega = round(operating_speed_rpm * 2 * math.pi / 60, PRECISION)   # the value the solver can reach
```

Every downstream quantity then derives from `omega`, not from `freq_ratio_r`. The sampler may still *choose* `freq_ratio_r` to place the item in a useful regime — it simply may not leak it into the solution.

### 1.3 Two cases needing a real design decision

**`damping_classification` — critically damped is not reachable by rounding.** After P3, ζ recomputed from a 2-dp `c` will essentially never equal exactly 1.0, so the "Critically Damped" class becomes unreachable and `elif zeta == 1` is dead.

**Decision: constrain the sample (option 1). Feasibility verified — see below.** The two alternatives (a stated band `|ζ − 1| ≤ 0.005`, or dropping the class) are recorded as rejected: the first changes the item, the second removes ~1/3 of its answer space and its most interesting case.

*Construction.* With `M = 100·m` (integer, since `m` is 2 dp) and `k` integer, `c_c = 2√(km) = √(kM)/5`. If `kM` is a perfect square `N²` then `c_c = N/5 = 0.2N`, **exact at 1 dp**. Sample `k = q·a²`, `M = q·b²`. Within the *current* ranges (`m ∈ [2, 600]`, `k ∈ [2000, 300000]`) this yields **680,869 feasible (k, m) pairs — 8,627 distinct stiffnesses, 3,944 distinct masses** — with entirely natural-looking values (`k = 27455 N/m, m = 160.55 kg → c_c = 4199.0`). Diversity is not a concern.

*Refined design — this also removes the float-equality test entirely.* Invert the current construction: sample `c` **directly** as an exact 2-dp value rather than deriving it from a chosen ζ. Then

- classification is an **exact comparison of two exact decimals**, `c` vs `c_c` (scaled-integer compare — the industrial branch's technique). `elif zeta == 1` disappears.
- ζ = `c / c_c` is computed for **display only** and never gates the classification.
- all three classes stay reachable and exactly representable.

This satisfies P1–P3 without any pedagogical change: the question still states `m`, `k`, `c` and asks for ζ, the class, and ω_d.

**`logarithmic_decrement` — the mismatch is deliberate.** The template back-solves amplitudes from a hidden `zeta_actual`, rounds them to 1 dp, and the recovered ζ *intentionally* differs — that is the measurement-uncertainty lesson. **This is not a defect.** Required change is only that the tolerance be **declared** in the template so the verifier knows the intended slack, not that the behaviour change. Flagged here so Phase 1 does not "fix" a feature.

### 1.4 Deliverables

| # | Deliverable |
|---|---|
| D1.1 | 8 templates edited to satisfy P1–P3 (`logarithmic_decrement`: declared tolerance only) |
| D1.2 | T2 round-trip oracles for all 9 |
| D1.3 | Per-template `TOLERANCE` declaration (value + justification comment), for the verifier to consume later |
| D1.4 | Decision record for `damping_classification` (§1.3) with the option chosen and why |
| D1.5 | T6 distribution diff, before/after, all 9 |
| D1.6 | **Item-pool impact note**: which published results are invalidated by regenerating these templates |
| D1.7 | Phase report: defect rate before → after (target: 0) |

### 1.5 Independent review

**Two reviewers, run in parallel, neither seeing the implementer's reasoning.**

*Reviewer A — Correctness.* Given the diff and the spec: verify P1–P3 hold; independently sample 200 seeds per template and confirm the round-trip; attempt to construct a seed that still breaks closure; confirm the T2 oracles are genuinely independent re-derivations and not copies of the template's own arithmetic. **This last point is the main failure mode of this phase** — an oracle that reuses the template's expression proves nothing.

*Reviewer B — Physics & pedagogy.* Confirm each item still tests what it tested: governing equations unchanged, parameter ranges still physically sensible, difficulty label still justified, and that the `damping_classification` decision does not degrade the item. Explicitly asked to answer: *would a domain expert notice this item got easier or less interesting?*

Both file written verdicts with reproduction commands. **Any CONFIRMED finding blocks the merge.**

### 1.6 Exit gate

- [ ] T1 closure 100% on all 9 (no non-closing printed line at any of 1,000 seeds)
- [ ] T2 round-trip passes at 1,000 seeds each
- [ ] T3 determinism passes
- [ ] T6 within tolerance, or breach explicitly signed off
- [ ] T7 ≥2 asserts each
- [ ] Reviewer A and B both clear
- [ ] D1.6 filed
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 15–25 h implementation + 8–12 h review.**

---

## Phase 2 — Determinism and solver-in-the-loop (4 templates)

### 2.1 Scope

| Template | Defect |
|---|---|
| `template_levenspiel_plot_interpretation` | **Unseeded `np.random`** — the recorded seed cannot regenerate the item. Plus 6 blanket `try/except` that silently substitute a different formula while the trace still prints the original one; step numbering restarts `1,2,3,1,2,3,4,5,6` |
| `template_adiabatic_flame_temperature` | `scipy.fsolve` — milestone only tolerance-reproducible; blanket `except (…, Exception)`; `while/else` fallback returns a tuple with no steps and no answer |
| `template_pfr_volume_changing_rate` | Prints the integrator's own error estimate as a step value (machine-dependent); `n` rounds to exactly 2.0 in ~10% of draws, contradicting the printed "non-integer order" note |
| `template_vibration_isolator_design` | Quadratic root **selection** (`max` of two roots) not expressible as a formula; Step 4 **hardcodes** its verdict without testing; two answer-less error returns; `{:,.0f}` breaks the parser |

### 2.2 Required changes

**Blocking dependency:** `adiabatic_flame_temperature` must not be re-baselined until `CP_PARAMS` is corrected on the constants track — its current output (~1130–1493 K where methane/air should be ~2200 K) is wrong for a different reason, and fixing determinism against wrong data wastes the work. **Sequence the constants correction for this table before Phase 2 starts.**

1. **Seed all generators.** `np.random.seed()` alongside `random.seed()`, or replace numpy sampling with stdlib. Enforced by T3.
2. **Delete every blanket `except`.** Per P5, a failure raises. Where a fallback exists to handle a genuinely rare sampling failure, replace with a bounded resample loop that raises on exhaustion — the pattern already used throughout civil/industrial.
3. **Make the solver step reproducible.** Replace `fsolve` with a *stated* algorithm: fixed iteration count, stated initial guess, stated convergence rule, all printed. The trace becomes something a student can reproduce by hand. Where a closed form exists, prefer it.
4. **Remove machine-dependent printed values** — the integration error estimate goes.
5. **Fix `n == 2.0`** — resample to exclude, or reword the note.
6. **Renumber Levenspiel's steps** contiguously (T4).
7. **`vibration_isolator_design`:** make root selection an explicit predicate step that *states and evaluates* the `r > √2` condition rather than asserting "Yes … met" unconditionally; replace the two answer-less returns with raises; drop the `{:,.0f}` format.

### 2.3 Deliverables

D2.1 four templates edited · D2.2 T3 determinism proof, 200 seeds × 2 processes · D2.3 written justification for each removed `try/except` (what it was masking, why removal is safe) · D2.4 stated-algorithm spec for the replaced solver, with a hand-worked example · D2.5 T6 distribution diff — **expected to move**, since removing silent fallbacks changes the accepted sample set; movement must be explained and approved · D2.6 confirmation that `CP_PARAMS` landed first

### 2.4 Independent review

*Reviewer A — Determinism.* Re-run T3 on a **different machine and Python version**. Grep for every non-deterministic source (unseeded RNG, `set`/`dict` iteration order affecting output, `time`, `id()`, floating-point-dependent branching). Confirm no blanket excepts remain.

*Reviewer C — Numerical methods.* Given the replaced solver: confirm the stated algorithm converges over the full sampled parameter range, that the stated iteration count is sufficient at the extremes, and that the printed intermediate values are reproducible by hand. Verify `adiabatic_flame_temperature` now produces physically plausible flame temperatures against an external source.

### 2.5 Exit gate

- [ ] T3 passes on two machines / two Python versions
- [ ] Zero blanket `except` in scope; zero answer-less return paths
- [ ] T4 contract conformance (Levenspiel numbering)
- [ ] Solver replacement converges across the full range (Reviewer C)
- [ ] `adiabatic_flame_temperature` physically plausible post-`CP_PARAMS`
- [ ] Reviewers A and C clear
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 12–24 h implementation + 8–12 h review.** Gated on the constants track.

---

## Phase 3 — Trace shape: iteration and search (2 templates)

The only phase where the **schema**, not the template, may be what changes.

### 3.1 `template_normal_depth_iteration`

A secant solve inside the trace. `n_steps` is constant 4, but Step 3 holds **8–28 quantities**, and `y_curr`/`g_curr` are rebound every pass — no stable `{id, symbol}` for "the third trial depth."

**Recommended: schema route.** Iterative solution *is* the skill this item tests; flattening it to three fixed updates would gut it. Introduce an `iteration` node type: an ordered list of homogeneous sub-traces, each with stable within-iteration symbols (`y_k`, `g_k`, `A_k`, `P_k`) and an explicit convergence predicate.

Template-side changes are then minimal but non-zero: name the loop variables per iteration rather than rebinding; guard the unguarded `g_curr - g_prev` division (a `ZeroDivisionError` is reachable in principle, 0 occurrences in 6,000 seeds); remove the `fmt3()` display-only `-0.000` fixer, which currently lets the printed operand differ in sign from the stored value (a direct P2 violation).

*Alternative if the schema route is rejected:* fix the iteration count at 3 with stable names — cheaper (~4 h), makes it class A, but pins the pedagogy to "do exactly three secant updates." **This is a P6 trade and needs sign-off.**

### 3.2 `template_line_balancing_heuristic`

The trace is a **search log**: milestones are task *sets* (`Station 1 contents: {a, b}`), a running remaining-time counter, and boolean fit tests, with one step per opened workstation (5 or 6).

**Recommended: schema route, same as §3.1.** *(This revises an earlier recommendation of "redesign to a fixed-shape decision table" — that option is now rejected.)* The reason is the same objection that rules out fixing the iteration count in §3.1: **the station count `n` is a computed result, not an input.** Line efficiency is a function of `n`, so fixing the station count to make the trace a fixed shape gives away part of the answer and materially weakens the item.

Instead, add a `decision` node type: an ordered list of per-station records `{station_id, assigned_set, remaining_time, rejected_candidates}`, count not fixed. Verification is set equality plus a replay of the greedy rule — deterministic.

Template-side changes are small: name the per-station accumulators stably rather than rebinding, and emit the eligibility snapshot as structured data rather than only as prose.

**Note that §3.1 and §3.2 now share one solution family** — `iteration` and `decision` are the same shape (an ordered list of homogeneous sub-traces with stable within-element symbols and a termination predicate). Specify them together; that is cheaper than two bespoke redesigns and it generalises to `linear_reservoir_routing_step` (a repeated sub-chain unrolled into the trace) and `qr_policy_one_iteration` (iterative in principle, one iteration emitted).

> **AMENDED by SPEC-CHANGE 8 (D-038).** The structural claim above is correct and the **equivalence is not**. The two types differ on whether the sequence's cardinality is an observable of the answer — an update count is incidental and is not even stable under the numerical slack the comparator tolerates, while the station count **is** the answer and is exact — so a comparator built from one merged type is wrong on one of the two templates. They share a base structure, a `carry` mechanism and a verification algorithm; they do **not** share a comparator. The specification delivered is [`phase3_node_types.md`](phase3_node_types.md), schema 1.5. The generalisation to `linear_reservoir_routing_step` and `qr_policy_one_iteration` remains **argued, not measured**, and is now a Phase 4 deliverable (D4.7).

### 3.3 Deliverables

D3.1 decision record: schema route vs. redesign for each, with the P6 trade stated · D3.2 template edits · D3.3 `iteration` and `decision` node type specifications, with worked examples from both templates — **this is the primary deliverable and feeds directly into the milestone-model design** · D3.4 comparator semantics: how a model's iteration count differing from gold is scored (a correct answer reached in 4 iterations instead of 3 must not be marked wrong) · D3.5 T6 distribution diff

### 3.4 Independent review

*Reviewer B — Pedagogy.* Does either item still test what it tested? Specifically: does the `line_balancing` redesign still require the solver to *run* the heuristic, or has it become a lookup?

*Reviewer D — Schema.* Given only the D3.3 node specs and a set of generated traces: attempt to write a verifier against the spec. Report every ambiguity. **The spec passes only if a reviewer who did not design it can implement against it** — this is the real test of D3.3, and it is the deliverable most likely to fail its gate on first attempt.

### 3.5 Exit gate — **PASSED 2026-09-06**, merged `27054c9`

Close-out in [`phase3_summary.md`](phase3_summary.md); the gate table there is authoritative.

- [x] D3.1 decisions recorded and signed off — schema route for both, argued separately
- [~] T1–T7 pass on both — **T1–T5 and T7 pass; T6 does not and cannot**, the committed baseline being stale corpus-wide (142/150 on `master`). Not regenerated: a baseline refreshed by the phase it gates is not a gate (D-043). Replaced by a direct before/after distribution diff. Phase 6 owns the regeneration
- [x] `g_curr - g_prev` guarded; `fmt3()` sign discrepancy removed — both were **latent**, and are recorded as latent rather than claimed as live fixes (D-042)
- [x] Reviewer D implemented a working verifier from D3.3 alone — twice, independently: 80/80 + 33/33, then 80/80 + 32/32 from a reviewer that never saw the first
- [x] D3.4 comparator semantics cover the iteration-count-mismatch case
- [x] **Independent review filed (R2) and every §5 suggestion triaged (R4)** — 28 findings and 27 §5 suggestions across **six** rounds, all disposed

**Carried out of scope, recorded not fixed:** `line_balancing_heuristic` is ~90% shortcuttable from its question text (D-046), pre-existing and measured identically on `master`.

**Effort: 8–16 h implementation + 10–14 h review** (higher review share — the deliverable is a specification).

---

## Phase 4 — Non-numeric verifiable templates (4 templates, no code edits)

`template_system_properties_memory_causality` · `template_system_property_linearity` · `template_signal_operations` · `template_incompressible_continuity`

### 4.1 Position

**These templates are not defective and must not be redesigned.** They have zero numeric content, but they are all **deterministically verifiable by a non-arithmetic comparator**. Adding a numeric payload to make them fit a numeric checker would damage four sound conceptual items to satisfy an implementation convenience. The schema adapts, not the items.

| Template | Comparator | Determinism |
|---|---|---|
| `system_properties_memory_causality` | Canonical label tuple (`Memoryless: No`, `Causal: Yes`) | Exact match |
| `system_property_linearity` | Canonical label (`linear` / `not linear`) | Exact match |
| `signal_operations` | Sequence-with-origin equality — values **and** the `n=0` index | Exact match |
| `incompressible_continuity` | Symbolic equivalence (sympy), including the discarded arbitrary `f(x)` | CAS-deterministic |

### 4.2 Required work — all specification

1. Define milestone `kind ∈ {numeric, categorical, sequence, symbolic, narrative, check}` and the comparator contract for each.
2. **Label normalisation** for categorical: case, synonyms (`non-linear`/`nonlinear`/`not linear`), hedging ("appears to be linear"). Must be specified from real model outputs, not invented — draw the vocabulary from the 2,200 archived traces in `error_analysis_annotation/samples/`.
3. **Sequence-with-origin format** and its tolerance for presentation variation (braces vs. brackets, asterisk vs. arrow origin marker).
4. **Symbolic equivalence policy**: what counts as equivalent, how the arbitrary integration function is handled, sympy timeout/failure behaviour.
5. Fix the cosmetic `New Value y[k')` header bug in `signal_operations` (mismatched bracket, wrong variable) — the one code change in this phase.

These four comparators also serve the 51 class-C templates. **Phase 4 is where most of class C's design cost is actually paid** — treat it as leverage, not a footnote.

### 4.3 Deliverables

D4.1 comparator specification for all six `kind` values · D4.2 categorical label-normalisation vocabulary, derived from archived model outputs with frequency counts · D4.3 reference implementation + unit tests per comparator · D4.4 adversarial test set: ≥20 hand-written near-miss model answers per comparator (correct-but-phrased-differently, and wrong-but-similar) · D4.5 the one cosmetic fix

**Added by Phase 3's R4 triage** (`phase3_summary.md` §11.4, §11.5, §12.5). Each was a reviewer suggestion dispositioned `ADOPT-PHASE-4`, so it is a named deliverable here rather than a note:

- **D4.6 — a conformance corpus for D3.4 §7.** *The largest gap Phase 3 left*, raised independently by both schema reviewers. Every verifier built in Phase 3 implements §6 (gold checking); the comparator rules that actually decide a model's score are **unexercised, with no corpus at all**. Needs candidate traces: wrong count with right answer, right count with different packing, early stop, fabricated frames.
- **D4.7 — fit a third template to the `iteration` node type by declaring a binding only**, no verifier edit permitted (`linear_reservoir_routing_step` or `qr_policy_one_iteration`). Turns §3.2's "argued, not measured" into measured. A *synthetic* renamed node already passes; a real one has never been tried.
- **D4.8 — fit a second `decision`-shaped template.** Until then `decision` is a type by construction and by the rename test, but on one instance and one precedence DAG.
- **D4.9 — classify trace-length semantics corpus-wide, once.** Apply D-038's `incidental` / `answer_bearing` test as a survey rather than rediscovering it per template.
- **D4.10 — operationalise "difficulty unchanged"** so a pedagogy reviewer has something to measure rather than reason about (a required-inference-count or step-count proxy). Phase 3's Reviewer B could only argue this.
- **D4.11 — generate §8-style reject lists from the verifier** rather than writing them alongside it. A hand-written reject list is a floor: Phase 3's was "a list of remembered failure modes" and five guessed invariant holes all hit.

### 4.4 Independent review

*Reviewer E — Comparator adversary.* Attempts to break each comparator: find a **correct** model answer it rejects, and an **incorrect** one it accepts. Draws from the real archived traces, not invented examples. Reports precision/recall on D4.4 plus their own additions.

*Reviewer B — Pedagogy.* Confirms no item was weakened, and that the normalisation vocabulary does not accept answers a human grader would reject (e.g. accepting "linear" when the model reached it by invalid reasoning is out of scope here, but accepting a hedge that never commits is not).

### 4.5 Exit gate

- [x] All six `kind` comparators specified and implemented, and the `multipart`
      question answered (SPEC-CHANGE 12 / D-047)
- [x] ≥95% precision and recall on D4.4 — **97.2% / 98.6%**, with two declared
      false accepts (`num-16b`, D-052; `cat-16`, D-056)
- [x] Reviewer E found no false accept **on real archived traces** — 0 across
      four rounds. **Read this with E's F0 and R2-F10**: the archive holds
      16 wrong answers, 15 decided, all in two of the six kinds, and the
      commitment machinery has no archive instances at all. The gate item is
      met and it is uninformative for four kinds; that is stated, not hidden.
- [x] Normalisation vocabulary traced to observed outputs, with **per-template**
      coverage and resolution limits, never averaged (D-048)
- [x] **Independent review filed (R2) and every §5 suggestion triaged (R4)** —
      two reviewers, **four rounds each**, 34 findings, all dispositioned

**Effort: 12–20 h + 8–10 h review.** No template redesign.

> **Actual: four review rounds per reviewer, not one.** The comparator contract is a natural-language surface, and Reviewer E's round-4 characterisation is that it is *unbounded* — every round's findings were one- or two-token mutations of cases the suite already passed. The phase closed by **reducing the mechanism's scope** (SPEC-CHANGE 11), not by exhausting the findings. Budget later natural-language work the same way: the stopping rule has to be an evidence threshold, not a clean round.

---

## Phase 5 — Output-contract hygiene, and comparator bindings (two tracks)

### Phase 5, Track A — Output-contract hygiene (11 templates, adjacent scope)

Not class-D, but cheap, and **every structural migration trips over these**. Fold in here or run in parallel with Phase 1.

Counts below were verified against source, not inferred from the archive.

| Defect | n | Templates | Fix |
|---|---:|---|---|
| Malformed `**Step N:**` (`**Step 2: **`, `**Step 3: **…*`, colon inside the bold) | 3 | `cd_dc_system_analysis`, `euclidean_distance_binary`, `finite_convolution` | Normalise the marker |
| No `**Answer:**` — uses `**Final Answer**` / `**Final Answers:**` | 5 | `batch_moles_vs_conversion`, `flow_system_molar_flow_rates`, `gas_phase_concentration`, `limiting_reactant`, `levenspiel_plot_interpretation` | Normalise, **or** widen the accepted marker set — decide once and apply corpus-wide |
| `+ j-51.22` malformed complex | 2 | `time_to_phasor`, `phasor_addition` | Sign formatting |
| `omega_a = 0*pi` from an unreachable `elif` | 1 | `decimation_aliasing_analysis` | Reorder the guard |

**11 distinct templates.** `levenspiel_plot_interpretation` is also in Phase 2 (step renumbering) — coordinate the two edits or sequence Phase 5 after Phase 2 for that file.

`cd_dc_system_analysis` is the urgent one: a strict marker parser drops **Step 3 in 100% of instances**, and Step 3 computes the answer.

**Deliverables:** D5.1 edits · D5.2 corpus-wide marker-conformance scan, all 150, proving zero remaining violations · D5.3 the marker-set decision record.
**Review:** *Reviewer A* re-runs T4 across all 150 templates, not just the 8 touched.
**Exit gate (Track A):** T4 passes on 150/150; T6 unchanged (these are formatting-only, so any distribution movement indicates an unintended change).
**Effort: 6–10 h + 3–4 h review.**

---

### Phase 5, Track B — Comparator bindings and the instruments that validate them

**Added after Phase 4 merged**, from D-058 and Phase 4's residual-risk register.
**This is a separate track with a separate gate and a separate reviewer.** It
shares a phase number with Track A and nothing else: Track A edits templates and
is gated by T4; Track B edits no template at all and is gated by cross-pairing.
Running them under one gate would be the R6 kitchen-sink anti-pattern, and the
brief must not merge them.

**The finding that creates the track.** Phase 4 delivered the comparator contract
and bound **4 of 150 templates** to it. Gold×gold cross-pairing over all 150
found three real defects the phase's own corpora could not see — because
`numeric` and `check` had *zero* archived negative instances (Reviewer E, F0):

| # | Defect | Evidence |
|---|---|---|
| **N1** | `parse_number` reads the **first** number in the answer span, which is often a grade, a temperature or a formula subscript | confirmed false accept on `hagen_poiseuille_flowrate`; 3 scalar templates in a 12-instance sample; **58 of 150 exposed** |
| **N2** | an **empty parse compares equal to an empty parse** — `{} == {}` is a `MATCH` | 128 of 132 pairs on `autocorrelation_rect_pulse` |
| **N3** | **no bindings exist for 146 templates** — no label set, no unit, no answer shape | `euclidean_distance_binary` matches 132/132 |

**N3 is larger than N1 and N2 together.** The comparator does not need more
rules; it needs bindings, and each binding needs evidence that it discriminates.

| # | Deliverable |
|---|---|
| **D5.4** | **Gold×gold cross-pairing as a standing check.** Every kind, no archive needed, truth known without labels. 19,668 pairs over 150 templates today. |
| **D5.5** | **Archive×gold cross-pairing**, Reviewer E's round-1 instrument promoted from review artefact to standing check. **22,982 real-text pairs** available. |
| **D5.6** | **Fix N1 by measurement, not by example.** Implement ≥3 candidate extraction rules (last number; number adjacent to a declared unit; `UNRESOLVED` on ambiguity) and **score each against D5.4/D5.5**, then adopt the winner. Choosing on one example is the error shape Phase 4 committed six times. |
| **D5.7** | **Fix N2**: an unrecoverable parse is `UNRESOLVED`, never `MATCH`, on every code path — including the ones that never reach the CAS. |
| **D5.8** | **Bindings for the corpus, validated per binding.** Declare `kind`, label sets, units and answer shape per template, and require each to clear a cross-pairing threshold before it counts as bound. Start with the 5 `classification`, 9 `symbolic`, 15 `vector`/`array` and 32 `multipart`; the 89 `scalar` need only a unit. |
| **D5.9** | **Per-kind negative-instance counts reported beside precision** (Reviewer E, §5: *"precision without it is unreadable"*). |
| **D5.10** | **Declare units for the 87 templates whose gold answer carries one**, closing D-052 for the bulk of the corpus and leaving the residual named. |
| **D5.11** | Phase 4 residuals: **R4-11** (`require_origin` is route-3, 0 of 16 archived traces reach its branch), **RB4-1** (the bare-comment proxy treats expletive `it`/`there` as anaphors), **R4-12** (the negative frames' nearest neighbours). |

**Review (Track B):** *Reviewer E — comparator adversary*, whose single mandatory
task is the gate property below. **Do not re-commission Track A's T4 sweep from
this reviewer.**

**Exit gate (Track B):**
- [ ] D5.4 and D5.5 land as runnable standing checks, with their pair counts reported
- [ ] **Zero false accepts on gold×gold across all 150 templates**, or every remaining one named with a written reason
- [ ] N1 fixed by a rule **chosen on measured evidence**, with the losing candidates' numbers recorded
- [ ] N2 fixed: no code path returns `MATCH` from an unrecoverable parse
- [ ] Per-kind negative counts published beside every precision figure
- [ ] Bindings declared and cross-pair-validated for every `classification`, `symbolic`, `vector`, `array` and `multipart` template
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 15–25 h + 8–10 h review.** No template redesign; Track B edits
`tests/comparators/` and `template_inventory.csv` only.

---

## Phase 6 — Consolidation and re-audit

### 6.1 Deliverables

D6.1 **Full re-run of the original audit harness** across all 150 templates; regenerate `template_inventory.csv` with post-change classifications · D6.2 class-migration table (D → A/B/C, with the 4 Phase-4 templates recorded as *verifiable, non-numeric* rather than resolved-by-edit) · D6.3 updated `template_audit_report.md` with post-change figures · D6.4 **item-pool impact statement**: every published number invalidated by regeneration, consolidated across phases · D6.5 regression suite promoted to CI on `master` · D6.6 residual-risk register: everything knowingly not fixed, and why

### 6.2 Independent review

*Reviewer F — Audit replication.* Given only the original audit report and the changed repo, **independently re-runs the audit** and compares. Must confirm: the class-D count fell as claimed, no new class-D templates were introduced elsewhere, and no metric regressed silently. Explicitly tasked with looking for **collateral damage** — a Phase-1 rounding change that broke a Phase-5 template's contract, or similar.

### 6.3 Exit gate

- [ ] Class D reduced from 16 to ≤4 (the Phase-4 four, reclassified not fixed)
- [ ] No template regressed in class
- [ ] T1–T7 green corpus-wide
- [ ] Reviewer F's independent re-audit matches D6.3 within stated tolerance
- [ ] D6.4 filed and circulated before any paper revision uses the regenerated items
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 10–16 h + 8–10 h review.**

---

---

# Track B — Constants re-grounding

Everything above is **Track A** (template integrity). Constants re-grounding is **Track B**: it runs in parallel, not as a Track A phase, because only a small slice of it is on Track A's critical path. Two hard sync points bind them.

**Why parallel and not folded in.** Track B is ~90–150 h — comparable to all of Track A — and only one file's worth of it blocks a template phase. Sequencing all of Track A behind all of Track B would delay Phase 1 (which fixes *demonstrably wrong gold traces*) by months for no benefit.

**Why it cannot simply be deferred either.** Changing a constant changes generated items, exactly as changing a template does. If Track B lands after the item pool is regenerated, you regenerate twice. Sync point **S2** exists to prevent that.

### The gap being closed

| Branch | Citation-bearing comments | Top-level constant tables |
|---|---:|---:|
| chemical | 2 | 19 |
| electrical | 1 | 14 |
| mechanical | 1 | 10 |
| civil | 34 | 32 |
| industrial | 39 | 39 |

Civil and industrial already carry a workable provenance vocabulary — `[ON-DISK]` (transcribed verbatim, page cited), `[VERIFY: <source>]`, `[DERIVABLE]`, `[REALISM]`, `[POLICY: sampling-only]` — plus a **given-values rule**: a value restated in the question text is correctness-neutral and needs a plausibility window, not a citation. **Adopt both rather than inventing a new scheme.**

---

## Phase C1 — Sources and convention

**Deliverables.** C1.1 acquired reference set, with edition and page ranges recorded · C1.2 provenance-tag convention documented and applied to one pilot table end-to-end · C1.3 **classification of all 43 tables** in the three original branches into `[ON-DISK]` (needs a citation) vs `[POLICY: sampling-only]` (needs only a plausibility window), applying the given-values rule · C1.4 a machine-readable unit declaration per table — the audit found **no unit field anywhere** in any branch; units live in identifier suffixes and comments.

C1.3 is the deliverable that sizes the rest of the track, and it will shrink it substantially. Provisional split from the audit:

| Branch | Needs real sourcing | Sampling-only / trivial |
|---|---|---|
| **chemical** | ~250 values in 10 tables | reactant/product name lists |
| **mechanical** | ~200 values in 6 tables | `OBJECT_SHAPES`, `OBJECT_MATERIALS` (60 flavour strings) |
| **electrical** | **~23 values** — `C0`, `EPSILON_0` (CODATA), `MEDIA_VELOCITIES` (21) | 11 of 14 tables are sampling ranges |

**Review.** *Reviewer G — Provenance.* Confirms the sampling-only classification is not being used to avoid work: spot-checks 10 tables classified `[POLICY]` and verifies each really is restated in the question text for every template that consumes it.

**Exit gate.** All 43 tables classified; pilot table complete; unit declarations specified. **Independent review filed (R2) and every §5 suggestion triaged (R4).**
**Effort: 12–20 h.**

---

## Phase C2 — Chemical thermochemistry ← **critical path**

**This is the slice that gates Track A Phase 2**, and it is far smaller than the full track.

| Table | Values | Status |
|---|---:|---|
| `CP_PARAMS` | 33 | **Known defective** — B and C coefficients 10× too large (`E-2`/`E-5` where Smith–Van Ness gives `E-3`/`E-6`); Cp(N₂, 300 K) = 42.4 vs 29.1 J/mol·K |
| `HEATS_OF_FORMATION` | 21 | Unverified |
| `COMBUSTION_REACTIONS` | 11 | Unverified |

**65 values, three tables, and all three are consumed by exactly one file** — `chemical_engineering/thermodynamics/heat_effects.py`, 5 templates. Tight, well-bounded, and independently verifiable.

**Sources.** Smith–Van Ness–Abbott (already the implied source for `CP_PARAMS`), cross-checked against the NIST Chemistry WebBook.

**Deliverables.** C2.1 the three tables re-derived with `[ON-DISK]` citations to edition and page · C2.2 **physical-plausibility test suite** — Cp at reference temperatures against literature, adiabatic flame temperatures against published values for each of the 11 reactions · C2.3 before/after impact statement for the 5 affected templates · C2.4 confirmation that `template_adiabatic_flame_temperature` now lands near ~2200 K for methane/air (currently 1130–1493 K)

**Review.** *Reviewer H — Domain (thermochemistry).* Independently re-derives a 10-value sample directly from the cited source pages, without reference to the implementer's transcription. Confirms C2.2 targets are correct physics, not just self-consistent.

**Exit gate.**
- [ ] All 65 values carry a citation to edition + page
- [ ] C2.2 plausibility suite green
- [ ] Methane/air flame temperature within the literature range
- [ ] Reviewer H's independent sample matches
- [ ] **Independent review filed (R2) and every §5 suggestion triaged (R4)**

**Effort: 15–25 h.** → **Sync point S1: this must land before Track A Phase 2 begins.**

---

## Phase C3 — Remaining tables

Everything else, in parallel with Track A Phases 1, 4 and 5.

- **Chemical (~185 values, 7 tables):** `CRITICAL_PROPERTIES`, `REAL_FLUID_DATA`, `THERMO_SUBSTANCES`, `GAS_MOLECULAR_PARAMS` (Lennard-Jones σ and ε/k), `COMMON_LIQUIDS`, `COMMON_GASES`, `POWER_LAW_FLUIDS`. Sources: Perry's, Bird–Stewart–Lightfoot App. E, Poling–Prausnitz–O'Connell.
- **Mechanical (~200 values, 6 tables):** `MATERIAL_PROPERTIES`, `SHEAR_MODULUS_VALUES`, `MATERIAL_DENSITIES`, `FLUID_DENSITIES`, `PIPE_FLUIDS`, `MANOMETER_FLUIDS`. Sources: Hibbeler / Beer–Johnston appendices, ASM Metals Handbook, White or Çengel–Cimbala.
- **Electrical (~23 values):** CODATA 2022 for `C0` and `EPSILON_0`; Balanis or Ulaby for `MEDIA_VELOCITIES`.
- **Sweep civil and industrial too.** Industrial's `constants.py` carries an explicit *"PROVENANCE CAVEAT — UNVERIFIED"* block and a logged withdrawn citation; the audit also found the `SCS_CURVE_NUMBERS` A/B/C/D ordering exists only in a comment and an inline dict literal. Being tagged is not the same as being verified.

**Two known data-integrity defects to resolve here**, both found in the audit: `template_two_phase_specific_volume` names a real substance from `THERMO_SUBSTANCES` but **invents** V_l/V_v randomly, ignoring `REAL_FLUID_DATA` (Ammonia v_g 0.754 vs true 0.1284 m³/kg); and `template_floating_object_submersion_depth` has a hardcoded fallback injecting "Pine Wood" at 500 kg/m³, which is not a `MATERIAL_DENSITIES` entry.

**Deliverables.** C3.1 re-derived tables with citations · C3.2 plausibility suite per table (order-of-magnitude and cross-property consistency checks) · C3.3 the two data-integrity fixes · C3.4 impact statement — which templates' answer distributions move, and by how much · C3.5 residual `[UNVERIFIED]` register for anything that could not be sourced, **left explicitly tagged rather than quietly accepted**

**Review.** *Reviewer H* (domain, per branch) on an independent sample per table; *Reviewer G* confirms every table now carries a tag and a unit declaration, and that C3.5 is honest.

**Exit gate.** Every `[ON-DISK]` table cited to page; plausibility suites green; C3.5 filed; zero silently-unverified tables. **Independent review filed (R2) and every §5 suggestion triaged (R4).**

**Effort: 60–105 h.** → **Sync point S2: must land before Track A Phase 6**, so the item pool is regenerated exactly once.

---

## Review agent protocol

The value of these gates depends entirely on the reviewers being genuinely independent.

### R0 — Mandatory phase-completion review

**No phase is complete until an independent reviewer has reviewed it and filed a report.** This is unconditional and applies to every phase in both tracks, including phases that produce no code (Phase 4, C1) and phases whose automated checks are all green. **A green test suite is not a completed phase.** Automated checks confirm the things we already knew to look for; the reviewer's job is to find the things we did not.

The review is **not** only a pass/fail gate. Every review must produce two kinds of output:

- **Backward-looking** — findings against what was built (defects, unmet criteria, weak evidence).
- **Forward-looking** — *suggestions for further probing, and improvements to carry into later phases.* A phase that passes every check but teaches nothing about where to look next has been under-reviewed.

A reviewer who reports no findings and no suggestions must state explicitly what they attempted and why it failed to surface anything. That is an acceptable outcome; silence is not.

### R1 — Rules

1. Reviewers receive: this spec, the phase diff, read access to the repo, and the phase deliverables. They **do not** receive the implementer's reasoning, notes, or self-assessment.
2. Reviewers **re-derive** every acceptance number independently. Confirming the implementer's number by reading it is not review.
3. Reviewers are instructed to **falsify**, not to confirm.
4. Every finding carries a **reproduction command**.
5. Findings are `CONFIRMED` (reproduced) or `PLAUSIBLE` (reasoned). **CONFIRMED blocks the merge.** `PLAUSIBLE` requires a written response, not necessarily a change.
6. A reviewer who cannot reproduce a claimed acceptance number **escalates rather than assuming they are wrong**.
7. Where a phase names two or more reviewers, they work **in parallel and in isolation from one another**. Convergent independent findings are a strong signal; a reviewer who has seen another's report cannot provide one.

### R2 — Review report structure

Filed as `docs/re-implementation-sep/reviews/phase<N>_<reviewer>.md`. Every report carries all five sections; an empty section must say why it is empty.

| § | Section | Content |
|---|---|---|
| 1 | **Verdict** | `PASS` / `PASS WITH FINDINGS` / `BLOCKED`, one line |
| 2 | **Independent re-derivation** | Each acceptance number, re-derived, with the command used. States agreement or divergence per number. |
| 3 | **Findings** | `CONFIRMED` / `PLAUSIBLE`, each with a reproduction command and an assessed impact |
| 4 | **Falsification attempts that failed** | What was tried to break the work and did not. **This is what makes a clean verdict meaningful** — without it, "no findings" is indistinguishable from "did not look". |
| 5 | **Further probing and improvements** | Forward-looking. See R3. |

### R3 — Section 5: further probing and improvements

The deliverable the rest of the protocol exists to enable. Reviewers are explicitly asked to answer:

- **Where is the evidence thinnest?** Which acceptance criterion is technically met but weakly supported, and what additional measurement would settle it?
- **What would you probe with more time?** Concrete next checks, ranked, with the reasoning for the ranking.
- **What generalises?** Does a defect or weakness found here likely exist in templates outside this phase's scope? Name them. *(This is how the corpus-wide back-solving pattern was found in the original audit — a defect seen in one template, hypothesised across others, then confirmed.)*
- **What should later phases do differently?** Changes to method, harness, or spec that would make the remaining phases more effective.
- **What in the spec was wrong or ambiguous?** Including this document. The spec is an artefact under review, not a fixed constraint.

### R4 — Disposition of suggestions

Suggestions are worthless if they are read and forgotten. Every item in §5 is triaged **before the phase closes**, and each receives one of:

| Disposition | Meaning |
|---|---|
| `ADOPT-NOW` | Actioned within this phase; the phase does not close until done |
| `ADOPT-PHASE-<N>` | Added as a named deliverable to a later phase — **and that phase's deliverable list is edited to include it** |
| `SPEC-CHANGE` | This document is amended; note the revision |
| `BACKLOG` | Recorded in the residual-risk register (D6.6) with reasoning |
| `REJECT` | Declined, with a written reason |

`BACKLOG` and `REJECT` require a one-line justification. **An untriaged suggestion blocks the phase gate**, exactly as a `CONFIRMED` finding does — a reviewer's forward-looking work is only as good as the disposition it receives.

### R6 — Scoping a review

**This is the rule most likely to be violated, and violating it costs both time and rigour at once.**

The instinct when commissioning a review is to ask for everything, on the theory that more asked for
is more assurance obtained. It is not. A review's value comes from **covering distinct properties**,
never from **repeating the same measurement**. An over-scoped review stalls, produces nothing, and
the gate it was guarding silently degrades into a formality — so bloat does not trade time against
robustness, it spends time *and* loses robustness.

There is also a second-order cost. Review is ~30% of Track A effort and is the first thing anyone
proposes cutting when a schedule slips. Every wasted reviewer-hour is ammunition for that argument.
**Keeping reviews lean is how the review protocol survives contact with a deadline.**

#### R6.1 — One gate property per review, mandatory; everything else is time-boxed

Every review has exactly one **mandatory** task: the single property the phase's exit gate depends
on. For Phase 0 that is *"the harness detects a planted defect"* — nothing else certifies the gate.
Forward-looking probing (R3) stays valuable but is explicitly time-boxed, and the brief must tell the
reviewer to **file what it has when the box expires** rather than continuing. State the box in the
brief as a number, not as "be efficient".

#### R6.2 — Every measurement has exactly one owner

Two reviewers must never be assigned the same number. Before dispatching a batch, write a
**measurement ownership register** — one line per number, one owner each — and check that nothing
appears twice:

| Measurement | Owner | Not re-measured by |
|---|---|---|
| audit defect rates (all 9 claims) | D0.5 agent | adversary, oracle agents |
| planted-defect detection | adversary | — |
| per-template round-trip failure rate | the oracle for that template | D0.5, adversary |

Convergence is only evidence when it comes from **different questions**. The two independent readings
of `template_damping_classification` — one agent applying a strict `zeta == 1` test, another
propagating the stated precision of `c` into a band — were worth having precisely because the briefs
differed; their disagreement located a real ambiguity in the item. The same brief handed to two
agents produces agreement that proves nothing, and costs twice.

#### R6.3 — Supply the tooling; a reviewer must not build its own instrument

If a review needs a driver, a scratch harness, a fixture, or a non-obvious API, hand it over **as
working code in the brief**. Time spent reverse-engineering the instrument is time not spent
reviewing, and it is the most common way a review stalls. If you cannot write that snippet quickly
yourself, the reviewer certainly cannot, and the task is mis-scoped.

#### R6.4 — Name what is already settled, and forbid re-deriving it

Briefs must state explicitly what is **out of scope because it is already known**, and cite where the
answer lives. "Do NOT re-measure defect rates — they are in `phase0_baseline.md`" is a required line,
not a courtesy. Reviewers default to thoroughness and will re-derive anything not explicitly fenced
off.

#### R6.5 — If a review cannot land inside its box, the scope was wrong

Re-cut to R6.1 and re-dispatch. **Do not extend the box.** A review that needs more time than its box
is telling you the brief bundled several reviews together.

#### Anti-patterns, named so they are recognisable at dispatch time

| Anti-pattern | What it looks like | Cost |
|---|---|---|
| **Duplicate commissioning** | two briefs in one batch ask for the same number | 2× cost, 0 extra assurance; and it starves the gate task |
| **Kitchen-sink brief** | "also analyse all seven checks for gaps" | unbounded; the mandatory task competes with it and loses |
| **Instrument-building** | reviewer must construct fixtures before it can start | stall risk; often the whole budget |
| **Re-deriving settled numbers** | not fencing off what is already measured | silent duplication the author never sees |
| **Unbounded probing** | R3 questions with no time box | a review that never files |

#### Pre-dispatch checklist

Run this before launching any review or batch of reviews:

- [ ] Exactly **one** mandatory task, and it is the phase's gate property
- [ ] Measurement ownership register written; **no number appears twice**
- [ ] Everything already settled is named, with a pointer to where it lives
- [ ] Required tooling supplied as working code in the brief
- [ ] A time box stated as a number, with "file what you have" instructions
- [ ] Optional work explicitly marked optional and separately boxed

> **Why this rule exists.** The first Phase 0 Harness Adversary brief carried three tasks: plant
> defects and confirm detection (the gate), independently re-derive three audit claims, and analyse
> all seven checks for gaps. The second task **duplicated the D0.5 agent's brief**, dispatched in the
> same batch by the same author, which was re-measuring all nine claims — the same numbers
> commissioned twice. The review ran 45 minutes past its last write and produced nothing; D0.5
> delivered the superset in 23. **The only mandatory task went unanswered because it was bundled with
> work already being done elsewhere.** Re-dispatched at one task with the plumbing supplied and a
> 15-minute box, it was tractable.

### R5 — Standing deliverable and gate item

Appended implicitly to every phase in both tracks:

- **D\<N\>.R** — independent review report(s) per R2, plus the R4 triage table.
- Exit-gate item: **☐ Independent review filed (R2) and every §5 suggestion triaged (R4)**.

For Phase 6, D6.R additionally carries a **retrospective**: which suggestions across all phases were adopted, which deferred, and what the review process itself should do differently next time.

### Roster

| ID | Role | Phases |
|---|---|---|
| A | Correctness / determinism | 1, 2, 5 |
| B | Physics & pedagogy (P6 guard) | 1, 3, 4 |
| C | Numerical methods | 2 |
| D | Schema implementability | 3 |
| E | Comparator adversary | 4 |
| F | Audit replication | 6 |
| G | Provenance | C1, C3 |
| H | Domain (thermochemistry / materials / fluids) | C2, C3 |
| — | Harness adversary | 0 |

**Anti-pattern to watch for, phase 1 especially:** a T2 round-trip oracle that reuses the template's own expression proves nothing — it will agree with any arithmetic, correct or not. Reviewer A's primary task is to confirm each oracle is an independent re-derivation from the governing equation, not a copy.

---

## Sequencing and effort

```
TRACK A (templates)                          TRACK B (constants)

Phase 0 ──────────────────►                  C1 Sources & convention
   │  (blocks all of Track A)                    │   (starts immediately;
   │                                             │    independent of Phase 0)
   ├──► Phase 1  round-trip ───────►            │
   │                                             ├──► C2 chem thermochem
   ├──► Phase 5  contract ─────────►            │        (critical path)
   │                                             │             │
   ├──► Phase 4  non-numeric ──────►            │             │
   │        │                                    │        ┌────┘
   │        ▼                                    │        │  S1
   ├──► Phase 3  trace shape ──────►            │        ▼
   │                                             │   Phase 2  determinism
   │                                             │
   │                                             └──► C3 remaining tables
   │                                                          │
   │                                                     S2   │
   └──────────────────────────────────────────────────────────┤
                                                              ▼
                                                          Phase 6
                                                   (single item-pool
                                                    regeneration)
```

**S1 — C2 before Phase 2.** `template_adiabatic_flame_temperature` must not be re-baselined against `CP_PARAMS` values that are 10× wrong.
**S2 — C3 before Phase 6.** So the item pool is regenerated exactly once, covering both tracks.

| Phase | Implementation | Review | Total |
|---|---:|---:|---:|
| **Track A** | | | |
| 0 — Infrastructure | 25–40 h | 6–8 h | **31–48 h** |
| 1 — Round-trip | 15–25 h | 8–12 h | **23–37 h** |
| 2 — Determinism | 12–24 h | 8–12 h | **20–36 h** |
| 3 — Trace shape | 8–16 h | 10–14 h | **18–30 h** |
| 4 — Non-numeric | 12–20 h | 8–10 h | **20–30 h** |
| 5 — Contract hygiene | 6–10 h | 3–4 h | **9–14 h** |
| 6 — Consolidation | 10–16 h | 8–10 h | **18–26 h** |
| *Track A subtotal* | | | **139–221 h** |
| **Track B** | | | |
| C1 — Sources & convention | 10–16 h | 2–4 h | **12–20 h** |
| C2 — Chem thermochemistry *(S1)* | 12–20 h | 3–5 h | **15–25 h** |
| C3 — Remaining tables *(S2)* | 50–90 h | 10–15 h | **60–105 h** |
| *Track B subtotal* | | | **87–150 h** |
| | | | **226–371 h** |

Review is ~30% of Track A and ~15% of Track B (transcription verifies more cheaply than design does). It should not be trimmed — **Phase 0**, the **Phase 3 schema review**, and **C2's independent re-derivation** are the three places where cutting corners silently invalidates everything downstream.

**Critical path** is Track A: `Phase 0 → 1 → 4 → 3 → 6`. Track B's C2 is short enough (15–25 h) that if C1 starts on day one it will never be what holds up Phase 2. C3 is the larger risk to the end date — it is the one place worth adding a second person, since table transcription parallelises cleanly across branches and design work does not.

---

## Execution model — the table above is EFFORT, not wall-clock

Read the hours as effort. Wall-clock is a fraction of them, and the difference is entirely in how the work is scheduled. Four rules.

**E1 — Fan out within a phase, not just across phases.** Phase 1's 9 templates are independent of each other; so are C3's tables, and Phase 4's four comparators. Work them concurrently rather than in sequence. A phase's wall-clock should approach its *slowest single item*, not the sum of its items.

**E2 — Review gates the merge, not the start of the next phase.** Reviewer A can be reviewing Phase 1 while Phase 4 implementation is already under way. Since review is ~30% of Track A effort, serialising it would add ~30% to the schedule for no gain in rigour. The only hard rule stays: **nothing merges before its own gate passes.**

**E3 — Track A and Track B run concurrently from day one.** C1 has no dependency on Phase 0. Only S1 and S2 bind the tracks.

**E4 — Independence is a reason to parallelise, not a cost of it.** The T2 oracles are the clearest case: an oracle written by the same person who just read the template's computation will tend to reproduce it, which is the single most likely way this plan fails silently (see Risks). Writing each oracle in a separate context, from the governing equations with the computation code withheld, is both faster and *more* rigorous than writing them serially.

**E5 — Speed comes from not repeating work, never from checking less.** The three legitimate sources of speed in this plan are: doing independent things concurrently (E1–E3), refusing to commission the same measurement twice (**R6.2**), and building tooling once so verification is seconds rather than hours (Phase 0 — the full suite runs the corpus in ~17s). Everything else that looks like speed is scope reduction wearing a disguise. Before cutting anything to save time, check it is not actually duplicated work that R6 would have removed for free.

**What must not be compressed:** Phase 0 itself (it is what makes every later verification cheap), the independent reviews (run them concurrently and keep them lean per R6, but do not skip them), and the T2 oracles' independence (parallelise the writing; do not thin the oracles).

### Revision log

| # | Change | Reason |
|---|---|---|
| SPEC-CHANGE 1 | T5 extended from "value binding" to "value binding **and rounding discipline**" (T5a/T5b) | T1 cannot see the round-then-recompute (P2) defect class: in `template_cantilever_double_integration` every printed line closes within its display tolerance while the trace is still wrong. T5b catches it, and caught `template_impulse_response_from_lccde` too. |
| SPEC-CHANGE 2 | T3 batched into two child processes for the whole corpus rather than two per template | 300 interpreter starts → 2. ~500s → ~8s on 150 templates, identical comparison. |
| SPEC-CHANGE 3 | Execution model E1–E4 added above | The phase table was being read as a serial schedule. |
| SPEC-CHANGE 5 | **P2 amended: break the rounding tie in decimal** (D-012) | P2's own worked exemplar applies it correctly and still fails closure on 5.89% of instances. Applied as written, Phase 1 would have reproduced the residue eight more times. |
| SPEC-CHANGE 6 | **Phase 1 re-scoped 9 → 12 templates, split into chain breaks vs display defects** (D-011, D-013) | Phase 0 measurement: three templates are not defective, three more have the same defect shape, and the two categories need different fixes. |
| SPEC-CHANGE 7 | Worst-case defect rates to be quoted as distributions, not sample maxima | "7.77%" and "0.87%" are extreme values of a 200-seed sample; at 20,000 seeds they are 14.27% and 1.21%. |
| SPEC-CHANGE 4 | **R6 review-scoping rules added** | The first Phase 0 adversary brief bundled the gate task with work already assigned to the D0.5 agent. It stalled and produced nothing; D0.5 delivered the superset. Over-scoping a review does not make it more thorough — it makes it not happen. |
| SPEC-CHANGE 8 | **§3.2's claim that `iteration` and `decision` "are the same shape" is amended** (D-038) | The structure is shared; the equivalence is not. They differ on whether the sequence's cardinality is an observable of the answer, so a comparator built from one merged type is wrong on one of the two templates. Delivered as two types over one base in [`phase3_node_types.md`](phase3_node_types.md). |
| SPEC-CHANGE 9 | **A screen's *rejected slice profile* is now a required recorded measurement** in every item-pool-impact note — the marginal of every sampled parameter over the rejected set (Phase 3 Reviewer B, §5) | Both Phase 3 screens reject **clustered**, not scattered, instances — one at a single slope value, removing 25% of it. Both were described as removing "ill-posed instances"; both descriptions were true and incomplete. The profile would have surfaced it with no reviewer, and my own "scattered" claim was measured over a *union* of screens, which cannot see a single-valued component (D-045). |
| SPEC-CHANGE 10 | **A pedagogy lookup-check must fit a model, not enumerate rules**, and must report *lift over a blind-guess floor* against a stated threshold (Phase 3 Reviewer B, §5) | B enumerated every shortcut it could construct and reached 64.97% against a 53.05% floor, concluding the item was a search. A depth-2 decision tree reaches 90.18% and reduces to one line. **"The best rule I could think of" is a floor on shortcuttability, never a ceiling** (D-046). A bare rate is also unreadable without its floor: 58.6% sounds alarming until the floor is 53.1%. |
| SPEC-CHANGE 11 | **§4.2's hedge requirement is demoted from enforced to reported** (D-056) | Four review rounds and an ablation. The hedge layer fires on **2 of 2,200** archived answer spans and **0** in a kind it gates; ablating it changes no positive recall and no archive verdict, and it produced 13 of Reviewer E's 20 findings while repeatedly marking *correct* answers wrong. A hedge is now detected, annotated and counted, never scored. `ENGTRACE_HEDGE_POLICY=enforce` restores the old behaviour. Same treatment `narrative` gets (D-051). |
| SPEC-CHANGE 12 | **`multipart` is not a seventh comparator `kind`** (D-047) | It is a property of the *answer*: an ordered list of `parts`, each carrying one of the six kinds, with a `mode` of `all` or `any`. `system_properties_memory_causality` is the proof inside Phase 4 — typed `classification`, described by §4.1 as a "label tuple", actually two categorical parts. A seventh kind would stop the six composing. Governs **24 of the 51 class-C templates**. |
| SPEC-CHANGE 13 | **D3.4 §7 amended in four places** (D-054) | §7 shipped unexercised and both Phase 3 schema reviewers named the missing corpus as the largest gap. Building it (D4.6) found four defects in the *specification*: §7.2's headline disposition has no instances, §7 never says whether §6 step 8 binds a candidate, §7.4's direct solver cannot be an `iteration` node at all, and §7.3 row 3's set-versus-ordered distinction is vacuous. |

---

## Risks

| Risk | Impact | Mitigation |
|---|---|---|
| **T2 oracles reuse template arithmetic** | Phase 1 gate passes while defects remain — the single most likely failure of this plan | Reviewer A's explicit primary task |
| **Audit defect rates don't reproduce** | Phase assignments wrong; work misdirected | Phase 0 D0.5 re-measures before any edit |
| Phase 1 changes item difficulty | Benchmark statistics shift; paper claims affected | T6 gate + Reviewer B; P6 sign-off |
| Removing fallbacks changes the sample distribution | Items no longer comparable to published runs | Expected in Phase 2; D2.5 documents and approves |
| ~~`damping_classification` option 1 infeasible~~ | — | **Retired.** Feasibility verified: 680,869 valid (k, m) pairs in the current ranges (§1.3) |
| Raw `inference_results/` generations were not retained | Parser fix needs full re-inference (11 models × 1,350 items) instead of a free re-score | **Check before promising a corrected table** — they are gitignored and absent from this working copy (decision 5) |
| Phase 3 schema under-specified | Class C work later inherits an ambiguous model | Reviewer D must *implement* from the spec alone |
| Regenerated items invalidate published results | Paper revision uses stale numbers | D1.6 + D6.4, circulated before revision |
| Phase 2 starts before `CP_PARAMS` fixed | Determinism verified against wrong physics | Hard dependency in the sequencing graph |

---

## Open decisions — recommendations

These are benchmark-ownership calls, not implementation details. **Decisions 1–3 should be settled before their phases start**, not during. A recommendation is given for each, with the condition that would change it.

**1. `damping_classification` → constrain the sample, and invert the construction.** (§1.3)
Sample `(k, m)` from the exact-`c_c` family and sample `c` **directly** as a 2-dp value; classify by exact decimal comparison `c` vs `c_c`; derive ζ for display only. Feasibility verified: 680,869 pairs in the current ranges. All three classes preserved, the float-equality test disappears, no pedagogical change.
*Changes if:* Reviewer B judges the constrained `(k, m)` values implausible to a domain expert — measured samples suggest they are not.

**2. `normal_depth_iteration` → schema `iteration` node.** (§3.1)
Fixing the iteration count would convert "iterate until converged" into "perform three updates," which is a materially easier item; the convergence count is genuinely data-dependent (2–6 evaluations measured).
*Changes if:* the milestone model is frozen before Phase 3 and cannot accept a new node type — then fix the count and record the P6 trade.

**3. `line_balancing_heuristic` → schema `decision` node.** (§3.2)
Same objection as #2: the station count is a computed result and fixing it leaks part of the answer. Shares a solution family with #2, so the marginal cost over #2 alone is small.
*Changes if:* #2 is decided the other way — the two should not diverge.

**4. Answer marker → normalise the 5 chemical templates to `**Answer:**`.** (§5)
One canonical gold marker. A widened accepted set has to be known by every downstream consumer and will drift as templates are added; the 5 files are being touched in Phase 5 regardless.
*Important asymmetry:* strictness applies to **gold traces only**. The value extractor reading **model** output must stay permissive — models emit "Final answer:", "Therefore," headings, or nothing at all. Do not let this decision leak into the extractor's tolerance.

**5. Item-pool regeneration → split by cost, and check a dependency first.**
The two fixes have very different costs and should not be bundled:
- **Parser fix → re-score only, immediately.** It changes scoring of *existing* generations; no items change and no inference is re-run. **Blocking dependency: do the raw `inference_results/*.jsonl` still exist?** They are gitignored and **absent from this working copy**. If retained elsewhere, this is a near-free correction of a known error in published numbers and should happen this week. If lost, it needs full re-inference across 11 models × 1,350 items, which is a budget decision — **establish which, before promising a corrected table.**
- **Template changes → regenerate once, after Phase 6.** Per-phase regeneration would invalidate comparisons repeatedly; one regeneration gives a single clean before/after.

**6. Phase 5 scoping → separate parallel PR.**
Keep diffs semantically homogeneous: mixing cosmetic marker fixes into Phase 1's semantic rounding changes makes Reviewer A's job harder and weakens the strongest gate in the plan. Phase 5 has its own trivially clear gate (T4 green on 150/150) and its own reviewer.
*Exception:* `levenspiel_plot_interpretation` appears in both — sequence it after Phase 2.
