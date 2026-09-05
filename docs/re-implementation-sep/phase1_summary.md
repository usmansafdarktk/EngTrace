# Phase 1 — Summary and close-out

**Status: COMPLETE WITH FINDINGS. Exit gate met; see §6 for what is open.**
**Date:** 2026-09-06 · **Branch:** `redesign/phase1-round-trip`

Phase 1 fixed the templates whose gold reasoning traces did not reproduce their
own answers. These were **defects in published results**, not design
limitations: a solver who read the question, followed the printed steps and used
the printed operands did not arrive at the printed answer.

---

## 1. What was actually in scope: **ten** templates, not twelve

The phase is called "12 templates" in its title, in D-011 and in the
implementation brief. **The correct number is ten**, and ten is what was edited.

D-011 grew the phase "from 9 to 12" by adding three templates, without
subtracting the three that left in the same table — D-013 found `poissons_ratio`,
`logarithmic_decrement` and `system_properties` not defective. 9 + 3 − 3 = 10,
which is exactly what the authoritative scope table (spec §1.1a) enumerates.
Recorded as **D-023**; found by Reviewer A. A third number, "all 9", is still in
the exit gate at §1.6 and is stale.

| Category | Templates | Fix |
|---|---|---|
| **Chain break** (3) | `mean_variance`, `rotating_unbalance`, `vibration_transmissibility` | P3 — state, then recompute forward |
| **Display defect** (6) | `beam_deflection_formula` (first), `cantilever_double_integration`, `annulus_flowrate`, `statically_indeterminate_shaft`, `shaft_design_power`, `composite_shafts_series` | P2 as amended → **D-016** |
| **Ill-posed, not wrong** (1) | `damping_classification` | state `c` as the critical value; drop the float-equality test |

---

## 2. Defect rate, before → after

Measured at **20,000 seeds per template**, both checks, after the reviews. The
"before" column is the like-for-like measurement on `master`, not the figure in
`phase0_summary.md` — that one was taken at the runner's default 25 seeds.

| Template | T1 closure failures | T2 round-trip failures | worst rel. error |
|---|---|---|---|
| `mean_variance` | — → **0** | 56.8% → **0** | 4.09e-4 |
| `rotating_unbalance` | 597/200 seeds → **0** | 4.6% → **0** | 4.63e-4 |
| `vibration_transmissibility` | 212/200 seeds → **0** | 4.5% → **0** | 1.81e-3 |
| `damping_classification` | 0 → **0** | 0 → **0** | 2.53e-5 |
| `beam_deflection_formula` | 5.89% of SI → **0** | (new oracle) **0** | **0.00e+00** |
| `cantilever_double_integration` | 5.20% → **0** | 0 → **0** | **0.00e+00** |
| `annulus_flowrate` | 94% of Step-5 → **0** | 0 → **0** | 5.33e-5 |
| `statically_indeterminate_shaft` | 47% of lines → **0** | (new oracle) **0** | 1.35e-3 |
| `shaft_design_power` | 30% of lines → **0** | (new oracle) **0** | **0.00e+00** |
| `composite_shafts_series` | 89.9% verbatim → **0** | (new oracle) **0** | 9.93e-5 |

**Target was 0. Target met**, at a seed count 20× the one the written gate names
(see D-024 for why that mattered).

Corpus-wide, nothing outside scope regressed. At a matched 200 seeds on both
sides: T1 56 → 48 failing templates, T2 3 → 0, T5 73 → 71, T7 95 → 87, T4
unchanged at 4. **Zero newly-failing templates anywhere.**

---

## 3. The three findings that justified the phase

**P2 as amended was still wrong, and applying it verbatim was worse than the
defect.** D-012 said to break the rounding tie in decimal. Measured, that turns a
T1 MARGINAL into a T1 **FAIL** — 359 failures per 5,000 seeds on
`beam_deflection_formula`, confirmed independently by Reviewer A. At an exact
display tie `|evaluated − printed|` equals the tolerance exactly, and float
representation decides FAIL versus MARGINAL **in either rounding direction**. No
rounding convention can pass. A tie is not a rounding problem; it is an
**ill-posed instance**, and the fix is to remove it rather than resolve it
(**D-016**). This is the second time the recipe for this phase was wrong and was
caught only by measuring it.

**Rounding to display precision is not enough; the value must be bound
*through* its display.** `34.4 * 1e-6` is `3.4399999999999996e-05`, prints as
`3.4400e-05`, and reparses to `3.44e-05` — a different float. One ulp puts the
template and the reader on opposite sides of a tie. Reviewer A found four such
seeds in 55,000, **concentrated on particular (section, span) pairs rather than
spread uniformly**, so they would have recurred systematically in a regenerated
pool. Fixed by `_as_printed(x, spec)` and by quarantining a tie *band* rather
than a point.

**Two "dead" code paths were live, and one of them was rewriting the physics.**
`annulus_flowrate`'s `if pressure_drop_Pa == 0: pressure_drop_Pa = 10.0` fired
on **17% of instances**: the sampler's Reynolds targeting was discarded and the
stated pressure drop bore no relation to the flow described. Reviewer B
independently established that **17% of instances were also non-laminar**
(Re up to 2149) while asserting the laminar solution. Both were found by writing
the invariant asserts T7 requires — not by any check that existed before.

> **Carry this forward:** every one of these was found by measuring something
> the phase was not asked to measure. The gate as written would have passed all
> three.

---

## 4. Deliverables

| # | Deliverable | Where | Status |
|---|---|---|---|
| D1.1 | Templates edited, satisfying P1–P3 | 6 files under `data/templates/branches/` | ✅ **10**, per §1.1a (D-023) |
| D1.2 | T2 oracles for the newly-added templates | `tests/template_integrity/oracles/` | ✅ **4 written** — the 3 required, plus `composite_shafts_series`, which had none |
| D1.3 | Per-template `TOLERANCE` + justification | `TEMPLATE_TOLERANCES` in each edited file | ✅ + measured `AGREEMENT_FLOOR` / `DETECTION_FLOOR` in the new oracles |
| D1.4 | `damping_classification` decision record | D-010, **D-020** | ✅ |
| D1.5 | T6 before/after distribution diff | §5 below; `tests/template_integrity/t6_report.py` | ✅ at 5,000 seeds |
| D1.6 | Item-pool impact note | [`phase1_item_pool_impact.md`](phase1_item_pool_impact.md) | ✅ |
| D1.7 | Phase report | this file | ✅ |
| D1.R | Two independent reviews + R4 triage | §7, §8 | ✅ |

**Every oracle written this phase is mutation-tested**: injecting a 1% error
into the gold answer fails 200/200 sampled instances, and the unmutated text
passes 200/200. Reviewer A additionally planted **14 structural defects** —
branch swap, dropped unit conversion, radius/diameter confusion, dropped π/2,
inverted length ratio, a whole segment omitted — and every one turned its oracle
red on 98–100% of reachable instances.

---

## 5. T6 distribution diff, before → after (D1.5)

**Both sides at 5,000 seeds.** The committed 1,000-seed baseline is not usable
for this: a distinct-answer count cannot exceed the seed count, and at N=1,000 a
template with ~3,200 reachable answers reads ~900 whatever its true answer space
is. That saturation hid a real 29% regression (**D-021**).

| Template | distinct answers | median | steps | verdict |
|---|---:|---:|---|---|
| `mean_variance` | 2,975 → 4,132 (**+38.9%**) | +0.63% | same | pass |
| `rotating_unbalance` | 3,242 → 4,099 (**+26.4%**) | +0.00% | same | pass |
| `vibration_transmissibility` | 1,544 → 3,614 (**+134%**) | −0.62% | same | pass |
| `damping_classification` | 1,600 → 1,600 (+0.0%) | +0.00% | same | pass |
| `beam_deflection_formula` | 1,828 → 1,827 (−0.1%) | +0.00% | same | pass |
| `cantilever_double_integration` | 395 → 395 (+0.0%) | +0.00% | same | pass |
| `annulus_flowrate` | 4,929 → 4,940 (+0.2%) | **−28.2%** | same | **BREACH — signed off, D-022** |
| `statically_indeterminate_shaft` | 4,886 → 4,684 (−4.1%) | −0.02% | same | pass |
| `shaft_design_power` | 3,243 → 3,234 (−0.3%) | −0.04% | same | pass |
| `composite_shafts_series` | 2,730 → 4,899 (**+79.5%**) | +0.04% | same | pass |

**Step-count distributions are identical before and after for all ten**
(verified independently by Reviewer B). The large distinct-answer gains are
where an answer had been quoted to too few significant figures — a
transmissibility of 0.056 printed to 3 dp is two significant figures, and
ungradeable.

**One breach, explicitly signed off, tolerance not widened:** `annulus_flowrate`'s
median moves −28% because the fix removes an artificial floor that had pinned
17% of instances at ΔP = 10 Pa. The shift is entirely a downward extension of
the low-Q tail. Reviewer B, the P6 guard, signed it off independently (D-022).

Two apparent breaches at 1,000 seeds were **noise, not signal**, and are recorded
rather than silently dropped: `vibration_transmissibility` read −2.3% at N=1,000
and −0.62% at N=5,000, and the *unmodified* template's own median differs by
3.14% between disjoint halves of the same 1,000-seed range (**D-019**).

---

## 6. Exit gate

| Criterion | Status |
|---|---|
| T1 closure: no non-closing printed line, all in scope | ✅ **0 at 20,000 seeds each** (the gate says 1,000; D-024 restates it) |
| T2 round-trip passes for every template with an oracle | ✅ 0 failures at 20,000 seeds; 10 of 10 now have oracles |
| T5b: zero **confirmed** rounding violations | ⚠️ **see below** |
| T3 determinism | ✅ |
| T6 within tolerance, or breach signed off | ✅ one breach, signed off (D-022) |
| T7 ≥2 invariant asserts each | ✅ 3–4 output asserts each |
| Reviewer A and Reviewer B both clear | ✅ after remediation — see §7 |
| Every §5 suggestion triaged | ✅ §8 |
| Item-pool impact note filed | ✅ |

**T5b — the one criterion not cleanly green.** T5 reports one finding each
against `Ix` (`beam_deflection_formula`, `cantilever_double_integration`),
`c_cubed` (`shaft_design_power`), `j1`/`j2` (`composite_shafts_series`) and
`term1`/`q_flow_rate` (`annulus_flowrate`). **None is a real violation.** Each
value is bound through its own display with `_as_printed`, and is verified equal
to it on 4,000 instances. T5b flags them because its precision regex reads
`.3e`/`.4e`/`.5e` as "N decimal places" when the digits are **mantissa**
decimals — `{Ix:.4e}` on `8.49e-05` displays five significant figures and loses
nothing, while T5b computes `round(8.49e-05, 4) = 0.0001` and reports a loss of
5e-05.

This is a **harness defect (D-017), raised and deliberately not fixed**: the
spec forbids modifying the harness to make a fix pass, and Phase 0 found two
checks that were green while measuring nothing. Fixing it is a `SPEC-CHANGE` for
the harness owner. The gate criterion says "zero *confirmed* rounding
violations", and zero are confirmed.

---

## 7. Review outcome

Two reviewers, in parallel, in isolation from each other, neither given the
implementer's reasoning.

**Reviewer A — Correctness: PASS WITH FINDINGS**, two CONFIRMED blockers, both
now fixed.

| Finding | Disposition |
|---|---|
| **F-1** `annulus_flowrate` still failed T1 — 6 exact-display-ties per 20,000 seeds on the `kappa` line; the tie guard had been applied to two templates and not this one | **Fixed.** Tie guard extended to `kappa` and the shape factor; radii bound through their displays. 0 failures in 20,000. |
| **F-2** `shaft_design_power` TOLERANCE 1e-3 sat below its own quantisation floor; seed 11230 gave 9.13 mm through the stated intermediates and 9.12 mm from the givens | **Fixed, by removing the ambiguity rather than tolerating it** — the template now resamples when the givens-path and intermediates-path answers differ. Worst relative error 0.00e+00 at 20,000 seeds. Widening the tolerance instead would have raised the smallest detectable error from 0.2% to 0.5%. |
| **F-3** the `beam_deflection_formula` oracle docstring described a 5-dp intermediate the merged code no longer has | **Fixed.** Docstring corrected. |
| **F-4** every oracle's declared TOLERANCE is argued, and in several cases is not the number that binds — display quantisation sets the real sensitivity | **Adopted.** `AGREEMENT_FLOOR` and `DETECTION_FLOOR` are now **measured** and declared in the four new oracles, with `tests/template_integrity/oracle_floors.py` to reproduce them. Two template-level declarations were **below their agreement floor** and were corrected. |
| **F-5** oracle sensitivity is coupled to constants tables Track B is about to change | **Recorded as D-025**, assigned to Track B's sync points. |

**Reviewer B — Physics & pedagogy: BLOCKED**, three CONFIRMED findings, all now
addressed.

| Finding | Disposition |
|---|---|
| **F-1** `beam_deflection_formula`'s US Step 2, titled "Assemble consistent US customary units", no longer converted anything — the kip/ft→kip/in step had become a no-op, and the SI branch still converted numerically, so the same item tested different things by unit system | **Fixed.** Step 2 performs the conversion and states its value to 6 dp; the substitution still carries the exact ratio. |
| **F-2** the `annulus_flowrate` T6 breach, `rotating_unbalance`'s question change and F-1 were P6 scoping decisions with no written record | **Fixed — D-022**, with Reviewer B's sign-off recorded. |
| **F-3** `mean_variance` stated its variance before showing the working, then repeated the total | **Fixed.** Deviations moved to their own line, main line runs deviations → products → total once. The obvious repair — one chained line — was measured and **fails 1,990 of 2,000 instances**, because a segment beginning `(-5.890)^2*(…)` reads as a terminal value of −5.890. |

Reviewer B's verdict on the P6 question, per template: **nine of ten "no change
or better"**; the tenth was F-1, now fixed. Reviewer B independently confirmed
that the six edited files contain 24 templates and **exactly 10 changed output**
— no collateral damage — and that governing equations are unchanged in all ten.

Reviewer B also **accepted D-020** on `damping_classification` after reproducing
the separation independently, and recommended D-004 stay closed.

---

## 8. R4 triage — every §5 suggestion

| # | Suggestion | Disposition | Action |
|---|---|---|---|
| A-1 | Verify closure at the seed count D-016 mandates, not the gate's 1,000 | `ADOPT-NOW` | Done: all ten at 20,000 seeds. **D-024** restates the gate. |
| A-2 | Put a **floor** on T1 coverage; the spec's own 60% rule makes most of these findings today | `SPEC-CHANGE` | Measured: coverage is 0.13–0.50, but **0 unparseable lines in 39,314 `=` lines** — every skipped line is legitimately symbolic (givens restatements, `phi_AC = phi_CB`). The 60% rule conflates *unparseable* with *symbolic*. T1 should gate on `unparseable == 0`, and report coverage as `evaluated / (evaluated + unparseable)`. |
| A-3 | Independently re-derive the T1 marginal-band density that D-016's tie argument rests on | `ADOPT-PHASE-2` | Named deliverable for Phase 2's review. |
| A-4 | Declare a **measured** detection floor per oracle | `ADOPT-NOW` | Done for the four new oracles; `oracle_floors.py` committed. The nine pre-existing oracles are `ADOPT-PHASE-2`. |
| A-5 | Gate the constants track against templates' stated precisions | `ADOPT-PHASE-C1` | **D-025**. |
| A-6 | Reconcile the 12 / 10 / 9 template counts | `SPEC-CHANGE` | **D-023**. |
| A-7 | Put T2's display-quantisation practice, and its consequence, in the spec | `SPEC-CHANGE` | §0.2 T2 to state that both sides quantise and that this, not TOLERANCE, sets sensitivity. |
| A-8 | Don't cite `poissons_ratio` / `logarithmic_decrement` / `system_properties` as *verified* — D-013 judged them before F-4 was understood | `BACKLOG` | Residual-risk register, D6.6. Not reopened. |
| B-1 | Systematic physical-plausibility screening as a harness check, not a reviewer's eye | `ADOPT-PHASE-6` | Would have caught the 11.7 GPa annulus pressure drop and the 185° composite twist, **both of which predate this phase**. |
| B-2 | Check whether raised answer precision degrades grading tolerance | `ADOPT-PHASE-2` | A verifier-tolerance question; T6 structurally cannot see it. |
| B-3 | `damping_classification`'s 1/3 critically-damped population is implausible for a classification task | `BACKLOG` | Unchanged from `master`; not this phase's business. |
| B-4 | Make the P6 record a deliverable, not a footnote | `ADOPT-NOW` | **D-022**, and carried into the Phase 2 brief. |
| B-5 | Give T6 a **question-text hash** — it profiles only the solution, and P6 is about the question | `SPEC-CHANGE` | **D-022**. `rotating_unbalance`'s wording change was invisible to every gate in this phase. |
| B-6 | Brief the P6 review against **before/after instance pairs at matched seeds**, not the source diff | `ADOPT-NOW` | Adopted for Phase 2's brief. Five of B's findings came from the instance diff in ten minutes; none is legible in the 4,000-line source diff. |
| B-7 | Name the sign-off *artefact* and *approver* for a T6 breach; P6 gives no threshold for "alters what an item tests" | `SPEC-CHANGE` | Proposed: the artefact is a DECISIONS entry, the approver is the phase's P6 reviewer. |
| B-8 | ΔP stated to 11 significant figures on viscous fluids | `BACKLOG` | Pre-existing in kind; folded into B-1. |
| B-9 | `Pa.s` reads oddly; corpus mixes `N-s/m` and `N.s/m` | `ADOPT-PHASE-5` | Unit-symbol consistency is Phase 5's contract work. |
| B-10 | The annulus distinct-answer count was passed on a saturated statistic | `ADOPT-NOW` | Re-measured at 5,000 seeds (4,929 → 4,940). Generalised as **D-021**. |

**No suggestion is untriaged.**

---

## 9. Decisions recorded this phase

**D-016** P2 amended again — remove the tie, do not resolve it ·
**D-017** T5b misreads `%e` as decimal places (harness, not fixed) ·
**D-018** two templates had zero T1 coverage ·
**D-019** T6's median gate is inside its own noise at 1,000 seeds ·
**D-020** `damping_classification`'s strict flip rate cannot change without D-004 ·
**D-021** T6's distinct-answer gate saturates and hid a 29% regression ·
**D-022** three P6 scoping decisions, with Reviewer B's sign-off ·
**D-023** the phase edits ten templates, not twelve ·
**D-024** the exit gate's seed count is superseded ·
**D-025** oracle sensitivity is coupled to Track B's constants tables.

---

## 10. What Phase 2 inherits

1. **The fix pattern, in three parts** (D-016): one rounding at the answer's
   precision with upstream displays matched to it; every displayed-and-consumed
   operand bound *through* its display; and residual ties resampled, as a band
   not a point. Parts 2 and 3 exist only because reviewers found parts 1 and 2
   insufficient — assume the same of this pattern.
2. **Size the acceptance run to the defect rate before running it.** This phase
   wrote that lesson down in D-016 and then shipped a template verified at 1,000
   seeds that failed at 20,000. Phase 2's gate should name 20,000.
3. **Four harness `SPEC-CHANGE`s are open and none is fixed**: T5b's `%e`
   misparse (D-017), T1 coverage gating on `unparseable` rather than a
   percentage (A-2), T6's saturating distinct-answer count (D-021), and T6's
   blindness to the question text (D-022/B-5). Phase 2 inherits a harness that
   is known-wrong in four specific ways.
4. **Two invariant asserts are worth more than they look.** T7 was treated as a
   box to tick; writing the asserts found a live fallback rewriting 17% of one
   item's physics, and a non-laminar regime asserted as laminar. Write them
   first, not last.
5. **The item pool is not yet regenerated** and should not be until after
   Phase 6 (D-008). 120 of 1,350 published items are affected; ~65% need only a
   re-score. Whether the archived generations still exist is **D-003, open**.
