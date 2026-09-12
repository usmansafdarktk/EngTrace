# Phase C3 corrections — item-pool impact

The C3 phase **corrected nothing**: it changed comments, recorded 19 measured-wrong values
in `phaseC3_residual_register.md` §1, and left every one of them standing, because
correcting a constant moves the item pool and that was the repo owner's call. That call
has now been made, and this document measures what the corrections did.

**P6 event: YES, on 13 of 150 templates. No template gained a generation error.**

| | templates | instances |
|---|---:|---:|
| question changed | **13** | 824 |
| answer body changed (**re-score and re-inference**) | 11 | 713 |
| solution changed | 12 | 765 |
| generation errors | 0 | 0 before, 0 after |
| unchanged in every respect | 137 | — |

## The instrument

`tests/constants_integrity/c3_instance_dump.py`, 150 templates × 300 seeds, dumped from
two git worktrees in **two separate processes** — the script refuses two dumps sharing a
root, which is the in-process-reload trap that reports every instance identical:

```
git -c core.longpaths=true worktree add --detach <wt_c3head> 4eb6488
python -m tests.constants_integrity.c3_instance_dump <wt_c3head> c3_before.json 300
python -m tests.constants_integrity.c3_instance_dump .          c3_after.json  300
python -m tests.constants_integrity.c3_instance_dump --diff c3_before.json c3_after.json
```

`before` is **4eb6488**, the Phase C3 merge; `after` is the working tree.
Totals from the diff: `{'q': 824, 'ans': 713, 'sol': 765, 'err_before': 0, 'err_after': 0}`.

## What moved

| template | question | answer | solution | driven by |
|---|---:|---:|---:|---|
| `gas_viscosity_kinetic_theory` | 285/300 | 285/300 | 285/300 | `GAS_MOLECULAR_PARAMS` σ re-sourced to Svehla |
| `multi_segment_rod` | 129/300 | 125/300 | 129/300 | `MATERIAL_PROPERTIES` E |
| `composite_shafts_series` | 102/300 | 102/300 | 102/300 | `MATERIAL_PROPERTIES` E, G |
| `poissons_ratio` | 64/300 | 60/300 | 64/300 | `MATERIAL_PROPERTIES` ν |
| `statically_indeterminate_shaft` | 60/300 | **0** | **0** | restated, not used |
| `axial_deformation` | 43/300 | 41/300 | 43/300 | `MATERIAL_PROPERTIES` E |
| `statically_indeterminate` | 41/300 | **0** | 41/300 | restated, not used in the answer |
| `angle_of_twist` | 32/300 | 32/300 | 32/300 | `MATERIAL_PROPERTIES`, `SHEAR_MODULUS_VALUES` |
| `floating_object_submersion_depth` | 22/300 | 21/300 | 22/300 | `FLUID_DENSITIES` |
| `wave_parameters_basic` | 13/300 | 14/300 | 14/300 | `MEDIA_VELOCITIES` Helium, GaP |
| `basic_buoyant_force` | 11/300 | 11/300 | 11/300 | `FLUID_DENSITIES` |
| `hydrostatic_force_on_plane` | 11/300 | 11/300 | 11/300 | `FLUID_DENSITIES` |
| `hydrostatic_pressure_at_depth` | 11/300 | 11/300 | 11/300 | `FLUID_DENSITIES` |

## 1. `gas_viscosity_kinetic_theory` — the intended event, and its own check

285 of 300 instances moved. This is the σ column of `GAS_MOLECULAR_PARAMS` re-sourced to
Svehla (NASA TR R-132, Table I(a)), and the item pool reproduces the physics that
justified the correction:

```
seed 0: 'dynamic viscosity of Butane (C₄H₁₀) at 465 K is **1.537e-05 Pa·s**'
     -> 'dynamic viscosity of Butane (C₄H₁₀) at 465 K is **2.094e-05 Pa·s**'
```

n-butane's σ went 5.47 → 4.687 Å. Chapman–Enskog has μ ∝ σ⁻², so the ratio predicted from
the constants alone is (5.47/4.687)² = **1.362**, and 1.537 × 1.362 = **2.093**. The
emitted answer is 2.094. The correction was argued from Chapman–Enskog with the Neufeld
collision integral against the NIST isobars, and the emitted items agree with that argument
to the last digit — computed by the template, not by the patch that made the change.

This template is why the σ column was `[KNOWN-DEFECTIVE]` rather than a citation nit: it
prints σ into the question stem and computes a viscosity from it. The old items asked
about a named gas using a σ that was not that gas's, from a compilation the table could not
name.

## 2. The two that need re-inference and **not** re-scoring

`statically_indeterminate` (41 questions, **0 answers**) and
`statically_indeterminate_shaft` (60 questions, **0 answers, 0 solutions**).

A corrected constant is **restated in the question and not used in the arithmetic** on
those seeds. A model saw a question containing one modulus and will now see another; the
correct answer is the number it always was. This is the distinction the two-process dump
exists to make, and collapsing it would either force a needless re-score or hide a real
one.

## 3. What did not move, and why that is the result

137 templates are byte-identical in question, solution and answer. In particular:

- **No generation error appeared anywhere** (0 → 0 across 45,000 instances). A corrected
  constant that put a template outside its own domain would show here as an exception.
- The civil `WATER_KINEMATIC_VISCOSITY_M2_S` correction (1.004e-6 → 1.003e-6) moved
  **nothing**: no template consumes it. It is corrected because it was wrong, not because
  anything depended on it.
- The `MANOMETER_FLUIDS` and `MATERIAL_DENSITIES` rows were not touched in this tranche, so
  the manometer and density templates are unchanged.

## 4. Regression

Identical to the recorded baseline on every gate:

| gate | baseline | now |
|---|---|---|
| T1 / T2 / T4 / T5 / T7 / T8 | 29 / 0 / 0 / 66 / 83 / 0 | 29 / 0 / 0 / 66 / 83 / 0 |
| T6 | 142 failing | 142 failing |
| `phase5_contract_scan` | 150/150 clean | 150/150 clean |
| `cross_pair` | PASS | PASS |
| `audit_3_8` | 80/80 clean | 80/80 clean |
| `constants_integrity` | all pass | all pass |

**T6's 142 is weak evidence and is recorded as such.** Its baseline is stale corpus-wide
(D-043) and it fails 142 templates regardless, so "T6 did not move" does not mean the
distributions did not move — 13 of them demonstrably did. The baseline was **not**
regenerated; `run.py --baseline` was never invoked. The instance dump above is the
measurement; T6 is only the statement that its own gate is where it was.

---

# Tranche 2 — Decision 7, and 23 tags that moved nothing

**P6 event: YES, on 1 of 150 templates.** `before` is **7f63348**, the tranche-1 head;
`after` is the working tree. Same instrument, same two-process discipline.

Totals: `{'q': 31, 'ans': 31, 'sol': 31, 'err_before': 0, 'err_after': 0}`

| template | question | answer | solution |
|---|---:|---:|---:|
| `sensible_heat_temp_dependent_cp` | 31/300 | 31/300 | 31/300 |

```
seed 2: 'The heat required is **3.74 kJ**.'  ->  'The heat required is **3.32 kJ**.'
```

`CP_PARAMS` declares `@domain T_lo=298..1500 K`, and the template's liquid branch drew
`T1 = uniform(280.0, 300.0)`, reaching 281.67 K — the polynomial integrated 16 K below the
interval it was fitted over, on 27 of 300 seeds. Raised to `uniform(298.15, 318.0)`:
minimum T1 **280.75 → 298.49 K**, seeds below the floor **27 → 0**. Less heat on those
seeds, because the interval no longer extends below 298 K.

**Question, answer and solution move together here**, unlike tranche 1's
`statically_indeterminate` pair where the question moved and the answer did not. A changed
integration bound is used by the arithmetic; a restated modulus need not be. These 31
instances need **both** re-inference and re-scoring.

**The span was kept ~20 K wide** rather than clipped to [298.15, 300.0]. Distinct
questions **300 → 300**. Collapsing a 20 K draw to 1.85 K would have shrunk the answer
space, which `checks/t6_distribution` calls a downgrade rather than a fix — so the obvious
minimal edit was the wrong one.

## The 23 tags that moved nothing

18 PubChem citations and 5 wood residuals were applied in the same tranche and appear
**nowhere** in the diff. That is the evidence they are comment-only: had any of them
disturbed a value, this dump is where it would show, across 45,000 instances.

## Regression

Identical to baseline on every gate: T1 29, T2 0, T4 0, T5 66, T7 83, T8 0, T6 142,
`phase5_contract_scan` 150/150, `cross_pair` PASS, `audit_3_8` 80/80.

The consumer-domain ratchet moved **2 excursions / 2 registered → 1 / 1**. It ratchets in
both directions — a listed excursion that no longer occurs is a failure — so it was the
suite, after the fix, that said to remove the `CP_PARAMS` line from `domain_findings.txt`.
The removal followed the evidence rather than accompanying the change.

---

# Tranche 3 — a correction that moved nothing, and why that is a result

**P6 event: NO.** `before` is **8c90537**; `after` is the working tree.
Totals: `{'q': 0, 'ans': 0, 'sol': 0, 'err_before': 0, 'err_after': 0}` — **0 templates
moved**.

`MATERIAL_DENSITIES['Tungsten']` 19600 → 19300 is a real correction: five independent
PubChem entries print 19.3 g/cm³ (CAMEO, ILO-WHO, OSHA, PAC @25 °C, NIOSH) and the
committed value was +1.55% against every one of them. Its tag had said "no source on disk
for this material", which the acquisition had made false.

**A value was wrong, and correcting it changed no item. Those are two different facts**,
and collapsing them would misreport either the defect or the impact. The dump is the
evidence; the reason is in the code:

- `template_floating_object_submersion_depth` draws from `MATERIAL_DENSITIES` and then
  **rejects any draw that does not float**, re-drawing up to 20 times. At 19300 kg/m³
  tungsten sinks in every fluid in the table, so it is drawn and always discarded; the
  fallback hard-names "Pine Wood".
- `template_basic_buoyant_force` uses its material as a **word in the prose** ("A solid
  tungsten cylinder…") and computes from `FLUID_DENSITIES`. The name can be emitted while
  the density never is.

A null result is reported here with its mechanism rather than as a bare zero, because
"nothing moved" and "nothing was wrong" look identical in a totals line and are not.

**A claim made earlier in this work and withdrawn:** that `OBJECT_MATERIALS` does not
contain tungsten. It does — it is the last entry. The reason tungsten's density is
unreachable is the float rejection and the prose-only use, not absence from that list.

## The 20 tags that moved nothing

19 fuel/oil rows and the aluminium row are comment-only and appear nowhere in the diff.
Nonane, decane and dodecane were acquired to **re-point** those rows, and reading them
settled that plan against itself: at 293.15 K they are 718.03, 730.41 and 749.44 kg/m³
against Kerosene 810, Diesel 850, SAE 30 917. Pointing kerosene at decane would move it
810 → 730 — an 11% error — because a real cut carries aromatics and cycloalkanes a pure
n-alkane does not. **Sourcing a row must not change what the row names.** The tags now
carry the measured bracket and record that the search was made.

Gasoline (726) lands *inside* that span, and its tag says at length that this is
coincidence and not support: a petrol cut is C4–C12 with most of its mass in C5–C8,
lighter than every alkane in the bracket.

**Aluminium is deliberately not corrected.** The row holds 2710 and names no alloy;
PubChem prints 2.7 (pure, 2700), MIL-HDBK-5J prints 2713 for 6061 and 2768 for 2024.
Replacing 2710 with the pure-metal figure would be choosing which aluminium the row means
and calling it sourcing.

## Regression

Unchanged on every gate: T1 29, T2 0, T4 0, T5 66, T7 83, T8 0, T6 142, contract scan
150/150, `cross_pair` PASS, `audit_3_8` 80/80, resolver 224 resolved of 485 tags with
0 LEGACY, plausibility 72 checks.
