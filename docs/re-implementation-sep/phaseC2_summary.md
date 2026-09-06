# Phase C2 — Summary and close-out

**Status: COMPLETE WITH FINDINGS.** · **Date:** 2026-09-06
**Branch:** `redesign/phaseC2-thermochemistry-constants`

C2 re-grounded the three chemical thermochemistry tables — `CP_PARAMS` (33
rows), `HEATS_OF_FORMATION` (21) and `COMBUSTION_REACTIONS` (11). It is
**sync point S1**: Track A Phase 2 could not correctly re-baseline
`template_adiabatic_flame_temperature` until it landed.

> **Headline.** The adiabatic flame temperature moves from **1130–1493 K** to
> **2101–2844 K**. Methane/air now computes **2311 K** against a correct
> reference of **2326 K** — within 0.6%.

---

## 1. The defect was four defects

The spec records one: *"B and C coefficients 10× too large (`E-2`/`E-5` where
Smith–Van Ness gives `E-3`/`E-6`)"*. That is real and affects most rows. It is
also **one of four**, and correcting only it leaves seven rows still wrong.

| Class | What it is | Rows |
|---|---|---|
| **1 — exponent slip** | `B` written `E-2` where the source column heading is `10³B`; `C` written `E-5` where the heading is `10⁶C`. A column heading transcribed as part of the value. | most |
| **2 — row simply wrong** | coefficients that match no source | `O2(g)`, `C2H2(g)`, `C2H5OH(g)`, `C3H6O`, `C6H14` |
| **3 — column shift** | `D`'s mantissa −0.792 written into the `C` column, an unrelated value left in `D` | `NO2(g)` |
| **4 — wrong phase** | gas-phase coefficients under a liquid key | `C6H6(l)`, `C7H8(l)`, `C6H14(l)`, `C3H6O(l)` |

Class 1 alone took Cp(N₂, 300 K) from 42.4 to 29.12 J/mol·K against NIST's
29.12. But `O2` still read 34.71 against 29.38, and `NO2` 40.27 against 36.97,
until each was replaced individually.

---

## 2. Method

Every row is checked against the **NIST Chemistry WebBook (SRD 69)**, retrieved
to `docs/references/nist_webbook/` so each `[ON-DISK]` citation resolves to a
file a reviewer can open without network access or a library.

NIST publishes the **Shomate** form; the constants use the **Smith–Van Ness**
form. They are different fits to the same underlying thermochemistry, so
agreement between them is evidence rather than a tautology. Two good fits agree
to about 3%; the acceptance gate is 8%, which passes everything correct and
still catches every defect found — the smallest was `NO2` at 14.8%.

A second, independent check was available for some rows: Smith–Van Ness
publishes a `Cp₂₉₈/R` column alongside A, B, C, D, so the table self-checks.
Where available it is decisive — methane computes **4.2171** against a published
4.217, ethane **6.3686** against 6.369, propane **9.0109** against 9.011, and
`O2` computes **4.175** against a published **3.535**, which is how that row was
condemned.

---

## 3. What changed

| Row | Action | Evidence |
|---|---|---|
| most rows | exponent corrected | e.g. N₂ 42.4 → 29.12 vs NIST 29.12 |
| `O2(g)` | **replaced** — A, B, C and D all wrong | 34.71 → 29.39 vs NIST 29.38 |
| `NO2(g)` | **replaced** — column shift | 40.27 → 36.98 vs NIST 36.97 |
| `C2H2(g)` | **refitted** to NIST Shomate | was 68.39 vs NIST 44.04; fit worst error 3.10% |
| `C2H5OH(g)` | **refitted** to TRC 1997 | was 74.40 vs 65.21; worst 1.79% |
| `C3H6O` | **refitted** to Chao 1986, re-keyed `(l)`→`(g)` | was 81.32 vs 75.02; worst 1.11% |
| `C6H14` | **refitted** to Scott 1974, re-keyed `(l)`→`(g)` | was 125.70 vs 142.60; worst 2.79% |
| `C6H6`, `C7H8` | re-keyed `(l)`→`(g)` | held gas values; benzene computed 85.3 where NIST *liquid* is 135.69 |
| `NaCl(s)` | **kept** | its 5.5% gap to NIST is a genuine disagreement between two fits, not a transcription defect |
| `H2SO4(l)`, `CaCO3(s)` | **tagged `[UNVERIFIED]`** | NIST's free tier paywalls both |

Coefficients are now written with the exponent attached to each value — `B` as
`E-3`, `C` as `E-6`, `D` as `E5` — precisely because the source tabulates them
as the column headings `10³B`, `10⁶C`, `10⁻⁵D`. Writing them this way makes
Class 1 impossible to repeat silently.

**`HEATS_OF_FORMATION`: all 21 verified, none outside its stated uncertainty.**
The spec listed the table as "Unverified"; it is sound. Largest gap is methanol
gas at 4.3 kJ/mol against a NIST value carrying ±10.

**`COMBUSTION_REACTIONS`: all 11 balance element-by-element** and all carry the
theoretical-air ratio N₂/O₂ = 3.760. `REACTIONS`: all 4 balance. Every species
a reaction names has a heat of formation.

---

## 4. The flame temperature, and a correction to the spec

The spec's exit gate asks for *"methane/air flame temperature within the
literature range"* and quotes **~2200 K**. **That is the wrong reference for
this item.**

`template_adiabatic_flame_temperature` models a single balanced reaction with
**no dissociation**. For air–methane, stoichiometric, reactants at 298.15 K:

| model | T_ad |
|---|---|
| **complete combustion, no dissociation** ← what the template computes | **2326.35 K** |
| chemical equilibrium with dissociation (GRI-Mech 3.0) | 2224.25 K |

*(ETASR / arXiv:2503.11826, saved to `docs/references/README.md`.)*

The template computes **2311 K** — within **0.6%** of the correct reference. The
~100 K gap to the spec's figure is the dissociation effect the item's model
deliberately excludes, and every fuel shows the same consistent offset, which is
the signature of a modelling assumption rather than a data defect.

**`SPEC-CHANGE`:** the C2.4 gate should read 2326 K, or state that ~2200 K is
the equilibrium value and not comparable to this item's model.

---

## 5. Deliverables

| # | Deliverable | Where | Status |
|---|---|---|---|
| C2.1 | Three tables re-derived with `[ON-DISK]` citations | `constants.py`, per row | ✅ 31 species on disk; 2 tagged `[UNVERIFIED]` |
| C2.2 | Physical-plausibility test suite | `tests/constants_integrity/test_chemical_thermochemistry.py` | ✅ **208 checks, all passing** |
| C2.3 | Before/after impact for the 5 affected templates | §6 | ✅ |
| C2.4 | Methane/air near the literature value | §4 | ✅ 2311 K vs 2326 K |
| C2.R | Independent review (Reviewer H) + R4 triage | `reviews/` | see §8 |

---

## 6. Item-pool impact (C2.3)

300 seeds per template, `master` vs this branch:

| Template | identical | gold-only (re-score) | question changed (re-inference) |
|---|---:|---:|---:|
| `adiabatic_flame_temperature` | 0/300 | **300/300** | 0/300 |
| `sensible_heat_temp_dependent_cp` | 26/300 | 0/300 | **274/300** |
| `heat_of_reaction_formation` | 300/300 | — | — |
| `latent_heat_vaporization` | 300/300 | — | — |
| `sensible_heat_constant_cp` | 300/300 | — | — |

The two are affected differently, and the second is worth stating plainly.
`sensible_heat_temp_dependent_cp` **prints the coefficients in its question**:

```
before:  Cp/R = 3.156 + 6.230e-03*T + 1.510e+04*T^-2
after:   Cp/R = 3.156 + 6.230e-04*T + 1.510e+04*T^-2
```

Those 274 items were not merely mis-answered. They **stated a heat capacity ten
times too large and were then self-consistent about it**, so a model could be
scored correct for doing sound arithmetic on physics that does not exist. The 26
unchanged instances are the monatomic gases, whose row was already exact.

`adiabatic_flame_temperature` keeps its question and gets a corrected answer, so
it is a **re-score**, not re-inference — assuming the archived generations still
exist, which is **D-003, still open**.

Per D-008, the pool is regenerated **once, after Phase 6**. This note counts the
cost in advance.

---

## 7. Source deviation, recorded not hidden

The spec names **Smith–Van Ness–Abbott** as primary with NIST as cross-check.
**No fetchable, citable copy of SVN Table C.1 could be obtained** — the
accessible copies are Scribd and SlideShare, which cannot be downloaded and are
not legitimate "edition + page" citations. Fabricating a page number would be
exactly the failure this project exists to prevent.

So **NIST is the primary on-disk source** and the SVN functional form is
retained. Four rows are least-squares fits to NIST tables rather than
transcriptions; each reports its residual. Those four no longer match a textbook
table a student might hold — a real cost, accepted deliberately and reversible
if a citable SVN copy is obtained.

---

## 8. What the phase's own suite caught — including my own error

The C2.2 suite refuses to pass a row that has no cited reference, rather than
skipping it. That refusal earned itself immediately, finding three defects the
manual pass had missed:

- **`H2O(l)` and `H2SO4(l)` — 10× errors I introduced while rewriting the
  table.** I transcribed `B`'s mantissa as 12.5 where the original read 1.25.
  That is the Class-1 defect, committed by the person removing it. `H2O(l)`
  computed 102.09 against NIST's 75.30.
- **`C6H14`** — 11.8% below the NIST *gas* value, so the row was wrong as a gas
  as well as mis-keyed as a liquid.

An audit now compares every non-replaced row's mantissa against the original,
confirming only the exponent moved: 3 mismatches found, 3 fixed, 0 remaining.

> **Carry this forward:** a table that has just been "corrected" is exactly when
> a fresh transcription error is most likely, and least likely to be looked for.

---

## 9. Residual risk register

| Item | Status |
|---|---|
| `H2SO4(l)`, `CaCO3(s)` | `[UNVERIFIED]` — NIST's free tier paywalls both. Tagged, not quietly accepted. |
| `NaCl(s)` | 5.5% from NIST; a genuine fit disagreement, kept as the correct source row |
| four refitted rows | no longer match a textbook table; residuals 1.1–3.1% |
| SVN page citations | unobtainable; NIST used instead (§7) |
| liquid branch | now only `H2O(l)` and `CH3OH(l)`, both verified as genuine liquids |

---

## 10. What Phase 2 inherits

1. **S1 is satisfied.** `template_adiabatic_flame_temperature` can now be
   re-baselined; its flame temperatures are physically right.
2. **Phase 2's other three templates were never blocked** by C2 —
   `levenspiel_plot_interpretation`, `pfr_volume_changing_rate` and
   `vibration_isolator_design` reference none of these tables.
3. **The C2.2 suite is a template for C3.** C3 covers ~400 more values across
   four branches; the same two-independent-checks method, and the same refusal
   to pass an uncited row, should carry over.
4. **Expect more than one defect class per table.** C2's spec entry described
   one and the table had four. C3's entries are one line each.
