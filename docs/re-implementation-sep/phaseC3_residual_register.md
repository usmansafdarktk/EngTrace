# Phase C3.5 — residual register

**What this is.** Spec §C3.5: every constant C3 could not resolve to an artefact, or
resolved and found wrong, with the decision it needs from the repo owner.
**Listing is not acceptance.** Every line is an open defect or an open question.

**Rewritten after the C3 review.** Reviewer G's F-2 found that §C1.2 states register
membership as the check for both residual classes and `R6` never performed it — so 20
residuals existed here only inside a remainder integer ("not covered by any phrase above:
24"). This version **itemises every one**: all 19 defects individually, and all 198
unverified rows by the table that carries them, so membership is checkable rather than
asserted. Three reviewer findings also overturned entries that were here; §5 records what
I got wrong.

**Re-derive, don't trust these lines:**

```
python -m tests.constants_integrity.test_citations_resolve
python -m tests.constants_integrity.test_consumer_domains
grep -c "^\s*# \[KNOWN-DEFECTIVE\]" data/templates/branches/*/constants.py   # 19
grep -c "^\s*# \[UNVERIFIED\]"      data/templates/branches/*/constants.py   # 198
```

(The greps anchor on a line that STARTS a tag. An earlier version of this file counted
bare occurrences, and my own explanatory prose — sentences saying a row was *no longer*
a class — inflated it by three, one of which the resolver was parsing as a real tag.)

**At the review fixes:** 452 tags — **156 resolved, of which only 109 are value
comparisons**; the other 47 are C2-form CAS tags whose check is the CAS number (G F-3).
19 locator-only, 277 stated, **0 LEGACY**, 0 unresolvable from a clone. 108 tables all
declare `@kind`, `@units`, `@domain`. **15 `[DERIVED]` tags are never recomputed** — the
`UNEXECUTED` class §C1.2 promises does not exist (G F-4).

| branch | tags | resolved | KNOWN-DEFECTIVE | UNVERIFIED |
|---|---|---|---|---|
| mechanical | 208 | 21 | 15 | 172 |
| chemical | 138 | 105 | 1 | 10 |
| civil | 34 | 12 | 1 | 8 |
| electrical | 34 | 16 | 2 | 4 |
| industrial | 38 | 2 | 0 | 4 |
| **total** | **452** | **156** | **19** | **198** |

---

## 1. Measured wrong — all 19, itemised

Each was compared with a named artefact at the table's own conditions and is not a
rounding of it. **No value was changed**: correcting one moves the item pool (P6).

| # | branch.table[row] | field | committed | artefact | error |
|---|---|---|---|---|---|
| 1 | chemical `GAS_MOLECULAR_PARAMS` | σ column | see below | Svehla TR R-132 | **n-butane −26.1%, SF₆ −12.7%** |
| 2 | electrical `MEDIA_VELOCITIES[Helium]` | n | 1.000036 | 1.0000349 (Ermolov, Mansfield) | +1.1e-06 on n, **+3.20% on n−1** |
| 3 | electrical `MEDIA_VELOCITIES[GaP]` | n | 3.5 | 3.3616–3.4068, five datasets | +2.74% to +4.12% |
| 4 | mechanical `MATERIAL_PROPERTIES[Steel]` | nu | 0.30 | p.62 AISI 1025, µ 0.32 | **−6.25%** |
| 5 | mechanical `MATERIAL_PROPERTIES[Stainless Steel]` | E_ksi | 27500 | p.277, 29.0e3 ksi | −5.2% |
| 6 | " | E_GPa | 190 | = 199.9 GPa | −5.0% |
| 7 | " | nu | 0.30 | mu 0.27 | +11.1% |
| 8 | mechanical `MATERIAL_PROPERTIES[Aluminum 6061-T6]` | E_ksi | 10000 | p.566, 9.9e3 ksi | +1.0% |
| 9 | " | E_GPa | 68.9 | = 68.26 GPa | +0.9% |
| 10 | mechanical `MATERIAL_PROPERTIES[Titanium]` | E_ksi | 16800 | p.899, 15.5e3 ksi | +8.4% |
| 11 | " | E_GPa | 116 | = 106.9 GPa | +8.5% |
| 12 | mechanical `MATERIAL_PROPERTIES[Ti-6Al-4V]` | E_ksi | 16500 | p.945, 16.0e3 ksi | +3.1% |
| 13 | " | E_GPa | 114 | = 110.316 GPa | **+3.3395%** |
| 14 | " | nu | 0.34 | mu 0.31 | +9.7% |
| 15 | mechanical `MATERIAL_PROPERTIES[Magnesium]` | nu | 0.29 | p.841, mu 0.35 | **−17.1%** |
| 16 | mechanical `SHEAR_MODULUS_VALUES[Ti-6Al-4V]` | value | 41.4 | p.945, 6.2e3 ksi = 42.75 GPa | −3.2% |
| 17 | mechanical `FLUID_DENSITIES[Benzene]` | value | 876 | NIST 878.92 @ 293.15 K | −0.33% |
| 18 | mechanical `FLUID_DENSITIES[R-134a]` | value | 1206 | NIST 1206.7 @ 298.15 K | −0.06% |
| 19 | civil `WATER_KINEMATIC_VISCOSITY_M2_S` | value | 1.004e-6 | mu/rho = 1.003396e-6 | +0.060% |

**Rows 4–15 are one decision, not twelve.** E and nu for five alloy families disagree with
the MIL-HDBK-5J this repo ships. AZ31B's nu at −17.1% is the largest measured error in the
corpus; Steel's nu at −6.25% is in its most-used material.

**Row 1 is the newest and the sharpest.** `GAS_MOLECULAR_PARAMS` pairs Svehla's ε/k with a
σ that is not Svehla's. Chapman–Enskog with the Neufeld collision integral at 293.15 K
against the NIST isobars: n-butane **−26.1%**, SF₆ **−12.7%**, propane +3.5%, methane
−0.7%, ethane +2.9%. Substituting Svehla's σ brings **all five inside 2%**. Since
μ ∝ σ⁻² and the template prints σ into the question stem, this is a wrong viscosity, not a
citation problem.

---

## 2. Unchecked — all 198, by the table that carries them

Not known to be wrong; known to be **unchecked**. 27 tables.

| branch.table | n | why |
|---|---:|---|
| mechanical `MATERIAL_DENSITIES` | 54 | one number for a material whose density depends on species/alloy/condition |
| mechanical `PIPE_FLUIDS` | 32 | no source on disk |
| mechanical `FLUID_DENSITIES` | 30 | mixtures, solutions and commercial products with no single composition |
| mechanical `MATERIAL_PROPERTIES` | 22 | no MIL-HDBK-5J design table under the names searched |
| mechanical `MANOMETER_FLUIDS` | 20 | not in the 36-fluid NIST set; WebBook pages carry no saturated volumes |
| mechanical `SHEAR_MODULUS_VALUES` | 14 | no handbook table for the named grade |
| chemical `REAL_FLUID_DATA` | 3 | R-410A (a blend), ethanol, acetone — absent from the NIST set |
| chemical `AIR_COMPOSITION` | 1 | US Standard Atmosphere 1976 not on disk |
| chemical `COMMON_GASES` | 1 | textbook values; 10 (row, field) pairs do equal NIST |
| chemical `COMMON_LIQUIDS` | 1 | textbook values; 3 pairs equal NIST — **diagnosis open, §5** |
| chemical `CRITICAL_PROPERTIES` | 1 | 3 of 27 rows resolve to a WebBook page; 24 have no artefact |
| chemical `POWER_LAW_FLUIDS` | 1 | K and n for ketchup, mayonnaise, toothpaste — no primary source |
| chemical `SUBSTANCES_FOR_HEATING` | 1 | Cp in J/(g·K) with no temperature stated |
| chemical `SUBSTANCES_FOR_VAPORIZATION` | 1 | 4 rows cite a WebBook page — **count corrected, §5** |
| civil `STEEL_E_KSI`, `STEEL_E_GPA`, `STEEL_FY_KSI`, `STEEL_FY_MPA` | 4 | AISC Manual and ASTM not on disk |
| civil `CONCRETE_EC_COEFF_PSI`, `CONCRETE_EC_COEFF_MPA` | 2 | ACI 318-19 not on disk |
| civil `LIVE_LOADS_PSF`, `LIVE_LOADS_KPA` | 2 | ASCE 7-22 not on disk; the SI column is a mixed soft conversion |
| electrical `MEDIA_VELOCITIES` | 4 | PTFE and polyethylene absent from the archive; "amorphous glass" names a class; the tissue row is a microwave √ε_r, not an index |
| industrial `IEC60063_E24_5PCT`, `IEC60063_E12_10PCT`, `RESISTOR_SERIES_BY_TOLERANCE` | 3 | IEC 60063 not on disk |
| industrial `HARD_ANODIZE_THICKNESS` | 1 | MIL-A-8625 not on disk |

**The mechanical 172 are ~87% of everything unchecked in the corpus**, and six tables
carry all of it. That is one decision taken six times, not 172.

**Named standards absent:** AISC Manual, ASTM A992/A36, ACI 318-19, ASCE 7-22, IEC 60063,
MIL-A-8625, US Standard Atmosphere 1976. MIL-STD-105E was in this list until C3.7 found a
hashed public copy already under `docs/references/` (D-076).

---

## 3. Substance defects — the row names the wrong thing (D-033)

| row | branch | finding |
|---|---|---|
| `R-410A` | chemical | a blend, absent from the 36 pure fluids; needs a mixture source or a replacement |
| `Tellurium Mercury` | mechanical | WebBook returns "Name Not Found"; whether it means mercury telluride is not established |
| `Tungsten Hexafluoride` | mechanical | listed as a manometer liquid; the on-disk Antoine parameters put P = 1.074 bar at 290.4 K, so it is a **gas** at room temperature |
| `2024-T4` | mechanical | MIL-HDBK-5J **does** carry the table (p.374); its text layer does not parse |
| `Mercury` ΔHvap | chemical | the WebBook phase-change page has no vaporization section |
| `Cork`, `Cork Board`, `Bamboo` | mechanical | **not wood** — a bark tissue and a grass; the wood reason never applied (H-mech F4) |

---

## 4. Tables that look sourced and are not

| table | finding |
|---|---|
| chemical `GAS_MOLECULAR_PARAMS` | §1 row 1 — now `[KNOWN-DEFECTIVE]` |
| chemical `COMMON_LIQUIDS` | organic viscosities 4.62–7.16% low at the declared 20 °C; **diagnosis open**, §5 |
| chemical `COMMON_GASES` | the steam row is **+25.00%** against NIST at the 373.15 K its own comment names, and within 0.25% of the ideal gas at 293.15 K — 20 °C vapour mislabelled as 100 °C steam |
| chemical `HEATS_OF_FORMATION[CH4(g)]` | the tag cites −74.6 ± 0.3 and the row uses −74.8 (inside the uncertainty, but the tag names a measurement the value is not) |
| civil `LIVE_LOADS_KPA` | 4 of 5 rows are the soft conversion 1 psf = 0.048 kPa exactly; the corridor row is neither that nor SP 811 |
| chemical `SUBSTANCES_FOR_VAPORIZATION` | **11 of 11** derivable rows agree with a NIST saturation derivation to within 1.454% — near-derivation, not derivation (D-077) |
| chemical `CRITICAL_PROPERTIES` | 3 of 27 rows resolve (ethanol 513.9 K, acetone 508.2 K, p-xylene 616.2 K, exact) |
| mechanical `MATERIAL_PROPERTIES` (Nickel) | the "no design table" reason cited a **copper-base** search set; the handbook carries 22 nickel-base captions (A-286 p.1045, Hastelloy X p.1061, Inconel 600 pp.1067–1070), none commercially pure — **conclusion open** |

---

## 5. What the review overturned in this file

Three entries here were wrong, and a fourth was unsafe. Recorded because a register that
quietly absorbs its corrections is worth no more than the claims it corrected.

- **The two civil unit-weight rows are gone from §1.** They were listed as defects on the
  claim that 9.81 kN/m³ "implies ρ = 1000.34 kg/m³ — water near 4 °C". **That is
  physically impossible**: liquid water at 1 atm peaks at 999.975 kg/m³, at 3.98 °C. They
  are the ρ = 1000 exactly convention, which Das states on PDF p.92 — one page before the
  Table 3.1 this file already cites. They now resolve as `[ON-DISK:LOCAL-ONLY]`.
- **`SUBSTANCES_FOR_VAPORIZATION` said "7 of 12".** Measured: 11 of the 11 rows with a
  saturation file agree within 1.454%, eight within 0.2%. Mercury is the twelfth and has
  no file — so *no row misses*, and "the rows that miss" was the stated reason for the
  class.
- **`COMMON_LIQUIDS`' diagnosis is withdrawn.** I called the viscosity gap "a source
  difference, not scatter". The review proposes a 25 °C column. **Neither is
  established**, and I could not test theirs: no on-disk grid carries 298.15 K. n-hexane,
  the one species with a row near 25 °C, moves −7.27% → −1.10%; water, the control, is
  +0.04% at 20 °C. The remedy I prescribed would have discarded rows that may be right at
  a temperature the table does not state.
- **Ti-6Al-4V E_GPa was printed as +3.4%** where the tag says +3.3% and the truth is
  +3.3395%. A wrong number inside a defect record is worse than an untagged value.

---

## 6. Open questions for the repo owner

1. The twelve mechanical alloy defects (§1 rows 4–15) — correct, or restate the grades.
2. `GAS_MOLECULAR_PARAMS`' σ column — re-source from Svehla (moves items; μ ∝ σ⁻²).
3. Six mechanical tables carry 172 unchecked rows — acquire sources, or accept.
4. `COMMON_LIQUIDS`' conditions — fetch a 298.15 K grid and settle §5's open diagnosis.
5. **`tol=` has no rule for its SIZE.** D-074 says when a tolerance may be used, never how
   large. 4 of the 10 `tol=` tags are their own residual rounded up (three electrical rows
   and `FLUID_DENSITIES['Liquid Propane']`); the other six are half a unit in the
   handbook's last digit, as intended. Each of the four now states what does and does not
   set it; none was re-sized. A resolver rule requiring a stated basis is a SPEC-CHANGE.
6. **15 `[DERIVED]` tags are never recomputed** (G F-4) — the `UNEXECUTED` class is unbuilt.
7. `CP_PARAMS` is integrated from **281.67 K**, below the 298 K floor its fit declares —
   registered in `domain_findings.txt`; raise the draw, refit, or restate the floor.

## 7. Registered elsewhere

| register | holds |
|---|---|
| `tests/constants_integrity/domain_findings.txt` | 2 consumer-domain excursions: `MEDIA_VELOCITIES` optical indices used for 50–500 MHz radio, and `CP_PARAMS` integrated below its floor |
| `phaseC3_literal_copies.md` | C3.8 — declared `@copied-in` copies and `_SCS_COMBOS` |
| `phaseC3_inline_windows.md` | C3.9 — named-entity facts written as template literals |

## 8. What C3 did not do

Nothing in §1 was corrected. C3 changed comments, not values, with three exceptions, all
recorded as P6 events and all measured: `REAL_FLUID_DATA` re-derived from NIST,
`CONTROL_CHART_FACTORS`' 29 + 2 corrected last digits, and two template edits. The
item-pool consequence is measured in `phaseC3_item_pool_impact.md`: **2 of 150 templates
moved**.
