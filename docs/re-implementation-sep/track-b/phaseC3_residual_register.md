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
python -m tests.constants_integrity.test_registers_reconcile            # this file vs the tree
python -m tests.constants_integrity.test_citations_resolve
python -m tests.constants_integrity.test_consumer_domains
grep -c "^\s*# \[KNOWN-DEFECTIVE\]" data/templates/branches/*/constants.py   # 0  (19 at the review)
grep -c "^\s*# \[UNVERIFIED\]"      data/templates/branches/*/constants.py   # 170 (198 at the review)
```

The first command is the one that keeps the rest honest: it parses the figures **this
document states** and re-derives each from the branch files, so a drift on either side
fails and the message names which side moved. It exists because every count here was
right when written and several were wrong when read.

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

**Now** (2026-09-16, from the resolver and the greps above) — the corrections tranche and
`269fc76`'s seven deletions have moved every column:

| branch | tags | resolved | KNOWN-DEFECTIVE | UNVERIFIED |
|---|---|---|---|---|
| mechanical | 202 | 60 | **0** | 142 |
| chemical | 175 | 140 | **0** | 12 |
| civil | 34 | 12 | **0** | 8 |
| electrical | 34 | 18 | **0** | 4 |
| industrial | 38 | 2 | **0** | 4 |
| **total** | **483** | **232** | **0** | **170** |

Of the 232 resolved, **47 are still C2-form CAS tags** whose check is the CAS number rather
than a value, so **185 are value comparisons** — G F-3's distinction survives the corrections
and is now printed by the resolver on every run. 20 locator-only, 231 stated, **0 LEGACY**,
0 unresolvable from a clone. `[DERIVED]`: 8 recomputed here, 1 by a named checker,
**7 UNEXECUTED** — the `UNEXECUTED` class G F-4 asked for now exists and reports itself.

---

## 1. Measured wrong — all 19, itemised, and **all 19 now corrected**

Each was compared with a named artefact at the table's own conditions and was not a
rounding of it. C3 left every one standing, because correcting a constant moves the item
pool (P6) and that was the repo owner's call.

**That call has been made, and all 19 are corrected.** The table below is kept as the
record of what was wrong; the branch files now carry the artefact's value. `[KNOWN-DEFECTIVE]`
is **zero corpus-wide**. The item-pool consequence is measured in
`phaseC3_corrections_item_pool_impact.md`: 13 of 150 templates moved, 824 questions and 713
answers, and **no template gained a generation error**.

Row 1 is corrected in part, and says so. Its σ column re-sourced **16 of 17** rows to
Svehla, which now carry 31 per-row `[ON-DISK]` citations the resolver re-reads from the
scan. Two elements are deliberately not written:

- **Xenon's σ** — it has no token in the OCR, and the positional read the C3 review
  recorded is not reproducible: indexing p.26's 20th number tuple returns 2.608 / 10.22,
  which is *helium's* second determination. Anchoring on its committed ε/k found no unique
  tuple either. It keeps its value under a row-level `[UNVERIFIED]`.
- **Ammonia's ε/k** — the scan prints it "55& 3" and the parser returns 55. By eye it is
  558.3, but by eye is not a parse.

A twentieth defect was found while closing these and is corrected with them:
`SHEAR_MODULUS_VALUES['Aluminum 2024-T4']`, 28.0 GPa against 27.58 (see §3).

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

## 2. Unchecked — **170**, by the table that carries them

Not known to be wrong; known to be **unchecked**.

**Measured now** (tag-starting lines, per branch): mechanical 142, chemical 12, civil 8,
electrical 4, industrial 4 — **170**. The per-table breakdown below is the **C3-era
snapshot** and is kept as the record of what the phase inherited; it no longer sums to the
current total.

**177 → 170, and the seven are deletions, not sourcing.** `269fc76` removed the seven
substance-defect rows of §3, owner-directed: `Cork`, `Cork Board`, `Bamboo`,
`Tellurium Mercury`, `Tungsten Hexafluoride` and `Tetrabromoethane` from mechanical, and
`R-410A` from chemical. Each had been deferred twice because deleting a row shifts every
key after it and so moves the item pool; with the pool being rebuilt that cost went to
zero. Mechanical 148 − 6 = 142, chemical 13 − 1 = 12.

The last change was **176 → 177**, and it is an improvement rather than a regression:
`SUBSTANCES_FOR_VAPORIZATION`'s Mercury row carried no tag of its own and was covered only
by the table-level one. It now states its own reason — that no artefact on disk gives a
ΔHvap for mercury and none can be derived, there being no mercury saturation file. A row
that states why it is unverified is more specific than a row sheltering under a blanket.

The count moved twice and in both directions, which is worth keeping visible:

- **198 → 200.** `GAS_MOLECULAR_PARAMS` entered this section for Xenon's σ and Ammonia's
  ε/k (§1). A table that *gains* a residual by declining to write a value it cannot read
  is in better condition than one that quietly writes it.
- **200 → 180.** 18 rows were sourced to PubChem and Tungsten was corrected, each moving
  from `[UNVERIFIED]` to `[ON-DISK]`. A further 24 rows stay here but no longer say "no
  on-disk source, and none was searched for": the fuel and oil rows carry a measured
  n-alkane bracket, and the wood rows carry FPL's actual species range (§6 item 3).

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

> **Seven of these rows are gone.** `269fc76` deleted them, owner-directed: `Cork`,
> `Cork Board`, `Bamboo`, `Tellurium Mercury`, `Tungsten Hexafluoride` and
> `Tetrabromoethane` from mechanical, and `R-410A` from chemical. The table below is kept
> as the record of *why* each was a defect and what was tried before deleting it — that is
> the evidence for the deletion, not a list of open questions. Each had been deferred
> twice on the ground that removing a row shifts every key after it and so invalidates the
> published pool; rebuilding the pool removed that cost. **The blocker was never the edit,
> it was the consequence.** `Acetylene Tetrabromide` was kept over `Tetrabromoethane` as
> the unambiguous name for the same substance. Two entries below are **not** deletions and
> remain open: the annotation fork, and Mercury's ΔHvap.

| row | branch | finding |
|---|---|---|
| `R-410A` | chemical | a blend, absent from the **39** pure fluids on disk (36 before this phase added nonane, decane and dodecane) — **and so are both its components**, R-32 and R-125, so it cannot be reconstructed either. **NOT substituted, for two reasons:** R-22 is the nearest on-disk refrigerant (298.15 K: P 1.0439 MPa, v_f 0.00083987, v_g 0.022608) and this row is **+9.78%** from its v_f and **−33.21%** from its v_g, so substituting would replace the row's *values* rather than relabel it; and `REAL_FLUID_DATA` **already carries an R-22 row**, so it would put one fluid in the table twice under two keys — the very defect recorded one row below. The committed volumes are unsourced: nothing on disk states or contradicts them |
| `REAL_FLUID_DATA`, the annotation copy | chemical | **NEW — a stale fork, flagged and deliberately not synced.** `templates_annotation/annotation_app/data/templates/branches/chemical_engineering/constants.py` carries the same table with **pre-C3 values**: R-12 `v_g 0.0268` where the main file has 0.0272, and R-22 `v_f 0.000845` where the main file has 0.000840 — both corrected during C3 against the NIST TSVs. Verified by reading that file directly. **Not synced on purpose:** that directory backs the annotation pilot, and a frozen snapshot may be intentional so annotations stay reproducible against what annotators actually saw. Syncing it silently could invalidate the pilot's provenance, so it is the repo owner's call. Only `constants.py` is forked — `fetch_references.py` is not duplicated, so `NIST_FLUIDS` has a single authoritative definition |
| `Tellurium Mercury` | mechanical | WebBook returns "Name Not Found"; whether it means mercury telluride is not established |
| `Tungsten Hexafluoride` | mechanical | listed as a manometer liquid; the on-disk Antoine parameters put P = 1.074 bar at 290.4 K, so it is a **gas** at room temperature |
| `2024-T4` | mechanical | **CLOSED — and this entry was wrong on both counts.** p.374's text layer parses perfectly well; it prints "See Table 3.2.3.0(d)" in all four elastic cells, so there is no value on that page to read and the reader is right to refuse. p.373 (Table 3.2.3.0(b1)) parses cleanly and prints G 4.0 ×10³ ksi in all four columns, and following the deferral finds Table 3.2.3.0(d) on p.376, captioned "…Sheet and Plate, **All Tempers**" — which is what licenses applying it to T4, the question the row's own tag left open. `SHEAR_MODULUS_VALUES['Aluminum 2024-T4']` was 28.0 GPa against 27.58 (+1.5%, outside the 1.25% half-unit bound); corrected to 27.6 and cited to p.373, which is machine-readable where p.376 is not |
| `Mercury` ΔHvap | chemical | the WebBook phase-change page has no vaporization section |
| `Cork`, `Cork Board`, `Bamboo` | mechanical | **not wood** — a bark tissue and a grass; the wood reason never applied (H-mech F4). **Now searched exhaustively and confirmed absent:** a full 546-page scan of FPL-GTR-282 for *cork*, *bamboo*, *Bambusa*, *Quercus suber* returns three hits and all three are false positives — "cork elm" (*Ulmus thomasii*, p.34), "Corkbark fir" (*A. lasiocarpa*, p.378) and "bamboo fiber" in a wood-plastic-composite narrative with no property value (p.302); *Bambusa* and *Quercus suber* return zero. No PubChem record (a tissue is not a compound with a CID — the polymer failure mode again), and nothing under any of the 13 `docs/references/` directories, all of which are structurally unrelated to a bark tissue or a grass. **No substitute value was read, so none is proposed.** These rows need a materials handbook |
| `Tetrabromoethane` **and** `Acetylene Tetrabromide` | mechanical | **NEW — one substance under two names.** `MANOMETER_FLUIDS` carries both, each at **2960**, and acetylene tetrabromide *is* 1,1,2,2-tetrabromoethane. The identical value is strong evidence the two rows were meant as one substance. Stated with its limit: **no artefact on disk distinguishes them**, and "tetrabromoethane" could in principle name the 1,1,1,2 isomer, so this is a duplicate row to resolve rather than an identification to act on. A template drawing `random.choice` over this table draws the same fluid twice under two labels |

---

## 4. Tables that look sourced and are not

| table | finding |
|---|---|
| chemical `GAS_MOLECULAR_PARAMS` | §1 row 1 — **resolved.** 16 of 17 rows re-sourced to Svehla, carrying 31 per-row `[ON-DISK]` citations that the resolver re-reads from the 1962 scan through `tests/constants_integrity/svehla.py`. Xenon's σ and Ammonia's ε/k stay `[UNVERIFIED]`, each with its reason |
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

1. ~~The twelve mechanical alloy defects (§1 rows 4–15) — correct, or restate the grades.~~
   **Answered: corrected**, to the MIL-HDBK-5J values this repo already ships. P6 measured.
2. ~~`GAS_MOLECULAR_PARAMS`' σ column — re-source from Svehla (moves items; μ ∝ σ⁻²).~~
   **Answered: re-sourced.** The emitted items then confirmed the argument independently —
   butane's viscosity moved 1.537e-05 → 2.094e-05, and (5.47/4.687)² = 1.362 predicts
   exactly that ratio from the constants alone.
3. ~~Six mechanical tables carry 172 unchecked rows — acquire sources, or accept.~~
   **Partly answered: 172 → 148, and the rest are measured rather than merely unchecked.**

   **The PubChem route is now exhausted, not merely paused.** A final sweep fetched the
   remaining rows that are single named substances and tagged **none** of them. Quartz
   returned the wrong *form* — PubChem's silicon dioxide record gives 2200/2300/2334/2600
   (amorphous silica, fumed silica, silica gel), none of them crystalline α-quartz's 2650 —
   and all six polymers returned `PUGREST.NotFound`, because a polymer is not a single
   compound with a CID, which is the same reason its density depends on grade rather than
   on a formula. Those rows need a **materials handbook**, and no further fetching will
   reach them.

   | what | rows | |
   |---|---:|---|
   | sourced to PubChem, cited to one named entry each | **18** | `[UNVERIFIED]` → `[ON-DISK]` |
   | Tungsten corrected | 1 | 19600 → 19300, five sources at 19.3 |
   | fuel and oil rows given a measured n-alkane bracket | 19 | still `[UNVERIFIED]`, and now say why |
   | wood rows given FPL's actual species range | 5 | still `[UNVERIFIED]`, and now say why |

   The remaining 142 are not one decision taken 142 times. **Counted by the reason each
   row states** (tag-starting lines, sums to 142; the 148 column is the C3-era figure,
   before `269fc76` deleted five of these rows and closed one):

   | rows | then | reason |
   |---:|---:|---|
   | 39 | 39 | no source on disk *for this material* — the bulk-material rows: granite, marble, brick, concrete, glass, pumice, asphalt, coal, paper, leather, bone |
   | 19 | 19 | a mixture or commercial grade with no single composition |
   | 15 | 16 | no source on disk |
   | 12 | 12 | no MIL-HDBK-5J design table is this material, and no other source |
   | 12 | 12 | a mixture, solution or commercial product with no single composition |
   | 10 | 10 | the table states no temperature or pressure, so no row of an artefact *is* the value |
   | 8 | 8 | **searched, and found unreachable by a chemical database** (PubChem: no CID for a polymer) |
   | 5 | 5 | the row names a material FPL indexes by **species** and states no moisture condition |
   | 4 | 4 | no MIL-HDBK-5J design table under the names searched |
   | — | 3 | ~~**not wood** — cork, cork board, bamboo~~ — **rows deleted** (§3) |
   | 3 | 3 | no MIL-HDBK-5J design table (its copper-base tables listed) |
   | 2 | 2 | no temperature or pressure stated, for a cryogen |
   | — | 2 | ~~**one substance under two names**~~ — `Tetrabromoethane` **deleted**, `Acetylene Tetrabromide` kept as the unambiguous name (§3) |
   | 13 | 13 | singletons, each naming its own specific miss. The membership changed with the deletions: `Tellurium Mercury` and WF₆'s phase left, and two rows that now state their own reason arrived — `Acetylene Tetrabromide`'s unsourced 2960, and the row whose tag reads *"a source IS on disk, and it disagrees"*. The rest stand: Xylene, the unnamed aluminium alloy, 1100/3003, quartz's wrong form, CP titanium's unprinted µ, isopropanol, ethylene glycol, chloroform, bromine, ASTM A36, air |

   **The shape of that table is the answer to "why are so many still unverified".** Roughly
   86 rows (39 + 19 + 16 + 12) are bulk materials and mixtures that no chemical database
   indexes — no pure substance *is* the row, and the fuel acquisition proved it:
   nonane/decane/dodecane bracket 718–749 kg/m³ against Kerosene 810, so re-pointing would
   change what the row names. 19 more need a **mechanical-properties handbook**, which
   PubChem cannot supply at any effort. And **17 rows now say "searched, and here is what
   the search showed"** rather than "no source" — which is the difference between a row
   nobody has looked at and a row that has been looked at and found genuinely unreachable.

   What changed is that "no on-disk source, and none was searched for" has become, for 43
   of them, a statement of what *was* searched and what it showed.
4. ~~`COMMON_LIQUIDS`' conditions — fetch a 298.15 K grid and settle §5's open diagnosis.~~
   **Answered: the grid was fetched, and the answer is neither reading offered.**

   The diagnosis could not be tested because the 1-atm isobar steps 10 K from 253.15 and
   so carries no row at 298.15 K — a structural gap, not an oversight. Five grids landing
   exactly on 298.15 K were acquired (the four organics plus **water as the control**).

   | viscosity | at 20 °C | at 25 °C |
   |---|---:|---:|
   | methanol | −7.06% | **+0.05%** |
   | benzene | −7.16% | **−0.33%** |
   | n-hexane | −6.12% | **−1.33%** |
   | toluene | −4.62% | **+1.41%** |
   | **water (control)** | **+0.04%** | +12.58% |

   The four organic **viscosities fit 25 °C**; **water does not**, and the **density
   column fits 20 °C throughout**. So the column is **MIXED** — the organic viscosities
   are at 25 °C while water and every density are at 20 °C.

   The review was right about the organics and wrong about the table. My own reading
   ("a source difference, not scatter") stays **withdrawn**: nothing here tests it, and it
   is not revived by this outcome.

   This is worse than either proposal and different in kind. It **cannot** be repaired by
   editing `@domain`, because no single temperature is true of the table; `@domain` keeps
   293.15 K, which is true of the densities and of water.

   **CORRECTED.** The four organic viscosities were re-sourced to the 293.15 K the table
   declares — methanol 0.000544 → 0.000585, benzene 0.000601 → 0.000647, toluene 0.000560
   → 0.000587, n-hexane 0.000294 → 0.000313 — each now carrying an `[ON-DISK]` tag the
   resolver checks. The densities were not touched.

   P6: **5 templates moved**, q 90 / ans 122 / sol 123, errors 0 → 0. One of them is the
   first **hidden-constant** case of the whole phase: `newtons_law_shear_stress` moved
   **0 questions and 32 answers** — the viscosity is used in the arithmetic and never
   printed in the question, so those instances need **re-scoring but not re-inference**,
   the exact inverse of tranche 1's `statically_indeterminate`. That is the more dangerous
   direction: the question a model saw is byte-identical while the gold answer moved, and
   nothing in the text signals it.
5. ~~**`tol=` has no rule for its SIZE.**~~ **Answered: the rule is built and enforced.**
   D-074 licensed *when* a tolerance may be used and never *how large*, so one could be
   fitted to its own residual and never fail — a check that cannot fail is not a check.
   Measured over **all 17** `tol=` tags (the count in the earlier note, 10, predated the
   mechanical corrections), not only the four whose own prose confessed it:

   | basis | n | |
   |---|---:|---|
   | half a unit in the artefact's last printed digit | **13** | recomputable, and every one matches |
   | its own residual, rounded up | **4** | circular |

   `tol=` now requires `basis=half-unit` — the resolver recomputes half a unit from the
   artefact's printed digits and the stated value must match — or `basis=condition`, where
   the resolver **cannot** size it, counts it, and says so in its disclosure line. A
   `tol=` with no basis **fails**: all 17 failed until each declared one.

   The four are not alike, and the rule respects that. Three cite **dispersion formulas**,
   which print no last digit at all (half a unit of one is ~1e-15), so `half-unit` is
   impossible for them by construction. **`FLUID_DENSITIES['Liquid Propane']` is the sharp
   one**: its artefact *does* print digits, half a unit is **0.0010%**, and the tag states
   **0.1300%**. It is *not* re-sized to 0.0010% and left to fail — half a unit of the
   row's own last digit (0.5/493 = 0.101%) would not admit the value either, so the
   warrant is a pressure the row never states ("At 25 °C, under pressure", while the file
   is saturated liquid at 0.95207 MPa). **A tolerance cannot be sized from a condition that
   is not named.**

   `MEDIA_VELOCITIES['Glycerine']` is labelled and **not thereby fixed**: its 0.12% is
   absorbing a *composition* ambiguity (Gupta states no conditions; glycerol is
   hygroscopic), which its own note already calls the weakest warrant in the branch. The
   label makes it countable, not sound.
6. ~~**16 `[DERIVED]` tags are never recomputed** (G F-4) — the `UNEXECUTED` class is
   unbuilt.~~ **Answered: built, and eight of them now run on every resolver pass.**
   (15 in the earlier note; the civil kinematic-viscosity reclass made it 16.)

   The single number was folding together three different situations — the resolver
   *could* recompute it and doesn't; something else recomputes it and the resolver doesn't
   know; **nothing** recomputes it. Only the third is UNEXECUTED in any useful sense, and
   it is the one a reader needs to find. `[DERIVED]` now takes `recompute=<name>`, naming
   an entry in a registry the resolver **calls**, or `by=<module>`, naming a checker that
   does it elsewhere:

   > `[DERIVED]: 8 RECOMPUTED here, 1 recomputed by a named checker, 7 UNEXECUTED`

   | recomputation | covers |
   |---|---|
   | `nu_from_tsv` | ν = μ/ρ, both operands one row of one NIST isobar |
   | `g_ft_from_codata` | gₙ / 0.3048, CODATA and the SP 811 foot |
   | `normal_quantile` | all 6 `Z_QUANTILES` against `NormalDist().inv_cdf` at 4 dp |
   | `atom_balance` | all 15 reactions balanced from the species keys themselves |
   | `monatomic_cp` | Cp/R = 2.5 exactly with B = C = D = 0, three rows |

   These had been asserted since C1 and tested by nothing. A wrong one now **fails**.

   `CONTROL_CHART_FACTORS` declares `by=control_chart_factors`, which the resolver imports
   to prove the reference is live and does **not** re-run — that module already derives all
   384 cells by quadrature every run, and doing it twice would cost the work to learn
   nothing.

   **Seven remain UNEXECUTED deliberately**, and that is the finding rather than an
   omission: the four `CP_PARAMS_COMBUSTION` rows are least-squares refits of NIST Shomate
   over 400 points, `CP_PARAMS['Air(g)']` is a mixture with no NIST entry, and the two
   `VALID_T_MAX` constants declare the interval that was fitted. Re-running a fit inside
   the resolver would **re-derive** the constant rather than check it — the fit would agree
   with itself by construction, which is the same circularity item 5 was built to refuse.
7. `CP_PARAMS` is integrated from **281.67 K**, below the 298 K floor its fit declares —
   registered in `domain_findings.txt`; raise the draw, refit, or restate the floor.

## 7. Registered elsewhere

| register | holds |
|---|---|
| `tests/constants_integrity/domain_findings.txt` | **1** consumer-domain excursion: `MEDIA_VELOCITIES` optical indices used for 50–500 MHz radio. The second, `CP_PARAMS` integrated below its fit floor, was closed in corrections tranche 2 by raising the liquid branch's draw to `uniform(298.15, 318.0)`; the ratchet fails a listed excursion that no longer occurs, so the suite itself required the line's removal |
| `phaseC3_literal_copies.md` | C3.8 — declared `@copied-in` copies and `_SCS_COMBOS` |
| `phaseC3_inline_windows.md` | C3.9 — named-entity facts written as template literals |

## 8. What C3 did not do

**As C3 shipped**, nothing in §1 was corrected. C3 changed comments, not values, with three
exceptions, all recorded as P6 events and all measured: `REAL_FLUID_DATA` re-derived from
NIST, `CONTROL_CHART_FACTORS`' 29 + 2 corrected last digits, and two template edits. The
item-pool consequence is measured in `phaseC3_item_pool_impact.md`: **2 of 150 templates
moved**.

**That is no longer this file's final state**, and the sentence above is kept in the past
tense rather than deleted, for the reason §5 gives: a register that quietly absorbs its own
corrections is worth no more than the claims it corrected. All 19 of §1 are now corrected,
plus a twentieth found in the process (§3, 2024-T4), and the consequence is measured
separately in `phaseC3_corrections_item_pool_impact.md`: **13 of 150 templates moved**.
What remains open is §2's **170** unchecked rows and §6's items 3–7. (177 when that sentence
was written; `269fc76` then deleted the seven substance-defect rows of §3, owner-directed.)
