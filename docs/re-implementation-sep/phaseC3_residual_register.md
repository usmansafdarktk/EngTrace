# Phase C3.5 — residual register

**What this is.** Spec §C3.5: every constant C3 could not resolve to an artefact on disk,
or resolved and found wrong, with the decision it needs from the repo owner.
**Listing is not acceptance.** Every line below is an open defect or an open question,
recorded so it cannot be silent and so the next phase starts from a count rather than a
search.

**The tags are the authority, not this file.** Each entry corresponds to a
`[KNOWN-DEFECTIVE]` or `[UNVERIFIED]` tag in a branch `constants.py`, which the resolver
checks on every run (R6: a tag of either class must state its reason). Re-derive every
count below rather than trusting these lines:

```
grep -c "\[KNOWN-DEFECTIVE\]" data/templates/branches/*/constants.py   # 19
grep -c "\[UNVERIFIED\]"      data/templates/branches/*/constants.py   # 200
python -m tests.constants_integrity.test_citations_resolve
```

**At the C3.7 commits** (branch `redesign/phaseC3-remaining-tables`): 219 tagged
residuals — **19 [KNOWN-DEFECTIVE]** and **200 [UNVERIFIED]** — over 450 tags, of which
152 resolve to an artefact, 19 are locator-only and 279 state a class. No branch carries
a LEGACY tag; no numeric table is missing `@domain`.

| branch | tags | resolved | KNOWN-DEFECTIVE | UNVERIFIED |
|---|---|---|---|---|
| chemical | 138 | 105 | 0 | 11 |
| mechanical | 206 | 19 | 14 | 173 |
| civil | 34 | 10 | 3 | 8 |
| electrical | 34 | 16 | 2 | 4 |
| industrial | 38 | 2 | 0 | 4 |

---

## 1. Measured wrong — the 19 [KNOWN-DEFECTIVE]

Each was compared with a named artefact **at the table's own declared conditions** and is
not a rounding of it. **No value was changed**: correcting a constant moves the item pool,
which is a P6 event and the owner's call (D-033).

### Mechanical — 14, against MIL-HDBK-5J design tables and NIST

| constant | table | handbook / NIST | error |
|---|---|---|---|
| AISI 301 stainless `E_ksi` | 27500 | p.277 Table 2.7.1.0(b), E:L 29.0e3 ksi | -5.2% |
| AISI 301 stainless `E_GPa` | 190 | = 199.9 GPa | -5.0% |
| AISI 301 stainless `nu` | 0.30 | mu 0.27 | +11.1% |
| 6061 aluminium `E_ksi` | 10000 | p.566 Table 3.6.2.0(b1), E 9.9e3 ksi | +1.0% |
| 6061 aluminium `E_GPa` | 68.9 | = 68.26 GPa | +0.9% |
| CP titanium `E_ksi` | 16800 | p.899 Table 5.2.1.0(b), E 15.5e3 ksi | +8.4% |
| CP titanium `E_GPa` | 116 | = 106.9 GPa | +8.5% |
| Ti-6Al-4V `E_ksi` | 16500 | p.945 Table 5.4.1.0(b), E 16.0e3 ksi | +3.1% |
| Ti-6Al-4V `E_GPa` | 114 | = 110.3 GPa | +3.4% |
| Ti-6Al-4V `nu` | 0.34 | mu 0.31 | +9.7% |
| Ti-6Al-4V `G_GPa` | 41.4 | G 6.2e3 ksi = 42.75 GPa | -3.2% |
| AZ31B magnesium `nu` | 0.29 | p.841 Table 4.2.1.0(b), mu 0.35 | -17.1% |
| benzene density | 876 | NIST 878.92 at 293.15 K | -0.33% |
| R-134a density | 1206 | NIST 1206.7 at 298.15 K | -0.06% |

The elastic-property cluster is **one decision, not twelve**: E and nu for four alloy
families disagree with the handbook this repo ships. AZ31B's `nu` at -17.1% is the largest
measured error in the corpus.

### Civil — 3, all water at 20 C (C3.7)

| constant | table | artefact | error |
|---|---|---|---|
| `UNIT_WEIGHT_WATER_KN_M3` | 9.81 | rho*g_n = 9.7891 kN/m^3 | +0.21% |
| `UNIT_WEIGHT_WATER_PCF` | 62.4 | 62.316 lb/ft^3 | +0.13% |
| `WATER_KINEMATIC_VISCOSITY_M2_S` | 1.004e-6 | mu/rho = 1.003396e-6 | +0.060% |

(rho = 998.21 kg/m^3 and mu = 0.0010016 Pa*s at 293.15 K from the NIST isobar; g_n from
CODATA; the lb/ft^3 factor from SP 811 p.66.)

One defect with one cause: 9.81 kN/m^3 implies rho = 1000.34 kg/m^3 and 62.4 lbf/ft^3
implies 999.55 — **water near 4 C**, not the "~15-20 C" the file claims. The viscosity
follows from a 4-digit mu of 1.002e-3 Pa*s rather than the 5-digit mu NIST tabulates.
**Decision:** correct all three to the declared conditions, or restate the conditions to
the ones the values describe.

### Electrical — 2

| constant | table | artefact |
|---|---|---|
| `Helium` refractive index | 1.000036 | 1.0000349 (Ermolov) and 1.0000349 (Mansfield), both at 273.15 K / 101325 Pa — the group comment's own conditions |
| `Gallium Phosphide (GaP)` | 3.5 | 3.36-3.41 from all five on-disk datasets covering 589 nm |

Helium is +1.1e-06 on n — below every display in the corpus — but **+3.2% on the
refractivity n-1**, which is the quantity that does work. GaP is +2.7% to +4.1% high and
every dataset rounds to 3.4 at the row's own 1-dp precision.

---

## 2. No artefact on disk — the 200 [UNVERIFIED]

Not known to be wrong; known to be **unchecked**. Counted by the phrase each tag uses
(`grep -c` on the strings below), so the classes are measured, not estimated:

| class | count | branch |
|---|---|---|
| "no source on disk for this material" / "no source on disk. Residual register" | **95** | mechanical |
| "no MIL-HDBK-5J design table is this material" | 19 | mechanical |
| "a mixture, solution or commercial product with no single composition" | 19 | mechanical |
| "not in the NIST fluid database on disk" (36 pure fluids) | 15 | mechanical 12, chemical 3 |
| "the table states no temperature or pressure, so no row of an artefact is THE row" | 12 | mechanical |
| "names no species and no moisture content" (wood) | 8 | mechanical |
| "…is not on disk" (a named standard) | 8 | civil 4, industrial 3, chemical 1 |
| not covered by any phrase above | 24 | mechanical 8, chemical 7, electrical 4, civil 4, industrial 1 |

The 95-entry class is the single largest fact in this register: **almost half of every
tagged residual in the corpus is one mechanical table's worth of materials with no
primary source at all.** It is one decision, taken once, not 95.

**The named standards (8).** AISC Steel Construction Manual and ACI 318-19 and ASCE 7-22
(civil), IEC 60063 and MIL-A-8625 (industrial), US Standard Atmosphere 1976 (chemical).
Acquire, or accept as unverifiable. Note that MIL-STD-105E, long cited as if it were in
the copyrighted-books tree, **is** a public on-disk artefact and now resolves.

**Electrical's 4**, which no class covers, each needing its own decision:

| row | why |
|---|---|
| `Teflon (PTFE)` | no PTFE dataset in the archive (searched: PTFE, teflon, tetrafluoroethylene, (C2F4)n) |
| `Polyethylene` | both (C2H4)n datasets have no rows bracketing 589 nm |
| `Glass (amorphous semiconductor)` | names a class, not a material: 18 amorphous datasets span n = 1.74 to 3.88 at 589 nm |
| human tissue / muscle | 7.14 is a **microwave** value (sqrt of the relative permittivity ~51 its own comment states), not an optical index — a different quantity from the row's |

---

## 3. Tables that mix sources

Resolved far enough to show the table is not what it says it is. These matter most,
because each *looks* sourced.

| table | finding | decision |
|---|---|---|
| chemical `GAS_MOLECULAR_PARAMS` | Svehla NASA TR R-132 Table I(a) **is** on disk and is the classic source — but this is not it. Parsed from the scan over 15 of 17 rows: eps/k agrees for 11 of 14, sigma for only 3 (CO2, Ar, He). n-butane is 5.47 against Svehla 4.687. Air, N2 and H2 differ in **both** columns, so they come from a third compilation | re-source sigma, or name the compilation the table follows |
| chemical `COMMON_LIQUIDS` | every organic liquid with a NIST isobar on disk has its **viscosity low** at the conditions the table declares — methanol -7.06%, benzene -7.16%, toluene -4.62%, n-hexane -6.12% — while densities agree to 0.001-0.7%. One-sided across four independent species | re-source the viscosity column |
| chemical `COMMON_GASES` | 10 (row, field) pairs equal NIST; the rest sit 0.1-6% off (worst SF6 viscosity +6.06%, butane +3.45%) | re-source or accept |
| civil `LIVE_LOADS_KPA` | not a conversion of `LIVE_LOADS_PSF`: 4 of 5 rows are the **soft** conversion 1 psf = 0.048 kPa exactly (+0.250% against SP 811); "corridor" is neither — 4.79, where soft gives 4.80 and SP 811 gives 4.788 | it is ASCE's own SI column; acquire ASCE 7-22 or accept |
| chemical `SUBSTANCES_FOR_VAPORIZATION` | **nearly** derivable: dHvap = H(v) - H(l) at 1 atm from the NIST saturation tables reproduces 7 of the 12 unsourced rows to within 0.04-1.45%, but only methanol, benzene, ammonia and oxygen are roundings of it — so the table is UNVERIFIED, not DERIVED | decide whether ~1% agreement is a derivation |
| chemical `CRITICAL_PROPERTIES` | 3 of 27 rows resolve to a WebBook page (ethanol 513.9 K, acetone 508.2 K, p-xylene 616.2 K, exact); 24 have no artefact | fetch 24 species pages, or accept |

---

## 4. Substance defects — the row names the wrong thing (D-033)

Not a wrong number: a row that names something the repo cannot source **because the
substance itself is the problem**. Replacing a substance is explicitly the owner's call.

| row | branch | finding |
|---|---|---|
| `R-410A` | chemical | a blend, absent from the 36 pure fluids on disk; needs a mixture source or a replacement species |
| `Tellurium Mercury` | mechanical | the WebBook name search returns "Name Not Found"; whether the row means **mercury telluride** is not established |
| `Tungsten Hexafluoride` | mechanical | listed as a manometer liquid, but [DERIVED] from the on-disk WebBook Antoine parameters (201.5-290.4 K, A = 4.55569, B = 1021.208, C = -64.7): at 290.4 K, P = 1.074 bar, already above 1 atm and rising — it is a **gas** at room temperature |
| `2024-T4` | mechanical | MIL-HDBK-5J **does** carry the table (p.374); its text layer does not parse. A page-image or OCR pass would resolve it |
| `Mercury` dHvap | chemical | 59.11 kJ/mol; the WebBook phase-change page on disk has no vaporization section at all |

---

## 5. Data integrity found while sourcing

- **civil `PERMEABILITY_RANGES_CM_S`** — four of five rows are Das Table 7.1 endpoint for
  endpoint. The fifth is not: Das gives clay as "<0.000001", an upper bound with **no
  floor**, where the table writes `(1e-8, 1e-6)`. The 1e-8 is the authoring round's.
- **civil `SPECIFIC_GRAVITY_RANGES`** — was cited to Das, but §2.6 gives only "most of the
  values fall within a range of 2.6 to 2.9" and its Table 2.4 is *Specific Gravity of
  Common **Minerals***, not soil types. Retagged `[POLICY: sampling-only]`.

## 6. Referrals that are not defects

`grep -c "C3.5"` finds **216** mentions across four files, against 219 KD/UNVERIFIED
tags — the two sets are not the same. Some rows **resolve** to an artefact and still refer
a question, because the row names a family and the implementer picked a member:

- electrical `Crown Glass (typical)` resolves against Schott **N-BK7**, `Flint Glass
  (dense)` against **N-SF2** — chosen by the implementer, not by the row.
- mechanical density rows that name no grade resolve against AISI 1025 and commercially
  pure titanium, each inside half a unit in the handbook's last digit, with the grade
  still the implementer's choice.
- `Sapphire` resolves against the **ordinary** ray; sapphire is birefringent and the row
  names no ray (the extraordinary ray gives 1.76, not 1.77).

## 7. Referred elsewhere, not repeated here

| register | holds |
|---|---|
| `tests/constants_integrity/domain_findings.txt` | consumers measured outside a table's `@domain` — currently 1: `MEDIA_VELOCITIES` optical indices (0.589 um) used for 50-500 MHz radio waves |
| `phaseC3_literal_copies.md` | C3.8 — 4 declared `@copied-in` copies, plus `_SCS_COMBOS`: ten hard-coded curve-number triples a template draws from while reading the table only in an assert |
| `phaseC3_inline_windows.md` | C3.9 — named-entity facts written as template literals, including chart-pair subgroup windows that **disagree** with `XBAR_R/S_SUBGROUP_N`, and nine soil unit-weight windows backed by no table |

## 8. What C3 did not do

Nothing in this register was fixed. C3 changed comments, not values — every patch
re-executed its module before and after and asserted the namespace was repr-identical
with unchanged dict key order. The two exceptions are recorded as P6 events and were made
**against artefacts**: `REAL_FLUID_DATA` re-derived from NIST, and `CONTROL_CHART_FACTORS`'
29 corrected last digits.
