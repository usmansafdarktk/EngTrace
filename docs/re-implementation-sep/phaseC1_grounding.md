# Phase C1.1 — which on-disk artefact could ground each table

**Status: a map of candidate routes, not a verification.** C1.1 asks, for every
table C1.3 classifies as needing a citation, *which on-disk artefact grounds it -
or that none does*. This file answers at table level. **No value in it has been
checked against its artefact**; that is C3's work, row by row, under the resolver.

Every coverage statement below cites the probe that established it (§1). A
statement no probe established is marked **unprobed**, and a belief held from
general knowledge rather than from a file is marked **unverified** - the grounded-
retrieval rule applies to what a planning document claims, not only to values.

Classification from [`phaseC1_census.md`](phaseC1_census.md) (36 CITATION tables).
Artefact paths are relative to `docs/references/`.

## 1. What was probed (2026-09-11)

| # | probe | result |
|---|---|---|
| p1 | column header of `nist_fluid_properties/nitrogen_C7727379_saturation.tsv` | liquid **and** vapour density, specific volume, enthalpy, Cp, viscosity, thermal conductivity; last row is the critical point (126.19 K, 3.3958 MPa, 313.30 kg/m³) |
| p2 | `T<sub>c</sub>`, `P<sub>c</sub>`, `vap</sub>H` in `nist_webbook_species/*_phase_change.html` | Tc on 11 of 12 pages (not radon); Pc on 10 (not mercury, radon); ΔvapH on 9 (not chlorine, mercury, radon) |
| p3 | `C<sub>p,liquid</sub>` in `nist_webbook_species/*_condensed_phase.html` | 8 of 12 (not acetylene, chlorine, mercury, radon) |
| p4 | PyMuPDF text layer, 40 sampled pages per PDF | MIL-HDBK-5J 1733 pp, median 993 chars/page; Wood Handbook 546 pp, median 4586; NASA TR R-132 122 pp - cover is OCR noise, table headings (`TABLE III`) read; NAVFAC DM-7.01 has pages with **no** text layer; MIL-STD-105E 73 pp with text |
| p5 | MIL-HDBK-5J text, regex per page | "E, 10³ ksi" 224 pp, "G, 10³ ksi" 225, "lb/in.³" 222; "ALUMINUM ALLOYS" 41; "MAGNESIUM ALLOYS" 15; "TITANIUM" 73; "HEAT-RESISTANT ALLOYS" 6; "nickel-base" 17; "COPPER" 30 (a word, not established as a chapter); "BERYLLIUM" 26; "tungsten" 11; **"cast iron" 0**; a steel-chapter regex matched nothing, so **steel coverage is unprobed** |
| p6 | Wood Handbook text | "specific gravity" 78 pp; "specific heat" 1 (p. 116); "teak" 16; "ebony" 2; "balsa" 12 |
| p7 | NAVFAC text | DM-7.02 "bearing capacity factor" 7 pp (first 12, 121, 149); DM-7.01 "compression index" 5 pp, "0.009" 1 p (243) |
| p8 | folder names in the refractiveindex.info archive | `main/` holds H2O, SiO2, Al2O3, C, GaP, CS2, CO2, He among others |
| p9 | `MANIFEST.json` | 120 sources present: 13 binary (gitignored), 107 text (committed) |

**No NIST saturation grid was checked for a specific temperature row** (77 K, 90 K,
a normal boiling point). The acquisition's lesson is that a nearby row is not the
row, and that check belongs to the transcription, per row.

## 2. Chemical (10)

| table | candidate artefacts | known gaps |
|---|---|---|
| `CRITICAL_PROPERTIES` | Tc, Pc, ρc: last row of `nist_fluid_properties/*_saturation.tsv` (p1) for the 23 of 27 species the fluid set carries; Tc, Pc: `nist_webbook_species/` (p2) for p-xylene, ethanol, acetone, chlorine. **Vc = M/ρc and Zc = PcVc/RTc are `[DERIVED]`; ω = −log10(Psat/Pc at Tr = 0.7) − 1 is `[DERIVED]`** from a saturation curve | **molar masses: no atomic-weights artefact is on disk**; ω for the four WebBook-only species has no saturation curve on disk |
| `SUBSTANCES_FOR_HEATING` | solids: `nist_janaf/` Fe, Cu, Al, Pb, W, Si, C; liquid Hg: `nist_janaf/Hg_ref.txt`; fluids: NIST isobar Cp (p1's column set) for water, methanol, N₂, O₂, H₂, He, Ar, CO₂, CH₄, NH₃, steam; liquids: condensed-phase pages (p3); wood: Wood Handbook "specific heat" (p6, one page) | per-gram Cp needs molar masses (above); **Au, Ag** (NIST-JANAF does not carry them); glass, concrete, polyethylene, olive oil, engine oil, sulfuric acid, ice, air, chlorine gas - unprobed or absent |
| `SUBSTANCES_FOR_VAPORIZATION` | ΔHvap = H_v − H_l at 1 atm from `*_saturation.tsv` - **`[DERIVED]`, an interpolation**, since the grid is in T; ΔvapH from WebBook pages (p2) for ethanol, isopropanol, acetone, CCl₄ | mercury (p2: no ΔvapH) |
| `HEATS_OF_FORMATION`, `CP_PARAMS`, `CP_PARAMS_COMBUSTION` | `nist_webbook/shomate_coefficients.json` - **closed by C2**, not reopened | - |
| `COMMON_LIQUIDS` | 293.15 K isobar density + viscosity for water, methanol, benzene, toluene, hexane; saturation curves for liquid N₂ and O₂ | no viscosity on disk for ethanol, IPA, acetone, diethyl ether, glycerol, ethylene glycol, mercury (the WebBook pages are phase-change and condensed-phase thermochemistry only); **mixtures unsearched**: seawater, oils, fuels, honey, milk, blood, syrup, 98 % H₂SO₄ |
| `COMMON_GASES` | 293.15 K isobar density + viscosity for N₂, O₂, CO₂, Ar, He, Ne, Kr, Xe, CH₄, C₂H₆, C₃H₈, C₄H₁₀, H₂, NH₃, SF₆; steam from the water isobar's vapour rows | air and natural gas (mixtures, no fluid file); acetylene and chlorine (no transport data on disk); radon (p2: no data) |
| `GAS_MOLECULAR_PARAMS` | `nasa_tr_r132/` (p4) - **transcribe from page images**, and grep the transcription for scaled-heading literals first (Reviewer H) | per-species coverage unprobed |
| `POWER_LAW_FLUIDS` | **none on disk; unsearched** | all 22 rows |

## 3. Electrical (2)

| table | candidate artefacts | known gaps |
|---|---|---|
| `EPSILON_0` | `codata_2022/allascii.txt` - **done in C1.7** | - |
| `MEDIA_VELOCITIES` | refractive index read at the stated wavelength from `refractiveindex_info/` YAML (p8) for water, fused silica, sapphire, diamond, GaP, CS₂, CO₂, He; ethanol, glycerol, benzene, air unprobed | "Crown Glass (typical)", "Flint Glass (dense)" and "Glass (amorphous semiconductor)" name families, not one material; "Human Body Tissue (muscle, ~3 GHz)" is a microwave permittivity, not an optical index; ice, PTFE, polyethylene, polystyrene unprobed |

## 4. Mechanical (6)

| table | candidate artefacts | known gaps |
|---|---|---|
| `MATERIAL_PROPERTIES`, `SHEAR_MODULUS_VALUES` | `mil_hdbk_5j/` E, G tables (p5) for its aluminium, magnesium, titanium and heat-resistant/nickel-base alloys; **temper and product form are part of the value** (D-032) | steel unprobed; cast iron absent from the text (p5); copper alloys, lead, molybdenum unprobed; polymers, composites, concrete, glass, alumina - no candidate on disk |
| `FLUID_DENSITIES`, `PIPE_FLUIDS`, `MANOMETER_FLUIDS` | 293.15 K isobars (water, methanol, benzene, toluene, hexane, He, CO₂); saturation curves for liquid N₂, O₂, H₂, propane at 25 °C, R-134a at 25 °C | ethanol, IPA, acetone, glycerol, CCl₄, chloroform, xylene, mercury, bromine, acids, oils, fuels, milk, blood, glycol, syrups; "Tungsten Hexafluoride" and "Tellurium Mercury" need an established phase before any density (brief, §Entries to check) |
| `MATERIAL_DENSITIES` | `mil_hdbk_5j/` "lb/in.³" (p5) for the alloys above; `usda_wood_handbook/` specific gravity (p6) for woods including teak, ebony, balsa | pure elements (NIST-JANAF carries no densities), plastics, stone, biological materials |

## 5. Civil (11)

| table | candidate artefacts | known gaps |
|---|---|---|
| `AISC_W_SHAPES` | `civil/aisc_shapes_database_v16.xlsx` - the resolver needs an xlsx locator (C3.7) | - |
| `MANNINGS_N_CHANNELS` | `civil/fhwa_hds4_highway_hydraulics.pdf` Table B.2 (the branch's own citation) | - |
| `RATIONAL_C` | `civil/fhwa_hec22_urban_drainage.pdf` Table 3-1 (the branch's own citation) | - |
| `SCS_CURVE_NUMBERS` | `civil/nrcs_tr55_urban_hydrology.pdf` Tables 2-2a/b/c (the branch's own citation) | - |
| `TERZAGHI_BEARING_FACTORS`, `TERZAGHI_MODIFIED_FACTORS` | Das (**local-only**); `civil/navfac_dm7_02_foundations.pdf` "bearing capacity factor" pages (p7) as a public cross-check | whether DM-7.02 tabulates the same factors is unprobed |
| `SKEMPTON_CC_COEFF`, `SKEMPTON_CC_OFFSET` | Das (**local-only**); `civil/navfac_dm7_01_soil_mechanics.pdf` p. 243 carries "0.009" and five pages carry "compression index" (p7) | whether p. 243 is Skempton's correlation is unprobed |
| `STEEL_E_KSI`, `STEEL_E_GPA` | the AISC Manual the branch cites is **not on disk**; `mil_hdbk_5j/` as an independent public source, steel unprobed (p5) | - |
| `UNIT_WEIGHT_WATER_KN_M3` | **`[DERIVED]`** from the water isobar density and CODATA's standard acceleration of gravity | the row's "~15-20 C" is not one condition |

## 6. Industrial (7)

| table | candidate artefacts | known gaps |
|---|---|---|
| `MIL_STD_105E_*` (3 tables) | `industrial/mil_std_105e_sampling.pdf`, PDF pp. 18-19 (the branch's own citation) | - |
| `IEC60063_E24_5PCT`, `IEC60063_E12_10PCT`, `RESISTOR_SERIES_BY_TOLERANCE` | **none on disk** | all; whether a public-domain copy of the E-series exists is unsearched |
| `HARD_ANODIZE_THICKNESS` | **none on disk**; the branch cites MIL-A-8625, which was not acquired | whether MIL-A-8625 is publicly released is **unverified** - check its distribution statement before adding it to `fetch_references.py` |

## 7. The other classes, for completeness

| class | tables | warrant route |
|---|---|---|
| DERIVATION | `REACTIONS`, `COMBUSTION_REACTIONS` | element balance, already asserted by the C2.2 suite |
| DERIVATION | `Z_QUANTILES` | `statistics.NormalDist().inv_cdf` |
| DERIVATION | `CONTROL_CHART_FACTORS` | c₄ in closed form (`industrial/nist_sematech_pmc32.htm`); A₂, D₃, D₄ from d₂, d₃ (`industrial/nist_sematech_pmc321.htm`). **d₂ and d₃ are not tabulated in any clone-recoverable artefact** (`README.md` records that pmc321 does not tabulate them; Montgomery, local-only, does), so they need their integral definitions |
| DEFINITION | `C0` | done in C1.7 |
| DEFINITION | `GRAVITY`, `GRAVITY_M_S2`, `GRAVITATIONAL_ACCELERATION` | CODATA's standard acceleration of gravity, at the 3 s.f. the tables use |
| DEFINITION | `SHEWHART_K_SIGMA` | the three-sigma convention - `[BY-DEFINITION]` |
| UNCONSUMED | `REAL_FLUID_DATA` and 31 others | **`REAL_FLUID_DATA` stops being unconsumed at C3.3**, when `template_two_phase_specific_volume` reads it; its route is the saturation curves (p1), and R-410A, a blend, has no fluid file |
