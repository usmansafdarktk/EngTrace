# Phase C3 — Reviewer H: Domain (thermochemistry / materials / fluids), CIVIL branch

**Filed:** 2026-09-12 · **Frozen ref:** `2fe50c72153892d0b1fa67e7ed486ebb36abca0f`
**Roster:** H — Domain, phases C2 and C3 · **Scope:** civil branch only (34 tags)
**Brief:** one mandatory gate task — *independently re-derive at least 8 civil values
directly from the primary artefacts on disk, without reference to the implementer's
transcription, and judge whether each row's provenance CLASS is the honest one.*
40-minute box on the mandatory task; §5 separately boxed to 15 minutes.

> Worked without sight of the implementer's reasoning. Code comments, commit messages,
> `phaseC3_residual_register.md` and `phaseC3_literal_copies.md` were treated as **claims
> under test**. Tag syntax, `@domain` presence and register completeness are Reviewer G's
> and are not judged here; I own whether the NUMBERS are right and whether each row's
> class is honest to a domain reader.

**Local-only material.** Nothing from the Das books is reproduced below. Where a row's
source is one of them I record *whether* my reading agrees, not the book's cells.

---

## 1. Verdict

**PASS WITH FINDINGS** — every leaf I could reach reproduces (196/196 AISC cells exact;
all 12 + 5 Manning's rows, all 16 rational-C rows and all 24 curve-number rows exact from
the public PDFs; all six local-only Das citations carry the claimed table) — but **two of
the three `[KNOWN-DEFECTIVE]` water rows are misclassified: the diagnosis names a cause
that is physically impossible, and the value is a discipline convention stated verbatim on
a page of the very book six other rows in this file already cite.**

---

## 2. Independent re-derivation

Method: PDFs read with `pypdf` page by page; the AISC workbook read with `openpyxl` using
column indices I located myself from the header row (**not** via
`tests/constants_integrity/xlsx_cells.py`, which is itself a claim under test); NIST TSV
and CODATA parsed directly. Scratch scripts under my temp dir; no tracked file modified
except this report.

| # | Value | Artefact and page consulted | My number | Committed | Agree? |
|---|---|---|---|---|---|
| 1 | `WATER_DENSITY_KG_M3` | `water_C7732185_isobar_1atm.tsv` T=293.15, col "Density (kg/m3)" | 998.21 → **998.2** (4 s.f.) | 998.2 | **Yes** |
| 2 | `UNIT_WEIGHT_WATER_KN_M3` | same TSV + `codata_2022/allascii.txt` "standard acceleration of gravity" = 9.806 65 (exact) | ρ·g = **9.78910** kN/m³ | 9.81 | **No**, +0.2135% — *class disputed, see F-1* |
| 3 | `UNIT_WEIGHT_WATER_PCF` | same TSV + SP 811 p.66 `lb/ft³ → kg/m³ = 1.601 846 E+01` | 998.21/16.018463 = **62.3162** | 62.4 | **No**, +0.1344% — *class disputed, see F-1* |
| 4 | `WATER_KINEMATIC_VISCOSITY_M2_S` | same TSV, T=293.15, μ = 0.0010016 Pa·s | μ/ρ = **1.0033961e-6** | 1.004e-6 | **No**, +0.0602% — class defensible, see F-3 |
| 5 | `GRAVITY_M_S2` | `codata_2022/allascii.txt` | 9.806 65 → **9.81** (3 s.f.) | 9.81 | **Yes** |
| 6 | `AISC_W_SHAPES` — **all 196 leaves** | `aisc_shapes_database_v16.xlsx`, "Database v16.0", US block label col C / SI block label col CH | 196 exact matches, **0 mismatches** | — | **Yes** |
| 7 | `AISC_W_SHAPES["W8X24"]["us"]["Ix"]` | same, row `W8X24` | **82.7** in⁴ | 82.7 | **Yes** |
| 8 | `AISC_W_SHAPES["W36X135"]["si"]["Sx"]` | same, row `W920X201` | **7190** ×10³ mm³ | 7190 | **Yes** |
| 9 | `MANNINGS_N_CHANNELS["minor stream, sluggish and weedy"]` | HDS-4 **p.199** (footer B-2), "Table B.2" | **0.050–0.080** | (0.050, 0.080) | **Yes** |
| 10 | `MANNINGS_N_CHANNELS` — other 11 rows | HDS-4 p.199 | all 11 reproduce | — | **Yes** |
| 11 | `MANNINGS_N_CONDUITS["steel pipe"]` | HDS-4 **p.201** (footer B-4), "Table B.3" | **0.009–0.013** | (0.009, 0.013) | **Yes** |
| 12 | `MANNINGS_N_CONDUITS` — other 4 rows | HDS-4 p.201 | all 4 reproduce | — | **Yes** |
| 13 | `RATIONAL_C["streets: concrete"]` | HEC-22 **p.56** (footer 3-6), "Table 3-1" | **0.80–0.95** | (0.80, 0.95) | **Yes** — and it *discriminates*: HDS-4 Table B.1 p.198 gives 0.70–0.95 |
| 14 | `RATIONAL_C["roofs"]` | HEC-22 p.56 | **0.75–0.95** | (0.75, 0.95) | **Yes** — HDS-4 B.1 gives 0.70–0.95 |
| 15 | `RATIONAL_C` — other 14 rows | HEC-22 p.56 | all 14 reproduce | — | **Yes** |
| 16 | `SCS_CURVE_NUMBERS["residential: 1/4 acre lots (38% impervious)"]` | TR-55 **p.14** (footer 2-5), "Table 2-2a" | **61, 75, 83, 87**; impervious 38 | (61, 75, 83, 87) | **Yes** |
| 17 | `SCS_CURVE_NUMBERS` — 16 other 2-2a rows | TR-55 p.14 | all 16 reproduce, impervious percentages included | — | **Yes** |
| 18 | `SCS_CURVE_NUMBERS` — 4 rows from 2-2b | TR-55 **p.15** (footer 2-6) | all 4 reproduce | — | **Yes** |
| 19 | `SCS_CURVE_NUMBERS["woods, good condition"]` | TR-55 **p.16** (footer 2-7), "Table 2-2c" | **30, 55, 70, 77** (the "30" carries footnote 4 in the text layer as `430`) | (30, 55, 70, 77) | **Yes** |
| 20 | `SCS_CURVE_NUMBERS` — 3 other 2-2c rows | TR-55 p.16 | all 3 reproduce | — | **Yes** |
| 21 | `SCS_IA_RATIO` | TR-55 **p.10** — Eq. 2-2 present verbatim | **0.2** | 0.2 | **Yes** (but the prose page label is wrong — F-4) |
| 22 | `@domain: AMC=II, Ia_ratio=0.2` on the CN table | TR-55 pp.14, 15, 16 | each page footnotes "Average runoff condition, and Ia = 0.2S" | — | **Yes** |
| 23 | `PERMEABILITY_RANGES_CM_S` — 4 non-clay rows | Das PGE **p.242**, "Table 7.1" present | endpoint for endpoint | — | **Yes** |
| 24 | `PERMEABILITY_RANGES_CM_S["clay"]` floor `1e-8` | Das PGE p.242 | book gives an **upper bound only, no floor** | (1e-8, 1e-6) | **Claim confirmed** — the floor is not the book's |
| 25 | `DAS_NATURAL_STATE_SOILS` — all 9 rows × 4 fields | Das PGE **p.93**, "Table 3.1" present | all 36 reproduce | — | **Yes** |
| 26 | `FRICTION_ANGLE_RANGES_DEG` — all 6 rows | Das PGE **p.495**, "Table 12.1" present | all 6 reproduce endpoint for endpoint | — | **Yes** |
| 27 | `TERZAGHI_BEARING_FACTORS` — 9 rows × 3 | Das PGE **p.742**, "Table 16.1" present | all 27 reproduce | — | **Yes** |
| 28 | `TERZAGHI_MODIFIED_FACTORS` — 4 rows × 3 | Das PFE **p.160**, caption present | all 12 reproduce | — | **Yes** |
| 29 | `SKEMPTON_CC_COEFF` / `_OFFSET` | Das PGE **p.447**, Eq. 11.39 present | 0.009 and 10 | 0.009 / 10 | **Yes** |
| 30 | `LIVE_LOADS_KPA` — 4 rows | SP 811 p.65 `lbf/ft² → Pa = 4.788 026 E+01` | hard: 2.394013 / 1.915210 / 0.957605; each committed row is **exactly +0.25013%** | 2.40 / 1.92 / 0.96 | **Claim confirmed** — soft 0.048 kPa/psf |
| 31 | `LIVE_LOADS_KPA["corridor (first floor)"]` | same | hard 4.788026, +0.0412% to 4.79; soft would give 4.80 | 4.79 | **Claim confirmed** |

Tag accounting cross-checked by running the repo's own resolver:
`python -m tests.constants_integrity.test_citations_resolve` → `civil_engineering tags 34
resolved 10 locator-only 8 stated 16 … all pass`, matching the brief's 10 / 8 / 3+8.

---

## 3. Findings

### F-1 — `CONFIRMED`. The two unit-weight rows are misclassified: the stated cause is physically impossible, and the value is a sourced convention, not a defect

**Three separate claims fail.**

*(a) The register's causal diagnosis cannot be true.* `phaseC3_residual_register.md` §1
says 9.81 kN/m³ "implies rho = 1000.34 kg/m³ and 62.4 lbf/ft³ implies 999.55 — **water
near 4 C**, not the '~15-20 C' the file claims." No liquid water at 1 atm has a density of
1000.34 kg/m³ at any temperature: water's maximum density is 999.975 kg/m³ at 3.98 °C. The
implied density is *above the physical maximum*, so "water near 4 °C" is not an available
explanation. (999.55 is a real water density, but at ~11 °C, not 4 °C.)

*(b) The actual provenance is a convention, and it is on disk.* Both literals follow from
one rule — ρ taken as **1000 kg/m³ exactly** — which is what geotechnical practice does:

```
1000 * 9.81 / 1000       = 9.81      kN/m^3   exactly
1000 / (0.45359237/0.3048**3) = 62.4280 lb/ft^3  -> 62.4 at 3 s.f.
```

Das & Sobhan, *Principles of Geotechnical Engineering*, **PDF p.92** — one page before
Table 3.1, which this same file transcribes as `DAS_NATURAL_STATE_SOILS` — states in its
own text that the unit weight of water is 9.81 kN/m³, or 62.4 lb/ft³, or 1000 kgf/m³, and
defines γ(kN/m³) = gρ(kg/m³)/1000 with g = 9.81 m/s². The trio 9.81 / 62.4 / 1000 kgf makes
the convention explicit. So this is not an unsourced number: it resolves
`[ON-DISK:LOCAL-ONLY] … @ page=92` on exactly the footing the branch's six other Das rows
already stand on.

*(c) The repo's own rule was applied to a domain the row invented.* SPEC-CHANGE 23 (D-074)
makes `[KNOWN-DEFECTIVE]` fire "when the **conditions match** and the literal is not a
rounding". The conditions here are `@domain: T=288.15..293.15 K` — written by this row
about itself, sourced from nothing. Against a *conventional* γw there is no temperature to
match. The defect, if the row has one, is in the `@domain` line, not in the number.

*The internal inconsistency this creates.* `GRAVITY_M_S2 = 9.81`, **three lines above**, is
tagged `@kind: defined`, `@domain: none (a defined standard value, not a measurement at a
condition)`, `[ON-DISK] … precision=3sf` — the same numeral, reached by the same 3-s.f.
rounding step, accepted without comment. γw = ρg with ρ conventional is the same kind of
object as g itself. The file calls one 9.81 *defined* and the next 9.81 *defective*.

**Reproduction**

```bash
python -c "print(998.21*9.80665/1000, 9.81*1000/9.80665, 1000*9.81/1000)"
python -c "print(998.21/(0.45359237/0.3048**3), 1000/(0.45359237/0.3048**3))"
python -c "from pypdf import PdfReader; p='pilot/references/public/full_books_civil_engineering/Das, Sobhan — Principles of Geotechnical Engineering.pdf'; t=' '.join((PdfReader(p).pages[91].extract_text() or '').split()); i=t.find('9.81 kN/m'); print(t[i-260:i+120])"
```

**Assessed impact.** *On classification: HIGH.* Two of civil's three `[KNOWN-DEFECTIVE]`
rows — two of the corpus's nineteen — are not defects. The register's §1 "Decision:
correct all three to the declared conditions, or restate the conditions" is half right: for
these two only the second branch is correct, and the register does not say so.

*On items if left alone: ZERO.* All seven consumers print γw verbatim in the question text
(`"Taking the unit weight of water as {gamma_w:.2f} kN/m^3, determine …"`,
`phase_and_index_properties.py:91`, `:183`, `:313`; `permeability_seepage_effective_stress.py:169`,
`:266`; `stress_distribution_and_consolidation.py:128`), so it is a stated given under the
repo's own given-values rule and the graded answer follows from it.

*On items if "corrected": NEGATIVE.* Replacing 9.81 with 9.79 would (i) put a number no
geotechnical text uses into the question stem, (ii) desynchronise this file from
`DAS_NATURAL_STATE_SOILS`, whose dry-unit-weight column Das derives on p.92 using g = 9.81,
so that `gamma_sat = gamma_d + (e/(1+e))*gamma_w` — the relation this file's own header
prescribes — would mix two conventions, and (iii) move the item pool (a P6 event) for no
gain in truth.

**Recommended class:** `[ON-DISK:LOCAL-ONLY] … Das, Sobhan … @ page=92 text="9.81 kN/m3"`
with `@domain: none (conventional value; rho taken as 1000 kg/m^3 exactly)`, and a note
that the *measured* γw at 20 °C is 9.789 kN/m³. Equivalently `[BY-DEFINITION]`. Not
`[KNOWN-DEFECTIVE]`, and **not a value change**.

### F-2 — `CONFIRMED`. The cited artefact cannot speak to half the declared domain

`docs/references/nist_fluid_properties/water_C7732185_isobar_1atm.tsv` holds **19 rows,
T = 293.15 … 453.15 K in 10 K steps**. There is no row at 288.15 K, and none below 293.15 K
at all.

Both unit-weight rows declare `@domain: T=288.15..293.15 K`, and the comment reasons "at
the ~15-20 C this line claims, gamma_w = rho*g_n = 9.7891 kN/m^3" — a number obtained at
293.15 K only, because the other endpoint is not in the file. The spec's own NIST-TSV
locator rule is explicit that "the row must be **on the grid**, not near it"; here half the
declared window is off the artefact entirely.

The conclusion's *direction* survives (at 288.15 K, ρ ≈ 999.1 kg/m³ gives γw ≈ 9.798 →
9.80 at 3 s.f., still not 9.81) — but that check cannot be run from anything on disk, and
the tag quotes its error to four figures as though it could.

**Reproduction**

```bash
python -c "L=open('docs/references/nist_fluid_properties/water_C7732185_isobar_1atm.tsv',encoding='utf-8').read().splitlines(); print(len(L)-1, L[1].split(chr(9))[0], L[-1].split(chr(9))[0])"
```

**Assessed impact.** MEDIUM. A `[KNOWN-DEFECTIVE]` verdict whose error is stated to 4
significant figures rests on a one-endpoint evaluation of a two-endpoint domain. Compounds
F-1: the domain that drove the verdict is both unsourced *and* unverifiable from the cited
file.

### F-3 — `PLAUSIBLE`. The viscosity row is correctly classed, but does not belong in the same bundle as the other two

`WATER_KINEMATIC_VISCOSITY_M2_S` is the only one of the three whose `@domain: T=293.15 K`
sits **on** the grid. μ/ρ = 0.0010016/998.21 = 1.0033961e-6; the literal 1.004e-6 is not
its 4-s.f. rounding (1.003e-6 is). By SPEC-CHANGE 23 the class is honest, and I do not
dispute it.

What I dispute is the framing. The register's heading is "Civil — 3, all water at 20 C" and
its text says "**One defect with one cause**". That is not so: the two unit weights follow
from ρ = 1000 exactly (F-1), while the viscosity follows from a 4-digit μ of 1.002e-3 Pa·s
— the implementer's own comment says so. Different cause, different class, different
decision. The +0.060% is also below every display precision in the corpus, and
`git grep` finds **no consumer of this constant anywhere under `data/`** — it is defined
and never read.

**Reproduction**

```bash
git grep -n "WATER_KINEMATIC_VISCOSITY_M2_S\|WATER_DENSITY_KG_M3\|UNIT_WEIGHT_WATER_PCF" 2fe50c7 -- data/   # constants.py only
```

**Assessed impact.** LOW on items (nil — unconsumed), MEDIUM on the register, which asks the
owner for one decision where three rows need two different ones.

### F-4 — `CONFIRMED`. `SCS_IA_RATIO`'s prose page reference is off by one

The comment reads `TR-55 p.2-2: "the following empirical equation: Ia = 0.2S. [Eq. 2-2]"`.
The machine locator (`page=10 text="Ia = 0.2S"`) is **correct and resolves** — I read the
page. But that page's own footer is **2-1**, not 2-2; Eq. 2-2 sits on document page 2-1.
(The adjacent `SCS_CURVE_NUMBERS` doc-page labels 2-5 / 2-6 / 2-7 are all correct.)

**Reproduction**

```bash
python -c "from pypdf import PdfReader; print(' '.join((PdfReader('docs/references/civil/nrcs_tr55_urban_hydrology.pdf').pages[9].extract_text() or '').split())[-60:])"
```

**Assessed impact.** Cosmetic. The resolvable locator is right; only the human-readable
page label is wrong.

### F-5 — `CONFIRMED` (not a defect; recorded so it is not later "fixed"). The AISC SI block is not a conversion of the US block

I verified all 196 leaves independently and they are exact. I then tried the test that was
applied to `LIVE_LOADS_KPA` — is the SI column a conversion of the US column? — and it
fails, by up to **−1.44%** (`W14X30` → `W360X44`: 30 lb/ft converts to 44.64 kg/m, the
database says 44). Eleven of the fourteen shapes disagree by more than 0.2% in at least one
field.

That is correct and expected: `W360X44` is AISC's **nominal designation mass**, not a
conversion of 30 lb/ft, and the SI section properties are the database's own rounded
column. Both blocks are read from the workbook in their own units, and the file transcribes
both faithfully.

**Assessed impact.** None today. Recorded because a future sweep that generalises the
soft/hard-conversion test from `LIVE_LOADS_KPA` to every dual-unit table in the corpus
would report up to 14 spurious civil defects here.

---

## 4. Falsification attempts that failed

- **`AISC_W_SHAPES` (the largest claim in the branch).** Read the workbook my own way —
  `openpyxl`, header row parsed by me, US fields at the first occurrence of each header and
  SI at the second (label columns 3 and 86, 1-based) — and compared all 196 leaves:
  **0 mismatches.** Then tried three ways to break the check itself: (i) looked for
  duplicate `AISC_Manual_Label` values that would let a wrong row satisfy a lookup — none
  in 2301 rows, for either block; (ii) confirmed the US and SI label lookups land on the
  *same worksheet row object* for each shape, so the `si_label` pairing is genuinely tested
  rather than assumed; (iii) checked the declared SI multipliers (Ix in 10⁶ mm⁴, Sx/Zx in
  10³ mm³) against the US columns — consistent to ~0.1% where they should be. The claim
  holds as written.
- **`LIVE_LOADS_KPA`.** Tried to show the "soft conversion 0.048" story was a coincidence.
  It is not: all four rows sit at **+0.25013%** over the SP 811 hard conversion — the
  identical ratio 0.048/0.04788026, to five figures, across four independent rows. That is a
  signature. The corridor row is +0.0412%. I verified the SP 811 factor two ways: read
  `4.788 026 E+01` off p.65 (`pound-force per square foot → pascal`) and recomputed it
  exactly as 4.4482216152605/0.09290304 = 47.88025898. Claim confirmed.
  *One refinement:* the note calls the corridor row "neither that nor the SP 811 factor",
  which is true as written, but 4.79 = 100 × 0.0479 is the same factor rounded to 3 s.f.
  where the other four use 2 s.f. The column is one factor at two precisions, not two
  unrelated rules — a cleaner statement of the same fact, and it strengthens the
  `[UNVERIFIED]` class rather than weakening it.
- **`PERMEABILITY_RANGES_CM_S`.** Read Das PGE p.242 myself. Table 7.1 is on that page under
  the claimed caption; the four non-clay rows match endpoint for endpoint; the clay row is
  given as an upper bound with **no floor**, exactly as the comment says. Could not break it.
  Also checked consumption: the only row any template reads is `["coarse sand"]`
  (`permeability_seepage_effective_stress.py:67`), so the unsourced `1e-8` floor reaches
  **zero items** — the honest tag has no live blast radius.
- **The public-document page locators.** All five pages carry what they claim, confirmed by
  the pages' own footers: HDS-4 p.199 = Table B.2 (footer B-2), p.201 = Table B.3 (B-4);
  HEC-22 p.56 = Table 3-1 (3-6); TR-55 pp.14/15/16 = Tables 2-2a/b/c (2-5/2-6/2-7), each
  carrying the `Ia = 0.2S` / average-runoff footnote that the `@domain` line claims. Every
  transcribed value reproduces: 12 channel rows, 5 conduit rows, 16 rational-C rows, 24
  curve-number quadruples including the impervious-area percentages.
- **Tried to catch the wrong source for `RATIONAL_C`.** HDS-4 Table B.1 (p.198) is *also* on
  disk and is nearly the same table — so I looked for a row that separates them. Two do:
  `streets: concrete` (HEC-22 0.80–0.95, HDS-4 0.70–0.95) and `roofs` (HEC-22 0.75–0.95,
  HDS-4 0.70–0.95). The committed values follow **HEC-22**, which is the cited source. The
  citation is the right one of two candidates, not the convenient one.
- **Tried to reopen the Terzaghi factor tables** — the file's own comment admits an earlier
  version mixed factor families at φ ≤ 25 and φ = 40. All 27 entries of
  `TERZAGHI_BEARING_FACTORS` and all 12 of `TERZAGHI_MODIFIED_FACTORS` reproduce from Das
  PGE p.742 and Das PFE p.160. No mixing survives. Both text-layer traps the comments flag
  are real and I hit them independently (PGE Table 16.1 prints the φ = 15 row's label as
  "16"; PFE spaces its caption as "T able 3.2Terzaghi's…") — the anchors chosen work around
  both correctly.
- **Tried to find an unflagged unsourced literal** in `DAS_NATURAL_STATE_SOILS` (9 rows × 4
  fields), `FRICTION_ANGLE_RANGES_DEG` (6 rows) and the Skempton pair. All 44 reproduce.
- **Ran the repo's own resolver** (`test_citations_resolve`) to see whether it disagreed with
  my hand checks: it does not — all pass, civil 34/10/8/16.

---

## 5. Further probing and improvements

*Separately time-boxed to 15 minutes. Forward-looking; none of this is a gate condition.*

1. **Retag, do not correct, the two unit weights.** `@domain: none (conventional value; rho
   taken as 1000 kg/m^3 exactly)` plus `[ON-DISK:LOCAL-ONLY] … Das, Sobhan … @ page=92
   text="9.81 kN/m3"`, keeping a one-line note that the measured γw at 20 °C is 9.789. This
   is a comment-only change — the C3 discipline — and it *removes* two rows from the
   nineteen-defect count honestly rather than by widening a tolerance.
2. **Add a precondition to SPEC-CHANGE 23.** `[KNOWN-DEFECTIVE]` should require that the
   `@domain` it tests against is itself **sourced**. Where the domain came from the row's
   own comment, the defect must first be tested against the hypothesis that *the domain* is
   wrong. Two of civil's three KD rows fail that test; it is worth re-running against the
   other sixteen in the corpus, and mechanical's fourteen are the obvious place to look.
3. **Guard the artefact grid against the declared domain.** Refuse (or flag) a `@domain`
   whose endpoints are not both on the cited artefact's grid. The water isobar begins at
   293.15 K; anything claiming 288.15 K from it is unverifiable by construction. Cheap to
   implement inside the existing TSV locator.
4. **Make `UNIT_WEIGHT_WATER_KN_M3` and `WATER_DENSITY_KG_M3` one object, not two.**
   γw = ρg is a `[DERIVED]` relation the resolver can execute. Today the pair is internally
   inconsistent by 0.21% and nothing notices, because neither is ever checked against the
   other. Whichever convention the owner picks, the two rows should be forced to agree.
5. **Split the register's civil §1.** "One defect with one cause" is wrong: two rows are a
   convention question and one is a 4-digit-μ question. As written, the owner is asked for
   one decision where two are needed, and the wrong decision on the pair is the expensive one.
6. **Two cheap wins.** Fix `p.2-2` → `p.2-1` in the `SCS_IA_RATIO` comment (F-4). And decide
   the clay permeability row: with no consumer, deleting it or demoting it to
   `[POLICY: sampling-only]` both beat shipping an invented floor behind an honest tag.
7. **One ambiguity worth a later pass, not a finding.** `MANNINGS_N_CHANNELS` flattens three
   distinct blocks of HDS-4 Table B.2 (rigid-boundary, minor streams, floodplains) into one
   dict. `"mountain stream, rocky bed": (0.040, 0.050)` is the rigid-boundary block's row;
   the same page's *Mountain Streams* block gives 0.030–0.050 for gravel/cobble beds. Both
   are in the cited table and the committed value is correct — but the key does not say
   which block it came from, and a future re-derivation will have to rediscover that.
8. **Not attempted, and worth someone's time:** the eight civil `[UNVERIFIED]` rows all
   reduce to three absent standards (AISC Manual, ACI 318-19, ASCE 7-22). Of these, ASCE
   7-22 is the one with a measurable consequence — `LIVE_LOADS_KPA`'s five rows are the only
   `[UNVERIFIED]` civil values that are *numbers in a table* rather than single coefficients,
   and the 0.048 signature I confirmed means the SI column can be reconstructed to the digit
   from the psf column once one row of the real table is seen.

---

*Reviewer H, civil branch. No tracked file modified other than this report.*
