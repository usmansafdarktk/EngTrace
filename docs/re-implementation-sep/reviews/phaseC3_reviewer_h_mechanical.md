# Phase C3 — Reviewer H (Domain: thermochemistry / materials / fluids), MECHANICAL branch

Frozen ref `2fe50c72153892d0b1fa67e7ed486ebb36abca0f`. Working tree clean at that commit.
Every number below was re-derived from the primary artefact by reading the PDF page text or
the NIST TSV directly, before consulting the implementer's transcription.

Branch under review: `data/templates/branches/mechanical_engineering/constants.py`
(729 lines; reproduced counts: 14 `[KNOWN-DEFECTIVE]`, 173 `[UNVERIFIED]`, 19 `[ON-DISK]`).

---

## 1. Verdict

**PASS WITH FINDINGS** — all 14 `[KNOWN-DEFECTIVE]` measurements are numerically correct
against the artefacts, but one register percentage disagrees with its own tag, and four
`[UNVERIFIED]` reason-strings state grounds that the artefacts on disk falsify.

---

## 2. Independent re-derivation

Conversion constants used, taken from the tags' own cited authority (NIST SP 811):
`1e3 ksi = 6.894757 GPa`, `1 lb/in^3 = 27679.9 kg/m^3`.

### 2a. The `[KNOWN-DEFECTIVE]` rows (all 14 re-derived)

Page captions were confirmed by reading the page's own first line; every cited page does
carry the cited table.

| # | value | artefact / page consulted | my number | committed | agree? |
|---|---|---|---|---|---|
| 1 | AISI 301 `E_ksi` 27500 | p.277 Table 2.7.1.0(b), E:L = 29.0e3 ksi | −5.172% → −5.2% | −5.2% | yes |
| 2 | AISI 301 `E_GPa` 190 | p.277, 29.0×6.894757 = 199.948 GPa | −4.975% → −5.0% | −5.0% (=199.9) | yes |
| 3 | AISI 301 `nu` 0.30 | p.277, µ = 0.27 | +11.111% → +11.1% | +11.1% | yes |
| 4 | 6061 `E_ksi` 10000 | p.566 Table 3.6.2.0(b1), E = 9.9e3 ksi | +1.010% → +1.0% | +1.0% | yes |
| 5 | 6061 `E_GPa` 68.9 | p.566, 9.9×6.894757 = 68.258 GPa | +0.941% → +0.9% | +0.9% (=68.26) | yes |
| 6 | CP Ti `E_ksi` 16800 | p.899 Table 5.2.1.0(b), E = 15.5e3 ksi | +8.387% → +8.4% | +8.4% | yes |
| 7 | CP Ti `E_GPa` 116 | p.899, 15.5×6.894757 = 106.869 GPa | +8.544% → +8.5% | +8.5% (=106.9) | yes |
| 8 | Ti-6Al-4V `E_ksi` 16500 | p.945 Table 5.4.1.0(b), E = 16.0e3 ksi | +3.125% → +3.1% | +3.1% | yes |
| 9 | Ti-6Al-4V `E_GPa` 114 | p.945, 16.0×6.894757 = 110.316 GPa | +3.340% → **+3.3%** | tag +3.3% / **register +3.4%** | **tag yes, register NO** |
| 10 | Ti-6Al-4V `nu` 0.34 | p.945, µ = 0.31 | +9.677% → +9.7% | +9.7% | yes |
| 11 | Ti-6Al-4V `G_GPa` 41.4 | p.945, G 6.2×6.894757 = 42.7475 GPa | −3.152% → −3.2% | −3.2% (=42.75) | yes |
| 12 | AZ31B `nu` 0.29 | p.841 Table 4.2.1.0(b), µ = 0.35 | −17.143% → −17.1% | −17.1% | yes |
| 13 | Benzene density 876 | `benzene_C71432_isobar_1atm.tsv` @ 293.15 K = **878.92** | −0.332% → −0.33% | −0.33% | yes |
| 14 | R-134a density 1206 | `r134a_C811972_saturation_298.15K.tsv` @ 298.15 K = **1206.7** | −0.058% → −0.06% | −0.06% | yes |

Page-caption check, read from the page text itself — all five match the citation exactly:

- p.277 → "Table 2.7.1.0(b). … AISI 301and Relateda,b,c Stainless Steels" (5 columns; the
  annealed column, which the tags use, is column 1)
- p.566 → "Table 3.6.2.0(b1). … 6061 Aluminum Alloy Sheet"
- p.841 → "Table 4.2.1.0(b). … AZ31B Magnesium Alloy Sheet and Plate"
- p.899 → "Table 5.2.1.0(b). … Commercially Pure Titanium"
- p.945 → "Table 5.4.1.0(b). … Ti-6Al-4V Sheet, Strip, and Plate"

### 2b. Rows that resolve (`[ON-DISK]`), re-derived from the artefact

| value | artefact / page | my number | committed | inside stated tol? |
|---|---|---|---|---|
| `MATERIAL_PROPERTIES['Aluminum 6061-T6'].nu` | p.566, µ 0.33 | 0.33 | 0.33 (`precision=exact`) | yes |
| `MATERIAL_PROPERTIES['Magnesium'].E_ksi` | p.841, E 6.5e3 | 6500 | 6500 (`precision=exact`) | yes |
| `MATERIAL_PROPERTIES['Magnesium'].E_GPa` | p.841, 6.5×6.894757 | 44.816 GPa; row +0.411% | 45 (`tol=0.77%`) | yes — ½ unit of 6.5 = 0.769% |
| `SHEAR_MODULUS_VALUES['Stainless Steel (304)']` | p.277, G 11.2e3 | 77.221 GPa; row −0.286% | 77.0 (`tol=0.45%`) | yes — ½ unit of 11.2 = 0.446% |
| `SHEAR_MODULUS_VALUES['Aluminum 6061-T6']` | p.566, G 3.8e3 | 26.200 GPa; row −0.763% | 26.0 (`tol=1.32%`) | yes — ½ unit of 3.8 = 1.316% |
| `SHEAR_MODULUS_VALUES['Magnesium Alloy (AZ31B)']` | p.841, G 2.4e3 | 16.547 GPa; row −0.286% | 16.5 (`tol=2.09%`) | yes — ½ unit of 2.4 = 2.083% |
| `MATERIAL_DENSITIES['Titanium']` | p.899, ρ 0.163 lb/in³ | 4511.8 kg/m³; row −0.262% | 4500 (`tol=0.31%`) | yes — ½ unit = 0.307% |
| `MATERIAL_DENSITIES['Steel (Carbon)']` | p.62 Table 2.2.1.0(b) AISI 1025, ρ 0.284 | 7861.1 kg/m³; row −0.141% | 7850 (`tol=0.18%`) | yes — ½ unit = 0.176% |
| `FLUID_DENSITIES['Fresh Water']` | `water…isobar_1atm` @ 293.15 | 998.21 | 998 (3sf) | yes |
| `FLUID_DENSITIES` methanol | `methanol…isobar_1atm` @ 293.15 | 791.01 | 791 (3sf) | yes |
| `FLUID_DENSITIES['Toluene']` | `toluene…isobar_1atm` @ 293.15 | 866.89 | 867 (3sf) | yes |
| `FLUID_DENSITIES['Liquid Nitrogen…']` | `nitrogen…sat_77.15K` @ 77.15 | 807.01 | 807 (3sf) | yes |
| `FLUID_DENSITIES['Liquid Oxygen…']` | `oxygen…sat_90.15K` @ 90.15 | 1141.4 | 1141 (4sf) | yes |
| `FLUID_DENSITIES['Liquid Hydrogen…']` | `hydrogen…sat_20.15K` @ 20.15 | 71.096 | 71 (2sf) | yes |
| `FLUID_DENSITIES['Liquid Propane']` | `propane…sat_298.15K` @ 298.15 | 492.36; row +0.130% | 493 (`tol=0.13%`) | yes, at the boundary |
| `FLUID_DENSITIES['Helium (at 20°C)']` | `helium…isobar_1atm` @ 293.15 | 0.16631 | 0.166 (3sf) | yes |

Every tolerance in the branch is exactly "half a unit in the handbook's last printed
digit", computed correctly in all eight MIL-HDBK cases. I could not find a tolerance that
had been loosened to admit a row that would otherwise fail.

### 2c. Handbook table inventory (used to test the "no design table" claims)

I scanned all 1733 pages of the handbook for `Design Mechanical and Physical Properties of
…` and recovered **223** design-table captions
(`scratchpad/captions.txt`). This is the evidence base for Findings 3 and 4.

---

## 3. Findings

### F1 — CONFIRMED (low impact, but it is an error inside a defect record)

The register prints the Ti-6Al-4V `E_GPa` error as **+3.4%**; the tag in `constants.py`
prints **+3.3%**; the true value is **+3.3395%**, which rounds to +3.3%. The register is
wrong and the code is right — the one place in the branch where the two disagree.

```
python -c "print(114/(16.0*6.894757)-1)"      # 0.0333951...
grep -n 'E_GPa.*114' data/templates/branches/mechanical_engineering/constants.py
grep -n 'Ti-6Al-4V .E_GPa' docs/re-implementation-sep/phaseC3_residual_register.md
```

Impact: the register is the artefact the owner will read when deciding this cluster.
A defect record that misstates its own measured error by a digit is the class of error
this review exists to catch, even when the magnitude is small.

### F2 — CONFIRMED (highest impact finding)

**`MATERIAL_PROPERTIES['Steel'].nu = 0.30` is tagged `[UNVERIFIED]` while its own tag
quotes the handbook value that contradicts it.** The tag reads:

```
# [UNVERIFIED] names no grade. MIL-HDBK-5J's carbon-steel table (p.62, AISI 1025) prints E
#   29.0, G 11.0 x10^3 ksi, mu 0.32; choosing a grade is the repo owner's call (C3.5).
'Steel': {'E_GPa': 200, 'E_ksi': 29000, 'nu': 0.30},
```

I re-derived p.62 independently: Table 2.2.1.0(b), AISI 1025 Carbon Steel, E 29.0, Ec 29.0,
G 11.0 ×10³ ksi, **µ 0.32**, ρ 0.284 lb/in³. So against that table the row is:

- `E_ksi` 29000 vs 29.0e3 — **exact**
- `E_GPa` 200 vs 199.948 — **+0.026%**, inside half a unit (0.17%)
- `nu` 0.30 vs 0.32 — **−6.25%**

Two of three fields resolve exactly; the third is off by −6.25%, which is **larger than
five of the fourteen errors the branch does report as `[KNOWN-DEFECTIVE]`** (6061 E_ksi
+1.0%, 6061 E_GPa +0.9%, benzene −0.33%, R-134a −0.06%, Ti-6Al-4V E_ksi +3.1%).

The class is not honest, because the file does not apply its stated policy consistently:

- The **same page 62 / AISI 1025 table** is used as a *resolving* `[ON-DISK]` source for
  `MATERIAL_DENSITIES['Steel (Carbon)']`, where the tag explicitly says "the row names no
  grade, so the grade is the implementer's choice — referred (C3.5)". Naming no grade was
  not an obstacle there.
- `MATERIAL_PROPERTIES['Stainless Steel']` names no grade either (comment says "304"), and
  there the implementer *did* pick a table — p.277, whose own footnote b is what extends it
  to AISI 304 — and filed three `[KNOWN-DEFECTIVE]` rows off it.

So grade-choice is treated as permissible when it makes a row resolve (steel density) and
when it makes a row defective (stainless), but as a blocker precisely where it would
expose a defect in the most-used material in the branch. The net effect is that the single
most common engineering material in the corpus carries a −6.25% Poisson's ratio that the
register never counts.

```
python -c "
from pypdf import PdfReader
t=' '.join((PdfReader('docs/references/mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf').pages[61].extract_text() or '').split())
import re; print(re.search(r'E, 103.{0,120}', t).group(0))"
python -c "print(0.30/0.32-1)"    # -0.0625
```

Impact: either `Steel.nu` is a 15th `[KNOWN-DEFECTIVE]` row, or the three
`Stainless Steel` rows are not defective either — the branch cannot have it both ways. This
is a decision the owner is currently not being asked to make.

### F3 — CONFIRMED (medium impact)

**The "no MIL-HDBK-5J design table is this material" reason carries a copy-pasted
parenthetical that is about the wrong alloy family for 4 of the 7 rows using it.** All
seven rows state the searched set as "(its copper-base tables: C86500 Manganese Bronze;
C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium Rod and Bar;
C17200 Copper Beryllium Mechanical Tubing)". That set is apt for `Copper`, `Brass` and
`Bronze`. It is applied unchanged to **`Tungsten`, `Cast Iron`, `Nickel` and `Lead`**.

For `Nickel` this is materially misleading: my caption scan finds the handbook carries a
large nickel-base chapter that the tag's stated search never touches — A-286 (p.1046),
Hastelloy X (p.1062), Inconel 600 (pp.1068–1070), 625 (pp.1075–1076), 706 (p.1086),
718 (pp.1092–1094), X-750 (p.1118), Waspaloy (p.1131), HAYNES 230 (pp.1137–1138),
HAYNES HR-120 (p.1150) — 15 tables in all.

```
grep -c "its copper-base tables" data/templates/branches/mechanical_engineering/constants.py   # 7
grep -A4 "its copper-base tables" data/templates/branches/mechanical_engineering/constants.py | grep -E "^\s+'"
# then see scratchpad/captions.txt, or re-run the caption scan in §2c
```

Impact: the conclusion ("no table *is* pure nickel / tungsten / gray cast iron / lead")
is, as far as I can tell, still true — the handbook is an aerospace alloy document and has
no pure-element tables. But the *evidence offered for it is false*, and a later reader
auditing "was the handbook searched for nickel?" is told the searched set was five copper
alloys. A reason string in a residual register is a claim; this one does not survive
contact with the artefact.

### F4 — CONFIRMED (medium impact; this is the "could have been sourced" case)

**The 8 rows carrying "names no species and no moisture content, and a wood density is not
one number without both. The USDA Wood Handbook is on disk but was not read to row level"
are misdescribed on both halves of the reason.**

The 8 rows are: `Pine Wood` 500, `Oak Wood` 750, `Cork` 240, `Teak Wood` 630,
`Maple Wood` 740, `Ebony Wood` 1200, `Bamboo` 300, `Cork Board` 240.

1. **Three of the eight are not wood.** `Cork` and `Cork Board` are bark (the phellem of
   *Quercus suber*), and `Bamboo` is a grass. Neither appears in the Wood Handbook's
   species tables at all, so "the USDA Wood Handbook … was not read to row level" is not
   why they are unsourced — they are unsourceable *from that artefact*, which is a
   different class with a different decision.
2. **"Names no species" is false for at least Teak.** *Tectona grandis* L. f. is a single
   species, and it is on disk: FPL-GTR-282 p.72 lists "Tectona grandis L. f. … Teak" in
   its imported-woods nomenclature table.
3. **The artefact does give exactly the two things the tag says are missing.**
   FPL-GTR-282 Table 5–3a/5–3b ("Strength properties of some commercially important woods
   grown in the United States", metric on pp.117–123, inch–pound from p.126) is indexed by
   *common species name* × *moisture content* (Green and 12%) with a **Specific gravity**
   column — e.g. Oak red / Northern red 0.56 green and 0.63 at 12%; Maple sugar 0.56 and
   0.63; Maple red 0.49 and 0.54.

```
python -c "
from pypdf import PdfReader; import re
r=PdfReader('docs/references/usda_wood_handbook/FPL-GTR-282_2021.pdf')
t=' '.join((r.pages[120].extract_text() or '').split())
print(re.search(r'Table 5.3a.{0,200}', t).group(0))
print(re.search(r'Northern red.{0,80}', t).group(0))"
```

Impact: this is the finding the brief asked for — `[UNVERIFIED]` rows that *could* have
been sourced from an artefact already present. I am deliberately **not** asserting the
values are wrong: the handbook's specific gravity is ovendry-mass-over-volume at a stated
volume basis, so turning 0.63 into a kg/m³ at 12% MC needs the handbook's own conversion
and is a genuine piece of work. But the honest class for `Pine/Oak/Teak/Maple/Ebony Wood`
is "artefact on disk, not consulted" (which the tag's last clause admits) rather than
"the row cannot be sourced because it names no species", which the artefact falsifies.

### F5 — PLAUSIBLE (low impact; inflates the register's headline number)

The same material gets two different `[UNVERIFIED]` reasons in the same file depending on
which table it sits in. `Brass`, `Copper`, `Nickel` and `Cast Iron` in `MATERIAL_PROPERTIES`
carry the specific "no MIL-HDBK-5J design table …" reason; `Brass` 8600, `Copper` 8940,
`Nickel` 8900 and `Iron (Wrought)` 7750 in `MATERIAL_DENSITIES` carry the generic
"no source on disk for this material".

```
grep -c "no source on disk" data/templates/branches/mechanical_engineering/constants.py   # 95
grep -c "no MIL-HDBK-5J design table is this material" data/templates/branches/mechanical_engineering/constants.py
```

Impact: the register's headline — "**almost half of every tagged residual in the corpus is
one mechanical table's worth of materials with no primary source at all**" (the 95-row
class) — is partly an artefact of the generic string absorbing rows for which a more
specific and more informative finding had already been established a few hundred lines
earlier. The 95 is a true grep count, but it is not 95 *distinct* unsourceable materials.

---

## 4. Falsification attempts that failed

Everything below I tried to break and could not:

1. **All 14 `[KNOWN-DEFECTIVE]` percentages.** Recomputed from the artefact; 13 match the
   tag to the stated rounding, and the 14th matches the tag (only the register disagrees —
   F1). No defect record overstates its error, and none is a rounding artefact dressed up
   as a defect.
2. **Every cited page really carries the cited table.** I checked the caption line of
   pp.62, 277, 566, 841, 899, 945 against the citation text; all six match verbatim,
   including the handbook's own mangled "AISI 301and Relateda,b,c".
3. **Column selection.** p.277 is a 5-column table and the tags use column 1 (annealed).
   I read the block manually and independently confirmed with the helper: col 1 is
   E:L 29.0 / G 11.2 / µ 0.27. Had the tags silently used a harder column (E:L 26.0) the
   stainless error would flip sign; they did not.
4. **The `milhdbk.py` helper, itself a claim under test.** Ran `design_values` on all six
   pages; its output agrees with my by-hand reading of the page text in every case,
   including the interleaved/block layout split (p.841 and p.62 interleaved; 277/566/899/945
   block) and the density field.
5. **The "`2024-T4` text layer does not parse" claim (register §4).** Reproduced exactly:
   p.374 raises `ValueError: 4 labels and values [0, 0, 0, 0] fit neither layout`. The
   claim is honest and the page is genuinely unreadable from its text layer.
6. **The `Titanium.nu` `[UNVERIFIED]`.** p.899 really does print `µ ... ...` — the helper
   returns `mu: [None]`. Declining to call 0.34 defective there is correct, not evasive.
7. **All 8 MIL-HDBK tolerances.** Each `tol=` equals half a unit in the last printed digit
   of the handbook value, and every row is genuinely inside it. No tolerance was fitted to
   the row.
8. **All 9 NIST locators** resolve to the exact tabulated value at the exact stated
   temperature and column, including the saturated-liquid-vs-isobar distinction (the
   cryogens and R-134a correctly use saturation files; benzene/toluene/methanol correctly
   use the 1-atm isobar). The benzene defect specifically is measured against the isobar,
   which is the right file for a table declaring 20 °C and 1 atm.
9. **Count claims.** 14 KD, 173 UNVERIFIED, 19 ON-DISK and the 95-row generic class all
   reproduce by grep on the branch file exactly as the register states.

---

## 5. Further probing and improvements (separately time-boxed, 15 min)

1. **Resolve F2 before anything else.** Either add `Steel.nu` as a 15th `[KNOWN-DEFECTIVE]`
   at −6.25% against p.62, or downgrade the three `Stainless Steel` rows. As it stands the
   branch's most-used material hides an error larger than five it reports.
2. **Make reason strings falsifiable per row, not per class.** F3 and F4 are both the same
   failure mode: a reason written once for a group and pasted onto members it does not
   describe. A cheap mitigation is a check that any tag naming a searched set
   ("its copper-base tables: …") may only appear on rows whose material is in that family.
3. **Split the 95-row generic class.** At minimum separate "no artefact of this kind exists
   on disk" (polymers, ceramics, foods, rocks) from "an artefact exists but was not read"
   (all the wood rows, and `Aluminum 2024-T4` p.374). They need opposite decisions: the
   first is an acquisition, the second is an afternoon's reading.
4. **The woods are the cheapest win in the branch.** FPL-GTR-282 Table 5–3a is on disk,
   machine-readable from the text layer (I extracted it in seconds), and covers Pine, Oak
   and Maple by species and moisture content. Pick a species and a moisture content, derive
   the density with the handbook's own relation, and five `[UNVERIFIED]` rows become either
   `[ON-DISK]` or honestly `[KNOWN-DEFECTIVE]`. Note my quick look suggests some rows will
   land in the second category (Oak red, northern, SG 0.63 at 12% vs the row's 750 kg/m³).
5. **`Cork`, `Cork Board`, `Bamboo`** should be moved out of the wood class entirely — they
   are a substance-class question (D-033), not a sourcing question.
6. **`p.374` deserves the OCR pass** the register suggests: it is the only place where the
   branch is blocked by a tooling limit rather than by a missing artefact, and it would
   settle `Aluminum 2024-T4` G = 28.0 against a handbook value that is almost certainly
   4.0e3 ksi = 27.58 GPa (i.e. +1.5%, a probable 15th/16th defect).

---

*Reviewer H — mechanical branch. Scope observed: numbers and provenance class only; tag
syntax, unit declarations, `@domain` presence and register completeness belong to
Reviewer G, and the other four branches to the other H reviewers.*
