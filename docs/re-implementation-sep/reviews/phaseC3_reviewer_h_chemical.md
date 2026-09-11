# Phase C3 — Reviewer H (Domain: thermochemistry / materials / fluids), CHEMICAL branch

Frozen ref `f747fb6845b330958e0abfe7c2d7eac9849942bc`. Working tree at that commit, clean.
Scope: whether the NUMBERS in `data/templates/branches/chemical_engineering/constants.py`
are right and whether each row's provenance CLASS is honest. Tag syntax, units, `@domain`
and register completeness are Reviewer G's; other branches are the other H reviewers'.

Every number below was re-derived from the artefact on disk without reading the
implementer's transcription first. All 138 chemical tags were read; 26 values were
independently re-derived.

---

## 1. Verdict

**PASS WITH FINDINGS** — every value I could re-derive is right, and `REAL_FLUID_DATA`'s
re-derivation is sound, but three of the branch's four headline provenance CLAIMS are
themselves wrong in a way that misdirects the decision they ask the owner to take.

---

## 2. Independent re-derivation

Twenty-six values. "Mine" is computed from the artefact named, with no reference to the
committed number.

| # | value | artefact consulted | my number | committed | agree? |
|---|---|---|---|---|---|
| 1 | `REAL_FLUID_DATA` Water `v_g` @ 100 °C | `water_C7732185_saturation_373.15K.tsv`, T=373.15 exact | 1.6718 | 1.6718 | yes |
| 2 | Water `v_f` @ 100 °C | same | 0.0010435 | 0.001044 | yes (4sf) |
| 3 | Carbon Dioxide `v_g` @ 20 °C | `carbon_dioxide_C124389_saturation_293.15K.tsv` | 0.0051493 | 0.00515 | yes (3sf) |
| 4 | Methane `v_g` @ −161 °C | `methane_C74828_saturation_112.15K.tsv` | 0.53089 | 0.53089 | yes |
| 5 | Methanol `v_f` @ 65 °C | `methanol_C67561_saturation_338.15K.tsv` | 0.0013372 | 0.001337 | yes (4sf) |
| 6 | R-12 `v_g` @ 25 °C | `r12_C75718_saturation_298.15K.tsv` | 0.027153 | 0.0272 | yes (3sf) |
| 7 | n-Octane `v_f` @ 126 °C | `octane_C111659_saturation_399.15K.tsv` | 0.0016370 | 0.001637 | yes |
| 8 | *(all 22 sourced `REAL_FLUID_DATA` rows, 44 fields)* | 22 targeted saturation TSVs | max deviation **0.173 %** (R-12 `v_g`, a 3sf rounding) | — | yes, all |
| 9 | `SUBSTANCES_FOR_VAPORIZATION` Methanol ΔH_vap | `methanol_..._saturation.tsv`, H(v)−H(l) at 0.101325 MPa | 35.283 kJ/mol | 35.3 | yes (rounding) |
| 10 | Propane ΔH_vap | `propane_..._saturation.tsv` | 18.767 | 19.04 | +1.45 % |
| 11 | Toluene ΔH_vap | `toluene_..._saturation.tsv` | 33.235 | 33.48 | +0.74 % |
| 12 | Water ΔH_vap | `water_..._saturation.tsv` | 40.666 | 40.66 | −0.01 % |
| 13 | *(all 11 derivable rows)* | 11 saturation TSVs | **all within 1.45 %**, 8 within 0.2 % | — | see F-2 |
| 14 | `COMMON_LIQUIDS` Methanol μ @ 293.15 K | `methanol_C67561_isobar_1atm.tsv` | 5.853e-4 Pa·s | 5.44e-4 | **no, −7.06 %** |
| 15 | Methanol μ @ 298.15 K | same, interpolated | 5.4576e-4 | 5.44e-4 | **yes, −0.32 %** |
| 16 | n-Hexane ρ @ 293.15 K | `hexane_C110543_isobar_1atm.tsv` | 659.38 kg/m³ | 654.8 | **no, −0.69 %** |
| 17 | n-Hexane ρ @ 298.15 K | same, interpolated | 654.84 | 654.8 | **yes, −0.01 %** |
| 18 | Water μ @ 293.15 K (control) | `water_C7732185_isobar_1atm.tsv` | 1.0016e-3 | 1.002e-3 | yes, +0.04 % |
| 19 | `COMMON_GASES` N₂ ρ, μ @ 293.15 K | `nitrogen_..._isobar_1atm.tsv` | 1.1648, 1.7573e-5 | 1.165, 1.76e-5 | yes |
| 20 | `COMMON_GASES` SF₆ μ @ 293.15 K | `sulfur_hexafluoride_..._isobar_1atm.tsv` | 1.4992e-5 | 1.59e-5 | +6.06 % |
| 21 | `COMMON_GASES` Steam ρ | NIST water isobar @ 373.15 K = 0.59761; ideal gas @ 293.15 K = 0.7489 | — | 0.747 | **no, see F-4** |
| 22 | `GAS_MOLECULAR_PARAMS` Cl₂ σ, ε/k | Svehla TR R-132 Table I(a), PDF p.24 | 4.217 Å, 316.0 K | 4.40, 316.0 | ε/k yes, σ no |
| 23 | Xe σ, ε/k | Svehla Table I(a), PDF p.26 | 4.047 Å, 231.0 K | 4.10, 231.0 | ε/k yes, σ no |
| 24 | SF₆ σ, ε/k | Svehla Table I(a), PDF p.26 (`SF_`) | 5.128 Å, 222.1 K | 5.51, 222.1 | ε/k yes, σ no |
| 25 | `CRITICAL_PROPERTIES` ethanol / acetone / p-xylene `Tc` | the three WebBook `_phase_change.html` pages | 513.9 / 508.2 / 616.2 all present in the page | same | yes |
| 26 | `GRAVITATIONAL_ACCELERATION` | `codata_2022/allascii.txt` | 9.80665 (exact) → 9.81 | 9.81 | yes |

Also spot-checked against `nist_webbook/shomate_coefficients.json`: `NO2(g)` −33.1 → 33.1 ✓,
`C3H8(g)` −103.8 = the Prosen and Rossini 1945 entry the comment names ✓, `CH3OH(l)` −238.6
does sit between the Baroody (−238.4) and Green (−238.9 ± 3.6) entries the comment names ✓.

---

## 3. Findings

### F-1 · CONFIRMED · `COMMON_LIQUIDS` viscosity is a **25 °C column**, not "a source difference"

The register (§3) and the table header both state the finding as *"for every organic liquid
that has a NIST isobar on disk, the VISCOSITY here is low by 4.6-7.2 % at the very conditions
this table declares"*, conclude *"a one-sided gap across four independent species is a source
difference, not scatter"*, and prescribe *"re-source the viscosity column"*.

The measured deviations reproduce exactly (−7.06 / −7.16 / −4.62 / −6.12 %). The **diagnosis
does not**. Re-evaluate the same four rows at 298.15 K and the gap collapses:

| row | μ dev @ 293.15 K (declared) | μ dev @ 298.15 K | ρ dev @ 293.15 K | ρ dev @ 298.15 K |
|---|---|---|---|---|
| Methanol | −7.06 % | **−0.32 %** | +0.04 % | +0.63 % |
| Benzene | −7.16 % | **−0.73 %** | −0.28 % | +0.33 % |
| Toluene | −4.62 % | **+1.12 %** | +0.00 % | +0.54 % |
| n-Hexane | −6.12 % | **−1.50 %** | −0.69 % | **−0.01 %** |
| Water *(control)* | **+0.04 %** | +11.41 % | −0.00 % | +0.13 % |

Water is the control and it settles it: water's viscosity is a 20 °C value (+0.04 %) and is
11.4 % wrong at 25 °C, while all four organics are the reverse. These are the standard
handbook 25 °C entries — toluene 0.560 mPa·s and n-hexane 0.294 mPa·s are the CRC 25 °C
values to three figures. n-Hexane's *density* is at 25 °C too (−0.01 %), while toluene's and
methanol's are at 20 °C.

So the table is not one source read at the wrong temperature; it is **mixed row by row and
field by field between 20 °C and 25 °C**, under a single `@domain: T=293.15 K` declaration.

*Why the class is still wrong in substance.* `[UNVERIFIED]` is defensible for a table whose
majority (honey, blood, engine oils) has no source. But the tag's stated FINDING asserts a
conclusion ("a source difference") that the artefact on disk falsifies, and the remedy it
prescribes — re-source the viscosity column — would not fix n-hexane's density and would
throw away four values that are *correct at 25 °C*. The cheaper and more honest fix is the
one the civil branch was told to consider for water: restate the conditions, or move the five
rows to one temperature.

*Impact: real and graded.* `COMMON_LIQUIDS` is consumed by
`transport_phenomena/shell_momentum_balances.py` (lines 84, 179, 287) and
`transport_phenomena/viscosity_and_momentum_transport.py` (lines 29, 105, 280) — falling-film,
pipe-flow and Reynolds-number items. A 7 % viscosity error propagates linearly into Re and
into every laminar/turbulent regime call made at the margin.

```bash
python - <<'EOF'
D='docs/references/nist_fluid_properties/'
def interp(fn,T,c):
    L=open(D+fn,encoding='utf-8').read().splitlines(); h=L[0].split('\t'); j=h.index(c)
    p=[(float(r.split('\t')[0]),float(r.split('\t')[j])) for r in L[1:] if r.split('\t')[j]]
    for k in range(len(p)-1):
        if p[k][0]<=T<=p[k+1][0]:
            f=(T-p[k][0])/(p[k+1][0]-p[k][0]); return p[k][1]+f*(p[k+1][1]-p[k][1])
for n,fn,rho,mu in [('Methanol','methanol_C67561',791.3,0.544e-3),('Benzene','benzene_C71432',876.5,0.601e-3),
                    ('Toluene','toluene_C108883',866.9,0.560e-3),('n-Hexane','hexane_C110543',654.8,0.294e-3),
                    ('Water','water_C7732185',998.2,1.002e-3)]:
    f=fn+'_isobar_1atm.tsv'
    for T in (293.15,298.15):
        print(n,T,'rho%+.2f%%'%(100*(rho-interp(f,T,'Density (kg/m3)'))/interp(f,T,'Density (kg/m3)')),
                  'mu%+.2f%%'%(100*(mu-interp(f,T,'Viscosity (Pa*s)'))/interp(f,T,'Viscosity (Pa*s)')))
EOF
```

---

### F-2 · CONFIRMED · `SUBSTANCES_FOR_VAPORIZATION` — the "7 of 12" count is wrong; **11 of 11** derivable rows agree

The header says the NIST-saturation derivation *"reproduces 7 of them to within 0.04-1.45 %"*
and concludes *"Claiming [DERIVED] would be false for the rows that miss"*. I derived
ΔH_vap = H(v) − H(l) at exactly 0.101325 MPa from each species' full saturation TSV:

| row | my ΔH_vap (kJ/mol) | committed | dev | rounding of mine? |
|---|---|---|---|---|
| Water | 40.666 | 40.66 | −0.01 % | no |
| Methanol | 35.283 | 35.3 | +0.05 % | **yes** |
| Propane | 18.767 | 19.04 | +1.45 % | no |
| n-Butane | 22.419 | 22.44 | +0.09 % | no |
| n-Hexane | 28.881 | 28.85 | −0.11 % | no |
| Benzene | 30.753 | 30.8 | +0.15 % | **yes** |
| Toluene | 33.235 | 33.48 | +0.74 % | no |
| Ammonia | 23.327 | 23.3 | −0.11 % | **yes** |
| Nitrogen | 5.580 | 5.57 | −0.17 % | no |
| Oxygen | 6.818 | 6.82 | +0.04 % | **yes** |
| Argon | 6.437 | 6.43 | −0.11 % | no |

**Every one of the eleven agrees**, eight of them within 0.2 %. The twelfth unsourced row,
Mercury, has no NIST fluid file at all, so it is not "a row that misses" — it is a row that
cannot be tried. **There are no rows that miss.**

Two further defects in how this is written down:

1. **The two documents contradict each other.** `constants.py` names 4 roundings *plus* 7
   non-roundings — eleven species — while saying "7 of them". The register (§3) says
   *"reproduces 7 of the 12 unsourced rows … but only methanol, benzene, ammonia and oxygen
   are roundings of it"*, which reads the 4 as a **subset** of the 7. One of the two is
   wrong regardless of my derivation; the arithmetic 7 + 4 = 11 says both understate the
   denominator.
2. **The stated reason for the class does not exist.** The table stays `[UNVERIFIED]`
   because "[DERIVED] would be false for the rows that miss". No row misses. The honest
   statement of the residual is *"a 1-atm NIST derivation reproduces all eleven derivable
   rows to within 1.45 %, eight within 0.2 %; only Mercury has no artefact"* — and the
   question put to the owner ("decide whether ~1 % agreement is a derivation") is a
   materially different question when the alternative is 11/11 rather than 7/12.

I do **not** claim the class should flip to `[DERIVED]`: propane at +1.45 % and toluene at
+0.74 % are not roundings of the on-disk derivation, and the table is consumed verbatim by
`thermodynamics/heat_effects.py:128`, which prints ΔH_vap into the question stem. The class
is arguably right; **its stated reason is false**, which is the thing a domain reviewer is
here to catch.

```bash
# see §2 script; bracket each saturation TSV at P=0.101325 MPa and take H(v)-H(l)
```

---

### F-3 · CONFIRMED · `GAS_MOLECULAR_PARAMS` — all 17 rows resolve against Svehla, and the defect is **within** rows, not between tables

The tag claims 15 of 17 rows were parsed, that *"Xenon and Chlorine have no unambiguous token
in the OCR, so they are not claimed either way"*, that ε/k agrees for **11 of 14**, and that
σ agrees for **3 of 15**. Reading Table I(a) myself off PDF pp.22-26:

| row | Svehla σ | committed σ | Svehla ε/k | committed ε/k | PDF p. |
|---|---|---|---|---|---|
| Air | 3.711 | 3.62 | 78.6 | 97.0 | 22 |
| N₂ (`N_`) | 3.798 | 3.70 | 71.4 | 95.05 | 25 |
| O₂ (`02`) | 3.467 | 3.46 | 106.7 | **106.7** | 25 |
| CO₂ | 3.941 | **3.94** | 195.2 | **195.2** | 23 |
| Ar | 3.542 | **3.54** | 93.3 | **93.3** | 22 |
| He | 2.551 | **2.55** | 10.22 | **10.22** | 22 |
| Ne | 2.820 | 2.92 | 32.8 | **32.8** | 25 |
| Kr | 3.655 | 3.69 | 178.9 | **178.9** | 22 |
| **Xe** | 4.047 | 4.10 | 231.0 | **231.0** | 26 |
| CH₄ | 3.758 | 3.78 | 148.6 | **148.6** | 23 |
| C₂H₆ | 4.443 | 4.42 | 215.7 | **215.7** | 23 |
| C₃H₈ (`C3Hs`) | 5.118 | 5.06 | 237.1 | **237.1** | 23 |
| n-C₄H₁₀ | 4.687 | 5.47 | 531.4 | **531.4** | 23 |
| H₂ | 2.827 | 2.93 | 59.7 | 33.3 | 24 |
| NH₃ | 2.900 | 2.92 | 558.3 (`55& 3`) | **558.3** | 25 |
| **Cl₂** (`C_`) | 4.217 | 4.40 | 316.0 | **316.0** | 24 |
| **SF₆** (`SF_`) | 5.128 | 5.51 | 222.1 | **222.1** | 26 |

Three corrections to the claim:

1. **Xenon and Chlorine and SF₆ do resolve.** Their OCR tokens are mangled (`C_` for Cl₂,
   a column-shifted `Xe`, `SF_` for SF₆) but each one's ε/k matches the committed value to
   the last digit — 316.0, 231.0, 222.1 — at exactly the alphabetical position the species
   belongs in. That agreement *is* the unambiguous identification. The correct denominators
   are **17 of 17 rows parsed**, ε/k agreeing for **14 of 17**, not 11 of 14.
2. **Krypton's σ is Svehla too.** Svehla gives a second, thermal-conductivity-derived set at
   the foot of p.26 (Ar 3.408/119.9, He 2.608/10.22, Kr **3.690**/164.7, Ne 2.764/40.2,
   Xe 4.082/206.9). Committed Kr σ = **3.69** is that value exactly, while committed Kr ε/k
   = 178.9 is the viscosity-derived one. So the row pairs one Svehla column with the other
   Svehla column. I flag this as a strong coincidence rather than a certainty — but it means
   "sigma agrees for only 3" is not safe as stated.
3. **The real defect is within-row, and the register does not say it.** Fourteen rows carry
   **Svehla's ε/k with a σ that is not Svehla's**. Lennard-Jones σ and ε/k are a jointly
   fitted *pair*; splitting them produces a parameter set no source ever published and that
   reproduces nothing. The register frames this as "the table mixes sources … re-source
   sigma, or name the compilation the table follows", which reads as a citation problem. It
   is a physics problem, and σ enters viscosity **squared**:

| gas | μ_CE with committed σ | with Svehla σ | NIST @ 293.15 K | committed dev | Svehla dev |
|---|---|---|---|---|---|
| **n-Butane** | 5.354e-6 | 7.293e-6 | 7.250e-6 | **−26.2 %** | +0.6 % |
| **SF₆** | 1.309e-5 | 1.511e-5 | 1.499e-5 | **−12.7 %** | +0.8 % |
| Neon | 2.870e-5 | 3.077e-5 | 3.136e-5 | −8.5 % | −1.9 % |
| Xenon | 2.200e-5 | 2.258e-5 | 2.267e-5 | −2.9 % | −0.4 % |
| Krypton | 2.425e-5 | 2.472e-5 | 2.496e-5 | −2.8 % | −1.0 % |
| Methane | 1.084e-5 | 1.097e-5 | 1.091e-5 | −0.7 % | +0.5 % |

Restoring Svehla's own σ brings every one of these inside ~2 % of the NIST isobar; the
committed σ puts n-butane 26 % low. The consumer,
`viscosity_and_momentum_transport.py:192`, hands σ straight into
μ = 2.6693e-6·√(MT)/(σ²Ω) and **prints σ into the question stem**, so a student is asked to
compute a viscosity from a diameter that is 17 % too large for butane. `[UNVERIFIED]` is the
right class; the register entry understates the severity by a full order of magnitude and
mis-names the remedy.

```bash
python -c "
from pypdf import PdfReader
r=PdfReader('docs/references/nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf')
for i in (21,23,24,25):
    print([l for l in (r.pages[i].extract_text() or '').splitlines()
           if l.strip().startswith(('C_','SF_','Xe','Kr','Ne','N_','02','Air','Ar ','He '))])"
```

---

### F-4 · PLAUSIBLE · `COMMON_GASES` "Steam (Water Vapor)" is +25 % against its own stated condition

The row is `"Steam (Water Vapor)": (0.747, 1.02e-5)` with the inline comment
`# At 100°C (373 K), 1 atm`.

- NIST water isobar at 373.15 K, 1 atm: **0.59761 kg/m³**. Committed 0.747 is **+25.0 %**.
- Ideal-gas H₂O at 293.15 K, 1 atm: **0.7489 kg/m³**. Committed 0.747 is −0.25 % of that.

So the value is water vapour at the table's *declared* 20 °C, computed as an ideal gas —
not steam at the 100 °C its own comment names. Either way the row is wrong: at 20 °C and
1 atm water is a liquid, so a 20 °C vapour density is unphysical, and at 100 °C the number
is a quarter high.

This also puts a hole in the register's characterisation of the table ("the rest sit 0.1-6 %
off it, worst SF6 viscosity +6.06 %"). A +25 % row is outside that band; it was missed
because water at 20 °C is a liquid in the isobar file, so an automated 293.15 K sweep
silently skips it. Untagged, so the row currently makes no claim at all — but the table's
blanket `[UNVERIFIED]` inherits a stated error range that does not cover it.

```bash
python -c "
R=8.314462618
print('ideal H2O @293.15K:',18.01528e-3*101325/(R*293.15))
L=open('docs/references/nist_fluid_properties/water_C7732185_isobar_1atm.tsv',encoding='utf-8').read().splitlines()
h=L[0].split('\t'); j=h.index('Density (kg/m3)')
print([r.split('\t')[j] for r in L[1:] if abs(float(r.split('\t')[0])-373.15)<1e-9])"
```

---

### F-5 · PLAUSIBLE (minor) · `CH4(g)` heat of formation cites one measurement and uses another

`# [ON-DISK] NIST 74-82-8, -74.6+/-0.3 (Manion 2002, adopting Gurvich 1991)` sits above
`"CH4(g)": -74.8`. The JSON's alternative, `dHf_alternative: -74.87 (Chase, 1998)`, is the
nearer of the two to the committed value. The row is inside the cited ±0.3 (0.2 away), so it
is not the G-2 defect Phase C2 corrected in C₂H₆/C₃H₈/CH₃OH(l) — those sat *outside*. But it
is the same shape, and this table's header promises each citation "names the MEASUREMENT the
row uses". Raised as a probe, not a defect; C2 owns this table and I found no evidence C3
disturbed it.

---

## 4. Falsification attempts that failed

- **`REAL_FLUID_DATA`'s re-derivation claim, attacked hardest.** All 22 sourced rows, 44
  fields, re-read from the targeted saturation TSVs. Every row is **exactly on the grid** —
  each targeted file carries three rows and the row read is the stated temperature to the
  last digit (`|T − T_stated| = 0`), so there is no "near it" interpolation hiding anywhere.
  Every committed value is the NIST value at the stated `precision=Nsf`; largest deviation
  0.173 % (R-12 `v_g`, 0.027153 → 0.0272 at 3sf). Each row's `temp_C` converts to its file's
  Kelvin exactly (100 °C → 373.15 K and so on). The C3 corrections the comments claim are
  real and large: CO₂ `v_g` was +265 %, methane `v_g` +242 %, ethane `v_g` +43 % before this
  phase. I could not break this table.
- **The three `[UNVERIFIED]` rows inside `REAL_FLUID_DATA`.** Ethanol, Acetone and R-410A —
  I confirmed there is genuinely no ethanol, acetone or R-410A file among the 36 fluids, and
  that the ethanol/acetone WebBook species pages carry no saturated volumes. The class is
  honest and the reason given is exactly right. Iso-pentane's identification as NIST's
  2-methylbutane (CAS 78-78-4 = the file's `C78784`) also checks out.
- **`CRITICAL_PROPERTIES`' three "locator-only" rows.** I looked for a reason to call these
  overclaimed, because the pages' *headline* values are averages that differ from the rows
  (ethanol `Tc 514. ± 7. AVG`, acetone `508. ± 2. AVG`, p-xylene `617. ± 3. AVG`). But
  513.9, 508.2 and 616.2 each do appear in their page's individual data points. The tag's
  own wording — "Locator-only: the value sits in a positional HTML table the resolver does
  not parse, so it pins the page, not the number" — is precisely the right class, and the
  blanket `[UNVERIFIED]` over the other 24 rows is honest.
- **`COMMON_GASES`' ten tagged pairs.** All ten equal the NIST isobar at 293.15 K within
  their stated precision (N₂ +0.02 %/+0.15 %, Ar −0.05 %/−0.03 %, He −0.19 %/−0.09 %,
  CH₄ −0.02 %/−0.13 %, CO₂ μ +0.17 %, O₂ ρ −0.02 %, H₂ ρ +0.06 %). The untagged rows sit
  where the register says (SF₆ μ +6.06 %, butane μ +3.45 %, Kr μ +2.16 %) — except the steam
  row, F-4.
- **`REACTIONS` and `COMBUSTION_REACTIONS`' `[DERIVED]` claims.** Counted atoms on both sides
  of all 4 + 11 reactions from the species keys by hand. All balance, including the ammonia
  one where the N₂ product is written 13.28 = 2 + 11.28 with the arithmetic in a comment.
  Theoretical air is 3.76 N₂/O₂ throughout. Honest.
- **`CP_VALID_T_MAX` and `CP_PARAMS_COMBUSTION`.** Classes are honest: the former is
  `[DERIVED]` from this repo's own fit ranges and says so; the latter's per-row residuals
  are the ones a 298-3000 K refit produces. `CP_PARAMS`/`HEATS_OF_FORMATION` spot checks
  (§2) found nothing to suggest C3 disturbed C2's work.
- **`AIR_COMPOSITION`, `POWER_LAW_FLUIDS`, `SUBSTANCES_FOR_HEATING`.** `[UNVERIFIED]` with
  correct reasons; the three air fractions do sum to 0.99964 as claimed, and the point that
  K for a shear-thinning fluid is not a substance constant is the right domain objection.
- **`GRAVITATIONAL_ACCELERATION`** resolves to CODATA's exact 9.806 65 at 3sf.
- **Resolver.** `python -m tests.constants_integrity.test_citations_resolve` passes:
  chemical 138 tags, 105 resolved, 7 locator-only, 26 stated, 0 unresolvable, 0 LEGACY —
  matching the register's table exactly.

---

## 5. Further probing and improvements

Forward-looking; separately time-boxed to 15 minutes and not part of the verdict.

1. **Settle `COMMON_LIQUIDS` by declaring 25 °C, not by re-sourcing.** Four of the five
   NIST-checkable rows are already correct at 298.15 K. Changing `@domain` to
   `T=298.15 K` fixes methanol, benzene, toluene and n-hexane's density and viscosity at a
   stroke and moves only water (μ 1.002e-3 → 0.8901e-3) and toluene/methanol density by
   <0.7 %. That is a far smaller item-pool disturbance than re-sourcing the column, and it
   is the option the register never put on the table. Worth costing before C4 decides.
2. **Extend the isobar sweep to species that are liquid only outside 293.15 K.** The steam
   row (F-4) was missed because an automated 293.15 K probe finds water in the liquid phase
   and skips. A sweep that reads each row's *own* inline condition comment — `# At 100°C`,
   `# Temperature is 77 K`, `# Temperature is 90 K` — rather than the table's `@domain`
   would also catch the two cryogen rows in `COMMON_LIQUIDS`, which I did not have time to
   check against `nitrogen_..._saturation.tsv` at 77.15 K and `oxygen_...` at 90.15 K. Those
   two targeted saturation files exist on disk and are currently unused by any tag.
3. **A pair-integrity check for Lennard-Jones rows.** F-3's defect class — a row whose two
   fields come from different fits of the same source — is invisible to a per-field resolver,
   because each field individually "agrees with an artefact". A cheap rule: when a table's
   `@kind: property` row carries fields that a source publishes as a *set*, require the tag
   to name one row of one table, not one column. This would also catch the Krypton case.
4. **Re-derive `SUBSTANCES_FOR_VAPORIZATION` at the NIST T_sat rather than 1 atm.** My
   bracket at 0.101325 MPa lands at T_sat = 372.79 K for water, 0.36 K below the true normal
   boiling point, because the saturation grid is 1 K and I interpolated linearly in P. A
   targeted saturation fetch at each species' normal boiling point (the same mechanism
   `NIST_SATURATION_POINTS` already uses for `REAL_FLUID_DATA`) would remove that last
   ambiguity and probably move propane's +1.45 % — the single worst row and the one the
   `[DERIVED]`-vs-`[UNVERIFIED]` decision turns on.
5. **Mercury.** The only `SUBSTANCES_FOR_VAPORIZATION` row with no path to an artefact: no
   NIST fluid file, and I confirmed its WebBook phase-change page has no vaporization section
   at all (and a corrupt `T c 0. K` entry). The register records this; it is the one row where
   "acquire or replace" is genuinely the only move.
6. **`CRITICAL_PROPERTIES` is the largest cheap win left.** 24 of 27 rows have no artefact,
   and the `fetch_references.py` machinery that produced the 28 species pages already on disk
   would cover most of them. `Vc`, `Zc` and `omega` will not come from WebBook phase-change
   pages, so scope it to `Tc`/`Pc` and expect the acentric factors to stay unverified.

---

*Filed by Reviewer H (chemical). Mandatory task ran to its 40-minute box; §5 to its own 15.*
