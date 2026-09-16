# Phase C2 — interim findings: `CP_PARAMS` is worse than recorded

**Date:** 2026-09-06 · **Branch:** `redesign/phaseC2-thermochemistry-constants`
**Status:** IN PROGRESS — gases complete, condensed phases outstanding

The spec records one defect in `CP_PARAMS`: *"B and C coefficients 10× too large
(`E-2`/`E-5` where Smith–Van Ness gives `E-3`/`E-6`)"*. That is real, and it is
the largest defect by count — but it is **one of four**, and correcting only it
leaves six or more rows still wrong.

---

## Method

Every row is checked against the **NIST Chemistry WebBook (SRD 69)**, retrieved
to `docs/references/nist_webbook/` so each citation resolves to a file on disk.
NIST publishes Shomate coefficients — a *different* functional form fitted to
the same underlying thermochemistry — so agreement between the two is real
evidence rather than a tautology.

Two independent checks are used:

1. **Cp(T) against NIST** at 298/500/800/1200 K. Two good fits of the same data
   agree to about 3%; anything above 5% is a defective row.
2. **The Smith–Van Ness `Cp₂₉₈/R` column**, which the table publishes alongside
   A, B, C, D and which therefore self-checks each row. Where that column is
   available it is decisive: methane computes 4.2171 against a published 4.217,
   ethane 6.3686 against 6.369, propane 9.0109 against 9.011.

---

## The four defect classes

### Class 1 — the recorded exponent error (systematic, most rows)

`B` is written `E-2` where the source column is `10³B`, and `C` is written
`E-5` where the column is `10⁶C`. Correcting it takes N₂ from Cp(300 K) = 42.4
to **29.12**, against NIST's 29.12.

### Class 2 — individually corrupted rows

| Row | Repo (exponents corrected) | Correct | Evidence |
|---|---|---|---|
| `O2(g)` | A 3.630, B 1.794e-3, C −0.658e-6, D +0.061e5 → Cp 34.71 | A 3.639, B 0.506e-3, C 0, D −0.227e5 → Cp **29.39** | NIST 29.38; SVN Cp₂₉₈/R published 3.535, repo row computes 4.175 |
| `C2H2(g)` | Cp₂₉₈ = 68.39 | — | NIST 44.04 → **+55% wrong** |
| `C2H5OH(g)` | Cp₂₉₈ = 74.40 | — | NIST 65.21 → **+14% wrong** |

### Class 3 — column shift

| Row | Repo | Correct | Evidence |
|---|---|---|---|
| `NO2(g)` | C −0.792e-6, D −0.377e5 → Cp 40.27 | C **0**, D **−0.792e5** → Cp **36.98** | NIST 36.97. The `D` mantissa −0.792 was written into the `C` column and an unrelated value put in `D`. |

### Class 4 — gas-phase coefficients under a liquid label

`C6H6(l)` computes Cp₂₉₈ = 85.3. NIST gives **135.69** for liquid benzene — but
~85 is the **gas-phase** value. The row carries gas coefficients under a liquid
key.

Three more liquid rows are low by the same ~35%, consistent with the same
cause: `C7H8(l)` −31.6%, `C3H6O(l)` −35.2%, `C6H14(l)` −35.7%. `H2SO4(l)` is
−59% and `CaCO3(s)` −10.5%. All six need individual sourcing.

---

## Status by row (33 total)

**Verified correct after the Class-1 fix — 21**
`CH4(g)` `CO(g)` `CO2(g)` `H2(g)` `H2O(g)` `N2(g)` `SO2(g)` `NH3(g)` `NO(g)`
`Cl2(g)` `HCl(g)` `H2S(g)` `C2H4(g)` `C3H8(g)` `C4H10(g)` `C2H6(g)` `Air(g)`
`He(g)` `Ar(g)` `Ne(g)` `H2O(l)`

- The 16 gases agree with NIST to within 3.3% across 298–1200 K.
- `He/Ar/Ne` are exact by theory: monatomic ideal gas, Cp = 5R/2 = 20.786, and
  the rows are A = 2.5, B = C = D = 0. No source needed.
- `Air(g)` is a mixture with no NIST entry; it verifies internally against the
  mole-weighted average of N₂/O₂/Ar to within 0.36%.

**Confirmed defective, replacement derived and verified — 2**
`O2(g)` `NO2(g)` — both now within 0.03% of NIST.

**Confirmed defective, replacement not yet derived — 2**
`C2H2(g)` `C2H5OH(g)`

**Outstanding, need individual sourcing — 8**
`C6H6(l)` `C7H8(l)` `C3H6O(l)` `CH3OH(l)` `C6H14(l)` `H2SO4(l)` `NaCl(s)`
`CaCO3(s)` — plus `CH3OH(l)` and `NaCl(s)`, which screen as plausible (−1.0%,
−5.5%) but have not been checked against a cited source.

---

## Source deviation, recorded not hidden

The spec names **Smith–Van Ness–Abbott** as the primary source with NIST as the
cross-check. No fetchable, citable copy of SVN Table C.1 could be obtained: the
accessible copies are Scribd and SlideShare, which cannot be downloaded and are
not appropriate to cite as "edition + page". Fabricating a page citation would
be precisely the failure this project exists to prevent.

**NIST is therefore used as the primary on-disk source** and the SVN functional
form is retained. Where a row must be re-derived rather than corrected, the
coefficients will be fitted to NIST data over the template's temperature range
with the fit residuals reported, so the value is *derived from a citable
artefact by a documented procedure* rather than transcribed from a table nobody
can check.

---

## What this means for the effort estimate

The spec sizes C2 at 15–25 h on the assumption of one mechanical defect. Four
defect classes, with at least 10 rows needing individual sourcing and some
needing a refit, is a materially larger job. `HEATS_OF_FORMATION` by contrast
is looking sound so far — benzene 49.0 against NIST 49.0, CO₂ −393.5 against
−393.51, H₂O(g) −241.8 against −241.826 — with one candidate discrepancy,
propane at −103.8 against NIST's −104.7.
