# Phase C2 — Reviewer H: Domain (thermochemistry)

**Filed:** 2026-09-06 · **Branch reviewed:** `redesign/phaseC2-thermochemistry-constants`
**Roster:** H — Domain (thermochemistry / materials / fluids), phases C2 and C3
**Brief:** one mandatory gate task — *independently re-derive a 10-value sample
directly from primary sources, without reference to the implementer's
transcription.* 35-minute box.

> Worked without sight of the implementer's reasoning (R1.1). Code comments,
> commit messages and the interim findings document were supplied as **claims
> under test**.

**Remediation is recorded in the right-hand column of §3 and triaged in
[`../phaseC2_summary.md`](../phaseC2_summary.md) §11. All CONFIRMED findings are
closed.**

---

## 1. Verdict

**PASS WITH FINDINGS** — every value in my independent 10-sample reproduces from
live NIST within the phase's own tolerance, and the disputed 2326 K
flame-temperature reference is confirmed by my own re-solve; but the polynomials
are extrapolated ~800 K past their validity in the one template that matters
most, and 32 of the 65 in-scope values were never touched.

---

## 2. Independent re-derivation

All NIST values fetched live from `webbook.nist.gov` and evaluated with the
Shomate form myself. I did **not** read the implementer's
`shomate_coefficients.json` until after deriving.

| # | Value | Source consulted (live) | My value | Committed | Agree? |
|---|---|---|---|---|---|
| 1 | `O2(g)` Cp(298.15) — *replaced* | `C7782447&Mask=1`, Shomate 100–700 K | 29.383 | 29.387 | **Yes** (+0.02%) |
| 2 | `NO2(g)` Cp(298.15) — *replaced* | `C10102440&Mask=1`, 298–1200 K | 36.975 | 36.977 | **Yes** (+0.01%) |
| 3 | `C2H2(g)` Cp(298.15) — *refit* | `C74862&Mask=1`, Shomate + tabulated 44.04 | 44.084 | 45.452 | **Marginal** (+3.10%) |
| 4 | `C6H14(g)` Cp — *refit* | `C110543&Mask=1`, Scott 1974, 14 points | 142.60 / 217.28 / 331.37 | 146.59 / 214.24 / 334.04 | **Yes**, max +2.79% |
| 5 | `N2(g)` Cp(298.15) — *exponent fix* | `C7727379&Mask=1`, 100–500 K | 29.124 | 29.116 | **Yes** (−0.03%) |
| 6 | `CO2(g)` Cp(298.15) — *exponent fix* | `C124389&Mask=1`, 298–1200 K | 37.130 | 37.141 | **Yes** (+0.03%) |
| 7 | `SO2(g)` Cp(298.15) — *exponent fix* | `C7446095&Mask=1`, 298–1200 K | 39.874 | 39.876 | **Yes** (+0.01%) |
| 8 | ΔHf `CO2(g)` | `C124389&Mask=1` — CODATA −393.51 ± 0.13 | −393.51 | −393.5 | **Yes** |
| 9 | ΔHf `C3H8(g)` | `C74986&Mask=1` — NIST lists **two**: −104.7 (Pittam & Pilcher 1972) and **−103.8** (Prosen & Rossini 1945) | −103.8 | −103.8 | **Yes** — this **retires** the interim doc's "propane discrepancy": a source choice, not an error |
| 10 | `H2SO4(l)` — *tagged* | `C7664939&Mask=2` — **no condensed-phase Cp section exists**, free or paid-linked | not derivable | Cp₂₉₈ = 56.9 | Tag honest as to NIST; see F-1 |

Bonus: ΔHf `C2H2(g)` NIST 226.73 vs committed 226.7 — agree. `CaCO3(s)`
confirms no free-tier Cp — tag honest.

**Transcription audit.** Having derived independently, I diffed
`shomate_coefficients.json` against my live fetches for O2, NO2, N2, CO2, SO2,
C2H2 and the C6H14 Scott-1974 table. **Every coefficient is digit-for-digit
identical.** The transcription is faithful; I could not find a defect in it.

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **F-1** | **CONFIRMED, moderate.** `[UNVERIFIED]` is the wrong tag for `H2SO4(l)`; it is *known-defective*. The row gives Cp₂₉₈ = 56.91 against a literature ~138.9 — **~59% low**, which the implementer's own interim findings state. It then ships behind a tag reading "we couldn't check this." We *did* check it; it failed. A row demonstrated wrong needs a `[KNOWN-DEFECTIVE]` class, not the same tag as a row nobody looked at. `CaCO3(s)` (−9%) is the milder same case. | **FIXED** — both retagged `[KNOWN-DEFECTIVE]` with the measured error |
| **F-2** | **CONFIRMED, moderate.** The polynomials are extrapolated ~800 K beyond validity in the flame template, and the suite does not test that range. Header declares 298–1500 K; the suite checks 298–1200 K; `adiabatic_flame_temperature` integrates to 2300–2850 K. Against live NIST, CO₂ is +8.9% at 2500 K and breaches the suite's own 8% tolerance at a temperature Acetylene (2844 K) and CO (2623 K) both exceed. **The suite is green only because it never looks there.** | **FIXED + HANDED ON** — `CP_VALID_T_MAX` declared (no validity range existed before, in any form); suite now measures the temperature the template actually reaches and reports the overshoot; the refit-vs-restrict trade handed to Phase 2 as **D-032** |
| **F-3** | **CONFIRMED, low–moderate.** 32 of 65 values not re-derived and carrying no row citation. The diff is one hunk, entirely inside `CP_PARAMS`; `HEATS_OF_FORMATION` and `COMBUSTION_REACTIONS` are byte-identical to master with no `[NIST …]` comment. The gate "all 65 values carry a citation" is **not met on the face of the table** — the evidence exists only in the suite and the JSON. | **FIXED** — all 21 heats of formation now carry an `[ON-DISK]` citation; 54 cited values total |
| **F-4** | **PLAUSIBLE, low.** The re-key is the right call but silently drops four liquid-phase items. `heat_effects.py:259` branches on `"(l)" in substance_name`; those four now take the gas branch, so the liquid sensible-heat pool collapses from six substances to two. The comment preserves *key order* for pool identity — but changing the *keys* moves the same seeds onto a different phase branch, which C2.3 needs to state. | **FIXED** — stated in the summary §11 and §6 |

---

## 4. Falsification attempts that failed

- **I tried to break the transcription.** Six species fetched live, every Shomate
  coefficient compared. Zero discrepancies. The JSON is not doing hidden work.
- **I tried to break the 2326 K argument.** I re-solved methane/air AFT from
  scratch — CH₄ + 2O₂ + 7.52N₂ → CO₂ + 2H₂O(g) + 7.52N₂, ΔHc = −802.31 kJ/mol,
  bisection on the product-enthalpy integral — twice: once with the committed
  table, once with **live NIST Shomate only**. Result: **repo 2311.1 K, NIST
  2327.2 K.** The claimed reference reproduces to within 1 K. **The argument is
  sound, not a rationalisation.** The spec's C2.4 is what is wrong here, not the
  implementation. The residual 16 K is *fully explained* by F-2 — the repo's Cp
  runs high at flame temperature, which depresses T. Fix F-2 and the number
  moves onto 2326 K by itself.
- **I tried to break the refits.** Both hold across the range the templates use:
  `C6H14(g)` max 2.79% over 298–1500 K (≤1.8% above 400 K); `C2H2(g)` ≤0.9% over
  400–1000 K. Only the 298 K endpoints are weak (+2.8%, +3.1%), and 298 K is
  precisely where these tables are most used.
- **I tried to make the propane ΔHf a defect.** It is not — NIST lists −103.8 as
  one of two measured values.
- **The suite runs green**: 208 checks, all pass.

---

## 5. Further probing and improvements

Triage in [`../phaseC2_summary.md`](../phaseC2_summary.md) §11.

**Thinnest evidence.** (a) `Air(g)` is self-verified against a mole-weighted
average — internal consistency, not an external source, and the assumed
composition is not recorded. (b) `NaCl(s)` at −5.5% is called "a genuine source
disagreement" with no argument offered; NIST *does* publish a solid-phase
Shomate, so it is resolvable and should not be waved through. (c) Five of the
eight condensed-phase rows the interim doc said needed sourcing were resolved by
re-keying rather than by sourcing.

**With more time I would** extend C2.2 to each template's true maximum
temperature and re-tune the tolerance per range (F-2 is the highest-value fix and
probably shifts several flame temperatures by 10–30 K); check the other seven
combustion reactions' flame temperatures the way I checked methane — Acetylene at
2844 K sits deepest into the extrapolated region; and verify the 11 reaction
stoichiometries balance atomically, which nothing in this phase appears to test.

**Does this generalise to C3?** Yes, two ways. First, the root cause was
*transcribing a column heading as part of the value*. Any table lifted from a
handbook with scaled column headings is exposed — Bird–Stewart–Lightfoot App. E
and Poling's Lennard-Jones tables (`GAS_MOLECULAR_PARAMS`, ε/k) use exactly that
convention, as do Perry's viscosity correlations. **C3 should grep for
`E-2`/`E-5`/`E5` literals as a first pass.** Second, and more important: **the
failure mode this phase actually exhibits is not bad values but unrecorded
validity ranges.** F-2 is a correct table used outside its domain.
`POWER_LAW_FLUIDS`, `REAL_FLUID_DATA` and the mechanical `MATERIAL_PROPERTIES`
all have implicit ranges that a value-by-value check will pass and a template
will then violate. **C3.2's consistency checks will not catch that class at
all** — it needs a per-table declared validity domain and an assertion that
consuming templates sample inside it.

**What the spec got wrong.** (1) C2.4's "~2200 K" is the *with-dissociation*
flame temperature and the wrong criterion for a complete-combustion template; my
independent solve puts the right target at 2327 K. The implementer was right to
deviate but changed a written exit gate on their own authority — that needs R4
triage, not a code comment. (2) The spec sized C2 at 15–25 h on one defect class;
there are four. (3) The spec assumed Smith–Van Ness would be citable; it is not
fetchable, and the substitution of NIST is well-reasoned and openly recorded —
but the exit-gate wording "citation to **edition + page**" is now unsatisfiable
by construction and should be reworded to "citation to a retrievable on-disk
artefact" before C3 inherits the same impossible gate.

---

*Note on §4: the reviewer's `COMBUSTION_REACTIONS` stoichiometry remark is the
one place the review is mistaken — the C2.2 suite does balance all 11 reactions
element-by-element and check the N₂/O₂ ratio. That check existed before the
review and passes. Recorded here rather than silently dropped.*
