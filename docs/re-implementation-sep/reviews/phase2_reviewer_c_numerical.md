# Phase 2 — Reviewer C: Numerical methods

**Filed:** 2026-09-06 · **Frozen ref reviewed:** `76a66ff`
**Roster:** C — Numerical methods, phases 2 and 3
**Brief:** one mandatory gate task — *does the replaced solver converge to the
right answer across the full sampled parameter range, and is every printed
intermediate reproducible by hand?* 45-minute box, plus a post-box addendum.

> Worked read-only. NIST Shomate coefficients were fetched **live** from
> `webbook.nist.gov` and parsed from raw HTML rather than read from the
> implementer's on-disk JSON, so the verification is independent of the
> transcription it is checking.

**Remediation is triaged in [`../phase2_summary.md`](../phase2_summary.md) §12.
All findings are closed.**

---

## 1. Verdict

**PASS WITH FINDINGS.** Both numerical replacements are sound. The
mean-heat-capacity iteration converges to the correct root for all eleven
reactions; six passes is comfortably enough (worst deviation from the exact
fixed point **0.035 K**, against my claimed 0.06 K); rounding the divisor each
pass does **not** produce a limit cycle; the map is a strong local contraction
(|g′| ≈ 0.09–0.135) so the 2000 K guess is safe across a 400–10000 K basin; and
**no reaction lands on the wrong side of a 0.5 K rounding boundary** — tightest
margin 0.116 K against a 0.035 K error, a 3.3× safety factor. The methane
instance reproduces by hand from the printed numbers alone. The closed form is
correct in sign, has the right n→1 logarithmic limit, and agrees with C's own
quadrature to **2.673e-15**, matching my 2.7e-15 to three digits.

---

## 2. Independent re-derivation

C re-derived rather than checked, which is what makes this worth having:

| Quantity | Committed | C derived | Source |
|---|---|---|---|
| Shomate A–H, 4 species, 10 ranges | on-disk JSON | **identical, zero mismatches** | live `webbook.nist.gov` |
| T_ad, all 11, from NIST | "within 1.9 K" | worst **1.73 K** | own multi-range solve + brentq |
| 6-pass vs exact fixed point | "within 0.06 K" | worst **0.035 K** | own brentq |
| Cp refit coefficients | 16 values | reproduced to **≤0.21%** | own `lstsq`, refit from scratch |
| Worst Cp residual at 298 K | −5.43% | **−5.45%** | own residuals |
| Cp residual above 1000 K | ≤1.30% | **1.30%** exactly | own residuals |
| PFR closed form vs quad | 2.7e-15 | **2.673e-15** | own `quad` |
| n→1 limit | asserted valid | → −ln(1−X) to 8 digits | own sweep |

C's solve passes its own zero-at-298.15 K and range-continuity checks
(≤3.2e-3 kJ/mol across every join) before being used as an arbiter.

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **C-1** | **CONFIRMED, low.** `constants.py` claimed agreement "within **1.7 K**" while the commit message said 1.9 K. Two numbers for one claim, and the stricter is false: carbon monoxide is 1.93 K out. | **FIXED** — re-measured template-vs-NIST across all 11; worst is **1.93 K** (CO). The 1.7 K was the *refit*-vs-NIST residual, a different and smaller quantity, and the file now says which is which |
| **C-2** | **CONFIRMED, medium.** The answer is quoted to the nearest kelvin but the Cp model is only good to ~2 K: the nearest-kelvin value differs from a live-NIST solve on **6 of 11** reactions. "A grader marking to ±1 K would fail a student who used NIST data." | **FIXED** — Step 6 now states the model tolerance in the trace itself: *"the model is good to about ±2 K against NIST reference data, so the answer should be read as N ± 2 K."* The gold trace declares its own precision rather than implying four significant figures it cannot support |
| **C-3** | **CONFIRMED, medium — the one C would block a merge on.** `f"k = {k} L^{n-1}/..."` interpolates `n-1` raw, printing `L^0.6000000000000001` on **38.6% of instances**, twice each. Step 3 twenty lines below already rounds it. | **FIXED** — bound as `n_minus_1`/`one_minus_n`. Reproduced at 771/2000 seeds before fixing |
| **C-4** | **CONFIRMED, low–medium.** `PFR_DP = 4` is **decimal** places on a quantity spanning 2.5 decades, leaving two significant figures at the bottom. C's post-box **exhaustive** grid — all 26,967,996 sampled points — puts the worst answer-key error at **0.5122%**, with **1.17% of the space above 0.1%**. Not an outlier: a population. | **FIXED** — switched to six **significant figures**. Worst answer error falls to **0.00123%**, a 417× improvement, and the printed answers now match the exact closed form on every seed that draws the same parameters |
| **C-5** | **CONFIRMED, low.** The docstring said the template "integrates to 2844 K"; the true maximum is **acetylene at 2908 K**. | **FIXED** — 2844 K was measured with the *old* table, so correcting the table raised the temperatures and the stale figure understated the very margin it described. Now 2908 K, with the 92 K of headroom stated |
| **C-6** | **CONFIRMED, low.** The docstring said the fixed count "needs no convergence test at run time" — but line 456 is exactly that assert. | **FIXED** — the docstring stopped denying it. C is right that the assert should stay, and right that `python -O` strips it, so it guards development rather than production. Both now said |
| **C-7** | **CONFIRMED, low.** The `H2O(g)` `[DERIVED]` tag claims a refit "over 298–3000 K", but NIST's lowest gas-phase Shomate range for water **starts at 500 K**, so the 298–500 K part of the grid is extrapolated Shomate and the quoted +0.87% is measured against an extrapolation. | **FIXED** — caveat added to the tag. C checked the damage itself and found it negligible (extrapolated Cp(298.15) = 33.590 vs JANAF 33.58, 0.03%), so the value stands and only the claim needed narrowing |

---

## 4. Falsification attempts that failed

Every one of these is a claim I made that C tried to break and could not:

- **"Six iterations is not enough at the extremes."** All 11 run to a 400-pass
  rounded attractor: worst 6-pass deviation **0.035 K**.
- **"Rounding the divisor produces a limit cycle."** All eleven reach a
  **period-1** attractor by pass 7 at the latest.
- **"The fixed count lands on the wrong side of a 0.5 K boundary."** Zero flips
  in 11/11; thinnest margin propane at **0.116 K**.
- **"The 2000 K guess is unsafe."** Swept 400 → 50000 K; every guess from 400 to
  10000 K gives the identical answer for all eleven. C also found the map is
  **not** a global contraction (|g′| reaches 1.083 for CO somewhere in
  400–6000 K) — the basin merely contains the stated guess comfortably. That is
  a sharper and more honest statement than mine.
- **"The closed form has a sign error or breaks near n=1."** Correct, positive
  throughout, → −ln(1−X) verified to 8 digits at n = 1+10⁻⁸.
- **"The 2191 K contamination survived into the final numbers."** Could not find
  it. C's own multi-range solve gives methane 2327.19 K. *This was the specific
  thing I asked C to look for, having nearly decided D-036 on that wrong number.*
- **Float-display artefacts in the flame trace.** Scanned all 11: none. The
  defect is confined to the PFR (C-3).

---

## 5. Further probing

**Thinnest evidence.** C's own arbiter shares `HEATS_OF_FORMATION` with the
template by construction — correctly, since that isolates the Cp question and
heats of formation were fenced off as settled by Phase C2. A genuinely external
check would be Cantera or a JANAF equilibrium solve. Second, the 0.116 K
rounding-boundary margin is adequate but is a property of these eleven
reactions, not a guarantee.

**With more time C would** solve all eleven in Cantera as a third independent
path, and check whether the 92 K of headroom survives any *future* reaction
added to `COMBUSTION_REACTIONS` — acetylene is already at 2908 K and a hotter
fuel would trip the assert. **Recorded as an open item:** the assert is the
guard, and it is correct that it would fire rather than silently extrapolate.

**Does the failure mode generalise? C's answer is the most valuable thing in
this review, and it partly contradicts my own §10.**

> This template is safe because its parameter space is a **closed set of eleven
> discrete items** — no continuous random inputs — so the rounding-boundary
> question could be verified *exhaustively*. `normal_depth_iteration` will not
> have that property.

Normal depth is a solve over continuously sampled geometry, slope and roughness,
so a fixed iteration count **cannot** be exhaustively validated and there will be
a measure-zero-but-nonempty set of parameters where the fixed count lands on the
wrong side of a display boundary. It needs either a proven a-priori error bound
over its sampled box, or a run-time convergence assert **not stripped by `-O`**.

And `line_balancing_heuristic` is a different problem entirely: **a heuristic has
no fixed point, so "six passes" reasoning does not transfer at all.** The
question there is whether the heuristic is deterministic under tie-breaking —
Reviewer A's territory, not a convergence question.

My §10 offered Phase 3 the criterion *"fix the iteration count when the count is
not part of what the item tests."* C's finding sharpens it: that criterion is
necessary but **not sufficient**. It also requires a parameter space you can
verify over. Phase 3 should treat its two templates as two different problems,
and neither of them is this one.

C-3 and C-4 generalise immediately: any template printing a raw arithmetic
expression like `{n-1}` in an f-string, or applying fixed decimal places to a
quantity spanning decades, has the same defect. C suggests a corpus-wide regex
for `\{[a-z_]+ *[-+*/] *[0-9.]+\}` inside f-strings. **Carried to Phase 3.**
