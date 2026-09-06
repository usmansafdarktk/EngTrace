# Phase 2 — item-pool impact (D2.3, D2.5)

**Date:** 2026-09-06 · **Branch:** `redesign/phase2-determinism`
**Method:** `(question, solution)` dumped for seeds 0–299 from a `master`
worktree and from the branch, **each in its own process**. An in-process reload
resolves both sides to the same already-imported modules and reports everything
identical — the trap `tests/template_integrity/instance_dump.py` documents, and
which has already caught someone on this project.

---

## 1. Summary

| Template | Question changes | Answer changes | Median shift | Max shift |
|---|---:|---:|---:|---:|
| `levenspiel_plot_interpretation` | **300/300** | **300/300** | — (whole pool) | — |
| `pfr_volume_changing_rate` | 300/300 (wording) | **24/300** — all of them a *different draw*, not a moved answer | — | — |
| `adiabatic_flame_temperature` | 0/300 | **300/300** | 0.845% | 2.256% |
| `vibration_isolator_design` | 0/300 | **300/300** | 0.012% | 0.128% |

Zero generation errors on either side, all four templates, 300 seeds.

**The PFR row changed after review** and is worth reading carefully. It first
measured 79 answers moved by a median 0.041% — that movement was the **4-decimal
rounding error Reviewer C found in C-4**, not something intended. Switching to
six significant figures removed it: of the 276 seeds that draw the same
parameters as master, **the answer moved on none**. The remaining 24 (8.0%) draw
*different parameters* entirely, because Reviewer A's A-1 tie guard resamples
0.057 draws per instance and a rejected draw consumes randomness. Those are
different items, not different answers to the same item.

The lesson is worth keeping: **the first measurement was of my own defect, and
it looked small enough to accept.** A 0.041% median reads as rounding noise. It
was a 0.5% worst-case error in the answer key, and only an exhaustive grid over
all 26,967,996 sampled points showed that 1.17% of the space exceeded 0.1%.

---

## 2. Why each pool moved

### `levenspiel_plot_interpretation` — the whole pool, unavoidably

This template sampled its measurement noise from `np.random`, whose generator
`random.seed()` does not touch. **The recorded seed never regenerated the item
in the first place**: the production path (`main()`) seeds only `random`, so the
published items for this template cannot be reproduced from their seeds at all.
Any fix therefore changes the pool — there is no "before" to preserve, because
"before" was not a function of the seed.

Sampling is now stdlib `random` throughout and numpy is gone from the file, so
one seed determines the item. The table's *character* is unchanged: same
conversion grid construction, same rising-rate shape, same 74%-of-draws
monotonicity enforcement, same distributions for `F_A0` and target conversion.
What changed is which draw a given seed produces.

### `pfr_volume_changing_rate` — every question, and 8% a different draw

The **question** changed on all 300 by wording alone: *"An {n}-order …"* became
*"An order-{n} …"*, and the closing note claiming analytical solutions are
"complex" was dropped. That note was the reason `n = 2.0` (9.6% of draws) was a
defect: the question promised a *non-integer* order and then sometimes gave an
integer one. Rewording rather than resampling keeps the item pool's *parameters*
identical — `n` still draws from the same 11 values — which is the cheaper fix
under P6 (D-016's principle: change the claim, not the sample, when the claim is
what is wrong).

The **answer** is unchanged on every seed that draws the same parameters —
276 of 300, moved on none. Intermediates are bound through their display, so
`V = (F_A0 / 0.227294) × 2.46757` is computed from the numbers the reader can
see; at six significant figures that binding costs at most 0.00123% of the
answer, which is invisible at the 2 dp the answer is quoted to. The closed form
moved nothing either: it agrees with `quad` to 2.7e-15 over 3000 seeds.

The other **24 of 300 (8.0%) draw different parameters**, because the A-1 tie
guard rejects 0.057 draws per instance and a rejected draw consumes randomness.
Those seeds produce a different — equally valid — item, not a different answer
to the same one. That is the cost of removing the display ties, and it is the
same cost Phase 1 paid wherever it resampled.

### `adiabatic_flame_temperature` — every answer, and every one of them closer to reality

Questions are untouched. Every answer moves, by a median 0.845% and at most
2.256%, for two compounding reasons:

1. **`CP_PARAMS_COMBUSTION`** (D-036). The old answers were computed from a
   polynomial fitted to 1500 K and evaluated to 2908 K. They were all low.
2. **The stated algorithm.** The answer is now the sixth iterate of a printed
   iteration rather than a root found by `fsolve` to its own tolerance, and it
   is quoted to the nearest kelvin rather than 2 dp — spurious precision on a
   model that excludes dissociation.

The movement is a **correction, not a drift**. Against a direct NIST Shomate
solve of the same energy balance:

| fuel | before | after | NIST | after − NIST |
|---|---:|---:|---:|---:|
| Methane | 2311.2 | **2327** | 2327.2 | −0.2 |
| Propane | 2373.9 | **2394** | 2394.1 | −0.1 |
| Hydrogen | 2494.2 | **2524** | 2524.5 | −0.5 |
| Carbon monoxide | 2623.2 | **2662** | 2663.9 | −1.9 |
| Acetylene | 2843.8 | **2908** | 2908.6 | −0.6 |
| Ethane | 2361.4 | **2380** | 2380.7 | −0.7 |
| Butane | 2378.9 | **2399** | 2399.3 | −0.3 |
| Octane | 2388.5 | **2409** | 2409.6 | −0.6 |
| Ammonia | 2101.5 | **2108** | 2106.8 | +1.2 |
| Methanol | 2315.8 | **2332** | 2332.3 | −0.3 |
| Ethanol | 2335.8 | **2353** | 2353.4 | −0.4 |

Methane against the independent complete-combustion reference of **2326.35 K**
(`docs/references/README.md`): **2327 K, 0.03%.** Before: 2311 K, 0.65%.

**These items were previously scored against answers that were wrong by 5–65 K.**

### `vibration_isolator_design` — every answer, mostly in the last digit

Questions are untouched. Every answer moves by a median 0.012% — the
round-then-recompute binding again. This template failed T1 on **272 of 300**
instances: `omega_n = omega / r` was computed at full precision and printed from
rounded operands, so a reader following the printed numbers did not reach the
printed answer. The answer is now what the trace derives.

Separately, the answer **format** changed on all 300: `{:,.0f}` emitted
`3,733,620 N/m`, and the thousands separators break the evaluation parser (the
D-003 defect that reads `**4,921**` as `4.0`). It now prints `3734321 N/m`.
Every one of these 300 items was previously unparseable by the scorer.

---

## 3. T6 (distribution non-regression) — cannot be measured against the baseline

**T6 fails on 142 of 150 templates at `master`.** Verified in a worktree, not
assumed. The committed baseline in `tests/template_integrity/baseline/` is stale
corpus-wide and predates work that has since landed, so it cannot distinguish
this phase's four templates from the 138 others that fail against it.

D2.5 asks for a T6 distribution diff. It is delivered above as a **direct
before/after instance diff** instead, which is strictly more informative for
these four: it reports what changed per seed rather than whether a summary
statistic drifted.

**`SPEC-CHANGE` / finding, raised not fixed:** the T6 baseline needs
regenerating, and regenerating it is not this phase's to do — a baseline
refreshed by the phase it is meant to gate is not a gate. It should be
regenerated deliberately, on `master`, as its own change, with the movement
since it was last written examined rather than absorbed. Recorded here and in
the phase summary.

---

## 4. Removed fallbacks (D2.3) — measured before removal

Rule from D-024/D-026: a fallback you remove is either **provably unreachable**
(show the argument) or it **was changing answers** (report the rate). Measured
over **5000 seeds**, which resolves a rate down to ~0.06%.

| Template | Fallback | Fired | Why removal is safe |
|---|---|---:|---|
| levenspiel | `except Exception` around the CSTR volume | 0.00% | Guarded `V_CSTR <= 0`; `F_A0`, `y`, `X` are all positive by construction and asserted |
| levenspiel | `except Exception` around the PFR volume | 0.00% | Guarded a length mismatch, non-positive rates, non-monotonic X, and `V_PFR > 10·V_CSTR` — the first three are now asserts, the fourth is impossible (see below) |
| levenspiel | `except Exception` around the trapezoid loop | 0.00% | Guarded `dx <= 0 or avg_height <= 0`; both are asserted positive |
| levenspiel | `V_CSTR <= 0 → abs(V_CSTR) + 1.0` | 0.00% | A "fix" that would have silently invented a volume |
| levenspiel | `V_PFR <= 0 → abs(V_PFR) + 1.0` | 0.00% | Same |
| levenspiel | `unusual_case` (PFR > 1.5·CSTR) | 0.00% | **Provably unreachable.** `1/(-r_A)` is enforced non-decreasing, so the exit rectangle contains the area under the curve and `V_PFR ≤ V_CSTR` always. Now an assert with that argument in its message |
| flame | `except (KeyError, ValueError, IndexError, Exception)` in a retry loop | never retried in 2000 seeds | It caught the solver's own failures and silently asked a **different question**. The data preconditions are now explicit asserts |
| flame | `while/else` returning an answer-less tuple | 0.00% | Returned `("Error: …", "Possible causes: …")` — not a benchmark item at all |
| isolator | `discriminant < 0 →` answer-less return | 0.00% | **Provably unreachable.** `TR ≤ 0.20 < 1`, so `C = TR² − 1 < 0` and `A = TR² > 0`, hence `−4AC > 0` and `D = B² − 4AC > 0` always |
| isolator | `r² < 0 →` answer-less return | 0.00% | **Provably unreachable.** The roots' product is `C/A < 0`, so exactly one root is positive. This is now the *stated selection rule* rather than `max()` |

Every one was dead. That is the point worth keeping: **none of them was load
bearing, and each would have produced a confidently wrong gold trace on the day
it fired.** The isolator's `max()` root selection is the sharpest case — it
happened to pick the right root for the wrong reason, and "take the bigger one"
is not a rule a solver could state or a verifier could check.
