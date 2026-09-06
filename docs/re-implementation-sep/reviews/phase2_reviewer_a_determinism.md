# Phase 2 — Reviewer A: Determinism

**Filed:** 2026-09-06 · **Frozen ref reviewed:** `76a66ff`
**Roster:** A — Determinism, phases 1 and 2
**Brief:** one mandatory gate task — *can a recorded seed regenerate every one
of these four items, byte for byte, under conditions the implementer did not
test?* 40-minute box.

> Worked read-only against a worktree pinned to `76a66ff`, never the live tree.
> The branch advanced to `3be62e3` during the review window; the diff was
> docs-only, and A noted that the frozen-ref discipline is what made that a
> footnote rather than a lost afternoon. That discipline exists because Phase C2
> lost a quarter of a reviewer's findings to a moving branch.

**Remediation is triaged in [`../phase2_summary.md`](../phase2_summary.md) §11.
All findings are closed.**

---

## 1. Verdict

**PASS WITH FINDINGS.** The gate passes: a recorded seed regenerates all four
templates byte for byte across every condition A could construct — two CPython
versions, three `PYTHONHASHSEED` settings, reversed seed order, absent numpy
seeding, and a seed range 500× beyond anything tested. **2,000 seeds per
template × 9 process-isolated configurations, zero diffs, zero errors.**

---

## 2. The gate, and the four conditions I did not test

| # | Interpreter | `PYTHONHASHSEED` | Seeds | Result |
|---|---|---|---|---|
| baseline | 3.13.13 | 0 | 0–999 | — |
| 1–3 | 3.13.13 / 3.12.13 | 0 and 1 | 0–999 | identical |
| 4 | 3.13.13 | **unset (randomised)** | 0–999 | identical |
| 5 | 3.13.13 | 0 | **999→0 reversed** | identical |
| 6 | 3.13.13 | 0, **no numpy seeding** | 0–999 | identical |
| 7/8 | 3.13 hs0 vs 3.12 hs1 | | **500000–500999** | identical |

Four of those are conditions I did not run and should have: **reversed seed
order** (rules out cross-seed generator-state leakage), **omitting the
`np.random.seed` my own dump script performs** (a template quietly using numpy
would look deterministic under my harness and fail under the real one),
a **randomised hash seed**, and a **seed range far outside 0–199**. My 200-seed
claim was true but thin; A's 1000-seed, 9-configuration result is what actually
establishes the gate.

A also confirmed: zero `except` of any kind in the four templates or the
constants modules they import; zero answer-less returns; numpy genuinely absent
from `conversion_and_reactor_sizing.py`; the one `set` reaching output is
`sorted()` first; T3 corpus-wide 0 failing.

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **A-1** | **CONFIRMED exposure, medium.** `pfr_volume_changing_rate` was the only one of the four with **no display-tie guard**, and its `rate_coefficient = _as_printed(k * (C_A0 ** n), spec)` rounds an **un-pre-rounded libm `pow()`**. C99 does not require `pow` to be correctly rounded. **1,170 of 749,111 reachable cells (0.156%) put `k·C_A0ⁿ` exactly on a 4-dp half-way point** — at `k=0.05, C_A0=1.21, n=1.5` the double sits *on* the tie and one ulp flips the printed `0.0665` to `0.0666`, taking the answer with it. A was explicit that it could **not** demonstrate a divergence: both interpreters here link MSVC v1944, so the 9-configuration byte-identity result is silent on this. Exposure, not an observed failure — and said so. | **FIXED** — bounded resample loop with a new `_display_is_fragile` guard. Reproduced A's 1170/749,111 exactly before fixing |
| **A-2** | **CONFIRMED, low (bookkeeping).** "Eight fallbacks removed" does not reconcile: diffing the four files gives **four `except` handlers**. | **FIXED** — a real inconsistency, and mine. The commit message said "eight"; the item-pool document's own table lists **ten** and is the correct count. A counted `except` statements; the table counts *fallback branches*, which includes two `abs(x)+1.0` repairs, a dead `unusual_case` branch, and two answer-less returns. The number is now stated once, as ten, with the accounting shown |
| **A-3** | **CONFIRMED (favourable).** libm exposure is confined to the PFR template; the other three are provably portable. `heat_effects` uses only `**2`/`**3` on values already bound through `_as_printed`; the isolator uses `math.sqrt`, which IEEE-754 *requires* to be correctly rounded. | **RECORDED** — it is what scopes A-1 to one template |

**A corrected itself mid-review**, and the correction is the useful part: it
first flagged `integral_term` as the finding, then traced the rounding chain and
found that pre-rounding `power_term` through `_as_printed` insulates it — the
derived value is a division of two clean doubles and is bit-identical
everywhere, **even on the 6 of 396 cells that sit exactly on a tie**.
`rate_coefficient` is the one place that pre-rounding was missing. That is
precisely why A-1 is the finding and the other is not.

---

## 4. Falsification attempts that failed

The load-bearing section. Each is something A expected might break:

- **Cross-version identity at 5× the claimed sample**, then again at seeds
  500000–500999. Identical.
- **Cross-seed state leakage** — seeds run 999→0. Identical to the forward run.
- **Hidden numpy dependence** — a run seeding `random` only. Identical, which
  also confirms the harness convention in `instance_dump.py` is safe here.
- **The isolator algebra, including the part my prose skipped.** My argument was
  over *exact* values; the code operates on `_hu(..., 4)` half-up-rounded ones,
  so the proof does not automatically transfer. A checked that it survives:
  `A = TR² ∈ [0.0025, 0.04]` and `C ∈ [−0.9975, −0.96]` never round to zero, so
  `−4AC ≥ 0.0096` and the discriminant is **~96 units in the last printed place**
  from zero — not marginally positive. `|sol2| ≈ 1` across the range, so rounding
  cannot collapse it and trip `len(admissible) != 1`. **The claim holds for the
  whole sampled range including the rounding the code performs** — which is more
  than I had established.
- **Isolator resample exhaustion** — over 20,000 seeds the worst case used
  **2 passes of 200**. 100× headroom.
- **Flame fixed iteration count** — a literal, no convergence test in the loop,
  every iterate bound before feeding the next.

---

## 5. Further probing

**Thinnest evidence.** A-1 is untestable on this machine: both interpreters link
the same MSVC runtime. The decisive experiment is one run on Linux/glibc — ten
minutes in a container. *Since the fix removes the fragile draws entirely rather
than arguing about which libm is right, that experiment is no longer on the
critical path, but it remains the way to confirm the exposure was real.*

**The "5000 seeds at 0.00%" measurement A took entirely on trust** — there is no
artefact in the frozen commit to check it against. Fair. The measurement script
is `scratchpad`-only and was not committed; the rates are recorded in
`phase2_item_pool_impact.md` §4 but the *evidence* for them is not reproducible
from the repo. Noted as a gap, not fixed here.

**Generalisation.** A's sharpest structural observation: **two rounding
conventions now coexist in the corpus** — `_hu` (Decimal `ROUND_HALF_UP`, per
D-012) in the vibrations file, and `_as_printed`/`format` (binary
`ROUND_HALF_EVEN`) in the three chemical-engineering files — and only the
vibrations file carried a tie guard. That unevenness is what A-1 is an instance
of. A recommends running its tie-margin enumeration over every template
combining a transcendental with a fixed-decimal print, and expects the PFR
template not to be the only one. **Carried forward to Phase 3** as an open item;
it is a corpus-wide sweep, not a Phase 2 fix.
