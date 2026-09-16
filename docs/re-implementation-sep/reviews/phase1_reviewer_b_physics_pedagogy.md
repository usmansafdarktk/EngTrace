# Phase 1 — Reviewer B: Physics and Pedagogy (the P6 guard)

**Filed:** 2026-09-06 · **Branch reviewed:** `redesign/phase1-round-trip`
**Roster:** B — Physics & pedagogy (spec §Roster)
**Brief:** one mandatory gate task — *for each edited template, would a domain
expert notice this item got easier, less interesting, or less physically
sensible?* 35-minute box.

> The reviewer worked in isolation from Reviewer A and without sight of the
> implementer's reasoning, per R1.1. Code comments and commit messages were
> supplied as **claims under test**, not as evidence. Arithmetic closure,
> round-trip correctness, oracle independence and determinism were explicitly
> out of scope — Reviewer A owned all of it, in parallel.

**Implementer's remediation is recorded in the right-hand column of §3 and in
[`../phase1_summary.md`](../track-a/phase1_summary.md) §7. All CONFIRMED findings are
closed.**

---

## 1. Verdict

**BLOCKED** — narrowly: one CONFIRMED pedagogical regression
(`beam_deflection_formula`, US-customary branch) and two unrecorded item-content
changes, including the annulus sign-off the brief says exists. No physics error,
no difficulty collapse, and no distribution regression that could be
substantiated; the remedies are small.

---

## 2. Independent re-derivation

**Scope check.** Six files are edited, containing 24 templates. Every template
in every edited file was generated at seeds 0–5 under `master` and under the
branch and compared as tuples. **Exactly 10 templates change output.** The other
14 — `shear_stress_torsion`, `angle_of_twist`, `system_properties`,
`vibration_isolator_design`, `logarithmic_decrement`,
`equivalent_stiffness_frequency`, the four signal templates and the rest — are
byte-identical. **No collateral damage**, which was a real risk given the files
were largely re-emitted.

**T6, all ten, 1,000 seeds vs the committed baseline.** Step-count distributions
are **identical before and after for all ten** (`beam` 520/480 3-step/4-step;
`damping` 678/322; the rest single-valued). Distinct-answer counts rose or held
everywhere. Two breaches: `vibration_transmissibility` −2.3% median,
`annulus_flowrate` −28.7% median.

**`annulus_flowrate` — the 28.7% breach.** Both samplers' draw sequences were
replicated over 1,000 seeds:

- `master`'s silent fallback (`if pressure_drop_Pa == 0: pressure_drop_Pa = 10.0`)
  fired on **172/1000 (17.2%)**, matching the in-code comment's 17.8% but **not**
  the "284 of 1,000" in the brief — that figure is the separate 283/300 Step-5
  non-closure count from `phase0_oracle_findings.md`. *The claim as briefed does
  not reproduce; the claim as coded does.*
- Independently and more seriously: on `master`, **171/1000 instances had
  Re ≥ 2100 computed from their own printed answer**. The item asserts the
  laminar annulus solution over transitional flow on ~17% of instances. The new
  `Re < 2100` guard removes that.
- New rejection cost: **25 rejected draws per 1,025** (18 for ΔP rounding to 0,
  7 for Re ≥ 2100) = 2.4%.
- The shift is **not** sampling noise (contrast D-019). At **5,000 seeds** per
  side: p5 2.68e-5→1.30e-5, **p50 2.40e-4→1.72e-4 (−28.1%)**, p95
  2.78e-1→2.85e-1 (+2.7%), distinct answers **4929→4940**. Entirely a downward
  extension of the low-Q tail — the signature of removing a fallback that pinned
  17% of instances to ΔP = 10 Pa when the sampler asked for 1–4 Pa. The answer
  space is not narrowed and Q stays in scientific notation, so nothing about the
  arithmetic burden changes.

**`damping_classification` — D-020's separation claim.** ζ recomputed in 50-digit
decimal from the *stated* m, k, c over 5,000 seeds, both branches:

```
OLD critical n=1676 max|z-1|=2.054e-05 median=2.172e-07 | other n=3324 min|z-1|=0.1500
NEW critical n=1676 max|z-1|=2.054e-05 median=2.172e-07 | other n=3324 min|z-1|=0.1500
```

D-020's numbers reproduce to three digits and — the point D-020 does not make —
**they are identical on `master`.** Label mix is 1652/1676/1672
(Under/Critical/Over) on both sides. The branch changes nothing about
answerability here; it merges two display lines and quotes `c_c` at 2 dp and ζ
at 4 dp.

**Beam/cantilever precision.** The metre display drops 5 dp → 4 dp
(`0.01120`→`0.0112`). 0.0001 m **is** 0.1 mm exactly, so the mm answer is
unchanged in every instance checked and no resolution is lost.

**Physical plausibility**, 300 instances each: RU amplitude 0.15–20.4 mm; shaft
diameter 9.7–88.2 mm; SIS T_A 167–3814 N·m; TR 0.056–4.90. All
engineering-sensible. Two absurdities exist — annulus ΔP up to 11.7 GPa and
composite-shaft twist to 185.6° — and **both reproduce identically on `master`**
(CSS old max 185.5944° vs new 185.5942°; annulus seed 19 old 11686663.280 kPa vs
new 11686787.580 kPa). Pre-existing.

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **F-1** | **CONFIRMED, blocking.** `beam_deflection_formula`'s US-customary Step 2 no longer converts anything. `w = 1.48 kip/ft = 1.48/12 kip/in`, with the substitution carrying the exact ratio. As hand-calc notation this is unimpeachable; as a **gold reasoning trace** it is a degradation. (a) The step is titled "Assemble consistent US customary units" and its other three lines are restatements of givens, so `w` was the only quantity it produced — the step is now a no-op. (b) kip/ft→kip/in is the single most error-prone operation in the US branch and the trace no longer demonstrates it. (c) It is asymmetric: the SI branch still evaluates every conversion numerically, so the same item tests different things depending on which unit system it samples. **This is a correctness fix that removed pedagogy as a side effect, the exact shape P6 prohibits.** | **FIXED** — Step 2 now performs the conversion and states its value to 6 dp (`w = 1.48 kip/ft = 1.48/12 = 0.123333 kip/in`); the substitution still carries the exact ratio. |
| **F-2** | **CONFIRMED, blocking.** The `annulus_flowrate` T6 breach has **no written sign-off**, and neither do two other item-content changes. `grep -n "annulus\|28.7" DECISIONS.md` returns nothing in D-016…D-021. D-019 signs off the *smaller, non-real* `vibration_transmissibility` move and D-021 signs off `shaft_design_power`; the 28.7% one — the only genuine distributional relocation in the phase — is unrecorded. Also unrecorded: `rotating_unbalance`'s question gained "(this total includes the rotating unbalanced mass)", which is physically correct for the Rao formulation used and removes a real ambiguity, **but removes a modelling judgement from the solver's task** — a P6 scoping decision. **If this review is the intended sign-off vehicle: I sign off the annulus distribution change on the evidence in §2.** It must still be written down. | **FIXED — D-022**, recording all three with this reviewer's sign-off. |
| **F-3** | **CONFIRMED, minor.** `mean_variance`'s variance section states the answer before the working: line 1 jumps substitution→answer, line 2 backfills the intermediates and repeats the same total. Precision is also internally inconsistent — mean terms at 3 dp, variance terms at 6 dp. Content and difficulty are preserved; the trace reads as sloppy and is a poor exemplar of "show your work". The underlying numeric change (13.889→13.890, from consuming the *stated* 3-dp probabilities) is a correct P3 fix and improves the item. | **FIXED** — deviations moved to their own line, main line runs deviations → products → total once. The precision difference is retained: each quantity is shown at the precision it exactly has. |
| **F-4** | **PLAUSIBLE, minor.** The 1 Pa grid pushes the stated ΔP to 11 significant figures on viscous fluids (`11686787.580 kPa`). `master` already read 10 s.f., so this is a one-digit worsening of a pre-existing absurdity. At the other end the grid yields `0.001 kPa` givens, which would read better as `1 Pa`. Neither affects difficulty. | `BACKLOG` — folded into the plausibility-screening item. |
| **F-5** | **PLAUSIBLE, cosmetic.** `Pa·s` → `Pa.s`. Defensible against the cp1252 mojibake, but `Pa-s` or `Pa*s` would read better. The sibling vibrations templates already mix `N-s/m` and `N.s/m`. | `ADOPT-PHASE-5` — unit-symbol consistency is Phase 5's contract work. |
| **F-6** | **PLAUSIBLE, process.** D-021 establishes that a distinct-answer count above ~50% of N is an instrument artefact — and `annulus_flowrate` reads 998→1000 of 1,000, i.e. **99.8% saturated**, and was passed on that number. Re-measured at 5,000 (4929→4940) there is no hidden regression, but the gate as run did not establish that. | **ADOPTED** — re-measured at 5,000 seeds; generalised as D-021. |

---

## 4. Falsification attempts that failed

1. **"The resample loops remove a physically meaningful region."** They do not.
   Annulus rejects 25/1,025 draws: 18 where ΔP rounds below 1 Pa
   (arithmetically degenerate — the question can only state 3 dp of kPa) and 7
   where the *answer's* Reynolds number reaches 2100, which the governing
   equation forbids. **The loop removes strictly less physics than `master`'s
   fallback destroyed**, and it repairs 171/1,000 instances that were
   transitional-flow items wearing a laminar solution. Beam/cantilever
   tie-resample rates (0.02%/0.115%) are too small to shape anything.
2. **"D-020's separation is an artefact of the fix."** No — it reproduces to
   three digits *and is identical on `master`*, as is the label mix.
   **D-020 accepted.** A 7,300× gap between the worst critical instance
   (|ζ−1| = 2.05e-5) and the nearest non-critical one (0.150) means no solver
   applying any precision-aware rule can land on the wrong label. **D-004 should
   stay closed**: constraining k·m to a perfect square buys an exactness nothing
   consumes, at the cost of parameter values that look manufactured — a worse
   trade than the one it fixes. Reservation recorded, not raised: 33.5% of
   instances being exactly critically damped is an implausible population for a
   classification task, and it is unchanged from `master`.
3. **"The significant-figure moves made items easier or revealed structure."**
   The opposite. `composite_shafts_series` J at 6 s.f. and φ at 6 s.f.,
   `shaft_design_power` c at 6 dp, annulus κ and shape factor at 6 dp all
   *increase* the arithmetic burden and none reveals a shortcut.
   `composite_shafts_series` is strictly better: `master` quoted the same answer
   as 0.0901 rad (3 s.f.) and 5.1609° (5 s.f.) — two precisions for one quantity.
   `rotating_unbalance` moves the *other* way (6 s.f. → 4 s.f.) and its distinct
   count still rose 901→964, so no answer-space loss there either.
4. **"Dropping metres from 5 dp to 4 dp costs resolution."** It cannot: 1e-4 m
   is exactly 0.1 mm, the answer's own precision. Every mm answer compared is
   byte-identical.
5. **"The branch introduced physical absurdities."** The 11.7 GPa annulus and
   the 185.6° composite twist both reproduce on `master` at the same seeds.
6. **"Governing equations changed."** They did not, in any of the ten. Annulus
   keeps `Q = πΔP·R⁴/(8μL)·[(1−κ⁴) − (1−κ²)²/ln(1/κ)]`; RU keeps
   `X = mₑeω²/√((k−mω²)² + (cω)²)`; TR keeps
   `√[(1+(2ζr)²)/((1−r²)²+(2ζr)²)]`; beam keeps `5wL⁴/384EI` and `wL⁴/8EI`;
   torsion keeps `φ = TL/JG`, `c³ = 2T/πτ`, `Tᴀ(1+L₁J₂/L₂J₁) = T`. Where a
   *number* moved (RU ω 106.55926→106.60471, VT ω 187.336→187.3646, SIS Tᴀ
   3278.75→3279.94) **the new value is the one derivable from the stated givens
   and the old one was not** — these are P3 repairs, and every one makes the item
   *more* answerable.

---

## 5. Further probing and improvements

Triage in [`../phase1_summary.md`](../track-a/phase1_summary.md) §8.

**Thinnest evidence, ranked.**
1. Instances were read at ~8 seeds per template. **Systematic absurdity
   screening** (velocity, stress, deflection ratio, Re) across the full 1,000 is
   the next measurement to buy, and it belongs in the harness as a per-template
   plausibility band rather than in a reviewer's head. It is what would have
   caught the 11.7 GPa annulus and the 185° twist years ago.
2. Whether **raising displayed precision degrades grading** is unverified.
   `composite_shafts_series` now asks for 0.090074 rad; a solver who rounds
   intermediates sensibly will differ in the 5th s.f. That is a T2/verifier-
   tolerance question — **T6 cannot see it.**
3. `damping_classification`'s 1/3 critically-damped population deserves a look
   on its own merits in a later phase.

**Where the same weaknesses live outside this phase.** The unevaluated-fraction
pattern (F-1) will recur wherever a US-customary conversion sits on a display
tie — check `virtual_work_truss_deflection` and the rest of `structural_analysis`,
plus any civil template with kip/ft or psf loads. The saturated distinct-answer
count (F-6) affects **every** template whose answer space exceeds ~500, which
from the T6 profiles is most of the corpus. The `Pa·s` / `N-s/m` / `N.s/m`
unit-symbol inconsistency is corpus-wide.

**What later phases should do differently.**
- Make the P6 record a **deliverable, not a footnote**: any diff that changes
  question *text*, a step's content, or a displayed precision should require a
  named DECISIONS line before the phase can close. **Two of my three blocking
  findings are things nobody wrote down, and one of them the brief asserts
  *was* written down.**
- Give T6 a **question-text hash** alongside the answer profile.
  `rotating_unbalance`'s wording change is invisible to every gate in the phase,
  because T6 profiles answers and step counts, not what the solver is told.
- Run the P6 review against **before/after instance diffs at matched seeds**,
  not the source diff. Five of these findings came from the instance diff in
  under ten minutes; none is legible in the 4,000-line source diff.

**Spec defects.**
- **§0.2 T6 is ambiguous about what "sign-off" is.** D-019 signs off inside
  DECISIONS, D-021 inside DECISIONS, the annulus breach apparently inside a
  review brief. Name the artefact and the approver.
- **§0.2 T6's gate set does not include the question.** It gates distinct
  answers, quantiles, step counts, branch proportions and difficulty label —
  every one a property of the *solution*. **P6 is about what the item tests,
  which lives in the question.** This is a genuine gap in the P6 guard.
- **P6 says "recorded and approved" but names no approver and no threshold** for
  what counts as altering "what an item tests".

---

## P6 answer, per template

| Template | Would a domain expert notice it got easier, less interesting, or less physically sensible? |
|---|---|
| `mean_variance` | **No easier, slightly worse-presented.** Same 2 steps, same task. Numbers corrected (P3). F-3. |
| `rotating_unbalance` | **No — marginally easier, and better physics.** ω is now derivable from the stated RPM (it was not). The mass clarification removes a modelling ambiguity: a small easing, unrecorded — F-2. |
| `vibration_transmissibility` | **No.** Same 4 steps, same formula. ω now derivable from stated f. TR to 4 s.f. instead of 3; distinct answers 676→928. The 2.3% median move is sampling noise. |
| `damping_classification` | **No.** Purely presentational. Question, parameters, label mix and the four-orders-of-magnitude separation are all bit-identical to `master`. D-020 accepted. |
| `beam_deflection_formula` | **Yes, in the US-customary uniform-load branch** — Step 2 no longer performs the conversion. F-1, blocking. SI branch: no change. |
| `cantilever_double_integration` | **No.** Same 6 steps, same development, same mm answers. |
| `annulus_flowrate` | **More physically sensible, differently distributed.** 17% of `master`'s instances asserted a laminar solution at Re ≥ 2100 and 17% had their stated ΔP silently overridden; both gone. Median −28% is the low tail being restored, not the item changing. F-4. |
| `statically_indeterminate_shaft` | **No.** Same 3 steps, same method. Tᴀ moves ~0.04% because the trace now consumes the printed ratio (P3 fix). |
| `shaft_design_power` | **No.** Same 4 steps. The 0.002 mm diameter grid *restores* answer space rather than reducing it. |
| `composite_shafts_series` | **No — better.** Same 4 steps, same physics. `master` quoted one quantity at two precisions; the branch makes them consistent. |
