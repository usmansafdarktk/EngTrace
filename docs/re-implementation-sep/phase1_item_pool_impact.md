# Phase 1 — Item-pool impact note (D1.6)

**Which published results does regenerating these 12 templates invalidate?**

**Date:** 2026-09-06 · **Branch:** `redesign/phase1-round-trip`
**Companion to:** [`phase1_summary.md`](phase1_summary.md), [`DECISIONS.md`](DECISIONS.md) D-008

---

## 1. Headline

| | |
|---|---|
| Published item pool | **1,350 items** = 90 templates × 15 seeded instances, per model |
| Templates this phase edits | **12** |
| …of which are **in** the published pool | **8** |
| Published items affected | **120 of 1,350 = 8.9%** |
| …needing only a **re-score** (question unchanged, gold trace corrected) | **~74 of 120** |
| …needing **re-inference** (the item itself changed) | **~44 of 120** |
| Model runs affected | all of them — the paper evaluates 27 LLMs |

**Four of the 12 templates carry no published-results cost at all.** The paper's
90 templates span **three** branches — chemical, electrical and mechanical
(§3.1, "Template Selection"; the domain-ranking prompt at line 1281 names them).
Civil and industrial were authored later and are not in the 1,350. So
`beam_deflection_formula` and `cantilever_double_integration` — both civil — are
free to change. `statically_indeterminate_shaft` and `shaft_design_power` are
mechanical and *are* in the pool.

---

## 2. Per template

Measured by generating each template before (at `master`) and after, in
**separate processes** so the two trees cannot share imported modules, and
comparing the `(question, solution)` pair per seed.

Three outcomes are distinguished, because they cost very different amounts:

- **identical** — nothing to redo.
- **gold-only change** — the question is byte-identical and only the solution
  changed. **The archived model generations remain valid**; what changed is the
  gold trace they are scored against. This is a re-score, not re-inference —
  the same split D-008 draws for the parser fix.
- **item changed** — the question text changed, so the item a model was asked
  is no longer the item in the pool. This needs re-inference.

### Over the first 15 seeds (the pool's own shape: 15 instances per template)

| Template | Branch | In the 1,350? | identical | gold-only | item changed |
|---|---|---|---:|---:|---:|
| `mean_variance` | electrical | yes | 0/15 | 7/15 | **8/15** |
| `annulus_flowrate` | chemical | yes | 0/15 | 0/15 | **15/15** |
| `rotating_unbalance` | mechanical | yes | 0/15 | 0/15 | **15/15** |
| `vibration_transmissibility` | mechanical | yes | 0/15 | 12/15 | **3/15** |
| `damping_classification` | mechanical | yes | 0/15 | 12/15 | **3/15** |
| `statically_indeterminate_shaft` | mechanical | yes | 2/15 | 13/15 | 0/15 |
| `shaft_design_power` | mechanical | yes | 0/15 | 15/15 | 0/15 |
| `composite_shafts_series` | mechanical | yes | 0/15 | 15/15 | 0/15 |
| `beam_deflection_formula` | civil | **no** | 6/15 | 9/15 | 0/15 |
| `cantilever_double_integration` | civil | **no** | 0/15 | 15/15 | 0/15 |
| **published-pool total** | | | **2/120** | **74/120** | **44/120** |

Of those 44 item changes, **43 state genuinely different values**; 1 differs only
in number formatting.

### Rate over 1,000 seeds

The 15 seeds above are a proxy: **the seeds actually used to build the published
pool cannot be recovered from this working copy** — `inference_results/` and
`evaluation_results/` are gitignored and absent (D-003). The rates below are the
robust version of the same measurement and are what to plan against.

| Template | identical | gold-only | item changed |
|---|---:|---:|---:|
| `mean_variance` | 0.0% | 57.0% | 43.0% |
| `annulus_flowrate` | 0.0% | 0.0% | **100.0%** |
| `rotating_unbalance` | 0.0% | 0.0% | **100.0%** |
| `vibration_transmissibility` | 0.0% | 97.4% | 2.6% |
| `damping_classification` | 0.0% | 74.1% | 25.9% |
| `statically_indeterminate_shaft` | 6.9% | 92.9% | 0.2% |
| `shaft_design_power` | 0.0% | 100.0% | 0.0% |
| `composite_shafts_series` | 0.0% | 100.0% | 0.0% |
| **published-branch aggregate** | **0.9%** | **65.2%** | **34.0%** |

Reproduce with `_dump.py` in a `master` worktree and in this tree, then diff.

---

## 3. Why each template's questions changed, where they did

Only three templates change the **question**, and in each case the question was
part of the defect:

- **`rotating_unbalance` (100%)** — the question now states that the total mass
  *includes* the rotating unbalanced mass. It did not before, and the solution
  assumed it did, so a solver who subtracted it was marked wrong. The stated
  damping coefficient also moved to the value the chain consumes. **The old
  question was ambiguous; every archived generation against it was scored
  against an under-specified item.**
- **`annulus_flowrate` (100%)** — the stated pressure drop moved from a 10 Pa
  grid to a 1 Pa grid. On 284 of 1,000 instances the old grid rounded the
  required drop to zero and the template silently substituted 10 Pa, discarding
  its own Reynolds targeting. **Those items stated a pressure drop unrelated to
  the flow they describe.**
- **`mean_variance` (43%)** — the stated probabilities are renormalised to sum
  to exactly 1.000. They previously summed to 0.999–1.002, so the question did
  not describe a probability distribution.
- **`damping_classification` (26%)** — the stated damping coefficient is now the
  half-up rounding rather than `round()`'s, moving the last displayed digit on
  some instances.

Everywhere else the question is untouched and only the gold trace is corrected.

---

## 4. What this means for the published tables

1. **No published number is safe to leave as-is for these 8 templates.** They
   are 8.9% of the pool, spread across all 27 evaluated models and both the
   final-answer and process-supervision metrics.
2. **About two-thirds of the affected items (65%) need only a re-score.** The
   question is unchanged, so archived generations — *if they still exist* — can
   be re-scored against the corrected gold trace with no new inference. Whether
   they exist is **D-003, still open**.
3. **About one-third (34%) need re-inference** because the item itself changed.
   At 8 templates × 15 instances × 34% ≈ **41 items × 27 models ≈ 1,100
   generations** — small, and cheap enough that it should not gate anything.
4. **The direction of the correction is not neutral.** For
   `rotating_unbalance` and `annulus_flowrate` the *old* items were defective as
   questions, not merely as traces. Any re-scored comparison against the old
   numbers is comparing against a partly ill-posed baseline, and the write-up
   should say so rather than presenting a like-for-like delta.

**Recommendation — unchanged from D-008:** do not regenerate now. Regenerate the
pool **once**, after Phase 6, so comparisons are not invalidated repeatedly.
This note exists so that when that happens, the cost and the caveats are already
counted.

---

## 5. What is not established here

- **Which seeds built the published pool.** Not recoverable from this working
  copy (D-003). The 15-seed table assumes seeds 0–14; the 1,000-seed rates do
  not depend on that assumption and should be preferred.
- **Whether the archived generations still exist.** D-003 is open, and it
  decides whether the 65% re-score slice is free or has to be re-inferred too.
- **Any effect on the AI-Tribunal process scores** beyond the arithmetic. A
  corrected gold trace changes what the tribunal compares against; whether that
  moves the process-quality numbers more or less than the final-answer numbers
  is not measured here and is a Phase 6 question.
