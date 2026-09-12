# Phase 6 — consolidated item-pool impact (D6.4)

**Which published numbers does the redesign invalidate, and what does correcting them cost?**

**Date:** 2026-09-13 · **Branch:** `redesign/phase6-consolidation`
**Consolidates:** [`phase1_`](phase1_item_pool_impact.md), [`phase2_`](phase2_item_pool_impact.md),
[`phase3_`](phase3_item_pool_impact.md), [`phase5_`](phase5_item_pool_impact.md),
[`phaseC3_`](phaseC3_item_pool_impact.md),
[`phaseC3_corrections_`](phaseC3_corrections_item_pool_impact.md), and the two tranches
measured in this phase.

---

## 1. The published pool, and why the question changed shape

| | |
|---|---|
| Published item pool | **1,350 items** = 90 templates × 15 seeded instances, per model |
| Models evaluated | 11 in `inference_results`; the paper reports 27 LLMs |
| Templates in the corpus **today** | **150** — civil and industrial were authored after the published run |

**The 90 is not a subset choice; it is all there was.** Only 20 template modules ever carried a
generator, covering chemical, electrical and mechanical — exactly 90 template functions. Civil and
industrial have never appeared in a published result.

## 2. Every measured event, in order

Each was measured from **two git worktrees in two separate processes**. An in-process reload
resolves both sides to the same already-imported modules and reports everything identical — the
trap `instance_dump.py` documents and which has caught four people.

| phase | templates moved | measured |
|---|---:|---|
| **1** round-trip integrity | 10 (8 in the published pool) | 120 of 1,350 items = **8.9%**; ~74 re-score, ~44 re-inference |
| **2** determinism | 4 | `levenspiel` 300/300; `adiabatic_flame_temperature` 300/300 answers (median 0.845%); `vibration_isolator_design` 300/300 (0.012%); `pfr_volume_changing_rate` 24/300, all a *different draw* |
| **3** trace shape | 2 | `normal_depth_iteration` **2.30%**, `line_balancing_heuristic` **0.65%**; answer *space* unchanged on both |
| **5** contract hygiene | 4 of 11 | emitted text only — **no question content and no answer value moved**, 22,000 instances |
| **C3** constants | 1 + 1 | 1 answer-body (re-score *and* re-inference), 1 question-only (re-inference only) |
| **C3 corrections** | 13 | q 824, ans 713, sol 765; 0 generation errors |
| **Track B** deletions | 3 | 576 instances; q = ans = sol = 576; 0 errors |
| **D6.7** doubled sign | 14 | q 1,008, ans 381, sol 2,234; 0 errors |

**No template gained a generation error in any tranche.**

## 3. The distinction that actually costs money

Three outcomes, and they are not interchangeable:

- **Solution-only** — question and answer byte-identical, only the gold *trace* differs. No
  re-inference, no re-scoring. Ten of D6.7's fourteen are this.
- **Question changed** — the item itself differs, so the model must be re-run. **Re-inference.**
- **Answer changed, question byte-identical** — a **hidden constant**. Nothing about the item
  looks different, and a reader diffing questions sees nothing. **Re-scoring**, and it is the
  direction that gets missed.

### Hidden constants found

| tranche | hidden-constant instances |
|---|---|
| Track B | **0** — every change was a visibly different draw |
| **D6.7** | **19**, all in `template_lorentz_force` |

D6.7's aggregate diff reports `lorentz_force` at q 279 / ans 223. Those totals cannot tell you
whether the two sets coincide, and they do not: checked per instance, **204** changed question
*and* answer, **76** solution-only, and **19** changed the answer while the question stayed
byte-identical. The mechanism is plain once seen — the question prints the **input** vectors, the
answer prints the **computed** cross product, so all-positive inputs can still yield a negative
output component.

**This was found only by a per-instance check.** The aggregate would have reported it as ordinary
re-inference.

## 4. D-003 — closed

Open since 2026-09-05 and carried by prompts 05, 06 and 07: *do the raw generations still exist,
so the parser fix is a free re-score?*

Measured: `inference_results/` **absent**, `evaluation_results/` **absent**, and `testset/`
**absent** — the input corpus the published numbers were generated from is not in this working
copy either. The only `*results*.json` anywhere belongs to the annotation pilot.

**So there is no free re-score.** But the conclusion is better than the spec's "budget decision":
the repo owner has directed that the testset be rebuilt and inference and evaluation re-run, so
there is no stale artefact to correct in place. **Recovery is not attempted, by decision.**

## 5. What supersedes all of it

D6.12 rebuilds the pool from one seeded generator: **150 templates × 15 instances = 2,250 items**,
against the published 1,350. Every tranche above is therefore *historical* — it explains why the
published numbers changed, and no longer describes work anyone has to do.

**The per-tranche figures are exact; the union is not computable from the summaries.** Templates
appear in more than one tranche and several tranches record counts without naming their templates.
Named lists exist for Phase 3 (2), Track B (3) and D6.7 (14). Anyone needing the union must
recompute it from the dumps, not add these rows.

## 6. What a paper revision must not say

- Not *"the benchmark was corrected"* — 60 templates in the current corpus have **never** been
  evaluated.
- Not *"results were re-scored"* — the generations are gone; everything is re-inferred.
- Not a comparison of old and new aggregate scores as though they measure the same instrument.
  The pool changes size (1,350 → 2,250), composition (3 branches → 5) and content.
