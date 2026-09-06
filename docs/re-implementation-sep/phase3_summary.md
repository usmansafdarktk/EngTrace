# Phase 3 — Summary and close-out

**Trace shape: iteration and search** · two templates · branch
`redesign/phase3-trace-shape` off `master` at `00bd0b0`
**Date:** 2026-09-06
**Companions:** [`phase3_node_types.md`](phase3_node_types.md) (D3.3, the primary
deliverable) · [`phase3_item_pool_impact.md`](phase3_item_pool_impact.md) (D3.5)
· [`DECISIONS.md`](DECISIONS.md) D-038 – D-044 ·
[`reviews/`](reviews/) phase3_reviewer_b_pedagogy, phase3_reviewer_d_schema

---

## 1. What this phase was for

Every other template in the corpus emits a **flat milestone list** — a fixed
sequence of steps with a stable `{id, symbol}` per step, known before the
template runs. These two do not. `normal_depth_iteration` runs a secant solve
whose update count varies with the data; `line_balancing_heuristic` runs a greedy
heuristic whose *number of workstations is part of the answer*. Neither fits a
fixed-shape milestone model, and closing that gap is what the milestone model —
and therefore the deterministic verification the whole effort exists to reach —
is built on.

So this is the phase where the **schema**, not the template, was allowed to be
what changes, and where the primary deliverable is a specification rather than
code.

**The framing that turned out to matter.** Both templates come from the corpus's
later, spec-driven authoring era: provenance-tagged constants, textbook
grounding, stated physical bounds, ten `assert`s each, bounded resample loops.
The defects here are not carelessness, and none of them were found by looking for
sloppiness. Every one is a structural consequence of a trace whose shape depends
on its data.

---

## 2. D3.1 — the two route decisions, argued separately

The spec recommends the schema route for both. Both decisions landed there, and
**the two arguments are different**; neither would carry the other's template.

### 2.1 `normal_depth_iteration` — **schema route**

The alternative on offer was the cheap one Phase 2 proved out on
`adiabatic_flame_temperature`: fix the iteration count, print every pass, ~4 h,
makes it class A. Three independent reasons it does not transfer:

1. **The termination predicate is part of the question.** The question says
   *"update the depth by linear interpolation between successive trials until the
   depth changes by less than 0.002 m."* Fixing the count means rewriting that to
   "do exactly three updates" — changing what is asked, not how it is recorded.
2. **Three updates is not enough for a third of the pool.** Measured over 4,000
   seeds the update count is `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}`:
   **32.9% of instances need more than three.** Truncating at three changes the
   gold answer on those, and screening the sample space so three always suffices
   biases the pool toward channels whose `g` is near-linear — a quiet difficulty
   reduction.
3. **Phase 2 Reviewer C's limit.** The flame template is verifiable at a fixed
   count only because its parameter space is a closed set of eleven discrete
   reactions, so the rounding-boundary question could be checked *exhaustively*.
   Normal depth samples geometry, slope and roughness continuously; a fixed count
   cannot be exhaustively validated, and would need either a proven a-priori error
   bound over the sampled box or a run-time convergence check not stripped by
   `python -O`.

Phase 2 §10 offered the criterion *"fix the iteration count when the count is not
part of what the item tests"*, and §12 corrected it to *necessary but not
sufficient — it also requires a parameter space you can verify over*. **Both
halves of the corrected criterion fail here**, and reason (2) is a third the
criterion does not cover at all.

**P6 trade: none taken.** The item continues to test what it tested — running an
iterative scheme to a stated tolerance. The cost is paid by the schema, which
must carry a variable-length node.

### 2.2 `line_balancing_heuristic` — **schema route, for an unrelated reason**

The cheap route is not merely worse here, it is **not the right shape of
argument**. A heuristic has no fixed point, so "N passes suffices" cannot be
established at all — Phase 2 Reviewer C's point, and it is decisive.

The binding reason is different again and is the whole of P6 for this template:
**the station count `n` is a computed result that appears in the answer.**
Balance delay is `(n·CT − Σt)/(n·CT)·100`. Fixing `n` to make the trace a fixed
shape publishes part of the answer in the question.

The specific thing the brief asked to check — **is the heuristic deterministic
under tie-breaking?** — resolves cleanly and by construction:

- The five task durations are `random.sample(range(15, 111), 5)`, so they are
  **distinct integers**. Measured: **0 duplicate-duration instances in 4,000
  seeds**, so `max(fitting, key=duration)` is always unique.
- Eligibility is built by scanning the **tuple** `_T19_ORDER`, not a set, so no
  implementation-defined iteration order reaches the output. T3 (determinism)
  passes, before and after.

So there was no non-determinism to fix. The tie-break rule is nonetheless now
**stated explicitly** in the node's `rule` field, because a verifier must be
total over traces a *model* might produce, and a model may well emit a tie.

**P6 trade: none taken.**

**Why these were argued separately.** They reach the same route, but if the
`normal_depth` argument were applied to `line_balancing` it would conclude that
the count is incidental — and §3 below is the finding that this is exactly
backwards.

---

## 3. The finding: the two node types are not interchangeable (D-038)

`template_redesign_spec.md` §3.2 says the two node types "are the same shape (an
ordered list of homogeneous sub-traces with stable within-element symbols and a
termination predicate)" and recommends specifying them together as the cheaper
route. **The structural claim is right; the equivalence is wrong**, and building
one node type from it yields a verifier that is wrong on one of the two
templates.

They differ on the property that decides how a candidate trace is *scored*:

> **Is the cardinality of the sequence an observable of the answer?**

| | `normal_depth_iteration` | `line_balancing_heuristic` |
|---|---|---|
| Count enters the answer? | No — the answer is the converged depth | **Yes — `n` is in the delay formula** |
| Count stable under permitted slack? | **No** — the test is `\|Δy\| < 0.002` on 4-dp values, so a solver carrying more digits converges in a different number of updates | Yes — every fit decision is an integer comparison |
| Comparator disposition | count mismatch **tolerated** | count mismatch is an **answer failure** |

D3.4's comparator must therefore tolerate on one and fail on the other. One node
type cannot express both dispositions. The discriminating rule, written so it
applies to a template this phase never saw:

> Cardinality is `incidental` when it is sensitive to the numerical slack the
> comparator already tolerates, and `answer_bearing` when it is invariant under
> that slack and appears in the answer.

Filed as a `SPEC-CHANGE` against §3.2. What survives of §3.2 is "specify them
together" — they share one base structure, one `carry` mechanism, one
verification algorithm. What is struck is the implication that they share a
comparator.

---

## 4. D3.3 — the node specification, and how it is made testable

[`phase3_node_types.md`](phase3_node_types.md) specifies both types: a common
`sequence` base (§3), the two types (§4, §5), a **normative verification
algorithm** (§6), the **D3.4 comparator semantics** (§7), and — the section that
makes it a specification rather than a description — **§8, the traces a
conforming verifier must reject**.

The mechanism that makes a variable-length sequence checkable is `carry`: a map
from each symbol of element *k+1* to a path expression over element *k*. It lets
a verifier check element *k* against *k−1* in isolation, without knowing the
count in advance. That is what replaces the flat model's per-step `{id, symbol}`.

**One structure, two presentations (D-039).** Both templates now *build* their
trace as a structured node and **render the printed prose from it**. This is
D-034's rule applied to traces: if a template built a node and formatted prose
independently, the two could disagree and nothing would notice — the same shape
as a test carrying its own answer key. `tests/trace_schema/extract.py` therefore
holds no reference values of its own; it lifts the node the prose came from.

The node is a frame local, not a return value, deliberately: changing the return
contract is a corpus-wide decision, not a two-template one. The cost is stated —
nothing outside the extractor can consume it yet, and no check gates its shape.
Promoting it belongs to the milestone-model work.

**Evidence the restructure was behaviour-preserving:** across 4,000 seeds, on
every seed the screens did not resample, **every arithmetic line of both
solutions is byte-identical** to `master`. The only changed lines anywhere are
two deliberately added prose sentences in the civil item.

---

## 5. D3.2 — what changed in the templates

### 5.1 `normal_depth_iteration`

| Change | Why | Measured |
|---|---|---|
| `fmt3()` removed; the **stored** `g` is normalised instead of the printed string | The helper rewrote `"-0.000"` to `"0.000"`, so the printed operand could differ in sign from the stored one — a direct P2 violation | **0 occurrences in 4,000 seeds.** Latent, not active (D-042) |
| `g_curr - g_prev` guarded by an explicit `raise` | Unguarded division | **0 in 20,000 seeds**; smallest `\|g_k − g_{k−1}\|` is 0.004, four steps of the 3-dp residual grid |
| Non-convergence raises instead of relying on an `assert` | Step 4 told the reader the tolerance was met while only an assert enforced it, and `python -O` strips asserts | 0 non-convergences in 20,000 seeds; guard **verified to fire** under `-O` (§7.3) |
| Step 4 now **prints** the last change rather than asserting it is small | The reader can check the claim instead of taking it — the same move as Phase 2's C-2 fix | — |
| Display-tie resample screen, **two** populations | D-016: at a tie no rounding convention closes in both directions | 2.30% of seeds redraw |

### 5.2 `line_balancing_heuristic`

| Change | Why | Measured |
|---|---|---|
| Display-tie screen on the exact rational balance delay | D-016 | 0.65% of seeds redraw |
| `n*CT - total` and the duration sum bound out of the f-string | T5 result-inline; Phase 2 Reviewer C's C-3 pattern | 2 T5 findings → 0 |
| `_t19_assign` returns the `decision` node; prose rendered from it | D-039 | 0 output change |
| Tie-break rule stated explicitly | A verifier must be total over model traces | never exercised |

### 5.3 The defect the brief did not list (D-040)

The brief lists one T1 defect for the civil template: *"1 hard failure in 200
seeds (seed 11, the `K = Q*n/S^(1/2)` line)"*. Accurate, and not the whole
defect. Measuring the **exact rational** value of every printed expression against
its display precision over 20,000 seeds, rather than counting T1 hard failures,
finds two tie populations of comparable size:

| Line | Display | Exact ties | T1 hard failures / 20,000 |
|---|---|---|---|
| Step 1, `K = Q*n/S^(1/2)` | 3 dp | **1.14%** | 22 |
| Step 3, the secant update | 4 dp | **1.02%** | 5 |

**Counting T1 hard failures understates a tie population by roughly 10×** — 27
hard failures against 427 ties in the same seeds. D-016 explains why: at an exact
tie `|evaluated − printed| == tol` and float representation error alone decides
FAIL versus MARGINAL. A tie reading MARGINAL is the same ill-posed instance as
one reading FAIL. **Any phase that sizes a display-tie fix by counting T1
failures will under-fix by an order of magnitude.** That generalises beyond these
two templates and is the most transferable thing this phase found.

D-037's cheap escape was tested before resampling and does not apply — `K`
terminates at its 3-dp display on only 10.9% of instances and at 6 dp on 12.9%,
so rounding persists at any length and the ties would relocate rather than
vanish. Establishing that properly produced **D-044**: D-037's test must be read
as *exact for every instance*, never *exact for the tie instances*, since the
latter is circular and always answers yes.

---

## 6. Verification

All commands run with `PYTHONIOENCODING=utf-8`.

### 6.1 The two templates

| Check | Before (`00bd0b0`) | After | Note |
|---|---|---|---|
| T1 closure, 200 seeds | 2 failing | **0** | |
| T1 closure, **20,000 seeds** | 27 + 204 hard failures | **0 + 0** | the acceptance run |
| T2 round-trip | pass (no oracle) | pass | neither has an oracle |
| T3 determinism | pass | pass | |
| T4 output contract | pass | pass | |
| T5 binding | 2 failing (3 + 2 findings) | **0** | |
| T6 distribution | fail | fail | stale baseline, corpus-wide — §6.3 |
| T7 asserts | pass | pass | assert count unchanged; guards **added** |

### 6.2 Corpus non-regression, measured not assumed

`master` was re-measured in a worktree at `00bd0b0` rather than trusted from the
brief; the brief's table reproduced exactly (T1 30, T2 0, T3 0, T4 3, T5 68,
T6 142, T7 83). Comparing full `--checks all` runs per template across all 150:

```
templates in before/after: 150 / 150
check-status changes: 3
  template_line_balancing_heuristic   T5: False -> True
  template_normal_depth_iteration     T1: False -> True
  template_normal_depth_iteration     T5: False -> True
```

**Three changes, all improvements, all in scope. 147 templates unchanged on
every check.** (`line_balancing`'s T1 shows no change here because at the
runner's default 25 seeds it passed both before and after; its 204 failures per
20,000 seeds appear only at scale.)

Constants suites: `test_chemical_thermochemistry` 222 checks all pass,
`test_citations_resolve` 55 checks all pass.

### 6.3 The T6 wall, and why the baseline was not regenerated (D-043)

T6 fails 142/150 on `master`; the committed baseline is stale corpus-wide, so it
cannot gate anything here. Regenerating it would turn T6 green and is the wrong
move: **a baseline refreshed by the phase it is meant to gate is not a gate**,
and it would silently absorb corpus-wide movement unrelated to Phase 3. Phase 2
hit this and delivered the diff directly; Phase 3 does the same, via a
before/after instance dump of 4,000 seeds per template **run as a separate
process per tree** (`tests/template_integrity/phase3_instance_dump.py`, which
refuses to run if a template resolves outside the tree root it was given — the
in-process-reload trap has now caught three people).

Regenerating the baseline corpus-wide remains right **as its own deliberate
change on `master`**, and is left to Phase 6.

### 6.4 Resolution limits, stated

Per the standing rule (D-024, D-026): N seeds cannot resolve a rate below ~3/N.

| Run | Seeds | Smallest resolvable rate | Rate being excluded |
|---|---|---|---|
| T1 acceptance | 20,000 | 0.015% | was 0.135% / 0.510% — resolved |
| tie census | 20,000 | 0.015% | 2.14% / 0.65% — resolved |
| `fmt3` firing | 4,000 | 0.075% | 0 observed — **excluded only above 0.075%** |
| zero denominator | 20,000 | 0.015% | 0 observed — **excluded only above 0.015%** |
| non-convergence | 20,000 | 0.015% | 0 observed — **excluded only above 0.015%** |
| item-pool churn | 4,000 | 0.075% | 2.30% / 0.65% — resolved |

The three zero-observation rows are the honest weak points: those defects are
excluded above their stated rate and above nothing below it. Two of them are now
guarded by an explicit `raise`, which is what covers the remainder.

---

## 7. Item-pool impact (D3.5)

Full numbers in [`phase3_item_pool_impact.md`](phase3_item_pool_impact.md).
Headline:

| | `normal_depth` | `line_balancing` |
|---|---|---|
| Instances resampled | 2.30% | 0.65% |
| Answers changed | 2.30% | 0.62% |
| …of which kept the same question | **0** | **0** |
| Distinct answers | 267 → 267 | 335 → 335 |
| p5 / p50 | identical | identical |
| p95 | −0.518% (one display bin) | identical |
| Branch mix | 0.10 pts | 0.22 pts |

No instance kept its question and changed its answer — that is the number that
matters, because it means no solver who saw a `master` instance would now be
graded against a different answer for the same problem.

**The industrial screen is not uniform, and that is recorded rather than
smoothed.** Of its 26 rejections, 22 have `n = 4` and 4 have `n = 3` — a 6×
skew, because `n·CT` at `n = 4` carries two extra powers of two and so lands on a
terminating tie more readily. Net movement on the station split is 0.22 points
against a 5-point tolerance, so it is immaterial here; it is written down because
"the screen rejects at random" was the comfortable claim and is not the true one,
and a future screen on a template with a finer-grained answer-bearing branch
could be moved by the same mechanism.

---

## 8. Errors I made in this phase

Every phase so far has recorded at least one. This one has five, and the
fifth is a process error rather than a numerical one.

1. **A 2× transcription slip in a docstring I wrote.** The civil template's new
   docstring cited the update-count distribution `{1: 16, 2: 163, 3: 1182, ...}`
   as "measured over 4,000 seeds". It is the **2,000**-seed figure. Caught by the
   4,000-seed item-pool dump disagreeing with it. The correct figure is
   `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}` — and correcting it *strengthened*
   the D3.1 argument, since the fraction of instances needing more than three
   updates is **32.9%, not the 16% I had first written**. An argument built on a
   number I had halved.
2. **My own worked example violated my own specification.** The `decision`
   example in `phase3_node_types.md` §9.2 omitted `carry`, a field the same
   document's §3 marks required. Caught by machine-comparing every JSON block in
   the spec against the committed conformance corpus — not by re-reading it,
   which I had already done twice.
3. **An arithmetic reading off in the fourth decimal.** The §9.1 walkthrough said
   the secant step gives `1.37855…`; it gives `1.378646…`. Same cause, same
   catch. Both corrected before the reviewers were dispatched.

4. **A justification that was true in spirit and false as written.** Both
   templates' resample-screen comments said the quantity "does not terminate at
   any fixed number of places". A denominator of the form 2^a·5^b does terminate,
   so the claim is false as stated and would not have survived a reviewer with a
   calculator. Caught by going back to check it before the reviewers did.
   Measuring it properly (D-044) turned up something better than a fix: the
   *naive* reading of D-037 is itself a trap, because a value on a half-way
   boundary at *k* places is by definition exact at *k+1* places, so "can I
   lengthen the display?" always answers yes and always relocates the tie
   population rather than removing it. The decision was right; only the reason
   was wrong, and the corrected reason generalises.

5. **I moved the branch under a reviewer, having been warned not to.** The brief
   for this phase says explicitly: *"Give each reviewer a frozen ref — a commit
   SHA, not a branch name"*, because in Phase C2 a branch moved under a reviewer
   and roughly a quarter of that review was wasted on stale findings. I did give
   both reviewers a SHA — and then committed D-044 to the same branch while they
   were still working. Reviewer B caught it as its **first** action, diffed
   `953adc9` against the working tree before trusting a single number, and
   established that both drifts were comment-only so its measurements held.
   Nothing was lost, and only because the reviewer was more careful than the
   implementer. The rule I had internalised was "hand out a SHA"; the rule that
   matters is **"do not commit to a branch a review is in flight against"** — a
   frozen ref is a promise about the repository, not a string in a brief. Round 2
   was dispatched against `9398ccc` with no commits made until it filed.

The pattern in 2 and 3: **a document that quotes a computed artefact must be
checked against that artefact by machine.** Re-reading finds neither. That check
is four lines and should be standard for any deliverable that embeds generated
data. The pattern in 4 is different and worth separating: **a justification that
is qualitatively right can be quantitatively false**, and the way to find out is
to measure the thing you asserted rather than to re-read the sentence.

---

## 9. Residual risk register

| # | Risk | Severity | Disposition |
|---|---|---|---|
| R3-1 | The `iteration`/`decision` generalisation to `linear_reservoir_routing_step` and `qr_policy_one_iteration` is **argued, not measured**. Neither has been fitted. | medium | BACKLOG → apply the `incidental`/`answer_bearing` test to both before the milestone model freezes |
| R3-2 | The structured node is a frame local with no interface and no check gating its shape. It can silently drift out of the spec. | medium | D-039; ADOPT in the milestone-model phase — promote to a return value and add a conformance check |
| R3-3 | Three defects are excluded only above 0.015–0.075%. `fmt3` in particular was never observed firing. | low | Two are now guarded by an explicit `raise`; the third is removed by construction |
| R3-4 | T6 cannot gate either template until the corpus baseline is regenerated. | medium | D-043; Phase 6 owns it |
| R3-5 | The corpus-wide sweep Phase 2 Reviewer A recommended — templates combining a transcendental with a fixed-decimal print — was **not** run. `normal_depth` is exactly that shape and was fixed; the other candidates were not looked for. | medium | BACKLOG, carried from Phase 2 §11 |
| R3-6 | Phase 2 Reviewer C's corpus-wide regex for `\{[a-z_]+ *[-+*/] *[0-9.]+\}` in f-strings was applied to these two templates only. | low–medium | BACKLOG, carried from Phase 2 §12 |
| R3-7 | Two rounding conventions still coexist corpus-wide (`_hu` decimal half-up vs `_as_printed` binary). | medium | BACKLOG, carried from Phase 2 §11 |
| R3-8 | Prose/node agreement is specified as a non-goal, not a check — vacuous by construction today because the prose is rendered from the node, and a real gap the day a template stops doing that. | low | Recorded in the spec's §10 |
| R3-9 | A solver who solves the civil item **directly** rather than by the prescribed secant gets a different gold answer on ~3.7% of instances (always by 0.001 m). Pre-existing on `master`, surfaced by Phase 3's Step 4 change. D3.4 §7.4 now says the comparator must accept within the method's own tolerance, but no extractor enforces it yet. | medium | Reviewer B F-B3; ADOPT-NOW in the spec, enforcement deferred to the extractor spike |
| R3-10 | Both display-tie screens reject **clustered**, not scattered, instances — the civil one at a single slope value. Immaterial here (0.10 and 0.22 pts of branch movement) but the mechanism is general. | low–medium | Reviewer B F-B1/F-B2; recorded as **D-045** |

---

## 10. Exit gate

| Gate item | Status |
|---|---|
| D3.1 decisions recorded, **argued separately per template**, P6 trade stated | ✅ §2 — schema route both, on unrelated arguments, no P6 trade taken on either |
| T1–T7 pass on both templates | ⚠️ T1–T5 and T7 pass; **T6 blocked corpus-wide by a stale baseline** (D-043), replaced by a direct distribution diff |
| `g_curr - g_prev` guarded; the `fmt3()` sign discrepancy removed | ✅ §5.1 |
| **Reviewer D implemented a working verifier from D3.3 alone** | *see §11* |
| D3.4 comparator semantics cover the iteration-count-mismatch case | ✅ `phase3_node_types.md` §7.2 |
| No corpus regression against the master baseline, measured not assumed | ✅ §6.2 — 3 status changes across 150 templates, all improvements, all in scope |
| Reviewers B and D filed, every finding and every §5 suggestion triaged | *see §11, §12* |
| Item-pool impact note filed | ✅ [`phase3_item_pool_impact.md`](phase3_item_pool_impact.md) |

---

## 11. R4 triage — Reviewer D (schema implementability)

*Pending — the review is in flight. This section is completed before the phase
closes; an untriaged suggestion blocks the gate exactly as a CONFIRMED finding
does (R4).*

---

## 12. R4 triage — Reviewer B (pedagogy)

[`reviews/phase3_reviewer_b_pedagogy.md`](reviews/phase3_reviewer_b_pedagogy.md).
**PASS WITH FINDINGS.** Five findings, all closed.

### 12.1 The gate question, answered

**Is `line_balancing_heuristic` a lookup?** No, and B settled it with a number I
had not measured and deliberately did not commission from myself:

| Predictor of `n` from the question text alone | Accuracy, 8,000 seeds |
|---|---|
| `n = N_min = ceil(Σt/CT)` — the trivially computable bound | 58.33% |
| `n = 3` — best blind constant | 53.05% |
| `n = max(N_min, 3)` — **best rule B could find** | **64.97%** |
| per-`N_min` majority vote — the ceiling of any `N_min`-based rule | 65.55% |

So the best shortcut buys **~12 points over guessing**, and the greedy rule must
actually be executed on at least 35% of instances. The structure behind it is
that `N_min` determines `n` only in its tails (`N_min = 2 → n = 3`, 277/277;
`N_min = 4 → n = 4`, 497/497), while **81% of the pool sits at `N_min = 3`, where
the split is 1848/1378 — a 57/43 coin.** B also searched exhaustively over
`n = N_min + [cheap predicate]` for all 15 constructible predicates; nothing beat
65.0%.

Master measures 58.27% / 64.84% on the same seed span, against the branch's
58.33% / 64.97%. **Phase 3 did not make this item more guessable.**

This is the single most valuable thing either review produced, because it
converts the D3.1 argument for `line_balancing` from a structural claim into a
measured one.

### 12.2 Findings

| # | Finding | Disposition | Action |
|---|---|---|---|
| **F-B1** | **CONFIRMED.** The civil `K` display-tie screen is **single-valued in slope**: all 718 rejections in 60,000 draws sit at `S = 0.0016` exactly, because it is the only 4-dp grid value in the window whose square root is exact at 5 dp (`0.04`), making `K = 25·Q·n` terminating and a 3-dp tie reachable at all. It depletes that slope 4.83% → 3.50%, ~28% relative. | **ACCEPTED — and it corrects a claim of mine** | Reproduced exactly: 709 K-rejections in 60,000, **1 distinct slope**, 25.4% of that slope removed, against 23 distinct slopes for the update screen. My item-pool note §2.4 said the rejections were "scattered", measured on the **combined** rejection set where the update screen's spread masks the K screen's concentration. Corrected there, recorded in the docstring, and generalised as **D-045** |
| **F-B2** | **CONFIRMED.** The industrial screen is also clustered: 913/913 rejections have `n·CT` divisible by 8, and 85.5% have `n = 4` against 46.6% in the pool. Net effect 0.23 pp. | **ACCEPTED (convergent)** | I had found this independently and recorded it in `phase3_item_pool_impact.md` §3.4 before the review landed — 22/26 = 84.6% against B's 85.5%. **Two different questions reaching the same number is the kind of convergence R6.2 exists to produce**; B derived it from `n·CT` divisibility, I from the tie census |
| **F-B3** | **PLAUSIBLE.** The new Step 4 sentence makes a pre-existing precision mismatch legible: the item converges to ±2 mm but reports `yn` to 3 dp, so a solver who solves Manning's equation **directly** rather than by the prescribed secant gets a different gold answer on ~3% of instances, always by exactly 0.001 m. | **ACCEPTED as real; NOT a Phase 3 regression; comparator fix ADOPT-NOW, question-text fix BACKLOG** | Verified on **both trees**, which is the claim that decides the disposition: master **3.70%**, branch **3.65%** — statistically identical, so Phase 3 surfaced it and did not cause it. (B reports 3.08% from a differently-filtered probe; the divergence is in the accepted-instance filter, not the phenomenon.) The real fix belongs in D3.4, not the template: see §12.3 |
| **F-B4** | **PLAUSIBLE.** Step 3 now discloses the update count before the trace, removing a little "have I finished?" work. | **ACCEPTED as a deliberate trade, recorded** | It is not an answer leak — the count is incidental (§3) — and the sentence exists to make exactly the point that the stopping rule is a tolerance and not a recipe. But B is right that it is the one change pointing toward "easier", and it was not recorded as a trade before B named it. Recorded now rather than argued away |
| **F-B5** | No finding on engineering content: both items remain faithful to Sturm Ch. 4 and Nahmias §9.10, both correctly Advanced. | **RECORDED (favourable)** | This is the P6 gate itself, and it is the section that makes the verdict mean something |

### 12.3 F-B3's disposition, in full

B offers two fixes: state in the question that the answer is the depth *the
scheme returns*, or tighten the convergence tolerance to 0.0005 m so the secant
result and the true root agree at 3 dp.

**Both are rejected for this phase, and a third is adopted.** The first rewrites
the question text on 100% of instances to resolve a 3% grading ambiguity, in a
phase scoped to trace shape; the second changes the item's stated tolerance,
which is a P6 change to what the item asks for. Neither is a trace-shape fix.

The defect is really in the **comparator**, which is this phase's own deliverable:
when a question *prescribes a method*, the gold answer is that method's output,
and the answer tolerance must be the looser of the display tolerance and the
method's own convergence tolerance. Marking a direct solver wrong for 0.001 m is
the same error as marking a four-update solver wrong for using four updates —
§7.2 already refuses the second and should refuse the first. **Added to D3.4**
(`phase3_node_types.md` §7.4) and carried on the risk register as **R3-9**, since
the extractor that will enforce it does not exist yet.

### 12.4 What B bought, and what it cost me

B ran the one measurement I could not run for myself without answering my own
gate question, and it ran it six ways plus an exhaustive predicate search before
concluding. It also **checked the ref for drift before trusting a single number**
— and found drift, because I committed D-044 while the review was in flight
(§8.5). Both drifts were comment-only and B verified that before proceeding, so
nothing was wasted; but a reviewer who had not checked would have been silently
reviewing something else, which is the exact failure the frozen-ref rule exists
to prevent and which I reintroduced by hand.

---

## 13. What Phase 4 inherits

Phase 4 is the four non-numeric verifiable templates, and the spec is explicit
that they **must not be redesigned**: the schema adapts, not the items. Phase 3
is the precedent for that, and hands over three things.

1. **The node model, and the axis it turns on.** Phase 4's templates have zero
   numeric content and are verified by a non-arithmetic comparator. The question
   to ask of each is the D-038 question in its general form: *what property of
   the trace is answer-bearing, and what is incidental?* For an iteration it is
   the count; for a classification item it will be something else. The
   `incidental` / `answer_bearing` split is the reusable part, not the two
   concrete node types.
2. **A warning about the comparator.** D3.4 §7 specifies two dispositions for two
   node types. Phase 4 adds a third kind of item, and the temptation will be to
   reuse whichever existing disposition is closest. This phase's central finding
   is that doing exactly that — merging two node types because they look alike —
   produces a verifier that is wrong on one of them.
3. **A measurement discipline that transfers.** Count the *ill-posed instances*,
   not the check failures. T1's FAIL/MARGINAL split is decided by float
   representation, not by the defect, so it understates a display-tie population
   by about 10× (§5.3, D-040). Phase 4 has no arithmetic and so no display ties —
   but the general form, *the check you have may be measuring a proxy for the
   defect rather than the defect*, is the reason Phase 0 found two checks that
   were green while measuring nothing.

**Still open and not Phase 4's to fix**, but carried on the register: the
corpus-wide sweeps from Phase 2 (R3-5, R3-6, R3-7), and the T6 baseline
regeneration (R3-4), which Phase 6 owns.
