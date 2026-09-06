# Decision and pivot log

Running record of every decision, reversal and scope change on the template
redesign work. **Append-only.** When a decision is superseded, mark the old
entry `SUPERSEDED` and add a new one — never edit history, because the reason a
decision changed is usually more useful later than the decision itself.

Companion to [`template_redesign_spec.md`](template_redesign_spec.md) (the plan),
[`template_audit_report.md`](template_audit_report.md) (the findings),
[`phase0_baseline.md`](phase0_baseline.md) (the measurements) and
[`reviews/`](reviews/) (independent reviews).

**Status values:** `OPEN` (needs sign-off) · `DECIDED` · `SUPERSEDED` · `REVERSED`.

---

## D-001 — Proceed with milestone verification, with a modified design

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** audit report §1

The hypothesis holds: 93% of f-string interpolations are already bound
variables, and 94.5% of step-result numbers are recoverable from the generator's
locals. Instrumentation is *emission*, not rewriting.

Three framing assumptions dropped:

1. **Static milestone declaration is impossible.** 48% of templates change the
   governing equation, step count, or milestone set with the sample. The trace
   must be emitted per instance.
2. **"Milestone = scalar with a unit" is insufficient.** 61 of 150 templates
   carry non-scalar content.
3. **"No LLM in the critical path" cannot be claimed universally.** ~13% of the
   corpus cannot be verified numerically in principle. The defensible claim is
   "deterministic verification for N of 150, remainder explicitly scoped".

---

## D-002 — Do the extractor before instrumenting templates

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** audit report §9, §11.3

The binding constraint is the *model* side, not the template side.
`engineering_parser.py` reduces a solution to one unitless float, reads
`**4,921**` as `4.0`, and cannot see the `**Final Answer**` marker five
templates use. 74.7% of gold answer blocks hold more than one distinct number.

Instrumenting 150 templates while the model side still yields one float
produces a structured gold trace with nothing to compare against. The extractor
spike is the gate for the whole project.

---

## D-003 — Parser fix and re-score are urgent and independent

**Date:** 2026-09-05 · **Status:** OPEN (blocked on a dependency check)

The thousands-separator bug corrupts the gold reference for two templates in
100% of instances; ≥43 archived rows state the correct value and are scored 0.
Two-line fix.

**Blocking dependency:** `inference_results/` and `evaluation_results/` are
gitignored and **absent from this working copy**. If the raw generations were
retained elsewhere this is a free re-score; if not it needs full re-inference
across 11 models × 1,350 items, which is a budget decision.
**Action: establish which, before promising a corrected table.**

---

## D-004 — `damping_classification`: constrain the sample so `c_c` is exact

**Date:** 2026-09-05 · **Status:** SUPERSEDED by D-010 · **Source:** spec §1.3

Original reasoning: after applying P3, ζ recomputed from a 2-dp `c` would never
equal exactly 1.0, so "Critically Damped" becomes unreachable. Verified 680,869
feasible (k, m) pairs exist within current ranges, so constraining the sample
was feasible and preserved all three classes.

---

## D-005 — `normal_depth_iteration`: schema `iteration` node, not a fixed count

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** spec §3.1

Fixing the iteration count converts "iterate until converged" into "perform
three updates" — a materially easier item. The convergence count is genuinely
data-dependent (2–6 evaluations measured).

*Reverses if:* the milestone model is frozen before Phase 3 and cannot take a
new node type. Then fix the count and record the P6 trade.

---

## D-006 — `line_balancing_heuristic`: schema `decision` node (REVERSAL)

**Date:** 2026-09-05 · **Status:** DECIDED — reverses an earlier recommendation
**Source:** spec §3.2

**Originally recommended** redesigning to a fixed-shape decision table.
**Reversed** on the same objection that rules out fixing the iteration count in
D-005: the **station count is a computed result**, and line efficiency is a
function of it. Fixing it leaks part of the answer.

`iteration` and `decision` turn out to be the same shape — an ordered list of
homogeneous sub-traces with stable within-element symbols and a termination
predicate — so specifying them together is cheaper than two bespoke redesigns,
and generalises to `linear_reservoir_routing_step` and `qr_policy_one_iteration`.

---

## D-007 — Answer marker: normalise gold, keep the extractor permissive

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** spec §5

Normalise the five chemical templates to `**Answer:**` rather than widening the
accepted set, which would have to be known by every downstream consumer and
would drift.

**The asymmetry is the point:** strictness applies to *gold traces only*. The
extractor reading **model** output must stay permissive — models emit "Final
answer:", bare headings, or nothing.

---

## D-008 — Item-pool regeneration: split by cost

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** spec §5 decision 5

- **Parser fix → re-score only, immediately** (no items change, no re-inference).
- **Template changes → regenerate once, after Phase 6.** Per-phase regeneration
  would invalidate comparisons repeatedly.

---

## D-009 — Constants re-grounding runs as a parallel Track B, not a Track A phase

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** spec Track B

Only a small slice gates Track A: `CP_PARAMS`, `HEATS_OF_FORMATION` and
`COMBUSTION_REACTIONS` — 65 values in three tables, all consumed by one file
(`heat_effects.py`, 5 templates). That is a 15–25 h job, not the 87–150 h track.

Two hard sync points: **S1** (C2 before Phase 2, so determinism is not verified
against 10×-wrong Cp data) and **S2** (C3 before Phase 6, so the item pool is
regenerated exactly once).

---

## D-010 — `damping_classification` is ill-posed, not wrong (SUPERSEDES D-004)

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** two independent oracle
agents; corroborated against D0.5 §4c

Two agents, working from different briefs, independently concluded the gold
label is **correct**. Both reproduced the audit's ~33.7–35% flip rate under a
strict `zeta == 1` test, then showed **every flip lies within 1.0× the half-unit
uncertainty of the 2-dp damping coefficient** (max ratio 0.99). The stated data
is consistent with `c = c_c` to its own printed precision.

The audit's 35% figure assumes exact float equality as the classification rule.
No solver reading `c` to 2 dp could conclude "Overdamped" from ζ = 1.000000854527.

**Revised fix: state `c` to more digits** (or state it as equal to the critical
value) and drop the fragile `elif zeta == 1` float test. Constraining the
sample (D-004) is no longer needed.

**Also recorded:** D0.5 found audit claim 4c **false as stated** — the flips are
~50/50 Overdamped/Underdamped (89/86), not all Overdamped, and the Underdamped
half is worse than described because gold then omits the ω_d step those
instances require.

---

## D-011 — Phase 1 grows from 9 templates to 12

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** `phase0_baseline.md`
NEW-1 and §6

D0.5 found three more templates with the same defect shape:

| Template | Closure failure |
|---|---|
| `template_beam_deflection_formula` | 5.89% of SI instances |
| `template_statically_indeterminate_shaft` | 47% |
| `template_shaft_design_power` | 30% |

`beam_deflection_formula` should be fixed **first**, because it is the pattern
the others copy (see D-012).

---

## D-012 — P2 amended: break the rounding tie in decimal, not on the binary float

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** `phase0_baseline.md`
NEW-1 / SPEC-2 · **Severity: this invalidated the original Phase 1 recipe**

P2's own worked exemplar, `template_beam_deflection_formula`, applies P2
correctly — `delta = round(delta, 5)` before use — and **still** fails printed-
arithmetic closure on 5.89% of SI instances:

```
0.01185 * 1000  ->  binary 11.85  ->  round() gives 11.8;  decimal half-up gives 11.9
```

Rounding to display precision removes the *double* rounding but leaves an exact
half-way tie that Python resolves on the binary value. **P2 as written is
necessary but not sufficient.**

Applied as originally specified to the other eight Phase 1 templates, the
transformation would have reproduced this residue eight more times. Caught by
the Phase 0 gate before any template was edited — which is the clearest
justification in this project for not compressing Phase 0.

**Amendment:** quantise in decimal — `Decimal(f"{x:.5f}") * 1000` with
`ROUND_HALF_UP` — not `round(x * 1000, 1)`.

---

## D-013 — Four Phase 1 templates are not chain-break defects

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** independent T2 oracles,
reconciled against D0.5

Round-trip oracles (question → answer, computation code withheld) split the
original nine Phase 1 templates into three categories rather than one:

| Category | Templates | Fix |
|---|---|---|
| **Chain break** — the answer is wrong from the stated givens | `mean_variance` (56.8%), `rotating_unbalance` (4.6% @ 0.5% tol, worst 8.43%), `vibration_transmissibility` | P3: state, then recompute forward |
| **Display defect** — the answer round-trips, an intermediate line does not | `cantilever_double_integration` (3.7% of Step-6 lines), `annulus_flowrate` (94% of Step-5 products) | P2 as amended by D-012 |
| **Not defective** | `poissons_ratio`, `logarithmic_decrement`, `system_properties` | none — see below |

**Reconciliation with D0.5.** The oracles and D0.5 appear to disagree; they do
not. D0.5 asked *"does a strict recomputation flip the label / fail to close?"*
and found the audit's rates reproduce. The oracles asked *"can a solver reach
the gold answer from the stated givens?"* and often found yes. Both are correct.
What is contested is **interpretation** — whether a line-level defect makes a
template a chain break — not measurement. An earlier summary of mine framed
this as the audit's claims failing to survive; that was wrong and is corrected
here.

**On the three "not defective":** `poissons_ratio`'s sign leak is a
*presentation* inconsistency (the question states |P| with "(compression)" in
prose while the trace prints a negative), and it affects US-customary instances
too, so the audit's "76/300 SI" undercounts the phenomenon.
`logarithmic_decrement`'s input back-solving is deliberate — the
measurement-uncertainty framing of the exercise. `system_properties`' residue is
two orders below its print precision.

---

## D-014 — R6 review scoping added after a review stalled

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** SPEC-CHANGE 4

The first Phase 0 adversary brief bundled the gate task with re-deriving three
audit claims — work already assigned, in the same batch by the same author, to
the D0.5 agent. It ran 45 minutes past its last write and produced nothing;
D0.5 delivered the superset in 23 minutes. **The only mandatory task went
unanswered because it was bundled with duplicated work.**

R6 now requires: one mandatory gate task, a stated time box, a measurement
ownership register with no number appearing twice, tooling supplied as working
code, and settled questions explicitly fenced off. Re-dispatched under those
rules the review was tractable and found a genuine blind spot (see D-015).

---

## D-015 — Two harness defects that made a check measure nothing

**Date:** 2026-09-05 · **Status:** DECIDED (both fixed) · **Source:** self-caught
and adversary review F1

Recorded because both are the same failure mode — **a check that reports green
while testing nothing** — and both would have silently invalidated a phase gate.

1. **T3 seeded numpy.** `seed_all()` seeded numpy as well as `random`, so both
   of T3's processes reproduced the unseeded `np.random` draw and T3 **passed**
   `template_levenspiel_plot_interpretation` — the one template it exists to
   catch. Fixed: T3's children seed only `random`. Verified 40/40 seeds now
   differ.
2. **T1 was blind to prose-labelled lines.** A colon-prefixed label made
   `_strip_units` treat "Row sum" as an operand symbol, poisoning the segment,
   which then landed in the same bucket as a legitimate formula-only step. The
   adversary planted `Row sum: 0.4347 + 0.2952 + 0.2701 = 1.5000 = 1` — a gross
   non-closure — and the whole suite passed it. Fixed: `:` is a clause
   separator. Verified the defect now fails, the correct version passes, and
   corpus failures are unchanged at 38.

**Lesson to carry:** a green check is not evidence until something known-broken
has been shown to turn it red.

---

---

## D-016 — P2 amended again: remove the tie, do not resolve it (SUPERSEDES the fix half of D-012)

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** Phase 1 implementation,
confirmed independently by Phase 1 Reviewer A (pattern review, finding F-1)
**Severity: this invalidated the amended Phase 1 recipe, as D-012 invalidated
the original one**

D-012 amended P2 to break the rounding tie in decimal:

```python
delta_mm = float((Decimal(f"{delta:.5f}") * 1000).quantize(Decimal("0.1"),
                                                           rounding=ROUND_HALF_UP))
```

**Applied verbatim this is worse than the defect it replaces.** Measured on
`template_beam_deflection_formula`:

```
line      : delta = 0.01185 * 1000 = 11.9 mm      <- D-012's output
evaluated : 11.85                                  <- T1, in binary floats
printed   : 11.9   (tol 0.05, delta 0.050000000000000711)   -> FAIL
```

The pre-D-012 `round()` gives 11.8 and lands at `delta 0.049999999999998934`
-> MARGINAL. So the amendment converts a T1 MARGINAL into a T1 FAIL. Reviewer A
reproduced this at scale — 359 T1 failures per 5,000 seeds on
`beam_deflection_formula` and 301 on `cantilever_double_integration` under the
D-012 recipe — and established the general statement:

> At an exact display tie, `|evaluated - printed| == tol` exactly, and float
> representation error alone decides FAIL versus MARGINAL, **in either rounding
> direction**. No rounding convention can pass T1 on a tie.

That is the point. A half-way tie is not a rounding-convention problem; it is
an **ill-posed instance**. A reader doing decimal arithmetic and applying
half-up reads 0.02325 m as 23.3 mm, a reader using binary floats and `round()`
reads it as 23.2 mm, and the printed line closes for exactly one of them
whichever the template picks. This is the same shape as D-010's finding about
`damping_classification`: an instance sitting on a measure-zero boundary has no
defensible gold answer.

**Amendment — three parts, in order of how much they remove:**

1. **One rounding, at the precision the answer is quoted to, with the upstream
   display precision matched to it** so a decimal unit conversion is exact.
   4 dp in metres IS 1 dp in mm, exactly; 5 dp in metres against 1 dp in mm is
   a 10x precision drop that lands on a tie for ~10% of instances. This removes
   the systematic tie class outright rather than resolving it.

2. **Bind every displayed-and-consumed operand through its own display**
   (`_as_printed(x, spec)` — `float(format(x, spec))`), so the stored value IS
   the value a reader recovers from the printed text. Rounding alone does not
   achieve this: `34.4 * 1e-6` is `3.4399999999999996e-05`, prints as
   `3.4400e-05`, and reparses to `3.44e-05` — a different float. One ulp is
   enough to put the template's value and the reader's value on opposite sides
   of a tie. **This was Reviewer A's finding F-2**, and without it the fix left
   a residual non-closure at roughly 1 instance in 14,000, concentrated on
   particular (section, span) pairs rather than spread uniformly — so it would
   have recurred systematically in any regenerated item pool.

3. **Resample the residual instances that still land exactly on a tie**
   (`_is_display_tie`). Measured rejection rate after (1) and (2): 0.02% on
   `beam_deflection_formula`, 0.115% on `cantilever_double_integration`.
   Reviewer A checked this is not hiding hard cases — the T1 marginal band is
   at the density rounding alone predicts, the last-digit distribution of the
   gold answer stays flat, and rejected values are scattered across the whole
   deflection range. No attainable gold answer becomes unreachable; only the
   ambiguous half-step boundary is removed.

**What this does NOT change.** D-012's diagnosis stands and its `Decimal`
machinery is still used — `_hu()` is decimal half-up throughout, because where
a single rounding is genuinely required it must still not be resolved on the
binary float. What is superseded is the claim that decimal tie-breaking is
*sufficient*.

**Carried forward:** the acceptance evidence for closure must be sized to the
defect rate. A 1-in-14,000 residual has roughly a 7% chance of appearing in a
1,000-seed run (Reviewer A, F-3), so Phase 1 verifies closure at 60,000 seeds,
not the 1,000 the exit gate names.

---

## D-017 — T5b misreads an `%e` format spec as decimal places (harness finding, NOT fixed in Phase 1)

**Date:** 2026-09-05 · **Status:** DECIDED (recorded, deferred) · **Source:**
Phase 1 implementation

`t5_binding.py`'s `_PRECISION_SPEC` is `\.(\d+)[feEg]`, so it reads `.4e` as
"4 decimal places". In an `%e` format the digits are **mantissa** decimals:
`{Ix:.4e}` on `8.49e-05` displays five significant figures, losing nothing,
while T5b computes `round(8.49e-05, 4) = 0.0001` and reports that the display
"discards up to 5e-05".

This makes T5b report a rounding violation against
`template_beam_deflection_formula` and `template_cantilever_double_integration`
that does not exist: `.4e` is lossless in decimal for every `Ix` in
`AISC_W_SHAPES` (verified over the whole table).

**Not fixed here.** The spec forbids modifying the harness to make a fix pass,
and Phase 0 found two checks that were green while measuring nothing (D-015), so
harness edits during an implementation phase are treated with suspicion. Recorded
as a `SPEC-CHANGE` for the harness owner: `_PRECISION_SPEC` must branch on the
presentation type, taking `.Ne`/`.Ng` as significant figures and only `.Nf` as
decimal places.

**Caveat worth keeping.** Reviewer A's F-6 notes that `.4e` is lossless *in
decimal* but not *in binary* — the reparsed float differs by an ulp — so the
harness is complaining about the right operand for the wrong reason. D-016 part
2 fixes the underlying hazard in the templates; the T5b misparse is a separate
harness defect and both are real.

---

## D-018 — Two Phase 1 templates had zero T1 coverage; traces restructured so closure is checkable

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** Phase 1 implementation

`template_mean_variance` reported T1 **coverage 0.0 with 0 checks performed** —
the check passed because it examined nothing. Two causes, both presentational:

1. Each stage was split across lines, so no single line carried both an
   expression and its result (`mu_X = <substitution>` on one line,
   `mu_X = <value>` on the next). T1 cannot link across lines, by design.
2. Operands were joined by juxtaposition, `(-9)(0.172)`, which no evaluator
   reads as a product.

Fixed by putting each stage on one line as a chained equation and using an
explicit `*`. T1 now performs 5,000 checks over 1,000 seeds at coverage 0.375,
all closing. `template_rotating_unbalance` gained the same treatment: 8,000
checks at coverage 0.381.

**This is the D-015 lesson recurring in a new place:** a green check is not
evidence until something known-broken has been shown to turn it red. `run.py`
already reports `coverage` per template; **it should fail, not pass, a template
whose closure coverage is 0** — a template T1 cannot read is unverified, not
verified. Raised as a `SPEC-CHANGE` against the harness for Phase 2.

The remaining uncovered `=` lines in both templates are genuine prose and
formula statements (`mu_X = sum(x_i * P(x_i))`, `Values: X = {...}`), which T1
is correct to skip.

---

## D-019 — T6's median-answer gate is inside its own sampling noise at 1,000 seeds

**Date:** 2026-09-05 · **Status:** DECIDED · **Source:** Phase 1, the only T6
breach raised in the phase

`template_vibration_transmissibility` tripped T6 after its P3 fix:

```
median answer moved 0.411 -> 0.40165 (2.3% > 2%)
```

**This is not a distributional change; it is noise in the statistic.** Measured
on the **unmodified** template, the median of two disjoint halves of the same
1,000-seed range differs by **3.14%** — larger than the "breach" itself. The
true shift converges as the sample grows:

| Seeds | old median | new median | shift |
|---:|---:|---:|---:|
| 1,000 | 0.41100 | 0.40165 | **-2.27%** |
| 5,000 | 0.39400 | 0.39155 | -0.62% |
| 20,000 | 0.38600 | 0.38495 | **-0.27%** |

Per-seed, the fix moves the answer by a median of 6.8e-4 relative, p95 3.9e-3,
with only 26 of 1,000 instances changing by more than 1% — consistent with the
0.48% worst genuine magnitude the Phase 0 oracles measured for this defect. The
distinct-answer count **rose**, 676 -> 928, because the answer is now quoted to
four significant figures instead of three.

**Disposition: not a breach; no sign-off required, and the tolerance is NOT
widened.** The measurement is simply taken at a sample size where the statistic
is stable.

**`SPEC-CHANGE` for the harness.** Phase 0's NEW-7 already required T6's
*distinct-answer* count to be taken at >= 1,000 seeds. The same argument applies
with more force to the **median**, and 1,000 is not enough for it: a template
whose answers span two orders of magnitude has an order statistic that moves
several percent between samples. T6 should either take the median at >= 20,000
seeds, or report a confidence interval and gate on that rather than on a point
estimate. As written, T6's median gate will fire on noise and pass real
regressions of the same size — the failure mode D-015 warns about, in a
different check.


---

## D-020 — `damping_classification`: the strict flip rate is unchanged by design, and cannot be changed without D-004

**Date:** 2026-09-05 · **Status:** DECIDED, but **flagged for Reviewer B**
**Source:** Phase 1 implementation · **Refines D-010**

D-010 prescribed: *state `c` to more digits (or state it as equal to the
critical value) and drop the fragile `elif zeta == 1` float test.* Both were
done. The float-equality test is gone - classification is now an exact
comparison of two scaled integers, and `zeta` is computed for display only and
gates nothing - and for a critically damped instance `c` **is** the critical
value: both are stated to 2 dp and the trace prints
`zeta = 3399.61 / 3399.61 = 1.0000`.

**What did NOT change, and will not: the strict recomputation flip rate is
still 1,676 in 5,000 (33.5%), identical to the pre-fix figure.**

This is not a failed fix. It is a property of the item that no amount of stated
precision can remove: `c_c = 2*sqrt(k*m)` is irrational for essentially every
(k, m) the sampler draws, so **no finite decimal `c` can equal it exactly**.
Stating `c` to 4 or 10 decimals shrinks the gap but never closes it, and a
strict `c > c_c` test flips every critically damped instance regardless.

The only construction that closes it is the one D-010 superseded: constrain the
sample so `k*m` is a perfect square and `c_c` is exactly representable (D-004).
The Phase 1 brief rules that out explicitly.

**Why the item is nevertheless answerable — the measurement that settles it.**
Over 5,000 seeds, computing `zeta` in 50-digit decimal from the stated `m`, `k`
and `c`:

| population | \|zeta - 1\| |
|---|---|
| critically damped (1,676 instances) | max **2.05e-5**, median 2.17e-7 |
| everything else (3,324 instances) | min **0.150** |

The two populations are separated by **more than four orders of magnitude**.
Every critically damped instance has zeta = 1 to at least 4.7 decimal digits;
no other instance comes within 0.15 of 1. Any solver applying any sane rule -
"is zeta equal to 1 to the precision the data supports?" - lands on the gold
label. Only exact float equality fails, and that rule is unanswerable in
principle for this item, which is precisely what "ill-posed, not wrong" meant.

**The 33.5% figure should therefore be retired as an acceptance metric for this
template.** It measures the strictness of the comparison rule, not a property
of the trace, and it will read 33.5% forever. The metric that means something
is the separation above, and the T2 oracle's band test, which passes 1,000/1,000.

**Flagged for Reviewer B (P6 guard).** The judgement that a 7,500x separation
makes the item answerable is a pedagogical call, not a numerical one. If a
domain reviewer disagrees, the remedy is to reopen D-004 and constrain the
sample - which costs the natural-looking parameter values and needs its own
sign-off.


---

## D-021 — T6's distinct-answer gate SATURATES at 1,000 seeds and hid a 29% regression

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 1
**Companion to D-019, which found the same check's median gate is noisy at the
same sample size.**

`template_shaft_design_power` was measured against the committed 1,000-seed
baseline and read **913 -> 822 distinct answers, -10.0%** — passing the "must
not fall by more than 10%" gate by a hair, and reported as a pass.

Re-measured with **both sides at 5,000 seeds** it reads **3,243 -> 2,293,
-29.3%**: a clear breach, and a real one. Rounding the shaft radius to 5 dp and
then doubling it put the diameter on a 0.02 mm grid, so only EVEN hundredths of
a millimetre were reachable and half the item's answer space disappeared.

**A distinct-answer count cannot exceed the seed count.** At N=1,000 a template
with ~3,200 reachable answers reads ~900 whatever its true answer space is, so
the statistic is saturated and the gate measures the sample size rather than
the item. Phase 0's NEW-7 already required this count to be taken at ">= 1,000
seeds"; that floor is too low by a factor of at least five for the templates
that matter.

**Fixed in the template** by carrying the radius at 6 dp, which steps the
diameter by 0.002 mm and makes every hundredth reachable again: 3,243 -> 3,169,
**-2.3%**, comfortably inside the gate.

**`SPEC-CHANGE` for the harness.** T6 must either
1. take the distinct-answer count at a seed count well above the expected
   answer-space size — 5,000 is enough for this corpus, and the committed
   baseline should be regenerated at that N — or
2. report saturation explicitly, e.g. flag any template whose distinct-answer
   count exceeds 50% of the seed count as "count not resolved at this N".

Option 2 is cheaper and safer: it makes the failure visible rather than
relying on everyone remembering to raise N.

**Lesson, and it is the same one as D-015 and D-018:** a check that reports
green while measuring nothing is worse than no check. Here the check reported a
*number*, and the number was an artefact of the instrument's range.


---

## D-022 — Three P6 scoping decisions Phase 1 made and had not written down

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 1 Reviewer B
(physics and pedagogy), finding F-2

Reviewer B's finding is procedural and correct: three changes altered **what
the item states or teaches**, which P6 requires be *recorded and approved*
rather than absorbed into a correctness fix. Two were argued only in commit
messages and one only in a review brief. They are recorded here, and Reviewer B
— the P6 guard, working without sight of the implementer's reasoning — has
signed off all three.

### 1. `annulus_flowrate` — T6 distribution breach, median answer -28%

The only genuine distributional relocation in the phase. Measured at 5,000
seeds per side: p5 2.68e-5 -> 1.30e-5, **p50 2.40e-4 -> 1.72e-4 (-28.1%)**,
p95 2.78e-1 -> 2.85e-1 (+2.7%), distinct answers 4,929 -> 4,940.

**Signed off.** The cause is the removal of an artificial floor, not a change to
the item. On `master`, `if pressure_drop_Pa == 0: pressure_drop_Pa = 10.0` fired
on **17.2%** of instances (Reviewer B's independent replication of the sampler's
draw sequence; the implementer measured 17.8% per draw and 28.4% of instances
including legitimate 10 Pa values). On those instances the sampler's Reynolds
targeting was discarded and the stated pressure drop bore no relation to the
flow described. The shift is entirely a downward extension of the low-Q tail,
which is the signature of un-pinning values the fallback had held at 10 Pa when
the physics asked for 1-4 Pa. Steps, governing equation, answer format and
answer-space size are all unchanged.

Reviewer B additionally established, independently, that **17.1% of `master`'s
instances had Re >= 2100 computed from their own printed answer** — the item
asserted the laminar annulus solution over transitional flow on about one
instance in six. The new guard removes that. Rejection cost: 2.4% of draws.

### 2. `rotating_unbalance` — the question now states that the total mass includes the unbalance mass

Added because the solution assumed it and the question never said so, and a
solver who subtracted `m_e` was marked wrong (Phase 0 oracle finding).

**Signed off, with the trade named:** this is physically correct for the Rao
formulation the template uses, and it removes a genuine ambiguity — but it also
removes a modelling judgement from the solver's task, so the item is
**marginally easier**. That is a P6 trade, accepted because the alternative is
an item whose gold answer depends on an assumption the question does not state.

### 3. `beam_deflection_formula` — the US-customary unit conversion

Reviewer B found (F-1, CONFIRMED, blocking) that the first version of this fix
carried the exact ratio `1.48/12` into the substitution and **deleted the
numeric conversion from Step 2 entirely**, leaving a step titled "Assemble
consistent US customary units" that converted nothing. kip/ft -> kip/in is the
most error-prone operation in that branch, and the SI branch still evaluated
every conversion numerically, so the same item tested different things
depending on which unit system it sampled.

**Not signed off — fixed.** Step 2 now performs the conversion and states its
value to 6 dp (`w = 1.48 kip/ft = 1.48/12 = 0.123333 kip/in`) while the
substitution still carries the exact ratio, so the step demonstrates the
operation and nothing rounded sits between the stated load and the answer.

### The general point, which is worth more than the three entries

**T6 cannot see any of this.** It profiles distinct answers, magnitude
quantiles, step counts, branch proportions and the difficulty label — every one
a property of the *solution*. P6 is about what the item **tests**, which lives
in the **question**. `rotating_unbalance`'s wording change was invisible to
every automated gate in this phase.

`SPEC-CHANGE`, carried to Phase 2: T6 must profile a **question-text hash**
alongside the answer distribution, and any diff that changes question text, a
step's content, or a displayed precision must carry a named DECISIONS entry
before its phase can close. Reviewer B's remark that "five of my findings came
from the instance diff in under ten minutes; none is legible in the 4,000-line
source diff" should also change how the P6 review is briefed: run it against
**before/after instance pairs at matched seeds**, not against the source diff.


---

## D-023 — Phase 1 edits TEN templates, not twelve. The spec says 12, 10 and 9 in three places

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 1 Reviewer A,
spec finding

Three different counts are in circulation for the same phase:

| Where | Count | How it is reached |
|---|---:|---|
| Phase title, D-011, the Phase 1 brief | **12** | 9 original + 3 added by D-011 |
| Spec §1.1a, the authoritative scope table | **10** | 3 chain break + 6 display defect + 1 ill-posed |
| Spec §1.6, the exit gate | **9** | the original nine, never updated |

**Ten is right.** D-011 says Phase 1 "grows from 9 to 12" by adding
`beam_deflection_formula`, `statically_indeterminate_shaft` and
`shaft_design_power` — but it adds them without subtracting the **three that
left in the same table**: D-013 found `poissons_ratio`,
`logarithmic_decrement` and `system_properties` not defective. 9 + 3 - 3 = 10,
which is exactly what §1.1a enumerates and exactly what this phase edited.

The "12" is therefore an arithmetic slip that propagated into the phase title,
the effort estimate and the implementation brief. The "9" in §1.6 is simply
stale.

**Action:** §1.1a stands as authoritative and the phase reports **10 templates
edited**. The deliverable "12 templates edited" is met by the 10 the scope table
names; nothing in scope was skipped. §1.6's "all 9" is superseded.

`SPEC-CHANGE`: amend the Phase 1 heading, §1.4 D1.1, §1.6 and D-011 to read 10,
with a note that the 9 -> 12 arithmetic omitted the three departures.

---

## D-024 — The exit gate's seed count is superseded by D-016 and must be restated

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 1 Reviewer A,
finding F-1 and spec finding

Spec §1.6 gate item 1 reads *"no non-closing printed line at any of 1,000
seeds"*. D-016 supersedes it with 60,000, on the ground that a 1-in-14,000
residual has roughly a 7% chance of appearing in a 1,000-seed run.

**Phase 1 then shipped `annulus_flowrate` verified at 1,000 seeds, and Reviewer
A found six T1 failures in 20,000** — an exact-display-tie on the `kappa` line,
the same class D-016 exists to remove, at a rate of 3e-4. At the gate's own
seed count it had roughly a 26% chance of being seen. It was not seen.

This is the second time in one phase that the acceptance evidence was smaller
than the defect rate it had to resolve (the first was Reviewer A of the pattern
review, finding F-3), and the third time overall that a check reported green
while not resolving the thing it was pointed at (D-015, D-018, D-021).

**Action:** the gate is restated for this phase and carried forward — **closure
is verified at 20,000 seeds minimum, and every template in scope is verified at
that count, not a sample of them.** All ten Phase 1 templates now pass T1 and T2
at 20,000 seeds with zero failures.

`SPEC-CHANGE`: §1.6 gate item 1 must carry the amended number, or a phase can
pass its own written gate while failing the decision that governs it — which is
literally what happened here.

---

## D-025 — The Phase 1 oracles' sensitivity is coupled to the constants tables Track B is about to change

**Date:** 2026-09-06 · **Status:** DECIDED (recorded; action assigned to Track B)
**Source:** Phase 1 Reviewer A, finding F-5

Reviewer A demonstrated two latent P3 leaks that today's constants happen to
hide, by patching the tables in memory:

| Template | If the constants gain one significant figure | Result |
|---|---|---|
| `composite_shafts_series` | `SHEAR_MODULUS_VALUES` at 4 s.f. (77.15 GPa) | worst relative error **1.977e-3** against a declared TOLERANCE of 2e-3 — passes with 1% margin, undetected |
| `beam_deflection_formula` | AISC SI `Ix` at 2 dp | 1.8% of instances fail |

Neither is a defect today: every `SHEAR_MODULUS_VALUES` entry above 20 GPa is
at most 3 s.f., so the template's `_as_printed(G*1e9, ".2e")` is lossless, and
all 14 AISC SI `Ix` values are exactly 1 dp, so the question's `.1f` is
lossless. Both verified over the whole tables.

**It is a coupling, not a bug** — and it is the D-015/D-018 failure mode
arriving through the constants track rather than the harness: a green check
that is green because of an accident of the input data.

**Action, assigned to Track B (C1-C3):** before any constants table lands, add
an assertion that **every constant consumed by a Phase 1 template is exactly
representable at the precision its question states it to**. Track B's sync
points S1 and S2 should carry this as an explicit gate item. If a re-grounded
value needs more digits than the question shows, the template's display
precision must move with it.


---

## D-026 — A third acceptance run, a third defect the previous seed count could not see

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 1 close-out, the
60,000-seed verification D-016 mandates
**Companion to D-024, which recorded the same failure mode one run earlier.**

`template_mean_variance` passed T2 with zero failures at 1,000 seeds and at
20,000. **At 60,000 it fails 3 times**, worst relative error 1.19e-3 against a
declared `TOLERANCE` of 1e-3 (seeds 23265, 35279, 45463).

**The trace is correct in all three.** Seed 23265: the exact variance is
0.358624 and the trace prints 0.359, which is the correct 3-dp rounding. The
failure is the *check's*, not the template's.

**Root cause: the Phase 0 oracle's tolerance was argued from the wrong range.**
Its comment reads *"sigma_X^2 is typically 10-90 … that is at most ~1e-4
relative"*. Measured over 20,000 seeds the variance actually spans **0.284 to
214**, and at the low end one 3-dp display step is **1.76e-3 relative** — above
the declared tolerance. 15 instances in 20,000 fall below 0.5, where 3 dp is
fewer than four significant figures.

This is exactly Reviewer A's finding F-4 — *the declared tolerance is argued,
and the number that actually binds is the display quantisation* — arriving on a
template whose oracle predates the phase.

**Fix: the template, not the oracle.** The spec forbids editing an oracle to
make a fix pass, and the oracle is not wrong to complain: an answer quoted to
three significant figures cannot be round-tripped to 1e-3. The variance is now
quoted to **at least four significant figures** (`var_dp = max(3, 3 -
exponent)`), the same significant-figure treatment already applied to
`rotating_unbalance`, `vibration_transmissibility`, `composite_shafts_series`
and `shaft_design_power` in this phase. Large variances are unchanged; 0.359
becomes 0.3586.

Result at 60,000 seeds: **T1 0, T2 0, worst relative error 4.30e-4** (was
1.19e-3). T6 at 5,000 seeds: distinct answers 2,975 → 4,133 (+38.9%), median
unchanged.

**The mean needs no such treatment** and is left alone: integer values times
exact 3-dp probabilities sum to an exact 3-dp decimal, so the mean is stated
exactly whatever its magnitude, and its measured relative error is 0.

### The lesson, now recorded three times

| Run | Seeds | What it found that the previous run could not |
|---|---:|---|
| Pattern review (Reviewer A, F-2) | 55,000 | 4 non-closing instances after a fix verified at 1,000 |
| Phase close-out (Reviewer A, F-1) | 20,000 | 6 non-closing `annulus_flowrate` instances after a fix verified at 1,000 |
| This one | 60,000 | 3 `mean_variance` round-trip failures after a fix verified at 20,000 |

Each time the acceptance evidence was smaller than the defect rate it had to
resolve, and each time the gate reported green. **The rule that follows is not
"use a bigger number" but "state the smallest rate the run can resolve, and
check it against the rate you are trying to exclude."** A run of N seeds cannot
resolve a defect rarer than about 3/N. The Phase 1 gate is therefore stated as
**60,000 seeds, resolving rates down to ~5e-5**, and any phase quoting a
different number must say what rate it resolves.

`SPEC-CHANGE`: `run.py` should print the resolvable rate alongside its pass
line, so a green result carries its own limit of detection.


---

## D-027 — Branching model: one branch PER PHASE, merged to `master` when the phase completes

**Date:** 2026-09-06 · **Status:** DECIDED (repo owner) · **Source:** raised
during Phase 1 close-out

The spec's working conventions (§0) say:

> All work on branch `redesign/template-integrity`, off `master`. One PR per
> phase.

**That convention was already broken by Phase 0's own merge, before Phase 1
started, and it is not recoverable as written.** `redesign/template-integrity`
sits at `9733a7e` and was merged into `master` at `23a520e`; `master` has since
advanced to `f3688fb`. Branching Phase 1 from `redesign/template-integrity`
would therefore have *discarded* the Phase 0 close-out (`4ccf8b9`) and both
prompt commits. The single-branch model and "one PR per phase" were never
compatible once a phase merged.

The Phase 1 brief resolved this in practice by naming a phase-specific branch —
`redesign/phase1-round-trip` off `master` — which is what Phase 1 used.

**Decision, from the repo owner: merge to `master` after every phase
completes.** The branching model is therefore:

- one branch per phase, named `redesign/phase<N>-<topic>`, taken off `master`;
- merged back to `master` with a `--no-ff` merge commit when the phase's exit
  gate passes, matching how Phase 0 was merged;
- `master` is the integration point every subsequent phase branches from, so
  each phase inherits the last one's close-out.

**Pushing remains a separate decision.** `master` is ahead of `origin/master`
and Phase 1 does not push.

`SPEC-CHANGE`: amend §0 Working conventions to describe the per-phase branch
model. As written it names a branch that is behind `master` and would silently
lose completed work — the third spec inconsistency this phase surfaced,
alongside the template count (D-023) and the exit-gate seed count (D-024,
D-026).


---

## D-028 — `CP_PARAMS` had four defect classes, not one; and the fix introduced a fifth

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2

The spec records one defect: B and C coefficients 10x too large. Correcting only
it left seven rows still wrong.

| Class | Mechanism | Rows |
|---|---|---|
| 1 | column heading (`10^3 B`, `10^6 C`) transcribed as part of the value | most |
| 2 | coefficients matching no source | `O2`, `C2H2`, `C2H5OH`, `C3H6O`, `C6H14` |
| 3 | column shift — `D`'s mantissa written into `C` | `NO2` |
| 4 | gas-phase coefficients under a liquid key | `C6H6(l)`, `C7H8(l)`, `C6H14(l)`, `C3H6O(l)` |

**And a fifth, mine.** Rewriting the table I transcribed `B`'s mantissa as 12.5
where the original read 1.25, for `H2O(l)` and `H2SO4(l)` — the Class-1 defect
committed by the person removing it. Caught by the C2.2 suite, not by review.
An audit now compares every non-replaced row's mantissa against the original.

**Method that worked, and should carry into C3:** two independent checks per
row — NIST's Shomate fit (a *different* functional form, so agreement is
evidence not tautology) and, where published, Smith-Van Ness's own `Cp298/R`
self-check column. The second condemned `O2` outright: the row computes 4.175
against a published 3.535.

**Lesson:** a table that has just been corrected is exactly when a fresh
transcription error is most likely and least likely to be looked for. The suite
must refuse to pass an uncited row rather than skip it.

---

## D-029 — The C2.4 flame-temperature gate quotes the wrong reference

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2

The spec's exit gate asks for methane/air "within the literature range" and
quotes **~2200 K**. `template_adiabatic_flame_temperature` models a single
balanced reaction with **no dissociation**, and for that model the published
value is **2326.35 K**; 2224.25 K is the chemical-equilibrium value *with*
dissociation (ETASR / arXiv:2503.11826, on disk).

The template computes **2311 K** — within **0.6%** of the correct reference.
Every fuel shows the same consistent positive offset against equilibrium
figures, which is the signature of a modelling assumption, not a data defect.

`SPEC-CHANGE`: the C2.4 gate should read 2326 K, or state that ~2200 K is the
equilibrium value and not comparable to this item's model. Judging the item
against 2200 K would have pushed someone to "fix" a correct template.

---

## D-030 — NIST replaces Smith-Van Ness as the primary on-disk source for C2

**Date:** 2026-09-06 · **Status:** DECIDED (reversible) · **Source:** Phase C2

The spec names Smith-Van Ness-Abbott primary with NIST as cross-check. **No
fetchable, citable copy of SVN Table C.1 could be obtained** — accessible copies
are Scribd and SlideShare, which cannot be downloaded and are not legitimate
"edition + page" citations. Fabricating a page number would be the exact failure
this project exists to prevent.

**Decision: NIST Chemistry WebBook (SRD 69) is the primary on-disk source; the
SVN functional form is retained.** 31 species are stored under
`docs/references/nist_webbook/` so every citation resolves to a file.

**Cost, accepted deliberately:** four rows are least-squares fits to NIST tables
rather than transcriptions, and no longer match a textbook table a student might
hold. Residuals 1.1-3.1%, each reported in the row's comment.

**Reversible.** If a citable SVN copy is obtained, those four rows should be
re-transcribed and the fits retired.

---

## D-031 — A dict's key ORDER can be part of the item pool's identity

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2

`template_sensible_heat_temp_dependent_cp` draws its substance with
`random.choice(list(CP_PARAMS.keys()))`. Regrouping `CP_PARAMS` by chemical
family — purely cosmetic, and more readable — changed which substance **274 of
300 seeds drew**.

Order restored, with the reason recorded in the file so the next editor does not
undo it.

**Generalise before C3 and Phase 6:** any constants table consumed by
`random.choice(list(...))` has its ORDER baked into the item pool. C3 touches
~400 values across four branches and will be tempted to tidy exactly this way.
Either freeze key order, or make the templates sample from an explicitly sorted
list so order stops mattering — the second is better, and is a Phase 6 candidate.


---

## D-032 — The real C2 failure mode was an unrecorded validity range, not a bad value

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2 Reviewer H,
finding F-2 · **Hands an item to Track A Phase 2**

Reviewer H's sharpest finding is not about a wrong number. Every value now
reproduces from live NIST. But `CP_PARAMS` is fitted to **1500 K**, and
`template_adiabatic_flame_temperature` integrates Cp to **2844 K** — 90% past
validity — while the C2.2 suite only checked to 1200 K. **The gate was green
because it never looked where the data is actually used.**

Measured against NIST:

| | 1500 K | 2000 K | 2500 K | 2900 K |
|---|---|---|---|---|
| `N2(g)` | −0.5% | +3.2% | +8.2% | +12.5% |
| `CO2(g)` | −0.7% | +3.6% | +8.9% | +13.5% |
| `H2O(g)` | −0.3% | +3.5% | +9.5% | +15.2% |

Reviewer H independently re-solved methane/air AFT twice — once with the
committed table, once from live NIST Shomate alone — getting **2311.1 K** and
**2327.2 K**. The 16 K residual is *entirely* this extrapolation: the table's
Cp runs high at flame temperature, which depresses T.

**Decisions.**

1. **Declare the range.** `CP_VALID_T_MAX` is added to `constants.py`. It did
   not exist before, in any form. A value without its domain is half a fact.
2. **Make the suite look where the data is used.** C2.2 now measures the
   temperature the flame template actually reaches and reports how far past
   validity it runs, with the worst Cp error there.
3. **Do NOT refit the products to 298-3000 K in C2.** A wide-range refit fixes
   the high end (N2 12.5% → 0.7%) but costs low-temperature accuracy, where
   `sensible_heat_temp_dependent_cp` lives (CO2 3.6% → 4.7% at 298 K). That is
   a trade between two consuming templates and belongs with the template owner,
   not the constants table.
4. **Handed to Phase 2**, which owns `adiabatic_flame_temperature`. Its options:
   accept the extrapolation and state it in the trace (textbooks do exactly
   this with Smith-Van Ness tables); restrict the sampled flame temperature to
   the validity range; or carry a separate high-temperature product table. This
   is a real deliverable, not a note — **Phase 2's D2.6 should require it.**

**Generalise before C3.** Reviewer H's words, and they are the most valuable
thing in the review: *"the failure mode this phase actually exhibits is not bad
values but unrecorded validity ranges."* `POWER_LAW_FLUIDS`, `REAL_FLUID_DATA`
and the mechanical `MATERIAL_PROPERTIES` all carry implicit domains — shear
rate, reduced temperature, temper — that a value-by-value check passes and a
template then violates. **C3.2's "order-of-magnitude and cross-property
consistency" checks cannot catch this class at all.** C3 needs a declared
validity domain per table, and an assertion that consuming templates sample
inside it.


---

## D-033 — Two unsourceable species REPLACED rather than shipped behind a tag

**Date:** 2026-09-06 · **Status:** DECIDED (repo owner directed) · **Source:**
Phase C2 completion, after Reviewer H finding F-1

C2 closed with two rows tagged `[KNOWN-DEFECTIVE]`: `H2SO4(l)` (~59% low) and
`CaCO3(s)` (~9% low). Both were measurably wrong and neither could be sourced —
NIST carries no condensed-phase Cp for either, verified against the free and
paid-linked sections of the sulfuric-acid page and both the calcite and
calcium-carbonate pages.

Tagging is honest, but it still ships two wrong numbers into a benchmark. The
repo owner directed: find a replacement source, **or replace the species with
one that has an authoritative citable source.**

| Removed | Replaced by | Source | Fit |
|---|---|---|---|
| `H2SO4(l)` | **`C6H6(l)`** benzene, liquid | four NIST condensed-phase measurements over 293–322 K (Kalali 1987; Grolier & Roux-Desgranges 1993; Reddy 1986; Naziev & Bashirov 1986) | worst 0.49%; Cp₂₉₈ 135.83 vs NIST 135.69 |
| `CaCO3(s)` | **`Al2O3(s)`** corundum, solid | NIST Shomate, α phase, 298–2327 K, CAS 1344-28-1 | worst 0.06% over 298–1200 K; Cp₂₉₈ 78.76 vs NIST 78.80 |

**Result: every one of the 33 `CP_PARAMS` rows now carries a citation, and none
is known-defective.** The C2.2 suite's exclusion list is empty — it went from
209 checks with two rows excluded to **215 checks with none**, which is a
strictly stronger gate.

**Why these two.** Benzene liquid has genuine temperature dependence from real
measurements (Cp rises 134.6 → 139.9 over 293–322 K), so the item stays a
*temperature-dependent* heat-capacity problem rather than becoming a constant-Cp
one; and benzene boils at 353 K, so the template's 280–350 K liquid window stays
in phase. Corundum is a standard refractory with NIST coverage to 2327 K, so it
is comfortably inside validity across the template's whole solid range.

**P6 cost, accepted.** Two substances leave the item pool and two enter.
Sulfuric acid and calcite are arguably more evocative than benzene and alumina —
but a benchmark item built on a number that is 59% wrong tests nothing, and
neither species could be rescued. This also partly answers Reviewer H's F-4: the
liquid sensible-heat pool goes back from two substances to three.

**Reversible.** If a citable source for either original is obtained — Robie &
Hemingway (USGS Bulletin 2131) covers calcite and is freely available, but is
461 pages and could not be resolved to a specific page here — the species can be
restored.


## D-034 — A test may not carry its own answer key

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2 Reviewer G,
findings G-1, G-2, G-3, G-12 (and G's own revision, which found the root cause)

C2 merged with a green 215-check suite that certified three heats of formation
against **nothing on disk at all**. The mechanism, which Reviewer G found only
on a second pass:

`test_chemical_thermochemistry.py` held a hardcoded `DHF_REF` dict of reference
values. `CP_PARAMS` was checked against the on-disk artefact; `HEATS_OF_FORMATION`
was checked against this dict. For `C8H18(g)`, `CH3OH(g)` and `C2H5OH(l)` the
reference value existed **only** in that dict — so the citation said `[ON-DISK]`,
the row said −208.4, the test said −208.4, the deviation was 0.000, and no
artefact anywhere backed either number. **The suite was measuring its own
agreement with itself and reporting it as verification.**

Two aggravating details in the same dict:

- `DHF_TOL_FLOOR = 2.0` floored every tolerance, so three rows passed while
  sitting outside the uncertainty they themselves cited.
- `C2H6(g)`'s uncertainty had been widened from NIST's 0.4 to **0.7** — exactly
  enough to cover its 0.7 deviation.

Neither was deception; both are what happens when the reference values are
editable in the same file as the assertion. That is the point.

**The rule.** A test that verifies a constant against a source must read the
source from the artefact the constant cites. Reference values may not be
transcribed into the test file. `test_citations_resolve.py` property P4
enforces this by refusing any `*REF`/`*_VALUES`/`*_TABLE` dict literal in the
constants-integrity tests, so the failure cannot return quietly.

**What the values turned out to be.** Nothing was actually wrong, which is the
uncomfortable part — the *evidence* was wrong, not the numbers. Every disputed
row matched a specific NIST measurement once all of them were on disk:

| Row | value | matches | the row had cited |
|---|---|---|---|
| `C2H6(g)` | −84.7 | −84.67 ± 0.49 Prosen & Rossini 1945 | −84.0 ± 0.4 (Manion, recommended) → looked 1.75σ out |
| `C3H8(g)` | −103.8 | −103.8 ± 0.59 Prosen & Rossini 1945 | Pittam & Pilcher only was on disk |
| `C8H18(g)` | −208.4 | −208.4 ± 0.67 Prosen & Rossini 1945 | nothing on disk |
| `CH3OH(g)` | −200.7 | −205 ± 10 NIST average of 9 | nothing on disk |
| `C2H5OH(l)` | −277.7 | −276 ± 2 NIST average of 6 | nothing on disk |
| `CH3OH(l)` | −238.6 | bracketed by Baroody −238.4 and Green −238.9 ± 3.6 | Chao & Rossini −239.5 ± 0.2 → 4.5σ out |

The table follows **Prosen & Rossini (1945)** for hydrocarbons wherever NIST
lists it — a self-consistent set from one laboratory. That was a real and
defensible source choice that nobody had written down, so it read as three
errors. It is now stated in the table header.

**One value was genuinely wrong.** `NO2(g)` read 33.2 against the single value
NIST lists, 33.10. Corrected to 33.1. It is consumed by no reaction and no
template, so the item pool is unchanged — but it is the only defect in the
table, and it was found by requiring each row to match a *measurement* rather
than merely name a *species*.

**Generalisation for C3.** "Does the value agree with the reference?" and "does
the reference exist?" are different questions, and C2 shows a suite can answer
the first convincingly while never asking the second. C3 carries ~400 values
across four branches; `test_citations_resolve.py` is a hard entry condition for
it, per Reviewer G, whose exact words were that this is "the only finding here
I would call blocking for C3."


## D-035 — Four provenance classes, because a value can be warranted four ways

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase C2 Reviewer G,
findings G-5, G-6, G-8; spec C1.3

Spec C1.3 defines two tags, `[ON-DISK]` and `[POLICY: sampling-only]`. Phase C2
met at least four distinct evidential situations, and squeezing them into two
tags is what let rows be tagged stronger than their evidence (G-8):

| Class | Warrant | Enforced by |
|---|---|---|
| `[ON-DISK]` | an entry in the cited file whose CAS matches the tag | P2/P3 of `test_citations_resolve.py` |
| `[DERIVED]` | fitted here from on-disk data; range stated | P3 |
| `[BY-DEFINITION]` | true by construction of the scale | P3 |
| `[KNOWN-DEFECTIVE]` | checked and failed; a record, not a claim | not a correctness claim |

`[BY-DEFINITION]` was added when the new test flagged the three elements in
their standard state as `[ON-DISK]` citations naming no CAS. They are 0 because
the enthalpy scale is *defined* that way; pointing at a NIST page for them
claims an artefact is the warrant when it cannot be. Small, but it is the same
error class as the rest of this phase: a tag asserting more than it has.

**The spec's vocabulary should be fixed before C3 tags 400 more values against
it** (SPEC-CHANGE, raised not applied — the spec is not this phase's to edit).

**Also closed here (G-8): three rows advertised weaker evidence than the
artefact already held.** `C2H6(g)` recorded a Smith–Van Ness `Cp₂₉₈/R`
self-check — a check against the source's own column, the exact tautology the
table header disclaims — while a 5-point Gurvich 1989 table sat unused on disk.
`C6H6(g)` and `C7H8(g)` cited *liquid* values on *gas* rows, which justifies the
re-key but verifies nothing about the coefficients. The suite consulted only two
of the six multi-point tables in the artefact; it now consults all of them.
**215 checks → 222, all passing**, and the three rows now rest on real
independent data (worst deviations −2.6%, −4.5%, +3.6%).

**Left open, reported not fixed:** `C3H8(g)` and `C4H10(g)` are verified at
298.15 K alone — the only Cp NIST publishes for them — yet `CP_VALID_T_MAX`
licenses 1500 K. One point licensing a 1200 K extrapolation. The range is the
*source's* claim, not this repo's verification, and the suite now says so on
every run rather than letting the ceiling pass as verified.


## D-036 — D-032 resolved: split the heat-capacity table by consumer

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 2, resolving the
trade Phase C2 handed on (D-032, Reviewer H finding F-2)

Phase C2 found that `CP_PARAMS` is fitted to 1500 K while
`template_adiabatic_flame_temperature` integrates to **2844 K**, and declined to
decide the fix: it is a trade between two consuming templates, and C2 did not
own either of them. Phase 2 owns both. Decided here.

**The cost of the extrapolation, measured.** Solving the same energy balance
with the repo's polynomial against a direct NIST Shomate solve:

| fuel | `CP_PARAMS` | NIST | error |
|---|---:|---:|---:|
| Methane | 2311.2 | 2327.2 | −16.0 |
| Acetylene | 2843.8 | 2908.6 | **−64.7** |
| Carbon monoxide | 2623.2 | 2663.9 | −40.7 |
| Ammonia | 2101.5 | 2106.8 | −5.3 |

Every flame temperature was low, by 0.2% to 2.2%. The reference for methane/air
with no dissociation is **2326.35 K**; the repo computed 2311 K.

**The three options, and why two fail.**

1. **Restrict the sampled range** — *not available*. The flame temperature is an
   **output**, not a sampled input. There is no knob to turn. This is worth
   stating because "restrict the range" is the standard answer to an
   extrapolation problem and it simply does not apply here.
2. **Refit `CP_PARAMS` wide** — fixes this template and **damages the other
   consumer**. `template_sensible_heat_temp_dependent_cp` lives at 298–1000 K,
   where a 298–3000 K fit is materially worse: CO₂'s residual goes from −0.4%
   to −5.4% at 298 K. Trading a correct template for a broken one is not a fix.
3. **Split the table by consumer** — chosen.

**What was added.** `CP_PARAMS_COMBUSTION`: the same Smith–Van Ness functional
form, least-squares refitted to NIST Shomate Cp over 298–3000 K (400 points),
for the four combustion products only — CO₂, H₂O, N₂, O₂. Read **only** by the
flame template. `CP_PARAMS` is byte-identical to before, so
`sensible_heat_temp_dependent_cp` is untouched and cannot regress.

**Result.** All eleven reactions now agree with a direct NIST Shomate solve to
within **1.9 K** (was 5–65 K low). Methane computes **2327 K** against the
2326.35 K reference — **0.03%**.

**The cost, stated.** The wide-range fit is worse at the low end: CO₂ −5.4% at
298 K, O₂ −3.6%. That is the price of the wide range and it is why this is a
*second* table rather than a replacement. It is acceptable here because the
flame integral has almost none of its mass near 298 K — above 1000 K the worst
residual is 1.3% — and the eleven flame temperatures landing within 1.9 K of
NIST is the direct evidence that it does not matter.

**Two tables for four species is a real cost.** A future editor can change one
and not the other. `CP_PARAMS_COMBUSTION` carries a header saying exactly why it
exists and who reads it, and the C2.2 suite checks the flame template stays
inside `CP_COMBUSTION_VALID_T_MAX`. The alternative — one table serving two
consumers with incompatible requirements — is what produced this defect.

**A near miss worth recording.** My first attempt at the evidence above
integrated a *single* Shomate range across a span covering several, and reported
methane at **2191 K** — which would have said the extrapolation made the answer
*better*, and that the honest fix was to leave it alone. The error surfaced only
because Reviewer H had independently derived 2327.2 K in Phase C2 and the
numbers disagreed. NIST publishes the sensible enthalpy in closed form
(`H(T) − H(298.15) = A·t + B·t²/2 + C·t³/3 + D·t⁴/4 − E/t + F − H`, with `F` and
`H` calibrated per range); using it removes the piecewise-continuity trap
entirely. **A decision is only as good as the script behind it, and that script
had no independent check until a prior phase's reviewer supplied one.**


## D-037 — Display ties can be removed by construction, not only by resampling

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 2,
`template_levenspiel_plot_interpretation`

D-016 established that a half-way display tie is **removed, not resolved** — at
a tie no rounding convention closes in both directions, so the instance is
rejected. Phase 1 implemented that as a bounded resample loop, and it worked
there because ties were rare.

It does not work here. Levenspiel's trapezoid terms are built from a table
quoted to 2 dp, and the average height `(y_i + y_i+1)/2` of two 2-dp values
lands on a **half-cent whenever the sum is odd in its last digit** — by
construction, 50% of intervals. With 6–9 intervals per instance, **99.53% of
instances carry at least one tie** (measured, 3000 seeds). There is nothing left
to resample to.

**The tie was an artefact of the display, not of the data.** A half-sum of two
2-dp values is *exact* in 3 dp; that times a 2-dp interval width is *exact* in
5 dp. Printing them at 2 dp and 3 dp was throwing away digits that existed, and
then rounding what remained. Displaying at the precision the values actually
have means **nothing rounds, so no tie can arise** — the tie is removed, as
D-016 requires, but by construction rather than by rejection.

Evidence: the template now reports **zero T1 marginals**, where it previously
carried them on 38% of lines.

**The general rule this adds to D-016.** Before resampling a display tie, ask
whether the quantity is *exactly representable* at a slightly longer display. If
it is, lengthen the display: it removes the tie with no rejection, no change to
the sampled distribution, and no loss of information. Resampling is for
quantities that are genuinely irrational at any finite precision — a square
root, a logarithm, a ratio — where no display makes the rounding go away.

**P6 cost, accepted.** The trace prints `0.09 × 4.250 = 0.38250` where it used
to print `0.09 × 4.25 = 0.383`. Marginally heavier to read; exactly right
instead of approximately right, and a reader who checks the arithmetic now finds
it closes.

---

## D-038 — `iteration` and `decision` share a structure but NOT a comparator; the spec's equivalence is wrong

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, D3.3
**This is a `SPEC-CHANGE` against `template_redesign_spec.md` §3.2.**

The spec says the two node types "are the same shape (an ordered list of
homogeneous sub-traces with stable within-element symbols and a termination
predicate)" and recommends specifying them together because that "is cheaper than
two bespoke redesigns".

**The structural claim is right; the equivalence is wrong, and specifying them as
one type yields a verifier that is wrong on one of the two templates.** They
differ on the property that decides how a candidate trace is *scored*:

> Is the cardinality of the sequence an observable of the answer?

- `normal_depth_iteration` — the answer is the converged depth. The update count
  does not enter it, and is not even **stable under the numerical slack the
  comparator already tolerates**: the termination test is `|Δy| < 0.002 m` on
  4-dp values, so a solver carrying more digits can converge in a different
  number of updates to the same depth. Measured over 4,000 seeds the count is
  `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}` — **32.9% of instances need more
  than three updates**, so this is not a corner case.
- `line_balancing_heuristic` — the answer is
  `(n·CT − Σt)/(n·CT)·100` and `n` **is** the element count. It is also exactly
  determined: every fit decision is an integer comparison, so no rounding slack
  can move it.

So D3.4's comparator must **tolerate** a count mismatch on the first and **fail**
it on the second. One node type cannot express both dispositions.

**The discriminating rule**, stated so it applies to a template this phase never
saw:

> Cardinality is `incidental` when it is sensitive to the numerical slack the
> comparator already tolerates, and `answer_bearing` when it is invariant under
> that slack and appears in the answer. A quantity that is both sensitive to
> slack *and* present in the answer is an ill-posed item, not a node-type
> question.

**Amendment.** §3.2's "specify them together" stands — they share one base
structure, one `carry` mechanism, one verification algorithm. What is struck is
the implication that they share a comparator. The full specification is
[`phase3_node_types.md`](phase3_node_types.md); its §2 carries this reasoning and
its §7 the two comparator dispositions.

**Carried forward.** The spec names `linear_reservoir_routing_step` and
`qr_policy_one_iteration` as generalisation targets. Neither has been fitted, so
the claim that this generalises is **argued, not measured**. Apply the
`incidental`/`answer_bearing` test to them first.

---

## D-039 — The structured trace is a frame local, not a return value; deliberately, and provisionally

**Date:** 2026-09-06 · **Status:** DECIDED (provisional; revisit in the milestone
model) · **Source:** Phase 3, D3.2/D3.3

Both Phase 3 templates now build their trace as a structured node and **render
the printed prose from it**. The node is bound to a local named `trace_nodes`;
the templates' public contract is still `(question, solution)`, and
`tests/trace_schema/extract.py` lifts the node out with `sys.settrace`.

**Why not return it.** Changing the return contract is a corpus-wide decision
affecting all 150 templates and every consumer, not a two-template one. Phase 3's
scope is two templates. Making the change here would have meant either a special
case for two templates or an out-of-scope corpus edit.

**Why render the prose from the node rather than beside it.** This is D-034's
rule applied to traces. If the template built a node *and* formatted the prose
independently, the two could disagree and nothing would notice — the same shape
as a test carrying its own answer key. Rendering from the node makes divergence
impossible by construction, and it is why `extract.py` holds no reference values
of its own. Evidence that the refactor was behaviour-preserving: over 4,000 seeds,
on every seed the screens did not resample, **every arithmetic line of both
solutions is byte-identical** to `master`; the only prose changes are two
deliberately added sentences in the civil item's Steps 3 and 4.

**The cost, stated.** The node is not part of any interface, so nothing outside
`extract.py` can consume it and no check gates its shape. **Promoting it to a
return value belongs to the milestone-model work**, and until then D3.3 is a
specification with two reference implementations rather than a live contract.

---

## D-040 — `normal_depth_iteration` had TWO display-tie populations, not the one the spec listed

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, own measurement

The Phase 3 brief and the spec list one T1 defect for this template: *"T1 has 1
hard failure in 200 seeds (seed 11, the `K = Q*n/S^(1/2)` line)"*. That is
accurate as far as it goes and it is not the whole defect.

Measuring the **exact rational** value of every printed expression against its
display precision over 20,000 seeds — rather than counting T1 hard failures —
finds two tie populations of comparable size:

| Line | Display | Exact ties per instance | T1 hard failures / 20,000 |
|---|---|---|---|
| Step 1, `K = Q*n/S^(1/2)` | 3 dp | **1.14%** | 22 |
| Step 3, the secant update `y_next = …` | 4 dp | **1.02%** | 5 |

**Why counting T1 failures understates it by roughly 10×.** D-016 established
that at an exact tie `|evaluated − printed| == tol` and float representation error
alone decides FAIL versus MARGINAL, *in either rounding direction*. So only a
minority of ties surface as hard failures: 27 hard failures against 427 ties in
the same 20,000 seeds. **A tie that currently reads MARGINAL is the same
ill-posed instance as one that reads FAIL** — it has no defensible gold reading
either way — so the screen removes both populations, and the acceptance evidence
must be the tie census, not the failure count.

**The lesson, generalisable.** T1's FAIL/MARGINAL split is a property of float
representation, not of the defect. Any phase that sizes a display-tie fix by
counting T1 hard failures will under-fix by about an order of magnitude. Measure
the exact-rational tie rate.

**Resolution.** D-037's cheap escape was tested first and does not apply: `K` is a
quotient by a 5-dp square root and the secant update divides by a difference of
3-dp residuals, so neither is exactly representable at any longer display.
Removed by resampling (D-016 part 3). Measured rejection rate **2.30% of seeds**
over 4,000. T1 hard failures over 20,000 seeds: **27 → 0**.

---

## D-041 — Convergence guards that the PROSE relies on must raise, not assert

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, actioning Phase 2
Reviewer C's C-6 and its `-O` observation

Phase 2's Reviewer C observed that `python -O` strips `assert`, so the flame
template's convergence assert "guards development, not production". Phase 2
recorded that and left it, correctly — nothing in that template's prose depended
on the assert holding.

**Here it does.** `normal_depth_iteration` Step 4 told the reader *"The last
change is below the tolerance, so yn = …"* unconditionally, while the only thing
enforcing it was `assert updates <= 5 and …`. Under `-O` a non-converged
iteration would have emitted a trace stating something false about itself. The
budget loop `for _ in range(5)` would simply fall through with an unconverged
depth.

**The rule.** A guard whose failure would make the emitted prose *false* is part
of the output contract and must be an explicit `raise`. A guard that only records
a developer's expectation may stay an `assert`. In this template convergence is
now a `raise`; the sampling-envelope checks (depth within 0.015 m of target,
residual, Froude number) stay asserts.

Step 4 also now **prints the last change** rather than asserting it is small, so
the reader can check the claim rather than take it. That is the same move as
Phase 2's C-2 fix, which put the model tolerance into the flame trace.

**Measured:** 0 non-convergences in 20,000 seeds, so this was latent. The
smallest rate a 20,000-seed run can resolve is ~0.015%; a defect rarer than that
would not have shown, and the guard is what covers the difference.

---

## D-042 — Two Phase 3 defects were latent, not active, and are recorded as such

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, own measurement

The brief listed two `normal_depth_iteration` defects that turn out **not to
fire in the sampled parameter space**. Both were fixed anyway; both are recorded
with their measured rate so a later reader does not mistake a latent defect for a
demonstrated one.

| Defect | Brief's description | Measured | Resolution |
|---|---|---|---|
| `fmt3()` rewrote `"-0.000"` to `"0.000"`, letting the printed operand differ in sign from the stored value (P2) | "A direct P2 violation" | **0 occurrences in 4,000 seeds.** `g` is only near zero at convergence, and the loop exits on the depth change *before* re-evaluating, so `g ∈ (−0.0005, 0)` is essentially unreachable | Helper removed; the **stored** value is normalised instead of the string, so printed and stored agree by construction |
| `g_curr - g_prev` unguarded division | "Reachable in principle; 0 occurrences in 6,000 seeds" | **0 in 20,000 seeds**; smallest `\|g_k − g_{k−1}\|` observed is **0.004**, four steps of the 3-dp residual grid on which `g` lives | Explicit `raise`, not a resample: a vanishing denominator means the sampling box moved, and that should surface rather than be redrawn around |

**Why this is worth an entry.** The brief described the `fmt3` sign discrepancy
as "a direct P2 violation… T5 flags it three times". T5's three flags are real —
they are static findings about a *call in result position* — but they are not
evidence the sign discrepancy ever occurred. It never did. Fixing it was right;
claiming it was an active defect would not have been. The distinction matters
because the phase's acceptance evidence is quoted in the item-pool impact note,
and a latent defect changes no published result.

**Resolution limit, stated per the standing rule (D-024, D-026):** 20,000 seeds
cannot resolve a rate below ~0.015%. Both defects are excluded above that rate
and above nothing below it.

---

## D-043 — Phase 3 does NOT regenerate the T6 baseline; the diff is delivered directly

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, D3.5

T6 fails 142/150 on `master` — re-measured this phase in a `git worktree` at
`00bd0b0` rather than trusted from the brief, and the brief's table reproduced
exactly (T1 30, T2 0, T3 0, T4 3, T5 68, T6 142, T7 83). The committed baseline
is stale corpus-wide, so T6 cannot gate either Phase 3 template.

**Two routes were available, and the cheap one was rejected.** Regenerating the
baseline would turn T6 green, but *a baseline refreshed by the phase it is meant
to gate is not a gate* — and the corpus-wide movement it would absorb has nothing
to do with Phase 3 and would be examined by nobody. Phase 2 faced this and
delivered a direct before/after instance dump instead; Phase 3 does the same.

**What replaces it:** 4,000 seeds per template dumped from a `master` worktree
and from the branch, **each in its own process** (`tests/template_integrity/
phase3_instance_dump.py`, which refuses to run if the module resolves outside the
tree root it was given — the in-process-reload trap has now caught three people).
Full numbers in [`phase3_item_pool_impact.md`](phase3_item_pool_impact.md).

**Regenerating the baseline remains the right corpus-wide action**, as its own
deliberate change on `master` with the movement examined, and is left to Phase 6
where the re-audit owns it. Recorded there rather than done here.

---

## D-044 — D-037's test is "exact for EVERY instance", not "exact for the tie instances"

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3, caught while
justifying two resample screens

D-037 says: *before resampling a display tie, ask whether the quantity is exactly
representable at a slightly longer display. If it is, lengthen the display.*

**Applied naively that test always says yes, and is always wrong.** A value that
lands exactly on a half-way boundary at *k* decimal places is, by definition,
exactly representable at *k+1* places — that is what a tie **is**. So checking
"are the tie instances representable one digit longer?" is circular: the answer
is unconditionally yes, and acting on it relocates the tie population to *k+1*
rather than removing it.

**The test D-037 actually intends** is the one its own worked example satisfies.
Levenspiel's half-sums are exact at 3 dp for **every** instance, so lengthening
the display means **nothing rounds at all** and no tie can arise anywhere. The
question is not about the ties; it is about the whole population.

**Restated, unambiguously:**

> Lengthen the display only if the quantity is exactly representable at that
> display for **every instance the template can draw**, so that no rounding
> occurs at all. If rounding still occurs for some instances, ties are still
> possible and lengthening moves them rather than removing them — then resample.

**Measured on Phase 3's two quantities**, which is what turned this up:

| Quantity | Exactly representable at its display | one digit longer | three longer |
|---|---|---|---|
| `K = Q·n/√S` (normal depth, 3 dp) | 10.9% | 11.3% | 12.9% (6 dp) |
| `100·idle/(n·CT)` (balance delay, 1 dp) | 3.97% (1 dp) | 3.97% (2 dp) | 4.90% (4 dp) |

3,000 seeds each. Both are far from 100% at any reachable display, so both
resample. Note the delay is no more representable at 2 dp than at 1 — lengthening
buys literally nothing there.

**A second reason, specific to answers.** The balance delay's precision is
**stated in the question** ("in percent to one decimal, round half up").
Lengthening an *answer's* display changes what the item asks for, which is a P6
change requiring sign-off; lengthening an *intermediate* quantity's display, as
D-037 did, does not. D-037's rule is safe for intermediates and needs this
caveat for answers.

**What was nearly shipped.** The first draft of both templates' screen comments
justified resampling with "the quotient does not terminate at any fixed number of
places" — true in general, false as stated (a denominator of the form 2^a·5^b
does terminate), and it would not have survived a reviewer with a calculator. The
justification is now the measured population figure above.

---

## D-045 — A display-tie screen rejects the TERMINATING pre-images, so it clusters rather than scatters

**Date:** 2026-09-06 · **Status:** DECIDED · **Source:** Phase 3 Reviewer B,
findings F-B1 and F-B2 (F-B2 convergent with the implementer's own measurement)

D-016 removes an instance that lands exactly on a half-way display boundary,
because at such a tie no rounding convention closes in both directions. The
natural mental picture is that such instances are *arithmetic accidents*,
scattered thinly across the parameter space. **For a quantity built by division,
that picture is wrong, and the error is systematic rather than random.**

A value can land exactly on a half-way boundary at *k* decimal places only if it
**terminates** at *k+1* places. So a tie screen does not sample the parameter
space uniformly — it selects, with probability zero elsewhere, precisely those
pre-images that make the quotient terminating.

**Measured, and it is stark.** `template_normal_depth_iteration` screens
`K = Q·n/√S`, with `S` on a 4-dp grid over `[0.0008, 0.003]`:

| | rejections | distinct slopes |
|---|---|---|
| the `K` screen | 709 in 60,000 draws | **1** — every one at `S = 0.0016` |
| the secant-update screen | 579 in 60,000 draws | 23, spread over the whole window |

`S = 0.0016` is the only value on that grid whose square root is exactly
representable at 5 dp (`0.04`), so `K` reduces to `25·Q·n`, which terminates;
everywhere else `K` is a quotient by a non-terminating root and a tie is
unreachable. The screen removes **25.4% of the instances at that one slope** and
none anywhere else, depleting it from 4.83% of the pool to 3.50%.

*(Note that `√0.0009 = 0.03` is also exact, yet produces no rejections: dividing
by `0.03` multiplies by `100/3` and reintroduces a factor of three, so the
quotient does not terminate. Exact-square-root is necessary, not sufficient.)*

The industrial template shows the same mechanism from the other side: all 913
rejections of its balance-delay screen have `n·CT` divisible by 8, and 85.5% have
`n = 4`, because `n = 4` supplies two extra powers of two and so makes the
percentage terminate more often.

**Why this is worth a decision entry rather than a footnote.** The Phase 3
item-pool note first described the civil rejections as "scattered", supported by
min/median/max over `Q`, `b` and `S` that did span the full box. That measurement
was taken on the **combined** rejection set, where the update screen's genuine
spread across 23 slopes masks the `K` screen's concentration on one. Aggregate
scatter statistics over a union of screens cannot see a single-valued component.
**Decompose a rejection set by screen before claiming it is scattered.**

**Disposition here: accepted, not fixed.** `S = 0.0016` is not an
engineering-distinguished slope — not a steep/mild boundary, not a Froude
threshold — so no class of channel disappears and nothing the item tests changes;
branch movement is 0.10 points against a 5-point tolerance. The alternatives are
worse: excluding `0.0016` from the grid removes a legitimate slope entirely
rather than a quarter of it, and lengthening the display is ruled out by D-044.

**The general rule.** When a display-tie screen guards a quotient, expect the
rejections to concentrate on the divisors that make it terminate, and **check
whether the concentrated parameter is answer-bearing**. Here it is not. On a
template where it were — a screen that removed a quarter of one material, one
support condition, or one flow regime — the same mechanism would be a P6 breach
wearing the costume of an arithmetic clean-up.

## Open decisions

| # | Decision | Needed before |
|---|---|---|
| D-003 | Do the raw `inference_results/` generations still exist? | promising any corrected results table |
| — | Phase 5 scoping: fold into Phase 1 or run as a parallel PR | Phase 1 start |
| — | Whether the 9 self-inconsistent templates are fixed or replaced | Phase 1 start (item-pool ownership call) |
