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


## Open decisions

| # | Decision | Needed before |
|---|---|---|
| D-003 | Do the raw `inference_results/` generations still exist? | promising any corrected results table |
| — | Phase 5 scoping: fold into Phase 1 or run as a parallel PR | Phase 1 start |
| — | Whether the 9 self-inconsistent templates are fixed or replaced | Phase 1 start (item-pool ownership call) |
