# Decision and pivot log

Running record of every decision, reversal and scope change on the template
redesign work. **Append-only.** When a decision is superseded, mark the old
entry `SUPERSEDED` and add a new one — never edit history, because the reason a
decision changed is usually more useful later than the decision itself.

Companion to [`template_redesign_spec.md`](template_redesign_spec.md) (the plan),
[`template_audit_report.md`](audit/template_audit_report.md) (the findings),
[`phase0_baseline.md`](track-a/phase0_baseline.md) (the measurements) and
[`reviews/`](reviews) (independent reviews).

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
sits at `4feb9d3` and was merged into `master` at `d83c794`; `master` has since
advanced to `7b2655c`. Branching Phase 1 from `redesign/template-integrity`
would therefore have *discarded* the Phase 0 close-out (`115cdaf`) and both
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
[`phase3_node_types.md`](track-a/phase3_node_types.md); its §2 carries this reasoning and
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
`1cacc59` rather than trusted from the brief, and the brief's table reproduced
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
Full numbers in [`phase3_item_pool_impact.md`](track-a/phase3_item_pool_impact.md).

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

---

## D-046 — `line_balancing_heuristic` is ~90% shortcuttable, pre-existing, and Reviewer B's own suggestion is what found it

**Date:** 2026-09-06 · **Status:** DECIDED (recorded; the fix is not Phase 3's)
**Source:** Phase 3, actioning Reviewer B's §5 item 1 — *"I scored shortcut rules
I could think of, not a learned upper bound … that is the number I would most
want checked."*

Reviewer B's mandatory gate task was *does the `line_balancing` result still
require the solver to run the heuristic, or has it become a lookup?* B answered
**it is a search, not a lookup**, on strong evidence: the best rule it could
construct from question text reached **64.97%** against a **53.05%**
blind-constant floor, and it proved no better `N_min`-based rule exists by
computing the per-`N_min` majority ceiling (65.55%). It then said explicitly that
a *learned* bound was the one number it could not produce and most wanted
checked, and set its own flip threshold at ~85%.

**Actioning that suggestion flips the answer.** A depth-2 decision tree over
question-text features reaches **90.18%** on seeds held out from training
(8,000 train / 4,000 disjoint test). It is not an opaque model — it reduces to a
rule anyone can state in one line:

> **Count the pairs of tasks that cannot share a station (`t_i + t_j > CT`). If
> five or fewer, answer three stations; if six or more, answer four.**

That bare rule alone scores **90.41%**. Since the answer is
`(n·CT − Σt)/(n·CT)·100` and both `CT` and `Σt` are given, predicting `n` *is*
answering the item. **A solver can score ~90% without executing the greedy rule,
without using the precedence network, and without constructing a single
station.**

**Measured on both trees, and this is what decides the disposition:** `master`
**90.44%**, the Phase 3 branch **90.41%**. Statistically identical.
**Pre-existing. Phase 3 neither caused it nor made it worse.**

**Disposition: recorded, not fixed, and it does not block the Phase 3 gate.**
Phase 3's scope is trace shape. The item's answer-space weakness is an
item-design question owned by whoever owns the item pool, and fixing it means
changing what the item samples — widening the `n` range beyond {3, 4}, or
sampling durations so the pair-count feature stops separating — which is a P6
change requiring sign-off and a regenerated pool.

**Three things this does NOT overturn.**

1. **D3.1's route decision stands, and is strengthened.** The argument for the
   schema route was that `n` is answer-bearing, so fixing it to a constant
   publishes part of the answer. Still true — and an item already 90%
   shortcuttable is one that could least afford it.
2. **D-038's `answer_bearing` classification stands.** `n` being *predictable*
   is unrelated to `n` being *the answer*; the comparator disposition in
   D3.4 §7.3 is driven by the second, not the first.
3. **Phase 3's own deliverable is the partial mitigation.** D3.4 §7.3 separates
   answer credit from process credit for exactly this node type: a candidate
   whose element count matches but whose `committed` sets differ takes the
   answer and **fails the process**. A model shortcutting to the right `n` is
   therefore visible as a shortcut rather than indistinguishable from a solver.
   That is the structured trace doing the job the flat answer cannot.

**What this says about review scoping.** B's finding and B's §5 suggestion
pointed in opposite directions, and B was right about both: right that no
*hand-built* rule clears 65%, right that a learned bound was the missing number.
A reviewer that reports a clean verdict *and* names the measurement that would
overturn it is doing the job R3 exists for. **The lesson for later phases is that
"the best rule I could think of" is a floor on shortcuttability, never a
ceiling** — any future P6 lookup-check should fit a model, not enumerate rules.
Recorded as a `SPEC-CHANGE` against the review protocol's pedagogy brief.

## D-047 — `multipart` is NOT a seventh comparator `kind`; it is a property of the ANSWER

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, Decision 1
**This is a `SPEC-CHANGE` against `template_redesign_spec.md` §4.2/§4.3.**

The spec defines `kind ∈ {numeric, categorical, sequence, symbolic, narrative,
check}` and separately records `answer_type` for all 150 templates, where
`multipart` is **32 of them — nearly a quarter of the corpus, and 24 of the 51
class-C templates.** The six `kind` values do not obviously contain it, and the
brief poses three readings: a seventh kind, a composition of kinds, or a property
of a milestone rather than a kind.

**Decision: the second.** An answer carries an ordered list of `parts`; each part
carries exactly one of the six kinds. The six stand and the **answer schema**
grows.

**The evidence is inside this phase, which is why it is decidable rather than a
preference.** `system_properties_memory_causality` is typed `classification` by
the inventory, described by the spec's own §4.1 table as a *"canonical label
tuple"*, and is in fact **a two-part answer whose parts are both categorical**.
Under a seventh-kind reading it would have to be typed `multipart` — and its
comparator would then have no way to say *"each part is categorical, normalise
each against the categorical vocabulary"*. **The six kinds would stop composing
at exactly the point they are most useful.** Under the composition reading it is
`parts=[categorical, categorical]` and reuses the categorical vocabulary
unchanged, which is what `tests/comparators/kinds.py::compare_label_tuple` does.

**Two structures wear the word, and conflating them is a scoring error in both
directions.** Sampled across the 32:

| `mode` | meaning | example |
|---|---|---|
| `all` | every part required, every part must match | `volumetric_flow_rate` — flow rate **and** average velocity |
| `any` | alternative renderings of one quantity | `undamped_natural_frequency_torsional` — "51.184 rad/s **or** 8.146 Hz" |

Read `any` as `all` and a model answering in rad/s alone is marked wrong; read
`all` as `any` and a model that gets the flow rate right and the velocity wrong
takes full credit.

**Partial credit is reported, never scored.** A partially correct `all` answer is
not a correct answer. A comparator that says otherwise weakens 32 items silently,
which is a larger P6 surface than any template edit this project has made. The
per-part outcomes are always in `observations`.

**Consequence for the class-C claim.** The spec says these comparators "also
serve the 51 class-C templates". The 51 reproduces (the inventory *does* carry an
`instrumentation_class` column, contrary to the brief), and **24 of the 51 are
`multipart`** — so this decision governs nearly half of the class Phase 4 is
supposed to be paying for.

---

## D-048 — Decision 2: route 1 taken for the vocabulary, and its argument FAILS the test

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, Decision 2

The archive holds 2,200 traces but **only 61 for the four Phase 4 templates**:
`memory_causality` 24, `signal_operations` 16, `linearity` 15,
**`incompressible_continuity` 6.** By D-024/D-026's rule — N samples cannot
resolve a variant rarer than ~3/N — the four have resolution limits of 12%, 19%,
20% and **50%**. Six traces are blind to anything occurring in under half of
model outputs.

**Route 1 is taken** (draw from all 2,200) **and the brief's requirement to
*test* the argument rather than assert it was met — the test does not support the
argument as stated.**

**The test.** If these forms are *model habits*, a rule's rate should be driven
more by which model wrote a trace than by which template it answers. Measured
over 27 rules, comparing per-model spread with per-template spread:

| verdict | rules |
|---|---:|
| model habit | **2** |
| item-driven | **15** |
| both / rare everywhere | 10 |

Only `latex-boxed` and `The answer is` are cleanly habits. Every Unicode rule and
most LaTeX rules are item-driven. **A model does not write a superscript because
it is that model; it writes one because the item has an exponent.**

**What that invalidates and what it does not.** It invalidates any *coverage*
claim borrowed across templates, so every coverage number is stated per template
at its own limit and none is averaged. It does **not** invalidate the *rules*: a
normalisation rule is a claim about meaning, not frequency, and one occurrence
anywhere fixes what `\frac{3x^2}{2}` denotes. That narrower claim is the part of
route 1 that survives and the part the rules use.

**A confound, named because it weakens my own test.** The statistic cannot
separate "the model does not use this form" from "the item gave it no
opportunity". Superscripts look item-driven partly because
`incompressible_continuity` is the only Phase 4 template with an exponent.
The 15 is an **upper bound**.

**Residual risk, carried not solved.** The archive cannot say which forms are
*missing* for `incompressible_continuity`. Route 2 needs D-003. For the symbolic
comparator the honest position is close to **route 3**: specified, adversarially
exercised against 22 cases, **not validated against a representative sample**.
Five rules have zero Phase 4 support at all (`unicode-operator`, `latex-text`,
`latex-exponent-brace`, `hedge`, `thousands-separator`) and are named in
`phase4_vocabulary.md` §3.1. **`hedge` — the rule spec §4.4 asks Reviewer B to
guard — rests on a single observation in 2,200.**

---

## D-049 — Decision 3: the comparator is three-valued, and the bias is toward refusing to decide

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, Decision 3

A comparator's failure modes are asymmetric. A **false accept** credits reasoning
that did not happen — which is the AI Tribunal criticism in a new costume, and
the criticism this whole effort exists to answer. A **false reject** penalises a
correct solver and makes the benchmark measure phrasing rather than reasoning.

**Decision: neither is acceptable, so the comparator stops guessing.** Verdicts
are `MATCH` / `MISMATCH` / **`UNRESOLVED`**, and `UNRESOLVED` can never
contribute to a pass. Modelled directly on D3.4 §8B.10's `unchecked` channel: an
unchecked relation may never contribute to a pass, so a node carrying one cannot
report a bare `PASS`.

**Stated bias:** toward rejection over acceptance, and toward `UNRESOLVED` over
both.

**What stops `UNRESOLVED` being a free escape**, which is the obvious objection:

| | definition | effect of `UNRESOLVED` |
|---|---|---|
| precision | `MATCH ∧ correct / MATCH` | none — not in the denominator |
| recall | `MATCH ∧ correct / correct` | **costs recall**, exactly as a wrong `MISMATCH` does |
| decided | `(MATCH + MISMATCH) / all` | reported alongside, never traded against the other two |

A comparator that declines everything scores 100% precision, **0% recall**, 0%
decided. The gate requires ≥95% on both, so the exit gate carries the decision
rather than hiding it, which is what the brief asked for.

**Quantified.** D4.4, 157 hand-written near misses: precision 98.6%, recall
100%, decided 85.4%. Archive, 61 real traces: 100% / 100% / 95.1% — **and that
number is weak evidence, because the archive is the corpus the vocabulary was
derived from.**

**The one accepted false accept is named and kept**: `num-16b`, `7.65 mL`
against a `7.65 L` gold with no declared unit. See D-052.

---

## D-050 — `signal_operations`' gold answer omits the origin on 16.4% of instances

**Date:** 2026-09-07 · **Status:** DECIDED (recorded; not fixed) · **Source:**
Phase 4, D4.1 §4.3

The template renders its answer sequence with an asterisk on the `n = 0` element
— but **only when the result's support contains `n = 0`**. A shift that moves the
support off the origin prints a bare value list, `y[n] = {-1, 5, -4, 3}`, from
which the origin **cannot be recovered**.

Measured over 4,000 seeds: gold states the origin on **3,344 (83.6%)** and is
silent on **656 (16.4%)**.

**This is not a defect being fixed.** The *item* is well posed: the question
states `x[n]` with its origin marked and the transformation is given, so a
solver can derive the indices. What is under-determined is the **printed gold
answer string**, and only for a comparator that reads it in isolation.

**How the comparator handles it**, stated because the handling is the interesting
part: when gold states no origin the answer is the value list alone and is
compared as such; a candidate that pins an origin gold does not is checked for
consistency rather than punished for saying more. When gold *does* pin the origin
and the candidate does not, the comparator first asks whether **any** placement
would match — if none does, the values are wrong whatever the origin and the case
is a decidable `MISMATCH`. That converts 3 of 16 archived traces from undecided
to a defensible verdict and converts none the other way.

**Carried forward.** If a later phase ever wants the answer string to be
self-describing, the fix is to widen the printed support to include `n = 0`. That
changes emitted text on 16.4% of instances and is a P6 event, so it is a Phase 5
or 6 scoping decision, not a Phase 4 edit.

---

## D-051 — `narrative` is specified as undecidable, and its gate is a census

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, D4.1 §4.5

One of the six `kind` values names milestones whose content is prose: a
justification, an interpretation, a physical explanation. There are two ways to
specify a comparator for it and only one of them is compatible with this
project's purpose.

**Rejected: "compare with an LLM."** That puts the AI Tribunal back on the
critical path for exactly the milestones where it is least accountable — the ones
with no checkable answer. The thing five reviewers objected to would return,
wearing a `kind` name.

**Decision: `narrative` always returns `UNRESOLVED`.** Declaring a milestone
`narrative` **removes it from the automatic score** and routes it to whatever
human or model process the benchmark chooses, with that routing **visible in the
results** rather than hidden inside an accuracy number.

**Its gate is therefore not precision but census.** The fraction of milestones
declared `narrative` is a reported quantity, and **a rising fraction is a
benchmark quietly returning to LLM judging.** That number should appear in every
results table.

Conformance: 20 adversarial cases, including one where the candidate restates
gold *verbatim*, and all 20 must return `UNRESOLVED`. A single `MATCH` would mean
the kind had started deciding prose. Currently 20/20.

---

## D-052 — Right number, wrong unit is an accepted false accept; units are checked only when DECLARED

**Date:** 2026-09-07 · **Status:** DECIDED (residual risk) · **Source:** Phase 4,
D4.4 case `num-16b`

D3.4 §10 records "no dimensional checking" as a non-goal: a dimensional
comparator needs units on every symbol, and that is a milestone-model decision.
The residue is a real false accept — **`7.65 mL` matches a `7.65 L` gold** on the
number.

**Two bad options, and the reason neither was taken.** Inferring gold's unit from
its string and requiring the candidate to match it rejects `7.65 litres`, which is
correct — a false reject bought with a false accept. Ignoring units entirely
leaves the hole open with nothing recording it.

**Decision: `compare_numeric` takes an optional declared `unit`.** When gold
declares one, a mismatched unit is a `MISMATCH` and a missing one is
`UNRESOLVED`; when gold declares none, the number alone decides. This closes the
hole **per item, by declaration** — the same "declaration beats inference"
principle that governs the label sets — and leaves it open, **visibly**, for
every item that has not declared.

**The open case is kept in the adversarial corpus as `num-16b` and it is the one
false accept the exit gate carries.** Deleting it would raise D4.4 precision from
98.6% to 100% and the gate would then be concealing a known hole rather than
carrying it. `num-16` (the same answer, unit declared) and `num-16c`
(`7.65 litres`, unit declared, correct) are its companions and show both sides of
the trade.

---

## D-053 — `sympy` is a new dependency, and the symbolic comparator is designed not to need it

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, D4.3

Spec §4.1 names "symbolic equivalence (sympy)" and §4.2.4 asks for its
"timeout/failure behaviour". **`sympy` was not installed and is not in
`requirements.txt`.** It has been installed for this work; adding it to
`requirements.txt` is a repository-level decision left to whoever owns that file.

**The comparator is built so that the dependency does not gate this phase's
template.** Every answer `incompressible_continuity` produces is a bivariate
polynomial with rational coefficients and degree ≤ 2, and equality there is
decidable by expanding to a coefficient map: exact, local, no timeout, and **no
`UNRESOLVED` outcome**. That matters more than the convenience — this is the
template with **six** archived traces, and a comparator whose verdict depended on
whether an optional package happened to be installed would make the thinnest
evidence in the phase thinner still.

A CAS is genuinely needed for the general kind: the corpus's other eight
`symbolic` templates carry sinc, exp and Q-functions. There, an import failure, a
parse failure, a timeout, or a symbol outside the declared alphabet all return
`UNRESOLVED` — **never `MATCH`**. Nothing about a comparator failing to parse an
answer is evidence the answer is right.

---

## D-054 — D3.4 §7 has four defects, all found by exercising it for the first time

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, D4.6
**This is a `SPEC-CHANGE` against `phase3_node_types.md` §7.**

Phase 3 shipped §7 unexercised and said so; both schema reviewers named the
missing corpus as the largest gap. Building it (16 candidate traces + 9 tolerance
assertions, all 15 dispositions covered) found four defects **in the
specification**.

**F7-1 — §7.2's headline disposition has no instances.** *"A correct answer
reached in four iterations instead of three is correct"* cannot occur for
`normal_depth_iteration`. §4.3.3 ties `converged` to the node's own tolerance,
§8A.4 forbids `converged: true` before the last element, §4.3.1 recomputes every
update, §4.6 recomputes every frame, and the preamble is fixed by the question.
Together **these make the conforming node for a given question unique.** Verified
both ways: clearing `converged` on the terminating element fails 4.3.3 and 8A.4;
appending after it fails 8A.3 and 8A.4.

D-038 measured the count varying 1..5 over 4,000 seeds and that measurement
stands — **but the variation is between QUESTIONS, not between solvers of one
question**, and §7.2 conflates the two.

> **This is the six review rounds' bill arriving.** Each round closed a way for a
> candidate to lie. Together they also closed every way for a candidate to be
> *differently right*. Nobody recorded that trade, because §7 was never run.

**F7-2 — §7 never says whether §6 step 8 binds a candidate.** It does in every
verifier built, and it should: a model whose stated answer contradicts its own
steps has failed the process. The consequence is that answer credit and process
credit are **not independent** for `iteration`, so "wrong answer, clean process"
is not constructible. §7 is written as though they vary freely.

**F7-3 — §7.4's direct solver cannot be an `iteration` node at all**, because a
direct solve has no iteration. Its tolerance is correct and testable only at the
function level, where the corpus now tests it (9 assertions).

**F7-4 — §7.3 row 3's set-versus-ordered distinction is vacuous.** §5.4.6 already
requires `committed` to equal the ordered chosen values, so a candidate whose set
matches but whose order differs is rejected before §7.3 runs.

**D4.9 reaches F7-1 from the other direction and on most of the corpus:** 141 of
150 templates emit a trace whose length is a **constant**, so `cardinality` has
one possible value and both §7.2 and §7.3 dispositions are vacuous there too.

---

## D-055 — The Phase 3 node types are not under-applied; they fit the only two templates that have their shape

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4, D4.7 and D4.8

Both deliverables were written on the assumption that `iteration` and `decision`
generalise and simply had not been tried. **Measured, they do not.**

**D4.7 — a third `iteration` template, binding only, no verifier edit.** Both
templates the spec names were attempted against an unmodified
`reviewer_d2_verifier_15`. Both **rejected**, for different reasons, which is what
makes the pair informative:

- `linear_reservoir_routing_step`: §4 fixes `termination.kind` to `"convergence"`
  and routing **exhausts a two-interval hydrograph** instead. Its element count
  is neither `incidental` (a third interval gives a different answer) nor
  `answer_bearing` (the count is not the answer, `O3` is), so **D-038's
  discriminating rule has no verdict for it**. And `FrameEval` admits only
  constants, the iterate, and earlier frame symbols — there is nowhere to put the
  per-element inflow pair. The verifier says so itself: *"'coeff' is not a role of
  this node"*.
- `qr_policy_one_iteration`: iterates on a **pair** coupled through `n(R)`, and
  §4.1 declares exactly one iterate triple.

**D4.8 — a second `decision` template.** All 150 scanned; 19 carry
decision-shaped vocabulary and **none has the shape**. A `decision` node needs
three things at once — an ordered list of commitments, a shared budget they
consume, and a count that is part of the answer — and `line_balancing_heuristic`
is the only template with all three.

**Consequence.** D3.3 §10's two limitations **stand and cannot be closed by this
corpus**: `decision` remains exercised on one instance and one precedence DAG.
Closing them requires a **new item**, which is an item-design decision and not
Phase 4's to take.

**The reframing, which is the useful part.** `iteration` is not "a repeated
sub-chain"; it is *"a convergence-terminated refinement of one quantity driven by
its own residual"*. D3.3 §10 should say so, and should retire
`linear_reservoir_routing_step` and `qr_policy_one_iteration` as named
generalisation targets rather than leaving them as unfinished work.

---

## D-056 — The hedge policy is ADVISORY: detected, annotated, never scored

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 4 Reviewer E,
round 4, recommendation 3c
**This is a `SPEC-CHANGE` against `template_redesign_spec.md` §4.2 and §4.4.**

§4.2 requires label normalisation to handle *"hedging ('appears to be linear')"*
and §4.4 names *"accepting a hedge that never commits"* as the specific case
Reviewer B must guard. Versions 1–3 of `tests/comparators/commitment.py`
implemented that literally: a hedge made the answer `UNRESOLVED`.

**Four review rounds say the enforcement is not worth what it costs.**

Reviewer E's round-4 census splits a number the phase had been treating as one:

| | count over 2,200 archived answer spans |
|---|---:|
| hedge markers fired | **2** (and **0** in a commitment-gated answer type) |
| subordinator openers | **43** (concessive 31, hypothetical 11) |

Two mechanisms inside one module with **opposite evidence**. E then ablated the
hedge layer, leaving segmentation untouched:

| | full | ablated |
|---|---|---|
| recall-corpus positive frames | 351/351 | **351/351** |
| recall-corpus total | 507/507 | 477/507 |

**The entire hedge-governance layer buys 30 synthetic negative controls that the
reviewers themselves wrote, and zero archive verdicts and zero positive recall.**
It is ~250 lines and **13 of E's 20 findings across four rounds**, and it
repeatedly marked *correct* answers wrong — which D4.1 §1 states is as
unacceptable as a false accept.

**Decision.** `HEDGE_POLICY = "advisory"`. A hedge is still detected, still
named in `observations`, still counted by the census — and never converts a
`MATCH` into an `UNRESOLVED`. `ENGTRACE_HEDGE_POLICY=enforce` restores the old
behaviour, and flipping that constant is the whole change.

**This is the same move D-051 makes for `narrative`**: a property the comparator
cannot decide reliably becomes a *reported quantity* rather than a silent
verdict. If the hedge rate ever becomes non-trivial, the census makes it visible
and the evidence for enforcing will exist.

**The cost, stated rather than absorbed.** `cat-16` — a hedge naming the right
label — is a **false accept** under the shipped default, and it is exactly the
case Reviewer B was commissioned to guard. It stays in the adversarial corpus
beside `num-16b` so the exit gate carries the trade. D4.4 precision moves
98.6% → 97.2% against a 95% gate.

**Why this is a P6 decision and not an implementation detail.** It changes what
counts as a correct answer, corpus-wide, for `categorical`, `categorical[tuple]`
and `check`. Reviewer B's position (a non-committal answer earns no credit) and
Reviewer E's measurement (the layer enforcing it has no observed instances and
produces false rejects) are both right, and this resolves in E's favour **on the
evidence**, not on preference. B did not object: its round-4 verdict is that the
phase can close.

---

## D-057 — Both Phase 4 classification items are 100% shortcuttable; recorded as a P6 scoping decision, not fixed

**Date:** 2026-09-07 · **Status:** DECIDED (recorded; item redesign out of
scope) · **Source:** Phase 4 Reviewer B, rounds 1 and 4, under SPEC-CHANGE 10

SPEC-CHANGE 10 requires a pedagogy lookup-check to **fit a model, not enumerate
rules**, and to report lift over a blind-guess floor. Reviewer B did:

| item | blind-guess floor | depth-2 held-out | lift |
|---|---:|---:|---:|
| `system_property_linearity` | 50.02% | **100.00%** | **+49.98** |
| `system_properties_memory_causality` | 34.02% | **100.00%** | **+65.98** |

Not overfitting: the generators emit **5** and **7** distinct right-hand-side
shapes and **no shape ever carries two labels**, so the form→label map is a
total function. 100% is the ceiling and the floor at once. The linearity rule
reduces to *"is there a `*`, or an `x[n -`?"*, and the template's docstring
claims it tests additivity and homogeneity.

**Why this is Phase 4's to record and not Phase 4's to fix.** Spec §4.1 forbids
redesigning these templates — *"the schema adapts, not the items"*. But Phase 4
is the phase that decides they are scored **on their answer alone**, and P6
requires a change to what an item tests to be recorded and approved rather than
arrive as a side effect. So it is recorded.

**Reviewer B withdrew its own proposed remedy, on measurement.** B had suggested
extending `check`'s "score the supporting quantity too" pattern to `categorical`
— require a linearity answer to name *which* property failed. B then measured
that only **4 of 15** archived linearity answer spans name one, so requiring it
would drop recall to ~27% and breach the §4.5 gate; it is also asymmetric, since
only `not linear` can carry a failing property, which *adds* a cue. **Filed as a
Phase 5 recommendation**, where the item can change with its gold.

**Actioned now:** `blind_guess_floor` and `surface_model_heldout` are two new
columns in [`template_inventory.csv`](audit/template_inventory.csv), computed for
every template whose answer space is small enough for the statistic to mean
anything (5 of 150; the rest are numeric and have hundreds of distinct answers,
which is the correct reason to leave them blank). B's threshold: flag at **lift
≥ 40 points or held-out ≥ 95% regardless of lift**, because a high floor masks a
total lookup. **Three templates are over it.** A flag means "this needs a P6
register entry", not "block the template".

**A first attempt at these columns was wrong and is worth recording.** It
labelled each instance by its digit-masked answer *shape*, which collapses every
scalar template to one class and reported **112 of 145** templates as fully
shortcuttable. The unmasked answer is the right label, and a numeric template
then drops out on its own. Caught by disbelieving the number, not by a check.


---

## D-058 — Cross-pairing finds three comparator defects Phase 4's corpora could not, and one binding gap that is larger than all three

**Date:** 2026-09-07 · **Status:** DECIDED (recorded; fixes assigned to Phase 5)
**Source:** post-merge experiment answering Reviewer E's F0 and R2-F10

Phase 4 closed with two gate items **met and uninformative**: the archive holds
no wrong answer at all for `categorical` or `categorical[tuple]`, and no trace of
any kind for `numeric`, `check` or `narrative`. Reviewer E said precision without
a per-kind negative count is unreadable. This is the measurement that fills it.

**The instrument.** Gold is by definition the correct answer to its own question,
so pairing gold *A* as a **candidate** against gold *B* as the **gold** has a
truth known with no label at all: match iff the two answers are identical. That
manufactures negatives for every kind, at arbitrary scale, from the templates
alone. Run over all 150 templates × 12 instances: **19,668 pairs.**

### Three real defects

**N1 — `parse_number` reads the first number in the answer span, which is often
not the answer.** Confirmed false accept:

```
cand: The volumetric flow rate of Engine Oil (SAE 50) ... is 0.001183 m^3/s.
gold: The volumetric flow rate of Engine Oil (SAE 50) ... is 0.002738 m^3/s.
-> MATCH        both canonicalise to 50
```

The incidental number is a fluid grade. `gas_viscosity_kinetic_theory` shows the
same defect reading a **temperature** (541) and a **subscript digit from a
chemical formula** (the 4 of C₄H₁₀). **3 scalar templates produce a false accept
in a 12-instance sample**; 58 of 150 carry a leading incidental number and are
therefore exposed to it.

**This is D-003's defect class reappearing inside the comparator built to replace
it**, and no corpus Phase 4 shipped could see it, because `numeric` had zero
archived negative instances.

**N2 — an empty parse compares equal to an empty parse.** `_as_polynomial`
returns an empty coefficient map when it recovers no polynomial content, and
`{} == {}` is a `MATCH`. On `autocorrelation_rect_pulse` that is **128 of 132
pairs**: every expression the parser fails on matches every other one it fails
on. **An unrecoverable parse must be `UNRESOLVED`** — this is D4.1 §4.4's own
stated rule ("failure is `UNRESOLVED`, never `MATCH`") violated by a code path
that never reaches the CAS.

**N3 — the comparator has bindings for 4 of 150 templates.** D4.1 specifies the
contract and `PHASE4_BINDINGS` implements it for the four Phase 4 items only.
Everything else has no declared label set, no declared unit, no declared answer
shape. Scored with defaults, `euclidean_distance_binary` matches **132 of 132**
pairs because a multipart answer read by a single-part comparator returns the
`1` from `s1`.

**N3 is larger than N1 and N2 together, and it reframes Phase 5.** The
comparator does not need more rules; it needs **bindings**, and each binding
needs evidence that it discriminates. Cross-pairing is that evidence.

### What is NOT a defect, stated because the same run reports it

The sweep's 68 `categorical` and 136 `sequence` false *rejects*, and most of its
560 false accepts, are **artefacts of my harness**, not comparator defects: it
scored every template with *default* options, so `damping_classification` was
judged against the **linearity** label set and every `multipart` answer was
judged by the `numeric` comparator. Those numbers say only that N3 is real. The
three defects above are the ones that survive per-kind binding.

**An earlier framing of mine was wrong and is corrected here.** I first reported
"58 of 150 templates" as the size of N1. 58 is the count *exposed* to it — having
a leading incidental number is necessary and not sufficient. The measured count
that produces a false accept in a 12-instance sample is **3**. Both numbers are
useful and they are not the same number.

### Standing instruments, not one-off experiments

Reviewer E said of its own cross-pairing that it *"should be run as a standing
check rather than as a one-off review artefact"*. Two complementary corpora, both
assigned to Phase 5:

| instrument | scale | tests |
|---|---:|---|
| **gold × gold** | 19,668 pairs, all 150 templates, no archive needed | equality and discrimination, every kind |
| **archive × gold** (E's) | **22,982** real-text pairs available — 11,384 `scalar`, 5,511 `multipart` | the same, in real model surface forms |

Gold×gold cannot test normalisation, because gold text lacks model surface
variation; archive×gold cannot cover kinds the archive lacks. Neither replaces
the other and Phase 5 lands both.

## D-059 — One canonical answer marker in gold; the candidate-side list stays wide

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track A, D5.3

Five templates terminated with `**Final Answer**` (four) or `**Final Answers:**`
(one) instead of the canonical `**Answer:**`. The brief offered two remedies —
normalise the templates, or widen the accepted marker set corpus-wide — and
required the decision to be made once and applied everywhere.

**Decision: normalise the five templates. Gold emits exactly one answer marker.**

**The two lists are different lists, and conflating them is the mistake this
entry exists to prevent.** `tests/template_integrity/core.py::ANSWER_MARKERS` is
what *gold* may emit; `tests/comparators/normalize.py::ANSWER_MARKERS` is what a
*candidate* may emit. Gold is a corpus we control and can make uniform. Model
output is not, and never will be: the archive shows `## Final Answer`,
`Final Answer:`, `The answer is:` and a bare `Answer:` across 2,200 traces. So
the candidate-side list stays wide and keeps its priority order; the gold-side
one narrows to a single marker.

**"Agree" is a predicate, not a sentiment, and it is now checked.** For every
gold solution in the corpus the candidate-side parser must recover the canonical
marker with no marker debris glued to the front of the span. Run:

```
python -m tests.template_integrity.phase5_contract_scan --markers
```

**Measured, before and after, over the eleven Track A templates × 200 seeds:**

| | spans recovering `**Answer:**` | spans carrying debris |
|---|---:|---:|
| before | 1,200 / 2,200 | **1,000** |
| after | 2,200 / 2,200 | **0** |

Corpus-wide after the edit: **1,200 of 1,200** gold spans (150 templates × 8
seeds) recover `**Answer:**`, none empty.

**This was not cosmetic and the "widen the set" option would not have fixed it.**
The candidate-side priority-2 pattern `Final\s+Answer:?` matched *inside*
`**Final Answer**`, so the recovered marker was `Final Answer` and the answer
span began `**\nAfter reaching a conversion of...`. On the plural
`**Final Answers:**` it began **`s:**\na) CSTR Volume = 13...`**. A comparator
reading those spans was reading `s:**` as part of the answer. Widening the
accepted set would have licensed the collision permanently; the priority order
is load-bearing (Reviewer E's Phase 4 F4 was this same mechanism firing inside a
caveat) and every marker added to it is another way for `answer_span` to peel
into the wrong place.

**P6.** The five templates' answer *bodies* are byte-identical before and after
on all 2,000 seeds each; only the marker line differs. Questions unchanged.
Evidence in [`phase5_item_pool_impact.md`](track-a/phase5_item_pool_impact.md).

**Reviewer A owns checking this**, and is asked to treat `--markers` as a claim
under test rather than as supplied tooling.

---

## D-060 — T4's non-canonical-marker finding becomes a FAILURE (SPEC-CHANGE 14)

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track A
**This is a `SPEC-CHANGE` against `template_redesign_spec.md` Phase 5 Track A's
exit gate.**

T4 printed `non-canonical marker {'**Final Answer**': 25}` for five templates
**and passed them**. The corpus could therefore be reported "T4 150/150 green"
with a known output-contract defect standing on five items — the third instance
in this project of a check that is green while measuring something it can see.

**The distinction that decides the remedy, because an earlier draft of the
Phase 5 brief got it wrong.** *Ungated* means the check sees the defect and does
not fail on it: the fix is one line of severity. *Invisible* means the check
cannot see it: only then is a new instrument warranted. Five templates were
ungated and three were invisible; treating the first group as the second would
have built a duplicate scanner.

**Decision.** `ContractResult.passed` now includes `non_canonical_answer`.
The *recognised* marker set stays wide deliberately, so a regression is reported
as `non-canonical marker {'**Final Answer**': 25}` rather than as the far less
useful `no answer marker on 25 seeds`. Narrowing the recognised set would gate
the same defect with a worse diagnostic.

**The severity change is verified by a planted defect, not by a green corpus.**
After Track A's edits T4 is 150/150 green *whether or not this change was made* —
the corpus is clean because the templates were fixed. A green suite is therefore
no evidence at all for this line. `--selftest` plants each of the four defect
classes into real generated output and requires detection:

```
python -m tests.template_integrity.phase5_contract_scan --selftest
```

| planted class | scan | T4 |
|---|---|---|
| `step_marker` | CAUGHT | **FAILS** |
| `answer_marker` | CAUGHT | **FAILS** ← the line this entry adds |
| `complex_sign` | CAUGHT | passes (expected: blind) |
| `degenerate_product` | CAUGHT | passes (expected: blind) |
| `degenerate_product` planted in a *derivation* step | silent (correct) | — |

---

## D-061 — The doubled-sign defect is 16 templates, not 2; Track A fixes its 2 and names the rest

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track A,
corpus-wide sweep

The spec records the malformed-complex defect as `+ j-51.22` on two templates.
The underlying fault is more general: **a hard-coded `+` in a format string
followed by an interpolated signed value.** The `j` is incidental.

Swept corpus-wide, 150 templates × 120 seeds, pattern
`[-+]\s+[-+]\s*\d` or `[-+]\s*j\s*[-+]\s*\d`:

| | templates |
|---|---:|
| emit a doubled sign anywhere in the solution | **16** |
| emit one **inside the answer span** | **4** |
| of those, in Track A's scope | **2** (`time_to_phasor`, `phasor_addition`) |

The other two answer-span carriers are `lorentz_force` (118 answer-span hits per
120 seeds) and `continuous_to_discrete_conversion` (66). The twelve
derivation-only carriers are `multi_segment_rod`, `vorticity_check`,
`nyquist_rate_determination`, `impulse_response_from_lccde`,
`work_isothermal_virial`, `mean_variance`, `coulombs_law`,
`wave_equation_interpretation`, `system_property_linearity`,
`signal_operations`, `sensible_heat_constant_cp`, `pitzer_correlation_z`.

**Decision: fix the two in scope, completely; record the other fourteen with the
measurement and assign them.** Editing fourteen unassigned templates inside a
phase whose named risk is *"changing an item pool by accident"* would trade the
thing the phase is for. Disposition: **`ADOPT-PHASE-6`**, added to Phase 6's
deliverable list.

**"Completely" was larger than the spec's row.** Fixing only the `j` sites in
the two templates would have left `cos(361*t + -146.39 deg)` standing in
`phasor_addition`'s **answer**. The fix therefore covers every doubled sign in
both templates, and a third defect found while doing it:

> `A_total = sqrt(-41.2^2 + 86.45^2) = 95.77` — evaluated as printed this is
> `sqrt(-1697.4 + 7473.6) = 76.00`, because `-41.2^2` is `-(41.2^2)` under
> ordinary precedence. **A P1 violation**, printed on every instance with a
> negative component. Parenthesising the operands cut `phasor_addition`'s T1
> closure failures from **9 to 5** over 25 seeds and removed the large-delta
> class entirely — the surviving five are ordinary rounding-boundary misses
> (delta ≈ 0.006 against a 0.005 tolerance). T1 coverage on the template rose
> 0.27 → 0.29 and lines checked 142 → 152.

That P1 defect was *not* in the brief, the spec, or the audit. It was found
because the doubled-sign sweep put the line in front of me.

---

## D-062 — Reviewer A's triage: every detector was narrower than the class it was named after

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Reviewer A
(correctness), `reviews/phase5_reviewer_a_correctness.md`

**Verdict `PASS WITH FINDINGS`, and the shape of the findings is the result.**
Nine findings, all CONFIRMED, and **not one is about the corpus.** The reviewer
rewrote all four detectors from the spec's defect table, swept 150 × 400 before
and after, and found the corpus genuinely clean. What it broke was the *gate*:

> The pattern this review keeps finding is *the detector is narrower than the
> class it is named after* — F1, F3, F5 and F6 are all instances.

That is worth more than the individual fixes. The phase had the right instinct —
`--selftest` plants a defect per class — and the instinct was defeated by
**writing the plant and the regex with the same hand**, so every plant was a
shape its regex already matched. All four narrow detectors passed the self-test.

### Findings

| # | Finding | Disposition | Where actioned |
|---|---|---|---|
| **F1** | `complex_sign` requires a literal `j`, so it gates only half the class D-061 defines — the half that did *not* change 2,560 question strings | `ADOPT-NOW` | `_DOUBLED_SIGN` added; gated on Track A's eleven, census corpus-wide (14 templates) |
| **F2** | The scan never reads `inst.question`; nothing in T1–T7 does either except for emptiness | `ADOPT-NOW` | `scan_solution(sol, res, question)`; classes 3–4 read the question, classes 1–2 deliberately do not |
| **F3** | `_DEGENERATE_PRODUCT` cannot match `0.0*pi`, and its census was short by one template | `ADOPT-NOW` | regex takes `0(?:\.0+)?`; census 2 → 3 templates |
| **F4** | `phase5_contract_scan` is in no standing gate — not in `ALL_CHECKS`, not imported by `run.py`, no CI in the repo | `ADOPT-NOW` | **T8** (`checks/t8_emission.py`), in `ALL_CHECKS` *and* `DEFAULT_CHECKS`, at 400 seeds |
| **F5** | A step heading that loses its bold entirely is invisible, and T4's contiguity check passes when the lost heading is the **last** one | `ADOPT-NOW` | `_STEP_ANY` is line-oriented, not marker-oriented |
| **F6** | A non-canonical marker *added beside* the canonical one passes T4 and the scan — and is priority 0 on the candidate side | `ADOPT-NOW` | `## Final Answer` recognised; T4 counts every marker, not only when the canonical one is absent |
| **F7** | `levenspiel_plot_interpretation`'s answer span swallows a 363-char `**Note:**` block whose last number is `31.8` | `ADOPT-NOW` | the Note moves **above** the answer marker; span 363 → 48 chars |
| **F8** | "7 of 150 templates give a different T6 report" is not reproducible | `ADOPT-NOW` | corrected in the item-pool note; the claim now rests on the order-insensitive result alone |
| **F9** | "T1–T7 unchanged per template" was a verdict claim; T5's `operand_restatements` moved 0→5 / 0→1, recorded nowhere | `ADOPT-NOW` | recorded in the item-pool note |

### §5 suggestions

| # | Suggestion | Disposition |
|---|---|---|
| 1 | Wire the detectors into `run.py` as T8 at 400 seeds | `ADOPT-NOW` — done; also added to `DEFAULT_CHECKS`, which the reviewer did not ask for and which is where it will actually run |
| 2 | Widen `_COMPLEX_SIGN` to the class, census the residual | `ADOPT-NOW` — the reviewer's own regex, adopted verbatim in substance |
| 3 | Scan the question; decide explicitly about classes 1–2 | `ADOPT-NOW` — decided: questions carry no steps and state no answer, so classes 1–2 are solution-only, and that is now a comment rather than an omission |
| 4 | Fix `_DEGENERATE_PRODUCT` for the float-zero form | `ADOPT-NOW` |
| 5 | Line-oriented step probe | `ADOPT-NOW` |
| 6 | Phase 6 worklist of 14 doubled-sign templates; promote `_signed_term`/`_rect_str` to a shared emission helper before nine templates are fixed by hand | `ADOPT-PHASE-6` — added to Phase 6's deliverable list as **D6.7**; the census now prints the list on every run |
| 7 | Decide what an answer span may contain, and gate it | `ADOPT-NOW` for option (a), the template edit. Option (b) — a `**Note:**` terminator in `normalize.answer_span` — is **REJECTED**: it adds a rule to the *candidate*-side parser to accommodate a shape only *gold* produced, which is the wrong side of the contract and the "add a rule to make one template pass" move the phase brief warns against |
| 8 | Evidence is thinnest on question text; the marker predicate tests the marker, not the span | `ADOPT-NOW` for the question (F2/T8). The span-shape assertion is `ADOPT-PHASE-6` (**D6.8**): after F7 the corpus's longest span is 199 characters, so a bound can be set from measurement rather than guessed |
| 9 | Two plants of materially different surface form per detector, written from the class definition | `ADOPT-NOW` — `PLANTS` now carries 10 plants over 4 classes, including one in the **question**, and every pair differs in surface form: malformed bold vs bold lost; marker swapped vs marker added; `+ j-5` vs `+ -5`; `0*pi` vs `0.0*pi`. `SPEC-CHANGE 17` carries the rule forward |
| 10 | The gate says "clean across all four classes" and the corpus is not clean across class 3 | `SPEC-CHANGE 16` — the gate now says what it means |
| — | `stoichiometry.py`'s `SyntaxWarning: invalid escape sequence '\%'` | `ADOPT-NOW` — one character, in a file this phase already edits |

### What this changed about the phase's own claims

**The scan's `complex_sign 0 templates` line was a corpus claim that D-061
contradicted three paragraphs earlier**, and neither the implementer nor the
self-test caught it. The class is present on 14 templates; the detector could
see two of its shapes. Reporting `0` for that is worse than reporting `14`, and
the census now does the second.

**T8 is in `DEFAULT_CHECKS`, not only `ALL_CHECKS`.** The reviewer asked for
`ALL_CHECKS`. `ALL_CHECKS` is the opt-in list; `DEFAULT_CHECKS` is what runs when
someone types the command with no arguments, which is the only invocation that
happens by habit. A gate nobody types is the thing F4 is about.

---

## D-063 — N1's replacement is chosen by measurement, and the rule is PER-KIND

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track B, D5.6

`parse_number` read the **first** number in an answer span. The first number in
an answer sentence is very often a fluid grade (`Engine Oil (SAE 50)`), a
temperature (`at 541 K`) or a chemical-formula subscript (the 4 of C4H10).

**Eleven candidate rules, two axes, four corpora, held-out slice frozen before
any candidate was written.** `tests/comparators/n1_candidates.py` regenerates
the whole table, losers included; the split is a pure function of the template
id and a fixed salt, defined above the candidates in the file.

Held-out slice, the one that counts:

| rule | gold×gold FA | gold×gold decided | archive XI | archive decided |
|---|---:|---:|---:|---:|
| `first` (incumbent) | **16** | 100.0% | 144 | 99.1% |
| `last` | **424** | 100.0% | 1,812 | 99.1% |
| `unit_adjacent` | 0 | 70.6% | 1 | 70.0% |
| `unique_or_unresolved` | 0 | 65.2% | 0 | 3.6% |
| **shipped (per-kind)** | **0** | **100.0%** | 6 | **99.1%** |

### Three findings, none of which a single example would have produced

**1. "Last number" is far worse than "first", not better** — 424 held-out false
accepts against 16. The brief warned this; the measurement quantifies it.

**2. The tokeniser matters more than the choice of number.** A digit inside a
unit exponent (`m^2`) or a formula subscript (`C4H10`) is not a quantity at all.
Removing those two classes is orthogonal to *which* of the remaining numbers a
rule picks, so it applies under every rule equally — and it is what turns "last
number" from 424 false accepts into 0.

**3. The rule must be per-`kind`, and this is the reframing.** `check` answers
state the quantity first and the threshold second — *"the deflection is 18.4 mm,
less than the 25 mm limit"* — so every last-ward rule picks the **limit**. On
D4.4's 21 `check` cases:

| | decided | correct | false rejects |
|---|---:|---:|---:|
| first-number | 13 | **13** | 0 |
| every last-ward rule | 13 | **5** | 8 |

So `check` keeps first-number and `numeric` does not. **The comparator did not
need another rule; it needed the binding to say which rule applies.**

### Two confounds separated before the table meant anything

**A unit false accept is not an extraction error.** `400.0 MHz` against
`400.0 kHz` is the same number and a different unit; no extraction rule can fix
it and only a declared unit can (D-052). Counted together they would have
credited and blamed every rule for something outside its reach, so they are
separate columns.

**The incumbent's higher decided rate WAS the defect.** It decided more often
because it picked a non-exponent incidental number whose precision was
computable, while the real answer was frequently in exponent form and
`_decimals` returned `None` for those. It was deciding *by reading the wrong
number*. Implementing §7.1's display tolerance for scientific notation
(`extract.displayed_decimals`) lifted the archive decided rate under **every**
rule including the incumbent — which is what shows it to be an independent fix
and not a way of paying for this one.

---

## D-064 — N2's root cause was the isolation, and it is that function's THIRD wrong positional rule

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track B, D5.7

D-058 diagnosed N2 as *"an empty parse compares equal to an empty parse"* —
`_as_polynomial` returning `{}` and `{} == {}` being a `MATCH`. That is true and
it is not the root cause. Fixing only it left `autocorrelation_rect_pulse`
matching across instances, because the CAS then received `"0"` from **both**
sides and correctly agreed.

**The actual defect is upstream.** `_isolate_expression` took the **last**
`=`-bearing line, and that template's answer span ends with a sentence of prose:

```
R_g(tau) = 64*(8 - |tau|), for |tau| <= 8, and 0 otherwise.
This triangle has a peak value of 512 at tau = 0.
```

The last `=`-bearing line is the prose, so the isolated "expression" was the
string `0`, on every instance.

**This is the third wrong positional rule in one function, and the first two are
recorded as errors in `phase4_summary.md` §8** (#3: *"split on the **last** `=`
…, then on `.split("\n")[0]`. Two wrong rules from two examples, in one
function"*). So the fix is deliberately **not a different position**:

- an answer line **assigns to a symbol** and prose does not, so lines whose
  left-hand side is a bare identifier (optionally with arguments) are preferred;
- `=` is split on only when bare, never inside `<=`, `>=`, `!=`;
- and a **piecewise** answer — `…, for |tau| <= 8, and 0 otherwise` — is refused
  outright rather than having one branch picked, because a single-expression
  comparator cannot represent it. D4.1 §4.4: failure is `UNRESOLVED`, never
  `MATCH`.

`autocorrelation_rect_pulse` is now `UNRESOLVED` **in both directions**,
including for a pair of identical instances. That is the honest outcome and it
is why the template is named unbound rather than counted as fixed.

---

## D-065 — A binding that never decides is not a binding

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Track B, D5.8

D5.8's criterion as written is *"zero false accepts over N ≥ 50 instances
cross-paired within its own template"*. **Eleven bindings met it by producing
zero verdicts** — 8 of the 9 `symbolic` templates and 3 `multipart` ones decide
**0.0%** of their 2,450 pairs each.

That is the degenerate pass the phase brief names in its own D5.6 section —
*"the degenerate winner is candidate 3, which achieves zero false accepts by
deciding nothing"* — committed one deliverable later, on the axis where the
brief did not repeat the warning.

**Decision: a binding counts as bound only if it also decides.**

**The threshold is measured, not chosen.** The decided-rate distribution over
the 132 candidate bindings is bimodal with an empty middle:

| decided rate | bindings |
|---|---:|
| 0% | **11** |
| 0–25% | 0 |
| 25–50% | 0 |
| 50–90% | 1 (87.4%) |
| 90–100% | **120** |

Nothing lies between 0% and 87.4%, so **every threshold in that range gives the
same partition** and the result does not depend on where the line is drawn.

**This floor alone took the count from 132 to 121.** Reviewer E then found
four more reasons a binding can be wrong while passing an accept-only gate, and the
final figure is **118 of 150 bound, 32 named unbound** (D-067). Each successive
number is smaller and truer than the one before it, and that is Track B's version
of the Phase 4 stopping rule — *bind fewer templates and name the rest*, rather
than add rules.

**Why the eight `symbolic` templates cannot be rescued cheaply, measured.** The
first diagnosis was that `_to_sympy` splits function names (`cos` → `c*o*s`), so
`c`, `o`, `s` read as symbols outside the declared alphabet. Masking the
function names changed the decided rate by **nothing at all** — 11.1% before and
after. The real blockers are several and none is small:

| template | isolated expression | blocker |
|---|---|---|
| `phasor_addition` | `73.5 * cos(361*t - 146.39 deg)` | the unit word `deg` inside the expression |
| `bpsk_energy_basis` | `b) The basis function (psi_1(t)) is 16.24 * cos(...)` | the isolation takes prose |
| `standing_wave_formation` | `0.0, 0.69, and 1.38` | a two-equation answer; the isolation takes the node list |

Diagnosing from an error message and shipping the fix would have added a
mechanism that buys zero verdicts. **Assigned to Phase 6 (D6.10)** with these
three shapes named.

---

## D-066 — D-057's supporting-quantity remedy is declined by this phase, with a reason

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5, D-057 /
`phase4_summary.md` §12

`phase4_summary.md` §12 assigns Phase 5 the supporting-quantity remedy for the
two 100%-shortcuttable classification items. **Declined here and reassigned to
Phase 6**, and the disposition is recorded in §12 itself as the brief requires,
not only in this register.

**Reason.** The remedy changes what the item *asks* and what its gold *answers*.
That is item design. Track B's charter is `tests/comparators/` and
`template_inventory.csv` **only** — the spec says so in the same paragraph that
sets its effort box — and Track A's scope is the eleven output-contract
templates. Taking it here would be scope creep on the two items where a mistake
is least recoverable, in a phase whose named risk is *changing an item pool by
accident*.

**What this phase contributes instead, and it is not nothing.** Both items are
now bound and cross-paired at N=50 with **zero false accepts** — and both remain
**100% shortcuttable**:

| template | blind-guess floor | held-out surface model | gold×gold false accepts |
|---|---:|---:|---:|
| `system_property_linearity` | 0.5008 | **1.0000** | 0 |
| `system_properties_memory_causality` | 0.3450 | **1.0000** | 0 |

Those two facts are independent and the point is to hold them side by side:
**a comparator that scores an answer correctly cannot tell you whether the
answer required the reasoning.** A bound template is not a hard one, and the
binding work must not be read as evidence about difficulty.

One further measurement worth recording: `levenspiel_plot_interpretation`'s
blind-guess floor is **1.0000** — a blind guess scores 100%, so the statistic is
degenerate on that item and its difficulty label is unsupported by it.

---
## D-067 — Reviewer E's triage: the two defect classes were masking each other

**Date:** 2026-09-07 · **Status:** DECIDED · **Source:** Phase 5 Reviewer E
(comparator adversary), `reviews/phase5_reviewer_e_comparator.md`

**Verdict `BLOCKED`, and it was the right verdict.**

> **75 of the 132 bound templates do not return `MATCH` when the candidate is a
> verbatim copy of gold.** `gold_gold`'s loop skips `a == b`, so this case was
> not merely unmeasured — it was **excluded by construction**, and no gate saw
> it.

```
Counter({'MATCH': 57, 'UNRESOLVED': 45, 'MISMATCH': 30})
```

**And the load-bearing point, which is worth more than any single finding:**
the reported *"zero false accepts over 340,550 pairs"* was **produced by** an
over-rejection defect. E found a real over-acceptance mechanism — a part
comparing a different part's number — and showed that both of its realised
instances return `MISMATCH` **because of the unit defect**, not because the
comparator noticed anything. Repairing one defect uncovers the other. An
accept-only reading of that gate was measuring the interaction of two bugs.

This is the Phase 4 lesson arriving on a new surface: *a layer that marks
correct answers wrong is as unacceptable as a false accept* (D4.1 §1), and it
hides false accepts while it does it.

### Findings

| # | Finding | Disposition |
|---|---|---|
| **E-1** | 75 of 132 bindings reject a verbatim copy of gold; the identity case is excluded by the loop's own `a == b` skip | `ADOPT-NOW` — identity is a first-class term in **both** the binding gate and `cross_pair`'s gate |
| **E-2** | `_resolve_unit` knows 11 canonical units; **40 of the 48 declared values are outside it**, covering 72 of 103 templates — and it matches substrings, so `N` is found inside `N/m` → `MISMATCH` | `ADOPT-NOW` — the derived unit is no longer passed to the comparator |
| **E-3** | `DECLARED_UNITS` is a trailing-token heuristic, not units: `otherwise`, `e-05`, `units`, `percent`, `dollars`, `subgroups`. It re-creates exactly the rejection D4.1 §4.1 says the opt-in design exists to prevent | `ADOPT-NOW` — same fix; the **census** stays, the comparator does not consume it |
| **E-4** | `numeric_parts(n, unit)` applies one unit to all `n` parts of an answer whose parts carry different units by construction; 13 templates MISMATCH themselves | `ADOPT-NOW` — parts carry no unit; per-part units need a per-part declaration (`ADOPT-PHASE-6`, D6.11) |
| **E-5** | The unit check costs **19 of the 82** matches on real archived model answers, and catches 2 | `ADOPT-NOW` — this is the ablation the brief prescribes, and it says delete rather than patch |
| **E-6** | 11 templates are "bound" while deciding **nothing** — 2,450/2,450 UNRESOLVED | `ADOPT-NOW` — a decided-rate floor (D-065); found independently before the report arrived, which does not make it less E's finding |
| **E-7** | The shipped tool **already printed** `false rejects 56` and it never reached the claim table | `ADOPT-NOW` — false rejects gate now, and the count is a headline |
| **E-8** | `nth_quantity`'s 16-character tail makes part *i* select part *i+1*'s number on **6 of 31** multipart bindings; two real gold pairs are accepted for each other | `ADOPT-NOW` — the slice ends at the number |
| **E-9** | Three multipart bindings contain parts that read a **constant**, so those parts compare nothing | `ADOPT-NOW` — a constant-part check in the binding gate |
| **E-10** | `partial` cannot fire on a numeric or symbolic binding, and two symbolic bindings are narrower than their answers | `ADOPT-NOW` for the flag's scope; the two symbolic ones are unbound anyway under D-065 |

### The three claims that did not survive

**E8 diverged on all three of its parts and the divergence is instructive.**
The commit claimed N1/N2/N4 were fixed and demonstrated it with
`compare_kind(...)`, which uses *kind defaults*. Under the shipped **bindings**:

- `hagen_poiseuille_flowrate` is **unbound**, so `compare_template` raises
  `KeyError` rather than returning `MISMATCH`. The claim was true of the kind
  and not of the template.
- `autocorrelation_rect_pulse` stops matching partly because the pseudo-unit
  `otherwise` suppressed it — a defect doing the work of a fix.
- `continuous_to_discrete_conversion` traded crashing for refusing, which is
  progress and is not the same as being fixed.

**A demonstration must run through the shipped path.** Demonstrating a fix with
the generic entry point while shipping a bound one is the same class of error as
verifying a claim by re-reading it.

**E also noted that `phase4_comparators.md` §8 at that ref still prints
98.6%/100%** where the current measurement is 95.8%/98.6% — a Phase 4 document
inconsistency, `ADOPT-PHASE-6`.

---

## D-068 — Reviewer E round 2: the block lifts, and two of my own numbers were wrong

**Date:** 2026-09-08 · **Status:** DECIDED · **Source:** Phase 5 Reviewer E,
round 2, `reviews/phase5_reviewer_e_comparator.md`

**`PASS WITH FINDINGS` — the block lifts.** E answered its own block question on
its own terms rather than accepting mine:

> Neither colliding pair is a false accept, and both `MISMATCH` **for the right
> reason**. `hydrostatic_pressure_at_depth` part 1 now selects `13.29` — the
> depth — instead of `101.347`, the pressure it had been comparing twice;
> `rackett_equation_volume` part 1 selects the temperature instead of the
> volume. *The masking is gone and nothing was traded for it.*

And it closed the gap it had declined to claim in its first filing: gold×gold at
**N = 100 over all 33 structural-kind templates — 326,700 pairs, 0 false
accepts, 0 false rejects, 0 errors, 0 UNRESOLVED**, with round 1's N = 250
collision search re-run byte-identical and returning **0** collisions on all
five formerly mis-slicing templates. Identity, re-derived independently rather
than read from my tool: **5,900 / 5,900**.

### Findings

| # | Finding | Disposition |
|---|---|---|
| **R2-F1** | `vector_components` reads ASCII `x_hat`, which is what **gold** writes; every archived model answer writes `x̂`. The comparator accepts **0 of 83** real answers across the 3 vector templates — and this is **structurally invisible to the identity gate**, because gold matches gold. | `ADOPT-NOW` |
| **R2-F2** | The gold×gold truth predicate was **textual identity of the span**, so `5.0 seconds` against `4.99 seconds` read as a false accept when §7.1's boundary-inclusive tolerance makes the MATCH correct. **10 of the 32 unbound entries rested on that predicate, 9 on it alone.** | `ADOPT-NOW` |
| **R2-F3** | For `symbolic`, "bind fewer and name the rest" *is* laundering: `derive_kind` routes to `symbolic` on a regex matching exactly the construct the comparator cannot read. "1 of 9" measures the comparator, not the corpus. | `ADOPT-NOW` (the characterisation) + `ADOPT-PHASE-6` (the fix, D6.10) |

### R2-F1 is the sharper of the two mechanical findings

A missing **surface declaration** — the mechanism D4.1 §4.2 already uses for
categorical labels — and not undecidability. `x̂` is a combining circumflex
(U+0302) or a precomposed `ŷ`/`ẑ` (U+0177, U+1E91); models also write the
component list as `[a, b]` where gold writes `<a, b>`. Both are now folded.

**And it names a real limitation of the identity gate, one round after that gate
was added to catch exactly this class.** Identity compares gold to gold, so a
defect that lives entirely in *model* surface forms passes it by construction.
The identity term is necessary and it is not sufficient; archive×gold is what
covers the other half, and its per-kind decided rate is what would have shown
this — `vector` at **12.0%** was in the table and I did not read it. That is the
same error as E-7 one round earlier: **the number was printed and not read.**

### R2-F3, and a case where E and I were each right about different things

E diagnosed the `symbolic` failures as the parser splitting function names into
free symbols — `['c','o','s']` is `cos`. I had already implemented that fix and
measured **no change at all** (11.1% decided, before and after), and recorded it
in D-065 as a known-bad fix.

Both are correct, and resolving it required a third measurement. My first mask
was itself broken — `_to_sympy`'s `x*(` rule re-inserted a `*` inside the
placeholder — so *that* measurement was worthless. Repaired, the mask works on
the probe (`cos(361*t - 146.39)` survives intact) and **still rescues 0 of 8**,
because function-name splitting is one of **five** independent blockers:

| template | isolated expression | blocker |
|---|---|---|
| `phasor_addition` | `38.93 * cos(420*t + 7.48 deg)` | the unit word `deg` → symbols `d`,`e`,`g` |
| `ft_esd_rect_pulse` | `2304.0 * sinc^2(4.0*f)` | `sinc^2(` — a function with an exponent, which the mask's `name(` pattern misses |
| `cd_dc_system_analysis` | `27.35 * cos(657*pi*(t - 3.81e-03))` | `pi` → `p*i`, and `3.81e-03` → `3.81*e-03` |
| `undamped_response_initial_conditions` | `-0.006*cos(34.6384*t) (m)` | a trailing unit annotation → symbol `m` |
| `standing_wave_formation` | `0.0, 0.99, and 1.98` | the isolation takes the node list, not the equation |

So E's characterisation stands — the routing regex selects for what the parser
cannot read — and D-065's "the obvious fix is known not to work" also stands,
with the correction that it is one of five and that my first attempt to measure
it was broken. **Recorded as an error (§7 #15): I measured a fix, saw no change,
and concluded the diagnosis was wrong, when the fix was.**

Not fixed here. Five surface rules added at the close of a phase, to a parser
whose *last* three positional rules were all wrong (D-064), is exactly the move
the stopping rule forbids. **Phase 6 (D6.10) now has a specification instead of
a symptom**, and one ruled-out fix.

---

---


## D-069 — One provenance vocabulary: seven classes, a local-only qualifier, and a stated relation

**Date:** 2026-09-11 · **Status:** DECIDED · **Source:** Phase C1.2; applies the
SPEC-CHANGE C2 raised and did not apply (D-035, Reviewer G §7) ·
**SPEC-CHANGE 20, 21**

Two vocabularies were live. C2's four classes (`[ON-DISK]`, `[DERIVED]`,
`[BY-DEFINITION]`, `[KNOWN-DEFECTIVE]`) were off-spec; civil and industrial had
grown `[VERIFY: X]`, `[POLICY: sampling-only]`, `[REALISM]`, `[DERIVABLE]` and
`[UNVERIFIED]`, with `[ON-DISK]` meaning three different things across them - a
CAS in a JSON, a page in a book on one machine, and a page image read by eye.

**Decision.** Seven classes with one meaning each, in spec §C1.2:
`[ON-DISK]`, `[ON-DISK:LOCAL-ONLY]`, `[DERIVED]`, `[BY-DEFINITION]`,
`[POLICY: sampling-only]`, `[KNOWN-DEFECTIVE]`, `[UNVERIFIED]`. Three choices in
it are worth recording:

1. **Local-only is a qualifier on the class, not a comment.** Two citations to the
   same page of Das are equally honest and only one can be checked from a clone;
   resolvability is a property of the evidence, so it lives in the tag the
   resolver reads. The resolver counts a local-only citation it cannot open as
   `UNRESOLVABLE-FROM-CLONE` - never passes it, never fails it.
2. **`[ON-DISK]` states its relation.** `precision=` means the constant *is* the
   artefact's value (a transcription, checked exactly); `tol=` means it *agrees
   with* an independent artefact (a verification). This is Reviewer G's G-7 -
   *"the record documents verification, not origin"* - made part of the tag, so
   the two can no longer read alike. A citation stating neither is counted
   `LOCATOR-ONLY`.
3. **Deprecated tags are counted, not failed.** `[VERIFY]`, `[REALISM]`,
   `[DERIVABLE]` and the unlocated `[ON-DISK]` forms are `LEGACY` until C3.7
   retags them. At C1 the resolver lists **47** (civil 22, industrial 25). C3's
   gate is zero.

**The citation gate** (SPEC-CHANGE 21) reads "a retrievable on-disk artefact and a
locator within it" in C2's and C3's exit gates - C2 had in fact been gated on that
form, and C3 was about to inherit the unsatisfiable "edition + page".

**What enforces it.** `test_citations_resolve.py` keeps C2's P1-P4 and adds R1-R7,
one per clause; `--selftest` plants ten failure defects and one attachment
fixture, each judged by the failure it *adds* over a clean fixture.

## D-070 — The given-values rule excuses correctness, not truth; the classification is measured

**Date:** 2026-09-11 · **Status:** DECIDED · **Source:** Phase C1.3 · **SPEC-CHANGE 20**

The spec's C1.3 classifies each table as *needs a citation* or *plausibility
window only* "applying the given-values rule". Applied literally, the rule excuses
every value a question restates - and `phaseC1_census.md` measures that
`MATERIAL_DENSITIES`, `FLUID_DENSITIES`, `CRITICAL_PROPERTIES` and 20 other
property tables are restated in every consumer. They would all have been
"sampling policy". But C2 §6 found 274 items that stated a heat capacity ten times
too large and were then self-consistent about it: **restating a value makes an
item self-consistent, not true.**

**Decision.** The class is computed from two things, never one:

- a declared **`@kind`** per table (the judgement, written where it can be
  challenged), and
- a **measured** given-values status, by perturbation: nudge one field, re-run
  every consumer on the same seed, and see whether the question changed or only
  the answer did (`census.py` P-GIVEN).

Only a `range` - a window a value is drawn from - can be a plausibility window,
and only when every consumer restates it or reads it as a guard. At C1 the 107
numeric tables classify as **CITATION 36 · PLAUSIBILITY 29 · DERIVATION 4 ·
DEFINITION 5 · DOMAIN 1 · UNCONSUMED 32**, regenerated by
`python -m tests.constants_integrity.census --seeds 40 --markdown …`.

**"43 tables" was a definition, not a count.** Under the published predicate the
three original branches hold 38 and all five hold 107.

**A correction inside this decision.** The first draft of `@kind` drew the line
between `property` and `range` as "real material vs. abstract" and declared four
civil soil-class windows `property`. R7 - a tag must agree with its `@kind` -
flagged the contradiction with civil's own `[POLICY]` tag on `CV_RANGES_M2_YR`. The
right line is **single value vs. window**: mercury's density at 20 °C is a fact
every sample shares; a clay's coefficient of consolidation is drawn for a
hypothetical specimen and asserts nothing universal. `SPECIFIC_GRAVITY_RANGES`,
`PERMEABILITY_RANGES_CM_S`, `FRICTION_ANGLE_RANGES_DEG` and `CV_RANGES_M2_YR` are
`range`; three moved from CITATION to PLAUSIBILITY. Recorded because a
classification rule that moved three tables on the implementer's re-reading is
exactly what Reviewer G is asked to test.

**Limits the census states about itself**, in its docstring: a coarse display
rounding reads HIDDEN (the conservative direction); guard consumption reads ERROR
or NO-EFFECT; `REACTIONS` restates its coefficients through a parallel `equation`
string the probe does not nudge; a constants-level copy (`MEDIA_VELOCITIES
["Vacuum"]` from `C0`) is seen statically and not dynamically.

## D-071 — C1.7 pilot: `C0` is a transcription at `precision=exact`; `EPSILON_0` is a rounding, and a rounding is a transcription with a stated precision

**Date:** 2026-09-11 · **Status:** DECIDED · **Source:** Phase C1.7

Measured against `codata_2022/allascii.txt`: `C0 = 299792458` equals CODATA's
`299 792 458 (exact)`; `EPSILON_0 = 8.854e-12` is CODATA's `8.854 187 8188 e-12`
rounded to four significant figures, 2.12e-5 relative below it.

**`C0` → `[ON-DISK] … precision=exact`, not `[BY-DEFINITION]`.** The 2019 SI does
fix c by definition, and `@kind: defined` says so. But `[BY-DEFINITION]` is for a
value no artefact can be the warrant for (an element's ΔHf° = 0). A defined
*number* can be mistyped, and an artefact catches that; a definition cannot.

**`EPSILON_0` → `[ON-DISK] … precision=4sf`, not `[DERIVED]`.** The brief left this
open. A rounding is fully specified by its precision, which the resolver checks
exactly - `round(CODATA, 4 s.f.) == 8.854e-12` - so nothing is computed that a
reader would need to re-derive. C2's heats of formation are the precedent: rounded
to 1 dp and tagged `[ON-DISK]`. `[DERIVED]` stays for fits, identities and
conversions.

**D-025, applied.** No consuming question states ε₀: the census measures it HIDDEN
in all four electrostatics templates, which print it only in the solution, at
`.4e` and `.3e` - both lossless for a four-significant-figure value. The cost of the
rounding is therefore on the solver's side, not in the trace: a solver using the
full CODATA value lands 2.1e-5 relative from gold. None of the four declares a
`TOLERANCE`; D-025's representability assertion (C3.6) is the place to hold it.

## D-072 — The acquisition's record was wrong in two ways `--verify` could not see

**Date:** 2026-09-11 · **Status:** DECIDED · **Source:** Phase C1.1

Both commits are `9b45baf` and `fc1e501`.

1. **The manifest recorded a request that did not produce its file.** A present
   file was re-recorded under the first-choice URL its caller rebuilds. For the
   water and cyclohexane isobars that is the request NIST clamps to the triple
   point - the response the acquisition rejected - while the files start at
   293.15 K. `--verify` passed because it compared hashes, and a hash says the
   file is unchanged, not that the record describing it is true. Fixed; the URL
   is now re-derived from content, and `--verify` fails an isobar whose recorded
   URL cannot have produced its file (`--selftest`, four cases).
2. **A clone could not have verified at all.** git here runs with
   `core.autocrlf=true`, every committed source is LF-only, and `MANIFEST.json`
   pins bytes. A Windows checkout would have converted all 107 committed text
   sources (`MANIFEST.json`: 120 present, 13 binary) to CRLF and `--verify` would
   have failed every one. The first draft of this entry said 106, from memory;
   the manifest says 107. `.gitattributes`
   now marks `docs/references/**` `-text`. **Measured on a clean clone:** 13
   MISSING (exactly the gitignored binaries), 0 CHANGED, 0 URL-NOT-FILE.

And one crash: `os.replace` on the manifest raised `PermissionError` on the 78th of
~150 saves, a transient Windows file lock. Bounded retry. The handed-over
`MANIFEST.json` was overwritten by that run before it crashed and was never
committed, so whether it already carried the wrong URLs cannot be established.

## D-073 — Reviewer G (C1): four of my census predicates were wrong, and fixing one exposed eleven tables whose class nobody had checked

**Date:** 2026-09-12 · **Status:** DECIDED · **Source:** C1 Reviewer G,
`reviews/phaseC1_reviewer_g_provenance.md` (filed at `7fc5407`, committed
unmodified as `75df94d`) · **SPEC-CHANGE 22**

**Verdict PASS WITH FINDINGS.** G checked 14 PLAUSIBILITY tables from template
source and runs - all agreed - and 10 UNCONSUMED tables individually, with all 32
swept by name, dynamic access and literal copy. Three findings are CONFIRMED and
blocked the merge; all three are fixed, each with plants. **The fixes were not
re-reviewed by G**; C3's Reviewer G reads the same census.

### Findings

| # | Finding | Disposition | Action |
|---|---|---|---|
| **G-1** | `SCS_IA_RATIO` is consumed through a copied literal `0.2` (and `0.8`, `0.4`, `0.04` derived from it) that name-based detection cannot see; the census called it UNCONSUMED; a correction to the table would never reach the item | `ADOPT-NOW` + `ADOPT-PHASE-C3` | P-COPY: a copy is declared in the table header (`@copied-in`) and the census VERIFIES the literal is in that template's source, counting it as a consumer - COPIED → CITATION. Making the template read the table is C3.8, a P6 event with a before/after dump; it is not byte-identical by construction (`0.2**2 != 0.04` in binary) |
| **G-2** | An all-crash probe rolled up as NO-EFFECT, and `classify` then granted PLAUSIBILITY "as a guard" (`SERVICE_LEVELS`) | `ADOPT-NOW` + `SPEC-CHANGE` | ERROR is its own measurement; a `range` measured ERROR or INDETERMINATE is REVIEW until its header declares `@given: stated\|guard (evidence)`; `census --check` fails on REVIEW. **The first fix was incomplete one level down** - below |
| **G-3** | A graded chart-selection rule (Montgomery's "n > 10 or 12") rests only on two UNCONSUMED tables; the question never states it | `PLAUSIBLE` → `ADOPT-PHASE-C3` | Written response: UNCONSUMED is right by the letter, and the rule is still an unsourced warrant for a hidden, answer-deciding choice. C3.9's inline-window register names it; C3 decides whether the template states the rule or reads and cites the table |
| **G-4** | `PHASE_RANGE_RAD = (-math.pi, math.pi)` is read by a template but has no numeric literal, so it sat outside the census and the metadata check | `ADOPT-NOW` | P-TABLE-LIVE: a table whose *value* holds a number is classified and held to `@kind`/`@units`; the coverage predicate is unchanged, so the brief's figures still reproduce. Measured: the only such table in five branches |

### §5 suggestions

| Suggestion | Disposition | Reason / action |
|---|---|---|
| A mechanical literal-copy sweep over every table leaf | `ADOPT-PHASE-C3` (C3.8) | G-1 came from doing this by hand for ~12 values; it also catches CITATION values drifting into inline copies |
| An inline-window register for named-entity facts that never reach a table | `ADOPT-PHASE-C3` (C3.9) | Includes G-3 and `two_phase_specific_volume` |
| `SPECIFIC_GRAVITY_RANGES`' sand row is nearly a single mineral value | `ADOPT-PHASE-C3` (C3.2) | Stays `range` - stated in all five consumers - and C3.2 checks the window against Das |
| No kind for a one-sided screening threshold | `SPEC-CHANGE` (22) | `range` includes a one-sided screening bound on drawn values |
| Define "guard" statically rather than from a probe that cannot tell a guard from a crash | `SPEC-CHANGE` (22) + `ADOPT-PHASE-C3` (C3.10) | C1 requires a declared, evidenced verdict; the static detector is C3's |
| Consumer counts look inflated (`QUEUE_SCENARIOS` 4 vs 2) | `ADOPT-NOW` | **Confirmed, and the class was wider than G saw.** P-CONSUMER was flow-insensitive: a helper name bound from a table in one module-level loop and rebound from unrelated data in the next "stood for" the table in both. Rewritten flow-sensitively with a name-reuse plant; `QUEUE_SCENARIOS` now has 2 consumers |
| Look closer at `continuous_to_discrete_conversion`'s phase | `ADOPT-NOW` | `given_evidence.py` measures the degree phase printed on 24/24 of the seeds that draw it; `PHASE_RANGE_RAD` measures ALL-RESTATED |
| The optional part was not attempted | - | Nothing to triage |

### The first fix for G-2 was incomplete, and a plant found it

The first repair added an ERROR rollup. The census self-test's new `WINDOW2` plant -
a window restated on some seeds and rewording the question on others - still read
ALL-RESTATED. `field_verdict` returned the FIRST of HIDDEN, RESTATED, STRUCTURAL,
ERROR present, so a field restated on one seed and crashing on another reported
RESTATED, and G-2's escape stayed open. Precedence now puts what the probe cannot
settle above RESTATED.

**That exposed eleven `range` tables classified PLAUSIBILITY at `7fc5407` on a
RESTATED that masked a crash or a rewording.** Each now declares `@given: stated`
with committed evidence:

- five cite G's committed scripts, which measured exactly this (G §2 rows 1,
  10-13): `SPECIFIC_GRAVITY_RANGES`, `QUEUE_SCENARIOS`, `NEWSVENDOR_ITEMS`,
  `COMPONENT_RELIABILITY_CLASSES`, `SPC_CHARACTERISTICS`;
- six cite `tests/constants_integrity/given_evidence.py`, which captures each
  consumer's locals at return and requires every drawn value printed in the
  question: `FREQUENCY_RANGE_HZ`, `PHASE_RANGE_DEG`, `DECIMATION_FACTOR_M_RANGE`,
  `OMEGA_DENOMINATOR_RANGE`, `HOLDING_RATE_PER_YR`, `INVENTORY_ITEMS` - 19 checks,
  all on every seed that draws the value. `OMEGA_DENOMINATOR_RANGE`'s draw is not
  printed itself; the question states the reduced fraction it produces, exactly.

Fifteen `@given` declarations in all (these eleven and the four G's report settled
directly); the ratchet fails one that cites no committed evidence.

### And my rewrite broke something G had not

The flow-sensitive rewrite replaced a table's own entry with the tables it is built
from, so `RESISTOR_SERIES_BY_TOLERANCE` (built from the IEC lists) lost its consumer
and dropped to UNCONSUMED, and `MEDIA_VELOCITIES` (built from `C0`) lost its
order-dependent draw; and it reported the scalar `SHEWHART_K_SIGMA` as drawn by
order. All three were found by diffing every table's class, consumers and draws
against the pre-review census - not by a plant, because there was none for a table
that is both a table and an alias. There are two now.

### Result

**108 tables** (107 P-TABLE + 1 P-TABLE-LIVE): **CITATION 37 · PLAUSIBILITY 30 ·
DERIVATION 4 · DEFINITION 5 · DOMAIN 1 · UNCONSUMED 31**, `census --check` clean.
Against `7fc5407` the only class moves are `SCS_IA_RATIO` (UNCONSUMED → CITATION)
and the new `PHASE_RANGE_RAD`. The eleven were right at `7fc5407` by accident: the
measurement under them could not have said otherwise.

### A correction to D-070

D-070 says the literal given-values rule would have excused "`MATERIAL_DENSITIES`,
`FLUID_DENSITIES`, `CRITICAL_PROPERTIES` and 20 other property tables" - 23, a
number I did not count. **Counted from the census:** at `7fc5407`, **25** tables of
kind `property`, `standard` or `measured-constant` measured ALL-RESTATED, **15** of
them `property`; under the corrected precedence, **16** (10 `property`), the other
nine restating on most seeds and rewording the question on some. The argument
stands; the figure was wrong. DECISIONS is append-only, so it is corrected here.

## D-074 — `precision=`, `tol=` and `[KNOWN-DEFECTIVE]`: when a disagreement is a tolerance and when it is a defect

**Date:** 2026-09-12 · **Status:** DECIDED · **Source:** C3.1/C3.7 · **SPEC-CHANGE 23**

C1.2 gave `[ON-DISK]` two relations and **no rule for choosing between them**. That is an
open door: any constant that disagrees with its artefact can be made to pass by widening
`tol=`, and the tag still reads as evidence. Electrical's Helium went in through it - first
tagged `tol=0.0002%` (9380721), a tolerance wide enough to swallow a value that is **+3.2%
on the refractivity** it actually determines.

**The rule, now normative:**

- **`precision=`** when the literal is a rounding of the artefact **at the table's own
  declared conditions**. This is the default and the only one that says *this value came
  from here*.
- **`tol=`** only when the conditions differ, when none are stated, or when the value
  crosses a non-decimal unit conversion - and then `tol` is **half a unit in the artefact's
  last printed digit**, derived, never a number chosen so the row passes.
- **`[KNOWN-DEFECTIVE]`** when the conditions match and the literal is **not** a rounding.

**What it cost to apply, which is the point.** Helium became `[KNOWN-DEFECTIVE]`. So did
three civil water constants that had read as ordinary: `UNIT_WEIGHT_WATER_KN_M3` (+0.21%),
`_PCF` (+0.13%) and `WATER_KINEMATIC_VISCOSITY_M2_S` (+0.060%) - each small enough to hide
under a tolerance and each, measured, describing water at a different temperature from the
one its own comment claims. A rule that never reclassifies anything is not a rule.

## D-075 — Six locator forms and three relations the spec did not have; and the `xlsx` form it did specify was never built

**Date:** 2026-09-12 · **Status:** DECIDED · **Source:** C1.7, C3.7 · **SPEC-CHANGE 23**

C1.2's locator table specified `xlsx` as `sheet="<name>" row="<key>" col="<header>"` - one
tag per leaf. `AISC_W_SHAPES` has **196 leaves** (14 shapes x 7 US and 7 SI fields), so that
form buries the table it documents. Built instead: a **whole-table** relation,
`sheet= label_col= rows=all blocks="us:1:key,si:2:si_label" precision=exact`, where `blocks`
names each row sub-dict's header block and where its label comes from; every leaf is
compared with its cell, a missing header or label is R3, a disagreeing leaf is R4.

Also added, none of them in the spec: `member=` + `wavelength=` (a refractiveindex.info
dataset evaluated at a wavelength, refusing one outside its range), `page=` + `mil=`
[+ `col=`] (a MIL-HDBK-5J design table, the caption pinning the table and `mil=` the row),
and the relations `field=[key]`, `via="NAME/x"` and `scale=<f>`.

**`image-only` was specified and never implemented.** A bare `page=` already resolves and
counts `LOCATOR-ONLY`, which says the same thing without a second vocabulary. MIL-STD-105E
is cited that way deliberately: its OCR text layer extracts as `":: .... m n mO ::,i:,o Vtm
Loi o, batcb ah:e ."`, and **a text anchor that cannot be trusted is worse than none**.

## D-076 — MIL-STD-105E was a public on-disk artefact, cited for two phases as if it were a copyrighted book

**Date:** 2026-09-12 · **Status:** DECIDED · **Source:** C3.7 industrial

Its header cited `pilot/references/public/mil_std_105e_sampling.pdf` - the gitignored
copyrighted-books tree - and its three tables carried `[ON-DISK: PDF p. 18 (document
p. 13)]`, a page number in prose. `MANIFEST.json` vouches for the same document at
`docs/references/industrial/mil_std_105e_sampling.pdf`, hash and all. All three now cite it
directly and resolve.

**The general lesson is about what a checker cannot see.** The resolver validates a path
only once a tag is in `artefact @ locator` form; a LEGACY bracket with the path in prose is
counted but never resolved. So a document can sit in `docs/references/`, vouched for and
hashed, while the tables that depend on it claim to be unresolvable from a clone. Every
`[ON-DISK: …]` prose variant was a place this could hide, which is why C3.7's exit gate is
**zero LEGACY tags**, not "fewer".

## D-077 — A near-derivation is not a derivation

**Date:** 2026-09-12 · **Status:** DECIDED · **Source:** C3.7 chemical and civil

Three tables came within ~1% of a clean derivation and are **not** tagged `[DERIVED]`:

- `SUBSTANCES_FOR_VAPORIZATION`: dHvap = Enthalpy(v) - Enthalpy(l) at 1 atm from the NIST
  saturation tables reproduces 7 of the 12 unsourced rows to within 0.04-1.45%, but only
  methanol, benzene, ammonia and oxygen are **roundings** of it.
- `LIVE_LOADS_KPA`: 4 of 5 rows are the soft conversion 1 psf = 0.048 kPa exactly; the
  corridor row is neither that nor the SP 811 factor.
- `GRAVITY_FT_S2`: 32.2 is a rounding of g_n/0.3048 taken **after** the conversion, where
  `scale=` compares in the artefact's own unit - so it is `[DERIVED]` with the arithmetic
  stated, not an `[ON-DISK]` relation that would have to fail.

`[DERIVED]` asserts that *this* value follows from *those* inputs. A derivation that
reproduces most of a table is evidence about the table's provenance and is recorded as
such - in the C3.5 register, where the owner can act on it - but tagging it `[DERIVED]`
would make a claim that is false for the rows that miss, and those are exactly the rows a
reader would most want flagged.

---

# The evaluator pilot (`evaluator_pilot_17092026/`)

The workstream-01/04/05 pilot from `EngTrace_Suggested_Actions.pdf`: compare
evaluator candidates E0-E5 against expert labels on one frozen slice. Its record is
[`evaluator_pilot_17092026/`](../../evaluator_pilot_17092026/) — README (setup and
stage 1), FINDINGS (defects in the published framework), RESULTS_E0 / RESULTS_E3,
and `analysis/` (the script behind every reported number). The decisions below are
the ones a later reader would otherwise have to reconstruct from commit messages.

## D-078 — The pilot slice spans all five branches and is pinned by text

**Date:** 2026-09-17 · **Status:** DECIDED · **Source:** Suggested Actions §01

The Suggested Actions sketch stratified the slice across the three original
branches. It spans all five: civil and industrial were built to be in the next
revision, and a slice without them validates an evaluator on a corpus the paper
will not report. 15 cells (5 branches x 3 levels), one template per cell, four
instances each: 60 items, 300 traces at five models.

Items are pinned by **question and gold text plus SHA-256**, not by seed. A seed
does not pin an item across template versions — Phase 1 moved 34% of its
instances and C3 moved 824 questions — so the freeze records the text the
annotators will read, and `freeze.py --verify` fails the moment it moves.

## D-079 — Trace sources, routing, and three substitutions forced by availability

**Date:** 2026-09-17/18 · **Status:** DECIDED · **Source:** pilot stage 1

The five trace sources are the Suggested Actions' own (GPT-5, Claude Opus 4.7,
Gemini 3.1 Pro, DeepSeek R1, Llama 3.1 70B): three judge families, two controls.
Three things had to change, each recorded in `models.json`:

- **`OPENAI_API_KEY` and `ANTHROPIC_API_KEY` return 401.** GPT-5 and Claude go
  through OpenRouter. `openai/gpt-5` is priced there identically to the direct API.
- **`deepseek/deepseek-r1` has one provider with a hard 16,000-token output cap**,
  and R1 truncated mid-reasoning on 28 of 60 items. Moved to `deepseek-r1-0528`
  (same model, four providers, all >= 32,000) and re-ran all 60, not the 28:
  the survivors were a different checkpoint, skewed to the easy items.
- **E0's published Gemini judge, `gemini-3-pro-preview`, returns 404.** Replaced
  by Google's named successor. (Largely moot — see D-080, E0-F6.)

## D-080 — E0 is run unmodified; its defects are recorded, not fixed

**Date:** 2026-09-17 · **Status:** DECIDED · **Source:** pilot stage 2

E0 is the anchor every candidate is compared against, so it must be the framework
the paper describes. `evaluators/e0_tribunal.py` imports the published file and
calls its own `evaluate_entry`; the file's SHA-256 is in the config hash. Seven
defects were found in it (FINDINGS E0-F1 to E0-F7), among them that the Tribunal
is **two judges, not three** (`genai.get_model_info` does not exist and is
swallowed) and that a random 20% sample meant for the error taxonomy moves the
headline score. None is fixed in E0. `e0_3j` is the one counterfactual, and it
restores the third judge by supplying the missing function, not by editing code.

Deviations that could not be avoided are written on every scored row: D1 Gemini
judge substituted, D2 OpenRouter transport, D3 seeded sampling.

## D-081 — Scoring libraries are pinned, and their versions are config

**Date:** 2026-09-17 · **Status:** DECIDED · **Source:** FINDINGS E0-F3

Under transformers 5.x, BERTScore raises on every entry and the framework returns
0.0 — a column of zeros that looks like data. `requirements.txt` is unpinned, so a
fresh install reproduces it. The pilot venv pins transformers 4.57.3,
sentence-transformers 5.1.x and bert-score 0.3.13 (the line current at the
published run), and the versions are in the evaluator config hash so a cached score
is never served across library stacks. torch is excluded from the hash and
recorded per row instead, so GPU (Kaggle) and CPU results can share a cache —
after `--import-kaggle` proves they agree (159 traces: Tier 1 identical, BERTScore
within 1.8e-7).

## D-082 — Every judged evaluator samples the same wrong answers

**Date:** 2026-09-18 · **Status:** DECIDED · **Supersedes:** the original seed scheme

The harness first seeded E0's 20% wrong-answer draw from (evaluator, item, model),
so each candidate judged a different sample: e0 judged 175 traces, e0_3j 179, 160
in common, and column deltas looked like mechanism differences. The seed is now
(item, model) alone, and the scheme is hashed into the config so rows under the
two schemes never mix. E0 was re-run on it ($6.53).

## D-083 — Small open-weight roster models join as an unlabelled robustness cohort

**Date:** 2026-09-18 · **Status:** DECIDED · **Source:** roster planning

If the main benchmark moves to smaller open-weight models (the supervisor's budget
steer), an evaluator validated on the pilot's five must still function there.
Qwen3.8-27B and Gemma-4-31B-it were added as a `robustness` cohort over the same
frozen items. They are **not** in the labelled 300: the question is mechanical
(step markers, parseable answers, milestones reachable) and needs no expert labels,
while adding them would cost 40% more annotation. Result (FINDINGS R-F1, R-F2):
structure does not degrade; runaway reasoning does.

## D-084 — E3's milestones are derived per instance, by rule, from the template's own values

**Date:** 2026-09-18 · **Status:** DECIDED · **Implements:** D-001

Per D-001, milestones are emitted per instance. No template is edited: each frozen
item is regenerated at its seed with the repo's frame-local capture, and must
reproduce **byte-identically** (60 of 60 do) before a value is used. A milestone is
a value the template **computed**, the gold **states**, and the question does
**not** give. Iteration trajectories are excluded, because a correct trace from a
different starting guess cannot reach them.

Tolerance is 0.5%, the rule's own display tolerance. The first version used E0's
2% and failed its null baseline: traces "reached" 23% of a sibling item's
milestones. At 0.5% that coincidence rate is 4%, real coverage 0.749.

## D-085 — E4 checks the arithmetic a trace states, not symbolic equivalence to the gold formula

**Date:** 2026-09-19 · **Status:** DECIDED

The Suggested Actions phrase E4 as algebraic equivalence against the gold formula.
That needs aligning a trace's symbols with the gold's (`V_sat`, `V_{sat}`, `V_L`),
which free text does not support reliably. E4 instead checks what catches a right
number reached for the wrong reason: whether the arithmetic a trace **shows**
produces the number it **states**. The checker (`evaluators/arith.py`) was
validated on the 60 gold solutions before any trace — gold arithmetic is correct by
construction, so every inconsistency there is a checker bug. It started at 61.5%
consistent and reached 100% (232 of 232 claims) after ten parser fixes, each with a
plant; they are listed in the module docstring.

## D-086 — Compute goes to Kaggle; API keys never do

**Date:** 2026-09-17 · **Status:** DECIDED

The Tier 1 scorer stack runs on a Kaggle GPU (300 traces in 769 s on a T4, against
hours on the laptop CPU). The paid judge calls stay local: a CLI-pushed kernel
cannot attach Kaggle Secrets, and sending the keys to a third party is not a
decision to take implicitly. The bundle is staged only after the full freeze
rebuild and T1-T7 pass, and is scanned for every key value in `.env` before upload.

## D-087 — A checker is validated on real output, not only on gold

**Date:** 2026-09-19 · **Status:** DECIDED · **Extends:** D-085

At 100% consistency on the 60 gold solutions, E4's arithmetic checker still flagged
15 milestones in real traces as contradicted, and reading them showed 13 were
checker bugs. Gold is formatted uniformly by the templates; model traces are not
(glued variables like `16t`, Greek `μ`, `\mathrm{m}^2`, superscript runs, clause
fragments). Gold validation catches the checker's bugs on gold-shaped text only.

The rule, for E4 and for any later evaluator: validate on gold, **then read every
flag the checker raises on real traces before reporting a number it produces**,
and turn each genuine catch into a plant. Applied here it left 2 contradicted
milestones in 1,494, both genuine, and measured the unread side metric
(arithmetic consistency) at about two-thirds precision on a fixed-seed sample,
which is reported with that caveat rather than as a validated per-model measure.

## D-088 — Judges for E1 and E5

**Date:** 2026-09-19 · **Status:** DECIDED (E1 and E5 both run the same day with these
judges; RESULTS_E1, RESULTS_E5) · **Evidence:**
`evaluator_pilot_17092026/JUDGE_SELECTION.md`

Table 13's 27 models span seven families once backbones are counted (OpenAI,
Anthropic, Google incl. Gemma, DeepSeek, Meta incl. MetaMath, Qwen, Mistral incl.
WizardMath). Every capable remaining family has its own base but documented exposure
to an evaluated family's outputs, so independence is a degree, not a yes/no. A
21-step probe with known labels ruled out Nemotron 3 Ultra and Seed 2.1 (half their
replies unparseable) and showed no candidate rubber-stamps; it cannot separate the
top five.

Recommended: E1 panel **MiniMax M3 + Xiaomi MiMo-V2.5-Pro + xAI Grok 4.6** (three
families, residual exposure spread across Claude and OpenAI rather than stacked;
~$3.50). E5: **MiMo-V2.5-Pro** alone - revised from MiniMax M3 the same day: a
single judge carries its exposure alone, and MiMo matched MiniMax on the probe with
the smallest documented Claude exposure (0.4M vs 13M+ exchanges). Kimi K3 and GLM-5.3 deliberately excluded to keep
the strongest open families available for the next roster. The prior decision is the
roster's families: a family cannot be both judge and judged.

## D-089 — E1's judges get uniform call settings, and one provider is excluded

**Date:** 2026-09-19 · **Status:** DECIDED · **Source:** E1 smoke test and first replay

The framework gives each judge slot its own settings, tuned to E0's judges; the
Anthropic slot's 2,048-token cap would have truncated all three E1 judges (4, 9 and 16
of 21 probe replies over it), and a truncated reply is silently dropped. E1's judges
all get JSON mode, temperature 0 and 16,384 tokens (D6). Empty or malformed replies are
re-requested (D7). OpenRouter provider ModelRun is excluded for MiniMax M3 after 11 of
its 17 replies came back as malformed JSON against 0 of 161 from its other providers
(D8). All are written on every E1 row.

## D-090 — E2's PRMs, and a self-check on verdicts rather than exact rewards

**Date:** 2026-09-19 · **Status:** DECIDED (user's choice after the diagnostic) ·
**Evidence:** `evaluator_pilot_17092026/hpg/README.md`

E2 runs three open PRMs on HiPerGator (RTX PRO 6000 Blackwell):
- **Qwen2.5-Math-PRM-72B** is the primary.
- **VersaPRM** is the multi-domain PRM (a LoRA adapter on Llama-PRM800K).
- **Qwen2.5-Math-PRM-7B** is a size cross-check.

Every repo is pinned to a commit, VersaPRM's base included. Each job first scores its
card's example and refuses to go on unless the check passes.

The check was first designed as a 0.02 tolerance on every reward. On this hardware the
72B misses on the card example's two borderline steps (0.16 / 0.56 against 0.32 / 0.82)
while every step keeps the card's verdict at 0.5. A diagnostic ruled out each
explanation in turn:
- **The transformers version.** 4.46.3 and 4.57.3 gave bit-identical rewards.
- **Softmax precision.** Taking it in bf16, as the card does, changed nothing.
- **The pipeline.** The 7B lands within 0.033 of its own card with the same code.

What remains is a sensitivity to the attention kernel: sdpa vs eager moves the 72B by
about 0.05. flash-attn is not in the container, and no non-Blackwell GPU here can hold
the 72B.

The self-check now passes on per-step verdict agreement at 0.5. The continuous gap is
recorded in every run's metadata (`max_abs_diff`, `within_tol`). Qwen runs with eager
attention, which lands closest to the card. E2 is therefore reported primarily on
thresholded verdicts and on ranking of the probe's known-label steps, and its
continuous rewards are not compared with published numbers.

**Outcome (RESULTS_E2).** All 360 traces scored. On the probe's known-label steps the
72B separates slips from clean steps (AUROC 0.935), the 7B less well (0.824), and
**VersaPRM fails validation**: it passed 11 of 12 known arithmetic slips and rates
95-99% of every model's steps correct. E2's score is therefore the 72B's; VersaPRM is
reported as a negative result, not as a score.

## D-091 — Commit messages are one line, and the history was rewritten to match

**Date:** 2026-09-20 · **Status:** DECIDED (owner's standing preference) ·
**Evidence:** `docs/tasks/rewrite-commit-messages.md`, `tools/remap_commit_hashes.py`,
`tools/verify_commit_messages.py`

A commit message is its subject line and nothing else: no body, no
`Co-Authored-By` trailer. This is the convention for new commits, and the
existing history was rewritten so the whole log reads the same way.

`git filter-repo` rewrote all 17 local refs in one pass over 301 commits. The
message callback keeps the *subject* in git's sense — the first paragraph,
wrapped lines joined with a space — rather than the first physical line. Four
commits had a subject that wrapped onto a second line, and taking the first
line alone would have truncated them mid-sentence. Subjects of 100-135
characters already existed in this history, so a joined subject is in keeping
with it.

Nothing else moved: trees, authors, committers and author dates are
byte-identical through the commit map for all 301 commits, and the branch
topology is unchanged (`pilot/llm-annotation-pilot` keeps its 1 own commit,
`pilot/phase1-template-authoring` its 5, and the 14 `redesign/*` branches
remain ancestors of `master`). 37 commits kept their hash; the rest changed.

Changed hashes dangle wherever a file cites one, so `tools/remap_commit_hashes.py`
rewrote **199 citations of 76 distinct commits across 49 tracked files** to the
new hash of the same commit, at the same abbreviation length, as one follow-up
commit rather than by rewriting blobs in history. Two of those citations are
more than prose: `AUDIT_REV` in `tests/template_integrity/regen_inventory.py`
(now `67d41f4`), which still reproduces the Phase 0 calibration 65/37/89 with
+0 deltas, and the frozen slice's provenance in
`evaluator_pilot_17092026/slice/FREEZE.json` (now `d430246`). SHA-256 content
digests are untouched: the scanner matches hex runs of 7-40 characters with
hex-aware boundaries, so a 64-character digest can never match.

The old history is preserved on the branch
`backup/pre-commit-message-rewrite-20260920` (old `master` tip `2046e26`), on
`origin` as well as locally, so every pre-rewrite hash still resolves.
`origin/gh-pages` is separate history and was left alone.

## D-091 — The hard case is scored on the step labels, not on a second labelling round

**Date:** 2026-09-23 · **Status:** DECIDED · **Evidence:**
`evaluator_pilot_17092026/analysis/hard_case_pool.py`

A reasoning evaluator earns its keep on traces whose final answer is right but whose
reasoning is not. X1 could not run that comparison: of the 228 correct-answer traces
the experts' holistic verdict calls only 3 unsound. The plan was a second round —
generate traces from weaker models, mine them, have the experts label them.

The labels already hold the set. The experts' own step labels mark at least one step
incorrect in **56** of those 228 traces (102 steps), and call 53 of those same traces
sound overall. The trace-level target, not the data, was hiding the hard case. So the
hard case is scored as "does this trace contain an incorrect step, given the answer is
correct", and no new traces are generated or labelled.

The result on that target is negative and worth reporting as such: no evaluator beats
E0 (best is E2's 72B minimum reward, 0.601 against 0.545, CI of the difference
-0.041 to +0.166), and E3, E4 and E5 are significantly *worse* than E0 — they score a
trace by milestones the answer already implies. Every AUROC sits near chance.

A second round was also costed and rejected on its own terms. The signals that would
select candidates barely enrich: the 72B's minimum reward under 0.20 yields 28%
against a 25% base rate (1.16x), a failed arithmetic claim yields 33% at n=12. At
those rates ~170 new traces must be labelled to add ~50 hard cases, and the resulting
intervals (±0.06 instead of ±0.085) would not change the conclusion. 101 of the 102
incorrect steps are calculation slips, one is conceptual, so the population being
bought is mostly arithmetic noise that leaves the answer intact.

The candidate miner stays in the script for the record: 87 unlabelled robustness-cohort
traces (gemma-4-31b, qwen3.8-27b) pass a deterministic answer check, 36 of them fire at
least one signal. Re-run it if a later round is ever funded.

## D-092 — Layer 0: a deterministic gate admits templates to certification, sized to the defect rate

**Date:** 2026-09-23 · **Status:** DECIDED · **Source:** `template_annotation_23092026/layer0/`
(`gate.py`, `check_limits.json`, `closure_fixes.md`, `tie_census.py`)

The December 2025 certification pipeline was an AI Tribunal followed by human review,
and neither stage saw the defects the September audit found in the same 90
certified templates (chains that do not reproduce their own answers, unseeded
generators, doubled signs, wrong constants). The deterministic integrity suite
does see them, so it becomes **Layer 0**: no template reaches a judge or an
expert until it passes T1 closure, T3 determinism, T4 contract and T8 emission
with no generation error. T5 binding and T7 asserts stay advisory, as
`phase6_residual_register.md` already classifies them; T2 is reported where an
oracle exists.

**The gate runs at 500 seeds.** A 25-seed snapshot showed 29 templates failing
closure; 100 seeds exposed 6 more and 500 seeds 16 more, all the same
product-of-short-decimals tie class at rates of 0.1–3% per instance. That is
D-016's carried-forward point: acceptance evidence must be sized to the defect
rate, and a gate that is green at 25 seeds and red at 500 is not green.

**Resolution of the closure failures.** D-016 and D-037 applied unchanged: a
displayed-and-consumed operand is bound through its display; a half-way tie is
removed, by lengthening an exact display first and by a bounded redraw only where
the quantity is a quotient exact at no display; a tie in a final answer is redrawn
rather than lengthened (D-044). Every edit, its cause, and how much of the item
pool it moved is in `closure_fixes.md`; the corpus-wide movement is measured by
`c3_instance_dump.py --diff` against the pre-edit tree.

**Four lines are excused, not fixed**, in `check_limits.json`, each keyed to its
template and a line pattern so anything else in that template still fails: a
mixed-unit ratio rewritten across `=` (consolidation), hours-to-minutes
conversions across `=` (server selection), stated round-down and round-up lines
(Cpk), and the omega line of `wave_parameters_basic`, where T1 sizes a `%e`
result's tolerance from the mantissa alone — `core.printed_precision()` drops the
exponent — which is the closure counterpart of D-017 and is deferred to the
harness owner on the same reasoning.

**Finding, recorded for the harness owner:** T1 files an exact tie as FAIL or
MARGINAL by floating-point luck, so its failure count undercounts ties several-
fold. An exact-decimal census (`tie_census.py`, adapted from the civil round-2
agent's instrument) now reports residual ties per template beside the gate; it
gates nothing. Also for the harness owner: for negative exponents the same
`printed_precision()` defect makes T1 far too lenient (a `1.3200e-05` result gets
a tolerance of 5e-5, not 5e-10), so closure failures in small scientific-notation
results are invisible today.

**Sign-offs this creates (D-044, the owner's):** the answer display was
lengthened in `hydrostatic_pressure_at_depth` (kPa 3 → 4 dp) and
`max_hump_height_no_choking` (m 3 → 4 dp); it is recommended but NOT done in
`beam_internal_moment` (2 → 3 dp would replace a 10.5% redraw that halves every
x = k.5 section), `terzaghi_strip_footing_bearing` (1 → 2 dp would replace a
one-third depletion of c' = 15 general shear) and `effective_stress_profile`
(1 → 3 dp would replace a 1.3% redraw). Gold answers move on unchanged questions
in `absorbing_chain_time_to_failure` (8.3%, ±0.01 week) and
`effective_stress_profile` (24%, ±0.1 kPa), in both cases because a double
rounding was removed; the alternative in each is a tie screen with no movement.

---

## D-093 — Layer 1: the template screen keeps the published prompt, changes the judges, and is capped at two passes

**Date:** 2026-09-23 · **Status:** DECIDED · **Source:** `template_annotation_23092026/screen/`
(`run_screen.py`, `analyze_screen.py`); analysis of 2026-09-22 in the session record

The AI Tribunal is kept as a **reader**, not a gate, for three reasons: it is
the only exhaustive plausibility-and-prose read of all 150 templates and that is
exactly what it caught in December (rubber loaded to 365 kN, a 55 MPa pipe
pressure drop, a "titanium plastic container"); the July 2026 rebuttal promised
the aggregate sigma-max and pairwise Gwet's AC1 among the template judges, and
the December per-judge outputs are not in the repository (only the 91-row
`tribunal_summary.csv` survives, without per-judge flags), so the statistic needs
a re-run; and one method should cover all five branches.

**Judges:** the pilot's non-suite panel (D-088) — xAI Grok 4.6, MiniMax M3,
Xiaomi MiMo-V2.5-Pro — through OpenRouter with D-089's uniform settings (JSON
mode, temperature 0, 16,384 output tokens, three attempts on an empty or
malformed reply, provider ModelRun excluded for MiniMax). No judge shares a
family with an evaluated model, which removes Reviewer yAYU's judge-overlap
objection from certification; Kimi K3 and GLM-5.3 stay excluded so the strongest
open families remain available for the roster (D-088). The served model id,
serving provider, finish reason, tokens and OpenRouter's own reported cost are
written on every row, because the published run recorded none of them and its
Google judge has since been retired (E0-F4); the `.env` of the December run names
a Flash model where the paper says Gemini 3, which a recorded id would have
settled.

**The prompt is Appendix H's, verbatim,** read by AST out of
`ai_assisted_quality_assurance/run_ai_tribunal.py`, with three instances from
fixed seeds through the integrity suite's generator, so a pass is reproducible.

**Two passes, no more:** one before human certification and one after, so the
report can state the panel's view of the certified corpus. Templates are never
iterated against the judges — that made the Tribunal the oracle in December and
hid that the human stage rejected nothing. Layer 0 is the gate. `--pass 3` is
refused in code. Replies are committed (unlike December's), so the statistic
can be recomputed later.

**Cost:** about $4.60 per pass on measured reply sizes (each judge spends 1,000–
3,000 billed reasoning tokens behind a 150-token JSON row; a flat 200-token
assumption underestimated the pass threefold), so $9.20 for both passes. Spend
so far: about $0.04 (a three-judge probe and a three-row smoke test). **Pass 1
runs only on the user's explicit approval** (user rule of 2026-09-23: no paid
API work without prior approval).

---

## Open decisions

| # | Decision | Needed before |
|---|---|---|
| D-092 | ~~Answer-display lengthening and gold-movement sign-offs~~ **Decided by the owner 2026-09-23:** the two lengthened answers stay; `beam_internal_moment`, `terzaghi_strip_footing_bearing` and `effective_stress_profile` answers are lengthened too; the gold movements are accepted; and a scoped third round removes the census ties at 3% and above (eight templates), with the frozen pool to be censused before inference and stragglers fixed then | — |
| D-093 | Approval to run screening pass 1 (~$4.60) once the Layer 0 branch is merged | Layer 2 (human certification) starts |
| — | Which families the next roster will evaluate (Kimi, GLM stay available as long as they are not judges) | the next benchmark run |
| — | Expert annotation of the frozen 300 (stage 3): annotators, protocol, the ~100-trace triple-labelled overlap | any X1 agreement number |
| D-003 | Do the raw `inference_results/` generations still exist? | promising any corrected results table |
| — | Phase 5 scoping: fold into Phase 1 or run as a parallel PR | Phase 1 start |
| — | Whether the 9 self-inconsistent templates are fixed or replaced | Phase 1 start (item-pool ownership call) |
