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

*(2026-09-27, D-112: a second probe step first labelled clean carried a real slip. Re-scored,
Kimi K3 and Grok 4.6 lead the probe at 0.96 and MiMo and MiniMax tie at 0.92, so "matched
MiniMax" still holds and the probe still cannot separate the top. The choice stands: see
D-112 and JUDGE_SELECTION, "Does the relabel change the choice?")*

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
reported as a negative result, not as a score. *(Re-scored 2026-09-27 after the D-112
relabel, 13 slips and 8 clean steps: the 72B 0.923, the 7B 0.779, VersaPRM 0.548, passing 12
of 13 slips. The conclusions do not change.)*

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

## D-096 — The hard case is scored on the step labels, not on a second labelling round

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
before the pass: about $0.04 (a three-judge probe and a three-row smoke test).
**Pass 1 runs only on the user's explicit approval** (user rule of 2026-09-23: no
paid API work without prior approval).

**Pass 1 ran 2026-09-23 on the user's approval**, judge by judge under a spend
cap (`--judge`, `--max-usd`, added for the purpose): 450 of 450 rows parsed on the
first or second attempt, no failures, $4.15 in OpenRouter-reported cost (Grok
$3.29, MiniMax $0.39, MiMo $0.47). Outcome under the paper's rule: 126 pass, 13
controversial, 11 critical failure. Inter-judge agreement on the review flag:
Gwet's AC1 0.836 across the three, 0.80–0.87 pairwise; Fleiss' kappa 0.287 beside
it, the same prevalence artefact the paper's Appendix K discusses. Replies,
config and analysis are committed under `template_annotation_23092026/screen/pass1/`.
The 24 flags are read in the folder README: presentation defects the gate cannot
see, physical-range concerns for an expert, six rounding-chain claims on lines T1
cannot parse, one claimed logic error, and two judge artefacts of the prompt's
design. Pass 2 runs after human certification.

---

## D-094 — The screen's flags are verified claims, not verdicts; 34 of 45 were real

**Date:** 2026-09-24 · **Status:** DECIDED (owner: "apply fixes properly based on the judges'
reviews") · **Source:** `template_annotation_23092026/screen/pass1_fixes.md`

Pass 1 flagged 24 templates. Each judge sentence was treated as a claim and verified
against the code and generated instances before any edit: 34 claims confirmed and
fixed, 9 rejected with the evidence (among them the one-judge claim that a wave's
direction convention was inverted, and a "calculation discrepancy" that was a
correctly rounded quotient), 2 artefacts of the published prompt showing a function
without its imports. Two findings were physics, not presentation: the virial-work
template mixed pressure-explicit and volume-explicit forms so its printed non-ideal
deviation was an artefact, and the gas-phase concentration template sampled ε
independently of the stoichiometry it stated. Fixes follow D-016/D-037 for any
numeric change; validity conditions (Vr ≥ 2 for the two-term virial, GM > 0 for
upright floating, the elastic range for Poisson's ratio) were added as sampling
constraints with their rejection rates measured and recorded. D-050's deferred
origin-marker fix is applied. The item-pool consequence is large for four templates
and is in `pass1_fixes.md` for the owner.

**Pass 2 is targeted, with carry-forward.** A judge's row is a verdict on a specific
prompt text, and each row records the prompt's hash; where a template's prompt is
byte-identical to pass 1, the pass-1 row is carried into pass 2 marked as carried,
cost zeroed, and only templates whose prompt changed are re-judged. That gives a
complete pass-2 table over all 150 on the shipped corpus for about a sixth of a
full pass's cost, and it is what a false-positive rate against the experts needs:
the panel's verdict on the same version the experts will see. A full third pass is
not needed unless human certification reworks a broad share of the corpus.

**Pass 2 ran 2026-09-24:** 24 templates re-judged (72 calls, $0.70, all parsed on the
first attempt), 378 rows carried. Outcome on the shipped corpus: 147 pass, 3
controversial, 0 critical failure; AC1 on the flag 0.929 (up from 0.836). The three
controversial templates are single-judge flags, two of them repeating claims already
rejected with evidence; they route to the experts as the paper's rule says, and are
not iterated against the judges.

---

## D-095 — Layer 2: certification that leaves evidence — hand checks, planted defects, timestamps

**Date:** 2026-09-24 · **Status:** DECIDED (protocol; roster and dates open) · **Source:**
`template_annotation_23092026/layer2/` (README, guide.md, build_tasks.py, app.py, plants/CONTRACT.md, score.py)

The December 2025 certification (Appendix K) approved 270 of 270 rows with every
mathematical-correctness score at 5, at 12–30 seconds per template, by the same
three reviewers across all three branches. Reviewer gFWV called its perfect kappa
unconvincing and the record could not rebut him: it held no dwell time, no hand
check, no rejection and no negative to detect. Layer 2 is designed so that the
record can.

- **Own-branch experts, three per template**, the pilot's 15 by default (5 branches
  × 3), 34 items each: the branch's 30 templates plus 4 planted defects.
- **A hand check before the solution is shown.** The expert enters their own
  answer to one instance; the app records it, compares within 1%, then reveals the
  solution. The same instance for the three experts of a branch.
- **Planted defects**, four per branch, one of each class — a wrong constant, a
  wrong unit conversion, a flipped sign, a printed step that does not follow from
  its operands — as mutations of real templates verified on every build
  (`plants/CONTRACT.md`). Experts are told quality-control items exist. The
  detection rate is the sensitivity figure a set of approvals cannot give; it is
  reported per expert and overall.
- **Opaque codes, shuffled order, no screen verdict shown, timestamps** on open,
  hand check and submit; a rejection requires a defect type and a note that goes to
  the author.
- **Two routes**, app or workbook, producing identical label rows; no template code
  ships to an expert (instances precomputed at build).
- **Reported by `score.py`:** plant detection, hand-check agreement, Fleiss κ and
  Gwet AC1 on Approve/Reject, AC2 on the scores, the screening panel's false-
  positive rate against the experts and the MAD between their medians, dwell time,
  and the fix list. Rejected templates are fixed and re-judged individually by the
  screen (D-094), not in a new pass.

**Why plants rather than a larger sample.** A perfect approval rate on 150
templates says nothing about the reviewer; a perfect approval rate on 150
templates *and* 16 of 20 planted defects rejected says the review was real. The
plant set is small because each is a hand-written mutation with a measured,
detectable defect, and four per branch is enough to distinguish a reviewer who
reads from one who clicks.

**Open:** the roster (names in `annotators.json` when confirmed), the dates, whether
experts are compensated, and adjudication of split verdicts (proposal: the
majority decides Approve/Reject; a 1-of-3 rejection with a substantive note is
still fixed).

---

## D-097 - The hard case is answered by a digit-level arithmetic check, not by a model

**Date:** 2026-09-24 · **Status:** DECIDED · **Evidence:**
`evaluator_pilot_17092026/analysis/digit_rule.py`, RESULTS_X1 Finding 5

D-096 recorded that no evaluator detects a flawed step behind a correct final answer,
and that buying more labelled hard cases would not change that. It also recorded why
the flaws are there: 164 of the 167 are calculation slips.

The annotation guide's rule for a slip is mechanical - rounding is not an error, a wrong
digit is - so it was implemented: recompute every arithmetic claim from the numbers the
trace itself shows, and reject a displayed value that is not a correct rounding at the
precision shown. On the 228 correct-answer traces it reaches step precision 0.488 and
recall 0.485 (the 72B PRM: 0.240 / 0.266) and trace-level AUROC 0.669, CI 0.604-0.732 -
the only result in the pilot whose hard-case interval clears chance. It calls no model,
so it costs nothing, and every flag names the claim, the value shown and the value
recomputed.

This is E4's own checker with its tolerance changed. E4 ships at 1% relative tolerance,
which is right for catching a fabricated number and blind to every slip the experts
marked: at 1% it reaches recall 0.036 on the same steps. The defect was the tolerance,
not the design, and E4's null result in RESULTS_E4 should be read that way.

Scope of the claim, for the paper: arithmetic flaws behind a correct answer are
deterministically detectable and need no judge. Conceptual flaws behind a correct answer
remain unmeasured - the corpus holds 3 - and no method can be validated on 3 cases.

## D-098 - The final-answer check is corrected and reported offline; E0 is not re-scored

**Date:** 2026-09-24 - **Status:** DECIDED (user's call) - **Evidence:**
`evaluator_pilot_17092026/evaluators/answer.py`, `analysis/answer_check.py`,
`evaluator_pilot_17092026/E0_RERUN.md`

The published framework's final-answer check disagrees with the experts on 72 of 300
traces, 68 of them traces the experts call correct. On the pilot slice it understates
accuracy by about 21 points overall, 0 to 33 per model (corrected 2026-09-27, D-111: this
first read "every model"), and ranks GPT-5 fourth where the experts rank it first.

`evaluators/answer.py` replaces it: the answer segment is read whole, each target is a
quantity the gold COMPUTED rather than any number it prints, every part of a multi-part
answer is scored, and the verdict is correct / partial / incorrect because 19 traces are
genuinely partial. It agrees with the experts on **0.947** of the 281 non-partial traces
against E0's 0.747, and 0.893 three ways over all 300. The relative tolerance is the one
fitted parameter, chosen on one half of the traces and reported on the other (0.876,
0.905), and a self-test pins each fixed defect.

The defect that mattered most was not either of the two on record: the gold value was read
as the last number in the solution, which for `manning_rectangular_discharge` is the 3 in
`m^3/s`. A trace scored correct if it wrote its unit in ASCII and wrong if it wrote the
unicode exponent.

**E0 is not re-scored with it.** The corrected accuracy is a property of the traces and is
computed offline for nothing; the paid re-run would only change which traces reach the
Tribunal, and no finding in the pilot turns on that. Re-running E0 alone would also put it
on a different check from E0-3J and E1, breaking a controlled comparison: "swapping the
judges barely moves the score" holds only while all three arms share one check. So the
pilot reports the corrected accuracy from `analysis/answer_check.py`, and E0/E0-3J/E1 stay
as measured, each labelled as scored under the published framework's own check.

If it is ever wanted, the re-run now costs about **$8.50, not $6.53**: the corrected check
calls 75 more traces correct, and with 288 of 300 traces below the Tier-1 match ratio
threshold, nearly every newly-correct trace reaches the judges - about 233 judged traces
against 178. The harness will run it in 12-50 minutes rather than 3.9 hours: the per-pair
cross-encoder cache never fired (Tier 1's matrix cache sits above it) and was not
output-identical when it did, and judge calls are now fetched concurrently while scoring
stays serial, because the wrong-answer sample reads a process-wide seed.

## D-099 - Electrical's third label set was re-annotated, and the ground truth did not move

**Date:** 2026-09-24 - **Status:** DECIDED - **Evidence:**
`evaluator_pilot_17092026/annotation/rater_diagnostics.py`, RESULTS_X1

The verification round put one electrical expert's self-agreement at 0.433, the lowest of
the fifteen, and the diagnostics showed the signature: "not a claim" used four times as
often as by the branch's other two experts, fewer steps marked incorrect, the odd label in
47 of 77 split steps, and the adjudication going against that set in 66 of 77. Electrical's
between-rater step kappa was 0.584.

That expert re-annotated their 60 solutions (`ele-3.1`; `ele-3` kept under `superseded/`).
The request named the two patterns rather than the numbers. Both moved: "not a claim" 6.1%
-> 0.4%, steps marked incorrect 10.3% -> 22.2%, every one carrying a written reason.
Electrical's kappa is now **0.762** and the pooled figure **0.781** (was 0.750); the branch's
split steps fell from 77 to 47, and the adjudication now goes against that set 4 times in 47.

**The ground truth is byte-identical before and after.** All 2,091 step labels and all 300
trace verdicts are unchanged. Every step where the replaced set could have swung a majority
had already been settled by the blind three-way adjudication, and the new disagreements are
cases where one expert differs from a unanimous pair, which the majority rule absorbs. So no
result in RESULTS_X1, RESULTS_E2 or the hard-case analysis changes; what changed is the
reliability the labels are held to.

This is worth reporting as a robustness property rather than a footnote: an entire expert's
set - a fifteenth of the annotation effort, and the least self-consistent one - was replaced,
and the truth the evaluators are scored against did not move.

**Both loose ends were then closed (2026-09-24, same day).** `ele-3.1` was verified on the
same seven traces as the branch's other two experts - 0.945 step agreement, kappa 0.845,
which lifts electrical's intra-rater figure above the pooled average and covers all fifteen
sets. And electrical was re-adjudicated against `ele-3.1`: all 47 of its current splits were
reviewed blind, 51 rows unanimous and 7 by majority, with the superseded round kept under
`adjudication/superseded/`.

That second round is what finally moved the truth, by 12 step labels, all correct ->
incorrect: the 12 steps that became splits only once `ele-3.1` was labelled, which the
reviewers judged real errors. The truth now holds **388 incorrect steps rather than 376**.
Verdicts and final answers are unchanged, so every trace-level result stands - E0 0.850,
E5 0.886, the expert answer verdict 0.974 - and only the step-level numbers move: the 72B
PRM to precision 0.539 / recall 0.515, the hard-case pool from 87 to 93 traces, and the
digit rule to precision 0.506 / recall 0.472 (trace AUROC 0.655, still the only result
clearing chance on that target).

So the honest form of the robustness claim is two-stage: replacing an entire expert's set
changed no ground-truth label, and re-adjudicating the disputes it created changed 12 of
2,091, none of them a verdict.

## D-100 - E2 keeps the stock 0.5 threshold, and that is measured rather than assumed

**Date:** 2026-09-24 - **Status:** DECIDED - **Evidence:**
`evaluator_pilot_17092026/analysis/prm_threshold.py`, RESULTS_E2

E2 flags a step when its process reward falls below 0.5, the PRM cards' own default. Every
precision and recall in RESULTS_E2 and RESULTS_X1 Finding 3 rests on it, so it had to be
shown that the number was not fitted to the data it is reported on.

The threshold was calibrated the only way that answers the question: maximise step F1
against the expert labels on one deterministic half of the 300 traces, report on the other,
both directions. For the 72B - the PRM E2 reports - the held-out gain is **-0.018 and
+0.006**, both intervals straddling zero, and F1 moves 0.02-0.03 over a +-0.10 band around
the fitted value. 0.5 sits just below that plateau. In-sample optimism, fitting and
reporting on all 300 at once, is +0.019.

So E2 keeps 0.5, and the claim in the paper is not "we used the default" but "the default
was tested against a calibrated alternative and lost nothing".

Two findings fall out of the same analysis:

- **VersaPRM's rewards are mis-scaled as well as undiscriminating.** It rates nearly
  everything above 0.9, so at 0.5 it flags almost nothing (recall 0.057-0.080). Calibrated
  to ~0.9 it recovers recall and still reaches only F1 0.345 against the 72B's 0.520.
- **The hard case cannot be tuned.** Inside correct-answer traces no PRM reaches precision
  0.80 on held-out data at any threshold, and the 72B's two halves choose 0.740 and 0.928 on
  sharp peaks rather than a plateau. The over-flagging RESULTS_X1 Finding 5 reports is a
  property of the model, not of the cut-off.

## D-101 - E4 scores its arithmetic on the digit rule; the 1% tolerance stays, as a second reading

**Date:** 2026-09-24 - **Status:** DECIDED - **Evidence:**
`evaluator_pilot_17092026/evaluators/arith.py`, `evaluators/e4_arith.py`,
`analysis/arith_gold_validation.py`, `analysis/digit_rule.py`, RESULTS_E4 (re-run)

D-097 measured the experts' rule and said E4's null result should be read as a defect in
its tolerance. E4 has now been re-run with the rule in the evaluator (300 gold traces plus
the Gemma column, $0.00, config `0b81b063cab8`; the old rows are kept). Every claim carries
both verdicts - `ok` at 1% ("is this number fabricated") and `ok_digit` ("is this digit
wrong") - and E4 scores on `ok_digit` while reporting the 1% rate beside it.

**The gold validation had to be re-run, and it found two bugs in the rule, not in the gold.**
At 1% the checker flags nothing on the 227 gold claims. The digit rule as D-097 measured it
flags **36** (15.9%): 24 are unit conversions compared before their unit factor
(`2.47 mm = 2.47e-03 m`), and 12 are `lorentz_force` cross products where the gold states an
operand to three figures and computes with the unrounded value, so the numbers it *shows*
cannot pin the result's last digit. Both were fixed rather than excluded - the comparison
happens after the unit factor, and `arith.shown_uncertainty` widens the tolerance by what
the displayed operands leave undetermined - and gold is back to **0 flagged, 227 of 227**.

**What the fixes cost, measured.** Inside correct-answer traces: bare rule precision 0.506 /
recall 0.472, shipped rule **0.750 / 0.320**; trace-level AUROC 0.655 (0.594, 0.717) against
0.639 (0.587, 0.692) - indistinguishable. A rule that flags 16% of gold cannot ship, so the
evaluator takes precision; `digit_rule.py` keeps both and prints the difference.

**What it changes in E4.** Contradicted milestones 2 -> 10 of 1,494, `e4_coverage` down at
most 0.016 per model, milestone precision/recall against the experts unchanged (0.926 /
0.915 - the experts' milestone label asks E3's question). What moves is the arithmetic
score: on the hard case `arith_consistency` goes from AUROC 0.494 (0.459, 0.532) to **0.661
(0.599, 0.726)**, +0.149 over E0 with the difference excluding zero, and as a filter it
selects 43 of 228 traces at precision **0.767** where the 1% rule selected 12 at 0.417.
E3 is untouched: E4's `e3_coverage` equals E3's `milestone_coverage` on all 360 rows.

So RESULTS_E4's "E4 adds nothing over E3" stands for E4's *milestone* arithmetic and fails
for its *per-claim* arithmetic, which is the pilot's only signal on the hard case.

**Re-running E3/E4 needs `pinned_templates.py`.** Template work on 2026-09-23/24 (`3fad887`,
`fc1a6dc`) rewrote how five of the pilot's templates display a derivation, so 17 of the 60
items no longer regenerate byte-identically and `milestones.py` refuses - correctly. The pin
reads those five files back out of git at `26f9048` and installs them for the run; it
relaxes no check, and `--check` reports 60 of 60 reproducing.

## D-102 - Planted defects, so the hard case is measured against a truth we set

**Date:** 2026-09-24 - **Status:** DECIDED - **Evidence:**
`evaluator_pilot_17092026/analysis/planted.py`, RESULTS_X1 Finding 7

Two objections sit on Findings 5 and 6. The digit rule implements the rule the annotation
guide gave the experts, so their agreement is partly built in. And conceptual error behind
a correct answer cannot be studied on this corpus: 3 such steps against 175 calculation
slips.

A planted set answers both. From the 129 traces the experts called clean, 60 get one digit
of one displayed intermediate changed, 60 get the REASONING of one step corrupted with no
digit anywhere in the trace moved, and 60 are left alone. Truth is by construction, every
plant is verified without the checker under test, and the build is seeded and reproduces
byte-identically.

  arithmetic   digit rule 0.750 (as Finding 5 runs it), 0.683 as E4 ships it, 0.533 at 1%
  conceptual   0.000, every evaluator, both readings
  controls     false alarms 0.267 / 0.117 / 0.083

**Nothing detects a conceptual defect.** A trace that computes correctly, misstates the
rule it is applying and lands the right answer passes every deterministic evaluator the
pilot has. That is the open problem, stated now without reference to the guide.

Two things the split settles. The digit rule scores **45 of 45** inside a claim arith.py
parses and **0 of 15** on a value with no parseable working, so Finding 5's recall of 0.472
is a parse ceiling and extending coverage - not changing the rule - is the remaining gain.
And E4's old tolerance is confirmed as the defect by construction: 1.000 on planted errors
of 1% or more, 0.000 on the 14 below it.

The set is a diagnostic, not evidence about model behaviour: it says what an evaluator
catches when a defect of a given shape is present, not how often such defects occur. The
conceptual rules vary in kind but not in phrasing, so it is not a held-out test.

An unrelated defect surfaced while building it, and it should not be lost: **17 of the 60
frozen items no longer reproduce byte-identically from the repo's templates**, so
`milestones.build_all` now raises on the pilot manifest and the build fell back to the
milestone values frozen into the annotation tasks - which is what the pilot's E3 run and
the experts both used, so no number here is affected. But the templates have drifted from
the September freeze, and `pinned_templates.py` (added with D-101) exists to work around
it. Fixing the drift, or pinning the templates properly, is unowned.

## D-103 - The judges do catch conceptual defects; the framework does not ask them

**Date:** 2026-09-24 - **Status:** DECIDED (user approved the spend) - **Evidence:**
`evaluator_pilot_17092026/analysis/planted_judges.py`, RESULTS_X1 Finding 8

D-102 reported that no free evaluator detects a planted conceptual defect, 0 of 60. Every
evaluator it tested reads numbers rather than prose, so the claim could not be general, and
the evaluators that might catch a misstated rule - the judges - had not been asked.

They were, on the framework's own Tribunal prompt with one step under review, matched: each
of the 120 planted defects judged both in the planted trace and in the same step of the
untouched original. 480 calls, **$4.67**, no failures.

  conceptual   GPT-5 0.333, Opus 4.5 0.133 (0.483 counting its `Other` verdicts),
               either judge 0.350 - against 0.000 for every deterministic evaluator
  arithmetic   0.717 each, 0.800 either
  originals    120 of 120 `Alternative Correct`, both judges: no false alarm at all

Two things follow.

**The capability exists and only the judges have it.** A third of conceptual defects is not
a solution, but it is the difference between "no evaluator can do this" and "only a judge
can, a third of the time". The paper's claim changes accordingly.

**Judges and the digit rule are complementary, exactly.** On arithmetic defects inside a
claim `arith.py` parses the digit rule scores 1.000 and the judges 0.667; on defects in a
stated value with no parseable working the digit rule scores 0.000 and both judges 0.867.
The checker is perfect where it can parse and blind where it cannot; the judges are
strongest where it is blind. **The evaluator this argues for is a router**: verify what can
be verified deterministically, and spend a judge call only on the residue - which is what
E5 already does for milestones, applied to steps.

What is NOT established: whether E0 would ever show a judge these steps. Every planted trace
has a correct final answer and E0 samples wrong-answer traces at 0.20, so in a real run most
would never reach the Tribunal. The capability is measured; the routing is the open question,
and it is now a design question rather than a research one.

## D-104 - E0's routing is indifferent to a conceptual defect, so the judges rarely see one

**Date:** 2026-09-25 - **Status:** DECIDED - **Evidence:**
`evaluator_pilot_17092026/analysis/planted_routing.py`, RESULTS_X1 Finding 8

D-103 established that E0's judges catch a third of planted conceptual defects and four
fifths of the arithmetic ones when the step is put in front of them, with no false alarms.
It could not say whether E0 puts it there. Tier 1 was therefore run for real on all 120
planted defects and their unmodified originals, with the Tribunal replaced by a recorder:
163 minutes of local compute, **$0**, no judge called.

                    triggered   step shown   same step, original   end to end
  conceptual          0.650        0.500          0.500              0.175
  arithmetic          0.833        0.767          0.683              0.613

**For a conceptual defect the routing carries no signal at all.** The corrupted step reaches
a judge exactly as often as the untouched one. Tier 1 forwards it because it cannot match
the step to the gold, not because anything about it is wrong, so whether E0 catches a
misstated rule is decided by a draw it was already making. Multiplying through D-103's
detection rates, E0 catches about **18%** of conceptual defects of this shape end to end and
**61%** of arithmetic ones. (Corrected 2026-09-27, D-111: the end-to-end column above is that
product. Counted defect by defect - shown AND caught - it is 9 of 60, **0.150**, and 36 of
60, **0.600**, because the judges catch fewer of the defects E0 shows them;
`analysis/router_planted.py`.)

This is the evidence for the router that RESULTS_X1's recommendation 5 argues for. E0 spends
on judges without aiming them; a checker that verified what it can and sent only the residue
to a judge would aim the same spend at the steps that need it.

A second finding falls out. E0's answer check calls the final answer wrong on **32 of the
120** planted traces, although the plants never touch the answer and every original was
correct. Those traces enter the wrong-answer sample at probability 0.20, so part of E0's
routing today is driven by the answer-check defect D-098 measured and deliberately did not
re-run E0 against. It does not change any reported number; it does mean E0's judge spend is
partly directed by a parser error.

The measurement is a diagnostic: it describes what E0 does with defects of this shape, not
how often models produce them.

## D-105 - The full run's evaluator stack, and no E0 at scale

**Date:** 2026-09-25 - **Status:** DECIDED (user's call) - **Evidence:** the evaluator
pilot, RESULTS_X1 Findings 1-8, `evaluator_pilot_17092026/analysis/judge_cost.py`

What the pilot was for. The full benchmark - 12 models x 2,250 problems = 27,000 traces -
is scored with:

  answer correctness   evaluators/answer.py, three-way correct / partial / incorrect.
                       No model. 0.947 against the experts where E0 scores 0.747, and since
                       the answer nearly determines the trace verdict (0.974) this is the
                       headline metric rather than a component of one.
  milestone coverage   E3, derived from the templates. No model.
  arithmetic           E4 on the digit rule (D-101). No model. Three real flags in four
                       inside correct-answer traces, and clean on gold.
  residual judging     E5: MiMo-V2.5-Pro on the 23.6% of milestones E3 cannot settle.
                       ~$79 over 27,000 traces (18,000 small-model traces at $0.00390,
                       9,000 frontier at $0.00098).

Total evaluation cost **~$79**. The optional step router would add $100-200 and is not
budgeted until MiMo's conceptual detection rate is measured (below). E2's PRMs can run on
HiPerGator for $0 as a secondary signal; they are not a headline number, because they
over-flag inside correct-answer traces (precision 0.246) and their threshold is not the
cause (D-100).

**E0 is NOT run on the full benchmark.** It would cost ~$390 in judges (18,000 x $0.00919
plus 9,000 x $0.02492) and produce numbers with no ground truth to validate them against.
The evaluator comparison belongs to the labelled slice, where 300 traces carry three expert
labels each: that is where "the new stack beats E0" is established, and RESULTS_X1 is the
record. The full run is the benchmark measurement using the evaluator the pilot chose.

Two consequences to carry into the paper, neither of them a cost:

1. **Table 1 will move.** The published numbers came from E0's answer check, which
   understates accuracy by about 21 points and ranked GPT-5 fourth of five where the experts
   and the corrected check both rank it first (D-098). Regenerating the table with the
   corrected check raises every model and reorders the top. A reader comparing versions will
   see that jump, so the paper has to say why.
2. **The design is 150 templates, not 15.** The pilot could detect AUROC differences above
   0.12-0.19 once clustering by template was accounted for (Finding 6). The full benchmark
   has ten times the clusters, so intervals should again be cluster-robust by template, and
   the comparisons the design can support should be settled before the run rather than after.

Open before generation starts:

- **The template drift (D-102).** E3, E4 and E5 all derive milestones from the repo's
  templates, and 17 of the 60 frozen items no longer reproduce byte-identically. At 2,250
  items this has to be fixed or formally pinned.
- **MiMo's conceptual rate.** Finding 8 measured GPT-5 and Claude Opus 4.5, which share
  families with the evaluated roster - the objection E1 exists to answer. The judge that
  survives it is MiMo, and its rate on the planted defects is unknown. ~$3.20 and an hour
  to close, and it decides whether the step router is worth building at all.

**Superseded in part, 2026-09-27.** The roster and the costs by D-110: 11 models, 24,750
traces, E5 about $77, and the router in the stack at $79 batched (untested) or $371 one step
per call. MiMo's rate was measured: 16 of 51 conceptual defects. Point 1's "understates
accuracy by about 21 points" is a pilot-slice figure, and the full-benchmark shift has to be
measured when the table is regenerated (D-111).

## D-106 — Layer 2 round 1: the experts' rejections are verified claims too; 20 templates change, 2 do not

**Date:** 2026-09-25 · **Status:** DECIDED · **Source:**
`template_annotation_23092026/layer2/RESULTS.md`, `layer2/fixes_round1.md`, `layer0/gate_report.md`,
`layer0/item_pool_impact.md`

Round 1 of Layer 2 (15 own-branch experts, 2026-09-24/25) caught all 60 planted defects,
matched the template's answer on 459 of 504 hand checks, agreed at Gwet AC1 0.91, and rejected
22 templates (15 by majority). The screening panel had passed all 22, so its false-positive
rate against the experts is 10.2% (15 of 147). The rejections are the pipeline working: every
one named a defect the gate cannot read and the judges did not see.

**Decision.** Each rejection note was treated as the screen's flags were (D-094): verified
against the source and the seeds the expert saw before any edit. The rule for the fix follows
the kind of defect the experts found, and it is the same rule for all 22:

- **Loads sampled independently of material strength** (five mechanical templates): the load
  or torque is bounded by an allowable stress from two new sampling-only tables in
  `mechanical_engineering/constants.py` (`[POLICY: sampling-only]`, never printed), and the
  materials with no usable design strength are redrawn. Bounds act by redrawing the load
  alone so the material and section shares stay flat (D-045); where a plain redraw would have
  rejected 89% of draws (`multi_segment_rod`) the loads are scaled together instead.
- **Absurd operating conditions** (three chemical transport templates, one vibration template):
  a physical cap with a named source (3 m/s pumped-liquid velocity, Sinnott §5.4.3; 500 kPa;
  base acceleration 1 g), applied as a bounded redraw; no fluid is deleted.
- **Sampled values that should be computed** (`gas_viscosity_kinetic_theory`): Ω_μ from the
  Neufeld–Janzen–Aziz correlation at the drawn T*, so the question's given value is right.
- **Wrong table rows**: molten metals out of `MANOMETER_FLUIDS`; blood plasma corrected to a
  Newtonian 1.2 mPa·s.
- **Presentation** (fixed-decimal displays that lose figures, raw floats, an unstated
  viscosity, an unstated origin convention, "leads" for "lags", "mm/m" for "mm/mm"): displays
  lengthened to the exact value or to four significant figures and bound through the display
  (D-016, D-037); the missing statements added.
- **Two claims not adopted**, with the reason in the docstring: the thermal-stability cap on
  `work_isothermal_virial` and `pitzer_correlation_z`. Both are confirmed as observations
  (organics at 900–1188 K), but there is no on-disk stability source, and a 750 K cap would
  remove 10 of 27 substances from the virial template because at Tr ≤ 1.5 no sampled P2
  satisfies the Vr ≥ 2 validity filter. That is the P2/Pc sampling-range decision already open
  under D-094, now with a second reason to take it.

**Evidence after the edits.** Only the 22 assigned functions differ from HEAD. Gate green,
150 of 150 at 500 seeds; the tie census is unchanged; the item pool moved in 20 templates
(question 2,929 / answer 3,850 / solution 4,718 of 45,000), cumulatively 77 templates since
`b1445ad`. Every number is in `fixes_round1.md` with the command that produced it.

**What follows.** Round 2 re-certifies the 22 templates alone, without plants and with fresh
hand-check instances, by the same experts (`build_tasks --only … --round 2`, `make_kits
--round 2`). The unchanged pair go back with the docstring argument so the rejecting expert
can answer it. Whether the screen re-judges the 20 changed templates first is a few cents and
the owner's call; the runner refuses a third pass, so it would be the targeted re-judge D-094
prescribes.

**Open.** The virial P2/Pc range; a per-row provenance tag for the plasma row (it sits under
the table-level UNVERIFIED tag because a new tag moves the register's chemical count).

## D-107 — Layer 2 round 2: 17 of 22 approved by all three; five objections, each confirmed on the instance cited

**Date:** 2026-09-26 · **Status:** RECORDED (results); the fixes and any round 3 are open · **Source:**
`template_annotation_23092026/layer2/RESULTS_round2.md` (`score.py --round 2 --prev-labels`)

Round 2 (2026-09-25) sent the 22 templates changed after round 1 (D-106) to the same nine experts
of the three branches concerned: 66 reviews, no planted defects, and a hand-check instance none of
them had seen. Seventeen templates are approved by all three experts, two by majority
(`finite_convolution`, `newtons_law_shear_stress`, one rejection each), and three are rejected by
majority (`basic_stress_strain` and `utube_manometer` 3 of 3, `multi_segment_rod` 2 of 3). Of the
43 round-1 rejections of these templates, 36 became approvals and 7 stayed rejections, most of
them for a new reason; 3 verdicts went from approve to reject. Every hand check matched (66 of
66), and Gwet AC1 on the verdict is 0.88.

What the rejections say, each checked against the instance the expert cited (the numbers reproduce):

- **Round 1's objections are resolved, with one exception.** Every defect named in round 1 is gone
  except `finite_convolution`'s origin. The app renders a question as Markdown, so the asterisks
  that mark n = 0 show as italics, and the fallback "x[0] = v" does not locate the origin when v
  repeats (instance 1, x = {2, *-2*, -2}). The models receive the raw text with its asterisks, so
  the item is unambiguous to them; stating the start index removes the dependence on rendering,
  and the app should show a question as typed.
- **The two claims not adopted in round 1 are accepted.** The expert who asked for a
  thermal-stability cap in `work_isothermal_virial` and `pitzer_correlation_z` approves both after
  reading the docstring argument.
- **Four objections are new.** `basic_stress_strain` draws load, diameter and elongation
  independently with no material, so instance 2 implies E = σ/ε = 1423 GPa (439 kN on a 36 mm bar,
  1.17 mm over 3.86 m). `utube_manometer` still pairs room-temperature manometer liquids with
  cryogenic pipe fluids: `PIPE_FLUIDS` holds liquid oxygen and liquid nitrogen, and only this
  template reads it. `newtons_law_shear_stress` has no laminar check, so water at 1.59 m/s across
  1.96 cm is at Re ≈ 3.1e4, where the linear Couette profile does not hold. `multi_segment_rod`'s
  objection was caused by the round-1 fix: scaling the loads to the allowable stress left US-unit
  deformations near 1e-3 in, and the fixed 4-dp display then keeps one or two figures, so instance
  3 sums 0.0036 + 0.0027 − 0.0054 = 0.0009 in where the exact sum is 0.000800 in.
- **The reviews were fast.** Median 1.1 to 1.5 minutes per item and 64 of 66 under two minutes,
  against 1.9 to 3.3 minutes for the same experts in round 1. They knew these templates and the
  notes are specific, but round 1 is the stronger record of review time and the one to quote.

**Open.** Whether a template approved 2 of 3 with a confirmed minority objection counts as
certified (the D-095 adjudication rule, now with two concrete cases); fixing the five and a round
3 for them; the screen re-judge (D-106). None of the five is fixed yet.

## D-108 — Layer 2 round 2 fixes: five templates fixed; the review app's rendering is a finding for rounds 1 and 2

**Date:** 2026-09-26 · **Status:** DECIDED (the five fixes); OPEN (the app, `signal_operations`) · **Source:**
`template_annotation_23092026/layer2/fixes_round2.md`, `layer2/markdown_scan.md`, `layer0/gate_report.md`,
`layer0/item_pool_impact.md`

The owner asked for the five templates D-107 lists to be fixed (2026-09-26). Each objection was
checked against the instance cited, then fixed by the rule its kind of defect calls for:

- **A scenario that no material could produce** (`basic_stress_strain`): the material is drawn and
  named, the load is bounded by its allowable stress, and the elongation follows from its modulus,
  so σ/ε reproduces a real E.
- **Physically impossible pairings** (`utube_manometer`): a stated rule, after Çengel & Cimbala,
  that the manometer liquid is denser than and immiscible with the pipe fluid, applied by
  redrawing the manometer liquid alone (D-031, D-045). Cryogenic pipe fluids and two rows that are
  not manometer liquids (bromine, sodium polysulfide) are never drawn. At HEAD 55.2% of instances
  used a pair the rule now excludes, so the manometer-liquid shares move (mercury 8.7% to 30.9%).
- **A display that lost figures after round 1's fix** (`multi_segment_rod`): deformations at four
  significant figures, the total as the exact sum of what is printed, and a 0.5% guard.
- **An assumption outside its validity** (`newtons_law_shear_stress`): a laminar cap, Re ≤ 1000,
  about 30% under the plane-Couette transition (Tillmark & Alfredsson 1992), applied by redrawing
  the speed below the cap so no fluid is lost and no speed is pinned at the range floor.
- **Meaning carried by rendering** (`finite_convolution`): every sequence states its index list.

**Evidence.** Only the five template functions changed. Gate green, 150 of 150 at 500 seeds; the
tie census drops to 12 templates, with none of the five among them; the item pool moved in exactly these five
(question 1,053 / answer 1,327 / solution 1,355 of 45,000). Answer displays changed in all five;
the owner's instruction to fix them is taken as the D-044 sign-off.

**The finding.** The review app renders questions and solutions as Markdown, and the models read
the raw text. `markdown_scan.py` finds constructs that remove characters in 31 questions and 63
solutions, 65 templates in either. Nearly all are unspaced products such as `2*pi*f*t`, whose
paired asterisks become italics, and dollar amounts, which become inline math. The experts judged
those templates in rounds 1 and 2 from that view. No verdict is known to be wrong because of it,
but none was made on the text the models read. `signal_operations` has `finite_convolution`'s
defect exactly: its origin is marked only by asterisks, in the question and the answer. Round 3's
five templates carry no lossy construct beyond the now-redundant convolution marker, so round 3
can go out as built.

**Open.** Whether the app shows questions and solutions as plain text from now on, and whether
the templates judged through the rendered view get a plain-text look; fixing `signal_operations`.

## D-109 — Layer 2 complete: 150 of 150 templates certified, each by a unanimous latest round on the code as it stands

**Date:** 2026-09-26 · **Status:** RECORDED · **Source:** `template_annotation_23092026/layer2/RESULTS_round3.md`,
`layer2/CERTIFICATION.md` (`certification.py`)

Round 3 (2026-09-26) sent the five templates fixed under D-108 to the nine experts of their
branches: 15 reviews, no planted defects, and a hand-check instance none of them had seen. All five
are approved by all three experts, and the ten round-2 rejections of these templates all became
approvals. Two hand checks were scored as mismatches; both answered the manometer in pascals
(6831 Pa) where the template answers in kilopascals (6.831 kPa). That is the same value, and the
app compares the numbers without their units.

`certification.py` combines the three rounds. A template counts as certified when the latest round
that reviewed it is a unanimous approval by its three own-branch experts and the current code
regenerates, byte for byte, the five instances they were shown. A later edit to a template, or to a
table it reads, therefore voids its certification until it is reviewed again. All 150 are
certified: 128 last reviewed in round 1, 17 in round 2 and 5 in round 3. No split verdict is in
force, so the pool needs no adjudication rule (D-095).

**What the certification rests on.** Fifteen own-branch experts gave 531 verdicts on real templates
over three rounds, plus 60 on planted defects, of which they caught 60 in round 1. A hand check came
before every solution. The rounds produced 53 rejections of real templates, 43 in round 1 and 10 in
round 2; each was verified against its instance, then fixed or, in two cases, answered with an
argument the expert accepted.

**What it does not yet rest on.** The experts judged 65 templates through the app's Markdown
rendering, which dropped multiplication and dollar signs from what they saw, and
`signal_operations` still marks its origin only by asterisks. Both stay open under D-108.

## D-110 — The full run's roster is the pricing document's eleven, and the router is in the stack

**Date:** 2026-09-27 · **Status:** DECIDED (the owner's call, "for now") · **Supersedes in part:**
D-105 (roster size and costs) · **Evidence:** `docs/inference_pricing/build_pricing_doc.js`,
`evaluator_pilot_17092026/PILOT_SUMMARY.md` section 3, `analysis/judge_cost.py`,
`analysis/router_residue.py`

The roster is the eleven models of the inference pricing document as it stands on 2026-09-27,
which the pilot summary reproduces:

- open weight, $212.81: gpt-oss-20b, gemma-4-26b-a4b-it, deepseek-v4.1-flash, qwen3-235b-a22b,
  glm-5.3-flash, glm-5.3, muse-glimmer-30b, kimi-k3
- closed weight, $190.18: gpt-5.4-mini, gemini-3.1-flash-lite, claude-sonnet-5

The pricing list dropped gpt-5.4-nano and claude-haiku-4.5 the same day. NEXT_CYCLE_REVIEW
section 9.5 had proposed reaching eleven by dropping Kimi K3 instead ($288 of inference); this
keeps Kimi K3. Generation is $402.99 on the pricing document's basis - every model writes as
much per problem as GPT-5 did on the pilot, 5,232 output tokens with reasoning - and a
reasoning-heavy model can exceed it (DeepSeek R1's median completion on the pilot was 9,658).

The rule the roster follows is the pricing document's (2026-09-22): no model the pilot generated
traces with, and no model used as a judge. The judge-and-judged objection needs only the second.
The rule leaves out three things NEXT_CYCLE_REVIEW names as anchors: Llama 3.1 70B, which the
expert labels anchor on; the Qwen2.5-Math-7B and Qwen2.5-7B pair behind the math-pretraining
claim (9W1B 5); and any flagship (section 9.3 item 9). They stay open questions for the paper.

**Evaluation, 11 x 2,250 = 24,750 traces, judge calls only:**

| layer | judge | cost | source |
|---|---|---|---|
| answer check, E3 milestones, E4 digit rule | none | $0 | — |
| E5, residual milestones | MiMo-V2.5-Pro | $76.82 | `judge_cost.py`, ROSTER column: 8 open at the pilot's weak-model rate, 3 closed at its frontier rate |
| step router, residue batched per trace | MiMo-V2.5-Pro | $79.10 | `router_residue.py`; the batched prompt is untested |
| step router, one call per residue step | MiMo-V2.5-Pro | $371.38 | `router_residue.py`; the design the planted probe measured |

With a batched router that is $155.92, or $448.20 one step per call. The published framework
on the same basis is $333.57 ($227 to $617 across the two rates) and is not run at scale
(D-105).

**The budget, counted.** About $47 is spent this round (NEXT_CYCLE_REVIEW section 7). With
$402.99 of generation and $155.92 of evaluation the plan comes to about $606 against the ~$500
round, and about $898 with a per-step router, before any paraphrase, tool or flagship
condition. Closing that gap - a cut, a batched router, or a larger round - is the supervisor's
call and is not made here.

**The router** is in the stack as the pilot summary presents it: it is the only component aimed
at conceptual error behind a correct answer, and on planted defects it lifts end-to-end
detection from E0's 0.150 to 0.294 (D-111). It is not built, and its batched prompt has not been
measured. D-105's condition for budgeting it - measure MiMo's conceptual rate - is met: 16 of
51.

## D-111 — The pilot summary said more than its data; each claim corrected against the script that owns it

**Date:** 2026-09-27 · **Status:** DECIDED · **Evidence:** `analysis/router_planted.py`,
`analysis/summary_numbers.py`, `analysis/cluster_bootstrap.py` (section 1b),
`analysis/answer_check.py`, `analysis/judge_cost.py`, `analysis/router_residue.py`,
`analysis/planted_judges.py`

The pilot summary (PILOT_SUMMARY.md and its PDF) was checked number by number against the
analyses. The claims below said more than the data, or something other than it. Each is
corrected in the summary, its figures and RESULTS_X1, and the number that replaces it is printed
by the script named. Overstating a result is not acceptable in this project in any document, and
this entry is the record of what was wrong.

1. **"The final answer beats the best evaluator by +0.124, with a clustered interval excluding
   zero."** +0.124 is the margin over E0. Over E5, the best evaluator, it is +0.088 with a
   template-level interval of −0.002 to +0.181, which includes zero. Over every other evaluator
   the interval excludes zero (`cluster_bootstrap.py`, section 1b).
2. **The router "shows the judge the flawed step every time" and lifts conceptual detection "to
   about 0.35", "doubling" it.** The first was assumed and never measured. The router forwards 51
   of the 60 conceptual plants: the other 9 sit in steps holding a claim that verifies, which it
   settles without a judge. With MiMo it catches 15 of 51 end to end, 0.294. The 0.35 was E0's
   two-judge union, not the router's judge (`router_planted.py`).
3. **E0's end-to-end detection, 0.175 conceptual and 0.613 arithmetic, was a product of two
   rates.** Counted defect by defect it is 9 of 60 (0.150) and 36 of 60 (0.600): the judges catch
   fewer of the defects E0 shows them. This corrects D-104's table too.
4. **"The benchmark understates every model by about 21 points"; "every model's accuracy rises
   by roughly 21 points".** 21 is the pooled gap on the slice; per model it is 0 to 33
   (`answer_check.py`). The slice was built to over-represent the answer types the comparator
   reads worst - 20% of its items are scalar, against 59% of the pool's templates - so the size of
   the full-benchmark shift has to be measured. This corrects D-098 and D-105 too.
5. **E5 is "more accurate" than E0, "the most accurate reasoning metric".** E0 has no milestone
   score to compare; at the trace level E5 does not beat E0; and on correct-answer traces with a
   flawed step E3 to E5 score AUROC 0.39 to 0.43, below chance. The summary had left that result
   out, and now reports it.
6. **The digit rule "scores 1.000 where it can parse the claim".** True of the bare rule, 45 of
   45. The rule E4 ships catches 40 of 45, and the five it gives up run from 1e-7 to 3e-5
   relative - not all "below one part in 100,000" (`router_planted.py`).
7. **The false-alarm row set a per-trace rate beside per-step ones.** The digit rule's 0.117 is
   any flag anywhere in a clean trace; the judges' 0.000 is on one step each. On the same 120
   untouched steps the digit rule flags 3.
8. **The router was priced at $79 and credited with the probe's detection rate.** The $79 is for
   a batched prompt nobody has tested; the detection comes from single-step prompts, and that
   design costs $371 (`router_residue.py`).
9. **The opening paragraph and section 2.3.** "A deterministic stack scores the same traces more
   accurately" holds only for the answer check. "A judge is worth paying for in exactly one
   place" contradicted the stack, which pays for two. "Judge identity does not matter, so the
   panel is not buying accuracy": the first half is measured; the second is not, and E0 does beat
   its own answer check at the trace level (0.850 against 0.812, significant when traces are
   resampled; RESULTS_X1 Finding 1). "This answers the independence objection": it is first
   evidence, not a per-judge bias test.
10. **Smaller ones.** The design effect was "about 3.7" (it is 2.6 to 4.2 across evaluators); the
    detectable difference "0.12 to 0.19 and no smaller" (0.06 for E0-3J and 0.004 for E1); "no
    evaluator separates from E0" (E1 does, by a negligible −0.002); "a third of the templates
    answer with a word" (2 of 15; 8 of 15 put no number on the answer line); "each of the 129
    clean traces receives one defect" (120 do; the 60 controls are untouched traces, 58 of them
    sources of a plant); VersaPRM "7%" (6%); "the framework already had this checker" (E4 had
    it; the published framework has no arithmetic check); "blind to every slip" at 1% (it catches
    3.4%); "45 of the 68 have the correct value in the answer line" (no committed script prints
    it; `answer_check.py` now prints that the corrected check calls 56 of the 68 correct); "about
    13 hours of expert time" (13.4 hours recorded by the app with a trace open, not counting the
    two later rounds); self-consistency above between-expert agreement (pooled only: civil and
    industrial invert it); MiMo's conceptual rate "of 60" (of 51); the published framework "judges
    every trace" (it judges the steps Tier 1 cannot match); the closed-weight roster models headed
    "frontier" (they are the small closed tiers).
11. **The figures typed their values in.** `make_summary_figures.py` now reads
    `figures/summary_numbers.json`, which `analysis/summary_numbers.py` computes from the
    analyses, so a figure cannot disagree with the scripts. The routing figure's legend no
    longer covers a bar.

## D-112 — The two label-review follow-ups closed: one probe step relabelled, one guide change withdrawn, and a gap in the router's rule

**Date:** 2026-09-27 · **Status:** DECIDED · **Evidence:** `analysis/judge_probe.py` (`RELABEL`;
`report --offline`), `analysis/e2_analysis.py`, `analysis/router_residue.py` (last section),
the adjudicated labels in `experts_filled_labels/version_2/` (local only)

RESULTS_X1 closed with two follow-ups from reading the experts' labels. Both are closed here,
and closing the first exposed something larger than either.

**1. Claude's `aoq_ati_rectifying#3` step is a slip, and is relabelled.** The judge probe had
it CLEAN. All three experts mark it incorrect: it writes 50 · 0.070 · (0.930)^49 = 50 · 0.070 ·
0.02857 = 0.1000, and 0.93^49 is 0.0285538, so P(X=1) is 0.0999. `judge_probe.RELABEL` now
carries it, and the probe (13 slips, 8 clean) and E2's validation are re-scored:

| | as first scored | after the relabel |
|---|---|---|
| probe leaders (balanced) | MiniMax, MiMo 0.96 | Kimi K3, Grok 4.6 0.96 |
| MiMo, MiniMax | 11/12 caught, 0/9 false, 0.96 | 11/13, 0/8, 0.92 |
| GPT-5, GLM-5.3 | 10/12, 1/9, 0.86 | 11/13, 0/8, 0.92 |
| Opus 4.5 | 10/12, 0/9, 0.92 | 10/13, 0/8, 0.88 |
| E2 probe AUROC, 72B / 7B / VersaPRM | 0.935 / 0.824 / 0.565 | 0.923 / 0.779 / 0.548 |

The five judges charged with a false alarm on this step were right, and MiMo, MiniMax and Opus
missed it. The 72B passed it at 0.959, so it misses all three non-Llama slips in the probe.

**It does not change the choice of judge.** The probe was never the deciding criterion, and it
still cannot separate the judges that parse (one item moves "balanced" by about 0.05, and the
relabel is one item). E5's judge was chosen for independence, and the two judges the relabel
lifts to the top cannot take its place: Kimi K3 is on the roster (D-110), and Grok 4.6 has
closed weights and documented exposure to OpenAI, whose gpt-5.4-mini and gpt-oss-20b are on the
roster. MiMo's documented exposure, to Claude, remains the smallest of any candidate. On the
larger matched test since, 120 planted defects, MiMo detects about as much as GPT-5 with no
false alarm (D-103); Grok was not run on it, and that run (about $3.40) is the one to buy if a
reviewer presses, not a re-reading of 21 steps.

**2. The guide change is withdrawn, because its case was misread.** The follow-up said GPT-5's
`normal_depth_iteration#1` step has a notational slip with a right value, which the experts
called correct. It has a wrong sign in its displayed formula, and a wrong value as well: it
writes 1.86921 where the formula gives 1.86942, because 0.046272/2.8172 is 0.016425, not the
0.01621 it uses. After adjudication all three experts mark it a calculation error with that
arithmetic as the reason; only one expert's first pass had called it correct. The guide already
covers it - "a wrong digit is a calculation error, even if the final answer survives it" - the
digit rule flags it, and the pilot holds no step with a notational slip and a right value. The
guide is unchanged. RESULTS_E4's and JUDGE_SELECTION's descriptions of the step, which said the
value is right, are corrected.

**3. What the relabelled step exposed: a gap in the router's rule.** Neither digit rule flags the
`aoq_ati_rectifying#3` slip, because it sits inside a chained equality the checker does not
compare, and the router, as specified in `router_residue.py`, would settle the step without a
judge, because another claim in it verifies. Placing every step the experts call incorrect:
inside correct-answer traces, 57 of 178 are flagged by the digit rule, 63 are residue a judge
sees, and **58 (32.6%) are never flagged and never shown to a judge** - 31 whose error the checker
does not reach, and 27 the bare digit rule flags but the shipped rule's operand-uncertainty
widening passes. Over all 300 traces it is 103 of 388 (26.5%). A router that settles a step on
any verifiable claim would never show a judge a third of the slips behind a correct answer, so
its rule has to forward every step the checker has not positively cleared. That enlarges the
residue beyond the 63% the $79 and $371 were priced on, and it goes into the router's design.

**What this changes in the evaluation stack.** The components stand: the answer check, E3
milestones, E4's digit rule, E5 with MiMo, and a step router with MiMo. The judge stands. What
changes is the router's routing rule, still to be designed, and two documented descriptions -
the probe's ranking and the GPT-5 step - which are corrected in JUDGE_SELECTION, RESULTS_E2,
RESULTS_E4, FINDINGS E0-F5, RESULTS_X1 and the pilot summary.

## D-113 — The router is batched or it is not run; batched, it holds up, and costs about $114

**Date:** 2026-09-27 · **Status:** DECIDED (the owner's call on the design; the smoke checks run
with approval, $2.30) · **Supersedes in part:** D-110 (the router's cost) · **Evidence:**
`analysis/router_batched.py`, RESULTS_X1 "The router, batched"

**The design.** The full run can afford the router only batched, one judge call per trace. One call
per residue step, about $371, is not coverable, so the choice is the batched router or none. The
routing rule is rule C (D-112): every step the digit rule, as E4 ships it, does not flag goes to
MiMo, a trace's steps together, on the framework's own Tribunal prompt.

**The two smoke checks**, called as the single-step probe was called, every call answered after
retries (a laptop shutdown cut four calls mid-run; `router_batched.py` now gives up on a call after
900 s and asks it again):

- *Planted defects*, 240 calls, $0.92. MiMo catches 19 of the 58 conceptual defects it is sent, 19 of
  60 end to end against E0's 9; on the 49 also judged one step at a time it catches 16 either way,
  and 14 of the 19 arithmetic defects sent against 12 one at a time. 2 false alarms on 842 clean
  steps. **Batching costs no detection.**
- *The 300 labelled traces*, 299 calls, $1.38. The router - the digit rule plus the batched judge -
  finds 234 of the 388 steps the experts call incorrect, precision 0.707 and recall 0.603, against
  the digit rule's 0.825 and 0.255. Almost all of the gain is in traces whose answer is not
  correct; inside correct-answer traces recall moves from 0.320 to 0.360 and the trace AUROC from
  0.639 to 0.675 (E0 0.542).

**The cost.** A batched call cost $0.0038 on the planted traces and $0.0046 on the labelled ones, not
the $0.0034 D-110 assumed, so the router over the full run is about $114 at the labelled rate ($95 at
the planted one), and evaluation with it about $190. The budget of D-110 becomes, with this round's
spend at about $49: $49 + $403 generation + $77 E5 + $114 router = about $643 against the ~$500 round.
Whether the router is run at that price is the supervisor's call; what the pilot now says is what it
buys - step-error recall more than doubled, mostly outside the hard case, and the only conceptual
signal in the stack.

**Still open:** building the router in the harness, and what it reports in the results table.

## D-114 — Certification closed; the full run's pool is 150 x 15 from a private seed, and only hashes are committed

**Date:** 2026-09-28 · **Status:** DECIDED (the owner's calls on the certification and the pool size; the seed rule as recommended) · **Evidence:** `full_run_28092026/freeze.py`, `full_run_28092026/FREEZE.json`, `template_annotation_23092026/layer0/gate_report.md`, `template_annotation_23092026/layer2/CERTIFICATION.md`

**The templates are final.** The owner, 2026-09-28: "template annotation has been done over three
rounds and properly closed". No template changes from here. D-108's open items (a plain-text look
for the 65 templates judged through the Markdown rendering, `signal_operations`'s origin marker),
D-106's targeted re-screen and D-094's residual redraw rates are closed without further action.
What the experts saw is recorded in `layer2/markdown_scan.md`, and the paper's account of the
certification should match it. At HEAD (`1acc7cf`) the gate passes 150 of 150 and all 150
templates are certified, each regenerating byte for byte the instances its experts reviewed.

**Why a private seed.** The templates and the generator's default master seed are public, so a
committed per-item seed publishes the item: rerun the template with it. NEXT_CYCLE_REVIEW 9.1
item 3 had said to commit the seed, which is corrected there. The default seed also produced the
pilot slice, whose items are public with their gold and on whose traces the answer check was
tuned. The pool is therefore drawn from a random 128-bit seed held in `SEED.secret`, gitignored
with `pool/`. Committed: `manifest.jsonl` (per-item SHA-256 over question + NUL + solution, no
text, no seed) and `FREEZE.json`, with the seed's SHA-256 as a commitment. Reviewers get the pool
in the ARR supplementary archive. At publication the seed is revealed and the pool released.
Publishing the pool now could not contaminate the models evaluated now, which were trained
earlier; the rule protects the pool after release and lets the paper say the evaluated
items were not public when the models ran.

**Selection.** Indices 0, 1, 2, ... per template; an index is replaced by the next when its
question repeats one already kept, matches a pilot-slice question, carries an exact display tie
on a line T1 can parse (the tie census's test, D-016), or fails to generate. This carries out
D-092's instruction to census the frozen pool and fix stragglers without editing a certified
template. 56 indices were replaced: 51 repeated questions and 5 display ties, 4 of them in
`server_configuration_selection`. None matched a pilot question.

**15 per template, the owner's rule.** The owner, 2026-09-28: "150 x 15 = 2250, pool should have
2250 instances in total". Two templates cannot supply 15 distinct questions in 75 indices:
`adiabatic_flame_temperature` has 11 and `heat_of_reaction_formation` 4. Each keeps every distinct
question once and fills to 15 with repeats in index order, 4 and 11 of them, each marked
`repeat_of` in the pool and the manifest. The pool is **2,250 instances and 2,235 distinct
questions**: 870 Easy, 870 Intermediate, 510 Advanced, 450 per branch. The paper should state both
numbers; "2,250 unique problems" would be wrong.

**Verified.** `freeze.py --verify` in a separate process regenerates all 2,250 byte-identically,
and `--check-files` matches `pool/` to the manifest. Manifest SHA-256 `f0ffa108e79312e6…`, seed
commitment `cd53376fb0fe90ac…`, both in full in `FREEZE.json`.

**Open.** A private backup of `pool/` and `SEED.secret`; the tag, on the commit inference runs at.

## D-115 — Two chemical templates widened so each can supply 15 distinct problems; round 4 for the chemical experts

**Date:** 2026-09-28 · **Status:** DECIDED (the owner's call to fix and re-certify); round 4 OPEN · **Evidence:** `template_annotation_23092026/layer2/fixes_round4.md`, `round4_checks.py` / `round4_checks.md`, `full_run_28092026/DIVERSITY.md`

The diversity analysis found two templates that could not produce 15 distinct questions from any seed:
`heat_of_reaction_formation` had 4 and `adiabatic_flame_temperature` 11, so D-114's freeze filled them
with repeats. The owner, 2026-09-28: "i think we should fix these templates and we can ask the chemical
annotators to look these up again ... 2 templates are doable". This reopens the certification for these
two templates only.

- **`heat_of_reaction_formation`** draws from a new table, `HESS_REACTIONS`: the four `REACTIONS` rows
  plus 21 reactions built only from species already priced in `HEATS_OF_FORMATION`. The question, the
  reasoning and the level are unchanged. `REACTIONS` itself is untouched, because three stoichiometry
  templates also read it: widening it was tried first, and the gate caught generation errors in all
  three. The constants suite now checks atom balance over `HESS_REACTIONS`.
- **`adiabatic_flame_temperature`** samples 0-100% excess air in whole percent, 1,111 fuel-and-level
  cases, with exact decimal coefficients; at 0% it reads as before. Its documented margins were measured
  at theoretical air, so all were re-measured: six passes within 0.102 K of the fixed point, contraction
  0.089-0.135, 1430-2908 K, and at most 2.41 K from a NIST Shomate solve, so Step 6 now says 2.5 K. The
  rounding margin did not carry over, so a guard redraws the excess level in the 20 cases where the
  answer's last kelvin would depend on the iteration count or a half-kelvin display.

Gate 150 of 150; the tie census and the Markdown scan are unchanged, and neither template appears in
either; certification is 148 of 150 until round 4 returns. Round 4 is built: the two templates for the
three chemical experts, no plants, fresh hand-check instances.

**Open.** Round 4's verdicts; the pool is re-frozen with these two templates as they now stand (D-116).

## D-116 — The pool is re-frozen with coverage selection: 15 per template, spread across its reasoning paths and answer forms

**Date:** 2026-09-28 · **Status:** DECIDED (the owner approved "coverage selection, then re-freeze") · **Supersedes in part:** D-114 (the selection rule) · **Evidence:** `full_run_28092026/freeze.py`, `FREEZE.json`, `DIVERSITY.md`

D-114 took each template's first 15 acceptable draws, so a template's rarer branches and labels
entered the pool only by chance: Reynolds' regime came out 12 turbulent, 3 laminar, 0 transitional.
The rule now examines the first 100 draws, groups the acceptable ones by reasoning path and answer
form (the lower reading of `diversity.py`: equation lines and answer segment with numbers and
question-varying words masked), and takes the 15 round-robin across the groups, lowest index first.
The exclusions are D-114's (repeated question, pilot question, display tie, generation error). A
template with one group keeps exactly the instances D-114 froze.

**What it changed.** Re-frozen from the same seed, after the two widened templates (D-115): 2,250
items and 2,250 distinct questions, no repeats needed; 332 items in 93 templates differ from D-114's
pool, 29 of them in the two widened templates. Measured by `diversity.py` on the new pool:

- templates with one reasoning path in all 15 instances: 54 to 58, where D-114's pool had 56 to 61
  and 500 public draws show 54 to 58. The pool now contains a second path for every template that
  can produce one within 500 draws.
- classification labels: Reynolds 8 turbulent, 6 laminar, 1 transitional (was 12, 3, 0); damping
  5, 5, 5 (was 7, 6, 2); linearity 7 linear, 8 not (was 10, 5); memory and causality 6, 5, 4 across
  its three combinations (was 7, 6, 2).
- near-duplicate pairs at 5%: 61 (was 71).

**What it costs, and the paper must say.** Rare branches are over-represented relative to how often a
template produces them, so any per-branch or per-label statistic describes the pool's design, not the
template's natural mix. `critical_depth_froude_classification` still shows one answer form: its answer
line states only the Froude number, so the class goes unscored. That is a template defect, left open.

`--verify` in a separate process regenerates all 2,250 byte-identically; the seed commitment is
unchanged. If round 4 changes either widened template, its 15 items are regenerated and re-frozen.

## D-117 — The analysis plan is written before the data

**Date:** 2026-09-28 · **Status:** PROPOSED as a whole (binds when the owner confirms, and in any case before the first inference call); its two scoring rules DECIDED by the owner the same day · **Evidence:** `full_run_28092026/ANALYSIS_PLAN.md`

Reviewers in both cycles asked for variance, significance and fixed thresholds (yAYU 2 and 4, 9W1B 3
and 4); a plan written after the results cannot answer them. `ANALYSIS_PLAN.md` fixes, before any
inference: the template as the unit of analysis, with every interval resampling templates; five
confirmatory questions with their tests (model accuracy and Holm-corrected pairwise McNemar tests; the
Easy-to-Advanced cliff with its detectable gap, 12 to 18 points for a template standard deviation of 0.20
to 0.30; what process scores add on wrong- and right-answer traces; consistency within a template, split
by single- and multi-path templates; the paraphrase test on 450 items with its paired test and ranking
stability); the sensitivity analyses; and what is reported but not tested. The owner decided its two
scoring rules: an unusable trace scores 0, defined by whether the check can read a final answer rather
than by the marker, with persistent service failures reported as missing and never scored; and a partial
answer scores 0.5, so the headline is an answer score, reported beside the fully-solved rate. The May
version's math-pretraining claim is dropped, as the roster has no math-specialised
model (D-110).

## D-118 — Correction to D-116: the Froude template's answer line is not a defect

**Date:** 2026-09-28 · **Status:** RECORDED · **Evidence:** `data/templates/branches/civil_engineering/water_resources/energy_and_rapidly_varied_flow.py`

D-116 said `critical_depth_froude_classification`'s answer line "states only the Froude number, so the
class goes unscored" and called it a template defect. It is not. The question asks to "determine the
Froude number of the flow" and, "in your solution", to compute the critical depth and use it to justify
whether the flow is subcritical or supercritical. The answer line states the Froude number, which is what
the question asks for as its answer; the regime is a justification step inside the solution, the kind of
step process evaluation exists to read. The claim came from the diversity report listing the template
under classification with a single answer form, read without opening the template.

What is off is metadata: the audit inventory labels the template's answer type `classification`, for a
template that answers with a number, so it is counted among the classification templates in the pool's
answer-type tally and the diversity report. Relabelling it `scalar` needs no template change and no
expert round, but the inventory is read by code, so what reads that label is checked first, in the gold
validation. The analysis plan's sensitivity item for this template is withdrawn.

## D-119 — Round 4: both widened templates approved by all three chemical experts; 150 of 150 certified

**Date:** 2026-09-28 · **Status:** RECORDED · **Source:** `template_annotation_23092026/layer2/RESULTS_round4.md` (`score.py --round 4`), `layer2/CERTIFICATION.md` (`certification.py`, rounds 1-4)

Round 4 sent the two templates widened under D-115, `heat_of_reaction_formation` and
`adiabatic_flame_temperature`, to the three chemical experts: no planted defects, and a hand-check
instance none of them had seen (seed 2401). All six verdicts approve, all six hand checks match the
template within 1%, and Gwet AC1 on the verdict is 1.000. Their previous verdicts, all approvals, are
round 1's, so round 4 is scored without a previous-round comparison, which `score.py` would have
headed "round 3".

The reviews were quick: a median of 1.8 to 2.0 minutes per template, 4 of the 6 under two minutes, for
templates that had changed substantively. The hand checks, computed by each expert before the solution
is shown, all matched, and they are the evidence the templates were worked rather than skimmed.

**Certified: 150 of 150**, each by a unanimous latest round on the code as it stands: 126 last reviewed
in round 1, 17 in round 2, 5 in round 3 and 2 in round 4. The round-4 tasks were built from the working
tree at `fb88cf1` carrying the fixes, which were committed straight after as `03b1dd7`;
`certification.py` confirms the current code regenerates every reviewed instance byte for byte. The pool
(D-116) was frozen from the same code, and `freeze.py --verify` regenerates all 2,250 items
byte-identically, so no re-freeze follows.

**Not done:** the screening panel has not re-judged the two templates; a targeted re-judge costs a few
cents if wanted (D-106's precedent).

## D-120 — The full pool's gold validation, and the three evaluator gaps it closed

**Date:** 2026-09-28 · **Status:** DECIDED · **Evidence:** `full_run_28092026/gold_validation.py`, `GOLD_VALIDATION.md`, `evaluator_pilot_17092026/evaluators/answer.py` and `arith.py`, the pilot's analyses re-run under `pinned_templates`

The deterministic evaluators - the answer check, E3 and E4's digit rule - were run over the 2,250
gold solutions, whose right answer is known, before any trace exists. All 2,250 items regenerate
byte-identically from their templates, and E3 finds every milestone in its own gold. Two evaluators
had gaps the pilot's 15 templates never exercised:

- **The answer check called 20 gold answers wrong, in three templates.** Two are the classification
  templates the pilot excluded as surface-predictable (D-057): "linear" and "not linear" were outside
  its vocabulary, and a labelled yes/no answer ("Memoryless: No, Causal: Yes") could not be read,
  since its word pattern skips two-letter words and its last-verdict-word rule cannot score two
  labels. The third is one item answering with a bare pi. Fixed: linear and nonlinear form their own
  family; a labelled yes/no is scored label by label, reading negated forms ("non-causal", "has
  memory"); pi is read as a value.
- **The digit rule flagged 103 of 6,970 gold claims, in seven templates.** Sentences listing
  assignments - "34 kN at a = 3.3 m", "D3 = 0 for n = 4", "from t = -2 to t = 2" - were read as the
  chain 34 = 3.3; a correct 58 ft = 696 in was flagged; an '=' inside a subscript label split a
  claim. Fixed: a connector word that opens a new assignment ends the claim; ft and in convert; a
  label's '=' is protected.

After the fixes: 2,250 of 2,250 gold answers correct in all six answer types, and 0 of 6,917 gold
claims flagged. Each rule fires only on the new forms, and the pilot's own analyses, re-run under the
pin, print what they printed before - `answer_check`, `arith_gold_validation` and `planted`
identical, every expert-agreement figure and the digit rule as shipped unchanged - except two
secondary readings in `digit_rule`: at the 1% and 0.1% tolerances the trace AUROC moves from 0.480 to
0.479 and from 0.527 to 0.526, corrected in RESULTS_X1 and RESULTS_E4.

**Found and not fixed.** Three templates have no milestones (`gauss_law_symmetric`,
`system_properties_memory_causality`, `system_property_linearity`), so milestone coverage is
undefined for their 45 items, which are left out of milestone aggregates rather than scored 0
(ANALYSIS_PLAN Q3); 52 templates have some item with a single milestone. E3's raw null, a sibling
item's gold scored against an item's milestones, is 0.123 over the pool against the pilot's 0.072,
and highest in `signal_operations` at 0.70.

**D-118's relabel is not needed.** The answer check scores `critical_depth_froude_classification`
the same under either label, because its answer line holds one number and no verdict word, and the
comparator bindings carry their own copy of the label, so the inventory is left as it is.

## D-121 — The full run's harness, and what its dry run found

**Date:** 2026-09-28 · **Status:** DECIDED (the harness); OPEN (`qwen3-235b-a22b`; approval of the check and the calibration run) · **Evidence:** `full_run_28092026/run_traces.py`, `models.json`

`run_traces.py` is the pilot's runner generalised: the same prompt, whose hash it checks against the
pilot traces', and the same call, so the stack validated on those traces reads the same kind of
output. It reads the frozen pool and refuses to run unless the pool matches the committed manifest;
ends each (item, model) as answered, empty or service failure, the D-117 rules; reads the billed cost
from every response rather than estimating it; routes open weights to the cheapest endpoint serving
fp8 or better and closed models by price; leaves decoding at each provider's default, as the pilot
did; and sets a 32,768-token ceiling. Every mode that bills refuses to start without `--yes`.

The dry run reads OpenRouter's public endpoint lists, with no key and no model call. On the pricing
document's basis, 270 input and 5,232 output tokens per item, the run is $402.98 for 24,750 calls;
at today's cheapest eligible endpoints it is $319.09 for ten models. Two endpoints changed since the
pricing document:

- `qwen3-235b-a22b` has one endpoint left, Alibaba's, with no declared quantization, output capped at
  8,192 tokens, at $0.455 and $1.82 per million against the document's $0.087 and $0.35. No endpoint
  meets the routing rule, so the paid modes skip it until the owner decides.
- `muse-glimmer-30b` has one eligible endpoint, DeepInfra at bf16, capping output at 16,384, which is
  its ceiling.

On the same basis `--check` bills about $0.01 and `--calibrate 20`, 220 calls, about $3.58; the
calibration replaces the assumed length with measured ones, and its traces count toward the run.

## D-122 — Inference starts on ten models, at provider-default decoding, with the analysis plan confirmed

**Date:** 2026-09-28 · **Status:** DECIDED (Qwen skipped for now; decoding; the analysis plan; the check and the calibration approved; the tag); OPEN (the full run's approval) · **Evidence:** `full_run_28092026/run_traces.py`, its dry run

The harness as committed in D-121 could not be imported: the `\n` escape opening the paid modes'
skip message (`run_traces.py`, line 351) had been written as a literal line break, a `SyntaxError`
that stopped every mode, the dry run included. Rejoining the line is the only change. The dry run
then passes: 2,250 items in 150 templates match the manifest, the prompt is the pilot's (sha256
`c2bcb87984c4e50b`), and the run is $402.98 on the pricing document's basis ($4.17 of it Qwen's) and
$334.65 for the ten models at today's cheapest eligible endpoints (D-121's dry run gave $319.09).

The owner's decisions, 2026-09-28:

- **`qwen3-235b-a22b` is skipped for now.** The run goes ahead with the other ten. Keeping it on
  Alibaba's endpoint or replacing it stays open; either adds its traces without touching the other
  models'.
- **Decoding stays at each provider's default**, temperature and reasoning effort both, as for the
  pilot traces the evaluator stack was validated on. The inference guide had listed temperature as
  open.
- **The analysis plan (D-117) is confirmed as a whole.**
- **The check and the calibration run are approved** on D-121's estimate. With Qwen skipped,
  calibration is 200 calls.
- **Inference runs at the commit tagged `full-run-inference`** (D-114).

A correction to D-121: `--check` does not skip a model with no eligible endpoint; it calls every
model it is given. A Qwen call would carry the fp8-or-better filter that no endpoint meets, so the
check runs on the ten by name.

## D-123 — Calibration: 200 of 200 answered, $1.80 billed; the full run re-estimated at about $202

**Date:** 2026-09-28 · **Status:** DECIDED (the record); OPEN (the full run's approval) · **Evidence:** `full_run_28092026/calibration_estimate.py` (`--providers` for the endpoints), `run_traces.py --status`, the traces (local)

The check, one tiny call per model at the tagged commit, answered for all ten and billed about
$0.006. The calibration then ran 20 items per model, 200 calls, one model at a time with 16 workers:
200 answered, none empty, none stopped at the output cap, no service failure and no call needing a
retry, so the rows' $1.797 is the calibration's whole bill. Its traces count toward the run.

`calibration_estimate.py` projects each model from those bills: what is billed plus the mean bill per
item times the 2,230 items left. The interval is a bootstrap over the 20 items, so it assumes they
represent the pool. They are one item from each of 20 templates, 15% Advanced against the pool's 23%;
reweighting the same bills by level gives $198.15 instead of $202.15.

| model | output tokens per item | billed $ | full run $ | 95% interval | by level $ | assumed $ | hours at 15 workers |
|---|---|---|---|---|---|---|---|
| gpt-oss-20b | 2,284 | 0.004 | 0.48 | 0.33-0.65 | 0.49 | 1.07 | 1.3 |
| gemma-4-26b-a4b | 958 | 0.007 | 0.74 | 0.57-0.96 | 0.73 | 3.59 | 0.4 |
| deepseek-v4.1-flash | 3,581 | 0.067 | 7.57 | 4.66-11.37 | 7.37 | 5.03 | 1.0 |
| glm-5.3-flash | 3,509 | 0.025 | 2.78 | 1.82-3.91 | 2.74 | 2.99 | 2.6 |
| glm-5.3 | 5,264 | 0.419 | 47.09 | 29.70-66.42 | 45.18 | 52.65 | 1.4 |
| muse-glimmer-30b | 3,379 | 0.083 | 9.28 | 7.80-10.83 | 9.32 | 13.13 | 0.8 |
| kimi-k3 | 3,365 | 0.715 | 80.49 | 49.40-126.95 | 77.70 | 130.18 | 3.5 |
| gpt-5.4-mini | 673 | 0.064 | 7.24 | 6.13-8.43 | 7.36 | 53.43 | 0.2 |
| gemini-3.1-flash-lite | 646 | 0.021 | 2.33 | 1.99-2.69 | 2.34 | 17.81 | 0.1 |
| claude-sonnet-5 | 1,889 | 0.393 | 44.16 | 36.09-52.49 | 44.92 | 118.93 | 0.7 |
| **ten models** | | **1.797** | **202.15** | **163.47-251.32** | **198.15** | **398.81** | |

The full run for the ten comes to about $202, about half the $398.81 the same ten cost on the
pricing document's basis, which is the basis of the plan's $403 generation line for all eleven
(D-110, D-113). Nine of the ten wrote less than
the 5,232 output tokens assumed, the closed models and Gemma far less; GLM-5.3, at 5,264, was close
to it. DeepSeek V4.1 Flash is the one model over its assumed cost: price-sorted routing with
fallbacks served its 20 rows from four endpoints, 3 of them from Morph, the cheapest. GLM-5.3's came
19 from Sail Research and 1 from Morph, not from Novita, which the dry run listed as cheapest. The
requests carry the fp8-or-better filter, so the endpoint moves the price, not the rule. With ten
processes in flight at once in the full run, the spread over endpoints, and so the price, may differ
from the calibration's. The hours assume each call takes as long as the calibration's did.

**Open.** The full run's approval, on this estimate.

## D-124 — The full run starts on seven models; Kimi K3, GLM-5.3 and Claude Sonnet 5 wait

**Date:** 2026-09-28 · **Status:** DECIDED (the seven); OPEN (the three) · **Evidence:** D-123; `calibration_estimate.py --model KEY ...`

The owner approved the full run for seven of the ten: `gpt-oss-20b`, `gemma-4-26b-a4b`,
`deepseek-v4.1-flash`, `glm-5.3-flash`, `muse-glimmer-30b`, `gpt-5.4-mini` and
`gemini-3.1-flash-lite`. On the calibration's bills the seven come to $30.41, with a 95% bootstrap
interval of $26.67 to $34.66; the three held back, `kimi-k3`, `glm-5.3` and `claude-sonnet-5`, come
to $171.74, interval $134.07 to $221.32 (the script with each set's `--model` keys). The three keep
their 20 calibration traces, which count toward the run whenever it resumes for them.

The seven started 2026-09-28 at 17:27 local time, one detached process per model with 15 calls in
flight, the harness unchanged since the tag `full-run-inference`. A spend line sits at the seven's
upper bound, $34.66: if the bill passes it, or the projection does once a quarter of the calls are
in, the run stops for the owner.

**Open.** Whether and when the three run.

## D-125 — The run stopped at the spend line, and resumed with the line at $45.06

**Date:** 2026-09-28 · **Status:** DECIDED · **Evidence:** `calibration_estimate.py --empties --until 2026-09-28T12:34:24Z` with the seven `--model` keys; OpenRouter's account usage

At 17:33 local time, with a quarter of the seven's calls in, their projected cost passed D-124's line
of $34.66, and the six processes still running were stopped at 17:34. At the stop, $9.76 was billed
across the seven, calibration included, and the projection stood at $41.42, 95% interval $37.98 to
$45.06.

- Gemini 3.1 Flash-Lite had finished all 2,250 items, at $2.48, inside its calibration interval of
  $1.99 to $2.69.
- gpt-oss-20b, Gemma, Muse and GPT-5.4 mini projected inside their calibration intervals.
- DeepSeek V4.1 Flash projected $15.48 against its calibration's $7.57, and GLM-5.3 Flash $4.88
  against $2.78. The harness submits items in the pool's sorted order, so a slow model's first rows
  come from the first templates alphabetically, not from the calibration's spread. Among them is
  `adiabatic_flame_temperature`, where each of the two had 3 items end empty at the output cap.
- Muse's 8 empty rows are all `aoq_ati_rectifying`, at its 16,384-token ceiling. The analysis plan
  scores an empty row 0 (D-117).
- OpenRouter's account usage rose $9.67 between a read before the launch and one after the stop,
  against $9.49 in the rows the pass wrote. The difference is at most the cost of the calls in flight
  when the processes were stopped.

The owner resumed all six at 17:45 with the line moved to the projection's upper bound, $45.06: the
two models' projections rest on their first templates, and only more of their items can correct them.

## D-126 — The seven are complete: 15,750 of 15,750 items, $33.23 in the rows

**Date:** 2026-09-28 · **Status:** DECIDED (the record); OPEN (the traces' Kaggle copy) · **Evidence:** `calibration_estimate.py --empties` with the seven `--model` keys; `run_traces.py --status`; OpenRouter's account usage; the traces (local)

The six resumed processes ended by 22:18 local time. Each of the seven has a row for every one of
its 2,250 items, answered or empty, and none is left as a service failure, so no re-run was needed.

| model | answered | empty | output tokens per item | billed $ | calibration's estimate $ (D-123) | attempts beyond the first |
|---|---|---|---|---|---|---|
| gpt-oss-20b | 2,197 | 53 | 3,343 | 0.690 | 0.48 | 20 |
| gemma-4-26b-a4b | 2,250 | 0 | 963 | 0.726 | 0.74 | 4 |
| deepseek-v4.1-flash | 2,236 | 14 | 4,753 | 7.041 | 7.57 | 7 |
| glm-5.3-flash | 2,208 | 42 | 4,901 | 4.512 | 2.78 | 27 |
| muse-glimmer-30b | 2,221 | 29 | 3,658 | 10.045 | 9.28 | 10 |
| gpt-5.4-mini | 2,250 | 0 | 721 | 7.737 | 7.24 | 0 |
| gemini-3.1-flash-lite | 2,250 | 0 | 690 | 2.480 | 2.33 | 0 |
| **seven** | **15,612** | **138** | | **33.230** | **30.41** | **68** |

The seven's $33.23 is inside the approval's first interval, $26.67 to $34.66 (D-124), and under the
$45.06 line (D-125). Two models ended above their calibration intervals: gpt-oss-20b, which wrote
3,343 output tokens per item against the calibration's 2,284, and GLM-5.3 Flash, $4.51 against an
upper end of $3.91. DeepSeek V4.1 Flash, whose early projection stopped the run, ended under its
calibration estimate.

The rows record the attempt each row kept. OpenRouter's account usage rose $35.75 over the run, from
a read before the launch to one after the end, against $32.96 in the run's rows. The $2.79 the rows
do not show covers the 68 attempts beyond the first, the calls in flight when the run was stopped
(D-125), and any other use of the account in those hours, which cannot be told apart from here. The
run itself therefore cost between $32.96 and $35.75, on top of the calibration's $1.80 and the
check's $0.006.

138 rows are empty. 136 of them ended at the output cap with no answer text; the other 2
(gpt-oss-20b, `finite_convolution`) ended with the model stopping without any. The analysis plan
scores an empty row 0 (D-117). They spread over 25 templates for gpt-oss-20b, 10 for GLM-5.3 Flash,
7 for Muse and 4 for DeepSeek; Gemma, GPT-5.4 mini and Gemini have none. Seven more rows stopped at
the cap after some answer text, and are scored on what they state.

The traces, the seven's and the three held models' calibration rows, are archived locally beside the
pool's backup with their checksum, every member checked against its source by SHA-256. The private
Kaggle copy waits for the owner, because the upload needs manual mode.

**Open.** The traces' Kaggle copy; the three held models (D-124); Qwen (D-121).

## D-127 — The three held models run: Kimi K3, GLM-5.3 and Claude Sonnet 5

**Date:** 2026-09-28 · **Status:** DECIDED · **Evidence:** D-123, D-124

The owner approved the full run for the three held back in D-124, on their calibration estimate of
$171.74, 95% interval $134.07 to $221.32, their 20 calibration items each included. They started
together at 23:27 local time, one detached process per model with 15 calls in flight, the harness
unchanged since the tag `full-run-inference`, once the laptop was back on mains power. The run stops
for the owner only if the three's recorded bill passes $221.32, the interval's upper end. A
projection is not a trigger this time: the pool's sorted order makes an early projection run high
(D-125).

## D-128 — Two Qwen candidates are calibrated: Qwen3-235B-A22B-2507 and Qwen3.8-27B

**Date:** 2026-09-29 · **Status:** DECIDED (the calibration); OPEN (the Qwen slot, D-121) · **Evidence:** `full_run_28092026/models.json`, `run_traces.py --dry-run` with the two `--model` keys

`qwen3-235b-a22b` has no endpoint meeting the routing rule (D-121). OpenRouter's public list, read on
2026-09-29 with the harness's own routing function, shows two candidates the owner chose to
calibrate:

- `qwen/qwen3-235b-a22b-2507`, the same model's July 2025 release, without thinking: 6 eligible
  endpoints, the cheapest GMICloud at fp8 for $0.087/$0.35 per M, the price the pricing document gives
  the original `qwen/qwen3-235b-a22b`.
- `qwen/qwen3.8-27b`: 8 eligible endpoints, the cheapest DeepInfra at bf16 for $0.15/$1.875 per M. It
  made traces in the evaluator pilot's robustness cohort, not as an evaluator.

Both endpoints allow more than the 32,768-token ceiling. `models.json` gains the two entries; the
existing entries, `run_traces.py`, the prompt and the pool are unchanged since `full-run-inference`,
and the two run at the commit tagged `full-run-inference-qwen`. The dry run puts their calibration,
20 items each, the same items as the others', at about $0.23 on the pricing document's basis. The
owner approved the check and the calibration. Neither model is part of the roster until the owner
decides the Qwen slot, and their full runs need their own approval.

## D-129 — The two Qwen candidates calibrated; their full runs started under the owner's $320 ceiling

**Date:** 2026-09-29 · **Status:** DECIDED (the calibration; the start); OPEN (the Qwen slot, D-121) · **Evidence:** `calibration_estimate.py` with the twelve `--model` keys, at 02:35; `run_traces.py --status`

After a check that answered for both, the calibration ran 20 items each at the tag
`full-run-inference-qwen`: 40 of 40 answered, none empty, none at the output cap, no retry, and
$0.443 billed.

| model | output tokens per item | billed $ | full run $ | 95% interval | hours at 15 workers |
|---|---|---|---|---|---|
| qwen3-235b-a22b-2507 | 1,304 | 0.010 | 1.07 | 0.75-1.52 | 2.4 |
| qwen3.8-27b | 8,480 | 0.433 | 48.72 | 32.18-68.20 | 8.3 |

The owner approved their full runs on one condition: that the twelve models' total inference cost
stay at or under $320. At 02:35, with Kimi K3 at 926 of its 2,250 items, the script put the twelve at
$286.60, 95% interval $269.66 to $306.30: the finished models' bills, Kimi's projection from its rows,
and the two calibrations. That interval treats Kimi's rows so far as representative, which the pool's
sorted order does not guarantee (D-125), and the rows leave out retried and stopped attempts (D-126).
Both runs started at 02:36, one detached process each with 15 calls in flight. The run stops for the
owner if the twelve's recorded bill passes $320.

## D-130 — The three held models and Qwen3-235B-A22B-2507 are complete

**Date:** 2026-09-29 · **Status:** DECIDED (the record) · **Evidence:** `calibration_estimate.py --empties` with the four `--model` keys; `run_traces.py --status`; the traces (local)

Claude Sonnet 5 finished at 00:21, GLM-5.3 at 01:48 and Kimi K3 at 08:54, local time; Qwen3-235B-A22B-2507
finished at 04:15. Each has all 2,250 items recorded, answered or empty, and none is left as a
service failure. Every row that stopped at the output cap used the full 32,768 tokens.

| model | answered | empty | output tokens per item | billed $ | calibration's estimate $ | attempts beyond the first |
|---|---|---|---|---|---|---|
| claude-sonnet-5 | 2,250 | 0 | 2,566 | 59.447 | 44.16 | 2 |
| glm-5.3 | 2,150 | 100 | 6,515 | 55.064 | 47.09 | 7 |
| kimi-k3 | 2,245 | 5 | 3,866 | 99.000 | 80.49 | 17 |
| qwen3-235b-a22b-2507 | 2,250 | 0 | 1,400 | 1.177 | 1.07 | 0 |

The three held models came to $213.51, under their $221.32 line (D-127) and above their $171.74
estimate. Claude Sonnet 5 ended above its interval of $36.09 to $52.49, writing 2,566 output tokens
per item against the calibration's 1,889; GLM-5.3 and Kimi K3 ended inside theirs, and so did
Qwen3-235B-A22B-2507.

GLM-5.3's 100 empty rows cluster on a few templates: every item of `qr_policy_one_iteration` and
`work_isothermal_virial`, 14 of `vdw_solve_for_volume`, 12 of `adiabatic_flame_temperature`, and 10
each of `aoq_ati_rectifying` and `normal_depth_iteration`. These are largely the templates on which
the other models' empty rows fall too (`--empties`).

## D-131 — Inference is complete for twelve models: $304.62 in the rows, under the $320 ceiling

**Date:** 2026-09-29 · **Status:** DECIDED (the record); OPEN (the Qwen slot; the traces' Kaggle copy) · **Evidence:** `run_traces.py --status`; `calibration_estimate.py --empties` with the twelve `--model` keys; OpenRouter's account usage; the traces (local, backed up)

Qwen3.8-27B finished at 10:13 local time: 2,161 answered and 89 empty of its 2,250, $56.705 billed,
9,148 output tokens per item, 28 attempts beyond the first, inside its calibration interval of
$32.18 to $68.20 (D-129). All 98 of its rows that stopped at the output cap used the full 32,768
tokens. Its empty rows fall on the same templates as the other models': 15 of `qr_policy_one_iteration`,
14 each of `adiabatic_flame_temperature` and `work_isothermal_virial`, 12 of `normal_depth_iteration`
and 10 of `vdw_solve_for_volume` among them.

Every one of the twelve models has all 2,250 items recorded, answered or empty, with no service
failure left: 27,000 rows, 332 of them empty. The original `qwen3-235b-a22b` stays skipped (D-121).
The rows record $304.62 across the twelve, their calibrations included, under the owner's $320 ceiling
(D-129); the checks add under a cent. OpenRouter's account usage rose $316.96 from the read before the
seven started (D-126) to the end, on top of the $1.80 of check and calibration billed before that read,
so $318.76 by the account. The $14.14 the rows do not show covers the attempts beyond the first, the
calls in flight when the run was stopped (D-125), the checks, and any other use of the account in those
hours, which cannot be told apart from here.

The twelve models' traces are archived locally beside the pool's backup with their checksum, every
member checked against its source by SHA-256. The private Kaggle copy waits for the owner, because the
upload needs manual mode.

**Open.** Which Qwen candidate fills the Qwen slot (D-121); the traces' Kaggle copy.

## D-132 — The roster rule stands: Qwen3-235B-A22B-2507 fills the Qwen slot, and Qwen3.8-27B is set aside

**Date:** 2026-09-29 · **Status:** DECIDED · **Evidence:** D-110; `docs/inference_pricing/build_pricing_doc.js`, the roster's exclusion list; D-128 to D-131

The owner keeps D-110's roster rule: no model the pilot generated traces with, and no model used as a
judge.

- **`qwen3-235b-a22b-2507` fills the Qwen slot**, closing D-121. It is the same model's July 2025
  release, without thinking, at the price the pricing document gives the original; the pilot did not
  generate with it, and it is not a judge. It replaces the original like for like, because the
  original's endpoints no longer meet the routing rule. Its full run is complete (D-130).
- **`qwen3.8-27b` is set aside.** The pilot generated traces with it: it is one half of the
  robustness pair that the pricing document's exclusion list names, with `gemma-4-31b-it`, and the
  pilot's evaluators were checked on those two models' traces (R-F1, R-F2). Its 2,250 traces are kept,
  labelled, outside the main results, and can serve only as a robustness check.
- **No other Gemma is added.** `gemma-4-31b-it`, the only newer Gemma on OpenRouter, falls under the
  same rule, and the Gemma 3 models are older and smaller.

The roster the analysis scores is D-110's eleven with the 2507 release in the Qwen slot: 11 × 2,250 =
24,750 traces, the size the evaluation was costed for (D-113).

**Correction to D-128 and D-129.** D-128 described Qwen3.8-27B's part in the pilot, but not that
D-110's rule excludes it: the candidates were read against the routing rule only, not the roster
rule. Its calibration and full run, $56.705 in the rows, were therefore spent on a model the rule
leaves out. The twelve-model totals in D-129 and D-131 stand as the record of what was billed.

## D-133 — The traces reviewed before scoring: every row the frozen item, one per item, one model per key

**Date:** 2026-09-29 · **Status:** RECORDED · **Evidence:** `full_run_28092026/trace_review.py`, `TRACE_REVIEW.md`

An independent check of the twelve trace files, before anything is scored. Every model has a final row
for each of its 2,250 items and no item has two; every row's item hash matches the manifest; no line is
malformed; every row carries the pilot's prompt hash; and each model key has one set of request
parameters and one served model id, so no model changed under the run. The newest local archive holds
the twelve files byte for byte.

The roster's eleven hold 24,750 final rows, 243 of them empty, and $247.92 in the rows; the
set-aside Qwen3.8-27B adds $56.705, D-131's total.

What the scoring has to carry:

- **Empty rows concentrate.** 153 of the 243 fall on six templates, out of 165 rows each:
  `qr_policy_one_iteration` 41, `work_isothermal_virial` 28,
  `normal_depth_iteration` 23, `adiabatic_flame_temperature` 22,
  `aoq_ati_rectifying` 21 and `vdw_solve_for_volume` 18. Four of them are iterative by
  construction. An empty row scores 0 (D-117), so on these templates the score measures finishing within
  the output ceiling as well as solving, and the paper should say so where it reports them.
- **Muse Glimmer ran at a 16,384-token ceiling**, its only eligible endpoint's cap, where the others had
  32,768 (D-121), so its empty rows are not strictly comparable with theirs.
- **Among answered rows**, 18 stopped at the cap after some answer text and are scored on what they
  state; 29 carry no answer marker the check looks for; 3 carry neither a marker nor a
  number at the end, the candidates for the plan's "never states one", which the answer check decides.
- **No `<think>` block** leaks into any answer text.
- **Open-weight models were served by several endpoints**, up to 11 for
  `glm-5.3-flash`, all under the fp8-or-better filter. A quantization change can move outputs, so the
  paper should say the endpoint varied within the rule.

## D-134 — The full-run scorer: every evaluator's raw output in one store, clean on gold, matching the pilot's expert agreement

**Date:** 2026-09-29 · **Status:** DECIDED · **Evidence:** `full_run_28092026/score.py`, `validate_scorer.py`, `SCORER_VALIDATION.md`

The deterministic stack runs once per trace and keeps what each evaluator returned, so every later
analysis reads the store instead of re-scoring: `scores/<variant>/<model>.jsonl`, gitignored with the
traces. A row holds:
- the item's metadata;
- the trace's metadata and the SHA-256 of its text;
- the answer check's label at the fitted tolerance and at half and double it, with the targets it looked for;
- E3's reached flag per milestone;
- per step, the claims the digit rule checked and flagged.

Steps are split by `e2_prm.steps_of`, the split the experts labelled, so a later judge's or router's
step scores refer to the same steps. The evaluators are imported unmodified, and each run records their
hashes and the commit in `CONFIG.json`.

- **Variants.** `main` reads `traces/<model>.jsonl`; any other variant reads `traces/<variant>/<model>.jsonl`.
  A paraphrased trace is scored against its original item, so its gold, milestones and question are the
  original's. The two arms of Q5 then differ only in what the model wrote.
- **Unusable** follows D-117. It is an empty row, or an answered one from which the answer check can read
  no final answer: no answer marker, and no number or verdict word in the last 700 characters.
- **Gold.** All 2,250 gold solutions score correct at all three tolerances, and none is unusable. E3 finds
  every milestone in the 2,180 items that have one. The digit rule flags none of 6,625 claims in 9,391
  steps. That count is below GOLD_VALIDATION.md's 6,917 because here a claim is read only within a step;
  neither check flags one.
- **The pilot's 300 traces** were scored by the same function as the full run and compared with the
  experts' labels. All ten published agreement figures reproduce:
  - the answer check's 0.947 and 0.893;
  - E3's precision, recall and F1;
  - the digit rule's precision and recall, on all traces and in the hard case.

  The answer labels equal the pilot's own call on 300 of 300. The digit rule's tp/fp/fn equal
  `digit_rule.py`'s: 99/21/289, and 57/19/121 in the hard case.
- **The free main pass** has scored all twelve models, 2,250 final rows each. E5, the paid judge on the
  milestones E3 does not find, is the next stage. It reads this store and writes beside it.

## D-135 — Correction to D-120 and ANALYSIS_PLAN Q3: 70 items have no milestones, not 45

**Date:** 2026-09-29 · **Status:** RECORDED · **Evidence:** `full_run_28092026/GOLD_VALIDATION.md` (its per-template counts), `score.py --gold`, `analyze.py`

D-120 and the plan's Q3 say that "45 items, the whole of three templates" have no milestones. Those three
templates are `gauss_law_symmetric`, `system_properties_memory_causality` and `system_property_linearity`,
45 items in all. But 25 more items, in 12 other templates, also have none, so the total is 70.
GOLD_VALIDATION.md's own per-template list shows a 0 in those templates' counts.

Both documents also say that 52 templates have "some item with a single milestone". GOLD_VALIDATION.md
counts items with at most one milestone: 52 templates have such an item, and 46 of them an item with
exactly one. The rule is unchanged: these items are left out of milestone aggregates, not scored 0. The
plan carries a dated note.

## D-136 — The analysis script: where the plan is silent, the method is fixed before any result was read

**Date:** 2026-09-29 · **Status:** DECIDED · **Evidence:** `full_run_28092026/analyze.py` (its docstring and `--selftest`)

`analyze.py` computes Q1 to Q5 from the score store alone, together with the sensitivity analyses and the
tables the plan reports without testing. It stops unless the store holds the pool the plan describes: 150
templates of 15 items, 58 Easy, 34 Advanced and 58 single-path. It is committed before its first run on
the scores. Where the plan names a quantity but not its method, the script fixes it:

- **Q2's p-value** comes from permuting the tier labels among the 92 Easy and Advanced templates. It sits
  beside the plan's within-tier bootstrap interval.
- **Q3's groups.** A wrong-answer trace scores 0 (incorrect or unusable), and a correct-answer trace is
  fully solved; a partial answer is in neither. A trace counts as flagged by the digit rule when any of its
  steps is. Until E5 runs, coverage is E3's, which is E5's deterministic part, and E5's columns are empty.
- **Q1's fully-solved check** agrees with the answer score on a pair when both tests hold at Holm-adjusted
  0.05, or neither does. When both hold, they must also point in the same direction.
- **Q5.** The paired difference is the item mean, and the sign flips act on each template's summed
  difference. Kendall's τ's interval resamples templates, as the plan asks of every interval.
- **Resampled p-values** are (count + 1) / (draws + 1). Each test has its own seed, so a re-run prints the
  same numbers.
- **Per-label accuracy** on classification templates takes the gold's label from its answer targets. The
  Froude template's label is the regime its Froude number implies (D-118).
- **The plan's fourth sensitivity**, without the two round-4 templates, does not arise, because round 4
  returned and certified both (D-119).

The self-test checks each statistic on synthetic data with a known answer. Among its checks, it plants a
paraphrase loss in three of eleven synthetic models and finds it in those three and in no other.

## D-137 — The answer check and E3 misread LaTeX numbers; fixed, and the pilot's agreement rises to 0.982

**Date:** 2026-09-29 · **Status:** DECIDED · **Evidence:** `full_run_28092026/parser_fix.py`, `PARSER_FIX.md`, `validate_scorer.py`, `SCORER_VALIDATION.md`; `evaluator_pilot_17092026/analysis/answer_check.py`; the self-test in `evaluators/answer.py`

The first run of `analyze.py` on the free scores showed an anomaly: DeepSeek V4.1 Flash solved 3 of 8
turbulent Reynolds items. Its traces were right. It wrote `\(Re \approx 3.04 \times 10^5\)`, and the
answer check read that as three numbers, 3.04, 10 and 5.

- **The answer check** (`answer.values`) had an ordering bug. Its exponent rule ran before `\times` and
  `\cdot` were rewritten, so the LaTeX forms never reached it. The docstring says the forms are handled.
  The order has been this way since `cb84f32`, the commit that produced the published 0.947.
- **E3's reader** (`milestones.numbers`) never rewrote `\times` or `\cdot` at all.
- **Both readers** split thousands written `28\,570`, `11{,}003` or with a narrow space, and the answer
  check did not read `\tfrac`.

The fix is minimal. Both readers now read these forms as one number, and 14 new self-test cases pin
them. A plain space still separates numbers, because `\log_2 256` is not 2,256. That leaves one residual
form, "14 056": it appears in 20 of Muse Glimmer's answers and in one each of three other models'.

**Checks.**
- **Gold.** All 2,250 gold answers still score correct at all three tolerances. None of the 2,250
  milestone sets changes.
- **The pilot's 300 traces.** 10 answer verdicts change, and all 10 now equal the experts' verdict.
  - The answer check's agreement goes from 0.947 to **0.982** on the traces the experts judged fully right
    or wrong, and from 0.893 to **0.927** three ways.
  - E3's precision, recall and F1 go from 0.926/0.915/0.921 to 0.927/0.920/0.923.
  - The split-half tolerance fit still gives 0.0015 and 0.0020, so REL stays at 0.002.
  - The pilot's DeepSeek R1 moves from 0.817 to 0.950 against the experts' 0.917. RESULTS_X1 carries a
    dated addendum.

  `validate_scorer.py` reproduces every published figure with the code it was published with. It reports
  the corrected figures beside them, and the paper should cite the corrected ones.
- **The full run** (`PARSER_FIX.md`). Changed verdicts range from 1 to 116 per model. Examples of the
  score changes:
  - DeepSeek V4.1 Flash: 0.941 → 0.974.
  - GPT-5.4 mini: 0.811 → 0.838.
  - Kimi K3: 0.942 → 0.963.
  - gpt-oss-20b: 0.807 → 0.826.
  - Claude Sonnet 5 and GLM-5.3: unchanged at three decimals.

  Some verdicts move from correct: 18 for GPT-5.4 mini, 10 for gpt-oss-20b, 8 for Kimi K3. Every such
  case read had matched through a stray exponent digit, and the model's answer was wrong. For example,
  the "10" of `1.10 \times 10^4` passed for 10,739 at the ×1000 scale, where the stated 11,000 is 2.4% off.

The pilot analyses that also call `milestones.numbers` (`e3_grid.py`, `e3_null.py`, `hard_case_pool.py`,
`judge_probe.py`) would read the corrected numbers if re-run.

## D-138 — No credit from a subscript or a bare 0 or 1: the answer check's match rule tightened

**Date:** 2026-09-29 · **Status:** DECIDED (the owner may reverse it: one commit) · **Evidence:** `full_run_28092026/match_audit.py`, `MATCH_AUDIT.md`; the self-test in `evaluators/answer.py`

Reading D-137's changed verdicts turned up an older weakness. `answer.match` accepts a stated number
within one unit of its own last digit, and it applies that after scaling by a unit factor of up to
10^12. A number whose one-unit window reaches zero therefore matches every gold at some scale. A bare
`1`, for example, is 1,000,000 ± 1,000,000 at 10^6. `answer.values` also read subscripts as numbers,
so `p_1`, `x_{1}` and `N_0` put exactly such digits into answers. In the full run, a formula with no
value, `P_b = Q(\sqrt{2E_b/N_0})`, "matched" 15.15 through the 0 of `N_0`.

**Two rules**, measured alone and together against the D-137 check before either was adopted:
- **S.** A digit written as a subscript is not a value. E3's reader already skips such digits.
- **R.** A number vouches for a gold at its own precision only when one unit of its last digit is
  smaller than the number itself. The relative tolerance and the gold's own precision are unchanged.

| | gold correct | pilot: verdicts changed | full run, 12 models: verdicts changed |
|---|---:|---:|---:|
| S | 2,250 | 0 | 66 |
| R | 2,250 | 0 | 243 |
| S+R, adopted | 2,250 | 0 | 252, all away from correct |

The experts cannot arbitrate: no pilot verdict moves under either rule, so the pilot's agreement stays
at D-137's 0.982 and 0.927. The full-run changes were read instead.
- **123 are numeric answers** (scalar, multipart, array, vector). Every one read was a wrong value that
  had been credited through a stray digit. Examples: 0.238 m for 0.151 m via `y_1`; 1,700 K for
  1,594 K via the 1 of a space-split `1 700`; a queue's W via the 0 of `P_0`.
- **129 are symbolic answers**, mostly `ber_estimation_mary`, `impulse_response_from_lccde` and
  `autocorrelation_rect_pulse`. Number matching cannot verify a formula under either rule. Some formulas
  are right, such as `36(8-|τ|)` for a gold that states 288 = 36 × 8, and they drop from correct to
  partial. Their credit had come from stray digits, not from their numbers. The paper should say that
  symbolic answers are scored by the numbers they state.

Answer scores fall by 0.001 to 0.028 per model. gpt-oss-20b falls most, 0.826 → 0.798, and Claude
Sonnet 5 falls 0.967 → 0.962. A one-significant-figure rounding still counts: "3.4 x 10^3" for 3,364
and "0.5" for 0.4987 are among the new self-test cases. What R still allows is a single digit from 2
to 9 at its own precision, such as a formula's `3` for 3.79; that residual is left as found.

## D-139 — A verdict word is decided within its own family, a negated mention skipped

**Date:** 2026-09-29 · **Status:** DECIDED · **Evidence:** `full_run_28092026/word_audit.py`, `WORD_AUDIT.md`; the self-test in `evaluators/answer.py`

Reading the first results also showed a defect in how the answer check reads the word of a
classification answer. The last mention decides, so that an answer restating the criteria is not
credited for every regime it names. But every verdict word competed with every other, and a negated
mention counted. Claude Sonnet 5's correct "ζ ≈ 1.00 (Critically Damped); ω_d = 0 rad/s (no oscillation
occurs)" ended on `no`, and "... since the system is not underdamped" ended on `underdamped`.

The rule adopted has two parts:
- **Families.** Only words of the target's own family compete: damping, flow regime, yes/no and so on.
  A word outside every family competes with all, as before, and `linear`/`nonlinear` keep their own
  handling.
- **Negation.** A mention negated just before it ("not", "non-", "n't", "never", "neither", "nor") is
  skipped.

Measured before adoption, alone and together:
- The gold is unchanged at 2,250.
- No pilot verdict moves.
- In the full run, 25 verdicts change, all partial → correct and all in `damping_classification`. There
  are 0 to 5 per model, out of its 75 classification items.

The self-test pins both defects, and also a wrong regime that is still wrong: "not critically damped
but overdamped" against a critically damped gold.

## D-140 — The deterministic stack's results on the free pass, and what the evaluator fixes changed

**Date:** 2026-09-29 · **Status:** RECORDED · **Evidence:** `full_run_28092026/results/RESULTS.md` and `results.json` (`analyze.py`); `results/main_pre_d137/` (the same, before the fixes)

`analyze.py` ran on the eleven models' 24,750 traces, scored with the evaluators of D-137 to D-139.
These are the answer check, E3 and the digit rule; E5 has not run, and neither has Q5.
- **Q1.** Answer scores run from 0.799 (gpt-oss-20b) to 0.969 (DeepSeek V4.1 Flash).
  - 32 of the 55 pairwise differences hold at a Holm-adjusted p below 0.05 on the template-level test.
  - None holds among the top five: DeepSeek V4.1 Flash, Claude Sonnet 5, Kimi K3, GLM-5.3 Flash and
    Muse Glimmer, which lie within 0.012 of each other.
  - McNemar's check agrees on 46 pairs. On the other 9 it holds and the template-level test does not,
    because it treats items as independent (D-111).
- **Q2.** The Easy-to-Advanced gap holds for 7 of 11 models, with gaps from +0.047 to +0.194.
- **Sensitivity.** The model ordering barely moves with the tolerance: Kendall's τ is 0.93 at half
  tolerance and 0.96 at double. It barely moves with the fully-solved rate either (0.89), and does not
  move without the shortcut templates (1.00). It does move when unusable traces are excluded (0.60),
  because GLM-5.3's 100 empty rows are then no longer counted.

**What the fixes changed.** Under the evaluators before D-137, the Q2 gap held for 2 of 11 models. The
fixes moved five models' verdicts into it. For all eleven, the Advanced score fell and the Easy score
rose or held. The fixes were found by reading these results, so both versions are kept: the one before
the fixes as a labelled record, `results/main_pre_d137/`. The paper should say that the answer check was
corrected after the first scoring pass, and why (D-137 to D-139). It should also report the Q2 result
under both versions.

**Added after the first run**, descriptive and outside the plan's tests (D-136 did not have them):
- Q3's claims checked per trace, beside each digit-rule rate. The rule reads 1.68 claims per answered
  trace for Kimi K3 and 8.35 for Qwen3-235B, so a low flag rate can mean little was read.
- E3 coverage on the readable wrong answers alone, since 100 of GLM-5.3's 123 wrong-answer traces are
  unusable.
- The note on McNemar's nine disagreements.

## D-141 — Variant runs, the paraphrase pipeline and the expert check: built and tested, nothing billed

**Date:** 2026-09-29 · **Status:** DECIDED (the runs are the owner's to start, by `EVALUATION_GUIDE.md`) · **Evidence:** `full_run_28092026/subsamples.py`, `run_traces.py --variant`, `paraphrase.py`, `paraphrase_kit.py`, `paraphrase_app.py`, `paraphrase_guide.md`; the self-tests named below

**The subsamples** are fixed in one module that the harness, the pipeline and the analysis all read.
- **Paraphrase:** the plan's 1st, 6th and 11th item of each template in manifest order, 450 items,
  90 per branch.
- **Repeat:** the 1st and 8th of each template, 300 items. The plan says only "2 per template", so
  this choice is fixed now, before any repeat runs.

The two share 150 items, each template's first. Positions are counted, not indices, because the
indices have gaps (D-116).

**The harness.** `run_traces.py --variant paraphrase | repeat1..3` runs the same models, prompt,
settings and row states into `traces/<variant>/`.
- The paraphrase run refuses unless the paraphrase pool matches its committed manifest. Its rows
  record the paraphrase's hash and the original's.
- A repeat needs a named model.
- A variant's dry run estimates from each model's own main-run bills on the same items:
  - a repeat of `gemma-4-26b-a4b` costs $0.093, so the three repeats cost $0.28, against the plan's
    $1.40 on pricing-document rates;
  - the paraphrase run over the eleven models costs about $49.2. Kimi K3 is $19.6, Claude Sonnet 5
    $11.8 and GLM-5.3 $10.9.

**The writer** is Mistral Large 3 (`mistralai/mistral-large-2512`), from a family on neither the
roster nor the judge's side. Only Mistral serves it, with no declared quantization, so it is routed
as the closed-weight models are. It runs at temperature 0.7, with the prompt hashed into every row.
- **Checks.** An attempt passes only if all of these hold:
  - the numbers are the same, as a multiset;
  - every technical token is kept;
  - the part labels are in the same order;
  - word similarity is at most 0.75;
  - the length is 0.7 to 1.5 times the original's;
  - there is no preamble.

  The writer gets up to three attempts per item.
- **Dry run.** $0.15 if every item passes at the first attempt, $0.45 at most.
- **Self-test.** Constructed cases give the verdicts they should. Each of the 450 originals, checked
  against itself, fails only the near-copy check.

**The expert check.** One own-branch expert judges each paraphrase.
- **Assignment.** An expert gets whole templates, 10 templates and 30 items each.
- **Questions.** Three: same problem? same answer? a new ambiguity, error or hint? A pair is kept on
  yes, yes, no, and any other answer needs a note.
- **The app.** Plain text in read-only boxes, not rendered Markdown, which was round 1's trap
  (D-108). It records timestamps. Kits carry only the expert's id and opaque codes.
- **Self-tests.**
  - The kits are built and scored on a stand-in pool: every paraphrase is assigned once, to its own
    branch, and no item id is revealed.
  - The app is driven headless with Streamlit's AppTest. Submitting is blocked until all three
    questions are answered, a rejection is blocked without a note, and each row is saved.
- **Q5.** `analyze.py` drops a rejected pair from both arms once the check returns, and labels Q5
  provisional until then.

**Timing.** The plan runs the paraphrase arm in the same week as the main run (28–29 September), so
that the served models match: by about 5 October.

## D-142 — E5 and the step router built for the full run, each reproducing the pilot exactly; nothing billed

**Date:** 2026-09-29 · **Status:** DECIDED (the runs are the owner's to start, by `EVALUATION_GUIDE.md`) · **Evidence:** `full_run_28092026/judge.py`, `router.py`, `judge_calls.py`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, `EVALUATION_GUIDE.md`

Both stages read score.py's store and write beside it: `scores/<variant>/e5/` and
`scores/<variant>/router/`. Their replies are kept in stores keyed by the SHA-256 of (model, settings,
prompt), so no reply is bought twice. A run refuses to start without `--yes` and stops starting
calls at a spend cap. Each call has a wall-clock deadline of 900 s.
- **E5** is the pilot's (D-105): the prompt and parser are `e5_hybrid`'s and the call is E1's, with
  JSON mode, temperature 0, 16,384 tokens and three attempts. It goes out at OpenRouter's default
  routing, as the pilot's did.
  - **One call per trace** whose item has milestones E3 did not reach, holding the question, the
    trace and those milestones.
  - **The score is E5-strict.** E5-lenient is kept as a diagnostic only (RESULTS_E5).
- **The router** is D-113's design.
  - **Rule C:** every step the shipped digit rule does not flag, a trace's steps in one call.
  - **The prompt** is the framework's Tribunal prompt, evaluated from the f-string in
    `_tier2_tribunal_batch` as it stands in the framework's source. The framework's model libraries
    are not needed.
  - **The call:** 8,192 tokens at the provider's defaults, three attempts, and a retry at 16,384
    when a reply is cut off.
  - **A step is flagged** by the digit rule, or when the judge's category contains "error".

**Validated on the pilot**, free: each stage replays the pilot's 300 labelled traces through the
functions the full run uses, reading the pilot's stored replies.
- **E5:** E3's flags equal the pilot's on 300 of 300 traces. The same 142 traces are sent, and every
  prompt is found under the pilot's key. Every milestone's source and every E5-strict score equal the
  pilot's on 300 of 300. The published totals reproduce: 951 of 1,245 milestones by E3; 73 REACHED,
  36 NOT_NEEDED and 185 MISSING; the per-model scores from 0.962 to 0.321.
- **The router:** all 299 prompts are byte-identical to those the pilot sent. The steps under review
  are the same. The published figures reproduce: 8 steps unjudged; precision and recall 0.707 and
  0.603 over all steps, and 0.703 and 0.360 inside correct-answer traces; AUROC 0.675 (0.622 to 0.730)
  on 228 traces.

**Priced on the full run** (`--dry-run`, the eleven models; tokens and call times are the pilot's,
prices OpenRouter's public list):

| | calls | input tokens | at Xiaomi's endpoint | across the six endpoints | time |
|---|---:|---:|---:|---:|---|
| E5 | 8,032 | 12.5 M | $26.56 | $18.59 to $49.71 | about 7.3 h at 16 workers |
| router | 24,506, 7.0 steps each | 44.8 M | $103.32 (output per step), $98.41 (per call) | $72.32 to $194.94 | about 33 h at 8 workers |

E5 comes in well under D-110's $76.82 because only 32% of the full run's traces leave a milestone to
judge, against 47% of the pilot's. The router's estimate agrees with D-113's $114.

`analyze.py` fills Q3's E5 and router columns when the stages have written their rows, and adds the
table of what points at a wrong answer (digit rule, E5 MISSING, the router's judge). Its self-test
checks these on a case worked by hand. The evaluator code gets its tag at the commit the paid runs
start from.

## D-143 — Two independent reviews of the full run's evaluation machinery: what they found, and every item acted on

**Date:** 2026-09-30 · **Status:** DECIDED (the owner, 2026-09-29: "act on all of them, all the action items"; nothing paid runs without notice) · **Evidence:** D-144 to D-149; `full_run_28092026/BOUNDARY_AUDIT.md`, `results/RESULTS.md` (regenerated), `EVALUATION_GUIDE.md`, `ANALYSIS_PLAN.md` (dated notes)

On 2026-09-29 two reviewers with no shared context, one on each of two models, each read every script and
document under `full_run_28092026/` and the evaluators it imports, ran every free check, and then assessed a
list of nine proposed actions from a first review. Their reports were compared and every finding re-verified
on the store before being acted on. What they found, in order of severity, and where each is closed:

1. *Would change a reported number or claim.* The complexity cliff's "7 of 11" moved with the seed and rested on
   a test that is liberal under unequal variances (D-146). glm-5.3's cliff was an output-ceiling effect (D-146).
   The digit rule's wrong-answer rate counted empty traces (D-149). A quantity the question asks for, stated in
   the body and left off the Answer line, scores partial, four times as often for Qwen3-235B-2507 as for any
   other model (D-147). The plan's "only multipart items can be partial" was false (D-145). A partly failed
   judge stage would have printed an understated E5 number, and a failed router call was recorded as a clean
   trace (D-148). The answer check's last-digit windows were decided at their edge by binary rounding (D-147).
   One of the 32 pairwise claims rested on partial credit and the seed (D-146).
2. *Would fail a reproducibility check.* The score store's recorded commit did not contain the code that scored
   it, and its hashes were line-ending dependent (D-144).
3. *Robustness of the paid stages.* `--max-usd` reset on every resume; a call past the deadline was bought
   again; a failing prompt was retried every run; stage rows were not tied to the store they read; a bare
   variant run would have billed two models outside the roster; seven rows with a provider fault inside a
   200 were scored (D-148).
4. *Completeness.* The per-template SD the July rebuttal promised, the first-flagged-step position, failure
   against milestone count, the 1% rate, the unusable split, E3's chance floor per model, the endpoint
   table matched on template, a noise floor for Kendall's tau, Q5's coverage delta, a committed per-template
   table, and the results' own provenance (D-149).

Of the first review's nine proposals, the reviewers found seven correct or correct-but-incomplete and two in
need of correction, and both corrections were taken: `arith.py` is not edited for a console-encoding
cosmetic, because the edit would change a recorded evaluator hash (the guide sets `PYTHONIOENCODING`
instead); and the paraphrase inference runs beside E5, not after it, since the same-week rule binds the
inference and the two use different providers. One reviewer suggestion was not taken: pinning the
paraphrase arm to each model's majority endpoint would break the plan's "same settings as the main run";
the arms keep the same routing, each row records its endpoint, and Q5 reports the same-endpoint pairs beside
the whole (D-149).

Everything free was done and is recorded below; the paid steps (E5, the router, the paraphrase arm, the
repeats) still wait for the owner's approval, each with its dry-run estimate.

## D-144 — Provenance: the store's recorded commit was wrong and its hashes line-ending dependent; fixed, and the store re-scored from a clean, recorded commit

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `full_run_28092026/score.py` (`provenance`, `--replace`), `scores/main/CONFIG.json` (local), `results/results.json` (`provenance`), the archived store `scores/_replaced/main_20260929T194154Z/` (local)

`scores/main/CONFIG.json` named commit `b73c711`, but the `answer.py` that scored the store carried the
D-138 and D-139 rules committed 25 minutes later as `1363e80`: the store was scored from a dirty tree, and a
reader checking out the recorded commit would have got an answer check that gives gpt-oss-20b 0.826 where
RESULTS.md printed 0.799. `results/main_pre_d137/` named `9915400`, a commit at which `score.py` did not yet
exist. The recorded evaluator hashes were over working-tree bytes under `core.autocrlf=true`, so `arith.py`,
whose code had not changed, hashed differently from its blob, and no other checkout could match all seven.
The rows themselves were right: re-scored at HEAD they were byte-identical (both reviewers).

**What changed.** `score.py` now records, for each evaluator file, its SHA-256 over LF-normalised bytes and
its blob at HEAD; the commit, its tag if any, and whether any evaluator or input file was dirty; the hashes
of the manifest and `diversity.json`; and per model the trace file's hash and when it was scored. CONFIG is
written after each model, never before. A store scored with other code or inputs is refused; `--replace`
archives it under `scores/_replaced/` first, and a model already scored with the same code on the same
traces is skipped. Rows are written atomically. The milestone cache carries a sidecar naming the manifest
and the `milestones.py` it was built from. `judge.py` and `router.py` record their commit, dirty flag, tag,
reply-store summary and the digest of the store CONFIG they read; `analyze.py` refuses a stage built from
another store, and `results.json` records its own commit and LF hash, the store CONFIG and every stage
CONFIG.

**Re-scored** from the clean commit `bf4a43b` (dirty false): 12 models, 27,000 rows. Against the archived
store, every old field is identical except the 81 verdicts D-147 changes (and the half- and double-tolerance
sensitivity labels, which move on 134 and 31 rows under the same rule); four fields are new. The tag
`full-run-evaluation` marks the commit the paid stages start from; its evaluator files are those the store's
CONFIG names. `results/main_pre_d137/` is regenerated at the same commit with a note saying what its store's
CONFIG cannot: the evaluators were `answer.py` and `milestones.py` at `3a7f247`, the rest as at `a72399b`.

## D-145 — Correction to ANALYSIS_PLAN: partial credit applies to every item with more than one target, not to multipart items only

**Date:** 2026-09-30 · **Status:** RECORDED · **Evidence:** `analyze.py`, the store (local); ANALYSIS_PLAN.md carries a dated note

The plan said "only multipart items can be partial, 480 of 2,250". `answer.verdict` gives `partial` whenever
some but not all of an item's targets match, and 531 items carry more than one target: 352 multipart, 83
symbolic, 51 vector and 45 classification. In the store 278 of the roster's 514 partial verdicts fall on
non-multipart items (symbolic 180, vector 57, classification 41). Seven multipart templates have a single
target on every item (`aoq_ati_rectifying`, `arl_beta_mean_shift`, `chase_vs_level_aggregate`,
`coaxial_capacitance`, `poissons_ratio`, `server_configuration_selection`, `signal_energy_power`), and 84
multipart items in 8 templates are checked on fewer quantities than the question's labelled parts, so
"correct" means every quantity the check verifies. The scoring is what the experts validated three-way and is
unchanged; the paper must describe it as it is and list the templates scored on only some of their parts.

## D-146 — Correction to D-136 and D-140: the cliff count rests on a test valid under unequal variances, and it is 3 of 11, not 7

**Date:** 2026-09-30 · **Status:** DECIDED · **Supersedes in part:** D-136 (Q2's test), D-140 (the cliff count) · **Evidence:** `analyze.py` (`welch`, `q2`, `--selftest`), `results/RESULTS.md` Q2

D-136 set Q2's p-value as a permutation of the tier labels on the raw Easy-minus-Advanced difference, and
D-140 reported the gap holding for 7 of 11 models. The second reviewer found, and the re-check confirmed,
that the count moved with the seed (7, 7 and 5 of 11 at 10,000 draws) and that the test is liberal here:
Advanced template means spread two to four times as widely as Easy ones (deepseek-v4.1-flash 0.063 against
0.209, glm-5.3 0.082 against 0.299, claude-sonnet-5 0.102 against 0.217), and a raw-difference permutation
over-rejects when the smaller group has the larger variance; the self-test now shows it on synthetic nulls.

**Adopted.** The count rests on Welch's t-test, Holm across the eleven models: **3 of 11** (gpt-oss-20b,
gpt-5.4-mini, gemini-3.1-flash-lite). The planned permutation is printed beside it and gives 6 of 11 at
100,000 draws. Every resampling test now draws 100,000 permutations, so the Holm floor over 55 pairs is
0.0006 rather than 0.0055. The detectable gap is given as the plan defined it and at the strictest Holm step
from the Welch standard error (gpt-oss-20b 0.156 becomes 0.221, deepseek 0.082 becomes 0.135). The plan's
within-tier bootstrap intervals are unchanged and remain what the paper reports per model.

**Two variations, added.** With unusable rows left out of the template means, glm-5.3's gap falls from +0.133
to +0.009, glm-5.3-flash's from +0.072 to +0.036 and muse's from +0.050 to +0.029: 145 of the roster's 246
unusable rows sit on Advanced templates, so those cliffs measure finishing within the output ceiling. Without
the nine symbolic templates the count is unchanged. On the pre-D137 record the corrected count is 0 of 11
against 2 under the planned test.

**Q1.** The template-level test on the fully-solved rate now stands beside McNemar's item-level one. The pair
gpt-oss-20b against qwen3-235b-a22b-2507, the 32nd claim, holds on the answer score at 0.0530 after D-147 and
not on the fully-solved rate (0.51), so **31 of 55** pairs hold; the two template-level tests agree on all 55.

## D-147 — The answer check's last-digit windows: the inclusive boundary adopted; the half-unit and whole-trace readings reported as sensitivities

**Date:** 2026-09-30 · **Status:** DECIDED (the owner may reverse it: one commit) · **Evidence:** `full_run_28092026/boundary_audit.py`, `BOUNDARY_AUDIT.md`; the self-test in `evaluators/answer.py`; `SCORER_VALIDATION.md`

`answer.match` accepts a value within one unit of its own last digit, or of the gold's. A value exactly one
unit off sits on the edge, and in binary the edge is not exact: 0.063 - 0.062 is 0.0010000000000000009, so
whether such a value passed depended on how the difference rounded. The windows are documented as inclusive
("within one unit"), so the comparison now carries a relative slack of 1e-9 (`answer.SLACK`) and the verdict
follows the rule. Measured as D-138 was, against the check before the change:

| reading | gold | pilot verdicts moved | full run, 12 models |
|---|---|---|---|
| inclusive, adopted | 2,250 of 2,250 | 0; agreement 0.982 / 0.927 unchanged | 81 change: 67 incorrect to correct, 14 partial to correct; scores +0.001 to +0.008 (gpt-5.4-mini 0.829 to 0.837, gpt-oss-20b 0.799 to 0.805); 27 templates, led by `hydraulic_jump_energy_loss` 9 and `rational_method_peak_flow` 8 |
| half-unit, a sensitivity | 2,250 | 0 | 275 change, 201 correct to incorrect; scores -0.002 to -0.026 |
| whole trace, a sensitivity | 2,250 | 2: 1 to the experts, 1 away; non-partial agreement 0.986 | 425 against the check before, 358 of them partial to correct; qwen3-235b-a22b-2507 0.873 to 0.893, the others +0.002 to +0.011 |

The half-unit reading requires a correct rounding at the precision shown; the experts validated the one-unit
rule and cannot arbitrate (no pilot verdict moves), so it is reported, not scored. The whole-trace reading
credits a numeric part the Answer line leaves out when the trace states it anywhere, for a trace whose
Answer line already matches at least one part: it answers the question the format-only partials raised (92
of Qwen3-235B-2507's 115 partial verdicts, 6 to 48 for the other models), and the pilot holds only two such
traces, which the experts split one each way. It is reported, not scored, because the prompt asks for the
final result on the Answer line and the experts validated the check there. A first version of that reading
re-read every part from the whole trace; it credited wrong Answer lines whose working held the right number,
turned 1,089 incorrect verdicts correct and moved 27 pilot verdicts away from the experts, and was rejected
before it was adopted. E5 and the router depend on E3 and the digit rule, not on the answer check, so none of
this touches a judge prompt.

## D-148 — Guards for the paid stages, the harness and the roster file

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `full_run_28092026/judge_calls.py`, `judge.py`, `router.py`, `analyze.py` (`stage_q3`, `--selftest`), `run_traces.py`, `models.json`

- **An unanswered call leaves a rate.** `judge.summarise` gave a trace whose call got no reply an E5 score
  equal to E3's fraction, and `analyze` counted it; `router.summarise` recorded a failed call as 0 unjudged
  and no flags, a clean trace. Now the router leaves every sent step unjudged, and `analyze` computes every
  judged rate over the answered calls, prints the count without a reply beside it, and marks the stage
  incomplete in the header. The self-test covers both.
- **The cap is cumulative.** `judge_calls.run` started its spend count at zero on every invocation, so two
  resumes at `--max-usd 50` could bill $100; the guide read as a cumulative cap. The count now starts from
  what the reply store records over every line, replies and failures alike.
- **Late replies are kept.** A worker stores its own result on return, so a call that outlives the deadline
  is stored when it arrives and read on the next run instead of bought again. The client makes no SDK
  retries and its timeout is 300 s, so one fetch of three attempts stays inside the deadline (960 s); the
  pilot's 600 s with two SDK retries could outlive it. A key that has failed in three runs is left alone and
  counted. Jobs are interleaved across models, so a stop at the cap leaves every model partly judged.
- **Stages are tied to the store.** Both stages refuse to start if a store row's recorded trace hash no
  longer matches the trace on disk, record the digest of the store CONFIG they read, and get `--status`.
- **The roster file.** `qwen3-235b-a22b` (never run; no endpoint met the routing rule) and the set-aside
  `qwen3.8-27b` carry `"run": false` with the reason, and no mode calls an inert entry unless `--model` names
  it; a bare `--variant paraphrase --yes` would otherwise have billed the set-aside model about $11.
- **A provider fault inside a 200** (`finish_reason` `error`) is a service failure and is retried. Seven such
  rows in the main run (4 correct, 2 incorrect, 1 unusable) are kept and reported as a count.

## D-149 — Reporting additions to the analysis, made after the first results were read and labelled as such

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `analyze.py` (docstring, `--selftest`), `score.py` (the new row fields), `results/RESULTS.md`, `results/per_template.csv`

Each of these answers a reviewer's ask (NEXT_CYCLE_REVIEW sections 3.1, 3.3 and 9.3, the July rebuttal) or
one of the two reviews, reads the existing store, and is labelled in RESULTS.md as added after the data.
Where a number is quoted, it is from the regenerated results.

- Beside every headline score, the per-template SD the July rebuttal promised, in both senses: within a
  template over its 15 items (0.041 for deepseek-v4.1-flash to 0.213 for gpt-oss-20b) and between template
  means (0.117 to 0.276).
- The digit rule's wrong-answer rate over the answered wrong answers (glm-5.3 0.111 where the empties gave
  0.016); the 1% rate beside it; the position of the first flagged step in a fully solved trace; and the
  wrong-answer rate against the item's milestone count.
- E3's chance floor per model, the trace scored against a sibling item's milestones (0.107 to 0.196 on the
  readable wrong answers), stored per row so E5's coverage can be read against it.
- "Unusable" split into empty and unreadable (243 and 3 on the roster); the answered rows cut at the output
  cap and scored (18; 8 of them correct); the seven odd finish reasons.
- The score by serving endpoint matched on template, because dispatch order confounds a raw per-endpoint
  mean with the templates each endpoint served: every endpoint with more than 36 templates in common with
  the others differs from them by at most 0.016; gpt-oss-20b's two minor endpoints served 8 rows in all and
  their differences are noise.
- A noise floor for Kendall's tau, the ordering on one half of each template's items against the other:
  median 0.881 over 200 splits, because the top five models lie within 0.012 of each other.
- Q5: the paired E3 coverage difference, the E5-strict difference once E5 has run on both arms (the coverage
  delta section 6 asks for; without it the $5 for E5 on the paraphrase arm would feed nothing), and the
  same-endpoint pairs alone.
- The sensitivity rows of D-147 and the pool without the nine symbolic templates (D-138).
- `results/per_template.csv`: one row per template and model, aggregates only, so figures and a reader's own
  template bootstrap can be redone without the private store.
- The results carry their provenance (D-144), the header says which stages have run and whether any is
  incomplete, and the caption states that the Holm floor is 0.0006 and that a rate over a subset of traces
  resamples the templates with a qualifying trace.
- The Wilcoxon test on milestone coverage the July rebuttal promised is not added: it was tied to the
  decoupling claim the paper no longer makes (NEXT_CYCLE_REVIEW 9.3 item 13); if that claim returns, the
  test returns with it.

## D-150 — Stream 1 readiness: the judged-stage paths dry-tested end to end, a sample of the digit rule's flags drawn, the variant review, and one runbook per stream

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `full_run_28092026/EVALUATION_GUIDE.md` (stream 1), `PARAPHRASE_RUNBOOK.md` (stream 2), `flag_sample.py`, `trace_review.py --variant`, `paraphrase_kit.py`, `analyze.py`

The owner will run the two remaining streams in separate sessions: the evaluation of the eleven models'
traces (E5, the router, the repeats, the flag reading, the analysis) and the paraphrase experiment. Before
that, a last check that the machinery is in place for both, done by exercising the paths rather than by
reading them:

- **The judged-stage paths, free.** `judge --score` and `router --score` were run against empty reply
  stores, which is the "every call unanswered" case: both wrote their rows and a CONFIG carrying the commit
  (`7ebc348`, tag `full-run-evaluation`, dirty false) and the store's digest; `analyze` marked every model
  incomplete in its header and printed the count without a reply beside each E5 rate; a store CONFIG
  rewritten under the stage made `analyze` refuse with the message that names `judge --score`; the dry rows
  were then removed and the committed results regenerated identically, apart from the provenance line, which
  names the commit the analysis ran at (the parent of the commit holding the results). No reply store was
  created: nothing was called.
- **The digit rule's flags.** `flag_sample.py --draw` takes up to 20 flagged claims per model, round-robin
  across the model's templates, recomputed from the traces at the store's step index: 220 claims from
  gpt-oss-20b 607, gemma-4-26b-a4b 796, deepseek-v4.1-flash 24, qwen3-235b-a22b-2507 1,547, glm-5.3-flash
  185, glm-5.3 191, muse-glimmer-30b 464, kimi-k3 51, gpt-5.4-mini 357, gemini-3.1-flash-lite 416 and
  claude-sonnet-5 375 flags in all. The sample and its CSV stay local under `scores/flag_review/`; an
  author fills the verdicts (slip, checker, unsure) and `--score` writes `FLAG_REVIEW.md`, counts only,
  with a per-model precision and its Wilson interval. A checker verdict with a pattern is fixed the D-137
  way. The reading itself is still to be done.
- **The variant review.** `trace_review.py --variant <name>` applies the main run's integrity checks to a
  variant's traces: a paraphrase row against `paraphrase/manifest.jsonl` and the original's hash, a repeat's
  against the pool manifest, and each model's served id against the main run's ("as main"), since the arms
  must be answered by the same checkpoint. `judge --status` and `router --status` print the billed total
  and the serving-provider mix per model, for the record the paper needs of the judge's endpoints.
- **The experts' check, partial returns.** `accepted.json` now lists every assigned pair; one not yet
  returned is kept provisionally and counted as outstanding, and only a rejected pair leaves both arms.
  The first version dropped every unreturned pair as if rejected.
- **The runbooks.** `EVALUATION_GUIDE.md` is rewritten as the stream 1 runbook and `PARAPHRASE_RUNBOOK.md`
  written for stream 2: session-start checks with their expected output, each step's command, cost, time
  and finish condition, what to commit and what stays local, what to record in DECISIONS, and what "done"
  means. Every self-test passes after the changes.

Nothing paid has run. The costs stand as D-142 and D-143 give them.

## D-151 — The decoding repeats over four cheap models, not one, and their run

**Date:** 2026-09-30 · **Status:** DECIDED (the owner, before the runs) · DONE (the runs) · **Evidence:** `run_traces.py --variant repeatN --model … --dry-run` and `--status`; `trace_review.py --variant repeatN` (`TRACE_REVIEW_repeat1.md` to `_repeat3.md`); `score.py --variant repeatN`; `analyze.py` (`results/RESULTS.md`, Decoding repeats)

The plan asks for one cheap model (ANALYSIS_PLAN, "Also reported, not tested"). One model's spread says
nothing about another's, so the owner extended the repeats to the four cheapest models of the roster by
the variant dry runs: `gemma-4-26b-a4b`, `gpt-oss-20b`, `qwen3-235b-a22b-2507` and
`gemini-3.1-flash-lite`, from three developers. The dry runs priced a repeat at $0.674 for the four,
$2.02 for three, against $0.28 for Gemma alone and $97.53 for all eleven (Kimi K3 $13.05, Claude Sonnet 5
$7.96 and GLM-5.3 $6.94 a repeat). The other seven models have no decoding figure, and the paper says so.
`analyze.py` needed no change: it prints a row for every model with repeat rows.

- **The runs.** 2026-09-29, 22:13 to 23:46 UTC, at `047c79f`: one process per model, `repeat1` to
  `repeat3` in turn, 15 workers, the main run's settings. Every variant has all 300 items for every model.
  gpt-oss-20b left 5, 7 and 6 empty; the others none. The rows record $0.625, $0.615 and $0.623, $1.863
  in all. The OpenRouter key's usage rose $1.471 from a read at 22:11:51 UTC, before the launch, to one at
  22:23:59, by when Gemini and Gemma had finished and the other two had not. From 22:24 the repeats
  shared the account with E5 (D-152), so from there on the rows carry the per-stage figure.
- **The review.** Every variant: the pilot's prompt, one request set and one served model per key, the
  same served model as the main run for all four ("as main" yes), no row outside the frozen items, none
  malformed. The first run of the review counted a repeat against the whole pool, 1,950 "missing" per
  model; `variant_manifest` now expects a repeat's 300 items (`fb0ea02`). The main run's review
  regenerates unchanged.
- **The scores.** Scored at `1c16b6b`, the store's CONFIG clean. `RESULTS.md` prints, per model, the three
  repeat scores on the 300 items, the main run's on the same items, the SD and range over the three
  repeats, and the share of items with the same verdict in all three: SDs of 0.003 to 0.016, ranges of
  0.005 to 0.032, and 0.850 to 0.937 of items with the same verdict every time.
- **The backup.** The twelve repeat trace files are archived locally, `full_run_traces_repeats_2026-09-30.zip`
  with its checksum, every member checked against its source by `backup_archive.py`. PowerShell's
  `Compress-Archive` had refused one of them as held open by another process, which Python's zipfile does
  not. A private Kaggle copy was made the same morning in manual mode: Kaggle unpacked the zip into its
  twelve files, and a download-back matched all thirteen, the checksum file with them, byte for byte.

## D-152 — E5's first night: calls that never return, the harness fix, and the cost basis

**Date:** 2026-09-30 · **Status:** DECIDED (the fix) · OPEN (E5 running) · **Evidence:** `judge_calls.py --selftest`; `judge.py --status`; `judge.py --dry-run`; `router.py --status`; `scores/e5.log` (local)

The owner approved E5 over the eleven models on the dry run's $26.56 at Xiaomi's prices ($18.59 to
$49.71 across the six endpoints), with a $55 cumulative cap. It was launched 2026-09-29 at 22:24 UTC at
`efa8fe9`, 16 workers, through `keepawake_run.ps1` started from WMI, so it holds a keep-awake request and
outlives the session that launched it. The laptop was on mains power.

- **Calls that never return.** Nine calls of the first run had no reply at the 960 s deadline and never
  returned. Each kept its thread, so the fixed pool of 16 lost a worker each time, and within half an hour
  the run was returning calls at less than half its first rate. The client's 300 s timeout
  is between bytes, not over the call, and OpenRouter keeps a waiting connection open, so it never fires
  (D-148 assumed it bounds a fetch). A restart at 22:58 UTC freed the workers; the nine prompts were asked
  again and answered normally. Four more calls stuck in the next 350, about one in a hundred.
- **The fix** (`1c16b6b`). A call past the deadline gives back its worker: the pool holds the `workers`
  live calls plus room for `STUCK_THREADS` (256) abandoned ones, and a new call starts whenever a live one
  returns or passes the deadline. Its reply is still stored if it arrives. What is sent and how a reply is
  read do not change, and `judge --validate` and `router --validate` reproduce the pilot as before.
  `judge_calls.py --selftest` runs `run()` offline with simulated hangs, eight checks. Run against the old
  code, it fails the two hang checks, because three calls hanging on three workers stop the run until they
  return: the self-test took 62 s there against 4.5 s on the fix. E5 was restarted on the fix at 23:47 UTC. A restart drops the calls in
  flight, at most 16, and they are asked again; anything they billed shows in the account, not in the rows.
- **The cost basis.** The dry runs take MiMo's output per call from the pilot's replies, and
  `judge --status` and `router --status` now print the per-reply figures beside that basis. E5's first 508
  replies average 4,209 completion tokens (median 3,122.5) against the pilot's 3,024 (2,469), 1.39 times;
  the median call takes 70.8 s against 52.6. The mean billed per reply is $0.00501. Over 8,032 calls that
  is about $40, against the $26.56 estimate and inside the cap. The first replies ran longer than these, and
  the calls are taken one model at a time in turn, each model's in the order of its score store, so the
  figure moves as the run reaches other templates.
- **The router, by the same basis.** Its dry run's $98.41 and $103.32 at Xiaomi's prices assume the pilot's
  output: 3,702 completion tokens per call for the first, 562 per step under review for the second. At E5's
  ratio of 1.39 on output they become about $129 and $136. That is an extrapolation, not a measurement:
  only the router's own replies will say how long they run.

**Open.** E5's finish, `without_reply`, the billed total by the rows and by the account, and the providers,
recorded when the run ends.

## D-153 — The flag reading: one reader for all 220, a workbook for a reader outside the code, the model hidden

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `flag_sample.py` (`--reader-copy`, `--merge`, `--score --read-by`); a round trip on copies (below)

The owner may give D-150's reading to someone other than an author. So that it reads the same either way:

- **One reader for all 220, not one per branch.** The sample is balanced by model, not by branch, and by
  branch it falls 65 industrial (27 templates), 51 civil (23), 35 mechanical (20), 35 electrical (20) and
  34 chemical (19). The output is a precision per model, and the models' flags sit unevenly across
  branches: 13 of DeepSeek's 20 are industrial, 7 of gpt-oss-20b's chemical. A reader per branch would tie
  a model's figure to that branch's reader. No branch expertise is needed, because each row is arithmetic
  and rounding at the displayed precision.
- **A workbook for a reader outside the code.** `sample.csv` has no step text, which a reader needs to
  judge a rounding chain, a unit or a wrongly split clause, and it names the model on every row, which can
  bias the reading. `--reader-copy` writes `scores/flag_review/flag_review_reader.xlsx`: the claims in the
  sample's shuffled order, with the step's full text, the units the checker read, a slip / checker / unsure
  list, numbers kept as text, and no model column. The instructions are a separate file, committed because
  it holds no trace text, so the record shows what the reader was told: `FLAG_READER_INSTRUCTIONS.md`,
  copied beside the workbook as `flag_review_instructions.md`. `--merge` writes the
  returned verdicts and notes into `sample.csv` by `code`. It refuses a code or a verdict it does not know
  and reports rows missing, verdicts changed and `checker` verdicts without a note.
- **The encoding.** `--score` and `--merge` read a CSV with or without a byte-order mark, which Excel adds
  when it saves "CSV UTF-8".
- **Who read it.** `--score --read-by` names the reader in `FLAG_REVIEW.md`'s title (the default is "an
  author"), and the paper says the same.

A round trip on copies filled the workbook, merged it, scored it, read a CSV with a byte-order mark,
and refused an unknown code and an unknown verdict: all eight checks passed. The real `sample.csv` was
unchanged, and no `FLAG_REVIEW.md` was written. The workbook and any filled copy stay local, like the
experts' labels.

## D-154 — The flag reading: half of this roster's digit-rule flags are the checker's

**Date:** 2026-09-30 · **Status:** MEASURED · OPEN (what to do about it) · **Evidence:** `flag_sample.py --merge` and `--score --read-by "a domain expert"` (`FLAG_REVIEW.md`); the filled workbook (local)

A domain expert read the 220 sampled flags in the reader's workbook (D-153) and decided every one: 111
slips and 109 checker misreadings, none unsure, with a note on every checker verdict.

- **Precision 0.505** among the decided (Wilson 0.439 to 0.570), against the 0.750 the pilot measured
  (SCORER_VALIDATION.md). By model it runs from 0.100 (glm-5.3, 2 of 20) to 1.000 (gpt-5.4-mini, 20 of
  20); at 20 claims a model, every interval is wide.
- **The notes** give the correct computation each time and name what the checker misread.
  `FLAG_REVIEW.md` groups them by their first words: "split at" 34, "thousands separator" 14 in three
  spellings, "missing brackets" 7, "unit conversion" 6, and single notes for the rest. Read in full, they
  fall into five kinds:
  - where a claim starts and ends: two claims read as one, split at '=>', 'Then', a period or
    `\qquad`, or cut inside a chain of partial products;
  - variables read as units (E I, D S, μ L, d, P, t, V);
  - number formats: thousands separators written as spaces or `\,`, scientific notation in a
    denominator read without its brackets, a trailing '…' that means truncation, ranges written with a
    dash, a percentage of a quantity;
  - units the checker does not convert (h and min, J and kJ, m and mm, ha, ksi and psi), and degrees
    read as radians;
  - relations other than equality ('>', '≈' as a comparison, '≠') read as '='.
  
  None of them is about the models' arithmetic.
- **What follows.**
  1. Q3 cannot print the pilot's 0.750 beside the digit rule's rates. This roster's is 0.505, and by
     model the share of real slips among the flags varies so much that the raw flag rates do not compare
     the models' arithmetic.
  2. The misreadings are the kind D-137 fixed. `arith.py`, the rule's parser, shares no code with the
     answer check or E3: neither imports it. So a fix confined to it moves no answer score, no E3 result
     and no E5 prompt. The re-score, the rebuild of E5's rows from the stored replies and the analysis
     are free, and no reply is bought again.
  3. These 220 would be what a fix is built from, so their precision after it would flatter it. An
     honest figure needs a fresh sample, read after the fix.
  4. The router waits, since the rule's flags are part of its prompts.

## D-155 — E5's run: every call answered, $36.30 in the rows

**Date:** 2026-09-29 to 30 · **Status:** DONE · **Evidence:** `judge.py --status`; `scores/e5.log` and the reply store (local); `analyze.py` (`results/RESULTS.md`, the analysis at `774cfc0`)

- **The run.** The owner approved it on the dry run's $26.56 ($18.59 to $49.71 across the six endpoints),
  with a $55 cumulative cap and 16 workers.
  - It was launched 2026-09-29 at 22:24 UTC at `efa8fe9`.
  - It was restarted twice (D-152): at 22:58 UTC to free the workers held by calls that never returned,
    and at 23:47 UTC on the fix at `1c16b6b`.
  - The main pass ended at 11:00 UTC on the 30th with 19 calls without a reply. The guide's re-run of the
    same command asked those 19 at `774cfc0`, and all 19 were answered.
  - No call is without a reply for any of the eleven models. 53 keys failed at least once and were answered
    later, and none was given up. The watchdog never saw the log fall silent.
- **The bill.** $36.3016 over the store's 8,090 lines: 8,032 replies and 58 failure lines. The mean billed
  per reply is $0.00451. Mean completion tokens were 3,354 (median 2,473), against the pilot's 3,024
  (2,469), and the median call took 56.0 s.
- **By the account.** The key's usage rose $38.325, from $356.587 at 22:11:51 UTC on the 29th, before the
  repeats, to $394.912 at 11:10:56 UTC after the retry. The rows of the two stages hold $38.165, the
  repeats' $1.863 and E5's $36.302. The $0.160 the rows do not show covers the calls in flight at the two
  restarts and any other use of the key in those hours, which cannot be told apart from here.
- **The providers.** At OpenRouter's default routing, `judge --status` shows each model's calls served
  mostly by Xiaomi, DigitalOcean, Novita and StreamLake. AtlasCloud and GMICloud served a few of the first
  replies.

`analyze.py` now fills Q3's E5 columns, and its header names the four deterministic and judged stages
present, with none incomplete. The digit rule's caption still carries the pilot's 0.750, which D-154
replaces.

## D-156 — The digit rule's parser fixed from the expert's notes, adopted with a documented exception, the store re-scored

**Date:** 2026-09-30 · **Status:** DECIDED (the owner: fix, and adopt with the exception) · DONE (re-score, E5's rows, the analysis, round 2 drawn) · **Evidence:** `evaluators/arith.py` (its docstring and self-test); `digit_fix.py` (`DIGIT_FIX.md`); `validate_scorer.py` (`SCORER_VALIDATION.md`); `gold_validation.py` (`GOLD_VALIDATION.md`); `router.py --validate` and `--dry-run`; `judge.py --validate` and `--score`; `score.py --status`; `analyze.py` (`results/RESULTS.md`); `flag_sample.py --round 2`

The owner chose the fix (D-154).

- **The fix.** It is confined to `arith.py`; its rules are in the docstring, under READING THE FULL RUN'S
  FLAGS. 27 new self-test cases pin them: the forms the notes name, and three slips the expert confirmed,
  which must stay flagged. Three rules were changed or dropped on evidence:
  - The unit factor 1e4 (ha and m²) passed a real slip, 4.036e-5 m³/mol "=" 0.4036 cm³/mol, and was left
    out.
  - Scientific notation was at first grouped for any mantissa. That flagged 15 gold claims in one
    template: `(92 - 20)/92 * 10^6` is a conversion to ppm. It is now grouped only for a decimal mantissa.
  - A unit tail holding a function (`2 sqrt(m k)`) is no longer a unit. This bug was older; the parser
    only reached it once it read `\(…\)` lines.
- **Before against after** (`DIGIT_FIX.md`: `arith.py` at `3a7f247` against the fix).
  - **Gold.** Claims read went from 6,917 to 7,137; claims flagged stayed at 0.
  - **The pilot**, against the experts' step labels:
    - All traces: tp/fp/fn went from 99/21/289 to 122/22/266; precision from 0.825 to 0.847; recall
      from 0.255 to 0.314.
    - Correct-answer traces: from 57/19/121 to 76/19/102; precision from 0.750 to 0.800; recall from
      0.320 to 0.427; F1 from 0.449 to 0.557.
    - 32 steps change. 27 move toward the experts' label: 23 flags added to steps they call incorrect
      and 4 removed from steps they call correct. 5 move away.
  - **The 220 read in round 1**, the set the fix was built from. Of the 109 checker verdicts, 97 steps
    are no longer flagged, 11 are still flagged by the same claim, and 1 is flagged by a real false
    equation on another line. Of the 111 slips, 109 are still flagged by the same claim and 2 by another;
    none is lost.
  - **The full run**, 12 models and 26,668 answered traces. Claims read went from 92,148 to 102,682,
    claims flagged from 5,029 to 4,384, and traces with a flag from 2,634 to 2,375. A flag disappears
    from 946 steps and appears on 520.
- **The exception.** D-137's rule allows no pilot step to move away from the experts, and five do. All
  five are DeepSeek-R1 steps on lines the old parser could not read: some sit inside `\(…\)` delimiters,
  and scientific notation was judged at its mantissa's precision. Each is a real wrong digit at the
  precision the trace displays, and the experts' step labels, which judge a step's engineering, accept
  it. For example, √(1.49486×10⁻⁹) is written 3.867×10⁻⁵ for 3.8663×10⁻⁵, and 3.252354/417 is written
  0.00779939 for 0.00779941. The owner adopted the fix with these five as a documented exception
  (`digit_fix.EXCEPTIONS`); any other step that moves away fails the audit again.
- **Left as they are.** Eleven misreadings in the design set remain, because a rule for them would pass
  real slips or change the rule's definition:
  - three running computations written as one chain;
  - three units the model did not state (m→mm, J→mJ, ksi→psi);
  - one each of: ≈ used as a comparison, V as an unknown volume, a carried rounding, "22.5%" of a
    quantity, and the hectare.
- **The published figures stand.**
  - `validate_scorer.py` reproduces the pilot's digit-rule figures with `arith.py` at `3a7f247` and
    prints the fixed rule's figures beside them.
  - `router.py --validate` replays the pilot with the pre-fix flags, so its prompts are rebuilt as sent,
    and it reproduces.
  - `judge.py --validate` is unaffected, and gold is clean.
- **The re-score.**
  - The main store was re-scored with `score --variant main --replace` at `bc7dfa4`; the old store is
    kept in `scores/_replaced/`. The three repeat stores were re-scored the same way, at `ac134b9`, whose
    only change is the analysis caption.
  - `judge --score` rebuilt E5's rows from the stored replies. No call was made, the reply store is
    unchanged at 8,090 lines and $36.3016, and no call is without a reply.
  - The analysis at `ac134b9` changes only the digit rule's figures: its flag rates and intervals, claims
    per trace, first-flag positions, the 1% rule's rates and the attribution column. No answer score, E3
    figure or E5 figure moves, since neither the answer check nor E3 imports `arith.py`.
  - The digit-rule caption now carries the fixed rule's pilot figures and the roster's round-1 reading,
    labelled as an addition.
- **The router, re-priced** on the new flags: 24,506 calls and 172,025 steps under review, $98.41 to
  $103.53 at Xiaomi's prices on the pilot's output basis. D-152's ratio would make it about $129 to $136.
- **Round 2.** The fix was built from round 1's notes, so round 1 cannot measure the fixed rule.
  - `flag_sample --round 2` drew 206 claims from the fixed rule's flags, with seed 1. That is 20 per
    model, except DeepSeek V4.1 Flash: it has only 6 flagged claims outside round 1's steps.
  - No step is shared with round 1.
  - The workbook and its instructions are in `scores/flag_review/round2/`. `FLAG_READER_INSTRUCTIONS.md`
    no longer states a count.

**Open.** Round 2's reading, and the router's run, which the owner set for after the flag reading.

## D-158 — The flag reading, round 2: the fixed rule's precision on this roster is 0.752

**Date:** 2026-09-30 · **Status:** MEASURED · OPEN (a second fix: the owner's call) · **Evidence:** `flag_sample.py --round 2 --merge` and `--score --read-by "the same domain expert as round 1"` (`FLAG_REVIEW_2.md`); the filled workbook (local)

The paraphrase stream's commits use D-157, so this entry is D-158.

- **The reading.** The domain expert who read round 1 read round 2's 206 claims (D-156) and decided
  every one: 155 slips and 51 checker misreadings, none unsure, with a note on every checker verdict.
- **Precision 0.752** among the decided (Wilson 0.689 to 0.806). Round 1 measured 0.505 with the rule
  before the fix, on other steps of the same traces. The pilot's 0.750, per step against its experts'
  labels on its own five models, is a different measure, printed beside it for scale.
- **By model:**
  - gpt-5.4-mini: 1.000 (20 of 20);
  - gpt-oss-20b, gemma-4-26b-a4b, gemini-3.1-flash-lite and claude-sonnet-5: 0.900 each;
  - qwen3-235b-a22b-2507: 0.800; muse-glimmer-30b: 0.750;
  - glm-5.3-flash: 0.550; kimi-k3: 0.500; glm-5.3: 0.450;
  - deepseek-v4.1-flash: 0.333, on its 6 claims.
  
  At 20 claims a model every interval is wide.
- **The notes:**
  - prose read as equations: "accepted if X = 0 or X = 1", "since Cp = 1.14 < 1.67", "using Q = 1090";
  - markdown table cells run together, all six in Kimi's traces;
  - `ẋ(0)` read as a function, four notes;
  - units the rule does not know (Stokes, bar, ppm, µV), and units the model left unstated (m→mm,
    Pa→MPa);
  - a `…` after a correctly rounded number, where the expert suggests accepting rounding as well as
    truncation;
  - `a/2π` read as (a/2)π, and `\mu L` written in LaTeX;
  - single cases: degrees in a conversion, a quadrant, ratio notation, carried roundings, a sign
    convention and a running product.
- **Housekeeping.** The owner moved round 1's files into `scores/flag_review/round1/`. `flag_sample.py`
  and `digit_fix.py` now read them there, and round 1's report regenerates byte-identical.

**Open.** Whether to fix the rule a second time, with a round 3 to measure it and no round 4, or to
report 0.752. Either way the decision comes before the router runs, since the rule's flags are part of
its prompts.

## D-159 — The digit rule's second fix, from round 2's notes: adopted with no pilot step moving away; round 3 drawn

**Date:** 2026-09-30 · **Status:** DECIDED (the owner: a second fix, measured by a round 3, no round 4) · DONE (the fix, the re-score, E5's rows, the analysis, round 3 drawn) · **Evidence:** `evaluators/arith.py` (THE SECOND READING, self-test); `digit_fix.py --fix 2` (`DIGIT_FIX_2.md`); `validate_scorer.py`; `gold_validation.py`; `router.py --validate` and `--dry-run`; `judge.py --validate` and `--score`; `score.py --status`; `analyze.py`; `flag_sample.py --round 3`

- **The fix.** It is confined to `arith.py`, like D-156's, and its rules are in the docstring under
  THE SECOND READING:
  - table cells;
  - prose connectors, and a comparison's threshold;
  - variables written in any script, and Greek letter names;
  - the unit a side implies when it writes none;
  - pair factors for Stokes, bar, atm and hectares;
  - `...` read as a rounding or a truncation;
  - `a/2π`, a degree sign as a unit, and a sum written a term a line.
  
  28 new self-test cases pin them. They include three slips the rules must go on flagging: a wrong
  power of ten, a zero result, and a claim after a comparison.
- **Five rules were narrowed on evidence while the fix was built:**
  - The Greek-names rule first held `psi`, which turned the unit psi into a variable; it was found
    tracing round 1's `21.8 × 10⁶ psi`. psi stays a unit. With it, a formula written only in Greek
    letters (`σ/ε`) now counts as a formula rather than as an unreadable segment, so it links the chain.
  - The implied unit first allowed any standard factor. That passed a wrong power of ten,
    (338.86×10⁶)/(122×10⁶) written 2.7775×10⁻³ m, and a zero result, 131 − 79 = 0 s through 10⁻¹². It now
    allows only the factor the written unit implies, and no conversion makes a zero.
  - The comparison rule first skipped every right side, which lost the slip W_L ≥ 824/179 = 4.61. It now
    skips only a bare threshold.
  - The zero guard first used a threshold, which flagged 1 mm⁴ = 10⁻¹² m⁴ on a pilot step the experts
    call correct. It now takes only an exact zero.
  - A truncated operand is widened by a whole unit; applied to a truncated result as well, that passed
    the slip 0.03461552… It now applies only inside an expression.
- **Before against after** (`DIGIT_FIX_2.md`: `arith.py` at `bc7dfa4` against `e783962`).
  - **Gold.** Claims read went from 7,137 to 7,167; claims flagged stayed at 0.
  - **The pilot**, against the experts' step labels:
    - All traces: tp/fp/fn went from 122/22/266 to 123/18/265; precision from 0.847 to 0.872; recall
      from 0.314 to 0.317.
    - Correct-answer traces: from 76/19/102 to 76/17/102; precision from 0.800 to 0.817; recall stayed
      at 0.427; F1 is 0.561.
    - Five steps change, all toward the experts: four flags removed from steps they call correct, and
      one added to a step they call incorrect. None moves away, so no exception is needed.
  - **Round 2's 206 flags**, the design set. Of the 51 misreadings, 42 steps are no longer flagged and 9
    still are. Of the 155 slips, 154 are still flagged. The other is `-109 000 × 0.841 = -91.6 kN`. The
    old rule flagged it only because it did not know kN from N. With the unit read, the rule's own
    propagation accepts it, since 0.841 shown to three decimals can be 0.8405, which gives −91.6. The
    expert, computing with 0.841 exactly, calls it a slip. This is a limit of the rule's definition, not
    a misreading.
  - **Round 1's flags**, whose slips must stay flagged. All 111 slips' steps stay flagged. Of the 11
    misreadings D-156 left, 8 remain, the ksi→psi conversion among them, and 1 step is flagged by a real
    false equation on another line.
  - **The full run**, 12 models. Claims read went from 102,682 to 103,834, claims flagged from 4,384 to
    4,141, and traces with a flag from 2,375 to 2,226. A flag disappears from 230 steps and appears on 30.
- **The published figures stand.** `validate_scorer.py` reproduces them all with the published code, and
  `router.py --validate` and `judge.py --validate` reproduce the pilot. The first fix's report is now
  pinned to its commits: `digit_fix.py` measures each fix between two fixed versions of `arith.py`, and
  `DIGIT_FIX.md` regenerates with the same figures.
- **Left as they are:**
  - running computations written as one chain: four;
  - carried roundings: three;
  - a unit no written unit implies, such as ksi→psi, and hours→minutes through a symbolic middle;
  - ratio notation, a sign convention, a quadrant, and `wL` and `Zc` as variables;
  - from round 1: ≈ used as a comparison, V as an unknown volume, and "22.5%" of a quantity.
- **The re-score.**
  - The main store and the three repeat stores were re-scored with `--replace` at `c221a9a`. A first
    re-score at `dfcef12` was superseded by the psi correction. The old stores are kept in
    `scores/_replaced/`.
  - `judge --score` rebuilt E5's rows at `c221a9a`, clean. No call was made, the reply store is
    unchanged at $36.3016, and no call is without a reply.
  - The analysis at `c221a9a` changes only the provenance and Q3's digit-rule fields; no answer score,
    E3 figure or E5 figure moves. Its caption carries the fixed rule's pilot figures and both readings,
    labelled.
- **The router, re-priced:** $98.41 to $103.63 at Xiaomi's prices on the pilot's output basis ($72.54 to
  $195.58 across the endpoints), about 33 hours at 8 workers. D-152's ratio would make it about $129 to
  $136.
- **Round 3.** `flag_sample --round 3` drew 181 claims from the second fix's flags, with seed 2. That is
  20 per model for nine models, 1 for DeepSeek V4.1 Flash and none for Kimi K3, whose remaining flags all
  sit on steps rounds 1 and 2 drew. No step is shared with an earlier round. It was drawn once before
  the psi correction and drawn again after it, before anyone read it; the counts are the same. Round 3
  measures the second fix; by the owner's rule, no round 4 follows.

**Open.** Round 3's reading, and then the router, which needs the owner's approval of the estimate above.
*Superseded 2026-09-30: round 3 was redrawn after D-160, before anyone read it.*

## D-160 — Two review agents checked both digit-rule fixes: their corrections adopted, the LaTeX control space read, round 3 redrawn

**Date:** 2026-09-30 · **Status:** DECIDED (the owner: two independent agents check the fixes before round 3, which stays the last; of the two older gaps they found, the control space is closed and `≈`-only lines are left) · DONE (the review, the fix, the re-score, E5's rows, the analysis, round 3 redrawn) · **Evidence:** `digit_fix.py --review-sample` and `--fix 3` (`DIGIT_FIX_3.md`); `gap_check.py` (`GAP_CHECK.md`); `evaluators/arith.py` (THE REVIEW, self-test); `validate_scorer.py`; `gold_validation.py`; `score.py`; `judge.py --score`; `analyze.py`; `flag_sample.py --round 3`

- **The review.** `digit_fix --review-sample` drew 80 full-run steps whose flag the two fixes changed
  between `3a7f247` and `e783962`: 40 of the 924 that lost a flag and 40 of the 339 that gained one,
  round-robin across ten models (seed 7). None is on a step an expert round drew, and none is from the
  set-aside model. Two agents, each on its own and read-only, checked both fixes' code against the
  docstring's rules, hunted regressions with probe lines, ran the self-test and gold, and judged the 80
  steps. Their verdicts are local, beside the sample. No API was called.
  - They agreed on 77 of the 80. Of the removed flags, A called 37 correct and 3 hid slips; B called 39
    correct and 1 hid a slip. Of the added flags, both called 37 real slips; A called 3 new misreadings,
    B called 2 new misreadings and 1 unsure.
  - The three they split on turn on older rules, not on the fixes: the shown-rounding allowance given to
    an exact datum (RV-3a6f2917), a right side `√N` judged at N's precision (RV-6caff853), and a carried
    rounding of an integer (RV-1f8325a8).
- **What they found in the fixes.** B reported nine regressions and two older misreadings the first fix
  exposed. A reported sixteen regressions and five older gaps. Each input was reproduced on the rule as
  it stood before anything changed. The corrections, each with a plant in the self-test, are in the
  docstring under THE REVIEW:
  - `a/Nπ` read both ways. D-159's a/(Nπ) had turned `4/3 π r³` into 4/(3π) r³.
  - The unit a right side implies, read as written (`kN T`, `m s⁻²`, `\mu m`) and to its power (`mm²`,
    `h⁻¹`).
  - The continuation rules:
    - A markdown closer, a rule or a blank line is not an operator. The first fix had hidden gemma's
      `2.6693 × 24.4694 = 65.314`, where the product is 65.3162.
    - A sum is written a term a line only by lines without an `=`.
    - `+(` continues a sum only when it opens a clause. 124 claims had been lost, among them glm-5.3's
      three Lorentz-force products, each off by ×1000.
  - The `...` of a series; an operand cut short, which is widened upward only; an exponent not taken
    for a mantissa.
  - Whole-number ranges only as a clause's last segment.
  - LaTeX letters (`\psi`, `\Pi`, `\ell`, `\hbar`, accents) as variables; T alone as a variable and a
    tesla only inside a unit; `°R`; ångströms.
- **Three rules narrowed on evidence while being built:**
  - Dropping T from the unit letters, as B proposed, hid muse's `47(0.09873) − 22(0.01994) = 4.20063 T m
    s⁻¹` (4.20163). So only T alone is a variable.
  - Reading a range only as a chain's last segment, for A's `Δ = 3.25 − 3.26 = 0.01`, flagged again a
    round-1 claim the expert called the checker's (`2π/18.67 = 0.336–0.337 m ≈ 0.34 m`). Decimal ranges
    stay as D-156 read them, and A's form is left.
  - The thousands fix below, alone, hid `2385.76 × 70 = 166,\!003.2 \ \text{cm}^3·\text{s}`, where the
    product is 167,003.2. The old rule flagged it only through a misreading, and without the control
    space the corrected right side is unread. The control-space fix restores it.
- **Found by listing the flags the review pointed to.** A claim between two written units is a slip only
  if no conversion makes it hold, so every such flag was listed:
  - correct conversions with no factor: poise and Pa·s, inches and millimetres;
  - one thousands spelling, `8,\!753,\!000`, which split its numbers into claims such as `000 = 776`. It
    raised the flag on 64 full-run steps, mostly Qwen3-235B's, none of which any round had drawn.

  Both are read now, and `\ll` and `\gg` are read as comparisons.
- **The two older gaps: the owner's choice.** Both agents measured two gaps older than either fix as far
  larger than the regressions. `gap_check.py` measures each as a one-line switch on the adopted rule
  (`GAP_CHECK.md`); the samples are of the steps each flagged against the candidate rule the owner chose
  from:
  - **A LaTeX control space before a unit** (`\ \text{m}`) left its segment unread. Reading it as a space
    flags 358 more steps. The samples read 18 of 19 decided (mine), 23 of 25 (B) and 12 of 14 (A) as real
    wrong digits, and it changes nothing on gold, the pilot or either round. Closed.
  - **A line whose only relation is `≈`** is never read. Reading it would flag 495 more steps, 16 of 20
    real in my sample. It also brings misreadings of kinds no round has read (`≈ 17120` for 17,122,
    `515 ≈ 514.9`), and it flags two pilot steps the experts call correct, both real wrong digits of
    D-156's exception class. Left, as a recall gap.
- **Before against after** (`DIGIT_FIX_3.md`: `e783962` against `5c83916`).
  - **Gold.** 7,167 claims, 0 flagged, unchanged.
  - **The pilot.** No step changes. Precision and recall are 0.872 and 0.317 over all traces, and 0.817
    and 0.427 inside correct-answer traces.
  - **The agents' 80 steps.** The two steps both called misreadings lose the flag, and the 37 real slips
    keep it. RV-3a6f2917's `√52,407.228 ≈ 229.0` (228.93) is now flagged.
  - **Round 1.** All 111 slips' steps stay flagged. 8 misreadings stay flagged, and 1 step is flagged on
    another line.
  - **Round 2.** 154 of its 155 slips stay flagged, as after D-159, and 9 misreadings stay flagged.
  - **The full run.** Claims read went from 103,834 to 119,573, claims flagged from 4,141 to 4,383, and
    traces with a flag from 2,226 to 2,470. A flag leaves 82 steps and reaches 370. The control space
    accounts for most of the gain: gpt-5.4-mini gains 123 steps, claude-sonnet-5 92 and gpt-oss-20b 81.
    The thousands spelling accounts for most of the loss: Qwen3-235B loses 66.
- **Left as they are:**
  - a threshold followed by its own derivation, and the zero guard (`2 mm = 0 m`);
  - a decimal subtraction written like a range, and a running computation over a count;
  - the shown rounding given to exact data. Both agents called RV-15fc9f0c a hidden slip: 244/125.7 is
    1.94113, written 1.9403…, where 125.7 = 3 × 41.9 is exact;
  - a right side `√N` at N's precision, `e^{…}`, and a temperature's offset;
  - a power of ten between two written units that the generic factors pass;
  - lines whose only relation is `≈`.
- **The re-score.**
  - The main store was re-scored with `--replace` at `eea846e`, and the three repeat stores at
    `937c842`. Only `flag_sample.py`'s report title changed between the two commits. `judge --score`
    rebuilt E5's rows at `937c842` with no call, and every call sent has its reply. The old stores are
    kept in `scores/_replaced/`.
  - The analysis at `937c842` moves only the provenance and Q3's digit-rule fields; no answer score, E3
    figure or E5 figure moves.
  - The paraphrase store, scored by its own stream under D-159's rule, is not re-scored here: Q5 reads no
    digit field, and `score --variant paraphrase --replace` redoes it free.
- **Round 3, redrawn.** The draw made before the review (181 claims, never sent) is kept aside in
  `scores/flag_review/_superseded/round3_before_d160/`. `flag_sample --round 3` drew 190 claims on 188
  steps, with seed 2. That is 20 per model for nine models, plus all that DeepSeek V4.1 Flash (6) and
  Kimi K3 (4) have left. None is on a step rounds 1 and 2 or the agents read; 496 steps were left out.
  Round 3 measures the rule as it now stands, and by the owner's rule no round 4 follows.

**Open.** Round 3's reading; then the router, which needs the owner's approval of its re-priced estimate.

## D-157 — The paraphrases written: four prompts, a throttled writer, and the notation restored by script

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `paraphrase.py` (docstring, `restore`, `--selftest`, `--check`, `--limit`), `PARAPHRASE.md`, `paraphrase/manifest.jsonl`; the set-aside attempt files `paraphrase/attempts_v{1,2,3}_*.jsonl` (local)

Stream 2's first step (`PARAPHRASE_RUNBOOK.md` section 1) was approved at $0.45 and expected to pass most
of the 450 items at the first attempt. It took four prompts and one change to the pipeline, all in one
day, before the writer's output held the items' notation. What was found, what was changed and why, in
order:

- **Prompt 1 (`6091f248`, the D-141 prompt).** 62 real attempts, 13 passed, 10 items exhausted: the
  writer prettified notation, `1/(-r_A)` as `\frac{1}{-r_A}`, `CH4(g)` as `$CH4(g)$` or `CO₂`, `m^3/s` as
  `m³/s`, `*` as `·`, a Unicode minus for `-`. The tokens and numbers checks rejected every one, rightly:
  the two arms must differ in wording only. $0.015.
- **Prompt 2 (`1ae64428`)** forbade reformatting in the abstract. 83 attempts, 29 passed, 14 exhausted:
  acronyms spelt out (`PFR`, `CSTR`), reaction orders turned into words (`order-2.1` as "second-order", a
  real error the numbers check also caught), `=` as "equals", `$` delimiters dropped where the original
  had them, near-copies. $0.022.
- **Prompt 3 (`bd06fbd6`)** named each of those. A 40-item pilot, all chemical, passed 64% of attempts and
  kept 91% of resolved items; the run then reached the civil items and fell to 30%: caret exponents in
  units became superscripts (`kN/m^3` as `kN/m³` 52 times in 138 attempts), `=`, `>` and `%` became words,
  "in m" became "in meters". 225 attempts, 91 passed, 31 exhausted. $0.068. Chemical kept 59 of 68
  resolved items, civil 29 of 51.
- **Prompt 4 (`aa5648ef`)** gives every rule a worked example, and `--limit` now pilots the first items of
  each branch in turn rather than the first N of the list. The pilot, 10 items per branch: 22 of 66
  attempts passed, 11 of 33 resolved items exhausted; chemical 5 of 5 attempts, civil and electrical 38%,
  mechanical 23%, industrial 16%. New habits with each branch: `<=` as "at most", `2D` as
  "two-dimensional", `x0` as `x₀`, `mu_X` as `μ_X`, `EOQ` spelt out, and `kN/m³` still 9 times. $0.019.

**The decision (the owner, on the options laid out with these numbers):** the prompt is not the lever.
Instead:

1. **`restore()`**: before the checks, the original's ASCII notation is put back wherever the writer used a
   Unicode form the original itself does not use (superscripts to `^`, subscript digits to plain or `_`
   digits as the original writes them, `≤ ≥ ≠`, the Unicode minus, `·` and `×`, Greek letters to the names
   the original uses). It is deterministic, the writer's raw text stays in `attempts.jsonl`, the served text
   is what `pool.jsonl` and the manifest carry (flagged `restored`), and the six checks judge it like any
   other text. Words written for symbols ("at most", "equals", "meters", an acronym spelt out) are prose
   and are not mapped back; such an attempt fails as before. Measured on the stored attempts before it
   was adopted: it recovered 12 of 44 failed attempts of the prompt-4 pilot and 37 of 134 of prompt 3's
   run, and turned no pass into a failure (0 of 113).
2. **Em and en dashes are punctuation**, not technical symbols: the tokens check had counted every
   non-ASCII character as a symbol to keep, so a paraphrase that wrote a comma for an em dash failed. 21 of
   the 450 originals contain one; 4 attempts on 2 items were failing for that alone.
3. **The checks themselves are unchanged**, and the copy threshold in particular: 27 copy-only failures on
   10 items were examined; their similarity over prose words alone was 0.67 to 0.76 against 0.39 to 0.62
   for every passing paraphrase, so they are near-copies in prose too, not data-heavy items the measure
   misjudged.
4. **The part-label regex was a false positive**: `\b[Pp]art\s+\w+` read the prose "part of" (28 times in
   the attempts), "part per" and "part rod" as part labels, so a paraphrase that used those words failed
   the parts check against an original with no labels at all; no original in the pool writes a "Part N"
   label. The pattern now takes a number or a single letter after "Part". On the stored attempts it
   regained 12 items (7 of them exhausted) and broke no pass.
5. **Examined and not adopted:** reading all-caps emphasis words (`WHOLE`, `UP`, `DOWN` in the industrial
   templates) case-insensitively would regain 2 of 51 exhausted items, the rest failing on other counts
   too; the writer's corrections of the originals (`-253°c` to `-253°C`, `N.s/m` to `N·s/m`) fail as
   written and are left to the expert check's account, since the paper cannot claim character-for-character
   fidelity and quietly accept a corrected original.
6. **From scratch under prompt 4**, so that every item has three attempts under one prompt: the passes of
   prompts 1 to 3 were set aside with their attempt files, not carried over ($0.06 and about 40 minutes
   redone), and the manifest records the prompt hash per item.

**The writer's endpoint.** Mistral alone serves `mistral-large-2512`, and through OpenRouter it was
throttled all day ("temporarily rate-limited upstream"): the admitted rate swung between 0.7 and 5 calls
a minute with the upstream's state, at two, eight and sixteen workers alike, so no worker count buys
throughput when the endpoint is closed. A refused call bills nothing; `write_one` retries each call
through `RETRY_SLEEPS` (fourteen tries over about three and a half minutes, jittered) before recording a
service failure, and `--until-done` passes again over the deferred items until every item is resolved,
waiting five minutes after a pass that admitted nothing. The run went at sixteen workers.
`resolved()` judges stored rows through the current `restore()` and checks, so the to-do list and
`--status` agree with `--check` after a change to the checks (before it, both read the flags stored at
write time, and one pass re-asked six items that already passed under restore).

**The run under prompt 4 with restore** (`PARAPHRASE.md`): 450 items written in six passes from 12:47 to
15:32 UTC, 817 admitted calls against 416 refused retry cycles. **316 pass** (chemical 71, civil 71,
electrical 64, industrial 49, mechanical 61 of 90 each), 264 at the first attempt, 39 at the second, 13
at the third; 134 have no paraphrase after three; 35 of the passing ones carry restored notation; **29 of
the 150 templates lose all three items**. Word similarity of the passing ones 0.198 to 0.75, median
0.624; the failed attempts by check: tokens 362, numbers 163, copy 75, parts 20, length 13. Billed
$0.312 including the pilot. Step 1 in all, prompts 1 to 4: **$0.416**, against the $0.69 ceiling approved
for prompt 4 ($0.45, then $0.54 and $0.64 for prompts 2 and 3, each approved on its dry run); no cap was
raised without an approval.

**What this means for Q5.** The arm has 316 pairs instead of 450, over 121 templates instead of 150; the
losses are the items the writer turns into words (`=`, `2D`, `%`, units, acronyms) or copies nearly
verbatim, and they fall unevenly: industrial keeps 49 of 90 items, mechanical 61, electrical 64, chemical
and civil 71 each. The template-level interval loses the 29 templates; the paper says so and reports the
survival per branch from `PARAPHRASE.md`, and Q5's ranking stability is read over the 121 templates
both arms share. The experts' check (`paraphrase_kit`) sees the restored text, as the models will, and
adjudicates what the script cannot, including the passing paraphrases with the lowest word similarity.
The alternative not taken today, a different writer family (Cohere, Amazon Nova), stays open if the
experts reject many pairs.

## D-161 — The paraphrase arm run over the eleven models the same day as its writing; scored and analysed, provisional

**Date:** 2026-09-30 · **Status:** DECIDED · **Evidence:** `TRACE_REVIEW_paraphrase.md`, `trace_review_paraphrase.json`, `results/RESULTS.md` (Q5), `traces/paraphrase/*.log` (local)

**The order.** The owner asked whether the experts' check (step 2) had to finish before the inference
(step 3); it does not, and the runbook has the inference first so that the served checkpoints are the
main run's. Step 3 was approved on its dry run at 16:5x UTC: 316 items over 121 templates, $33.05 on the
main run's own per-item bills, under a $40 ceiling; the kits went out in parallel. A pair an expert later
rejects leaves both arms at analysis; no trace is re-run for a verdict.

**The run.** Eleven processes at 15 workers under the keep-awake launcher, 15:52 UTC. Ten complete by
16:23; Kimi K3 ended its first pass at 16:13 with 55 Moonshot service failures ("temporarily
rate-limited upstream", 64 of its 77 failure rows over the passes) and needed three more passes for those
alone, complete at 16:53 with 315 answered and 1 empty. **$34.17** billed in the rows against the $33.05
estimate and the $40 ceiling: Kimi $18.60 against $13.08, about 40% above its main-run per-item cost;
GLM-5.3 $3.03 against $7.45, the same served model routed over four providers at the cheapest eligible
price; the other nine within cents of their estimates.

**The review** (`trace_review --variant paraphrase`): every model 316 items, **"as main" yes for all
eleven**, prompt "pilot", one request set and one served model each, 0 duplicate final rows, 0 items
outside the paraphrase pool, 3,476 final rows, 22 empty over 8 templates, the main run's usual ones
(`vdw_solve_for_volume` 5, `work_isothermal_virial` 5, `adiabatic_flame_temperature` 4); an empty row
scores 0 as in the main run.

**Backup.** `~/EngTrace_private_backup/full_run_paraphrase_2026-09-30.zip`: `paraphrase/` (the attempts
of all four prompts, the pool, the keyfile, the kits) and `traces/paraphrase/`, 64 files checked member
by member, 8.4 MB, sha256 `0d13e18d…0845f7f` beside it. The private Kaggle copy is still to be made
(manual mode, the owner's account), as for the main run's traces (D-126, D-131).

**Scoring and analysis.** `score --variant paraphrase` (5½ minutes) against the original items' gold,
milestones and question; answer means on the arm from 0.815 (`gpt-oss-20b`) to 0.979 (`muse-glimmer-30b`),
unusable rows equal to the empties. `analyze` fills Q5: paraphrase minus original, paired by item over the
316, lies between -0.032 (`qwen3-235b-a22b-2507`) and +0.014 (`muse-glimmer-30b`) for every model, no
difference survives Holm (the smallest, `glm-5.3-flash` -0.027, CI -0.049 to -0.008, p 0.21), and the E3
coverage differences are of the same size. Kendall's tau between the two arms' orderings of the models is
0.550 (95% CI 0.449 to 0.849) against the roster's noise floor of 0.881: the top five models sit within
0.012 of one another, so their order is not stable under any resampling, which is the reading, not
template exploitation. Q5 is **provisional** until the experts' files are all scored
(`paraphrase_kit --score`); the rejected pairs then leave both arms and the analysis is re-run.

**Committed as `f8b63fc` under the number D-160 by mistake; D-160 is the review sample's (`afdc05b`, `3847fb5`).**

**Open.** *(Both closed in D-162: E5 on the arm ran for $4.56, the experts returned in full.)* E5 on the paraphrase arm (runbook section 4, about $3 at 316 items) once approved; the experts'
returns (a later entry when they are in).

## D-162 — The experts' check of the paraphrases returned in full: 39 of 316 pairs rejected, Q5 final; a tau that had ignored the rejections fixed

**Date:** 2026-10-01 · **Status:** DECIDED · **Evidence:** `PARAPHRASE_REVIEW.md`, `paraphrase_kit.py --score`, `analyze.py` (`q5`), `results/RESULTS.md` (Q5); the experts' files and `paraphrase/accepted.json` (local)

**The returns.** All fifteen experts returned their files within the night of 30 September (the last
at 00:53 UTC on 1 October), 316 of 316 pairs, none outstanding; the folder they landed in,
`paraphrase/experts_filled_paraphrases/`, was added to `.gitignore` before scoring. Median 41 seconds per
item (32 to 300). **277 kept, 39 rejected** (12%): chemical 65 / 6, civil 60 / 11, electrical 61 / 3,
industrial 43 / 6, mechanical 48 / 13. The answers behind the rejections: "same problem" no 31, "adds
an ambiguity, an error or a hint" yes 24, "the answer still answers it" no 4; by pattern, 15 pairs
differ in the problem but not in the answer, 12 differ and add something, 8 are the same problem with
something added, 4 change the answer.

**What the notes say** (they stay local; the reading is the author's). The rejections are of one kind:
the writer drops or alters a qualifier that a script cannot weigh. "Solid" rod or shaft dropped
(`statically_indeterminate` ×3, `statically_indeterminate_shaft` ×2, `angle_of_twist`); "uniform flow"
made "steady, uniform" (`manning_trapezoidal_velocity` ×2); "approximately normally distributed" made
"normal" (`newsvendor_normal_demand` ×3); "the magnitude of the force" made "the force", which admits a
signed answer (`truss_method_of_joints` ×3, the 3 of the 4 answer-changing cases); "valued at" made
"priced at" where the holding rate applies to the value (`epq_finite_production` ×3); "moisture content"
made "contains 13.5% moisture", which can be read on the total mass (`phase_relations_degree_of_saturation`
×3, `borrow_pit_fill_volume`); "normal boiling point" made "standard" (`latent_heat_vaporization` ×2);
"submerged in diesel" made "underwater" (`hydrostatic_force_on_plane` ×2); "rigidly connected at B" made
"a rigid composite shaft" (`composite_shafts_series` ×2); "absolute amplitude" made "amplitude"
(`vibration_transmissibility` ×2). Five templates lose all three of their pairs to one expert's consistent
reading, which is the check working as designed. Two `signal_operations` pairs were rejected for a
leftover preamble, "Here's the rephrased version of the problem:", that the script's `clean` check let
through: its regex expects an ASCII apostrophe and the writer wrote a curly one. The experts caught
both; the regex is to be fixed before any further writing run, not now, since a `--check` after the fix
would change the committed manifest under kits already judged (runbook, section 6).

**The arm.** 277 pairs over 115 templates (of 150): 134 items never got a paraphrase (D-157) and 39 lost
theirs here. The paper reports both losses per branch and names the check as the reason the arm is
smaller than the 450 planned, not the models.

**Q5, final** (`analyze`, the rejected pairs dropped from both arms): paraphrase minus original, paired by
item over the 277, lies between -0.031 (`qwen3-235b-a22b-2507`, CI -0.074 to 0.011) and +0.016
(`muse-glimmer-30b`) for every model; no difference survives Holm (the two nearest, `claude-sonnet-5`
-0.018, CI -0.036 to -0.004, p 0.51, and `glm-5.3-flash` -0.020, CI -0.042 to 0.000, p 1.0); the E3
coverage differences are of the same size (`gemma-4-26b-a4b` +0.022, `glm-5.3-flash` -0.031,
`gpt-5.4-mini` -0.027 with intervals off zero, no adjustment applied to that column). The reading: no
model's score depends on the templates' wording to any degree the paired test can see, which answers the
reviewers' template-exploitation objection and the contamination question with the same measurement.

**The tau that ignored the rejections.** `q5` applied the kept set to every per-model row but built the
item set for Kendall's tau between the arms without it, so the provisional and the first final results
both said "over the 316 items". Fixed: the tau's items are intersected with the kept pairs when the check
has returned (the self-test, which passes no check, is unchanged). Over the 277: **tau 0.636, 95% CI
0.457 to 0.871** (0.550 over the 316 before the fix), against the roster's noise floor of 0.881; as before, the top five models lie within 0.012 of one
another, so the ordering between arms is not stable under any resampling, which is the reading, not
exploitation.

**Open.** *(Closed below: E5 on the arm ran for $4.56.)* E5 on the paraphrase arm (runbook section 4, about $3 at 277 kept pairs plus the 39 the judge
would see anyway, since it runs on what is on disk) once approved.

**E5 on the arm (2026-09-30, approved at $10 of new spend).** 1,227 judge calls over the 316 paraphrase traces' milestones, all returned after one re-run for 3 failures; **$4.56** against the $4.06 estimate. The first launch passed `--max-usd 10`, which the stage reads as a cumulative cap over the store (already $36.30 from the main arm), so it made no call; it ran at `--max-usd 46.30`, the same $10, and the runbook now says so. Q5's E5-strict difference, paraphrase minus original over the kept pairs, lies between -0.014 (`gpt-oss-20b`) and +0.018 (`gemma-4-26b-a4b`, CI 0.000 to 0.037) for every model, every other interval spanning zero: the reasoning the judge reads holds under the rewording as the answers do.

**Backup (2026-10-01).** The arm archived again with the experts' returns: `full_run_paraphrase_2026-10-01.zip`, 80 files checked member by member, sha256 `24d9f614…af9293`; uploaded in manual mode as the private Kaggle dataset `ayeshaiq/engtrace-full-run-paraphrase`, downloaded back and every member matched by path and hash (the one extra file is the checksum itself). The 30 September archive D-161 cites stays beside it.

## D-163 — The flag reading, round 3: the corrected rule's precision on this roster is 0.905; no round 4

**Date:** 2026-10-01 · **Status:** DONE (read by the same domain expert as rounds 1 and 2; by the owner's rule, the last round) · **Evidence:** `flag_sample.py --round 3 --merge` and `--score` (`FLAG_REVIEW_3.md`); `analyze.py`

- **The reading.** The same domain expert read round 3's 190 claims. They were drawn after D-160 from
  the rule as it stands, and none is on a step an earlier round or the review agents read.
  - 171 are slips, 18 the checker's and 1 unsure. Precision among the 189 decided is 0.905 (95% Wilson
    0.854 to 0.939). It was 0.505 before the first fix and 0.752 after it.
  - Per model it runs from 0.632 for GLM-5.3-Flash (7 checker verdicts of 19 decided) to 1.000 for
    GPT-5.4-mini and Gemini 3.1 Flash-Lite (20 each), DeepSeek V4.1 Flash (6) and Kimi K3 (4).
- **The 18 misreadings, by the expert's notes:**
  - five letters that are variables, read as units: `wL`, `SE`, and the unknowns `h`, `g` and `A`;
  - three unit conversions the rule lacks: kbaud, √J to √pJ, and hours to minutes through a product;
  - three expressions whose brackets the model left out but meant;
  - two roundings to significant figures that end in zeros (`≈ 2700 K`, `2600 m³`);
  - two cases of prose or list items set side by side;
  - a range written to unequal decimals (`0.0795 – 0.080`), a chain of partial products, and
    `14 → 14 = 0`, a change of zero.

  The unsure claim sits on a step cut off in the middle of a number.
- **No round 4.** The owner set round 3 as the last, so the rule is left as it stands. 0.905 is the
  precision of the rule the results report. The misreadings above are recorded as its known limits,
  with those D-159 and D-160 left.
- **The caption.** Q3's caption carries the reading, labelled. The analysis was regenerated at
  `a11c130`, and only the caption and the provenance move.
- **Local only.** The verdicts and the filled workbook stay in `scores/flag_review/round3/`, gitignored
  with the store. `FLAG_REVIEW_3.md` holds counts only.

**Open.** The router, which needs the owner's approval of a fresh dry-run estimate; and the Kaggle copy of
`scores/` (manual mode, without `scores/flag_review/`).

## D-164 — The router run approved at a $140 cap and launched on the final digit rule

**Date:** 2026-09-30 · **Status:** DECIDED (the owner: approve, cap $140) · DONE 2026-10-02 · **Evidence:** `router.py --validate`, `--dry-run` and `--status`; `scores/router.log`

- **Before the launch, all free:**
  - `router --validate` reproduces the pilot. All 299 prompts are byte-identical to those sent, and every
    figure matches its published value.
  - `router --dry-run` on the store re-scored under the final rule (D-160) gives 24,506 calls and 171,941
    steps under review. The estimate is $98.41 to $103.49 at Xiaomi's prices, $72.44 to $195.30 across
    the six endpoints, and about 33 hours at 8 workers. The count is unchanged from D-159's, since the digit flags
    change which steps a call sends, not the calls.
  - The reply store was empty.
- **The approval.** The owner approved a cumulative cap of $140. It covers the estimate at Xiaomi's prices
  with E5's longer replies (about $129 to $136, D-152). If routing lands on dearer endpoints, the run stops
  at the cap with every answer kept, and the owner decides then.
- **The launch.** It ran through `keepawake_run.ps1`, detached, with `router --yes --max-usd 140
  --workers 8`, at 2026-09-30T21:21:29Z. The paraphrase stream's own run went on beside it, untouched.
- **The first 53 replies** (`router --status`): no failure; 3,482 completion tokens on the mean against
  the basis's 3,702; $0.00392 a reply, about $96 over the run. The median call takes 59 s against the
  pilot's 39, so the run looks nearer 50 hours than 33.

- **The run** (`router --status`, `scores/router.log`). The owner asked for more speed, so the router
  was restarted at 16 workers on 2026-10-01 and at 24 later the same day. A laptop restart came
  between, and the run resumed from the store each time. The main pass ended at 2026-10-02T05:55Z with
  783 calls without a reply, 732 of them on dropped connections overnight. A retry at 16 workers made
  the 782 not yet given up and answered 780.
  - 24,503 of the 24,506 calls have a reply. The 3 without one are one each for Gemma, Qwen3-235B and
    Gemini; one key was given up after three failed runs.
  - Billed: $104.16 cumulative, under the cap. The mean is $0.00423 a reply, and 3,508 completion
    tokens against the basis's 3,702.
- **The results** (`analyze` at `142fa7f`). Q3's router columns are filled, and the attribution table
  gains its router column. Only those fields and the provenance move.
- **Backup (2026-10-02).** `backup_archive.py scores --exclude scores/flag_review scores/_replaced
  scores/paraphrase` wrote `full_run_scores_2026-10-02.zip`, 85 files, each checked against its source,
  sha256 `6e29439b…79c74919`. It holds the main and repeat stores, the stages' rows, E5's and the router's
  reply stores, and the audit files. The expert's flag-reading labels are not in it. Superseded stores
  can be regenerated, and the paraphrase arm has its own backup (D-162). A private Kaggle copy, uploaded
  in manual mode, was downloaded back, and all 85 members matched by path and hash.

## D-165 — Q5 read as a bound: the equivalence margin, the noise floor at the arm's size, the repeats beside it; the stream audited and the paper notes written

**Date:** 2026-10-01 · **Status:** DECIDED · **Evidence:** `analyze.py` (`paired`, `tau_noise_arm`, `vs_repeats`, `EQUIV_MARGIN`, `--selftest`), `results/RESULTS.md` Q5, `paraphrase_audit.py` and `PARAPHRASE_AUDIT.md`, `PARAPHRASE_PAPER_NOTES.md`

The owner asked whether Q5 is publishable. It is, as a robustness section answering a named reviewer
objection, once the null result is stated as a bound rather than as "nothing significant", once the
ordering's stability is read against a noise floor computed at the arm's own size, and once the
wording effect is set beside the run-to-run effect the decoding repeats measured. Three additions to
`analyze.py`, each reusing Q5's definitions (the template bootstrap, its seeds, the item-weighted
ordering), and nothing in the pre-registered test changed:

1. **Bounds.** `paired()` now returns, from the same bootstrap draws, the 90% interval beside the 95%
   one, read against `EQUIV_MARGIN = 0.05`: a model is within the margin when the whole 90% interval is
   (the two-one-sided-tests rule at 5%). The margin is the plan's largest detectable paired difference
   (5.1 points at 15% discordance, `ANALYSIS_PLAN.md` Q5); it was fixed after the point estimates were
   known (all within ±3.2 points) and before the intervals were computed, and the results and the paper
   notes say so. **Ten of eleven models are within ±5 points**; `qwen3-235b-a22b-2507` is not (−0.031,
   90% CI −0.067 to +0.005). Two intervals lie wholly below zero, `glm-5.3-flash` (−0.020, −0.038 to
   −0.004) and `claude-sonnet-5` (−0.018, −0.032 to −0.005): small decreases the data favour, within the
   bound, not significant after Holm.
2. **The noise floor at the arm's size.** The roster-wide floor (0.881) splits each template's 15 items
   in half and is not the comparison for an arm of 277 items over 115 templates. `tau_noise_arm` draws,
   for each template with k kept pairs, 2k of the main run's items and splits them into two halves of k,
   orders the models on each by the item-weighted mean as `kendall_boot` orders the arms, and takes tau,
   200 times: **median 0.783, quartiles 0.722 to 0.844, 5th to 95th percentile 0.636 to 0.917**. The
   arm's tau, 0.636 (95% CI 0.457 to 0.871), sits at that distribution's 5th percentile: consistent with
   sampling noise, at its lower edge. The reading for the paper: the tiers hold; the order within the top
   tier, five models within 0.012 of one another, is not resolved by this arm.
3. **Against run-to-run noise.** `vs_repeats` pairs each decoding repeat against the main run on the 300
   repeat items with the same `paired()`, and also both arms on the 92 kept pairs the two subsamples
   share. The paraphrase difference lies within the repeats' spread for `gemma-4-26b-a4b` and
   `gemini-3.1-flash-lite`; not for `qwen3-235b-a22b-2507` (−0.031 against repeats within ±0.008) and
   marginally not for `gpt-oss-20b` (+0.011 against +0.008). The sentence "rewording moves a model less
   than re-running it" is therefore not available for the roster; the four rows are reported as they
   are.

The self-test gained two assertions: a zero difference sits within the margin and a nine-point loss
outside it; the arm-size floor is finite. The results were regenerated from the clean tree.

**The audit.** `paraphrase_audit.py` recomputes every number the stream reports from the artifacts on
disk and compares it with what is committed: the manifest against the local pool (hashes, one prompt,
the original's hash the frozen pool's, every served text passing the six checks as the code stands,
the restored flag), the kits against the keyfile and the layer-2 roster (whole templates, own branch,
at most ten, the kit text the restored text), the returns against the assignment (one submission per
code from its expert, the kept rule recomputed, a note on every rejection), the traces (one final row
per item and model, both hashes the manifest's, one served model each as the main run's, the spend),
the scores (hashed to their traces, unusable rows the empty rows), the judge (a reply for every judged
trace), the results (eleven rows over the kept pairs, a clean provenance) and the backups (the archive's
checksum and members, the Kaggle listing). **28 checks, all pass**; `PARAPHRASE_AUDIT.md` is
the record. Its first run failed three checks on the script's own field names and none on the data.

**The paper notes.** `PARAPHRASE_PAPER_NOTES.md`: the figures to report with their sources, what to
state (the bound, what it answers, the ordering, the two small decreases, Qwen, the smaller arm and
why, the expert step's necessity, provenance), what not to state ("robust", "no contamination", "the
ranking is stable", "less than re-running", 450, the roster-wide floor, the E3 columns as findings,
per-branch conclusions), the margin's history, suggested wording for the section and the limitations,
and the figures worth drawing.

## D-166 — Related Work revised for the October submission: the May citations audited, the new literature added, the positioning claims the combination

**Date:** 2026-10-02 · **Status:** DECIDED (the text); OPEN (the authors' choices listed below) · **Evidence:** `docs/related_work_oct2026/` (`README.md`, `CITED_SET_AUDIT.md`, `NEW_PAPERS_AUDIT.md`, `CHANGES.md`, `notes/fact_check.md`)

Section 2 of the May submission was audited and rewritten. Each of the 46 works it cites, or names but
attributes to another reference, was read by three independent reviewer agents. A search of six topics
found 156 newer candidates. The revision cites 99 works: 67 in the main text and 32 in an appendix with
a comparison table. Every statement it makes about another paper was checked against that paper's text
(`verify_facts.py`, 177 of 177 claims found). What the audit changed:

1. **The section's central claim is withdrawn.** "Existing benchmarks remain limited to outcome matching"
   and the Introduction's "no existing benchmark verifies this process" are contradicted by seven of the
   May citations, among them PhysReason, FinChain, FEABench and TransportBench, and by newer engineering
   benchmarks: ThermoQA, TPS-CalcBench and EngVQA. The revision instead claims the combination, hedged
   with "to our knowledge": generated instances across five branches, gold traces from the template code,
   deterministic checks of intermediate values, and an evaluator validated against experts' step labels.
2. **Four citations pointed at the wrong paper.** These were "Glue", "BBH", "BIG-Bench" and the "Liu et
   al." CIRCUIT entry. Three characterisations were the opposite of the paper: ABench-Physics, LLM-SRBench
   and FEABench. LLM-SRBench and the two storytelling self-citations are dropped, and the reference list
   is regenerated from source records.
3. **FinChain and CIRCUIT become named precedents.** FinChain is described in the third person, for
   anonymity. A paragraph on process supervision and evaluator validity is added, as NEXT_CYCLE_REVIEW
   §4.2 asked.
4. **Templates are no longer presented as a remedy for saturation.** Akhtar et al. (2026) find such
   safeguards have limited effect. The private seed is presented as the contamination answer.

The reviewer panel was stopped on 2026-10-02 at the owner's request, to save cost. The remaining
candidates were assessed through the fact check (`notes/review_protocol.md`).

## D-167 — Related Work: the 18 unreviewed main-text papers read, the positioning holds; ten citing sentences or rows corrected, OPT-Engine added, two statements about EngTrace fixed

**Date:** 2026-10-02 · **Status:** DECIDED (the text) · **Evidence:** `docs/related_work_oct2026/notes/followup_review.md`, its 18 reviews in `reviews/candidate/`, `notes/fact_check.md` (200 of 200), `notes/venue_check.json`

D-166's section cited 18 new papers in its main text that no reviewer had read beyond the sentence
citing them. Each was reviewed once, 12 in the session and 6 by one reviewer agent whose findings were
re-checked against the papers' text. The claim that no benchmark combines generated instances across
engineering branches, gold traces, deterministic checks of intermediate values and an evaluator
validated against experts' step labels holds against Table A, the 100 cited works and the 86 uncited
candidates. What changed in substance:

1. **Saturation.** Akhtar et al. find that templated benchmarks saturate no differently from others
   (14 against 46, p = 0.10, their H6), and private test sets did not help either. The text now says
   that neither has measurably slowed saturation. D-166's conclusion stands: neither templates nor the
   private seed is presented as a remedy. The seed answers contamination.
2. **EngTrace's own description.** The digit rule recomputes each claim from the model's own numbers;
   it does not check arithmetic against the gold trace, as the text had said. And the judge decides
   "what these checks leave open", which covers the step router (D-110, D-113) as well as E5.
3. **Precedents.**
   - PRISM-Physics' step check is named as rule-based, with no LLM. It is the nearest physics
     precedent for EngTrace's deterministic checks; its validation is two experts' scores of 70
     solutions.
   - OPT-Engine (ICML 2026) is added. It generates operations-research instances, inventory and
     production among them, and is a closer reference for the industrial branch than the
     multiple-choice ORQA.
   - Table A's AtmosSci-Bench row is corrected: its generated items are multiple choice, and its LLM
     grader agreed with human graders on 92.79% and 93.02% of answers.
4. **Scope.** Seven citations now claim no more than their papers:
   - PRIME is college-level STEM, 18% of it engineering;
   - Długosz et al. re-analyse GSM-Symbolic only;
   - Lunardi et al. studied multiple-choice benchmarks;
   - Krumdick et al. require a correct reference;
   - on JudgeBench, many judges are near chance, not all;
   - MATH-Perturb probes memorised methods as well as instances;
   - ControlBench's experts examined the reasoning but did not grade it.
5. **Length.** The [cut first] sentence on physics re-grading is cut, and its paper is cited for
   rule-based grader errors. The main text is 704 words (it was 723). The section cites 100 works:
   68 in the main text and 32 in the appendix.

**Venues** (`check_venues.py`, `notes/venue_check.json`).

- ERI is confirmed: the DOI on its arXiv page is Computers & Industrial Engineering 221 in Crossref.
- PRIME, PRISM-Physics, HiPhO and MATH-Perturb are confirmed in OpenReview or the ACL Anthology.
- Two papers listed as preprints have appeared, and are now cited at their venues: SymPyBench (EACL 2026
  Industry Track) and AtmosSci-Bench (NeurIPS 2025 Datasets and Benchmarks).
- EMNLP 2026 for Huang et al. and Długosz et al. still rests on arXiv comments.

The Introduction sentence in `CHANGES.md` §2 no longer cites Akhtar et al. for BIG-Bench Hard, which
they count as unsaturated.

## D-168 — The full run reviewed end to end before the paper: the record reproduces; the answer check misreads pi-fractions and one-quantity-two-units lines, measured and not yet adopted

**Date:** 2026-10-02 · **Status:** MEASURED · OPEN (adoption is the owner's call, as D-138 and D-147 were) · **Evidence:** `full_run_28092026/answer_form_audit.py`, `ANSWER_FORM_AUDIT.md`; `scores/answer_form_audit_changes.jsonl` (local); the free checks named below

**What was checked, all free.** `freeze --check-files` (FILES OK); `run_traces --status` (11 models, missing 0,
$247.917 in the rows); `score --status` (main at `eea846e`, clean, 12 models; the evaluator blobs its CONFIG
names are HEAD's); `judge --status` (8,032 main-arm calls answered, none without a reply, 1,227 on the
paraphrase arm; $40.86 over both); `router --status` (24,503 of 24,506 answered, one key given up; $104.16);
`analyze --selftest` (all pass); `judge --validate` and `router --validate` (both reproduce the pilot);
`validate_scorer` (gold clean; all ten published figures reproduce); `trace_review` (identical to the committed
review apart from its archive line, because the script now picks the newer repeats archive, a cosmetic
defect); and `analyze` regenerated `results/` byte-identical apart from the provenance commit. Q1's table was
also recomputed from the store by a separate script: every count and score matches, including the 7 odd
finish reasons of D-148.

**What a read of the verdicts found.** Twelve "incorrect" verdicts drawn at random from the top five models'
answered traces were read by hand: six were right answers in another form. The forms were measured by
`answer_form_audit.py` the D-137 way, on the gold, on the pilot's 300 labelled traces and on every full-run
trace, with `answer.py` unchanged:

- **Pi-fractions (reading P).** `answer.PI_EXPR` reads `(841*pi)/2447` as 841*pi, which the template also
  computes, so the target itself is wrong, and reads `\dfrac{841\pi}{2447}` as pi. Normalising the text first
  so the existing rule reads a*pi/b: gold 2,250 of 2,250; no pilot verdict moves, because the form does not
  occur in the 15 pilot templates, so the experts cannot arbitrate (as under D-138 and D-147); in the full run
  142 verdicts change, 138 to correct, 3 from correct to partial and 1 from incorrect to partial (three-part
  items answered in one part), all on `continuous_to_discrete_conversion` (75) and
  `decimation_aliasing_analysis` (67), both signals_and_systems, the domain every model scores lowest on.
  Answer scores rise by 0.004 to 0.006 per model; 13 of DeepSeek V4.1 Flash's 40 incorrect verdicts and 11
  of Claude Sonnet 5's 65 are of this kind.
- **One quantity in two units (reading A, on top of P).** A scalar gold line "0.088960 radians, or 5.0970
  degrees" yields one target, the degrees, and a trace in radians is incorrect. Accepting another computed
  value on the line that is the target in another unit (ratio 180/pi, 2*pi or a power of ten): 56 more
  verdicts to correct, 55 on `composite_shafts_series` and 1 on `angle_of_twist`; Kimi K3 gains 10. A first
  version that accepted any computed value on the line was rejected before adoption: it credited traces right
  on one of two different quantities (`vdw_solve_for_volume`, `levenspiel_plot_interpretation`).
- **Two-quantity lines typed scalar or array.** `answer.targets` keeps the last computed value for `scalar`
  and `array` items, so the templates the audit lists (`vdw_solve_for_volume`, `levenspiel_plot_interpretation`,
  `standing_wave_formation`, `finite_convolution`, `gas_viscosity_kinetic_theory` on one item) are verified on
  one of the two or more values their gold line states: a trace right on the last alone is correct, one right
  on the other alone is incorrect. A property of the type, not of a reading; it goes beside D-145's list, or
  the type changes.
- **Symbolic answers written as their value (reading S, counted only).** On `ber_estimation_mary`, 73 of the
  roster's 117 incorrect or partial verdicts state the numeric value of the gold's `a * Q(b)`. D-138 scores a
  symbolic answer by the numbers it states, and the nine symbolic templates are a sensitivity row in RESULTS.md;
  this count says what that row absorbs.

**If P and A are adopted.** The change is confined to `answer.py`; E5 and the router read E3 and the digit rule,
so no reply is bought again. `score --replace` on the main, paraphrase and repeat stores, `judge --score`,
`router --score` and `analyze` are free, and the entry records it as D-137 to D-139 did. Under P+A the top
five lie between 0.965 and 0.976, still within 0.011 of one another.

**Also found, no action taken.** Four roster models wrote no reasoning tokens at provider-default decoding:
`gemma-4-26b-a4b`, `qwen3-235b-a22b-2507`, `gpt-5.4-mini` and `gemini-3.1-flash-lite`; the other seven did, at
a median of 718 (Claude Sonnet 5) to 2,473 (GLM-5.3) per trace. D-122 chose the providers' defaults; Appendix
P must state what each model ran with, and the closed tier's scores are read against it.

## D-169 — The answer check corrected after the run was read: pi-fractions, a scalar gold line's second unit and num/den fractions; every changed verdict read, the change measured between pinned commits, the stores re-scored

**Date:** 2026-10-02 · **Status:** DECIDED (the owner: adopt readings P, A and N; C declined) · DONE (the change, the gate, the measurement, the re-score, the stage rows, the analysis, the archive) · **Evidence:** `evaluator_pilot_17092026/evaluators/answer.py` (docstring and self-test, 64 cases); `full_run_28092026/answer_form_audit.py` (`ANSWER_FORM_AUDIT.md`, the readings; `--fix`, `ANSWER_FORM_FIX.md`, the shipped change between `e116b4f` and `5415f61`); `SCORER_VALIDATION.md`; `GOLD_VALIDATION.md`; `results/RESULTS.md` (`analyze.py` at `2950875`); `scores/answer_form_fix_changes.jsonl` and `scores/rescore_d169.log` (local)

**The gate.** Before anything changed, every verdict the readings move was read by hand with the gold's true
value beside the trace's answer segment: 198 under P and A, 238 once reading N was added. 233 of the 234
credits are right answers. One is not: `gpt-5.4-mini` on `continuous_to_discrete_conversion#1` writes
"ω ≈ 0.7995π rad/sample (≈ 2.511 rad/sample)", a wrong answer, and is credited because the bare coefficient
0.7995 lies within 0.2% of the true 0.7987 rad/sample. The four other moves are three-part aliasing items
answered in one part, scored partial as D-145 scores every multi-target item; one of them, Gemma's
`cos((π/2)n)` for `cos((3π/2)n)`, is equivalent for integer n and is a limitation. Two refinements were
measured against the wrong credit:

- Reading A as first written, matched with the check's unit factors, credited a cut-off Gemma trace of 832
  characters (`composite_shafts_series#16`, finish reason `error`, one of D-148's seven) holding `πd⁴/32` and
  no answer, through 32 ± 1 at the 1/60 factor. Matched at unit scale it does not, and no right credit is lost.
- Reading C, a π-coefficient is not a standalone value, removes the wrong credit and two accidental credits of
  wrong answers on `bpsk_energy_basis`, but strips the accidental credit of four right answers on
  `cd_dc_system_analysis` (a phase written −0.75π for a 0.0075 s shift) and of two on
  `continuous_to_discrete_conversion#10`: six outcomes worse for three better. Not adopted. The one wrong
  credit is recorded as a known limit of a check that credits any number within 0.2% of the target on the
  Answer line, and is pinned in the self-test so that a change to it is noticed.

The owner chose P + A + N on 2026-10-02, with the stop condition of the gate met and explained.

**What changed in `answer.py`, and nothing else.**
- **P.** A pi-fraction is one value in every shape, `(841*pi)/2447`, `\dfrac{841\pi}{2447}`, `841\pi/2447`,
  `(3/17)\pi`, `0.4\pi`: `PI_EXPR` allows the parentheses, `_pi_text` resolves LaTeX's `\pi`, `\frac`, `\left`,
  `\right` and `$`, and `(a/b)π` is rewritten `a*pi/b`; `pi*n` and `pi(` are still not values. Until now the
  parenthesis stopped the rule at 841π, a value the template also computes, so the gold's own target was
  wrong on the seven `continuous_to_discrete_conversion` items with no milestones and on
  `decimation_aliasing_analysis`.
- **N.** A pi-fraction on the gold's answer line whose numerator and denominator are both computed, the
  reduced fraction's `num` and `den` that the other eight `continuous_to_discrete_conversion` items expose as
  their only milestones, is one quantity: the pair of integer targets becomes the fraction's value, so a trace
  stating ω as a decimal or unreduced is right.
- **A.** A scalar gold line stating one quantity in two units ("0.088960 radians, or 5.0970 degrees") accepts
  either: the second is a rendering, related to the target by 180/π or 2π either way, matched at unit scale
  (`match(scales=(1.0,))`). Powers of ten are not renderings: two different quantities can differ by one, and
  the unit factors already serve the target.

`milestones.py` is untouched, so E3, E5's prompts and the router's are what they were and no reply was bought
again. 24 new self-test cases pin the forms, the counter-cases and the known limit; 64 of 64 pass.

**Measured between pinned commits** (`answer_form_audit.py --fix`, `ANSWER_FORM_FIX.md`: `answer.py` at
`e116b4f` against `5415f61`, both from git, on the same milestones). Gold 2,250 of 2,250 before and after.
The pilot's 300 traces: 0 verdicts move and agreement stays 0.982 / 0.927, so the published figures stand;
the forms occur in none of the 15 pilot templates, so the experts cannot vouch for them, and the paper says so
as it does for D-138 and D-147. The full run: 238 verdicts change, the 238 the gate read and no other, with
the same new verdicts: 234 to correct, 3 correct → partial, 1 incorrect → partial, on four templates,
`continuous_to_discrete_conversion` 116, `decimation_aliasing_analysis` 67, `composite_shafts_series` 54,
`angle_of_twist` 1. Per model 12 (GLM-5.3) to 32 (Gemma 4 26B, Qwen3-235B-2507) verdicts; answer scores
rise by 0.005 (GLM-5.3-Flash) to 0.014 (Gemini 3.1 Flash-Lite).

**Re-scored, all free** (`scores/rescore_d169.log`, 2026-10-02 12:32 to 13:16 UTC, at `2950875`, tree
clean): `score --replace` on `main` (12 models), `paraphrase` and `repeat1` to `repeat3`, the old stores
archived under `scores/_replaced/`; `judge --score` and `router --score` for `main` and `judge --variant
paraphrase --score`, every sent call still with its reply (the router's 3 without one as in D-164);
`analyze`; `gold_validation`, unchanged. `validate_scorer` regenerated with the change named; no figure moves.

**What moved in RESULTS.md** (`analyze.py` at `2950875`):
- **Q1.** DeepSeek V4.1 Flash 0.976, Kimi K3 0.974, Claude Sonnet 5 0.972, GLM-5.3-Flash 0.969, Muse Glimmer
  0.967, GLM-5.3 0.948, Qwen3-235B-2507 0.886, Gemini 3.1 Flash-Lite 0.872, Gemma 4 26B 0.864, GPT-5.4 mini
  0.847, gpt-oss-20b 0.814. Kimi K3 and Claude Sonnet 5 change places; the top five lie within 0.009. 32 of
  55 pairs hold (31 before); the fully-solved test agrees on 54 of 55 (55), McNemar on 47 of 55 (46).
- **Q2.** 4 of 11 under Welch with Holm (3 before; Gemma 4 26B joins gpt-oss-20b, GPT-5.4 mini and Gemini 3.1
  Flash-Lite at p 0.043); 9 of 11 under the planned permutation (6).
- **Q3.** Only the group memberships move, since 238 traces change group; no E3, E5 or digit-rule figure does.
- **Q5.** Paired differences −0.031 (Qwen3-235B-2507) to +0.020 (Muse Glimmer); none survives Holm; 10 of 11
  within ±5 points as before, Qwen not; Muse's 90% interval now lies wholly above zero (+0.005 to +0.036),
  beside GLM-5.3-Flash's and Claude Sonnet 5's wholly below; Kendall's τ 0.673 (CI 0.455 to 0.881) against a
  roster-wide floor of 0.891 (0.636 and 0.881 before).
- **Domains.** signals_and_systems moves from 0.777–0.920 to 0.883–1.000 and is no longer every model's
  weakest domain; thermodynamics, unchanged at 0.517–0.922, is the lowest for nine of the eleven. Electrical
  for the top five: 0.947–0.981 (0.912–0.954 before).
- **Sensitivity.** The ordering's τ with the headline: 0.927 at half and at double tolerance (0.891 and 0.964
  before), 0.964 on the fully-solved rate (0.891), 0.991 under the half-unit window (0.917), 1.000 without the
  shortcut templates and under the whole-trace reading, 0.964 without the symbolic templates, 0.636 with
  unusable rows excluded (0.600).
- **Repeats.** SDs 0.004 to 0.013 (0.003 to 0.016). The set-aside Qwen3.8-27B: 0.937 (0.924).

**Backups.** `~/EngTrace_private_backup/full_run_scores_2026-10-02.zip` was rewritten with the re-scored
stores and the rebuilt stage rows: 88 files, each checked against its source, sha256 `1a3df7ec…dd3e07`. It
replaces the morning's archive of the same name (D-164), whose stores are the superseded ones under
`scores/_replaced/` and regenerate from `e116b4f`. Uploaded the same day, with the owner's go-ahead, as a new
version of the private Kaggle dataset `ayeshaiq/engtrace-full-run-scores` (`kaggle datasets version`, never
`--public`), downloaded back and compared: 87 of the 88 members match their local source by path and SHA-256,
and the 88th, `scores/rescore_d169.log`, is the live log's first 5,371 bytes, because the chain appended its
archive-step report to the log after the archive had copied it; the archive's own `.sha256` matches. The
check folder was removed afterwards.

**For the paper.** The check was corrected after the run on a reading of its verdicts, as D-137 to D-139 and
D-147 were, and the paper reports it that way. The `array` type and the two-quantity scalar lines of D-168
stay as they are and are listed beside D-145's multipart templates. `trace_review.py` now reads the dated
main-run archive only (D-168's cosmetic defect); the committed review regenerates identically.

## D-170 — Reporting addition and documentation pass after D-169: the reasoning score gets a headline table; the paper notes and the next-steps file are written; NEXT_CYCLE_REVIEW is marked stale

**Date:** 2026-10-02 · **Status:** DECIDED · **Evidence:** `full_run_28092026/analyze.py` (`q3_overall`, docstring), `results/RESULTS.md` (Q3, "Milestone coverage over all traces"), `RESULTS_PAPER_NOTES.md`, `docs/EVALUATION_NEXT_STEPS.md`, `ANALYSIS_PLAN.md` (dated note), `EVALUATION_GUIDE.md`, `README.md`, `answer_form_audit.py` (`--readings` guard), `docs/NEXT_CYCLE_REVIEW.md` (status note)

- **The reasoning score had no headline.** Q3 reported E5 coverage on the wrong-answer traces only, as the
  plan asked ("what process scores add beyond the answer"). The paper's emphasis is the reasoning trace, so
  `analyze.py` now also prints, per model with template-level intervals, E3 and E5-strict coverage over every
  trace whose item has milestones (2,180 per model), an unusable trace scoring what it reached and an empty one
  nothing, and over the readable traces alone, with the digit rule's and the router judge's flag rates on fully
  solved traces repeated beside them. Descriptive, outside the plan's tests, labelled as added after the data.
  E5-strict runs from 0.816 (gpt-oss-20b) to 0.923 (Claude Sonnet 5); the intervals of the top eight overlap;
  the order differs from the answer score's (Claude first on coverage, DeepSeek V4.1 Flash first on answers),
  and the caption says what coverage measures: progress through the gold derivation, not the absence of error
  (RESULTS_X1 Finding 5). Results regenerated at `2950875`; only the new table and the header move.
- **`RESULTS_PAPER_NOTES.md`** says, for the whole run as `PARAPHRASE_PAPER_NOTES.md` does for Q5, which figures
  to report with their sources, what to state, what not to state, suggested wording for the results, the
  evaluator section and the limitations, and the figures worth drawing. Its rule is the project's: tiers, not
  an order within the top tier; coverage is coverage; flag rates are flags; four models ran without reasoning;
  no claim the data do not carry.
- **`docs/EVALUATION_NEXT_STEPS.md`** carries each remaining analysis and experiment with what it answers, the
  method, the unit, the cost and where its output lands: the Wilcoxon and bootstrap comparison of coverage across
  models the July rebuttal promised, branch and level intervals, leave-one-judge-out with controls, the X2 bias
  estimate, the threshold appendix, attribution by error type, the shortcut audit, the decoding table, power
  statements, the expert readings (the answer check on this roster, the error sample, the judge), and the paid
  conditions (reasoning-on, judge swap, flagship anchor, tool use). It supersedes section 10 of
  `NEXT_CYCLE_REVIEW.md`, which now carries a status note saying it is the September plan and no longer
  maintained.
- **Documentation made consistent with D-169:** the plan's dated note; the guide's store commit (`2950875`),
  its "Done means" list closed for the router, the results, the records and the archive, and a closing
  paragraph pointing to the review, the notes and the next steps; the README's account of the audit as a
  frozen design record and of the fix record; `answer_form_audit.py` refuses its default mode without
  `--readings`, since re-running the readings against the corrected code would overwrite the record with
  nothing; `validate_scorer.py`'s sentence names D-169 and `SCORER_VALIDATION.md` is regenerated.

## D-171 — The remaining incorrect verdicts measured: two Advanced chemical templates under-specify the answer, 17 templates require the digits; what it means for Q2, Q3 and the experts' reading

**Date:** 2026-10-02 · **Status:** MEASURED · OPEN (what the paper states, and whether a chemical expert reads the two templates, is the owner's call) · **Evidence:** `full_run_28092026/residual_incorrect.py`, `RESIDUAL_INCORRECT.md`; `docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` section 5

**Why.** D-168 read twelve of the top five models' "incorrect" verdicts and found six to be right answers in another
form; D-169 fixed those forms. Nothing had measured what the verdicts that remain look like, and the top models'
wrong-answer subsets (26 to 50 traces) carry Q3's process-score analyses for those models.

**What was measured, all free.** For every answered, usable trace the store labels `incorrect` (1,558 across the
roster at `2950875`), the closest approach of any number on the Answer segment to any numeric target under the check's
own unit scales, bucketed by relative distance; symbolic and word-only items counted apart. Roster-wide: 120 within
0.2%, 406 at 0.2% to 1%, 336 at 1% to 5%, 228 at 5% to 20%, 365 beyond 20%, 90 symbolic, 11 word-only, 2 with no
number. Every one of the 120 within the tolerance sits on one of the 17 templates whose question prescribes the
rounding (`exact`), where a last-digit miss scores 0 by the experts' rule; the (Q, R) question, for one, prescribes
every intermediate rounding, and the top models compute at full precision and land 1 to 4 units off. Those are
legitimate zeros and a finding about the models. The top five's 171: 38 symbolic (34 on `ber_estimation_mary`,
D-168's reading S), 27 within 0.2% (all exact-digit items), 31 at 0.2% to 1%, 40 at 1% to 5%, 33 beyond 5%, 2
word-only; by template `work_isothermal_virial` 38, `ber_estimation_mary` 34, `adiabatic_flame_temperature` 19,
`qr_policy_one_iteration` 13.

**Two Advanced chemical templates do not pin the answer to the check's tolerance.** `work_isothermal_virial` asks for
the work to compress one mole isothermally and reversibly "based on the virial equation of state truncated to two
terms" and names neither the closed-system form (W = -∫P dV with Z = 1 + B/V, the gold's) nor the flow form
(W = ∫V dP with Z = 1 + BP/RT, which is RT ln(P2/P1) + B(P2 - P1)); computing the flow-work value from the gold's own
stated B and ideal-gas work and matching it with `answer.match`, **53 of its 97 incorrect traces state the flow work**
(0 the ideal-gas work, 44 neither); three top models give 11,867 J/mol to five figures where the gold says 11,590.
`adiabatic_flame_temperature` gives no heat-capacity data and asks for an estimate; the gold iterates a polynomial Cp
with uncited coefficients; 48 of its 89 wrong answers lie at 0.2% to 5%. Both were certified by three chemical experts
(round 2 and round 4); the wording was not the objection either time. Both are Advanced: the roster loses 237 of 330
points on them, and without them the top five's Easy-to-Advanced gap is +0.029 to +0.047 instead of +0.055 to
+0.073, GLM-5.3's +0.087 instead of +0.134, the four significant models' +0.105 to +0.157 instead of +0.151 to
+0.199 (descriptive, no interval).

**What changes and what does not.** No score, store or result changes: the pool is frozen (D-114), inference is done,
and the sensitivity table already carries the symbolic templates. Proposed for the paper: the two templates as a
stated limitation with a Q2 sensitivity row (A2 computes its interval); the 17 exact-digit templates named with the
rule; the top models' wrong-answer subset sizes stated beside every Q3 row that uses them, or those rows restricted to
models with more than 100 wrong answers. For the experts: one chemical expert reads the two questions and golds, three
traces each, in the same request as B1 to B3 (next steps B4). The per-template near-miss table names the next
templates to look at (`vdw_solve_for_volume`, `pfr_volume_changing_rate`, `best_hydraulic_rectangular_section`,
`manning_rectangular_discharge`, `annulus_flowrate`, `pitzer_correlation_z`). The measurement judges no verdict; the
reading of the three templates above is this session's and is recorded as such.

## D-172 — The experts' reading request built: B1 to B4 as one kit per expert, an app and a guide; sending is the owner's step

**Date:** 2026-10-02 · **Status:** BUILT · OPEN (the owner sends the kits and sets the return date; the sizes are flags and can be cut before sending) · **Evidence:** `full_run_28092026/expert_kits.py` (`--build`, `--score`, `--selftest`), `reading_app.py`, `EXPERT_READING_GUIDE.md`, `EXPERT_REQUEST.md` (the composition, counts only); `expert_request/` (local: tasks, keyfile, kits)

**Why now.** B1 to B4 read the frozen store and traces as they stand at `2950875` and depend on none of the free
analyses (section A), while the experts' time is the long pole and D-171 made B1 and B4 urgent. The owner asked for the
request to go out first, with the instructions in a markdown file and an app, as for the paraphrase check.

**What was drawn** (seed 20261002; `EXPERT_REQUEST.md` has the tables): 418 items, 1,004 readings.
- *B1, 150 final answers* from the top five models: 75 the check called incorrect (49) or partial (26) and 75 it called
  correct, each half round-robin over templates so no template dominates (103 templates; symbolic 17, multipart 30,
  scalar 79). The expert sees the problem, the reference answer, the model's final answer and both in full; not the
  model or the verdict. Two readers each, rotating pairs.
- *B3, 100 milestones*: 50 the judge ruled REACHED and 50 MISSING, round-robin over the eleven models and their
  templates, one per trace (53 templates). The expert sees the quantity's name and value, the reference and the trace;
  not the verdict. Two readers each.
- *B2, 160 wrong answers*: 40 each for `claude-sonnet-5`, `gpt-5.4-mini`, `gpt-oss-20b` and `gemma-4-26b-a4b`, split
  over levels in proportion to each model's wrong answers (Claude's are mostly Advanced: 29 / 10 / 1), round-robin
  over templates (79). The May version's six categories by its stop-at-first-yes hierarchy, plus "No error: the answer
  is correct, or the question admits it" and "Incomplete", with a required excerpt for an error. Three readers each,
  for the agreement figure May reported.
- *B4, 8 templates*, chemical and civil: the two D-171 templates with one problem and three top-model answers each,
  and six near-miss templates (`vdw_solve_for_volume`, `pfr_volume_changing_rate`, `annulus_flowrate`,
  `pitzer_correlation_z`; `best_hydraulic_rectangular_section`, `manning_rectangular_discharge`) with one near-miss
  answer each; two to four questions on whether the wording pins a single answer at 0.2%. All three experts of the
  branch.

**Who reads what.** Each item goes to experts of its own branch (the layer-2 roster of 15); queues run template, answer,
milestone, error, shuffled within each block per expert. Load per expert: 53 to 81 items (chemical 79 to 81, civil 62
to 64, electrical 70 to 71, industrial 67 to 69, mechanical 53 to 54); B2 is the bulk (`--b2-per-model 30` would cut
about 15 items from each kit). Codes are opaque; the keyfile stays local.

**What the app and the guide are.** `reading_app.py` is `paraphrase_app.py` generalised to the four kinds: plain-text
display, one radio per question, a required note where a verdict needs one, a required excerpt for an error category,
and `<id>.jsonl` written next to it. `EXPERT_READING_GUIDE.md` is the one-page instruction per kind; the kit carries it
as `guide.md` with a `README.txt` on running it. `--selftest` builds a small kit in a temp folder, checks that no kit
holds a model name, a verdict or an item id, drives the app through every kind with Streamlit's test harness, submits
one item and scores the return.

**What `--score` will report**, counts only, in `EXPERT_REQUEST.md`: B1 agreement with the check three-way and on each
side, per branch and between readers (Cohen's kappa); B3 the share of REACHED the experts confirm and of MISSING they
confirm, with the "route does not need it" share; B2 the category distribution per model and level, item majorities,
the "No error" share and Fleiss' kappa; B4 the three experts' answers per question. The experts' notes stay local.

**Not decided here.** The send itself, the return date, and whether B2 is cut. No score, store or result changes.

*Amended 2026-10-02, at the owner's request:* the kits are delivered as one bundle per expert with a sub-folder per
task (`B1_final_answers`, `B3_milestones`, `B2_wrong_answers`, `B4_templates`), each with its own copy of the app, which
takes its title from the one kind it finds, a guide made of the general part plus that task's section, and its own return
file; the tasks can be sent or dropped separately, and `--score` reads a folder tree. The samples and readers are unchanged.

## Open decisions

| # | Decision | Needed before |
|---|---|---|
| D-172 | Send the fifteen kits in `full_run_28092026/expert_request/dist/` (app.py, README.txt, guide.md and the kit_<id> folder each) with a return date; cut B2 first if the load is too much (`--b2-per-model`); when the files come back, `expert_kits.py --score <folder>` | the paper's evaluator section, Q3 and the D-171 limitation |
| D-171 | What the paper states about the two under-specified chemical templates and the 17 exact-digit templates (a limitation and a Q2 sensitivity row, as proposed), and whether one chemical expert reads `work_isothermal_virial` and `adiabatic_flame_temperature` with B1 (next steps B4) | the paper's sections 5 and 6; the experts' request |
| D-170 | The analyses of `docs/EVALUATION_NEXT_STEPS.md` section A (free) in the order given; the expert request of section B, batched; a decision on each paid condition of section C after its dry run | the paper's results section; the experts' availability; the owner's approval per condition |
| D-169 | ~~Upload the rewritten `full_run_scores_2026-10-02.zip` as a new version of the private Kaggle scores dataset (manual mode)~~ **Done 2026-10-02: new version of `ayeshaiq/engtrace-full-run-scores`, downloaded back, 87 of 88 members identical and the log a prefix of the live one (D-169, Backups).** In the paper: the `array` type and the two-quantity scalar lines are listed with D-145's multipart templates; the one wrong credit is a stated limit | the paper's scoring statements |
| D-168 | ~~Adopt readings P and A of the answer check (`ANSWER_FORM_AUDIT.md`), then `score --replace` on main, paraphrase and repeats, `judge --score`, `router --score`, `analyze`, and a DECISIONS entry; or report both as sensitivities~~ **Adopted 2026-10-02 with reading N, after every changed verdict was read (D-169); re-scored and regenerated.** ~~Decide how `array` items and the two-quantity scalar lines are scored, or list those templates with D-145's~~ **Listed with D-145's (D-169).** ~~Fix `trace_review.py`'s archive glob~~ **Fixed (D-169)** | — |
| D-167 | Confirm EMNLP 2026 for Huang et al. and Długosz et al. once the proceedings appear (ERI is confirmed). The authors' call, outside the Related Work folder: `JUDGE_SELECTION.md` should say self-preference "can be" more than 50%, and the cliff should be reported as a difference between the templates in each tier (`notes/followup_review.md`) | the reference list; the paper's sections 3.3 and 5 |
| D-166 | The authors' calls on the Related Work revision (`docs/related_work_oct2026/CHANGES.md`): adopt the replacement Introduction sentences (§2); whether the Limitations acknowledge Mondorf et al.'s finding on randomly sampled instances (§4); confirm three single-source venues (§6) | the paper's sections 1 and 2 |
| D-161 | ~~Run: E5 on the paraphrase arm, about $3 at 316 items (`PARAPHRASE_RUNBOOK.md` section 4), for Q5's E5 column~~ **Done 2026-09-30 (D-162): $4.56, every call answered.** ~~The private Kaggle copy of the arm's archive~~ **Done 2026-10-01: `ayeshaiq/engtrace-full-run-paraphrase`, private, from `full_run_paraphrase_2026-10-01.zip` (80 files, the experts' returns included), downloaded back and matched member by member** | Q5's E5 column in the paper |
| D-160 | ~~Round 3 of the flag reading: the 190 claims in `scores/flag_review/round3/` (workbook and instructions), then `flag_sample --round 3 --merge <file>` and `--score --read-by "<who>"`, giving the corrected rule's precision on this roster; no round 4~~ **Read 2026-10-01 by the same domain expert: precision 0.905 (D-163, `FLAG_REVIEW_3.md`)** | — |
| D-159 | ~~Round 3 of the flag reading: the 181 claims in `scores/flag_review/round3/`~~ **Redrawn 2026-09-30, unread, after two review agents' corrections (D-160)** | — |
| D-158 | ~~A second fix of the digit rule's parser from round 2's notes, measured by a round 3 (no round 4), or 0.752 reported as measured~~ **Decided 2026-09-30: the second fix, adopted with no pilot step moving away (D-159); round 3 drawn** | — |
| D-156 | ~~Round 2 of the flag reading: the 206 claims in `scores/flag_review/round2/` (workbook and instructions), then `flag_sample --round 2 --merge <file>` and `--score --read-by "<who>"`, giving the fixed rule's precision on this roster~~ **Read 2026-09-30 by the same domain expert: precision 0.752 (D-158, `FLAG_REVIEW_2.md`)** | — |
| D-154 | ~~Fix the digit rule's parser the D-137 way, then read a fresh sample; or keep the rule and report its measured 0.505 (per model 0.10 to 1.00) in Q3~~ **Decided 2026-09-30: fixed, and adopted with a documented exception (D-156); round 2 drawn** | — |
| D-150 | ~~The author's reading of the 220 sampled digit-rule flags (`scores/flag_review/sample.csv`), then `flag_sample --score` and a checker fix if one is needed~~ **Read 2026-09-30 by a domain expert, in the D-153 workbook: precision 0.505 against the pilot's 0.750 (D-154, `FLAG_REVIEW.md`)** | — |
| D-147 | The inclusive last-digit boundary is adopted; the owner may reverse it (one commit). Whether the paper reports the half-unit or whole-trace reading as anything more than a sensitivity would need an expert spot-check of the format-only partials | the paper's tables |
| D-146 | The paper reports the cliff per model with its interval and, if it states a count, the Welch count with the planned one beside it (as RESULTS.md now does) | the paper's section 5 |
| — | The router's funding: with everything else run the round lands at about $450 on the account's basis (about $480 at E5's dearest endpoint); the router adds $72 to $195 and does not fit in $500. **2026-09-30 (D-151, D-152): with E5 at its projected $40 and the four-model repeats, about $465; the router's estimate rests on the pilot's output lengths, about $129 to $136 at Xiaomi's prices at E5's ratio (an extrapolation). The owner: the router waits for the author's flag reading** | the router's run |
| — | A stratified read of the digit rule's flags on this roster (Qwen3-235B-2507 has 1,547 on 18,785 claims; 16 of deepseek's 24 sit in one template), the pilot's own rule before Q3 is reported; author time, no code | Q3 in the paper |
| D-141 | ~~Run: writing the 450 paraphrases with Mistral Large 3, $0.15 to $0.45 (`EVALUATION_GUIDE.md`, the paraphrase arm)~~ **Done 2026-09-30 (D-157): 316 of 450 pass after four prompts and the notation restore, 29 templates lose all three items, $0.416 in all; the 15 kits are built** | — |
| D-141 | ~~Run: the paraphrase run over the eleven models, about $49.2 on the main run's bills for the same items~~ **Done 2026-09-30 (D-161): 316 items, $34.17, every served model as the main run's; scored and analysed, Q5 provisional until the experts return** | — |
| D-141 | ~~Run: the three decoding repeats of `gemma-4-26b-a4b`, about $0.28~~ **Done 2026-09-30 over four models, by the owner's decision (D-151): $1.863 in the rows; `RESULTS.md` prints the table** | — |
| D-138 | The stricter match rule removes 252 credits, 129 of them on symbolic answers the check cannot verify either way: keep it (as decided) or reverse it | the paper's tables |
| D-142 | ~~Run: E5 over the eleven models, 8,032 calls, about $26.56 at Xiaomi's prices ($18.59 to $49.71), about 7 hours (`EVALUATION_GUIDE.md`, step 3)~~ **Done 2026-09-30 (D-155): every call answered, $36.30 in the rows; `RESULTS.md` carries the E5 columns** | — |
| D-142 | Run: the step router over the eleven models, 24,506 calls, about $98 to $103 at Xiaomi's prices ($72 to $195), about 33 hours at 8 workers (`EVALUATION_GUIDE.md`, step 4). **The owner, 2026-09-30: after the author's flag reading; the estimate rests on the pilot's output lengths (D-152). Re-priced on the fixed flags: $98.41 to $103.53 at Xiaomi's prices, about $129 to $136 at E5's ratio (D-156); on the second fix's, $98.41 to $103.63 (D-159). Approved 2026-09-30 at a $140 cap on the final rule's flags ($98.41 to $103.49) and run: 24,503 of 24,506 calls answered, $104.16 (D-164)** | — |
| D-108 | ~~Whether the review app shows questions and solutions as plain text, as the models read them, and whether the 65 templates judged through its Markdown rendering get a plain-text look; fixing `signal_operations`'s origin marker~~ **Closed by the owner 2026-09-28 (D-114): the certification is closed and no template changes** | — |
| D-107 | ~~Fix the five templates round 2 objected to~~ **Fixed 2026-09-26 (D-108) and re-certified in round 3 (D-109): all five approved by all three** | — |
| D-106 | ~~Whether the screen re-judges the changed templates (a few cents, targeted; not run before round 2); the plasma row's per-row tag~~ **Closed by the owner 2026-09-28 (D-114): no re-judge** | — |
| D-095 | ~~Layer 2 roster and dates~~ **Round 1 run 2026-09-24/25 with the pilot's 15 experts (D-106).** **Round 2 run 2026-09-25 (D-107).** **Round 3 run 2026-09-26 (D-109).** ~~Adjudication rule for split verdicts~~ **Not needed for this pool: after round 3 no split verdict is in force; all 150 templates are certified unanimously (D-109)** | — |
| D-094 | ~~Whether to narrow the P2/Pc range in `work_isothermal_virial` (now 50% redraw; an expert's thermal-stability objection, D-106, turns on the same range) and accept the 41% stability redraw in `floating_object_submersion_depth`; the residuals in `pass1_fixes.md`~~ **Closed by the owner 2026-09-28 (D-114): the redraw rates stand** | — |
| D-092 | ~~Answer-display lengthening and gold-movement sign-offs~~ **Decided by the owner 2026-09-23:** the two lengthened answers stay; `beam_internal_moment`, `terzaghi_strip_footing_bearing` and `effective_stress_profile` answers are lengthened too; the gold movements are accepted; and a scoped third round removes the census ties at 3% and above (eight templates), with the frozen pool to be censused before inference and stragglers fixed then | — |
| D-093 | ~~Approval to run screening pass 1~~ **Both passes run (2026-09-23/24, $4.85 total); the 24 flags resolved (D-094).** **The Layer 2 protocol was decided (D-095) and run to completion (D-109).** | — |
| — | ~~Whether the full run adds the step router~~ **In the stack (D-110), batched or not at all (D-113).** ~~Its routing rule~~ **rule C (D-112).** ~~Measuring the batched prompt~~ **Measured: no detection lost, about $114 at full scale (D-113).** ~~Building it, its reported score~~ **Built and validated; Q3 reports its flags (D-142).** Still open: the approval of its run (above) | before evaluation starts |
| — | The repo's templates no longer reproduce 17 of the 60 frozen items byte-identically (D-102); `milestones.build_all` raises on the pilot manifest; `pinned_templates.py` works around it for the pilot | any re-derivation of milestones from templates |
| — | ~~Which families the next roster will evaluate~~ **Decided 2026-09-27 (D-110): the pricing document's eleven** | — |
| D-115 | ~~Round 4: the three chemical experts review `heat_of_reaction_formation` and `adiabatic_flame_temperature` as widened~~ **Done (D-119): both approved by all three; 150 of 150 certified** | — |
| D-118 | ~~Relabel `critical_depth_froude_classification`'s answer type~~ **Not needed (D-120): scoring is the same under either label** | — |
| D-121 | ~~`qwen3-235b-a22b` has no endpoint meeting the routing rule: keep it on Alibaba's endpoint (8,192-token cap, undeclared quantization) or replace it~~ **Replaced by `qwen3-235b-a22b-2507`; `qwen3.8-27b` set aside under D-110's rule (D-132)** | — |
| D-121 | ~~Approve the harness check (about $0.01) and the calibration run (about $3.58, 220 calls)~~ **Approved by the owner 2026-09-28 (D-122)** | — |
| D-122 | ~~Approve the full run: about $202 for the ten, 95% interval $163 to $251, on the calibration's bills (D-123)~~ **Approved for seven by the owner 2026-09-28 (D-124): $30.41, interval $26.67 to $34.66; the spend line moved to $45.06 after the stop (D-125)** | — |
| D-124 | ~~Whether and when to run `kimi-k3`, `glm-5.3` and `claude-sonnet-5`: $171.74, interval $134.07 to $221.32~~ **Approved by the owner and started 2026-09-28 (D-127)** | — |
| D-126 | ~~A private Kaggle copy of the traces~~ **Done 2026-09-29: the twelve models' archive as a private Kaggle dataset in the owner's account, downloaded back and matched file for file; the local archive is made and verified (D-131)** | — |
| D-117 | ~~Confirm the analysis plan as a whole; its two scoring rules are decided~~ **Confirmed by the owner 2026-09-28 (D-122)** | — |
| D-114 | ~~A private backup of `full_run_28092026/pool/` and `SEED.secret`~~ **Backed up 2026-09-28: a private Kaggle dataset in the owner's account, downloaded back and matched file for file, and a local archive with its checksum.** ~~The tag on the commit inference runs at~~ **`full-run-inference` (D-122)** | — |
| — | The budget: the plan comes to about $606 against the ~$500 round with a batched router (D-110). Calibration puts generation at about $202 for the ten against the plan's $403 (D-123); the seven's run cost $32.96 to $35.75 (D-126) | generation starts |
| — | ~~Expert annotation of the frozen 300 (stage 3)~~ **Done: 15 experts, every trace labelled three times, with verification and adjudication rounds (RESULTS_X1)** | — |
| D-003 | Do the raw `inference_results/` generations still exist? | promising any corrected results table |
| — | Phase 5 scoping: fold into Phase 1 or run as a parallel PR | Phase 1 start |
| — | Whether the 9 self-inconsistent templates are fixed or replaced | Phase 1 start (item-pool ownership call) |
