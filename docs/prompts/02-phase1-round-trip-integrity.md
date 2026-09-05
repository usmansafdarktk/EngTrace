# Prompt 02 — Phase 1: round-trip integrity

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

This is the **first stage that modifies templates.** Phase 0 (the verification harness and
baseline) is complete and merged; this stage depends on it.

---

I need you to implement **Phase 1** of the template redesign: fixing 12 templates whose gold
traces do not reproduce their own answers. This is an implementation task with a hard
verification gate, not a design task — the design decisions were made and are recorded.

## What EngTrace is, and why this matters

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from
**150 parameterized Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each, in 47 files). Each template is
a function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of natural-language strings. The solution is the **gold
reasoning trace**.

**The paper has been rejected twice from ACL ARR** (Jan 2026 and May 2026). The most damaging
criticism, raised independently by five reviewers, is that the evaluation framework verifies
model reasoning using an "AI Tribunal" of three frontier LLMs (GPT-5, Claude Opus 4.5, Gemini
3) while simultaneously evaluating models from those same families — GPT-5 is both judge and
evaluated model. The paper is titled *"Verifiable Process Supervision"*, but verification is a
majority vote among LLMs.

The long-term fix is deterministic verification: templates emit a structured trace, and
checking a model's reasoning becomes a numeric and dimensional check with no LLM in the
critical path. **Phase 1 is a prerequisite for that.** You cannot build a deterministic
verifier against gold traces that contradict themselves — any tolerance tight enough to catch
a real reasoning error would reject the gold. That is what these 12 templates do today.

For the full picture (optional, but read it if a judgement call turns on what the benchmark is
claiming):

- `docs/_ARR_May__EngTrace.txt` — the current paper. **Section 4** is the evaluation
  framework, **Appendix G** has three full template implementations, **Appendix M** is the
  framework validation.
- `docs/EngTrace_Rebuttal_Jul2026.txt` — the most recent reviewer objections and our
  responses. These promises are the real to-do list behind this work.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not
  flag it.
- `evaluation/` is a **separate track** and out of scope here. Its parser has known defects
  (it reads the gold answer `**4,921**` as `4.0`), recorded as D-003. Do not fix them in this
  phase.

## The current state of the templates

**Structure.** Every template follows one shape: sample parameters (often inside a rejection
loop) → compute into local variables → concatenate one large f-string → return. The solution
string is a `**Given:**` block, then `**Step N:**` blocks, then a final `**Answer:**`. Step
counts run 2–9, median 4.

**Health, as measured by Phase 0** across all 150 (this is the baseline you are working
against, not a to-do list):

| Check | Failing | Note |
|---|---:|---|
| T1 printed-arithmetic closure | 38/150 | 95 more sit on the rounding boundary (marginals) |
| T3 determinism | 1/150 | `levenspiel_plot_interpretation`, unseeded `np.random` |
| T4 output contract | 4/150 | malformed step markers; 5 more use `**Final Answer**` |
| T5 binding / rounding | 73/150 | 39 are confirmed P2 violations |
| T7 invariant asserts | 95/150 | expected — see below |

All 150 import and run cleanly: **zero exceptions across 30,000 generations.** The corpus is
healthy; the defects are specific, not systemic.

**Two authoring eras, and it matters.** Civil and industrial were written later under a
spec-driven process: provenance-tagged constants, textbook citations, stated physical bounds,
**60 templates carrying 142 `assert` statements between them, and a "round-then-recompute"
convention**. Chemical, electrical and mechanical have **zero asserts** and no provenance.

That is why T7 fails 95/150 — it is advisory corpus-wide and a gate only on templates a phase
edits. It is also why the civil convention is the model P2 is built on. But note the trap
Phase 0 found: civil holds round-then-recompute in 29/30 templates and its **exemplar still
fails closure on 5.89% of instances**, because the convention as written is incomplete. Do not
assume a civil template is correct because it looks disciplined.

**Useful structural facts** (from the audit): 93% of f-string interpolations are already bound
variables rather than inline computations; 48% of templates change their governing equation,
step count or milestone set based on a sampled parameter, so **several of the 12 have multiple
structural branches you must exercise** — generate enough seeds to hit each one.

## Read these first, in this order

1. **[`docs/re-implementation-sep/phase0_summary.md`](../re-implementation-sep/phase0_summary.md)** —
   start here. §8 "What Phase 1 inherits" is your brief in five points. §3 explains why the
   original fix recipe was wrong.
2. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   the plan. Read §0 governing principles (**P1–P6**, especially the amended **P2**), all of
   **Phase 1**, and the review protocol **R0–R6**. Note that **§1.1a is authoritative** and
   §1.1 is retained only for provenance.
3. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   D-010 through D-013 are the decisions that set this phase's scope and fix recipes.
   Append-only: add entries, never edit existing ones.
4. **[`docs/re-implementation-sep/reviews/phase0_oracle_findings.md`](../re-implementation-sep/reviews/phase0_oracle_findings.md)** —
   per-template measurements, the tolerance rationale, and §5's list of question-text
   fragilities you will run into.
5. **[`docs/re-implementation-sep/phase0_baseline.md`](../re-implementation-sep/phase0_baseline.md)** —
   the independent re-measurement. Consult for specific numbers; **do not re-measure what is
   already in here.**

Background on the whole effort: `docs/re-implementation-sep/template_audit_report.md`.

## Why this phase exists

These templates produce **incorrect gold traces that are already in published results**. A
solver who reads the question, follows the printed steps and uses the printed operands does
not arrive at the printed answer. Measured worst cases: `mean_variance` 56.8% of instances,
`rotating_unbalance` 33.7% of answers not the correct rounding of what the question implies,
`annulus_flowrate` 94% of Step-5 products not closing.

This is not cosmetic. Any verification tolerance tight enough to catch a real reasoning error
will reject these gold traces.

## Scope — 12 templates, three categories, different fixes

**Do not treat these as one problem.** The categories need different remedies, and applying
the wrong one is how Phase 0's original recipe would have failed.

### Category A — chain breaks (3). Fix with **P3**.

`template_mean_variance` · `template_rotating_unbalance` · `template_vibration_transmissibility`

The template samples a hidden target, derives the "given", states it **rounded**, then
computes everything downstream from the **unrounded** pre-image. Fix: round the stated given
first, then recompute forward from the rounded value. The sampler may still *choose* the
hidden target to place the item in a useful regime — it may not leak it into the solution.

### Category B — display defects (6). Fix with **P2 as amended**.

`template_beam_deflection_formula` ← **do this one first** · `template_cantilever_double_integration` ·
`template_annulus_flowrate` · `template_statically_indeterminate_shaft` ·
`template_shaft_design_power` · `template_composite_shafts_series`

The answer round-trips correctly; an intermediate printed line does not close. Rounding to
display precision is **necessary but not sufficient** — it removes the double rounding but
leaves an exact half-way tie that `round()` resolves on the binary float. Break the tie in
decimal (see P2 in the spec for worked code).

**Fix `beam_deflection_formula` first and get it reviewed before touching the other five** —
it is the pattern they copy, and it is the template whose apparent fix misled the original
spec.

### Category C — ill-posed, not wrong (1).

`template_damping_classification` — the gold label is correct; the item is ill-posed. A third
of instances sit on a measure-zero boundary and the 2-dp rounding of `c` destroys the digit
that would settle it. Fix: **state `c` to more digits** (or state it as equal to the critical
value) and drop the fragile `elif zeta == 1` float-equality test. Do **not** constrain the
sampling — that approach was superseded (D-004 → D-010).

### Explicitly out of scope — do not "fix" these

`template_poissons_ratio`, `template_logarithmic_decrement`, `template_system_properties`
were in the original Phase 1 and left it: independent oracles found them **not defective**
(D-013). `logarithmic_decrement`'s input back-solving is deliberate — it is the
measurement-uncertainty framing of the exercise. If you believe one of these is defective,
say so and stop; do not change it.

## Three defects to fix alongside

Found during Phase 0, not in the original brief:

- `mean_variance` — the probabilities are **never renormalised** (they sum to 0.999–1.002).
- `rotating_unbalance` — the question never states whether the unbalance mass is included in
  the stated total mass. The solution assumes it is; a solver who subtracts it is marked wrong.
- `rotating_unbalance` — some answers print to **one significant figure** (`0.002 mm`),
  ungradeable regardless of the chain fix.

## How to verify — the harness already exists

`tests/template_integrity/` runs all 150 templates in ~17s. Do not write new measurement code
before checking whether a check already covers it.

```bash
python -m tests.template_integrity.run --checks T1,T2,T4,T5,T7            # 9s
python -m tests.template_integrity.run --checks T3 --det-seeds 200        # 8s
python -m tests.template_integrity.run --checks T6                        # vs baseline
python -m tests.template_integrity.run --templates template_rotating_unbalance --checks all
python -m tests.template_integrity.run --checks T2 --seeds 500            # the 9 oracles
```

**Verify fixes against T5b, not T1 alone.** T1's ±0.5-ulp tolerance rule structurally cannot
see the half-way-tie class — that is the whole reason T5b exists (adversary finding F3). A
Category B fix that turns T1 green while T5b stays red is not fixed.

**Nine round-trip oracles already exist** in `tests/template_integrity/oracles/`, written
independently of the templates' computation code. They are your acceptance test: a Category A
fix must take its oracle from failing to passing. **Do not modify an oracle to make a fix
pass** — if an oracle looks wrong, that is a finding to raise, not a file to edit.

**T6 is the P6 guard.** Every fix must keep the distribution within tolerance: distinct
answers must not fall by more than 10%, branch proportions stay within ±5 points, median
answer magnitude within ±2%. A breach needs explicit sign-off, **not a widened tolerance**.

**T7:** every template you touch must carry ≥2 physical invariant asserts in the
civil/industrial style before the phase closes.

## Reviews — read R0–R6 before dispatching any

Phase 1 requires **two independent reviewers working in parallel and in isolation from each
other**:

- **Reviewer A — correctness.** Verifies P1–P3 hold, samples ≥200 seeds per template
  independently, and confirms each T2 oracle is a genuine independent re-derivation rather
  than a transcription of the template's arithmetic. **This is the single most likely way
  this phase passes while still broken**, and it is Reviewer A's primary task.
- **Reviewer B — physics and pedagogy (the P6 guard).** Confirms each item still tests what
  it tested: governing equations unchanged, parameter ranges still sensible, difficulty label
  still justified. Explicitly asked: *would a domain expert notice this item got easier or
  less interesting?*

**Scope every review to R6 before dispatching it.** One mandatory gate task, a stated time
box, tooling supplied as working code, everything already settled fenced off explicitly, and
a measurement ownership register with **no number commissioned twice**. A Phase 0 review
stalled and produced nothing because its brief duplicated another agent's — do not repeat
that. In particular: **do not ask a reviewer to re-measure defect rates**; they are in
`phase0_baseline.md`.

Reports go to `docs/re-implementation-sep/reviews/phase1_<reviewer>.md` in the R2 five-section
structure. Every §5 suggestion must be triaged (R4) before the phase closes; an untriaged
suggestion blocks the gate exactly as a CONFIRMED finding does.

## Git

- Branch **`redesign/phase1-round-trip`** off `master`. Do not work directly on `master`.
- Logical commits, not one lump — I suggest one per category, plus the alongside-defects.
- **Do not push.** `master` is ahead of `origin/master`; publishing is a separate decision.
- Do not modify `evaluation/` — the parser fix is a separate track (D-003).
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`
- Merge to `master` only after the exit gate passes, or record explicitly what is open.

## Deliverables

1. **12 templates edited**, satisfying P1–P3 and passing their checks.
2. **T2 oracles for the 3 newly-added templates** (`beam_deflection_formula`,
   `statically_indeterminate_shaft`, `shaft_design_power`) — written from the governing
   equations with the computation code withheld, per the contract in
   `tests/template_integrity/checks/t2_roundtrip.py`.
3. **A per-template `TOLERANCE` declaration** with a justification comment, for the verifier
   to consume later.
4. **T6 before/after distribution diff** for all 12.
5. **Item-pool impact note** — which published results are invalidated by regenerating these
   templates.
6. **`docs/re-implementation-sep/phase1_summary.md`** — close-out in the shape of
   `phase0_summary.md`: defect rate before → after (target 0), R4 triage table, exit-gate
   status, and what Phase 2 inherits.
7. **New entries in `DECISIONS.md`** for every decision or reversal you make.
8. **Two review reports** under `reviews/`.

## Exit gate

- [ ] T1 closure: no non-closing printed line at 1,000 seeds, all 12
- [ ] T2 round-trip passes at 1,000 seeds for every template with an oracle
- [ ] T5b: zero confirmed rounding violations in the 12
- [ ] T3 determinism passes
- [ ] T6 within tolerance, or breach explicitly signed off
- [ ] T7 ≥2 invariant asserts each
- [ ] Reviewer A and Reviewer B both clear; every §5 suggestion triaged
- [ ] Item-pool impact note filed

## Constraints and cautions

- **Do not touch `data/templates/` outside the 12 templates in scope.** If you find a defect
  elsewhere, record it in `DECISIONS.md` or the phase summary — do not fix it.
- **Do not modify the harness to make a fix pass.** If a check is wrong, that is a finding
  and a `SPEC-CHANGE`, raised explicitly. Phase 0 found two checks that were green while
  measuring nothing; treat harness changes with suspicion.
- **P6 is a real constraint.** These are benchmark items. If a correctness fix would change
  what an item tests or how hard it is, stop and record the trade rather than absorbing it.
- Several templates have fragile question text (split scientific notation, signs carried only
  in prose, duplicated noun phrases, inconsistent thousands separators). See
  `reviews/phase0_oracle_findings.md` §5 — anything that makes an oracle's parser fragile
  makes a model's reading fragile too, so these are worth fixing where cheap, but record them
  either way.
- Where you are uncertain, say so and explain what evidence would settle it. Do not smooth
  over a gap to make the gate look clean — Phase 0's most valuable finding was that its own
  recommended fix was insufficient, and it only surfaced because it was measured rather than
  assumed.
