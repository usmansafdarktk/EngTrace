# Prompt 03 — Phase 2: determinism and solver-in-the-loop

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

This stage was **blocked on the constants track** and is now unblocked: Phase C2 landed the
corrected chemical thermochemistry constants (sync point S1, merge `2e9c8f2`). Phase 1
(round-trip integrity) is complete and merged.

---

I need you to implement **Phase 2** of the template redesign: making four templates
deterministic and replacing an opaque numerical solver with a stated algorithm. This is an
implementation task with a hard verification gate.

## What EngTrace is, and why this matters

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from
**150 parameterized Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each, in 47 files). Each template is
a function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of strings. The solution is the **gold reasoning trace**.

**The paper has been rejected twice from ACL ARR.** The most damaging criticism, raised
independently by five reviewers, is that the framework verifies model reasoning using an "AI
Tribunal" of three frontier LLMs while evaluating models from those same families. The paper
is titled *"Verifiable Process Supervision"*, but verification is a majority vote among LLMs.

The long-term fix is deterministic verification: templates emit a structured trace, and
checking a model's reasoning becomes a numeric and dimensional check with no LLM in the
critical path.

**Phase 2 is where "deterministic" stops being a figure of speech.** Phase 1 fixed traces that
contradicted their own answers. This phase fixes something more basic: four templates where
*the recorded seed does not regenerate the item*, or where the printed reasoning is not
something any reader could reproduce. A benchmark whose items cannot be regenerated from their
seeds cannot be audited, replicated, or corrected after publication — and one of these four is
already in published results.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy.
- `evaluation/` is a **separate track** and out of scope. Its parser has known defects
  (D-003). Do not fix them here.

## Read these first, in this order

1. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 governing principles **P1–P6** (especially **P4 determinism** and **P5 no silent
   fallbacks**), all of **Phase 2**, and the review protocol **R0–R6**.
2. **[`docs/re-implementation-sep/track-a/phase1_summary.md`](../re-implementation-sep/track-a/phase1_summary.md)** —
   how the previous implementation phase was run, gated and reviewed. §"What Phase 2
   inherits" is part of your brief.
3. **[`docs/re-implementation-sep/track-b/phaseC2_summary.md`](../re-implementation-sep/track-b/phaseC2_summary.md)** —
   the constants this phase now depends on. **§12 "What Phase 2 inherits"** and **D-032** are
   a required deliverable of this phase, not background.
4. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   append-only. **D-016** (why display ties are removed, not resolved), **D-024/D-026** (the
   evidence-rate rule), **D-032** (the Cp extrapolation trade you must decide), **D-034** (a
   test may not carry its own answer key).

## Scope — four templates, and nothing else

| Template | File |
|---|---|
| `template_levenspiel_plot_interpretation` | `chemical_engineering/reaction_kinetics/conversion_and_reactor_sizing.py` |
| `template_pfr_volume_changing_rate` | `chemical_engineering/reaction_kinetics/mole_balances.py` |
| `template_adiabatic_flame_temperature` | `chemical_engineering/thermodynamics/heat_effects.py` |
| `template_vibration_isolator_design` | `mechanical_engineering/vibrations_and_acoustics/harmonically_excited_vibrations.py` |

**Do not touch `data/templates/` outside these four.** If you find a defect elsewhere, record
it in `DECISIONS.md` or the phase summary — do not fix it.

**Do not modify `evaluation/`.**

**Do not modify the harness to make a fix pass.** If a check is wrong, that is a finding and a
`SPEC-CHANGE`, raised explicitly. The same applies to oracles: if an oracle looks wrong, that
is a finding to raise, not a file to edit.

**Do not push.** `master` is ahead of `origin/master`; publishing is a separate decision.

## The defects

Spec §2.1 lists these. The spec's list is not complete — Phase 1 and C2 both found the
recorded defect was one of several, and you should expect the same here. Verify each claim
before fixing it, and report what the spec missed.

**`levenspiel_plot_interpretation`** — unseeded `np.random.uniform` (the only T3 failure in
the whole corpus: the recorded seed cannot regenerate the item). Six blanket `try/except
Exception` that silently substitute a different formula while the trace still prints the
original one. Step numbering restarts `1,2,3` then `1,2,3,4,5,6`, which breaks T4.

**`pfr_volume_changing_rate`** — prints `scipy.integrate.quad`'s own error estimate as a step
value, which is machine-dependent. `n = round(uniform(1.5, 2.5), 1)` lands on exactly `2.0`,
contradicting the printed "non-integer order" note. Note that the integral
`∫₀ˣ(1−X)⁻ⁿdX` **has a closed form** for `n ≠ 1`; the spec says prefer it where one exists.

**`adiabatic_flame_temperature`** — `scipy.fsolve`, so the milestone is only
tolerance-reproducible. Blanket `except (…, Exception)`. A `while/else` fallback that returns
a tuple with no steps and no answer. Beyond the spec's list: **Step 4 says "we use a numerical
solver" and shows no arithmetic at all**, so the trace is not verifiable even in principle.

**`vibration_isolator_design`** — quadratic root *selection* (`max` of two roots) not
expressible as a formula. Step 4 **hardcodes** its verdict (`"Yes, … the condition is met"`)
without testing it. Two answer-less error returns. `{:,.0f}` emits thousands separators that
break the parser.

## Required changes

Spec §2.2 is authoritative. In brief:

1. **Seed all generators.** `np.random.seed()` alongside `random.seed()`, or replace numpy
   sampling with stdlib. Enforced by T3.
2. **Delete every blanket `except`.** Per P5, a failure raises. Where a fallback handles a
   genuinely rare sampling failure, replace it with a bounded resample loop that raises on
   exhaustion — the pattern already used throughout civil/industrial.
3. **Make the solver step reproducible.** Replace `fsolve` with a *stated* algorithm: fixed
   iteration count, stated initial guess, stated convergence rule, **all printed**. The trace
   must become something a student can reproduce by hand. Where a closed form exists, prefer
   it.
4. **Remove machine-dependent printed values.**
5. **Fix `n == 2.0`** — resample to exclude, or reword the note.
6. **Renumber Levenspiel's steps** contiguously (T4).
7. **`vibration_isolator_design`:** make root selection an explicit predicate step that
   *states and evaluates* the `r > √2` condition; replace both answer-less returns with
   raises; drop the `{:,.0f}` format.

**Before deleting any fallback, measure how often it currently fires.** Phase 1 found two
"dead" paths that were live on 17% of instances; deleting one silently would have changed the
item pool without anyone noticing. A fallback you remove is either provably unreachable
(show the argument) or it was changing answers (report the rate).

## The extra deliverable: D-032

Phase C2 handed this phase an open decision, and it is **not optional**.

`CP_PARAMS` is fitted to 1500 K. `adiabatic_flame_temperature` integrates to **2844 K**, where
the worst Cp error reaches **14.3%**. The C2.2 suite reports this on every run. A wide-range
refit fixes the high end (N₂: 12.5% → 0.7%) but costs the low end, where
`sensible_heat_temp_dependent_cp` lives — so it is a trade between two consuming templates,
and C2 deliberately left it to the template owner. **You are the template owner.**

Decide and record it: refit, restrict the sampled range, split the table by consumer, or
accept and document. State the P6 pedagogy cost either way. Whatever you choose, the flame
temperature must stay physically plausible — methane/air against the correct
no-dissociation reference of **2326.35 K** (`docs/references/README.md`; the spec's ~2200 K is
the *with-dissociation* figure and is the wrong target — D-029).

## Verification

The harness is `tests/template_integrity/`. Checks: **T1** printed-arithmetic closure, **T2**
round-trip oracle, **T3** determinism (two processes), **T4** output contract, **T5** value
binding and rounding, **T6** distribution non-regression, **T7** invariant asserts.

```
python -m tests.template_integrity.run                     # T1,T4,T5,T7 corpus-wide
python -m tests.template_integrity.run --checks all        # add T3 and T6
python -m tests.template_integrity.run --templates a,b
python -m tests.template_integrity.gate_report <N>         # exit-gate detail at N seeds
python -m tests.template_integrity.t6_report               # distribution profile
python -m tests.template_integrity.instance_dump           # (question, solution) per seed
```

**`instance_dump` must be run as separate processes per tree.** An in-process reload silently
reported "120 identical" once. Corpus baseline at master, for the non-regression comparison —
**verify it yourself against a `git worktree` rather than trusting this number**:

| T1 | T2 | T4 | T5 | T7 | failing |
|---|---|---|---|---|---|
| 32 | 0 | 4 | 71 | 87 | 85 |

*(Use a short worktree path such as `C:/wtm` — the repo has paths long enough to fail a
checkout under a temp directory.)*

**State the smallest rate your run can resolve, and check it against the rate you need to
exclude.** N seeds cannot resolve a defect rarer than about 3/N. This rule exists because
acceptance evidence was smaller than the defect rate three separate times (D-024, D-026).

**T6 is expected to move.** Removing silent fallbacks changes the accepted sample set. Movement
must be measured, explained and approved (D2.5) — not waved through, and not suppressed.

## Deliverables

D2.1 four templates edited · D2.2 T3 determinism proof, 200 seeds × 2 processes · D2.3 written
justification for each removed `try/except` — what it was masking, why removal is safe, and
**how often it fired** · D2.4 stated-algorithm spec for the replaced solver, with a
hand-worked example · D2.5 T6 distribution diff with explanation · D2.6 confirmation that
`CP_PARAMS` landed first · **D2.7 the D-032 decision, with its P6 trade stated**

Write the phase summary to `docs/re-implementation-sep/track-a/phase2_summary.md` and the item-pool
impact to `phase2_item_pool_impact.md`, following the Phase 1 and C2 documents' shape.

## Independent review — mandatory, two reviewers

Scope each to **R6** before dispatching: **one** mandatory gate task, a stated time box,
tooling supplied, settled things fenced off, and **no measurement commissioned twice**. An
over-scoped review stalls and leaves its gate unchecked.

*Reviewer A — Determinism.* Re-run T3 on a different Python version. Grep for every
non-deterministic source: unseeded RNG, `set`/`dict` iteration order affecting output, `time`,
`id()`, floating-point-dependent branching. Confirm no blanket excepts remain.

*Reviewer C — Numerical methods.* Confirm the stated algorithm converges over the **full**
sampled parameter range, that the stated iteration count suffices at the extremes, and that
printed intermediates are hand-reproducible. Verify the flame temperatures against an external
source.

**Give each reviewer a frozen ref.** In C2 the branch moved while a reviewer worked, and three
of its twelve findings were stale on arrival — about a quarter of the review wasted. That was
the implementer's error, not the reviewer's.

File each report to `docs/re-implementation-sep/reviews/`, triage every finding **and** every
§5 suggestion in the summary (R4), and state honestly which findings you rejected and why.

## Exit gate

- [ ] T3 passes across two processes and two Python versions
- [ ] Zero blanket `except` in scope; zero answer-less return paths
- [ ] T4 contract conformance (Levenspiel numbering)
- [ ] Solver replacement converges across the full sampled range
- [ ] `adiabatic_flame_temperature` physically plausible against 2326.35 K
- [ ] No corpus regression against the master baseline, measured not assumed
- [ ] D-032 decided and recorded
- [ ] Reviewers A and C filed, every finding and suggestion triaged (R2, R4)

Work on a branch. Merge to `master` with `--no-ff` when the gate is green.
