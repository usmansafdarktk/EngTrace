# Prompt 04 — Phase 3: trace shape, iteration and search

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Phases 0, 1, C2 and 2 are complete and merged (`3cb147c`). This stage depends on all of
them, and on Phase 2 in particular — it inherits a criterion that Phase 2's own review
corrected.

---

I need you to implement **Phase 3** of the template redesign: the two templates whose traces
are *shaped* by iteration and search. **This is the one phase where the right answer may be to
change the schema rather than the template**, and where the primary deliverable is a
specification rather than code.

## What EngTrace is, and why this matters

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from
**150 parameterized Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each, in 47 files). Each template is
a function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of natural-language strings. The solution is the **gold
reasoning trace**.

**The paper has been rejected twice from ACL ARR** (January 2026 and May 2026). The most
damaging criticism, raised independently by five reviewers, is that the evaluation framework
verifies model reasoning using an "AI Tribunal" of three frontier LLMs (GPT-5, Claude Opus 4.5,
Gemini 3) while simultaneously evaluating models from those same families — GPT-5 is both judge
and evaluated model. The paper is titled *"Verifiable Process Supervision"*, but verification is
a majority vote among LLMs.

The long-term fix is deterministic verification: templates emit a **structured trace**, and
checking a model's reasoning becomes a numeric and dimensional check with no LLM in the
critical path. Phases 1, C2 and 2 made the gold traces trustworthy — they now reproduce their
own answers, rest on constants verified against NIST, and regenerate from their recorded seeds.

**Phase 3 is where the trace stops being a flat list of steps.** Every other template in the
corpus emits a fixed sequence of milestones. These two do not: one runs a secant solve whose
iteration count varies with the data, and the other runs a greedy heuristic whose *number of
workstations is part of the answer*. A milestone model with a fixed `{id, symbol}` per step
cannot represent either. That is the gap this phase closes, and **D3.3 feeds directly into the
milestone-model design** — it is the deliverable everything after Phase 3 is built on.

For the full picture (optional, but read it if a judgement call turns on what the benchmark is
claiming):

- `docs/_ARR_May__EngTrace.txt` — the current paper. **Section 4** is the evaluation framework,
  **Appendix G** has three full template implementations, **Appendix M** is the framework
  validation.
- `docs/EngTrace_Rebuttal_Jul2026.txt` — the most recent reviewer objections and our responses.
  **These promises are the real to-do list behind this work**, and a structured trace for
  iterative items is one of them.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not
  flag it.
- `evaluation/` is a **separate track** and out of scope here. Its parser has known defects
  (it reads the gold answer `**4,921**` as `4.0`), recorded as D-003. Do not fix them.

## The corpus these two come from, and why it matters

**There are two authoring eras in this repo, and both of your templates are from the good
one.** Civil and industrial were written later under a spec-driven process: provenance-tagged
constants, textbook citations, stated physical bounds, invariant `assert`s, and a
round-then-recompute convention. Chemical, electrical and mechanical came earlier with none of
that.

`normal_depth_iteration` is civil. `line_balancing_heuristic` is industrial. **Each carries 10
`assert` statements, textbook grounding (Sturm's *Open Channel Hydraulics* Ch. 4), stated
physical bounds, and a bounded resample loop.** They are among the best-authored templates in
the corpus.

That reframes the phase, and you should hold onto it: **the defects here are not carelessness.**
They are structural consequences of a trace whose shape depends on its data, in templates whose
authors did everything else right. You will not fix them by tightening discipline, because the
discipline is already there. Do not go looking for sloppiness — look at the shape.

Corpus health at `master`, for context (this is the baseline you work against, not a to-do
list):

| Check | Failing | Note |
|---|---:|---|
| T1 printed-arithmetic closure | 30/150 | 97 more sit on the rounding boundary (marginals) |
| T2 round-trip oracle | 0/150 | only templates with an oracle are checked |
| T3 determinism | 0/150 | was 1 until Phase 2 |
| T4 output contract | 3/150 | malformed step markers |
| T5 binding / rounding | 68/150 | |
| T6 distribution | 142/150 | **the committed baseline is stale corpus-wide — see Verification** |
| T7 invariant asserts | 83/150 | advisory corpus-wide; a gate on templates a phase edits |

## Read these first, in this order

1. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 governing principles **P1–P6**, all of **Phase 3**, and the review protocol **R0–R6**.
   Read R0–R6 properly before you dispatch anything; it is short and it is the part most often
   skipped.
2. **[`docs/re-implementation-sep/track-a/phase2_summary.md`](../re-implementation-sep/track-a/phase2_summary.md)** —
   **§10 is written for you, and §12 corrects it.** Read both, and believe §12.
3. **[`docs/re-implementation-sep/reviews/phase2_reviewer_c_numerical.md`](../re-implementation-sep/reviews/phase2_reviewer_c_numerical.md)** —
   §5 is the sharpest statement of what makes these two templates hard, and it is not what the
   spec says.
4. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   append-only; add entries, never edit existing ones. **D-016** and **D-037** (display ties:
   removed, and the two ways to remove them), **D-024/D-026** (the evidence-rate rule),
   **D-031** (a dict's key order can be part of the item pool's identity), **D-034** (a test
   may not carry its own answer key).
5. **[`docs/re-implementation-sep/track-a/phase1_summary.md`](../re-implementation-sep/track-a/phase1_summary.md)**
   and **[`phaseC2_summary.md`](../re-implementation-sep/track-b/phaseC2_summary.md)** — the shape your
   own close-out should take, and the R4 triage tables to model.

## Scope — two templates, and nothing else

| Template | File |
|---|---|
| `template_normal_depth_iteration` | `civil_engineering/water_resources/uniform_flow.py:287` |
| `template_line_balancing_heuristic` | `industrial_engineering/production_and_inventory/production_planning.py:158` |

**Do not touch `data/templates/` outside these two.** A defect found elsewhere is recorded in
`DECISIONS.md` or the phase summary — not fixed.

**Do not modify the harness or an oracle to make a fix pass.** If a check looks wrong, that is
a finding and a `SPEC-CHANGE`, raised explicitly. Phase 0 found two checks that were green
while measuring nothing; treat harness changes with suspicion.

## What is actually wrong, measured

Verify these yourself before acting — **every prior phase found the spec's defect list was a
starting point rather than an inventory.** Phase 1 found its prescribed fix was insufficient,
Phase C2 found one recorded defect was four, and Phase 2 found a template failing T1 on 91% of
instances that the spec did not mention at all.

**`template_normal_depth_iteration`** — a secant solve on Manning's equation. Measured over
400 seeds: the **step count is constant at 4**, but Step 3 carries **1 to 5 interpolation
updates** (distribution `{1: 5, 2: 31, 3: 230, 4: 128, 5: 6}`), each with its own `A`, `P`,
`A*R^(2/3)`, `g`. So there is no stable `{id, symbol}` for "the third trial depth" —
`y_curr`/`g_curr` are rebound every pass. Also:

- `g_curr - g_prev` is an **unguarded division**. Reachable in principle; 0 occurrences in
  6,000 seeds. Measure, then decide — do not assert it away.
- `fmt3()` rewrites `"-0.000"` to `"0.000"`, so **the printed operand can differ in sign from
  the stored value.** A direct P2 violation. T5 flags it three times as
  `[result-inline] fmt3(g_prev)` / `fmt3(g_curr)`.
- T1 has **1 hard failure in 200 seeds** (seed 11, the `K = Q*n/S^(1/2)` line).

**`template_line_balancing_heuristic`** — a greedy longest-eligible-task rule. The trace is a
**search log**: milestones are task *sets* (`Station 1 contents: {a, b}`), a running remaining
time, and boolean fit tests. Measured over 400 seeds: **the step count itself varies, 5 or 6**,
because the station count is 3 or 4 — and **`n` is a computed result, not an input.** Balance
delay is a function of `n`, so fixing the station count to make the trace a fixed shape
**gives away part of the answer.** T1 has 4 hard failures in 200 seeds; T5 flags 2
`[result-inline]` expressions.

## The decision you are here to make (D3.1)

For each template, **independently**: schema route or redesign?

The spec recommends the schema route for both — an `iteration` node type and a `decision` node
type, which it observes are the same shape (an ordered list of homogeneous sub-traces with
stable within-element symbols and a termination predicate). It offers a cheaper alternative for
§3.1: fix the iteration count at 3 with stable names, class A for ~4 h, at the cost of pinning
the pedagogy to "do exactly three secant updates."

**Phase 2 supplies a worked precedent for the cheap route, and its review supplies the limit.**
`template_adiabatic_flame_temperature` *was* an iterative solve, and Phase 2 resolved it with no
schema change — fixed 6 passes, every pass printed, answer to the nearest kelvin. Phase 2's
summary §10 offered you the criterion:

> *fix the iteration count when the count is not part of what the item tests.*

**Reviewer C showed that is necessary but not sufficient. Treat §12 as authoritative over §10.**
The flame template is verifiable only because its parameter space is **a closed set of eleven
discrete reactions** — the rounding-boundary question could be checked *exhaustively*. Normal
depth samples channel geometry, slope and roughness **continuously**, so a fixed iteration count
cannot be exhaustively validated, and there will be a measure-zero-but-nonempty set of
parameters where the fixed count lands on the wrong side of a display boundary. Taking the cheap
route there needs **either a proven a-priori error bound over the sampled box, or a run-time
convergence assert that is not stripped by `python -O`** — note that `assert` is, which is why
the flame template's existing one guards development only.

For `line_balancing_heuristic` the reasoning does not transfer *at all*: **a heuristic has no
fixed point**, so "N passes suffices" is not even the right shape of argument. The question
there is whether the heuristic is deterministic under tie-breaking — check what happens when two
eligible tasks have equal duration, and whether `_T19_ORDER` or a `set` iteration order reaches
the output.

**Two templates, two decisions, argued separately. Do not pick one route for both because they
look alike.**

## Deliverables

1. **D3.1** — decision record per template, with the **P6 trade stated explicitly**.
2. **D3.2** — template edits.
3. **D3.3** — the **`iteration` and `decision` node-type specifications**, with worked examples
   from both templates. **This is the primary deliverable.**
4. **D3.4** — comparator semantics: how a model's iteration count differing from gold is scored.
   **A correct answer reached in 4 iterations instead of 3 must not be marked wrong.**
5. **D3.5** — T6 distribution diff (see the wall below).
6. **Item-pool impact note** — `phase3_item_pool_impact.md`: which published results are
   invalidated by regenerating these templates, and by how much.
7. **`phase3_summary.md`** — close-out in the shape of `phase1_summary.md` / `phase2_summary.md`:
   what changed, the R4 triage tables, exit-gate status, residual risk register, and what
   Phase 4 inherits.
8. **New entries in `DECISIONS.md` for every decision or reversal you make** — including the
   ones you reverse mid-phase. Phase 1's most valuable record is the one where its own
   prescribed fix turned out to be wrong.
9. **Two review reports** under `docs/re-implementation-sep/reviews/`.

## Verification

```
python -m tests.template_integrity.run --checks all          # corpus, ~60s
python -m tests.template_integrity.run --templates a,b --checks all
python -m tests.template_integrity.gate_report <N>
python -m tests.constants_integrity.test_chemical_thermochemistry
python -m tests.constants_integrity.test_citations_resolve
```

**Verify the corpus baseline yourself against a `git worktree` at `master` rather than trusting
the table above** — Phase 2 caught a stale recalled number this way. Use a short worktree path
such as `C:/wtm`; the repo contains paths long enough to fail a checkout under a temp directory.

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and will
crash on `→`/`Σ`/`²`.

**The T6 wall.** T6 fails 142/150 **on `master`**, because the committed baseline is stale
corpus-wide. D3.5 will hit this exactly as Phase 2 did. Either deliver the distribution diff as
a direct before/after instance dump (Phase 2's approach — see `phase2_item_pool_impact.md` §3),
or regenerate the baseline **as its own deliberate change on `master`**, with the corpus-wide
movement examined rather than absorbed. A baseline refreshed by the phase it is meant to gate is
not a gate. If you regenerate it, that decision needs its own record.

**State the smallest rate your run can resolve and check it against the rate you need to
exclude.** N seeds cannot resolve a defect rarer than ~3/N. This rule exists because acceptance
evidence was smaller than the defect rate three separate times (D-024, D-026).

**Instance dumps must run as separate processes per tree.** An in-process reload resolves both
sides to the same already-imported modules and reports everything identical — documented in
`tests/template_integrity/instance_dump.py`, and it has caught two people.

## Reviews — read R0–R6 before dispatching any

Phase 3 requires **two independent reviewers, working in parallel and in isolation from each
other.**

- **Reviewer B — pedagogy (the P6 guard).** Does either item still test what it tested?
  Specifically: **does the `line_balancing` result still require the solver to *run* the
  heuristic, or has it become a lookup?** Also: would a domain expert notice either item got
  easier or less interesting?
- **Reviewer D — schema.** **Given only the D3.3 node specs and a set of generated traces,
  write a working verifier.** Report every ambiguity. **The spec passes only if a reviewer who
  did not design it can implement against it.** This is the real test of D3.3 and the
  deliverable most likely to fail its gate on first attempt — budget for a second round rather
  than treating one failure as the end.

**Scope every review to R6 before dispatching it.** One mandatory gate task, a stated time box,
tooling supplied as working code, everything already settled fenced off explicitly, and a
measurement-ownership register with **no number commissioned twice**. A Phase 0 review stalled
and produced nothing because its brief duplicated another agent's. Do not ask a reviewer to
re-measure something a prior phase already established — say so explicitly in the brief.

**Independence, concretely:**

- **R1.1 — work without sight of the implementer's reasoning.** Supply code comments, commit
  messages and summaries as **claims under test**, not as background. Phase C2's Reviewer G was
  fenced off from values entirely and still found the phase's deepest defect *because* it was
  asked whether the evidence was real rather than whether the numbers were right.
- **Give each reviewer a frozen ref — a commit SHA, not a branch name.** In Phase C2 the branch
  moved under a reviewer and three of its twelve findings were stale on arrival, roughly a
  quarter of the review wasted. Phase 2 fixed this and it worked; Reviewer A's report says so.
- **Do not let the two reviewers see each other's briefs or reports.** Scope them
  orthogonally — B asks whether the item still teaches, D asks whether the spec is
  implementable. Neither should be re-deriving the other's numbers.

**Reports go to `docs/re-implementation-sep/reviews/phase3_<reviewer>.md` in the R2 five-section
structure** (verdict · gate result · findings · falsification attempts that failed · further
probing). **Every §5 suggestion must be triaged (R4) before the phase closes; an untriaged
suggestion blocks the gate exactly as a CONFIRMED finding does.** State honestly which findings
you rejected and why.

`docs/re-implementation-sep/reviews/` is now tracked — a bare `reviews` entry in `.gitignore`
had silently untracked all eight prior reports until Phase 2 caught it. **Confirm your reports
are actually committed, not merely written**; an ignored file shows as neither modified nor
untracked.

## Git

- Branch **`redesign/phase3-trace-shape`** off `master`. Do not work directly on `master`.
- Logical commits, not one lump — the node-type spec, the two templates, and the reviews are
  naturally separate.
- **Do not push.** `master` is 50 commits ahead of `origin/master`; publishing is a separate
  decision.
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
- Merge to `master` with `--no-ff` only after the exit gate passes, or record explicitly what
  is left open.

## Exit gate

- [ ] D3.1 decisions recorded, **argued separately per template**, P6 trade stated
- [ ] T1–T7 pass on both templates
- [ ] `g_curr - g_prev` guarded; the `fmt3()` sign discrepancy removed
- [ ] **Reviewer D implemented a working verifier from D3.3 alone**
- [ ] D3.4 comparator semantics cover the iteration-count-mismatch case
- [ ] No corpus regression against the master baseline, measured not assumed
- [ ] Reviewers B and D filed, **every finding and every §5 suggestion triaged (R2, R4)**
- [ ] Item-pool impact note filed

## Constraints and cautions

- **P6 is a real constraint.** These are benchmark items. If a fix would change what an item
  tests or how hard it is, stop and record the trade rather than absorbing it. For
  `line_balancing` this is the whole phase: the station count is part of the answer.
- **Where you are uncertain, say so and explain what evidence would settle it.** Do not smooth
  over a gap to make the gate look clean. Phase 0's most valuable finding was that its own
  recommended fix was insufficient, and it surfaced only because it was measured rather than
  assumed. Phase 2's most valuable finding was a 0.5% error in an answer key that its own
  impact measurement had already seen and dismissed as rounding noise.
- **Report your own errors in the summary.** Every phase so far has recorded at least one — a
  10× transcription slip caught by its own test suite, a defect introduced and caught by T1
  within a minute, an evidence script that was 135 K wrong and nearly drove a decision the
  wrong way. That record is worth more than a clean-looking close-out.

## Open items inherited from Phase 2 — carry forward, fix only if they block your gate

- **Two rounding conventions coexist in the corpus** — `_hu` (Decimal `ROUND_HALF_UP`, D-012)
  in the vibrations file, `_as_printed`/`format` (binary `ROUND_HALF_EVEN`) in the chemical
  files — and only some carry a display-tie guard. Reviewer A recommends a corpus-wide sweep of
  every template combining a transcendental with a fixed-decimal print. **`normal_depth_iteration`
  is exactly that shape** (`(A/P)**(2/3)` printed at 3 dp), so this may land in your lap whether
  or not you go looking for it.
- **Reviewer C's two generalisable defects:** a raw arithmetic expression interpolated into an
  f-string (`{n-1}` printing `0.6000000000000001`), findable corpus-wide with a regex for
  `\{[a-z_]+ *[-+*/] *[0-9.]+\}`; and fixed *decimal* places applied to a quantity spanning
  decades, which cost 0.5% of an answer key in Phase 2.
- **Fallback firing-rate evidence lives in scratch scripts, not the repo**, so a reviewer cannot
  check it. If you remove a fallback, consider committing the measurement.
