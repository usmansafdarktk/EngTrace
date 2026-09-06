# Prompt 04 — Phase 3: trace shape, iteration and search

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Phases 0, 1, C2 and 2 are complete and merged (`103fa8b`). This stage depends on all of
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
returns a `(question, solution)` pair of strings. The solution is the **gold reasoning trace**.

**The paper has been rejected twice from ACL ARR.** The most damaging criticism, raised
independently by five reviewers, is that the framework verifies model reasoning using an "AI
Tribunal" of three frontier LLMs while evaluating models from those same families. The paper
is titled *"Verifiable Process Supervision"*, but verification is a majority vote among LLMs.

The fix is deterministic verification: templates emit a **structured trace**, and checking a
model's reasoning becomes a numeric and dimensional check with no LLM in the critical path.
Phases 1, C2 and 2 made the gold traces trustworthy — they now reproduce their own answers,
rest on correct constants, and regenerate from their seeds.

**Phase 3 is where the trace stops being a flat list of steps.** Every other template in the
corpus emits a fixed sequence of milestones. These two do not: one runs a secant solve whose
iteration count varies with the data, and the other runs a greedy heuristic whose *number of
workstations is part of the answer*. A milestone model with a fixed `{id, symbol}` per step
cannot represent either. That is the gap this phase exists to close, and **D3.3 feeds directly
into the milestone-model design** — it is the deliverable everything after Phase 3 is built on.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy.
- `evaluation/` is a **separate track** and out of scope. Its parser has known defects
  (D-003). Do not fix them here.

## Read these first, in this order

1. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 governing principles **P1–P6**, all of **Phase 3**, and the review protocol **R0–R6**.
2. **[`docs/re-implementation-sep/phase2_summary.md`](../re-implementation-sep/phase2_summary.md)** —
   **§10 is written for you**, and **§12 corrects it.** Read both, and believe §12. Details
   below.
3. **[`docs/re-implementation-sep/reviews/phase2_reviewer_c_numerical.md`](../re-implementation-sep/reviews/phase2_reviewer_c_numerical.md)** —
   §5 is the sharpest statement of what makes these two templates hard, and it is not what the
   spec says.
4. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   append-only. **D-016** and **D-037** (display ties: removed, and the two ways to remove
   them), **D-024/D-026** (the evidence-rate rule), **D-031** (a dict's key order can be part
   of the item pool's identity), **D-034** (a test may not carry its own answer key).

## Scope — two templates, and nothing else

| Template | File |
|---|---|
| `template_normal_depth_iteration` | `civil_engineering/water_resources/uniform_flow.py:287` |
| `template_line_balancing_heuristic` | `industrial_engineering/production_and_inventory/production_planning.py:158` |

**Do not touch `data/templates/` outside these two.** A defect found elsewhere is recorded in
`DECISIONS.md` or the phase summary — not fixed.

**Do not modify the harness or an oracle to make a fix pass.** If a check looks wrong, that is
a finding and a `SPEC-CHANGE`, raised explicitly.

**Do not push.** `master` is 50 commits ahead of `origin/master`; publishing is a separate
decision.

## What is actually wrong, measured

Verify these yourself before acting — every prior phase found the spec's defect list was a
starting point rather than an inventory.

**`template_normal_depth_iteration`** — a secant solve on Manning's equation. Measured over
400 seeds: the **step count is constant at 4**, but Step 3 carries **1 to 5 interpolation
updates** (distribution `{1: 5, 2: 31, 3: 230, 4: 128, 5: 6}`), each with its own `A`, `P`,
`A*R^(2/3)`, `g`. So there is no stable `{id, symbol}` for "the third trial depth" —
`y_curr`/`g_curr` are rebound every pass. Also:

- `g_curr - g_prev` is an **unguarded division**. Reachable in principle; 0 occurrences in
  6,000 seeds. Do not remove the risk by asserting it cannot happen — measure, then decide.
- `fmt3()` rewrites `"-0.000"` to `"0.000"`, so **the printed operand can differ in sign from
  the stored value.** That is a direct P2 violation, and T5 flags it three times as
  `[result-inline] fmt3(g_prev)` / `fmt3(g_curr)`.
- T1 currently has **1 hard failure in 200 seeds** (seed 11, the `K = Q*n/S^(1/2)` line).

**`template_line_balancing_heuristic`** — a greedy longest-eligible-task rule. The trace is a
**search log**: milestones are task *sets* (`Station 1 contents: {a, b}`), a running remaining
time, and boolean fit tests. Measured over 400 seeds: **the step count itself varies, 5 or 6**,
because the station count is 3 or 4 — and **`n` is a computed result, not an input.** Balance
delay is a function of `n`, so fixing the station count to make the trace a fixed shape
**gives away part of the answer**. T1 has 4 hard failures in 200 seeds; T5 flags 2
`[result-inline]` expressions.

## The decision you are here to make (D3.1)

For each template, independently: **schema route or redesign?**

The spec recommends the schema route for both — an `iteration` node type and a `decision` node
type, which it observes are the same shape (an ordered list of homogeneous sub-traces with
stable within-element symbols and a termination predicate). It offers a cheaper alternative for
§3.1: fix the iteration count at 3 with stable names, making it class A for ~4 h, at the cost
of pinning the pedagogy to "do exactly three secant updates."

**Phase 2 supplies a worked precedent for the cheap route, and its review supplies the limit.**
`template_adiabatic_flame_temperature` *was* an iterative solve, and Phase 2 resolved it
without any schema change — fixed 6 passes, every pass printed, answer to the nearest kelvin.
Phase 2's summary §10 offered you the criterion:

> *fix the iteration count when the count is not part of what the item tests.*

**Reviewer C showed that criterion is necessary but not sufficient, and you should treat §12 as
authoritative over §10.** The flame template is verifiable only because its parameter space is
**a closed set of eleven discrete reactions** — the rounding-boundary question could be checked
*exhaustively*. Normal depth samples channel geometry, slope and roughness **continuously**, so
a fixed iteration count cannot be exhaustively validated, and there will be a
measure-zero-but-nonempty set of parameters where the fixed count lands on the wrong side of a
display boundary. If you take the cheap route there, it needs **either a proven a-priori error
bound over the sampled box, or a run-time convergence assert that is not stripped by `python
-O`** — note that `assert` is, which is why the existing one guards development only.

And for `line_balancing_heuristic` the reasoning does not transfer *at all*: **a heuristic has
no fixed point**, so "N passes suffices" is not even the right shape of argument. The question
there is whether the heuristic is deterministic under tie-breaking — check what happens when
two eligible tasks have equal duration, and whether `_T19_ORDER` or a `set` iteration order
reaches the output.

**So: two templates, two separate decisions, argued separately.** Do not pick one route for
both because they look alike.

## Deliverables

D3.1 decision record per template, with the **P6 trade stated** · D3.2 template edits · **D3.3
the `iteration` and `decision` node-type specifications, with worked examples from both
templates — the primary deliverable** · D3.4 comparator semantics: how a model's iteration
count differing from gold is scored (**a correct answer reached in 4 iterations instead of 3
must not be marked wrong**) · D3.5 T6 distribution diff

Write the phase summary to `docs/re-implementation-sep/phase3_summary.md` and the item-pool
impact to `phase3_item_pool_impact.md`, following the Phase 1/C2/2 documents' shape.

## Verification

```
python -m tests.template_integrity.run --checks all          # corpus, ~60s
python -m tests.template_integrity.run --templates a,b --checks all
python -m tests.template_integrity.gate_report <N>
python -m tests.constants_integrity.test_chemical_thermochemistry
python -m tests.constants_integrity.test_citations_resolve
```

Corpus baseline at `master` — **verify it yourself against a `git worktree` rather than
trusting this table**, and use a short path such as `C:/wtm`, because the repo contains paths
long enough to fail a checkout under a temp directory:

| T1 | T2 | T3 | T4 | T5 | T6 | T7 | failing |
|---|---|---|---|---|---|---|---|
| 30 | 0 | 0 | 3 | 68 | 142 | 83 | 146 |

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and will
crash on `→`/`Σ`/`²`.

**T6 fails 142/150 on `master`.** The committed baseline is stale corpus-wide; this is a known
open item, raised as a `SPEC-CHANGE` in Phase 2 and deliberately not fixed there, because a
baseline refreshed by the phase it gates is not a gate. **D3.5 will hit the same wall.** Either
deliver the distribution diff as a direct before/after instance dump (Phase 2's approach — see
`phase2_item_pool_impact.md` §3), or regenerate the baseline **as its own deliberate change on
`master`** with the corpus-wide movement examined rather than absorbed. If you do regenerate
it, that is a decision worth its own record.

**State the smallest rate your run can resolve and check it against the rate you need to
exclude.** N seeds cannot resolve a defect rarer than ~3/N. This rule exists because acceptance
evidence was smaller than the defect rate three separate times (D-024, D-026).

**Instance dumps must run as separate processes per tree.** An in-process reload resolves both
sides to the same already-imported modules and reports everything identical. This is documented
in `tests/template_integrity/instance_dump.py` and has caught two people.

## Independent review — mandatory, two reviewers

Scope each to **R6**: one mandatory gate task, a stated time box, tooling supplied, settled
things fenced off, no measurement commissioned twice. **Give each a frozen ref** — a commit
SHA, not a branch name. In Phase C2 a branch moved under a reviewer and a quarter of its
findings were stale on arrival; Phase 2 fixed this and it worked.

*Reviewer B — Pedagogy.* Does either item still test what it tested? Specifically: **does the
`line_balancing` result still require the solver to *run* the heuristic, or has it become a
lookup?**

*Reviewer D — Schema.* **Given only the D3.3 node specs and a set of generated traces, write a
working verifier.** Report every ambiguity. **The spec passes only if a reviewer who did not
design it can implement against it.** This is the real test of D3.3 and the deliverable most
likely to fail its gate on first attempt — budget for a second round rather than treating one
failure as the end.

File both reports to `docs/re-implementation-sep/reviews/`, triage every finding **and** every
§5 suggestion (R4), and state honestly which you rejected and why. Note that
`docs/re-implementation-sep/reviews/` is now tracked — a bare `reviews` entry in `.gitignore`
had silently untracked all eight prior reports until Phase 2 caught it. Confirm your reports
are actually committed, not merely written.

## Exit gate

- [ ] D3.1 decisions recorded, **argued separately per template**, P6 trade stated
- [ ] T1–T7 pass on both templates
- [ ] `g_curr - g_prev` guarded; the `fmt3()` sign discrepancy removed
- [ ] **Reviewer D implemented a working verifier from D3.3 alone**
- [ ] D3.4 comparator semantics cover the iteration-count-mismatch case
- [ ] No corpus regression against the master baseline, measured not assumed
- [ ] Reviewers B and D filed, every finding and suggestion triaged (R2, R4)

Work on a branch. Merge to `master` with `--no-ff` when the gate is green.

## Open items inherited from Phase 2, for the summary rather than the code

These are recorded, not assigned. Fix them only if they block your gate; otherwise carry them
forward.

- **Two rounding conventions coexist in the corpus** — `_hu` (Decimal `ROUND_HALF_UP`, D-012)
  in the vibrations file, `_as_printed`/`format` (binary `ROUND_HALF_EVEN`) in the chemical
  files — and only some carry a display-tie guard. Reviewer A recommends a corpus-wide sweep of
  every template combining a transcendental with a fixed-decimal print. **`normal_depth_iteration`
  is exactly that shape** (`(A/P)**(2/3)` printed at 3 dp), so this one may land in your lap
  whether or not you go looking.
- **Reviewer C's two generalisable defects:** a raw arithmetic expression interpolated into an
  f-string (`{n-1}` printing `0.6000000000000001`), findable corpus-wide with a regex for
  `\{[a-z_]+ *[-+*/] *[0-9.]+\}`; and fixed *decimal* places applied to a quantity spanning
  decades, which cost 0.5% of an answer key in Phase 2.
- **Fallback firing-rate evidence lives in scratch scripts, not the repo**, so a reviewer cannot
  check it. If you remove a fallback, consider committing the measurement.
