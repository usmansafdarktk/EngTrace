# Session prompts — evaluation framework rework

Paste-ready prompts for running the evaluation-framework rework as a series of focused
sessions, one stage per session. Each prompt is self-contained: a fresh session needs no
prior conversation, only the repository.

Background for the whole effort is in
[`docs/re-implementation-sep/EngTrace_Suggested_Actions.pdf`](../re-implementation-sep/EngTrace_Suggested_Actions.pdf).

## The sequence

| # | Stage | Prompt | Status |
| --- | --- | --- | --- |
| 01 | **Template structure audit** — can the 150 templates emit structured traces at all, and at what cost? | [`01-template-structure-audit.md`](01-template-structure-audit.md) | ready |
| 02 | **Phase 1: round-trip integrity** — fix 12 templates whose gold traces do not reproduce their own answers | [`02-phase1-round-trip-integrity.md`](02-phase1-round-trip-integrity.md) | ready |
| 03 | Phase 2: determinism and solver-in-the-loop (4 templates) — gated on constants phase C2 | not written yet | |
| 04 | Phase 3/4: trace-shape schema and non-numeric comparators | not written yet | |
| 05 | Value-extractor spike — the gate for the whole project (D-002) | not written yet | |

The audit (01) reshaped what follows. It found the binding constraint is the *model* side,
not the template side, so the sequence now runs template-integrity phases first (they fix
gold traces that are already wrong in published results) with the extractor spike as the
gate for the milestone work proper. The full ten-phase plan across two tracks is in
[`template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md); these
prompts cover the phases as they become ready.

Prompts 03–05 are deliberately unwritten. Each depends on what the stage before it finds —
Phase 0 is the case in point: it discovered that the fix recipe Phase 1 was originally given
was insufficient, and would have propagated a defect eight more times. Write each prompt at
the end of the preceding session, while its findings are fresh.

## Conventions

- One stage per session. Stages are long enough that mixing them loses fidelity.
- Every stage writes its findings to `docs/re-implementation-sep/` so the next session can
  pick them up from disk rather than from conversation history.
- Audit and review stages are **read-only**. Implementation stages modify templates only on
  a branch, never on `master` directly.
- Every stage ends with a **mandatory independent review** (protocol R0–R6 in the spec) that
  files findings *and* forward-looking suggestions, each triaged before the stage closes.
- Scope every review to **R6** before dispatching: one mandatory gate task, a stated time
  box, tooling supplied, and no measurement commissioned twice. An over-scoped review stalls
  and leaves its gate unchecked — this has already happened once.
- Decisions and pivots go in
  [`DECISIONS.md`](../re-implementation-sep/DECISIONS.md), append-only.
