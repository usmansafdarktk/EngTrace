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
| 03 | **Phase 2: determinism and solver-in-the-loop** — 4 templates: seed every generator, delete silent fallbacks, replace `fsolve`/`quad` with stated algorithms | [`03-phase2-determinism.md`](03-phase2-determinism.md) | ready |
| 04 | **Phase 3: trace shape, iteration and search** — 2 templates whose traces are shaped by iteration; the one phase where the *schema* may be what changes, and where the primary deliverable is a specification | [`04-phase3-trace-shape.md`](04-phase3-trace-shape.md) | **done** — [`phase3_summary.md`](../re-implementation-sep/phase3_summary.md), spec at [`phase3_node_types.md`](../re-implementation-sep/phase3_node_types.md) v1.5 after five review rounds |
| 05 | **Phase 4: non-numeric templates and the comparator contract** — 4 templates, one cosmetic fix; specifies the six comparator `kind` values that all 150 templates are scored by, and delivers the conformance corpus Phase 3 left open | [`05-phase4-non-numeric-comparators.md`](05-phase4-non-numeric-comparators.md) | **done** — [`phase4_summary.md`](../re-implementation-sep/phase4_summary.md), spec at [`phase4_comparators.md`](../re-implementation-sep/phase4_comparators.md) and [`phase4_vocabulary.md`](../re-implementation-sep/phase4_vocabulary.md), after **four review rounds per reviewer** |
| 06 | **Phase 5: output-contract hygiene, and comparator bindings** — two SEQUENCED tracks. Track A fixes 11 templates' output contract; Track B binds the comparator to the corpus and lands the cross-pairing instruments that found four defects Phase 4 could not see (D-058) | [`06-phase5-contract-hygiene-and-bindings.md`](06-phase5-contract-hygiene-and-bindings.md) | **ready** — reviewed by two independent agents (facts, scoping); both reports in [`reviews/`](../re-implementation-sep/reviews/) |
| 07 | **Phases C1 + C3: constants sourcing and the remaining tables** — two SEQUENCED stages. C1 classifies every constants table and fixes the provenance vocabulary; C3 re-grounds the ~40 untagged tables in chemical, mechanical and electrical against sources acquired to `docs/references/`, and makes civil/industrial citations resolvable off one machine. Gates Phase 6 (sync point S2) | [`07-phaseC1-C3-constants-sourcing.md`](07-phaseC1-C3-constants-sourcing.md) | **ready** — sources pre-acquired by `docs/references/fetch_references.py`; not yet reviewed |
| 08 | Value-extractor spike — the gate for the whole project (D-002) | not written yet | |

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
