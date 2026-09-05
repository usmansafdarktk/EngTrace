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
| 02 | Literature review — prior art on process supervision, deterministic verification, judge bias | not written yet | |
| 03 | Independent design review — parallel agents propose and critique instrumentation designs | not written yet | |
| 04 | Implementation — instrument the pilot slice against the chosen design | not written yet | |
| 05 | Tests and experiments — correctness of the instrumentation, then the evaluator bake-off | not written yet | |

Prompts 02–05 are deliberately unwritten. Each depends on what the stage before it finds:
the literature review should be aimed by the obstacles the audit surfaces, the design review
should respond to real constraints rather than imagined ones, and there is no point
specifying implementation or tests before a design exists. Write each one at the end of the
preceding session, while its findings are fresh.

## Conventions

- One stage per session. Stages are long enough that mixing them loses fidelity.
- Every stage writes its findings to `docs/re-implementation-sep/` so the next session can
  pick them up from disk rather than from conversation history.
- Audit and review stages are **read-only**. Only stage 04 modifies templates, and only on
  a branch.
