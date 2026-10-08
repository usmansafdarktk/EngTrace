# Mock-review workstreams, October 2026 submission

One brief per isolated session. Each brief is self-contained: it carries the context it needs, the rules, the
files it may edit, its steps, its acceptance criteria and the report it must leave in `reports/`. The session
that wrote these files (the orchestrator) integrates the reports and runs Phase 2 in order.

Deadline: 12 October 2026 (ARR, anywhere on Earth). Today is 6 October. Phase 1 runs 7 to 9 October in parallel;
Phase 2 runs 10 to 12 October.

## Files

| File | Phase | Session | What it does |
|---|---|---|---|
| `PLAN_CONTEXT.md` | context | read by every session | Why these changes: the two mock reviews, what the record holds, decisions, the full work-package list and the traceability to every review item. Read sections 1 to 3 before any brief. |
| `00_ORCHESTRATION.md` | both | the orchestrator | Dependency graph, start order, handoff signals, decisions log, budget ledger, file-ownership matrix, label registry, Phase 2 sequence, master checklist, risks. |
| `WS-A_template_repair_and_rerun.md` | 1 | A | Repair the two chemical templates, re-certify (round 5), re-draw their 30 instances, re-run every model and every arm on them, re-score. |
| `WS-B_reasoning_on_inference.md` | 1 | B | Run every model that returned no reasoning tokens with reasoning on, on the full set; Gemma 4 calibration; the Qwen Thinking sibling; the paraphrase and repeat cascade. |
| `WS-C_analysis_code.md` | 1 | C (shared context) | Every new analysis as code: the single-path table, MC variants and their tests, scoring-clause counts, carried precision, the depth model, the matched-settings family, Table 1 redesign, per-area table, worked example. No tex edits. |
| `WS-C1_matched_family_and_exports.md`, `WS-C2_coverage_variants.md`, `WS-C3_scoring_variants_precision_depth.md`, `WS-C4_paper_machinery.md`, `WS-C5_statistics_tables_and_worked_example.md` | 1 | C1 to C5 | The WS-C split (2026-10-07): five parallel sessions of three to four hours, each owning its own scripts and result files; read the shared context first. |
| `WS-E_expert_tasks_and_symbolic_check.md` | 1 | E | Expert kits and scoring: difficulty ratings (E1), symbolic grading (E2), the error-analysis top-up (E3b); the symbolic equivalence check in `answer.py`; label adoption script. |
| `WS-F_figures.md` | 1 | F | Move the drawing code out of `paper_results.py`; redraw the four figures to the reviewers' points; approval gate; captions in the report. |
| `WS-G_integration.md` | 2, first | the orchestrator (or one session) | Final re-score with provenance, stages for new rows, `analyze.py`, `paper_results.py --write` on the real tree, every `--check`, the numbers sheet for the writers. |
| `WS-D1_main_text.md` | 2, after G | D1 | Title, abstract, introduction, related work, design, evaluation framework, results prose, conclusion, Limitations, Ethics. |
| `WS-D2_appendices_and_bibliography.md` | 2, after G | D2 | Every appendix's prose, the release appendix, the extended related work, the listings swap, the bibliography, the small corrections. |
| `reports/` | both | every session | `WS-<X>_report.md` per session, written from the template at the end of each brief. The orchestrator reads only these. |

## How to run a session

1. Open the brief. Read `PLAN_CONTEXT.md` sections 1 to 3 (ten minutes).
2. Check the brief's "Before you start" list. If an input is missing, ask the owner in that session; do not guess.
3. Work through the steps in order. Write new code freely within the files the brief assigns; do not touch files
   assigned to another stream (the matrix is in `00_ORCHESTRATION.md`, section 6).
4. When a step says "signal", write the signal line at the top of your report file at once, so the orchestrator
   and the dependent stream can proceed.
5. Finish with the report. Every number in it names the file or command that printed it.

## Start order (Phase 1, morning of 7 October)

1. **WS-F** first (30 minutes to extract the drawing code; signals FIGURES EXTRACTED).
2. **WS-B** at the same time (the first hour adds the harness variant and the item filter; signals HARNESS READY).
3. **WS-A**, **WS-C**, **WS-E** at once. WS-A's first hours need no harness; WS-C's first hours need no
   `paper_results.py`; WS-E's first hours are kits.
4. Phase 2 starts when the five reports are in: WS-G, then WS-D1 and WS-D2 in parallel, then the orchestrator's
   final checks and the second mock review.
