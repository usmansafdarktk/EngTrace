# Evaluation-framework rework — the record

Everything here derives from the first section of
[`EngTrace_Suggested_Actions.pdf`](EngTrace_Suggested_Actions.pdf), which asked whether the 150
templates can support trace-level evaluation at all. The answer turned into a **ten-phase plan
across two parallel tracks**, and this directory is that plan's working record.

Read in this order if you are new to it:

1. **[`template_redesign_spec.md`](template_redesign_spec.md)** — the plan. Governing principles,
   every phase's deliverables and exit gate, the review protocol, sequencing and risks.
2. **[`audit/template_audit_report.md`](audit/template_audit_report.md)** — the 2026-09-05 audit
   that started it, and the baseline every later measurement is scored against.
3. The track you care about, below.
4. **[`DECISIONS.md`](DECISIONS.md)** — append-only, 77 entries, spanning both tracks. When a
   document and a decision disagree, the decision is why.

---

## Layout

```
├── EngTrace_Suggested_Actions.pdf   the source document
├── template_redesign_spec.md        the plan derived from it
├── DECISIONS.md                     append-only decision log, both tracks
├── audit/                           the 2026-09-05 baseline, and its re-measurement
├── track-a/                         template integrity — Phases 0-6
├── track-b/                         constants re-grounding — Phases C1-C3
└── reviews/                         every independent review, as filed
```

The rule is mechanical, so it stays true as files are added: **`phase0`–`phase6*` belong to
track-a, `phaseC*` to track-b, the audit artefacts to `audit/`, and reviews to `reviews/`.**

## The two tracks, and why they are separate

**Track A — template integrity (Phases 0–6).** Fixes what the templates *emit*: gold traces that
did not reproduce their own answers, unseeded generators, trace shape, the comparator contract,
output-contract hygiene, and finally consolidation and re-audit.

**Track B — constants re-grounding (Phases C1–C3).** Fixes what the templates *draw from*: the
provenance of every constants table, against sources acquired to `docs/references/`.

They run in parallel because only a small slice of B is on A's critical path, and they are bound by
**two hard sync points** — S1 before Phase 2, S2 before Phase 6 — so the item pool is regenerated
exactly once. The spec's "Sequencing and effort" section is the authority on this.

## Track A — [`track-a/`](track-a/)

| phase | what it did | summary |
|---|---|---|
| 0 | verification infrastructure and an independent baseline | [`phase0_summary.md`](track-a/phase0_summary.md) · [`phase0_baseline.md`](track-a/phase0_baseline.md) |
| 1 | round-trip integrity — traces that did not reproduce their own answers | [`phase1_summary.md`](track-a/phase1_summary.md) |
| 2 | determinism; solver-in-the-loop removed | [`phase2_summary.md`](track-a/phase2_summary.md) |
| 3 | trace shape for iteration and search; the node-type spec | [`phase3_summary.md`](track-a/phase3_summary.md) · [`phase3_node_types.md`](track-a/phase3_node_types.md) |
| 4 | the comparator contract — six `kind` values all 150 are scored by | [`phase4_summary.md`](track-a/phase4_summary.md) · [`phase4_comparators.md`](track-a/phase4_comparators.md) · [`phase4_vocabulary.md`](track-a/phase4_vocabulary.md) |
| 5 | output-contract hygiene, and binding the comparator to the corpus | [`phase5_summary.md`](track-a/phase5_summary.md) |
| 6 | consolidation and re-audit | [`phase6_summary.md`](track-a/phase6_summary.md) · [`phase6_residual_register.md`](track-a/phase6_residual_register.md) |

**Item-pool impact** is recorded per phase (`phase*_item_pool_impact.md`) and consolidated in
[`phase6_item_pool_impact.md`](track-a/phase6_item_pool_impact.md) — which is the one to read if
you need to know what the redesign invalidates and what correcting it costs.

## Track B — [`track-b/`](track-b/)

| phase | what it did | summary |
|---|---|---|
| C1 | classify every constants table; fix the provenance vocabulary | [`phaseC1_summary.md`](track-b/phaseC1_summary.md) · [`phaseC1_census.md`](track-b/phaseC1_census.md) |
| C2 | chemical thermochemistry — the critical path | [`phaseC2_summary.md`](track-b/phaseC2_summary.md) |
| C3 | the remaining tables, across all five branches | [`phaseC3_summary.md`](track-b/phaseC3_summary.md) · [`phaseC3_residual_register.md`](track-b/phaseC3_residual_register.md) |
| C3+ | the corrections tranche, owner-directed | [`phaseC3_corrections_summary.md`](track-b/phaseC3_corrections_summary.md) |

## audit/ — the baseline

[`template_audit_report.md`](audit/template_audit_report.md) and
[`template_inventory.csv`](audit/template_inventory.csv) are the 2026-09-05 artefacts. **The
inventory is read by code** — `generate_testset.py`, `derive_bindings.py`, `n1_candidates.py`,
`derive_vocabulary.py` and `regen_inventory.py` — so it is not free to move again without updating
them.

[`template_inventory_regen.csv`](audit/template_inventory_regen.csv) is Phase 6's re-measurement.
It sits **beside** the original rather than replacing it, so the two can be diffed: the original
remains the reference the reproduction is checked against.

## reviews/ — [`reviews/`](reviews/)

Every independent review, committed **unmodified and before any fix**. That ordering is the point:
a report edited after the fact is not an independent review. Reviewers are lettered by role —
A correctness, B pedagogy, C numerical, D schema, E comparator, F audit replication, G provenance,
H domain.

## Conventions worth knowing before you add to this

- **Every count ships with its predicate or a committed script.** A number without the rule that
  produced it cannot be checked, and several in this record turned out to mean different things
  under different rules.
- **Measure two trees in two separate processes.** An in-process reload resolves both sides to the
  same already-imported modules and reports everything identical. That trap has caught four people.
- **A summary that records no mistakes is not a clean phase; it is an unexamined one.** Each phase
  summary carries its own errors.
