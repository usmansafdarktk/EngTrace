# R4 triage — the two reviews of the Phase 5 brief

**Date:** 2026-09-07 · **Artefact reviewed:**
[`docs/prompts/06-phase5-contract-hygiene-and-bindings.md`](../../prompts/06-phase5-contract-hygiene-and-bindings.md)
at frozen ref `ecff2d1`
**Reports:** [`phase5_prompt_facts.md`](phase5_prompt_facts.md) ·
[`phase5_prompt_scoping.md`](phase5_prompt_scoping.md)

R4 requires every §5 suggestion to be dispositioned before the artefact it governs is used,
and **an untriaged suggestion blocks the gate exactly as a CONFIRMED finding does**. That rule
has been enforced on every phase in this project; it applies to a review of a *brief* too, and
this table existed only after the omission was noticed on a readiness check.

Dispositions: `ADOPT-NOW` · `ADOPT-PHASE-<N>` · `SPEC-CHANGE` · `BACKLOG` · `REJECT`.
`BACKLOG` and `REJECT` carry a one-line reason.

---

## Facts reviewer — findings (§3)

| # | Finding | Disposition | Action |
|---|---|---|---|
| F1 | `19,668` is 149 templates, not 150; `continuous_to_discrete_conversion` raises and **satisfies the gate by crashing** | `ADOPT-NOW` | Named **N4** in the brief and in spec Phase 5 Track B; **D5.7b** added; both exit gates now count errors separately from verdicts |
| F2 | "invisible to T4" false for the five `**Final Answer**` templates — T4 *reports* and *passes* them | `ADOPT-NOW` | Corrected; the ungated/invisible distinction is now stated as a general rule |
| F3 | `22,982` archive pairs not reproducible (second derivation 23,292, splits diverging in opposite directions) | `ADOPT-NOW` | Number removed; replaced with "derive it and record the construction". Archive coverage (90/150 templates, 259 orphan traces) now stated |
| F4 | `58` and `87` have no operational definition | `ADOPT-NOW` | Both marked definition-dependent with ranges; brief now requires the predicate or a script |
| F5 | `Round 4 — Task 3` does not exist | `ADOPT-NOW` | Re-cited as the round-4 §3 `RECOMMENDATION`; **verified to resolve** (report line 1155) |
| F6 | `origin/master` **is** fetched; master 87 ahead, 0 behind | `ADOPT-NOW` | Corrected. My "~86" was nearly right and my *correction* to "not fetched" was the error |
| F7–F9 | "three people" → four; "eight prior reports" unsourced; spec says "the 8 touched" in an 11-template section | `ADOPT-NOW` (F7, F8) · `SPEC-CHANGE` (F9) | F7/F8 corrected in the brief. **F9 is against the spec, not the brief** — carried to Phase 5 to fix in `template_redesign_spec.md` where it lives |
| — | **Methodological note:** `compare_kind(kind, gold, candidate)`; wrong order reports **zero false accepts** across all 150 | `ADOPT-NOW` | Signature stated in §Verification **plus a required assertion**: the `hagen_poiseuille` false accept must reproduce before any sweep is trusted. Worth more than several findings |

## Facts reviewer — §5 suggestions

| # | Suggestion | Disposition | Action |
|---|---|---|---|
| 5.1 | Never state a count without its predicate; better, commit a script beside it | `ADOPT-NOW` | Stated as a standing rule in the brief: every count ships with its predicate or a regenerating script. **This is D-034 applied to prose, and the brief was the one artefact in the chain exempt from it** |
| 5.2 | A pair count must state its denominator in *templates* | `ADOPT-NOW` | D5.4/D5.5 now require `pairs / templates-scored / templates-skipped`, skips named |
| 5.3 | Distinguish "not gated by X" from "invisible to X" | `ADOPT-NOW` | Stated as a general rule beside the T4 passage, with the remedy each implies |
| 5.4 | Cite sections by literal heading text | `ADOPT-NOW` | Done for the one reference that failed; verified to resolve |
| 5.5 | Claims about repository history decay silently — source them or drop them | `ADOPT-NOW` | Three corrected. **The asymmetry is the lesson: everything from a run was right, everything from memory was not** |
| 5.6 | Say what the archive does *not* cover | `ADOPT-NOW` | 90/150 templates and 259 orphan traces now stated as D5.5's coverage ceiling |

## Scoping reviewer — §5 recommendations

| # | Recommendation | Disposition | Action / reason |
|---|---|---|---|
| 1 | **Split the phase: Track A becomes Phase 5, Track B becomes Phase 7** | **`REJECT`, with the reviewer's own fallback adopted in full** | The repo owner asked for these to be filed as *Phase 5* deliverables, so renumbering is not mine to take. **The technical core of the recommendation is adopted**: the tracks are explicitly sequenced (Track A merges first, Track B branches from that merge), *"they share nothing else"* is deleted as false, and the coupling is named with the five templates that carry it. This is the reviewer's stated minimum for keeping them in one phase. **Re-raise it if Track B overruns its box** — the split remains the cleaner design and this is a scoping call, not a disagreement about the facts |
| 2 | Give D5.6 an objective, a held-out corpus, a tie-break | `ADOPT-NOW` | Two axes, an adoption rule, and a 20% held-out slice chosen before any candidate is examined. Without the second axis the degenerate winner **decides nothing and passes** |
| 3 | Restore the false-reject half of Reviewer E's task | `ADOPT-NOW` | E's mandate is two-directional again; Track B's gate gains a decided-rate non-regression check. 13 of E's 20 Phase 4 findings were over-rejection |
| 4 | Amend Reviewer A's commission to the gate the brief repairs | `ADOPT-NOW` | A checks all four defect classes and owns the `ANSWER_MARKERS` seam; the Track A gate leads with the four-class item, not with T4 150/150 |
| 5 | State a threshold and a lever for D5.8 | `ADOPT-NOW` | N ≥ 50 (12 resolves nothing rarer than 25%), binding order restored, "named unbound with a reason is a result" |
| 6 | Make the stopping rule operational, or drop it | `ADOPT-NOW` | Replaced the maxim with the instrument that produced it: findings per reviewer-hour, and **stub-then-re-run before adding any rule**. Track B's "reduce scope" is *bind fewer templates and name the rest* |
| 7 | State the boxes as numbers | `ADOPT-NOW` | Both tracks' figures stated, with Phase 4's overrun as the calibration and an explicit stop trigger at 25 h |
| 8 | Carry the two dropped inheritances | `ADOPT-NOW` | R4-1 restored to D5.11; D-057 gets "decide and record" rather than a deferral |
| 9 | Cut, in this order | `ADOPT-NOW` (partial) | N1–N3 transcript reduced to one line each pointing at D-058; T6 collapsed from five mentions to two; three doubled warnings reduced to one each. **The brief grew from 408 to 485 lines despite the cuts**, because the corrections added more than the cuts removed. Recorded rather than hidden: if a later reader finds it long, the cut list in the report is where to start |

---

## What this triage changed about the readiness answer

The brief was **not** ready when it was first called done. Three things were outstanding and
only a readiness check surfaced them:

1. the SHA references had gone stale when `master` moved past `9e288d9`;
2. three §5 suggestions were partly actioned and none was dispositioned;
3. **this table did not exist**, which is the same omission — an untriaged suggestion — that
   the project's own review protocol treats as blocking.

The last one is the one worth carrying forward. **The rule was applied to every phase and not
to the reviews of the artefact that governs the next one.**
