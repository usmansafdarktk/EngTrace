# Phase 4 — Summary and close-out

**Phase:** 4 — non-numeric templates and the comparator contract
**Branch:** `redesign/phase4-comparators`, off `master` at `51a4aef`
**Date:** 2026-09-07
**Deliverables:** [`phase4_comparators.md`](phase4_comparators.md) (D4.1) ·
[`phase4_vocabulary.md`](phase4_vocabulary.md) (D4.2) ·
[`tests/comparators/`](../../../tests/comparators) (D4.3, D4.4) ·
[`phase4_conformance/candidates.json`](phase4_conformance/candidates.json) (D4.6) ·
[`reviews/phase4_reviewer_e_comparator.md`](../reviews/phase4_reviewer_e_comparator.md) ·
[`reviews/phase4_reviewer_b_pedagogy.md`](../reviews/phase4_reviewer_b_pedagogy.md)

---

## 1. What this phase was for

Phases 1, C2, 2 and 3 made the gold traces trustworthy. **None of that scores a
model.** The comparator turns a gold trace into a verdict, and it did not exist —
so the AI Tribunal was still doing the verifying and the objection five reviewers
raised across two ARR rejections still stood.

One measurement frames everything below. On the 61 archived traces for the four
templates:

| | |
|---|---|
| scored correct by the deployed parser | **0 / 61** |
| actually correct, independently labelled | **45 / 61** |

The parser is not merely wrong on these items — it compares numbers that do not
exist. Gold `a) Memoryless: **No**` parses to `1.0`; `z[n] = {-9, 6, *-5*, …}` to
`-9.0`, the sequence's *first element*; `u = -2x^2 + 6xy` to `-2.0`, the *leading
coefficient*.

---

## 2. The three decisions

**Decision 1 — `multipart` is not a seventh `kind`** (D-047, SPEC-CHANGE 12). It
is a property of the *answer*: an ordered list of `parts`, each carrying one of
the six kinds, combined under `mode` `all` or `any`.
`system_properties_memory_causality` is the proof inside the phase — typed
`classification`, called a "label tuple" by the spec, actually two categorical
parts. **Governs 24 of the 51 class-C templates.**

**Decision 2 — route 1, and its argument failed its own test** (D-048). §5.

**Decision 3 — three-valued verdicts, biased toward refusing to decide**
(D-049). `UNRESOLVED` never contributes to a pass, and **costs recall exactly as
a wrong `MISMATCH` does**, so a comparator that declines everything scores 100%
precision and 0% recall. The gate carries the decision instead of hiding it.

---

## 3. What was delivered

| # | Deliverable | Status |
|---|---|---|
| D4.1 | Comparator specification, six kinds + `multipart` | `phase4_comparators.md` |
| D4.2 | Vocabulary, per-template counts, machine-checked | `phase4_vocabulary.md` |
| D4.3 | Reference implementation + tests per comparator | `tests/comparators/` |
| D4.4 | Adversarial set, ≥20 near misses per comparator | 159 cases |
| D4.5 | The cosmetic fix | applied; item-pool impact measured |
| **D4.6** | **Conformance corpus for D3.4 §7** | 16 candidates + 9 assertions, 15/15 dispositions, **4 findings against §7** |
| D4.7 | Third `iteration` template, binding only | attempted, **negative**, measured |
| D4.8 | Second `decision` template | surveyed, **none exists** |
| D4.9 | Trace-length semantics corpus-wide | surveyed, reported **with its own failure** |
| D4.10 | "Difficulty unchanged" operationalised | five proxies + control; superseded in part by D-057 |
| D4.11 | Reject lists generated, not remembered | three mechanisms |
| — | **Reviewer battery**, 81 cases | every case both reviewers used, all four rounds |
| — | **Recall corpus**, 507 cases | 351 positive + 156 negative — the artefact D4.4 was missing; **100% under both hedge policies** |

**D4.6 was the one that must not slip. It did not, and it produced the most.**

---

## 4. D4.6 — §7 exercised for the first time

Four defects **in the specification** (D-054, SPEC-CHANGE 13; actioned in
`phase3_node_types.md` §7).

**F7-1 — §7.2's headline disposition has no instances.** *"A correct answer
reached in four iterations instead of three is correct"* cannot occur for
`normal_depth_iteration`: §4.3.3, §8A.4, §4.3.1, §4.6 and the fixed preamble
together make the conforming node for a question **unique**. Verified both ways.
D-038's measured 1–5 spread is **between questions**, not between solvers of one.

> This is the six review rounds' bill arriving. Each closed a way for a candidate
> to lie; together they closed every way for a candidate to be *differently
> right*. Nobody recorded that trade, because §7 was never run.

**F7-2** — §7 never says whether §6 step 8 binds a candidate. It does, so answer
and process credit are **not independent** and "wrong answer, clean process" is
not constructible. **F7-3** — §7.4's direct solver cannot be an `iteration` node
at all. **F7-4** — §7.3 row 3's set-vs-ordered distinction is **vacuous**.

D4.9 reaches F7-1 from the other direction and on most of the corpus: **141 of
150 templates emit a trace whose length is a constant.**

---

## 5. D4.2 — an argument that did not survive its own test

61 traces across the four templates, **6 for the symbolic comparator**.
Resolution limits, per template and never averaged: 12%, 19%, 20%, **50%**.

Route 1 was taken and its argument *tested*, as the brief required. Over 27
rules, comparing per-model spread with per-template spread: **2 model habits, 15
item-driven, 10 both**. So coverage cannot be borrowed across templates — but the
*rules* can be, because a normalisation rule is a claim about meaning, not
frequency. The confound is named: the test cannot separate "the model does not
use this form" from "the item gave it no opportunity", so 15 is an upper bound.

**Residual, carried not solved:** for the symbolic comparator the honest position
is the brief's **route 3** — specified, adversarially exercised, **not validated
against a representative sample**. Five rules have zero Phase 4 support.

---

## 6. Verification

**Corpus baseline, re-measured against a `master` worktree, before and after:**

| Check | master | branch |
|---|---:|---:|
| T1 / T2 / T3 / T4 | 29 / 0 / 0 / 3 | 29 / 0 / 0 / 3 |
| T5 / T6 / T7 | 66 / 142 / 83 | 66 / 142 / 83 |
| T1 marginals | 97 | 97 |

**Zero cells moved.** `audit_3_8`: 80 nodes clean.

**Two divergences from the brief's own table, found by re-deriving it.**

1. The brief says `incompressible_continuity` "also fails T6". **T6 fails on two
   of the four** — `signal_operations` as well.
2. **T6 cannot see these templates' answers at all.** Its extractor returns the
   sequence's first element and the polynomial's leading coefficient; its breach
   line reads *"median answer moved 1.5 → -0.5"*, a comparison between two parse
   artifacts. This is D-003's defect class arriving through the **gate** rather
   than the evaluation track, and it will affect T6 on ~61 of 150 templates by
   `answer_type`. **Not fixed** — the brief forbids modifying a check to make
   something pass — and it is a stronger reason than staleness for T6 being
   uninformative here.

**Does T7 bind?** The brief asks explicitly. **No.** T7 is advisory corpus-wide
and a gate on templates a phase *edits*; the only edit is one column header's
presentation, which cannot add or remove an invariant assert. All four failed T7
on `master` with `0 asserts` and still do.

**D4.5 before/after**, 4,000 seeds, **each tree its own process**:

| | |
|---|---:|
| questions changed | **0 / 4000** |
| answers changed | **0 / 4000** |
| solutions changed | 4000 / 4000 — one table column header |
| — bracket only (shift) / + variable `y`→`z` (reversal) | 2024 / 1976 |

**Item pool unchanged.** Reviewer B read the diff and judged it pedagogically
neutral.

**A brief/spec divergence, resolved toward the spec and recorded.** The brief
calls D4.5 "the one-character-class fix"; spec §4.2.5 says "mismatched bracket,
**wrong variable**". The spec governs, so both halves are fixed and the dump
classifies them separately so either can be reverted alone.

---

## 7. Item-pool impact

No item pool changed: zero questions and zero answers differ on any of 4,000
seeds for any of the four templates.

**The larger P6 surface is not the edit.** This phase changed what counts as a
correct answer, corpus-wide. That is why Decisions 1 and 3, the hedge policy
(D-056) and the shortcuttability record (D-057) are all filed as pedagogy
decisions rather than implementation details.

---

## 8. Errors I made in this phase

Nineteen. **Three shapes account for thirteen of them**, and the shapes matter
more than the list.

### Shape 1 — I derived a rule from one example and it broke on the second (6)

| # | Error |
|---|---|
| 1 | The comparator double-negated fused negatives: `nonlinear` → `linear` on 3 of 15 archived linearity traces. |
| 2 | The sequence parser took the **last** brace group, which on a mixed-form answer is the *index list*. |
| 3 | `_isolate_expression` split on the **last** `=` (losing the answer to `= C(x)`), then on `.split("\n")[0]` (returning prose). Two wrong rules from two examples, in one function. |
| 4 | `_BRACED_RE`'s bracket alternative matched an index **subscript**, so `y[4]` was read as the value group. |
| 5 | The hedge list contained `it depends` and `either` — **domain vocabulary here** — producing `UNRESOLVED` on a correct archived answer. |
| 6 | Relaxing `_HYPO_RE`'s anchor for *all* hypotheticals matched ordinary past tense in "conditions **were** checked". The fix's own new surface, caught one run later. |

Every one was in a parser, every one in a template with fewer than 20 traces.
**Thin evidence does not only fail to validate a rule — it fails to generate
one**, and the second failure is quieter.

### Shape 2 — a claim about my own work verified by re-reading, not by running (4)

| # | Error |
|---|---|
| 7 | `normalize.py`'s support counts were written **from impression**: "44/61" and "47/61" against measured 57 and 55. D-034's failure mode, in the file arguing for D-034. |
| 8 | The probe that caught #7 was itself wrong first — it ran against the *answer span*, from which the marker has already been consumed, reporting 0/61. |
| 9 | I wrote that the presentation rules are model habits; **my own test contradicted it** (2 of 27). |
| 10 | Having rewritten #9 I overstated it the other way: "**every** Unicode rule and **most** LaTeX rules are item-driven". Measured: 2 of 3, and **2 of 7**. |

Phase 3 recorded ten errors and three shared exactly this shape. **It recurred
here at the same rate despite being written down.** All four were caught by
something that *ran* and disagreed. Writing the lesson down did not prevent it;
building a check that runs did.

### Shape 3 — I built a corpus that could not fail (3)

| # | Error |
|---|---|
| 11 | Two D4.4 labels were mine and wrong (`cat-16`, `tup-15`); a third (the E-F4 control) contradicted its own note. |
| 12 | **The recall corpus shipped with no negative control.** Reviewer B supplied one and 16 of 39 answers still matched under a hypothetical frame — 41% of the corpus was inert while reporting 100%. |
| 13 | The D4.10 shortcuttability proxy labelled instances by their **digit-masked** answer shape, collapsing every scalar template to one class and reporting **112 of 145** as fully shortcuttable. Caught by disbelieving the number. |

This is D-034 one level up, and **Reviewer B diagnosed it in round 1 and I
committed it again in round 2**: *"the corpus samples the inside of the list it
certifies."*

### The rest

| # | Error |
|---|---|
| 14 | The label-tuple comparator read only gold's encoding (`Memoryless: No`), scoring **21 of 24** archived traces wrong. |
| 15 | The §7 conformance candidates `C1`/`C8` were **structurally invalid nodes** — they tested the mutation, not the specification. |
| 16 | I deferred Reviewer E's F1 (the false trace-8 claim) after round 1 and **forgot it**; E had to file it again in round 2 and round 3. |
| 17 | My round-2 fix for the recall corpus (ignore non-enumerated prefixes) **relocated a defect** into the tuple path — a hypothetical frame was then dropped unread. Found by B's negative control. |
| 18 | `\b` written in a **bash heredoc** became a literal backspace (0x08) inside three compiled regexes. They compiled, matched nothing, and looked like a design error. Found by printing the pattern. It is the fourth time heredoc escaping cost me a fix in this session. |
| 19 | The trailing-question rule shipped guarded by `not _label_hit`, so `"Or is it?"` was caught and `"Or is it nonlinear?"` **credited** — the *more explicit* withdrawal passing. |

---

## 9. The reviews — four rounds each, and what they cost and bought

Two reviewers, **four rounds each**, **34 findings**, all dispositioned. Full
tables in the two reports; §11 and §12 carry the R4 triage.

**Four convergent findings**, reached independently by isolated reviewers — the
strong signal R1.7 describes:

| finding | E | B |
|---|---|---|
| label selection ranked by length across the whole span | R1-F2 | R1-F4 |
| concessive openers treated as hypotheticals | R2-F9a | R2-F1.1 |
| the missing marker→label governance relation | R2-F7 | R2-F1.2 |
| the undirected neighbour rule | R2-F9c | R2-F1.3 |

**The pattern, and it is Phase 3's exactly.** Every round found defects in the
layer added to fix the previous round's. Rounds 1→3 were relocations. Reviewer E
characterised it in round 4 and the characterisation is the phase's most useful
output:

> All five round-4 findings are **one- or two-token mutations of cases the suite
> already contains and passes**. The residual is *holding*; what collapsed is the
> cost of drawing from it — 6 findings per session in round 2, 10 in round 3,
> **5 in 25 minutes** in round 4 by mechanically editing the suite's own cases.
> 588 cases green, zero robustness to editing them.

**And E then made the measurement that ended the phase.** The commitment census
splits one number into two mechanisms with **opposite evidence** — hedge markers
2/2200 (0 in a gated kind), subordinator openers 43. Ablating the hedge layer:

| | full | ablated |
|---|---|---|
| recall-corpus positive frames | 351/351 | **351/351** |
| recall-corpus total | 507/507 | 477/507 |

**The whole layer buys 30 synthetic negative controls the reviewers wrote, and
zero archive verdicts and zero positive recall** — for ~250 lines and 13 of E's
20 findings, while repeatedly marking *correct* answers wrong.

So the phase closed by **reducing the mechanism's scope**, not by exhausting the
findings: hedges are now **advisory** — detected, annotated, counted, never
scored (D-056, SPEC-CHANGE 11). E declined "ship with the residual stated" and
said why: it would sign that for segmentation alone, but not while a layer with
no observed instances can silently mark a correct answer wrong.

**Reviewer B's round 4: the phase can close**, with one required edit (RB4-2 —
the shortcuttability caveat in D4.1 §8 plus the two inventory columns), both
actioned. B blocked three times on sound grounds and did not block a fourth,
which is what makes the close meaningful.

**What the reviews cost:** eight review rounds against a plan that budgeted one.
**What they bought:** every false accept in this document, the recall corpus, the
negative controls, the ablation, and the decision to stop.

---

## 10. Exit gate

- [x] All six `kind` comparators specified and implemented; `multipart` answered (D-047)
- [x] ≥95% precision and recall on D4.4 — **95.8% / 98.6%**
- [x] **Reviewer E found no false accept on real archived traces** — 0, across four rounds
- [x] Vocabulary traced to observed outputs, **coverage per template, not averaged** (D-048)
- [x] **D4.6 delivered**: §7 has a conformance corpus and every disposition is exercised (15/15)
- [x] D4.5 applied; no other change to `data/templates/`
- [x] `audit_3_8` clean; **no corpus regression, measured** (T1–T7 identical)
- [x] Reviewers E and B filed; every finding and §5 suggestion triaged; every `SPEC-CHANGE` actioned in the document it names
- [x] **Item pool unchanged**, with the cosmetic fix's before/after as evidence

**Two gate items are met and uninformative, and that is stated rather than
hidden.** E's F0: the archive contains **no wrong answer at all** for
`categorical` or `categorical[tuple]` and no trace of any kind for `numeric`,
`check` or `narrative`, so "no false accept on real archived traces" rests on 15
negative instances in two of six kinds. E's R2-F10: the commitment machinery has
**zero** archive instances in either direction.

**Three false accepts are declared, not concealed:** `num-16b` (right number,
wrong unit, none declared — D-052), and `cat-16` and `tup-18` (hedges naming the
right label, credited under the advisory policy — D-056). All three stay in the
corpus so the gate carries the trades rather than hiding them. Under
`ENGTRACE_HEDGE_POLICY=enforce` the last two are `UNRESOLVED` and precision is
97.2%; the difference between those two numbers **is** the cost of D-056, and it
is meant to be visible.

---

## 11. Residual risk register

| # | Risk | Severity | Owner |
|---|---|---|---|
| R4-1 | **The symbolic comparator is not validated against a representative sample.** 6 traces, ~50% resolution. Route 2 needs D-003. | high | Phase 5/6 |
| R4-10 | **The commitment machinery has no archive support** — 2 markers in 2,200 spans, 0 in a gated kind. A policy validated on constructed text. Mitigated by demoting it to advisory (D-056). | high | Phase 6 |
| R4-11 | **`require_origin` has no archive support for its discriminating power** — 0 of 16 traces reach that branch (E's F1). A route-3 rule. | medium | Phase 5 |
| R4-3 | Right number, wrong unit is accepted where no unit is declared (D-052). | medium | milestone model |
| R4-4 | `check` and `numeric` have **zero** Phase 4 archive support; `check`'s 21 verdict surfaces were written from the domain. | medium | first `check` template fitted |
| R4-5 | **T6 cannot see a non-scalar answer** and misreports on ~61 of 150. Not fixed; forbidden by the brief. | medium | Phase 6 (D-043) |
| R4-6 | `decision` remains exercised on one instance and one DAG; cannot be closed by this corpus (D-055). | medium | item design |
| RB4-1 | The bare-comment proxy treats **expletive** `it`/`there` and topic-shifting `we` as anaphors (B, round 4). Lands on `UNRESOLVED`, not a wrong score; 1-in-2,200 evidential standing. | medium | Phase 5 |
| R4-12 | **B's negative frames are right in kind, wrong in coverage**: each one's nearest linguistic neighbour fails, and *anchoring rather than vocabulary* decides. Two were fixed; the class is open. | medium | Phase 5 |
| R4-7 | 100%/97.8% on the archive is measured **on the corpus the vocabulary came from**. | stated | — |
| R4-8 | `sympy` added as a dependency, not yet in `requirements.txt` (D-053). | low | repo owner |
| R4-9 | D4.10's difficulty proxies are computed from the **solution string**, so a 100%-shortcuttable item and a hard one look identical (B, round 4). Superseded for classification items by D-057's two columns. | low | Phase 6 |

---

## 12. What Phase 5 and 6 inherit

**The stopping rule, which is the most transferable thing here.** A
natural-language surface cannot be finished by a rule set, and four rounds of
evidence say so quantitatively. The rule that follows: **when a mechanism's
findings become mutations of cases it already passes, and its census shows no
observed instances, reduce its scope rather than add a rule.** Budget later
natural-language work against an evidence threshold, not against a clean round.

**D-055's reframing.** `iteration` is not "a repeated sub-chain" — it is *a
convergence-terminated refinement of one quantity driven by its own residual*.
Both named generalisation targets are retired in `phase3_node_types.md` §10.

**D-057's two columns** (`blind_guess_floor`, `surface_model_heldout`) are in
`template_inventory.csv` for the 5 templates whose answer space is small enough
for the statistic to mean anything. **Three are over B's threshold.** Phase 5
owns the supporting-quantity remedy, where the item can change with its gold.

> **Phase 5's disposition, recorded here rather than only in its own summary
> (D-066): DECLINED, and reassigned to Phase 6.**
>
> The remedy changes what the item *asks* and what its gold *answers* — item
> design, not comparator work. Track B's charter is `tests/comparators/` and
> `template_inventory.csv` only, stated in the same spec paragraph that sets its
> effort box, and Track A's scope is the eleven output-contract templates.
> Taking it in Phase 5 would have been scope creep on the two items where a
> mistake is least recoverable, in a phase whose named risk is *changing an item
> pool by accident*.
>
> **What Phase 5 contributes instead.** Both items are now bound to the
> comparator, cross-paired at N = 50, credit a verbatim copy of their own gold,
> and produce **zero false accepts** — while remaining **100% shortcuttable**:
>
> | template | blind-guess floor | held-out surface model | gold×gold false accepts |
> |---|---:|---:|---:|
> | `system_property_linearity` | 0.5008 | **1.0000** | 0 |
> | `system_properties_memory_causality` | 0.3450 | **1.0000** | 0 |
>
> Those two facts are independent, and the reason to hold them side by side is
> that **a comparator which scores an answer correctly cannot tell you whether
> the answer required the reasoning**. A bound template is not a hard one, and
> the Phase 5 binding work must not be read as evidence about difficulty.
>
> One further measurement: `levenspiel_plot_interpretation`'s blind-guess floor
> is **1.0000** — a blind guess scores 100%, so the statistic is degenerate on
> that item and its difficulty label is unsupported by it.

**The rule this phase adds to D-034:** *a number written beside the rule it
describes is not evidence; a number a script recomputes from the artefact is.*
Four of my nineteen errors were numbers typed from impression, and every one was
caught by a check that ran.

**Still open and not this phase's:** D-003; the Phase 2 corpus-wide sweeps;
`_froude_capped_slope` across three templates (D-045); the T6 baseline
regeneration (D-043).
