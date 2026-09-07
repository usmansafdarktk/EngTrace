# Phase 5 — Summary and close-out

**Two tracks, sequenced.** Track A fixed eleven templates' output contract and
merged (`4dcc829`); Track B branched from that merge and bound the corpus to the
Phase 4 comparator. They share a phase number and nothing else: Track A is gated
by T4/T8 and Reviewer A, Track B by cross-pairing and Reviewer E.

**Both reviewers found the gate rather than the work.** Reviewer A: *"the
detector is narrower than the class it is named after."* Reviewer E: the
reported zero false accepts *"is produced by"* an over-rejection defect. Neither
found a corpus that was dirty; both found instruments that could not see.

---

## 1. What this phase was for

Phase 4 built the comparator that removes the LLM from the critical path — the
defect five ARR reviewers raised independently — and bound **4 of 150**
templates to it. Cross-pairing then found four defects Phase 4's own corpora
could not see, because `numeric` and `check` had *zero* archived negative
instances (Reviewer E, F0). Track B is the rest of that work.

Track A is adjacent and cheap: eleven templates whose emitted text breaks a
mechanical parser. `cd_dc_system_analysis` was urgent — a strict marker parser
dropped **Step 3 on 100% of instances**, and Step 3 computes the answer.

---

## 2. The decisions

| # | Decision |
|---|---|
| **D-059** | Gold emits **one** answer marker; the candidate-side list stays wide. They are different lists. |
| **D-060** | T4's non-canonical-marker finding becomes a **failure**, verified by a planted defect (SPEC-CHANGE 14). |
| **D-061** | The doubled-sign class is **16 templates, not 2** (SPEC-CHANGE 15). |
| **D-062** | Reviewer A's triage: every detector was narrower than its class (SPEC-CHANGE 16, 17). |
| **D-063** | N1 fixed by measurement over two axes; the rule is **per-kind**. |
| **D-064** | N2's root cause was the *isolation* — that function's third wrong positional rule. |
| **D-065** | A binding that never decides is not a binding. |
| **D-066** | D-057's supporting-quantity remedy declined, with a reason, and reassigned. |
| **D-067** | Reviewer E's triage: the two defect classes were **masking each other** (SPEC-CHANGE 18, 19). |

---

## 3. What was delivered

**Track A** — D5.1 the eleven edits · D5.2 `phase5_contract_scan.py`, four
detectors with a planted-defect self-test that reports its own false-positive
rate every run · D5.3 the marker decision with agreement stated as a checkable
predicate · **T8** (`checks/t8_emission.py`) in `DEFAULT_CHECKS`, so the two
classes T4 cannot see are gated by something other than memory ·
`phase5_instance_dump.py` and the item-pool note.

**Track B** — D5.4 gold×gold and D5.5 archive×gold as standing checks · D5.6 N1
fixed by measurement with twelve rules' numbers committed · D5.7 N2 fixed at its
root · D5.7b N4 fixed and errors counted separately from verdicts · D5.8 **118
of 150 bound and validated on four terms** · D5.9 per-kind negative counts ·
D5.10 the unit census with a published predicate · D5.11 the four Phase 4
residuals, each with a measurement (§8).

---

## 4. The numbers

**Track A**, 150 templates × 400 seeds = 60,000 instances:

| class | before | after |
|---|---:|---:|
| malformed `**Step N:**` | 3 templates | **0** |
| non-canonical answer marker | 5 | **0** |
| doubled sign, inside Track A | 2 | **0** |
| degenerate product in an answer span | 1 | **0** |
| generation errors | 0 | **0** |

Corpus: **T4 3 → 0 failing**, T8 0; T1/T2/T3/T5/T6/T7 identical, verified per
template against a `master` worktree under an order-insensitive comparison.

**Track B.** Three columns, because the middle one is the point:

| | Phase 4 | first Track B build | shipped |
|---|---:|---:|---:|
| templates bound | 4 / 150 | 132 / 150 | **118 / 150** |
| gold×gold false accepts | 560 | 0 | **0** |
| gold×gold **false rejects** | — | 56 | **0** |
| gold×gold **identity failures** | — | **75 of 132** | **0 of 1,416** |
| gold×gold UNRESOLVED | — | 1,514 | **0** |
| gold×gold errors | 132 | 0 | **0** |
| archive decided rate | — | 56.8% | **67.8%** |

**The middle column is what an accept-only gate reports.** It was green on the
one term anybody was looking at, and wrong on three that were not being counted
— one of which the tool was *already printing*.

D4.4 precision/recall are **unchanged from Phase 4 at 95.8% / 98.6%**, archive
100% / 97.8%, 0 archive false accepts, under both hedge policies.

**The 32 unbound, by measured reason:**

| reason | n |
|---|---:|
| the answer's shape varies across instances | 10 |
| false accepts at N = 50 | 10 |
| rejects a verbatim copy of its own gold | 9 |
| a part reads a constant, so it compares nothing | 3 |

---

## 5. Three numbers in the brief that did not survive derivation

The brief flagged three of its own figures as soft and told me to derive them.
All three moved.

**"68–87 templates carry a unit"** — under a published predicate the readings are
**seed 0 115 / any 116 / always 113 / invariant 103**. None is 87. And the
number turned out to matter less than expected, because **the census must not be
wired into the comparator at all** (§7, E-2/E-3).

**"22,982 / 23,292 archive×gold pairs"** — neither reproduces and neither had a
recorded construction. Mine is **18,115** (1,273 positive + 16,842 negative),
and the construction is in `cross_pair.py`'s docstring.

**"19,668 gold×gold pairs"** — arithmetically right and epistemically wrong: it
is 19,800 minus the 132 pairs of the template that *crashed*. Reporting a crash
as a smaller denominator is exactly how N4 survived the sweep that found N1–N3.

---

## 6. Item-pool impact

Full note in [`phase5_item_pool_impact.md`](phase5_item_pool_impact.md).
**4 of 11 templates changed emitted text; no question's content and no answer's
value changed on any of 22,000 instances.** The largest single change is a
correctness fix that was in no brief, spec or audit:
`sqrt(-41.2^2 + 86.45^2) = 95.77` evaluates **as printed** to 76.00 — a P1
violation on every instance with a negative component.

---

## 7. Errors I made in this phase

Fourteen. **Three shapes account for nine of them**, and two of the three are
shapes `phase4_summary.md` §8 already records.

### Shape A — a mechanism that looked green because nothing measured the thing it got wrong (5)

| # | Error |
|---|---|
| 1 | **Eleven bindings "passed" by deciding 0.0% of 2,450 pairs each.** The brief warns about this exact degeneracy in its D5.6 section; I committed it one deliverable later, on the axis where the warning was not repeated. |
| 2 | **75 of 132 bindings rejected a verbatim copy of their own gold**, and the cross-pairing loop skips `a == b`, so the case was excluded by construction. Reviewer E found it with one line. |
| 3 | **`cross_pair` printed `false rejects 56` on every run and I never put it in the claim table.** The number was on the screen. This is §8 Shape 2 — a claim about my own work verified by re-reading rather than by running — and the "running" had already happened. |
| 4 | **Four detectors were narrower than the class they were named after** (Reviewer A, F1/F3/F5/F6), and my planted-defect self-test passed anyway, because I wrote the plants and the regexes with the same hand. |
| 5 | **I claimed in a decision record that the census "prints the list on every run"** when it printed a count. Caught by reading my own claim against the output. |

### Shape B — a diagnosis from an error message rather than from a measurement (3)

| # | Error |
|---|---|
| 6 | **"T6 is nondeterministic."** It was the *ordering* of a list. Corrected in the item-pool note before review; Reviewer A then found my replacement number ("7 of 150") was not reproducible either. |
| 7 | **A "zero-collapsed decimal" corpus finding that was a regex backtracking artefact** — `0.000083` matching `0\.0{2,}(?![1-9])`. I nearly recorded it as a defect on 8 templates. |
| 8 | **"92.8% of archived gold is stale."** 83.6% of that is simply a different seed. Real drift is 9.2%, and 126 of those 179 rows were my *own* Track A marker change. |

### Shape C — inferring what should be declared (2)

| # | Error |
|---|---|
| 9 | **`DECLARED_UNITS` was a trailing-token heuristic**, producing `otherwise`, `e-05`, `percent`, `dollars` as "units" — and D4.1 §4.1 says in terms why the unit check is opt-in and declared. I re-created the exact rejection that paragraph exists to prevent, and it cost **19 of 82** matches on real archived answers. |
| 10 | **One declared unit applied to every part of a multipart answer**, whose parts carry different units by construction. 13 templates rejected their own gold. |

### The rest

| # | Error |
|---|---|
| 11 | **Shell and CRLF quoting corrupted source four times** — `\n` collapsing inside a *quoted* heredoc, two CRLF anchor failures, and backticks inside a double-quoted `bash -c` eating three comment fragments. `phase4_summary.md` §8 #18 records this as "the fourth time heredoc escaping cost me a fix"; for me it was four times in one phase. Fixed durably only after the fourth, by writing every patch to a file. |
| 12 | **My first class-4 plant landed in a derivation step**, where the detector is designed not to fire, and I read the correct refusal as a miss. |
| 13 | **I ran the corpus scan at 25 seeds and it reported this phase's own headline defect as absent** — class 4 fires at 1.95% and 25 seeds resolves ~12% (D-024/D-026). On the first run. |
| 14 | **The identity check I added for E-1 called `compare_template` before the binding was installed**, so all 73 candidates failed with the same `KeyError` and the gate rejected them for the wrong reason. Caught because 73 identical error messages is not what a real finding looks like. |

**What actually caught them.** Not one of these was found by reasoning about the
code. #1–#5 were found by a reviewer or by an instrument reporting a number I
had not asked for; #6–#8 by re-running with a different question; #9–#10 by a
reviewer; #11–#14 by an assertion that ran. **The two errors I caught myself
(#5, #14) were both caught by disbelieving a result, not by inspecting a
diff.** Phase 4 §8's conclusion — *"writing the lesson down did not prevent it;
building a check that runs did"* — reproduces here exactly, including for the
lessons this project has already written down twice.

---

## 8. D5.11 — the Phase 4 residuals, each with a measurement

**R4-1 — the symbolic comparator is unvalidated.** Now measured, and the answer
is worse than Phase 4 knew: **8 of the 9 `symbolic` templates decide 0.0%** of
their 2,450 cross-paired instances. They were passing "zero false accepts" by
never accepting anything. All eight are now **named unbound** (D-065). The
obvious fix — masking function names in `_to_sympy`, since `cos` was becoming
`c*o*s` — was implemented and changed the decided rate by **nothing at all**,
11.1% before and after; the real blockers are a unit word inside the expression,
an isolation that takes prose, and two-equation answers. **Assigned to Phase 6
(D6.10) with the three shapes named and one known-bad fix ruled out.**

**R4-11 — `require_origin` has no archive support.** Phase 4: 0 of 16 traces
reached the branch. Answered by ablation:

| corpus | pairs / traces | verdicts changed by the rule |
|---|---:|---:|
| gold×gold, both sequence templates | 1,104 | **0** |
| archive, both templates | 25 | **6** (5 `MATCH`→`UNRESOLVED`, 1 →`MISMATCH`) |

No effect on gold — gold states the origin consistently, so both sides agree —
and **6 of 25** real archived verdicts changed, every one in the direction of
declining to credit. **R4-11 is closed**: the branch has observed instances and
they are the discriminating ones.

**RB4-1 and R4-12 — the bare-comment proxy and the negative frames.** Answered
by the ablation the brief prescribes — stub the mechanism, re-run every corpus:

| corpus | advisory (shipped) | enforce |
|---|---|---|
| `score` (D4.4 + archive) | **identical** | **identical** |
| `recall_corpus` | **identical** | differs |
| `reviewer_battery` | differs | differs |

Ablating `is_bare_comment` changes nothing on D4.4, nothing on the archive and
nothing on the recall corpus under the shipped policy. Its entire observable
footprint is the reviewer battery — cases the reviewers themselves wrote. Both
are **`BACKLOG`**, and the justification is a measurement rather than a shrug.

---

## 9. D-057 — disposition, recorded in the document that names it

**Declined by this phase, with a reason, and reassigned to Phase 6.** See D-066
and the note added to `phase4_summary.md` §12. The remedy changes what the item
asks and what its gold answers; Track B's charter is `tests/comparators/` and
`template_inventory.csv` only.

**What this phase adds instead:** both items are now bound, cross-paired at
N = 50, credit a verbatim copy of gold, and produce **zero false accepts** —
while remaining **100% shortcuttable** (held-out surface model 1.0000 against
floors of 0.5008 and 0.3450). Those facts are independent, and holding them
side by side is the point: **a comparator that scores an answer correctly cannot
tell you whether the answer required the reasoning.** A bound template is not a
hard one.

---

## 10. Verification

All green, under both hedge policies (`advisory` default and `enforce`):

```
python -m tests.template_integrity.run --checks all        # T1-T8, 150 templates
python -m tests.template_integrity.phase5_contract_scan    # D5.2, 60,000 instances
python -m tests.template_integrity.phase5_contract_scan --selftest   # planted defects
python -m tests.template_integrity.phase5_contract_scan --markers    # D5.3 agreement
python -m tests.comparators.derive_bindings                # D5.8/D5.10, checks for drift
python -m tests.comparators.cross_pair                     # D5.4 + D5.5
python -m tests.comparators.n1_candidates --check          # D5.6, twelve rules
python -m tests.comparators.score                          # archive + D4.4
python -m tests.comparators.reviewer_battery               # 81 reviewer cases
python -m tests.comparators.recall_corpus                  # 507 cases
python -m tests.comparators.derive_vocabulary              # SUPPORT census
python -m tests.trace_schema.audit_3_8
python -m tests.trace_schema.candidate_7 docs/re-implementation-sep/phase4_conformance/candidates.json
```

Corpus, unchanged by Track B and moved only by Track A's T4:
T1 29 · T2 0 · T3 0 · **T4 0** · T5 66 · T6 142 · T7 83 · **T8 0**.

---

## 11. Exit gates

### Track A

- [x] Acceptance evidence for all four defect classes, by **planted defect**;
      the gate amended so a later phase cannot repeat it (SPEC-CHANGE 14, 16, 17)
- [x] T4 passes 150/150, the five `**Final Answer**` templates no longer merely reported
- [x] D5.2 clean within Track A's eleven, corpus-wide residual as a **named census** of 14
- [x] Marker decision recorded and the predicate checked — **verified by Reviewer A**
- [x] T6 movement explained by the before/after dump; baseline **not** regenerated
- [x] Item-pool impact stated, separating a marker rename from an answer-body change
- [x] The two classes T4 cannot see are in a **standing** gate (T8, `DEFAULT_CHECKS`)
- [x] Reviewer A filed; every finding and §5 suggestion triaged (D-062)

### Track B

- [x] D5.4 and D5.5 are runnable standing checks reporting **pairs /
      templates-scored / templates-skipped**, every skip named, with the construction
- [x] **All four terms zero on gold×gold**: 0 false accepts, 0 false rejects,
      0 identity failures of 1,416, 0 errors
- [x] N1 fixed on measured evidence over two axes, twelve rules' numbers
      committed, winner reported on a **held-out slice frozen first**
- [x] The decided rate did not fall — and the incumbent's higher archive rate was
      itself the defect, which is measured rather than asserted
- [x] N2 and N4 fixed: no path returns `MATCH` from an unrecoverable parse, none raises,
      errors counted separately
- [x] Per-kind negative counts beside every rate
- [x] **118 bound and validated at N = 50 on four terms; 32 named unbound with a
      measured reason**
- [x] No unit inferred from gold's string (SPEC-CHANGE 19)
- [x] All suites green under both hedge policies; no corpus regression, measured
- [x] Reviewer E filed; every finding and §5 suggestion triaged (D-067)

**Two gate items are met and carry a stated trade.** `symbolic` is bound on
**1 of 9** templates, and `vector`'s archive decided rate is **12.0%**. Both are
reported rather than smoothed: the alternative was eight bindings that decide
nothing, which is what the first build shipped.

---

## 12. Residual risk register

| # | Risk | Severity | Owner |
|---|---|---|---|
| R5-1 | **`symbolic` is bound on 1 of 9 templates.** The comparator cannot represent a unit word inside an expression, a piecewise answer, or a two-equation answer. | high | Phase 6 (D6.10) |
| R5-2 | **113 templates carry a unit and none is declared to the comparator.** Right-number-wrong-unit is accepted corpus-wide — D-052's residue, now measured and larger than Phase 4 knew. | high | Phase 6 (D6.11) |
| R5-3 | **The gold×gold truth predicate is textual identity of the answer span.** Two instances carrying the same answer written differently would read as a false accept. Zero observed, but the proxy is a proxy. | medium | Phase 6 |
| R5-4 | **N = 50 resolves ~0.12% per template.** E found a 1-in-31,000 collision by sweeping 250 instances; the shipped gate cannot see one. | medium | Phase 6 |
| R5-5 | **9.2% of archived gold no longer matches what the template emits** (9 templates), so D5.5 pairs some model answers against a retired gold. 126 of the 179 rows are Track A's own marker change. | medium | Phase 6 |
| R5-6 | **14 templates still emit a doubled sign**, 2 in the answer span, 5 in a question. Censused and named. | medium | Phase 6 (D6.7) |
| R5-7 | **`multipart` archive decided rate is 52.7%.** Models frequently state fewer parts than gold. Not a comparator defect; an open question about what partial answers should score. | medium | milestone model |
| R5-8 | The two Phase 4 classification items remain **100% shortcuttable** (D-057/D-066). | stated | Phase 6 |
| R5-9 | `sympy` is still not in `requirements.txt` (D-053). | low | repo owner |

---

## 13. What Phase 6 inherits

**The rule this phase adds, and it is the one that would have saved it:**
*a gate with one term is passed by a mechanism that fails on the others, and
failure modes mask each other.* "Zero false accepts" was green while 75 of 132
bindings rejected their own gold — and the zero was **produced by** those
rejections. Both reviewers arrived at the same shape from opposite directions:
A found detectors narrower than their class, E found a gate narrower than its
property.

**The corollary for review scoping:** the plant and the detector must not be
written by the same hand. Four of Reviewer A's findings survived a
planted-defect suite that passed, for exactly that reason (SPEC-CHANGE 17), and
E's decisive finding is a case the loop *excluded by construction* — the kind of
blind spot only someone who did not write the loop will look for.

**And the one that keeps recurring:** three of this phase's fourteen errors are
shapes `phase4_summary.md` §8 already names, including one it explicitly says
was "the fourth time". Writing a lesson down does not prevent it. **Only a check
that runs does** — and this phase's evidence for that is that every error it
caught itself was caught by disbelieving a number, never by reading a diff.

New deliverables filed against Phase 6: **D6.7–D6.11**.

**Still open and not this phase's:** D-003; the Phase 2 corpus-wide sweeps;
`_froude_capped_slope` across three templates (D-045); the T6 baseline
regeneration (D-043); T6's blindness to non-scalar answers (R4-5).
