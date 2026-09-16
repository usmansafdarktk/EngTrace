# Prompt 06 — Phase 5: output-contract hygiene, and comparator bindings

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Phases 0, 1, C2, 2, 3 and 4 are complete and merged. Phase 4's merge commit is `9e288d9`;
`master` has since taken **documentation-only** commits (this brief, D-058, the spec
amendments and two reviews of this brief). **Branch from `master`, not from `9e288d9`** — and
re-derive the table below at whatever `master` is when you start, rather than trusting that
those commits stayed documentation-only.

**This phase has two tracks and they are SEQUENCED, not parallel.** Track A edits eleven
templates. Track B binds the comparator to the corpus. **Track A lands and merges first**,
because it rewrites gold text that Track B's instruments consume — read
§"How the two tracks relate" before you plan anything.

---

I need you to implement **Phase 5**: **Track A**, output-contract hygiene across eleven
templates, and **Track B**, the comparator bindings Phase 4 specified but did not deliver,
plus the two instruments that make a binding checkable.

## What EngTrace is, and why this matters

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from
**150 parameterized Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each, in 47 files). Each template is
a function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of natural-language strings. The solution is the **gold
reasoning trace**.

**The paper has been rejected twice from ACL ARR** (January 2026 and May 2026). The most
damaging criticism, raised independently by five reviewers, is that the evaluation framework
verifies model reasoning using an "AI Tribunal" of three frontier LLMs while simultaneously
evaluating models from those same families. **Phase 4 built the comparator that removes the
LLM from the critical path — and bound 4 of 150 templates to it.** Track B is the rest.

Optional, but read it if a judgement call turns on what the benchmark is claiming:
`docs/_ARR_May__EngTrace.txt` (§4 the evaluation framework, Appendix M its validation) and
`docs/EngTrace_Rebuttal_Jul2026.txt` — **the rebuttal's promises are the real to-do list.**

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — an **LLM
  annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not flag it.
  (You *will* use `error_analysis_annotation/samples/`, which is different and is real
  archived model output.)
- `evaluation/` is a **separate track** with known defects (D-003). **Do not fix them.** You
  are extending its replacement.

## How the two tracks relate

| | Track A | Track B |
|---|---|---|
| edits | 11 templates in `data/templates/` | `tests/comparators/`, `template_inventory.csv` |
| gate | acceptance evidence for four defect classes | zero false accepts on cross-pairing |
| reviewer | **A** — correctness | **E** — comparator adversary |
| risk | changing an item pool by accident | binding a template wrongly and not noticing |

**They are not independent, and an earlier draft of this brief wrongly said they were.**
Track A's marker decision (D5.3) constrains `tests/comparators/normalize.py::ANSWER_MARKERS`,
which is Track B's edit surface. And Track A rewrites the gold text of eleven templates —
which is the corpus Track B's gold×gold cross-pairing consumes. **Five of Track A's eleven
carry gold×gold false accepts today** (`euclidean_distance_binary`, `limiting_reactant`,
`flow_system_molar_flow_rates`, `batch_moles_vs_conversion`, `decimation_aliasing_analysis`),
so they are in both tracks at once.

**Therefore: Track A lands and merges, then Track B builds its corpora.** If you run them
concurrently, every D5.4/D5.6 number is provisional until Track A merges and you must **re-run
rather than re-read** them.

**The seam has an owner, and it is Reviewer A.** `ANSWER_MARKERS`-agreement is the one
comparator file A may read; without that assignment the coupling is an exit-gate item that
only the implementer checks, which is the Phase 0 "green while measuring nothing" pattern
rebuilt inside this brief's own gate.

**Boxes.** Track A: 6–10 h implementation, 3–4 h review. Track B: 15–25 h, 8–10 h review —
and **Phase 4 budgeted 12–20 h + 8–10 h for four templates and spent eight review rounds**, so
treat Track B's figure as a trigger, not a forecast. **When Track B passes 25 h with D5.8
incomplete, ship the bound templates, name the unbound, and close.** If the phase cannot carry
both, Track A ships and Track B re-scopes.

## Read these first, in this order

1. **[`template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 principles **P1–P6**, all of **Phase 5** (both tracks), and the review protocol
   **R0–R6**. Read R0–R6 properly before you dispatch anything. Note **SPEC-CHANGE 11, 12
   and 13** — all three are Phase 4's.
2. **[`phase4_summary.md`](../re-implementation-sep/track-a/phase4_summary.md)** — **§8, §9 and §12
   are written for you.** §8 groups nineteen errors into three shapes and you are exposed to
   all three; §9 is what eight review rounds cost and bought; §12 states the stopping rule.
3. **[`phase4_comparators.md`](../re-implementation-sep/track-a/phase4_comparators.md)** — D4.1, the
   contract you are binding to. **§2 is the pair of rules every binding must satisfy**, §3 is
   `multipart`, §5 is what it does not settle.
4. **[`DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** — append-only. **D-058** creates
   Track B. Also **D-047** (`multipart` composes), **D-049** (three-valued verdict), **D-052**
   (units by declaration), **D-056** (hedges advisory), **D-034** (a test may not carry its own
   answer key), **D-024/D-026** (the evidence-rate rule).
5. **[`reviews/phase4_reviewer_e_comparator.md`](../re-implementation-sep/reviews/phase4_reviewer_e_comparator.md)** —
   four rounds. **Read the round-4 section's verdict and its §3 `RECOMMENDATION` before you
   write a line of comparator code**, and read **F0** before you quote any precision figure.

## What is actually here, measured

**Verify these yourself before acting — every prior phase found the brief's list was a
starting point rather than an inventory**, and *this brief was itself corrected twice by
review before you received it.* Phase 1 found its prescribed fix insufficient, C2 found one
recorded defect was four, Phase 2 found a template failing T1 on 91% of instances, Phase 3
found two display-tie populations where the brief listed one, Phase 4 found T6 failing on two
of four where the brief said one.

Corpus health, **re-measured and confirmed exact by an independent reviewer** at `9e288d9`
and unchanged at `c69a54d` — the first brief in this series whose corpus table needed no correction. Still
re-derive it in a worktree:

| Check | Failing | Note |
|---|---:|---|
| T1 printed-arithmetic closure | 29/150 | 97 more sit on the rounding boundary (marginals) |
| T2 round-trip oracle | 0/150 | only templates with an oracle are checked |
| T3 determinism | 0/150 | |
| T4 output contract | **3/150** | see below — it does not gate most of Track A |
| T5 binding / rounding | 66/150 | |
| T6 distribution | 142/150 | baseline stale corpus-wide (D-043); Phase 6 owns regeneration |
| T7 invariant asserts | 83/150 | advisory corpus-wide; a gate on templates a phase edits |

### Track A's gate does not gate Track A, and the fix is smaller than it looks

T4 **fails** on exactly `euclidean_distance_binary`, `cd_dc_system_analysis` and
`finite_convolution` — the spec's *malformed `**Step N:**`* row. Of the other eight:

- the five `**Final Answer**` templates **are already reported by T4** — it prints
  `non-canonical marker {'**Final Answer**': 25}` and **passes them anyway**. The remedy is a
  **severity change, not a new scanner.** An earlier draft of this brief said they were
  invisible to T4 and would have had you build a duplicate.
- the two malformed complex numbers and the unreachable `elif` **are** invisible to T4.

So "T4 passes 150/150" would be true with eight of the eleven edits unmade.

**Note the distinction, because it decides the remedy and this brief got it wrong once.**
*Ungated* means the check sees the defect and does not fail on it — a **severity** change.
*Invisible* means the check cannot see it — a **new instrument**. Five templates are ungated
and three are invisible; treating the first group as the second builds a duplicate scanner.
Apply the same care anywhere else you find a check that is green on a known defect.

**Fix the gate, then satisfy it** — extend T4's severity, add a check for the remaining three, or state
per-class acceptance evidence — and record it as a `SPEC-CHANGE` against Phase 5 Track A's
exit gate. Phase 0 found two checks that were green while measuring nothing; this is a third,
caught in advance.

### Track B's starting point (D-058 — re-derive all of it)

| | |
|---|---|
| templates bound to the comparator contract | **4 of 150** |
| gold×gold pairs, 12 instances each | **19,668** = **149** × 132 — *not* 150; see N4 |
| gold answers carrying a unit token | **68–87, definition-dependent** — see below |
| templates exposed to N1 | **56–64, definition-dependent** — see below |
| templates *measured* to produce an N1 false accept, 12 instances | **3** |
| archive×gold pairs available | **derive it yourself; two derivations disagreed** — see below |

**Three numbers in that table are soft and you must not treat them as hard.**

- **58/87 have no operational definition.** Two independent derivations gave 56–64 exposed and
  68–87 unit-carrying, depending on whether "carries a unit" means *at seed 0*, *on all 12
  instances*, or *on any*. **Define it, publish the definition beside the number, then use it.**
  D5.10 says "declare units for the 87" and that count is only meaningful once you have.
- **The archive×gold pair count is not reproducible.** One derivation gave 22,982 (11,384
  `scalar`, 5,511 `multipart`), another 23,292 (9,900 / 6,744) — diverging in *opposite*
  directions, so it is not an off-by-one and no construction is recorded anywhere. **Derive it,
  record the construction, and treat any prior figure as unverified.** Also unstated until now:
  **the archive covers only 90 of 150 templates and 259 traces map to no template at all.**
- **Never quote an exposure count as a defect count.** 56–64 templates are *exposed* to N1;
  **3** are measured to fail. An earlier draft of this brief reported the exposure count as the
  defect size. Necessary is not sufficient.

**And a rule for every count you produce, which this brief failed three times before you got
it.** D-034 says a number a script recomputes from the artefact is evidence and a number typed
beside it is not. **That applies to your prose as much as to your code.** Every count you state
— in a summary, a decision record or a review brief — ships with either its predicate written
out (*"87 = templates where all 12 gold spans contain a token from `<vocabulary>`"*) or a
committed script that regenerates it. The three numbers in this brief that failed review
(`58`, `87`, `22,982`) failed the same way: the number was recorded and the question it answers
was not. The measured tables, which came from runs, were flawless.

## Track A — output-contract hygiene

Eleven templates, four defect classes, in spec §Phase 5 Track A. `cd_dc_system_analysis` is
urgent: a strict marker parser drops **Step 3 on 25/25 seeds**, and Step 3 computes the answer.

**One decision to make and record (D5.3).** Five templates use `**Final Answer**` instead of
`**Answer:**`. Normalise the templates, **or** widen the accepted marker set corpus-wide —
decide once, apply everywhere, record it. `normalize.py::ANSWER_MARKERS` already accepts
several variants and **its priority order is load-bearing**; whichever way you decide, the two
must agree, and **Reviewer A owns checking that they do.**

`levenspiel_plot_interpretation` is also in Phase 2's scope — check whether that edit landed
before touching the file.

## Track B — bindings, and the instruments that make a binding checkable

Phase 4 specified the contract and bound four templates. Cross-pairing over the corpus found
four defects Phase 4's own corpora could not see, because `numeric` and `check` had **zero
archived negative instances** (Reviewer E, F0). **Full detail and reproductions are in D-058
and spec §Phase 5 Track B — read them there rather than here.** In one line each:

- **N1** — `parse_number` reads the **first** number in the answer span, often a fluid grade, a
  temperature or a chemical-formula subscript. Confirmed false accept on
  `hagen_poiseuille_flowrate`.
- **N2** — an **empty parse compares equal to an empty parse**: `{} == {}` is a `MATCH`,
  128/132 pairs on `autocorrelation_rect_pulse`.
- **N3** — **no bindings exist for 146 templates.** Larger than N1 and N2 together.
- **N4 — `continuous_to_discrete_conversion` RAISES** `AttributeError: 'tuple' object has no
  attribute 'free_symbols'` on every pair, so it contributes zero and **satisfies "zero false
  accepts" by crashing.** Found by the reviewer of this brief, not by the sweep that produced
  N1–N3. **A crash is not a pass; the gate below counts errors separately.**

**The reframing that matters: the comparator does not need more rules, it needs bindings, and
each binding needs evidence that it discriminates.** Phase 4 spent eight review rounds adding
rules to a natural-language surface and closed by *reducing* the mechanism's scope. If you find
yourself adding a rule to make one template pass, stop and ask whether a binding would do it.

### D5.6 — fix N1 by measurement, and the objective is stated

"First number" is wrong. **"Last number" is also wrong** — it breaks `check` answers like
*"18.4 mm, less than the 25 mm limit"*. Picking one from a single example is the error shape
`phase4_summary.md` §8 records **six times**.

Implement at least three candidates — the last number; the number adjacent to a declared unit;
`UNRESOLVED` when more than one number is present and nothing disambiguates — and add your own.
Then:

> Score every candidate on **two** axes and publish both for all of them:
> **(a) false accepts** on gold×gold and archive×gold;
> **(b) decided rate on answers that ought to be decided** — D4.4's `numeric` cases plus a
> positive set you build from real archive spans for templates with a declared unit.
>
> **Adopt the rule with the highest decided rate among those with zero false accepts.** If none
> reaches zero, say so and adopt on the stated trade rather than silently.
>
> **Hold out 20% of both corpora, chosen before you look at any candidate, and report the
> winner's numbers on the held-out slice.**

Without axis (b) the degenerate winner is candidate 3, which achieves zero false accepts by
deciding nothing; without the held-out split you have chosen and certified a rule on one
corpus, which is Shape 3 in `phase4_summary.md` §8. This brief's constraint that *every corpus
you add gets a negative control* has a mirror here: **a corpus of negatives certifying an
extraction rule needs positives.**

### D5.8 — bindings, with a threshold

> A binding counts as **bound** when it produces **zero false accepts over N ≥ 50 instances**
> cross-paired within its own template. Twelve instances resolves nothing rarer than 25%
> (D-024/D-026) and is not enough. Generation is deterministic and cheap, so raise N rather
> than lowering the bar.
>
> Bind in this order and **ship what is bound**: 5 `classification`, 9 `symbolic`, 15
> `vector`/`array`, 32 `multipart`. The 89 `scalar` need only a unit (D5.10).
>
> **A named unbound template with a written reason is a result. An unvalidated binding is a
> false accept waiting to happen.**

## Deliverables

**Track A:** D5.1 the eleven edits · D5.2 corpus-wide marker-conformance scan across all 150
proving zero remaining violations · D5.3 the marker-set decision record, agreed with
`ANSWER_MARKERS`.

**Track B:** D5.4 gold×gold cross-pairing as a standing check · D5.5 archive×gold
cross-pairing, Reviewer E's instrument promoted from review artefact ·

> **Both instruments report `pairs / templates-scored / templates-skipped`, with every skip
> named.** `19,668` reads as complete and is not — it is 149 of 150, and the missing one is
> N4. A pair count without its denominator in *templates* hides a crash as a pass, which is
> exactly how N4 survived the sweep that found N1–N3.

 D5.6 N1 fixed by
measurement, losing candidates' numbers recorded · D5.7 N2 fixed on every path · **D5.7b N4
fixed, and errors counted separately from verdicts in every corpus you build** · D5.8 bindings
declared and cross-pair-validated · D5.9 per-kind negative counts beside every precision
figure · D5.10 units declared, with the definition of "carries a unit" published · D5.11 the
Phase 4 residuals **R4-1** (the symbolic comparator is unvalidated — high severity, and the 9
`symbolic` templates are directly in your scope), R4-11, RB4-1, R4-12.

**D-057 needs a disposition, not a deferral.** `phase4_summary.md` §12 assigns Phase 5 the
supporting-quantity remedy for the two 100%-shortcuttable classification items; decide whether
you take it, and **record the disposition in the document that names it** either way.

**Both:** `phase5_summary.md` in the shape of `phase4_summary.md` — what changed, the R4 triage
tables, **both** exit gates, residual risk, what Phase 6 inherits, **and your own errors**. New
`DECISIONS.md` entries for every decision and reversal, including the ones you reverse
mid-phase. Review reports under `docs/re-implementation-sep/reviews/`.

## Verification

```
python -m tests.template_integrity.run --checks all          # corpus, ~25s
python -m tests.template_integrity.run --templates a,b --checks all
python -m tests.trace_schema.audit_3_8                       # must stay clean
python -m tests.comparators.score                            # archive + D4.4
python -m tests.comparators.reviewer_battery                 # 81 reviewer cases
python -m tests.comparators.recall_corpus                    # 507 cases
python -m tests.comparators.derive_vocabulary                # SUPPORT + commitment census
python -m tests.trace_schema.candidate_7 docs/re-implementation-sep/track-a/phase4_conformance/candidates.json
```

All eight exit 0 today, under both hedge policies. **`reviewer_battery` and `recall_corpus`
run under two policies** and the default is `advisory` (D-056); `ENGTRACE_HEDGE_POLICY=enforce`
checks the reviewers' cases against the policy they were argued for. **Both must stay green
under both**; green under one and red under the other is a finding, not a pass.

**The signature is `compare_kind(kind, gold, candidate)`.** Getting the order wrong raises
`KeyError` on all 150 templates and reports **zero false accepts** — a clean sweep that looks
exactly like N1–N4 being fixed. D5.4 and D5.5 are built on this call. **Assert a known
false accept before trusting any sweep**: `hagen_poiseuille_flowrate` seeds 3 and 9 must
report `MATCH` before your fix and `MISMATCH` after.

**Verify the corpus baseline against a `git worktree` at `master`** rather than trusting any
table. Use a short path such as `C:/wtm`; the repo contains paths long enough to fail a
checkout under a temp directory. **Remove the worktree before you merge** — it holds `master`
and the merge will refuse.

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and will
crash on `→`/`Σ`/`²`/`·`.

**Track A owes a before/after instance dump.**
`tests/template_integrity/phase4_instance_dump.py` is the pattern to copy: it takes the tree
root as an argument and refuses to run if a template resolves outside it. **Run each tree as
its own process.** An in-process reload resolves both sides to the same already-imported
`data.templates.*` modules and reports every instance identical — the trap documented in
`instance_dump.py`, which has now caught four people.

**T6 will move for Track A's marker edits and that is expected.** Its baseline is stale
corpus-wide (D-043), Phase 6 owns regeneration, and **regenerating it to make a gate green is
forbidden** — the before/after dump replaces it, as in Phases 3 and 4. Separately, T6's answer
extractor cannot see a non-scalar answer at all (R4-5, `phase4_summary.md` §6); its breach
lines on those templates compare parse artifacts. That is Phase 6's, not yours.

**State the smallest rate your run can resolve and check it against the rate you need to
exclude.** N samples cannot resolve a defect rarer than ~3/N (D-024, D-026). It binds hardest
on D5.8, which is why that threshold is N ≥ 50 rather than 12.

## Reviews — read R0–R6 before dispatching any

Two independent reviewers, one per track, **working in parallel and in isolation.**

- **Reviewer A — correctness (Track A).** *The gate.* Re-runs T4 across all 150, **and
  independently checks the acceptance evidence for the three defect classes T4 cannot see and
  the five it reports but passes.** You will hand A a check or a per-class evidence statement;
  **A treats it as a claim under test, not as tooling.** A also owns the
  `ANSWER_MARKERS`-agreement item — the one comparator file A may read. **Do not ask A about
  the rest of the comparator.**
- **Reviewer E — comparator adversary (Track B).** *The gate.* Given the bindings and nothing
  else: **find a binding that accepts an answer to a different instance of the same template,
  AND find a correct answer that the N1 rule leaves `UNRESOLVED`.** Both directions — 13 of
  E's 20 Phase 4 findings were over-rejection, and a comparator that resolves less than the one
  it replaced passes every accept-only gate. Draw from gold×gold and archive×gold, not from
  invented examples. **Do not ask E to re-run T4.**

**Scope every review to R6 before dispatching.** One mandatory gate task, **a time box stated
as a number**, tooling supplied as working code, everything settled fenced off explicitly, and
a measurement-ownership register with **no number commissioned twice**.

**Independence, concretely:**

- **R1.1 — work without sight of the implementer's reasoning.** Supply code comments, commit
  messages and summaries as **claims under test**, not as background.
- **Give each reviewer a frozen ref — a SHA, not a branch name.** Then **do not run
  `git commit` between dispatching a review and its filing.** Phase 3 broke this three times;
  Phase 4 kept it with exactly that mechanical rule — **and also had to stop *writing to the
  working tree*, not merely committing**, because a reviewer reading the tree at a frozen SHA
  sees uncommitted edits too.
- Tell reviewers to read with `git show <sha>:<path>` and **never `git show <sha> --stat`** —
  the latter prints the commit body and contaminated a Phase 3 reviewer.
- **Do not let the two reviewers see each other's briefs or reports.**

**Budget for more than one round, and know what "done" looks like.** Phase 3's specification
took six rounds; **Phase 4 took four per reviewer and closed not by exhausting findings but by
reducing the mechanism's scope.** Track B is a *binding* problem rather than a
natural-language one, so it should converge faster. If it does not, the stopping rule is a
**budget and an ablation, not a judgement**:

> **(a)** Record **findings per reviewer-hour** each round. When it stops falling, the surface
> is not converging and no further round will close it.
> **(b)** Before adding any rule to fix a finding, **stub the mechanism it lives in and re-run
> every corpus.** If deleting it costs zero decided verdicts on real archive text, delete it
> instead of patching it. **That ablation, not the maxim, is what ended Phase 4.**
> **(c)** Track B's version of "reduce scope" is **bind fewer templates and name the rest** —
> not "add fewer rules".

Reports go to `docs/re-implementation-sep/reviews/phase5_<reviewer>.md` in the R2 five-section
structure. **Every §5 suggestion must be triaged (R4) before the phase closes; an untriaged
suggestion blocks the gate exactly as a CONFIRMED finding does.** A `SPEC-CHANGE` **amends
`template_redesign_spec.md`**; an `ADOPT-PHASE-<N>` **edits that phase's deliverable list**.
Writing the disposition in your summary is not discharging it.

**State honestly which findings you rejected and why.** A `REJECT` with a written reason is a
legitimate and respected outcome; a finding quietly absent from the triage table is not.

**Confirm your reports are actually committed, not merely written** — a bare `reviews` entry in
`.gitignore` once silently untracked every prior report.

## Git

- Branch **`redesign/phase5-hygiene-bindings`** off `master`. Do not work directly on `master`.
- Logical commits, not one lump. **Commit each reviewer's report unmodified and before any
  fix**, so the record is what they wrote rather than what survived your response.
- **Do not push.** `origin/master` is at `1736447` and `master` is **87 ahead, 0 behind**;
  publishing is a separate decision and not this phase's.
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
- Merge Track A to `master` with `--no-ff` when its gate passes; **Track B branches from that
  merge**, not from the commit you started Track A at.

## Exit gate

**Track A** — the four-class item is first *because T4 alone is not the gate*.

- [ ] **Acceptance evidence exists for all four defect classes**, and the gate is amended so a
      later phase cannot repeat the mistake
- [ ] T4 passes **150/150**, with the five `**Final Answer**` templates no longer merely
      *reported*
- [ ] Marker-set decision recorded, and `normalize.py::ANSWER_MARKERS` agrees — **verified by
      Reviewer A, not by the implementer**
- [ ] T6 movement explained by the before/after dump; **baseline not regenerated**
- [ ] Item-pool impact stated, with the dump as evidence

**Track B**

- [ ] D5.4 and D5.5 land as runnable standing checks, with pair counts **and the construction
      that produced them** recorded
- [ ] **Zero false accepts on gold×gold across all 150**, with **errors counted and reported
      separately from verdicts** — a template that raises has not passed (N4)
- [ ] N1 fixed by a rule chosen on **measured evidence over two axes**, losing candidates'
      numbers recorded, winner reported on a **held-out slice**
- [ ] **The decided rate for `numeric` did not fall against the Phase 4 baseline, measured**
- [ ] N2 and N4 fixed: no code path returns `MATCH` from an unrecoverable parse, and none
      raises
- [ ] Per-kind negative-instance counts published beside every precision figure
- [ ] Every `classification`, `symbolic`, `vector`, `array` and `multipart` template either
      **bound and validated at N ≥ 50**, or **named unbound with a written reason**
- [ ] `reviewer_battery`, `recall_corpus`, `score`, `candidate_7` green **under both hedge
      policies**; `audit_3_8` clean; no corpus regression, measured not assumed

**Both**

- [ ] Reviewers A and E filed, **every finding and every §5 suggestion triaged (R2, R4)**, and
      every `SPEC-CHANGE` / `ADOPT-PHASE-N` actioned in the documents they name

## Constraints and cautions

- **Do not modify a check or an oracle to make something pass.** If a check looks wrong that is
  a finding and a `SPEC-CHANGE`, raised explicitly.
- **Do not build a binding to agree with the existing parser** (D-003, D-034). Reproducing its
  defects is the same error as a test carrying its own answer key — and N1 *is* that error,
  already committed once inside the comparator built to replace it.
- **A corpus that certifies a rule set must contain cases from outside it.** Reviewer B
  diagnosed this in Phase 4 round 1 — *"the corpus samples the inside of the list it
  certifies"* — and the implementer committed it again in round 2 by shipping a recall corpus
  with no negative control. **Every corpus you add gets a negative control before you quote a
  number from it.**
- **P6 is a real constraint.** Track A changes emitted text; Track B changes what counts as a
  correct answer for up to 146 templates, silently and corpus-wide. The second is the larger P6
  surface even though it edits no item.
- **Where you are uncertain, say so and explain what evidence would settle it.** Do not smooth
  over a gap to make a gate look clean.
- **Report your own errors in the summary.** Phase 4 recorded nineteen; thirteen shared three
  shapes. If your summary records none, that is not a clean phase, it is an unexamined one.

## Open items inherited — carry forward, fix only if they block your gate

- **D-003 is still open**: do the raw `inference_results/` generations still exist? It gates any
  corrected results table.
- **The corpus-wide sweeps from Phase 2** — templates combining a transcendental with a
  fixed-decimal print; the f-string regex `\{[a-z_]+ *[-+*/] *[0-9.]+\}`; two coexisting
  rounding conventions. None swept.
- **`_froude_capped_slope`'s exactness pathology affects three templates** in
  `civil_engineering/water_resources/uniform_flow.py`, not one (D-045).
- **`line_balancing_heuristic` is ~90% shortcuttable** (D-046) and **both Phase 4 classification
  items are 100% shortcuttable** (D-057). `blind_guess_floor` and `surface_model_heldout` are
  columns in `template_inventory.csv` for the five templates where the statistic means
  anything. **Do not build a binding on the assumption that an item tests what its difficulty
  label claims.**
- **The T6 baseline regeneration** (D-043) and **T6's blindness to non-scalar answers** (R4-5) —
  Phase 6.
- **`sympy` is a dependency and is not in `requirements.txt`** (D-053).
