# Prompt 05 — Phase 4: non-numeric templates and the comparator contract

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Phases 0, 1, C2, 2 and 3 are complete and merged; Phase 3's last content commit is
`1a71b70` and this prompt sits on top of it. This stage depends on Phase 3
in particular: it inherits six named deliverables from Phase 3's review triage, and it
inherits one gap Phase 3 could not close.

---

I need you to implement **Phase 4** of the template redesign: the four templates with **zero
numeric content**, and the **comparator contract** that all 150 templates will eventually be
scored by. **This is a specification-and-implementation phase with almost no template
editing** — one cosmetic fix is the only code change to `data/templates/`.

## What EngTrace is, and why this matters

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from
**150 parameterized Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each, in 47 files). Each template is
a function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of natural-language strings. The solution is the **gold
reasoning trace**.

**The paper has been rejected twice from ACL ARR** (January 2026 and May 2026). The most
damaging criticism, raised independently by five reviewers, is that the evaluation framework
verifies model reasoning using an "AI Tribunal" of three frontier LLMs (GPT-5, Claude Opus 4.5,
Gemini 3) while simultaneously evaluating models from those same families — GPT-5 is both judge
and evaluated model. The paper is titled *"Verifiable Process Supervision"*, but verification is
a majority vote among LLMs.

**Phase 4 is the phase that removes the LLM from the critical path.** Phases 1, C2, 2 and 3
made the gold traces trustworthy — they reproduce their own answers, rest on constants verified
against NIST, regenerate from their recorded seeds, and (for the two iterative items) emit a
structured trace. **None of that scores a model.** The comparator is what turns a gold trace
into a verdict, and it does not exist yet. Until it does, the AI Tribunal is still the thing
doing the verifying, and the reviewers' central objection stands.

That is the actual stake in this phase. The four templates are the vehicle; the comparator
contract is the deliverable.

For the full picture (optional, but read it if a judgement call turns on what the benchmark is
claiming):

- `docs/_ARR_May__EngTrace.txt` — the current paper. **Section 4** is the evaluation framework,
  **Appendix G** has three full template implementations, **Appendix M** is the framework
  validation.
- `docs/EngTrace_Rebuttal_Jul2026.txt` — the most recent reviewer objections and our responses.
  **These promises are the real to-do list behind this work.**

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not
  flag it. (You *will* be working in `error_analysis_annotation/samples/`, which is different
  and is real archived model output.)
- `evaluation/` is a **separate track**. Its parser has known defects (it reads the gold answer
  `**4,921**` as `4.0`), recorded as D-003. **Do not fix them** — but do read the parser before
  you design a comparator, because you are specifying its replacement and its defects are the
  clearest available statement of what not to do.

## The position, which is unusual and load-bearing

**These four templates are not defective and must not be redesigned.** They have zero numeric
content and they are all deterministically verifiable by a non-arithmetic comparator. Adding a
numeric payload to make them fit a numeric checker would damage four sound conceptual items to
satisfy an implementation convenience. **The schema adapts, not the items.**

Phase 3 is the precedent and it went the other way — there, the schema changed to fit two
templates whose trace shape was data-dependent. The principle is the same in both directions:
the item states what it tests, and the machinery accommodates it.

| Template | File | Comparator kind |
|---|---|---|
| `template_signal_operations` | `electrical_engineering/signals_and_systems/discrete_time_signals.py:6` | sequence-with-origin |
| `template_system_properties_memory_causality` | `.../discrete_time_signals.py:152` | categorical (label tuple) |
| `template_system_property_linearity` | `.../discrete_time_signals.py:426` | categorical (single label) |
| `template_incompressible_continuity` | `mechanical_engineering/fluid_mechanics/fluid_kinematics.py:546` | symbolic (CAS) |

## Read these first, in this order

1. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 governing principles **P1–P6**, all of **Phase 4** (§4.1–4.5, including the six
   deliverables Phase 3's triage added), and the review protocol **R0–R6**. Read R0–R6
   properly before you dispatch anything; it is short and it is the part most often skipped.
   Note **SPEC-CHANGE 8, 9 and 10** in the revision log — all three are Phase 3's, and 9 and 10
   change what your own reviews must do.
2. **[`docs/re-implementation-sep/phase3_node_types.md`](../re-implementation-sep/phase3_node_types.md)** —
   the trace node specification, schema 1.5. **§7 is the part you are extending**, and §3.8 is
   the rule you will be judged against. Read §10's non-goals: two of them are now your
   deliverables.
3. **[`docs/re-implementation-sep/phase3_summary.md`](../re-implementation-sep/phase3_summary.md)** —
   **§11.6, §12.6 and §13 are written for you.** §13.3 states the gap you are inheriting;
   §12.6 is a worked example of a reviewer overturning its own verdict; §11.6 carries the one
   sentence that predicts how this phase will fail if it fails.
4. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   append-only; add entries, never edit existing ones. **D-038** (cardinality: incidental vs
   answer-bearing), **D-039** (the node is a frame local — you may change this), **D-046**
   (an item 90% shortcuttable, and how it was found), **D-034** (a test may not carry its own
   answer key), **D-024/D-026** (the evidence-rate rule).
5. **[`docs/re-implementation-sep/reviews/phase3_reviewer_d2_schema.md`](../re-implementation-sep/reviews/phase3_reviewer_d2_schema.md)** —
   six rounds against one specification. Read §9's closing paragraph before you write a line of
   comparator code.

## What is actually here, measured

**Verify these yourself before acting — every prior phase found the spec's list was a starting
point rather than an inventory.** Phase 1 found its prescribed fix insufficient, Phase C2 found
one recorded defect was four, Phase 2 found a template failing T1 on 91% of instances the spec
did not mention, and Phase 3 found two display-tie populations where the brief listed one.

Measurements I have already taken, which you should re-derive rather than trust:

- **The four templates pass T1–T5 and fail T7** (advisory corpus-wide, but a gate on templates a
  phase edits). `incompressible_continuity` also fails T6. Since this phase makes one cosmetic
  edit, decide explicitly whether T7 binds here and record the answer.
- **The archive holds 2,200 traces** across 11 models in
  `error_analysis_annotation/samples/*.jsonl`, 200 per model, with `problem_statement`,
  `gold_answer`, `model_reasoning`, `final_answer_acc` and a `question_id` whose stem is the
  template name. This is the evidence base for D4.2 and the only source of real model
  phrasing you have.
- **But only 61 of those 2,200 are for your four templates**: `memory_causality` 24,
  `signal_operations` 16, `linearity` 15, **`incompressible_continuity` 6.**

**That last number is the single most important fact in this brief, and it is a problem you
must solve before D4.2 is meaningful.** The spec requires the normalisation vocabulary to be
"derived from archived model outputs, not invented". Six traces cannot establish a vocabulary
for the symbolic comparator: by the standing rule (D-024, D-026), N samples cannot resolve a
variant rarer than ~3/N, so six traces are blind to anything occurring in under half of model
outputs. You have three honest routes and must pick one explicitly:

1. Draw the vocabulary from the **whole 2,200**, on the argument that hedging and label
   phrasing are model habits rather than template-specific — then *test* that argument rather
   than assert it.
2. Generate new traces for the thin templates, which needs D-003 resolved (do the raw
   generations still exist?) and is probably out of scope.
3. Specify the symbolic comparator **without** an observational vocabulary and state plainly
   that it is unvalidated against real output, carrying it as a residual risk.

Do not paper over this by counting the 2,200 and implying the coverage is uniform.

Corpus health at `master` (`1a71b70`), for context — **this is the baseline you work against,
not a to-do list.** Re-measure it yourself in a worktree rather than trusting the table:

| Check | Failing | Note |
|---|---:|---|
| T1 printed-arithmetic closure | 29/150 | 97 more sit on the rounding boundary (marginals) |
| T2 round-trip oracle | 0/150 | only templates with an oracle are checked |
| T3 determinism | 0/150 | |
| T4 output contract | 3/150 | malformed step markers |
| T5 binding / rounding | 66/150 | |
| T6 distribution | 142/150 | **the committed baseline is stale corpus-wide — see Verification** |
| T7 invariant asserts | 83/150 | advisory corpus-wide; a gate on templates a phase edits — **including all four of yours** |

Phase 3 moved exactly three of these cells and nothing else. **Your phase should move none of
them except through D4.5**, and if it moves one you have edited something you should not have.

## The leverage, and the reason this phase is not small

The four comparators do not serve four templates. `template_inventory.csv` records an
`answer_type` for all 150:

| `answer_type` | count | comparator kind |
|---|---:|---|
| scalar | 89 | numeric |
| multipart | 32 | composite — **not in the spec's list of six; decide what it is** |
| symbolic | 9 | symbolic |
| vector | 8 | sequence |
| array | 7 | sequence |
| classification | 5 | categorical |

The spec says these comparators "also serve the 51 class-C templates" and that **Phase 4 is
where most of class C's design cost is actually paid**. Verify that 51 — the inventory CSV has
no `class` column, so the number comes from the audit report and I could not confirm it. But
the shape of the claim is right and the table above is the reason: **`multipart` is 32
templates, nearly a quarter of the corpus, and the spec's six `kind` values do not obviously
contain it.** Decide whether `multipart` is a seventh kind, a composition of kinds, or a
property of a milestone rather than a kind — and record the decision.

## Scope — four templates, one cosmetic edit, and a contract that outlives them

**Do not edit `data/templates/` except for D4.5**, the one-character-class fix at
`discrete_time_signals.py:113`. A defect found elsewhere is recorded in `DECISIONS.md` or the
phase summary — not fixed. That includes defects in the four templates themselves: they are
sound, and the position of this phase is that the schema adapts to them.

**Do not modify a check or an oracle to make something pass.** If a check looks wrong, that is
a finding and a `SPEC-CHANGE`, raised explicitly. Phase 0 found two checks that were green
while measuring nothing.

**Do not "fix" `evaluation/`.** You are specifying its replacement, not repairing it.

What you *are* free to change, and should consider changing: **`phase3_node_types.md` §7**,
which is the comparator half of the trace-node spec and is the thing D4.6 will exercise for the
first time. Phase 3 shipped it unexercised and said so. If exercising it shows it is wrong, say
so and amend it — that is the deliverable working, not a regression. **D-039 is likewise open
to you**: the structured trace is a frame local rather than a return value because changing the
return contract was out of Phase 3's two-template scope. It may be in yours.

## The three decisions you are here to make

Make each explicitly, record each in `DECISIONS.md`, and state the trade rather than absorbing
it. **Do not let any of the three be settled implicitly by an implementation choice.**

**Decision 1 — is `multipart` a seventh `kind`?** It is 32 of 150 templates, nearly a quarter
of the corpus, and the spec's six `kind` values do not obviously contain it. Three readings:
a seventh kind; an ordered composition of other kinds; or a property of a *milestone* rather
than a *kind*, in which case the six stand and the answer schema grows instead. The choice
determines how a partially-correct multipart answer scores, which is a P6 question with real
consequences for 32 items.

**Decision 2 — where does the D4.2 vocabulary come from, given the archive is thin?** 61 traces
across the four templates, six of them for the symbolic comparator. The three routes are in
"What is actually here". Pick one, state the resolution limit it leaves you with, and do not
average the coverage across templates to make it look uniform.

**Decision 3 — which way is the comparator biased, and why?** A false accept credits reasoning
that did not happen; a false reject penalises a correct solver. The first is the AI Tribunal
criticism in a new costume, and this whole effort exists to answer that criticism — which is an
argument, not a proof, for biasing toward rejection. But a comparator that rejects correct
answers makes the benchmark measure phrasing rather than reasoning. Decide, quantify the
asymmetry you are accepting on D4.4, and make the exit gate's "≥95% precision and recall" carry
your decision rather than hide it.

## Deliverables

From the spec (§4.3), unchanged:

1. **D4.1** — comparator specification for all six `kind` values
   (`numeric`, `categorical`, `sequence`, `symbolic`, `narrative`, `check`), plus whatever you
   conclude about `multipart`.
2. **D4.2** — categorical label-normalisation vocabulary, **derived from archived output with
   frequency counts**, not invented. Read the §"What is actually here" warning first.
3. **D4.3** — reference implementation plus unit tests per comparator.
4. **D4.4** — adversarial test set: **≥20 hand-written near-miss model answers per
   comparator**, both correct-but-phrased-differently and wrong-but-similar.
5. **D4.5** — the one cosmetic fix: `discrete_time_signals.py:113` prints a table header
   `New Value y[k')` — it opens a bracket and closes a parenthesis, and the sibling column
   reads `New Index (k')`. This is the only edit to `data/templates/` this phase.

Added by Phase 3's R4 triage, each a reviewer suggestion dispositioned `ADOPT-PHASE-4`:

6. **D4.6** — **a conformance corpus for D3.4 §7.** *The largest gap Phase 3 left.* Every
   verifier built in Phase 3 implements §6, gold-checking. The comparator rules that decide a
   model's score have **no corpus at all**. Needs candidate traces: wrong count with right
   answer, right count with different packing, early stop, fabricated frames. **This is the
   deliverable that connects Phase 3's work to a score, and both Phase 3 schema reviewers named
   it independently as the top gap.**
7. **D4.7** — fit a third template to the `iteration` node type **by declaring a binding only**,
   no verifier edit permitted (`linear_reservoir_routing_step` or `qr_policy_one_iteration`).
   Turns "argued, not measured" into measured. A *synthetic* renamed node already passes; a real
   one has never been tried.
8. **D4.8** — fit a second `decision`-shaped template. Until then `decision` is a type by
   construction and by the rename test, on one instance and one precedence DAG.
9. **D4.9** — classify trace-length semantics corpus-wide, once, applying D-038's
   `incidental` / `answer_bearing` test as a survey.
10. **D4.10** — operationalise "difficulty unchanged" so a pedagogy reviewer has something to
    measure rather than reason about.
11. **D4.11** — generate reject lists from the verifier rather than writing them alongside it.
12. **`phase4_summary.md`** — close-out in the shape of `phase3_summary.md`: what changed, the
    R4 triage tables, exit-gate status, residual risk register, what Phase 5 inherits, **and
    your own errors**.
13. **New entries in `DECISIONS.md`** for every decision or reversal, including the ones you
    reverse mid-phase.
14. **Two review reports** under `docs/re-implementation-sep/reviews/`.

**That is a lot, and the six inherited items are not all equally urgent.** If the phase cannot
carry all of them, **D4.6 is the one that must not slip** — say so explicitly and re-scope the
rest rather than delivering all of them thinly.

## Verification

```
python -m tests.template_integrity.run --checks all          # corpus, ~35s
python -m tests.template_integrity.run --templates a,b --checks all
python -m tests.trace_schema.audit_3_8                       # 3.8 audit, must stay clean
python -m tests.trace_schema.extract <out.json> <n>          # regenerate the trace corpus
python -m tests.constants_integrity.test_chemical_thermochemistry
python -m tests.constants_integrity.test_citations_resolve
```

**Verify the corpus baseline yourself against a `git worktree` at `master`** rather than
trusting any table. Use a short worktree path such as `C:/wtm`; the repo contains paths long
enough to fail a checkout under a temp directory. **Remove the worktree before you merge** —
it holds the `master` branch and the merge will refuse.

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and will
crash on `→`/`Σ`/`²`.

**The one before/after you do owe** is the cosmetic fix's. The exit gate asks you to state
that no item pool changed, and D4.5 changes emitted text, so you need a dump of the affected
template from a `master` worktree and from your branch. **Run each tree as its own process.**
An in-process reload resolves both sides to the same already-imported `data.templates.*`
modules and reports every instance identical — documented in
`tests/template_integrity/instance_dump.py`, and it has now caught three people.
`tests/template_integrity/phase3_instance_dump.py` is the pattern to copy: it takes the tree
root as an argument and refuses to run if a template resolves outside it, so the mistake fails
loudly instead of producing a clean-looking null result.

**T6 and your four templates.** T6 fails 142/150 corpus-wide because the committed baseline is
stale, and `incompressible_continuity` is one of the 142. **Do not regenerate the baseline** to
make it green — a baseline refreshed by the phase it gates is not a gate (D-043), and Phase 6
owns the regeneration. Since you edit only one template's *formatting*, the honest statement is
that T6 is uninformative here and the before/after dump replaces it, exactly as Phase 3 did.

**State the smallest rate your run can resolve and check it against the rate you need to
exclude.** N samples cannot resolve a defect rarer than ~3/N. This rule exists because
acceptance evidence was smaller than the defect rate three separate times (D-024, D-026), and
in this phase it bites hardest on D4.2, where one template has six traces.

**T6 fails 142/150 on `master`** — the committed baseline is stale corpus-wide. Do not
regenerate it to make your gate green; a baseline refreshed by the phase it gates is not a gate
(D-043). Phase 6 owns the regeneration.

## Reviews — read R0–R6 before dispatching any

Phase 4 requires **two independent reviewers, working in parallel and in isolation from each
other.**

- **Reviewer E — comparator adversary.** *The gate.* Given the comparators and nothing else:
  **find a correct model answer they reject, and an incorrect one they accept.** Draw from the
  real archived traces, not invented examples. Report precision and recall on D4.4 plus its own
  additions. **A comparator that has not been attacked has not been verified** — Phase 3's gold
  corpus passed 80/80 against a verifier that hardcoded one template's symbol names, because
  gold traces do not lie. Yours will not either.
- **Reviewer B — pedagogy (the P6 guard).** Does the normalisation vocabulary accept answers a
  human grader would reject? A hedge that never commits is the specific case named in the spec.
  Also: **SPEC-CHANGE 10 now binds you** — a pedagogy lookup-check must **fit a model, not
  enumerate rules**, and report lift over a blind-guess floor against a stated threshold. In
  Phase 3, enumerating rules gave 65% and fitting a model gave 90%, which overturned the
  reviewer's own verdict (D-046).

**Scope every review to R6 before dispatching it.** One mandatory gate task, a stated time box,
tooling supplied as working code, everything already settled fenced off explicitly, and a
measurement-ownership register with **no number commissioned twice**. A Phase 0 review stalled
and produced nothing because its brief duplicated another agent's — the only mandatory task
went unanswered because it was bundled with work already being done elsewhere. Do not ask a
reviewer to re-measure something a prior phase established; say so explicitly in the brief, and
run R6's pre-dispatch checklist every time.

**Independence, concretely:**

- **R1.1 — work without sight of the implementer's reasoning.** Supply code comments, commit
  messages and summaries as **claims under test**, not as background.
- **Give each reviewer a frozen ref — a commit SHA, not a branch name.** And then **do not
  commit to that branch until the review files.** I broke this three times in Phase 3, twice
  after recording it as an error; each time the reviewer caught it. A frozen ref is a promise
  about the repository, not a string in a brief. The rule that works is mechanical: do not run
  `git commit` at all between dispatching a review and its filing.
- Tell reviewers to read files with `git show <sha>:<path>` and **never `git show <sha> --stat`**
  — the latter prints the commit body and contaminated a Phase 3 reviewer that was meant to be
  blind to it.
- **Do not let the two reviewers see each other's briefs or reports.**

**Budget for more than one round.** Phase 3's brief anticipated a second round on its
specification; it took **six**, and five of the six found a defect *in the layer added to fix
the previous round's defect*. The fix was sound every time; the surface it introduced had not
been reviewed. **Whatever you add to fix a finding, review the mechanism it is built on, not
the fix.** That is the single most useful sentence Phase 3 produced and it was written by a
reviewer, in round six.

Reports go to `docs/re-implementation-sep/reviews/phase4_<reviewer>.md` in the R2 five-section
structure. **Every §5 suggestion must be triaged (R4) before the phase closes; an untriaged
suggestion blocks the gate exactly as a CONFIRMED finding does.** And note what R4 actually
requires, which Phase 3 recorded and initially failed to do: a `SPEC-CHANGE` **amends
`template_redesign_spec.md`**, and an `ADOPT-PHASE-<N>` **edits that phase's deliverable list**.
Writing the disposition in your summary is not discharging it.

**State honestly which findings you rejected and why.** A `REJECT` with a written reason is a
legitimate outcome and a respected one; a finding quietly absent from the triage table is not.

**Confirm your reports are actually committed, not merely written** — a bare `reviews` entry in
`.gitignore` once silently untracked all eight prior reports.

## Git

- Branch **`redesign/phase4-comparators`** off `master`. Do not work directly on `master`.
- Logical commits, not one lump.
- **Do not push.** `master` is ~70 commits ahead of `origin/master`; publishing is a separate
  decision.
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
- Merge to `master` with `--no-ff` only after the exit gate passes, or record explicitly what
  is left open.

## Exit gate

- [ ] All six `kind` comparators specified and implemented, and the `multipart` question answered
- [ ] ≥95% precision and recall on D4.4
- [ ] **Reviewer E found no false accept on real archived traces**
- [ ] Normalisation vocabulary traced to observed outputs, not invented — **with its coverage
      per template stated, not averaged**
- [ ] **D4.6 delivered: D3.4 §7 has a conformance corpus and its dispositions are exercised**
- [ ] D4.5 cosmetic fix applied; no other change to `data/templates/`
- [ ] `audit_3_8.py` still clean; no corpus regression, measured not assumed
- [ ] Reviewers E and B filed, **every finding and every §5 suggestion triaged (R2, R4)**, and
      every `SPEC-CHANGE` / `ADOPT-PHASE-N` actioned in the documents they name
- [ ] Item-pool impact: **state explicitly that no item pool changed**, with the cosmetic fix's
      before/after as the evidence

## Constraints and cautions

- **P6 is a real constraint, and this phase's version of it is unusual.** You are not changing
  items — you are changing what counts as a correct answer to them. A comparator that accepts a
  hedge weakens every item it scores, silently and corpus-wide. That is a larger P6 surface
  than any template edit.
- **A comparator's failure mode is asymmetric and you should say which way you biased it.**
  A false accept credits reasoning that did not happen — which is precisely the AI Tribunal
  criticism this whole effort exists to answer. A false reject penalises a correct solver. They
  are not equally bad here; decide, state, and justify.
- **Do not build the comparator to agree with the existing parser.** It has known defects
  (D-003) and reproducing them would be the same error as a test carrying its own answer key
  (D-034).
- **Where you are uncertain, say so and explain what evidence would settle it.** Do not smooth
  over a gap to make the gate look clean.
- **Report your own errors in the summary.** Phase 3 recorded ten; three shared one shape — a
  claim about its own work verified by re-reading rather than by running something. If your
  summary records none, that is not a clean phase, it is an unexamined one.

## Open items inherited — carry forward, fix only if they block your gate

- **D-003 is still open**: do the raw `inference_results/` generations still exist? It gates any
  corrected results table, and it gates route 2 for the thin-archive problem above.
- **The corpus-wide sweeps from Phase 2** — templates combining a transcendental with a
  fixed-decimal print; the f-string regex `\{[a-z_]+ *[-+*/] *[0-9.]+\}`; two coexisting
  rounding conventions. None swept.
- **`_froude_capped_slope`'s exactness pathology affects three templates** in
  `civil_engineering/water_resources/uniform_flow.py`, not one (D-045). Phase 5 territory.
- **`line_balancing_heuristic` is ~90% shortcuttable** from its question text (D-046),
  pre-existing. Not yours to fix, but do not build the `decision` comparator on the assumption
  that the item tests what its difficulty label claims.
- **The T6 baseline regeneration** (D-043) — Phase 6.
