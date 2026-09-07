# Prompt 06 — Phase 5: output-contract hygiene, and comparator bindings

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Phases 0, 1, C2, 2, 3 and 4 are complete and merged; Phase 4's merge commit is `9e288d9`
and this prompt sits on top of it. **This phase has two tracks that share a number and
nothing else.** Track A edits eleven templates and is gated by T4. Track B edits no
template at all and is gated by cross-pairing. **Read §"How the two tracks relate" before
you plan anything** — merging them is the specific way this phase fails.

---

I need you to implement **Phase 5** of the template redesign: **Track A**, output-contract
hygiene across eleven templates, and **Track B**, the comparator bindings Phase 4 specified
but did not deliver, plus the two instruments that make a binding checkable.

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

For the full picture (optional, but read it if a judgement call turns on what the benchmark
is claiming):

- `docs/_ARR_May__EngTrace.txt` — the current paper. **Section 4** is the evaluation
  framework, **Appendix M** is the framework validation.
- `docs/EngTrace_Rebuttal_Jul2026.txt` — the most recent reviewer objections and our
  responses. **These promises are the real to-do list behind this work.**

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — these are
  an **LLM annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not
  flag it. (You *will* be working in `error_analysis_annotation/samples/`, which is different
  and is real archived model output.)
- `evaluation/` is a **separate track** with known defects (D-003). **Do not fix them.** You
  are extending its replacement.

## How the two tracks relate — read this before planning

They share a phase number because they land together. They share nothing else:

| | Track A | Track B |
|---|---|---|
| edits | 11 templates in `data/templates/` | `tests/comparators/`, `template_inventory.csv` |
| gate | T4 passes 150/150, T6 unchanged | zero false accepts on gold×gold |
| reviewer | **A** — correctness | **E** — comparator adversary |
| risk | changing an item pool by accident | binding a template wrongly and not noticing |

**Do not commission one reviewer for both.** R6.2 requires a measurement-ownership register
with no number owned twice, and the Phase 0 precedent is that a bundled brief stalls and the
gate it guards silently becomes a formality. If you find yourself writing a single review
brief covering T4 and cross-pairing, you have mis-scoped it.

**Track A is small and Track B is not.** If the phase cannot carry both, **Track A ships and
Track B re-scopes** — Track A unblocks structural migration for every later phase, and its
absence is felt corpus-wide. Say so explicitly rather than delivering both thinly.

## Read these first, in this order

1. **[`docs/re-implementation-sep/template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** —
   §0 governing principles **P1–P6**, all of **Phase 5** (both tracks), and the review
   protocol **R0–R6**. Read R0–R6 properly before you dispatch anything; it is short and it
   is the part most often skipped. Note **SPEC-CHANGE 11, 12 and 13** in the revision log —
   all three are Phase 4's.
2. **[`docs/re-implementation-sep/phase4_summary.md`](../re-implementation-sep/phase4_summary.md)** —
   **§8, §9 and §12 are written for you.** §8 groups nineteen errors into three shapes and
   two of the three are shapes you are about to be exposed to; §9 is what eight review rounds
   cost and bought; §12 states the stopping rule.
3. **[`docs/re-implementation-sep/phase4_comparators.md`](../re-implementation-sep/phase4_comparators.md)** —
   D4.1, the contract you are binding templates to. **§2 is the pair of rules that must hold
   for every binding you add**, §3 is `multipart`, §5 is what it does *not* settle.
4. **[`docs/re-implementation-sep/DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** —
   append-only; add entries, never edit existing ones. **D-058** is the finding that creates
   Track B. Also **D-047** (`multipart` composes, it is not a kind), **D-049** (the
   three-valued verdict), **D-052** (units by declaration), **D-056** (hedges advisory),
   **D-034** (a test may not carry its own answer key), **D-024/D-026** (the evidence-rate
   rule).
5. **[`docs/re-implementation-sep/reviews/phase4_reviewer_e_comparator.md`](../re-implementation-sep/reviews/phase4_reviewer_e_comparator.md)** —
   four rounds. **Read §"Round 4 — Task 3" before you write a line of comparator code**, and
   read F0 before you quote any precision figure.

## What is actually here, measured

**Verify these yourself before acting — every prior phase found the spec's list was a
starting point rather than an inventory.** Phase 1 found its prescribed fix insufficient,
Phase C2 found one recorded defect was four, Phase 2 found a template failing T1 on 91% of
instances the spec did not mention, Phase 3 found two display-tie populations where the brief
listed one, and Phase 4 found T6 failing on two of its four templates where the brief said
one.

Corpus health at `master` (`9e288d9`) — **this is the baseline you work against, not a
to-do list.** Re-measure it in a worktree rather than trusting the table:

| Check | Failing | Note |
|---|---:|---|
| T1 printed-arithmetic closure | 29/150 | 97 more sit on the rounding boundary (marginals) |
| T2 round-trip oracle | 0/150 | only templates with an oracle are checked |
| T3 determinism | 0/150 | |
| T4 output contract | **3/150** | **Track A's stated gate — and it does not gate most of Track A. See below.** |
| T5 binding / rounding | 66/150 | |
| T6 distribution | 142/150 | committed baseline is stale corpus-wide (D-043); Phase 6 owns regeneration |
| T7 invariant asserts | 83/150 | advisory corpus-wide; a gate on templates a phase edits |

**Track A's gate covers 3 of its 11 templates, and this is the first thing to fix.**
T4 fails on exactly `euclidean_distance_binary`, `cd_dc_system_analysis` and
`finite_convolution` — the spec's *malformed `**Step N:**`* row. The other three defect
classes in that table (5 templates using `**Final Answer**`, 2 with malformed complex
numbers, 1 with an unreachable `elif`) **are invisible to T4**, so "T4 passes 150/150" would
be true with 8 of the 11 edits unmade or made wrongly. The spec's 11 and T4's 3 do not
contradict each other; the *gate* is the thing that is wrong.

Decide how to close that — extend T4, add a check, or state per-defect-class acceptance
evidence — and record it as a `SPEC-CHANGE` against Phase 5 Track A's exit gate. **Do not
simply satisfy the gate as written.** Phase 0 found two checks that were green while
measuring nothing, and this is a third in advance.

**Track B's measured starting point** (D-058, re-derive it):

| | |
|---|---|
| templates bound to the comparator contract | **4 of 150** |
| gold×gold pairs available, all 150 × 12 instances | **19,668** |
| archive×gold pairs available from the 2,200-trace archive | **22,982** (11,384 `scalar`, 5,511 `multipart`) |
| gold answers carrying a unit token | **87 of 150** |
| templates exposed to N1 (leading incidental number) | **58 of 150** |
| templates measured to produce an N1 false accept, 12 instances | **3** |

**That last pair of rows is the shape of every number in this phase**, and getting it wrong
is the easiest mistake available to you. 58 is the count *exposed*; 3 is the count *measured*
to fail. I first reported 58 as the defect size and it is not. **Never quote an exposure count
as a defect count.**

## Track A — output-contract hygiene

Eleven templates, four defect classes, listed in spec §Phase 5 Track A. `cd_dc_system_analysis`
is the urgent one: a strict marker parser drops **Step 3 in 100% of instances**, and Step 3
computes the answer.

**One decision to make and record.** Five templates use `**Final Answer**` / `**Final
Answers:**` instead of `**Answer:**`. Normalise the templates, **or** widen the accepted
marker set corpus-wide — decide once, apply everywhere, and record it. Note that
`tests/comparators/normalize.py::ANSWER_MARKERS` already accepts several variants and its
priority order is load-bearing; whichever way you decide, the two must agree.

`levenspiel_plot_interpretation` is also in Phase 2's scope — check whether that edit landed
before touching the file.

## Track B — bindings, and the instruments that make a binding checkable

**The finding that creates this track.** Phase 4 specified the contract and bound four
templates. Gold×gold cross-pairing over all 150 found three defects Phase 4's own corpora
could not see, because `numeric` and `check` had **zero archived negative instances**
(Reviewer E, F0):

**N1 — `parse_number` reads the first number in the answer span, which is often not the
answer.** Confirmed:

```
cand: The volumetric flow rate of Engine Oil (SAE 50) ... is 0.001183 m^3/s.
gold: The volumetric flow rate of Engine Oil (SAE 50) ... is 0.002738 m^3/s.
-> MATCH        both canonicalise to 50
```

The same defect reads a temperature (541) and a chemical-formula subscript (the 4 of C₄H₁₀)
on `gas_viscosity_kinetic_theory`.

**N2 — an empty parse compares equal to an empty parse.** `_as_polynomial` returns an empty
coefficient map when it recovers nothing, and `{} == {}` is a `MATCH`. On
`autocorrelation_rect_pulse` that is 128 of 132 pairs. D4.1 §4.4 already states the rule this
violates: *failure is `UNRESOLVED`, never `MATCH`*.

**N3 — no bindings exist for 146 templates**, and it is larger than N1 and N2 together.
`euclidean_distance_binary` matches 132/132 because a `multipart` answer read by a single-part
comparator returns the `1` from `s1`.

**The reframing that matters: the comparator does not need more rules, it needs bindings, and
each binding needs evidence that it discriminates.** Phase 4 spent eight review rounds adding
rules to a natural-language surface and closed by *reducing* the mechanism's scope. Do not
repeat that. If you find yourself adding a rule to make one template pass, stop and ask
whether a binding would do it instead.

### D5.6 is the one to get right, and it has a named failure mode

Fixing N1 needs an extraction rule. "First number" is wrong. **"Last number" is also wrong** —
it breaks `check` answers like *"18.4 mm, less than the 25 mm limit"*, where the last number
is the limit. Neither is obviously right, and picking one from a single example is the error
shape `phase4_summary.md` §8 records **six times**.

**So: implement at least three candidate rules and score each against D5.4 and D5.5, then
adopt the winner and record the losers' numbers.** Candidates worth including — you may add
others, and you should:

1. the last number in the span;
2. the number adjacent to a declared unit (D5.10 gives you 87 of these);
3. `UNRESOLVED` when the span holds more than one number and nothing disambiguates.

Candidate 3 will lose on decided-rate and may still be right; report all three rather than
only the winner.

## Deliverables

**Track A:** D5.1 the eleven template edits · D5.2 corpus-wide marker-conformance scan across
all 150 proving zero remaining violations · D5.3 the marker-set decision record.

**Track B:** D5.4 gold×gold cross-pairing as a standing check · D5.5 archive×gold
cross-pairing, Reviewer E's instrument promoted from review artefact · D5.6 N1 fixed **by
measurement**, losing candidates recorded · D5.7 N2 fixed on every path · D5.8 bindings
declared and cross-pair-validated · D5.9 per-kind negative counts beside every precision
figure · D5.10 units declared for the 87 · D5.11 the Phase 4 residuals R4-11, RB4-1, R4-12.

**Both:** `phase5_summary.md` in the shape of `phase4_summary.md` — what changed, the R4
triage tables, **both** exit gates, residual risk, what Phase 6 inherits, **and your own
errors**. New `DECISIONS.md` entries for every decision and reversal, including the ones you
reverse mid-phase. Review reports under `docs/re-implementation-sep/reviews/`.

**If the phase cannot carry both tracks, Track A ships and Track B re-scopes.** Say so and
re-scope explicitly rather than delivering both thinly.

## Verification

```
python -m tests.template_integrity.run --checks all          # corpus, ~25s
python -m tests.template_integrity.run --templates a,b --checks all
python -m tests.trace_schema.audit_3_8                       # must stay clean
python -m tests.comparators.score                            # archive + D4.4
python -m tests.comparators.reviewer_battery                 # 81 reviewer cases
python -m tests.comparators.recall_corpus                    # 507 cases
python -m tests.comparators.derive_vocabulary                # SUPPORT + commitment census
python -m tests.trace_schema.candidate_7 docs/re-implementation-sep/phase4_conformance/candidates.json
```

**`reviewer_battery` and `recall_corpus` run under two hedge policies** and the default is
`advisory` (D-056). `ENGTRACE_HEDGE_POLICY=enforce` checks the reviewers' cases against the
policy they were argued for. **Both must stay green under both policies**; a change that is
green under one and red under the other is a finding, not a pass.

**Verify the corpus baseline yourself against a `git worktree` at `master`** rather than
trusting any table. Use a short worktree path such as `C:/wtm`; the repo contains paths long
enough to fail a checkout under a temp directory. **Remove the worktree before you merge** —
it holds the `master` branch and the merge will refuse.

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and will
crash on `→`/`Σ`/`²`/`·`.

**Track A owes a before/after instance dump.** `tests/template_integrity/phase4_instance_dump.py`
is the pattern to copy: it takes the tree root as an argument and refuses to run if a template
resolves outside it. **Run each tree as its own process.** An in-process reload resolves both
sides to the same already-imported `data.templates.*` modules and reports every instance
identical — documented in `instance_dump.py`, and it has now caught four people.

**Track A changes emitted text, so T6 will move and that is expected for the marker edits.**
T6's committed baseline is stale corpus-wide (D-043) and Phase 6 owns regeneration — **do not
regenerate it to make a gate green.** The before/after dump replaces it, exactly as Phases 3
and 4 did.

**A caution about T6 you will hit and should not chase.** T6's answer extractor cannot see a
non-scalar answer at all: it returns a sequence's first element and a polynomial's leading
coefficient (`phase4_summary.md` §6). Its breach lines on those templates compare parse
artifacts. That is recorded as R4-5, it is Phase 6's, and it is *not* something to fix here.

**State the smallest rate your run can resolve and check it against the rate you need to
exclude.** N samples cannot resolve a defect rarer than ~3/N (D-024, D-026). In this phase it
binds hardest on D5.8: a binding validated on 12 instances resolves nothing rarer than 25%.

## Reviews — read R0–R6 before dispatching any

Phase 5 requires **two independent reviewers, one per track, working in parallel and in
isolation from each other.**

- **Reviewer A — correctness (Track A).** *The gate.* Re-runs T4 across all 150 templates,
  not just the ones touched, and confirms the item pool did not move. **Do not ask A about
  the comparator.**
- **Reviewer E — comparator adversary (Track B).** *The gate.* Given the bindings and nothing
  else: **find a template whose binding accepts an answer to a different instance of the same
  template.** Draw from gold×gold and archive×gold, not from invented examples. **Do not ask
  E to re-run T4.**

**Scope every review to R6 before dispatching it.** One mandatory gate task, a stated time
box, tooling supplied as working code, everything already settled fenced off explicitly, and a
measurement-ownership register with **no number commissioned twice**.

**Independence, concretely:**

- **R1.1 — work without sight of the implementer's reasoning.** Supply code comments, commit
  messages and summaries as **claims under test**, not as background.
- **Give each reviewer a frozen ref — a commit SHA, not a branch name.** Then **do not run
  `git commit` at all between dispatching a review and its filing.** Phase 3 broke this three
  times and Phase 4 kept it by adopting exactly that mechanical rule. **Note that Phase 4 also
  had to stop *writing to the working tree*, not merely committing** — a reviewer reading the
  tree at a frozen SHA sees uncommitted edits too.
- Tell reviewers to read files with `git show <sha>:<path>` and **never `git show <sha> --stat`**
  — the latter prints the commit body and contaminated a Phase 3 reviewer.
- **Do not let the two reviewers see each other's briefs or reports.**

**Budget for more than one round, and know what "done" looks like.** Phase 3's specification
took six rounds. **Phase 4 took four rounds per reviewer and did not close by exhausting the
findings — it closed when Reviewer E characterised the surface as unbounded and recommended
reducing the mechanism's scope.** Track B is a *binding* problem rather than a
natural-language one, so it should converge faster; if it does not, that is itself the
finding, and the stopping rule is the same:

> **When a mechanism's findings become mutations of cases it already passes, and its census
> shows no observed instances, reduce its scope rather than add a rule.**

Reports go to `docs/re-implementation-sep/reviews/phase5_<reviewer>.md` in the R2 five-section
structure. **Every §5 suggestion must be triaged (R4) before the phase closes; an untriaged
suggestion blocks the gate exactly as a CONFIRMED finding does.** A `SPEC-CHANGE` **amends
`template_redesign_spec.md`**, and an `ADOPT-PHASE-<N>` **edits that phase's deliverable
list**. Writing the disposition in your summary is not discharging it.

**State honestly which findings you rejected and why.** A `REJECT` with a written reason is a
legitimate and respected outcome; a finding quietly absent from the triage table is not.

**Confirm your reports are actually committed, not merely written** — a bare `reviews` entry
in `.gitignore` once silently untracked all eight prior reports.

## Git

- Branch **`redesign/phase5-hygiene-bindings`** off `master`. Do not work directly on `master`.
- Logical commits, not one lump. **Commit each reviewer's report unmodified and before any
  fix**, so the record is what they wrote rather than what survived your response.
- **Do not push.** `master` has never been published and `origin/master` is not fetched in
  this working copy, so *check the actual divergence yourself before believing any number for
  it* — including this one. Publishing is a separate decision and not this phase's.
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
- Merge to `master` with `--no-ff` only after both exit gates pass, or record explicitly what
  is left open.

## Exit gate

**Track A**

- [ ] T4 passes **150/150**
- [ ] **Acceptance evidence exists for all four defect classes, not just the marker class
      T4 sees** — and the gate is amended so a later phase cannot repeat the mistake
- [ ] Marker-set decision recorded, and `normalize.py::ANSWER_MARKERS` agrees with it
- [ ] T6 movement explained by the before/after dump; **baseline not regenerated**
- [ ] Item-pool impact stated, with the dump as evidence

**Track B**

- [ ] D5.4 and D5.5 land as runnable standing checks, pair counts reported
- [ ] **Zero false accepts on gold×gold across all 150**, or every remainder named with a reason
- [ ] N1 fixed by a rule **chosen on measured evidence**, losing candidates' numbers recorded
- [ ] N2 fixed: no code path returns `MATCH` from an unrecoverable parse
- [ ] **Per-kind negative-instance counts published beside every precision figure**
- [ ] Bindings declared and cross-pair-validated for every `classification`, `symbolic`,
      `vector`, `array` and `multipart` template
- [ ] `reviewer_battery`, `recall_corpus`, `score` and `candidate_7` green **under both hedge
      policies**; `audit_3_8` clean; no corpus regression, measured not assumed

**Both**

- [ ] Reviewers A and E filed, **every finding and every §5 suggestion triaged (R2, R4)**, and
      every `SPEC-CHANGE` / `ADOPT-PHASE-N` actioned in the documents they name

## Constraints and cautions

- **Do not modify a check or an oracle to make something pass.** If a check looks wrong, that
  is a finding and a `SPEC-CHANGE`, raised explicitly. Phase 0 found two checks that were
  green while measuring nothing.
- **Do not build a binding to agree with the existing parser** (D-003, D-034). Reproducing its
  defects is the same error as a test carrying its own answer key — and N1 is that error
  already committed once, inside the comparator built to replace it.
- **A corpus that certifies a rule set must contain cases from outside it.** Reviewer B
  diagnosed this in Phase 4 round 1 — *"the corpus samples the inside of the list it
  certifies"* — and the implementer committed it again in round 2 by shipping a recall corpus
  with no negative control. **Every corpus you add gets a negative control before you quote a
  number from it.**
- **P6 is a real constraint.** Track A changes emitted text; Track B changes what counts as a
  correct answer for up to 146 templates, silently and corpus-wide. The second is the larger
  P6 surface even though it edits no item.
- **Never quote an exposure count as a defect count.** See §"What is actually here".
- **Where you are uncertain, say so and explain what evidence would settle it.** Do not smooth
  over a gap to make a gate look clean.
- **Report your own errors in the summary.** Phase 4 recorded nineteen; thirteen shared three
  shapes. If your summary records none, that is not a clean phase, it is an unexamined one.

## Open items inherited — carry forward, fix only if they block your gate

- **D-003 is still open**: do the raw `inference_results/` generations still exist? It gates
  any corrected results table.
- **The corpus-wide sweeps from Phase 2** — templates combining a transcendental with a
  fixed-decimal print; the f-string regex `\{[a-z_]+ *[-+*/] *[0-9.]+\}`; two coexisting
  rounding conventions. None swept.
- **`_froude_capped_slope`'s exactness pathology affects three templates** in
  `civil_engineering/water_resources/uniform_flow.py`, not one (D-045).
- **`line_balancing_heuristic` is ~90% shortcuttable** (D-046) and **both Phase 4
  classification items are 100% shortcuttable** (D-057). Pre-existing; `blind_guess_floor` and
  `surface_model_heldout` are columns in `template_inventory.csv` for the five templates where
  the statistic means anything. **Not yours to fix**, but do not build a binding on the
  assumption that an item tests what its difficulty label claims.
- **The T6 baseline regeneration** (D-043) and **T6's blindness to non-scalar answers**
  (R4-5) — Phase 6.
- **`sympy` is a dependency and is not in `requirements.txt`** (D-053).
