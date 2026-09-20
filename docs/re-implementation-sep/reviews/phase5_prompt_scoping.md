# Phase 5 prompt — scoping and design review

**Artefact:** `docs/prompts/06-phase5-contract-hygiene-and-bindings.md`
**Frozen ref:** `bbccae510b8eac64f45c3e2b5e92a5be8b3f067b`
**Scope:** scoping and design only. Every count, SHA, path and command in the brief is
another reviewer's; where one is mentioned below it is because a *structure* depends on it,
never because I re-derived it.
**Box:** 35 minutes, filed at the box.

---

## 1. Verdict

**It will produce Track A and a partial, over-rejecting Track B.** The two-track split is
right in principle but the brief's central claim — that the tracks "share a phase number and
nothing else" — is false by its own text, and the one place they *do* couple is fenced off
from both reviewers; Track B has no stated budget, no threshold, no objective function for its
one measured decision, and a stopping rule imported from a problem it is not solving.

---

## 2. Assessment against R6

| Anti-pattern | Verdict | Evidence |
|---|---|---|
| **Duplicate commissioning** | **Clean.** The best-handled rule in the brief. | Explicit negative fences in both directions — *"Do not ask A about the comparator"*, *"Do not ask E to re-run T4"* — plus a `reviewer` row in the relate-table and the R6.2 register demanded by name. No number is visibly commissioned twice. |
| **Kitchen-sink brief** | **Violated, twice over.** | (a) The *session* carries eleven deliverables across two unrelated surfaces under two gates. R6.1's "one mandatory gate property per review" is honoured per reviewer, but the session itself has two gates and no arbitration rule beyond a contingency. (b) Reviewer A's mandatory task is *two* properties with two instruments: "re-runs T4 across all 150" **and** "confirms the item pool did not move" (brief §Reviews). |
| **Instrument-building** | **Mostly clean, one hole.** | Reviewer E inherits D5.4/D5.5 as standing checks — the right move, and the brief promotes E's own round-1 instrument rather than asking E to rebuild it. Reviewer A gets T4, which exists. The hole: nobody is given an instrument for the four-defect-class acceptance evidence, because that instrument does not exist yet and the brief does not say who builds it or what it is. |
| **Re-deriving settled numbers** | **Inverted deliberately, and it is defensible.** | R6.4 says fence off what is settled; the brief instead says *"Verify these yourself before acting"* and *"Re-measure it in a worktree rather than trusting the table"* for the whole baseline. Given five prior phases where the spec's list was a starting point, this is a considered trade, not an oversight. But it is a trade: nothing at all is fenced off as settled, so R6.4's required line is absent for the *implementer*. For the reviewers it is not supplied either — the brief tells the session to write the fences without listing a single one. |
| **Unbounded probing** | **Violated.** | R6.1: *"State the box in the brief as a number, not as 'be efficient'."* The brief tells the next session to state a time box (§Reviews) and never states one anywhere — not for the phase, not for either track, not for either review. The word "hours" does not appear. Prompt 05 had the same gap and Phase 4 ran eight review rounds against a plan that budgeted one. |
| **Pre-dispatch checklist** | Recited, not instantiated. | The brief reproduces the checklist as prose (one mandatory task, box, tooling, fences, register) but supplies zero of the six items in concrete form. See F6. |

---

## 3. Findings

### F1 — CONFIRMED. "A phase number and nothing else" is false; the tracks share `normalize.py` and they share the gold corpus.

The brief's own Track A section: *"Note that `tests/comparators/normalize.py::ANSWER_MARKERS`
already accepts several variants and its priority order is load-bearing; whichever way you
decide, the two must agree."* `tests/comparators/` is, by the same brief's table, **Track B's
edit surface**. So D5.3 (Track A's marker decision) either edits Track B's file or constrains
it.

It is worse than a shared file. `ANSWER_MARKERS` is what selects the answer span; the answer
span is what `parse_number` reads; `parse_number` is N1. And Track A rewrites the emitted gold
text of eleven templates, which is the text D5.4's gold×gold cross-pairing consumes for those
eleven. **Track A's output is an input to Track B's measurement.**

*Consequence:* run in parallel as the brief invites ("they land together"), D5.6's candidate
bake-off is scored on a corpus that changes under it, and the losing candidates' recorded
numbers — an explicit deliverable — are recorded against a corpus that no longer exists at
merge. Re-running after Track A lands is cheap; *not noticing* is the failure, and nothing in
the brief prompts the check.

### F2 — CONFIRMED. The seam between the tracks is unowned by both reviewers.

A is told not to look at the comparator. E is told not to re-run T4. The marker decision's
agreement with `ANSWER_MARKERS` — the one artefact spanning both tracks — sits exactly in the
gap, and it is an exit-gate item ("Marker-set decision recorded, and `normalize.py::ANSWER_MARKERS`
agrees with it") with no owner in the measurement-ownership register the brief demands.

*Consequence:* the gate item is self-certified by the implementer. That is the Phase 0 pattern
the brief itself cites — a check that is green while measuring nothing — reproduced in the
brief's own gate.

### F3 — CONFIRMED. Reviewer A's commission re-instates the gate the brief spends a page diagnosing as broken.

§"What is actually here" is the brief's strongest passage: T4 sees 3 of the 11 edits, so
*"T4 passes 150/150 would be true with 8 of the 11 edits unmade or made wrongly"*, and the
session is told to amend the gate and record a `SPEC-CHANGE`. Two pages later, Reviewer A's
mandatory task is *"Re-runs T4 across all 150 templates"* — the unamended gate — and the Track A
exit gate's **first** line is still `T4 passes 150/150`, with the amendment demoted to the
second bullet.

*Consequence:* the reviewer certifies 3 of 11 edits. The other 8 are certified by whatever the
implementer invents and then reviews itself. The diagnosis is excellent and the dispatch does
not act on it.

### F4 — CONFIRMED. Reviewer E's mandatory task drops the false-reject half that Phase 4 gave it, and D5.6 has no objective function.

Prompt 05: *"find a correct model answer they reject, **and** an incorrect one they accept."*
Prompt 06: *"find a template whose binding accepts an answer to a different instance of the
same template."* Accepts only. The Track B exit gate is likewise single-sided: *zero false
accepts on gold×gold*.

Now read D5.6 against that. The brief says implement ≥3 candidate rules and *"score each against
D5.4 and D5.5, then adopt the winner"* — and never says what is maximised. The two obvious
objectives are degenerate in opposite directions:

- Minimise false accepts on cross-pairs → **candidate 3 wins trivially at its limit**: a rule
  that returns `UNRESOLVED` whenever a span holds more than one number has a false-accept rate
  of zero and decides nothing.
- Maximise decided rate → **"first number" wins**, which is N1 itself.

The brief half-sees this (*"Candidate 3 will lose on decided-rate and may still be right"*) but
supplies no combination, no tie-break, no threshold, and — critically — **no corpus that
measures the decided rate on answers that ought to be decided.** Gold×gold same-instance pairs
are byte-identical text; decided-rate measured there is uninformative. The only corpus in the
repo that constrains numeric over-rejection is D4.4's 24 hand-written `numeric` cases, held
green by the `score` gate; D5.6's procedure does not name it, and 24 cases resolve nothing
rarer than ~12% under the brief's own D-024/D-026 rule.

*Consequence:* a session can follow D5.6 exactly, pick candidate 3, publish three candidates'
numbers, pass a single-sided gate with an adversary commissioned only for accepts — and ship a
comparator that resolves less than the one it replaced. Every Phase 4 reviewer-E finding class
that mattered (13 of 20) was **over-rejection**; the brief removes the only party who was
looking for it.

### F5 — CONFIRMED. D5.8 is an open job wearing a bounded one's clothes.

The spec requires each binding to *"clear a cross-pairing threshold before it counts as bound"*.
The brief carries the requirement (*"cross-pair-validated"*, exit gate: for every
`classification`, `symbolic`, `vector`, `array` and `multipart` template) and **never states the
threshold**. It states the obstacle instead — *"a binding validated on 12 instances resolves
nothing rarer than 25%"* — and stops there.

Two further omissions make it open rather than merely under-specified:

- **The lever is never named.** Instances per template is a free parameter; gold generation is
  deterministic and cheap. Raising it is the only way to move the 25% floor, and the brief does
  not mention that it can be raised.
- **The spec's ordering is dropped.** Spec D5.8: *"Start with the 5 `classification`, 9
  `symbolic`, 15 `vector`/`array` and 32 `multipart`."* The brief's exit gate flattens these to
  a single conjunction with no priority, so a session that runs short has no principled way to
  choose what to ship.

*Consequence:* ~61 bindings × an unstated evidential bar. This is the deliverable most likely to
be delivered thinly and reported as complete.

### F6 — CONFIRMED. The brief demands the R6 apparatus from the next session and supplies none of it.

*"Scope every review to R6 before dispatching it. One mandatory gate task, a stated time box,
tooling supplied as working code, everything already settled fenced off explicitly, and a
measurement-ownership register with no number commissioned twice."* Against that list the brief
provides: the two mandatory tasks (yes, and well fenced — see §4), and nothing else. No box
number, no register rows, no settled-list, no tooling snippet. Every one of these is harder
here than in Phase 4 because the corpus baseline is explicitly *not* trusted, so "what is
settled" is genuinely non-obvious and the brief passes that judgement downstream unaided.

*Consequence:* PLAUSIBLE-level, but the Phase 0 precedent the brief itself cites is a review
that stalled and produced nothing because its brief was written under time pressure at the end
of an implementation phase. That is exactly when these two briefs will be written.

### F7 — CONFIRMED. The stopping rule is a conclusion quoted without its instrument.

> *"When a mechanism's findings become mutations of cases it already passes, and its census
> shows no observed instances, reduce its scope rather than add a rule."*

Two clauses, neither operational as written, and the round-4 report shows why. Reviewer E did
not apply this rule; E **manufactured the evidence for it**, in three moves the brief does not
transmit:

1. **A discovery-cost series across rounds** — "6 findings per session in round 2, 10 in round
   3, 5 in 25 minutes in round 4" — which is what "mutations of cases it already passes"
   actually cashes out to. Nobody is told to record findings-per-reviewer-hour, so the series
   cannot exist in Phase 5.
2. **A census that separated one object into two mechanisms with opposite evidence** (hedge
   markers 2/2200 vs subordinator openers 43). E states the split *"was invisible while the
   commitment machinery was discussed as one object."* The brief says "its census" as though a
   census is a thing one has.
3. **An ablation** — stub the layer, re-run the corpora, price the deletion at 0 positive
   recall / 0 archive verdicts / 30 synthetic negatives. This is the measurement that ended the
   phase, and **the brief does not mention ablation at all.**

Worse, the rule is scoped to *a mechanism with a census*, and the brief says in the same
paragraph that Track B *"is a binding problem rather than a natural-language one"*. A binding
has no census of observed instances; scope reduction for 61 bindings means "bind fewer
templates", which the rule does not say and the exit gate forbids.

*Consequence:* the brief has quoted the lesson and not internalised it — the specific charge
the review task raised, and it is upheld. Phase 4 §8 records the same failure at one level up:
*"Writing the lesson down did not prevent it; building a check that runs did."*

### F8 — CONFIRMED. There is no stop trigger, and no budget for it to be measured against.

The only stop signal is *"If the phase cannot carry both tracks, Track A ships and Track B
re-scopes"* — stated twice, with no condition attached to "cannot". No hours, no round count, no
deliverable-completion trigger. The spec's 15–25 h + 8–10 h for Track B is never repeated in
the brief, so the session's only budget signal is a number it must go and find in a 891-line
spec.

On credibility: Phase 4 budgeted 12–20 h + 8–10 h review for **four templates and no code
edits**, and took **eight review rounds, 34 findings, 19 self-recorded errors**. Track B is
61 bindings + 87 unit declarations + two new standing instruments + two parser defects + three
Phase 4 residuals + a three-way rule bake-off, and its review is a *comparator adversary that
found five blocking findings in 25 minutes* against a mature suite. The brief's argument for
faster convergence — bindings, not natural language — is contradicted by its own D5.6, which is
a rule over prose, is the deliverable the brief singles out as *"the one to get right"*, and is
where every Phase 4 error of Shape 1 lived.

*Consequence:* PLAUSIBLE but strongly evidenced — Track B overruns, and the brief's fallback
fires late, after Track A has already been merged and the sunk cost argues for finishing.

### F9 — CONFIRMED. Two Phase-5-owned inheritances are dropped, one of them contradicted.

- `phase4_summary.md` §12: *"**Phase 5 owns the supporting-quantity remedy**, where the item can
  change with its gold"* (D-057, the two 100%-shortcuttable classification items). The brief's
  Open-items section says of the same D-057: *"**Not yours to fix**."* A reassignment of an
  inherited deliverable, with no `SPEC-CHANGE` and no R4 disposition — the exact discharge
  failure the brief's own §Reviews warns about ("Writing the disposition in your summary is not
  discharging it").
- Residual **R4-1** — *"the symbolic comparator is not validated against a representative
  sample", severity **high**, owner **Phase 5/6*** — is absent from the brief. D5.11 carries
  R4-11, RB4-1 and R4-12 (all medium) and not the high one. Phase 5 binds 9 `symbolic`
  templates, so this is the residual most directly in scope.

*Consequence:* the highest-severity open item touching this phase's work is invisible to the
session, and a deliverable is silently cancelled.

### F10 — CONFIRMED. The Track A gate contradicts itself on T6, and only one of the two amendments is commissioned.

The relate-table says Track A's gate is *"T4 passes 150/150, **T6 unchanged**"* (following the
spec: *"these are formatting-only, so any distribution movement indicates an unintended
change"*). §Verification then says *"Track A changes emitted text, so **T6 will move** and that
is expected"*, and the exit gate reads *"T6 movement explained by the before/after dump"*.

The brief is right on the substance and it is amending a second spec gate clause — but it
demands a `SPEC-CHANGE` only for the T4 clause. Meanwhile the first thing a planning session
reads is the table, which states the superseded rule.

*Consequence:* low; likely caught. But it is a `SPEC-CHANGE` that will not be filed because the
brief does not ask for it, and the spec keeps a gate clause that is now known to be wrong.

### F11 — PLAUSIBLE. Reviewer E is rostered for Phase 4 only.

Spec §Roster: A → phases 1, 2, 5; E → phase 4. Using E for Phase 5 Track B is obviously the
right call, and it is what the spec's own Phase 5 Track B section assumes — but the roster table
was not updated when Track B was added, and the brief neither notes the discrepancy nor asks for
the amendment. Cheap to fix, and the brief is otherwise scrupulous about this class of thing.

### F12 — PLAUSIBLE. Nothing tells the session to fan out.

Spec E1: *"A phase's wall-clock should approach its slowest single item, not the sum of its
items."* Sixty-one bindings across five branches is the most parallelisable unit of work in the
entire plan — more so than Phase 1's nine templates, which the spec names as the exemplar — and
the brief never mentions E1, concurrency, or per-branch division. Combined with F8's absent
budget, the default is serial grind.

---

## 4. What the brief gets right

Specifically, because most of it is right.

- **The T4 gate diagnosis (§"What is actually here") is the best thing in the document and has
  no precedent in prompt 05.** It works out that the phase's own stated gate covers 3 of 11
  edits, says *"the spec's 11 and T4's 3 do not contradict each other; the gate is the thing
  that is wrong"*, forbids the easy path (*"Do not simply satisfy the gate as written"*), and
  routes the fix through `SPEC-CHANGE`. Finding a broken gate *before* dispatch, and saying so
  in the brief, is the review protocol working upstream of the reviewers.
- **The exposure/defect distinction, with the author's own error attached.** *"58 is the count
  exposed; 3 is the count measured to fail. I first reported 58 as the defect size and it is
  not."* Naming one's own live error in the brief that inherits it is worth more than the
  warning.
- **Reviewer isolation is handled better than in any prior brief in the series.** Frozen SHA not
  branch name; no `git commit` between dispatch and filing; and the Phase 4 upgrade the earlier
  brief lacked — *"Phase 4 also had to stop **writing to the working tree**, not merely
  committing"*. Plus `git show <sha>:<path>` and the `--stat` prohibition, and mutual blindness
  between A and E.
- **The negative fences on the two reviewers are exactly R6.2.** Two lines — *"Do not ask A
  about the comparator"*, *"Do not ask E to re-run T4"* — do more anti-duplication work than the
  paragraph of protocol around them.
- **Do-not-fix boundaries are drawn cleanly and repeatedly**: `evaluation/`, the T6 baseline,
  R4-5's non-scalar blindness, the annotation-pilot false alarm. A session will not waste a day
  on any of them.
- **The instance-dump trap is transmitted as a mechanical rule with a failure count attached**
  — *"An in-process reload resolves both sides to the same already-imported modules... it has
  now caught four people."* Run-each-tree-as-its-own-process is actionable in a way that "be
  careful with imports" is not.
- **D5.6 is framed as a bake-off rather than a fix.** Even without an objective function (F4),
  "implement three, publish the losers' numbers" is structurally right and directly targets the
  error shape Phase 4 committed six times. *"Report all three rather than only the winner"* is
  the correct instinct.
- **The N1-as-D-034 observation is sharp**: building a binding to agree with the existing parser
  is *"the same error as a test carrying its own answer key — and N1 is that error already
  committed once, inside the comparator built to replace it."*
- **The re-measure-the-baseline instruction is justified with five phases of evidence** rather
  than asserted, which is what makes it a defensible inversion of R6.4 rather than a lapse.

---

## 5. Rewrite recommendations, ranked

**1. Split the phase. Track A becomes Phase 5 and merges; Track B becomes Phase 7.**
The brief already knows this — *"Track A ships and Track B re-scopes"* is its contingency, F1
shows Track A's output is Track B's input, and F2 shows the shared seam has no reviewer. Make
the order the plan. Replace the relate-table's framing with:

> Track A lands and merges first. Its marker decision (D5.3) changes the gold text that Track
> B's cross-pairing consumes and constrains `ANSWER_MARKERS`, so Track B's corpora are built
> **after** Track A merges. If you run them concurrently, every D5.4/D5.6 number is provisional
> until Track A's merge, and you must re-run rather than re-read them.

If they must share a session, at minimum delete *"They share nothing else"* — it is false and it
is load-bearing for the reviewer split.

**2. Give D5.6 an objective, a held-out corpus, and a tie-break.** Replace *"score each against
D5.4 and D5.5, then adopt the winner"* with:

> Score every candidate on **two** axes and publish both for all three:
> (a) **false accepts** on gold×gold and archive×gold cross-pairs (the negatives);
> (b) **decided rate on answers that ought to be decided** — D4.4's 24 `numeric` cases plus a
> new positive set you build from real archive spans for templates whose gold answer carries a
> declared unit (D5.10 gives you 87).
> Adopt the rule with the highest decided rate among those with **zero** false accepts; if none
> reaches zero, say so and adopt on the stated trade rather than silently. **Hold out 20% of
> both corpora, chosen before you look at any candidate, and report the winner's numbers on the
> held-out slice.** A rule chosen and certified on one corpus is Shape 3 from `phase4_summary.md`
> §8, and this brief's own constraint — *every corpus you add gets a negative control* — has a
> mirror image here: a corpus of negatives certifying an extraction rule needs positives.

**3. Restore the false-reject half of Reviewer E's mandatory task, or give it an owner.**
Change E's commission to:

> Given the bindings and nothing else: **find a binding that accepts an answer to a different
> instance of the same template, and find a correct answer that the N1 rule leaves
> `UNRESOLVED`.** Draw both from gold×gold and archive×gold.

and add to the Track B exit gate: *"the decided rate for `numeric` did not fall against the
Phase 4 baseline, measured"*. Without this, F4's degenerate winner passes every gate in the
document. (If one reviewer with two tasks offends R6.1 — it does — that is another argument
for recommendation 1: Track B alone can afford two reviewers.)

**4. Amend Reviewer A's commission to the gate the brief itself repairs.** Replace *"Re-runs T4
across all 150 templates"* with:

> Re-runs T4 across all 150, **and independently checks the acceptance evidence for the three
> defect classes T4 cannot see** — the five `**Final Answer**` templates, the two malformed
> complex numbers, the unreachable `elif`. The implementer will hand you a check or a per-class
> evidence statement; **treat it as a claim under test, not as tooling.** This is the gate.

and reorder the Track A exit gate so the four-class item is first and `T4 150/150` second.
Assign the `ANSWER_MARKERS`-agreement item to A explicitly, as the one comparator file A may
read, or it has no owner (F2).

**5. State a threshold and a lever for D5.8, and restore the spec's ordering.** Add:

> A binding counts as bound when it produces **zero false accepts over N ≥ 50 instances**
> cross-paired within its template — 12 instances resolves nothing rarer than 25% and is not
> enough. Instances per template is yours to raise; generation is deterministic and cheap.
> Bind in this order and ship what is bound: **5 `classification`, 9 `symbolic`, 15
> `vector`/`array`, 32 `multipart`.** A named unbound template with a reason is a result; an
> unvalidated binding is a false accept waiting to happen.

(N is illustrative — pick it from the generation cost, but pick one.)

**6. Make the stopping rule operational, or drop it.** Replace the block quote with the
instrument that produced it:

> Track B's stopping rule is a **budget and an ablation**, not a judgement. Concretely:
> (a) record **findings per reviewer-hour** each round; when it stops falling, the surface is
> not converging and no further round will close it.
> (b) before adding any rule to fix a finding, **stub the mechanism it lives in and re-run every
> corpus.** If deleting it costs zero decided verdicts on real archive text, delete it instead
> of patching it. That ablation, not the quotation below, is what ended Phase 4.
> (c) Track B's version of "reduce scope" is **bind fewer templates, and name the rest** — not
> "add fewer rules". Say which templates you did not bind and why.

**7. State the boxes as numbers.** One line in §"How the two tracks relate": *"Track A: 6–10 h
implementation, 3–4 h review. Track B: 15–25 h, 8–10 h review — and Phase 4 budgeted 12–20 h +
8–10 h for four templates and spent eight review rounds, so treat Track B's figure as a trigger,
not a forecast. **When Track B passes 25 h with D5.8 incomplete, ship the bound templates, name
the unbound, and close.**"* Then state a number for each review brief's box, since the brief
requires the next session to do so and models it nowhere.

**8. Carry the two dropped inheritances.** Add R4-1 (high, symbolic comparator unvalidated —
directly in scope, 9 `symbolic` templates) to D5.11. And resolve D-057: `phase4_summary.md` §12
assigns Phase 5 the supporting-quantity remedy and the brief says "not yours to fix" — pick one
and record the disposition, in the document that names it.

**9. Cut, in this order.** The brief is 408 lines and its best passage is buried at line 116.

- **§"Track B — bindings…" lines 162–192 (~30 lines): cut the N1/N2/N3 prose to five.** It
  restates the spec's Phase 5 Track B table, which is required reading as item 1 of §"Read these
  first", and the confirmed-defect transcript adds no decision the session must make. Keep the
  reframing sentence (*"the comparator does not need more rules, it needs bindings"*) and the
  pointer.
- **T6 appears in five places** — the baseline table, two consecutive §Verification paragraphs,
  an exit-gate line, and Open items. Two of the three prose paragraphs say the same thing
  (baseline is stale, do not regenerate). Collapse to one.
- **Two warnings are stated twice, verbatim in substance,** and will be tuned out on the second
  reading rather than the first: *"Track A ships and Track B re-scopes"* (lines 66–68, 228–229)
  and *"never quote an exposure count as a defect count"* (143, 386). Also *"Phase 0 found two
  checks that were green while measuring nothing"* (127, 374). Keep each once, at the point of
  use.
- **Do not cut** §"What is actually here", the instance-dump paragraph, or the independence
  bullets. They are the parts that are doing work.

---

## Appendix — what I tried that failed to break the brief

Recorded because a review that only lists hits is indistinguishable from one that only looked
for them.

- **I expected `decimation_aliasing_analysis`'s "unreachable `elif`" fix to change computed
  answers**, which would have made Track A's "item pool did not move" gate unsatisfiable and the
  "formatting-only" framing wrong. It does not. The dead branch is `elif nf == 0` behind
  `if df == 1`, and reaching it changes the printed string `0*pi` to `0` — text, not value. The
  brief's classification of all four defect classes as output-contract hygiene holds.
- **I expected the two tracks to duplicate a measurement** (the R6.2 anti-pattern). They do not;
  the fences are clean. The defect is the opposite one — a gap, not an overlap (F2).
- **I expected D5.5 to be unusable for N1** on the strength of Reviewer E's F0 (`numeric` and
  `check` have zero archived negative instances). Cross-pairing manufactures negatives from real
  text without labels, which is precisely what F0 said was missing. That design is sound and it
  is the best idea in Track B.
- **I checked whether the brief's `ANSWER_MARKERS` claim was decorative.** It is not: the list is
  priority-ordered with `## Final Answer` above `**Answer:**`, and the comment block in
  `normalize.py` records a live defect caused by marker ordering. The brief is right that the two
  must agree — which is what makes F1 a coupling rather than a nit.

---

*Filed at the box. Not committed.*
