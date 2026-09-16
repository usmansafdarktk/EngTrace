# D4.1 — The comparator contract

**Version:** 1.0 · **Status:** normative for the milestone model · **Phase:** 4
**Date:** 2026-09-07
**Companion:** [`template_redesign_spec.md`](../template_redesign_spec.md) §4 ·
[`phase4_vocabulary.md`](phase4_vocabulary.md) (D4.2) ·
[`phase3_node_types.md`](phase3_node_types.md) §7 (which this extends)
**Reference implementation:** [`tests/comparators/`](../../../tests/comparators)
**Conformance:** `python -m tests.comparators.score` ·
`python -m tests.trace_schema.candidate_7 docs/re-implementation-sep/phase4_conformance/candidates.json`

---

## 0. What this is for

The paper is titled *Verifiable Process Supervision* and the verification is
currently a majority vote among three frontier LLMs, two of which come from
families the paper also evaluates. Five reviewers raised that independently
across two ARR rejections. Phases 1, C2, 2 and 3 made the gold traces
trustworthy — they reproduce their own answers, rest on constants verified
against NIST, regenerate from their recorded seeds, and emit a structured trace
where the trace has shape. **None of that scores a model.** The comparator is
what turns a gold trace into a verdict, and until it exists the tribunal is
still doing the verifying.

This document specifies it.

### 0.1 What it replaces, and what not to copy from it

`evaluation/engineering_parser.py` extracts **one float** from an answer and
compares it. It is the entire scoring surface today, and on the four Phase 4
templates it produces this:

```
python -m tests.comparators.score
```

| | archived traces for the four templates |
|---|---|
| scored correct by the deployed parser | **0 / 61** |
| actually correct (independent labels, `tests/comparators/ground_truth.py`) | **45 / 61** |

The parser is not merely wrong on these items; it is comparing numbers that do
not exist. Gold `a) Memoryless: **No** / b) Causal: **Yes**` parses to `1.0`.
Gold `z[n] = {-9, 6, *-5*, -4, 5, -3, 3}` parses to `-9.0` — the first element
of the sequence. Gold `u = -2x^2 + 6xy` parses to `-2.0` — the leading
coefficient. Every comparison it then makes is between two parse artifacts.

**Do not build this comparator to agree with that parser** (D-034, D-003). The
grammar, not the incumbent, is the reference: `parse_number` reads `4,921` as
four thousand nine hundred and twenty-one because that is what it means, not
because a test said so.

---

## 1. The three-valued verdict

**A comparator returns `MATCH`, `MISMATCH` or `UNRESOLVED`, and `UNRESOLVED`
can never contribute to a pass.**

This is the whole of Decision 3 and it is load-bearing (D-049). A two-valued
comparator must guess whenever a candidate does not commit, and every guess is
one of two errors:

- a **false accept** credits reasoning that did not happen — which is the AI
  Tribunal criticism in a new costume, and this entire effort exists to answer
  that criticism;
- a **false reject** penalises a correct solver and makes the benchmark measure
  phrasing rather than reasoning.

They are not equally bad *and neither is acceptable*, so the resolution is not
to pick a threshold between them. It is to stop guessing. `UNRESOLVED` is the
refusal, and it is modelled directly on D3.4 §8B.10's `unchecked` channel: an
unchecked relation may never contribute to a pass, so a node carrying one cannot
report a bare `PASS`.

**The bias, stated as the brief requires:** toward rejection over acceptance,
and toward `UNRESOLVED` over both. Quantified in §8.

**What stops `UNRESOLVED` becoming a free escape** is the definition of the two
gate numbers:

| | definition | what `UNRESOLVED` does to it |
|---|---|---|
| precision | `MATCH ∧ correct / MATCH` | nothing — it is not in the denominator |
| recall | `MATCH ∧ correct / correct` | **costs recall**, exactly as a `MISMATCH` on a correct answer does |
| decided | `(MATCH + MISMATCH) / all` | reported alongside, never traded against the other two |

A comparator that declines to decide everything scores 100% precision, 0%
recall, and 0% decided. The gate requires ≥95% on *both*.

---

## 2. Two rules that run through all six kinds

### 2.1 The gold side owns every parameter of the comparison

Precision, the label set, the origin convention, the symbol alphabet, the unit:
all read from **gold**, never from the candidate.

This is D3.4 §7.1's rule — *"always the GOLD node's `symbol_precision`, never the
candidate's"* — generalised. Reviewer D2 asked in round 6 for it to be stated
explicitly rather than left implicit in "the last place gold displays it",
because a careless implementer reads that as the candidate's rendering. Then a
candidate declaring coarser precisions widens the tolerance it is judged against.
**A candidate that can set a parameter of its own comparison is a candidate that
grades itself.** Exercised as `C9-candidate-widens-tolerance` in the §7
conformance corpus.

### 2.2 A value the comparator consumes is declared once, or recomputed

D3.3 §3.8, applied one layer up. The comparators hold **no reference answers**:
a gold canonical form is derived from the gold string by the same code path that
derives the candidate's. There is no `EXPECTED = {...}` dict anywhere in
`tests/comparators/`, and that is not stylistic — it is D-034, which is on the
register because a green 215-check suite once certified three constants against
a dict in the test file and nothing else.

The vocabulary tables are the one place counts *are* written down, and
`derive_vocabulary.py` recomputes every one from the archive and fails on
disagreement (§D4.2).

---

## 3. `multipart` is not a seventh kind — Decision 1

**`multipart` is a property of the *answer*, not a `kind`.** An answer carries an
ordered list of `parts`; each part carries exactly one of the six kinds. The six
stand and the answer schema grows. Recorded as **D-047**.

The evidence is inside this phase. `system_properties_memory_causality` is typed
`classification` by the inventory, described by the spec as a *"canonical label
tuple"*, and is in fact **a two-part answer whose parts are both categorical**.
Under a seventh-kind reading it would have to be typed `multipart` — and its
comparator would then have no way to say *"each part is categorical, normalise
each against the categorical vocabulary"*. The six kinds would stop composing at
exactly the point they are most useful. Under the composition reading it is
`parts=[categorical, categorical]` and reuses the categorical vocabulary
unchanged.

Sampling the 32 `multipart` templates shows the same structure, and shows
something the word hides — **two different structures wear it**:

| `mode` | meaning | example |
|---|---|---|
| `all` | every part required, every part must match | `volumetric_flow_rate`: flow rate **and** average velocity |
| `any` | the parts are alternative renderings of one quantity | `undamped_natural_frequency_torsional`: "51.184 rad/s **or** 8.146 Hz" |

Conflating them is a scoring error in both directions: read `any` as `all` and a
model answering in rad/s alone is marked wrong; read `all` as `any` and a model
that gets the flow rate right and the velocity wrong takes full credit.

**Partial credit is reported, never folded into the verdict.** A partially
correct `all` answer is not a correct answer; a comparator that says otherwise
weakens 32 items silently, which is a larger P6 surface than any template edit
this project has made. The per-part outcomes are always in `observations`, so
partial correctness is legible in the results without being scored as a pass.

**Precedence for the composite verdict:**

- `all`: any `UNRESOLVED` part makes the whole `UNRESOLVED`. A part we could not
  decide might have been the one that was wrong; resolving it as `MISMATCH`
  reports a model error nobody observed.
- `any`: a single `MATCH` is a `MATCH` even beside an `UNRESOLVED` sibling,
  because one alternative rendering has already settled it.

**Leverage.** `multipart` is 32 of 150 and **24 of the 51 class-C templates** —
so this decision governs nearly half of class C. Verified: the inventory does
carry an `instrumentation_class` column (the brief says it does not), and the
51 reproduces.

---

## 4. The six kinds

### 4.1 `numeric`

Half a unit in the last place **gold** displays, boundary **inclusive**:
`|c − g| ≤ 0.5 × 10^(−p)`, `p` from gold's declared precision or, absent one,
from gold's own rendering — never the candidate's.

Widened by **§7.4** when the question prescribes a method: the answer tolerance
is the looser of the display tolerance and the method's own stated tolerance. An
item whose question prescribes *how* to solve it has two defensible answers, and
marking a direct solver wrong for one display unit is the same error as marking
a four-update solver wrong for taking four updates.

**Thousands separators are parsed, not stripped by accident.** `4,921` is 4921.
D-003's defect is that the deployed parser applies its comma strip to a
different string from the one its regex matched.

**Units are checked only when declared.** D3.4 §10 records "no dimensional
checking" as a non-goal because a dimensional comparator needs units on every
symbol, which is a milestone-model decision. The residue is a real false accept:
`7.65 mL` against a gold of `7.65 L` matches on the number. `compare_numeric`
takes an optional declared `unit`, closing it **per item by declaration** and
leaving it open, visibly, for every item that does not declare one. Inferring
the unit from gold's string instead would reject `7.65 litres`, which is
correct. That trade is why it is opt-in, and the open residue is **D-052**, kept
in the adversarial corpus as `num-16b` so the gate carries it rather than hides
it.

### 4.2 `categorical`

A single label from a **closed, declared** set. `labels` maps each canonical
label to its accepted surfaces; it is gold-side problem data and is the only
vocabulary in play.

**Surfaces match longest-first.** `nonlinear` beats the `linear` inside it.
Matching the short one scores all twelve archived "not linear" answers as
"linear" — which the first draft of the implementation did.

**Negation is scoped, not global.** A negator governs a label within ~40
characters, and a clause boundary (`.`, `;`, `but`, `however`, `whereas`,
`while`) ends its scope. `neither A nor B` negates both labels and is handled
separately because its scope spans two of them.

**A fused negative is a declaration, not an inference.** `nonlinear` and
`non-causal` are declared as surfaces of the negative label. Inferring a
negation *on top of* the declaration double-negates and returns the positive.

**Negation of a two-valued set is its other member; of a larger set it is
nothing.** "not laminar" leaves transitional and turbulent open, so it is
`UNRESOLVED` rather than silently the complement.

**A hedge is `UNRESOLVED`, and hedges are scoped like negation.** A
non-committal answer has not answered — this is the specific case spec §4.4
names for Reviewer B. Two words are *deliberately absent* from the hedge list:
`it depends` and `either`, because both are domain vocabulary here (`a) Not
memoryless (it depends on the past input x[n-2])` is a committed, correct
answer, and `fails at least one of the tests (additivity or homogeneity)` is
gold's own phrasing). With them in the list the comparator returned `UNRESOLVED`
on a correct archived answer. **A hedge vocabulary built by asking "what sounds
non-committal?" rather than by testing it against the corpus produces exactly
that, and produces it as a false reject — invisible in an accuracy number.**

**A label-tuple** (`categorical[tuple]`) is a composite of categorical parts,
per §3, and is the worked proof of D-047. Two encodings are in play and **gold
uses the rarer one**: gold writes `a) Memoryless: **No**` (slot named, answered
yes/no) while models overwhelmingly write `a) Not memoryless` (property asserted
or denied). A comparator reading only gold's encoding scores 21 of 24 archived
traces wrong. Slots are located by their own vocabulary when present and **by
position** otherwise — position is not a convenience fallback: `a) Yes, b) Yes`
names neither property and is a complete answer.

### 4.3 `sequence`

Ordered values **and** the index of `n = 0`. Presentation is not part of the
answer; the origin is.

**The origin is answer-bearing — and the archive does NOT prove it.** Version
1.0 of this section claimed signal trace 8 emitted "exactly gold's multiset of
values in exactly gold's order" and called it "1 of 16 archived traces". Both
halves were false, and Reviewer E filed it as F1: trace 8's answer is a
**rotation**, which a values-only comparator also rejects, and
`ground_truth.py`'s own label for that trace said so three files away.
Measured: **0 of 16 archived traces reach the origin-discriminating branch.**

What survives is an argument from the **item**, not from observed model error.
The question states `x[n]` with its origin marked and the transformation maps
indices, so a value list at the wrong indices is a different signal — and traces
1, 6, 10 and 13 *are* labelled wrong for exactly that, they simply have wrong
values too, so the values check catches them first. `require_origin` is
therefore a **route-3** rule in D4.2 §2.4's sense: specified and adversarially
exercised, **not validated against real output**. It is on the residual-risk
register as R4-11.

**Five packings are parsed and two archived traces use more than one at once**,
so the parser produces a single index→value map and merges rather than
dispatching on a form:

| packing | example | n |
|---|---|---:|
| braces, asterisk origin | `{2, 3, *-8*, 6}` | 1 |
| braces, no origin stated | `{-7, 7, 5, -1}` | 9 |
| braces + explicit index list | `{…} for n = {-6, -5, …}` | 2 |
| per-element assignments | `z[-1] = -8, z[-2] = -3` | 1 |
| braces **and** assignments together | `{6,0,0,0,0} for n={0..4} and z[-1] = -8, …` | 1 |

Leading and trailing zeros are padding and are trimmed; **interior zeros are
structural and are never trimmed**.

**Gold itself omits the origin on 16.4% of instances** (measured, 4,000 seeds):
the template marks `n = 0` only when the result's support contains it. When gold
is silent the answer is the value list alone and is compared as such.

**A candidate that pins an origin where gold cannot is `UNRESOLVED`, not a
match.** Version 1.0 said such a candidate was "checked for consistency rather
than punished for saying more", and Reviewer E's F3 showed what that accepted:
three candidate forms taken from real archived traces, all pinning *wrong*
origins, all scored `MATCH`. On these instances the item cannot decide, so the
comparator does not either. The cost is visible and is reported rather than
absorbed: archived `signal_operations` recall falls from 100% to 50% and its
decided rate to 75%, because **one of its two matches was a match on an item
that cannot discriminate** — which is Reviewer E's F6, settled by measurement.

Recorded as **D-050**, because a printed gold answer from which the origin
cannot be recovered is a property of the item worth knowing, and E's §5.2 notes
the better fix is in the template — print the support explicitly, `y[n] = {…}
for n = 1…4`, even when `n = 0` lies outside it. That changes emitted text on
16.4% of instances, so it is a Phase 5/6 item-design decision, not a Phase 4
edit.

**When gold pins the origin and the candidate does not**, the comparator first
asks whether *any* placement would match. If none does, the values are wrong
whatever the origin and the case is a decidable `MISMATCH`; only if some
placement would match is it `UNRESOLVED`. This converts 3 of 16 archived traces
from undecided to a defensible verdict and converts none the other way.

### 4.4 `symbolic`

Two stages.

**Stage 1 — the polynomial fragment, decided exactly and locally.** Every answer
this kind sees for `incompressible_continuity` is a bivariate polynomial with
rational coefficients and degree ≤ 2, and equality there is decidable by
expanding to a coefficient map: no CAS, no timeout, no dependency, and **no
`UNRESOLVED` outcome**. That matters more than the convenience — this is the
template with **six** archived traces, and a comparator whose verdict depended on
whether an optional package happened to be installed would make the thinnest
evidence in the phase thinner still.

**Stage 2 — a CAS outside the fragment.** The corpus's other eight `symbolic`
templates carry sinc, exp and Q-functions, so the fragment is not enough in
general. `sympy` is used there and **is a new dependency** (it was not installed
and is not in `requirements.txt`); recorded as **D-053**.

**The arbitrary-function policy.** The item asks for the *simplest* expression
and gold sets the arbitrary function to zero. A candidate carrying it explicitly
— `3x²/2 − 6xy + C(y)` — has answered the question and additionally shown it
understands why the function is there. **Accepted**, with the retained term
recorded as an observation. This is not lenience: the general solution and the
simplest particular solution differ by exactly that term, and rejecting it marks
a *more* complete answer wrong. A candidate whose answer is *only* the arbitrary
term (`u(x, y) = C`) is a `MISMATCH` — it dropped the solution rather than
annotating it. Both cases are in the archive.

**Failure is `UNRESOLVED`, never `MATCH`.** A CAS timeout, a parse failure, or a
symbol outside the declared alphabet all decline. Nothing about a comparator
failing to parse an answer is evidence the answer is right, and a `sympy`
exception must not become a silent pass — §8B.10's rule for this kind.

### 4.5 `narrative`

**Always `UNRESOLVED`. That is the kind's entire specification** (D-051).

`narrative` names a milestone whose content is prose no comparator here can
decide — a justification, an interpretation, a physical explanation. Specifying
it as "compare with an LLM" would put the AI Tribunal back on the critical path
for exactly the milestones where it is least accountable, which is the thing this
phase exists to remove.

So the kind is **undecidable by construction**: declaring a milestone
`narrative` *removes it from the automatic score* and routes it to whatever
human or model process the benchmark chooses, with that routing visible in the
results rather than hidden inside an accuracy number.

**Its gate is a census, not precision.** The fraction of milestones declared
`narrative` is a reported quantity, and a rising fraction is a benchmark quietly
returning to LLM judging. Conformance: 20 adversarial cases, including one where
the candidate restates gold *verbatim*, and all 20 must return `UNRESOLVED`. A
single `MATCH` would mean the kind had started deciding prose.

### 4.6 `check`

A verdict-bearing assertion — "verify the beam is safe", "confirm the flow is
laminar". Its answer has two parts and **both are scored**: the boolean, and,
when gold states one, the quantity the verdict rests on, at §4.1 tolerance.

**Scoring only the boolean makes the kind a coin flip.** A candidate that
computes a deflection of 40 mm against a 25 mm limit and concludes "not
acceptable" agrees with a gold that computed 30 mm. Crediting it credits
reasoning that did not happen, on a two-valued answer where guessing is free.
This is where Decision 3's asymmetry bites hardest, and the supporting quantity
is what removes it. Adversarial case `chk-11`.

A candidate that states the verdict but not gold's supporting quantity is
`UNRESOLVED`, not `MATCH`.

---

## 5. What this specification does **not** settle

Stated as non-goals so an implementer meets the edge here rather than in the
middle of a scoring run.

- **No dimensional checking by default.** §4.1. Inherited from D3.4 §10 and only
  partly closed. `num-16b` in the adversarial corpus is the live residue and it
  is the one false accept the gate carries (D-052).
- **`narrative` is not scored at all**, by design. What routes it, and how that
  routing is reported, is a milestone-model decision (D-051).
- **The comparator scores an answer, not a derivation.** Whether the *reasoning*
  reaching a correct label was valid is D3.4 §7's process credit for the two
  node types that have one, and is out of scope for the other four kinds. A
  model that reaches "linear" by invalid reasoning takes answer credit here;
  spec §4.4 puts that explicitly out of Reviewer B's scope too.
- **The vocabulary is not validated against a representative sample for
  `symbolic`.** Six traces. §D4.2 and the residual-risk register.
- **`check`'s truth/false surface lists are not derived from the archive.** No
  Phase 4 template is a `check`, so those 21 surfaces are the one part of this
  document written from the domain rather than from observation. Flagged rather
  than smoothed over: the first `check` template to be fitted should re-derive
  them.

---

## 6. Amendments to `phase3_node_types.md` §7

Building the §7 conformance corpus (D4.6) found four defects in §7 itself. They
are defects in the *specification*, which had never been exercised — the gap
both Phase 3 schema reviewers named as the largest.

| # | Finding | Disposition |
|---|---|---|
| **F7-1** | §7.2's headline disposition — *"a correct answer reached in four iterations instead of three is correct"* — **has no instances**. §4.3.3 ties `converged` to the node's own tolerance, §8A.4 forbids `converged: true` before the last element, §4.3.1 recomputes each update, §4.6 recomputes each frame, and the preamble is fixed by the question. Together these make the conforming node for a question **unique**. Verified both ways. D-038's measured 1..5 spread is **between questions**, not between solvers of one, and §7.2 conflates them. | Amend §7.2 |
| **F7-2** | §7 never says whether §6 step 8 binds a candidate. It does in every verifier built, and it should — but the consequence is that answer credit and process credit are **not independent** for `iteration`, so "wrong answer, clean process" is not constructible. §7 is written as though they vary freely. | Amend §7 |
| **F7-3** | §7.4's direct solver cannot be represented as an `iteration` node at all, because a direct solve has no iteration. Its tolerance is correct and testable only at the function level. | Note in §7.4 |
| **F7-4** | §7.3 row 3's set-versus-ordered distinction is **vacuous** for a §6-conforming node: §5.4.6 already requires `committed` to equal the ordered chosen values. | Amend §7.3 |

**F7-1 is the one that matters, and it is the six review rounds' bill arriving.**
Each round closed a way for a candidate to lie. Together they also closed every
way for a candidate to be *differently right*. Nobody recorded that trade
because §7 was never run.

D4.9 finds the same thing from the other direction and on most of the corpus:
**141 of 150 templates emit a trace whose length is a constant**, so
`cardinality` has one possible value and both §7.2 and §7.3 dispositions are
vacuous there too.

---

## 7. Serving class C

The spec says these comparators "also serve the 51 class-C templates" and that
Phase 4 is where most of class C's design cost is paid. Verified: **51
reproduces** from `instrumentation_class`, and the composition is what makes the
claim concrete.

| `answer_type` | all 150 | class C | kind |
|---|---:|---:|---|
| scalar | 89 | 6 | `numeric` |
| multipart | 32 | **24** | composite — §3 |
| symbolic | 9 | 8 | `symbolic` |
| vector | 8 | 8 | `sequence` |
| array | 7 | 3 | `sequence` |
| classification | 5 | 2 | `categorical` |

**Nearly half of class C is `multipart`**, so Decision 1 is not a bookkeeping
question — it is the single decision that governs most of the class the phase is
supposed to be paying for.

**A caution for whoever fits class C.** Two of the six kinds are specified
against evidence and four are not, in different degrees:

| kind | Phase 4 evidence | status |
|---|---|---|
| `categorical` | 39 archived traces | derived |
| `sequence` | 16 | derived |
| `symbolic` | **6** | specified, adversarially exercised, **not validated** |
| `numeric` | 0 Phase 4 traces (rules from the wider archive and the grammar) | specified |
| `check` | 0 | specified from the domain |
| `narrative` | 0 | undecidable by construction, so no evidence is needed |

---

## 8. Measured

```
python -m tests.comparators.score
```

**Against the 61 real archived traces**, scored against independent labels in
`tests/comparators/ground_truth.py`, which does not import the comparator:

| template | n | correct | accepted | precision | recall | decided |
|---|---:|---:|---:|---:|---:|---:|
| `signal_operations` | 16 | 2 | 2 | 100% | 100% | 81.2% |
| `system_properties_memory_causality` | 24 | 24 | 24 | 100% | 100% | 100% |
| `system_property_linearity` | 15 | 15 | 15 | 100% | 100% | 100% |
| `incompressible_continuity` | 6 | 4 | 4 | 100% | 100% | 100% |
| **total** | **61** | **45** | **45** | **100%** | **100%** | **95.1%** |

**This number is weak evidence and is reported as such.** The archive is where
the vocabulary came from, so scoring 100% on it measures memorisation. It is the
same shape as Phase 3's gold corpus passing 80/80 against a verifier that
hardcoded one template's symbol names — gold traces do not lie, and neither does
a training set.

> **And a 100% here says nothing about whether the item measures reasoning.**
> Two of these four templates are **fully predictable from the question surface**
> — `system_property_linearity` at 100% held-out against a 50.02% blind-guess
> floor, `system_properties_memory_causality` at 100% against 34.02%, with the
> form→label map a total function over 5 and 7 distinct right-hand sides
> (Reviewer B, D-057). The comparator scores those answers perfectly and a
> solver who has never heard of homogeneity scores them perfectly too. The two
> numbers now live per template as `blind_guess_floor` and
> `surface_model_heldout` in [`template_inventory.csv`](../audit/template_inventory.csv);
> flag at lift ≥ 40 points **or** held-out ≥ 95% regardless of lift, since a
> high floor masks a total lookup. **Three of the 150 templates are over that
> threshold.**

**Against D4.4**, 157 hand-written near misses, ≥20 per comparator:

| kind | n | correct | accepted | precision | recall | decided |
|---|---:|---:|---:|---:|---:|---:|
| `categorical` | 24 | 13 | 13 | 100% | 100% | 75.0% |
| `categorical[tuple]` | 22 | 13 | 13 | 100% | 100% | 86.4% |
| `sequence` | 24 | 11 | 11 | 100% | 100% | 91.7% |
| `symbolic` | 22 | 11 | 11 | 100% | 100% | 95.5% |
| `numeric` | 24 | 13 | 14 | 92.9% | 100% | 100% |
| `check` | 21 | 9 | 9 | 100% | 100% | 61.9% |
| **total** | **137** | **70** | **71** | **98.6%** | **100%** | **85.4%** |
| `narrative` | 20 | — | — | — | — | 20/20 `UNRESOLVED` as required |

**Gate: ≥95% precision and ≥95% recall. Met.**

**The single remaining false accept is deliberate and named:** `num-16b`,
`7.65 mL` against a `7.65 L` gold with no unit declared. It is kept in the
corpus so the gate carries D3.4 §10's open non-goal rather than concealing it.

**The adversarial set earned its place**: it found six things the archive could
not. Three were comparator defects (a bare `No` in "No memory" read as the slot
*value*; `nonlinearity` matching no label at all; `x^{2}` surviving LaTeX
stripping as a brace group). Two were errors in the labels I wrote. One is the
residual above. All six are in `phase4_summary.md` §8.

---

## 9. Conformance

A conforming implementation must:

1. Return one of the three outcomes, and never let `UNRESOLVED` reach a pass.
2. Read every comparison parameter from gold (§2.1).
3. Carry no reference answers (§2.2).
4. Score `narrative` as `UNRESOLVED` unconditionally, including against a
   verbatim restatement of gold.
5. Score `check` on the quantity as well as the boolean whenever gold states one.
6. Reach ≥95% precision **and** ≥95% recall on `tests/comparators/adversarial.json`.
7. Reproduce the 16 dispositions of
   `docs/re-implementation-sep/phase4_conformance/candidates.json` and the 9
   tolerance assertions, with all 15 §7 dispositions exercised.

Items 6 and 7 are runnable; items 1–5 are behavioural obligations in the sense of
D3.4 §8B, and item 4 is the one that is testable only as a behaviour.
