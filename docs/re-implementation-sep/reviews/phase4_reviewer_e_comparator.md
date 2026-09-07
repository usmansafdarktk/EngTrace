# Phase 4 — Reviewer E: comparator adversary

Frozen ref: `f6bb653745785a2be164ce2756cd748f32fb4a51` (branch `redesign/phase4-comparators`).
All commands below are run from the repository root with `PYTHONIOENCODING=utf-8` set.

---

## 1. Verdict

**BLOCKED** — the phase's exit gate ("Reviewer E found no false accept on real archived traces") does not hold: four false accepts are reproduced below, one of them built entirely from real archived answer forms and a real archived gold, and the archive is structurally incapable of testing four of the six kinds.

---

## 2. Independent re-derivation

| number | as reported | as re-derived | agreement |
|---|---|---|---|
| archive precision / recall / decided | 100% / 100% / 95.1% over 61 | reproduced exactly by `python -m tests.comparators.score` | **agrees**, but see F0 — the denominator is the problem, not the arithmetic |
| D4.4 precision / recall | 98.6% / 100% over 137, 1 false accept (`num-16b`) | reproduced exactly | **agrees** |
| narrative conformance | 20/20 UNRESOLVED | reproduced | **agrees** |
| ground-truth census | 45 correct / 16 wrong of 61 | `python -m tests.comparators.ground_truth` → 45 / 16 | **agrees** as arithmetic |
| `signal_operations` label self-check | clean | reproduced clean | **agrees** |
| gold states the origin on 83.6% of instances (silent on 16.4%) | 3344 / 4000 | reproduced 3344 / 4000 | **agrees** |
| `normalize.SUPPORT`, 27 entries | machine-checked against the archive | `python -m tests.comparators.derive_vocabulary` → 27/27 agree | **agrees** |
| rules with zero Phase 4 support | 5 | reproduced (unicode-operator, latex-text, latex-exponent-brace, hedge, thousands-separator) | **agrees** |
| archived signal golds that omit the origin | not reported anywhere | **7 of 16 (43.8%)** | **new number**, see F3 |
| archived traces exercising the "values agree, origin differs" branch | claimed as trace 8, "1 of the 16 archived traces" | **0 of 16** | **DIVERGES**, see F1 |
| wrong answers in the archive, per template | not reported anywhere | signal 14, continuity 2, memory-causality **0**, linearity **0** | **new number**, and it is F0 |

Additional instrument built (not in the brief's toolkit): **cross-pairing**, in which every archived model answer labelled correct is scored against the *real gold of every other archived item in the same template*. Truth for each pair is settled by an independent gold parser that does not import `tests.comparators`. Result: **782 real-text pairs, 0 false accepts, 0 false rejects, 0 UNRESOLVED** (210 linearity, 552 memory-causality, 20 continuity). Reported in §4 — it is the strongest positive evidence in this review, and it is why the §3 failures matter: every one of them lies outside the shapes the archive contains.

---

## 3. Findings

### F0 — CONFIRMED (blocking, structural): the archive gate is vacuous for four of the six kinds

"100% precision on 61 real archived traces" is computed over an archive that contains **no incorrect answer at all** for `categorical` (0 of 15 linearity traces wrong) or `categorical[tuple]` (0 of 24 memory-causality traces wrong), and **no trace of any kind** for `numeric`, `check` or `narrative`. Precision measured over an empty set of negatives cannot fail. The entire archive-side false-accept gate rests on **13 decided sequence rejections plus 2 symbolic rejections — 15 negative instances**, all in two kinds.

Reproduction:

```bash
python -c "import sys; sys.path.insert(0,'.'); \
from tests.comparators.ground_truth import LABELS; \
[print(k, 'wrong answers in archive =', sum(1 for ok,_ in v.values() if not ok)) for k,v in LABELS.items()]"
```

→ `signal_operations 14` / `system_properties_memory_causality 0` / `system_property_linearity 0` / `incompressible_continuity 2`

**Impact.** The gate wording must not be read as evidence about `categorical`, `categorical[tuple]`, `numeric`, `check` or `narrative`. As written it is true and uninformative for 39 of the 61 traces. This is not a dispute with any label; it is a dispute with what the number licenses.

---

### F1 — CONFIRMED (blocking): the archive does not support the sequence comparator's central claim; "1 of the 16 archived traces" is wrong

`kinds.py:682` and `docs/re-implementation-sep/phase4_comparators.md:251` both assert: *"Signal trace 8 emits exactly gold's multiset of values in exactly gold's order and is wrong, because it indexes the reversal from n = 0 … A values-only comparator scores it a match. That is the false accept this kind exists to prevent, and it is not hypothetical — it is 1 of the 16 archived traces."*

Both halves are false.

* Trace 8's answer is `{3, -9, 6, -5, -4, 5, -3}`; gold's value list is `{-9, 6, -5, -4, 5, -3, 3}`. Same multiset, **rotated by one** — not "exactly gold's order". `ground_truth.py`'s own label for trace 8 says "is a rotation of the truth", contradicting the docstring three files away.
* A values-only comparator therefore **also rejects it**. The verdict actually taken is `"no placement of n=0 makes the candidate's values equal gold's"` — the values-only placement search, reached before any origin comparison.
* **No archived trace reaches the origin-discriminating branch at all.**

Reproduction:

```bash
python - <<'PY'
import sys, json, glob; sys.path.insert(0, '.')
from collections import defaultdict
from tests.comparators.answer import compare_template
rows = []
for p in sorted(glob.glob('error_analysis_annotation/samples/*.jsonl')):
    for line in open(p, encoding='utf-8'): rows.append(json.loads(line))
G = defaultdict(list)
for r in rows: G[r['question_id'].rsplit('__', 1)[0]].append(r)
n = 0
for i, r in enumerate(G['signal_operations']):
    v = compare_template('template_signal_operations', r['gold_answer'], r['model_reasoning'])
    if 'origin differs' in (v.reason or ''): n += 1
    print(f"[{i:2}] {v.outcome:10} {v.reason[:80]}")
print("traces taking the 'values agree but the origin differs' branch:", n)
PY
```

→ `traces taking the 'values agree but the origin differs' branch: 0`

**Impact.** At the one place the spec claims the archive *proves* a design decision, the archive proves nothing. `require_origin` has zero archive support for its discriminating power; it is a route-3 rule exactly like the symbolic ones, and D4.2 §5 should say so. The prose must be corrected before merge — it is the load-bearing empirical claim of the `sequence` kind, and it is the sentence a later reviewer will cite.

---

### F2 — CONFIRMED (blocking): **false accept.** `_label_hit` prioritises by surface *length*, overriding position; the docstring claims the opposite

`_resolve_label`'s docstring: *"a set that declares both polarities explicitly is resolved by position, latest wins."* It is not. `_label_hit` sorts surfaces longest-first and `break`s on the first surface that hits *anywhere*, so a single mention of `nonlinear` anywhere in the answer span outranks the committed answer `linear`, regardless of where each sits.

Reproduction:

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:**\nThe system is **not linear**."
c = "## Final Answer\n**Answer:** The system is linear. Unlike a nonlinear system, it satisfies both tests."
v = compare_kind('categorical', g, c)
print(v.outcome, v.gold_canonical, v.cand_canonical)
PY
```

→ `MATCH not linear not linear`

A candidate that commits to **linear** — the wrong answer — is credited as `not linear`. The defect also fires in reverse: `"Both tests pass, so the system is linear (it is not nonlinear)."` resolves to `linear`, so the sign of the error depends on incidental phrasing rather than on the answer.

**Realism.** Contrastive mention of the opposite class is ordinary answer prose for this item family, and three of the fifteen archived linearity golds are `linear`, so both polarities are live. The archive happens not to contain such a span — because (F0) it contains no wrong linearity answer at all — which is precisely why this survived to review.

**Impact.** A false accept and a false reject in one defect, in the kind with the largest archived trace count and zero archived negatives. Fix: rank hits by position first, use length only to break ties at the same position; or scope resolution to the final clause of the span.

---

### F3 — CONFIRMED (blocking): **false accept.** When gold omits the origin, `require_origin=True` is silently inert and a candidate's *wrong* origin is accepted

`compare_sequence` documents the origin as answer-bearing, then carries a branch `g.origin is None and c.origin is not None` that compares **trimmed values only** and records "candidate stated an origin gold does not" as an *observation* on a `MATCH`. Gold omits the origin on **7 of the 16 archived signal traces (43.8%)** and on **16.4% of the 4,000-seed item pool**, so the flag is inert on a large minority of instances.

Reproduction — gold below is the real archived gold of signal traces 0, 7 and 14 (truth: n = 1…4):

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:**\ny[n] = {-1, 5, -4, 3}"
for c in ["**Answer:** y[n] = {*-1*, 5, -4, 3}",
          "**Answer:** y[n] = {-1, 5, -4, 3} for n = -3, -2, -1, 0",
          "**Answer:** y[0]=-1, y[1]=5, y[2]=-4, y[3]=3"]:
    print(compare_kind('sequence', g, c).outcome, '|', c)
PY
```

→ `MATCH` three times: the first pins n = 0 at the wrong element, the second places the support four indices off, the third is the off-by-one placement.

The candidate forms are not invented. The asterisk-origin form is archived signal trace 0; the `for n = {...}` index-list form is archived traces 1 and 3; per-element assignment is archived trace 1. The error itself — a value list placed at the wrong indices — is the labelled error of archived traces 1, 6, 10 and 13, four of sixteen. All that is combined here is a real archived error mode with a real archived gold from the same template.

**Impact.** The comparator's stated purpose ("the origin is part of the answer") is unenforced on 16.4% of the item pool, and what it accepts there are wrong signals, not merely under-specified ones. The observation string is not a defence: it is emitted on a `MATCH` and is invisible in precision. This is an *item* defect as much as a comparator one (D-050 records that gold cannot print an origin when the support excludes n = 0) — but the comparator converts an undecidable item into a silent accept rather than an `UNRESOLVED`, which is the one outcome the evidence does not support.

---

### F4 — CONFIRMED (blocking): **false accept.** `answer_span` peels past the committed answer into a trailing caveat

`answer_span` peels the *last* occurrence of the highest-priority marker, repeatedly, up to four times. A caveat sentence containing `Answer:` after the committed answer therefore becomes the answer.

Reproduction:

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
from tests.comparators.normalize import answer_span
g = "**Answer:**\nThe system is **not linear**."
c = "## Final Answer\n**Answer:** The system is linear.\n\nNote: Answer: not linear if the offset were nonzero."
print(repr(answer_span(c)[0]))
print(compare_kind('categorical', g, c).outcome)
PY
```

→ `'not linear if the offset were nonzero.'` then `MATCH`

The `Answer:?` marker (priority 5) is bare enough that any "Answer: …" or "The answer is …" in a post-hoc remark outranks the model's actual commitment. `marker-the-answer-is` has support 4 in 61 and 202 in 2,200, so the low-priority markers do fire on real output; peel-to-last is what makes them dangerous.

**Impact.** A false accept whose sign is arbitrary. The docstring justifies last-wins by "models restate the answer and the final statement is the committed one" — true for restatements, false for conditionals and caveats, and the rule cannot tell them apart.

---

### F5 — PLAUSIBLE: the `neither … nor` window (60 chars) inverts both slots when a clause intervenes

`NEITHER_NOR_RE` allows at most 60 characters between `neither` and `nor`. Both archived uses have a **12-character** gap, so the archive never probes the bound by a factor of five. A longer, ordinary interposed clause escapes it and *both* labels then read positive — the exact inversion the constant's own comment calls "the single most dangerous construction in the archive".

Reproduction:

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:**\na) Memoryless: **Yes**\nb) Causal: **Yes**"
c = ("## Final Answer\n**Answer:** The system is neither memoryless, because the output "
     "depends on past inputs as shown in step one, nor causal.")
v = compare_kind('categorical[tuple]', g, c)
print(v.outcome, v.gold_canonical, v.cand_canonical)
PY
```

→ `MATCH ['Yes', 'Yes'] ['Yes', 'Yes']`

A candidate asserting (No, No) is read as (Yes, Yes) and, against a (Yes, Yes) gold, accepted. Against a (No, No) gold the same defect is a false reject. Rated PLAUSIBLE rather than CONFIRMED only because the interposed causal clause is my construction; it is stylistically unremarkable, and the two archived instances come from different model families (gemma, meta-llama), so the form is not one model's tic.

---

### F6 — PLAUSIBLE: ground-truth label `signal_operations[15]` records as correct an answer the trace's own reasoning contradicts

Label: `15: (True, "{-10,-10,-7,-4,6,9,-8} equals gold exactly; gold states no origin either")`. The reading of the *printed* answer is right. But the trace reaches it by asserting `x[n-(-2)] = {x[2], x[3], … x[8]}` — a **left** shift in 1-based positional indices, ignoring the `*9*` origin the question states — and then emits the input unchanged. The printed answer coincides with gold only because gold's support excludes n = 0 and so gold prints no origin either (F3).

I am not asking for the label to be flipped: `correct` is defined in `ground_truth.py` as a property of the model's *final answer*, and as printed that answer is indistinguishable from gold. The finding is that this trace is the archive's clearest instance of F3's degeneracy — one of only two archived `MATCH`es in `signal_operations` is a match on an item that cannot discriminate — and the justification text should say so, because as written it reads as positive evidence for the comparator when it is evidence about the item.

I attempted to dispute the other 44 labels and could not (§4).

---

## 4. Falsification attempts that failed

These are the attacks that did **not** break the comparators. They are what makes F0–F6 a bounded list rather than an impression.

1. **Cross-pairing, 782 real-text pairs.** Every archived model answer labelled correct, scored against the real gold of every *other* archived item in the same template. Truth for each pair is settled by an independent gold parser that does not import `tests.comparators` (regex over gold's `**Answer:**` block for the two categorical templates; `sympy` over gold's `u = …` line for continuity). Pairs whose golds differ must not match; pairs whose golds agree must match.
   → **210 linearity + 552 memory-causality + 20 continuity = 782 pairs; 0 false accepts, 0 false rejects, 0 UNRESOLVED.**
   The categorical, label-tuple and symbolic comparators are genuinely robust across every real answer form the archive contains, at every real gold in the template — a far stronger statement than the 61-trace diagonal, and it holds.

   ```bash
   # sketch; the full script is ~40 lines and lives in the review's working notes
   python - <<'PY'
   import sys, json, glob, re; sys.path.insert(0, '.')
   from collections import defaultdict
   from tests.comparators.answer import compare_template
   from tests.comparators.ground_truth import LABELS
   rows = []
   for p in sorted(glob.glob('error_analysis_annotation/samples/*.jsonl')):
       for line in open(p, encoding='utf-8'): rows.append(json.loads(line))
   G = defaultdict(list)
   for r in rows: G[r['question_id'].rsplit('__', 1)[0]].append(r)
   def ans(g): i = g.rfind('**Answer:**'); return g[i:] if i >= 0 else g
   def lin(r):
       b = ans(r['gold_answer']).lower().replace('*', '')
       return 'not linear' if ('not linear' in b or 'nonlinear' in b) else 'linear'
   rs = G['system_property_linearity']; truths = [lin(r) for r in rs]
   ok = [i for i, (c, _) in LABELS['system_property_linearity'].items() if c]
   bad = [(i, j) for i in ok for j in range(len(rs)) if i != j
          and compare_template('template_system_property_linearity',
                               rs[j]['gold_answer'], rs[i]['model_reasoning']).is_match
          != (truths[i] == truths[j])]
   print('linearity cross-pair disagreements:', len(bad))
   PY
   ```

2. **Negation-window escape (40 chars).** `"The system does not, on the evidence of both tests above, qualify as linear."` places the negator 55 characters upstream. Correctly `MISMATCH`: the label resolves as `linear` against a `not linear` gold. I could not turn the window overrun into an *accept* on this binding, because with a two-member set an escaped negation yields the wrong label rather than the right one. The window is a false-reject risk here, not a false-accept one.
3. **Clause-boundary cut in `_negated`.** `"It is not additive; linear."` and `"The system is not homogeneous but linear."` both resolve to `linear` (the cut does what it claims) and correctly `MISMATCH` against a `not linear` gold.
4. **Slot mislocation in `compare_label_tuple`.** I could not make content-location and position-location disagree in a way that produced an accept. `"a) Not memoryless\nb) Causal, since it does not depend on future inputs"` correctly resolves (No, Yes) — the `not` in the trailing clause does not reach `Causal` — and `"a) Memoryless: Yes b) Causal: Yes"` against a (No, Yes) gold correctly mismatches on the memoryless slot. `_resolve_property`'s stated precedence order survived every case I built.
5. **Symbolic chained equality.** `"u(x, y) = -2x^2 + 6xy = -2x^2 + 6xy + C(y)"` — `_isolate_expression` takes the first non-LHS segment and `_drop_arbitrary` removes the trailing term; `MATCH`, correct. The multi-line preamble form (`"… is given by:\nu(x, y) = 10y^2 + 4xy"`) is handled by the `with_eq[-1]` line selection as documented. I could not construct a wrong polynomial the fragment parser accepted.
6. **Symbolic arbitrary-function policy.** `+ C(y)`, `+ C`, and `u(x,y) = C` alone all behave as specified — accept, accept, and mismatch (not UNRESOLVED) respectively.
7. **Superscript / NFKC ordering.** `normalize_chars` folds `x²` to `x^2` before NFKC as claimed; I confirmed the failure it prevents (NFKC alone yields `x2`) and found no ordering that broke a real archived continuity answer.
8. **`SUPPORT` table drift.** All 27 entries recompute from the archive. The self-critical note about the first draft's 44/61 and 47/61 figures is borne out: the measured values are 57 and 55.

---

## 5. Further probing and improvements

**Where the evidence is thinnest, ranked.**

1. **Four of six kinds cannot be falsified by this archive, and the fifth turns out not to test the property it was built for (F0, F1).** With more time I would spend all of it outside these four templates: the corpus holds 2,200 traces and only 61 were used, and the `numeric` / `multipart` families are where real wrong answers actually live. A per-kind count of archived *negative* instances should be reported next to precision — precision without it is unreadable.
2. **The origin degeneracy (F3) is an item problem the comparator is papering over.** The better fix is probably not in `compare_sequence`: make the template *always* print the origin (state the support explicitly, `y[n] = {…} for n = 1…4`, even when n = 0 is outside it). That removes 16.4% of the item pool from the undecidable class and makes `require_origin` mean what it says. Failing that, the `g.origin is None and c.origin is not None` branch must return `UNRESOLVED`, not `MATCH`.
3. **Surface priority (F2) generalises badly.** Length-first-then-break is inherited by every future label set, and it is wrong wherever one surface is lexically inside another — which is the common case for negated technical adjectives. Later phases adding class-C bindings will hit it on their first item.
4. **`answer_span` needs a stop rule (F4).** Peel-to-last is right for restatement and wrong for qualification. Cheapest general improvement: refuse to peel across a sentence opening with a discourse marker (`Note`, `However`, `Unless`, `If`), and record the peel count as an observation so the effect is visible in results instead of silent.
5. **Three magic windows with no sensitivity analysis** (`_NEG_WINDOW = 40`, `HEDGE_WINDOW = 45`, `neither…nor` 60). The archive's two `neither…nor` gaps are both 12 characters. I would sweep each bound over the full 2,200-trace archive and report the verdict-flip rate as a function of the bound: a bound whose verdicts are stable over 30–120 characters is defensible; one that is not should be replaced by clause-structural scoping rather than a character count.

**What the spec and the two documents get wrong or leave ambiguous.**

* `phase4_comparators.md:251` and `kinds.py:682` carry the trace-8 claim that F1 falsifies. It must be corrected, not softened: it is the sole empirical justification offered for the `sequence` kind's design, and `ground_truth.py`'s own label for the same trace already contradicts it.
* `_resolve_label`'s docstring claims resolution "by position, latest wins"; the code resolves by surface length. One of the two must change, and F2 argues it should be the code.
* `phase4_comparators.md:273` reports "gold itself omits the origin on 16.4% of instances" as a property of the item, and the surrounding prose describes the resulting comparison as "checked for consistency". It should state plainly that on those instances the origin is **not scored at all** and `require_origin=True` has no effect. That is the operative fact for anyone reading the gate.
* D4.2 §5's residual-risk section correctly places the symbolic rules on "route 3" — stated, adversarially exercised, not validated against real output. After F1, `require_origin` belongs in the same paragraph.

**Process note.** The brief's framing — "the archive is also the corpus the vocabulary was derived from, so 100% on it may mean nothing at all" — is right, and F0 supplies the number: the archive holds 16 wrong answers, 15 of them decided, all in two of the six kinds. The cross-pairing instrument in §4 is the cheapest way I found to get more signal out of the same 61 traces without inventing text, and it should be run as a standing check rather than as a one-off review artefact.

---
---

# Round 2 — comparator adversary

**Verdict: BLOCKED — the unreviewed module the fixes are built on reproduces two of the defects it was written to close, one clause boundary to the left, and it does so with no archive support in either direction.**

Frozen ref: `335ad7a3c2a7ec621db62180d74e286d2f59e0f1`.

---

## 1. Verdict

`tests/comparators/commitment.py` is a sound idea implemented as three magic
constants wearing the vocabulary of linguistics. It replaced a 45-character
hedge window with a **±1-clause hedge window**, and a longest-surface rule with
a last-clause rule that still resolves *within* the clause by last position —
so both of the defects it names in its own docstring survive, reachable by
changing one period to a comma (R2-F7) or by adding one clause of distance
(R2-F8). Against that it now rejects ordinary correct engineering prose
(R2-F9), most severely in `check`, the one gated kind with a numeric register.
And the whole apparatus has **essentially zero exercise on real archived text**
(R2-F10): across all 2,200 archived answer spans the hedge machinery fires
**once**, and never in any of the three kinds it gates.

Round-1 status: **F2 relocated. F3 addressed. F4 addressed. F5 addressed. F1
not addressed, and now contradicted by the F3 fix as well as by the archive.
F6 not addressed; the parser fix behind it is correct and settles F6 anyway.**

---

## 2. Independent re-derivation (round 2)

| claim under test | as re-derived | agreement |
|---|---|---|
| archive: 100% precision, 0 false accepts | reproduced (`score`); recall now **97.8%**, decided **93.4%** | **agrees** |
| D4.4: 98.6% / 98.6%, 1 false accept (`num-16b`), 1 false reject (`seq-08`) | reproduced | **agrees** |
| reviewer battery 43/43 | reproduced | **agrees** — but 43/43 is a fixed set, and R2-F7/F8 live one edit outside it |
| **F2 fix**: position first, length only at a tie | reproduced in `_label_hit`; **but resolution is still last-by-position *within* a clause, and the splitter does not cut on commas** | **DIVERGES**, R2-F7 |
| **F3 fix**: silent gold + stated origin → UNRESOLVED / MISMATCH | reproduced; archived signal trace 0 and `seq-08` are now UNRESOLVED | **agrees** — what I asked for; the cost is visible, signal recall 100% → 50% |
| **F4 fix**: peel refuses to cross a discourse marker | reproduced; F4 battery cases pass | **agrees** |
| **F5 fix**: `neither…nor` scoped to the sentence | reproduced; my 60-char escape now correctly resolves (No, No) | **agrees** |
| **F1 fix** | `kinds.py:732-737` and `phase4_comparators.md:251` are **unchanged from round 1** | **not addressed** |
| `_BRACED_RE` subscript bug | real, and the fix is right — `y[4]` is not a one-element sequence | **agrees** |
| **cross-pairing, re-run against the new clause machinery** | **762 real-text pairs, 0 false accepts, 0 false rejects, 0 UNRESOLVED** (210 linearity, 552 memory-causality) | **agrees — no regression**; the round's strongest positive result |
| hedge markers fired across all 2,200 archived answer spans | **1** (`modal: would be`, once) | **new number**, R2-F10 |
| clauses opening with a `NON_ASSERTING_OPENER`, all 2,200 spans | **20**, of which 9 are a span's final clause — **all 9 in numeric templates, none in a commitment-gated kind, none a committed answer** | **new number**, R2-F10 |

The cross-pairing re-run answers the brief's worry that the clause machinery
could break the instrument: it does not. Every real archived answer form still
scores correctly against every real gold in its template. The failures below are
all outside the shapes the archive contains — the same sentence I wrote in
round 1, now the finding rather than a caveat.

(762 rather than round 1's 782: the 20 continuity pairs are omitted here, since
`symbolic` is not commitment-gated and was out of the round's scope.)

---

## 3. Findings

### R2-F7 — CONFIRMED (blocking): **false accept.** F2 is *relocated*, not closed — the clause splitter does not cut on commas, and within a clause resolution is still last-by-position

`_label_hit` was correctly changed to rank by position. But `find_commitment`
only chooses *which clause*; `_resolve_label` is then run on that clause and
still takes the **last** surface in it. `_CLAUSE_SPLIT` cuts on `.;!?`, dashes
and enumerators — **not on commas**, and not on `unlike` — so a contrastive
mention inside the committed clause outranks the commitment exactly as before.

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:**\nThe system is **not linear**."
for c in ["## Final Answer\n**Answer:** The system is linear, unlike a nonlinear system, since both tests pass.",
          "**Answer:** The system is linear (a nonlinear system would fail additivity)."]:
    print(compare_kind('categorical', g, c).outcome, '|', c[-60:])
PY
```

→ `MATCH` twice. A candidate committing to **linear** — the wrong answer — is
credited as `not linear`.

This is my round-1 F2 case with **one period changed to a comma**. Battery case
`E-F2` uses the period form and passes; the comma form is the more natural
English of the two. `unlike` *is* in `NON_ASSERTING_OPENERS`, but `_OPENER_RE`
anchors with `^`, so it fires only when the contrast is punctuated into a clause
of its own. The fix is conditioned on punctuation the writer chooses freely.

**Impact.** The defect the module was written to close is reachable in the kind
with the largest archived trace count and zero archived negatives, by ordinary
prose. Blocking for the same reason F2 was.

---

### R2-F8 — CONFIRMED (blocking): **false accept.** `_governing_neighbour` is a fixed-width window in clauses instead of characters — a hedge one clause further away escapes it

The module's stated advance is *"scope is the clause, bidirectionally — not N
characters backwards."* `_governing_neighbour` iterates `for j in (i - 1, i +
1)`. That is not clause scope; it is **N = 1 clause** — a magic constant of the
same shape as `_NEG_WINDOW = 40` and `HEDGE_WINDOW = 45`, chosen with no
sensitivity analysis and unmeasured on real text.

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:**\nThe system is **linear**."
for c in ["**Answer:** I am not sure. Let me redo the algebra. The system is linear.",
          "**Answer:** I cannot really tell. The two tests are in step 4. The system is linear.",
          "**Answer:** It seems additive. The homogeneity check is in step 3. The system is linear.",
          "**Answer:** The system is linear. Step 4 shows the algebra. But I am not sure."]:
    print(compare_kind('categorical', g, c).outcome, '|', c[13:])
PY
```

→ `MATCH` four times. `I am not sure` and `I cannot really tell` are Reviewer
B's F1/F3 vocabulary; `UNCERTAINTY` lists both; `hedge_markers` detects both
correctly. They escape because one intervening clause puts them at distance 2.
The interposed clause is label-free and hedge-free, so nothing blocks the
governance chain — it is simply never looked for. The fourth case shows the same
escape in the trailing direction, which is the direction B's F2 was about.

Interaction with R2-F11: `clauses()` splits on decimal points and abbreviations,
so a hedge that *was* adjacent can be pushed to distance 2 by a number appearing
between it and the label.

**Impact.** B's F1/F3 are relocated rather than closed: the answer is still
"unhedged" if the hedge is far enough away, and "far enough" is now one clause
rather than 45 characters. This is exactly the Phase 3 lesson — the fix is
sound, the surface it introduced was not reviewed.

---

### R2-F9 — CONFIRMED (blocking): **over-rejection.** `is_assertion` and the hedge classes refuse ordinary correct answers; `check` is worst affected

Four of five naturally-phrased **correct** `check` answers are now UNRESOLVED:

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_kind
g = "**Answer:** The deflection is 18.0 mm, which is acceptable."
for c in ["**Answer:** The deflection is 18.0 mm. Given that the limit is 25 mm, the design is acceptable.",
          "**Answer:** The deflection is 18.0 mm, roughly 72% of the 25 mm limit, so the design is acceptable.",
          "**Answer:** The deflection is 18.0 mm. Note that the limit is 25 mm, so it is acceptable.",
          "**Answer:** The deflection is 18.0 mm, so the design is acceptable. The margin seems comfortable."]:
    v = compare_kind('check', g, c); print(v.outcome, '|', v.reason[:72])
PY
```

→ `UNRESOLVED` × 4. The same forms fail on `categorical`
(`Given that both tests pass, the system is linear` /
`Although the offset complicates matters, the system is linear` /
`Note that both tests pass, so the system is linear`) and on
`categorical[tuple]` (`a) Memoryless: Yes b) Causal: Yes. Both look immediate.`
→ UNRESOLVED on the causal slot).

Three distinct causes, each a design error rather than a coverage gap:

1. **A fronted subordinate clause is not a non-assertion.** `_OPENER_RE` matches
   the clause opener, but `Given that X, Y` and `Although X, Y` *assert Y*. The
   rule confuses the opener of a subordinate clause with the mood of the
   sentence. `note that` has the same problem when it *precedes* the answer
   rather than following it: position, not vocabulary, distinguishes a caveat
   from a preface, and `is_assertion` sees no position.
2. **`roughly` is in `EVIDENTIALS` but is not epistemic in engineering
   register.** It is a quantity approximator (`roughly 72%`, `roughly 3 mm`),
   and `check` is the one gated kind whose answers carry quantities.
   `nominally` has the same problem (`nominally 25 mm`).
3. **`_governing_neighbour` is bidirectional over appearance copulas**, so a
   trailing remark about *anything else* (`The margin seems comfortable`,
   `Both conditions look satisfied`, `Roughly 3 lines of algebra confirm it`)
   retracts a commitment it does not modify. The docstring's justification is
   `Linear. It seems.` — an anaphoric hedge — but the rule cannot tell an
   anaphoric hedge from an unrelated sentence containing `seems`.

**Impact.** `check` already has the lowest decided rate on D4.4 (61.9%). These
are false *rejects*, so they do not move precision and are invisible in the
false-accept half of the gate — but the gate is precision **and** recall, and
recall is measured on 137 hand-written cases whose phrasing was chosen by the
same author. Ownership note: this is not the pedagogy question Reviewer B owns
(whether the vocabulary is too generous to a bad answer); it is mechanism
robustness, which is mine.

---

### R2-F10 — CONFIRMED (blocking, structural): the new mechanism has **no archive support at all** — this is F0 one layer up, and it is the module's own charge against its predecessor

The module's docstring convicts the hedge blocklist of D-034: *"every hedged
adversarial case used a phrase already in the list — the corpus samples the
inside of the list it certifies."* The replacement is in the same position, and
by a wider margin.

```bash
python - <<'PY'
import sys, json, glob; sys.path.insert(0, '.')
from collections import Counter
from tests.comparators.normalize import answer_span, prepare
from tests.comparators.commitment import clauses, hedge_markers, _OPENER_RE
rows = []
for p in sorted(glob.glob('error_analysis_annotation/samples/*.jsonl')):
    for line in open(p, encoding='utf-8'): rows.append(json.loads(line))
ev, op = Counter(), Counter()
for r in rows:
    for c in clauses(prepare(answer_span(r['model_reasoning'])[0])):
        m = _OPENER_RE.match(c)
        if m: op[m.group(0).strip().lower()] += 1
        for h in hedge_markers(c): ev[h] += 1
print('rows:', len(rows), '| hedge markers fired:', sum(ev.values()), ev)
print('non-asserting openers:', sum(op.values()), op.most_common())
PY
```

→ **2,200 answer spans. Hedge markers fired: 1** (`modal: would be`, once).
**Non-asserting openers: 20**, of which 9 land on a span's final clause.

The distribution is the finding, not the total. All 9 final-clause openers are
in **numeric** templates (`flow_rates_vs_conversion`,
`gauss_law_symmetric_charge`, `hagen_poiseuille_flowrate`,
`finding_limiting_reactant`, `levenspiel_plot_interpretation`, …) — a kind
`commitment.py` does not gate. In `categorical`, `categorical[tuple]` and
`check` — the three kinds it *does* gate — the hedge probes and the opener rule
fire **zero times across every archived trace**.

So: ~250 lines of natural-language heuristics, gating three of six kinds,
validated entirely against 43 battery cases and 137 D4.4 cases **hand-written by
the two people whose findings it answers**. It cannot be falsified by the
archive in either direction — the archive contains nothing it rejects and
nothing it should reject. Every number in §2 that speaks well of it comes from
text the module's author or its reviewers wrote.

Reported as honestly as F0: I am not claiming the archive contradicts the
module. I am claiming it says nothing about it, and the gate wording must not be
read as though it did.

---

### R2-F11 — CONFIRMED (non-blocking alone; the enabler for R2-F8): `clauses()` splits inside decimals, abbreviations and initials

```python
clauses("The gain is 2.5 and the system is linear.")
# ['The gain is 2', '5 and the system is linear']
clauses("The system is linear i.e. additive and homogeneous")
# ['The system is linear i', 'e', 'additive and homogeneous']
clauses("By Prof. Smith's test the system is linear.")
# ['By Prof', "Smith's test the system is linear"]
clauses("Part b. The system is linear.")   # the enumerator rule
# ['Part', 'The system is linear']
```

No false accept follows from the split *alone* — the label usually stays with
its predicate, and `_governing_neighbour` covers the immediate spill. It matters
because it is a **distance amplifier for R2-F8**: any decimal between a hedge
and its label pushes the hedge from distance 1 to distance 2, outside the
governing window. `check` answers contain decimals by construction. The fix is
the standard guard (`(?<!\d)\.(?!\d)` plus an abbreviation exception list), and
it should land before the window question is settled, because the two interact.

---

### R2-F12 — PLAUSIBLE (minor): `_POST_NEG_RE` fires inside a parenthetical

`_POST_NEG_RE` opens `^\W*(?:is|was|would\s+be|seems)?\s*(?:the\s+)?` — `\W*`
steps over an opening bracket and the optional `the` absorbs the article, so a
parenthetical *about* wrongness negates the label it follows:

`"The system is linear (the wrong answer would be nonlinear)."` → `UNRESOLVED`
against a `not linear` gold.

The rule itself is right and its intended cases work
(`"linear is the wrong description here."` → `not linear`, correct;
`"Calling it nonlinear is incorrect; the system is linear."` → `linear`,
correct). PLAUSIBLE because the trigger phrasing is mine. The cheap tightening
is to forbid `\W*` from crossing an unclosed bracket.

---

### F1 — **NOT ADDRESSED**, and now wrong in a second way

`kinds.py:732-737` and `phase4_comparators.md:251` are unchanged from round 1
and still assert that signal trace 8 "emits exactly gold's multiset of values in
exactly gold's order" and that this is "1 of the 16 archived traces". Both
halves remain false: the answer is a **rotation**, `ground_truth.py`'s own label
says so, and 0 of 16 archived traces reach the origin-discriminating branch.

The F3 fix has now made a *second* sentence in the same passage false.
`phase4_comparators.md:273-277` still reads: *"a candidate that pins an origin
gold does not is **checked for consistency rather than punished for saying
more**."* That describes the branch as it was **before** F3 was actioned. The
code now returns `UNRESOLVED` when the values agree and `MISMATCH` when they do
not — archived signal trace 0 and `seq-08` are both demonstrations, and both
appear in the `score` output as false rejects. The document describes behaviour
the build no longer has, at the exact place a reader goes to understand what the
gate measured.

---

### F6 — not addressed; the parser fix behind it is correct, and it settles F6 anyway

`_BRACED_RE`'s bracket alternative was matching `y[4]` and `val_groups[-1]` was
taking the subscript as the answer. The guard (*a bracket group is a sequence
only if it contains a separator*) is right, and the resulting reclassification
of archived signal trace 0 from `MATCH` to `UNRESOLVED` is **correct and is
exactly what F3 asked for**: gold omits the origin, the candidate states one via
the `*-1*` form, the item cannot discriminate, and `UNRESOLVED` is the only
outcome the evidence supports.

It also settles F6 in my favour without the label being touched. Of the
archive's two `signal_operations` `MATCH`es, one is now gone *because it was
never decidable* — which is what F6 said the justification text should record.
Signal recall is 100% → **50%** and decided 100% → **75%**. That is the honest
number and it should be reported beside the gate, not absorbed into a total.

---

## 4. Falsification attempts that failed

1. **Cross-pairing, re-run against the clause machinery — 762 real-text pairs,
   0 false accepts, 0 false rejects, 0 UNRESOLVED** (210 linearity, 552
   memory-causality). The brief was right that this was the change most likely
   to break it, and it did not: `clauses`, `is_assertion`, `hedge_markers` and
   `find_commitment` are transparent on every real archived answer form, at
   every real gold in its template. This is the strongest positive evidence in
   either round, and R2-F10 is why it is not sufficient.
2. **F5 (`neither…nor`) — genuinely closed.** My 60-char escape case now
   resolves (No, No), and the sentence-scoped rule survived every longer
   interposition I built.
3. **F4 (peel into a caveat) — genuinely closed.** `answer_span` no longer peels
   into `Note: Answer: …`, and `is_assertion` catches the residue that remains
   inside the span, as the module's comment claims.
4. **`must` / `can` excluded from `MODALS` — correct.** `The system must be
   linear` and `The system can be shown to be linear` both commit, as they
   should. I could not build a hedge needing `must` or `can` that is not already
   caught by an evidential or an appearance copula.
5. **`"no"` removed from `NEGATORS` — correct.** `No, the system is linear`
   resolves to `linear` and `No, the system is not linear` to `not linear`; the
   sentential `No` no longer inverts the label it precedes.
6. **`_POST_NEG_RE` on its intended cases — correct** (see R2-F12 for the
   residual).
7. **`DISCOURSE_MARKERS` refusing a legitimate marker.** I could not find one.
   The 9 real final-clause caveats in the corpus are all genuine trailing
   commentary, correctly declined. `if` and `but` in the list worried me; both
   behave.
8. **Decimal / abbreviation splitting as a direct false accept.** R2-F11: I
   could not turn the split into an accept on its own, only into extra distance
   for R2-F8.
9. **`TASK_RESTATEMENTS` misfiring on a real answer.** Zero hits across the
   2,200 spans, and I could not construct a natural committed answer carrying
   `determine whether` / `verify whether` in the *committed* clause.

---

## 5. Further probing and improvements

**Ranked by what would change the verdict.**

1. **Close R2-F7 at the resolution step, not the segmentation step.** Adding
   commas to `_CLAUSE_SPLIT` would break `a) Memoryless: Yes, b) Causal: No`.
   The right fix is that `_resolve_label`, applied to a committed clause, must
   ignore surfaces inside a **contrastive or parenthetical constituent** — a
   mention after `unlike`, `rather than`, `as opposed to`, `not`, or inside
   brackets, is not the commitment. Small, testable, and it closes the
   parenthetical half of R2-F12 too.
2. **Replace the ±1-clause window with a real scope rule (R2-F8).** Either let
   *any* preceding label-free clause in the span carrying `UNCERTAINTY` govern
   (stated uncertainty is span-scoped, not clause-scoped — a model that says "I
   am not sure" anywhere has not committed), or scan outward until a
   label-bearing clause blocks. Whichever is chosen, **report the flip rate as a
   function of the radius** — the sensitivity analysis round 1 asked for on
   `_NEG_WINDOW`/`HEDGE_WINDOW`, which this change should not have been allowed
   to skip.
3. **Split `is_assertion` by position (R2-F9).** A `NON_ASSERTING_OPENER` at the
   *start of the span* prefaces an answer; the same opener on the *last clause*
   qualifies one. The corpus supports exactly this: all 9 real openers are final
   clauses and all 9 are caveats. Restricting the opener rule to non-initial
   clauses costs nothing on real text and recovers `Given that …, the system is
   linear`. Separately, drop `roughly` and `nominally` from `EVIDENTIALS`, or
   require them not to be adjacent to a numeral.
4. **Correct the F1 prose — both sentences — before merge.** It has now survived
   a round in which everything around it changed, and the F3 fix has made a
   second sentence in the same passage describe behaviour the code no longer
   has. It is the sentence a later reviewer will cite.
5. **The gate needs a per-kind negative count *and* a per-kind
   commitment-exercise count printed next to precision.** R2-F10 is F0 with a
   new subject, and the reason both were findable is that the score table
   reports a denominator including instances the mechanism never touched. One
   extra column — "archived spans on which any commitment probe fired" — would
   have made this visible without a reviewer.

**Process note, and it is the whole round.** Every fix from round 1 is
individually correct. Three of them are correct *and* the surface they were
built on reproduces the defect one edit away. The module states the Phase 3
lesson in its own docstring — "the corpus samples the inside of the list it
certifies" — and was then certified against 43 cases written by the two people
it was answering. A new mechanism introduced to fix review findings should carry
its own evidence base, measured on text nobody in the review wrote, before it is
allowed to gate three kinds.

---

# Round 3 — comparator adversary

**Frozen ref `3a00875b5cc5c8baa440fdd55fa131d5a96abfe6`.**

## 1. Verdict

**BLOCKED.** `commitment.py` v2 is a better mechanism than v1 and I withdraw
nothing I said about the direction of travel — but the round-2 pattern repeats a
third time: **five of the six v2 mechanisms fail on their own first
generalisation**, and four of my earlier findings are *relocated* rather than
closed. I have both directions the round asked for: correct answers the module
rejects (R3-F13, R3-F16, R3-F17, R3-F18, R3-F19) and non-commitments it accepts
(R3-F14, R3-F15). The four worst take a *correct* answer to **MISMATCH**, not
UNRESOLVED — a model error reported that did not occur, which is the failure
mode F4 and `answer_span` exist to prevent.

## 2. Independent re-derivation (round 3)

Baseline reproduced at the frozen ref before attacking:

```
reviewer_battery   64/64   every case both reviewers used
recall_corpus     351/351  9 frames x 39 archived answers, 0 false rejects
derive_vocabulary  27/27 SUPPORT entries agree with the archive
```

The **commitment census** (R2-F10) is honest and I accept it as filed. It
reports 2 hedge firings in 2,200 spans, **0 in a gated kind**, and states in its
own output that "the archive cannot falsify the commitment machinery in either
direction". That is the correct framing and it is what I asked for. One number
in it is load-bearing and is not called out: `factive: 1`. The FACTIVE class is
the *least*-supported of the three new subordinator classes — one occurrence in
2,200 spans — and it is the class R3-F16 breaks.

Both green suites are, again, closed sets. All 64 battery cases and all 351
recall frames were written by the implementer or by B and me. §4 records what I
could not break; everything in §3 was reachable from the frozen tree in under
thirty minutes with no archive access.

## 3. Findings

### R3-F13 — CONFIRMED (blocking): **over-rejection.** `is_bare_comment`'s subjectless fallback classifies ordinary technical sentences as bare epistemic comments — R2-F9c relocated

`is_bare_comment` ends with a *determiner* blocklist:

```python
return not re.match(r"^\W*(?:the|a|an|both|each|every|all|its|his|her|their)\b", stripped, re.I)
```

Anything <= 6 words that is label-free and does not open with one of those
eleven determiners is a "bare comment". The subject vocabulary of the three
gated kinds is overwhelmingly **bare abstract nouns**, which take no determiner:

| candidate (gold `linear` / `not linear`) | verdict | should be |
|---|---|---|
| `The system is linear. Superposition seems to hold.` | UNRESOLVED | MATCH |
| `The system is linear. Scaling appears satisfied.` | UNRESOLVED | MATCH |
| `The system is linear. Additivity seems fine.` | UNRESOLVED | MATCH |
| `The system is nonlinear. Homogeneity appears violated.` | UNRESOLVED | MATCH |

`is_bare_comment("Superposition seems to hold", surfaces)` -> `True`.

This is **exactly R2-F9c**, in the same register, one construction over. My
round-2 case was `... so the design is acceptable. The margin seems comfortable.`
— refused. That case now passes *only because "The margin" opens with `the`*.
Swap the determiner-headed NP for the bare noun the domain actually uses
(superposition, additivity, homogeneity, scaling, causality, stability) and the
refusal is back. The fix distinguished the two sentences by the one feature that
does not distinguish them.

`recall_corpus`'s `post-quantified` frame is 39/39, so the corpus does not see
this: its frame sentence evidently carries a determiner.

### R3-F14 — CONFIRMED (blocking): **false accept.** `_strip_label_parens` runs *before* hedge detection, so a hedge inside a label-naming parenthetical is deleted along with it — R2-F12 relocated

`find_commitment` computes `bare = _strip_label_parens(matrix, surfaces)` and
then runs `hedge_markers(bare)`. Any parenthetical that names a label is removed
**entirely**, including whatever epistemic material it carries:

| candidate (gold `linear`) | verdict | should be |
|---|---|---|
| `The system is linear (or possibly nonlinear).` | **MATCH** | UNRESOLVED |
| `The system is linear (or nonlinear, I am not sure).` | **MATCH** | UNRESOLVED |
| `The system is linear (though it may be nonlinear).` | **MATCH** | UNRESOLVED |

Row 2 contains a verbatim `INABILITY` marker and is scored a clean commitment.

The R2-F12 fix is right for the case it was built on — `linear (a nonlinear
system would fail additivity)` is a gloss whose `would` is about a hypothetical
answer — and wrong for the case one word away: `linear (or possibly nonlinear)`
offers an *alternative*, and the parenthetical is the whole hedge. The mechanism
has no way to tell a gloss about the other label from an alternative between
both, because it looks only at whether a label surface is present, which is true
of both. The docstring's stated invariant — "a parenthetical naming **no** label
is a genuine qualifier and stays" — is the wrong test: `(or possibly nonlinear)`
names a label *and* is a genuine qualifier.

Minimum fix: strip the parenthetical for **label resolution** only, and run
`hedge_markers` over the *unstripped* matrix.

### R3-F15 — CONFIRMED (blocking): `_MAX_COMMENT_WORDS = 6` is the magic constant `(i-1, i+1)` was, and one word escapes it

The distance bound left `_governing_comment` and reappeared in
`is_bare_comment` as a **word count**. It is escapable exactly the way the old
one was — by making the qualification one word longer:

| candidate (gold `linear`) | words | verdict |
|---|---|---|
| `The system is linear. I am not sure.` | 4 | UNRESOLVED (correct) |
| `The system is linear. I am not entirely sure about that.` | 7 | **MATCH** (wrong) |
| `The system is linear. Honestly I cannot really tell for certain.` | 8 | **MATCH** (wrong) |

`not entirely sure` and `cannot really tell` are both **declared** in
`INABILITY`; the module tokenises them, then discards the clause before it
looks. The docstring says the magic constant "leaves rather than shrinking".
It did not leave; it changed units, from clauses to words, and there are now
**three** unjustified integers on this path — `_MAX_COMMENT_WORDS = 6`,
`range(4)` in the excision loop, and `_NEG_WINDOW = 40` in `kinds.py`
(untouched, still character-scoped, and now operating on text `segment()` has
spliced — see R3-F19). None carries a sensitivity analysis. That is the same
request round 1 made about `_NEG_WINDOW`/`HEDGE_WINDOW` and it has now been
skipped three times.

### R3-F16 — CONFIRMED (blocking): **MISMATCH on a correct answer.** The new FACTIVE class re-opens round-1 F4, from inside `segment()`, and contradicts `normalize.DISCOURSE_MARKERS`

Two hand-built opener lists in two modules read the same opener in **opposite
directions**:

* `normalize.DISCOURSE_MARKERS` contains `note`, `notes` — a `Note ...` sentence
  is a *qualification*, so `answer_span` must not peel into it. That list is the
  round-1 F4 fix.
* `commitment.FACTIVE` contains `note that`, `notes that`, `recall that` — a
  `Note that X` clause **asserts X**, so `find_commitment` may select it.

So `answer_span` correctly declines to peel into the caveat, and `segment()`
then hands the caveat to the comparator as the asserting clause. Because
`find_commitment` takes the **last** labelled clause, the caveat outranks the
answer:

```
## Final Answer
**Answer:** The system is linear.

Recall that a nonlinear map fails additivity.
```

gold `linear` -> **MISMATCH `not linear`**. `segment("Recall that a nonlinear
map fails additivity")` -> `'a nonlinear map fails additivity'`.

Likewise `Note that the answer is nonlinear with an offset` after a committed
`linear` -> **MISMATCH**.

The battery's two `E-F4` cases both use `Note:` (colon), which `_BARE_QUAL_RE`
catches. They cover precisely the form v2 suppresses and miss the form v2 newly
*asserts*. The docstring claims the two are "distinguished by what follows,
which is the distinction E's R2-F9a said version 1 could not make" — the
distinction exists in `commitment.py` and does not exist in `normalize.py`, and
the two modules are on the same path.

The underlying error is scope: a factive asserts its complement **about the
world**, but `Recall that a nonlinear map fails additivity` is a general remark
about the *other label*, not a verdict on *this system*. FACTIVE has archive
support of 1 (the census, §2). This class should not gate three kinds on one
observation and a contradicted opener list.

### R3-F17 — CONFIRMED (blocking): **over-rejection.** A fronted concessive without a comma has "no matrix to assert" — B's R2-F1.1 relocated to an orthographic condition

`segment()` locates the matrix of a fronted concessive with `c.find(",")` and
returns `None` when there is none. The comma is punctuation, not grammar:

| candidate (gold `linear`) | verdict |
|---|---|
| `As expected the system is linear.` | UNRESOLVED |
| `Since both tests pass the system is linear.` | UNRESOLVED |
| `Given that additivity holds the system is linear.` | UNRESOLVED |
| `Because superposition holds the system is linear.` | UNRESOLVED |
| `Although the gain varies the system is linear.` | UNRESOLVED |

All five refuse with `the clause is a bare subordinate clause with no matrix to
assert` — v1's refusal message, on v1's construction, gated now by whether the
model typed a comma. `recall_corpus`'s `fronted-concessive` and `fronted-given`
frames are 39/39 each because both frames insert the comma; the corpus tests the
punctuation it supplies.

The first comma is also the wrong boundary whenever the subordinate clause
contains one of its own, which is common in this register.

### R3-F18 — CONFIRMED (blocking): **MISMATCH on a correct answer.** `answer_span`'s qualifier guard has both a short opener list *and* a live 80-character distance bound

The `Answer\s*:` colon fix is correct and I confirm it (`The system is linear
(the wrong answer would be nonlinear)` no longer peels). The guard *behind* it
did not get the same attention. `_QUALIFIER_RE` is:

```python
r"(?:^|[.;!?\n])\s*(?:\*+\s*)?(?:" + DISCOURSE_MARKERS + r")\b[^.;!?\n]{0,80}$"
```

Two independent defects, both reachable:

**(a) The opener list is missing openers `commitment.py` itself declares.**
`given`, `when`, `for`, `strictly`, `whether` are absent from
`DISCOURSE_MARKERS`, while `given` is in `commitment.CONCESSIVE` and `given
that` is in `commitment.FACTIVE`. Two lists of the same linguistic object,
maintained separately, diverging:

```
"**Answer:** The system is linear. Given a nonzero offset, the answer is nonlinear."
  answer_span -> 'nonlinear.'   compare -> MISMATCH   (gold: linear)
"**Answer:** The system is linear. When the gain varies, the answer is nonlinear."
  answer_span -> 'nonlinear.'   compare -> MISMATCH
"**Answer:** The system is linear. For a squaring element, the answer is nonlinear."
  answer_span -> 'nonlinear.'   compare -> MISMATCH
```

Note the marker doing the peeling here is `The\s+answer\s+is:?\s*`, whose colon
is **still optional** and which outranks the bare `Answer\s*:` in priority. The
round's stated fix hardened the low-priority marker and left the higher-priority
one matching the ordinary English phrase "the answer is".

**(b) `{0,80}` is a character distance bound on the guard for my own F4.**
Push the qualifying opener more than 80 characters from the marker and the guard
silently stops applying, even with a listed opener:

```
"**Answer:** The system is linear. Note that for a system with a constant
 additive offset term applied to the output, the answer is nonlinear."
  answer_span -> peels to 'nonlinear.'   compare -> MISMATCH   (gold: linear)
```

R2-F8 said a fixed-width window is not scope. That claim was accepted and the
window was removed from `commitment.py`. **The identical construct survives, in
character units, in the module that runs first and decides what `commitment.py`
ever sees.** Removing a distance bound from the consumer while leaving one on
the producer does not close the finding.

### R3-F19 — CONFIRMED (blocking): **MISMATCH on a correct answer.** `_TRAILING_SUB_RE`'s excision closes on the wrong comma and splices a fragment of the subordinate clause into the matrix

```python
close = c.find(",", m.end())
c = (c[:m.start()] + c[close:]).strip() if close != -1 else c[:m.start()].strip()
```

`m.end()` is the end of the *subordinator*, so `close` is the first comma
**anywhere after it** — which is the subordinate clause's own internal comma
whenever it has one, not its closing comma. The excision then cuts mid-clause
and splices the remainder into the matrix:

```
"The system is nonlinear, because the squaring term, which we verified, dominates"
  -> "The system is nonlinear, which we verified, dominates"
"The system is linear, since the map, a sum of scaled inputs, is additive"
  -> "The system is linear, a sum of scaled inputs, is additive"
```

Both are ungrammatical splices. When the spliced-in fragment carries the other
label, the splice **flips the verdict**, because `_resolve_label` is still
last-by-position within the matrix:

```
"The system is linear, since the check, a nonlinear probe, passed."
  segment -> "The system is linear, a nonlinear probe, passed"
  gold linear -> MISMATCH 'not linear'
```

A correct answer scored wrong. Note the interaction with `_NEG_WINDOW = 40`
(R3-F15): `_negated` now runs its 40-character window over text `segment()` has
**cut and rejoined**, so the window's contents no longer correspond to anything
the model wrote. Character-scoped negation over surgically spliced text is not a
scope test at all.

The loop bound is the third magic integer. Five comma-subordinators exhaust it
and the fifth survives into the matrix:

```
"The system is linear, since additivity holds, because scaling holds, given that
 both pass, as the map is a sum, unlike a nonlinear map, and that settles it"
  segment -> "The system is linear, unlike a nonlinear map, and that settles it"
  gold linear -> MISMATCH
```

Contrived as written, but `range(4)` is unjustified and the failure is silent —
the loop exits, the caller cannot tell an exhausted excision from a completed
one, and there is no `UNRESOLVED` on the path.

### R3-F20 — CONFIRMED (non-blocking): three of the four lookbehinds added for R2-F11 are **no-ops**; the fix is carried entirely by `_ABBREV`

```
(?<!\d)\.(?!\d)(?<!\bi\.e)(?<!\be\.g)(?<!\bcf)(?<![A-Z])
```

The decimal guard is correct — `(?<!\d)` precedes the `\.`. The other four
lookbehinds are placed **after** `\.` has been consumed, so each inspects text
ending at the dot, one character past what it was written to test:

* `(?<![A-Z])` examines the `.` itself, which is never `[A-Z]`. **Always true.**
  Demonstrated: `clauses("The gain is set by R. The system is linear.")` ->
  `['The gain is set by R', 'The system is linear']`. The docstring's comment
  `# . but not an initial` describes behaviour the regex does not have.
* `(?<!\bcf)` examines `f.`, never `cf`. **Always true.** Masked only because
  `_ABBREV` stashes `cf.` first.
* `(?<!\bi\.e)` / `(?<!\be\.g)` examine `.e.` / `.g.`. Same, same mask.

R2-F11 is therefore **half not addressed**: the abbreviation half works and does
so entirely through `_ABBREV`, a closed hand list of eleven items (missing at
least `approx.`, `Sec.`, `Ref.`, `Tab.`, `resp.`, `Ch.`, `al.`, `w.r.t.`); the
initials half does not exist. Either delete the three dead lookbehinds and say
so, or move them before the `\.`.

### R3-F21 — CONFIRMED (minor): `label_segment` is dead code whose `_strip_label_parens` call has different semantics from the live one

`label_segment` (`commitment.py:414`) has no callers anywhere in `tests/` or
`docs/`. It calls `_strip_label_parens(matrix)` with **no surfaces**, which
takes the `if not surfaces` branch and strips *every* parenthetical — including
the label-free ones that `find_commitment`'s call deliberately preserves to keep
`Linear (unclear)` non-committal (B, round-1 F2). A dead function that
contradicts a live invariant is a trap for the next editor. Delete it, or give
it the same call.

### R3-F22 — CONFIRMED (minor): `_TASK_RE` is an unanchored `search` over the whole clause

`determine whether`, `check whether`, `test whether`, `verify whether` are
matched anywhere in the clause, so a clause that *performs* the check is refused
as a *restatement* of it:

```
"Both tests verify whether the system is linear and both pass."  -> UNRESOLVED
"We check whether the system is linear: it is."                  -> UNRESOLVED
```

Damage is limited because `find_commitment` falls back to earlier labelled
clauses, so this only kills single-clause answers — but it does kill them.
Anchor it to the clause opener, as `_HYPO_RE`/`_CONC_RE`/`_FACTIVE_RE` are.

## Status of every round-1 and round-2 finding

| finding | claimed | my assessment |
|---|---|---|
| F0 | (open, structural) | **not addressed**; correctly carried as risk. Superseded by R2-F10 / R4-10. |
| F1 | now addressed | **addressed — complete and correct.** Both false halves are retracted in `kinds.py:768-791`, the measured `0 of 16` replaces the fabricated `1 of 16`, and `require_origin` is relabelled a route-3 rule in the docstring *and* in `phase4_comparators.md:263-264` citing D4.2 §2.4. The second sentence I flagged in round 2 is gone. I close F1. |
| F2 | closed r1 | **addressed** (position-first ranking holds). |
| F3 | closed r1 | **addressed.** |
| F4 | closed r1 | **RELOCATED — see R3-F16.** The `answer_span` guard holds; the caveat now reaches the comparator through `segment()`'s FACTIVE class instead, and through the omissions and the 80-char bound in `_QUALIFIER_RE` (R3-F18). |
| F5 | closed r1 | **addressed** (`neither ... nor` survives the new excision; verified §4). |
| F6 | settled by `_BRACED_RE` | **addressed**, as I argued. |
| R2-F7 | closed | **addressed.** Comma-level segmentation genuinely closes it; my r2 case is 2/2 and I could not revive it. |
| R2-F8 | closed | **RELOCATED.** The clause bound is gone from `_governing_comment` — that part is real. But the bound reappears as `_MAX_COMMENT_WORDS = 6` (R3-F15) and `range(4)` (R3-F19), and a *character* distance bound on the same path survives untouched in `normalize._QUALIFIER_RE` (R3-F18b). Three bounds where there was one; none has a sensitivity analysis. |
| R2-F9a | closed | **addressed** for the comma'd form; **RELOCATED to the comma-less form** (R3-F17). |
| R2-F9b | closed | **addressed.** `roughly 72% of the 25 mm limit` no longer hedges; `roughly speaking` still does. Verified. |
| R2-F9c | closed | **RELOCATED — see R3-F13.** Directionality and the bare-comment gate fix my exact sentence; the same sentence with a bare-noun subject still fails. |
| R2-F10 | not fixable, machine-reported | **addressed as far as it can be.** The census is honest, states its own limits in its own output, and is carried as R4-10. Accepted. One remark: `factive: 1` is the weakest support in the census and it gates R3-F16. |
| R2-F11 | closed | **half not addressed — see R3-F20.** Three of four lookbehinds are no-ops; the abbreviation list is doing all the work. |
| R2-F12 | closed | **RELOCATED — see R3-F14.** Stripping the parenthetical stops it *outranking*, and starts it *deleting a hedge*. |

## 4. Falsification attempts that failed

Recorded so the clean parts of the verdict mean something. Each of these I
expected to break and could not:

1. **`neither ... nor` under the new excision.** R2-F5's fix is the case the
   excision loop was explicitly built not to break, and it holds:
   `neither memoryless, because the output depends on past inputs, nor causal`
   -> MATCH `(No, No)`. The excise-don't-truncate design is correct and the
   docstring's account of why is accurate.
2. **Defeating HYPOTHETICAL through a concessive wrapper.**
   `Although X holds, if Y then the system is linear` -> refused, `the matrix
   clause is itself hypothetical`. The re-check after the concessive branch is
   there and works.
3. **Reviving R2-F7.** Every comma variant of my round-1 F2 case I tried now
   resolves correctly. Comma-level segmentation genuinely closed it.
4. **Reviving R2-F9c's derivation opener.** `Imagine a scaled input a*x[n]; the
   output scales, so the system is linear` -> MATCH. The direction rule and the
   `NEIGHBOUR_CLASSES` exclusion of `verb of opinion` are both correct.
5. **`must` / `can` as hedges.** Still correctly excluded.
6. **Deleting a negator by excision.** I tried several shapes where excising a
   comma'd subordinate should have stranded a label outside its negator's scope
   (`The system is not, given the squaring, linear`; `The system is not linear,
   unlike an LTI system, since it squares the input`). All resolved correctly —
   in one case by an accidental double negation, which is luck, not design, and
   is why R3-F19's interaction with `_NEG_WINDOW` still worries me.
7. **The `Answer\s*:` colon fix itself.** `The system is linear (the wrong
   answer would be nonlinear)` no longer peels. The fix is correct; my objection
   in R3-F18 is to the *unchanged* `The\s+answer\s+is:?` beside it.
8. **The `Note:` form of round-1 F4.** Still correctly suppressed by
   `_BARE_QUAL_RE`. It is only `Note that` / `Recall that` that broke.
9. **The census.** I tried to find a number in it that flattered the mechanism.
   There isn't one; it is the most honest artefact in the module.

## 5. Further probing and improvements

1. **Stop hand-listing openers in two places.** `normalize.DISCOURSE_MARKERS`
   and `commitment.{HYPOTHETICAL,CONCESSIVE,FACTIVE,BARE_QUALIFIERS}` are five
   lists of one linguistic object, maintained independently, and R3-F16 and
   R3-F18a are both *divergences between them*. Derive `DISCOURSE_MARKERS` from
   the commitment classes, or make `answer_span` call `segment()`. One list, one
   semantics.
2. **A factive must assert about the *subject under test*.** `Note that the
   system is linear` and `Recall that a nonlinear map fails additivity` are not
   the same act. Requiring the factive complement to name the item's subject —
   or, cheaper and defensible, refusing to let a factive clause **outrank** an
   earlier unhedged commitment — closes R3-F16 without a new list.
3. **Order the parenthetical strip after hedge detection** (R3-F14). Two lines,
   and it is the only fix in this report that costs nothing elsewhere.
4. **Replace the determiner blocklist with a positive test** (R3-F13). A bare
   comment is one whose subject is anaphoric *or* absent; "has a subject that is
   not a determiner-headed NP" is not that test. If a positive test is too
   expensive, invert the burden: govern only from `_ANAPHORIC_SUBJECT` matches
   and drop the subjectless fallback, which is where every false positive I
   found came from.
5. **Print the integers, or justify them.** `_MAX_COMMENT_WORDS = 6`,
   `range(4)`, `_NEG_WINDOW = 40`, plus `{0,80}` in `normalize`. This is the
   fourth time across three rounds that a magic bound has been replaced by
   another magic bound and the sensitivity analysis has been skipped. **Report
   the flip rate as a function of each**, on the 351-frame corpus, as a table in
   the module. If a bound is flat over its range, say so and keep it; if it is
   not, it is a tuning parameter fitted to reviewer-authored text and it belongs
   on R4-10.
6. **The recall corpus tests its own frames.** Nine frames, 39 answers, 351
   cases, 100% — and R3-F13 and R3-F17 both live *inside* frames the corpus
   already has, differing by a determiner and by a comma. Vary the frames
   mechanically: for each frame, emit a comma-less variant, a bare-noun-subject
   variant, and a parenthetical variant. That is a one-afternoon change and it
   would have caught five of the ten findings above.

**Process note.** Round 2's note said a new mechanism should carry its own
evidence base before gating three kinds. Version 2 answered that by *measuring
the absence* honestly (the census) rather than by building the evidence — which
is the right thing to say and not the same as the right thing to do. The result
is a third round in which every individual fix is defensible and the surface
they are built on reproduces the defect one determiner, one comma, one word, or
one module away. The module is now 512 lines of hand-built English grammar
gating three answer kinds, validated on 415 sentences all written by its three
critics, with three unjustified integers and five divergent opener lists. That
is not a comparator; it is a parser, and it should either be scoped down to what
the archive can falsify or moved to R4-10 in its entirety.

---

# Round 4 — is the mechanism converging?

## 1. Verdict

**UNBOUNDED — reduce the mechanism's scope. Do not spend a round 5 on these
five findings.** I found both halves of the gate again (3 false rejects, 2
false accepts, plus a cross-module contradiction), all five in 25 minutes, and
**all five are one-word or one-clause mutations of cases already in the passing
suite.** The defect *classes* are not closing; the *cost of finding an instance*
is falling. That is the signature of an unbounded surface, not a converging one.

The decisive number is not a finding. It is this: **deleting the hedge-governance
layer costs 0 of 351 positive recall cases, 0 of 2,200 archived spans, and 30
synthetic negative controls the reviewers wrote.** That layer is ~250 lines,
5 probe classes, 5 word lists, and **13 of my 20 findings across four rounds.**

## 2. Independent re-derivation (round 4)

Working tree at the frozen ref. Baselines reproduce: `reviewer_battery` 81/81,
`recall_corpus` 507/507 (351 positive + 156 negative), `score` gate PASS
(archive precision 100.0% / recall 97.8%; adversarial 98.6% / 98.6%).

Census, re-read rather than re-derived: **2 hedge markers in 2,200 spans**
(1 modal, 1 stated uncertainty), **0 in a commitment-gated answer type**.
One number I had not separated before, and it changes my recommendation:
**subordinator openers fire 43 times** (concessive 31, hypothetical 11,
factive 1). So the module is **two mechanisms with opposite evidence**:
segmentation is observed 43 times in the archive and is falsifiable; hedge
governance is observed **zero** times in the kinds it gates and is not.
That asymmetry is the basis of §3's recommendation, and it was invisible while
"the commitment machinery" was discussed as one object.

**Ablation, the measurement this phase was missing.** I stubbed
`_hedges_governing` and `_governing_comment` to return nothing, leaving
segmentation untouched, and re-ran the corpus:

| | full | hedge layer ablated |
|---|---|---|
| positive frames (9 x 39) | 351/351 | **351/351** |
| `neg-hypothetical`, `neg-task-restatement` | 78/78 | **78/78** |
| `neg-explicit-hedge` | 39/39 | 24/39 |
| `neg-question` | 39/39 | 24/39 |
| **total** | 507/507 | **477/507 (94.1%)** |

Ablating costs **30 cases, all in two negative frames, all reviewer-authored,
and none on the archive.** It costs **zero** false rejects. The two battery
cases it loses are `E-R3-F15` and `B-R3-F2` — mine and B's, written to specify
this layer. The layer is validated exclusively by the cases written to specify
it, which is the circularity of D4.4 and of round-1 F1 for the third time.

## 3. RECOMMENDATION — read this before the findings

> ### The mechanism is **unbounded**. Reduce its scope; do not patch it again.
>
> Natural-language commitment detection is not finishable by a rule set at this
> budget, and four rounds are enough evidence to stop. Concretely:
>
> **3a. Delete the hedge-governance layer from the gate.** Remove
> `_hedges_governing`, `_governing_comment`, `is_bare_comment`,
> `_ANAPHORIC_SUBJECT`, `_SUBJECTLESS`, `NEIGHBOUR_CLASSES`, `_COMMENT_CLASSES`
> and the five `PROBES` **from the decision path**. Cost, measured, not
> estimated: 0 archive verdicts, 0 positive recall, 30 synthetic negatives.
> Benefit: 13 of my 20 findings and 4 of this round's 5 become unreachable,
> because the code they live in is gone.
>
> **3b. Keep segmentation, which has 43 archive instances and does real work.**
> `clauses()`, the hypothetical/concessive/factive distinction, `TASK_RESTATEMENTS`
> and `conjuncts()` stay. They are falsifiable against the archive; the hedge
> layer is not. **R4-F23 below is a live segmentation defect and it should be
> fixed** — that is a two-line change, and it is the only fix I would spend a
> round on.
>
> **3c. Make hedges advisory, not gating.** Keep `hedge_markers()` as a
> reporting function. If a marker fires anywhere in an answer span, attach it to
> the verdict as an annotation and let it route to human review. **It must never
> convert a MATCH into UNRESOLVED.** Every one of my false rejects across four
> rounds — R2-F9a/b/c, R3-F13, R3-F17, R4-F23, R4-F24, R4-F25, R4-F27 — is a
> hedge or a segment silently overturning a correct verdict. Removing the
> silence removes the harm without removing the signal. Given 2 markers in
> 2,200 spans, the human-review queue this creates is **two items.**
>
> **3d. State the residual honestly and close the phase.** Spec §4.4 requires a
> hedged answer to be non-committal. After 3a-3c that requirement is *reported*
> rather than *enforced*, and the phase should say so on R4-10 in those words,
> rather than claim an enforcement it has never once exercised on real output.
>
> **What I am explicitly NOT recommending.** Not "one or two more rounds"
> (§4 says why), and not "ship it with the residual stated" — I would sign that
> for segmentation alone, but not while a layer with zero observed instances can
> silently mark a correct answer wrong. 3c is what makes "ship it" honest.

### R4-F23 — CONFIRMED (blocking): **over-rejection.** R3-F17's comma fix takes the *first* comma unconditionally, so a comma-less fronted concessive with any later comma loses its matrix

`segment()` does `comma = c.find(","); c = c[comma+1:] if comma != -1 else c[cm.end():]`.
The R3-F17 fix made the comma *optional*, but it did not make it *the concessive's
own* comma. Any later comma in the clause is taken as the boundary:

```
"Since both tests pass the system is linear."                  -> MATCH   (R3-F17 case, fixed)
"Since both tests pass the system is linear, as expected."     -> UNRESOLVED
"Although the algebra is fiddly the system is linear,
                          which the definitions require."      -> UNRESOLVED
"Since both tests pass, the system is linear, as expected."    -> MATCH   (control)
```

`segment()` returns `('as expected', '')` and `('which the definitions require', '')`.
The label is gone; the reason given is "the label is inside a backgrounded or
parenthetical part". **Two words — `, as expected` — flip a correct answer.**
This is R3-F17 relocated: the orthographic condition moved from "a comma must
exist" to "the first comma must be the right one". Note the second case combines
two frames the corpus already ships (`fronted-concessive` and
`appositive-relative`) — the corpus tests them **only separately**, and that is
why 507/507 does not see it.

**Fix (the one I would spend a round on):** when there is no comma before the
matrix, the comma-less branch must be taken. Search for the boundary rather
than `find(",")` — e.g. take the first comma **only if** the text before it
contains no label surface, else fall back to `c[cm.end():]`.

### R4-F24 — CONFIRMED (blocking): **over-rejection.** Existential `there` is in the anaphoric-subject allowlist, so a remark about the *algebra* retracts the verdict

`_ANAPHORIC_SUBJECT` lists `it|this|that|these|those|i|we|there`. Existential
`there` is an **expletive**, not an anaphor: in "there seems to be an issue with
my algebra" the logical subject is *an issue with my algebra*, a full NP, and the
clause is a claim about something else — exactly what `is_bare_comment` exists to
exclude.

```
"There seems to be an issue with my algebra, but the system is linear."  -> MATCH
"The system is linear, but there seems to be an issue with my algebra."  -> UNRESOLVED
"The system is linear, but the algebra seems fiddly."                    -> MATCH   (control)
```

The **same proposition in the same sentence flips on clause order**, and the
non-expletive paraphrase is credited while the expletive one is not. R3-F13
relocated: the blocklist became an allowlist, and the allowlist admits one word
it should not. Removing `there` fixes these three cases and nothing else; it does
not fix R4-F25, which is why 3a is the real answer.

### R4-F25 — CONFIRMED (blocking): **over-rejection.** The bare-comment test scopes by *subject*, never by the *proposition the hedge is about*

This is the through-line, and it is why R4-F24 is not worth fixing on its own.
`_hedges_governing`'s stated discriminator is "*whose* confidence is qualified".
The implementation tests only the **syntactic subject of the hedge-bearing
unit** — never what the hedge is predicated of. Any first-person or anaphoric
clause carrying an INABILITY marker retracts the verdict, whatever it is about:

```
"The system is linear. We cannot determine the ROC without more information."    -> UNRESOLVED
"The system is linear. I am not sure what the professor means by 'input space'." -> UNRESOLVED
```

Both commit to *linear* and disclaim something else — a region of convergence,
a piece of the prompt's wording. A grader credits both. The docstring's own
contrast pair ("the speaker's" vs "the reasoning's") is a distinction about
**propositional content**, and subject-form cannot express it: *I am not sure*
and *I am not sure what "input space" means* have the same subject and opposite
scope. No allowlist over subjects closes this, because the distinguishing
information is not in the subject. **This is R2-F9c and R3-F13 at their third
location, and it is the finding that makes the surface unbounded.**

Reviewer B reached the same layer independently this round (B's R4-F2: the hedge
rule's proxy refuses committed answers at 0/39 on four ordinary frames). Two
reviewers converging on *the same layer* from different directions, in the round
that asked whether the mechanism converges, is itself evidence for 3a.

### R4-F26 — CONFIRMED (blocking): **false accept.** The trailing-question withdrawal is guarded by `not _label_hit`, so the *more explicit* withdrawal is credited

The new rule fires only when the trailing question is **label-free**:

```python
if follows and text.rstrip().endswith("?") and not _label_hit(text, surfaces):
    marks.append("withdrawn by a following question")
```

So the shipped negative control is caught and its stronger form is not:

```
"The system is linear. Or is it?"                        -> UNRESOLVED  (neg-question, passing)
"The system is linear. Or is it nonlinear?"              -> MATCH       false accept
"The system is linear. But is the system really linear?" -> MATCH       false accept
```

Naming the alternative makes the retraction *less* ambiguous and the comparator
*more* confident. `neg-question` is a shipped corpus frame and case 2 is that
frame **plus one word**. The guard exists because `segment()` already refuses a
label-bearing question — but that only stops the question from *committing*; it
does nothing to stop the earlier clause committing, so the withdrawal is
dropped on the floor. R3-F16's shape: one rule's precondition is another rule's
blind spot.

### R4-F27 — CONFIRMED (blocking): **over-rejection**, and a *new* cross-module contradiction in the function added this round

`suspends_what_follows` returns `"the answer is a question"` for any preface
ending in `?`. Restating the question before answering is among the most common
things a model does:

```
"a) No, b) No"                                       -> MATCH
"Is the system memoryless and causal? a) No, b) No"  -> UNRESOLVED  "the answer is a question"
"So which is it? a) No, b) No"                       -> UNRESOLVED
"The two tests are below. a) No, b) No"              -> MATCH       (control)
```

And the two paths **disagree about the same construction**: on the categorical
path `"Is the system linear? Yes, it is linear."` is a **MATCH**, because
`clauses()` splits at the `?` and the following clause commits. The enumerated
path reads the same preface and refuses. This is R3-F16 exactly — two code paths
contradicting each other about one construction — reproduced inside the function
written **this round** to fix a different problem. `suspends_what_follows`
already exempts a preface closed by `.`; a preface closed by `?` is the same
case and needs the same exemption.

## 4. Falsification attempts that failed — and the convergence measurement

I tried to break segmentation in five further ways and could not: nested
concessives (`"Although X, although Y, Z"`), the `neither … nor` excision with
an interpolated comma (R3-F19's fix holds), `as opposed to` / `rather than`
trailing subordinators, factive-after-commitment ordering (R3-F16's fix holds),
and the abbreviation guard (`i.e.`, `2.5`, initials). **R3-F19, R3-F16 and
R3-F20 are genuinely closed.** That is real progress and I want it on the record:
segmentation minus R4-F23 is in good shape, which is why 3b keeps it.

**Now the question the round was asked.** Are these findings (a) the same kind,
one small edit from an existing case, or (b) new classes?

| round-4 finding | layer | relocation of | distance from a **passing** shipped case |
|---|---|---|---|
| R4-F23 | segmentation | R3-F17 | battery case `E-R3-F17` **+ `, as expected`** (2 words) |
| R4-F24 | hedge gov. | R3-F13 | frame `post-quantified`, trailing clause swapped |
| R4-F25 | hedge gov. | R2-F9c, R3-F13 | frame `post-quantified`, trailing clause swapped |
| R4-F26 | hedge gov. | R3-F16 | frame `neg-question` **+ one word** |
| R4-F27 | hedge gov. / cross-module | R3-F16 | frame `neg-question`, moved to the preface |

**Answer: (a), unanimously, and worse than in round 3.** Every finding is a
relocation of a named earlier finding, and **5 of 5 are one-or-two-token
mutations of cases the suite already contains and passes.** Not one is a new
class.

**Is the residual shrinking, holding, or growing?** The residual is *holding*;
what has changed is the **cost of drawing from it**, and it has collapsed:

| round | findings | how they were found |
|---|---|---|
| 2 | 6 | full session, reading the module cold |
| 3 | 10 | full session, targeted constructions |
| 4 | **5** | **25 minutes, mechanically mutating the suite's own cases by one word** |

A surface where a fixed budget yields a roughly constant number of defects, and
where the search that finds them gets *cheaper* each round because the previous
round's fixes leave one-token neighbours exposed, is not converging. Three rounds
of "every CONFIRMED finding actioned" produced a suite that is 100% green on 588
cases and 0% robust to editing those cases. **The suite measures its cases, not
the mechanism** — which is the R2-F10 charge, one level up, for the third time.

The strongest single argument against a round 5 is R4-F25. R4-F23, F24, F26 and
F27 all have two-line fixes and I could write them. R4-F25 does not: it needs to
know *what a hedge is about*, which is propositional content, and no list of
subjects, determiners, classes, positions or distances encodes that. Rounds 2, 3
and 4 have each tried a different **syntactic proxy** for a **semantic**
relation — character distance, clause distance, segment membership, conjunct
membership, subject form — and each proxy has failed on a sentence one word from
one it handles. That is not a sequence converging on an answer; it is a sequence
of equivalent-cost approximations to something outside the method's reach.

## 5. Further probing and improvements

1. **Adopt 3a-3d.** Everything else in this list is subordinate to it.
2. **If 3a is rejected, fix R4-F23 first and alone.** It is the only round-4
   finding in the layer worth keeping, it is two lines, and it costs nothing
   elsewhere. Leave F24-F27 unfixed and state them, rather than growing the
   surface by four more patches.
3. **Make the corpus a mutation harness, not a frame list.** My round-3 §5.6
   asked for mechanical frame variation and it was not built; every round-4
   finding is a case that harness would have generated. Minimum viable version:
   for each of the 13 frames, emit variants with (i) the comma deleted, (ii)
   `, as expected` appended, (iii) the trailing clause replaced by a first-person
   disclaimer about a *different* object, (iv) label surfaces inserted into
   label-free frames. That is ~50 lines and it would have found F23, F24, F25 and
   F26 without me.
4. **Publish the ablation table in the module docstring.** "This layer changes 0
   of 2,200 archived verdicts and 30 of 507 corpus cases, all of them negative
   controls we wrote" is the most useful sentence anyone can write about
   `commitment.py`, and it belongs where the next reader will see it.
5. **Retire the phrase "commitment machinery" from the spec.** It names two
   mechanisms with opposite evidence (43 archive instances vs 0) and opposite
   correct dispositions. Naming them separately — *clause segmentation* and
   *hedge governance* — is what made this round's recommendation visible, and it
   should have been visible in round 2.

**Reproduction.** All figures above come from the working tree at the frozen ref
with `PYTHONPATH` set to the repo root:
`reviewer_battery`, `recall_corpus`, `score`, `derive_vocabulary` for the
baselines and the census; for the ablation, stub
`commitment._hedges_governing` to return `[]` and `commitment._governing_comment`
to return `""` **before** importing `recall_corpus`, then re-run its `build()`
through `answer.compare_template`. Segmentation is left untouched by the stub,
which is what makes the 351/351 positive column meaningful.

**Process note, and my last one.** Rounds 1-3 asked "is this defect fixed?" and
the answer was always yes; the phase kept moving because each answer was true.
Round 4 asked a different question and got a different kind of answer: the fixes
were real and the trend they belong to is flat. I have filed 20 findings against
this module and I would file a 21st in another 25 minutes, and that fact — not
any one of the 20 — is the finding. The phase does not need a sixteenth patch.
It needs to stop enforcing a policy it has never once observed being violated,
report it instead, and close.
