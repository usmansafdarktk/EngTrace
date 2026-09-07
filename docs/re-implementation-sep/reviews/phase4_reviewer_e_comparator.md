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
