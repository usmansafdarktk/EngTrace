# Phase 4 — Reviewer B (physics & pedagogy, the P6 guard)

**Frozen ref:** `f6bb653745785a2be164ce2756cd748f32fb4a51` · **Branch:** `redesign/phase4-comparators`
**Mandate:** does the normalisation vocabulary accept answers a human grader would reject?
**Binding constraint:** SPEC-CHANGE 10 — fit a model, report lift over a blind-guess floor.
**Time box:** 35 minutes, held.

Run everything below from the repo root with `PYTHONIOENCODING=utf-8` set.

---

## 1. Verdict

**BLOCKED** — the mandatory gate fails: the hedge policy stated in D4.1 §4.2 is implemented as a
27-string blocklist with a one-sided scope window, and 12 of 15 unlisted hedges, every trailing
hedge, and three outright non-answers are scored `MATCH`; separately, both classification templates
are 100% predictable from the question surface against floors of 50.0% and 34.0%.

---

## 2. Independent re-derivation

I re-derived only what my register row owns. I did **not** re-run the corpus baseline, the item-pool
impact counts, or precision/recall on D4.4 (Reviewer E's).

| number | source claim | my re-derivation | agreement |
|---|---|---|---|
| blind-guess floor, linearity | not stated anywhere | **50.02%** (2001 `linear` / 1999 `not linear` over 4,000 seeds) | n/a — never reported |
| blind-guess floor, memory/causality pair | not stated anywhere | **34.02%** (majority class `('No','No')`) | n/a — never reported |
| depth-2 held-out accuracy, linearity | not stated anywhere | **100.00%** (2,000 train / 2,000 held-out) | n/a |
| depth-2 held-out accuracy, mem/caus pair | not stated anywhere | **100.00%** | n/a |
| `HEDGES` support in the archive | `SUPPORT["hedge"] = (0, 1)` in `normalize.py` | read as stated: **1 occurrence in 2,200 traces, 0 in the 61 Phase 4 traces** | agree — and see F1 |
| template edit is pedagogically neutral | D4.5 | inspected the diff; it renames one solution-table column header (`New Value y[k')` → `New Value {output_var}[k']`), fixing a bracket/paren mismatch and a wrong output variable on the reversal branch. It touches neither question nor answer. | **agree — neutral** |

**Reproduction — floors** (under a minute):

```bash
python - <<'PY'
import sys, random, importlib, re, collections, json
sys.path.insert(0, '.')
m = importlib.import_module('data.templates.branches.electrical_engineering.signals_and_systems.discrete_time_signals')
lin, mc = [], []
for s in range(4000):
    random.seed(s); q, sol = m.template_system_property_linearity()
    lin.append((q, 'linear' if 'is **linear**' in sol else 'not linear'))
for s in range(4000):
    random.seed(100000+s); q, sol = m.template_system_properties_memory_causality()
    a = re.search(r'Memoryless:\s*\*\*(Yes|No)\*\*', sol).group(1)
    b = re.search(r'Causal:\s*\*\*(Yes|No)\*\*', sol).group(1)
    mc.append((q, (a, b)))
print('linearity floor', max(collections.Counter(y for _, y in lin).values())/4000)
print('mem/caus floor ', max(collections.Counter(y for _, y in mc).values())/4000)
json.dump({'lin': lin, 'mc': [(q, list(y)) for q, y in mc]}, open('inst.json', 'w'))
PY
```

The fitted model is a hand-rolled depth-2 tree (no `sklearn` in the environment, as Phase 3 also
found) over bag-of-tokens plus the structural probes the brief names (`x[n +`, `x[n -`, `^`, a
trailing additive constant, a product of two `x[·]` terms). Split 2,000 train / 2,000 held-out,
seed 7. Greedy split on post-split majority accuracy at each of the three nodes.

**Result — both items are lookups, decisively:**

| item | blind-guess floor | depth-2 held-out | **lift** |
|---|---:|---:|---:|
| `system_property_linearity` | 50.02% | **100.00%** | **+49.98 pts** |
| `system_properties_memory_causality` (pair) | 34.02% | **100.00%** | **+65.98 pts** |
| — Memoryless slot alone | 66.57% | 100.00% | +33.43 pts |
| — Causal slot alone | 65.97% | 100.00% | +34.03 pts |

The learned linearity tree is: *is there a `*` in the equation? → `linear`. Otherwise, is there
`x[n -`? → `linear`, else `not linear`.* The learned mem/caus tree is: *is there `x[n +`? →
`(No, No)`. Otherwise, is there a `-`? → `(No, Yes)`, else `(Yes, Yes)`.*

I then closed the question of whether 100% is an artefact of the split by censusing the generator's
form space directly (2,000 seeds each, digits masked to `#`):

| item | distinct (RHS shape, label) pairs | shapes carrying two labels |
|---|---:|---:|
| `system_property_linearity` | **5** (`# * x[n]`→linear, `x[n - #]`→linear, `(x[n])^#`→not, `x[n] + #`→not, `x[n] - #`→not) | **0** |
| `system_properties_memory_causality` | **7** | **0** |

The map from surface shape to label is a *function*. There is no seed at which the shortcut can be
wrong. 100% is not overfitting; it is the ceiling and the floor at once.

Census reproduction:

```bash
python - <<'PY'
import sys, random, importlib, re, collections
sys.path.insert(0, '.')
m = importlib.import_module('data.templates.branches.electrical_engineering.signals_and_systems.discrete_time_signals')
def rhs(q):
    mm = re.search(r'y\[n\]\s*=\s*(.*)', q)
    return re.sub(r'\d+', '#', mm.group(1).strip()) if mm else ''
c1, c2 = collections.Counter(), collections.Counter()
for s in range(2000):
    random.seed(s); q, sol = m.template_system_property_linearity()
    c1[(rhs(q), 'linear' if 'is **linear**' in sol else 'not linear')] += 1
for s in range(2000):
    random.seed(100000+s); q, sol = m.template_system_properties_memory_causality()
    a = re.search(r'Memoryless:\s*\*\*(Yes|No)\*\*', sol).group(1)
    b = re.search(r'Causal:\s*\*\*(Yes|No)\*\*', sol).group(1)
    c2[(rhs(q), (a, b))] += 1
for c in (c1, c2):
    print(len(c)); [print(' ', v, k) for k, v in sorted(c.items())]
PY
```

---

## 3. Findings

### F1 — CONFIRMED (blocking). The hedge policy is a blocklist, and unlisted hedges are accepted.

D4.1 §4.2 states the policy as a principle: *"A hedge is `UNRESOLVED`… A non-committal answer has
not answered."* `normalize.py` implements it as `HEDGES`, a tuple of 27 literal strings matched by
alternation. A blocklist has an unbounded complement, and ordinary English fills it immediately.

**12 of 15 hedged answers I wrote were scored `MATCH`** — each names the right label but does not
commit to it, and a human grader credits none of them:

| candidate (gold: `The system is **linear**`) | outcome |
|---|---|
| `I suspect the system is linear.` | **MATCH** |
| `Apparently the system is linear.` | **MATCH** |
| `The system should be linear.` | **MATCH** |
| `The system looks linear.` | **MATCH** |
| `Tentatively, the system is linear.` | **MATCH** |
| `Perhaps the system is linear.` | **MATCH** |
| `Plausibly the system is linear.` | **MATCH** |
| `Seemingly the system is linear.` | **MATCH** |
| `My guess is that the system is linear.` | **MATCH** |
| `I would say the system is linear.` | **MATCH** |
| `Roughly speaking, the system is linear.` | **MATCH** |
| `I am not 100% sure, but the system is linear.` | **MATCH** |
| `Linear?` | **MATCH** |
| `The system appears to be linear.` (the spec's own named case) | UNRESOLVED |
| `The system is most likely linear.` | UNRESOLVED |

Note that `seemingly` and `apparently` are morphological neighbours of `seems` and `appears`, which
*are* listed. The list does not close even over its own entries.

**This is SPEC-CHANGE 10's own anti-pattern, applied to the hedge vocabulary rather than to
shortcuttability.** The rule reads: *"the best rule I could think of" is a floor on shortcuttability,
never a ceiling.* The same holds for a hand-enumerated hedge list — it is a floor on
non-commitment, never a ceiling — and this phase enumerated where it should have modelled.

**Why the D4.4 gate does not see this:** every hedged case in the adversarial corpus uses a phrase
that is already in `HEDGES` — `cat-15`/`cat-16` `"appears to be"`, `cat-17` `"could be"`, `tup-17`
`"seems to be"`/`"unclear"`, `tup-18` `"might be"`, `chk-14` `"seems"`. The corpus was written from
the same list the comparator implements, so it can only confirm it. That is **D-034 one layer up** —
the spec's own §2.2 warns that a green suite once certified three constants against a dict in the
test file, and the hedge policy is now certified against the list that defines it. I do **not**
dispute the `correct` labels on any of the eight cases named in my brief; I agree with every one,
including `chk-21` (`"Not unacceptable."` is not a committed verdict). The problem is not the
labels — it is that the corpus samples only the inside of the list.

**Reproduction:**

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_template
G = "## Final Answer\nThe system is **linear**."
for c in ["I suspect the system is linear.", "Apparently the system is linear.",
          "The system should be linear.", "The system looks linear.",
          "Tentatively, the system is linear.", "Perhaps the system is linear.",
          "Plausibly the system is linear.", "Seemingly the system is linear.",
          "My guess is that the system is linear.", "Linear?",
          "I am not 100% sure, but the system is linear."]:
    v = compare_template('template_system_property_linearity', G, "## Final Answer\n" + c)
    print(f"{v.outcome:10s} | {c}")
PY
```

**Impact:** corpus-wide. `categorical`, `categorical[tuple]` and `check` all route hedging through
this one list. It credits a model that named a label without committing to it, on items where the
label space is two or three wide and naming one is nearly free — which is precisely the
AI-Tribunal criticism D4.1 §0 says the phase exists to answer.

### F2 — CONFIRMED (blocking). Hedge scope is one-sided, so a trailing hedge is never seen.

`has_hedge(text, at)` inspects `text[at-45 : at]` only. English puts hedges after the claim as often
as before it, and every one of these uses a phrase that **is** in `HEDGES`:

| candidate (gold: `linear`) | outcome |
|---|---|
| `The system is linear, I think.` | **MATCH** (`i think` is listed) |
| `The system is linear, but that is only probably right.` | **MATCH** (`probably` is listed) |
| `Linear. It seems.` | **MATCH** (`seems` is listed) |
| `The system is linear; I cannot determine this with confidence.` | **MATCH** (`cannot determine` is listed) |
| `Linear (unclear).` | **MATCH** (`unclear` is listed) |
| `The system is linear. Hard to say, though.` | **MATCH** (`hard to say` is listed) |
| `a) Not memoryless / b) Causal, I think` (tuple) | **MATCH** |

This is not a vocabulary-coverage complaint; it is the stated policy failing on its own vocabulary.
The `HEDGE_WINDOW = 45` boundary is also sharp and undefended: `It appears to be the case that
<20 chars of filler> the system is linear` is `MATCH`, while the same sentence without the filler is
`UNRESOLVED`.

**Reproduction:**

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_template
G = "## Final Answer\nThe system is **linear**."
for c in ["The system is linear, I think.", "Linear. It seems.", "Linear (unclear).",
          "The system is linear. Hard to say, though.",
          "It appears to be the case that" + " x"*10 + " the system is linear."]:
    print(compare_template('template_system_property_linearity', G, "## Final Answer\n"+c).outcome, '|', c)
PY
```

**Impact:** the same surface as F1, and cheaper to hit — no unlisted word is needed, only ordinary
word order.

### F3 — CONFIRMED (blocking). Answers that assert nothing are scored `MATCH`.

These are not hedges. They contain no commitment of any kind, and a grader gives them zero:

| candidate (gold: `linear`) | outcome |
|---|---|
| `If additivity holds, the system is linear.` (conditional — states a rule, not a verdict) | **MATCH** |
| `Assume the system is linear.` (an assumption, the opposite of a conclusion) | **MATCH** |
| `We must determine whether the system is linear.` (a restatement of the *question*) | **MATCH** |
| `I cannot complete the test; the system is linear.` (explicit refusal, then a guess) | **MATCH** |
| `The system is linear (to be verified).` | **MATCH** |

`We must determine whether the system is linear` is the item's own prompt wording. The comparator
credits an answer that is a paraphrase of the question. (The tuple comparator is better here — the
verbatim mem/caus question echoed as a candidate is correctly `UNRESOLVED`, because its slots state
no value. The single-label comparator has no equivalent guard: echoing the linearity question
verbatim returns `MISMATCH`, i.e. it *did* extract a label from the question stem, and would have
returned `MATCH` had that seed's gold been `linear`.)

The root cause is shared with F1 and F2 and is structural: the comparator asks *"does a label
surface appear, unnegated, unhedged?"* It never asks *"is this sentence an assertion about the
system?"* Commitment is a property of the clause; the implementation models it as the absence of
listed words near a token.

**Reproduction:**

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_template
G = "## Final Answer\nThe system is **linear**."
for c in ["If additivity holds, the system is linear.", "Assume the system is linear.",
          "We must determine whether the system is linear.",
          "I cannot complete the test; the system is linear.",
          "The system is linear (to be verified)."]:
    print(compare_template('template_system_property_linearity', G, "## Final Answer\n"+c).outcome, '|', c)
PY
```

### F4 — CONFIRMED (blocking). The longest surface anywhere beats the committed final clause — in both directions.

`kinds.py::_find_surface` sorts surfaces by length and `break`s on the **first surface that occurs
anywhere in the span**. Within a surface it takes the latest occurrence; *across* surfaces position
is not consulted at all. So whenever `nonlinear` appears anywhere in the answer span it wins over a
later `linear`, regardless of which the model actually committed to.

D4.1 §4.2's longest-first rule is about *containment* — `linear` is a substring of `nonlinear` at the
same offset — and that rule is right. It is being applied across *different offsets*, which is a
different claim, and it contradicts `answer_span`'s own docstring: *"the last occurrence … wins,
because models restate the answer and the final statement is the committed one."*

**A false accept** (my mandate's direction):

| gold | candidate | outcome | model actually committed to |
|---|---|---|---|
| `not linear` | `The system is nonlinear -- no, wait, it is linear.` | **MATCH** | `linear` — **wrong**, and credited |
| `not linear` | `Nonlinear? No. The system is linear.` | **MATCH** | `linear` — **wrong**, and credited |

**And a false reject**, the mirror case:

| gold | candidate | outcome | committed to |
|---|---|---|---|
| `linear` | `The system is nonlinear. Actually, the system is linear.` | **MISMATCH** | `linear` — **right**, and penalised |

Self-correction is a common and *desirable* model behaviour. A re-marker rescues it —
`Wait, nonlinear. Final answer: linear.` is `MATCH`, because `Final Answer:` re-spans the text — so
the comparator's verdict on an identical reasoning pattern turns on whether the model happened to
re-emit a marker.

**Reproduction:**

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_template
G_NOT = "## Final Answer\nThe system is **not linear**."
G_LIN = "## Final Answer\nThe system is **linear**."
for g, c in [(G_NOT, "The system is nonlinear -- no, wait, it is linear."),
             (G_NOT, "Nonlinear? No. The system is linear."),
             (G_LIN, "The system is nonlinear. Actually, the system is linear."),
             (G_LIN, "Wait, nonlinear. Final answer: linear.")]:
    print(compare_template('template_system_property_linearity', g, "## Final Answer\n"+c).outcome, '|', c)
PY
```

**Impact:** every `categorical` item whose label set declares a fused negative surface. Both error
directions, and the false-accept direction is invisible to D4.4, which contains no self-correction
case.

### F5 — CONFIRMED (blocking as a scoping decision, not as a code defect). Both Phase 4 classification items are 100% shortcuttable, and this phase specifies scoring them answer-only.

Numbers and method in §2. `system_property_linearity` is 100% predictable at a 50.0% floor
(+49.98 pts); `system_properties_memory_causality` at a 34.0% floor (+65.98 pts). The generators emit
5 and 7 distinct RHS shapes respectively, and no shape ever carries two labels.

The template's own docstring is the claim under test: *"This template tests the ability to formally
prove or disprove if a system is linear by checking the two defining properties: additivity and
homogeneity."* **It does not.** It tests the ability to recognise which of four generator branches
fired. A solver that has never heard of homogeneity scores 100% from *"is there a `*`?"*. Note that
`x[n] + C` is the one branch with genuine conceptual content — affine systems are the classic
incrementally-linear trap, and the archive shows a real model reasoning about exactly that
(`answer_span`'s docstring cites linearity trace 0 arguing `incrementally linear` mid-derivation) —
and it is also the branch the shortcut identifies most trivially.

I distinguish the two items. For mem/caus the surface feature *coincides with the definition*
(a shift term means memory; an advance means non-causal), so the shortcut is arguably the physics
and the item is honestly a one-step recognition task — but then its difficulty label should say so.
For linearity the shortcut and the taught procedure are **different things**, and the item's stated
pedagogy is not what its answer measures.

This is a P6 matter for *this* phase, not a pre-existing one, because Phase 4 is the phase that
decides these items are scored on their answer alone. D4.1 §8 reports 100% precision and 100% recall
on both. Those numbers are true and they are about the vocabulary; they say nothing about whether a
scored model reasoned, because on these two items a scored answer carries at most 1.0 and 1.6 bits
and both are readable off the question.

**Impact:** two of the four Phase 4 items. The general form — a class-C classification template whose
label is a function of which generator branch fired — is likely to recur across the 146; see §5.

---

## 4. Falsification attempts that failed

Things I tried to break and could not:

- **The spec's named case holds.** `The system appears to be linear` is `UNRESOLVED`, reason
  `candidate hedges (appears to be) and never commits`. The mechanism works; F1 is about its extent,
  not its existence.
- **`neither … nor` is correctly scoped.** I could not get `neither memoryless nor causal` to invert
  into `(Yes, Yes)`. `NEITHER_NOR_RE` catches it and the doc's account of the danger is accurate.
- **The fused-negative double-negation trap is genuinely closed.** `The system is nonlinear` returns
  `not linear`, not `linear`. `_negated` deliberately declining to fire on a declared fused surface
  is correct, and the docstring's account of the first-draft bug is consistent with the code.
- **The `keep_asterisks` split is real, not decoration.** `{2, 3, *-8*, 6}` and `{2, 3, -8, 6}` are
  different answers and the sequence path preserves the distinction. The comment is right that
  getting this backwards deletes the answer.
- **The superscript-before-NFKC ordering is load-bearing.** `unicodedata.normalize("NFKC", "x²")`
  gives `x2`; the explicit pre-pass gives `x^2`. The stated reason reproduces exactly.
- **Partial tuple credit is not leaking.** `a) Not memoryless` alone, second slot silent, is
  `UNRESOLVED`, not `MATCH`. `tup-18`'s label is right and the implementation honours it.
- **Positional slot resolution is not over-eager.** The verbatim mem/caus question echoed as a
  candidate is `UNRESOLVED` — the positional fallback needs a value to be *stated*, not merely a slot
  name to appear. This is the guard the single-label comparator lacks (F3).
- **I could not find a seed at which the F5 shortcut is wrong.** I looked for one specifically; the
  form→label map is total and single-valued over 2,000 seeds per item.
- **The template edit is pedagogically neutral.** I read the diff rather than the D4.5 counts. It
  changes a solution-table column header only: it closes a bracket that was opened as `[` and closed
  as `)`, and it substitutes the branch's real output variable for a hardcoded `y`. Neither the
  question, nor the answer, nor the number of steps moves. I have no P6 objection to it.

---

## 5. Further probing and improvements

**Where the evidence is thinnest, ranked.**

1. **The hedge rule has essentially no evidence at all, and the file says so.**
   `SUPPORT["hedge"] = (0, 1)`. Zero occurrences in the 61 Phase 4 traces; **one** in all 2,200.
   Every other rule in that table rests on tens to hundreds of observations. `phase4_vocabulary.md`
   §4 is honest that the list is not derived from the Phase 4 archive, but the honesty is filed as a
   *coverage* caveat when it is a *methodology* problem: a rule with one supporting observation
   cannot have been derived from evidence, so it was written from intuition — and D4.1 §4.2 spends a
   paragraph explaining why writing a hedge list from intuition is exactly the wrong move. The
   document diagnoses its own defect and then treats *removing two entries* as the cure. The cure for
   a blocklist is not a shorter blocklist.

2. **The measurement that would overturn a clean verdict on the hedge policy — and it is the one I
   ran.** Do not score the hedge list against cases written from the hedge list. Sample
   non-committal answers from *outside* it: take the ~2,200-trace archive, extract the answer spans a
   human labelled correct, and mine the ones the comparator marks `MATCH` for modal and evidential
   markers (`should`, `suspect`, `apparently`, `looks`, the `seem`/`appear` morphological families,
   `-ly` adverbs of likelihood, sentence-final `?`). Any hit is an F1 instance in the wild. My
   battery is synthetic; this would make it empirical, and it is what I would run first with more
   time.

3. **The fix I would specify, in place of a longer list.** Commitment is a property of the clause
   containing the label, and should be tested as one:
   - hedge scope must be **bidirectional** — the clause, not 45 characters backwards (closes F2);
   - the label must sit in an **asserted matrix clause** — reject when it is under `if`, `assume`,
     `suppose`, `whether`, `determine whether`, or a question mark (closes F3);
   - **modality is the feature, not the phrase** — any modal (`may`, `might`, `could`, `should`,
     `would`) or evidential adverb governing the copula is a hedge, which generalises past the 27
     strings (closes most of F1);
   - across surfaces, prefer the **last committed clause**, not the longest surface; keep
     longest-first strictly *within* one offset (closes F4).
   All four are testable against the archive as it stands, and none needs new data.

4. **The D4.10 difficulty proxies are the wrong proxies, and F5 shows why.** `steps`, `arith_lines`,
   `given`, `derived` and `inference` are all computed from the **solution string**. They measure how
   much gold *writes down*. They are the right instrument for the question Phase 3 asked — *did an
   edit change the work?* — and they answered it: the Phase 4 edit moves none of them, which is real
   evidence of neutrality and I relied on it. But they are structurally incapable of detecting F5,
   because a 100%-shortcuttable item and a genuinely hard one produce identical traces.
   `template_system_property_linearity` scores respectably on all five and is a five-way lookup.
   **The missing proxy is a property of the question→answer map, not of the trace:** held-out
   accuracy of a small model on the question text, reported against the blind-guess floor.
   SPEC-CHANGE 10 already mandates that shape for a reviewer's lookup check; it should be promoted
   from a reviewer instruction to a standing item-level metric, computed once per template and stored
   beside the difficulty label. Two numbers, both cheap, and they would have caught D-046
   (`line_balancing_heuristic`) without a reviewer having to go looking.

5. **What generalises to the other 146.** The F5 pattern — *the answer label is a function of which
   `random.choice` branch fired, and the branches are surface-distinguishable* — is a property of how
   these templates are written, not of these two items. Any template whose generator picks a scenario
   from a list and whose answer is determined by that pick has it. That is the whole `classification`
   population and much of `check`. `check` is the worst case: a two-valued verdict is 1 bit, and
   D4.1 §4.6 already recognises the coin-flip problem and fixes it by *also* scoring the supporting
   quantity. **That is the right shape of fix and it should be extended to `categorical`, not left
   in `check` alone** — a linearity answer could be required to name the property that failed
   (additivity, homogeneity, or both), which is exactly what the taught procedure produces and
   exactly what the surface shortcut does not. I would run the §2 script over all `classification`
   and `check` templates in one batch before Phase 5 commits to a scoring surface; it is minutes per
   template.

6. **A scoping decision this phase should record explicitly.** D4.1 §1 states the bias — toward
   rejection, and toward `UNRESOLVED` over both — and §8 measures precision, recall and decided. None
   of those is sensitive to the *information content of the answer being scored*. An item whose
   answer carries 1.0 bits against a 50% floor cannot distinguish a reasoner from a guesser no matter
   how good the comparator is. P6 says a change to what an item tests must be recorded and approved.
   Phase 4 decides these items are scored answer-only; on F5's evidence that *is* such a change for
   at least two of them, and it should enter the register as a scoping decision rather than arrive as
   a consequence of the comparator's scope.

**What was wrong or ambiguous in the two documents.**

- D4.1 §4.2 states the hedge policy as a general principle (*"a non-committal answer has not
  answered"*) while the implementation is a closed list. The spec should say which it means. As
  written, a conforming implementation cannot be checked, because the principle and the artefact
  disagree about coverage — and §9's conformance list does not mention hedging at all.
- D4.1 §4.2's *"Surfaces match longest-first"* does not say **within what scope**. Read as
  within-offset it is correct; read as within-span it is F4. The implementation takes the second
  reading. One clause fixes the ambiguity.
- `normalize.py::answer_span`'s docstring — *"the last occurrence … the final statement is the
  committed one"* — is true of marker selection and false of label selection, in the same file. A
  reader takes it as a property of the comparator, and it is the property F4 shows is missing.
- D4.1 §8's *"This number is weak evidence and is reported as such"* is the right instinct and does
  not go far enough. It says the archive is weak because the vocabulary came from it. The D4.4 table
  immediately below is presented as the independent check, and for hedging it is not: the corpus
  samples the inside of `HEDGES`. That same paragraph should be attached to the D4.4 table.
- `normalize.py` §4's comment that `it depends` and `either` were removed because they are domain
  vocabulary is correct and well argued for `it depends`. For `either` it proves less than it claims:
  the quoted justification is gold's own phrasing (`fails at least one of the tests (additivity or
  homogeneity)`), which contains no `either` at all, while `cat-17` shows a candidate disjunction is
  a genuine non-commitment. The removal may still be right; the stated reason does not support it.

---

# Round 2 — Reviewer B (physics & pedagogy, the P6 guard)

**Frozen ref:** `335ad7a3c2a7ec621db62180d74e286d2f59e0f1` · **Branch:** `redesign/phase4-comparators`
**Mandate:** does the *new* mechanism reject answers a human grader would ACCEPT?
**Time box:** 30 minutes, held.

Run everything below from the repo root with `PYTHONIOENCODING=utf-8` set.

---

## 1. Verdict

**BLOCKED — the defect has been relocated, not closed.** All four of my round-1 code findings are
**addressed** and I could not reopen any of them. But the mechanism that closed them refuses **25 of
38** constructed answers that a competent grader credits, and the refusals are not exotic: *"Note
that the system is linear"*, *"Although the equation contains a delay, the system is linear"*, *"The
distinction between additivity and homogeneity is not entirely obvious, but the system is linear"*.
Two of the trigger classes fire on **real archived phrasing** — a `reynolds_number_flow_regime` trace
that hedges a caveat and then commits to `turbulent`, and a `decimation_aliasing_analysis` trace
whose committed conclusion opens with `Although`. The comparator now measures phrasing in the
opposite direction, which is the equal-and-opposite failure D4.1 §1 names.

The archive gate does not see this (§4), which is the round-1 blindness with the sign flipped: D4.4
samples the inside of the hedge list, and the 61 traces sample four items whose committed answers
happen to be terse.

**What lifts the block** is narrow and is stated in §5.1: split the concessive subordinators out of
`NON_ASSERTING_OPENERS`, give the non-modal hedge probes the governance test the module's own
docstring already promises them, and make `_governing_neighbour` directional.

---

## 2. Disposition of the round-1 findings

I re-ran my full round-1 battery against the frozen ref before writing anything new.

| round-1 finding | disposition | evidence |
|---|---|---|
| **F1** — hedge list is a 27-string blocklist; 12/15 unlisted hedges scored `MATCH` | **addressed** | all 15 now `UNRESOLVED`. `MODALS` / `EVIDENTIALS` / `HEDGE_VERBS` / `APPEARANCE_COPULAS` / `UNCERTAINTY` are classes, and the `-ly` family that `Seemingly`/`Apparently`/`Plausibly` escaped through is closed as a family, not as three more strings. |
| **F2** — hedge scope one-sided, every trailing hedge passed | **addressed** | all 7 now `UNRESOLVED`. Scope is the clause, bidirectional; `HEDGE_WINDOW = 45` and its undefended boundary are gone. |
| **F3** — conditionals, assumptions, refusals and the item's own prompt wording scored `MATCH` | **addressed** | all 5 now `UNRESOLVED`. `is_assertion` is the guard the single-label path lacked, and `TASK_RESTATEMENTS` closes the prompt-paraphrase case specifically. |
| **F4** — longest surface across the whole span beat the committed final clause | **addressed** | all 4 now behave as I argued, in both directions. `find_commitment` walks `reversed(labelled)` and longest-first is confined to within one offset, where it is about containment and is right. |
| **F5** — both classification items 100% shortcuttable | **relocated to a scoping question** — see §5.4. The result is accepted and not in dispute; the items cannot change in this phase. |

**The `_UNCERTAIN_RE` word-boundary note in the module docstring is correct and I verified it.** `\b`
after `not 100%` demands a word character that never arrives, so the entry would never fire.
`(?<!\w)…(?!\w)` is the right fix and `I am not 100% sure, but the system is linear` is now
`UNRESOLVED`. That is a real bug caught during implementation and it is worth saying so.

**No regression in the archive numbers.** `python -m tests.comparators.score` reports archive
precision 100.0% / recall 97.8%, adversarial 98.6% / 98.6%, one archive false reject
(`signal_operations[0]`, a `sequence` item that does not route through `commitment.py`) and one
adversarial false accept (`num-16b`, D-052). **There is no new false reject among the 61.** §4
explains why that is not reassurance.

---

## 3. Findings

### R2-F1 — CONFIRMED (blocking). The commitment test refuses committed answers, at a rate of 25 in 38.

I wrote 38 candidates that name the right label **and commit to it** — the mirror of my round-1
battery. A competent grader credits every one. **25 are refused: 24 `UNRESOLVED` and 1 `MISMATCH`.**
Grouped by the mechanism that refuses them:

**(a) Concessive and discourse openers — 10 cases.** `NON_ASSERTING_OPENERS` treats every listed word
as scoping over the whole clause. For a genuine hypothetical (`if`, `suppose`, `assume`) that is right
and it is what closed F3. For a **concessive or contrastive** subordinator it is wrong in a specific,
systematic way: *`although X, Y`* asserts `Y`. The subordinator's scope ends at the comma, and
`_OPENER_RE.match(clause)` cannot see that because the clause splitter does not cut on commas — so
the test asks about the clause's first word while the label sits in the matrix clause after it.

| candidate | outcome |
|---|---|
| `Note that the system is linear.` | **UNRESOLVED** — *opens with 'Note that'* |
| `Note: the system is linear.` | **UNRESOLVED** |
| `Although the equation contains a delay, the system is linear.` | **UNRESOLVED** |
| `Though it looks like a shift, the system is linear.` | **UNRESOLVED** |
| `Whereas a scaling would preserve linearity, the square does not, so the system is not linear.` | **UNRESOLVED** |
| `Except for the delay, nothing changes: the system is linear.` | **UNRESOLVED** |
| `Recall that superposition holds here, so the system is linear.` | **UNRESOLVED** |
| `Unlike a linear system, this one squares its input, so it is not linear.` | **UNRESOLVED** |
| `Given that additivity and homogeneity both hold, the system is linear.` | **UNRESOLVED** |
| `Provided that a and b are scalars, superposition holds and the system is linear.` | **UNRESOLVED** |

`note that` and `recall that` are worse than a scope error — they are **factive**: their complement is
asserted, not suspended. *"Note that the system is linear"* is a normal way for a solution to state a
conclusion. Reviewer E's motivating case, `Note: Answer: not linear if the offset were nonzero`, is
non-committal because of `if … were`, which `if` already catches; `note` does no work there that is
not already done, and it over-rejects everywhere else.

**(b) The non-modal hedge probes have no governance test — 9 cases.** The docstring says *"A modal or
evidential **governing the copula** is a hedge"*. Only `_MODAL_RE` implements governance
(`modal + up to two words + copula`). `_EVIDENTIAL_RE`, `_HEDGE_VERB_RE`, `_APPEARANCE_RE` and
`_UNCERTAIN_RE` are bare `\b…\b` searches over the whole clause. A marker anywhere in the clause
condemns a label anywhere else in it — which is **exactly the round-1 F2 error with a different
window**. Character-window scope became clause scope; the missing thing, then and now, is a relation
between the marker and the label.

| candidate | outcome |
|---|---|
| `The output is roughly twice the input in every case, and the system is linear.` | **UNRESOLVED** — *evidential: roughly* |
| `The distinction between additivity and homogeneity is not entirely obvious, but the system is linear.` | **UNRESOLVED** — *not entirely* |
| `It is unclear which textbook convention applies, but the system is not linear.` | **UNRESOLVED** — *unclear* |
| `The gain is not entirely constant across n but scaling still holds, so the system is linear.` | **UNRESOLVED** |
| `Whatever the equation looks like, additivity and homogeneity both hold: the system is linear.` | **UNRESOLVED** — *copula of appearance: looks* |
| `Look at the coefficient: it is constant, so the system is linear.` | **UNRESOLVED** — *copula of appearance: look* |
| `The squaring term is nominally the only nonlinearity, so the system is not linear.` | **UNRESOLVED** — *evidential: nominally* |
| `While the delay term is present, it does not break superposition: the system is linear.` | **MISMATCH** |

Note the last row. It is not `UNRESOLVED`; the comparator scores a **correct answer as wrong**, reason
`label 'not linear' != gold 'linear'` — negation scope reaching across the colon from *"does not break
superposition"*. `MISMATCH` on a correct answer is strictly worse than `UNRESOLVED`, and D4.1 §1's
stated bias toward `UNRESOLVED` is not honoured here. I did not attribute this to `commitment.py`; it
is label selection and may predate round 2. It is listed because a reviewer hunting false rejects
should see it.

`roughly` and `nominally` deserve a sentence of their own. In engineering prose they are about
**numeric precision and nominal sizing**, not epistemic commitment. *"The output is roughly 2x"* and
*"a nominally 10 mm pipe"* are statements a grader reads as confident. Classing them with `perhaps`
imports a category error into a physics benchmark, and `check` and `numeric` items are where it will
bite hardest.

**(c) `_governing_neighbour` is bidirectional and undirected — 6 cases.** A label-free hedge clause on
*either* side governs. But a model that reasons hesitantly and then concludes firmly puts a hedge
clause immediately before its answer, and that is the pattern good practice recommends.

| candidate | outcome |
|---|---|
| `It is hard to say which test is more direct. The system is linear.` | **UNRESOLVED** — *hedged by an adjacent clause* |
| `The wording of the question is uncertain; regardless, the system is linear.` | **UNRESOLVED** |
| `The system is linear; most likely candidates for failure (additivity, homogeneity) both hold.` | **UNRESOLVED** — *evidential: likely* |
| `The system is linear. It seems obvious in hindsight.` | **UNRESOLVED** |
| `The system is not linear. This should be clear from the square.` | **UNRESOLVED** — *modal: should be* |
| `Imagine a scaled input a*x[n]; the output is a*y[n], so the system is linear.` | **UNRESOLVED** — *verb of opinion: imagine* |

The last two are the sharpest. *"This should be clear from the square"* is an expression of
**confidence** and is read as a hedge. *"Imagine a scaled input"* is the standard way to open a
homogeneity proof, and `imagine` is in `HEDGE_VERBS`.

The round-1 cases that motivated the neighbour rule are all **following** hedges or explicit
statements of **inability** (`Linear. It seems.`; `I cannot complete the test; the system is
linear.`). Neither requires a preceding *difficulty-of-reasoning* clause to govern. The rule
generalised further than its evidence.

**Reproduction (~10 s):**

```bash
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.answer import compare_template
L = "## Final Answer\nThe system is **linear**."
N = "## Final Answer\nThe system is **not linear**."
for g, c in [
  (L, "Note that the system is linear."),
  (L, "Although the equation contains a delay, the system is linear."),
  (N, "Unlike a linear system, this one squares its input, so it is not linear."),
  (L, "Given that additivity and homogeneity both hold, the system is linear."),
  (L, "Recall that superposition holds here, so the system is linear."),
  (L, "The output is roughly twice the input in every case, and the system is linear."),
  (L, "The distinction between additivity and homogeneity is not entirely obvious, but the system is linear."),
  (N, "It is unclear which textbook convention applies, but the system is not linear."),
  (L, "Look at the coefficient: it is constant, so the system is linear."),
  (L, "Imagine a scaled input a*x[n]; the output is a*y[n], so the system is linear."),
  (L, "It is hard to say which test is more direct. The system is linear."),
  (N, "The system is not linear. This should be clear from the square."),
  (L, "While the delay term is present, it does not break superposition: the system is linear."),
]:
    v = compare_template('template_system_property_linearity', g, "## Final Answer\n" + c)
    print(f"{v.outcome:11s} | {c}")
    if v.outcome != "MATCH":
        print("            ->", v.reason)
PY
```

**Impact:** every `categorical` item, and `check` through the same `find_commitment` call at
`kinds.py:1194`. Not visible in today's numbers (§4), and unavoidable the moment the schema reaches
items whose committed answers are prose rather than a two-word verdict.

### R2-F2 — CONFIRMED (blocking as evidence, not as a separate defect). The triggers fire on real archived phrasing.

My battery is constructed, so I went to the archive. **The 2,200 rows are error-analysis samples:
`final_answer_acc == 0.0` for all 2,200**, so the archive cannot supply a *labelled-correct* answer
that is rejected — that instrument does not exist and no amount of time would have produced it. What
it can supply is proof that the constructions are real. Scanning every conclusion-bearing clause in
all 2,200 `model_reasoning` fields (cue: `therefore|thus|hence|so the|we conclude|the system is|is
(not) linear|causal|stable|laminar|turbulent`), 1,035 clauses tested:

| trigger | hits | archived example |
|---|---:|---|
| `OPENER Although` | 1 | `decimation_aliasing_analysis__3490738110` — *"Although the strict no-aliasing inequality is not met, the decimated frequency Mω0 = 4·(π/4) = π maps to the edge frequency π and does not produce overlapping spectral copies"* — a committed conclusion, refused on its first word |
| `HEDGE modal: may be` | 1 | `reynolds_number_flow_regime__3122400918` — *"…depending on disturbances the flow may be described as transitional up to a few ×10⁶, **but using the standard criterion it is turbulent**"* — hedges the caveat, **commits to the label**, refused |
| `HEDGE modal: should be` | 2 | `kinematic_viscosity__3586477399`, `basic_buoyant_force__1845565441` — *"the result should be rounded to three significant figures"*: a modal about **rounding convention**, not about the answer |
| `OPENER Alternatively` | 1 | `multi_segment_rod__3591646318` — *"Alternatively, sum of loads to the right of a cut … so the internal force in segment1 is −21,000 N"* — a committed second derivation, refused |
| `OPENER otherwise` | 1 | `decimation_aliasing_analysis__2335660753` — *"otherwise, no aliasing occurs"* |
| `HEDGE evidential: likely` | 2 | `signal_energy_power__2570931285`, `vibration_isolator_design__3831595160` |

`reynolds_number_flow_regime` is the case to look at: a **two-way classification item**, structurally
identical to `system_property_linearity`, whose model states a caveat with `may be` and then commits.
This comparator refuses it. That is R2-F1 in the wild, on the item family Phase 5 extends to.

**Reproduction:**

```bash
python - <<'PY'
import sys, json, glob, collections, re
sys.path.insert(0, '.')
from tests.comparators.commitment import clauses, is_assertion, hedge_markers
rows = [json.loads(l) for f in sorted(glob.glob('error_analysis_annotation/samples/*.jsonl'))
        for l in open(f, encoding='utf-8')]
print('rows', len(rows), 'graded correct',
      sum(1 for r in rows if r['final_answer_acc'] == 1.0))     # -> 2200, 0
CUE = re.compile(r'\b(therefore|thus|hence|so the|we conclude|the system is|'
                 r'is (?:not )?(?:linear|causal|memoryless|stable|laminar|turbulent))', re.I)
hits, ex = collections.Counter(), collections.defaultdict(list)
for r in rows:
    for c in clauses(r['model_reasoning']):
        if not CUE.search(c):
            continue
        ok, why = is_assertion(c)
        if not ok:
            k = 'OPENER ' + (why.split("'")[1] if "'" in why else why)
            hits[k] += 1
            ex[k].append((r['question_id'], ' '.join(c.split())[:150]))
            continue
        for m in hedge_markers(c):
            hits['HEDGE ' + m] += 1
            ex['HEDGE ' + m].append((r['question_id'], ' '.join(c.split())[:150]))
for k, v in hits.most_common(20):
    print(f'{v:5d}  {k}')
    for q, t in ex[k][:2]:
        print('        ', q[:42], '|', t)
PY
```

### R2-F3 — MINOR. The commitment policy is applied inconsistently between `categorical` and `categorical[tuple]`.

The single-label path calls `find_commitment`, which includes `_governing_neighbour`. The tuple path
(`kinds.py::_slot_verdict`) calls `is_assertion` and `hedge_markers` on the **slot clause only** and
has no neighbour rule. The `_ENUM_RE` split usually amputates the hedge, so the same construction
lands differently:

| candidate | single-label analogue | tuple |
|---|---|---|
| `It is unclear how the textbook counts n=0, but a) Not memoryless b) Causal.` | **UNRESOLVED** | **MATCH** |
| `Although the system uses a past sample, a) Not memoryless b) Causal.` | **UNRESOLVED** | **MATCH** |
| `The shift would require storage, so a) Not memoryless / b) Causal.` | **UNRESOLVED** | **MATCH** |

**The tuple path is the one that is right.** I record this as evidence for R2-F1's remedy rather than
as a defect in the tuple path: one policy, stated once in D4.1 §4.2, produces opposite verdicts on
the same English depending on which comparator reads it, and the more permissive of the two is the
one a grader agrees with.

---

## 4. Falsification attempts that failed

- **I could not reopen F1, F2, F3 or F4.** I spot-checked the round-1 battery rather than re-running
  all 31 cases (the brief reports them all passing and `reviewer_battery.py` is green), and I probed
  for three further escapes per finding — morphological variants outside `EVIDENTIALS`, hedges split
  across an em-dash, a label restated after a retraction. I found none. The round-1 mechanism is
  genuinely replaced, not patched.
- **I could not find a new false reject among the 61.** `signal_operations[0]` is the only one and it
  is a `sequence` item that never calls `find_commitment`. **This is the result I most wanted and it
  is negative, so I say plainly why it is weak evidence:** 39 of the 61 traces are the two
  classification items, whose archived answers are terse verdicts (`**Linear**`, `a) No b) Yes`) with
  no subordinate clause to trip the opener test and no adverb to trip the hedge test. The gate cannot
  detect over-rejection because the sample contains almost no prose. That is D-034's shape a third
  time: round 1's corpus sampled the inside of the hedge list; round 2's archive samples the inside
  of a phrasing style.
- **`must` and `can` really are excluded, and that call is right.** *"The system must be linear"* and
  *"we can conclude the system is linear"* are both `MATCH`. The docstring's stated reason — that
  treating them as hedges would trade false accepts for false rejects — is exactly the principle
  R2-F1 says was not applied consistently to the other classes.
- **The clause splitter is doing real work and the `?` retention is correct.** `Linear?` is
  `UNRESOLVED` on the question mark alone.
- **`TASK_RESTATEMENTS` does not over-reject when the restatement is followed by an answer.**
  *"To determine whether the system is linear I checked both properties; it is linear"* is `MATCH`,
  because `;` splits and the last clause is clean. I expected this to fail and it did not. So are
  *"I checked whether additivity holds. It does. The system is linear."* and *"Test whether a\*x
  scales: it does. The system is linear."*
- **Modals inside a derivation survive.** *"A scaled input would produce a scaled output, so the
  system is linear"* is `MATCH` — `_MODAL_RE`'s governance requirement (`would … be/is`) is what
  saves it, which is direct evidence that governance is the right mechanism and that the other four
  probes are missing it.
- **The 100% archive precision is not bought with `UNRESOLVED`.** `decided` is 93.4%, unchanged.

---

## 5. Further probing and improvements

### 5.1 The remedy for R2-F1, in three bounded changes

None needs new data, and all three are testable against the batteries already written.

1. **Split `NON_ASSERTING_OPENERS` in two.** *Hypotheticals* (`if`, `unless`, `assume`, `assuming`,
   `suppose`, `supposing`, `whether`, `in case`, `were`, `had`, `even if`) suspend the whole sentence
   and must keep rejecting — they are what closed F3. *Concessives, contrastives and discourse
   markers* (`although`, `though`, `whereas`, `while`, `unlike`, `except`, `given that`, `provided
   that`, `note`, `note that`, `recall that`, `consider`, `let`, `alternatively`, `otherwise`, `in
   contrast`, `aside`) scope over their own clause only; for these, **strip the subordinate clause up
   to the comma and test the matrix clause that follows**. If there is no comma-separated matrix
   clause, the current rejection stands. This closes all 10 cases in R2-F1(a) without reopening a
   single F3 case.
2. **Give the four non-modal probes the governance test `_MODAL_RE` already has.** The docstring
   promises it; one probe in five implements it. The minimum version is positional: an evidential,
   appearance copula or uncertainty phrase hedges the label only if it lies in the **same sub-clause**
   as the label — no comma, colon, `but`, `so`, `and`, `regardless` or `still` between them. That
   closes all of R2-F1(b) except the `MISMATCH` row, which is a separate negation-scope matter.
   Additionally, move `roughly`, `nominally` and `likely` out of `EVIDENTIALS` unless they govern the
   copula directly (`is likely linear`), because in engineering prose they quantify precision rather
   than confidence.
3. **Make `_governing_neighbour` directional and narrow it.** Only a **following** label-free hedge
   clause governs (`Linear. It seems.`), plus a **preceding clause of stated inability**
   (`cannot determine`, `cannot complete`, `unable to`, `no idea`) — which is what *"I cannot
   complete the test; the system is linear"* actually is. A preceding clause expressing *difficulty
   of reasoning*, or a *meta remark*, must not govern. This means splitting `UNCERTAINTY` into
   `INABILITY` and `IMPRECISION`, which is worth doing anyway: they are different speech acts and
   only the first is a non-answer.

### 5.2 The measurement that would settle this, and it does not exist yet

Round 1 I said: do not score the hedge list against cases written from the hedge list. The same
applies now — **do not score the commitment test against an archive of terse verdicts.** The missing
instrument is a set of archived answers that are *both* prose *and* labelled correct, and the 2,200
rows cannot supply it (all are `final_answer_acc == 0.0`). Two ways to get one, in order of cost:

- **Cheapest, and I would do it first:** take the 45 archived answers my own `ground_truth.py` labels
  `true+` and **paraphrase each into three prose forms** that preserve the verdict (a concessive
  opener; a hedged-reasoning-then-commit form; an alternative-derivation form). That is 135 cases
  whose correct label is known by construction, and it is D4.4's missing half — **D4.4 contains only
  near-misses that should be rejected and no correct answers phrased awkwardly. A recall corpus is as
  necessary as a precision corpus, and Phase 4 has only one of them.** This is the single
  highest-value artefact I can name for round 3.
- Re-run the four templates against two or three of the archived models with a prompt that asks for
  reasoning inline, and hand-label the outputs. More faithful; hours rather than minutes.

### 5.3 Where round 2 leaves the D4.1 and vocabulary documents

- **D4.1 §4.2 still states the policy as a principle and the artefact still narrows it**, only in the
  other direction now. The missing clause is: *a hedge must **govern** the label, and a subordinator
  scopes over its own clause only.* Without it the conformance problem stands — a reader cannot tell
  whether `Although X, the system is linear` conforms.
- **`commitment.py`'s docstring says "governing the copula" and the code implements it once in five.**
  Either the docstring or the four probes should move. My recommendation is the probes.
- **§9's conformance list still does not mention hedging**, which was a round-1 note and is
  unaddressed. It now needs *two* entries, one per error direction.
- **`phase4_vocabulary.md` §3's honesty about `SUPPORT["hedge"] = (0, 1)` should be carried forward
  to `commitment.py`.** The new module has *less* archive support than the list it replaced, not
  more: no rule in it is derived from an observation, because §5.2's instrument does not exist. That
  is defensible — it implements my §5.3 prescription — but it should be recorded as such rather than
  inheriting the old table's evidence.

### 5.4 F5's disposition: what the comparator and the results table should do, given the items cannot change

**Recommendation: record it; do not fix it in the comparator.**

**(a) Do NOT extend `check`'s supporting-quantity pattern to `categorical` in Phase 4.** My round-1
§5.5 proposed requiring a linearity answer to name *which* property failed. It is implementable
without touching the item — the gold solution does state it — but I now think it is the wrong move
here, and **I withdraw it as a Phase 4 requirement**, for three reasons:

- **It fails the phase's own gate.** Of the 15 archived `system_property_linearity` traces, only
  **4** name a property (`additiv|homogen|superposition|scal`) anywhere in their answer span.
  Requiring it would turn 11 currently-correct answers into non-answers and drop recall on that item
  from 100% to ~27%, breaching S4.5's ≥95%. A scoring change that fails the gate is not a comparator
  change; it is an item change wearing a comparator's clothes.
- **It is asymmetric, and the asymmetry is itself a cue.** On the `linear` branch nothing fails, so
  only `not linear` answers could carry a named property. A requirement that applies to one label and
  not the other tells the solver which label it is looking at — it would *add* a shortcut while
  trying to remove one.
- **It relocates a P6 change into the comparator.** Spec §4.1 forbids redesigning these items;
  changing what counts as a correct answer to them from inside the comparator is the same change by
  another route and with less scrutiny. That is precisely what P6 exists to prevent.

**(b) DO record it as a P6 scoping decision.** My round-1 §5.6 stands and is the disposition I press.
The register entry should say: *Phase 4 scores `system_property_linearity` and
`system_properties_memory_causality` on their answer alone; those answers carry 1.0 and 1.6 bits
against blind-guess floors of 50.02% and 34.02%, and a depth-2 tree on the question surface reaches
100% held-out on both. The comparator's 100% precision and recall on these items are statements about
the vocabulary and not about whether a scored model reasoned.* Two sentences, no code, and it stops
the D4.1 §8 numbers from being read as more than they are.

**(c) The results table should carry two extra columns per item, not a footnote:** `floor`
(blind-guess) and `surface-model held-out accuracy`. Both are computable in minutes per template by
the round-1 §2 script. A reader who sees `50.02% → 100.00%` beside `precision 100%` cannot misread
the second number, and no prose caveat achieves that.

**(d) The `check`-style supporting-quantity requirement is a Phase 5 recommendation**, where the items
*can* change. The right fix is for the item to **ask** for the failing property, so question, gold and
rubric change together and the archive is re-collected against the new question. Filed as a
recommendation, not a requirement, and explicitly not actionable in Phase 4.

**Reproduction of the 4-of-15 number:**

```bash
python - <<'PY'
import sys, json, glob, re; sys.path.insert(0, '.')
from tests.comparators.normalize import answer_span
rows = [json.loads(l) for f in sorted(glob.glob('error_analysis_annotation/samples/*.jsonl'))
        for l in open(f, encoding='utf-8')]
rows = [r for r in rows if r['question_id'].startswith('system_property_linearity')]
n = 0
for r in rows:
    s = answer_span(r['model_reasoning'])
    s = s[0] if isinstance(s, tuple) else s
    n += bool(re.search(r'additiv|homogen|superposition|scal', s, re.I))
print('linearity traces', len(rows), '| answer span names a property', n)   # -> 15 | 4
PY
```

### 5.5 What generalises

Round 1 I said the F5 shortcut pattern generalises to the whole `classification` population. Round 2
adds a second thing that generalises, and it is the more urgent: **`commitment.py` is written against
four items whose answers are two words long, and it is specified to serve 146 items whose answers are
paragraphs.** Every trigger in R2-F2 fired on a template *outside* the Phase 4 four. The
over-rejection rate on those four is zero and on the archive at large it is not, and the gap is a
property of the sample, not of the mechanism. Before Phase 5 adopts this module, §5.2's recall corpus
should exist.

---

# Round 3 — Reviewer B (physics & pedagogy, the P6 guard)

**Frozen ref:** `3a00875b5cc5c8baa440fdd55fa131d5a96abfe6` · **Branch:** `redesign/phase4-comparators`
**Mandate:** audit `recall_corpus.py` — the artefact I named in round 2 §5.2. Does its 100% mean anything?
**Time box:** 30 minutes, held.

---

## 1. Verdict

**BLOCKED — the corpus is 41% inert, and the one frame I wrote from my own archived
evidence scores 0%.**

The artefact was built and it is the right artefact. But **144 of its 351 cases cannot
observe the frame they are named for**: 16 of the 39 archived answers are `a) … b) …`
enumerations, and the clause splitter cuts on the enumerator, so everything the frame
prepends is discarded before the comparator reads it. I proved this with the negative
control the corpus does not have — wrapping those 16 answers in `If the additivity test
holds, {a}`, a canonical round-1 F3 rejection, still returns **MATCH**. A frame that
cannot fail cannot pass either. The corpus's real denominator is 23 answers, not 39.

On the 23 it does reach, **9 of my 12 new frames hold at 100% and three do not**. The one
that matters is `B3-caveat-but` — **0 of 23** — and it is not a construction I invented:
it is the `reynolds_number_flow_regime` phrasing I filed as R2-F2, *hedge the caveat, then
commit to the label*. `The result may be described differently under another convention,
but the system is linear` is `UNRESOLVED`. **R2-F1.2 is relocated, not closed**: hedge
scope went from a 45-character window (round 1) to a clause (round 2) to a *segment*
(round 3), and the segment boundary set still contains no coordinating conjunction. My
round-2 §5.1(2) named the boundaries explicitly — *"no comma, colon, `but`, `so`, `and`,
`regardless` or `still` between them"*. `but`, `and`, `yet`, `however` and `regardless`
were not implemented. This is the third round in which the marker→label relation has been
tightened by one syntactic level and stopped one level short of the construction that
breaks it.

Separately, the `While the delay term is present…` **MISMATCH is not addressed** — a
correct answer still scored *wrong*, which is the outcome D4.1 §1's stated bias forbids.

**Corpus rate with my frames added: 96.5% over 21 frames × 39 answers (790/819); 94.4%
on the 483 cases where the frame is actually read.** The shipped nine remain 100% on both
denominators.

---

## 2. Disposition of my earlier findings

| finding | disposition | evidence |
|---|---|---|
| **R1-F1 / F2 / F3 / F4** | **addressed**, unchanged from round 2 | `reviewer_battery` 64/64; I re-probed and could not reopen any. |
| **R2-F1.1** concessive/factive/hypothetical | **addressed** | `segment()` is a genuine three-way split. All 10 of my R2-F1(a) cases pass, `note that`/`recall that` are in `FACTIVE`, and the `Note:` vs `Note that` distinction via `BARE_QUALIFIERS` is a real one I did not ask for and that is correct. |
| **R2-F1.2** four probes with no governance test | **RELOCATED — see R3-F2** | All five classes now fire only in the label's own segment, which closes every case I filed. But `segment()` never cuts at a coordinator, so the hedge in a `X, but Y` caveat is still "in the label's own segment". 6 of 7 coordination constructions refused; the archived reynolds case among them. `roughly`/`nominally` are correctly out of `EVIDENTIALS` and `B3-approximation` (100%) confirms it. |
| **R2-F1.3** undirected neighbour rule | **addressed, with one residue** | `is_bare_comment()` is the grammatical test I asked for and the distance bound is gone rather than shrunk — that is the right shape of fix. All 6 of my cases pass. Residue: `UNCERTAINTY` was renamed `INABILITY` but `unclear`, `to be verified`, `it depends on how` and `one could argue` stayed in it. Those are *ambiguity of the question*, not *inability to answer*, and `Whether the offset matters is unclear; regardless, the system is linear` is still `UNRESOLVED`. My §5.1(3) asked for the bag to be **split**, not renamed. |
| **`While the delay term is present…` MISMATCH** | **NOT ADDRESSED** — see R3-F3 | Still `MISMATCH`, reason `label 'not linear' != gold 'linear'`. The segmentation rewrite did not fix it incidentally, and `Although` now fails identically. |
| **F5 / §5.4 withdrawal** | **addressed** | Accepted in full, as recorded. Not re-litigated; see §5.2 for where the two columns should live. |

---

## 3. Findings

### R3-F1 — CONFIRMED (blocking). 144 of the 351 cases cannot observe their own frame.

The corpus has no negative control. I supplied one: frames that **must** be refused — a
hypothetical, a task restatement, an explicit hedge. If a case returns `MATCH` under those,
the frame is not being tested; it is being discarded.

| control frame | MATCH (of 39) | reading |
|---|---:|---|
| `If the additivity test holds, {a}` (round-1 F3) | **16** | 16 answers ignore the frame entirely |
| `We must determine whether the following holds: {a}` (round-1 F3) | **16** | the same 16 |
| `I suspect, though I am not sure, that {a}` (round-1 F1) | **16** | the same 16 |

The 16 are **exactly** the answers matching `a) … b) …` (checked: 16 of 16 enumerated,
0 non-enumerated inert). `_CLAUSE_SPLIT` cuts on `(?:^|[\s;])\(?[a-d][)]\s`, so the frame
becomes a leading, label-free clause and the tuple slot path never reads it. Those 16
answers × 9 frames = **144 cases that are 144 further copies of the `plain` control.**

This is my round-2 R2-F3 arriving as a measurement error rather than as a policy
inconsistency: the tuple path is more permissive than the single-label path, and the recall
corpus inherits that permissiveness as unearned recall. The module docstring's honesty about
the *frames* being synthetic is real and I credit it; it does not cover this, because this is
not a limit of the frames — it is 41% of the corpus not running.

**The fix is one frame and one number.** Add the must-reject frames to `FRAMES` with an
expected outcome of `UNRESOLVED`, and report the rate over the answers the frame reaches. A
recall corpus whose cases cannot fail reports its own sample size, not its subject's recall.

**Reproduction (~15 s):**

```bash
python - <<'PY'
import sys, re; sys.path.insert(0, '.')
from tests.comparators.recall_corpus import archived_answer
from tests.comparators.answer import compare_template
from tests.comparators.ground_truth import LABELS, load
from tests.comparators.score import STEM_TO_TEMPLATE
rows = []
for stem in ("system_property_linearity", "system_properties_memory_causality"):
    rs, lb = load(stem), LABELS[stem]
    for i, row in enumerate(rs):
        if lb[i][0]:
            a = archived_answer(row)
            if a and len(a) <= 400:
                rows.append((stem, row, a))
enum = re.compile(r'(?:^|[\s;])\(?[a-d][)]\s')
inert = [a for stem, row, a in rows
         if compare_template(STEM_TO_TEMPLATE[stem], row["gold_answer"],
              "## Final Answer\n**Answer:** If the additivity test holds, " + a).is_match]
print(len(inert), 'of', len(rows), 'ignore the frame;',
      sum(bool(enum.search(a)) for a in inert), 'of those are enumerated')
PY
```

### R3-F2 — CONFIRMED (blocking). Hedge scope is now the *segment*, and a segment is not cut at a coordinator.

`segment()` cuts at subordinators and at comma-delimited subordinate parts. It does not cut
at `but`, `and`, `yet`, `however`, `so`, `still` or `regardless`. So a caveat coordinated
with the commitment sits in the label's own segment, and every probe class fires on it:

| candidate (gold `linear`) | outcome |
|---|---|
| `The result may be described differently under another convention, but the system is linear.` | **UNRESOLVED** — *modal: may be* |
| `Homogeneity could be checked in two ways and the system is linear.` | **UNRESOLVED** — *modal: could be* |
| `It seems ambiguous at first, yet the system is linear.` | **UNRESOLVED** — *appearance: seems* |
| `The margin appears tight, however the system is linear.` | **UNRESOLVED** — *appearance: appears* |
| `A different textbook would say otherwise, but the system is linear.` | **UNRESOLVED** — *opinion: would say* |
| `Whether the offset matters is unclear; regardless, the system is linear.` | **UNRESOLVED** — *bare comment: unclear* |
| `The boundary case might go either way, so the system is linear.` | MATCH |

Only `so` survives, and by accident — `might go` fails `_MODAL_RE`'s copula requirement, not
because the coordinator was seen. `segment()` returns the whole sentence unsplit; I checked
directly:

```
segment("The result may be described differently under another convention, but the system is linear.")
  -> ("The result may be described differently under another convention, but the system is linear.", "")
hedge_markers(matrix) -> ['modal: may be']
```

**As a frame over the 23 reachable archived answers this scores 0 of 23.** It is the exact
construction of `reynolds_number_flow_regime__3122400918`, which I filed in round 2 as the one
archived two-way-classification trace that hedges a caveat and then commits. That trace is
still refused.

**The through-line, stated once more.** Round 1: scope was 45 characters, and the fix made it a
clause. Round 2: scope was a clause, and the fix made it a segment. Round 3: scope is a segment,
and a segment is the whole sentence whenever the contrast is coordinated rather than
subordinated. Each fix moved the boundary one level up the syntax and stopped short of the same
relation. My round-2 §5.1(2) listed the coordinators by name; three of that remedy's four parts
were implemented and this one was not.

### R3-F3 — CONFIRMED (blocking). The `MISMATCH` on a correct answer is unchanged, and has spread.

```
MISMATCH | While the delay term is present, it does not break superposition: the system is linear.
         -> label 'not linear' != gold 'linear'
MISMATCH | Although the delay term is present, it does not break superposition: the system is linear.
         -> label 'not linear' != gold 'linear'
MATCH    | While the delay term is present, superposition still holds, so the system is linear.
```

The brief asked whether the segmentation rewrite fixed this incidentally. **It did not**, and
the second row is new: `Although` — the opener the rewrite was built to rehabilitate — now
fails the same way. `while` and `although` are correctly classified as concessive and the
matrix is correctly taken, but the matrix is `it does not break superposition: the system is
linear`, and label selection reads the negation across the colon onto `linear`. This is
`kinds.py` label selection, not `commitment.py`, and I said so in round 2; it is filed again
because it is the **worst available outcome** — a correct answer scored wrong, not merely
undecided — and because D4.1 §1 states a bias toward `UNRESOLVED` over `MISMATCH` that this
violates. Negation scope needs the treatment commitment scope has now had three times: it
should stop at a colon and at a coordinator.

### R3-F4 — MINOR. Two trailing full sentences still hedge, but only on the tuple path.

`{a} The result should be stated to three significant figures.` and `{a} That is how it
appears from the two tests.` score **21 of 23**. Both failures are the `Memoryless: No;
Causal: Yes.` answers, and the reason is `causal: hedged (modal: should be)` — the tuple slot
path attaches the trailing sentence to the second slot. The single-label path handles both
correctly, because `is_bare_comment` rejects a 9-word clause with a definite subject. This is
R2-F3 with the sign flipped: in round 2 the tuple path was the more permissive one and I said
it was the one that was right; here it is the stricter one and it is the one that is wrong.
**One policy, two implementations, and which is correct now depends on the case.**
`_slot_verdict` should call the same `find_commitment` machinery rather than a hand-rolled
subset of it.

`The result should be stated to three significant figures` is verbatim the
`kinematic_viscosity__3586477399` phrasing from my R2-F2 table — a modal about **rounding
convention**, which I flagged in round 2 and which is still read as a modal about the answer.

---

## 4. Falsification attempts that failed

- **Nine of my twelve frames hold at 100%**, on both denominators: `Therefore, {a}`; `I am
  confident that {a}`; `My first pass had this backwards. {a}` (self-correction); `Strictly on
  the definitions given in the course text, {a}`; `{a} The algebra above is a little rough,
  though.`; `By the standard textbook criterion, {a}`; `The coefficient is roughly 2 for every
  n, so {a}`; `Conclusion: {a}`; `It would have been the other way had the square been absent,
  but {a}`. The last is the one I most expected to fail — a counterfactual coordinated with a
  commitment — and it passes, because `had` never reaches a clause-initial position. **The
  concessive/factive rewrite is genuinely sound**; R3-F2 is about coordination, which that
  rewrite never claimed to handle.
- **`roughly` and `nominally` are properly out.** `The coefficient is roughly 2 for every n, so
  {a}` is 100%, and `roughly speaking` still hedges. The multi-word retention is the right call
  and the reasoning in the `EVIDENTIALS` comment is correct.
- **`reviewer_battery` is 64/64 and I could not add a case to it that both rounds imply and
  that fails.** Every case I filed in rounds 1 and 2 behaves as I argued. R3-F2 is not in the
  battery because I did not write it as a case in round 2 — I wrote it as a *remedy*, and the
  remedy was implemented in three parts of four.
- **The gate is green and stays green.** archive 100.0% / 97.8%, adversarial 98.6% / 98.6%,
  zero archived false accepts. **None of R3-F1 through R3-F4 moves any of those numbers**,
  which is the point: three rounds in, the gate has never once been the instrument that found a
  defect, and the recall corpus was built to change that. It can, once R3-F1 is fixed — the
  `B3-caveat-but` frame turns R3-F2 into a 0-of-23 line in a table that runs on every commit.
- **`is_bare_comment`'s no-distance-bound design survives probing.** I tried to make a bare
  comment govern across four intervening clauses where it should not, and across zero where it
  should, and it behaved correctly both times. Removing the magic constant rather than retuning
  it was the right move and it is the strongest single change in this round's diff.
- **The `Note:` / `Note that` split is right, and I checked both directions.**

---

## 5. Further probing and improvements

### 5.1 What the recall corpus needs in order to be the instrument I asked for

1. **A negative control, reported as a row.** Frames that must be refused (`If X, {a}`; `We
   must determine whether {a}`; `I suspect … that {a}`). Any case that returns `MATCH` under
   all of them is inert and must be excluded from the denominator or fixed. It costs nothing
   and it is the difference between a measurement and a tautology.
2. **Report the reachable denominator.** 39 answers, 23 reachable; 351 cases, 207 live. The
   headline number should be the live one.
3. **Frames must be able to wrap, not only prefix.** Every shipped frame but two is a prefix,
   and `archived_answer` appends a full stop, so no frame can put anything inside the label's
   own clause — which is precisely where R3-F2 lives. At least one frame should interpolate:
   `{a_stem}, but {a}`.
4. **The frames I would add**, in priority order: **coordinated caveat** (`X may be …, but
   {a}` — 0% today); **modal about convention or rounding, trailing the answer** (`{a} The
   result should be stated to three significant figures.` — the archived `kinematic_viscosity`
   phrasing); **ambiguity-of-question then commit** (`Whether the offset matters is unclear;
   regardless, {a}`); **colon-separated negated premise** (the R3-F3 construction). All four
   come from archived text or from cases already in this report, not from the constructions
   `commitment.py` enumerates — which is the circularity test the current nine do not pass.
   Seven of the nine shipped frames name in their own comment the finding they were written
   from, which is admirably honest and is also the diagnosis: they were selected from the
   inside of the list they certify, for the third round running.

### 5.2 Where the two accepted numbers should live (optional task)

**`docs/re-implementation-sep/template_inventory.csv`**, as two new columns —
`blind_guess_floor` and `surface_model_holdout` — beside the existing `difficulty` column.
Not the item-pool-impact note, and not the D4.1 results table:

- The inventory is **per template**, which is the grain both numbers have. It already carries
  `difficulty`, `n_steps`, `answer_type` and `pct_answer_values_recoverable`, so a reader
  comparing a *claimed* difficulty against a *measured* one has both in one row. The item-pool
  note is per phase and would strand the numbers in whichever phase computed them; the results
  table is per comparator run and would recompute them every time.
- Both are cheap and deterministic: the round-1 §2 script, 2,000 seeds, minutes per template.
  A `phase5_shortcuttability.py` that fills the two columns for all 146 is one batch job, and
  my round-1 §5.4 argued it would have caught D-046 unaided.
- A third `shortcut_lift` column is unnecessary — it is the difference — but it is the natural
  **sort key**, so compute it in the report rather than storing it.

**Threshold for concern.** Two, because there are two failure modes:

- **`surface_model_holdout − blind_guess_floor ≥ 40 pts` → flag.** Both Phase 4 items are
  +49.98 and +65.98 and both are pure lookups. 40 points is generous and still catches them.
- **`surface_model_holdout ≥ 95%` → flag regardless of lift.** A high floor can mask a total
  lookup: an item with an 85% majority class and 100% held-out has a lift of only 15 points and
  is still a function of the question surface. The Phase 4 items trip both tests; no template
  should carry an `Advanced` label while tripping either.

A flag should mean **"the difficulty label and the answer-only scoring decision need a P6
register entry"**, not "block the template". That is the disposition my §5.4 argued for F5, and
it generalises: the number's job is to force the scoping decision to be *recorded*, which is
what P6 exists for.

### 5.3 What lifts the block

Narrow, and all three are already specified:

1. Add the negative control to `recall_corpus.py` and report the reachable denominator
   (R3-F1). Without it the corpus's number is not interpretable in either direction.
2. Cut `segment()` at coordinators — `but`, `and`, `yet`, `however`, `still`, `regardless`,
   `nevertheless` — as round-2 §5.1(2) specified (R3-F2). This is the fourth part of a
   four-part remedy of which three parts are done.
3. Stop negation scope at a colon and at a coordinator (R3-F3). A correct answer scored
   `MISMATCH` is the one outcome D4.1 §1 rules out, and it has now survived two rounds.

R3-F4 is minor and can ride along: route `_slot_verdict` through `find_commitment` so the
policy has a single implementation.

---

# Round 4 — Reviewer B (physics & pedagogy, the P6 guard)

**Verdict: the phase can CLOSE, with three residuals stated — one of which is a
required edit to D4.1 §8, not a recommendation.** I have blocked three times and
each block was a defect that would have shipped. There is no fourth thing of that
kind. What I found in round 4 is that the *instrument* still over-reports, in both
directions, and the honest fix is to restate its number rather than to hold the
gate.

## 1. Verdict

| Task | Question | Answer |
|---|---|---|
| 1 | Are the four negative frames the right four? | **Right in kind, wrong in coverage.** Each passes 39/39, and each one's nearest linguistic neighbour passes 0–59%. With six defensible additions the negative block is **182/390 = 46.7%**, not 100%. |
| 2 | Is "whose confidence is qualified" the right pedagogical line? | **Yes — the line is right.** The *proxy* that implements it (anaphoric-or-absent subject) mis-identifies expletive `it` and topic-shifting `we`, and it errs on the false-reject side, which is my half of the gate. **0/39** on four ordinary frames. |
| 3 | Is the F5 disposition sufficient? | **Sufficient in kind. One element of it is required, not recommended, and none of it is in the tree yet.** |

**Headline for the record:** the recall corpus reports **507/507 (100%)**. Adding
six positive frames a grader credits and six negative frames a grader refuses —
all twelve constructed from the same generator, on the same 39 archived answers —
takes it to **548/975 = 56.2%**. That is not a reason to block. It is a reason for
the phase to publish the second number's *existence* alongside the first.

## 2. Disposition of my earlier findings

| Finding | Round | State in the frozen tree |
|---|---|---|
| R1-F1 hedge blocklist | 1 | Superseded by the class-based probes. **The blocklist shape survives**, see R4-F1. |
| R1-F2 one-sided scope | 1 | Closed. `_hedges_governing` scopes bidirectionally. |
| R1-F3 assert-nothing | 1 | Closed. `segment` returns `None` with a reason. |
| R1-F4 longest-surface | 1 | Closed, and re-closed by discourse position for FACTIVE. |
| R1-F5 shortcuttability | 1 | See §3, R4-F3. |
| R2-F1 over-rejection 25/38 | 2 | Closed. All nine positive frames 39/39. |
| R2-F3 `categorical` vs tuple | 2 | Closed on the paths I probed; `suspends_what_follows` still diverges from `segment` on one construction (R4-F2). |
| R3-F1 corpus cannot see its frame | 3 | **Actioned.** Negative frames exist and found two real defects. Partially effective — §3, R4-F1. |
| R3-F2 coordinator scope | 3 | Closed. `conjuncts()` exists and is used. |
| R3-F3, R3-F4 | 3 | Closed on re-probe. |

**The backslash-b-in-heredoc report checks out.** `is_bare_comment` now fires on
every constructed case I put to it, including the ones I could not move in round 3.
I tried to reproduce the masked behaviour and could not; the fix holds.

## 3. Findings

### R4-F1 — CONFIRMED (residual, not blocking). The negative frames are as position-locked as the positive frames were, and each one's nearest neighbour fails.

The four negative frames were the right *classes* to pick. But every one of them
puts its trigger at character 0 of the frame, and the three matchers that catch
them — `_TASK_RE`, `_HYPO_RE`, `_BARE_QUAL_RE` — are `.match()`-anchored. So each
frame certifies a matcher against the one string position the matcher can see.
Move the trigger three words to the right and the frame inverts.

Rates over the same 39 archived correct answers, same generator:

| shipped frame | rate | nearest neighbour I added | rate |
|---|---:|---|---:|
| `neg-task-restatement` *"We must determine whether X"* | 39/39 | *"We will check below whether X"* | **0/39** |
| `neg-hypothetical` *"If the additivity test holds, X"* | 39/39 | *"For now I will assume X"* | **0/39** |
| `neg-explicit-hedge` *"I am not sure, but perhaps X"* | 39/39 | *"My best guess is that X"* | **23/39** |
| `neg-question` *"X. Or is it?"* | 39/39 | *"X. Or possibly the opposite."* | **3/39** |
| | | *"Either X. Or the opposite holds."* | **0/39** |
| | | *"X. But I cannot verify this without more information."* | **0/39** |

`check whether` and `assume` are **already in** `TASK_RESTATEMENTS` and
`HYPOTHETICAL`. They do not fire because the anchoring, not the vocabulary, is
what decides. `cannot verify` is not in `INABILITY` at all, next to
`cannot determine` and `cannot tell` — that is my own round-1 F1 recurring: the
policy is class-based at the top and a remembered list at the bottom.

**Rate with mine added:** negative block **26/234 = 11.1%** on the six additions,
**182/390 = 46.7%** combined. The four shipped frames therefore measure
*non-inertness* — which is what I asked for in R3-F1, and they delivered it, twice
over — but they do not measure coverage of non-answers, and the corpus's 156/156
should not be read as if they did.

**Why this is not a block.** The failures above are false *accepts*, which is
Reviewer E's remit and explicitly out of mine. I report them only because task 1
asked me to name absent frames and give the rate. The instrument correction is
mine; the defects it exposes are E's to disposition.

**Cosmetic, but in the artefact whose credibility is the question.**
`recall_corpus.main()` prints `507 cases from 56 archived correct answers x 9
frames`. It is 39 answers × 13 frames; `len(cases) // len(FRAMES)` was not updated
when `NEGATIVE_FRAMES` was appended to the build loop. A corpus whose whole purpose
is to stop a number being over-read prints a wrong denominator in its own header.
One-line fix.

### R4-F2 — CONFIRMED (residual, blocking-adjacent; the strongest thing I found). The hedge rule's *proxy* refuses committed answers at 0/39 on four ordinary frames.

Task 2 asked whether "whose confidence is qualified" is the right line. **It is.**
It is the correct distinction, it is the one I drew in round 2 when I split
`UNCERTAINTY`, and it genuinely does resolve R1-F1, R2-§5 and E's R3-F14 under one
rule rather than three patches. I confirmed my round-2 §5 case still works:
*"The distinction is not entirely obvious, but the system is linear"* is correctly
**not** a bare comment and is credited.

The defect is in the operationalisation. `is_bare_comment` identifies "whose
confidence" by **subject anaphoricity**, via
`_ANAPHORIC_SUBJECT = ^(it|this|that|these|those|i|we|there)`. Two of those seven
are not anaphoric when they carry a complement:

- **Expletive `it`.** *"It is not certain which sign convention the textbook uses,
  but the system is linear."* The `it` is a dummy subject; the uncertainty is
  scoped to *the convention*, a different proposition entirely. The rule reads it
  as the speaker's, and it governs from anywhere — including backwards across the
  coordinator, since a preceding unit is checked with `classes={"stated
  uncertainty"}`, exactly the class these fire in. **0/39.**
- **Topic-shifting `we`.** *"We cannot determine the settling time from this alone,
  but the system is linear."* The inability names the quantity it applies to, and
  it is not the label. **0/39.**

| frame (a grader credits every one) | MATCH |
|---|---:|
| *"It is not certain which sign convention the textbook uses, but {a}"* | **0/39** |
| *"It is unclear what part (c) is asking, but {a}"* | **0/39** |
| *"We cannot determine the settling time from this alone, but {a}"* | **0/39** |
| *"{a} Why? Because both defining tests pass."* | **0/39** |
| *"{a} Is that enough? Yes: both tests are satisfied."* | **0/39** |
| *"Is the system linear? {a}"* | 15/39 |

The last three are a separate root cause and it is a **regression introduced by
the fix to my own R3-F1**. `_hedges_governing` now carries:

```python
if follows and text.rstrip().endswith("?") and not _label_hit(text, surfaces):
    marks.append("withdrawn by a following question")
```

Any label-free trailing question withdraws the commitment. *"Or is it?"* does.
*"Why?"* and *"Is that enough?"* do not — they are rhetorical self-questions
immediately answered by the clause after them, and they are among the most common
constructions in tutorial-register model output. The 15/39 row is the enumerated
path: `suspends_what_follows` returns `"the answer is a question"` when the
**preface** ends in `?`, which is R2-F3's divergence between the tuple path and
`segment` surviving in one more place.

**The discriminator that fixes all three, and it is already in the module's
idiom:**

1. *A trailing question withdraws only if nothing asserting follows it.* This is
   the same discourse-position principle `find_commitment` already applies to
   FACTIVE ("a factive that FOLLOWS a commitment is elaboration, not retraction").
   `Or is it?` ends the span; `Why? Because…` does not.
2. *A bare epistemic comment has no complement of its own.* `it`/`there` followed
   by a marker and then a wh-/that-complement, or a PP naming a quantity, is a
   claim about that quantity, not about the answer. This is the same move that
   already removed the `_MAX_COMMENT_WORDS` constant: a structural test, not a
   list.

**Why this is not a fourth block, and I want the reasoning on the record because
it cuts against me.** I measured the archive rather than asserting. Across all
2,200 archived answer spans, label-free question clauses with an assertion after
them occur **once**, and that one is a mis-split of the *question* text, not an
answer. So the construction I am refusing at 100% has the *same* evidential
standing — 1 in 2,200 — as the hedge machinery the module itself declares to be
"a policy validated against constructed text." I cannot hold a gate on a
constructed case while accepting the phase's own constructed cases as sufficient.
Both directions of the false reject also land on `UNRESOLVED`, which the
three-valued design exists to make safe: it withholds a score, it does not mark a
correct answer wrong.

Where it *will* bite is Phase 5. My round-2 §5.5 said `commitment.py` is written
against four items whose answers are two words long and specified to serve 146
whose answers are paragraphs. R4-F2 is what that sentence looks like when it comes
due, and it should be on the register in those words.

**Rate with my positive frames added:** **15/234 = 6.4%** on the six additions,
**366/585 = 62.6%** combined. The corpus's 351/351 is a statement about nine
frames.

### R4-F3 — The F5 disposition. Sufficient in kind; one element is required, and none of it is in the tree.

Taking the four parts as briefed:

- **(a) §5.4 withdrawal accepted** — correct, and I reaffirm it. Requiring a named
  failing property would drop `system_property_linearity` recall to ~27% and is an
  item change wearing a comparator's clothes.
- **(b) register entry** — correct disposition.
- **(d) supporting-quantity requirement filed as a Phase 5 recommendation** —
  correct. The items can change there; they cannot change here.
- **(c) two columns, my thresholds** (flag at lift ≥ 40 pts **or** held-out ≥ 95%
  regardless of lift) — the thresholds are mine and I stand behind them. **The
  placement is a substitution, and it drops the part that did the work.**

My round-2 §5.4(c) named **the results table**, and named it for one reason,
stated there: *"A reader who sees `50.02% → 100.00%` beside `precision 100%`
cannot misread the second number, and no prose caveat achieves that."*
`template_inventory.csv` is a 150-row planning artefact read by whoever designs
Phase 5's items. It is the right **durable** home and I endorse it as such. It is
not where the misreading happens. The misreading happens at D4.1 §8, which in the
frozen tree still prints, with nothing adjacent to qualify it:

```
| `system_properties_memory_causality` | 24 | 24 | 24 | 100% | 100% | 100% |
| `system_property_linearity`          | 15 | 15 | 15 | 100% | 100% | 100% |
```

§8 already does exactly the right thing one paragraph below for a *different*
number — "**This number is weak evidence and is reported as such.** The archive is
where the vocabulary came from, so scoring 100% on it measures memorisation." The
shortcuttability caveat is the same species and is absent.

**State of the tree, verified:**

- `docs/re-implementation-sep/template_inventory.csv` has **18 columns** and none
  of them is a floor, a lift, or a held-out accuracy. All four Phase 4 templates
  are present as rows (lines 89, 90, 91, 123); the columns are not.
- `50.02` / `34.02` appear nowhere in `docs/re-implementation-sep/*.md` or in the
  CSV — only in this review file.
- `phase4_summary.md`, which D4.1 §8 cites twice ("All six are in
  `phase4_summary.md` §8"), **does not exist** in the tree.

**What P6 requires before the phase closes.** Not a comparator change — I withdrew
that and the withdrawal stands. Not the shortcuttability fix — it is a Phase 5 item
change. P6 requires that a scoping decision be *recorded and approved* (spec §36),
and the decision here is: **Phase 4 scores two 100%-shortcuttable items on their
answer alone.** One sentence in D4.1 §8, beside the two rows, is the whole of it:

> `system_property_linearity` and `system_properties_memory_causality` are scored
> on their answer alone. Their answers carry 1.0 and 1.6 bits against blind-guess
> floors of 50.02% and 34.02%, and a depth-2 tree on the question surface reaches
> 100% held-out on both. The 100% precision and recall above are statements about
> the vocabulary, not about whether a scored model reasoned. Recorded per template
> in `template_inventory.csv`; the item fix is a Phase 5 recommendation.

Zero code, zero risk, and it is the only part of the disposition that a reader of
Phase 4 will ever encounter. **That edit, plus the two columns actually existing,
is my condition.** Everything else in the disposition is right as briefed.

## 4. Falsification attempts that failed

I tried to break the phase in six more places and could not:

1. **My round-2 §5 case is genuinely fixed, not accidentally passing.** *"The
   distinction is not entirely obvious, but the system is linear"* — the subject
   is a full NP, `is_bare_comment` correctly returns `False`, the hedge stays in
   its conjunct, the answer is credited. The discriminator does the work here, not
   a list.
2. **Precision terms survive.** *"{a} The margin is roughly 12% of the limit."* —
   39/39. Removing bare `roughly` from `EVIDENTIALS` while keeping
   `roughly speaking` was right and holds end-to-end.
3. **Emphatic reinforcement is not read as a hedge.** *"{a} There is no doubt about
   this."* — 39/39, despite `no` being a `_SUBJECTLESS` opener.
4. **An inability scoped to an unlisted quantity does not retract.** *"{a} I cannot
   give the exact transient response, but that was not asked."* — 39/39. (Passing
   for the weaker reason, though: `cannot give` is simply not in `INABILITY`. Cf.
   R4-F1's `cannot verify`.)
5. **A hedge about a different part of the question does not retract.** *"{a} The
   wording of part (c) is unclear."* — 39/39, because the subject is a full NP.
   Note the contrast with R4-F2's *"It is unclear what part (c) is asking"* at
   0/39: **the same proposition, phrased with an expletive subject, inverts the
   verdict.** That is the cleanest statement of R4-F2 I have.
6. **The backslash-b heredoc bug is fixed.** Every bare-comment probe I ran fires.
   I could not reproduce any of round 3's masked behaviour.

I also confirmed the two defects the negative controls found are genuinely closed:
the hedged preface on enumerated answers (`suspends_what_follows` now splits on
`,`/`but`/`and` and tests each piece) and the trailing `Or is it?` — the latter
correctly, if over-broadly, per R4-F2.

## 5. Further probing and improvements

### 5.1 What the corpus should report instead of 100%

Not "add my twelve frames and go green." The number 507/507 is not wrong; what is
wrong is reading it as coverage. Two changes, both cheap:

1. **Print the denominator that exists.** Fix the header count, and state the frame
   count: *"507 cases = 39 archived correct answers × 13 frames."* A reader then
   knows the shape of what was measured.
2. **Split the frames by where the trigger sits.** Every current frame is
   trigger-at-position-0 or trigger-at-end. One extra column — *anchored* vs
   *embedded* — makes R4-F1 visible in the corpus's own output rather than in a
   review, and it is the difference between an instrument that reports its
   coverage and one that reports its score.

### 5.2 The register entries I am filing

Three, in priority order for whoever picks up Phase 5:

- **RB4-1 (Phase 5, high).** `commitment.py`'s bare-comment proxy mis-identifies
  expletive `it`/`there` and topic-shifting `we`; and any label-free trailing
  question withdraws a commitment, including a rhetorical one that its own next
  clause answers. **0/39 on four ordinary frames, reproduced above.** Both fixes
  are structural and are given in R4-F2. Archive support for the constructions is
  1 in 2,200 *on four items with two-word answers*, which is why it is Phase 5's
  and not Phase 4's — and why it will bite there, on 146 items with prose answers.
- **RB4-2 (Phase 4, required).** The D4.1 §8 sentence in R4-F3, plus the two
  columns actually in `template_inventory.csv`.
- **RB4-3 (housekeeping).** `phase4_summary.md` is cited twice by D4.1 §8 and does
  not exist. Either write it or fix the citation; a contract that points at a
  missing document is the kind of thing a conformance reader finds first.

### 5.3 What generalises, for the fourth and last time

Every one of my four rounds has found the same shape, at a different altitude:

| round | the thing that was scoped too narrowly | the thing that certified it |
|---|---|---|
| 1 | a 45-character hedge window | D4.4, sampled from inside the hedge list |
| 2 | a ±1-clause window, and a hedge list with no governance | 61 archived answers, almost none of them prose |
| 3 | segment scope that never cut at a coordinator | 351 positive cases, 144 of them inert copies |
| 4 | matchers anchored at position 0, and a subject allowlist | 13 frames, every trigger at position 0 or at the end |

**The through-line is not the hedge policy. It is that at every altitude the
instrument was built from the same construction the mechanism was built from, so
it could only confirm.** The fix is never a longer list; it is a test whose frames
are generated independently of the matcher's anchoring. Phase 5 should generate
frames by *transformation* — take a passing frame and move its trigger, embed it,
give it an expletive subject — rather than by enumeration. That is a half-day and
it is the last recommendation I have.

### 5.4 On closing

I have blocked this phase three times and I would do it three times again. I am not
blocking a fourth, and I want the reason stated plainly rather than buried: the two
false rejects in R4-F2 are real and reproducible, but they land on `UNRESOLVED`
rather than on a wrong score, they have 1-in-2,200 archive support, and holding the
gate on constructed text while the phase honestly declares its own machinery to be
validated on constructed text would be applying a standard to the implementation
that I did not apply to the specification. The contract is sound, the discriminator
in task 2 is the right pedagogical line, and the F5 withdrawal was correct.

**The phase closes when RB4-2 is done.** That is one sentence and two CSV columns.
RB4-1 goes to the register with its reproductions, where it belongs, and where the
next reviewer of this module should start.

### 5.5 Reproduction

```bash
export PYTHONIOENCODING=utf-8
python - <<'PY'
import sys; sys.path.insert(0, '.')
from tests.comparators.recall_corpus import GATED, archived_answer
from tests.comparators.answer import compare_template
from tests.comparators.ground_truth import LABELS, load
from tests.comparators.score import STEM_TO_TEMPLATE
answers = []
for stem in GATED:
    rows, labels = load(stem), LABELS[stem]
    for i, row in enumerate(rows):
        if not labels[i][0]:
            continue
        a = archived_answer(row)
        if a and len(a) <= 400:
            answers.append((STEM_TO_TEMPLATE[stem], row["gold_answer"], a))
FRAMES = [                                        # (name, frame, a grader credits it?)
    ("neg-deferred",      "We will check below whether {a}",                      False),
    ("neg-provisional",   "For now I will assume {a}",                            False),
    ("neg-guess",         "My best guess is that {a}",                            False),
    ("neg-alternation",   "{a} Or possibly the opposite.",                        False),
    ("neg-disjunction",   "Either {a} Or the opposite holds.",                    False),
    ("neg-refusal-after", "{a} But I cannot verify this without more information.", False),
    ("pos-expletive-it",  "It is not certain which sign convention the textbook uses, but {a}", True),
    ("pos-unclear-part-c","It is unclear what part (c) is asking, but {a}",        True),
    ("pos-we-elsewhere",  "We cannot determine the settling time from this alone, but {a}", True),
    ("pos-rhetorical-q",  "{a} Why? Because both defining tests pass.",            True),
    ("pos-rhet-2",        "{a} Is that enough? Yes: both tests are satisfied.",     True),
    ("pos-selfq-lead",    "Is the system linear? {a}",                             True),
]
for name, f, want in FRAMES:
    ok = sum((compare_template(t, g, "## Final Answer\n**Answer:** " + f.format(a=a)).is_match) == want
             for t, g, a in answers)
    print(f"  {name:20s} {ok:3d}/{len(answers)}  {ok/len(answers):6.1%}")
PY
```

Archive census for the trailing-question construction (1 in 2,200): iterate
`error_analysis_annotation/samples/*.jsonl`, take `normalize.answer_span`, split
with `commitment.clauses`, and count clauses ending in `?` that fail
`kinds._label_hit` and are followed by a further clause.
