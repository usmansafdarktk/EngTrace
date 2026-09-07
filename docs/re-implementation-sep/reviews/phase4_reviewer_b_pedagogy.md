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
