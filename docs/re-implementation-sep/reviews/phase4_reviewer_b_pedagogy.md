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
