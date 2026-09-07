# D4.2 — The normalisation vocabulary, and what it rests on

**Version:** 1.0 · **Phase:** 4 · **Date:** 2026-09-07
**Companion:** [`phase4_comparators.md`](phase4_comparators.md) (D4.1)
**Source of truth:** `tests/comparators/normalize.py` — the rules and the
`SUPPORT` table live in the module, not here.
**Regenerate every number in this document:**
`python -m tests.comparators.derive_vocabulary`

---

## 1. The rule this document exists to obey, and where it bites

Spec §4.2 requires the vocabulary to be *"specified from real model outputs, not
invented"*. The brief adds the constraint that makes that hard:

> **only 61 of those 2,200 are for your four templates**: `memory_causality` 24,
> `signal_operations` 16, `linearity` 15, **`incompressible_continuity` 6.**
> That last number is the single most important fact in this brief.

It reproduces exactly. And by the standing rule (D-024, D-026) — *N samples
cannot resolve a variant rarer than about 3/N* — the four templates have four
different limits of detection, which is why **every coverage number in this
document is stated per template and none is averaged**:

| template | traces | resolves down to |
|---|---:|---:|
| `system_properties_memory_causality` | 24 | ~12% |
| `signal_operations` | 16 | ~19% |
| `system_property_linearity` | 15 | ~20% |
| **`incompressible_continuity`** | **6** | **~50%** |
| whole archive | 2,200 | ~0.14% |

**Six traces are blind to anything occurring in under half of model outputs.**

---

## 2. Decision 2 — where the vocabulary comes from

The brief offers three routes and requires one to be picked explicitly.

**Route 1 is taken** — draw from the whole 2,200 — **and its argument was tested
rather than asserted, and the test does not support the argument as stated.**

### 2.1 The test

If these surface forms are *model habits*, a rule's rate should be driven far
more by **which model** wrote a trace than by **which template** it answers. So:
for each of 27 rules, compare the spread of its per-model rates with the spread
of its per-template rates across the whole archive.

### 2.2 The result, reported as measured

| verdict | rules |
|---|---:|
| model habit (model spread ≥ 2× template spread) | **2** |
| item-driven (template spread ≥ 2× model spread) | **15** |
| both, or rare everywhere | 10 |

Only `latex-boxed` and `The answer is` are cleanly model habits. **Every Unicode
rule and most LaTeX rules are item-driven.** A model does not write a superscript
because it is that model; it writes one because the item has an exponent in it.

### 2.3 What that invalidates, and what it does not

**It invalidates any coverage claim borrowed across templates.** "The vocabulary
handles 95% of observed forms" cannot be argued from the 2,200 and applied to the
6. Hence §1's per-template table, and hence no averaged number anywhere here.

**It does not invalidate the rules.** A normalisation rule is a claim about
*meaning*, not frequency: `\frac{3x^2}{2}` denotes 1.5x² whoever writes it and
however often. The archive's job for a rule is to show the form exists and to fix
what it means, and one occurrence anywhere in 2,200 traces does both. **That is
the part of route 1 that survives, and it is the part the rules actually use.**

**A confound, named because it weakens the test itself.** The statistic cannot
separate *"the model does not use this form"* from *"the item gave it no
opportunity"*. Superscripts look item-driven partly because
`incompressible_continuity` is the only Phase 4 template with an exponent at all.
The better statistic is the rate **given the opportunity**, which needs an
opportunity detector per rule and is not built. So **15 is an upper bound** on
how item-driven the vocabulary really is.

### 2.4 Residual risk, carried not solved

**The archive cannot tell me which forms are *missing* for
`incompressible_continuity`.** Six traces, item-driven forms so the other 2,194
do not fill the gap, and route 2 (generating more) needs D-003 resolved — the raw
`inference_results/` generations are gitignored and absent from this working
copy.

For the **symbolic comparator specifically** the honest position is close to the
brief's **route 3**: the rules are stated, they are exercised against 22
hand-written adversarial cases, and they are **not validated against a
representative sample of real output**. That is on the residual-risk register, not
a solved problem. Recorded as **D-048**.

---

## 3. The table

Machine-checked. `tests/comparators/normalize.py` carries `SUPPORT` as
`rule -> (count in the 61, count in all 2,200)`, and
`derive_vocabulary.py` recomputes every pair from the archive and **fails on
disagreement**. 27/27 agree today.

Read the two columns together and the second is not reassurance: because rates
are item-driven (§2.2), a large whole-archive count is evidence that a form
*exists*, not evidence about how often it occurs in these four templates.

| rule | sig_ops (16) | mem_caus (24) | linearity (15) | continuity (6) | all 2,200 |
|---|---:|---:|---:|---:|---:|
| `marker-final-answer-heading` | 14 | 24 | 15 | 4 | 1824 |
| `marker-bold-answer` | 14 | 22 | 15 | 4 | 1736 |
| `marker-the-answer-is` | 2 | 0 | 0 | 2 | 202 |
| `negator-not` | 0 | 16 | 9 | 0 | 34 |
| `negator-fused-non` | 0 | 0 | 3 | 0 | 3 |
| `negator-neither-nor` | 0 | **2** | 0 | 0 | 2 |
| `negator-fails` | 0 | 0 | 3 | 0 | 3 |
| `trailing-annotation` | 2 | 10 | 4 | 0 | 86 |
| `brace-sequence` | 16 | 0 | 0 | 1 | 447 |
| `origin-asterisk` | **1** | 0 | 0 | 0 | 19 |
| `explicit-index-list` | 2 | 0 | 0 | 0 | 6 |
| `per-element-assignment` | 2 | 0 | 0 | 0 | 2 |
| `latex-escaped-brace` | 5 | 0 | 0 | 0 | 6 |
| `latex-inline-math` | 2 | 1 | 2 | 0 | 332 |
| `latex-display-math` | 0 | 0 | 0 | 1 | 66 |
| `latex-frac` | 0 | 0 | 0 | **1** | 53 |
| `latex-boxed` | 0 | 0 | 0 | **1** | 104 |
| `arbitrary-function` | 0 | 0 | 0 | **2** | 30 |
| `unicode-minus` | 0 | 2 | 0 | 2 | 71 |
| `unicode-superscript` | 0 | 0 | 0 | **3** | 236 |
| `markdown-bold` | 0 | 2 | 1 | 0 | 239 |
| `bracket-sequence` | 0 | 0 | 0 | 1 | 64 |
| `unicode-operator` | **0** | **0** | **0** | **0** | 330 |
| `latex-text` | **0** | **0** | **0** | **0** | 274 |
| `latex-exponent-brace` | **0** | **0** | **0** | **0** | 165 |
| `hedge` | **0** | **0** | **0** | **0** | 1 |
| `thousands-separator` | **0** | **0** | **0** | **0** | 32 |

### 3.1 The five rules with zero Phase 4 support

`unicode-operator` · `latex-text` · `latex-exponent-brace` · `hedge` ·
`thousands-separator`

Each is carried on the wider archive alone. Each is a place the comparator is
specified against forms these four templates **have never been observed to
produce**. This is the concrete content of §2.4's residual risk, and it is listed
by name rather than left as a caveat.

Two are worth singling out:

- **`hedge` occurs once in 2,200 traces.** The hedge rule — the specific thing
  spec §4.4 asks Reviewer B to guard — is supported by a single observation. It
  is a *policy* rather than an observation-derived rule, and it is fair to
  attack on exactly that ground.
- **`latex-exponent-brace` occurs 165 times in the archive and zero times in
  Phase 4.** The rule that handles it was added because adversarial case
  `sym-04` broke without it — evidence from D4.4, not from the archive, which is
  the intended division of labour between the two corpora.

---

## 4. The rules the archive genuinely establishes

Six rules exist because the archive would be misread without them. Each is here
with the trace that forces it.

**1. Longest surface first.** `nonlinear` must beat the `linear` inside it. All
12 archived "not linear" answers score as "linear" otherwise — and the first
draft of the implementation did exactly that on 3 of 15.

**2. A fused negative is declared, not inferred.** `nonlinear` is a surface of
"not linear" in its own right; inferring a negation *on top of* the declaration
double-negates back to "linear". Support: 3 linearity traces.

**3. `neither A nor B` negates both slots.** Support: 2 of 24 memory-causality
traces, and it is the most dangerous construction in the archive — a matcher
seeing `memoryless` and `causal` without scoping the negation reads *"neither
memoryless nor causal"* as (Yes, Yes), the exact inversion, scored as a match.

**4. Two encodings for a label tuple, and gold uses the rarer one.** Gold writes
`Memoryless: **No**`; models write `Not memoryless` on 11 of 24 and
`Memoryless: No` on 6. Reading only gold's encoding scores 21 of 24 wrong.

**5. The origin marker is answer-bearing and usually absent.** 1 of 16 archived
traces emits `*value*`; 2 give an explicit index list; 1 gives per-element
assignments; the rest state no origin at all. And **gold omits it on 16.4% of
instances** (4,000 seeds) because the template marks `n = 0` only when the
result's support contains it (D-050).

**6. The arbitrary function of integration is retained by real models.** 2 of 6
continuity traces carry `C(y)` or `C(x)`; one of those two is a correct answer
that has additionally shown its working, and the other has dropped the solution
and kept only the constant. **The same surface form, opposite verdicts** — which
is why the policy distinguishes "annotated" from "instead of" rather than keying
on the presence of the symbol.

---

## 5. Two things I got wrong here

**The support counts in the module were invented.** The first draft wrote
"44/61" and "47/61" beside the two answer markers, from impression. The measured
counts are **57** and **55**. Numbers written from memory next to the rule they
describe are D-034's failure mode in miniature, committed in the file that argues
for D-034. They now live in a machine-checked table and
`derive_vocabulary.py` fails if they drift.

**The probe that caught it was itself wrong first.** The marker probes ran
against the *answer span* — from which `answer_span` has already consumed the
marker — so they reported 0/61 for markers the same file documented as appearing
on most traces. A self-contradiction visible in the output, which is the only
reason it was caught at all.

Both are in `phase4_summary.md` §8.
