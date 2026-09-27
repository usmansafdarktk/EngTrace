# The EngTrace evaluator pilot — what we did and what we found

September 2026. The pilot asked one question: **how should EngTrace evaluate reasoning at
full scale, and is the published framework the right instrument?** It ran six candidate
evaluators over a frozen slice of the benchmark, scored them against expert labels, and
produced the design the full run will use.

---

## 1. How it was run

**The slice.** 60 problems — five engineering branches, three difficulty levels, 15
templates with four instances each — frozen with a recorded seed and hash so every later
run scores exactly the same bytes. Five models generated solutions, giving **300 traces**,
plus two smaller models as a robustness cohort.

**The candidates.** Six evaluators, each scoring the same 300 traces:

| | what it does |
|---|---|
| **E0** | the published framework: step matching, then a panel of LLM judges on what it cannot match |
| **E0-3J** | E0 with its third judge connected (it had never actually been called) |
| **E1** | E0's method with all three judges swapped for families outside the evaluated suite |
| **E2** | open process reward models scoring every step, run on local GPUs |
| **E3** | deterministic milestone verification — did the trace reach the quantities a correct derivation must reach |
| **E4** | E3 plus a check of the arithmetic each trace displays |
| **E5** | E3 first, then a judge only on the milestones E3 cannot settle |

**The expert labels.** 15 domain experts, three per branch. Every trace was labelled by
three experts of its own branch: a label for each step, each milestone, the final answer,
and an overall verdict on whether the reasoning was sound. Traces were blinded and a shared
calibration set was labelled across branches. That is 1,020 submissions, about 13 hours of
expert time.

**Three quality rounds, which turned out to matter more than the labelling itself:**

- **Reasons.** Every step called incorrect carries a written reason. All 1,042 of them.
- **Verification.** Each expert re-labelled seven of their own traces, blind to their first
  pass, so self-consistency is measured rather than assumed.
- **Adjudication.** Every step where a branch's three experts split 2–1 went back to all
  three, blind and shuffled, and the majority of that second pass became the label. 272
  split steps, reviewed in 134 traces.

Two of the fifteen label sets were re-annotated after the verification round identified them
as least self-consistent; both originals were kept.

**The analysis.** Bootstrap confidence intervals on every comparison, two baselines, and —
after we found that 300 traces are really 15 templates — intervals that resample templates
rather than traces. Plus a set of **planted defects**: known errors injected into traces the
experts had called clean, so detection could be measured against a truth we set rather than
against the annotation guide.

---

## 2. What we found

### The final answer decides almost everything

The experts' own final-answer verdict predicts their soundness verdict better than any
evaluator does (AUROC 0.974 against the best evaluator's 0.886). Of 228 traces with a
correct final answer, their holistic verdict calls only **three** unsound.

The consequence is structural: **a trace-level score cannot rank reasoning evaluators**,
because it mostly measures answer correctness. Everything interesting happens at the step and
milestone level, and that is where the pilot did its work.

### The published framework's answer check is wrong on a quarter of traces

It disagrees with the experts on **72 of 300** traces, and almost always in one direction —
68 are traces the experts call correct and it calls wrong. In 45 of those, the correct value
is sitting in the trace's own answer line.

Two defects cause it: the reference answer line often restates the question's inputs, and a
third of the templates answer with a word rather than a number. A third, found later, was
worse: the reference value was read as "the last number in the solution", which for one
template is the exponent inside the unit `m^3/s`.

| answer accuracy | experts | published framework | corrected check |
|---|---|---|---|
| best model | 0.950 | **0.617** | 0.950 |
| all 300 traces | 0.760 | **0.547** | 0.753 |

**The benchmark understates every model by about 21 points and mis-ranks the top.** A
deterministic replacement — read the answer segment, take each target from a quantity the
reference actually computed, score every part of a multi-part answer, return correct /
partial / incorrect — agrees with the experts on **0.947** of traces where the framework
manages 0.747. Its one fitted parameter was chosen on one half of the traces and reported on
the other.

### Swapping the judges changes nothing; the judges themselves are not the cost driver

Replacing all three judges with families outside the evaluated suite moved the pooled score
from 0.379 to 0.378. Connecting the never-called third judge moved it +0.013, most of which
turned out to be a sampling artifact. **Judge identity does not buy accuracy** — which
answers the independence objection and also means the panel's cost is not buying
discrimination.

### Milestone coverage with a residual judge is the best reasoning metric, and the cheapest

| milestones, against the experts | precision | recall | F1 |
|---|---|---|---|
| deterministic matching | 0.926 | 0.915 | 0.921 |
| **deterministic, then a judge on what it misses** | **0.930** | **0.989** | **0.958** |

Three quarters of milestones are settled with no model involved. Measured cost over the same
300 traces: **$0.47 against the published framework's $6.53** — a thirteenth, and more
accurate.

### Process reward models are not the answer, and their threshold is not why

The best of three open PRMs ranks steps well (AUROC 0.825) but only about half its flags are
real errors, and inside correct-answer traces — where a step evaluator would earn its keep —
that falls to **0.246**. A multi-domain PRM found 7% of the errors and was dropped.
Calibrating the flag threshold on half the labels and reporting on the other half gains
**−0.006**: the default sits on a flat optimum, and over-flagging cannot be tuned away.

### Arithmetic errors are deterministically detectable, and the framework was looking for the wrong thing

The experts' own rule is mechanical — *rounding is not an error, a wrong digit is*. Applied
by machine, recomputing each claim from the numbers the trace itself shows:

| inside correct-answer traces | precision | recall |
|---|---|---|
| the arithmetic check as the framework shipped it (1% tolerance) | 0.154 | 0.034 |
| **the digit rule** | **0.750** | 0.320 |
| best process reward model, for comparison | 0.246 | 0.255 |

**Three flags in four are real errors, for no cost.** The shipped 1% tolerance answers "is
this number fabricated", which on this corpus is almost never; it is blind to every slip the
experts actually marked. On planted arithmetic defects the digit rule scores **1.000 where it
can parse the claim and 0.000 where it cannot** — so its recall ceiling is a parsing limit,
not a rule limit.

### Nothing automatic catches a wrong idea behind a right answer — except a judge, a third of the time

This is the pilot's central negative result, and the planted defects made it measurable. We
injected 60 arithmetic defects, 60 conceptual defects (a criterion flipped, a rule misstated,
a stated formula contradicting the one used — with no digit in the trace changed) and kept 60
traces untouched as controls.

| | deterministic evaluators | judges, asked directly |
|---|---|---|
| conceptual defects caught | **0 of 60** | 0.13 – 0.33 (two of three ≈ 0.32) |
| arithmetic defects caught | 0.68 – 0.75 | 0.72 |
| false alarms on untouched traces | low | **0.000** |

Three judges were probed, each defect judged both in the planted trace and in the same step
of the original, so only the difference counts. **Not one false alarm on any untouched
step, from any of the three.**
The judge that satisfies the independence constraint performs like the frontier judges, so
the full run does not have to choose between independence and sensitivity.

But the framework rarely shows them the flawed step. Running its routing for real: the
corrupted step reaches a judge **0.500** of the time for conceptual defects — *exactly as
often as the untouched step*. It is forwarded because the matcher cannot match it, not
because anything is wrong with it. End to end, the framework catches about **18%** of
conceptual defects and 61% of arithmetic ones.

### The labels are reliable, and robust to replacing an expert

| | |
|---|---|
| between experts, step labels (Fleiss κ) | **0.781** |
| within an expert, blind re-label (Cohen κ) | **0.828** |
| final-answer verdict, between experts | 0.966 |

Self-consistency slightly exceeds between-expert agreement, which is the relationship you
want: the labels are close to the ceiling the experts' own consistency sets. Adjudication
settled every split and changed 85 step labels, raising the count of steps called incorrect
from 322 to 388.

One stress test is worth repeating: **an entire expert's set — a fifteenth of the annotation
effort, and the least self-consistent one — was replaced, and not one ground-truth label
moved.** Re-adjudicating the new disputes then changed 12 of 2,091. The adjudication design
is what absorbs a weak rater.

### 300 traces are really 15 templates

Instances of a template are the same problem with different numbers. Resampling templates
instead of traces, the design effect is about 3.7: the slice is worth roughly **80
independent traces, not 300**, and can detect AUROC differences of about **0.12–0.19** and no
smaller. Every null result here means "no difference this large", not "no difference" — and
the one finding that clears the bar comfortably is the dominance of the final answer.

---

## 3. What the full run will do

| layer | how | cost over ~25,000 traces |
|---|---|---|
| answer correctness | deterministic, three-way correct / partial / incorrect | **$0** |
| milestone coverage | deterministic, derived from the templates | **$0** |
| arithmetic integrity | the digit rule | **$0** |
| residual judging | one independent judge, only on milestones the deterministic pass cannot settle | **~$77** |

**The judge is MiMo-V2.5-Pro.** It was chosen for independence — it is not from a family on
the evaluated roster, which is what the judge/judged objection demands — and the pilot then
showed it does not cost anything to get that: on planted defects it matches the frontier
judges (0.314 against 0.333 on conceptual defects, identical 0.717 on arithmetic) with no
false alarms. It is also the cheapest of the three at about a third of a cent per call. The
deterministic layers use no model at all.

**The benchmark.** Five branches x 30 templates x 15 instances = **2,250 problems** (870
easy, 870 intermediate, 510 advanced), answered once by each model on the roster.

**The roster — 11 models, about $403 to generate:**

| open-weight ($212.81) | frontier ($190.18) |
|---|---|
| gpt-oss-20b · gemma-4-26b-a4b-it · deepseek-v4.1-flash · qwen3-235b-a22b · glm-5.3-flash · glm-5.3 · muse-glimmer-30b · kimi-k3 | gpt-5.4-mini · gemini-3.1-flash-lite · claude-sonnet-5 |

Every model the pilot generated traces with, or used as a judge, is deliberately excluded
from the roster, so nothing on it has been either a subject or an instrument of the
evaluation. Costs assume one pass per problem at the pilot's measured token profile.

**About $77 to evaluate the full benchmark**, against roughly $334 for the published
framework — which we are not running at scale, because the evaluator comparison belongs to
the 300 labelled traces where ground truth exists, and unlabelled full-run numbers would add
cost without adding evidence.

**Under consideration (~$79):** a router — verify deterministically what can be verified, and
send only the residue to a judge. 63% of steps carry nothing a checker can recompute, and 14%
of that residue is a genuinely incorrect step. It would roughly double end-to-end conceptual
detection, from 18% to about 31%.

**Not running:** the published tribunal as the reasoning metric (its score is dominated by
answer correctness and its parser is broken); the third judge (adds nothing); the judge swap
as a main arm (already shown equivalent); the multi-domain PRM (failed); threshold tuning
(gains nothing). The strongest process reward model, Qwen2.5-Math-PRM-72B, can still be run
on local GPUs for no API cost as a secondary per-step signal — about 15 GPU-hours for the
full benchmark — but not as a headline number, since it flags one real error in four inside
correct-answer traces.

**Two design notes.** The full benchmark has 150 templates against the pilot's 15, so
intervals should again be cluster-robust by template, with the supportable comparisons
settled before the run rather than after. And when the results table is regenerated with the
corrected answer check, **every model's accuracy will rise by roughly 21 points and the
ranking will change** — that needs a sentence in the paper, since readers will compare
against the published version.

---

## 4. Limits worth stating plainly

- **The slice is 15 templates**, five models and one weak model. It separates frontier from
  weak; it cannot rank frontier models, and no evaluator claim should be read as one.
- **The hard-case result is about arithmetic.** 175 of the 178 flawed steps behind a correct
  answer are calculation slips. Three are conceptual — too few to study, which is why the
  planted set exists.
- **Planted defects are a diagnostic.** They say what an evaluator catches when a defect of a
  given shape is present, not how often models produce one.
- **The digit rule shares its rule with the annotation guide**, so its agreement with the
  experts is partly by construction. The planted set is the independent check, and there it
  scores 1.000 inside claims it can parse.
- **The judges were asked about one step at a time.** A router would batch several into one
  prompt; that is untested.
- **Expert labels are held outside the repository** by decision, so the analyses are
  reproducible from the committed scripts only with the label files supplied.

---

## 5. Open items

1. **Template drift.** Three of the four metrics derive milestones from the benchmark's
   templates, and a subset of the frozen items no longer reproduces byte-identically. This
   must be fixed or formally pinned before the full run.
2. **The router decision** — measured, costed, not yet taken.
3. **An anonymised ground-truth file** would make the analyses reproducible without exposing
   annotator data. Pending a compliance decision.
4. **Two small follow-ups** from the label review: one probe step needs relabelling, and the
   annotation guide should state explicitly how to treat a notational slip whose computed
   result is right.

---

*Every number in this summary is produced by a committed script, and every decision it
refers to is recorded with its evidence in the project's decision log.*
