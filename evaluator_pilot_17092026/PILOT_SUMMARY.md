# The EngTrace evaluator pilot

**September 2026. What we ran, what we learned, and how the full benchmark will be scored.**

The pilot asked one question: *how should EngTrace evaluate reasoning at full scale, and is
the published framework the right instrument?* Six candidate evaluators, and a variant of the
published one, were run over a frozen slice of the benchmark, scored against labels from 15
domain experts, and stress tested with defects we planted ourselves. This is the result.

**In short.** The published framework's final-answer check is wrong on a quarter of the traces,
and its reasoning score tracks that check closely (correlation 0.74). A deterministic
replacement agrees with the experts on 0.947 of the traces where the published check manages
0.747. For the reasoning itself, no candidate beats the published framework at the trace level
on this slice. Below the trace level, arithmetic slips behind a correct answer can be flagged
without a model: a digit-level check finds a third of them, and three of its flags in four are
real. A misstated rule behind a correct answer is caught only by a judge, and on planted
defects the best judges catch about a third. So the recommended stack is deterministic
wherever a check exists, and spends judge calls only on what no check can settle.

---

## 1. How the pilot was built

### 1.1 The frozen slice

60 problems drawn from five engineering branches and three difficulty levels: 15 templates
with four instances each. The slice is frozen with a recorded master seed and a manifest hash,
and every scoring run checks the problems against the manifest before it is allowed to score, so
two runs months apart are known to have read the same problems. Since the template
certification edits of 23 and 24 September, 17 of the 60 items no longer regenerate
byte-identically from the current templates, so the evaluators that derive milestones run
against the templates as they were at the freeze.

The slice was chosen to over-represent the answer types the published parser handles worst:
12 of its 60 items (20%) have a single scalar answer, against 89 of the benchmark's 150
templates (59%). That makes it a hard test of answer checking, and it means the size of any
answer-accuracy correction measured here does not carry over to the full benchmark.

Five models generated solutions, giving **300 traces**. Two smaller models were added as a
robustness cohort, checked mechanically but not labelled by the experts.

### 1.2 The candidates

| | what it computes | headline score |
|---|---|---|
| **E0** | the published framework: semantic step matching, then a panel of LLM judges on the steps it cannot match | share of reference steps recovered |
| **E0-3J** | a variant of E0 with its third judge connected, which in the published code was never actually called | same |
| **E1** | E0's method with all three judges replaced by families outside the evaluated suite | same |
| **E2** | three open process reward models scoring every step, run on university GPUs with no API cost | fraction of steps above the reward threshold, and the minimum step reward |
| **E3** | deterministic milestone verification: did the trace state the intermediate quantities the reference solution computes, order free and unit aware | milestone coverage |
| **E4** | E3 plus a check of the arithmetic claims the trace displays, where the checker can parse them (about a third of the equations) | coverage, plus arithmetic consistency |
| **E5** | E3 first, then one judge only on the milestones E3 cannot settle | strict milestone coverage |

All of them ran on the same 300 traces (the two Qwen PRMs skip 2 that exceed their context), resumable
and keyed by a configuration hash, so a change in scoring cannot be confused with a change in
data.

### 1.3 The expert labels

15 domain experts, three per branch. Every trace was labelled by three experts **of its own
branch**, each producing:

* a label for every step: correct, alternative correct, incorrect, or not a claim, with an
  error type and a written reason when incorrect;
* a status for every milestone the template defines;
* a verdict on the final answer: correct, partial, incorrect, or not stated;
* an overall verdict on whether the reasoning was sound.

Traces were blinded behind opaque codes, assignment was fixed in advance, and a shared
calibration set was labelled across branches so that agreement *between* branches could be
distinguished from agreement *within* one. 1,020 submissions in total, 900 of them forming the
ground truth; the annotation app recorded 13.4 hours with a trace open across them, not counting
the two later rounds.

### 1.4 Three quality rounds

Three rounds followed the first pass, and they changed the ground truth materially.

1. **A reason on every step called incorrect.** All 1,042 of them. Without this a
   disagreement cannot be adjudicated, only counted.
2. **Verification.** Each expert re-labelled seven of their own traces, blind to their first
   pass. This measures self-consistency, the ceiling against which between-expert agreement
   should be read.
3. **Adjudication.** Every step where a branch's three experts split 2 to 1 went back to all
   three, blind, with the three original labels shown anonymised and shuffled so no reviewer
   could tell which was their own. The majority of that second pass became the label. 272
   split steps across 134 traces.

Two of the fifteen label sets were re-annotated after verification identified them as least
self-consistent. Both originals were kept.

### 1.5 How the comparisons were computed

* **Intervals.** 2,000 bootstrap resamples on every comparison, with each evaluator's
  difference from the published framework computed on the same resamples, so a paired
  comparison is not reported as two independent ones.
* **Clustering.** The 300 traces are 15 templates with four instances each, and instances of
  a template are the same problem with different numbers. Every headline comparison is
  therefore also reported with **templates** resampled rather than traces.
* **Power.** From the clustered standard error, the smallest difference detectable at 80%
  power and 95% confidence, which is what turns a null result into a statement about what the
  design can see.
* **Baselines.** Two: the framework's own answer check, and the experts' answer verdict,
  which is the best any answer check could do.

### 1.6 Planted defects

Two objections follow any result drawn from the labels. The arithmetic rule we tested is the
same rule the annotation guide gave the experts, so their agreement is partly built in. And
conceptual error behind a correct answer cannot be studied on this corpus at all: the experts
found three such steps against 175 calculation slips.

So we built a set whose ground truth is true **by construction**. From the 129 traces the
experts called clean, 120 each receive exactly one defect and 60 are kept untouched as
controls (58 of those are also the source of a planted trace):

| family | n | what changes | what is held fixed |
|---|---|---|---|
| arithmetic | 60 | one digit of one displayed intermediate value | later steps and the final answer |
| conceptual | 60 | the reasoning a step states: a criterion flipped, a rule misstated, a stated formula contradicting the one actually used | **every digit in the trace** |
| control | 60 | nothing | everything |

Each plant is verified without the checker under test: an arithmetic plant is proved wrong by
recomputing the claim at full precision or against the template's own milestone value, and a
conceptual plant must leave the trace's digit string byte identical. The build is seeded and
reproduces exactly.

---

## 2. What we found

### 2.1 The final answer decides the trace-level verdict

![trace](figures/trace_level.png)

As a point estimate, the experts' own final-answer verdict predicts their soundness verdict
better than any evaluator: AUROC 0.974, against 0.850 for the published framework and 0.886 for
the best evaluator, E5. Of 228 traces with a correct final answer, their holistic verdict calls
only **three** unsound.

No evaluator beats the published framework at the trace level: E5 leads it by +0.036, well
inside the 0.19 this slice could detect (section 2.11), and E1 trails it by a statistically
distinguishable but negligible 0.002. On this slice, then, **the holistic verdict cannot rank
reasoning evaluators**: it is almost entirely decided by whether the answer is right. What
distinguishes the candidates has to be measured at the step and milestone level, which is where
the rest of the pilot works.

### 2.2 The published answer check is wrong on a quarter of traces

It disagrees with the experts on **72 of 300**, and almost always in one direction: 68 are
traces the experts call correct and it calls wrong. The corrected check below calls 56 of those
68 correct.

Five defects cause it, each measured. Where the reference answer line restates the question's
inputs, the first number on it is read as the answer: for one template that is the
temperature, not the molar volume. For another, the reference value read is the exponent inside
the unit `m^3/s`, so a trace was scored correct or wrong according to whether it wrote its units
in ASCII. Multi-part answers were cut to their last part. Superscripts, LaTeX and thousands
separators were read as the wrong numbers, or none. And a single relative tolerance cannot tell
a correct rounding from a wrong digit. Behind several of these, eight of the fifteen references
put no number on the answer line at all, two of them because the answer is a word.

| answer accuracy | experts | published framework | corrected check |
|---|---|---|---|
| strongest model (GPT-5) | 0.950 | **0.617** | 0.950 |
| all 300 traces | 0.760 | **0.547** | 0.753 |

**On this slice the published check understates answer accuracy by 21 points overall, and it
mis-ranks the top: it puts GPT-5 fourth of five, where the experts and the corrected check put
it first.** Per model the understatement is 33 points for GPT-5, 28 for DeepSeek R1, 23 for
Gemini 3.1 Pro, 22 for Claude Opus 4.7 and none for Llama 3.1 70B. The replacement is
deterministic and needs no model:

1. read the trace's answer segment from its final-answer heading rather than its last
   "Answer" marker, so a multi-part answer is not truncated to its last part;
2. take each target from a quantity the reference **computed**, not from any number it
   prints, which removes restated inputs and formula constants without per-template rules;
3. score every part separately and return **correct, partial or incorrect**, because 19 of
   the 300 traces are genuinely partial and a binary check cannot represent them;
4. accept a value within a relative tolerance **or** within one unit of the last digit it
   displays, or of the reference's, since `114` for 113.55 is a correct rounding while
   `50.20` for 50.05 claims two decimals and gets them wrong;
5. allow unit factors, and read superscripts, LaTeX and fractions as the numbers they are.

It agrees with the experts on **0.947** of the 281 traces they did not call partial, where the
published check manages 0.747, and on 0.893 of all 300 with partial scored as a verdict of its
own. Its one fitted parameter, the relative tolerance, was also fitted on each half of the
traces and scored on the other: held out, the three-way agreement is 0.876 and 0.905. 32
disagreements remain, 15 of them on one template whose question asks for six quantities while
its reference states one.

### 2.3 Which judges sit on the panel does not move the score

Replacing all three judges with families outside the evaluated suite moved the pooled score from
0.379 to 0.378, and no model's score by more than 0.004. Connecting the never-called third judge
moved it by +0.013, most of which proved to be a sampling artifact rather than the judge. So the
judge-and-judged objection finds no aggregate self-preference to point at: GPT-5's traces score
0.474 with GPT-5 on the panel and 0.470 without it. This is not a per-judge bias test against the
expert labels, and both comparisons were run under the published framework's own answer check,
so they compare panels with each other, not with the experts.

### 2.4 Milestone coverage with a residual judge is the most accurate milestone evaluator

| milestones, against the experts | precision | recall | F1 |
|---|---|---|---|
| deterministic matching (E3) | 0.926 | 0.915 | 0.921 |
| **deterministic, then a judge on the residue (E5)** | **0.930** | **0.989** | **0.958** |

Three quarters of milestones (76.4%) are settled with no model in the loop, and the judge's
"reached" verdict was validated before it was used: it credited none of 88 values the trace
never states. E5 cost $0.47 in judge calls over the 300 traces, plus $0.43 once to validate its
judge, where one run of the published framework cost $6.28 to $6.53. The published framework
produces no milestone score, so there is no like-for-like accuracy comparison at this level,
and at the trace level E5 does not beat it (section 2.1).

### 2.5 Milestone coverage cannot see a flawed step behind a correct answer

This is the case a reasoning evaluator exists for. The experts' step labels mark at least one
incorrect step in 93 of the 228 correct-answer traces (178 steps), although their holistic
verdict calls 90 of those 93 sound. Asked to find those 93:

| correct-answer traces, 93 of 228 flawed | AUROC | minus the published framework, templates resampled |
|---|---|---|
| published framework (E0) | 0.542 | — |
| best PRM, lowest step reward (E2) | 0.583 | +0.044 (−0.130 to +0.217) |
| milestones (E3) | 0.392 | −0.150 (−0.320 to +0.024) |
| milestones and a residual judge (E5) | 0.426 | −0.116 (−0.259 to +0.053) |
| **E4's arithmetic score, on the digit rule** | **0.661** | +0.149 (+0.000 to +0.300) |

The milestone evaluators score below chance here, because a trace with a slip usually still
reaches every milestone; with 15 templates that deficit is not significant, but its direction is
consistent. Only the arithmetic score separates from the published framework when traces are
resampled (it is scored on the 178 traces where the checker finds something to check), and with
templates resampled its interval touches zero. The stack's milestone layer therefore measures
progress through a derivation and says nothing about a flawed step behind a right answer; that
is what its arithmetic layer and its judge are for.

### 2.6 Process reward models over-flag, and the threshold is not the reason

The strongest of three open PRMs ranks steps well (AUROC 0.825), but only about half its flags
are real errors (precision 0.539), and inside correct-answer traces that falls to **0.246**. A
multi-domain PRM found 6% of the errors and was dropped.

Calibrating the flag threshold on one half of the traces and reporting on the other changes the
held-out F1 by **−0.006** (−0.018 and +0.006 in the two directions, both intervals straddling
zero): the default sits just below a flat plateau. Inside correct-answer traces no threshold
reaches precision 0.80 on held-out data, so the over-flagging is a property of the model rather
than of the cut-off.

### 2.7 Arithmetic errors can be flagged without a model

![steps](figures/step_level.png)

The experts' own rule is mechanical: *rounding is not an error, a wrong digit is*. Applied by
machine, it recomputes each claim from the numbers the trace itself shows and asks whether the
displayed value is a correct rounding at the precision shown.

E4 already had this checker and was asking it the wrong question. Its first tolerance, 1%
relative, answers "is this number fabricated", which on this corpus is rare, and it catches
almost none of the slips the experts marked: recall **0.034** on the steps inside correct-answer
traces.

Two corrections were forced by re-validating on the reference solutions, whose arithmetic is
correct by construction: compare a unit conversion after its unit factor, and do not hold a
result to a precision its own displayed operands cannot pin down. Together they took the
false-flag rate on reference solutions from 15.9% to **zero**. They cost recall: against the
experts' labels, precision rose from 0.506 to 0.750 and recall fell from 0.472 to 0.320.

On planted arithmetic defects inside a claim the checker can parse, the rule before those
corrections catches 45 of 45 and the rule as E4 ships it 40 of 45; the five it gives up are
last-digit slips of about 1e-7 to 3e-5 relative. Outside a parseable claim it catches 1 of 15.
So recall is bounded mainly by what the checker can parse, and partly by the precision the
corrections buy.

### 2.8 No deterministic check catches a wrong idea behind a right answer; a judge catches about a third

![planted](figures/planted.png)

This is the pilot's central negative result, and the planted defects are what made it
measurable rather than merely plausible.

Three judges were probed on the framework's own tribunal prompt, with exactly one step under
review. Each defect was judged **twice**: once in the planted trace and once in the same step of
the untouched original. A judge that flags both is flagging the step, not detecting the defect,
so only the difference counts.

| | digit rule (as E4 ships it) | GPT-5 | Claude Opus 4.5 | MiMo-V2.5-Pro |
|---|---|---|---|---|
| conceptual defects | **0 of 60** | 20 of 60 (0.333) | 8 of 60 (0.133) | 16 of 51 (0.314) |
| arithmetic defects | 41 of 60 (0.683) | 43 of 60 (0.717) | 43 of 60 (0.717) | 43 of 60 (0.717) |
| the same steps untouched, flagged | 3 of 120 | **0 of 120** | **0 of 120** | **0 of 111** |

Every deterministic evaluator the pilot has scores 0 of 60 on the conceptual defects, not only
the digit rule. MiMo returned a verdict on both arms for 51 of the 60 conceptual defects and for
111 of the 120 untouched steps: 17 of its 240 calls never returned, after retries.

Two things follow. A trace that computes correctly, misstates the rule it is applying, and lands
the right answer **passes every deterministic check we have**. And on this probe the judge that
satisfies the independence constraint detects about as much as GPT-5, the stronger of the two
frontier judges, so the full run does not have to trade sensitivity for independence. With about
60 defects per cell these rates are coarse: MiMo's 16 of 51 has a 95% interval of 0.20 to 0.45.

### 2.9 The framework shows a judge the flawed step as often as a clean one

![routing](figures/routing.png)

Being able to catch a defect is not the same as being given the chance. The framework's routing
was run for real on all 120 planted defects, with the judges stubbed so nothing was spent, and a
defect counts as caught end to end only when it was both shown and flagged:

| | conceptual | arithmetic |
|---|---|---|
| tribunal triggered | 0.650 | 0.833 |
| **the flawed step is shown to a judge** | **0.500** | 0.767 |
| *the same step, unmodified, is shown* | *0.500* | *0.683* |
| end to end: shown, and caught by either of its two judges | **0.150** | 0.600 |

For a conceptual defect the routing carries **no signal**: the corrupted step reaches a judge
exactly as often as the untouched one. It is forwarded because the matcher cannot match it, not
because anything is wrong with it. And the judges catch fewer of the defects they are shown (9 of
30) than of the whole set (21 of 60).

A verify-first router would send a judge every step no checker can verify. On the same planted
set it forwards 51 of the 60 conceptual defects; the other 9 sit in steps whose arithmetic the
checker verifies, which the router settles without a judge. With MiMo as its judge it catches 15
of 51 end to end (0.294; 0.267 with GPT-5 as its judge), about twice the published framework's
0.150. Those rates come from asking about one step per call.

A second finding fell out of the same run. The framework's answer check calls the final answer
wrong on 32 of the 120 planted traces, although the plants never touch the answer. Those traces
become eligible for its 20% wrong-answer sample, so part of its judge budget today is aimed by a
parser error.

### 2.10 The labels are reliable, and robust to replacing an expert

| | value |
|---|---|
| between experts, step labels (Fleiss kappa) | **0.781** |
| within an expert, blind re-label (Cohen kappa) | **0.828** |
| between experts, milestone status | 0.880 |
| between experts, final-answer verdict | 0.966 |

Pooled, self-consistency slightly exceeds between-expert agreement, which is the relationship
you want to see: the labels sit close to the ceiling the experts' own consistency sets. By branch
the between-expert step kappa runs from 0.645 (civil) to 0.867 (industrial), and in two branches,
the two added in September, the relationship inverts: civil's experts agree with their own first
pass at kappa 0.519, below the 0.645 they reach with each other, and industrial's at 0.820 against
0.867. Each within-expert figure rests on seven re-labelled traces per expert. Adjudication
settled every split and changed 85 step labels, raising the count of steps called incorrect from
322 to 388.

One stress test is worth repeating. **An entire expert's set, a fifteenth of the annotation
effort and the least self-consistent one, was replaced.** With the adjudication held fixed, no
ground-truth label moved; re-adjudicating the 12 disputes the replacement created then moved 12
of electrical's 406 step labels, all from correct to incorrect, and no trace verdict. The
adjudication design absorbed the weak set.

### 2.11 300 traces are worth far fewer independent ones

Resampling templates instead of traces gives a design effect of 2.6 to 4.2 depending on the
evaluator (3.7 for the published framework): the slice is worth roughly 70 to 115 independent
traces, not 300. That treats every instance of a template as the same problem; clustering by item
instead, the weaker assumption, gives a design effect of 1.5 to 1.8, and the true penalty lies
between the two. The smallest AUROC difference from the published framework this design can
detect is 0.12 to 0.19 for the evaluators built differently from it; only its own variants, which
track it trace by trace, are measured more finely (0.06 for E0-3J, 0.004 for E1).

Every null result in this pilot therefore means "no difference this large", not "no difference".
The dominance of the final answer clears the bar against every evaluator but the best: the
experts' answer verdict beats the published framework by +0.124 (template interval +0.043 to
+0.218), and beats every other evaluator with an interval that excludes zero except E5, where the
margin is +0.088 and the interval, −0.002 to +0.181, includes zero.

---

## 3. How the full benchmark will be scored

![cost](figures/cost.png)

### 3.1 The stack

| layer | how it works | model | judge cost, full run |
|---|---|---|---|
| answer correctness | three-way check, targets taken from computed quantities | none | **$0** |
| milestone coverage | deterministic, order free, unit aware | none | **$0** |
| arithmetic integrity | the digit rule, at displayed precision | none | **$0** |
| residual milestone judging | one judge, only on the milestones the deterministic pass cannot settle | MiMo-V2.5-Pro | **~$77** |
| step routing | verify what can be verified, send only the residue to a judge | MiMo-V2.5-Pro | **~$79 batched, ~$371 one step per call** |

With a batched router that is **about $156** in judge calls to evaluate the full benchmark, or
about $448 if every residue step is its own call. The published framework would cost about $334
on the same basis, and we are not running it at scale. The milestone judge and the published
framework are scaled from their measured per-trace judge spend on the pilot, pricing the eight
open-weight roster models at the pilot's weak model's rate and the three closed ones at its
frontier models' rate: for E5 that is the expensive end, since a weak model leaves more milestones
to the judge, and for the published framework it is the cheap end ($227 to $617 across the two
rates). The router is priced at the $0.0034 MiMo averaged per call in E5.

The router is in the stack because it is the only component aimed at the pilot's central
weakness. 63% of steps carry nothing a checker can recompute, 94% of traces hold at least one
such step, and 14% of that residue is a step the experts call incorrect. On the planted set it
lifts end-to-end conceptual detection from the published framework's 0.150 to 0.294 (section
2.9). Two things are not established. That detection was measured asking about one step per
call, which is the $371 design; the $79 figure assumes a trace's residue steps are batched into
one prompt, as the Tribunal batches its mismatched steps, and a batched prompt has not been
tested. And the router is not built.

Its routing rule also needs changing before it is built. As specified, a step counts as
verified as soon as any one of its claims can be recomputed, so a claim that passes can settle
a step whose error sits elsewhere. Of the 178 steps the experts call incorrect inside
correct-answer traces, 58 (a third) are never flagged by the digit rule and never reach a judge:
31 hold their error where the checker does not reach, such as a link inside a chained equality,
and 27 are slips the shipped digit rule declines to flag because the displayed operands cannot
pin the result down. The router has to send a judge every step the checker has not positively
cleared, which makes its residue larger than the 63% priced above.

### 3.2 The models

**The benchmark.** Five branches × 30 templates × 15 instances = **2,250 problems** (870 easy,
870 intermediate, 510 advanced), answered once by each model on the roster.

**The roster, 11 models, about $403 to generate**, on the inference pricing document's basis:
every model writes as much per problem as GPT-5 did on the pilot, 5,232 output tokens with
reasoning included. Reasoning-heavy models can exceed that; DeepSeek R1's median completion on the
pilot was 9,658 tokens.

| open weight ($212.81) | closed weight ($190.18) |
|---|---|
| gpt-oss-20b, gemma-4-26b-a4b-it, deepseek-v4.1-flash, qwen3-235b-a22b, glm-5.3-flash, glm-5.3, muse-glimmer-30b, kimi-k3 | gpt-5.4-mini, gemini-3.1-flash-lite, claude-sonnet-5 |

Generation and judging together come to about $559 with a batched router.

No model that judged in any of the pilot's evaluators is on the roster, so nothing being evaluated
has also been an instrument of the evaluation. The five pilot trace models and the robustness pair
are not on it either. Kimi K3 and GLM-5.3 were screened as judge candidates and deliberately not
chosen, to keep them free to be evaluated.

**The judge is MiMo-V2.5-Pro**, for two reasons the pilot separated. It is independent: no family
on the roster, which is what the judge-and-judged objection requires. And on planted defects it
detects about as much as GPT-5 (0.314 against 0.333 on conceptual defects, the same 0.717 on
arithmetic) with no false alarm, at about a third of a cent per call. It is not perfectly
reliable: 17 of its 240 probe calls never returned, after retries, so the full run needs a retry
and fallback policy.

### 3.3 What is not being run, and why

| | why not |
|---|---|
| the published tribunal as the reasoning metric | its answer check is wrong on a quarter of traces, and its reasoning score tracks that check (correlation 0.74) |
| the third judge | moves the score by +0.013, most of it a sampling artifact |
| the judge swap as a main arm | gave the same pooled score, 0.379 against 0.378; keep it as a robustness check on a sample |
| the multi-domain PRM | finds 6% of the errors |
| threshold tuning for the PRM | held-out change in F1 of −0.006 |
| the published framework at full scale | about $334, and unlabelled full-run numbers would add cost without adding evidence, since the evaluator comparison belongs to the 300 labelled traces |

The strongest process reward model can still be run on HiPerGator at no API cost as a secondary
per-step signal, about 14 GPU-hours for the full benchmark, but not as a headline number.

### 3.4 Two notes on reporting

**Intervals.** The full benchmark has 150 templates against the pilot's 15, so the same
cluster-robust treatment applies with ten times the clusters, and the comparisons the design can
support should be settled before the run rather than after.

**The table will move.** On this slice the corrected answer check raises accuracy by 3 to 33
points per model (21 overall) and puts GPT-5 first instead of fourth. The slice over-represents
the answer types the published parser misreads, so the size of the shift on the full benchmark
has to be measured when the table is regenerated, not carried over from here. Readers will
compare against the published version, so the paper needs a sentence explaining the change.

---

## 4. Limits worth stating plainly

* **The slice is 15 templates**, five models, one of them weak. It separates frontier from
  weak; it cannot rank frontier models, and no evaluator claim should be read as one.
* **The stack was validated on the pilot's five models, four of them frontier.** The roster is
  mostly smaller open models. Trace structure holds on two such models (Gemma-4-31B and
  Qwen3.8-27B), but the stack's accuracy against experts on them is not measured.
* **The hard-case finding is about arithmetic.** 175 of the 178 flawed steps behind a correct
  answer are calculation slips. Three are conceptual, which is why the planted set exists.
* **Planted defects are a diagnostic.** They say what an evaluator catches when a defect of a
  given shape is present, not how often models produce one. The conceptual defects vary in kind
  but not in phrasing, so they are not a held-out test.
* **The digit rule shares its rule with the annotation guide**, so its agreement with the
  experts is partly by construction. The planted set is the independent check.
* **The judges were asked about one step at a time.** A batched router would put several into
  one prompt: more context per call, but also more steps competing for attention, and that is not
  measured.
* **Expert labels, traces and scores are not in the repository.** The labels are kept outside by
  decision and the traces and evaluator scores are not committed, so the committed scripts
  reproduce the numbers only with those files supplied. Regenerating the traces would not
  reproduce them, since model outputs vary between calls.

---

## 5. Open items

1. **The stack has run on 15 of the 150 templates.** Before any full-run trace exists, run the
   answer check, E3 and E4 on the reference solutions of all 150, and check milestone counts per
   template: a template whose reference states one intermediate reduces E3 to an answer check.
2. **Template versions.** 17 of the 60 frozen items no longer regenerate byte-identically from
   the current templates; `pinned_templates.py` runs the pilot against the templates as they were
   at the freeze. The full run needs the same guarantee: freeze the 2,250-item pool with per-item
   hashes, and score against the template commit it was generated from.
3. **The router needs designing and building.** Its routing rule has to forward the steps
   the checker has not positively cleared (section 3.1), its output has to be defined as a
   reported score, and its batched prompt has to be measured, since the detection rate above
   comes from single-step prompts.
4. **MiMo's non-returns.** The full run needs a retry count, a fallback (an unjudged milestone
   counts as missing under strict coverage), and the unjudged rate reported per model.
5. **An anonymised ground-truth file**, together with the traces and evaluator scores (which
   hold no annotator data), would let the analyses reproduce without exposing annotator data.
   Pending a compliance decision.

---

*Every number here is printed by a committed script under `analysis/`, whose README lists which
script reports what; the figures are drawn from `figures/summary_numbers.json`, which
`analysis/summary_numbers.py` computes. Every decision the summary refers to is recorded with its
evidence in the project's decision log.*
