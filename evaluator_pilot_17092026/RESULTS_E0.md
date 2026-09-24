# E0 baseline over the frozen pilot slice

Run 2026-09-18. 300 traces (60 frozen items x 5 models), the published framework
called unmodified through `evaluators/e0_tribunal.py`. **$6.28, 300 of 300 scored,
0 judge failures, 0 truncations, 0 errors.** Scores live in `scores/e0/`.

This is the anchor, not a result about the models: it is what every other evaluator
candidate is compared against once the expert labels exist. Read it with the defects
in [FINDINGS.md](FINDINGS.md) in view - particularly E0-F6, that the Tribunal is two
judges rather than three.

## Per trace model

| Model | Final-answer acc | Reasoning F1 | recall | precision | BERTScore | reached judges |
|---|---|---|---|---|---|---|
| claude-opus-4.7 | **0.700** | **0.465** | 0.599 | 0.397 | 0.880 | 47 of 60 |
| gpt-5 | 0.617 | 0.430 | 0.528 | 0.379 | 0.882 | 36 of 60 |
| gemini-3.1-pro | 0.633 | 0.418 | 0.576 | 0.346 | 0.873 | 37 of 60 |
| deepseek-r1 | 0.633 | 0.398 | 0.572 | 0.342 | 0.866 | 38 of 60 |
| llama-3.1-70b | 0.150 | 0.145 | 0.186 | 0.126 | 0.882 | 17 of 60 |

BERTScore separates nothing: 0.866 to 0.882 across models whose accuracy runs from
0.15 to 0.70. Llama scores 0.882, the joint highest, while answering 15% of items
correctly.

## Per answer type

| Answer type | Reasoning F1 | Final-answer acc | n |
|---|---|---|---|
| symbolic | 0.675 | 0.800 | 40 |
| array | 0.407 | 0.800 | 20 |
| classification | 0.391 | 0.475 | 40 |
| scalar | 0.368 | 0.500 | 60 |
| multipart | 0.280 | 0.425 | 80 |
| vector | 0.268 | 0.550 | 60 |

Multipart and vector answers score lowest, and they are 140 of the 300 traces. The
frozen slice was built to over-represent exactly these kinds because the comparator
is known to be weak on them; this is that choice paying off, and it is where E3/E4
have the most room to beat E0.

## What the judges actually said

2,310 votes across 175 judged traces:

| Category | Votes | Share |
|---|---|---|
| Alternative Correct | 2,156 | 93.3% |
| Calculation Error | 101 | 4.4% |
| Conceptual Error | 45 | 1.9% |
| Other | 8 | 0.3% |

The Tribunal says "Alternative Correct" to 93% of the steps Tier 1 rejected. Given
E0-F5 - that Tier 1 matches almost nothing and 93% of correct traces reach the
judges - E0's reasoning score is close to: *Tier 1 rejects nearly every step, then
the judges reinstate nearly every step.* Whether that is a sound mechanism or an
expensive way to approximate "did it get the answer right" is precisely what the
expert labels will settle.

**A hypothesis this data refutes.** If the judges reinstate nearly everything, F1
should be driven by the ratio of trace steps to gold steps rather than by
correctness. Measured across all 300 traces, the correlation between Reasoning-F1
and step-count agreement is **0.069** - nothing. Correlation with final-answer
correctness is **0.765**. Among the 175 judged traces alone, step-count agreement
does matter moderately (0.386), but the metric is not a step-count artefact.

## Cost, measured

$6.28 for 300 traces, $0.0359 per judged trace. Two judges: GPT-5 $0.0311 per call
(2,861 output tokens, mostly reasoning) and Claude $0.0191 (304 output tokens).
The pre-run estimate of $17 assumed three judges; the third was never billed
because it was never called (E0-F6).

Wall clock 3.9 hours, of which 53% was serial waiting on judge APIs and 47% was
the framework's post-judgement recovery path recomputing cross-encoder pairs one at
a time - pairs its own Tier 1 had already scored in a batch, and that the Kaggle GPU
dry run had already computed. Both are addressed before the next candidate runs:
cache `CROSS_ENCODER.predict` per pair (output-identical) and give the harness
`--workers`.

---

# E0-3J: the same framework with the third judge connected

Run 2026-09-18, 300 of 300 scored, **$7.82**, 0 failures. Every judged trace
reports three judges. `evaluators/e0_3j_tribunal.py` restores Google by supplying
the missing `get_model_info` the framework calls (E0-F6) - the framework's own
check then passes and appends `'google'` by its own logic. Nothing else differs.

## The third judge barely moves anything

Compared on the **160 traces both panels judged**:

| | |
|---|---|
| Traces whose Reasoning-F1 changed | **3 of 160** |
| Direction | 3 up, 0 down |
| Mean change on those 3 | +0.127 |

So the defect in E0-F6 is serious as a mechanism - the paper describes a panel that
never existed - but its effect on the scores is small. That is itself evidence for
E0-F5: a panel voting "Alternative Correct" 93% of the time does not discriminate
much, so a third vote of the same kind rarely changes the outcome.

**Directionally it can only help, and did.** With two judges the majority test
`count > total/2` demands unanimity, so any split falls to the conservative
`min()`. With three, a 2-1 split carries. None of the 3 changes went down.

## A caveat that matters more than the result

The whole-column averages look bigger than the like-for-like comparison:

| Model | e0 F1 | e0_3j F1 | delta |
|---|---|---|---|
| llama-3.1-70b | 0.145 | 0.176 | +0.030 |
| deepseek-r1 | 0.398 | 0.417 | +0.019 |
| gemini-3.1-pro | 0.418 | 0.432 | +0.014 |
| gpt-5 | 0.430 | 0.434 | +0.004 |
| claude-opus-4.7 | 0.465 | 0.461 | -0.004 |
| **all 300** | **0.371** | **0.384** | **+0.013** |

**Most of that +0.013 is not the third judge.** The harness seeds the framework's
20% wrong-answer sample from `(evaluator, item, model)`, so changing the evaluator
id changed which wrong answers were sampled: e0 judged 175 traces, e0_3j judged
179, and only 160 are common. 19 traces were judged only by the three-judge panel
and 15 only by the two-judge one. Reading the column deltas as the third judge's
effect would attribute a sampling difference to a mechanism.

**Fix for the remaining candidates:** seed the sample from `(item, model)` alone,
so every evaluator that uses this trigger sees the same sampled traces. E1 and E5
share E0's architecture and would otherwise each draw their own sample, making
their comparisons against E0 partly a comparison of samples. E3 and E4 are
deterministic and never sample, so they are unaffected.

## One more judge behaviour worth recording

Across e0_3j's 3,528 votes the Gemini judge returned a category the prompt never
defines - `"Standard processing applied"`, 3 times. The framework's mapping is
`"alternative" -> 1.0`, `"calculation" -> 0.5`, **everything else -> 0.0**, so an
invented category scores identically to "Conceptual Error". A judge that fails to
follow the label set is silently counted as judging the step wrong.

# The answer check, corrected (2026-09-24)

Every number above comes from the published framework's own final-answer check, which the
expert labels show is wrong on 72 of 300 traces - almost always calling a correct answer
wrong (68 of the 72). Its accuracy column therefore understates every model by about 21
points and ranks GPT-5 fourth of five.

`evaluators/answer.py` replaces it and agrees with the experts on 0.947 of the non-partial
traces against E0's 0.747; `analysis/answer_check.py` is the measurement, RESULTS_X1
Finding 1b the result. **E0 has not been re-scored with it** (D-098), so the figures in this
document stand as the published framework produces them, and the corrected accuracy is
reported separately. Re-running E0 would cost about $8.50 rather than the $6.53 spent here,
because a correct check sends roughly 233 traces to the judges instead of 178.
