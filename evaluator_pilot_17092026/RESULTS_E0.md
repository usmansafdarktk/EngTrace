# E0 baseline over the frozen pilot slice

Run 2026-09-18. 300 traces (60 frozen items x 5 models), the published framework
called unmodified through `evaluators/e0_tribunal.py`. **$6.28, 300 of 300 scored,
0 judge failures, 0 truncations, 0 errors.** Scores live in `scores/e0/`.

This is the anchor, not a result about the models: it is what every other evaluator
candidate is compared against once the human labels exist. Read it with the defects
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
human labels will settle.

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
