# Analysis plan for the full run

Written 2026-09-28, before any inference on the pool (D-117). It fixes which claims the paper will
test, how, and at what unit, so the paper can say the analysis was set before the data. Items marked
**proposed** are the owner's to confirm; the plan binds once confirmed, and in any case before the first
inference call. Nothing here approves spend: each paid step still needs its own estimate and approval.

## What is analysed

| | |
|---|---|
| items | the frozen pool, 2,250 = 150 templates x 15 (D-116); its manifest SHA-256 is in `FREEZE.json` |
| models | the eleven of D-110 |
| completions | one per item and model, the Appendix P prompt byte-identical and hashed; the harness records the served model id and the decoding settings of every call |
| evaluation | the answer check (correct / partial / incorrect), E3 milestones, E4 digit rule, E5 with MiMo-V2.5-Pro on residual milestones (D-105, D-110); the batched step router if funded (D-113) |

**Unit of analysis: the template.** Every template has exactly 15 items, so the item mean equals the
mean of template means. Every interval resamples templates, not items (B = 10,000, percentile), because
instances of a template share a derivation: the pilot's design effect was 2.6 to 4.2 (D-111).

**Unusable traces** (empty, stopped at the token ceiling, or no answer segment): **proposed**, counted as
incorrect in the headline, since not answering within the budget is the model's failure, and reported
per model beside it. Sensitivity: the same figures with unusable traces excluded.

**Partial answers** (a multipart answer with some parts right): **proposed**, not correct in the headline
and reported separately; sensitivity with a partial counted as 0.5.

## Confirmatory questions

Each question names its test before the data. Within a family of comparisons, significance claims use
Holm's correction; everything is also reported as an estimate with its interval, and a null is stated as
"no difference as large as the detectable one", never as equality.

**Q1. How accurate is each model, and which differences hold?** Accuracy per model with a template-level
interval. Pairwise, 55 pairs: McNemar's exact test on the paired correct/incorrect verdicts and a
template-level bootstrap interval of the difference; a model is said to beat another only at a Holm-
adjusted p below 0.05.

**Q2. The complexity cliff.** Per model, accuracy on the 58 Easy templates minus the 34 Advanced ones,
templates resampled within each tier; Holm across the eleven models. With σ the standard deviation of
template-mean accuracy within a tier, the smallest gap this detects at 80% power and two-sided 0.05 is
2.80 x σ x sqrt(1/58 + 1/34) = 0.605 σ:

| σ | 0.20 | 0.25 | 0.30 |
|---|---|---|---|
| detectable Easy-to-Advanced gap | 12 points | 15 points | 18 points |

The paper reports each model's gap with its interval and this table; the May version's cliff claim is not
carried over unless it survives this test.

**Q3. What process scores add beyond the answer.** Descriptive, with template-level intervals and no
ranking claims: E5 milestone coverage on wrong-answer traces (how far a model gets before it fails); the
E4 digit rule's flag rate on correct-answer traces (arithmetic slips behind a right answer); the router's
step-error flags if it runs. Per model, the judged fraction and the unjudged-milestone rate, so the
reader sees how much of each score a judge decided.

**Q4. Consistency within a template.** Per model, the share of templates solved on all 15 instances, on
some, and on none, reported separately for the single-path templates (one reasoning path across their 15
pool instances by `diversity.py`'s lower reading: 58 templates) and the rest. A single-path template
answered on some instances but not all was missed on the numbers, not the method. Descriptive, with
template-level intervals.

**Q5. Paraphrase robustness** (run after the main inference; its spend needs approval). Subsample: 450
items, the 1st, 6th and 11th of each template's 15 in manifest order. One paraphrase per item by a model
from a family on neither the roster nor the judge's side, keeping every number, unit, symbol, technical
term and the order of the parts; a script rejects a paraphrase whose numbers or units changed and a
near-copy; one own-branch expert confirms each is the same problem with the same answer, and a rejected
pair is dropped from both arms. The same models, prompt, settings and scoring as the main run, run in the
same week so the served models are the same.
Per model: the paired difference, paraphrase minus original, with McNemar's exact test and a template-
level interval, Holm across the models tested. Across models: Kendall's τ between their accuracies on the
originals and on the paraphrases, with a bootstrap interval. Per model the detectable difference is
2.80 x sqrt(d / 450) for a discordance rate d, before clustering widens it:

| discordance d | 5% | 10% | 15% |
|---|---|---|---|
| detectable paired difference | 3.0 points | 4.2 points | 5.1 points |

## Sensitivity analyses

- The answer check's relative tolerance at half and at double the fitted value.
- Headline accuracy without the four templates the pilot excluded for shortcuts: two predictable from the
  question surface (D-057), one about 90% shortcuttable (D-046), one with a blind-guess floor of 1.0 (D-066).
- Without `critical_depth_froude_classification`, whose class goes unscored (D-116).
- Without the two templates widened for round 4, if round 4 has not returned when the paper is written.

## Also reported, not tested

Branch and domain accuracy; accuracy by answer type; token use against accuracy; the error attribution
from E4, E5 and the router; per-label accuracy on classification templates, noting that the pool
over-represents rare branches by design (D-116). A decoding-variance repeat: one cheap model, 300 items
(2 per template) three times at the main-run settings, reporting the spread of accuracy across repeats
(on the pricing document's basis about $1.40 for gemma-4-26b-a4b-it; needs approval).

## Not tested, and why

The math-pretraining claim of the May version: the roster has no math-specialised model (D-110), so the
paper does not make it.

## Rules for reporting

Every number in the paper is printed by a committed script reading the run's outputs. No claim beyond the
questions above; an unplanned analysis is labelled exploratory where it appears.
