# Analysis plan for the full run

Written 2026-09-28, before any inference on the pool (D-117). It fixes which claims the paper will
test, how, and at what unit, so the paper can say the analysis was set before the data. The owner
decided its two scoring rules on 2026-09-28; the plan as a whole binds once the owner confirms it, and in
any case before the first inference call. Nothing here approves spend: each paid step still needs its own
estimate and approval.

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

**An unusable trace scores 0** (the owner, 2026-09-28). Not answering within the budget is the model's
failure, and excluding such traces would score each model on a different set of items. Three rules keep
that fair:

- **Defined by the answer, not the format.** A trace is unusable when the answer check can read no final
  answer: it is empty, cut off before answering, or never states one. A missing `**Answer:**` marker alone
  does not make a trace unusable: in the pilot GPT-5 used the exact marker in 87% of traces and Llama 3.1
  70B in 70% (FINDINGS R-F1). The gold validation confirms how the check reads such traces.
- **Service failures are not the model's.** The harness retries HTTP errors, timeouts and provider
  faults; one that persists is reported as missing and never scored.
- **Generous, stated token ceilings.** Each model's ceiling is set high and recorded in Appendix P, so
  running out means the model's reasoning ran away, not that the budget was tight.

Per model, the unusable rate is reported beside the score, and an appendix gives the score with unusable
traces excluded.
*Note, 2026-09-30 (D-148):* seven answered rows in the main run ended with a finish reason of `error` or
none, a provider fault reported inside a 200, and were scored on what they state (4 correct, 2 incorrect,
1 unusable); they are reported as a count, and the harness now retries such a reply as a service failure.
The unusable count is reported as empty plus unreadable.

**A partial answer scores 0.5** (the owner, 2026-09-28: not fully correct, but awarded some points). A
correct answer scores 1 and an incorrect or unusable one 0, so a model's **answer score** is its mean over
the items. Half a point uses only the three-way verdict, correct, partial or incorrect, which is what the
expert validation measured: agreement 0.893 three-way and 0.947 on the traces the experts judged fully
right or fully wrong (RESULTS_X1 Finding 1b). A fraction of parts correct would be finer, but it was never
compared with the experts. Only multipart items can be partial, 480 of 2,250. Reported beside the answer
score: the **fully-solved rate**, in which a partial scores 0, and each model's three-way split.
*Correction, 2026-09-30 (D-145):* the sentence "only multipart items can be partial" was wrong. The check
gives `partial` whenever some but not all of an item's targets match, and 531 items carry more than one
target: 352 multipart, 83 symbolic, 51 vector and 45 classification. In the store 278 of the roster's 514
partial verdicts fall on non-multipart items. Seven multipart templates have a single target on every
item, and 84 multipart items are checked on fewer quantities than the question's labelled parts; "correct"
means every quantity the check verifies. The rule is unchanged.

## Confirmatory questions

Each question names its test before the data. Within a family of comparisons, significance claims use
Holm's correction; everything is also reported as an estimate with its interval, and a null is stated as
"no difference as large as the detectable one", never as equality.

**Q1. How well does each model score, and which differences hold?** The answer score per model with a
template-level interval. Pairwise, 55 pairs: a template-level paired bootstrap interval of the score
difference and a sign-flip permutation test on the 150 per-template mean differences (10,000
permutations); a model is said to beat another only at a Holm-adjusted p below 0.05. The same comparison
on the fully-solved rate, with McNemar's exact test on the paired verdicts, checks that partial credit
does not decide an ordering.

**Q2. The complexity cliff.** Per model, the answer score on the 58 Easy templates minus the 34 Advanced
ones, templates resampled within each tier; Holm across the eleven models. With σ the standard deviation
of the template-mean score within a tier, the smallest gap this detects at 80% power and two-sided 0.05 is
2.80 x σ x sqrt(1/58 + 1/34) = 0.605 σ:

| σ | 0.20 | 0.25 | 0.30 |
|---|---|---|---|
| detectable Easy-to-Advanced gap | 12 points | 15 points | 18 points |

The paper reports each model's gap with its interval and this table; the May version's cliff claim is not
carried over unless it survives this test.
*Correction, 2026-09-30 (D-146):* the permutation test `analyze.py` added beside this interval (D-136)
shuffles tier labels on the raw difference, which is liberal when the smaller tier has the larger spread;
Advanced template means spread two to four times as widely as Easy ones, and the test's count of models
moved with its seed. The count now rests on Welch's t-test with Holm across the eleven models, with the
planned permutation printed beside it, and the detectable gap is also given at the strictest Holm step
from the Welch standard error. The gap is also shown with unusable rows left out of the template means
and without the nine symbolic templates. The intervals above are unchanged.

**Q3. What process scores add beyond the answer.** Descriptive, with template-level intervals and no
ranking claims: E5 milestone coverage on wrong-answer traces (how far a model gets before it fails); the
E4 digit rule's flag rate on correct-answer traces (arithmetic slips behind a right answer); the router's
step-error flags if it runs. Per model, the judged fraction and the unjudged-milestone rate, so the
reader sees how much of each score a judge decided.
Milestone coverage is not defined for an item with no milestones: 45 items, the whole of three
templates (`GOLD_VALIDATION.md`), are left out of milestone aggregates rather than scored 0, and 52
templates have some item with a single milestone, where coverage is close to an answer check.
*Correction, 2026-09-29 (D-135):* 70 items have no milestones: the 45 of those three templates and 25
more in 12 other templates. The 52 templates are those with an item of at most one milestone; 46 have an
item with exactly one. The rule is unchanged.
*Corrections and additions, 2026-09-30 (D-148, D-149):* the digit rule's flag rate on wrong answers is
over the answered wrong answers, because an empty trace has no step to flag; the first version counted
the empties. Reported beside the rates, descriptive: the 1% flag rate E4 shipped with; E3's chance floor
per model, the trace scored against a sibling item's milestones; where the first flagged step falls in a
fully solved trace; and the wrong-answer rate against the item's milestone count. A judged stage with
calls left unanswered is incomplete: its rates are over the answered calls and the count is printed.

**Q4. Consistency within a template.** Per model, the share of templates solved, meaning fully correct,
on all 15 instances, on some, and on none, reported separately for the single-path templates (one reasoning path across their 15
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
Per model: the paired difference in answer score, paraphrase minus original, with a sign-flip
permutation test over templates and a template-level interval, Holm across the models tested; McNemar's
exact test on the fully-solved verdicts as a check. Across models: Kendall's τ between their answer scores
on the originals and on the paraphrases, with a bootstrap interval. On the fully-solved verdict, the
detectable difference per model is 2.80 x sqrt(d / 450) for a discordance rate d, before clustering
widens it:

| discordance d | 5% | 10% | 15% |
|---|---|---|---|
| detectable paired difference | 3.0 points | 4.2 points | 5.1 points |

*Additions, 2026-09-30 (D-149):* beside the answer score, the paired difference in E3 milestone coverage
on the items with milestones, and in E5-strict once E5 has run on both arms (the coverage delta
NEXT_CYCLE_REVIEW section 6 asks for); and the answer-score difference on the pairs both arms served from
the same endpoint, since the arms run the same routing and each row records its provider. Kendall's tau is
read against its noise floor (below).

## Sensitivity analyses

- The answer check's relative tolerance at half and at double the fitted value.
- The fully-solved rate in place of the answer score, and unusable traces excluded rather than scored 0.
- Headline accuracy without the four templates the pilot excluded for shortcuts: two predictable from the
  question surface (D-057), one about 90% shortcuttable (D-046), one with a blind-guess floor of 1.0 (D-066).
- Without the two templates widened for round 4, if round 4 has not returned when the paper is written.
- *Added 2026-09-30 (D-147, D-149), after the first results were read:* without the nine symbolic
  templates, whose answers the check scores by the numbers they state (D-138); the answer check's
  half-unit window, a correct rounding at the precision shown where the rule accepts one unit either way;
  the whole-trace reading, which credits a quantity the question asks for when it is computed in the body
  and left off the Answer line; and a noise floor for Kendall's tau, the ordering on one half of each
  template's items against the other. The experts could not arbitrate the first two readings (no pilot
  verdict moves) and split the third one each way, so none is the score.

## Also reported, not tested

Branch and domain scores; scores by answer type; token use against score; the error attribution
from E4, E5 and the router; per-label accuracy on classification templates, noting that the pool
over-represents rare branches by design (D-116). A decoding-variance repeat: one cheap model, 300 items
(2 per template) three times at the main-run settings, reporting the spread of the score across repeats
(on the pricing document's basis about $1.40 for gemma-4-26b-a4b-it; needs approval).

## Not tested, and why

The math-pretraining claim of the May version: the roster has no math-specialised model (D-110), so the
paper does not make it.

## Rules for reporting

Every number in the paper is printed by a committed script reading the run's outputs. No claim beyond the
questions above; an unplanned analysis is labelled exploratory where it appears.
