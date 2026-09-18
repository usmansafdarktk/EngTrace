# Findings from the evaluator pilot

Things the pilot surfaced about the **published** framework and its environment, as
opposed to defects in the pilot's own tooling (those are in the commit history).
Each one is measured, with how to reproduce it. None has been fixed in E0: E0 is
the baseline, and correcting it would move the anchor every candidate is compared
against. They are recorded so the comparison can be read correctly.

---

## E0-F1 · The final-answer parser reads a restated input as the gold answer

**Template:** `rackett_equation_volume` — 4 of the 60 pilot items.

The gold line is

    **Answer:** The estimated molar volume of saturated liquid n-Pentane at 283.81 K is **113.55 cm³/mol**.

`evaluation/engineering_parser.py::extract_final_answer_eng` takes the **first**
number after `**Answer:**`, so the gold final answer is read as `283.81` — the
temperature, an input — instead of `113.55`. On `rackett_equation_volume#0` Claude
Opus 4.7 answers `113.55 cm³/mol`, exactly right, and E0 scores
`final_answer_acc = 0`.

**This is not a rework regression.** The answer wording is identical at the
template's first commit (`645fb54`, 2025-09-15), at `6286022` (2025-10-14) and at
HEAD. So in the published results, final-answer accuracy on this template was 0
for every model whether it was right or wrong, unless a model's own answer line
happened to lead with the temperature too.

Reproduce: `extract_steps(item['solution'])[2]` on any `rackett_equation_volume`
item returns the temperature.

## E0-F2 · 8 of 15 pilot templates have no number on the Answer line

`reynolds_number_flow_regime`, `damping_classification` (classification),
`phasor_addition`, `incompressible_continuity` (symbolic), `lorentz_force`,
`fluid_particle_acceleration` (vector), `gas_phase_concentration`,
`coaxial_capacitance` (multipart).

With no number on the Answer line the parser falls back to the last `= number` in
the whole solution, so E0's "final answer" for these items is whichever quantity
happens to be computed last. For a classification or a symbolic answer that number
is not the answer at all. Final-answer accuracy is a scalar-only mechanism applied
to a corpus where scalars are 12 of the 60 pilot items; this is the gap E3/E4 are
meant to be tested against, and it should be read into E0's column rather than
discovered in it.

## E0-F3 · Under transformers 5.x, BERTScore is silently 0.0 on every entry

`bert_score` raises `OverflowError: int too big to convert` because the Longformer
tokenizer reports `model_max_length = 1e30`. `safe_bert_score` catches it, prints a
warning and returns `0.0`, so the column fills with zeros that look like scores.
`requirements.txt` leaves `transformers` unpinned, so a fresh `pip install -r` today
gets 5.x and reproduces this.

The pilot venv pins `transformers==4.57.3`, `sentence-transformers` 5.1.x and
`bert-score==0.3.13` — the line current when the published run was made, which its
Opus 4.5 judge dates to about December 2025. Verified afterwards: BERTScore 0.948 on
a paraphrase pair, and the cross-encoder separates same-meaning (0.923) from
unrelated (0.009). The harness now flags a BERTScore failure on the row instead of
accepting the zero.

## E0-F4 · One of the three published judges no longer exists

`gemini-3-pro-preview` returns 404 from Google: "no longer available, use
gemini-3.1-pro-preview" (probed 2026-09-17). It is also absent from OpenRouter's
catalogue. The published E0 cannot be re-run exactly by anyone; the pilot uses
Google's named successor and records the substitution on every scored row.
Precisely: `models.get` on the id still returns its metadata, `models.list` no
longer includes it, and `generateContent` returns 404. A get-by-name check is
therefore not an availability check. Re-checked 2026-09-17 16:33 UTC: Gemini 3.1 Pro
Preview is the only Gemini 3.x Pro text model Google lists.

## E0-F5 · Tier 1 matches almost no steps, so E0's reasoning score is the Tribunal's

From the full E0 dry run over all 300 traces (Tier 1 real, judges not called):

| Trace model | final-answer acc | Tier 1 F1, before the Tribunal | correct answers sent to judges |
|---|---|---|---|
| claude-opus-4.7 | 0.70 | 0.102 | 42 of 42 |
| gemini-3.1-pro | 0.63 | 0.175 | 34 of 38 |
| deepseek-r1 | 0.63 | 0.147 | 36 of 38 |
| gpt-5 | 0.62 | 0.170 | 34 of 37 |
| llama-3.1-70b | 0.15 | 0.058 | 6 of 9 |

Tier 1 accepts a step only if its number is within 2% **and** the cross-encoder
scores it ≥ 0.70 against a gold step. On this corpus that conjunction almost never
holds: Tier 1 F1 is 0.06-0.18. So 152 of 164 correct answers (93%) fall below the
80% alignment trigger and go to the judges, and for those traces the reasoning
score E0 reports is, in effect, the Tribunal's. Tier 1 is described as the
automated backbone with the LLMs as a fallback; measured here, the fallback is the
main path. That matters for every comparison against E0, and it is the strongest
argument that E1 (judge-family removal) and E3 (deterministic milestones) are
testing the part of E0 that actually decides scores.

Also confirmed at scale: **E0-F1 holds for all five models** - `rackett_equation_volume`
final-answer accuracy is 0 of 20 across every model, right or wrong. The
wrong-answer Tribunal sample realised at 24 of 136 (17.6%) against the framework's
20% rate, with the draw seeded per (item, model).

## Provenance · GPU and CPU Tier 1 agree exactly

The dry run was computed on a Kaggle Tesla T4 (torch 2.10.0+cu128, 769 s for all
300) and imported only after comparison with 159 traces independently scored on
the laptop CPU (torch 2.14.0+cpu), under identical pinned scoring libraries: Tier 1
metrics and ROUGE identical (max difference 0.0), BERTScore max difference
1.8 × 10⁻⁷, Tribunal trigger identical on all 159. The CPU rows are kept in
`scores/e0_dry_cpu_reference/`.

---

# Roster robustness cohort (2026-09-18)

Two open-weight models from the proposed main roster, run over the same 60 frozen
items to ask one question: if the main benchmark moves to smaller open-weight
models, does an evaluator validated on the pilot's five still work there? They are
NOT in the labelled 300 - the question is mechanical, so it needs no human labels,
and 420 traces would have been 40% more expert annotation.

## R-F1 · Structure does not degrade. The formatting worry was unfounded

| Model | Traces | Avg steps | No step marker | Has `**Answer:**` | Final unparseable | Numbers per step |
|---|---|---|---|---|---|---|
| gemma-4-31b-it | 60 | 6.8 | 0% | **100%** | 0% | **0.99** |
| qwen3.8-27b | 46 clean | 7.3 | 0% | **100%** | 0% | 0.96 |
| claude-opus-4.7 | 60 | 6.8 | 0% | 100% | 0% | 0.89 |
| gemini-3.1-pro | 60 | 7.2 | 0% | 100% | 0% | 0.96 |
| deepseek-r1 | 60 | 8.5 | 0% | 95% | 0% | 0.98 |
| gpt-5 | 60 | 6.0 | 0% | 87% | 0% | 0.94 |
| llama-3.1-70b | 60 | 6.3 | 0% | 70% | 0% | 0.94 |

Not one trace in 416, from any model, lacks step markers, and every trace yields a
parseable final answer. The two open-weight models hold the `**Answer:**` line in
100% of traces - better than GPT-5 at 87% and Llama at 70% - and Gemma carries a
number on 99% of its steps, the highest of any model here, which is the property
E3's milestone extraction depends on.

So the risk that motivated this cohort is not real: an evaluator that reads step
structure will function on a small open-weight roster. The prompt does the work,
and a 31B model follows it as well as a frontier one.

## R-F2 · The real small-model failure is runaway reasoning, not formatting

| Model | Usable | Unusable | Cost, 60 traces |
|---|---|---|---|
| gemma-4-31b-it (no reasoning) | **60 of 60** | 0 | **$0.02** |
| qwen3.8-27b (reasoning) | 46 clean, 10 truncated | 14 of 60 | $2.13 |

Qwen3.8-27B never answered 4 items even with a 32,768-token budget: on
`normal_depth_iteration` it wrote 50,613 characters of reasoning and stopped. Ten
more answers end mid-derivation. The same failure DeepSeek R1, GPT-5 and Gemini all
showed, and the reason three ceilings were raised during this pilot.

Two consequences for roster planning. A reasoning model's cost is not predictable
from its price per token: Qwen is 100x cheaper than Gemma per token and cost 100x
more here, because it thinks. And a small reasoning model can simply fail to answer
the hardest items, which is a coverage hole in a benchmark table, not a low score -
worth reporting as such rather than averaging away.

Gemma-4-31B-it is the strongest anchor this cohort found: 60 of 60, the cleanest
structure measured, $0.02 for the set.
