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
Google's named successor and records the substitution on every scored row. Listing
the model through the Models API still succeeds, which is why a listing check is
not an availability check.
