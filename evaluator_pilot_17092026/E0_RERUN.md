# E0 re-run: fix the final-answer check, then score E0 once more

One task, one run. Everything here is about E0's answer check and nothing else.

## Why

E0 decides whether a trace's final answer is right, and it disagrees with the experts on
**72 of 300 traces** — 68 of them traces the experts call **correct** and E0 calls wrong.
In 45 of those the gold value is sitting in the trace's own answer segment
(`analysis/answer_check.py`). Two known defects cause it: the gold's answer line restates
the inputs (E0-F1), and 8 of 15 templates answer with a word, not a number (E0-F2).

What that does to the headline number:

| | experts | E0 reports |
|---|---|---|
| gpt-5 | 0.950 | **0.617** |
| all 300 traces | 0.760 | **0.547** |

The benchmark understates every model by about 21 points, and the ranking is wrong: E0
puts claude-opus-4.7 first and **GPT-5 fourth**, where the experts put GPT-5 first. A
results table regenerated with the current check reproduces that error. The check also
gates the judges (wrong-answer traces are sampled to the tribunal at 0.20), so the 68
mislabelled traces move the reasoning score and the bill too.

## Scope

**In:** E0's answer check, and one re-score of E0 over the existing 300 traces.

**Out, deliberately:**
- No re-annotation. The expert labels are the yardstick, not the thing being fixed.
- No new traces, no new models.
- No E0-3J or E1 re-run. Their comparison with E0 ("swapping judges barely moves the
  score") is controlled — it needs both arms on the *same* check, not a correct one. It
  stays as measured, reported as "under the published framework's answer check". Revisit
  only if a reviewer asks, $12.20.
- Not the digit rule (D-092), not the PRM threshold, not parse coverage. Separate tasks.

## How

The check must match the task, not patch the parser:

1. Read the trace's **answer segment** — after its final Answer marker.
2. **Split it into parts.** These items ask for (a) a value and (b) a regime.
3. Score each part; return **correct / partial / incorrect**. 19 of 300 traces are
   genuinely partial, and a binary check cannot represent them.
4. **Accept a correct rounding** (56.27 for 56.29). The opposite of the digit rule, on
   purpose: there a wrong digit is an error, here a rounded display is not.
5. Take the gold value from the gold's **last computed number**, never its answer line
   (E0-F1).
6. Require the **word** when the gold's answer is qualitative (E0-F2).

Every change is measured against the 300 expert verdicts before it is kept.

## Action items

| # | item | cost | done when |
|---|---|---|---|
| 1 | Build the per-part, three-way check | $0 | `analysis/answer_check.py` reports it |
| 2 | Fix the six templates that carry the residual errors: `aoq_ati_rectifying` (15), `rackett_equation_volume` (11), `manning_rectangular_discharge` (11), `reynolds_number_flow_regime` (9), `fluid_particle_acceleration` (8), `coaxial_capacitance` (6) | $0 | agreement ≥ 0.90 on the 281 non-partial traces (now 0.797, E0 0.747) |
| 3 | Cache `CROSS_ENCODER.predict` per pair and add `--workers`, verified output-identical | $0 | a 5-trace re-score matches the current rows exactly |
| 4 | **Re-score E0 over the 300 traces with the corrected check** | **$6.53** | new rows in `scores/e0/`, new `config_sha256` |
| 5 | Regenerate the affected tables and documents | $0 | RESULTS_E0, RESULTS_X1 Finding 1, FINDINGS E0-F1/F2 |
| 6 | Record the decision | $0 | a D-09x entry in DECISIONS.md |

Item 4 is the only paid step and needs explicit approval before it runs. Items 1–3 are
independent and can run in parallel; item 4 waits on 1–2 being validated, because
re-scoring with a half-finished check wastes the run.

**Total: $6.53 and roughly 3–4 hours**, most of it unattended.

## Checks that keep this honest

- Agreement against the expert verdicts is reported before and after, on the same 300
  traces, with the per-model table.
- The corrected check's residual disagreements are listed by template, not summarised.
- The old E0 rows are kept; the re-run adds rows under a new config hash rather than
  replacing them, so both are reportable.

## Status, 2026-09-24

**Items 1-3 are done** (`cb84f32`, `2a7900d`, and the harness commit above).

- The check is `evaluators/answer.py`, measured by `analysis/answer_check.py`:
  **0.947 against the experts on the 281 non-partial traces, where E0 scores 0.747**;
  0.893 three-way over all 300. The tolerance is chosen on one half of the traces and
  reported on the other (0.876 / 0.905). A 10-case self-test pins every defect fixed.
- 32 disagreements remain, 15 of them `aoq_ati_rectifying` (its question asks six
  quantities and the gold states one on its answer line) and 7 `lorentz_force`.
- The harness is faster and its cache is now correct: the old per-pair cache never fired,
  because Tier 1's matrix cache sits above it, and three defects meant it was not
  output-identical when it did. Judge calls are now fetched concurrently while scoring
  stays serial, because the wrong-answer sample reads a process-wide seed. Verified
  offline: 300/300 D3 draws identical, 21/21 rows identical at 1 and 16 workers, 14/14
  traces identical with the pair cache on and off. **3.9 h becomes 12-50 min.**

**Item 4 has not run, and its cost was understated.** Two things changed:

1. `e0_tribunal.py` does not import the new check yet. Wiring it changes E0's config
   hash, which is what makes the re-run score anything at all.
2. The corrected check calls 75 more traces correct, and 288 of 300 traces have a Tier-1
   match ratio below 0.80, so nearly every newly-correct trace now reaches the judges:
   roughly **233 judged traces against today's 178, about +30%**, so the re-run is
   **~$8.50, not $6.53**.

**The open decision.** The corrected accuracy needs no paid run - `analysis/answer_check.py`
computes it offline from the existing traces for $0. The $8.50 buys E0's *reasoning* score
under the corrected gating, which differs on the traces whose judging changes. Report the
corrected accuracy offline, or re-run E0 and have both under one configuration.
