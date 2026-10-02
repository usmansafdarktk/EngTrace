# analysis/ — the scripts behind every reported number

Each finding and results table in this pilot was computed by a script here, and every
script here is listed below. Run them from the repository root:

```bash
python evaluator_pilot_17092026/analysis/<script>.py
```

Most need only the system Python and make no network call. The exceptions are marked in
the table:

- **paid**: the mode named bills API tokens and must not run without the owner's approval
  for that run. The script's other modes are free.
- **net**: reads OpenRouter's public catalogue or live prices. No tokens, no cost.
- **venv**: imports the published framework or builds its judge prompts, so it needs the
  pinned scorer stack in `.venv` (D-081).

They read `slice/`, `traces/`, `scores/`, `experts_filled_labels/` and `analysis/out/`. All
but `slice/` are gitignored, and the expert labels are kept out of the repository by
decision and never committed. On a fresh clone the runs have to be regenerated first
(`run_traces.py`, then `run_evaluator.py`) and the label files supplied by the authors;
regenerated traces will not reproduce these numbers exactly, because model outputs vary
between calls.

| Script | Reproduces | Where it is reported |
|---|---|---|
| `milestone_probe.py` | all 60 frozen items regenerate **byte-identically** under frame-local capture | RESULTS_E3 §How E3 decides |
| `gold_parse.py` | how E0's parser reads each template's gold answer; `rackett` reads the temperature; 8 of 15 put no number on the answer line | FINDINGS E0-F1, E0-F2 |
| `trace_shape.py` | step markers, `**Answer:**` presence, parseability across all seven models | FINDINGS R-F1 |
| `e0_retest.py` | E0 run-to-run change split into Tier 1 / judges / sampling | FINDINGS E0-F7 |
| `compare_e0_e3.py` | E3 vs E0 per model, correlations, largest disagreements — **current configs only** | RESULTS_E3 §Results, §Where they disagree |
| `e3_analysis.py` | E3 per model and answer type, the `rackett` check, raw null | RESULTS_E3 |
| `e3_null.py` | the null split into shared values (3.8%) and coincidence | RESULTS_E3 §The null baseline |
| `e3_grid.py` | E3 real / null / separation across tolerance x unit scaling | RESULTS_E3 table, `e3_milestones.py` docstring |
| `e4_analysis.py` | E4 vs E3 per model, milestone statuses, parse coverage, every contradicted milestone printed for reading | RESULTS_E4 §Results, §The finding |
| `e4_claim_audit.py [N]` | a fixed-seed sample of E4's inconsistent claims, for hand classification | RESULTS_E4 §The side metric (19 real / 10 FP / 1 unclear of 30) |
| `arith_gold_validation.py [--rule]` | the arithmetic checker's gold validation under the 1% reading and the digit rule, with the two fixes it forced | RESULTS_E4 §The gold validation, D-101 |
| `digit_rule.py` | the experts' digit rule by machine, four rules side by side, per step and per trace inside correct-answer traces | RESULTS_X1 Finding 5, RESULTS_E4, D-097, D-101 |
| `judge_candidates.py` | every OpenRouter family outside the 27-model suite, with price, openness and JSON support (**net**) | JUDGE_SELECTION §The constraint |
| `judge_probe.py build/run/report [--offline]` | 21 known-label steps sent to 9 judges with the framework's own prompt, two labels corrected in `RELABEL` (`run` is **paid**, ~$3; `report` reads live prices unless `--offline`; **venv**) | JUDGE_SELECTION §The probe; FINDINGS E0-F5 update; D-112 |
| `e1_panel_check.py` | the E1 panel's availability on OpenRouter (listing, endpoints, a live JSON call that bills a few tokens) and E1's cost from probe token usage (**net**) | JUDGE_SELECTION; E1 cost estimate |
| `e1_analysis.py` | E1 run health, E1 vs E0 and E0-3J, inter-judge agreement (Cohen's and Fleiss' kappa) for both panels | RESULTS_E1; FINDINGS E0-F8 |
| `e5_analysis.py validate / report` | the E5 judge's validity on known answers (`validate` is **paid**, ~$0.43), then E5 against E3, E4 and E0 | RESULTS_E5 |
| `e2_analysis.py` | the three PRMs on the probe's known-label steps, per model, their agreement, and against E0/E3/E5 | RESULTS_E2 |
| `prm_threshold.py` | E2's 0.5 threshold fitted on one half of the labelled traces and reported on the other | RESULTS_E2 §The 0.5 threshold, D-100 |
| `x1_analysis.py` | every evaluator against the expert labels: trace-level AUROC with bootstrap intervals and baselines, the correct-answer hard case, step level, per model, without chemical engineering | RESULTS_X1 |
| `answer_check.py` | the corrected final-answer check against the expert verdicts: agreement, E0's errors by direction, the split-half tolerance, per-model accuracy, every residual disagreement | RESULTS_X1 Finding 1b, E0_RERUN, D-098 |
| `hard_case_pool.py` | the hard case on the step labels (93 of 228 correct-answer traces flawed), and ranked candidates for a second labelling round | RESULTS_X1 Finding 5, D-096 |
| `cluster_bootstrap.py` | design effect, template-level intervals and the smallest detectable difference; section 1b, the experts' answer verdict minus every evaluator | RESULTS_X1 Finding 6, D-111 |
| `planted.py` | builds the planted-defect set (seeded, byte-identical) and scores the deterministic evaluators on it | RESULTS_X1 Finding 7, D-102 |
| `planted_judges.py build/run/report` | the matched judge probe on the planted set: GPT-5, Opus 4.5, MiMo, each defect judged planted and untouched (`run` is **paid**; `report` is free; **venv**) | RESULTS_X1 Finding 8, D-103 |
| `planted_routing.py` | E0's Tier 1 routing on the planted defects, the Tribunal replaced by a recorder, $0 (**venv**) | RESULTS_X1 Finding 8, D-104 |
| `router_residue.py [--models]` | the steps no checker can verify, what the experts said about them, a router's full-run cost batched and per step, and where every step the experts call incorrect lands (flagged, shown to a judge, or neither) | RESULTS_X1 Finding 8, PILOT_SUMMARY §3.1, D-110, D-112 |
| `router_planted.py` | a verify-first router on the planted defects: forwarded, caught when asked, and both jointly; E0's routing on the same basis; the digit rule's false alarms on the judged steps and the slips it gives up (**venv**) | PILOT_SUMMARY §2.7-2.9, RESULTS_X1 Findings 7-8, D-111 |
| `router_batched.py build/run/report` | the batched router's smoke checks: rule C, a trace's steps in one MiMo call, on the planted defects and on the 300 labelled traces; `build` is free and prints the estimate, `run CHECK --budget` is **paid** and capped (**venv**) | RESULTS_X1 "The router, batched", PILOT_SUMMARY §2.9 and §3.1, D-113 |
| `lojo.py` | leave-one-judge-out on E0-3J, replayed offline from the cached Tier 1 matrices and the stored replies (reproduces all 179 judged F1s), with placebo and control drops; the panel against the experts' step labels; X2, each judge's bias per trace model with the family effect (**venv**) | RESULTS_LOJO; D-174 |
| `attribution.py` | each component's flags (the digit rule, the router's judge by category, E5 MISSING) against the experts' step error types: precision and recall by type, steps and traces (**venv**) | RESULTS_ATTRIBUTION; D-176 |
| `summary_numbers.py` | every value the pilot summary's figures draw, and the prose numbers no other script prints, written to `figures/summary_numbers.json` (**venv**) | PILOT_SUMMARY; `make_summary_figures.py` |
| `inference_cost.py` | per-model inference cost of the full benchmark from the pilot's measured token counts (**net**: live prices) | the inference pricing document |
| `judge_cost.py [--open N --closed M]` | the judge cost of E0, E0-3J, E1 and E5 on the full run, per trace and at the roster mix | D-105, D-110, PILOT_SUMMARY §3 |
| `roster_cost.py` | per-model cost of the full 2,250-item benchmark under E0 | conversation with supervisor; roster planning |
| `roster_candidates.py` | open-weight roster candidates from OpenRouter's catalogue, with the families a judge or a PRM blocks (**net**) | roster planning |
| `model_health.py` | every pilot model listed, un-deprecated, served, answering (makes live calls, ~$0.01: **paid**) | README §The five models |

## Two null baselines, deliberately different

`e3_analysis.py` prints a **raw** null — a trace scored against a sibling item's
milestones as they stand: **0.072** at the 0.5% tolerance.

`e3_grid.py` and `e3_null.py` first remove the sibling milestones that equal one of
this item's own values or givens, because reaching a genuinely shared value is not
coincidence. That **coincidence-only** null is **0.040**, and it is the one
RESULTS_E3 reports, because it is the rate at which E3 credits a trace for something
it did not do.

Both are correct. The gap between them (~0.03) is the shared-value effect.

## Rule for anything added here

A result goes into FINDINGS, RESULTS or the summary only with the script that computed it,
the script is added to the table above, and a script reads **the evaluator's current
config only**. One of these scripts used to keep the last row per trace regardless of
config; after the E0 re-seed it silently mixed two runs into one table. That is fixed, and
the filter is the rule.

A number is quoted with its denominator and its basis, and a figure draws from computed
values rather than typed ones (D-111: the summary's router bars were typed in, and one of
them had never been measured).
