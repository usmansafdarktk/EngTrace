# analysis/ — the scripts behind every reported number

Each finding and results table in this pilot was computed by a script here. Run any
of them from the repository root with the system Python; none makes an API call
except `model_health.py`, and none needs the `.venv`.

```bash
python evaluator_pilot_17092026/analysis/<script>.py
```

They read `slice/`, `traces/` and `scores/`. `traces/` and `scores/` are gitignored,
so on a fresh clone these scripts need the runs regenerated first
(`run_traces.py`, then `run_evaluator.py`).

| Script | Reproduces | Where it is reported |
|---|---|---|
| `milestone_probe.py` | all 60 frozen items regenerate **byte-identically** under frame-local capture | RESULTS_E3 §How E3 decides |
| `gold_parse.py` | how E0's parser reads each template's gold answer; `rackett` reads the temperature | FINDINGS E0-F1, E0-F2 |
| `trace_shape.py` | step markers, `**Answer:**` presence, parseability across all seven models | FINDINGS R-F1 |
| `e0_retest.py` | E0 run-to-run change split into Tier 1 / judges / sampling | FINDINGS E0-F7 |
| `compare_e0_e3.py` | E3 vs E0 per model, correlations, largest disagreements — **current configs only** | RESULTS_E3 §Results, §Where they disagree |
| `e3_analysis.py` | E3 per model and answer type, the `rackett` check, raw null | RESULTS_E3 |
| `e3_null.py` | the null split into shared values (3.8%) and coincidence | RESULTS_E3 §The null baseline |
| `e3_grid.py` | E3 real / null / separation across tolerance x unit scaling | RESULTS_E3 table, `e3_milestones.py` docstring |
| `roster_cost.py` | per-model cost of the full 2,250-item benchmark under E0 | conversation with supervisor; roster planning |
| `model_health.py` | every pilot model listed, un-deprecated, served, answering (makes live calls, ~$0.01) | README §The five models |

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

A result goes into FINDINGS or RESULTS only with the script that computed it, and a
script reads **the evaluator's current config only**. One of these scripts used to
keep the last row per trace regardless of config; after the E0 re-seed it silently
mixed two runs into one table. That is fixed, and the filter is the rule.
