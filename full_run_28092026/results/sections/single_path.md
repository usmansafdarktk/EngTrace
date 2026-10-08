# Consistency within a template: single-path and other templates

**QUICK**: development values, 1,000 bootstrap draws and 10,000 sign flips or permutations (the full pass draws 10,000 and 100,000); the integration pass replaces every number here.

Printed by `analyze.py` (WS-C1) from Q4, default configuration (`scores/main`). Per model, the share of the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) fully solved on all 15 instances, on some and on none, and the same for the 92 others; 95% intervals, templates resampled.

| model | single: all | 95% CI | some | 95% CI | none | 95% CI | others: all | 95% CI | some | 95% CI | none | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `glm-5.3-flash` | 0.914 | 0.828 to 0.983 | 0.086 | 0.017 to 0.155 | 0.000 | 0.000 to 0.000 | 0.826 | 0.750 to 0.902 | 0.174 | 0.098 to 0.261 | 0.000 | 0.000 to 0.000 |
| `deepseek-v4.1-flash` | 0.897 | 0.810 to 0.966 | 0.103 | 0.034 to 0.173 | 0.000 | 0.000 to 0.000 | 0.848 | 0.772 to 0.913 | 0.152 | 0.076 to 0.239 | 0.000 | 0.000 to 0.000 |
| `claude-sonnet-5` | 0.897 | 0.810 to 0.966 | 0.103 | 0.034 to 0.172 | 0.000 | 0.000 to 0.000 | 0.826 | 0.761 to 0.902 | 0.174 | 0.109 to 0.250 | 0.000 | 0.000 to 0.000 |
| `muse-glimmer-30b` | 0.862 | 0.775 to 0.931 | 0.138 | 0.052 to 0.241 | 0.000 | 0.000 to 0.000 | 0.826 | 0.750 to 0.902 | 0.174 | 0.109 to 0.250 | 0.000 | 0.000 to 0.000 |
| `kimi-k3` | 0.862 | 0.775 to 0.948 | 0.138 | 0.052 to 0.224 | 0.000 | 0.000 to 0.000 | 0.793 | 0.707 to 0.870 | 0.207 | 0.130 to 0.294 | 0.000 | 0.000 to 0.000 |
| `glm-5.3` | 0.845 | 0.741 to 0.931 | 0.138 | 0.052 to 0.224 | 0.017 | 0.000 to 0.052 | 0.772 | 0.685 to 0.848 | 0.217 | 0.141 to 0.304 | 0.011 | 0.000 to 0.033 |
| `gemini-3.1-flash-lite` | 0.672 | 0.552 to 0.793 | 0.276 | 0.155 to 0.397 | 0.052 | 0.000 to 0.121 | 0.598 | 0.500 to 0.696 | 0.402 | 0.304 to 0.500 | 0.000 | 0.000 to 0.000 |
| `gpt-5.4-mini` | 0.483 | 0.345 to 0.621 | 0.483 | 0.362 to 0.621 | 0.034 | 0.000 to 0.086 | 0.467 | 0.359 to 0.576 | 0.511 | 0.402 to 0.609 | 0.022 | 0.000 to 0.054 |
| `gemma-4-26b-a4b` | 0.466 | 0.328 to 0.587 | 0.500 | 0.379 to 0.638 | 0.034 | 0.000 to 0.086 | 0.511 | 0.413 to 0.609 | 0.489 | 0.380 to 0.598 | 0.000 | 0.000 to 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.328 to 0.569 | 0.534 | 0.414 to 0.672 | 0.017 | 0.000 to 0.052 | 0.543 | 0.435 to 0.641 | 0.435 | 0.337 to 0.533 | 0.022 | 0.000 to 0.054 |
| `gpt-oss-20b` | 0.362 | 0.241 to 0.466 | 0.569 | 0.448 to 0.707 | 0.069 | 0.017 to 0.138 | 0.337 | 0.239 to 0.435 | 0.652 | 0.565 to 0.739 | 0.011 | 0.000 to 0.033 |

The "none" column (templates no instance of which the model fully solves) ranges from 0.000 to 0.069 over the models and both groups.
