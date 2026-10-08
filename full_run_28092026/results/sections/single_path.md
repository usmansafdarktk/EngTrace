# Consistency within a template: single-path and other templates

Printed by `analyze.py` (WS-C1) from Q4, default configuration (`scores/main`). Per model, the share of the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) fully solved on all 15 instances, on some and on none, and the same for the 92 others; 95% intervals, templates resampled.

| model | single: all | 95% CI | some | 95% CI | none | 95% CI | others: all | 95% CI | some | 95% CI | none | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.931 | 0.862 to 0.983 | 0.069 | 0.017 to 0.138 | 0.000 | 0.000 to 0.000 | 0.880 | 0.815 to 0.946 | 0.120 | 0.054 to 0.185 | 0.000 | 0.000 to 0.000 |
| `claude-sonnet-5` | 0.914 | 0.844 to 0.983 | 0.086 | 0.017 to 0.155 | 0.000 | 0.000 to 0.000 | 0.837 | 0.761 to 0.913 | 0.163 | 0.087 to 0.239 | 0.000 | 0.000 to 0.000 |
| `glm-5.3-flash` | 0.897 | 0.810 to 0.966 | 0.103 | 0.034 to 0.190 | 0.000 | 0.000 to 0.000 | 0.826 | 0.739 to 0.902 | 0.174 | 0.098 to 0.250 | 0.000 | 0.000 to 0.000 |
| `kimi-k3` | 0.862 | 0.776 to 0.948 | 0.138 | 0.052 to 0.225 | 0.000 | 0.000 to 0.000 | 0.815 | 0.728 to 0.891 | 0.185 | 0.109 to 0.261 | 0.000 | 0.000 to 0.000 |
| `glm-5.3` | 0.845 | 0.741 to 0.931 | 0.138 | 0.052 to 0.224 | 0.017 | 0.000 to 0.052 | 0.793 | 0.707 to 0.880 | 0.196 | 0.120 to 0.283 | 0.011 | 0.000 to 0.033 |
| `muse-glimmer-30b` | 0.776 | 0.672 to 0.879 | 0.224 | 0.121 to 0.328 | 0.000 | 0.000 to 0.000 | 0.826 | 0.750 to 0.902 | 0.174 | 0.098 to 0.261 | 0.000 | 0.000 to 0.000 |
| `gemini-3.1-flash-lite` | 0.672 | 0.552 to 0.793 | 0.276 | 0.172 to 0.397 | 0.052 | 0.000 to 0.121 | 0.598 | 0.500 to 0.696 | 0.402 | 0.304 to 0.500 | 0.000 | 0.000 to 0.000 |
| `gpt-5.4-mini` | 0.500 | 0.379 to 0.638 | 0.431 | 0.310 to 0.552 | 0.069 | 0.017 to 0.138 | 0.478 | 0.380 to 0.587 | 0.500 | 0.391 to 0.598 | 0.022 | 0.000 to 0.054 |
| `gemma-4-26b-a4b` | 0.466 | 0.345 to 0.603 | 0.500 | 0.362 to 0.621 | 0.034 | 0.000 to 0.086 | 0.522 | 0.424 to 0.620 | 0.478 | 0.380 to 0.576 | 0.000 | 0.000 to 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.328 to 0.586 | 0.534 | 0.414 to 0.655 | 0.017 | 0.000 to 0.052 | 0.554 | 0.457 to 0.652 | 0.424 | 0.326 to 0.522 | 0.022 | 0.000 to 0.054 |
| `gpt-oss-20b` | 0.345 | 0.224 to 0.466 | 0.586 | 0.466 to 0.707 | 0.069 | 0.017 to 0.138 | 0.293 | 0.207 to 0.391 | 0.696 | 0.598 to 0.783 | 0.011 | 0.000 to 0.033 |

The "none" column (templates no instance of which the model fully solves) ranges from 0.000 to 0.069 over the models and both groups.
