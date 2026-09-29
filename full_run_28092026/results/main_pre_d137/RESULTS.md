# Results of the full run: the deterministic stack

**A record, not the results:** computed from `scores/main_pre_d137`. The same traces scored by the evaluators as they were before D-137, D-138 and D-139 (score.py at a72399b); kept to show what those fixes changed.

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring. The answer check, E3 and E4's digit rule only: E5's judge has not run, so its columns are empty. Every interval is 95% and resamples the 150 templates (B = 10,000).

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0.

| model | answer score | 95% CI | fully solved | 95% CI | correct | partial | incorrect | unusable |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `claude-sonnet-5` | 0.967 | 0.947 to 0.983 | 0.961 | 0.940 to 0.979 | 2163 | 26 | 61 | 0 |
| `muse-glimmer-30b` | 0.959 | 0.936 to 0.978 | 0.954 | 0.929 to 0.975 | 2146 | 24 | 51 | 29 |
| `glm-5.3-flash` | 0.954 | 0.931 to 0.974 | 0.950 | 0.926 to 0.971 | 2138 | 18 | 52 | 42 |
| `kimi-k3` | 0.942 | 0.919 to 0.963 | 0.933 | 0.906 to 0.956 | 2099 | 43 | 103 | 5 |
| `glm-5.3` | 0.942 | 0.912 to 0.968 | 0.937 | 0.906 to 0.964 | 2109 | 22 | 19 | 100 |
| `deepseek-v4.1-flash` | 0.941 | 0.912 to 0.965 | 0.929 | 0.896 to 0.957 | 2091 | 51 | 94 | 14 |
| `qwen3-235b-a22b-2507` | 0.873 | 0.836 to 0.906 | 0.851 | 0.809 to 0.890 | 1914 | 99 | 236 | 1 |
| `gemini-3.1-flash-lite` | 0.855 | 0.808 to 0.897 | 0.845 | 0.798 to 0.888 | 1901 | 44 | 305 | 0 |
| `gemma-4-26b-a4b` | 0.847 | 0.804 to 0.886 | 0.831 | 0.787 to 0.872 | 1869 | 75 | 306 | 0 |
| `gpt-5.4-mini` | 0.811 | 0.766 to 0.853 | 0.792 | 0.745 to 0.836 | 1783 | 83 | 384 | 0 |
| `gpt-oss-20b` | 0.807 | 0.764 to 0.847 | 0.789 | 0.745 to 0.832 | 1776 | 79 | 340 | 55 |

Of the 55 pairs, 32 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The fully-solved check (McNemar exact, Holm) gives the same verdict on 43 of 55: both tests or neither hold, and when both hold, in the same direction. On 12 of the 12 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level test.

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | McNemar p (Holm) | same verdict |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.160 | -0.202 to -0.119 | 0.0055 | -0.172 | -0.215 to -0.130 | 0.0000 | yes |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.156 | -0.198 to -0.117 | 0.0055 | -0.169 | -0.214 to -0.127 | 0.0000 | yes |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.152 | -0.190 to -0.117 | 0.0055 | -0.164 | -0.204 to -0.128 | 0.0000 | yes |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.148 | 0.113 to 0.186 | 0.0055 | +0.161 | 0.125 to 0.200 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.147 | -0.189 to -0.108 | 0.0055 | -0.161 | -0.204 to -0.119 | 0.0000 | yes |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.143 | 0.109 to 0.181 | 0.0055 | +0.158 | 0.121 to 0.197 | 0.0000 | yes |
| `gpt-oss-20b` | `kimi-k3` | -0.136 | -0.179 to -0.094 | 0.0055 | -0.144 | -0.189 to -0.100 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3` | -0.135 | -0.181 to -0.090 | 0.0055 | -0.148 | -0.194 to -0.100 | 0.0000 | yes |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.134 | -0.178 to -0.090 | 0.0055 | -0.140 | -0.185 to -0.094 | 0.0000 | yes |
| `kimi-k3` | `gpt-5.4-mini` | +0.132 | 0.096 to 0.168 | 0.0055 | +0.140 | 0.104 to 0.180 | 0.0000 | yes |
| `glm-5.3` | `gpt-5.4-mini` | +0.131 | 0.094 to 0.170 | 0.0055 | +0.145 | 0.105 to 0.187 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.130 | 0.093 to 0.168 | 0.0055 | +0.137 | 0.098 to 0.176 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.120 | -0.156 to -0.086 | 0.0055 | -0.131 | -0.169 to -0.095 | 0.0000 | yes |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.112 | -0.152 to -0.077 | 0.0055 | -0.116 | -0.156 to -0.080 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.112 | -0.148 to -0.079 | 0.0055 | -0.123 | -0.163 to -0.086 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.107 | -0.141 to -0.074 | 0.0055 | -0.120 | -0.158 to -0.084 | 0.0000 | yes |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.104 | 0.068 to 0.143 | 0.0055 | +0.109 | 0.071 to 0.149 | 0.0000 | yes |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.100 | 0.067 to 0.135 | 0.0055 | +0.105 | 0.070 to 0.144 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.095 | -0.133 to -0.060 | 0.0055 | -0.102 | -0.141 to -0.064 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.095 | -0.131 to -0.061 | 0.0055 | -0.107 | -0.145 to -0.072 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.094 | -0.125 to -0.066 | 0.0055 | -0.111 | -0.148 to -0.076 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.093 | -0.133 to -0.056 | 0.0055 | -0.099 | -0.140 to -0.060 | 0.0000 | yes |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.088 | 0.054 to 0.125 | 0.0055 | +0.088 | 0.052 to 0.128 | 0.0000 | yes |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.088 | 0.055 to 0.122 | 0.0055 | +0.092 | 0.060 to 0.129 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.086 | -0.117 to -0.058 | 0.0055 | -0.103 | -0.141 to -0.069 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.086 | 0.047 to 0.127 | 0.0055 | +0.084 | 0.044 to 0.127 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.082 | -0.111 to -0.054 | 0.0055 | -0.100 | -0.135 to -0.067 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.070 | -0.100 to -0.042 | 0.0055 | -0.082 | -0.120 to -0.048 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.070 | -0.101 to -0.040 | 0.0055 | -0.087 | -0.125 to -0.052 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.068 | 0.038 to 0.100 | 0.0055 | +0.079 | 0.040 to 0.118 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.062 | 0.031 to 0.094 | 0.0055 | +0.058 | 0.018 to 0.097 | 0.0000 | yes |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.044 | -0.071 to -0.017 | 0.0360 | -0.052 | -0.084 to -0.022 | 0.0000 | yes |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.066 | -0.111 to -0.023 | 0.0828 | -0.061 | -0.110 to -0.012 | 0.0000 | no |
| `glm-5.3` | `claude-sonnet-5` | -0.025 | -0.049 to -0.006 | 0.4114 | -0.024 | -0.048 to -0.005 | 0.0000 | no |
| `kimi-k3` | `claude-sonnet-5` | -0.025 | -0.046 to -0.005 | 0.4683 | -0.028 | -0.054 to -0.006 | 0.0000 | no |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.036 | 0.005 to 0.068 | 0.5139 | +0.038 | 0.005 to 0.071 | 0.0003 | no |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.026 | -0.052 to -0.003 | 0.6725 | -0.032 | -0.063 to -0.004 | 0.0000 | no |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.048 | -0.097 to 0.002 | 1.0000 | -0.056 | -0.106 to -0.004 | 0.0000 | no |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.040 | -0.084 to 0.003 | 1.0000 | -0.041 | -0.087 to 0.003 | 0.0006 | no |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.025 | -0.055 to 0.004 | 1.0000 | -0.020 | -0.057 to 0.017 | 0.2556 | yes |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | -0.018 | -0.045 to 0.006 | 1.0000 | -0.024 | -0.053 to 0.002 | 0.0006 | no |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.018 | -0.010 to 0.048 | 1.0000 | +0.006 | -0.032 to 0.042 | 1.0000 | yes |
| `glm-5.3` | `muse-glimmer-30b` | -0.017 | -0.044 to 0.008 | 1.0000 | -0.016 | -0.046 to 0.010 | 0.0499 | no |
| `muse-glimmer-30b` | `kimi-k3` | +0.017 | -0.007 to 0.040 | 1.0000 | +0.021 | -0.003 to 0.046 | 0.0066 | no |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.014 | -0.033 to 0.005 | 1.0000 | -0.021 | -0.045 to 0.001 | 0.0024 | no |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.013 | -0.032 to 0.004 | 1.0000 | -0.011 | -0.032 to 0.008 | 0.3462 | yes |
| `glm-5.3-flash` | `glm-5.3` | +0.012 | -0.007 to 0.031 | 1.0000 | +0.013 | -0.007 to 0.033 | 0.1553 | yes |
| `glm-5.3-flash` | `kimi-k3` | +0.012 | -0.007 to 0.031 | 1.0000 | +0.017 | -0.004 to 0.039 | 0.0364 | no |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.008 | -0.031 to 0.015 | 1.0000 | -0.008 | -0.033 to 0.018 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.007 | -0.035 to 0.020 | 1.0000 | -0.014 | -0.043 to 0.015 | 0.5711 | yes |
| `glm-5.3-flash` | `muse-glimmer-30b` | -0.005 | -0.023 to 0.013 | 1.0000 | -0.004 | -0.022 to 0.015 | 1.0000 | yes |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.004 | -0.049 to 0.041 | 1.0000 | -0.003 | -0.049 to 0.042 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.002 | -0.016 to 0.014 | 1.0000 | -0.004 | -0.020 to 0.015 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.002 | -0.032 to 0.029 | 1.0000 | -0.008 | -0.042 to 0.025 | 1.0000 | yes |
| `glm-5.3` | `kimi-k3` | -0.000 | -0.030 to 0.028 | 1.0000 | +0.004 | -0.028 to 0.035 | 1.0000 | yes |

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier; p from permuting the tier labels, Holm across the eleven models. Sigma is the within-tier SD of the template-mean score; the detectable gap, 2.80 x sigma x sqrt(1/58 + 1/34), is the smallest this test finds at 80% power.

| model | Easy | Advanced | gap | 95% CI | p (Holm) | sigma | detectable gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.859 | 0.760 | +0.099 | -0.012 to 0.217 | 0.5795 | 0.262 | 0.158 |
| `gemma-4-26b-a4b` | 0.897 | 0.771 | +0.127 | 0.025 to 0.233 | 0.1359 | 0.238 | 0.144 |
| `deepseek-v4.1-flash` | 0.945 | 0.932 | +0.013 | -0.053 to 0.090 | 1.0000 | 0.160 | 0.097 |
| `qwen3-235b-a22b-2507` | 0.900 | 0.839 | +0.061 | -0.033 to 0.157 | 0.9929 | 0.212 | 0.128 |
| `glm-5.3-flash` | 0.966 | 0.932 | +0.033 | -0.018 to 0.090 | 0.9929 | 0.123 | 0.074 |
| `glm-5.3` | 0.975 | 0.851 | +0.124 | 0.033 to 0.230 | 0.0297 | 0.190 | 0.115 |
| `muse-glimmer-30b` | 0.954 | 0.959 | -0.005 | -0.058 to 0.048 | 1.0000 | 0.136 | 0.082 |
| `kimi-k3` | 0.955 | 0.929 | +0.026 | -0.031 to 0.086 | 1.0000 | 0.137 | 0.083 |
| `gpt-5.4-mini` | 0.857 | 0.740 | +0.117 | 0.001 to 0.239 | 0.3392 | 0.264 | 0.160 |
| `gemini-3.1-flash-lite` | 0.911 | 0.734 | +0.177 | 0.058 to 0.306 | 0.0360 | 0.266 | 0.161 |
| `claude-sonnet-5` | 0.979 | 0.928 | +0.051 | -0.009 to 0.122 | 0.5795 | 0.136 | 0.082 |

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 70 items with no milestones (all 15 of 3 templates and some items of 12 more); 52 templates have an item with at most one milestone (46 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The digit rule's flag counts a trace when any step is flagged. Against the experts, on fully solved traces, it had precision 0.750 and recall 0.320 (SCORER_VALIDATION.md), so its rate is not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 395 | 55 | 0.455 | 0.388 to 0.529 | 0.529 | 0.465 to 0.597 | pending |
| `gemma-4-26b-a4b` | 306 | 0 | 0.591 | 0.519 to 0.669 | 0.591 | 0.518 to 0.667 | pending |
| `deepseek-v4.1-flash` | 108 | 14 | 0.443 | 0.321 to 0.582 | 0.514 | 0.380 to 0.655 | pending |
| `qwen3-235b-a22b-2507` | 237 | 1 | 0.581 | 0.511 to 0.656 | 0.584 | 0.515 to 0.658 | pending |
| `glm-5.3-flash` | 94 | 42 | 0.279 | 0.171 to 0.433 | 0.540 | 0.371 to 0.745 | pending |
| `glm-5.3` | 119 | 100 | 0.084 | 0.032 to 0.192 | 0.785 | 0.631 to 0.950 | pending |
| `muse-glimmer-30b` | 80 | 29 | 0.440 | 0.260 to 0.699 | 0.737 | 0.608 to 0.886 | pending |
| `kimi-k3` | 108 | 5 | 0.674 | 0.547 to 0.784 | 0.710 | 0.601 to 0.812 | pending |
| `gpt-5.4-mini` | 384 | 0 | 0.475 | 0.403 to 0.551 | 0.475 | 0.404 to 0.553 | pending |
| `gemini-3.1-flash-lite` | 305 | 0 | 0.495 | 0.438 to 0.555 | 0.495 | 0.435 to 0.557 | pending |
| `claude-sonnet-5` | 61 | 0 | 0.691 | 0.597 to 0.857 | 0.691 | 0.598 to 0.854 | pending |

**The digit rule.**

| model | flag rate, wrong answers | 95% CI | fully solved traces | flag rate, fully solved | 95% CI | claims checked per trace | traces with a claim |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.167 | 0.104 to 0.244 | 1776 | 0.153 | 0.118 to 0.192 | 2.25 | 0.604 |
| `gemma-4-26b-a4b` | 0.389 | 0.278 to 0.504 | 1869 | 0.135 | 0.102 to 0.173 | 4.12 | 0.673 |
| `deepseek-v4.1-flash` | 0.000 | 0.000 to 0.000 | 2091 | 0.005 | 0.001 to 0.010 | 1.82 | 0.527 |
| `qwen3-235b-a22b-2507` | 0.380 | 0.247 to 0.511 | 1914 | 0.216 | 0.174 to 0.261 | 8.35 | 0.855 |
| `glm-5.3-flash` | 0.053 | 0.014 to 0.103 | 2138 | 0.063 | 0.046 to 0.082 | 4.10 | 0.788 |
| `glm-5.3` | 0.017 | 0.000 to 0.057 | 2109 | 0.066 | 0.051 to 0.083 | 4.44 | 0.843 |
| `muse-glimmer-30b` | 0.125 | 0.067 to 0.208 | 2146 | 0.122 | 0.095 to 0.151 | 3.10 | 0.679 |
| `kimi-k3` | 0.000 | 0.000 to 0.000 | 2099 | 0.018 | 0.012 to 0.024 | 1.68 | 0.471 |
| `gpt-5.4-mini` | 0.172 | 0.087 to 0.262 | 1783 | 0.096 | 0.069 to 0.127 | 2.56 | 0.611 |
| `gemini-3.1-flash-lite` | 0.308 | 0.187 to 0.434 | 1901 | 0.079 | 0.056 to 0.107 | 3.65 | 0.740 |
| `claude-sonnet-5` | 0.164 | 0.031 to 0.291 | 2163 | 0.116 | 0.089 to 0.145 | 4.09 | 0.785 |

## Q4. Consistency within a template

The share of templates fully solved on all 15 instances, on some, and on none: the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) and the 92 others. Intervals are in `results.json`.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.345 | 0.638 | 0.017 | 0.293 | 0.707 | 0.000 |
| `gemma-4-26b-a4b` | 0.466 | 0.500 | 0.034 | 0.500 | 0.478 | 0.022 |
| `deepseek-v4.1-flash` | 0.724 | 0.276 | 0.000 | 0.783 | 0.196 | 0.022 |
| `qwen3-235b-a22b-2507` | 0.466 | 0.517 | 0.017 | 0.500 | 0.467 | 0.033 |
| `glm-5.3-flash` | 0.828 | 0.172 | 0.000 | 0.750 | 0.250 | 0.000 |
| `glm-5.3` | 0.810 | 0.172 | 0.017 | 0.739 | 0.250 | 0.011 |
| `muse-glimmer-30b` | 0.793 | 0.207 | 0.000 | 0.772 | 0.228 | 0.000 |
| `kimi-k3` | 0.707 | 0.293 | 0.000 | 0.717 | 0.283 | 0.000 |
| `gpt-5.4-mini` | 0.397 | 0.586 | 0.017 | 0.370 | 0.576 | 0.054 |
| `gemini-3.1-flash-lite` | 0.672 | 0.276 | 0.052 | 0.576 | 0.402 | 0.022 |
| `claude-sonnet-5` | 0.897 | 0.103 | 0.000 | 0.783 | 0.217 | 0.000 |

## Q5. Paraphrase robustness

Not run yet.

## Sensitivity

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.800 | 0.807 | 0.815 | 0.789 | 0.827 | 0.805 |
| `gemma-4-26b-a4b` | 0.826 | 0.847 | 0.878 | 0.831 | 0.847 | 0.844 |
| `deepseek-v4.1-flash` | 0.937 | 0.941 | 0.942 | 0.929 | 0.947 | 0.939 |
| `qwen3-235b-a22b-2507` | 0.859 | 0.873 | 0.885 | 0.851 | 0.873 | 0.870 |
| `glm-5.3-flash` | 0.950 | 0.954 | 0.956 | 0.950 | 0.972 | 0.953 |
| `glm-5.3` | 0.939 | 0.942 | 0.944 | 0.937 | 0.986 | 0.941 |
| `muse-glimmer-30b` | 0.957 | 0.959 | 0.962 | 0.954 | 0.972 | 0.958 |
| `kimi-k3` | 0.939 | 0.942 | 0.945 | 0.933 | 0.945 | 0.942 |
| `gpt-5.4-mini` | 0.799 | 0.811 | 0.822 | 0.792 | 0.811 | 0.807 |
| `gemini-3.1-flash-lite` | 0.845 | 0.855 | 0.868 | 0.845 | 0.855 | 0.852 |
| `claude-sonnet-5` | 0.960 | 0.967 | 0.973 | 0.961 | 0.967 | 0.966 |
| tau with the headline | 0.954 | 1 | 0.964 | 0.964 | 0.673 | 1.000 |

The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).

## Also reported, not tested

### By branch and level

| model | chemical | civil | electrical | industrial | mechanical | Easy | Intermediate | Advanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.857 | 0.687 | 0.853 | 0.804 | 0.833 | 0.859 | 0.782 | 0.760 |
| `gemma-4-26b-a4b` | 0.736 | 0.876 | 0.858 | 0.873 | 0.894 | 0.897 | 0.843 | 0.771 |
| `deepseek-v4.1-flash` | 0.867 | 0.989 | 0.890 | 0.978 | 0.980 | 0.945 | 0.941 | 0.932 |
| `qwen3-235b-a22b-2507` | 0.824 | 0.931 | 0.868 | 0.869 | 0.871 | 0.900 | 0.865 | 0.839 |
| `glm-5.3-flash` | 0.931 | 0.978 | 0.941 | 0.942 | 0.979 | 0.966 | 0.956 | 0.932 |
| `glm-5.3` | 0.898 | 0.964 | 0.938 | 0.922 | 0.989 | 0.975 | 0.963 | 0.851 |
| `muse-glimmer-30b` | 0.970 | 0.958 | 0.939 | 0.940 | 0.989 | 0.954 | 0.964 | 0.959 |
| `kimi-k3` | 0.910 | 0.980 | 0.896 | 0.980 | 0.947 | 0.955 | 0.937 | 0.929 |
| `gpt-5.4-mini` | 0.733 | 0.831 | 0.833 | 0.836 | 0.821 | 0.857 | 0.806 | 0.740 |
| `gemini-3.1-flash-lite` | 0.774 | 0.864 | 0.840 | 0.858 | 0.937 | 0.911 | 0.869 | 0.734 |
| `claude-sonnet-5` | 0.949 | 0.987 | 0.936 | 0.973 | 0.991 | 0.979 | 0.978 | 0.928 |

### By domain

| domain | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| digital_communications | 0.983 | 0.967 | 0.990 | 0.893 | 0.987 | 0.947 | 0.997 | 0.937 | 0.977 | 0.893 | 0.953 |
| electromagnetics_and_waves | 0.790 | 0.777 | 0.803 | 0.863 | 0.897 | 0.943 | 0.900 | 0.857 | 0.680 | 0.807 | 0.950 |
| fluid_mechanics | 0.847 | 0.893 | 0.983 | 0.913 | 0.990 | 0.987 | 1.000 | 0.973 | 0.880 | 0.950 | 0.987 |
| geotechnical_engineering | 0.527 | 0.913 | 0.987 | 0.947 | 0.987 | 0.980 | 0.920 | 0.993 | 0.913 | 0.960 | 0.980 |
| mechanics_of_materials | 0.807 | 0.827 | 0.973 | 0.843 | 0.977 | 1.000 | 0.990 | 0.927 | 0.753 | 0.867 | 1.000 |
| production_and_inventory | 0.707 | 0.873 | 0.933 | 0.847 | 0.887 | 0.887 | 0.887 | 0.967 | 0.780 | 0.853 | 0.987 |
| quality_and_reliability_control | 0.773 | 0.800 | 1.000 | 0.787 | 0.947 | 0.907 | 0.933 | 0.980 | 0.760 | 0.800 | 0.940 |
| reaction_kinetics | 0.970 | 0.913 | 1.000 | 0.973 | 1.000 | 0.993 | 1.000 | 0.993 | 0.900 | 0.913 | 1.000 |
| signals_and_systems | 0.787 | 0.830 | 0.877 | 0.847 | 0.940 | 0.923 | 0.920 | 0.893 | 0.843 | 0.820 | 0.903 |
| stochastic_operations | 0.933 | 0.947 | 1.000 | 0.973 | 0.993 | 0.973 | 1.000 | 0.993 | 0.967 | 0.920 | 0.993 |
| structural_analysis | 0.727 | 0.907 | 1.000 | 0.947 | 1.000 | 1.000 | 0.980 | 0.987 | 0.827 | 0.913 | 0.993 |
| thermodynamics | 0.783 | 0.522 | 0.906 | 0.728 | 0.900 | 0.761 | 0.944 | 0.944 | 0.694 | 0.678 | 0.878 |
| transport_phenomena | 0.825 | 0.833 | 0.642 | 0.783 | 0.892 | 0.983 | 0.971 | 0.754 | 0.583 | 0.746 | 0.992 |
| vibrations_and_acoustics | 0.847 | 0.963 | 0.983 | 0.857 | 0.970 | 0.980 | 0.977 | 0.940 | 0.830 | 0.993 | 0.987 |
| water_resources | 0.807 | 0.807 | 0.980 | 0.900 | 0.947 | 0.913 | 0.973 | 0.960 | 0.753 | 0.720 | 0.987 |

### By answer type

| answer type | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| array | 0.924 | 0.867 | 0.981 | 0.933 | 0.914 | 0.771 | 0.981 | 0.962 | 0.762 | 0.771 | 1.000 |
| classification | 0.620 | 0.953 | 0.907 | 0.880 | 0.967 | 0.973 | 0.940 | 0.907 | 0.927 | 0.973 | 0.980 |
| multipart | 0.863 | 0.890 | 0.950 | 0.883 | 0.980 | 0.952 | 0.972 | 0.973 | 0.860 | 0.902 | 0.985 |
| scalar | 0.779 | 0.823 | 0.948 | 0.871 | 0.956 | 0.953 | 0.963 | 0.948 | 0.793 | 0.849 | 0.969 |
| symbolic | 0.856 | 0.833 | 0.907 | 0.841 | 0.919 | 0.878 | 0.922 | 0.841 | 0.819 | 0.730 | 0.837 |
| vector | 0.854 | 0.879 | 0.850 | 0.825 | 0.900 | 0.988 | 0.904 | 0.875 | 0.771 | 0.863 | 0.979 |

### Tokens against score

Median completion tokens as billed.

| model | answer score | median tokens | on fully solved | on the rest |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.807 | 1795 | 1612 | 2834 |
| `gemma-4-26b-a4b` | 0.847 | 797 | 768 | 970 |
| `deepseek-v4.1-flash` | 0.941 | 3060 | 2986 | 4537 |
| `qwen3-235b-a22b-2507` | 0.873 | 1013 | 993 | 1272 |
| `glm-5.3-flash` | 0.954 | 2500 | 2370 | 8284 |
| `glm-5.3` | 0.942 | 3270 | 3082 | 32768 |
| `muse-glimmer-30b` | 0.959 | 2988 | 2948 | 4288 |
| `kimi-k3` | 0.942 | 2192 | 2141 | 3726 |
| `gpt-5.4-mini` | 0.811 | 632 | 615 | 761 |
| `gemini-3.1-flash-lite` | 0.855 | 653 | 632 | 796 |
| `claude-sonnet-5` | 0.967 | 1804 | 1783 | 2747 |

### Classification templates: fully solved by the gold label

The pool over-represents the rare labels by design (D-116), so these rates are per label, not a pooled accuracy.

| template | gold label | items | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `critical_depth_froude_classification` | subcritical | 7 | 0.143 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `critical_depth_froude_classification` | supercritical | 8 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `damping_classification` | critically | 5 | 0.200 | 0.600 | 0.600 | 0.400 | 0.200 | 0.800 | 0.400 | 0.200 | 0.400 | 1.000 | 0.400 |
| `damping_classification` | overdamped | 5 | 0.400 | 1.000 | 0.800 | 0.800 | 0.800 | 0.400 | 0.400 | 0.600 | 1.000 | 1.000 | 1.000 |
| `damping_classification` | underdamped | 5 | 0.800 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | laminar | 6 | 0.500 | 0.667 | 0.500 | 0.667 | 1.000 | 1.000 | 1.000 | 0.833 | 0.500 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | transitional | 1 | 0.000 | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | turbulent | 8 | 0.500 | 0.750 | 0.375 | 0.375 | 1.000 | 1.000 | 1.000 | 0.625 | 0.500 | 1.000 | 1.000 |
| `system_properties_memory_causality` | memoryless no, causal no | 4 | 0.750 | 1.000 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 | 0.750 | 1.000 | 0.500 | 1.000 |
| `system_properties_memory_causality` | memoryless no, causal yes | 6 | 1.000 | 1.000 | 0.667 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_properties_memory_causality` | memoryless yes, causal yes | 5 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_property_linearity` | linear | 7 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.857 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_property_linearity` | nonlinear | 8 | 0.375 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.875 | 1.000 | 1.000 | 1.000 |

## Set aside (D-132), outside every comparison

- `qwen3.8-27b`: answer score 0.915 (95% CI 0.878 to 0.947), fully solved 0.906, unusable 89
