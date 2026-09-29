# Results of the full run: the deterministic stack

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring. The answer check, E3 and E4's digit rule only: E5's judge has not run, so its columns are empty. Every interval is 95% and resamples the 150 templates (B = 10,000).

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0.

| model | answer score | 95% CI | fully solved | 95% CI | correct | partial | incorrect | unusable |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.969 | 0.948 to 0.985 | 0.962 | 0.939 to 0.981 | 2165 | 29 | 42 | 14 |
| `claude-sonnet-5` | 0.963 | 0.941 to 0.982 | 0.957 | 0.933 to 0.977 | 2153 | 28 | 69 | 0 |
| `kimi-k3` | 0.962 | 0.941 to 0.979 | 0.956 | 0.934 to 0.976 | 2152 | 26 | 67 | 5 |
| `glm-5.3-flash` | 0.962 | 0.939 to 0.980 | 0.958 | 0.934 to 0.977 | 2155 | 18 | 35 | 42 |
| `muse-glimmer-30b` | 0.957 | 0.934 to 0.977 | 0.952 | 0.928 to 0.973 | 2143 | 21 | 57 | 29 |
| `glm-5.3` | 0.940 | 0.910 to 0.966 | 0.935 | 0.904 to 0.963 | 2104 | 23 | 23 | 100 |
| `qwen3-235b-a22b-2507` | 0.870 | 0.833 to 0.903 | 0.844 | 0.801 to 0.885 | 1899 | 115 | 235 | 1 |
| `gemini-3.1-flash-lite` | 0.857 | 0.811 to 0.898 | 0.848 | 0.802 to 0.891 | 1907 | 41 | 302 | 0 |
| `gemma-4-26b-a4b` | 0.848 | 0.806 to 0.886 | 0.832 | 0.790 to 0.872 | 1871 | 75 | 304 | 0 |
| `gpt-5.4-mini` | 0.829 | 0.786 to 0.869 | 0.814 | 0.771 to 0.856 | 1832 | 67 | 351 | 0 |
| `gpt-oss-20b` | 0.799 | 0.755 to 0.842 | 0.784 | 0.738 to 0.827 | 1763 | 71 | 361 | 55 |

Of the 55 pairs, 32 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The fully-solved check (McNemar exact, Holm) gives the same verdict on 46 of 55: both tests or neither hold, and when both hold, in the same direction. On 9 of the 9 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level test.

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | McNemar p (Holm) | same verdict |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.169 | -0.212 to -0.129 | 0.0055 | -0.179 | -0.222 to -0.137 | 0.0000 | yes |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.164 | -0.206 to -0.126 | 0.0055 | -0.173 | -0.215 to -0.134 | 0.0000 | yes |
| `gpt-oss-20b` | `kimi-k3` | -0.163 | -0.206 to -0.122 | 0.0055 | -0.173 | -0.217 to -0.131 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.162 | -0.204 to -0.123 | 0.0055 | -0.174 | -0.217 to -0.134 | 0.0000 | yes |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.158 | -0.196 to -0.121 | 0.0055 | -0.169 | -0.210 to -0.131 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3` | -0.141 | -0.184 to -0.099 | 0.0055 | -0.152 | -0.195 to -0.109 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.140 | 0.101 to 0.180 | 0.0055 | +0.148 | 0.108 to 0.190 | 0.0000 | yes |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.134 | -0.174 to -0.097 | 0.0055 | -0.143 | -0.183 to -0.105 | 0.0000 | yes |
| `kimi-k3` | `gpt-5.4-mini` | +0.133 | 0.095 to 0.173 | 0.0055 | +0.142 | 0.104 to 0.183 | 0.0000 | yes |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.133 | 0.099 to 0.171 | 0.0055 | +0.144 | 0.108 to 0.181 | 0.0000 | yes |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.128 | 0.094 to 0.165 | 0.0055 | +0.138 | 0.103 to 0.175 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.120 | -0.157 to -0.087 | 0.0055 | -0.131 | -0.169 to -0.096 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.115 | -0.150 to -0.082 | 0.0055 | -0.125 | -0.162 to -0.092 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.114 | -0.150 to -0.081 | 0.0055 | -0.125 | -0.161 to -0.091 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.114 | -0.147 to -0.082 | 0.0055 | -0.126 | -0.163 to -0.093 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.112 | 0.075 to 0.151 | 0.0055 | +0.115 | 0.078 to 0.155 | 0.0000 | yes |
| `glm-5.3` | `gpt-5.4-mini` | +0.111 | 0.079 to 0.146 | 0.0055 | +0.121 | 0.087 to 0.157 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.109 | -0.144 to -0.076 | 0.0055 | -0.121 | -0.160 to -0.085 | 0.0000 | yes |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.106 | -0.145 to -0.072 | 0.0055 | -0.109 | -0.148 to -0.074 | 0.0000 | yes |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.106 | 0.072 to 0.143 | 0.0055 | +0.109 | 0.073 to 0.148 | 0.0000 | yes |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.105 | 0.072 to 0.141 | 0.0055 | +0.110 | 0.076 to 0.148 | 0.0000 | yes |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.100 | 0.066 to 0.138 | 0.0055 | +0.105 | 0.068 to 0.144 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.099 | 0.071 to 0.131 | 0.0055 | +0.118 | 0.084 to 0.157 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.094 | -0.124 to -0.066 | 0.0055 | -0.113 | -0.151 to -0.078 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.093 | -0.123 to -0.065 | 0.0055 | -0.112 | -0.149 to -0.079 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.092 | -0.121 to -0.065 | 0.0055 | -0.114 | -0.151 to -0.080 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.092 | -0.126 to -0.061 | 0.0055 | -0.104 | -0.138 to -0.072 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.088 | -0.117 to -0.060 | 0.0055 | -0.108 | -0.146 to -0.074 | 0.0000 | yes |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.084 | 0.052 to 0.118 | 0.0055 | +0.088 | 0.056 to 0.122 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.071 | -0.101 to -0.042 | 0.0055 | -0.091 | -0.129 to -0.056 | 0.0000 | yes |
| `glm-5.3-flash` | `glm-5.3` | +0.022 | 0.009 to 0.037 | 0.0300 | +0.023 | 0.009 to 0.038 | 0.0000 | yes |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.070 | -0.113 to -0.028 | 0.0408 | -0.060 | -0.108 to -0.012 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.028 | 0.010 to 0.051 | 0.0897 | +0.027 | 0.007 to 0.050 | 0.0000 | no |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.040 | 0.009 to 0.073 | 0.3058 | +0.030 | -0.011 to 0.068 | 0.0178 | no |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.057 | -0.103 to -0.011 | 0.3591 | -0.064 | -0.110 to -0.018 | 0.0000 | no |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.049 | -0.089 to -0.009 | 0.3600 | -0.048 | -0.090 to -0.008 | 0.0000 | no |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.028 | -0.055 to -0.002 | 0.7295 | -0.033 | -0.061 to -0.007 | 0.0002 | no |
| `glm-5.3` | `claude-sonnet-5` | -0.023 | -0.048 to -0.003 | 0.7343 | -0.022 | -0.046 to -0.002 | 0.0002 | no |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.030 | -0.074 to 0.014 | 1.0000 | -0.031 | -0.077 to 0.013 | 0.0352 | no |
| `glm-5.3` | `kimi-k3` | -0.022 | -0.048 to 0.001 | 1.0000 | -0.021 | -0.047 to 0.002 | 0.0014 | no |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.021 | -0.051 to 0.008 | 1.0000 | -0.012 | -0.048 to 0.026 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.019 | -0.011 to 0.049 | 1.0000 | +0.017 | -0.015 to 0.048 | 0.6023 | yes |
| `glm-5.3` | `muse-glimmer-30b` | -0.017 | -0.041 to 0.005 | 1.0000 | -0.017 | -0.044 to 0.005 | 0.0206 | no |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.013 | -0.015 to 0.043 | 1.0000 | -0.004 | -0.042 to 0.034 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.012 | -0.008 to 0.030 | 1.0000 | +0.010 | -0.011 to 0.029 | 0.5819 | yes |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.008 | -0.036 to 0.018 | 1.0000 | -0.016 | -0.043 to 0.010 | 0.5650 | yes |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.007 | -0.004 to 0.019 | 1.0000 | +0.004 | -0.009 to 0.018 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.006 | -0.004 to 0.018 | 1.0000 | +0.006 | -0.005 to 0.019 | 1.0000 | yes |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.006 | -0.026 to 0.014 | 1.0000 | -0.004 | -0.027 to 0.017 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.006 | -0.006 to 0.018 | 1.0000 | +0.005 | -0.008 to 0.020 | 1.0000 | yes |
| `muse-glimmer-30b` | `kimi-k3` | -0.005 | -0.027 to 0.016 | 1.0000 | -0.004 | -0.026 to 0.018 | 1.0000 | yes |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.005 | -0.014 to 0.023 | 1.0000 | +0.005 | -0.015 to 0.024 | 1.0000 | yes |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.001 | -0.018 to 0.013 | 1.0000 | +0.001 | -0.016 to 0.016 | 1.0000 | yes |
| `kimi-k3` | `claude-sonnet-5` | -0.001 | -0.015 to 0.012 | 1.0000 | -0.000 | -0.016 to 0.014 | 1.0000 | yes |
| `glm-5.3-flash` | `kimi-k3` | -0.000 | -0.016 to 0.016 | 1.0000 | +0.001 | -0.015 to 0.017 | 1.0000 | yes |

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier; p from permuting the tier labels, Holm across the eleven models. Sigma is the within-tier SD of the template-mean score; the detectable gap, 2.80 x sigma x sqrt(1/58 + 1/34), is the smallest this test finds at 80% power.

| model | Easy | Advanced | gap | 95% CI | p (Holm) | sigma | detectable gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.885 | 0.702 | +0.183 | 0.068 to 0.301 | 0.0154 | 0.258 | 0.156 |
| `gemma-4-26b-a4b` | 0.906 | 0.758 | +0.148 | 0.049 to 0.255 | 0.0161 | 0.228 | 0.138 |
| `deepseek-v4.1-flash` | 0.987 | 0.915 | +0.072 | 0.009 to 0.153 | 0.0474 | 0.136 | 0.082 |
| `qwen3-235b-a22b-2507` | 0.905 | 0.825 | +0.080 | -0.016 to 0.179 | 0.1856 | 0.214 | 0.129 |
| `glm-5.3-flash` | 0.983 | 0.914 | +0.070 | 0.011 to 0.144 | 0.0474 | 0.130 | 0.079 |
| `glm-5.3` | 0.975 | 0.842 | +0.133 | 0.040 to 0.239 | 0.0154 | 0.192 | 0.116 |
| `muse-glimmer-30b` | 0.966 | 0.919 | +0.047 | -0.013 to 0.111 | 0.1856 | 0.141 | 0.085 |
| `kimi-k3` | 0.983 | 0.916 | +0.068 | 0.008 to 0.142 | 0.0720 | 0.134 | 0.081 |
| `gpt-5.4-mini` | 0.897 | 0.706 | +0.191 | 0.073 to 0.317 | 0.0154 | 0.254 | 0.153 |
| `gemini-3.1-flash-lite` | 0.921 | 0.726 | +0.194 | 0.075 to 0.325 | 0.0154 | 0.262 | 0.159 |
| `claude-sonnet-5` | 0.979 | 0.912 | +0.067 | -0.003 to 0.150 | 0.1044 | 0.155 | 0.094 |

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 70 items with no milestones (all 15 of 3 templates and some items of 12 more); 52 templates have an item with at most one milestone (46 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The digit rule's flag counts a trace when any step is flagged. Against the experts, on fully solved traces, it had precision 0.750 and recall 0.320 (SCORER_VALIDATION.md), so its rate is not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 416 | 55 | 0.406 | 0.341 to 0.475 | 0.467 | 0.401 to 0.538 | pending |
| `gemma-4-26b-a4b` | 304 | 0 | 0.591 | 0.520 to 0.670 | 0.591 | 0.519 to 0.669 | pending |
| `deepseek-v4.1-flash` | 56 | 14 | 0.404 | 0.269 to 0.591 | 0.565 | 0.410 to 0.721 | pending |
| `qwen3-235b-a22b-2507` | 236 | 1 | 0.591 | 0.520 to 0.669 | 0.594 | 0.521 to 0.673 | pending |
| `glm-5.3-flash` | 77 | 42 | 0.265 | 0.137 to 0.444 | 0.662 | 0.518 to 0.841 | pending |
| `glm-5.3` | 123 | 100 | 0.104 | 0.046 to 0.216 | 0.755 | 0.614 to 0.936 | pending |
| `muse-glimmer-30b` | 86 | 29 | 0.410 | 0.262 to 0.612 | 0.653 | 0.571 to 0.778 | pending |
| `kimi-k3` | 72 | 5 | 0.586 | 0.454 to 0.710 | 0.635 | 0.535 to 0.742 | pending |
| `gpt-5.4-mini` | 351 | 0 | 0.438 | 0.368 to 0.518 | 0.438 | 0.367 to 0.518 | pending |
| `gemini-3.1-flash-lite` | 302 | 0 | 0.514 | 0.448 to 0.584 | 0.514 | 0.449 to 0.581 | pending |
| `claude-sonnet-5` | 69 | 0 | 0.682 | 0.593 to 0.816 | 0.682 | 0.595 to 0.811 | pending |

**The digit rule.**

| model | flag rate, wrong answers | 95% CI | fully solved traces | flag rate, fully solved | 95% CI | claims checked per trace | traces with a claim |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.178 | 0.122 to 0.242 | 1763 | 0.151 | 0.114 to 0.190 | 2.25 | 0.604 |
| `gemma-4-26b-a4b` | 0.408 | 0.297 to 0.525 | 1871 | 0.132 | 0.099 to 0.169 | 4.12 | 0.673 |
| `deepseek-v4.1-flash` | 0.018 | 0.000 to 0.057 | 2165 | 0.005 | 0.001 to 0.010 | 1.82 | 0.527 |
| `qwen3-235b-a22b-2507` | 0.386 | 0.255 to 0.517 | 1899 | 0.217 | 0.174 to 0.262 | 8.35 | 0.855 |
| `glm-5.3-flash` | 0.078 | 0.030 to 0.146 | 2155 | 0.062 | 0.045 to 0.081 | 4.10 | 0.788 |
| `glm-5.3` | 0.016 | 0.000 to 0.055 | 2104 | 0.065 | 0.050 to 0.082 | 4.44 | 0.843 |
| `muse-glimmer-30b` | 0.116 | 0.055 to 0.202 | 2143 | 0.123 | 0.096 to 0.152 | 3.10 | 0.679 |
| `kimi-k3` | 0.000 | 0.000 to 0.000 | 2152 | 0.017 | 0.012 to 0.023 | 1.68 | 0.471 |
| `gpt-5.4-mini` | 0.182 | 0.093 to 0.284 | 1832 | 0.095 | 0.069 to 0.125 | 2.56 | 0.611 |
| `gemini-3.1-flash-lite` | 0.308 | 0.187 to 0.435 | 1907 | 0.079 | 0.055 to 0.106 | 3.65 | 0.740 |
| `claude-sonnet-5` | 0.145 | 0.028 to 0.264 | 2153 | 0.116 | 0.089 to 0.145 | 4.09 | 0.785 |

## Q4. Consistency within a template

The share of templates fully solved on all 15 instances, on some, and on none: the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) and the 92 others. Intervals are in `results.json`.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.328 | 0.603 | 0.069 | 0.283 | 0.707 | 0.011 |
| `gemma-4-26b-a4b` | 0.466 | 0.500 | 0.034 | 0.478 | 0.500 | 0.022 |
| `deepseek-v4.1-flash` | 0.862 | 0.138 | 0.000 | 0.826 | 0.174 | 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.534 | 0.017 | 0.489 | 0.478 | 0.033 |
| `glm-5.3-flash` | 0.879 | 0.121 | 0.000 | 0.772 | 0.228 | 0.000 |
| `glm-5.3` | 0.793 | 0.190 | 0.017 | 0.728 | 0.261 | 0.011 |
| `muse-glimmer-30b` | 0.828 | 0.172 | 0.000 | 0.772 | 0.228 | 0.000 |
| `kimi-k3` | 0.810 | 0.190 | 0.000 | 0.750 | 0.250 | 0.000 |
| `gpt-5.4-mini` | 0.431 | 0.534 | 0.034 | 0.391 | 0.576 | 0.033 |
| `gemini-3.1-flash-lite` | 0.672 | 0.276 | 0.052 | 0.576 | 0.402 | 0.022 |
| `claude-sonnet-5` | 0.879 | 0.121 | 0.000 | 0.783 | 0.217 | 0.000 |

## Q5. Paraphrase robustness

Not run yet.

## Sensitivity

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.789 | 0.799 | 0.811 | 0.784 | 0.819 | 0.797 |
| `gemma-4-26b-a4b` | 0.826 | 0.848 | 0.881 | 0.832 | 0.848 | 0.845 |
| `deepseek-v4.1-flash` | 0.964 | 0.969 | 0.971 | 0.962 | 0.975 | 0.969 |
| `qwen3-235b-a22b-2507` | 0.856 | 0.870 | 0.882 | 0.844 | 0.870 | 0.866 |
| `glm-5.3-flash` | 0.958 | 0.962 | 0.964 | 0.958 | 0.980 | 0.961 |
| `glm-5.3` | 0.937 | 0.940 | 0.942 | 0.935 | 0.984 | 0.939 |
| `muse-glimmer-30b` | 0.954 | 0.957 | 0.960 | 0.952 | 0.970 | 0.956 |
| `kimi-k3` | 0.959 | 0.962 | 0.965 | 0.956 | 0.964 | 0.962 |
| `gpt-5.4-mini` | 0.815 | 0.829 | 0.844 | 0.814 | 0.829 | 0.825 |
| `gemini-3.1-flash-lite` | 0.846 | 0.857 | 0.871 | 0.848 | 0.857 | 0.854 |
| `claude-sonnet-5` | 0.956 | 0.963 | 0.969 | 0.957 | 0.963 | 0.962 |
| tau with the headline | 0.927 | 1 | 0.964 | 0.891 | 0.600 | 1.000 |

The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).

## Also reported, not tested

### By branch and level

| model | chemical | civil | electrical | industrial | mechanical | Easy | Intermediate | Advanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.807 | 0.662 | 0.869 | 0.816 | 0.843 | 0.885 | 0.771 | 0.702 |
| `gemma-4-26b-a4b` | 0.734 | 0.873 | 0.866 | 0.873 | 0.894 | 0.906 | 0.843 | 0.758 |
| `deepseek-v4.1-flash` | 0.956 | 0.989 | 0.924 | 0.980 | 0.994 | 0.987 | 0.982 | 0.915 |
| `qwen3-235b-a22b-2507` | 0.824 | 0.931 | 0.847 | 0.869 | 0.877 | 0.905 | 0.861 | 0.825 |
| `glm-5.3-flash` | 0.947 | 0.976 | 0.954 | 0.942 | 0.990 | 0.983 | 0.968 | 0.914 |
| `glm-5.3` | 0.898 | 0.962 | 0.926 | 0.922 | 0.993 | 0.975 | 0.963 | 0.842 |
| `muse-glimmer-30b` | 0.963 | 0.958 | 0.937 | 0.940 | 0.988 | 0.966 | 0.971 | 0.919 |
| `kimi-k3` | 0.962 | 0.980 | 0.920 | 0.980 | 0.969 | 0.983 | 0.968 | 0.916 |
| `gpt-5.4-mini` | 0.796 | 0.820 | 0.858 | 0.824 | 0.848 | 0.897 | 0.834 | 0.706 |
| `gemini-3.1-flash-lite` | 0.776 | 0.862 | 0.851 | 0.858 | 0.937 | 0.921 | 0.869 | 0.726 |
| `claude-sonnet-5` | 0.949 | 0.987 | 0.912 | 0.973 | 0.994 | 0.979 | 0.978 | 0.912 |

### By domain

| domain | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| digital_communications | 0.930 | 0.883 | 0.923 | 0.840 | 0.967 | 0.923 | 0.947 | 0.897 | 0.923 | 0.877 | 0.913 |
| electromagnetics_and_waves | 0.900 | 0.887 | 0.990 | 0.877 | 0.977 | 0.947 | 0.997 | 0.980 | 0.817 | 0.870 | 0.953 |
| fluid_mechanics | 0.880 | 0.893 | 0.997 | 0.913 | 0.987 | 0.987 | 0.987 | 0.987 | 0.947 | 0.950 | 0.987 |
| geotechnical_engineering | 0.513 | 0.913 | 0.987 | 0.947 | 0.987 | 0.980 | 0.920 | 0.993 | 0.913 | 0.960 | 0.980 |
| mechanics_of_materials | 0.800 | 0.827 | 0.993 | 0.857 | 0.987 | 1.000 | 0.990 | 0.927 | 0.803 | 0.867 | 1.000 |
| production_and_inventory | 0.767 | 0.873 | 0.940 | 0.847 | 0.887 | 0.887 | 0.887 | 0.967 | 0.780 | 0.853 | 0.987 |
| quality_and_reliability_control | 0.773 | 0.800 | 1.000 | 0.787 | 0.947 | 0.907 | 0.933 | 0.980 | 0.760 | 0.800 | 0.940 |
| reaction_kinetics | 0.960 | 0.913 | 0.993 | 0.973 | 1.000 | 0.993 | 1.000 | 0.993 | 0.900 | 0.913 | 1.000 |
| signals_and_systems | 0.777 | 0.827 | 0.860 | 0.823 | 0.920 | 0.907 | 0.867 | 0.883 | 0.833 | 0.807 | 0.870 |
| stochastic_operations | 0.907 | 0.947 | 1.000 | 0.973 | 0.993 | 0.973 | 1.000 | 0.993 | 0.933 | 0.920 | 0.993 |
| structural_analysis | 0.720 | 0.907 | 1.000 | 0.947 | 1.000 | 1.000 | 0.980 | 0.987 | 0.827 | 0.913 | 0.993 |
| thermodynamics | 0.656 | 0.517 | 0.900 | 0.728 | 0.872 | 0.761 | 0.911 | 0.922 | 0.667 | 0.678 | 0.878 |
| transport_phenomena | 0.842 | 0.838 | 0.992 | 0.783 | 0.992 | 0.983 | 0.996 | 0.983 | 0.858 | 0.750 | 0.992 |
| vibrations_and_acoustics | 0.850 | 0.963 | 0.993 | 0.860 | 0.997 | 0.993 | 0.987 | 0.993 | 0.793 | 0.993 | 0.997 |
| water_resources | 0.753 | 0.800 | 0.980 | 0.900 | 0.940 | 0.907 | 0.973 | 0.960 | 0.720 | 0.713 | 0.987 |

### By answer type

| answer type | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| array | 0.886 | 0.867 | 0.962 | 0.924 | 0.914 | 0.762 | 0.971 | 0.962 | 0.752 | 0.771 | 0.990 |
| classification | 0.680 | 0.953 | 0.960 | 0.887 | 0.993 | 1.000 | 0.960 | 0.960 | 0.973 | 0.973 | 1.000 |
| multipart | 0.870 | 0.891 | 0.980 | 0.892 | 0.978 | 0.950 | 0.972 | 0.980 | 0.887 | 0.906 | 0.979 |
| scalar | 0.762 | 0.822 | 0.978 | 0.871 | 0.963 | 0.952 | 0.960 | 0.971 | 0.810 | 0.849 | 0.969 |
| symbolic | 0.789 | 0.778 | 0.822 | 0.756 | 0.893 | 0.852 | 0.822 | 0.793 | 0.767 | 0.700 | 0.785 |
| vector | 0.946 | 0.971 | 1.000 | 0.833 | 0.988 | 0.988 | 1.000 | 0.988 | 0.858 | 0.925 | 0.983 |

### Tokens against score

Median completion tokens as billed.

| model | answer score | median tokens | on fully solved | on the rest |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.799 | 1795 | 1588 | 2931 |
| `gemma-4-26b-a4b` | 0.848 | 797 | 763 | 987 |
| `deepseek-v4.1-flash` | 0.969 | 3060 | 2994 | 7705 |
| `qwen3-235b-a22b-2507` | 0.870 | 1013 | 980 | 1300 |
| `glm-5.3-flash` | 0.962 | 2500 | 2396 | 25561 |
| `glm-5.3` | 0.940 | 3270 | 3076 | 32768 |
| `muse-glimmer-30b` | 0.957 | 2988 | 2958 | 4300 |
| `kimi-k3` | 0.962 | 2192 | 2148 | 4513 |
| `gpt-5.4-mini` | 0.829 | 632 | 609 | 811 |
| `gemini-3.1-flash-lite` | 0.857 | 653 | 632 | 801 |
| `claude-sonnet-5` | 0.963 | 1804 | 1783 | 2775 |

### Classification templates: fully solved by the gold label

The pool over-represents the rare labels by design (D-116), so these rates are per label, not a pooled accuracy.

| template | gold label | items | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `critical_depth_froude_classification` | subcritical | 7 | 0.143 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `critical_depth_froude_classification` | supercritical | 8 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `damping_classification` | critically | 5 | 0.400 | 0.600 | 0.600 | 0.400 | 0.800 | 1.000 | 0.400 | 0.800 | 0.400 | 1.000 | 1.000 |
| `damping_classification` | overdamped | 5 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `damping_classification` | underdamped | 5 | 0.800 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | laminar | 6 | 1.000 | 0.667 | 1.000 | 0.667 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | transitional | 1 | 0.000 | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 1.000 | 1.000 |
| `reynolds_number_flow_regime` | turbulent | 8 | 0.625 | 0.750 | 0.750 | 0.375 | 1.000 | 1.000 | 1.000 | 0.875 | 0.875 | 1.000 | 1.000 |
| `system_properties_memory_causality` | memoryless no, causal no | 4 | 0.750 | 1.000 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 | 0.750 | 1.000 | 0.500 | 1.000 |
| `system_properties_memory_causality` | memoryless no, causal yes | 6 | 1.000 | 1.000 | 0.667 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_properties_memory_causality` | memoryless yes, causal yes | 5 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_property_linearity` | linear | 7 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.857 | 1.000 | 1.000 | 1.000 | 1.000 |
| `system_property_linearity` | nonlinear | 8 | 0.375 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.875 | 1.000 | 1.000 | 1.000 |

## Set aside (D-132), outside every comparison

- `qwen3.8-27b`: answer score 0.922 (95% CI 0.887 to 0.953), fully solved 0.917, unusable 89
