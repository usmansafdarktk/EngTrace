# Results of the full run: the deterministic stack

**A record, not the results:** computed from `scores/main_pre_d137`. The same traces scored by the evaluators as they were before D-137, D-138 and D-139 (score.py at a72399b); kept to show what those fixes changed. Regenerated at bf4a43b with the corrected analysis (D-146 to D-149): its CONFIG names commit 9915400, at which score.py did not yet exist and answer.py and milestones.py stood at 3a7f247; the rows reproduce with that code (D-144). Its store predates the sensitivity fields, so those columns are blank.

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring, and the corrections and additions made after the first results were read are labelled where they appear (D-146 to D-149). The answer check, E3 and E4's digit rule; E5's judge has not run, so its columns are empty. Every interval is 95% and resamples templates (B = 10,000): all 150 for a model's score, and the templates with a qualifying trace for a rate over a subset of traces. Every resampling test draws 100,000 permutations, so its Holm-adjusted floor over 55 pairs is 0.0006.

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0. SD within: the mean over templates of the SD of the score across a template's 15 items; SD between: the SD of the 150 template means (the per-template variability the July rebuttal promised). Unusable is empty plus unreadable; capped: answered rows that stopped at the output cap and are scored on what they state (D-117).

| model | answer score | 95% CI | SD within | SD between | fully solved | 95% CI | correct | partial | incorrect | unusable (empty + unreadable) | capped, scored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `claude-sonnet-5` | 0.967 | 0.947 to 0.983 | 0.052 | 0.113 | 0.961 | 0.940 to 0.979 | 2163 | 26 | 61 | 0 (0 + 0) | 0 |
| `muse-glimmer-30b` | 0.959 | 0.936 to 0.978 | 0.064 | 0.129 | 0.954 | 0.929 to 0.975 | 2146 | 24 | 51 | 29 (29 + 0) | 0 |
| `glm-5.3-flash` | 0.954 | 0.931 to 0.974 | 0.071 | 0.134 | 0.950 | 0.926 to 0.971 | 2138 | 18 | 52 | 42 (42 + 0) | 1 |
| `kimi-k3` | 0.942 | 0.919 to 0.963 | 0.091 | 0.139 | 0.933 | 0.906 to 0.956 | 2099 | 43 | 103 | 5 (5 + 0) | 4 |
| `glm-5.3` | 0.942 | 0.912 to 0.968 | 0.066 | 0.177 | 0.937 | 0.906 to 0.964 | 2109 | 22 | 19 | 100 (100 + 0) | 6 |
| `deepseek-v4.1-flash` | 0.941 | 0.912 to 0.965 | 0.071 | 0.166 | 0.929 | 0.896 to 0.957 | 2091 | 51 | 94 | 14 (14 + 0) | 1 |
| `qwen3-235b-a22b-2507` | 0.873 | 0.836 to 0.906 | 0.158 | 0.220 | 0.851 | 0.809 to 0.890 | 1914 | 99 | 236 | 1 (0 + 1) | 1 |
| `gemini-3.1-flash-lite` | 0.855 | 0.808 to 0.897 | 0.125 | 0.274 | 0.845 | 0.798 to 0.888 | 1901 | 44 | 305 | 0 (0 + 0) | 0 |
| `gemma-4-26b-a4b` | 0.847 | 0.804 to 0.886 | 0.162 | 0.253 | 0.831 | 0.787 to 0.872 | 1869 | 75 | 306 | 0 (0 + 0) | 1 |
| `gpt-5.4-mini` | 0.811 | 0.766 to 0.853 | 0.198 | 0.273 | 0.792 | 0.745 to 0.836 | 1783 | 83 | 384 | 0 (0 + 0) | 0 |
| `gpt-oss-20b` | 0.807 | 0.764 to 0.847 | 0.228 | 0.262 | 0.789 | 0.745 to 0.832 | 1776 | 79 | 340 | 55 (53 + 2) | 4 |

Of the 55 pairs, 31 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The same template-level test on the fully-solved rate gives the same verdict on 53 of 55. McNemar's exact test on the paired item verdicts (Holm) gives the same verdict on 42 of 55; on 13 of the 13 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level tests. 7 answered rows across the roster ended with a finish reason other than stop or length (a provider fault inside a 200) and are scored on what they state; the harness now retries such a reply (D-148).

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | template p (Holm) | McNemar p (Holm) | same verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.160 | -0.202 to -0.119 | 0.0005 | -0.172 | -0.215 to -0.130 | 0.0005 | 0.0000 | yes |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.156 | -0.198 to -0.117 | 0.0005 | -0.169 | -0.214 to -0.127 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.152 | -0.190 to -0.117 | 0.0005 | -0.164 | -0.204 to -0.128 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.148 | 0.113 to 0.186 | 0.0005 | +0.161 | 0.125 to 0.200 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.147 | -0.189 to -0.108 | 0.0005 | -0.161 | -0.204 to -0.119 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.143 | 0.109 to 0.181 | 0.0005 | +0.158 | 0.121 to 0.197 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `kimi-k3` | -0.136 | -0.179 to -0.094 | 0.0005 | -0.144 | -0.189 to -0.100 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3` | -0.135 | -0.181 to -0.090 | 0.0005 | -0.148 | -0.194 to -0.100 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.134 | -0.178 to -0.090 | 0.0005 | -0.140 | -0.185 to -0.094 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gpt-5.4-mini` | +0.132 | 0.096 to 0.168 | 0.0005 | +0.140 | 0.104 to 0.180 | 0.0005 | 0.0000 | yes |
| `glm-5.3` | `gpt-5.4-mini` | +0.131 | 0.094 to 0.170 | 0.0005 | +0.145 | 0.105 to 0.187 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.130 | 0.093 to 0.168 | 0.0005 | +0.137 | 0.098 to 0.176 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.120 | -0.156 to -0.086 | 0.0005 | -0.131 | -0.169 to -0.095 | 0.0005 | 0.0000 | yes |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.112 | -0.152 to -0.077 | 0.0005 | -0.116 | -0.156 to -0.080 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.112 | -0.148 to -0.079 | 0.0005 | -0.123 | -0.163 to -0.086 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.107 | -0.141 to -0.074 | 0.0005 | -0.120 | -0.158 to -0.084 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.104 | 0.068 to 0.143 | 0.0005 | +0.109 | 0.071 to 0.149 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.100 | 0.067 to 0.135 | 0.0005 | +0.105 | 0.070 to 0.144 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.095 | -0.133 to -0.060 | 0.0005 | -0.102 | -0.141 to -0.064 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.095 | -0.131 to -0.061 | 0.0005 | -0.107 | -0.145 to -0.072 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.094 | -0.125 to -0.066 | 0.0005 | -0.111 | -0.148 to -0.076 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.093 | -0.133 to -0.056 | 0.0005 | -0.099 | -0.140 to -0.060 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.088 | 0.054 to 0.125 | 0.0005 | +0.088 | 0.052 to 0.128 | 0.0006 | 0.0000 | yes |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.088 | 0.055 to 0.122 | 0.0005 | +0.092 | 0.060 to 0.129 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.086 | -0.117 to -0.058 | 0.0005 | -0.103 | -0.141 to -0.069 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.086 | 0.047 to 0.127 | 0.0005 | +0.084 | 0.044 to 0.127 | 0.0016 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.082 | -0.111 to -0.054 | 0.0005 | -0.100 | -0.135 to -0.067 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.070 | -0.100 to -0.042 | 0.0005 | -0.082 | -0.120 to -0.048 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.070 | -0.101 to -0.040 | 0.0005 | -0.087 | -0.125 to -0.052 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.068 | 0.038 to 0.100 | 0.0008 | +0.079 | 0.040 to 0.118 | 0.0016 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.062 | 0.031 to 0.094 | 0.0045 | +0.058 | 0.018 to 0.097 | 0.1075 | 0.0000 | no |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.044 | -0.071 to -0.017 | 0.0516 | -0.052 | -0.084 to -0.022 | 0.0262 | 0.0000 | no |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.066 | -0.111 to -0.023 | 0.0828 | -0.061 | -0.110 to -0.012 | 0.3678 | 0.0000 | yes |
| `glm-5.3` | `claude-sonnet-5` | -0.025 | -0.049 to -0.006 | 0.4466 | -0.024 | -0.048 to -0.005 | 0.6390 | 0.0000 | yes |
| `kimi-k3` | `claude-sonnet-5` | -0.025 | -0.046 to -0.005 | 0.4466 | -0.028 | -0.054 to -0.006 | 0.4336 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.036 | 0.005 to 0.068 | 0.5040 | +0.038 | 0.005 to 0.071 | 0.5462 | 0.0003 | yes |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.026 | -0.052 to -0.003 | 0.6646 | -0.032 | -0.063 to -0.004 | 0.6286 | 0.0000 | yes |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.048 | -0.097 to 0.002 | 1.0000 | -0.056 | -0.106 to -0.004 | 0.6572 | 0.0000 | yes |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.040 | -0.084 to 0.003 | 1.0000 | -0.041 | -0.087 to 0.003 | 1.0000 | 0.0006 | yes |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.025 | -0.055 to 0.004 | 1.0000 | -0.020 | -0.057 to 0.017 | 1.0000 | 0.2556 | yes |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | -0.018 | -0.045 to 0.006 | 1.0000 | -0.024 | -0.053 to 0.002 | 1.0000 | 0.0006 | yes |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.018 | -0.010 to 0.048 | 1.0000 | +0.006 | -0.032 to 0.042 | 1.0000 | 1.0000 | yes |
| `glm-5.3` | `muse-glimmer-30b` | -0.017 | -0.044 to 0.008 | 1.0000 | -0.016 | -0.046 to 0.010 | 1.0000 | 0.0499 | yes |
| `muse-glimmer-30b` | `kimi-k3` | +0.017 | -0.007 to 0.040 | 1.0000 | +0.021 | -0.003 to 0.046 | 1.0000 | 0.0066 | yes |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.014 | -0.033 to 0.005 | 1.0000 | -0.021 | -0.045 to 0.001 | 1.0000 | 0.0024 | yes |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.013 | -0.032 to 0.004 | 1.0000 | -0.011 | -0.032 to 0.008 | 1.0000 | 0.3462 | yes |
| `glm-5.3-flash` | `glm-5.3` | +0.012 | -0.007 to 0.031 | 1.0000 | +0.013 | -0.007 to 0.033 | 1.0000 | 0.1553 | yes |
| `glm-5.3-flash` | `kimi-k3` | +0.012 | -0.007 to 0.031 | 1.0000 | +0.017 | -0.004 to 0.039 | 1.0000 | 0.0364 | yes |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.008 | -0.031 to 0.015 | 1.0000 | -0.008 | -0.033 to 0.018 | 1.0000 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.007 | -0.035 to 0.020 | 1.0000 | -0.014 | -0.043 to 0.015 | 1.0000 | 0.5711 | yes |
| `glm-5.3-flash` | `muse-glimmer-30b` | -0.005 | -0.023 to 0.013 | 1.0000 | -0.004 | -0.022 to 0.015 | 1.0000 | 1.0000 | yes |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.004 | -0.049 to 0.041 | 1.0000 | -0.003 | -0.049 to 0.042 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.002 | -0.016 to 0.014 | 1.0000 | -0.004 | -0.020 to 0.015 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.002 | -0.032 to 0.029 | 1.0000 | -0.008 | -0.042 to 0.025 | 1.0000 | 1.0000 | yes |
| `glm-5.3` | `kimi-k3` | -0.000 | -0.030 to 0.028 | 1.0000 | +0.004 | -0.028 to 0.035 | 1.0000 | 1.0000 | yes |

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier. **Corrected after the first results (D-146):** the planned test permuted the tier labels on the raw difference, which is liberal when the smaller tier has the larger spread, and Advanced template means spread two to four times as widely as Easy ones (the two SD columns); its count of models also moved with the seed. The count now rests on Welch's t-test, Holm across the eleven models: **0 of 11** hold, against 2 under the planned permutation, printed beside it. The detectable gap at 80% power is given as the plan defined it (2.80 x pooled sigma x sqrt(1/58 + 1/34)) and at the strictest Holm step from the Welch standard error.

| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | planned p (Holm) | SD Easy | SD Adv | detectable, planned | detectable, Holm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.859 | 0.760 | +0.099 | -0.012 to 0.217 | 0.7061 | 0.5798 | 0.240 | 0.295 | 0.158 | 0.219 |
| `gemma-4-26b-a4b` | 0.897 | 0.771 | +0.127 | 0.025 to 0.233 | 0.2115 | 0.1328 | 0.217 | 0.270 | 0.144 | 0.200 |
| `deepseek-v4.1-flash` | 0.945 | 0.932 | +0.013 | -0.053 to 0.090 | 1.0000 | 1.0000 | 0.150 | 0.177 | 0.097 | 0.133 |
| `qwen3-235b-a22b-2507` | 0.900 | 0.839 | +0.061 | -0.033 to 0.157 | 1.0000 | 0.9619 | 0.189 | 0.247 | 0.128 | 0.181 |
| `glm-5.3-flash` | 0.966 | 0.932 | +0.033 | -0.018 to 0.090 | 1.0000 | 0.9619 | 0.112 | 0.140 | 0.074 | 0.104 |
| `glm-5.3` | 0.975 | 0.851 | +0.124 | 0.033 to 0.230 | 0.2115 | 0.0232 | 0.082 | 0.294 | 0.115 | 0.190 |
| `muse-glimmer-30b` | 0.954 | 0.959 | -0.005 | -0.058 to 0.048 | 1.0000 | 1.0000 | 0.147 | 0.113 | 0.082 | 0.101 |
| `kimi-k3` | 0.955 | 0.929 | +0.026 | -0.031 to 0.086 | 1.0000 | 1.0000 | 0.134 | 0.141 | 0.083 | 0.110 |
| `gpt-5.4-mini` | 0.857 | 0.740 | +0.117 | 0.001 to 0.239 | 0.5049 | 0.3366 | 0.232 | 0.312 | 0.160 | 0.226 |
| `gemini-3.1-flash-lite` | 0.911 | 0.734 | +0.177 | 0.058 to 0.306 | 0.0955 | 0.0254 | 0.212 | 0.340 | 0.161 | 0.238 |
| `claude-sonnet-5` | 0.979 | 0.928 | +0.051 | -0.009 to 0.122 | 0.8717 | 0.5798 | 0.095 | 0.186 | 0.082 | 0.126 |

**The cliff under two variations.** Unusable rows left out of the template means, because an empty row at the output cap measures finishing within the ceiling as well as solving, and most such rows fall on Advanced templates; and without the nine symbolic templates (D-138). Welch p, Holm across models.

| model | gap, as scored | gap, unusable left out | 95% CI | Welch p (Holm) | gap, no symbolic | 95% CI | Welch p (Holm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | +0.099 | +0.087 | -0.021 to 0.200 | 1.0000 | +0.128 | 0.010 to 0.251 | 0.3184 |
| `gemma-4-26b-a4b` | +0.127 | +0.127 | 0.025 to 0.233 | 0.2334 | +0.152 | 0.045 to 0.266 | 0.0964 |
| `deepseek-v4.1-flash` | +0.013 | +0.001 | -0.062 to 0.075 | 1.0000 | +0.019 | -0.048 to 0.104 | 1.0000 |
| `qwen3-235b-a22b-2507` | +0.061 | +0.061 | -0.033 to 0.157 | 1.0000 | +0.076 | -0.019 to 0.179 | 0.7953 |
| `glm-5.3-flash` | +0.033 | -0.009 | -0.047 to 0.029 | 1.0000 | +0.042 | -0.010 to 0.102 | 0.7953 |
| `glm-5.3` | +0.124 | -0.001 | -0.027 to 0.021 | 1.0000 | +0.134 | 0.036 to 0.249 | 0.1949 |
| `muse-glimmer-30b` | -0.005 | -0.025 | -0.070 to 0.015 | 1.0000 | +0.009 | -0.041 to 0.062 | 1.0000 |
| `kimi-k3` | +0.026 | +0.016 | -0.039 to 0.076 | 1.0000 | +0.017 | -0.032 to 0.068 | 1.0000 |
| `gpt-5.4-mini` | +0.117 | +0.117 | 0.001 to 0.239 | 0.5680 | +0.146 | 0.029 to 0.272 | 0.2189 |
| `gemini-3.1-flash-lite` | +0.177 | +0.177 | 0.058 to 0.306 | 0.0955 | +0.186 | 0.063 to 0.313 | 0.0781 |
| `claude-sonnet-5` | +0.051 | +0.051 | -0.009 to 0.122 | 1.0000 | +0.052 | -0.001 to 0.126 | 0.7953 |

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 70 items with no milestones (all 15 of 3 templates and some items of 12 more); 52 templates have an item with at most one milestone (46 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The floor is the same trace scored against a sibling item's milestones (D-149): what coverage a trace reaches by chance, per model, on the readable wrong answers. The digit rule's flag counts a trace when any step is flagged. Against the experts, on fully solved traces, it had precision 0.750 and recall 0.320 (SCORER_VALIDATION.md), so its rate is not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | floor, readable | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 395 | 55 | 0.455 | 0.388 to 0.529 | 0.529 | 0.465 to 0.597 |  | pending |
| `gemma-4-26b-a4b` | 306 | 0 | 0.591 | 0.519 to 0.669 | 0.591 | 0.518 to 0.667 |  | pending |
| `deepseek-v4.1-flash` | 108 | 14 | 0.443 | 0.321 to 0.582 | 0.514 | 0.380 to 0.655 |  | pending |
| `qwen3-235b-a22b-2507` | 237 | 1 | 0.581 | 0.511 to 0.656 | 0.584 | 0.515 to 0.658 |  | pending |
| `glm-5.3-flash` | 94 | 42 | 0.279 | 0.171 to 0.433 | 0.540 | 0.371 to 0.745 |  | pending |
| `glm-5.3` | 119 | 100 | 0.084 | 0.032 to 0.192 | 0.785 | 0.631 to 0.950 |  | pending |
| `muse-glimmer-30b` | 80 | 29 | 0.440 | 0.260 to 0.699 | 0.737 | 0.608 to 0.886 |  | pending |
| `kimi-k3` | 108 | 5 | 0.674 | 0.547 to 0.784 | 0.710 | 0.601 to 0.812 |  | pending |
| `gpt-5.4-mini` | 384 | 0 | 0.475 | 0.403 to 0.551 | 0.475 | 0.404 to 0.553 |  | pending |
| `gemini-3.1-flash-lite` | 305 | 0 | 0.495 | 0.438 to 0.555 | 0.495 | 0.435 to 0.557 |  | pending |
| `claude-sonnet-5` | 61 | 0 | 0.691 | 0.597 to 0.857 | 0.691 | 0.598 to 0.854 |  | pending |

**The digit rule.** The wrong-answer rate is over the answered wrong answers (D-149: an empty trace has no step to flag). Beside it, the 1% tolerance E4 shipped with, blind to most slips (RESULTS_X1), and where the first flag falls in a fully solved trace (0 = first step, 1 = last).

| model | wrong answers, answered | flag rate | 95% CI | fully solved traces | flag rate | 95% CI | at 1% | claims checked per trace | traces with a claim | first flag: traces, median position |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 342 | 0.193 | 0.124 to 0.278 | 1776 | 0.153 | 0.118 to 0.192 | 0.082 | 2.25 | 0.604 | 272, 0.600 |
| `gemma-4-26b-a4b` | 306 | 0.389 | 0.278 to 0.504 | 1869 | 0.135 | 0.102 to 0.173 | 0.035 | 4.12 | 0.673 | 253, 0.571 |
| `deepseek-v4.1-flash` | 94 | 0.000 | 0.000 to 0.000 | 2091 | 0.005 | 0.001 to 0.010 | 0.008 | 1.82 | 0.527 | 11, 0.500 |
| `qwen3-235b-a22b-2507` | 237 | 0.380 | 0.247 to 0.511 | 1914 | 0.216 | 0.174 to 0.261 | 0.096 | 8.35 | 0.855 | 414, 0.577 |
| `glm-5.3-flash` | 52 | 0.096 | 0.029 to 0.200 | 2138 | 0.063 | 0.046 to 0.082 | 0.061 | 4.10 | 0.788 | 134, 0.840 |
| `glm-5.3` | 19 | 0.105 | 0.000 to 0.207 | 2109 | 0.066 | 0.051 to 0.083 | 0.070 | 4.44 | 0.843 | 139, 0.778 |
| `muse-glimmer-30b` | 51 | 0.196 | 0.108 to 0.326 | 2146 | 0.122 | 0.095 to 0.151 | 0.103 | 3.10 | 0.679 | 262, 0.414 |
| `kimi-k3` | 103 | 0.000 | 0.000 to 0.000 | 2099 | 0.018 | 0.012 to 0.024 | 0.019 | 1.68 | 0.471 | 37, 0.429 |
| `gpt-5.4-mini` | 384 | 0.172 | 0.087 to 0.262 | 1783 | 0.096 | 0.069 to 0.127 | 0.016 | 2.56 | 0.611 | 172, 0.652 |
| `gemini-3.1-flash-lite` | 305 | 0.308 | 0.187 to 0.434 | 1901 | 0.079 | 0.056 to 0.107 | 0.034 | 3.65 | 0.740 | 151, 0.500 |
| `claude-sonnet-5` | 61 | 0.164 | 0.031 to 0.291 | 2163 | 0.116 | 0.089 to 0.145 | 0.037 | 4.09 | 0.785 | 251, 0.615 |

**Wrong-answer rate against the item's milestone count** (D-149): how failure grows with the depth of the gold derivation.

| model | 0 (70 items) | 1 (290 items) | 2 (452 items) | 3 (367 items) | 4-5 (532 items) | 6+ (539 items) |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.214 | 0.066 | 0.126 | 0.180 | 0.195 | 0.249 |
| `gemma-4-26b-a4b` | 0.100 | 0.045 | 0.100 | 0.093 | 0.086 | 0.299 |
| `deepseek-v4.1-flash` | 0.100 | 0.052 | 0.024 | 0.011 | 0.038 | 0.095 |
| `qwen3-235b-a22b-2507` | 0.157 | 0.034 | 0.084 | 0.084 | 0.079 | 0.195 |
| `glm-5.3-flash` | 0.100 | 0.021 | 0.024 | 0.014 | 0.038 | 0.083 |
| `glm-5.3` | 0.100 | 0.003 | 0.015 | 0.027 | 0.051 | 0.124 |
| `muse-glimmer-30b` | 0.114 | 0.003 | 0.027 | 0.030 | 0.028 | 0.061 |
| `kimi-k3` | 0.114 | 0.045 | 0.024 | 0.014 | 0.039 | 0.093 |
| `gpt-5.4-mini` | 0.114 | 0.100 | 0.111 | 0.106 | 0.173 | 0.308 |
| `gemini-3.1-flash-lite` | 0.114 | 0.024 | 0.080 | 0.074 | 0.141 | 0.282 |
| `claude-sonnet-5` | 0.086 | 0.000 | 0.027 | 0.011 | 0.021 | 0.052 |

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

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one. The three added readings (D-147, D-149) the experts could not arbitrate: the half-unit window requires a correct rounding at the precision shown where the rule accepts one unit either way; the whole-trace reading credits a quantity the question asks for when it is stated in the body and left off the Answer line (the prompt asks for it there); the pool without the nine symbolic templates, whose answers the check scores by the numbers they state (D-138). "Unusable excluded" is an item mean over the usable rows, not a mean of template means. Tau's noise floor: the ordering on one random half of each template's items against the other, median 0.881 over 200 splits (quartiles 0.855 to 0.891); the top five models lie within 0.012 of each other, so tau falls below 1 from sampling noise alone.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates | without the 9 symbolic templates | half-unit window | whole trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.800 | 0.807 | 0.815 | 0.789 | 0.827 | 0.805 | 0.804 |  |  |
| `gemma-4-26b-a4b` | 0.826 | 0.847 | 0.878 | 0.831 | 0.847 | 0.844 | 0.848 |  |  |
| `deepseek-v4.1-flash` | 0.937 | 0.941 | 0.942 | 0.929 | 0.947 | 0.939 | 0.943 |  |  |
| `qwen3-235b-a22b-2507` | 0.859 | 0.873 | 0.885 | 0.851 | 0.873 | 0.870 | 0.875 |  |  |
| `glm-5.3-flash` | 0.950 | 0.954 | 0.956 | 0.950 | 0.972 | 0.953 | 0.957 |  |  |
| `glm-5.3` | 0.939 | 0.942 | 0.944 | 0.937 | 0.986 | 0.941 | 0.946 |  |  |
| `muse-glimmer-30b` | 0.957 | 0.959 | 0.962 | 0.954 | 0.972 | 0.958 | 0.961 |  |  |
| `kimi-k3` | 0.939 | 0.942 | 0.945 | 0.933 | 0.945 | 0.942 | 0.949 |  |  |
| `gpt-5.4-mini` | 0.799 | 0.811 | 0.822 | 0.792 | 0.811 | 0.807 | 0.810 |  |  |
| `gemini-3.1-flash-lite` | 0.845 | 0.855 | 0.868 | 0.845 | 0.855 | 0.852 | 0.863 |  |  |
| `claude-sonnet-5` | 0.960 | 0.967 | 0.973 | 0.961 | 0.967 | 0.966 | 0.975 |  |  |
| tau with the headline | 0.954 | 1 | 0.964 | 0.964 | 0.673 | 1.000 | 1.000 |  |  |

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

### By serving endpoint

Open-weight models were served by several endpoints under the fp8-or-better rule (D-133). The harness dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is confounded with the templates each endpoint happened to serve; the matched difference compares an endpoint with the other endpoints on the same templates (D-149).

| model | endpoint | rows | raw score | unusable | matched difference | templates matched |
|---|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | Darkbloom | 2242 | 0.807 | 0.025 | -0.179 | 8 |
| `gpt-oss-20b` | DekaLLM | 7 | 0.857 | 0.000 | +0.184 | 7 |
| `gpt-oss-20b` | DeepInfra | 1 | 1.000 | 0.000 | +0.143 | 1 |
| `gemma-4-26b-a4b` | DekaLLM | 1263 | 0.850 | 0.000 | +0.006 | 150 |
| `gemma-4-26b-a4b` | NextBit | 987 | 0.843 | 0.000 | -0.006 | 150 |
| `deepseek-v4.1-flash` | CoreWeave | 793 | 0.929 | 0.006 | -0.006 | 143 |
| `deepseek-v4.1-flash` | Parasail | 472 | 0.947 | 0.000 | +0.002 | 140 |
| `deepseek-v4.1-flash` | Morph | 458 | 0.968 | 0.000 | +0.007 | 135 |
| `deepseek-v4.1-flash` | DeepInfra | 419 | 0.934 | 0.019 | -0.001 | 79 |
| `deepseek-v4.1-flash` | Makora | 89 | 0.899 | 0.000 | -0.004 | 44 |
| `deepseek-v4.1-flash` | Novita | 17 | 0.912 | 0.059 | -0.048 | 8 |
| `deepseek-v4.1-flash` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `qwen3-235b-a22b-2507` | GMICloud | 2200 | 0.874 | 0.000 | -0.016 | 36 |
| `qwen3-235b-a22b-2507` | Parasail | 50 | 0.820 | 0.020 | +0.016 | 36 |
| `glm-5.3-flash` | GMICloud | 990 | 0.951 | 0.017 | -0.004 | 146 |
| `glm-5.3-flash` | Novita | 398 | 0.975 | 0.010 | +0.002 | 120 |
| `glm-5.3-flash` | Parasail | 224 | 0.929 | 0.045 | -0.006 | 76 |
| `glm-5.3-flash` | Sail Research | 205 | 0.946 | 0.039 | -0.014 | 116 |
| `glm-5.3-flash` | Phala | 167 | 0.973 | 0.006 | +0.012 | 75 |
| `glm-5.3-flash` | Z.AI | 84 | 0.952 | 0.012 | +0.005 | 49 |
| `glm-5.3-flash` | StreamLake | 77 | 0.974 | 0.000 | +0.000 | 34 |
| `glm-5.3-flash` | NextBit | 43 | 0.849 | 0.023 | -0.050 | 13 |
| `glm-5.3-flash` | AtlasCloud | 42 | 1.000 | 0.000 | +0.060 | 15 |
| `glm-5.3-flash` | SiliconFlow | 19 | 1.000 | 0.000 | +0.119 | 6 |
| `glm-5.3-flash` | Morph | 1 | 1.000 | 0.000 | +0.000 | 1 |
| `glm-5.3` | Baidu | 918 | 0.953 | 0.036 | -0.012 | 86 |
| `glm-5.3` | Sail Research | 748 | 0.944 | 0.037 | +0.010 | 77 |
| `glm-5.3` | Morph | 245 | 0.961 | 0.020 | +0.007 | 122 |
| `glm-5.3` | Novita | 172 | 0.910 | 0.081 | -0.007 | 42 |
| `glm-5.3` | AtlasCloud | 147 | 0.864 | 0.136 | -0.006 | 55 |
| `glm-5.3` | BaseTen | 15 | 1.000 | 0.000 | +0.000 | 5 |
| `glm-5.3` | AkashML | 3 | 1.000 | 0.000 | +0.024 | 3 |
| `glm-5.3` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `kimi-k3` | Morph | 2231 | 0.943 | 0.002 | -0.023 | 7 |
| `kimi-k3` | BaseTen | 19 | 0.921 | 0.000 | +0.023 | 7 |
| `gpt-5.4-mini` | Azure | 2249 | 0.811 | 0.000 | -0.107 | 1 |
| `gpt-5.4-mini` | OpenAI | 1 | 1.000 | 0.000 | +0.107 | 1 |

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

## Provenance

`analyze.py` at commit `bf4a43b`; the score store `main_pre_d137` scored at commit `9915400` on 2026-09-29T10:23:42+00:00. The evaluator hashes and the per-model trace hashes are in `results.json`.
