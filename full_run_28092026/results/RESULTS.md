# Results of the full run: the deterministic stack

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring, and the corrections and additions made after the first results were read are labelled where they appear (D-146 to D-149). The answer check, E3 and E4's digit rule; E5's judge has not run, so its columns are empty. Every interval is 95% and resamples templates (B = 10,000): all 150 for a model's score, and the templates with a qualifying trace for a rate over a subset of traces. Every resampling test draws 100,000 permutations, so its Holm-adjusted floor over 55 pairs is 0.0006.

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0. SD within: the mean over templates of the SD of the score across a template's 15 items; SD between: the SD of the 150 template means (the per-template variability the July rebuttal promised). Unusable is empty plus unreadable; capped: answered rows that stopped at the output cap and are scored on what they state (D-117).

| model | answer score | 95% CI | SD within | SD between | fully solved | 95% CI | correct | partial | incorrect | unusable (empty + unreadable) | capped, scored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.970 | 0.949 to 0.986 | 0.041 | 0.117 | 0.963 | 0.940 to 0.982 | 2167 | 29 | 40 | 14 (14 + 0) | 1 |
| `claude-sonnet-5` | 0.965 | 0.942 to 0.983 | 0.047 | 0.128 | 0.959 | 0.935 to 0.979 | 2157 | 28 | 65 | 0 (0 + 0) | 0 |
| `kimi-k3` | 0.964 | 0.942 to 0.981 | 0.056 | 0.122 | 0.958 | 0.936 to 0.977 | 2156 | 26 | 63 | 5 (5 + 0) | 4 |
| `glm-5.3-flash` | 0.964 | 0.941 to 0.982 | 0.049 | 0.131 | 0.960 | 0.935 to 0.980 | 2159 | 18 | 31 | 42 (42 + 0) | 1 |
| `muse-glimmer-30b` | 0.960 | 0.937 to 0.979 | 0.057 | 0.131 | 0.956 | 0.932 to 0.976 | 2150 | 19 | 52 | 29 (29 + 0) | 0 |
| `glm-5.3` | 0.942 | 0.912 to 0.969 | 0.061 | 0.180 | 0.937 | 0.906 to 0.965 | 2109 | 23 | 18 | 100 (100 + 0) | 6 |
| `qwen3-235b-a22b-2507` | 0.873 | 0.836 to 0.906 | 0.152 | 0.222 | 0.847 | 0.804 to 0.889 | 1906 | 115 | 228 | 1 (0 + 1) | 1 |
| `gemini-3.1-flash-lite` | 0.858 | 0.812 to 0.900 | 0.123 | 0.274 | 0.849 | 0.803 to 0.892 | 1910 | 41 | 299 | 0 (0 + 0) | 0 |
| `gemma-4-26b-a4b` | 0.851 | 0.809 to 0.888 | 0.162 | 0.247 | 0.834 | 0.793 to 0.874 | 1877 | 75 | 298 | 0 (0 + 0) | 1 |
| `gpt-5.4-mini` | 0.837 | 0.794 to 0.877 | 0.172 | 0.261 | 0.824 | 0.780 to 0.865 | 1853 | 61 | 336 | 0 (0 + 0) | 0 |
| `gpt-oss-20b` | 0.805 | 0.760 to 0.847 | 0.213 | 0.276 | 0.790 | 0.744 to 0.833 | 1778 | 65 | 352 | 55 (53 + 2) | 4 |

Of the 55 pairs, 31 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The same template-level test on the fully-solved rate gives the same verdict on 55 of 55. McNemar's exact test on the paired item verdicts (Holm) gives the same verdict on 46 of 55; on 9 of the 9 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level tests. 7 answered rows across the roster ended with a finish reason other than stop or length (a provider fault inside a 200) and are scored on what they state; the harness now retries such a reply (D-148).

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | template p (Holm) | McNemar p (Holm) | same verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.165 | -0.208 to -0.124 | 0.0005 | -0.173 | -0.217 to -0.131 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.160 | -0.202 to -0.121 | 0.0005 | -0.168 | -0.210 to -0.128 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `kimi-k3` | -0.159 | -0.203 to -0.118 | 0.0005 | -0.168 | -0.212 to -0.125 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.159 | -0.201 to -0.119 | 0.0005 | -0.169 | -0.212 to -0.129 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.155 | -0.194 to -0.118 | 0.0005 | -0.165 | -0.206 to -0.127 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3` | -0.138 | -0.182 to -0.096 | 0.0005 | -0.147 | -0.191 to -0.105 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.132 | 0.094 to 0.173 | 0.0005 | +0.140 | 0.100 to 0.181 | 0.0005 | 0.0000 | yes |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.128 | -0.167 to -0.091 | 0.0005 | -0.135 | -0.175 to -0.098 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gpt-5.4-mini` | +0.127 | 0.089 to 0.166 | 0.0005 | +0.135 | 0.096 to 0.175 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.126 | 0.093 to 0.164 | 0.0005 | +0.136 | 0.101 to 0.174 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.123 | 0.090 to 0.158 | 0.0005 | +0.132 | 0.098 to 0.168 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.119 | -0.155 to -0.086 | 0.0005 | -0.129 | -0.167 to -0.095 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.114 | -0.149 to -0.082 | 0.0005 | -0.124 | -0.160 to -0.092 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.113 | -0.149 to -0.080 | 0.0005 | -0.124 | -0.160 to -0.090 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.113 | -0.146 to -0.082 | 0.0005 | -0.125 | -0.161 to -0.093 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.112 | 0.074 to 0.151 | 0.0005 | +0.114 | 0.078 to 0.154 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.109 | -0.144 to -0.077 | 0.0005 | -0.121 | -0.160 to -0.086 | 0.0005 | 0.0000 | yes |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.107 | -0.146 to -0.072 | 0.0005 | -0.110 | -0.149 to -0.074 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.106 | 0.072 to 0.143 | 0.0005 | +0.109 | 0.073 to 0.148 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.106 | 0.073 to 0.141 | 0.0005 | +0.111 | 0.077 to 0.148 | 0.0005 | 0.0000 | yes |
| `glm-5.3` | `gpt-5.4-mini` | +0.105 | 0.074 to 0.140 | 0.0005 | +0.114 | 0.080 to 0.150 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.102 | 0.068 to 0.139 | 0.0005 | +0.107 | 0.071 to 0.145 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.097 | 0.069 to 0.128 | 0.0005 | +0.116 | 0.081 to 0.155 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.092 | -0.122 to -0.065 | 0.0005 | -0.112 | -0.149 to -0.076 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.092 | -0.126 to -0.061 | 0.0005 | -0.103 | -0.138 to -0.072 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.091 | -0.121 to -0.064 | 0.0005 | -0.111 | -0.148 to -0.077 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.091 | -0.120 to -0.064 | 0.0005 | -0.112 | -0.149 to -0.079 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.087 | -0.116 to -0.060 | 0.0005 | -0.108 | -0.146 to -0.074 | 0.0005 | 0.0000 | yes |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.084 | 0.053 to 0.118 | 0.0005 | +0.088 | 0.057 to 0.123 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.070 | -0.101 to -0.041 | 0.0005 | -0.090 | -0.129 to -0.055 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `glm-5.3` | +0.021 | 0.008 to 0.037 | 0.0285 | +0.022 | 0.008 to 0.038 | 0.0345 | 0.0000 | yes |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.068 | -0.111 to -0.026 | 0.0530 | -0.057 | -0.105 to -0.008 | 0.5100 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.027 | 0.009 to 0.049 | 0.1575 | +0.026 | 0.005 to 0.049 | 0.4425 | 0.0000 | yes |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.053 | -0.099 to -0.007 | 0.5509 | -0.059 | -0.105 to -0.012 | 0.3456 | 0.0000 | yes |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.046 | -0.086 to -0.006 | 0.5509 | -0.044 | -0.086 to -0.004 | 0.8047 | 0.0001 | yes |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.036 | 0.005 to 0.068 | 0.5840 | +0.024 | -0.016 to 0.061 | 1.0000 | 0.1320 | yes |
| `glm-5.3` | `claude-sonnet-5` | -0.022 | -0.047 to -0.003 | 0.9048 | -0.021 | -0.046 to -0.002 | 1.0000 | 0.0002 | yes |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.032 | -0.076 to 0.011 | 1.0000 | -0.033 | -0.079 to 0.009 | 1.0000 | 0.0136 | yes |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.022 | -0.052 to 0.008 | 1.0000 | -0.013 | -0.048 to 0.025 | 1.0000 | 1.0000 | yes |
| `glm-5.3` | `kimi-k3` | -0.022 | -0.047 to 0.001 | 1.0000 | -0.021 | -0.047 to 0.002 | 1.0000 | 0.0016 | yes |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.021 | -0.048 to 0.005 | 1.0000 | -0.025 | -0.053 to 0.001 | 1.0000 | 0.0136 | yes |
| `glm-5.3` | `muse-glimmer-30b` | -0.017 | -0.042 to 0.005 | 1.0000 | -0.018 | -0.044 to 0.004 | 1.0000 | 0.0105 | yes |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.015 | -0.014 to 0.045 | 1.0000 | -0.002 | -0.040 to 0.036 | 1.0000 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.014 | -0.015 to 0.043 | 1.0000 | +0.011 | -0.020 to 0.041 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.010 | -0.010 to 0.029 | 1.0000 | +0.008 | -0.013 to 0.027 | 1.0000 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.007 | -0.034 to 0.019 | 1.0000 | -0.015 | -0.041 to 0.012 | 1.0000 | 0.8220 | yes |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.006 | -0.005 to 0.018 | 1.0000 | +0.004 | -0.010 to 0.017 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.006 | -0.004 to 0.017 | 1.0000 | +0.005 | -0.006 to 0.017 | 1.0000 | 1.0000 | yes |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.005 | -0.025 to 0.015 | 1.0000 | -0.003 | -0.025 to 0.019 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.005 | -0.007 to 0.017 | 1.0000 | +0.004 | -0.009 to 0.019 | 1.0000 | 1.0000 | yes |
| `muse-glimmer-30b` | `kimi-k3` | -0.004 | -0.026 to 0.017 | 1.0000 | -0.003 | -0.024 to 0.019 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.004 | -0.015 to 0.022 | 1.0000 | +0.004 | -0.016 to 0.022 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.001 | -0.018 to 0.013 | 1.0000 | +0.001 | -0.016 to 0.016 | 1.0000 | 1.0000 | yes |
| `kimi-k3` | `claude-sonnet-5` | -0.001 | -0.015 to 0.012 | 1.0000 | -0.000 | -0.016 to 0.013 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `kimi-k3` | -0.000 | -0.016 to 0.016 | 1.0000 | +0.001 | -0.015 to 0.017 | 1.0000 | 1.0000 | yes |

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier. **Corrected after the first results (D-146):** the planned test permuted the tier labels on the raw difference, which is liberal when the smaller tier has the larger spread, and Advanced template means spread two to four times as widely as Easy ones (the two SD columns); its count of models also moved with the seed. The count now rests on Welch's t-test, Holm across the eleven models: **3 of 11** hold, against 6 under the planned permutation, printed beside it. The detectable gap at 80% power is given as the plan defined it (2.80 x pooled sigma x sqrt(1/58 + 1/34)) and at the strictest Holm step from the Welch standard error.

| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | planned p (Holm) | SD Easy | SD Adv | detectable, planned | detectable, Holm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.891 | 0.706 | +0.185 | 0.070 to 0.304 | 0.0353 | 0.0115 | 0.227 | 0.305 | 0.156 | 0.221 |
| `gemma-4-26b-a4b` | 0.907 | 0.760 | +0.148 | 0.049 to 0.254 | 0.0593 | 0.0203 | 0.199 | 0.270 | 0.137 | 0.195 |
| `deepseek-v4.1-flash` | 0.988 | 0.917 | +0.071 | 0.009 to 0.153 | 0.2999 | 0.0546 | 0.063 | 0.209 | 0.082 | 0.135 |
| `qwen3-235b-a22b-2507` | 0.906 | 0.826 | +0.079 | -0.016 to 0.179 | 0.3240 | 0.1871 | 0.184 | 0.259 | 0.130 | 0.186 |
| `glm-5.3-flash` | 0.986 | 0.914 | +0.072 | 0.014 to 0.146 | 0.2546 | 0.0381 | 0.072 | 0.192 | 0.078 | 0.126 |
| `glm-5.3` | 0.978 | 0.844 | +0.133 | 0.040 to 0.240 | 0.1076 | 0.0115 | 0.082 | 0.299 | 0.116 | 0.193 |
| `muse-glimmer-30b` | 0.970 | 0.920 | +0.050 | -0.009 to 0.113 | 0.3240 | 0.1871 | 0.127 | 0.155 | 0.084 | 0.115 |
| `kimi-k3` | 0.983 | 0.920 | +0.064 | 0.004 to 0.138 | 0.3240 | 0.0924 | 0.076 | 0.199 | 0.082 | 0.131 |
| `gpt-5.4-mini` | 0.901 | 0.721 | +0.181 | 0.065 to 0.305 | 0.0448 | 0.0115 | 0.197 | 0.324 | 0.152 | 0.226 |
| `gemini-3.1-flash-lite` | 0.921 | 0.726 | +0.194 | 0.075 to 0.325 | 0.0433 | 0.0090 | 0.199 | 0.345 | 0.159 | 0.238 |
| `claude-sonnet-5` | 0.980 | 0.914 | +0.066 | -0.004 to 0.149 | 0.3240 | 0.1078 | 0.102 | 0.217 | 0.093 | 0.146 |

**The cliff under two variations.** Unusable rows left out of the template means, because an empty row at the output cap measures finishing within the ceiling as well as solving, and most such rows fall on Advanced templates; and without the nine symbolic templates (D-138). Welch p, Holm across models.

| model | gap, as scored | gap, unusable left out | 95% CI | Welch p (Holm) | gap, no symbolic | 95% CI | Welch p (Holm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | +0.185 | +0.175 | 0.063 to 0.290 | 0.0464 | +0.201 | 0.083 to 0.324 | 0.0270 |
| `gemma-4-26b-a4b` | +0.148 | +0.148 | 0.049 to 0.254 | 0.0593 | +0.161 | 0.056 to 0.271 | 0.0436 |
| `deepseek-v4.1-flash` | +0.071 | +0.059 | 0.000 to 0.138 | 0.6992 | +0.056 | 0.005 to 0.134 | 0.4933 |
| `qwen3-235b-a22b-2507` | +0.079 | +0.079 | -0.016 to 0.179 | 0.6992 | +0.079 | -0.014 to 0.181 | 0.4933 |
| `glm-5.3-flash` | +0.072 | +0.036 | -0.012 to 0.101 | 0.7050 | +0.075 | 0.019 to 0.151 | 0.2426 |
| `glm-5.3` | +0.133 | +0.009 | -0.021 to 0.042 | 0.7050 | +0.137 | 0.037 to 0.252 | 0.1410 |
| `muse-glimmer-30b` | +0.050 | +0.029 | -0.023 to 0.084 | 0.7050 | +0.043 | -0.006 to 0.099 | 0.4933 |
| `kimi-k3` | +0.064 | +0.054 | -0.004 to 0.129 | 0.6992 | +0.046 | 0.003 to 0.105 | 0.4933 |
| `gpt-5.4-mini` | +0.181 | +0.181 | 0.065 to 0.305 | 0.0464 | +0.203 | 0.085 to 0.330 | 0.0313 |
| `gemini-3.1-flash-lite` | +0.194 | +0.194 | 0.075 to 0.325 | 0.0464 | +0.196 | 0.074 to 0.323 | 0.0399 |
| `claude-sonnet-5` | +0.066 | +0.066 | -0.004 to 0.149 | 0.6992 | +0.058 | 0.002 to 0.132 | 0.4933 |

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 70 items with no milestones (all 15 of 3 templates and some items of 12 more); 52 templates have an item with at most one milestone (46 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The floor is the same trace scored against a sibling item's milestones (D-149): what coverage a trace reaches by chance, per model, on the readable wrong answers. The digit rule's flag counts a trace when any step is flagged. Against the experts, on fully solved traces, it had precision 0.750 and recall 0.320 (SCORER_VALIDATION.md), so its rate is not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | floor, readable | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 407 | 55 | 0.394 | 0.329 to 0.464 | 0.456 | 0.389 to 0.527 | 0.113 | pending |
| `gemma-4-26b-a4b` | 298 | 0 | 0.588 | 0.517 to 0.666 | 0.588 | 0.515 to 0.667 | 0.142 | pending |
| `deepseek-v4.1-flash` | 54 | 14 | 0.385 | 0.241 to 0.570 | 0.549 | 0.379 to 0.708 | 0.107 | pending |
| `qwen3-235b-a22b-2507` | 229 | 1 | 0.586 | 0.511 to 0.664 | 0.589 | 0.516 to 0.667 | 0.196 | pending |
| `glm-5.3-flash` | 73 | 42 | 0.225 | 0.094 to 0.393 | 0.619 | 0.441 to 0.819 | 0.168 | pending |
| `glm-5.3` | 118 | 100 | 0.069 | 0.021 to 0.158 | 0.697 | 0.535 to 0.938 | 0.117 | pending |
| `muse-glimmer-30b` | 81 | 29 | 0.370 | 0.224 to 0.564 | 0.613 | 0.536 to 0.731 | 0.154 | pending |
| `kimi-k3` | 68 | 5 | 0.567 | 0.419 to 0.690 | 0.618 | 0.507 to 0.727 | 0.135 | pending |
| `gpt-5.4-mini` | 336 | 0 | 0.422 | 0.352 to 0.498 | 0.422 | 0.351 to 0.500 | 0.114 | pending |
| `gemini-3.1-flash-lite` | 299 | 0 | 0.511 | 0.445 to 0.581 | 0.511 | 0.445 to 0.579 | 0.107 | pending |
| `claude-sonnet-5` | 65 | 0 | 0.660 | 0.569 to 0.789 | 0.660 | 0.569 to 0.791 | 0.158 | pending |

**The digit rule.** The wrong-answer rate is over the answered wrong answers (D-149: an empty trace has no step to flag). Beside it, the 1% tolerance E4 shipped with, blind to most slips (RESULTS_X1), and where the first flag falls in a fully solved trace (0 = first step, 1 = last).

| model | wrong answers, answered | flag rate | 95% CI | fully solved traces | flag rate | 95% CI | at 1% | claims checked per trace | traces with a claim | first flag: traces, median position |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 354 | 0.206 | 0.142 to 0.278 | 1778 | 0.151 | 0.115 to 0.190 | 0.084 | 2.25 | 0.604 | 269, 0.600 |
| `gemma-4-26b-a4b` | 298 | 0.409 | 0.295 to 0.524 | 1877 | 0.133 | 0.100 to 0.170 | 0.034 | 4.12 | 0.673 | 249, 0.571 |
| `deepseek-v4.1-flash` | 40 | 0.025 | 0.000 to 0.071 | 2167 | 0.005 | 0.001 to 0.010 | 0.007 | 1.82 | 0.527 | 10, 0.500 |
| `qwen3-235b-a22b-2507` | 229 | 0.393 | 0.259 to 0.524 | 1906 | 0.217 | 0.174 to 0.262 | 0.096 | 8.35 | 0.855 | 413, 0.571 |
| `glm-5.3-flash` | 31 | 0.194 | 0.104 to 0.385 | 2159 | 0.062 | 0.045 to 0.081 | 0.060 | 4.10 | 0.788 | 133, 0.833 |
| `glm-5.3` | 18 | 0.111 | 0.000 to 0.222 | 2109 | 0.065 | 0.050 to 0.082 | 0.069 | 4.44 | 0.843 | 137, 0.778 |
| `muse-glimmer-30b` | 52 | 0.173 | 0.089 to 0.294 | 2150 | 0.123 | 0.096 to 0.152 | 0.103 | 3.10 | 0.679 | 265, 0.400 |
| `kimi-k3` | 63 | 0.000 | 0.000 to 0.000 | 2156 | 0.017 | 0.012 to 0.023 | 0.018 | 1.68 | 0.471 | 37, 0.429 |
| `gpt-5.4-mini` | 336 | 0.176 | 0.085 to 0.277 | 1853 | 0.097 | 0.070 to 0.126 | 0.016 | 2.56 | 0.611 | 179, 0.600 |
| `gemini-3.1-flash-lite` | 299 | 0.308 | 0.185 to 0.436 | 1910 | 0.080 | 0.056 to 0.107 | 0.034 | 3.65 | 0.740 | 152, 0.500 |
| `claude-sonnet-5` | 65 | 0.154 | 0.034 to 0.283 | 2157 | 0.116 | 0.089 to 0.145 | 0.037 | 4.09 | 0.785 | 250, 0.615 |

**Wrong-answer rate against the item's milestone count** (D-149): how failure grows with the depth of the gold derivation.

| model | 0 (70 items) | 1 (290 items) | 2 (452 items) | 3 (367 items) | 4-5 (532 items) | 6+ (539 items) |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.229 | 0.034 | 0.108 | 0.191 | 0.224 | 0.265 |
| `gemma-4-26b-a4b` | 0.100 | 0.045 | 0.104 | 0.104 | 0.085 | 0.275 |
| `deepseek-v4.1-flash` | 0.100 | 0.000 | 0.011 | 0.014 | 0.015 | 0.054 |
| `qwen3-235b-a22b-2507` | 0.157 | 0.041 | 0.075 | 0.084 | 0.081 | 0.182 |
| `glm-5.3-flash` | 0.100 | 0.003 | 0.024 | 0.008 | 0.019 | 0.076 |
| `glm-5.3` | 0.100 | 0.003 | 0.020 | 0.025 | 0.049 | 0.122 |
| `muse-glimmer-30b` | 0.114 | 0.000 | 0.033 | 0.033 | 0.030 | 0.056 |
| `kimi-k3` | 0.114 | 0.000 | 0.018 | 0.019 | 0.028 | 0.056 |
| `gpt-5.4-mini` | 0.100 | 0.028 | 0.104 | 0.093 | 0.156 | 0.291 |
| `gemini-3.1-flash-lite` | 0.114 | 0.021 | 0.077 | 0.079 | 0.143 | 0.269 |
| `claude-sonnet-5` | 0.100 | 0.000 | 0.035 | 0.014 | 0.023 | 0.046 |

## Q4. Consistency within a template

The share of templates fully solved on all 15 instances, on some, and on none: the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) and the 92 others. Intervals are in `results.json`.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.362 | 0.569 | 0.069 | 0.326 | 0.663 | 0.011 |
| `gemma-4-26b-a4b` | 0.466 | 0.500 | 0.034 | 0.500 | 0.478 | 0.022 |
| `deepseek-v4.1-flash` | 0.897 | 0.103 | 0.000 | 0.826 | 0.174 | 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.534 | 0.017 | 0.533 | 0.435 | 0.033 |
| `glm-5.3-flash` | 0.914 | 0.086 | 0.000 | 0.793 | 0.207 | 0.000 |
| `glm-5.3` | 0.845 | 0.138 | 0.017 | 0.750 | 0.239 | 0.011 |
| `muse-glimmer-30b` | 0.862 | 0.138 | 0.000 | 0.793 | 0.207 | 0.000 |
| `kimi-k3` | 0.862 | 0.138 | 0.000 | 0.761 | 0.239 | 0.000 |
| `gpt-5.4-mini` | 0.483 | 0.483 | 0.034 | 0.446 | 0.522 | 0.033 |
| `gemini-3.1-flash-lite` | 0.672 | 0.276 | 0.052 | 0.576 | 0.402 | 0.022 |
| `claude-sonnet-5` | 0.897 | 0.103 | 0.000 | 0.804 | 0.196 | 0.000 |

## Q5. Paraphrase robustness

Not run yet.

## Sensitivity

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one. The three added readings (D-147, D-149) the experts could not arbitrate: the half-unit window requires a correct rounding at the precision shown where the rule accepts one unit either way; the whole-trace reading credits a quantity the question asks for when it is stated in the body and left off the Answer line (the prompt asks for it there); the pool without the nine symbolic templates, whose answers the check scores by the numbers they state (D-138). "Unusable excluded" is an item mean over the usable rows, not a mean of template means. Tau's noise floor: the ordering on one random half of each template's items against the other, median 0.881 over 200 splits (quartiles 0.844 to 0.891); the top five models lie within 0.012 of each other, so tau falls below 1 from sampling noise alone.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates | without the 9 symbolic templates | half-unit window | whole trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.797 | 0.805 | 0.812 | 0.790 | 0.825 | 0.803 | 0.806 | 0.774 | 0.813 |
| `gemma-4-26b-a4b` | 0.831 | 0.851 | 0.882 | 0.834 | 0.851 | 0.847 | 0.856 | 0.836 | 0.862 |
| `deepseek-v4.1-flash` | 0.965 | 0.970 | 0.971 | 0.963 | 0.976 | 0.970 | 0.979 | 0.966 | 0.972 |
| `qwen3-235b-a22b-2507` | 0.862 | 0.873 | 0.884 | 0.847 | 0.873 | 0.870 | 0.880 | 0.864 | 0.893 |
| `glm-5.3-flash` | 0.962 | 0.964 | 0.965 | 0.960 | 0.982 | 0.963 | 0.968 | 0.960 | 0.966 |
| `glm-5.3` | 0.941 | 0.942 | 0.943 | 0.937 | 0.986 | 0.941 | 0.948 | 0.938 | 0.946 |
| `muse-glimmer-30b` | 0.958 | 0.960 | 0.961 | 0.956 | 0.972 | 0.959 | 0.968 | 0.947 | 0.961 |
| `kimi-k3` | 0.962 | 0.964 | 0.965 | 0.958 | 0.966 | 0.964 | 0.975 | 0.960 | 0.968 |
| `gpt-5.4-mini` | 0.826 | 0.837 | 0.848 | 0.824 | 0.837 | 0.833 | 0.841 | 0.803 | 0.843 |
| `gemini-3.1-flash-lite` | 0.849 | 0.858 | 0.871 | 0.849 | 0.858 | 0.855 | 0.868 | 0.852 | 0.865 |
| `claude-sonnet-5` | 0.960 | 0.965 | 0.970 | 0.959 | 0.965 | 0.964 | 0.976 | 0.958 | 0.969 |
| tau with the headline | 0.891 | 1 | 0.964 | 0.891 | 0.600 | 1.000 | 0.964 | 0.917 | 1.000 |

The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).

## Also reported, not tested

### By branch and level

| model | chemical | civil | electrical | industrial | mechanical | Easy | Intermediate | Advanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.809 | 0.673 | 0.872 | 0.820 | 0.849 | 0.891 | 0.776 | 0.706 |
| `gemma-4-26b-a4b` | 0.734 | 0.882 | 0.866 | 0.878 | 0.894 | 0.907 | 0.848 | 0.760 |
| `deepseek-v4.1-flash` | 0.956 | 0.993 | 0.924 | 0.980 | 0.994 | 0.988 | 0.982 | 0.917 |
| `qwen3-235b-a22b-2507` | 0.827 | 0.940 | 0.849 | 0.871 | 0.877 | 0.906 | 0.867 | 0.826 |
| `glm-5.3-flash` | 0.947 | 0.984 | 0.954 | 0.942 | 0.990 | 0.986 | 0.971 | 0.914 |
| `glm-5.3` | 0.898 | 0.973 | 0.926 | 0.922 | 0.993 | 0.978 | 0.965 | 0.844 |
| `muse-glimmer-30b` | 0.963 | 0.969 | 0.938 | 0.940 | 0.989 | 0.970 | 0.974 | 0.920 |
| `kimi-k3` | 0.964 | 0.987 | 0.920 | 0.980 | 0.969 | 0.983 | 0.971 | 0.920 |
| `gpt-5.4-mini` | 0.800 | 0.838 | 0.862 | 0.831 | 0.854 | 0.901 | 0.841 | 0.721 |
| `gemini-3.1-flash-lite` | 0.776 | 0.862 | 0.856 | 0.860 | 0.937 | 0.921 | 0.872 | 0.726 |
| `claude-sonnet-5` | 0.949 | 0.996 | 0.912 | 0.973 | 0.994 | 0.980 | 0.980 | 0.914 |

### By domain

| domain | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| digital_communications | 0.930 | 0.883 | 0.923 | 0.840 | 0.967 | 0.923 | 0.950 | 0.897 | 0.923 | 0.877 | 0.913 |
| electromagnetics_and_waves | 0.910 | 0.887 | 0.990 | 0.883 | 0.977 | 0.947 | 0.997 | 0.980 | 0.830 | 0.883 | 0.953 |
| fluid_mechanics | 0.880 | 0.893 | 0.997 | 0.913 | 0.987 | 0.987 | 0.987 | 0.987 | 0.950 | 0.950 | 0.987 |
| geotechnical_engineering | 0.520 | 0.913 | 1.000 | 0.947 | 0.993 | 0.993 | 0.933 | 1.000 | 0.913 | 0.960 | 0.987 |
| mechanics_of_materials | 0.800 | 0.827 | 0.993 | 0.857 | 0.987 | 1.000 | 0.990 | 0.927 | 0.803 | 0.867 | 1.000 |
| production_and_inventory | 0.767 | 0.873 | 0.940 | 0.847 | 0.887 | 0.887 | 0.887 | 0.967 | 0.780 | 0.853 | 0.987 |
| quality_and_reliability_control | 0.773 | 0.800 | 1.000 | 0.787 | 0.947 | 0.907 | 0.933 | 0.980 | 0.760 | 0.800 | 0.940 |
| reaction_kinetics | 0.967 | 0.913 | 0.993 | 0.980 | 1.000 | 0.993 | 1.000 | 1.000 | 0.913 | 0.913 | 1.000 |
| signals_and_systems | 0.777 | 0.827 | 0.860 | 0.823 | 0.920 | 0.907 | 0.867 | 0.883 | 0.833 | 0.807 | 0.870 |
| stochastic_operations | 0.920 | 0.960 | 1.000 | 0.980 | 0.993 | 0.973 | 1.000 | 0.993 | 0.953 | 0.927 | 0.993 |
| structural_analysis | 0.733 | 0.913 | 1.000 | 0.947 | 1.000 | 1.000 | 0.987 | 0.993 | 0.827 | 0.913 | 1.000 |
| thermodynamics | 0.656 | 0.517 | 0.900 | 0.728 | 0.872 | 0.761 | 0.911 | 0.922 | 0.667 | 0.678 | 0.878 |
| transport_phenomena | 0.842 | 0.838 | 0.992 | 0.783 | 0.992 | 0.983 | 0.996 | 0.983 | 0.858 | 0.750 | 0.992 |
| vibrations_and_acoustics | 0.867 | 0.963 | 0.993 | 0.860 | 0.997 | 0.993 | 0.990 | 0.993 | 0.810 | 0.993 | 0.997 |
| water_resources | 0.767 | 0.820 | 0.980 | 0.927 | 0.960 | 0.927 | 0.987 | 0.967 | 0.773 | 0.713 | 1.000 |

### By answer type

| answer type | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| array | 0.886 | 0.867 | 0.962 | 0.924 | 0.914 | 0.762 | 0.971 | 0.962 | 0.800 | 0.771 | 0.990 |
| classification | 0.680 | 0.953 | 0.960 | 0.887 | 0.993 | 1.000 | 0.960 | 0.960 | 0.973 | 0.973 | 1.000 |
| multipart | 0.876 | 0.893 | 0.980 | 0.892 | 0.978 | 0.950 | 0.973 | 0.980 | 0.892 | 0.908 | 0.979 |
| scalar | 0.768 | 0.825 | 0.979 | 0.876 | 0.966 | 0.956 | 0.964 | 0.974 | 0.816 | 0.849 | 0.972 |
| symbolic | 0.789 | 0.778 | 0.822 | 0.756 | 0.893 | 0.852 | 0.826 | 0.793 | 0.770 | 0.707 | 0.785 |
| vector | 0.954 | 0.971 | 1.000 | 0.842 | 0.988 | 0.988 | 1.000 | 0.988 | 0.871 | 0.933 | 0.983 |

### Tokens against score

Median completion tokens as billed.

| model | answer score | median tokens | on fully solved | on the rest |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.805 | 1795 | 1588 | 2943 |
| `gemma-4-26b-a4b` | 0.851 | 797 | 765 | 986 |
| `deepseek-v4.1-flash` | 0.970 | 3060 | 2994 | 7719 |
| `qwen3-235b-a22b-2507` | 0.873 | 1013 | 982 | 1302 |
| `glm-5.3-flash` | 0.964 | 2500 | 2396 | 26629 |
| `glm-5.3` | 0.942 | 3270 | 3080 | 32768 |
| `muse-glimmer-30b` | 0.960 | 2988 | 2958 | 4490 |
| `kimi-k3` | 0.964 | 2192 | 2151 | 4513 |
| `gpt-5.4-mini` | 0.837 | 632 | 609 | 811 |
| `gemini-3.1-flash-lite` | 0.858 | 653 | 632 | 792 |
| `claude-sonnet-5` | 0.965 | 1804 | 1783 | 2775 |

### By serving endpoint

Open-weight models were served by several endpoints under the fp8-or-better rule (D-133). The harness dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is confounded with the templates each endpoint happened to serve; the matched difference compares an endpoint with the other endpoints on the same templates (D-149).

| model | endpoint | rows | raw score | unusable | matched difference | templates matched |
|---|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | Darkbloom | 2242 | 0.804 | 0.025 | -0.188 | 8 |
| `gpt-oss-20b` | DekaLLM | 7 | 0.857 | 0.000 | +0.194 | 7 |
| `gpt-oss-20b` | DeepInfra | 1 | 1.000 | 0.000 | +0.143 | 1 |
| `gemma-4-26b-a4b` | DekaLLM | 1263 | 0.849 | 0.000 | -0.002 | 150 |
| `gemma-4-26b-a4b` | NextBit | 987 | 0.854 | 0.000 | +0.002 | 150 |
| `deepseek-v4.1-flash` | CoreWeave | 793 | 0.965 | 0.006 | -0.004 | 143 |
| `deepseek-v4.1-flash` | Parasail | 472 | 0.976 | 0.000 | +0.004 | 140 |
| `deepseek-v4.1-flash` | Morph | 458 | 0.984 | 0.000 | +0.003 | 135 |
| `deepseek-v4.1-flash` | DeepInfra | 419 | 0.951 | 0.019 | -0.010 | 79 |
| `deepseek-v4.1-flash` | Makora | 89 | 1.000 | 0.000 | +0.010 | 44 |
| `deepseek-v4.1-flash` | Novita | 17 | 0.941 | 0.059 | +0.005 | 8 |
| `deepseek-v4.1-flash` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `qwen3-235b-a22b-2507` | GMICloud | 2200 | 0.874 | 0.000 | -0.016 | 36 |
| `qwen3-235b-a22b-2507` | Parasail | 50 | 0.820 | 0.020 | +0.016 | 36 |
| `glm-5.3-flash` | GMICloud | 990 | 0.967 | 0.017 | +0.000 | 146 |
| `glm-5.3-flash` | Novita | 398 | 0.969 | 0.010 | -0.002 | 120 |
| `glm-5.3-flash` | Parasail | 224 | 0.942 | 0.045 | -0.003 | 76 |
| `glm-5.3-flash` | Sail Research | 205 | 0.937 | 0.039 | -0.028 | 116 |
| `glm-5.3-flash` | Phala | 167 | 0.982 | 0.006 | +0.009 | 75 |
| `glm-5.3-flash` | Z.AI | 84 | 0.946 | 0.012 | +0.007 | 49 |
| `glm-5.3-flash` | StreamLake | 77 | 0.987 | 0.000 | -0.006 | 34 |
| `glm-5.3-flash` | NextBit | 43 | 0.965 | 0.023 | +0.023 | 13 |
| `glm-5.3-flash` | AtlasCloud | 42 | 1.000 | 0.000 | +0.014 | 15 |
| `glm-5.3-flash` | SiliconFlow | 19 | 0.947 | 0.000 | -0.012 | 6 |
| `glm-5.3-flash` | Morph | 1 | 1.000 | 0.000 | +0.000 | 1 |
| `glm-5.3` | Baidu | 918 | 0.956 | 0.036 | -0.012 | 86 |
| `glm-5.3` | Sail Research | 748 | 0.937 | 0.037 | +0.004 | 77 |
| `glm-5.3` | Morph | 245 | 0.969 | 0.020 | +0.018 | 122 |
| `glm-5.3` | Novita | 172 | 0.916 | 0.081 | -0.006 | 42 |
| `glm-5.3` | AtlasCloud | 147 | 0.864 | 0.136 | -0.007 | 55 |
| `glm-5.3` | BaseTen | 15 | 1.000 | 0.000 | +0.000 | 5 |
| `glm-5.3` | AkashML | 3 | 1.000 | 0.000 | +0.000 | 3 |
| `glm-5.3` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `kimi-k3` | Morph | 2231 | 0.965 | 0.002 | +0.010 | 7 |
| `kimi-k3` | BaseTen | 19 | 0.895 | 0.000 | -0.010 | 7 |
| `gpt-5.4-mini` | Azure | 2249 | 0.837 | 0.000 | +0.750 | 1 |
| `gpt-5.4-mini` | OpenAI | 1 | 0.000 | 0.000 | -0.750 | 1 |

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

- `qwen3.8-27b`: answer score 0.924 (95% CI 0.888 to 0.954), fully solved 0.918, unusable 89

## Provenance

`analyze.py` at commit `7ebc348` (tag `full-run-evaluation`), dirty; the score store `main` scored at commit `bf4a43b` on 2026-09-29T19:49:30+00:00. The evaluator hashes and the per-model trace hashes are in `results.json`.
