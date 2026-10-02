# Results of the full run

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring, and the corrections and additions made after the first results were read are labelled where they appear (D-146 to D-149, D-170). The answer check is the one corrected after the run was read, D-137 to D-139, D-147 and D-169 (the run's README and `ANSWER_FORM_FIX.md`). The answer check, E3 and E4's digit rule, E5 and the step router. Every interval is 95% and resamples templates (B = 10,000): all 150 for a model's score, and the templates with a qualifying trace for a rate over a subset of traces. Every resampling test draws 100,000 permutations, so its Holm-adjusted floor over 55 pairs is 0.0006. **Incomplete judged stages:** `gemma-4-26b-a4b`, `qwen3-235b-a22b-2507`, `gemini-3.1-flash-lite` have calls without a reply; their judged rates are over the answered calls only.

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0. SD within: the mean over templates of the SD of the score across a template's 15 items; SD between: the SD of the 150 template means (the per-template variability the July rebuttal promised). Unusable is empty plus unreadable; capped: answered rows that stopped at the output cap and are scored on what they state (D-117).

| model | answer score | 95% CI | SD within | SD between | fully solved | 95% CI | correct | partial | incorrect | unusable (empty + unreadable) | capped, scored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.976 | 0.957 to 0.990 | 0.034 | 0.105 | 0.970 | 0.948 to 0.987 | 2182 | 28 | 26 | 14 (14 + 0) | 1 |
| `kimi-k3` | 0.974 | 0.956 to 0.988 | 0.046 | 0.099 | 0.969 | 0.950 to 0.984 | 2180 | 24 | 41 | 5 (5 + 0) | 4 |
| `claude-sonnet-5` | 0.972 | 0.952 to 0.988 | 0.041 | 0.112 | 0.966 | 0.944 to 0.984 | 2173 | 27 | 50 | 0 (0 + 0) | 0 |
| `glm-5.3-flash` | 0.969 | 0.947 to 0.986 | 0.041 | 0.123 | 0.965 | 0.943 to 0.984 | 2172 | 18 | 18 | 42 (42 + 0) | 1 |
| `muse-glimmer-30b` | 0.967 | 0.946 to 0.984 | 0.050 | 0.116 | 0.963 | 0.941 to 0.981 | 2166 | 19 | 36 | 29 (29 + 0) | 0 |
| `glm-5.3` | 0.948 | 0.918 to 0.973 | 0.054 | 0.175 | 0.943 | 0.911 to 0.970 | 2121 | 22 | 7 | 100 (100 + 0) | 6 |
| `qwen3-235b-a22b-2507` | 0.886 | 0.852 to 0.916 | 0.150 | 0.202 | 0.860 | 0.820 to 0.898 | 1935 | 116 | 198 | 1 (0 + 1) | 1 |
| `gemini-3.1-flash-lite` | 0.872 | 0.829 to 0.911 | 0.121 | 0.256 | 0.863 | 0.819 to 0.903 | 1941 | 41 | 268 | 0 (0 + 0) | 0 |
| `gemma-4-26b-a4b` | 0.864 | 0.825 to 0.899 | 0.160 | 0.231 | 0.848 | 0.808 to 0.885 | 1907 | 74 | 269 | 0 (0 + 0) | 1 |
| `gpt-5.4-mini` | 0.847 | 0.806 to 0.885 | 0.168 | 0.252 | 0.834 | 0.792 to 0.873 | 1876 | 59 | 315 | 0 (0 + 0) | 0 |
| `gpt-oss-20b` | 0.814 | 0.770 to 0.856 | 0.207 | 0.270 | 0.800 | 0.755 to 0.842 | 1800 | 65 | 330 | 55 (53 + 2) | 4 |

Of the 55 pairs, 32 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The same template-level test on the fully-solved rate gives the same verdict on 54 of 55. McNemar's exact test on the paired item verdicts (Holm) gives the same verdict on 47 of 55; on 8 of the 8 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level tests. 7 answered rows across the roster ended with a finish reason other than stop or length (a provider fault inside a 200) and are scored on what they state; the harness now retries such a reply (D-148).

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | template p (Holm) | McNemar p (Holm) | same verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.162 | -0.205 to -0.121 | 0.0005 | -0.170 | -0.214 to -0.128 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `kimi-k3` | -0.160 | -0.204 to -0.119 | 0.0005 | -0.169 | -0.212 to -0.127 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.157 | -0.199 to -0.119 | 0.0005 | -0.166 | -0.208 to -0.126 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.155 | -0.197 to -0.116 | 0.0005 | -0.165 | -0.208 to -0.125 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.152 | -0.191 to -0.116 | 0.0005 | -0.163 | -0.204 to -0.124 | 0.0005 | 0.0000 | yes |
| `gpt-oss-20b` | `glm-5.3` | -0.133 | -0.177 to -0.091 | 0.0005 | -0.143 | -0.186 to -0.100 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.129 | 0.092 to 0.169 | 0.0005 | +0.136 | 0.097 to 0.177 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gpt-5.4-mini` | +0.127 | 0.088 to 0.169 | 0.0005 | +0.135 | 0.096 to 0.176 | 0.0005 | 0.0000 | yes |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.125 | -0.164 to -0.088 | 0.0005 | -0.132 | -0.172 to -0.095 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.122 | 0.089 to 0.160 | 0.0005 | +0.132 | 0.096 to 0.169 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.120 | 0.087 to 0.156 | 0.0005 | +0.129 | 0.094 to 0.165 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.112 | -0.146 to -0.081 | 0.0005 | -0.122 | -0.158 to -0.090 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.110 | -0.146 to -0.077 | 0.0005 | -0.121 | -0.157 to -0.088 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.108 | -0.140 to -0.077 | 0.0005 | -0.118 | -0.152 to -0.087 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.105 | -0.136 to -0.077 | 0.0005 | -0.118 | -0.152 to -0.087 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.104 | 0.069 to 0.142 | 0.0005 | +0.107 | 0.073 to 0.145 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.103 | -0.137 to -0.073 | 0.0005 | -0.115 | -0.151 to -0.082 | 0.0005 | 0.0000 | yes |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.102 | 0.068 to 0.140 | 0.0005 | +0.106 | 0.071 to 0.146 | 0.0005 | 0.0000 | yes |
| `glm-5.3` | `gpt-5.4-mini` | +0.101 | 0.069 to 0.135 | 0.0005 | +0.109 | 0.076 to 0.145 | 0.0005 | 0.0000 | yes |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.100 | -0.137 to -0.067 | 0.0005 | -0.103 | -0.141 to -0.070 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.098 | 0.067 to 0.131 | 0.0005 | +0.103 | 0.070 to 0.139 | 0.0005 | 0.0000 | yes |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.095 | 0.063 to 0.130 | 0.0005 | +0.100 | 0.066 to 0.137 | 0.0005 | 0.0000 | yes |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.090 | 0.064 to 0.120 | 0.0005 | +0.110 | 0.076 to 0.147 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.088 | -0.118 to -0.061 | 0.0005 | -0.109 | -0.145 to -0.075 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.086 | -0.113 to -0.061 | 0.0005 | -0.106 | -0.141 to -0.073 | 0.0005 | 0.0000 | yes |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.084 | -0.115 to -0.055 | 0.0005 | -0.095 | -0.128 to -0.066 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.084 | -0.110 to -0.059 | 0.0005 | -0.105 | -0.140 to -0.073 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.081 | -0.107 to -0.057 | 0.0005 | -0.103 | -0.138 to -0.070 | 0.0005 | 0.0000 | yes |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.076 | 0.048 to 0.107 | 0.0005 | +0.080 | 0.052 to 0.112 | 0.0005 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.062 | -0.089 to -0.036 | 0.0005 | -0.083 | -0.119 to -0.049 | 0.0005 | 0.0000 | yes |
| `glm-5.3-flash` | `glm-5.3` | +0.022 | 0.009 to 0.037 | 0.0135 | +0.023 | 0.009 to 0.038 | 0.0230 | 0.0000 | yes |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.071 | -0.113 to -0.030 | 0.0254 | -0.060 | -0.108 to -0.013 | 0.3199 | 0.0000 | no |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.028 | 0.010 to 0.050 | 0.0750 | +0.027 | 0.007 to 0.050 | 0.2551 | 0.0000 | yes |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.057 | -0.103 to -0.012 | 0.3067 | -0.063 | -0.108 to -0.017 | 0.1982 | 0.0000 | yes |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.050 | -0.089 to -0.010 | 0.3067 | -0.048 | -0.088 to -0.009 | 0.4763 | 0.0000 | yes |
| `glm-5.3` | `kimi-k3` | -0.027 | -0.051 to -0.006 | 0.3322 | -0.026 | -0.051 to -0.006 | 0.4763 | 0.0000 | yes |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.039 | 0.007 to 0.073 | 0.4429 | +0.026 | -0.015 to 0.066 | 1.0000 | 0.0585 | yes |
| `glm-5.3` | `claude-sonnet-5` | -0.024 | -0.048 to -0.005 | 0.4808 | -0.023 | -0.047 to -0.004 | 0.8316 | 0.0000 | yes |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.032 | -0.077 to 0.012 | 1.0000 | -0.034 | -0.080 to 0.011 | 1.0000 | 0.0113 | yes |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.025 | -0.054 to 0.002 | 1.0000 | -0.029 | -0.058 to -0.001 | 0.9938 | 0.0027 | yes |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.022 | -0.052 to 0.008 | 1.0000 | -0.012 | -0.048 to 0.025 | 1.0000 | 1.0000 | yes |
| `glm-5.3` | `muse-glimmer-30b` | -0.019 | -0.044 to 0.002 | 1.0000 | -0.020 | -0.046 to 0.002 | 1.0000 | 0.0019 | yes |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.017 | -0.013 to 0.048 | 1.0000 | +0.014 | -0.019 to 0.047 | 1.0000 | 1.0000 | yes |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.014 | -0.015 to 0.045 | 1.0000 | -0.003 | -0.041 to 0.035 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.009 | -0.011 to 0.028 | 1.0000 | +0.007 | -0.013 to 0.026 | 1.0000 | 1.0000 | yes |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.008 | -0.035 to 0.018 | 1.0000 | -0.015 | -0.042 to 0.011 | 1.0000 | 0.7116 | yes |
| `muse-glimmer-30b` | `kimi-k3` | -0.007 | -0.028 to 0.012 | 1.0000 | -0.006 | -0.027 to 0.014 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.007 | -0.004 to 0.018 | 1.0000 | +0.004 | -0.008 to 0.018 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `kimi-k3` | -0.005 | -0.019 to 0.008 | 1.0000 | -0.004 | -0.018 to 0.009 | 1.0000 | 1.0000 | yes |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.005 | -0.025 to 0.015 | 1.0000 | -0.003 | -0.025 to 0.019 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.004 | -0.006 to 0.016 | 1.0000 | +0.004 | -0.009 to 0.018 | 1.0000 | 1.0000 | yes |
| `kimi-k3` | `claude-sonnet-5` | +0.002 | -0.006 to 0.013 | 1.0000 | +0.003 | -0.008 to 0.015 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.002 | -0.016 to 0.020 | 1.0000 | +0.003 | -0.017 to 0.021 | 1.0000 | 1.0000 | yes |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.002 | -0.019 to 0.012 | 1.0000 | -0.000 | -0.017 to 0.014 | 1.0000 | 1.0000 | yes |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.002 | -0.006 to 0.009 | 1.0000 | +0.001 | -0.008 to 0.009 | 1.0000 | 1.0000 | yes |

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier. **Corrected after the first results (D-146):** the planned test permuted the tier labels on the raw difference, which is liberal when the smaller tier has the larger spread, and Advanced template means spread two to four times as widely as Easy ones (the two SD columns); its count of models also moved with the seed. The count now rests on Welch's t-test, Holm across the eleven models: **4 of 11** hold, against 9 under the planned permutation, printed beside it. The detectable gap at 80% power is given as the plan defined it (2.80 x pooled sigma x sqrt(1/58 + 1/34)) and at the strictest Holm step from the Welch standard error.

| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | planned p (Holm) | SD Easy | SD Adv | detectable, planned | detectable, Holm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.905 | 0.716 | +0.189 | 0.079 to 0.307 | 0.0234 | 0.0049 | 0.207 | 0.302 | 0.149 | 0.215 |
| `gemma-4-26b-a4b` | 0.921 | 0.771 | +0.151 | 0.054 to 0.255 | 0.0429 | 0.0076 | 0.175 | 0.270 | 0.130 | 0.190 |
| `deepseek-v4.1-flash` | 0.996 | 0.929 | +0.067 | 0.011 to 0.144 | 0.2906 | 0.0081 | 0.015 | 0.197 | 0.073 | 0.125 |
| `qwen3-235b-a22b-2507` | 0.917 | 0.839 | +0.078 | -0.011 to 0.174 | 0.2906 | 0.0756 | 0.168 | 0.247 | 0.121 | 0.176 |
| `glm-5.3-flash` | 0.995 | 0.922 | +0.073 | 0.020 to 0.145 | 0.1898 | 0.0071 | 0.021 | 0.190 | 0.070 | 0.120 |
| `glm-5.3` | 0.986 | 0.852 | +0.134 | 0.042 to 0.240 | 0.0984 | 0.0053 | 0.051 | 0.300 | 0.112 | 0.191 |
| `muse-glimmer-30b` | 0.982 | 0.927 | +0.055 | 0.004 to 0.112 | 0.2906 | 0.0502 | 0.086 | 0.152 | 0.069 | 0.104 |
| `kimi-k3` | 0.993 | 0.928 | +0.065 | 0.009 to 0.137 | 0.2906 | 0.0168 | 0.020 | 0.196 | 0.072 | 0.124 |
| `gpt-5.4-mini` | 0.918 | 0.734 | +0.184 | 0.076 to 0.306 | 0.0316 | 0.0047 | 0.157 | 0.326 | 0.141 | 0.219 |
| `gemini-3.1-flash-lite` | 0.937 | 0.738 | +0.199 | 0.083 to 0.327 | 0.0315 | 0.0029 | 0.164 | 0.347 | 0.150 | 0.233 |
| `claude-sonnet-5` | 0.993 | 0.923 | +0.070 | 0.006 to 0.151 | 0.2906 | 0.0433 | 0.037 | 0.214 | 0.080 | 0.136 |

**The cliff under two variations.** Unusable rows left out of the template means, because an empty row at the output cap measures finishing within the ceiling as well as solving, and most such rows fall on Advanced templates; and without the nine symbolic templates (D-138). Welch p, Holm across models.

| model | gap, as scored | gap, unusable left out | 95% CI | Welch p (Holm) | gap, no symbolic | 95% CI | Welch p (Holm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | +0.189 | +0.179 | 0.073 to 0.293 | 0.0310 | +0.192 | 0.074 to 0.314 | 0.0389 |
| `gemma-4-26b-a4b` | +0.151 | +0.151 | 0.054 to 0.255 | 0.0429 | +0.149 | 0.045 to 0.259 | 0.0780 |
| `deepseek-v4.1-flash` | +0.067 | +0.055 | 0.003 to 0.128 | 0.5869 | +0.042 | 0.001 to 0.114 | 0.7893 |
| `qwen3-235b-a22b-2507` | +0.078 | +0.078 | -0.011 to 0.174 | 0.5869 | +0.065 | -0.024 to 0.163 | 0.7893 |
| `glm-5.3-flash` | +0.073 | +0.037 | -0.001 to 0.098 | 0.5869 | +0.067 | 0.012 to 0.141 | 0.3842 |
| `glm-5.3` | +0.134 | +0.010 | -0.009 to 0.036 | 0.5869 | +0.128 | 0.028 to 0.242 | 0.2017 |
| `muse-glimmer-30b` | +0.055 | +0.034 | -0.007 to 0.083 | 0.5869 | +0.035 | -0.013 to 0.090 | 0.7893 |
| `kimi-k3` | +0.065 | +0.055 | 0.003 to 0.127 | 0.5869 | +0.036 | -0.003 to 0.094 | 0.7893 |
| `gpt-5.4-mini` | +0.184 | +0.184 | 0.076 to 0.306 | 0.0316 | +0.188 | 0.070 to 0.315 | 0.0609 |
| `gemini-3.1-flash-lite` | +0.199 | +0.199 | 0.083 to 0.327 | 0.0315 | +0.183 | 0.061 to 0.311 | 0.0700 |
| `claude-sonnet-5` | +0.070 | +0.070 | 0.006 to 0.151 | 0.4714 | +0.048 | -0.004 to 0.122 | 0.7893 |

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 70 items with no milestones (all 15 of 3 templates and some items of 12 more); 52 templates have an item with at most one milestone (46 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The floor is the same trace scored against a sibling item's milestones (D-149): what coverage a trace reaches by chance, per model, on the readable wrong answers. The digit rule's flag counts a trace when any step is flagged. Against the experts' step labels on the pilot's fully solved traces it has precision 0.817 and recall 0.427 (0.750 and 0.320 before D-156; 0.800 and 0.427 between D-156 and D-159, unchanged by D-160; SCORER_VALIDATION.md). *Added 2026-09-30 and 2026-10-01 (D-154 to D-165), after the first results were read:* on this roster a domain expert found 111 of 220 sampled flags real before the first fix, 0.505 (FLAG_REVIEW.md), and 155 of 206 real after it, 0.752 (FLAG_REVIEW_2.md). The second fix, built from those notes, took the flag off 42 of their 51 misread steps and left it on 154 of the 155 slips' steps (DIGIT_FIX_2.md). Two review agents then checked both fixes; the corrections they led to, with LaTeX control spaces now read, took the flag off 82 full-run steps and put it on 370, and kept every slip's step the second fix had kept (D-160, DIGIT_FIX_3.md). The same expert then read 190 flags of the rule as it now stands: 171 of the 189 decided are real, 0.905 (95% Wilson 0.854 to 0.939; per model 0.632 to 1.000; FLAG_REVIEW_3.md), and by the owner's rule no fourth reading follows. So its rate is still not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestone coverage over all traces** (*added 2026-10-02, D-170, descriptive and outside the plan's tests*). Per model, the share of the gold's milestones a trace states (E3), or states or the judge rules REACHED (E5-strict), averaged over every trace whose item has milestones, an unusable trace scoring what it reached and an empty one nothing, and over the readable traces alone; template-level intervals. Coverage measures progress through the gold derivation, not the absence of error: on the pilot's correct-answer traces with a flawed step the milestone evaluators scored below chance (RESULTS_X1 Finding 5). The last two columns repeat, from the tables below, the share of fully solved traces the digit rule flags and the share the router's judge flags.

| model | traces with milestones | E3, all | 95% CI | E3, readable | 95% CI | E5-strict, all | 95% CI | E5-strict, readable | 95% CI | digit rule, fully solved | router judge, fully solved |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 2180 | 0.793 | 0.755 to 0.829 | 0.813 | 0.777 to 0.846 | 0.816 | 0.778 to 0.851 | 0.836 | 0.802 to 0.867 | 0.148 | 0.211 |
| `gemma-4-26b-a4b` | 2180 | 0.855 | 0.823 to 0.884 | 0.855 | 0.824 to 0.884 | 0.875 | 0.847 to 0.902 | 0.875 | 0.847 to 0.902 | 0.129 | 0.143 |
| `deepseek-v4.1-flash` | 2180 | 0.860 | 0.830 to 0.888 | 0.865 | 0.836 to 0.893 | 0.896 | 0.870 to 0.921 | 0.902 | 0.876 to 0.925 | 0.005 | 0.011 |
| `qwen3-235b-a22b-2507` | 2180 | 0.869 | 0.840 to 0.896 | 0.869 | 0.841 to 0.896 | 0.893 | 0.868 to 0.917 | 0.893 | 0.868 to 0.917 | 0.225 | 0.109 |
| `glm-5.3-flash` | 2180 | 0.882 | 0.854 to 0.908 | 0.899 | 0.875 to 0.922 | 0.911 | 0.885 to 0.935 | 0.929 | 0.907 to 0.949 | 0.019 | 0.033 |
| `glm-5.3` | 2180 | 0.875 | 0.840 to 0.908 | 0.917 | 0.893 to 0.939 | 0.900 | 0.866 to 0.931 | 0.943 | 0.923 to 0.962 | 0.018 | 0.026 |
| `muse-glimmer-30b` | 2180 | 0.859 | 0.830 to 0.888 | 0.871 | 0.843 to 0.897 | 0.898 | 0.871 to 0.923 | 0.910 | 0.888 to 0.932 | 0.056 | 0.059 |
| `kimi-k3` | 2180 | 0.878 | 0.852 to 0.902 | 0.880 | 0.854 to 0.903 | 0.915 | 0.892 to 0.935 | 0.917 | 0.895 to 0.937 | 0.006 | 0.019 |
| `gpt-5.4-mini` | 2180 | 0.800 | 0.761 to 0.837 | 0.800 | 0.759 to 0.836 | 0.835 | 0.797 to 0.869 | 0.835 | 0.798 to 0.870 | 0.141 | 0.220 |
| `gemini-3.1-flash-lite` | 2180 | 0.840 | 0.808 to 0.870 | 0.840 | 0.807 to 0.871 | 0.865 | 0.835 to 0.892 | 0.865 | 0.835 to 0.893 | 0.070 | 0.122 |
| `claude-sonnet-5` | 2180 | 0.902 | 0.878 to 0.924 | 0.902 | 0.879 to 0.924 | 0.923 | 0.902 to 0.943 | 0.923 | 0.902 to 0.943 | 0.139 | 0.057 |

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | floor, readable | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 385 | 55 | 0.389 | 0.325 to 0.461 | 0.453 | 0.385 to 0.527 | 0.112 | 0.401 (0.336 to 0.474); readable 0.464 |
| `gemma-4-26b-a4b` | 269 | 0 | 0.587 | 0.510 to 0.670 | 0.587 | 0.511 to 0.669 | 0.147 | 0.601 (0.523 to 0.683); readable 0.601 |
| `deepseek-v4.1-flash` | 40 | 14 | 0.398 | 0.226 to 0.602 | 0.612 | 0.479 to 0.739 | 0.123 | 0.440 (0.254 to 0.674); readable 0.676 |
| `qwen3-235b-a22b-2507` | 199 | 1 | 0.580 | 0.501 to 0.668 | 0.583 | 0.506 to 0.670 | 0.200 | 0.636 (0.559 to 0.717); readable 0.636 |
| `glm-5.3-flash` | 60 | 42 | 0.154 | 0.049 to 0.265 | 0.513 | 0.308 to 0.667 | 0.220 | 0.162 (0.052 to 0.277); readable 0.540 |
| `glm-5.3` | 107 | 100 | 0.053 | 0.013 to 0.129 | 0.810 | 0.571 to 1.000 | 0.167 | 0.055 (0.014 to 0.137); readable 0.845 |
| `muse-glimmer-30b` | 65 | 29 | 0.347 | 0.198 to 0.571 | 0.635 | 0.540 to 0.791 | 0.184 | 0.380 (0.219 to 0.615); readable 0.694 |
| `kimi-k3` | 46 | 5 | 0.525 | 0.366 to 0.662 | 0.591 | 0.484 to 0.713 | 0.104 | 0.608 (0.396 to 0.778); readable 0.684 |
| `gpt-5.4-mini` | 315 | 0 | 0.436 | 0.366 to 0.517 | 0.436 | 0.364 to 0.517 | 0.118 | 0.473 (0.398 to 0.559); readable 0.473 |
| `gemini-3.1-flash-lite` | 268 | 0 | 0.510 | 0.441 to 0.581 | 0.510 | 0.440 to 0.581 | 0.110 | 0.537 (0.465 to 0.610); readable 0.537 |
| `claude-sonnet-5` | 50 | 0 | 0.706 | 0.620 to 0.868 | 0.706 | 0.619 to 0.873 | 0.177 | 0.753 (0.661 to 0.883); readable 0.753 |

**E5's judge**, MiMo-V2.5-Pro on the milestones E3 did not find (`judge.py`). E5 coverage is E5-strict: E3's milestones plus those the judge rules REACHED. On the pilot the judge never called a fabricated value REACHED and found 76% of true ones, so the score is conservative (RESULTS_E5). Judged fraction: milestones sent to the judge over those required, on answered traces; unjudged: milestones the reply did not name, over those sent. A trace whose call got no reply is left out of every rate and counted (D-148).

| model | calls | without a reply | judged fraction | of the judged, REACHED | unjudged |
|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 852 | 0 | 0.239 | 0.104 | 0.002 |
| `gemma-4-26b-a4b` | 743 | 0 | 0.187 | 0.125 | 0.000 |
| `deepseek-v4.1-flash` | 775 | 0 | 0.168 | 0.244 | 0.004 |
| `qwen3-235b-a22b-2507` | 717 | 0 | 0.163 | 0.167 | 0.000 |
| `glm-5.3-flash` | 610 | 0 | 0.138 | 0.255 | 0.001 |
| `glm-5.3` | 499 | 0 | 0.112 | 0.277 | 0.000 |
| `muse-glimmer-30b` | 736 | 0 | 0.169 | 0.287 | 0.002 |
| `kimi-k3` | 732 | 0 | 0.156 | 0.269 | 0.001 |
| `gpt-5.4-mini` | 923 | 0 | 0.252 | 0.133 | 0.002 |
| `gemini-3.1-flash-lite` | 824 | 0 | 0.209 | 0.149 | 0.002 |
| `claude-sonnet-5` | 621 | 0 | 0.131 | 0.217 | 0.004 |

**The digit rule.** The wrong-answer rate is over the answered wrong answers (D-149: an empty trace has no step to flag). Beside it, the 1% tolerance E4 shipped with, blind to most slips (RESULTS_X1), and where the first flag falls in a fully solved trace (0 = first step, 1 = last).

| model | wrong answers, answered | flag rate | 95% CI | fully solved traces | flag rate | 95% CI | at 1% | claims checked per trace | traces with a claim | first flag: traces, median position |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 332 | 0.238 | 0.173 to 0.306 | 1800 | 0.148 | 0.117 to 0.181 | 0.021 | 3.74 | 0.782 | 266, 0.667 |
| `gemma-4-26b-a4b` | 269 | 0.450 | 0.330 to 0.569 | 1907 | 0.129 | 0.096 to 0.166 | 0.015 | 4.39 | 0.689 | 246, 0.539 |
| `deepseek-v4.1-flash` | 26 | 0.038 | 0.000 to 0.095 | 2182 | 0.005 | 0.002 to 0.009 | 0.006 | 3.27 | 0.804 | 12, 0.643 |
| `qwen3-235b-a22b-2507` | 199 | 0.477 | 0.343 to 0.606 | 1935 | 0.225 | 0.181 to 0.272 | 0.053 | 9.83 | 0.906 | 436, 0.600 |
| `glm-5.3-flash` | 18 | 0.111 | 0.000 to 0.444 | 2172 | 0.019 | 0.012 to 0.027 | 0.015 | 4.75 | 0.851 | 42, 0.750 |
| `glm-5.3` | 7 | 0.000 | 0.000 to 0.000 | 2121 | 0.018 | 0.012 to 0.026 | 0.017 | 5.72 | 0.915 | 39, 0.727 |
| `muse-glimmer-30b` | 36 | 0.028 | 0.000 to 0.086 | 2166 | 0.056 | 0.042 to 0.071 | 0.020 | 3.26 | 0.706 | 122, 0.500 |
| `kimi-k3` | 41 | 0.000 | 0.000 to 0.000 | 2180 | 0.006 | 0.003 to 0.010 | 0.004 | 2.80 | 0.722 | 14, 0.563 |
| `gpt-5.4-mini` | 315 | 0.276 | 0.182 to 0.377 | 1876 | 0.141 | 0.108 to 0.178 | 0.012 | 3.79 | 0.797 | 265, 0.600 |
| `gemini-3.1-flash-lite` | 268 | 0.321 | 0.183 to 0.458 | 1941 | 0.070 | 0.046 to 0.095 | 0.014 | 3.86 | 0.757 | 135, 0.500 |
| `claude-sonnet-5` | 50 | 0.300 | 0.125 to 0.491 | 2173 | 0.139 | 0.107 to 0.171 | 0.015 | 5.30 | 0.857 | 301, 0.667 |

**Wrong-answer rate against the item's milestone count** (D-149): how failure grows with the depth of the gold derivation.

| model | 0 (70 items) | 1 (290 items) | 2 (452 items) | 3 (367 items) | 4-5 (532 items) | 6+ (539 items) |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.143 | 0.034 | 0.088 | 0.188 | 0.222 | 0.256 |
| `gemma-4-26b-a4b` | 0.000 | 0.045 | 0.082 | 0.104 | 0.085 | 0.252 |
| `deepseek-v4.1-flash` | 0.000 | 0.000 | 0.000 | 0.011 | 0.015 | 0.052 |
| `qwen3-235b-a22b-2507` | 0.057 | 0.041 | 0.058 | 0.082 | 0.081 | 0.156 |
| `glm-5.3-flash` | 0.000 | 0.003 | 0.013 | 0.008 | 0.019 | 0.074 |
| `glm-5.3` | 0.000 | 0.003 | 0.011 | 0.025 | 0.049 | 0.122 |
| `muse-glimmer-30b` | 0.014 | 0.000 | 0.015 | 0.033 | 0.030 | 0.054 |
| `kimi-k3` | 0.014 | 0.000 | 0.007 | 0.019 | 0.028 | 0.037 |
| `gpt-5.4-mini` | 0.000 | 0.028 | 0.075 | 0.090 | 0.156 | 0.291 |
| `gemini-3.1-flash-lite` | 0.029 | 0.021 | 0.049 | 0.076 | 0.143 | 0.249 |
| `claude-sonnet-5` | 0.000 | 0.000 | 0.018 | 0.014 | 0.023 | 0.046 |

**The step router** (`router.py`): the digit rule's flags, and MiMo-V2.5-Pro on every other step in one batched call per trace. On the pilot's labelled traces it had precision 0.707 and recall 0.603 over all steps, and 0.703 and 0.360 inside correct-answer traces (ROUTER_VALIDATION.md), so its rates are flags, not counts of errors. A trace counts as flagged when any step is; a trace whose call got no reply is left out and counted.

| model | calls | without a reply | flagged, fully solved | 95% CI | by the judge | 95% CI | flagged, wrong answers | 95% CI | steps flagged per trace | unjudged steps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 2197 | 0 | 0.315 | 0.278 to 0.353 | 0.211 | 0.182 to 0.241 | 0.931 | 0.887 to 0.965 | 1.00 | 0.001 |
| `gemma-4-26b-a4b` | 2249 | 1 | 0.244 | 0.200 to 0.291 | 0.143 | 0.110 to 0.179 | 0.885 | 0.806 to 0.947 | 0.66 | 0.000 |
| `deepseek-v4.1-flash` | 2236 | 0 | 0.016 | 0.010 to 0.022 | 0.011 | 0.006 to 0.016 | 0.423 | 0.053 to 0.733 | 0.03 | 0.000 |
| `qwen3-235b-a22b-2507` | 2250 | 1 | 0.292 | 0.246 to 0.340 | 0.109 | 0.085 to 0.134 | 0.849 | 0.766 to 0.917 | 0.84 | 0.001 |
| `glm-5.3-flash` | 2208 | 0 | 0.051 | 0.039 to 0.063 | 0.033 | 0.024 to 0.043 | 0.611 | 0.100 to 0.824 | 0.08 | 0.000 |
| `glm-5.3` | 2150 | 0 | 0.042 | 0.030 to 0.055 | 0.026 | 0.016 to 0.037 | 0.143 | 0.000 to 0.500 | 0.06 | 0.000 |
| `muse-glimmer-30b` | 2221 | 0 | 0.110 | 0.092 to 0.130 | 0.059 | 0.045 to 0.073 | 0.500 | 0.261 to 0.722 | 0.15 | 0.000 |
| `kimi-k3` | 2245 | 0 | 0.026 | 0.018 to 0.034 | 0.019 | 0.012 to 0.027 | 0.488 | 0.167 to 0.828 | 0.06 | 0.000 |
| `gpt-5.4-mini` | 2250 | 0 | 0.314 | 0.264 to 0.365 | 0.220 | 0.177 to 0.264 | 0.886 | 0.815 to 0.937 | 0.78 | 0.000 |
| `gemini-3.1-flash-lite` | 2250 | 1 | 0.176 | 0.137 to 0.217 | 0.122 | 0.090 to 0.157 | 0.861 | 0.743 to 0.951 | 0.49 | 0.000 |
| `claude-sonnet-5` | 2250 | 0 | 0.185 | 0.153 to 0.218 | 0.057 | 0.043 to 0.071 | 0.640 | 0.375 to 0.875 | 0.28 | 0.000 |

**Also reported, not tested: what points at a wrong answer.** On the answered traces that score 0, the share with a digit-rule flag, with a milestone E5 rules MISSING, and with a step the router's judge flags. A trace can be in several columns, or in none.

| model | answered wrong-answer traces | digit rule | E5 MISSING | router judge |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 332 | 0.238 | 0.786 | 0.895 |
| `gemma-4-26b-a4b` | 269 | 0.450 | 0.643 | 0.743 |
| `deepseek-v4.1-flash` | 26 | 0.038 | 0.462 | 0.385 |
| `qwen3-235b-a22b-2507` | 199 | 0.477 | 0.593 | 0.668 |
| `glm-5.3-flash` | 18 | 0.111 | 0.444 | 0.611 |
| `glm-5.3` | 7 | 0.000 | 0.143 | 0.143 |
| `muse-glimmer-30b` | 36 | 0.028 | 0.611 | 0.472 |
| `kimi-k3` | 41 | 0.000 | 0.366 | 0.488 |
| `gpt-5.4-mini` | 315 | 0.276 | 0.771 | 0.848 |
| `gemini-3.1-flash-lite` | 267 | 0.322 | 0.730 | 0.753 |
| `claude-sonnet-5` | 50 | 0.300 | 0.240 | 0.400 |

## Q4. Consistency within a template

The share of templates fully solved on all 15 instances, on some, and on none: the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) and the 92 others. Intervals are in `results.json`.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.362 | 0.569 | 0.069 | 0.337 | 0.652 | 0.011 |
| `gemma-4-26b-a4b` | 0.466 | 0.500 | 0.034 | 0.511 | 0.478 | 0.011 |
| `deepseek-v4.1-flash` | 0.897 | 0.103 | 0.000 | 0.848 | 0.152 | 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.534 | 0.017 | 0.543 | 0.435 | 0.022 |
| `glm-5.3-flash` | 0.914 | 0.086 | 0.000 | 0.826 | 0.174 | 0.000 |
| `glm-5.3` | 0.845 | 0.138 | 0.017 | 0.772 | 0.217 | 0.011 |
| `muse-glimmer-30b` | 0.862 | 0.138 | 0.000 | 0.826 | 0.174 | 0.000 |
| `kimi-k3` | 0.862 | 0.138 | 0.000 | 0.793 | 0.207 | 0.000 |
| `gpt-5.4-mini` | 0.483 | 0.483 | 0.034 | 0.467 | 0.511 | 0.022 |
| `gemini-3.1-flash-lite` | 0.672 | 0.276 | 0.052 | 0.598 | 0.391 | 0.011 |
| `claude-sonnet-5` | 0.897 | 0.103 | 0.000 | 0.826 | 0.174 | 0.000 |

## Q5. Paraphrase robustness

The experts' check: 316 pairs returned, 39 rejected and dropped from both arms, 0 not yet returned and kept provisionally.

Paraphrase minus original, paired by item; besides the answer score, E3 coverage on the items with milestones, E5-strict where both arms carry it, and the answer score on the pairs both arms served from the same endpoint (D-149).

| model | items | answer score diff | 95% CI | p (Holm) | McNemar p (Holm) | E3 coverage diff | 95% CI | E5 diff | 95% CI | same provider: pairs, diff |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 277 | +0.011 | -0.036 to 0.057 | 1.0000 | 1.0000 | -0.009 | -0.042 to 0.024 | -0.014 | -0.046 to 0.017 | 275, +0.011 |
| `gemma-4-26b-a4b` | 277 | +0.007 | -0.035 to 0.047 | 1.0000 | 1.0000 | +0.022 | 0.002 to 0.045 | +0.018 | 0.000 to 0.037 | 166, +0.015 |
| `deepseek-v4.1-flash` | 277 | -0.007 | -0.029 to 0.015 | 1.0000 | 1.0000 | +0.007 | -0.014 to 0.029 | +0.000 | -0.019 to 0.022 | 31, +0.000 |
| `qwen3-235b-a22b-2507` | 277 | -0.031 | -0.074 to 0.011 | 1.0000 | 1.0000 | -0.005 | -0.027 to 0.017 | -0.010 | -0.033 to 0.012 | 186, -0.030 |
| `glm-5.3-flash` | 277 | -0.022 | -0.043 to -0.002 | 0.6473 | 0.9844 | -0.031 | -0.053 to -0.009 | -0.009 | -0.029 to 0.009 | 23, +0.000 |
| `glm-5.3` | 277 | -0.007 | -0.032 to 0.016 | 1.0000 | 1.0000 | +0.000 | -0.021 to 0.021 | -0.009 | -0.030 to 0.011 | 85, +0.012 |
| `muse-glimmer-30b` | 277 | +0.020 | 0.003 to 0.040 | 0.6473 | 0.7031 | -0.014 | -0.041 to 0.011 | -0.010 | -0.037 to 0.014 | 277, +0.020 |
| `kimi-k3` | 277 | +0.000 | -0.009 to 0.011 | 1.0000 | 1.0000 | -0.008 | -0.029 to 0.012 | -0.004 | -0.024 to 0.016 | 3, +0.000 |
| `gpt-5.4-mini` | 277 | -0.018 | -0.053 to 0.018 | 1.0000 | 1.0000 | -0.027 | -0.053 to -0.003 | -0.011 | -0.038 to 0.015 | 277, -0.018 |
| `gemini-3.1-flash-lite` | 277 | +0.007 | -0.032 to 0.046 | 1.0000 | 1.0000 | +0.006 | -0.012 to 0.025 | +0.003 | -0.014 to 0.021 | 277, +0.007 |
| `claude-sonnet-5` | 277 | -0.020 | -0.037 to -0.006 | 0.3414 | 0.6875 | -0.008 | -0.025 to 0.008 | -0.006 | -0.023 to 0.012 | 277, -0.020 |

Kendall's tau between the models' answer scores on the originals and on the paraphrases, over the 277 items every tested model holds: 0.673, 95% CI 0.455 to 0.881; the noise floor for tau on this roster is 0.891 (below).

The noise floor at the arm's own size and template mix (D-165): two disjoint draws of 277 main-run items with the arm's per-template counts over its 115 templates, the models ordered on each, 200 draws: median tau 0.782, quartiles 0.721 to 0.855, 5th to 95th percentile 0.624 to 0.927. This, not the roster-wide floor (halves of 7 or 8 items per template), is the comparison for the arm's tau.

**Bounds (D-165).** The same bootstrap's 90% interval per model, read against a margin of ±0.05 in answer score: a model is within the margin when the whole interval is (the two-one-sided-tests rule at 5%). The margin is the plan's largest detectable paired difference (5.1 points at 15% discordance); it was fixed after the point estimates were known and before these intervals were computed.

| model | answer score diff | 90% CI | within the margin |
|---|---:|---:|---|
| `gpt-oss-20b` | +0.011 | -0.028 to 0.050 | yes |
| `gemma-4-26b-a4b` | +0.007 | -0.028 to 0.040 | yes |
| `deepseek-v4.1-flash` | -0.007 | -0.025 to 0.011 | yes |
| `qwen3-235b-a22b-2507` | -0.031 | -0.067 to 0.005 | no |
| `glm-5.3-flash` | -0.022 | -0.039 to -0.005 | yes |
| `glm-5.3` | -0.007 | -0.028 to 0.012 | yes |
| `muse-glimmer-30b` | +0.020 | 0.005 to 0.036 | yes |
| `kimi-k3` | +0.000 | -0.008 to 0.009 | yes |
| `gpt-5.4-mini` | -0.018 | -0.047 to 0.011 | yes |
| `gemini-3.1-flash-lite` | +0.007 | -0.025 to 0.040 | yes |
| `claude-sonnet-5` | -0.020 | -0.034 to -0.007 | yes |

**Against run-to-run noise (D-165).** For the models with decoding repeats: the paraphrase difference beside each repeat minus the main run on the 300 repeat items, each paired by item with the same template bootstrap; and both on the kept pairs the two subsamples share.

| model | paraphrase − original (pairs) | repeat1 − main (300) | repeat2 − main (300) | repeat3 − main (300) | shared items | paraphrase on them | repeats on them | |paraphrase| within the repeats' spread |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | +0.011 -0.036 to 0.057 (277) | +0.003 -0.030 to 0.035 | +0.012 -0.022 to 0.043 | +0.007 -0.028 to 0.040 | 92 | +0.005 -0.065 to 0.076 | -0.027, +0.000, -0.005 | yes |
| `gemma-4-26b-a4b` | +0.007 -0.035 to 0.047 (277) | -0.020 -0.055 to 0.015 | -0.003 -0.037 to 0.028 | +0.005 -0.020 to 0.030 | 92 | -0.038 -0.103 to 0.027 | -0.022, +0.005, -0.027 | yes |
| `qwen3-235b-a22b-2507` | -0.031 -0.074 to 0.011 (277) | -0.002 -0.030 to 0.028 | +0.005 -0.023 to 0.035 | +0.005 -0.025 to 0.033 | 92 | -0.027 -0.082 to 0.027 | -0.005, +0.000, +0.027 | no |
| `gemini-3.1-flash-lite` | +0.007 -0.032 to 0.046 (277) | -0.005 -0.027 to 0.017 | -0.012 -0.035 to 0.012 | +0.007 -0.015 to 0.030 | 92 | -0.027 -0.076 to 0.016 | -0.027, -0.033, -0.033 | yes |

## Sensitivity

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one. The three added readings (D-147, D-149) the experts could not arbitrate: the half-unit window requires a correct rounding at the precision shown where the rule accepts one unit either way; the whole-trace reading credits a quantity the question asks for when it is stated in the body and left off the Answer line (the prompt asks for it there); the pool without the nine symbolic templates, whose answers the check scores by the numbers they state (D-138). "Unusable excluded" is an item mean over the usable rows, not a mean of template means. Tau's noise floor: the ordering on one random half of each template's items against the other, median 0.891 over 200 splits (quartiles 0.855 to 0.927); the top five models lie within 0.012 of each other, so tau falls below 1 from sampling noise alone.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates | without the 9 symbolic templates | half-unit window | whole trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.807 | 0.814 | 0.822 | 0.800 | 0.835 | 0.813 | 0.811 | 0.785 | 0.823 |
| `gemma-4-26b-a4b` | 0.844 | 0.864 | 0.895 | 0.848 | 0.864 | 0.861 | 0.864 | 0.849 | 0.876 |
| `deepseek-v4.1-flash` | 0.972 | 0.976 | 0.977 | 0.970 | 0.982 | 0.976 | 0.983 | 0.973 | 0.979 |
| `qwen3-235b-a22b-2507` | 0.876 | 0.886 | 0.897 | 0.860 | 0.886 | 0.883 | 0.889 | 0.876 | 0.907 |
| `glm-5.3-flash` | 0.968 | 0.969 | 0.970 | 0.965 | 0.988 | 0.968 | 0.970 | 0.965 | 0.972 |
| `glm-5.3` | 0.946 | 0.948 | 0.948 | 0.943 | 0.992 | 0.946 | 0.950 | 0.943 | 0.951 |
| `muse-glimmer-30b` | 0.965 | 0.967 | 0.968 | 0.963 | 0.980 | 0.966 | 0.971 | 0.955 | 0.968 |
| `kimi-k3` | 0.972 | 0.974 | 0.976 | 0.969 | 0.976 | 0.974 | 0.982 | 0.970 | 0.977 |
| `gpt-5.4-mini` | 0.836 | 0.847 | 0.859 | 0.834 | 0.847 | 0.843 | 0.845 | 0.813 | 0.853 |
| `gemini-3.1-flash-lite` | 0.863 | 0.872 | 0.885 | 0.863 | 0.872 | 0.869 | 0.876 | 0.865 | 0.879 |
| `claude-sonnet-5` | 0.966 | 0.972 | 0.977 | 0.966 | 0.972 | 0.971 | 0.978 | 0.965 | 0.976 |
| tau with the headline | 0.927 | 1 | 0.927 | 0.964 | 0.636 | 1.000 | 0.964 | 0.991 | 1.000 |

The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).

## Also reported, not tested

### By branch and level

| model | chemical | civil | electrical | industrial | mechanical | Easy | Intermediate | Advanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.809 | 0.673 | 0.908 | 0.820 | 0.862 | 0.905 | 0.782 | 0.716 |
| `gemma-4-26b-a4b` | 0.734 | 0.882 | 0.904 | 0.878 | 0.921 | 0.921 | 0.861 | 0.771 |
| `deepseek-v4.1-flash` | 0.956 | 0.993 | 0.954 | 0.980 | 0.997 | 0.996 | 0.983 | 0.929 |
| `qwen3-235b-a22b-2507` | 0.827 | 0.940 | 0.886 | 0.871 | 0.906 | 0.917 | 0.882 | 0.839 |
| `glm-5.3-flash` | 0.947 | 0.984 | 0.981 | 0.942 | 0.992 | 0.995 | 0.972 | 0.922 |
| `glm-5.3` | 0.898 | 0.973 | 0.951 | 0.922 | 0.993 | 0.986 | 0.965 | 0.852 |
| `muse-glimmer-30b` | 0.963 | 0.969 | 0.971 | 0.940 | 0.991 | 0.982 | 0.975 | 0.927 |
| `kimi-k3` | 0.964 | 0.987 | 0.949 | 0.980 | 0.991 | 0.993 | 0.982 | 0.928 |
| `gpt-5.4-mini` | 0.800 | 0.838 | 0.911 | 0.831 | 0.854 | 0.918 | 0.841 | 0.734 |
| `gemini-3.1-flash-lite` | 0.776 | 0.862 | 0.900 | 0.860 | 0.961 | 0.937 | 0.885 | 0.738 |
| `claude-sonnet-5` | 0.949 | 0.996 | 0.947 | 0.973 | 0.994 | 0.993 | 0.980 | 0.923 |

### By domain

| domain | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| digital_communications | 0.930 | 0.883 | 0.923 | 0.840 | 0.967 | 0.923 | 0.950 | 0.897 | 0.923 | 0.877 | 0.913 |
| electromagnetics_and_waves | 0.910 | 0.887 | 0.990 | 0.883 | 0.977 | 0.947 | 0.997 | 0.980 | 0.830 | 0.883 | 0.953 |
| fluid_mechanics | 0.880 | 0.893 | 0.997 | 0.913 | 0.987 | 0.987 | 0.987 | 0.987 | 0.950 | 0.950 | 0.987 |
| geotechnical_engineering | 0.520 | 0.913 | 1.000 | 0.947 | 0.993 | 0.993 | 0.933 | 1.000 | 0.913 | 0.960 | 0.987 |
| mechanics_of_materials | 0.840 | 0.907 | 1.000 | 0.943 | 0.993 | 1.000 | 0.997 | 0.993 | 0.803 | 0.940 | 1.000 |
| production_and_inventory | 0.767 | 0.873 | 0.940 | 0.847 | 0.887 | 0.887 | 0.887 | 0.967 | 0.780 | 0.853 | 0.987 |
| quality_and_reliability_control | 0.773 | 0.800 | 1.000 | 0.787 | 0.947 | 0.907 | 0.933 | 0.980 | 0.760 | 0.800 | 0.940 |
| reaction_kinetics | 0.967 | 0.913 | 0.993 | 0.980 | 1.000 | 0.993 | 1.000 | 1.000 | 0.913 | 0.913 | 1.000 |
| signals_and_systems | 0.883 | 0.943 | 0.950 | 0.933 | 1.000 | 0.983 | 0.967 | 0.970 | 0.980 | 0.940 | 0.973 |
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
| multipart | 0.886 | 0.904 | 0.994 | 0.905 | 0.986 | 0.958 | 0.981 | 0.990 | 0.906 | 0.921 | 0.989 |
| scalar | 0.772 | 0.834 | 0.980 | 0.885 | 0.966 | 0.956 | 0.965 | 0.981 | 0.816 | 0.857 | 0.972 |
| symbolic | 0.870 | 0.867 | 0.874 | 0.830 | 0.952 | 0.907 | 0.907 | 0.856 | 0.881 | 0.811 | 0.867 |
| vector | 0.954 | 0.971 | 1.000 | 0.842 | 0.988 | 0.988 | 1.000 | 0.988 | 0.871 | 0.933 | 0.983 |

### Tokens against score

Median completion tokens as billed.

| model | answer score | median tokens | on fully solved | on the rest |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.814 | 1795 | 1588 | 3060 |
| `gemma-4-26b-a4b` | 0.864 | 797 | 767 | 1044 |
| `deepseek-v4.1-flash` | 0.976 | 3060 | 2991 | 12171 |
| `qwen3-235b-a22b-2507` | 0.886 | 1013 | 986 | 1334 |
| `glm-5.3-flash` | 0.969 | 2500 | 2402 | 32768 |
| `glm-5.3` | 0.948 | 3270 | 3072 | 32768 |
| `muse-glimmer-30b` | 0.967 | 2988 | 2952 | 6134 |
| `kimi-k3` | 0.974 | 2192 | 2164 | 6033 |
| `gpt-5.4-mini` | 0.847 | 632 | 610 | 836 |
| `gemini-3.1-flash-lite` | 0.872 | 653 | 632 | 810 |
| `claude-sonnet-5` | 0.972 | 1804 | 1780 | 3346 |

### By serving endpoint

Open-weight models were served by several endpoints under the fp8-or-better rule (D-133). The harness dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is confounded with the templates each endpoint happened to serve; the matched difference compares an endpoint with the other endpoints on the same templates (D-149).

| model | endpoint | rows | raw score | unusable | matched difference | templates matched |
|---|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | Darkbloom | 2242 | 0.814 | 0.025 | -0.223 | 8 |
| `gpt-oss-20b` | DekaLLM | 7 | 1.000 | 0.000 | +0.235 | 7 |
| `gpt-oss-20b` | DeepInfra | 1 | 1.000 | 0.000 | +0.143 | 1 |
| `gemma-4-26b-a4b` | DekaLLM | 1263 | 0.865 | 0.000 | +0.000 | 150 |
| `gemma-4-26b-a4b` | NextBit | 987 | 0.863 | 0.000 | -0.000 | 150 |
| `deepseek-v4.1-flash` | CoreWeave | 793 | 0.977 | 0.006 | -0.003 | 143 |
| `deepseek-v4.1-flash` | Parasail | 472 | 0.981 | 0.000 | +0.001 | 140 |
| `deepseek-v4.1-flash` | Morph | 458 | 0.986 | 0.000 | +0.003 | 135 |
| `deepseek-v4.1-flash` | DeepInfra | 419 | 0.953 | 0.019 | -0.003 | 79 |
| `deepseek-v4.1-flash` | Makora | 89 | 1.000 | 0.000 | +0.010 | 44 |
| `deepseek-v4.1-flash` | Novita | 17 | 0.941 | 0.059 | +0.005 | 8 |
| `deepseek-v4.1-flash` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `qwen3-235b-a22b-2507` | GMICloud | 2200 | 0.885 | 0.000 | -0.011 | 36 |
| `qwen3-235b-a22b-2507` | Parasail | 50 | 0.900 | 0.020 | +0.011 | 36 |
| `glm-5.3-flash` | GMICloud | 990 | 0.974 | 0.017 | +0.002 | 146 |
| `glm-5.3-flash` | Novita | 398 | 0.969 | 0.010 | -0.005 | 120 |
| `glm-5.3-flash` | Parasail | 224 | 0.955 | 0.045 | -0.001 | 76 |
| `glm-5.3-flash` | Sail Research | 205 | 0.946 | 0.039 | -0.021 | 116 |
| `glm-5.3-flash` | Phala | 167 | 0.988 | 0.006 | +0.005 | 75 |
| `glm-5.3-flash` | Z.AI | 84 | 0.946 | 0.012 | +0.007 | 49 |
| `glm-5.3-flash` | StreamLake | 77 | 0.987 | 0.000 | -0.017 | 34 |
| `glm-5.3-flash` | NextBit | 43 | 0.965 | 0.023 | +0.023 | 13 |
| `glm-5.3-flash` | AtlasCloud | 42 | 1.000 | 0.000 | +0.014 | 15 |
| `glm-5.3-flash` | SiliconFlow | 19 | 0.947 | 0.000 | -0.012 | 6 |
| `glm-5.3-flash` | Morph | 1 | 1.000 | 0.000 | +0.000 | 1 |
| `glm-5.3` | Baidu | 918 | 0.957 | 0.036 | -0.006 | 86 |
| `glm-5.3` | Sail Research | 748 | 0.951 | 0.037 | +0.001 | 77 |
| `glm-5.3` | Morph | 245 | 0.969 | 0.020 | +0.016 | 122 |
| `glm-5.3` | Novita | 172 | 0.916 | 0.081 | -0.006 | 42 |
| `glm-5.3` | AtlasCloud | 147 | 0.864 | 0.136 | -0.007 | 55 |
| `glm-5.3` | BaseTen | 15 | 1.000 | 0.000 | +0.000 | 5 |
| `glm-5.3` | AkashML | 3 | 1.000 | 0.000 | +0.000 | 3 |
| `glm-5.3` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `kimi-k3` | Morph | 2231 | 0.975 | 0.002 | +0.010 | 7 |
| `kimi-k3` | BaseTen | 19 | 0.895 | 0.000 | -0.010 | 7 |
| `gpt-5.4-mini` | Azure | 2249 | 0.847 | 0.000 | +0.750 | 1 |
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

## Decoding repeats

| model | items | repeat1 | repeat2 | repeat3 | main run, same items | SD | range | same verdict in every repeat |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 300 | 0.833 | 0.842 | 0.837 | 0.830 | 0.004 | 0.008 | 0.920 |
| `gemma-4-26b-a4b` | 300 | 0.857 | 0.873 | 0.882 | 0.877 | 0.013 | 0.025 | 0.857 |
| `qwen3-235b-a22b-2507` | 300 | 0.893 | 0.900 | 0.900 | 0.895 | 0.004 | 0.007 | 0.887 |
| `gemini-3.1-flash-lite` | 300 | 0.890 | 0.883 | 0.902 | 0.895 | 0.009 | 0.018 | 0.937 |

## Set aside (D-132), outside every comparison

- `qwen3.8-27b`: answer score 0.937 (95% CI 0.906 to 0.963), fully solved 0.931, unusable 89

## Provenance

`analyze.py` at commit `2f388c3`, dirty; the score store `main` scored at commit `2950875` on 2026-10-02T13:09:33+00:00; stages: e5 at `2950875`, router at `2950875`. The evaluator hashes and the per-model trace hashes are in `results.json`.
