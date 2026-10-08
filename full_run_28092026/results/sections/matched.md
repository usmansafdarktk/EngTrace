# Matched settings

**QUICK**: development values, 1,000 bootstrap draws and 10,000 sign flips or permutations (the full pass draws 10,000 and 100,000); the integration pass replaces every number here.

Printed by `analyze.py` (WS-C1) from `results/matched_config.json` (written 2026-10-07T13:20:31Z, sha256 `e1e4d06b307e`): each model at its reasoning-on store where the file names one, else at its default store. With reasoning on: `gpt-5.4-mini` (effort=medium, `scores/reasoning-medium-full`), `gemma-4-26b-a4b` (effort=medium, `scores/reasoning-medium-full`), `gemini-3.1-flash-lite` (effort=medium, `scores/reasoning-medium-full`); the others at `scores/main`. Skipped: `qwen3-235b-a22b-thinking-2507` (run: false). The functions and seeds are those of Q1, Q2 and Q3, so a model whose rows did not change shows the numbers of the default configuration. FAC is the answer score (correct 1, partial 0.5, incorrect or unusable 0), the mean of the 150 template means; every interval is 95% and resamples templates. The letters are the compact letter display of the sign-flip pairs at a Holm-adjusted p below 0.05: models that share a letter are not separated. No readable answer: empty, or answered without an answer the check can read. MC is Milestone Coverage, the mean of the template means over the templates every model has coverage on (147 with the judge, 147 by matching alone): with the judge (E5-strict) and by matching alone (E3). Judge-decided share: the milestones sent to the judge over those required, on answered responses. The flag rates are on fully solved responses: the arithmetic check's flags and the step check's judge flags, template intervals.

| model | store | FAC | 95% CI | letter | letter, default | no readable answer | MC with the judge | 95% CI | letter | MC by matching alone | 95% CI | judge-decided share | Easy | Intermediate | Advanced | gap | 95% CI | Welch p (Holm) | arithmetic flags | 95% CI | calculations parsed per response | judged step flags | 95% CI |
|---|---|---:|---:|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | main | 0.981 | 0.967 to 0.991 | a | a | 0.004 (8) | 0.900 | 0.876 to 0.923 | abc | 0.862 | 0.834 to 0.889 | 0.167 | 0.996 | 0.983 | 0.951 | +0.045 | 0.009 to 0.096 | 0.3399 | 0.005 | 0.002 to 0.009 | 3.27 | 0.010 | 0.005 to 0.015 |
| `claude-sonnet-5` | main | 0.981 | 0.967 to 0.991 | ab | ab | 0.001 (2) | 0.924 | 0.904 to 0.943 | a | 0.903 | 0.881 to 0.924 | 0.125 | 0.993 | 0.980 | 0.962 | +0.031 | -0.004 to 0.084 | 0.5247 | 0.143 | 0.112 to 0.180 | 5.32 | 0.055 | 0.042 to 0.068 |
| `kimi-k3` | main | 0.974 | 0.957 to 0.988 | abcd | ab | 0.007 (16) | 0.913 | 0.889 to 0.934 | ab | 0.876 | 0.849 to 0.900 | 0.150 | 0.993 | 0.982 | 0.926 | +0.067 | 0.010 to 0.140 | 0.3399 | 0.007 | 0.003 to 0.011 | 2.79 | 0.017 | 0.011 to 0.024 |
| `glm-5.3-flash` | main | 0.972 | 0.954 to 0.988 | ac | a | 0.020 (45) | 0.912 | 0.882 to 0.938 | abc | 0.882 | 0.849 to 0.911 | 0.135 | 0.995 | 0.972 | 0.935 | +0.060 | 0.017 to 0.114 | 0.1957 | 0.019 | 0.012 to 0.027 | 4.74 | 0.032 | 0.023 to 0.042 |
| `muse-glimmer-30b` | main | 0.969 | 0.950 to 0.985 | abcd | ab | 0.014 (31) | 0.898 | 0.872 to 0.924 | abc | 0.860 | 0.831 to 0.889 | 0.168 | 0.982 | 0.975 | 0.935 | +0.047 | -0.001 to 0.099 | 0.3412 | 0.059 | 0.044 to 0.075 | 3.27 | 0.060 | 0.046 to 0.075 |
| `gpt-5.4-mini` (reasoning on) | reasoning-medium-full | 0.960 | 0.944 to 0.976 | cde | cd | 0.000 (0) | 0.892 | 0.869 to 0.918 | bc | 0.848 | 0.820 to 0.879 | 0.183 | 0.982 | 0.946 | 0.947 | +0.035 | -0.003 to 0.087 | 0.5247 | 0.076 | 0.058 to 0.096 | 3.27 | 0.115 | 0.092 to 0.139 |
| `glm-5.3` | main | 0.948 | 0.918 to 0.972 | bdef | b | 0.044 (99) | 0.901 | 0.868 to 0.931 | abc | 0.876 | 0.842 to 0.908 | 0.113 | 0.986 | 0.965 | 0.854 | +0.132 | 0.042 to 0.241 | 0.1271 | 0.018 | 0.012 to 0.026 | 5.72 | 0.026 | 0.017 to 0.037 |
| `gemma-4-26b-a4b` (reasoning on) | reasoning-medium-full | 0.932 | 0.905 to 0.956 | efg | cd | 0.001 (3) | 0.885 | 0.857 to 0.912 | c | 0.852 | 0.824 to 0.882 | 0.185 | 0.969 | 0.923 | 0.886 | +0.083 | 0.016 to 0.154 | 0.2082 | 0.100 | 0.064 to 0.138 | 3.30 | 0.075 | 0.053 to 0.098 |
| `gemini-3.1-flash-lite` (reasoning on) | reasoning-medium-full | 0.910 | 0.875 to 0.938 | fg | cd | 0.000 (0) | 0.880 | 0.848 to 0.907 | c | 0.856 | 0.823 to 0.887 | 0.186 | 0.960 | 0.922 | 0.803 | +0.157 | 0.062 to 0.268 | 0.0395 | 0.070 | 0.049 to 0.097 | 3.54 | 0.094 | 0.068 to 0.121 |
| `qwen3-235b-a22b-2507` | main | 0.894 | 0.865 to 0.921 | g | c | 0.000 (0) | 0.894 | 0.871 to 0.918 | bc | 0.871 | 0.844 to 0.897 | 0.157 | 0.917 | 0.882 | 0.876 | +0.041 | -0.027 to 0.111 | 0.5247 | 0.233 | 0.190 to 0.283 | 10.21 | 0.107 | 0.083 to 0.132 |
| `gpt-oss-20b` | main | 0.814 | 0.773 to 0.855 | h | d | 0.026 (59) | 0.815 | 0.780 to 0.852 | d | 0.793 | 0.754 to 0.832 | 0.237 | 0.905 | 0.782 | 0.716 | +0.189 | 0.085 to 0.302 | 0.0234 | 0.148 | 0.115 to 0.182 | 3.77 | 0.210 | 0.182 to 0.240 |

**Tiers.** Of the 55 pairs on FAC, 32 differ at a Holm-adjusted p below 0.05 under the matched settings (33 under the defaults). Matched settings: a: `deepseek-v4.1-flash`, `claude-sonnet-5`, `kimi-k3`, `glm-5.3-flash`, `muse-glimmer-30b`; b: `claude-sonnet-5`, `kimi-k3`, `muse-glimmer-30b`, `glm-5.3`; c: `kimi-k3`, `glm-5.3-flash`, `muse-glimmer-30b`, `gpt-5.4-mini`; d: `kimi-k3`, `muse-glimmer-30b`, `gpt-5.4-mini`, `glm-5.3`; e: `gpt-5.4-mini`, `glm-5.3`, `gemma-4-26b-a4b`; f: `glm-5.3`, `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`; g: `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`, `qwen3-235b-a22b-2507`; h: `gpt-oss-20b`. Defaults: a: `deepseek-v4.1-flash`, `claude-sonnet-5`, `kimi-k3`, `glm-5.3-flash`, `muse-glimmer-30b`; b: `claude-sonnet-5`, `kimi-k3`, `muse-glimmer-30b`, `glm-5.3`; c: `gpt-5.4-mini`, `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`, `qwen3-235b-a22b-2507`; d: `gpt-5.4-mini`, `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`, `gpt-oss-20b`. On MC with the judge, 16 of 55 pairs differ (27 under the defaults); by matching alone 16 (24 under the defaults). Kendall's tau between the default and the matched orderings on FAC: 0.745 (95% CI 0.673 to 0.855, templates resampled, 11 models).

## Paired change per re-run model

Reasoning-on rows minus default rows, paired by item over every item both stores hold: the item mean, the template bootstrap at 95% and 90%, the sign-flip test over templates with Holm over the re-run models, the smallest change the design detects at 80% power (two-sided 0.05), and the 90% interval read against the ±0.05 margin of D-165 ("yes": the change is bounded inside it). MC changes are over the items with milestones where both responses have a coverage value.

| model | setting | items | FAC default | FAC reasoning on | change | 95% CI | p (Holm) | detectable | 90% CI | within ±0.05 | MC change with the judge (items) | 95% CI | p (Holm) | MC change by matching alone | 95% CI | empty, default / reasoning on | no readable answer, default / reasoning on |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---|
| `gemma-4-26b-a4b` | effort=medium | 2250 | 0.868 | 0.932 | +0.064 | 0.042 to 0.088 | 0.0003 | 0.032 | 0.046 to 0.084 | no | +0.008 (2180) | -0.003 to 0.020 | 0.1770 | -0.003 | -0.019 to 0.010 | 0 / 3 | 0 / 3 |
| `gpt-5.4-mini` | effort=medium | 2250 | 0.850 | 0.960 | +0.111 | 0.078 to 0.149 | 0.0003 | 0.053 | 0.082 to 0.142 | no | +0.057 (2180) | 0.035 to 0.082 | 0.0003 | +0.049 | 0.029 to 0.073 | 0 / 0 | 0 / 0 |
| `gemini-3.1-flash-lite` | effort=medium | 2250 | 0.874 | 0.910 | +0.036 | 0.020 to 0.053 | 0.0003 | 0.023 | 0.023 to 0.051 | no | +0.016 (2180) | 0.006 to 0.025 | 0.0014 | +0.016 | 0.006 to 0.026 | 0 / 0 | 0 / 0 |

## FAC pairs, matched settings

The sign-flip permutation test on the 11 models' per-template differences, Holm over the 55 pairs.

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.166 | -0.208 to -0.122 | 0.0055 | 0.060 / 0.090 |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.166 | -0.212 to -0.128 | 0.0055 | 0.060 / 0.089 |
| `gpt-oss-20b` | `kimi-k3` | -0.159 | -0.202 to -0.119 | 0.0055 | 0.060 / 0.089 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.158 | -0.201 to -0.120 | 0.0055 | 0.059 / 0.088 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.154 | -0.194 to -0.117 | 0.0055 | 0.056 / 0.083 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.146 | -0.187 to -0.105 | 0.0055 | 0.059 / 0.087 |
| `gpt-oss-20b` | `glm-5.3` | -0.134 | -0.180 to -0.091 | 0.0055 | 0.061 / 0.091 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.118 | -0.156 to -0.080 | 0.0055 | 0.056 / 0.083 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.095 | -0.141 to -0.047 | 0.0055 | 0.064 / 0.094 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.087 | 0.063 to 0.114 | 0.0055 | 0.038 / 0.057 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.086 | -0.113 to -0.062 | 0.0055 | 0.038 / 0.056 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.080 | -0.110 to -0.051 | 0.0055 | 0.040 / 0.060 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.078 | -0.104 to -0.054 | 0.0055 | 0.037 / 0.055 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.074 | -0.098 to -0.053 | 0.0055 | 0.033 / 0.050 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.071 | 0.045 to 0.100 | 0.0055 | 0.040 / 0.060 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.071 | -0.100 to -0.044 | 0.0055 | 0.041 / 0.061 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | -0.066 | -0.094 to -0.040 | 0.0055 | 0.038 / 0.057 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.064 | 0.037 to 0.093 | 0.0055 | 0.041 / 0.061 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.063 | 0.039 to 0.091 | 0.0055 | 0.038 / 0.057 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.059 | 0.035 to 0.085 | 0.0055 | 0.039 / 0.057 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.048 | -0.071 to -0.029 | 0.0055 | 0.031 / 0.046 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.048 | -0.073 to -0.028 | 0.0055 | 0.032 / 0.047 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.036 | -0.056 to -0.018 | 0.0066 | 0.027 / 0.041 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.041 | -0.068 to -0.020 | 0.0096 | 0.034 / 0.051 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.040 | -0.062 to -0.019 | 0.0124 | 0.031 / 0.047 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.080 | -0.125 to -0.036 | 0.0150 | 0.061 / 0.091 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | +0.050 | 0.022 to 0.082 | 0.0150 | 0.042 / 0.063 |
| `glm-5.3-flash` | `glm-5.3` | +0.024 | 0.011 to 0.041 | 0.0150 | 0.022 / 0.033 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.054 | -0.084 to -0.025 | 0.0189 | 0.043 / 0.063 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.021 | 0.008 to 0.033 | 0.0189 | 0.017 / 0.025 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.020 | -0.032 to -0.009 | 0.0250 | 0.017 / 0.026 |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.033 | 0.011 to 0.059 | 0.0480 | 0.032 / 0.047 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | +0.038 | 0.014 to 0.063 | 0.0575 | 0.035 / 0.052 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.038 | 0.012 to 0.068 | 0.1034 | 0.038 / 0.057 |
| `glm-5.3` | `claude-sonnet-5` | -0.033 | -0.060 to -0.007 | 0.2079 | 0.037 / 0.055 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | -0.028 | -0.051 to -0.008 | 0.3280 | 0.033 / 0.048 |
| `glm-5.3` | `kimi-k3` | -0.026 | -0.050 to -0.005 | 0.4579 | 0.033 / 0.048 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.023 | 0.005 to 0.044 | 0.4734 | 0.028 / 0.042 |
| `kimi-k3` | `gpt-5.4-mini` | +0.014 | 0.001 to 0.027 | 0.9587 | 0.020 / 0.029 |
| `glm-5.3` | `muse-glimmer-30b` | -0.021 | -0.045 to 0.001 | 1.0000 | 0.034 / 0.050 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | -0.016 | -0.040 to 0.013 | 1.0000 | 0.040 / 0.060 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.016 | -0.043 to 0.012 | 1.0000 | 0.039 / 0.058 |
| `glm-5.3` | `gpt-5.4-mini` | -0.012 | -0.042 to 0.013 | 1.0000 | 0.038 / 0.057 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.012 | -0.001 to 0.028 | 1.0000 | 0.021 / 0.030 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.012 | -0.008 to 0.029 | 1.0000 | 0.026 / 0.038 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.012 | -0.030 to 0.003 | 1.0000 | 0.024 / 0.036 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.008 | -0.002 to 0.021 | 1.0000 | 0.017 / 0.025 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.008 | -0.010 to 0.024 | 1.0000 | 0.025 / 0.036 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.008 | -0.025 to 0.007 | 1.0000 | 0.024 / 0.035 |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.007 | -0.003 to 0.019 | 1.0000 | 0.016 / 0.023 |
| `kimi-k3` | `claude-sonnet-5` | -0.007 | -0.019 to 0.004 | 1.0000 | 0.017 / 0.025 |
| `muse-glimmer-30b` | `kimi-k3` | -0.005 | -0.025 to 0.013 | 1.0000 | 0.028 / 0.041 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.004 | -0.010 to 0.020 | 1.0000 | 0.022 / 0.033 |
| `glm-5.3-flash` | `kimi-k3` | -0.001 | -0.017 to 0.015 | 1.0000 | 0.023 / 0.034 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.000 | -0.010 to 0.011 | 1.0000 | 0.014 / 0.021 |

## MC pairs with the judge, matched settings

As Q3's coverage pairs: E5-strict template means over 147 templates, sign flips, Holm over the pairs.

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.108 | -0.141 to -0.079 | 0.0055 | 0.043 / 0.064 |
| `gpt-oss-20b` | `kimi-k3` | -0.097 | -0.130 to -0.066 | 0.0055 | 0.044 / 0.066 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.097 | -0.131 to -0.066 | 0.0055 | 0.047 / 0.069 |
| `gpt-oss-20b` | `glm-5.3` | -0.086 | -0.123 to -0.051 | 0.0055 | 0.053 / 0.078 |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.085 | -0.114 to -0.056 | 0.0055 | 0.043 / 0.063 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.083 | -0.112 to -0.056 | 0.0055 | 0.040 / 0.060 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.079 | -0.111 to -0.051 | 0.0055 | 0.044 / 0.065 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.077 | -0.106 to -0.048 | 0.0055 | 0.041 / 0.061 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.069 | -0.096 to -0.043 | 0.0055 | 0.041 / 0.060 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.065 | -0.094 to -0.037 | 0.0055 | 0.043 / 0.064 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.044 | -0.064 to -0.023 | 0.0055 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.039 | -0.058 to -0.020 | 0.0055 | 0.027 / 0.040 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.031 | -0.045 to -0.019 | 0.0055 | 0.020 / 0.030 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.033 | 0.016 to 0.052 | 0.0210 | 0.026 / 0.039 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.028 | -0.045 to -0.011 | 0.0246 | 0.025 / 0.036 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.029 | -0.047 to -0.012 | 0.0280 | 0.026 / 0.038 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.023 | -0.040 to -0.009 | 0.1053 | 0.021 / 0.031 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.032 | 0.010 to 0.055 | 0.1748 | 0.031 / 0.047 |
| `kimi-k3` | `gpt-5.4-mini` | +0.020 | 0.006 to 0.036 | 0.2146 | 0.021 / 0.031 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.027 | -0.048 to -0.008 | 0.2988 | 0.029 / 0.043 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.026 | -0.046 to -0.008 | 0.2988 | 0.028 / 0.042 |
| `glm-5.3` | `claude-sonnet-5` | -0.022 | -0.052 to 0.005 | 1.0000 | 0.041 / 0.062 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.021 | -0.007 to 0.050 | 1.0000 | 0.041 / 0.061 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.020 | 0.001 to 0.040 | 1.0000 | 0.027 / 0.040 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.020 | -0.002 to 0.041 | 1.0000 | 0.030 / 0.045 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.018 | -0.036 to -0.000 | 1.0000 | 0.026 / 0.038 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.018 | -0.000 to 0.035 | 1.0000 | 0.027 / 0.040 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.018 | -0.038 to 0.003 | 1.0000 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.017 | -0.044 to 0.011 | 1.0000 | 0.040 / 0.059 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.016 | -0.033 to -0.000 | 1.0000 | 0.024 / 0.035 |
| `muse-glimmer-30b` | `kimi-k3` | -0.015 | -0.035 to 0.002 | 1.0000 | 0.027 / 0.041 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.014 | 0.001 to 0.030 | 1.0000 | 0.022 / 0.032 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.014 | -0.003 to 0.032 | 1.0000 | 0.024 / 0.036 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.013 | -0.030 to 0.004 | 1.0000 | 0.025 / 0.038 |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.013 | -0.025 to -0.000 | 1.0000 | 0.018 / 0.026 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | +0.012 | -0.004 to 0.029 | 1.0000 | 0.024 / 0.036 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.012 | -0.028 to 0.005 | 1.0000 | 0.026 / 0.038 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.012 | -0.034 to 0.007 | 1.0000 | 0.029 / 0.043 |
| `glm-5.3` | `kimi-k3` | -0.011 | -0.039 to 0.012 | 1.0000 | 0.037 / 0.055 |
| `kimi-k3` | `claude-sonnet-5` | -0.011 | -0.025 to 0.004 | 1.0000 | 0.020 / 0.030 |
| `glm-5.3-flash` | `glm-5.3` | +0.011 | -0.005 to 0.027 | 1.0000 | 0.023 / 0.034 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.010 | -0.023 to 0.004 | 1.0000 | 0.020 / 0.030 |
| `glm-5.3` | `gpt-5.4-mini` | +0.009 | -0.020 to 0.039 | 1.0000 | 0.043 / 0.064 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | -0.008 | -0.024 to 0.008 | 1.0000 | 0.023 / 0.034 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.008 | -0.005 to 0.023 | 1.0000 | 0.020 / 0.030 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.007 | -0.036 to 0.023 | 1.0000 | 0.042 / 0.063 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.006 | -0.011 to 0.024 | 1.0000 | 0.027 / 0.040 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.005 | -0.014 to 0.022 | 1.0000 | 0.027 / 0.040 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.004 | -0.007 to 0.016 | 1.0000 | 0.016 / 0.024 |
| `glm-5.3` | `muse-glimmer-30b` | +0.004 | -0.020 to 0.027 | 1.0000 | 0.032 / 0.047 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.003 | -0.020 to 0.014 | 1.0000 | 0.025 / 0.037 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.002 | -0.014 to 0.022 | 1.0000 | 0.026 / 0.039 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.002 | -0.014 to 0.019 | 1.0000 | 0.024 / 0.036 |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.001 | -0.026 to 0.027 | 1.0000 | 0.038 / 0.056 |
| `glm-5.3-flash` | `kimi-k3` | -0.001 | -0.020 to 0.016 | 1.0000 | 0.026 / 0.039 |

## MC pairs by matching alone, default configuration

E3 coverage (matching alone, the judge left out) on the default configuration, template means over 147 templates, sign flips, Holm over the pairs: the family Q3 runs on E5-strict. Per model: `claude-sonnet-5` 0.903 (0.881 to 0.924); `glm-5.3-flash` 0.882 (0.849 to 0.911); `glm-5.3` 0.876 (0.842 to 0.908); `kimi-k3` 0.876 (0.849 to 0.900); `qwen3-235b-a22b-2507` 0.871 (0.844 to 0.897); `deepseek-v4.1-flash` 0.862 (0.834 to 0.889); `muse-glimmer-30b` 0.860 (0.831 to 0.889); `gemma-4-26b-a4b` 0.855 (0.825 to 0.884); `gemini-3.1-flash-lite` 0.840 (0.808 to 0.870); `gpt-5.4-mini` 0.798 (0.764 to 0.837); `gpt-oss-20b` 0.793 (0.754 to 0.832).

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.110 | -0.143 to -0.082 | 0.0055 | 0.042 / 0.062 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.105 | -0.134 to -0.076 | 0.0055 | 0.041 / 0.061 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.090 | -0.124 to -0.057 | 0.0055 | 0.048 / 0.072 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.084 | 0.057 to 0.117 | 0.0055 | 0.045 / 0.067 |
| `gpt-oss-20b` | `glm-5.3` | -0.084 | -0.119 to -0.048 | 0.0055 | 0.053 / 0.079 |
| `gpt-oss-20b` | `kimi-k3` | -0.083 | -0.115 to -0.050 | 0.0055 | 0.045 / 0.067 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.079 | -0.109 to -0.049 | 0.0055 | 0.045 / 0.067 |
| `glm-5.3` | `gpt-5.4-mini` | +0.078 | 0.045 to 0.113 | 0.0055 | 0.049 / 0.073 |
| `kimi-k3` | `gpt-5.4-mini` | +0.078 | 0.049 to 0.107 | 0.0055 | 0.041 / 0.060 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.073 | 0.047 to 0.102 | 0.0055 | 0.040 / 0.060 |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.070 | -0.100 to -0.041 | 0.0055 | 0.044 / 0.066 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.067 | -0.095 to -0.043 | 0.0055 | 0.039 / 0.058 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.064 | 0.040 to 0.094 | 0.0055 | 0.039 / 0.058 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.063 | -0.085 to -0.039 | 0.0055 | 0.032 / 0.048 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.063 | -0.090 to -0.036 | 0.0055 | 0.039 / 0.058 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.061 | 0.035 to 0.089 | 0.0055 | 0.038 / 0.057 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.057 | 0.032 to 0.084 | 0.0055 | 0.036 / 0.054 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.048 | -0.069 to -0.026 | 0.0055 | 0.031 / 0.046 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.040 | -0.059 to -0.023 | 0.0055 | 0.026 / 0.038 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.043 | -0.064 to -0.023 | 0.0072 | 0.031 / 0.045 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.042 | -0.066 to -0.019 | 0.0105 | 0.034 / 0.050 |
| `kimi-k3` | `claude-sonnet-5` | -0.027 | -0.042 to -0.011 | 0.0105 | 0.022 / 0.033 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.042 | 0.016 to 0.069 | 0.0297 | 0.035 / 0.052 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.036 | 0.015 to 0.057 | 0.0448 | 0.031 / 0.046 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.047 | -0.076 to -0.019 | 0.0527 | 0.043 / 0.064 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.031 | -0.053 to -0.012 | 0.0660 | 0.029 / 0.043 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.031 | 0.012 to 0.051 | 0.0986 | 0.028 / 0.041 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.036 | 0.010 to 0.064 | 0.3640 | 0.040 / 0.060 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.027 | -0.050 to -0.003 | 0.8207 | 0.035 / 0.052 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.023 | 0.003 to 0.044 | 0.8735 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.015 | 0.001 to 0.029 | 0.9424 | 0.020 / 0.030 |
| `glm-5.3` | `claude-sonnet-5` | -0.027 | -0.053 to 0.000 | 1.0000 | 0.040 / 0.060 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.022 | 0.000 to 0.047 | 1.0000 | 0.032 / 0.047 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.021 | -0.047 to 0.009 | 1.0000 | 0.042 / 0.063 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.021 | -0.043 to 0.001 | 1.0000 | 0.031 / 0.046 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.020 | -0.043 to -0.001 | 1.0000 | 0.030 / 0.044 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.020 | -0.042 to 0.001 | 1.0000 | 0.032 / 0.047 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.020 | 0.001 to 0.040 | 1.0000 | 0.030 / 0.045 |
| `glm-5.3` | `muse-glimmer-30b` | +0.017 | -0.007 to 0.042 | 1.0000 | 0.034 / 0.051 |
| `muse-glimmer-30b` | `kimi-k3` | -0.016 | -0.038 to 0.003 | 1.0000 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.016 | -0.032 to 0.001 | 1.0000 | 0.025 / 0.038 |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.014 | -0.041 to 0.017 | 1.0000 | 0.042 / 0.062 |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.013 | -0.027 to 0.000 | 1.0000 | 0.020 / 0.029 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | +0.012 | -0.010 to 0.032 | 1.0000 | 0.031 / 0.047 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.011 | -0.034 to 0.012 | 1.0000 | 0.033 / 0.049 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | -0.009 | -0.032 to 0.014 | 1.0000 | 0.034 / 0.050 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.007 | -0.030 to 0.014 | 1.0000 | 0.032 / 0.047 |
| `glm-5.3-flash` | `kimi-k3` | +0.006 | -0.014 to 0.026 | 1.0000 | 0.029 / 0.043 |
| `glm-5.3-flash` | `glm-5.3` | +0.006 | -0.010 to 0.023 | 1.0000 | 0.024 / 0.035 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.005 | -0.036 to 0.027 | 1.0000 | 0.047 / 0.069 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.005 | -0.034 to 0.027 | 1.0000 | 0.044 / 0.065 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.004 | -0.026 to 0.016 | 1.0000 | 0.029 / 0.043 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.004 | -0.025 to 0.017 | 1.0000 | 0.030 / 0.045 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.003 | -0.018 to 0.024 | 1.0000 | 0.032 / 0.047 |
| `glm-5.3` | `kimi-k3` | +0.000 | -0.027 to 0.025 | 1.0000 | 0.038 / 0.056 |

## Serving endpoints of the reasoning-on stores

As the default configuration's endpoint table (`providers.json`), for the re-run models served by more than one endpoint.

| model | store | endpoint | rows | raw score | unusable | matched difference | templates matched | templates served |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `gemma-4-26b-a4b` | reasoning-medium-full | Io Net | 2194 | 0.934 | 0.001 | -0.010 | 36 | 150 |
| `gemma-4-26b-a4b` | reasoning-medium-full | DekaLLM | 56 | 0.884 | 0.018 | +0.010 | 36 | 36 |
| `gemini-3.1-flash-lite` | reasoning-medium-full | Google AI Studio | 2132 | 0.906 | 0.000 | +0.001 | 18 | 150 |
| `gemini-3.1-flash-lite` | reasoning-medium-full | Google | 118 | 0.975 | 0.000 | -0.001 | 18 | 18 |
