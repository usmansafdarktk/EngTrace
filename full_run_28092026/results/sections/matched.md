# Matched settings

Printed by `analyze.py` (WS-C1) from `results/matched_config.json` (written 2026-10-07T13:20:31Z, sha256 `e1e4d06b307e`): each model at its reasoning-on store where the file names one, else at its default store. With reasoning on: `gpt-5.4-mini` (effort=medium, `scores/reasoning-medium-full`), `gemma-4-26b-a4b` (effort=medium, `scores/reasoning-medium-full`), `gemini-3.1-flash-lite` (effort=medium, `scores/reasoning-medium-full`); the others at `scores/main`. Skipped: `qwen3-235b-a22b-thinking-2507` (run: false). The functions and seeds are those of Q1, Q2 and Q3, so a model whose rows did not change shows the numbers of the default configuration. FAC is the answer score (correct 1, partial 0.5, incorrect or unusable 0), the mean of the 150 template means; every interval is 95% and resamples templates. The letters are the compact letter display of the sign-flip pairs at a Holm-adjusted p below 0.05: models that share a letter are not separated. No readable answer: empty, or answered without an answer the check can read. MC is Milestone Coverage, the mean of the template means over the templates every model has coverage on (148 with the judge, 148 by matching alone): with the judge (E5-strict) and by matching alone (E3). Judge-decided share: the milestones sent to the judge over those required, on answered responses. The flag rates are on fully solved responses: the arithmetic check's flags and the step check's judge flags, template intervals.

| model | store | FAC | 95% CI | letter | letter, default | no readable answer | MC with the judge | 95% CI | letter | MC by matching alone | 95% CI | judge-decided share | Easy | Intermediate | Advanced | gap | 95% CI | Welch p (Holm) | arithmetic flags | 95% CI | calculations parsed per response | judged step flags | 95% CI |
|---|---|---:|---:|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | main | 0.988 | 0.979 to 0.995 | a | a | 0.004 (8) | 0.899 | 0.874 to 0.922 | abc | 0.861 | 0.831 to 0.890 | 0.169 | 0.995 | 0.986 | 0.979 | +0.016 | -0.001 to 0.041 | 0.7911 | 0.006 | 0.003 to 0.010 | 3.27 | 0.010 | 0.005 to 0.015 |
| `claude-sonnet-5` | main | 0.987 | 0.980 to 0.993 | ab | ab | 0.001 (2) | 0.922 | 0.902 to 0.942 | a | 0.903 | 0.880 to 0.925 | 0.125 | 0.993 | 0.983 | 0.984 | +0.008 | -0.007 to 0.024 | 1.0000 | 0.143 | 0.111 to 0.177 | 5.32 | 0.054 | 0.041 to 0.068 |
| `kimi-k3` | main | 0.980 | 0.967 to 0.990 | abc | ab | 0.007 (16) | 0.912 | 0.888 to 0.934 | ab | 0.877 | 0.850 to 0.902 | 0.150 | 0.992 | 0.985 | 0.953 | +0.039 | -0.002 to 0.094 | 0.7911 | 0.007 | 0.003 to 0.011 | 2.79 | 0.017 | 0.011 to 0.024 |
| `glm-5.3-flash` | main | 0.972 | 0.953 to 0.987 | bc | b | 0.020 (45) | 0.912 | 0.886 to 0.936 | abc | 0.891 | 0.864 to 0.916 | 0.126 | 0.994 | 0.970 | 0.939 | +0.055 | 0.012 to 0.113 | 0.3679 | 0.019 | 0.012 to 0.026 | 4.74 | 0.032 | 0.023 to 0.042 |
| `gpt-5.4-mini` (reasoning on) | reasoning-medium-full | 0.969 | 0.958 to 0.980 | cde | de | 0.000 (0) | 0.892 | 0.866 to 0.917 | bc | 0.847 | 0.815 to 0.877 | 0.184 | 0.982 | 0.953 | 0.975 | +0.006 | -0.013 to 0.028 | 1.0000 | 0.075 | 0.057 to 0.096 | 3.27 | 0.115 | 0.091 to 0.140 |
| `muse-glimmer-30b` | main | 0.969 | 0.949 to 0.985 | abcd | abc | 0.014 (31) | 0.897 | 0.871 to 0.922 | abc | 0.862 | 0.832 to 0.890 | 0.167 | 0.976 | 0.968 | 0.958 | +0.019 | -0.026 to 0.067 | 1.0000 | 0.058 | 0.044 to 0.074 | 3.27 | 0.059 | 0.046 to 0.074 |
| `glm-5.3` | main | 0.950 | 0.921 to 0.975 | def | c | 0.044 (99) | 0.900 | 0.866 to 0.930 | abc | 0.881 | 0.846 to 0.912 | 0.107 | 0.984 | 0.966 | 0.865 | +0.120 | 0.030 to 0.223 | 0.2164 | 0.018 | 0.012 to 0.025 | 5.72 | 0.026 | 0.016 to 0.036 |
| `gemma-4-26b-a4b` (reasoning on) | reasoning-medium-full | 0.942 | 0.916 to 0.964 | ef | d | 0.001 (3) | 0.884 | 0.856 to 0.910 | bc | 0.850 | 0.819 to 0.880 | 0.186 | 0.967 | 0.934 | 0.912 | +0.055 | -0.003 to 0.116 | 0.5284 | 0.099 | 0.066 to 0.135 | 3.30 | 0.074 | 0.053 to 0.097 |
| `gemini-3.1-flash-lite` (reasoning on) | reasoning-medium-full | 0.917 | 0.884 to 0.945 | fg | d | 0.000 (0) | 0.881 | 0.851 to 0.908 | c | 0.855 | 0.824 to 0.886 | 0.187 | 0.959 | 0.928 | 0.826 | +0.132 | 0.049 to 0.225 | 0.0617 | 0.069 | 0.047 to 0.094 | 3.54 | 0.093 | 0.068 to 0.121 |
| `qwen3-235b-a22b-2507` | main | 0.899 | 0.869 to 0.926 | g | d | 0.000 (0) | 0.892 | 0.867 to 0.916 | bc | 0.866 | 0.836 to 0.894 | 0.158 | 0.917 | 0.890 | 0.885 | +0.031 | -0.035 to 0.097 | 1.0000 | 0.232 | 0.187 to 0.279 | 10.21 | 0.107 | 0.084 to 0.131 |
| `gpt-oss-20b` | main | 0.808 | 0.764 to 0.849 | h | e | 0.026 (59) | 0.812 | 0.775 to 0.847 | d | 0.790 | 0.752 to 0.826 | 0.238 | 0.900 | 0.768 | 0.719 | +0.181 | 0.070 to 0.299 | 0.0353 | 0.147 | 0.116 to 0.181 | 3.77 | 0.206 | 0.177 to 0.237 |

**Tiers.** Of the 55 pairs on FAC, 36 differ at a Holm-adjusted p below 0.05 under the matched settings (38 under the defaults). Matched settings: a: `deepseek-v4.1-flash`, `claude-sonnet-5`, `kimi-k3`, `muse-glimmer-30b`; b: `claude-sonnet-5`, `kimi-k3`, `glm-5.3-flash`, `muse-glimmer-30b`; c: `kimi-k3`, `glm-5.3-flash`, `gpt-5.4-mini`, `muse-glimmer-30b`; d: `gpt-5.4-mini`, `muse-glimmer-30b`, `glm-5.3`; e: `gpt-5.4-mini`, `glm-5.3`, `gemma-4-26b-a4b`; f: `glm-5.3`, `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`; g: `gemini-3.1-flash-lite`, `qwen3-235b-a22b-2507`; h: `gpt-oss-20b`. Defaults: a: `deepseek-v4.1-flash`, `claude-sonnet-5`, `kimi-k3`, `muse-glimmer-30b`; b: `claude-sonnet-5`, `kimi-k3`, `glm-5.3-flash`, `muse-glimmer-30b`; c: `muse-glimmer-30b`, `glm-5.3`; d: `gpt-5.4-mini`, `gemma-4-26b-a4b`, `gemini-3.1-flash-lite`, `qwen3-235b-a22b-2507`; e: `gpt-5.4-mini`, `gpt-oss-20b`. On MC with the judge, 15 of 55 pairs differ (28 under the defaults); by matching alone 20 (26 under the defaults). Kendall's tau between the default and the matched orderings on FAC: 0.709 (95% CI 0.636 to 0.818, templates resampled, 11 models).

## Paired change per re-run model

Reasoning-on rows minus default rows, paired by item over every item both stores hold: the item mean, the template bootstrap at 95% and 90%, the sign-flip test over templates with Holm over the re-run models, the smallest change the design detects at 80% power (two-sided 0.05), and the 90% interval read against the ±0.05 margin of D-165 ("yes": the change is bounded inside it). MC changes are over the items with milestones where both responses have a coverage value.

| model | setting | items | FAC default | FAC reasoning on | change | 95% CI | p (Holm) | detectable | 90% CI | within ±0.05 | MC change with the judge (items) | 95% CI | p (Holm) | MC change by matching alone | 95% CI | empty, default / reasoning on | no readable answer, default / reasoning on |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---|
| `gemma-4-26b-a4b` | effort=medium | 2250 | 0.872 | 0.942 | +0.070 | 0.049 to 0.092 | <0.0001 | 0.030 | 0.053 to 0.088 | no | +0.008 (2193) | -0.004 to 0.021 | 0.1786 | -0.005 | -0.020 to 0.009 | 0 / 3 | 0 / 3 |
| `gpt-5.4-mini` | effort=medium | 2250 | 0.852 | 0.969 | +0.118 | 0.082 to 0.156 | <0.0001 | 0.053 | 0.088 to 0.149 | no | +0.058 (2193) | 0.036 to 0.082 | <0.0001 | +0.050 | 0.029 to 0.072 | 0 / 0 | 0 / 0 |
| `gemini-3.1-flash-lite` | effort=medium | 2250 | 0.878 | 0.917 | +0.038 | 0.023 to 0.056 | <0.0001 | 0.024 | 0.025 to 0.053 | no | +0.016 (2193) | 0.007 to 0.025 | 0.0010 | +0.016 | 0.006 to 0.025 | 0 / 0 | 0 / 0 |

## FAC pairs, matched settings

The sign-flip permutation test on the 11 models' per-template differences, Holm over the 55 pairs.

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.180 | -0.223 to -0.140 | 0.0005 | 0.060 / 0.089 |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.179 | -0.222 to -0.140 | 0.0005 | 0.059 / 0.088 |
| `gpt-oss-20b` | `kimi-k3` | -0.172 | -0.215 to -0.133 | 0.0005 | 0.059 / 0.088 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.164 | -0.207 to -0.124 | 0.0005 | 0.059 / 0.088 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.161 | -0.201 to -0.123 | 0.0005 | 0.057 / 0.085 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.161 | -0.200 to -0.124 | 0.0005 | 0.055 / 0.081 |
| `gpt-oss-20b` | `glm-5.3` | -0.142 | -0.186 to -0.101 | 0.0005 | 0.061 / 0.090 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.134 | -0.173 to -0.097 | 0.0005 | 0.054 / 0.081 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.109 | -0.152 to -0.065 | 0.0005 | 0.062 / 0.093 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.089 | 0.064 to 0.118 | 0.0005 | 0.038 / 0.057 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.088 | -0.114 to -0.064 | 0.0005 | 0.037 / 0.055 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.081 | -0.110 to -0.054 | 0.0005 | 0.040 / 0.060 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.073 | -0.100 to -0.048 | 0.0005 | 0.037 / 0.055 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.072 | 0.044 to 0.101 | 0.0005 | 0.041 / 0.060 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.070 | -0.102 to -0.044 | 0.0005 | 0.041 / 0.061 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | -0.070 | -0.098 to -0.045 | 0.0005 | 0.038 / 0.056 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.070 | -0.093 to -0.048 | 0.0005 | 0.033 / 0.049 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.064 | 0.037 to 0.094 | 0.0005 | 0.041 / 0.061 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.056 | 0.032 to 0.083 | 0.0005 | 0.037 / 0.054 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.046 | -0.069 to -0.026 | 0.0005 | 0.030 / 0.045 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.052 | 0.028 to 0.080 | 0.0007 | 0.038 / 0.056 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.045 | -0.068 to -0.025 | 0.0007 | 0.031 / 0.046 |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.038 | 0.018 to 0.062 | 0.0007 | 0.031 / 0.047 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.091 | -0.134 to -0.049 | 0.0013 | 0.061 / 0.091 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | +0.053 | 0.025 to 0.083 | 0.0043 | 0.041 / 0.061 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | +0.043 | 0.018 to 0.068 | 0.0141 | 0.035 / 0.052 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.019 | 0.008 to 0.030 | 0.0141 | 0.016 / 0.023 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.051 | -0.081 to -0.021 | 0.0246 | 0.043 / 0.063 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.018 | -0.029 to -0.008 | 0.0256 | 0.015 / 0.022 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.038 | -0.063 to -0.016 | 0.0257 | 0.034 / 0.051 |
| `glm-5.3` | `claude-sonnet-5` | -0.037 | -0.065 to -0.014 | 0.0495 | 0.037 / 0.055 |
| `glm-5.3` | `kimi-k3` | -0.030 | -0.055 to -0.010 | 0.0495 | 0.031 / 0.047 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.027 | -0.046 to -0.010 | 0.0495 | 0.026 / 0.038 |
| `glm-5.3-flash` | `glm-5.3` | +0.022 | 0.008 to 0.039 | 0.0495 | 0.022 / 0.033 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.016 | 0.006 to 0.028 | 0.0495 | 0.016 / 0.023 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.030 | -0.052 to -0.011 | 0.0500 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.025 | 0.007 to 0.045 | 0.1718 | 0.027 / 0.041 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.019 | 0.005 to 0.035 | 0.1820 | 0.022 / 0.033 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.033 | 0.009 to 0.061 | 0.1945 | 0.037 / 0.056 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | -0.027 | -0.050 to -0.006 | 0.2805 | 0.032 / 0.048 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.018 | -0.036 to -0.002 | 0.5979 | 0.025 / 0.037 |
| `glm-5.3` | `gpt-5.4-mini` | -0.019 | -0.048 to 0.004 | 1.0000 | 0.037 / 0.055 |
| `glm-5.3` | `muse-glimmer-30b` | -0.019 | -0.044 to 0.004 | 1.0000 | 0.035 / 0.053 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | -0.018 | -0.044 to 0.012 | 1.0000 | 0.039 / 0.059 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.015 | -0.033 to -0.001 | 1.0000 | 0.023 / 0.034 |
| `muse-glimmer-30b` | `kimi-k3` | -0.011 | -0.032 to 0.008 | 1.0000 | 0.029 / 0.043 |
| `kimi-k3` | `gpt-5.4-mini` | +0.011 | -0.003 to 0.025 | 1.0000 | 0.020 / 0.029 |
| `glm-5.3-flash` | `kimi-k3` | -0.008 | -0.024 to 0.006 | 1.0000 | 0.021 / 0.031 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.008 | -0.034 to 0.018 | 1.0000 | 0.038 / 0.056 |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.008 | -0.002 to 0.020 | 1.0000 | 0.015 / 0.023 |
| `kimi-k3` | `claude-sonnet-5` | -0.007 | -0.020 to 0.004 | 1.0000 | 0.017 / 0.026 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.003 | -0.014 to 0.021 | 1.0000 | 0.025 / 0.036 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.003 | -0.015 to 0.018 | 1.0000 | 0.024 / 0.035 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.001 | -0.008 to 0.010 | 1.0000 | 0.013 / 0.019 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | -0.000 | -0.018 to 0.016 | 1.0000 | 0.024 / 0.035 |

## MC pairs with the judge, matched settings

As Q3's coverage pairs: E5-strict template means over 148 templates, sign flips, Holm over the pairs.

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.110 | -0.141 to -0.080 | 0.0005 | 0.043 / 0.064 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.100 | -0.132 to -0.068 | 0.0005 | 0.047 / 0.069 |
| `gpt-oss-20b` | `kimi-k3` | -0.100 | -0.130 to -0.070 | 0.0005 | 0.045 / 0.066 |
| `gpt-oss-20b` | `glm-5.3` | -0.088 | -0.126 to -0.050 | 0.0005 | 0.053 / 0.078 |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.087 | -0.117 to -0.057 | 0.0005 | 0.043 / 0.064 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.085 | -0.114 to -0.057 | 0.0005 | 0.041 / 0.061 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.080 | -0.110 to -0.051 | 0.0005 | 0.042 / 0.062 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.080 | -0.111 to -0.050 | 0.0005 | 0.044 / 0.065 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.072 | -0.101 to -0.044 | 0.0005 | 0.041 / 0.061 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.068 | -0.099 to -0.039 | 0.0009 | 0.043 / 0.064 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.030 | -0.044 to -0.016 | 0.0013 | 0.020 / 0.029 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.042 | -0.063 to -0.021 | 0.0026 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.038 | -0.057 to -0.020 | 0.0030 | 0.027 / 0.040 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.030 | -0.049 to -0.013 | 0.0323 | 0.025 / 0.038 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.031 | 0.013 to 0.050 | 0.0381 | 0.027 / 0.040 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.028 | -0.045 to -0.011 | 0.0548 | 0.025 / 0.037 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.023 | -0.038 to -0.009 | 0.0686 | 0.021 / 0.031 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.031 | 0.010 to 0.053 | 0.1592 | 0.031 / 0.046 |
| `kimi-k3` | `gpt-5.4-mini` | +0.020 | 0.005 to 0.034 | 0.2431 | 0.021 / 0.030 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.028 | -0.047 to -0.008 | 0.2617 | 0.029 / 0.043 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.025 | -0.045 to -0.007 | 0.3759 | 0.028 / 0.041 |
| `glm-5.3` | `claude-sonnet-5` | -0.022 | -0.052 to 0.005 | 1.0000 | 0.041 / 0.061 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.020 | -0.041 to 0.002 | 1.0000 | 0.030 / 0.045 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.020 | -0.002 to 0.041 | 1.0000 | 0.030 / 0.045 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.020 | -0.038 to -0.002 | 1.0000 | 0.026 / 0.038 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.019 | -0.009 to 0.047 | 1.0000 | 0.041 / 0.061 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.018 | 0.000 to 0.038 | 1.0000 | 0.027 / 0.040 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.017 | -0.002 to 0.035 | 1.0000 | 0.027 / 0.040 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.016 | -0.043 to 0.012 | 1.0000 | 0.040 / 0.059 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.015 | -0.031 to 0.001 | 1.0000 | 0.024 / 0.035 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.015 | -0.002 to 0.031 | 1.0000 | 0.024 / 0.035 |
| `muse-glimmer-30b` | `kimi-k3` | -0.014 | -0.035 to 0.004 | 1.0000 | 0.027 / 0.041 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.013 | -0.030 to 0.005 | 1.0000 | 0.025 / 0.037 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.013 | -0.030 to 0.006 | 1.0000 | 0.026 / 0.038 |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.013 | -0.025 to -0.001 | 1.0000 | 0.017 / 0.026 |
| `glm-5.3-flash` | `glm-5.3` | +0.012 | -0.003 to 0.029 | 1.0000 | 0.023 / 0.034 |
| `glm-5.3` | `kimi-k3` | -0.012 | -0.039 to 0.012 | 1.0000 | 0.037 / 0.055 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | +0.012 | -0.005 to 0.029 | 1.0000 | 0.024 / 0.036 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.012 | -0.003 to 0.027 | 1.0000 | 0.022 / 0.033 |
| `kimi-k3` | `claude-sonnet-5` | -0.010 | -0.025 to 0.003 | 1.0000 | 0.020 / 0.030 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.010 | -0.031 to 0.009 | 1.0000 | 0.029 / 0.043 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | -0.008 | -0.024 to 0.008 | 1.0000 | 0.023 / 0.034 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.008 | -0.022 to 0.006 | 1.0000 | 0.020 / 0.030 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.008 | -0.035 to 0.022 | 1.0000 | 0.042 / 0.063 |
| `glm-5.3` | `gpt-5.4-mini` | +0.008 | -0.024 to 0.037 | 1.0000 | 0.043 / 0.063 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.007 | -0.011 to 0.026 | 1.0000 | 0.027 / 0.040 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.007 | -0.006 to 0.021 | 1.0000 | 0.020 / 0.030 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.005 | -0.022 to 0.013 | 1.0000 | 0.025 / 0.037 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.005 | -0.014 to 0.023 | 1.0000 | 0.027 / 0.040 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.004 | -0.008 to 0.015 | 1.0000 | 0.016 / 0.024 |
| `glm-5.3` | `muse-glimmer-30b` | +0.003 | -0.021 to 0.024 | 1.0000 | 0.032 / 0.047 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.002 | -0.015 to 0.021 | 1.0000 | 0.026 / 0.039 |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.001 | -0.027 to 0.026 | 1.0000 | 0.038 / 0.056 |
| `glm-5.3-flash` | `kimi-k3` | +0.000 | -0.019 to 0.018 | 1.0000 | 0.027 / 0.039 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | -0.000 | -0.017 to 0.016 | 1.0000 | 0.024 / 0.035 |

## MC pairs by matching alone, default configuration

E3 coverage (matching alone, the judge left out) on the default configuration, template means over 148 templates, sign flips, Holm over the pairs: the family Q3 runs on E5-strict. Per model: `claude-sonnet-5` 0.903 (0.880 to 0.925); `glm-5.3-flash` 0.891 (0.864 to 0.916); `glm-5.3` 0.881 (0.846 to 0.912); `kimi-k3` 0.877 (0.850 to 0.902); `qwen3-235b-a22b-2507` 0.866 (0.836 to 0.894); `muse-glimmer-30b` 0.862 (0.832 to 0.890); `deepseek-v4.1-flash` 0.861 (0.831 to 0.890); `gemma-4-26b-a4b` 0.855 (0.825 to 0.884); `gemini-3.1-flash-lite` 0.840 (0.808 to 0.870); `gpt-5.4-mini` 0.797 (0.758 to 0.834); `gpt-oss-20b` 0.790 (0.752 to 0.826).

| a | b | a - b | 95% CI | p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.113 | -0.143 to -0.084 | 0.0005 | 0.042 / 0.063 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.106 | -0.135 to -0.078 | 0.0005 | 0.040 / 0.060 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.102 | -0.135 to -0.069 | 0.0005 | 0.047 / 0.070 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.095 | 0.065 to 0.127 | 0.0005 | 0.044 / 0.065 |
| `gpt-oss-20b` | `glm-5.3` | -0.092 | -0.129 to -0.055 | 0.0005 | 0.052 / 0.078 |
| `gpt-oss-20b` | `kimi-k3` | -0.087 | -0.119 to -0.057 | 0.0005 | 0.045 / 0.067 |
| `glm-5.3` | `gpt-5.4-mini` | +0.084 | 0.052 to 0.119 | 0.0005 | 0.048 / 0.071 |
| `kimi-k3` | `gpt-5.4-mini` | +0.080 | 0.053 to 0.109 | 0.0005 | 0.040 / 0.059 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.076 | -0.111 to -0.043 | 0.0005 | 0.049 / 0.072 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.072 | -0.100 to -0.046 | 0.0005 | 0.039 / 0.058 |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.071 | -0.103 to -0.041 | 0.0005 | 0.045 / 0.067 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.066 | -0.094 to -0.039 | 0.0005 | 0.039 / 0.058 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.064 | 0.039 to 0.092 | 0.0005 | 0.038 / 0.057 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.063 | -0.086 to -0.041 | 0.0005 | 0.032 / 0.048 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.059 | 0.035 to 0.086 | 0.0005 | 0.036 / 0.054 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.042 | -0.060 to -0.024 | 0.0005 | 0.026 / 0.038 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.069 | 0.043 to 0.097 | 0.0008 | 0.038 / 0.057 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.065 | 0.039 to 0.092 | 0.0008 | 0.038 / 0.056 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.051 | 0.029 to 0.075 | 0.0008 | 0.033 / 0.049 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.047 | -0.069 to -0.027 | 0.0008 | 0.030 / 0.045 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.041 | -0.062 to -0.021 | 0.0031 | 0.030 / 0.045 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.043 | -0.068 to -0.020 | 0.0082 | 0.034 / 0.050 |
| `kimi-k3` | `claude-sonnet-5` | -0.026 | -0.042 to -0.011 | 0.0290 | 0.022 / 0.033 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.050 | -0.080 to -0.021 | 0.0352 | 0.043 / 0.064 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.037 | 0.015 to 0.059 | 0.0356 | 0.031 / 0.046 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.036 | -0.061 to -0.015 | 0.0356 | 0.033 / 0.049 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.036 | -0.060 to -0.013 | 0.0536 | 0.033 / 0.049 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.030 | 0.010 to 0.050 | 0.0832 | 0.028 / 0.042 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.041 | 0.014 to 0.069 | 0.0910 | 0.039 / 0.058 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.030 | -0.052 to -0.008 | 0.1576 | 0.031 / 0.047 |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.016 | -0.029 to -0.003 | 0.5512 | 0.019 / 0.029 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.026 | 0.002 to 0.048 | 0.5935 | 0.033 / 0.049 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.015 | 0.001 to 0.030 | 0.8515 | 0.020 / 0.030 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.022 | 0.001 to 0.043 | 0.9082 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.026 | -0.055 to 0.003 | 1.0000 | 0.041 / 0.061 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.025 | -0.051 to 0.000 | 1.0000 | 0.036 / 0.054 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.021 | -0.044 to -0.001 | 1.0000 | 0.031 / 0.046 |
| `glm-5.3` | `claude-sonnet-5` | -0.021 | -0.050 to 0.004 | 1.0000 | 0.039 / 0.058 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.021 | -0.002 to 0.044 | 1.0000 | 0.032 / 0.048 |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.020 | -0.048 to 0.010 | 1.0000 | 0.041 / 0.061 |
| `glm-5.3` | `muse-glimmer-30b` | +0.019 | -0.004 to 0.043 | 1.0000 | 0.033 / 0.050 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.015 | -0.047 to 0.018 | 1.0000 | 0.046 / 0.069 |
| `muse-glimmer-30b` | `kimi-k3` | -0.015 | -0.038 to 0.005 | 1.0000 | 0.030 / 0.044 |
| `glm-5.3-flash` | `kimi-k3` | +0.015 | -0.006 to 0.034 | 1.0000 | 0.029 / 0.042 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.011 | -0.033 to 0.008 | 1.0000 | 0.029 / 0.043 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.011 | -0.030 to 0.011 | 1.0000 | 0.029 / 0.043 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.011 | -0.035 to 0.012 | 1.0000 | 0.034 / 0.050 |
| `glm-5.3-flash` | `glm-5.3` | +0.010 | -0.005 to 0.026 | 1.0000 | 0.023 / 0.033 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.007 | -0.040 to 0.027 | 1.0000 | 0.047 / 0.070 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.007 | -0.027 to 0.014 | 1.0000 | 0.030 / 0.044 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.006 | -0.029 to 0.016 | 1.0000 | 0.032 / 0.048 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | -0.005 | -0.031 to 0.021 | 1.0000 | 0.037 / 0.056 |
| `glm-5.3` | `kimi-k3` | +0.004 | -0.023 to 0.029 | 1.0000 | 0.037 / 0.056 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | +0.004 | -0.022 to 0.029 | 1.0000 | 0.036 / 0.053 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | -0.001 | -0.021 to 0.022 | 1.0000 | 0.031 / 0.047 |

## Serving endpoints of the reasoning-on stores

As the default configuration's endpoint table (`providers.json`), for the re-run models served by more than one endpoint.

| model | store | endpoint | rows | raw score | unusable | matched difference | templates matched | templates served |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `gemma-4-26b-a4b` | reasoning-medium-full | Io Net | 2194 | 0.943 | 0.001 | -0.003 | 36 | 150 |
| `gemma-4-26b-a4b` | reasoning-medium-full | DekaLLM | 56 | 0.911 | 0.018 | +0.003 | 36 | 36 |
| `gemini-3.1-flash-lite` | reasoning-medium-full | Google AI Studio | 2132 | 0.913 | 0.000 | +0.001 | 18 | 150 |
| `gemini-3.1-flash-lite` | reasoning-medium-full | Google | 118 | 0.975 | 0.000 | -0.001 | 18 | 18 |
