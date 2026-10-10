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

## Pairwise comparisons, matched settings

Of the 55 pairs, these differ at a Holm-adjusted p below 0.05: on FAC under the sign-flip test over templates 36; on strict FAC (fully solved, sign flips over templates) 29, which agrees with FAC on 48; under McNemar's exact test on the instances 45; on MC with the judge under the sign-flip test 15 and under Wilcoxon's signed-rank test 17 (the two agree on 51); on MC by matching alone 20.

## Answers per model (Q1), matched settings

Verdict counts over the instances; the fully solved rate (strict FAC) with its template interval; the SD of the score over a template's 15 instances (median and quartiles over the templates) and between the template means; the templates with no instance variance, all solved and none solved.

| model | correct | partial | incorrect | unusable | empty | strict FAC | 95% CI | within-template SD, median (quartiles) | between-template SD | no instance variance | all solved | none solved |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 2217 | 13 | 12 | 8 | 8 | 0.985 | 0.976 to 0.993 | 0.000 (0.000 to 0.000) | 0.051 | 135 | 135 | 0 |
| `claude-sonnet-5` | 2214 | 14 | 20 | 2 | 2 | 0.984 | 0.974 to 0.992 | 0.000 (0.000 to 0.000) | 0.042 | 130 | 130 | 0 |
| `kimi-k3` | 2199 | 14 | 21 | 16 | 15 | 0.977 | 0.964 to 0.988 | 0.000 (0.000 to 0.000) | 0.075 | 125 | 125 | 0 |
| `glm-5.3-flash` | 2182 | 11 | 12 | 45 | 45 | 0.970 | 0.950 to 0.986 | 0.000 (0.000 to 0.000) | 0.110 | 128 | 128 | 0 |
| `gpt-5.4-mini` (reasoning on) | 2173 | 16 | 61 | 0 | 0 | 0.966 | 0.954 to 0.976 | 0.000 (0.000 to 0.129) | 0.068 | 108 | 108 | 0 |
| `muse-glimmer-30b` | 2176 | 9 | 34 | 31 | 31 | 0.967 | 0.948 to 0.983 | 0.000 (0.000 to 0.000) | 0.112 | 121 | 121 | 0 |
| `glm-5.3` | 2132 | 11 | 8 | 99 | 99 | 0.948 | 0.917 to 0.973 | 0.000 (0.000 to 0.000) | 0.171 | 124 | 122 | 2 |
| `gemma-4-26b-a4b` (reasoning on) | 2108 | 23 | 116 | 3 | 3 | 0.937 | 0.911 to 0.960 | 0.000 (0.000 to 0.129) | 0.149 | 110 | 110 | 0 |
| `gemini-3.1-flash-lite` (reasoning on) | 2051 | 23 | 176 | 0 | 0 | 0.912 | 0.879 to 0.941 | 0.000 (0.000 to 0.251) | 0.192 | 105 | 104 | 1 |
| `qwen3-235b-a22b-2507` | 1972 | 102 | 176 | 0 | 0 | 0.876 | 0.840 to 0.911 | 0.000 (0.000 to 0.258) | 0.178 | 79 | 77 | 0 |
| `gpt-oss-20b` | 1787 | 62 | 342 | 59 | 56 | 0.794 | 0.749 to 0.836 | 0.258 (0.000 to 0.352) | 0.268 | 52 | 47 | 5 |

Answer points lost by the first five on FAC: `deepseek-v4.1-flash` 26.5; `claude-sonnet-5` 29; `kimi-k3` 44; `glm-5.3-flash` 62.5; `gpt-5.4-mini` 69; 231 in all (no readable answer 71, partial 68 at half a point, incorrect 126). Their FAC spread: 0.019.

## Level gap (Q2), matched settings

Easy minus Advanced on the template means, templates resampled within each tier; Welch's t-test and the planned tier-label permutation, each with Holm over the 11 models, and the gap Welch detects at the strictest Holm step; then the gap with its Welch p (Holm) with unusable rows left out, without the symbolic templates and without the two chemical templates. The gap holds under Welch for 1 of 11 (`gpt-oss-20b`) and under the permutation for 4 (`glm-5.3-flash`, `glm-5.3`, `gemini-3.1-flash-lite`, `gpt-oss-20b`).

| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | permutation p (Holm) | detectable (Holm) | unusable excluded | without symbolic | without two chemical |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.995 | 0.979 | +0.016 | -0.001 to 0.041 | 0.7911 | 0.3412 | 0.040 | +0.002 (1.0000) | +0.018 (0.7540) | +0.005 (1.0000) |
| `claude-sonnet-5` | 0.993 | 0.984 | +0.008 | -0.007 to 0.024 | 1.0000 | 1.0000 | 0.029 | +0.004 (1.0000) | +0.005 (1.0000) | +0.005 (1.0000) |
| `kimi-k3` | 0.992 | 0.953 | +0.039 | -0.002 to 0.094 | 0.7911 | 0.2380 | 0.093 | +0.016 (1.0000) | +0.043 (0.7540) | +0.011 (1.0000) |
| `glm-5.3-flash` | 0.994 | 0.939 | +0.055 | 0.012 to 0.113 | 0.3679 | 0.0210 | 0.098 | +0.007 (1.0000) | +0.055 (0.5092) | +0.030 (0.7732) |
| `gpt-5.4-mini` (reasoning on) | 0.982 | 0.975 | +0.006 | -0.013 to 0.028 | 1.0000 | 1.0000 | 0.039 | +0.006 (1.0000) | +0.002 (1.0000) | +0.001 (1.0000) |
| `muse-glimmer-30b` | 0.976 | 0.958 | +0.019 | -0.026 to 0.067 | 1.0000 | 1.0000 | 0.088 | -0.007 (1.0000) | +0.020 (1.0000) | +0.017 (1.0000) |
| `glm-5.3` | 0.984 | 0.865 | +0.120 | 0.030 to 0.223 | 0.2164 | 0.0201 | 0.187 | -0.004 (1.0000) | +0.126 (0.2644) | +0.074 (0.7521) |
| `gemma-4-26b-a4b` (reasoning on) | 0.967 | 0.912 | +0.055 | -0.003 to 0.116 | 0.5284 | 0.3412 | 0.113 | +0.054 (0.7405) | +0.059 (0.5310) | +0.036 (1.0000) |
| `gemini-3.1-flash-lite` (reasoning on) | 0.959 | 0.826 | +0.132 | 0.049 to 0.225 | 0.0617 | 0.0091 | 0.169 | +0.132 (0.0617) | +0.133 (0.1022) | +0.108 (0.1943) |
| `qwen3-235b-a22b-2507` | 0.917 | 0.885 | +0.031 | -0.035 to 0.097 | 1.0000 | 1.0000 | 0.125 | +0.031 (1.0000) | +0.025 (1.0000) | +0.024 (1.0000) |
| `gpt-oss-20b` | 0.900 | 0.719 | +0.181 | 0.070 to 0.299 | 0.0353 | 0.0094 | 0.216 | +0.170 (0.0471) | +0.193 (0.0375) | +0.147 (0.1245) |

## Branches and domains, matched settings

Branch means of the template means with template intervals; the spread between the highest and the lowest branch; the smallest branch difference the design detects; the within-model branch pairs at a Holm-adjusted Welch p below 0.05; the lowest and the highest domain (item means, as RESULTS.md reports them).

| model | chemical engineering | civil engineering | electrical engineering | industrial engineering | mechanical engineering | spread | lowest branch | detectable | pairs at Holm 0.05 | lowest domain | highest domain |
|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|---|
| `deepseek-v4.1-flash` | 0.980 (0.953 to 0.996) | 0.993 (0.982 to 1.000) | 0.990 (0.980 to 0.998) | 0.980 (0.944 to 1.000) | 0.998 (0.993 to 1.000) | 0.018 | industrial_engineering | 0.037 | none | production_and_inventory 0.940 | fluid_mechanics 1.000 |
| `claude-sonnet-5` | 0.993 (0.987 to 1.000) | 0.996 (0.987 to 1.000) | 0.974 (0.952 to 0.991) | 0.973 (0.949 to 0.993) | 0.999 (0.997 to 1.000) | 0.026 | industrial_engineering | 0.030 | none | quality_and_reliability_control 0.940 | fluid_mechanics 1.000 |
| `kimi-k3` | 0.960 (0.900 to 0.996) | 0.982 (0.951 to 1.000) | 0.984 (0.973 to 0.993) | 0.980 (0.962 to 0.993) | 0.996 (0.990 to 1.000) | 0.036 | chemical_engineering | 0.055 | none | thermodynamics 0.917 | fluid_mechanics 1.000 |
| `glm-5.3-flash` | 0.960 (0.902 to 0.998) | 0.982 (0.949 to 1.000) | 0.981 (0.963 to 0.996) | 0.942 (0.869 to 0.989) | 0.996 (0.990 to 1.000) | 0.053 | industrial_engineering | 0.080 | none | production_and_inventory 0.887 | reaction_kinetics 1.000 |
| `gpt-5.4-mini` (reasoning on) | 0.961 (0.929 to 0.987) | 0.971 (0.949 to 0.989) | 0.981 (0.963 to 0.996) | 0.958 (0.922 to 0.984) | 0.976 (0.954 to 0.992) | 0.023 | industrial_engineering | 0.049 | none | production_and_inventory 0.933 | stochastic_operations 1.000 |
| `muse-glimmer-30b` | 0.963 (0.909 to 0.996) | 0.962 (0.924 to 0.991) | 0.989 (0.980 to 0.997) | 0.938 (0.864 to 0.989) | 0.993 (0.986 to 0.999) | 0.056 | industrial_engineering | 0.081 | none | production_and_inventory 0.880 | digital_communications 1.000 |
| `glm-5.3` | 0.900 (0.800 to 0.987) | 0.973 (0.927 to 1.000) | 0.961 (0.931 to 0.987) | 0.918 (0.831 to 0.984) | 0.998 (0.993 to 1.000) | 0.098 | chemical_engineering | 0.123 | none | thermodynamics 0.767 | fluid_mechanics 1.000 |
| `gemma-4-26b-a4b` (reasoning on) | 0.896 (0.804 to 0.967) | 0.960 (0.924 to 0.987) | 0.960 (0.930 to 0.984) | 0.920 (0.840 to 0.982) | 0.974 (0.956 to 0.991) | 0.079 | chemical_engineering | 0.107 | none | thermodynamics 0.806 | signals_and_systems 1.000 |
| `gemini-3.1-flash-lite` (reasoning on) | 0.846 (0.749 to 0.929) | 0.918 (0.833 to 0.978) | 0.947 (0.911 to 0.977) | 0.893 (0.800 to 0.967) | 0.980 (0.960 to 0.996) | 0.134 | chemical_engineering | 0.136 | none | thermodynamics 0.778 | structural_analysis 0.987 |
| `qwen3-235b-a22b-2507` | 0.869 (0.782 to 0.940) | 0.940 (0.911 to 0.964) | 0.911 (0.860 to 0.954) | 0.871 (0.773 to 0.951) | 0.904 (0.851 to 0.951) | 0.071 | chemical_engineering | 0.129 | none | transport_phenomena 0.783 | reaction_kinetics 0.980 |
| `gpt-oss-20b` | 0.798 (0.664 to 0.907) | 0.671 (0.540 to 0.796) | 0.897 (0.852 to 0.937) | 0.816 (0.738 to 0.882) | 0.859 (0.791 to 0.921) | 0.226 | civil_engineering | 0.188 | civil_engineering - electrical_engineering -0.226 | geotechnical_engineering 0.513 | reaction_kinetics 0.953 |

Domain means (item means) per model:

| model | digital communications | electromagnetics and waves | fluid mechanics | geotechnical engineering | mechanics of materials | production and inventory | quality and reliability control | reaction kinetics | signals and systems | stochastic operations | structural analysis | thermodynamics | transport phenomena | vibrations and acoustics | water resources |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.990 | 0.990 | 1.000 | 1.000 | 1.000 | 0.940 | 1.000 | 0.993 | 0.990 | 1.000 | 1.000 | 0.961 | 0.992 | 0.993 | 0.980 |
| `claude-sonnet-5` | 0.980 | 0.953 | 1.000 | 0.987 | 1.000 | 0.987 | 0.940 | 1.000 | 0.990 | 0.993 | 1.000 | 0.989 | 0.992 | 0.997 | 1.000 |
| `kimi-k3` | 0.987 | 0.977 | 1.000 | 1.000 | 0.993 | 0.967 | 0.980 | 0.993 | 0.990 | 0.993 | 0.993 | 0.917 | 0.983 | 0.993 | 0.953 |
| `glm-5.3-flash` | 0.980 | 0.977 | 0.997 | 0.993 | 0.993 | 0.887 | 0.947 | 1.000 | 0.987 | 0.993 | 1.000 | 0.906 | 0.992 | 0.997 | 0.953 |
| `gpt-5.4-mini` (reasoning on) | 0.967 | 0.980 | 0.990 | 0.960 | 0.953 | 0.933 | 0.940 | 0.993 | 0.997 | 1.000 | 0.987 | 0.944 | 0.946 | 0.983 | 0.967 |
| `muse-glimmer-30b` | 1.000 | 0.987 | 0.997 | 0.927 | 0.997 | 0.880 | 0.933 | 1.000 | 0.980 | 1.000 | 0.987 | 0.922 | 0.979 | 0.987 | 0.973 |
| `glm-5.3` | 0.947 | 0.943 | 1.000 | 0.993 | 1.000 | 0.880 | 0.907 | 0.993 | 0.993 | 0.967 | 1.000 | 0.767 | 0.983 | 0.993 | 0.927 |
| `gemma-4-26b-a4b` (reasoning on) | 0.970 | 0.910 | 0.973 | 0.953 | 0.963 | 0.907 | 0.860 | 0.960 | 1.000 | 0.993 | 0.973 | 0.806 | 0.950 | 0.987 | 0.953 |
| `gemini-3.1-flash-lite` (reasoning on) | 0.947 | 0.927 | 0.977 | 0.973 | 0.977 | 0.920 | 0.813 | 0.903 | 0.967 | 0.947 | 0.987 | 0.778 | 0.875 | 0.987 | 0.793 |
| `qwen3-235b-a22b-2507` | 0.897 | 0.890 | 0.913 | 0.947 | 0.943 | 0.847 | 0.787 | 0.980 | 0.947 | 0.980 | 0.947 | 0.833 | 0.783 | 0.857 | 0.927 |
| `gpt-oss-20b` | 0.933 | 0.867 | 0.887 | 0.513 | 0.833 | 0.767 | 0.773 | 0.953 | 0.890 | 0.907 | 0.733 | 0.650 | 0.825 | 0.857 | 0.767 |

## Wrong answers, coverage and flags (Q3), matched settings

RESULTS.md's Q3 on each model's matched store. Wrong answers score 0; coverage on the readable ones with milestones, by matching alone (E3) beside its chance floor (a sibling instance's milestones) and with the judge (E5-strict), template intervals; the judge's reached share and the share of milestones it decided; the digit rule on answered wrong answers; the router's judge flags on fully solved responses with their template interval, its steps flagged per answered response, the router's flags on fully solved and on wrong responses; the share of answered wrong answers each component points at (digit rule / E5 missing a milestone / router judge, responses in brackets); the first flagged step's position on fully solved responses (median, 0 the first step; responses in brackets); arithmetic claims per answered response.

| model | wrong | with milestones | readable | E3 on readable wrong | 95% CI | E3 floor | E5-strict on readable wrong | 95% CI | E5 responses | reached share | judge-decided share | digit rule on wrong | router judge on fully solved | 95% CI | router steps flagged per response | router on fully solved | router on wrong | attribution | first flag, median | claims per response |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| `deepseek-v4.1-flash` | 20 | 20 | 12 | 0.667 | 0.293 to 0.807 | 0.141 | 0.667 | 0.293 to 0.811 | 20 | 0.253 | 0.169 | 0.000 | 0.010 | 0.005 to 0.015 | 0.026 | 0.015 | 0.417 | 0.000 / 0.833 / 0.417 (12) | 0.667 (13) | 3.27 |
| `claude-sonnet-5` | 22 | 22 | 20 | 0.850 | 0.738 to 0.961 | 0.150 | 0.869 | 0.774 to 0.969 | 22 | 0.218 | 0.125 | 0.200 | 0.054 | 0.041 to 0.068 | 0.275 | 0.186 | 0.650 | 0.200 / 0.200 / 0.500 (20) | 0.667 (316) | 5.32 |
| `kimi-k3` | 37 | 36 | 20 | 0.584 | 0.411 to 0.752 | 0.113 | 0.584 | 0.409 to 0.748 | 36 | 0.276 | 0.150 | 0.000 | 0.017 | 0.011 to 0.024 | 0.049 | 0.024 | 0.545 | 0.000 / 0.636 / 0.545 (22) | 0.571 (15) | 2.79 |
| `glm-5.3-flash` | 57 | 57 | 12 | 0.378 | 0.156 to 0.606 | 0.134 | 0.458 | 0.241 to 0.665 | 57 | 0.215 | 0.126 | 0.250 | 0.032 | 0.023 to 0.042 | 0.072 | 0.049 | 0.333 | 0.250 / 0.750 / 0.250 (12) | 0.750 (41) | 4.74 |
| `gpt-5.4-mini` (reasoning on) | 61 | 61 | 61 | 0.640 | 0.526 to 0.750 | 0.137 | 0.689 | 0.572 to 0.800 | 61 | 0.255 | 0.184 | 0.246 | 0.115 | 0.091 to 0.140 | 0.268 | 0.178 | 0.770 | 0.246 / 0.459 / 0.672 (61) | 0.600 (163) | 3.27 |
| `muse-glimmer-30b` | 65 | 64 | 33 | 0.583 | 0.479 to 0.744 | 0.178 | 0.653 | 0.546 to 0.858 | 64 | 0.275 | 0.167 | 0.029 | 0.059 | 0.046 to 0.074 | 0.153 | 0.113 | 0.529 | 0.029 / 0.588 / 0.500 (34) | 0.571 (127) | 3.27 |
| `glm-5.3` | 107 | 107 | 8 | 0.486 | 0.188 to 0.781 | 0.146 | 0.596 | 0.292 to 0.875 | 107 | 0.235 | 0.107 | 0.125 | 0.026 | 0.016 to 0.036 | 0.055 | 0.042 | 0.250 | 0.125 / 0.500 / 0.125 (8) | 0.727 (39) | 5.72 |
| `gemma-4-26b-a4b` (reasoning on) | 119 | 119 | 116 | 0.548 | 0.468 to 0.635 | 0.099 | 0.569 | 0.492 to 0.658 | 119 | 0.184 | 0.186 | 0.448 | 0.074 | 0.053 to 0.097 | 0.350 | 0.159 | 0.853 | 0.448 / 0.724 / 0.707 (116) | 0.500 (208) | 3.30 |
| `gemini-3.1-flash-lite` (reasoning on) | 176 | 174 | 174 | 0.515 | 0.436 to 0.600 | 0.116 | 0.540 | 0.466 to 0.626 | 174 | 0.169 | 0.187 | 0.358 | 0.093 | 0.068 to 0.121 | 0.367 | 0.153 | 0.875 | 0.358 / 0.750 / 0.773 (176) | 0.400 (142) | 3.54 |
| `qwen3-235b-a22b-2507` | 176 | 173 | 173 | 0.586 | 0.504 to 0.676 | 0.175 | 0.616 | 0.538 to 0.697 | 173 | 0.163 | 0.158 | 0.449 | 0.107 | 0.084 to 0.131 | 0.864 | 0.297 | 0.864 | 0.449 / 0.631 / 0.722 (176) | 0.600 (457) | 10.21 |
| `gpt-oss-20b` | 401 | 392 | 335 | 0.455 | 0.389 to 0.531 | 0.117 | 0.465 | 0.399 to 0.537 | 392 | 0.106 | 0.238 | 0.229 | 0.206 | 0.177 to 0.237 | 1.006 | 0.311 | 0.916 | 0.229 / 0.777 / 0.881 (345) | 0.667 (263) | 3.77 |

Wrong-answer rate by the instance's milestone count (instances in brackets):

| model | 0 milestones | 1 milestones | 2 milestones | 3 milestones | 4-5 milestones | 6+ milestones |
|---|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.000 (57) | 0.000 (293) | 0.000 (462) | 0.000 (368) | 0.006 (542) | 0.032 (528) |
| `claude-sonnet-5` | 0.000 (57) | 0.000 (293) | 0.017 (462) | 0.003 (368) | 0.013 (542) | 0.011 (528) |
| `kimi-k3` | 0.018 (57) | 0.000 (293) | 0.006 (462) | 0.014 (368) | 0.028 (542) | 0.025 (528) |
| `glm-5.3-flash` | 0.000 (57) | 0.003 (293) | 0.017 (462) | 0.011 (368) | 0.022 (542) | 0.061 (528) |
| `gpt-5.4-mini` (reasoning on) | 0.000 (57) | 0.007 (293) | 0.022 (462) | 0.014 (368) | 0.048 (542) | 0.034 (528) |
| `muse-glimmer-30b` | 0.018 (57) | 0.003 (293) | 0.017 (462) | 0.038 (368) | 0.030 (542) | 0.047 (528) |
| `glm-5.3` | 0.000 (57) | 0.003 (293) | 0.013 (462) | 0.027 (368) | 0.057 (542) | 0.112 (528) |
| `gemma-4-26b-a4b` (reasoning on) | 0.000 (57) | 0.017 (293) | 0.041 (462) | 0.046 (368) | 0.042 (542) | 0.104 (528) |
| `gemini-3.1-flash-lite` (reasoning on) | 0.035 (57) | 0.007 (293) | 0.039 (462) | 0.057 (368) | 0.098 (542) | 0.152 (528) |
| `qwen3-235b-a22b-2507` | 0.053 (57) | 0.034 (293) | 0.056 (462) | 0.073 (368) | 0.085 (542) | 0.121 (528) |
| `gpt-oss-20b` | 0.158 (57) | 0.061 (293) | 0.091 (462) | 0.190 (368) | 0.236 (542) | 0.254 (528) |

## Coverage over every response with milestones (q3_overall), matched settings

E3 and E5-strict over every response whose instance has milestones (an unusable response scores what it reached, an empty one nothing), and over the readable ones; template intervals.

| model | responses | E3, all | 95% CI | E3, readable | 95% CI | E5-strict, all | 95% CI | E5-strict, readable | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 2193 | 0.860 | 0.830 to 0.889 | 0.864 | 0.835 to 0.892 | 0.898 | 0.873 to 0.922 | 0.901 | 0.876 to 0.924 |
| `claude-sonnet-5` | 2193 | 0.903 | 0.880 to 0.925 | 0.904 | 0.880 to 0.925 | 0.923 | 0.902 to 0.942 | 0.923 | 0.903 to 0.942 |
| `kimi-k3` | 2193 | 0.877 | 0.849 to 0.902 | 0.883 | 0.858 to 0.906 | 0.911 | 0.887 to 0.933 | 0.918 | 0.896 to 0.938 |
| `glm-5.3-flash` | 2193 | 0.891 | 0.863 to 0.917 | 0.909 | 0.886 to 0.930 | 0.911 | 0.885 to 0.935 | 0.930 | 0.909 to 0.950 |
| `gpt-5.4-mini` (reasoning on) | 2193 | 0.848 | 0.817 to 0.878 | 0.848 | 0.817 to 0.878 | 0.893 | 0.868 to 0.918 | 0.893 | 0.867 to 0.918 |
| `muse-glimmer-30b` | 2193 | 0.862 | 0.832 to 0.890 | 0.874 | 0.848 to 0.900 | 0.897 | 0.870 to 0.922 | 0.910 | 0.887 to 0.932 |
| `glm-5.3` | 2193 | 0.880 | 0.845 to 0.912 | 0.922 | 0.899 to 0.942 | 0.899 | 0.865 to 0.930 | 0.941 | 0.921 to 0.960 |
| `gemma-4-26b-a4b` (reasoning on) | 2193 | 0.850 | 0.818 to 0.880 | 0.851 | 0.820 to 0.881 | 0.884 | 0.855 to 0.909 | 0.885 | 0.858 to 0.910 |
| `gemini-3.1-flash-lite` (reasoning on) | 2193 | 0.856 | 0.824 to 0.886 | 0.856 | 0.825 to 0.886 | 0.881 | 0.852 to 0.909 | 0.881 | 0.853 to 0.908 |
| `qwen3-235b-a22b-2507` | 2193 | 0.866 | 0.835 to 0.894 | 0.866 | 0.836 to 0.894 | 0.892 | 0.866 to 0.915 | 0.892 | 0.867 to 0.916 |
| `gpt-oss-20b` | 2193 | 0.790 | 0.751 to 0.827 | 0.811 | 0.775 to 0.844 | 0.812 | 0.775 to 0.848 | 0.833 | 0.800 to 0.865 |

## Coverage across models (q3_coverage), matched settings

As Q3's coverage table (E5-strict, 148 templates): the answer and coverage ranks, Spearman's rho across templates between coverage and the steps and claims per readable response, coverage on fully solved and on wrong readable responses (responses in brackets), the share of fully solved responses below 0.5 coverage and of wrong ones at full coverage, template intervals. Kendall's tau between the answer and the coverage orderings: 0.709 (95% CI 0.345 to 0.818).

| model | MC | 95% CI | answer rank | coverage rank | rho, steps | rho, claims | median steps | on fully solved | on wrong | solved below 0.5 | 95% CI | wrong at 1.0 | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.899 | 0.874 to 0.922 | 1 | 5 | -0.239 | -0.143 | 6 | 0.903 (2164) | 0.667 (12) | 0.034 | 0.017 to 0.054 | 0.083 | 0.000 to 0.429 |
| `claude-sonnet-5` | 0.922 | 0.902 to 0.942 | 2 | 1 | -0.183 | -0.233 | 7 | 0.926 (2157) | 0.869 (20) | 0.025 | 0.012 to 0.040 | 0.550 | 0.250 to 0.875 |
| `kimi-k3` | 0.912 | 0.888 to 0.934 | 3 | 3 | -0.213 | -0.129 | 6 | 0.921 (2145) | 0.584 (20) | 0.028 | 0.012 to 0.047 | 0.200 | 0.000 to 0.444 |
| `glm-5.3-flash` | 0.912 | 0.886 to 0.936 | 4 | 2 | -0.095 | -0.211 | 7 | 0.933 (2126) | 0.458 (12) | 0.023 | 0.008 to 0.042 | 0.083 | 0.000 to 0.273 |
| `gpt-5.4-mini` (reasoning on) | 0.892 | 0.866 to 0.917 | 5 | 7 | -0.279 | -0.160 | 6 | 0.899 (2118) | 0.689 (61) | 0.040 | 0.021 to 0.063 | 0.311 | 0.159 to 0.485 |
| `muse-glimmer-30b` | 0.897 | 0.871 to 0.922 | 6 | 6 | -0.208 | -0.148 | 5 | 0.914 (2123) | 0.653 (33) | 0.030 | 0.016 to 0.047 | 0.364 | 0.167 to 0.750 |
| `glm-5.3` | 0.900 | 0.866 to 0.930 | 7 | 4 | -0.147 | -0.219 | 8 | 0.943 (2075) | 0.596 (8) | 0.020 | 0.004 to 0.040 | 0.500 | 0.125 to 0.875 |
| `gemma-4-26b-a4b` (reasoning on) | 0.884 | 0.856 to 0.910 | 8 | 9 | -0.233 | -0.212 | 5 | 0.903 (2055) | 0.569 (116) | 0.041 | 0.020 to 0.067 | 0.181 | 0.073 to 0.325 |
| `gemini-3.1-flash-lite` (reasoning on) | 0.881 | 0.851 to 0.908 | 9 | 10 | -0.184 | -0.238 | 5 | 0.911 (2000) | 0.540 (174) | 0.036 | 0.017 to 0.059 | 0.161 | 0.074 to 0.272 |
| `qwen3-235b-a22b-2507` | 0.892 | 0.867 to 0.916 | 10 | 8 | -0.189 | -0.183 | 10 | 0.914 (1919) | 0.616 (173) | 0.036 | 0.017 to 0.060 | 0.249 | 0.155 to 0.368 |
| `gpt-oss-20b` | 0.812 | 0.775 to 0.847 | 11 | 11 | -0.262 | -0.137 | 9 | 0.907 (1742) | 0.465 (335) | 0.033 | 0.017 to 0.052 | 0.146 | 0.088 to 0.221 |

## Consistency within a template (Q4), matched settings

The share of the 58 single-path templates and of the 92 others fully solved on all 15 instances, on some and on none; 95% intervals, templates resampled.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.931 (0.862 to 0.983) | 0.069 (0.017 to 0.138) | 0.000 (0.000 to 0.000) | 0.880 (0.815 to 0.946) | 0.120 (0.054 to 0.185) | 0.000 (0.000 to 0.000) |
| `claude-sonnet-5` | 0.914 (0.844 to 0.983) | 0.086 (0.017 to 0.155) | 0.000 (0.000 to 0.000) | 0.837 (0.761 to 0.913) | 0.163 (0.087 to 0.239) | 0.000 (0.000 to 0.000) |
| `kimi-k3` | 0.862 (0.776 to 0.948) | 0.138 (0.052 to 0.225) | 0.000 (0.000 to 0.000) | 0.815 (0.728 to 0.891) | 0.185 (0.109 to 0.261) | 0.000 (0.000 to 0.000) |
| `glm-5.3-flash` | 0.897 (0.810 to 0.966) | 0.103 (0.034 to 0.190) | 0.000 (0.000 to 0.000) | 0.826 (0.739 to 0.902) | 0.174 (0.098 to 0.250) | 0.000 (0.000 to 0.000) |
| `gpt-5.4-mini` (reasoning on) | 0.655 (0.534 to 0.776) | 0.345 (0.224 to 0.466) | 0.000 (0.000 to 0.000) | 0.761 (0.674 to 0.848) | 0.239 (0.152 to 0.326) | 0.000 (0.000 to 0.000) |
| `muse-glimmer-30b` | 0.776 (0.672 to 0.879) | 0.224 (0.121 to 0.328) | 0.000 (0.000 to 0.000) | 0.826 (0.750 to 0.902) | 0.174 (0.098 to 0.261) | 0.000 (0.000 to 0.000) |
| `glm-5.3` | 0.845 (0.741 to 0.931) | 0.138 (0.052 to 0.224) | 0.017 (0.000 to 0.052) | 0.793 (0.707 to 0.880) | 0.196 (0.120 to 0.283) | 0.011 (0.000 to 0.033) |
| `gemma-4-26b-a4b` (reasoning on) | 0.759 (0.655 to 0.862) | 0.241 (0.138 to 0.345) | 0.000 (0.000 to 0.000) | 0.717 (0.620 to 0.804) | 0.283 (0.196 to 0.380) | 0.000 (0.000 to 0.000) |
| `gemini-3.1-flash-lite` (reasoning on) | 0.776 (0.655 to 0.879) | 0.207 (0.103 to 0.310) | 0.017 (0.000 to 0.052) | 0.641 (0.543 to 0.739) | 0.359 (0.261 to 0.457) | 0.000 (0.000 to 0.000) |
| `qwen3-235b-a22b-2507` | 0.448 (0.328 to 0.586) | 0.534 (0.414 to 0.655) | 0.017 (0.000 to 0.052) | 0.554 (0.457 to 0.652) | 0.424 (0.326 to 0.522) | 0.022 (0.000 to 0.054) |
| `gpt-oss-20b` | 0.345 (0.224 to 0.466) | 0.586 (0.466 to 0.707) | 0.069 (0.017 to 0.138) | 0.293 (0.207 to 0.391) | 0.696 (0.598 to 0.783) | 0.011 (0.000 to 0.033) |

## Paraphrase (Q5), matched settings, one family

Paraphrase minus original, paired by instance on the 275 pairs the experts kept, each model at its matched store (a re-run model's reasoning-on rows against its reasoning-on paraphrase store, the others' default rows against `scores/paraphrase`): the item mean, the template bootstrap at 95% and 90%, the sign-flip test over templates and McNemar's test on the fully solved verdict, each with Holm over the 11 models; the 90% interval read against the margin of ±0.05 (D-165). Within the margin: 9 of 11. Kendall's tau between the paraphrase and the original orderings: 0.722 (95% CI 0.449 to 0.844, 275 instances); its noise floor at the arm's size and template mix, median 0.807 (quartiles 0.748 to 0.860).

| model | paraphrase store | items | templates | change | 95% CI | 90% CI | within ±0.05 | p (Holm) | McNemar p (Holm) | E3 change | E5-strict change | same endpoint: change (items) |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | paraphrase | 275 | 114 | -0.013 | -0.029 to 0.000 | -0.026 to -0.002 | yes | 1.0000 | 1.0000 | -0.000 | -0.009 | -0.030 (33) |
| `claude-sonnet-5` | paraphrase | 275 | 114 | -0.011 | -0.025 to 0.000 | -0.022 to -0.003 | yes | 1.0000 | 1.0000 | -0.007 | -0.005 | -0.011 (275) |
| `kimi-k3` | paraphrase | 275 | 114 | +0.002 | -0.005 to 0.011 | -0.004 to 0.009 | yes | 1.0000 | 1.0000 | -0.007 | -0.005 | +0.167 (6) |
| `glm-5.3-flash` | paraphrase | 275 | 114 | -0.011 | -0.029 to 0.005 | -0.026 to 0.004 | yes | 1.0000 | 1.0000 | -0.028 | -0.015 | +0.000 (25) |
| `gpt-5.4-mini` (reasoning on) | paraphrase-reasoning-medium | 275 | 114 | +0.005 | -0.012 to 0.024 | -0.009 to 0.021 | yes | 1.0000 | 1.0000 | +0.016 | +0.023 | +0.005 (275) |
| `muse-glimmer-30b` | paraphrase | 275 | 114 | +0.025 | 0.004 to 0.051 | 0.007 to 0.046 | yes | 0.5273 | 0.7197 | -0.009 | -0.007 | +0.025 (275) |
| `glm-5.3` | paraphrase | 275 | 114 | +0.002 | -0.018 to 0.022 | -0.015 to 0.019 | yes | 1.0000 | 1.0000 | +0.007 | +0.001 | +0.011 (88) |
| `gemma-4-26b-a4b` (reasoning on) | paraphrase-reasoning-medium | 275 | 114 | +0.000 | -0.033 to 0.028 | -0.026 to 0.023 | yes | 1.0000 | 1.0000 | +0.003 | -0.009 | +0.009 (235) |
| `gemini-3.1-flash-lite` (reasoning on) | paraphrase-reasoning-medium | 275 | 114 | -0.013 | -0.031 to 0.004 | -0.028 to 0.000 | yes | 1.0000 | 1.0000 | -0.011 | -0.006 | -0.008 (182) |
| `qwen3-235b-a22b-2507` | paraphrase | 275 | 114 | -0.036 | -0.080 to 0.005 | -0.073 to -0.002 | no | 1.0000 | 1.0000 | -0.005 | -0.015 | -0.033 (184) |
| `gpt-oss-20b` | paraphrase | 275 | 114 | +0.018 | -0.031 to 0.067 | -0.023 to 0.059 | no | 1.0000 | 1.0000 | -0.013 | -0.021 | +0.018 (273) |

Beside the run-to-run differences (each repeat minus the matched rows on the repeat instances): `gpt-oss-20b` paraphrase +0.018, repeats +0.008, +0.012, +0.007 (beyond the largest absolute repeat difference); `gemma-4-26b-a4b` paraphrase +0.000, repeats +0.017, +0.010, +0.002 (within the largest absolute repeat difference); `qwen3-235b-a22b-2507` paraphrase -0.036, repeats +0.003, +0.000, +0.010 (beyond the largest absolute repeat difference); `gemini-3.1-flash-lite` paraphrase -0.013, repeats +0.012, +0.005, +0.018 (within the largest absolute repeat difference).

## Decoding repeats, matched settings

Each model's repeat stores at its matched setting (a re-run model's reasoning-on repeats, the others' default repeats), on the instances every repeat holds: the score per repeat, the matched rows on the same instances, the SD and the range over the repeats, and the share of instances with the same verdict in every repeat.

| model | items | repeats | matched rows on the items | SD | range | same verdict every repeat |
|---|---:|---|---:|---:|---:|---:|
| `gemma-4-26b-a4b` (reasoning on) | 300 | repeat1-reasoning-medium 0.955, repeat2-reasoning-medium 0.948, repeat3-reasoning-medium 0.940 | 0.938 | 0.008 | 0.015 | 0.917 |
| `gemini-3.1-flash-lite` (reasoning on) | 300 | repeat1-reasoning-medium 0.927, repeat2-reasoning-medium 0.920, repeat3-reasoning-medium 0.933 | 0.915 | 0.007 | 0.013 | 0.920 |
| `qwen3-235b-a22b-2507` | 300 | repeat1 0.910, repeat2 0.907, repeat3 0.917 | 0.907 | 0.005 | 0.010 | 0.893 |
| `gpt-oss-20b` | 300 | repeat1 0.835, repeat2 0.838, repeat3 0.833 | 0.827 | 0.003 | 0.005 | 0.910 |

## Sensitivity, matched settings

Each model's item mean under each reading of the answer check and each pool, and Kendall's tau of each ordering with the headline (as scored); the tau noise floor (two halves of each template's instances): median 0.917, quartiles 0.891 to 0.927.

| model | as scored | half tol | double tol | fully solved | unusable excluded | without shortcut templates | without symbolic templates | half unit | whole trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.988 | 0.986 | 0.990 | 0.985 | 0.992 | 0.989 | 0.987 | 0.983 | 0.989 |
| `claude-sonnet-5` | 0.987 | 0.983 | 0.989 | 0.984 | 0.988 | 0.987 | 0.988 | 0.980 | 0.988 |
| `kimi-k3` | 0.980 | 0.979 | 0.982 | 0.977 | 0.987 | 0.981 | 0.979 | 0.975 | 0.982 |
| `glm-5.3-flash` | 0.972 | 0.971 | 0.973 | 0.970 | 0.992 | 0.971 | 0.972 | 0.968 | 0.974 |
| `gpt-5.4-mini` (reasoning on) | 0.969 | 0.965 | 0.978 | 0.966 | 0.969 | 0.968 | 0.969 | 0.950 | 0.970 |
| `muse-glimmer-30b` | 0.969 | 0.967 | 0.970 | 0.967 | 0.983 | 0.969 | 0.968 | 0.952 | 0.970 |
| `glm-5.3` | 0.950 | 0.949 | 0.951 | 0.948 | 0.994 | 0.949 | 0.949 | 0.944 | 0.951 |
| `gemma-4-26b-a4b` (reasoning on) | 0.942 | 0.927 | 0.952 | 0.937 | 0.943 | 0.941 | 0.939 | 0.933 | 0.944 |
| `gemini-3.1-flash-lite` (reasoning on) | 0.917 | 0.906 | 0.924 | 0.912 | 0.917 | 0.916 | 0.916 | 0.909 | 0.919 |
| `qwen3-235b-a22b-2507` | 0.899 | 0.888 | 0.910 | 0.876 | 0.899 | 0.897 | 0.898 | 0.889 | 0.917 |
| `gpt-oss-20b` | 0.808 | 0.801 | 0.815 | 0.794 | 0.830 | 0.806 | 0.802 | 0.780 | 0.811 |
| tau with the headline | | 0.964 | 0.927 | 0.964 | 0.636 | 0.964 | 0.964 | 0.964 | 1.000 |
