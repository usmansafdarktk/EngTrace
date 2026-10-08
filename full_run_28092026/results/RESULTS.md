# Results of the full run

Printed by `analyze.py` from the score store (`score.py`) under ANALYSIS_PLAN.md (D-117); the method choices the plan leaves open are fixed in the script's docstring, and the corrections and additions made after the first results were read are labelled where they appear (D-146 to D-149, D-170). The answer check is the one corrected after the run was read, D-137 to D-139, D-147 and D-169 (the run's README and `ANSWER_FORM_FIX.md`). The answer check, E3 and E4's digit rule, E5 and the step router. Every interval is 95% and resamples templates (B = 10,000): all 150 for a model's score, and the templates with a qualifying trace for a rate over a subset of traces. Every resampling test draws 100,000 permutations, so its Holm-adjusted floor over 55 pairs is 0.0006. **Incomplete judged stages:** `qwen3-235b-a22b-2507` have calls without a reply; their judged rates are over the answered calls only.

## Q1. Answer score per model

Correct 1, partial 0.5, incorrect or unusable 0. Fully solved: a partial scores 0. SD within: the mean over templates of the SD of the score across a template's 15 items; SD between: the SD of the 150 template means (the per-template variability the July rebuttal promised). Unusable is empty plus unreadable; capped: answered rows that stopped at the output cap and are scored on what they state (D-117).

| model | answer score | 95% CI | SD within | SD between | fully solved | 95% CI | correct | partial | incorrect | unusable (empty + unreadable) | capped, scored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `deepseek-v4.1-flash` | 0.988 | 0.979 to 0.995 | 0.026 | 0.051 | 0.985 | 0.976 to 0.993 | 2217 | 13 | 12 | 8 (8 + 0) | 1 |
| `claude-sonnet-5` | 0.987 | 0.980 to 0.993 | 0.035 | 0.042 | 0.984 | 0.974 to 0.992 | 2214 | 14 | 20 | 2 (2 + 0) | 0 |
| `kimi-k3` | 0.980 | 0.967 to 0.990 | 0.043 | 0.075 | 0.977 | 0.964 to 0.988 | 2199 | 14 | 21 | 16 (15 + 1) | 2 |
| `glm-5.3-flash` | 0.972 | 0.953 to 0.987 | 0.043 | 0.110 | 0.970 | 0.950 to 0.986 | 2182 | 11 | 12 | 45 (45 + 0) | 2 |
| `muse-glimmer-30b` | 0.969 | 0.949 to 0.985 | 0.054 | 0.112 | 0.967 | 0.948 to 0.983 | 2176 | 9 | 34 | 31 (31 + 0) | 0 |
| `glm-5.3` | 0.950 | 0.921 to 0.975 | 0.053 | 0.171 | 0.948 | 0.917 to 0.973 | 2132 | 11 | 8 | 99 (99 + 0) | 5 |
| `qwen3-235b-a22b-2507` | 0.899 | 0.869 to 0.926 | 0.148 | 0.178 | 0.876 | 0.840 to 0.911 | 1972 | 102 | 176 | 0 (0 + 0) | 0 |
| `gemini-3.1-flash-lite` | 0.878 | 0.837 to 0.915 | 0.124 | 0.244 | 0.872 | 0.830 to 0.909 | 1961 | 30 | 259 | 0 (0 + 0) | 0 |
| `gemma-4-26b-a4b` | 0.872 | 0.835 to 0.905 | 0.164 | 0.216 | 0.860 | 0.824 to 0.894 | 1934 | 56 | 260 | 0 (0 + 0) | 1 |
| `gpt-5.4-mini` | 0.852 | 0.810 to 0.890 | 0.159 | 0.254 | 0.842 | 0.800 to 0.882 | 1894 | 45 | 311 | 0 (0 + 0) | 0 |
| `gpt-oss-20b` | 0.808 | 0.764 to 0.849 | 0.221 | 0.268 | 0.794 | 0.749 to 0.836 | 1787 | 62 | 342 | 59 (56 + 3) | 5 |

**Instance variance within a template** (*added 2026-10-03, D-178*). The SD-within column above is the mean over templates of the SD of the score across a template's 15 instances; here its distribution: the quartiles over the 150 templates, the templates with no instance variance at all (the same verdict on all 15 instances, split into all solved and none solved), and the three templates with the largest instance variance. `results/per_template.csv` carries every template's SD (`answer_score_sd`).

| model | per-template SD: lower quartile / median / upper quartile | templates with no instance variance (all solved / none solved) | highest instance variance |
|---|---|---|---|
| `deepseek-v4.1-flash` | 0.000 / 0.000 / 0.000 | 135 (135 / 0) | `qr_policy_one_iteration` 0.52; `work_isothermal_virial` 0.49; `newsvendor_normal_demand` 0.35 |
| `claude-sonnet-5` | 0.000 / 0.000 / 0.000 | 130 (130 / 0) | `single_sampling_oc_point` 0.46; `chart_pair_selection` 0.41; `ber_estimation_mary` 0.35 |
| `kimi-k3` | 0.000 / 0.000 / 0.000 | 125 (125 / 0) | `normal_depth_iteration` 0.51; `adiabatic_flame_temperature` 0.46; `work_isothermal_virial` 0.46 |
| `glm-5.3-flash` | 0.000 / 0.000 / 0.000 | 128 (128 / 0) | `normal_depth_iteration` 0.52; `work_isothermal_virial` 0.46; `single_sampling_oc_point` 0.46 |
| `muse-glimmer-30b` | 0.000 / 0.000 / 0.000 | 121 (121 / 0) | `aoq_ati_rectifying` 0.52; `primary_consolidation_settlement` 0.51; `relative_density_of_sand` 0.49 |
| `glm-5.3` | 0.000 / 0.000 / 0.000 | 124 (122 / 2) | `aoq_ati_rectifying` 0.49; `mm1k_finite_capacity` 0.49; `normal_depth_iteration` 0.49 |
| `qwen3-235b-a22b-2507` | 0.000 / 0.000 / 0.258 | 79 (77 / 0) | `annulus_flowrate` 0.52; `aoq_ati_rectifying` 0.51; `chart_pair_selection` 0.51 |
| `gemini-3.1-flash-lite` | 0.000 / 0.000 / 0.258 | 97 (94 / 3) | `ber_estimation_mary` 0.52; `manning_rectangular_discharge` 0.52; `power_law_fluid_shear` 0.52 |
| `gemma-4-26b-a4b` | 0.000 / 0.000 / 0.306 | 77 (75 / 2) | `aoq_ati_rectifying` 0.52; `ideal_gas_volume` 0.52; `utube_manometer` 0.52 |
| `gpt-5.4-mini` | 0.000 / 0.000 / 0.352 | 79 (73 / 6) | `annulus_flowrate` 0.52; `chart_pair_selection` 0.52; `rotating_unbalance` 0.52 |
| `gpt-oss-20b` | 0.000 / 0.258 / 0.352 | 52 (47 / 5) | `aoq_ati_rectifying` 0.52; `phase_relations_degree_of_saturation` 0.52; `statically_indeterminate_shaft` 0.52 |

Of the 55 pairs, 38 differ at a Holm-adjusted p below 0.05 on the answer score (sign-flip permutation over the 150 per-template differences). The same template-level test on the fully-solved rate gives the same verdict on 51 of 55. McNemar's exact test on the paired item verdicts (Holm) gives the same verdict on 48 of 55; on 7 of the 7 others McNemar holds and the template-level test does not: McNemar treats the 2,250 items as independent, and instances of a template are not (D-111), so a claim rests on the template-level tests. 6 answered rows across the roster ended with a finish reason other than stop or length (a provider fault inside a 200) and are scored on what they state; the harness now retries such a reply (D-148).

| a | b | a - b | 95% CI | p (Holm) | fully solved a - b | 95% CI | template p (Holm) | McNemar p (Holm) | same verdict | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.180 | -0.223 to -0.140 | 0.0005 | -0.191 | -0.235 to -0.150 | 0.0005 | 0.0000 | yes | 0.060 / 0.089 |
| `gpt-oss-20b` | `claude-sonnet-5` | -0.179 | -0.222 to -0.140 | 0.0005 | -0.190 | -0.232 to -0.150 | 0.0005 | 0.0000 | yes | 0.059 / 0.088 |
| `gpt-oss-20b` | `kimi-k3` | -0.172 | -0.215 to -0.133 | 0.0005 | -0.183 | -0.227 to -0.142 | 0.0005 | 0.0000 | yes | 0.059 / 0.088 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.164 | -0.207 to -0.124 | 0.0005 | -0.176 | -0.219 to -0.135 | 0.0005 | 0.0000 | yes | 0.059 / 0.088 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.161 | -0.200 to -0.124 | 0.0005 | -0.173 | -0.214 to -0.135 | 0.0005 | 0.0000 | yes | 0.055 / 0.081 |
| `gpt-oss-20b` | `glm-5.3` | -0.142 | -0.186 to -0.101 | 0.0005 | -0.153 | -0.196 to -0.112 | 0.0005 | 0.0000 | yes | 0.061 / 0.090 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.136 | 0.099 to 0.177 | 0.0005 | +0.144 | 0.105 to 0.185 | 0.0005 | 0.0000 | yes | 0.056 / 0.083 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.135 | -0.175 to -0.098 | 0.0005 | -0.142 | -0.183 to -0.105 | 0.0005 | 0.0000 | yes | 0.055 / 0.082 |
| `kimi-k3` | `gpt-5.4-mini` | +0.129 | 0.091 to 0.170 | 0.0005 | +0.136 | 0.097 to 0.176 | 0.0005 | 0.0000 | yes | 0.055 / 0.081 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.120 | 0.085 to 0.161 | 0.0005 | +0.128 | 0.091 to 0.167 | 0.0005 | 0.0000 | yes | 0.054 / 0.080 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.117 | 0.084 to 0.154 | 0.0005 | +0.125 | 0.091 to 0.164 | 0.0005 | 0.0000 | yes | 0.051 / 0.075 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.116 | -0.150 to -0.086 | 0.0005 | -0.126 | -0.161 to -0.094 | 0.0005 | 0.0000 | yes | 0.046 / 0.068 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.115 | -0.148 to -0.084 | 0.0005 | -0.124 | -0.160 to -0.093 | 0.0005 | 0.0000 | yes | 0.046 / 0.069 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.110 | 0.074 to 0.148 | 0.0005 | +0.114 | 0.079 to 0.152 | 0.0005 | 0.0000 | yes | 0.052 / 0.077 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.109 | -0.148 to -0.074 | 0.0005 | -0.112 | -0.151 to -0.078 | 0.0005 | 0.0000 | yes | 0.052 / 0.078 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.108 | -0.143 to -0.077 | 0.0005 | -0.118 | -0.152 to -0.085 | 0.0005 | 0.0000 | yes | 0.047 / 0.070 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.102 | 0.069 to 0.138 | 0.0005 | +0.106 | 0.071 to 0.144 | 0.0005 | 0.0000 | yes | 0.051 / 0.075 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.100 | -0.131 to -0.072 | 0.0005 | -0.110 | -0.143 to -0.080 | 0.0005 | 0.0000 | yes | 0.044 / 0.065 |
| `glm-5.3` | `gpt-5.4-mini` | +0.098 | 0.065 to 0.134 | 0.0005 | +0.106 | 0.071 to 0.143 | 0.0005 | 0.0000 | yes | 0.051 / 0.075 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.097 | -0.129 to -0.068 | 0.0005 | -0.108 | -0.141 to -0.076 | 0.0005 | 0.0000 | yes | 0.044 / 0.065 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.094 | 0.064 to 0.128 | 0.0005 | +0.098 | 0.067 to 0.134 | 0.0005 | 0.0000 | yes | 0.047 / 0.069 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.091 | 0.059 to 0.126 | 0.0005 | +0.096 | 0.062 to 0.132 | 0.0005 | 0.0000 | yes | 0.048 / 0.072 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.089 | 0.064 to 0.118 | 0.0005 | +0.109 | 0.076 to 0.145 | 0.0005 | 0.0000 | yes | 0.038 / 0.057 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.088 | -0.114 to -0.064 | 0.0005 | -0.108 | -0.142 to -0.075 | 0.0005 | 0.0000 | yes | 0.037 / 0.055 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.081 | -0.110 to -0.054 | 0.0005 | -0.101 | -0.137 to -0.067 | 0.0005 | 0.0000 | yes | 0.040 / 0.060 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.078 | -0.110 to -0.049 | 0.0005 | -0.088 | -0.120 to -0.059 | 0.0005 | 0.0000 | yes | 0.043 / 0.064 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.073 | -0.100 to -0.048 | 0.0005 | -0.093 | -0.128 to -0.061 | 0.0005 | 0.0000 | yes | 0.037 / 0.055 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.072 | 0.044 to 0.103 | 0.0005 | +0.076 | 0.048 to 0.108 | 0.0005 | 0.0000 | yes | 0.042 / 0.063 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.070 | -0.093 to -0.048 | 0.0005 | -0.091 | -0.124 to -0.060 | 0.0005 | 0.0000 | yes | 0.033 / 0.049 |
| `deepseek-v4.1-flash` | `glm-5.3` | +0.038 | 0.018 to 0.062 | 0.0005 | +0.038 | 0.017 to 0.062 | 0.0037 | 0.0000 | yes | 0.031 / 0.047 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.091 | -0.134 to -0.049 | 0.0010 | -0.082 | -0.130 to -0.035 | 0.0223 | 0.0000 | yes | 0.061 / 0.091 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.051 | -0.081 to -0.021 | 0.0211 | -0.071 | -0.108 to -0.035 | 0.0029 | 0.0000 | yes | 0.043 / 0.063 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.064 | -0.104 to -0.025 | 0.0315 | -0.065 | -0.106 to -0.027 | 0.0337 | 0.0000 | yes | 0.056 / 0.083 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.070 | -0.115 to -0.025 | 0.0436 | -0.077 | -0.123 to -0.032 | 0.0244 | 0.0000 | yes | 0.064 / 0.095 |
| `glm-5.3` | `claude-sonnet-5` | -0.037 | -0.065 to -0.014 | 0.0436 | -0.036 | -0.064 to -0.013 | 0.0722 | 0.0000 | no | 0.037 / 0.055 |
| `glm-5.3` | `kimi-k3` | -0.030 | -0.055 to -0.010 | 0.0436 | -0.030 | -0.054 to -0.009 | 0.1184 | 0.0000 | no | 0.031 / 0.047 |
| `glm-5.3-flash` | `glm-5.3` | +0.022 | 0.008 to 0.039 | 0.0436 | +0.022 | 0.008 to 0.039 | 0.0858 | 0.0000 | no | 0.022 / 0.033 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | +0.016 | 0.006 to 0.028 | 0.0436 | +0.016 | 0.005 to 0.028 | 0.1663 | 0.0005 | no | 0.016 / 0.023 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.047 | 0.014 to 0.082 | 0.1132 | +0.035 | -0.007 to 0.074 | 1.0000 | 0.0018 | yes | 0.049 / 0.072 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.019 | 0.005 to 0.035 | 0.1618 | +0.018 | 0.004 to 0.035 | 0.3988 | 0.0001 | yes | 0.022 / 0.033 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.018 | -0.036 to -0.002 | 0.5979 | -0.017 | -0.036 to 0.000 | 1.0000 | 0.0018 | yes | 0.025 / 0.037 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.044 | -0.090 to 0.001 | 0.8128 | -0.048 | -0.094 to -0.003 | 0.6720 | 0.0000 | yes | 0.064 / 0.095 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.026 | -0.055 to -0.000 | 0.8549 | -0.030 | -0.058 to -0.003 | 0.6563 | 0.0010 | yes | 0.040 / 0.059 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.027 | -0.058 to 0.003 | 0.9257 | -0.017 | -0.054 to 0.022 | 1.0000 | 0.4399 | yes | 0.043 / 0.064 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.015 | -0.033 to -0.001 | 0.9257 | -0.014 | -0.032 to 0.000 | 1.0000 | 0.0093 | yes | 0.023 / 0.034 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.021 | -0.009 to 0.053 | 1.0000 | +0.005 | -0.035 to 0.044 | 1.0000 | 1.0000 | yes | 0.044 / 0.065 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.020 | -0.011 to 0.052 | 1.0000 | +0.018 | -0.014 to 0.050 | 1.0000 | 0.3268 | yes | 0.044 / 0.066 |
| `glm-5.3` | `muse-glimmer-30b` | -0.019 | -0.044 to 0.004 | 1.0000 | -0.020 | -0.047 to 0.004 | 1.0000 | 0.0027 | yes | 0.035 / 0.053 |
| `muse-glimmer-30b` | `kimi-k3` | -0.011 | -0.032 to 0.008 | 1.0000 | -0.010 | -0.031 to 0.010 | 1.0000 | 0.3268 | yes | 0.029 / 0.043 |
| `glm-5.3-flash` | `kimi-k3` | -0.008 | -0.024 to 0.006 | 1.0000 | -0.008 | -0.023 to 0.006 | 1.0000 | 0.4487 | yes | 0.021 / 0.031 |
| `deepseek-v4.1-flash` | `kimi-k3` | +0.008 | -0.002 to 0.020 | 1.0000 | +0.008 | -0.002 to 0.020 | 1.0000 | 0.2734 | yes | 0.015 / 0.023 |
| `kimi-k3` | `claude-sonnet-5` | -0.007 | -0.020 to 0.004 | 1.0000 | -0.007 | -0.021 to 0.006 | 1.0000 | 0.5032 | yes | 0.017 / 0.026 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | -0.006 | -0.033 to 0.020 | 1.0000 | -0.012 | -0.039 to 0.015 | 1.0000 | 0.5354 | yes | 0.038 / 0.056 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.003 | -0.014 to 0.021 | 1.0000 | +0.003 | -0.014 to 0.021 | 1.0000 | 1.0000 | yes | 0.025 / 0.036 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | +0.001 | -0.008 to 0.010 | 1.0000 | +0.001 | -0.010 to 0.012 | 1.0000 | 1.0000 | yes | 0.013 / 0.019 |

*Added 2026-10-03 (D-173, next steps A9):* the last column is the smallest mean per-template difference a paired test on these 150 differences detects at 80% power, at two-sided 0.05 (2.80 x the SD of the differences / sqrt(150)) and at the strictest Holm step over the 55 pairs. A pair whose difference is below its detectable value is one this design could not have separated: 19 of the 55 pairs, 17 of them not significant.

## Q2. The complexity cliff: Easy minus Advanced

Answer score on the 58 Easy templates minus the 34 Advanced, templates resampled within each tier. **Corrected after the first results (D-146):** the planned test permuted the tier labels on the raw difference, which is liberal when the smaller tier has the larger spread, and Advanced template means spread two to four times as widely as Easy ones (the two SD columns); its count of models also moved with the seed. The count now rests on Welch's t-test, Holm across the eleven models: **2 of 11** hold, against 6 under the planned permutation, printed beside it. The detectable gap at 80% power is given as the plan defined it (2.80 x pooled sigma x sqrt(1/58 + 1/34)) and at the strictest Holm step from the Welch standard error.

| model | Easy | Advanced | gap | 95% CI | Welch p (Holm) | planned p (Holm) | SD Easy | SD Adv | detectable, planned | detectable, Holm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.900 | 0.719 | +0.181 | 0.070 to 0.299 | 0.0353 | 0.0094 | 0.207 | 0.303 | 0.149 | 0.216 |
| `gemma-4-26b-a4b` | 0.919 | 0.802 | +0.117 | 0.032 to 0.205 | 0.1029 | 0.0361 | 0.176 | 0.224 | 0.118 | 0.165 |
| `deepseek-v4.1-flash` | 0.995 | 0.979 | +0.016 | -0.001 to 0.041 | 0.6593 | 0.2327 | 0.018 | 0.063 | 0.024 | 0.040 |
| `qwen3-235b-a22b-2507` | 0.917 | 0.885 | +0.031 | -0.035 to 0.097 | 0.9200 | 1.0000 | 0.167 | 0.152 | 0.098 | 0.125 |
| `glm-5.3-flash` | 0.994 | 0.939 | +0.055 | 0.012 to 0.113 | 0.2759 | 0.0198 | 0.021 | 0.154 | 0.057 | 0.098 |
| `glm-5.3` | 0.984 | 0.865 | +0.120 | 0.030 to 0.223 | 0.1683 | 0.0198 | 0.054 | 0.293 | 0.110 | 0.187 |
| `muse-glimmer-30b` | 0.976 | 0.958 | +0.019 | -0.026 to 0.067 | 0.9200 | 1.0000 | 0.105 | 0.115 | 0.066 | 0.088 |
| `kimi-k3` | 0.992 | 0.953 | +0.039 | -0.002 to 0.094 | 0.6593 | 0.1700 | 0.022 | 0.146 | 0.055 | 0.093 |
| `gpt-5.4-mini` | 0.916 | 0.756 | +0.160 | 0.049 to 0.284 | 0.1029 | 0.0198 | 0.169 | 0.329 | 0.145 | 0.223 |
| `gemini-3.1-flash-lite` | 0.937 | 0.762 | +0.175 | 0.068 to 0.292 | 0.0458 | 0.0059 | 0.164 | 0.317 | 0.140 | 0.215 |
| `claude-sonnet-5` | 0.993 | 0.984 | +0.008 | -0.007 to 0.024 | 0.9200 | 1.0000 | 0.037 | 0.037 | 0.022 | 0.029 |

**The cliff under three variations.** Unusable rows left out of the template means, because an empty row at the output cap measures finishing within the ceiling as well as solving, and most such rows fall on Advanced templates; without the nine symbolic templates (D-138); and (*added 2026-10-03, D-171 and D-173*) without the two Advanced chemical templates whose question does not pin the answer to the check's tolerance, `work_isothermal_virial` and `adiabatic_flame_temperature`. Welch p, Holm across models.

| model | gap, as scored | gap, unusable left out | 95% CI | Welch p (Holm) | gap, no symbolic | 95% CI | Welch p (Holm) | gap, without the two chemical | 95% CI | Welch p (Holm) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | +0.181 | +0.170 | 0.063 to 0.284 | 0.0471 | +0.193 | 0.075 to 0.314 | 0.0375 | +0.147 | 0.041 to 0.258 | 0.1245 |
| `gemma-4-26b-a4b` | +0.117 | +0.117 | 0.032 to 0.205 | 0.1029 | +0.127 | 0.035 to 0.218 | 0.0912 | +0.088 | 0.007 to 0.167 | 0.3148 |
| `deepseek-v4.1-flash` | +0.016 | +0.002 | -0.005 to 0.011 | 1.0000 | +0.018 | -0.001 to 0.045 | 0.6283 | +0.005 | -0.005 to 0.016 | 1.0000 |
| `qwen3-235b-a22b-2507` | +0.031 | +0.031 | -0.035 to 0.097 | 1.0000 | +0.025 | -0.043 to 0.094 | 1.0000 | +0.024 | -0.043 to 0.092 | 1.0000 |
| `glm-5.3-flash` | +0.055 | +0.007 | -0.005 to 0.024 | 1.0000 | +0.055 | 0.009 to 0.117 | 0.3819 | +0.030 | 0.001 to 0.067 | 0.5850 |
| `glm-5.3` | +0.120 | -0.004 | -0.016 to 0.007 | 1.0000 | +0.126 | 0.028 to 0.240 | 0.2057 | +0.074 | 0.004 to 0.161 | 0.5850 |
| `muse-glimmer-30b` | +0.019 | -0.007 | -0.040 to 0.027 | 1.0000 | +0.020 | -0.027 to 0.073 | 1.0000 | +0.017 | -0.028 to 0.067 | 1.0000 |
| `kimi-k3` | +0.039 | +0.016 | -0.006 to 0.046 | 1.0000 | +0.043 | -0.001 to 0.103 | 0.6283 | +0.011 | -0.008 to 0.040 | 1.0000 |
| `gpt-5.4-mini` | +0.160 | +0.160 | 0.049 to 0.284 | 0.1029 | +0.175 | 0.055 to 0.303 | 0.0912 | +0.138 | 0.031 to 0.261 | 0.2221 |
| `gemini-3.1-flash-lite` | +0.175 | +0.175 | 0.068 to 0.292 | 0.0471 | +0.174 | 0.057 to 0.297 | 0.0816 | +0.140 | 0.037 to 0.250 | 0.1593 |
| `claude-sonnet-5` | +0.008 | +0.004 | -0.011 to 0.019 | 1.0000 | +0.005 | -0.010 to 0.020 | 1.0000 | +0.005 | -0.010 to 0.020 | 1.0000 |

Without the two chemical templates the gap holds for 0 of 11 models under Welch with Holm.

## Q3. What the process scores add beyond the answer

Wrong-answer traces score 0 (incorrect or unusable); fully solved ones score 1; a partial answer is in neither. E3 coverage leaves out the 57 items with no milestones (all 15 of 2 templates and some items of 13 more); 52 templates have an item with at most one milestone (47 with exactly one), where coverage is close to an answer check. E3 is the deterministic part of E5, which adds the judge's verdict on the milestones E3 does not find. The floor is the same trace scored against a sibling item's milestones (D-149): what coverage a trace reaches by chance, per model, on the readable wrong answers. The digit rule's flag counts a trace when any step is flagged. Against the experts' step labels on the pilot's fully solved traces it has precision 0.817 and recall 0.427 (0.750 and 0.320 before D-156; 0.800 and 0.427 between D-156 and D-159, unchanged by D-160; SCORER_VALIDATION.md). *Added 2026-09-30 and 2026-10-01 (D-154 to D-165), after the first results were read:* on this roster a domain expert found 111 of 220 sampled flags real before the first fix, 0.505 (FLAG_REVIEW.md), and 155 of 206 real after it, 0.752 (FLAG_REVIEW_2.md). The second fix, built from those notes, took the flag off 42 of their 51 misread steps and left it on 154 of the 155 slips' steps (DIGIT_FIX_2.md). Two review agents then checked both fixes; the corrections they led to, with LaTeX control spaces now read, took the flag off 82 full-run steps and put it on 370, and kept every slip's step the second fix had kept (D-160, DIGIT_FIX_3.md). The same expert then read 190 flags of the rule as it now stands: 171 of the 189 decided are real, 0.905 (95% Wilson 0.854 to 0.939; per model 0.632 to 1.000; FLAG_REVIEW_3.md), and by the owner's rule no fourth reading follows. So its rate is still not a count of slips; and it reads a different amount of arithmetic in each model's traces (claims checked per answered trace, shown), so a low rate can mean little was read.

**Milestone coverage over all traces** (*added 2026-10-02, D-170, descriptive and outside the plan's tests*). Per model, the share of the gold's milestones a trace states (E3), or states or the judge rules REACHED (E5-strict), averaged over every trace whose item has milestones, an unusable trace scoring what it reached and an empty one nothing, and over the readable traces alone; template-level intervals. Coverage measures progress through the gold derivation, not the absence of error: on the pilot's correct-answer traces with a flawed step the milestone evaluators scored below chance (RESULTS_X1 Finding 5). The last two columns repeat, from the tables below, the share of fully solved traces the digit rule flags and the share the router's judge flags.

| model | traces with milestones | E3, all | 95% CI | E3, readable | 95% CI | E5-strict, all | 95% CI | E5-strict, readable | 95% CI | digit rule, fully solved | router judge, fully solved |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 2193 | 0.790 | 0.751 to 0.827 | 0.811 | 0.775 to 0.844 | 0.812 | 0.775 to 0.848 | 0.833 | 0.800 to 0.865 | 0.147 | 0.206 |
| `gemma-4-26b-a4b` | 2193 | 0.855 | 0.824 to 0.884 | 0.855 | 0.824 to 0.883 | 0.875 | 0.846 to 0.902 | 0.875 | 0.847 to 0.902 | 0.132 | 0.146 |
| `deepseek-v4.1-flash` | 2193 | 0.860 | 0.830 to 0.889 | 0.864 | 0.835 to 0.892 | 0.898 | 0.873 to 0.922 | 0.901 | 0.876 to 0.924 | 0.006 | 0.010 |
| `qwen3-235b-a22b-2507` | 2193 | 0.866 | 0.835 to 0.894 | 0.866 | 0.836 to 0.894 | 0.892 | 0.866 to 0.915 | 0.892 | 0.867 to 0.916 | 0.232 | 0.107 |
| `glm-5.3-flash` | 2193 | 0.891 | 0.863 to 0.917 | 0.909 | 0.886 to 0.930 | 0.911 | 0.885 to 0.935 | 0.930 | 0.909 to 0.950 | 0.019 | 0.032 |
| `glm-5.3` | 2193 | 0.880 | 0.845 to 0.912 | 0.922 | 0.899 to 0.942 | 0.899 | 0.865 to 0.930 | 0.941 | 0.921 to 0.960 | 0.018 | 0.026 |
| `muse-glimmer-30b` | 2193 | 0.862 | 0.832 to 0.890 | 0.874 | 0.848 to 0.900 | 0.897 | 0.870 to 0.922 | 0.910 | 0.887 to 0.932 | 0.058 | 0.059 |
| `kimi-k3` | 2193 | 0.877 | 0.849 to 0.902 | 0.883 | 0.858 to 0.906 | 0.911 | 0.887 to 0.933 | 0.918 | 0.896 to 0.938 | 0.007 | 0.017 |
| `gpt-5.4-mini` | 2193 | 0.799 | 0.760 to 0.836 | 0.799 | 0.760 to 0.835 | 0.835 | 0.800 to 0.870 | 0.835 | 0.799 to 0.869 | 0.141 | 0.219 |
| `gemini-3.1-flash-lite` | 2193 | 0.840 | 0.809 to 0.870 | 0.840 | 0.809 to 0.871 | 0.865 | 0.836 to 0.893 | 0.865 | 0.836 to 0.893 | 0.070 | 0.122 |
| `claude-sonnet-5` | 2193 | 0.903 | 0.880 to 0.925 | 0.904 | 0.880 to 0.925 | 0.923 | 0.902 to 0.942 | 0.923 | 0.903 to 0.942 | 0.143 | 0.054 |

**Coverage compared across models** (*exploratory; added 2026-10-03, D-173, next steps A1 and A9: the comparison the July rebuttal promised, Wilcoxon on the continuous reasoning score beside McNemar on the answer*). Per model, E5-strict coverage as the mean of its template means over the 148 templates every model has coverage on (an unusable trace scores what it reached), with a template bootstrap, the per-template SD within and between as Q1 has them, and the model's rank on the answer score beside its rank on coverage. Below, the 55 pairs: the paired template bootstrap of the difference, the sign-flip permutation test and the Wilcoxon signed-rank test on the per-template differences, each Holm-corrected over the 55 pairs, and the smallest difference the design detects at 80% power. Coverage credits stated intermediates, so a terser model scores lower without reasoning worse: the comparison ranks what traces state, not how soundly they reason.

| model | coverage | 95% CI | SD within | SD between | rank on answers | rank on coverage |
|---|---:|---:|---:|---:|---:|---:|
| `claude-sonnet-5` | 0.922 | 0.902 to 0.942 | 0.067 | 0.125 | 2 | 1 |
| `glm-5.3-flash` | 0.912 | 0.886 to 0.936 | 0.077 | 0.158 | 4 | 2 |
| `kimi-k3` | 0.912 | 0.888 to 0.934 | 0.081 | 0.143 | 3 | 3 |
| `glm-5.3` | 0.900 | 0.866 to 0.930 | 0.076 | 0.199 | 6 | 4 |
| `deepseek-v4.1-flash` | 0.899 | 0.874 to 0.922 | 0.076 | 0.150 | 1 | 5 |
| `muse-glimmer-30b` | 0.897 | 0.871 to 0.922 | 0.082 | 0.160 | 5 | 6 |
| `qwen3-235b-a22b-2507` | 0.892 | 0.867 to 0.916 | 0.100 | 0.152 | 7 | 7 |
| `gemma-4-26b-a4b` | 0.876 | 0.847 to 0.903 | 0.106 | 0.172 | 9 | 8 |
| `gemini-3.1-flash-lite` | 0.864 | 0.835 to 0.892 | 0.104 | 0.179 | 8 | 9 |
| `gpt-5.4-mini` | 0.834 | 0.798 to 0.868 | 0.111 | 0.223 | 10 | 10 |
| `gpt-oss-20b` | 0.812 | 0.775 to 0.847 | 0.146 | 0.226 | 11 | 11 |

Of the 55 pairs, 28 differ at a Holm-adjusted p below 0.05 under the sign-flip test and 29 under Wilcoxon. Kendall's tau between the models' answer-score order and their coverage order: 0.782 (95% CI 0.564 to 0.891, templates resampled).

| a | b | a - b | 95% CI | sign-flip p (Holm) | Wilcoxon p (Holm) | detectable, 0.05 / Holm |
|---|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | `claude-sonnet-5` | -0.110 | -0.141 to -0.080 | 0.0005 | 0.0000 | 0.043 / 0.064 |
| `gpt-oss-20b` | `glm-5.3-flash` | -0.100 | -0.132 to -0.068 | 0.0005 | 0.0000 | 0.047 / 0.069 |
| `gpt-oss-20b` | `kimi-k3` | -0.100 | -0.130 to -0.070 | 0.0005 | 0.0000 | 0.045 / 0.066 |
| `gpt-5.4-mini` | `claude-sonnet-5` | -0.088 | -0.117 to -0.062 | 0.0005 | 0.0000 | 0.039 / 0.058 |
| `gpt-oss-20b` | `glm-5.3` | -0.088 | -0.126 to -0.050 | 0.0005 | 0.0000 | 0.053 / 0.078 |
| `gpt-oss-20b` | `deepseek-v4.1-flash` | -0.087 | -0.117 to -0.057 | 0.0005 | 0.0000 | 0.043 / 0.064 |
| `gpt-oss-20b` | `muse-glimmer-30b` | -0.085 | -0.114 to -0.057 | 0.0005 | 0.0000 | 0.041 / 0.061 |
| `gpt-oss-20b` | `qwen3-235b-a22b-2507` | -0.080 | -0.111 to -0.050 | 0.0005 | 0.0001 | 0.044 / 0.065 |
| `glm-5.3-flash` | `gpt-5.4-mini` | +0.078 | 0.051 to 0.109 | 0.0005 | 0.0000 | 0.042 / 0.062 |
| `kimi-k3` | `gpt-5.4-mini` | +0.078 | 0.053 to 0.107 | 0.0005 | 0.0000 | 0.038 / 0.057 |
| `deepseek-v4.1-flash` | `gpt-5.4-mini` | +0.065 | 0.041 to 0.092 | 0.0005 | 0.0000 | 0.037 / 0.055 |
| `muse-glimmer-30b` | `gpt-5.4-mini` | +0.064 | 0.040 to 0.088 | 0.0005 | 0.0000 | 0.034 / 0.051 |
| `gpt-oss-20b` | `gemma-4-26b-a4b` | -0.064 | -0.092 to -0.037 | 0.0005 | 0.0027 | 0.040 / 0.059 |
| `gemini-3.1-flash-lite` | `claude-sonnet-5` | -0.058 | -0.080 to -0.036 | 0.0005 | 0.0000 | 0.031 / 0.046 |
| `kimi-k3` | `gemini-3.1-flash-lite` | +0.047 | 0.028 to 0.068 | 0.0005 | 0.0004 | 0.028 / 0.042 |
| `qwen3-235b-a22b-2507` | `gpt-5.4-mini` | +0.059 | 0.035 to 0.084 | 0.0008 | 0.0002 | 0.036 / 0.053 |
| `glm-5.3-flash` | `gemini-3.1-flash-lite` | +0.047 | 0.027 to 0.069 | 0.0008 | 0.0002 | 0.030 / 0.045 |
| `gemma-4-26b-a4b` | `claude-sonnet-5` | -0.046 | -0.068 to -0.026 | 0.0008 | 0.0004 | 0.030 / 0.045 |
| `glm-5.3` | `gpt-5.4-mini` | +0.066 | 0.035 to 0.100 | 0.0018 | 0.0000 | 0.046 / 0.069 |
| `deepseek-v4.1-flash` | `gemini-3.1-flash-lite` | +0.034 | 0.016 to 0.055 | 0.0122 | 0.0247 | 0.028 / 0.041 |
| `gemma-4-26b-a4b` | `gpt-5.4-mini` | +0.042 | 0.020 to 0.068 | 0.0133 | 0.0093 | 0.034 / 0.051 |
| `gemma-4-26b-a4b` | `kimi-k3` | -0.036 | -0.058 to -0.016 | 0.0136 | 0.0213 | 0.030 / 0.044 |
| `gpt-oss-20b` | `gemini-3.1-flash-lite` | -0.052 | -0.082 to -0.023 | 0.0231 | 0.0183 | 0.043 / 0.063 |
| `gemma-4-26b-a4b` | `glm-5.3-flash` | -0.036 | -0.058 to -0.015 | 0.0246 | 0.0084 | 0.031 / 0.046 |
| `qwen3-235b-a22b-2507` | `claude-sonnet-5` | -0.030 | -0.049 to -0.013 | 0.0246 | 0.0216 | 0.025 / 0.038 |
| `muse-glimmer-30b` | `gemini-3.1-flash-lite` | +0.033 | 0.014 to 0.052 | 0.0330 | 0.0216 | 0.028 / 0.042 |
| `qwen3-235b-a22b-2507` | `gemini-3.1-flash-lite` | +0.028 | 0.011 to 0.045 | 0.0380 | 0.0783 | 0.024 / 0.036 |
| `deepseek-v4.1-flash` | `claude-sonnet-5` | -0.023 | -0.038 to -0.009 | 0.0493 | 0.0234 | 0.021 / 0.031 |
| `gpt-5.4-mini` | `gemini-3.1-flash-lite` | -0.031 | -0.053 to -0.010 | 0.1320 | 0.1823 | 0.031 / 0.046 |
| `glm-5.3` | `gemini-3.1-flash-lite` | +0.036 | 0.008 to 0.063 | 0.2792 | 0.0123 | 0.039 / 0.058 |
| `muse-glimmer-30b` | `claude-sonnet-5` | -0.025 | -0.045 to -0.007 | 0.2792 | 0.4366 | 0.028 / 0.041 |
| `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | -0.023 | -0.043 to -0.005 | 0.4716 | 1.0000 | 0.028 / 0.041 |
| `gemma-4-26b-a4b` | `muse-glimmer-30b` | -0.021 | -0.041 to -0.002 | 0.7516 | 0.2718 | 0.028 / 0.042 |
| `qwen3-235b-a22b-2507` | `kimi-k3` | -0.020 | -0.038 to -0.002 | 0.7516 | 0.0726 | 0.026 / 0.038 |
| `deepseek-v4.1-flash` | `kimi-k3` | -0.013 | -0.025 to -0.001 | 0.8413 | 0.4198 | 0.017 / 0.026 |
| `gemma-4-26b-a4b` | `qwen3-235b-a22b-2507` | -0.016 | -0.032 to -0.001 | 0.8812 | 0.4406 | 0.023 / 0.033 |
| `gemma-4-26b-a4b` | `glm-5.3` | -0.024 | -0.052 to 0.005 | 1.0000 | 0.0221 | 0.041 / 0.060 |
| `glm-5.3` | `claude-sonnet-5` | -0.022 | -0.052 to 0.005 | 1.0000 | 1.0000 | 0.041 / 0.061 |
| `gpt-oss-20b` | `gpt-5.4-mini` | -0.021 | -0.054 to 0.010 | 1.0000 | 1.0000 | 0.045 / 0.067 |
| `qwen3-235b-a22b-2507` | `glm-5.3-flash` | -0.020 | -0.041 to 0.002 | 1.0000 | 0.0726 | 0.030 / 0.045 |
| `glm-5.3-flash` | `muse-glimmer-30b` | +0.015 | -0.002 to 0.031 | 1.0000 | 1.0000 | 0.024 / 0.035 |
| `muse-glimmer-30b` | `kimi-k3` | -0.014 | -0.035 to 0.004 | 1.0000 | 1.0000 | 0.027 / 0.041 |
| `deepseek-v4.1-flash` | `glm-5.3-flash` | -0.013 | -0.030 to 0.006 | 1.0000 | 0.3688 | 0.026 / 0.038 |
| `glm-5.3-flash` | `glm-5.3` | +0.012 | -0.003 to 0.029 | 1.0000 | 1.0000 | 0.023 / 0.034 |
| `glm-5.3` | `kimi-k3` | -0.012 | -0.039 to 0.012 | 1.0000 | 1.0000 | 0.037 / 0.055 |
| `gemma-4-26b-a4b` | `gemini-3.1-flash-lite` | +0.011 | -0.002 to 0.025 | 1.0000 | 1.0000 | 0.019 / 0.029 |
| `kimi-k3` | `claude-sonnet-5` | -0.010 | -0.025 to 0.003 | 1.0000 | 1.0000 | 0.020 / 0.030 |
| `glm-5.3-flash` | `claude-sonnet-5` | -0.010 | -0.031 to 0.009 | 1.0000 | 1.0000 | 0.029 / 0.043 |
| `qwen3-235b-a22b-2507` | `glm-5.3` | -0.008 | -0.035 to 0.022 | 1.0000 | 0.1271 | 0.042 / 0.063 |
| `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | +0.007 | -0.011 to 0.026 | 1.0000 | 1.0000 | 0.027 / 0.040 |
| `qwen3-235b-a22b-2507` | `muse-glimmer-30b` | -0.005 | -0.022 to 0.013 | 1.0000 | 1.0000 | 0.025 / 0.037 |
| `glm-5.3` | `muse-glimmer-30b` | +0.003 | -0.021 to 0.024 | 1.0000 | 1.0000 | 0.032 / 0.047 |
| `deepseek-v4.1-flash` | `muse-glimmer-30b` | +0.002 | -0.015 to 0.021 | 1.0000 | 1.0000 | 0.026 / 0.039 |
| `deepseek-v4.1-flash` | `glm-5.3` | -0.001 | -0.027 to 0.026 | 1.0000 | 1.0000 | 0.038 / 0.056 |
| `glm-5.3-flash` | `kimi-k3` | +0.000 | -0.019 to 0.018 | 1.0000 | 1.0000 | 0.027 / 0.039 |

**Does coverage track verbosity?** (*exploratory; D-173, next steps A10*). Per model, Spearman's rho across templates between the template's mean coverage and its mean number of steps and of arithmetic claims per readable trace; and coverage on the fully solved traces against the answered wrong answers, readable traces, template intervals.

| model | rho(coverage, steps) | rho(coverage, claims) | median steps | coverage, fully solved | 95% CI | coverage, wrong answers | 95% CI | wrong answers |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | -0.262 | -0.137 | 9 | 0.907 | 0.883 to 0.929 | 0.465 | 0.400 to 0.540 | 335 |
| `gemma-4-26b-a4b` | -0.210 | -0.294 | 5 | 0.915 | 0.891 to 0.937 | 0.602 | 0.524 to 0.686 | 260 |
| `deepseek-v4.1-flash` | -0.239 | -0.143 | 6 | 0.903 | 0.878 to 0.925 | 0.667 | 0.293 to 0.811 | 12 |
| `qwen3-235b-a22b-2507` | -0.189 | -0.183 | 10 | 0.914 | 0.889 to 0.937 | 0.616 | 0.538 to 0.699 | 173 |
| `glm-5.3-flash` | -0.095 | -0.211 | 7 | 0.933 | 0.913 to 0.952 | 0.458 | 0.247 to 0.667 | 12 |
| `glm-5.3` | -0.147 | -0.219 | 8 | 0.943 | 0.922 to 0.962 | 0.596 | 0.275 to 0.896 | 8 |
| `muse-glimmer-30b` | -0.208 | -0.148 | 5 | 0.914 | 0.892 to 0.935 | 0.653 | 0.546 to 0.858 | 33 |
| `kimi-k3` | -0.213 | -0.129 | 6 | 0.921 | 0.899 to 0.942 | 0.584 | 0.411 to 0.745 | 20 |
| `gpt-5.4-mini` | -0.222 | -0.168 | 6 | 0.897 | 0.869 to 0.922 | 0.485 | 0.410 to 0.572 | 311 |
| `gemini-3.1-flash-lite` | -0.206 | -0.248 | 5 | 0.912 | 0.888 to 0.933 | 0.537 | 0.467 to 0.609 | 257 |
| `claude-sonnet-5` | -0.183 | -0.233 | 7 | 0.926 | 0.905 to 0.945 | 0.869 | 0.774 to 0.966 | 20 |

**Verdict against coverage at the trace level** (*exploratory; D-173, next steps A11*). The share of fully solved traces with coverage below 0.5, and of answered wrong answers with coverage 1.0, template intervals. The second is what a process score adds on a wrong answer: a complete derivation to a wrong value.

| model | fully solved traces (items with milestones) | with coverage < 0.5 | 95% CI | answered wrong answers (items with milestones) | with coverage 1.0 | 95% CI |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 1742 | 0.033 | 0.017 to 0.052 | 335 | 0.146 | 0.088 to 0.221 |
| `gemma-4-26b-a4b` | 1879 | 0.037 | 0.019 to 0.060 | 260 | 0.258 | 0.169 to 0.364 |
| `deepseek-v4.1-flash` | 2164 | 0.034 | 0.017 to 0.054 | 12 | 0.083 | 0.000 to 0.429 |
| `qwen3-235b-a22b-2507` | 1919 | 0.036 | 0.017 to 0.060 | 173 | 0.249 | 0.155 to 0.368 |
| `glm-5.3-flash` | 2126 | 0.023 | 0.008 to 0.042 | 12 | 0.083 | 0.000 to 0.273 |
| `glm-5.3` | 2075 | 0.020 | 0.004 to 0.040 | 8 | 0.500 | 0.125 to 0.875 |
| `muse-glimmer-30b` | 2123 | 0.030 | 0.016 to 0.047 | 33 | 0.364 | 0.167 to 0.750 |
| `kimi-k3` | 2145 | 0.028 | 0.012 to 0.047 | 20 | 0.200 | 0.000 to 0.444 |
| `gpt-5.4-mini` | 1840 | 0.049 | 0.028 to 0.074 | 311 | 0.164 | 0.096 to 0.253 |
| `gemini-3.1-flash-lite` | 1906 | 0.029 | 0.014 to 0.048 | 257 | 0.140 | 0.076 to 0.217 |
| `claude-sonnet-5` | 2157 | 0.025 | 0.012 to 0.040 | 20 | 0.550 | 0.250 to 0.875 |

**Milestones on the wrong-answer traces.** An unusable trace reaches only what it wrote before it stopped, and an empty one nothing, so coverage is also shown on the readable wrong answers alone.

| model | wrong-answer traces | of them unusable | E3 coverage | 95% CI | readable only | 95% CI | floor, readable | E5 coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 401 | 59 | 0.390 | 0.326 to 0.463 | 0.455 | 0.389 to 0.531 | 0.117 | 0.402 (0.340 to 0.474); readable 0.465 |
| `gemma-4-26b-a4b` | 260 | 0 | 0.588 | 0.510 to 0.672 | 0.588 | 0.510 to 0.669 | 0.145 | 0.602 (0.524 to 0.683); readable 0.602 |
| `deepseek-v4.1-flash` | 20 | 8 | 0.400 | 0.083 to 0.661 | 0.667 | 0.293 to 0.807 | 0.141 | 0.400 (0.083 to 0.661); readable 0.667 |
| `qwen3-235b-a22b-2507` | 176 | 0 | 0.586 | 0.504 to 0.676 | 0.586 | 0.504 to 0.676 | 0.175 | 0.616 (0.538 to 0.699); readable 0.616 |
| `glm-5.3-flash` | 57 | 45 | 0.080 | 0.030 to 0.164 | 0.378 | 0.156 to 0.606 | 0.134 | 0.096 (0.047 to 0.180); readable 0.458 |
| `glm-5.3` | 107 | 99 | 0.036 | 0.008 to 0.096 | 0.486 | 0.188 to 0.781 | 0.146 | 0.045 (0.012 to 0.109); readable 0.596 |
| `muse-glimmer-30b` | 65 | 31 | 0.300 | 0.170 to 0.489 | 0.583 | 0.479 to 0.744 | 0.178 | 0.337 (0.202 to 0.549); readable 0.653 |
| `kimi-k3` | 37 | 16 | 0.329 | 0.180 to 0.633 | 0.584 | 0.411 to 0.752 | 0.113 | 0.329 (0.179 to 0.644); readable 0.584 |
| `gpt-5.4-mini` | 311 | 0 | 0.442 | 0.368 to 0.522 | 0.442 | 0.369 to 0.524 | 0.117 | 0.485 (0.409 to 0.573); readable 0.485 |
| `gemini-3.1-flash-lite` | 259 | 0 | 0.511 | 0.443 to 0.580 | 0.511 | 0.441 to 0.582 | 0.112 | 0.537 (0.466 to 0.610); readable 0.537 |
| `claude-sonnet-5` | 22 | 2 | 0.773 | 0.603 to 0.906 | 0.850 | 0.738 to 0.961 | 0.150 | 0.790 (0.631 to 0.918); readable 0.869 |

**E5's judge**, MiMo-V2.5-Pro on the milestones E3 did not find (`judge.py`). E5 coverage is E5-strict: E3's milestones plus those the judge rules REACHED. On the pilot the judge never called a fabricated value REACHED and found 76% of true ones, so the score is conservative (RESULTS_E5). Judged fraction: milestones sent to the judge over those required, on answered traces; unjudged: milestones the reply did not name, over those sent. A trace whose call got no reply is left out of every rate and counted (D-148).

| model | calls | without a reply | judged fraction | of the judged, REACHED | unjudged |
|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 863 | 0 | 0.238 | 0.106 | 0.002 |
| `gemma-4-26b-a4b` | 747 | 0 | 0.187 | 0.124 | 0.000 |
| `deepseek-v4.1-flash` | 791 | 0 | 0.169 | 0.253 | 0.004 |
| `qwen3-235b-a22b-2507` | 729 | 0 | 0.158 | 0.163 | 0.000 |
| `glm-5.3-flash` | 582 | 0 | 0.126 | 0.215 | 0.001 |
| `glm-5.3` | 491 | 0 | 0.107 | 0.235 | 0.000 |
| `muse-glimmer-30b` | 730 | 0 | 0.167 | 0.275 | 0.002 |
| `kimi-k3` | 721 | 0 | 0.150 | 0.276 | 0.001 |
| `gpt-5.4-mini` | 934 | 0 | 0.250 | 0.140 | 0.000 |
| `gemini-3.1-flash-lite` | 827 | 0 | 0.207 | 0.150 | 0.002 |
| `claude-sonnet-5` | 620 | 0 | 0.125 | 0.218 | 0.005 |

**The digit rule.** The wrong-answer rate is over the answered wrong answers (D-149: an empty trace has no step to flag). Beside it, the 1% tolerance E4 shipped with, blind to most slips (RESULTS_X1), and where the first flag falls in a fully solved trace (0 = first step, 1 = last).

| model | wrong answers, answered | flag rate | 95% CI | fully solved traces | flag rate | 95% CI | at 1% | claims checked per trace | traces with a claim | first flag: traces, median position |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 345 | 0.229 | 0.158 to 0.305 | 1787 | 0.147 | 0.116 to 0.181 | 0.021 | 3.77 | 0.782 | 263, 0.667 |
| `gemma-4-26b-a4b` | 260 | 0.473 | 0.349 to 0.597 | 1934 | 0.132 | 0.099 to 0.168 | 0.017 | 4.56 | 0.689 | 256, 0.558 |
| `deepseek-v4.1-flash` | 12 | 0.000 | 0.000 to 0.000 | 2217 | 0.006 | 0.003 to 0.010 | 0.006 | 3.27 | 0.805 | 13, 0.667 |
| `qwen3-235b-a22b-2507` | 176 | 0.449 | 0.318 to 0.589 | 1972 | 0.232 | 0.187 to 0.279 | 0.054 | 10.21 | 0.907 | 457, 0.600 |
| `glm-5.3-flash` | 12 | 0.250 | 0.000 to 0.545 | 2182 | 0.019 | 0.012 to 0.026 | 0.014 | 4.74 | 0.852 | 41, 0.750 |
| `glm-5.3` | 8 | 0.125 | 0.000 to 0.375 | 2132 | 0.018 | 0.012 to 0.025 | 0.017 | 5.72 | 0.915 | 39, 0.727 |
| `muse-glimmer-30b` | 34 | 0.029 | 0.000 to 0.075 | 2176 | 0.058 | 0.044 to 0.074 | 0.019 | 3.27 | 0.708 | 127, 0.571 |
| `kimi-k3` | 22 | 0.000 | 0.000 to 0.000 | 2199 | 0.007 | 0.003 to 0.011 | 0.004 | 2.79 | 0.721 | 15, 0.571 |
| `gpt-5.4-mini` | 311 | 0.267 | 0.172 to 0.370 | 1894 | 0.141 | 0.108 to 0.178 | 0.012 | 3.78 | 0.797 | 268, 0.620 |
| `gemini-3.1-flash-lite` | 259 | 0.363 | 0.211 to 0.504 | 1961 | 0.070 | 0.048 to 0.095 | 0.016 | 3.87 | 0.757 | 137, 0.500 |
| `claude-sonnet-5` | 20 | 0.200 | 0.045 to 0.412 | 2214 | 0.143 | 0.111 to 0.177 | 0.019 | 5.32 | 0.857 | 316, 0.667 |

**Wrong-answer rate against the item's milestone count** (D-149): how failure grows with the depth of the gold derivation.

| model | 0 (57 items) | 1 (293 items) | 2 (462 items) | 3 (368 items) | 4-5 (542 items) | 6+ (528 items) |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.158 | 0.061 | 0.091 | 0.190 | 0.236 | 0.254 |
| `gemma-4-26b-a4b` | 0.000 | 0.044 | 0.080 | 0.106 | 0.094 | 0.227 |
| `deepseek-v4.1-flash` | 0.000 | 0.000 | 0.000 | 0.000 | 0.006 | 0.032 |
| `qwen3-235b-a22b-2507` | 0.053 | 0.034 | 0.056 | 0.073 | 0.085 | 0.121 |
| `glm-5.3-flash` | 0.000 | 0.003 | 0.017 | 0.011 | 0.022 | 0.061 |
| `glm-5.3` | 0.000 | 0.003 | 0.013 | 0.027 | 0.057 | 0.112 |
| `muse-glimmer-30b` | 0.018 | 0.003 | 0.017 | 0.038 | 0.030 | 0.047 |
| `kimi-k3` | 0.018 | 0.000 | 0.006 | 0.014 | 0.028 | 0.025 |
| `gpt-5.4-mini` | 0.000 | 0.020 | 0.076 | 0.095 | 0.170 | 0.271 |
| `gemini-3.1-flash-lite` | 0.035 | 0.024 | 0.048 | 0.073 | 0.149 | 0.227 |
| `claude-sonnet-5` | 0.000 | 0.000 | 0.017 | 0.003 | 0.013 | 0.011 |

**The step router** (`router.py`): the digit rule's flags, and MiMo-V2.5-Pro on every other step in one batched call per trace. On the pilot's labelled traces it had precision 0.707 and recall 0.603 over all steps, and 0.703 and 0.360 inside correct-answer traces (ROUTER_VALIDATION.md), so its rates are flags, not counts of errors. A trace counts as flagged when any step is; a trace whose call got no reply is left out and counted.

| model | calls | without a reply | flagged, fully solved | 95% CI | by the judge | 95% CI | flagged, wrong answers | 95% CI | steps flagged per trace | unjudged steps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 2194 | 0 | 0.311 | 0.272 to 0.349 | 0.206 | 0.177 to 0.237 | 0.916 | 0.871 to 0.950 | 1.01 | 0.001 |
| `gemma-4-26b-a4b` | 2249 | 0 | 0.248 | 0.204 to 0.295 | 0.146 | 0.112 to 0.182 | 0.908 | 0.834 to 0.960 | 0.66 | 0.000 |
| `deepseek-v4.1-flash` | 2242 | 0 | 0.015 | 0.010 to 0.022 | 0.010 | 0.005 to 0.015 | 0.417 | 0.000 to 0.636 | 0.03 | 0.000 |
| `qwen3-235b-a22b-2507` | 2250 | 1 | 0.297 | 0.251 to 0.346 | 0.107 | 0.084 to 0.131 | 0.864 | 0.776 to 0.930 | 0.86 | 0.001 |
| `glm-5.3-flash` | 2205 | 0 | 0.049 | 0.038 to 0.062 | 0.032 | 0.023 to 0.042 | 0.333 | 0.077 to 0.667 | 0.07 | 0.000 |
| `glm-5.3` | 2151 | 0 | 0.042 | 0.030 to 0.054 | 0.026 | 0.016 to 0.036 | 0.250 | 0.000 to 0.625 | 0.06 | 0.000 |
| `muse-glimmer-30b` | 2219 | 0 | 0.113 | 0.093 to 0.135 | 0.059 | 0.046 to 0.074 | 0.529 | 0.208 to 0.743 | 0.15 | 0.000 |
| `kimi-k3` | 2235 | 0 | 0.024 | 0.016 to 0.031 | 0.017 | 0.011 to 0.024 | 0.545 | 0.350 to 0.765 | 0.05 | 0.000 |
| `gpt-5.4-mini` | 2250 | 0 | 0.313 | 0.265 to 0.364 | 0.219 | 0.176 to 0.265 | 0.923 | 0.883 to 0.953 | 0.80 | 0.000 |
| `gemini-3.1-flash-lite` | 2250 | 0 | 0.176 | 0.138 to 0.217 | 0.122 | 0.091 to 0.158 | 0.907 | 0.831 to 0.960 | 0.50 | 0.000 |
| `claude-sonnet-5` | 2248 | 0 | 0.186 | 0.154 to 0.220 | 0.054 | 0.041 to 0.068 | 0.650 | 0.412 to 0.875 | 0.28 | 0.000 |

**Also reported, not tested: what points at a wrong answer.** On the answered traces that score 0, the share with a digit-rule flag, with a milestone E5 rules MISSING, and with a step the router's judge flags. A trace can be in several columns, or in none.

| model | answered wrong-answer traces | digit rule | E5 MISSING | router judge |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 345 | 0.229 | 0.777 | 0.881 |
| `gemma-4-26b-a4b` | 260 | 0.473 | 0.635 | 0.762 |
| `deepseek-v4.1-flash` | 12 | 0.000 | 0.833 | 0.417 |
| `qwen3-235b-a22b-2507` | 176 | 0.449 | 0.631 | 0.722 |
| `glm-5.3-flash` | 12 | 0.250 | 0.750 | 0.250 |
| `glm-5.3` | 8 | 0.125 | 0.500 | 0.125 |
| `muse-glimmer-30b` | 34 | 0.029 | 0.588 | 0.500 |
| `kimi-k3` | 22 | 0.000 | 0.636 | 0.545 |
| `gpt-5.4-mini` | 311 | 0.267 | 0.778 | 0.887 |
| `gemini-3.1-flash-lite` | 259 | 0.363 | 0.726 | 0.795 |
| `claude-sonnet-5` | 20 | 0.200 | 0.200 | 0.500 |

## Q4. Consistency within a template

The share of templates fully solved on all 15 instances, on some, and on none: the 58 single-path templates (one reasoning path across their instances, `diversity.py`'s lower reading) and the 92 others. Intervals are in `results.json`.

| model | single: all | some | none | others: all | some | none |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.345 | 0.586 | 0.069 | 0.293 | 0.696 | 0.011 |
| `gemma-4-26b-a4b` | 0.466 | 0.500 | 0.034 | 0.522 | 0.478 | 0.000 |
| `deepseek-v4.1-flash` | 0.931 | 0.069 | 0.000 | 0.880 | 0.120 | 0.000 |
| `qwen3-235b-a22b-2507` | 0.448 | 0.534 | 0.017 | 0.554 | 0.424 | 0.022 |
| `glm-5.3-flash` | 0.897 | 0.103 | 0.000 | 0.826 | 0.174 | 0.000 |
| `glm-5.3` | 0.845 | 0.138 | 0.017 | 0.793 | 0.196 | 0.011 |
| `muse-glimmer-30b` | 0.776 | 0.224 | 0.000 | 0.826 | 0.174 | 0.000 |
| `kimi-k3` | 0.862 | 0.138 | 0.000 | 0.815 | 0.185 | 0.000 |
| `gpt-5.4-mini` | 0.500 | 0.431 | 0.069 | 0.478 | 0.500 | 0.022 |
| `gemini-3.1-flash-lite` | 0.672 | 0.276 | 0.052 | 0.598 | 0.402 | 0.000 |
| `claude-sonnet-5` | 0.914 | 0.086 | 0.000 | 0.837 | 0.163 | 0.000 |

## Q5. Paraphrase robustness

The experts' check: 314 pairs returned, 39 rejected and dropped from both arms, 0 not yet returned and kept provisionally.

Paraphrase minus original, paired by item; besides the answer score, E3 coverage on the items with milestones, E5-strict where both arms carry it, and the answer score on the pairs both arms served from the same endpoint (D-149).

| model | items | answer score diff | 95% CI | p (Holm) | McNemar p (Holm) | E3 coverage diff | 95% CI | E5 diff | 95% CI | same provider: pairs, diff |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 275 | +0.018 | -0.031 to 0.067 | 1.0000 | 1.0000 | -0.013 | -0.047 to 0.018 | -0.021 | -0.054 to 0.012 | 273, +0.018 |
| `gemma-4-26b-a4b` | 275 | +0.004 | -0.039 to 0.044 | 1.0000 | 1.0000 | +0.020 | -0.001 to 0.044 | +0.019 | 0.002 to 0.038 | 166, +0.012 |
| `deepseek-v4.1-flash` | 275 | -0.013 | -0.029 to 0.000 | 1.0000 | 1.0000 | -0.000 | -0.019 to 0.018 | -0.009 | -0.026 to 0.008 | 33, -0.030 |
| `qwen3-235b-a22b-2507` | 275 | -0.036 | -0.080 to 0.005 | 1.0000 | 1.0000 | -0.005 | -0.032 to 0.019 | -0.015 | -0.041 to 0.009 | 184, -0.033 |
| `glm-5.3-flash` | 275 | -0.011 | -0.029 to 0.005 | 1.0000 | 1.0000 | -0.028 | -0.049 to -0.008 | -0.015 | -0.034 to 0.003 | 25, +0.000 |
| `glm-5.3` | 275 | +0.002 | -0.018 to 0.022 | 1.0000 | 1.0000 | +0.007 | -0.013 to 0.026 | +0.001 | -0.019 to 0.019 | 88, +0.011 |
| `muse-glimmer-30b` | 275 | +0.025 | 0.004 to 0.051 | 0.5273 | 0.7197 | -0.009 | -0.038 to 0.017 | -0.007 | -0.035 to 0.017 | 275, +0.025 |
| `kimi-k3` | 275 | +0.002 | -0.005 to 0.011 | 1.0000 | 1.0000 | -0.007 | -0.028 to 0.013 | -0.005 | -0.025 to 0.015 | 6, +0.167 |
| `gpt-5.4-mini` | 275 | +0.000 | -0.035 to 0.035 | 1.0000 | 1.0000 | -0.028 | -0.055 to -0.003 | -0.013 | -0.041 to 0.013 | 275, +0.000 |
| `gemini-3.1-flash-lite` | 275 | -0.009 | -0.051 to 0.031 | 1.0000 | 1.0000 | +0.003 | -0.016 to 0.023 | +0.002 | -0.015 to 0.019 | 275, -0.009 |
| `claude-sonnet-5` | 275 | -0.011 | -0.025 to 0.000 | 1.0000 | 1.0000 | -0.007 | -0.024 to 0.008 | -0.005 | -0.023 to 0.013 | 275, -0.011 |

Kendall's tau between the models' answer scores on the originals and on the paraphrases, over the 275 items every tested model holds: 0.807, 95% CI 0.547 to 0.844; the noise floor for tau on this roster is 0.964 (below).

The noise floor at the arm's own size and template mix (D-165): two disjoint draws of 275 main-run items with the arm's per-template counts over its 114 templates, the models ordered on each, 200 draws: median tau 0.807, quartiles 0.754 to 0.860, 5th to 95th percentile 0.660 to 0.925. This, not the roster-wide floor (halves of 7 or 8 items per template), is the comparison for the arm's tau.

**Bounds (D-165).** The same bootstrap's 90% interval per model, read against a margin of ±0.05 in answer score: a model is within the margin when the whole interval is (the two-one-sided-tests rule at 5%). The margin is the plan's largest detectable paired difference (5.1 points at 15% discordance); it was fixed after the point estimates were known and before these intervals were computed.

| model | answer score diff | 90% CI | within the margin | detectable at 0.05 (A9) |
|---|---:|---:|---|---:|
| `gpt-oss-20b` | +0.018 | -0.023 to 0.059 | no | 0.079 |
| `gemma-4-26b-a4b` | +0.004 | -0.031 to 0.038 | yes | 0.064 |
| `deepseek-v4.1-flash` | -0.013 | -0.026 to -0.002 | yes | 0.017 |
| `qwen3-235b-a22b-2507` | -0.036 | -0.073 to -0.002 | no | 0.052 |
| `glm-5.3-flash` | -0.011 | -0.026 to 0.004 | yes | 0.020 |
| `glm-5.3` | +0.002 | -0.015 to 0.019 | yes | 0.034 |
| `muse-glimmer-30b` | +0.025 | 0.007 to 0.046 | yes | 0.044 |
| `kimi-k3` | +0.002 | -0.004 to 0.009 | yes | 0.009 |
| `gpt-5.4-mini` | +0.000 | -0.029 to 0.029 | yes | 0.054 |
| `gemini-3.1-flash-lite` | -0.009 | -0.043 to 0.025 | yes | 0.067 |
| `claude-sonnet-5` | -0.011 | -0.022 to -0.003 | yes | 0.014 |

**Against run-to-run noise (D-165).** For the models with decoding repeats: the paraphrase difference beside each repeat minus the main run on the 300 repeat items, each paired by item with the same template bootstrap; and both on the kept pairs the two subsamples share.

| model | paraphrase − original (pairs) | repeat1 − main (300) | repeat2 − main (300) | repeat3 − main (300) | shared items | paraphrase on them | repeats on them | |paraphrase| within the repeats' spread |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `gpt-oss-20b` | +0.018 -0.031 to 0.067 (275) | +0.008 -0.027 to 0.042 | +0.012 -0.023 to 0.047 | +0.007 -0.030 to 0.042 | 92 | +0.000 -0.076 to 0.076 | -0.033, -0.016, -0.011 | no |
| `gemma-4-26b-a4b` | +0.004 -0.039 to 0.044 (275) | -0.023 -0.058 to 0.012 | -0.002 -0.035 to 0.030 | +0.008 -0.017 to 0.033 | 92 | -0.043 -0.109 to 0.022 | -0.038, +0.000, -0.027 | yes |
| `qwen3-235b-a22b-2507` | -0.036 -0.080 to 0.005 (275) | +0.003 -0.025 to 0.033 | +0.000 -0.028 to 0.028 | +0.010 -0.020 to 0.040 | 92 | -0.022 -0.071 to 0.027 | +0.000, -0.005, +0.033 | no |
| `gemini-3.1-flash-lite` | -0.009 -0.051 to 0.031 (275) | +0.002 -0.022 to 0.025 | -0.008 -0.032 to 0.015 | +0.010 -0.012 to 0.033 | 92 | -0.022 -0.076 to 0.027 | -0.005, -0.022, -0.022 | yes |

## C1 and C4. Arms run on the same items against a base run: reasoning on, open book, and the tool

*Exploratory; added 2026-10-03 (D-179, D-180, D-182, D-183, D-184).* Each arm keeps the base run's prompt, ceiling, routing and scoring on the originals of the 450-item subsample (three per template) and changes one thing. `reasoning-<effort>`: OpenRouter's reasoning parameter at that effort, for the closed roster models whose endpoints reported no reasoning tokens at the provider's default (base: the main run). `openbook`: the template's governing equations appended to the question (`openbook.py`; 429 items, the seven templates without a symbolic equation left out), for one model from each tier (base: the main run). `tool`: a Python tool offered in the request, the model's scripts run in an isolated interpreter and their output returned (`run_traces.py`, D-184), for the same three (base: the main run; the tool-use columns say how much it was used). `flagship-reasoning-<effort>`: the closed anchor of C3 with the parameter, against the `flagship` arm, the same anchor at the provider's default. Arm minus base, paired by item, with Q5's machinery: the item mean, a template bootstrap, a sign-flip test over templates with Holm across the models in the arm, the smallest paired difference the arm detects at 80% power, and McNemar on the fully-solved verdict; the level means are descriptive. The base run's scores remain the headline; an arm says what the change cost or bought.

| arm | model | items | base on these items | arm | change | 95% CI | p (Holm, arm family) | detectable | 90% CI, within ±0.05 | change on items usable in both | run-to-run noise (repeats minus main) | fully solved change | sign-flip p (Holm) | McNemar p (Holm) | by level, change with 95% CI: Easy / Intermediate / Advanced | E3 coverage change (items) | 95% CI | E5-strict change (items) | 95% CI | unusable, base / arm | digit flags on fully solved, base / arm (unpaired) | router judge flags on fully solved, base / arm (unpaired) | tool use: share of traces, calls per trace |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---|---|---|---|
| reasoning-medium | `gpt-5.4-mini` | 450 | 0.862 | 0.964 | +0.102 | 0.063 to 0.143 | <0.0001 | 0.058 | 0.070 to 0.137, no | +0.102 (450) | no repeats | +0.107 | <0.0001 (<0.0001) | <0.0001 (<0.0001) | +0.066 (0.026 to 0.118; 58 templates) / +0.089 (0.029 to 0.155; 58 templates) / +0.186 (0.083 to 0.304; 34 templates) | +0.055 (440) | 0.028 to 0.083 | +0.057 (440) | 0.031 to 0.085 | 0 / 0 | 0.125 / 0.067 | 0.229 / 0.130 |  |
| reasoning-medium | `gemini-3.1-flash-lite` | 450 | 0.884 | 0.898 | +0.013 | -0.014 to 0.042 | 0.4072 | 0.040 | -0.010 to 0.038, yes | +0.013 (450) | -0.008 to +0.010 (3 repeats) | +0.013 | 0.4599 (0.4599) | 0.3915 (0.3915) | +0.006 (-0.029 to 0.034; 58 templates) / +0.020 (-0.029 to 0.072; 58 templates) / +0.015 (-0.054 to 0.088; 34 templates) | +0.013 (440) | -0.007 to 0.033 | +0.017 (440) | -0.001 to 0.036 | 0 / 0 | 0.086 / 0.075 | 0.114 / 0.067 |  |
| openbook | `gpt-oss-20b` | 423 | 0.835 | 0.878 | +0.044 | 0.004 to 0.084 | 0.1114 | 0.057 | 0.009 to 0.077, no | +0.022 (408) | +0.007 to +0.012 (3 repeats) | +0.052 | 0.0184 (0.0553) | 0.0092 (0.0276) | +0.021 (-0.024 to 0.076; 55 templates) / +0.100 (0.036 to 0.173; 55 templates) / -0.016 (-0.113 to 0.070; 31 templates) | +0.064 (414) | 0.032 to 0.097 | +0.064 (414) | 0.034 to 0.096 | 12 / 3 | 0.156 / 0.182 | not run |  |
| openbook | `gpt-5.4-mini` | 423 | 0.877 | 0.870 | -0.007 | -0.034 to 0.020 | 0.8900 | 0.039 | -0.030 to 0.015, yes | -0.007 (423) | no repeats | -0.007 | 0.7638 (1.0000) | 0.7660 (1.0000) | +0.000 (-0.033 to 0.036; 55 templates) / -0.024 (-0.076 to 0.027; 55 templates) / +0.011 (-0.048 to 0.070; 31 templates) | +0.038 (414) | 0.012 to 0.065 | +0.032 (414) | 0.005 to 0.060 | 0 / 0 | 0.114 / 0.157 | not run |  |
| openbook | `claude-sonnet-5` | 423 | 0.981 | 0.988 | +0.007 | -0.006 to 0.021 | 0.8900 | 0.020 | -0.005 to 0.019, yes | +0.007 (423) | no repeats | +0.007 | 0.5610 (1.0000) | 0.5488 (1.0000) | +0.000 (-0.009 to 0.009; 55 templates) / +0.006 (-0.018 to 0.030; 55 templates) / +0.022 (0.000 to 0.065; 31 templates) | +0.033 (414) | 0.015 to 0.053 | +0.031 (414) | 0.013 to 0.051 | 0 / 0 | 0.130 / 0.122 | not run |  |
| openbook2 | `gpt-oss-20b` | 405 | 0.810 | 0.878 | +0.068 | 0.031 to 0.106 | 0.0013 | 0.054 | 0.037 to 0.100, no | +0.044 (388) | +0.007 to +0.012 (3 repeats) | +0.074 | 0.0003 (0.0010) | 0.0002 (0.0005) | +0.028 (-0.022 to 0.085; 53 templates) / +0.103 (0.037 to 0.177; 50 templates) / +0.078 (0.016 to 0.146; 32 templates) | +0.072 (399) | 0.042 to 0.106 | +0.068 (399) | 0.037 to 0.102 | 13 / 5 | 0.165 / 0.179 | not run |  |
| openbook2 | `gpt-5.4-mini` | 405 | 0.853 | 0.837 | -0.016 | -0.047 to 0.016 | 0.6292 | 0.044 | -0.042 to 0.010, yes | -0.016 (405) | no repeats | -0.017 | 0.3773 (0.7547) | 0.3489 (0.6978) | +0.006 (-0.025 to 0.044; 53 templates) / -0.010 (-0.070 to 0.053; 50 templates) / -0.062 (-0.135 to 0.000; 32 templates) | +0.037 (399) | 0.010 to 0.064 | +0.032 (399) | 0.007 to 0.059 | 0 / 0 | 0.123 / 0.161 | not run |  |
| openbook2 | `claude-sonnet-5` | 405 | 0.980 | 0.989 | +0.009 | -0.002 to 0.021 | 0.6292 | 0.017 | 0.000 to 0.020, yes | +0.009 (405) | no repeats | +0.005 | 0.7524 (0.7547) | 0.7266 (0.7266) | +0.000 (0.000 to 0.000; 53 templates) / +0.013 (-0.013 to 0.040; 50 templates) / +0.016 (0.000 to 0.047; 32 templates) | +0.028 (399) | 0.009 to 0.048 | +0.028 (399) | 0.010 to 0.048 | 0 / 0 | 0.146 / 0.133 | not run |  |
| tool | `gpt-5.4-mini` | 450 | 0.862 | 0.859 | -0.003 | -0.040 to 0.036 | 0.9088 | 0.054 | -0.034 to 0.029, yes | -0.003 (450) | no repeats | -0.002 | 1.0000 (1.0000) | 1.0000 (1.0000) | -0.017 (-0.057 to 0.017; 58 templates) / -0.003 (-0.063 to 0.063; 58 templates) / +0.020 (-0.078 to 0.132; 34 templates) | +0.005 (440) | -0.017 to 0.028 | +0.010 (440) | -0.012 to 0.032 | 0 / 0 | 0.125 / 0.094 | not run | 0.18, 0.4 |
| tool | `claude-sonnet-5` | 450 | 0.982 | 0.992 | +0.010 | 0.000 to 0.023 | 0.2521 | 0.017 | 0.001 to 0.021, yes | +0.010 (450) | no repeats | +0.009 | 0.3133 (0.6267) | 0.2891 (0.5781) | +0.000 (0.000 to 0.000; 58 templates) / +0.014 (-0.003 to 0.034; 58 templates) / +0.020 (0.000 to 0.059; 34 templates) | +0.002 (440) | -0.011 to 0.015 | +0.002 (440) | -0.011 to 0.015 | 0 / 0 | 0.136 / 0.128 | not run | 0.66, 0.9 |
| flagship-reasoning-medium vs flagship | `gpt-5.4` | 450 | 0.951 | 0.997 | +0.046 | 0.021 to 0.074 | 0.0002 | 0.039 | 0.024 to 0.069, no | +0.046 (450) | no repeats | +0.047 | 0.0003 (0.0003) | <0.0001 (<0.0001) | +0.009 (0.000 to 0.023; 58 templates) / +0.046 (0.011 to 0.092; 58 templates) / +0.108 (0.029 to 0.206; 34 templates) | +0.015 (440) | 0.003 to 0.027 | +0.012 (440) | -0.002 to 0.025 | 0 / 0 | 0.117 / 0.027 | not run |  |

Reading the columns: a claim rests on the template-level tests (the sign-flip p with Holm over the arm's models); McNemar is item-level and shown as Q5 shows it. The 90% interval is the two-one-sided-tests reading against the ±0.05 equivalence margin of D-165: "yes" means the change is bounded inside the margin. "Change on items usable in both" leaves out items unusable in either arm, so a change that comes from fewer empty traces shows as the gap between the two. The run-to-run noise is the model's decoding repeats minus the main run where repeats exist (D-151). The flag rates are on each arm's own fully solved traces, unpaired and without intervals: descriptive.

Tool arm, the answer score by whether the trace called the tool (descriptive: whether to call it is the model's choice and depends on the item): `gpt-5.4-mini` with calls 0.872, without 0.856, 0 traces at the call limit; `claude-sonnet-5` with calls 0.988, without 1.000, 0 traces at the call limit.

## C3. Flagship anchors on the 450-item subsample: `flagship`

*Exploratory; added 2026-10-03 (D-182).* Flagships that pass the roster rule (neither a pilot generator nor a judge), on the originals of the 450-item subsample, scored by the same stack; the roster's eleven on the same items from the main run beside them. An anchor is a reference point outside the pairwise family: no test is run against it. The `flagship` arm ran at each provider's default, as the main run did; a `flagship-reasoning-<effort>` arm carries the reasoning parameter. Score and coverage carry template intervals (150 templates of three items), as do the level means. GPT-5.4's reasoning arm is paired against its default arm in the block above, the one paired test among the anchors.

| model | arm | items | answer score | 95% CI | fully solved | unusable (empty / unreadable) | Easy, 95% CI | Intermediate, 95% CI | Advanced, 95% CI | E3 coverage | 95% CI | coverage (stage) | 95% CI | digit flags on fully solved | median completion tokens |
|---|---|---:|---:|---:|---:|---:|---|---|---|---:|---:|---|---:|---:|---:|
| `deepseek-v4-pro` | flagship | 450 | 0.968 | 0.947 to 0.984 | 0.962 | 4 (4 / 0) | 0.974 (0.951 to 0.991) | 0.971 (0.945 to 0.991) | 0.951 (0.873 to 1.000) | 0.873 | 0.839 to 0.904 | 0.907 (E5-strict) | 0.877 to 0.933 | 0.104 | 2596 |
| `gpt-5.4` | flagship | 450 | 0.951 | 0.920 to 0.977 | 0.949 | 0 (0 / 0) | 0.989 (0.974 to 1.000) | 0.948 (0.897 to 0.989) | 0.892 (0.794 to 0.971) | 0.862 | 0.826 to 0.894 | 0.894 (E5-strict) | 0.866 to 0.921 | 0.117 | 738 |
| `deepseek-v4.1-flash` | main | 450 | 0.990 | 0.981 to 0.998 | 0.987 | 1 (1 / 0) | 0.997 (0.991 to 1.000) | 0.986 (0.968 to 1.000) | 0.985 (0.961 to 1.000) | 0.860 | 0.828 to 0.891 | 0.895 (E5-strict) | 0.868 to 0.920 | 0.011 | 3178 |
| `claude-sonnet-5` | main | 450 | 0.982 | 0.968 to 0.993 | 0.980 | 0 (0 / 0) | 0.994 (0.986 to 1.000) | 0.971 (0.943 to 0.994) | 0.980 (0.941 to 1.000) | 0.907 | 0.883 to 0.929 | 0.924 (E5-strict) | 0.902 to 0.944 | 0.136 | 1773 |
| `kimi-k3` | main | 450 | 0.979 | 0.959 to 0.994 | 0.976 | 5 (4 / 1) | 0.991 (0.980 to 1.000) | 0.989 (0.971 to 1.000) | 0.941 (0.863 to 1.000) | 0.875 | 0.844 to 0.903 | 0.912 (E5-strict) | 0.885 to 0.937 | 0.007 | 2294 |
| `glm-5.3-flash` | main | 450 | 0.977 | 0.958 to 0.991 | 0.973 | 7 (7 / 0) | 0.989 (0.974 to 1.000) | 0.971 (0.931 to 1.000) | 0.966 (0.931 to 0.995) | 0.889 | 0.858 to 0.918 | 0.909 (E5-strict) | 0.880 to 0.935 | 0.011 | 2600 |
| `muse-glimmer-30b` | main | 450 | 0.959 | 0.934 to 0.981 | 0.956 | 6 (6 / 0) | 0.971 (0.937 to 0.997) | 0.957 (0.908 to 0.994) | 0.941 (0.882 to 0.990) | 0.861 | 0.829 to 0.892 | 0.898 (E5-strict) | 0.870 to 0.924 | 0.072 | 2972 |
| `glm-5.3` | main | 450 | 0.949 | 0.921 to 0.973 | 0.947 | 17 (17 / 0) | 0.989 (0.974 to 1.000) | 0.948 (0.902 to 0.983) | 0.882 (0.794 to 0.961) | 0.879 | 0.844 to 0.911 | 0.900 (E5-strict) | 0.867 to 0.929 | 0.023 | 3155 |
| `qwen3-235b-a22b-2507` | main | 450 | 0.898 | 0.862 to 0.931 | 0.873 | 0 (0 / 0) | 0.931 (0.885 to 0.971) | 0.871 (0.796 to 0.934) | 0.887 (0.824 to 0.946) | 0.863 | 0.830 to 0.895 | 0.894 (E5-strict) | 0.866 to 0.921 | 0.188 | 1020 |
| `gemini-3.1-flash-lite` | main | 450 | 0.884 | 0.840 to 0.924 | 0.878 | 0 (0 / 0) | 0.940 (0.891 to 0.977) | 0.888 (0.810 to 0.951) | 0.784 (0.672 to 0.887) | 0.837 | 0.802 to 0.871 | 0.862 (E5-strict) | 0.829 to 0.893 | 0.086 | 659 |
| `gemma-4-26b-a4b` | main | 450 | 0.868 | 0.827 to 0.906 | 0.856 | 0 (0 / 0) | 0.908 (0.845 to 0.960) | 0.853 (0.787 to 0.914) | 0.824 (0.735 to 0.902) | 0.848 | 0.814 to 0.880 | 0.874 (E5-strict) | 0.845 to 0.901 | 0.119 | 806 |
| `gpt-5.4-mini` | main | 450 | 0.862 | 0.817 to 0.904 | 0.853 | 0 (0 / 0) | 0.922 (0.871 to 0.966) | 0.853 (0.770 to 0.925) | 0.775 (0.657 to 0.877) | 0.799 | 0.759 to 0.839 | 0.835 (E5-strict) | 0.796 to 0.872 | 0.125 | 643 |
| `gpt-oss-20b` | main | 450 | 0.821 | 0.771 to 0.868 | 0.804 | 13 (11 / 2) | 0.920 (0.862 to 0.971) | 0.782 (0.695 to 0.859) | 0.721 (0.603 to 0.828) | 0.802 | 0.763 to 0.839 | 0.826 (E5-strict) | 0.788 to 0.861 | 0.157 | 1816 |

Level means carry template intervals (58, 58, 34 templates for Easy, Intermediate and Advanced on the subsample); the anchors' Advanced intervals are wide and overlap the roster's, so a level difference between an anchor and a roster model is not a finding.

## C3. Flagship anchors on the 450-item subsample: `flagship-reasoning-medium`

*Exploratory; added 2026-10-03 (D-182).* Flagships that pass the roster rule (neither a pilot generator nor a judge), on the originals of the 450-item subsample, scored by the same stack; the roster's eleven on the same items from the main run beside them. An anchor is a reference point outside the pairwise family: no test is run against it. The `flagship` arm ran at each provider's default, as the main run did; a `flagship-reasoning-<effort>` arm carries the reasoning parameter. Score and coverage carry template intervals (150 templates of three items), as do the level means. GPT-5.4's reasoning arm is paired against its default arm in the block above, the one paired test among the anchors.

| model | arm | items | answer score | 95% CI | fully solved | unusable (empty / unreadable) | Easy, 95% CI | Intermediate, 95% CI | Advanced, 95% CI | E3 coverage | 95% CI | coverage (stage) | 95% CI | digit flags on fully solved | median completion tokens |
|---|---|---:|---:|---:|---:|---:|---|---|---|---:|---:|---|---:|---:|---:|
| `gpt-5.4` | flagship-reasoning-medium | 450 | 0.997 | 0.991 to 1.000 | 0.996 | 0 (0 / 0) | 0.997 (0.991 to 1.000) | 0.994 (0.983 to 1.000) | 1.000 (1.000 to 1.000) | 0.877 | 0.845 to 0.906 | 0.906 (E5-strict) | 0.877 to 0.931 | 0.027 | 1239 |
| `deepseek-v4.1-flash` | main | 450 | 0.990 | 0.981 to 0.998 | 0.987 | 1 (1 / 0) | 0.997 (0.991 to 1.000) | 0.986 (0.968 to 1.000) | 0.985 (0.961 to 1.000) | 0.860 | 0.828 to 0.891 | 0.895 (E5-strict) | 0.868 to 0.920 | 0.011 | 3178 |
| `claude-sonnet-5` | main | 450 | 0.982 | 0.968 to 0.993 | 0.980 | 0 (0 / 0) | 0.994 (0.986 to 1.000) | 0.971 (0.943 to 0.994) | 0.980 (0.941 to 1.000) | 0.907 | 0.883 to 0.929 | 0.924 (E5-strict) | 0.902 to 0.944 | 0.136 | 1773 |
| `kimi-k3` | main | 450 | 0.979 | 0.959 to 0.994 | 0.976 | 5 (4 / 1) | 0.991 (0.980 to 1.000) | 0.989 (0.971 to 1.000) | 0.941 (0.863 to 1.000) | 0.875 | 0.844 to 0.903 | 0.912 (E5-strict) | 0.885 to 0.937 | 0.007 | 2294 |
| `glm-5.3-flash` | main | 450 | 0.977 | 0.958 to 0.991 | 0.973 | 7 (7 / 0) | 0.989 (0.974 to 1.000) | 0.971 (0.931 to 1.000) | 0.966 (0.931 to 0.995) | 0.889 | 0.858 to 0.918 | 0.909 (E5-strict) | 0.880 to 0.935 | 0.011 | 2600 |
| `muse-glimmer-30b` | main | 450 | 0.959 | 0.934 to 0.981 | 0.956 | 6 (6 / 0) | 0.971 (0.937 to 0.997) | 0.957 (0.908 to 0.994) | 0.941 (0.882 to 0.990) | 0.861 | 0.829 to 0.892 | 0.898 (E5-strict) | 0.870 to 0.924 | 0.072 | 2972 |
| `glm-5.3` | main | 450 | 0.949 | 0.921 to 0.973 | 0.947 | 17 (17 / 0) | 0.989 (0.974 to 1.000) | 0.948 (0.902 to 0.983) | 0.882 (0.794 to 0.961) | 0.879 | 0.844 to 0.911 | 0.900 (E5-strict) | 0.867 to 0.929 | 0.023 | 3155 |
| `qwen3-235b-a22b-2507` | main | 450 | 0.898 | 0.862 to 0.931 | 0.873 | 0 (0 / 0) | 0.931 (0.885 to 0.971) | 0.871 (0.796 to 0.934) | 0.887 (0.824 to 0.946) | 0.863 | 0.830 to 0.895 | 0.894 (E5-strict) | 0.866 to 0.921 | 0.188 | 1020 |
| `gemini-3.1-flash-lite` | main | 450 | 0.884 | 0.840 to 0.924 | 0.878 | 0 (0 / 0) | 0.940 (0.891 to 0.977) | 0.888 (0.810 to 0.951) | 0.784 (0.672 to 0.887) | 0.837 | 0.802 to 0.871 | 0.862 (E5-strict) | 0.829 to 0.893 | 0.086 | 659 |
| `gemma-4-26b-a4b` | main | 450 | 0.868 | 0.827 to 0.906 | 0.856 | 0 (0 / 0) | 0.908 (0.845 to 0.960) | 0.853 (0.787 to 0.914) | 0.824 (0.735 to 0.902) | 0.848 | 0.814 to 0.880 | 0.874 (E5-strict) | 0.845 to 0.901 | 0.119 | 806 |
| `gpt-5.4-mini` | main | 450 | 0.862 | 0.817 to 0.904 | 0.853 | 0 (0 / 0) | 0.922 (0.871 to 0.966) | 0.853 (0.770 to 0.925) | 0.775 (0.657 to 0.877) | 0.799 | 0.759 to 0.839 | 0.835 (E5-strict) | 0.796 to 0.872 | 0.125 | 643 |
| `gpt-oss-20b` | main | 450 | 0.821 | 0.771 to 0.868 | 0.804 | 13 (11 / 2) | 0.920 (0.862 to 0.971) | 0.782 (0.695 to 0.859) | 0.721 (0.603 to 0.828) | 0.802 | 0.763 to 0.839 | 0.826 (E5-strict) | 0.788 to 0.861 | 0.157 | 1816 |

Level means carry template intervals (58, 58, 34 templates for Easy, Intermediate and Advanced on the subsample); the anchors' Advanced intervals are wide and overlap the roster's, so a level difference between an anchor and a roster model is not a finding.

## Sensitivity

Answer score under each variation; the last row is Kendall's tau between that ordering of the models and the headline one. The three added readings (D-147, D-149) the experts could not arbitrate: the half-unit window requires a correct rounding at the precision shown where the rule accepts one unit either way; the whole-trace reading credits a quantity the question asks for when it is stated in the body and left off the Answer line (the prompt asks for it there); the pool without the nine symbolic templates, whose answers the check scores by the numbers they state (D-138). "Unusable excluded" is an item mean over the usable rows, not a mean of template means. Tau's noise floor: the ordering on one random half of each template's items against the other, median 0.964 over 200 splits (quartiles 0.927 to 0.964); the top five models lie within 0.012 of each other, so tau falls below 1 from sampling noise alone.

| model | tolerance half | fitted (headline) | tolerance double | fully solved | unusable excluded | without the 4 shortcut templates | without the 9 symbolic templates | half-unit window | whole trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.801 | 0.808 | 0.815 | 0.794 | 0.830 | 0.806 | 0.802 | 0.780 | 0.811 |
| `gemma-4-26b-a4b` | 0.850 | 0.872 | 0.904 | 0.860 | 0.872 | 0.869 | 0.867 | 0.855 | 0.876 |
| `deepseek-v4.1-flash` | 0.986 | 0.988 | 0.990 | 0.985 | 0.992 | 0.989 | 0.987 | 0.983 | 0.989 |
| `qwen3-235b-a22b-2507` | 0.888 | 0.899 | 0.910 | 0.876 | 0.899 | 0.897 | 0.898 | 0.889 | 0.917 |
| `glm-5.3-flash` | 0.971 | 0.972 | 0.973 | 0.970 | 0.992 | 0.971 | 0.972 | 0.968 | 0.974 |
| `glm-5.3` | 0.949 | 0.950 | 0.951 | 0.948 | 0.994 | 0.949 | 0.949 | 0.944 | 0.951 |
| `muse-glimmer-30b` | 0.967 | 0.969 | 0.970 | 0.967 | 0.983 | 0.969 | 0.968 | 0.952 | 0.970 |
| `kimi-k3` | 0.979 | 0.980 | 0.982 | 0.977 | 0.987 | 0.981 | 0.979 | 0.975 | 0.982 |
| `gpt-5.4-mini` | 0.840 | 0.852 | 0.864 | 0.842 | 0.852 | 0.848 | 0.846 | 0.819 | 0.854 |
| `gemini-3.1-flash-lite` | 0.869 | 0.878 | 0.892 | 0.872 | 0.878 | 0.876 | 0.877 | 0.870 | 0.880 |
| `claude-sonnet-5` | 0.983 | 0.987 | 0.989 | 0.984 | 0.988 | 0.987 | 0.988 | 0.980 | 0.988 |
| tau with the headline | 1.000 | 1 | 0.964 | 1.000 | 0.709 | 1.000 | 0.964 | 1.000 | 1.000 |

The plan's fourth sensitivity, without the two templates widened for round 4, applies only if round 4 had not returned; it returned and certified both (template_annotation_23092026/layer2/CERTIFICATION.md).

## Also reported, not tested

### By branch and level

| model | chemical | civil | electrical | industrial | mechanical | Easy | Intermediate | Advanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.798 | 0.671 | 0.897 | 0.816 | 0.859 | 0.900 | 0.768 | 0.719 |
| `gemma-4-26b-a4b` | 0.754 | 0.882 | 0.927 | 0.878 | 0.919 | 0.919 | 0.866 | 0.802 |
| `deepseek-v4.1-flash` | 0.980 | 0.993 | 0.990 | 0.980 | 0.998 | 0.995 | 0.986 | 0.979 |
| `qwen3-235b-a22b-2507` | 0.869 | 0.940 | 0.911 | 0.871 | 0.904 | 0.917 | 0.890 | 0.885 |
| `glm-5.3-flash` | 0.960 | 0.982 | 0.981 | 0.942 | 0.996 | 0.994 | 0.970 | 0.939 |
| `glm-5.3` | 0.900 | 0.973 | 0.961 | 0.918 | 0.998 | 0.984 | 0.966 | 0.865 |
| `muse-glimmer-30b` | 0.963 | 0.962 | 0.989 | 0.938 | 0.993 | 0.976 | 0.968 | 0.958 |
| `kimi-k3` | 0.960 | 0.982 | 0.984 | 0.980 | 0.996 | 0.992 | 0.985 | 0.953 |
| `gpt-5.4-mini` | 0.807 | 0.838 | 0.926 | 0.831 | 0.858 | 0.916 | 0.844 | 0.756 |
| `gemini-3.1-flash-lite` | 0.784 | 0.862 | 0.922 | 0.860 | 0.962 | 0.937 | 0.888 | 0.762 |
| `claude-sonnet-5` | 0.993 | 0.996 | 0.974 | 0.973 | 0.999 | 0.993 | 0.983 | 0.984 |

### By branch and level, with template intervals

*Added 2026-10-03 (D-173, next steps A2).* The mean of the branch's template means with its template bootstrap (30 templates per branch); the last column is the smallest difference between two branches a model's within-branch spread lets 30 templates detect at 80% power (2.80 x pooled SD x sqrt(2/30)). Then the same by level.

| model | chemical | civil | electrical | industrial | mechanical | detectable branch difference |
|---|---:|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.798 (0.664 to 0.907) | 0.671 (0.540 to 0.796) | 0.897 (0.852 to 0.937) | 0.816 (0.738 to 0.882) | 0.859 (0.791 to 0.921) | 0.188 |
| `gemma-4-26b-a4b` | 0.754 (0.640 to 0.857) | 0.882 (0.816 to 0.940) | 0.927 (0.887 to 0.961) | 0.878 (0.773 to 0.962) | 0.919 (0.874 to 0.958) | 0.152 |
| `deepseek-v4.1-flash` | 0.980 (0.953 to 0.996) | 0.993 (0.982 to 1.000) | 0.990 (0.980 to 0.998) | 0.980 (0.944 to 1.000) | 0.998 (0.993 to 1.000) | 0.037 |
| `qwen3-235b-a22b-2507` | 0.869 (0.782 to 0.940) | 0.940 (0.911 to 0.964) | 0.911 (0.860 to 0.954) | 0.871 (0.773 to 0.951) | 0.904 (0.851 to 0.951) | 0.129 |
| `glm-5.3-flash` | 0.960 (0.902 to 0.998) | 0.982 (0.949 to 1.000) | 0.981 (0.963 to 0.996) | 0.942 (0.869 to 0.989) | 0.996 (0.990 to 1.000) | 0.080 |
| `glm-5.3` | 0.900 (0.800 to 0.987) | 0.973 (0.927 to 1.000) | 0.961 (0.931 to 0.987) | 0.918 (0.831 to 0.984) | 0.998 (0.993 to 1.000) | 0.123 |
| `muse-glimmer-30b` | 0.963 (0.909 to 0.996) | 0.962 (0.924 to 0.991) | 0.989 (0.980 to 0.997) | 0.938 (0.864 to 0.989) | 0.993 (0.986 to 0.999) | 0.081 |
| `kimi-k3` | 0.960 (0.900 to 0.996) | 0.982 (0.951 to 1.000) | 0.984 (0.973 to 0.993) | 0.980 (0.962 to 0.993) | 0.996 (0.990 to 1.000) | 0.055 |
| `gpt-5.4-mini` | 0.807 (0.684 to 0.913) | 0.838 (0.758 to 0.907) | 0.926 (0.878 to 0.967) | 0.831 (0.724 to 0.922) | 0.858 (0.756 to 0.942) | 0.184 |
| `gemini-3.1-flash-lite` | 0.784 (0.658 to 0.898) | 0.862 (0.764 to 0.942) | 0.922 (0.876 to 0.962) | 0.860 (0.751 to 0.949) | 0.962 (0.927 to 0.989) | 0.173 |
| `claude-sonnet-5` | 0.993 (0.987 to 1.000) | 0.996 (0.987 to 1.000) | 0.974 (0.952 to 0.991) | 0.973 (0.949 to 0.993) | 0.999 (0.997 to 1.000) | 0.030 |

| model | Easy (58) | Intermediate (58) | Advanced (34) |
|---|---:|---:|---:|
| `gpt-oss-20b` | 0.900 (0.843 to 0.948) | 0.768 (0.694 to 0.837) | 0.719 (0.614 to 0.812) |
| `gemma-4-26b-a4b` | 0.919 (0.870 to 0.959) | 0.866 (0.802 to 0.923) | 0.802 (0.724 to 0.873) |
| `deepseek-v4.1-flash` | 0.995 (0.990 to 0.999) | 0.986 (0.967 to 0.998) | 0.979 (0.956 to 0.996) |
| `qwen3-235b-a22b-2507` | 0.917 (0.870 to 0.955) | 0.890 (0.834 to 0.937) | 0.885 (0.831 to 0.930) |
| `glm-5.3-flash` | 0.994 (0.989 to 0.999) | 0.970 (0.930 to 0.993) | 0.939 (0.883 to 0.982) |
| `glm-5.3` | 0.984 (0.968 to 0.996) | 0.966 (0.925 to 0.990) | 0.865 (0.761 to 0.955) |
| `muse-glimmer-30b` | 0.976 (0.945 to 0.997) | 0.968 (0.933 to 0.990) | 0.958 (0.915 to 0.989) |
| `kimi-k3` | 0.992 (0.986 to 0.997) | 0.985 (0.975 to 0.994) | 0.953 (0.898 to 0.992) |
| `gpt-5.4-mini` | 0.916 (0.868 to 0.955) | 0.844 (0.771 to 0.905) | 0.756 (0.642 to 0.858) |
| `gemini-3.1-flash-lite` | 0.937 (0.890 to 0.973) | 0.888 (0.821 to 0.945) | 0.762 (0.652 to 0.861) |
| `claude-sonnet-5` | 0.993 (0.982 to 0.999) | 0.983 (0.970 to 0.994) | 0.984 (0.971 to 0.994) |

**Branch pairs within a model.** Welch's t-test on the two branches' template means, Holm over the pairs within a model; the pairs that hold at 0.05 are listed as the higher branch, the lower, and the difference. A sentence of the form "branch X is hardest" needs that branch below every other at this test; "X is harder than Y" needs the pair listed.

| model | pairs that hold | which |
|---|---:|---|
| `gpt-oss-20b` | 1 of 10 | electrical > civil (0.226) |
| `gemma-4-26b-a4b` | 0 of 10 | none |
| `deepseek-v4.1-flash` | 0 of 10 | none |
| `qwen3-235b-a22b-2507` | 0 of 10 | none |
| `glm-5.3-flash` | 0 of 10 | none |
| `glm-5.3` | 0 of 10 | none |
| `muse-glimmer-30b` | 0 of 10 | none |
| `kimi-k3` | 0 of 10 | none |
| `gpt-5.4-mini` | 0 of 10 | none |
| `gemini-3.1-flash-lite` | 0 of 10 | none |
| `claude-sonnet-5` | 0 of 10 | none |

### By domain

| domain | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| digital_communications | 0.933 | 0.917 | 0.990 | 0.897 | 0.980 | 0.947 | 1.000 | 0.987 | 0.967 | 0.917 | 0.980 |
| electromagnetics_and_waves | 0.867 | 0.887 | 0.990 | 0.890 | 0.977 | 0.943 | 0.987 | 0.977 | 0.830 | 0.877 | 0.953 |
| fluid_mechanics | 0.887 | 0.893 | 1.000 | 0.913 | 0.997 | 1.000 | 0.997 | 1.000 | 0.963 | 0.953 | 1.000 |
| geotechnical_engineering | 0.513 | 0.913 | 1.000 | 0.947 | 0.993 | 0.993 | 0.927 | 1.000 | 0.913 | 0.960 | 0.987 |
| mechanics_of_materials | 0.833 | 0.907 | 1.000 | 0.943 | 0.993 | 1.000 | 0.997 | 0.993 | 0.803 | 0.940 | 1.000 |
| production_and_inventory | 0.767 | 0.873 | 0.940 | 0.847 | 0.887 | 0.880 | 0.880 | 0.967 | 0.780 | 0.853 | 0.987 |
| quality_and_reliability_control | 0.773 | 0.800 | 1.000 | 0.787 | 0.947 | 0.907 | 0.933 | 0.980 | 0.760 | 0.800 | 0.940 |
| reaction_kinetics | 0.953 | 0.913 | 0.993 | 0.980 | 1.000 | 0.993 | 1.000 | 0.993 | 0.913 | 0.913 | 1.000 |
| signals_and_systems | 0.890 | 0.977 | 0.990 | 0.947 | 0.987 | 0.993 | 0.980 | 0.990 | 0.980 | 0.973 | 0.990 |
| stochastic_operations | 0.907 | 0.960 | 1.000 | 0.980 | 0.993 | 0.967 | 1.000 | 0.993 | 0.953 | 0.927 | 0.993 |
| structural_analysis | 0.733 | 0.913 | 1.000 | 0.947 | 1.000 | 1.000 | 0.987 | 0.993 | 0.827 | 0.913 | 1.000 |
| thermodynamics | 0.650 | 0.567 | 0.961 | 0.833 | 0.906 | 0.767 | 0.922 | 0.917 | 0.683 | 0.700 | 0.989 |
| transport_phenomena | 0.825 | 0.838 | 0.992 | 0.783 | 0.992 | 0.983 | 0.979 | 0.983 | 0.858 | 0.750 | 0.992 |
| vibrations_and_acoustics | 0.857 | 0.957 | 0.993 | 0.857 | 0.997 | 0.993 | 0.987 | 0.993 | 0.807 | 0.993 | 0.997 |
| water_resources | 0.767 | 0.820 | 0.980 | 0.927 | 0.953 | 0.927 | 0.973 | 0.953 | 0.773 | 0.713 | 1.000 |

### By answer type

| answer type | `gpt-oss-20b` | `gemma-4-26b-a4b` | `deepseek-v4.1-flash` | `qwen3-235b-a22b-2507` | `glm-5.3-flash` | `glm-5.3` | `muse-glimmer-30b` | `kimi-k3` | `gpt-5.4-mini` | `gemini-3.1-flash-lite` | `claude-sonnet-5` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| array | 0.886 | 0.867 | 0.962 | 0.924 | 0.905 | 0.762 | 0.971 | 0.943 | 0.790 | 0.771 | 0.990 |
| classification | 0.680 | 0.953 | 0.960 | 0.887 | 0.993 | 1.000 | 0.960 | 0.960 | 0.973 | 0.973 | 1.000 |
| multipart | 0.877 | 0.901 | 0.993 | 0.903 | 0.985 | 0.956 | 0.980 | 0.988 | 0.905 | 0.921 | 0.989 |
| scalar | 0.766 | 0.841 | 0.988 | 0.900 | 0.969 | 0.955 | 0.962 | 0.980 | 0.819 | 0.860 | 0.987 |
| symbolic | 0.900 | 0.948 | 1.000 | 0.911 | 0.981 | 0.963 | 0.993 | 0.996 | 0.944 | 0.896 | 0.974 |
| vector | 0.904 | 0.967 | 1.000 | 0.850 | 0.988 | 0.988 | 0.983 | 0.988 | 0.871 | 0.925 | 0.983 |

### Tokens against score

Median completion tokens as billed.

| model | answer score | median tokens | on fully solved | on the rest |
|---|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.808 | 1795 | 1583 | 3138 |
| `gemma-4-26b-a4b` | 0.872 | 797 | 772 | 1016 |
| `deepseek-v4.1-flash` | 0.988 | 3060 | 3025 | 16700 |
| `qwen3-235b-a22b-2507` | 0.899 | 1013 | 998 | 1218 |
| `glm-5.3-flash` | 0.972 | 2500 | 2406 | 32768 |
| `glm-5.3` | 0.950 | 3270 | 3073 | 32768 |
| `muse-glimmer-30b` | 0.969 | 2988 | 2952 | 7353 |
| `kimi-k3` | 0.980 | 2192 | 2177 | 6598 |
| `gpt-5.4-mini` | 0.852 | 632 | 610 | 836 |
| `gemini-3.1-flash-lite` | 0.878 | 653 | 635 | 819 |
| `claude-sonnet-5` | 0.987 | 1804 | 1796 | 2860 |

### By serving endpoint

Open-weight models were served by several endpoints under the fp8-or-better rule (D-133). The harness dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is confounded with the templates each endpoint happened to serve; the matched difference compares an endpoint with the other endpoints on the same templates (D-149).

| model | endpoint | rows | raw score | unusable | matched difference | templates matched |
|---|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | Darkbloom | 2240 | 0.808 | 0.026 | -0.228 | 9 |
| `gpt-oss-20b` | DekaLLM | 9 | 0.889 | 0.000 | +0.239 | 8 |
| `gpt-oss-20b` | DeepInfra | 1 | 1.000 | 0.000 | +0.143 | 1 |
| `gemma-4-26b-a4b` | DekaLLM | 1245 | 0.883 | 0.000 | +0.002 | 148 |
| `gemma-4-26b-a4b` | NextBit | 975 | 0.875 | 0.000 | -0.002 | 148 |
| `gemma-4-26b-a4b` | Io Net | 30 | 0.333 | 0.000 |  | 0 |
| `deepseek-v4.1-flash` | CoreWeave | 782 | 0.992 | 0.001 | -0.000 | 141 |
| `deepseek-v4.1-flash` | Parasail | 466 | 0.992 | 0.000 | -0.002 | 138 |
| `deepseek-v4.1-flash` | Morph | 456 | 0.991 | 0.000 | +0.002 | 133 |
| `deepseek-v4.1-flash` | DeepInfra | 411 | 0.984 | 0.005 | +0.000 | 78 |
| `deepseek-v4.1-flash` | Makora | 87 | 1.000 | 0.000 | +0.003 | 43 |
| `deepseek-v4.1-flash` | InferenceNet | 30 | 0.800 | 0.167 |  | 0 |
| `deepseek-v4.1-flash` | Novita | 16 | 1.000 | 0.000 | +0.016 | 7 |
| `deepseek-v4.1-flash` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `qwen3-235b-a22b-2507` | GMICloud | 2192 | 0.899 | 0.000 | -0.010 | 35 |
| `qwen3-235b-a22b-2507` | Parasail | 54 | 0.907 | 0.000 | +0.006 | 35 |
| `qwen3-235b-a22b-2507` | Nebius | 4 | 1.000 | 0.000 | +0.182 | 1 |
| `glm-5.3-flash` | GMICloud | 981 | 0.976 | 0.015 | -0.003 | 145 |
| `glm-5.3-flash` | Novita | 419 | 0.962 | 0.031 | +0.003 | 120 |
| `glm-5.3-flash` | Parasail | 220 | 0.964 | 0.036 | +0.002 | 75 |
| `glm-5.3-flash` | Sail Research | 201 | 0.965 | 0.020 | -0.022 | 114 |
| `glm-5.3-flash` | Phala | 166 | 0.991 | 0.006 | +0.004 | 74 |
| `glm-5.3-flash` | StreamLake | 80 | 0.950 | 0.037 | -0.025 | 35 |
| `glm-5.3-flash` | Z.AI | 79 | 0.994 | 0.000 | +0.004 | 48 |
| `glm-5.3-flash` | NextBit | 43 | 0.965 | 0.023 | +0.031 | 13 |
| `glm-5.3-flash` | AtlasCloud | 42 | 1.000 | 0.000 | +0.014 | 15 |
| `glm-5.3-flash` | SiliconFlow | 18 | 1.000 | 0.000 | +0.000 | 5 |
| `glm-5.3-flash` | Morph | 1 | 1.000 | 0.000 | +0.000 | 1 |
| `glm-5.3` | Baidu | 918 | 0.956 | 0.036 | -0.006 | 86 |
| `glm-5.3` | Sail Research | 739 | 0.964 | 0.030 | -0.004 | 77 |
| `glm-5.3` | Morph | 257 | 0.940 | 0.054 | +0.021 | 123 |
| `glm-5.3` | Novita | 167 | 0.943 | 0.054 | -0.004 | 41 |
| `glm-5.3` | AtlasCloud | 137 | 0.927 | 0.073 | -0.010 | 54 |
| `glm-5.3` | AkashML | 15 | 0.267 | 0.733 | -0.020 | 5 |
| `glm-5.3` | BaseTen | 15 | 1.000 | 0.000 | +0.000 | 5 |
| `glm-5.3` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 |
| `kimi-k3` | Morph | 2231 | 0.980 | 0.007 | -0.017 | 7 |
| `kimi-k3` | BaseTen | 19 | 1.000 | 0.000 | +0.017 | 7 |
| `gpt-5.4-mini` | Azure | 2249 | 0.852 | 0.000 | -0.036 | 1 |
| `gpt-5.4-mini` | OpenAI | 1 | 1.000 | 0.000 | +0.036 | 1 |

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
| `gpt-oss-20b` | 300 | 0.835 | 0.838 | 0.833 | 0.827 | 0.003 | 0.005 | 0.910 |
| `gemma-4-26b-a4b` | 300 | 0.862 | 0.883 | 0.893 | 0.885 | 0.016 | 0.032 | 0.857 |
| `qwen3-235b-a22b-2507` | 300 | 0.910 | 0.907 | 0.917 | 0.907 | 0.005 | 0.010 | 0.893 |
| `gemini-3.1-flash-lite` | 300 | 0.900 | 0.890 | 0.908 | 0.898 | 0.009 | 0.018 | 0.933 |

## Set aside (D-132), outside every comparison

- `qwen3.8-27b`: answer score 0.956 (95% CI 0.932 to 0.976), fully solved 0.953, unusable 61

## Provenance

`analyze.py` at commit `b8b3358`; the score store `main` scored at commit `b8b3358`, dirty on 2026-10-08T15:11:25+00:00; stages: e5 at `b8b3358`, router at `b8b3358`. The evaluator hashes and the per-model trace hashes are in `results.json`.
