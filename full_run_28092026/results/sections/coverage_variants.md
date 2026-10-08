## Coverage variants

Printed by `coverage_variants.py`; the method is in its docstring. MC as scored is Table 1's Milestone Coverage, recomputed from the judge's per-milestone sources and checked equal to the stored `e5_strict` for every response. Each model's value is the mean of its template means over the templates on which every model has a value under that reading, with a 95% template-bootstrap interval; "wrong" is the item mean over readable responses scoring 0.

**Target rule.** A milestone is an answer target when a number in the answer row's targets.numbers equals its value under one of the answer check's unit factors (answer.SCALES) within the milestone display tolerance (0.5%), the test answer.targets uses to make a gold answer value a target; milestones.py records no answer line. Of 2,250 instances, 70 have no milestones and 440 have only target milestones, so 510 are left without milestones under intermediate-only MC (23 templates have none left). 2,736 of 8,754 milestones are targets; on 159 instances the gold's answer value is not a milestone, so all their milestones stay. Templates in each reading's comparison: as scored (E3 + judge) 147, matching alone (E3) 147, route-adjusted 147, intermediate-only 127.

**Route-adjusted.** (e3 + REACHED) / (n - NOT_NEEDED); a response whose milestones the judge rules all not needed has no value and is left out (count per model): gpt-oss-20b 7, gemma-4-26b-a4b 5, deepseek-v4.1-flash 3, qwen3-235b-a22b-2507 3, glm-5.3-flash 1, glm-5.3 0, muse-glimmer-30b 3, kimi-k3 2, gpt-5.4-mini 10, gemini-3.1-flash-lite 7, claude-sonnet-5 4. Of the milestones E3 leaves unmatched, the share the judge rules not needed (and the count unmatched): gpt-oss-20b 0.314 (2,290), gemma-4-26b-a4b 0.462 (1,633), deepseek-v4.1-flash 0.608 (1,535), qwen3-235b-a22b-2507 0.546 (1,371), glm-5.3-flash 0.466 (1,475), glm-5.3 0.354 (1,595), muse-glimmer-30b 0.479 (1,642), kimi-k3 0.520 (1,398), gpt-5.4-mini 0.364 (2,182), gemini-3.1-flash-lite 0.445 (1,815), claude-sonnet-5 0.634 (1,110).

| model | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong: matching | wrong: route-adj. | wrong: interm. | wrong n |
|---|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | 0.815 [0.777, 0.851] | 0.793 [0.754, 0.830] | 0.875 [0.841, 0.905] | 0.761 [0.713, 0.807] | 0.461 | 0.451 | 0.506 | 0.527 | 318 |
| gemma-4-26b-a4b | 0.876 [0.848, 0.903] | 0.855 [0.825, 0.885] | 0.942 [0.922, 0.959] | 0.820 [0.776, 0.860] | 0.602 | 0.590 | 0.669 | 0.629 | 260 |
| deepseek-v4.1-flash | 0.900 [0.875, 0.924] | 0.862 [0.832, 0.891] | 0.982 [0.974, 0.988] | 0.826 [0.783, 0.868] | 0.687 | 0.607 | 0.799 | 0.678 | 21 |
| qwen3-235b-a22b-2507 | 0.894 [0.869, 0.918] | 0.871 [0.844, 0.897] | 0.957 [0.942, 0.969] | 0.837 [0.795, 0.876] | 0.629 | 0.592 | 0.681 | 0.678 | 176 |
| glm-5.3-flash | 0.912 [0.885, 0.936] | 0.882 [0.853, 0.909] | 0.966 [0.947, 0.981] | 0.861 [0.820, 0.900] | 0.477 | 0.477 | 0.513 | 0.535 | 8 |
| glm-5.3 | 0.901 [0.867, 0.931] | 0.876 [0.840, 0.908] | 0.947 [0.917, 0.972] | 0.854 [0.809, 0.895] | 0.845 | 0.810 | 0.881 | 0.783 | 7 |
| muse-glimmer-30b | 0.898 [0.872, 0.923] | 0.860 [0.830, 0.888] | 0.962 [0.944, 0.977] | 0.832 [0.789, 0.873] | 0.711 | 0.669 | 0.735 | 0.810 | 29 |
| kimi-k3 | 0.913 [0.889, 0.934] | 0.876 [0.850, 0.900] | 0.974 [0.961, 0.983] | 0.851 [0.811, 0.888] | 0.734 | 0.631 | 0.796 | 0.711 | 30 |
| gpt-5.4-mini | 0.835 [0.798, 0.870] | 0.798 [0.759, 0.837] | 0.900 [0.869, 0.929] | 0.762 [0.710, 0.811] | 0.485 | 0.442 | 0.539 | 0.530 | 309 |
| gemini-3.1-flash-lite | 0.865 [0.835, 0.893] | 0.840 [0.808, 0.871] | 0.932 [0.911, 0.951] | 0.802 [0.758, 0.843] | 0.542 | 0.515 | 0.613 | 0.590 | 262 |
| claude-sonnet-5 | 0.924 [0.903, 0.943] | 0.903 [0.880, 0.924] | 0.985 [0.978, 0.991] | 0.865 [0.824, 0.902] | 0.844 | 0.821 | 0.924 | 0.756 | 28 |

At the reasoning store `matched_config.json` names (same templates):

| model | store | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong n |
|---|---|---|---|---|---|---|---|
| gemma-4-26b-a4b | reasoning-medium-full | 0.885 [0.857, 0.911] | 0.852 [0.821, 0.881] | 0.959 [0.944, 0.973] | 0.815 [0.772, 0.854] | 0.583 | 124 |
| gpt-5.4-mini | reasoning-medium-full | 0.892 [0.866, 0.917] | 0.848 [0.816, 0.878] | 0.972 [0.961, 0.980] | 0.818 [0.772, 0.861] | 0.692 | 72 |
| gemini-3.1-flash-lite | reasoning-medium-full | 0.880 [0.850, 0.908] | 0.856 [0.824, 0.886] | 0.953 [0.934, 0.970] | 0.815 [0.770, 0.856] | 0.558 | 181 |

### The separation rule

Claude Sonnet 5 minus DeepSeek V4.1 Flash, Holm-adjusted p over each reading's 55 pairs: as scored (E3 + judge) +0.023, p = 0.0459; matching alone (E3) +0.040, p = 0.0005; route-adjusted +0.003, p = 1.0000; intermediate-only +0.038, p = 0.0251. 
**The separation does not hold**: the pair separates after Holm under as scored (E3 + judge) (p = 0.0459) and matching alone (E3) (p = 0.0005) but not under route-adjusted (p = 1.0000), so the results sentence does not claim it.

### The 55 pairs under each reading

Difference a minus b in the mean of template means, Holm-adjusted sign-flip p within each reading's family. Pairs separated at 0.05: as scored (E3 + judge) 28, matching alone (E3) 24, route-adjusted 29, intermediate-only 19.

| a | b | as scored (E3 + judge): diff | p Holm | matching alone (E3): diff | p Holm | route-adjusted: diff | p Holm | intermediate-only: diff | p Holm |
|---|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | gemma-4-26b-a4b | -0.061 | 0.0008* | -0.063 | 0.0005* | -0.067 | 0.0005* | -0.058 | 0.0619 |
| gpt-oss-20b | deepseek-v4.1-flash | -0.085 | 0.0005* | -0.070 | 0.0015* | -0.107 | 0.0005* | -0.065 | 0.0300* |
| gpt-oss-20b | qwen3-235b-a22b-2507 | -0.079 | 0.0005* | -0.079 | 0.0005* | -0.082 | 0.0005* | -0.075 | 0.0088* |
| gpt-oss-20b | glm-5.3-flash | -0.097 | 0.0005* | -0.090 | 0.0005* | -0.091 | 0.0005* | -0.100 | 0.0005* |
| gpt-oss-20b | glm-5.3 | -0.086 | 0.0008* | -0.084 | 0.0008* | -0.072 | 0.0042* | -0.092 | 0.0040* |
| gpt-oss-20b | muse-glimmer-30b | -0.083 | 0.0005* | -0.067 | 0.0005* | -0.087 | 0.0005* | -0.070 | 0.0078* |
| gpt-oss-20b | kimi-k3 | -0.097 | 0.0005* | -0.083 | 0.0005* | -0.099 | 0.0005* | -0.089 | 0.0014* |
| gpt-oss-20b | gpt-5.4-mini | -0.020 | 1.0000 | -0.005 | 1.0000 | -0.026 | 1.0000 | -0.000 | 1.0000 |
| gpt-oss-20b | gemini-3.1-flash-lite | -0.049 | 0.0321* | -0.047 | 0.0583 | -0.058 | 0.0066* | -0.040 | 0.9384 |
| gpt-oss-20b | claude-sonnet-5 | -0.108 | 0.0005* | -0.110 | 0.0005* | -0.110 | 0.0005* | -0.103 | 0.0005* |
| gemma-4-26b-a4b | deepseek-v4.1-flash | -0.024 | 0.3348 | -0.007 | 1.0000 | -0.040 | 0.0005* | -0.007 | 1.0000 |
| gemma-4-26b-a4b | qwen3-235b-a22b-2507 | -0.018 | 0.5260 | -0.016 | 1.0000 | -0.015 | 0.6257 | -0.017 | 1.0000 |
| gemma-4-26b-a4b | glm-5.3-flash | -0.036 | 0.0263* | -0.027 | 0.8103 | -0.024 | 0.2127 | -0.042 | 0.2749 |
| gemma-4-26b-a4b | glm-5.3 | -0.025 | 1.0000 | -0.021 | 1.0000 | -0.005 | 1.0000 | -0.034 | 1.0000 |
| gemma-4-26b-a4b | muse-glimmer-30b | -0.022 | 0.7132 | -0.004 | 1.0000 | -0.020 | 0.4779 | -0.012 | 1.0000 |
| gemma-4-26b-a4b | kimi-k3 | -0.037 | 0.0098* | -0.021 | 1.0000 | -0.032 | 0.0165* | -0.031 | 0.6507 |
| gemma-4-26b-a4b | gpt-5.4-mini | +0.041 | 0.0181* | +0.057 | 0.0008* | +0.041 | 0.0149* | +0.058 | 0.0071* |
| gemma-4-26b-a4b | gemini-3.1-flash-lite | +0.011 | 1.0000 | +0.015 | 0.9895 | +0.010 | 1.0000 | +0.018 | 1.0000 |
| gemma-4-26b-a4b | claude-sonnet-5 | -0.047 | 0.0008* | -0.048 | 0.0008* | -0.043 | 0.0005* | -0.045 | 0.1099 |
| deepseek-v4.1-flash | qwen3-235b-a22b-2507 | +0.006 | 1.0000 | -0.009 | 1.0000 | +0.025 | 0.0040* | -0.010 | 1.0000 |
| deepseek-v4.1-flash | glm-5.3-flash | -0.012 | 1.0000 | -0.020 | 1.0000 | +0.016 | 0.1065 | -0.035 | 0.6709 |
| deepseek-v4.1-flash | glm-5.3 | -0.001 | 1.0000 | -0.014 | 1.0000 | +0.035 | 0.0445* | -0.027 | 1.0000 |
| deepseek-v4.1-flash | muse-glimmer-30b | +0.002 | 1.0000 | +0.003 | 1.0000 | +0.020 | 0.1230 | -0.005 | 1.0000 |
| deepseek-v4.1-flash | kimi-k3 | -0.013 | 0.9526 | -0.013 | 1.0000 | +0.008 | 0.6093 | -0.024 | 0.6270 |
| deepseek-v4.1-flash | gpt-5.4-mini | +0.065 | 0.0005* | +0.064 | 0.0005* | +0.081 | 0.0005* | +0.065 | 0.0031* |
| deepseek-v4.1-flash | gemini-3.1-flash-lite | +0.035 | 0.0083* | +0.022 | 1.0000 | +0.050 | 0.0005* | +0.025 | 1.0000 |
| deepseek-v4.1-flash | claude-sonnet-5 | -0.023 | 0.0459* | -0.040 | 0.0005* | -0.003 | 1.0000 | -0.038 | 0.0251* |
| qwen3-235b-a22b-2507 | glm-5.3-flash | -0.018 | 1.0000 | -0.011 | 1.0000 | -0.009 | 1.0000 | -0.024 | 1.0000 |
| qwen3-235b-a22b-2507 | glm-5.3 | -0.007 | 1.0000 | -0.005 | 1.0000 | +0.010 | 1.0000 | -0.017 | 1.0000 |
| qwen3-235b-a22b-2507 | muse-glimmer-30b | -0.003 | 1.0000 | +0.012 | 1.0000 | -0.005 | 1.0000 | +0.005 | 1.0000 |
| qwen3-235b-a22b-2507 | kimi-k3 | -0.018 | 0.9526 | -0.004 | 1.0000 | -0.017 | 0.6257 | -0.014 | 1.0000 |
| qwen3-235b-a22b-2507 | gpt-5.4-mini | +0.060 | 0.0005* | +0.073 | 0.0005* | +0.056 | 0.0005* | +0.075 | 0.0005* |
| qwen3-235b-a22b-2507 | gemini-3.1-flash-lite | +0.030 | 0.0167* | +0.031 | 0.0603 | +0.025 | 0.0249* | +0.035 | 0.1595 |
| qwen3-235b-a22b-2507 | claude-sonnet-5 | -0.029 | 0.0321* | -0.031 | 0.0710 | -0.028 | 0.0007* | -0.028 | 0.9206 |
| glm-5.3-flash | glm-5.3 | +0.011 | 1.0000 | +0.006 | 1.0000 | +0.019 | 0.3080 | +0.007 | 1.0000 |
| glm-5.3-flash | muse-glimmer-30b | +0.014 | 1.0000 | +0.023 | 0.8468 | +0.004 | 1.0000 | +0.029 | 0.7917 |
| glm-5.3-flash | kimi-k3 | -0.001 | 1.0000 | +0.006 | 1.0000 | -0.008 | 1.0000 | +0.011 | 1.0000 |
| glm-5.3-flash | gpt-5.4-mini | +0.077 | 0.0005* | +0.084 | 0.0005* | +0.065 | 0.0005* | +0.100 | 0.0005* |
| glm-5.3-flash | gemini-3.1-flash-lite | +0.047 | 0.0005* | +0.042 | 0.0224* | +0.034 | 0.0048* | +0.059 | 0.0064* |
| glm-5.3-flash | claude-sonnet-5 | -0.012 | 1.0000 | -0.020 | 1.0000 | -0.019 | 0.3559 | -0.004 | 1.0000 |
| glm-5.3 | muse-glimmer-30b | +0.004 | 1.0000 | +0.017 | 1.0000 | -0.015 | 1.0000 | +0.022 | 1.0000 |
| glm-5.3 | kimi-k3 | -0.011 | 1.0000 | +0.000 | 1.0000 | -0.027 | 0.4264 | +0.003 | 1.0000 |
| glm-5.3 | gpt-5.4-mini | +0.067 | 0.0015* | +0.078 | 0.0005* | +0.047 | 0.0593 | +0.092 | 0.0028* |
| glm-5.3 | gemini-3.1-flash-lite | +0.037 | 0.2407 | +0.036 | 0.3716 | +0.015 | 1.0000 | +0.052 | 0.2749 |
| glm-5.3 | claude-sonnet-5 | -0.022 | 1.0000 | -0.027 | 1.0000 | -0.037 | 0.0900 | -0.011 | 1.0000 |
| muse-glimmer-30b | kimi-k3 | -0.015 | 1.0000 | -0.016 | 1.0000 | -0.012 | 1.0000 | -0.019 | 1.0000 |
| muse-glimmer-30b | gpt-5.4-mini | +0.063 | 0.0005* | +0.061 | 0.0005* | +0.061 | 0.0005* | +0.070 | 0.0005* |
| muse-glimmer-30b | gemini-3.1-flash-lite | +0.033 | 0.0259* | +0.020 | 1.0000 | +0.030 | 0.0141* | +0.030 | 0.7917 |
| muse-glimmer-30b | claude-sonnet-5 | -0.026 | 0.2163 | -0.043 | 0.0022* | -0.023 | 0.0745 | -0.033 | 0.4362 |
| kimi-k3 | gpt-5.4-mini | +0.078 | 0.0005* | +0.078 | 0.0005* | +0.073 | 0.0005* | +0.089 | 0.0005* |
| kimi-k3 | gemini-3.1-flash-lite | +0.048 | 0.0005* | +0.036 | 0.0467* | +0.042 | 0.0005* | +0.049 | 0.0097* |
| kimi-k3 | claude-sonnet-5 | -0.011 | 1.0000 | -0.027 | 0.0187* | -0.011 | 0.6093 | -0.014 | 1.0000 |
| gpt-5.4-mini | gemini-3.1-flash-lite | -0.030 | 0.1687 | -0.042 | 0.0147* | -0.032 | 0.0445* | -0.040 | 0.3081 |
| gpt-5.4-mini | claude-sonnet-5 | -0.089 | 0.0005* | -0.105 | 0.0005* | -0.084 | 0.0005* | -0.103 | 0.0005* |
| gemini-3.1-flash-lite | claude-sonnet-5 | -0.059 | 0.0005* | -0.063 | 0.0005* | -0.052 | 0.0005* | -0.063 | 0.0014* |

Smallest detectable difference (80% power, at the strictest Holm step), median over the 55 pairs: as scored (E3 + judge) 0.046, matching alone (E3) 0.050, route-adjusted 0.040, intermediate-only 0.065.

### Verbosity

Response-level least squares of MC as scored on the numeric values shown and the visible completion tokens, template fixed effects, pooled over 11 models (23,722 readable responses, 147 templates); 95% intervals resample templates. Per 10 numeric values: +0.0014 [-0.0002, +0.0033]; per 1,000 visible tokens: +0.0042 [-0.0066, +0.0152]. Within model and template (template-by-model fixed effects): +0.0021 [+0.0003, +0.0041] and -0.0080 [-0.0206, +0.0072]. Interquartile range of the count: 35 to 74; of visible tokens: 559 to 1,048. Spearman across the models' means (MC against the mean count): +0.309. Rows with more reasoning than completion tokens, floored at zero visible tokens: glm-5.3-flash 3, glm-5.3 18.

The within-model Spearman is over a model's readable responses with no template control, so it mixes in the templates' difficulty (a harder template draws a longer response and a lower MC); the regression's fixed effects remove that.

| model | MC (readable) | numeric values, mean | median | visible tokens, median | Spearman MC vs count | vs tokens |
|---|---|---|---|---|---|---|
| gpt-oss-20b | 0.836 | 56.6 | 49 | 736 | -0.190 | -0.243 |
| gemma-4-26b-a4b | 0.876 | 66.3 | 48 | 802 | -0.220 | -0.228 |
| deepseek-v4.1-flash | 0.903 | 52.5 | 45 | 626 | -0.109 | -0.134 |
| qwen3-235b-a22b-2507 | 0.895 | 101.3 | 63 | 1,022 | -0.146 | -0.164 |
| glm-5.3-flash | 0.930 | 60.4 | 52 | 672 | -0.067 | -0.113 |
| glm-5.3 | 0.943 | 64.9 | 58 | 741 | -0.126 | -0.174 |
| muse-glimmer-30b | 0.911 | 56.4 | 48 | 1,114 | -0.138 | -0.283 |
| kimi-k3 | 0.919 | 54.9 | 46 | 670 | -0.117 | -0.124 |
| gpt-5.4-mini | 0.837 | 54.0 | 44 | 636 | -0.207 | -0.212 |
| gemini-3.1-flash-lite | 0.865 | 49.9 | 43 | 654 | -0.235 | -0.233 |
| claude-sonnet-5 | 0.925 | 61.7 | 54 | 1,048 | -0.111 | -0.135 |

### MC against reasoning tokens (descriptive)

Quartiles of reasoning tokens within each model; MC as scored per quartile (item mean), the median reasoning tokens and the share of empty responses per quartile.

| model | store | Q1 | Q2 | Q3 | Q4 | median tokens Q1 to Q4 | empty Q1 to Q4 | Spearman |
|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | main | 0.930 | 0.872 | 0.827 | 0.632 | 298, 708, 1,423, 4,150 | 0.004, 0.000, 0.000, 0.095 | -0.401 |
| gemma-4-26b-a4b | reasoning-medium-full | 0.934 | 0.887 | 0.895 | 0.820 | 2,370, 4,052, 5,680, 9,127 | 0.000, 0.000, 0.000, 0.006 | -0.203 |
| deepseek-v4.1-flash | main | 0.930 | 0.903 | 0.916 | 0.848 | 636, 1,712, 3,351, 8,261 | 0.000, 0.000, 0.000, 0.015 | -0.197 |
| qwen3-235b-a22b-2507 | | no reasoning tokens in any store | | | | | | |
| glm-5.3-flash | main | 0.941 | 0.948 | 0.926 | 0.830 | 514, 1,248, 2,786, 8,746 | 0.000, 0.000, 0.000, 0.083 | -0.215 |
| glm-5.3 | main | 0.956 | 0.957 | 0.939 | 0.749 | 530, 1,652, 3,837, 13,955 | 0.000, 0.000, 0.000, 0.182 | -0.310 |
| muse-glimmer-30b | main | 0.941 | 0.931 | 0.915 | 0.804 | 1,082, 1,586, 2,152, 3,498 | 0.000, 0.000, 0.000, 0.057 | -0.259 |
| kimi-k3 | main | 0.941 | 0.926 | 0.915 | 0.867 | 484, 1,034, 2,195, 6,483 | 0.000, 0.000, 0.000, 0.028 | -0.194 |
| gpt-5.4-mini | reasoning-medium-full | 0.933 | 0.889 | 0.898 | 0.854 | 195, 437, 963, 2,967 | 0.000, 0.000, 0.000, 0.000 | -0.193 |
| gemini-3.1-flash-lite | reasoning-medium-full | 0.941 | 0.904 | 0.858 | 0.820 | 524, 732, 931, 1,310 | 0.000, 0.000, 0.000, 0.000 | -0.258 |
| claude-sonnet-5 | main | 0.945 | 0.936 | 0.910 | 0.907 | 160, 472, 1,048, 2,921 | 0.000, 0.000, 0.000, 0.004 | -0.169 |
