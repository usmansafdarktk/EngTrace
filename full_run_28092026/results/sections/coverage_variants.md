## Coverage variants

Printed by `coverage_variants.py`; the method is in its docstring. MC as scored is Table 1's Milestone Coverage, recomputed from the judge's per-milestone sources and checked equal to the stored `e5_strict` for every response. Each model's value is the mean of its template means over the templates on which every model has a value under that reading, with a 95% template-bootstrap interval; "wrong" is the item mean over readable responses scoring 0.

**Target rule.** A milestone is an answer target when a number in the answer row's targets.numbers equals its value under one of the answer check's unit factors (answer.SCALES) within the milestone display tolerance (0.5%), the test answer.targets uses to make a gold answer value a target; milestones.py records no answer line. Of 2,250 instances, 57 have no milestones and 445 have only target milestones, so 502 are left without milestones under intermediate-only MC (23 templates have none left). 2,750 of 8,779 milestones are targets; on 159 instances the gold's answer value is not a milestone, so all their milestones stay. Templates in each reading's comparison: as scored (E3 + judge) 148, matching alone (E3) 148, route-adjusted 148, intermediate-only 127.

**Route-adjusted.** (e3 + REACHED) / (n - NOT_NEEDED); a response whose milestones the judge rules all not needed has no value and is left out (count per model): gpt-oss-20b 7, gemma-4-26b-a4b 5, deepseek-v4.1-flash 3, qwen3-235b-a22b-2507 4, glm-5.3-flash 1, glm-5.3 0, muse-glimmer-30b 4, kimi-k3 2, gpt-5.4-mini 12, gemini-3.1-flash-lite 7, claude-sonnet-5 4. Of the milestones E3 leaves unmatched, the share the judge rules not needed (and the count unmatched): gpt-oss-20b 0.314 (2,307), gemma-4-26b-a4b 0.464 (1,641), deepseek-v4.1-flash 0.608 (1,549), qwen3-235b-a22b-2507 0.547 (1,390), glm-5.3-flash 0.493 (1,401), glm-5.3 0.363 (1,556), muse-glimmer-30b 0.483 (1,637), kimi-k3 0.522 (1,396), gpt-5.4-mini 0.366 (2,197), gemini-3.1-flash-lite 0.446 (1,821), claude-sonnet-5 0.643 (1,110).

| model | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong: matching | wrong: route-adj. | wrong: interm. | wrong n |
|---|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | 0.812 [0.775, 0.847] | 0.790 [0.752, 0.826] | 0.872 [0.840, 0.902] | 0.759 [0.711, 0.805] | 0.465 | 0.455 | 0.509 | 0.534 | 335 |
| gemma-4-26b-a4b | 0.876 [0.847, 0.903] | 0.855 [0.825, 0.884] | 0.942 [0.922, 0.959] | 0.817 [0.773, 0.858] | 0.602 | 0.588 | 0.667 | 0.627 | 260 |
| deepseek-v4.1-flash | 0.899 [0.874, 0.922] | 0.861 [0.831, 0.890] | 0.982 [0.974, 0.988] | 0.822 [0.778, 0.864] | 0.667 | 0.667 | 0.711 | 0.644 | 12 |
| qwen3-235b-a22b-2507 | 0.892 [0.867, 0.916] | 0.866 [0.836, 0.894] | 0.956 [0.942, 0.969] | 0.834 [0.792, 0.873] | 0.616 | 0.586 | 0.675 | 0.673 | 173 |
| glm-5.3-flash | 0.912 [0.886, 0.936] | 0.891 [0.864, 0.916] | 0.966 [0.947, 0.981] | 0.859 [0.818, 0.898] | 0.458 | 0.378 | 0.524 | 0.462 | 12 |
| glm-5.3 | 0.900 [0.866, 0.930] | 0.881 [0.846, 0.912] | 0.946 [0.916, 0.971] | 0.852 [0.807, 0.893] | 0.596 | 0.486 | 0.596 | 0.825 | 8 |
| muse-glimmer-30b | 0.897 [0.871, 0.922] | 0.862 [0.832, 0.890] | 0.962 [0.944, 0.977] | 0.830 [0.787, 0.871] | 0.653 | 0.583 | 0.681 | 0.788 | 33 |
| kimi-k3 | 0.912 [0.888, 0.934] | 0.877 [0.850, 0.902] | 0.973 [0.961, 0.983] | 0.847 [0.806, 0.885] | 0.584 | 0.584 | 0.627 | 0.646 | 20 |
| gpt-5.4-mini | 0.834 [0.798, 0.868] | 0.797 [0.758, 0.834] | 0.901 [0.870, 0.929] | 0.760 [0.708, 0.810] | 0.485 | 0.442 | 0.537 | 0.531 | 311 |
| gemini-3.1-flash-lite | 0.864 [0.835, 0.892] | 0.840 [0.808, 0.870] | 0.932 [0.911, 0.951] | 0.800 [0.757, 0.841] | 0.537 | 0.511 | 0.605 | 0.588 | 257 |
| claude-sonnet-5 | 0.922 [0.902, 0.942] | 0.903 [0.880, 0.925] | 0.984 [0.978, 0.990] | 0.859 [0.817, 0.897] | 0.869 | 0.850 | 0.935 | 0.834 | 20 |

At the reasoning store `matched_config.json` names (same templates):

| model | store | as scored | matching alone | route-adjusted | intermediate-only | wrong: as scored | wrong n |
|---|---|---|---|---|---|---|---|
| gemma-4-26b-a4b | reasoning-medium-full | 0.884 [0.856, 0.910] | 0.850 [0.819, 0.880] | 0.959 [0.944, 0.973] | 0.813 [0.770, 0.853] | 0.569 | 116 |
| gpt-5.4-mini | reasoning-medium-full | 0.892 [0.866, 0.917] | 0.847 [0.815, 0.877] | 0.972 [0.962, 0.981] | 0.817 [0.771, 0.860] | 0.689 | 61 |
| gemini-3.1-flash-lite | reasoning-medium-full | 0.881 [0.851, 0.908] | 0.855 [0.824, 0.886] | 0.953 [0.935, 0.969] | 0.815 [0.770, 0.856] | 0.540 | 174 |

### The separation rule

Claude Sonnet 5 minus DeepSeek V4.1 Flash, Holm-adjusted p over each reading's 55 pairs: as scored (E3 + judge) +0.023, p = 0.0493; matching alone (E3) +0.042, p = 0.0005; route-adjusted +0.003, p = 1.0000; intermediate-only +0.037, p = 0.0269. 
**The separation does not hold**: the pair separates after Holm under as scored (E3 + judge) (p = 0.0493) and matching alone (E3) (p = 0.0005) but not under route-adjusted (p = 1.0000), so the results sentence does not claim it.

### The 55 pairs under each reading

Difference a minus b in the mean of template means, Holm-adjusted sign-flip p within each reading's family. Pairs separated at 0.05: as scored (E3 + judge) 28, matching alone (E3) 26, route-adjusted 29, intermediate-only 18.

| a | b | as scored (E3 + judge): diff | p Holm | matching alone (E3): diff | p Holm | route-adjusted: diff | p Holm | intermediate-only: diff | p Holm |
|---|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | gemma-4-26b-a4b | -0.064 | 0.0005* | -0.066 | 0.0005* | -0.070 | 0.0005* | -0.058 | 0.0850 |
| gpt-oss-20b | deepseek-v4.1-flash | -0.087 | 0.0005* | -0.071 | 0.0005* | -0.110 | 0.0005* | -0.063 | 0.0673 |
| gpt-oss-20b | qwen3-235b-a22b-2507 | -0.080 | 0.0005* | -0.076 | 0.0005* | -0.084 | 0.0005* | -0.075 | 0.0090* |
| gpt-oss-20b | glm-5.3-flash | -0.100 | 0.0005* | -0.102 | 0.0005* | -0.094 | 0.0005* | -0.100 | 0.0005* |
| gpt-oss-20b | glm-5.3 | -0.088 | 0.0005* | -0.092 | 0.0005* | -0.073 | 0.0037* | -0.093 | 0.0037* |
| gpt-oss-20b | muse-glimmer-30b | -0.085 | 0.0005* | -0.072 | 0.0005* | -0.090 | 0.0005* | -0.071 | 0.0066* |
| gpt-oss-20b | kimi-k3 | -0.100 | 0.0005* | -0.087 | 0.0005* | -0.101 | 0.0005* | -0.088 | 0.0014* |
| gpt-oss-20b | gpt-5.4-mini | -0.021 | 1.0000 | -0.007 | 1.0000 | -0.028 | 1.0000 | -0.001 | 1.0000 |
| gpt-oss-20b | gemini-3.1-flash-lite | -0.052 | 0.0231* | -0.050 | 0.0352* | -0.060 | 0.0053* | -0.041 | 0.9053 |
| gpt-oss-20b | claude-sonnet-5 | -0.110 | 0.0005* | -0.113 | 0.0005* | -0.112 | 0.0005* | -0.100 | 0.0005* |
| gemma-4-26b-a4b | deepseek-v4.1-flash | -0.023 | 0.4716 | -0.006 | 1.0000 | -0.040 | 0.0005* | -0.005 | 1.0000 |
| gemma-4-26b-a4b | qwen3-235b-a22b-2507 | -0.016 | 0.8812 | -0.011 | 1.0000 | -0.014 | 0.7124 | -0.017 | 1.0000 |
| gemma-4-26b-a4b | glm-5.3-flash | -0.036 | 0.0246* | -0.036 | 0.0536 | -0.024 | 0.2124 | -0.042 | 0.2330 |
| gemma-4-26b-a4b | glm-5.3 | -0.024 | 1.0000 | -0.026 | 1.0000 | -0.004 | 1.0000 | -0.035 | 1.0000 |
| gemma-4-26b-a4b | muse-glimmer-30b | -0.021 | 0.7516 | -0.007 | 1.0000 | -0.020 | 0.4864 | -0.013 | 1.0000 |
| gemma-4-26b-a4b | kimi-k3 | -0.036 | 0.0136* | -0.021 | 1.0000 | -0.031 | 0.0195* | -0.030 | 0.9053 |
| gemma-4-26b-a4b | gpt-5.4-mini | +0.042 | 0.0133* | +0.059 | 0.0005* | +0.042 | 0.0138* | +0.057 | 0.0132* |
| gemma-4-26b-a4b | gemini-3.1-flash-lite | +0.011 | 1.0000 | +0.015 | 0.8515 | +0.010 | 1.0000 | +0.017 | 1.0000 |
| gemma-4-26b-a4b | claude-sonnet-5 | -0.046 | 0.0008* | -0.047 | 0.0008* | -0.042 | 0.0005* | -0.042 | 0.2183 |
| deepseek-v4.1-flash | qwen3-235b-a22b-2507 | +0.007 | 1.0000 | -0.005 | 1.0000 | +0.025 | 0.0025* | -0.013 | 1.0000 |
| deepseek-v4.1-flash | glm-5.3-flash | -0.013 | 1.0000 | -0.030 | 0.1576 | +0.016 | 0.1003 | -0.038 | 0.4920 |
| deepseek-v4.1-flash | glm-5.3 | -0.001 | 1.0000 | -0.020 | 1.0000 | +0.036 | 0.0264* | -0.030 | 1.0000 |
| deepseek-v4.1-flash | muse-glimmer-30b | +0.002 | 1.0000 | -0.001 | 1.0000 | +0.020 | 0.1379 | -0.008 | 1.0000 |
| deepseek-v4.1-flash | kimi-k3 | -0.013 | 0.8413 | -0.016 | 0.5512 | +0.008 | 0.5461 | -0.025 | 0.4921 |
| deepseek-v4.1-flash | gpt-5.4-mini | +0.065 | 0.0005* | +0.064 | 0.0005* | +0.081 | 0.0005* | +0.061 | 0.0054* |
| deepseek-v4.1-flash | gemini-3.1-flash-lite | +0.034 | 0.0122* | +0.021 | 1.0000 | +0.050 | 0.0005* | +0.022 | 1.0000 |
| deepseek-v4.1-flash | claude-sonnet-5 | -0.023 | 0.0493* | -0.042 | 0.0005* | -0.003 | 1.0000 | -0.037 | 0.0269* |
| qwen3-235b-a22b-2507 | glm-5.3-flash | -0.020 | 1.0000 | -0.025 | 1.0000 | -0.010 | 1.0000 | -0.025 | 1.0000 |
| qwen3-235b-a22b-2507 | glm-5.3 | -0.008 | 1.0000 | -0.015 | 1.0000 | +0.011 | 1.0000 | -0.018 | 1.0000 |
| qwen3-235b-a22b-2507 | muse-glimmer-30b | -0.005 | 1.0000 | +0.004 | 1.0000 | -0.006 | 1.0000 | +0.005 | 1.0000 |
| qwen3-235b-a22b-2507 | kimi-k3 | -0.020 | 0.7516 | -0.011 | 1.0000 | -0.017 | 0.6009 | -0.013 | 1.0000 |
| qwen3-235b-a22b-2507 | gpt-5.4-mini | +0.059 | 0.0008* | +0.069 | 0.0008* | +0.056 | 0.0007* | +0.074 | 0.0005* |
| qwen3-235b-a22b-2507 | gemini-3.1-flash-lite | +0.028 | 0.0380* | +0.026 | 0.5935 | +0.025 | 0.0361* | +0.034 | 0.2005 |
| qwen3-235b-a22b-2507 | claude-sonnet-5 | -0.030 | 0.0246* | -0.036 | 0.0356* | -0.028 | 0.0005* | -0.024 | 1.0000 |
| glm-5.3-flash | glm-5.3 | +0.012 | 1.0000 | +0.010 | 1.0000 | +0.020 | 0.1476 | +0.007 | 1.0000 |
| glm-5.3-flash | muse-glimmer-30b | +0.015 | 1.0000 | +0.030 | 0.0832 | +0.004 | 1.0000 | +0.030 | 0.7868 |
| glm-5.3-flash | kimi-k3 | +0.000 | 1.0000 | +0.015 | 1.0000 | -0.008 | 1.0000 | +0.012 | 1.0000 |
| glm-5.3-flash | gpt-5.4-mini | +0.078 | 0.0005* | +0.095 | 0.0005* | +0.065 | 0.0005* | +0.099 | 0.0005* |
| glm-5.3-flash | gemini-3.1-flash-lite | +0.047 | 0.0008* | +0.051 | 0.0008* | +0.034 | 0.0028* | +0.059 | 0.0073* |
| glm-5.3-flash | claude-sonnet-5 | -0.010 | 1.0000 | -0.011 | 1.0000 | -0.018 | 0.3613 | +0.001 | 1.0000 |
| glm-5.3 | muse-glimmer-30b | +0.003 | 1.0000 | +0.019 | 1.0000 | -0.016 | 1.0000 | +0.022 | 1.0000 |
| glm-5.3 | kimi-k3 | -0.012 | 1.0000 | +0.004 | 1.0000 | -0.028 | 0.2937 | +0.005 | 1.0000 |
| glm-5.3 | gpt-5.4-mini | +0.066 | 0.0018* | +0.084 | 0.0005* | +0.045 | 0.0667 | +0.092 | 0.0028* |
| glm-5.3 | gemini-3.1-flash-lite | +0.036 | 0.2792 | +0.041 | 0.0910 | +0.014 | 1.0000 | +0.052 | 0.2742 |
| glm-5.3 | claude-sonnet-5 | -0.022 | 1.0000 | -0.021 | 1.0000 | -0.039 | 0.0624 | -0.007 | 1.0000 |
| muse-glimmer-30b | kimi-k3 | -0.014 | 1.0000 | -0.015 | 1.0000 | -0.011 | 1.0000 | -0.017 | 1.0000 |
| muse-glimmer-30b | gpt-5.4-mini | +0.064 | 0.0005* | +0.065 | 0.0008* | +0.061 | 0.0005* | +0.069 | 0.0005* |
| muse-glimmer-30b | gemini-3.1-flash-lite | +0.033 | 0.0330* | +0.022 | 0.9082 | +0.030 | 0.0158* | +0.030 | 0.9053 |
| muse-glimmer-30b | claude-sonnet-5 | -0.025 | 0.2792 | -0.041 | 0.0031* | -0.022 | 0.0811 | -0.029 | 0.9053 |
| kimi-k3 | gpt-5.4-mini | +0.078 | 0.0005* | +0.080 | 0.0005* | +0.073 | 0.0005* | +0.087 | 0.0005* |
| kimi-k3 | gemini-3.1-flash-lite | +0.047 | 0.0005* | +0.037 | 0.0356* | +0.042 | 0.0005* | +0.047 | 0.0285* |
| kimi-k3 | claude-sonnet-5 | -0.010 | 1.0000 | -0.026 | 0.0290* | -0.011 | 0.6009 | -0.012 | 1.0000 |
| gpt-5.4-mini | gemini-3.1-flash-lite | -0.031 | 0.1320 | -0.043 | 0.0082* | -0.031 | 0.0435* | -0.040 | 0.3227 |
| gpt-5.4-mini | claude-sonnet-5 | -0.088 | 0.0005* | -0.106 | 0.0005* | -0.084 | 0.0005* | -0.098 | 0.0005* |
| gemini-3.1-flash-lite | claude-sonnet-5 | -0.058 | 0.0005* | -0.063 | 0.0005* | -0.052 | 0.0005* | -0.059 | 0.0073* |

Smallest detectable difference (80% power, at the strictest Holm step), median over the 55 pairs: as scored (E3 + judge) 0.046, matching alone (E3) 0.050, route-adjusted 0.040, intermediate-only 0.066.

### Verbosity

Response-level least squares of MC as scored on the numeric values shown and the visible completion tokens, template fixed effects, pooled over 11 models (23,865 readable responses, 148 templates); 95% intervals resample templates. Per 10 numeric values: +0.0015 [-0.0001, +0.0034]; per 1,000 visible tokens: +0.0036 [-0.0067, +0.0145]. Within model and template (template-by-model fixed effects): +0.0023 [+0.0005, +0.0043] and -0.0094 [-0.0220, +0.0049]. Interquartile range of the count: 35 to 74; of visible tokens: 560 to 1,048. Spearman across the models' means (MC against the mean count): +0.309. Rows with more reasoning than completion tokens, floored at zero visible tokens: glm-5.3-flash 3, glm-5.3 18.

The within-model Spearman is over a model's readable responses with no template control, so it mixes in the templates' difficulty (a harder template draws a longer response and a lower MC); the regression's fixed effects remove that.

| model | MC (readable) | numeric values, mean | median | visible tokens, median | Spearman MC vs count | vs tokens |
|---|---|---|---|---|---|---|
| gpt-oss-20b | 0.833 | 56.5 | 49 | 739 | -0.190 | -0.249 |
| gemma-4-26b-a4b | 0.875 | 66.1 | 48 | 802 | -0.222 | -0.228 |
| deepseek-v4.1-flash | 0.901 | 52.5 | 45 | 628 | -0.114 | -0.139 |
| qwen3-235b-a22b-2507 | 0.892 | 101.3 | 63 | 1,026 | -0.147 | -0.168 |
| glm-5.3-flash | 0.930 | 60.3 | 52 | 673 | -0.069 | -0.110 |
| glm-5.3 | 0.941 | 64.8 | 58 | 744 | -0.133 | -0.178 |
| muse-glimmer-30b | 0.910 | 56.3 | 48 | 1,120 | -0.140 | -0.285 |
| kimi-k3 | 0.918 | 54.8 | 46 | 670 | -0.125 | -0.132 |
| gpt-5.4-mini | 0.835 | 53.9 | 44 | 639 | -0.207 | -0.211 |
| gemini-3.1-flash-lite | 0.865 | 49.8 | 43 | 655 | -0.239 | -0.235 |
| claude-sonnet-5 | 0.923 | 61.6 | 53 | 1,049 | -0.113 | -0.139 |

### MC against reasoning tokens (descriptive)

Quartiles of reasoning tokens within each model; MC as scored per quartile (item mean), the median reasoning tokens and the share of empty responses per quartile.

| model | store | Q1 | Q2 | Q3 | Q4 | median tokens Q1 to Q4 | empty Q1 to Q4 | Spearman |
|---|---|---|---|---|---|---|---|---|
| gpt-oss-20b | main | 0.930 | 0.862 | 0.829 | 0.629 | 298, 710, 1,425, 4,150 | 0.004, 0.000, 0.000, 0.095 | -0.402 |
| gemma-4-26b-a4b | reasoning-medium-full | 0.934 | 0.886 | 0.894 | 0.821 | 2,389, 4,064, 5,680, 9,114 | 0.000, 0.000, 0.000, 0.005 | -0.201 |
| deepseek-v4.1-flash | main | 0.931 | 0.902 | 0.913 | 0.846 | 640, 1,722, 3,338, 8,239 | 0.000, 0.000, 0.000, 0.015 | -0.200 |
| qwen3-235b-a22b-2507 | | no reasoning tokens in any store | | | | | | |
| glm-5.3-flash | main | 0.942 | 0.946 | 0.927 | 0.829 | 516, 1,251, 2,784, 8,684 | 0.000, 0.000, 0.000, 0.082 | -0.215 |
| glm-5.3 | main | 0.957 | 0.955 | 0.936 | 0.747 | 535, 1,656, 3,822, 13,932 | 0.000, 0.000, 0.000, 0.181 | -0.315 |
| muse-glimmer-30b | main | 0.942 | 0.929 | 0.915 | 0.803 | 1,084, 1,587, 2,150, 3,498 | 0.000, 0.000, 0.000, 0.057 | -0.261 |
| kimi-k3 | main | 0.942 | 0.926 | 0.913 | 0.865 | 485, 1,034, 2,180, 6,478 | 0.000, 0.000, 0.000, 0.027 | -0.199 |
| gpt-5.4-mini | reasoning-medium-full | 0.929 | 0.892 | 0.898 | 0.855 | 196, 437, 960, 2,952 | 0.000, 0.000, 0.000, 0.000 | -0.193 |
| gemini-3.1-flash-lite | reasoning-medium-full | 0.941 | 0.906 | 0.858 | 0.821 | 525, 733, 930, 1,308 | 0.000, 0.000, 0.000, 0.000 | -0.257 |
| claude-sonnet-5 | main | 0.945 | 0.936 | 0.905 | 0.905 | 160, 473, 1,044, 2,908 | 0.000, 0.000, 0.000, 0.004 | -0.174 |
