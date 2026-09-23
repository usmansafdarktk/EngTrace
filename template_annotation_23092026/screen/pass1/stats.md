# Pass 1 - the numbers

Computed by `analyze_screen.py` from `replies.jsonl` (450 rows, 450 ok). Judges: grok-4.6, minimax-m3, mimo-v2.5-pro. Corpus at git `2e3480ca43`, prompt template `93a395de04`, instance seeds [1001, 1002, 1003].

## Certification outcome (paper rule, section 3.3)

| Category | Templates |
|---|---:|
| Pass | 126 |
| Controversial | 13 |
| Critical Failure | 11 |
| judged by all 3 | 150 |
| incomplete (a judge failed) | 0 |

| Branch | Pass | Controversial | Critical Failure |
|---|---:|---:|---:|
| chemical_engineering | 21 | 4 | 5 |
| civil_engineering | 29 | 1 | 0 |
| electrical_engineering | 23 | 3 | 4 |
| industrial_engineering | 29 | 0 | 1 |
| mechanical_engineering | 24 | 5 | 1 |

## Disagreement among the judges (sigma_max, equation 2)

| Statistic | Value |
|---|---:|
| mean sigma_max | 0.293 |
| median sigma_max | 0.000 |
| max sigma_max | 0.943 |
| templates with sigma_max <= 0.5 | 128 of 150 |

## Inter-judge agreement on the review flag (binary)

All 3 judges: Gwet AC1 0.836, Fleiss kappa 0.287, percent agreement 80.0%.

| Pair | Gwet AC1 | Cohen/Fleiss kappa | Percent agreement |
|---|---:|---:|---:|
| grok-4.6 vs minimax-m3 | 0.799 | 0.350 | 84.7% |
| grok-4.6 vs mimo-v2.5-pro | 0.840 | 0.212 | 86.7% |
| minimax-m3 vs mimo-v2.5-pro | 0.866 | 0.258 | 88.7% |

Per judge: flag rate and mean scores.

| Judge | Flag rate | mean phys | mean math | mean ped | mean confidence |
|---|---:|---:|---:|---:|---:|
| grok-4.6 | 14.7% | 4.83 | 4.71 | 4.69 | 4.91 |
| minimax-m3 | 12.7% | 4.91 | 4.89 | 4.67 | 4.88 |
| mimo-v2.5-pro | 4.0% | 4.98 | 4.99 | 4.89 | 4.98 |

## Agreement on the 1-5 scores

| Dimension | All three identical | Gwet AC1 (5 categories, unweighted) |
|---|---:|---:|
| phys | 86.0% | 0.892 |
| math | 76.0% | 0.820 |
| ped | 61.3% | 0.710 |

## The run itself

Cost: $4.147 for 450 judged rows (450 priced by OpenRouter's own usage.cost, the rest from the catalogue).

| Judge | Served model id | Rows |
|---|---|---:|
| grok-4.6 | x-ai/grok-4.6 | 150 |
| mimo-v2.5-pro | xiaomi/mimo-v2.5-pro | 150 |
| minimax-m3 | minimax/minimax-m3 | 150 |

| Judge | Attempts needed | Rows |
|---|---:|---:|
| grok-4.6 | 1 | 150 |
| mimo-v2.5-pro | 1 | 150 |
| minimax-m3 | 1 | 145 |
| minimax-m3 | 2 | 5 |
