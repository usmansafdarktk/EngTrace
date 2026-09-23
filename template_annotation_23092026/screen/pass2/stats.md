# Pass 2 - the numbers

Computed by `analyze_screen.py` from `replies.jsonl` (450 rows, 450 ok). Judges: grok-4.6, minimax-m3, mimo-v2.5-pro. Corpus at git `fc1a6dcb35`, prompt template `93a395de04`, instance seeds [1001, 1002, 1003].

## Certification outcome (paper rule, section 3.3)

| Category | Templates |
|---|---:|
| Pass | 147 |
| Controversial | 3 |
| Critical Failure | 0 |
| judged by all 3 | 150 |
| incomplete (a judge failed) | 0 |

| Branch | Pass | Controversial | Critical Failure |
|---|---:|---:|---:|
| chemical_engineering | 30 | 0 | 0 |
| civil_engineering | 29 | 1 | 0 |
| electrical_engineering | 29 | 1 | 0 |
| industrial_engineering | 30 | 0 | 0 |
| mechanical_engineering | 29 | 1 | 0 |

## Disagreement among the judges (sigma_max, equation 2)

| Statistic | Value |
|---|---:|
| mean sigma_max | 0.219 |
| median sigma_max | 0.000 |
| max sigma_max | 0.943 |
| templates with sigma_max <= 0.5 | 147 of 150 |

## Inter-judge agreement on the review flag (binary)

All 3 judges: Gwet AC1 0.929, Fleiss kappa -0.034, percent agreement 90.0%.

| Pair | Gwet AC1 | Cohen/Fleiss kappa | Percent agreement |
|---|---:|---:|---:|
| grok-4.6 vs minimax-m3 | 0.890 | -0.053 | 90.0% |
| grok-4.6 vs mimo-v2.5-pro | 0.980 | -0.010 | 98.0% |
| minimax-m3 vs mimo-v2.5-pro | 0.913 | -0.042 | 92.0% |

Per judge: flag rate and mean scores.

| Judge | Flag rate | mean phys | mean math | mean ped | mean confidence |
|---|---:|---:|---:|---:|---:|
| grok-4.6 | 2.0% | 4.90 | 4.84 | 4.81 | 4.92 |
| minimax-m3 | 8.0% | 4.94 | 4.91 | 4.75 | 4.91 |
| mimo-v2.5-pro | 0.0% | 4.99 | 4.99 | 4.97 | 4.99 |

## Agreement on the 1-5 scores

| Dimension | All three identical | Gwet AC1 (5 categories, unweighted) |
|---|---:|---:|
| phys | 88.0% | 0.918 |
| math | 80.7% | 0.866 |
| ped | 66.0% | 0.755 |

## The run itself

Cost: $0.704 for 72 judged rows (72 priced by OpenRouter's own usage.cost, the rest from the catalogue). 378 rows were carried forward from pass 1 because the template's prompt text is byte-identical there; their original cost was $3.402.

| Judge | Served model id | Rows |
|---|---|---:|
| grok-4.6 | x-ai/grok-4.6 | 150 |
| mimo-v2.5-pro | xiaomi/mimo-v2.5-pro | 150 |
| minimax-m3 | minimax/minimax-m3 | 150 |

| Judge | Attempts needed | Rows |
|---|---:|---:|
| grok-4.6 | 1 | 150 |
| mimo-v2.5-pro | 1 | 150 |
| minimax-m3 | 1 | 147 |
| minimax-m3 | 2 | 3 |
