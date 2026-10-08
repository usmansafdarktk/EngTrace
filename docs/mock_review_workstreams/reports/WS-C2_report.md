# WS-C2 report: Milestone Coverage variants

No signal (C2 posts none). Written 2026-10-07. Every number below is from `full_run_28092026/results/coverage_variants.json`
at full resolution (10,000 bootstrap draws, 100,000 sign flips) on the current stores (main, and `reasoning-medium-full`
for the three re-run models); the section `results/sections/coverage_variants.md` prints the same with intervals.

## Script and flags

| script | what it computes | output | check |
|---|---|---|---|
| `full_run_28092026/coverage_variants.py` (new) | four readings of MC per response; per-model means (all, wrong answers) at the default and the reasoning store; the 55-pair family under each reading; the separation rule; the verbosity model; MC by reasoning-token quartile | `results/coverage_variants.json`, `results/sections/coverage_variants.md` | `--selftest` (alias `--check`) passes |

- `python -m full_run_28092026.coverage_variants` is the full pass (about 7 minutes). `--quick` runs 1,000 draws and 10,000 flips (about 4 minutes), and both files then say QUICK. `--render` rewrites the section from the JSON. `--selftest` takes about 1 minute.
- What the self-test checks:
  - MC as scored equals the stored `e5_strict` and `analyze.coverage_value` on all 30,520 responses with milestones, across 14 model stores (11 default, 3 reasoning).
  - The four readings on a hand-made response, and its edge cases: every milestone not needed, every milestone a target, no reply.
  - The number reader.
  - The fixed-effects fit recovers known slopes, and a bootstrap resample solved from cross-products equals a direct fit.
  - The quartiles.
- **Cross-check.** As scored, the family reproduces `analyze.q3_coverage` run on the current main store with no changes: Table 1's MC and its interval for every model, all 28 separated pairs, and the Claude–DeepSeek pair's difference, interval, p and Holm p, bit for bit.
- **Imports from `analyze.py`:** `load`, `load_stage`, `per_template`, `coverage_value`, `boot_mean`, `cluster_mean`, `sign_flip_p(d, seed, draws=)`, `holm`, `detectable_paired`, `pct`, and the globals `B` and `B_TEST`. Quick mode sets `analyze.B = 1000` in its own process.

## The rules

- **As scored:** (e3 + REACHED) / n, judge.py's `e5_strict`.
- **Matching alone:** E3's coverage.
- **Route-adjusted:** (e3 + REACHED) / (n − NOT_NEEDED). MISSING and UNJUDGED stay in the denominator. A response whose milestones the judge rules all not needed has no value and is left out: 0 to 10 per model, 45 in all.
- **Target rule (intermediate-only).** `milestones.py` records no answer line, so the rule inverts `answer.targets`. A milestone is an answer target when a number in the answer row's `targets.numbers` equals its value:
  - under one of the answer check's unit factors (`answer.SCALES`);
  - within the milestone display tolerance (0.5%), which is the test `answer.targets` uses to make a gold answer value a target.

  The targets are the same in every row of an item (asserted).
- **Count left without milestones: 510 of 2,250 instances.** 70 have no milestones at all, and 440 have only target milestones. 23 templates have none left, so the intermediate-only comparison runs on 127 templates; the other three readings run on 147.
  - 2,736 of 8,754 milestones are targets.
  - On 159 instances (40 templates) the gold's answer value is not a milestone, so all their milestones stay.
- **No failed replies.** No response in any store has a judge call without a reply, so all four readings use the same responses (2,180 per model).

## Findings on the current stores

Per model: the mean of template means with the template bootstrap. As scored is Table 1's MC. Wrong answers are readable responses that score 0.

| model | as scored | matching alone | route-adjusted | intermediate-only | wrong answers, as scored (n) |
|---|---|---|---|---|---|
| gpt-oss-20b | 0.815 | 0.793 | 0.875 | 0.761 | 0.461 (318) |
| gemma-4-26b-a4b | 0.876 | 0.855 | 0.942 | 0.820 | 0.602 (260) |
| deepseek-v4.1-flash | 0.900 | 0.862 | 0.982 | 0.826 | 0.687 (21) |
| qwen3-235b-a22b-2507 | 0.894 | 0.871 | 0.957 | 0.837 | 0.629 (176) |
| glm-5.3-flash | 0.912 | 0.882 | 0.966 | 0.861 | 0.477 (8) |
| glm-5.3 | 0.901 | 0.876 | 0.947 | 0.854 | 0.845 (7) |
| muse-glimmer-30b | 0.898 | 0.860 | 0.962 | 0.832 | 0.711 (29) |
| kimi-k3 | 0.913 | 0.876 | 0.974 | 0.851 | 0.734 (30) |
| gpt-5.4-mini | 0.835 | 0.798 | 0.900 | 0.762 | 0.485 (309) |
| gemini-3.1-flash-lite | 0.865 | 0.840 | 0.932 | 0.802 | 0.542 (262) |
| claude-sonnet-5 | 0.924 | 0.903 | 0.985 | 0.865 | 0.844 (28) |

At the reasoning store, over the same templates:

| model | as scored | matching alone | route-adjusted | intermediate-only |
|---|---|---|---|---|
| gemma-4-26b-a4b | 0.885 | 0.852 | 0.959 | 0.815 |
| gpt-5.4-mini | 0.892 | 0.848 | 0.972 | 0.818 |
| gemini-3.1-flash-lite | 0.880 | 0.856 | 0.953 | 0.815 |

**The separation rule does not hold.** Claude Sonnet 5 minus DeepSeek V4.1 Flash, Holm over each reading's 55 pairs:

| reading | difference | Holm p | separates? |
|---|---|---|---|
| as scored | +0.023 | 0.0459 | yes |
| matching alone | +0.040 | 0.0005 | yes |
| route-adjusted | +0.003 | 1.0000 | no |
| intermediate-only | +0.038 | 0.0251 | yes |

So the pair separates as scored and by matching alone, but not route-adjusted. The results sentence should not claim the separation as a difference in reaching the milestones.

- **Why it fails.** Route adjustment raises DeepSeek most: 0.900 to 0.982, against Claude's 0.924 to 0.985.
  - DeepSeek leaves more milestones unmatched by E3: 1,535 against Claude's 1,110, over the same 2,180 responses each.
  - The judge rules a similar share of each model's unmatched milestones not needed: 0.608 for DeepSeek (933) and 0.634 for Claude (704).
  - Taking those out of the denominator closes the gap.
  - Under route adjustment DeepSeek newly separates from Gemma 4, Qwen3-235B and GLM-5.3.

  The counts are `models[key].sources` in the JSON.
- **Pairs separated after Holm:**
  - as scored: 28;
  - matching alone: 24 (10 verdicts differ from as scored);
  - route-adjusted: 29 (7 differ);
  - intermediate-only: 19.
- **The median smallest detectable difference** at the strictest Holm step is 0.046, 0.050, 0.040 and 0.065 for the four readings.

**Verbosity.** Pooled over the 11 default models, with template fixed effects: 23,722 readable responses on 147 templates, 95% intervals resampling templates.

| slope | numeric values shown (per 10) | visible completion tokens (per 1,000) |
|---|---|---|
| pooled | +0.0014 [−0.0002, +0.0033] | +0.0042 [−0.0066, +0.0152] |
| within model and template (template-by-model fixed effects) | +0.0021 [+0.0003, +0.0041] | −0.0080 [−0.0206, +0.0072] |

- Over the pooled interquartile range of the count (35 to 74 values), the within slope comes to about +0.008 MC.
- The within-model Spearman of MC against the count runs from −0.235 to −0.067. It has no template control, so harder templates drawing longer responses push it negative.
- Across the 11 models' means, the Spearman is +0.309.
- 21 GLM rows (GLM-5.3 18, GLM-5.3 Flash 3) report more reasoning than completion tokens and are floored at zero visible tokens.

**Reasoning tokens (descriptive).** In every model that reasons, MC is lower in the highest reasoning-token quartile than in the lowest. It does not always fall monotonically from quartile to quartile:

- Q1 is 0.930 to 0.956, and Q4 is 0.632 to 0.907.
- The within-model Spearman runs from −0.401 to −0.169.
- Q4 holds the empty responses; GLM-5.3 is 18.2% empty in Q4.
- Gemma 4, GPT-5.4 mini and Gemini 3.1 Flash-Lite are read from `reasoning-medium-full`. Qwen3-235B has no reasoning tokens in any store.

## For other streams

1. **C1, C4, D1: `results/results.json` is stale for MC.**
   - On the current main store (after the 30-item re-run), Table 1's MC pair Claude Sonnet 5 vs DeepSeek V4.1 Flash has Holm p **0.0459**. `results.json` still has 0.0238 from before the re-run.
   - I confirmed 0.0459 with `analyze.q3_coverage` itself at 100,000 flips. The count of separated MC pairs is 28 in both.
   - The as-scored separation is now marginal. Whatever sentence states it should read the regenerated value.
2. **WS-G and the orchestrator: run C2 without `--quick`.** Quick mode cannot decide this pair: with 10,000 flips the as-scored Holm p is 0.076 against 0.046 at full resolution, and the pair sits at rank 28 of 55. Re-run C2 after the one re-score as well. If `answer.py` changes the stored `targets` (the symbolic check), the intermediate-only reading moves with them.
3. **C1: keep the imported names and signatures listed above.** If `B` stops being a module global read at call time, quick mode must change; tell C2.
4. **C4: schema.** As the brief gives it, plus these extra keys:
   - per reading: `wrong_ci`, `wrong_n`, `responses`;
   - per model: `store`, `route_adjusted_left_out` and `sources` (the judge's source counts, `unmatched` and `not_needed_share_of_unmatched`), and `matched_store.store`;
   - `pairs.as_scored`, and each pair's `p` and `detectable_holm`;
   - `separation_rule.claude_vs_deepseek.diff` (per reading);
   - top level: `draws`, `templates`, `keys` and `provenance`;
   - verbosity: `within_model`, `spearman_within_tokens`, `spearman_between_models`, `per_model`, `iqr_*` and `tokens_floored`;
   - per reasoning model: `store`, `tokens_median`, `empty_share`, `responses` and `spearman`.

   Models are keyed by store key; `reasoning_tokens[key]` is null for Qwen3-235B.

## Open items

- None blocking. The full-resolution JSON and section on disk are from the final script. `provenance.script_sha256_lf` matches the file.
