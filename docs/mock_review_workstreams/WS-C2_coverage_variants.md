# WS-C2: Milestone Coverage variants (Phase 1, parallel split of WS-C)

Read first: `WS-C_analysis_code.md` (context, rules, amendments; its step C2 has the method) and
`00_ORCHESTRATION.md` section 7. This brief adds only what is specific to C2. Target: three to four hours.

## Mission

Compute the Milestone Coverage variants and their tests from the stored scores, so the paper can say what MC
measures and whether the one top-tier separation survives its own controls. Closes the analysis side of
review 1 M1 (a to f) and review 2 W3, Q3, Q4.

## Files you own

New: `full_run_28092026/coverage_variants.py`, `results/coverage_variants.json`, `results/sections/coverage_variants.md`.
Nothing else. Import the pairwise functions from `analyze.py` (C1 exposes them; until then, copy the sign-flip,
Holm and detectable-difference logic into a private helper and replace it with the import when C1 reports the
signatures). Do not edit `analyze.py`, `answer.py` or any store.

## Quick mode

`--quick`: 1,000 resamples, 10,000 permutations; print "QUICK" in the outputs.

## Steps

1. Per response, from `scores/<store>/<model>.jsonl` (E3 `reached` per milestone) and `scores/<store>/e5/<model>.jsonl`
   (per-milestone `sources`): MC as scored (assert it equals `e5_strict`), MC by matching alone, route-adjusted MC
   (milestones ruled `NOT_NEEDED` leave the denominator), intermediate-only MC (the answer-target milestones leave
   the set). For the target rule read `evaluation/milestones.py` first; if it records which milestones the gold
   answer line states, use that; otherwise match milestone values to the answer row's `targets.numbers` under the
   check's unit factors. State the rule and the count of instances left without milestones.
2. Per model: means with template bootstrap for all responses and for wrong answers, for the default stores and,
   where `results/matched_config.json` names one, the reasoning store.
3. The 55-pair tests (sign-flip, Holm, detectable) for matching alone, route-adjusted and intermediate-only, on
   the default configuration.
4. The cross-model verbosity model: response-level regression of MC on the number of numeric values the response
   displays (count with the arithmetic module's reader) and on completion tokens, with template fixed effects,
   pooled over models; slope with a template-cluster bootstrap interval; the within-model Spearman beside it.
5. MC against reasoning tokens per model, descriptive (means by within-model quartile of `reasoning_tokens`).
6. The decision rule output, stated explicitly in the section file: whether Claude Sonnet 5 and DeepSeek V4.1
   Flash separate after Holm under (a) matching alone and (b) route-adjusted MC.

## Output schema (shared with C4)

```
{"quick": bool, "rule": {"targets": "...", "instances_without_milestones": n},
 "models": {key: {"as_scored": {"all","ci","wrong"}, "matching_only": {...}, "route_adjusted": {...},
                  "intermediate_only": {...}, "matched_store": {... same four, or null}}},
 "pairs": {"matching_only": [{"a","b","diff","ci","p_holm","detectable"}], "route_adjusted": [...], "intermediate_only": [...]},
 "separation_rule": {"claude_vs_deepseek": {"as_scored": p, "matching_only": p, "route_adjusted": p, "holds": bool}},
 "verbosity": {"slope_numbers","ci","slope_tokens","ci_tokens","spearman_within": {key: rho}},
 "reasoning_tokens": {key: {"q1","q2","q3","q4"}}}
```

## Acceptance

`coverage_variants.py --quick` runs on the current stores and writes both files; `--selftest` reproduces
`e5_strict` for every response and the variants on a hand-made example. Report: the per-model variant means, the
separation-rule outcome, the verbosity slope, the target rule used.
