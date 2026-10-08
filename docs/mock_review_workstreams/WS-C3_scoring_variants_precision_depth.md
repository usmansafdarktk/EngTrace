# WS-C3: Scoring-rule variants, carried precision, depth (Phase 1, parallel split of WS-C)

Read first: `WS-C_analysis_code.md` (context, rules, amendments; its steps C3, C4 and C5 have the methods) and
`00_ORCHESTRATION.md` section 7. This brief adds only what is specific to C3. Target: three to four hours. If
two sessions are available, split it: C3a (step 1) and C3b (steps 2 and 3), which own different files.

## Mission

Three offline analyses from the stored scores: how many verdicts depend on each clause of the final-answer rule
and on prescribed digits, with proportional partial credit as a variant; which arithmetic flags are carried
precision; and a controlled depth model. Closes the analysis side of review 1 M3, M6 (depth), Q4 and review 2
W8, W10, Q5, Q13.

## Files you own

New: `full_run_28092026/clause_variants.py` with `results/sensitivity_variants.json` and
`results/sections/sensitivity_variants.md`; `full_run_28092026/flag_precision.py` with
`results/flag_precision.json` and `results/sections/flag_precision.md`; `full_run_28092026/depth_model.py` with
`results/depth_model.json` and `results/sections/depth_model.md`. Nothing else. `answer.py` and `arith.py` are
imported, never edited; if a clause cannot be switched through their parameters, implement the variant on the
stored `targets` and the response text with a local copy of the rule and say so in the docstring.

## Quick mode

`--quick`: 1,000 resamples where intervals are drawn; print "QUICK".

## Steps

1. **Clause variants** (`clause_variants.py`), over the main store and each reasoning store named in
   `results/matched_config.json`: FAC under (a) the absolute-value clause off, (b) the last-digit term bounded by
   min(u(ŷ), 0.01|y|), (c) the exact-digit requirement of the prescribing templates relaxed to the tolerance,
   (d) proportional partial credit (`matched / of`); per model the number of verdicts that change and in which
   direction; the distribution of relative error among accepted answers (≤0.2%, 0.2 to 1%, 1 to 5%, >5%); Kendall's
   tau of each variant's ordering with the headline.
2. **Carried precision** (`flag_precision.py`): classify every arithmetic flag on a correct-answer response as
   carried precision (the printed result is a correct rounding of the value recomputed from the unrounded upstream
   values that appear earlier in the response) or other, per model; on the flags the domain expert read
   (`FLAG_REVIEW_3.md`; the item list is local), how many of the 171 confirmed slips are carried precision.
3. **Depth model** (`depth_model.py`): wrong-answer rate against the gold trace's milestone count over readable
   responses only, per model: per-bin means with template intervals and a template-clustered logistic model of
   wrong answer on milestone count with answer kind as a covariate; the slope per model, its interval, and whether
   it holds after Holm over models.

## Output schemas (shared with C4)

`results/sensitivity_variants.json`:
```
{"quick": bool, "stores": {store: {"models": {key: {"headline","abs_clause_off","last_digit_bounded",
 "prescribed_relaxed","per_part_credit", "changed": {"abs_clause_off": {"up","down"}, ...}}},
 "tau_with_headline": {...}}}, "relative_error_bins": {key: {"le_0.2","0.2_1","1_5","gt_5"}}}
```
`results/flag_precision.json`: `{"models": {key: {"flags","carried_precision","other"}}, "expert_confirmed": {"slips": 171, "carried_precision": n, "other": n}}`
`results/depth_model.json`: `{"quick": bool, "models": {key: {"bins": {"1": {"rate","ci","n"}, "2": ..., "3": ..., "4-5": ..., "6+": ...}, "slope","ci","p","p_holm"}}}`

## Acceptance

Each script runs in quick mode on the current stores and writes its two files; each has a `--selftest` on a
hand-made example (a clause flip, a carried-precision chain, a depth table). Report: the clause counts per model,
the carried-precision split overall and on the 171 confirmed slips, the depth slopes and their Holm outcome.
