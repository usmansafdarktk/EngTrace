# WS-C1: Matched-settings family and exports (Phase 1, parallel split of WS-C)

Read first: `WS-C_analysis_code.md` (context, rules, the amendments of 2026-10-07) and `00_ORCHESTRATION.md`
section 7 (the registry, including the output-file schema added for the split). This brief adds only what is
specific to C1. Target: three to four hours.

## Mission

Make `analyze.py` produce the matched-settings family and the exports the paper needs, each into its own output
file so that four other sessions can work at the same time. Closes the analysis side of review 1 M4, M9a, M9d
and review 2 W6, W12, W13, Q7.

## Files you own

`full_run_28092026/analyze.py`; new outputs `full_run_28092026/results/matched.json`, `results/single_path.json`,
`results/providers.json`, `results/sections/matched.md`, `results/sections/single_path.md`,
`results/sections/providers.md`; `results/judge_swap_main.json` (re-run `judge_swap.py` unchanged). Nothing
else. Do not write `results/results.json` or `results/RESULTS.md` except through `analyze.py`'s existing path,
and do not run `paper_results.py`.

## Quick mode

Add `--quick` to `analyze.py` (or honour the environment variable `ENGTRACE_QUICK=1`): 1,000 bootstrap
resamples and 10,000 permutations instead of 10,000 and 100,000. Develop and test in quick mode; print "QUICK"
in every output written under it. The orchestrator runs the full-resolution pass at integration.

## Steps

1. **Matched family.** From `results/matched_config.json`: skip any entry with `run: false` or no store; the
   configuration "each model at its reasoning store where one exists, else its default store". Compute per
   model: FAC with template interval, the letter display, the share with no readable answer, level means and the
   Easy-minus-Advanced gap with its Holm-adjusted Welch p, MC with and without the judge with intervals, the
   arithmetic and judged-step flag rates on correct answers; the 55 pairwise sign-flip tests with Holm and the
   detectable difference; the MC pairwise family; Kendall's tau between the default and matched orderings. Per
   re-run model, the paired change against its default rows over all 2,250 items (change, 95% interval,
   Holm-adjusted sign-flip p over templates, detectable change, 90% interval against ±0.05, MC change, empty
   responses in each configuration). Write `results/matched.json` and `results/sections/matched.md`.
2. **Single path.** Export Q4 with intervals, including the "none" column and the 92-others columns, to
   `results/single_path.json` and `results/sections/single_path.md`.
3. **Providers.** Export "By serving endpoint" with the matched difference, templates matched, and a flag for
   endpoints that served fewer than 20 templates, to `results/providers.json` and `results/sections/providers.md`.
4. **E3 pairwise family** (coverage by matching alone) as Q3 does for E5-strict, into `results/matched.json`
   under `e3_pairs_default`, so C2 does not have to duplicate the pairwise machinery; expose the pairwise
   functions (sign-flip, Holm, detectable, letters) as importable functions with a stable signature, and say so
   in your report: C2 and C3 import them.
5. **Judge swap.** Re-run `judge_swap.py`; confirm 211 responses remain; note the count in your report.
6. **Twelfth model.** Model lists come from the stores present plus `matched_config.json`, never a hard-coded
   eleven; the Qwen Thinking entry is skipped.

## Output schema (shared with C4)

`results/matched.json`:
```
{"quick": bool, "config": {...}, "models": {key: {"store", "fac", "fac_ci", "letter", "unreadable", "mc_strict",
 "mc_strict_ci", "mc_e3", "mc_e3_ci", "levels": {"Easy","Intermediate","Advanced"}, "gap", "gap_ci", "gap_p_holm",
 "digit_flag_rate", "router_flag_rate"}},
 "pairs": [{"a","b","diff","ci","p_holm","detectable"}], "mc_pairs": [...], "e3_pairs_default": [...],
 "paired_change": {key: {"fac_default","fac_reasoning","change","ci","p_holm","detectable","ci90","mc_change",
 "empty_default","empty_reasoning"}}, "tau_default_vs_matched": {"tau","ci"}}
```
`results/single_path.json`: `{key: {"single": {"all","some","none", "ci": {...}}, "others": {...}}}` plus
`"n_single": 58, "n_others": 92`.
`results/providers.json`: `[{"model","endpoint","rows","raw_score","unusable","matched_diff","templates_matched",
"few_templates"}]`.

## Acceptance

`analyze.py --quick` runs end to end on the current stores and writes the three files; the existing `RESULTS.md`
and `results.json` outputs are unchanged in content apart from the swap count; a self-test (`--selftest` or a
small assert block) checks the schema. Report: the matched-family tiers in quick mode, the paired changes, the
function signatures C2 and C3 import, anything another stream must do.
