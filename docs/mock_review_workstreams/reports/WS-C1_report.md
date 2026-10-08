SIGNAL: C1 PAIRWISE API 2026-10-07 15:05 UTC (C2 and C3 may import the functions below)
SIGNAL: C1 EXPORTS WRITTEN (quick) 2026-10-07 15:05 UTC (results/matched.json, single_path.json, providers.json and their sections)

# WS-C1 report: matched-settings family and exports

Done in quick mode. Every number below is from the quick pass (1,000 resamples, 10,000 permutations) on the stores as
they are on 2026-10-07: main re-scored after round 5, and WS-B's `reasoning-medium-full` with its judge and step-check
rows. The orchestrator's full pass replaces them. No commits; no paid calls; no `.tex` touched.

## Scripts and flags

| script | what it computes | output | check |
|---|---|---|---|
| `analyze.py` | everything it computed before, unchanged | `results.json`, `RESULTS.md`, `per_template.csv` (full pass, or `--quick --out DIR`) | `--selftest`; HEAD comparison below |
| `analyze.py` | matched-settings family, paired change per re-run model, tau between orderings, E3 pairs on defaults | `results/matched.json`, `results/sections/matched.md` | schema check before every write; `--selftest` |
| `analyze.py` | Q4 single-path table with intervals | `results/single_path.json`, `results/sections/single_path.md` | the same |
| `analyze.py` | endpoint table, templates served, two few-template flags | `results/providers.json`, `results/sections/providers.md` | the same |
| `judge_swap.py` (run unchanged) | judge swap on the remaining sample | `results/judge_swap_main.json`; also rewrites `JUDGE_SWAP.md` | none |

Commands:
- Development: `python -m full_run_28092026.analyze --quick`, or set `ENGTRACE_QUICK=1`.
  - Writes only the three files and their sections.
  - Never overwrites `results.json`, `RESULTS.md` or `per_template.csv` in `results/`, because C4's `--check` reads them.
  - `--out DIR` writes every output, the old three included, under DIR.
- Integration (orchestrator): `python -m full_run_28092026.analyze`, without `--quick`.
  - Writes all six outputs and the three sections.
  - Not timed; the quick pass takes about 4 minutes.
- Self-test: `python -m full_run_28092026.analyze --selftest` passes. It takes 5 to 6.5 minutes at full resolution, mostly the existing tests' 100,000-permutation draws.

What else changed in `analyze.py`:
- **Draw counts.** `B` and `B_TEST` are read when a test runs. `sign_flip_p` and `tier_perm_p` take `draws=None`, so quick mode reaches them.
- **Roster.** `ROSTER` is the eleven, plus any model `matched_config.json` names with rows in the main store (a twelfth row). An entry with `run: false` or no store is skipped, as the Qwen Thinking sibling is.
- **Store per model.** `q3`, `q3_overall` and `q3_coverage` accept a mapping from model to store. `q3_coverage(judge=False)` gives coverage by matching alone.
- **Arm block.** The C1/C4 arm block of `RESULTS.md` leaves out the full-set reasoning stores `matched_config.json` names. HEAD would have added three `reasoning-medium-full` rows there; the matched family reports that change instead, over all 2,250 items.
- **Output order.** `main()` writes every file before it prints.
  - HEAD's `print(text)` raises `UnicodeEncodeError` on a redirected Windows stdout (cp1252 has no "−"); stdout now uses `errors='replace'`.
  - In HEAD that crash came after every file was written. With the new exports it would have come first, and they would have been lost.

## Function signatures C2 and C3 import (stable)

```python
from full_run_28092026 import analyze
analyze.set_quick(on=True)                       # B = 1,000, B_TEST = 10,000 (or run under ENGTRACE_QUICK=1)
analyze.sign_flip_p(d, seed, draws=None)         # draws defaults to analyze.B_TEST at call time
analyze.holm(ps) -> list[float]
analyze.detectable_paired(d, m=1) -> {'detectable', 'detectable_holm', 'sd_paired'}
analyze.boot_mean(v, seed) -> [lo, hi]           # reads analyze.B at call time
analyze.pairwise_family(M, keys, seed_ci, seed_p, wilcoxon=False) -> list[dict]
    # M: models x templates (the caller keeps the templates every model has a value on); per pair a, b, diff, ci
    # (boot_mean, seed_ci + n), p (sign flips, seed_p + n), detectable, detectable_holm, sd_paired, p_holm, and
    # p_wilcoxon, p_wilcoxon_holm when asked. q3_coverage is this with (6000, 7000, wilcoxon=True).
analyze.separated(pairs, key='p_holm', alpha=0.05) -> set[frozenset]
analyze.letters(order, sep) -> {model: letters}  # paper_results.py's insert-and-absorb display, verbatim
analyze.q3_coverage(runs, templates, keys, store='main' | {model: store}, judge=True)
```

C2's current `coverage_variants.py` already calls `load`, `load_stage`, `cluster_mean`, `per_template`, `boot_mean`,
`sign_flip_p(..., draws=)`, `holm` and `detectable_paired`. All keep their signatures, and setting `analyze.B` directly
still works.

## Checks

- **HEAD comparison.**
  - Setup: HEAD's `analyze.py` (a saved copy, forced to the quick draw counts) and the new one (`--quick --out`) on the same stores.
  - `per_template.csv` is byte-identical.
  - These `results.json` keys are identical value for value: `q1`, `q2`, `q3`, `q3_overall`, `q3_coverage`, `branches_levels`, `q4`, `q5`, `sensitivity`, `reported`, `anchors`, `repeats`, `set_aside`, `milestones` and `symbolic_templates`.
  - `reasoning_arms` is identical once HEAD's three `reasoning-medium-full` rows are removed.
  - Provenance: store and stages identical.
  - `RESULTS.md` differs only by the QUICK banner and those three arm rows.
- **Same numbers for unchanged models.** For the 8 models whose rows did not change, every matched-family number equals its default-configuration number: FAC and interval, MC and interval, flag rates and intervals, gap, gap interval, Welch p, level means and intervals, no-readable-answer count, calculations per response, judge-decided share. The paired change's default FAC equals Q1's for the 3 re-run models.
- **Old outputs untouched.** The real quick pass left `results/results.json`, `RESULTS.md` and `per_template.csv` byte-identical (sha256 checked before and after).
- **Schemas.** Each of the three files passes its schema check before it is written. The self-test also checks the checker catches a missing letter and a broken interval.

## Findings on the current stores (quick mode; source `results/matched.json` unless named)

Matched configuration, from `matched_config.json` (sha256 `e1e4d06b307e…`):
- Reasoning at effort medium (`scores/reasoning-medium-full`): GPT-5.4 mini, Gemma 4 26B, Gemini 3.1 Flash-Lite.
- The other eight at `scores/main`.
- Skipped: `qwen3-235b-a22b-thinking-2507` (`run: false`). No twelfth model.

**Tiers on FAC** (letters; models sharing a letter are not separated at Holm 0.05):

| model | FAC, matched | letter, matched | FAC, default | letter, default |
|---|---:|---|---:|---|
| DeepSeek V4.1 Flash | 0.981 | a | 0.981 | a |
| Claude Sonnet 5 | 0.981 | ab | 0.981 | ab |
| Kimi K3 | 0.974 | abcd | 0.974 | ab |
| GLM-5.3-Flash | 0.972 | ac | 0.972 | a |
| Muse Glimmer 30B | 0.969 | abcd | 0.969 | ab |
| GPT-5.4 mini | **0.960** (reasoning on) | cde | 0.850 | cd |
| GLM-5.3 | 0.948 | bdef | 0.948 | b |
| Gemma 4 26B | **0.932** (reasoning on) | efg | 0.868 | cd |
| Gemini 3.1 Flash-Lite | **0.910** (reasoning on) | fg | 0.874 | cd |
| Qwen3-235B-2507 | 0.894 | g | 0.894 | c |
| gpt-oss-20b | 0.814 | h | 0.814 | d |

- **The top five stay the top five.** They are the same five models under both configurations.
- **GPT-5.4 mini.** With reasoning on it is sixth. It is not separated from Kimi K3, GLM-5.3-Flash or Muse Glimmer 30B, and it is separated from DeepSeek V4.1 Flash and Claude Sonnet 5.
- **Pairs separated** (of 55, Holm 0.05; matched against default):

  | family | matched | default |
  |---|---:|---:|
  | FAC | 32 | 33 |
  | MC with the judge | 16 | 27 |
  | MC by matching alone | 16 | 24 |

  Only the 27 pairs that involve a re-run model can change. The re-run models moved into the middle of the MC range.
- **Kendall's tau** between the default and matched FAC orderings: 0.745 (95% CI 0.673 to 0.855, templates resampled).
- **Easy-minus-Advanced gap.** Under matched settings it holds after Holm for gpt-oss-20b and Gemini 3.1 Flash-Lite.

**Paired change per re-run model** (reasoning on minus default, all 2,250 items; Holm over the three):

| model | FAC default → on | change (95% CI) | p (Holm) | detectable | 90% CI | within ±0.05 | MC change with judge (2,180 items) | MC change, matching alone | empty, default / on |
|---|---|---|---:|---:|---|---|---|---|---|
| GPT-5.4 mini | 0.850 → 0.960 | +0.111 (0.078 to 0.149) | 0.0003 | 0.053 | 0.082 to 0.142 | no | +0.057 (0.035 to 0.082), Holm p 0.0003 | +0.049 (0.029 to 0.073) | 0 / 0 |
| Gemma 4 26B | 0.868 → 0.932 | +0.064 (0.042 to 0.088) | 0.0003 | 0.032 | 0.046 to 0.084 | no | +0.008 (−0.003 to 0.020), Holm p 0.177 | −0.003 (−0.019 to 0.010) | 0 / 3 |
| Gemini 3.1 Flash-Lite | 0.874 → 0.910 | +0.036 (0.020 to 0.053) | 0.0003 | 0.023 | 0.023 to 0.051 | no | +0.016 (0.006 to 0.025), Holm p 0.0014 | +0.016 (0.006 to 0.026) | 0 / 0 |

0.0003 is the quick-mode floor (3 × 1/10,001).

**Claude Sonnet 5 against DeepSeek V4.1 Flash on MC, default configuration.** This is C2's separation rule; the numbers are here because the store update moved them.
- With the judge (E5-strict): the difference is −0.023 (DeepSeek minus Claude), raw sign-flip p 0.0027, Holm rank 28 of 55, Holm p 0.076. That is **not separated** at 0.05 (Wilcoxon Holm 0.026).
  - In the committed `results.json` of 4 October (store before round 5, full draws) it was −0.0255, Holm p 0.024, separated.
  - The full pass decides. At rank 28, a raw p of 0.0027 would need to fall below 0.0018.
- By matching alone (E3): −0.040 (−0.059 to −0.023), Holm p at the quick floor (0.0055): **separated** (`e3_pairs_default`).

**Single-path table** (`results/single_path.json`, 58 and 92 templates): the "none" column ranges from 0.000 to 0.069
across models and both groups. The highest is gpt-oss-20b on the single-path templates (4 of 58).

**Endpoints** (`results/providers.json`):
- 40 endpoint rows over 8 models. Claude Sonnet 5, Muse Glimmer 30B and Gemini 3.1 Flash-Lite each had one endpoint in the main run, so they are not listed.
- 16 endpoints served fewer than 20 templates (`few_templates`, as the brief defines it).
- 19 have a matched difference resting on fewer than 20 templates (`few_matched`). Example: Azure served 150 templates for GPT-5.4 mini, but its +0.750 rests on 1 template against one OpenAI row.
- Over the 21 endpoints matched on at least 20 templates, the matched difference lies between −0.028 and +0.018. The largest is GLM-5.3-Flash on StreamLake, −0.028 over 35 templates. So "within ±0.02" (plan WP4.5) is not quite true; the sentence should say ±0.03, or name the exception.
- Reasoning-on stores (`matched.json` → `providers_reasoning`):
  - Gemma 4: Io Net / DekaLLM, ∓0.010 over 36 templates.
  - Gemini: Google AI Studio / Google, ±0.001 over 18 templates.

**Judge swap:** **211 responses remain** (18 to 20 per model, 13 to 19 templates per model), 424 milestones judged by
both (466 before). Pooled agreement 0.844 three-way (0.815 before); kappa on REACHED against not 0.778 (0.790 before).

## Notes for C4 (schema as written)

- **`matched.json`.** It has the brief's keys, plus extra ones. Per model:
  - `configuration` (`reasoning`/`default`), `setting`, `items`, `fully_solved`, `fully_ci`;
  - `all15_share`: Table 1's "All 15 instances solved" (templates solved on all 15 / 150);
  - `no_readable_answer_n`, `empty_n`, `answered_unreadable_n`;
  - `mc_stage`, `mc_letter`, `mc_e3_letter`, `judge_decided_share`;
  - `levels_ci`, `gap_p`, `detectable_gap_holm`;
  - `digit_flag_ci`, `claims_per_trace` ("calculations parsed per response"), `router_flag_ci`, `router_any_flag_rate`;
  - `flag_rates` ({field: {rate, ci, traces}} over every step-flag field present: `digit_flags` today, `formula_flags` if WS-H lands);
  - `single_path`.

  At the top level: `stand_in: false`, `order` (FAC descending) and `default`. `default` holds `order`, `letters`, `mc_letters`, `e3_letters`, the significant-pair counts, `fac` and `mc_strict`, so `--headline matched` can print both tiers. Also: `mc_e3_pairs`, `e3_default`, `providers_reasoning` and `provenance`.
- **`unreadable`** is the share with no readable answer: empty plus answered without a readable answer, Q1's "unusable".
- **Pair rows** add `p` and `detectable_holm`, and `p_wilcoxon_holm` on the MC families.
- **`single_path.json`.** Top-level keys are the model keys plus `quick`, `n_single` and `n_others`. Each group also carries `templates`.
- **`providers.json`** is a list, as specified, so `quick` sits on every row. The rows add `templates_served` and `few_matched`; flag the matched-difference column on `few_matched`, since `few_templates` misses the dominant endpoints.
- **Letters.** `letters` in `analyze.py` is a verbatim copy of `paper_results.py`'s. MC letters are drawn in FAC order, as `CLD_MC` is.

## Open items

- **`judge_swap.py` (not this stream's; nobody owns it in Phase 1).** `JUDGE_SWAP.md` now has two stale phrases.
  - The header says "on 220 sampled traces": `out['sample']` counts the `sample` list in `scores/main/e5_grok-4-6/CONFIG.json`, which still has 220 entries. It should count the traces judged: 211.
  - The closing paragraph says "20 traces per model over 14 to 19 templates": a literal in the script. The current figures are 18 to 20 traces over 13 to 19 templates.
  - `docs/appendix_evaluation.py` (line 164) parses that literal with `(\d+) traces per model over (\d+) to (\d+) templates`, and `docs/check_plan_claims.py` reads the swap table. Fix in Phase 2 (WS-G or D2): generate the sentence and the header count from `out['models']`.
- **`docs/appendix_evaluation.py --check`** will see the new `JUDGE_SWAP.md` figures until Phase 2 updates the appendix.
- **WS-G / orchestrator.**
  - The `results.json` and `RESULTS.md` on disk are still from 4 October, before the round-5 re-run. The full pass changes many numbers because the stores changed, not because of this code; the Claude–DeepSeek MC pair is one.
  - Run the full pass once, after the one re-score at ANSWER FINAL.
- **D1 / WS-G.** The paired-change rows supersede the 450-item reasoning arm's figures for the three models; the arm stays in `RESULTS.md`'s arm block.
- **For C2 (cross-check).** On the default stores, `matched.json`'s `e3_pairs_default` and C2's `matching_only` pairs should agree number for number. Both use q3_coverage's estimand and seeds 6000/7000; compare them at integration.

## Files

- Changed: `full_run_28092026/analyze.py`.
- New:
  - `full_run_28092026/results/matched.json`, `single_path.json`, `providers.json`;
  - `results/sections/matched.md`, `single_path.md`, `providers.md`;
  - this report.
- Regenerated by `judge_swap.py`: `results/judge_swap_main.json` and `full_run_28092026/JUDGE_SWAP.md`.
- Untouched: `results/results.json`, `RESULTS.md`, `per_template.csv`.
