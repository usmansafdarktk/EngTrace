# SWITCH W1A: matched.json carries every per-model family

Only `full_run_28092026/analyze.py` was edited. Nothing was committed.

## Keys added to `results/matched.json`
Every existing key is unchanged. All 22 differences from the previous file are additions.
- `models[k].q3`: q3()'s full row at the model's matched store. Same shape as an entry of `results.json['q3']`.
- `q1`: `{models, pairs}`, as in `results.json['q1']`. The pairs carry `fully_p_holm` and `mcnemar_p_holm`.
- `q3_coverage`: `{stage, templates, models, pairs, tau}`. The pairs carry `p_wilcoxon_holm`.
- `pairs_separated`: counts at Holm 0.05, in these fields: `pairs`, `fac`, `fac_strict`, `fac_strict_agree`, `fac_mcnemar`, `mc`, `mc_wilcoxon`, `mc_wilcoxon_agree`, `mc_e3`.
- `q2`, `q3_overall`, `branches_levels`, `q4`: lists with one row per model, as in results.json.
- `q5`: `{models, tau, margin, expert_check, vs_repeats, expert_stats, stores}`. `stores` gives each model's paraphrase store.
- `repeats`: `{model: {items, scores{store: score}, main_on_same_items, sd, range, same_verdict_every_repeat}}`.
- `sensitivity`: `{models, tau_with_headline, tau_noise}`.
- `reported`: `{models, classification, providers}`.

`sections/matched.md` prints each family as one table or list.

## Matched headline (full resolution)
- **FAC pairs separated, of 55:**
  - sign-flip test: 36
  - strict FAC: 29 (agrees with FAC on 48)
  - McNemar: 45
- **MC pairs separated:**
  - sign-flip test: 15
  - Wilcoxon: 17 (the two agree on 51)
  - matching alone: 20
- **Kendall tau, default vs matched FAC ordering:** 0.709 (0.636 to 0.818).
- **Level gap (Easy minus Advanced):**
  - Welch (Holm): holds only for gpt-oss-20b (+0.181).
  - Permutation (Holm): holds for glm-5.3-flash (+0.055), glm-5.3 (+0.120), gemini (+0.132) and gpt-oss (+0.181).
  - Without the two chemical templates: holds for none.
  - Other gaps: deepseek +0.016, claude +0.008, kimi +0.039, gpt-5.4-mini +0.006, muse +0.019, gemma +0.055, qwen +0.031.
- **Paraphrase (275 kept pairs, Holm over 11):**
  - No model's change is significant; the lowest adjusted p is 0.527 (muse).
  - Within ±0.05 for 9 of 11. The two outside are gpt-oss +0.018 (90% CI −0.023 to 0.059) and qwen −0.036 (−0.073 to −0.002).
  - Re-run models: gemma +0.000, gpt-5.4-mini +0.005, gemini −0.013.
  - Others: deepseek −0.013, glm-5.3-flash −0.011, glm-5.3 +0.002, muse +0.025, kimi +0.002, claude −0.011.
  - Paraphrase-vs-original tau: 0.722 (0.449 to 0.844), against an arm-size noise median of 0.807.
- **Repeats SD (300 instances, reasoning-on repeats):** gemma 0.0075, gemini 0.0067. The default entries stay: gpt-oss 0.0025, qwen 0.0051.
- **Sensitivity taus:**
  - 0.964: half tolerance, fully solved, without shortcut, without symbolic, half unit
  - 0.927: double tolerance
  - 0.636: unusable excluded
  - 1.000: whole trace
- **First five on FAC:** deepseek, claude, kimi, glm-5.3-flash and gpt-5.4-mini. gpt-5.4-mini (0.9693) is now fifth, ahead of muse (0.9691).
  - FAC spread across the five: 0.019.
  - Lowest branch and branch spread: industrial for deepseek (0.018), claude (0.026), glm-5.3-flash (0.053) and gpt-5.4-mini (0.023); chemical for kimi (0.036).
  - Lowest domain: production_and_inventory for deepseek (0.940), glm-5.3-flash (0.887) and gpt-5.4-mini (0.933); quality_and_reliability_control for claude (0.940); thermodynamics for kimi (0.917).
  - Answer points lost: 231 in total: 71 with no readable answer, 68 partial answers at half a point, 126 incorrect.
  - One branch pair still holds: gpt-oss-20b, civil below electrical (−0.226).

## Default results unchanged
- `results.json` differs only in `provenance.analyze`: `git`, `dirty` and `sha256_lf`, which record the script's commit, its dirty flag and its own hash. No timestamp changed.
- `RESULTS.md` differs only in its provenance line.
- `per_template.csv`, `single_path.json` and `providers.json` are identical to before.

## Runtime and self-test
- Full pass: 3 min 54 s. Quick pass: 47 s.
- `--selftest` passes. It now feeds stand-in paraphrase and repeat stores and checks two things: a re-run model is paired against its reasoning-on stores, and an unchanged model keeps its default Q2, branch and Q5 numbers.

## None or not computed
- No value is None for the three re-run models. The judge and router stages are complete for them: `reasoning-medium-full` has e5 and router rows at 2,250 each, and `paraphrase-reasoning-medium/e5` has 275 each.
- gpt-5.4-mini has no reasoning-on repeat runs (its cascade lists none), so it has no `repeats` or `vs_repeats` entry.
- qwen3's router stage is still marked incomplete, as it is at the default.
