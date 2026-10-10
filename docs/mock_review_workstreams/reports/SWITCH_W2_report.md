# SWITCH W2: paper_results.py at the matched headline

Only `full_run_28092026/paper_results.py` was edited. Nothing was committed.

## Blocks under `--headline matched`
- **tab:main_results:** the rows are unchanged. In the caption, the judge share is now 11–24% (was 25%) and ∗ reads "No reasoning setting offered".
- **tab:matched:** the default columns always come from results.json. Rows follow the matched order.
- **tab:judged_steps:** full rows for all eleven models, and the bottom block is dropped. GPT-5.4 mini: 0.115, with 0.770 on wrong answers.
- **tab:single_path:** from `matched.json['q4']`. GPT-5.4 mini's single-path "some" share is 0.345.
- **tab:coverage_variants:** Holm p is 0.069 as scored, 0.001 by matching alone, 1.000 route-adjusted and 0.032 intermediate only. The slope is +0.0020 (0.0005 to 0.0041); the interval now prints at 4 decimals.
- **tab:depth_model:** the slope holds for 3 models (was 5). The pooled odds ratio is 1.21 (1.09–1.33).
- **tab:flag_precision:**
  - The three re-run models now read `reasoning-medium-full` (GPT-5.4 mini 384→187, Gemma 487→347).
  - All models: 3,077 flags, 12 carried.
  - The caption now says the expert's slips were read on responses at the providers' default settings.
- **tab:branch_domain and tab:level_gap:** GPT-5.4 mini's gap is +0.006 (was +0.160). The Welch test holds for gpt-oss-20b alone.
- **tab:coverage:** "--" in the flag-precision column for the three re-run models, with a caption clause.
- **tab:paraphrase:** from `matched.json['q5']`.
- **tab:experiments:** the rows are unchanged and the caption has one added sentence.
- **tab:errors and tab:providers:** unchanged. models.tex is byte-identical.

## Phrases whose claim changed
- **Tiers:** two tiers → "a graded order of eight overlapping groups, FAC 0.808–0.988; 36 of 55 pairs differ". Only gpt-oss-20b differs from every other model.
- **First five:** "only three of ten pairs differ", detectable at one to two points. The GLM-5.3 clause is dropped.
- **Solved on all instances:** first five 72–90%, other six 31–81%.
- **"Four of the lower five…"** → "only GPT-5.4 mini of the three re-run models reaches the first five".
- **Claude vs DeepSeek MC:** separated as scored → not separated (0.069). Matching alone still separates them (0.001).
- **Readable-only:** reorders only the first seven; GLM-5.3 moves from seventh to first; τ 0.636.
- **"Two tiers give way…":** now "rises to 0.969, level with Muse, no longer differs from Kimi K3 and GLM-5.3-Flash".
- **Wording:** "lower tier" → "last four" (wrong answers) and "last five" (branch, depth). Muse is the model "between".
- **Depth:** the slope holds for three of eleven.
- **MC pairs:** "nearly the same (28 vs 26)" → "matching alone separates more (20 vs 15)".
- **Lost points (residual file):** 231 points lost; 126 incorrect, of which 38 are near and 5 symbolic. Exact-digit templates: 240 incorrect verdicts.
- **Experts' sample:** "contrast the tiers" → "span the FAC order".
- **Captions:** GPT-5.4 mini is described as "(closed, fifth on FAC)".

## Figures
- **Copied to `overleaf_source_04102026/figs/`:** branch-bars, level-bars, domain-radar-labeled and -unlabeled, level-gap, coverage-wrong and paraphrase.
- **PNGs:** in `<scratchpad>/w2/png/`.
- **error-categories.pdf:** byte-identical.
- **Axes:** no axis change, and paper_figures.py is untouched.
- **New drawing:** `write_blocks` now also draws the three figures not placed in the paper.

## `--check --out`
- Block failures: 0. Claim and drift failures: 0.
- Prose: 63 items pending (35 phrases, 28 numbers), plus 6 prose lines over 100 characters in results.tex.

## Caveats
- **Anchors sentence:** GPT-5.4 mini's subset score for "the first five" (0.964) comes from its separate reasoning-medium experiment on the same 450 instances. That run has no interval of its own.
- **Not switched:**
  - The judge-swap shift range is still measured on default responses.
  - The expert-reading phrases are still on default responses.
