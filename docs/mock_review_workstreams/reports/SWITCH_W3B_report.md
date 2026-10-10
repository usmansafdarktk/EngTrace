# SWITCH W3B: appendix prose at matched settings

Five files edited, prose only. The generated blocks were not touched (checked by script). Nothing was committed.

## results.tex
- Pairwise: 38 → 36. The list of separated pairs is rewritten for the new first five (it ends "gpt-oss-20b from every other model"). Detectable range 0.064 → 0.039.
- MC: "28 differ ... Claude and DeepSeek ... do" → "15 differ, each of them also on FAC ... neither do Claude and DeepSeek (Holm p 0.069)". τ 0.782 (0.564–0.891) → 0.709 (0.345–0.818).
- Strict 34/51 → 29/48. Wilcoxon 29 → 17.
- Matched Settings: matched is now "the configuration reported throughout". The defaults sentence now reads: 0.878, 0.872, 0.852; the default first five include Muse; τ 0.709. Changes +0.038/+0.070/+0.118 added.
- Consistency: "lower-tier vs upper-tier" → "the two lowest, Qwen3 and gpt-oss-20b, ... than any other model". First five against the rest does not hold, because GPT-5.4 mini has 20 such templates and Muse has 13.
- Scoring: halving swaps GPT-5.4 mini/Muse and doubling swaps two (τ 0.964, 0.927). Bound 0.033 → 0.028. Readable-only τ 0.709 → 0.636.
- Coverage:
  - Slope bound 0.003 → 0.004.
  - ρ steps −0.26 → −0.28.
  - Slope +0.0020 (0.001 to 0.004).
  - Within-model ρ −0.24 → −0.21.
  - Claude–DeepSeek: separated by matching alone (0.001) and intermediate-only (0.032), not as scored (0.069) or route-adjusted.
  - Judge adds 0.019 to 0.045. "nearly the same (28 against 26)" → "with it fewer pairs separate (15 against 20)".
  - Carried precision 11/3,409 (0.3%) → 12/3,077 (0.4%).
- Difficulty:
  - Permutation: six → four.
  - Odds ratio 1.23 (1.13–1.34) → 1.21 (1.09–1.33).
  - Depth slope: "five ... four weakest" → "three: DeepSeek, Gemini, gpt-oss-20b; other eight".

## paraphrase.tex
- τ 0.807 (0.547–0.844) "is what noise gives" → 0.722 (0.449–0.844) "within what noise gives". Quartiles 0.748–0.860, and the 5th percentile 0.647 is added.
- Repeat SD 0.003–0.016 → 0.003–0.008.

## further_experiments.tex
- Intro: one sentence added. The experiments and their evaluation column run at the providers' defaults, so GPT-5.4 mini and Gemini are without reasoning there. The pointer to `tab:matched` stays.
- Anchors: "upper tier" → "first five". The comparison now reads "first five at the providers' defaults".
- Open Book: "top-tier/lower-tier" → "one of the first five and two lower-scoring ones at the providers' defaults".

## error_analysis.tex
- "contrast the tiers" → "span the range".
- Incorrect verdicts: 99 → 126. Answer points 231.5/102/30.5/61 → 231/71/34/68.
- Shortcut-template bound 0.004 → 0.002.

## branch_domain.tex
- Rationale: GPT-5.4 mini is now "closed and fifth on FAC".
- Lowest branch: chemical 6 / industrial 4 → 5 / 5.
- Bars: "the three highest of the four stay close (GPT-5.4 mini lowest in industrial)".
- Radar: DeepSeek and GPT-5.4 mini are lowest on production and inventory. "Lower-tier" wording removed.

No heading was changed.

## Prose word counts (outside generated blocks)
| file | before | after |
|---|---:|---:|
| results.tex | 1154 | 1219 |
| paraphrase.tex | 481 | 483 |
| further_experiments.tex | 803 | 862 |
| error_analysis.tex | 794 | 794 |
| branch_domain.tex | 168 | 149 |

## Needs the owner's eye
1. **error_analysis, verdict breakdown:** "28 on three templates (13/8/7)", "27 within 0.2%", "3 symbolic", and "311 ... 120 on 11 templates" are still the default numbers, with Muse where GPT-5.4 mini now belongs.
   - `matched.json['reported']['classification']` holds the accuracy of the classification templates, not this breakdown.
   - `RESIDUAL_INCORRECT.md` covers `main` only, so `residual_incorrect.py` needs a matched run.
2. **Anchors:** 0.959–0.990 and Advanced 0.94–0.99 are the default first five. The `paper_results.py` claim "top five on the subset" will fail under the matched order.
3. **Carried precision:** the count is 12, not the brief's 11. It is summed over stores, because no script prints the matched total.
4. **Matching-alone p:** kept at 0.001. That is how `pv()` prints 0.00055.
5. **Second-judge shift:** the −0.018 to +0.019 figure is from `judge_swap_main.json`, the default store. Unchanged.
