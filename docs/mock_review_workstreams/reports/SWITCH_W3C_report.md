# SWITCH W3C: numbers sheet at matched settings

Only `numbers_sheet.py` edited; `NUMBERS_SHEET.md` regenerated (518 lines). Nothing committed.

## Sections
- **Changed, matched first:** Setup, Table 1, pairs, first five, stability, single-path, levels, depth, MC readings, flags, scoring variants, error bins, repeats, paraphrases. The re-run models' default rows follow each table.
- **"Defaults by design" notes:** judge swap, further conditions, error analysis, answer-kind audit, endpoints.
- **Added:** default-family tiers; strict-FAC, McNemar, Wilcoxon counts; sixth-place FAC; tau noise floor; permutation and no-chemical gap holders; branch pairs; paraphrase margin count and tau; matched carried precision.

## Headline at matched settings
- **First five:** DeepSeek, Claude, Kimi, GLM-5.3-Flash, GPT-5.4 mini (0.9693; Muse sixth at 0.9691).
- **Pairs separated (of 55):**
  - FAC 36 (strict 29, McNemar 45)
  - MC 15 (Wilcoxon 17)
  - matching alone 20
  - default: 38 / 28 / 26
- **Tau, default vs matched:** 0.709 [0.636, 0.818].
- **Level gap holders (Holm):**
  - Welch: gpt-oss only
  - permutation: GLM-5.3-Flash, GLM-5.3, Gemini, gpt-oss
  - without the chemical pair: none
- **Depth:** the slope holds for 3 of 11 (DeepSeek, Gemini, gpt-oss); 5 at default. Pooled +0.187 [0.088, 0.285].
- **Paraphrase:** 9 of 11 within ±0.05 (Qwen3 and gpt-oss outside). Tau 0.722 against a noise median of 0.807.

## Missing from the result files
- No matched totals for the relative-error bins or `flag_precision`. The sheet sums them per model and labels them "summed here".
- `answer_audit.json` exists only for the main run.
- GPT-5.4 mini has no repeat runs.
