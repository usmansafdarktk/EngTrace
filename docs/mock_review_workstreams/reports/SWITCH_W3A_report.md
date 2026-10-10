# SWITCH W3A: 6_results.tex prose at matched settings

## Sentences changed
- 5.3.1 S1: two tiers, 0.950–0.988 / 0.808–0.899 → graded order 0.808–0.988; 36 of 55 pairs differ; only gpt-oss-20b differs from all
- S2: GLM-5.3 sixth; 1 of 10 pairs; 1–3 points → Muse sixth, GLM-5.3 seventh; 3 of 10 (DeepSeek > GLM-5.3-Flash, GPT-5.4 mini; Claude > GPT-5.4 mini); 1–2 points
- S4: >6 templates; 81–90% / 31–63% → >5; first five 72–90%, two lowest 51%, 31%
- S5: four of the lower five lack reasoning tokens → only Qwen3-235B-2507 does
- Claude–DeepSeek: 0.049 / 0.001 / 1.000 → matching 0.001, intermediate only 0.032, as scored 0.069, route-adjusted 1.000
- Judged step flags: 21.9% → 20.6%
- Wrong answers: five lower-tier models, 44–59% → four lowest, 46–59%; 26% → 25%; missing 63–78% / 50–83% → 63–78% / 46–83%
- Single-path: ≤59% / ≤22% → two lowest 53%, 59%; other nine ≤34%
- Stability S1: four lower-tier models, across the tiers → four lowest models, past one it differs from
- Repeats SD: ≤0.016 → ≤0.008 (models in table order)
- Readable only: upper tier, τ 0.709, sixth → first seven, τ 0.636, seventh
- Closing: defaults, tiers survive → matched settings, the pairs that differ survive
- Reasoning Setting: defaults 0.878, 0.872, 0.852; +0.038, +0.070, +0.118; two tiers at defaults, GPT-5.4 mini tenth; with reasoning it is in the first five and level with Muse, and no longer differs from Kimi, GLM-5.3-Flash, Muse, GLM-5.3
- Branch: 0.018–0.056 → 0.018–0.053; lower tier <0.79, GLM-5.3 between → four lowest 0.071–0.226, <0.81; Muse, GLM-5.3 between (0.056, 0.098; 0.880, 0.767); thermodynamics 5 → 4
- Level gap: 0.008–0.181, 5 include zero → 0.006–0.181, 7; Welch: Gemini and gpt-oss → gpt-oss alone; permutation: four
- Depth: lower tier 0.121–0.271 → four lowest 0.104–0.254; OR 1.23 (1.13–1.34), 5 → 1.21 (1.09–1.33), 3
- Experiments: clause added (defaults, paired with default responses); first five on subset 0.959–0.990 → 0.967–0.990; "upper tier" → "first five"
- Error analysis: "contrast the tiers" → "span the range"; largest share: no readable answer → incorrect verdicts; 231.5 → 231; 102 → 71; 30.5/61 → 34/68; 99/27/3/69 → 126/38/5/83

## Headings
- "Two Tiers and a Ceiling" → "A Graded Order and a Ceiling"; "Stability of the Tiers" → "Stability of the Order". Both changed because there are no tiers at matched settings.

## Words (outside generated blocks)
2,070 → 2,110.

## Needs the owner's eye
- Stability S1: for paraphrase, the claim holds at the point estimates only; at the 90% bounds some pairs that differ could swap.
- No results file holds 0.967 (GPT-5.4 mini's reasoning-on store on the 450-item subset) or 38/5/83 (GPT-5.4 mini at reasoning: 14 near misses, 3 symbolic), so I computed both separately. `reported.classification` holds answer-kind scores, not this breakdown. Both numbers need a committed script.
- "No model of the first five" includes GPT-5.4 mini, which was tested at its default.
- `appendix:error_analysis` and `tab:coverage_variants` may still show default figures.
