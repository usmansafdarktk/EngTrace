# SWITCH W4: aligning prose and generated phrases (matched headline)

## Phrases changed in `paper_results.py`
- **Anchors (first five on the subset):** GPT-5.4 mini's subset-arm score was replaced by its `reasoning-medium-full` responses on the 450 `subsamples.paraphrase_ids()` items. The score is the mean of template means, as `analyze.anchors` computes it, and it has no interval. Two drift checks were added. Result: 0.964 → **0.967** to 0.990; Advanced 0.94 to 0.99 is unchanged. "Inside the interval of" is computed over the four models that have intervals (DeepSeek V4.1 Flash).
- **Paraphrase τ:** "below what sampling noise alone gives" → "within the range sampling noise alone gives, below its lower quartile". "; 5th percentile 0.647" was added from `q5.tau.noise_arm.p5`.
- **MC, main text:** the misplaced p value was moved. Old: "nor does the final answer (Holm-adjusted $p$ 0.069)". New: "after Holm's correction (Holm-adjusted $p$ 0.069), nor does the final answer". The 0.069 is MC's own p value.
- **Level gap:** the clause "and for four models under the permutation of level labels" is now generated.
- **MC pairs, appendix:** "each of them also on FAC" is added when SEP_MC ⊆ SEP_FAC, which holds.

## Phrases added (matched only)
- "GPT-5.4 mini, Gemma 4 26B, and Gemini 3.1 Flash-Lite return no reasoning tokens and score 0.852, 0.872, and 0.878" (main text and results.tex). A claim checks the no-reasoning-tokens part.
- "At these defaults, … an upper six each of which differs from each of a lower five, and GPT-5.4 mini ranks tenth"

## Numbers removed from the prose
- Main text: 0.032; 51 and 53; 0.098, 0.767 and 0.112 (GLM-5.3 is now one of the last five).
- results.tex: +0.038/+0.070/+0.118; "other eight".
- further_experiments: 0.959.
- error_analysis: 28, 13, 8, 27, 3, 311, 120. These are replaced by 34, 18, 9, 38, 5, 240, 108.
- branch_domain: "three highest".

The six long lines in results.tex were wrapped, and so was one in paraphrase.tex.

## Checks (final)
| check | result |
|---|---|
| `paper_results.py --check --headline matched --repaired` | 0 failures |
| `paper_setup.py --check` | pass |
| `appendix_evaluation.py --check` | pass |
| `appendix_statistics.py --check` | pass |
| `appendix_certification.py --check` | pass |
| `appendix_listings.py --check` | pass |

## For the owner's eye
- **results.tex, pairwise:** the generated list replaces the prose agent's grouped list. The generated list is longer and says only "27 of the 30" for the pairs across groups. The phrase could be improved.
- **Details dropped because the phrases don't hold them:**
  - the named tolerance-swap pairs (results.tex);
  - "GPT-5.4 mini is lowest in industrial" (branch_domain);
  - the main text's "ahead of Muse… 99 empty… seventh" (repeated in the readable-only sentence).
- **Main text:** GPT-5.4 mini "no longer differs from Kimi K3 and GLM-5.3-Flash" follows the phrase, which looks only at the first five. The prose agent also listed Muse Glimmer 30B and GLM-5.3.
- **Number words:** the check tests them only as a set, per file. Words it never flagged, such as "nine of the eleven" in results.tex, were not verified sentence by sentence.
