SIGNAL: FIGURES REFRESHED 2026-10-08 20:20 UTC
SIGNAL: FIGURES EXTRACTED 2026-10-07 14:35 UTC (step F1, done by WS-C4; see `WS-C4_report.md`)

# WS-F report

Owner's amendment of 2026-10-08: the data of every figure in `figs/` is brought to the final results, with no redesign
and no new figure. Drawn by `paper_results.py --out <scratch> --headline default --repaired`. The three unplaced figures
are drawn from the same `figure_data()`.

## Design check

Today's `paper_figures.py`, fed the inputs the committed PDFs were drawn from (`results.json` at `081373e`,
`expert_request/scored.json`), reproduces all 8 committed PDFs byte for byte. For `error-categories`, this holds with
the committed bar order. The design is therefore unchanged, and only the data and the three approved fixes below differ.
Every figure's data changed, so none stays byte-identical.

## Fixes in `paper_figures.py` (approved by the owner, 2026-10-08)

- `fig_error_categories`: the bars keep the committed order (Claude Sonnet 5, GPT-5.4 mini, Gemma 4 26B, GPT OSS 20B),
  which is also the order of `tab:errors`. Since `b8b3358`, the pipeline passes the models in table order.
- `fig_level_gap`: the x-axis starts at −0.04, not −0.03. Qwen3-235B-2507's interval starts at −0.035.
- `fig_level_gap`: the "without chemical pair" crosses and their legend entry are removed (the two chemical templates are
  repaired; owner, 2026-10-09).
- `fig_paraphrase`: the "No change" label moves to the top gap, because it sat on Kimi K3's marker. The x-axis runs from
  −0.075 to 0.061, so that the intervals of Qwen3-235B-2507 (−0.073) and gpt-oss-20b (0.059) fit.

## Figures

| label | file | design identical | what changed in the data | approved |
|---|---|---|---|---|
| `fig:level_bars` | `level-bars.pdf` | yes | Advanced: DeepSeek V4.1 Flash 0.93→0.98, Claude Sonnet 5 0.92→0.98, GPT-5.4 mini 0.73→0.76; GPT OSS 20B Easy 0.91→0.90, Intermediate 0.78→0.77 | 2026-10-08 |
| `fig:error_categories` | `error-categories.pdf` | yes (bar order kept by the fix) | Claude Sonnet 5: 120→57 readings, no error 64%→11%, calculation 26%→63%, setup (16) and formula (11) now printed; the other three models shift by at most 4.4 points | 2026-10-08 |
| `fig:branch_bars` | `branch-bars.pdf` | yes | DeepSeek chemical 0.96→0.98, electrical 0.95→0.99, FAC 0.98→0.99; Claude chemical 0.95→0.99, electrical 0.95→0.97, FAC 0.97→0.99; the lower two models move by at most 0.015 | 2026-10-08 |
| `fig:domain_radar` | `domain_radar/domain-radar-labeled.pdf`, `-unlabeled.pdf` | yes | DeepSeek and Claude rise on thermodynamics (+0.06, +0.11) and digital communications (+0.07 each); GPT OSS 20B moves on 11 domains by at most 0.04; every value stays above 0.4 | 2026-10-08 |
| not placed | `level-gap.pdf` | axis from −0.04; chemical-pair crosses removed | every gap shrinks; Gemma 4 26B and GPT-5.4 mini no longer hold after Holm | 2026-10-08 |
| not placed | `coverage-wrong.pdf` | yes | n for the top three 26→12, 50→20, 40→20; with the judge, Claude 0.753→0.869, Kimi K3 0.684→0.584, GLM-5.3 0.845→0.596 | 2026-10-08 |
| not placed | `paraphrase.pdf` | label and axis fixed | gpt-oss-20b no longer within the margin (its 90% interval reaches 0.059) | 2026-10-08 |

In the three unplaced figures, Claude Sonnet 5 and Kimi K3 swap rows (FAC order).

## Captions and prose (for WS-D2)

No caption sentence fails to fit. Claude Sonnet 5 is first on MC and DeepSeek V4.1 Flash first on FAC; the B2 range is
57 to 120 readings; the gpt-oss-20b electrical-above-civil pair holds; every radar value is above 0.4. The prose that
cites the four figures (`6_results.tex`, `appendices/branch_domain.tex`) also fits.
