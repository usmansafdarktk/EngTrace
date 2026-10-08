SIGNAL: INTEGRATED 2026-10-08 17:20 UTC (generated blocks written into the real tree with `--headline default --repaired --text-only`; WS-F's figures and captions left for later, on the owner's instruction)
SIGNAL: NUMBERS SHEET 2026-10-08 16:35 UTC (G2, G2b, G3 and G5 done on the final stores, the milestone matcher's floor fixed)

# WS-G report

Status 2026-10-08 17:20 UTC. Committed and pushed: `e83084c` (this stream's scripts and records), `412656a` (the paper-text edits pending in the tree, committed first on the owner's instruction), then the G4 writes. The numbers sheet is `reports/NUMBERS_SHEET.md` (`python -m full_run_28092026.numbers_sheet`).

D2 (owner, 2026-10-08): **provider defaults** are the headline. Table 1 holds the eleven at their providers' defaults with a marked block for the three re-run with reasoning on, and the matched family is in `tab:matched`. This is the fallback `00_ORCHESTRATION.md` names, since the error analysis was not re-read for the re-run models. Its section 3 row is left for the orchestrator: another session has that file open.

## The final evaluator

- `answer.py` `8fbe1072e1be` (ANSWER FINAL, commit `c7894ed`), and `milestones.py` `9788a9ac6f41`: ANSWER FINAL's file with one change, made on the owner's decision (2026-10-08). This change is not committed: every store's CONFIG names commit `b8b3358`, dirty, with these hashes. Commit `milestones.py` with the rest.
- **The change.** `milestones.close` treated any two numbers below 1e-12 as equal. At the 1e-9 unit factor, a milestone such as 6.69e-4 was therefore matched by any small number in a response. The floor now applies only to a target of exactly zero.
- **Checks:** 10 of 10 targeted cases; `answer.py` self-test 74/74; symbolic equivalence 69/69.
  - Gold (`GOLD_VALIDATION.md`): every gold answer is correct and every milestone is found in its own gold. The chance floor is 0.117 (was 0.123).
  - Expert study (`SCORER_VALIDATION.md`): unchanged, at answer 0.986 / 0.930 and milestones 0.927 / 0.924 / 0.926.
- **What it moved:**
  - Milestone lists: 25 items on 6 templates (26 milestones added, 1 removed). These are small computed values the builder had taken for restated givens. Items without milestones fall from 70 to 57, and MC now uses 148 templates.
  - E3 changed on 332 main-store responses, and 1 verdict moved (Qwen3-235B, incorrect to correct, through a restored answer target).
  - Record: `RESCORE_DIFF_floor_fix.md`.

## Stores re-scored (G2)

- All 16 live stores were re-scored twice with `repair_round5 --rescore --variant V`: once at ANSWER FINAL and once after the floor fix (logs `scores/rescore_answer_final.log`, `scores/rescore_floor_fix.log`). `score.py` exited 0 every time.
  - Only `answer`, `score`, `e3` and (after the fix) `milestones_required` differ: no step, readability or usable flag changed, and no row exists on one side only.
- From the store before ANSWER FINAL to the final store, over the eleven models in `main`: **276 verdicts moved** (`RESCORE_DIFF.md`, `rescore_diff.py --since 20261008T000000Z`).
  - **198** were raised by the symbolic step, on 5 templates.
  - **77** were lowered and **1** raised by the number rule.
  - These are the preview's 275 (`symbolic/RESCORE_PREVIEW.md`, matching per model and per template) plus the floor fix's 1. The set-aside `qwen3.8-27b` adds 28.
- Judge calls (approved, sent in one process with `judge_stores.py`; 0 failed):

  | run | calls | billed |
  |---|---:|---:|
  | after ANSWER FINAL | 152 | $0.40 |
  | after the floor fix | 196 | **$0.68**, above the $0.20–0.50 estimate (the $0.60 cap was passed by calls already in flight) |
  | total | 348 | $1.08 |

  - The reply store was backed up before each run. Afterwards the old bytes were an exact prefix, with only well-formed lines added.
  - The router has no new job; its rows were rebuilt from its store.
- Judge swap, on the owner's choice "same prompts only": **201** responses compared. Of the 220 drawn, 9 belong to re-run items and 10 have a changed prompt; both groups are left out.

## Integration fixes (code)

- `milestones.py`: the floor (above).
- `clause_variants.py`: the headline copy of the rule takes `answer.OWN_DIGIT_CAP`; the variant is `last_digit_unbounded`, which raises verdicts.
- `paper_results.py`:
  - reads the renamed key and caption;
  - the top five's lost points carry a half point (the partial count is odd);
  - the tolerance phrase handles zero swapped pairs.
- `decoding_table.py`: the matched-settings variants are headed as such.
- `judge_swap.py`: same-prompt pairs only; generated counts and per-model range.
- `paraphrase_kit.py --summary`: writes `PARAPHRASE_REVIEW.md` from the merged `accepted.json` (275 kept). `--score` would have dropped round 5's verdicts.
- `paper_setup.py`:
  - reads the current expert-study sentence;
  - reports its Advanced-versus-Easy claim instead of asserting it;
  - gains the matched configuration's data and two phrases.
- `worked_example.py`: docstring; its self-test passes, and the account is C5's (3 of 7, m(r) = 3/7).
- `docs/check_plan_claims.py`:
  - expects the validation figures 0.986, 0.930 and 0.926;
  - the depth table's milestone count is no longer typed.
- `docs/appendix_certification.py`: rounds 5 and 6; kappa "--" for civil and industrial.
- `docs/appendix_evaluation.py`: current B1 and B3 figures, the swap phrase, the symbolic step and 1% cap; claims reported.
- `docs/appendix_statistics.py`: `tab:levels_agreement` is written only once ratings exist (D7).
- New scripts: `rescore_diff.py`, `judge_stores.py`, `numbers_sheet.py`.

## Scripts run (G3, full resolution, final stores)

- `analyze` (B 10,000, B_TEST 100,000), `coverage_variants`, `clause_variants` (every readable verdict reproduced), `flag_precision`, `depth_model`, `judge_swap`, `residual_incorrect`, `threshold_appendix`, `gold_validation`, `validate_scorer`.
- `expert_kits --score`: B2 137 items; B1 144, scored against the current verdicts; B3 97.
- Decoding tables (12) and trace reviews (16), run on the unchanged traces: clean except the intended gaps (`qwen3.8-27b`, `openbook` v1).
- `layer2.markdown_scan`.
- `paper_results.py --out` (scratch), both headlines, `--repaired`: no STAND-IN or QUICK caption and no SOURCE DRIFT.

## For the writers

1. **Claims the final data no longer supports.** Each script prints these as CLAIM FAILS: 20 from `paper_results.py`, 1 from `paper_setup.py` (Advanced varies more than Easy, not for Qwen3-235B or Claude) and 1 from `appendix_evaluation.py` (most partial verdicts called correct: 18 of 36). In addition, `docs/check_plan_claims.py` fails 47 of the old plan's 123 claims.
2. **Results that changed** (numbers sheet):
   - matched settings: the first five are DeepSeek V4.1 Flash, Claude Sonnet 5, Kimi K3, GLM-5.3-Flash and GPT-5.4 mini (Muse Glimmer 30B is sixth);
   - tau between the orderings is 0.709 [0.636, 0.818];
   - Claude Sonnet 5 against DeepSeek V4.1 Flash on MC: Holm p 0.049 as scored, 0.0005 by matching alone, 1.0 route-adjusted, 0.027 intermediate-only. The pair is separated in every reading but route adjustment, so the separation is route conformity;
   - the depth slope holds for 5 of 11 models (DeepSeek's rests on 12 wrong answers);
   - the level gap holds after Holm (Welch) for gpt-oss-20b and Gemini 3.1 Flash-Lite (default), and for gpt-oss-20b (matched);
   - two paraphrase changes are not within the margin.
3. B3: 5 of its 97 readings are of rulings now settled by matching (genuine matches at factor 1); the current-figures rule keeps them.
4. The paper-text edits that were pending in the tree are committed (`412656a`), and the generated blocks were written on top of them. The writers continue from the committed tree and re-run `--write --text-only` after changing a phrase.
5. `paraphrase_audit.py` predates round 5 and was not re-run.

## Blocks written (G4)

- **Markers placed** (the writers may move a pair; the block follows its markers):
  - `appendices/results.tex`: `tab:matched` and `tab:single_path` after the pairwise paragraph; `tab:scoring_variants` before the coverage paragraph; `tab:coverage_variants`, `tab:judged_steps` and `tab:flag_precision` after `tab:coverage`; `tab:depth_model` after `tab:level_gap`.
  - `appendices/models.tex`: `tab:providers` after `tab:decoding`.
  - `appendices/taxonomy_content.tex`: `tab:area` before the difficulty paragraph.
  - `appendices/scoring.tex`: `\input{sections/appendices/worked_example}` at the end.
- **`paper_results.py --write --text-only --headline default --repaired`:** all 19 blocks.
  - New: `tab:matched`, `tab:single_path`, `tab:coverage_variants`, `tab:scoring_variants`, `tab:providers`, `tab:depth_model`, `tab:flag_precision`, `tab:judged_steps`.
  - Rewritten: `tab:main_results`, `tab:coverage`, `tab:level_gap`, `tab:branch_domain`, `tab:paraphrase`, `tab:experiments`, `tab:errors` and the figure blocks `fig:level_bars`, `fig:error_categories`, `fig:branch_bars`, `fig:domain_radar`.
- **`docs/appendix_statistics.py --write`:** `tab:area`. **`worked_example.py --write`:** `appendices/worked_example.tex` (62 lines, self-test passes).
- **Left for WS-F:**
  - no figure PDF redrawn: `figs/error-categories.pdf` still shows the earlier B2 readings while its caption block now gives the current ones;
  - `fig:domain_heatmap` and `fig:level_bars_all` are not placed;
  - WS-F's captions replace the figure blocks' text later.

## --check outcomes (real tree)

| check | outcome |
|---|---|
| `docs/appendix_statistics.py --check` | **passes** (28 of 28 rows; `tab:area` current) |
| `full_run_28092026/paper_results.py --check --headline default --repaired` | every block placed and current. 162 prose failures: 61 phrases to place, 81 old numbers no phrase generates, 20 claims |
| `full_run_28092026/paper_setup.py --check` | waits for prose: `6_experiments.tex` 7 of 9 phrases, 2 old numbers; `models.tex` 31 of 44, 3 old numbers; 1 claim |
| `docs/appendix_certification.py --check` | waits for prose: 30 of 39 (the six rounds, rounds 5 and 6, the hand checks, the undefined kappa) |
| `docs/appendix_evaluation.py --check` | waits for prose: `5_evaluation.tex` 4 of 5, `scoring.tex` 13 of 24 (the symbolic step, the 1% cap), `validation.tex` 39 of 52; 1 claim |

## Regeneration commands for the writers (from the repository root)

```
python full_run_28092026/paper_results.py --write --text-only --headline default --repaired   # after editing a phrase
python full_run_28092026/paper_results.py --check --headline default --repaired
python docs/appendix_statistics.py --write
python docs/appendix_statistics.py --check
python -m full_run_28092026.worked_example --write
python docs/appendix_certification.py --check
python docs/appendix_evaluation.py --check
python full_run_28092026/paper_setup.py --check
python -m full_run_28092026.numbers_sheet
```

Each script prints the rows and phrases it generates (run it without `--check`). A sentence that waits on a number carries `%% FILL: <what>` until it is written.
