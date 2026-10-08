SIGNAL: NUMBERS SHEET 2026-10-08 16:35 UTC (G2, G2b, G3 and G5 done on the final stores, the milestone matcher's floor fixed; G4 waits for D2 and WS-F's captions)

# WS-G report

Status 2026-10-08 16:35 UTC. INTEGRATED is not posted: G4 has not started, and no `.tex` file in the real tree has been written. No commit (none asked). The numbers sheet is `reports/NUMBERS_SHEET.md` (`python -m full_run_28092026.numbers_sheet`).

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
4. `appendices/validation.tex` holds uncommitted edits from another session. The writers pause while G4 runs `--write`.
5. `paraphrase_audit.py` predates round 5 and was not re-run.

## --check outcomes now

| check | outcome |
|---|---|
| `paper_results.py --check` (scratch copies, both headlines) | fails: markers and phrases not yet placed, and 20 claims |
| `docs/appendix_statistics.py --check` | passes (28 of 28 rows; `tab:area` markers to place) |
| `full_run_28092026/paper_setup.py --check` | 1 claim; `6_experiments.tex`: 2 sentences, stale 246 and 1.0; `models.tex`: 13 rows and phrases, stale 147 |
| `docs/appendix_certification.py --check` | 30 of 39; 9 rows and phrases wait for the certification prose |
| `docs/appendix_evaluation.py --check` | 1 claim; 1 + 11 + 13 rows and phrases wait (`5_evaluation.tex`, `scoring.tex`, `validation.tex`) |

## G4: what remains (once D2 and WS-F's captions are in)

Steps:
1. Place the registry's markers.
2. Paste WS-F's captions into the `figure()` calls.
3. Run `paper_results.py --write --headline <D2> --repaired`.
4. Run `worked_example.py --write`.
5. Run `docs/appendix_statistics.py --write`.
6. Run every `--check`.

About 30 minutes. The writers pause meanwhile.

## Regeneration commands (from the repository root)

```
python -m full_run_28092026.numbers_sheet
python full_run_28092026/paper_results.py --write --text-only --headline <D2> --repaired
python full_run_28092026/paper_results.py --check --headline <D2> --repaired
python docs/appendix_statistics.py --write
python docs/appendix_statistics.py --check
python docs/appendix_certification.py --check
python docs/appendix_evaluation.py --check
python full_run_28092026/paper_setup.py --check
```
