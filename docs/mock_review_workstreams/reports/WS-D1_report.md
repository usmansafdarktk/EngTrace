SIGNAL: WS-D1 DONE 2026-10-08 23:45 (+05:00)

# WS-D1 report

Status: every file I own is rewritten. `paper_results.py --check --headline default --repaired` and `paper_setup.py --check` both exit 0, and so do the four `docs/appendix_*.py --check`. Two FILL markers remain, both waiting on E4. Nothing is committed.

## Files rewritten

- **`current_overleaf_project/main.tex`**
  - Title: "\ourdataset: Verifiable Process Evaluation of Engineering Reasoning from Symbolic Templates".
  - A three-sentence Conclusion: the ceiling on the final answer, what the process measures add, next steps.
  - Limitations: nine scope sentences.
  - Ethics: synthetic data; experts (E4 FILL); labels and names withheld; no deployment on this evidence; MIT; `appendix:release`.
  - The commented acknowledgements are removed.
- **`0_abstract.tex`** (199 words), in this order:
  - deterministic checks against judged ones, with the symbolic rule (equivalence on 5 templates, stated numbers on 4);
  - a misstated rule escapes the checks;
  - both configurations;
  - the five strongest within two points, near the ceiling;
  - reasoning lifts three models by 0.04 to 0.12;
  - the error-type sentence, with its sample (137) and levels.
- **`1_intro.tex`**:
  - Reframed premise: in the expert study the answer nearly settles soundness, and 175 of 178 flawed steps behind correct answers are slips.
  - What step checks add; "process evaluation"; the scope sentence.
  - Three contributions: validation 0.930 / 0.958; what the evaluator cannot see; near the ceiling at matched settings.
- **`2_relatedwork.tex`**:
  - Added: dynamic and functional benchmarks, reference-free evaluators, PRMBench, the five science benchmarks.
  - The ThermoQA and FinChain deltas, one sentence each.
  - A pointer to `appendix:related_work`.
- **`3_4_design.tex`**:
  - "follows" curricula; the selection rule in WS-D2's appendix words (ABET, NCEES FE specifications, ASCE, IISE; three-LLM panel; majority vote, three per branch).
  - A `tab:area` pointer and the uncovered core areas.
  - Difficulty: the experts' original labels and the protocol; no agreement figure and no re-rating.
  - Planted-defect scope; six rounds, unanimous; "physically validated".
- **`5_evaluation.tex`**:
  - Deterministic-first opening; the validation figures, which are in-sample except the cross-fitted ε.
  - A worked-example pointer and two-sentence targets ("label").
  - Equation (1) now carries the 1% own-digit cap, min(c·u(ŷ), 0.01|y|), and the symbolic rule is stated.
  - Visible text only; a valid route that bypasses a milestone is not credited; the judged step check runs at provider defaults and is judged, not verified.
- **`6_experiments.tex`**: both configurations (6,750 matched responses; Qwen3-235B-2507 offers no setting); 260 unreadable (1.1%); the unit-of-analysis sentence corrected; Holm families.
- **`6_results.tex`**:
  - Tiers:
    - Defaults: two tiers, 6 against 5, with every cross pair separated. The first five lie within 0.019, and one pair among them separates (DeepSeek V4.1 Flash above GLM-5.3-Flash).
    - The ceiling, with the Advanced means beside it.
    - Matched settings: a graded order. Only gpt-oss-20b is separated from all; the first four are still above GLM-5.3 and below.
  - MC: the Claude–DeepSeek separation is route conformity (p 0.049 as scored, 0.001 by matching, 1.000 route-adjusted).
  - Flags with their precision (0.905; judged 7 of 15); carried precision; the single-path pointer; the repeat models named; the readable-only reordering.
  - Level gap: 2 models (Welch), 6 (permutation), 1 (matched). Depth model: odds ratio 1.23 pooled, 5 of 11 per model.
  - Experiments: four experiments; anchors with intervals; the tool and equations results as bounds.
  - Error analysis: 137 readings at defaults; the ceiling close. No "headroom"; the open-versus-closed sentence is removed.

## `paper_results.py` and `paper_setup.py`

- **Main-text phrases:** all 45 rewritten to the new prose.
- **The 20 CLAIM FAILS:** each one is either replaced by a claim that guards the new wording or removed together with its sentence. New guards cover the tiers under both configurations, the readable-only reordering, the paraphrase bounds, the separation rule, the depth model and the tool use. Still guarded: 0 claims fail.
- **Appendix phrase generators fixed:**
  - FAC pairs, strict FAC, the tolerance swaps;
  - Claude–DeepSeek (the old wording called a p of 0.001 "no longer differ");
  - the radar's lowest domains; the paraphrase margin, τ position and below/above list; the anchors;
  - the near-miss wording; the error-sample sentence (it said GPT-5.4 mini "has 39 wrong answers in all");
  - new results.tex phrases for the judge-swap gap and carried precision (11 of 3,409, 0.3%).
- **Captions:**
  - Table 1 is shorter.
  - `tab:judged_steps` now gives the judged check's own precision (7 of 15) instead of the combined 0.70 / 0.36.
  - All four figure captions carry the D9 note ("GPT OSS 20B in the figure").
  - The error figure says the readings are at the providers' default settings.
  - `ci()` never prints −0.000.
- **`paper_setup.py`:**
  - Its phrases match the new §5.1 and §5.2.
  - The Advanced-variance phrase is now "for nine of the eleven models".
  - The endpoint sentence is generated from `providers.json`.
  - Generated blocks are left out of its number count.
- **Regeneration:** yes (`--write --text-only --headline default --repaired`). Both `--check` runs pass.

## FILL markers left

| File | Line | Needed |
|---|---|---|
| `overleaf_source_04102026/3_4_design.tex` | 110 | E4: who wrote the templates and who certified them (independence sentence) |
| `current_overleaf_project/main.tex` | 278 | E4: who the 15 experts are, recruitment, compensation, consent |

## Bib entries added under `% --- added by WS-D1`

13 keys, each checked against its proceedings page:
- `zhu2024dyval`, `prasad2023receval`, `chen2023theoremqa`
- `mishra2024mathcamps`: arXiv v1; no archival version exists, and the arXiv id now shows a later paper
- `srivastava2024functionalbench`: arXiv-only
- `li2024gsmplus`, `golovneva2022roscoe`, `xia2024reasoneval`, `song2025prmbench`
- `wang2023scibench`, `he2024olympiadbench`, `xu2025ugphysics`, `qiu2025phybench`

## Sentences for WS-D2

- **Already aligned:** your appendix prose now holds every phrase, and both checks pass.
- **If you edit again:** run both checks. The judge-swap and carried-precision sentences are checked in `results.tex`, and the endpoint sentence in `models.tex`.
- **Matching wording:**
  - §3.1 uses your taxonomy-appendix words, so keep them in step.
  - Your difficulty paragraph points to §3.1, which now states the protocol.

## "Process supervision" in files I do not own

None anywhere in `overleaf_source_04102026/`.

## Open items

- **Recorded decisions:** D1 (I used the recommended title), D8 (the Limitations and Ethics follow the recommendation) and D10 (§3.1 follows WS-D2's appendix wording) have no "Chosen" entry in `00_ORCHESTRATION.md`.
- **Page count:**
  - The main text ends on page 8, measured locally with the committed overview figure (aspect 0.30).
  - `figs/engtrace-overview-5branch.pdf` is missing from `overleaf_source_04102026/figs/`.
  - The taller `-large` redraw (aspect 0.45) would push about 0.15 column over the limit.
- **`main.tex` preamble** (outside my sections):
  - `\usepackage[preprint]{acl}` prints names and links; review needs `[review]`.
  - `todonotes` is loaded twice, which causes an option-clash error.
- **Not applicable:** `docs/check_plan_claims.py` checks the old plan, not the paper.
- **WS-F:** `figs/error-categories.pdf` still shows the earlier readings.
