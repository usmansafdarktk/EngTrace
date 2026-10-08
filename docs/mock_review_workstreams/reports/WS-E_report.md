SIGNAL: KAPPA IN: none to report. E1 is set aside; D7 keeps the domain experts' original labels, and the paper describes the labelling protocol without agreement figures until the experts' count files arrive
SIGNAL: GRADES IN 2026-10-08 (362 of 362 readings; GRADES.md; five templates enabled, SYMBOLIC_CHECK.md)
SIGNAL: ANSWER FINAL 2026-10-08 12:34 +0500 (answer.py, milestones.py and symbolic_equivalence.py are final; WS-G re-scores every store once)

# WS-E report

Phase 1 done, 2026-10-08. No paid call. The experts' files and the per-item grade files stay local; the reports
carry counts only. One `.tex` file was edited, at the owner's request: `appendices/taxonomy_content.tex`.

## E1 difficulty ratings
- E1 is set aside and not used.
- D7: the domain experts' original labels stay (58 / 58 / 34), so nothing is relabelled or re-frozen.
- The paper describes the procedure in `template_annotation_23092026/levels/difficulty_labelling_protocol.md`. No
  agreement figure until the experts' count files arrive; a committed script computes it then.
- Domain and area validation: `appendices/taxonomy_content.tex` now names `Grok 4.6`, `MiniMax M3` and
  `MiMo-V2.5-Pro`.

## E2 symbolic grading
- 302 items on the nine symbolic templates: 90 incorrect, 176 partial and 36 correct verdicts as controls.
- 362 readings by 6 experts, all returned. 3 items split. On the 60 items read twice the grades agree on 0.950
  (Cohen's kappa 0.882).
- Of the 264 graded incorrect or partial verdicts, 200 are equivalent to the reference. Controls: 27 of 35.
- Per template and model: `symbolic/GRADES.md`. Per item: `symbolic/grades.json` (local).

## The symbolic check
- Rule (full statement in `symbolic/SYMBOLIC_CHECK.md`): an answer is equivalent when it is the reference's quantity
  at the precision the reference shows. On an enabled template the check raises an incorrect or partial verdict to
  correct, and it never lowers one.
- Against the grades, on the verdicts it can move:

  | template | graded | precision | recall |
  |---|---:|---:|---:|
  | `autocorrelation_rect_pulse` | 18 | 1.000 | 1.000 |
  | `ber_estimation_mary` | 116 | 1.000 | 0.979 |
  | `cd_dc_system_analysis` | 40 | 1.000 | 1.000 |
  | `impulse_response_from_lccde` | 39 | 1.000 | 1.000 |
  | `incompressible_continuity` | 30 | 1.000 | 1.000 |
  | all templates with a spec | 264 | 1.000 | 0.990 |

- 4 disagreements, all explained. The rule was settled on these grades. The earlier B1 and B2 readings give an
  independent check: 30 readings, precision 1.000, recall 0.840.
- These five templates are enabled. The rest stay by the numbers:
  - `phasor_addition` and `undamped_response_initial_conditions`: no graded verdict is equivalent;
  - `continuous_to_discrete_conversion` and `ft_esd_rect_pulse`: there is nothing to move.

## ANSWER FINAL (D5, D11a to D11d)
- Files: `evaluator_pilot_17092026/evaluators/answer.py`, `milestones.py`, and `symbolic_equivalence.py` (new).
  `answer.py` pins the new module's SHA-256, so the provenance hash `score.py` records covers it.
- Changes:
  - D5: the symbolic step;
  - D11a: superscript exponents read by the milestone reader;
  - D11b: an exponent after a bracket is no longer read as a value;
  - D11c: documented, not changed;
  - D11d: an answer's own last digit vouches for at most 1% of the target.
- Checks: self-tests 74 of 74 and 69 of 69; all 2,250 golds correct; `validate_scorer.py` passes.
- Expert-study figures (`SCORER_VALIDATION.md`):
  - non-partial agreement 0.982 → 0.986;
  - three-way agreement 0.927 → 0.930;
  - milestone recall 0.920 → 0.924;
  - milestone F1 0.923 → 0.926.
- What the main store's re-score will change (`symbolic/RESCORE_PREVIEW.md`):
  - 275 verdicts move: 198 raised by the symbolic step, 77 lowered by the number rule;
  - model means change by −0.006 to +0.007;
  - D11a changes the matched milestones of 216 responses, and 123 judge prompts.

## E3b top-up (E3c not built)
- Dropped, because their questions changed in WS-A's repair: Claude Sonnet 5 17, GPT-5.4 mini 2, gpt-oss-20b 1.
- Added: Claude Sonnet 5 5 (all it had left), GPT-5.4 mini 2, gpt-oss-20b 1. All 24 readings returned.
- B2 now: 148 items, 444 readings, Fleiss' kappa 0.925.
  - Per model: gpt-oss-20b 0.835, Gemma 4 0.965, GPT-5.4 mini 0.924, Claude Sonnet 5 1.000.
  - Claude Sonnet 5's "No error" share: 0.357.
- B1 now: 144 items, 288 readings. B3 now: 97 items, 194 readings.
- After the re-score, 11 answers leave B2 (10 of them had a "No error" majority). That leaves 137 items, 19 of them
  Claude's.
- `scored.json` is unchanged. The current figures are in `scored_current.json` and `EXPERT_REQUEST.md`.
- E3c not built (no further expert reading). The code is ready.

## Files other streams read
- `symbolic/grades.json` (local): `items["<item_id>|<model>"]`, `adjudicated_credit_by_model`, `agreement`. WS-C's
  expert-adjudicated sensitivity reads it.
- `symbolic/validation.json` (local): `golds`, `grades`, `readings`, `reach` (the number rule's verdict computed live;
  the step moves 221 verdicts) and `no_spec`.
- `symbolic/SYMBOLIC_CHECK.md`, `symbolic/RESCORE_PREVIEW.md`, `SCORER_VALIDATION.md`, `EXPERT_REQUEST.md`, and
  `expert_request/scored_current.json` (local).
- In `answer.py`: `SYMBOLIC_EQUIVALENCE_TEMPLATES`, `OWN_DIGIT_CAP` and
  `symbolic_equivalence.apply(label, item, text, enabled, tol, unit, segment)`.

## Open items
WS-G:
1. No `score.py` run on any store before the re-score.
2. Re-score every store (`repair_round5 --rescore --variant V`).
3. Judge: 123 changed prompts in the main store; the dry run counts the other stores. The calls are paid and need
   the owner's approval.
4. Re-run `python -m full_run_28092026.expert_kits --score full_run_28092026/expert_request/expert_requests_filled`.

WS-C:
- C3, `clause_variants.py`: put the 1% cap (`answer.OWN_DIGIT_CAP`) in the headline copy, or its reproduction check
  stops after the re-score.
- C5, `worked_example.py`: its item gains two matched milestones under D11a. Re-check it after the re-score.
- `docs/check_plan_claims.py` lines 154 to 156: update to 0.986, 0.930 and 0.926.

WS-D:
- §3.1: describe the labelling procedure from the protocol file, with no agreement figure.
- `3_4_design.tex` line 31 says "four LLMs"; the appendix now names three.
- `5_evaluation.tex` line 9 and `appendices/validation.tex` lines 74 to 76: the new expert-study figures.
- Scoring appendix: state the 1% bound (D11d) and the symbolic equivalence rule.

Owner:
- E2: confirm how the grades were collected before the paper cites the graders' agreement.
