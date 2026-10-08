SIGNAL: KAPPA IN (pending: no E1 ratings have come back; `levels/score_levels.py` posts it when they do)
SIGNAL: GRADES IN 2026-10-08 (362 of 362 readings; grades.json (local) and GRADES.md; the check validated, SYMBOLIC_CHECK.md: five templates enabled)
SIGNAL: ANSWER FINAL 2026-10-08 12:34 +0500 (answer.py, milestones.py and the new symbolic_equivalence.py are final; WS-G re-scores every store once. Re-posted: the 12:10 version quoted pool items in a self-test; the data is now synthetic and every verdict is unchanged)

# WS-E report

Phase 1 is done as of 2026-10-08, except E1, whose ratings have not come back. Nothing is committed (not asked) and no
paid call was made. No `.tex` file and no other stream's file was edited. The experts' files stay local and
gitignored, and this report and every generated file carry counts only.

## E1 difficulty ratings
- Kits built 2026-10-07 (`levels/build_kit.py`, seed 20261007) from the current manifest (sha256 `f2c1dd2f1693`,
  after WS-A's repair): 150 templates, 15 experts, 450 ratings. Each expert rates the 30 templates of their branch.
  - What each expert sees: one evaluation-set problem per template with its reference solution, in their own order,
    under opaque codes, with the §3.1 rubric. The current label is hidden.
  - What each expert answers: the level, and whether the problem states the governing formula or names the method.
  - Composition: `levels/RESULTS.md`.
- Returned: none of 15 as of 2026-10-08 (`levels/experts_filled_levels/` does not exist). Whether the kits went out
  is the owner's record.
- Ready for the returns: `python -m template_annotation_23092026.levels.score_levels` writes `agreement.json` and
  `RESULTS.md`, and KAPPA IN follows. It computes:
  - Fleiss' kappa per branch and overall, Gwet's linear and quadratic weighted forms, and bootstrap 95% intervals
    over templates;
  - the majority against the current label (by one step or two), with Cohen's kappa;
  - the formula-stated counts by branch and level.

  Its self-test passes and reproduces Fleiss' 1971 example (0.210).
- `levels/apply_labels.py` runs only on D7.
  - Rule: a template whose two-of-three majority differs from its label takes the majority; a three-way split
    keeps its label.
  - It changes only the `difficulty` cell of `docs/re-implementation-sep/audit/template_inventory.csv` (what
    `difficulty_map()` reads).
  - It refuses to run if `agreement.json` was scored against other labels.
  - It prints WS-G's follow-up: re-run `freeze.py` for the manifest's `level` and `FREEZE.json`'s `by_level`.

  Its self-test on a copy passes.
- D7 recommendation: after the returns.

## E2 symbolic grading
- Sent: 302 items on the nine symbolic templates. These are every incorrect (90) and partial (176) verdict of the
  eleven models, plus 36 correct ones as controls. That makes 362 readings: one expert of the branch per item, and a
  second on 60 items.
- Returned 2026-10-08: 362 of 362, from 6 of 6 experts. All 302 items are graded, and 3 are split between their two
  readers. On the 60 items read twice, the grades agree on 0.950 (Cohen's kappa 0.882).
- Of the 266 incorrect or partial verdicts, 264 have an agreed grade: 200 equivalent to the reference and 64 not
  (one of the 64 unreadable).
- Controls: 27 of 35 graded equivalent. The other 8 are wrong answers that the number rule credits.
- Per template and per model: `symbolic/GRADES.md`. Per item: `symbolic/grades.json`, which stays local like the
  experts' files.

## The symbolic check
- Enabled (`answer.SYMBOLIC_EQUIVALENCE_TEMPLATES`): `autocorrelation_rect_pulse`, `ber_estimation_mary`,
  `cd_dc_system_analysis`, `impulse_response_from_lccde`, `incompressible_continuity`.
- The rule (stated in full in `symbolic/SYMBOLIC_CHECK.md`):
  - An answer is equivalent when it is the reference's quantity to the precision the reference shows: within half
    a unit of the last digit of every decimal the gold writes, carried through the gold's expression.
  - For the bit error rate, a single value is also equivalent when the reference's value rounds to it.
  - An answer that also states a different value is not equivalent.
  - On an enabled template, the check raises an incorrect or partial verdict to correct. It never lowers one.
- Validation against the grades, on the verdicts the check can move (the answer check's incorrect and partial):

  | template | graded | precision | recall |
  |---|---:|---:|---:|
  | `autocorrelation_rect_pulse` | 18 | 1.000 | 1.000 |
  | `ber_estimation_mary` | 116 | 1.000 | 0.979 |
  | `cd_dc_system_analysis` | 40 | 1.000 | 1.000 |
  | `impulse_response_from_lccde` | 39 | 1.000 | 1.000 |
  | `incompressible_continuity` | 30 | 1.000 | 1.000 |
  | all seven templates with a spec | 264 | 1.000 | 0.990 |

  - Over all seven: TP 198, FP 0, FN 2, TN 64. With the controls added: 291 graded, precision 1.000, recall 0.982.
  - The enabling bar: precision and recall both at least 0.9, at least one graded verdict moved, and every
    disagreement explained.
- Disagreements: 4, all explained.
  - Two Claude Sonnet 5 BER answers (the two FN). The expert accepted roundings beyond the reference's printed
    precision.
  - Two phasor-addition controls, on a template that is not enabled.
- The rule was settled on these same grades, so the figures above are not an estimate on unseen answers.
  - The first version had precision 0.909 and recall 1.000. All 22 of its disagreements were answers it accepted
    and the experts did not.
  - Reading those 22 gave the definition above and fixed five reader defects.
  - The B1 and B2 readings, made before the check existed, give an independent view: 30 readings, precision 1.000,
    recall 0.840.
- Sensitivity: with a whole unit of the gold's last digit instead of half a unit, precision is 0.976 and recall
  1.000 over all.
- Golds: on all seven templates with a spec, each of the 15 golds is equivalent to itself. Of the 1,470
  cross-instance pairs, one is accepted: two BER instances, two modulations whose values coincide within the
  precision shown.
- These templates stay by the numbers:
  - `phasor_addition`: the experts graded none of its 18 incorrect or partial verdicts equivalent.
  - `undamped_response_initial_conditions`: none of its 3 graded equivalent.
  - `continuous_to_discrete_conversion` and `ft_esd_rect_pulse`: no spec. The number rule scores none of their 164
    and 165 usable responses incorrect or partial, so there is nothing to move.
- answer.py changed: yes. WS-G must re-score every store.

## ANSWER FINAL: the evaluator changes (D5, D11a to D11d)
- Files:
  - `evaluator_pilot_17092026/evaluators/answer.py`: the file `score.py` imports and hashes. The brief's
    `full_run_28092026/answer.py` does not exist.
  - `milestones.py`, in the same folder.
  - `symbolic_equivalence.py` (new), in the same folder. `answer.py` pins its LF-normalised SHA-256 (`d39cf431…`)
    and refuses to import a different version. The hash of `answer.py` that `score.py` records therefore covers
    the check, and `score.py` is unchanged. `full_run_28092026/symbolic/equivalence.py` re-exports the module.
  - The self-tests of both files use synthetic items only: no question or gold of the pool. Moving them off pool
    items at 12:34 changed the pin and the hash of `answer.py` and no verdict: the gold check, `SCORER_VALIDATION.md`
    and `RESCORE_PREVIEW.md` came out identical. No store had been re-scored with the 12:10 hash.
- D5: the symbolic step above, called at the end of `answer.verdict`.
- D11a: `milestones.numbers` now reads `1.152 × 10⁻¹⁹` written with superscript digits as one number. The milestone
  cache is rebuilt: its values are byte-identical, and its sidecar names the new `milestones.py`.
- D11b: `answer.values` no longer reads an exponent after a bracket as a value (the 2 of `(signal units)^2`).
- D11c: documented, not changed. A response with no answer heading is read from its last `WINDOW` characters, and
  the D11d bound keeps a coarse intermediate there from matching at a unit factor. WS-C3's example (0.5 × 10³ for
  540) is now an incorrect case in `answer.py`'s self-test.
- D11d: `OWN_DIGIT_CAP = 0.01`. The answer's own last digit vouches for at most 1% of the target, min(u(ŷ), 0.01|y|).
  The scoring appendix's match rule must say so (WS-D2).
- Checks, 2026-10-08:
  - `answer.py` self-test: 74 of 74.
  - `symbolic_equivalence` self-test: 69 of 69.
  - `score --gold`: all 2,250 golds correct at all three tolerances, at half unit and on the whole trace; none
    unusable; no digit flag.
  - `validate_scorer.py` passes, and `symbolic/validate.py --selftest` passes.
- Re-validated on the expert study (`SCORER_VALIDATION.md`, "now" column):

  | figure | before | now |
  |---|---:|---:|
  | answer check, non-partial agreement | 0.982 | 0.986 |
  | answer check, three-way agreement | 0.927 | 0.930 |
  | milestone matching, precision | 0.927 | 0.927 |
  | milestone matching, recall | 0.920 | 0.924 |
  | milestone matching, F1 | 0.923 | 0.926 |

  Its new symbolic section: every enabled gold is equivalent to itself, and graded-equivalent answers score correct
  22/22, 95/97, 36/36, 31/31 and 32/32.
- What the re-score changes in the main store (`symbolic/RESCORE_PREVIEW.md`, written by
  `symbolic/rescore_preview.py`; no store is written). The store was scored with the `answer.py` these changes edit:
  its CONFIG.json records that file's hash.
  - 275 verdicts move:
    - 198 raised by the symbolic step: 139 partial → correct, 59 incorrect → correct;
    - 77 lowered by the number rule: 38 correct → incorrect, 26 correct → partial, 13 partial → incorrect.
  - Each change on its own, applied in turn:
    - D11b lowers 14;
    - D11d lowers 107 (121 together with D11b, the count WS-C3 gave for the bound);
    - the symbolic step raises 221. Of these, 23 are verdicts D11b or D11d had lowered.
  - Mean score per model moves by between −0.0064 (gpt-oss-20b) and +0.0073 (deepseek-v4.1-flash).
  - D11a changes the reached milestones of 216 responses: 221 milestones newly reached and 39 no longer reached. Of
    these responses, GLM-5.3-Flash has 70 and Claude Sonnet 5 has 17, the amendment's upper bounds. 123 judge
    prompts change and are sent again at the re-score.

## E3b top-up (E3c not built)
- Dropped because WS-A's repair changed the question (`adiabatic_flame_temperature`, `work_isothermal_virial`):
  Claude Sonnet 5 17, GPT-5.4 mini 2, gpt-oss-20b 1, Gemma 4 none.
- Added from the wrong answers not yet read, by the level rule: Claude Sonnet 5 5 (all it had left; it has 28 wrong
  answers in all), GPT-5.4 mini 2, gpt-oss-20b 1. That is 8 items and 24 readings (che 2, ele 2 and ind 4 per
  expert).
- Returned 2026-10-08: 24 of 24, from 9 experts, merged by `--score`.
- B2 now (`EXPERT_REQUEST.md`, "The current figures"):
  - 148 items, 444 readings, Fleiss' kappa 0.925.
  - Per model: gpt-oss-20b 0.835, Gemma 4 0.965, GPT-5.4 mini 0.924, Claude Sonnet 5 1.000.
  - Claude Sonnet 5: 28 items, "No error" share 0.357 (0.642 on the original 40).
- New B1 and B3 denominators (the left-out items are those whose question changed):
  - B1: 144 items, 288 readings (6 out); three-way agreement with the check's current verdict 0.778.
  - B3: 97 items, 194 readings (3 out).
- After WS-G's re-score, 11 answers in the sample are no longer scored incorrect. They leave it when `--score` is
  re-run: Claude Sonnet 5 9, Gemma 4 1, GPT-5.4 mini 1.
  - For 10 of them the readers' majority was "No error", so the new check now agrees with the experts.
  - For 1 the majority was an error: a Claude BER answer with a calculation slip whose final value still lies
    within the reference's printed precision. The E2 grader also called it equivalent.
  - The sample then has 137 items, 19 of them Claude's.
- No further top-up: the owner asked for no more expert work. `--b2-topup` builds a round 2 if that changes.
- `scored.json` and `EXPERT_REQUEST.md`'s first two sections stay as returned, so `paper_results.py`'s asserts on 480
  readings and 40 per model still hold. The current figures are in `scored_current.json` and the new section.
- E3c: not built, because the owner declined further expert reading (2026-10-08). If D2 makes the matched
  configuration the headline, the error analysis remains a reading of the default-setting responses and the paper
  says so, as it does for B3. The code is written and self-tested (`--b2-matched VARIANT`, `--score-matched VARIANT
  DIR`).

## Files other streams read
The three JSON files below carry the experts' judgements per template or per item, so they stay local and
gitignored, like `expert_request/scored.json`. Other streams read them in this working tree. The count reports are
committed.

- `template_annotation_23092026/levels/agreement.json` (written after the E1 returns):
  - `agreement.{overall,by_branch}.{fleiss,fleiss_ci95,weighted_linear,weighted_linear_ci95,weighted_quadratic,formula_fleiss}`;
  - `majority.{overall,by_branch}`, `formula_stated`;
  - `templates.<id>.{current,majority,ratings,formula_majority}`.
- `full_run_28092026/symbolic/grades.json`:
  - `items["<item_id>|<model>"].{grade,check_label,group,template_id,model,level,branch,counts,readings}`;
  - `adjudicated_credit_by_model`, `agreement`, `by_template`, `by_model`, `controls`.

  WS-C's expert-adjudicated sensitivity reads this file.
- `full_run_28092026/symbolic/validation.json`:
  - `golds.<id>`;
  - `grades.{acting,all_graded,acting_one_unit,all_graded_one_unit,enabled,disagreements,unexplained}`;
  - `readings`;
  - `reach.<id>`: the number rule's verdict against the check, with `would_move_by_model`;
  - `no_spec.<id>`.

  `reach` now computes the number rule's verdict live instead of reading the stored label, so it is the same before
  and after the re-score. Measured that way, the step moves 221 verdicts (198 when measured against the stored
  verdicts).
- `full_run_28092026/symbolic/SYMBOLIC_CHECK.md`: the rule, the validation, and what stays by the numbers.
- `full_run_28092026/symbolic/RESCORE_PREVIEW.md`: the change at ANSWER FINAL against the stored verdicts. Re-run
  after the re-score, it would show nothing moving.
- In `answer.py`: `SYMBOLIC_EQUIVALENCE_TEMPLATES`, `OWN_DIGIT_CAP`, and
  `symbolic_equivalence.apply(label, item, text, enabled, tol, unit, segment)` → `(label, detail or None)`. These are
  the names `clause_variants.py`'s `symbolic_hook()` expects.
- `full_run_28092026/SCORER_VALIDATION.md`: the figures table (the "now" column moved, see above) and the new section
  "The symbolic equivalence check".
- `full_run_28092026/expert_request/scored_current.json` (local): `kinds.{answer,milestone,error}`,
  `kinds.error.composition`, `exclusions`, `topup`.
- `full_run_28092026/EXPERT_REQUEST.md`: the sections "The B2 top-up: what was sent" and "The current figures".

## Open items
For WS-G, in order:
1. No stream runs `score.py` on a store before WS-G's re-score: every store's CONFIG names the old hash of
   `answer.py`.
2. Re-score every store once (`repair_round5 --rescore --variant V`).
3. Judge: D11a changes 123 judge prompts in the main store; the judge's dry run counts the other stores. The new
   calls are paid: dry run first, then the owner's approval.
4. Run `python -m full_run_28092026.expert_kits --score full_run_28092026/expert_request/expert_requests_filled`.
   It finds the top-up's returns in their own folder. `EXPERT_REQUEST.md`'s current figures then show B2 at 137
   items and B1 scored against the new verdicts.
5. Not needed again: the milestone cache (already rebuilt), `validate_scorer.py` and `symbolic/validate.py`. The
   last two read the new code, not the stores' verdicts.

For WS-C:
- C3, `clause_variants.py`:
  - Its copy of `answer.match` applies the 1% cap only in the `last_digit_bounded` variant.
  - After the re-score, the stored verdicts carry the cap, so its reproduction check stops on the first row the cap
    moved.
  - Fix: the cap (`answer.OWN_DIGIT_CAP`) goes into the headline copy, which makes `last_digit_bounded` equal to the
    headline. The old unbounded window can replace it as the variant.
  - No change is needed for D11b (it reaches the script through `answer.values`), and the symbolic hook's names
    match.
- C5, `worked_example.py`: its default item gains two matched milestones under D11a: 3 of 7 reached instead of 1.
  The verdict stays incorrect, but the judge prompt changes, so the account of what was matched and what was judged
  needs re-checking after the re-score. This was seen with a one-off read, not a committed script.
- C4:
  - "40 from each of four models" stays true of the as-returned B2 sample. The current sample (148 items, 137 after
    the re-score) is in `scored_current.json` if the paper moves to it.
  - The claims about Claude's majority labels change again after the re-score; re-run `--score` (item 4 above).
- `docs/check_plan_claims.py` lines 154 to 156 check the old "now" figures (0.982, 0.927, 0.923). Those checks fail
  until updated to 0.986, 0.930 and 0.926.

For Phase 2 (WS-D; no `.tex` was edited here):
- `5_evaluation.tex` line 9 and `appendices/validation.tex` lines 74 to 76 cite 0.927, 0.982 and 0.927 / 0.920 /
  0.923. The figures are now 0.930, 0.986 and 0.927 / 0.924 / 0.926.
- The scoring appendix needs two statements:
  - the match rule's own-digit term is min(u(ŷ), 0.01|y|) (D11d);
  - symbolic answers on the five enabled templates are scored by equivalence to the reference at its printed
    precision, and the other four by the numbers they state (`SYMBOLIC_CHECK.md`).
- `docs/appendix_evaluation.py` still hard-codes 150 and 100 (already on WS-G's list).

For the owner:
- E1: the returns go to `template_annotation_23092026/levels/experts_filled_levels/`. I then score them and post
  KAPPA IN, and D7 follows.
- D5 is implemented as its default: the check is validated and enabled where the grades support it. D11a, D11b and
  D11d are fixed and D11c is documented, as decided on 2026-10-08.
- The optional E6 (the absolute-value clause sample) is not built; it was not asked for.
- Committed and pushed on 2026-10-08, at the owner's request:
  - `template_annotation_23092026/levels/`: scripts, app, guide, `RESULTS.md`, `.gitignore`;
  - `full_run_28092026/symbolic/`: scripts, app, guide, `GRADES.md`, `SYMBOLIC_CHECK.md`, `RESCORE_PREVIEW.md`,
    `.gitignore`;
  - `evaluator_pilot_17092026/evaluators/answer.py`, `milestones.py` and `symbolic_equivalence.py`;
  - `full_run_28092026/expert_kits.py`, `EXPERT_REQUEST.md`, `validate_scorer.py`, `SCORER_VALIDATION.md`;
  - this report.

  Local and gitignored: kits, keyfiles, the experts' returns, and `grades.json`, `validation.json` and
  `agreement.json`. Before committing, a scan of the committed files against the pool found:
  - no question sentence, and no gold line or item id of an instance, apart from `answer.py`'s older self-test
    cases, which were already in the repository;
  - template sentences that every instance shares, which the template code in `data/templates/` already holds;
  - modulation coefficients that every instance of a modulation shares (0.750, 0.583).
