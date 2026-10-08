# WS-A: Template repair and re-run (Phase 1)

## Mission

Repair the two chemical templates whose questions do not pin the answer, have the three chemical experts
re-certify them (round 5), redraw their 30 evaluation-set instances with the same private seed, re-run the eleven
models and every experiment arm on those instances, re-score every store, and leave the stores ready for
integration. No tex is edited; the facts for the certification appendix go into your report.

Closes: review 1 F1 (the artifact part), M2 (the two templates); review 2 W1 (part), W2, Q2.

## Context (read `PLAN_CONTEXT.md` sections 1 to 3 first)

- The three chemical experts' reading of 2026-10-03 (`full_run_28092026/EXPERT_REQUEST.md`, section B4):
  `work_isothermal_virial`: the wording decides neither the closed-system nor the flow reading, nor which form of
  the truncated virial equation; no single correct answer; every response shown is correct under one reading.
  `adiabatic_flame_temperature`: the method is standard (two of three), but standard heat-capacity sources differ
  by more than the 0.2% tolerance (two of three); some responses shown are correct under some reading.
- History: round 4 (2026-09-28, D-115, `template_annotation_23092026/layer2/fixes_round4.md`) widened the flame
  template (excess air) and the heat-of-reaction template so each yields 15 distinct questions; both were
  re-certified. The virial template was approved in round 2 unchanged; its round-1 rejection (high temperatures)
  was not adopted (D-094). Neither file has changed since 2026-09-28. So this defect is open.
- Today's wrong answers on them (`full_run_28092026/RESIDUAL_INCORRECT.md`): virial 97 across the eleven models
  (38 in the top five; 53 of the 97 state the flow-work reading), flame 89 (19 in the top five; 48 within 5%).
- The evaluation set is private (local zips, a private Kaggle dataset); `FREEZE.json` holds the SHA-256 of the
  seed, not the seed. The seed does not change, so the commitment stands; only the two templates' code changes.
- What the same rows cost today, from the stored bills: main run $15.31 (330 rows; these templates cost 2 to 6
  times an average instance), arms $3.40, paraphrase pairs $2.99, judge and step-check calls $4.61.

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, branch `master` only. Never create, switch to, list or inspect
  another branch.
- Other sessions work in this same working tree at the same time. Edit only the files under "Files you own". If a
  step needs another stream's file, do not touch it: write the need into your report and continue.
- Paid API calls: the owner approved a cap of **$30** for this stream. Before every billed launch run the dry
  run, print the estimate, and launch only if the running total stays within the cap; record every bill in the
  report. Dry runs, checks and re-scoring are free.
- No commits or pushes unless the owner asks in this session. When asked: one short line, no body, push at once.
- The experts' filled files stay local and gitignored; the report carries counts only; no expert names anywhere.
- Scratch files go to the session's scratchpad; scripts whose outputs the paper will report live in the
  repository under the paths named below, each with a `--check` or self-test.
- Write new code freely within your files. No `.tex` file is edited in Phase 1.
- Finish by writing `docs/mock_review_workstreams/reports/WS-A_report.md` from the template at the end. Post the
  signals at its top as soon as they are true.

## Before you start

- The private seed is available locally (find how `freeze.py` reads it; if it lives only with the owner, ask).
- The three chemical experts can take a round-5 kit this week (the owner dispatches; you build).
- Read: `template_annotation_23092026/README.md` (gate commands), `layer2/README.md` (kit and scoring commands),
  `layer2/RESULTS.md` (round 1, the che-3 rejection of virial), `layer2/fixes_round1.md` (virial lines 57 and
  87 to 92), `layer2/fixes_round4.md`, `layer2/round4_checks.md`, `layer2/CERTIFICATION.md`;
  `full_run_28092026/README.md`, `INFERENCE_GUIDE.md`, `EVALUATION_GUIDE.md`, `FREEZE.json` (the selection rule
  text), `freeze.py`, `subsamples.py`, `score.py` (its docstring and `run()`: the provenance rules), `judge.py`,
  `router.py`, `openbook.py`, `paraphrase.py`, `paraphrase_kit.py`, `decoding_table.py`, `backup_archive.py`;
  `docs/EVALUATION_NEXT_STEPS.md` section B4.
- HARNESS READY from WS-B (the `--only-items` filter in `run_traces.py`) is needed at step A8. Until then do
  steps A1 to A7, which need no harness.

## Files you own

- `data/templates/branches/chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.py`: the
  function `template_work_isothermal_virial` and helpers only it uses.
- `data/templates/branches/chemical_engineering/thermodynamics/heat_effects.py`: the function
  `template_adiabatic_flame_temperature` and helpers only it uses. Do not touch shared constants or any other
  template (round 4's lesson: a shared table widened once changed instances other experts had certified).
- `template_annotation_23092026/layer0/*` outputs; `layer2/tasks_round5/`, `dist_round5/`,
  `experts_filled_annotations_round_5/` (local), `RESULTS_round5.md`, `fixes_round5.md`, `CERTIFICATION.md`.
- `full_run_28092026/freeze.py` (a restricted mode if needed), `FREEZE.json`, `manifest.jsonl`, `pool/`,
  `diversity.py` outputs (`diversity.json`, `DIVERSITY.md`).
- `full_run_28092026/repair_round5.py` (new).
- The rows of the 30 repaired items in every store: `traces/<model>.jsonl`, `traces/<variant>/`,
  `scores/<variant>/<model>.jsonl`, `scores/<variant>/e5/`, `scores/<variant>/router/`,
  `scores/main/e5_grok-4-6/`; the archive `scores/_replaced/round5_*/` and `traces/_replaced_round5/`.
- `full_run_28092026/openbook/` items of the two templates; `paraphrase/` for the 6 new instances.
- `DECODING_TABLE*.md`, `TRACE_REVIEW*.md` regenerated; the local backup zips.

Not yours: `run_traces.py`, `models.json` (WS-B); `score.py`, `judge.py`, `router.py` (run only, no edits);
`analyze.py`, `paper_results.py` (WS-C); `answer.py`, `expert_kits.py` (WS-E); every `.tex` and the
`docs/appendix_*.py` scripts (Phase 2).

## Steps

### A1. Repair `work_isothermal_virial`

1. Read the function and its gold trace. Establish from the code which reading the gold computes (closed-system
   reversible isothermal work on the gas, or flow work) and which form of the truncated virial equation it uses
   (pressure-explicit Z = 1 + BP/RT, or the volume form), and how B is obtained (the Pitzer correlation).
2. Change the question so that it states exactly that: the system and process (one reading, named), the form
   of the equation to use (one form, written out), and the quantity asked (work per mole on the gas, with sign
   convention). Change the gold trace's first step so it names the same reading and form. Do not change the
   sampled parameters, their ranges, the constants, or the arithmetic: the certified physics stays.
3. Keep the output format the integrity checks require (numbered steps, one answer marker), the display rules
   (every displayed value bound through its display before use), and the answer line's form.
4. Run the template at 20 seeds and read five instances end to end.

### A2. Repair `adiabatic_flame_temperature`

1. Read the function: which heat-capacity representation the gold uses (the coefficient table), the heats of
   formation, the reference temperature, the iteration and its rounding (to the kelvin), and the excess-air
   sampling from round 4.
2. Add a data block to the question: for every species in the gold's energy balance, the heat-capacity
   coefficients exactly as the code uses them (with the form of the polynomial and its units stated), the
   standard heats of formation, the reference state, and the convention "report to the nearest kelvin". The
   question then pins the data, and 0.2% holds a solver to the given data, not to a textbook's.
3. Keep the sampling (0 to 100% excess air), the format, and the answer line. Check the question length is
   reasonable (it carries a table now) and that the integrity checks still parse every printed calculation.
4. Run at 20 seeds and read five instances end to end, including one at 0% excess air.

### A3. Integrity checks

Run the layer-0 gate at 500 seeds (closure of every printed calculation, determinism across processes, output
format, every seed generates), restricted to the two templates if the gate supports it, else for the corpus;
then the display-tie census (`template_annotation_23092026/layer0/tie_census.py`). Both must be green. Record the
commands and outputs in `fixes_round5.md`.

### A4. LLM screen (optional, cents)

If `template_annotation_23092026/screen/run_screen.py` can run on two templates, run pass 1 on them (about
$0.10) and record the verdicts; if it cannot without a code change, skip and say so.

### A5. Round-5 certification kits

1. Build the tasks the way round 4 did (`fixes_round4.md`, line 82 shows the command):
   `python -m template_annotation_23092026.layer2.build_tasks --only template_work_isothermal_virial,template_adiabatic_flame_temperature --round 5`
   then `make_kits --round 5` (check the exact flags in `build_tasks.py` and `make_kits.py`). The kit: one
   instance to hand-solve before the solution is shown, then the solution, the code and the other instances;
   approve or reject with a reason. Deliver `dist_round5/` to the owner for the three chemical experts.
2. When the three return (into `experts_filled_annotations_round_5/`), run the layer-2 scorer for round 5 to
   `RESULTS_round5.md`; update `CERTIFICATION.md` (the two rows: round 5, the three verdicts, hand check matched).
3. If any expert rejects: revise per the reason, repeat A3 and A5, and tell the orchestrator at once, because
   the freeze (A6) and the re-run (A8) may have to be repeated (about $25).
4. Write `fixes_round5.md`: the defect as the experts stated it on 2026-10-03, the repair, the checks, the
   round-5 verdicts.

Do not wait for the experts before A6 to A8: the freeze and the re-run proceed on the repaired code, and are
repeated only if round 5 asks for a change.

### A6. Re-freeze the two templates

1. Run `freeze.py` to a scratch output first (`--verify` if it has it) and diff the result against the current
   `manifest.jsonl`: only the 30 rows of the two templates may differ (item ids keep their form
   `<template>#<index>`; `sha256` and the questions change). If any other row differs, stop: something else
   changed, report it.
2. If `freeze.py` cannot restrict itself to two templates, add `--only-templates` (you own it), keeping the
   selection rule unchanged (candidates from indices 0 to 99, grouped by reasoning path and answer form,
   round-robin).
3. Write for real: `manifest.jsonl`, `FREEZE.json` (add a `round5` entry: the two templates, the old and new
   sha256 of each of the 30 items, the date, the commit of the repaired code), `pool/` for the two templates.
4. `python -m full_run_28092026.diversity` to update `diversity.json` and `DIVERSITY.md` (the two templates'
   reasoning-path counts may change; note them).
5. Post **FREEZE DONE** in your report (the orchestrator tells WS-B and WS-C).

### A7. Archive the old rows: `repair_round5.py` (new)

Write `full_run_28092026/repair_round5.py` with these modes, each idempotent and logged:

- `--list`: the 30 item ids, from the manifest.
- `--archive [--variant V ...]`: for every store (default: `main`, `reasoning-medium`, `flagship`,
  `flagship-reasoning-medium`, `openbook2`, `openbook`, `tool`, `paraphrase`, `repeat1`, `repeat2`, `repeat3`,
  and any store WS-B names), move the rows whose `item_id` is one of the 30 out of `traces/...jsonl`,
  `scores/<variant>/<model>.jsonl`, `scores/<variant>/e5/<model>.jsonl`, `scores/<variant>/router/<model>.jsonl`
  and `scores/main/e5_grok-4-6/*` into `scores/_replaced/round5_<timestamp>/` and `traces/_replaced_round5/`,
  with a manifest of what moved (store, model, item id, old sha256). The trace rows are what `run_traces.py`
  resumes on: read its resume logic first and confirm it keys rows by `item_id` within `traces/<...>.jsonl`.
- `--verify-manifest`: the diff of A6.1 as a check.
- `--rescore --variant V`: the provenance procedure of A9, scripted so that WS-G can repeat it.
- `--verify-rescore --variant V`: every row of an unchanged item in the re-scored store equals the archived row
  except provenance fields.
- `--check`: self-test on a copy.

### A8. Re-run inference on the 30 items (needs HARNESS READY)

Dry run first, every time; the cap is $30 for the stream.

1. Main run: `python -m full_run_28092026.run_traces --variant main --dry-run` should list exactly 330 missing
   rows (30 per model); then `--yes`. Resumable; keep-awake launch per `INFERENCE_GUIDE.md`.
2. Arms, each on the two templates' subset instances (`subsamples.py`: positions 0, 5, 10 of each template for
   the 450-item subset; 0 and 7 for the repeats):
   - `reasoning-medium` (GPT-5.4 mini, Gemini 3.1 Flash-Lite): 12 rows.
   - `flagship` (DeepSeek V4 Pro, GPT-5.4): 12 rows; `flagship-reasoning-medium` (GPT-5.4): 6 rows.
   - `openbook2`: first rebuild the open-book items for the two templates (`openbook.py --build --version 2`),
     diff the item file so that only their equation blocks changed, read the two blocks; then 18 rows (Claude
     Sonnet 5, GPT-5.4 mini, gpt-oss-20b). `openbook` (version 1) is obsolete: archive only.
   - `tool` (Claude Sonnet 5, GPT-5.4 mini): 12 rows, if the sandbox runner works unattended; if not, report
     and leave the tool experiment at 444 instances.
   - `repeat1` to `repeat3` (gpt-oss-20b, Gemma 4 26B, Qwen3-235B-2507, Gemini 3.1 Flash-Lite): 4 items × 4
     models × 3 = 48 rows.
   - `paraphrase`: write paraphrases for the 6 new subset instances with `paraphrase.py` (Mistral Large 3, the
     scripted checks, three attempts); build the expert kit (`paraphrase_kit.py`) for one chemical expert; after
     the return, run the eleven models on the kept pairs (up to 66 rows). If the expert cannot return by day 3,
     drop the pairs and report the new kept count (277 minus the dropped).
3. WS-B's reasoning-on stores: not yours. WS-B re-runs the 30 items there after FREEZE DONE.

### A9. Re-score with provenance

`score.py` ties every store to the hashes of `manifest.jsonl` and `diversity.json` (`INPUT_FILES`) and of the
evaluator files; after A6 every store's `CONFIG.json` is stale and `score.py` refuses to score in place. It also
refuses `--replace` while `e5/` exists under the store. Script this in `repair_round5.py --rescore`:

1. Move `scores/<variant>/e5/`, `router/` and (main) `e5_grok-4-6/` aside into the round-5 archive.
2. `python -m full_run_28092026.score --variant <variant> --replace` (archives the old store, re-scores every
   row from the traces; deterministic and free).
3. Restore the moved stage files (they now hold only unchanged items, since the 30 items' rows were archived in
   A7); the judge and step-check rows of unchanged items remain valid because their inputs did not change.
4. `--verify-rescore`: every unchanged item's new row equals its archived row except provenance.
5. Repeat for every arm store that holds the two templates' items.

WS-G repeats this once more at the end, after WS-E's `answer.py` is final (its hash is in the provenance too).

### A10. Judge and step-check stages for the new rows

`python -m full_run_28092026.judge --variant <variant>` and `... router --variant <variant>` resume on the items
without rows (30 per model in `main`, the handful per arm). Expected: about $5 in all. Record `--status` before
and after. For `e5_grok-4-6` (the judge-swap sample of 220 responses): count the sampled responses that sit on
the two templates; do not re-judge them; report the count so WS-C drops them from the swap table.

### A11. Tables, reviews, backups

`decoding_table.py` for `main` and each arm (`--variant`), `trace_review.py` likewise; `backup_archive.py` for
the local zips; the Kaggle copy is the owner's manual step (say so in the report).

### A12. Post STORES UPDATED

When A8 to A10 are complete for `main` (the arms may still be finishing), post the signal: WS-E builds the
error-analysis top-up kit from the re-scored main store.

## Acceptance

- Gate green at 500 seeds for both templates; tie census clean; the five read instances of each are correct and
  unambiguous (state the system, the equation form, the data).
- Round 5: three approvals per template, hand checks matched, `CERTIFICATION.md` updated, `fixes_round5.md`
  written.
- `manifest.jsonl` differs from the archived copy in exactly 30 rows; `FREEZE.json` carries the round-5 entry;
  `diversity.json` regenerated.
- `run_traces.py --status`: missing 0 for `main` and for every arm you re-ran; every re-run row answered.
- Every store you touched: `CONFIG.json` current, `--verify-rescore` passes, judge and step-check rows present
  for every answered row of the 30 items.
- Bills recorded; total within $30.

## Signals

`SIGNAL: FREEZE DONE <date time>` after A6; `SIGNAL: STORES UPDATED <date time>` after A12.

## Report template (`reports/WS-A_report.md`)

```
SIGNAL: FREEZE DONE <date time>
SIGNAL: STORES UPDATED <date time>

# WS-A report

## The repair
- virial: what the question now states (reading, equation form, quantity); what changed in the gold trace; nothing else changed (confirm).
- flame: the data block added (species, representation, reference state, rounding); nothing else changed (confirm).
- Gate and tie census: commands and outcomes. LLM screen: run or skipped, verdicts.

## Round 5
- Kits built (date), returned (date); verdicts per template (A/R × 3), hand checks matched; revisions if any.
- Facts for the certification appendix (Phase 2 writes the sentences): rounds now 5; what round 5 re-certified and why, in one sentence without history words.

## The freeze
- Rows changed: 30 (list the item ids); old and new sha256 recorded in FREEZE.json under `round5`.
- diversity.json changes for the two templates (paths before / after).

## Inference and stages
| store | rows re-run | bill | judge calls | step-check calls | notes |
- Paraphrase: pairs written / kept by the expert / dropped; new kept count.
- Tool arm: re-run or not; instance count for the experiment.
- Judge-swap sample: responses on the two templates (count) to drop.

## Stores
- For each store: CONFIG current (yes/no), verify-rescore passed (yes/no), stage rows complete (yes/no).
- repair_round5.py: modes implemented; the exact command WS-G runs to repeat the re-score.

## Open items
- Anything not done, and why; anything another stream must do (file, action).
```
