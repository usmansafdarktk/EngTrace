# WS-B: Reasoning-on inference, full set (Phase 1)

## Mission

Evaluate every model that ran without reasoning tokens with reasoning on, over the full evaluation set, so the
paper can report every model at matched settings and the configuration confound disappears. Models: GPT-5.4 mini
and Gemini 3.1 Flash-Lite (the serving API's reasoning-effort parameter, medium, as the 450-instance arm used);
Gemma 4 26B if its endpoint returns reasoning tokens with the parameter (calibrate first); the Qwen Thinking
sibling of Qwen3-235B-2507 as a twelfth model if the roster rule admits it. Then the cascade for the re-run
models: the paraphrase pairs and the repeats, so that the headline can be the matched configuration without a
hole. No tex is edited.

Closes: review 1 M4 (the reasoning part), M9e; review 2 W6, Q7.

## Context (read `PLAN_CONTEXT.md` sections 1 to 3 first)

- At their providers' defaults, four of the eleven models returned no reasoning tokens: Gemma 4 26B,
  Qwen3-235B-2507 (the Instruct variant, which has no thinking mode), Gemini 3.1 Flash-Lite and GPT-5.4 mini.
  Four of the lower five are these models. Both reviews call the tiers a configuration artifact.
- The 450-instance arm `reasoning-medium` (D-180, `full_run_28092026/RESULTS.md`, "C1 and C4") lifted GPT-5.4
  mini by +0.093 FAC and Gemini 3.1 Flash-Lite by +0.016; both endpoints returned reasoning tokens on every row.
  It cost $9.77 with the judge and step-check stages; inference alone $4.44 and $1.10
  (`DECODING_TABLE_reasoning-medium.md`).
- Main-run inference bills: GPT-5.4 mini $7.74, Gemini $2.48, Gemma $0.73, Qwen Instruct $1.18. Stage costs on
  the main run: judge $36.30 for 8,032 calls (one call per response with a milestone matching missed),
  step check $104.16 for 24,506 calls (one per response).
- The owner's instruction of 2026-10-06: "I want inference also re-run on those models on which thinking was not
  enabled properly before; I don't want a configuration mistake to become a big issue on paper."

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, branch `master` only. Never create, switch to, list or inspect
  another branch.
- Other sessions work in this same working tree at the same time. Edit only the files under "Files you own". If a
  step needs another stream's file, write the need into your report and continue.
- Paid API calls: the owner approved a cap of **$90** for this stream. Before every billed launch run the dry
  run, print the estimate, and launch only if the running total stays within the cap; record every bill. Priority
  if the cap binds: the two closed models, then Gemma 4, then the cascade, then Qwen Thinking.
- No commits or pushes unless the owner asks in this session. When asked: one short line, no body, push at once.
- Scratch files go to the session's scratchpad; code lives in the repository.
- Write new code freely within your files. No `.tex` file is edited in Phase 1.
- Finish by writing `docs/mock_review_workstreams/reports/WS-B_report.md` from the template at the end. Post the
  signals at its top as soon as they are true.

## Before you start

- Read: `full_run_28092026/run_traces.py` (the variants, `items_of()`, `--dry-run`, `--calibrate`, pricing,
  resume), `models.json` (the roster, routing, the `anchor` flag), `subsamples.py`, `INFERENCE_GUIDE.md`,
  `EVALUATION_GUIDE.md` (the arms section, the stage commands and their costs), `DECODING_TABLE_reasoning-medium.md`,
  `score.py` (provenance: `INPUT_FILES`, `EVALUATOR_FILES`), `judge.py`, `router.py`, `decoding_table.py`,
  `trace_review.py`, `keepawake_run.ps1`; `docs/inference_pricing/build_pricing_doc.js` (the roster rule D-110 and
  its exclusion list, lines 21 to 30).
- The decision D3 in `00_ORCHESTRATION.md` (which models), and D2 (whether the matched configuration becomes
  the headline, which decides whether the cascade runs; default: run it, it is cheap).
- WS-A will post FREEZE DONE when the two repaired chemical templates' 30 instances are redrawn. If it arrives
  before you launch, you run on the new pool; if after, you re-run those 30 items in each of your stores (B9).

## Files you own

- `full_run_28092026/run_traces.py`, `models.json`, `subsamples.py` (if a cascade subsample helper is needed).
- `full_run_28092026/traces/<your variants>/`, `scores/<your variants>/` and their stage folders; for the Qwen
  Thinking row, `traces/qwen3-235b-a22b-thinking-2507.jsonl` and `scores/main/qwen3-235b-a22b-thinking-2507.jsonl`
  with its stage rows (the main store is per-model files, so this does not collide with WS-A).
- `DECODING_TABLE_<variant>.md`, `TRACE_REVIEW_<variant>.md`, `results/decoding_table_<variant>.json`.
- `full_run_28092026/results/matched_config.json` (new).
- `EVALUATION_GUIDE.md`: append a section on the full-set reasoning variants.

Not yours: `score.py`, `judge.py`, `router.py`, `decoding_table.py`, `trace_review.py` (run only); `freeze.py`,
`manifest.jsonl`, `pool/` (WS-A); `analyze.py`, `paper_results.py` (WS-C); `answer.py` (WS-E); every `.tex`.

## Steps

### B1. Roster check and prices

1. In `docs/inference_pricing/build_pricing_doc.js`, confirm the exclusion list: pilot generators (including
   `gemma-4-31b-it` and `qwen3.8-27b`) and judges. `qwen/qwen3-235b-a22b-thinking-2507` is eligible if it is on
   neither list and was never used to generate or judge anything in `evaluator_pilot_17092026/` (grep it).
2. Record its OpenRouter identifier, the Hugging Face repository, the providers at 8-bit precision or higher, and
   the price, the way the pricing document does.

### B2. Harness code

1. In `run_traces.py` add the variant `reasoning-medium-full`: all 2,250 pool items with their original
   questions, the reasoning-effort parameter at medium, the main run's prompt (assert the prompt hash equals the
   main run's in `--check`), the main run's output ceiling (32,768) and routing rule, for the models given by
   `--model`. Pricing in `--dry-run` from the 450-item arm's measured token statistics (its decoding table), not
   the main run's. Stores: `traces/reasoning-medium-full/<model>.jsonl`, `scores/reasoning-medium-full/`.
2. Add `--only-items PATH` (a text file of item ids, one per line), honored by every variant, for the resume of a
   subset: WS-A needs it for the 30 repaired items, and you need it in B9.
3. Add the cascade variants: `paraphrase-reasoning-medium` (the kept paraphrase pairs with reasoning on) and
   `repeat1-reasoning-medium` to `repeat3-reasoning-medium` (the 300 repeat items), each restricted to the re-run
   models.
4. For the Qwen Thinking model, if eligible: a roster entry in `models.json` (open-weights, the routing rule,
   `added_for: "matched settings"`, the date), run as `--variant main --model qwen3-235b-a22b-thinking-2507`, so it
   is scored and staged like any roster model.
5. `python -m full_run_28092026.run_traces --check` passes. Post **HARNESS READY**.

### B3. Gemma 4 calibration (cents)

`--variant reasoning-medium-full --model gemma-4-26b-a4b --calibrate 20 --yes`. Read the 20 rows: if the
endpoint reports reasoning tokens on most rows (or returns reasoning text), Gemma is in. If it ignores the
parameter, try the documented alternatives for this model family (a provider flag or a thinking mode the model
card names) once; if none works, Gemma stays at its default. Record every parameter tried and the rows' fields,
so the setup paragraph can say "the endpoint offers no reasoning setting" with evidence.

### B4. Dry runs and launch

1. `--dry-run` per model; print the estimate. Expected inference: GPT-5.4 mini about $22, Gemini about $6, Gemma
   $2 to $4, Qwen Thinking $6 to $12.
2. If FREEZE DONE has not arrived and it is past 14:00 on day 1, launch anyway (B9 covers the 30 items later).
3. Launch each model detached with the keep-awake script, one quiet background poller per run (one line per
   event). `Re-run the same command to resume` on any failure. Watch only your own runs.

### B5. Validate each run

`--status` missing 0; finish reasons; empty responses at the ceiling (reasoning models may hit it more often:
count them; the ceiling stays 32,768 for comparability, do not raise it); share of rows with reasoning tokens
(should be about 100%); providers that served; the prompt hash equals the main run's.
`decoding_table.py --variant reasoning-medium-full`, `trace_review.py --variant reasoning-medium-full`.

### B6. Score and stage

1. `score.py --variant reasoning-medium-full` (free). If FREEZE DONE has arrived, the manifest is the new one and
   your store's provenance is current; if not, WS-G re-scores everything at the end (deterministic).
2. `judge.py --variant reasoning-medium-full` and `router.py --variant reasoning-medium-full` (about $13 per
   model). Dry run (`--status` prints the calls and the estimate) before `--yes`.
3. The same for the Qwen Thinking row under `main`, and for the cascade variants.

### B7. Cascade (default: run)

`paraphrase-reasoning-medium` for the re-run models on the kept pairs (277, or WS-A's new count); the repeats
only for the re-run models that are in the repeat set (Gemini 3.1 Flash-Lite, Gemma 4 if it runs). Score and
stage them. Expected: paraphrase about $3, repeats about $2.

### B8. `results/matched_config.json`

For every roster model: `{"model": key, "default_store": "main", "reasoning_store": "reasoning-medium-full" | null,
"reasoning_setting": "effort=medium" | "none offered" | "reasons by default", "evidence": "..."}`, plus the Qwen
Thinking entry with `sibling_of: qwen3-235b-a22b-2507`. WS-C's matched family reads this file.

### B9. After FREEZE DONE: the 30 repaired items

If any of your runs started before FREEZE DONE, archive the rows of the 30 items in your stores (WS-A's
`repair_round5.py --archive --variant reasoning-medium-full` once it exists, or your own equivalent that moves
the rows to `scores/_replaced/round5_*` with a manifest), then `--only-items <the 30 ids> --yes`, then the
stages for those rows. Report the extra bill (a dollar or two).

### B10. Post RUNS SCORED

When every store you own is scored and staged and `matched_config.json` is written.

## Acceptance

- `run_traces.py --check` passes with the new variants and `--only-items`.
- For each run model: `--status` missing 0 on the full set; reasoning tokens on about every row; the prompt hash
  equals the main run's; decoding table and trace review written.
- Stores scored; judge and step-check rows complete; `matched_config.json` written and self-consistent.
- Bills recorded; total within $90.

## Signals

`SIGNAL: HARNESS READY <date time>` after B2; `SIGNAL: RUNS SCORED <date time>` after B10.

## Report template (`reports/WS-B_report.md`)

```
SIGNAL: HARNESS READY <date time>
SIGNAL: RUNS SCORED <date time>

# WS-B report

## Harness
- Variants added, flags added, --check output.

## Roster check
- Qwen Thinking: eligible yes/no and why; identifier, weights, providers, price.

## Gemma 4 calibration
- Parameters tried; rows with reasoning tokens / text out of 20; decision.

## Runs
| model | variant | rows | empty at ceiling | reasoning-token share | providers | inference $ | judge calls / $ | step-check calls / $ |
- Prompt hash check: equal to main (yes/no).
- Cascade: paraphrase rows and bill; repeats rows and bill.
- The 30 repaired items: re-run in your stores (yes/no, bill).

## Stores
- For each: path, CONFIG current (yes/no, against which manifest), stage rows complete.
- matched_config.json: contents summary.

## Facts for the setup paragraph (Phase 2 writes the sentences)
- Which models were run with reasoning on and at what setting; which endpoints offer no setting; the twelfth model if any.

## Open items
```
