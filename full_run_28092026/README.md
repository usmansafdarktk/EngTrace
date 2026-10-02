# full_run_28092026 — the full benchmark run

Started 2026-09-28. This directory holds the frozen item pool for the full run: 150 templates
× 15 instances, 2,250 items (D-114). The roster (D-110) and the evaluation stack (D-105, D-110,
D-113) are decided in DECISIONS.md. The inference harness will live here too.

## What is committed, and what is not

| File | Committed | Holds |
|---|---|---|
| `freeze.py` | yes | builds, verifies and checks the pool |
| `manifest.jsonl` | yes | one row per item: id, template, index, branch, domain, area, level, answer type, SHA-256 over question + NUL + solution, and `repeat_of` on a repeat. No text, no seed |
| `FREEZE.json` | yes | the counts, the selection rule, every replaced index and why, the short templates, the SHA-256 of the seed, the manifest's hash |
| `pool/` | **no** | the question and gold text, one `.jsonl` per module, in `testset/`'s schema plus `item_id`, `instance_index` and `repeat_of` |
| `SEED.secret` | **no** | the 128-bit master seed |

The templates are public, so publishing the seed or the text would publish the items. Both
stay local until publication. At publication the seed is revealed, and anyone can regenerate the
pool and check it against the manifest and the seed commitment. Reviewers get `pool/` in the
ARR supplementary archive and can check it with `--check-files`, which needs no seed.

**Back up `pool/` and `SEED.secret` privately.** Without the seed the pool cannot be
regenerated, and `freeze.py` refuses to draw a new seed while a manifest exists. Backed up
2026-09-28 as a private dataset in the owner's Kaggle account, downloaded back and matched file
for file, and as a local archive with its checksum. After any re-freeze, upload a new version.

## Commands

```bash
python -m full_run_28092026.freeze --verify       # regenerate from the seed: must print VERIFY OK and FILES OK
python -m full_run_28092026.freeze --check-files  # pool/ against the manifest; no seed needed
python -m full_run_28092026.freeze                # rebuild. Only if a template changes, which voids the freeze
```

## The pool as frozen, 2026-09-28 (D-116)

Re-frozen after the two widened chemical templates (D-115) and with the coverage selection
(D-116), from the same seed. The first freeze (D-114) took each template's first 15 acceptable
draws; this one takes 15 round-robin across the reasoning paths and answer forms of the first 100.

| | |
|---|---|
| items | 2,250: 870 Easy, 870 Intermediate, 510 Advanced; 450 per branch |
| distinct questions | 2,250; no template needs a repeat |
| templates with more than one path-and-answer group in their first 100 draws | 101 |
| indices not taken as candidates | 585: 549 repeated questions, 35 display ties, 1 pilot question |
| items changed from the D-114 freeze | 332, in 93 templates: 29 in the two widened templates, 303 from the coverage rule |
| checks at the freeze, commit `03b1dd7` | Layer 0 gate 150 of 150; certification 148 of 150, the two widened templates awaiting round 4; `--verify` byte-identical in a separate process |
| certification after round 4 | 150 of 150 (D-119): round 4 approved both widened templates, and the frozen items are the code they approved |
| manifest SHA-256 | `028c7a637eb4895ed061bcbf79a69c656c57f4adc327a4c81cc4faf6d898e2b3` |
| seed commitment | `cd53376fb0fe90ac807ff369a845a5a45e47c1bc48aa4fd470b709a191ddd64e`, unchanged |

Every figure here is printed by `freeze.py` and recorded in `FREEZE.json`.

**How varied the items are** is measured by `diversity.py`, which writes `DIVERSITY.md` and
`diversity.json`: question skeletons, reasoning paths, answer variants and near-duplicates, on the
pool and over 500 public draws per template. It writes counts, answer labels and item ids only.

## Gold validation (D-120)

`gold_validation.py` runs the deterministic evaluators over the 2,250 gold solutions: every item
regenerates byte-identically, the answer check scores every gold answer correct, E3 finds every
milestone, and the digit rule flags no gold claim (`GOLD_VALIDATION.md`). It writes counts and ids;
the text of anything that fails goes to `pool/_gold_validation_details.txt`, local like the pool.

## The harness (D-121)

Step by step, from the dry run to saving the results locally and to Kaggle: `INFERENCE_GUIDE.md`.

`run_traces.py` runs the roster over the pool: `--dry-run` and `--status` are free; `--check`,
`--calibrate N` and the run bill, and refuse to start without `--yes`. The roster, routing and
ceilings are in `models.json`, with the reasons. Traces go to `traces/`, gitignored.

Inference runs at the commit tagged `full-run-inference`, on ten models: `qwen3-235b-a22b` is skipped
for now, and decoding stays at each provider's default (D-122). Two Qwen candidates added to
`models.json` afterwards, `qwen3-235b-a22b-2507` and `qwen3.8-27b`, run at the tag
`full-run-inference-qwen` (D-128). The 2507 release fills the Qwen slot; `qwen3.8-27b` falls under the
roster rule, as a model the pilot generated with, and its traces are set aside (D-132). The roster
the analysis scores is eleven models, 24,750 traces.

`calibration_estimate.py` re-estimates the run from the traces recorded so far: each model's billed
cost per item, a bootstrap interval, a level-reweighted figure and the hours left; `--providers` lists
the endpoints that served the rows, `--empties` the templates of the empty rows, and `--until` limits
it to the rows written by a given time. It is free and calls nothing. The calibration's figures: D-123;
the stop at the spend line: D-125.

## Scoring and analysis (D-134 to D-136, corrected in D-143 to D-149)

Step by step, from scoring to the judged stages, the analysis, the paraphrase and repeat arms and
the backups: `EVALUATION_GUIDE.md`.

`score.py` runs the deterministic evaluators over the traces once and keeps every raw output, one row
per trace, in `scores/<variant>/<model>.jsonl`, gitignored like the traces: the answer check at three
tolerances and under two sensitivity readings, E3 per milestone with the chance floor against a
sibling item, the digit rule per step. `--gold` runs it on the gold solutions. `validate_scorer.py`
checks it on the gold and on the pilot's 300 expert-labelled traces, where it reproduces all ten
published agreement figures (`SCORER_VALIDATION.md`). A paraphrase or repeat run is scored as its own
variant, against the original items. `CONFIG.json` beside the rows records what scored them: the
commit and tag, a dirty flag, every evaluator file's LF-normalised hash and blob, the inputs' hashes
and each model's trace hash (D-144); a store scored with other code is archived by `--replace`, never
overwritten.

`analyze.py` computes ANALYSIS_PLAN.md's questions from the store and writes `results/RESULTS.md`,
`results/results.json` and `results/per_template.csv`, aggregates only; `--selftest` checks its
statistics on synthetic data. The method choices the plan leaves open are fixed in its docstring
(D-136), and the corrections and additions made after the first results were read are labelled where
they appear (D-146, D-148, D-149). `--store main_pre_d137` writes the same tables from the scores
before the evaluator fixes, as a labelled record (D-140).

The first results showed three defects in the answer check, each fixed after being measured on the
gold, on the pilot's expert labels and on the full run: LaTeX numbers (D-137, `parser_fix.py`,
`PARSER_FIX.md`), credit from stray digits (D-138, `match_audit.py`, `MATCH_AUDIT.md`), and verdict
words (D-139, `word_audit.py`, `WORD_AUDIT.md`). Two independent reviews then found the last-digit
windows decided by binary rounding at their edge; the inclusive rule was adopted the same way, and
the half-unit and whole-trace readings are reported as sensitivities (D-147, `boundary_audit.py`,
`BOUNDARY_AUDIT.md`). A review of the finished run read a sample of its verdicts and measured three more
misreadings the pilot's templates never exercised: pi-fractions, one quantity stated in two units on a
scalar gold line, and a pi-fraction whose numerator and denominator are the computed milestones (D-168,
`answer_form_audit.py`, `ANSWER_FORM_AUDIT.md`). All three were adopted after every changed verdict was
read, measured between pinned commits (`--fix`, `ANSWER_FORM_FIX.md`), and the stores re-scored (D-169).
`ANSWER_FORM_AUDIT.md` is the design measurement as made at `e116b4f`, before the change, against the store as
it then stood; it is kept as that record and not regenerated (the script refuses without `--readings`).
`ANSWER_FORM_FIX.md` is the record of the change itself. One wrong credit the change introduces is recorded there
and pinned in `answer.py`'s self-test as a known limit.

`residual_incorrect.py` measures how far the verdicts that remain incorrect after D-169 sit from the gold, per model
and per template, and what the wrong answers on `work_isothermal_virial` state (`RESIDUAL_INCORRECT.md`, D-171; the
reading is in `docs/PILOT_AND_FULL_RUN_ASSESSMENT.md`). `--sample N` prints answer segments for reading and writes nothing.

`expert_kits.py` builds the experts' reading request, B1 to B4 of `docs/EVALUATION_NEXT_STEPS.md`, as one kit per
expert under `expert_request/dist/` (gitignored), with `reading_app.py` as the app and `EXPERT_READING_GUIDE.md` as
the instructions, and scores the returned files into `EXPERT_REQUEST.md` (counts only; D-172). `--selftest` builds a
small kit in a temp folder and drives the app through it.

**For the paper:** `RESULTS_PAPER_NOTES.md` says which figures to report, with their sources, what they
support and what they do not (D-170); `PARAPHRASE_PAPER_NOTES.md` does the same for Q5. The analyses still to
do, in order, are in `docs/EVALUATION_NEXT_STEPS.md`.

```bash
python -m full_run_28092026.score --variant main --workers 8   # free
python -m full_run_28092026.validate_scorer                    # free; needs the experts' labels, local
python -m full_run_28092026.analyze                            # free
```

## The judged stages: E5 and the step router (D-142)

`judge.py` is E5: MiMo-V2.5-Pro on the milestones E3 did not find. `router.py` is the step router:
the digit rule's flags, and MiMo on every other step, batched per trace, on the framework's
Tribunal prompt. Both read the score store and write `scores/<variant>/e5/` and
`scores/<variant>/router/`. Their replies are cached by prompt, so none is bought twice, and the
shared call code, with its deadline and spend cap, is in `judge_calls.py`. `--validate` replays
the pilot's 300 labelled traces from its stored replies and reproduces its published figures
(`E5_VALIDATION.md`, `ROUTER_VALIDATION.md`); `--dry-run` prices the full run. Both are free; the
runs bill and need `--yes`.

```bash
python -m full_run_28092026.judge --validate        # free
python -m full_run_28092026.judge --dry-run         # free
python -m full_run_28092026.judge --yes --max-usd 50            # bills, once approved
python -m full_run_28092026.router --validate       # free
python -m full_run_28092026.router --dry-run        # free
python -m full_run_28092026.router --yes --max-usd 200 --workers 8   # bills, once approved
python -m full_run_28092026.analyze                 # free: Q3's E5 and router columns
```

## The two remaining streams, each with its runbook (D-150)

- **Stream 1, the eleven models' traces:** E5, the router, the decoding repeats, an author's reading
  of the digit rule's flags (`flag_sample.py`, `FLAG_REVIEW.md`), the analysis, the backups and the
  record: `EVALUATION_GUIDE.md`.
- **Stream 2, the paraphrase experiment (Q5):** writing and checking the paraphrases, the experts'
  check, the inference on the paraphrases, scoring and analysis: `PARAPHRASE_RUNBOOK.md`. **Done
  2026-10-01** (D-157, D-161, D-162): 277 expert-kept pairs over 115 templates, Q5 final, $39.15.

`trace_review.py --variant <name>` reviews a variant's traces as it reviewed the main run's, checking a
paraphrase row against `paraphrase/manifest.jsonl` and the original's hash, and each model's served id
against the main run's (`TRACE_REVIEW_<variant>.md`).

## Variant runs, paraphrases and the expert check (D-141)

`subsamples.py` fixes the paraphrase subsample (the plan's 450) and the repeat subsample (300).
`run_traces.py --variant paraphrase | repeat1..3` runs them into `traces/<variant>/`, with a dry run
priced on each model's main-run bills. `paraphrase.py` writes and checks the paraphrases (`--dry-run`
and `--selftest` are free; writing bills and needs `--yes`). `paraphrase_kit.py` builds the experts'
kits and scores their returns; `paraphrase_app.py` is the app shipped in them, and
`paraphrase_guide.md` the guide. The paraphrase text, the kits and the returns stay local; the
committed record is `paraphrase/manifest.jsonl` (hashes only), `PARAPHRASE.md` and
`PARAPHRASE_REVIEW.md`. A pair the experts have not yet returned is kept provisionally and counted as
outstanding; only a rejected pair leaves both arms. Before its checks, `paraphrase.py` restores the
original's ASCII notation where the writer used Unicode the original lacks (D-157).

```bash
python -m full_run_28092026.paraphrase --dry-run                                  # free
python -m full_run_28092026.paraphrase --yes --limit 50                           # bills: a pilot, 10 items per branch
python -m full_run_28092026.paraphrase --yes --until-done --workers 16            # bills, once approved
python -m full_run_28092026.paraphrase --check                                    # free: re-checks the stored attempts
python -m full_run_28092026.paraphrase_kit --build                                # free: the kits
python -m full_run_28092026.run_traces --variant paraphrase --dry-run             # free
python -m full_run_28092026.score --variant paraphrase                            # free, after the run
python -m full_run_28092026.paraphrase_kit --score <returned folder>              # free
```

`testset/` is not this pool. It is the 2026-09-17 generation from the default seed, and the
pilot's `freeze.py --verify` rebuilds its slice from it, so it is left as it is.
