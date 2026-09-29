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

## Scoring and analysis (D-134 to D-136)

`score.py` runs the deterministic evaluators over the traces once and keeps every raw output, one row
per trace, in `scores/<variant>/<model>.jsonl`, gitignored like the traces: the answer check at three
tolerances, E3 per milestone, the digit rule per step. `--gold` runs it on the gold solutions.
`validate_scorer.py` checks it on the gold and on the pilot's 300 expert-labelled traces, where it
reproduces all ten published agreement figures (`SCORER_VALIDATION.md`). A paraphrase or repeat run
is scored as its own variant, against the original items.

`analyze.py` computes ANALYSIS_PLAN.md's questions from the store and writes `results/RESULTS.md` and
`results/results.json`, aggregates only; `--selftest` checks its statistics on synthetic data. The
method choices the plan leaves open are fixed in its docstring (D-136). `--store main_pre_d137` writes
the same tables from the scores before the evaluator fixes, as a labelled record (D-140).

The first results showed three defects in the answer check, each fixed after being measured on the
gold, on the pilot's expert labels and on the full run: LaTeX numbers (D-137, `parser_fix.py`,
`PARSER_FIX.md`), credit from stray digits (D-138, `match_audit.py`, `MATCH_AUDIT.md`), and verdict
words (D-139, `word_audit.py`, `WORD_AUDIT.md`).

```bash
python -m full_run_28092026.score --variant main --workers 8   # free
python -m full_run_28092026.validate_scorer                    # free; needs the experts' labels, local
python -m full_run_28092026.analyze                            # free
```

## Variant runs, paraphrases and the expert check (D-141)

`subsamples.py` fixes the paraphrase subsample (the plan's 450) and the repeat subsample (300).
`run_traces.py --variant paraphrase | repeat1..3` runs them into `traces/<variant>/`, with a dry run
priced on each model's main-run bills. `paraphrase.py` writes and checks the paraphrases (`--dry-run`
and `--selftest` are free; writing bills and needs `--yes`). `paraphrase_kit.py` builds the experts'
kits and scores their returns; `paraphrase_app.py` is the app shipped in them, and
`paraphrase_guide.md` the guide. The paraphrase text, the kits and the returns stay local; the
committed record is `paraphrase/manifest.jsonl` (hashes only), `PARAPHRASE.md` and
`PARAPHRASE_REVIEW.md`.

```bash
python -m full_run_28092026.paraphrase --dry-run                                  # free
python -m full_run_28092026.paraphrase --yes                                      # bills, once approved
python -m full_run_28092026.paraphrase_kit --build                                # free: the kits
python -m full_run_28092026.run_traces --variant paraphrase --dry-run             # free
python -m full_run_28092026.score --variant paraphrase                            # free, after the run
python -m full_run_28092026.paraphrase_kit --score <returned folder>              # free
```

`testset/` is not this pool. It is the 2026-09-17 generation from the default seed, and the
pilot's `freeze.py --verify` rebuilds its slice from it, so it is left as it is.
