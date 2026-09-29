# Evaluation guide: scoring the full run

How to score the full run's traces, run the two judged stages and compute the results, and how to
run the paraphrase and repeat arms. The reasons behind each rule are in `README.md`,
`ANALYSIS_PLAN.md` and DECISIONS D-134 to D-142.

## The stages

| Stage | Script | What it does | Cost | Writes |
|---|---|---|---|---|
| Answer check, E3 milestones, E4 digit rule | `score.py` | Scores every trace once and keeps each evaluator's raw output | free | `scores/<variant>/<model>.jsonl` |
| E5, the milestone judge | `judge.py` | MiMo-V2.5-Pro on the milestones E3 did not find, one call per trace | about $27 ($19 to $50), 8,032 calls | `scores/<variant>/e5/` |
| The step router | `router.py` | The digit rule's flags, and MiMo on every other step, one batched call per trace | about $100 ($72 to $195), 24,506 calls | `scores/<variant>/router/` |
| Analysis | `analyze.py` | The plan's Q1 to Q5, the sensitivity analyses and the reported tables | free | `results/RESULTS.md`, `results/results.json` |

The costs are the free dry runs' estimates for the main run's eleven models, at 29 September's prices.
The judge is called at OpenRouter's default routing, as in the pilot, so its price depends on the
provider that serves each call. Every reply records its provider and billed cost.

## Before you start

- [ ] `full_run_28092026/pool/` and `traces/` are present. If they are not, restore them from the private backups.
- [ ] The experts' labels are in `evaluator_pilot_17092026/experts_filled_labels/version_2/`, and the pilot's stored rows and judge replies are in `evaluator_pilot_17092026/scores/`. The validations replay them. Both stay local.
- [ ] `OPENROUTER_API_KEY` is set in the repo's `.env`. Only the paid stages need it.
- [ ] The spend for the paid stage you are about to run is approved.
- [ ] The laptop is plugged in and set not to sleep. E5 takes about 7 hours and the router about 33.
- [ ] No evaluator changes between a stage's validation and the end of its paid run.

Run every command from the repository root.

## Steps

**1. Deterministic scoring. Free.** This has been run for the main run, and its results are D-140.
Run it again only if the traces or the evaluators change.

```bash
python -m full_run_28092026.score --gold                    # all 2,250 gold solutions
python -m full_run_28092026.score --variant main --workers 8
python -m full_run_28092026.score --status                  # 2,250 rows per model
python -m full_run_28092026.validate_scorer                 # writes SCORER_VALIDATION.md
```

`--gold` must print 2,250 correct at all three tolerances, 0 unusable and 0 digit flags.
`SCORER_VALIDATION.md` must end with "The gold is clean and the published code reproduces every figure".

**2. Tag the evaluator code**, at the commit the paid stages start from, as was done for inference:

```bash
git tag full-run-evaluation
git push origin full-run-evaluation
```

**3. E5, the milestone judge.**

```bash
python -m full_run_28092026.judge --validate                     # free
python -m full_run_28092026.judge --dry-run                      # free: calls, cost, time
python -m full_run_28092026.judge --yes --max-usd 50 --workers 16 > full_run_28092026/scores/e5.log 2>&1
```

- `--validate` replays the pilot's 300 traces from its stored replies and writes `E5_VALIDATION.md`.
  It must end with "The stage reproduces the pilot". If it does not, do not run the paid stage.
- The paid run writes `scores/main/e5/` when it finishes, and prints `without_reply` per model. That
  count should be 0.

**4. The step router.**

```bash
python -m full_run_28092026.router --validate                    # free
python -m full_run_28092026.router --dry-run                     # free
python -m full_run_28092026.router --yes --max-usd 200 --workers 8 > full_run_28092026/scores/router.log 2>&1
```

- `--validate` writes `ROUTER_VALIDATION.md`, which must end with "The stage reproduces the pilot".
- Run the router after E5, not beside it: both call the same judge, and its providers refuse traffic
  beyond a point. The pilot had to drop to 4 workers. If many calls fail, re-run with `--workers 4`.

**Resuming either stage.** Re-run the same `--yes` command. Only the calls not already in the
reply store are made, so nothing answered is paid for twice. `--score` rewrites a stage's rows from
the stored replies without calling anything:

```bash
python -m full_run_28092026.judge --score
python -m full_run_28092026.router --score
```

**5. The analysis. Free.**

```bash
python -m full_run_28092026.analyze --selftest
python -m full_run_28092026.analyze
```

It writes `results/RESULTS.md` and `results/results.json`, aggregates only, and those two files are
committed. Q3's E5 and router columns fill once steps 3 and 4 have written their rows; until then
they show "pending".

**6. Save the paid results locally and to Kaggle.** The judge's replies in `scores/_judge/` are the
paid-for result, so keep two private copies of `scores/`, as for the traces (`INFERENCE_GUIDE.md`,
step 6). First, a local archive with its checksum:

```powershell
$d = "$env:USERPROFILE\EngTrace_private_backup"; $z = "$d\full_run_scores_$(Get-Date -Format yyyy-MM-dd).zip"
Compress-Archive -Path full_run_28092026\scores\* -DestinationPath $z -Force
(Get-FileHash $z -Algorithm SHA256).Hash | Out-File -Encoding ascii "$z.sha256"
```

Then a private Kaggle dataset, in manual mode, staged in a short path and never with `--public`:

```powershell
$s = "$d\kaggle_scores"; New-Item -ItemType Directory -Force $s | Out-Null; Copy-Item "$z*" $s
[IO.File]::WriteAllText("$s\dataset-metadata.json",
  '{"title": "engtrace-full-run-scores", "id": "<your-kaggle-user>/engtrace-full-run-scores", "licenses": [{"name": "other"}]}')
cd $s; kaggle datasets create -p .
```

Download it back and compare every file's hash, as the inference guide does for the traces.

## The paraphrase arm (Q5)

The plan runs it in the same week as the main run (28 to 29 September), so that the served models
are the same: by about 5 October. The experts' check can run while the paraphrase traces are
generated, because a pair the experts reject is dropped from both arms at analysis.

```bash
python -m full_run_28092026.paraphrase --selftest       # free
python -m full_run_28092026.paraphrase --dry-run        # free: $0.15 to $0.45
python -m full_run_28092026.paraphrase --yes            # BILLS: writes and checks the 450 paraphrases
python -m full_run_28092026.paraphrase_kit --build      # free: one kit per expert
```

**The experts' check.** From `paraphrase/dist/`, send each expert `app.py`, `README.txt` and
`guide.md`, together with their own `kit_<id>/` folder. Each expert runs `streamlit run app.py` and
sends back `<id>.jsonl`. Put the returned files in `paraphrase/returned/`.

**The paraphrase run** is launched like the main run (`INFERENCE_GUIDE.md`, step 4), with the
variant added:

```bash
python -m full_run_28092026.run_traces --variant paraphrase --dry-run           # free: about $49
python -m full_run_28092026.run_traces --variant paraphrase --model <key> --workers 15 --yes
python -m full_run_28092026.score --variant paraphrase
python -m full_run_28092026.paraphrase_kit --score full_run_28092026/paraphrase/returned
python -m full_run_28092026.analyze
```

Q5 is marked provisional until the experts' returns are scored.

## The decoding repeats

One cheap model, 300 items, three times, at about $0.09 a repeat:

```bash
python -m full_run_28092026.run_traces --variant repeat1 --model gemma-4-26b-a4b --dry-run
python -m full_run_28092026.run_traces --variant repeat1 --model gemma-4-26b-a4b --workers 15 --yes
python -m full_run_28092026.score --variant repeat1
```

Then the same for `repeat2` and `repeat3`, and `analyze` prints the decoding-repeats table.

## What is committed and what stays local

| Committed | Local only, gitignored |
|---|---|
| The scripts | `pool/`, `SEED.secret`, `traces/` |
| `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md` | `scores/`: the store, the E5 and router rows, the judge's replies |
| `PARSER_FIX.md`, `MATCH_AUDIT.md`, `WORD_AUDIT.md` | `paraphrase/pool.jsonl`, `attempts.jsonl`, `tasks/`, `dist/`, `returned/`, `accepted.json` |
| `results/RESULTS.md`, `results/results.json` | the experts' labels and returns |
| `paraphrase/manifest.jsonl`, `PARAPHRASE.md`, `PARAPHRASE_REVIEW.md` | |

## If something goes wrong

- **A validation does not reproduce.** An evaluator or an input has changed. Do not run the paid stage until the change is found and recorded in DECISIONS.md.
- **`analyze.py` stops with "the store is not the pool the plan describes".** A model's scores are missing or incomplete. Re-run step 1 for that model.
- **The run stops at the cap.** The calls already running finish, so the spend can pass the cap by about one call per worker. Raise `--max-usd` only for an approved new total, then re-run.
- **A call hangs.** After 900 s it is written as a failure. The process writes its rows and ends, and the next run asks that call again.
- **Error 402.** The account is out of credit. Top it up and re-run.
- **Error 429, or many failed calls.** The judge's providers are refusing traffic. Re-run with fewer workers.
- **"without a reply" is above 0 in RESULTS.md.** Some calls failed. Re-run the stage's `--yes` command, then `analyze`.

## Rules

- **Nothing changes while scores exist.** An evaluator change means re-validating, re-scoring and a DECISIONS entry, as D-137 to D-139 were.
- **Every reported number** is printed by a committed script reading these outputs.
- **The store, the replies and anything holding pool text stay local.** Only aggregates are committed.
- **Every paid step needs its own approval**, after its dry-run estimate.
