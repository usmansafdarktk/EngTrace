# Stream 2 runbook: the paraphrase experiment (Q5)

What it tests: whether the models exploit the templates' fixed wording. Each of 450 pool items (the
1st, 6th and 11th of every template) gets one paraphrase that keeps every number, unit, symbol and
part in place; an expert of the item's branch confirms it is the same problem with the same answer;
the eleven models answer the paraphrases under the same prompt, settings and scoring as the main run;
and Q5 reports, per model, the paired difference paraphrase minus original with a template-level
interval, plus the ranking's stability across models. A stable score answers the reviewers' "template
exploitation" objection; a drop is a finding about the models. It is also the paper's contamination
test, since the templates have been public since January (NEXT_CYCLE_REVIEW section 6, 9.3 item 6).

Three steps, in this order; 2 and 3 run in parallel once 1 is done:

1. **Write and check the paraphrases** (Mistral Large 3, $0.18 to $0.54; minutes when Mistral's endpoint is open, an hour or more when it throttles).
2. **The experts' check**: kits out, verdicts back (free; the long pole).
3. **Inference on the paraphrases** (about $49) and scoring; then the analysis.

The reasons behind each rule are in `ANALYSIS_PLAN.md` Q5, D-141, D-149 and D-150. Every paid step
needs the owner's approval on its dry-run estimate.

**Timing.** The plan runs this arm in the same week as the main run (28 to 29 September) so the served
models are the same: by about 5 October. `trace_review --variant paraphrase` checks that each model's
served id equals the main run's; a mismatch is a different checkpoint and must be recorded.

## 0. Session start

Run every command from the repository root; on Windows set `PYTHONIOENCODING=utf-8`.

```bash
git status --short                                   # nothing but the untracked figures_oct_12/
git log --oneline -1                                 # at or after 7ebc348
python -m full_run_28092026.freeze --check-files     # FILES OK - 2250 items in pool/ match manifest.jsonl
python -m full_run_28092026.run_traces --status      # 11 models, missing 0: the main run this arm is paired with
python -m full_run_28092026.paraphrase --selftest    # selftest: all pass
python -m full_run_28092026.paraphrase_kit --selftest  # selftest: all pass
python -m full_run_28092026.subsamples               # 450 paraphrase items; 300 repeat items
```

- [ ] `OPENROUTER_API_KEY` is set in `.env`; the writer and the roster both go through OpenRouter.
- [ ] `full_run_28092026/pool/` and `traces/` are present (restore from the private backups if not).
- [ ] The experts' roster is the layer-2 one (`template_annotation_23092026/layer2/build_tasks.roster`, the pilot's fifteen unless `layer2/annotators.json` overrides it): three per branch, so each expert gets ten templates, 30 items.
- [ ] The spend for the step you are about to run is approved.

## 1. Write and check the paraphrases

The writer is `mistralai/mistral-large-2512`, a family on neither the roster nor the judge's side,
routed like the closed models (its only endpoint is Mistral's). The prompt (hashed into every row)
asks for a rewrite in new words that keeps every number, unit, symbol, variable, formula and technical
term as written and the parts in order, adds nothing, and hints at nothing. An attempt passes only if
all six script checks hold: the same numbers as a multiset, every technical token kept, the part
labels in order, word similarity at most 0.75 (a near-copy tests nothing), length 0.7 to 1.5 times the
original, no preamble. Three attempts per item; an item with no passing attempt leaves both arms.

```bash
python -m full_run_28092026.paraphrase --dry-run     # free: selection, writer's endpoint, $0.18 if all pass first time, $0.54 at most
python -m full_run_28092026.paraphrase --yes --limit 40   # BILLS: a pilot of the prompt on the first 40 unresolved items; read --status before going on
python -m full_run_28092026.paraphrase --yes         # BILLS: writes and checks; resumable; retries each call through the upstream throttle
python -m full_run_28092026.paraphrase --status      # items attempted, passing, billed
```

What it writes, under `paraphrase/`: `attempts.jsonl` (every attempt, with text; local),
`pool.jsonl` (the passing paraphrase per item; local), `manifest.jsonl` (per item its id, the
original's and the paraphrase's SHA-256, the attempt that passed, the hash of the prompt that wrote it
and each check's result; no text; **committed**), and `PARAPHRASE.md` (counts; **committed**). `--check` re-runs every check over the
stored attempts and rewrites the manifest and pool, free.

Expected: most of the 450 pass at the first attempt. Read `PARAPHRASE.md`: the number passing, the
failed attempts by check, the similarity range. If more than a few items have no paraphrase after three
attempts, look at which check they fail (`attempts.jsonl` holds each attempt's check results) before
spending on more attempts; a template whose questions the writer cannot rephrase within the rules is
worth a sentence in the paper, not a fourth attempt.

Commit `paraphrase/manifest.jsonl` and `PARAPHRASE.md`. Record the run in DECISIONS: items written,
passing, billed.

## 2. The experts' check

```bash
python -m full_run_28092026.paraphrase_kit --build   # free: paraphrase/tasks/ and paraphrase/dist/
```

`--build` deals each branch's templates to its three experts in turn, so an expert judges whole
templates: ten templates, up to 30 items, in an order shuffled per expert. It writes `paraphrase/dist/`
with `app.py`, `README.txt`, `guide.md` and one `kit_<id>/` folder per expert holding only opaque codes
and the text they need; `paraphrase/tasks/keyfile.json` maps codes to items and stays local.

**Send each expert** `app.py`, `README.txt`, `guide.md` and their own `kit_<id>/` folder, and nothing
else. A message that has worked for the earlier rounds:

> Thank you again for the certification rounds. One more short task, about an hour: 30 problems
> from your branch, each shown as the original and a reworded version. For each, three questions:
> is it the same problem (same givens, same quantities asked, same conditions); does the original's
> answer still answer it exactly; does the rewording add an ambiguity, an error or a hint. The app
> shows plain text, saves after every item, and can be paused. Unzip the folder, `pip install
> streamlit`, `streamlit run app.py`, pick your id in the sidebar. When it says you are done, send me
> the file `<your id>.jsonl` from the same folder. Please do not discuss items with the other
> reviewers. If possible by <date>.

The app records timestamps, blocks a submission until all three questions are answered, and requires
a note for any answer that rejects the pair.

**When files come back**, put each `<id>.jsonl` in `paraphrase/returned/` and score:

```bash
python -m full_run_28092026.paraphrase_kit --score full_run_28092026/paraphrase/returned
```

It refuses a code that is not in the build's keyfile or a file whose annotator was not assigned that
code, writes `paraphrase/accepted.json` (local: every assigned pair, returned or not, and whether it
was kept) and `PARAPHRASE_REVIEW.md` (**committed**: assigned, returned, kept, rejected, not returned,
per branch, the answers behind the rejections, the median seconds per item). A pair is kept when
"same problem" is yes, "same answer" is yes and "adds ambiguity" is no; any other answer rejects it and
it leaves both arms at analysis. A pair not yet returned is kept provisionally and counted as
outstanding; Q5 stays marked provisional until none is outstanding. Re-run `--score` as more files
arrive; the last submission per item counts.

Commit `PARAPHRASE_REVIEW.md` each time it changes. Record in DECISIONS the returns and the rejections;
the rejected items' ids are in `accepted.json`, and a rejection concentrated in one template is worth
reading (the expert's note says why).

## 3. Inference on the paraphrases, and scoring

The harness runs the same models, prompt, decoding settings, routing and row states as the main run,
with the paraphrase in place of the question; a row records the paraphrase's hash and the original's,
and the run refuses unless `paraphrase/pool.jsonl` matches the committed manifest. The two inert
entries in `models.json` (the original `qwen3-235b-a22b`, the set-aside `qwen3.8-27b`) are never
called by a bare run.

```bash
python -m full_run_28092026.run_traces --variant paraphrase --dry-run     # free: about $49 on the main run's bills for these items
```

Then one process per model, as the main run was launched (`INFERENCE_GUIDE.md` step 4):

```powershell
New-Item -ItemType Directory -Force full_run_28092026\traces\paraphrase | Out-Null
$models = "gpt-oss-20b","gemma-4-26b-a4b","deepseek-v4.1-flash","qwen3-235b-a22b-2507","glm-5.3-flash",
          "glm-5.3","muse-glimmer-30b","kimi-k3","gpt-5.4-mini","gemini-3.1-flash-lite","claude-sonnet-5"
foreach ($m in $models) {
  Start-Process python -ArgumentList "-m full_run_28092026.run_traces --variant paraphrase --model $m --workers 15 --yes" `
    -RedirectStandardOutput "full_run_28092026\traces\paraphrase\$m.log" `
    -RedirectStandardError "full_run_28092026\traces\paraphrase\$m.err" -NoNewWindow
}
```

Watch and finish:

```bash
python -m full_run_28092026.run_traces --variant paraphrase --status   # complete when every model shows missing 0
python -m full_run_28092026.trace_review --variant paraphrase          # writes TRACE_REVIEW_paraphrase.md
```

The review must show, per model: items = the number of passing paraphrases, 2 final rows 0, not a
frozen item 0, prompt "pilot", one request set, one served model and **"as main" yes**. A model served
by a different checkpoint than the main run breaks the pairing for that model; record it and report
that model's Q5 row with the caveat. The dearest models are Kimi K3 (about $19.6), Claude Sonnet 5
($11.8) and GLM-5.3 ($10.9); an empty row scores 0 as in the main run. Only unrun items and service
failures are called again on a re-run.

**Back up the traces** as the main run's were (`INFERENCE_GUIDE.md` step 6), with the subfolder:

```powershell
$d = "$env:USERPROFILE\EngTrace_private_backup"; $z = "$d\full_run_traces_paraphrase_$(Get-Date -Format yyyy-MM-dd).zip"
Compress-Archive -Path full_run_28092026\traces\paraphrase\*.jsonl -DestinationPath $z -Force
(Get-FileHash $z -Algorithm SHA256).Hash | Out-File -Encoding ascii "$z.sha256"
```

Then score and analyse (both free):

```bash
python -m full_run_28092026.score --variant paraphrase      # against the ORIGINAL items: gold, milestones and question are the original's
python -m full_run_28092026.analyze                          # Q5 fills; provisional until the experts' check is complete
```

Q5 prints, per model, the paired answer-score difference with its template-level interval, the
sign-flip test (Holm across models), McNemar as a check, the paired E3 coverage difference on the items
with milestones, the E5-strict difference once E5 has run on both arms, and the difference on the
pairs both arms served from the same endpoint; across models, Kendall's tau between the two arms with
its interval, read against the noise floor in the sensitivity section.

Commit `TRACE_REVIEW_paraphrase.md`, `trace_review_paraphrase.json` and the three `results/` files.
Record the run in DECISIONS as the main run's was (D-126, D-130): rows, empties, billed in the rows
and by the account, served models against the main run's.

## 4. E5 on the paraphrase arm (optional, about $5)

Only after stream 1's E5 has run on the main traces, since Q5's E5 column is a paired difference and
needs both arms:

```bash
python -m full_run_28092026.judge --variant paraphrase --dry-run
python -m full_run_28092026.judge --variant paraphrase --yes --max-usd 10 --workers 16
python -m full_run_28092026.judge --variant paraphrase --status
python -m full_run_28092026.analyze
```

The judge sees the original question with the paraphrase's trace, as the scorer does, so the two arms
differ only in what the model wrote. Its replies join the same store (`scores/_judge/e5_replies.jsonl`)
and the same cumulative cap.

## 5. What is committed and what stays local

| Committed | Local only, gitignored |
|---|---|
| `paraphrase/manifest.jsonl` (hashes and check results, no text), `PARAPHRASE.md`, `PARAPHRASE_REVIEW.md` | `paraphrase/attempts.jsonl`, `pool.jsonl`, `tasks/`, `dist/`, `returned/`, `accepted.json` |
| `TRACE_REVIEW_paraphrase.md`, `trace_review_paraphrase.json` | `traces/paraphrase/`, `scores/paraphrase/` |
| `results/RESULTS.md`, `results.json`, `per_template.csv` | the experts' returned files |

## 6. If something goes wrong

- **`paraphrase --yes` reports service failures.** They are not counted as attempts; re-run the same
  command and only those items are asked again. "Rate-limited upstream" is Mistral throttling
  OpenRouter's shared capacity, not the worker count: on 30 September most calls were refused at eight
  workers and at two alike, and the throughput scaled with the workers. A refused call bills nothing, so
  the writer retries each call up to nine times with a jittered backoff (`RETRY_SLEEPS`) before writing
  a service failure; if a run still defers many items, re-run it or wait for the throttle to ease.
- **Most attempts fail the tokens or numbers check.** Look at `tokens_missing` and `numbers_changed` in
  `attempts.jsonl` before anything else. On 30 September the first prompt lost most attempts to
  reformatted notation (LaTeX, `$` delimiters, Unicode superscripts and minus signs, `·` for `*`) and the
  second to acronyms spelt out, "order-2.1" turned into "second-order", `=` written as "equals" and `$`
  delimiters dropped; the prompt now names each. The checks are not to be loosened for this: the two arms
  must differ in wording only. Try a changed prompt on the first items with `--limit` before the rest,
  and restart from scratch (set `attempts.jsonl` aside under a v-name; the gitignore covers
  `attempts*.jsonl`) so that every item has its three attempts under the one prompt the manifest
  records.
- **An item has no paraphrase after three attempts.** It leaves both arms; `PARAPHRASE.md` says how
  many and which check failed. Do not edit a paraphrase by hand: the manifest's hash would no longer
  match and the run would refuse it.
- **`run_traces --variant paraphrase` refuses: "the paraphrase pool does not match its manifest".** The
  pool or the manifest was edited after the other. `paraphrase --check` rebuilds both from the
  attempts; commit the manifest again if it changed.
- **A model's endpoint has gone, or a different checkpoint is served.** The harness skips a model with
  no eligible endpoint and says so; a served id different from the main run's shows as "as main NO"
  in the trace review. Record it; that model's Q5 row is then a comparison across checkpoints and the
  paper must say so or drop it.
- **A returned file names a code not in the keyfile, or an expert not assigned to it.** The kit was
  built again after it was sent (`--build` rebuilds the codes' assignment with the same seed, so the
  codes are stable unless the pool changed). Score against the build the experts received.
- **The experts are slow.** Q5 is reported provisional with the outstanding count; the inference and
  scoring do not wait for them.

## 7. Done means

- [ ] `PARAPHRASE.md` and `paraphrase/manifest.jsonl` committed; the pool and seed of the arm backed up with the traces.
- [ ] All fifteen experts' files returned and scored; `PARAPHRASE_REVIEW.md` committed with 0 outstanding.
- [ ] Traces for all eleven models, `missing` 0, `TRACE_REVIEW_paraphrase.md` with "as main" yes throughout, backed up.
- [ ] `scores/paraphrase/` scored; `results/RESULTS.md` Q5 no longer marked provisional; committed and pushed.
- [ ] DECISIONS entries for the writing run, the experts' returns and the inference run.
