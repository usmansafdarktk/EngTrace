# Stream 1 runbook: evaluating the eleven models' traces

The full run's inference is complete and its deterministic scoring is done, validated and committed.
This runbook is everything that remains on the traces of the eleven roster models: the two judged
stages (E5, the router), the decoding repeats, an author's reading of the digit rule's flags, the
analysis, the backups and the record. It is written for a fresh session, in the order things should
happen, with the output each command must give before the next is run.

The paraphrase arm (Q5) is a separate stream with its own runbook, `PARAPHRASE_RUNBOOK.md`; the only
point where the two meet is E5 on the paraphrase traces, which needs this stream's E5 first. That
stream is done, E5 on its traces included (2026-10-01, D-157, D-161, D-162); its E5 replies are in
the same store, `scores/_judge/e5_replies.jsonl`, and count toward the same cumulative cap.

The reasons behind each rule are in `README.md`, `ANALYSIS_PLAN.md` and DECISIONS D-134 to D-150.
Every paid step needs the owner's approval on its dry-run estimate; nothing here starts a paid call
on its own.

## 0. What is already done, and what must not change

| Done | Where | Record |
|---|---|---|
| Pool frozen, 150 x 15 from a private seed | `pool/` (local), `manifest.jsonl`, `FREEZE.json` | D-114 to D-119 |
| Inference, eleven roster models plus the set-aside Qwen3.8-27B, 27,000 rows | `traces/` (local), backed up locally and on Kaggle | D-122 to D-133 |
| Deterministic scoring: answer check, E3, digit rule | `scores/main/` (local), re-scored at `2950875` after the D-169 answer-check change, CONFIG clean | D-134 to D-140, D-144, D-147, D-156 to D-160, D-169 |
| Validations: gold, pilot expert labels, E5 and router replay | `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, `BOUNDARY_AUDIT.md` | D-134, D-142, D-147 |
| Analysis of the deterministic stack | `results/RESULTS.md`, `results.json`, `per_template.csv` | D-140, D-146, D-149 |
| Tag `full-run-evaluation` on `7ebc348`, the evaluator code the paid stages start from | `git tag` | D-144 |

**Nothing below changes a template, the pool, a trace file or an evaluator.** An evaluator change means
re-validating, re-scoring with `--replace`, rebuilding the stage rows with `--score`, and a DECISIONS
entry, as D-137 to D-139 and D-147 were done. The judged stages depend on E3 and the digit rule, not
on the answer check, so an answer-check change after them costs no reply; a change to E3, the milestone
derivation or the digit rule changes the prompts and the unanswered ones are bought again.

## 1. Session start: check before anything else

Run every command from the repository root. On Windows set `PYTHONIOENCODING=utf-8` first: two
self-tests print unicode. Expected output is in the comment.

```bash
git status --short                                   # nothing but the untracked figures_oct_12/
git log --oneline -1                                 # at or after 7ebc348
git tag                                              # full-run-inference, full-run-inference-qwen, full-run-evaluation
python -m full_run_28092026.freeze --check-files     # FILES OK - 2250 items in pool/ match manifest.jsonl
python -m full_run_28092026.run_traces --status      # 11 models, missing 0, $247.917; two inert entries named
python -m full_run_28092026.score --status           # main: commit 2950875, 12 models in CONFIG, 2250 rows each
python -m full_run_28092026.analyze --selftest       # selftest: all pass
python -m full_run_28092026.judge --validate         # ends "The stage reproduces the pilot."
python -m full_run_28092026.router --validate        # ends "The stage reproduces the pilot."
```

Also:

- [ ] `OPENROUTER_API_KEY` is set in the repo's `.env`. It is the only key the stages use; the OpenAI and Anthropic keys there are dead (D-121).
- [ ] The experts' labels are in `evaluator_pilot_17092026/experts_filled_labels/version_2/` and the pilot's stored replies in `evaluator_pilot_17092026/scores/`; the validations above read them. Both stay local.
- [ ] The local backups exist: `~/EngTrace_private_backup/full_run_pool_and_seed_2026-09-28.zip` and `full_run_traces_2026-09-29.zip` with their `.sha256` files.
- [ ] The laptop is on mains power with the lid open, and cannot sleep for the length of the stage (section 2). The stages hold no keep-awake request themselves; launch them through `keepawake_run.ps1` (section 3) or set sleep to never.
- [ ] The spend for the step you are about to run is approved by the owner, on the estimate the dry run prints today, not on the figures in this file.

If any check fails, stop: section 9.

## 2. Costs and times, as of 29 September

| Step | Calls | Estimate | Range over MiMo's endpoints | Time |
|---|---:|---:|---|---|
| E5 | 8,032 | $26.56 | $18.59 to $49.71 | 7.3 h at 16 workers; about 29 h at 4 |
| Router | 24,506 | $103.32 | $72.32 to $194.94 | 32.9 h at 8 workers; about 66 h at 4 |
| Repeats, Gemma, 3 x 300 items | 900 | $0.28 | | minutes |

Spent before this stream: about $354 on the rows' basis, $368 by the account (D-131, D-143). The
router is the one step that does not fit the ~$500 round; whether it runs is the owner's call
(Open decisions).

The estimates take MiMo's output per call from the pilot's replies. E5's first 508 replies ran 1.39 times
as long on the mean, which puts E5 near $40 (D-152). `judge --status` and `router --status` print the
per-reply figures beside that basis, so a stage's estimate can be checked in its first hour.

MiMo is called at OpenRouter's default routing, as in the pilot, so the price depends on the
provider that serves each call; every reply records its provider and billed cost, and `--status`
prints the mix. The pilot's providers held 4 workers, not 16; the estimates at 16 and 8 are the
optimistic case.

## 3. E5, the milestone judge

What it does: for every answered trace whose item has milestones E3 did not find, one call asks
MiMo-V2.5-Pro whether each missed milestone is REACHED, NOT_NEEDED or MISSING; E5-strict credits E3's
milestones plus REACHED. It fills Q3's E5 columns and the attribution table.

```bash
python -m full_run_28092026.judge --dry-run                      # free: calls, cost at today's prices, time
python -m full_run_28092026.judge --yes --max-usd 55 --workers 16 > full_run_28092026/scores/e5.log 2>&1
python -m full_run_28092026.judge --status                       # free, any time: the reply store and the rows
```

- **Launching on Windows.** `keepawake_run.ps1` runs the same command under a keep-awake request, which
  lasts as long as the process and changes no power setting, and appends to the log with a header and an
  exit line. Started through WMI, the process outlives the session that launched it:

  ```powershell
  $cl = "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Minimized -File $PWD\full_run_28092026\keepawake_run.ps1 -Log full_run_28092026\scores\e5.log -m full_run_28092026.judge --yes --max-usd 55 --workers 16"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cl; CurrentDirectory = "$PWD" }
  ```

  To stop it, end the `python.exe` whose command line holds `judge --yes`
  (`Get-CimInstance Win32_Process -Filter "Name='python.exe'"`); every returned call is already in the store.
- **The cap.** `--max-usd` is cumulative (D-148): what the reply store already records, replies and
  failures alike, counts toward it, so it bounds the stage across resumes. Set it to the approved
  amount. The dearest endpoint's estimate is $49.71, so a $50 cap may stop a few calls short if
  routing lands there; the calls in flight when the cap is passed still finish, so the spend can pass
  it by about one call per worker.
- **While it runs**, `judge --status` in another shell shows lines in the store, replies, failed keys,
  the billed total over every line, and, once the rows are written, calls without a reply per model.
  The log prints every 50 returns with the cumulative bill.
- **Finish.** The run writes `scores/main/e5/` and prints `without_reply` per model. It should be 0.
  If it is not, re-run the same `--yes` command: only the unanswered calls are made. A key that has
  failed in three runs is left alone and counted (`keys_given_up` in `--status`); a handful is
  acceptable and is reported, hundreds means a provider problem (section 9).
- **Rebuild without calling.** `judge --score` rewrites the rows from the stored replies and is free;
  it is what to run after any re-score of the store.
- **Then:** `analyze` (section 7) fills Q3's E5 columns. Until `without_reply` is 0 the header says the
  stage is incomplete and the rates are over the answered calls only.

## 4. The step router, if funded *(funded and run, 2026-09-30 to 10-02: D-164)*

What it does: the digit rule's flags, plus MiMo on every other step of a trace in one batched call on
the framework's own Tribunal prompt; a step is flagged when the rule flags it or the judge's category
contains "error". It fills Q3's router columns. It is the only component with any signal on
conceptual error behind a correct answer (D-113), and it costs the most.

```bash
python -m full_run_28092026.router --dry-run                     # free
python -m full_run_28092026.router --yes --max-usd 200 --workers 8 > full_run_28092026/scores/router.log 2>&1
python -m full_run_28092026.router --status                      # free
```

- Run it after E5, not beside it: both call MiMo, and its providers refuse traffic beyond a point. If
  many calls fail, re-run with `--workers 4`; expect about 66 hours then.
- Everything said of the launch, the cap, `--status`, the finish and `--score` for E5 holds here, with
  `scores/router.log` for the log; the reply store is `scores/_judge/router_replies.jsonl`.
- A reply cut off at 8,192 tokens is asked again at 16,384; three attempts per call.

## 5. The decoding repeats

**Done 2026-09-30 (D-151).** Four cheap models, by the owner's decision (the plan names one): Gemma 4
26B, gpt-oss-20b, Qwen3-235B-2507 and Gemini 3.1 Flash-Lite, on the 300-item repeat subsample (the 1st
and 8th item of each template, `subsamples.py`), three times at the main run's settings, so the paper can
say what part of a score is sampling noise at the providers' default decoding. About $0.62 a repeat for
the four; the runs billed $1.863 in the rows.

```bash
M="--model gemma-4-26b-a4b --model gpt-oss-20b --model qwen3-235b-a22b-2507 --model gemini-3.1-flash-lite"
python -m full_run_28092026.run_traces --variant repeat1 $M --dry-run
python -m full_run_28092026.run_traces --variant repeat1 $M --workers 15 --yes   # D-151 ran one process per model
python -m full_run_28092026.run_traces --variant repeat1 $M --status
python -m full_run_28092026.trace_review --variant repeat1        # writes TRACE_REVIEW_repeat1.md; "as main" must be yes
python -m full_run_28092026.score --variant repeat1
```

Then the same with `repeat2` and `repeat3`. `--status` must show missing 0 before scoring; a re-run
of the same command fills the gaps. `analyze` then prints the decoding-repeats table: the three
scores on the same items, their SD and range, and the share of items with the same verdict in every
repeat. A repeat needs `--model`: the harness refuses a bare repeat run or status.

## 6. The digit rule's flags: an author reads a sample

Why: the rule's precision (three flags in four real) was measured on the pilot's five models; this
roster writes differently, and its flags are uneven (Qwen3-235B-2507 has 1,547 on 18,785 claims; 16 of
DeepSeek's 24 sit in one template). The pilot's own rule for the checker is "validate on gold, then
read every flag raised on real traces", and reading is what found its last four parser defects. Q3
prints the rule's rates with the pilot's precision beside them; if this roster's is lower, the paper
must carry this figure instead.

```bash
python -m full_run_28092026.flag_sample --draw          # free: 220 flagged claims, 20 per model, seed 0
```

The sample is already drawn and lives in `scores/flag_review/round1/sample.csv` (local, it holds trace
text); `--draw` refuses to draw a round again, so its verdicts are kept. Open the CSV, and for each row recompute
the left side from the numbers shown and ask whether the right side is a correct rounding at the
precision it displays. Fill `verdict` with one of:

- `slip`: the trace's arithmetic is wrong at the displayed digit; a real flag;
- `checker`: the arithmetic is right and the checker misread it (a unit it did not know, a rounding
  chain it should tolerate, a clause split wrongly, a value inherited from the wrong side); say why in
  `note`, since the note is what a fix is built from;
- `unsure`.

**A reader outside the code (D-153)** gets a workbook instead of the CSV: each claim with its step's full
text and the units the checker read, a verdict list, and no model names, so they cannot bias the reading.
It goes with the instructions in `FLAG_READER_INSTRUCTIONS.md`, which `--reader-copy` copies beside it
as `flag_review_instructions.md`. The workbook holds trace text: send it privately, and never commit it
or the filled copy. Use one reader for all 220, not one per branch. The result is a precision per model, and the
models' flags sit unevenly across branches (13 of DeepSeek's 20 are industrial), so a reader per branch
would tie a model's figure to one reader. No branch expertise is needed: each row is arithmetic and
rounding.

```bash
python -m full_run_28092026.flag_sample --reader-copy                  # flag_review_reader.xlsx and flag_review_instructions.md in scores/flag_review/round1/
python -m full_run_28092026.flag_sample --merge <the returned .xlsx>   # its verdicts into sample.csv, by code
```

Then:

```bash
python -m full_run_28092026.flag_sample --score --read-by "an author"  # writes FLAG_REVIEW.md: counts only, committed
```

`--read-by` names who read the flags in the report's title, and the paper must say the same.

**Rounds 1 and 2 are done, and so are the two fixes they forced (D-154, D-156, D-158, D-159).** The same
domain expert read both. In round 1, 109 of the 220 flags were the checker's; the first fix removed 97.
Round 2 measured that fixed rule at 0.752, and the second fix removed 42 of its 51 misreadings. Before
round 3 went out, two review agents checked both fixes, each on its own. What they found wrong is
corrected, and LaTeX control spaces are now read (D-160, `DIGIT_FIX_3.md`). Round 3 is drawn from that
rule: 190 claims, none on a step an earlier round or the agents read, with the reader's workbook and
instructions beside them in `scores/flag_review/round3/`. It measures the rule as it now stands, and by
the owner's rule no round 4 follows. Every command takes `--round 3`:

```bash
python -m full_run_28092026.flag_sample --round 3 --merge <the returned .xlsx>          # scores/flag_review/round3/
python -m full_run_28092026.flag_sample --round 3 --score --read-by "<who read it>"     # FLAG_REVIEW_3.md
```

If the checker verdicts show a pattern, the fix follows D-137's procedure: an audit script that
measures the change on the gold (must stay 2,250 correct), the pilot's 300 traces against the experts
and every full-run trace; adopt it only if the gold stays clean and no pilot verdict moves away from
the experts; then `score --variant main --replace`, `judge --score` and `router --score` if those
stages exist, `analyze`, and a DECISIONS entry. The digit rule feeds the router's prompts (rule C), so
a change to it after the router has run means the changed prompts are bought again.

## 7. The analysis, and what to commit

```bash
python -m full_run_28092026.analyze --selftest          # all pass
python -m full_run_28092026.analyze                     # writes results/RESULTS.md, results.json, per_template.csv
```

- It refuses if a stage's rows were built from another store than the one on disk ("run `judge
  --score`"), or if the store is not the pool the plan describes.
- The header names the stages present and, if any has calls without a reply, marks it incomplete.
- Its provenance block records the commit the script ran at, the store's CONFIG and each stage's
  CONFIG. The commit it names is the one the analysis ran at, so it is the parent of the commit that
  holds the results; that is expected.
- Commit `results/RESULTS.md`, `results/results.json` and `results/per_template.csv`, and any
  regenerated `TRACE_REVIEW*.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, `FLAG_REVIEW.md`. Never
  anything under `scores/`, `traces/` or `pool/`. One short commit line, then push.

## 8. Backups and the record

**Back up the replies after each paid stage.** The judge's replies in `scores/_judge/` are the
paid-for result. First a local archive with its checksum:

```bash
python -m full_run_28092026.backup_archive scores --exclude scores/flag_review scores/_replaced scores/paraphrase
# ~/EngTrace_private_backup/full_run_scores_<date>.zip, its .sha256, every member checked; the expert's labels never leave
```

Run it after the stage ends: a store still being written is reported as differing. `Compress-Archive`
refuses a file another process holds open, which is why the script uses Python's zipfile (D-151).

Then a private Kaggle dataset, in manual mode (auto mode blocks the upload), staged in a short path,
never `--public`:

```powershell
$d = "$env:USERPROFILE\EngTrace_private_backup"; $z = "$d\full_run_scores_$(Get-Date -Format yyyy-MM-dd).zip"
$s = "$d\kaggle_scores"; New-Item -ItemType Directory -Force $s | Out-Null; Copy-Item "$z*" $s
[IO.File]::WriteAllText("$s\dataset-metadata.json",
  '{"title": "engtrace-full-run-scores", "id": "<your-kaggle-user>/engtrace-full-run-scores", "licenses": [{"name": "other"}]}')
cd $s; kaggle datasets create -p .     # later versions: kaggle datasets version -p . -m "<what changed>"
```

Download it back and compare every file's hash, as `INFERENCE_GUIDE.md` step 6 does for the traces.
The repeat traces under `traces/repeat*/` go into a new traces archive the same way:
`backup_archive traces/repeat1 traces/repeat2 traces/repeat3 --name full_run_traces_repeats` (made and
checked 2026-09-30, locally and on Kaggle, D-151).

**Record each paid run in DECISIONS.md** the way D-123 to D-131 record inference: the date and
commit, the dry-run estimate it was approved on, calls made, billed in the rows and by the account
(the difference is retries and stopped calls), `without_reply` and `keys_given_up`, the providers that
served the calls (`--status` prints them), and anything that stopped or resumed the run. Every number
the paper will quote from these stages must be printed by `analyze.py` or `--status`, never typed.

## 9. If something goes wrong

- **A validation does not reproduce.** An evaluator or an input has changed since the tag. Do not run
  a paid stage; find the change (`git status`, `git diff full-run-evaluation`) and record it.
- **`score --status` shows a commit other than `2950875` or DIRTY.** The store was re-scored. Fine if
  recorded in DECISIONS; otherwise find out why before paying for anything built on it.
- **A stage refuses: "store rows no longer match their traces".** The traces changed after scoring.
  Restore them from the backup or re-score the variant.
- **`analyze` refuses: "was built from another score store".** The store was re-scored after a stage
  wrote its rows. Run the stage's `--score` and try again.
- **The run stops at the cap.** Raise `--max-usd` only for an approved new total; the spend already
  in the store counts toward it.
- **A call hangs.** After 960 s it is no longer waited for and a new call takes its worker (D-152); its
  reply is stored if it arrives, and the next run asks again otherwise. About one call in a hundred did
  this on E5's first night. The process writes its rows and ends.
- **Error 402.** The account is out of credit.
- **Error 429, or failures in the hundreds.** MiMo's providers are refusing traffic. Re-run with fewer
  workers (4 held in the pilot). If one provider returns malformed replies, the pilot's precedent is
  to exclude it by name (`e1_panel.PROVIDER_IGNORE`), measured first.
- **A laptop shutdown mid-run.** Nothing is lost: re-run the same `--yes` command; the store holds
  every answered call.
- **`without a reply` stays above 0 after three runs.** Those keys are given up and counted; the
  analysis leaves their traces out of the judged rates and prints the count. Record it.

## 10. Done means

- [x] E5 rows for all eleven models, `without_reply` 0 or a recorded handful; replies backed up. *(2026-09-30, `without_reply` 0, $36.30, D-155; archived locally and on Kaggle with the router's replies, D-164 and D-169)*
- [x] The router rows likewise, or a recorded decision not to run it. *(2026-10-02, 24,503 of 24,506 answered, $104.16, D-164)*
- [x] Three repeat variants scored; `TRACE_REVIEW_repeat*.md` say "as main" yes. *(2026-09-30, four models, D-151)*
- [x] `FLAG_REVIEW.md` committed, and any checker fix it forced done the D-137 way. *(2026-09-30, D-154, D-156)*
- [x] Round 2 read, `FLAG_REVIEW_2.md` committed: the fixed rule's precision on this roster. *(0.752, D-158; the second fix followed, D-159)*
- [x] Both fixes checked by two independent review agents before round 3; their corrections adopted, the control space read. *(D-160, `DIGIT_FIX_3.md`)*
- [x] Round 3 read, `FLAG_REVIEW_3.md` committed: the corrected rule's precision on this roster; no round 4. *(0.905, D-163)*
- [x] `results/RESULTS.md` regenerated with no "incomplete" in its header, committed and pushed. *(2026-10-02; regenerated again after the D-169 re-score and with the D-170 coverage table)*
- [x] A DECISIONS entry per paid run, and the Open decisions table updated. *(D-155, D-164)*
- [x] `scores/` archived locally and on Kaggle, hashes checked. *(2026-10-02, D-164; rewritten and re-uploaded after the re-score, D-169)*

**After the run was closed (D-168 to D-170).** An end-to-end review read a sample of the verdicts and found three
answer-check misreadings the pilot templates never exercised; they were measured, read verdict by verdict, adopted,
measured between pinned commits and the five stores re-scored (`ANSWER_FORM_AUDIT.md`, `ANSWER_FORM_FIX.md`,
`scores/rescore_d169.log`). The results the paper cites are those regenerated at `2950875` or later; the notes for
the paper are `RESULTS_PAPER_NOTES.md`, and the analyses still to do are in `docs/EVALUATION_NEXT_STEPS.md`.

**After the next steps' section A (D-173 to D-178, 2026-10-03).** `analyze.py` gained the coverage comparison, the branch
and level intervals, the detectable-difference columns and the per-template SD distribution; `decoding_table.py`,
`threshold_appendix.py` and `shortcut_audit.py` write their own files; the pilot gained `analysis/lojo.py` and
`analysis/attribution.py`. All free; `analyze --selftest` covers the additions. The results were regenerated from a
clean tree afterwards, so the provenance line names the commit without "dirty".

**After the next steps' section C (D-180 to D-184, 2026-10-03).** `run_traces.py` gained the arms `reasoning-<effort>`, `flagship`, `flagship-reasoning-<effort>`, `openbook` and `openbook2` (with `openbook.py`, two filter versions) and `tool`, all on the originals of the 450-item subsample; the scorer, `trace_review.py`, `decoding_table.py --variant`, E5 and the router take an arm's name and its models from the arm's own files; `analyze.py` reports the paired arms against their base run and the anchors beside the roster; `judge_swap.py` reports the second judge. Every paid arm was approved on its dry run and recorded with its bill; the tool arm is built and priced, not run. The arms' traces, scores and the judge reply stores are archived locally (`full_run_traces_c_arms_2026-10-03.zip`, `full_run_scores_c_arms_and_replies_2026-10-03.zip`, `backup_archive.py`); the Kaggle copy is the owner's manual-mode step.

## 11. The full-set reasoning variants: the matched configuration (2026-10-06)

**Why.** At their providers' defaults, four roster models returned no reasoning tokens on any row of the main run
(`DECODING_TABLE.md`): Gemma 4 26B, Qwen3-235B-2507, Gemini 3.1 Flash-Lite and GPT-5.4 mini. Every model whose
endpoint offers a reasoning setting is run again with it on, on all 2,250 items, so each model can be reported at
matched settings. The setting is OpenRouter's unified reasoning parameter at medium effort, as the 450-item arm
`reasoning-medium` sent it (D-180). Everything else is the main run's: the prompt (`--check` and every launch refuse
unless its hash equals the one on the main run's rows), the 32,768 ceiling and the routing.

| model | setting | evidence |
|---|---|---|
| GPT-5.4 mini, Gemini 3.1 Flash-Lite | effort medium (`reasoning-medium-full`) | the 450-item arm returned reasoning tokens on every row |
| Gemma 4 26B | effort medium (`reasoning-medium-full`) | the calibration below: 20 of 20 rows reasoned |
| Qwen3-235B-2507 | none offered | the Instruct release "supports only non-thinking mode" (its Hugging Face card); OpenRouter lists no reasoning parameter for it |
| the other seven | reasons by default | reasoning tokens on 99.6% to 100% of their main-run rows |

Qwen3-235B-A22B-Thinking-2507 is the Instruct release's Thinking sibling, a different model. It passes D-110's roster
rule and is in `models.json` as an inert entry (`run: false`, `added_for: "matched settings"`). It runs only when named:
`--variant main --model qwen3-235b-a22b-thinking-2507`, as a labelled twelfth row outside the eleven. It was calibrated
on 20 items on 2026-10-07 ($0.298; a median 191 s and 3,930 completion tokens an item, all answered, served by Novita;
the rows are in `traces/_calibration/`) and not run, by the owner's decision: the full set would have billed about
$33.60 plus about $10 of judged stages, and taken about 3 h of inference at 64 workers.

**The variants** (`run_traces.py`):

- `reasoning-medium-full`: the 2,250 items with the parameter on, for `MATCHED_MODELS` or the models `--model` names.
- `paraphrase-reasoning-medium`: the 275 kept paraphrase pairs (the 314 passing paraphrases less the 39 the expert
  rejected, the pairs `analyze.py` keeps), for `MATCHED_MODELS`.
- `repeat1-reasoning-medium` to `repeat3-reasoning-medium`: the 300 repeat items, for the `MATCHED_MODELS` in the
  decoding repeats' set (Gemini 3.1 Flash-Lite and Gemma 4).

**Commands, in order** (one process per model and variant, launched through `keepawake_run.ps1` as in section 3):

```bash
python -m full_run_28092026.run_traces --variant reasoning-medium-full --model gemma-4-26b-a4b --calibrate 20 --yes  # cents
python -m full_run_28092026.run_traces --variant reasoning-medium-full --dry-run        # free: priced from measured bills
python -m full_run_28092026.run_traces --variant reasoning-medium-full --model <key> --yes
python -m full_run_28092026.run_traces --variant reasoning-medium-full --review         # free: exit 0 when complete and clean
python -m full_run_28092026.run_traces --variant reasoning-medium-full --trace-review   # free: TRACE_REVIEW_<variant>.md
python -m full_run_28092026.decoding_table --variant reasoning-medium-full              # DECODING_TABLE_<variant>.md
python -m full_run_28092026.score --variant reasoning-medium-full                       # free
python -m full_run_28092026.judge --variant reasoning-medium-full --dry-run             # free
python -m full_run_28092026.judge --variant reasoning-medium-full --yes --max-usd <store total + approved> --workers 4
python -m full_run_28092026.router --variant reasoning-medium-full --dry-run            # free
python -m full_run_28092026.router --variant reasoning-medium-full --yes --max-usd <store total + approved> --workers 4
# the cascade: the same run, review and score for paraphrase-reasoning-medium and repeat1..3-reasoning-medium;
# E5 on the paraphrase variant only, as the provider-default paraphrase arm has it; the repeats are scored only, as
# the provider-default repeats are
python -m full_run_28092026.run_traces --matched-config                                 # free: results/matched_config.json
```

- `--only-items PATH` (a file of item ids, one per line) restricts every mode of every variant to those items: the
  resume of a subset whose rows were archived.
- `--trace-review` exists because `trace_review.py` reads any variant name that begins `reasoning-` as the 450-item
  subsample and an unknown name as the whole pool. It passes the variant's own item set and runs that script unchanged.
- The judged stages share the judge model with every other stream. Run one judged process at a time across the
  machine, check the other streams' reports first, and remember that `--max-usd` is cumulative over the reply store.
- `--matched-config` refuses to write while a store it names is incomplete or its rows contradict the setting.

**Costs and times** (the bills are the rows' `billed_usd`; `run_traces --status` prints them):

| run | items | inference billed |
|---|---:|---:|
| `--check` calls, five | | $0.008 |
| GPT-5.4 mini, `reasoning-medium-full` | 2,250 | $22.650 |
| Gemini 3.1 Flash-Lite, `reasoning-medium-full` | 2,250 | $5.526 |
| Gemma 4 26B, `reasoning-medium-full` (its 20 calibration rows, $0.032, included) | 2,250 | $3.795 |
| GPT-5.4 mini, `paraphrase-reasoning-medium` | 275 | $2.600 |
| Gemini 3.1 Flash-Lite, `paraphrase-reasoning-medium` | 275 | $0.660 |
| Gemma 4 26B, `paraphrase-reasoning-medium` | 275 | $0.481 |
| Gemini 3.1 Flash-Lite, `repeat1..3-reasoning-medium` | 3 x 300 | $2.196 |
| Gemma 4 26B, `repeat1..3-reasoning-medium` | 3 x 300 | $1.566 |
| **Inference, all** | | **$39.48** |

GPT-5.4 mini took about 35 minutes at 16 workers, Gemini about 11. Gemma took about 3 hours 10 minutes, at about a
minute per item, while its four cascade runs shared its endpoints. Gemma's figures include 26 rows moved to
`traces/_unfinished/` (below) and the calls that replaced them.

**Rows without a finish reason.** Io Net ended about 1% of Gemma 4's completions with no finish reason, the text cut off
mid-sentence far below the ceiling. In these variants `call()` treats that as a provider fault and asks again, as it
does `finish_reason=error` (D-148). Rows written before that rule went to `traces/_unfinished/<variant>/`
(`run_traces --requeue-unfinished`, which refuses a trace file changed in the last two minutes), and the next resume
asked those items again.

**The judged stages** (2026-10-07, 02:07 to 12:32 UTC; every call answered, none failed):

| stage | calls | billed | dry run |
|---|---:|---:|---|
| E5, `reasoning-medium-full` | 2,391 | $7.285 | $7.79 at Xiaomi's endpoint ($5.45 to $14.67) |
| E5, `paraphrase-reasoning-medium` | 281 | $0.833 | $0.91 ($0.64 to $1.72) |
| router, `reasoning-medium-full` | 6,747 | $22.764 | $22.01 to $26.78 ($15.41 to $40.65) |
| **Judged stages, all** | | **$30.88** | |

The repeats are scored only. Inference and stages together: $70.36; with the Qwen Thinking calibration below, $70.66. At 4 workers the judge model answered about
2.6 calls a minute (a call is generation-bound, a median 74 s at about 41 tokens a second, spread over several
providers). It ran at 16, 24 and, for the router, 32 workers without a failure: E5 at about 17 calls a minute at 24
workers, the router at 18 to 24 a minute at 32. A stage stopped mid-run resumes from its reply store: relaunch the
same command, or the chain, and only the unanswered calls are made.

The traces (with `_unfinished` and `_calibration`), the five score stores and the reply stores are archived locally
(`full_run_traces_reasoning_on_2026-10-07.zip`, `full_run_scores_reasoning_on_2026-10-07.zip`, `backup_archive.py`);
a private Kaggle copy (a new dataset, manual mode) was checked by download, 69 of 69 files matching by path and hash.
