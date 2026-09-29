# Inference guide: running the full run

How to run the roster over the frozen pool with `run_traces.py`. The reasons behind each rule are
in `README.md`, `models.json` and DECISIONS D-114 to D-121.

## Before the first paid call

- [ ] `OPENROUTER_API_KEY` is set in the repo's `.env`. It is the only key the harness uses.
- [ ] `full_run_28092026/pool/` and `SEED.secret` are present. If they are not, restore them from the private backup.
- [ ] The open decisions are settled: Qwen3-235B, the temperature and the analysis plan (DECISIONS, "Open decisions").
- [ ] The spend for the step you are about to run is approved.
- [ ] The laptop is plugged in and set not to sleep. A shutdown cuts calls mid-run.

Run every command from the repository root.

## Steps

**1. Dry run. Free: it makes no model call.**

```bash
python -m full_run_28092026.run_traces --dry-run
```

It confirms that the pool matches the committed manifest and that the prompt matches the pilot's. It
also shows each model's routed endpoint, its output cap and the estimate. A model with no eligible
endpoint is skipped by every paid mode.

**2. Check. About $0.01: one tiny call per model.**

```bash
python -m full_run_28092026.run_traces --check --yes
```

Every model should print `answered`, with its served model id and provider.

**3. Calibration. About $3.50: 20 items per model.**

```bash
python -m full_run_28092026.run_traces --calibrate 20 --yes
python -m full_run_28092026.run_traces --status
```

`--status` shows the billed cost and the measured output tokens per item. Use them to re-estimate
the full run before asking for its approval. Calibration traces count toward the run.

**4. The full run. One process per model, 15 calls in flight each.**

```powershell
New-Item -ItemType Directory -Force full_run_28092026\traces | Out-Null
$models = "gpt-oss-20b","gemma-4-26b-a4b","deepseek-v4.1-flash","glm-5.3-flash","glm-5.3",
          "muse-glimmer-30b","kimi-k3","gpt-5.4-mini","gemini-3.1-flash-lite","claude-sonnet-5"
foreach ($m in $models) {
  Start-Process python -ArgumentList "-m full_run_28092026.run_traces --model $m --workers 15 --yes" `
    -RedirectStandardOutput "full_run_28092026\traces\$m.log" `
    -RedirectStandardError "full_run_28092026\traces\$m.err" -NoNewWindow
}
```

Each process writes its own `traces/<model>.jsonl`, so the ten cannot collide.

**5. Watch and finish.**

```bash
python -m full_run_28092026.run_traces --status
```

The run is complete when every model shows `missing` at 0. To fill the gaps, re-run the same
command. Only unrun items and service failures are called again, and nothing already answered is
paid for twice.

**6. Save the results locally and to Kaggle.**

Once the run is complete, the trace files are the paid-for result, so keep two copies of them. The
same applies later to the evaluation's score files. Both copies stay private: the traces restate the
pool's questions.

First, a local archive with its checksum, beside the pool backup:

```powershell
$d = "$env:USERPROFILE\EngTrace_private_backup"; $z = "$d\full_run_traces_$(Get-Date -Format yyyy-MM-dd).zip"
Compress-Archive -Path full_run_28092026\traces\*.jsonl -DestinationPath $z -Force
(Get-FileHash $z -Algorithm SHA256).Hash | Out-File -Encoding ascii "$z.sha256"
```

Then a private Kaggle dataset. Switch Claude Code to manual mode for this step, because auto mode
blocks the upload as data exfiltration. Stage it in a short path, because the Kaggle CLI fails on long
temporary-folder paths. Never pass `--public`.

```powershell
$s = "$d\kaggle_traces"; New-Item -ItemType Directory -Force $s | Out-Null; Copy-Item "$z*" $s
[IO.File]::WriteAllText("$s\dataset-metadata.json",
  '{"title": "engtrace-full-run-traces", "id": "<your-kaggle-user>/engtrace-full-run-traces", "licenses": [{"name": "other"}]}')
cd $s; kaggle datasets create -p .
```

Kaggle unpacks the archive. To check the upload, download it back and compare every file's hash with
the local traces. The comparison must print nothing:

```powershell
kaggle datasets download <your-kaggle-user>/engtrace-full-run-traces -p "$d\check" --unzip
$a = Get-ChildItem "$d\check" -Recurse -Filter *.jsonl | Get-FileHash
$b = Get-ChildItem "$env:USERPROFILE\EngTrace\full_run_28092026\traces\*.jsonl" | Get-FileHash
Compare-Object $a.Hash $b.Hash
```

Once it matches, delete `kaggle_traces` and `check`. When traces or scores change, upload a new
version with `kaggle datasets version -p . -m "<what changed>"`, and check it again.

## What a row means

| State | Meaning | Called again? | Scored |
|---|---|---|---|
| `answered` | Text came back, even if truncated. | No | Yes, by the answer check |
| `empty` | The model returned no text. | No | 0 |
| `service_failure` | The HTTP error, timeout or provider fault persisted through 4 attempts. A reply whose finish reason is `error` counts as a fault since D-148; seven such rows in the main run were kept and are reported as a count. | Yes | Never |

An entry in `models.json` with `"run": false` is inert (D-148): no mode calls it unless `--model`
names it, and the dry run and the status say so. The original `qwen3-235b-a22b` and the set-aside
`qwen3.8-27b` are inert, so a bare variant run cannot bill them.

## If something goes wrong

- **The pool on disk does not match the manifest.** The harness refuses to run. Restore `pool/` from the backup and never edit it by hand.
- **The prompt differs from the pilot's.** The harness refuses. `evaluation/run_inference.py` must not change during the run.
- **Error 402.** The account is out of credit. Top it up and re-run.
- **Error 429.** A provider is rate-limiting. The retries absorb short bursts. If failures persist, re-run that model with fewer workers, such as `--workers 8`.

## Rules for the whole run

- **Nothing changes while traces exist.** No template, pool or evaluator edits. A template change means a re-freeze, and every trace of that template must be run again.
- **Traces stay local.** `traces/` is gitignored, because traces restate the pool's questions.
- **Every paid step needs its own approval**, after its estimate.
