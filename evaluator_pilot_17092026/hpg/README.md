# E2 on HiPerGator: open process reward models

E2 scores every step of every trace with an open process reward model (PRM), run locally
on HiPerGator GPUs. No API key, token or credential is involved: the models are public
and ungated, and inference happens on the node.

| File | What it does |
|---|---|
| `build_e2_inputs.py` | Writes `e2_inputs.jsonl`: 360 traces (gold 300 + gemma 60), 2,496 steps, split with E0's own `extract_steps` so step k is the same step in every evaluator. |
| `e2_score.py` | Loads one PRM exactly as its model card specifies, runs the self-check, then scores. Records per trace: repo + pinned revision, step probabilities, token count vs context, `over_length`, `step_mismatch`, hashes. |
| `e2.sbatch` | One PRM over the inputs (`5` = smoke, `all` = full). Weights go to node-local `/tmp` in a per-job directory that a trap deletes however the job ends (`/blue` is 99% full). |
| `e2_diag.sbatch` | Self-check only, under each attention implementation; used to diagnose the Qwen-72B card mismatch below. |

Files reach the cluster as base64 heredocs piped through the already-open SSH control
connection (`ssh -O proxy hpg`), and are checked by md5 on both ends.

## PRMs and pins

| PRM | Repo @ revision | transformers |
|---|---|---|
| qwen72 | `Qwen/Qwen2.5-Math-PRM-72B` @ `9df429b0` (2 GPUs) | 4.46.3 |
| qwen7 | `Qwen/Qwen2.5-Math-PRM-7B` @ `06107400` (size cross-check, 1 GPU) | 4.46.3 |
| versa | `UW-Madison-Lee-Lab/VersaPRM` @ `cd537f77`, a LoRA adapter on `UW-Madison-Lee-Lab/Llama-PRM800K` @ `1973a85d` | 4.57.3 + peft 0.17.1 |

Hardware: RTX PRO 6000 Blackwell (97.9 GB), torch 2.13.0+cu130, container
`vllm-openai.sif`. flash-attn is not installed in the container.

## Smoke log (2026-09-19)

| Attempt | Outcome | Cause / fix |
|---|---|---|
| 1 | both failed at install | the container's tokenizers 0.23.2 is outside transformers 4.57.3's range; pinned 0.22.2. `/tmp` turned out to be the node's shared disk, so a per-job directory and a trap cleanup were added. |
| 2 | both failed at start | binding a `/tmp` subdirectory into the container collides with a file the image ships (FATAL mount hook); removed the bind, since Apptainer mounts `/tmp` by default. |
| 3 | versa: base not found; qwen72: crash in forward | versa: loading the adapter through `AutoModelForCausalLM` passed the adapter's revision on to the base repo. The base is now loaded explicitly at its own pinned commit. qwen72: the card's `modeling_qwen2_rm.py` calls a cache API that transformers 4.57 removed; Qwen now runs on 4.46.3. |
| 4 | **versa passed**; qwen72 self-check failed | see below |

**VersaPRM passed.** Its card prints no reference rewards, so the self-check requires one
score per step, each in [0, 1], and that the same example scored with the adapter
switched off gives different scores. It did: the scores shifted by up to 0.34, so the
adapter is really loaded. The 5 smoke traces scored with no over-length traces and no
step mismatches, and the scores were bit-identical across two runs.

**Qwen-72B does not reproduce its card's example rewards** (tolerance 0.02). The
scorer therefore refused to score, as designed. The diagnostic job (`e2_diag.sbatch`,
job 42701952) scored the card example under every available condition:

| Card example (454 tokens) | step 1 | step 2 | step 3 | step 4 |
|---|---|---|---|---|
| Qwen-7B, card | 1.0 | 0.190 | 0.977 | 1.0 |
| Qwen-7B, here (sdpa) | 0.9995 | 0.157 | 0.974 | 0.9997 |
| Qwen-72B, card | 0.992 | 0.0048 | 0.324 | 0.820 |
| Qwen-72B, here (sdpa) | 0.993 | 0.0037 | 0.159 | 0.561 |
| Qwen-72B, here (eager) | 0.993 | 0.0036 | 0.184 | 0.608 |

What this establishes:

- **Not the library version.** transformers 4.57.3 and 4.46.3 gave bit-identical
  rewards. The rotary position buffer was fp32 under 4.46.3, which rules out the
  first hypothesis, a bf16 position table.
- **Not the softmax precision.** Taking the softmax in bf16, as the card does, moves
  nothing by more than 0.002.
- **Not the pipeline.** The Qwen-7B, with the same code, template, step token and
  extraction, lands within 0.033 of its card with the identical pattern.
- **The attention kernel moves the 72B's borderline steps** by about 0.05 (sdpa vs
  eager). The remaining gap to the card is most plausibly the GPU/kernel stack the
  card's numbers were produced on. This cannot be tested directly here: the only
  other GPU partitions are B200 (also Blackwell) and L4 (24 GB, too small for the 72B).
- **The verdicts agree.** At the 0.5 threshold PRMs are used with, every step gets
  the same correct/incorrect verdict as the card, for both models.

**Decision (D-090).** The self-check now passes on per-step verdict agreement at 0.5.
The continuous gap is kept in every run's `.meta.json` (`max_abs_diff`, `within_tol`).
Qwen runs with eager attention, the closest to its card. The 7B joins as a size
cross-check. E2 is reported on thresholded verdicts and on the probe's known-label
steps, not by comparing continuous rewards with published numbers.

**With eager attention the 7B passes even the original 0.02 check** (max gap 0.003,
against 0.033 under sdpa). That confirms the attention kernel as the cause. The 72B
under eager passes on verdicts, with its 0.21 gap recorded (smoke job 42710089).

Every run's `.meta.json` also records the scorer's own md5, the Slurm job id and the
host, so a result file names the exact code that produced it.

## Storage

On 2026-09-19 `/blue/fire-finai` hit 100% (10 TB of 10 TB, of which 9.3 TB is the
group's data, not E2's). The first full VersaPRM run scored all 360 traces, then lost
its output file to `Disk quota exceeded`. Results and logs (a few MB) now go to
`/home/zhuohan.xie/engtrace-e2/{out,logs}`, and the batch scripts live there too.
The scorer and the inputs stay on `/blue` and are only read. Weights stay on node-local
`/tmp`.

## Job names

Jobs carry no model names, since the account is shared: `e2_<n>_<smoke|full>`, where
**1 = versa, 2 = qwen72, 3 = qwen7**, plus `e2_diag`. Logs are at
`logs/<name>_<jobid>.out`. The smoke attempts before this scheme were named
`e2-versa-smoke`, `e2-qwen72-smoke` and `e2-diag` (jobs 42698932-42701952); finished
jobs cannot be renamed.
