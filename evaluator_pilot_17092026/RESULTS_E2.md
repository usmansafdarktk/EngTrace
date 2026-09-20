# E2 — an open process reward model scores every step

Run 2026-09-19 on HiPerGator. **360 of 360 traces** (gold 300 + gemma 60), **three
PRMs**, **no API calls, $0**. Traces over a PRM's context limit are not scored:
2 for each Qwen PRM, both DeepSeek-R1, and 0 for VersaPRM. There were 0 step-count
mismatches. The GPU runs are `hpg/` (setup, smoke log, diagnosis: [hpg/README.md](hpg/README.md)).
The harness evaluator is `evaluators/e2_prm.py`, and every figure is printed by
`analysis/e2_analysis.py`. The PRM outputs are committed in `hpg/results/` and the
job logs in `hpg/logs/`.

## How it works

Each PRM reads the question and the trace's steps, split exactly as E0 splits them.
It gives each step the probability that the step is correct, using its model card's
own prompt format and reward extraction. A step is **flagged** when that probability
is below 0.5.

| PRM | Role | Run as |
|---|---|---|
| **Qwen2.5-Math-PRM-72B** | primary | 2 GPUs, eager attention, transformers 4.46.3 |
| VersaPRM (LoRA on Llama-PRM800K) | the multi-domain PRM | 1 GPU, transformers 4.57.3 + peft |
| Qwen2.5-Math-PRM-7B | size cross-check | 1 GPU, eager attention |

Per trace, E2 records:
- `frac_ok`: the share of steps not flagged. This is **the E2 score**, taken from the 72B.
- `min`: the lowest step reward.
- `clean`: 1 when no step is flagged.

Every repo is pinned to a commit, VersaPRM's base included. **Every job first scores its
card's example and refuses to continue unless the check passes.**
- The **7B** reproduces its card to within 0.003.
- The **72B** keeps every per-step verdict, but misses two borderline rewards by up
  to 0.21. That gap was traced to attention-kernel numerics on this GPU, not to the
  pipeline (D-090).
- **VersaPRM's** card prints no reference numbers. Its check proves the adapter is
  really applied.

## Validation: do they catch known errors?

Tested on the judge probe's steps, whose status is known without asking any model
(D-087):
- **SLIP** (12 steps): the shown arithmetic is wrong, confirmed by reading.
- **CLEAN** (9 steps): a step from a right-answer trace whose checkable claims all hold.

| PRM | mean reward, SLIP | mean reward, CLEAN | slips flagged | clean passed | AUROC |
|---|---|---|---|---|---|
| **Qwen-72B** | **0.251** | 0.894 | **9 / 12** | 8 / 9 | **0.935** |
| Qwen-7B | 0.448 | 0.841 | 7 / 12 | 8 / 9 | 0.824 |
| VersaPRM | 0.831 | 0.896 | **1 / 12** | 9 / 9 | 0.565 |

- **VersaPRM does not discriminate.** It passed 11 of 12 known arithmetic slips and
  rates 95–99% of all steps correct for every model. It barely correlates with anything
  (below). As an error detector on this benchmark it **fails validation**, and it is
  reported only as that finding. Its training data (MMLU-Pro reasoning across domains)
  seems not to have taught it to check arithmetic. Its base model with the adapter
  switched off scored the card's wrong step lower (0.59 vs 0.93), so the adapter made
  it more lenient.
- **The 72B separates the two sets well, and the 7B less so.** Scale helps.
- **Caveat on the 72B's AUROC: the probe is confounded by model.** 10 of the 12 slips
  are Llama's, and 8 of the 9 clean steps come from frontier models. The 72B caught 9
  of the 10 Llama slips but **neither non-Llama slip**: 0.969 on Gemma's
  `aoq_ati_rectifying#0`, 0.739 on GPT-5's `normal_depth_iteration#1`. So part of that
  0.935 may be "Llama's steps look worse", not "a wrong step is caught". The one clean
  Llama step scored 0.992, which argues against a pure model effect, but one step is
  not evidence. A clean test needs slips and clean steps from the same models, which
  is what the expert annotation will give.

## Results

| Model | Qwen-72B `frac_ok` | 72B `min` | 72B clean traces | Qwen-7B `frac_ok` | VersaPRM `frac_ok` | E5 strict | E0 F1 |
|---|---|---|---|---|---|---|---|
| gpt-5 | 0.821 | 0.396 | 33% | 0.845 | 0.991 | 0.962 | 0.474 |
| claude-opus-4.7 | 0.864 | 0.431 | 42% | 0.852 | 0.975 | 0.946 | 0.449 |
| gemini-3.1-pro | 0.876 | 0.455 | 48% | 0.809 | 0.969 | 0.929 | 0.454 |
| deepseek-r1 | **0.921** | 0.542 | 57% | 0.850 | 0.969 | 0.923 | 0.402 |
| llama-3.1-70b | **0.560** | 0.122 | 12% | 0.668 | 0.952 | 0.321 | 0.115 |
| *gemma-4-31b* | *0.865* | *0.433* | *43%* | *0.823* | *0.983* | — | — |

What the table shows:

1. **The 72B separates Llama from the rest**, like every other evaluator: 0.560
   against 0.82–0.92.
2. **It reorders the frontier.** DeepSeek-R1 ranks first and GPT-5 last, the reverse
   of E0 and E5. GPT-5 writes terse steps. The one GPT-5 clean probe step the 72B
   flagged (0.243, `coaxial_capacitance#2`) is an example of what this may cost.
   Whether this is a style effect or real error-finding is exactly what the expert
   labels must decide.
3. **It flags much of the frontier.** Only 33–57% of frontier traces pass with no step
   flagged, although most of them reach the right answer. Each flag is either an error
   the other evaluators miss or a false alarm from a math PRM reading engineering
   derivations (units, table look-ups, design choices). E2 cannot tell which, and
   neither can we without labels.

## Agreement between PRMs

| | step verdicts agree | Cohen's κ | trace `frac_ok` Spearman |
|---|---|---|---|
| Qwen-72B vs Qwen-7B | 0.870 | **0.532** | 0.547 |
| Qwen-72B vs VersaPRM | 0.830 | 0.074 | 0.106 |
| Qwen-7B vs VersaPRM | 0.843 | 0.155 | 0.116 |

The two Qwen sizes agree moderately. VersaPRM agrees with neither beyond chance: its
raw agreement of 83–84% is the base rate of calling almost everything correct.

## Against the other evaluators (gold five, per trace)

| | E0 final-answer correct | E3 coverage | E5 strict |
|---|---|---|---|
| Qwen-72B `frac_ok` | 0.467 | 0.477 | **0.550** |
| Qwen-72B `min` | 0.380 | 0.302 | 0.361 |
| Qwen-7B `frac_ok` | 0.251 | 0.279 | 0.379 |
| VersaPRM `frac_ok` | 0.021 | 0.106 | 0.160 |

The 72B lines up best with E5-strict. That is some independent support: E5 settles
76% of its milestones by deterministic number-matching and asks an LLM judge
(MiMo-V2.5-Pro) only about the rest, so the two share no machinery. The caveat in RESULTS_E5
applies unchanged. None of these correlations ranks an evaluator, and E0's
correctness check is flawed (E0-F1, E0-F2). Only the expert annotation can referee.

## Cost and run

$0 in API spend. GPU time on HiPerGator's RTX PRO 6000 Blackwell:
- the 72B: 6 min on 2 GPUs
- VersaPRM: 1.3 min on 1 GPU
- the 7B: 1.2 min on 1 GPU

The three full runs took 15 GPU-minutes. The smoke attempts, the diagnostic and one
run whose output was lost to the full `/blue` disk took another 55, so 70
GPU-minutes in all, measured from Slurm's accounting (`sacct`). Weights were
fetched to node-local disk and deleted by each job's exit trap. Results were copied
down with md5 checks, and no E2 job is left on the cluster.
