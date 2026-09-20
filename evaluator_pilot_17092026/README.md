# evaluator_pilot_17092026 — the evaluation pilot

Frozen 2026-09-17. Five models answer the same 60 problems, giving **300 traces**
that domain experts then annotate step by step. Those annotations are the ground
truth the six evaluator candidates get scored against, which is the whole point:
until now the framework has been compared with itself.

Not to be confused with [`pilot_new_branches/`](../pilot_new_branches/), which was called
`pilot/` until today and is staging from *building* the civil and industrial
branches.

## Where things stand

| Stage | State | Spend |
|---|---|---|
| 0 · freeze 60 items | done | — |
| 1 · 300 traces, five models | done | $6.53 |
| 1b · 120 traces, two open-weight roster models ([D-083](../docs/re-implementation-sep/DECISIONS.md)) | done | $2.15 |
| 2 · E0, the published framework | done, twice ([E0-F7](FINDINGS.md)) | $12.81 |
| 2 · E0-3J, E0 with its third judge connected | done | $7.82 |
| 2 · E3, deterministic milestones | done | $0.00 |
| 2 · E4, E3 + stated-arithmetic checking | done ([RESULTS_E4](RESULTS_E4.md)) | $0.00 |
| 2 · E1, non-suite judge panel | done ([RESULTS_E1](RESULTS_E1.md)) | $4.38 (+ probe $3) |
| 2 · E5, E3 first, a judge on the residuals | done ([RESULTS_E5](RESULTS_E5.md)) | $0.47 (+ $0.43 judge validation) |
| 2 · E2, open process reward models on GPUs | done ([RESULTS_E2](RESULTS_E2.md)) | $0.00 (70 GPU-min, HiPerGator) |
| 3 · expert annotation of the 300 | package ready ([annotation/](annotation/README.md)): app, guide, scoring; awaiting 15 experts | — |

**Read these first:** [FINDINGS.md](FINDINGS.md) (seven defects in the published
framework, plus the robustness cohort), [RESULTS_E0.md](RESULTS_E0.md),
[RESULTS_E3.md](RESULTS_E3.md), [RESULTS_E4.md](RESULTS_E4.md), [RESULTS_E1.md](RESULTS_E1.md), [RESULTS_E5.md](RESULTS_E5.md), [RESULTS_E2.md](RESULTS_E2.md), and the pilot's decisions D-078 to D-090 in
[DECISIONS.md](../docs/re-implementation-sep/DECISIONS.md). Every number in them is
reproduced by a script in [analysis/](analysis/README.md).

## Layout

```
evaluator_pilot_17092026/
  README.md          this file
  FINDINGS.md        defects in the published framework; robustness cohort results
  RESULTS_E0.md      E0 and E0-3J
  RESULTS_E2.md      E2, the three open PRMs (GPU runs and their smoke log: hpg/)
  annotation/        stage 3: the expert annotation app, guide and scoring
  RESULTS_E3.md      E3, its null baseline, where it disagrees with E0
  freeze.py          cuts and pins the slice; --verify fails if anything drifts
  run_traces.py      models over the slice; resumable, --check first
  verify_traces.py   checks produced traces against the slice; T1-T7 + plants
  run_evaluator.py   one harness for every evaluator: (evaluator, traces) -> scores
  models.json        trace models, cohorts, routes, ceilings, prices - and why
  evaluators/
    e0_tribunal.py     E0: the published framework, imported unmodified
    e0_3j_tribunal.py  E0 with the missing get_model_info supplied (third judge)
    milestones.py      per-instance milestones, derived by rule (D-084)
    e3_milestones.py   E3: milestone coverage, order-free, unit-aware
    arith.py           sympy checker for stated arithmetic, validated on gold AND traces (D-085, D-087)
    e4_arith.py        E4: E3 milestones classed verified / contradicted / stated
  analysis/          the script behind every reported number
  kaggle/            GPU offload for the scorer stack (D-086)
  slice/             FROZEN. Do not edit by hand.
    manifest.jsonl     60 items: question, gold solution, seed, SHA-256
    FREEZE.json        the rule, the audit, the counts, the manifest hash
  traces/            PRODUCED, one JSONL per model. Not committed.
  scores/            PRODUCED, one directory per evaluator. Not committed.
  .venv/             the pinned scorer stack for E0 (D-081). Not committed.
```

`slice/` is committed; `traces/` and `scores/` are not. The freeze is the record
that the annotations were made against a specific set of items; untracked, it could
not be verified later, which is the one thing a freeze is for.

## Order of operations

```bash
python -m evaluator_pilot_17092026.freeze --verify      # free. Confirms the slice has not moved.
python -m evaluator_pilot_17092026.run_traces --check   # ~$0.01. Keys, model ids, live prices.
python -m evaluator_pilot_17092026.run_traces --dry-run # free. The plan and the estimated spend.
python -m evaluator_pilot_17092026.run_traces           # the paid run. Resumable.
python -m evaluator_pilot_17092026.run_traces --status  # what exists, and what it actually cost.
python -m evaluator_pilot_17092026.verify_traces        # free. T1-T7 over the produced traces.
```

`--status` counts rows; `verify_traces` is what says the rows mean anything. Add
`--redo-truncated` to a run to re-call only the rows that ended mid-derivation.

`--check` is not a formality. It is what settles whether OpenAI accepts `gpt-5`
verbatim and whether Google's OpenAI-compatible endpoint accepts
`gemini-3.1-pro-preview`; `models.json` carries an OpenRouter fallback for Gemini
if it does not.

### Stage 2 — the evaluators

E0 and E0-3J need the pinned scorer stack; E3 and E4 run on the system Python.

```bash
PY=evaluator_pilot_17092026/.venv/Scripts/python
$PY -m evaluator_pilot_17092026.run_evaluator e0 --dry-run    # free: Tier 1 real, judges recorded not called
$PY -m evaluator_pilot_17092026.run_evaluator e0 --smoke 5    # a few real judged traces, pennies
$PY -m evaluator_pilot_17092026.run_evaluator e0              # the paid run, resumable
$PY -m evaluator_pilot_17092026.run_evaluator e0 --status
python -m evaluator_pilot_17092026.run_evaluator e3           # deterministic, free
python -m evaluator_pilot_17092026.run_evaluator e3 --cohort robustness --model gemma-4-31b
```

The harness refuses to score unless the freeze verifies and the requested columns
pass T1-T7, seeds E0's wrong-answer sample from (item, model) so every judged
evaluator samples the same traces (D-082), and resumes only when both the trace hash
and the evaluator's config hash match.

**Going faster.** The dry run's Tier 1 can go to a Kaggle GPU
(`kaggle/stage_bundle.py`, then push `kaggle/e0-dryrun`, then
`run_evaluator e0 --import-kaggle DIR`, which refuses unless GPU and CPU agree on
reference rows). Paid runs parallelise by running one process per model column
with `ENGTRACE_TORCH_THREADS=2`; the framework's globals make threads unsafe.

The experts' annotation is a separate step and does not wait on any of this.

## The slice

15 templates, one per branch × level cell, four instances each. Chosen by the rule
in `freeze.py`'s docstring and audited in `FREEZE.json`: within a cell, take the
template whose answer type is least represented so far, breaking ties on
BLAKE2b(master_seed, branch, level, template_id). That preference is deliberate —
it pulls symbolic, vector and multipart answers in instead of leaving a slice of
fifteen scalars, and those are precisely the comparator kinds known to be weak.

| Branch | Level | Template | Answer type | Candidates | Excluded from the cell |
|---|---|---|---|---|---|
| chemical | Easy | `rackett_equation_volume` | scalar | 13 | - |
| chemical | Intermediate | `reynolds_number_flow_regime` | classification | 9 | - |
| chemical | Advanced | `gas_phase_concentration` | multipart | 8 | levenspiel_plot_interpretation |
| civil | Easy | `manning_rectangular_discharge` | scalar | 12 | - |
| civil | Intermediate | `cantilever_double_integration` | scalar | 12 | - |
| civil | Advanced | `normal_depth_iteration` | array | 6 | - |
| electrical | Easy | `lorentz_force` | vector | 11 | system_properties_memory_causality |
| electrical | Intermediate | `phasor_addition` | symbolic | 12 | system_property_linearity |
| electrical | Advanced | `coaxial_capacitance` | multipart | 7 | - |
| industrial | Easy | `two_state_steady_state` | vector | 12 | - |
| industrial | Intermediate | `server_configuration_selection` | multipart | 12 | - |
| industrial | Advanced | `aoq_ati_rectifying` | multipart | 6 | line_balancing_heuristic |
| mechanical | Easy | `fluid_particle_acceleration` | vector | 10 | - |
| mechanical | Intermediate | `damping_classification` | classification | 13 | - |
| mechanical | Advanced | `incompressible_continuity` | symbolic | 7 | - |

**All five branches, evenly.** 12 items each, 20 per difficulty level. Answer types:
multipart 16, scalar 12, vector 12, classification 8, symbolic 8, array 4.

Four templates are excluded everywhere, each for a measured reason rather than a
hunch: `system_property_linearity` and `system_properties_memory_causality` are
100% predictable from the question surface (D-057), `line_balancing_heuristic` is
about 90% shortcuttable by counting task pairs (D-046), and
`levenspiel_plot_interpretation` has a blind-guess floor of 1.0000 (D-066). A trace
whose answer is reachable without the reasoning teaches a step-level annotation
nothing.

One instance was dropped: `incompressible_continuity` emits identical text at
indices 1 and 2, so index 2 is skipped and indices 0, 1, 3, 4 are used. A duplicate
costs twice — an expert annotates the same problem again, and the item is
double-weighted in every aggregate. The generator reports 38 duplicate questions
across the full pool, so this was not a one-off.

## Why the text is pinned, not the seed

A seed does not pin an item across template versions, and that is measured rather
than feared: Phase 1 moved 34% of its templates' instances, Phase 2 moved a
template's whole pool, and the C3 corrections moved 824 questions across 13
templates. `seed 7` means one thing today and another after the next constants fix.
So the manifest carries the question and gold text plus a SHA-256 per item, and
`--verify` reports `TEXT MOVED` per item rather than silently relabelling what the
experts read.

The slice was frozen at commit `d430246` against master seed `20260912`; manifest
SHA-256 begins `925b6da7d72c77f6`. `run_traces.py` refuses to make a single call
unless `--verify` passes.

## The five models

| Key | Route | Model id | Role |
|---|---|---|---|
| `gpt-5` | OpenRouter | `openai/gpt-5` | judge family |
| `claude-opus-4.7` | OpenRouter | `anthropic/claude-opus-4.7` | judge family |
| `gemini-3.1-pro` | Google | `gemini-3.1-pro-preview` | judge family |
| `deepseek-r1` | OpenRouter | `deepseek/deepseek-r1-0528` | non-judge control |
| `llama-3.1-70b` | OpenRouter | `meta-llama/llama-3.1-70b-instruct` | non-judge control |

Two axes on purpose: capability, and whether the model's family also sits on the
tribunal. Without the two controls there is no way to tell a judge's competence
from a judge's preference for its own family's output.

All five go through the one `openai` SDK at three base URLs, so the keys stay
scoped exactly as they are in `.env`: `OPENAI_API_KEY` for OpenAI,
`GEMINI_API_KEY` for Google's compatible endpoint, `OPENROUTER_API_KEY` for the
rest. `ANTHROPIC_API_KEY` is not used, and neither is `OPENAI_API_KEY`: both
return 401, so Claude and GPT-5 are both reached through OpenRouter. `openai/gpt-5`
is priced there exactly as the direct API is, so nothing is lost by it. The direct
route stays in `models.json` under `fallback`, to switch back and confirm with
`--check` when a working key exists.

## What the runner records, and why

One row per (item, model), written as it completes:

- the item's SHA-256 from the manifest, so a trace cannot drift from the frozen text
- the prompt's SHA-256, read live out of `evaluation/run_inference.py` rather than
  copied, so a pilot trace and an archive trace were produced by the same
  instructions or the run stops
- the model id as configured **and** the id the API says it served, because a
  `-preview` alias can repoint underneath you
- token usage including reasoning tokens, computed cost, attempts, finish reason

Three behaviours are deliberate:

**Resume, never restart.** `evaluation/run_inference.py` opens its output with
`"w"`, so a rate-limit stall loses everything before it. Here each completed
(item, model) is skipped on a re-run, so an interruption costs only the calls it
had not yet made.

**An empty completion is a failure, even at HTTP 200.** A reasoning model that
spends its whole budget thinking returns `content=""` with tokens billed. Counted
as a success that is a trace nobody can annotate and a row that silently lowers
every score. This repository's recurring defect is a check that is green while
measuring nothing.

**The token ceiling is per model.** Reasoning is billed inside the completion
budget, so the three reasoning models get 8192 and Llama gets 2048. OpenAI's
reasoning models also reject `max_tokens` and want `max_completion_tokens`, while
the other two endpoints want `max_tokens`; the route picks, and a 400 naming the
other one switches.

## Stage 1 result, 2026-09-17

**300 of 300 traces. `verify_traces.py` passes T1-T7, all eight plants fire.**

| Model | served | traces | reported cost | median trace | median completion tokens |
|---|---|---|---|---|---|
| gpt-5 | `openai/gpt-5` | 60 | $2.96 | 1,449 chars | 4,491 |
| claude-opus-4.7 | `anthropic/claude-opus-4.7` | 60 | $1.47 | 1,321 chars | 833 |
| gemini-3.1-pro | `gemini-3.1-pro-preview` | 60 | $0.77 (floor) | 2,376 chars | 1,047 |
| deepseek-r1 | `deepseek/deepseek-r1-0528` | 60 | $1.31 | 3,105 chars | 9,658 |
| llama-3.1-70b | `meta-llama/llama-3.1-70b-instruct` | 60 | $0.02 | 1,406 chars | 583 |
| **total** | | **300** | **$6.53** | 1,646 chars | |

Getting there took three corrections, each one found by a check rather than by
reading output.

**DeepSeek could not fit its own answer.** `deepseek/deepseek-r1` has a single
provider whose hard output cap is 16,000 tokens, and R1 writes 30,000-34,000
characters of reasoning on the iterative items. 28 of 60 truncated mid-reasoning,
and raising `max_tokens` did nothing because the number never reached the model.
Moved to `deepseek/deepseek-r1-0528`: same model, updated checkpoint, four
providers, none capped below 32,000, and cheaper. All 60 re-run rather than the 28
that failed, because the 32 survivors were a different checkpoint and the items
they failed on were the hard ones - topping them up would have left the column
standing on its easy half. The old rows are in `traces/superseded/`.

**19 answers stopped mid-derivation.** 15 Gemini, 4 GPT-5, at
`finish_reason=length`. They had text, so the empty-completion guard passed them
and `--status` counted them as answers; one ended mid-LaTeX-fraction. A quarter of
the Gemini column. T7 exists because of this, and it is a failure now rather than
a footnote.

**Google does not report thinking tokens.** The 15 Gemini truncations reported
323-1,455 completion tokens against an 8,192 ceiling, which cannot be that ceiling
being hit. `reasoning_tokens` is `None` on every Gemini row: the OpenAI-compatible
endpoint omits thinking from `usage`, so thinking was consuming the budget
invisibly and cutting the visible answer. Hence the ceilings at 32,768 against
provider maxima of 65,536 and 128,000 - ours was the binding constraint, not
theirs.

## Cost, and what the first pass corrected

$6.53 against a $1.36 estimate, 4.8x. The whole gap is reasoning tokens, which the
estimator does not model because it is calibrated on archive completions from
non-reasoning models. GPT-5 alone is $2.96 against $0.26 predicted: its median
completion is 4,491 tokens of which only about 362 are the visible answer.

**Gemini's figure is a floor, not a bill.** Google bills thinking tokens that its
OpenAI-compatible endpoint declines to report, so anything computed from returned
usage understates it. It is the only model here with that gap, and any budget
built on these numbers should carry it.

**Stage 2 is not affected, and it is worth being precise about why.** Judges read
the visible trace, not the reasoning. Median visible trace across all 300 is 411
tokens, against the 401 the archive predicted - so the trace-length assumption was
right all along, and only the billing assumption was wrong. E0 over these 300
traces is three judges x 300 calls at about 2,170 input tokens each, roughly
**$10.31**. The pilot phase envelope of $70-150 holds.
