# evaluator_pilot_17092026 — the evaluation pilot

Frozen 2026-09-17. Five models answer the same 60 problems, giving **300 traces**
that human experts then annotate step by step. Those annotations are the ground
truth the six evaluator candidates get scored against, which is the whole point:
until now the framework has been compared with itself.

Not to be confused with [`pilot_new_branches/`](../pilot_new_branches/), which was called
`pilot/` until today and is staging from *building* the civil and industrial
branches.

## Layout

```
evaluator_pilot_17092026/
  README.md        this file
  freeze.py        cuts and pins the slice; --verify fails if anything drifts
  run_traces.py    the five models over the slice; resumable, --check first
  models.json      the five models, their routes, keys, ceilings and prices
  slice/           FROZEN. Do not edit by hand.
    manifest.jsonl   60 items: question, gold solution, seed, SHA-256
    FREEZE.json      the rule, the audit, the counts, the manifest hash
  traces/          PRODUCED by run_traces.py, one JSONL per model. Not committed.
```

`slice/` is committed, `traces/` is not. The freeze is the record that 300 paid
annotations were made against a specific set of items; untracked, it cannot be
verified later, which is the one thing a freeze is for.

## Order of operations

```bash
python -m evaluator_pilot_17092026.freeze --verify      # free. Confirms the slice has not moved.
python -m evaluator_pilot_17092026.run_traces --check   # ~$0.01. Keys, model ids, live prices.
python -m evaluator_pilot_17092026.run_traces --dry-run # free. The plan and the estimated spend.
python -m evaluator_pilot_17092026.run_traces           # the paid run. Resumable.
python -m evaluator_pilot_17092026.run_traces --status  # what exists, and what it actually cost.
```

`--check` is not a formality. It is what settles whether OpenAI accepts `gpt-5`
verbatim and whether Google's OpenAI-compatible endpoint accepts
`gemini-3.1-pro-preview`; `models.json` carries an OpenRouter fallback for Gemini
if it does not.

After the traces exist: the evaluator candidates run over them (E0, the current
three-LLM tribunal, is the baseline), and the experts annotate them. Those are
separate steps. This directory only produces the material.

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

The slice was frozen at commit `0f86bc7` against master seed `20260912`; manifest
SHA-256 begins `925b6da7d72c77f6`. `run_traces.py` refuses to make a single call
unless `--verify` passes.

## The five models

| Key | Route | Model id | Role |
|---|---|---|---|
| `gpt-5` | OpenAI | `gpt-5` | judge family |
| `claude-opus-4.7` | OpenRouter | `anthropic/claude-opus-4.7` | judge family |
| `gemini-3.1-pro` | Google | `gemini-3.1-pro-preview` | judge family |
| `deepseek-r1` | OpenRouter | `deepseek/deepseek-r1` | non-judge control |
| `llama-3.1-70b` | OpenRouter | `meta-llama/llama-3.1-70b-instruct` | non-judge control |

Two axes on purpose: capability, and whether the model's family also sits on the
tribunal. Without the two controls there is no way to tell a judge's competence
from a judge's preference for its own family's output.

All five go through the one `openai` SDK at three base URLs, so the keys stay
scoped exactly as they are in `.env`: `OPENAI_API_KEY` for OpenAI,
`GEMINI_API_KEY` for Google's compatible endpoint, `OPENROUTER_API_KEY` for the
rest. `ANTHROPIC_API_KEY` is not used — Claude is reached through OpenRouter.

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

## Cost

`--check` re-reads the live OpenRouter catalogue, so the estimate is never a stale
snapshot. At archive-mean lengths `--dry-run` puts the 300 traces at **$1.36**
(gpt-5 $0.26, Claude $0.69, Gemini $0.32, DeepSeek $0.07, Llama $0.02), and the 300
traces plus one clean E0 tribunal pass came to **$11.40**, $6.41 of it Claude.
Reasoning tokens push the real figure above the mean-based estimate, which is why
`--status` reports what was actually spent from returned usage. The pilot phase as
a whole was quoted at $70-150, covering the other evaluator candidates, X2 and
retries.
