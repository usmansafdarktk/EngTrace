# E1 — E0's Tribunal, judged by families outside the evaluated suite

Run 2026-09-19. **300 of 300 scored, all 178 judged traces decided by all three
judges, 0 judge failures, 0 truncations, $4.38.** Evaluator
`evaluators/e1_panel.py`; every figure below is printed by `analysis/e1_analysis.py`.

The panel (D-088): **Grok 4.6** (xAI), **MiniMax M3**, **MiMo-V2.5-Pro** (Xiaomi) —
none from a family in the paper's 27-model suite. The framework's prompt, parsing,
majority vote and recovery run unmodified; only the judges change.

## Result 1: swapping the judges barely moves the scores

| Trace model | E0 (published judges) | **E1 (non-suite panel)** |
|---|---|---|
| gpt-5 | 0.474 | 0.470 |
| gemini-3.1-pro | 0.454 | 0.452 |
| claude-opus-4.7 | 0.449 | 0.448 |
| deepseek-r1 | 0.402 | 0.401 |
| llama-3.1-70b | 0.115 | 0.117 |
| **all 300** | **0.379** | **0.378** |

E0 and E1 share the same wrong-answer sample (D-082), so they judged exactly the same
178 traces — this is a like-for-like comparison. Against E0-3J (the original three
judges, all connected) on the 154 traces both judged, 3 changed, all slightly down.

**What this says about the judge/judged objection (yAYU #1).** The concern is that
GPT-5 and Gemini inflate their own traces when judging them. Under E0, GPT-5's traces
are judged partly by GPT-5 and score 0.474; judged entirely by models from other
families they score 0.470. Gemini's go 0.454 to 0.452. **No aggregate self-preference
is visible.** This is not yet X2 — a proper test is a per-judge, per-family model
against the expert labels — but it is the first direct evidence, and it points the
same way: removing the overlap does not change the scores.

## Result 2: the new panel agrees with itself less than the original did

Inter-judge agreement, which reviewer 9W1B asked for and E0 never reported. Each
judge's verdict per (trace, step) is parsed from its raw reply with the framework's
own logic; κ corrects for the ~93% base rate of "Alternative Correct".

| Panel | Steps | Fleiss κ, 4 categories | Fleiss κ, correct vs error |
|---|---|---|---|
| **E0-3J** — GPT-5, Opus 4.5, Gemini 3.1 Pro | 1,174 | **0.725** | **0.776** |
| **E1** — Grok 4.6, MiniMax M3, MiMo-V2.5-Pro | 1,165 | **0.563** | **0.622** |

| E1 pair | κ, 4-way | κ, binary |
|---|---|---|
| MiniMax / Grok | 0.567 | 0.637 |
| MiniMax / MiMo | 0.564 | 0.642 |
| Grok / MiMo | 0.561 | 0.590 |

| E0-3J pair | κ, 4-way | κ, binary |
|---|---|---|
| Opus 4.5 / GPT-5 | 0.756 | 0.817 |
| Opus 4.5 / Gemini | 0.730 | 0.761 |
| Gemini / GPT-5 | 0.691 | 0.751 |

The original panel's agreement is **substantial**; the new panel's is **moderate**.
Yet the two produce almost identical scores, because with a 93% base rate and a
majority vote, most disagreements do not change the outcome.

**Read this carefully: agreement is not accuracy.** Higher agreement can mean more
reliable judges, or judges that share training lineage and therefore share errors.
The original three are frontier models from labs that also appear, via distillation,
in each other's orbit; the new three were chosen to spread that. Which panel is
closer to the truth is exactly what the expert labels (X1) will say.

## What the run exposed

**A published-framework robustness defect.** The framework iterates a judge's
`results` calling `.get()` on each item. A reply whose `results` is a string — valid
JSON, wrong shape — does not just lose that judge's vote: it raises and crashes the
whole trace. E0's original judges never produced such replies, so it never showed.
In the first E1 replay it crashed 6 traces (FINDINGS E0-F8).

**One bad provider.** All 11 malformed replies came from one of MiniMax's 13
OpenRouter providers, ModelRun (11 of its 17 replies), against 0 of 161 from the rest.
ModelRun is excluded for MiniMax (D8) and all 17 of its replies were re-fetched.

## Deviations, on every row

| | |
|---|---|
| D2 | all three judges through OpenRouter |
| D3 | wrong-answer sample seeded from (item, model), shared with E0 |
| D4 | `get_model_info` supplied so the third slot is really called |
| D6 | uniform call settings in every slot — JSON mode, temperature 0, 16,384 tokens — replacing per-slot settings tuned to E0's judges (the 2,048 cap would have truncated all three) |
| D7 | empty or malformed replies re-requested, up to 3 attempts |
| D8 | provider ModelRun excluded for MiniMax |

## How it ran: ~20 minutes instead of ~10 hours

The framework calls its judges one after another, and the laptop had room for one
scorer process (3.6 GB free). So E1 ran in three places, and the API key never left
the laptop:

| Step | Where | Time |
|---|---|---|
| capture every judge prompt (framework path, judges stubbed) | Kaggle T4 | 153 s |
| fetch all 534 replies, 48 at a time (HTTP only) | laptop | 11.6 min |
| replay: the real E1 run, replies served from the store | Kaggle T4 | 111 s |

Replies are keyed by (model, settings, prompt), so a reply can never be served to the
wrong prompt, and a missing one is an error on its row, never a silently dropped
judge. The replay's config hash (`f49fa5304dcf`) matches the laptop's exactly.
