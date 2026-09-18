# Choosing the judges for E1 and E5

2026-09-19. Evidence: `analysis/judge_candidates.py` (catalogue screen),
`analysis/judge_probe.py` (discrimination probe, ~$3), and the public sources
linked below. **Recommendation, pending sign-off — D-088 is OPEN.**

## The constraint

E1 is "the same architecture, judges drawn only from families not in the evaluation
suite". The suite is the paper's Table 13, 27 models, which spans **seven families
once backbones are counted**:

| Family | In the suite as |
|---|---|
| OpenAI | GPT-5, GPT-5 Mini, GPT-4.1, GPT-4.1 Mini |
| Anthropic | Claude Opus 4.7, Sonnet 4.5, Sonnet 4, Sonnet 3.7 |
| Google | Gemini 3.1 Pro, 3 Pro, 2.5 Pro, 2.5 Flash — **and Gemma 3 27B, Gemma 2 9B** |
| DeepSeek | V4 Pro, V3, R1 |
| Meta | Llama 3.1 70B, 8B — **and MetaMath (Llemma, a Llama derivative)** |
| Alibaba / Qwen | Qwen 2.5 72B/14B/7B, Qwen 3 8B, Qwen2.5-Math |
| Mistral | Mathstral — **and WizardMath (Mistral-7B backbone)** |

That rules out every fine-tune of those bases too (Hermes, Magnum, Dolphin,
Euryale, WizardLM-2), and `meta/muse-spark`, which is Meta.

## Independence is a matter of degree

No capable candidate is clean. Each has its own pretrained base, and nearly every
one has **documented exposure to an evaluated family's outputs**:

| Candidate | Own base | Documented exposure |
|---|---|---|
| Moonshot — Kimi K3 | yes | Claude: 23M exchanges, and forwarding its own customers' requests to Claude ([Anthropic, Feb](https://x.com/AnthropicAI/status/2025997928242811253); [Sept](https://thehackernews.com/2026/09/anthropic-says-seven-china-based-ai.html)) |
| MiniMax — M3 | yes | Claude: 13M+ exchanges; a proxy service to Anthropic and OpenAI models ([Sept](https://thehackernews.com/2026/09/anthropic-says-seven-china-based-ai.html)) |
| Z.ai — GLM-5.3 | yes | Claude: 3.4M exchanges ([Sept](https://thehackernews.com/2026/09/anthropic-says-seven-china-based-ai.html)) |
| Xiaomi — MiMo-V2.5-Pro | yes | Claude: 0.4M exchanges, the smallest named ([Sept](https://thehackernews.com/2026/09/anthropic-says-seven-china-based-ai.html)) |
| xAI — Grok 4.6 | yes | OpenAI: Musk testified Grok was "partly" trained on OpenAI models ([TechCrunch](https://techcrunch.com/2026/04/30/elon-musk-testifies-that-xai-trained-grok-on-openai-models/)) |
| NVIDIA — Nemotron 3 Ultra | yes, from scratch | its own report: SFT math data from **DeepSeek-V3.2**, plus GPT-OSS and Qwen ([report](https://arxiv.org/html/2606.15007v1)) |
| ByteDance — Seed 2.1 | yes | none reported; its founder ruled distillation out ([AI Weekly](https://aiweekly.co/alerts/bytedances-zhang-yiming-rules-out-distillation-for-ai-models)) |
| StepFun | yes | named by US agencies ([THN](https://thehackernews.com/2026/09/us-agencies-accuse-china-ai-firms-of.html)) |

The concern is well founded: judges are over 50% more likely to wrongly pass a
rubric item when the output is their own family's, even under objective criteria
([Self-Preference Bias in Rubric-Based Evaluation](https://arxiv.org/abs/2604.06996)).
Distillation exposure is weaker than shared family, but it is the same mechanism
in miniature, and it is why the panel below spreads it across different evaluated
families instead of concentrating it.

## The probe: do they discriminate?

Desk research cannot say whether a judge actually notices a wrong step. E0's judges
answered "Alternative Correct" to 93% of the steps sent to them, so that was tested
directly: **21 steps from real pilot traces whose status is known without asking a
model** — 12 containing an arithmetic error (found by E4's sympy checker, then read),
9 clean. Each prompt is the framework's own Tribunal prompt, built by its own code,
with one step under review; replies are parsed with the framework's own logic.

| Judge | Family | Parsed | Errors caught | Clean flagged | **Balanced** | $/call | Median s |
|---|---|---|---|---|---|---|---|
| minimax-m3 | MiniMax | 21/21 | 11/12 | 0/9 | **0.96** | **0.0021** | **5.5** |
| mimo-v2.5-pro | Xiaomi | 21/21 | 11/12 | 0/9 | **0.96** | 0.0034 | 45.9 |
| *opus-4.5* | *E0 judge* | 21/21 | 10/12 | 0/9 | 0.92 | 0.0198 | 8.7 |
| kimi-k3 | Moonshot | 21/21 | 11/12 | 1/9 | 0.90 | 0.0200 | 19.3 |
| grok-4.6 | xAI | 21/21 | 11/12 | 1/9 | 0.90 | 0.0141 | 26.4 |
| *gpt-5* | *E0 judge* | 21/21 | 10/12 | 1/9 | 0.86 | 0.0124 | 13.9 |
| glm-5.3 | Z.ai | 20/21 | 10/12 | 1/9 | 0.86 | 0.0069 | 7.5 |
| seed-2.1 | ByteDance | **11/21** | 5/12 | 0/9 | 0.71 | 0.0154 | 106.2 |
| nemotron-ultra | NVIDIA | **12/21** | 5/12 | 1/9 | 0.65 | 0.0083 | 15.1 |

**How far to trust it.** With 21 items, one item moves "balanced" by about 0.05. The
top five (0.86-0.96) are **not separable** by this probe. What it does settle:

- **Nemotron 3 Ultra and Seed 2.1 are out.** Nearly half their verdicts could not be
  parsed; Seed also truncated 10 of 21 replies at an 8,192-token ceiling, taking a
  median 106 s per step. A judge the framework cannot read contributes nothing.
- **No candidate rubber-stamps.** Every judge that parses caught at least 10 of 12
  known errors.

## Two things the probe showed that are bigger than the choice

**1. E0's judges are not rubber stamps either.** Asked about one known-wrong step,
GPT-5 and Opus 4.5 each catch 10 of 12. So E0's 93% "Alternative Correct" (RESULTS_E0)
is better explained by *what gets sent to them*: Tier 1's strict AND rejects almost
every step, including correct ones (E0-F5), so the judges mostly see correct steps
and mostly say so. That is a fairer reading than the "rubber stamp" framing earlier
in this pilot, and FINDINGS E0-F5 is updated to say so.

**2. The deterministic evaluators have a blind spot the judges do not.** One step
selected as clean — Llama computing √(124343 × 582.96) as √72,311,151 (the product is
72,486,995) — is a real error of **0.24%**. It sits inside E4's 1% arithmetic
tolerance and E3's 0.5% milestone tolerance, so both passed it. **Every judge flagged
it.** The label was corrected (`judge_probe.RELABEL`, with the arithmetic). In the
other direction, GPT-5's sign slip — right value, wrong written formula — got past
six of the seven judges that returned a verdict; only Nemotron flagged it, and E4 catches it
deterministically. Each mechanism misses what the other catches. That is the case
for E5.

## Recommendation

### E1: a three-judge panel from three families

| Judge | Why | Exposure it carries |
|---|---|---|
| **MiniMax M3** | joint-best discrimination, cheapest, fastest, open weights, 21/21 parsed | Claude, OpenAI |
| **Xiaomi MiMo-V2.5-Pro** | joint-best discrimination, cheap, open weights, the **smallest** documented exposure | Claude (0.4M) |
| **xAI Grok 4.6** | top-tier discrimination; a non-Chinese lab whose documented exposure is **OpenAI, not Claude** | OpenAI |

Three judges keeps the framework's majority vote meaningful (a 2-1 split carries),
and E1 exists partly to report inter-judge agreement, which needs at least two. The
panel spreads its residual exposure across two evaluated families rather than
stacking three Claude-distilled judges — and X2 can then measure whether any judge
favours the family it is exposed to, because traces from all five generator families
exist.

**Cost:** about $0.020 per judged trace for the panel, so ~$3.50 over the ~178 traces
E0 sends to judges.

**Caveat on Grok:** closed weights, so it cannot be pinned the way the open two can.
If reproducibility outranks lineage spread for you, swap it for **GLM-5.3** (0.86,
open, $0.007) — but see the roster point below.

### E5: MiniMax M3 alone, on the residuals

E5 runs E3 first and asks a judge only about milestones E3 could not match, so volume
is small and one judge suffices. MiniMax M3 is joint-best on discrimination, cheapest,
fastest (5.5 s median) and parsed every reply. At most ~$0.75 for the whole pilot.
Second choice: MiMo-V2.5-Pro (same discrimination, eight times slower).

### Deliberately not recommended: Kimi K3 and GLM-5.3

Both discriminate well. They are left out for two reasons:

1. **They are the strongest open-weights models** (Artificial Analysis index 43.8 and
   44.9 — [source](https://benchlm.ai/benchmarks/artificialanalysis)), exactly what the
   next EngTrace roster is likely to *evaluate* if it moves to open models. A family
   cannot be both judge and judged. Choosing them as judges now forecloses them from
   the roster.
2. Kimi carries the heaviest documented exposure of any candidate (23M Claude
   exchanges, plus forwarding customer traffic to Claude).

**The decision that actually needs making first:** fix the next roster's families,
then the judges are whatever is left. If the roster will not include MiniMax, Xiaomi
or xAI models, the recommendation above stands.
