# Choosing the judges for E1 and E5

2026-09-19. Evidence: `analysis/judge_candidates.py` (catalogue screen),
`analysis/judge_probe.py` (discrimination probe, ~$3), and the public sources
linked below. **Recommendation, pending sign-off — D-088 is OPEN.** *(D-088 was decided the
same day, and E1 and E5 ran with these judges.)*

**Corrected 2026-09-27 (D-112).** A second step selected as clean, Claude's
`aoq_ati_rectifying#3`, carries an arithmetic slip that all three experts found (0.93^49 is
0.0285538, not 0.02857, so P(X=1) is 0.0999, not 0.1000). It is relabelled in
`judge_probe.RELABEL`, and the probe below is re-scored: 13 slips and 8 clean steps. Before the
relabel MiniMax and MiMo led at 0.96; after it Kimi K3 and Grok 4.6 do, with MiMo, MiniMax,
GPT-5 and GLM-5.3 at 0.92. The recommendation does not change, and the section "Does the
relabel change the choice?" below says why.

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
model** — 13 containing an arithmetic error (found by E4's sympy checker and read, or by
the experts), 8 clean. Each prompt is the framework's own Tribunal prompt, built by its own
code, with one step under review; replies are parsed with the framework's own logic. The
table is `judge_probe.py report` after both relabels (the second, D-112, on 2026-09-27; the
dollar column is the live price at the time of the probe).

| Judge | Family | Parsed | Errors caught | Clean flagged | **Balanced** | $/call | Median s |
|---|---|---|---|---|---|---|---|
| kimi-k3 | Moonshot | 21/21 | 12/13 | 0/8 | **0.96** | 0.0200 | 19.3 |
| grok-4.6 | xAI | 21/21 | 12/13 | 0/8 | **0.96** | 0.0141 | 26.4 |
| mimo-v2.5-pro | Xiaomi | 21/21 | 11/13 | 0/8 | 0.92 | 0.0034 | 45.9 |
| minimax-m3 | MiniMax | 21/21 | 11/13 | 0/8 | 0.92 | **0.0021** | **5.5** |
| *gpt-5* | *E0 judge* | 21/21 | 11/13 | 0/8 | 0.92 | 0.0124 | 13.9 |
| glm-5.3 | Z.ai | 20/21 | 11/13 | 0/8 | 0.92 | 0.0069 | 7.5 |
| *opus-4.5* | *E0 judge* | 21/21 | 10/13 | 0/8 | 0.88 | 0.0198 | 8.7 |
| nemotron-ultra | NVIDIA | **12/21** | 6/13 | 0/8 | 0.73 | 0.0083 | 15.1 |
| seed-2.1 | ByteDance | **11/21** | 5/13 | 0/8 | 0.69 | 0.0154 | 106.2 |

As first scored, with Claude's `aoq_ati_rectifying#3` counted as clean, MiniMax and MiMo led
at 0.96 (11/12, 0/9) and the five judges that flagged that step were charged a false alarm.

**How far to trust it.** With 21 items, one item moves "balanced" by about 0.05, and the
relabel moved the order by exactly one item. The seven judges that parse 20 or 21 of the 21
replies (0.88-0.96) are **not separable** by this probe. What it does settle:

- **Nemotron 3 Ultra and Seed 2.1 are out.** Nearly half their verdicts could not be
  parsed; Seed also truncated 10 of 21 replies at an 8,192-token ceiling, taking a
  median 106 s per step. A judge the framework cannot read contributes nothing.
- **No candidate rubber-stamps.** Every judge that parses caught at least 10 of 13
  known errors.

## Two things the probe showed that are bigger than the choice

**1. E0's judges are not rubber stamps either.** Asked about one known-wrong step,
GPT-5 catches 11 of 13 and Opus 4.5 10 of 13. So E0's 93% "Alternative Correct" (RESULTS_E0)
is better explained by *what gets sent to them*: Tier 1's strict AND rejects almost
every step, including correct ones (E0-F5), so the judges mostly see correct steps
and mostly say so. That is a fairer reading than the "rubber stamp" framing earlier
in this pilot, and FINDINGS E0-F5 is updated to say so.

**2. The deterministic evaluators have a blind spot the judges do not.** One step
selected as clean — Llama computing √(124343 × 582.96) as √72,311,151 (the product is
72,486,995) — is a real error of **0.24%**. It sits inside E4's 1% arithmetic
tolerance and E3's 0.5% milestone tolerance, so both passed it. **Every judge flagged
it.** The label was corrected (`judge_probe.RELABEL`, with the arithmetic). A second
step selected as clean, Claude's `aoq_ati_rectifying#3`, writes 0.93^49 as 0.02857 (it is
0.0285538) inside a chained equality the checker does not compare; all three experts mark
it incorrect, GPT-5, GLM-5.3, Kimi K3, Grok 4.6 and Nemotron flagged it, and MiMo, MiniMax
and Opus did not (relabelled 2026-09-27, D-112). In the other direction, GPT-5's slip in
`normal_depth_iteration#1` — a wrong sign in the written formula, and a value off in its
fourth decimal (1.86921 where the formula gives 1.86942) — got past seven of the eight
judges that returned a verdict, both E0 judges included; only Nemotron flagged it, and E4
catches it deterministically. Each mechanism misses what the other catches. That is the
case for E5.

## Recommendation

### E1: a three-judge panel from three families

| Judge | Why | Exposure it carries |
|---|---|---|
| **MiniMax M3** | 0.92 on the probe, one item below the best (joint-best at 0.96 as first scored); cheapest, fastest, open weights, 21/21 parsed | Claude, OpenAI |
| **Xiaomi MiMo-V2.5-Pro** | 0.92 on the probe, one item below the best (joint-best as first scored); cheap, open weights, the **smallest** documented exposure | Claude (0.4M) |
| **xAI Grok 4.6** | joint-best on the probe after the D-112 relabel (0.96); a non-Chinese lab whose documented exposure is **OpenAI, not Claude** | OpenAI |

Three judges keeps the framework's majority vote meaningful (a 2-1 split carries),
and E1 exists partly to report inter-judge agreement, which needs at least two. The
panel spreads its residual exposure across two evaluated families rather than
stacking three Claude-distilled judges — and X2 can then measure whether any judge
favours the family it is exposed to, because traces from all five generator families
exist.

**Cost:** about $0.020 per judged trace for the panel, so ~$3.50 over the ~178 traces
E0 sends to judges.

**Caveat on Grok:** closed weights, so it cannot be pinned the way the open two can.
If reproducibility outranks lineage spread for you, swap it for **GLM-5.3** (0.92 after
the relabel, open, $0.007) — but see the roster point below: GLM-5.3 is on the roster
D-110 fixed, so it can no longer judge.

### E5: MiMo-V2.5-Pro alone, on the residuals

*Revised 2026-09-19. The first version recommended MiniMax M3 here, on cost and
speed. On defensibility that was the wrong basis; see "How defensible" below.*

E5 runs E3 first and asks a judge only about milestones E3 could not match, so volume
is small and one judge suffices. A single judge carries its exposure alone, so the
deciding criterion is independence, not price: MiMo-V2.5-Pro matched MiniMax exactly
on the probe (11/12 and 0/9 as first scored, 11/13 and 0/8 after the D-112 relabel, 21/21
parsed), is open-weights, and has the **smallest**
documented exposure of any candidate (0.4M Claude exchanges, against MiniMax's 13M+).
MiniMax's advantages — $0.002 vs $0.003 a call, 5.5 s vs 46 s median — are irrelevant
when all of E5 costs under a dollar and runs in parallel.

## How defensible is this in the paper?

**The probe cannot carry a claim, and should not be asked to.** After the relabel, MiMo's
11/13 caught has a 95% Wilson interval of [0.58, 0.96], the same as GPT-5's, and the
leaders' 12/13 is [0.67, 0.99]. "0/8 false flags" is compatible with a true rate up to
0.32. And 10 of the 13 errors are Llama's — the obvious kind. Of the two frontier-model
errors, GPT-5's in `normal_depth_iteration#1` got past nearly every judge, MiniMax and MiMo
included, and Claude's in `aoq_ati_rectifying#3` got past MiMo, MiniMax and Opus. The probe
screened out two unusable judges; it did not validate the rest.

## Does the relabel change the choice? (2026-09-27)

No. The probe was never the deciding criterion, and it still cannot separate the judges
that parse. What decided E5's judge was independence, and the constraints on that have only
tightened:

- **Kimi K3 and GLM-5.3**, which the relabel lifts, are on the roster D-110 fixed, so they
  cannot judge.
- **Grok 4.6**, now joint-best on the probe, has closed weights, so it cannot be pinned, and
  its documented exposure is to OpenAI, whose gpt-5.4-mini and gpt-oss-20b are on the
  roster. MiMo's documented exposure, the smallest of any candidate, is to Claude, and
  claude-sonnet-5 is on the roster too: no candidate is free of exposure, and MiMo's remains
  the least.
- **The larger, matched test run since** — 120 planted defects, each judged planted and
  untouched (RESULTS_X1 Finding 8) — has MiMo detecting about as much as GPT-5 (16 of 51
  conceptual defects against 20 of 60; 43 of 60 arithmetic each) with no false alarm. Grok
  was not run on it. If a reviewer presses the point, that is the run to buy (240 calls,
  about $3.40 at Grok's probe rate, which needs approval), not a re-reading of 21 steps.

What the relabel does change is a description: MiMo and MiniMax are no longer the probe's
best. And the slip MiMo missed is arithmetic inside a checkable step, which the stack assigns
to the digit rule rather than to the judge; D-112 records that the digit rule and the router
miss this step too.

**The most attackable point is MiniMax's independence.** Anthropic named it the
largest distiller of Claude (13M+ of 16M exchanges) and as running a proxy to Anthropic
and OpenAI models. The reply a reviewer raising judge/judged overlap (yAYU #1) will
reach for is: *a Claude judge was replaced by a model trained on Claude's outputs,
judging Claude's traces.* That is why MiniMax is one vote of three in E1 and is not
E5's sole judge.

**What makes the choice defensible — to do, not to argue:**

1. **Panel, not single judge** (E1): exposure diluted, and a 2-1 vote means no single
   judge decides.
2. **Disclose the exposure table** in the paper.
3. **Validate against the expert labels** (X1), with confidence intervals. That, not
   the selection procedure, is what makes a judge defensible.
4. **Measure the bias** (X2) *(done 2026-10-03: RESULTS_LOJO.md, D-174)*: whether each judge favours the family it is exposed to,
   against the expert labels. The pilot has traces from all five generator families;
   a null result settles the objection as well as a positive one.
5. **Show the conclusion survives a swap**: re-run E1 with a different judge in each
   seat (~$3 each) and report whether the ranking of evaluators changes.

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
or xAI models, the recommendation above stands. *(Fixed 2026-09-27, D-110: the roster
includes none of the three, and it includes both Kimi K3 and GLM-5.3.)*
