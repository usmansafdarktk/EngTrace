SIGNAL: HARNESS READY 2026-10-06 21:40 UTC
SIGNAL: RUNS SCORED 2026-10-07 12:45 UTC (every WS-B store run, scored and staged; results/matched_config.json written)

# WS-B report

Done. The sections below are current as of 2026-10-07 13:30 UTC. Spend: $70.66 of the $90 cap ($39.48 inference,
$30.88 judged stages, $0.30 the Qwen Thinking calibration).

- **Judged stages**, launched 2026-10-07 02:07 UTC after WS-A posted STAGES DONE (02:03 UTC), finished 12:32 UTC. One
  detached chain ran one judge-model process at a time: E5 for the closed models, Gemma and the paraphrase pairs, then
  the router for the closed models and Gemma, then one more pass over every step. The reply stores held $52.0955 (E5)
  and $109.3397 (router) at launch, so the cumulative caps were $61.10, $66.10 and $68.10 (E5) and $131.34 and
  $142.34 (router). No cap was reached. The chain log is `scores/reasoning-medium-full/stages_chain.log`, and each step
  logged to `scores/<variant>/e5.log` or `router.log`.
- Workers were raised on the owner's instruction, 4 to 16 (04:28 UTC) to 24 (04:37) and, for the router, to 32
  (07:16); the main run's stages ran at 16 to 24. Each time the running stage process was stopped, the chain stopped
  itself, the reply store ended cleanly (every line parses) and the relaunched chain resumed from the store, so no
  call was bought twice. No call failed at any setting. Throughput was about 2.6 calls a minute at 4 workers; E5 ran
  at about 17 at 24 workers; the router at about 10 at 24 workers and 18 to 24 at 32.
- Every call of every stage has a reply (`without_reply` 0). Two calls passed the 960 s deadline; both replies
  arrived and are stored.

## Harness
- `full_run_28092026/run_traces.py` (and `models.json`) installed 2026-10-06 21:36 UTC, after WS-A's
  `traces/_round5_rerun/DONE.json` (21:31 UTC). The existing variants' free modes print exactly what the committed
  version prints (status and dry runs of main, paraphrase, the repeats, the reasoning, flagship, openbook and tool arms,
  and `--selftest`, compared output for output against `HEAD` on the same `models.json`).
- Variants added: `reasoning-medium-full` (all 2,250 items, the reasoning parameter at effort medium, the main run's
  prompt, ceiling and routing; models from `MATCHED_MODELS` or `--model`); the cascade `paraphrase-reasoning-medium`
  (the 275 kept pairs: the 314 passing paraphrases less the 39 the expert rejected, the rule of `analyze.accepted_pairs`)
  and `repeat1-reasoning-medium` to `repeat3-reasoning-medium` (the 300 repeat items), restricted to `MATCHED_MODELS`
  (the repeats further to the decoding repeats' set of D-151). `MATCHED_MODELS` = GPT-5.4 mini, Gemini 3.1 Flash-Lite,
  Gemma 4 26B.
- Flags added, all free except where `--yes` is needed anyway:
  - `--only-items PATH` (every mode of every variant). An id that is no pool item stops the run; an id outside the
    variant is skipped and counted.
  - `--review`: missing items, finish reasons, empty rows at the ceiling, reasoning-token and reasoning-text shares,
    request sets, providers, the prompt hash against the main run's rows, and the bill of any rows moved to
    `traces/_unfinished/`. It exits 1 if an item is missing or a prompt differs.
  - `--trace-review`: runs `trace_review.py` on a matched variant, given that variant's item set (see Open items).
  - `--requeue-unfinished`: moves a matched variant's rows that ended without a finish reason to
    `traces/_unfinished/` (see Runs).
  - `--matched-config`: writes `results/matched_config.json`.
- Dry runs of the new variants price from bills measured with the parameter on: the 450-item arm's rows, or the model's
  own calibration rows; the cascade from the model's `reasoning-medium-full` rows on the same items.
- `--check` (billed, $0.0077 in all) passes:

  ```
  --variant reasoning-medium-full --model gpt-5.4-mini --model gemini-3.1-flash-lite --model gemma-4-26b-a4b --only-items <the 30 ids> --check --yes
  --only-items ...: 30 of its 30 ids are items of this variant
    prompt sha256 c2bcb87984c4e50b; the main run's rows carry c2bcb87984c4e50b (26970 rows): the same
    gemma-4-26b-a4b        answered  served=google/gemma-4-26b-a4b-it provider=Io Net billed=$0.00127 reasoning_tokens=4849
    gpt-5.4-mini           answered  served=openai/gpt-5.4-mini provider=Azure billed=$0.00038 reasoning_tokens=52
    gemini-3.1-flash-lite  answered  served=google/gemini-3.1-flash-lite provider=Google AI Studio billed=$0.00101 reasoning_tokens=568
  --variant main --model qwen3-235b-a22b-thinking-2507 --check --yes
    qwen3-235b-a22b-thinking-2507 answered  served=qwen/qwen3-235b-a22b-thinking-2507 provider=Novita billed=$0.00396
  --variant paraphrase-reasoning-medium --model gemini-3.1-flash-lite --check --yes
    gemini-3.1-flash-lite  answered  served=google/gemini-3.1-flash-lite provider=Google AI Studio billed=$0.00105 reasoning_tokens=592
  ```

## Roster check
- Qwen Thinking (`qwen3-235b-a22b-thinking-2507`): eligible.
  - It is on neither exclusion list in `docs/inference_pricing/build_pricing_doc.js`. The pilot generators are gpt-5,
    claude-opus-4.7, gemini-3.1-pro, deepseek-r1, llama-3.1-70b, gemma-4-31b-it and qwen3.8-27b; the judges are gpt-5,
    claude-opus-4.5, gemini-3.1-pro, grok, minimax and mimo.
  - In `evaluator_pilot_17092026/`, ignored files included, its id appears only in `scores/_cache/openrouter_prices.json`,
    the price snapshot `run_evaluator.py` reads (458 models). No trace, judge row or configuration names it.
- Identifiers: OpenRouter `qwen/qwen3-235b-a22b-thinking-2507`; Hugging Face `Qwen/Qwen3-235B-A22B-Thinking-2507`
  (Apache 2.0). Its card says it "supports only thinking mode" and recommends 32,768 output tokens.
- Endpoints, read 2026-10-06:
  - Novita fp8, $0.30/$3.00 per M, output cap 32,768: the one routing picks.
  - Venice fp8, $0.45/$3.50, cap 16,384.
  - Alibaba, $0.23/$2.30, quantization undeclared: not eligible under the fp8-or-better rule.
- `models.json` entry: inert (`run: false`), `added_for: "matched settings"`, `added: 2026-10-06`,
  `sibling_of: qwen3-235b-a22b-2507`. It runs only as `--variant main --model qwen3-235b-a22b-thinking-2507`, so a plain
  resume of the main run by another stream never calls it.
- Price: the dry run on the pricing basis (270 in, 5,232 out per item) gave $35.50 for the 2,250 items, before about
  $13 of judged stages, against the brief's expected $6 to $12 (which matches the Instruct release's price).
- Calibration, 2026-10-07 13:03 UTC, at the owner's request (`--variant main --calibrate 20 --workers 20 --yes`):
  - 20 of 20 answered, all `stop`, all served by Novita; none empty or at the ceiling.
  - Reasoning tokens on 20 of 20 (median 2,530); reasoning text on 16; no `<think>` tag in the answer text.
  - Completion tokens median 3,930, max 13,153. Time per item median 191 s, mean 301 s, max 897 s.
  - Billed $0.298, $0.0149 an item. The full set would bill about $33.60, plus about $10 of judged stages. It would take
    about 3 h of inference at 64 workers, or 6 h at 32, before the stages.
  - The 20 rows are in `traces/_calibration/qwen3-235b-a22b-thinking-2507.jsonl`, out of the main traces folder.
- **Not run, by the owner's decision (2026-10-07): too little time before the deadline.** By the brief it was to run
  as a twelfth row once it passed the roster check, unless the cap bound. When it was priced (2026-10-06 21:15 UTC) it
  did not fit the cap, which WS-B applied without asking. Raised then, the run could have finished overnight.

## Gemma 4 calibration
- Parameter tried: OpenRouter's unified `reasoning: {"effort": "medium"}`, the request the closed models get
  (`--variant reasoning-medium-full --model gemma-4-26b-a4b --calibrate 20 --yes`, 2026-10-06 21:37 UTC). No alternative
  was needed. The model card (Hugging Face `google/gemma-4-26B-A4B-it`) names "configurable thinking modes", switched on
  by `enable_thinking`, and OpenRouter lists the `reasoning` parameter for every Gemma 4 endpoint (read 2026-10-06).
- Rows:
  - 20 of 20 answered, finish `stop` on all 20.
  - Reasoning tokens on 20 of 20 (median 4,494, max 15,134); completion median 5,339, max 16,835.
  - Reasoning text on 20 of 20.
  - Served by Io Net (bf16) on all 20; billed $0.032 ($0.00159 per row).
  - In the main run, at the provider's default, Gemma returned reasoning tokens on 0 of 2,250 rows.
- Decision: Gemma 4 is in. It is in `MATCHED_MODELS`, and its 20 calibration rows count toward its full run, as the
  arm's calibration rows did (D-180).

## Runs
| model | variant | rows | empty at ceiling | reasoning-token share | providers | inference $ | judge calls / $ | step-check calls / $ |
|---|---|---:|---:|---:|---|---:|---|---|
| GPT-5.4 mini | `reasoning-medium-full` | 2,250 | 0 | 1.000 | Azure 2,250 | 22.650 | 831 / 2.442 | 2,250 / 7.988 |
| Gemini 3.1 Flash-Lite | `reasoning-medium-full` | 2,250 | 0 | 1.000 | Google AI Studio 2,132, Google 118 | 5.526 | 768 / 2.364 | 2,250 / 7.666 |
| Gemma 4 26B | `reasoning-medium-full` | 2,250 | 3 | 1.000 | Io Net 2,194, DekaLLM 56 | 3.795 | 792 / 2.479 | 2,247 / 7.110 |
| GPT-5.4 mini | `paraphrase-reasoning-medium` | 275 | 0 | 1.000 | Azure 275 | 2.600 | 97 / 0.270 | not staged |
| Gemini 3.1 Flash-Lite | `paraphrase-reasoning-medium` | 275 | 0 | 1.000 | Google AI Studio 184, Google 91 | 0.660 | 91 / 0.283 | not staged |
| Gemma 4 26B | `paraphrase-reasoning-medium` | 275 | 0 | 1.000 | Io Net 239, DekaLLM 34, NextBit 2 | 0.481 | 93 / 0.280 | not staged |
| Gemini 3.1 Flash-Lite | `repeat1..3-reasoning-medium` | 3 x 300 | 0 | 1.000 | Google AI Studio and Google | 2.196 | not staged | not staged |
| Gemma 4 26B | `repeat1..3-reasoning-medium` | 3 x 300 | 1 | 1.000 | Io Net, DekaLLM, NextBit | 1.566 | not staged | not staged |

- Inference total $39.48: $39.47 in the rows (final rows, plus $0.044 in the 26 rows moved to `traces/_unfinished/`),
  plus $0.008 of `--check` calls. Gemma's figures include those moved rows and its $0.032 calibration. Estimates
  were $27.91 for the two closed models and $3.55 for Gemma. GPT-5.4 mini billed $22.65 against $22.36.
- Every store: `run_traces --review` exits 0, with missing 0, one request set, every row on the main run's prompt
  (`c2bcb87984c4e50b`) and reasoning tokens on every row. Reasoning text is on 87.6% of GPT-5.4 mini's full-set rows
  (88.4% in the 450-item arm) and on 99.6% to 100% of the others'.
- Finish reasons: `stop` on every row except Gemma's: 4 rows at the 32,768 ceiling in the full set, 3 of them empty,
  and 1 empty in repeat1. The ceiling stays 32,768.
- Prompt hash check: equal to main, yes, on every row of every store.
- Gemini was served by Google AI Studio and by Google (Vertex); the main run and the 450-item arm used Google AI Studio
  alone. The served model id is the same.
- Cascade bills: paraphrase pairs $3.741 (825 rows); repeats $3.762 (1,800 rows).
- **Rows without a finish reason are asked again: decided in the stream, confirmed by the owner on 2026-10-07.**
  - What happened: Io Net ended 26 of Gemma's completions (about 1%) with no finish reason. The text was cut off
    mid-sentence at 3,111 to 12,006 completion tokens and 40 to 332 s, far below the ceiling. DekaLLM did so on none of
    its rows, and the main run had one such row in 26,970.
  - Why it matters: scored as they stood, these rows lowered Gemma's paraphrase score from 0.933 to 0.920.
  - The rule: the harness already treats `finish_reason=error` inside a 200 as a provider fault and asks again (D-148).
    `call()` now does the same for a missing finish reason, in the matched-configuration variants only; the other
    variants keep the main run's reading.
  - Handling: the rows written before the rule were moved, not deleted, to `traces/_unfinished/<variant>/`
    (`run_traces --requeue-unfinished`), and their items asked again: 15 full-set items, 5 paraphrase pairs, and 3, 1
    and 2 repeat items. All answered on the next call, finishing `stop`, for $0.048.
  - Verified 2026-10-07: each of the 26 items has a newer final row that finished `stop`, and its scored row was built
    from that row; no final row of any WS-B store is left without a finish reason.
  - Alternative: scoring those rows as they stand can be rebuilt from the archive.
- The 30 repaired items: nothing to re-run. Every WS-B run started after FREEZE DONE, on the final pool (manifest
  `f2c1dd2f…`).
- B5 files written: `DECODING_TABLE_<variant>.md`, `results/decoding_table_<variant>.json`, `TRACE_REVIEW_<variant>.md`
  and `trace_review_<variant>.json` for the five variants. The trace reviews are clean: no missing item, no row outside
  the variant's items, no second final row, no malformed line, every model "as main".

## Stores
- `scores/reasoning-medium-full/`, `scores/paraphrase-reasoning-medium/`, `scores/repeat1..3-reasoning-medium/`: every
  model scored, all rows.
  - CONFIG is current against the final manifest (`f2c1dd2f…`) and `diversity.json` (`c71061f4…`), the inputs WS-A's
    re-scored `main` store carries. It names commit `740ee25` and is DIRTY because those two inputs are uncommitted
    (WS-A's re-freeze).
  - Mean scores: full set GPT-5.4 mini 0.960, Gemma 4 0.932 (3 unusable), Gemini 0.910; paraphrase pairs 0.964, 0.933,
    0.931; repeats Gemini 0.922, 0.912, 0.927 and Gemma 0.943 (1 unusable), 0.937, 0.930.
- Stage rows, complete and current:
  - `scores/reasoning-medium-full/e5/` and `router/`, and `scores/paraphrase-reasoning-medium/e5/`. Each stage CONFIG
    lists all three models and was built from the store's current CONFIG (`store_config_sha256` equal).
  - E5, full set: 2,391 calls (GPT-5.4 mini 831, Gemini 768, Gemma 792), $7.285. Inside complete replies the judge left
    14 milestones unjudged: Gemma 13, Gemini 1.
  - E5, paraphrase pairs: 281 calls (97, 91, 93), $0.833.
  - Router, full set: 6,747 calls (2,250, 2,250, 2,247; Gemma's 3 empty rows have nothing to check), $22.764. Inside
    complete replies it left 2 steps unjudged, both GPT-5.4 mini's.
  - Served by GMICloud, Xiaomi, Novita and DigitalOcean (`judge --status`, `router --status`).
  - Dry runs, for comparison: E5 $7.79 and $0.91 at Xiaomi's endpoint; router $22.01 to $26.78.
  - The repeats are scored only, as the provider-default repeats are. The paraphrase pairs get E5 only, as the
    provider-default paraphrase arm has.
- `results/matched_config.json` (`run_traces --matched-config`; it refuses to write while a store it names is
  incomplete or the rows contradict the setting): 12 entries.
  - GPT-5.4 mini, Gemini 3.1 Flash-Lite and Gemma 4: `effort=medium`, `reasoning_store: reasoning-medium-full`, with
    their cascade stores (repeats for Gemini and Gemma only). Evidence: 0 of 2,250 main-run rows with reasoning tokens,
    2,250 of 2,250 with the parameter.
  - Qwen3-235B-2507: `none offered`, with the card's "supports only non-thinking mode" and no reasoning parameter
    listed on OpenRouter.
  - The seven others: `reasons by default`, with reasoning tokens on 2,242 to 2,250 of their 2,250 main-run rows.
  - Qwen Thinking: `sibling_of: qwen3-235b-a22b-2507`, `run: false`, `default_store: null`, with its card's "supports
    only thinking mode" and the reason it was not run.

## Facts for the setup paragraph (Phase 2 writes the sentences)
- Run with reasoning on at matched settings: GPT-5.4 mini, Gemini 3.1 Flash-Lite and Gemma 4 26B. Each used OpenRouter's
  reasoning parameter at medium effort, on all 2,250 items, with the main run's prompt, 32,768-token ceiling and
  routing. At the providers' defaults, none of the three returned reasoning tokens on any of their 2,250 main-run rows;
  with the parameter on, all three did on every row.
- Endpoint with no setting: Qwen3-235B-A22B-Instruct-2507, whose card says it "supports only non-thinking mode" and for
  which OpenRouter lists no reasoning parameter.
- Reasoning by default: the other seven, with reasoning tokens on 99.6% to 100% of their main-run rows.
- Twelfth model: none. Qwen3-235B-A22B-Thinking-2507 passes the roster rule and was calibrated on 20 items, but it was
  not run, by the owner's decision. Suggested scope sentence: Qwen3-235B-2507 is evaluated in its Instruct release; its
  Thinking release is a separate model and was not evaluated.
- For the cascade: the paraphrase pairs (275) and, for the two re-run models in the repeat set (Gemini 3.1 Flash-Lite,
  Gemma 4), the three decoding repeats (300 items each) were re-run at the same setting.
- Provider faults: Io Net cut about 1% of Gemma's reasoning-on completions short without a finish reason. Those items
  were asked again (the owner confirmed the rule).

## Open items
- **`trace_review.py`** (not WS-B's file) misreads the new variants. `variant_manifest()` takes any name starting with
  `reasoning-` as the 450-item subsample, so `reasoning-medium-full` would show 1,800 rows as "not a frozen item". The
  names `paraphrase-reasoning-medium` and `repeatN-reasoning-medium` fall through to the whole pool. Fix for the
  orchestrator, at the top of `variant_manifest`, before the existing branches:

  ```python
  from full_run_28092026 import run_traces
  if variant == run_traces.FULL_REASONING:
      return man
  if variant in run_traces.CASCADE_REPEATS:
      return {i: man[i] for i in subsamples.repeat_ids()}
  if variant == run_traces.CASCADE_PARAPHRASE:                 # the kept pairs, as the harness runs them
      keep = run_traces.kept_pairs()
      return {i: r for i, r in variant_manifest('paraphrase', man).items() if i in keep}
  ```

  Until then, `run_traces --variant <v> --trace-review` passes the right item set and runs the script unchanged. The
  `TRACE_REVIEW_*reasoning-medium*.md` files were made this way.
- **`decoding_table.py`** heads any `reasoning-` variant "(C1)" with "(D-179, D-180)". That is cosmetic for these
  variants: the request line under it is correct.
- **For WS-G**: after any re-score of these stores, run `judge --score` and `router --score` on
  `reasoning-medium-full` and `judge --score` on `paraphrase-reasoning-medium` (free). E5 and the router now hold rows
  from WS-B, so each store's cumulative `--max-usd` total grows by WS-B's spend.
- **A decision entry** (next free D-number) for the owner to write when asked: the matched configuration, the Gemma
  calibration, the no-finish-reason rule, and Qwen Thinking calibrated and not run.
