SIGNAL: C3 FILES WRITTEN 2026-10-07 15:39 UTC (quick mode: sensitivity_variants.json, flag_precision.json, depth_model.json and their sections)

# WS-C3 report

Done in quick mode, 2026-10-07. Three new scripts, each with a self-test, each run on the current stores. Nothing
committed (none asked). No paid call. No `.tex` file and no other stream's file edited; `answer.py` and `arith.py`
imported only.

## Scripts and flags

| script | what it computes | section | result file | check |
|---|---|---|---|---|
| `clause_variants.py` | FAC under four one-clause variants of the final-answer rule; verdicts moved up and down; Kendall's tau with the headline; relative error of accepted answers | `results/sections/sensitivity_variants.md` | `results/sensitivity_variants.json` | `--selftest` 32/32; the run stops unless the local copy of the rule reproduces every stored readable verdict (it did, on all 14 model-store pairs) |
| `flag_precision.py` | every digit-rule flag on a correct answer: carried precision, or other (truncated, no upstream, not reproduced); the expert's 171 confirmed slips | `results/sections/flag_precision.md` | `results/flag_precision.json` | `--selftest` 16/16; every step it reads must give the stored claim and flag counts (0 unverified of the 3,415 + 720 flags) |
| `depth_model.py` | wrong-answer rate by milestone count over readable responses; per-model logistic slope with answer kind, clustered by template; Holm | `results/sections/depth_model.md` | `results/depth_model.json` | `--selftest` 10/10 |

- All three take `--quick` (or `ENGTRACE_QUICK=1`): 1,000 resamples for the intervals, and the files carry
  `"quick": true`. Counts, slopes and p-values do not depend on resampling; only the intervals change at full
  resolution.
- Full-resolution commands for the integration pass (times on this machine in quick mode):
  `python -m full_run_28092026.clause_variants` (about 2 min),
  `python -m full_run_28092026.flag_precision` (about 15 min),
  `python -m full_run_28092026.depth_model` (under 1 min).
- Stores: `main` (the eleven roster models) and each reasoning store in `results/matched_config.json`
  (`reasoning-medium-full`: GPT-5.4 mini, Gemini 3.1 Flash-Lite, Gemma 4 26B). The `run: false` sibling is skipped,
  and the cascade stores (paraphrase, repeats) are not included. `clause_variants` and `depth_model` also report a
  `matched` configuration: each model at its reasoning store where it has one.
- Schemas: the brief's keys, plus extras. C4's loader reads them as written (C4 report, "Schema"). Extras:
  - `sensitivity_variants.json`: per model `ci`, `delta_ci`, `changed_templates`, `raised_by_symbolic`; per store
    `where_changed`; per bin set `accepted`, `numeric`, `no_numeric`, `symbolic`, `share`, `over_1pct`.
  - `flag_precision.json`: `other_truncated`, `other_no_upstream`, `other_not_reproduced`, `share_carried`,
    `share_ci`; `expert_confirmed.step_changed` and its other counters.
  - `depth_model.json`: `configurations.{main,matched}`, `pooled`, `odds_ratio`, `kinds_dropped`, `rows`,
    `wrong_rows`.

## Findings on the current store

### Clause counts (`sensitivity_variants.json`, store `main`; verdicts up / down)

| model | FAC | absolute-value clause off | last digit bounded | prescribed digits relaxed | per-part credit |
|---|---:|---:|---:|---:|---:|
| `gpt-oss-20b` | 0.814 | 0 / 50 | 0 / 33 | 38 / 0 | 8 / 3 |
| `gemma-4-26b-a4b` | 0.868 | 0 / 23 | 0 / 12 | 23 / 0 | 14 / 3 |
| `deepseek-v4.1-flash` | 0.981 | 0 / 16 | 0 / 11 | 8 / 0 | 3 / 2 |
| `qwen3-235b-a22b-2507` | 0.894 | 0 / 10 | 0 / 8 | 18 / 0 | 6 / 26 |
| `glm-5.3-flash` | 0.972 | 0 / 93 | 0 / 9 | 4 / 0 | 2 / 2 |
| `glm-5.3` | 0.948 | 0 / 80 | 0 / 5 | 0 / 0 | 5 / 2 |
| `muse-glimmer-30b` | 0.969 | 0 / 13 | 0 / 16 | 5 / 0 | 1 / 1 |
| `kimi-k3` | 0.974 | 0 / 62 | 0 / 9 | 6 / 0 | 1 / 5 |
| `gpt-5.4-mini` | 0.850 | 0 / 33 | 0 / 10 | 25 / 0 | 8 / 6 |
| `gemini-3.1-flash-lite` | 0.874 | 0 / 21 | 0 / 4 | 17 / 0 | 7 / 4 |
| `claude-sonnet-5` | 0.981 | 0 / 33 | 0 / 4 | 11 / 0 | 3 / 2 |
| all | | 0 / 434 | 0 / 121 | 155 / 0 | 58 / 56 |
| tau with the headline | | 0.927 | 0.964 | 0.964 | 1.000 |

- Matched configuration: tau 0.818 (absolute-value clause off), 0.964, 0.964 and 1.000. The reasoning store's
  three models move 66, 19 and 47 verdicts, and 18 up / 8 down under per-part credit; tau 1.000 for each variant.
- Where the changes sit (main):
  - Absolute-value clause: 434 verdicts on 22 templates. Most are on `incompressible_continuity` (48),
    `cd_dc_system_analysis` (39), `vorticity_check` (39), `heat_of_reaction_formation` (36) and `time_to_phasor` (36).
  - Last-digit bound: 121 verdicts on 38 templates. Most are on three symbolic templates: `impulse_response_from_lccde`
    25, `ber_estimation_mary` 18, `autocorrelation_rect_pulse` 17.
  - Prescribed digits: 155 verdicts on 15 of the 17 prescribing templates.
  - Per-part credit: 114 verdicts on 13 templates.
- Relative error of accepted answers (main, 22,533 accepted, 316 with word or label targets only). Of the 22,217
  with a numeric target: 92.4% within 0.2%, 6.9% from 0.2 to 1%, 0.3% (64) from 1 to 5%, 0.5% (103) above 5%.
  Of the 167 above 1%:
  - 80 stay accepted with the last-digit window capped (through the gold's own last digit or another reading);
  - 50 are on symbolic templates.
- What drives the accepted answers above 5% (read on a sample of them):
  - `answer.values` drops a unit exponent only when a letter precedes it. So the `2` of `(signal units)^2` or
    `\text{(...)}^2` is read as a number, and `match` then credits the autocorrelation target 288 as 2 x 100
    (the percent factor), with a one-unit window of ±100.
  - On responses with no final-answer heading, a small stated number on a segment that is not the final line is
    matched at a unit factor (an aliasing item: 0.5 x 10^3 ± 100 against 540).
  - The last-digit bound removes this credit. See the open items for WS-E.

### Carried precision (`flag_precision.json`)

- Main run: 11 of 3,415 flags on correct answers are carried precision (0.3%). By model: Qwen3-235B 4, Muse
  Glimmer 3, Claude Sonnet 5 2, GLM-5.3 1, GPT-5.4 mini 1, the other six 0.
- The other 3,404 split three ways:
  - truncated: 398 (11.7%), the printed result cut instead of rounded;
  - no upstream: 1,991 (58.3%), no operand has an earlier, more precise value;
  - not reproduced: 1,015 (29.7%), one has, and the recomputation still misses the printed digits.
- Reasoning store: 2 of 720 carried precision.
- The 171 slips the expert confirmed (FLAG_REVIEW_3.md), re-found in the current traces:
  - 0 are carried precision;
  - 168 are other: 17 truncated, 118 no upstream, 33 not reproduced;
  - 3 sit in steps whose text has changed (two on the virial template, one on the flame template: items re-run in
    the round-5 repair);
  - none is unflagged by the current rule, and none failed to be re-found.
- Why carried precision is rare here: the digit rule already lets each shown operand move half a unit of its last
  digit (`arith.shown_uncertainty`). So a correctly rounded intermediate that is carried forward is not flagged in
  the first place.
- How the classifier links an operand to an upstream value: the value lies within one unit of the operand's last
  digit (a rounding, a truncation or a double rounding). It also takes whole numbers of 100 or more, which
  `arith.py` treats as exact, and the values of earlier expressions.
- I read 10 of the 11 main-run cases. Each is a carried value, for example `52677 / 1758 = 29.97`, where 1758 is the
  EOQ rounded from the 1757.59 computed earlier.

### Depth model (`depth_model.json`)

Main run: 4 of 11 slopes hold after Holm.

| model | wrong / rows | slope per milestone (95% CI) | odds ratio | p Holm | holds |
|---|---:|---:|---:|---:|---|
| `gpt-oss-20b` | 318 / 2,123 | +0.227 (+0.135 to +0.319) | 1.25 | 1.4e-05 | yes |
| `gemma-4-26b-a4b` | 260 / 2,180 | +0.221 (+0.107 to +0.336) | 1.25 | 0.0011 | yes |
| `deepseek-v4.1-flash` | 21 / 2,172 | +0.277 (-0.013 to +0.568) | 1.32 | 0.41 | no |
| `qwen3-235b-a22b-2507` | 176 / 2,180 | +0.122 (-0.005 to +0.248) | 1.13 | 0.41 | no |
| `glm-5.3-flash` | 8 / 2,135 | +0.204 (-0.107 to +0.515) | 1.23 | 0.99 | no |
| `glm-5.3` | 7 / 2,081 | -0.104 (-0.359 to +0.152) | 0.90 | 0.99 | no |
| `muse-glimmer-30b` | 29 / 2,149 | +0.115 (-0.106 to +0.336) | 1.12 | 0.99 | no |
| `kimi-k3` | 30 / 2,164 | +0.125 (-0.071 to +0.321) | 1.13 | 0.99 | no |
| `gpt-5.4-mini` | 309 / 2,180 | +0.247 (+0.129 to +0.364) | 1.28 | 0.00038 | yes |
| `gemini-3.1-flash-lite` | 262 / 2,180 | +0.239 (+0.126 to +0.353) | 1.27 | 0.00038 | yes |
| `claude-sonnet-5` | 28 / 2,178 | +0.023 (-0.158 to +0.204) | 1.02 | 0.99 | no |

- The four slopes that hold belong to the four models with the most wrong answers. Qwen3-235B, fifth with 176, has
  p = 0.059 before correction. Ten of the eleven point estimates are positive.
- Pooled, with a model fixed effect: +0.204 (+0.118 to +0.290), odds ratio 1.23 per milestone, p = 3.4e-06, over
  23,722 rows on 147 templates.
- Matched configuration: 2 of 11 hold, gpt-oss-20b and Gemini 3.1 Flash-Lite at medium effort (+0.228, p Holm
  0.0056). Gemma 4 at medium effort: +0.218, p Holm 0.074. GPT-5.4 mini at medium effort: +0.115, p Holm 0.53, with
  72 wrong answers against 309. Pooled: +0.181 (+0.079 to +0.283).
- Where a model has no wrong answer on an answer kind, that kind has no finite coefficient, so it is left out of
  that model's fit. This gives the slope the full likelihood tends to. The kinds left out are listed per model;
  classification is left out for most models.

## Open items

- WS-E (owner of `answer.py`):
  - The unit-exponent reading above (`)^2` read as a 2) credits wrong or unstated values through the percent
    factor, mostly on symbolic templates. Your symbolic check raises verdicts and never lowers them, so it will not
    remove these credits.
  - Whether to change `values` is your call. Any change needs the one re-score at ANSWER FINAL.
  - The last-digit bound removes these credits, along with every other acceptance that rests on a coarse last
    digit of the response: 121 main-run verdicts down in all.
- WS-E, at ANSWER FINAL: `clause_variants.py` reads the symbolic step from `answer.py` under the names your report
  gives: `answer.symbolic_equivalence.apply(label, item, text, enabled, tol, unit, segment)` and
  `answer.SYMBOLIC_EQUIVALENCE_TEMPLATES`. It applies the step after each variant's label. If you use other names,
  its reproduction check stops on the first raised row. Tell WS-G, or me, so `symbolic_hook()` can be adapted.
- WS-D1 or WS-D2:
  - The absolute-value clause carries 434 main-run verdicts on 22 templates, and its tau is the lowest (0.927, or
    0.818 matched). The script counts these verdicts but does not judge them.
  - Some will be sign conventions (a deflection stated downward, heat released). Others may be wrong signs
    (velocity components, vorticity).
  - If the paper calls the clause a convention allowance, a reading of a sample is needed first.
- WS-C4: `flag_precision.json` is now written, so `tab:flag_precision` and its phrase can drop STAND-IN.
- Orchestrator: before my fix, a run of `flag_precision.py` printed the right-side result as an upstream value (a
  bug). I stopped it in its expert phase, and it never wrote. The file on disk carries `other_truncated` and
  `written_at_utc` 15:18:52 UTC, from the fixed run. If that stopped process were still alive and wrote late, the
  file would lack `other_truncated`; the integration pass rewrites it anyway.
