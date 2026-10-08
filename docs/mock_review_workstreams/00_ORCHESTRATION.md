# Orchestration

For the orchestrating session. Written 2026-10-06. Everything the isolated sessions need is in their briefs; this
file holds what only the orchestrator needs: the graph, the signals, the decisions, the budget, the ownership
matrix, the registry of names both phases share, the Phase 2 sequence, the master checklist and the fallbacks.

## 1. Phases and the dependency graph

**Phase 1 (parallel, 7 to 9 October): code, data, experts, figures. No tex file is edited.**

```
WS-F (extract drawing code) --FIGURES EXTRACTED--> WS-C (may now edit paper_results.py)
WS-B (harness variant + --only-items) --HARNESS READY--> WS-A (re-run of the 30 items)
WS-A (repair, gate, round 5, re-freeze) --FREEZE DONE--> WS-B (launch, or re-run 30 items later), WS-C (manifest-based tables)
WS-A (re-run + re-score) --STORES UPDATED--> WS-E (error-analysis top-up kit)
WS-B (runs scored) --RUNS SCORED + matched_config.json--> WS-C (matched family on real stores)
WS-E (E1 kappa, E2 grades, answer.py final) --KAPPA IN / GRADES IN / ANSWER FINAL--> WS-G
WS-F (PDFs approved) --CAPTIONS--> WS-D2 (figure captions are generated blocks; D2 passes them to paper_results.py's figure() calls)
```

Nothing in Phase 1 blocks on the writing. Every Phase 1 stream can start on day 1 morning.

**Phase 2 (sequential, 10 to 12 October): integration, then writing.**

```
WS-G (one re-score with provenance, stages for new rows, analyze, --write on the real tree, numbers sheet)
   -> WS-D1 (main text, back matter)  ||  WS-D2 (appendices, bibliography, listings)   [parallel, disjoint files]
   -> orchestrator: FILL markers resolved, vocabulary grep, every --check, Overleaf sync, compile,
      second mock review, revision-letter notes, Responsible NLP checklist, submit.
```

## 2. Start order and signals

| Signal | From | To | Means | Expected |
|---|---|---|---|---|
| FIGURES EXTRACTED | F | C | `paper_figures.py` exists, `paper_results.py` imports it, the current PDFs regenerate identically | day 1, first hour |
| HARNESS READY | B | A | `run_traces.py` has the full-set reasoning variant and `--only-items`; `--check` passes | day 1, first two hours |
| FREEZE DONE | A | B, C | `manifest.jsonl`, `FREEZE.json`, `pool/` carry the two templates' new instances; only 30 rows differ | day 1 afternoon |
| STORES UPDATED | A | E, C | the 30 items are re-run in every store, re-scored, judge and step-check rows present | day 2 to 3 |
| RUNS SCORED | B | C | every reasoning-on store scored with its stages; `results/matched_config.json` written | day 2 to 3 |
| KAPPA IN | E | G, D1 | `template_annotation_23092026/levels/agreement.json` written | day 3 |
| GRADES IN | E | G, D2 | symbolic grades scored; the check's validation written | day 3 |
| ANSWER FINAL | E | G | `answer.py` will not change again (its hash is in every store's provenance) | day 3 |
| CAPTIONS | F | D2 | approved PDFs in `figs/`, captions in F's report | day 3 |

A signal is a line at the top of the stream's report file: `SIGNAL: <name> <date time>`. The orchestrator polls the
reports folder; dependent sessions are told by the owner.

## 3. Decisions (record the choice in the last column)

| # | Decision | Recommendation | Chosen |
|---|---|---|---|
| D1 | Title and term | "Verifiable Process Evaluation of Engineering Reasoning from Symbolic Templates"; "process evaluation" throughout | |
| D2 | Headline configuration | Matched settings as the headline: each model at its reasoning-on setting where its endpoint offers one, the provider-default rows in the appendix with the paired change. Consequences the owner accepts: the error analysis is re-read for the re-run models that are in it (GPT-5.4 mini; Gemma 4 if it runs), the paraphrase pairs and the repeats are re-run for the re-run models (WS-B's cascade), the judge-ruling readings keep their default-run source and say so. Fallback if the re-read cannot happen: defaults stay the headline, the matched family sits beside it, the tiers are stated under both. WS-C builds both layouts behind `--headline`. | |
| D3 | Models to re-run | GPT-5.4 mini and Gemini 3.1 Flash-Lite on all 2,250; Gemma 4 26B if the calibration shows reasoning tokens; Qwen3-235B-A22B-Thinking-2507 as a twelfth row if the roster check passes | |
| D4 | The two chemical templates | Path B, decided by the owner on 2026-10-06: repair, re-certify, re-run | Path B |
| D5 | Symbolic answers | Expert grading now; a validated equivalence check for `ber_estimation_mary` and any template the grades validate; the by-numbers rule stated for the rest | |
| D6 | Prescribed digits | Keep the exact-digit requirement as the construct; report the relaxed score as a sensitivity | |
| D7 | Difficulty labels after E1 | Report agreement and keep the labels if the majority differs on a handful; adopt the majority and regenerate if it differs widely (WS-E's `apply_labels.py`, then WS-G) | The authors' labels stay (58 / 58 / 34): E1 set aside, no agreement figure in the paper, nothing relabelled or re-frozen (owner, 2026-10-08) |
| D8 | Release and contamination policy | Release templates, generator, evaluator code, all responses and scores, judge rulings and hashes under MIT; the seed at publication; a held-out seed kept; names and per-rater files withheld | |
| D9 | Figure model labels | Keep "GPT OSS 20B" in figures with one note per caption, or switch to `gpt-oss-20b`; consistent either way | |
| D11 | **Evaluator fixes before ANSWER FINAL** (WS-E amendments of 2026-10-07 evening): (a) fold unicode superscript exponents in the milestone reader; (b) stop reading a unit exponent as a number in `answer.values`; (c) the non-final-segment fallback window; (d) bound the last-digit term by min(u(y_hat), 0.01|y|) as the rule | (a) fix; (b) fix; (c) fix or document after reading C3's sample; (d) adopt, state in the scoring appendix, re-validate on the expert study | |
| D12 | **The absolute-value clause** (434 main-run verdicts on 22 templates depend on it) | Keep the clause and report the count as a sensitivity; have a 40-item sample read (E6) only if the text is to call it a convention allowance | |
| D10 | Taxonomy facts | The rule that cut the prompted 5 to 7 domains to 3 per branch; what the LLM panel did (cross-check or vote); stated in the same words in §3.1 and the appendix | |

The budget was approved on 2026-10-06 at $75 to $110 for all paid work. Each launch still prints its dry-run
estimate and stays within its stream's cap (section 4).

## 4. Budget ledger

| Stream | Item | Estimate | Cap | Actual |
|---|---|---|---|---|
| A | 30 items × 11 models, main run | $15 | | $21.24 |
| A | arms on the 6 subset items (effort, anchors, open book, tool, repeats) | $3 | | $2.96 |
| A | paraphrase pairs re-run (5 × 11), if the expert keeps them | $3 | | $2.25 (3 virial pairs × 11) |
| A | judge and step-check calls on the new rows | $5 | | $4.36 |
| A | LLM screen of the two templates | $0.10 | **$35** (raised 2026-10-06) | skipped; **A total $30.82** |
| B | GPT-5.4 mini, inference + stages | $33 | | inference $22.65; stages running |
| B | Gemini 3.1 Flash-Lite, inference + stages | $16 | | inference $5.53; stages running |
| B | Gemma 4 26B with reasoning, if supported | $13 | | supported (20 of 20 reasoned); inference $3.80; stages running |
| B | Qwen3-235B Thinking, if eligible | $20 | | eligible, priced $35.50 + $13 stages; **not run, cap binds; owner's call** |
| B | cascade: paraphrase pairs and repeats for the re-run models | $6 | **$90** | $7.50; **B inference total $39.48**, stages projected about $40 |
| E | symbolic check: none (offline); error-analysis top-up: none (expert time) | $0 | $0 | |
| C | the step-check stability sample (optional) | $2 | $3 | |
| | **Total** | **$116 if everything runs** | **$123** | |

Priority if the cap binds: A first, then B's two closed models, Gemma, the cascade, Qwen Thinking last.

## 5. Expert dispatch log

| Ask | Kit from | Sent | Returned | Scored by |
|---|---|---|---|---|
| E1 difficulty ratings, all 15 experts | WS-E | | | WS-E |
| E2 symbolic grading, electrical (and any branch with symbolic templates) | WS-E | | | WS-E |
| E3a round-5 re-certification, three chemical experts | WS-A | 2026-10-06 | 2026-10-06 (round 5: flame A A A, virial A R A; round 6: virial A A A) | WS-A, done |
| Paraphrase re-check of up to 5 pairs, one chemical expert | WS-A | 2026-10-06 | 2026-10-06 (3 virial pairs kept 3 of 3; flame pairs could not be written) | WS-A, done; kept pairs now 275 |
| E3b error-analysis top-up, three readers per affected branch | WS-E (after STORES UPDATED) | | | WS-E |
| E3c error-analysis re-read for the re-run models (only under D2 = matched headline) | WS-E (after RUNS SCORED) | | | WS-E |
| E4 facts (experts, recruitment, pay, consent, independence, study-template selection, taxonomy rule) | the authors | | | WS-D1, WS-D2 |

## 6. File-ownership matrix (Phase 1)

| Path | Owner | Others |
|---|---|---|
| `data/templates/branches/chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.py` (the virial function), `heat_effects.py` (the flame function) | A | read only |
| `template_annotation_23092026/layer0/*` outputs, `layer2/*round5*`, `layer2/CERTIFICATION.md`, `layer2/fixes_round5.md` | A | read only |
| `full_run_28092026/freeze.py`, `FREEZE.json`, `manifest.jsonl`, `pool/`, `diversity.json`, `DIVERSITY.md` | A | read only (C reads the manifest) |
| `full_run_28092026/traces/*` and `scores/*` rows of the 30 repaired items, in every store; `scores/_replaced/round5_*` | A | B writes its own variant stores; E does not write stores |
| `full_run_28092026/paraphrase/*` (the 5 pairs), `openbook/` items for the two templates | A | |
| `full_run_28092026/repair_round5.py` (new) | A | G reuses it |
| `full_run_28092026/run_traces.py`, `models.json`, `subsamples.py` | B | A uses `--only-items` after HARNESS READY |
| `full_run_28092026/traces/<reasoning variants>/`, `scores/<reasoning variants>/`, `DECODING_TABLE_<variant>.md`, `TRACE_REVIEW_<variant>.md`, `results/matched_config.json` | B | |
| `full_run_28092026/analyze.py`, `results/results.json`, `results/RESULTS.md`, `results/matched.json`, `results/single_path.json`, `results/providers.json`, `results/judge_swap_main.json`, `results/sections/{matched,single_path,providers}.md` | C1 | C2, C3 import analyze.py's pairwise functions |
| `full_run_28092026/coverage_variants.py` (new), `results/coverage_variants.json`, `results/sections/coverage_variants.md` | C2 | |
| `full_run_28092026/clause_variants.py`, `flag_precision.py`, `depth_model.py` (new), `results/sensitivity_variants.json`, `flag_precision.json`, `depth_model.json` and their sections | C3 | |
| `full_run_28092026/paper_results.py` (after FIGURES EXTRACTED; or C4 extracts first) | C4 | F edits it once, before the signal; D1 in Phase 2 |
| `docs/appendix_statistics.py`, `full_run_28092026/worked_example.py` (new), `appendices/worked_example.tex` (under `--out` only) | C5 | |
| `full_run_28092026/paper_figures.py` (new), `overleaf_source_04102026/figs/*.pdf` | F | C's `--write` calls F's module |
| `full_run_28092026/answer.py`, `validate_scorer.py`, `symbolic/` (new), `expert_kits.py`, `EXPERT_REQUEST.md`, `template_annotation_23092026/levels/` (new) | E | A and B run `score.py`, which imports `answer.py`: E announces ANSWER FINAL before G's re-score |
| `full_run_28092026/score.py`, `judge.py`, `router.py`, `decoding_table.py`, `trace_review.py` | nobody edits in Phase 1 (run only); fixes needed go to the orchestrator | |
| every `.tex`, `custom.bib`, `WRITING_RULES.md`, `docs/appendix_evaluation.py`, `docs/appendix_certification.py`, `docs/appendix_listings.py` | nobody in Phase 1 (scripts may write to `--out` scratch copies); D1, D2, G in Phase 2 | |

Phase 2: G owns every generated block and runs the scripts; D1 owns `main.tex`, `0_abstract.tex`, `1_intro.tex`,
`2_relatedwork.tex`, `3_4_design.tex`, `5_evaluation.tex`, `6_results.tex` (prose) and the phrases section of
`paper_results.py`; D2 owns `appendices/*.tex`, `7_appendix.tex`, `custom.bib`, `docs/appendix_listings.py`,
`docs/appendix_evaluation.py`, `docs/appendix_certification.py`, `WRITING_RULES.md`. Both append to `custom.bib`
under their own marker and never edit the other's block.

## 7. Registry of names shared across streams

Variants (WS-B): `reasoning-medium-full` (all 2,250 items; models by `--model`), `paraphrase-reasoning-medium`,
`repeat1-reasoning-medium` to `repeat3-reasoning-medium`; the Qwen Thinking model as a `main` roster row with key
`qwen3-235b-a22b-thinking-2507`. The 450-item arm `reasoning-medium` stays as it is.

Script flags (WS-C): `paper_results.py --out DIR` (write blocks and figures into a copy of the tex tree),
`--headline default|matched`, `--repaired` (the two templates are repaired: no "with and without" clauses),
`--text-only`; `docs/appendix_statistics.py --write` (fills marked blocks) and `--check`.

Generated block labels (WS-C writes, D1/D2 place the markers): `tab:main_results` (redesigned), `tab:single_path`,
`tab:matched`, `tab:coverage_variants`, `tab:scoring_variants` (the clause and prescribed-digit sensitivities),
`tab:providers`, `tab:depth_model`, `tab:flag_precision`, `tab:judged_steps` (moved out of Table 1), `tab:area`,
`tab:levels_agreement`, `tab:errors` (gains the no-error-removed shares), `box:worked_example`
(`appendices/worked_example.tex`). Figures keep their labels: `fig:level_bars`, `fig:error_categories`,
`fig:branch_bars`, `fig:domain_heatmap` (replaces `fig:domain_radar`).

Result files of the WS-C split (each script writes only its own; schemas in the C1 to C3 briefs; every file
carries `"quick": true` when written in quick mode and C4 prints "STAND-IN" in a caption that reads a stand-in):
`results/matched.json`, `single_path.json`, `providers.json` (C1); `coverage_variants.json` (C2);
`sensitivity_variants.json`, `flag_precision.json`, `depth_model.json` (C3); Markdown sections under
`results/sections/<name>.md`. Quick mode: `--quick` or `ENGTRACE_QUICK=1` gives 1,000 resamples and 10,000
permutations; the orchestrator's integration pass runs without it.

FILL markers (Phase 2): a sentence that waits on a number is written with `%% FILL: <what>` on the line above;
the final grep `grep -rn "FILL:" overleaf_source_04102026 current_overleaf_project` must return nothing.

## 8. Phase 2 sequence

1. All five Phase 1 reports in; decisions D2, D5, D7 recorded.
2. **WS-G** (`WS-G_integration.md`): the one re-score of every store (manifest and `answer.py` final), stages for
   the new rows only, `analyze.py`, `paper_results.py --write` on the real tree with the chosen flags, the appendix
   scripts' `--write` and `--check`, `paper_setup.py --check`, the numbers sheet `reports/NUMBERS_SHEET.md`.
3. **WS-D1** and **WS-D2** in parallel, from the numbers sheet.
4. Orchestrator: resolve FILL markers; vocabulary grep (PAPER_PLAN section 8); every `--check`; sync to Overleaf
   per `docs/OVERLEAF_WORKFLOW.md`; compile; eight pages; every `\autoref` resolves; Limitations and Ethics present
   and unnumbered; no author information.
5. Second mock review: the same prompt as the first cycle to the same two models, with `main.tex`, the compiled
   PDF and every figure in the bundle; save under `ACL Mock Reviews/second cycle/`; read for regressions only.
6. `docs/REVISION_LETTER_NOTES.md`: round 4, the B4 reading, the repair and re-run, the reasoning-on run, which
   real-reviewer points each answers. `docs/re-implementation-sep/DECISIONS.md`: D-197 onward, when the owner asks.
7. Responsible NLP checklist; commit and push when the owner says; submit.

## 9. Master checklist (review item → stream; tick when the report confirms it)

| Review item | Stream(s) | Done |
|---|---|---|
| R1 F1 / R2 W1: saturation, residual artifacts, ceiling stated | A, E (symbolic), C (ceiling numbers), D1 | |
| R1 F2 / R2 W4, W5: process-supervision framing, premise | D1 | |
| R1 C1 / R2 W11: Limitations, Ethics, Conclusion, experts, release, dates | D1, D2, E4 | |
| R1 M1 / R2 W3, Q3, Q4: MC variants, verbosity, the separation rule, tau wording | C, D1 | |
| R1 M2 / R2 W2, Q2: defective templates, symbolic scoring, prescribed digits, read-off templates, certification's blind spot | A, E, C, D2 | |
| R1 M3 / R2 Q13: sign and precision clauses, partial credit per part | C, D2 | |
| R1 M4 / R2 W6, W12, Q7, Q9: matched settings, providers, ceiling effects, open-vs-closed dropped | B, C, D1, D2 | |
| R1 M5 / R2 W14, Q5: validation wording, per-category readings, small samples | D2, C (no-error-removed shares) | |
| R1 M6 / R2 Q6: difficulty protocol and agreement, depth, listings, method hints | E, C, D1, D2 | |
| R1 M7 / R2 W9, Q12: taxonomy rule, curricular citations, per-area table, coverage wording | C (table), D1, D2 | |
| R1 M8 / R2 Q10: contamination after release | D1 (Limitations), D2 (release appendix) | |
| R1 M9a to g / R2 W13, Q8: missing or inconsistent statements | C (single-path, names, counts), D1 | |
| R1 M10 / R2 W15: related work, ChainEval baseline named, deltas | D1, D2 | |
| R1 M11 / R2 W16: worked example, dense sentences, main-text mentions, conclusion, naming | C (box), D1, D2, F | |
| R1 m1 to m13 | m1, m2, m3, m8, m12: D2; m4, m5, m9: D2; m6: C and D1; m7, m13: C and F; m10: E and D1; m11: D1 | |
| R2 W7, W8: judged-step column, carried precision | C, D1 | |
| R2 W10: scoring policy sensitivity | C, D1 | |
| R1 §5 and R2 §5 comments | F (figures), D1 and D2 (wording, citations, typography) | |

## 10. Risks and fallbacks

| Risk | Fallback |
|---|---|
| The chemical experts reject the repaired wording | A revises and re-runs the 30 items again (about $25); the paper waits for round 5 before the freeze is final |
| The experts' returns (E1, E2, E3b) are late | The paper ships without the agreement sentence and the symbolic sensitivity, and says so in Limitations; the top-up reports the reduced error-analysis sample |
| Gemma 4's endpoint ignores the reasoning parameter | Gemma stays at its default; the setup paragraph says the endpoint offers no reasoning setting |
| The Qwen Thinking model fails the roster check or is unavailable | Dropped; the setup paragraph says Qwen3-235B-2507 Instruct has no thinking mode |
| The symbolic parser cannot be validated in time | The expert-adjudicated score is a sensitivity row; the by-numbers rule is stated in the abstract |
| The tool arm cannot be re-run on the 6 repaired items | The tool experiment's instance count drops by six; the table says so |
| A reasoning-on run hits the output ceiling often | The ceiling stays 32,768 for comparability; the empty count is reported per configuration |
| Two sessions edit one file | The matrix in section 6 forbids it; a session that needs another stream's file writes the need in its report and stops that step |

## 11. Status log

**2026-10-07, 04:40 UTC.** Reports in: WS-A (final), WS-B (in progress).

- Signals received: FREEZE DONE (A, 2026-10-06 19:05 UTC, final after round 6 at 19:44), STORES UPDATED (A,
  2026-10-07 01:35 UTC), WS-A STAGES DONE (02:03 UTC), HARNESS READY (B, 2026-10-06 21:40 UTC). RUNS SCORED is
  pending: WS-B's inference is complete, validated and scored; its judge and step-check stages run since 02:07 UTC,
  restarted at 16 workers at 04:28 UTC (100 of 1,198 closed-model judge calls returned in the first six minutes; at
  that pace every stage ends on 7 October, so the parallel router launch WS-B asked about is not needed; re-check at
  noon UTC).
- Facts that change other briefs (registry additions):
  - Certification has six rounds: round 5 (flame A A A; virial A R A on physical plausibility), a 700 K cap on
    organics, round 6 (virial A A A). 150 of 150 certified; kappa 1.000 printed where every verdict approves is
    undefined (say so, as for civil and industrial).
  - The virial template's 30 items: 27 keep their ids; `work_isothermal_virial#12`, `#19`, `#24` left and `#14`,
    `#15`, `#16` came in. Subset positions keep their ids.
  - Paraphrase kept pairs: 275 (the 5 old pairs of the two templates dropped; 3 new virial pairs kept; no flame
    pair could be written, since its data block must be copied verbatim).
  - Judge-swap sample: 9 of the 220 responses sat on the two templates and are archived; 211 remain (WS-C drops them).
  - The set-aside `qwen3.8-27b` store covers 2,220 items (not re-run).
  - Matched models: GPT-5.4 mini, Gemini 3.1 Flash-Lite, Gemma 4 26B at effort medium, every row with reasoning
    tokens; Qwen3-235B-2507 offers no thinking mode; the Qwen Thinking sibling is eligible and priced but not run.
  - WS-B's full-set means (scored, stages pending): GPT-5.4 mini 0.960, Gemma 4 0.932, Gemini 3.1 Flash-Lite 0.910.
- Orchestrator actions taken: `full_run_28092026/trace_review.py` patched as WS-B specified (the three new
  variant families) and re-run on them.
- Decisions now pending for the owner (add to section 3): (a) confirm WS-B's rule that a reasoning-on row with no
  finish reason is asked again (mirrors D-148; 26 Gemma rows, archived, $0.05); (b) Qwen Thinking: run (about $50
  beyond the cap) or state the sibling in the setup paragraph; (c) D2, D5, D7 still to record.
- Unblocked: WS-E (E3b top-up kit from the re-scored main store), WS-C (manifest final; paraphrase count and the
  swap drop known), WS-F (no dependency). Not yet started per the reports folder: C, E, F.
- Open for WS-G: `judge --score` and `router --score` on `reasoning-medium-full`, `judge --score` on
  `paraphrase-reasoning-medium` after any re-score; `decoding_table.py`'s cosmetic heading for `reasoning-` variants;
  `template_annotation_23092026/layer2/markdown_scan.md` stale for the two templates (regenerate).
- Phase 2 facts to carry: the certification appendix (six rounds), the revision letter (B4, the repairs, rounds 5 and
  6, the caps, the matched configuration, Gemma's calibration, the no-finish-reason rule, Qwen Thinking priced and
  not run), DECISIONS.md entries D-197 onward when the owner asks.

**2026-10-07, evening.** WS-C split into C1 to C5 (owner's request: a few hours, not two days). Five sessions in
parallel, three to four hours each; C4 reads the others' result files by schema and uses marked stand-ins until
they land; after the five reports, the orchestrator runs every script without `--quick`, then
`paper_results.py --out` once more, and checks that no caption still says STAND-IN. WS-C's amendments of
2026-10-07 apply to every split session.

**2026-10-07, 21:30 local.** WS-C1 to C5 reports in and verified: the seven result files exist (coverage variants
at full resolution, the rest quick), every self-test passes (analyze, coverage_variants, clause_variants,
flag_precision, depth_model, worked_example, appendix_statistics), no `.tex` in the real tree was written by any
stream (the modified list equals the session-start list; `6_results.tex`'s timestamp is 6 October), FIGURES
EXTRACTED was done by C4 on WS-F's behalf.

Findings that change the paper (sources: the C1 to C3 reports and result files):
- Matched settings: the top five are the same five under both configurations; GPT-5.4 mini with reasoning on is
  sixth (0.960) and not separated from three of the top five; Gemma 0.932; Gemini 0.910; tau 0.745 between the
  orderings; paired changes +0.111, +0.064, +0.036, none within the 0.05 margin.
- Claude Sonnet 5 against DeepSeek V4.1 Flash on MC, after the re-run: as scored Holm p 0.046 (marginal), by
  matching alone separated, route-adjusted not separated (p 1.0), intermediate-only separated. The rule in the C2
  brief fails: the separation is route conformity, not progress. WS-D1 writes it that way.
- Verbosity: within-template slope +0.0021 MC per 10 numeric values shown (about +0.008 over the interquartile
  range): MC does not reward display. Reasoning tokens: MC lower in the highest quartile (harder items).
- Carried precision: 11 of 3,415 arithmetic flags (0.3%); 0 of the 171 expert-confirmed slips. Review 2's W8 is
  refuted by the data.
- Depth: the per-milestone slope holds after Holm for the four weakest models only; pooled over models the odds of
  a wrong answer rise by 1.23 per milestone. "Depth lowers every model" is not supported per model.
- Clause counts: absolute-value clause 434 verdicts, last-digit bound 121, prescribed digits 155, per-part credit
  58 up and 56 down; 92.4% of accepted numeric answers within 0.2%, 0.5% above 5%, the latter driven by two reader
  defects (D11b, D11c).
- Providers: matched differences within 0.03 over endpoints matched on at least 20 templates (the WP4.5 sentence
  says 0.03, or names GLM-5.3-Flash on StreamLake at -0.028).
- Judge swap: 211 responses, pooled agreement 0.844, kappa 0.778 on reached against not.
- Single path: the "none" column ranges 0.000 to 0.069 (gpt-oss-20b, 4 of 58).
- C4's checks: two source drifts (PARAPHRASE_REVIEW.md still says 316 returned and 277 kept; results.json is the
  4 October pass) and two failing claims (Claude's largest majority label is now calculation 10 against no error
  9; the Advanced no-error share is no longer "largely" Claude's).

For WS-G, added to its list: regenerate `PARAPHRASE_REVIEW.md` (275 kept; `repair_round5.py --paraphrase-score`
or `paraphrase_kit.py --score`); `JUDGE_SWAP.md`'s two stale phrases (220, "20 traces per model") become
generated, and `docs/appendix_evaluation.py` line 164 parses them; `docs/appendix_evaluation.py` still
hard-codes 150 and 100; the full-resolution passes after ANSWER FINAL in this order: analyze (the main pass plus
the matched family), coverage_variants (7 minutes), clause_variants (2 minutes), flag_precision (15 minutes),
depth_model, then `paper_results.py --out` once more and the STAND-IN grep; the judge jobs D11a creates (about
183 calls, about $1); optionally C1 adds the router fields to `matched.json` so `tab:judged_steps` prints full
rows for the re-run models.

**2026-10-08, morning. Where things stand and what comes next (the orchestration continues in a new session).**

Streams: A done and verified; B done and verified (RUNS SCORED); C1 to C5 done and verified (quick-mode result
files, every self-test passing, FIGURES EXTRACTED done by C4); E in progress (three kits built, not yet sent; the
symbolic check built and pre-validated; the evening amendments of 2026-10-07 list four evaluator fixes D11a to d to
do before ANSWER FINAL); F not started (only the redraws remain, with the owner's approval gate); G, D1, D2 not
started.

Order of work from here:
1. WS-E: the D11 fixes, re-validate on the expert study, dispatch E2 then E1 then E3b (the owner sends), build E3c
   if D2 is the matched headline, then ANSWER FINAL (when the E2 grades are scored, or by the evening of 9 October
   with the check left disabled).
2. WS-F beside it (no dependency).
3. WS-G after ANSWER FINAL: the one re-score of every store (`repair_round5.py --rescore --variant V` per store),
   the new judge jobs D11a creates (about 183 calls, about $1, approval first), the full-resolution passes in the
   order listed in the 2026-10-07 21:30 entry, `paper_results.py --write` into the real tree with the chosen
   flags, the numbers sheet.
4. WS-D1 and WS-D2 may start their number-independent parts before WS-G (title, abstract framing, introduction,
   related work, Limitations, Ethics, the release appendix, the certification account, the bibliography, the
   listing swap), with `%% FILL` markers for every number; they must pause while the orchestrator runs any
   `--write` into the tex tree, and resume from the numbers sheet afterwards.
5. The orchestrator's final checks, the second mock review, the letter notes, submission.

Decisions still open for the owner: D2 (headline configuration; the recommendation is matched), D5 and D7 (after
the E2 and E1 returns), D11a to d (the evaluator fixes; recommended: fix, fix, fix or document, adopt the bound),
D12 (the absolute-value clause; recommended: keep and report). The expert kits E1, E2, E3b are built and waiting
to be sent.

**2026-10-08, afternoon. WS-E done.** ANSWER FINAL posted at 12:34 (+05:00): `evaluator_pilot_17092026/evaluators/
answer.py`, `milestones.py` and the new `symbolic_equivalence.py` are final and committed (the owner asked for the
commit). GRADES IN: 362 of 362 readings from six experts, reader agreement 0.950 (kappa 0.882); 200 of 264 agreed
incorrect-or-partial symbolic verdicts are equivalent to the reference. The symbolic check is enabled on five
templates (`autocorrelation_rect_pulse`, `ber_estimation_mary`, `cd_dc_system_analysis`,
`impulse_response_from_lccde`, `incompressible_continuity`) with precision 1.000 and recall 0.990 on the graded
verdicts; the rule was settled on these grades, so those figures are not an estimate on unseen answers (the
earlier B1 and B2 readings give precision 1.000, recall 0.840). D11a, D11b and D11d are fixed, D11c documented
(`OWN_DIGIT_CAP = 0.01`). Re-validation on the expert study: answer check 0.986 non-partial (was 0.982), 0.930
three-way (0.927); milestone F1 0.926 (0.923). Preview of the re-score on the main store: 275 verdicts move (198
raised by the symbolic step, 77 lowered by the number rule), mean score per model moves by at most 0.0073; D11a
changes 123 judge prompts (paid re-calls at the re-score). E3b top-up returned (24 of 24): B2 is 148 items,
kappa 0.925; after the re-score 11 sampled answers are no longer incorrect and the sample becomes 137 (Claude 19).
E1 ratings not returned; KAPPA IN pending. E3c not built: the owner declined further expert reading, so the
error analysis stays a reading of the default-setting responses and the paper says so. E6 not built.

Consequences:
- WS-G is unblocked and should run now in the new session: every store's provenance names the old `answer.py`,
  so `score.py` refuses until the one re-score; then the judge's new calls (123 prompts in main plus the other
  stores' dry-run count; approval first), `expert_kits.py --score` (B2 to 137 items, B1 against the new verdicts),
  the full-resolution passes, the blocks.
- Before the full pass, two small fixes: `clause_variants.py`'s headline copy of the match rule must carry
  `answer.OWN_DIGIT_CAP` (its reproduction check otherwise stops after the re-score; the old unbounded window
  becomes the variant); `worked_example.py`'s default item gains two matched milestones under D11a, so its account
  must be re-checked after the re-score.
- `docs/check_plan_claims.py` lines 154 to 156 and the prose in `5_evaluation.tex` (line 9) and
  `appendices/validation.tex` (lines 74 to 76) carry the old validation figures: 0.930, 0.986 and
  0.927 / 0.924 / 0.926 now (WS-D1, WS-D2).
- The scoring appendix must state the own-digit cap and the five enabled templates (WS-D2).
- `tab:levels_agreement` stays a stand-in until E1 returns; WS-G re-runs `docs/appendix_statistics.py --write`
  then, and D7 is decided then.

**Correction, 2026-10-08 afternoon.** WS-E's report says its files were committed and pushed on the 8th at the
owner's request. They were not: no commit exists after 2026-10-07, master equals origin/master, and the WS-E files
(answer.py, milestones.py, symbolic_equivalence.py, symbolic/, levels/, expert_kits.py, EXPERT_REQUEST.md,
validate_scorer.py, SCORER_VALIDATION.md, the report) are modified or untracked in the working tree, as is all
other Phase 1 work (74 entries). The owner decides whether to commit the Phase 1 baseline before WS-G, so that
every store's provenance after the re-score names a commit holding the final evaluator rather than a dirty tree.
