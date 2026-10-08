# WS-G: Integration (Phase 2, first)

## Mission

Bring every Phase 1 output together once: the final re-score of every store under the final manifest, labels
and `answer.py`; the judge and step-check stages for any row still missing; every analysis script; the
generated blocks and figures written into the real tex tree with the chosen flags; every `--check`; and a
numbers sheet from which the two writing sessions work. This is the one session that writes generated blocks
into the real tree before the writers start.

## Preconditions (all five Phase 1 reports in `reports/`)

- WS-A: FREEZE DONE and STORES UPDATED; `repair_round5.py --rescore` exists.
- WS-B: RUNS SCORED; `results/matched_config.json` written.
- WS-C: scripts and flags in place; `paper_results.py --out` tested.
- WS-E: ANSWER FINAL (2026-10-08 12:34: `answer.py`, `milestones.py` and `symbolic_equivalence.py` as committed in
  `c7894ed`); GRADES IN; the top-up scored. KAPPA IN is not needed: E1 is set aside and D7 keeps the authors' labels.
- WS-F: CAPTIONS; approved PDFs in `figs/`.
- Decisions in `00_ORCHESTRATION.md` section 3: D7 recorded (the authors' labels stay); D5 as WS-E implemented it
  (the check enabled on five templates); D2 (`--headline`) recorded before G4.

## The final evaluator (reviewed 2026-10-08)

Two independent reviews of WS-E's changes to `answer.py` and the orchestrator's check found it final: no code change
before the re-score.
- The `answer.py` the stores were scored with (`1e81131`) reproduces every stored main-run verdict. The committed
  one moves exactly the verdicts `symbolic/RESCORE_PREVIEW.md` lists for the eleven models (275: 198 raised by the
  symbolic step, 77 lowered by the 1% bound, 9 of these with D11b), each explained by one change. D11a moves no
  answer verdict.
- `qwen3.8-27b` (set aside, not in the paper) also has a file in `scores/main`, and its verdicts move too, so the
  re-score's own count for `main` is above 275. The paper's figure is the eleven models'.
- Three bracket-exponent forms that D11b does not drop (`(V)²`, `(V)^(2)`, `(x)^{2.5}`) change no verdict in any
  store, under any of the five labels `score.py` records. They are a reader limit, not a fix.
- `score.py` does not store the symbolic step's per-item detail, so a failure inside the check would be silent. None
  occurs on the five enabled templates in any store.
- Run the re-score on this machine's environment (Python 3.13.13, sympy 1.14.0): the symbolic step uses sympy, and
  the pin covers the module, not the library.
- The milestone cache is current (its sidecar names the committed `milestones.py`); do not rebuild it.

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, `master` only. No other session edits code or stores while this
  runs; the writers start after INTEGRATED.
- Paid calls only for the judge and router jobs the re-score creates (D11a changes 123 judge prompts in `main`;
  the dry run counts the other stores): dry run first, the owner's approval of its estimate, cap $5.
- No commits or pushes unless the owner asks. When asked: one short line, no body, push at once.
- Every number in the numbers sheet names the file that printed it.

## Files you own in this phase

Everything under `full_run_28092026/` (run scripts; edit only to fix an integration defect, and record the
fix), `docs/appendix_certification.py` and `docs/appendix_evaluation.py` (teach them round 5 and the new
expert counts: these are code changes, so they belong here, not with the writers), the generated blocks in every
`.tex` (written by the scripts), `reports/NUMBERS_SHEET.md`, `reports/WS-G_report.md`.

## Steps

### G1. Labels: not run

D7 (owner, 2026-10-08): the authors' labels stay (58 Easy, 58 Intermediate, 34 Advanced). Nothing is relabelled or
re-frozen: `apply_labels.py` is not run, and the manifest's `level`, `FREEZE.json` and `template_inventory.csv` stay
as they are. The paper describes the labelling procedure
(`template_annotation_23092026/levels/difficulty_labelling_protocol.md`) with no agreement figure, so
`tab:levels_agreement` is not placed. `docs/appendix_statistics.py` still lists that block and reads
`levels/agreement.json`, which does not exist: if its `--check` asks for the block, remove `tab:levels_agreement`
from the script's block list rather than place a stand-in.

### G2. The one re-score of every store

For each of the 16 live stores (`main`, `flagship`, `flagship-reasoning-medium`, `openbook`, `openbook2`, `paraphrase`,
`paraphrase-reasoning-medium`, `reasoning-medium`, `reasoning-medium-full`, `repeat1` to `repeat3`,
`repeat1-reasoning-medium` to `repeat3-reasoning-medium`, `tool`; not the archive `main_pre_d137`; the Qwen Thinking
row was not run): `python -m full_run_28092026.repair_round5 --rescore --variant <v>` (moves the stage folders aside,
`score.py --replace`, restores the stage rows of unchanged items, verifies unchanged rows identical except
provenance). No stream runs `score.py` on a store before this.

Then `judge.py --status` and `router.py --status` per store. D11a changes the matched milestones of 216 `main`
responses and so 123 judge prompts: before any call, confirm by reading `judge.py`'s job selection that these
responses are re-queued and their old rulings not restored. Dry run, the owner's approval, then the calls. Run
`judge --score` and `router --score` on `reasoning-medium-full`, and `judge --score` on
`paraphrase-reasoning-medium`. Record bills.

Confirm the changed verdicts against `symbolic/RESCORE_PREVIEW.md` (275 over the eleven models: 198 raised, 77
lowered; `qwen3.8-27b`'s on top) and record the counts.

### G2b. Analysis-script fixes the re-score requires (no evaluator change)

- `clause_variants.py`: its headline copy of the match rule takes `answer.OWN_DIGIT_CAP`; the old unbounded window
  becomes the variant. Otherwise its reproduction check stops on the first row the cap moved.
- `worked_example.py`: its item (Qwen3-235B, `gas_viscosity_kinetic_theory#0`) gains two matched milestones under
  D11a (3 of 7); re-check its account of what was matched and what was judged.
- `docs/check_plan_claims.py` lines 154 to 156: the validation figures become 0.986, 0.930 and 0.926.

### G3. Analyses

In this order, each printing into `RESULTS.md` and `results.json`, all without `--quick` (full resolution):
`analyze.py` (the default family, the matched family from `matched_config.json`, sensitivity, endpoints, Q4);
`coverage_variants.py`; `clause_variants.py`; `flag_precision.py`; `depth_model.py`; `judge_swap.py` (with WS-A's
dropped responses: 211 remain); `residual_incorrect.py`; `threshold_appendix.py` if it reads the store;
`python -m full_run_28092026.expert_kits --score full_run_28092026/expert_request/expert_requests_filled`, so that
`EXPERT_REQUEST.md` is current (B2 becomes 137 items; B1 is scored against the new verdicts); `decoding_table.py`
for every variant (with the heading fix for `reasoning-` variants); `trace_review.py` likewise. Then regenerate
`PARAPHRASE_REVIEW.md` (275 kept pairs) and `template_annotation_23092026/layer2/markdown_scan.md`, run
`paper_results.py --out` once more, and check that no caption still says STAND-IN.

### G4. Generated blocks and figures into the real tree

1. Place the block markers the registry names (`00_ORCHESTRATION.md` section 7) in the tex files where the writers
   will want them (the writers may move a block later; the markers travel with it): new tables in
   `appendices/results.tex` (`tab:single_path`, `tab:matched`, `tab:coverage_variants`, `tab:scoring_variants`,
   `tab:depth_model`, `tab:flag_precision`, `tab:judged_steps`), `appendices/models.tex` (`tab:providers`),
   `appendices/taxonomy_content.tex` (`tab:area`; not `tab:levels_agreement`, see G1), `appendices/scoring.tex`
   (`\input{sections/appendices/worked_example}` with the file written by `worked_example.py`),
   `appendices/branch_domain.tex` (`fig:domain_heatmap` replacing the radar block), `6_results.tex` (the
   redesigned `tab:main_results`, the figure blocks with WS-F's captions), `appendices/results.tex`
   (`fig:level_bars_all`).
2. `python full_run_28092026/paper_results.py --write --headline <D2> [--repaired]` on the real tree (figures
   drawn by `paper_figures.py` from the final `results.json`); paste WS-F's captions into the `figure()` calls
   first.
3. `python docs/appendix_statistics.py --write` then `--check`.
4. `python docs/appendix_certification.py`: teach it rounds 5 and 6 from WS-A's `RESULTS_round5.md` and
   `CERTIFICATION.md`; `--check` will fail until WS-D2 writes the prose; leave it failing with a note.
5. `python docs/appendix_evaluation.py`: update the expert-reading counts (B1, B3 denominators; B2 top-up) and
   the symbolic rule's numbers from WS-E's files; replace its hard-coded 150 and 100; make the two phrases of
   `JUDGE_SWAP.md` that its line 164 parses (220 responses, "20 traces per model") generated; `--check` likewise
   waits for WS-D2.
6. `python full_run_28092026/paper_setup.py --check`: the setup numbers (responses in the main run, the
   configurations, the twelfth model, the output ceiling) change; update its data and leave the prose to WS-D1,
   noting which sentences it checks.
7. `python full_run_28092026/paper_results.py --check` passes on the real tree.

### G5. The numbers sheet (`reports/NUMBERS_SHEET.md`)

Write `full_run_28092026/numbers_sheet.py` (new) that prints, from `results.json`, `matched_config.json`,
`agreement.json`, `grades.json`, `SYMBOLIC_CHECK.md`, `EXPERT_REQUEST.md`, `RESULTS_round5.md`,
`CERTIFICATION.md`, the decoding tables and the Phase 1 reports, every number the writers need, each with its
source, grouped by section of the paper:

- Setup: models and configurations (which ran with reasoning on, at what setting; which endpoints offer none;
  the twelfth model), responses in the main run and in the matched configuration, the output ceiling, empty
  responses per model and configuration, query dates.
- Results: Table 1 under the chosen headline (and the other layout); tiers and letters under both
  configurations; the paired change per re-run model on the full set; the readable-only reordering; the
  single-path table's ranges; the branch spread at three decimals; the domain lows after the repair; the level
  gaps and their Holm outcome; the depth model's slopes; MC as scored, by matching, route-adjusted,
  intermediate-only, and the pairwise outcomes (Claude vs DeepSeek under each: the rule's result); the verbosity
  slope; tau between orderings with its interval; judge-decided shares; judge-swap shifts with the dropped
  responses; arithmetic flags with carried-precision split; judged step flags and precision; clause counts;
  prescribed-digit relaxed scores; per-part credit; the four repeat models and their 300 instances; paraphrase
  kept pairs (new count) and the per-model bounds; the experiments table (with the six repaired instances in
  or out).
- Error analysis: counts by model and level after the top-up (E3c did not run: the readings are of the
  default-setting responses, and the paper says so), no-error shares, kappa; what remains of the "templates whose wording does not pin the answer"
  paragraph (the six near-miss templates only).
- Certification: six rounds, the round-5 and round-6 verdicts and hand checks, what rounds 5 and 6 re-certified
  (one sentence without
  history words), the hand-check "comparable number" explanation (504 of 510), the undefined kappa cells.
- Difficulty: the level counts per branch and overall (58 / 58 / 34; `python -m full_run_28092026.export_levels`
  prints them) and the labelling procedure (`levels/difficulty_labelling_protocol.md`); no agreement figure.
- Symbolic: the enabled templates, validation precision and recall (in-sample, with the independent B1 and B2
  readings beside them), verdicts changed over the eleven models, what stays by the numbers.
- Validation: the expert-study figures under the final evaluator (`SCORER_VALIDATION.md`, "now" column), the
  judge's unreturned calls (17 of 240), the two F1 measurements, the per-category readings with
  their new denominators.
- Providers: the matched differences and the endpoints with few templates.

### G6. Post INTEGRATED

`SIGNAL: INTEGRATED <date time>` at the top of `reports/WS-G_report.md`, with the list of `--check` outcomes
(which pass, which wait for prose) and the exact regeneration command lines for the writers to re-run after
they edit phrases.

## Report template (`reports/WS-G_report.md`)

```
SIGNAL: INTEGRATED <date time>

# WS-G report
- Stores re-scored: list, verify-rescore outcomes, stage rows run (count, bill).
- Verdicts changed by the symbolic check: count per template and model.
- Scripts run and their RESULTS.md sections.
- Blocks written (labels, files); figures drawn (files).
- --check outcomes: paper_results, appendix_statistics, paper_setup (pass); appendix_certification, appendix_evaluation (waiting for prose: which sentences).
- Regeneration commands for the writers.
- Open items.
```
