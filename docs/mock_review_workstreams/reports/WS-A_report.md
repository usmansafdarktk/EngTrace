SIGNAL: FREEZE DONE 2026-10-06 19:05 UTC (final: round 6 approved work_isothermal_virial 3 of 3 at 19:44 UTC; adiabatic_flame_temperature final since 17:53 UTC)
SIGNAL: STORES UPDATED 2026-10-07 01:35 UTC (main: re-run, re-scored with verification, judge and step-check rows for every answered row of the 30 repaired items)
SIGNAL: WS-A STAGES DONE 2026-10-07 02:03 UTC (every judge and step-check run of WS-A has ended; the shared reply stores scores/_judge/e5_replies.jsonl and router_replies.jsonl are free for WS-B)

# WS-A report

Final, 2026-10-07 02:10 UTC.

Changed from the brief, by the owner's decisions in session:
- The virial template needed a second certification round, round 6: round 5's rejection led to the
  700 K cap on organics.
- The re-run used the current harness by plain resume rather than WS-B's `--only-items`, launched
  before HARNESS READY.
- The stream's cap was raised from $30 to $35 after a cost drift (below).

## The repair

- virial (`work_isothermal_virial`):
  - The question now names one reading. It is a closed system (a piston-cylinder device),
    compressed isothermally and mechanically reversibly, and the quantity is the work done on the
    gas, W = -∫P dV from the initial to the final molar volume, positive for a compression. It
    writes out the form of the equation, the two-term virial equation in its volume form,
    Z = PV/(RT) = 1 + B/V, and names B's source, the Pitzer correlation
    B·Pc/(R·Tc) = B0 + ω·B1 with Abbott's equations for B0 and B1.
  - Step 1 of the gold names the same reading and form.
  - Two behaviours changed besides:
    - (1) A draw whose printed work lies more than 0.1% (half the answer tolerance) from the exact
      work for its stated data is redrawn. Re-solving every question from its own text had found the
      gold outside the 0.2% tolerance at 13 of 500 gate seeds, all helium; the cause is display
      rounding.
    - (2) After round 5, organic compounds are capped at 700 K. 11 of the 16 organics leave the
      template (n-butane to n-octane, benzene, toluene, p-xylene, methanol, ethanol, acetone);
      methane, ethane, ethylene, propane, propylene and the 11 inorganics remain. The substituted
      lines print `·` and `^` where they printed `*` and `**`, so the review app stops italicising.
  - Ranges, constants and arithmetic are unchanged. At 193 of 500 seeds the instance states the
    same numbers and its solution from Step 2 on is identical but for the glyphs. The other 307 are
    redraws, every one explained: 292 by the cap and 15 by (1).
- flame (`adiabatic_flame_temperature`): a data block is added to the question.
  - Heats of formation: every species in the equation, at 298.15 K, in kJ/mol, with water as vapour.
  - Heat capacities: the products' Cp/R coefficients exactly as `CP_PARAMS_COMBUSTION` holds them,
    with the form Cp/R = A + B·T + C·T² + D·T⁻², its units and R = 8.314 J/(mol·K).
  - Dry air at 3.76 mol N2 per mol O2, products without dissociation, and the answer to the nearest
    kelvin.

  Nothing else changed: the solution is identical at 500 of 500 seeds and each question is the old
  one with the block appended. Shared files (`constants.py`, `_emission.py`) are untouched.
- Evidence (`python -m template_annotation_23092026.layer2.round5_checks --seeds 500 --kit 6 --pool`,
  record `round5_checks.md`): every question, solved again from its own text by independent code,
  reaches its gold.
  - Virial: within 0.1% at 500 of 500 seeds. The other readings the experts named miss: the
    pressure-explicit form at 53 of 500, steady-flow work at 469 of 500.
  - Flame: rounds to the gold kelvin at 500 of 500.
  - The 30 frozen items and the 5 round-6 kit instances pass too.
- Gate: `python -m template_annotation_23092026.layer0.gate --seeds 500` gives 150 of 150 passing.
  T1 has no failures on the two templates (4,000 and 4,982 checks); the register (1,494 lines) and
  the advisories (T5 60, T7 78) are as before. Tie census
  (`python -m template_annotation_23092026.layer0.tie_census`): 12 templates with a tie, neither of
  these two, as before.
- LLM screen: skipped, nothing billed. Its passes are tied to the corpus they started on and refuse
  once a prompt hash moves, and a third pass is refused by design, so it cannot run on two templates
  without a code change.
- Five instances of each template were read end to end, flame including one at 0% excess air. The
  full record is `template_annotation_23092026/layer2/fixes_round5.md`.

## Round 5 and round 6

- Round 5: kits built 2026-10-06 (`layer2/dist_round5/`, seeds 2501-2505, hand check 2501),
  returned 2026-10-06 (`RESULTS_round5.md`).
  - flame: A A A.
  - virial: A R A. The rejection was on physical plausibility: organics were compressed far above
    their decomposition range, the objection of round 1 (D-094). A minor note concerned the app's
    italics.
  - Hand checks: 6 of 6 matched within 1%.
- Revision: the owner chose a 700 K cap on organics, within the expert's 700-750 K, together with
  the glyph fix.
- Round 6 (virial only): kits built 2026-10-06 (`layer2/dist_round6/`, seeds 2601-2605, hand check
  2601), returned 2026-10-06 (`RESULTS_round6.md`).
  - virial: A A A.
  - Hand checks: 3 of 3 matched.
- `CERTIFICATION.md` (`certification.py` over rounds 1-6): 150 of 150 certified. 126 were last
  reviewed in round 1, 16 in round 2, 5 in round 3, 1 in round 4, 1 in round 5 (flame) and 1 in
  round 6 (virial). Every template's current code regenerates the instances its experts saw.
- Facts for the certification appendix (Phase 2 writes the sentences; sources `CERTIFICATION.md`,
  `RESULTS_round5.md`, `RESULTS_round6.md`):
  - Rounds: 6. Round 5 reviewed 2 templates (6 verdicts) and round 6 reviewed 1 (3 verdicts).
  - What round 5 re-certified and why: the two templates whose questions now state their reading and
    their data.
  - What round 6 re-certified and why: the virial template with organic compounds held below their
    decomposition temperature.
  - Every template is certified by a unanimous latest round.
  - Caveat: `RESULTS_round6.md` prints Fleiss' kappa 1.000 where every verdict approves. Kappa is
    undefined there (0/0), the issue the review raised for civil and industrial; use percent
    agreement or AC1.

## The freeze

- Writes (record `round5` in `FREEZE.json`, with each id's sha256 before round 5 and now, the
  template files' sha256, the previous manifest and pool hashes, and the history in `refreezes`):
  - 2026-10-06 17:53 UTC:
    `python -m full_run_28092026.freeze --only-templates template_work_isothermal_virial,template_adiabatic_flame_temperature --write --record round5`.
  - 19:05 UTC, after the cap:
    `python -m full_run_28092026.freeze --only-templates template_work_isothermal_virial --write --record round5 --amend`.
- Rows (`python -m full_run_28092026.repair_round5 --verify-manifest` against the pre-round-5 copy):
  - all 2,220 rows outside the two templates are identical;
  - 27 rows changed only their sha256: `adiabatic_flame_temperature#0`-`#14` and
    `work_isothermal_virial#0`-`#11`;
  - 3 ids left the selection (`work_isothermal_virial#12`, `#19`, `#24`) and 3 came in
    (`work_isothermal_virial#14`, `#15`, `#16`), by the unchanged selection rule over the re-drawn
    walk;
  - the subsample positions (0, 5, 10; 0, 7) keep their ids.
- `python -m full_run_28092026.freeze --verify`: VERIFY OK (2,250 items regenerate byte-identically), FILES OK.
- `diversity.json`: only the virial row changed. Pool: question skeletons 12 to 9, reasoning paths 2
  (unchanged). 500-draw reach: question skeletons 27 to 16, paths (upper and lower) 3 to 2.

## Inference and stages

Archive first (`repair_round5 --archive`, 2026-10-06 19:57 UTC). It moved the 1,896 rows written
for the old questions into `traces/_replaced_round5/` and `scores/_replaced/round5_20261006T195726Z/`,
with `MOVED.jsonl`. Re-run: `repair_round5 --rerun`, detached under keep-awake from 20:04 to 21:31
UTC. It ran 13 per-model chains in parallel (main first, then the model's arms), 44 steps, all exit 0.

| store | rows re-run | bill | judge calls (bill) | step-check calls (bill) | notes |
|---|---:|---:|---:|---:|---|
| main | 330 | $21.24 | 252 ($1.72) | 271 ($2.01) | 11 roster models; 63 of 330 empty answers (GLM-5.3 26, GLM-5.3 Flash 13, Kimi K3 11, DeepSeek V4.1 Flash 5, gpt-oss-20b 4, Claude 2, Muse 2); set-aside `qwen3.8-27b` archived, not re-run |
| reasoning-medium | 12 | $0.34 | 12 ($0.08) | 12 ($0.08) | |
| flagship | 12 | $0.30 | 9 ($0.05) | - | DeepSeek V4 Pro 3 empty at the ceiling |
| flagship-reasoning-medium | 6 | $0.61 | 6 ($0.03) | - | |
| openbook2 | 18 | $0.96 | 16 ($0.09) | - | items rebuilt (`openbook.py --build --version 2`); equation blocks unchanged |
| tool | 12 | $0.70 | 13 ($0.09) | - | sandbox self-test passed; ran unattended |
| paraphrase | 33 | $2.24 | 25 ($0.16) | - | the 3 virial pairs × 11 models |
| repeat1-3 | 48 | $0.09 | - | - | never judged, as before |
| openbook (v1) | 0 | - | 0 | - | obsolete: archived only, judge rows rebuilt from stored replies |

- Bills:
  - inference $26.45 (`bill_breakdown` of the new rows' `billed_usd`);
  - judge $2.28 (reply store $49.82 to $52.10);
  - step-check $2.08 (reply store $107.26 to $109.34);
  - paraphrase writing $0.01.
  - **Total $30.82 within the $35 cap.**
- Cost drift: at 20:17 UTC the projection reached $30-32 against the original $30 cap. Claude Sonnet
  5 writes about twice as much on the repaired questions (main $5.34 against $2.85 before), and
  GLM-5.3 runs to the 32,768-token ceiling. The owner raised the cap to $35.
- Paraphrase:
  - Written for the 6 subset items: virial 3 passed, flame 0 of 3 (its data block must be copied
    verbatim, which fails the near-copy check).
  - The assigned expert kept 3 of 3; the other two chemical experts also kept 3 of 3 (agreement
    only).
  - Dropped: the 5 old pairs of the two templates.
  - New kept count: **275** (`analyze.accepted_pairs`: 275 kept, 39 rejected, 0 outstanding). The
    manifest passes 314.
- Tool arm: re-run (12 rows); the experiment stays at 450 instances for Claude Sonnet 5 and GPT-5.4
  mini.
- Judge-swap sample (`scores/main/e5_grok-4-6`): **9 of its 220 responses sit on the two templates**
  (flame 6, virial 3; DeepSeek V4.1 Flash 1, Gemma 4 2, GPT-5.4 mini 2, gpt-oss-20b 2, Kimi K3 1,
  Qwen 1). They are archived, not re-judged, and 211 remain. WS-C drops them from the swap table.
  The folder's `CONFIG.json` still names the pre-round-5 store.
- Step-check gaps: two step-check rows of other items lack a reply, as they did after the original
  run (`rotating_unbalance#8` Gemma 4, `annulus_flowrate#6` Qwen). None of the repaired items does.
  The main judge run logged one "abandoned" call, yet every sent job has a reply.

## Stores

| store | CONFIG current | verify-rescore | stage rows complete |
|---|---|---|---|
| main | yes | yes, 26,640 of 26,640 | yes: judge 267 of 267 answered repaired rows, step-check 267 of 267 |
| reasoning-medium | yes | yes, 888 of 888 | yes (judge and step-check) |
| flagship | yes | yes, 888 of 888 | yes (judge) |
| flagship-reasoning-medium | yes | yes, 444 of 444 | yes (judge) |
| openbook2 | yes | yes, 1,197 of 1,197 | yes (judge) |
| tool | yes | yes, 888 of 888 | yes (judge) |
| paraphrase | yes | yes, 3,421 of 3,421 | yes (judge) |
| repeat1, repeat2, repeat3 | yes | yes, 1,184 of 1,184 each | no stages, as before |
| openbook (v1) | yes | yes, 1,269 of 1,269 | yes (judge rows rebuilt); its 6 repaired items are not re-run |

- verify-rescore: every row of an item the round does not cover equals its archived row. The only
  difference is format: four tool-arm keys (`turns`, `tool_calls`, `tool_limit`, `tool_errors`)
  present as None in stores scored before D-184.
- `answer.py` differs from the stores' previous CONFIG (D-187) and moved no row.
- Trace reviews (`trace_review.py`, every store): nothing missing, duplicated, outside the frozen set
  or malformed, and the pilot prompt everywhere. Two gaps are intended: `qwen3.8-27b` and `openbook`
  v1 lack the repaired items. Decoding tables regenerated.
- `full_run_28092026/repair_round5.py` modes: `--list`, `--verify-manifest`,
  `--archive [--variant V] [--dry-run]`, `--archive-paraphrases`, `--paraphrase-kit`,
  `--paraphrase-score DIR`, `--rerun`, `--stages --variant V`, `--rescore --variant V`,
  `--verify-rescore --variant V`, and `--check` (a self-test on a temporary tree; passes).
- The exact command WS-G runs to repeat the re-score, for each store in the table above:
  `python -m full_run_28092026.repair_round5 --rescore --variant <V>`. It moves the stage folders
  aside, runs `score.py --replace`, restores them and verifies. After it, `--stages --variant <V>`
  rebuilds the judge and step-check rows from the stored replies, at no cost unless a row is new.
- Local backups, checked file by file, with `.sha256` beside each, in `~/EngTrace_private_backup/`:
  - `full_run_pool_round5_2026-10-07.zip` (43 files);
  - `full_run_traces_round5_2026-10-07.zip` (169 files; WS-B's stores excluded);
  - `full_run_scores_round5_2026-10-07.zip` (170 files; expert labels, `_replaced` and WS-B's stores
    excluded).

  The seed is unchanged, and its 2026-09-28 backup stands.
- Off-machine copies, made 2026-10-07 on the owner's go-ahead: new versions of the three private
  Kaggle datasets.
  - The pool copy is a 45-file pool-and-seed zip in the 2026-09-28 layout (pool, seed, manifest,
    FREEZE.json).
  - The traces and scores copies are the two zips above.
  - Each was downloaded back and compared file by file: 45 of 45, 169 of 169 and 170 of 170
    identical, with the `.sha256` files identical too.

## Open items

- Commits: none made (none asked). This stream's files are as follows; the experts' files stay
  local and gitignored.
  - Templates: `data/templates/.../heat_effects.py` and `volumetric_properties_pure_fluids.py`.
  - `template_annotation_23092026/`:
    - `layer0/`: `gate_report.*`, `tie_census.*`;
    - `layer2/`: `round5_checks.py`/`.md`/`.json`, `fixes_round5.md`, `RESULTS_round5.md`,
      `RESULTS_round6.md`, `CERTIFICATION.md`.
  - `full_run_28092026/`:
    - freeze: `freeze.py`, `FREEZE.json`, `manifest.jsonl`, `diversity.json`, `DIVERSITY.md`;
    - `repair_round5.py`;
    - arms: `paraphrase/manifest.jsonl`, `PARAPHRASE.md`, `openbook/manifest2.jsonl`;
    - tables and reviews: `DECODING_TABLE{,_flagship,_flagship-reasoning-medium,_openbook,_openbook2,_reasoning-medium,_tool}.md`
      with `results/decoding_table*.json`, and `TRACE_REVIEW*.md` with `trace_review*.json` for
      main and the ten arms above.
  - `.gitignore`: one line, so the experts' paraphrase returns are never committed.
  - This report.
- Not this stream's: `EVALUATION_GUIDE.md`, `models.json`, `run_traces.py` and the
  `*-reasoning-medium*` tables are WS-B's working changes.
- WS-C:
  - three virial item ids changed (above);
  - the paraphrase kept count is 275;
  - drop the 9 judge-swap responses;
  - the set-aside `qwen3.8-27b` store covers 2,220 items;
  - `template_annotation_23092026/layer2/markdown_scan.md` (not this stream's) is stale for the two
    templates. Its owner regenerates it with `python -m template_annotation_23092026.layer2.markdown_scan`.
- WS-E: STORES UPDATED is posted; the error-analysis top-up kit can be built from the re-scored main
  store.
- WS-B: free to use the judge and step-check reply stores (STAGES DONE above). The re-drawn items
  are in the frozen pool, so its stores need no B9 archive for runs launched after FREEZE DONE.
- Phase 2:
  - the certification appendix facts and the kappa caveat (above);
  - the revision letter: the B4 reading, the two repairs, rounds 5 and 6, the cap;
  - DECISIONS.md entries (D-197 onward) when the owner asks.
- `PARAPHRASE.md`'s "billed" line now excludes the archived attempts ($0.0031).
