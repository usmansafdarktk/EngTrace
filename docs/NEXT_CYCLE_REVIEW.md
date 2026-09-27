# EngTrace, third submission: where it stands and what is left

Written 2026-09-27 against the May 2026 submission (`docs/_ARR_May__EngTrace.pdf`), both
rebuttals (`docs/EngTrace Rebuttals Jan 2026.docx`, `docs/EngTrace_Rebuttal_Jul2026 (1).docx`),
the revision letter, `docs/re-implementation-sep/EngTrace_Suggested_Actions.pdf`, and the
September record: `docs/re-implementation-sep/` (template rework, D-001 to D-109),
`template_annotation_23092026/` (certification), `evaluator_pilot_17092026/` (the evaluator
pilot and X1). Every number below is quoted from one of those files.

The short version: the September work answers the two objections that sank both cycles,
the evaluator's validity and the certification's credibility, and it does so with evidence
rather than argument. What remains is (1) a dozen concrete items before the pool is frozen
and inference starts, most of them free; (2) a results section that reports the statistics
the July rebuttal promised, on the new run; (3) two experiments that answer the objections
every reviewer repeated and that the current plan does not fund: a paraphrase robustness
test and a tool-use condition; and (4) an honest account in the paper of what changed and
why, including that the new branches were AI-drafted and expert-certified.

---

## 1. Every objection from both cycles, and what now answers it

| Objection (reviewer) | Answered by | Still missing |
|---|---|---|
| Judges share families with evaluated models; framework validated by the kind of system it uses (yAYU 1, 9W1B 2, gFWV 4) | E1: non-suite panel moves per-model scores by at most 0.004 (RESULTS_E1). New stack's only judge is MiMo-V2.5-Pro, outside every roster family (D-088, D-105). Template screen uses non-suite judges (D-093). Exposure table in JUDGE_SELECTION. | The literal promise: leave-one-judge-out on the Tribunal with a placebo and untouched-model controls. Computable offline on the 300 pilot traces from the stored per-judge votes, $0 (section 3.4). |
| rho = 0.632 on n = 100, no CI, "only moderate" (yAYU 3, 9W1B 2, gFWV 2) | X1: 300 traces, 15 experts, 3 per trace, verification and adjudication rounds. Human ceiling: kappa 0.78 between raters, 0.83 within. Answer check 0.947 vs experts; milestone F1 0.958; trace AUROC 0.886 with cluster-robust intervals and a stated detectable difference (RESULTS_X1 Findings 1, 2, 6). | Nothing on the slice. The paper must report the ceiling and the power, not only the point estimates. |
| Single-run point scores, no variance or significance (yAYU 2, 9W1B 4) | Nothing yet: the new run has not happened. | Per-template SD column, template-level bootstrap CIs, McNemar and Wilcoxon paired tests, decoding-variance repeat (section 3.1). Plan it before the run. |
| Thresholds 2% / 0.7 / 0.8 never varied (yAYU 4, 9W1B 3) | The new stack has different thresholds and each one is already justified: E3's 0.5% sits on a plateau against a null baseline (RESULTS_E3 grid), E4 uses the displayed precision (D-097, D-101), E2's 0.5 is calibrated held-out (D-100), the answer check's tolerance is fitted on one half and reported on the other (D-098). | Report all four grids in one appendix table. Re-run the answer-check tolerance split on the full run. |
| Step scoring coarse, 1.0 / 0.5 / 0.0 (9W1B 3) | E3/E5 continuous milestone coverage; E4 per-claim arithmetic. | Nothing further. |
| Inter-judge agreement never reported (9W1B) | E0-3J kappa 0.725 (4-way) / 0.776 (binary); E1 0.563 / 0.622 (RESULTS_E1). Screen AC1 0.929 (pass 2). | Report both. |
| Certification IAA weak, perfect kappa unconvincing (gFWV 1) | Layer 2: 15 own-branch experts, 60 of 60 planted defects caught, hand checks 459 of 504, AC1 0.914, 53 rejections over three rounds, 25 templates fixed, all 150 certified unanimously (D-106 to D-109). | D-108: 65 templates were judged through a Markdown renderer that dropped `*` and `$` (section 2, item 1). |
| Instance-level verification shallow (gFWV 3) | Layer 0 gate at 500 seeds: closure, determinism, contract, emission; 54 templates edited to pass (D-092). | Tie census on the frozen pool (section 2, item 4). |
| Template authorship ambiguous (gFWV 5) | The original 90 were expert-authored. | The 60 new templates were AI-agent-drafted under `docs/pilot_template_authoring_spec.md` (commit 63bb6ac) and then expert-certified. The paper must say so (section 4.3). |
| Only three branches (ynoK 1, AC) | Five branches, 150 templates, 2,250 items. | Taxonomy figure, Appendix C textbooks, Appendix L statistics for five branches. The new domains came from a working ABET/ASCE/IISE taxonomy, not the paper's four-model panel (Appendix A/B); either run the panel procedure on the two branches (cheap) or describe what was actually done. |
| 90 templates are 90 problems; linguistic diversity; template exploitation (nWW3 2, cqGs 2, 9W1B 1, ynoK 2, AC) | Structural diversity argument, Appendix G. Shortcut audit exists for four templates (D-046, D-057, D-066). | A paraphrase robustness experiment and a corpus-wide surface-shortcut audit (section 6). |
| Tool-free setting conflates arithmetic with reasoning; no RAG or tool baseline (cqGs 1, 4; ynoK 3; AC) | Error taxonomy showing frontier failures are arithmetic. | A tool-use condition on a subset (section 5). RAG stays future work, with the reason stated. |
| Error analysis coarse; insights descriptive (nWW3 3, 4; cqGs 3) | Six-category human taxonomy on 2,200 traces (May version). | Redo on the new run (section 3.6). The pilot adds a genuine finding: what verification can and cannot catch (section 4.2). |
| Math-pretraining claim confounded by scale and family (9W1B 5) | The Qwen2.5-Math-7B vs Qwen2.5-7B pair in the old Table 1. | Keep that pair in the roster (local, $0) and report branch-level FAC for it. |
| Validation buried in appendices (9W1B 6) | — | Validation paragraph in section 4 of the paper; demote BERTScore and ROUGE (section 4.1). |
| Datasets and software unavailable (9W1B) | — | Anonymised code and data archive with the frozen pool manifest; README is stale (still says 90 templates, three branches). |
| DeepSeek V4 Pro grouping, Claude non-monotonic scaling, `??` ref, `**` typos | — | Moot if the roster changes; proofreading pass regardless. |

The Suggested Actions document's decision criterion was: do E3/E4/E5 beat E0 and E1 on
agreement with the experts, with non-overlapping intervals? The measured answer is no at
the trace level (E5 0.886 vs E0 0.850, intervals overlap, detectable difference 0.19) and yes
where it matters for the title: E5 is the most accurate milestone evaluator (F1 0.958 vs
experts), the corrected answer check is the largest single improvement (0.947 vs 0.747), and
the deterministic digit rule is the only signal on flawed-reasoning-behind-a-correct-answer
(AUROC 0.661 vs 0.542). D-105 took the "instrument broadly" branch: full run on the
deterministic stack with a residual judge, no E0 at scale. That is the right call and the
paper should present it exactly that way: not "our evaluator beats the old one" but "the old
one was measured against experts, found wrong on a quarter of answers, and replaced by one
whose every component is validated against the same experts".

---

## 2. Before inference: the checklist

Ordered. Items 1 to 9 cost nothing in API spend. Nothing paid starts without approval
(D-093's rule).

1. **Close D-108, the Markdown rendering.** The review app rendered questions and solutions
   as Markdown, so in 65 templates the experts saw text with `*` and `$` removed
   (`layer2/markdown_scan.md`). The certification claim is "experts approved the instance
   the models will read"; for 65 templates that is not literally true. Two options: a short
   round 4 where each branch's three experts re-confirm their templates' instances in plain
   text (about 13 templates per branch, no plants, an hour each), or a stated limitation.
   Recommendation: the round. The experts have shown they read (60 of 60 plants), and
   "re-confirmed as plain text" is a sentence a reviewer cannot argue with. Fix the
   `signal_operations` origin marker at the same time.
2. **Close D-094's residuals**: `work_isothermal_virial`'s P2/Pc range (50% redraw) and
   `floating_object_submersion_depth`'s 41% stability redraw. A high redraw rate narrows the
   sampled space and shows up in T6. Decide, edit, re-gate. Any edit that moves the seeds the
   experts saw (2101 to 2105) voids that template's certification (`certification.py`), so
   fold these into the same round 4 as item 1.
3. **Re-screen the 25 templates edited after screen pass 2** (D-106): the panel's
   false-positive rate against the experts (10.2%) and the MAD table must be on the final
   corpus. Targeted re-judge, a few cents, approval needed.
4. **Regenerate the pool, then census it.** `testset/` on disk is untracked and dated
   2026-09-17; templates were edited on 09-24, 09-25 and 09-26 and 77 templates moved
   (`layer0/item_pool_impact.md`). It is stale. Run `generate_testset.py` on the final
   templates, then run `tie_census.py` on the 2,250 items actually in the pool (D-092 said
   this would be done at freeze): 12 templates still carry ties, `server_configuration_selection`
   at 15.4% of instances. A tied instance has no gold value a decimal and a binary reader
   agree on. Redraw tied instances by index, do not edit the template. Also resolve the 38
   duplicate questions the generator reports (`heat_of_reaction_formation` draws from a
   4-entry table): a duplicate double-weights an item. Redraw or accept and say so.
5. **Freeze and tag.** Commit the pool manifest with a SHA-256 per item, the master seed,
   the template commit hash and the per-record `unit` field, the way `freeze.py` does for the
   pilot. Tag the commit. Milestones for E3/E4/E5 are derived by regenerating each item from
   the templates (`milestones.py`), so the evaluation must run at the same commit the pool
   was generated at, or the D-102 drift recurs at 2,250 items. A tag makes that checkable.
6. **Validate the evaluator stack on gold, all 150 templates.** E3, E4 and the answer check
   have only ever run on the 15 pilot templates. Before any trace exists: run
   `milestones.build_all` over all 2,250 items and audit the milestone counts per template
   (templates with zero or one milestone reduce E3 to an answer check; RESULTS_E3 already
   notes `coaxial_capacitance` and `reynolds_number_flow_regime`); run the answer check gold
   against gold (must be 100% correct, all six answer kinds); run E4 on the gold solutions
   (must flag nothing, as the 227-claim validation did); run E3 on gold (coverage 1.0) and
   against a sibling item (the null must stay near 0.04). This is the single largest free
   check and it will find template-specific parsing gaps, because the 60 new templates were
   never exercised by any evaluator.
7. **Port the inference harness.** `evaluation/run_inference.py` is the January script: one
   provider, opens its output with `"w"`, no resume, no empty-completion guard, 4,096-token
   cap. The pilot's `run_traces.py` has resume, T1 to T7 verification, per-model token
   ceilings, served-model-id and prompt-hash recording, and the reasoning-token cost lesson
   (GPT-5 cost 4.8x the estimate). Generalise it to the 2,250-item pool and the roster; keep
   the Appendix P prompt byte-identical and hash it. Decide the coverage-hole policy up front:
   a truncated or empty completion is reported as unusable, not averaged (FINDINGS R-F2:
   Qwen3.8-27B never answered 4 items at 32,768 tokens).
8. **Decide the roster against the claims it must support** (section 7). The pricing
   document's 12 models cost $417.82 for inference; with E5's ~$79 that is $497 against a
   ~$500 round, leaving nothing for the router ($86), the paraphrase test, the tool condition
   or decoding repeats. This is the decision that gates everything else.
9. **Write the analysis plan before the run.** Which claims, which comparisons, which tests,
   with template-level clustering, and the smallest effect 150 templates can detect. D-105
   says to settle this before generation; the pilot's `cluster_bootstrap.py` is the tool.
   Record it as a D-entry so the paper can say the analysis was fixed before the data.
10. **Recruit the experts now for the post-run work** (section 3.6 and section 6): the
    error-analysis sample and the paraphrase spot-check. Expert time is the long pole; the
    Suggested Actions document said so and it was true in September.
11. **Estimate, then ask.** Dry-run the harness (`--dry-run`, `--check`), report the cost
    with its basis including the Gemini thinking-token floor, then request approval for the
    run.

Not blocking but do before submission: fix `README.md` (still 90 templates, three branches,
1,350 cases), commit `figures_oct_12/` (untracked, holds the five-branch overview figure and
its generator), and decide on the anonymised truth file (RESULTS_X1 Limitations: codes and
majority labels, no per-rater rows) so X1 can be reproduced from a clone. That last one is a
compliance decision for the annotators' side; the per-annotator files stay out of git either
way.

---

## 3. What the results section must contain

### 3.1 Variance and significance, as promised in July

- A per-template standard deviation column beside every headline score, computed over the 15
  instances (the rebuttal quoted 20.1 pp on FAC and 15.6 pp on F1 for the old run and
  promised the column).
- Bootstrap intervals that resample templates, not items, for every aggregate: the pilot's
  design effect was 3.7 to 4.2 (RESULTS_X1 Finding 6), and at 150 templates the honest unit
  is still the template.
- The complexity cliff with template-level bootstrap on 34 advanced templates (the old run
  had 22 and the intervals included zero for three of four models). Report the gap with its
  interval per model, not "significant / not".
- Paired tests for any decoupling claim: McNemar on the binary answer verdict, Wilcoxon on
  milestone coverage, as promised.
- A decoding-variance repeat: one cheap open model, three repeats over a stratified 300-item
  subsample (about $1), so the paper can say what part of a score is sampling noise at the
  chosen temperature. "No seeds" was a stated objection and this closes it for a dollar.
- Coverage holes reported per model (unusable traces), separately from accuracy.

### 3.2 The evaluator validation, in the main text

A validation paragraph in section 4 (9W1B 6 asked for exactly this) with: 300 traces, 15
experts, 3 per trace, kappa 0.78 between and 0.83 within; answer check 0.947 vs 0.747 for
the published check; milestone F1 0.958; trace AUROC 0.886 (template-level interval) against
E0's 0.850; the detectable difference; the hard case (section 4.2). Table 12's rho = 0.632
disappears, and with it the n = 100 study; if it is kept at all it is one sentence of history.

### 3.3 Threshold sensitivity, in one appendix table

Every threshold in the new stack with its grid and its justification: E3 0.5% (real / null /
separation at 2%, 1%, 0.5%, 0.2%); E4 digit rule vs 1% vs 0.1%; E2 0.5 (held-out
calibration, gain -0.006); the answer check's relative tolerance (split-half 0.876 / 0.905,
re-run on the full pool). The story is that no threshold was tuned on what it reports, and
that is a better answer than a sensitivity grid on arbitrary constants.

### 3.4 Judge overlap, discharged three ways

1. E1 vs E0 per model (RESULTS_E1 Result 1) and both panels' inter-judge kappa.
2. Leave-one-judge-out on E0-3J over the 300 pilot traces: for GPT-5 traces drop the GPT-5
   judge, for Gemini traces drop the Gemini judge, majority of the remaining two with the
   conservative tie-break, delta F1 against a placebo (drop a non-family judge) and against
   the DeepSeek and Llama traces. The per-judge votes are already parsed
   (`e1_analysis.py`), so this is a $0 analysis. It discharges the literal promise.
3. The X2 bias estimate: per generator family, mean (evaluator score minus expert score)
   for E0 and E1, against the expert labels. Also $0. A null result settles the objection as
   well as a positive one.
4. In the new stack the point is structural: the residual judge (MiMo) shares no family with
   any roster model, and 76% of milestones never reach a judge (RESULTS_E5). Say the judged
   fraction per model in the results table.

### 3.5 What to put in the main table instead of BERTScore

BERTScore ranged 0.866 to 0.882 across models whose accuracy ran 0.15 to 0.70, with Llama
joint highest (RESULTS_E0). gFWV said it cannot discriminate; the pilot proved it. Replace
the textual columns with: answer verdict (correct / partial / incorrect, three-way), milestone
coverage (E5), arithmetic consistency (E4 digit rule), judged fraction, unusable traces.
Keep ROUGE and BERTScore in an appendix if at all.

### 3.6 Error analysis on the new run

The May section 5.4 (2,200 traces, 11 models) is about a pool and a roster that no longer
exist. Two layers for the new run:

- **Automated attribution over all traces**: E4's digit rule flags calculation slips
  deterministically (three flags in four are real); the judge flags conceptual and
  unsupported steps. Validate the attribution on the pilot's 300 expert-labelled traces,
  whose step labels already carry calculation / conceptual / unsupported types
  (`annotation/guide.md`).
- **A human-annotated stratified sample** for continuity with the six-category taxonomy:
  100 failure traces per model for a handful of contrastive models, three annotators,
  the same decision hierarchy and agreement reporting as Appendices R to T. Smaller than
  2,200 and it can be, because the automated layer covers the population.

State plainly, as RESULTS_X1 does, that inside correct-answer traces almost every flawed
step is arithmetic (175 of 178) and conceptual error behind a correct answer is rare in the
wild and hard to detect.

---

## 4. Paper changes

### 4.1 Rewrite section 4 around the mechanism the title promises

The Suggested Actions diagnosis stands: the title says verifiable, the May mechanism was a
vote. The new section 4: answer check (deterministic, three-way, six answer kinds); milestone
coverage derived per instance from the template's own computed values, order-free
(D-084); per-claim arithmetic on the displayed precision (D-097, D-101); a judge only on the
residue, with the judged fraction reported. Then the validation paragraph (3.2), then the
limits (4.2). The comparator contract (Phase 4: scalar, vector, array, symbolic,
classification, multipart) belongs here too: the old parser reduced every answer to one
float and 74.7% of gold answers held more than one number (audit report section 1).

### 4.2 Add the finding about the limits of verification

This is the paper's best new material and it answers "lack of profound insights" better than
more tables would. From RESULTS_X1 Findings 5, 7 and 8 and D-102 to D-104:

- The final answer nearly determines the expert's trace verdict (AUROC 0.974).
- Arithmetic flaws behind a correct answer are deterministically detectable: the digit rule
  scores 45 of 45 on planted slips inside a parseable claim and 0 of 15 outside one; three
  flags in four are real on the experts' labels; it costs nothing.
- No deterministic evaluator detects a conceptual defect behind a correct answer (0 of 60
  planted). Judges catch about a third (GPT-5 0.333, MiMo 0.314), with no false alarms in
  240 clean steps. The published Tribunal's routing showed the judge the corrupted step
  exactly as often as the clean one (0.500 vs 0.500), so end to end it caught 18% of them.
- An off-the-shelf math PRM ranks steps well overall (72B AUROC 0.825) but three of its four
  flags inside correct-answer traces are false; VersaPRM finds 7% of incorrect steps;
  thresholds are not the cause (D-100).

This positions the paper in the process-supervision literature (ProcessBench, Math-Shepherd,
PRM800K, the self-preference results) and Related Work needs a paragraph on PRMs, LLM-judge
validity and meta-evaluation, which the May version lacks. If the router is funded, its
measured contribution (0.175 to about 0.31 on conceptual defects, $86) goes here.

### 4.3 Say what changed, and be exact about authorship

- **Authorship.** The paper says "the 90 symbolic Python templates were fully authored by
  domain experts", added in response to gFWV 5. The 60 civil and industrial templates were
  drafted by AI agents under a textbook-grounded authoring spec with multiple blind AI review
  cycles (`docs/pilot_template_authoring_spec.md`, commit 63bb6ac), then passed the
  deterministic gate, the non-suite screen and three rounds of own-branch expert certification
  that rejected 22 of them in round 1 and required 25 fixes. That is a defensible pipeline
  and it must be described as what it was. Writing "authored by domain experts" over 150
  templates would be false, and the reviewer who raised authorship once will look again.
- **The taxonomy for the new branches.** Section 3.1 and Appendices A, B and D describe a
  four-model panel with majority voting and a pedagogical significance score. The new
  domains came from the authoring spec's working taxonomy (ABET, ASCE, IISE). Either run the
  Appendix A/B/D prompts on the two branches now (minutes, negligible cost) and report
  agreement with the chosen domains, or describe the procedure that was used.
- **Certification, honestly.** The December pipeline approved 270 of 270 rows at 12 to 30
  seconds per template (D-095). The new one has hand checks before the solution is shown,
  planted defects, timestamps, three rounds, 53 rejections, and the experts found real physics
  problems (a nylon shaft at 743 MPa, honey needing 2.3 GPa of pressure drop). Report the
  rejections as a strength. Appendix K becomes the Layer 0 / 1 / 2 description with the
  plant-detection rate, hand-check agreement, AC1 and kappa, and the panel's false-positive
  rate against the experts.
- **Why the numbers moved.** The May Table 1 came from an answer check that understated
  every model by about 21 points and ranked GPT-5 fourth where the experts rank it first
  (D-098, D-105). The gold traces of eight published templates did not reproduce their own
  answers (Phase 1). The published Tribunal was two judges, not three (E0-F6), and a 20%
  sampling step reordered the ranking between two runs of identical code (E0-F7). None of
  this needs to be a confession in the abstract, but the framework validation section should
  state that the earlier check was measured against experts and found wrong, because a reader
  comparing versions will see the jump. Do not compare new numbers with old ones: the old
  pool is not reproducible (`inference_results/` is gone, phase6_item_pool_impact section 4).
- **The roster changed and shrank.** 27 models became about 12, and the flagship tier is out
  by the supervisor's steer. Frame the roster by capability tier and by what each claim
  needs, not as a regression from 27. Keep the anchors that claims depend on (section 7).
- **Scope sentences.** The abstract's "stress-test generalization across diverse physical
  scenarios" becomes physical and structural diversity, as promised to 9W1B 1; the
  safety-critical framing gets the promised distinction; the Limitations section is rewritten
  for five branches and for what the paraphrase test shows.

### 4.4 Appendix housekeeping

Figure 2 taxonomy for five branches and 15 domains (figures_oct_12 has the overview figure);
Table 2 textbooks for civil and industrial (from `docs/references/`); Table 9 dataset
statistics for 150 templates (870 / 870 / 510 items); Appendix O for the new roster with
served model ids as the harness records them; Appendix P with the inference config actually
used (token ceilings per model, not 4,096); the pinned scoring-library versions (D-081: the
BERTScore column silently reads zero under transformers 5.x); a note that the published
Gemini judge no longer exists (E0-F4), which is the reproducibility argument for an
open-weights judge.

---

## 5. RAG baseline: no. Tool-use condition: yes, on a subset

Every reviewer in both cycles and the January AC asked for tool or retrieval baselines, and
both rebuttals answered "future work". A third "future work" will read as a pattern.

**RAG: do not add it.** A retrieval baseline over textbooks needs a corpus the paper cannot
release, a retrieval design that becomes its own attack surface (which chunks, which
reranker, which k), and it would be answering the wrong question: the error taxonomy shows
formula recall is a minor failure mode at the frontier (Formula/Principle errors 0 to 3.5%
for frontier models, Table 18), and Hallucination of constants 1.5 to 3.5%. Say this in one
paragraph with those numbers, name RAG as future work with the reason, and move on.

**Tool use: add a code-execution condition on a subset.** It tests the paper's own central
claim, that frontier failures are arithmetic execution. If accuracy rises sharply with a
calculator, the diagnosis is confirmed and the paper has a second finding; if it does not,
that is a finding too. Design: a Python tool the model may call, same prompt otherwise, on
three roster models spanning the tiers (or all models on a stratified 750-item subsample),
scored with the same stack. E4's arithmetic check will flag less because the arithmetic is
delegated, and that difference is itself the measurement. Cost is roughly 1.5 to 2 times the
inference cost of the models chosen, so about $100 to $150 for three models on the full pool,
less on a subsample. It needs a roster or budget decision (section 7).

---

## 6. Linguistic diversity: run the paraphrase test, and the shortcut audit

The January AC's third suggested revision was paraphrased or naturally sourced formulations.
The team has answered "structural diversity" twice. The structural argument is right and
insufficient, because the objection is empirical: can a model exploit the fixed wording?
That is testable for very little money.

**Paraphrase robustness experiment.**

1. Generate one paraphrase of each question in the pool (2,250) with a model that is neither
   on the roster nor a judge, under a constrained instruction: reword the prose, keep every
   number, unit, symbol and technical term unchanged, change nothing about the scenario.
   Cost: dollars.
2. Verify deterministically that the multiset of numbers and units in the paraphrase equals
   the original's, and reject and regenerate otherwise. Gold solutions are untouched, so
   every evaluator runs unchanged.
3. Expert spot-check a stratified sample (two instances per template, 300 items, split across
   the 15 experts: 20 each) for "same problem, same answer". The experts' own hand check
   format works for this.
4. Run the same models as the tool condition on the paraphrased pool and report the delta in
   answer accuracy and milestone coverage per model, with template-level intervals. A stable
   score is the direct answer to "template-level pattern exploitation"; a drop is a finding
   about the models, not a flaw in the benchmark, and either way the paper reports it.

**Surface-shortcut audit, corpus-wide.** D-057 measured two classification templates as
100% predictable from the question surface, D-046 one at about 90%, and D-066 a blind-guess
floor of 1.0; four templates are excluded from the pilot for it. Extend the audit to all 150
(the depth-2 held-out rule fit from D-057 is the method), report the blind-guess floor and
the surface-model lift per template in an appendix, and exclude or flag the shortcuttable
ones in the headline numbers. This is the rigorous version of "instances are not just
different numbers".

**Wording.** The abstract and section 3.2 keep the distinction between physical, structural
and linguistic diversity, and the Limitations section says what the paraphrase test did and
did not show. Multi-modal and real-artifact seeds stay future work.

---

## 7. Budget and the decisions only you can make

| Item | Cost | Status |
|---|---|---|
| Spent so far this round (pilot, screens, planted-defect judges) | about $47 | done |
| Inference, 12-model roster (pricing doc; Gemini figures are floors) | $417.82 | needs roster decision and approval |
| E5 residual judge over 27,000 traces (D-105) | about $79 | approved in principle by D-105 |
| Step router (D-105, `router_residue.py`) | $86 | open |
| Paraphrase generation plus paraphrased-pool inference for 3 models | about $5 plus roughly $80 to $130 | open |
| Tool-use condition, 3 models | about $100 to $150 | open |
| Decoding repeats, LOJO, X2 bias, shortcut audit, gold validation | $0 to $1 | open, all free |
| Post-run expert error-analysis sample, paraphrase spot-check | expert time | recruit now |

The base plan alone reaches the ~$500 round. Three ways through, for the supervisor:

1. Trim the roster. Kimi K3 is $130 of the $418 and carries the heaviest documented
   Claude-distillation exposure of any candidate (JUDGE_SELECTION); GLM-5.3 is $53. Dropping
   both funds the paraphrase test and the tool condition on three models.
2. Ask for the increase with this table: each optional line answers an objection raised in
   both cycles.
3. Run the tool and paraphrase conditions on a 750-item stratified subsample for all models
   instead of the full pool for three; same money, wider coverage, wider intervals.

Roster anchors the claims depend on, whichever way the budget goes:

- Qwen2.5-Math-7B and Qwen2.5-7B, the same-family same-scale pair, run locally on
  HiPerGator at $0. It is the only controlled answer to 9W1B 5, and it is free. E2 is
  secondary (D-105), so the Qwen-as-PRM block in `roster_candidates.py` does not bar it.
- Llama 3.1 70B: under a dollar for 2,250 items, and it is the model the expert labels and
  the pilot anchor on.
- At least two closed models from different families so "capability tier" means something;
  the pricing roster has four.
- No Xiaomi, MiniMax or xAI model, because they judge (D-088).

Decisions to record as D-entries once made: D-108 (round 4 or limitation), D-094 residuals,
the roster and its rationale, the router, the paraphrase and tool conditions, the
analysis plan, the anonymised truth file.

---

## 8. Risks

- **Expert availability.** Round 4 (item 2.1), the paraphrase spot-check and the error
  analysis all need the 15 experts again. Batch the asks so each expert is approached once.
- **Milestone derivation on the 60 new templates.** E3 was built on 15 templates. Templates
  whose gold states one intermediate reduce E3 to an answer check; templates with iteration
  or search shapes (Phase 3 node types) need the trajectory exclusions checked. Item 2.6 is
  where this surfaces, before money is spent.
- **Reasoning-model cost and truncation.** Every reasoning model in the pilot needed 32,768
  tokens and still left unusable traces; GPT-5 cost 4.8x the estimate. Budget from measured
  token profiles (`inference_cost.py`), not list prices, and decide the unusable-trace policy
  before the run.
- **A reviewer reads the September record.** The repository is public and the DECISIONS log
  says, correctly, that the published framework had eight defects and eight published
  templates had wrong gold. That is a strength only if the paper says it first. A paper that
  presents the new numbers without explaining the old ones invites the reader to find the
  explanation in the repo.
- **Power.** 150 templates support model-level claims; they do not support fine evaluator
  rankings, and the pilot's 15 could not either. Report the detectable difference and let
  every null result say "no difference this large".
