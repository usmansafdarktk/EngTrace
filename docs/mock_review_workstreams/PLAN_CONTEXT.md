# Mock-review fix plan, October 2026 submission

For the authors. Written 2026-10-06, after the two mock ARR reviews in `ACL Mock Reviews/first cycle/`
(`ARR_Review_EngTrace_1.md`, `ARR_Review_EngTrace_2.md`). Deadline: 12 October 2026 (ARR, anywhere on Earth).
This plan amends `docs/PAPER_PLAN_OCT2026.md`; that plan's writing rules (section 8) and
`overleaf_source_04102026/WRITING_RULES.md` still govern every sentence. Nothing below has been done yet.

Contents: 1 why; 2 rules; 3 decisions the authors must make; 4 what the experts are asked for; 5 work packages
with their steps; 6 schedule; 7 out of scope, with the Limitations sentence for each; 8 traceability from every
review item to its work package.

---

## 1. Why these changes

### 1.1 What the two reviews say

Two AI models reviewed the LaTeX sources with the ARR form and a strict brief. Review 1 recommends "resubmit next
cycle" (overall 2.0; soundness 2.5, excitement 2.5, reproducibility 3, confidence 4). Review 2 recommends Findings
(overall 3; soundness 3, excitement 3, reproducibility 3, confidence 4). Both read every file and re-derived
numbers; both found the descriptive numbers accurate and internally consistent. Neither saw `main.tex` or
Figures 1 and 2, which were not in the bundle, so "no Limitations section" is partly a bundling artifact. But
`current_overleaf_project/main.tex` still carries May's Conclusion, Limitations and Ethics text (27 LLMs, three
branches, a "complexity cliff"), which contradicts the new paper and must be rewritten regardless.

The objections the two reviews share, most serious first:

1. **Saturation presented as a finding.** The five strongest models lie within 0.009 on the final answer, two
   flagships do no better, and the paper's own accounting shows that most of the top five's remaining incorrect
   verdicts fall on two ambiguous chemical templates, one symbolic template scored by its numbers, and prescribed
   digits. The abstract says "final-answer accuracy alone cannot separate the strongest models" and the results
   section ends on "headroom". Both reviewers read that as spin; both want the ceiling stated plainly, the
   headline numbers reported on a cleaned set, and the defective templates repaired or removed.
2. **"Verifiable process supervision" outruns the evidence.** Conceptual step errors are not verified (no
   deterministic check catches any of 60 planted misstated rules; judges catch about a third), the step checks
   have recall below one half, and "process supervision" conventionally means a training signal. The motivating
   premise, a correct answer behind flawed reasoning, is rare on this benchmark (AUROC 0.974; 175 of the 178
   flawed steps behind correct answers are arithmetic slips).
3. **Milestone Coverage is confounded.** The final answer is itself a milestone (290 instances have exactly one);
   MC credits what a model displays; a valid alternative route is not credited (a third of the confirmed
   "missing" rulings); the one top-tier separation (Claude Sonnet 5 against DeepSeek V4.1 Flash, 0.025) is of the
   size of judge noise (the judge decides 13% to 17% of their milestones at precision 0.85; re-judging moved
   coverage by up to 0.020); Kendall's tau 0.709 is moderate agreement, not a different ranking.
4. **Unequal inference conditions.** Four models return no reasoning tokens at their providers' defaults, output
   ceilings differ, open-weights models are served by up to eleven providers, and "four of the top five are
   open-weights" is a configuration artifact. Reasoning effort alone moves GPT-5.4 mini by +0.093.
5. **Known-defective templates stay in the scored set**, and the planted defects (transposed digits, factors of
   100, flipped signs) would not have caught the under-specification that got through. 36 templates show 6 to 15
   reasoning paths, yet each expert hand-solved one instance.
6. **The evaluator was validated mostly on the data used to design it**, on 15 templates of which 3 have a scalar
   answer; agreement on the evaluated models is lower (85 of 100 "reached" rulings, 83 of 98 incorrect verdicts,
   34 of 52 partial verdicts called fully correct), and no per-model scorer accuracy is reported.
7. **Difficulty labels have no protocol or agreement**, the gap survives correction for none of the models once
   the two chemical templates are removed, the showcased Advanced listing is a two-formula problem, and the
   showcased Easy question prints its formula.
8. **The taxonomy rests on LLM votes** and the "curricular standards" cited for AIChE, ASME and IEEE are a
   constitution, a mission statement and a standards homepage; the cut from the 5 to 7 prompted domains to 3 per
   branch is unexplained; core areas are missing (circuit analysis; thermodynamics and heat transfer in
   mechanical; transportation and environmental in civil).
9. **Unsupported or inconsistent statements.** The single-path consistency claim in §5.3.1 cites an appendix that
   does not hold it; the four re-decoded models and their 300 instances are unnamed; "almost no template defeats
   a model on all 15 instances" is unquantified; "Three experiments" against the appendix's four; "0.04 to 0.05"
   against a spread of 0.051; the unit-of-analysis justification ("15 instances share one derivation")
   contradicts the 92 multi-path templates; 504 against 510 hand checks; 51 of 60 for one judge; kappa 1.000
   where it is undefined.
10. **Procedural gaps.** Nothing on the experts (qualifications, recruitment, pay, consent, independence from the
    template authors), no release statement or licence in the bundle, no query dates, no Ethics section.
11. **Novelty is incremental and the related work has gaps**: dynamic and functional benchmarks (DyVal, Functional
    Benchmarks, MathCAMPS, GSM-Plus), reference-free step evaluators (ROSCOE, ReCEval, ReasonEval), science
    benchmarks (SciBench, TheoremQA, OlympiadBench, PHYBench, UGPhysics), PRMBench; no ChainEval-style baseline
    is named; the delta from ThermoQA and FinChain is not stated precisely.
12. **Clarity.** Sentences carrying four to six statistics; no worked scoring example; the defective templates
    are mentioned three times in the main text and explained only in the appendix; Table 1 bolds a "best" inside
    a tie; the branch, radar and error-category figures are hard to read; no Conclusion.

### 1.2 What the record says (checked 2026-10-06 against the scripts and result files)

| Objection | What the code and results hold |
|---|---|
| Single-path claim has no table | `full_run_28092026/results/results.json` key `q4` and `RESULTS.md` Q4 hold the numbers; `paper_results.py` reads them for the prose (lines 363, 586, 676) and writes no table. |
| MC by matching alone, all responses | `RESULTS.md`, "Milestone coverage over all traces" (D-170): E3 and E5-strict per model with template intervals, all and readable. Only the wrong-answer version is printed (`tab:coverage`). |
| Per-provider effects | `RESULTS.md`, "By serving endpoint" (D-133): raw score, unusable share and a template-matched difference per endpoint. Not in the paper. |
| Empty responses per model | `RESULTS.md` Q1, "unusable (empty + unreadable)" per model. The paper marks only GLM-5.3's 100. |
| Readable-only reordering | `RESULTS.md` Sensitivity, tau 0.636; in the appendix's scoring-rules paragraph, not in "Stability of the Tiers". |
| Is reasoning text searched? | `score.py` scores `row['text']`; `run_traces.py` stores the endpoint's reasoning text in a separate field. Matching and the arithmetic check read the visible response only. The paper does not say so. |
| Can MC variants be computed offline? | Yes. E3 rows carry a per-milestone `reached` list; `judge.py` rows carry a per-milestone source (`e3`, `REACHED`, `NOT_NEEDED`, `MISSING`, `UNJUDGED`) and `e5_lenient`; answer rows carry the numeric targets and `matched` / `of` counts. |
| Were the two chemical templates fixed? | Partly, and not for this defect. Round 4 (2026-09-28, D-115, `layer2/fixes_round4.md`) widened `adiabatic_flame_temperature` (samples excess air) and `heat_of_reaction_formation` (25 reactions) so each yields 15 distinct questions, and the chemical experts re-certified both before the freeze. `work_isothermal_virial` was approved in round 2 unchanged (its round-1 rejection concerned high temperatures, not adopted, D-094). The ambiguity the paper reports (wording decides neither the closed-system nor the flow reading nor the virial form; flame data sources differ) was found by the B4 reading on 2026-10-03 (D-171, D-185), after the run; the decision then was a stated limitation and a sensitivity row, no re-score. No template file changed after 2026-09-28 (`git log`). The repair the reviewers ask for is open. |
| Cost of reasoning on | The `reasoning-medium` arm on the 450-instance subset cost $9.77 with judge and router (D-180); inference alone $4.44 (GPT-5.4 mini) and $1.10 (Gemini 3.1 Flash-Lite). Both endpoints returned reasoning tokens on every row. Gemma 4 26B was never tried with the parameter. Qwen3-235B-2507 is the Instruct variant and has no thinking mode; its Thinking sibling is a different model. |
| 504 against 510 hand checks | `layer2/RESULTS.md` line 52: "504 hand checks with a comparable number"; six items had no number to compare. |
| 51 of 60 for MiMo | Pilot `RESULTS_X1.md`: 17 of its 240 calls never returned after retries. |
| Kappa 1.000 for civil and industrial | All 90 verdicts approve, so Fleiss' kappa is 0/0; `layer2/RESULTS.md` prints 1.000. |
| Related-work gaps | `docs/related_work_oct2026/related_work_sources.bib` already holds ROSCOE, ReasonEval, Functional Benchmarks (Srivastava), SciBench, OlympiadBench, PHYBench, UGPhysics, PRMBench and GSM-Plus; DyVal, ReCEval, TheoremQA and MathCAMPS are not in it. The "step matching with in-family judges" design in the validation appendix is the FinChain-style step-alignment baseline in all but name; confirm against the pilot's design record before saying so. |
| Letter display for tiers | `paper_results.py` already has `letters()` (line 259); Table 1 uses bold instead. |

### 1.3 Goals, in order

1. Every claim matches its evidence: title, abstract, contributions and the results prose say what is verified
   and what is not, state the ceiling, and keep the MC separation only if it survives its own controls.
2. The configuration confound is removed: every model that can reason is also evaluated with reasoning on, on the
   full set, and the tiers are stated under both configurations.
3. Every analysis the record holds and a reviewer will ask for is printed, and the two MC variants and the
   scoring-clause counts are added from the stored scores.
4. The procedural gaps are closed: Conclusion, Limitations, Ethics, the expert paragraph, the release statement.
5. The template defects are repaired where the experts can turn it around in time; otherwise the cleaned set is
   reported beside the full set and the defects are named in the main text.

---

## 2. Rules for this work

- The paper states the final state only. Round 4, the B4 reading, the repair and the re-run are history: they go
  to `docs/REVISION_LETTER_NOTES.md`, not the paper.
- Every number comes from a committed script whose `--check` passes: `full_run_28092026/paper_results.py`,
  `paper_setup.py`, `docs/appendix_statistics.py`, `docs/appendix_certification.py`, `docs/appendix_evaluation.py`,
  `docs/appendix_listings.py`. No typed numbers, no numbers copied from Markdown tables.
- Paid runs: a dry run with an estimate first; nothing is billed without the owner's explicit approval of that
  run; the bill is recorded with the run.
- The experts' filled files stay local and gitignored; the paper reports counts and agreement.
- Commits only when the owner says; each commit pushed at once; one-line messages.
- Figures: draw, render, show, wait for approval; no tex edits and no other figure touched until then.
- Decisions go into `docs/re-implementation-sep/DECISIONS.md` (next numbers D-197 onward) when the owner asks.
- No internal vocabulary in the paper (PAPER_PLAN section 8): the proofreading pass greps for it.

---

## 3. Decisions the authors must make

Each with a recommendation. The schedule assumes they are made on the morning of day 1.

| # | Decision | Recommendation |
|---|---|---|
| D1 | **Title and the term.** Keep "Verifiable Process Supervision" or switch to the plan's alternative, "Verifiable Process Evaluation of Engineering Reasoning from Symbolic Templates". | Switch. Both reviews object that "supervision" names a training signal, and the paper runs no supervision experiment. Use "process evaluation" and "step-level checks" throughout. |
| D2 | **Headline configuration.** Keep the pre-specified evaluation at provider defaults as the headline, with the reasoning-on configuration beside it; or make reasoning-on the headline for the models that support it. | Keep the defaults as the headline (the expert readings, the error analysis, the paraphrase test and the repeats are all on that run), add the reasoning-on rows to Table 1 marked as such, and state the tiers under both. The matched configuration gets its own pairwise family in the appendix. |
| D3 | **Which models to re-run with reasoning on.** GPT-5.4 mini and Gemini 3.1 Flash-Lite (the parameter works); Gemma 4 26B (unknown; calibrate 20 instances first); Qwen3-235B-A22B-Thinking-2507 as a twelfth model (a sibling model, not a setting; eligible under the roster rule if it is neither a pilot generator nor a judge, to be checked in `docs/inference_pricing/build_pricing_doc.js`). | Run the two closed models on all 2,250 instances; calibrate Gemma 4 and run it if its endpoint returns reasoning tokens; add the Qwen Thinking model only if the roster check passes and the dry run is cheap, as a labelled extra row. |
| D4 | **The two chemical templates.** Path B: repair the wording and data, re-certify, re-run their 30 instances, regenerate. Path A: keep them, report the cleaned set beside the full set, name them in the main text. | Path B if the three chemical experts can return the re-certification within 48 hours and the readers can top up the error analysis within the window; otherwise Path A. Decide on day 1 after asking the experts. |
| D5 | **Symbolic answers.** Grade the roughly 90 incorrect symbolic verdicts by expert, then (a) implement symbolic equivalence for `ber_estimation_mary` and re-score, or (b) report the expert-adjudicated score as a sensitivity. | Do the expert grading now (it is needed either way); attempt (a) for `ber_estimation_mary`, validated against the grades; fall back to (b) for the rest. |
| D6 | **Prescribed digits (17 templates).** Keep the exact-digit requirement as the construct, or relax to the tolerance. | Keep, and add the relaxed score as a sensitivity row with one sentence in the main text; state the policy in the scoring appendix. |
| D7 | **Difficulty labels.** After the experts' independent ratings: report agreement only, or adopt the majority label where it differs. | Report agreement and keep the labels if the majority differs on a handful; adopt the majority and regenerate if it differs widely. Decide when the ratings return. |
| D8 | **Release and contamination policy.** What is released (templates, evaluator code, all responses and scores, judge rulings, hashes), the licence (main.tex says MIT), and the plan after the seed is public. | Release all of it under MIT; state that the generator draws fresh sets from a new seed and that the authors keep a held-out seed for a hosted or later evaluation; withhold names and the experts' per-rater files. |
| D9 | **Figure labels.** Figures say "GPT OSS 20B", the text `gpt-oss-20b`. | Keep the figure label (the authors' choice) and add "GPT OSS 20B is `gpt-oss-20b`" once in each caption, or switch the figures; either way, consistent. |
| D10 | **Taxonomy facts** (not a decision but needed): the rule that cut the prompted 5 to 7 domains to 3 per branch, and whether the LLMs "cross-checked" the authors' lists (main text) or "finalized" them by vote (appendix). | State the actual process in both places in the same words. |

---

## 4. What the experts are asked for

Nothing new from the experts is strictly required for the largest fixes (back matter, reframing, the printed
analyses, the MC variants, the reasoning re-run: the evaluator's validation carries over). Four asks improve the
answers to specific objections; three are small. Send them on day 1 so that they return by day 3.

| # | Ask | Who | Effort | What it unlocks |
|---|---|---|---|---|
| E1 | **Independent difficulty ratings.** Each expert assigns Easy, Intermediate or Advanced to the 30 templates of their branch, with the three dimensions of §3.1 as the rubric, without seeing the current labels (one instance per template in the kit). | all 15 experts | 30 to 45 minutes each | Fleiss' kappa on the levels per branch and overall, and a one-sentence protocol in §3.1. Answers the "no labelling protocol or agreement" objection in both reviews. |
| E2 | **Symbolic-answer grading.** For every incorrect verdict on the nine symbolic templates (about 90 across the eleven models, `RESIDUAL_INCORRECT.md`, "symbolic" column), the gold expression beside the response's final expression: equivalent, not equivalent, or unreadable. | one expert per branch with symbolic templates (mostly electrical) | one to two hours | Either a validated symbolic check for `ber_estimation_mary` (34 of the top five's incorrect verdicts) or an expert-adjudicated sensitivity row; the abstract's "six kinds" sentence can then be stated exactly. |
| E3 | **Path B only.** (a) Re-certify the repaired `work_isothermal_virial` and `adiabatic_flame_temperature`: hand-solve one instance, read the solution, approve or reject, as in rounds 1 to 4. (b) Read the wrong answers that the re-run adds to the error-analysis sample, with the same six-category hierarchy: mostly Claude Sonnet 5's (up to 25 of its 40 read items fall on these templates), a few for the other three models. | (a) the three chemical experts; (b) the three readers of each affected branch | (a) about 20 minutes each; (b) roughly 25 to 40 items, three readers | (a) "certified by a unanimous vote" stays true after the repair; (b) the error analysis stays a 40-per-model sample. |
| E4 | **Facts, not work.** Qualifications and roles of the 15 experts; recruitment; compensation; consent; whether the template authors certified their own templates; whether the same 15 labelled the expert study, vetted the paraphrases and read the error analysis; whether the three chemical experts who called the two templates ambiguous are their certifiers; how the 15 expert-study templates were chosen; the taxonomy facts of D10. | the authors | an hour | The expert paragraph of the Ethics statement, the independence sentence in §3.3, and the study-sampling sentence in the validation appendix. |

Not asked this cycle: a fresh expert-labelled validation sample on the evaluated models (per-model scorer
accuracy), and readings of the reasoning-on run's flags. Both go to the Limitations as scope and to the next
cycle.

---

## 5. Work packages

Each package names the review items it answers (R1 = review 1, R2 = review 2), where it lands in the paper, the
steps, the scripts, the checks and the effort. Dependencies on decisions (D) and expert asks (E) are marked.

### WP1. Back matter: Conclusion, Limitations, Ethics, release (R1 C1, M11; R2 W11, W16)

Where: `current_overleaf_project/main.tex` (the three sections live there, after the inputs) and a new appendix
file `overleaf_source_04102026/appendices/release.tex` (the plan's Appendix N).

1. Replace May's Conclusion with one paragraph (0.3 page): what EngTrace is, what the evaluation found (the ceiling
   on final answers for the strongest models, what the process measures add on wrong answers and behind right
   ones, what reasoning effort changes), what comes next (method memorisation, multi-modal artefacts, the
   open-weights tier under tool use), in the paper's own words and without "first".
2. Write Limitations as one paragraph of scope sentences, each a single sentence, from PAPER_PLAN section 5
   ("Limitations") plus the sentences of section 7 below. No stacked numbers, no process detail.
3. Write the Ethics statement: synthetic data from referenced constants; no private data; the experts (E4: who
   they are, how recruited, compensated, consent); their labels and names stay out of the release; models should
   not be deployed in safety-critical systems on this evidence; the licence.
4. Write the release appendix (D8): what is released and when (templates and generator, evaluator code, the
   evaluation set with its seed at publication, the 24,750 responses with scores and judge rulings, the paraphrase
   pairs' hashes), what is withheld and why, query dates of the run (28 September 2026 and the arm dates from the
   decoding tables), the script behind each table, the contamination policy after release. Pointer sentences from
   Limitations and Ethics.
5. Remove the commented acknowledgement block from the review version.

Checks: compile; the Limitations and Ethics sections are unnumbered (`\section*`) as the ACL template requires.
Effort: half a day. Depends on E4, D8.

### WP2. Framing: title, abstract, introduction, contributions, results prose (R1 F1, F2, M9e, M9f; R2 W1, W4, W5, W6)

Where: `main.tex` (title), `0_abstract.tex`, `1_intro.tex`, `5_evaluation.tex` (first paragraph), `6_results.tex`.

1. Title per D1. Replace "process supervision" everywhere with "process evaluation" (grep the sources).
2. Abstract: say concretely what is deterministic (the final value, the presence of the gold trace's intermediate
   values, displayed arithmetic) and what is judged; say that conceptual step errors are not detected
   deterministically; replace "the five strongest within 0.01 of one another" with "the five strongest are not
   separable at this size: the evaluation set is near ceiling for them on the final answer"; say that symbolic
   answers are scored by the numbers they state (unless D5(a) lands); qualify the Easy-slips sentence per WP8.
3. Introduction: reframe the motivation. Not "a correct answer can follow flawed reasoning" as the premise, but:
   on examination-style derivations the answer nearly settles soundness, and what step-level checks add is where a
   wrong derivation breaks, arithmetic slips behind right answers, and a comparison when answers tie. Keep the
   ProcessBench citation as the general observation and state the finding beside it.
4. Contribution 3: "we show that final-answer accuracy no longer separates the strongest models on these problems,
   and what the process measures add" rather than "cannot separate".
5. Results prose: drop "Four of the top five are open-weights"; rewrite the "headroom" close of §5.4 as the ceiling
   it is ("on the final answer the evaluation set is solved at the top by the models of September 2026; what
   separates them is the derivation"); make the tool and equations conclusion a bound ("neither an offered
   interpreter nor the governing equations moves the two closed models beyond ±0.05, which for Claude Sonnet 5
   is the ceiling"); say that each experiment-reading agreement in §5.4 rests on one model; keep the Claude
   against DeepSeek sentence only under WP5's rule.
6. Say once, in §4, that matching reads the visible response and not the reasoning text (WP4.7).

Checks: `docs/check_plan_claims.py` if it still applies; the "do not state" list of PAPER_PLAN 3.2. Effort: one day
of careful text. Depends on D1, WP5.

### WP3. Re-inference with reasoning on, full set (R1 M4, M9e; R2 W6, Q7)

Where: `full_run_28092026/run_traces.py`, `score.py`, `judge.py`, `router.py`, `decoding_table.py`, `analyze.py`,
`paper_results.py`; `6_experiments.tex`, `6_results.tex`, `appendices/models.tex`, `appendices/results.tex`,
`appendices/further_experiments.tex`.

1. **Roster check** (D3): confirm in `docs/inference_pricing/build_pricing_doc.js` that the Qwen Thinking sibling is
   neither a pilot generator nor a judge; record the identifier and price.
2. **Variant.** Add a variant family `reasoning-<effort>-full` to `run_traces.py` (`VARIANTS`, `items_of()`,
   pricing): all 2,250 pool items with their original questions, the reasoning-effort parameter at medium (the
   effort the 450-instance arm used), the main run's prompt, ceiling and routing, for the models named on the
   command line. Its `--check` asserts the item count and that the prompt hash equals the main run's.
3. **Gemma 4 calibration.** `--calibrate 20 --yes` with the parameter (cents). If the endpoint reports reasoning
   tokens on the 20 rows, Gemma joins the run; if not, Gemma stays at its default and the setup paragraph says the
   endpoint offers no reasoning setting.
4. **Dry run and approval.** `--dry-run` prices the full set from the arm's measured token statistics. Expected
   inference: about $22 (GPT-5.4 mini), $6 (Gemini 3.1 Flash-Lite), Gemma and Qwen Thinking to be priced; judge
   and router stages about $5 per model. Record the estimate; launch only on the owner's approval of that
   estimate and the cap.
5. **Launch** detached with the keep-awake script, one quiet background poller; expect a few hours per model.
6. **Score.** `score.py` on the variant (its own store with `CONFIG.json` provenance); `judge.py` and `router.py`
   stages on it (approval of their calls); `decoding_table.py --variant`; `trace_review.py`.
7. **Analysis.** In `analyze.py`, add a *matched configuration* family: the eleven models with each re-run model
   replaced by its reasoning-on rows; the 55 pairs with the sign-flip test and Holm within the family; level gaps;
   MC with and without the judge; the share of responses with no readable answer; the paired change per model
   against its default rows (the full-set version of the reasoning-effort experiment). Write it to `RESULTS.md`
   under a new heading and to `results.json`.
8. **Paper.** `paper_results.py`: Table 1 gains one marked row per re-run model ("with reasoning at medium effort")
   or a second block, per D2; a new appendix table "Matched reasoning settings" with the pairwise outcome; the
   generated tier sentences carry both configurations ("the tiers hold when the four models reason" or whatever
   the data say); the further-experiments appendix's reasoning-effort paragraph reports the full-set paired change
   and the 450-instance arm is retired or kept as the first estimate. `6_experiments.tex`: the setup paragraph
   states both configurations and which models have no reasoning setting. The Limitations sentence on unequal
   effort shrinks to what remains.
9. Say in the error-analysis section that the readings are of the default-setting responses.

Checks: `paper_results.py --check`, `paper_setup.py --check` (the setup numbers change: 24,750 becomes the main
run's count plus the re-run's; say both). Effort: inference day 1 to 2, scoring and analysis day 2 to 3, text
day 3. Depends on D2, D3, the owner's approval of each bill.

### WP4. Analyses the record already holds, printed (R1 M4, M9a, M9c, M9d, M9g, m13, M3; R2 W12, W13, W10, Q9)

Where: `paper_results.py` (new generated blocks), `appendices/results.tex`, `appendices/models.tex`,
`6_results.tex`.

1. **Single-path consistency table** from `results.json` `q4`: per model, the share of single-path templates
   solved on all, some and no instances, and the same for the 92 others, with intervals. Fixes the dangling
   pointer in §5.3.1 and quantifies "almost no template defeats a model on all 15 instances" (the "none" column).
2. **MC by matching alone for all responses** (D-170 table) as a Table 1 column or a column of `tab:coverage`, with
   the judge-decided share beside it. Run the 55-pair family on it (WP5.3).
3. **Cleaned-set headline**: FAC and MC for every model without the two chemical templates and without
   `ber_estimation_mary` (from `per_template.csv` for FAC and E3; from the judge store for E5-strict), with
   intervals, as an appendix table and one sentence in §5.3.1 (Path A) or as the headline itself (Path B, where
   the table is unnecessary).
4. **Responses with no readable answer per model** as a Table 1 column (or footnote per model), from Q1's
   "unusable"; the readable-only FAC and its reordering (tau 0.636) moved into "Stability of the Tiers" with a
   sentence that GLM-5.3's position is an output-ceiling effect.
5. **Per-provider table** in the decoding appendix from "By serving endpoint": endpoint, rows, raw score,
   template-matched difference; one sentence that the matched differences lie within ±0.02 for every endpoint
   that served more than a handful of templates (verify the number before writing it).
6. **Name the four re-decoded models and the 300 instances** (the repeat subsample: `subsamples.py`) in the
   paraphrase appendix and in §5.3.1.
7. **Visible text only**: one sentence in §4 and in the prompt appendix that matching and the arithmetic check read
   the response text and not the reasoning the endpoint returns separately; a descriptive table of MC against
   reasoning tokens per model (from the stored `reasoning_tokens`) in the coverage appendix.
8. **Scoring-clause counts** (R1 M3, Q4): re-score the main store offline with (a) the absolute-value clause off and
   (b) the last-digit term bounded, min(u(ŷ), 0.01|y|); count the correct verdicts that depend on each; report the
   distribution of relative error among accepted answers (a `residual_incorrect.py`-style table for accepted
   answers). One sentence in the scoring appendix; the sign rule stays if the count is small, with the count.
9. **Prescribed digits relaxed** (D6): a sensitivity column with the exact-digit requirement relaxed to the
   tolerance; one sentence in §5.4.
10. **Credit per part** (R1 M3): from `matched` / `of` in the answer rows, FAC with proportional partial credit as
    a sensitivity column.
11. **Depth, controlled** (R1 M6): the share scored 0 by milestone count over readable responses only, and a
    template-clustered model of wrong answer against milestone count with answer kind as a covariate; replace the
    untested "depth lowers every model" sentence with what the model supports.

Checks: `paper_results.py --check` after `--write`; every new block between its BEGIN and END markers. Effort: one
and a half days of script work.

### WP5. Milestone Coverage variants and the separation rule (R1 M1; R2 W3, Q3, Q4)

Where: a new script `full_run_28092026/coverage_variants.py` (reads the E3 and judge stores), `analyze.py`,
`paper_results.py`; `5_evaluation.tex`, `6_results.tex`, `appendices/results.tex`.

1. **Route-adjusted MC**: per response, milestones the judge rules "not needed" leave the denominator; from the
   judge rows' per-milestone sources. Report per model with intervals beside MC as scored.
2. **Intermediate-only MC**: the answer targets leave the milestone set. Identify a target milestone by matching
   its value to the answer row's numeric targets under the check's unit factors (or by the milestone builder's
   line index for the final-answer line, if `evaluation/milestones.py` records it; check first). Report per model;
   state how many instances then have no milestone.
3. **Pairwise tests** on MC by matching alone, route-adjusted MC and intermediate-only MC: the 55 pairs with Holm
   within each family.
4. **Verbosity, cross-model** (R1 M1f; R2 W3): a response-level model of coverage on the number of numeric values
   the response displays (or visible tokens) with template fixed effects, pooled over models; report the slope.
   Keep the within-model Spearman beside it.
5. **The rule for the Claude against DeepSeek sentence**: keep "MC separates them" only if the pair separates after
   Holm under matching alone and under the route-adjusted variant. Otherwise the sentence becomes "as scored, MC
   separates Claude Sonnet 5 from DeepSeek V4.1 Flash; by matching alone it does not", or is dropped.
6. **Judge noise beside the gap**: in the coverage appendix, put the judge-swap shifts (−0.017 to +0.020) in the
   same paragraph as the top-tier MC gaps, so the reader sees them together.
7. **Kendall's tau wording**: "orders the models differently" becomes "agrees with the final-answer ordering at
   tau 0.709", with the interval.

Checks: the variants reproduce MC as scored when no milestone is removed (a unit check in the script). Effort: one
day. Feeds WP2.

### WP6. Template defects: repair or report (R1 F1, M2, m8, m11; R2 W1, W2, Q2)

Where: `data/templates/branches/chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.py`
(`work_isothermal_virial`) and `heat_effects.py` (`adiabatic_flame_temperature`); `template_annotation_23092026/`
(layer 0 gate, layer 2 round 5); `full_run_28092026/freeze.py`, `run_traces.py`, `score.py`; `3_4_design.tex`,
`6_results.tex`, `appendices/certification.tex`, `appendices/error_analysis.tex`.

**Path B (D4), in order:**

1. **Repair.** `work_isothermal_virial`: the question states the system (closed, reversible, isothermal, or the
   steady-flow reading, one of them) and the virial form to use (Z = 1 + BP/RT or the volume form, one of them),
   and the gold trace follows that reading. `adiabatic_flame_temperature`: the question prints the heat-capacity
   data the gold uses (the coefficients or the values at the reference states) and names the reference; the
   tolerance then holds the response to the data given. The template authors make the edits; no other template
   changes (the round-4 lesson: a shared table widened once changed certified instances elsewhere).
2. **Integrity checks** at 500 seeds (the layer-0 gate: closure, determinism, format); the display-tie census for
   the two templates; `docs/appendix_listings.py --check` is unaffected unless a listing changes.
3. **Round 5 certification**: `layer2/build_tasks --only template_work_isothermal_virial,template_adiabatic_flame_temperature`,
   `make_kits --round 5`, kits to the three chemical experts (E3a), `score.py` to `RESULTS_round5.md`,
   `CERTIFICATION.md` updated; `docs/appendix_certification.py` and the certification appendix then say five rounds
   and the round-5 counts.
4. **Re-freeze the two templates** with the same private seed and selection rule (`freeze.py` restricted to the
   two templates): 30 new instances; `FREEZE.json`, `manifest.jsonl` and `pool/` updated for them; the replaced
   hashes recorded; the Kaggle backup refreshed (the owner's manual step).
5. **Re-run** the 30 instances for the eleven models, plus the two templates' three subset instances in every arm
   that includes them (flagship, flagship-reasoning, reasoning-medium, openbook2 if their equations are stated,
   tool, paraphrase if their pairs were kept, repeats if in the 300): a `main-repair` variant of `run_traces.py`
   over those items, rows merged into the main store with provenance. Cost: a few dollars in all. Approval per
   the rules.
6. **Re-score** the affected rows (the store's `CONFIG.json` ties scores to the inputs' hashes, so either re-score
   the whole store or replace the rows with recorded provenance, per `score.py`'s own rules), then the judge and
   router stages for those rows.
7. **Expert readings.** B4 becomes letter material; B1 readings on the two templates' items leave the counts
   (`docs/appendix_evaluation.py` recomputes); the error-analysis sample loses the wrong answers on the two
   templates: redraw replacements per affected model from the re-scored store (`expert_kits.py --build` for the
   top-up), send E3b, and until it returns report the reduced sample with its counts.
8. **Regenerate**: `paper_results.py --write`, `docs/appendix_statistics.py --check`, `docs/appendix_evaluation.py
   --check`, `paper_setup.py --check`. Remove every "with and without the two chemical templates" clause (level
   gap, Thermodynamics domain, §5.4) and the "Templates Whose Wording Does Not Pin the Answer" paragraph's first two
   items; the six near-miss templates' finding stays.

**Path A (fallback):** keep the two templates; WP4.3's cleaned set beside the full set; introduce the two
templates in the main text at their first mention (one sentence saying what the three chemical experts found)
rather than in the appendix only; the Limitations sentence of section 7.

**Symbolic answers (D5, E2):**

1. Build the grading kit from `RESIDUAL_INCORRECT.md`'s symbolic column (the response's final expression and the
   gold's, template by template), as a reading sheet like the B kits; send E2.
2. Implement an equivalence check for `ber_estimation_mary` in `answer.py`: parse the final expression (LaTeX or
   plain) into SymPy, compare by simplification and by numeric sampling over the free symbols; validate against
   the expert grades (precision and recall on the graded items); apply only where validated.
3. Re-score the answer stage for the templates the check covers; regenerate. State the rule in the scoring
   appendix ("symbolic answers are checked for equivalence on N templates and by the numbers they state on the
   others").
4. If the parser cannot be validated by day 3: the expert-adjudicated FAC for the nine templates as a sensitivity
   row, and the by-numbers rule stated in the abstract.

**Read-off templates (4):** keep in the set (the sensitivity moves FAC by at most 0.004), name them in the
appendix, flag them in the release metadata.

**Certification's blind spot (R1 M2 last bullet; R2 W2):** one sentence in the certification appendix that the
planted defects test numeric correctness and not under-specification, and that the B4-type reading is what caught
the ambiguity; under Path B, a sentence that the repaired templates were re-certified. Report how many reasoning
paths the experts inspected per template (the five instances in the kit) as a fact.

Effort: Path B two days spread over days 1 to 4, gated by the experts; symbolic one day. Depends on D4, D5, E2, E3.

### WP7. Difficulty labels and the showcased templates (R1 M6, m10; R2 W1, Q6)

Where: `3_4_design.tex`, `appendices/taxonomy_content.tex`, `appendices/template_examples.tex`,
`docs/appendix_listings.py`, `docs/appendix_statistics.py`, `6_results.tex`.

1. **Rating kit** (E1): per branch, the 30 templates with one instance each and the three-dimension rubric; the
   current labels hidden. Score: Fleiss' kappa per branch and overall; majority label per template.
2. **Protocol sentence** in §3.1: who assigned the levels originally, and that three experts per branch rated them
   independently with the agreement; per D7, adopt the majority if it differs widely (then the level field in the
   template metadata and the manifest changes, every level table regenerates, and `FREEZE.json`'s by-level counts
   are recomputed).
3. **Depth as validation**: report the correlation of the level with the gold trace's milestone count and with
   the top tier's accuracy (offline, from the manifest and the scores); WP4.11 supplies the controlled depth
   result.
4. **Swap the listings**: `docs/appendix_listings.py` takes a template whose Advanced label meets §3.1's criteria
   (several coupled quantities, an iterative or implicit solve, six or more milestones; candidates from the
   manifest's milestone counts and levels) and an Easy template whose question does not state the formula;
   regenerate the listings; keep one listing per level from three branches.
5. **Method hints**: count the templates whose question states the governing formula or names the method (an
   audit by reading the 150 question skeletons, or by the experts in E1's kit with one extra tick box); report
   the count in the statistics appendix and qualify the Easy-slips sentence accordingly (R1 m10).
6. **Figure 3**: either all eleven models in a compact appendix figure or the four with the non-significant gaps
   marked (WP12).

Effort: half a day of scripting, the rest waits on E1. Depends on E1, D7.

### WP8. Evaluator validation: what the numbers mean (R1 M5, m4, m5, m9, m12; R2 W7, W8, W14, Q5)

Where: `5_evaluation.tex`, `appendices/validation.tex`, `appendices/scoring.tex`, `6_results.tex` (Table 1).

1. State how the 15 expert-study templates were chosen (E4) and that only 3 have a scalar answer, as now.
2. Say which settings were chosen on the expert labels and which were cross-fitted (ε on half the responses),
   in one sentence beside the headline agreement, so "validated against" is read correctly.
3. Keep the per-verdict-category readings on the evaluated models (they exist) and add a sentence that per-model
   scorer accuracy is not measured (Limitations).
4. **Judged step flags**: move the column to the appendix table or keep it with the caption saying its precision
   on correct answers (0.467 on 15 flags); either way the main text stops citing "up to 22.0%" without the
   precision beside it.
5. **Arithmetic flags and carried precision** (R2 W8): classify the 171 confirmed slips into cases where the
   result is consistent with unrounded upstream values present in the response and true slips (an offline check
   with the arithmetic module on the stored responses: recompute each flagged result from the unrounded values
   stated earlier); report both counts in the coverage appendix; if carried precision dominates for Claude
   Sonnet 5, say so in §5.3.1 instead of "slips that the answer passes".
6. **The judged step check is not deterministic**: say in §4 that it runs at the provider's default settings and
   that its flags are judged, not verified; optional, a stability sample (the router on 220 responses twice,
   about $3) with the flag agreement in the appendix.
7. Footnotes: the judge's 17 unreturned calls (51 of 60, 111 of 120); the two F1 values are two measurements
   (recorded matches against live matching), said in the caption in plain words.
8. PRM thresholds (R1 m9): one sentence that the 0.5 threshold is the models' default and that the pilot found the
   threshold is not the cause of the false flags (RESULTS_X1), or report the AUROC alone.
9. Error analysis by level with "no error" readings removed (R2 Q14): add the column to `tab:errors` and qualify
   the abstract's sentence on error types by difficulty (the Easy count is 26 wrong answers, 78 readings).

Effort: one day. Depends on E4.

### WP9. Taxonomy, sources and coverage (R1 M7, Q20; R2 W9, Q12)

Where: `3_4_design.tex`, `appendices/taxonomy_content.tex`, `custom.bib`, `docs/appendix_statistics.py`.

1. Replace `aiche2025constitution`, `asme2025vision` and `ieee2025standards` with curricular documents: the NCEES FE
   exam specifications for Chemical, Civil, Electrical and Computer, Industrial and Systems, and Mechanical (one
   entry each, current edition), or the ABET program criteria per discipline; keep ABET, ASCE BOK3 and IISE BOK.
   Check every entry by hand per the bibliography rules.
2. State the selection rule (D10) in §3.1 and the appendix in the same words: how the authors' lists were drawn
   from the standards and textbooks, what the LLM panel did (cross-check or vote), and how 5 to 7 prompted domains
   became 3 per branch.
3. Add a table of templates per area and level from `manifest.jsonl` (it carries `area`): 42 rows or grouped by
   domain, with the pedagogical-significance score if the authors can supply it (the outputs are held by the
   experts, PAPER_PLAN section 5); `docs/appendix_statistics.py --check` covers it.
4. Tone down "chosen to mirror formal engineering curricula" to "follows"; add one scope sentence naming the core
   areas the benchmark does not cover (circuit analysis; thermodynamics and heat transfer in mechanical;
   transportation and environmental in civil) and that each branch holds three domains by design.
5. Caption of the statistics table: "30 templates per branch" instead of "perfectly balanced distribution".
6. "ART-DEIT Database" needs a citation or a plain name; say which parameter values come from IEEE 100 (a
   dictionary of terms) or drop it.

Effort: half a day. Depends on D10 and the authors' facts.

### WP10. Related work and the novelty case (R1 M10; R2 W15)

Where: `2_relatedwork.tex`, `custom.bib`, `appendices/validation.tex`; optionally a new
`appendices/related_work.tex` (the plan's Appendix M, from `RELATED_WORK_v3.md`'s appendix part).

1. Add one sentence each, with citations, on dynamic and functional benchmarks (Functional Benchmarks, DyVal,
   MathCAMPS, GSM-Plus) in the symbolic paragraph, and on reference-free step evaluators (ROSCOE, ReCEval,
   ReasonEval) and PRMBench in the step-level paragraph; name SciBench, TheoremQA, OlympiadBench, PHYBench and
   UGPhysics among science benchmarks with step- or expression-level scoring. Add the missing bib entries (DyVal,
   ReCEval, TheoremQA, MathCAMPS) in their published versions.
2. State the delta from ThermoQA and FinChain in one sentence each: what they check (steps against computed
   values on a fixed set; step alignment with LLM judging) and what EngTrace does differently (deterministic
   recomputation of displayed arithmetic, a judge only on residual milestones, the certification protocol,
   instance generation from a private seed).
3. In the validation appendix, name the step-matching design as the FinChain-style step-alignment baseline if the
   pilot's design record confirms the correspondence (check `evaluator_pilot_17092026/` before writing it); then
   §2's "FinChain scores by step alignment" sentence can point to the comparison.
4. If the extended related work (the fifteen closest benchmarks) is added as an appendix, the main text ends with
   its pointer sentence.

Effort: half a day.

### WP11. Clarity: a worked example, the dense sentences, the main-text mentions (R1 M11; R2 W16)

Where: `5_evaluation.tex` (a box or figure), `6_results.tex`, `appendices/scoring.tex`.

1. **Worked example**: one real response from the store (a lower-tier wrong answer with a judge ruling and an
   arithmetic flag reads best): the question in brief, the gold milestones with values, which matched
   deterministically, which the judge ruled reached, not needed or missing, the arithmetic flag, then v(r) and
   m(r). As a `tcolorbox` in §4 if the page budget allows, else in the scoring appendix with a pointer. Built by a
   small script from the E3, judge and step rows so that every value is from the store.
2. **Dense-sentence pass** on §4 "Preliminaries" and "Milestone Coverage", §5.3.1 "Derivations Behind Wrong
   Answers" and the last paragraph of §5.4: one or two statistics per sentence, the definition of targets in two
   sentences.
3. Under Path A, the two chemical templates are introduced at their first main-text mention (§5.3.2) in one
   sentence; under Path B the mentions disappear.
4. "Verdict" is used for classification answer words and for scorer outcomes: use "label" for the former.
5. Conclusion (WP1).

Effort: one day.

### WP12. Table 1 and the figures (R1 m7, m13, §5 Tables and Figures; R2 W7, §5)

Where: `paper_results.py` (Table 1 block), the figure scripts under `full_run_28092026/` and `figures_oct_12/`.

1. **Table 1**: tier letters from `letters()` instead of bold "best"; the no-readable-answer column (WP4.4); the
   MC-by-matching column (WP4.2); "Calculations parsed per response"; the judged step flags moved or caveated
   (WP8.4); the reasoning-on rows (WP3.8). Caption rewritten accordingly.
2. **Figures**, each drawn, rendered and shown for approval before any tex changes: the error-category figure with
   distinct hues (and the empty "Unit" category merged into "Other"); the branch figure as grouped bars or a dot
   plot instead of fifths of a stacked bar; the radar replaced by a heatmap or dropped (its table carries the
   information); the level figure with all eleven models in the appendix or with non-significant gaps marked;
   the model label note per D9.
3. Greyscale check on every figure.

Effort: one day, after approval of the drawings.

### WP13. Small corrections (R1 §5, m1, m2, m3, m8, m12; R2 W13, W16, §5)

All in one pass, then `--check` on every script:

- "Three experiments" in §5.3.4 becomes four (or the paragraph names the four).
- The branch-spread sentence prints three decimals or says "about 0.05" (`paper_results.py` line 683, `f2`).
- The unit-of-analysis justification: "because a template's instances are generated by one procedure and
  scored against one gold derivation", not "share one derivation".
- 504 of 510 hand checks: "504 of the 510 hand checks compared a number"; split the 45 mismatches into real and
  planted items if the record allows (`layer2/RESULTS.md`).
- Round 3 approvals despite two mismatched hand checks, and the round-4 "extension": one sentence each from
  `fixes_round4.md` and `RESULTS_round3.md` (what the extension was; why the two mismatches did not block).
- Kappa for civil and industrial: print "undefined (all approve)" in the agreement table and say why.
- The judge's 17 unreturned calls; the two F1 values (WP8.7).
- `Mistral Large 3` gets a citation; `piepho2004letter` is cited or removed from `custom.bib`.
- "physically-validated" becomes "physically validated"; the quotation marks inside `\ttfamily` boxes use
  matching straight or curly quotes consistently; captions in sentence case where the rules say so.
- The listing precision policy: one sentence on which templates prescribe rounding and why (WP6's prescribed
  digits).
- Table 13's "Without two chemical" p-values: say whether Holm applies within that column.
- Query dates for every model (WP1.4).
- The abstract's "deterministic wherever a check exists" rewritten per WP2.2.

Effort: two to three hours.

### WP14. Final checks, the second mock review, submission

1. `paper_results.py --write` then `--check`; `paper_setup.py --check`; the four `docs/appendix_*.py --check`;
   `docs/check_plan_claims.py` if applicable.
2. The vocabulary grep of PAPER_PLAN section 8 (E3, E5, router, pilot, roster, run, arm, D-xxx, file names).
3. Sync to Overleaf per `docs/OVERLEAF_WORKFLOW.md`; compile; check every `\autoref` resolves; eight pages of
   content; the Limitations and Ethics sections present and unnumbered; no author information.
4. Proofread once end to end for the ARR guide's rules (American English, active voice, abbreviations defined).
5. **Second mock review**: the same prompt to the same two AI models on the revised bundle, this time with
   `main.tex`, the compiled PDF and every figure included, saved under `ACL Mock Reviews/second cycle/`. Read for
   regressions only; no new scope.
6. The Responsible NLP checklist (annotator details, licence, compute).
7. `docs/REVISION_LETTER_NOTES.md`: add the round-4 history, the B4 reading, the repair and re-run (Path B), the
   reasoning-on full-set run, and which real-reviewer points each answers.

---

## 6. Schedule

Day 1 is Tuesday 7 October. The deadline is Sunday 12 October, anywhere on Earth.

| Day | Owner (decisions, approvals, experts) | Writing | Scripts and runs |
|---|---|---|---|
| **1, Tue 7** | Decisions D1 to D10 in the morning. Expert asks E1, E2, E4 sent; E3 asked for if Path B. Approve the dry-run estimates for WP3 (and WP6.5 under Path B). | WP1 back matter drafted. WP13 small corrections. | WP3.1 to 3.5: roster check, variant, Gemma calibration, dry run, launch. WP4.1, 4.4, 4.5, 4.6 (single-path table, unreadable column, per-provider, model names). Path B: WP6.1 to 6.3 (repairs, gate, kits). WP6 symbolic kit built and sent. |
| **2, Wed 8** | Path B: chemical experts return the re-certification. | WP2 framing drafted (holding the MC sentence for WP5). WP9, WP10 text. | WP5 MC variants and tests. WP4.2, 4.3, 4.7 to 4.11. WP3.6 scoring stages as inference completes. Path B: WP6.4 to 6.6 (re-freeze, re-run, re-score). Symbolic parser attempt. |
| **3, Thu 9** | E1 ratings and E2 grades return. D7 decided. Figures approved. | WP2 finalized with WP5's outcome. WP8 text. WP11 worked example and dense-sentence pass. | WP3.7 to 3.8 (matched family, Table 1 rows, text). WP7 kappa and listings. Symbolic re-score or sensitivity. Path B: WP6.7 top-up kit sent; regenerate. WP12 drawings rendered and shown. |
| **4, Fri 10** | Path B: readers return the top-up. | Full read-through; Limitations final. | Path B: regenerate with the top-up. Every `--write` and `--check`. WP12 figures placed. Overleaf sync and compile. |
| **5 to 6, Sat 11 to Sun 12** | Final approval. | Proofread; second mock review read for regressions; checklist. | Buffer for anything late. Submit. |

Critical path: the owner's approvals and the expert asks on the morning of day 1; inference complete by the
morning of day 2; the experts' returns by the evening of day 3. If E1 or E2 is late, the paper ships with the
agreement sentence and the symbolic sensitivity left out, and the Limitations say so.

---

## 7. Out of scope this cycle, and the Limitations sentence for each

| Not done this cycle | Sentence for the Limitations (one each, measured tone) |
|---|---|
| A harder template tier (coupled systems, design problems, multi-domain chains) | The templates are examination-style problems of undergraduate curricula; the strongest models solve nearly all of them on the final answer, and the benchmark separates those models on the derivation rather than the answer. |
| A fresh expert-labelled validation sample on the evaluated models; per-model scorer accuracy | The evaluator's agreement figures come from an expert study on five models outside the evaluated eleven; on the evaluated models, experts read samples of the final-answer verdicts, the judge's rulings and the arithmetic flags, reported per category, and scorer accuracy per model is not measured. |
| A required-tool condition for every model | The tool experiment offers the interpreter to two closed models; it bounds what an offered tool adds, not what required use would, and the open-weights model could not be served with a tool. |
| One pinned provider per open-weights model | Open-weights models were served by several providers at 8-bit precision or higher; the per-provider differences, matched on template, are reported in the appendix. |
| A full second-judge pass | The judge decides 11% to 25% of milestones; a second judge from a third family moves no model's coverage by more than 0.02 on a 220-response sample. |
| Contamination after release (D8) | The evaluation set becomes public with its seed at publication; the generator draws fresh sets from a new seed, and a held-out seed is kept. |
| MC's nature | Milestone Coverage credits the intermediate values a response states and the judge rules reached; a valid route that bypasses a milestone is not credited, and matching reads the visible response, not the reasoning text. |
| Unequal effort for the models with no reasoning setting (after WP3) | Models that offer no reasoning setting are compared at their providers' defaults. |
| The two chemical templates (Path A only) | Two Advanced chemical templates admit more than one textbook reading or depend on the data source; the headline results are reported with and without them. |
| The LLM screen's contribution | The LLM screen guides revision and does not certify. |

---

## 8. Traceability: every review item and where it is answered

Review 1 (`ARR_Review_EngTrace_1.md`):

| Item | In a few words | Answered by |
|---|---|---|
| F1 | Saturated at the top; residual errors are artifacts | WP2 (ceiling stated), WP4.3 (cleaned set), WP6 (repair or report), section 7 |
| F2 | "Verifiable process supervision" unsupported; premise undercut | WP2 (title, abstract, intro, contributions), WP8.4 to 8.6 |
| C1 | No Limitations, Ethics, Conclusion, expert details, release | WP1, E4 |
| M1 | MC confounded; separation rests on one pair | WP5 (variants, verbosity model, rule), WP4.2, WP4.7 |
| M2 | Defective templates and scorer artifacts in the set | WP6 (Path B or A; symbolic; prescribed digits; read-off), WP4.3, WP4.9 |
| M3 | Sign and precision clauses; crude partial credit | WP4.8, WP4.10 |
| M4 | Unequal conditions; providers; frontier exclusion | WP3, WP4.4, WP4.5, WP2.5 (open-vs-closed dropped); frontier models stay anchors (stated) |
| M5 | Validated on design data; unrepresentative templates | WP8.1 to 8.3, section 7 |
| M6 | Difficulty labels unvalidated; depth confounded | WP7, WP4.11 |
| M7 | Taxonomy by LLM votes; sources not curricula; coverage narrow | WP9 |
| M8 | Contamination control narrow; no plan after release | WP1.4 (D8), section 7; the paraphrase margin's basis stated in the appendix |
| M9a to g | Missing or selective statements | WP4.1 (a), WP3 and the appendix (b: the one model named), WP4.6 (c), WP4.1 (d), WP2.5 (e, f), WP4.4 (g) |
| M10 | Novelty incremental; related work gaps; no ChainEval baseline | WP10 |
| M11 | Dense; no worked example; no conclusion | WP11, WP1 |
| m1 to m13 | Minor | m1, m2, m3, m8, m12: WP13; m4, m5, m9: WP8; m6: WP5.6 (judge noise stated beside the gap); m7, m13: WP12; m10: WP7.5, WP8.9; m11: section 7 |
| §4 Q1 to Q20 | Questions | Q1, Q10, Q11, Q19, Q20: E4 and WP13; Q2, Q15: WP7, E1; Q3, Q4: WP4.8; Q5, Q6: WP5, WP4.7; Q7, Q8: WP4.6, WP4.1; Q9: WP6; Q12: section 7; Q13: WP4.5; Q14: WP3; Q16, Q17: WP1.4; Q18: section 7 |
| §5 | Suggested analysis, wording, tables, figures, citations, typography | Carried-precision analysis WP8.5; human solve time not done (section 7); wording WP2 and WP13; tables WP12; figures WP12; citations WP9, WP13; typography WP13 |

Review 2 (`ARR_Review_EngTrace_2.md`):

| Item | In a few words | Answered by |
|---|---|---|
| W1 | Near ceiling; Advanced not advanced | WP2, WP7, WP4.3, section 7 |
| W2 | Certification missed defects; defective templates kept | WP6, the certification sentence in WP6 |
| W3 | MC measures route conformity and display | WP5, WP4.7 |
| W4 | Overclaim in the title and abstract | WP2 (D1) |
| W5 | Premise undercut by own data | WP2.3 |
| W6 | Configuration confound; open-vs-closed misleading | WP3, WP2.5 |
| W7 | Step-flag columns over-presented | WP8.4, WP12.1 |
| W8 | Slips may be carried precision | WP8.5 |
| W9 | LLM-decided composition; citations not curricula; no per-area table | WP9 |
| W10 | Scoring policy moves FAC more than the spread | WP4.9, WP6 (prescribed digits), WP4.3 |
| W11 | Required ARR information missing | WP1, E4 |
| W12 | Serving heterogeneity; output ceiling | WP4.5, WP4.4 |
| W13 | Unsupported or inconsistent numbers | WP4.1, WP13 |
| W14 | Small validation samples | WP8.3, WP8.9, section 7 |
| W15 | Novelty; related-work gaps | WP10 |
| W16 | Clarity; undefined macro; figures missing; naming | WP11, WP13, WP14.5 (bundle with main.tex and figures), D9 |
| Q1 to Q14 | Questions | Q1, Q2, Q11: E4 and WP6; Q3, Q4: WP4.7, WP5; Q5: WP8.5; Q6: WP7; Q7: WP3; Q8: WP4.1, WP13; Q9: WP1.4, WP4.5; Q10: WP1.4; Q12: WP9.3; Q13: WP13 (one sentence on the 1% hand-check tolerance against the 0.2% answer tolerance); Q14: WP8.9 |
| §5 | Comments | Macro and figures: WP14.5; abstract: WP2.2; curricular citations and per-area table: WP9; Eq. 1 clarifications: WP8.2, WP11.2; Table 1: WP12.1; unit of analysis: WP13; branch spread: WP13; Fig. 3: WP7.6; Fig. 4: WP12.2; radar: WP12.2; "perfectly balanced": WP9.5; pre- or post-screen versions in round one: WP13; Opus 4.5 against 4.7: one sentence in the validation appendix; script release: WP1.4; naming: D9; electrical kappa 0.365: one sentence in the certification appendix |
