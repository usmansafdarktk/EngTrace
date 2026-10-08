# WS-D1: Main text and back matter (Phase 2, after INTEGRATED)

## Mission

Rewrite the main text so that every claim matches its evidence: the title and the term, the abstract, the
introduction and contributions, the design section's taxonomy and difficulty paragraphs, the evaluation
framework, the results prose around the regenerated tables, and the back matter (Conclusion, Limitations,
Ethics). Every number comes from `reports/NUMBERS_SHEET.md`; a sentence that waits on a number carries a FILL
marker.

Closes: review 1 F1, F2, C1 (Limitations, Ethics, Conclusion), M4 (text), M7 (§3.1 text), M8 (Limitations), M9e,
M9f, M10 (§1, §2), M11 (dense sentences, main-text mentions), m6, m8, m10, m11, §5 wording; review 2 W1, W4, W5,
W6 (text), W9 (§3.1 text), W11, W15, W16, Q1 (the facts), Q10 (the Limitations sentence), §5 abstract and
wording comments.

## Context (read `PLAN_CONTEXT.md` sections 1 to 3, then `reports/NUMBERS_SHEET.md`)

The two reviews' framing objections, in one paragraph each:
- **Saturation.** The strongest models are at ceiling on the final answer and most of their remaining incorrect
  verdicts were artifacts; after Phase 1 the two chemical templates are repaired and the symbolic answers are
  checked for equivalence where validated, so the residual is smaller, but the ceiling remains. State it as a
  finding, not as "headroom".
- **The term.** "Verifiable process supervision" promises step-level verification and a training signal; what
  is verified deterministically is the final value, the presence of the gold trace's intermediate values and
  displayed arithmetic; conceptual step errors are not detected deterministically (0 of 60 planted) and judges
  catch about a third. The premise "a correct answer can follow flawed reasoning" is rare on this benchmark
  (AUROC 0.974; 175 of 178 flawed steps behind correct answers are slips).
- **Configuration.** With the matched-settings evaluation done, say both configurations and state the tiers
  under each; never "four of the top five are open-weights".
- **Milestone Coverage.** Say what it credits and does not; keep the Claude-versus-DeepSeek separation sentence
  only if the numbers sheet says it holds under matching alone and route-adjusted MC; tau 0.709 is agreement, not
  a different ranking.

## Rules in this session

- Repository `C:\Users\ayesha.gull01\EngTrace`, `master` only.
- WS-D2 works at the same time on the appendices and the bibliography. You own the files listed below and
  nothing else. `custom.bib`: add entries only under a `% --- added by WS-D1` marker at the end; never edit
  WS-D2's block. `paper_results.py`: you alone edit it in Phase 2 (the phrases section and the captions of
  main-text and appendix blocks; WS-D2 sends you appendix caption changes through its report or the owner);
  after an edit run `python full_run_28092026/paper_results.py --write --text-only --headline <D2> [--repaired]`
  and then `--check`.
- Writing rules: `overleaf_source_04102026/WRITING_RULES.md` and `docs/PAPER_PLAN_OCT2026.md` section 8 (plain
  language; no internal vocabulary; the "do not state" list in section 3.2; final state only, no history, no
  cost figures; `\ourdataset` never "EngTrace"; citations tied with `~`; American English; active voice;
  "measure" not "metric"; "significant" only with a test).
- The main text is eight pages of content; cut evidence to the appendix before cutting a limit.
- A sentence that waits on a number: write it with the number's name in brackets and a `%% FILL: <what>` line
  above it. The orchestrator resolves every FILL at the end.
- No commits or pushes unless the owner asks. When asked: one short line, no body, push at once.
- Finish by writing `docs/mock_review_workstreams/reports/WS-D1_report.md` from the template at the end.

## Before you start

- `reports/NUMBERS_SHEET.md` (WS-G) and the decisions D1, D2, D8, D9, D10 and the E4 facts in
  `00_ORCHESTRATION.md` (sections 3 and 5). If an E4 fact is missing, write the sentence with a FILL marker.
- Read the current files end to end: `current_overleaf_project/main.tex` (title, back matter),
  `overleaf_source_04102026/0_abstract.tex`, `1_intro.tex`, `2_relatedwork.tex`, `3_4_design.tex`,
  `5_evaluation.tex`, `6_experiments.tex`, `6_results.tex` (with its regenerated blocks and phrases), and
  `docs/related_work_oct2026/RELATED_WORK_v3.md` (the revised related work and its appendix part).
- The registry of labels (`00_ORCHESTRATION.md` section 7), so every `\autoref` points at a block that exists.

## Files you own

`current_overleaf_project/main.tex` (title, the Conclusion, Limitations and Ethics sections, the
`\input` list if a new appendix file is added by WS-D2: coordinate through the owner);
`overleaf_source_04102026/0_abstract.tex`, `1_intro.tex`, `2_relatedwork.tex`, `3_4_design.tex`,
`5_evaluation.tex`, `6_experiments.tex`, `6_results.tex` (prose outside the generated blocks; the blocks through
`paper_results.py`); `full_run_28092026/paper_results.py` (phrases and captions); `full_run_28092026/paper_setup.py`
(the sentences it checks in `6_experiments.tex` and `appendices/models.tex` prose it owns: coordinate the
appendix prose with WS-D2); `custom.bib` under your marker.

## Steps

### D1-1. Title and the term

`\title{}` per D1. Replace "process supervision" throughout your files with "process evaluation" (grep the whole
`overleaf_source_04102026/` tree and report any occurrence in WS-D2's files to the owner rather than editing
them). "Verifiable" stays where what follows is verified: final values, intermediate values, displayed arithmetic.

### D1-2. Abstract (at most 200 words; no citations, no section references)

In order: the gap (engineering derivations have intermediate quantities that can be checked; most engineering
benchmarks score only final answers; those that assess steps rely on LLM judges or cover one domain); the
benchmark (150 expert-certified templates, five branches, 15 domains, 42 areas; one computation emits the
answer and the gold trace); what is checked deterministically (the final value for six answer kinds, with the
symbolic rule stated as the numbers sheet gives it; coverage of the gold trace's intermediate quantities;
displayed arithmetic) and what is judged (residual milestones, by a judge from outside the evaluated families);
what is not detected deterministically (a misstated rule behind a correct answer); the evaluation (the number of
models and configurations, 2,250 instances from a private seed); the finding stated as a ceiling ("the five
strongest are not separable at this size on the final answer"); what the process measures add (where wrong
derivations break; arithmetic slips behind right answers; a comparison when answers tie); the
matched-settings sentence (FILL from the sheet); the error-type sentence qualified by level and sample size.

### D1-3. Introduction

- Motivation reframed: on examination-style derivations the answer nearly settles soundness; what step-level
  checks add is where a wrong derivation breaks, the slips behind right answers, and a comparison when answers
  tie. Keep the ProcessBench citation as the general observation, and state the finding beside it.
- The benchmark paragraph as now, with the certification sentence and the private seed.
- The evaluator paragraph: deterministic wherever a check exists, said concretely; the judge only on residual
  milestones; why from outside the evaluated families.
- Contributions, three bullets: (1) the certified benchmark; (2) the evaluator validated against step-level
  expert labels, with what verification can and cannot see; (3) the evaluation at matched settings: the
  final answer no longer separates the strongest models on these problems, and what milestone coverage, the
  step checks and the expert reading show instead.
- The scope sentence: a depth-oriented subset of engineering; examination-style by construction.
- Do not: "first", "no existing benchmark verifies this process", "cannot separate", "safety-critical" beyond
  one clause.

### D1-4. Related work (0.75 page)

Add one sentence with citations on dynamic and functional benchmarks (Functional Benchmarks, DyVal, MathCAMPS,
GSM-Plus) in the symbolic paragraph; one on reference-free step evaluators (ROSCOE, ReCEval, ReasonEval) and
PRMBench in the step-level paragraph; name SciBench, TheoremQA, OlympiadBench, PHYBench and UGPhysics among
science benchmarks with step- or expression-level scoring. State the delta from ThermoQA and FinChain in one
sentence each (what they check; what this benchmark does differently: recomputation of displayed arithmetic, a
judge only on residual milestones, the certification protocol, instances from a private seed). Point to the
extended related work appendix if WS-D2 adds it. Bib keys: use those in
`docs/related_work_oct2026/related_work_sources.bib` where they exist; add the missing entries (DyVal, ReCEval,
TheoremQA, MathCAMPS) in published form under your marker in `custom.bib`.

### D1-5. Design section

- §3.1: the taxonomy sentence "chosen to mirror formal engineering curricula" becomes "follows"; the standards
  citations as WS-D2 replaces them (use the new keys from its report; FILL until then); the selection rule per
  D10 in the same words the appendix uses; a pointer to `tab:area`; one scope sentence naming the core areas not
  covered (circuit analysis; thermodynamics and heat transfer in mechanical; transportation and environmental in
  civil).
- Difficulty (D7, owner 2026-10-08): the protocol sentence from
  `template_annotation_23092026/levels/difficulty_labelling_protocol.md`: the domain experts' original labels,
  unchanged since the templates were written. No agreement figure, no `tab:levels_agreement`, no statistic from
  the protocol file, and no re-rating is mentioned.
- §3.3: the planted defects test numeric correctness, not under-specification (one sentence); the independence
  sentence from E4 (who wrote the templates and who certified them); five rounds; "all 150 certified by a
  unanimous vote" stays if true after round 5.

### D1-6. Evaluation framework (§4)

- The first sentence: "we check deterministically what can be checked: the final value, the presence of the
  gold trace's intermediate values, and displayed arithmetic; a judge decides what these checks leave open".
- Preliminaries and targets in two sentences each (the dense-sentence pass).
- Milestone Coverage: say that matching reads the visible response and not the reasoning text an endpoint
  returns separately; that a valid route bypassing a milestone is not credited; the judge-decided share.
- Step-level diagnostics: the precision beside each flag rate; the judged step check runs at the provider's
  default settings and is judged, not verified; one sentence on carried precision if the sheet says it matters.
- The validation sentence: which settings were chosen on the expert labels and which were cross-fitted.
- A pointer to the worked example (`box:worked_example`, in the scoring appendix unless the page budget admits
  it here).

### D1-7. Experiments and results (§5)

- §5.1 and §5.2 (`6_experiments.tex`): the models and both configurations (the provider-default run and the
  matched-settings run; which endpoints offer no reasoning setting; the twelfth model if any); the response
  counts; the unit-of-analysis sentence corrected ("a template's instances are generated by one procedure and
  scored against one gold derivation"); `paper_setup.py --check` must pass after your edit.
- §5.3 (`6_results.tex` prose and the phrases in `paper_results.py`): under the chosen headline, read Table 1
  first; "two tiers" stated under both configurations; the open-versus-closed sentence removed; the
  Claude-versus-DeepSeek sentence per the rule; "orders the models differently" replaced by the tau sentence;
  the single-path sentence points to `tab:single_path`; the four repeat models named; "Stability of the Tiers"
  gains the readable-only reordering and the matched-settings sentence; "four experiments"; the branch spread
  at three decimals; the tool and equations conclusion as a bound; the "with and without the two chemical
  templates" clauses gone (`--repaired`); Thermodynamics sentences as regenerated.
- §5.4: the error reading as regenerated (the top-up; the matched-configuration readings if E3c ran, with the
  configuration said); the close rewritten as the ceiling: on the final answer the evaluation set is solved at
  the top by the models evaluated, what separates them is the derivation, and a response that misstates the
  rule it applies but computes the right answer passes every deterministic check, so the process measures flag
  where to read and the expert reading certifies. No "headroom".
- The dense-sentence pass on §5.3.1 "Derivations Behind Wrong Answers" and the last paragraph of §5.4: one or two
  statistics per sentence.

### D1-8. Back matter (`main.tex`)

- Conclusion (0.3 page): what the benchmark is, what the evaluation found at matched settings, what the process
  measures add, what comes next (method memorisation, multi-modal artefacts, the open-weights tier under tool
  use). No "first".
- Limitations (unnumbered, required): one scope sentence each, measured tone, from `PLAN_CONTEXT.md` section 7
  minus the ones Phase 1 resolved (the two chemical templates; unequal effort, now reduced to the models whose
  endpoints offer no setting), plus: examination-style problems and the ceiling; the evaluator's validation
  population and the absence of per-model scorer accuracy; the offered tool; several providers; the single
  judge with the swap result; contamination after release (D8); what Milestone Coverage credits; the paraphrase
  test's coverage; the subset experiments.
- Ethics (unnumbered): synthetic data from referenced constants; the experts (E4: who, recruitment,
  compensation, consent; FILL where missing); their labels and names withheld; models not to be deployed in
  safety-critical systems on this evidence; the licence; a pointer to the release appendix (`appendix:release`,
  WS-D2).
- Remove the commented acknowledgements block from the review version.

### D1-9. Consistency pass on your files

The vocabulary grep (PAPER_PLAN section 8: E3, E5, router, pilot, roster, run, arm, D-xxx, file names);
"verdict" only for scorer outcomes ("label" for classification words); one spelling of every model name; every
`\autoref` points at an existing label; abbreviations defined at first use in the abstract and the body;
`paper_results.py --check`, `paper_setup.py --check`, `docs/check_plan_claims.py` (if it still applies) pass.

## Acceptance

- Every file you own rewritten per the steps; eight pages of content when compiled with WS-D2's appendices;
  no FILL marker left that you can resolve from the numbers sheet; the remaining FILL markers listed in the
  report.
- `paper_results.py --check` and `paper_setup.py --check` pass; the vocabulary grep returns nothing in your
  files.

## Report template (`reports/WS-D1_report.md`)

```
# WS-D1 report
- Files rewritten; what each now says that it did not (one line per file).
- Phrases and captions changed in paper_results.py; regeneration run (yes/no), --check (pass/fail).
- FILL markers left (file, line, what is needed).
- Bib entries added under the WS-D1 marker.
- Sentences for WS-D2 (appendix prose that must match the main text: list).
- Occurrences of "process supervision" or other must-change terms found in files you do not own.
- Open items.
```
