# Supervisor comments of 7 October 2026 and the changes made in response

The four comments Zhuohan Xie left on the draft on 7 October 2026, verbatim, each followed by what the paper does
about it and where. The comments predate the mock-review workstreams (WS-A to WS-G, 7 to 9 October) and the
review pass of 9 October (D-197); the responses describe the paper as it stands on 9 October.

---

## 1. Contribution and Main Story

> **Zhuohan Xie, 7 October, 3:44 pm. [Zhuohan: Contribution and Main Story]**
> The motivation is useful, but please align the title, abstract and contributions with what the evaluator
> actually establishes. Milestone coverage measures the intermediate values a response states, and arithmetic
> checks inspect parsable displayed calculations. These do not certify the correctness of the complete reasoning
> process: the deterministic checks miss all 60 planted conceptual defects, as your own study reports. Define the
> scope of verifiable process supervision explicitly, including its dependence on LLM judging and expert review,
> and avoid implying that high coverage guarantees sound reasoning. The strongest story is what process
> diagnostics reveal when final-answer scores are nearly tied, together with their measured blind spots. Bring
> that story into the introduction and conclusion, and condense repeated tier descriptions and robustness
> narration. This does not require adding a training experiment just to justify the title.

**Response.**

- Title: "Verifiable Process Evaluation of Engineering Reasoning from Symbolic Templates" (decision D1); "process
  supervision" appears nowhere in the paper.
- Abstract (`overleaf_source_04102026/0_abstract.tex`): says what is deterministic (final answers of six kinds,
  coverage of the gold trace's intermediate quantities, displayed arithmetic), that residual milestones go to an
  LLM judge outside the evaluated families, and that "a misstated rule behind a correct answer escapes these
  checks"; the five strongest are "within 0.02 and near the ceiling".
- Introduction (`1_intro.tex`): the expert-study premise names its population (300 responses from five other
  LLMs); the sentence on what step-level checks add states the tied-models story (how much of the gold derivation
  a wrong answer keeps and which quantities it misses, the slips behind right answers, a comparison of models
  whose answers tie) and that the blind spot is measured; contribution 2 gives the validation figures with their
  basis (in-sample; milestone F1 0.926 by matching alone, 0.958 with the judge) and the blind spot.
- Section 4 (`5_evaluation.tex`): MC "measures how much of the gold route a derivation states, not whether it is
  free of error"; matching compares values only; the flags are "flags with measured precision, not scores"; the
  0 of 60 planted misstated rules are stated.
- Section 5 (`6_results.tex`): "MC credits the intermediate values a response states, so it ranks how much of the
  gold route a model shows, not how soundly it reasons"; the Claude Sonnet 5 against DeepSeek V4.1 Flash separation
  is read as route conformity, with the basis of that reading (the judge's not-needed rulings, which the experts
  did not read); the tier description appears once, with robustness in one paragraph and its detail in
  `appendices/results.tex` and `appendices/paraphrase.tex`.
- Conclusion (`current_overleaf_project/main.tex`): the ceiling on the final answer, what the process measures
  show, the blind spot, the reasoning gains and the error mix; no rankings.
- No training experiment was added.

## 2. Completion and Presentation

> **Zhuohan Xie, 7 October, 3:52 pm. [Zhuohan: Completion and Presentation]**
> Please complete the currently empty Conclusion and Limitations. The conclusion should state the contribution
> and strongest supported findings, rather than repeat model rankings. Limitations should cover the bounded
> curriculum/template scope, ambiguous questions, answer-checking simplifications, valid alternative derivations,
> incomplete detection of conceptual errors, and provider-default reasoning settings. Condense the main text
> around the key findings, with detailed robustness results in the appendix. Update the old EngChain name inside
> Figure 2 to EngTrace, make figure text readable, and position figures after their first textual citation,
> preferably on the same or next page. Shorten unnecessarily long headings and balance section lengths. Resolve the
> four duplicate BibTeX entries reported by the current build: kendall1938, spearman1904, wilcoxon1945 and
> mcnemar1947. Prioritize these revisions and analysis of saved outputs before proposing additional model or judge
> runs.

**Response.**

- Conclusion and Limitations are written (`current_overleaf_project/main.tex`). Limitations cover, one sentence
  each: the examination-style curriculum scope; under-specified wording (certification catches numeric defects,
  and the experts' near-miss reading covers six templates); the answer-rule simplifications (symbolic answers on
  four templates scored by stated numbers, half credit for a partly correct multipart answer, prescribed digits)
  with their measured effect (at most 0.030 per clause, Kendall's tau 0.927 or above); valid alternative routes
  and "high coverage does not certify sound reasoning"; conceptual errors (0 of 60 planted misstated rules, recall
  about a third of flawed steps inside correct answers, judged-flag precision 0.467); and the provider-default
  settings with the matched configuration beside them.
- The main text was rewritten around the findings (WS-D1, D-195, D-197); robustness detail sits in the appendix.
- Figure 2 (`figs/template-generation.pdf`) carries no benchmark name; the overview figure says EngTrace. The
  old name survives only in the Overleaf mirror's unused file `current_overleaf_project/figs/engchain-overview.pdf`.
  Figures are placed after their first citation in the sources; readability is checked on the compiled build.
- Headings: the experiments appendix heading names its four experiments by the authors' rule
  (`docs/PAPER_PLAN_OCT2026.md`); no other heading is long.
- BibTeX: `overleaf_source_04102026/custom.bib` holds one entry per key (86 entries, every one cited). The
  duplicates came from an earlier sync that concatenated the old and new sets into the mirror's bibliography.
- No further model or judge runs were proposed; the remaining work was analysis of saved outputs.

## 3. Final Benchmark Quality

> **Zhuohan Xie, 7 October, 3:46 pm. [Zhuohan: Final Benchmark Quality]**
> The expert certification work is substantive, but the final dataset status must be reconciled with the later
> audit. Section 5.4 and Appendix Q identify two Advanced chemical templates whose wording does not determine a
> unique answer, and these templates affect the reported domain and difficulty patterns. Please state their
> current status explicitly. Correct their assumptions/data sources and re-certify them, or exclude them from the
> primary valid set while retaining them as documented diagnostic cases. Use the saved predictions to report the
> effect on the main results where rescoring is sufficient; a footnote or an optional exclusion analysis does not
> by itself establish that the primary benchmark has unambiguous gold answers. Clarify the version and timing of
> certification relative to these later findings, and make the final valid template/instance counts consistent
> throughout.

**Response.** The first option was taken (decision D4, Path B; WS-A).

- The two templates (`work_isothermal_virial`, `adiabatic_flame_temperature`) were repaired: the question now
  states the system, the form of the equation, and the property data on which the answer depends. They were
  re-certified by the three chemical experts in rounds five and six (round five: one rejection, on physical
  plausibility; round six: unanimous), their 30 instances were re-drawn from the same seed and re-run for every
  model and every experiment, and every result was regenerated from the saved outputs.
- The paper states the final state: `3_4_design.tex` ("after six rounds, all 150 templates are certified by a
  unanimous vote"); `appendices/certification.tex` describes rounds five and six and ends "Every template is thus
  approved unanimously in the latest round that reviewed it, its current code regenerates, byte for byte, the
  instances its experts reviewed, and the evaluation set holds instances of these certified versions only." The
  counts are 150 templates and 2,250 instances throughout, with no "with and without" clauses left.
- The audit that found the wording defects, the repair and the re-run are recorded for the revision letter in
  `docs/REVISION_LETTER_NOTES.md` (section 7) and in `docs/mock_review_workstreams/reports/WS-A_report.md`; the
  paper itself describes the benchmark as it is.
- The near-miss reading (`appendices/error_analysis.tex`) covers the six other templates on which many wrong
  answers land close to the gold value: all three experts find one correct answer and the response shown wrong.

## 4. Coverage and Alternative Derivations

> **Zhuohan Xie, 7 October, 3:47 pm. [Zhuohan: Coverage and Alternative Derivations]**
> Please make the interpretation of MC explicit in the main text and model comparisons. It is coverage of the
> reference intermediate quantities, not verification that each step is valid. Appendix J reports that some
> missing milestones are unnecessary for a valid alternative route, while 15 of 100 readings of reached rulings
> were not confirmed. Show one concrete alternative correct derivation and one case with matched values but flawed
> reasoning. Explain which quantity correspondence, units and dependencies are checked and which are not.
> Distinguish deterministic matches from judge-assisted coverage and explain the observed judge contribution. A
> lower MC score can reflect a different route or less displayed detail, so do not equate the MC ranking with a
> ranking of reasoning soundness. Keep the measured limitations prominent when explaining what the process
> diagnostics add beyond final-answer scores.

**Response.**

- Interpretation: Section 4 defines MC as "how much of the gold route a derivation states, not whether it is free
  of error"; Section 5 adds "it ranks how much of the gold route a model shows, not how soundly it reasons";
  Limitations: "high coverage does not certify sound reasoning".
- What is checked: Section 4 and `appendices/scoring.tex` now say that matching compares values only, not the
  quantity's name, the unit written beside it, the step it appears in, or the order and dependencies among steps,
  so a coincidental value counts, which the chance floor of `tab:coverage` measures; the judge sees each
  milestone's name and value and rules on the response's own route.
- Two concrete cases, `box:coverage_cases` (`appendices/coverage_cases.tex`, generated by
  `full_run_28092026/coverage_cases.py` from the saved outputs): A, a correct vorticity derivation by Kimi K3 that
  differentiates symbolically and never states the two partial derivatives the gold trace evaluates, so the judge
  rules them not needed and m(r) = 1/3 for a valid derivation; B, a planted defect of the validation study in which
  "R = A / P" became "R = P / A" with no digit changed, which matching, the arithmetic check and the answer check
  all pass, and which two of the four judges label a conceptual error and the paper's judge does not. Pointed to
  from Section 4 and the scoring appendix.
- Deterministic against judge-assisted coverage: Table 1 states that the judge decides 11% to 25% of a model's
  milestones and points to `tab:coverage_variants`, which prints MC as scored and by matching alone;
  `appendices/results.tex` now states the observed contribution: the judge adds 0.019 to 0.038 to a model's
  coverage over matching alone, and the pairs separated are nearly the same (28 against 26 of the 55); the
  Claude Sonnet 5 against DeepSeek V4.1 Flash separation is shown not to rest on the judge (matching alone
  separates them more strongly) and to disappear route-adjusted.
- The readings: `appendices/validation.tex` keeps both figures (a third of the confirmed missing rulings are on a
  route that does not need the milestone; 15 of the 98 readings of reached rulings do not confirm them), and
  Section 4 now cites the out-of-sample readings on the evaluated models.
- Measured limitations beside the diagnostics: every flag rate carries its precision (0.905 for arithmetic flags on
  189 decided flags; 0.467 for judged step flags inside correct answers) and the recall of both checks together
  (0.36 inside correct answers), in Section 5 and in Limitations.
