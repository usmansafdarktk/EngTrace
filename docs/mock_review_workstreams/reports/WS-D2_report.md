SIGNAL: WS-D2 DONE 2026-10-09 (all six checks pass on the final tree, after WS-D1's phrase updates)

# WS-D2 report

2026-10-09. Appendices, bibliography and listings written against `reports/NUMBERS_SHEET.md`. Nothing committed.

## Files rewritten

- `appendices/certification.tex`: six rounds as a list (round 3's two mismatches are the gold value entered in pascals
  against an answer in kilopascals; round 4 widens two templates to 15 distinct instances each; round 5 re-certifies
  the two chemical templates whose questions state system, equation form and property data, one rejection on hot
  organics; round 6 holds organics below decomposition, unanimous); 504 of the 510 hand checks compare a number;
  planted defects test wrong numbers, not under-specification; κ "--" (undefined) for civil and industrial, with the
  caption note; electrical κ 0.365 explained (5 of 90 decisions are rejections); "decisions" for the experts'
  approve/reject; sentence-case captions.
- `appendices/scoring.tex`: the 1% cap on the last-digit term; new paragraphs Symbolic Answers (five of nine enabled;
  1.000 / 0.990 on the 264 graded verdicts, which also settled the rule; 1.000 / 0.840 on 30 readings the rule did not
  see), Prescribed Digits (17 templates; kept; relaxed: 154 verdicts on 15 templates, τ 0.964), Sensitivity of the Rule
  (absolute-value clause 325 on 20, τ 0.927; uncapped last digit 68 on 36; proportional credit 109 on 14; accepted
  answers by relative error); 57 instances without milestones; `tab:scoring_settings` rows current; every check reads
  the response, not the separate reasoning text; judged step flags are rulings, with pointers; Worked Example
  subsection introducing `box:worked_example`; "label word" for classification answers.
- `appendices/validation.tex`: study-template rule (one per branch and level, least-represented answer kind, seeded
  ties, the four read-off templates left out); only ε cross-fitted, so the agreement figures are in-sample; Opus 4.5
  in the panel against Opus 4.7 in the study; PRM threshold sentence (−0.006; no threshold reaches 0.80); table rows
  current; ‡ note on the judge's 17 of 240 calls without a reply (51 and 111); † note: the two milestone F1 values are
  two measurements; judge swap 201 responses, −0.018 to +0.019; new subsection `appendix:readings` (B1 144 items, B3 97,
  partial verdicts 18 of 36 called correct and 9 incorrect, per-model scorer accuracy not measured);
  `tab:judged_steps` introduced with precision 0.640 / 0.467 (7 of 15).
- `appendices/models.tex`: decoding rows from `paper_setup.py`; Matched Settings paragraph (three models at effort
  medium, every response with reasoning tokens, Gemma 3 empty; Qwen3-235B-2507 has no setting; its Thinking sibling not
  evaluated; dates in `appendix:release`); Serving Endpoints paragraph (21 endpoints matched on ≥ 20 templates, −0.025
  to +0.021); protocol: 148 templates with milestones, "nine of the eleven", a matched-settings bullet.
- `appendices/results.tex`: pairwise paragraph (38 pairs; strict FAC 34, agreeing on 51 of 55; McNemar 45); Matched
  Settings paragraph (first five at matched settings, τ 0.709, readings are of default responses); Consistency Within
  a Template (`tab:single_path`); Scoring Rules (halving swaps none, doubling one non-separated pair; readable-only τ
  0.709); coverage: verbosity, the Claude–DeepSeek separation under four readings (route conformity), the judge-swap
  shifts beside their +0.023 MC gap, `tab:judged_steps`, `tab:flag_precision` (11 of 3,409 carried precision, none of
  171 confirmed slips); level gap: permutation holds for six, Holm note for `tab:level_gap`, GLM-5.3 readable-only
  −0.004; `tab:depth_model` paragraph (five of eleven, the four weakest among them).
- `appendices/taxonomy_content.tex`: domain selection in the owner's current wording (candidates from curricular
  documents and textbooks; three-LLM panel prompted individually; 15 domains, three per branch, and the areas by
  majority vote); panel's role stated once; `tab:area` sentence; textbook and data-source paragraphs rewritten; ART-DEIT
  and IEEE 100 dropped from `tab:authoritative_sources`; statistics caption "each branch holds 30 templates";
  difficulty paragraph names the domain experts' levels (D7), no agreement figure; straight quotes in prompt boxes.
- `appendices/template_examples.tex`, `docs/appendix_listings.py`: listings swapped (below); introduction rewritten;
  precision-policy sentence.
- `appendices/error_analysis.tex`: sample 137 (19 / 39 / 39 / 40) on 77 templates, default settings only; the
  no-error-removed shares introduced; Near Misses paragraph with the six templates named; incorrect verdicts and the
  answer's form (99 incorrect, 28 on three templates; prescribed digits; four of nine symbolic by numbers; rule
  sensitivities); the four read-off templates named; the two chemical templates' paragraph removed.
- `appendices/paraphrase.tex`: Mistral Large 3 cited; 314 / 275 / 114; two models outside the margin; four intervals on
  one side of zero; τ 0.807 is what sampling noise gives; the four repeat models named, 300 instances.
- `appendices/further_experiments.tex`: full-set changes (+0.038, +0.070, +0.118, none bounded within the margin);
  anchors current; open book +0.044; tool 66% / 18%, the "harder instances" claim dropped.
- `appendices/branch_domain.tex`: numbers only (lowest branch counts; the radar's top-tier lowest domains).
- `appendices/7_appendix.tex`: inputs `related_work` and `release` after the error analysis.
- `WRITING_RULES.md`: label list completed.
- `docs/appendix_certification.py`, `docs/appendix_evaluation.py`: checks for every new sentence (sources added to
  the docstring); the failed partial-verdict claim replaced by "0.5 is conservative on average".

## New labels (for WS-D1's pointers)

- `appendix:release` (release, query dates, scripts, contamination): Limitations and Ethics.
- `appendix:related_work`, table `tab:related_work` (fifteen closest benchmarks): §2.
- `appendix:readings` (experts' readings of the evaluated models): §4 or Limitations, if wanted.

## Bibliography

- Removed: `aiche2025constitution`, `ieee2025standards`, `asme2025vision`, `piepho2004letter` (uncited).
- Added under `% --- added by WS-D2` (last block of the file): `ncees2020fechemical`, `ncees2020fecivil`,
  `ncees2020feelectrical`, `ncees2020feindustrial`, `ncees2020femechanical` (NCEES FE CBT exam specifications,
  effective with the July 2020 examinations, checked against the PDFs); `mistral2025mistral3`; 33 related-work keys
  (published venue, location and pages from the ACL Anthology where they exist; papers accepted at venues not yet held,
  EMNLP 2026 and NeurIPS 2026, as arXiv preprints).
- §3.1 cites the five NCEES keys in place of the removed ones (WS-D1, done).
- `related_work.tex` also cites `golovneva2022roscoe`, `xia2024reasoneval`, `he2024olympiadbench`, `li2024gsmplus`
  from WS-D1's block; they must stay there.

## Hand-offs

Done by WS-D1: every generated phrase my appendix prose needed (pairs, strict FAC, tolerance, the Claude–DeepSeek
separation, the judge-swap gap and carried precision in `results.tex`, the paraphrase list, margin and τ, the
near-miss wording, the error sample, the anchors, the radar), the endpoint sentence in `paper_setup.py`, generated
blocks left out of its number count, the `tab:judged_steps` caption, "physically validated".

For the orchestrator (generated captions: an edit in `paper_results.py`, then `--write --text-only --headline default
--repaired`):
- `tab:scoring_variants`: "re-applied offline over the main store" → "re-applied to every response of the
  evaluation"; "no longer capped" → "not capped".
- `tab:providers`: "since the router assigns instances unevenly" → "since OpenRouter assigns instances to endpoints
  unevenly".
- `tab:errors`: "3 domain experts of the answer's branch read each, and 137 have been read (411 readings)" → "three
  domain experts of the answer's branch read each (411 readings)"; "over the 3 experts" → "over the three experts";
  drop "A reading counts only while the store holds the question and response as read."
- Sentence case for the bold leads: "Milestone Coverage behind wrong answers.", "Difficulty level and trace depth.",
  "Change under paraphrase.", "The four experiments against the evaluation.", "Final Answer Accuracy by branch and by
  domain.", "Final Answer Accuracy by engineering branch for four representative models.", "Final Answer Accuracy by
  domain for four representative models."
- `current_overleaf_project/main.tex` inputs `sections/7_appendix`; the appendix is `sections/appendices/7_appendix`
  (WRITING_RULES).

## Listings

- Easy: `template_undamped_natural_frequency_torsional` (mechanical). The question gives J and k_t and asks for ω_n,
  f_n and the period without stating the formula.
- Intermediate: `template_exponential_mttf_topology` (industrial), unchanged; the sampled topology changes the
  derivation.
- Advanced: `template_vdw_solve_for_volume` (chemical). a and b from the critical properties, coupled to the cubic
  equation of state, solved numerically for three real roots, redrawn into the two-phase region. Every instance has
  6 or more milestones; the script asserts this when `scores/milestones.json` is present. Its non-ASCII characters
  print through the listing's `literate` option, so the code is shown unchanged.
- Not chosen: the civil candidates (constant-head permeability, relative density) carry development notes in their
  code comments ("R1, cycle 1", "AUTHOR_NOTES", "D-016").

## FILL markers left

None.

## Check outcomes

| check | outcome |
|---|---|
| `docs/appendix_certification.py --check` | pass, 40 of 40 |
| `docs/appendix_evaluation.py --check` | pass: `5_evaluation.tex` 5 of 5, `scoring.tex` 32 of 32, `validation.tex` 64 of 64, no claim fails |
| `docs/appendix_statistics.py --check` | pass, 28 of 28, `tab:area` current |
| `docs/appendix_listings.py --check` | pass, 3 of 3 |
| `paper_results.py --check --headline default --repaired` | pass (0 failures; my five files hold every phrase) |
| `paper_setup.py --check` | pass (`models.tex` 45 of 45, no number outside the phrases) |
| labels, citations, FILL, vocabulary | every `\autoref` resolves, every key resolves, no FILL; script names only in `appendix:release` |

## Open items

- D10 is not recorded as chosen. §3.1 and the taxonomy appendix now say the same thing (candidates from the
  curricular documents and textbooks; a three-LLM panel; the 15 domains, three per branch, and the 42 areas by majority
  vote). The authors confirm the rule.
- D8 is not recorded as chosen. `release.tex`, Limitations and Ethics follow the recommendation (MIT License; seed at
  publication; held-out seed; names and per-rater files withheld). The owner confirms.
- Query dates (main 28 September to 6 October 2026, matched 6 to 7 October) come from `numbers_sheet.py`, which reads
  the local traces; no tex check covers them; the further experiments' dates are not in the sheet.
- `tab:authoritative_sources` still lists handbooks the constants tables do not cite (Perry's, ASM Handbook, Marks',
  Shigley's, Standard Handbook for Electrical Engineers). The owner decides whether to keep them.
- The round-1 mismatch split into real and planted items is not stated: `layer2/RESULTS.md` does not give it.
