# Follow-up review: the 18 main-text papers without a panel review, and a final pass

2026-10-02. The brief is `notes/followup_review_prompt.md`. The revised section cited 18 new papers in its
main text that no reviewer had read for anything beyond the sentence that cites them. The fact check
(`verify_facts.py`) had confirmed what the text says about each, but not whether the paper changes how
EngTrace should be positioned. This note gives each paper's verdict, what was found and what changed. It
then lists the final pass over the whole section.

**How.** One reviewer per paper, under the schema of `notes/review_protocol.md`. Papers 1 to 12 were
reviewed in the coordinating session. Papers 13 to 18 were reviewed by one reviewer agent (Sonnet).
Every claim of theirs that changed the text was re-checked against the paper's text before it was used.
The reviews are `reviews/candidate/<key>.r1.json`. Page numbers refer to `text/<key>.txt`. Two figures
were read from the PDFs because their text is not extracted: PRIME's Figure 3 and PRISM-Physics'
Figure 2.

## Summary

- **The positioning claim holds.** It was checked against Table A, the 100 cited works and the 86
  candidates that are not cited. None of these combines generated instances across engineering branches
  with gold traces, deterministic checks of intermediate values and an evaluator validated against step
  labels. The nearest are three. PRISM-Physics has a rule-based step check but static physics problems
  and solution-level validation. ChemCoTBench-V2 has deterministic step rules and expert agreement, but
  static chemistry. AtmosSci-Bench has expert templates and a rule-first grader checked against humans,
  but multiple-choice instances and no steps.
- **Ten of the 18 led to a change in the sentence that cites them or in their Table A row.** For two
  more, the section changed around them. Six needed no change.
  - **Changed:** PRIME, PRISM-Physics, AtmosSci-Bench, Akhtar et al., Długosz et al., Lunardi et al.,
    Krumdick et al., JudgeBench, MATH-Perturb and ControlBench.
  - **Changed around them:** the re-grading sentence (Ansari et al.) was cut and its paper cited in
    another sentence; ORQA led to OPT-Engine being added beside it.
  - **No change:** HiPhO, SymPyBench, Huang et al. (one appendix citation added), Pombal et al., MATH()
    and GSM1k.
- **Two statements about EngTrace itself were corrected** in the Positioning paragraph (Part 2, item 3).
- **Length.** The main text went from 723 words to 704; it is 678 without the remaining [cut first]
  sentence. The fact check went from 177 of 177 claims to 200 of 200.

## Part 1: the 18 papers

| # | Paper | v3 sentence | Change |
|---:|---|---|---|
| 1 | PRIME, Wang et al. 2026 | imprecise | "in mathematics and engineering" becomes "in college-level STEM including engineering" |
| 2 | PRISM-Physics, Zhao et al. 2026 | accurate | "symbolic" becomes "rule-based"; Table A row sharpened |
| 3 | HiPhO, Yu et al. 2026 | accurate | none |
| 4 | SymPyBench, Imani et al. 2025 | accurate | none |
| 5 | AtmosSci-Bench, Li et al. 2025 | the sentence accurate; the Table A row imprecise, one cell wrong | Table A row corrected |
| 6 | Akhtar et al. 2026 | accurate | the clause now carries H6, at no extra length; Introduction citation corrected (CHANGES.md) |
| 7 | Długosz et al. 2026 | imprecise | "the drops they report" becomes "GSM-Symbolic's drops" |
| 8 | Lunardi et al. 2025 | imprecise | scoped "on multiple-choice benchmarks" |
| 9 | Huang et al. 2026 | accurate | cited also in the appendix, for the rule-first hybrid |
| 10 | Krumdick et al. 2026 | imprecise | "a reference answer" becomes "a correct reference" |
| 11 | Pombal et al. 2026 | accurate | none (a note on `JUDGE_SELECTION.md` below) |
| 12 | JudgeBench, Tan et al. 2025 | imprecise | "they do" becomes "many do" |
| 13 | Ansari et al. 2026 | imprecise | the [cut first] sentence cut; cited instead for rule-based grader errors |
| 14 | MATH-Perturb, Huang et al. 2025 | imprecise | "the original instance" becomes "the original instance or its method" |
| 15 | MATH(), Srivastava et al. 2024 | accurate | none |
| 16 | GSM1k, Zhang et al. 2024 | accurate (the reviewer said imprecise; see below) | none; a fact entry added |
| 17 | ControlBench, Kevian et al. 2024 | imprecise | "experts graded the reasoning" becomes "experts examined the reasoning" |
| 18 | ORQA, Mostajabdaveh et al. 2025 | accurate | OPT-Engine added beside it |

### The closest precedents

**1. PRIME.** A benchmark of verifiers, not of solvers. It has 2,530 model responses to textbook and
exam problems. Eighteen domain experts gave each response two binary labels: outcome correct, and
derivation sound and reaching the answer (p.2, p.5-6). Verifiers are scored on the second.

- **How much is engineering.** About a fifth. Figure 3 (p.5, an image) gives mathematics 36%, chemistry
  28%, biology 19% and engineering 18%. The paper contradicts itself here: its results table has a
  physics column and no engineering column (p.6).
- **The expert labels** are per response, not per step, and no agreement statistic is reported.
- **Table A.** PRIME does not belong there: it has static items, no gold steps, and it measures
  verifiers.
- **The change.** The sentence keeps its point (verifiers must catch answers the reasoning does not
  support) but now says "college-level STEM including engineering", the paper's own term (p.1, p.4).

**2. PRISM-Physics.** Reference solutions become graphs of formulas. A model's formulas are matched to
them by a rule-based equivalence test, constant substitution and then random numerical substitution
(p.6, p.16), and every matched node's ancestors are credited (p.18). The authors call it "our non-LLM
PRISM-DAG" (p.9). That makes it the nearest physics precedent for EngTrace's deterministic step checks,
which the word "symbolic" did not convey. The main text now says "rule-based formula matching".

Three findings go into Table A:

- **Source.** The problems come from a book of PhD qualifying-exam questions (p.19). The rewriting into
  graphs is done by an LLM pipeline with rule and LLM checks (p.6, p.19).
- **The human check.** It covers 70 DeepSeek-V3 solutions, each scored by two experts with a third
  adjudicating (p.9). It is a per-solution score, and the experts are mostly co-authors (p.31).
- **The 0.294 baseline** is an outcome-only LLM judge giving 0/1 scores (p.9).

Its size is not stated in the text: Figure 2 shows topics without counts. Ansari et al. (item 13) found
that grader errors dominate in their audit of PRISM-Physics, because of its rule-based evaluator (p.4).
That is the rule-based failure mode the evaluator paragraph cites.

**3. HiPhO.** It has 13 olympiad exams, 360 problems and 519 subquestions (p.4). Only 7 exams have
marking schemes, and step grading applies only to those. A problem scores the maximum of its answer and
step scores (p.6), so a correct final answer gets full credit whatever the steps.

The LLM grader (Gemini-2.5-Flash) was compared with human graders only in one illustrative table: IPhO
2024, two contestant models, one run, and no agreement statistic (Table 7, p.19). "Fully aligned with
human examiners" (p.1) refers to using the official schemes, not to a measured agreement. v3 claims no
validation for HiPhO, which is right, so there is no change.

**4. SymPyBench.** Every Table A cell holds. LLaMA-3.2-90B wrote the templates, paraphrases and Python
solvers from Creative Commons problems; about 88% got a solver that passed validation, and all problems
were reviewed by hand (p.3-5). The benchmark scores answers (partial and exact match) and three
cross-variant metrics. It does not score steps: the step-by-step reasoning is a reference no metric
uses (p.4-6). How free-form answers are matched is not described.

Neither the paper (a Meta report dated 8 Dec 2025) nor its arXiv page names a venue. The search note's
"EACL 2026 (Industry Track)" is nevertheless right: the ACL Anthology has the paper as 2026.eacl-industry.8,
pp. 105-118 (see the venue lookups). The bibliography now cites that version.

**5. AtmosSci-Bench.** The generation sentence is right, but the Table A row was not.

- **Its generated instances are multiple choice.** They come from 67 expert templates, with GPT-4o
  solvers that experts reviewed, and are scored by the option chosen (p.5-6).
- **The 391 open-ended questions are static.** A cascade scores them: a numeric check within 5%, then
  SymPy, then GPT-4o-mini (p.7).
- **"–" for human validation was wrong.** Appendix J reports 92.79% and 93.02% agreement between the LLM
  grader and human graders on the open-ended answers it decided (p.31).
- **No reference steps** are scored.

The row now says all of this. The rule-first cascade, checked against humans, is the closest
answer-level analogue of EngTrace's design in the cited set.

### Claims the framing leans on

**6. Akhtar et al., benchmark saturation (ICML 2026, per the arXiv comment).** Their hypothesis H6 is
answered on p.6: "templated benchmarks (N = 14) do not differ significantly from non-templated ones
(N = 46) in saturation behaviour (p = 0.10)". H1, on private test sets, is rejected, though on 4 private
benchmarks against 56 public (p.5-6).

- **What it decides.** EngTrace may not present templates or the private seed as a saturation
  safeguard. The private seed supports a contamination claim: no evaluated item was public when the
  models ran.
- **The change.** The clause now says "neither private test sets nor templating has measurably slowed
  saturation". It has the same length as before, and both halves are in `facts.py`.
- **BIG-Bench Hard.** The paper counts BBH among the benchmarks that "remain unsaturated despite
  prolonged exposure" (p.6). The Introduction replacement in CHANGES.md had cited it in the bracket for
  "as GLUE and BIG-Bench Hard did". That bracket is now split: Akhtar et al. for saturation in general,
  BBEH for BBH (CHANGES.md, section 2).

**7. Długosz et al., re-evaluation of GSM-Symbolic.** Using GLMMs with per-question random effects,
they find significant changes for only 8 of 20 models. The variants' integers run larger than GSM8K's,
which accounts for half of the remaining effects (p.1). They apply Holm across models (p.5).

- **The change.** The paper re-analyses GSM-Symbolic only, so "the drops they report" became
  "GSM-Symbolic's drops", which also saves two words.
- **Does the critique apply to EngTrace?** Mostly, EngTrace's analyses already meet it.
  - Every interval resamples templates, and model differences use a sign-flip permutation over
    per-template differences, because "instances of a template are not" independent
    (`full_run_28092026/results/RESULTS.md`, lines 3 and 23; D-111).
  - Holm runs across pairs and across models.
  - The paraphrase arm keeps the original's numbers as a multiset, so no magnitude shift can arise
    there.
- **What stays exposed** is the difficulty cliff. It compares different templates in different tiers,
  so tier is confounded with whatever else differs between them. RESULTS.md checks two such factors
  (line 101). The paper should not read the gap as an effect of difficulty alone.
- **Venue.** EMNLP 2026 rests on the arXiv comment ("Accepted to EMNLP 2026 (main conference), track:
  Resources and Evaluation"). The acknowledgements thank reviewers of "a 2026 ARR cycle" (p.10). There
  is no proceedings record yet.

**8. Lunardi et al., paraphrase robustness (ECAI 2025, per the arXiv comment).** GPT-4o mini wrote five
paraphrases of each of 53,966 questions from six multiple-choice benchmarks, with no human check (p.3-4).
Accuracy drops, while Kendall's τ between rankings stays above 0.9 (p.6-7).

- **The mismatch.** v3 stated this as a general fact, and EngTrace's own arm reproduces neither half. No
  model's paraphrase difference survives Holm, and τ is 0.636, the 5th percentile of the noise
  distribution at the arm's size (`full_run_28092026/PARAPHRASE_PAPER_NOTES.md`).
- **The change.** The clause is scoped "on multiple-choice benchmarks".
- **For the robustness section, within the notes' rules:**
  - On multiple-choice benchmarks with unchecked paraphrases, scores drop while rankings hold.
  - On EngTrace's expert-checked paraphrases, no drop survives correction. The tiers hold; the order
    within the top tier is not resolved.
  - The experts rejected 12% of the paraphrases that had passed every script check.
  - The banned phrasings stay banned: "robust to paraphrasing", "the ranking is stable" and "no
    contamination".

**9. Huang et al., rule- and model-based verifiers.** Rule-based verifiers have precision above 99% but
recall of about 86% (p.1, p.3), so the v3 sentence is accurate.

- **The hybrid.** The paper's rule-first hybrid asks a model only about responses the rules reject. It
  adds about 3 points of recall at over 98% precision (p.4). That is EngTrace's design: the
  deterministic milestone check first, MiMo only on what it missed (`full_run_28092026/judge.py`). The
  appendix sentence on escalation now cites it for this.
- **Reward hacking.** Model verifiers can be hacked, but that result concerns a policy trained against
  the verifier. EngTrace's models are not trained against its judge.
- **For the paper.** The judge's measured precision against experts belongs beside it (milestone
  precision 0.930, recall 0.989; RESULTS_X1 via the state brief).
- **Venue.** EMNLP 2026 rests on the arXiv comment alone.

**10. Krumdick et al., No Free Labels (COLM 2026, the paper's own header).**

- **The change.** The condition that matters is a correct reference. A reference the judge writes
  itself does not help, and "slightly incorrect references can be worse than providing no reference at
  all" (p.2, p.6, p.10). v3's "a reference answer" became "a correct reference".
- **The judge's inputs.** By analogy, the paper supports what `judge.py` (lines 9 to 12) sends: each
  missed milestone's id and value, not the gold solution. The values are code-computed and the gold is
  verified (2,250 of 2,250), so the judge never depends on solving the problem itself.
- **The limits of the analogy.** The paper's references are full answers, and its task is judging
  whole responses.

**11. Pombal et al., self-preference on rubrics (COLM 2026, the paper's header).** Accurate.
`evaluator_pilot_17092026/JUDGE_SELECTION.md` (lines 51 to 53) cites the paper in the right direction,
and with its family finding (p.2). But it turns "can be more than 50% more likely" (p.1-2), the largest
effect observed, into "judges are over 50% more likely". If the paper repeats the figure, write "can
be". That file is outside this folder and was not edited.

**12. JudgeBench (ICLR 2025, the paper's header).** GPT-4o scores 56.57% where chance is 50%. But in the
same table, o3-mini (high) reaches 80.86%, o1-preview 75.43% and DeepSeek-R1 73.14% (p.7-8). "They do
little better than chance" became "many do", the abstract's own scope ("many strong models", p.1).

### The rest

**13. Ansari et al., expert re-grading of physics benchmarks.** The cut sentence's numbers were right.

- **The figures.** 143 benchmark errors, 95 grader errors and 12 model errors make 250 (p.13). The
  re-graders were faculty and graduate researchers (p.1-2). The 250 are questions on which all of one
  model's attempts were marked wrong, on four audited benchmarks (p.3).
- **Why it was cut.** It was the first [cut first] sentence, and cutting it pays for the additions.
- **Where the paper is now cited.** Its finding that grader errors arise "mainly with rule-based
  evaluators, which can fail to recognize a correct answer written in an equivalent mathematical form"
  (p.3-4) now sits beside Huang et al. in the evaluator paragraph.
- **Venue.** None stated (arXiv 2609.13009v1).

**14. MATH-Perturb.** The simple perturbations test memorisation of the instance. The hard ones test
whether a model "memorizes the problem-solving techniques ... and blindly applies them" (p.2). The
sentence now reads "memorised the original instance or its method". That distinction matches
EngTrace's own paraphrase notes, which separate memorised wording from memorised methods.

The paper has 279 problems in each split (p.1). The text and the arXiv page name no venue, but
OpenReview lists it as an ICML 2025 poster.

**15. MATH().** Accurate: 2,060 of MATH's 5,000 problems are functionalised (p.7), and the reasoning gap
runs from 58.35% to 80.31% (p.3). No venue is stated (arXiv 2402.19450v1), so it stays a preprint.

**16. GSM1k (NeurIPS 2024 Datasets and Benchmarks, per p.1 and the arXiv comment).** The reviewer rated
"leak into training data" imprecise, because GSM1k measures overfitting (drops of up to 8%, p.1).

That is not adopted here. GSM1k itself reads the overfitting as leakage: "one important component of
this overfitting is that models have partially memorized examples from GSM8k" (p.3). It bases this on
a likelihood correlation (Spearman's r² = 0.36, as printed, p.1). The evidence is indirect, but the
claim is the paper's own, and Deng et al. carries the clause with it. No change; a fact entry now
records the passage.

**17. ControlBench.** It has 147 undergraduate control problems, and ControlBench-C is a multiple-choice
subset (p.2, p.5). The experts judged each answer right or wrong and sorted the failures into seven
error types (p.11). They did not score the reasoning.

"Experts graded the reasoning in TransportBench and ControlBench" therefore said too much. It now reads
"experts examined the reasoning", which fits both papers. The matching cell of CITED_SET_AUDIT.md for
TransportBench was updated. No venue is stated (arXiv 2404.03647v1).

**18. ORQA (AAAI 2025: the paper's copyright footer and arXiv comment).** It has 1,513 four-option
questions on formulating optimisation models (p.2, p.4), so it scores the option chosen. It is a weak
match for EngTrace's industrial branch, which covers production and inventory, quality and
reliability, and stochastic operations.

Among the candidates that were not cited, OPT-Engine is a better one (ICML 2026: "Proceedings of the
43rd International Conference on Machine Learning ... PMLR 306", p.1). Its programmatic generators
produce instances of ten operations-research problem classes, among them inventory and production,
with verifiable optimal solutions (p.1-4). It is now cited:

- beside ORQA in the single-branch list;
- in the sentence on benchmarks that generate their instances, where a reader looking for every
  generated engineering benchmark would expect it.

A2utoLPBench is optional, and SupChain-Bench (tool orchestration) is out of scope. ERI, already cited,
lists industrial engineering among its nine fields (p.1).

## Part 2: the final pass

1. **Accuracy.** `facts.py` changed in two ways.
   - Five claims the text no longer makes were removed: four from the cut sentence and AtmosSci-Bench's
     "final answers".
   - Twenty-eight were added, one or more for each new or changed statement.
   - `verify_facts.py --quiet` finds 200 of 200.
2. **Positioning.** Two things were checked against Table A, the 100 cited works and the 86 uncited
   candidates (`NEW_PAPERS_AUDIT.md`): the "To our knowledge, none combines ..." sentence, and the
   Positioning paragraph.
   - Six uncited candidates were read closely because their titles suggested they might combine some of
     the four elements: engineering-equation solving, VEHBench, the Infinite Problem Generator, VPRMs,
     SciR and EngIntervene. None combines generated instances across engineering branches with gold
     traces, deterministic intermediate checks and expert step labels.
   - The claim holds. The section makes no "first" or "novel" claim: every "first" in the body is the
     [cut first] marker or "first-error".
3. **EngTrace's own facts.** Two statements in the Positioning paragraph were wrong in detail.
   - **The arithmetic check.** "Intermediate values and arithmetic are checked deterministically against
     that trace" put the arithmetic check against the gold trace. The digit rule (E4) instead recomputes
     each claim from the numbers the model's own trace shows (state brief, section 2). It now reads
     "Intermediate values are checked deterministically against that trace and arithmetic by
     recomputation".
   - **The judge's scope.** "Decides only the milestones these checks leave open" left out the step
     router. The router is in the evaluator stack (D-110, D-113): steps the digit rule does not flag go
     to the same judge, and its run is in progress (D-164). The text now says "decides only what these
     checks leave open, given each milestone's expected value".
   - **The rest matches its sources:**
     - the private 128-bit seed: `full_run_28092026/README.md`;
     - the judge's inputs: `judge.py`, lines 9 to 12;
     - 15 experts and 300 traces, three own-branch experts per trace: `RESULTS_X1.md`, lines 3 and 4;
     - the PRM baselines: `RESULTS_E2.md`;
     - the expert-checked paraphrases: `PARAPHRASE_PAPER_NOTES.md`.
4. **Consistency across the folder.**
   - **CITED_SET_AUDIT.md:** TransportBench's "In v3" cell now says "examined". No other cell describes a
     changed sentence. The BBH and GLUE rows already cite BBEH and SuperGLUE, not Akhtar et al., for
     saturation.
   - **CHANGES.md:** updated for every change above, including the Introduction bracket and the counts:
     100 works, 68 in the main text, 32 in the appendix (31 from the search, plus NEWTON), and 15 Table A
     rows besides EngTrace's. It also corrects the venue list in section 6.
   - **README.md:** the counts, the status line, and a row for this review.
   - **NEW_PAPERS_AUDIT.md:** regenerated (39 in the main text, 31 in the appendix, 86 not cited). Its
     generator now mentions this review.
   - **`notes/review_tally.md`:** regenerated.
   - **DECISIONS.md:** D-167 records the substance. D-166 is unchanged.
5. **Bibliography.**
   - **Keys.** `render_related_work.py` reports all 100 cited keys in the .bib.
   - **Spot check.** Ten entries were checked against their first pages: Długosz, Düzkar, SuperGPQA,
     Krumdick, Pombal, JudgeBench, GSM1k, Lunardi, ORQA and PRIME, plus OPT-Engine. Authors and years
     match in all. Długosz's second author is "Arlindo Oliveira" in the record and "Arlindo L. Oliveira"
     on the paper; that is harmless. The venues come from:
     - the page header or footer for four;
     - the arXiv comment for three;
     - for PRIME, the ACL Anthology (2026.acl-long.683); for SuperGPQA, the catalogue (its NeurIPS 2025
       entry was audited with the May citations).
   - **Venues.**
     - **ERI** is now confirmed. The DOI on its arXiv page, 10.1016/j.cie.2026.112333, is Computers &
       Industrial Engineering volume 221 in Crossref (`check_venues.py`, `notes/venue_check.json`). It
       has left `make_bib.py`'s single-source list.
     - **Huang et al. and Długosz et al. (EMNLP 2026)** still rest on their arXiv comments, as no
       proceedings record exists yet.
     - **Two preprints had appeared.** SymPyBench is in the EACL 2026 Industry Track, and AtmosSci-Bench
       in the NeurIPS 2025 Datasets and Benchmarks Track. Both are now cited at those venues, through
       `make_bib.py`'s `VENUE_CONFIRMED`; SymPyBench's label becomes "Imani et al., 2026".
     - The other uncertain venues are confirmed (see the lookups below). DBLP answers scripts with a bot
       check, so OpenReview, Semantic Scholar, Crossref and the ACL Anthology were used.
   - **Characters.** The .tex is pure ASCII. The .bib holds twelve accented letters (é, í, á, ü, ł, Š,
     ă), all covered by LaTeX's UTF-8 input under the T1 encoding the ACL template loads. "M-A-P Team" is
     brace-protected and prints as "M-A-P Team et al., 2025". This is a static check: no TeX
     installation is available on this machine, so nothing was compiled.
6. **Readability and length.**
   - **Length.** The main text is 704 words, against 723 before, as `render_related_work.py` counts
     them. It is 678 without the remaining [cut first] sentence.
   - **Paragraphs.** Each opens with its point. The generation paragraph now opens with what does not
     slow saturation, which is the argument for the private seed and the step checks.
   - **Cutting further.** One [cut first] sentence remains (design, multiphysics, industrial practice).
     The header's cut order is unchanged.
7. **Anonymity.**
   - FinChain is cited in the third person only.
   - The .tex and the section body mention neither EngTrace's arXiv preprint (2511.01650, "EngChain")
     nor that ERI, Sci-ρ or LPDS cite it. Those facts appear only in the notes for the authors and in
     CHANGES.md.

## Venue lookups

`check_venues.py` (`notes/venue_check.json`) covers the twelve keys whose venue was uncertain. It
fetches all of Semantic Scholar's records in one batch request and searches OpenReview by title. It
also reads Crossref for ERI's DOI and the ACL Anthology for two ids found on its event pages.

| Paper | Venue in the bibliography | Evidence | Status |
|---|---|---|---|
| ERI | Computers & Industrial Engineering 2026 | Crossref: the DOI on its arXiv page is volume 221; Semantic Scholar agrees | confirmed |
| PRIME | ACL 2026 | ACL Anthology 2026.acl-long.683, pp. 14973-14985; Semantic Scholar | confirmed |
| PRISM-Physics | ICLR 2026 | OpenReview: ICLR 2026 Poster | confirmed |
| HiPhO | ICML 2026 | OpenReview: ICML 2026 (rejected at ICLR 2026 first) | confirmed |
| MATH-Perturb | ICML 2025 | OpenReview: ICML 2025 poster; Semantic Scholar | confirmed |
| ORQA | AAAI 2025 | arXiv comment; the PDF's AAAI copyright footer; OpenReview's import | confirmed |
| SymPyBench | EACL 2026 Industry Track (was a preprint) | ACL Anthology 2026.eacl-industry.8, pp. 105-118; Semantic Scholar | confirmed, bibliography corrected |
| AtmosSci-Bench | NeurIPS 2025 Datasets and Benchmarks (was a preprint) | OpenReview: D&B Track poster; Semantic Scholar | confirmed, bibliography corrected |
| Huang et al. 2026 | EMNLP 2026 | the arXiv comment only; OpenReview shows an ICLR 2026 withdrawn submission | single source |
| Długosz et al. 2026 | EMNLP 2026 | the arXiv comment only; OpenReview shows an ARR May 2026 submission | single source |
| MATH() | preprint | no venue in any index | preprint |
| ControlBench | preprint | no venue in any index | preprint |

AtmosSci-Bench's title in the bibliography comes from the arXiv record ("... the Recent Advance of
Large Language Model ..."). The paper's title page reads "Recent Advances of Large Language Models".
Use the proceedings title when the reference list is regenerated.

## Outside this folder: for the authors, not edited

- **`evaluator_pilot_17092026/JUDGE_SELECTION.md`, lines 51 to 53.** Write "can be more than 50% more
  likely", not "are over 50% more likely" (item 11).
- **The difficulty cliff (RESULTS.md Q2).** Present it as a difference between the templates in each
  tier, not as an effect of difficulty alone (item 7).
- **The robustness section.** For the contrast with Lunardi et al. (item 8), use the notes' wording.

## Not resolved

- **Venues.** EMNLP 2026 for Huang et al. and Długosz et al. still rests on their arXiv comments; no
  proceedings record exists yet. Every other uncertain venue is now confirmed (see the lookups above).
- **Missing details in the papers.**
  - PRISM-Physics' size is not stated in its text.
  - The names of ControlBench's seven error types are in a figure that is not extracted.
  - SymPyBench's rule for matching free-form answers is not described.
- **The LaTeX was not compiled.**
