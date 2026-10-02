We are revising the Related Work section of EngTrace, our benchmark paper for the ARR submission that closes on 12 October 2026. The repository is C:\Users\ayesha.gull01\EngTrace, on branch master; stay on master. All of this work lives in docs/related_work_oct2026/.

An earlier session audited every paper the May 2026 submission cites, using three independent full-text reviews per paper. It then searched for new work and rewrote the section. The revised section cites 69 new papers, 38 of them in the main text. Eighteen of those 38 were never properly reviewed. The only things behind them are a short description written by the search agent that found each one, plus a claim check: verify_facts.py matches each statement the text makes about a paper to a passage in that paper's text. That confirms what we say, but nobody has read these papers looking for anything that would change how EngTrace should be positioned. Your job is to review those 18, then make one final pass over the whole section and make sure it all makes sense.

## Read first

1. docs/related_work_oct2026/README.md: how the folder is organised and how each file is produced.
2. docs/related_work_oct2026/related_work_v3.src.md: the section's source. It is the only file you edit. RELATED_WORK_v3.md and related_work_v3.tex are generated from it.
3. docs/related_work_oct2026/CHANGES.md: what changed from the May text and why.
4. docs/related_work_oct2026/notes/engtrace_state_brief.md: what EngTrace now is. It has 150 templates in five branches and 2,250 items drawn from a private seed. Its evaluator combines a deterministic answer check, milestone coverage, a per-claim arithmetic check, and a residual LLM judge from outside the evaluated roster. That evaluator was validated against 15 experts on 300 traces.
5. docs/related_work_oct2026/notes/review_protocol.md: the review schema, in the section "Per-paper review file".
6. docs/related_work_oct2026/facts.py and notes/fact_check.md: every claim the text makes, with its evidence.
7. docs/related_work_oct2026/CITED_SET_AUDIT.md: an example of the depth and tone expected.

## Part 1: review the 18 papers

Each paper's text is in docs/related_work_oct2026/text/<key>.txt, with `===== PAGE n =====` markers for page references. The PDFs are in papers/. For each paper:

- Read the abstract, introduction, benchmark or method construction, and evaluation sections in full, and skim the rest. Use grep for section headings so you read what matters rather than the whole file.
- Write docs/related_work_oct2026/reviews/candidate/<key>.r1.json in the protocol's schema. No review exists yet for these keys; check with ls first. Put the v3 sentence that cites the paper in `engtrace_use.sentence`, and grade it on the protocol's scale: accurate, imprecise, misleading or wrong. Write the file with the Write tool, not a shell heredoc, which fails on apostrophes here. Keep the file ASCII, and validate it with python json.load.
- Look beyond the sentence. Does the paper score intermediate steps, generate instances, give gold traces, validate its scorer against experts, or cover engineering? If so, the positioning claim may need to change. That claim is "To our knowledge, none combines generated instances across engineering branches with gold traces, deterministic checks of intermediate values, and an evaluator validated against experts' step-level labels."

The 18, highest priority first. Each line gives where the paper is cited and what v3 says about it.

**Closest precedents:**

1. **wang2026prime** (PRIME, ACL 2026)
   - Cited in: the process paragraph.
   - v3 says: it "tests whether verifiers catch answers that the reasoning does not support, in mathematics and engineering."
   - Find out: how much of it is engineering (the search could not tell, and its results tables show math, physics, chemistry and biology); whether its expert labels are per solution or per step; and whether it belongs in Table A.
2. **zhao2026prismphysics** (PRISM-Physics, ICLR 2026)
   - Cited in: paragraph 1, and Table A.
   - v3 says: it "score[s] steps against reference solutions ... [by] symbolic formula matching". Table A gives: physics | static | formula graphs | steps, by symbolic formula matching | Kendall τb 0.346 with expert annotations (an LLM judge: 0.294).
   - Find out: whether the matching is fully deterministic; the size of the benchmark; how the expert annotations were made. It is the nearest physics precedent for deterministic step checks. Should v3 credit it more explicitly?
3. **yu2026hipho** (HiPhO, ICML 2026)
   - Cited in: paragraph 1.
   - v3 says: it scores steps with "official marking schemes applied by an LLM judge".
   - Find out: whether that LLM application of the marking schemes was checked against human graders.
4. **imani2025sympybench** (SymPyBench, arXiv 2512.05954)
   - Cited in: the generation paragraph, and Table A.
   - v3 says: it "instantiate[s] science problems from code templates". Table A gives: university physics | parameterised Python | step-by-step reasoning | answers, and consistency across variants | –.
   - Find out: whether every cell is right, especially whether it scores steps; and its venue (one search note says EACL 2026 Industry Track).
5. **li2025atmosscibench** (AtmosSci-Bench, arXiv 2502.01159)
   - Cited in: the generation paragraph, and Table A.
   - Table A gives: atmospheric science | 67 templates with Python solvers | – | final answers | –.
   - Find out: whether it gives reference steps; how its numeric, SymPy and LLM grading cascade works; and its venue.

**Claims the framing leans on:**

6. **akhtar2026benchmarksaturation** (ICML 2026)
   - Cited in: the generation paragraph, and the Introduction replacement in CHANGES.md.
   - v3 says: static test sets saturate, and "safeguards such as private test sets have had limited effect on saturation".
   - Find out: what its hypothesis H6 found for templated versus non-templated benchmarks. The earlier session could not confirm this, and it decides how EngTrace may talk about saturation.
7. **dlugosz2026gsmsymbolicreeval** (EMNLP 2026, per the arXiv comment only)
   - v3 says: the drops variant benchmarks report "do not always survive re-analysis".
   - Find out: whether its statistical critique also applies to EngTrace's own robustness and complexity analyses; and confirm the venue.
8. **lunardi2025paraphraserobustness** (ECAI 2025)
   - v3 says: "paraphrasing lowers absolute scores while leaving model rankings relatively stable".
   - Find out: how this sits beside EngTrace's own paraphrase result in full_run_28092026/PARAPHRASE_PAPER_NOTES.md. Respect that file's banned phrasings.
9. **huang2026verifierrobustness** (EMNLP 2026, per the arXiv comment only)
   - v3 says: "Rule-based checkers avoid these biases but can reject correct answers given in unexpected formats."
   - Find out: what it implies for EngTrace's deterministic answer check and its model-based residual judge; and confirm the venue.
10. **krumdick2026nofreelabels** (COLM 2026)
    - v3 says: judges agree "with experts mainly on questions they can answer themselves unless given a reference answer".
    - Find out: whether it supports giving the residual judge each milestone's expected value (full_run_28092026/judge.py, lines 9 to 12).
11. **pombal2026rubricselfpreference** (COLM 2026)
    - v3 says: self-preference persists "even on objective rubric items".
    - Find out: whether this agrees with how evaluator_pilot_17092026/JUDGE_SELECTION.md already cites the paper.
12. **tan2025judgebench** (JudgeBench, ICLR 2025)
    - v3 says: "on difficult objective comparisons they do little better than chance".

**The rest:**

13. **ansari2026physicsregrading** (arXiv 2609.13009)
    - v3 says: "physicists re-grading 250 answers marked wrong on four physics benchmarks attributed all but 12 of them to errors in the benchmark or the grader". This is a sentence marked [cut first].
14. **huang2025mathperturb** (MATH-Perturb, ICML 2025)
    - v3 says: functional, symbolic and perturbed variants "test whether a model has memorised the original instance".
15. **srivastava2024functionalbench** (MATH(), arXiv 2402.19450)
    - Cited in: the same sentence as MATH-Perturb.
    - Find out: its venue.
16. **zhang2024gsm1k** (GSM1k, NeurIPS 2024 Datasets and Benchmarks)
    - v3 says: static test sets "leak into training data".
17. **kevian2024controlbench** (ControlBench, arXiv 2404.03647)
    - v3 says: it is the control benchmark, and "experts graded the reasoning in TransportBench and ControlBench".
18. **mostajabdaveh2025orqa** (ORQA, AAAI 2025)
    - v3 says: it is the operations-research benchmark, and the only reference for the industrial branch.
    - Find out: its format (multiple choice or not); whether v3 should say what it scores; and whether a better industrial-engineering reference exists among the candidates in NEW_PAPERS_AUDIT.md.

**Cost.** Do these reviews yourself in this session. If you need to split the work, use at most three subagents, one reviewer per paper, on Sonnet. Never use Fable, and do not run multi-reviewer panels. The earlier session hit its usage limit on panels.

## Part 2: one final pass over the section

Check each of these, and fix at the source:

1. **Accuracy.** Every statement about another paper must be in facts.py and pass `python docs/related_work_oct2026/verify_facts.py --quiet`. Add or adjust entries for anything Part 1 changes.
2. **Positioning.** The "To our knowledge, none combines ..." sentence and the Positioning paragraph must hold against Table A, all 99 cited works, and the 87 candidates in NEW_PAPERS_AUDIT.md that are not cited. The section makes no "first" or "novel" claims, and must not start.
3. **EngTrace's own facts.** Every statement about EngTrace must match notes/engtrace_state_brief.md and the files it cites:
   - private seed: full_run_28092026/README.md
   - residual judge and its inputs: full_run_28092026/judge.py
   - the 15 experts and 300 traces: evaluator_pilot_17092026/RESULTS_X1.md
   - PRMs as baselines: evaluator_pilot_17092026/RESULTS_E2.md
   - expert-checked paraphrases: full_run_28092026/PARAPHRASE_PAPER_NOTES.md
4. **Consistency across the folder.** These must agree with the final text and with each other:
   - CITED_SET_AUDIT.md (the "In v3" column)
   - CHANGES.md (its counts: 99 works, 67 in the main text and 32 in the appendix, 15 Table A rows, the claim total; and the replacement sentences' keys)
   - README.md
   - D-166 in docs/re-implementation-sep/DECISIONS.md
5. **Bibliography.** `python render_related_work.py` must report every cited key present in the .bib.
   - Spot-check ten entries' authors, years and venues against the papers' first pages.
   - Try to confirm the three venues that rest on one source: Huang 2026 (EMNLP), Długosz 2026 (EMNLP) and ERI (Computers & Industrial Engineering).
   - Check that the LaTeX characters render: ł in Długosz, ü in Düzkar, the M-A-P Team author.
6. **Readability.** Each paragraph should open with its point and read as an argument, not a list. The main text is 723 words of prose. Do not let it grow; if you add, cut elsewhere, starting with the sentences marked [cut first].
7. **Anonymity.** FinChain shares authors with us and is cited only in the third person. The section must not reveal that EngTrace has a public arXiv preprint (2511.01650), or that ERI, Sci-ρ or LPDS cite it.

## How to edit

All commands run from the repository root.

- Edit only docs/related_work_oct2026/related_work_v3.src.md. Citations are `[@key]` (parenthetical) or a bare `@key` (textual). Then run `python docs/related_work_oct2026/render_related_work.py`, which rewrites RELATED_WORK_v3.md and related_work_v3.tex.
- To add a citation, add its key to related_work_keys.txt under "# Main text" or "# Appendix". Then run `python docs/related_work_oct2026/make_bib.py --keys docs/related_work_oct2026/related_work_keys.txt`. It caches raw records, so only new ones are fetched.
- After any change, run:
  - `python docs/related_work_oct2026/verify_facts.py --quiet` (all claims must be found)
  - `python docs/related_work_oct2026/make_new_papers_audit.py`
  - `python docs/related_work_oct2026/tally_reviews.py`

## Deliverables and rules

- **Your note.** Write docs/related_work_oct2026/notes/followup_review.md. Give each of the 18 papers its verdict, what you found, and any change to the text. Then list your final-pass findings and how you fixed each. Update the status line in README.md.
- **Decision log.** If the section's claims change in substance, add a new D-entry to docs/re-implementation-sep/DECISIONS.md. The log is append-only, so do not edit D-166. Check the latest D number first, because another session also writes there.
- **Commits.** Commit messages are one short line, with no body and no Co-Authored-By trailer. Push to origin master immediately after every commit.
- **Scope.** Stay inside docs/related_work_oct2026/, apart from that DECISIONS.md entry. Other sessions share this machine and a long evaluation run may be in progress, so do not inspect or stop processes you did not start. Do not publish artifacts, and make no paid API calls.
- **Final report.** End with a short message: which of the 18 needed changes, what changed in the section, and anything you could not resolve.
