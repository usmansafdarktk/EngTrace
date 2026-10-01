# Search notes: physics_science

Topic: physics and physical-science reasoning benchmarks, and general science or graduate-level suites, 2024 to October 2026.
Searched 2026-10-01 by the physics_science agent. Output catalogue: `candidates/physics_science.json` (30 entries: 8 must-cite, 15 should-cite, 7 optional). All 30 PDFs are on disk (27 fetched by this agent, 3 already fetched by symbolic_contamination under the same keys; see "Download record").

Facts below come from arXiv abstract pages, and for the must-cites also from the arXiv HTML full text. WebFetch summarised both with a small model, so the numbers here still need page references. The review panel should confirm every number against `text/<key>.txt` before the rewrite uses it. Where a fact came from a secondary page, the note says so.

## 1. The positioning question

### Which physics benchmarks evaluate intermediate steps?

- **Cited set:**
  - PhysReason scores steps (PSAS-S, alongside the answer-level PSAS-A).
  - PHYBench gives expression-level partial credit (the EED score) on the final expression, not on steps.
- **New, step-level or process-level:**
  - PRISM-Physics: DAG-of-formulas step credit with rule-based symbolic matching, validated against experts.
  - HiPhO: official olympiad marking schemes, applied by an LLM judge.
  - PhysicsArena: scores variable identification, process formulation and solution derivation separately.
  - PhysElite: process-level failure localisation.
  - ScienceArena: an LLM judge with process-credit rubrics.
  - FrontierScience: its research track is graded with rubrics.
  - CritPt: 190 checkpoints.
  - SciCode: subproblem unit tests.
  - SymPyBench: partial accuracy over subproblems.
  - Sci-Rho: LLM-judged step-level F1 against gold steps generated with each instance.
  - PhysicsEval: an LLM rubric that includes logical consistency.
  - In `process_supervision.json`: PhysPRM with the human-verified PhysProcessBench, and ChemCoTBench-V2 with deterministic intermediate chemical states.
- **Expression-level metrics:** PHYBench's EED and CMPhysBench's SEED.
- **Still final-answer or multiple-choice only:** UGPhysics, PHYSICS, OlympiadBench (it has step annotations but does not score them), SeePhys, PhysUniBench, PhyX, JEEBench, GPQA, MMLU-Pro, HLE, SciDA, AtmosSci-Bench, MatSciBench and ChemBench. TPBench uses auto-verification plus holistic grading.

### Which use generated instances?

- **Cited set:**
  - ABench-Physics Phy_B: 100 problems with an automatic variation engine.
  - NEWTON: templated QA from 2,800 object-attribute pairs.
  - GSM-Symbolic (mathematics).
- **New:**
  - SymPyBench: 15,045 problems, each with a parameterised executable-Python solver.
  - SciDA: parameters re-sampled every evaluation round.
  - AtmosSci-Bench: 67 templates with expert-verified Python solvers, instantiated 10 or 30 times.
  - Sci-Rho: Python templates that emit instances together with gold step-by-step solutions.
  - SciR: generated from formal structures.
  - Infinite Problem Generator: Formula-as-Code generation, for training data.
- **Rejected for scope:** PRiSM (multimodal, agent-generated and Python-verified) and PhysGym (simulated environments).
- **Contamination handled by recency or novelty instead:** HiPhO and ScienceArena use recent exams; CritPt and TPBench use unpublished or novel problems; PRL-Bench uses recent PRL papers.

### Does "these benchmarks remain confined to theoretical physics, lacking the integrative context of engineering" hold up?

Not as written.

1. **The papers the sentence summarises are not theoretical physics.**
   - UGPhysics, PhysReason, PHYBench and ABench-Physics are undergraduate and olympiad problem solving in mechanics, electromagnetism, thermodynamics and optics.
   - NEWTON is object-attribute physical commonsense.
   - LLM-SRBench is equation discovery across physics, chemistry, biology and materials science.
   - "Theoretical physics" fits only the research-level strand: TPBench, CritPt, CMT-Benchmark and PRL-Bench.
2. **The engineering half is still defensible.** None of the physics benchmarks found sets problems in engineering contexts (handbook data, design constraints, integration across engineering branches). Neighbouring benchmarks do edge towards engineering:
   - JEEBench: pre-engineering exams.
   - MatSciBench: materials.
   - AtmosSci-Bench: hydrology and geophysics.
   - ChemBench: chemistry.
   - MMLU-Pro: has an engineering category.
   - SimulCost: simulator tuning for fluids, solids and plasmas.
3. **The May paper's broader claim, that existing benchmarks are "limited to outcome matching", is contradicted within physics.** PhysReason itself (a cited paper) scores steps, as do PRISM-Physics, HiPhO, PhysicsArena, PhysElite, ScienceArena and Sci-Rho.
4. **Generation with code-computed answers is not unique to EngTrace in physical science.** SymPyBench, SciDA, AtmosSci-Bench, Sci-Rho and ABench-Physics all do it. Contamination resistance through generation cannot carry the novelty claim alone.

**What stays distinctive.** This is what the abstracts and the HTML pages I read support; the panel should confirm it against full texts. The distinctive part is the combination of five things:

- (a) engineering problems across five branches and 15 domains;
- (b) a gold trace produced by the same code as the answer, so intermediate values can be checked on every fresh instance;
- (c) deterministic milestone coverage and per-claim arithmetic checks inside free-form traces, with an out-of-roster LLM judge only for the residual;
- (d) validation of the step evaluation against 15 experts' labels on 300 traces;
- (e) paraphrase robustness and template-level bootstrap intervals.

Sci-Rho (school STEM, multimodal, LLM-judged step F1) and SymPyBench (university physics, subproblem partial credit) each share parts of (b) and (c). Neither is engineering.

**Dates, for priority.** The arXiv API returned EngTrace's own preprint, arXiv 2511.01650 (first posted 2025-11-03, "33 pages ... introduces the EngTrace benchmark"). It predates SymPyBench (2025-12-05) and Sci-Rho (2026-06-06). If the authors want to frame those two as concurrent or later work, the ARR anonymity rules apply.

**Draft replacement for the sentence.** Hedged, for the audit to check; the citation labels follow the May paper (Zhang et al. 2025a is PhysReason and 2025b is ABench-Physics):

> Physics benchmarks now go beyond final answers: some score intermediate steps against reference solutions or official marking schemes (Zhang et al., 2025a; Zhao et al., 2026; Yu et al., 2026), and some re-sample problem parameters to resist memorisation (Zhang et al., 2025b; Zhou et al., 2025; Imani et al., 2025). Expert audits, however, attribute most of the errors these benchmarks report to reference solutions and automatic graders rather than to the models (Ansari et al., 2026), and find correct final answers reached by invalid methods (Ren et al., 2026). Their problems come from textbooks, olympiads or physics research; none we know of pairs freshly generated engineering problems with code-derived gold intermediate values checked inside the model's own reasoning.

## 2. Queries

### 2.1 WebSearch (36 run, verbatim)

1. `SeePhys benchmark vision physics reasoning arXiv 2025`
2. `PhysicsEval benchmark LLM physics problems inference-time techniques`
3. `PHYSICS: Benchmarking Foundation Models on University-Level Physics Problem Solving Feng 2025`
4. `TPBench theoretical physics benchmark AI reasoning`
5. `PhysUniBench undergraduate-level physics multimodal benchmark`
6. `PhysicsArena multimodal physics reasoning benchmark variable process solution`
7. `"How Good Are Frontier Models at Physics" expert re-grading broken evaluations near-saturation`
8. `PRL-Bench frontier physics research benchmark LLM 2026`
9. `PhysElite Olympiad-level physics problems LLM benchmark`
10. `CritPt frontier physics research benchmark critical point AI reasoning`
11. `physics benchmark step-level evaluation intermediate steps LLM reasoning 2026`
12. `physics problems generated parameterized templates contamination-free benchmark LLM dynamic instances`
13. `OlympiadBench bilingual multimodal scientific problems ACL 2024 physics`
14. `CURIE scientific long-context understanding reasoning information extraction benchmark ICLR 2025`
15. `Humanity's Last Exam Nature 2025 publication benchmark`
16. `SciCode research coding benchmark curated by scientists NeurIPS 2024`
17. `ChemBench Nature Chemistry 2025 chemical knowledge reasoning large language models chemists`
18. `MatSciBench materials science reasoning benchmark LLM college-level`
19. `HiPhO high school physics olympiad benchmark LLM human-aligned evaluation medal`
20. `FrontierScience benchmark OpenAI expert-level scientific reasoning olympiad research physics chemistry biology`
21. `ChemIQ benchmark molecular comprehension chemical reasoning LLM`
22. `JEEBench "Have LLMs advanced enough" challenging problem solving benchmark EMNLP 2023`
23. `ScienceAgentBench rigorous assessment language agents data-driven scientific discovery ICLR 2025`
24. `multiphysics simulation reasoning LLM benchmark COMSOL OpenFOAM CFD 2025 2026`
25. `GPQA graduate-level Google-proof Q&A benchmark COLM 2024 Rein`
26. `MMLU-Pro more robust challenging multi-task language understanding NeurIPS 2024 datasets benchmarks`
27. `TheoremQA theorem-driven question answering dataset EMNLP 2023 physics EE`
28. `physics process reward model step verification physics reasoning LLM 2025 PRM`
29. `PhyX physical reasoning multimodal benchmark "Does your model have the wits"`
30. `NeurIPS 2025 datasets and benchmarks track physics reasoning benchmark LLM accepted`
31. `ICLR 2026 physics reasoning benchmark large language models`
32. `ACL 2026 physics benchmark language models findings main conference`
33. `EMNLP 2025 physics problem solving benchmark LLM evaluation`
34. `SciReas scientific reasoning suite knowledge versus reasoning probing LLMs 2025`
35. `ChemEval comprehensive multi-level chemical evaluation large language models`
36. `MaScQA materials science question answering benchmark LLM Zaki Digital Discovery`

Six more were sent but refused: the session's WebSearch budget ("200 of 200 WebSearch calls") was exhausted, the shared budget presumably spent mostly by the parallel agents. Their targets were checked instead through arXiv pages (sections 2.3 and 2.4).

37. `SciDA scientific dynamic assessor randomized numerical initialization olympiad problems`
38. `IsoSci isomorphic cross-domain science problems reasoning versus knowledge retrieval`
39. `ScienceArena benchmarking LLMs latest scientific olympiad competitions 2026`
40. `"Chain-of-Thought Physics Benchmark for Reasoning Models" openreview`
41. `PRISM-Physics causal DAG process evaluation physics reasoning ICLR 2026`
42. `HLE-Verified systematic verification structured revision Humanity's Last Exam errors`

### 2.2 arXiv API (export.arxiv.org/api/query)

Every query used `sortBy=submittedDate&sortOrder=descending&max_results=100` and kept entries first posted in 2024 to 2026. The script was `arxiv_query.py` in the scratchpad (curl, then `xml.etree`).

The first run, at 3.5 s spacing, completed 5 of the 14 queries. The other nine returned HTTP 429 (rate limit) and once 503. A retry job with 20 s spacing and backoff from 30 s to 150 s was still being refused on its first query after about 10 minutes, so I stopped it. The IP is shared with the other literature agents. Results are kept in the scratchpad (`arxiv_results_firstrun.json`). Because the queries sorted by date with a cap of 100, the three that hit the cap covered only their most recent 100 matches.

| Query (search_query, verbatim) | Result |
|---|---|
| `all:physics AND all:benchmark AND all:"language models"` | HTTP 429, not retrieved |
| `all:"physical reasoning" AND all:benchmark` | HTTP 429 |
| `all:"scientific reasoning" AND all:benchmark` | HTTP 429 |
| `all:"step-level" AND all:physics` | HTTP 429 |
| `all:symbolic AND all:physics AND all:problems AND all:LLM` | 51 total, 48 kept |
| `all:olympiad AND all:benchmark AND (all:physics OR all:science)` | 142 total, 100 kept |
| `ti:physics AND (ti:benchmark OR ti:benchmarking OR ti:evaluating)` | 439 total, 100 kept |
| `ti:physics AND ti:reasoning` | 261 total, 100 kept |
| `all:physics AND (all:"process evaluation" OR all:"intermediate steps" OR all:"process-level")` | HTTP 429 |
| `all:benchmark AND all:contamination AND (all:physics OR all:scientific) AND (all:dynamic OR all:generated OR all:templates)` | 64 total, 55 kept |
| `ti:chemistry AND (ti:benchmark OR ti:benchmarking) AND all:"language models"` | HTTP 429 |
| `all:"materials science" AND ti:benchmark AND all:"language models"` | HTTP 429 |
| `all:multiphysics AND all:"language models"` | HTTP 429 |
| `all:graduate AND all:science AND ti:benchmark AND all:LLM` | HTTP 429 |
| `abs:physics AND abs:LLMs AND abs:benchmark`, `abs:scientific AND abs:benchmark AND abs:step AND abs:LLMs`, `abs:physics AND abs:contamination AND abs:LLMs` | planned, never sent (job stopped) |

### 2.3 arXiv listing searches (arxiv.org/search via WebFetch, newest first)

These stood in for the refused API variants. WebFetch returns a summarised listing, and long listings were truncated: the first title search printed 107 of its 200 rows, and I fetched the rest with `start=100`.

- `query=physics benchmark "language models"`, all fields, size 100
- `query=physics benchmark`, title, size 200, then `start=100`
- `query=physics reasoning`, title, size 200
- `query=scientific reasoning benchmark`, title, size 200
- `query=physics problems LLM`, abstract, size 200
- `query=chemistry benchmark language models`, title, size 100 (4 hits)
- `query=multiphysics LLM`, all fields, size 100 (5 hits)
- `query=step-level physics`, abstract, size 100 (mostly off-topic)
- `query=CMT-Benchmark`, all fields, size 50 (found 2510.05228)

### 2.4 Pages opened for verification

- **arXiv abstract pages (WebFetch):**
  - Included candidates: 2510.03185, 2509.07894, 2512.05954, 2606.08034, 2506.12909, 2502.01159, 2609.13009, 2608.02442, 2608.25097, 2505.15472, 2509.26574, 2502.15815, 2503.21821, 2608.30517, 2601.21165, 2501.14249, 2311.12022, 2406.01574, 2402.14008, 2305.15074, 2407.13168, 2404.01475, 2510.12171, 2505.19099, 2602.13964, 2606.13020, 2508.18124, 2603.14486, 2305.12524, 2508.00079.
  - Rejected candidates: 2506.17667, 2604.15411, 2607.01431, 2508.19202, 2503.13517, 2505.07735, 2410.05080, 2507.15550, 2505.15929, 2505.24182, 2505.24823, 2607.05199, 2609.13152, 2608.22126, 2605.01203, 2607.13190, 2607.00276, 2608.17724, 2511.14366, 2605.14040, 2603.20253, 2308.09115, 2512.05930.
  - Cited-set updates: 2502.00334, 2502.12054, 2504.16074, 2507.04766, 2310.07018, 2502.14739.
- **arXiv HTML full text (WebFetch, summarised):** 2512.05954, 2506.12909, 2510.03185, 2608.02442, 2609.13009, 2509.07894, 2502.01159, 2606.08034, 2508.00079v2.
- **Other pages:**
  - The ACL Anthology page for PhysPRM (2026.findings-acl.458).
  - The SeePhys GitHub README, for its venue, size and judge.
  - Crossref API records (curl) for 10.1038/s41586-025-09962-4 (HLE in Nature), 10.1145/3770855.3818888 (MatSciBench at KDD 2026), 10.1088/2632-2153/adfcb0 (TPBench in MLST) and 10.1038/s41557-025-01815-x (ChemBench in Nature Chemistry).
  - nature.com redirected to a login page and dl.acm.org returned 403, which is why the Crossref records were used.

## 3. Included candidates

### Must-cite (8)

**zhao2026prismphysics: PRISM-Physics (Zhao et al., ICLR 2026; arXiv 2510.03185, Oct 2025).**
- Problems come from *Major American Universities Ph.D. Qualifying Questions and Solutions*, in seven areas and three difficulty levels. The items are static.
  - Ansari et al. (2026, appendix B.1.2) count 1,401 questions, 549 of them image-dependent; the PRISM-Physics pages I read do not state a total.
- Each reference solution is rewritten into a DAG of formulas by an LLM-assisted pipeline with rule-based and LLM checks. A response earns ancestor-closure credit for the formulas it matches, and matching is rule-based symbolic and numeric equivalence (constant substitution, unit conversion), not an LLM judge.
- Validation: two experts (an IPhO gold medallist and a physics PhD), with a third adjudicating, scored 70 problems. Kendall tau-b with the expert scores is 0.346 for PRISM-DAG, 0.294 for an LLM judge and 0.213 for PhysReason's PSAS-S.
- Bearing: the closest physics analogue to EngTrace's step evaluation. EngTrace's intermediate values come from the generating code on fresh instances, not from rewritten textbook solutions.
- Ansari et al. (2026) re-graded 74 retained PRISM-Physics items: grader and reference errors dominated, and the corrected mean@4 rose from 13.0% to 94.6%.

**yu2026hipho: HiPhO (Yu et al., ICML 2026; arXiv 2509.07894, Sep 2025).**
- 13 olympiad exams from 2024-2025 (IPhO, APhO, EuPhO and regional contests): 360 problems and 519 subquestions, text and diagram items. Static; recency is the contamination control.
- Grading applies the official marking schemes at answer and step level, with Gemini-2.5-Flash as judge. The paper reports judge-versus-human differences under one point on IPhO 2024 examples.
- Models are placed against human medal thresholds.
- Bearing: the venue-backed example of step-level physics grading. Its step credit depends on an LLM applying rubrics; EngTrace checks milestones deterministically and sends only the residual to a judge.

**imani2025sympybench: SymPyBench (Imani et al., arXiv 2512.05954, Dec 2025).**
- 15,045 university-physics problems (mechanics about 34%, electricity and magnetism about 27%, and other areas), built from open Creative-Commons problem sets.
- An LLM writes a parameterised Python solver for each problem. A solver is kept only if it reproduces the original answers, which held for about 88% of problems. The paper body says all problems were manually reviewed.
- Instances can be generated for any parameter setting. Question types are symbolic multiple-choice, numerical multiple-choice and free-form.
- Scoring: exact match, partial accuracy over the subproblems of a structured solution, and consistency, failure and confusion rates across variants. Nine models are evaluated.
- Bearing: the nearest prior art to EngTrace's generator and gold traces in a physical domain. It was posted a month after EngTrace's preprint, but the revision must still state the differences:
  - engineering contexts across five branches;
  - milestone matching and per-claim arithmetic checks inside free-form traces;
  - expert validation of the step evaluation.
  - Confirm SymPyBench's own step scoring against its text.
- Also in `symbolic_contamination.json`, under the same key.

**azmi2026scirho: Sci-Rho (Azmi et al., arXiv 2606.08034, Jun 2026; v2 Sep 2026).**
- A multilingual (seven languages), visually grounded STEM benchmark: 606 Python templates per language (4,242 in all; 42,420 instances; at least ten per template).
- Subjects: mathematics, physics, chemistry, biology and computer science at GCSE and A-level. 195 of the 606 base templates are physics (text p. 5).
- The templates were written by STEM graduates and olympiad medallists. A QC team then checked the generated steps and answers.
- Each instance comes with a dynamically generated step-by-step LaTeX solution.
- Evaluation: average and worst-case accuracy across variants, plus a step-level F1 in which an LLM judge matches model steps to gold steps. 17 VLMs are evaluated.
- Bearing: it already combines generated instances with code-generated gold steps and step-level scoring, so EngTrace cannot claim to be the first symbolic benchmark with gold reasoning traces. EngTrace differs in being:
  - text-only, university-level engineering;
  - checked by deterministic milestone and arithmetic checks validated against experts;
  - judged by an LLM only on the residual.
- Also in `symbolic_contamination.json`, under the same key.

**zhou2025scida: SciDA (Zhou et al., arXiv 2506.12909, Jun 2025).**
- Over 1,000 olympiad-level numerical problems in mathematics, physics, chemistry and biology.
- Sources: IMO, CMO, IPhO, CPhO, IChO, CChO and IBO problems, workbooks and textbooks. Over 20% were contributed privately by medallists, coaches and professors.
- Every modifiable parameter becomes a variable, re-sampled uniformly within its range on each evaluation round. The final number is scored within a tolerance; steps are not evaluated.
- Accuracy drops 10-20 points (20-60% relative) under random initialisation, which the authors read as memorised numeric patterns. They also hand-analysed the errors of 50 sampled problems per subject.
- Bearing: the most direct external evidence that re-sampled instances change measured performance. EngTrace adds gold traces and step checks.

**li2025atmosscibench: AtmosSci-Bench (Li et al., arXiv 2502.01159, Feb 2025; v3 Oct 2025).**
- Atmospheric science: hydrology, atmospheric dynamics, atmospheric physics, geophysics and physical oceanography, built from university course materials.
- Experts turned questions into 67 templates with variable ranges and constraints. GPT-4o drafted a Python solver for each, and experts verified and refined it.
- Templates are instantiated 10 times (MCQ10, 670 items) or 30 times (MCQ30, 2,010 items). There are also 391 open-ended questions.
- Scoring: multiple-choice by option match; open-ended answers by a cascade of a 5% numeric check, SymPy equivalence and a GPT-4o-mini rubric judge. Steps are not scored.
- Bearing: close to EngTrace's civil water-resources domain, it already does template generation with verified code solvers and a deterministic-then-LLM cascade. EngTrace's novelty must rest on gold traces and step checks, not generation.

**ansari2026physicsregrading: How Good Are Frontier Models at Physics? (Ansari et al., arXiv 2609.13009, Sep 2026; 51 authors).**
- Physicists (faculty and their graduate researchers) re-graded frontier-model answers marked wrong on six benchmarks: HLE-Physics, CMT-Benchmark, CritPt, UGPhysics, PRISM-Physics and PHYBench. The audit kept to text-only problems with deterministic solutions.
- The pooled audit covered 250 rejections across four benchmarks: 143 benchmark defects, 95 grader errors and 12 genuine model errors. CMT-Benchmark and CritPt were audited separately.
- Corrected scores rise sharply. On PHYBench (87 retained items) the mean@4 went from 26.5% to 90.2%, mostly because EED rejects equivalent forms. On UGPhysics (82 items) it went from 83.0% to 92.1%.
- Conclusion: closed-ended physics benchmarks understate frontier models and are near saturation.
- Bearing: it audits two cited benchmarks. It supports EngTrace's code-derived gold and expert-certified templates, and warns that answer keys and graders, not models, can dominate measured error.
- Also in `symbolic_contamination.json`, under the same key.

**ren2026solutionhacking: Right Answer, Wrong Method (Ren et al., arXiv 2608.02442, Aug 2026; comment "working in progress").**
- "Solution hacking" means a correct final answer reached by numerical search, enumeration, guessing or answer-first verification.
- Benchmark tiers: SciBench and MATH-500 (easy); IMO-Bench, PHYBench and SciOlympiad (medium); HLE (hard). Eleven models.
- 21 PhD experts (8 mathematicians, 7 physicists, 6 chemists) labelled 300 solutions. A majority vote of three LLMs then scaled up the labelling. The panel agrees with the experts 75.4% of the time, with 17% false positives and 39% false negatives, so its hack rates are lower bounds.
- Hack rates are 2.2% on the easy tier, 28.3% on the medium tier and 37.4% on HLE. Across models, 8.2%-44.1% of answers credited as correct were hacked.
- Bearing: direct evidence for EngTrace's premise that final-answer accuracy overstates reasoning. The panel's physics recall (0.58) is also a caution for any LLM judge, EngTrace's residual judge included.

### Should-cite (15)

**xu2026physelite: PhysElite (Xu et al., NeurIPS 2026; arXiv 2608.25097, Aug 2026).**
- A bilingual multimodal benchmark of 11,586 olympiad-tier physics problems, each with diagrams, bilingual step-by-step derivations and a verified final answer. Static.
- 18 models are scored on answers (best 33.7%), and a process-level evaluation locates where reasoning fails.
- The abstract does not detail contamination or human validation.
- Bearing: the newest and largest physics set, and another reason the revision cannot say physics benchmarks only score outcomes.

**dai2025physicsarena: PhysicsArena (Dai et al., Findings of EMNLP 2025; arXiv 2505.15472, May 2025).**
- Multimodal high-school physics problems, collected and AI-assisted-annotated. Static.
- Scores three stages: variable identification, physical-process formulation and solution derivation. Best derivation accuracy is 33.47%, and the first two stages predict the third.
- The abstract page did not give the size or per-dimension scoring.
- Bearing: explicit intermediate-stage scoring, coarser than EngTrace's milestones.

**zhu2025critpt: CritPt (Zhu et al., arXiv 2509.26574, Sep 2025; v4 May 2026; 65 authors).**
- 71 composite research challenges across modern physics, written by 50+ active researchers from their own work, and split into 190 checkpoints.
- Answers are unpublished, guess-resistant and machine-verifiable, graded by a physics-specific automatic pipeline.
- The best base model scores 5.7% on full challenges, about 10% with code tools. Static.
- Bearing: the research-level strand, where "theoretical physics" is accurate. Ansari et al. report much higher corrected scores after expert re-grading.

**chung2025tpbench: TPBench (Chung et al., Machine Learning: Science and Technology 6(3) 030505, 2025; arXiv 2502.15815, Feb 2025).**
- 57 novel high-energy-theory and cosmology problems across five difficulty levels, from undergraduate to research. They are deliberately absent from public collections, to limit leakage.
- Graded by auto-verification plus holistic grading. Static.
- Bearing: genuinely theoretical physics. Cite it if the revised sentence keeps that word, and only for this strand.

**feng2025physics: PHYSICS (Feng et al., Findings of ACL 2025; arXiv 2503.21821, Mar 2025).**
- 1,297 expert-annotated university-level problems in six areas: classical mechanics, quantum mechanics, thermodynamics and statistical mechanics, electromagnetism, atomic physics and optics.
- Automated answer evaluation (o3-mini best at 59.9%), plus error analysis, prompting and RAG studies. Static.
- Bearing: an ACL-venue, outcome-scored physics benchmark to list beside UGPhysics.

**zhao2026sciencearena: ScienceArena (Zhao et al., EMNLP 2026 main per arXiv comment; arXiv 2608.30517, Aug 2026).**
- Digitised official exams, figures, solutions and rubrics, verified by olympiad medallists. Sources: IPhO and IChO 2025-2026, IBO 2023, USAPhO 2026, USNCO 2025 and seven other public competitions.
- Static, with recency as the contamination control.
- An LLM judge calibrated against medallist ground truth applies process-credit rubrics. 14 models.
- Bearing: a current example of rubric-based process credit by an LLM judge, to contrast with EngTrace's deterministic milestones plus residual judge.

**wang2026frontierscience: FrontierScience (Wang et al., OpenAI; arXiv 2601.21165, Jan 2026).**
- Physics, chemistry and biology in two tracks:
  - an olympiad track by international medallists and coaches (short answers);
  - a research track of PhD-level open-ended sub-problems, written and verified by PhD scientists and graded with granular rubrics.
- Several hundred questions, 160 of them open-sourced as a gold set. Original, static items.
- Bearing: the frontier-lab reference for expert-level science. Its research rubrics are process evaluation by rubric, not by computed milestones.

**phan2026hle: Humanity's Last Exam (Phan et al.; Nature 649:1139-1146, 2026; arXiv 2501.14249, Jan 2025, v11 Jul 2026).**
- 2,500 questions across dozens of subjects, written by roughly a thousand experts and filtered against frontier models.
- Multiple-choice and exact-match short answers, graded automatically. Static.
- Its answer key has since been audited by HLE-Verified, by Ansari et al. (HLE-Physics) and by Ren et al.
- Nature published it under a different title ("A benchmark of expert-level academic questions to assess AI capabilities"), with the Center for AI Safety consortium as author.
- Bearing: the general frontier suite to cite with GPQA and the cited SuperGPQA.

**rein2024gpqa: GPQA (Rein et al., COLM 2024; arXiv 2311.12022, Nov 2023).**
- 448 expert-written multiple-choice questions in biology, physics and chemistry.
- Experts reach 65% (74% discounting clear mistakes). Skilled non-experts reach 34% even with web access.
- The 198-item Diamond subset is the usual reporting set. Static.
- Bearing: the canonical graduate-science suite, and the obvious companion to the cited SuperGPQA.

**wang2024mmlupro: MMLU-Pro (Wang et al., NeurIPS 2024 Datasets and Benchmarks Track, Spotlight; arXiv 2406.01574, Jun 2024).**
- 12K+ questions over 14 domains, including physics, chemistry and engineering. There are ten options instead of four, and trivial or noisy MMLU items were removed.
- Accuracy is 16-33 points lower than on MMLU. Prompt sensitivity falls from 4-5% to 2%, and chain-of-thought helps here, unlike on MMLU. Static.
- Bearing: multiple-choice suites have moved from recall towards reasoning, so the May label "factual recall" for generalist suites should be softened.

**he2024olympiadbench: OlympiadBench (He et al., ACL 2024; arXiv 2402.14008, Feb 2024).**
- 8,476 olympiad and Gaokao mathematics and physics problems, bilingual and partly multimodal, each with expert step-by-step solutions.
- Scored automatically on final answers: GPT-4V reaches 17.97% overall and 10.74% on physics. Static.
- Bearing: the reference olympiad-physics set. Its step annotations are not scored, which is the gap EngTrace's milestones fill.

**arora2023jeebench: JEEBench (Arora, Singh and Mausam, EMNLP 2023; arXiv 2305.15074, May 2023).**
- 515 pre-engineering mathematics, physics and chemistry problems from IIT JEE-Advanced: single- and multi-correct multiple-choice, integer and numeric answers. Static, scored on final answers.
- The best result was under 40% at the time.
- The error analysis names three failure types: algebra, translating concepts into equations, and retrieving domain concepts.
- Bearing: the exam benchmark nearest to engineering in level and content. It tests the physics EngTrace builds on, but not engineering contexts.

**tian2024scicode: SciCode (Tian et al., NeurIPS 2024 Datasets and Benchmarks Track; arXiv 2407.13168, Jul 2024).**
- 80 research problems from 16 natural-science subfields, decomposed into 338 subproblems. Each has scientist-written gold solutions and test cases.
- Models write code, which is scored by tests. The best model solves only a few percent of full problems in the most realistic setting. Static.
- Bearing: a different route to intermediate checking (decomposed subproblems with executable tests). EngTrace checks intermediate values inside a natural-language trace.

**mirza2025chembench: ChemBench (Mirza et al., Nature Chemistry 17:1027-1034, 2025; arXiv 2404.01475, Apr 2024, as "Are large language models superhuman chemists?").**
- 2,788 chemistry question-answer pairs (1,039 written by hand, 1,749 generated semi-automatically), tagged by topic, skill and difficulty. Multiple-choice and numeric formats.
- The best models beat the human chemists surveyed on average, but fail some basics and are overconfident. Static.
- Bearing: the reference chemistry benchmark beside EngTrace's chemical-engineering branch. It is outcome-scored, with a human-expert baseline.

**zhang2026matscibench: MatSciBench (Zhang et al., KDD 2026 per Crossref; arXiv 2510.12171, Oct 2025).**
- 1,340 college-level materials-science problems from ten textbooks, in 6 fields and 31 subfields, with three difficulty tiers by solution length. 946 have reference solutions and 315 include images.
- Scoring: rule-based checks with a 5% numeric tolerance, plus an LLM judge for formula answers. DeepSeek-R1 is best on text (75.22%) and GPT-5 on images (53.02%). Static.
- Bearing: the closest physical-science analogue to EngTrace's deterministic check followed by a judge, next to its mechanics-of-materials and geotechnical domains.

### Optional (7)

**xiang2025seephys: SeePhys (Xiang et al., NeurIPS 2025 per the authors' GitHub; arXiv 2505.19099, May 2025).**
- 2,000 physics questions, from middle school to PhD qualifying exams, in seven domains with 21 diagram types; 75% are vision-essential.
- Scored by an LLM judge, per the repository's evaluation script. Best models score under 60%. Static.
- Bearing: one citation to mark multimodal physics (with PhysUniBench and PhyX) as outside EngTrace's text-only scope.

**zhai2026hleverified: HLE-Verified (Zhai et al., arXiv 2602.13964, Feb 2026; v4 Aug 2026).**
- Experts and models cross-checked all of HLE in two stages: 668 items passed unchanged, and 1,143 were repaired by dual independent expert revisions with adjudication.
- Models gain 7-10 points on average, and 30-40 points on items whose statements or answers were flawed.
- Bearing: a second answer-key-validity audit, needed only if HLE is cited.

**beckmann2026scir: SciR (Beckmann, Valentino and Freitas, NeurIPS 2026 Evaluations and Datasets Track; arXiv 2606.13020, Jun 2026; v3 Sep 2026).**
- Generates deduction, induction and causal-abduction tasks from formal objects (deduction trees, rule hypotheses, causal graphs), rendered as multi-document scientific prose with verifiable answers.
- Two difficulty axes can be varied independently: information extraction and inference. Six models.
- Bearing: generated instances with controlled difficulty in science, but logical rather than quantitative reasoning.

**wang2025cmphysbench: CMPhysBench (Wang et al., arXiv 2508.18124, Aug 2025).**
- 520+ graduate-level condensed-matter calculation problems. Static.
- Introduces SEED, a tree-based expression edit distance that gives partial credit. Grok-4 is best, at 36 SEED and 28% accuracy.
- Bearing: the expression-metric line (EED, SEED) that EngTrace does not follow. Ansari et al. trace most PHYBench grader errors to EED's rigidity with equivalent forms.

**sharan2026infiniteproblemgenerator: Infinite Problem Generator (Sharan, Hebbale and Kumar, arXiv 2603.14486, Mar 2026).**
- An agentic Formula-as-Code pipeline writes each physics solution as executable Python, which guarantees solvability.
- It releases ClassicalMechanicsV1: 1,335 problems from 165 expert seeds over 102 formulas.
- Verification-code length tracks formula count (R^2 about 0.95).
- Bearing: code-verified physics generation, but for training data, not evaluation.

**chen2023theoremqa: TheoremQA (Chen et al., EMNLP 2023; arXiv 2305.12524, May 2023).**
- 800 expert-curated questions applying 350+ theorems in mathematics, physics, EE&CS and finance, scored on final answers.
- GPT-4 reaches 51% with program-of-thought. Static.
- Bearing: an early university-level science benchmark that already includes electrical-engineering theorems.

**siddique2025physicseval: PhysicsEval (Siddique et al., Findings of IJCNLP-AACL 2025; arXiv 2508.00079, Jul 2025).**
- 19,609 problems from 20 textbooks. Their forum and educational-site solutions were restructured into steps by Gemini 2.5 Pro; only a small sample was manually reviewed.
- Scored by a Gemini 2.5 Pro rubric with six sub-metrics (one is logical consistency), aggregated to 0-100. The paper also studies self-correction and multi-agent verification. Static.
- Bearing: a contrast case of LLM-restructured gold and LLM-rubric step scoring, against EngTrace's code-derived gold and deterministic checks.

## 4. Considered and rejected (one line each)

**Physics: text and problem solving**

- PhysUniBench (2506.17667): 3,304 multimodal undergraduate questions, answer-scored; SeePhys represents the multimodal line.
- PhyX (2505.15929): 3,000 visual-scenario physics problems, answer-scored, multimodal; same reason.
- Multi-Physics (2509.15839), OmniPhys (2608.25398, Findings of EMNLP 2026), SeePhys Pro (2605.09266), FeynmanBench (2604.03893), PhysAlign (2609.33319), mmJEE-Eval (2511.09339): multimodal or diagram benchmarks, or RLVR diagnostics; outside text-only scope.
- PRL-Bench (2604.15411): 100 research tasks from 2025-2026 PRL papers; research-level theory is already covered by CritPt and TPBench.
- CMT-Benchmark (2510.05228): expert-built condensed-matter theory, LLM-judged; covered through Ansari et al.'s audit.
- PhySense (2505.24823): principle-first reasoning benchmark, size not on the abstract page; peripheral.
- MVPBench (2505.24182): graph-based consistency for visual physical CoT; multimodal.
- Physics-R1 (2605.14040): single-author audit of physics evaluation infrastructure plus released verifiers. It is relevant to judge dependence, but Ansari et al. make the same case with expert evidence.
- PhysicsMinions (2509.24855), P1 (2511.13612), P1-VL (2602.09443), "Solving Physics Olympiad via RL on Physics Simulators" (2604.11805), "Mastering Olympiad-Level Physics with AI" (2511.10515, not opened): models and agents, not benchmarks.
- "Benchmarking Foundation Models with RAG in Olympic-Level Physics Problem Solving" (2510.00919, Findings of EMNLP 2025): a RAG study.
- PhysProver (2601.15737), Lean4Physics (2510.26094), Ax-Prover (2510.12787), AxQM (2609.05157): formal theorem proving, a different task.
- "Test-time Scaling Techniques ... on the TPBench Dataset" (2506.20729): a method paper.
- "Reason, Reward, Refine" (2607.05199), "Decoupled Physical Modeling and Execution" (2608.22126), "Scientific Logicality Enriched Methodology" (2605.17104), QuantumQA (2604.18176): training or method papers with no new benchmark.
- PhysMent (2609.13152; arXiv shows v1 2026-07-07): 105 interactive classical-mechanics scenes, an interactive-experimentation task.
- "Testing Frontier LLMs' Physics Literacy in Parallel Physical Worlds" (2607.00276): single author, three models, counterfactual laws; too small.
- "Reliable isomorphic physics problem generation" (2607.13190): an education tool tested on 13 multiple-choice questions.
- "Agentic Retrieval and RL Equation Chains ... Physics Word Problems" (2606.15591): a generation framework, not opened.
- IsoSci (2607.01431): isomorphic cross-domain pairs testing reasoning against knowledge retrieval; the abstract does not say how the pairs are built. Possible symbolic_contamination material.
- "Grading the Unspoken" (2604.14188, QFT and string theory) and "Assessing AI in Introductory Physics Problem Solving" (2607.14303): not opened; peripheral.
- "Chain-of-Thought Physics Benchmark for Reasoning Models" (OpenReview yHHOPIb645; Zenodo 21365628): seen in a search result, could not be verified (no arXiv id found; search budget gone).

**Embodied or visual physical reasoning (different modality and task)**

- PhysBench (2501.16411), PhysToolBench (2510.09507), DeepPHY (2508.05405), QuantiPhy (2512.19526), PhysicsMind (2601.16007), BilliardPhys-Bench (2605.30900), Causal Scaffolding (2606.05966, KDD 2026), MechReasoner (2609.34636), IntPhys 2 (2506.09849), plus video and world-model physics benchmarks.

**General science suites**

- SciReas (2508.19202, ICML 2026): a meta-suite of existing benchmarks with a knowledge-versus-reasoning probe. A defensible optional, left out for the 30-entry cap.
- CURIE (2503.13517, ICLR 2025): long-context extraction from papers; a different task.
- ScienceAgentBench (2410.05080, ICLR 2025) and AstaBench (2510.21652, ICLR 2026): agentic data-driven discovery and code tasks.
- ATLAS (2511.14366): about 800 expert problems judged by an LLM panel; FrontierScience and HLE cover the role.
- General365 (2604.11778), SciVQR (2605.10187, multimodal), Sci-MMR (2609.11243), TRACES (2608.11415), SPUR, PaperArena: off-topic or multimodal.
- PRiSM (2512.05930): the SymPyBench authors' multimodal sibling, 24,750+ agent-generated, Python-verified problems; multimodal. Cite with SymPyBench only if needed.
- LLM Olympiad: Why Model Evaluation Needs a Sealed Exam (2603.23292), and BeyondBench (2509.24210): contamination position papers or algorithmic generators; belong to symbolic_contamination.
- SCI-Verifier (2509.24285) and Sci-PRM (2606.04579): verifiers and PRMs; belong to process_supervision or llm_judge.

**Chemistry and materials**

- ChemIQ (2505.07735; J. Chem. Inf. Model. per the ACS listing): 816 organic-chemistry short-answer molecular-comprehension questions; ChemBench covers chemistry.
- ChemEval (2409.13989): a broad chemistry task suite; less related.
- ChemPro (2602.03108): not opened.
- ChemLLMBench (2305.18365, 2023): superseded.
- Chemistry olympiad exam evaluation (2512.14989, Communications Chemistry): multimodal exam study.
- TSBench (2609.08503): not opened.
- MaScQA (2308.09115; Digital Discovery 2024): 650 GATE-derived materials questions. MatSciBench is newer and broader, but the GATE source (an engineering entrance exam) may interest the engineering agents.
- MSQA (2505.23982): graduate materials; already in engineering_web.json.

**Multiphysics and simulation (left to the engineering catalogues; FEABench is in the cited set)**

- SimulCost (2603.20253; ICML 2026 per the arXiv comment): 2,643 single-round and 2,304 multi-round parameter-tuning tasks over 11 fluid, solid and plasma simulators. The best representative of simulation reasoning, cut at the cap; tool use, not closed-form reasoning.
- CFDLLMBench (2509.20374) and FEM-Bench (2512.20732): already in engineering_web.json.
- "Your Simulation Runs but Solves the Wrong Physics" (2605.09360), "False Summit and Silent Drift" (2606.21841), PHITSBench (2607.09789), GPUPhysBench (2609.35639), FireWorldBench (2609.23064), the Modelica Agent Workflow Benchmark (2608.23653), Deep Research in Physical Sciences (2606.18648): simulation and code agents, not opened in full.
- Engineering items seen in the listings, for the engineering agents: PolyBridgeBench (2609.21493), CADEngBench (2608.09296), EMRB (2608.24086), BuildArena (2510.16559), PCEval (2601.02404). MechReason (2609.16012) is already in engineering_web.json.

**Other**

- PhysGym (2507.15550; NeurIPS 2025 Datasets and Benchmarks per the proceedings URL): interactive law discovery in simulations with controlled priors; a different task.
- Gravity-Bench-v1 (2501.18411), Collider-Bench (2605.13950), HEPToolBench (2608.28232), gwBenchmarks (2605.11269), QC-Stark (2609.35581): agentic or tool tasks in specialised physics.
- VERaiPHY (2608.17724): a statistics-standards initiative, not a benchmark.
- "Physics as the label for measuring and correcting materials reasoning" (2609.12181), BioPhys-Bridge (2609.19180, EMNLP 2026), SFBench (2606.29630): not opened; peripheral.

## 5. Covered in other agents' catalogues (not duplicated here)

- `process_supervision.json` already holds five items I would otherwise have added:
  - dong2026physprm: PhysPRM, Findings of ACL 2026, with the human-verified PhysProcessBench.
  - sun2026lsrben: LSR-Ben, a PRM benchmark in scientific and logical reasoning.
  - guo2026chemcotbenchv2: ChemCoTBench-V2, deterministic intermediate-state checks in chemistry.
  - zi2026neurosymbolicprm.
  - yeadon2026answerargument: HLE physics traces whose correct answers hide errors.

  PhysPRM and ChemCoTBench-V2 also belong in the physics and science paragraph.
- `symbolic_contamination.json` chose the same keys and arXiv ids for imani2025sympybench, azmi2026scirho and ansari2026physicsregrading. The fetch script will list them as duplicate keys when run over all catalogues and keep a single copy.

## 6. Updates to existing citations

- **UGPhysics (xu2025ugphysics):**
  - Accepted to ICML 2025 (arXiv comment; v4 3 Jun 2025). papers.json says only "arXiv".
  - 5,520 problems in 13 subjects, in English and Chinese, judged by the MARJ pipeline (model-assistant rule-based judgement), so it is outcome-scored.
  - Ansari et al. 2026 audited 82 retained items: the pre-audit to corrected mean@4 went from 83.0% to 92.1%, and most rejections were benchmark defects.
- **PhysReason (zhang2025physreason):**
  - ACL 2025 main, long paper (ACL Anthology 2025.acl-long.811, seen as a search-result URL). arXiv v2 is dated 26 May 2025 and has no comment.
  - It scores steps (PSAS-S), so the May framing that existing benchmarks are limited to outcome matching does not apply to it.
  - PRISM-Physics reports PSAS-S agreement with expert scores at Kendall tau-b 0.213.
- **PHYBench (qiu2025phybench):**
  - NeurIPS 2025, Datasets and Benchmarks Track (neurips.cc poster page in a search result). arXiv v2 18 May 2025; 500 problems; EED score; human experts 61.9%.
  - Ansari et al.: 40 of 56 audited PHYBench rejections were grader errors from EED rigidity, and the corrected mean@4 went from 26.5% to 90.2%.
  - Ren et al. use PHYBench in their medium tier.
- **ABench-Physics (zhang2025abenchphysics):**
  - arXiv v1 only (7 Jul 2025); no venue found.
  - Phy_A has 400 static problems; Phy_B has 100 problems with an automatic variation engine. It is the cited set's own generated-instance physics benchmark, and the July rebuttal cites its 1% tolerance.
- **NEWTON (wang2023newton):**
  - Findings of EMNLP 2023 (arXiv journal-ref), which matches papers.json.
  - Its 160K QA items are generated from 2,800 object-attribute pairs. It is templated physical commonsense, not "theoretical physics".
- **LLM-SRBench (shojaee2025llmsrbench):** ICML 2025 Oral (arXiv comment); papers.json says ICML 2025, so add "Oral". It is scientific equation discovery, not theoretical-physics problem solving.
- **SciBench (wang2023scibench):** ICML 2024 (papers.json correct). Ren et al. use it as their easy tier, with a 2.2% hack rate.
- **SuperGPQA (du2025supergpqa):** arXiv v4 28 Mar 2025, with no venue in the arXiv comments; no acceptance verified. The May text's label "factual recall via multiple-choice" is loose for a 285-discipline suite that includes calculation items; the cited-set panel should check.

## 7. What I could not verify

- **Search tooling:**
  - The WebSearch budget ran out after 36 queries; six planned queries were refused (section 2.1).
  - The arXiv API returned HTTP 429 or 503 for 9 of 14 query variants. I used arxiv.org/search listings through WebFetch instead; these are summarised and may be incomplete.
  - No API pass covered materials science specifically; that coverage rests on WebSearch.
- **Venues taken from secondary pages, not opened:**
  - PRISM-Physics, ICLR 2026: proceedings.iclr.cc URL in a search result.
  - HiPhO, ICML 2026: icml.cc poster URL in a search result.
  - PHYBench, NeurIPS 2025 Datasets and Benchmarks: neurips.cc poster URL.
  - PhysGym, NeurIPS 2025 Datasets and Benchmarks: proceedings URL.
  - PhysReason, ACL 2025: anthology URL.
  - SeePhys, NeurIPS 2025: the authors' GitHub README; the arXiv comment says only "46 pages".
  - PhysElite: the "Evaluations and Datasets" track comes from GitHub; arXiv says only NeurIPS 2026.
  - SimulCost, ICML 2026: inferred from the arXiv comment wording.
  - ChemIQ: J. Chem. Inf. Model. 66(1), per an ACS URL.
- **HLE:** the Nature version (Crossref: 649(8099):1139-1146, online 2026-01-28) has a different title, with the consortium as first author. The key uses Long Phan, the arXiv first author, and the year 2026.
- **Already spot-checked against the downloaded texts** (`text/*.txt`, grep):
  - Ansari et al.: 143 benchmark, 95 grader and 12 model errors out of 250. Corrections: PHYBench 26.50% to 90.23% (87 items; 13 benchmark, 40 grader, 3 model errors among 56 rejections); PRISM-Physics 13.00% to 94.59% (74 items); UGPhysics 83.00% to 92.07% (82 items).
  - PRISM-Physics: tau-b values 0.346, 0.294 and 0.213.
  - HiPhO: 360 problems and 519 subquestions.
  - Ren et al.: 8.2%-44.1% hacked, and 75.4% panel-expert agreement.
  - SymPyBench: about 88% code success, manual review of all problems, and Partial Accuracy defined as the fraction of subproblems answered correctly.
  - AtmosSci-Bench: MCQ30, and 391 open-ended questions.
  - SciDA: uniform sampling; drops of 10%-20% absolute and 20%-60% relative.
  - Sci-Rho: 606 templates, 4,242 counting languages, and 42,420 instances.
  - FrontierScience: "several hundred questions (160 in the open-sourced gold set)".
- **Numbers from summarised full texts.** Everything else came through WebFetch summaries of the arXiv HTML and still needs page references from the panel:
  - SymPyBench: subproblem partial accuracy, the 88% code survival and the manual review.
  - SciDA: the 10-20 point drop.
  - PRISM-Physics: the tau-b values and the 70-problem study.
  - HiPhO: under one point of judge-versus-human difference.
  - Ansari et al.: the 143/95/12 breakdown and per-benchmark corrections.
  - Ren et al.: rates and judge errors.
  - AtmosSci-Bench: counts and cascade.
  - Sci-Rho: template counts and step-F1 judge.
  - PhysicsEval: rubric weights.
- **Counts I could not pin down:**
  - FrontierScience's total question count: the paper says "several hundred"; the "700+" from OpenAI's announcement is unverified.
  - PhysicsArena's size: not on the abstract page.
- **Contradictory source:** PhysicsEval's abstract says its solutions were scraped, while the body says Gemini 2.5 Pro restructured them into steps. I report both.
- **Dates:** PhysMent's arXiv first-posted date (2026-07-07) sits oddly with its 2609 identifier; this is irrelevant to the selection.

## 8. Download record

`python docs/related_work_oct2026/fetch_papers.py --file physics_science.json`, run from the repository root on 2026-10-01: "27 fetched, 0 failed, 0 unresolved, 131 files in MANIFEST.json".

- The other three keys (imani2025sympybench, azmi2026scirho, ansari2026physicsregrading) were already on disk from symbolic_contamination, with the same arXiv sources, and printed "present".
- All 30 PDFs have text layers (12 to 64 pages each).
- No entry needed an id or url fix, and the script printed no errors.

One observation, with no change made to the script. `fetch_papers.py` reads MANIFEST.json at start and rewrites it at the end. If two agents run it at the same time, the later save can drop the other's new records. Each PDF stays on disk, and a later run re-records any PDF that lacks a manifest entry ("present ... text extracted"), so a final full run (`python fetch_papers.py`) after all agents finish would reconcile the manifest.
