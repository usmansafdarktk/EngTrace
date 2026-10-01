# Literature search: process supervision and step-level verification (code name `process_supervision`)

Searched 2026-10-01 for the October 2026 revision of EngTrace's Related Work. The output catalogue is
`docs/related_work_oct2026/candidates/process_supervision.json`, with 30 entries: 9 must-cite, 18 should-cite and 3 optional.
All 30 PDFs were downloaded by `fetch_papers.py --file process_supervision.json` (30 fetched, 0 failed, 0 unresolved).

Conventions used in the catalogue:

- The year in each key and in the `year` field is the year the paper was first posted (arXiv v1). This matches the
  brief's example key `zheng2024processbench`. The publication venue is in `venue`. Where the two years differ, the
  paper should be cited by the venue year (for example, ROSCOE is `golovneva2022roscoe` but is cited as ICLR 2023).
- Each venue was checked against an official page: the arXiv comments field, the ACL Anthology, the AAAI OJS, the ICLR
  proceedings, mlanthology, or the "Published in" header of the arXiv PDF. Venues that rest on a weaker source are
  listed under "What I could not verify".
- The PRM findings EngTrace has to place in this literature are the ones in `docs/NEXT_CYCLE_REVIEW.md` §4.2. The
  off-the-shelf 72B math PRM ranks steps well overall (AUROC 0.825), but three of its four flags inside correct-answer
  traces are false. VersaPRM finds 6% of incorrect steps. No deterministic evaluator detects a conceptual defect behind
  a correct answer (0 of 60 planted), while the best judges catch about a third.

---

## 1. Queries (verbatim)

### 1.1 WebSearch: 63 queries run, in order

1. ProcessBench identifying process errors in mathematical reasoning Zheng 2024
2. PRMBench fine-grained challenging benchmark process-level reward models
3. The Lessons of Developing Process Reward Models in Mathematical Reasoning Qwen2.5-Math-PRM
4. VersaPRM multi-domain process reward model synthetic reasoning data
5. Math-Shepherd verify and reinforce LLMs step-by-step without human annotations
6. ThinkPRM process reward models that think verification chain-of-thought
7. GenPRM scaling test-time compute of process reward models via generative reasoning
8. Generative Verifiers: Reward Modeling as Next-Token Prediction GenRM
9. ReasonEval evaluating mathematical reasoning beyond accuracy validity redundancy steps
10. ROSCOE suite of metrics for scoring step-by-step reasoning Golovneva
11. ReCEval evaluating reasoning chains via correctness and informativeness Prasad
12. MR-Ben meta-reasoning benchmark evaluating system-2 thinking first error step
13. MR-GSM8K meta-reasoning benchmark large language model evaluation
14. CriticBench benchmarking LLMs for critique-correct reasoning
15. LLMs cannot find reasoning errors but can correct them BIG-Bench Mistake
16. REVEAL chain-of-thought is as strong as its weakest link benchmark verifiers of reasoning chains
17. process reward model benchmark 2026 arXiv step-level error detection
18. Hard2Verify step-level verification benchmark open-ended frontier math
19. DeltaBench can large language models detect errors in long chain-of-thought reasoning
20. process reward model physics chemistry science domain step-level evaluation 2025
21. Solving math word problems with process- and outcome-based feedback Uesato
22. correct final answer flawed reasoning steps LLM evaluation false positives process errors 2025
23. verify intermediate computations chain-of-thought calculator deterministic arithmetic check LLM reasoning steps
24. unfaithful chain of thought reasoning models don't always say what they think
25. "Evaluating Step-by-step Reasoning Traces: A Survey" arXiv
26. GR-Ben general reasoning benchmark evaluating process reward models
27. Rethinking reward models for multi-domain test-time scaling outcome versus process
28. PhysicsEval step-level grading LLM physics solutions human expert annotation process evaluation
29. engineering problem solving LLM step-level process evaluation intermediate values reference solution 2026
30. Do VLMs Reason Like Engineers benchmark stage-wise evaluation
31. Socratic-PRMBench benchmarking process reward models systematic reasoning patterns
32. process reward models fail out-of-distribution generalization non-math domains analysis 2025
33. ICLR 2026 process reward model step verification benchmark paper
34. ACL 2026 process reward models error identification benchmark findings
35. NeurIPS 2025 process reward model evaluation step-level verifier paper
36. EMNLP 2025 step-level error detection reasoning chain human annotated benchmark
37. "Verifiable Process Reward Models for Structured Reasoning" deterministic step-level verification
38. PRIME process-outcome alignment benchmark verifiable reasoning mathematics engineering
39. "From Answers to States" verifiable process-level evaluation chemical reasoning large language models
40. Neuro-symbolic PRM scientific reasoning structured traces symbolic verification 2026
41. "LLMs as a Jury" cross-model consensus outperform process reward models reasoning
42. "The Answer Is Not the Argument" arXiv 2026 reasoning
43. MedPRMBench fine-grained benchmark process reward models medical reasoning
44. LLMs cannot spot math errors even when allowed to peek into the solution
45. "Reasoning Jury" multi-model consensus evaluating reasoning traces
46. SPARE single-pass annotation reference-guided evaluation automatic process supervision reward modelling
47. Evaluating LLMs at detecting errors in LLM responses ReaLMistake Kamoi
48. Evaluating research-level math proofs via strict step-level verification 2026
49. Improve mathematical reasoning in language models by automated process supervision OmegaPRM
50. Rewarding progress scaling automated process verifiers for LLM reasoning Setlur
51. Deductive verification of chain-of-thought reasoning natural program Ling NeurIPS 2023
52. reasoning models don't always say what they think arXiv 2505 Chen Benton faithfulness
53. CoT-Pass@k reinforcement learning verifiable rewards implicitly incentivizes correct reasoning base LLMs
54. VerifyBench systematic benchmark evaluating reasoning verifiers across domains
55. Big-Math large-scale high-quality math dataset reinforcement learning verifiable answers Albalak
56. physics process reward model step-level verifier PhysPRM PhysicsPRM 2026
57. PhysPRM generative process reward model fine-grained diagnosis physics problem solving PhysProcessBench arXiv
58. MathCheck evaluating mathematical reasoning with checklist robustness
59. Fin-PRM domain-specialized process reward model financial reasoning
60. process reward models biased hidden bias length position step correct answer traces false alarms 2026
61. LLM-as-a-judge versus process reward model step error identification comparison critic models outperform PRMs ProcessBench
62. Chain-of-Thought reasoning in the wild is not always faithful Arcuschin
63. Stepwise verification and remediation of student reasoning errors LLM tutors first error step annotation

A 64th query, `"Free Process Rewards without Process Labels" implicit PRM Yuan`, was refused. The session's WebSearch
budget (200 calls, shared across agents) had been used up, so that paper is listed below without a verified id.

### 1.2 arXiv API: 16 queries

Every query used the same URL template, run with a 3 s pause between calls and parsed with `xml.etree`:
`http://export.arxiv.org/api/query?search_query=<Q>&start=0&max_results=100&sortBy=submittedDate&sortOrder=descending`.
Only entries first posted in 2024 to 2026 were kept. The number after each query is how many of its 100 newest results
were from 2024 to 2026.

| # | search_query | 2024-26 hits |
|---|---|---|
| 1 | `all:"process reward" AND all:benchmark` | 100 |
| 2 | `all:"step-level verification"` | 19 |
| 3 | `all:"process supervision" AND all:evaluation` | 95 |
| 4 | `all:"error detection" AND all:reasoning AND all:step` | 44 |
| 5 | `all:"process reward model" AND (all:science OR all:physics OR all:chemistry OR all:engineering)` | 100 |
| 6 | `all:"intermediate steps" AND all:verification AND all:reasoning` | 45 |
| 7 | `all:"first error" AND all:reasoning` | 21 |
| 8 | `all:"process reward" AND (all:domain OR all:generalization)` | 100 |
| 9 | `all:"process reward" AND all:engineering` | 34 |
| 10 | `all:"reasoning traces" AND all:evaluation AND (all:expert OR all:human) AND all:step` | 48 |
| 11 | `ti:faithful AND (ti:"chain-of-thought" OR ti:reasoning)` | 100 |
| 12 | `all:"step-level" AND all:"LLM-as-a-judge"` | 7 |
| 13 | `ti:process AND ti:benchmark AND (ti:reasoning OR ti:error OR ti:reward)` | 18 |
| 14 | `all:intermediate AND all:arithmetic AND all:verification AND all:"chain-of-thought"` | 3 |
| 15 | `ti:verifier AND ti:benchmark` | 77 |
| 16 | `all:"correct answer" AND all:"flawed reasoning"` | 23 |

Queries 1 to 8 returned 359 unique 2024 to 2026 papers. I screened all of their titles by hand, and did the same for
the new titles from queries 9 to 16. The 2026 candidates in the catalogue were found this way and then checked on the
web: LSR-Ben, PRIME, ChemCoTBench-V2, Neuro-symbolic PRM, The Answer Is Not the Argument, and the PhysPRM lead.

### 1.3 Metadata checks

- **arXiv abstract pages (WebFetch).** I opened the abstract page of every arXiv candidate to confirm the title, the
  full author list, the v1 date and the comments or journal-ref fields. I also opened 2305.20050 (Lightman et al.)
  and 2606.10833 (EngVQA).
- **Venue pages (WebFetch).** For venues I opened the ACL Anthology pages of PRIME, Hard2Verify, DeltaBench, VPRM,
  the Qwen PRM lessons paper, Math-Shepherd and PhysPRM. I also opened the AAAI OJS page of GenPRM and the
  mlanthology page of ThinkPRM.
- **ICLR 2024 proceedings index.** I fetched it with curl to confirm the venue of Lightman et al.
- **Lookups that failed.** The arXiv API `id_list` check returned HTTP 429 five times in a row. dblp and OpenReview
  showed bot walls, and the Semantic Scholar API returned 429. I worked around all of these.
- **Numbers from the PDFs.** After the download, I checked the numbers quoted below against the PDFs' extracted text
  in `docs/related_work_oct2026/text/`.

---

## 2. Included candidates (30)

Each entry gives what is labelled, who labelled it, the domain, the size, the headline finding, and how the paper
bears on EngTrace.

### A. Process versus outcome supervision; the PRMs EngTrace tests

**uesato2022processoutcome** (must-cite). Uesato et al., *Solving math word problems with process- and outcome-based
feedback*, arXiv 2211.14275, Nov 2022.

- **What it does.** This is the first systematic comparison of process-based and outcome-based supervision, on GSM8K.
  Human annotators marked whether each step so far was correct, and "trace error" is the share of correct-answer
  solutions with any reasoning mistake, as judged by humans.
- **Finding.** Outcome supervision matches process supervision on final-answer error with fewer labels. Low trace
  error, however, needs process feedback or a reward model that emulates it. With it, final-answer error fell from
  16.8% to 12.7% and trace error among answer-correct solutions from 14.0% to 3.4%.
- **For EngTrace.** This paper is the origin of the measurement EngTrace scales up: a correct answer does not
  certify the trace. It is math only and needs human labels for every trace; EngTrace replaces those labels with gold
  intermediate values generated from code.

**wang2023mathshepherd** (must-cite). Wang et al., *Math-Shepherd*, ACL 2024 (pp. 9426-9439), arXiv 2312.08935.

- **What it does.** It trains a PRM with no human labels. A step is labelled by whether rollouts from that step still
  reach the correct answer (Monte Carlo estimation). The PRM is used for reranking and for step-by-step PPO on
  GSM8K and MATH.
- **Finding.** With PPO, Mistral-7B rose from 77.9% to 84.1% on GSM8K and from 28.6% to 33.0% on MATH.
- **For EngTrace.** It defines the automatic-labelling paradigm behind most open PRMs. The plan of record names it.
  EngTrace's gold traces give exact step labels, where Monte Carlo labels can only estimate them.

**zhang2025prmlessons** (must-cite). Zhang et al., *The Lessons of Developing Process Reward Models in Mathematical
Reasoning*, Findings of ACL 2025 (pp. 10495-10516), arXiv 2501.07301.

- **What it does.** The Qwen team releases Qwen2.5-Math-PRM-7B and -72B, the PRMs EngTrace ran.
- **Findings.** Monte Carlo step labels yield worse PRMs than LLM-as-judge or human labels, because completions can
  reach correct answers from incorrect steps and incorrect answers from correct steps. Best-of-N evaluation is biased:
  policies produce correct answers with flawed processes, and PRMs tolerate them. Existing PRMs concentrate their
  minimum scores on the final-answer step, which is a drift from process to outcome assessment. The fix the paper
  proposes is consensus filtering of Monte Carlo and LLM-judge labels, plus step-level metrics.
- **For EngTrace.** It must be cited because EngTrace evaluates its model. Its own diagnosis, that PRMs tolerate
  correct-answer-flawed-process responses, anticipates EngTrace's finding that the 72B PRM's flags inside
  correct-answer traces are mostly false.

**zeng2025versaprm** (must-cite). Zeng et al., *VersaPRM*, ICML 2025 (PMLR 267), arXiv 2502.06737.

- **What it does.** A multi-domain PRM trained on synthetic chains of thought for MMLU-Pro's 14 categories, one of
  which is Engineering. Llama-3.1-8B wrote the chains and Llama-3.1-70B-Instruct auto-labelled the steps.
- **How it is evaluated.** Only by reranking accuracy (weighted majority voting), never against human step labels.
  In Law it gains 7.9% over majority voting, against Qwen2.5-Math-PRM's 1.3%.
- **For EngTrace.** It must be cited because EngTrace ran it. EngTrace supplies what VersaPRM's evaluation lacks: a
  step-level check against expert labels, where VersaPRM finds 6% of incorrect steps. LSR-Ben finds the same pattern
  of missed errors (see B).

**lee2025rethinkingrm** (should-cite). Lee et al., *Rethinking Reward Models for Multi-Domain Test-Time Scaling*,
TMLR (07/2026) per the PDF header, arXiv 2510.00492.

- **What it does.** A unified evaluation, across 14 domains, of discriminative and generative outcome reward models
  and PRMs (dORM, dPRM, gORM, gPRM).
- **Findings.** The dORM matches the dPRM, the gPRM is not competitive, and the gORM is the most robust. Stepwise
  scoring inherits label noise from LLM-based auto-labelling, and stepwise aggregation compounds errors as traces
  grow longer.
- **For EngTrace.** This is the outcome-versus-process evidence outside math. It supports EngTrace's decision not to
  rely on learned step scorers, and to use deterministic checks with gold values instead.

### B. Step-level verification benchmarks with human labels

**zheng2024processbench** (must-cite). Zheng et al., *ProcessBench*, ACL 2025, arXiv 2412.06559.

- **What is labelled.** 3,400 model-generated solutions to GSM8K, MATH, OlympiadBench and Omni-MATH problems are
  labelled with the earliest erroneous step. Annotators were doctoral-level math experts, three to five per solution
  until three agreed; about 30% of solutions were discarded.
- **Findings.** Existing PRMs fail to generalise beyond GSM8K and MATH and underperform prompted critic LLMs. The
  share of correct-answer solutions that contain a process error grows with difficulty: 3.5% on GSM8K, 18.8% on
  MATH, 32.2% on OlympiadBench and 51.8% on Omni-MATH (Table 2).
- **For EngTrace.** This is the reference first-error benchmark, and the plan of record names it. The contrast with
  EngTrace: ProcessBench is math only, its labels are human labels on model outputs, and it measures verifiers.
  EngTrace is engineering, its gold intermediate values come from code, and it validates its evaluators against 15
  experts. ProcessBench's rising rate of process errors behind correct answers is the math analogue of EngTrace's
  argument.

**song2025prmbench** (must-cite). Song et al., *PRMBench*, ACL 2025 main, arXiv 2501.03124.

- **What is labelled.** 6,216 problems and 83,456 step labels, built on PRM800K. GPT-4o modified correct solutions
  to inject errors, followed by LLM and human filtering. Errors fall into three categories (simplicity, soundness,
  sensitivity) with nine subcategories.
- **Findings.** Across 25 models, the best (Gemini-2-Thinking) scores 68.8 against a human score of 83.8, and the
  authors report a significant inconsistency between step-level and outcome-level evaluation.
- **For EngTrace.** The fine-grained PRM benchmark, math only, with constructed errors. EngTrace's planted conceptual
  defects are a counterpart in engineering, scored deterministically where possible.

**sun2026lsrben** (should-cite). Sun et al., *LSR-Ben: A Logical and Scientific Reasoning Benchmark for Evaluating
Process Reward Models*, arXiv 2605.01203 (v1 2 May 2026, titled *GR-Ben: A General Reasoning Benchmark ...*; three
versions).

- **What is labelled.** 3,595 solutions are labelled with the first erroneous step and its error type. Three trained
  annotators label each solution, with double cross-validation; about 15% were discarded. The domains are science
  (biology, physics, chemistry, computer science) and five logic types.
- **Findings, from v3 Table 3.** Qwen2.5-Math-PRM-72B, the best PRM, falls from 78.3 F1 on ProcessBench to 37.4
  average F1 on LSR-Ben, and VersaPRM scores 34.1. PRMs tend to overlook errors, while LLMs over-flag them: in Table
  4, Qwen2.5-Math-PRM-72B and VersaPRM miss the first error in 56.7% and 54.4% of flawed solutions (FN%).
- **For EngTrace.** The most direct external corroboration that math PRMs, and VersaPRM, transfer poorly to
  non-math technical reasoning. It is a preprint (see section 7).

**pandit2025hard2verify** (should-cite). Pandit et al., *Hard2Verify*, ACL 2026 (pp. 22502-22517), arXiv 2510.13744.

- **What is labelled.** 1,860 steps across 200 responses that GPT-5, Gemini 2.5 Pro and Claude Sonnet 4 wrote to 80
  recent Olympiad problems. PhD-level math experts labelled them over more than 500 hours, with three rounds of
  review.
- **Findings.** 29 verifiers were evaluated. On first-error identification, Qwen2.5-Math-PRM-72B drops from 78.3 on
  ProcessBench to 37.3 (Fig. 1).
- **For EngTrace.** It shows that a PRM's score on ProcessBench does not carry over to harder, naturally occurring
  errors. This is useful context for the 72B PRM result.

**he2025deltabench** (should-cite). He et al., *Can Large Language Models Detect Errors in Long Chain-of-Thought
Reasoning?* (DeltaBench), ACL 2025 (pp. 18468-18489), arXiv 2502.19361.

- **What is labelled.** 1,236 long chains of thought from o1-like models, divided into sections. Humans annotated
  usefulness, correctness, the first error and reflection. The domains are math, programming, PCB (physics,
  chemistry, biology) and general reasoning.
- **Finding.** The best critic, GPT-4-turbo-128k, reaches only 40.8% F1, and the Qwen2.5-Math PRMs were among the
  PRMs tested.
- **For EngTrace.** It extends step verification to reasoning-model traces and science domains.

**zeng2023mrgsm8k** (should-cite). Zeng et al., *MR-GSM8K*, ICLR 2025, arXiv 2312.17080.

- **What is labelled.** 2,999 model solutions (1,428 correct, 1,571 incorrect) to original, program-of-thought and
  reversed GSM8K questions. Trained annotators under quality control labelled correctness, the first error step and
  the error reason.
- **Finding.** Models that look alike on GSM8K separate by more than 20 points when asked to verify rather than
  solve.
- **For EngTrace.** Early evidence that solving and verifying are different skills, and an example of first-error
  labels with stated reasons.

**zeng2024mrben** (should-cite). Zeng et al., *MR-Ben*, NeurIPS 2024, arXiv 2406.13975.

- **What is labelled.** 5,975 questions with model solutions. Human experts annotated correctness, the first error
  step, the error reason and a correction.
- **Domains.** Math, physics, chemistry, biology, medicine, logic and coding.
- **Finding.** Most models fall far below the o1 series on the MR-Score.
- **For EngTrace.** The nearest multi-domain step-error benchmark built on human labels. It has no engineering and no
  gold intermediate values.

**tyen2023bigbenchmistake** (should-cite). Tyen et al., *LLMs cannot find reasoning errors, but can correct them
given the error location* (BIG-Bench Mistake), Findings of ACL 2024, arXiv 2311.08516.

- **What is labelled.** 2,186 chain-of-thought traces from PaLM 2 on five BIG-Bench tasks (word sorting, shuffled
  objects, logical deduction, multi-step arithmetic, Dyck languages), annotated with the first logical mistake. Humans
  annotated four of the tasks; agreement is reported with Krippendorff's alpha.
- **Findings.** LLMs fail to find mistakes even when the mistakes are objective and unambiguous, but they can correct
  them once given the location.
- **For EngTrace.** The classic evidence that LLM judges are weak mistake-finders. This motivates EngTrace's
  deterministic first tier and its reporting of judge recall.

**jacovi2024reveal** (should-cite). Jacovi et al., *A Chain-of-Thought Is as Strong as Its Weakest Link* (REVEAL),
ACL 2024, arXiv 2402.00559.

- **What is labelled.** 1,002 chain-of-thought answers (3,360 steps) to four open-domain QA datasets. Each step is
  labelled for relevance, its type (attribution or logic), attribution to evidence passages, and logical correctness.
- **Finding.** Verifiers struggle, especially with logical correctness.
- **For EngTrace.** The step-level verifier benchmark outside math. It separates factual attribution from logical
  validity, a split EngTrace mirrors with milestones (values) versus the judge (concepts).

### C. Reference-free and model-based step evaluators

**golovneva2022roscoe** (must-cite). Golovneva et al., *ROSCOE*, ICLR 2023, arXiv 2212.07919.

- **What it is.** Unsupervised, interpretable step-by-step metrics for semantic alignment, semantic similarity,
  logical inference and language coherence.
- **How it was validated.** On five human-annotated datasets and six programmatically perturbed diagnostic datasets.
- **For EngTrace.** The canonical reasoning-trace metric suite. EngTrace's Section 4.3 still scores solutions
  against the gold trace with BERTScore and ROUGE-2, and ROSCOE is the reasoning-specific comparison a reviewer will
  expect, in either direction.

**prasad2023receval** (should-cite). Prasad et al., *ReCEval*, EMNLP 2023, arXiv 2304.10703.

- **What it is.** Reasoning-chain evaluation by correctness (natural language inference over reasoning content units)
  and informativeness (pointwise V-information).
- **How it was validated.** Meta-evaluated on perturbed and human-annotated errors on EntailmentBank, GSM-8K and
  DROP, where it beats ROSCOE-style baselines.
- **For EngTrace.** It pairs with ROSCOE as the reference-free baseline that EngTrace's reference-based, deterministic
  design replaces.

**xia2024reasoneval** (must-cite). Xia et al., *Evaluating Mathematical Reasoning Beyond Accuracy* (ReasonEval), AAAI
2025, arXiv 2404.05692.

- **What it is.** LLM-based scorers of step validity and redundancy, trained on PRM800K.
- **How it was validated.** Meta-evaluated on MR-GSM8K and MR-MATH. MR-MATH's labels are undergraduate annotators'
  labels on Abel and WizardMath solutions: 83 solutions with incorrect steps and 76 without.
- **Finding.** Higher final-answer accuracy in math-specialised LLMs does not imply better reasoning steps on hard
  problems.
- **For EngTrace.** The closest model-based analogue of EngTrace's process scores, and the same motivating finding.
  EngTrace's difference is that most of its checks are deterministic against gold values, and that experts validate
  them in engineering.

**lee2025stepsurvey** (should-cite). Lee and Hockenmaier, *Evaluating Step-by-step Reasoning Traces: A Survey*,
Findings of EMNLP 2025, arXiv 2502.12289.

- **What it is.** A taxonomy of trace-evaluation criteria under four headings: factuality, validity, coherence and
  utility. (The v3 abstract says "factuality"; an earlier version's abstract said "groundedness".) The survey also reviews evaluators and
  meta-evaluation datasets.
- **For EngTrace.** One citation that locates the milestone check (factuality and validity of intermediate values)
  and the arithmetic check (validity) in the field.

### D. Generative verifiers

**zhang2024genrm** (should-cite). Zhang et al., *Generative Verifiers: Reward Modeling as Next-Token Prediction*, ICLR
2025, arXiv 2408.15240.

- **What it is.** Verification trained as next-token prediction, optionally with chain-of-thought (GenRM-CoT).
- **Finding.** It beats discriminative verifiers and LLM-as-judge in Best-of-N: from 5% to 45.3% on algorithmic tasks
  and from 73% to 93.4% on GSM8K.
- **For EngTrace.** The paradigm behind using an LLM to verify reasoning, which is the role EngTrace's residual judge
  plays.

**khalifa2025thinkprm** (should-cite). Khalifa et al., *Process Reward Models That Think*, TMLR 2026 (per mlanthology
and the authors' repository), arXiv 2504.16828.

- **What it is.** A long chain-of-thought verifier fine-tuned on synthetic verification chains, filtered with about
  1% of PRM800K's labels.
- **Findings.** It beats LLM-as-judge and discriminative PRMs on ProcessBench, MATH-500 and AIME '24. Out of domain,
  on GPQA-Diamond and LiveCodeBench subsets, it beats discriminative verifiers by 8% and 4.5%.
- **For EngTrace.** The strongest recent alternative to the discriminative PRMs EngTrace tested. A reviewer may ask
  why it was not run.

**zhao2025genprm** (should-cite). Zhao et al., *GenPRM*, AAAI 2026 (vol. 40, no. 41), arXiv 2504.00891.

- **What it is.** A generative PRM that reasons and executes code to verify each step before judging it. It is
  trained on 23K MATH examples, with labels from relative progress estimation.
- **Finding.** On ProcessBench, a 1.5B GenPRM beats GPT-4o and a 7B GenPRM beats Qwen2.5-Math-PRM-72B.
- **For EngTrace.** The program-aided step check closest to EngTrace's per-claim arithmetic rule. The difference is
  that GenPRM's code is model-written, while EngTrace re-computes stated arithmetic deterministically and checks it
  against gold values.

### E. Process evaluation in science and engineering; deterministic step checks

**wang2026prime** (must-cite). Wang et al., *PRIME: A Process-Outcome Alignment Benchmark for Verifiable Reasoning in
Mathematics and Engineering*, ACL 2026 (pp. 14973-14985), arXiv 2602.11570.

- **What is labelled.** 2,530 trajectories on college-level STEM problems drawn from textbooks and exams; the source
  pool spans mathematics, physics, chemistry and engineering, in 16 subdomains. Problems were filtered by consistency
  and difficulty. Domain experts labelled each trajectory with an outcome label and a binary label for process-outcome
  consistency, at the solution level rather than the step level.
- **Findings.** Verifiers often miss derivation flaws behind correct answers ("spurious correctness"). The best
  verifier, Gemini-2.5-Pro-thinking, reaches 88.56% overall accuracy. Verifier accuracy on PRIME predicts RLVR gains
  (R² > 0.92).
- **For EngTrace.** The closest competitor, since it makes the same correct-answer, flawed-process claim with
  engineering in the title. The differences: EngTrace labels at the step level with gold intermediate values from
  code, is engineering-specific across five branches, and evaluates the solvers as well as the verifiers. PRIME's
  results tables report Math, Physics, Chemistry and Biology, so how many of its items are engineering is not stated
  in the text I read (see section 7).

**dong2026physprm** (should-cite). Dong et al., *PhysPRM*, Findings of ACL 2026 (pp. 9410-9427), ACL Anthology PDF; no
arXiv version found.

- **What it is.** A generative physics PRM that outputs critiques, judgments and error types. It is trained on
  PhysPRM30K (12K problems from eight sources, solutions from five LLMs), auto-labelled by Gemini-2.5-Flash against
  ground truth, with Monte Carlo rollouts used to confirm incorrect steps.
- **What is labelled.** PhysProcessBench, built from five test sets, has every step verified by human experts.
- **Findings.** Many open-source LLMs score about 50 F1 on PhysProcessBench because they label most steps correct.
  For example, Qwen2.5-7B-Instruct scores 71.6 F1 on correct steps and 14.2 on incorrect ones. PhysPRM averages 89.0.
- **For EngTrace.** The nearest science-domain process verifier, and a further instance of the positive bias, where
  errors go unflagged, that EngTrace saw in VersaPRM.

**guo2026chemcotbenchv2** (should-cite). Guo et al., *From Answers to States: Verifiable Process-Level Evaluation of
Chemical Reasoning in Large Language Models* (ChemCoTBench-V2), arXiv 2606.03660.

- **What it is.** 5,620 samples across 18 tasks: molecular understanding, editing, optimisation and reaction
  prediction. Models must expose intermediate steps in expert-designed templates.
- **How steps are checked.** With deterministic chemistry rules and, for closed-answer tasks, reference traces, "rather
  than another LLM judge". It reports answer correctness, template adherence and step-wise verifier correctness as
  separate signals.
- **Finding.** Frontier models often answer correctly with weak supporting reasoning.
- **For EngTrace.** The closest methodological twin of EngTrace's deterministic milestone checks, in chemistry. It
  should be cited so the paper can claim the engineering version rather than the idea.

**zi2026neurosymbolicprm** (should-cite). Zi et al., *Neuro-symbolic PRM: Enhancing Scientific Reasoning via
Structured Traces and Symbolic Verification*, arXiv 2608.26329.

- **What it does.** It splits step checking into symbolic validity and semantic groundedness. Symbolic validity is
  guaranteed by a deterministic verifier: the step must execute and be unit-consistent. Semantic groundedness is
  scored by a PRM trained on hard negatives that pass the verifier but are wrong ("valid math, wrong variable").
- **How it is evaluated.** On ProcessBench, PRMBench, and math and science benchmarks; the PRM is trained from
  PRM800K queries.
- **For EngTrace.** It names the exact residual class EngTrace measured. Applying the right formula to the wrong
  variable passes every deterministic check, and EngTrace found that no deterministic evaluator catches such
  conceptual defects behind a correct answer. It supports EngTrace's tiered design.

**pronesti2026vprm** (optional). Pronesti, Belz and Hou, *Beyond Outcome Verification: Verifiable Process Reward
Models for Structured Reasoning*, Findings of ACL 2026 (pp. 32187-32202), arXiv 2601.17223.

- **What it does.** An RL framework in which each step is checked by a deterministic verifier grounded in guidelines,
  against gold labels. It is applied to risk-of-bias assessment for medical evidence synthesis.
- **Finding.** Up to 20% higher F1 than state-of-the-art models, and 6.5% higher than verifiable outcome rewards.
- **For EngTrace.** A training-side counterpart to EngTrace's deterministic step checks. Cite it if the paper
  mentions using EngTrace's checks as rewards.

**rizvi2025spare** (optional). Rizvi, Zhu and Gurevych, *SPARE*, AAAI 2026 (oral), arXiv 2506.15498.

- **What it does.** Single-pass step annotation that aligns each solution step to one or more reference-solution steps
  and judges it, on GSM8K, MATH, MuSiQue-Ans and SpaRP.
- **Findings.** It matches MCTS-based labelling at 2.3 times lower token cost and generalises to ProcessBench with
  about 16% of the training samples.
- **For EngTrace.** The LLM-based counterpart to EngTrace's gold-milestone coverage, since both are reference-guided
  step evaluation. EngTrace's gold steps are machine-generated and its alignment is deterministic.

### F. A correct answer over an unsound process; reference-guided judging

**yeadon2026answerargument** (should-cite). Yeadon et al., *The Answer Is Not the Argument*, arXiv 2609.00264 (31 Aug
2026).

- **What is labelled.** 237 step-numbered solutions from three frontier models to 79 Humanity's Last Exam physics
  questions, with no inserted errors. Each is labelled for answer correctness and first false step. The reference
  standard combines physicist annotation, an LLM debate and source-masked adjudication. 24 traces are "critical":
  the answer is correct but the trace contains a genuine error.
- **Findings.** Eight LLM monitors were tested. Giving a certified answer raised balanced accuracy from 0.637 to
  0.796, and raised recall on wrong-answer traces from 0.653 to 0.951. On critical traces it lowered recall from 0.521
  to 0.438, and the direction held for all eight monitors.
- **For EngTrace.** Independent evidence, in physics, for EngTrace's correct-answer blind spot. It also bears on
  EngTrace's finding that the expert verdict is nearly determined by the answer (AUROC 0.974), and it warns that a
  judge shown the gold answer checks conclusions rather than arguments.

**srivatsa2025peek** (should-cite). Srivatsa, Maurya and Kochmar, *LLMs cannot spot math errors, even when allowed to
peek into the solution*, EMNLP 2025, arXiv 2509.01395.

- **What it does.** First-error-step localisation on VtG (teacher-annotated student solutions; Daheim et al. 2024) and
  PRM800K.
- **Finding.** State-of-the-art LLMs struggle even when given the reference solution. Generating an intermediate
  corrected student solution helps.
- **For EngTrace.** Directly relevant to the residual judge, which sees gold material. Access to a reference does not
  make an LLM a reliable step verifier, which is why EngTrace reports judge recall on planted defects.

### G. Faithfulness of traces

**turpin2023unfaithful** (optional). Turpin et al., *Language Models Don't Always Say What They Think*, NeurIPS 2023,
arXiv 2305.04388.

- **Finding.** Chain-of-thought explanations can systematically misstate why a model answered. With biased contexts,
  accuracy fell by up to 36% on 13 BIG-Bench Hard tasks, and the bias went unmentioned.
- **For EngTrace.** A Limitations sentence: EngTrace scores the written trace, not the model's internal computation.
  Chen et al. 2025 is the reasoning-model alternate (see section 4).

---

## 3. Positioning in one place (for whoever drafts the paragraph)

- **Prior benchmarks.** Process supervision began as human step labels (Uesato; PRM800K, already cited) and
  automatic Monte Carlo labels (Math-Shepherd). Step-level verification benchmarks with human first-error labels are
  ProcessBench, PRMBench, MR-GSM8K, MR-Ben, BIG-Bench Mistake, REVEAL, DeltaBench and Hard2Verify. They are almost all
  math, or open-domain QA and logic, and they label model outputs after the fact.
- **The nearest technical-domain work.** PRIME labels process-outcome consistency per solution, in mathematics and
  engineering. LSR-Ben and PhysPRM label steps in science, and ChemCoTBench-V2 uses deterministic chemistry
  verifiers.
- **What EngTrace adds.**
  - Gold intermediate values generated by the same code as the problem, so steps can be checked deterministically
    rather than estimated (Monte Carlo) or judged (LLM).
  - Engineering coverage across five branches.
  - Evaluator validation against 15 experts' step labels.
  - Evidence that agrees with the external literature on limits:
    - Math PRMs transfer poorly (LSR-Ben, Hard2Verify, PhysPRM's positive bias).
    - Deterministic checks cannot see conceptual defects behind correct answers (the residual class the
      Neuro-symbolic PRM targets).
    - Judges shown the answer check conclusions rather than arguments (Yeadon et al.; Srivatsa et al.).
- **Two framing cautions.** Lee et al. (TMLR 2026) find that generative outcome verification beats process scoring
  across domains, so the paper should not claim that process scoring is better in general. Its claim is the narrower
  one about what answer checks miss. PRIME already uses the phrase "spurious correctness", so EngTrace should cite it
  rather than appear to coin the idea.

---

## 4. Considered and not included (one line each)

- OmegaPRM (Luo et al. 2024, arXiv 2406.06592): MCTS-automated step labels; training-side; Math-Shepherd and the Qwen lessons paper cover automatic labelling.
- Rewarding Progress / PAV (Setlur et al., arXiv 2410.08146, ICLR 2025): progress-based process verifiers for search and RL; training-side.
- AlphaMath Almost Zero (arXiv 2405.03553): process supervision without process labels; training-side.
- PRIME, *Process Reinforcement through Implicit Rewards* (Cui et al., arXiv 2502.01456): RL training method; its name collides with the PRIME benchmark (wang2026prime), so take care in the bibliography.
- Free Process Rewards without Process Labels (Yuan et al., implicit PRM): training-side; id not verified (search budget exhausted).
- Is PRM Necessary? (arXiv 2505.11227): RL implicitly induces PRM ability; training-side.
- CriticBench (Lin et al., arXiv 2402.14809, Findings ACL 2024): solution-level critique and correction over five domains; MR-Ben and BIG-Bench Mistake cover step-level critique.
- ReaLMistake (Kamoi et al., arXiv 2404.03602, COLM 2024): binary error detection at the response level, not step level.
- Daheim et al. 2024 (arXiv 2407.09136, EMNLP 2024): teacher-annotated first errors in about 1K student math solutions (VtG); included indirectly through Srivatsa et al., which uses VtG.
- Socratic-PRMBench (arXiv 2505.23474): math PRM benchmark organised by reasoning pattern; redundant with PRMBench.
- MedPRMBench (arXiv 2604.17282): medical PRM benchmark with 14 error types and 4-level severity grading; out of domain; first alternate if the paper discusses severity.
- Med-PRM (EMNLP 2025), Fin-PRM (arXiv 2508.15202): domain PRMs (medicine, finance); out of domain; Fin-PRM only if the FinChain sibling is discussed.
- Sci-PRM (arXiv 2606.04579): tool-aware scientific PRM; PhysPRM and the Neuro-symbolic PRM cover science-domain step verification.
- Process Reward Models Meet Planning (arXiv 2604.17957, ACL 2026), FoVer (arXiv 2505.15960): step labels from PDDL or formal verification; training-side.
- GroundedPRM (2510.14942), BiPRM (2508.01682, ACL 2026), R-PRM (2503.21295, EMNLP 2025), DG-PRM (2507.17849, ACL 2025), retrieval-augmented PRM (2502.14361), learned-reliability PRM (2605.15529), unsupervised PRMs (2605.10158), distributional PRMs (2605.06785), Cliff (2609.02817), progression-aware error localisation (2609.33297): PRM training or localisation methods with no new evaluation evidence relevant here.
- The Hidden Bias of PRMs / PRISM (arXiv 2606.09078): PRMs over-credit plausible incorrect steps because step labels are imbalanced (73.1% of PRM800K steps are correct); a candidate mechanism for low PRM recall; first alternate if the paper explains the PRM failure.
- EST-PRM (2606.00437), quality-diversity stress tests for PRMs (2608.08008): PRM stress-testing; cap.
- LLMs as a Jury (arXiv 2607.10139): cross-model answer consensus beats trained verifiers off-distribution; answer selection; belongs to the judge topic.
- Reasoning Jury (arXiv 2608.12585): a jury of open-weight LLMs with step-grounded verdicts beats frontier single judges at finding reasoning defects; belongs to the LLM-judge validity topic (PoLL is cited); flagged for that agent.
- Strict step-level verification of research-level proofs (arXiv 2606.10799), fine-grained evaluation of natural-language proofs (arXiv 2510.13888): proof grading, not numeric engineering derivations.
- MPBench (2503.12505), ViLBench (2503.20271), PRMBench-V (ACL 2026), PRISM-Bench (2510.23594), ErrorRadar (2410.04509), VisualPRM (2503.10291), AudioProcessBench (2606.09925): multimodal or audio process benchmarks; out of scope.
- ToolPRMBench (2601.12294), AgentProcessBench (2603.14465), ToolComp (2501.01290), MAS-ProVe (2602.03053), ParaRecover (2609.12345), CUARewardBench (2510.18596): agent and tool trajectories; out of scope.
- VerifyBench (Li et al., arXiv 2507.09884, AAAI 2026) and VerifyBench (Yan et al., arXiv 2505.15801): answer-level verifiers, with expert labels in the first case; better placed with the answer-checking discussion; the two papers share a name.
- RewardBench 2 (arXiv 2506.01937, ICLR 2026), SCI-Verifier (arXiv 2509.24285): outcome or answer verifiers.
- Big-Math (arXiv 2502.17387): 250K RL questions with verifiable final answers; outcome-only verification, a contrast rather than a comparator.
- CoT-Pass@k (Wen et al., arXiv 2506.14245) and *Does CoT-Pass@k Really Check the CoT?* (arXiv 2609.32622): an LLM-verified chain-of-thought correctness metric for RLVR; relevant but RL-focused.
- MathCheck (arXiv 2407.08733, ICLR 2025): checklist that includes a process-judging task on GSM8K and geometry; robustness focus.
- Deductive Verification of CoT (Ling et al., arXiv 2306.03872, NeurIPS 2023): an LLM checks each step given only its premises; first alternate for a "verify each claim in isolation" sentence.
- VeriCoT (arXiv 2511.04662), CoT verification via computational graphs (arXiv 2510.09312), SymCode (arXiv 2510.25975): logic-consistency, white-box or generation-side methods.
- Faithfulness: Chen et al. 2025 *Reasoning Models Don't Always Say What They Think* (arXiv 2505.05410; the alternate to Turpin for reasoning models); Arcuschin et al. 2025 (arXiv 2503.08679); FaithCoT-Bench (2510.04040); C2-Faith (2603.05167); RFEval (2602.17053); *Right Answer, Wrong Reason* (clinical, 2609.32817); *Legibility is Not Interpretability* (2609.04194). Faithfulness is tangential here, so one optional cite.
- Better Accuracies, Worse Reasoning (arXiv 2605.28301): step-level audit of medical chain-of-thought distillation (accuracy up, step quality down), parallel to ReasonEval; cap.
- The Correct Answer Trap (2606.23205), Right for Wrong Reasons in small models (2601.00513), Filtered Reasoning Score (2604.11996), verifying reasoning in self-improvement (2603.21558): the answer-process gap in tutoring, agent or training settings; cap.
- PhysicsEval (arXiv 2508.00079, Findings IJCNLP 2025), PhysicsArena (2505.15472): physics benchmarks with LLM-scored or multimodal process dimensions; not step-level verification.
- Do VLMs Reason Like Engineers? / EngVQA (arXiv 2606.10833): 696 multimodal engineering problems and an 8-stage automatic evaluator (Pearson 0.975 with human grades); belongs in the process paragraph too, but engineering_web.json already catalogues it under the same key, `wasiq2026engvqa`.
- StepORLM (arXiv 2509.22558): generative process supervision for operations-research modelling; training-side; possibly relevant to the industrial branch.
- RusFinChain (arXiv 2607.01388): a FinChain-style benchmark; belongs to the symbolic-benchmark topic.
- A Survey of PRMs (arXiv 2510.08049), Trust but Verify survey (arXiv 2508.16665): surveys; Lee and Hockenmaier was chosen as the single survey, with the PRM survey as alternate.
- OPV (arXiv 2512.10756), Reward Granularity in RLVR (arXiv 2607.02869): verification efficiency, and process versus outcome rewards in RL training; cap.

---

## 5. Overlaps with other catalogue files

- `wang2026prime` (2602.11570) was also in `engineering_web.json`, under the same key, earlier in the session. At my
  final check it was no longer there, so this file is now its only catalogue entry.
- `wasiq2026engvqa` (2606.10833) is in `engineering_web.json` and not in mine; see section 4.
- At the final check nothing in this file overlapped, by key or by arXiv id, with `papers.json`,
  `engineering_web.json`, `llm_judge.json`, `physics_science.json` or `symbolic_contamination.json`.

---

## 6. Updates to existing citations

**Lightman et al. 2023, *Let's Verify Step by Step* (cited set, key `lightman2023verify`).**

- **Venue.** The May reference says "ArXiv preprint, abs/2305.20050". The paper was published at **ICLR 2024** (The
  Twelfth International Conference on Learning Representations). It is listed in the ICLR 2024 proceedings at
  `proceedings.iclr.cc/paper_files/paper/2024/hash/aca97732e30bcf1303bc22ac3924fd16-Abstract-Conference.html`.
  `papers.json` already records ICLR 2024; the bibliography entry needs updating.
- **Author list (a correction, not just an update).** The May reference lists 30 authors. The arXiv page and the ICLR
  2024 proceedings list 10: Hunter Lightman, Vineet Kosaraju, Yura (Yuri) Burda, Harri (Harrison) Edwards, Bowen
  Baker, Teddy Lee, Jan Leike, John Schulman, Ilya Sutskever, Karl Cobbe.
  - Only 7 of the May list's 30 names are authors: Lightman, Kosaraju, Burda, Edwards, Leike, Schulman, Sutskever.
  - These 23 names are not on the paper: Lukas Mesnard, Tyna Wang, Farzad Khorrami, Nguyet Minh Nguyen, Shayne
    Mostyn, Max Miller, Chia Hsuan Wang, Sam Gelman, Denis Igor, Marina Polozov, Girish Sastry, Prachit Tara, Sandhini
    Agarwal, Nan Rosemary Sun, Shengjia Zhao, Jeffrey Wu, Szymon Sidor, Jiayi Weng, Yuan Cao, Adrià Puigdomènech Badia,
    Nikolas Tezak, Peter Welinder, Wojciech Zaremba.
  - These true authors are missing from it: Bowen Baker, Teddy Lee, Karl Cobbe.
  - A reviewer who checks a reference will see this, so the entry should be regenerated from the ICLR BibTeX.
- **Wording.** The in-text use ("Inspired by recent advances in process supervision (Lightman et al., 2023)") can stay,
  with the year 2024 if the paper cites by venue year. Under the brief's key convention the key stays
  `lightman2023verify`.

No other entry in `papers.json` belongs to this topic.

---

## 7. What I could not verify

- **PRIME's engineering share.** The title says "Mathematics and Engineering" and the source pool includes
  engineering, but the main results tables report Math, Physics, Chemistry and Biology columns. I did not find the
  number of engineering items among the 2,530 in the text I read; it may be in Figure 3, which is not in the text
  layer. This matters for the "first engineering process benchmark" framing: do not claim that PRIME has no
  engineering, and do not claim that it has a lot.
- **LSR-Ben's status.** The v3 PDF header reads "Published as a conference paper at ICLR 2027", which cannot reflect
  an acceptance on 2026-10-01. I treat it as an arXiv preprint. Model names changed between versions (v3 uses
  Gemini-3.8-flash and GPT-5.5), so the numbers quoted above are from v3, the downloaded version. The v1 title, GR-Ben,
  is taken from the arXiv v1 HTML page as returned by search; I did not open the v1 PDF.
- **ThinkPRM's venue (TMLR 2026)** rests on mlanthology and the authors' GitHub README. The arXiv PDF carries no TMLR
  header.
- **Rethinking RMs' venue (TMLR 07/2026)** rests on the arXiv PDF header ("Published in Transactions on Machine
  Learning Research (07/2026)"). I could not open OpenReview (bot wall).
- **PhysPRM** has no arXiv version that I could find; the PDF is the ACL Anthology one. The `date` 2026-07 is the ACL
  2026 Findings publication month, not a first-posting date.
- **Uesato et al. 2022.** Search results also show a NeurIPS 2022 MATH-AI workshop version (mathai2022.github.io). I
  did not open it, so the catalogue venue stays "arXiv 2211.14275".
- **Free Process Rewards without Process Labels** was not searched, because the web-search budget was exhausted.
- **Numbers from summaries.** The WebFetch summary of 2609.00264 got the direction of a recall change wrong ("increased
  from 0.521 to 0.438"). The PDF says recall on critical traces fell from 0.521 to 0.438, and I used the PDF.
  Srivatsa et al.'s exact-step accuracies (about 65% and 51%) appeared only in a search summary and are not used.
- **Search coverage.** The arXiv API caps each query at 100 results sorted by date, so the broad queries (1, 5, 8, 11)
  reach back only to early 2026 for some terms. Older work was covered by the targeted WebSearch queries.
- **The paper's own preprint is indexed.** arXiv 2511.01650, *EngTrace: A Symbolic Benchmark for Verifiable Process
  Supervision of Engineering Reasoning*, appeared in the API results (query 1). It is not a candidate. Whoever prepares
  the ARR version should keep the anonymity rules in mind when citing it.

---

## 8. Download status

`python docs/related_work_oct2026/fetch_papers.py --file process_supervision.json` reported 30 fetched, 0 failed and 0
unresolved. Every PDF validated, with a text layer of 46,641 to 284,863 characters. PhysPRM came from the ACL
Anthology; all the others came from arXiv.

There is a shared-file problem. `fetch_papers.py` reads `MANIFEST.json` at start and writes it at the end, so runs by
several agents at once overwrite each other's records:

- When my run finished (17:09:46Z), the manifest held my 30 records and 2 cited ones, while 27 PDFs from other runs
  were on disk with no record.
- Ten minutes later, another agent's run (manifest stamped 17:19:27Z, 131 files) dropped all 30 of my records.
- My 30 PDFs and texts are still on disk (`--list --file process_supervision.json` shows "have" for every entry).

Nothing needs re-downloading. Once every agent has finished, one `python docs/related_work_oct2026/fetch_papers.py`
run re-records every PDF on disk through its "present" path, then `--verify` confirms the result. On the
re-recorded entries, `retrieved_utc` will show the time of that run, not the time of the original download. I did not
edit the script or the manifest, and I did not re-run, because a re-run now would only start another round of
overwrites.
