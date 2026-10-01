# Literature search: symbolic, templated and dynamic benchmarks, contamination and saturation (code name `symbolic_contamination`)

Searched 2026-10-01 for the October 2026 revision of EngTrace's Related Work. The output catalogue is
`docs/related_work_oct2026/candidates/symbolic_contamination.json`, with 30 entries: 9 must-cite, 16 should-cite and
5 optional. All 30 PDFs and text files were fetched by `fetch_papers.py --file symbolic_contamination.json` (two runs:
30 fetched, then 1 fetched after a swap; 0 failed, 0 unresolved).

Conventions:

- The year in each key and in the `year` field is the year of arXiv v1, as in `mirzadeh2024gsmsymbolic`. The venue
  field holds the publication venue where one exists. Where the two years differ, cite by the venue year.
- `date` is the v1 submission month shown on the arXiv abstract page. Two ids carry a later month than their v1 date,
  probably because of moderation holds: VeRA (2602.13217, v1 23 Jan 2026) and Xiao and Cheng (2609.02899, v1 5 Jul 2026).
- Every included paper's arXiv abstract page was opened. Numbers quoted below come from those pages or from the
  downloaded full text in `docs/related_work_oct2026/text/`.
- PRIME (2602.11570, ACL 2026), ChemCoTBench-V2 (2606.03660) and OPT-Engine (2601.19924, ICML 2026) fall in this topic
  too, but they are already in `engineering_web.json` and/or `process_supervision.json` (ChemCoTBench-V2 under the
  same key, `guo2026chemcotbenchv2`). I left them out to avoid duplicates. They are discussed in section 3.
- Three of my entries were also chosen, under the identical key and arXiv id, by agents whose files appeared after mine:
  `imani2025sympybench` and `ansari2026physicsregrading` in `physics_science.json`, and `azmi2026scirho` in
  `engineering_arxiv.json`. They are the same papers, so a full `fetch_papers.py` run will report them as DUPLICATE and
  keep the earlier file's record. I kept them because they are central to this topic.

---

## 1. Answers to the four positioning questions

**Which generated benchmarks pair instances with gold intermediate traces, not just final answers?**

- **Finance:** FinChain (cited), RusFinChain (2026) and V-FiLLM (2026, not included).
- **Physics:** SymPyBench (December 2025; EACL 2026 Industry Track).
- **STEM, visual and multilingual:** Sci-Rho (June 2026).
- **Grade-school math training data:** TemplateGSM.
- **Final-answer scoring only:** DyVal and GSM-Infinite construct every intermediate node of a computation graph, but they
  score only the final answer. GSM-Symbolic, MATH(), Putnam-AXIOM, UGMathBench, VeRA, OMEGA, RE-IMAGINE, GSM-Plus and
  MATH-Perturb also score final answers. Some of them ship worked or reference solutions, but none checks the model's
  intermediate steps.

**Which evaluate process?**

- **Step alignment with the gold chain:** FinChain (ChainEval) and RusFinChain (fuzzy numeric alignment).
- **LLM-judged steps:** Sci-Rho (step precision, recall and F1, judged by GPT-5.4 against the gold steps).
- **Deterministic oracle:** Milestone Oracles uses a deterministic symbolic verifier, but on milestones that a teacher
  model writes.
- **Error separation:** GSM-Ranges separates logical from arithmetic errors by re-executing the model's solution as code.
- **Answer-level only:** SymPyBench ships step-by-step gold reasoning but scores only answers and consistency across
  variants.
- **Other agents' entries:** ChemCoTBench-V2 checks intermediate chemical states with deterministic rules, but it is not
  parameterised instance generation. PRIME evaluates verifiers on process-outcome alignment over curated STEM problems.

I found no generated benchmark that combines code-derived milestones, a deterministic per-claim arithmetic check and an
LLM judge validated against expert step labels. That is EngTrace's distinguishing combination.

**What does the latest evidence say about contamination and paraphrase sensitivity?**

- **Contamination is real but uneven.**
  - GSM1k found drops of up to 8%, mostly in some open families.
  - MathArena found "strong signs of contamination in AIME 2024".
  - Wu et al. found that Qwen2.5-Math-7B regenerates the rest of 54.60% of MATH-500 problems from their first 60%.
  - Schaeffer et al. (2026) show in controlled pretraining that one test-set replica lowers loss below the clean floor.
- **Rewriting items does not decontaminate them.** None of 20 mitigation strategies beats the unmodified benchmark (Sun
  et al., ICML 2025).
- **On public benchmarks, paraphrasing lowers absolute scores while rankings stay stable.**
  - Lunardi et al. (ECAI 2025): 34 models, six multiple-choice benchmarks.
  - Xiao and Cheng (2026): rank correlation 0.997 between the standard and the paraphrase-controlled leaderboard.
  - Semantic variants cut accuracy more, about 28% on average in GSM-SEM's strictest setting.
  - Sampling variations at random can understate fragility by up to 5 times (LPDS).
- **GSM-Symbolic's headline drops are partly statistical.** In a re-analysis with mixed models, 8 of 20 models show
  significant drops, and the variants shift integers to larger values (Długosz et al., EMNLP 2026).
- **Saturation.** Nearly half of 60 LLM benchmarks saturate, more with age. Once age is controlled, neither private test
  sets nor templating significantly slows saturation (Akhtar et al., ICML 2026). Leading physics benchmarks are near
  saturation once grader and reference errors are fixed (Ansari et al., September 2026).

**Does any 2026 work do "symbolic templates with verifiable traces" in a technical domain other than finance?**

Yes, so EngTrace can no longer claim to be the only one outside finance.

- **SymPyBench:** physics, 15,045 problems, executable Pint/SymPy code with step-by-step reasoning. It was posted in
  December 2025 and published at EACL 2026.
- **Sci-Rho:** 606 expert-written templates per language across mathematics, physics, chemistry, biology and computer
  science, with gold reasoning steps and step-level F1.
- **Neither covers engineering.** SymPyBench scores only answers. Sci-Rho's step scores come from a single LLM judge, and
  the full text I searched reports no agreement between that judge and human labels.

The defensible claim is narrower: the first such benchmark for engineering branches, with deterministic, code-derived
process checks and an expert-validated residual judge. Even that "first" should be hedged as "to our knowledge".

---

## 2. Included candidates (30)

Each entry gives what the paper generates or measures, its size, whether it has gold intermediate traces, how it
evaluates, its key finding, and how it bears on EngTrace.

### 2.1 Symbolic or templated generation with gold traces (closest analogues)

**`imani2025sympybench`: SymPyBench** (Imani et al.; arXiv 2512.05954, v1 Dec 2025; EACL 2026 Industry Track,
pp. 105-118). **must-cite.**

- **Construction.** SymPyBench builds 15,045 university-physics problems (90/10 train/test) from openly licensed
  (Creative Commons) problem sets. An LLM (Llama 3.2/3.3) restates each problem into a schema: question, step-by-step
  reasoning, input and output variables with units, and constants. It then writes a parameterised template for the
  question and the solution, three textual rephrasings and Python code using Pint and SymPy. The code is kept only if
  it reproduces the original answer, which about 88% did, and all problems were manually reviewed.
- **Gold traces.** Every instance carries templated step-by-step reasoning and executable code, so gold traces exist.
- **Evaluation.** Scoring is answer-level only: exact match over all sub-questions and partial accuracy. Across numeric,
  textual and format variants it adds a Consistency Score, a Complete Failure Rate and a Confusion Rate. The best
  Consistency Score was 42.42% (Sonnet-3.7).
- **Bearing.** This is the nearest technical-domain precedent: the same recipe, it predates EngTrace's October revision,
  and it even includes a rephrasing axis. EngTrace differs in four ways. It covers engineering branches. It scores
  process deterministically (milestone coverage, per-claim arithmetic). Its judge is validated against expert step
  labels. And its items come from a private seed, whereas SymPyBench's seed problems are public. The Related Work must
  name SymPyBench and drop any "first outside finance" wording.

**`azmi2026scirho`: Sci-Rho** (Azmi et al.; arXiv 2606.08034, v1 Jun 2026, v2 Sep 2026). **must-cite.**

- **Construction.** 4,242 Python templates (606 per language, seven languages) were written by domain experts,
  including Olympiad medallists. They cover mathematics, physics, chemistry, biology and computer science (199, 195,
  137, 23 and 52 per language). The templates vary numbers, visual patterns, shapes, colours and function types, giving
  42,420 instances, each paired with reasoning steps and a ground truth.
- **Evaluation.** 17 VLMs are scored on average versus worst-case accuracy (a template counts only if every variant is
  right) and on step recall, precision and F1. A GPT-5.4 judge scores the steps against the gold steps.
- **Findings.** Worst-case accuracy is well below average accuracy, and smaller models degrade across languages. Step
  scoring was run on a subset of five models. Even for GPT-5.4, worst-case step F1 is 13.0 to 20.2 points below average
  step F1, depending on the subject.
- **Bearing.** This is the 2026 paper closest to "symbolic templates with verifiable traces outside finance". Its step
  scoring rests on one LLM judge with no reported human validation. EngTrace's code-derived milestones and arithmetic
  check are deterministic, and its residual judge is validated on 300 expert-labelled traces. Its worst-case-over-variants
  metric is a natural comparison for EngTrace's per-template analysis.
- **Authorship caution.** Sci-Rho says its templates are expert-crafted. Per `engtrace_state_brief.md`, EngTrace's
  templates were agent-written from textbooks and then expert-certified, so the paper should not borrow "authored by
  domain experts".

**`arabov2026rusfinchain`: RusFinChain** (Arabov; arXiv 2607.01388, v1 Jul 2026, v2 Aug 2026). **should-cite.**

- **Construction.** A Russian-language FinChain counterpart: 17 domains, 172 topics and 5,280 parameterised examples
  from executable Python templates. Each example has a gold reasoning chain with intermediate numeric values.
- **Evaluation.** It adds Fuzzy Numeric Alignment and Soft-Attention Alignment to ChainEval and tests 8 open-weight
  models (8,100 responses).
- **Finding.** Step alignment is fairly high (Hard F1 about 0.65), but only about 29% of final answers are correct. The
  fuzzy metrics correlate better with answer correctness (ρ about 0.48) than ChainEval does (about 0.38 to 0.46).
- **Bearing.** It shows the FinChain design spreading, and it gives the metric design space in which EngTrace's
  order-free, tolerance-based milestone coverage sits. Caveats: single author, arXiv only.

**`zhang2024templategsm`: TemplateGSM / Template-based Data Generation** (Zhang; arXiv 2411.18104, v1 Nov 2024,
v6 May 2026; ICLR 2025 DATA-FM Workshop). **optional.**

- **Construction.** GPT-4 writes parameterised meta-templates that generate over 7 million grade-school problems, each
  with a programmatically verifiable solution in code and natural language.
- **Use.** It is aimed at SFT and RLVR data rather than evaluation.
- **Bearing.** It is the training-data counterpart of EngTrace's code-derived gold solutions. Cite it only where
  template-generated data with verifiable solutions is discussed.

### 2.2 Functional, dynamic and symbolic-variant benchmarks

**`srivastava2024functionalbench`: Functional Benchmarks / MATH()** (Srivastava et al.; arXiv 2402.19450, Feb 2024).
**must-cite.**

- **Construction.** The authors hand-rewrote the reasoning of 41.2% of MATH (2,060 of 5,000 problems) into code, so
  that seeded snapshots give new instances of each problem. Three snapshots were released.
- **Finding.** The "reasoning gap" between static and functional accuracy was 58.35% to 80.31% for the models tested.
- **Gold traces.** The functional code encodes each problem's reasoning but outputs only the answer, and only answers are scored.
- **Bearing.** This is where "functional variants" and the reasoning-gap metric come from. Putnam-AXIOM, which EngTrace
  cites, adopts both, and EngTrace's templates are functional in exactly this sense.

**`cheng2026vera`: VeRA** (Cheng et al.; arXiv 2602.13217, v1 Jan 2026, v2 Sep 2026). **must-cite.**

- **Construction.** Each seed item becomes an executable task family: a natural-language template, an input generator and
  a deterministic answer program. Execution checks, seed anchoring (the program must reproduce the seed answer), answer
  discrimination and independent human solving validate each family. VeRA-E draws fresh instances, while VeRA-H and
  H Pro make harder modifications.
- **Seeds.** GSM8K (1,319 seeds, 2,638 variants), AIME 2024 and 2025 (30 seeds and 60 variants each), Beyond-AIME and
  AMO-Bench.
- **Findings.** Across 16 models, AIME-2024 accuracy falls from 84.46% on seeds to 70.25% on fresh VeRA-E instances, and
  on AIME-2024-II to 58.57% for H Pro. Audits accept 75.4% of hardened candidates, rising to 95.1% after repair.
- **Gold traces.** Answer programs only.
- **Bearing.** VeRA is the 2026 general form of the executable-template idea, and its validation chain parallels
  EngTrace's Layer 0 to 2 certification. EngTrace adds code-derived intermediate milestones and process scoring in
  engineering.

**`xu2025ugmathbench`: UGMathBench** (Xu et al.; arXiv 2501.13766; ICLR 2025). **should-cite.**

- **Construction.** 5,062 undergraduate problems from a university online homework system, covering 16 subjects, 111
  topics and 10 answer types. Each problem has three randomised versions that differ in numbers or variables.
- **Metrics.** Effective accuracy (EAcc: correct on all versions) and the reasoning gap Δ (mean accuracy minus EAcc).
- **Finding.** Across 23 LLMs the best EAcc is 56.3% (OpenAI-o1-mini), with large Δ for every model.
- **Scoring.** Final answers only (I did not check whether worked solutions ship with the data).
- **Bearing.** This is the closest dynamic benchmark at EngTrace's undergraduate difficulty. EAcc and Δ are standard
  ways to report robustness across instances of the same template.

**`zhu2023dyval`: DyVal** (Zhu et al.; arXiv 2309.17167; ICLR 2024 spotlight). **should-cite.**

- **Construction.** Tree and general DAGs generate arithmetic, linear-equation, logic and algorithmic problems with
  controllable depth and width, to resist contamination. The non-leaf nodes are by construction the intermediate steps
  of inference.
- **Evaluation and finding.** Only the root (final) answer is scored automatically. Accuracy falls with complexity, and
  a manual analysis of 20 failures per task shows partial calculation errors.
- **Bearing.** It is an early generated benchmark whose intermediate values are known by construction. EngTrace checks
  such values automatically through milestone coverage.

**`zhou2025gsminfinite`: GSM-Infinite** (Zhou et al.; arXiv 2502.05252; ICML 2025, PMLR 267). **optional.**

- **Construction.** A fully synthetic grade-school problem generator built on computational graphs abstracted from
  GSM8K. Unnecessary nodes and edges add noise, giving unbounded reasoning complexity and context length.
- **Evaluation and findings.** Final answers are scored. Performance falls along a sigmoid as complexity grows, and
  exponentially more inference compute yields only linear gains.
- **Bearing.** A generator precedent with explicit intermediate quantities. It is optional because it is mainly a
  long-context study.

**`sun2025omega`: OMEGA** (Sun et al.; arXiv 2506.18880; NeurIPS 2025). **optional.**

- **Construction.** 40 templated problem generators across six math domains produce controlled train/test pairs. Answers
  are verified symbolically, numerically or graphically.
- **Axes.** Exploratory, compositional and transformative generalisation.
- **Findings.** Accuracy degrades sharply with complexity. RL fine-tuning helps the exploratory axis, while the
  compositional axis remains limited and the transformative axis improves little or not at all.
- **Scoring.** Final answers only; generators emit verified answers, and I found no step-level gold.
- **Bearing.** A recent templated-generator benchmark to cite with GSM-Symbolic when listing template-based math work.

**`xu2025reimagine`: RE-IMAGINE** (Xu et al.; arXiv 2506.15455; ICML 2025). **should-cite.**

- **Construction.** It translates problems into an intermediate symbolic (code) representation, mutates them at three
  levels of Pearl's ladder (association, intervention, counterfactual) and renders them back.
- **Data and finding.** Applied to GSM8K, CLadder, CRUXEval and Loop, it generates arbitrarily many variants, and
  performance drops on them.
- **Bearing.** It belongs to the symbolic-synthesis family EngTrace is in. Its framing of variants that cannot be solved
  by recall alone is useful for the parameter-variation argument.

### 2.3 Perturbation and robustness evidence

**`huang2025mathperturb`: MATH-Perturb** (Huang et al.; arXiv 2502.06453; ICML 2025). **must-cite.**

- **Construction.** MATH-P-Simple and MATH-P-Hard each hold 279 perturbations of level-5 MATH problems, written by
  annotators. Hard perturbations change the problem so the original solution path no longer applies.
- **Findings.** Accuracy falls sharply on MATH-P-Hard (o1-mini −16.49%, gemini-2.0-flash-thinking −12.9%). Models also
  apply memorised techniques where they no longer fit, and putting the original problems in context makes this worse.
- **Bearing.** EngTrace's parameter draws are "simple perturbations": the solution path is kept, and per the brief, 54 to
  58 templates follow one reasoning path in all 15 instances. The paper should present them as a test of numeric and
  structural robustness and contamination control, not of transfer to new solution paths.

**`li2024gsmplus`: GSM-Plus** (Li et al.; arXiv 2402.19255; ACL 2024). **should-cite.**

- **Construction.** Eight variants of each of the 1,319 GSM8K test questions (10,552 in all). The types are numerical
  substitution, digit expansion, integer-decimal-fraction conversion, adding and reversing operations, problem
  understanding (rephrasing), distractor insertion and critical thinking.
- **Finding.** Across 25 LLMs and 4 prompting methods, models fail on questions they had solved once statements are
  added or targets changed.
- **Bearing.** The standard perturbation suite. Its rephrasing type is the ancestor of EngTrace's paraphrase arm.

**`shrestha2025gsmranges`: GSM-Ranges** (Shrestha, Kim, Ross; arXiv 2502.08680). **should-cite.**

- **Construction.** It replaces GSM8K numbers at six levels of scale.
- **Evaluation.** A grader converts each model output to code and re-runs it. If the result matches the ground truth,
  the error was arithmetic (non-logical); otherwise it was logical.
- **Findings.** Logical error rates rise by up to 14 points as numbers grow. Models do well on isolated arithmetic but
  worse when the same arithmetic is embedded in word problems.
- **Bearing.** This is the nearest published analogue of EngTrace's per-claim digit rule (E4), which separates
  arithmetic slips from reasoning errors. It also shows that number magnitude is a confound in parameterised templates.

**`mondorf2026lpds`: LPDS** (Mondorf et al.; arXiv 2605.15393, May 2026). **should-cite.**

- **Method.** Logic-preserving difficulty scaling scores how hard each allowable variation is (names, numbers, context)
  and searches the variation space for the hardest ones.
- **Findings.** Testing a random subset of variations can overstate robustness: searched variations cause drops up to
  5 times larger. Reasoning-chain errors grow with difficulty, and fine-tuning on harder variations gives more consistent
  robustness gains.
- **Bearing.** EngTrace samples its parameters at random from a private seed. LPDS is the caveat to state: random
  instances estimate typical rather than worst-case robustness.

### 2.4 Paraphrase sensitivity

**`lunardi2025paraphraserobustness`: On Robustness and Reliability of Benchmark-Based Evaluation of LLMs** (Lunardi et
al.; arXiv 2509.04013; ECAI 2025). **should-cite.**

- **Method.** Each question in six multiple-choice benchmarks (ARC-C, HellaSwag, MMLU, OpenBookQA, RACE and SciQ) gets
  five automatic paraphrases. All 34 LLMs answer every version, zero-shot.
- **Finding.** Rankings stay relatively stable, but absolute accuracy drops significantly.
- **Bearing.** This is the main precedent for EngTrace's paraphrase arm and its rank-stability reading. Note the
  differences: these benchmarks are public and multiple-choice, while EngTrace's items are private open numeric problems.
  A drop on public items is evidence of memorised wording. EngTrace's null result on never-published items fits the
  brief's wording, "no evidence of wording-level memorisation".

**`xiao2026contaminationrankings`: Contamination Inflates Scores but Rarely Reorders LLM Leaderboards** (Xiao and Cheng;
arXiv 2609.02899, v1 Jul 2026). **should-cite.**

- **Method.** Contamination is measured as the item-level accuracy contrast between original and semantically
  equivalent paraphrased items, an "anchor-item invariance" test. The data are per-instance outputs released by
  Dekoninck et al. (2024) for ARC, GSM8K, HellaSwag and MMLU. The measure is calibrated on 74 fine-tuned models with a
  known dose of contamination (+0.187 for test-set leakage; −0.012 for the negative control).
- **Findings.** On 47 public models, the rank correlation between the standard and the paraphrase-controlled leaderboard
  is 0.997. Only 3 of 188 model-benchmark cases show corroborated differential contamination.
- **Bearing.** It gives a published methodology and vocabulary for EngTrace's paraphrase arm as an invariance test. It
  also recommends rankings with confidence intervals, which matches EngTrace's τ-against-noise-floor reporting. Caveat:
  arXiv only, two authors.

**`singh2026gsmsem`: GSM-SEM** (Singh et al., Oracle AI; arXiv 2605.07053, May 2026). **should-cite.**

- **Method.** A stochastic generate, validate and filter framework that alters entities, attributes and relations while
  keeping the original calculation and answer, so each run yields fresh semantic variants. It is applied to GSM8K,
  GSM-Symbolic and GSM-Plus, and the three variant sets are released fully human-validated.
- **Finding.** Across 14 models, performance drops consistently, by about 28% on average when semantic perturbation is
  layered on GSM-Symbolic or GSM-Plus at maximum strictness.
- **Bearing.** It contrasts with EngTrace's surface-paraphrase arm. GSM-SEM argues that surface rewording mostly
  preserves facts and is the weaker test, which helps scope what EngTrace's paraphrase result can show.

### 2.5 Contamination evidence

**`zhang2024gsm1k`: GSM1k** (Zhang et al.; arXiv 2405.00332; NeurIPS 2024 Datasets and Benchmarks). **must-cite.**

- **Construction.** 1,205 newly commissioned problems that match GSM8K on human solve rate, number of steps and answer
  magnitude. The set is withheld from public release.
- **Findings.** Accuracy drops by up to 8%, and several model families overfit systematically. A model's probability of
  generating GSM8K items relates to its GSM8K-to-GSM1k gap (Spearman r² = 0.36). Frontier models show minimal
  overfitting.
- **Bearing.** The canonical evidence for contamination, and a precedent for a privately held test set. Cite it with its
  frontier caveat so as not to overstate the threat.

**`balunovic2025matharena`: MathArena** (Balunović et al.; arXiv 2505.23281; NeurIPS 2025 Datasets and Benchmarks).
**should-cite.**

- **Method.** Models are evaluated on competition problems as soon as they are released (for example AIME, CMIMC and
  IMO 2025). The paper reports "strong signs of contamination in AIME 2024".
- **Scale.** Over 50 models on 162 problems from seven competitions (v3), plus proof grading.
- **Bearing.** It is the time-window alternative to EngTrace's generated, never-published instances, and it is direct
  evidence that a famous static set is contaminated.

**`wu2025reasoningmemorization`: Reasoning or Memorization?** (Wu et al.; arXiv 2507.10532; AAAI 2026).
**should-cite.**

- **Finding on rewards.** Gains from random or incorrect rewards in Qwen2.5 RL appear on MATH-500, AMC and AIME, but not
  for Llama.
- **Contamination probe.** Given the first 60% of each MATH-500 problem, Qwen2.5-Math-7B regenerates the remaining 40%
  exactly in 54.60% of cases, against 3.8% for Llama3.1-8B. On LiveMathBench (version 202505) Qwen's completion rate is
  0.0%.
- **Clean generator.** Their RandomCalculation generator, which creates fresh arithmetic problems free of leakage, shows
  that only correct rewards help.
- **Bearing.** Recent and widely cited evidence that conclusions drawn on public sets can be artefacts of contamination,
  and that generated instances fix this.

**`schaeffer2026contaminationgenerative`: Quantifying the Effect of Test Set Contamination on Generative Evaluations**
(Schaeffer et al.; arXiv 2601.04301, v1 Jan 2026, v3 Sep 2026). **optional.**

- **Method.** Language models are pretrained on web data mixed with MATH test-set replicas, varying model size and the
  number of replicas.
- **Findings.** Even one replica pushes loss below the irreducible error of the clean corpus. Higher sampling
  temperature reduces the effect, and longer solutions are harder to memorise than short ones.
- **Bearing.** Controlled, quantitative contamination evidence for generative evaluation, where EngTrace operates.

**`sun2025bdcmitigation`: The Emperor's New Clothes in Benchmarking?** (Sun et al.; arXiv 2503.16402; ICML 2025).
**should-cite.**

- **Method.** A controlled pipeline with two question-level metrics, fidelity and contamination resistance. It tests
  10 LLMs, 5 benchmarks, 20 mitigation strategies that rewrite or regenerate items, and 2 contamination scenarios.
- **Finding.** No strategy significantly beats the unmodified benchmark on resistance across all benchmarks, and none
  balances fidelity with resistance.
- **Bearing.** It supports presenting the paraphrase arm as a probe of wording sensitivity rather than as decontamination.
  The private pool is what handles contamination.

### 2.6 Statistical re-analysis, surveys, saturation and terminology

**`dlugosz2026gsmsymbolicreeval`: The Importance of Being Statistically Earnest** (Długosz, Oliveira, Díaz-Rodríguez;
arXiv 2605.28700, v1 May 2026, v3 Sep 2026; accepted at EMNLP 2026 main, Resources and Evaluation track, per arXiv
comments). **must-cite.**

- **Method.** It re-analyses 20 of GSM-Symbolic's open-weight models with bootstrapped generalised linear mixed models
  that include per-question random effects.
- **Findings** (v3 figures). Only 8 models show significant changes under the original prompt format. The GSM-Symbolic variants have
  integers shifted towards larger values than GSM8K's (K-S statistic 0.12, p < 0.001), and controlling for this accounts
  for significance in half of the remaining cases. The failures it does find are model-specific: variable binding,
  arithmetic and dual-task interference.
- **Bearing.** GSM-Symbolic is EngTrace's lead precedent. This paper backs EngTrace's template-clustered intervals and
  Holm-corrected tests, and it argues against blanket "lack of reasoning" claims. Number-magnitude confounds (see also
  GSM-Ranges) are worth checking across EngTrace's tiers.

**`chen2025staticdynamicsurvey`: Benchmarking LLMs Under Data Contamination: A Survey from Static to Dynamic Evaluation**
(Chen et al.; arXiv 2502.17521; EMNLP 2025, pp. 10080-10098). **must-cite.**

- **Content.** It surveys static mitigations and dynamic benchmarking: temporal cutoff, rule-based generation (including
  template-based, where it files GSM-Symbolic), LLM-based rewriting and hybrids. It proposes six criteria for dynamic
  benchmarks: correctness, scalability, collision, stable of complexity (sic), diversity and interpretability.
- **Title.** The arXiv abstract page still carries the v1 title, "Recent Advances in Large Langauge Model Benchmarks
  against Data Contamination: From Static to Dynamic Evaluation". Cite the EMNLP title.
- **Bearing.** It gives the taxonomy for placing EngTrace (rule-based, template-based, private seed). Its "collision"
  criterion matches the brief's finding that 38 templates have a single question skeleton across their 15 instances.

**`akhtar2026benchmarksaturation`: When AI Benchmarks Plateau** (Akhtar et al.; arXiv 2602.16763, v4 Aug 2026;
ICML 2026). **must-cite.**

- **Method.** It defines a saturation index and analyses 60 LLM benchmarks from developer technical reports against 14
  properties.
- **Findings.** Nearly half the benchmarks show saturation, and saturation rises with benchmark age. Once age is
  considered, the 4 private and 56 public benchmarks do not differ. Nor do the 14 templated and 46 non-templated
  benchmarks (p = 0.10), or closed and open-ended formats. Expert-curated benchmarks saturate less, as do larger test sets.
- **Bearing.** This is the LLM-era replacement for, or companion to, Ott et al. (2022). It cuts against any claim that
  templates or a private pool by themselves prevent saturation. EngTrace should claim contamination control and the
  ability to refresh (new draws, harder regimes), and point to expert certification and test-set size.

**`ansari2026physicsregrading`: How Good Are Frontier Models at Physics?** (Ansari et al.; arXiv 2609.13009,
v1 Sep 2026). **should-cite.**

- **Method.** Faculty and graduate researchers audit model responses on six physics benchmarks: HLE-Physics,
  CMT-Benchmark, CritPt, UGPhysics, PRISM-Physics and PHYBench.
- **Findings.** Most audited cases first graded wrong were grader errors, wrong reference solutions or ill-posed
  questions. After correction, GPT-5.6-Sol's mean@4 rises from 47.3% to 78.7% on HLE-Physics and from 61.0% to 87.2% on
  CMT-Benchmark, and its corrected pass@4 reaches 94.4% on 54 retained CritPt challenges.
- **Bearing.** EngTrace cites UGPhysics and PHYBench. This audit shows that static technical benchmarks are confounded by
  reference and grader errors and are close to saturation. That motivates code-derived gold (per the state brief, all 2,250
  gold answers pass EngTrace's deterministic answer check) and expert-validated grading. It is contemporaneous (September 2026).

**`allawati2026contaminationresistant`: LLM Benchmark Datasets Should Be Contamination-Resistant** (Al-Lawati et al.;
arXiv 2605.19999; ICML 2026 Position Paper Track). **optional.**

- **Content.** It defines a contamination-resistant dataset as one that stays usable for inference but is unlearnable
  (Definition 2.1), and grounds this in the asymmetry between Transformer training and inference.
- **Bearing.** The May paper calls EngTrace's instances "contamination-resistant" in the sense that they were private
  when the models ran. With an ICML 2026 position paper fixing a different meaning, EngTrace should define the term or
  use "private-pool" or "contamination-controlled".

### 2.7 Milestone-level diagnosis

**`wang2026milestoneoracles`: Solving Every Step Is Not Enough: Milestone Oracles** (Wang et al.; arXiv 2609.32235,
v1 26 Sep 2026; accepted at NeurIPS 2026 Evaluations and Datasets Track, per arXiv comments). **should-cite.**

- **Method.** In OracleLadder, a teacher model writes a fixed roadmap of milestones, and a deterministic symbolic verifier
  grades every answer. Comparing no help, the roadmap, the roadmap with milestone answers, and each milestone alone sorts
  every failure into one of five gaps.
- **Data.** 354 NuminaMath problems and six models from 8B to 671B parameters.
- **Findings.** The largest gap for every model is the composition gap: 33-48% of problems, or 24-37% after removing
  suspected grading errors. Accuracy and milestone-help recovery rank the two strongest models differently. The effect
  replicates on MATH500 and AIME 2024/25, and per-problem recovery agrees 83-87% under a second teacher.
- **Bearing.** It uses the same word as EngTrace's "milestone coverage" for a different construct. Theirs are
  teacher-written oracle hints; EngTrace's are code-derived values checked for presence in the model's own trace.
  The two should be told apart explicitly. It is contemporaneous (five days before this search), so cite it for
  completeness, not as a required comparison.

---

## 3. Considered and not included (one line each)

- **VarBench** (Qian et al., Findings EMNLP 2024, 2406.17681): verified; GSM-Symbolic and UGMathBench make its
  variable-perturbation point; use only if space allows.
- **DyVal 2 / Meta Probing Agents** (Zhu et al., ICML 2024, 2402.14865): LLM-agent rewriting of items; less relevant than
  DyVal.
- **LiveCodeBench** (Jain et al., ICLR 2025, 2403.07974): verified; a code benchmark; LiveBench and MathArena already
  represent the time-window approach.
- **MathCheck** (Zhou et al., ICLR 2025, 2407.08733): a checklist with process-judging tasks and robustness variants;
  better placed with the process-supervision agent's papers.
- **DynaMath** (Zou et al., ICLR 2025, 2411.00836): program-generated visual math variants; multimodal.
- **Reasoning Gym** (Stojanovski et al., NeurIPS 2025 D&B, 2505.24760): RLVR procedural generators with verifiers but no
  gold traces; training-oriented.
- **Physics of Language Models 2.1 / iGSM** (Ye et al., ICLR 2025, 2407.20311): a synthetic GSM generator with dependency
  graphs, used to train small models.
- **Correct Answers, Invalid Traces** (Puduppully et al., 2609.38107, 29 Sep 2026): on iGSM, 31.6% of correct answers
  on the hardest instances have invalid traces, but the models are trained on iGSM. Only two days old.
- **VAR-MATH** (Yao et al., 2507.12885) and **EFAGen / Executable Functional Abstractions** (Khan et al., 2504.09763):
  symbolic multi-instance and inferred generator programs for math; VeRA covers both ideas.
- **V-FiLLM** (Larsen et al., 2608.11047, Aug 2026): finance items from executable computation trees with verified CoT
  traces; dropped for space as a third finance analogue. It was downloaded in the first fetch run and its orphaned
  PDF and text were then deleted.
- **MatheMagic** (O'Brien et al., 2510.05962): counterfactual arithmetic generation; off EngTrace's axis.
- **LINGOLY-TOO** (Khouja et al., 2503.02972): templatised orthographic obfuscation for linguistics puzzles; out of domain.
- **GSM-Identity** (Negi, Puccetti, Esuli, *Machine Learning* 2026): equivalence transformations; no arXiv or open PDF
  found (see section 5).
- **GSM-Plus-BN** (2607.13248), **MGSM-Pro** (2601.21225) and **Lost in Cultural Translation** (2503.18018): language or
  culture extensions of GSM perturbation suites.
- **Robust Reasoning Benchmark** (Golikov et al., 2604.08571), **mathematically-equivalent transformations** (Hao et al.,
  2508.08833), **Toward Automated Robustness Evaluation of Mathematical Reasoning** (2506.05038, Findings ACL 2026),
  **RIDE** (2511.04120), **Numerical Sensitivity and Robustness** (2511.08022) and **numeric-remapping attacks** (Barker
  et al., 2606.03606): further perturbation studies on math; MATH-Perturb, GSM-Plus, LPDS and GSM-Ranges cover the point.
- **MathRobust-LV** (2510.06430), **Surface-Form Brittleness Under Paraphrase Stress Tests** (Navarro Carranza,
  2510.08616, NeurIPS 2025 workshop; 2 models on ARC) and **RoParQ** (2511.21568, training): paraphrase studies; Lunardi
  et al. and Xiao and Cheng are larger and cover the point.
- **DVD, variant contamination detection** (Liang et al., 2601.04895), **Zero-CoT Probe** (Lan et al., 2605.21856),
  **Controllable Contamination Detection with Statistical Guarantees** (ACL 2026), **When Benchmarks Leak** (ACL 2026,
  2601.19334), **Test of Time** (ACL 2026, 2509.00072) and **LeakScale** (2609.27176): detection and mitigation methods;
  EngTrace does not run detection.
- **Benchmark Contamination: A Taxonomy Organized by Defeated Mitigation** (Angulo et al., 2608.29463): useful framing
  ("a private test set addresses only the first type", i.e. direct contamination); low visibility; optional if the
  contamination paragraph grows.
- **A Survey on Data Contamination for LLMs** (Cheng et al., 2502.14425), **Benchmark Data Contamination of LLMs: A
  Survey** (Xu et al., 2406.04244), **Are LLM Benchmarks Already Contaminated?** (GEM 2026) and **Evaluation data
  contamination in LLMs** (Singh et al., 2411.03923): alternatives to the chosen EMNLP 2025 survey.
- **AntiLeakBench** (ACL 2025, 2412.13670), **LastingBench** (EMNLP 2025), **CapBencher** (2505.18102), **KUMO /
  Generative Evaluation** (2504.02810), **BeyondBench** (2509.24210), **LiveAoPSBench** (2501.14275) and **Mathador-LM**
  (2406.12572): other contamination-resistant or live constructions; not needed beyond LiveBench and MathArena.
- **Benchmarks Saturate When The Model Gets Smarter Than The Judge** (Ballon et al., 2601.19532): judge errors mask model
  differences on Omni-MATH; relevant to the judge-validity paragraph rather than to this topic.
- **The Ouroboros of Benchmarking** (2511.01365), **Measuring Five-Nines Reliability** (2605.11209) and **The Growing
  Pains of Frontier Models** (2605.18840, not opened): saturation essays and methods; Akhtar et al. is the systematic
  study.
- **Right Answer, Wrong Method** (Ren et al., 2608.02442): 8.2-44.1% of credited answers are "hacked"; strong process
  motivation, but in the process-supervision agent's scope.
- **PRISM-Physics** (2510.03185, ICLR 2026) and **Beyond the Answer Key** (2608.28725, ICMLA 2026): process evaluation in
  physics, and evaluator robustness to non-canonical traces; process-supervision scope. The second bears on whether
  milestone coverage penalises valid alternative paths.
- **PRIME** (2602.11570, ACL 2026), **ChemCoTBench-V2** (2606.03660) and **OPT-Engine** (2601.19924, ICML 2026): relevant,
  but already in `engineering_web.json` and/or `process_supervision.json`.
- **WirelessMathBench-XL** (2509.23219, NeurIPS 2026 E&D): a 13-gram contamination audit of a wireless-engineering
  benchmark; for the engineering agent.
- **SciR** (2606.13020, NeurIPS 2026 E&D), **Physics literacy in parallel physical worlds** (2607.00276), **Macaron**
  (ACL 2026), **PuzzleClone** (Findings ACL 2026), **ReTRE** (ACL 2026), **Think Like You Execute** (ACL 2026 Industry),
  **Loong** (2509.03059), **Reasoning Core** (2509.18083) and **SMDD-Bench / TSBench / oMeBench**: synthetic or verifiable
  data in other formats or domains, or training-data work.
- **"Revisiting GSM-Symbolic"** (LessWrong post, March 2026): reports that GSM-Symbolic's effects shrink once ambiguous
  items are audited out; non-archival, so it is context only.

---

## 4. Updates to existing citations

- **GSM-Symbolic** (`mirzadeh2024gsmsymbolic`).
  - **Venue.** ICLR 2025. The arXiv v2 (27 Aug 2025) comment reads "ICLR camera ready + additional discussion in the
    appendix", and Semantic Scholar lists it with DBLP key `conf/iclr/MirzadehASTBF25`. The May reference cites it as
    "arXiv preprint arXiv:2410.05229" (2024).
  - **Author names are wrong in the May reference list.** It prints "Kiarash Alizadeh" and "Hidetoshi Shahrokhi". The
    arXiv page lists Iman Mirzadeh, Keivan Alizadeh, Hooman Shahrokhi, Oncel Tuzel, Samy Bengio and Mehrdad Farajtabar.
  - **Context.** Długosz et al. (2026, above) re-analyses its statistics.
- **Putnam-AXIOM** (`gulati2024putnamaxiom`).
  - **Venue and title.** The May paper cites the NeurIPS 2024 MATH-AI workshop version. The archival version is arXiv
    2508.08292 (v1 5 Aug 2025, v2 27 Aug 2025), journal-ref "ICML 2025", DBLP key `conf/icml/GulatiMCXFDK25`. Its title
    adds "in LLMs": "Putnam-AXIOM: A Functional and Static Benchmark for Measuring Higher Level Mathematical Reasoning in
    LLMs".
  - **Content.** 522 original problems and 100 functional variants. o1-preview scores 41.9% on the originals and drops by
    19.6 points (46.8% relative) on the variants, and the paper proposes Teacher-Forced Accuracy.
  - **Action.** Set `arxiv` to 2508.08292 in papers.json and re-check any Putnam-AXIOM numbers the paper quotes. The May
    text also writes "Putnam-Axiom"; the paper's own spelling is "Putnam-AXIOM".
- **FinChain** (`xie2025finchain`).
  - **Venue.** ACL 2026, Long Papers, Anthology ID `2026.acl-long.662`, read from the ACL 2026 event page in the ACL
    Anthology, which lists the same 25 authors. The arXiv id is 2506.02515 (v1 3 Jun 2025, v4 30 Apr 2026).
  - **Published abstract.** 58 topics across 12 financial domains, parameterised symbolic templates with executable
    Python, ChainEval, and 26 LLMs.
  - **Action.** The May reference is "2025 ... arXiv preprint" with no id; update it to ACL 2026 and set `arxiv` to
    2506.02515. FinChain's GitHub README still says "EMNLP 2025 submission", which is out of date.
- **Deng et al. 2024** (`deng2024contamination`). Verified: arXiv 2311.09783, whose v2 (3 Apr 2024) is marked "NAACL 2024
  Version". The authors (Chunyuan Deng, Yilun Zhao, Xiangru Tang, Mark Gerstein, Arman Cohan) match the May reference.
  No change needed.
- **Ott et al. 2022** (`ott2022saturation`). Verified: Nature Communications 13, article 6793 (2022), per the arXiv
  journal-ref. It has an arXiv version, 2203.04592 (v4 7 Oct 2022), which can go into papers.json for fetching. Pair it
  with Akhtar et al. (ICML 2026) for current LLM saturation evidence.
- **Incidental, outside this topic.** EngiBench (Zhou et al., `zhou2025engibench`) appears in Findings of ACL 2026 as
  `2026.findings-acl.1810` (ACL 2026 event page). The May reference cites arXiv.
- **Observation.** EngTrace's own preprint, arXiv 2511.01650, is public and surfaces in searches for symbolic
  engineering benchmarks. v1 was titled "EngChain: A Symbolic Benchmark for Verifiable Multi-Step Reasoning in
  Engineering", and v3 is "EngTrace: ...", describing 90 templates, three branches, 1,350 items and 27 LLMs. Reviewers
  may find it, and its figures are the May ones.

---

## 5. What I could not verify

- **GSM-Identity** (K. Negi, G. Puccetti, A. Esuli, "GSM-Identity: Evaluating Mathematical Reasoning in LLMs via
  Equivalence Transformations", *Machine Learning* 115(4), 2026). I saw it only in a search-result snippet of the
  authors' publication page and in Semantic Scholar's citer lists (venue "Machine-mediated learning", 2026-03-10). I
  found no arXiv version or open PDF, so it is not included.
- **Venues from arXiv comments only, not checked against proceedings:**
  - Długosz et al.: EMNLP 2026; its proceedings are not yet out.
  - Wang et al.: NeurIPS 2026 Evaluations and Datasets Track.
  - Akhtar et al.: ICML 2026.
  - Al-Lawati et al.: ICML 2026 Position Paper Track.
  - Wu et al.: AAAI 2026.
  - Lunardi et al.: ECAI 2025.
  - TemplateGSM: ICLR 2025 DATA-FM Workshop.
  - GSM1k: NeurIPS 2024 Datasets and Benchmarks.
- **Venues from secondary metadata:**
  - OMEGA, NeurIPS 2025: from Semantic Scholar and DBLP metadata (`conf/nips/SunHZZHDS25`); track not confirmed.
  - MathArena, NeurIPS 2025 Datasets and Benchmarks: from the proceedings URL in a search result.
  - GSM-Infinite, ICML 2025: confirmed on its PMLR page, where the title reads "GSM-∞".
  - SymPyBench: confirmed on its ACL Anthology page.
- **Abstracts paraphrased by WebFetch.** For about a dozen papers, WebFetch returned paraphrased rather than verbatim
  abstracts (among them DyVal, UGMathBench, LiveBench, Lunardi, Reasoning or Memorization, GSM-Ranges, Schaeffer,
  Al-Lawati). I re-checked in the downloaded full texts the numbers used above for DyVal, UGMathBench, Lunardi, Wu et al.,
  GSM-Ranges and Al-Lawati. Schaeffer et al.'s findings and LiveBench's "below 70%" come from the abstract summaries only.
- **Sci-Rho's step judge.** My claim that it reports no human validation comes from a keyword search of the full text
  ("agreement", "human evaluat", "kappa", "validat"), not a full reading.
- **Failed arXiv API queries.** The arXiv API returned HTTP 429 repeatedly, because several agents share the IP.
  `all:"symbolic" AND all:benchmark AND all:contamination` failed in all three attempts, and `abs:"intermediate values"
  AND abs:benchmark AND abs:LLM` failed once before succeeding. The arXiv web search for "symbolic templates benchmark"
  covered similar ground.
- **Web search budget.** The session's shared cap of 200 web searches was reached after my 36th query. Four further
  queries were refused and are listed in section 6.1. Their targets were verified directly on arXiv abstract pages.
- **Semantic Scholar citers of SymPyBench.** It returned 0, so I could not citation-chase SymPyBench.
- **DBLP** pages were blocked by a bot check, so no venue rests on DBLP directly; DBLP keys were read via Semantic
  Scholar.
- **The fetch script's manifest.** `MANIFEST.json` is loaded at the start of each `fetch_papers.py` run and rewritten
  whole at the end. Concurrent runs by several agents therefore drop one another's records. My second run found my
  earlier entries missing from the manifest and re-extracted and re-recorded them ("present ... text extracted"). At one check there were 129 PDFs on disk but 62 manifest records. A final
  `python docs/related_work_oct2026/fetch_papers.py` after every agent has finished re-records every PDF on disk. This
  is a property of the script, which I did not edit, not an error message.

---

## 6. Queries (verbatim)

### 6.1 WebSearch (36 run, in order)

1. "GSM-Symbolic" follow-up benchmark 2026
2. "FinChain" symbolic benchmark verifiable chain-of-thought financial reasoning
3. Putnam-AXIOM functional variations arXiv 2025
4. symbolic templates benchmark engineering problems LLM verifiable reasoning traces 2026
5. benchmark data contamination survey LLM 2025 2026 arXiv
6. benchmark saturation evidence large language models 2026 arXiv
7. "GSM-Identity" Negi 2026 arXiv
8. MATH-Perturb benchmarking LLMs math reasoning hard perturbations
9. GSM-Infinite infinitely scalable increasing complexity long context reasoning benchmark
10. VarBench dynamic variable perturbation benchmark contamination
11. UGMathBench dynamic undergraduate mathematics benchmark randomized versions reasoning gap
12. RusFinChain Russian benchmark verifiable chain-of-thought
13. DyVal dynamic evaluation large language models reasoning tasks directed acyclic graphs; DyVal 2 meta probing agents
14. "Functional Benchmarks for Robust Evaluation of Reasoning Performance" reasoning gap MATH() Srivastava
15. TemplateGSM TemplateMath template-based data generation GPT-4 meta templates code solutions
16. "A careful examination of large language model performance on grade school arithmetic" GSM1k overfitting
17. MathCheck checklist evaluating mathematical reasoning robustness task generalization ICLR 2025
18. LiveBench challenging contamination-limited LLM benchmark ICLR 2025 objective ground truth
19. LiveCodeBench holistic contamination free evaluation code time window
20. GSM-Plus comprehensive benchmark robustness math word problem solvers perturbations ACL 2024
21. paraphrase robustness math word problems LLM memorization wording 2025 arXiv
22. template-based benchmark physics problems parameterized python generated instances step-by-step solutions LLM 2026
23. procedurally generated chemistry benchmark LLM contamination-free verifiable intermediate steps 2025 2026
24. "GSM-Ranges" logical and arithmetic errors wide numerical ranges
25. "Reasoning or Memorization" unreliable results reinforcement learning data contamination RandomCalculation Qwen2.5
26. MathArena evaluating LLMs on uncontaminated math competitions contamination evidence AIME
27. SymPyBench dynamic benchmark scientific reasoning executable Python code physics parameterized
28. "From Answers to States" verifiable process-level evaluation chemical reasoning ChemCoTBench
29. OMEGA can LLMs reason outside the box math programmatic problem generators exploratory compositional transformative
    generalization
30. Reasoning Gym procedural generators verifiable rewards reinforcement learning NeurIPS 2025
31. engineering benchmark parameterized problem generator LLM contamination thermodynamics circuits 2026 arXiv
32. MatheMagic generating dynamic mathematics benchmarks robust to memorization
33. "The Importance of Being Statistically Earnest" GSM-Symbolic re-evaluation
34. GSM-SEM semantically variant augmentations benchmark framework 2026
35. "On Robustness and Reliability of Benchmark-Based Evaluation of LLMs" paraphrased benchmark questions
36. "Robust Reasoning Benchmark" arXiv 2604.08571

Refused because the session's shared web-search cap was reached (targets then verified via arXiv abstract pages):
"LLMs Show Surface-Form Brittleness Under Paraphrase Stress Tests"; "Benchmarking Large Language Models Under Data
Contamination" survey static to dynamic evaluation Chen EMNLP 2025; "When AI Benchmarks Plateau" systematic study
benchmark saturation 60 benchmarks; "LLM Benchmark Datasets Should Be Contamination-Resistant" position.

### 6.2 arXiv API (export.arxiv.org, sorted by submission date, 100 results, 2024-2026 kept; 3-4 s pauses)

URL pattern: `https://export.arxiv.org/api/query?search_query=<query, URL-encoded>&start=0&max_results=100&sortBy=submittedDate&sortOrder=descending`.

First batch (stopped after three queries by a console-encoding error in my parser, then rerun):

- `all:"symbolic" AND all:benchmark AND all:contamination`: failed (HTTP 429).
- `abs:"symbolic templates" AND abs:benchmark`: 3 hits (FinChain, GSM-Symbolic, one pre-2024).
- `abs:"templated" AND abs:benchmark AND abs:LLM`: 100 hits, almost all off-topic.

Second batch:

1. `all:"symbolic" AND all:benchmark AND all:contamination`: failed (HTTP 429).
2. `abs:"dynamic benchmark" AND abs:contamination`: 20 hits.
3. `abs:"functional variants"`: 43 hits (only MATH() and Putnam-AXIOM relevant).
4. `ti:"data contamination"`: 81 hits.
5. `abs:perturbation AND abs:robustness AND abs:"mathematical reasoning"`: 26 hits.
6. `abs:"verifiable chain-of-thought"`: 6 hits.
7. `abs:parameterized AND abs:templates AND abs:benchmark AND abs:reasoning`: 9 hits.
8. `abs:"contamination-resistant"`: 42 hits.
9. `abs:"contamination-free" AND abs:benchmark AND abs:reasoning`: 25 hits.
10. `ti:saturation AND ti:benchmark`: 20 hits.
11. `abs:paraphrase AND abs:benchmark AND abs:memorization`: 19 hits.
12. `abs:"procedurally generated" AND abs:benchmark AND abs:LLM`: 41 hits.
13. `abs:symbolic AND abs:benchmark AND abs:engineering AND abs:"large language models"`: 62 hits.
14. `abs:"step-level" AND abs:symbolic AND abs:benchmark`: 19 hits.
15. `abs:executable AND abs:templates AND abs:physics AND abs:benchmark`: 9 hits.
16. `abs:"intermediate values" AND abs:benchmark AND abs:LLM`: failed (HTTP 429).

Third batch:

- `all:"symbolic" AND all:benchmark AND all:contamination`: failed again (HTTP 429).
- `abs:"intermediate values" AND abs:benchmark AND abs:LLM`: 3 hits, none relevant.
- `abs:"templated benchmark"`: 4 hits, none relevant.
- `abs:"verifiable" AND abs:"symbolic templates"`: 2 hits (FinChain).

### 6.3 arXiv web search (arxiv.org/search via WebFetch, newest first)

1. `query=symbolic+templates+benchmark&searchtype=abstract`: yielded Sci-Rho, the numeric-remapping attacks paper and
   the industrial-maintenance benchmark DiagnosticIQ.
2. `query=PRIME+process-outcome+alignment+verifiable+reasoning+mathematics+engineering&searchtype=all`: resolved PRIME to
   2602.11570.
3. `query="verifiable chain-of-thought"&searchtype=all`: yielded V-FiLLM.
4. `query="dynamic benchmark" contamination&searchtype=abstract`
5. `query="functional variants" OR "functional variant" benchmark&searchtype=abstract`: malformed; it returned all recent
   papers.
6. `query=symbolic+benchmark+parameterized+templates+executable+python&searchtype=abstract`: returned RusFinChain and
   FinChain only.
7. `query=contamination+benchmark+engineering+LLM+reasoning+generated&searchtype=abstract`
8. `query=paraphrase+benchmark+contamination+memorization+LLM&searchtype=abstract`: yielded Xiao and Cheng, and DVD.
9. `query=templates+physics+chemistry+benchmark+LLM+step-by-step+generated+instances&searchtype=abstract`: no results.
10. `query=benchmark+saturation+language+models&searchtype=title`
11. `query=perturbation+robustness+mathematical+reasoning+variants+LLM&searchtype=abstract`

### 6.4 Citation chasing and venue lookups

- **Semantic Scholar citations API** (`/graph/v1/paper/arXiv:<id>/citations`, all pages), filtered to 2025-2026 titles
  with topic keywords:
  - GSM-Symbolic 2410.05229: 679 citers, 233 kept. This is where Statistically Earnest, Milestone Oracles, VeRA, LPDS,
    the taxonomy paper, VAR-MATH, RE-IMAGINE, EFAGen and others came from.
  - FinChain 2506.02515: 19 (RusFinChain, LPDS, Loong).
  - Putnam-AXIOM 2508.08292: 15.
  - MATH-Perturb 2502.06453: 101.
  - MATH() 2402.19450: 87.
  - GSM1k 2405.00332: 240.
  - Chen et al. survey 2502.17521: 39.
  - SymPyBench 2512.05954: 0.
- **Semantic Scholar batch metadata** (`/graph/v1/paper/batch`) for 40 candidate ids: venues and DBLP keys.
- **OpenAlex** (`/works/doi:10.48550/arXiv.<id>`): GSM-Symbolic had 55 citations, and FinChain was not found by DOI.
  Abandoned for Semantic Scholar.
- **ACL Anthology:** the ACL 2026 event page `https://aclanthology.org/events/acl-2026/` was downloaded and scanned for
  titles, which gave the FinChain, EngiBench and PRIME anthology IDs. Also the paper page `2026.eacl-industry.8` for
  SymPyBench and `2025.emnlp-main.511` for the Chen et al. survey.
- **PMLR:** `proceedings.mlr.press/v267/zhou25m.html` for GSM-Infinite.
- **GitHub API:** the mbzuai-nlp/finchain README, which still says "EMNLP 2025 submission".
- **DBLP:** blocked by a bot check.
- **Abstract pages:** each of the 30 included papers, plus about 20 rejected candidates, was opened on arxiv.org/abs via
  WebFetch.
