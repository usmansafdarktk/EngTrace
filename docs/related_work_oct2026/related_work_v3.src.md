<!-- part: md-header -->
# Related Work, revised for the October 2026 submission

This is the proposed replacement for section 2 of the May 2026 submission (`docs/_ARR_May__EngTrace.pdf`),
with an appendix for what the main text cannot hold. It is rendered by `render_related_work.py` from
`related_work_v3.src.md`, which also writes the LaTeX (`related_work_v3.tex`); the bibliography is
`related_work_sources.bib`, and every change from the May text, with its reason, is in `CHANGES.md`.

Every statement about another paper below is checked against that paper's text: `facts.py` lists each
claim with the passage that supports it, and `verify_facts.py` finds all of them (`notes/fact_check.md`).
The choices rest on `CITED_SET_AUDIT.md` (the May citations) and `NEW_PAPERS_AUDIT.md` (the papers found
in the October search). Citation labels below are computed from the bibliography the way natbib prints
them, so they match the LaTeX.

Length, counted by `render_related_work.py` without citations: the main text has 704 words of prose, against
227 in the May version. It adds the paragraph on process supervision and evaluator validity that the plan of
record asks for, and the engineering paragraph now says what prior benchmarks actually score. Removing the
sentence marked [cut first] brings it to 678. If more must go, cut in this order: the coding benchmarks in
the first sentence, all but three of the single-branch benchmarks, and the PRIME and ChemCoTBench-V2 sentence.
The follow-up review (`notes/followup_review.md`) cut the other [cut first] sentence, on re-graded physics
answers, and cites that paper for its point about rule-based graders in the paragraph on evaluators.

---

<!-- part: body -->
## 2 Related Work

**Reasoning benchmarks and what they score.**
Mathematical and coding benchmarks such as GSM8K [@cobbe2021gsm8k], MATH [@hendrycks2021math], HumanEval [@chen2021humaneval], MBPP [@austin2021mbpp] and SWE-bench [@jimenez2024swebench] score a final answer or a test outcome, and broad suites such as MMLU [@hendrycks2021mmlu] and SuperGPQA [@du2025supergpqa] score a chosen option.
Yet a correct final answer can follow flawed reasoning [@cobbe2021gsm8k], increasingly so on harder problems: in ProcessBench, the share of correct-answer solutions that contain a process error rises from 3.5% on GSM8K to 51.8% on Omni-MATH [@zheng2024processbench].
Most physics benchmarks likewise grade final answers or expressions [@xu2025ugphysics; @qiu2025phybench; @zhang2025abenchphysics; @wang2023scibench], while PhysReason [@zhang2025physreason], PRISM-Physics [@zhao2026prismphysics] and HiPhO [@yu2026hipho] score steps against reference solutions, with an LLM judge, rule-based formula matching, and official marking schemes applied by an LLM judge, respectively.

**Generated instances and contamination.**
Static test sets saturate [@ott2022saturation; @akhtar2026benchmarksaturation] and leak into training data [@deng2024contamination; @zhang2024gsm1k], and neither private test sets nor templating has measurably slowed saturation [@akhtar2026benchmarksaturation].
Functional, symbolic and perturbed variants of existing problems test whether a model has memorised the original instance or its method [@srivastava2024functionalbench; @mirzadeh2024gsmsymbolic; @gulati2024putnamaxiom; @huang2025mathperturb], although GSM-Symbolic's drops do not always survive re-analysis [@dlugosz2026gsmsymbolicreeval], randomly sampled variants can overstate robustness [@mondorf2026lpds], and on multiple-choice benchmarks paraphrasing lowers absolute scores while leaving model rankings relatively stable [@lunardi2025paraphraserobustness].
Generation has also been paired with gold solutions: HARDMath produces applied-mathematics problems with numerically validated step-by-step solutions [@fan2024hardmath]; FinChain pairs each templated financial problem with an executable gold reasoning chain and scores a model's steps against it [@xie2025finchain]; and SymPyBench [@imani2025sympybench], AtmosSci-Bench [@li2025atmosscibench] and Sci-ρ [@azmi2026scirho] instantiate science problems from code templates.

**Engineering benchmarks.**
General suites test engineering mostly through multiple choice; SuperGPQA's 7,892 engineering questions often require calculation but are scored on the option chosen [@du2025supergpqa].
Dedicated benchmarks typically cover one branch: transportation [@syed2024transportbench], control [@kevian2024controlbench], power dispatch [@zhou2024elecbench], electrical and electronics engineering [@li2025eeebench], circuits [@skelic2025circuit], astrodynamics [@wu2025apbench], thermodynamics [@duzkar2026thermoqa], thermal-protection design [@zheng2026tpscalcbench] or operations research [@mostajabdaveh2025orqa; @chen2026optengine].
Multi-branch benchmarks grade answers with LLM judges: EngiBench [@zhou2025engibench] with two LLM evaluators, plus rubrics for its open-ended level, and ERI [@naser2026eri] against reference answers that an LLM wrote along with the questions.
<!-- cut first -->Other work evaluates design [@guo2025engdesign], multiphysics simulation [@mudur2025feabench] and questions from industrial practice [@heesch2025realworld], the last questioning evaluations built from examination materials whose correctness is easily verifiable.
A few benchmarks look past the answer: experts examined the reasoning in TransportBench and ControlBench, LLM judges score it in ElecBench, TPS-CalcBench and EngVQA [@wasiq2026engvqa], and ThermoQA scores property values and weighted solution steps against values computed with a thermodynamic property library.
Fewer generate their instances: CIRCUIT instantiates circuit templates with several numerical setups, OPT-Engine [@chen2026optengine] generates optimisation problems of controlled complexity, and a power-systems agent benchmark draws held-out cases from privately seeded generators [@trashchenkov2026psab].
To our knowledge, none combines generated instances across engineering branches with gold traces, deterministic checks of intermediate values, and an evaluator validated against experts' step-level labels.

**Process supervision and the validity of evaluators.**
Outcome supervision can match process supervision on final-answer error, but correct reasoning steps required process-based supervision or reward models that emulate it [@uesato2022processoutcome].
Process reward models (PRMs) learn from human step labels [@lightman2023verify] or from labels estimated by rollouts [@wang2023mathshepherd], which are noisy [@zhang2025prmlessons], and are evaluated on human-annotated first-error benchmarks in mathematics [@zheng2024processbench; @song2025prmbench]; PRMs trained on mathematics perform poorly in other domains [@zeng2025versaprm].
Beyond mathematics, PRIME tests whether verifiers catch answers that the reasoning does not support, in college-level STEM including engineering [@wang2026prime], and ChemCoTBench-V2 checks chemistry steps with deterministic rules [@guo2026chemcotbenchv2].
LLM judges are the common alternative, but they favour their own outputs and those of related models [@zheng2023mtbench; @panickssery2024selfrecognition; @spiliopoulou2025playfavorites; @li2026preferenceleakage], even on objective rubric items [@pombal2026rubricselfpreference]; their errors are correlated across models [@goel2025greatmodels; @kim2025correlatederrors], which limits what a panel of judges [@verga2024poll] gains by voting; and on difficult objective comparisons many do little better than chance [@tan2025judgebench], agreeing with experts mainly on questions they can answer themselves unless given a correct reference [@krumdick2026nofreelabels].
Rule-based checkers avoid these biases but can reject correct answers given in unexpected formats [@huang2026verifierrobustness; @ansari2026physicsregrading].

**Positioning.**
EngTrace brings these lines together.
As in GSM-Symbolic and FinChain, symbolic templates generate the instances, here across five engineering branches and from a private seed, and as in FinChain the gold trace comes from the code that computes the answer.
Intermediate values are checked deterministically against that trace and arithmetic by recomputation; an LLM judge from a model family outside the evaluated roster decides only what these checks leave open, given each milestone's expected value; and the full evaluator, with off-the-shelf PRMs as baselines, is validated against step labels from 15 domain experts.
Expert-checked paraphrases test whether models exploit template wording.

<!-- part: appendix -->
## Appendix: Extended related work

**Table A.** How EngTrace relates to the closest benchmarks. "–" means not provided, or not reported in the sections reviewed for this audit.

| Benchmark | Domain | Instances | Gold steps | What is scored, and how | Scorer checked against people |
|---|---|---|---|---|---|
| FinChain [@xie2025finchain] | finance | generated from symbolic templates | executable gold chain | steps, against the gold chain | Spearman 0.655 with experts' ratings of reasoning quality |
| HARDMath [@fan2024hardmath] | applied mathematics | generated by code | step-by-step, validated numerically | final answer, graded by an LLM with rubrics | – |
| SymPyBench [@imani2025sympybench] | university physics | parameterised Python | step-by-step reasoning | answers, and consistency across variants | – |
| AtmosSci-Bench [@li2025atmosscibench] | atmospheric science | 67 templates with Python solvers, as multiple choice; 391 static open-ended questions | – | the option chosen; open-ended answers by numeric, then symbolic, then LLM checks | its LLM grader agreed with human graders on 92.79% and 93.02% of the open-ended answers it decided (two models' outputs) |
| Sci-ρ [@azmi2026scirho] | school STEM | 606 Python templates | reasoning steps | steps, as an F1 from an LLM judge | – |
| PhysReason [@zhang2025physreason] | physics | static | annotated solution steps | steps, an LLM judge locating the first error | about 98% first-error accuracy on 1,000 annotated solutions |
| PRISM-Physics [@zhao2026prismphysics] | physics | static (PhD qualifying-exam problems) | formula graphs | steps, by rule-based formula matching | Kendall τb 0.346 with two experts' scores of 70 solutions (an outcome-only LLM judge: 0.294) |
| CIRCUIT [@skelic2025circuit] | circuits | 102 templates, several numerical setups each | – | final answers, per template | human review found about 5% false positives |
| ThermoQA [@duzkar2026thermoqa] | thermodynamics | static (293 problems) | states and steps computed with CoolProp | property values and weighted steps, within 2% | – |
| TPS-CalcBench [@zheng2026tpscalcbench] | thermal-protection design | static (420 items) | – | answer, and a reasoning rubric scored by an LLM judge | weighted κ 0.68 to 0.82 per rubric dimension |
| EngVQA [@wasiq2026engvqa] | engineering, multimodal | static (696 problems) | – | eight solution stages, by an LLM judge | reported against adjusted ("synthetic") human scores |
| EngiBench [@zhou2025engibench] | several engineering branches | static | – | answers by two LLM evaluators; rubrics for the open-ended level | humans resolve evaluator disagreements and review the open-ended scores |
| ERI [@naser2026eri] | nine engineering disciplines | generated by an LLM (57,750 records) | LLM-written references | a five-point score from three LLM judges | – |
| TransportBench [@syed2024transportbench] | transportation | static | – | answers and reasoning, graded by experts | the experts are the graders |
| ChemCoTBench-V2 [@guo2026chemcotbenchv2] | chemistry | static | reference traces | intermediate steps, by deterministic rules | 87.4% agreement with adjudicated expert labels (κ 0.74; experts with each other, 90.1%) |
| **EngTrace** | five engineering branches | 150 templates, instances from a private seed | from the template code | final answer; intermediate values and arithmetic, deterministically; a residual LLM judge from outside the evaluated roster | step labels from 15 domain experts on 300 traces |

**Further related work.**
*Engineering.* Recent benchmarks extend coverage to agents working in professional engineering software [@gao2026engiworld], strength of materials [@wan2025som1k], circuit analysis from schematics with synthetically generated problems [@akbari2025circuitsense], telecommunications mathematics [@colle2025telemath], executable structural analysis [@qin2026structureclaw], finite-element programming [@mohammadzadeh2025fembench], computational fluid dynamics [@somasekharan2026cfdllmbench], civil-engineering licensure questions [@silwal2026pecivilbench] and process systems engineering judged by several LLMs [@sukpancharoen2026psebench].
*Science.* OlympiadBench pairs olympiad problems with expert step-by-step solutions but scores final answers [@he2024olympiadbench]; SciDA randomises numerical parameters at each evaluation [@zhou2025scida]; NEWTON generates physical-commonsense questions from templates [@wang2023newton].
*Contamination.* VeRA turns benchmark items into executable task families [@cheng2026vera], GSM-Plus perturbs GSM8K problems [@li2024gsmplus], @chen2025staticdynamicsurvey survey the move from static to dynamic benchmarks, and @xiao2026contaminationrankings measure contamination as the gap between original and paraphrased items.
*Step verification.* Error-detection benchmarks extend beyond mathematics [@tyen2023bigbenchmistake; @zeng2024mrben; @he2025deltabench; @sun2026lsrben], physics now has a process reward model [@dong2026physprm], and GenPRM verifies each step with code [@zhao2025genprm]; reference-free step metrics include ROSCOE [@golovneva2022roscoe] and ReasonEval [@xia2024reasoneval].
Closest to EngTrace's milestones, @wang2026milestoneoracles grade mathematical reasoning against milestone roadmaps with a deterministic symbolic verifier, though their milestones are written by a teacher model rather than computed by the problem's own code.
*Judges.* @thakur2025judgingjudges and @ye2025justiceprejudice document the unreliability and biases of LLM judges; judges accept valid outputs but reject few invalid ones [@jain2025agreeableness], their agreement with human markers in physics depends on the task [@yeadon2026physicsjudge], and they should be validated against human judgments before use [@bavaresco2025llmsinsteadofhumans].
@roytburg2026narcissists re-examine self-preference findings. @jung2025trustescalate escalate only uncertain items to stronger judges, and @huang2026verifierrobustness ask a model only about answers a rule-based verifier rejects, the pattern EngTrace's residual judge follows.

<!-- part: md-notes -->
---

## Notes for the authors

1. **The Introduction cites several of the same works and needs the matching changes.** The "Glue" and "BBH"
   citations point to SuperGLUE and BIG-Bench Extra Hard. Zhou et al. and Mudur et al. are cited as
   relying "solely on outcome matching", which neither does. The opening sentence cites a storytelling
   study. "No existing benchmark verifies this process" is contradicted by PhysReason, FinChain,
   PRISM-Physics and others. Replacement sentences are in `CHANGES.md`, section 2.
2. **Anonymity.** FinChain shares four authors with EngTrace and thanks its first two, so it must be cited
   in the third person, as above. Three 2026 papers cite EngTrace by name and arXiv id: ERI, Sci-ρ and
   LPDS. LPDS also evaluates an earlier public release of EngTrace's templates. The text above cites all
   three without saying so. Whether to acknowledge LPDS's finding about EngTrace in the Limitations is the
   authors' call; the options are in `CHANGES.md`, section 4.
3. **The reference list.** Regenerate it from `related_work_sources.bib`, which takes authors and titles
   from arXiv, Crossref and the ACL Anthology. The May list has wrong first authors for ElecBench, APBench
   and the CIRCUIT entry, and cites the wrong paper for BIG-Bench, GLUE and BBH. It prints a 30-name author
   list for the 10-author Lightman et al. It also has wrong given names in ten entries, and lists as preprints
   papers that have since been published (`CITED_SET_AUDIT.md`).
