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

## 2 Related Work

**Reasoning benchmarks and what they score.**
Mathematical and coding benchmarks such as GSM8K (Cobbe et al., 2021), MATH (Hendrycks et al., 2021b), HumanEval (Chen et al., 2021), MBPP (Austin et al., 2021) and SWE-bench (Jimenez et al., 2024) score a final answer or a test outcome, and broad suites such as MMLU (Hendrycks et al., 2021a) and SuperGPQA (M-A-P Team et al., 2025) score a chosen option.
Yet a correct final answer can follow flawed reasoning (Cobbe et al., 2021), increasingly so on harder problems: in ProcessBench, the share of correct-answer solutions that contain a process error rises from 3.5% on GSM8K to 51.8% on Omni-MATH (Zheng et al., 2025).
Most physics benchmarks likewise grade final answers or expressions (Xu et al., 2025; Qiu et al., 2025; Zhang et al., 2025b; Wang et al., 2024b), while PhysReason (Zhang et al., 2025a), PRISM-Physics (Zhao et al., 2026b) and HiPhO (Yu et al., 2026) score steps against reference solutions, with an LLM judge, rule-based formula matching, and official marking schemes applied by an LLM judge, respectively.

**Generated instances and contamination.**
Static test sets saturate (Ott et al., 2022; Akhtar et al., 2026) and leak into training data (Deng et al., 2024; Zhang et al., 2024), and neither private test sets nor templating has measurably slowed saturation (Akhtar et al., 2026).
Functional, symbolic and perturbed variants of existing problems test whether a model has memorised the original instance or its method (Srivastava et al., 2024; Mirzadeh et al., 2025; Gulati et al., 2025; Huang et al., 2025), although GSM-Symbolic's drops do not always survive re-analysis (Długosz et al., 2026), randomly sampled variants can overstate robustness (Mondorf et al., 2026), and on multiple-choice benchmarks paraphrasing lowers absolute scores while leaving model rankings relatively stable (Lunardi et al., 2025).
Generation has also been paired with gold solutions: HARDMath produces applied-mathematics problems with numerically validated step-by-step solutions (Fan et al., 2025); FinChain pairs each templated financial problem with an executable gold reasoning chain and scores a model's steps against it (Xie et al., 2026); and SymPyBench (Imani et al., 2026), AtmosSci-Bench (Li et al., 2025a) and Sci-ρ (Azmi et al., 2026) instantiate science problems from code templates.

**Engineering benchmarks.**
General suites test engineering mostly through multiple choice; SuperGPQA's 7,892 engineering questions often require calculation but are scored on the option chosen (M-A-P Team et al., 2025).
Dedicated benchmarks typically cover one branch: transportation (Syed et al., 2024), control (Kevian et al., 2024), power dispatch (Zhou et al., 2024), electrical and electronics engineering (Li et al., 2025b), circuits (Skelic et al., 2025), astrodynamics (Wu et al., 2025), thermodynamics (Düzkar, 2026), thermal-protection design (Zheng et al., 2026) or operations research (Mostajabdaveh et al., 2025; Chen et al., 2026).
Multi-branch benchmarks grade answers with LLM judges: EngiBench (Zhou et al., 2026) with two LLM evaluators, plus rubrics for its open-ended level, and ERI (Naser et al., 2026) against reference answers that an LLM wrote along with the questions.
[cut first] Other work evaluates design (Guo et al., 2025), multiphysics simulation (Mudur et al., 2025) and questions from industrial practice (Heesch et al., 2025), the last questioning evaluations built from examination materials whose correctness is easily verifiable.
A few benchmarks look past the answer: experts examined the reasoning in TransportBench and ControlBench, LLM judges score it in ElecBench, TPS-CalcBench and EngVQA (Wasiq et al., 2026), and ThermoQA scores property values and weighted solution steps against values computed with a thermodynamic property library.
Fewer generate their instances: CIRCUIT instantiates circuit templates with several numerical setups, OPT-Engine (Chen et al., 2026) generates optimisation problems of controlled complexity, and a power-systems agent benchmark draws held-out cases from privately seeded generators (Trashchenkov, 2026).
To our knowledge, none combines generated instances across engineering branches with gold traces, deterministic checks of intermediate values, and an evaluator validated against experts' step-level labels.

**Process supervision and the validity of evaluators.**
Outcome supervision can match process supervision on final-answer error, but correct reasoning steps required process-based supervision or reward models that emulate it (Uesato et al., 2022).
Process reward models (PRMs) learn from human step labels (Lightman et al., 2024) or from labels estimated by rollouts (Wang et al., 2024a), which are noisy (Zhang et al., 2025c), and are evaluated on human-annotated first-error benchmarks in mathematics (Zheng et al., 2025; Song et al., 2025); PRMs trained on mathematics perform poorly in other domains (Zeng et al., 2025).
Beyond mathematics, PRIME tests whether verifiers catch answers that the reasoning does not support, in college-level STEM including engineering (Wang et al., 2026a), and ChemCoTBench-V2 checks chemistry steps with deterministic rules (Guo et al., 2026).
LLM judges are the common alternative, but they favour their own outputs and those of related models (Zheng et al., 2023; Panickssery et al., 2024; Spiliopoulou et al., 2025; Li et al., 2026), even on objective rubric items (Pombal et al., 2026); their errors are correlated across models (Goel et al., 2025; Kim et al., 2025), which limits what a panel of judges (Verga et al., 2024) gains by voting; and on difficult objective comparisons many do little better than chance (Tan et al., 2025), agreeing with experts mainly on questions they can answer themselves unless given a correct reference (Krumdick et al., 2026).
Rule-based checkers avoid these biases but can reject correct answers given in unexpected formats (Huang et al., 2026; Ansari et al., 2026).

**Positioning.**
EngTrace brings these lines together.
As in GSM-Symbolic and FinChain, symbolic templates generate the instances, here across five engineering branches and from a private seed, and as in FinChain the gold trace comes from the code that computes the answer.
Intermediate values are checked deterministically against that trace and arithmetic by recomputation; an LLM judge from a model family outside the evaluated roster decides only what these checks leave open, given each milestone's expected value; and the full evaluator, with off-the-shelf PRMs as baselines, is validated against step labels from 15 domain experts.
Expert-checked paraphrases test whether models exploit template wording.

## Appendix: Extended related work

**Table A.** How EngTrace relates to the closest benchmarks. "–" means not provided, or not reported in the sections reviewed for this audit.

| Benchmark | Domain | Instances | Gold steps | What is scored, and how | Scorer checked against people |
|---|---|---|---|---|---|
| FinChain (Xie et al., 2026) | finance | generated from symbolic templates | executable gold chain | steps, against the gold chain | Spearman 0.655 with experts' ratings of reasoning quality |
| HARDMath (Fan et al., 2025) | applied mathematics | generated by code | step-by-step, validated numerically | final answer, graded by an LLM with rubrics | – |
| SymPyBench (Imani et al., 2026) | university physics | parameterised Python | step-by-step reasoning | answers, and consistency across variants | – |
| AtmosSci-Bench (Li et al., 2025a) | atmospheric science | 67 templates with Python solvers, as multiple choice; 391 static open-ended questions | – | the option chosen; open-ended answers by numeric, then symbolic, then LLM checks | its LLM grader agreed with human graders on 92.79% and 93.02% of the open-ended answers it decided (two models' outputs) |
| Sci-ρ (Azmi et al., 2026) | school STEM | 606 Python templates | reasoning steps | steps, as an F1 from an LLM judge | – |
| PhysReason (Zhang et al., 2025a) | physics | static | annotated solution steps | steps, an LLM judge locating the first error | about 98% first-error accuracy on 1,000 annotated solutions |
| PRISM-Physics (Zhao et al., 2026b) | physics | static (PhD qualifying-exam problems) | formula graphs | steps, by rule-based formula matching | Kendall τb 0.346 with two experts' scores of 70 solutions (an outcome-only LLM judge: 0.294) |
| CIRCUIT (Skelic et al., 2025) | circuits | 102 templates, several numerical setups each | – | final answers, per template | human review found about 5% false positives |
| ThermoQA (Düzkar, 2026) | thermodynamics | static (293 problems) | states and steps computed with CoolProp | property values and weighted steps, within 2% | – |
| TPS-CalcBench (Zheng et al., 2026) | thermal-protection design | static (420 items) | – | answer, and a reasoning rubric scored by an LLM judge | weighted κ 0.68 to 0.82 per rubric dimension |
| EngVQA (Wasiq et al., 2026) | engineering, multimodal | static (696 problems) | – | eight solution stages, by an LLM judge | reported against adjusted ("synthetic") human scores |
| EngiBench (Zhou et al., 2026) | several engineering branches | static | – | answers by two LLM evaluators; rubrics for the open-ended level | humans resolve evaluator disagreements and review the open-ended scores |
| ERI (Naser et al., 2026) | nine engineering disciplines | generated by an LLM (57,750 records) | LLM-written references | a five-point score from three LLM judges | – |
| TransportBench (Syed et al., 2024) | transportation | static | – | answers and reasoning, graded by experts | the experts are the graders |
| ChemCoTBench-V2 (Guo et al., 2026) | chemistry | static | reference traces | intermediate steps, by deterministic rules | 87.4% agreement with adjudicated expert labels (κ 0.74; experts with each other, 90.1%) |
| **EngTrace** | five engineering branches | 150 templates, instances from a private seed | from the template code | final answer; intermediate values and arithmetic, deterministically; a residual LLM judge from outside the evaluated roster | step labels from 15 domain experts on 300 traces |

**Further related work.**
*Engineering.* Recent benchmarks extend coverage to agents working in professional engineering software (Gao et al., 2026), strength of materials (Wan et al., 2026), circuit analysis from schematics with synthetically generated problems (Akbari et al., 2026), telecommunications mathematics (Colle et al., 2025), executable structural analysis (Qin et al., 2026), finite-element programming (Mohammadzadeh et al., 2026), computational fluid dynamics (Somasekharan et al., 2026), civil-engineering licensure questions (Silwal et al., 2026) and process systems engineering judged by several LLMs (Sukpancharoen and Srinophakun, 2026).
*Science.* OlympiadBench pairs olympiad problems with expert step-by-step solutions but scores final answers (He et al., 2024); SciDA randomises numerical parameters at each evaluation (Zhou et al., 2025); NEWTON generates physical-commonsense questions from templates (Wang et al., 2023).
*Contamination.* VeRA turns benchmark items into executable task families (Cheng et al., 2026), GSM-Plus perturbs GSM8K problems (Li et al., 2024), Chen et al. (2025) survey the move from static to dynamic benchmarks, and Xiao and Cheng (2026) measure contamination as the gap between original and paraphrased items.
*Step verification.* Error-detection benchmarks extend beyond mathematics (Tyen et al., 2024; Zeng et al., 2024; He et al., 2025; Sun et al., 2026), physics now has a process reward model (Dong et al., 2026), and GenPRM verifies each step with code (Zhao et al., 2026a); reference-free step metrics include ROSCOE (Golovneva et al., 2023) and ReasonEval (Xia et al., 2025).
Closest to EngTrace's milestones, Wang et al. (2026b) grade mathematical reasoning against milestone roadmaps with a deterministic symbolic verifier, though their milestones are written by a teacher model rather than computed by the problem's own code.
*Judges.* Thakur et al. (2025) and Ye et al. (2025) document the unreliability and biases of LLM judges; judges accept valid outputs but reject few invalid ones (Jain et al., 2025), their agreement with human markers in physics depends on the task (Yeadon et al., 2026), and they should be validated against human judgments before use (Bavaresco et al., 2025).
Roytburg et al. (2026) re-examine self-preference findings. Jung et al. (2025) escalate only uncertain items to stronger judges, and Huang et al. (2026) ask a model only about answers a rule-based verifier rejects, the pattern EngTrace's residual judge follows.

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
