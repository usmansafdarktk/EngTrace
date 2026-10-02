# What changes from the May 2026 text, and why

Each change below points to the audit entry that justifies it. `CITED_SET_AUDIT.md` covers the May
citations, and `NEW_PAPERS_AUDIT.md` covers the October search. Replacement sentences are in LaTeX with
citation keys. The keys are in `related_work_sources.bib` or `intro_methods_sources.bib`. Every claim they
make about a cited paper is in `facts.py` and verified by `verify_facts.py`.

## 1. Related Work, paragraph by paragraph

**Mathematics and coding** becomes **Reasoning benchmarks and what they score.**
- **Framing:** "rigorously test abstract logical deduction and reasoning stability" fitted GSM8K and MATH
  loosely and HARDMath not at all. The paragraph now says what these benchmarks score: a final answer, a
  test outcome or a chosen option. It then cites the evidence that correct answers can rest on flawed
  steps: Cobbe et al., and ProcessBench's 3.5% to 51.8%.
- **Physics:** the physics benchmarks move here. The May "symbolic" and "qualitative" labels were wrong
  for ABench-Physics and LLM-SRBench. "Confined to theoretical physics" fitted none of the six. PhysReason
  is now credited with step scoring. PRISM-Physics and HiPhO are added as step-level physics evaluations.
  SciBench, which the July rebuttal used, is added.
- **Heesch et al.** leaves the transfer sentence, because the paper never tests transfer from math or
  code. It is cited in the engineering paragraph for what it is.
- **LLM-SRBench** and **NEWTON** leave the main text. LLM-SRBench is dropped and NEWTON moves to the
  appendix.

**New paragraph: Generated instances and contamination.** The May text treated symbolic templates as
EngTrace's mechanism and cited only GSM-Symbolic and FinChain for them.
- **What the paragraph covers:** the template and perturbation literature; the evidence on saturation
  (Ott et al.; Akhtar et al. 2026) and contamination (Deng et al.; GSM1k); the critiques (Długosz et al.
  2026; Mondorf et al. 2026); paraphrase sensitivity (Lunardi et al. 2025).
- **Precedents that already pair generation with gold solutions:** HARDMath, FinChain, SymPyBench,
  AtmosSci-Bench and Sci-ρ.
- **What the text no longer claims:** that templates counter saturation. Akhtar et al. find that
  safeguards such as private test sets have limited effect, and EngTrace's own scores already run high.

**Engineering benchmarks** is rewritten around what each benchmark scores.
- **Removed:** the "factual recall via multiple choice" sentence, which mis-cited BIG-Bench and
  mis-described MMLU and SuperGPQA.
- **Single-branch list:** it names each benchmark with its correct first author: Zhou et al. for
  ElecBench, Wu et al. for APBench, Skelic et al. for CIRCUIT. It adds ControlBench, ThermoQA,
  TPS-CalcBench and ORQA. The industrial branch had no reference before.
- **Multi-branch benchmarks:** EngiBench is described accurately. Two LLM evaluators score most problems,
  and rubrics apply only to the open-ended level. ERI is added.
- **Reasoning and generation:** the paragraph credits the benchmarks that look past the answer:
  TransportBench, ControlBench, ElecBench, TPS-CalcBench, EngVQA and ThermoQA. It also credits the two
  that generate instances, CIRCUIT and the power-systems agent benchmark.
- **The closing sentence:** "existing benchmarks remain limited to outcome matching" becomes a hedged
  claim about the combination. The old claim was contradicted by seven of the papers the May text itself
  cited.

**New paragraph: Process supervision and the validity of evaluators.** The plan of record
(`docs/NEXT_CYCLE_REVIEW.md` §4.2) asks for this paragraph.
- **Process supervision:** process against outcome supervision (Uesato et al.; Lightman et al.); PRMs and
  their benchmarks, which are mathematical; the failure of math-trained PRMs elsewhere (VersaPRM); and the
  nearest work outside mathematics (PRIME, ChemCoTBench-V2).
- **Judge validity:** self- and family preference, correlated errors, judges near chance on hard
  objective items, and judges needing references.
- **Rule-based checkers:** they have their own failure mode.
- **PoLL:** it is cited for panels of judges, not in support of the single residual judge. PoLL argues
  against single judges.

**New closing paragraph: Positioning.** It states the combination EngTrace adds, crediting FinChain and
GSM-Symbolic for the mechanism. It makes no "first" claim.

**New appendix.**
- **Table A:** compares EngTrace with the fifteen closest benchmarks. Every cell is checked against the
  paper.
- **Further related work:** cites 32 more papers from the search, grouped by topic.

## 2. The Introduction: sentences that cite the same works

These are outside section 2, but they rest on the same citations and two of them repeat claims the audit
found wrong.

**L45-48**, currently citing a general survey and a storytelling study (`xie2023storytelling`, dropped):
> As LLMs expand into high-stakes engineering workflows, rigorous evaluation of their reasoning has become
> essential, not least because a correct answer can follow an invalid reasoning path \citep{zhao2023survey}.

**L77-80**, "yet no existing benchmark verifies this process", which PhysReason, FinChain, PRISM-Physics,
ThermoQA, TPS-CalcBench and others contradict:
> Engineering reasoning requires physically grounded, multi-step problem solving where intermediate steps are
> as consequential as the final result, yet most engineering benchmarks score only final answers, and those
> that assess intermediate steps do so with LLM judges or human graders, or within a single domain
> (Section~\ref{sec:related}).

**L84-87**, on symbolic templates (keep the mechanism, credit the precedents precisely):
> Following symbolic-template benchmarks in mathematics and finance \citep{mirzadeh2024gsmsymbolic,xie2025finchain},
> templates serve as the generative mechanism, each instance paired with a gold reasoning trace computed by
> the same code;

**L97-105**, where "Glue" cited SuperGLUE and "BBH" cited BIG-Bench Extra Hard, and templates were said to
counter saturation:
> First, static benchmarks saturate, as GLUE \citep{wang2018glue} and BIG-Bench Hard \citep{suzgun2022bbh}
> did within a few years \citep{wang2019superglue,kazemi2025bbeh,ott2022saturation,akhtar2026benchmarksaturation},
> and their items leak into training data \citep{deng2024contamination}; EngTrace draws its instances from
> symbolic templates with a private seed, so no evaluated item was public when the models ran, and it tests
> sensitivity to wording with expert-checked paraphrases.

**L105-111**, where MMLU was "broad factual recall" and MATH "abstract logic":
> Second, most benchmarks evaluate skills in disciplinary silos: general benchmarks such as MMLU
> \citep{hendrycks2021mmlu} test knowledge and problem solving through multiple-choice questions, while
> specialized ones such as MATH \citep{hendrycks2021math} or HumanEval \citep{chen2021humaneval} test
> mathematical problem solving or code generation.

**L116-122**, where EngiBench and FEABench were said to "rely solely on outcome matching" (EngiBench uses
rubrics for its open-ended level; FEABench's main metrics are intermediate):
> EngTrace addresses both gaps through integrative reasoning templates whose every instance carries a gold
> reasoning trace, so that intermediate steps can be checked as well as the final answer, a crucial
> requirement in engineering, where a flawed process can lead to catastrophic failure \citep{leveson2012}.

(`leveson2012` is the May entry for Leveson's book; it was not part of this audit.)

## 3. Other sections that cite audited works

- **Section 3.3 (the AI Tribunal, now the Layer 1 screen).** Cite PoLL \citep{verga2024poll} for using
  a panel of judges from different families. Its evidence covers QA and chatbot preferences, not reasoning
  (p.6). Given correlated judge errors \citep{goel2025greatmodels,kim2025correlatederrors}, say that the
  screen's verdicts are checked by experts (Layer 2) rather than trusted as independent votes.
- **Section 4 (evaluation framework).** Replace "Inspired by recent advances in process supervision
  \citep{lightman2023verify} and cascaded model verification \citep{chen2023frugalgpt}" with: "Following
  work on process supervision \citep{uesato2022processoutcome,lightman2023verify}, the evaluator checks
  intermediate steps as well as the final answer; like cascaded evaluation
  \citep{chen2023frugalgpt,jung2025trustescalate}, it sends to an LLM judge only what cheaper,
  deterministic checks cannot settle." FrugalGPT is TMLR 2024, not ICML.
- **Section 4.3 (textual quality).** Drop `xie2023deltascore`. It is a reference-free story metric whose
  own results rank reference-based metrics lowest. If BERTScore and ROUGE stay in an appendix, cite
  \citet{zhang2020bertscore} and \citet{lin2004rouge}. The reason for not relying on them is in their own
  papers: BERTScore's authors note it has difficulty detecting factual errors (p.16). HumanEval's authors
  found BLEU does not track functional correctness (p.6).
- **Appendix K.** Keep \citet{plank2022hlv} only where the appendix discusses leftover disagreement on
  subjective ratings. The new Layer 2 certification has no such analysis yet.

## 4. Anonymity and authorship (decisions for the authors)

- **FinChain** shares four authors with EngTrace, and its acknowledgements thank EngTrace's first two. It
  is cited in the third person throughout, which ARR allows. Do not describe it as "our earlier work".
- **EngTrace's arXiv preprint** (2511.01650; v1 titled "EngChain"; v3 still describes 90 templates) is
  public. ARR permits preprints, but the submission must not cite it as the authors' own.
- **Three 2026 papers cite EngTrace by name:** ERI (Naser et al.), Sci-ρ (Azmi et al.) and LPDS
  (Mondorf et al.). LPDS also evaluates 100 templates derived from EngTrace's public repository and reports
  that randomly sampled instances "can overestimate robustness". The revised text cites LPDS for that
  general finding without naming EngTrace. For the Limitations, there are two options:
  - **A neutral sentence:** "Randomly sampled instances may overstate robustness
    \citep{mondorf2026lpds}; EngTrace samples 15 instances per template across its reasoning paths rather
    than searching for the hardest."
  - **Leave it out:** that is defensible because LPDS evaluated the May templates, not the certified 150.
  
  The authors decide.
- **A note for the authors:** LPDS's Figure 24 (p.34) reproduces an EngTrace spring-mass template with a
  wrong sine term. The repository's template has computed `sqrt(stiffness / mass)` since its first commit
  (108b0d6), as it does now. The error is LPDS's transcription, and nothing needs fixing.

## 5. Critiques the paper should answer, from the new literature

- **Exam-style items.** Heesch et al. criticise engineering evaluations "adapted from examination materials
  where correctness is easily verifiable" (p.1). EngTrace's items are exactly that. The Limitations should
  say why: verifiable gold traces need it. It should also say what that leaves out: ill-posed problems and
  real artefacts. The May Limitations already gestures at this.
- **Saturation.** EngTrace's answer scores on the new roster are high (state brief §4). The text should
  not present templates as a remedy for saturation (Akhtar et al. 2026). The private seed answers
  contamination, and the step-level checks answer the question that a saturated answer score no longer
  can.
- **Random instance sampling** (Mondorf et al. 2026): see section 4.

## 6. Things that could not be confirmed

- **Venues from one source:** EMNLP 2026 for Huang et al. and Długosz et al. (each from an arXiv comment)
  and Computers & Industrial Engineering for ERI. Each is flagged in `related_work_sources.bib`. Confirm
  them before submission.
- **LSR-Ben's venue:** its PDF header says "Published as a conference paper at ICLR 2027", which cannot be
  true in October 2026. It is cited as a preprint.
- **Abstract-only papers:** PE Civil Bench and PSE-Bench are on ScienceDirect, which blocks automated
  download. Only their abstracts were read (OpenAlex), and the appendix says nothing beyond them.
- **Gaps in the search:** the shared web-search quota ran out partway through. Later searches used the
  arXiv API, OpenAlex and direct page fetches, so journals that do not post to arXiv (IEEE, ASCE,
  Elsevier) are under-sampled.
- **Gaps in the panel:** the review panel covered every May citation (three reviews each) and the
  engineering and process-supervision must-cites (two or three each). It was stopped before reviewing the
  rest, at the owner's instruction, to save cost (`notes/review_protocol.md`). Every claim the new text
  makes about those papers was instead checked against their text (`notes/fact_check.md`).
