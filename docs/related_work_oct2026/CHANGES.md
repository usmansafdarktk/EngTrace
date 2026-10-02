# What changes from the May 2026 text, and why

Each change below points to the audit entry that justifies it. `CITED_SET_AUDIT.md` covers the May
citations, and `NEW_PAPERS_AUDIT.md` covers the October search. Replacement sentences are in LaTeX with
citation keys. The keys are in `related_work_sources.bib` or `intro_methods_sources.bib`. Every claim they
make about a cited paper is in `facts.py` and verified by `verify_facts.py`.

The follow-up review of 2026-10-02 (`notes/followup_review.md`) read the 18 main-text papers the panel
had not reviewed and made a final pass over the section. Its changes are folded into the entries below.

## 1. Related Work, paragraph by paragraph

**Mathematics and coding** becomes **Reasoning benchmarks and what they score.**
- **Framing:** "rigorously test abstract logical deduction and reasoning stability" fitted GSM8K and MATH
  loosely and HARDMath not at all. The paragraph now says what these benchmarks score: a final answer, a
  test outcome or a chosen option. It then cites the evidence that correct answers can rest on flawed
  steps: Cobbe et al., and ProcessBench's 3.5% to 51.8%.
- **Physics:** the physics benchmarks move here. The May "symbolic" and "qualitative" labels were wrong
  for ABench-Physics and LLM-SRBench. "Confined to theoretical physics" fitted none of the six. PhysReason
  is now credited with step scoring. PRISM-Physics and HiPhO are added as step-level physics evaluations.
  PRISM-Physics' step check is named as rule-based: it uses no LLM, which makes it the nearest physics
  precedent for EngTrace's deterministic checks. SciBench, which the July rebuttal used, is added.
- **Graders that fail:** the sentence on physicists re-grading 250 rejected answers (Ansari et al. 2026)
  was cut, as the first [cut first] sentence. The paper is now cited for its finding that rule-based
  graders reject correct answers in equivalent forms, in the paragraph on evaluators.
- **Heesch et al.** leaves the transfer sentence, because the paper never tests transfer from math or
  code. It is cited in the engineering paragraph for what it is.
- **LLM-SRBench** and **NEWTON** leave the main text. LLM-SRBench is dropped and NEWTON moves to the
  appendix.

**New paragraph: Generated instances and contamination.** The May text treated symbolic templates as
EngTrace's mechanism and cited only GSM-Symbolic and FinChain for them.
- **What the paragraph covers:** the template and perturbation literature; the evidence on saturation
  (Ott et al.; Akhtar et al. 2026) and contamination (Deng et al.; GSM1k); the critiques (Długosz et al.
  2026, who re-analyse GSM-Symbolic only; Mondorf et al. 2026); paraphrase sensitivity on multiple-choice
  benchmarks (Lunardi et al. 2025). Variants test memorisation of the original instance or, for
  MATH-Perturb's hard perturbations, of its solution method.
- **Precedents that already pair generation with gold solutions:** HARDMath, FinChain, SymPyBench,
  AtmosSci-Bench and Sci-ρ.
- **What the text no longer claims:** that templates counter saturation. Akhtar et al. find that neither
  private test sets (their H1) nor templating (H6: 14 templated against 46 other benchmarks, p = 0.10)
  measurably slows saturation, and EngTrace's own scores already run high.

**Engineering benchmarks** is rewritten around what each benchmark scores.
- **Removed:** the "factual recall via multiple choice" sentence, which mis-cited BIG-Bench and
  mis-described MMLU and SuperGPQA.
- **Single-branch list:** it names each benchmark with its correct first author: Zhou et al. for
  ElecBench, Wu et al. for APBench, Skelic et al. for CIRCUIT. It adds ControlBench, ThermoQA,
  TPS-CalcBench, ORQA and OPT-Engine. The industrial branch had no reference before. ORQA is multiple
  choice on formulating optimisation models; OPT-Engine generates inventory, production and other
  operations-research problems with verified optimal solutions, the closer match to that branch.
- **Multi-branch benchmarks:** EngiBench is described accurately. Two LLM evaluators score most problems,
  and rubrics apply only to the open-ended level. ERI is added.
- **Reasoning and generation:** the paragraph credits the benchmarks that look past the answer:
  TransportBench, ControlBench, ElecBench, TPS-CalcBench, EngVQA and ThermoQA. Its experts "examined"
  the reasoning in TransportBench and ControlBench: ControlBench's judged each answer right or wrong and
  sorted the failures into seven error types, so "graded" said too much. It also credits the three that
  generate instances: CIRCUIT, OPT-Engine and the power-systems agent benchmark.
- **The closing sentence:** "existing benchmarks remain limited to outcome matching" becomes a hedged
  claim about the combination. The old claim was contradicted by seven of the papers the May text itself
  cited.

**New paragraph: Process supervision and the validity of evaluators.** The plan of record
(`docs/NEXT_CYCLE_REVIEW.md` §4.2) asks for this paragraph.
- **Process supervision:** process against outcome supervision (Uesato et al.; Lightman et al.); PRMs and
  their benchmarks, which are mathematical; the failure of math-trained PRMs elsewhere (VersaPRM); and the
  nearest work outside mathematics (PRIME, whose problems are college-level STEM, 18% of them engineering by
  its Figure 3, and ChemCoTBench-V2).
- **Judge validity:** self- and family preference, correlated errors, many judges near chance on hard
  objective items (reasoning models do better on JudgeBench), and judges needing correct references.
- **Rule-based checkers:** they have their own failure mode (Huang et al. 2026; Ansari et al. 2026).
- **PoLL:** it is cited for panels of judges, not in support of the single residual judge. PoLL argues
  against single judges.

**New closing paragraph: Positioning.** It states the combination EngTrace adds, crediting FinChain and
GSM-Symbolic for the mechanism. It makes no "first" claim. It says the arithmetic check recomputes each
claim (it checks the model's own numbers, not the gold trace), and that the judge decides only what the
deterministic checks leave open, which covers the step router as well as the residual milestone judge.

**New appendix.**
- **Table A:** compares EngTrace with the fifteen closest benchmarks. Every cell is checked against the
  paper. The follow-up review corrected two rows: AtmosSci-Bench's generated instances are multiple
  choice, and its LLM grader was checked against human graders (92.79% and 93.02% agreement); PRISM-Physics'
  validation is two experts' scores of 70 solutions, against an outcome-only LLM judge.
- **Further related work:** cites 32 more papers, 31 of them from the search (NEWTON is the other),
  grouped by topic.

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
> First, static benchmarks saturate \citep{ott2022saturation,akhtar2026benchmarksaturation}, as GLUE
> \citep{wang2018glue} and BIG-Bench Hard \citep{suzgun2022bbh} did within a few years
> \citep{wang2019superglue,kazemi2025bbeh}, and their items leak into training data \citep{deng2024contamination};
> EngTrace draws its instances from symbolic templates with a private seed, so no evaluated item was public
> when the models ran, and it tests sensitivity to wording with expert-checked paraphrases.

Akhtar et al. are cited for saturation in general, not for the two examples. By their index, which
measures lost discrimination among the strongest models, BIG-Bench Hard is among the benchmarks that
"remain unsaturated despite prolonged exposure" (p.6). BBEH, cited for BBH, reports frontier scores above
90%. The earlier version of this sentence put the four citations in one bracket after both examples.

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

- **Venues from one source:** EMNLP 2026 for Huang et al. and Długosz et al. still rests on each paper's
  arXiv comment, and neither has a proceedings record yet. OpenReview is consistent with both: it shows
  Długosz et al. as an ARR May 2026 submission and Huang et al. as withdrawn from ICLR 2026
  (`notes/venue_check.json`). Both stay flagged in `related_work_sources.bib`; confirm them before
  submission.
- **ERI's venue is confirmed.** Its arXiv page links DOI 10.1016/j.cie.2026.112333. Crossref gives that
  DOI as Computers & Industrial Engineering, volume 221, and Semantic Scholar agrees.
- **Venues the papers do not state, checked in October** (`check_venues.py`):
  - ACL 2026 for PRIME (ACL Anthology 2026.acl-long.683);
  - ICLR 2026 for PRISM-Physics and ICML 2026 for HiPhO (OpenReview);
  - ICML 2025 for MATH-Perturb (OpenReview).
- **Two papers the catalogue listed as preprints have since appeared.** The bibliography now cites
  SymPyBench at the EACL 2026 Industry Track (2026.eacl-industry.8) and AtmosSci-Bench at the NeurIPS 2025
  Datasets and Benchmarks Track (OpenReview), through `make_bib.py`'s `VENUE_CONFIRMED`.
- **MATH() and ControlBench** have no venue in any index, and stay preprints.
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
