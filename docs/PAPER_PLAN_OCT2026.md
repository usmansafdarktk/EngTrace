# EngTrace, October 2026 submission: the paper as it stands, and how to write it

For the authors. Written 2026-10-03, revised 2026-10-04 after their review. Deadline: 12 October 2026 (ARR, anywhere
on Earth). Every number below names the file that prints it; nothing here is a new measurement.

**What this plan reads**

- The May 2026 submission (`docs/_ARR_May__EngTrace.pdf`), both rebuttals, both meta-reviews (`docs/meta-reviews.md`),
  and the current Overleaf sources (`current_overleaf_project/`).
- The template rework (`docs/re-implementation-sep/`) and the certification (`template_annotation_23092026/`).
- The evaluator study (`evaluator_pilot_17092026/`) and the full evaluation (`full_run_28092026/`; `results/RESULTS.md`
  at commit `451c111`, with the open-book and tool conditions complete).
- The revised Related Work (`docs/related_work_oct2026/`), the five-branch figure (`figures_oct_12/`) and the experts'
  readings (`full_run_28092026/EXPERT_REQUEST.md`).

**Companions**

- `docs/REVISION_LETTER_NOTES.md`: what changed since May and why. For the letter only.
- `docs/OVERLEAF_WORKFLOW.md`: how the Overleaf project is kept in sync.
- `full_run_28092026/RESULTS_PAPER_NOTES.md` and `PARAPHRASE_PAPER_NOTES.md`: number by number, what each result
  supports and does not.

**The rule that shapes everything.** The paper is written for readers who have never seen EngTrace. It describes the
benchmark, the evaluator and the results as they are, with their evidence and their limits. No "previously", no "in
the earlier version", no "as promised". And no internal vocabulary: section 8 gives the replacement for every
internal term.

---

## 1. The paper in one paragraph

EngTrace is a benchmark of 150 symbolic templates across five engineering branches (chemical, civil, electrical,
industrial, mechanical), 15 domains and 42 areas. Each template is executable code: it samples physically grounded
parameters from referenced constants, computes the answer, and writes the gold reasoning trace from the same
computation, so every intermediate quantity exists with its unit before any prose does.

Templates reach the benchmark through automated integrity checks, an LLM screen whose judges share no family with
any evaluated model, and certification by three experts of the template's own branch, who caught 60 of 60 planted
defects and rejected 53 templates across rounds until all 150 were approved unanimously. The evaluation set (2,250
instances, 15 per template) is drawn from a private seed, so no evaluated instance was public when the models ran.

A response is scored deterministically wherever a check exists: a three-way final-answer check over six answer
kinds, coverage of the gold derivation's intermediate quantities, and an arithmetic check of every displayed
calculation. A judge from a family outside the evaluated models decides only what those checks leave open. The
evaluator was chosen by comparing seven designs against 15 experts' step-level labels on 300 responses and 120
planted defects, which also fixed the limit of verification: no deterministic check detects a misstated rule behind
a correct answer, and LLM judges catch about a third.

Eleven open and closed models of September 2026 score 0.81 to 0.98 on the final answer; the top five are inseparable
at this size; every model scores lower on Advanced templates; milestone coverage orders the models differently from
the answer; on wrong answers models still reach most of the derivation; expert-checked paraphrases change no model's
score by more than five points for ten of the eleven; and neither the governing equations in the prompt nor a code
tool moves the strongest models' scores.

## 2. Title, abstract and contributions

**Title.** Keep *EngTrace: A Symbolic Benchmark for Verifiable Process Supervision of Engineering Reasoning*. The
mechanism now matches the title. A shorter form, if wanted: *EngTrace: Verifiable Process Evaluation of Engineering
Reasoning from Symbolic Templates*.

**Abstract, draft** (about 250 words; trim to taste).

> Engineering problems are solved by derivations whose intermediate quantities can be checked, yet most benchmarks
> score only the final answer, and those that assess steps rely on LLM judges or human graders. We present
> EngTrace, a benchmark of 150 symbolic templates spanning five engineering branches, 15 domains and 42 areas. Each
> template is executable code that samples physically grounded parameters, computes the answer and emits a gold
> reasoning trace from the same computation, so every intermediate quantity is known with its unit. Templates pass
> automated integrity checks and an LLM screen, and are certified by three experts of their own branch, who caught
> 60 of 60 planted defects. Evaluated instances are drawn from a private seed, so none was public when the models
> ran. EngTrace scores a response deterministically wherever a check exists: a three-way final-answer check over six
> answer kinds, coverage of the gold derivation's intermediate quantities, and an arithmetic check of every displayed
> calculation; a judge from a model family outside the evaluated models decides only what these checks leave open.
> The evaluator is validated against step-level labels from 15 experts on 300 responses (answer agreement 0.98,
> milestone F1 0.96) and against 120 planted defects, which also show the limit of verification: no deterministic
> check detects a misstated rule behind a correct answer, and judges catch about a third. On 2,250 instances, eleven
> models score 0.81 to 0.98 on the final answer, with the top five inseparable at this size; every model scores lower
> on Advanced templates; milestone coverage orders the models differently; and on wrong answers models still reach
> 46 to 85% of the derivation. Expert-checked paraphrases change no model's score by more than five points for ten
> of eleven models, and neither the governing equations in the prompt nor a code tool moves the strongest models.

**Contributions, draft.** Three, each owned by an active verb and each carrying its evidence.

1. **We build and certify EngTrace**, 150 symbolic templates across five engineering branches whose gold traces are
   computed by the template code, with an evaluation set drawn from a private seed and a three-layer certification
   whose evidence is reported: automated integrity checks, a screen by LLM judges from outside the evaluated
   families, and own-branch expert certification with hand checks and planted defects (60 of 60 caught), which
   rejected 53 templates before all 150 were approved unanimously.
2. **We design and validate a deterministic-first process evaluator**: a final-answer check over six answer kinds,
   milestone coverage derived from the template's own computed quantities, an arithmetic check of displayed
   calculations, and a judge from outside the evaluated families only where these checks leave a question open. We
   choose it by comparing seven evaluator designs, including a step-matching design with a panel of LLM judges and
   three process reward models, against 15 experts' step-level labels on 300 responses and 120 planted defects, and
   we report what verification can and cannot see.
3. **We evaluate eleven open and closed models** with confidence intervals and significance tests that treat the
   template as the unit, and report tiers rather than an order at the top, the drop on Advanced problems for every
   model, what milestone coverage adds on wrong answers, failure against derivation depth, a bound on the effect of
   rewording, and what reasoning effort, flagship models, the governing equations and a code tool change.

## 3. The claims the data support, and the claims not to make

Each claim is a sentence the paper may make, in the paper's own words, with its source. The "do not" list is binding;
it collects every overstatement the record has already caught (`RESULTS_PAPER_NOTES.md`, its "What not to state"
sections; `PARAPHRASE_PAPER_NOTES.md`; D-111, D-146, D-171, D-180 to D-184).

### 3.1 Claims

| # | Claim | Evidence (file) |
|---|---|---|
| C1 | 150 templates, 30 per branch, 15 domains, 42 areas, 58 / 58 / 34 templates by level; 2,250 instances (870 / 870 / 510; 450 per branch), all distinct questions, drawn from a private seed with the instances' hashes and a commitment to the seed published | `full_run_28092026/FREEZE.json`, `manifest.jsonl`, D-114 to D-116 |
| C2 | Every template passes automated integrity checks at 500 seeds (every printed calculation closes, generation is deterministic across processes, the output format holds); a screen by three LLM judges from outside the evaluated families passed 147 of 150 with 3 controversial and 0 critical; three own-branch experts certified all 150 after 53 rejections and fixes, catching 60 of 60 planted defects, with 459 of 504 hand checks within 1% | `template_annotation_23092026/layer0/gate_report.md`, `screen/pass2/stats.md`, `layer2/RESULTS.md`, `CERTIFICATION.md`, D-109 |
| C3 | The deterministic components agree with 15 experts' labels on 300 responses: the final-answer check 0.982 (correct or not) and 0.927 three-way; milestone F1 0.923 without the judge and 0.958 with it; the judge credited 0 of 88 fabricated values; the arithmetic check has precision 0.817 and recall 0.427 on slips inside correct-answer responses | `full_run_28092026/SCORER_VALIDATION.md`, `E5_VALIDATION.md`, pilot `RESULTS_X1.md`, `THRESHOLD_APPENDIX.md` |
| C4 | No deterministic check detects a misstated rule behind a correct answer (0 of 60 planted conceptual defects for every checker); LLM judges catch about a third (16 of 52, 20 of 60, 22 of 60 for three judges) with no false alarm on untouched steps; the final answer nearly decides the experts' verdict on a response (AUROC 0.974) | pilot `RESULTS_X1.md` Findings 1, 5, 7, 8; D-181 |
| C5 | Final Answer Accuracy 0.814 (gpt-oss-20b) to 0.976 (DeepSeek V4.1 Flash); the top five lie within 0.009 and no pair among them separates; 32 of 55 pairs differ after correction; 23 pairs lie below the smallest difference the design detects, all of them non-significant | `results/RESULTS.md` Q1 |
| C6 | Every model scores lower on Advanced than on Easy templates (+0.055 to +0.199); the gap holds after correction for the four lowest-scoring models and for none once two Advanced chemical templates whose wording does not pin the answer are removed | Q2 and its variations; D-171; `EXPERT_REQUEST.md` B4 |
| C7 | Milestone Coverage runs 0.816 to 0.923 and orders the models differently from the answer (Kendall's tau 0.709, 95% CI 0.514 to 0.855); 28 of 55 pairs differ; coverage does not rise with the number of steps written (rho −0.11 to −0.26) | Q3 "Coverage compared across models", "Does coverage track verbosity?" |
| C8 | On wrong answers models still reach 46% to 85% of the milestones (readable responses) against a chance floor of 10% to 22%; 13% to 26% of wrong answers are a complete derivation to a wrong value (the five models with more than 100 wrong answers) | Q3 "Milestones on the wrong-answer traces", "Verdict against coverage" |
| C9 | The arithmetic check flags 0.5% to 22.5% of correct-answer responses, reading 2.8 to 9.8 calculations per response, at precision 0.905 on these models (171 of 189 flags read by a domain expert); the judged step check flags 1.1% to 22.0% | Q3 "The digit rule", "The step router"; `FLAG_REVIEW_3.md` |
| C10 | The wrong-answer rate on instances whose gold derivation has six or more milestones is 0.037 to 0.291, against 0.000 to 0.045 on one-milestone instances, for every model | Q3 "Wrong-answer rate against the item's milestone count" |
| C11 | The strongest models solve all 15 instances of 85% to 91% of the single-path templates; the weakest 36% to 48% | Q4 |
| C12 | On 277 expert-kept paraphrase pairs over 115 templates, the paired change in Final Answer Accuracy is −0.031 to +0.020, none significant after correction; for ten of eleven models the 90% interval lies within ±5 points; the tiers hold and the order within the top tier is not resolved (tau 0.673 against a sampling-noise median of 0.782) | Q5; `PARAPHRASE_PAPER_NOTES.md` |
| C13 | Repeated decoding on four models (300 instances, three repeats): SD 0.004 to 0.013, the same verdict every time on 86% to 94% of instances | "Decoding repeats" |
| C14 | The ordering is unchanged at half and double the answer tolerance (tau 0.927), without the six templates whose answers can be read off the question (within 0.002) and without the nine symbolic-answer templates (+0.002 to +0.011) | Sensitivity; `THRESHOLD_APPENDIX.md`; `SHORTCUT_AUDIT.md` |
| C15 | Reasoning at medium effort moves GPT-5.4 mini from 0.858 to 0.951 on the 450-instance subset (+0.093, 95% CI 0.056 to 0.133) and Gemini 3.1 Flash-Lite by less than the design detects (+0.016, −0.010 to 0.042) | "C1 and C4"; D-180 |
| C16 | Two flagships that pass the selection rule do not exceed the top tier on the same 450 instances: DeepSeek V4 Pro 0.968 and GPT-5.4 with reasoning on 0.980 sit inside the top five's intervals (0.964 to 0.976); GPT-5.4 at its provider's default (no reasoning tokens) scores 0.941, below them | "C3"; D-182 |
| C17 | Supplying each template's governing equations with the question (405 instances) lifts gpt-oss-20b by +0.064 (95% CI 0.027 to 0.104; +0.043 on the instances answered in both conditions, since fewer responses run out of room) and changes GPT-5.4 mini (−0.019) and Claude Sonnet 5 (+0.007) by amounts bounded inside ±5 points, while coverage rises by 0.03 to 0.06 for all three; the condition supplies equations, not data tables | "C1 and C4", the `openbook2` rows; D-183 |
| C18 | Offered a Python tool with the prompt otherwise unchanged (450 instances), Claude Sonnet 5 calls it on 66% of instances and GPT-5.4 mini on 19%, and neither model's Final Answer Accuracy (+0.004; −0.004) or coverage changes beyond ±5 points; GPT-5.4 mini's arithmetic flags on correct-answer responses fall from 0.127 to 0.093; the open-weight model could not be served with a tool under the provider rule | "C1 and C4", the `tool` rows; D-184 |
| C19 | Judge independence: the judge shares no family with any evaluated model; on a 220-response sample a second judge from a third family moves no model's coverage beyond its interval (−0.017 to +0.020); on a panel-of-judges design replayed offline, dropping a family's own judge changes its score by at most 0.006, the same as a placebo drop, and every judge is lenient on every family alike | `JUDGE_SWAP.md`; pilot `RESULTS_LOJO.md`; D-174, D-181 |
| C20 | On these models, two experts read 150 of the top five models' final-answer verdicts (half incorrect or partial by design): correct verdicts confirmed in 142 of 150 readings, incorrect verdicts in 83 of 98 (12 expert-correct), partial verdicts called fully correct in 34 of 52, so the half-point is conservative; the judge's "reached" verdicts confirmed in 85 of 100 readings and its "missing" verdicts in 79 of 100; the experts agree with one another on 97% of double-read answers and 88% of milestones | `EXPERT_REQUEST.md` B1, B3 |
| C21 | Three readers labelled 160 wrong-answer responses of four contrastive models with the six-category taxonomy plus "No error" (Fleiss' kappa 0.930): for Claude Sonnet 5, 26 of 40 item majorities are "No error"; for gpt-oss-20b, Gemma 4 and GPT-5.4 mini, Calculation Error leads (17, 28 and 28 of 40), then Formula/Principle; Easy wrong answers are almost all calculation (69 of 78 readings) | `EXPERT_REQUEST.md` B2 |

### 3.2 Do not state

- Not "model X is the best": tiers, not an order, at the top; the top five are within 0.9 points and no pair separates.
- Not "trade-off", "decoupling", "trace fidelity" or "reasoning quality" for Milestone Coverage: it is coverage of
  stated intermediates; no test of a trade-off was planned and the intervals overlap.
- Not "x% of correct answers rest on flawed reasoning" from either diagnostic: flag rates with known precision
  (0.905; about 0.70) and recall below one half.
- Not a ranking of models by arithmetic reliability: the check reads 2.8 to 9.8 calculations per response depending
  on style.
- Not "the drop on Advanced problems holds" without both counts (4 of 11 as scored; 0 of 11 without the two chemical
  templates) and the per-model intervals.
- Not "closed models underperform open ones" or any comparison at equal effort: four models ran without reasoning
  tokens.
- Not "chemical engineering is hardest" or any branch order: one branch pair in the whole evaluation holds
  (gpt-oss-20b, electrical above civil); domain means are descriptive.
- Not "2,250 unique problems": 2,250 instances, all distinct questions, from 150 templates.
- Not "robust to paraphrasing" or "no contamination": a bound of ±5 points for ten models on wording memorisation;
  method memorisation is not tested.
- Not "the ranking is stable under rewording": tau 0.673 sits at the lower edge of the sampling-noise distribution;
  the tiers hold.
- Not "flagships do no better" and not "solved at the top" without the Advanced means beside it (0.87 to 0.95 for
  the anchors, 0.92 to 0.95 for the top five, wide intervals).
- Not "retrieval would not help" or "tools do not help": the open-book condition supplies equations only, no data
  tables; the tool condition covers the two closed models; neither says where the remaining errors lie, which is
  the experts' reading.
- Not gpt-oss-20b's open-book gain as a retrieval effect alone: part of it is responses that no longer run out of
  room (+0.043 on the instances answered in both conditions).
- Not "no self-preference" as evidence that LLM judges are accurate: they are uniformly lenient.
- Not the agreement of the final-answer check with the experts as one population rate: the 150-verdict sample is
  stratified by verdict; report it per verdict category.
- Not templates as a remedy for benchmark saturation (Akhtar et al. 2026 find none); the private seed answers
  contamination, the step-level checks answer the question a saturated answer score cannot.
- Not "no existing benchmark verifies the process"; the claim is the combination, hedged as in `RELATED_WORK_v3.md`.
- Not any number beside a May number; the paper describes the evaluator it uses.

## 4. Two metrics and two diagnostics

The paper names two scores and two diagnostics, and nothing else. Level gaps, depth, consistency and the paraphrase
change are presented as analyses of the two scores ("Final Answer Accuracy by level", "change in Final Answer
Accuracy under rewording"), not as metrics with names. Every score carries a 95% interval that resamples templates;
every null carries the smallest difference the design detects.

| Role | Name | Definition, as the paper gives it |
|---|---|---|
| Final answer | **Final Answer Accuracy (FAC)** | Mean per-instance score: 1 when every asked quantity matches the gold within the tolerance, 0.5 when some but not all do, 0 otherwise or when no final answer is readable. The tolerance is 0.002 relative or one unit of the last displayed digit, after unit conversion; targets are the quantities the template computed; six answer kinds (scalar, vector, array, symbolic, classification, multipart). |
| Reasoning | **Milestone Coverage (MC)** | Share of the gold derivation's intermediate quantities (milestones) that the response reaches: matched deterministically within 0.5% under unit conversion, order-free; the milestones the matcher does not find go to a judge from outside the evaluated families, and only its "reached" verdict is credited. Measures progress through the derivation, not the absence of error. |
| Diagnostic | **arithmetic flags** | Share of correct-answer responses with at least one displayed calculation that, recomputed from its printed operands, is not a correct rounding at the printed precision; reported with the calculations read per response and the measured precision (0.905). |
| Diagnostic | **judged step flags** | Share of correct-answer responses in which the judge flags a step the arithmetic check did not; reported with its precision and recall on the experts' labels. |

Two choices inside this scheme, adopted as recommendations and open to reversal:

- **FAC keeps partial credit.** It is the headline the analysis specified; the experts called the check's partial
  verdicts fully correct in 34 of 52 readings, so half a point is the conservative middle; and strict accuracy
  changes no conclusion (54 of 55 pairs agree) and sits in Appendix I.
- **"Milestone Coverage"** rather than "Derivation Coverage": the word milestone is defined once and then does work in
  every table. The choice is cosmetic.

**Table 1** is therefore: model, FAC with its interval, MC with its interval, and a compact letter display (models
sharing a letter do not differ after correction), so the tiers are visible without another column. A footnote marks
GLM-5.3's 100 responses with no readable answer. The strict FAC, the share of responses with no readable answer, the
per-template SD, the judged fraction and the two diagnostics go to Appendix I.

## 5. Sections

ARR long paper: eight pages of content, unlimited references and appendices; the review version must stand on its
own. Few numbered subsections; most structure is paragraphs with bold run-in headings (`\paragraph{}`), as the
current Overleaf sources already do. The budget sums to eight pages.

### 1 Introduction (1.0 page)

- Open with the gap the revised Related Work supports: engineering reasoning is a derivation whose intermediate
  quantities can be checked, yet most engineering benchmarks score only final answers, and those that assess steps
  do so with LLM judges or human graders, or within one domain.
- EngTrace in one paragraph (section 1 above), then Figure 1: a template, one instance, its gold trace, and the
  milestones and calculations the evaluator checks in it.
- The scope sentence: a depth-oriented evaluation of a principled subset of engineering; items are examination-style
  by construction, which is what verifiability needs and which leaves out ill-posed problems and real artefacts.
- The three contributions.
- `docs/related_work_oct2026/CHANGES.md` section 2 has the replacement sentences for every Introduction sentence
  that cites a prior benchmark; use them as written.
- Do not: "first"; "no existing benchmark verifies this process"; "safety-critical" beyond one motivating clause.

### 2 Related Work (0.75 page)

- The main text of `docs/related_work_oct2026/RELATED_WORK_v3.md` as five paragraphs (`related_work_v3.tex`;
  bibliography `related_work_sources.bib`, 100 entries). The file names the cut order from 704 words to about 0.75
  page. Table A and the further-related-work paragraphs go to Appendix M.
- Anonymity: FinChain in the third person; the EngTrace preprint not cited as the authors' own; ERI, Sci-rho and
  LPDS cited without saying they cite EngTrace. Whether the Limitations acknowledge LPDS's finding on randomly
  sampled instances is the authors' call (`CHANGES.md` section 4).
- Venues to confirm before submission: EMNLP 2026 for Huang et al. and Dlugosz et al.

### 3 The EngTrace Benchmark (1.6 pages)

**3.1 Taxonomy and Scope (0.3 page).** Five branches, 15 domains, 42 areas; 30 templates per branch, 10 per domain
except chemical (10 / 12 / 8). Figure 2 is the five-branch overview (`figures_oct_12/engtrace-overview-5branch.pdf`,
copied into `figs/`).

- *Domains and areas.* The lists follow curricular standards (ABET; the societies: ASME, IEEE and AIChE for the
  original three branches, ASCE and IISE for civil and industrial) and the textbooks of Appendix B. For the chemical,
  electrical and mechanical branches the lists were additionally cross-checked with an LLM panel prompted as a
  curriculum designer (Appendix A). For civil and industrial, the domains and areas were chosen by the colleague who
  wrote the templates, from the ABET, ASCE and IISE curricula and the branch textbooks. The paper says exactly that,
  and lets the own-branch expert certification stand as the validation of all 150 templates. No panel run on the new
  branches (decided 2026-10-04).
- *Difficulty.* Easy, Intermediate, Advanced by conceptual complexity, mathematical sophistication and procedural
  depth; 58 / 58 / 34 templates. One sentence on who assigned the levels, for the original branches and for the new
  ones (the 2026-09-05 inventory holds them).

**3.2 Templates and Instances (0.6 page).**

- *Template structure.* A seeded function that samples parameters, often inside a rejection loop enforcing physical
  validity (laminar regime, elastic range, upright floating); computes every quantity; binds each displayed value
  through its display before using it; asserts physical invariants; and emits the question and the gold trace with
  numbered steps and an answer line. Six answer kinds.
- *Grounding.* Every constant table carries a provenance class and resolves to a referenced source (NIST, CODATA,
  JANAF, NASA, military and federal handbooks, the textbooks; Appendix B).
- *Structural variation.* 72 of 150 templates change the governing equation, the step count or the computed
  quantities with a sampled parameter; the evaluation set holds a second reasoning path for every template that can
  produce one within 500 draws; 58 templates follow one path (`template_audit_report.md`, `DIVERSITY.md`).
- *The evaluation set.* 15 instances per template from a private 128-bit seed, chosen round-robin across each
  template's reasoning paths and answer forms, excluding repeated questions and display ties. The instances' hashes
  and a commitment to the seed are published now; the seed and the text at publication.

**3.3 Certification (0.7 page).**

- *Authorship.* All 150 templates were written by people with engineering training and grounded in textbooks: the
  chemical, electrical and mechanical templates by the team's domain experts, the civil and industrial templates by
  a colleague of the authors. The exact description of that colleague (role, expertise) is still to be supplied, as
  is a decision on whether any writing assistance needs disclosing under the ARR checklist. The textbooks and data
  sources the colleague used are listed under Appendix B below.
- *Automated integrity checks.* At 500 seeds per template: every printed calculation closes from its printed operands
  at the printed precision; generation is deterministic across processes; the output format holds; every seed
  generates. 150 of 150 pass; 54 templates were edited to pass; a census of display ties.
- *LLM screen.* The rubric of Appendix E; three judges from no evaluated family (Grok 4.6, MiniMax M3,
  MiMo-V2.5-Pro); two passes; 147 pass, 3 controversial, 0 critical on the shipped templates; agreement AC1 0.93 on
  the flag; every flag verified before any edit, 34 of 45 claims confirmed.
- *Expert certification.* Three experts of the template's own branch, 34 items each: the branch's 30 templates and 4
  planted defects, one of each kind (wrong constant, wrong unit conversion, flipped sign, a printed step that does
  not follow). A hand check before the solution is shown; opaque codes; timestamps.
- *Outcome.* 60 of 60 planted defects rejected; 459 of 504 hand checks within 1%, 44 of the 45 mismatches ending in a
  rejection; AC1 0.914 on approve/reject (Fleiss' kappa 0.589 beside it, deflated by prevalence); AC2 0.947 / 0.980 /
  0.954 on the three quality scores; 53 rejections over the rounds, each verified against its instance and fixed or,
  in two cases, answered; every civil and industrial expert approved all 30 of their templates in round 1, and the 22
  round-1 rejections fell on the original three branches; rounds 2 to 4 re-certified the changed templates; 150 of
  150 certified unanimously on code that still regenerates the reviewed instances byte for byte. The screen's
  false-positive rate against the experts (10.2%) is why the screen is a reader and the experts the gate.
- *State the limitation:* the review interface rendered text as Markdown, which removed multiplication and dollar
  signs from what the experts saw in 65 templates.

### 4 Evaluation Framework (1.5 pages)

Written in the register of the May section, with its equations.

**4.1 Scoring a Response (0.8 page).**

- *Preliminaries.* An instance is a question, a gold trace, a final answer and a set of milestones, the intermediate
  quantities the template computed with their units. A model returns a trace of steps and a final answer.
- *Final Answer Accuracy.* One target per asked quantity; the match rule as one equation (within 0.002 relative or
  one unit of the last displayed digit of either value, after unit conversion; LaTeX, fractions and superscripts
  read as numbers); the three-way verdict; FAC as the mean score. The two scoring rules in one sentence each: a
  response with no readable final answer scores 0, because not answering within the budget is the model's failure;
  a partial answer scores 0.5, the three-way verdict being what the experts validated.
- *Milestone Coverage.* Milestones derived per instance from the template's values; matched order-free within 0.5%
  under unit conversion. The milestones not found go to one judge with the question, the response and each
  milestone's name and value, never the gold; it answers reached, not needed or missing, and only "reached" is
  credited, because on validation the judge used "not needed" to excuse a quarter of fabricated values. MC as one
  equation. The judge decides 11% to 25% of required milestones. The judge is MiMo-V2.5-Pro, chosen because it
  shares no family with any evaluated model and detected planted defects as well as GPT-5 at a third of a cent per
  call.
- *Step-level diagnostics.* The arithmetic check and the judged step check in four sentences, as flags with
  precision; used in 5.2.

**4.2 Validation Against Expert Judgment (0.7 page).**

- *Expert study.* 60 problems from 15 templates across the five branches and three levels; five models; 300
  responses. 15 experts, three of the response's own branch per response, labelling every step (correct, alternative
  correct, incorrect with type and reason, or not a claim), every milestone, the final answer and the response's
  soundness; a reason on all 1,042 incorrect steps; a blind re-labelling round; a blind adjudication of the 272 split
  steps. Fleiss' kappa 0.781 between experts on steps, 0.828 within an expert, 0.880 on milestones, 0.966 on the
  answer.
- *Agreement with the experts* (Table 2). Seven designs compared: a step-matching evaluator with a three-judge panel,
  the same with a panel from outside the evaluated families, three process reward models, deterministic milestone
  coverage, coverage with the arithmetic check, coverage with the judge. The final-answer check 0.982 / 0.927 against
  the panel design's 0.760; milestone F1 0.923 and 0.958; the judge's validation; the arithmetic check 0.817 / 0.427
  inside correct-answer responses; the judged step check 0.707 / 0.603; the best reward model ranks steps well but
  three of four of its flags inside correct-answer responses are false.
- *Validation on the evaluated models.* Two experiments on this evaluation's own responses: (1) a domain expert read a
  fresh, stratified sample of the arithmetic check's flags, blind to the model, after the check had been refined on
  earlier samples: 171 of 189 decided flags are real slips, precision 0.905 (Wilson 0.854 to 0.939); (2) fifteen
  experts read 150 final-answer verdicts of the top five models, two readers each, and 100 judge decisions (C20).
  State the precisions per verdict category, never as one population rate.
- *Planted defects and the limit of verification.* 120 defects, 60 arithmetic and 60 conceptual, with 60 untouched
  controls; no deterministic check catches a misstated rule, judges catch about a third with no false alarm. The
  final answer nearly decides the experts' verdict (AUROC 0.974); 93 of 228 correct-answer responses carry an
  incorrect step, 175 of the 178 such steps arithmetic; coverage cannot see them. The design detects AUROC
  differences of 0.12 to 0.19 and no smaller.
- *Judge independence.* C19 in four sentences.

### 5 Experiments (2.85 pages)

**5.1 Setup (0.45 page).**

- *Models.* Eleven: eight open-weight (gpt-oss-20b, Gemma 4 26B-A4B, DeepSeek V4.1 Flash, Qwen3-235B-A22B-2507,
  GLM-5.3-Flash, GLM-5.3, Muse Glimmer 30B, Kimi K3) and three closed (GPT-5.4 mini, Gemini 3.1 Flash-Lite, Claude
  Sonnet 5), selected under one rule: no model that generated responses for the expert study and no judge. No
  flagship by a budget choice within the rule; two flagships run as anchors on a subset.
- *Decoding and prompt.* Each provider's default; a 32,768-token ceiling (16,384 for Muse); one response per instance,
  24,750 responses; the same zero-shot prompt throughout; closed book, no tools. Four models' endpoints returned no
  reasoning tokens (Gemma 4, Qwen3-235B-2507, GPT-5.4 mini, Gemini 3.1 Flash-Lite), so the models are not compared at
  equal reasoning effort; Appendix H states the setting per model. 246 responses (1.0%) ended without a readable
  answer, 243 of them empty at the ceiling, 153 on six templates, four of them iterative by construction.
- *Statistical protocol.* The template as the unit, because a template's 15 instances share a derivation; 95%
  intervals from resampling templates; within a family of comparisons Holm's correction, explained in one sentence;
  every null with the smallest difference the design detects at 80% power. Optionally one plain sentence that the
  confirmatory tests were specified before the models were run, which answers "which tests were planned"; keep or
  drop, never "pre-registered".

**5.2 Results (2.0 pages).**

- *Final answer accuracy.* Table 1 and Figure 3; the suggested wording of `RESULTS_PAPER_NOTES.md`: five models
  between 0.967 and 0.976 with no pair separating; GLM-5.3 at 0.948 owing most of its distance to 100 empty responses
  at the ceiling; the remaining five between 0.814 and 0.886; 32 of 55 pairs; 23 pairs below the detectable
  difference. A footnote with the two closed models' reasoning-on scores on the subset. Instance variance in one
  sentence.
- *Difficulty.* Figure 4. Every model lower on Advanced by 5.5 to 19.9 points; significant after correction for the
  four lowest-scoring models, three of which ran without reasoning tokens; for the top five 6 to 7 points with
  intervals mostly excluding zero but not surviving correction; the detectable gap 7 to 15 points; GLM-5.3's gap
  mostly an output-ceiling effect. Then the sensitivity: without the two Advanced chemical templates whose wording
  does not pin the answer (confirmed by three chemical experts), the gap holds for none. Both counts in one paragraph.
- *Milestone coverage.* The ordering with its intervals and both tests; coverage does not rise with the amount
  written; Figure 5, coverage on wrong answers against the chance floor; 13% to 26% of wrong answers are complete
  derivations to a wrong value; which signal points at a wrong answer, read through the diagnostics' precisions.
- *Step-level diagnostics.* The two flag rates with the calculations-per-response column and the precisions; why
  they do not rank models.
- *Consistency and depth.* C11 and C10; on a single-path template a partly solved template is a failure on the
  numbers, not the method.
- *Rewording.* The paraphrase result as a bound, in the suggested wording of `PARAPHRASE_PAPER_NOTES.md`: the 450 /
  316 / 277 funnel, ±5 points for ten of eleven, the one model whose interval reaches −6.7, tau against its noise
  floor, the experts' 12% rejection of script-passed pairs as evidence for the method; the decoding repeats; the
  sensitivity sentence on tolerance, shortcut and symbolic templates.
- *Reasoning effort, flagship anchors, the governing equations and a code tool.* C15, C16, C17 and C18 in four
  sentences, in the suggested wording of `RESULTS_PAPER_NOTES.md` ("For the retrieval objection"); the open-weight
  model's absence from the tool condition stated as a limit.

**5.3 Error Analysis (0.4 page).**

- *Expert reading of wrong answers.* C21: three readers, 160 wrong-answer responses of four contrastive models, the
  six categories plus "No error", kappa 0.930; what dominates per model and per level.
- *What the benchmark discriminates.* Stated as a finding: the evaluation set is solved at the top on final answers by
  models of September 2026; for the top five the remaining incorrect verdicts are mostly form (34 symbolic items on
  one template scored by their numbers), prescribed-precision misses on the 17 templates that require the digits the
  question prescribes, and the two under-specified chemical templates. The headroom is the Advanced tier, depth,
  prescribed-precision compliance and the process scores. The limit of verification as the methodological point: a
  response that computes correctly, misstates the rule it applies and lands the right answer passes every
  deterministic check, so process scores earn their place on the wrong answers and on the arithmetic behind right
  ones (`RESIDUAL_INCORRECT.md`; `docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` sections 5 and 6).

### 6 Conclusion (0.3 page)

One paragraph: what EngTrace is, what it found, what comes next (method memorisation, multi-modal artefacts, the
open-weight tier under tool use), without padding.

### Limitations (not counted; required)

From `RESULTS_PAPER_NOTES.md` "For the limitations" and `PARAPHRASE_PAPER_NOTES.md`:

- four models without extended reasoning at their providers' defaults;
- 1.0% of responses with no readable answer, concentrated on six templates where the score also measures finishing
  within the ceiling;
- coverage credits stated intermediates and is blind to a misstated rule; the arithmetic check's flags are 90% real
  but it finds under half of the slips the experts marked;
- the top five are not separable at this size (150 templates detect FAC differences of about 1 to 4 points and level
  gaps of 7 to 15);
- the evaluator was validated on five models not among the evaluated eleven and on 15 templates, with
  model-specific readings covering the final-answer check, the judge and the arithmetic check only;
- the two under-specified chemical templates and the 17 exact-digit templates; the certification's Markdown
  rendering;
- the paraphrase test covers 115 of 150 templates, tests wording not method, used one writer and a margin fixed
  after the point estimates;
- the open-book condition supplies equations only and the tool condition covers the two closed models;
- items are examination-style by construction; the experts' labels are not released (Appendix N).

### Ethics Statement

Synthetic data from referenced constants; no private data; models should not be deployed in safety-critical systems
on this evidence; experts' labels and personal data stay out of the release; the evaluation set and seed are released
at publication under the MIT licence; the experts were engineers who certified the templates and read the responses
(compensation, if any, stated).

## 6. The appendix, section by section

Letters follow the main text's order. "May" refers to the appendix of `docs/_ARR_May__EngTrace.pdf` (sections A to U,
Tables 2 to 18, Figures 6 to 8), which `current_overleaf_project/sections/7_appendix.tex` still holds. Every table is
generated by a committed script from the files named, never typed.

**A. Taxonomy and Content Selection.** Against May: A, B and D kept; one paragraph added.
- *Content and format:* the domain, area and pedagogical-significance prompts as listings, stated as used for the
  chemical, electrical and mechanical branches; a paragraph on how the civil and industrial domains and areas were
  chosen (by the templates' author, from the ABET, ASCE and IISE curricula and the branch textbooks); one sentence on
  who assigned difficulty levels.
- *Sources:* `7_appendix.tex` A, B, D; the authors.

**B. Source Texts and Reference Data.** Against May: C (Table 2) and E (Table 3) merged and extended to five branches.
- *Content and format:* a table of textbooks by domain, one primary text per domain as the May table does. For civil
  and industrial the colleague's source set (`civil_industrial_sources/references/public/`, local and gitignored;
  `MANIFEST.md` is the record) lists: Das and Sobhan, *Principles of Geotechnical Engineering*; Holtz, Kovacs and
  Sheahan, *An Introduction to Geotechnical Engineering*; Knappett and Craig, *Craig's Soil Mechanics*; Das and
  Sivakugan, *Principles of Foundation Engineering*; Hibbeler, *Structural Analysis*; Kassimali, *Structural
  Analysis*; Leet, Uang and Lanning, *Fundamentals of Structural Analysis*; Chow, Maidment and Mays, *Applied
  Hydrology*; Chin, *Water-Resources Engineering*; Sturm, *Open Channel Hydraulics*; Hillier and Lieberman,
  *Introduction to Operations Research*; Ross, *Introduction to Probability Models*; Taha, *Operations Research: An
  Introduction*; Nahmias and Olsen, *Production and Operations Analysis*; Silver, Pyke and Thomas, *Inventory and
  Production Management in Supply Chains*; Montgomery, *Introduction to Statistical Quality Control*; Grant and
  Leavenworth, *Statistical Quality Control*; Mitra, *Fundamentals of Quality Control and Improvement*. The colleague
  confirms the primary text for each of the six domains. A second table of reference data sources with their
  provenance classes: NIST WebBook and fluid properties, CODATA 2022, NIST-JANAF, NASA TR R-132, MIL-HDBK-5J with
  NIST SP 811, USDA Wood Handbook, refractiveindex.info, PubChem, NAVFAC DM-7.01 and 7.02, FHWA HDS-4 and HEC-22,
  NRCS TR-55, USGS WSP 2339, the AISC shapes database, MIL-STD-105E, the NIST/SEMATECH e-Handbook (the civil and
  industrial files are the same in the colleague's set and in `docs/references/`).
- *Sources:* `docs/references/README.md`, `MANIFEST.json`; the colleague's `MANIFEST.md`;
  `data/templates/branches/*/constants.py`.

**C. Template Construction and Examples.** Against May: F extended to five branches; G cut to one or two listings.
- *Content and format:* bullets per branch on how parameters are grounded and constrained (add civil and industrial);
  one or two template listings from different branches, shorter than May's three; a small table of structural
  variation (templates that change the governing equation, step count or computed quantities with a sampled
  parameter; single-path templates).
- *Sources:* `7_appendix.tex` F, G; `audit/template_audit_report.md`; `DIVERSITY.md`.

**D. Dataset Statistics and the Evaluation Set.** Against May: L (Table 9) replaced.
- *Content and format:* a table of templates by branch and level and by domain with area counts; instances by level;
  answer kinds with counts. A paragraph on how the 2,250 instances were drawn (private seed, round-robin over
  reasoning paths and answer forms, exclusions) and what is published now and at publication. A table or paragraph on
  question skeletons, reasoning paths and near-duplicate pairs.
- *Sources:* `manifest.jsonl`, `FREEZE.json`, `DIVERSITY.md`, D-114 to D-116.

**E. Template Certification.** Against May: H and I kept; J replaced; K (Tables 5 to 8) replaced.
- *Content and format:* E.1 automated integrity checks: what each check tests, 500 seeds, the outcome, the excused-line
  register, the display-tie census (prose and one table). E.2 the LLM screen: the prompt listing (May H, verbatim),
  the rubric table (May I), the judges and settings, the two passes with outcomes and agreement, the verified claims
  (table). E.3 the expert certification protocol: own-branch experts, hand check, the 20 planted defects (table:
  template, defect kind, description), opaque codes, timestamps, the interface (new screenshot or none). E.4 the
  outcome: plant detection, hand checks, agreement per branch (Fleiss, AC1, AC2), screen against experts
  (false-positive rate, MAD), round-by-round counts, the Markdown-rendering finding (tables).
- *Sources:* `template_annotation_23092026/README.md`, `layer0/gate_report.md`, `tie_census.md`, `screen/pass1/`,
  `pass2/stats.md`, `pass1_fixes.md`, `layer2/RESULTS*.md`, `CERTIFICATION.md`, `markdown_scan.md`.

**F. Scoring Details.** New; absorbs May N (prompt, adapted) and replaces May P.2.
- *Content and format:* F.1 the final-answer check: answer kinds (table with counts), how targets are derived, the
  match rule, partial and unreadable verdicts, and the scoring statements as lists (nine symbolic-answer templates
  scored by stated numbers; four templates verified on the last value; 84 multipart instances checked on fewer
  quantities than the question's labelled parts; the 17 exact-digit templates; the eleven symbolic, expression and
  classification templates on which the experts found the check conservative). F.2 milestones: derivation per
  instance, matching, the judge prompt (listing), verdict semantics, the judge's validation (table: 67 / 20 / 1 and
  0 / 66 / 22), judged fraction per model (table). F.3 the arithmetic check: the rule, worked examples, what it skips.
  F.4 the judged step check: the prompt (listing; May N's rubric adapted to a whole response), precision and recall.
  F.5 thresholds: the tables of `THRESHOLD_APPENDIX.md`.
- *Sources:* `full_run_28092026/answer.py`, `judge.py`, `router.py`, `THRESHOLD_APPENDIX.md`, D-138, D-145.

**G. Validation Against Expert Judgment.** Against May: M (Tables 10 to 12) replaced.
- *Content and format:* G.1 study design: problems, models, experts, the label scheme (table replacing Tables 10 and
  11), verification and adjudication rounds. G.2 expert agreement (table, overall and per branch). G.3 the evaluator
  comparison (table replacing Table 12: the seven designs; answer agreement, milestone precision, recall and F1,
  response-level AUROC with template intervals, cost). G.4 planted defects: design (table) and detection per checker
  and judge (table), the routing finding. G.5 judge independence: judge selection and the exposure table, the judge
  swap (table), the panel ablation and per-judge bias (tables). G.6 validation on the evaluated models: the
  arithmetic check's expert flag reading (protocol: stratified sample, model hidden, fresh sample after refinement;
  171 of 189 real), and the experts' readings of 150 final-answer verdicts and 100 judge decisions, per verdict
  category with inter-expert agreement (tables).
- *Sources:* pilot `PILOT_SUMMARY.md`, `RESULTS_X1.md`, `RESULTS_E1.md`, `RESULTS_E5.md`, `RESULTS_LOJO.md`,
  `RESULTS_ATTRIBUTION.md`, `JUDGE_SELECTION.md`; `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`,
  `JUDGE_SWAP.md`, `FLAG_REVIEW_3.md`, `FLAG_READER_INSTRUCTIONS.md`, `EXPERT_REQUEST.md`.

**H. Models and Inference Settings.** Against May: O (Table 13) replaced; P.1 kept.
- *Content and format:* a table from the decoding table: model, developer, open or closed, served identifier, ceiling,
  endpoints, median and 90th-percentile completion and reasoning tokens, share of responses with reasoning tokens,
  cost; the selection rule and the one model set aside; the prompt (May P.1, unchanged, as a listing); output reading
  in one sentence pointing to F.1; the statistical protocol in full (unit, resampling, correction, detectable
  differences, the confirmatory questions and when they were fixed, the labelled exploratory additions).
- *Sources:* `DECODING_TABLE.md`, `models.json`, `ANALYSIS_PLAN.md`, D-110, D-117, D-132.

**I. Full Results.** Against May: Q (Figures 7 and 8) replaced by tables with intervals.
- *Content and format:* I.1 Table 1 extended (FAC, strict FAC, unreadable share, SD within and between, MC with and
  without the judge, judged fraction, the two flag rates). I.2 pairwise comparisons (55 rows: difference, interval,
  corrected p, detectable). I.3 FAC and MC by branch, level and domain with intervals, and the branch-pair tests. I.4
  the level gap under its variations. I.5 by answer kind. I.6 the sensitivity table (and, if the authors adopt it,
  the sensitivity counting partial verdicts as correct on the eleven templates the experts flagged, D-185). I.7
  within-template consistency and instance variance. I.8 failure against depth. I.9 decoding repeats. I.10 tokens
  against score. I.11 coverage details (wrong-answer coverage against the floor, verdict against coverage, coverage
  against verbosity). I.12 diagnostics details (flag rates with intervals, first-flag position, what points at a
  wrong answer, precision by error type).
- *Sources:* `results/RESULTS.md`, `results.json`, `per_template.csv`, pilot `RESULTS_ATTRIBUTION.md`.

**J. Robustness to Rewording.** New.
- *Content and format:* design (subset, writer, the rewriting rules), the scripted checks, the expert check and its
  rejection reasons with counts, survival by branch and the list of lost templates, the results table (FAC, MC with
  and without the judge; intervals; corrected p), the bounds table (90% intervals against ±5), tau against its noise
  floor, the comparison with decoding repeats.
- *Sources:* `PARAPHRASE.md`, `PARAPHRASE_REVIEW.md`, `PARAPHRASE_AUDIT.md`, `RESULTS.md` Q5, D-157 to D-165.

**K. Additional Conditions.** New. Four conditions, all run on the 450-instance subset.
- *Content and format:* one design paragraph each for reasoning effort (two closed models), the flagship anchors (two
  flagships, one also with reasoning on), the governing equations in the prompt (the corrected build only: 135
  templates, 405 instances; how the equation blocks were built from the templates' own stated equations; the
  templates left out) and the Python tool (the tool's definition, the isolated interpreter, the call limit; the two
  closed models; why the open-weight model could not be served). The paired tables as `RESULTS.md` prints them, the
  anchors beside the eleven on the same instances, the per-condition decoding tables, and the tool-use columns
  (share of responses calling the tool, calls per response, score with and without a call).
- *Sources:* `RESULTS.md` "C1 and C4" (the `openbook2` and `tool` rows), "C3"; `DECODING_TABLE_*.md`;
  `OPENBOOK_SURVEY_2.md`; D-180 to D-184.

**L. Error Analysis.** Against May: R kept with two added categories; S, T, U replaced; Table 15 redrawn or dropped.
- *Content and format:* L.1 taxonomy (Table 14's definitions and the decision hierarchy, with "No error" and
  "Incomplete" added). L.2 sample and protocol (160 wrong-answer responses, four models, levels in proportion, three
  readers). L.3 agreement (Fleiss' kappa overall and per model). L.4 results by model and level (table). L.5
  representative responses, if new examples are drawn with the readers' consent. L.6 templates whose wording does not
  pin the answer: the two chemical templates with the experts' answers, the 17 exact-digit templates, the near-miss
  distance table, the six templates whose answers can be read off the question.
- *Sources:* `7_appendix.tex` R; `EXPERT_REQUEST.md` B2, B4; `RESIDUAL_INCORRECT.md`; `SHORTCUT_AUDIT.md`.

**M. Extended Related Work.** New.
- *Content and format:* Table A (the fifteen closest benchmarks) and the further-related-work paragraphs.
- *Sources:* `RELATED_WORK_v3.md` appendix.

**N. Reproducibility and Release.** New.
- *Content and format:* what is released (templates, code, the evaluation set with the seed at publication, responses
  and scores, the paraphrases' hashes); what is withheld and why (names; the experts' labels, unless the authors
  decide to release an anonymised file of codes and majority labels with no per-rater rows); the anonymised archive's
  contents for review (never the repository, which carries author names and the working record); the script behind
  each table.
- *Sources:* this plan; `full_run_28092026/README.md`.

**Dropped from May without replacement:** the textual-quality metrics (BERTScore, ROUGE) and their columns; the
ablation of evaluator configurations by rank correlation (Table 12); the Tribunal examples (Table 5), kappa and alpha
tables (Tables 6 and 7) and AI-human MAD table (Table 8) of the old certification; the annotation-interface
screenshot (Figure 6) unless replaced; the old sampling-strategy equations (S.1); the 11-model error tables (Tables
16 to 18); the word "Tribunal".

## 7. Figures and tables

| Item | Content | Source | State |
|---|---|---|---|
| Figure 1 | A template, one instance, its gold trace, and the milestones and calculations the evaluator checks | a certified template; `milestones.py` | revise `figs/engtrace-example-template.pdf` (overlay), or a civil or industrial template to show a new branch |
| Figure 2 | Five-branch taxonomy | `figures_oct_12/engtrace-overview-5branch.pdf` | exists; copy to `figs/` |
| Figure 3 | FAC and MC per model with intervals, tiers marked | `results/results.json` Q1, Q3 | script to write |
| Figure 4 | FAC on Easy, Intermediate and Advanced per model with intervals, the corrected verdict marked | `results.json` Q2 and the level table | script to write; may move to Appendix I |
| Figure 5 | Coverage on wrong answers against the chance floor, per model | `results.json` Q3 | script to write; may move to Appendix I |
| Figure 6 | Change in FAC under rewording, per model, 90% intervals and the ±5 band | `results.json` Q5 | script to write; may move to Appendix J |
| Table 1 | FAC, MC, compact letter display | `results.json` Q1, Q3 | script to write |
| Table 2 | Validation against the experts: each component's agreement, precision, recall; planted-defect detection | `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, RESULTS_X1 | collate by script |

Figures 1 to 3 and Tables 1 and 2 are the must-haves in the main text. Every figure and table is produced by a
committed script; figures are drawn to read in greyscale; the ACL column is 3.03 inches wide.

## 8. Writing rules

1. **Simple, straightforward language**, consistent and understandable as a whole. Short sentences; one idea each;
   the same word for the same thing throughout.
2. **No unnecessarily complex jargon or technical terms.** Where a statistical or engineering term is needed, define
   it once in plain words and use it the same way everywhere.
3. **Few terms, used consistently.** The paper invents nothing beyond the two metrics and the two diagnostics of
   section 4 and the terms it defines once (template, instance, gold trace, milestone, evaluation set, condition).
   Every new term is a cost for the reviewers.
4. **Do not overstate.** This is the most important rule. Say what the data show, with the interval, the unit and
   the limit beside it; section 3.2 lists the overstatements already caught, and every sentence is checked against it.
5. **No internal vocabulary in the paper. Strictly forbidden:** evaluator labels (E0 to E5, E5-strict), question and
   item labels (Q1 to Q5, A1 to A11, B1 to B4, C1 to C5), decision numbers, commit hashes, file, script and variable
   names, and the working vocabulary of the record. The table gives the replacement for each; the proofreading pass
   greps for the left column.

   | Internal term | Paper wording |
   |---|---|
   | E3, E4, E5, E5-strict, "coverage (stage)" | Milestone Coverage; "without the judge" and "with the judge" where both are shown |
   | digit rule, E4 flags, flag review | the arithmetic check; arithmetic flags; the expert reading of flags |
   | router, step router, rule C | the judged step check; judged step flags |
   | residual judge, judge calls | the judge (first use: "a judge on the milestones the matcher does not find") |
   | Tribunal, AI Tribunal, panel | a panel of LLM judges (only when describing the compared design) |
   | Layer 0 / gate / T1 closure, T3, T4, T8 | automated integrity checks (each described in words) |
   | Layer 1 / screen pass 1, pass 2 | the LLM screen; its two passes |
   | Layer 2 / plants / kits / app / workbook | expert certification; planted defects; the annotation interface |
   | pilot, slice, frozen slice | the expert study; the validation set |
   | pool, freeze, manifest, seed commitment | the evaluation set; the published hashes; a commitment to the seed |
   | roster, roster rule | the evaluated models; the selection rule |
   | run, main run, arm, variant, subsample, openbook2 | the evaluation; condition; the 450-instance subset; the open-book condition |
   | trace (where it means a model output) | response (or "trace" only for the reasoning text, defined once) |
   | unusable, empty, capped | no readable final answer; empty at the output ceiling |
   | fully solved, FSR, answer score | strict Final Answer Accuracy; Final Answer Accuracy |
   | REACHED / NOT_NEEDED / MISSING | reached, not needed, missing (in words) |
   | store, stage rows, score store, CONFIG | not mentioned |
   | owner, supervisor, agent, session, colleague's "Stage B/C" | not mentioned |
   | pre-registered | "specified before the models were run" (one optional sentence) |
   | D-xxx, `RESULTS.md`, `analyze.py`, any path | not mentioned; Appendix N names the scripts behind the tables once |

6. **One voice, present tense, no history.** The paper describes what EngTrace is; "we" build, certify, score,
   validate. No "previously", "originally", "in response to".
7. **Every number from a committed script**, regenerated from `results.json`, the manifest or the validation files;
   never copied from a Markdown table by hand.
8. **The unit of analysis is named once and kept**: the template; every interval resamples templates; every family of
   tests is corrected; every null carries its detectable difference.
9. **Tiers, not an order**, wherever intervals overlap; "not separable at this size", never "equal".
10. **Flags are flags; coverage is coverage.** The two naming rules of section 4.
11. **Say the setting with the number**: provider defaults, the four models without reasoning tokens, the ceiling,
    the unreadable-answer rule, wherever a closed-against-open or a level sentence appears.
12. **Name the limits where the claim is made**: the judged fraction beside coverage, the precision beside a flag
    rate, the two chemical templates beside the level gap, the 115 templates beside the paraphrase bound, the two
    closed models beside the tool result.
13. **Limits that must appear in the paper**: the Markdown rendering during certification; the two under-specified
    chemical templates; the 17 exact-digit templates; the four no-reasoning models; the validation population
    against the evaluated models; the open-weight tier's absence from the tool condition. The history of how the
    final-answer check was refined during the work (its corrected answer forms, the verdicts that moved, the one
    known wrong credit) is not paper material: the paper describes the check as it is and reports its validation,
    including the expert reading of the arithmetic flags; the history goes to the revision letter and stays in the
    repository record.
14. **Length discipline**: cut evidence to the appendix before cutting a limit from the main text.
