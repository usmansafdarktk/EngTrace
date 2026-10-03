# EngTrace, October 2026 submission: the paper as it stands, and how to write it

Written 2026-10-03 for the authors, nine days before the ARR October deadline (12 October 2026, anywhere on
Earth); revised the same evening after the authors' review (decision D-186). It reads the May 2026 submission
(`docs/_ARR_May__EngTrace.pdf`), both rebuttals and meta-reviews, the current Overleaf sources
(`current_overleaf_project/`), and the September and October record: the template rework
(`docs/re-implementation-sep/`), the certification (`template_annotation_23092026/`), the evaluator study
(`evaluator_pilot_17092026/`), the full evaluation (`full_run_28092026/`, `results/RESULTS.md` at commit `026210a`),
the revised Related Work (`docs/related_work_oct2026/`), the five-branch figure (`figures_oct_12/`) and the experts'
readings scored today (`full_run_28092026/EXPERT_REQUEST.md`). Every number names the file that prints it.

Two companions hold what this file no longer does: `docs/REVISION_LETTER_NOTES.md` (what changed since May and
why, for the letter only) and `docs/OVERLEAF_WORKFLOW.md` (how the Overleaf project is kept in sync).
`full_run_28092026/RESULTS_PAPER_NOTES.md` and `PARAPHRASE_PAPER_NOTES.md` say, number by number, what each result
supports; this file says what the paper is and where everything goes.

**The rule that shapes everything.** The paper is written for readers who have never seen EngTrace. It describes
the benchmark, the evaluator and the results as they are, with their evidence and their limits. No "previously",
no "in the earlier version", no "as promised". And no internal vocabulary: section 8 makes that rule binding and
gives the replacement for every internal term.

---

## 1. The paper in one paragraph

EngTrace is a benchmark of 150 symbolic templates across five engineering branches (chemical, civil, electrical,
industrial, mechanical), 15 domains and 42 areas. Each template is executable code: it samples physically grounded
parameters from referenced constants, computes the answer, and writes the gold reasoning trace from the same
computation, so every intermediate quantity exists with its unit before any prose does. Templates reach the
benchmark through automated integrity checks, an LLM screen whose judges share no family with any evaluated model,
and certification by three experts of the template's own branch, who caught 60 of 60 planted defects and rejected 53
templates across rounds until all 150 were approved unanimously. The evaluation set (2,250 instances, 15 per
template) is drawn from a private seed, so no evaluated instance was public when the models ran. A response is scored
deterministically wherever a check exists (a three-way final-answer check over six answer kinds, coverage of the gold
derivation's intermediate quantities, an arithmetic check of every displayed calculation), and a judge from a family
outside the evaluated models decides only what those checks leave open. The evaluator was chosen by comparing seven
designs against 15 experts' step-level labels on 300 responses and 120 planted defects, which also fixed the limit of
verification: no deterministic check detects a misstated rule behind a correct answer, and LLM judges catch about a
third. Eleven open and closed models of September 2026 score 0.81 to 0.98 on the final answer; the top five are
inseparable at this size; every model scores lower on Advanced templates; milestone coverage orders the models
differently from the answer; on wrong answers models still reach most of the derivation; and expert-checked
paraphrases change no model's score by more than five points for ten of the eleven.

## 2. Title, abstract and contributions

**Title.** Keep *EngTrace: A Symbolic Benchmark for Verifiable Process Supervision of Engineering Reasoning*. The
mechanism now matches the title. A shorter form, if wanted: *EngTrace: Verifiable Process Evaluation of Engineering
Reasoning from Symbolic Templates*.

**Abstract, draft (about 250 words; trim to taste).**

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
> of eleven models.

**Contributions, draft.** Three, each owned by an active verb and each carrying its evidence; nothing about how the
analysis was organised internally.

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
   rewording, and reasoning-effort, flagship and open-book conditions.

## 3. The claims the data support, and the claims not to make

Each claim is a sentence the paper may make, in the paper's own words, with its source. The "do not" list is
binding; it collects every overstatement the record has already caught (`RESULTS_PAPER_NOTES.md`, both "What not
to state" sections; `PARAPHRASE_PAPER_NOTES.md`; D-111, D-146, D-171, D-180 to D-183).

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
| C17 | Supplying each template's governing equations with the question (open book, as run) changed Final Answer Accuracy by +0.005, −0.007 and +0.041 for the three models tried, none holding after correction, while coverage rose by 0.03 to 0.06; the equation blocks had defects (D-183) and a corrected condition was running as this was written | "C1 and C4"; D-183 |
| C18 | Judge independence: the judge shares no family with any evaluated model; on a 220-response sample a second judge from a third family moves no model's coverage beyond its interval (−0.017 to +0.020); on a panel-of-judges design replayed offline, dropping a family's own judge changes its score by at most 0.006, the same as a placebo drop, and every judge is lenient on every family alike | `JUDGE_SWAP.md`; pilot `RESULTS_LOJO.md`; D-174, D-181 |
| C19 | On these models, two experts read 150 of the top five models' final-answer verdicts (half incorrect or partial by design): correct verdicts confirmed in 142 of 150 readings, incorrect verdicts in 83 of 98 (12 expert-correct), partial verdicts called fully correct in 34 of 52, so the half-point is conservative; the judge's REACHED verdicts confirmed in 85 of 100 readings and its MISSING verdicts in 79 of 100 | `EXPERT_REQUEST.md` B1, B3 |
| C20 | Three readers labelled 160 wrong-answer responses of four contrastive models with the six-category taxonomy plus "No error" (Fleiss' kappa 0.930): for Claude Sonnet 5, 26 of 40 item majorities are "No error"; for gpt-oss-20b, Gemma 4 and GPT-5.4 mini, Calculation Error leads (17, 28 and 28 of 40), then Formula/Principle; Easy wrong answers are almost all calculation (69 of 78 readings) | `EXPERT_REQUEST.md` B2 |

### 3.2 Do not state

- Not "model X is the best": tiers, not an order, at the top; the top five are within 0.9 points and no pair separates.
- Not "trade-off", "decoupling", "trace fidelity" or "reasoning quality" for Milestone Coverage: it is coverage of stated intermediates; no test of a trade-off was planned and the intervals overlap.
- Not "x% of correct answers rest on flawed reasoning" from either diagnostic: flag rates with known precision (0.905; about 0.70) and recall below one half.
- Not a ranking of models by arithmetic reliability: the check reads 2.8 to 9.8 calculations per response depending on style.
- Not "the drop on Advanced problems holds" without both counts (4 of 11 as scored; 0 of 11 without the two chemical templates) and the per-model intervals.
- Not "closed models underperform open ones" or any comparison at equal effort: four models ran without reasoning tokens.
- Not "chemical engineering is hardest" or any branch order: one branch pair in the whole evaluation holds (gpt-oss-20b, electrical above civil); domain means are descriptive.
- Not "2,250 unique problems": 2,250 instances, all distinct questions, from 150 templates.
- Not "robust to paraphrasing" or "no contamination": a bound of ±5 points for ten models on wording memorisation; method memorisation is not tested.
- Not "the ranking is stable under rewording": tau 0.673 sits at the lower edge of the sampling-noise distribution; the tiers hold.
- Not "flagships do no better" and not "solved at the top" without the Advanced means beside it (0.87 to 0.95 for the anchors, 0.92 to 0.95 for the top five, wide intervals).
- Not "retrieval would not help": the open-book condition supplies equations only, had defective blocks as run, and the tool condition is pending.
- Not "no self-preference" as evidence that LLM judges are accurate: they are uniformly lenient.
- Not templates as a remedy for benchmark saturation (Akhtar et al. 2026 find none); the private seed answers contamination, the step-level checks answer the question a saturated answer score cannot.
- Not the expert-study agreement figures as validation of answer forms the study's templates do not contain (pi-fractions, two-unit lines): say the experts did not arbitrate those.
- Not "no existing benchmark verifies the process"; the claim is the combination, hedged as in `RELATED_WORK_v3.md`.
- Not any number beside a May number; the paper describes the evaluator it uses.

## 4. Two metrics and two diagnostics

The paper names two scores and two diagnostics, and nothing else. Level gaps, depth, consistency and the
paraphrase change are presented as analyses of the two scores ("Final Answer Accuracy by level", "change in Final
Answer Accuracy under rewording"), not as metrics with names. Every score carries a 95% interval that resamples
templates; every null carries the smallest difference the design detects.

| Role | Name | Definition, as the paper gives it |
|---|---|---|
| Final answer | **Final Answer Accuracy (FAC)** | Mean per-instance score: 1 when every asked quantity matches the gold within the tolerance, 0.5 when some but not all do, 0 otherwise or when no final answer is readable. The tolerance is 0.002 relative or one unit of the last displayed digit, after unit conversion; targets are the quantities the template computed; six answer kinds (scalar, vector, array, symbolic, classification, multipart). |
| Reasoning | **Milestone Coverage (MC)** | Share of the gold derivation's intermediate quantities (milestones) that the response reaches: matched deterministically within 0.5% under unit conversion, order-free; the milestones the matcher does not find go to a judge from outside the evaluated families, and only its "reached" verdict is credited. Measures progress through the derivation, not the absence of error. |
| Diagnostic | **arithmetic flags** | Share of correct-answer responses with at least one displayed calculation that, recomputed from its printed operands, is not a correct rounding at the printed precision; reported with the calculations read per response and the measured precision (0.905). |
| Diagnostic | **judged step flags** | Share of correct-answer responses in which the judge flags a step the arithmetic check did not; reported with its precision and recall on the experts' labels. |

Two choices inside this scheme, adopted here as recommendations and open to reversal. FAC keeps partial credit
because it is the headline the analysis specified, because the experts called the check's partial verdicts fully
correct in 34 of 52 readings so half a point is the conservative middle, and because strict accuracy changes no
conclusion (54 of 55 pairs agree) and sits in Appendix I. "Milestone Coverage" is preferred over "Derivation
Coverage" because the word milestone is defined once and then does work in every table; the choice is cosmetic.

Table 1 is therefore: model, FAC with its interval, MC with its interval, and a compact letter display (models
sharing a letter do not differ after correction), so the tiers are visible without another column. A footnote marks
GLM-5.3's 100 responses with no readable answer. The strict FAC, the share of responses with no readable answer, the
per-template SD, the judged fraction and the two diagnostics go to Appendix I.

## 5. Sections

ARR long paper: eight pages of content, unlimited references and appendices; the review version must stand on its
own. Few numbered subsections; most structure is paragraphs with bold run-in headings (`\paragraph{}`), as the
current Overleaf sources already do. The budget sums to eight pages.

### 1 Introduction (1.0 page)

Open with the gap in the form the revised Related Work supports: engineering reasoning is a derivation whose
intermediate quantities can be checked, yet most engineering benchmarks score only final answers, and those that
assess steps do so with LLM judges or human graders, or within one domain. Then EngTrace in one paragraph (section 1
above), Figure 1 (a template, one instance, its gold trace, and the milestones and calculations the evaluator checks
in it), the scope sentence (a depth-oriented evaluation of a principled subset of engineering; items are
examination-style by construction, which is what verifiability needs and which leaves out ill-posed problems and
real artefacts), and the three contributions. `docs/related_work_oct2026/CHANGES.md` section 2 has the replacement
sentences for every Introduction sentence that cites a prior benchmark; use them as written. No "first", no
"no existing benchmark verifies this process", no "safety-critical" beyond one motivating clause.

### 2 Related Work (0.75 page)

The main text of `docs/related_work_oct2026/RELATED_WORK_v3.md` as five paragraphs (`related_work_v3.tex`;
bibliography `related_work_sources.bib`, 100 entries); the file names the cut order from 704 words to about 0.75
page. Table A and the further-related-work paragraphs go to Appendix M. Anonymity: FinChain in the third person;
the EngTrace preprint not cited as the authors' own; ERI, Sci-rho and LPDS cited without saying they cite EngTrace;
whether the Limitations acknowledge LPDS's finding is the authors' call (`CHANGES.md` section 4). Venues to confirm:
EMNLP 2026 for Huang et al. and Dlugosz et al.

### 3 The EngTrace Benchmark (1.6 pages)

**3.1 Taxonomy and Scope (0.3).** Five branches, 15 domains, 42 areas, 30 templates per branch, 10 per domain
except chemical (10 / 12 / 8); Figure 2, the five-branch overview (`figures_oct_12/engtrace-overview-5branch.pdf`,
to be committed and copied into `figs/`). Paragraphs: *Domains and areas* (curricular standards: ABET and the
societies, ASME, IEEE, AIChE for the original three, ASCE and IISE for civil and industrial; the textbooks of
Appendix B; the expert-persona check of the domain and area lists for the original three branches; the procedure
used for the two new branches, stated as it was, section 9 decision 3). *Difficulty* (Easy, Intermediate, Advanced
by conceptual complexity, mathematical sophistication and procedural depth; 58 / 58 / 34 templates).

**3.2 Templates and Instances (0.6).** Paragraphs: *Template structure* (a seeded function that samples parameters,
often inside a rejection loop enforcing physical validity, computes every quantity, binds each displayed value
through its display before using it, asserts physical invariants, and emits the question and the gold trace with
numbered steps and an answer line; six answer kinds). *Grounding* (every constant table carries a provenance class
and resolves to a referenced source: NIST, CODATA, JANAF, NASA, military and federal handbooks, the textbooks;
Appendix B). *Structural variation* (72 of 150 templates change the governing equation, the step count or the
computed quantities with a sampled parameter; the evaluation set holds a second reasoning path for every template
that can produce one within 500 draws; 58 templates follow one path; `template_audit_report.md`, `DIVERSITY.md`).
*The evaluation set* (15 instances per template from a private 128-bit seed, chosen round-robin across each
template's reasoning paths and answer forms, excluding repeated questions and display ties; the instances' hashes
and a commitment to the seed are published, the seed and the text at publication).

**3.3 Certification (0.7).** Paragraphs: *Authorship* (all 150 templates were written by people with engineering
training and grounded in textbooks: the chemical, electrical and mechanical templates by the team's domain experts,
the civil and industrial templates by a colleague of the authors; exact wording in section 9 decision 2). *Automated
integrity checks* (at 500 seeds per template: every printed calculation closes from its printed operands at the
printed precision, generation is deterministic across processes, the output format holds, every seed generates; 150
of 150 pass, 54 templates edited to pass; the display-tie census). *LLM screen* (the rubric of Appendix E; three
judges from no evaluated family, Grok 4.6, MiniMax M3 and MiMo-V2.5-Pro; two passes; 147 pass, 3 controversial, 0
critical on the shipped templates; agreement AC1 0.93 on the flag; every flag verified before any edit, 34 of 45
claims confirmed). *Expert certification* (three experts of the template's own branch, 34 items each: the branch's
30 templates and 4 planted defects, one of each kind: wrong constant, wrong unit conversion, flipped sign, a printed
step that does not follow; a hand check before the solution is shown; opaque codes, timestamps). *Outcome* (60 of 60
planted defects rejected; 459 of 504 hand checks within 1%, 44 of the 45 mismatches ending in a rejection; AC1 0.914
on approve/reject with Fleiss' kappa 0.589 beside it, deflated by prevalence; AC2 0.947 / 0.980 / 0.954 on the three
quality scores; 53 rejections over the rounds, each verified against its instance and fixed or, in two cases,
answered; every civil and industrial expert approved all 30 of their templates in round 1 and the 22 round-1
rejections fell on the original three branches; rounds 2 to 4 re-certified the changed templates; 150 of 150
certified unanimously on code that still regenerates the reviewed instances byte for byte; the screen's
false-positive rate against the experts, 10.2%, is why the screen is a reader and the experts the gate). State the
limitation: the review interface rendered text as Markdown, which removed multiplication and dollar signs from what
the experts saw in 65 templates.

### 4 Evaluation Framework (1.5 pages)

**4.1 Scoring a Response (0.8).** Paragraphs, in the register of the May section with its equations:
*Preliminaries* (an instance is a question, a gold trace, a final answer and a set of milestones, the intermediate
quantities the template computed with their units; a model returns a trace of steps and a final answer). *Final
Answer Accuracy* (one target per asked quantity; the match rule as one equation: within 0.002 relative or one unit
of the last displayed digit of either value, after unit conversion; LaTeX, fractions and superscripts read as
numbers; the three-way verdict; FAC as the mean score; the two scoring rules in one sentence each: a response with no
readable final answer scores 0 because not answering within the budget is the model's failure; a partial answer
scores 0.5, the three-way verdict being what the experts validated). *Milestone Coverage* (milestones derived per
instance from the template's values; matched order-free within 0.5% under unit conversion; the milestones not found
go to one judge with the question, the response and each milestone's name and value, never the gold, who answers
reached, not needed or missing; only "reached" is credited because on validation the judge used "not needed" to
excuse a quarter of fabricated values; MC as one equation; the judge decides 11% to 25% of required milestones; the
judge is MiMo-V2.5-Pro, chosen because it shares no family with any evaluated model and detected planted defects as
well as GPT-5 at a third of a cent per call). *Step-level diagnostics* (the arithmetic check and the judged step
check in four sentences, as flags with precision; used in 5.2).

**4.2 Validation Against Expert Judgment (0.7).** Paragraphs: *Expert study* (60 problems from 15 templates across
the five branches and three levels, five models, 300 responses; 15 experts, three of the response's own branch per
response, labelling every step as correct, alternative correct, incorrect with type and reason, or not a claim,
every milestone, the final answer and the response's soundness; a reason on all 1,042 incorrect steps, a blind
re-labelling round and a blind adjudication of the 272 split steps; Fleiss' kappa 0.781 between experts on steps,
0.828 within an expert, 0.880 on milestones, 0.966 on the answer). *Agreement with the experts* (Table 2: seven
designs compared, a step-matching evaluator with a three-judge panel, the same with a panel from outside the
evaluated families, three process reward models, deterministic milestone coverage, coverage with the arithmetic
check, coverage with the judge; the final-answer check 0.982 / 0.927 against the panel design's 0.760; milestone F1
0.923 and 0.958; the judge's validation; the arithmetic check 0.817 / 0.427 inside correct-answer responses; the
judged step check 0.707 / 0.603; the best reward model ranks steps well but three of four of its flags inside
correct-answer responses are false; and on these models, the arithmetic check's precision 0.905 on 189 flags read by
a domain expert, the experts' readings of the top models' verdicts, C19). *Planted defects and the limit of
verification* (120 defects, 60 arithmetic and 60 conceptual, with 60 untouched controls; no deterministic check
catches a misstated rule, judges catch about a third with no false alarm; the final answer nearly decides the
experts' verdict, AUROC 0.974; 93 of 228 correct-answer responses carry an incorrect step, 175 of the 178 such steps
arithmetic; coverage cannot see them; the design detects AUROC differences of 0.12 to 0.19 and no smaller). *Judge
independence* (C18 in four sentences). One sentence on corrections: the final-answer check was corrected after the
evaluation on answer forms the study's templates did not contain, each correction measured on the gold solutions,
the labelled responses and every response before adoption, none moving an expert-labelled verdict, the last moving
238 of 24,750 verdicts, 234 to correct, every one read; one wrong credit is known and stated (Appendix G).

### 5 Experiments (2.85 pages)

**5.1 Setup (0.45).** Paragraphs: *Models* (eleven: eight open-weight, gpt-oss-20b, Gemma 4 26B-A4B, DeepSeek V4.1
Flash, Qwen3-235B-A22B-2507, GLM-5.3-Flash, GLM-5.3, Muse Glimmer 30B, Kimi K3; three closed, GPT-5.4 mini, Gemini
3.1 Flash-Lite, Claude Sonnet 5; selected under one rule: no model that generated responses for the expert study and
no judge; no flagship by a budget choice within the rule, two run as anchors on a subset). *Decoding and prompt*
(each provider's default; a 32,768-token ceiling, 16,384 for Muse; one response per instance, 24,750 responses; the
same zero-shot prompt throughout, closed book, no tools; four models' endpoints returned no reasoning tokens, Gemma 4,
Qwen3-235B-2507, GPT-5.4 mini and Gemini 3.1 Flash-Lite, so the models are not compared at equal reasoning effort and
Appendix H states the setting per model; 246 responses, 1.0%, ended without a readable answer, 243 of them empty at
the ceiling, 153 on six templates, four of them iterative by construction). *Statistical protocol* (the template as
the unit, because a template's 15 instances share a derivation; 95% intervals from resampling templates; within a
family of comparisons Holm's correction, explained in one sentence; every null with the smallest difference the
design detects at 80% power; optionally, one plain sentence that the confirmatory tests were specified before the
models were run, which answers "which tests were planned"; keep or drop, never "pre-registered").

**5.2 Results (2.0).** Paragraphs: *Final answer accuracy* (Table 1 and Figure 3; the suggested wording of
`RESULTS_PAPER_NOTES.md`: five models between 0.967 and 0.976 with no pair separating; GLM-5.3 at 0.948 owing most of
its distance to 100 empty responses at the ceiling; the remaining five between 0.814 and 0.886; 32 of 55 pairs; 23
pairs below the detectable difference; a footnote with the two closed models' reasoning-on scores on the subset;
instance variance in one sentence). *Difficulty* (Figure 4; every model lower on Advanced by 5.5 to 19.9 points;
significant after correction for the four lowest-scoring models, three of which ran without reasoning tokens; for
the top five 6 to 7 points with intervals mostly excluding zero but not surviving correction; the detectable gap 7
to 15 points; GLM-5.3's gap mostly an output-ceiling effect; without the two Advanced chemical templates whose
wording does not pin the answer, confirmed by three chemical experts, the gap holds for none; both counts in one
paragraph). *Milestone coverage* (the ordering with its intervals and both tests; coverage does not rise with the
amount written; Figure 5, coverage on wrong answers against the chance floor; 13% to 26% of wrong answers are
complete derivations to a wrong value; which signal points at a wrong answer, read through the diagnostics'
precisions). *Step-level diagnostics* (the two flag rates with the calculations-per-response column and the
precisions; why they do not rank models). *Consistency and depth* (C11; C10; on a single-path template a partly
solved template is a failure on the numbers, not the method). *Rewording* (the paraphrase result as a bound, in the
suggested wording of `PARAPHRASE_PAPER_NOTES.md`: the 450 / 316 / 277 funnel, ±5 points for ten of eleven, the one
model whose interval reaches −6.7, tau against its noise floor, the experts' 12% rejection of script-passed pairs as
evidence for the method; the decoding repeats; the sensitivity sentence on tolerance, shortcut and symbolic
templates). *Reasoning effort, flagship anchors and open book* (C15, C16, C17 in three sentences, the open-book
sentence in the form the record allows at submission).

**5.3 Error Analysis (0.4).** Paragraphs: *Expert reading of wrong answers* (C20: three readers, 160 wrong-answer
responses of four contrastive models, the six categories plus "No error", kappa 0.930; what dominates per model and
per level). *What the benchmark discriminates* (stated as a finding: the evaluation set is solved at the top on
final answers by models of September 2026; for the top five the remaining incorrect verdicts are mostly form, 34
symbolic items on one template scored by their numbers, prescribed-precision misses on the 17 templates that
require the digits the question prescribes, and the two under-specified chemical templates; the headroom is the
Advanced tier, depth, prescribed-precision compliance and the process scores; the limit of verification as the
methodological point: a response that computes correctly, misstates the rule it applies and lands the right answer
passes every deterministic check, so process scores earn their place on the wrong answers and on the arithmetic
behind right ones; `RESIDUAL_INCORRECT.md`, `docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` sections 5 and 6).

### 6 Conclusion (0.3 page)

One paragraph: what EngTrace is, what it found, what comes next (the tool condition and the corrected open-book
condition if not yet reported, method memorisation, multi-modal artefacts), without padding.

### Limitations (not counted; required)

From `RESULTS_PAPER_NOTES.md` "For the limitations" and `PARAPHRASE_PAPER_NOTES.md`: four models without extended
reasoning at their providers' defaults; 1.0% of responses with no readable answer, concentrated on six templates where
the score also measures finishing within the ceiling; coverage credits stated intermediates and is blind to a
misstated rule; the arithmetic check's flags are 90% real but it finds under half of the slips the experts marked; the
top five are not separable at this size (150 templates detect FAC differences of about 1 to 4 points and level gaps of
7 to 15); the evaluator was validated on five models not among the evaluated eleven and on 15 templates, with
model-specific readings covering the final-answer check, the judge and the arithmetic check only; the two
under-specified chemical templates and the 17 exact-digit templates; the certification's Markdown rendering; the
paraphrase test covers 115 of 150 templates, tests wording not method, used one writer and a margin fixed after the
point estimates; items are examination-style by construction; the experts' labels are not released (section 9
decision 5); the open-book condition's defects and the tool condition's state, if still pending.

### Ethics Statement

Synthetic data from referenced constants; no private data; models should not be deployed in safety-critical systems
on this evidence; experts' labels and personal data stay out of the release; the evaluation set and seed are released
at publication under the MIT licence; the experts were engineers who certified the templates and read the responses
(compensation, if any, stated).

## 6. The appendix, section by section

Letters follow the main text's order. "May" refers to the appendix of `docs/_ARR_May__EngTrace.pdf` (sections A to U,
Tables 2 to 18, Figures 6 to 8), which `current_overleaf_project/sections/7_appendix.tex` still holds. Every table is
generated by a committed script from the files named, never typed.

| New | Title | Against the May appendix | Content and format | Sources |
|---|---|---|---|---|
| A | Taxonomy and Content Selection | May A, B and D kept; new paragraph added | The domain, area and pedagogical-significance prompts as listings (used for the original three branches); a paragraph stating how the civil and industrial domains and areas were chosen (section 9 decision 3); one sentence on who assigned difficulty levels | `7_appendix.tex` A, B, D; the authors |
| B | Source Texts and Reference Data | May C (Table 2) and E (Table 3) merged and extended | Table of textbooks by domain for all five branches (civil and industrial rows from the colleague's source set, section 9 decision 2, one primary text per domain); table of reference data sources with the provenance classes: NIST WebBook and fluid properties, CODATA 2022, NIST-JANAF, NASA TR R-132, MIL-HDBK-5J with NIST SP 811, USDA Wood Handbook, refractiveindex.info, PubChem, NAVFAC DM-7.01 and 7.02, FHWA HDS-4 and HEC-22, NRCS TR-55, USGS WSP 2339, the AISC shapes database, MIL-STD-105E, the NIST/SEMATECH e-Handbook | `docs/references/README.md`, `MANIFEST.json`; `civil_industrial_sources/references/public/MANIFEST.md` (local); `data/templates/branches/*/constants.py` |
| C | Template Construction and Examples | May F extended to five branches; May G cut to one or two listings | Bullets per branch on how parameters are grounded and constrained (add civil and industrial); one or two template listings from different branches, shorter than May's three; a small table of structural variation (templates that change the governing equation, step count or computed quantities with a sampled parameter; single-path templates) | `7_appendix.tex` F, G; `audit/template_audit_report.md`; `DIVERSITY.md` |
| D | Dataset Statistics and the Evaluation Set | May L (Table 9) replaced | Table: templates by branch and level, by domain with area counts; instances by level; answer kinds with counts. Paragraph: how the 2,250 instances were drawn (private seed, round-robin over reasoning paths and answer forms, exclusions), what is published now and at publication. Table or paragraph: question skeletons, reasoning paths, near-duplicate pairs | `manifest.jsonl`, `FREEZE.json`, `DIVERSITY.md`, D-114 to D-116 |
| E | Template Certification | May H and I kept; May J replaced; May K (Tables 5 to 8) replaced | E.1 Automated integrity checks: what each check tests, 500 seeds, the outcome, the excused-line register, the display-tie census (prose and one table). E.2 LLM screen: the prompt listing (May H, verbatim), the rubric table (May I), the judges and settings, the two passes with outcomes and agreement, the verified claims (table). E.3 Expert certification protocol: own-branch experts, hand check, the 20 planted defects (table: template, defect kind, description), opaque codes, timestamps, the interface (new screenshot or none). E.4 Outcome: plant detection, hand checks, agreement per branch (Fleiss, AC1, AC2), screen against experts (false-positive rate, MAD), round-by-round counts, the Markdown-rendering finding (tables) | `template_annotation_23092026/README.md`, `layer0/gate_report.md`, `tie_census.md`, `screen/pass1/`, `pass2/stats.md`, `pass1_fixes.md`, `layer2/RESULTS*.md`, `CERTIFICATION.md`, `markdown_scan.md` |
| F | Scoring Details | New; absorbs May N (prompt, adapted) and replaces May P.2 | F.1 Final-answer check: answer kinds (table with counts), how targets are derived, the match rule, partial and unreadable verdicts, the scoring statements as lists (nine symbolic-answer templates scored by stated numbers; four templates verified on the last value; 84 multipart instances checked on fewer quantities than the question's labelled parts; the 17 exact-digit templates). F.2 Milestones: derivation per instance, matching, the judge prompt (listing), verdict semantics, the judge's validation (table: 67 / 20 / 1 and 0 / 66 / 22), judged fraction per model (table). F.3 Arithmetic check: the rule, worked examples, what it skips. F.4 Judged step check: the prompt (listing; May N's rubric adapted to a whole response), precision and recall. F.5 Thresholds: the tables of `THRESHOLD_APPENDIX.md` | `full_run_28092026/answer.py`, `judge.py`, `router.py`, `THRESHOLD_APPENDIX.md`, D-138, D-145, D-169 |
| G | Validation Against Expert Judgment | May M (Tables 10 to 12) replaced | G.1 Study design: problems, models, experts, the label scheme (table replacing Tables 10 and 11), verification and adjudication rounds. G.2 Expert agreement (table, overall and per branch). G.3 Evaluator comparison (table replacing Table 12: the seven designs; answer agreement, milestone precision, recall and F1, response-level AUROC with template intervals, cost). G.4 Planted defects: design (table) and detection per checker and judge (table), the routing finding. G.5 Judge independence: judge selection and the exposure table, the judge swap (table), the panel ablation and per-judge bias (tables). G.6 Validation on the evaluated models: the arithmetic check's reading rounds (0.505, 0.752, 0.905) and the experts' readings of verdicts and judge decisions (tables). G.7 Corrections to the final-answer check after the evaluation: what, counts, the known wrong credit (list) | pilot `PILOT_SUMMARY.md`, `RESULTS_X1.md`, `RESULTS_E1.md`, `RESULTS_E5.md`, `RESULTS_LOJO.md`, `RESULTS_ATTRIBUTION.md`, `JUDGE_SELECTION.md`; `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, `JUDGE_SWAP.md`, `FLAG_REVIEW*.md`, `EXPERT_REQUEST.md`, `ANSWER_FORM_FIX.md` |
| H | Models and Inference Settings | May O (Table 13) and P.1 replaced and kept respectively | Table from the decoding table: model, developer, open or closed, served identifier, ceiling, endpoints, median and 90th-percentile completion and reasoning tokens, share of responses with reasoning tokens, cost; the selection rule and the one model set aside; the prompt (May P.1, unchanged, as a listing); output reading in one sentence pointing to F.1; the statistical protocol in full (unit, resampling, correction, detectable differences, the confirmatory questions and when they were fixed, the labelled exploratory additions) | `DECODING_TABLE.md`, `models.json`, `ANALYSIS_PLAN.md`, D-110, D-117, D-132 |
| I | Full Results | May Q (Figures 7 and 8) replaced by tables with intervals | I.1 Table 1 extended (FAC, strict FAC, unreadable share, SD within and between, MC with and without the judge, judged fraction, the two flag rates). I.2 Pairwise comparisons (55 rows: difference, interval, corrected p, detectable). I.3 FAC and MC by branch, level and domain with intervals, and the branch-pair tests. I.4 The level gap under its variations. I.5 By answer kind. I.6 Sensitivity table. I.7 Within-template consistency and instance variance. I.8 Failure against depth. I.9 Decoding repeats. I.10 Tokens against score. I.11 Coverage details (wrong-answer coverage against the floor, verdict against coverage, coverage against verbosity). I.12 Diagnostics details (flag rates with intervals, first-flag position, what points at a wrong answer, precision by error type) | `results/RESULTS.md`, `results.json`, `per_template.csv`, pilot `RESULTS_ATTRIBUTION.md` |
| J | Robustness to Rewording | New | Design (subset, writer, the rewriting rules), the scripted checks, the expert check and its rejection reasons with counts, survival by branch and the list of lost templates, the results table (FAC, MC with and without the judge; intervals; corrected p), the bounds table (90% intervals against ±5), tau against its noise floor, the comparison with decoding repeats | `PARAPHRASE.md`, `PARAPHRASE_REVIEW.md`, `PARAPHRASE_AUDIT.md`, `RESULTS.md` Q5, D-157 to D-165 |
| K | Additional Conditions | New | Reasoning effort, flagship anchors, open book, and the tool condition if run: one design paragraph each, the paired tables as `RESULTS.md` prints them, the anchors beside the eleven on the same instances, the per-condition decoding tables, the open-book block construction and, if reported as run, its defects | `RESULTS.md` "C1 and C4", "C3"; `DECODING_TABLE_*.md`; `OPENBOOK_SURVEY*.md`, `OPENBOOK_DIFF.md`; D-180 to D-184 |
| L | Error Analysis | May R kept with two added categories; May S, T, U replaced; May Table 15 redrawn or dropped | L.1 Taxonomy (Table 14's definitions and the decision hierarchy, with "No error" and "Incomplete" added). L.2 Sample and protocol (160 wrong-answer responses, four models, levels in proportion, three readers). L.3 Agreement (Fleiss' kappa overall and per model). L.4 Results by model and level (table). L.5 Representative responses, if new examples are drawn with the readers' consent. L.6 Templates whose wording does not pin the answer: the two chemical templates with the experts' answers, the 17 exact-digit templates, the near-miss distance table, the six templates whose answers can be read off the question | `7_appendix.tex` R; `EXPERT_REQUEST.md` B2, B4; `RESIDUAL_INCORRECT.md`; `SHORTCUT_AUDIT.md` |
| M | Extended Related Work | New | Table A (the fifteen closest benchmarks) and the further-related-work paragraphs | `RELATED_WORK_v3.md` appendix |
| N | Reproducibility and Release | New | What is released (templates, code, the evaluation set with the seed at publication, responses and scores, the paraphrases' hashes), what is withheld and why (the experts' labels, names), the anonymised archive's contents for review, and the script behind each table | this plan; `full_run_28092026/README.md`; section 9 decisions 1 and 5 |

Dropped from May without replacement: the textual-quality metrics (BERTScore, ROUGE) and their columns; the
ablation of evaluator configurations by rank correlation (Table 12); the Tribunal examples (Table 5), kappa and alpha
tables (Tables 6 and 7) and AI-human MAD table (Table 8) of the old certification; the annotation-interface
screenshot (Figure 6) unless replaced; the old sampling-strategy equations (S.1); the 11-model error tables (Tables
16 to 18); the word "Tribunal".

## 7. Figures and tables

| Item | Content | Source | State |
|---|---|---|---|
| Figure 1 | A template, one instance, its gold trace, and the milestones and calculations the evaluator checks | a certified template; `milestones.py` | revise `figs/engtrace-example-template.pdf` (overlay), or a civil or industrial template to show a new branch |
| Figure 2 | Five-branch taxonomy | `figures_oct_12/engtrace-overview-5branch.pdf` | exists, untracked; copy to `figs/` |
| Figure 3 | FAC and MC per model with intervals, tiers marked | `results/results.json` Q1, Q3 | script to write |
| Figure 4 | FAC on Easy, Intermediate and Advanced per model with intervals, the corrected verdict marked | `results.json` Q2 and the level table | script to write; may move to Appendix I if space is short |
| Figure 5 | Coverage on wrong answers against the chance floor, per model | `results.json` Q3 | script to write; may move to Appendix I |
| Figure 6 | Change in FAC under rewording, per model, 90% intervals and the ±5 band | `results.json` Q5 | script to write; may move to Appendix J |
| Table 1 | FAC, MC, compact letter display | `results.json` Q1, Q3 | script to write |
| Table 2 | Validation against the experts: each component's agreement, precision, recall; planted-defect detection | `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, RESULTS_X1 | collate by script |

Figures 1 to 3 and Tables 1 and 2 are the must-haves in the main text. Every figure and table is produced by a
committed script; figures are drawn to read in greyscale; the ACL column is 3.03 inches wide.

## 8. Writing rules

1. **No internal vocabulary in the paper. Strictly forbidden:** evaluator labels (E0 to E5, E5-strict), question
   and item labels (Q1 to Q5, A1 to A11, B1 to B4, C1 to C5), decision numbers, commit hashes, file, script and
   variable names, and the working vocabulary of the record. The table below gives the replacement for each term;
   the proofreading pass greps for the left column.

   | Internal term | Paper wording |
   |---|---|
   | E3, E4, E5, E5-strict, "coverage (stage)" | Milestone Coverage; "without the judge" and "with the judge" where both are shown |
   | digit rule, E4 flags | the arithmetic check; arithmetic flags |
   | router, step router, rule C | the judged step check; judged step flags |
   | residual judge, judge calls | the judge (first use: "a judge on the milestones the matcher does not find") |
   | Tribunal, AI Tribunal, panel | a panel of LLM judges (only when describing the compared design) |
   | Layer 0 / gate / T1 closure, T3, T4, T8 | automated integrity checks (each described in words) |
   | Layer 1 / screen pass 1, pass 2 | the LLM screen; its two passes |
   | Layer 2 / plants / kits / app / workbook | expert certification; planted defects; the annotation interface |
   | pilot, slice, frozen slice | the expert study; the validation set |
   | pool, freeze, manifest, seed commitment | the evaluation set; the published hashes; a commitment to the seed |
   | roster, roster rule | the evaluated models; the selection rule |
   | run, main run, arm, variant, subsample | the evaluation; condition; the 450-instance subset |
   | trace (where it means a model output) | response (or "trace" only for the reasoning text, defined once) |
   | unusable, empty, capped | no readable final answer; empty at the output ceiling |
   | fully solved, FSR, answer score | strict Final Answer Accuracy; Final Answer Accuracy |
   | REACHED / NOT_NEEDED / MISSING | reached, not needed, missing (in words) |
   | store, stage rows, score store, CONFIG | not mentioned |
   | owner, supervisor, agent, session | not mentioned |
   | pre-registered | "specified before the models were run" (one optional sentence) |
   | D-xxx, `RESULTS.md`, `analyze.py`, any path | not mentioned; Appendix N names the scripts behind the tables once |

2. **One voice, present tense, no history.** The paper describes what EngTrace is; "we" build, certify, score,
   validate. No "previously", "originally", "in response to".
3. **Every number from a committed script**, regenerated from `results.json`, the manifest or the validation files;
   never copied from a Markdown table by hand.
4. **The unit of analysis is named once and kept**: the template; every interval resamples templates; every family of
   tests is corrected; every null carries its detectable difference.
5. **Tiers, not an order**, wherever intervals overlap; "not separable at this size", never "equal".
6. **Flags are flags; coverage is coverage.** The two naming rules of section 4.
7. **Say the setting with the number**: provider defaults, the four models without reasoning tokens, the ceiling,
   the unreadable-answer rule, wherever a closed-against-open or a level sentence appears.
8. **Name the limits where the claim is made**: the judged fraction beside coverage, the precision beside a flag
   rate, the two chemical templates beside the level gap, the 115 templates beside the paraphrase bound.
9. **Honesty items that must appear**: the Markdown rendering during certification; the corrections to the
   final-answer check after the evaluation, with their counts and the one known wrong credit; the two under-specified
   chemical templates; the 17 exact-digit templates; the open-book condition's defects if reported as run; the four
   no-reasoning models; the validation population against the evaluated models.
10. **Length discipline**: cut evidence to the appendix before cutting a limit from the main text.

## 9. Decisions the authors must make before the text is final

1. **Anonymity.** Link an anonymised archive (ARR zip or anonymous.4open.science), never the repository, which
   carries author names, the decision log and the study record. FinChain in the third person; the preprint not cited
   as the authors' own; the project and code links removed from the title block of the review version.
2. **The authorship sentence for the civil and industrial templates.** They were written by a colleague of the
   authors and certified like the rest; the paper needs the exact description of that colleague (role and
   expertise) and a decision on whether any writing assistance needs disclosing under the ARR checklist. (The
   September plan's note that these templates were "AI-drafted" was a misreading, corrected on 2026-10-03 in
   `docs/NEXT_CYCLE_REVIEW.md`.) The sources the colleague used arrived on 2026-10-03 as
   `civil_industrial_sources/references/public/` (local, gitignored; its `MANIFEST.md` is the record): the
   public-domain data sources, which are the same files as `docs/references/civil` and `industrial` (NAVFAC DM-7.01
   and 7.02, FHWA HDS-4 and HEC-22, NRCS TR-55, USGS WSP 2339, the AISC shapes database; MIL-STD-105E and the
   NIST/SEMATECH e-Handbook), and the textbooks: for civil, Das and Sobhan, *Principles of Geotechnical
   Engineering*; Holtz, Kovacs and Sheahan, *An Introduction to Geotechnical Engineering*; Knappett and Craig,
   *Craig's Soil Mechanics*; Das and Sivakugan, *Principles of Foundation Engineering*; Hibbeler, *Structural
   Analysis*; Kassimali, *Structural Analysis*; Leet, Uang and Lanning, *Fundamentals of Structural Analysis*; Chow,
   Maidment and Mays, *Applied Hydrology*; Chin, *Water-Resources Engineering*; Sturm, *Open Channel Hydraulics*; for
   industrial, Hillier and Lieberman, *Introduction to Operations Research*; Ross, *Introduction to Probability
   Models*; Taha, *Operations Research: An Introduction*; Nahmias and Olsen, *Production and Operations Analysis*;
   Silver, Pyke and Thomas, *Inventory and Production Management in Supply Chains*; Montgomery, *Introduction to
   Statistical Quality Control*; Grant and Leavenworth, *Statistical Quality Control*; Mitra, *Fundamentals of
   Quality Control and Improvement*. Appendix B lists one primary text per domain, as the May table does; the
   colleague confirms which of these is primary for each of the six domains.
3. **The taxonomy and difficulty procedure for the new branches.** The May paper's section 3.1 says the domain and
   area lists were drawn from ABET and the professional societies and then cross-checked with an "expert-persona
   LLM": four models, prompted as a senior curriculum designer, ranked the cornerstone domains of each branch,
   listed the fundamental areas of each domain and scored each area's pedagogical significance, and the lists were
   settled by majority vote (May Appendices A, B and D). That check was run for the original three branches only.
   For civil and industrial, the domains and areas were chosen by the colleague who wrote the templates, from the
   ABET, ASCE and IISE curricula and the textbooks of decision 2. Two ways to write 3.1: (a) recommended, say
   exactly that, keep the LLM cross-check as a statement about the original three branches, and let the own-branch
   expert certification stand as the validation of all 150 templates; or (b) run the same three prompts on the two
   branches now with a current panel (about 120 short calls, one to three dollars, needs approval) and report how
   many of the chosen domains and areas the panel's lists contain, so one sentence covers all five branches. Who
   assigned the difficulty levels of the 60 new templates (the 2026-09-05 inventory holds them) is stated either way.
4. **The retrieval objection.** The corrected open-book condition and the tool condition appeared to be running in
   another session as this was written; once their rows are in `RESULTS.md`, the sentence in 5.2 is written from them
   in the form `RESULTS_PAPER_NOTES.md` gives. If either does not finish, the open-book condition is reported as run
   with its defects stated, or left out.
5. **The experts' labels.** Release an anonymised ground-truth file (codes and majority labels, no per-rater rows)
   with the responses and scores, or state that the labels are withheld and why.
6. **The two chemical templates.** The experts' answers support a stated limitation with the level-gap sensitivity
   (3 of 3 say the virial wording does not decide the form; 2 of 3 say standard heat-capacity sources differ by more
   than the tolerance). The evaluation set is fixed; no re-score.
7. **Strict or partial-credit accuracy as the headline, and the coverage metric's name** (section 4).
8. **LPDS.** Whether the Limitations carry the neutral sentence on randomly sampled instances (`CHANGES.md` section 4).
9. **Housekeeping that blocks claims**: `README.md` still says 90 templates and three branches; `figures_oct_12/` and
   `current_overleaf_project/` are untracked; Appendix H is generated from the decoding table; the pinned
   scoring-library versions are stated.
10. **Title**, as in section 2.

## 10. Order of work to 12 October

| Days | Work | Depends on |
|---|---|---|
| 4 to 5 Oct | Sections 3 and 4, the Related Work drop-in, Limitations and Ethics; the appendix skeleton with the tables that need no new data (certification, scoring details, thresholds, models and settings, dataset statistics) | nothing pending |
| 5 to 6 Oct | The figure and table scripts from `results.json` (Table 1, Figures 3 to 6, Appendix I); Figure 1 revised | nothing pending |
| 6 to 7 Oct | Sections 5.1 to 5.3; the conditions sentences in the form the decisions allow; Appendices J, K, L | decisions 4 and 6 |
| 7 to 8 Oct | Introduction, abstract, contributions, Conclusion; the authorship and taxonomy paragraphs | decisions 2 and 3 |
| 8 to 9 Oct | Appendices filled; cross-references; the bibliography regenerated from `related_work_sources.bib` plus the method citations; venue confirmations | — |
| 9 to 10 Oct | Two proofreading passes (numbers against `results.json`; the "do not" list, the naming rules and the forbidden-term grep); page budget; the anonymised archive; the responsible-NLP checklist | — |
| 11 Oct | Buffer; final compile in review mode; submit before 12 Oct AoE | — |
