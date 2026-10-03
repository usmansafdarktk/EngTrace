# EngTrace, October 2026 submission: the paper as it stands, and how to write it

Written 2026-10-03 for the authors, nine days before the ARR October deadline (12 October 2026, anywhere on
Earth). It reads the May 2026 submission (`docs/_ARR_May__EngTrace.pdf`), both rebuttals and both meta-reviews
(`docs/meta-reviews.md`), the current Overleaf sources (`current_overleaf_project/`), and the September and October
record: the template rework (`docs/re-implementation-sep/`), the certification (`template_annotation_23092026/`),
the evaluator pilot (`evaluator_pilot_17092026/`), the full run (`full_run_28092026/`, `results/RESULTS.md` at
commit `026210a`), the revised Related Work (`docs/related_work_oct2026/`), the five-branch overview figure
(`figures_oct_12/`) and the experts' readings returned today (`full_run_28092026/EXPERT_REQUEST.md`, uncommitted
as this is written). Every number quoted below names the file that prints it; nothing here is a new measurement.

**The one rule that shapes this document.** The paper is written for readers who have never seen EngTrace. It
describes the benchmark, the evaluator and the results as they are today, with their evidence and their limits.
Nothing in the paper says "previously", "in the earlier version" or "as promised in the rebuttal". Every sentence
of that kind belongs in the revision letter, for which section 11 keeps notes. The companion files
`full_run_28092026/RESULTS_PAPER_NOTES.md` and `PARAPHRASE_PAPER_NOTES.md` say, number by number, what each result
supports; this document says what the paper is and where everything goes.

---

## 1. The paper in one paragraph

EngTrace is a benchmark of 150 symbolic templates across five engineering branches (chemical, civil, electrical,
industrial, mechanical), 15 domains and 42 areas. Each template is executable code: it samples physically
grounded parameters from referenced constants, computes the answer, and writes the gold reasoning trace from the
same computation, so every intermediate quantity exists with its unit before any prose does. Templates reach the
benchmark through a deterministic integrity gate, an LLM screen whose judges share no family with any evaluated
model, and certification by three experts of the template's own branch, who caught 60 of 60 planted defects and
rejected 53 real templates across rounds until all 150 were approved unanimously. Evaluated instances (2,250: 15 per
template) are drawn from a private seed, so no evaluated item was public when the models ran. A trace is scored
deterministically wherever a check exists (a three-way answer check over six answer kinds, milestone coverage of
the gold derivation's intermediate quantities, an arithmetic check of every displayed claim) and a judge from a
family outside the roster decides only what those checks leave open. The evaluator was chosen by a meta-evaluation
of seven designs against 15 experts' step-level labels on 300 traces and 120 planted defects, which also fixed the
limit of verification: no deterministic check detects a misstated rule behind a correct answer, and LLM judges catch
about a third. Eleven open and closed models of September 2026 score 0.81 to 0.98 on the final answer; the top five
are inseparable at this size; every model scores lower on Advanced templates; milestone coverage orders the models
differently from the answer; on wrong answers models still reach most of the derivation; and expert-checked
paraphrases change no model's score by more than five points for ten of the eleven.

## 2. Title, abstract and contributions

**Title.** Keep *EngTrace: A Symbolic Benchmark for Verifiable Process Supervision of Engineering Reasoning*. The
mechanism now matches the title: intermediate values are verified against quantities the template computed, not
voted on. If a shorter form is wanted, *EngTrace: Verifiable Process Evaluation of Engineering Reasoning from
Symbolic Templates* says the same.

**Abstract, draft (about 250 words; trim to taste).**

> Engineering problems are solved by derivations whose intermediate quantities can be checked, yet most benchmarks
> score only the final answer, and those that assess steps rely on LLM judges or human graders. We present
> EngTrace, a benchmark of 150 symbolic templates spanning five engineering branches, 15 domains and 42 areas. Each
> template is executable code that samples physically grounded parameters, computes the answer and emits a gold
> reasoning trace from the same computation, so every intermediate quantity is known with its unit. Templates pass
> a deterministic integrity gate and an LLM screen, and are certified by three experts of their own branch, who
> caught 60 of 60 planted defects. Evaluated instances are drawn from a private seed, so no evaluated item was
> public when the models ran. EngTrace scores a trace deterministically wherever a check exists: a three-way answer
> check over six answer kinds, milestone coverage of the gold derivation's intermediate quantities, and an
> arithmetic check of every displayed claim; a judge from a model family outside the evaluated roster decides only
> what these checks leave open. The evaluator is validated against step-level labels from 15 experts on 300 traces
> (answer agreement 0.98, milestone F1 0.96) and against 120 planted defects, which also show the limit of
> verification: no deterministic check detects a misstated rule behind a correct answer, and judges catch about a
> third. On 2,250 instances, eleven models score 0.81 to 0.98 on the final answer, with the top five inseparable
> at this size; every model scores lower on Advanced templates; milestone coverage orders the models differently;
> and on wrong answers models still reach 46 to 85% of the derivation. Expert-checked paraphrases change no
> model's score by more than five points for ten of eleven models.

**Contributions, draft.** Three, each owned by an active verb and each carrying its evidence:

1. **We build and certify EngTrace**, 150 symbolic templates across five branches whose gold traces are computed
   by the template code, with a private-seed instance pool and a three-layer certification whose evidence is
   reported: an integrity gate at 500 seeds, a screen by judges from outside the roster, and own-branch expert
   certification with hand checks, planted defects (60 of 60 caught) and 53 rejections fixed before the
   unanimous approval of all 150.
2. **We design and validate a deterministic-first process evaluator**: an answer check over six answer kinds,
   milestone coverage derived from the template's own computed quantities, an arithmetic check of displayed claims,
   and a judge from outside the roster only on the residue. We choose it by a meta-evaluation of seven evaluator
   designs, including a matching-plus-judge-panel design and three process reward models, against 15 experts'
   step-level labels on 300 traces and 120 planted defects, and we report what verification can and cannot see.
3. **We evaluate eleven models under a pre-registered analysis** with the template as the unit, template-level
   intervals, Holm-corrected tests and the smallest difference the design detects: tiers rather than an order at
   the top, the level gap per model, what process scores add on wrong answers, failure against derivation depth, a
   paraphrase bound, decoding repeats, and reasoning-effort, flagship and open-book conditions.

## 3. The claims the data support, and the claims not to make

Each claim below is a sentence the paper may make, with its source. The "do not" list is binding: it collects
every overstatement the record has already caught (`RESULTS_PAPER_NOTES.md`, "What not to state", in both of its
sections; `PARAPHRASE_PAPER_NOTES.md`; D-111; D-146; D-171; D-180 to D-183).

### 3.1 Claims

| # | Claim | Evidence (file) |
|---|---|---|
| C1 | 150 templates, 30 per branch, 15 domains, 42 areas, 58 / 58 / 34 templates by level; 2,250 instances (870 / 870 / 510; 450 per branch), 2,250 distinct questions, drawn from a private seed with per-item hashes and a seed commitment committed | `full_run_28092026/FREEZE.json`, `manifest.jsonl`, D-114 to D-116 |
| C2 | Every template passes a deterministic gate at 500 seeds (printed-arithmetic closure, determinism across processes, output contract, emission); a screen by three judges from no roster family passed 147 of 150 on the shipped corpus with 3 controversial and 0 critical; three own-branch experts certified all 150, after 53 rejections and fixes, catching 60 of 60 planted defects, with 459 of 504 hand checks within 1% | `template_annotation_23092026/layer0/gate_report.md`, `screen/pass2/stats.md`, `layer2/RESULTS.md`, `CERTIFICATION.md`, D-109 |
| C3 | The evaluator's deterministic parts agree with 15 experts' labels on 300 traces: answer check 0.982 (correct-or-not) and 0.927 three-way; milestone F1 0.923 deterministic and 0.958 with the residual judge; the judge gave REACHED to 0 of 88 fabricated values; the arithmetic rule has precision 0.817 and recall 0.427 on slips inside correct-answer traces | `full_run_28092026/SCORER_VALIDATION.md`, `E5_VALIDATION.md`, pilot `RESULTS_X1.md`, `THRESHOLD_APPENDIX.md` |
| C4 | No deterministic check detects a misstated rule behind a correct answer (0 of 60 planted conceptual defects for every checker); LLM judges catch about a third (MiMo 16 of 52, GPT-5 20 of 60, Grok 22 of 60) with no false alarm on untouched steps; the final answer nearly decides the experts' trace verdict (AUROC 0.974) | pilot `RESULTS_X1.md` Findings 1, 5, 7, 8; D-181 |
| C5 | Answer score 0.814 (gpt-oss-20b) to 0.976 (DeepSeek V4.1 Flash); the top five lie within 0.009 and no pair among them separates; 32 of 55 pairs differ after Holm; 23 pairs lie below the detectable difference, all of them non-significant | `results/RESULTS.md` Q1 |
| C6 | Every model scores lower on Advanced than on Easy templates (+0.055 to +0.199); the gap holds after Holm for the four lowest-scoring models (Welch) and for none once two Advanced chemical templates whose wording does not pin the answer are removed | Q2 and its three variations; D-171; `EXPERT_REQUEST.md` B4 |
| C7 | Milestone coverage runs 0.816 to 0.923 and orders the models differently from the answer score (Kendall's tau 0.709, 95% CI 0.514 to 0.855); 28 of 55 pairs differ (sign-flip), 29 under Wilcoxon; coverage does not rise with the number of steps written (rho −0.11 to −0.26) | Q3 "Coverage compared across models", "Does coverage track verbosity?" |
| C8 | On wrong answers models still reach 46% to 85% of the milestones (readable traces) against a chance floor of 10% to 22%; 13% to 26% of wrong answers are a complete derivation to a wrong value (the five models with more than 100 wrong answers) | Q3 "Milestones on the wrong-answer traces", "Verdict against coverage" |
| C9 | The arithmetic rule flags 0.5% to 22.5% of fully solved traces, reading 2.8 to 9.8 claims per trace, at precision 0.905 on this roster (171 of 189 flags read by a domain expert); the router's judge flags 1.1% to 22.0% | Q3 "The digit rule", "The step router"; `FLAG_REVIEW_3.md` |
| C10 | The wrong-answer rate on items whose gold derivation has six or more milestones is 0.037 to 0.291, against 0.000 to 0.045 on one-milestone items, for every model | Q3 "Wrong-answer rate against the item's milestone count" |
| C11 | The strongest models solve all 15 instances of 85% to 91% of the single-path templates; the weakest 36% to 48% | Q4 |
| C12 | On 277 expert-kept paraphrase pairs over 115 templates, the paired change in answer score is −0.031 to +0.020, none significant after Holm; for ten of eleven models the 90% interval lies within ±5 points; the tiers hold and the order within the top tier is not resolved (tau 0.673 against a noise-floor median of 0.782 at this size) | Q5; `PARAPHRASE_PAPER_NOTES.md` |
| C13 | Decoding repeats on four models (300 items, three repeats): SD 0.004 to 0.013, the same verdict every time on 86% to 94% of items | "Decoding repeats" |
| C14 | The ordering is unchanged at half and double the answer tolerance (tau 0.927), without the six surface-shortcut templates (headline within 0.002) and without the nine symbolic templates (+0.002 to +0.011) | Sensitivity; `THRESHOLD_APPENDIX.md`; `SHORTCUT_AUDIT.md` |
| C15 | Reasoning at medium effort moves GPT-5.4 mini from 0.858 to 0.951 on the 450-item subsample (+0.093, 95% CI 0.056 to 0.133) and Gemini 3.1 Flash-Lite by less than the design detects (+0.016, −0.010 to 0.042) | "C1 and C4"; D-180 |
| C16 | Two flagships that pass the roster rule do not exceed the roster's top tier on the same 450 items: DeepSeek V4 Pro 0.968 and GPT-5.4 with reasoning on 0.980 sit inside the top five's intervals (0.964 to 0.976); GPT-5.4 at its provider's default (no reasoning tokens) scores 0.941, below them | "C3"; D-182 |
| C17 | Supplying each template's governing equations with the question (open book, as run) changed the answer score by +0.005, −0.007 and +0.041 for the three models tried, none holding after Holm, while stated coverage rose by 0.03 to 0.06; the equation blocks had defects (D-183) and a corrected arm awaits a decision | "C1 and C4"; D-183 |
| C18 | Judge independence: the full run's judge shares no family with the roster; on a 220-trace sample a second judge from a third family moves no model's coverage figure beyond its interval (−0.017 to +0.020); on a matching-plus-panel design replayed offline, dropping a family's own judge changes its score by at most 0.006, the same as a placebo drop, and every judge is lenient on every family alike | `JUDGE_SWAP.md`; pilot `RESULTS_LOJO.md`; D-174, D-181 |
| C19 | On this roster, two experts read 150 of the top five models' answer verdicts (half incorrect or partial by design): the check's correct verdicts were confirmed in 142 of 150 readings, its incorrect verdicts in 83 of 98 (12 expert-correct), and its partial verdicts were called fully correct in 34 of 52, so the half-point is conservative; the judge's REACHED verdicts were confirmed in 85 of 100 readings and its MISSING verdicts in 79 of 100 | `EXPERT_REQUEST.md`, "What came back" (B1, B3) |
| C20 | Three readers labelled 160 wrong-answer traces of four contrastive models with the six-category taxonomy plus "No error" (Fleiss' kappa 0.930): for Claude Sonnet 5, 26 of 40 item majorities are "No error" (the answer is correct or the question admits it); for gpt-oss-20b, Gemma 4 and GPT-5.4 mini, Calculation Error leads (17, 28 and 28 of 40), then Formula/Principle; Easy wrong answers are almost all calculation (69 of 78 readings) | `EXPERT_REQUEST.md` B2 |

### 3.2 Do not state

- Not "model X is the best": tiers, not an order, at the top; the top five are within 0.9 points and no pair separates.
- Not "trade-off", "decoupling" or "trace fidelity" for milestone coverage: it is coverage of stated intermediates; no test of a trade-off was planned and the intervals overlap.
- Not "x% of correct answers rest on flawed reasoning" from the arithmetic rule or the router: flag rates with known precision (0.905; about 0.70) and recall below one half.
- Not a ranking of models by arithmetic reliability: the rule reads 2.8 to 9.8 claims per trace depending on style.
- Not "the complexity cliff holds" without both counts (4 of 11 as scored; 0 of 11 without the two chemical templates) and the per-model intervals; not the planned permutation's 9 of 11 as the headline.
- Not "closed models underperform open ones" or any comparison at equal effort: four models ran without reasoning tokens.
- Not "chemical engineering is hardest" or any branch order: one branch pair in the roster holds (gpt-oss-20b, electrical above civil); domain means are descriptive.
- Not "2,250 unique problems": 2,250 instances and 2,250 distinct questions from 150 templates.
- Not "robust to paraphrasing" or "no contamination": a bound of ±5 points for ten models on wording memorisation; method memorisation is not tested.
- Not "the ranking is stable under rewording": tau 0.673 sits at the lower edge of the sampling-noise distribution; the tiers hold.
- Not "flagships do no better" and not "solved at the top" without the Advanced means beside it (0.87 to 0.95 for the anchors, 0.92 to 0.95 for the top five, wide intervals).
- Not "retrieval would not help": the open-book arm supplies equations only, had defective blocks as run, and the tool condition was not run.
- Not "no self-preference" as evidence that LLM judges are accurate: they are uniformly lenient; the ablation answers the overlap objection and nothing else.
- Not templates as a remedy for benchmark saturation (Akhtar et al. 2026 find no such effect); the private seed answers contamination, and the step-level checks answer the question a saturated answer score cannot.
- Not the pilot's agreement figures as validation of answer forms the pilot templates do not contain (pi-fractions, two-unit lines; D-169): say the experts did not arbitrate those.
- Not "no existing benchmark verifies the process" (PhysReason, FinChain, PRISM-Physics, ThermoQA contradict it); the claim is the combination, hedged as in `RELATED_WORK_v3.md`.
- Not any number beside a May number. The old pool, roster and check are not reproducible; the paper describes the evaluator it uses and the revision letter explains the change.

## 4. The metrics: names, definitions, and what each is reported with

The paper uses reader-facing names. The internal labels (E3, E4, E5, Q1 to Q5) stay in the repository and in one
appendix table that maps names to scripts. Every metric is reported with a 95% interval that resamples templates
(B = 10,000) unless marked descriptive, and every null with the smallest difference the design detects at 80% power.

| Paper name (abbreviation) | Definition | Reported with | Internal |
|---|---|---|---|
| **Answer Score (AS)** | Mean over items of 1 for a correct final answer, 0.5 for a partial one (some but not all verified targets match), 0 for incorrect or unusable. The check accepts a value within 0.002 relative of the gold or within one unit of the last displayed digit, after unit factors; targets are the quantities the template computed; six answer kinds (scalar, vector, array, symbolic, classification, multipart) | interval; per-template SD within and between; the fully-solved rate; unusable rate; pairwise Holm verdicts and detectable differences | `answer.py`; Q1 |
| **Fully-Solved Rate (FSR)** | Share of items on which every verified target matches; a partial scores 0. (The name "Final Answer Accuracy" may be kept for this column if the authors prefer; define it once) | interval | Q1 |
| **Unusable Rate** | Share of traces in which the check can read no final answer: empty at the output ceiling, cut off, or never stating one; scored 0 in AS and FSR | per model, beside AS; the six templates that hold most of them | Q1; `TRACE_REVIEW.md` |
| **Milestone Coverage (MC)** | Share of the gold derivation's milestones (the intermediate quantities the template computes and the gold trace states) that the trace states within 0.5% relative, unit-aware and order-free, or that the residual judge rules REACHED. Items with no milestones are left out. Measures progress through the derivation, not the absence of error | interval; the deterministic part alone (MC without the judge); the judged fraction and the share of judged milestones ruled REACHED; the unjudged rate; the chance floor (the trace scored against a sibling item's milestones) | E3 and E5-strict; Q3 |
| **Arithmetic Flag Rate (AFR)** | Share of fully solved traces with at least one displayed arithmetic claim whose recomputation from the printed operands is not a correct rounding at the printed precision (the digit rule) | claims checked per trace; the rule's precision on this roster (0.905) and its precision and recall on the expert labels (0.817 / 0.427); the rate at a 1% tolerance beside it | E4; Q3 |
| **Step Flag Rate (SFR)** | Share of fully solved traces in which the router's judge flags a step the digit rule did not (every unflagged step of a trace judged in one call) | its precision and recall on the expert labels (0.707 / 0.603 over all steps; 0.703 / 0.360 inside correct-answer traces); the share of flagged steps by category read through `RESULTS_ATTRIBUTION.md` | router; Q3 |
| **Level gap** | AS on the Easy templates minus AS on the Advanced templates, templates resampled within each tier; Welch's t-test with Holm across models | interval; the count under Welch and under the planned permutation; the detectable gap; the three variations (unusable rows left out, no symbolic templates, without the two chemical templates) | Q2 |
| **Within-template consistency** | Share of templates fully solved on all 15 instances, on some, and on none; single-path templates (one reasoning path across their instances, 58) reported apart from the other 92 | descriptive, intervals in `results.json` | Q4 |
| **Depth failure rate** | Wrong-answer rate by the item's milestone count (0, 1, 2, 3, 4 to 5, 6 or more) | descriptive | Q3 |
| **Paraphrase change** | Paired AS difference, paraphrase minus original, over the expert-kept pairs; sign-flip permutation over templates with Holm across models; the 90% interval read against a ±5-point margin (TOST) | interval; the detectable paired difference; the same for MC; Kendall's tau between arms against its noise floor at the arm's size | Q5 |
| **Repeat SD** | SD of AS across three decoding repeats on 300 items, four models | the share of items with the same verdict in every repeat | "Decoding repeats" |
| **Condition change** | Arm minus base on the same items, paired, with the paraphrase machinery (reasoning effort, open book, tool; the flagship anchors are outside the pairwise family) | interval; Holm within the arm; detectable; 90% interval against ±5 | "C1 and C4", "C3" |

Two naming rules. First, "coverage" never becomes "reasoning quality", "soundness" or "fidelity" anywhere in the
text, captions included. Second, flag rates are "flags", never "errors": the sentence "the rule flags x% of
traces, and 90% of its flags are real slips" is right; "x% of traces contain an arithmetic error" is not.

## 5. Sections, with page budget, content, sources and writing notes

ARR long paper: eight pages of content, unlimited references and appendices; the review version must stand on
its own without the appendix. The budget below sums to eight pages and is tight; the appendix carries everything
that is evidence rather than argument. The Overleaf project already uses the ACL style with `\usepackage[preprint]{acl}`;
switch to `review` for the submission.

### 5.1 Introduction (1.0 page)

**Write.** Open with the gap in the form the revised Related Work supports: engineering reasoning is a derivation
whose intermediate quantities can be checked, yet most engineering benchmarks score only final answers, and those
that assess steps do so with LLM judges or human graders, or within one domain. Then EngTrace in one paragraph
(section 1 above). Then a figure-referenced example (Figure 1: a template, one instance, its gold trace, and the
milestones the evaluator derives from it). Then the scope sentence: a depth-oriented evaluation of a principled
subset of engineering, verifiable because the gold trace is computed; items are examination-style by construction,
which is what verifiability needs and which leaves out ill-posed problems and real artefacts. Then the three
contributions of section 2.

**Sources.** `docs/related_work_oct2026/CHANGES.md` section 2 has the replacement sentences for every Introduction
sentence that cites a prior benchmark (GLUE and BIG-Bench Hard for saturation; MMLU, MATH and HumanEval; the
"outcome matching" claim about EngiBench and FEABench, which was wrong). Use them as written.

**Do not.** No "first", no "no existing benchmark verifies this process". No reviewer-facing language. Do not call
the benchmark "safety-critical"; the deployment context motivates it (one clause), and the Ethics section says
current models should not be deployed in safety-critical systems.

### 5.2 Related Work (0.75 page)

**Write.** Drop in the main text of `docs/related_work_oct2026/RELATED_WORK_v3.md` (`related_work_v3.tex` is the
LaTeX; `related_work_sources.bib` the bibliography, 100 entries, generated from arXiv, Crossref and the ACL
Anthology). Its five paragraphs: reasoning benchmarks and what they score; generated instances and contamination;
engineering benchmarks; process supervision and the validity of evaluators; positioning. At 704 words of prose it
is about 0.9 page; the file names the cut order to reach 0.75 (the coding benchmarks in the first sentence, all but
three single-branch benchmarks, the PRIME and ChemCoTBench-V2 sentence). Table A (fifteen closest benchmarks) and
the further-related-work paragraphs go to the appendix.

**Anonymity.** FinChain shares four authors with EngTrace: third person throughout, never "our earlier work". The
EngTrace preprint is public and is not cited as the authors' own. Three 2026 papers cite EngTrace by name (ERI,
Sci-rho, LPDS); the text cites them without saying so. Whether the Limitations acknowledge LPDS's finding on randomly
sampled instances is the authors' call (`CHANGES.md` section 4).

**Venues to confirm before submission** (`CHANGES.md` section 6): EMNLP 2026 for Huang et al. and Dlugosz et al.

### 5.3 The EngTrace benchmark (1.5 pages)

**3.1 Scope and taxonomy (0.3).** Five branches, 15 domains, 42 areas, 30 templates per branch, 10 templates per
domain except chemical (10 / 12 / 8). Figure 2 is the five-branch overview (`figures_oct_12/engtrace-overview-5branch.pdf`;
untracked, to be committed and copied into the Overleaf `figs/`). The procedure: curricular standards (ABET and the
societies: ASME, IEEE, AIChE for the original three; ASCE and IISE for civil and industrial), textbook grounding
(Appendix B lists the texts; `docs/references/` holds the on-disk sources every constant resolves to), and an
expert-persona model panel used to check the domain and area lists for the original three branches. **Say exactly
what was done for the new two branches** (section 9, decision 3); the May text's four-model panel was not run on
them.

**3.2 Templates and instances (0.5).** What a template is: a seeded Python function that samples parameters (often
inside a rejection loop that enforces physical validity: laminar regime, elastic range, upright floating), computes
every quantity, binds each displayed value through its display before it is consumed (round then recompute, with
decimal tie-breaking), asserts physical invariants, and emits the question and the gold trace with `**Step n:**`
markers and an `**Answer:**` line. Six answer kinds. Instances differ structurally as well as numerically: 72 of
150 templates change the governing equation, the step count or the computed quantities with a sampled parameter
(`docs/re-implementation-sep/audit/template_audit_report.md`), and the pool holds a second reasoning path for every
template that can produce one within 500 draws (`DIVERSITY.md`; 58 templates are single-path). Difficulty levels
(Easy, Intermediate, Advanced) by conceptual complexity, mathematical sophistication and procedural depth, with
the counts; who assigned them (section 9, decision 3). The pool: 15 instances per template drawn from a private
128-bit seed, chosen round-robin across each template's reasoning paths and answer forms from its first 100
acceptable draws, excluding repeated questions, pilot questions and display ties; per-item SHA-256 and a seed
commitment are committed, the seed and the text are released at publication. One sentence on constants: every
constant table carries a provenance class and resolves to a referenced source (NIST, CODATA, JANAF, MIL-HDBK-5J,
NAVFAC, FHWA, MIL-STD-105E, the NIST/SEMATECH handbook, and the textbooks).

**3.3 Authorship and certification (0.6).** Authorship, exactly: the 90 chemical, electrical and mechanical templates
were written by domain experts from the textbooks; the 60 civil and industrial templates were drafted with AI
assistance under a textbook-grounded authoring specification with blind AI review cycles and then put through the
same three layers as the rest, where the experts rejected 22 of the 150 in round 1 (section 9, decision 2, for the
record the paper must be able to point to). Then the three layers, each with its number: **Layer 0**, the
deterministic gate at 500 seeds (printed-arithmetic closure, determinism across processes, output contract,
emission), 150 of 150 pass, 54 templates edited to pass; **Layer 1**, the LLM screen with the published rubric and
three judges from no roster family (Grok 4.6, MiniMax M3, MiMo-V2.5-Pro), two passes capped, 147 pass / 3
controversial / 0 critical on the shipped corpus, AC1 0.93 on the flag, every flag verified before any edit (34 of
45 claims confirmed); **Layer 2**, three experts of the template's own branch per template, 34 items each (30
templates plus 4 planted defects: a wrong constant, a wrong unit conversion, a flipped sign, a printed step that
does not follow), a hand check before the solution is shown, opaque codes, timestamps; 60 of 60 plants rejected,
459 of 504 hand checks within 1% (44 of the 45 mismatches ended in a rejection), AC1 0.914 on approve/reject (Fleiss'
kappa 0.589 beside it, deflated by prevalence), AC2 0.947 / 0.980 / 0.954 on the three quality scores; 53 rejections
over the rounds (43 in round 1, 10 in round 2), each verified against its instance and fixed or, in two cases,
answered; rounds 2 to 4 re-certified the changed templates (17 of 22 approved by all three in round 2, 5 of 5 in
round 3, 2 of 2 in round 4); 150 of 150 certified, each by a unanimous latest round on code that still regenerates
the reviewed instances byte for byte. The screen's false-positive rate against the experts (10.2%) and the MAD
between panel and expert medians (0.59 / 0.19 / 0.79) show what an LLM screen does and does not see, and are the
reason the screen is a reader and the experts the gate. **State the limitation**: the review app rendered the text
as Markdown, which removed multiplication and dollar signs from what the experts saw in 65 templates
(`layer2/markdown_scan.md`).

**3.4 Dataset statistics (0.1).** One table (Appendix C holds the long form): templates per branch by level and by
domain, items per level, answer kinds (scalar 89 templates, multipart 32, symbolic 9, vector 8, array 7,
classification 5).

**Do not.** No "gold standard" in quotation marks, no "high-fidelity". Give the rejection counts as the evidence
that the review was real; a certification that rejected nothing would be the weaker claim.

### 5.4 The evaluation framework (1.5 pages)

**4.1 What a trace is scored on (0.5).** Figure 3, new: the stack as a pipeline. (a) The answer check: the answer
segment read from the final-answer heading; one target per quantity the template computed; correct / partial /
incorrect; tolerance 0.002 relative or one unit of the last displayed digit, after unit factors; LaTeX, fractions,
pi-fractions and superscripts read as numbers; six answer kinds. (b) Milestones: the intermediate quantities the
gold states, derived per instance from the template's own values, matched order-free within 0.5% and under unit
factors (the deterministic part); the milestones not found go to one judge with the question, the trace, and each
milestone's name and value (not the gold), who answers REACHED, NOT_NEEDED or MISSING; only REACHED is credited,
because on validation the judge used NOT_NEEDED to excuse a quarter of fabricated values. (c) The arithmetic
check: every displayed claim `expression = value` recomputed from the printed operands and compared at the printed
precision, with unit factors; the rule the experts' own guide uses (rounding is not an error, a wrong digit is).
(d) The step router: every step the digit rule does not flag is sent, a trace's steps in one call, to the same judge
on a four-way rubric (alternative correct / calculation error / conceptual error / other). The judge is MiMo-V2.5-Pro,
chosen because it shares no family with any evaluated model and detected planted defects as well as GPT-5 at a
third of a cent per call; everything else is code. State the judged fraction (11% to 25% of required milestones)
so the reader sees how much of the score a judge decided.

**4.2 Metrics (0.3).** The definitions of section 4, compressed; the unit of analysis (the template), intervals,
Holm, detectable differences, the two scoring rules (unusable scores 0; partial 0.5) with their reasons, and that
the analysis plan was fixed before the data (one sentence; Appendix H holds the plan).

**4.3 Validation against experts (0.5).** Present the pilot as a meta-evaluation, for a reader who has never seen
it: 60 problems (15 templates, five branches, three levels, four instances each), five models (GPT-5, Claude Opus
4.7, Gemini 3.1 Pro, DeepSeek R1, Llama 3.1 70B), 300 traces; 15 experts, three of the trace's own branch per
trace, labelling every step (correct / alternative correct / incorrect with type and reason / not a claim), every
milestone, the final answer and the trace's soundness; a reason on all 1,042 incorrect steps, a blind re-labelling
round and a blind adjudication of every 2-to-1 split (272 steps); Fleiss' kappa 0.781 between experts on steps and
0.828 within an expert, 0.880 on milestones, 0.966 on the answer. Seven evaluator designs were scored against those
labels: a step-matching evaluator with a three-judge panel (the design of several prior benchmarks and of an
earlier EngTrace evaluator; describe it generically in the anonymous paper), the same with a panel from outside the
roster, three process reward models, deterministic milestone coverage, coverage plus the arithmetic check, and
coverage plus a residual judge. Table 2: answer check 0.982 / 0.927 against the panel design's 0.760; milestone F1
0.923 and 0.958; the judge's validation (REACHED never given to a fabricated value, 67 of 88 true ones found);
digit rule 0.817 / 0.427 inside correct-answer traces, 0.872 precision over all steps; router 0.707 / 0.603; PRMs
rank steps well (AUROC 0.825 for the 72B model) but three of four flags inside correct-answer traces are false. Then
the finding (C4): the experts' own answer verdict predicts their soundness verdict at AUROC 0.974; 93 of 228
correct-answer traces carry an incorrect step, 175 of the 178 such steps arithmetic; milestone coverage cannot see
them (below chance); the digit rule finds about half of the slips it can parse; on 120 planted defects no
deterministic check catches a misstated rule (0 of 60) and judges catch about a third, with no false alarm. The
power statement: 300 traces over 15 templates detect AUROC differences of 0.12 to 0.19 and no smaller. Then the
roster-specific validation: the digit rule's precision 0.905 on 189 flags read by a domain expert after three
correction rounds, and the experts' readings of the top five models' verdicts (C19) and of the judge (C19). Then
one sentence on corrections: the answer check was corrected after the run on forms the pilot templates did not
contain (LaTeX numbers, stray digits, verdict words, inclusive last-digit windows, pi-fractions, two-unit lines),
each correction measured on the gold, the labelled traces and every trace before adoption, none moving a pilot
verdict against the experts, the last moving 238 of 24,750 verdicts, 234 to correct, every one read; one wrong
credit is known and stated.

**4.4 Judge independence (0.2).** C18 in one paragraph: the roster rule (no pilot generator, no judge); the judge
swap on a 220-trace sample; the leave-one-judge-out replay on the panel design with placebo and untouched controls
(no family effect at a detectable 0.003 to 0.025); the per-judge bias, positive and uniform across families, so a
panel's error is the task's, not self-preference.

**Do not.** No E0 to E5 in the main text. No "our evaluator beats the old one": the paper describes the evaluator
it uses and the meta-evaluation that chose it. Do not quote the threshold appendix's re-run pilot rows as the
published pilot figures (they differ in the third decimal); `SCORER_VALIDATION.md` reproduces the published ones.

### 5.5 Experimental setup (0.4 page)

**Write.** The roster: eleven models, eight open-weight (gpt-oss-20b, Gemma 4 26B-A4B, DeepSeek V4.1 Flash,
Qwen3-235B-A22B-2507, GLM-5.3-Flash, GLM-5.3, Muse Glimmer 30B, Kimi K3) and three closed (GPT-5.4 mini, Gemini
3.1 Flash-Lite, Claude Sonnet 5), under the roster rule; no flagship by a budget choice within the rule, with two
flagships run as anchors on a subsample (section 5.6). Decoding: each provider's default, a 32,768-token ceiling
(16,384 for Muse), one completion per item, the prompt byte-identical and hashed; **four models' endpoints
returned no reasoning tokens** (Gemma 4, Qwen3-235B-2507, GPT-5.4 mini, Gemini 3.1 Flash-Lite), so the roster is not
compared at equal reasoning effort, and Appendix I states the setting per model. Zero-shot, closed book, no tools,
with the open-book and tool conditions as labelled arms. Scoring: one trace per item, 24,750 traces; an unusable
trace scores 0 (246, 1.0%, 243 of them empty at the ceiling, 153 on six templates, four iterative by construction).
The analysis plan fixed before the data: the five confirmatory questions and their tests, the sensitivities, and
what is reported but not tested; exploratory additions labelled as such.

**Sources.** `DECODING_TABLE.md`, `TRACE_REVIEW.md`, `ANALYSIS_PLAN.md`, `models.json`, D-110, D-117, D-132.

### 5.6 Results (2.0 pages)

**6.1 Answer score and tiers (0.5).** Table 1 (AS with interval, FSR, unusable rate, SD within, MC with interval,
AFR with claims per trace, SFR; eleven rows) and Figure 4 (per-model intervals as dots, or the 55 Holm verdicts as
a matrix). The paragraph is the suggested wording of `RESULTS_PAPER_NOTES.md` ("For the results, first paragraph"):
five models between 0.967 and 0.976, no pair separating; GLM-5.3 at 0.948 owing most of its distance to 100
empty traces at the ceiling; the remaining five between 0.814 and 0.886; 32 of 55 pairs; 23 pairs below the
detectable difference. The footnote to Table 1 names the two closed models' reasoning-on scores on the subsample.
Instance variance in one sentence: most templates are solved on every instance or vary little; the variance
concentrates in a handful of templates per model, named in Appendix J.

**6.2 The level gap (0.3).** Figure 5 (Easy / Intermediate / Advanced per model with intervals, the Welch verdict
marked). Every model lower on Advanced by 5.5 to 19.9 points; significant after Holm for the four lowest-scoring
models, three of which ran without reasoning tokens; for the top five 6 to 7 points with intervals mostly
excluding zero but not surviving correction; the detectable gap 7 to 15 points; GLM-5.3's gap mostly an
output-ceiling effect; and the sensitivity: without the two Advanced chemical templates whose wording does not pin
the answer (the virial-work template admits two textbook readings, and the flame-temperature question gives no
heat-capacity data; three chemical experts confirmed both, Appendix L), the gap holds for none of the eleven.
Both counts in the same paragraph.

**6.3 What the process scores add (0.6).** Milestone coverage with its ordering (the suggested wording "For the
coverage comparison"); Figure 6 (coverage on wrong-answer traces against the chance floor, per model); the
sentence on complete derivations to a wrong value (13% to 26%); the attribution table read through its precisions
(a MISSING milestone or a router flag on most wrong answers of the weaker models, on fewer of the stronger ones,
whose wrong answers are mostly complete derivations with a wrong value); the arithmetic and step flag rates with
their precisions and the claims-per-trace column; and depth (Figure 7 or a sentence: 0.04 to 0.29 on items with six
or more milestones against 0.00 to 0.05 on one-milestone items). Within-template consistency in two sentences
(C11), with the point that on a single-path template a partly solved template is a failure on the numbers, not the
method.

**6.4 Robustness (0.3).** The paraphrase result as a bound, in the suggested wording of `PARAPHRASE_PAPER_NOTES.md`
(the 450 / 316 / 277 funnel, the ±5 bound for ten of eleven, Qwen3-235B-2507's −6.7 lower end, tau against its
noise floor, the experts' 12% rejection as evidence for the method), the decoding repeats, and the sensitivity
sentence (tolerance, shortcut templates, symbolic templates).

**6.5 Conditions (0.2).** The roster paragraph's reasoning-on sentence (C15); the flagship sentence (C16); the
open-book sentence (C17) in the form the record allows at writing time (section 9, decision 4). Appendix M holds
the arm tables.

**6.6 Error analysis (0.1 in this section; the rest in Discussion).** C20 in three sentences: three readers, 160
wrong-answer traces of four contrastive models, the six categories plus "No error", kappa 0.930; what dominates per
model and per level; and that for the top model most "wrong" answers are derivations the experts accept, which
is the saturation reading of 6.1 made concrete.

**Do not.** Every "do not" of section 3.2. Also: not the anchors in Table 1 or the pairwise tests; not the level
means of the subsample as tested differences; not the arms' coverage rise as better reasoning (handed equations
prompt traces to state them).

### 5.7 Discussion (0.5 page)

**Write.** What the benchmark discriminates now, stated as a finding: the pool is solved at the top on final
answers by models of September 2026, and for the top five the remaining incorrect verdicts are mostly form (34
symbolic items on one template scored by their numbers), prescribed-precision misses (17 templates require the
digits the question prescribes; the top models compute at full precision and land a unit off), and the two
under-specified chemical templates; the headroom is the Advanced tier, depth, prescribed-precision compliance and
the process scores (`docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` sections 5 and 6; `RESIDUAL_INCORRECT.md`). Then the
limit of verification as the paper's methodological point: a trace that computes correctly, misstates the rule it
applies and lands the right answer passes every deterministic check; only a judge sees it, at about a third; and
the final answer nearly decides the experts' trace verdict, so process scores earn their place on the wrong answers
and on the arithmetic behind right ones. Then what reasoning effort buys (nine points for one model, nothing
detectable for another) and what the flagship anchors say (the tier, not a jump). Then one paragraph on the
scoring statements the reader needs (nine symbolic templates scored by stated numbers; array items and two-quantity
lines verified on the last value; 84 multipart items checked on fewer quantities than the question's labelled parts;
the six surface-shortcut templates and the headline without them).

### 5.8 Conclusion (0.2 page)

One paragraph: what EngTrace is, what it found, and what comes next (the tool condition, the corrected open-book
arm, method memorisation, multi-modal artefacts), without "future work will explore" padding.

### 5.9 Limitations (not counted; required)

One paragraph each, from `RESULTS_PAPER_NOTES.md` "For the limitations" and `PARAPHRASE_PAPER_NOTES.md`: four
models without extended reasoning at their providers' defaults; 1.0% unusable traces concentrated on six templates
where the score also measures finishing within the ceiling; coverage credits stated intermediates and is blind to a
misstated rule; the digit rule's flags are 90% real on this roster but it finds under half of the slips the experts
marked; the top five are not separable at this size (150 templates detect AS differences of about 1 to 4 points and
level gaps of 7 to 15); the evaluator was validated on five models not on the roster and on 15 templates, with
roster-specific readings covering the answer check, the judge and the digit rule only; the two under-specified
chemical templates and the 17 exact-digit templates; the certification's Markdown rendering; the paraphrase test
covers 115 of 150 templates, tests wording not method, used one writer and a margin fixed after the point
estimates; items are examination-style by construction (what verifiability needs; ill-posed problems and real
artefacts are out of scope); the pilot's labels are not released (section 9, decision 5); the open-book arm as
run had defective equation blocks and the tool condition was not run (if that is still so at submission).

### 5.10 Ethics and data

Synthetic data from referenced constants; no private data; models should not be deployed in safety-critical
systems on this evidence; experts' labels and personal data stay out of the release; the pool and seed are
released at publication under the MIT licence; the annotators were engineers who certified the templates and read
the traces (compensation, if any, stated).

## 6. Main text against appendix

| Goes in the main text | Goes in the appendix |
|---|---|
| Figure 1 example; Figure 2 taxonomy; Figure 3 evaluation pipeline | Taxonomy procedure and prompts (A); textbooks and reference sources (B) |
| Section 3's certification numbers in one paragraph | Certification in full: gate checks and register, screen passes and agreement, Layer 2 protocol, plants, hand checks, agreement, panel FPR and MAD, round-by-round counts, the Markdown-rendering finding, fix-list summary (D) |
| Dataset statistics, one table | Templates by branch, domain, area, level, answer kind; pool freeze and selection rule; diversity (skeletons, paths, near-duplicates); the 58 single-path templates (C) |
| Authorship, one paragraph | The authoring specification summary for the 60 new templates and the review cycles they went through (E) |
| Section 4's evaluator description and Table 2 validation summary | Evaluator details: answer kinds and the comparator contract, milestone derivation, the digit rule, the judge prompt and its validation, the router rule and batching (F); the meta-evaluation in full: design, labels, agreement, every evaluator's figures with intervals, planted defects, routing, PRMs (G) |
| One sentence on thresholds | The threshold appendix as `THRESHOLD_APPENDIX.md` prints it (F) |
| One paragraph on judge independence | Judge selection and the exposure table; judge swap table; leave-one-judge-out and per-judge bias tables (G) |
| Analysis plan, one sentence | The plan with its tests and detectable differences; the exploratory additions, dated (H) |
| Roster and decoding, one paragraph | Appendix I from `DECODING_TABLE.md`: served ids, endpoints, ceilings, tokens, reasoning tokens, cost per model; the inference prompt and parsing (I) |
| Table 1; Figures 4 to 7; the results paragraphs | Pairwise table (55 rows) with detectable differences; instance variance; branch, level and domain tables with intervals and the branch pairs; answer types; sensitivity table; repeats; tokens against score; endpoints matched; classification by label (J) |
| Paraphrase paragraph | Paraphrase design, prompt, checks, survival by branch, the 35 lost templates, the experts' rejection reasons, the Q5 table (K) |
| Level-gap sensitivity sentence | The two chemical templates and the 17 exact-digit templates, with the experts' B4 answers; the residual-incorrect distance table (L) |
| Conditions in three sentences | C1, C3 and C4 tables as `RESULTS.md` prints them, with the arms' decoding tables and the open-book block defects (M) |
| Error analysis in three sentences | B1 to B3 in full: the answer-check reading with its branch and template breakdown, the judge reading, the error categories per model and level with kappa (N) |
| Scoring statements, one paragraph | The lists: nine symbolic, four array/two-quantity, 84 multipart items, six shortcut templates, the one known wrong credit; the answer-check corrections with their counts (O) |
| Related Work, five paragraphs | Table A and further related work (P) |
| — | Reproducibility: what is committed, what is released at publication, what is withheld and why; the anonymised archive's contents (Q) |

## 7. Figures and tables

| Item | Content | Source data | State |
|---|---|---|---|
| Figure 1 | A template (code excerpt), one instance, its gold trace, and the milestones and claims the evaluator derives | a certified template; `milestones.py` | the May CSTR figure exists (`figs/engtrace-example-template.pdf`); add the milestone and claim overlay, or pick a civil or industrial template to show a new branch |
| Figure 2 | Five-branch taxonomy | `figures_oct_12/engtrace-overview-5branch.pdf` | exists, untracked; copy to `figs/` |
| Figure 3 | The evaluation pipeline: answer check, milestones, digit rule, residual judge, router, with the judged fraction | this document, section 5.4 | to draw |
| Figure 4 | AS per model with intervals, tiers marked; or the 55 Holm verdicts as a matrix | `results/results.json` Q1 | script to write |
| Figure 5 | Easy / Intermediate / Advanced per model with intervals, Welch verdict marked | `results.json` Q2 and the level table | script to write |
| Figure 6 | Coverage on wrong-answer traces against the chance floor, per model | `results.json` Q3 | script to write |
| Figure 7 | Paraphrase change per model with 90% intervals and the ±5 band | `results.json` Q5 | script to write |
| Table 1 | AS, FSR, unusable, SD within, MC, AFR with claims per trace, SFR | `results.json` Q1, Q3 | script to write |
| Table 2 | Validation against experts: each component's agreement, precision, recall; planted-defect detection | `SCORER_VALIDATION.md`, `E5_VALIDATION.md`, `ROUTER_VALIDATION.md`, RESULTS_X1 | collate by script |
| Table 3 | Dataset statistics | `manifest.jsonl` | script to write (the counts are in this document's section 5.3) |

Every figure and table is produced by a committed script reading `results.json`, the manifest or the validation
files, and the script's commit is named in the appendix. Nothing is typed from a Markdown table. Figures are
drawn to read in greyscale; the ACL column is 3.03 inches wide.

## 8. Writing rules

1. **One voice, present tense, no history.** The paper describes what EngTrace is. "We" build, certify, score,
   validate. No "previously", "originally", "in response to".
2. **Every number from a committed script**, cited by file in the appendix; regenerate tables from `results.json`,
   never copy from `RESULTS.md` by hand.
3. **The unit of analysis is named once and kept**: the template; every interval resamples templates; every family
   of tests is Holm-corrected; every null carries its detectable difference.
4. **Tiers, not an order**, wherever intervals overlap; "not separable at this size", never "equal".
5. **Flags are flags; coverage is coverage.** The two naming rules of section 4.
6. **Say the setting with the number**: provider defaults, the four models without reasoning tokens, the ceiling,
   the unusable rule, wherever a closed-against-open or a level sentence appears.
7. **Name the limits where the claim is made**, not only in Limitations: the judged fraction beside coverage, the
   precision beside a flag rate, the two chemical templates beside the level gap, the 115 templates beside the
   paraphrase bound.
8. **Internal labels stay internal**: no E0 to E5, Q1 to Q5, D-numbers or file names in the text; one appendix table
   maps paper names to scripts.
9. **Honesty items that must appear** (D-114, D-169, D-171, D-183): the Markdown rendering during certification; the
   answer-check corrections after the run with their counts and the one known wrong credit; the two under-specified
   chemical templates; the 17 exact-digit templates; the open-book blocks' defects if that arm is reported as run;
   the four no-reasoning models; the validation population against the roster.
10. **Length discipline**: cut evidence to the appendix before cutting a limit from the main text.

## 9. Decisions the authors must make before the text is final

1. **Anonymity.** Link an anonymised archive (ARR zip or anonymous.4open.science), never the repository: it
   carries author names, the decision log and the pilot record. Cite FinChain in the third person; do not cite the
   preprint as the authors' own; remove the project and code links from the title block for the review version.
2. **The authoring record of the 60 civil and industrial templates.** The authoring specification is not on `master`
   (the pilot folder was removed on 2026-09-27; the templates arrived in commit `96448c7`, 2026-09-05). Commit a
   dated summary of how they were written (the specification, the AI drafting, the blind review cycles, who checked
   them) before submission, so that the paper's authorship paragraph points to a record. The paper must not say
   "fully authored by domain experts" of all 150.
3. **The taxonomy and difficulty procedure for the new branches.** Say what was done: the domains and areas came
   from an ABET / ASCE / IISE working taxonomy, and the difficulty tiers from the 2026-09-05 audit inventory; or run
   the expert-persona panel on the two branches now (minutes, negligible cost) and report its agreement. Either is
   fine; describing a procedure that was not run is not.
4. **The retrieval objection: the corrected open-book arm (`openbook2`, 249 items per model, $7.68 plus about $1 for
   the judge) and the tool condition (calibration about $1.50, then $20 to $54).** As this is written both arms
   appear to be running in another session (their trace reviews exist in the working tree, uncommitted). Once their
   rows are in `RESULTS.md`, the sentence in 6.5 is written from them in the form `RESULTS_PAPER_NOTES.md` gives
   ("For the retrieval objection"); if either arm does not finish, the open-book arm is reported as run with its
   defects stated, or left out.
5. **The experts' labels.** Decide whether an anonymised ground-truth file (codes and majority labels, no per-rater
   rows) is released with the traces and scores; otherwise the paper states that the labels are withheld and why,
   and that the validation figures reproduce only with them.
6. **The two chemical templates.** The experts' B4 answers support stating them as a limitation with the level-gap
   sensitivity row (3 of 3 say the virial wording does not decide the form; 2 of 3 say standard heat-capacity sources
   differ by more than the tolerance). The pool is frozen; no re-score.
7. **Metric names.** Confirm the names of section 4, and whether "Final Answer Accuracy" survives as the name of the
   fully-solved rate.
8. **LPDS.** Whether the Limitations carry the neutral sentence on randomly sampled instances (`CHANGES.md` section 4).
9. **Housekeeping that blocks the camera-ready claims**: `README.md` still says 90 templates and three branches;
   `figures_oct_12/` and `current_overleaf_project/` are untracked; Appendix I is generated from
   `DECODING_TABLE.md`; the pinned scoring-library versions are stated (D-081).
10. **Title**, as in section 2.

## 10. Order of work to 12 October

| Days | Work | Depends on |
|---|---|---|
| 4 to 5 Oct | Section 3 (benchmark, certification), Section 4 (evaluator, validation, independence), the Related Work drop-in, Limitations and Ethics; the appendix skeleton with the generated tables that need no new data (certification, thresholds, decoding, pool statistics) | nothing pending |
| 5 to 6 Oct | The figure and table scripts from `results.json` (Table 1, Figures 4 to 7, the appendix results tables); Figure 3 drawn; Figure 1 revised | nothing pending |
| 6 to 7 Oct | Section 5 (setup) and Section 6 (results) text; the conditions sentences in the form the decisions of section 9 allow; the error-analysis appendix from `EXPERT_REQUEST.md` | decisions 4 and 6 |
| 7 to 8 Oct | Introduction, abstract, contributions, Discussion, Conclusion; the authorship and taxonomy paragraphs | decisions 2 and 3 |
| 8 to 9 Oct | Appendices filled; cross-references; the bibliography regenerated from `related_work_sources.bib` plus the method citations; venue confirmations | — |
| 9 to 10 Oct | Two proofreading passes (one for numbers against `results.json`, one for the "do not" list and the naming rules); page budget; the anonymised archive; the responsible-NLP checklist | — |
| 11 Oct | Buffer; final compile in `review` mode; submit before 12 Oct AoE | — |

Nothing in the first two rows waits on a pending run. If `openbook2` or the tool arm runs, its figures drop into
Appendix M and one sentence in 6.5 changes.

## 11. Notes for the revision letter (to be written later)

The letter is where every "what changed and why" sentence lives. Its spine is a table: each reviewer point from
both cycles, where the paper now answers it, and what the answer is. The record for that table is
`docs/PILOT_AND_FULL_RUN_ASSESSMENT.md` section 4 and `docs/NEXT_CYCLE_REVIEW.md` section 1; the points below are
the ones the letter must handle with care.

1. **Why the numbers moved and cannot be compared.** The evaluator the May paper used was measured against 15
   experts' labels and found wrong on 72 of 300 traces (68 of them correct answers called wrong), ranking the
   strongest model fourth where the experts rank it first; its panel was two judges where the paper said three;
   its 20% wrong-answer sample reordered the models between identical runs; its routing showed a judge the flawed
   step as often as a clean one. It was replaced, not patched. The May pool is not reproducible: eight of the
   published templates' gold traces did not reproduce their own answers, 77 templates moved during the rework, and
   the roster changed. So no October number is placed beside a May number; the metric-continuity table of
   `RESULTS_PAPER_NOTES.md` says what became of each May column.
2. **Promises kept, in their new form.** Variance and significance (yAYU 2, 9W1B 4): per-template SD, template
   intervals, 55 Holm-corrected pairs, Welch on the level gap, detectable differences. Wilcoxon on the continuous
   reasoning score beside McNemar: the coverage comparison table. The judge-exclusion ablation with placebo and
   untouched controls (meta-review 1; yAYU 1; 9W1B 2; gFWV 4): replayed offline, no family effect; plus the per-judge
   bias, the judge swap on the new stack, and a judge from no roster family. Threshold sensitivity (yAYU 4, 9W1B 3):
   the threshold appendix; the cross-encoder and alignment-ratio thresholds no longer exist. Inter-judge agreement
   (9W1B): reported for the old panel (0.725 / 0.776) and the replacement panel (0.563 / 0.622). A fourth branch
   (ynoK 1; January AC): two added. Linguistic diversity (meta-review 3; 9W1B 1; cqGs 2; nWW3 2): the paraphrase
   arm, as a bound; the shortcut audit. Tool or retrieval baselines (cqGs 1 and 4; ynoK 3; January AC): the
   open-book arm and the tool condition, in whatever state they are at submission. Error analysis too coarse (nWW3
   3 and 4; cqGs 3): the automated attribution validated by error type on the experts' labels, and the human reading
   of 160 wrong answers. Validation in the main text (9W1B 6; meta-review 2): Section 4.3. Proofreading and the
   broken reference: moot, the text is new.
3. **Promises replaced with a reason.** Branch-level reporting for the math-specialised models and the
   math-pretraining claim (meta-review 2; yAYU 5; 9W1B 5): the claim is dropped; no math-specialised model is on the
   roster, which the roster rule and the budget set; the paper makes no such claim. An interval on rho = 0.632 and
   DeepSeek V4 Pro's category: those numbers no longer exist; the study was replaced by the 300-trace meta-evaluation.
4. **Objections answered by construction rather than by argument.** gFWV 1 (perfect kappa unconvincing): the new
   certification rejected 53 templates, caught 60 of 60 plants and records hand checks and timestamps. gFWV 3
   (instance-level verification shallow): the gate at 500 seeds, the gold validation of every pool item, the
   display-tie census. gFWV 5 (authorship): stated exactly for the 90 and the 60. gFWV 4 and 9W1B 2 (validated by
   the same kind of system): expert labels, deterministic components, a judge from outside the roster. nWW3 1 and 2
   (comprehensiveness, 90 problems): five branches, 150 templates, structural variation measured, the pool's
   second-path coverage. cqGs 3 (causes of the cliff): depth, attribution, the two-template sensitivity, the human
   reading by level.
5. **What the letter must disclose that a reader of the repository would otherwise find first**: the September
   audit's findings on the published templates, the evaluator's measured defects, the Markdown rendering during
   certification, the answer-check corrections after the run, and that the new numbers are high because the models
   of September 2026 solve the pool at the top on final answers.

## 12. Working with the Overleaf project

The local copy in `current_overleaf_project/` should be the source of truth, committed to this repository so every
edit is reviewable in git, and pushed to Overleaf rather than typed there. Overleaf offers three ways in:

1. **Git integration** (a premium feature; MBZUAI's institutional licence may provide it: Overleaf menu, "Sync",
   "Git"). Each project has a git URL (`https://git.overleaf.com/<project-id>`) and an authentication token
   (Account settings, "Git integration"). From this repository, `git subtree push --prefix=current_overleaf_project
   overleaf master` pushes the folder to the project's root, and `git subtree pull` brings co-authors' Overleaf
   edits back. One command each way, with history on both sides.
2. **GitHub sync** (also premium) needs the paper at the root of its own GitHub repository, so it would mean a
   second repository for the paper; not worth it beside option 1.
3. **Upload** (free accounts): Overleaf's "Upload" button accepts several files at once and overwrites same-named
   files, so the changed `.tex` and `figs/` files can be dropped in after each local commit; copy-paste is for a
   single small edit. A zip upload creates a new project, so it is not an update path. Co-authors' edits come back
   through "Download as zip" or the project history.

There is no TeX installation on this machine, so page counts and the final check are compiled on Overleaf
unless MiKTeX (pdfLaTeX, the engine Overleaf's ACL template assumes) is installed locally; the renderer in
`docs/render_md_pdf.py` is for documents like this one, not for the paper.
