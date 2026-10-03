# EngTrace as it stands: a brief for the Related Work rewrite

Written 2026-10-01 from `master` at `11e182c`, read-only. Each paragraph names its sources, and numbers are quoted, not rounded. "Counted" marks a count I made from the named files. The May 2026 paper's figures (90 templates, three branches, 1,350 items, 27 models) are history.

## 1. The benchmark

**Taxonomy** (counted from `data/templates/branches/` and `docs/re-implementation-sep/audit/template_inventory.csv`). Five branches, 15 domains, 150 templates and 42 area modules (one `.py` file per area):

| Branch | Domains (templates) | Areas | Templates E / I / A |
|---|---|---:|---|
| Chemical | reaction kinetics (10), thermodynamics (12), transport phenomena (8) | 7 | 13 / 9 / 8 |
| Civil (new) | geotechnical (10), structural analysis (10), water resources (10) | 11 | 12 / 12 / 6 |
| Electrical | digital communications (10), electromagnetics and waves (10), signals and systems (10) | 7 | 11 / 12 / 7 |
| Industrial (new) | production and inventory (10), quality and reliability control (10), stochastic operations (10) | 11 | 12 / 12 / 6 |
| Mechanical | fluid mechanics (10), mechanics of materials (10), vibrations and acoustics (10) | 6 | 10 / 13 / 7 |

By answer type, the 150 templates are scalar 89, multipart 32, symbolic 9, vector 8, array 7 and classification 5. `figures_oct_12/` (untracked) holds the three-branch overview figure, a regenerated five-branch version (`engtrace-overview-5branch.pdf`), its generator `make_overview5.py`, and 15 domain icons.

**Items** (`full_run_28092026/README.md`, D-114, D-116, `DIVERSITY.md`). Each template has 15 instances, giving 2,250 items: 870 Easy, 870 Intermediate, 510 Advanced, 450 per branch, all distinct questions. The 15 are taken round-robin across the reasoning-path and answer-form groups in the template's first 100 draws. Still, 54 to 58 templates follow one reasoning path in all 15 instances, and 38 have a single question skeleton (the wording with numbers masked).

**Private seed, hashed pool** (D-114, `full_run_28092026/README.md`). The pool is drawn from a private 128-bit seed. Only two files are committed: `manifest.jsonl`, with a per-item SHA-256 over question + NUL + solution and no text, and `FREEZE.json`, which commits to the seed's SHA-256. Reviewers get the pool in the ARR supplementary archive, and the seed is revealed at publication. The design "lets the paper say the evaluated items were not public when the models ran" (D-114).

**How the 60 new templates were written.** `docs/pilot_template_authoring_spec.md` is not on master. It exists only at commit `63bb6ac` (2026-08-18) on the local branch `pilot/phase1-template-authoring`, which has no remote-tracking ref. DECISIONS.md has no authoring entry, and the templates reached master on 2026-09-05 (`96448c7`). The spec (v1.0, "Draft for supervisor review") prescribes:

- **Scope.** AI agents work one branch at a time, from a taxonomy "derived from ABET curricula and professional-society standards: ASCE, IISE, AIAA". Each domain gets 10 templates: 4 Easy, 4 Intermediate, 2 Advanced.
- **Sources.** A Librarian agent picks the textbooks behind a human approval gate. Constants are cited, and an agent checks at least 20% of them.
- **Review.** Three fresh-context reviewer agents cover plausibility, an independent re-solve (within 0.1%) and pedagogy. A template passes only when every reviewer scores at least 4 on every dimension with no blocking flag; after 4 cycles it is discarded. An auditor re-labels difficulty blind.
- **Boundary.** Pilot templates "MUST NOT be reported in the paper unless they later pass the real Phase 2 (AI Tribunal) and Phase 3 (human certification) pipeline".

Branch reports there: civil "~50 template-review cycles" and 15 genuine gate failures; industrial "65 logged review cycles", one discard at the cap, two overruns, two cycles without a reviewer panel.

*Corrected 2026-10-03 (the owner):* the draft spec above describes a pilot that did not produce the shipped templates. The 60 civil and industrial templates were written by a colleague of the authors and certified by the own-branch experts like the other 90; the certification files cited in the "Disagreement" item below are consistent with that (every civil and industrial expert approved all 30 in round 1). The paper and the revision letter must not describe the new branches as AI-drafted or spec-generated; `docs/NEXT_CYCLE_REVIEW.md`'s passages saying so are corrected in place.

**Certification** (`template_annotation_23092026/README.md`, `layer2/RESULTS*.md`, `layer2/CERTIFICATION.md`, D-092 to D-119):

- **Layer 0**, a deterministic gate at 500 seeds: all 150 pass, after 54 were edited for closure.
- **Layer 1**, the Appendix H prompt with three judges from outside the roster (Grok 4.6, MiniMax M3, MiMo-V2.5-Pro). Pass 1 gave 126 pass, 13 controversial and 11 critical failure (AC1 0.836, Fleiss 0.287); 34 of the 45 claims were confirmed and fixed. Pass 2 gave 147, 3 and 0 (AC1 0.929).
- **Layer 2**, three own-branch experts per template, 15 in all, each solving one instance by hand first, with four planted defects per branch.
  - Round 1: 60 of 60 plants caught, 459 of 504 hand checks matched, AC1 0.914 (Fleiss 0.589), 22 templates rejected. The screening panel's false-positive rate against the experts: 10.2%.
  - Rounds 2 to 4 re-reviewed 22, 5 and 2 templates. All 150 are now certified unanimously, after 53 rejections (43 plus 10).
- **Caveat.** The experts judged 65 templates through a Markdown view that dropped `*` and `$` (D-108). Certification closed without a plain-text look (D-114).
- **Disagreement.** NEXT_CYCLE_REVIEW §4.3 says the 60 new templates went through "three rounds of own-branch expert certification that rejected 22 of them in round 1 and required 25 fixes". The certification files say otherwise:
  - `layer2/RESULTS.md` shows every civil and industrial expert approving 30 of 30 in round 1 (AC1 1.000).
  - `CERTIFICATION.md` shows all 60 last reviewed in round 1.
  - The 22 round-1 rejections fall in chemical (9), mechanical (10) and electrical (3) (counted).
- **Tiers** come from the 2026-09-05 inventory (§9.3 item 10). The Layer 2 guide and app never mention difficulty, so experts apparently did not see tiers.

## 2. The evaluation stack

Sources: D-105, D-110, D-113, D-142, `full_run_28092026/ANALYSIS_PLAN.md`, `results/RESULTS.md`, `GOLD_VALIDATION.md`.

| Component | Kind | What it does | State |
|---|---|---|---|
| Answer check | deterministic | Correct, partial or incorrect, against quantities the gold computed. Tolerance REL 0.002 (D-137) or one unit of the last digit. Partial scores 0.5 (D-117). Symbolic answers are scored by the numbers they state (D-138). | Gold: 2,250 of 2,250 correct, six answer types |
| E3 milestones | deterministic | Values the template computed, the gold states and the question does not give. Tolerance 0.5%, 13 unit factors, order-free (D-084). | 70 items have none (D-135) |
| E4 digit rule | deterministic | Recomputes each claim from the numbers the trace shows; flags a displayed value that is not a correct rounding (D-097, D-101). Parser fixed three times (D-156, D-159, D-160). | Gold: 0 of 7,167 claims flagged |
| E5 residual judge | LLM | MiMo-V2.5-Pro, only on milestones E3 missed. E5-strict credits REACHED. | 8,032 calls, $36.3016, all answered (D-155) |
| Step router | LLM | Steps the digit rule does not flag go to MiMo, batched per trace (D-113). | Running (D-164) |

**Judged fraction.** On the pilot, 23.6% of milestones went to the judge (`RESULTS_E5.md`). On the full run the judged fraction is 0.112 to 0.252 per model, and 0.000 to 0.004 is left unjudged (RESULTS.md Q3). E2's PRMs are not in the stack, and E0 is not run at scale (D-105).

**Independence.** The roster rule: "no model the pilot generated traces with, and no model used as a judge" (D-110, D-132). MiMo (Xiaomi) shares no family with the roster, but "Independence is a matter of degree" (`JUDGE_SELECTION.md`): MiMo's documented exposure is to Claude, "0.4M exchanges, the smallest named", and Claude Sonnet 5 is on the roster.

## 3. Validation evidence

**X1** (`evaluator_pilot_17092026/RESULTS_X1.md`, D-079, D-099, D-137):

- **Design.** 300 traces: 60 items × 5 models (GPT-5, Claude Opus 4.7, Gemini 3.1 Pro, DeepSeek R1, Llama 3.1 70B). 15 experts, three from the trace's own branch per trace.
- **Agreement.** Fleiss κ between raters: 0.781 on steps, 0.880 on milestones, 0.966 on the answer. Cohen's κ within raters: 0.828.
- **Answer check.** 0.947 agreement on the 281 non-partial traces (the published check: 0.747), and 0.893 three-way. After D-137: 0.982 and 0.927, the figures to cite.
- **Milestones.** E5 F1 0.958 (precision 0.930, recall 0.989); E3 F1 0.921.
- **Trace AUROC.** E5 0.886 (0.835 to 0.934 resampling traces, 0.784 to 0.978 resampling templates). E0 0.850 (0.804 to 0.892; 0.760 to 0.927). The experts' own answer verdict scores 0.974. E5 minus E0 is +0.036 (template interval −0.093 to +0.173), against a detectable difference of 0.192. Design effect: 2.6 to 4.2.
- **Hard case.** 93 of 228 correct-answer traces contain an incorrect step: 178 steps, 175 calculation slips and 3 conceptual. E4 scores AUROC 0.661 against E0's 0.542, a gain of +0.149 (+0.000 to +0.300; `PILOT_SUMMARY.md` §2.5). The fixed rule reaches precision 0.817 and recall 0.427 there (RESULTS.md Q3).

**Judge overlap** (`RESULTS_E1.md`). With non-suite judges (Grok 4.6, MiniMax M3, MiMo), the pooled score goes from 0.379 to 0.378, and GPT-5's traces from 0.474 to 0.470. Fleiss κ among the judges: E0-3J 0.725 four-way and 0.776 binary; E1 0.563 and 0.622.

**PRMs** (RESULTS_X1 Finding 3, `RESULTS_E2.md`, D-090, D-100). Qwen2.5-Math-PRM-72B reaches step AUROC 0.825 (0.798 to 0.850), with precision 0.539 and recall 0.515. Inside correct-answer traces those fall to 0.246 and 0.255. Qwen2.5-Math-PRM-7B reaches 0.767. VersaPRM, a LoRA on Llama-PRM800K, reaches 0.629 with recall 0.064, and passed 12 of 13 known slips. A held-out threshold moves F1 by only −0.018 and +0.006.

**Planted defects** (`PILOT_SUMMARY.md` §2.7 to §2.9, D-102 to D-104, D-111, D-113). The set is 60 arithmetic, 60 conceptual and 60 control traces.

- **No deterministic check** catches any of the 60 conceptual plants.
- **The digit rule.** The bare rule catches 45 of 45 inside a parseable claim and 0 of 15 outside one. As E4 shipped it, the figures were 40 of 45 and 1 of 15; these predate D-156.
- **Judges, one step per call:**

  | Judge | Conceptual | Arithmetic | Clean steps flagged |
  |---|---|---|---|
  | GPT-5 | 20 of 60 (0.333) | 43 of 60 | 0 of 120 |
  | Opus 4.5 | 8 of 60 (0.133) | 43 of 60 | 0 of 120 |
  | MiMo | 16 of 51 (0.314) | 43 of 60 | 0 of 111 |

- **The published routing** sends a judge the corrupted conceptual step as often as the clean one (0.500 against 0.500). End to end it catches 9 of 60 (0.150).
- **The router.** One step per call, it catches 15 of 51 (0.294). Batched, it catches 19 of 60 (0.317), with 2 false alarms on 842 clean steps. On the labelled traces it finds 234 of 388 incorrect steps (precision 0.707, recall 0.603), against the digit rule's 0.825 and 0.255.

**Flag reading** (D-154, D-158, D-163). One domain expert read samples of the digit rule's roster flags:

- **Before the fixes:** 111 of 220 real, precision 0.505 (Wilson 0.439 to 0.570).
- **After the first fix:** 155 of 206, 0.752 (0.689 to 0.806).
- **Final rule:** 171 of 189 decided, 0.905 (0.854 to 0.939); per model 0.632 to 1.000.

## 4. The full run

**Roster** (D-110, D-132, `models.json`). Eleven models, 24,750 traces; qwen3.8-27b was set aside.

- **Open weight:** gpt-oss-20b, gemma-4-26b-a4b-it, deepseek-v4.1-flash, qwen3-235b-a22b-2507, glm-5.3-flash, glm-5.3, muse-glimmer-30b, kimi-k3. The 2507 release replaced qwen3-235b-a22b, which no endpoint served under the routing rule.
- **Closed weight:** gpt-5.4-mini, gemini-3.1-flash-lite, claude-sonnet-5.

There is no flagship and no math-specialised model. D-111 calls the closed models "the small closed tiers", not "frontier". Decoding uses provider defaults with a 32,768-token ceiling, 16,384 for Muse (D-121, D-122, D-133).

**Headline** (`results/RESULTS.md` Q1). Answer scores run from 0.805 (gpt-oss-20b) to 0.970 (deepseek-v4.1-flash); D-140's 0.799 to 0.969 predates D-147. 31 of 55 pairwise differences hold at Holm p < 0.05, and the top five lie within 0.012 of each other. SD is 0.041 to 0.213 within templates and 0.117 to 0.276 between them.

**Cliff** (Q2, D-146). Easy-minus-Advanced gaps run from +0.050 to +0.194. By Welch with Holm, 3 of 11 hold (gpt-oss-20b, gpt-5.4-mini, gemini-3.1-flash-lite); the planned permutation gives 6 of 11, and D-140's "7 of 11" is superseded. Without unusable rows, glm-5.3's gap falls from +0.133 to +0.009.

**Paraphrase arm** (D-157, D-161, D-162, D-165, `PARAPHRASE.md`, `PARAPHRASE_PAPER_NOTES.md`):

- **Items:** 450, the 1st, 6th and 11th instance of each template.
- **Writer:** Mistral Large 3 (`mistralai/mistral-large-2512`), from a family on neither the roster nor the judge side, using the fourth prompt; notation was restored deterministically.
- **Script checks:** the same numbers as a multiset, technical tokens, part labels, word similarity at most 0.75, length 0.7 to 1.5 times the original, no preamble. 316 passed; 29 templates lost all three items.
- **Expert check:** own-branch experts kept 277 pairs over 115 templates and rejected 39 (12%), mostly over a dropped qualifier.
- **Results:** paraphrase minus original runs from −0.031 (qwen3-235b-a22b-2507) to +0.016 (muse-glimmer-30b), and none survives Holm. 10 of 11 models stay within ±0.05 at 90%, a margin fixed after the point estimates. E5-strict differences run from −0.014 to +0.018. Kendall's τ is 0.636 (0.457 to 0.871), against a noise floor at the arm's size with median 0.783 and 5th percentile 0.636. Cost: $39.15.

**Decoding repeats** (D-151, RESULTS.md). gemma-4-26b-a4b, gpt-oss-20b, qwen3-235b-a22b-2507 and gemini-3.1-flash-lite each ran three times on 300 items. SD was 0.003 to 0.016 and the range 0.005 to 0.032; 0.850 to 0.937 of items got the same verdict every time. The other seven models have no decoding figure.

**Step router** (D-164, the latest entry). It was approved at a $140 cap and launched 2026-09-30T21:21:29Z; status RUNNING. The router column in RESULTS.md is empty.

**Not run (no D-entry):** the tool condition, the open-book condition, the flagship anchor, LOJO/X2, and the Appendix A/B/D taxonomy panel for the new branches.

## 5. Claims dropped and kept, and the prescribed framing

**Dropped:**

- **Math pretraining.** "the roster has no math-specialised model (D-110), so the paper does not make it" (`ANALYSIS_PLAN.md`). The July meta-review had asked to "Strengthen the causal claim about mathematical pretraining with branch-level reporting" (`docs/meta-reviews.md`, untracked).
- **Decoupling.** "Do not carry the sentence over" (NEXT_CYCLE_REVIEW §9.3 item 13). The Wilcoxon test "was tied to the decoupling claim the paper no longer makes" (D-149).
- **BERTScore and ROUGE.** "Keep ROUGE and BERTScore in an appendix if at all" (§3.5). Pilot BERTScore spanned only 0.866 to 0.882 (`RESULTS_E0.md`).
- **The old human-alignment result.** "Table 12's rho = 0.632 disappears, and with it the n = 100 study" (§3.2).
- **Old-to-new comparisons.** "the old pool is not reproducible" (§4.3).

**Kept.** Section 4 covers mechanism, validation and the limits of verification (§4.1, §4.2), framed not as "our evaluator beats the old one" but as "the old one was measured against experts, found wrong on a quarter of answers, and replaced by one whose deterministic components and residual judge were validated against the same experts, on 15 templates and five models" (§1).

**Prescribed framing:**

- **Scope.** "The abstract's 'stress-test generalization across diverse physical scenarios' becomes physical and structural diversity" (§4.3). "The abstract and section 3.2 keep the distinction between physical, structural and linguistic diversity, and the Limitations section says what the paraphrase test did and did not show. Multi-modal and real-artifact seeds stay future work." (§6).
- **Authorship.** "Writing 'authored by domain experts' over 150 templates would be false" (§4.3).
- **Contamination.** "The paraphrase test (section 6) is also the direct test of wording memorisation, and the paper should say so; the private-pool rule in 9.1 item 3 is the other half." (§9.3 item 6). The paraphrase notes rule out "no contamination" (say "no evidence of wording-level memorisation"), "the models are robust to paraphrasing" and "the ranking is stable".
- **RAG.** "RAG: do not add it." (§5). Its reason rests on the May paper's Table 18 (old data).
- **Tools.** "Tool use: add a code-execution condition on a subset." (§5). §9.4 makes it the first cut; it has not been run.

## 6. Positioning implications for Related Work

**What NEXT_CYCLE_REVIEW.md asks for** (verbatim):

- §4.2: "This positions the paper in the process-supervision literature (ProcessBench, Math-Shepherd, PRM800K, the self-preference results) and Related Work needs a paragraph on PRMs, LLM-judge validity and meta-evaluation, which the May version lacks."
- §9.4: "Start writing sections 3.3 and 4, related work, limitations and appendices K, L, M, O and P now: none of that text depends on the run."
- §9.6: "No RAG; the tool condition if it fits, the open-book condition as the fallback answer to cqGs 4; the related-work paragraph on PRMs and judge validity; the honest account of authorship, the taxonomy procedure, the certification rejections and why the numbers moved."

**The May paper** (`docs/_ARR_May__EngTrace.txt`). Its §2 covers math and coding, physical-science and engineering benchmarks. It cites PoLL (Verga et al., 2024) for the certification Tribunal, now only the Layer 1 screen (D-093), and Lightman et al. (2023) for process supervision. It names no PRM.

**What the repository's `.md` files mention.** I searched 161 files (tracked, untracked and gitignored, excluding `node_modules`, `.venv`, `__pycache__` and `.pytest_cache`):

| Name | Hits | Where, with context |
|---|---:|---|
| ProcessBench, Math-Shepherd | 1 each | `NEXT_CYCLE_REVIEW.md:267`, the §4.2 sentence |
| PRM800K | 4 | `NEXT_CYCLE_REVIEW.md:268`; `DECISIONS.md:3297`, `RESULTS_E2.md:21` and `hpg/README.md:25`, as VersaPRM's base, `UW-Madison-Lee-Lab/Llama-PRM800K` |
| VersaPRM | 29 | `DECISIONS.md` (D-090, D-100, D-111, D-112), `RESULTS_E2.md`, `RESULTS_X1.md`, `hpg/README.md`, `NEXT_CYCLE_REVIEW.md:264` ("finds 6% of incorrect steps") |
| Qwen2.5-Math-PRM | 11 | `DECISIONS.md` (D-090), `hpg/README.md` (pinned repos), `RESULTS_E2.md`, `RESULTS_X1.md` |
| self-preference | 4 | `NEXT_CYCLE_REVIEW.md:268`; `JUDGE_SELECTION.md:53`; `PILOT_SUMMARY.md:197` and `RESULTS_E1.md:29` ("No aggregate self-preference is visible") |
| panel of | 1 | `PILOT_SUMMARY.md:46`: E0 is "a panel of LLM judges" |
| GSM-Symbolic, EngiBench, SciBench | 1, 1, 2 | Only in the writers' own files from today in `docs/related_work_oct2026/` |

Zero hits: PRMBench, Skywork, Panickssery, Verga, PoLL, G-Eval, MT-Bench, FinChain, Putnam, MATH-Perturb, UGPhysics, PhysReason, PHYBench, FEABench, TheoremQA, OlympiadBench, SeePhys, PhysicsEval, ReasonEval, ROSCOE and "LLM-as-a-judge". "judge" has 783 hits in 72 files, all about the team's own judges; most are in `DECISIONS.md` (197), `RESULTS_X1.md` (68) and `PILOT_SUMMARY.md` (59).

**What this means for the section:**

- **Judge bias.** The record cites one paper, at `JUDGE_SELECTION.md:53`: "judges are over 50% more likely to wrongly pass a rubric item when the output is their own family's, even under objective criteria ([Self-Preference Bias in Rubric-Based Evaluation](https://arxiv.org/abs/2604.06996))".
- **PRMs.** The PRM paragraph can rest on the measured Qwen2.5-Math-PRM and VersaPRM results (§3). ProcessBench, Math-Shepherd and PRM800K-as-dataset appear only in the plan's sentence.
- **Meta-evaluation.** Position the team's own evidence against judge-validity and process-error benchmarks: expert labels with measured agreement, template-clustered power, and planted defects.
