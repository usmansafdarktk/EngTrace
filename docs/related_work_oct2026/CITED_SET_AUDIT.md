# Audit of the May 2026 citations

The question for each paper the May 2026 submission cites in its Introduction and Related Work: does the
sentence that cites it say something the paper supports, and does the paper still earn its place?

**Method.** Each paper was downloaded (`papers.json`, `fetch_papers.py`) and read by three reviewer
agents working independently from the paper's text, under `notes/review_protocol.md`. Their 138 reviews
are in `reviews/cited/`; `tally_reviews.py` tabulates them (`notes/review_tally.md`). The three verdicts
agree for 36 of the 46 papers, and every disagreement is one step on the scale (for example, imprecise
against misleading). The consolidated verdict below takes the majority, and the note says where the
reviewers split. Verdict scale: **accurate** (the paper supports the sentence), **imprecise** (right in
gist, loose in detail), **misleading** (suggests something the paper does not support, or omits what
matters most for the comparison), **wrong** (the wrong paper, or a claim the paper contradicts).

## Summary

Of the 40 works the May text cites, 11 are cited accurately, 13 imprecisely, 9 misleadingly and 7 wrongly.

- **The claim the section rests on does not hold.** The closing sentence says "existing benchmarks remain
  limited to outcome matching", and the Introduction says "no existing benchmark verifies this process".
  Seven of the papers the May text itself cites contradict this. PhysReason scores steps against
  annotated solutions. FinChain scores steps against executable gold chains. FEABench's main metrics are
  intermediate ones. TransportBench's experts graded reasoning. ElecBench scores reasoning with LLM
  judges. EngiBench uses rubrics for its open-ended level. EngDesign checks designs by simulation with
  partial credit. The revision claims the combination instead, hedged with "to our knowledge".
- **Four citations point at the wrong paper.** "Glue" cites SuperGLUE, "BBH" cites BIG-Bench Extra Hard,
  "BIG-Bench" cites an unrelated text-to-image bias benchmark, and the "Liu et al., 2025" entry merges a
  homework-grading study with the CIRCUIT benchmark under an author who wrote neither.
- **Three characterisations are the opposite of the paper.** ABench-Physics is called "advanced symbolic"
  but scores only numerical answers. LLM-SRBench is called "qualitative" but is symbolic regression scored
  by equation recovery. FEABench is called outcome-only, but its main metrics are intermediate.
- **The closest precedents are hidden.** FinChain, cited once for "symbolic templates", already pairs
  generated instances with executable gold traces and scores steps against them. CIRCUIT, never named,
  already builds an engineering benchmark from templates. PhysReason already scores steps with a
  human-checked scorer. The revision names all three and says what EngTrace adds.
- **The "confined to theoretical physics" sentence fits none of the six physics papers.** Five are
  textbook, exam or olympiad physics. NEWTON is everyday physical commonsense. LLM-SRBench is mostly not
  physics.
- **Two self-citations do not support their sentences.** A storytelling study is cited for the importance
  of evaluating reasoning. A reference-free story metric is cited for "standard reference-based metrics",
  and its own results rank reference-based metrics lowest. Both are dropped.

## Paper by paper

`v3` is the revised text (`RELATED_WORK_v3.md`). Page numbers refer to the text files in `text/`.

| Paper (as cited) | Where | Verdicts r1 / r2 / r3 | Consolidated | In v3 |
|---|---|---|---|---|
| GSM8K, Cobbe et al. 2021 | RW math | imprecise ×3 | imprecise | kept; cited for final-answer scoring that credits flawed reasoning (p.7) |
| MATH, Hendrycks et al. 2021b | Intro, RW math | imprecise ×3 | imprecise | kept; it ships step-by-step solutions but scores the boxed answer |
| HARDMath, Fan et al. 2024 | RW math | misleading ×3 | misleading | kept as a generation-with-solutions precedent (applied, not "abstract" math; code-generated, numerically validated solutions) |
| Putnam-AXIOM, Gulati et al. 2024 | RW math | imprecise ×3 | imprecise | kept for functional variants; cite the ICML 2025 version (the paper contradicts itself on o1-preview's drop, so quote no figure) |
| GSM-Symbolic, Mirzadeh et al. 2024 | Intro, RW math | accurate ×3 | accurate | kept; cite ICLR 2025; its templates also emit a worked solution |
| HumanEval, Chen et al. 2021 | Intro, RW code | accurate ×3 | accurate | kept |
| MBPP, Austin et al. 2021 | RW code | accurate ×3 | accurate | kept |
| SWE-bench, Jimenez et al. 2024 | RW code | accurate ×3 | accurate | kept |
| Heesch et al. 2025 | RW math | misleading ×3 | misleading | moved to the engineering paragraph: it never tests transfer from math or code; it is itself an engineering evaluation and criticises exam-style items (p.1) |
| UGPhysics, Xu et al. 2025 | RW physics | imprecise ×3 | imprecise | kept as final-answer physics (rules plus a model judge); ICML 2025 |
| PhysReason, Zhang et al. 2025a | RW physics | misleading ×3 | misleading | kept as a step-scoring precedent (first-error scoring checked on 1,000 annotated solutions) |
| ABench-Physics, Zhang et al. 2025b | RW physics | wrong / wrong / misleading | wrong | kept as numerical-answer physics; "symbolic" removed |
| PHYBench, Qiu et al. 2025 | RW physics | imprecise ×3 | imprecise | kept as expression-scored physics |
| NEWTON, Wang et al. 2023 | RW physics | misleading ×3 | misleading | moved to the appendix (template-generated commonsense, not theoretical physics) |
| LLM-SRBench, Shojaee et al. 2025 | RW physics | wrong ×3 | wrong | dropped (2 of 3 reviewers; the third would keep it only for contamination) |
| MMLU, Hendrycks et al. 2020 / 2021a | Intro, RW eng. | imprecise ×3 | imprecise | kept; the two reference entries are one paper (ICLR 2021) and are merged |
| BIG-Bench, cited as Luo et al. 2024 | RW eng. | wrong ×3 | wrong | dropped; the cited paper is a text-to-image bias benchmark named BIGbench |
| SuperGPQA, Du et al. 2025 | RW eng. | misleading ×3 | misleading | kept, reworded: 7,892 engineering questions, many requiring calculation, scored on the option chosen |
| TransportBench, Syed et al. 2024 | RW eng. | accurate ×3 | accurate | kept; its experts graded the reasoning too |
| APBench, cited as Chen et al. 2025 | RW eng. | accurate ×3 | accurate | kept as Wu et al. 2025 (the first author is Di Wu) |
| ElecBench, cited as Cheng et al. 2024 | RW eng. | accurate / imprecise / imprecise | imprecise | kept as Zhou et al. 2024 (the first author is Xiyuan Zhou) and now named |
| EEE-Bench, Li et al. 2025 | RW eng. | imprecise ×3 | imprecise | kept; a whole branch in ten subdomains, so not "narrow" |
| "Liu et al. 2025" (CIRCUIT) | RW eng. | wrong ×3 | wrong | replaced by CIRCUIT, Skelic et al. 2025 |
| EngDesign, Guo et al. 2025b | RW eng. | imprecise ×3 | imprecise | kept; design is checked by simulation, not called "tangential"; NeurIPS 2025 |
| FEABench, Mudur et al. 2025 | Intro, RW eng. | wrong / misleading / wrong | wrong | kept as multiphysics simulation; no longer called outcome-only |
| EngiBench, Zhou et al. 2025 | Intro, RW eng. | misleading / misleading / wrong | misleading (RW), wrong (Intro) | kept, reworded: two LLM evaluators score most problems, rubrics only the open-ended level; Findings of ACL 2026 |
| FinChain, Xie et al. 2025 | Intro | misleading ×3 | misleading | promoted to the closest precedent and discussed in the third person; ACL 2026 |
| Zhao et al. 2023 (survey) | Intro | imprecise ×3 | imprecise | Intro only; better cited for correct answers via invalid reasoning (v19 p.61) or not at all |
| Xie et al. 2023a (storytelling) | Intro | misleading ×3 | misleading | dropped |
| Xie et al. 2023b (DeltaScore) | §4.3 | misleading / misleading / wrong | misleading | dropped |
| SuperGLUE, Wang et al. 2019 | Intro, as "Glue" | wrong ×3 | wrong | Intro: cite GLUE (Wang et al. 2018) for the name, SuperGLUE as evidence that GLUE saturated |
| BBEH, Kazemi et al. 2025 | Intro, as "BBH" | wrong ×3 | wrong | Intro: cite BBH (Suzgun et al. 2023) for the name, BBEH as evidence that BBH saturated |
| Ott et al. 2022 | Intro | accurate ×3 | accurate | kept (now with Akhtar et al. 2026 for LLM-era saturation); templates are not claimed to prevent saturation |
| Deng et al. 2024 | Intro | accurate ×3 | accurate | kept |
| Plank 2022 | §3.3, App. K | imprecise ×3 | imprecise | only if the appendix discusses leftover disagreement on subjective ratings |
| Let's Verify Step by Step, Lightman et al. 2023 | §4 | accurate ×3 | accurate | kept in the new process paragraph; ICLR 2024; author list regenerated (the May entry lists 30 names for 10 authors) |
| FrugalGPT, Chen et al. 2023 | §4 | imprecise ×3 | imprecise | method section only, as an LLM cascade; TMLR 2024, not ICML |
| PoLL, Verga et al. 2024 | §3.3 | imprecise / imprecise / accurate | imprecise | kept for panels of judges; not cited in support of the single residual judge, the arrangement it argues against |
| BERTScore, Zhang et al. 2020 | §4.3 | accurate ×3 | accurate | appendix only, if the similarity metrics are kept there |
| ROUGE, Lin 2004 | §4.3 | accurate ×3 | accurate | as BERTScore; the one citation covers ROUGE-L too |

Works the May text names but attributes to another reference, and two related works:

| Paper | Status | Verdicts | In v3 |
|---|---|---|---|
| GLUE, Wang et al. 2018 | named "Glue" in the Intro, not in the references | describes the named benchmark (2 of 3) | add to the Intro if GLUE stays named |
| BBH, Suzgun et al. 2023 | named "BBH" in the Intro, not in the references | describes the named benchmark (2 of 3) | add to the Intro if BBH stays named |
| BIG-bench, Srivastava et al. 2022 | named "BIG-Bench" in RW, not in the references | wrong / wrong / misleading: not multiple-choice recall, no engineering task | not cited |
| EngiBench, Felten et al. 2025 | named in a note under Zhou et al. | the note is accurate; it evaluates no language model | not cited; delete the note |
| SciBench, Wang et al. 2024 | cited in the July rebuttal only | accurate ×3 (5% relative tolerance on the textbook set) | added to the physics sentence |
| CIRCUIT, Skelic et al. 2025 | the arXiv id inside the "Liu et al." entry | misleading by omission: 102 templates, five numerical setups each | added; a templated engineering precedent |

## Reference-list corrections

`related_work_sources.bib` regenerates every entry the revision cites from arXiv, Crossref or the ACL
Anthology, so these are corrected there. Listed for the authors' check of the rest of the bibliography:

- **Wrong first author:** ElecBench (Zhou, not Cheng; the May entry prints the author list of a different
  Scientific Reports paper by Cheng et al.); APBench (Wu, not Chen; Jang and Lavezzi are not authors);
  the "Liu et al." entry (no author named Liu wrote either paper).
- **Wrong author list:** Lightman et al. (30 names printed for 10 authors; Baker, Lee and Cobbe missing).
- **Wrong given names:** ABench-Physics (12 of 12), BIGbench (9 of 9), UGPhysics (8 of 9), PhysReason
  (6 of 9), BBEH (6 of 20), LLM-SRBench (4 of 6), PHYBench (3 of 3), GSM-Symbolic (Keivan, not Kiarash;
  Hooman, not Hidetoshi), NEWTON (Jiafei, not Jia), SWE-bench (Shunyu, not Shun).
- **Venues now published** (cited as arXiv in May): GSM-Symbolic (ICLR 2025), HARDMath (ICLR 2025),
  Putnam-AXIOM (ICML 2025), UGPhysics (ICML 2025), PhysReason (ACL 2025), PHYBench (NeurIPS 2025
  Datasets and Benchmarks), SuperGPQA (NeurIPS 2025 Datasets and Benchmarks), EngDesign (NeurIPS 2025
  Datasets and Benchmarks), EngiBench (Findings of ACL 2026), FinChain (ACL 2026), Lightman et al. (ICLR
  2024), BBH (Findings of ACL 2023), SciBench (ICML 2024). FrugalGPT is TMLR, not ICML.
- **Duplicate:** MMLU appears as both Hendrycks et al. 2020 and 2021a.
- **Lost capitals:** about 20 titles print "Gpt-5", "llms", "Elecbench" and similar. The regenerated
  entries protect their titles.

## Reviewer disagreements, resolved

- **FEABench, EngiBench, ABench-Physics, DeltaScore:** two reviewers chose the more severe verdict, and
  that verdict is taken. In each case the facts were not in dispute, only the label.
- **ElecBench, PoLL:** two of three said imprecise, and imprecise is taken. ElecBench is never named in
  the May text, and its logicality metric scores reasoning. PoLL argues for small, diverse judges, and the
  May Tribunal used three frontier models.
- **GLUE, BBH:** one reviewer graded the May attribution (wrong) and two graded the named paper
  (accurate). Both views lead to the same action.
- **LLM-SRBench:** dropped on two of three recommendations.
- **CIRCUIT:** misleading by omission for two reviewers, imprecise for one. All three agree it must be
  named as a templated engineering precedent.
