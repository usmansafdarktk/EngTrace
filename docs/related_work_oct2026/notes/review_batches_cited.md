# Review batches: the cited set

Six batches, three independent reviewers each (r1, r2, r3). The sentences are quoted from
`docs/_ARR_May__EngTrace.txt` with its line numbers (spaces restored); `citation_context_map.md`
Part B has the full context. "Also used at" lists any other place the paper is cited.

## Batch M: mathematics and coding (Related Work paragraph 1)

Sentence S-M1 (L155-161): "In mathematics, benchmarks ranging from GSM8K (Cobbe et al., 2021) and
Math (Hendrycks et al., 2021b) to HardMath (Fan et al., 2024), Putnam-Axiom (Gulati et al., 2024),
and GSM-Symbolic (Mirzadeh et al., 2024) rigorously test abstract logical deduction and reasoning
stability."

Sentence S-M2 (L161-165): "Similarly, coding evaluations span from function-level generation in
HumanEval (Chen et al., 2021) and MBPP (Austin et al., 2021) to real-world software engineering in
SWE-Bench (Jimenez et al., 2024)."

Sentence S-M3 (L165-168): "Yet neither abstract logical deduction nor algorithmic proficiency
transfers to reliable reasoning under physical constraints (Heesch et al., 2025)."

Sentence S-I2 (L84-94, Introduction): "Symbolic templates (Mirzadeh et al., 2024; Xie et al., 2025)
serve as the generative mechanism, producing unique, contamination-resistant instances".

Sentence S-I3b (L105-111, Introduction): "specialized benchmarks (e.g., Math (Hendrycks et al.,
2021b) or HumanEval (Chen et al., 2021)) assess abstract logic and algorithmic translation."

| key | cited as | sentences |
|---|---|---|
| cobbe2021gsm8k | Cobbe et al., 2021 | S-M1 |
| hendrycks2021math | Hendrycks et al., 2021b | S-M1, S-I3b |
| fan2024hardmath | Fan et al., 2024 | S-M1 |
| gulati2024putnamaxiom | Gulati et al., 2024 | S-M1 (the downloaded text is the 2025 arXiv/ICML version; the paper cites the NeurIPS 2024 workshop version) |
| mirzadeh2024gsmsymbolic | Mirzadeh et al., 2024 | S-M1, S-I2 |
| chen2021humaneval | Chen et al., 2021 | S-M2, S-I3b |
| austin2021mbpp | Austin et al., 2021 | S-M2 |
| jimenez2024swebench | Jimenez et al., 2024 | S-M2 |
| heesch2025realworld | Heesch et al., 2025 | S-M3; note it is itself an engineering evaluation that the engineering paragraph never discusses |

## Batch P: physical sciences (Related Work paragraph 2)

Sentence S-P1 (L172-177): "While UGPhysics (Xu et al., 2025) and PhysReason (Zhang et al., 2025a)
effectively test the application of reasoning to natural laws, others target advanced symbolic
(Zhang et al., 2025b; Qiu et al., 2025) or qualitative (Wang et al., 2023; Shojaee et al., 2025)
understanding."

Sentence S-P2 (L177-180, no citation, characterises all six): "Despite this progress, these
benchmarks remain confined to theoretical physics, lacking the integrative context of engineering."

Rebuttal use R-1 (July 2026 rebuttal, not in the paper): "SciBench (Wang et al., 2023) allows a 5%
relative tolerance for grading numeric answers on college-level scientific problems, and
ABench-Physics (Zhang et al., 2025b) uses a 1% tolerance for physics reasoning."

| key | cited as | sentences |
|---|---|---|
| xu2025ugphysics | Xu et al., 2025 | S-P1 ("effectively test the application of reasoning to natural laws"), S-P2 |
| zhang2025physreason | Zhang et al., 2025a | S-P1 (same), S-P2 |
| zhang2025abenchphysics | Zhang et al., 2025b | S-P1 ("advanced symbolic"), S-P2, R-1 |
| qiu2025phybench | Qiu et al., 2025 | S-P1 ("advanced symbolic"), S-P2 |
| wang2023newton | Wang et al., 2023 | S-P1 ("qualitative"), S-P2 |
| shojaee2025llmsrbench | Shojaee et al., 2025 | S-P1 ("qualitative"), S-P2 |
| wang2023scibench | not in the paper; R-1 | R-1; also judge whether it belongs in the revised paragraph |

## Batch E1: engineering, generalist suites and siloed benchmarks (Related Work paragraph 3, first sentence)

Sentence S-E1 (L182-190): "Generalist suites like MMLU (Hendrycks et al., 2021a), BIG-Bench (Luo
et al., 2024), and SuperGPQA (Du et al., 2025) restrict assessment to factual recall via
multiple-choice questions, while specialized benchmarks like TransportBench, APBench, and EEE-Bench
(Syed et al., 2024; Chen et al., 2025; Cheng et al., 2024; Li et al., 2025; Liu et al., 2025) are
siloed within narrow sub-disciplines."

Sentence S-I3a (L105-111, Introduction): "general knowledge benchmarks (e.g., MMLU (Hendrycks et
al., 2020)) test broad factual recall".

| key | cited as | sentences and notes |
|---|---|---|
| hendrycks2021mmlu | Hendrycks et al., 2021a and 2020 (the same paper listed twice) | S-E1, S-I3a |
| luo2024bigbench_t2i | Luo et al., 2024 | S-E1, as the source for "BIG-Bench"; it is a text-to-image social-bias benchmark also named BIGbench |
| srivastava2022bigbench | not in the reference list | the BIG-bench the sentence names; judge whether S-E1's characterisation fits it |
| du2025supergpqa | Du et al., 2025 | S-E1 |
| syed2024transportbench | Syed et al., 2024 | S-E1 (TransportBench) |
| chen2025apbench | Chen et al., 2025 | S-E1 (APBench); the first author is Di Wu, not Chen |
| cheng2024elecbench | Cheng et al., 2024 | S-E1; ElecBench is cited but never named; the first author is Xiyuan Zhou |
| li2025eeebench | Li et al., 2025 | S-E1 (EEE-Bench) |
| liu2025circuit | Liu et al., 2025 | S-E1; the printed entry merges two papers; this file is the printed title (Chen et al., FIE 2025, homework assessment in circuit analysis, arXiv 2511.18221) |
| skelic2025circuit | not in the reference list | the CIRCUIT benchmark (arXiv 2502.07980) the printed entry's note names; judge which of the two the sentence should cite |

## Batch E2: engineering, broader efforts and the closest precedents (Related Work paragraph 3, second sentence, and the Introduction)

Sentence S-E2 (L190-198): "Others focus on tangential skills such as design generation (Guo et al.,
2025b) or software proficiency (Mudur et al., 2025), while broader efforts like EngiBench (Zhou et
al., 2025) rely on subjective rubric-based scoring."

Sentence S-I3c (L116-122, Introduction): "EngTrace addresses both gaps through integrative reasoning
templates that, unlike existing benchmarks (e.g., Zhou et al., 2025; Mudur et al., 2025) that rely
solely on outcome matching, demand holistic procedural reasoning".

Sentence S-E3 (L199-205, closing): "Evidently, existing benchmarks remain limited to outcome
matching and fail to validate the physically grounded reasoning required for engineering. EngTrace
addresses this limitation by introducing a symbolic benchmark that pairs unique problem instances
with gold-standard reasoning traces to enable verifiable process supervision."

Sentence S-I2 (L84-94, Introduction): "Symbolic templates (Mirzadeh et al., 2024; Xie et al., 2025)
serve as the generative mechanism, producing unique, contamination-resistant instances".

| key | cited as | sentences and notes |
|---|---|---|
| guo2025engdesign | Guo et al., 2025b | S-E2 ("design generation"), S-E3 |
| mudur2025feabench | Mudur et al., 2025 | S-E2 ("software proficiency"), S-I3c ("rely solely on outcome matching"), S-E3 |
| zhou2025engibench | Zhou et al., 2025 | S-E2 ("subjective rubric-based scoring"), S-I3c ("rely solely on outcome matching"), S-E3 |
| felten2025engibench | not in the reference list; a note under Zhou et al. says "Distinct from the design-focused EngiBench by Felten et al." | judge whether it needs citing to disambiguate |
| xie2025finchain | Xie et al., 2025 | S-I2 only; the citation map calls it the closest precedent (a symbolic benchmark with verifiable chain-of-thought, four shared authors). Judge what EngTrace must say it adds beyond the change of domain |

## Batch G: the Introduction's framing, saturation and contamination

Sentence S-I1 (L45-48): "As LLMs expand into high-stakes engineering workflows, rigorous evaluation
of their reasoning capabilities has become paramount (Zhao et al., 2023; Xie et al., 2023a)."

Sentence S-I3 (L97-105): "First, static benchmarks like Glue (Wang et al., 2019) and BBH (Kazemi et
al., 2025) are increasingly prone to model saturation (Ott et al., 2022) and data contamination
(Deng et al., 2024); EngTrace counters this through symbolic templates and domain-aware
parameterization, generating unique, physically grounded problems that resist rote memorization."

Sentence S-T (L466-468, section 4.3): "We evaluate the fluency and semantic fidelity of generated
solutions using standard reference-based metrics (Xie et al., 2023b)."

Use S-K (section 3.3 / Appendix K and the rebuttals): Plank (2022) is cited for treating residual
one-point disagreement between expert annotators as human label variation rather than noise.

| key | cited as | sentences and notes |
|---|---|---|
| zhao2023survey | Zhao et al., 2023 | S-I1 |
| xie2023storytelling | Xie et al., 2023a | S-I1 (a storytelling study, by an EngTrace author) |
| xie2023deltascore | Xie et al., 2023b | S-T (by an EngTrace author) |
| wang2019superglue | Wang et al., 2019 | S-I3, as the source for "Glue" |
| wang2018glue | not in the reference list | the GLUE the sentence names |
| kazemi2025bbeh | Kazemi et al., 2025 | S-I3, as the source for "BBH"; it is BIG-Bench Extra Hard |
| suzgun2022bbh | not in the reference list | the BBH the sentence names |
| ott2022saturation | Ott et al., 2022 | S-I3 (saturation) |
| deng2024contamination | Deng et al., 2024 | S-I3 (contamination) |
| plank2022hlv | Plank, 2022 | S-K |

## Batch V: the evaluation framework's method citations

Sentence S-V1 (L390-395, section 4): "Inspired by recent advances in process supervision (Lightman
et al., 2023) and cascaded model verification (Chen et al., 2023), we introduce a two-stage
evaluation framework that validates intermediate reasoning steps in addition to final-answer
accuracy through a tiered verification protocol."

Sentence S-V2 (L260-263, section 3.3): "Building upon the 'Panel of LLM' (PoLL) methodology (Verga
et al., 2024), we construct an 'AI Tribunal' comprising three frontier reasoning models".

Sentence S-V3 (L469-471, section 4.3): "We employ BERTScore (Zhang et al., 2020) to measure deep
semantic similarity via contextual embedding alignment, and ROUGE-2 (Lin, 2004) to quantify lexical
overlap."

For this batch also judge what each paper says that bears on the planned new paragraph on process
supervision and LLM-judge validity (the state brief, section 6), and what the revised text may cite
it for given that the AI Tribunal is now only a template screen and the evaluator is deterministic
with a residual judge.

| key | cited as | sentences |
|---|---|---|
| lightman2023verify | Lightman et al., 2023 | S-V1 (the reference list prints 30 authors; the paper has 10) |
| chen2023frugalgpt | Chen et al., 2023 | S-V1 ("cascaded model verification") |
| verga2024poll | Verga et al., 2024 | S-V2 |
| zhang2020bertscore | Zhang et al., 2020 | S-V3 |
| lin2004rouge | Lin, 2004 | S-V3 |
