# Citation-context map: EngTrace, ARR May 2026 submission

Source: `docs/_ARR_May__EngTrace.txt`, the 3,276-line text dump of `docs/_ARR_May__EngTrace.pdf`. Main text is L1–861, the references L862–1250 and the appendices L1254–3276. Every `L` number below is a line of that dump. Quotes restore the spaces the dump drops; they were checked against the PDF's own text layer, and paragraph breaks against the PDF's first-line indents. Prepared 2026-10-01 for the related-work panel. Nothing in the paper was changed.

## Summary

- **79 reference entries**, not about 95. Counting lines that look like entry starts gives 100, because continuation lines such as "Released November 24, 2025." also begin with a capital after a full stop or a page break; 21 of the 100 are such lines.
- **ORPHAN: none.** Every entry is cited at least once.
- **DANGLING: none.** Every author–year citation in the text resolves to an entry. Works the paper names without citing (the Appendix C and E sources, "EngiBench by Felten et al.", two models) are listed after Part A.
- **90 citation instances**: 81 in the main text, 9 in the appendices. Sections 1 and 2 hold 42 instances of 37 distinct entries.
- **Types**: 30 benchmark, 22 model or system card, 10 metric or statistics, 9 textbook or standard, 8 survey or other.
- **Worst mismatches (Part C)**: three benchmark names point at other works (Glue → SuperGLUE, BBH → BIG-Bench Extra Hard, BIG-Bench → a text-to-image social-bias benchmark); five keys hang on three names; MMLU is listed twice; EngiBench is described in two incompatible ways; the Liu et al. 2025 entry contradicts itself.
- **Supplement (C-S)**: the printed author lists were compared with the arXiv record already in `notes/arxiv_abs_check.json`. 13 of the 38 cited papers in that record have wrong author data, and for two of them the in-text key names the wrong first author.

### Legend for "Cited where"

| Label | Location |
|---|---|
| §1 ¶1, ¶2, ¶3 | Introduction, paragraphs 1–3 (L45–52, L77–94, L95–129) |
| §2 Math | Related Work, "Benchmarks in Mathematics and Coding" (L153–168) |
| §2 Phys | Related Work, "Benchmarks in the Physical Sciences" (L169–180) |
| §2 Eng | Related Work, "Benchmarks in Engineering" (L181–198) |
| §3.1 | Taxonomy and Content Selection |
| §3.3 | Template Validation |
| §4 | Evaluation Framework, opening paragraph |
| §4.3 | Textual Quality |
| §5.1 | Evaluation Model Suite |
| App. K.2 | Template Validation Results, Stage 2: Inter-Annotator Agreement |
| App. K.3 fn1 | AI-Human Alignment, footnote 1 (marker at L2338, text at L2352–2356) |
| App. T.1 | Inter-Annotator Agreement, Computed Metrics |

No citation occurs in the abstract, the contributions list (L130–151), the Related Work closing paragraph (L199–205), §3.2, §3.4, §4.1, §4.2, §4.4, §5.2–§6, Limitations, the ethics statement, or Appendices A–J, K.1, L–S, T.2 and U and the appendix tables. Appendix G's template code was scanned for author–year patterns and has none.

## Part A: Reference list, in printed order

The key is the in-text label (first-author surname, or both surnames for two-author works, plus year with the paper's a/b suffix), followed by the line where the entry starts. Titles and venues are as printed, with the dump's dropped spaces restored and URLs rejoined across line breaks; capitalisation is left as printed.

| # | Key (entry line) | Title as printed | Venue / id as printed | Type | Cited where | Part C |
|---|---|---|---|---|---|---|
| 1 | ABET 2024 (L863) | Criteria for accrediting engineering programs | Web page, `https://www.abet.org/accreditation/accreditation-criteria/criteria-for-accrediting-engineering-programs-2024-2025/` | textbook or standard (accreditation criteria) | §3.1 (L213). "ABET" also appears in the App. A, B and D prompt texts (L1266, L1322, L1417), not as a citation | — |
| 2 | AIChE 2025 (L867) | Constitution & bylaws | Web page, `https://www.aiche.org/about/governance/constitution-bylaws` | textbook or standard (society web page) | §3.1 (L215) | C10 |
| 3 | Anthropic 2025a (L870) | Claude 3.7 sonnet system card | `https://assets.anthropic.com/m/785e231869ea8b3b/original/claude-3-7-sonnet-system-card.pdf` | model or system card | §5.1 (L560) | C15h |
| 4 | Anthropic 2025b (L874) | Claude Opus 4.5 System Card | `https://assets.anthropic.com/m/64823ba7485345a7/Claude-Opus-4-5-System-Card.pdf`. Released November 24, 2025. | model or system card | §3.3 (L264–265; AI Tribunal judge) | — |
| 5 | Anthropic 2025c (L878) | Claude sonnet 4.5 system card | `https://assets.anthropic.com/m/12f214efcc2f457a/original/Claude-Sonnet-4-5-System-Card.pdf` | model or system card | §5.1 (L560) | C15h |
| 6 | Anthropic 2025d (L882) | System card: Claude opus 4 & claude sonnet 4 | `https://www-cdn.anthropic.com/6be99a52cb68eb70eb9572b4cafad13df32ed995.pdf` | model or system card | §5.1 (L560) | C15h |
| 7 | Anthropic 2026 (L885) | Claude Opus 4.7 model card | Technical report, Anthropic. Model identifier: claude-opus-4-7. Available via Anthropic API, Amazon Bedrock, Google Cloud Vertex AI, and Microsoft Foundry. (no URL) | model or system card | §5.1 (L560–561) | C15e |
| 8 | Artstein and Poesio 2008 (L890) | Inter-coder agreement for computational linguistics | Computational Linguistics, 34(4):555–596 | metric or statistics | App. K.3 fn1 (L2354) | C14b |
| 9 | ASME 2025 (L893) | Vision, mission, & core values | Web page, `https://www.asme.org/about-asme/vision-mission-core-values` | textbook or standard (society web page) | §3.1 (L214–215) | C10 |
| 10 | Austin 2021 (L896) | Program synthesis with large language models | arXiv preprint arXiv:2108.07732 | benchmark (MBPP) | §2 Math (L163–164) | — |
| 11 | Byrt 1993 (L901) | Bias, prevalence and kappa | Journal of Clinical Epidemiology, 46(5):423–429 | metric or statistics (PABAK) | App. T.1 (L2726) | — |
| 12 | Chen 2023 (L904) | FrugalGPT: How to use large language models while reducing cost and improving performance | In Proceedings of the 40th International Conference on Machine Learning (ICML) | survey or other (LLM cascade method) | §4 (L392) | C14a, C15f |
| 13 | Chen 2021 (L909) | Evaluating large language models trained on code | arXiv preprint arXiv:2107.03374 | benchmark (HumanEval) | §1 ¶3 (L110); §2 Math (L163) | — |
| 14 | Chen 2025 (L915) | Apbench and benchmarking large language model performance in fundamental astrodynamics problems for space engineering | Scientific Reports, 15 | benchmark (APBench) | §2 Eng (L188) | C4, C15f, C15h |
| 15 | Cheng 2024 (L921) | Elecbench: A power dispatch evaluation benchmark for large language models | ArXiv preprint, abs/2407.05365 | benchmark (ElecBench) | §2 Eng (L188–189) | C4, C15h, C-S |
| 16 | Cobbe 2021 (L930) | Training verifiers to solve math word problems | arXiv preprint arXiv:2110.14168 | benchmark (GSM8K) | §2 Math (L156) | — |
| 17 | Comanici 2025 (L933) | Gemini 2.5: Pushing the frontier with advanced reasoning, multimodality, long context, and next generation agentic capabilities | arXiv preprint arXiv:2507.06261 | model or system card | §5.1 (L563) | — |
| 18 | DeepSeek-AI 2026 (L940) | DeepSeek-V4: Towards highly efficient million-token context intelligence | Technical report, DeepSeek. Model identifier: deepseek-ai/DeepSeek-V4-Pro. 1.6T-parameter MoE model released under the MIT License. (no URL) | model or system card | §5.1 (L573–574) | C15e |
| 19 | Deng 2024 (L945) | Investigating data contamination in modern benchmarks for large language models | In Proceedings of the 2024 Conference of the North American Chapter of the Association for Computational Linguistics: Human Language Technologies (Volume 1: Long Papers), pages 8706–8719 | survey or other (contamination study) | §1 ¶3 (L101) | — |
| 20 | Du 2025 (L952) | Supergpqa: Scaling llm evaluation across 285 graduate disciplines | ArXiv preprint, abs/2502.14739 | benchmark (SuperGPQA) | §2 Eng (L184–185) | C15h |
| 21 | Dym 2012 (L956) | Engineering Design: A Project-Based Introduction, 4th edition | John Wiley & Sons | textbook or standard (textbook) | §1 ¶3 (L115) | — |
| 22 | Fan 2024 (L959) | HARDMath: A benchmark dataset for challenging problems in applied mathematics | arXiv preprint arXiv:2410.09988 | benchmark | §2 Math (L158) | — |
| 23 | Feinstein and Cicchetti 1990 (L964) | High agreement but low kappa: I. The problems of two paradoxes | Journal of Clinical Epidemiology, 43(6):543–549 | metric or statistics | §3.3 (L368); App. K.2 (L2238) | C14c |
| 24 | Fogler 2020 (L968) | Elements of Chemical Reaction Engineering, 6th edition | Pearson, Boston, MA | textbook or standard (textbook) | §3.1 (L224–225). Also named, without a key, in App. C Table 2 (L1367–1368) | — |
| 25 | Gemma Team and Google DeepMind 2025 (L970) | Gemma 3: Technical report | arXiv preprint, abs/2503.19786 | model or system card | §5.1 (L572–573) | — |
| 26 | Gemma Team 2024 (L972) | Gemma 2: Improving open language models at a practical size | arXiv preprint, abs/2408.00118 | model or system card | §5.1 (L571–572) | — |
| 27 | Google 2025 (L977) | Gemini 3: A new era of intelligence | `https://blog.google/products/gemini/gemini-3/`. Released November 18, 2025. Model documentation available at `https://ai.google.dev/gemini-api/docs/gemini-3` | model or system card (blog post) | §3.3 (L265; AI Tribunal judge); §5.1 (L562–563) | — |
| 28 | Google DeepMind 2026 (L982) | Gemini 3.1 Pro: A more capable reasoning model in the Gemini 3 series | Technical report, Google DeepMind. Model identifier: gemini-3.1-pro-preview. Available via Google AI Studio, Vertex AI, and Gemini API. (no URL) | model or system card | §5.1 (L562) | C15e |
| 29 | Grattafiori 2024 (L987) | The llama 3 herd of models | ArXiv preprint, abs/2407.21783 | model or system card | §5.1 (L569) | C15h |
| 30 | Gulati 2024 (L992) | Putnam-axiom: A functional and static benchmark for measuring higher level mathematical reasoning | In The 4th Workshop on Mathematical Reasoning and AI at NeurIPS'24 | benchmark | §2 Math (L158–159) | C15h, C-S |
| 31 | Guo 2025a (L998) | Deepseek-r1: Incentivizing reasoning capability in llms via reinforcement learning | arXiv preprint arXiv:2501.12948 | model or system card | §5.1 (L564) | C15h |
| 32 | Guo 2025b (L1003) | Toward engineering agi: Benchmarking the engineering design capabilities of llms | ArXiv preprint, abs/2509.16204 | benchmark (engineering design) | §2 Eng (L191–192) | C15h |
| 33 | Gwet 2008 (L1007) | Computing inter-rater reliability and its variance in the presence of high agreement | British Journal of Mathematical and Statistical Psychology, 61(1):29–48 | metric or statistics | §3.3 (L365); App. K.2 (L2232); App. T.1 (L2727) | C11 |
| 34 | Heesch 2025 (L1011) | Evaluating large language models for real-world engineering tasks | arXiv preprint arXiv:2505.13484 | benchmark (engineering evaluation) | §2 Math (L168) | C13 |
| 35 | Hendrycks 2020 (L1016) | Measuring massive multitask language understanding | arXiv preprint arXiv:2009.03300 | benchmark (MMLU) | §1 ¶3 (L107) | C5 |
| 36 | Hendrycks 2021a (L1020) | Measuring massive multitask language understanding | In International Conference on Learning Representations | benchmark (MMLU) | §2 Eng (L183) | C5 |
| 37 | Hendrycks 2021b (L1025) | Measuring mathematical problem solving with the MATH dataset | arXiv preprint arXiv:2103.03874 | benchmark (MATH) | §1 ¶3 (L109); §2 Math (L157) | C5 |
| 38 | IEEE 2025 (L1030) | IEEE standards association | Web page, `https://standards.ieee.org/` | textbook or standard (standards-body home page) | §3.1 (L215) | C10 |
| 39 | Jimenez 2024 (L1032) | SWE-bench: Can language models resolve real-world github issues? | In The Twelfth International Conference on Learning Representations | benchmark | §2 Math (L165) | C15h, C-S |
| 40 | Kazemi 2025 (L1041) | BIG-Bench Extra Hard | In Proceedings of the 63rd Annual Meeting of the Association for Computational Linguistics (Volume 1: Long Papers), pages 26473–26501, Vienna, Austria | benchmark (BBEH) | §1 ¶3 (L99) | C2, C-S |
| 41 | Landis and Koch 1977 (L1050) | The measurement of observer agreement for categorical data | Biometrics, 33(1):159–174 | metric or statistics | App. T.1 (L2735, narrative form "Landis and Koch (1977)") | — |
| 42 | Leveson 2012 (L1053) | Engineering a Safer World: Systems Thinking Applied to Safety | MIT Press, Cambridge, MA | textbook or standard (book) | §1 ¶3 (L122) | — |
| 43 | Li 2025 (L1056) | Eee-bench: A comprehensive multimodal electrical and electronics engineering benchmark | ArXiv preprint, abs/2411.01492. Accepted to CVPR 2025. | benchmark | §2 Eng (L189) | C4, C15h |
| 44 | Lightman 2023 (L1061) | Let's verify step by step | ArXiv preprint, abs/2305.20050 | survey or other (process-supervision method) | §4 (L391) | C-S |
| 45 | Lin 2004 (L1072) | ROUGE: A package for automatic evaluation of summaries | In Text Summarization Branches Out, pages 74–81, Barcelona, Spain. Association for Computational Linguistics | metric or statistics (metric) | §4.3 (L471) | — |
| 46 | Liu 2024 (L1076) | Deepseek-v3 technical report | arXiv preprint arXiv:2412.19437 | model or system card | §5.1 (L564) | C15h |
| 47 | Liu 2025 (L1080) | Enhancing large language models for automated homework assessment in undergraduate circuit analysis | ArXiv preprint, abs/2502.07980. Introduces the CIRCUIT benchmark. | benchmark (per its note) | §2 Eng (L189) | C4, C8, C-S |
| 48 | Luo 2023 (L1084) | Wizardmath: Empowering mathematical reasoning for large language models via reinforced evol-instruct | ArXiv preprint, abs/2308.09583 | model or system card | §5.1 (L580) | C15h |
| 49 | Luo 2024 (L1090) | BIGbench: A unified benchmark for evaluating multi-dimensional social biases in text-to-image models | arXiv preprint arXiv:2407.15240 | benchmark (text-to-image social bias) | §2 Eng (L183–184, as "BIG-Bench") | C3, C-S |
| 50 | Mirzadeh 2024 (L1096) | GSM-Symbolic: Understanding the limitations of mathematical reasoning in large language models | arXiv preprint arXiv:2410.05229 | benchmark | §1 ¶2 (L84); §2 Math (L159–160) | C-S |
| 51 | Mistral 2024 (L1101) | Mathstral | 7B parameter model for mathematical reasoning, released under Apache 2.0 license. (no URL; author printed "AI Mistral") | model or system card | §5.1 (L579) | C15c |
| 52 | Moaveni 2019 (L1104) | Engineering Fundamentals: An Introduction to Engineering, 6th edition | Cengage Learning | textbook or standard (textbook) | §1 ¶3 (L115) | — |
| 53 | Mudur 2025 (L1107) | FEABench: Evaluating language models on multiphysics reasoning ability | arXiv preprint arXiv:2504.06260 | benchmark | §1 ¶3 (L118); §2 Eng (L192–193) | C7 |
| 54 | OpenAI 2025a (L1112) | Gpt-5 system card | `https://cdn.openai.com/gpt-5-system-card.pdf`. System card describing GPT-5 variants (thinking/main/mini/nano). | model or system card | §3.3 (L264; AI Tribunal judge); §5.1 (L558–559) | C15h |
| 55 | OpenAI 2025b (L1116) | Introducing gpt-4.1 in the api | `https://openai.com/index/gpt-4-1/`. Product research post introducing GPT-4.1, 4.1 mini, and 4.1 nano. | model or system card (product post) | §5.1 (L558–559) | C15h |
| 56 | Ott 2022 (L1119) | Mapping global dynamics of benchmark creation and saturation in artificial intelligence | Nature Communications, 13(1):6793 | survey or other (benchmark meta-analysis) | §1 ¶3 (L100) | — |
| 57 | Plank 2022 (L1123) | The "problem" of human label variation: On ground truth in data, modeling and evaluation | In Proceedings of the 2022 Conference on Empirical Methods in Natural Language Processing, pages 10671–10682, Abu Dhabi, United Arab Emirates. Association for Computational Linguistics | survey or other (position paper) | App. K.2 (L2283) | — |
| 58 | Qiu 2025 (L1129) | PHYBench: Holistic evaluation of physical perception and reasoning in large language models | arXiv preprint arXiv:2504.16074 | benchmark | §2 Phys (L176) | C12, C-S |
| 59 | Qwen 2024 (L1133) | Qwen2.5 technical report | ArXiv preprint, abs/2412.15115 | model or system card | §5.1 (L570) | C15d |
| 60 | Qwen 2025 (L1135) | Qwen3 | (nothing printed) | model or system card | §5.1 (L570–571) | C15b, C15d |
| 61 | Shojaee 2025 (L1136) | LLM-SRBench: A new benchmark for scientific equation discovery with large language models | arXiv preprint arXiv:2504.10415 | benchmark | §2 Phys (L177) | C12, C-S |
| 62 | Syed 2024 (L1141) | Benchmarking the capabilities of large language models in transportation system engineering: Accuracy, consistency, and reasoning behaviors | ArXiv preprint, abs/2408.08302 | benchmark (cited as TransportBench) | §2 Eng (L188) | C4 |
| 63 | Ulaby and Ravaioli 2015 (L1151) | Fundamentals of Applied Electromagnetics, 7th edition | Pearson, Upper Saddle River, NJ | textbook or standard (textbook) | §3.1 (L225). Also named, without a key, in App. C Table 2 (L1378–1380) | — |
| 64 | Verga 2024 (L1154) | Replacing judges with juries: Evaluating LLM generations with a panel of diverse models | ArXiv preprint, abs/2404.18796 | survey or other (LLM-jury evaluation method) | §3.3 (L262) | — |
| 65 | Wang 2019 (L1160) | SuperGLUE: A stickier benchmark for general-purpose language understanding systems | arXiv preprint arXiv:1905.00537 | benchmark | §1 ¶3 (L98, as "Glue") | C1 |
| 66 | Wang 2023 (L1165) | NEWTON: Are large language models capable of physical reasoning? | arXiv preprint arXiv:2310.07018 | benchmark | §2 Phys (L176–177) | C12, C-S |
| 67 | Wongpakaran 2013 (L1169) | A comparison of Cohen's Kappa and Gwet's AC1 when calculating inter-rater reliability coefficients: a study conducted with personality disorder samples | BMC Medical Research Methodology, 13:61 | metric or statistics | §3.3 (L366); App. K.2 (L2232–2233) | C11 |
| 68 | Xie 2023a (L1175) | The next chapter: A study of large language models in storytelling | In Proceedings of the 16th International Natural Language Generation Conference, pages 323–351, Prague, Czechia. Association for Computational Linguistics | survey or other (storytelling study) | §1 ¶1 (L48) | C9 |
| 69 | Xie 2023b (L1181) | DeltaScore: Fine-grained story evaluation with perturbations | In Findings of the Association for Computational Linguistics: EMNLP 2023, pages 5317–5331, Singapore. Association for Computational Linguistics | metric or statistics (story-evaluation metric) | §4.3 (L468) | C9 |
| 70 | Xie 2025 (L1187) | FINCHAIN: A symbolic benchmark for verifiable chain-of-thought financial reasoning | arXiv preprint. (no identifier) | benchmark | §1 ¶2 (L85) | C15a; D.3 note 1 |
| 71 | Xu 2025 (L1197) | UGPhysics: A comprehensive benchmark for undergraduate physics reasoning with large language models | arXiv preprint arXiv:2502.00334 | benchmark | §2 Phys (L172–173) | C-S |
| 72 | Yang 2024 (L1202) | Qwen2.5-math technical report: Toward mathematical expert model via self-improvement | ArXiv preprint, abs/2409.12122 | model or system card | §5.1 (L578–579) | C15h |
| 73 | Yu 2023 (L1209) | Metamath: Bootstrap your own mathematical questions for large language models | ArXiv preprint, abs/2309.12284 | model or system card | §5.1 (L581) | C15h |
| 74 | Zapf 2016 (L1214) | Measuring inter-rater reliability for nominal data: which coefficients and confidence intervals are appropriate? | BMC Medical Research Methodology, 16(1):93 | metric or statistics | App. T.1 (L2734) | C14e |
| 75 | Zhang 2020 (L1219) | Bertscore: Evaluating text generation with BERT | In 8th International Conference on Learning Representations, ICLR 2020, Addis Ababa, Ethiopia, April 26-30, 2020. OpenReview.net | metric or statistics (metric) | §4.3 (L469) | C15h |
| 76 | Zhang 2025a (L1225) | PhysReason: A comprehensive benchmark towards physics-based reasoning | arXiv preprint arXiv:2502.12054 | benchmark | §2 Phys (L173) | C-S |
| 77 | Zhang 2025b (L1230) | ABench-Physics: Benchmarking physical reasoning in LLMs via high-difficulty and dynamic physics problems | arXiv preprint arXiv:2507.04766 | benchmark | §2 Phys (L175–176) | C12, C-S |
| 78 | Zhao 2023 (L1236) | A survey of large language models | ArXiv preprint, abs/2303.18223 | survey or other (survey) | §1 ¶1 (L47) | C14d |
| 79 | Zhou 2025 (L1244) | Engibench: A benchmark for evaluating large language models on engineering problem solving | ArXiv preprint, abs/2509.17677. Distinct from the design-focused EngiBench by Felten et al. | benchmark | §1 ¶3 (L118); §2 Eng (L193–198, split across the page break) | C6, C15g, C15h |

ORPHAN entries: none. DANGLING citations: none.

### Named in the paper without a citation or an entry

These are not counted as DANGLING, which the brief defines as author–year citations without an entry, but a reviewer may still ask for references.

- **Appendix C, Table 2 (L1366–1401).** §3.1 calls this table "the complete list of source texts" (L225–226). Seven of its nine textbooks have no entry: Smith, Van Ness and Abbott, *Introduction to Chemical Engineering Thermodynamics* (L1369–1372); Bird, Stewart and Lightfoot, *Transport Phenomena* (L1373–1375); Proakis and Salehi, *Digital Communications* (L1381–1384); Oppenheim and Schafer, *Discrete-Time Signal Processing* (L1385–1387); Beer, Johnston and DeWolf, *Mechanics of Materials* (L1388–1392); Munson, Young and Okiishi, *Fundamentals of Fluid Mechanics* (L1393–1395); Rao, *Mechanical Vibrations* (L1396–1399). Fogler and Ulaby–Ravaioli, in the same table, do have entries.
- **Appendix E, Table 3 (L1427–1451).** None of the nine data sources has an entry: Perry's Chemical Engineers' Handbook; NIST Chemistry WebBook; CRC Handbook of Chemistry and Physics; IEEE 100: The Authoritative Dictionary of IEEE Standards Terms; Standard Handbook for Electrical Engineers; ART-DEIT Database; ASM Handbook, Volume 1; Marks' Standard Handbook for Mechanical Engineers; Shigley's Mechanical Engineering Design.
- **"EngiBench by Felten et al."**, named inside the Zhou et al. 2025 entry (L1248–1249). The arXiv record identifies it as arXiv:2508.00831, "EngiBench: A Framework for Data-Driven Engineering Design Research".
- **GPT-4o and Claude 3.5 Sonnet**, on the Appendix A model panel (L1256). Gemini 2.5 Pro and DeepSeek V3, on the same panel, have entries but are not cited there.
- **Methods and tools named without a citation**: Fleiss' κ (first at L363), Krippendorff's α (L365), Cohen's κ (L2724), Spearman's ρ (L2336, L2392, L2411), the Hungarian matching algorithm (L543), the Cross-Encoder used in Tier 1 (L495; the model is not named), Chain-of-Thought (Figure 1, L64) and Streamlit (L361, L2133).
- **The abstract** names MMLU, MATH and HumanEval without citations (L18–19), as abstracts usually do. §1 cites all three.

## Part B: Every citation in §1 and §2, by paragraph

Each sentence that carries a citation is quoted once, with its lines, followed by one line per citation stating the claim it is used to support. A pointer to Part C marks where the entry does not obviously support that claim. Sentences without citations are omitted, except the Related Work closing paragraph, which the brief asks for. The Figure 1 caption and the contributions list (L130–151) carry no citations.

### B.1 Introduction, paragraph 1 (L45–52)

> "As LLMs expand into high-stakes engineering workflows, rigorous evaluation of their reasoning capabilities has become paramount (Zhao et al., 2023; Xie et al., 2023a)." (L45–48)

- **Zhao et al., 2023** is cited as support that rigorous evaluation of LLM reasoning has become paramount as LLMs move into high-stakes engineering workflows. The entry is a general survey of LLMs (C14d).
- **Xie et al., 2023a** is cited as support for the same claim. The entry is a study of LLMs in storytelling (C9).

### B.2 Introduction, paragraph 2 (L77–94)

> "Symbolic templates (Mirzadeh et al., 2024; Xie et al., 2025) serve as the generative mechanism, producing unique, contamination-resistant instances; Figure 1 illustrates one such instance, a CSTR volume calculation requiring synthesis of interdependent stoichiometry and reaction kinetics variables into a fully verifiable reasoning trace." (L84–94)

- **Mirzadeh et al., 2024** (GSM-Symbolic) is cited as precedent for symbolic templates as the mechanism that generates unique, contamination-resistant problem instances.
- **Xie et al., 2025** (FinChain) is cited as a second precedent for symbolic templates. Its title, "A symbolic benchmark for verifiable chain-of-thought financial reasoning", also claims verifiable reasoning, which the paragraph's first sentence (L77–80) says no existing benchmark offers for engineering; see D.3 note 1.

### B.3 Introduction, paragraph 3 (L95–129)

> "First, static benchmarks like Glue (Wang et al., 2019) and BBH (Kazemi et al., 2025) are increasingly prone to model saturation (Ott et al., 2022) and data contamination (Deng et al., 2024); EngTrace counters this through symbolic templates and domain-aware parameterization, generating unique, physically grounded problems that resist rote memorization." (L97–105)

- **Wang et al., 2019** is cited as the source for "Glue", an example of a static benchmark prone to saturation and contamination. The entry is SuperGLUE (C1).
- **Kazemi et al., 2025** is cited as the source for "BBH", a second example of a saturating static benchmark. The entry is BIG-Bench Extra Hard (C2).
- **Ott et al., 2022** is cited as evidence that static benchmarks become saturated.
- **Deng et al., 2024** is cited as evidence that static benchmarks suffer data contamination.

> "Second, most benchmarks evaluate skills in disciplinary silos: general knowledge benchmarks (e.g., MMLU (Hendrycks et al., 2020)) test broad factual recall, while specialized benchmarks (e.g., Math (Hendrycks et al., 2021b) or HumanEval (Chen et al., 2021)) assess abstract logic and algorithmic translation." (L105–111)

- **Hendrycks et al., 2020** is cited as the source for MMLU, the example of a general-knowledge benchmark that tests broad factual recall. It is the same paper as Hendrycks et al., 2021a (C5).
- **Hendrycks et al., 2021b** is cited as the source for MATH, an example of a specialized benchmark that assesses abstract logic.
- **Chen et al., 2021** is cited as the source for HumanEval, an example of a specialized benchmark that assesses algorithmic translation.

> "This fragmentation is ill-suited to engineering, which is fundamentally integrative, requiring the synthesis of scientific principles, mathematical modeling, and practical constraints (e.g., Moaveni, 2019; Dym et al., 2012)." (L111–115)

- **Moaveni, 2019** (an introductory engineering textbook) is cited for the claim that engineering is fundamentally integrative, combining scientific principles, mathematical modelling and practical constraints.
- **Dym et al., 2012** (an engineering-design textbook) is cited for the same claim.

> "EngTrace addresses both gaps through integrative reasoning templates that, unlike existing benchmarks (e.g., Zhou et al., 2025; Mudur et al., 2025) that rely solely on outcome matching, demand holistic procedural reasoning, a crucial requirement in engineering, where a flawed process can lead to catastrophic failure (Leveson, 2012)." (L116–122)

- **Zhou et al., 2025** (EngiBench) is cited as an example of an existing benchmark that relies solely on outcome matching. §2 describes the same work as relying on subjective rubric-based scoring (C6).
- **Mudur et al., 2025** (FEABench) is cited as a second example of a benchmark that relies solely on outcome matching. §2 describes it as targeting software proficiency (C7).
- **Leveson, 2012** is cited for the claim that in engineering a flawed process can lead to catastrophic failure.

### B.4 Related Work: Benchmarks in Mathematics and Coding (L153–168)

> "In mathematics, benchmarks ranging from GSM8K (Cobbe et al., 2021) and Math (Hendrycks et al., 2021b) to HardMath (Fan et al., 2024), Putnam-Axiom (Gulati et al., 2024), and GSM-Symbolic (Mirzadeh et al., 2024) rigorously test abstract logical deduction and reasoning stability." (L155–161)

- **Cobbe et al., 2021** is cited as the source for GSM8K, a math benchmark that tests abstract logical deduction. The name is not in the entry's title, but this is the GSM8K paper.
- **Hendrycks et al., 2021b** is cited as the source for MATH, in the same role.
- **Fan et al., 2024** is cited as the source for HARDMath, in the same role (applied mathematics, per its title).
- **Gulati et al., 2024** is cited as the source for Putnam-AXIOM, in the same role. Its title ("a functional and static benchmark") also fits the "reasoning stability" half of the claim.
- **Mirzadeh et al., 2024** is cited as the source for GSM-Symbolic, a math benchmark that tests reasoning stability.

> "Similarly, coding evaluations span from function-level generation in HumanEval (Chen et al., 2021) and MBPP (Austin et al., 2021) to real-world software engineering in SWE-Bench (Jimenez et al., 2024)." (L161–165)

- **Chen et al., 2021** is cited as the source for HumanEval, a function-level code-generation benchmark.
- **Austin et al., 2021** is cited as the source for MBPP, a function-level code-generation benchmark. The name is not in the entry's title, but this is the MBPP paper.
- **Jimenez et al., 2024** is cited as the source for SWE-bench, an evaluation on real-world software engineering.

> "Yet neither abstract logical deduction nor algorithmic proficiency transfers to reliable reasoning under physical constraints (Heesch et al., 2025)." (L165–168)

- **Heesch et al., 2025** is cited as evidence that mathematical and coding proficiency do not transfer to reliable reasoning under physical constraints. By its title it is itself an evaluation of LLMs on real-world engineering tasks, which "Benchmarks in Engineering" does not discuss (C13).

### B.5 Related Work: Benchmarks in the Physical Sciences (L169–180)

> "While UGPhysics (Xu et al., 2025) and PhysReason (Zhang et al., 2025a) effectively test the application of reasoning to natural laws, others target advanced symbolic (Zhang et al., 2025b; Qiu et al., 2025) or qualitative (Wang et al., 2023; Shojaee et al., 2025) understanding." (L172–177)

- **Xu et al., 2025** is cited as the source for UGPhysics, a physics benchmark that tests applying reasoning to natural laws.
- **Zhang et al., 2025a** is cited as the source for PhysReason, in the same role.
- **Zhang et al., 2025b** (ABench-Physics, not named in the text) is cited as a benchmark that targets advanced symbolic understanding.
- **Qiu et al., 2025** (PHYBench, not named in the text) is cited as a benchmark that targets advanced symbolic understanding.
- **Wang et al., 2023** (NEWTON, not named in the text) is cited as a benchmark that targets qualitative physical understanding.
- **Shojaee et al., 2025** (LLM-SRBench, not named in the text) is cited as a benchmark that targets qualitative understanding; its title describes scientific equation discovery (C12).

The next sentence carries no citation but characterises all six: "Despite this progress, these benchmarks remain confined to theoretical physics, lacking the integrative context of engineering." (L177–180)

### B.6 Related Work: Benchmarks in Engineering (L181–198)

> "Generalist suites like MMLU (Hendrycks et al., 2021a), BIG-Bench (Luo et al., 2024), and SuperGPQA (Du et al., 2025) restrict assessment to factual recall via multiple-choice questions, while specialized benchmarks like TransportBench, APBench, and EEE-Bench (Syed et al., 2024; Chen et al., 2025; Cheng et al., 2024; Li et al., 2025; Liu et al., 2025) are siloed within narrow sub-disciplines." (L182–190)

- **Hendrycks et al., 2021a** is cited as the source for MMLU, a generalist suite restricted to multiple-choice factual recall. Same paper as Hendrycks et al., 2020 (C5).
- **Luo et al., 2024** is cited as the source for BIG-Bench, a generalist multiple-choice suite. The entry is a text-to-image social-bias benchmark named "BIGbench" (C3).
- **Du et al., 2025** is cited as the source for SuperGPQA, a generalist multiple-choice suite.
- **Syed et al., 2024** is cited as the source for TransportBench, a specialized benchmark siloed in transportation engineering.
- **Chen et al., 2025** is cited as the source for APBench, siloed in astrodynamics for space engineering.
- **Cheng et al., 2024** is cited within the same group as a siloed specialized benchmark, but its benchmark (ElecBench, power dispatch) is not named in the sentence (C4).
- **Li et al., 2025** is cited as the source for EEE-Bench, siloed in electrical and electronics engineering.
- **Liu et al., 2025** is cited within the same group, but its work (circuit-analysis homework assessment, or the CIRCUIT benchmark per its note) is not named in the sentence (C4, C8).

> "Others focus on tangential skills such as design generation (Guo et al., 2025b) or software proficiency (Mudur et al., 2025), while broader efforts like EngiBench (Zhou et al., 2025) rely on subjective rubric-based scoring." (L190–198)

- **Guo et al., 2025b** is cited as an engineering benchmark that targets a tangential skill, design generation.
- **Mudur et al., 2025** is cited as an engineering benchmark that targets a tangential skill, software proficiency. Its title says "multiphysics reasoning ability" (C7).
- **Zhou et al., 2025** is cited as the source for EngiBench, a broader engineering benchmark that relies on subjective rubric-based scoring. The Introduction says the same work relies solely on outcome matching (C6).

### B.7 Related Work: closing paragraph (L199–205), no citations

> "Evidently, existing benchmarks remain limited to outcome matching and fail to validate the physically grounded reasoning required for engineering. EngTrace addresses this limitation by introducing a symbolic benchmark that pairs unique problem instances with gold-standard reasoning traces to enable verifiable process supervision." (L199–205)

The generalisation rests on the citations above. The sentence just before it (L193–198) describes EngiBench as rubric-scored, which is not outcome matching (C6). The PDF indents L199, so this is a separate paragraph, not the tail of "Benchmarks in Engineering".

## Part C: Mismatches visible from the reference list

No in-text key disagrees with its own entry, because the keys are generated from the entries. The problems lie between a name or claim in the text and the entry it points to, between two entries, or inside one entry. Severity: **High** means the citation points to the wrong work or contradicts the paper; **Medium** means the cited work does not obviously support the claim; **Low** means a judgement call or bibliographic hygiene. C1–C15 use only the paper and its reference list. C-S goes one step further and uses the arXiv record already in this folder.

**C1 (High). "Glue" is attributed to SuperGLUE.**
- Text, L98: "static benchmarks like Glue (Wang et al., 2019)".
- Entry, L1160–1164: Wang et al. 2019, "SuperGLUE: A stickier benchmark for general-purpose language understanding systems", arXiv:1905.00537.
- GLUE and SuperGLUE are different benchmarks with different papers. The repo's arXiv record lists the GLUE paper as `wang2018glue` (arXiv:1804.07461, ICLR 2019).

**C2 (High). "BBH" is attributed to BIG-Bench Extra Hard.**
- Text, L99: "and BBH (Kazemi et al., 2025) are increasingly prone to model saturation".
- Entry, L1041–1049: Kazemi et al. 2025, "BIG-Bench Extra Hard", ACL 2025.
- BBH is BIG-Bench Hard. BBEH is a later, separate benchmark whose title marks it as the harder successor, which makes it an odd example of a saturating benchmark. The arXiv record lists BBH as `suzgun2022bbh` (arXiv:2210.09261).

**C3 (High). "BIG-Bench" is attributed to a text-to-image social-bias benchmark.**
- Text, L183–184: "Generalist suites like MMLU (Hendrycks et al., 2021a), BIG-Bench (Luo et al., 2024), and SuperGPQA (Du et al., 2025) restrict assessment to factual recall via multiple-choice questions".
- Entry, L1090–1095: Luo et al. 2024, "BIGbench: A unified benchmark for evaluating multi-dimensional social biases in text-to-image models", arXiv:2407.15240.
- A name collision: the entry is neither the generalist BIG-bench suite nor a multiple-choice factual-recall test. The arXiv record lists BIG-bench as `srivastava2022bigbench` (arXiv:2206.04615, TMLR).

**C4 (High). Five keys hang on three benchmark names.**
- Text, L187–189: "specialized benchmarks like TransportBench, APBench, and EEE-Bench (Syed et al., 2024; Chen et al., 2025; Cheng et al., 2024; Li et al., 2025; Liu et al., 2025)".
- Entries: Syed 2024 (L1141–1150, transportation system engineering; the name TransportBench is not in its title); Chen 2025 (L915–920, APBench); Cheng 2024 (L921–925, "Elecbench: A power dispatch evaluation benchmark"); Li 2025 (L1056–1060, EEE-Bench); Liu 2025 (L1080–1083, circuit-analysis homework assessment, "Introduces the CIRCUIT benchmark").
- ElecBench and CIRCUIT are cited but never named, so a reader cannot tell which key goes with which name.

**C5 (Medium). MMLU is listed twice, under two keys.**
- Entries L1016–1019 (Hendrycks et al. 2020, arXiv:2009.03300) and L1020–1024 (Hendrycks et al. 2021a, ICLR) have the same title and the same seven authors.
- The text cites "MMLU (Hendrycks et al., 2020)" at L107 (§1 ¶3) and "MMLU (Hendrycks et al., 2021a)" at L183 (§2 Eng), so one benchmark appears with two years. Merging the entries also turns "2021b" (MATH, L109, L157) into "2021".

**C6 (High). EngiBench is described in two incompatible ways.**
- L116–119: "unlike existing benchmarks (e.g., Zhou et al., 2025; Mudur et al., 2025) that rely solely on outcome matching".
- L193–198: "broader efforts like EngiBench (Zhou et al., 2025) rely on subjective rubric-based scoring".
- Rubric-based scoring is not outcome matching. The entry's own note (L1248–1249) also mentions a second, same-named "EngiBench by Felten et al." that has no entry.

**C7 (Medium). FEABench is called "software proficiency", but its title says multiphysics reasoning.**
- L190–193: "Others focus on tangential skills such as design generation (Guo et al., 2025b) or software proficiency (Mudur et al., 2025)".
- Entry, L1107–1111: "FEABench: Evaluating language models on multiphysics reasoning ability".
- By its title, FEABench targets physical reasoning, the skill the paper claims as its own. It is also one of the two "outcome matching" examples at L118.

**C8 (High). The Liu et al. 2025 entry contradicts itself.**
- Entry, L1080–1083: the author is given only as "Jiaheng Liu et al."; the title is "Enhancing large language models for automated homework assessment in undergraduate circuit analysis"; the id is abs/2502.07980; a note adds "Introduces the CIRCUIT benchmark."
- The title (homework grading) and the note (a benchmark) describe different works. C-S shows that they are two different papers by different authors, neither of them named Liu.

**C9 (Medium). Two self-citations support claims their titles do not address.**
- Xie et al. 2023a (L1175–1180, "The next chapter: A study of large language models in storytelling", INLG 2023) supports "As LLMs expand into high-stakes engineering workflows, rigorous evaluation of their reasoning capabilities has become paramount" (L45–48).
- Xie et al. 2023b (L1181–1186, "DeltaScore: Fine-grained story evaluation with perturbations", Findings of EMNLP 2023) is the source given for "standard reference-based metrics" (L466–468). BERTScore and ROUGE get their own citations in the next sentence (L469–471).
- Zhuohan Xie is an author of EngTrace (L8) and of both papers.

**C10 (Medium). Three of the four "curricular standards" are not curricular standards.**
- L211–215: "we synthesized engineering curricular standards from ABET (ABET, 2024) and major professional societies (ASME (ASME, 2025), IEEE (IEEE, 2025), AIChE (AIChE, 2025))".
- Entries: ASME, "Vision, mission, & core values" (L893–895); IEEE, "IEEE standards association", the standards body's home page (L1030–1031); AIChE, "Constitution & bylaws" (L867–869). Only ABET's "Criteria for accrediting engineering programs" (L863–866) is a curricular standard.

**C11 (Medium). A paper about AC1 is cited for AC2.**
- L365–366 and L2231–2233: "Gwet's AC2 (Gwet, 2008; Wongpakaran et al., 2013)".
- Entry, L1169–1174: Wongpakaran et al. 2013, "A comparison of Cohen's Kappa and Gwet's AC1 when calculating inter-rater reliability coefficients". AC1 is the nominal coefficient and AC2 the weighted, ordinal one. App. T.1 (L2726–2730) itself explains why AC1 rather than AC2 suits nominal data.

**C12 (Low, judgement). The "symbolic" and "qualitative" labels sit oddly with the titles.**
- L175–177: "others target advanced symbolic (Zhang et al., 2025b; Qiu et al., 2025) or qualitative (Wang et al., 2023; Shojaee et al., 2025) understanding".
- Shojaee 2025 (L1136–1140) is "LLM-SRBench: A new benchmark for scientific equation discovery", by its title the most symbolic of the four, yet it is labelled qualitative. The two "symbolic" titles say "high-difficulty and dynamic physics problems" (ABench-Physics, L1233–1235) and "holistic evaluation of physical perception and reasoning" (PHYBench, L1129–1132).
- The next sentence (L177–180) says all six physics benchmarks are "confined to theoretical physics". "Scientific equation discovery" and NEWTON's "physical reasoning" (L1165–1168) do not read as theoretical physics.

**C13 (Low, positioning). An engineering evaluation is cited only as a supporting fact.**
- Heesch et al. 2025 (L1011–1015, "Evaluating large language models for real-world engineering tasks") appears only at L165–168, in the mathematics-and-coding paragraph, as evidence that math and coding skill does not transfer. By its title it is direct prior work on engineering evaluation, yet "Benchmarks in Engineering" (L181–198) does not mention it.

**C14 (Low to Medium, judgement). Weaker fits between the claim and the cited title.**
- a. Chen et al. 2023 (L904–908, "FrugalGPT: How to use large language models while reducing cost and improving performance") is cited for "cascaded model verification" (L391–392). The title is about cost-saving cascades.
- b. Artstein and Poesio 2008 (L890–892, "Inter-coder agreement for computational linguistics") is cited for the claim that, under range restriction, "rank correlation is dominated by tie-breaking noise" (L2352–2354). The title is about agreement coefficients, not rank correlation.
- c. Feinstein and Cicchetti 1990 (L964–967, "High agreement but low kappa: I. The problems of two paradoxes") is cited for AC2's "robustness to range restriction" (L367–368). The title documents kappa's paradox, the problem, rather than AC2's robustness, the remedy. At L2234–2239 it is cited for the behaviour of Krippendorff's α, while its title is about kappa.
- d. Zhao et al. 2023 (L1236–1243, "A survey of large language models") supports the engineering-specific claim at L45–48. General, not wrong.
- e. Zapf et al. 2016 (L1214–1218) is cited for "Fleiss' κ and Krippendorff's α (nominal variant), which are mathematically equivalent for three raters with complete nominal data" (L2731–2734). With complete data, the standard definitions give α = κ_F + (1 − κ_F)/n, where n is the total number of ratings (600 per model here, 6,600 overall). The two are close but not identical, and the gap does not depend on there being three raters. For example, the printed κ_F = 0.183 for Gemini 3.1 Pro (L3236) implies α of 0.184 or 0.185 to three decimals, against "identical for all models" (L2740–2742) and "equal to Krippendorff's α in all cases" (L3249–3250). Whether Zapf et al. state an exact equivalence should be checked; this point is derived from the definitions, not from reading their paper.

**C15 (Low, hygiene). Incomplete or malformed entries.**
- a. Xie et al. 2025, FinChain (L1187–1196): "arXiv preprint." with no identifier; the arXiv record gives 2506.02515.
- b. Team Qwen 2025 (L1135): the whole entry is "Team Qwen. 2025. Qwen3.", with no venue, URL or identifier.
- c. AI Mistral 2024 (L1101–1103): no URL, and the organisation is printed as "AI Mistral" (the text cites "Mistral, 2024").
- d. "Team Qwen" (L1133, L1135) inverts "Qwen Team" (the text cites "Qwen, 2024, 2025").
- e. Anthropic 2026 (L885–889), Google DeepMind 2026 (L982–986) and DeepSeek-AI 2026 (L940–944) are "Technical report" entries with a model identifier but no URL.
- f. Chen et al. 2023 (L904–908) gives the ICML proceedings without pages; Chen et al. 2025 (L915–920) gives "Scientific Reports, 15" without an article number.
- g. The Zhou et al. 2025 note (L1248–1249) names "EngiBench by Felten et al." without an entry.
- h. Acronyms and names lost their capitals (missing BibTeX braces): "Gpt-5 system card" (L1112), "Introducing gpt-4.1 in the api" (L1116), "Supergpqa: Scaling llm evaluation" (L952–953), "Toward engineering agi ... of llms" (L1004–1006), "Deepseek-r1 ... in llms" (L1000–1001), "Deepseek-v3" (L1078–1079), "Qwen2.5-math" (L1206), "Elecbench" (L923), "Eee-bench" (L1057), "Apbench" (L917), "Engibench" (L1245), "Bertscore" (L1220), "Wizardmath" (L1086–1087), "Metamath" (L1211), "Putnam-axiom" (L994), "The llama 3 herd" (L990), "Claude 3.7 sonnet" (L870), "Claude sonnet 4.5" (L878), "Claude opus 4 & claude sonnet 4" (L882–883), "github" (L1039).

### C-S Supplement: printed authors against the arXiv record

This goes beyond the reference list. `notes/arxiv_abs_check.json` is this folder's record of arXiv metadata, written by `check_arxiv_abs.py`. Other work in this folder regenerates the record as it goes. The results below were produced on the snapshot generated 2026-10-01T17:07:23Z (48 records) and came out identical on the one generated at 17:19:49Z (44 records). Both cover 38 papers of this reference list (39 entries, because MMLU has two) and the works the paper names but attributes elsewhere (GLUE, BBH, BIG-bench, CIRCUIT, Felten et al.'s EngiBench). Script 2 compares each printed author list with the arXiv one, position by position, looking records up by arXiv id: surnames must match, and so must the first given name, where a printed initial matches any name with that initial and middle names are ignored. The cited papers themselves were not read. A later run may cover more entries.

**Wrong paper, wrong first author or wrong author list**

| Entry | Printed | arXiv record | Effect |
|---|---|---|---|
| Liu 2025 (L1080) | "Jiaheng Liu et al."; title "Enhancing large language models for automated homework assessment in undergraduate circuit analysis"; id abs/2502.07980; note "Introduces the CIRCUIT benchmark" | The title is arXiv:2511.18221, by Liangliang Chen, Huiru Xie, Zhihao Qin, Yiming Guo, Jacqueline Rohde and Ying Zhang (FIE 2025). The id 2502.07980 is "CIRCUIT: A Benchmark for Circuit Interpretation and Reasoning Capabilities of LLMs", by Lejla Skelic, Yan Xu, Matthew Cox, Wenjie Lu, Tao Yu and Ruonan Han. Neither paper has an author named Liu | The entry merges two papers under an author who wrote neither. The key should be "Chen et al., 2025" or "Skelic et al., 2025", depending on which work is meant |
| Cheng 2024 (L921) | Yuheng Cheng, Huan Zhao, Xiyuan Zhou, Junhua Zhao, Yuji Cao, Chao Yang, Xinlei Cai (no "et al.") | Xiyuan Zhou, Huan Zhao, Yuheng Cheng, Yuji Cao, Gaoqi Liang, Guolong Liu, Wenxuan Liu, Yan Xu, Junhua Zhao | Wrong first author, so the key should read "Zhou et al., 2024". "Chao Yang" and "Xinlei Cai" are not authors, and four authors are missing |
| Lightman 2023 (L1061) | 30 names | 10 authors: Hunter Lightman, Vineet Kosaraju, Yura Burda, Harri Edwards, Bowen Baker, Teddy Lee, Jan Leike, John Schulman, Ilya Sutskever, Karl Cobbe | 23 of the 30 printed names are not authors (for example Lukas Mesnard, Tyna Wang, Farzad Khorrami); Baker, Lee and Cobbe are missing. The key is unaffected |

**Surnames right, given names wrong** (printed → arXiv)

| Entry | Wrong / compared | Names |
|---|---|---|
| Zhang 2025b (L1230) | 12 / 12 | Yifan → Yiming Zhang; Yidong → Yingfan Ma; Yuntian → Yanmei Gu; Ziyi → Zhengkai Yang; Yijia → Yihong Zhuang; Fanglei → Feng Wang; Zhaofeng → Zenan Huang; Yan → Yuanyuan Wang; Chaohui → Chao Huang; Biao → Bowen Song; Chang → Cheng Lin; Jun → Junbo Zhao |
| Luo 2024 (L1090) | 9 / 9 | Haoyang → Hanjun Luo; Hao → Haoyu Huang; Zhaowei → Ziye Deng; Xiaodong → Xinfeng Li; Haoran → Hewei Wang; Yaxin → Yingbin Jin; Yaling → Yang Liu; Weiming → Wenyuan Xu; Zenghao → Zuozhu Liu |
| Xu 2025 (L1197) | 8 / 9 | Xingyu → Xin Xu; Qian → Qiyun Xu; Tao → Tianhao Chen; Yufan → Yuchen Yan; Jifan → Jiaxin Zhang; Shuo → Shizhe Diao; Chen → Can Yang; Yizhong → Yang Wang |
| Zhang 2025a (L1225) | 6 / 9 | Xinyi → Xinyu Zhang; Yutao → Yanrui Wu; Jiahui → Jiaxing Huang; Chen → Chengyou Jia; Li → Lingling Zhang; Jielin → Jun Liu |
| Kazemi 2025 (L1041) | 6 / 20 | Mostafa → Mehran Kazemi; Behnam → Bahare Fatemi; Christos → Chrysovalantis Anastasiou; Shweta V. → Sanket Vaibhav Mehta; Lovish K. → Lalit K. Jain; Valentina → Virginia Aglietti. The 12 names printed as initials agree |
| Shojaee 2025 (L1136) | 4 / 6 | Payam → Parshin Shojaee; Nam-Huan → Ngoc-Hieu Nguyen; Kamyar → Kazem Meidani; Khanh Duy → Khoa D Doan |
| Qiu 2025 (L1129) | 3 / 3 printed | Siyu → Shi Qiu; Shuo → Shaoyang Guo; Ze-Yu → Zhuo-Yang Song |
| Mirzadeh 2024 (L1096) | 2 / 6 | Kiarash → Keivan Alizadeh; Hidetoshi → Hooman Shahrokhi |
| Wang 2023 (L1165) | 1 / 4 | Jia → Jiafei Duan |
| Jimenez 2024 (L1032) | 1 / 7 | Shun → Shunyu Yao |

**Version difference, not counted as an error.** Gulati 2024 (L992): the seven printed names match the arXiv list in order, except that the arXiv version (first posted August 2025) has an eighth author, Elyas Obbad, before Sanmi Koyejo. The entry cites the 2024 NeurIPS workshop version, which may have listed seven.

**Agree with the record** (24 papers, 25 entries): Austin 2021, Chen 2021, Chen 2023, Cobbe 2021, Deng 2024, Du 2025 (arXiv lists "M-A-P Team" first, so Script 2 shows each name one place later), Fan 2024, Guo 2025b, Heesch 2025, Hendrycks 2020 and 2021a, Hendrycks 2021b, Li 2025, Mudur 2025, Ott 2022, Plank 2022, Syed 2024, Verga 2024, Wang 2019, Xie 2023a, Xie 2023b, Xie 2025, Zhang 2020, Zhao 2023, Zhou 2025.

**Not in the record** (40 entries, not checked): the 22 model entries; the 9 textbooks and standards; 8 of the 10 metric and statistics entries (all but Zhang 2020 and Xie 2023b); and Chen 2025, APBench, which appeared in Scientific Reports and has no arXiv id in the folder's catalogue.

In ten entries the surnames match while the given names do not, and two more list people who are not authors. The pattern indicates these entries were not copied from the source records. Regenerating the entries from arXiv, the ACL Anthology or Crossref (this folder's `make_bib.py` builds such records) would fix C-S and much of C15; C5 still needs the duplicate removed, and C8 needs the authors to say which work they meant.

## Part D: The paper's self-description

### D.1 What the paper says is new (15 quotes)

1. **Symbolic benchmark; verifiable process supervision.** "EngTrace: A Symbolic Benchmark for Verifiable Process Supervision of Engineering Reasoning" (title, L5–6).
2. **Symbolic templates; contamination resistance.** "a symbolic benchmark built on 90 parameterized templates, each generating unique, contamination-resistant problem instances" (abstract, L25–28).
3. **Tiered protocol; AI Tribunal.** "Moving beyond outcome matching, we introduce a verifiable two-stage evaluation framework that uses a tiered protocol to validate intermediate reasoning traces alongside final answers through automated procedural checks and a heterogeneous AI Tribunal." (abstract, L32–37)
4. **Complexity cliff; integrative reasoning.** "Our evaluation of 27 leading LLMs reveals a distinct trade-off between numeric precision and trace fidelity, identifying a complexity cliff where abstract mathematical pre-training fails to translate into the integrative reasoning required for advanced engineering tasks." (abstract, L37–43)
5. **The gap.** "Engineering reasoning requires physically grounded, multi-step problem solving where intermediate steps are as consequential as the final result, yet no existing benchmark verifies this process." (§1 ¶2, L77–80)
6. **Gold traces; templates credited to prior work.** "We introduce EngTrace: a symbolic benchmark pairing each problem instance with a gold-standard reasoning trace to enable verifiable process supervision. Symbolic templates (Mirzadeh et al., 2024; Xie et al., 2025) serve as the generative mechanism, producing unique, contamination-resistant instances" (§1 ¶2, L81–87)
7. **Contamination resistance.** "EngTrace counters this through symbolic templates and domain-aware parameterization, generating unique, physically grounded problems that resist rote memorization." (§1 ¶3, L101–105)
8. **Integrative reasoning.** "EngTrace addresses both gaps through integrative reasoning templates that ... demand holistic procedural reasoning, a crucial requirement in engineering, where a flawed process can lead to catastrophic failure (Leveson, 2012)." (§1 ¶3, L116–122)
9. **Scope, tied to gold traces.** "Importantly, EngTrace is scoped as a depth-oriented evaluation of verifiable process supervision across a principled engineering subset rather than full disciplinary coverage or a production agent workflow; this scope is a deliberate design choice, as verifiability requires gold-standard reasoning traces that can only be guaranteed through controlled symbolic generation." (§1 ¶3, L122–129)
10. **Templates beyond variable substitution.** "creating a dynamic problem space that effectively resists data contamination. Critically, templates go beyond variable substitution: domain-aware parameter sampling can alter the physical regime, governing equation, and required reasoning path entirely (see Appendix G), though linguistic diversity through paraphrased or naturally sourced wording remains future work." (§3.2, L247–254)
11. **AI Tribunal, credited to PoLL.** "Building upon the “Panel of LLM” (PoLL) methodology (Verga et al., 2024), we construct an “AI Tribunal” comprising three frontier reasoning models" (§3.3, L260–263)
12. **Tiered verification protocol, credited.** "Inspired by recent advances in process supervision (Lightman et al., 2023) and cascaded model verification (Chen et al., 2023), we introduce a two-stage evaluation framework that validates intermediate reasoning steps in addition to final-answer accuracy through a tiered verification protocol." (§4, L390–395)
13. **Validation of the evaluator.** "We design and validate a two-stage evaluation framework combining automated symbolic verification with a heterogeneous AI Tribunal to assess reasoning process quality beyond final answer accuracy, achieving stronger expert alignment (ρ = 0.632) than standard NLP metrics (ρ = 0.325)." (§1, contribution 2, L138–144). §4 adds that the Tribunal-based approach "significantly outperforms standard metrics in capturing engineering reasoning quality" (L399–400).
14. **Error taxonomy; complexity cliff.** "We evaluate 27 LLMs across frontier, open-weights, and math-enhanced categories, and introduce a six-category error taxonomy applied to 2,200 annotated failure traces, revealing qualitatively distinct failure modes across model classes and a complexity cliff that abstract mathematical pre-training cannot bridge." (§1, contribution 3, L145–151). §5.4 adds: "this analysis employs a finer six-category taxonomy applied by human domain experts to diagnose the root causes of confirmed reasoning failures independently of the automated evaluation pipeline" (L761–765).
15. **The paper's own caveat.** "Second, EngTrace is entirely synthetic. While this ensures verifiability and contamination resistance, symbolic templates produce linguistically constrained problem statements that do not capture the full diversity or contextual ambiguity of real-world engineering." (Limitations, L829–835)

### D.2 Sentences that contrast EngTrace with named prior benchmarks

- Abstract, L17–23: "However, existing benchmarks such as MMLU, MATH, and HumanEval assess isolated cognitive skills, failing to capture the physically grounded reasoning central to engineering, where scientific principles, quantitative modeling, and practical constraints must converge."
- §1 ¶3, L97–105: "First, static benchmarks like Glue (Wang et al., 2019) and BBH (Kazemi et al., 2025) are increasingly prone to model saturation (Ott et al., 2022) and data contamination (Deng et al., 2024); EngTrace counters this through symbolic templates and domain-aware parameterization, generating unique, physically grounded problems that resist rote memorization."
- §1 ¶3, L105–111: "Second, most benchmarks evaluate skills in disciplinary silos: general knowledge benchmarks (e.g., MMLU (Hendrycks et al., 2020)) test broad factual recall, while specialized benchmarks (e.g., Math (Hendrycks et al., 2021b) or HumanEval (Chen et al., 2021)) assess abstract logic and algorithmic translation."
- §1 ¶3, L116–122: "EngTrace addresses both gaps through integrative reasoning templates that, unlike existing benchmarks (e.g., Zhou et al., 2025; Mudur et al., 2025) that rely solely on outcome matching, demand holistic procedural reasoning, a crucial requirement in engineering, where a flawed process can lead to catastrophic failure (Leveson, 2012)."
- §2 Math, L165–168: "Yet neither abstract logical deduction nor algorithmic proficiency transfers to reliable reasoning under physical constraints (Heesch et al., 2025)." (aimed at GSM8K, MATH, HARDMath, Putnam-AXIOM, GSM-Symbolic, HumanEval, MBPP and SWE-bench)
- §2 Phys, L177–180: "Despite this progress, these benchmarks remain confined to theoretical physics, lacking the integrative context of engineering." (aimed at UGPhysics, PhysReason, ABench-Physics, PHYBench, NEWTON and LLM-SRBench)
- §2 Eng, L182–190: "Generalist suites like MMLU (Hendrycks et al., 2021a), BIG-Bench (Luo et al., 2024), and SuperGPQA (Du et al., 2025) restrict assessment to factual recall via multiple-choice questions, while specialized benchmarks like TransportBench, APBench, and EEE-Bench (Syed et al., 2024; Chen et al., 2025; Cheng et al., 2024; Li et al., 2025; Liu et al., 2025) are siloed within narrow sub-disciplines."
- §2 Eng, L190–198: "Others focus on tangential skills such as design generation (Guo et al., 2025b) or software proficiency (Mudur et al., 2025), while broader efforts like EngiBench (Zhou et al., 2025) rely on subjective rubric-based scoring."
- §2 closing, L199–205: "Evidently, existing benchmarks remain limited to outcome matching and fail to validate the physically grounded reasoning required for engineering. EngTrace addresses this limitation by introducing a symbolic benchmark that pairs unique problem instances with gold-standard reasoning traces to enable verifiable process supervision."
- App. U.1, L2844–2847: "Two models with identical aggregate error rates thus represent fundamentally different failure modes, a distinction aggregate benchmarks cannot surface."

Implicit or method-level contrasts, with no benchmark named: "Critically, templates go beyond variable substitution" (§3.2, L249–250), which reads as a contrast with GSM-Symbolic-style templates; "Unlike standard string matching, this alignment must account for structural variations such as step merging, reordering, or alternative valid methodologies" (§4.1, L416–419); and "the poor performance of standard NLP metrics (ρ = 0.325) reinforces the need for the domain-specific numerical verification employed in EngTrace" (App. M, L2420–2423).

### D.3 Notes for the panel on positioning

1. **The closest precedent is cited only in passing.** L84–85 cites Xie et al. 2025 (FinChain) for symbolic templates and nothing else. Its title, "FINCHAIN: A symbolic benchmark for verifiable chain-of-thought financial reasoning" (L1194–1196), has the same structure as EngTrace's, and it shares four authors with EngTrace: Zhuohan Xie, Fan Zhang, Veselin Stoyanov and Preslav Nakov (L7–8; L1187–1194). Related Work does not discuss it. The claim at L77–80 ("no existing benchmark verifies this process") is scoped to engineering, so FinChain does not contradict it, but a reviewer will ask what EngTrace adds beyond the change of domain.
2. **Related Work has no paragraph on process supervision.** Its three paragraphs cover math and coding, physical-science and engineering benchmarks. The works on step-level verification and LLM judging (Lightman et al. 2023, Verga et al. 2024, Chen et al. 2023) appear only in §3.3 and §4, although "verifiable process supervision" is the claim in the title.
3. **"Outcome matching" is the dividing line, and two of its examples do not fit it cleanly.** L116–119 and L199–201 place all prior engineering benchmarks on the outcome-matching side, but L193–198 says EngiBench uses rubric scoring (C6), and FEABench's title claims multiphysics reasoning (C7).
4. **"Complexity cliff" is measured two ways.** Figure 4 (L730–734), the conclusion (L803–806) and App. S.3 (L2691–2695: "52.55% Easy to 23.33% Advanced, a 29-point drop") use the fall in final-answer accuracy from Easy to Advanced, contrasting open-weights with frontier models. §5.4 (L790–797: "cliff = −4.1 pp", "a 41.8 pp conceptual cliff") and App. U.1 (L2806–2810) define it as "the percentage-point change in conceptual error rate from Easy to Advanced difficulty". The abstract (L40–43) and contribution 3 (L150–151) tie the cliff to abstract mathematical pre-training.

## How this was produced

1. The whole dump was read: main text, reference list and appendices. Appendix G's template code was scanned for author–year patterns and has none.
2. **Entry count.** The 79 entries were listed by hand and cross-checked with a line-start heuristic (`awk` over L863–1249: a capitalised line that follows a full stop, a page break or the heading). The heuristic returns 100 candidates; 21 are continuation lines.
3. **Citation instances, orphans and danglings.** Script 1 below joins every non-reference line with all whitespace removed (so "Mirzadehetal.,2024" and "Mirzadeh et al., 2024" look the same), drops glyph lines, page markers and page numbers, maps each character back to its line, and searches one pattern per entry. It finds 90 instances and no entry with zero. Every author–year pattern in the same text was then listed and resolved to an entry by hand; none was left over.
4. **Quotes and paragraphs.** Sentences were taken from the PDF's text layer (PyMuPDF `page.get_text`), which keeps the spaces the dump loses. Paragraph boundaries were read from first-line indents: L77, L95, L130 and L199 are indented in the PDF.
5. **Author supplement (C-S).** Script 2 below compares the printed author lists, copied from the PDF text layer, with `abs_authors` in `notes/arxiv_abs_check.json`, looking each paper up by arXiv id. The numbers in C-S come from the snapshot generated 2026-10-01T17:07:23Z and were unchanged on the 17:19:49Z snapshot.

Both scripts run from the repository root with Python 3 and the standard library.

<details>
<summary>Script 1: citation instances per entry</summary>

```python
import re
L = [""] + open("docs/_ARR_May__EngTrace.txt", encoding="utf-8").read().split("\n")
SECTIONS = [(1,43,"Abstract"),(44,52,"§1 ¶1"),(53,76,"§1 Fig.1"),(77,94,"§1 ¶2"),(95,129,"§1 ¶3"),(130,151,"§1 contributions"),
 (152,168,"§2 Math"),(169,180,"§2 Phys"),(181,198,"§2 Eng"),(199,205,"§2 closing"),(206,235,"§3.1"),(236,254,"§3.2"),
 (255,378,"§3.3"),(379,388,"§3.4"),(389,401,"§4"),(402,448,"§4.1"),(449,464,"§4.2"),(465,474,"§4.3"),(475,549,"§4.4"),
 (550,583,"§5.1"),(584,599,"§5.2"),(600,756,"§5.3"),(757,799,"§5.4"),(800,816,"§6"),(817,847,"Limitations"),(848,861,"Ethics"),
 (862,1250,"References"),(1251,2226,"App. A–K.1"),(2227,2288,"App. K.2"),(2289,2338,"App. K.3"),(2339,2351,"App. L"),
 (2352,2356,"App. K.3 fn1"),(2357,2708,"App. L–S"),(2709,2735,"App. T.1"),(2736,3276,"App. T.2–U, tables")]
section = lambda n: next(s for a, b, s in SECTIONS if a <= n <= b)
buf, cmap = [], []
for n in range(1, len(L)):
    if 862 <= n <= 1250:
        buf.append("|"); cmap.append(n); continue
    if L[n].startswith("/uni") or L[n].startswith("===== PAGE") or re.fullmatch(r"\s*\d{1,2}\s*", L[n]):
        continue
    for ch in L[n]:
        if not ch.isspace():
            buf.append(ch); cmap.append(n)
S = "".join(buf)
E = r"etal\.,"
KEYS = [("ABET 2024","ABET,2024"),("AIChE 2025","AIChE,2025"),("Anthropic 2025a","Anthropic,2025c,d,a"),
 ("Anthropic 2025b","An-?thropic,2025b"),("Anthropic 2025c","Anthropic,2025c"),("Anthropic 2025d","Anthropic,2025c,d"),
 ("Anthropic 2026","Anthropic,2025c,d,a,2026"),("Artstein 2008","ArtsteinandPoesio,2008"),("ASME 2025","ASME,2025"),
 ("Austin 2021","Austin"+E+"2021"),("Byrt 1993","Byrt"+E+"1993"),("Chen 2023","Chen"+E+"2023"),("Chen 2021","Chen"+E+"2021"),
 ("Chen 2025","Chen"+E+"2025"),("Cheng 2024","Cheng"+E+"2024"),("Cobbe 2021","Cobbe"+E+"2021"),("Comanici 2025","Comanici"+E+"2025"),
 ("DeepSeek-AI 2026","DeepSeek-?AI,2026"),("Deng 2024","Deng"+E+"2024"),("Du 2025","Du"+E+"2025"),("Dym 2012","Dym"+E+"2012"),
 ("Fan 2024","Fan"+E+"2024"),("Feinstein 1990","FeinsteinandCicchetti,1990"),("Fogler 2020","Fogler,2020"),
 ("Gemma Team and GDM 2025","GemmaTeamandGoogleDeep-?Mind,2025"),("Gemma Team 2024","GemmaTeam"+E+"2024"),
 ("Google 2025","(?<!DeepMind,2026;)Google,2025|DeepMind,2026;Google,2025"),("Google DeepMind 2026","GoogleDeepMind,2026"),
 ("Grattafiori 2024","Grattafiori"+E+"2024"),("Gulati 2024","Gulati"+E+"2024"),("Guo 2025a","Guo"+E+"2025a"),("Guo 2025b","Guo"+E+"2025b"),
 ("Gwet 2008","Gwet,2008"),("Heesch 2025","Heesch"+E+"2025"),("Hendrycks 2020","Hendrycks"+E+"2020"),("Hendrycks 2021a","Hendrycks"+E+"2021a"),
 ("Hendrycks 2021b","Hendrycks"+E+"2021b"),("IEEE 2025","IEEE,2025"),("Jimenez 2024","Jimenez"+E+"2024"),("Kazemi 2025","Kazemi"+E+"2025"),
 ("Landis 1977",r"LandisandKoch\(1977\)"),("Leveson 2012","Leveson,2012"),("Li 2025","Li"+E+"2025"),("Lightman 2023","Lightman"+E+"2023"),
 ("Lin 2004","Lin,2004"),("Liu 2024","Liu"+E+"2024"),("Liu 2025","Liu"+E+"2025"),("Luo 2023","Luo"+E+"2023"),("Luo 2024","Luo"+E+"2024"),
 ("Mirzadeh 2024","Mirzadeh"+E+"2024"),("Mistral 2024","Mistral,2024"),("Moaveni 2019","Moaveni,2019"),("Mudur 2025","Mudur"+E+"2025"),
 ("OpenAI 2025a","Ope-?nAI,2025a"),("OpenAI 2025b","Ope-?nAI,2025a,b"),("Ott 2022","Ott"+E+"2022"),("Plank 2022","Plank,2022"),
 ("Qiu 2025","Qiu"+E+"2025"),("Qwen 2024","Qwen,2024"),("Qwen 2025","Qwen,2024,2025"),("Shojaee 2025","Shojaee"+E+"2025"),
 ("Syed 2024","Syed"+E+"2024"),("Ulaby 2015","UlabyandRavaioli,2015"),("Verga 2024","Verga"+E+"2024"),("Wang 2019","Wang"+E+"2019"),
 ("Wang 2023","Wang"+E+"2023"),("Wongpakaran 2013","Wongpakaran"+E+"2013"),("Xie 2023a","Xie"+E+"2023a"),("Xie 2023b","Xie"+E+"2023b"),
 ("Xie 2025","Xie"+E+"2025"),("Xu 2025","Xu"+E+"2025"),("Yang 2024","Yang"+E+"2024"),("Yu 2023","Yu"+E+"2023"),("Zapf 2016","Zapf"+E+"2016"),
 ("Zhang 2020","Zhang"+E+"2020"),("Zhang 2025a","Zhang"+E+"2025a"),("Zhang 2025b","Zhang"+E+"2025b"),("Zhao 2023","Zhao"+E+"2023"),
 ("Zhou 2025","Zhou"+E+"2025")]
total = 0
for label, rx in KEYS:
    hits = [(cmap[m.start()], cmap[m.end() - 1]) for m in re.finditer(rx, S)]
    total += len(hits)
    where = "; ".join(f"L{a}" + (f"-{b}" if b != a else "") + f" [{section(a)}]" for a, b in hits)
    print(f"{label:26s} n={len(hits)}  {where or 'ORPHAN'}")
print(len(KEYS), "entries;", total, "citation instances")
```

</details>

<details>
<summary>Script 2: printed authors against the arXiv record</summary>

```python
import json, re, unicodedata
DOC = json.load(open("docs/related_work_oct2026/notes/arxiv_abs_check.json", encoding="utf-8"))
BY_ID = {}
for rec in DOC["entries"].values():
    BY_ID.setdefault(rec.get("arxiv"), rec)
PRINTED = [  # (entry, arXiv id of the paper the entry names, authors as printed in the May reference list)
 ("Austin 2021", "2108.07732", "Jacob Austin, Augustus Odena, Maxwell Nye, Maarten Bosma, Henryk Michalewski, David Dohan, Ellen Jiang, Carrie Cai, Michael Terry, Quoc Le, et al."),
 ("Chen 2021", "2107.03374", "Mark Chen, Jerry Tworek, Heewoo Jun, Qiming Yuan, Henrique Ponde de Oliveira Pinto, Jared Kaplan, Harri Edwards, Yuri Burda, Nicholas Joseph, Greg Brockman, et al."),
 ("Chen 2023", "2305.05176", "Lingjiao Chen, Matei Zaharia, and James Zou"),
 ("Cheng 2024", "2407.05365", "Yuheng Cheng, Huan Zhao, Xiyuan Zhou, Junhua Zhao, Yuji Cao, Chao Yang, and Xinlei Cai"),
 ("Cobbe 2021", "2110.14168", "Karl Cobbe, Vineet Kosaraju, Mohammad Bavarian, et al."),
 ("Deng 2024", "2311.09783", "Chunyuan Deng, Yilun Zhao, Xiangru Tang, Mark Gerstein, and Arman Cohan"),
 ("Du 2025", "2502.14739", "Xinrun Du, Yifan Yao, Kaijing Ma, Bingli Wang, Tianyu Zheng, et al."),
 ("Fan 2024", "2410.09988", "J. Fan, S. Martinson, E. Y. Wang, K. Hausknecht, J. Brenner, D. Liu, N. Peng, C. Wang, and M. P. Brenner"),
 ("Gulati 2024", "2508.08292", "Aryan Gulati, Brando Miranda, Eric Chen, Emily Xia, Kai Fronsdal, Bruno de Moraes Dumont, and Sanmi Koyejo"),
 ("Guo 2025b", "2509.16204", "Xingang Guo, Yaxin Li, Xiangyi Kong, Yilan Jiang, Xiayu Zhao, et al."),
 ("Heesch 2025", "2505.13484", "Rene Heesch, Sebastian Eilermann, Alexander Windmann, Alexander Diedrich, Philipp Rosenthal, and Oliver Niggemann"),
 ("Hendrycks 2020/2021a", "2009.03300", "Dan Hendrycks, Collin Burns, Steven Basart, Andy Zou, Mantas Mazeika, Dawn Song, and Jacob Steinhardt"),
 ("Hendrycks 2021b", "2103.03874", "Dan Hendrycks, Collin Burns, Saurav Kadavath, Akul Arora, Steven Basart, Eric Tang, Dawn Song, and Jacob Steinhardt"),
 ("Jimenez 2024", "2310.06770", "Carlos E. Jimenez, John Yang, Alexander Wettig, Shun Yao, Kexin Pei, Ofir Press, and Karthik R. Narasimhan"),
 ("Kazemi 2025", "2502.19187", "Mostafa Kazemi, Behnam Fatemi, Hritik Bansal, John Palowitch, Christos Anastasiou, Shweta V. Mehta, Lovish K. Jain, Valentina Aglietti, D. Jindal, P. Chen, N. Dikkala, G. Tyen, X. Liu, U. Shalit, S. Chiappa, K. Olszewska, Y. Tay, V. Q. Tran, Q. V. Le, and O. Firat"),
 ("Li 2025", "2411.01492", "Ming Li, Jike Zhong, Tianle Chen, Yuxiang Lai, and Konstantinos Psounis"),
 ("Lightman 2023", "2305.20050", "Hunter Lightman, Vineet Kosaraju, Yura Burda, Harri Edwards, Lukas Mesnard, Tyna Wang, Farzad Khorrami, Nguyet Minh Nguyen, Shayne Mostyn, Max Miller, Chia Hsuan Wang, Sam Gelman, Denis Igor, Marina Polozov, Girish Sastry, Prachit Tara, Sandhini Agarwal, Nan Rosemary Sun, Shengjia Zhao, Jeffrey Wu, Szymon Sidor, Jiayi Weng, Yuan Cao, Adrià Puigdomènech Badia, Nikolas Tezak, Peter Welinder, Ilya Sutskever, John Schulman, Jan Leike, and Wojciech Zaremba"),
 ("Liu 2025 (printed title)", "2511.18221", "Jiaheng Liu et al."),
 ("Liu 2025 (printed id)", "2502.07980", "Jiaheng Liu et al."),
 ("Luo 2024", "2407.15240", "Haoyang Luo, Hao Huang, Zhaowei Deng, Xiaodong Li, Haoran Wang, Yaxin Jin, Yaling Liu, Weiming Xu, and Zenghao Liu"),
 ("Mirzadeh 2024", "2410.05229", "Iman Mirzadeh, Kiarash Alizadeh, Hidetoshi Shahrokhi, Oncel Tuzel, Samy Bengio, and Mehrdad Farajtabar"),
 ("Mudur 2025", "2504.06260", "Nayantara Mudur, Hao Cui, Subhashini Venugopalan, Paul Raccuglia, Michael P. Brenner, and Peter Norgaard"),
 ("Ott 2022", "2203.04592", "S. Ott, A. Barbosa-Silva, K. Blagec, J. Brauner, and M. Samwald"),
 ("Plank 2022", "2211.02570", "Barbara Plank"),
 ("Qiu 2025", "2504.16074", "Siyu Qiu, Shuo Guo, Ze-Yu Song, et al."),
 ("Shojaee 2025", "2504.10415", "Payam Shojaee, Nam-Huan Nguyen, Kamyar Meidani, Amir Barati Farimani, Khanh Duy Doan, and Chandan K. Reddy"),
 ("Syed 2024", "2408.08302", "Usman Syed, Ethan Light, Xingang Guo, Huan Zhang, Lianhui Qin, Yanfeng Ouyang, and Bin Hu"),
 ("Verga 2024", "2404.18796", "Pat Verga, Sebastian Hofstatter, Sophia Althammer, Yixuan Su, Aleksandra Piktus, Arkady Arkhangorodsky, Minjie Xu, Naomi White, and Patrick Lewis"),
 ("Wang 2019", "1905.00537", "Alex Wang, Yada Pruksachatkun, Nikita Nangia, Amanpreet Singh, Julian Michael, Felix Hill, Omer Levy, and Samuel R. Bowman"),
 ("Wang 2023", "2310.07018", "Yi R. Wang, Jia Duan, Dieter Fox, and Siddhartha S. Srinivasa"),
 ("Xie 2023a", "2301.09790", "Zhuohan Xie, Trevor Cohn, and Jey Han Lau"),
 ("Xie 2023b", "2303.08991", "Zhuohan Xie, Miao Li, Trevor Cohn, and Jey Lau"),
 ("Xie 2025", "2506.02515", "Zhuohan Xie, Daniil Orel, Rushil Thareja, Dhruv Sahnan, Hachem Madmoun, Fan Zhang, Debopriyo Banerjee, Georgi Georgiev, Xueqing Peng, Lingfei Qian, Jimin Huang, Jinyan Su, Aaryamonvikram Singh, Rui Xing, Rania Elbadry, Chen Xu, Haonan Li, Fajri Koto, Ivan Koychev, Tanmoy Chakraborty, Yuxia Wang, Salem Lahlou, Veselin Stoyanov, Sophia Ananiadou, and Preslav Nakov"),
 ("Xu 2025", "2502.00334", "Xingyu Xu, Qian Xu, Tong Xiao, Tao Chen, Yufan Yan, Jifan Zhang, Shuo Diao, Chen Yang, and Yizhong Wang"),
 ("Zhang 2020", "1904.09675", "Tianyi Zhang, Varsha Kishore, Felix Wu, Kilian Q. Weinberger, and Yoav Artzi"),
 ("Zhang 2025a", "2502.12054", "Xinyi Zhang, Yuxuan Dong, Yutao Wu, Jiahui Huang, Chen Jia, Basura Fernando, Mike Zheng Shou, Li Zhang, and Jielin Liu"),
 ("Zhang 2025b", "2507.04766", "Yifan Zhang, Yidong Ma, Yuntian Gu, Ziyi Yang, Yijia Zhuang, Fanglei Wang, Zhaofeng Huang, Yan Wang, Chaohui Huang, Biao Song, Chang Lin, and Jun Zhao"),
 ("Zhao 2023", "2303.18223", "Wayne Xin Zhao, Kun Zhou, Junyi Li, Tianyi Tang, Xiaolei Wang, Yupeng Hou, Yingqian Min, Beichen Zhang, Junjie Zhang, Zican Dong, Yifan Du, Chen Yang, Yushuo Chen, Zhipeng Chen, Jinhao Jiang, Ruiyang Ren, Yifan Li, Xinyu Tang, Zikang Liu, Peiyu Liu, Jian-Yun Nie, and Ji-Rong Wen"),
 ("Zhou 2025", "2509.17677", "Xiyuan Zhou, Xinlei Wang, Yirui He, Yang Wu, Ruixi Zou, et al."),
]
def toks(s):
    s = "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))
    return re.sub(r"[^a-z ]", "", s.lower().replace("-", " ")).split()
def same(printed, arxiv):
    """Same surname (last token) and same first given name; a printed initial matches any name with that initial."""
    p, a = toks(printed), toks(arxiv)
    if not p or not a or p[-1] != a[-1]:
        return False
    if len(p) == 1 or len(a) == 1:
        return True
    return p[0] == a[0] if len(p[0]) > 1 else p[0] == a[0][0]
print("record generated", DOC["generated_utc"], "with", len(DOC["entries"]), "records")
for label, aid, s in PRINTED:
    rec = BY_ID.get(aid)
    if rec is None:
        print(f"{label:26s} arXiv {aid}: not in the record"); continue
    etal = s.endswith("et al.")
    names = [x.strip() for x in re.sub(r",?\s*et al\.$", "", s).replace(", and ", ", ").replace(" and ", ", ").split(",") if x.strip()]
    arx = rec["abs_authors"]
    bad = [(i + 1, n, arx[i] if i < len(arx) else "(none)") for i, n in enumerate(names) if i >= len(arx) or not same(n, arx[i])]
    print(f"{label:26s} arXiv {aid}: printed {len(names)}{'+et al.' if etal else ''}, arXiv {len(arx)}: " + ("agree" if not bad else f"{len(bad)} differ"))
    for pos, n, a in bad:
        elsewhere = [j + 1 for j, x in enumerate(arx) if same(n, x)]
        note = f"  (printed name is arXiv author #{elsewhere[0]})" if elsewhere else ""
        print(f"    #{pos}: {n}  ->  {a}{note}")
```

</details>
