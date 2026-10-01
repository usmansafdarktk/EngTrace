# Literature search: LLM-as-a-judge validity and meta-evaluation (topic `llm_judge`)

Search run 2026-10-01 for the October 2026 Related Work revision. Catalogue:
`docs/related_work_oct2026/candidates/llm_judge.json` (30 entries: 11 must-cite, 12 should-cite,
7 optional). Download: `python docs/related_work_oct2026/fetch_papers.py --file llm_judge.json`
fetched 30 of 30, 0 failed, 0 unresolved.

Scope: self-preference and family bias, panels of judges, meta-evaluation against human experts
(technical domains first), judge calibration and reporting, the case for deterministic
verification, and evidence that reference-based overlap metrics do not measure reasoning quality.
Already cited and not re-added: Verga et al. 2024 (PoLL), Lightman et al. 2023, BERTScore, ROUGE,
FrugalGPT, Plank 2022.

## How the search was run, and what limited it

- **WebSearch was not available.** The session's shared budget of 200 WebSearch calls was already
  used up, presumably by the parallel literature agents, before this agent's first query. All six
  queries attempted (listed below) came back "not performed". The brief asked for at least 20; none
  ran. Discovery therefore used the arXiv API through curl (18 queries), OpenAlex search
  (5 queries) and targeted abstract-page checks. Venues came from arXiv comments, PDF headers, the
  OpenReview API and the Semantic Scholar batch API.
- **The arXiv API was heavily rate-limited** (HTTP 429/503, probably from the other agents). Each
  query was retried with backoff. q15, q16 and q18 succeeded only on late retries, after the
  selection was drafted. Their results were screened and changed nothing in it.
- Every included paper's arXiv abstract page was opened with WebFetch (44 pages in all). Its PDF
  was downloaded, and specific claims were checked against the extracted text in
  `docs/related_work_oct2026/text/`.

## Queries (verbatim)

### WebSearch (attempted; all returned "Web search was not performed: this session has used its web search budget (200 of 200)")

1. `"Judging LLM-as-a-judge with MT-Bench and Chatbot Arena" position verbosity self-enhancement bias agreement`
2. `"Preference Leakage" contamination problem in LLM-as-a-judge arXiv 2502.01534 venue`
3. `JudgeBench benchmark for evaluating LLM-based judges objective correctness ICLR 2025`
4. `RewardBench 2 advancing reward model evaluation arXiv 2025`
5. `JETTS benchmark LLM-as-judges test-time scaling evaluators ICML 2025`
6. `LLM-as-a-judge self-preference family bias judges favor models from same family 2026`

### arXiv API (curl)

URL template:
`http://export.arxiv.org/api/query?search_query=<Q>&start=0&max_results=100&sortBy=submittedDate&sortOrder=descending`.
Results were parsed with `xml.etree` and records first posted before 2024 were dropped. The query
strings below are exactly as sent (URL-encoded), with the number of entries returned in brackets.
Queries that hit 100 returned only the 100 most recent records, almost all from 2026. Older papers
for those queries were found by direct checks of known titles and IDs.

| id | search_query (sent) | returned |
|---|---|---|
| q01 | `all:%22LLM-as-a-judge%22+AND+all:bias` | 100 |
| q02 | `all:%22self-preference%22+AND+all:LLM` | 36 |
| q03 | `all:%22judge+reliability%22+OR+%28ti:judge+AND+abs:reliability+AND+abs:LLM%29` | 100 |
| q04 | `all:%22meta-evaluation%22+AND+all:judge+AND+all:LLM` | 51 |
| q05 | `%28all:%22panel+of+judges%22+OR+all:%22multiple+judges%22+OR+all:jury%29+AND+all:LLM` | 87 |
| q06 | `all:judge+AND+all:%22human+experts%22+AND+all:agreement+AND+all:LLM` | 21 |
| q07 | `all:%22LLM-as-a-judge%22+AND+%28all:math+OR+all:mathematical+OR+all:scientific%29` | 100 |
| q08 | `ti:%22LLM-as-a-judge%22` | 100 |
| q09 | `%28all:%22preference+leakage%22+OR+all:%22family+bias%22+OR+all:%22self-bias%22%29+AND+all:judge` | 5 |
| q10 | `all:verifier+AND+all:%22rule-based%22+AND+all:%22model-based%22+AND+all:reasoning` | 10 |
| q11 | `ti:%22LLM+judges%22+OR+ti:%22LLM+evaluators%22` | 100 |
| q12 | `all:%22LLM-as-a-judge%22+AND+all:calibration` | 83 |
| q13 | `all:%22reasoning+traces%22+AND+all:judge+AND+%28all:step+OR+all:process%29` | 37 |
| q14 | `ti:verifier+AND+%28ti:math+OR+ti:mathematical+OR+ti:answer%29` | 74 |
| q15 | `all:BERTScore+AND+all:reasoning+AND+all:%22human+judgments%22` | 7 (failed twice; succeeded on a final retry) |
| q16 | `%28all:%22LLM-as-a-judge%22+OR+all:%22LLM+judge%22%29+AND+%28all:physics+OR+all:engineering+OR+all:chemistry%29` | 100 (failed twice; succeeded on a final retry; mostly software-engineering agent papers) |
| q17 | `all:%22self-recognition%22+AND+all:LLM` | 12 |
| q18 | `%28all:%22LLM-as-a-judge%22+OR+all:%22LLM+judges%22%29+AND+all:family+AND+all:bias` | 34 (succeeded on a late retry) |

834 unique records from 2024 to 2026 were screened by title. About 80 abstracts were read in full.

### Other search and verification calls

- **OpenAlex** (`https://api.openalex.org/works?search=<q>&filter=from_publication_date:2024-01-01&per-page=25`):
  `LLM-as-a-judge engineering`; `LLM judge physics grading human markers`;
  `BERTScore reasoning evaluation chain-of-thought`; `self-preference bias LLM judge`;
  `panel of LLM judges correlated errors`. Results were noisy. They surfaced Kortemeyer et al. 2024
  on thermodynamics-exam grading and the SE-domain judge study (both listed below).
- **Semantic Scholar search**, five queries, all HTTP 429 with no results:
  `LLM-as-a-judge engineering problem solving agreement with experts`;
  `LLM judge physics chemistry grading agreement human experts`;
  `BERTScore ROUGE fail to evaluate reasoning chains`; `self-preference family bias LLM judge 2026`;
  `panel of LLM judges correlated errors`.
- **Semantic Scholar batch** (`POST /graph/v1/paper/batch`, 40 arXiv IDs; fields: venue,
  publicationVenue, externalIds, journal): used to cross-check venues.
- **OpenReview API** (`https://api2.openreview.net/notes/search?term=<title>&type=terms&content=title&source=forum`):
  used to check venues for 35 titles, including every included paper, PoLL, RewardBench 2,
  RM-Bench and ROSCOE.
- **DBLP** (`https://dblp.org/search/publ/api?q=<title>&format=json`), 20 title queries. Every
  response was an anti-bot challenge page, so DBLP could not be used.
- **WebFetch of arXiv abstract pages:** 2306.05685, 2303.16634, 2406.12624, 2404.13076, 2502.01534,
  2410.12784, 2410.02736, 2506.01937, 2504.15253, 2410.16184, 2305.17926, 2411.15594, 2506.07962,
  2505.22203 (v1 and v3), 2507.08794, 2407.18370, 2212.07919, 2311.08516, 2409.04168, 2503.21934,
  2410.13341, 2501.10970, 2406.18403, 2508.18076, 2508.06709, 2603.14732, 2604.06996, 2601.22548,
  2609.22512, 2510.11822, 2404.18796, 2602.20629, 2503.05061, 2504.03846, 2511.21140, 2502.04313,
  2410.21819, 2602.11570, 2608.20607, 2609.31857, 2410.20266, 2406.17859, 2605.29800.

## Included candidates (30)

### How these papers support the paragraph

1. **The paradigm and its biases.** Zheng 2023 and G-Eval established the paradigm. Ye 2025
   catalogues its biases.
2. **Self-preference.** Panickssery 2024 and Wataoka 2024 established it. Chen 2025 and Roytburg
   2026 qualify it: part of measured self-preference is judge error.
3. **Family bias.** Spiliopoulou 2025 measures family bias directly. Li 2026 (preference leakage)
   and Goel 2025 (similarity) show relatedness biases judges. Pombal 2026 shows self-preference
   survives objective, binary rubric items and is not removed by ensembling.
4. **Panels do not repair this.** Their errors are correlated (Kim 2025; Hossain 2026; JuryProbe),
   and judges are agreeable: they rarely reject invalid work (Jain 2025).
5. **Validation against people has to be done properly.** Use chance-corrected agreement
   (Thakur 2025), validate per task (Bavaresco 2025) and ground the judge in references
   (Krumdick 2026). Technical-domain evidence: physics marking (Yeadon 2026) and math proofs
   (QEDBench). Correct and report judge-based scores properly (Lee 2026; Dorner 2025; Li 2026
   sparse overlap).
6. **Where correctness can be computed, compute it.** JudgeBench shows judges near chance on
   objective pairs. Rule-based and model-based verifiers trade false negatives for false positives
   (Huang 2026), and LLM verifiers can be fooled (Zhao 2025). For routing only the residue to a
   judge, see Trust or Escalate and JuryProbe.
7. **Reference-based metrics.** ROSCOE, ReCEval and PRIME are in the process-supervision file (see
   "Overlaps" below), not this one.

### Must-cite (11)

**zheng2023mtbench**: Zheng et al., "Judging LLM-as-a-Judge with MT-Bench and Chatbot Arena",
NeurIPS 2023 Datasets and Benchmarks (arXiv 2306.05685, June 2023).
- **Measured:** how often strong LLM judges agree with human preferences. The paper also names
  position, verbosity and self-enhancement bias and the judges' limited reasoning ability.
- **Tasks:** the 80 multi-turn MT-Bench questions in eight categories (writing, roleplay,
  extraction, reasoning, math, coding, STEM and humanities knowledge), and Chatbot Arena
  crowd votes.
- **Headline:** GPT-4 as judge exceeds 80% agreement with human preferences, the same level as
  agreement between humans.
- **Technical coverage:** partial. There are math, coding and STEM categories, but the labels are
  preferences, not correctness.
- **Bearing on EngTrace:** this is the paradigm the May "AI Tribunal" inherited. Its own bias list,
  self-enhancement included, is the reviewers' objection. Agreement with chat preferences is not
  evidence of validity for step-level engineering correctness.

**panickssery2024selfrecognition**: Panickssery, Bowman, Feng, "LLM Evaluators Recognize and Favor
Their Own Generations", NeurIPS 2024 oral (OpenReview; the arXiv page states no venue) (2404.13076,
Apr 2024).
- **Measured:** self-recognition (can a model tell its own text from others') and self-preference,
  on summarisation (XSUM, CNN/DailyMail).
- **Headline:** GPT-4 and Llama 2 recognise their own outputs better than chance. Fine-tuning shows
  a linear correlation between self-recognition and the strength of self-preference, and the
  controlled experiments point to a causal link.
- **Technical coverage:** none.
- **Bearing on EngTrace:** the mechanism behind the circularity objection. A residual judge from a
  family outside the roster removes the self-recognition channel.

**thakur2025judgingjudges**: Thakur et al., "Judging the Judges: Evaluating Alignment and
Vulnerabilities in LLMs-as-Judges", GEM² Workshop 2025, ACL (the PDF carries the proceedings
header) (2406.12624, June 2024).
- **Measured:** how 13 judge models align with human judges when grading 9 exam-taker models on
  TriviaQA answers, a clean setting where humans agree with each other.
- **Headline:** only the best and largest judges align reasonably, and they remain well behind
  inter-human agreement. Their scores can differ from human scores by up to 5 points. Judges are
  lenient and sensitive to prompt complexity and length.
- **Metric caveat:** percent agreement can be high while scores differ greatly. The paper reports
  Scott's π; an earlier version used Cohen's κ.
- **Technical coverage:** none (factual QA).
- **Bearing on EngTrace:** answers the "thin human-alignment evidence" objection. Report
  chance-corrected agreement and per-class error rates against the 15 experts, not a single
  Spearman (the May paper's was 0.632).

**li2026preferenceleakage**: Li et al., "Preference Leakage: A Contamination Problem in
LLM-as-a-judge", ICLR 2026 (the PDF carries the ICLR 2026 header) (2502.01534, first posted
Feb 2025).
- **Measured:** bias that arises when the model generating synthetic training data and the judge
  are related, through being the same model, inheritance, or the same family.
- **Tasks:** Arena-Hard, AlpacaEval 2.0 and MT-Bench, scoring student models trained on the
  synthetic data.
- **Headline:** judges systematically favour related students. This "preference leakage" is
  harder to detect than previously known judge biases.
- **Technical coverage:** none.
- **Bearing on EngTrace:** relatedness through family is itself a contamination channel. This is
  the strongest published support for choosing a residual judge (Xiaomi MiMo) from outside every
  evaluated family. The brief called it "Li et al. 2025"; the venue year is 2026.

**spiliopoulou2025playfavorites**: Spiliopoulou et al., "Play Favorites: A Statistical Method to
Measure Self-Bias in LLM-as-a-Judge", arXiv 2508.06709 (Aug 2025; no venue found).
- **Measured:** self-bias and family bias, separated from genuine quality differences by a
  statistical model. Quality is anchored by an independent third-party judge, here human experts.
- **Data:** more than 5,000 prompt-completion pairs. Completions from nine LLMs (Llama 3, GPT,
  Mistral, Claude) on about 600 prompts, with expert ratings on six dimensions and verdicts from
  nine LLM judges.
- **Headline:** GPT-4o and Claude 3.5 Sonnet systematically score their own outputs higher, and
  also favour other models of their own family.
- **Technical coverage:** not domain-specific.
- **Bearing on EngTrace:** the reviewers' family-bias objection, measured. A GPT/Claude/Gemini
  panel scoring GPT/Claude/Gemini models is the confounded design this paper quantifies.

**goel2025greatmodels**: Goel et al., "Great Models Think Alike and this Undermines AI
Oversight", ICML 2025 spotlight (OpenReview; the arXiv page states no venue) (2502.04313,
Feb 2025).
- **Measured:** similarity between models via CAPA (chance-adjusted overlap of mistakes), and its
  effect on judge scores and on weak-to-strong training.
- **Tasks:** judges scored candidate models on 8,707 filtered MMLU-Pro questions across 14 domains.
- **Headline:** judge scores favour models similar to the judge (average Pearson r = 0.84 between
  score and similarity), which generalises self-preference. Models' mistakes grow more similar as
  capability rises.
- **Technical coverage:** yes (MMLU-Pro problem solving).
- **Bearing on EngTrace:** even a judge from another family can be biased through similarity, so
  family exclusion needs the expert-label check alongside it. It also explains why frontier
  panels are not independent.

**kim2025correlatederrors**: Kim, Garg, Peng, Garg, "Correlated Errors in Large Language Models",
ICML 2025 (arXiv comment and OpenReview) (2506.07962, June 2025).
- **Measured:** how correlated the errors of different LLMs are.
- **Data:** 349 LLMs on 12,032 MMLU questions (HuggingFace Open LLM Leaderboard), 71 LLMs on 14,042
  questions (HELM), and 20 LLMs on 1,800 résumé-job pairs.
- **Headline:** on one dataset, models agree 60% of the time when both err. Shared provider and
  architecture explain part of this, but larger, more accurate models are highly correlated even
  across providers. The paper shows downstream effects on LLM-as-judge evaluation and on hiring.
- **Technical coverage:** partial (MMLU includes STEM).
- **Bearing on EngTrace:** the PoLL rationale assumes diverse judges' errors cancel, and correlated
  errors break that assumption.

**tan2025judgebench**: Tan et al., "JudgeBench: A Benchmark for Evaluating LLM-based Judges",
ICLR 2025 (the PDF carries the ICLR 2025 header) (2410.12784, Oct 2024).
- **Measured:** judge accuracy on challenging response pairs whose labels reflect objective
  correctness.
- **Tasks:** knowledge, reasoning, math and coding.
- **Headline:** many strong models, GPT-4o among them, do only slightly better than random.
- **Technical coverage:** yes.
- **Bearing on EngTrace:** where correctness is objectively checkable, LLM judges are weak. That is
  the core argument for computing the answer check, milestone coverage and arithmetic
  deterministically, and leaving the judge only the remainder.

**krumdick2026nofreelabels**: Krumdick et al., "No Free Labels: Limitations of LLM-as-a-Judge
Without Human Grounding", COLM 2026 (the PDF carries the COLM 2026 header) (2503.05061, first
posted Mar 2025).
- **Measured:** how well automated graders, LLM judges among them, agree with experts' correctness
  labels.
- **Data:** BFF-Bench, 160 business and finance questions written by professionals, plus a hard
  MT-Bench subset. Experts labelled 1,200 responses.
- **Headline:** LLM judges beat the other graders. But without a correct reference, judges agree
  with the experts mainly on questions they can answer themselves. Expert-written references
  largely remove the problem.
- **Technical coverage:** partial (finance fundamentals, including quantitative questions).
- **Bearing on EngTrace:** the residual judge is given the gold milestone values from the
  template, and this paper is the evidence that such grounding makes verdicts trustworthy. It also
  implies that ungrounded Tribunal verdicts were capped by the judges' own engineering competence.

**pombal2026rubricselfpreference**: Pombal, Rei, Martins, "Self-Preference Bias in Rubric-Based
Evaluation of Large Language Models", COLM 2026 (the PDF carries the COLM 2026 header)
(2604.06996, Apr 2026).
- **Measured:** self-preference when judges issue binary verdicts on individual rubric criteria.
- **Tasks:** IFEval and LiveCodeBench, whose rubric items are programmatically verifiable, and
  HealthBench, whose rubrics are subjective.
- **Headline:** where the generator failed a rubric item, a judge can be more than 50% more likely
  to mark it satisfied when the output is its own. Ensembling judges mitigates this but does not
  eliminate it. On HealthBench, self-preference shifts scores by up to 10 points.
- **Technical coverage:** yes (code), plus medical.
- **Bearing on EngTrace:** the residual judge issues binary verdicts on milestones. Objective
  criteria do not immunise a judge, and ensembling (the Tribunal) is not a fix. Excluding the
  roster's families is warranted.

**huang2026verifierrobustness**: Huang et al., "From Accuracy to Robustness: A Study of Rule- and
Model-based Verifiers in Mathematical Reasoning", EMNLP 2026 (per the arXiv v3 comment of
2026-08-27) (2505.22203, May 2025; first titled "Pitfalls of Rule- and Model-based Verifiers -- A
Case Study on Mathematical Reasoning").
- **Measured:** how reliable rule-based and model-based answer verifiers are, both in static
  evaluation and inside RL training.
- **Headline:** rule-based verifiers miss equivalent answers written in other formats. These false
  negatives grow as the policy gets stronger. Model-based verifiers have higher static accuracy
  but are highly susceptible to reward hacking: they accept certain response patterns as correct,
  especially after fine-tuning.
- **Technical coverage:** yes (math).
- **Bearing on EngTrace:** frames the deterministic answer check honestly. A numeric check with
  relative tolerance cannot be hacked the way a model verifier can. Its known weakness, format
  false negatives, is narrow for numeric answers but should be stated.

### Should-cite (12)

**liu2023geval**: Liu et al., "G-Eval: NLG Evaluation using GPT-4 with Better Human Alignment",
EMNLP 2023 (2303.16634, Mar 2023).
- **Measured:** how well evaluation metrics correlate with human judgements.
- **Tasks:** summarisation and dialogue generation.
- **Headline:** G-Eval, which uses GPT-4 with chain-of-thought and form filling, reaches Spearman
  0.514 with humans on summarisation, well above previous methods. BLEU and ROUGE correlate
  weakly with humans. The paper also flags a possible bias of LLM evaluators toward LLM-written
  text.
- **Technical coverage:** none.
- **Bearing on EngTrace:** the classic citation for weak n-gram metrics and an early warning about
  LLM-favouring bias. It pairs with dropping ROUGE and BERTScore from the main table.

**wataoka2024selfpreference**: Wataoka, Takahashi, Ri, "Self-Preference Bias in LLM-as-a-Judge",
NeurIPS 2024 Safe Generative AI Workshop (rejected from ICLR 2025 per OpenReview) (2410.21819,
Oct 2024).
- **Measured:** self-preference, using a new fairness-style metric, on the Chatbot Arena dataset
  (33,000 dialogues).
- **Headline:** GPT-4 shows significant self-preference. Compared with human evaluators, LLMs give
  higher scores to lower-perplexity outputs whether or not they wrote them.
- **Technical coverage:** none.
- **Bearing on EngTrace:** a familiarity mechanism means same-family judges, which share data and
  style, can carry the bias even when they are not judging their own outputs.

**ye2025justiceprejudice**: Ye et al., "Justice or Prejudice? Quantifying Biases in
LLM-as-a-Judge", ICLR 2025 (OpenReview) (2410.02736, Oct 2024).
- **Measured:** 12 potential judge biases (position, verbosity, authority, self-enhancement and
  others), probed with the automated CALM framework.
- **Data:** the "fact-related" datasets include GSM8K, MATH and ScienceQA.
- **Headline:** significant biases persist on specific tasks even when overall judge performance
  is strong.
- **Technical coverage:** partial (math and science QA).
- **Bearing on EngTrace:** the single citation for the standard bias taxonomy.

**gu2024surveyjudge**: Gu et al., "A Survey on LLM-as-a-Judge", arXiv 2411.15594 (Nov 2024; CoRR
only).
- **Covers:** how to build reliable LLM-judge systems: consistency, bias mitigation, adaptation to
  different tasks, and how to evaluate the reliability of the judges themselves, with a benchmark
  for that.
- **Bearing on EngTrace:** the general pointer, standing in for three surveys. The alternatives are
  listed among the rejected candidates.

**chen2025selfpreferencereason**: Chen et al., "Do LLM Evaluators Prefer Themselves for a
Reason?", arXiv 2504.03846 (Apr 2025; rejected from ICLR 2026 per OpenReview).
- **Measured:** harmful self-preference (favouring objectively worse responses) separated from
  legitimate self-preference (favouring genuinely better ones). Verifiable math, factual and code
  benchmarks make the separation possible.
- **Scope:** seven model families, plus an LMArena extension.
- **Headline:** stronger models prefer themselves more, mostly legitimately. Harmful self-preference
  persists when the evaluator erred as a generator, and is stronger in stronger models. Long
  chain-of-thought before judging reduces it.
- **Technical coverage:** yes.
- **Bearing on EngTrace:** in technical tasks, harmful self-preference concentrates exactly where
  the judge is wrong itself, which is where error detection matters.

**roytburg2026narcissists**: Roytburg et al., "Are LLM Evaluators Really Narcissists? Sanity
Checking Self-Preference Evaluations", ICML 2026 (arXiv comment and OpenReview) (2601.22548,
Jan 2026).
- **Measured:** whether published self-preference findings survive an evaluator-quality baseline.
  The baseline compares a judge's voting when it evaluates itself with its voting when it
  evaluates another model.
- **Headline:** only 51% of the examples in earlier findings stay statistically significant,
  though these cover 89.6% of the total self-preference probability mass. An entropy analysis
  points to uncertainty-driven overlap.
- **Bearing on EngTrace:** the balanced citation. Part of measured self-preference is judge error,
  not narcissism. EngTrace avoids having to settle the question by design (no judge from a roster
  family) and checks the residual bias against expert labels.

**hossain2026agreementoverstates**: Hossain, Yousefi, Lim, "Agreement Overstates Evidence: Error
Dependence in LLM Judge Consensus", arXiv 2609.22512 (Sep 2026).
- **Measured:** how dependent judges' errors are, and what that does to consensus.
- **Data:** 1,195 stratified RewardBench v2 items (factuality, focus, reasoning, safety), plus
  RewardBench, UltraFeedback and PKU-SafeRLHF.
- **Headline:** in a bank of ten judges, the mean pairwise error correlation is 0.21, so the ten
  carry about as much information as 3.5 independent judges. Frontier judges are more dependent:
  mean error correlation 0.56, and they make the same error together 7.7 times as often as
  independence predicts, across providers too. Ignoring the dependence flips significance in up
  to 28% of system comparisons.
- **Technical coverage:** partial (a reasoning subset).
- **Bearing on EngTrace:** a three-judge frontier Tribunal gives far fewer than three independent
  votes, and its consensus is not evidence of correctness.

**jain2025agreeableness**: Jain et al., "Beyond Consensus: Mitigating the Agreeableness Bias in
LLM Judge Evaluations", arXiv 2510.11822 (Oct 2025; rejected from ICLR 2026 per OpenReview).
- **Measured:** judges' true-positive and true-negative rates, ensembles of judges, a minority-veto
  rule and a regression calibrated on human labels.
- **Task:** feedback on 366 high-school Python programs, judged by 14 LLMs.
- **Headline:** judges accept valid outputs well (TPR 96%) but reject few invalid ones (TNR under
  25%). With class imbalance, this inflates apparent reliability. Majority voting is not enough;
  the calibrated regression brings the maximum absolute error down to 1.2%.
- **Technical coverage:** yes (code).
- **Bearing on EngTrace:** the same pattern as EngTrace's planted-defect test, where judges rarely
  flag a flawed step. Report sensitivity and specificity separately, and do not lean on majority
  votes.

**bavaresco2025llmsinsteadofhumans**: Bavaresco et al., "LLMs instead of Human Judges? A Large
Scale Empirical Study across 20 NLP Evaluation Tasks", ACL 2025 (2406.18403, June 2024).
- **Measured:** how well 11 LLMs reproduce human annotations across JUDGE-BENCH, 20 NLP datasets.
- **Headline:** variance across models and datasets is substantial. Reliability depends on the
  property being judged, the expertise of the human annotators, and whether the text was written
  by people or models. The authors conclude LLMs must be validated against humans before being
  used as evaluators.
- **Technical coverage:** none.
- **Bearing on EngTrace:** makes expert validation a task-specific precondition, not an optional
  extra.

**yeadon2026physicsjudge**: Yeadon, Hardy, Mackay, Agra, "LLM-as-a-judge validity is strongly
task-dependent across physics assessment formats", arXiv 2603.14732 (Mar 2026; v3 dated
2026-09-29).
- **Measured:** LLM marking against human markers in physics, across three formats (structured
  questions, essays, scientific plots).
- **Conditions:** blind, solution-provided, false-solution and anchored.
- **Judges:** GPT-5.2, Grok 4.1, Claude Opus 4.5, DeepSeek-V3.2 and Gemini 3 Pro, plus committees
  of them.
- **Headline:**
  - Structured questions (n = 771 and n = 1,151): rank agreement with human markers above
    Spearman 0.6. Official solutions reduce errors.
  - False solutions: absolute accuracy drops, because models defer to the reference, but rank
    order survives.
  - Essays (275): rank agreement is near zero.
  - Plot elements (1,400): ρ > 0.83.
  - The reliability of the human markers bounds what can be claimed about AI-human agreement.
- **Technical coverage:** yes.
- **Bearing on EngTrace:** the closest technical-domain analogue. A Spearman of about 0.6 is what
  frontier judges reach on structured physics, so the May paper's 0.632 is ordinary rather than
  strong evidence. References matter (EngTrace's gold milestones), and inter-expert reliability
  should be reported as the ceiling.

**gonzalez2026qedbench**: Gonzalez et al., "QEDBench: Quantifying the Alignment Gap in Automated
Evaluation of University-Level Mathematical Proofs", ICML 2026 (OpenReview venue id
ICML.cc/2026/Conference; arXiv v3 states no venue) (2602.20629, Feb 2026).
- **Measured:** how closely LLM judges match human experts grading upper-undergraduate and early
  graduate proofs. Two rubrics: course-specific, and expert common-knowledge criteria.
- **Setup:** 7 judges by 5 solvers, against more than 1,000 hours of human evaluation.
- **Headline:** Claude Opus 4.5, DeepSeek-V3, Qwen 2.5 Max and Llama 4 Maverick inflate mean scores
  by up to +0.18, +0.20, +0.30 and +0.36 respectively. Several solvers degrade on discrete math and
  graph theory.
- **Technical coverage:** yes.
- **Bearing on EngTrace:** expert-anchored evidence of judge leniency in a technical domain.

**lee2026reportjudge**: Lee et al., "How to Correctly Report LLM-as-a-Judge Evaluations",
ICML 2026 (the PDF carries the PMLR 306 header) (2511.21140, Nov 2025).
- **Measured:** the bias that a judge's imperfect sensitivity and specificity put into naive
  judge-based accuracy.
- **Method:** a plug-in correction with confidence intervals that cover both test-set and
  calibration-set uncertainty, adaptive allocation of calibration samples, and unbiasedness under
  distribution shift.
- **Bearing on EngTrace:** the expert-labelled traces can serve as the calibration set to correct
  and bound whatever the residual judge contributes. Pairs with the plan to report the judged
  fraction per model.

### Optional (7)

**zhao2025onetoken**: Zhao et al., "One Token to Fool LLM-as-a-Judge", arXiv 2507.08794 (Jul 2025;
rejected from ICLR 2026 per OpenReview). In reference-based RLVR, generative reward models give
false-positive verdicts to "master keys": non-word symbols such as ":" or "." and openers such as
"Thought process:". This affects GPT-o1 and Claude-4, and adversarially trained "Master-RMs"
mitigate it. Technical: yes (reasoning verification). Bearing: an LLM checking an answer against a
reference can be exploited; a numeric check cannot.

**jung2025trustescalate**: Jung, Brahman, Choi, "Trust or Escalate: LLM Judges with Provable
Guarantees for Human Agreement", ICLR 2025 oral (2407.18370, Jul 2024). The judge evaluates
selectively, with a provable guarantee of human agreement: confidence comes from Simulated
Annotators, and Cascaded Selective Evaluation escalates from cheap to strong judges only when
needed. Starting from Mistral-7B, it guarantees over 80% human agreement on Chatbot Arena at high
coverage. Bearing: a principled precedent for EngTrace's routing, where the cheap, confident tier
is deterministic.

**dorner2025limits**: Dorner, Nastl, Hardt, "Limits to scalable evaluation at the frontier: LLM as
Judge won't beat twice the data", ICLR 2025 oral (2410.13341, Oct 2024). Theorem: when the judge
is no more accurate than the model it evaluates, no debiasing method can cut the ground-truth
labels needed by more than half. In practice the savings are smaller still. Bearing: at the
frontier (judges scoring GPT-5-class models), expert labels cannot be replaced. EngTrace keeps them
and minimises how much it relies on a judge.

**zhou2026juryprobe**: Zhou, Lin, "JuryProbe: An Empirical Consensus-Risk Diagnostic for Routing
Reference-Free Factuality Judge Panels to Grounded Verification", TMLR 2026 (2608.20607, Aug 2026).
On audited FEVER corruptions, reference-free panels share false negatives: false-negative-only
correlations are 0.402 and 0.368, and false-consensus lifts are 3.13x and 18.13x. In a best-case
diagnostic, routing majority accepts to the same judges with trusted references removes unanimous
false consensus. The paper offers no formal guarantee. Bearing: agreement among ungrounded judges
is not evidence; ground the check.

**li2026sparseoverlap**: Li, Mukherjee, Pal, "LLM Judge Validation Under Sparse Overlap: From
Inference to Design", NeurIPS 2026 (per arXiv comment) (2609.31857, Sep 2026). When few items carry
labels from more than one annotator, judge-validation decisions go wrong often. At 5% pairwise
overlap, wrong deployment decisions reach 25%, and the wrong best judge among ten is chosen 65% of
the time. The paper gives a minimum-overlap rule (0.25 suffices for judges that are not borderline)
and a stratified allocation that halves false rejections. Bearing: tells EngTrace how much
double-labelling across its 15 experts it needs to support claims about inter-expert reliability
and judge validity.

**szymanski2025expertknowledge**: Szymanski et al., "Limitations of the LLM-as-a-Judge Approach
for Evaluating LLM Outputs in Expert Knowledge Tasks", IUI 2025 (2410.20266, Oct 2024). Subject
experts (dietitians, clinical psychologists) and LLM judges made pairwise comparisons on domain
tasks. They agreed 68% of the time in dietetics and 64% in mental health, with variation by aspect.
Technical: expert domains, but not numeric. Bearing: the standard citation for keeping experts in
the loop on expert-knowledge tasks.

**zhou2025jetts**: Zhou et al., "Evaluating Judges as Evaluators: The JETTS Benchmark of
LLM-as-Judges as Test-Time Scaling Evaluators", ICML 2025 (2504.15253, Apr 2025). Ten judges on
math, code and instruction following, in three settings: response reranking, step-level beam
search and critique-driven refinement. Judges match outcome reward models at reranking but trail
process reward models in beam search, and their critiques do not help refinement. Technical: yes.
Bearing: judges are weaker at step-level than outcome-level judgement, so EngTrace's step-level
milestones should not depend on a judge wherever a deterministic check exists.

## Overlaps with other topic files (checked read-only; not added here to avoid duplicate keys)

- **golovneva2022roscoe** (process_supervision), ROSCOE (2212.07919): **published at ICLR 2023**
  (the PDF carries the "Published as a conference paper at ICLR 2023" header). This is the
  citation for reference-based metrics failing on reasoning. It compares its metrics with ROUGE-1/2/L,
  BLEURT, PRISM, BERTScore, BARTScore and CTC on five human-annotated and six perturbed diagnostic
  datasets, and reports that its own metrics consistently outperform those baselines. Its Appendix
  Table 19 repeats the comparison in the reference-based setting, against gold reasoning steps,
  with the same conclusion.
- **prasad2023receval** (process_supervision), ReCEval (2304.10703): the most direct numbers. On
  EntailmentBank-challenge error types, Somers' D for ROUGE-2 is -0.01, -0.02 and 0.14
  (hallucination, negation, swap), and for BERTScore 0.09, 0.02 and 0.07. ReCEval-correctness
  scores 0.89, 0.88 and 0.39. In that setup the baselines compare the chain with the input
  context, not with a gold chain.
- **wang2026prime** (process_supervision), PRIME (2602.11570): ACL 2026 main conference (OpenReview
  imports this as "ACL (1) 2026"). On 2,530 college-level math and engineering problems, verifiers
  frequently miss derivation flaws behind correct answers. The judge paragraph should cite it too.
- **tyen2023bigbenchmistake** (process_supervision), 2311.08516, ACL Findings 2024: LLMs cannot
  find reasoning errors but can correct them when told where the error is.
- **guo2026chemcotbenchv2** (process_supervision), 2606.03660: deterministic chemistry-rule checks
  of reasoning steps "rather than another LLM judge". The closest design precedent for EngTrace's
  deterministic step checks.
- **ansari2026physicsregrading** (physics_science and symbolic_contamination), 2609.13009: expert
  re-grading of physics evaluations. Relevant to grading validity; owned by those topics.

## Considered and not included (one line each)

- Wang et al., "Large Language Models are not Fair Evaluators" (2305.17926; ACL 2024 per Semantic Scholar): position bias; covered by MT-Bench and CALM, and the residual judge is pointwise.
- Lambert et al., RewardBench (2403.13787; NAACL Findings 2025 per OpenReview/DBLP): reward-model benchmark, not judge validity in technical domains.
- Malik et al., RewardBench 2 (2506.01937; ICLR 2026, verified): as above; JudgeBench is the better judge benchmark for objective correctness.
- Liu et al., RM-Bench (2410.16184; ICLR 2025, verified): style bias in reward models; peripheral.
- Li et al., "From Generation to Judgment" (2411.16594; EMNLP 2025 per Semantic Scholar): second survey; Gu et al. suffices.
- Li et al., "LLMs-as-Judges: A Comprehensive Survey" (2412.05579): third survey; redundant.
- Stephan et al., "From Calculation to Adjudication" (2409.04168; arXiv only): math judges track the candidate model's strength and 70-75% of verdicts are predictable from surface features; overlaps JudgeBench and Chen 2025. Good alternate if a math-specific judge study is wanted.
- Petrov et al., "Proof or Bluff?" (2503.21934): the abstract covers expert grading of USAMO proofs only, not LLM graders; QEDBench covers judges against experts.
- Calderon et al., "The Alternative Annotator Test" (2501.10970; ACL 2025 per Semantic Scholar): test for replacing annotators with an LLM; EngTrace does not replace experts. Alternate to Lee 2026.
- Chehbouni et al., "Neither Valid nor Reliable?" (2508.18076; NeurIPS 2025 per Semantic Scholar/DBLP): measurement-theory position paper; good framing; alternate. Same group's "LLJ Cards" (2609.24516) is a reporting checklist.
- Kortemeyer, Nöhl, Onishchuk, "Grading Assistance for a Handwritten Thermodynamics Exam using AI" (2406.17859; Phys. Rev. Phys. Educ. Res. 20, 020144, 2024; verified): engineering-thermodynamics exam grading by AI. Part-wise grading was more reliable, and failing exams still needed human graders. **Recommended swap-in** for an optional slot; left out only because the 30-entry cap was reached after downloads.
- Kohli, "Nine Judges, Two Effective Votes" (2605.29800; verified): 9 judges from 7 families equal about 2 independent votes on NLI. Same point as Hossain 2026 and Kim 2025; single author.
- Acharya et al., RoPoLL (2606.30931; ICML 2026 AgenticUQ workshop): robust aggregation for PoLL. Better aggregation does not remove shared bias.
- Wang et al., "Reasoning Jury" (2608.12585): deliberating open-weight juries beat single frontier judges at finding reasoning defects. Industry preprint; no expert validation in the abstract.
- Liu, "LLMs as a Jury: Cross-Model Consensus Can Outperform PRMs" (2607.10139): answer selection by consensus, not judging. Its shared-error floor is near zero on math, non-trivial on science.
- Kranti and Vajjala, "LLM Judges Can Be Too Generous When There Is No Reference Answer" (2607.12885): same finding as Krumdick 2026, which is peer-reviewed.
- Wahi, "LLM-as-a-Judge Is Not an Oracle: Why Self-Improving Agents Need Deterministic Guardrails" (2609.02246): supportive but an anecdotal single-author experience report.
- Norman et al., "Reliability without Validity" (2606.19544): 21 judges; exact match overstates agreement by 33-41 points relative to kappa. Strong alternate to Thakur.
- Yosef et al., "Rethinking Math Reasoning Evaluation: ... Beyond Symbolic Rigidity" (2604.22597): counterpoint that LLM answer checks beat symbolic comparison on varied answer formats. Covered by Huang 2026; EngTrace answers are numeric.
- Huang et al., "Reasoning Model Is Superior LLM-Judge, Yet Suffers from Biases" (2601.03630; ACL 2026 EvalEval workshop): reasoning-model judges remain biased; peripheral.
- Zhang et al., "Can We Trust LLM Judges: Capability-Dependent Biases ..." (2609.12002): more capable examinees get more lenient verdicts; preprint; alternate.
- Yang et al., "Quantifying and Mitigating Self-Preference Bias of LLM Judges" (2604.22891): equal-quality pairs; overlaps Pombal and Spiliopoulou.
- Guey and Bougault, "Self-Preference Is Weak or Absent in Verifiable Instruction-Following Revision" (2606.20093): small null result about revising, not judging.
- Lehr et al., "Extreme Self-Preference in Language Models" (2509.26464): self-preference in word association and hiring, not evaluation.
- Chen et al., "Beyond the Surface: Measuring Self-Preference" (2506.02592): DBG score against gold judgements; overlaps Spiliopoulou.
- "Breaking the Mirror" (2509.03647) and St. Amand et al., "Self-Generated Text Recognition" (2608.26159): mitigation and measurement detail; not needed.
- "Know Thyself? On the Incapability and Implications of AI Self-Recognition" (2510.03399): self-recognition counter-finding; peripheral.
- Feng et al., Sage (2512.16041; rejected from ICLR 2026): judge consistency without human labels; less relevant.
- Hwang et al., "Can You Trick the Grader?" (2508.07805): persuasive wording inflates math-judge scores of incorrect solutions by up to 8%. Technical and relevant, but adversarial; alternate.
- Findeis et al., "Can External Validation Tools Improve Annotation Quality for LLM-as-a-Judge?" (2507.17015; ACL 2025): code execution and web search help AI annotators on math and code; supports grounding; alternate.
- VerifyBench (2507.09884), xVerify (2504.10481), "Scaling Generative Verifiers ..." (2511.13027): RLVR verifier tooling, not validity studies.
- Puduppully et al., "Correct Answers, Invalid Traces" (2609.38107): programmatic trace checks on iGSM; on the hardest instances 31.6% of correct answers have invalid traces. Belongs to the process-supervision topic; flagged to the lead.
- Song et al., "Beyond the Illusion of Consensus" (2603.11027): model-level Spearman 0.99 hides sample-level r = 0.72; supports Hossain; alternate.
- Chen et al., "A Judge Should Know What Changed: Construct Validity ..." (2608.24419; under TMLR review): judges are invariant (S = 0.945) but insensitive to construct changes (R = 0.319); alternate.
- Zhu, "Three Ways Classical Test Theory Can Mislead About LLM Judges" (2609.29709); Usami et al., "LLM Judges Have Dark Current" (2606.15610): psychometric method papers.
- Li et al., "How Many Humans Are 32 LLM Judges Worth?" (2609.21277): ChaosNLI effective panel size; same point as Hossain.
- Fox et al., "LLM Judges Verify Presence, Not Absence" (2608.31016): judges miss clinical omissions, and per-fact checking recovers detection. An interesting analogue to milestone decomposition, but clinical.
- Gozel, "Commit-first LLM judging inherits the judge's own errors" (2609.00088): same mechanism as Krumdick; single author.
- Zhang et al., "When the Judge Should Not Decide" (2608.07813): judges inside reasoning pipelines; tangential.
- AdvancedMathBench (2607.11849): best proof verifier Balanced F1 65.1. Alternate to QEDBench, which has stronger expert anchoring and an ICML venue.
- Naik et al., "Do We Need Frontier Models to Verify Mathematical Proofs?" (2604.02450); Grayzel, "Cost-Effective Automated Judging of NL Math Proofs" (2608.00004): cost-focused.
- Sinhahajari et al., "On the Limits of LLM-as-Judge for Scientific Novelty Assessment" (2606.12071): novelty, not correctness.
- He et al., "Judging the Judges: Human Validation of Multi-LLM Evaluation for K-12 Science Materials" (2602.13243): two education experts; qualitative.
- Zhang, "Testing Frontier LLMs' Physics Literacy in Parallel Physical Worlds" (2607.00276): notes judge reliability does not transfer across frameworks; single author.
- Yagubyan, "The Coin Flip Judge?" (2606.13685): flip rates for two OpenAI judges only.
- Soumik, "Judging the Judges: ... Bias Mitigation Strategies" (2604.23178; TMLR 2026): compares debiasing strategies; not needed.
- Feng et al., "Noisy but Valid" (2601.20913; ICLR 2026): hypothesis testing with imperfect judges; alternate to Lee 2026.
- Dependence-Aware Label Aggregation via Ising Models (2601.22336; ICML 2026 per OpenReview): aggregation under dependent judge errors; alternate for the panel discussion.
- Kuai et al., "A Statistical Framework for Auditing Behavioral Dependence and Induced Bias in LLM Judges" (2604.07650; q18): 18 LLMs from six families; entanglement between models predicts judges over-endorsing on MMLU-Pro (rho 0.51-0.52) and MATH-500 (rho 0.44-0.46). Technical, but the same point as Goel 2025; alternate.
- Wang et al., "Chain-of-Models: Cross-Model Auditing for Bias-Robust LLM Judges" (2607.28636; q18; marked work in progress): whether a same-model, same-family or other-family auditor works best depends on the bias type; preliminary.
- Sunkavalli, "Anchor-Judge Error Correlation ..." (2609.08826; q18): estimator for panel common-mode error, with a family-block test; single author, validated mostly in simulation.
- Shi et al., position-bias study (AACL 2025); Chen et al., "Humans or LLMs as the Judge?" (EMNLP 2024); Koo et al., CoBBLEr: bias studies covered by CALM.
- Guerdan et al., "Validating LLM-as-a-Judge Systems under Rating Indeterminacy" (NeurIPS 2025 per OpenReview): relevant to human label variation (Plank 2022 is cited). Not opened, so not verified beyond the listing.
- SciArena (2507.01001) and HealthBench (2505.08775): scientific-literature and medical judge meta-evaluations; not opened. QEDBench and Yeadon cover technical-domain expert anchoring.
- "Can LLMs Replace Human Evaluators? An Empirical Study of LLM-as-a-Judge in Software Engineering" (PACMSE/FSE 2025, DOI 10.1145/3728963; from OpenAlex): an SE-domain judge study; no arXiv ID checked; alternate.

## Updates to existing citations

- **Verga et al. 2024, "Replacing Judges with Juries: Evaluating LLM Generations with a Panel of
  Diverse Models"** (`verga2024poll`, arXiv 2404.18796): **no peer-reviewed venue found as of
  2026-10-01.**
  - The arXiv page has v1 (29 Apr 2024) and v2 (1 May 2024), with no comments or journal-ref.
  - Semantic Scholar's venue is arXiv.org, and its DBLP key is the CoRR one
    (journals/corr/abs-2404-18796).
  - OpenReview lists only "CoRR 2024".
  - Cite it as an arXiv preprint (Verga et al., 2024).
  - Its central assumption, that a panel of diverse models reduces bias because their errors are
    independent, is now directly contested by Kim et al. 2025 (ICML), Goel et al. 2025 (ICML),
    Hossain et al. 2026 and Pombal et al. 2026 (ensembling reduces but does not remove
    self-preference). A 2026 robust-aggregation successor, RoPoLL (2606.30931), formalises PoLL's
    unbounded bias under contamination.
- Papers the brief named with a different year than their venue year: Li et al. "Preference
  Leakage" is ICLR 2026 (first posted Feb 2025). Thakur et al. "Judging the Judges" is GEM² 2025
  (first posted June 2024). Keys and the `year` field follow the cited set's venue-year
  convention; `date` gives the first-posted month.

## What I could not verify

- **WebSearch:** none of the 20+ queries the brief asked for could run (session budget exhausted).
  Recall rests on the arXiv API, OpenAlex and targeted checks of known titles and IDs, so
  non-arXiv venue papers (ACM, IEEE, journals) are under-sampled.
- **arXiv queries:** q15 (BERTScore AND reasoning AND "human judgments") returned only 7 records,
  none relevant. q16 (judge AND physics/engineering/chemistry) was swamped by software-engineering
  agent papers: "engineering" matched them. Neither surfaced a judge-versus-expert study in
  engineering; the closest are Yeadon 2026 (physics) and Kortemeyer 2024 (thermodynamics exam,
  found via OpenAlex). q18 ("family bias") confirmed the selection. Its new hits are listed among
  the rejected candidates.
- **Venues taken from one source only:**
  - Huang et al. as EMNLP 2026: arXiv v3 comment only; the PDF has no header, and OpenReview
    shows an ICLR 2026 withdrawal.
  - QEDBench as ICML 2026: OpenReview venue id only.
  - Li et al. (sparse overlap) as NeurIPS 2026: arXiv comment only.
  - Goel et al. as ICML 2025 spotlight and Panickssery et al. as NeurIPS 2024 oral: OpenReview and
    Semantic Scholar, not the arXiv pages.
- **No venue found:** Spiliopoulou 2025, Chen 2025, Jain 2025, Zhao 2025, Gu 2024, Hossain 2026
  and Yeadon 2026. Chen, Jain and Zhao appear on OpenReview as ICLR 2026 rejections.
- **Within-paper details:** taken from abstracts and spot-checks of the extracted text. Full
  results sections were not read.
- **DBLP:** unusable (anti-bot challenge), so DBLP keys were seen only through Semantic Scholar's
  external IDs.
