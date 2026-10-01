# Cited set: resolution and acquisition

2026-10-01. Every entry of `papers.json` (the May 2026 submission's cited set) was resolved to a
source, every arXiv id was checked against its abstract page, the current publication status of
every paper was recorded in a new field `venue_now` (with `arxiv_latest`, the date of the latest
arXiv version), and every paper was fetched with a text layer for the review panel.

## Result

| | |
|---|---|
| Entries in `papers.json` | 46: the file held 45, and `skelic2025circuit` was added (problems 1 and 2) |
| Resolved to a source | 46: 40 arXiv ids, 6 publisher or ACL Anthology PDFs |
| The six RESOLVE entries | all resolved: `chen2025apbench`, `gulati2024putnamaxiom`, `xie2023storytelling`, `xie2023deltascore`, `xie2025finchain`, `felten2025engibench` |
| arXiv ids corrected | 1: `liu2025circuit`, 2502.07980 to 2511.18221 |
| VERIFY entries | `deng2024contamination` and `ott2022saturation` confirmed as catalogued; `liu2025circuit` corrected |
| Fetched, with a text layer | 46 (44 in this run; `lin2004rouge` and `mirzadeh2024gsmsymbolic` were already on disk from the 16:40Z run) |
| Failed | 0 |

## How to re-check

    python docs/related_work_oct2026/check_arxiv_abs.py --id 2203.04592 --id 2211.02570 --id 2301.09790 --id 2303.08991
    python docs/related_work_oct2026/fetch_papers.py --file papers.json
    python docs/related_work_oct2026/fetch_papers.py --verify
    python docs/related_work_oct2026/check_text.py

- `check_arxiv_abs.py` fetches `https://arxiv.org/abs/<id>` for the 40 entries with an id, plus
  the arXiv preprints of four URL-sourced papers (Ott 2203.04592, Plank 2211.02570, Xie 2023a
  2301.09790, Xie 2023b 2303.08991). It compares the page title with `title` and the latest
  version's date with `arxiv_latest`, and writes `notes/arxiv_abs_check.json` (titles, authors,
  comments, journal-ref, every version and its date). Run of 2026-10-01: 39 titles match, 1 is a
  variant (`gulati2024putnamaxiom`, problem 3), 0 mismatch; all 40 `arxiv_latest` values match.
- `fetch_papers.py --verify` re-hashed every file in `MANIFEST.json` at 17:25Z: 161 files (the 46
  cited plus 115 candidate records written by other sessions), 0 problems.
- `check_text.py` looks in the first two pages of each text file for the entry's title and its
  first author: 46 of 46 show the full title and the first author, 0 flagged.
- `venue_now` comes from the arXiv comments and journal-ref (in the abs-page record), OpenReview
  venue records (including OpenReview's imports from DBLP), Crossref (the ACL Anthology, IEEE,
  Springer and Nature DOIs), Semantic Scholar, the PMLR volume 267 index, the ACL Anthology page
  of ROUGE (its workshop volume) and OpenAlex (GLUE at ICLR 2019). Every venue rests on at least
  two of these, except FEABench's workshops, which rest on its arXiv comment alone.
  "arXiv only as of 2026-10" means Semantic Scholar lists only the arXiv record, OpenReview shows
  no published version, and the arXiv page names no venue; Crossref was also searched for
  ElecBench, TransportBench, PoLL, ABench-Physics, BIGbench and CIRCUIT, with the same result.
- `authors` now holds each paper's author list as its source gives it (the arXiv page for the 40
  arXiv entries; the PDF, Crossref or the Anthology for the other six), cut at 25 names with the
  total for longer lists. The reference list's own lists are not kept where they differ
  (problem 6). The brief's rule, keep the reference list's list when it is more complete, never
  applied: the only printed list longer than the paper's (Lightman et al., 30 names for 10
  authors) is longer because 23 of its names are not authors.

## Table

Pages and text characters are from `MANIFEST.json`; "source" is what `fetch_papers.py` downloads.

| key | cited_as | source | venue_now | correction | status | pages | text chars |
|---|---|---|---|---|---|---:|---:|
| austin2021mbpp | Austin et al., 2021 | arXiv:2108.07732 | arXiv only as of 2026-10 | - | fetched | 34 | 126,717 |
| chen2021humaneval | Chen et al., 2021 | arXiv:2107.03374 | arXiv only as of 2026-10 | - | fetched | 35 | 155,072 |
| chen2025apbench | Chen et al., 2025 | https://www.nature.com/articles/s41598-025-91150-5.pdf | Scientific Reports 15, article 7944 (2025) | resolved (RESOLVE): Nature open-access PDF; author list replaced (the reference's is wrong) | fetched | 15 | 65,541 |
| cheng2024elecbench | Cheng et al., 2024 | arXiv:2407.05365 | arXiv only as of 2026-10 | - | fetched | 26 | 105,543 |
| cobbe2021gsm8k | Cobbe et al., 2021 | arXiv:2110.14168 | arXiv only as of 2026-10 | - | fetched | 22 | 45,479 |
| deng2024contamination | Deng et al., 2024 | arXiv:2311.09783 | NAACL 2024 (Long Papers), pp. 8706-8719 | VERIFY done: id is the NAACL 2024 paper | fetched | 14 | 62,027 |
| du2025supergpqa | Du et al., 2025 | arXiv:2502.14739 | NeurIPS 2025 Datasets and Benchmarks Track (poster) | - | fetched | 249 | 1,144,089 |
| fan2024hardmath | Fan et al., 2024 | arXiv:2410.09988 | ICLR 2025 (poster); earlier, NeurIPS 2024 MATH-AI workshop | - | fetched | 32 | 81,607 |
| gulati2024putnamaxiom | Gulati et al., 2024 | arXiv:2508.08292 | ICML 2025 (PMLR 267:20723-20747); the cited version is the NeurIPS 2024 MATH-AI workshop paper | resolved (RESOLVE): arXiv 2508.08292, the ICML 2025 version (title adds "in LLMs") | fetched | 27 | 82,636 |
| guo2025engdesign | Guo et al., 2025b | arXiv:2509.16204 | NeurIPS 2025 Datasets and Benchmarks Track (poster) | - | fetched | 61 | 175,343 |
| heesch2025realworld | Heesch et al., 2025 | arXiv:2505.13484 | AI 2025: Advances in Artificial Intelligence (Australasian Joint Conference on AI, per DBLP), Springer LNCS, pp. 54-66 | - | fetched | 23 | 57,769 |
| hendrycks2021mmlu | Hendrycks et al., 2020 / 2021a | arXiv:2009.03300 | ICLR 2021 | - | fetched | 27 | 83,915 |
| hendrycks2021math | Hendrycks et al., 2021b | arXiv:2103.03874 | NeurIPS 2021 Datasets and Benchmarks Track | - | fetched | 22 | 85,222 |
| jimenez2024swebench | Jimenez et al., 2024 | arXiv:2310.06770 | ICLR 2024 (oral) | - | fetched | 52 | 154,063 |
| kazemi2025bbeh | Kazemi et al., 2025 | arXiv:2502.19187 | ACL 2025 (Long Papers), pp. 26473-26501 | - | fetched | 36 | 128,782 |
| li2025eeebench | Li et al., 2025 | arXiv:2411.01492 | CVPR 2025, pp. 13337-13349 | - | fetched | 52 | 164,129 |
| lightman2023verify | Lightman et al., 2023 | arXiv:2305.20050 | ICLR 2024 (poster) | - | fetched | 29 | 53,786 |
| lin2004rouge | Lin, 2004 | https://aclanthology.org/W04-1013.pdf | Text Summarization Branches Out (ACL 2004 workshop), pp. 74-81 | - | fetched | 8 | 38,534 |
| liu2025circuit | Liu et al., 2025 | arXiv:2511.18221 | IEEE Frontiers in Education Conference (FIE 2025), pp. 1-9 | **id corrected** 2502.07980 -> 2511.18221 (the printed title's paper) | fetched | 9 | 47,412 |
| luo2024bigbench_t2i | Luo et al., 2024 | arXiv:2407.15240 | arXiv only as of 2026-10 | - | fetched | 23 | 82,629 |
| mirzadeh2024gsmsymbolic | Mirzadeh et al., 2024 | arXiv:2410.05229 | ICLR 2025 (poster) | - | fetched | 24 | 88,776 |
| mudur2025feabench | Mudur et al., 2025 | arXiv:2504.06260 | NeurIPS 2024 workshops (MATH-AI; Open-World Agents) per the arXiv comment; no archival venue found | - | fetched | 39 | 122,235 |
| ott2022saturation | Ott et al., 2022 | https://www.nature.com/articles/s41467-022-34591-0.pdf | Nature Communications 13, article 6793 (2022) | VERIFY done: DOI and PDF URL correct | fetched | 11 | 58,556 |
| plank2022hlv | Plank, 2022 | https://aclanthology.org/2022.emnlp-main.731.pdf | EMNLP 2022 (main), pp. 10671-10682 | - | fetched | 12 | 55,630 |
| qiu2025phybench | Qiu et al., 2025 | arXiv:2504.16074 | NeurIPS 2025 Datasets and Benchmarks Track (poster) | - | fetched | 34 | 98,849 |
| shojaee2025llmsrbench | Shojaee et al., 2025 | arXiv:2504.10415 | ICML 2025 (oral) | - | fetched | 35 | 114,925 |
| syed2024transportbench | Syed et al., 2024 | arXiv:2408.08302 | arXiv only as of 2026-10 (also an SSRN preprint, 10.2139/ssrn.4931447) | - | fetched | 23 | 89,511 |
| verga2024poll | Verga et al., 2024 | arXiv:2404.18796 | arXiv only as of 2026-10 | - | fetched | 17 | 78,180 |
| wang2019superglue | Wang et al., 2019 | arXiv:1905.00537 | NeurIPS 2019 | - | fetched | 29 | 92,841 |
| wang2018glue | (named as Glue in the intro, not in the reference list) | arXiv:1804.07461 | ICLR 2019 (also BlackboxNLP 2018 workshop at EMNLP) | - | fetched | 20 | 72,809 |
| wang2023newton | Wang et al., 2023 | arXiv:2310.07018 | Findings of EMNLP 2023, pp. 9743-9758 | - | fetched | 18 | 70,403 |
| xie2023storytelling | Xie et al., 2023a | https://aclanthology.org/2023.inlg-main.23.pdf | INLG 2023, pp. 323-351 | resolved (RESOLVE): ACL Anthology 2023.inlg-main.23 | fetched | 29 | 94,878 |
| xie2023deltascore | Xie et al., 2023b | https://aclanthology.org/2023.findings-emnlp.353.pdf | Findings of EMNLP 2023, pp. 5317-5331 | resolved (RESOLVE): ACL Anthology 2023.findings-emnlp.353 | fetched | 15 | 57,555 |
| xie2025finchain | Xie et al., 2025 | arXiv:2506.02515 | ACL 2026 (Long Papers), pp. 14529-14553 | resolved (RESOLVE): arXiv 2506.02515 | fetched | 25 | 101,022 |
| xu2025ugphysics | Xu et al., 2025 | arXiv:2502.00334 | ICML 2025 (poster) | - | fetched | 29 | 98,079 |
| zhang2020bertscore | Zhang et al., 2020 | arXiv:1904.09675 | ICLR 2020 | - | fetched | 43 | 130,564 |
| zhang2025physreason | Zhang et al., 2025a | arXiv:2502.12054 | ACL 2025 (Long Papers), pp. 16593-16615 | - | fetched | 23 | 86,141 |
| zhang2025abenchphysics | Zhang et al., 2025b | arXiv:2507.04766 | arXiv only as of 2026-10 | - | fetched | 7 | 28,272 |
| zhao2023survey | Zhao et al., 2023 | arXiv:2303.18223 | Frontiers of Computer Science 20(12), article 2012627 (2026); the arXiv version is still updated | - | fetched | 144 | 863,614 |
| zhou2025engibench | Zhou et al., 2025 | arXiv:2509.17677 | Findings of ACL 2026, pp. 36308-36334 | - | fetched | 28 | 116,792 |
| chen2023frugalgpt | Chen et al., 2023 | arXiv:2305.05176 | TMLR (published 2024-12) | - | fetched | 13 | 44,863 |
| suzgun2022bbh | (named as BBH in the intro; attributed to Kazemi et al. 2025) | arXiv:2210.09261 | Findings of ACL 2023, pp. 13003-13051 | - | fetched | 49 | 123,180 |
| srivastava2022bigbench | (named as BIG-Bench in Related Work; attributed to Luo et al. 2024) | arXiv:2206.04615 | TMLR (published 2023-05) | - | fetched | 95 | 400,029 |
| felten2025engibench | (named in a reference-list note: distinct from the design-focused EngiBench by Felten et al.) | arXiv:2508.00831 | NeurIPS 2025 Datasets and Benchmarks Track (poster) | resolved (RESOLVE): arXiv 2508.00831 | fetched | 43 | 117,474 |
| wang2023scibench | (named in the July rebuttal as SciBench, Wang et al. 2023; not in the paper) | arXiv:2307.10635 | ICML 2024 (poster) | - | fetched | 28 | 112,219 |
| skelic2025circuit | (the arXiv id abs/2502.07980 and the note 'Introduces the CIRCUIT benchmark' printed in the Liu et al., 2025 reference) | arXiv:2502.07980 | arXiv only as of 2026-10 (an ICLR 2025 submission that was not accepted, per OpenReview) | **added**: 2502.07980, the paper the printed id and note point to | fetched | 27 | 68,708 |

46 entries; 46 fetched; 0 failed; 0 without a text layer

## Spot-checks of the text files

The opening characters of 13 text files were read (`head -c 800`, or a little more). Each held
the paper named, with its title and authors, and none was an HTML or error page:

- `chen2025apbench`: the APBench title; "Di Wu, Raymond Zhang, Enrico M. Zucchelli, Yongchao
  Chen & Richard Linares".
- `gulati2024putnamaxiom`: "Putnam-AXIOM: A Functional & Static Benchmark for Measuring Higher
  Level Mathematical Reasoning in LLMs"; Gulati, Miranda, Chen, Xia, Fronsdal, Dumont, Obbad,
  Koyejo.
- `deng2024contamination`: the title; Deng, Zhao, Tang, Gerstein, Cohan.
- `liu2025circuit`: "Enhancing Large Language Models for Automated Homework Assessment in
  Undergraduate Circuit Analysis"; Liangliang Chen ... Ying Zhang (Georgia Tech).
- `skelic2025circuit`: "CIRCUIT: A Benchmark for Circuit Interpretation and Reasoning
  Capabilities of LLMs"; Lejla Skelic ... Ruonan Han (MIT, Analog Devices).
- `ott2022saturation`: the Nature Communications header with DOI 10.1038/s41467-022-34591-0; Ott,
  Barbosa-Silva, Blagec, Brauner, Samwald.
- `plank2022hlv`: the EMNLP 2022 proceedings header (pages 10671-10682); Barbara Plank.
- `xie2023storytelling`: the INLG 2023 proceedings header (pages 323-351); Xie, Cohn, Lau.
- `xie2023deltascore`: the Findings of EMNLP 2023 header (pages 5317-5331); Xie, Li, Cohn, Lau.
- `xie2025finchain`: the FinChain title; the 25 authors from Zhuohan Xie to Preslav Nakov.
- `felten2025engibench`: the EngiBench framework title; Felten ... Fuge (ETH Zurich, Maryland).
- `du2025supergpqa`: the SuperGPQA title; "M-A-P", "ByteDance Seed, 2077.AI".
- `zhang2025abenchphysics`: the arXiv stamp "2507.04766v1"; the title; Yiming Zhang ... Junbo
  Zhao. The file has 7 pages, which is the length of that version, not a truncation.

## Problems and judgement calls

1. **Entry count.** The brief and the README say 46 entries; `papers.json` held 45. It now holds
   46 because `skelic2025circuit` was added (problem 2).
2. **`liu2025circuit` is a conflated reference.** The reference prints "Jiaheng Liu et al. 2025.
   Enhancing large language models for automated homework assessment in undergraduate circuit
   analysis. abs/2502.07980. Introduces the CIRCUIT benchmark." arXiv 2502.07980 is "CIRCUIT: A
   Benchmark for Circuit Interpretation and Reasoning Capabilities of LLMs" (Skelic, Xu, Cox, Lu,
   Yu, Han; arXiv only, an ICLR 2025 submission that was not accepted). The printed title is a
   different paper, arXiv 2511.18221 (Chen, Xie, Qin, Guo, Rohde, Zhang; IEEE FIE 2025): an
   enhancement pipeline for GPT-4o's assessment of circuit-analysis homework, which introduces no
   benchmark (its data come from the authors' earlier benchmarking study) and cites CIRCUIT as its
   reference [16] with the id 2502.07980. Neither paper has an author "Jiaheng Liu". Per
   the brief (the id must match the title), the entry's id was corrected to 2511.18221, and
   CIRCUIT was added as `skelic2025circuit`, following the catalogue's convention for works the
   text names under another reference (`wang2018glue`, `suzgun2022bbh`, `srivastava2022bigbench`).
   The citing sentence lists specialized benchmarks siloed in narrow sub-disciplines, which fits
   CIRCUIT; if the authors meant CIRCUIT, the in-text key becomes "Skelic et al., 2025". Remove
   the added entry if only the printed title is wanted.
3. **Putnam-AXIOM: the ICML 2025 version, not the cited workshop version.** The reference cites
   the NeurIPS 2024 MATH-AI workshop paper. The catalogue now points to arXiv 2508.08292 (v2,
   2025-08-27), the ICML 2025 version (PMLR 267:20723-20747), because it is the complete and
   archival one. Its title adds "in LLMs" and it has an eighth author (Elyas Obbad); `title`
   keeps the cited wording, so the abs-page check reports a variant rather than a match. The
   workshop PDF is `https://openreview.net/pdf?id=YXnwlZe0yf` and the PMLR PDF is
   `https://raw.githubusercontent.com/mlresearch/v267/main/assets/gulati25a/gulati25a.pdf`, if
   the panel wants the cited text.
4. **The version of record over the arXiv preprint for the six URL-sourced papers.** The two Xie
   et al. 2023 papers point to their ACL Anthology PDFs, the versions the reference cites with
   page numbers (their arXiv preprints, 2301.09790 and 2303.08991, are in the abs-page record).
   `ott2022saturation` keeps the Nature Communications PDF, which is open access (CC BY 4.0), over
   its preprint 2203.04592. APBench (Scientific Reports, CC BY-NC-ND 4.0) has no arXiv version:
   Crossref links only a Research Square preprint (10.21203/rs.3.rs-5619028/v1) and Semantic
   Scholar lists no arXiv id. No journal paper was paywalled, so no preprint substitution was
   needed. For these six `arxiv_latest` is null.
5. **APBench's authors.** The reference prints Chen, Zucchelli, Jang, Wu, Zhang, Lavezzi,
   Linares. The PDF and Crossref give Di Wu, Raymond Zhang, Enrico M. Zucchelli, Yongchao Chen,
   Richard Linares: Daniel Jang and Giovanni Lavezzi are not authors, and the in-text "Chen et
   al., 2025" should read "Wu et al., 2025".
6. **Author lists in the reference list that do not match the papers** (each checked against the
   arXiv page record or the PDF):
   - Lightman et al. 2023: 30 names printed; the paper has 10 authors. Seven of the 30 are
     authors; the other 23 are not, and Bowen Baker, Teddy Lee and Karl Cobbe are missing.
   - Kazemi et al. 2025 (BBEH): six first names differ: Mostafa/Mehran Kazemi, Behnam/Bahare
     Fatemi, Christos/Chrysovalantis Anastasiou, Shweta V./Sanket Vaibhav Mehta, Lovish K./Lalit
     K. Jain, Valentina/Virginia Aglietti.
   - Zhang et al. 2025b (ABench-Physics): all 12 first names differ (surnames match), e.g.
     Yifan/Yiming Zhang, Jun/Junbo Zhao.
   - Luo et al. 2024 (BIGbench): all 9 first names differ (surnames match), e.g. Haoyang/Hanjun
     Luo, Zenghao/Zuozhu Liu.
   - Xu et al. 2025 (UGPhysics): 8 of 9 first names differ (only Tong Xiao matches), e.g.
     Xingyu/Xin Xu, Yizhong Wang/Yang Wang.
   - Zhang et al. 2025a (PhysReason): 6 of 9 differ, e.g. Xinyi/Xinyu Zhang, Jielin/Jun Liu.
   - Shojaee et al. 2025 (LLM-SRBench): 4 of 6 differ: Payam/Parshin Shojaee, Nam-Huan/Ngoc-Hieu
     Nguyen, Kamyar/Kazem Meidani, Khanh Duy/Khoa D. Doan.
   - Qiu et al. 2025 (PHYBench): Siyu/Shi Qiu, Shuo/Shaoyang Guo, Ze-Yu/Zhuo-Yang Song.
   - Mirzadeh et al. 2024 (GSM-Symbolic): Kiarash/Keivan Alizadeh, Hidetoshi/Hooman Shahrokhi.
   - Cheng et al. 2024 (ElecBench): the first author is Xiyuan Zhou in both arXiv versions (so
     "Zhou et al., 2024"); "Chao Yang" and "Xinlei Cai" are not authors of either version.
   - APBench (problem 5) and the CIRCUIT reference (problem 2).
   - Minor: "Jia Duan" for Jiafei Duan (NEWTON); "Shun Yao" for Shunyu Yao (SWE-bench).
7. **Venues the reference list misses.** FrugalGPT is printed as ICML 2023; it is TMLR (2024).
   Sixteen papers printed as arXiv preprints now have a venue: SuperGPQA, EngDesign, PHYBench
   (NeurIPS 2025 Datasets and Benchmarks), HARDMath and GSM-Symbolic (ICLR 2025), Let's Verify
   (ICLR 2024), LLM-SRBench and UGPhysics (ICML 2025), PhysReason (ACL 2025), FinChain (ACL 2026),
   EngiBench (Findings of ACL 2026), NEWTON (Findings of EMNLP 2023), MATH (NeurIPS 2021 Datasets
   and Benchmarks), SuperGLUE (NeurIPS 2019), Heesch et al. (AI 2025, Springer LNCS, whose
   version has five of the arXiv version's six authors) and the LLM survey (Frontiers of Computer
   Science, 2026). EEE-Bench is printed as arXiv with "Accepted to CVPR 2025"; Putnam-AXIOM is
   printed as the workshop paper and is now ICML 2025. One source record is itself wrong:
   BIG-bench's arXiv journal-ref reads "Transactions on Machine Learning Research, May/2022", but
   OpenReview dates its TMLR publication 2023-05-11 and DBLP lists it under TMLR 2023.
8. **The fetched arXiv file is the latest version.** Three are 2026 versions: `zhao2023survey`
   v19 (2026-03-18), `xie2025finchain` v4 (2026-04-30) and `zhou2025engibench` v2 (2026-05-02);
   every other latest version is from 2025 or earlier. Whether a 2026 version postdates what the
   May authors read is not checked here.
9. **The fetch shared its folder with another session.** The first run, as briefed, without
   `--file`, also picked up the candidate catalogues another session had begun writing (it
   printed "DUPLICATE key wang2026prime: kept engineering_web.json, dropped
   process_supervision.json"). It was stopped after four cited papers; those four files, not yet
   recorded in `MANIFEST.json`, were deleted and fetched again by a run limited to
   `--file papers.json`, so every cited record is a first-hand download. `fetch_papers.py`
   rewrites `MANIFEST.json` whole at the end of each run, so runs from different sessions that
   overlap can drop each other's records. At 17:24Z it held all 46 cited records (plus 115
   candidate records), and 30 candidate PDFs on disk had no record; `fetch_papers.py --file
   <catalogue>` records files that are already present. The script also reports candidate keys
   that appear in two candidate files (`wang2026prime`, `azmi2026scirho`,
   `ansari2026physicsregrading` at 17:25Z); the second copy is dropped.
10. **Lookup limits.** This session's web-search budget ran out after the first searches, the
    DBLP API sits behind a bot check, and the OpenReview API began requiring a challenge after
    about 50 queries; the venues therefore come from the sources listed above, all queried
    directly.
