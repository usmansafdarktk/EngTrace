# Related work, October 2026 revision

The Related Work section of the May 2026 submission (`docs/_ARR_May__EngTrace.pdf`, section 2) was
last revised for that cycle. This folder is the record of its revision for the 12 October 2026
submission: an audit of every paper it cites, a search for what has appeared since, a panel review
of both, and the rewritten section with its bibliography.

Status: in progress (started 2026-10-01). This README is completed when the work is.

## Layout

| Path | What |
|---|---|
| `papers.json` | The cited set: every benchmark, evaluation-method and positioning paper the May text cites, plus the four works it names but attributes to a different reference, and SciBench from the July rebuttal. 46 entries. |
| `candidates/<topic>.json` | New papers found in the October 2026 search, one file per search topic. |
| `fetch_papers.py` | Acquires the PDFs and extracts their text; `MANIFEST.json` records source, SHA-256, size, pages and retrieval time for every file. PDFs and text are gitignored; the script re-acquires them. |
| `make_bib.py` | Builds BibTeX records from arXiv, the ACL Anthology and Crossref for the keys the revised text cites. |
| `notes/` | The working record: the state-of-EngTrace brief the writers used, the citation-context map of the May paper, the search notes per topic, the resolution of the cited set, and the review protocol (written before the reviews). |
| `reviews/cited/`, `reviews/candidate/` | The three independent reviews of every paper, one JSON file each. |
| `CITED_SET_AUDIT.md`, `NEW_PAPERS_AUDIT.md` | The consolidated audits. |
| `RELATED_WORK_v3.md`, `related_work_v3.tex`, `related_work_sources.bib` | The rewritten section, its LaTeX, and the bibliography entries it needs. |
| `CHANGES.md` | Every change from the May text with the audit entry that justifies it. |

## How to re-acquire the papers

    python docs/related_work_oct2026/fetch_papers.py            # fetch what is missing
    python docs/related_work_oct2026/fetch_papers.py --verify   # re-hash what is on disk
