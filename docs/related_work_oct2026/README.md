# Related work, October 2026 revision

The Related Work section of the May 2026 submission (`docs/_ARR_May__EngTrace.pdf`, section 2) was last
revised for that cycle. This folder is the record of its revision for the 12 October 2026 submission. It
holds an audit of every paper the May text cites, a search for what has appeared since, panel reviews,
and the rewritten section with its bibliography.

Status: complete, 2026-10-02, including the follow-up review of the 18 main-text papers the panel had
not reached and a final pass over the section (`notes/followup_review.md`, D-167). Start with
`RELATED_WORK_v3.md`, then `CHANGES.md`.

## The deliverables

| File | What |
|---|---|
| `RELATED_WORK_v3.md` | The revised section 2 (four paragraphs and a positioning paragraph) and an appendix with a comparison table (Table A) and further related work; notes for the authors at the end. |
| `related_work_v3.tex` | The same in LaTeX (`\citep`/`\citet`, a `table*` for Table A). Generated: edit `related_work_v3.src.md` and run `render_related_work.py`. |
| `related_work_sources.bib` | The 100 entries the section cites, built from arXiv, Crossref and ACL Anthology records by `make_bib.py`, with published venues applied and title capitals protected. |
| `intro_methods_sources.bib` | Nine more entries used by the replacement sentences for the Introduction and sections 3.3, 4 and 4.3. |
| `CHANGES.md` | Every change from the May text with its reason; replacement sentences for the Introduction and other sections; anonymity decisions for the authors; what could not be confirmed. |
| `CITED_SET_AUDIT.md` | The audit of the May citations: three independent verdicts per paper, the consolidated verdict, what v3 does, and the reference-list corrections. |
| `NEW_PAPERS_AUDIT.md` | Every candidate from the October search, its panel reviews if any, and whether v3 cites it, with the reason if not. |

## How it was done, and the scripts that reproduce it

| Step | Record | Script |
|---|---|---|
| Catalogue the May citations (46 works, including those the text names but mis-attributes) and acquire their PDFs | `papers.json`, `MANIFEST.json`, `notes/cited_set_resolution.md` | `fetch_papers.py`, `check_arxiv_abs.py`, `check_text.py` |
| Map every citation in sections 1 and 2 to the claim it supports | `notes/citation_context_map.md` | — |
| Summarise EngTrace's current state for the writers | `notes/engtrace_state_brief.md` | — |
| Search for new work in six topics; download what was found | `candidates/*.json`, `notes/search_*.md` | `fetch_papers.py` |
| Panel review, three independent reviewers per paper (protocol written before the reviews) | `notes/review_protocol.md`, `notes/review_prompt_template.md`, `notes/review_batches_*.md`, `reviews/` | `make_candidate_batches.py`, `tally_reviews.py` → `notes/review_tally.md` |
| Check every claim the new text makes against the paper's own text | `facts.py` → `notes/fact_check.md` (200 of 200 found) | `verify_facts.py` |
| Build the bibliography; render the section | `related_work_sources.bib`, `related_work_v3.*` | `make_bib.py`, `render_related_work.py` |
| Generate the new-papers audit | `NEW_PAPERS_AUDIT.md` | `make_new_papers_audit.py` |
| Follow-up review: the 18 main-text papers without a panel review, one reviewer each; a final pass over the section; uncertain venues checked in OpenReview, Semantic Scholar, Crossref and the ACL Anthology | `reviews/candidate/*.r1.json`, `notes/followup_review.md`, `notes/venue_check.json` | `check_arxiv_abs.py`, `check_venues.py` |

The panel covered all 46 cited works (three reviews each) and the highest-priority candidates. It was
stopped on 2026-10-02 at the owner's request, to save cost. From then on, the fact check stood in for
reviewers on the remaining candidates: no statement in the new text rests on an unchecked description.
The 18 candidates the main text cited without a panel review were then each read by one reviewer, for
anything that would change the positioning as well as for the sentence that cites them
(`notes/followup_review.md`).

## How to re-acquire the papers

PDFs (`papers/`), extracted text (`text/`) and raw BibTeX (`bib_raw/`) are gitignored for size and
copyright. To rebuild them:

    python docs/related_work_oct2026/fetch_papers.py            # fetch what is missing
    python docs/related_work_oct2026/fetch_papers.py --verify   # re-hash what is on disk
    python docs/related_work_oct2026/make_bib.py --keys docs/related_work_oct2026/related_work_keys.txt

Two candidates (PE Civil Bench, PSE-Bench) are on ScienceDirect, which blocks scripts. Their `text/` files
hold the OpenAlex abstract only.
