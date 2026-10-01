# Search notes: engineering_arxiv (incomplete record)

The agent that ran this search (engineering benchmarks found through the arXiv API, citation
chasing and the Hugging Face papers index, as a second route beside the web search in
`search_engineering_web.md`) was terminated by the account's session limit on 2026-10-01 while
appending its last entries, before it wrote this note. Its catalogue, `candidates/engineering_arxiv.json`
(30 entries), is valid and complete as far as it got, and every entry's PDF and text were fetched.
The query log and the per-candidate descriptions it would have written here are lost; the
one-sentence rationale per candidate survives in each entry's `notes` field, and the panel reviews
in `reviews/candidate/` describe each paper from its text.

Entries this file shares with `engineering_web.json` (same key, same paper) are kept once by
`fetch_papers.py`, which takes the first catalogue file in alphabetical order, so the
`engineering_arxiv` entry is the one listed for those keys.
