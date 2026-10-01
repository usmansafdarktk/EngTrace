# The reviewer prompt, verbatim

Each reviewer agent in the panel (see `review_protocol.md`) receives this prompt with the
placeholders filled: `{SET}` is `cited` or `candidate`, `{REVIEWER}` is `r1`, `r2` or `r3`, and
`{PAPERS}` is the batch: for each paper its key, its text file, and for the cited set the May
sentence(s) that cite it (from `citation_context_map.md`) or, for a candidate, the search note's
one-paragraph description and the topic it was found under.

---

You are one reviewer on an independent panel auditing the literature behind the Related Work
section of EngTrace, an NLP benchmark paper being revised for ARR (deadline 12 October 2026).
Today is 2026-10-01. Work in C:\Users\ayesha.gull01\EngTrace; read files with Bash (cat, sed -n,
grep -n); run no git command that changes anything; write only the review files named below.

First read, in full: docs/related_work_oct2026/notes/engtrace_state_brief.md (what EngTrace now
is), docs/related_work_oct2026/notes/review_protocol.md (the schema you must follow), and the May
paper's Introduction and Related Work: sed -n 44,206p docs/_ARR_May__EngTrace.txt (the text dump
drops some spaces between words).

Then review each paper in your batch. For each one:

1. Read its text file docs/related_work_oct2026/text/<key>.txt: the abstract, introduction, the
   benchmark or method construction section and the evaluation section in full; skim the rest.
   The file has `===== PAGE n =====` markers; cite pages. If the file is not the paper named
   (wrong paper, HTML, no text layer), say so in the review and set "read" to "not the paper".
2. Establish what the paper actually is: task and domain; how the items were built (textbooks,
   exams, generated from templates, perturbed, expert-written); size; whether instances are static
   or generated; whether it provides gold intermediate steps or traces; how it evaluates (final
   answer, multiple choice, rubric, LLM judge, unit tests, step-level or process); what human
   validation it reports; any contamination measures; and the current venue if the text shows one.
3. For a CITED paper: judge the May sentence that cites it on the protocol's scale (accurate,
   imprecise, misleading, wrong) with page references, and recommend keep as is / keep, reword /
   drop / replace. For a CANDIDATE: judge whether a careful reviewer of EngTrace would expect it
   to be cited and recommend add (must) / add (should) / add (optional) / drop, saying which
   paragraph it belongs to (math and coding; physical sciences; engineering; process supervision
   and judge validity; or the introduction).
4. Write ONE sentence or clause that cites the paper accurately, in the register of the May text,
   and list the facts the revised section may state about it (each with a page reference).
5. Say what the paper shares with EngTrace and what it lacks, concretely (generated instances;
   gold traces; step-level evaluation; engineering domain; expert validation of the evaluator).

Write docs/related_work_oct2026/reviews/{SET}/<key>.{REVIEWER}.json for every paper, exactly in the
protocol's schema (valid JSON; validate each file with python -c "import json;json.load(open(...))").
Be exact and sceptical: do not credit a benchmark with step-level evaluation because it mentions
chain-of-thought, and do not call a benchmark static if it generates variants. Where the May
sentence groups several papers under one characterisation, judge the characterisation for THIS
paper only. Do not consult the web; the paper text and the files named here are the evidence.

Final message: one line per paper: key, verdict (cited) or recommendation (candidate), confidence,
and anything the consolidator must know (a wrong file, a naming clash, a surprise).

Your batch ({SET}, reviewer {REVIEWER}):
{PAPERS}
