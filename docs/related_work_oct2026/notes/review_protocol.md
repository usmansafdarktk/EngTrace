# Review protocol for the October 2026 related-work revision

Written 2026-10-01, before any review was run, so the procedure is on record ahead of its results.

## What is reviewed

Two sets of papers, catalogued in this folder and acquired by `fetch_papers.py` (PDF in `papers/`,
extracted text in `text/`, provenance in `MANIFEST.json`):

- **The cited set** (`papers.json`): every benchmark, evaluation-method and positioning paper the
  May 2026 submission cites in its Introduction and Related Work, plus the four works the text names
  but attributes to a different reference (GLUE, BBH, BIG-Bench, the design-focused EngiBench) and
  SciBench, which the July rebuttal cites. Model cards, textbooks, standards and the inter-rater
  statistics papers are not reviewed: the related work does not turn on them.
- **The candidate set** (`candidates/*.json`): papers found in the October 2026 search, by topic
  (`engineering_web`, `engineering_arxiv`, `physics_science`, `symbolic_contamination`,
  `process_supervision`, `llm_judge`). The search notes are `notes/search_<topic>.md`.

The question for a cited paper is whether the May text uses it accurately and whether it still earns
its place. The question for a candidate is whether a careful reviewer of EngTrace would expect it to
be cited, and what exactly it should be cited for.

## The panel

Each paper is read independently by **three reviewer agents**, each given the paper's text file,
the sentence the May paper uses it in (from `notes/citation_context_map.md`; null for candidates),
EngTrace's current state (`notes/engtrace_state_brief.md`), and this schema. Reviewers do not see
each other's output. A reviewer must read the paper's abstract, introduction, benchmark or method
construction and evaluation sections in full and skim the rest; every factual claim about the paper
must carry a page reference (the `===== PAGE n =====` markers in the text file).

One **consolidator agent per set** reads all three reviews of every paper, records where they agree
and disagree, re-reads the paper text wherever the verdicts differ, and writes the audit
(`CITED_SET_AUDIT.md`, `NEW_PAPERS_AUDIT.md`). The audits, not the individual reviews, are what the
rewrite uses; the individual reviews stay in `reviews/` as the record.

**Amendment, 2026-10-02.** The first launch of the cited-set panel (18 reviewers at once, on
2026-10-01) was cut off by the account's session limit before any batch finished, and the
candidate set turned out to hold 159 papers. The cited set keeps three reviewers per paper, run
in two waves. The candidate set is reviewed with a panel sized by priority: three reviewers for
every must-cite candidate, two for should-cite candidates in the engineering group (the paragraph
the paper's positioning turns on) and one elsewhere, and one for optional candidates. The
consolidator still re-reads any paper whose single review it doubts. `make_candidate_batches.py`
records the batches and counts in `review_batches_candidate.md`.

## Per-paper review file

`reviews/<set>/<key>.<reviewer>.json` (set is `cited` or `candidate`; reviewer is `r1`, `r2`, `r3`):

```json
{
 "key": "mirzadeh2024gsmsymbolic",
 "reviewer": "r1",
 "read": "full | partial: <sections read>",
 "title_confirmed": true,
 "venue_now": "ICLR 2025 (page 1 header)",
 "what_it_is": "Three to five sentences: task, how the items were built, size, what is measured.",
 "domain": "grade-school math | undergraduate physics | engineering: <branches> | ...",
 "size": "items, templates, models evaluated",
 "source_of_items": "textbooks | exams | competitions | generated from templates | perturbed from an existing set | expert-written",
 "instances": "static | generated | perturbed",
 "has_gold_traces": "yes | no | partial: <what intermediate supervision exists>",
 "evaluation": "final answer | multiple choice | rubric | LLM judge | unit tests | step-level or process | expression-level metric | ...",
 "human_validation": "what human checking of items or of the evaluator the paper reports, with numbers",
 "contamination_measures": "none | generated instances | perturbations | private test set | recency | ...",
 "engtrace_use": {
   "sentence": "the May sentence that cites it, or null for a candidate",
   "verdict": "accurate | imprecise | misleading | wrong | n/a",
   "explanation": "why, with page references"
 },
 "recommendation": "keep as is | keep, reword | drop | add (must) | add (should) | add (optional) | replace with <key>",
 "suggested_wording": "one clause or sentence that cites this paper accurately in the revised related work",
 "facts_to_cite": ["the specific facts, with page refs, that the revised text may state about this paper"],
 "overlap_with_engtrace": "what it shares (generated instances, gold traces, step-level evaluation, engineering domain, expert validation) and what it lacks",
 "confidence": "high | medium | low",
 "notes": "anything else: an error in the paper, a newer version, a naming clash"
}
```

Verdict scale for a cited paper's use in the May text:

- **accurate**: the sentence states something the paper supports.
- **imprecise**: the gist is right but the characterisation is loose (e.g. a benchmark called "multiple-choice" that is mostly but not only multiple-choice).
- **misleading**: the sentence suggests something the paper does not support, or omits what matters most for the comparison (e.g. calling a benchmark "outcome-only" when it also scores steps).
- **wrong**: the citation points at the wrong paper, or the claim contradicts the paper.

## Consolidation

For every paper, the audit records: the three verdicts and recommendations; whether they agree;
for a disagreement, what the consolidator found on re-reading and which view it adopts, with the
page reference; the consolidated description (what the paper is, in two or three sentences, for the
writers); the final recommendation; and the wording the revised text may use. A summary table at
the top lists every paper with its final verdict and recommendation. For candidates the table also
carries the paragraph the paper belongs to and its priority (must, should, optional).

## What the rewrite then does

The revised Related Work cites only papers that appear in an audit with a recommendation of keep or
add, and states about each only facts the audit lists. Every change from the May text is justified
by an audit entry; `CHANGES.md` lists them.

**Amendment, 2026-10-02 (second).** After the second wave was also cut off by the session limit, the owner
asked to stop spending on reviewer agents and move to synthesis. The panel had by then reviewed all 46
cited works (138 reviews, three per paper) and part of the candidate set: the eight engineering
must-cites (three reviews each), the process-supervision must-cites (one to three each) and part of the
LLM-judge must-cites (one each). No further reviewers were launched. For every candidate the revised text
cites, each statement it makes about the paper is listed in `facts.py` with the passage that supports it,
and `verify_facts.py` checks it against the paper's text (`fact_check.md`). That check replaces the
consolidator for the candidate set, and `NEW_PAPERS_AUDIT.md` records which candidates had panel reviews.
The cited-set consolidation was done by the orchestrator from the 138 reviews (`CITED_SET_AUDIT.md`).
