# Writing rules

These apply to every file in this folder. The rules on content (plain language, no overstatement, no internal
vocabulary) are in `docs/PAPER_PLAN_OCT2026.md`, section 8; this file covers form.

## The authors' rules

1. All text is written in LaTeX (`.tex`), ready to drop into the Overleaf project.
2. Benchmark names are set with `\textsc{}`, spelled as their authors spell them: `\textsc{MMLU}`, `\textsc{MATH}`,
   `\textsc{HumanEval}`, `\textsc{GSM8K}`, `\textsc{PhysReason}`, `\textsc{FinChain}`.
3. Never type "EngTrace"; always write `\ourdataset` (defined in `main.tex` with `\xspace`, so a following space or
   punctuation mark needs nothing extra). The one exception is text quoted verbatim in a code listing, where a macro
   does not expand.
4. A citation is tied to the word before it with a tilde, never a space: `templates~\citep{mirzadeh2024gsmsymbolic}`.
   The same holds for every cross-reference: `in~\autoref{fig:example}`, never `in \autoref{...}`.
5. "LLM" and "LLMs" are used without spelling them out. Experts are "domain experts" throughout.

## Conventions of the current Overleaf project

- `\citep{}` for parenthetical citations, several keys in one command (`\citep{a, b}`); `\citet{}` when the authors
  are part of the sentence (`\citet{xie2025finchain} build ...`).
- Cross-references with `\autoref{}`, tied with a tilde: `in~\autoref{sec:evaluation}` prints "§4", "Figure 2",
  "Table 1" or "Appendix A". Labels carry a prefix: `sec:`, `subsec:`, `fig:`, `tab:`, `eq:`, `appendix:`.
- Run-in headings with `\paragraph{}`; few numbered subsections.
- Numbers: thousands with a comma (`2,250`); percent as `\%`; ranges in tables as `0.81--0.98`; em dash as `---`.
- Quotation marks as ``` ``this'' ```, never straight double quotes; `e.g.,` with its comma; F1 written as `F1`.
- An abbreviation is spelled out at first use, in the abstract and again in the body: "large language models (LLMs)".

## From the ACL instructions (`current_overleaf_project/formatting.md`)

- Abstract: at most 200 words. By convention it has no citations and no references to sections, figures or tables.
- Long paper: eight pages of content; references and appendices do not count, but the paper must stand on its own.
- Review version is anonymous: no names, affiliations, project links or acknowledgements; own prior work is cited in
  the third person (`\citet{x} show`, never "we previously showed").
- Captions go below figures and tables; figures must read in greyscale.
- Cite the refereed version of a work when one exists; give DOIs where available.

## From the ARR writing guide (`ARR_ How to Write a Research Article.docx`)

- American English throughout (modeling, favor, artifacts); never mix spellings.
- Active voice ("we draw the instances", not "the instances are drawn"); no tense switching within a paragraph.
- Define every abbreviation before use, even LLM; mark a term at its definition with `\emph{}`, not quotes.
- One term per concept; "measure", not "metric"; "significant" only with a statistical test; "open-weights" LLMs;
  "use", not "employ".
- The Introduction mentions only the closest work, to contrast it with ours, and ends with bullet-point
  contributions that say why the work matters; the abstract and conclusion repeat the same points.
- Every figure and table is discussed in the text; captions are in sentence case and end with a full stop.
- Bibliography: cite the published version, not arXiv; booktitle "Proceedings of ...", a `series` field
  (e.g. `ICLR~'25`), the conference location as `address`, no `month`; an arXiv-only paper is `@article` with
  `journal={ArXiv preprint}` and `volume={arXiv:NNNN.NNNNN}`; protect capitals only where needed (`{LLM}`);
  check every entry by hand. `custom.bib` holds only the citations the written sections use;
  `custom_old.bib` is the previous bibliography, kept for reference.

## This folder

- One file per section, named as `main.tex` inputs it (`0_abstract.tex`, `1_intro.tex`, `2_relatedwork.tex`, ...);
  the bibliography is `custom.bib`.
- In the `.tex` files every sentence starts on a new line and wraps at 100 characters, never inside a citation,
  reference or `\textsc{}` name; changes then show sentence by sentence, and LaTeX still prints one paragraph.
  Prompt boxes copied verbatim from the May appendix keep their original source lines.
- The appendix lives in `appendices/`: `7_appendix.tex` inputs one file per appendix section
  (`taxonomy_content.tex`, ...), with paths relative to `main.tex` (`sections/appendices/...`). On Overleaf the folder
  is `sections/appendices/`, and `main.tex` inputs `sections/appendices/7_appendix` in place of `sections/7_appendix`.

## Cross-reference labels

Labels used before their sections exist. Each must be defined when its section is written, or the reference prints
"??" in the PDF.

- To define: `sec:evaluation` (Section 4), `sec:experiments` (Section 5), `appendix:templates` (template
  construction and examples), `appendix:statistics` (dataset statistics and the evaluation set),
  `appendix:certification` (template certification), `appendix:scoring` (scoring details).
- Defined: `sec:intro`, `sec:related`, `sec:benchmark`, `subsec:taxonomy`, `subsec:templates`,
  `subsec:certification`, `fig:engtrace-overview` (the five-branch figure, `figs/engtrace-overview-5branch.pdf`).
- Defined in `appendices/taxonomy_content.tex`, with the May appendix's labels: `sec:appendix_domain_validation`,
  `sec:appendix_area_validation`, `sec:appendix_books` and `sec:appendix_data_sources` (one section, two tables:
  `tab:foundational_textbooks`, `tab:authoritative_sources`), `sec:appendix_significance_scoring`.
