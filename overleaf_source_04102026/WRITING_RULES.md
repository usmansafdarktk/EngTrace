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

## This folder

- One file per section, named as `main.tex` inputs it (`0_abstract.tex`, `1_intro.tex`, `2_relatedwork.tex`, ...);
  the bibliography is `custom.bib`.
- One sentence per line in the `.tex` files, so changes show sentence by sentence; LaTeX prints them as one paragraph.

## Cross-reference labels

Labels used before their sections exist. Each must be defined when its section is written, or the reference prints
"??" in the PDF.

- To define: `sec:benchmark` (Section 3), `sec:evaluation` (Section 4), `sec:experiments` (Section 5),
  `appendix:templates` (template examples), `appendix:scoring` (scoring details).
- Defined: `sec:intro`, `fig:example`.
