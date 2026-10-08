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
- Headings follow the ACL instructions (`current_overleaf_project/formatting.md`), which describe two numbered
  levels: `\section` (12 pt bold, "1 Introduction") and `\subsection` (11 pt bold, "6.1 File Format"), numbered by
  `acl.sty`. No `\subsubsection`. Each heading stands on its own line, its `\label` on the next. Few subsections.
- Below them, run-in headings with `\paragraph{Title.}`: bold, in the text size, ending with a full stop, and run
  into the text, so the first sentence follows the heading on the same source line
  (`\paragraph{Milestones.} We rerun ...`).
- A new paragraph is a blank line; `acl.sty` indents it by 1em (the 0.4 cm the instructions ask for), except the first
  after a heading. No `\noindent`, `\\` or `\vspace` in running text.
- Headings in title case throughout (the ARR guide asks for consistent capitalization of headings).
- Numbers: thousands with a comma (`2,250`); percent as `\%`; ranges in tables as `0.81--0.98`; em dash as `---`.
- Quotation marks as ``` ``this'' ```, never straight double quotes; `e.g.,` with its comma; F1 written as `F1`.
- Model names in `\texttt{}`, spelled as their developers write them (`\texttt{MiMo-V2.5-Pro}`,
  `\texttt{GPT-5.4 mini}`), and never wrapped across lines.
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
  the bibliography is `custom.bib`. Section 5 is two files, `6_experiments.tex` (the section heading and the setup)
  and `6_results.tex` (results and error analysis), which `6_experiments_results.tex` inputs, so `main.tex` keeps
  `\input{sections/6_experiments_results}`; the wrapper compiles once `6_results.tex` exists.
- A `\paragraph{}` heading and the first sentence of its paragraph share a source line (the authors, 2026-10-05;
  Section 5's appendix so far, the earlier sections still start the text on the next line).
- The final state only (the authors, 2026-10-05): no history of the work in the paper, such as rounds of fixes,
  which analyses were planned or added, when a margin or threshold was fixed, or which test replaced which; and no
  cost figures. Appendices give implementation details, in lists where they read better than prose.
- In the `.tex` files every sentence starts on a new line and wraps at 100 characters, never inside a citation,
  reference or `\textsc{}` name; changes then show sentence by sentence, and LaTeX still prints one paragraph.
  Prompt boxes keep the prompt's own lines, which may exceed 100 characters.
- The appendix lives in `appendices/`: `7_appendix.tex` inputs one file per appendix section
  (`taxonomy_content.tex`, ...), with paths relative to `main.tex` (`sections/appendices/...`). On Overleaf the folder
  is `sections/appendices/`, and `main.tex` inputs `sections/appendices/7_appendix` in place of `sections/7_appendix`.

## Cross-reference labels

Labels used before their sections exist. Each must be defined when its section is written, or the reference prints
"??" in the PDF.

- Defined in `6_experiments.tex`: `sec:experiments`, `subsec:models` (the evaluated model suite) and `subsec:setup`
  (the experimental setup), one paragraph each with no paragraph headings (the authors, 2026-10-05); every number in
  it comes from `full_run_28092026/paper_setup.py`, whose `--check` must pass.
- Defined in `appendices/models.tex`: `appendix:models` (the selection rule; table `tab:models`, in the May layout)
  and `appendix:setup`, whose subsections are `appendix:decoding` (table `tab:decoding`), `appendix:prompt` (the
  inference prompt exactly as the run sent it, and output processing) and `appendix:protocol` (the statistical
  protocol in full). `full_run_28092026/paper_setup.py --check` checks every number, table row and identifier, and
  that the prompt box is the template every response used.
- Defined: `sec:intro`, `sec:related`, `sec:benchmark`, `subsec:taxonomy`, `subsec:templates`,
  `subsec:certification`, `fig:engtrace-overview` (the five-branch figure, `figs/engtrace-overview-5branch.pdf`),
  `fig:template-generation` (May's pipeline figure, `figs/template-generation.pdf`).
- Defined in `appendices/taxonomy_content.tex`, with the May appendix's labels: `sec:appendix_domain_validation`,
  `sec:appendix_area_validation`, `sec:appendix_books` and `sec:appendix_data_sources` (one section, two tables:
  `tab:foundational_textbooks`, `tab:authoritative_sources`), `sec:appendix_significance_scoring`, and
  `appendix:statistics` (dataset statistics; tables `tab:difficulty_stats`, `tab:template_variation`; every number
  from `docs/appendix_statistics.py`, whose `--check` must pass).
- Defined in `appendices/template_examples.tex`: `sec:appendix_param_details` (four run-in paragraphs by principle,
  with examples from all five branches; D-196) and `sec:appendix_template_examples` (three listings generated from the current template
  code by `docs/appendix_listings.py`, whose `--check` must pass after any template change).
- Defined in `appendices/certification.tex`: `appendix:certification` (integrity checks, LLM screen with May's
  prompt box, expert protocol with `tab:planted_defects`, results with `tab:expert_agreement` and the rounds
  in the prose (D-196); every number from `docs/appendix_certification.py`, whose `--check` must pass). Table rows may exceed 100 characters; prose may not.
- Defined in `5_evaluation.tex`: `sec:evaluation` (no subsections: the measures and the scoring process, with the
  expert study behind a pointer to `appendix:validation`), `eq:match` (the final-answer match rule) and
  `eq:coverage` (Milestone Coverage).
- Defined in `appendices/scoring.tex`: `appendix:scoring` (one table, `tab:scoring_settings`; the two judge prompts).
  Defined in `appendices/validation.tex`: `appendix:validation` (one table, `tab:validation_results`). Every number
  in these two files and in `5_evaluation.tex` comes from `docs/appendix_evaluation.py`, whose `--check` must pass.
