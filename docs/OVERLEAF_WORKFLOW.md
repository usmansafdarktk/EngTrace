# Working with the Overleaf project

Split out of `docs/PAPER_PLAN_OCT2026.md` on 2026-10-03; revised 2026-10-04.

## The rule

The paper lives on Overleaf. `current_overleaf_project/` is a local mirror of it for editing and review here, and it
is **not versioned in this repository** (gitignored on 2026-10-04). Changes made locally are pushed to Overleaf;
changes co-authors make on Overleaf are pulled into the mirror before the next local edit.

## Three ways to sync

1. **Git integration** (a premium feature; MBZUAI's institutional licence may provide it: Overleaf menu, "Sync",
   "Git"). Each project has a git URL, `https://git.overleaf.com/<project-id>`, and an authentication token (Account
   settings, "Git integration"; the username is `git`, the token is the password). Make the mirror a clone of the
   Overleaf repository:

   ```bash
   rm -rf current_overleaf_project
   git clone https://git.overleaf.com/<project-id> current_overleaf_project
   cd current_overleaf_project && git pull     # before editing
   git add -A && git commit -m "..." && git push   # after editing
   ```

   The clone is its own repository inside the ignored folder, so it does not interact with this repository's
   history. If the Git option is in the menu, send the project's git URL and the clone can be set up here.
2. **GitHub sync** (also premium) needs the paper at the root of its own GitHub repository; not worth it beside
   option 1.
3. **Upload** (free accounts): Overleaf's "Upload" button accepts several files at once and overwrites same-named
   files, so the changed `.tex` files and `figs/` can be dropped in after each local edit. Copy-paste is for a single
   small edit. A zip upload creates a new project, so it is not an update path. Co-authors' edits come back through
   "Download as zip" (replace the mirror) or the project history.

## Compiling

There is no TeX installation on this machine, so the page count and the final check are compiled on Overleaf unless
MiKTeX (pdfLaTeX, the engine the ACL template assumes; `\pdfoutput=1` is a pdfTeX primitive) is installed locally.
`docs/render_md_pdf.py` renders Markdown working documents through headless Edge and is not for the paper.

## Before submission

Switch `\usepackage[preprint]{acl}` to `review`, remove the author block and the project and code links, and point
the code link at the anonymised archive.
