# Working with the Overleaf project

Split out of `docs/PAPER_PLAN_OCT2026.md` on 2026-10-03.

## The rule

`current_overleaf_project/` in this repository is the source of truth for the paper. It is edited locally, every
change is a commit, and it is pushed to Overleaf rather than typed there. Co-authors who edit on Overleaf pull back
before the next push (below). The folder is untracked as this is written; commit it before the first edit.

## Three ways into Overleaf

1. **Git integration** (a premium feature; MBZUAI's institutional licence may provide it: Overleaf menu, "Sync",
   "Git"). Each project has a git URL, `https://git.overleaf.com/<project-id>`, and an authentication token (Account
   settings, "Git integration"; the username is `git`, the token is the password). From this repository:

   ```bash
   git remote add overleaf https://git.overleaf.com/<project-id>
   git subtree push --prefix=current_overleaf_project overleaf master      # local folder -> project root
   git subtree pull --prefix=current_overleaf_project overleaf master --squash   # Overleaf edits -> local folder
   ```

   One command each way, with history on both sides. If the Git option is in the menu, send the project's git URL
   and the remote can be set up here.
2. **GitHub sync** (also premium) needs the paper at the root of its own GitHub repository, so it would mean a
   second repository for the paper; not worth it beside option 1.
3. **Upload** (free accounts): Overleaf's "Upload" button accepts several files at once and overwrites same-named
   files, so the changed `.tex` files and `figs/` can be dropped in after each local commit. Copy-paste is for a
   single small edit. A zip upload creates a new project, so it is not an update path. Co-authors' edits come back
   through "Download as zip" (replace the local folder, review the diff, commit) or the project history.

## Compiling

There is no TeX installation on this machine, so the page count and the final check are compiled on Overleaf
unless MiKTeX (pdfLaTeX, the engine the ACL template assumes; `\pdfoutput=1` is a pdfTeX primitive) is installed
locally. `docs/render_md_pdf.py` renders Markdown working documents through headless Edge and is not for the paper.

## Before submission

Switch `\usepackage[preprint]{acl}` to `review`, remove the author block and the project and code links, and
point the code link at the anonymised archive.
