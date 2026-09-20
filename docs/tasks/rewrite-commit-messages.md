# Task: rewrite commit messages to one line, with no Co-Authored-By trailer

**Written:** 2026-09-20. **For:** a separate session. The evaluator-pilot session is
still committing to `master`, so read "Before you start" first.

## Goal

Every commit message on `master` becomes a **single short line**: the existing
subject line, with no body and no `Co-Authored-By: Claude …` trailer. Nothing else
changes: no file contents, authors, author dates or tree. The one exception is the
hash remapping below, which is needed so the repo does not point at commits that no
longer exist.

This is the owner's standing preference for new commits too. See the memory
`commit-message-style`: `git commit -m "<short subject>"`, no body, no trailer.

## Facts measured on 2026-09-20 (HEAD `bdee611`)

- `master` has **287 commits**. **240** carry the trailer, from 2026-09-05 onward.
  Almost all of those also have a multi-paragraph body.
- Authors: 280 by `ayeshag7`; 7 by Usman (`Usman`, `usmansafdarktk`,
  `Usmansafdarktk`). Rewrite messages only: authors, committers' identities and
  author dates stay as they are.
- **Remote:** `origin = https://github.com/usmansafdarktk/EngTrace.git`, a repo shared
  with Usman. The rewrite needs a force-push to `origin/master` (and to
  `origin/pilot/phase1-template-authoring`). Every existing clone will diverge.
  Usman must be told beforehand and must re-clone or hard-reset afterwards.
  `origin/gh-pages` is separate history; check whether it needs anything.
- **Branches based on `master`:** 16 local branches (`redesign/*` ×14, `pilot/*` ×2).
  14 have no commits of their own; `pilot/llm-annotation-pilot` has 1 and
  `pilot/phase1-template-authoring` has 5. After the rewrite they must point at the
  rewritten commits: rewrite all refs in one pass (`--refs --all`, or equivalent).
- **75 distinct commit hashes are cited inside 49 tracked files.** Once hashes change,
  each citation must be replaced by the new hash of the same commit, or it dangles.
  The files that are more than documentation:
  - `tests/template_integrity/regen_inventory.py`: **code** uses
    `AUDIT_REV = '67d41f4'` (lines 9, 174, 521). If this is not remapped, the
    script breaks.
  - `evaluator_pilot_17092026/slice/FREEZE.json`: `"commit":
    "d43024656e37c08520f7d86e9b5e5c8fca35e806"` and a `renamed_from` note citing
    `d43024656e37`. This is the provenance of the frozen evaluation slice. Remap
    both. The per-item SHA-256 values in that file are **content hashes, not
    commits**: do not touch them.
  - The rest are docs: `docs/re-implementation-sep/DECISIONS.md` (10),
    `reviews/*`, `track-a/*`, `track-b/*`, `docs/prompts/*`,
    `evaluator_pilot_17092026/FINDINGS.md` and `README.md`.

  The scan that produced this list matched every 7–40-hex-character token in tracked
  files against `git log` hashes. Re-run it: the list will have grown.

## Before you start

1. **Coordinate with the evaluator-pilot session.** It commits and pushes to `master`
   after every step. Do the rewrite only while nothing else is committing, and make
   sure that session re-fetches (`git fetch && git reset --hard origin/master` on a
   clean tree) before it commits again.
2. **Tell Usman** (repo owner) before force-pushing.
3. **Back up:** `git bundle create ../EngTrace-pre-rewrite.bundle --all`. Also note
   the old `master` SHA, so the old history can be restored exactly.
4. Work in a **fresh clone**, not in the pilot's working copy. It holds untracked and
   ignored pilot data (`evaluator_pilot_17092026/.venv`, scores, traces) that must not
   be disturbed.

## Suggested procedure

`git filter-repo` (pip install git-filter-repo) does the rewrite in one pass over all
refs:

1. **Messages.** A `--message-callback` keeps only the first line:
   `return message.split(b'\n', 1)[0].rstrip() + b'\n'`.
   filter-repo already rewrites abbreviated hashes that appear **inside commit
   messages**. It does not touch file contents.
2. **Hashes inside files.** File contents are not remapped automatically. Take the
   old→new commit map filter-repo writes (`.git/filter-repo/commit-map`). Then either:
   - (a) **Recommended.** Make one follow-up commit on the rewritten `master` that
     replaces every cited old hash with its new one, at the same abbreviation length.
     History stays message-only, and the fix is a single reviewable diff.
   - (b) Rewrite the blobs in history too. This is not recommended: changing a blob
     changes its commit's hash, so the map must be built incrementally in
     topological order.
3. **Verify** before pushing:
   - `git log --format=%B | grep -c Co-Authored-By` → 0, and every commit's `%b` is
     empty.
   - `git diff <old-master> <new-master> --stat` shows only the hash-remap commit's
     files. Per commit, `git diff-tree` trees are identical to the old ones (compare
     old vs new tree hashes through the commit map).
   - `python tests/template_integrity/regen_inventory.py` (or its test) still runs
     against the remapped `AUDIT_REV`.
   - Every hash cited in a tracked file resolves: `git cat-file -e <hash>^{commit}`.
   - Branch tips: each `pilot/*` branch still has its own 1 or 5 commits on top of
     the rewritten base.
4. **Push:** `git push --force-with-lease origin master` and the pilot branch. Then
   check GitHub shows the new history.

## Things to leave alone

- The evaluator pilot's data and caches, and anything under `.venv`.
- `Co-Authored-By` or trailer-like text **inside files** (docs quoting commit
  conventions): only commit messages change.
- Content hashes (`sha256` fields, `config_sha256`, `trace_sha256`, `item_sha256`):
  none of these are commit hashes.
