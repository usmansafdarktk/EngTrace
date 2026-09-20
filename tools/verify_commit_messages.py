#!/usr/bin/env python3
"""Check that every commit message is a single line with no trailer.

The owner's convention (memory ``commit-message-style``, and the rewrite in
docs/tasks/rewrite-commit-messages.md) is ``git commit -m "<short subject>"``:
one line, no body, no ``Co-Authored-By`` trailer.  This script walks the
commits reachable from the given revisions and reports any that break it.

Usage::

    python tools/verify_commit_messages.py [revisions...]   # default: --all

Exits non-zero when a commit has a body, a trailer or an empty message.

Name the branches explicitly rather than relying on the default when the
pre-rewrite history is still around: ``--all`` reaches
``backup/pre-commit-message-rewrite-20260920`` and the remote-tracking refs,
whose commits are *meant* to keep their bodies.  To check the live history::

    python tools/verify_commit_messages.py $(git for-each-ref refs/heads \\
        --format='%(refname:short)' | grep -v '^backup/')
"""

from __future__ import annotations

import subprocess
import sys

SEP = "\x01"
REC = "\x02"


def commit_messages(revs: list[str]) -> list[tuple[str, str]]:
    out = subprocess.run(
        ["git", "log", *revs, f"--format={REC}%H{SEP}%B"],
        check=True, capture_output=True, text=True, encoding="utf-8",
    ).stdout
    records = []
    for chunk in out.split(REC):
        if not chunk.strip():
            continue
        sha, _, message = chunk.partition(SEP)
        records.append((sha.strip(), message))
    return records


def main(argv: list[str]) -> int:
    revs = argv[1:] or ["--all"]
    records = commit_messages(revs)
    with_body, with_trailer, empty = [], [], []

    for sha, message in records:
        lines = [l for l in message.split("\n") if l.strip()]
        if not lines:
            empty.append(sha)
            continue
        if len(lines) > 1:
            with_body.append(sha)
        if "Co-Authored-By" in message or "Co-authored-by" in message:
            with_trailer.append(sha)

    print(f"commits checked:            {len(records)}")
    print(f"with a body (>1 line):      {len(with_body)}")
    print(f"with a Co-Authored-By line: {len(with_trailer)}")
    print(f"with an empty message:      {len(empty)}")

    for label, shas in (("body", with_body), ("trailer", with_trailer), ("empty", empty)):
        for sha in shas[:10]:
            print(f"  {label}: {sha}")

    longest = max(records, key=lambda r: len(r[1].strip()))
    print(f"longest subject: {len(longest[1].strip())} chars ({longest[0][:12]})")

    return 1 if (with_body or with_trailer or empty) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
