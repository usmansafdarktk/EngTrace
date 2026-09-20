#!/usr/bin/env python3
"""Remap commit hashes cited inside tracked files after a history rewrite.

The 2026-09-20 rewrite (docs/tasks/rewrite-commit-messages.md) reduced every
commit message on ``master`` to its subject line, which changed every commit
hash from the first rewritten commit onward.  Tracked files cite commit hashes
in prose, in JSON provenance and in one module constant; each citation has to
be moved to the new hash of the same commit or it dangles.

Input is the ``commit-map`` that ``git filter-repo`` writes
(``<repo>/filter-repo/commit-map``): two whitespace-separated 40-hex columns,
old hash then new hash, with a header line.

A token is rewritten only when it is a hex run of 7..40 characters that is an
unambiguous prefix of exactly one old commit hash.  The replacement keeps the
original abbreviation length, lengthened only if that prefix would be
ambiguous among the new hashes.  SHA-256 content digests (64 hex characters)
never match, because the hex run must be 40 characters or shorter and is
matched with hex-aware boundaries.

Usage::

    python tools/remap_commit_hashes.py --commit-map PATH [--apply]
    python tools/remap_commit_hashes.py --check [--commit-map PATH]

Without ``--apply`` it reports what it would change and writes nothing.
``--check`` instead reports the hash-shaped runs in tracked files that do not
resolve to a commit here; given the map as well, it fails only on the ones
that are pre-rewrite hashes, which are citations the remap missed.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# A hex run not adjacent to another hex character, so a 64-character SHA-256
# digest is never partially matched.
HEX_RUN = re.compile(rb"(?<![0-9a-fA-F])[0-9a-fA-F]{7,40}(?![0-9a-fA-F])")

MIN_LEN = 7
MAX_LEN = 40

#: Files that record the rewrite itself.  Their whole point is to hold the old
#: hashes next to the new ones, so remapping them would destroy the record.
SKIP = {"docs/tasks/commit-map-20260920.txt"}


def git(repo: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True, capture_output=True, text=True,
    )
    return out.stdout


def load_commit_map(path: Path) -> dict[str, str]:
    mapping = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 2 or len(parts[0]) != 40:
            continue  # header or blank
        old, new = parts
        if set(old) <= set("0123456789abcdef"):
            mapping[old] = new
    return mapping


def build_prefix_index(mapping: dict[str, str]) -> dict[str, str | None]:
    """prefix -> new hash, or None when the prefix is ambiguous."""
    index: dict[str, str | None] = {}
    for old, new in mapping.items():
        for length in range(MIN_LEN, MAX_LEN + 1):
            prefix = old[:length]
            if prefix in index and index[prefix] != new:
                index[prefix] = None
            else:
                index.setdefault(prefix, new)
    return index


def unambiguous_length(repo: Path, new_hash: str, wanted: int) -> int:
    """Shortest length >= wanted at which new_hash is unambiguous in repo."""
    for length in range(wanted, MAX_LEN + 1):
        probe = new_hash[:length]
        try:
            resolved = git(repo, "rev-parse", "--verify", f"{probe}^{{commit}}").strip()
        except subprocess.CalledProcessError:
            continue  # ambiguous or unknown at this length
        if resolved == new_hash:
            return length
    return MAX_LEN


def tracked_text_files(repo: Path) -> list[str]:
    names = git(repo, "ls-files", "-z").split("\0")
    return [n for n in names if n and n not in SKIP]


def commit_prefixes(repo: Path) -> set[str]:
    """Every prefix of length 7..40 of every commit reachable in repo.

    Built once from ``git rev-list --all``.  Asking git to resolve each token
    instead (``cat-file``, ``rev-parse``) is far too slow: disambiguating a
    7-character prefix searches the object database, and the corpus holds tens
    of thousands of hex runs that are not commits at all.
    """
    prefixes = set()
    for sha in git(repo, "rev-list", "--all").split():
        for length in range(MIN_LEN, MAX_LEN + 1):
            prefixes.add(sha[:length])
    return prefixes


def looks_like_a_hash(token: str) -> bool:
    """Whether a hex run is plausibly a commit id rather than a number.

    The corpus is full of long digit strings (JSONL sample values, numbers
    lifted out of PDFs) that are hex runs only by accident.  A real
    abbreviated commit mixes digits and letters; requiring both drops tens of
    thousands of those without hiding a genuine citation.
    """
    return (any(c.isdigit() for c in token)
            and any(c in "abcdef" for c in token))


def check_citations(repo: Path, mapping: dict[str, str] | None = None) -> int:
    """Report hex runs in tracked files that look like commits but do not resolve.

    Every hex run is a candidate; the ones that resolve are live citations and
    the rest are filtered by :func:`looks_like_a_hash` and listed so a human
    can judge them.  Plenty of those are legitimately unresolvable here --
    commits pinned in *other* repositories, and the ``e`` exponent of a float
    such as ``1.0033961e-05`` -- so the list on its own is advisory.

    Given a commit map, the check becomes decisive: any dangling run that is a
    prefix of a *pre-rewrite* hash is a citation the remap missed, and only
    those make this exit non-zero.
    """
    where: dict[str, list[str]] = defaultdict(list)
    for name in tracked_text_files(repo):
        path = repo / name
        try:
            data = path.read_bytes()
        except OSError:
            continue
        if b"\x00" in data:
            continue
        for match in HEX_RUN.finditer(data):
            where[match.group(0).decode("ascii").lower()].append(name)

    tokens = sorted(where)
    known = commit_prefixes(repo)
    live = {t for t in tokens if t in known}
    dangling = [t for t in tokens
                if t not in known and looks_like_a_hash(t)]

    print(f"distinct hex runs scanned:                {len(tokens)}")
    print(f"live commit citations (distinct hashes):  {len(live)}")
    print(f"hash-shaped runs that do not resolve:     {len(dangling)}")
    for token in dangling[:20]:
        files = sorted(set(where[token]))
        print(f"  {token}  in {len(files)} file(s): {files[0]}")
    if len(dangling) > 20:
        print(f"  ... and {len(dangling) - 20} more")

    if mapping is None:
        print("\nadvisory only: pass --commit-map to test these against the "
              "pre-rewrite hashes")
        return 0

    old_prefixes = {old[:n]: old for old in mapping for n in range(MIN_LEN, MAX_LEN + 1)}
    missed = [t for t in dangling if t in old_prefixes]
    print(f"\ndangling runs that are pre-rewrite commit hashes: {len(missed)}")
    for token in missed:
        print(f"  MISSED {token} -> should be "
              f"{mapping[old_prefixes[token]][:len(token)]}  in {sorted(set(where[token]))}")
    return 1 if missed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--commit-map", type=Path,
                    help="filter-repo commit-map; required unless --check")
    ap.add_argument("--repo", default=Path("."), type=Path)
    ap.add_argument("--apply", action="store_true",
                    help="write the files; otherwise report only")
    ap.add_argument("--check", action="store_true",
                    help="report dangling commit citations instead of remapping")
    args = ap.parse_args()

    repo = args.repo.resolve()
    if args.check:
        return check_citations(
            repo, load_commit_map(args.commit_map) if args.commit_map else None)
    if args.commit_map is None:
        ap.error("--commit-map is required unless --check is given")
    mapping = load_commit_map(args.commit_map)
    index = build_prefix_index(mapping)
    print(f"commit-map entries: {len(mapping)}")

    length_cache: dict[tuple[str, int], int] = {}
    per_file: dict[str, int] = {}
    cited: set[str] = set()
    ambiguous: dict[str, int] = defaultdict(int)
    unchanged_hits = 0
    changed_files = 0

    for name in tracked_text_files(repo):
        path = repo / name
        try:
            data = path.read_bytes()
        except OSError:
            continue
        if b"\0" in data:
            continue  # binary

        hits = 0

        def substitute(match: re.Match) -> bytes:
            nonlocal hits, unchanged_hits
            token = match.group(0)
            key = token.decode("ascii").lower()
            new = index.get(key)
            if new is None:
                if key in index:
                    ambiguous[key] += 1
                return token
            cited.add(key)
            if new.startswith(key):
                unchanged_hits += 1
                return token  # commit kept its hash
            cache_key = (new, len(key))
            if cache_key not in length_cache:
                length_cache[cache_key] = unambiguous_length(repo, new, len(key))
            hits += 1
            return new[: length_cache[cache_key]].encode("ascii")

        updated = HEX_RUN.sub(substitute, data)
        if hits:
            per_file[name] = hits
            changed_files += 1
            if args.apply:
                path.write_bytes(updated)

    for name in sorted(per_file):
        print(f"  {per_file[name]:4d}  {name}")
    print(f"distinct commits cited: {len(cited)}")
    print(f"citations rewritten:    {sum(per_file.values())} in {changed_files} files")
    print(f"citations already correct (hash unchanged by the rewrite): {unchanged_hits}")
    if ambiguous:
        print(f"ambiguous prefixes skipped: {dict(ambiguous)}")
    if not args.apply:
        print("\ndry run: no files written (pass --apply)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
