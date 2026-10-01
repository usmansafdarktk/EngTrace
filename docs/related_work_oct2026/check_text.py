#!/usr/bin/env python
"""Check that each extracted text file holds the paper its catalogue entry names.

For every entry of a catalogue file, read text/<key>.txt (written by fetch_papers.py) and look in
its first two pages for (a) the entry's title and (b) the surname of the entry's first author (for
a team author such as "M-A-P Team", the team's name). Text is compared after Unicode NFKC folding,
lower-casing, "&" read as "and", and removal of everything but letters and digits, so line breaks,
hyphenation, ligatures and punctuation do not matter. The title result is "full" when the whole
title is found, otherwise the longest leading part of the title that is found (in characters of
the folded title); a short leading part, or a missing first author, means the file may be a
different paper, an error page or a cover sheet, and is flagged. Nothing is written.

Usage
  python check_text.py                    papers.json (the May 2026 cited set)
  python check_text.py --file candidates/x.json
"""
import argparse
import json
import os
import re
import sys
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_PREFIX = 30  # folded title characters that must be found when the full title is not


def fold(s):
    s = unicodedata.normalize("NFKC", s or "").lower().replace("&", "and")
    return re.sub(r"[^a-z0-9]", "", s)


def first_pages(path, n=2):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    parts = re.split(r"\n\n===== PAGE \d+ =====\n\n", text)
    return " ".join(parts[1:n + 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", default="papers.json")
    args = ap.parse_args()
    path = args.file if os.path.isabs(args.file) else os.path.join(HERE, args.file)
    with open(path, encoding="utf-8") as fh:
        entries = json.load(fh)
    flagged = 0
    for e in entries:
        txt = os.path.join(HERE, "text", e["key"] + ".txt")
        if not os.path.exists(txt):
            print(f"NO TEXT  {e['key']}")
            flagged += 1
            continue
        page = fold(first_pages(txt))
        title = fold(e.get("title"))
        n = len(title)
        while n and title[:n] not in page:
            n -= 1
        found = "full" if n == len(title) else f"{n}/{len(title)}"
        first = (e.get("authors") or "").split(",")[0].strip()
        words = first.split()
        if len(words) > 1 and words[-1].lower() == "team":  # a team author, e.g. "M-A-P Team"
            words = words[:-1]
        surname = fold(words[-1]) if words else ""
        has_author = bool(surname) and surname in page
        bad = (n < len(title) and n < MIN_PREFIX) or not has_author
        flagged += bad
        print(f"{'FLAG    ' if bad else 'ok      '}{e['key']:28s} title {found:9s} "
              f"first author {first or '?'}: {'found' if has_author else 'NOT FOUND'}")
    print(f"\n{len(entries)} entries, {flagged} flagged")
    sys.exit(1 if flagged else 0)


if __name__ == "__main__":
    main()
