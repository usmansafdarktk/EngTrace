#!/usr/bin/env python
"""Check every arXiv id in a related-work catalogue against its arXiv abstract page.

For each entry with an "arxiv" id, fetch https://arxiv.org/abs/<id> and record the page's title,
authors, comments, journal-ref, related DOI and submission history (every version and its UTC
date). The page title is compared with the entry's "title" after normalisation (case, punctuation
and spacing ignored; a title that only extends the other is reported as a variant), and the latest
version's date with the entry's "arxiv_latest" field when the entry has one. Nothing in the
catalogue is changed: the record goes to notes/arxiv_abs_check.json, so every title check and
date in the notes can be re-derived.

Usage
  python check_arxiv_abs.py                     papers.json (the May 2026 cited set)
  python check_arxiv_abs.py --file candidates/x.json
  python check_arxiv_abs.py --only KEY          one entry (repeatable)
  python check_arxiv_abs.py --id 2303.08991     an id outside the catalogue (repeatable)
"""
import argparse
import difflib
import html
import json
import os
import re
import sys
import time
from datetime import datetime, timezone

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "notes", "arxiv_abs_check.json")
UA = {"User-Agent": "Mozilla/5.0 (compatible; EngTrace-related-work-fetch/1.0; research use)"}
PAUSE = 3.0  # seconds between arXiv requests


def clean(fragment):
    text = re.sub(r"<[^>]+>", " ", fragment or "")
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def norm(title):
    t = html.unescape(title or "").lower().replace("’", "'").replace("‘", "'")
    return re.sub(r"[^a-z0-9]", "", t)


def field(page, css):
    m = re.search(r'<td class="tablecell ' + css + r'[^"]*">(.*?)</td>', page, re.S)
    return clean(m.group(1)) if m else None


def parse(page):
    m = re.search(r'<h1 class="title mathjax">(.*?)</h1>', page, re.S)
    title = clean(re.sub(r'<span class="descriptor">Title:</span>', "", m.group(1))) if m else None
    m = re.search(r'<div class="authors">(.*?)</div>', page, re.S)
    authors = [clean(a) for a in re.findall(r"<a [^>]*>(.*?)</a>", m.group(1), re.S)] if m else []
    m = re.search(r'<div class="submission-history">(.*?)</div>', page, re.S)
    versions = []
    for v, stamp in re.findall(r"\[v(\d+)\](?:</a>)?</strong>\s*([A-Z][a-z]{2}, \d{1,2} [A-Z][a-z]{2} \d{4} "
                               r"\d{2}:\d{2}:\d{2}) UTC", m.group(1) if m else ""):
        day = datetime.strptime(stamp, "%a, %d %b %Y %H:%M:%S").strftime("%Y-%m-%d")
        versions.append({"version": int(v), "date": day})
    doi = re.findall(r"10\.\d{4,9}/[^\s\"<>]+", field(page, "doi") or "")
    return {"abs_title": title, "abs_authors": authors, "comments": field(page, "comments"),
            "journal_ref": field(page, "jref"), "related_doi": doi[0] if doi else None,
            "versions": versions,
            "latest_version": versions[-1]["version"] if versions else None,
            "latest_date": versions[-1]["date"] if versions else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", default="papers.json")
    ap.add_argument("--only", action="append")
    ap.add_argument("--id", action="append", default=[], help="an arXiv id outside the catalogue")
    args = ap.parse_args()

    path = args.file if os.path.isabs(args.file) else os.path.join(HERE, args.file)
    with open(path, encoding="utf-8") as fh:
        entries = json.load(fh)
    todo = [e for e in entries if e.get("arxiv") and (not args.only or e["key"] in args.only)]
    todo += [{"key": "id_" + i.replace(".", "_"), "arxiv": i, "title": None} for i in args.id]

    record = {}
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as fh:
            record = json.load(fh).get("entries", {})
    problems = 0
    for n, e in enumerate(todo):
        if n:
            time.sleep(PAUSE)
        aid = re.sub(r"^(arxiv:|https?://arxiv.org/(abs|pdf)/)", "", str(e["arxiv"]).strip(), flags=re.I)
        url = f"https://arxiv.org/abs/{aid}"
        try:
            r = requests.get(url, headers=UA, timeout=(20, 60))
            status = r.status_code
            info = parse(r.text) if status == 200 else {}
        except requests.RequestException as ex:
            status, info = f"{type(ex).__name__}: {ex}", {}
        rec = {"arxiv": aid, "url": url, "http": status, "catalogue_title": e.get("title"), **info,
               "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
        if e.get("title") is None:
            rec["title_check"] = "no catalogue title"
        elif not info.get("abs_title"):
            rec["title_check"] = "NO TITLE ON PAGE"
        elif norm(info["abs_title"]) == norm(e["title"]):
            rec["title_check"] = "match"
        elif norm(info["abs_title"]).startswith(norm(e["title"])) or norm(e["title"]).startswith(norm(info["abs_title"])):
            rec["title_check"] = "variant (one title extends the other; see the note)"
        else:
            ratio = difflib.SequenceMatcher(None, norm(info["abs_title"]), norm(e["title"])).ratio()
            rec["title_check"] = f"MISMATCH (similarity {ratio:.2f})"
        if "arxiv_latest" in e:
            rec["arxiv_latest_check"] = ("match" if e["arxiv_latest"] == info.get("latest_date")
                                         else f"MISMATCH (catalogue {e['arxiv_latest']})")
        bad = rec["title_check"].startswith(("MISMATCH", "NO TITLE")) or \
            str(rec.get("arxiv_latest_check", "")).startswith("MISMATCH")
        problems += bad
        record[e["key"]] = rec
        print(f"{'PROBLEM ' if bad else 'ok      '}{e['key']:28s} {aid:12s} v{rec.get('latest_version')} "
              f"{rec.get('latest_date')}  title {rec['title_check']}"
              f"{'  | page: ' + str(rec.get('abs_title')) if bad else ''}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                   "generator": "docs/related_work_oct2026/check_arxiv_abs.py",
                   "entries": dict(sorted(record.items()))}, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\n{len(todo)} checked, {problems} problems; record in {os.path.relpath(OUT, HERE)}")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
