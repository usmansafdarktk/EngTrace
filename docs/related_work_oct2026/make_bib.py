#!/usr/bin/env python
"""Build BibTeX entries for the catalogued papers from authoritative sources.

For every entry in papers.json and candidates/*.json (or only the keys given with --only), fetch the
BibTeX record from the source that holds it:
  arXiv            https://arxiv.org/bibtex/<id>                (the entry's "arxiv" field)
  ACL Anthology    https://aclanthology.org/<anthology id>.bib  (an aclanthology.org url)
  a DOI            https://doi.org/<doi> with Accept: application/x-bibtex (an "doi" field)
Each record is rewritten to carry the catalogue key as its BibTeX key and, when the catalogue entry
has a "venue_now" field that names a conference or journal, a "note" field with it, so the
bibliography can say where a preprint was published without inventing a booktitle. The raw records
are kept untouched under raw/ for checking.

Usage
  python make_bib.py                 all entries -> related_work_sources.bib
  python make_bib.py --only KEY ...  some entries
  python make_bib.py --keys FILE     the keys listed one per line in FILE (the ones the text cites)
  python make_bib.py --out FILE      write elsewhere
"""
import argparse
import glob
import json
import os
import re
import sys
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "bib_raw")
UA = {"User-Agent": "Mozilla/5.0 (compatible; EngTrace-related-work-bib/1.0; research use)"}


def load_entries():
    files = [os.path.join(HERE, "papers.json")] + sorted(glob.glob(os.path.join(HERE, "candidates", "*.json")))
    seen, out = set(), []
    for f in files:
        if not os.path.exists(f):
            continue
        with open(f, encoding="utf-8") as fh:
            for e in json.load(fh):
                if e["key"] in seen:
                    continue
                seen.add(e["key"])
                out.append(e)
    return out


def get(url, headers=None):
    h = dict(UA)
    if headers:
        h.update(headers)
    for attempt, delay in enumerate((0, 5, 15)):
        if delay:
            time.sleep(delay)
        try:
            r = requests.get(url, headers=h, timeout=(20, 60), allow_redirects=True)
        except requests.RequestException as ex:
            err = f"{type(ex).__name__}"
            continue
        if r.status_code == 200 and "@" in r.text:
            return r.text, None
        err = f"HTTP {r.status_code}"
        if r.status_code not in (429, 500, 502, 503, 504):
            break
    return None, err


def source_for(e):
    a = e.get("arxiv")
    if a:
        a = re.sub(r"^(arxiv:|https?://arxiv.org/(abs|pdf)/)", "", str(a).strip(), flags=re.I)
        return "arxiv", f"https://arxiv.org/bibtex/{a}", None
    if e.get("doi"):
        return "doi", f"https://doi.org/{e['doi']}", {"Accept": "application/x-bibtex"}
    u = e.get("url") or ""
    m = re.search(r"aclanthology\.org/([A-Za-z0-9.\-]+?)(?:\.pdf)?/?$", u)
    if m:
        return "acl", f"https://aclanthology.org/{m.group(1)}.bib", None
    return None, None, None


def rekey(bib, key, venue_now):
    bib = bib.strip()
    bib = re.sub(r"^@(\w+)\{[^,]*,", lambda m: f"@{m.group(1)}{{{key},", bib, count=1)
    if venue_now and not re.search(r"arxiv only", venue_now, re.I):
        bib = bib.rstrip()
        if bib.endswith("}"):
            bib = bib[:-1].rstrip().rstrip(",") + f",\n  note = {{{venue_now}}}\n}}"
    return bib + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--keys", help="file with one key per line")
    ap.add_argument("--out", default=os.path.join(HERE, "related_work_sources.bib"))
    args = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    entries = load_entries()
    wanted = None
    if args.only:
        wanted = set(args.only)
    if args.keys:
        with open(args.keys, encoding="utf-8") as fh:
            wanted = {line.strip() for line in fh if line.strip() and not line.startswith("#")}
    if wanted is not None:
        missing = wanted - {e["key"] for e in entries}
        if missing:
            sys.exit(f"keys not in any catalogue: {sorted(missing)}")
        entries = [e for e in entries if e["key"] in wanted]
    records, problems = [], []
    last_arxiv = 0.0
    for e in entries:
        kind, url, headers = source_for(e)
        if not url:
            problems.append((e["key"], "no arXiv id, DOI or ACL Anthology url; write the entry by hand"))
            continue
        raw_path = os.path.join(RAW, e["key"] + ".bib")
        if os.path.exists(raw_path):
            with open(raw_path, encoding="utf-8") as fh:
                text = fh.read()
        else:
            if kind == "arxiv":
                wait = 3.0 - (time.time() - last_arxiv)
                if wait > 0:
                    time.sleep(wait)
            text, err = get(url, headers)
            if kind == "arxiv":
                last_arxiv = time.time()
            if text is None:
                problems.append((e["key"], f"{err} from {url}"))
                continue
            with open(raw_path, "w", encoding="utf-8") as fh:
                fh.write(text)
        records.append(rekey(text, e["key"], e.get("venue_now")))
        print(f"ok       {e['key']}  ({kind})")
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("% Built by docs/related_work_oct2026/make_bib.py from arXiv, the ACL Anthology and Crossref records.\n")
        fh.write("% Keys are the catalogue keys in papers.json and candidates/*.json. A note field carries the\n")
        fh.write("% publication venue recorded in the catalogue (venue_now) when a preprint has since been published.\n\n")
        fh.write("\n".join(records))
    print(f"\n{len(records)} entries written to {args.out}")
    for k, p in problems:
        print(f"PROBLEM  {k}: {p}")


if __name__ == "__main__":
    main()
