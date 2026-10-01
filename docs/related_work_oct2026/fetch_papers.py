#!/usr/bin/env python
"""Fetch, verify and text-extract the papers catalogued for the October 2026 related-work revision.

Catalogue files: papers.json (the May 2026 submission's cited set) and candidates/*.json (new papers
found in the October 2026 search). Each is a JSON list of entries:

  {"key": "mirzadeh2024gsmsymbolic",          unique; the file stem for the PDF and the text
   "set": "cited" | "candidate",
   "cited_as": "Mirzadeh et al., 2024",       how the May paper cites it (cited set only)
   "title": "...", "authors": "...", "year": 2024, "venue": "...",
   "arxiv": "2410.05229",                     preferred source: https://arxiv.org/pdf/<id>
   "url": null,                               used when arxiv is null; must serve a PDF
   "notes": "..."}

Usage
  python fetch_papers.py              fetch every entry whose PDF is missing, then extract text
  python fetch_papers.py --only KEY   one entry (repeatable)
  python fetch_papers.py --file F     only the entries of one catalogue file (path or basename)
  python fetch_papers.py --verify     re-hash what is on disk against MANIFEST.json
  python fetch_papers.py --list       print the plan, fetch nothing
  python fetch_papers.py --text       re-extract text for every PDF on disk

A download is judged by its content, never by its HTTP status: it must start with %PDF, be larger
than 10 kB, open in PyMuPDF with at least one page, and yield at least 1,000 characters of text
(otherwise it is kept but flagged "no text layer"). PDFs (papers/) and text (text/) are gitignored;
MANIFEST.json records the source URL, SHA-256, bytes, pages, text characters and retrieval time of
every file, and the reason for every failure, so the acquisition is reproducible and checkable.
"""
import argparse
import glob
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime, timezone

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
PAPERS = os.path.join(HERE, "papers")
TEXT = os.path.join(HERE, "text")
MANIFEST = os.path.join(HERE, "MANIFEST.json")
UA = {"User-Agent": "Mozilla/5.0 (compatible; EngTrace-related-work-fetch/1.0; research use)",
      "Accept": "application/pdf,*/*;q=0.8"}
ARXIV_PAUSE = 3.0  # seconds between arXiv requests (their rate-limit guidance)
MIN_BYTES = 10_000
MIN_CHARS = 1_000


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def catalogue_files():
    files = [os.path.join(HERE, "papers.json")]
    files += sorted(glob.glob(os.path.join(HERE, "candidates", "*.json")))
    return [f for f in files if os.path.exists(f)]


def load_catalogue(only_file=None):
    entries, seen, duplicates = [], {}, []
    for f in catalogue_files():
        if only_file and os.path.basename(f) != os.path.basename(only_file):
            continue
        with open(f, encoding="utf-8") as fh:
            try:
                data = json.load(fh)
            except json.JSONDecodeError as e:
                sys.exit(f"{f}: invalid JSON: {e}")
        if not isinstance(data, list):
            sys.exit(f"{f}: expected a JSON list of entries")
        for e in data:
            key = e.get("key")
            if not key or not re.fullmatch(r"[a-z0-9_]+", key):
                sys.exit(f"{f}: bad or missing key in entry {e}")
            if key in seen:
                duplicates.append({"key": key, "kept": seen[key], "dropped": os.path.basename(f)})
                continue
            seen[key] = os.path.basename(f)
            e["_file"] = os.path.basename(f)
            entries.append(e)
    return entries, duplicates


def source_url(e):
    a = e.get("arxiv")
    if a:
        a = str(a).strip()
        a = re.sub(r"^(arxiv:|https?://arxiv.org/(abs|pdf)/)", "", a, flags=re.I)
        a = re.sub(r"(\.pdf)$", "", a)
        return f"https://arxiv.org/pdf/{a}"
    u = e.get("url")
    if u:
        return u.strip()
    return None


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url):
    """Return (bytes, final_url, error). Three attempts with backoff on network errors, 429 and 5xx."""
    delays = [5, 15, 45]
    last = None
    for attempt in range(3):
        try:
            r = requests.get(url, headers=UA, timeout=(20, 120), allow_redirects=True)
            if r.status_code == 200:
                return r.content, r.url, None
            last = f"HTTP {r.status_code}"
            if r.status_code not in (429, 500, 502, 503, 504):
                return None, r.url, last
        except requests.RequestException as ex:
            last = f"{type(ex).__name__}: {ex}"
        if attempt < 2:
            time.sleep(delays[attempt])
    return None, url, last


def validate_pdf(data):
    if not data.startswith(b"%PDF"):
        head = data[:200].decode("latin-1", "replace").replace("\n", " ")
        return f"not a PDF (starts with: {head[:80]!r})"
    if len(data) < MIN_BYTES:
        return f"too small ({len(data)} bytes)"
    return None


def extract_text(pdf_path, txt_path):
    import fitz  # PyMuPDF
    with fitz.open(pdf_path) as doc:
        pages = doc.page_count
        parts = []
        for i, page in enumerate(doc):
            parts.append(f"\n\n===== PAGE {i + 1} =====\n\n")
            parts.append(page.get_text())
    text = "".join(parts)
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return pages, len(text)


def load_manifest():
    if os.path.exists(MANIFEST):
        with open(MANIFEST, encoding="utf-8") as fh:
            return json.load(fh)
    return {"generated_utc": None, "generator": "docs/related_work_oct2026/fetch_papers.py",
            "files": {}, "failures": {}, "unresolved": [], "duplicates": []}


def save_manifest(m):
    m["generated_utc"] = utc_now()
    m["files"] = dict(sorted(m["files"].items()))
    with open(MANIFEST, "w", encoding="utf-8") as fh:
        json.dump(m, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


def record(manifest, e, url, final, data, pages, chars, retrieved):
    manifest["files"][e["key"]] = {
        "set": e.get("set"), "title": e.get("title"), "source_url": url, "final_url": final,
        "sha256": sha256(data), "bytes": len(data), "pages": pages, "text_chars": chars,
        "no_text_layer": chars < MIN_CHARS, "retrieved_utc": retrieved, "catalogue": e["_file"]}
    manifest["failures"].pop(e["key"], None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", action="append", help="fetch this key only (repeatable)")
    ap.add_argument("--file", help="only the entries of this catalogue file")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--text", action="store_true", help="re-extract text for PDFs on disk")
    ap.add_argument("--force", action="store_true", help="re-download even if the PDF exists")
    args = ap.parse_args()

    os.makedirs(PAPERS, exist_ok=True)
    os.makedirs(TEXT, exist_ok=True)
    entries, duplicates = load_catalogue(args.file)
    if args.only:
        entries = [e for e in entries if e["key"] in set(args.only)]
        missing = set(args.only) - {e["key"] for e in entries}
        if missing:
            sys.exit(f"unknown keys: {sorted(missing)}")
    manifest = load_manifest()
    manifest["duplicates"] = duplicates
    for d in duplicates:
        print(f"DUPLICATE key {d['key']}: kept {d['kept']}, dropped {d['dropped']}")

    if args.verify:
        problems = 0
        for key, rec in manifest["files"].items():
            p = os.path.join(PAPERS, key + ".pdf")
            if not os.path.exists(p):
                print(f"MISSING  {key}")
                problems += 1
                continue
            with open(p, "rb") as fh:
                h = sha256(fh.read())
            if h != rec["sha256"]:
                print(f"CHANGED  {key}: {h[:12]} != {rec['sha256'][:12]}")
                problems += 1
        print(f"{len(manifest['files'])} files, {problems} problems")
        return

    if args.list:
        for e in entries:
            have = "have" if os.path.exists(os.path.join(PAPERS, e["key"] + ".pdf")) else "----"
            print(f"{have}  {e['key']:40s} {e.get('set', '?'):9s} {source_url(e) or 'UNRESOLVED'}")
        return

    if args.text:
        for e in entries:
            p = os.path.join(PAPERS, e["key"] + ".pdf")
            if os.path.exists(p):
                pages, chars = extract_text(p, os.path.join(TEXT, e["key"] + ".txt"))
                rec = manifest["files"].setdefault(e["key"], {})
                rec.update({"pages": pages, "text_chars": chars, "no_text_layer": chars < MIN_CHARS})
                print(f"text     {e['key']}: {pages} pages, {chars} chars")
        save_manifest(manifest)
        return

    unresolved, fetched, failed = [], 0, 0
    last_arxiv = 0.0
    for e in entries:
        key = e["key"]
        pdf = os.path.join(PAPERS, key + ".pdf")
        txt = os.path.join(TEXT, key + ".txt")
        url = source_url(e)
        if not url:
            unresolved.append(key)
            print(f"UNRESOLVED {key}: no arxiv id and no url")
            continue
        if os.path.exists(pdf) and not args.force:
            if not os.path.exists(txt) or key not in manifest["files"]:
                pages, chars = extract_text(pdf, txt)
                with open(pdf, "rb") as fh:
                    data = fh.read()
                old = manifest["files"].get(key, {})
                record(manifest, e, url, old.get("final_url", url), data, pages, chars,
                       old.get("retrieved_utc", utc_now()))
                print(f"present  {key}: text extracted ({pages} pages, {chars} chars)")
            else:
                print(f"present  {key}")
            continue
        if "arxiv.org" in url:
            wait = ARXIV_PAUSE - (time.time() - last_arxiv)
            if wait > 0:
                time.sleep(wait)
        data, final, err = fetch(url)
        if "arxiv.org" in url:
            last_arxiv = time.time()
        if err is None:
            err = validate_pdf(data)
        if err:
            failed += 1
            manifest["failures"][key] = {"source_url": url, "error": err, "at_utc": utc_now()}
            print(f"FAILED   {key}: {err}  [{url}]")
            continue
        with open(pdf, "wb") as fh:
            fh.write(data)
        try:
            pages, chars = extract_text(pdf, txt)
        except Exception as ex:  # a corrupt file that still starts with %PDF
            os.remove(pdf)
            failed += 1
            manifest["failures"][key] = {"source_url": url, "error": f"PyMuPDF: {ex}", "at_utc": utc_now()}
            print(f"FAILED   {key}: PyMuPDF could not open it  [{url}]")
            continue
        record(manifest, e, url, final, data, pages, chars, utc_now())
        fetched += 1
        flag = "  (no text layer!)" if chars < MIN_CHARS else ""
        print(f"fetched  {key}: {len(data)} bytes, {pages} pages, {chars} chars{flag}")
    manifest["unresolved"] = unresolved
    if not args.only and not args.file:  # a full run: drop records of keys no catalogue lists any more
        listed = {e["key"] for e in entries}
        for key in [k for k in manifest["files"] if k not in listed]:
            print(f"pruned   {key}: no longer in any catalogue")
            manifest["files"].pop(key)
    save_manifest(manifest)
    print(f"\n{fetched} fetched, {failed} failed, {len(unresolved)} unresolved, "
          f"{len(manifest['files'])} files in MANIFEST.json")


if __name__ == "__main__":
    main()
