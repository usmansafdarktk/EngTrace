#!/usr/bin/env python
"""Look up the venue of record for papers whose venue rests on one source (an arXiv comment or the
search catalogue), in four free indexes. DBLP's search API was tried first; it answers scripts with a
bot check.

- Semantic Scholar: all records in one batch request (the unauthenticated pool answers 429 to a
  request per paper), by the arXiv ids that notes/arxiv_abs_check.json holds: venue, publication
  venue, journal, year, external ids.
- OpenReview: a title search per paper; the record keeps the venue and venue id of every note whose
  title matches the catalogue's after normalisation (an accepted paper's venue id is the conference's,
  a rejected one's ends in Rejected_Submission).
- Crossref: for a key with a DOI, the journal, volume, issue date and authors.
- ACL Anthology: for a key with an anthology id, the paper page's title, pages and first authors.

Nothing is changed: the record goes to notes/venue_check.json, and the follow-up review
(notes/followup_review.md) reads it.

Usage: python check_venues.py
"""
import glob
import json
import os
import re
import sys
import time
from datetime import datetime, timezone

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "notes", "venue_check.json")
UA = {"User-Agent": "EngTrace-related-work-venue-check/1.0 (research use)"}
S2_FIELDS = "title,venue,publicationVenue,journal,year,externalIds,publicationTypes"

KEYS = {
    # the three single-source venues of CHANGES.md section 6
    "huang2026verifierrobustness": None,
    "dlugosz2026gsmsymbolicreeval": None,
    "naser2026eri": "10.1016/j.cie.2026.112333",   # the related DOI on its arXiv page (notes/arxiv_abs_check.json)
    # venues that rest on the search catalogue alone (followup_review.md)
    "wang2026prime": None,
    "zhao2026prismphysics": None,
    "yu2026hipho": None,
    "imani2025sympybench": None,
    "huang2025mathperturb": None,
    "srivastava2024functionalbench": None,
    "kevian2024controlbench": None,
    "mostajabdaveh2025orqa": None,
    "li2025atmosscibench": None,
}
# ACL Anthology ids found on the event pages (aclanthology.org/events/acl-2026/, /events/eacl-2026/)
ANTHOLOGY = {"wang2026prime": "2026.acl-long.683", "imani2025sympybench": "2026.eacl-industry.8"}


def norm(t):
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def openreview(title):
    r = requests.get("https://api2.openreview.net/notes/search",
                     params={"term": title, "source": "forum", "limit": 10}, headers=UA, timeout=(20, 60))
    r.raise_for_status()
    out = []
    for n in r.json().get("notes", []):
        c = n.get("content", {})
        val = lambda f: (c.get(f) or {}).get("value") if isinstance(c.get(f), dict) else c.get(f)
        if norm(val("title"))[:60] == norm(title)[:60]:
            out.append({"venue": val("venue"), "venueid": val("venueid"), "forum": n.get("forum")})
    return out


def anthology(pid):
    r = requests.get(f"https://aclanthology.org/{pid}/", headers=UA, timeout=(20, 60))
    r.raise_for_status()
    title = re.search(r"<title>(.*?)</title>", r.text, re.S)
    pages = re.search(r"Pages:?</dt>\s*<dd[^>]*>(.*?)</dd>", r.text, re.S)
    return {"id": pid, "title": re.sub(r"\s+", " ", title.group(1)).replace(" - ACL Anthology", "") if title else None,
            "pages": re.sub(r"<[^>]+>|\s+", " ", pages.group(1)).strip() if pages else None,
            "first_authors": re.findall(r'<a href="?/people/[^>]*>([^<]+)</a>', r.text)[:4]}


def arxiv_ids():
    with open(os.path.join(HERE, "notes", "arxiv_abs_check.json"), encoding="utf-8") as fh:
        return {k: v.get("arxiv") for k, v in json.load(fh)["entries"].items()}


def catalogue_titles():
    titles = {}
    for path in glob.glob(os.path.join(HERE, "candidates", "*.json")) + [os.path.join(HERE, "papers.json")]:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        items = data if isinstance(data, list) else data.get("candidates", data.get("papers", []))
        for e in items:
            if isinstance(e, dict) and e.get("key") and e.get("title"):
                titles.setdefault(e["key"], e["title"])
    return titles


def s2_batch(arxiv_ids_in_order):
    """One POST for all papers; a busy pool (429) is retried with a growing wait."""
    for attempt in range(8):
        r = requests.post("https://api.semanticscholar.org/graph/v1/paper/batch",
                          params={"fields": S2_FIELDS},
                          json={"ids": [f"arXiv:{a}" for a in arxiv_ids_in_order]},
                          headers=UA, timeout=(20, 120))
        if r.status_code != 429:
            break
        time.sleep(10 * (attempt + 1))
    r.raise_for_status()
    out = []
    for m in r.json():
        if m is None:
            out.append(None)
            continue
        pv = m.get("publicationVenue") or {}
        out.append({"title": m.get("title"), "venue": m.get("venue"), "publication_venue": pv.get("name"),
                    "journal": m.get("journal"), "year": m.get("year"), "types": m.get("publicationTypes"),
                    "external_ids": m.get("externalIds")})
    return out


def crossref(doi):
    r = requests.get(f"https://api.crossref.org/works/{doi}", headers=UA, timeout=(20, 60))
    r.raise_for_status()
    m = r.json()["message"]
    return {"container": m.get("container-title"), "volume": m.get("volume"), "article": m.get("article-number"),
            "issued": m.get("issued", {}).get("date-parts"), "type": m.get("type"), "title": m.get("title"),
            "authors": [" ".join(x for x in (a.get("given"), a.get("family")) if x) for a in m.get("author", [])]}


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    titles = catalogue_titles()
    ids = arxiv_ids()
    keys = list(KEYS)
    record = {k: {"title": titles.get(k), "arxiv": ids.get(k)} for k in keys}
    try:
        for k, s2 in zip(keys, s2_batch([ids[k] for k in keys])):
            record[k]["semantic_scholar"] = s2 if s2 is not None else "no record"
    except Exception as ex:              # noqa: BLE001 - the error is the record
        for k in keys:
            record[k]["semantic_scholar_error"] = f"{type(ex).__name__}: {ex}"
    for k, doi in KEYS.items():
        time.sleep(1.5)
        try:
            record[k]["openreview"] = openreview(titles[k])
        except Exception as ex:          # noqa: BLE001
            record[k]["openreview_error"] = f"{type(ex).__name__}: {ex}"
        if doi:
            try:
                record[k]["crossref"] = {"doi": doi, **crossref(doi)}
            except Exception as ex:      # noqa: BLE001
                record[k]["crossref_error"] = f"{type(ex).__name__}: {ex}"
        if k in ANTHOLOGY:
            try:
                record[k]["acl_anthology"] = anthology(ANTHOLOGY[k])
            except Exception as ex:      # noqa: BLE001
                record[k]["acl_anthology_error"] = f"{type(ex).__name__}: {ex}"
    for k in keys:
        rec = record[k]
        ss = rec.get("semantic_scholar") if isinstance(rec.get("semantic_scholar"), dict) else {}
        journal = (ss.get("journal") or {}).get("name") if isinstance(ss.get("journal"), dict) else ss.get("journal")
        venue = " / ".join(str(x) for x in (ss.get("venue"), ss.get("publication_venue"), journal) if x)
        orv = "; ".join(f"{o['venue']} ({o['venueid']})" for o in rec.get("openreview", []))
        print(f"{k:32s} S2: {venue or rec.get('semantic_scholar_error') or 'no venue'}"
              + f" | OpenReview: {orv or rec.get('openreview_error') or 'no match'}"
              + (f" | Crossref: {rec['crossref']['container']} vol {rec['crossref']['volume']}" if "crossref" in rec else "")
              + (f" | ACL Anthology: {rec['acl_anthology']['id']} pp. {rec['acl_anthology']['pages']}" if "acl_anthology" in rec else ""))
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump({"generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                   "generator": "docs/related_work_oct2026/check_venues.py", "entries": record},
                  fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"record in {os.path.relpath(OUT, HERE)}")


if __name__ == "__main__":
    main()
