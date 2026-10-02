#!/usr/bin/env python
"""Build BibTeX entries for the catalogued papers from authoritative sources.

For every entry in papers.json and candidates/*.json (or only the keys given), fetch the record from
the source that holds it:
  arXiv            https://arxiv.org/bibtex/<id>                (the entry's "arxiv" field)
  a DOI            https://doi.org/<doi>, Accept: application/x-bibtex (the entry's "doi" field)
  ACL Anthology    https://aclanthology.org/<anthology id>.bib  (an aclanthology.org "url")
Authors and titles therefore come from the source, never from a hand-typed list. Then:
  - the BibTeX key is the catalogue key, or the entry's "bibkey" when the catalogue key misnames
    the first author (e.g. cheng2024elecbench -> zhou2024elecbench);
  - when the catalogue's venue ("venue_now" for the cited set, "venue" for candidates) names a
    conference or journal, an arXiv @misc record becomes @inproceedings/@article with that
    booktitle/journal and the venue's year; workshop-only venues stay @misc with a note;
  - titles from arXiv and Crossref are wrapped in double braces so styles keep their capitals.
Raw records are cached under bib_raw/ (gitignored) for checking.

Usage
  python make_bib.py --keys FILE [--out FILE]   the keys listed one per line in FILE
  python make_bib.py --only KEY ...             some entries
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

# (pattern, entry type, booktitle or journal template). Order matters: Findings before the main venue.
VENUES = [
    (r"Findings of (the )?ACL", "inproceedings", "Findings of the Association for Computational Linguistics: ACL {y}"),
    (r"Findings of (the )?EMNLP", "inproceedings", "Findings of the Association for Computational Linguistics: EMNLP {y}"),
    (r"Findings of IJCNLP-AACL", "inproceedings", "Findings of the Association for Computational Linguistics: IJCNLP-AACL {y}"),
    (r"\bNAACL\b", "inproceedings", "Proceedings of the {y} Conference of the North American Chapter of the Association for Computational Linguistics (NAACL)"),
    (r"\bEACL\b.{0,30}Industry", "inproceedings", "Proceedings of the {y} Conference of the European Chapter of the Association for Computational Linguistics (EACL): Industry Track"),
    (r"\bEACL\b", "inproceedings", "Proceedings of the {y} Conference of the European Chapter of the Association for Computational Linguistics (EACL)"),
    (r"\bEMNLP\b", "inproceedings", "Proceedings of the {y} Conference on Empirical Methods in Natural Language Processing (EMNLP)"),
    (r"\bACL\b", "inproceedings", "Proceedings of the {y} Annual Meeting of the Association for Computational Linguistics (ACL)"),
    (r"NeurIPS.{0,30}Datasets and Benchmarks", "inproceedings", "Advances in Neural Information Processing Systems (NeurIPS {y}), Datasets and Benchmarks Track"),
    (r"NeurIPS.{0,30}Evaluations", "inproceedings", "Advances in Neural Information Processing Systems (NeurIPS {y}), Evaluations and Datasets Track"),
    (r"\bNeurIPS\b", "inproceedings", "Advances in Neural Information Processing Systems (NeurIPS {y})"),
    (r"\bICLR\b", "inproceedings", "International Conference on Learning Representations (ICLR {y})"),
    (r"\bICML\b", "inproceedings", "Proceedings of the International Conference on Machine Learning (ICML {y})"),
    (r"\bCOLM\b", "inproceedings", "Conference on Language Modeling (COLM {y})"),
    (r"\bAAAI\b", "inproceedings", "Proceedings of the AAAI Conference on Artificial Intelligence ({y})"),
    (r"\bCVPR\b", "inproceedings", "Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR {y})"),
    (r"\bECAI\b", "inproceedings", "Proceedings of the European Conference on Artificial Intelligence (ECAI {y})"),
    (r"\bKDD\b", "inproceedings", "Proceedings of the ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD {y})"),
    (r"\bFIE\b|Frontiers in Education", "inproceedings", "IEEE Frontiers in Education Conference (FIE {y})"),
    (r"AI 20\d\d: Advances in Artificial Intelligence", "inproceedings", "AI {y}: Advances in Artificial Intelligence"),
    (r"\bTMLR\b", "article", "Transactions on Machine Learning Research"),
    (r"Scientific Reports", "article", "Scientific Reports"),
    (r"Nature Communications", "article", "Nature Communications"),
    (r"Frontiers of Computer Science", "article", "Frontiers of Computer Science"),
    (r"Computers & Industrial Engineering", "article", "Computers \\& Industrial Engineering"),
]
# venues that rest on one source (an arXiv comment or an OpenReview page): flagged in the output
# author strings the source garbles (arXiv prints SuperGPQA's "M-A-P Team" as "P Team")
AUTHOR_FIXES = {"du2025supergpqa": ("author={P Team and", "author={{M-A-P Team} and")}
# (naser2026eri left the list on 2026-10-02: Crossref gives the DOI on its arXiv page as Computers &
# Industrial Engineering 221, notes/venue_check.json)
SINGLE_SOURCE = {"huang2026verifierrobustness", "dlugosz2026gsmsymbolicreeval"}
# venues confirmed after the search, which the catalogue records as preprints (check_venues.py ->
# notes/venue_check.json; notes/followup_review.md). They take precedence over the catalogue's venue.
VENUE_CONFIRMED = {
    "imani2025sympybench": "EACL 2026 Industry Track (ACL Anthology 2026.eacl-industry.8)",
    "li2025atmosscibench": "NeurIPS 2025 Datasets and Benchmarks Track (OpenReview)",
}


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
    err = None
    for delay in (0, 5, 15):
        if delay:
            time.sleep(delay)
        try:
            r = requests.get(url, headers=h, timeout=(20, 60), allow_redirects=True)
        except requests.RequestException as ex:
            err = type(ex).__name__
            continue
        if r.status_code == 200 and "@" in r.text:
            r.encoding = "utf-8"
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


def venue_of(e):
    v = (VENUE_CONFIRMED.get(e["key"]) or e.get("venue_now") or e.get("venue") or "").strip()
    if not v or re.match(r"(?i)arxiv", v) or re.search(r"(?i)arxiv only", v):
        return None
    return v


def braced_value(bib, field):
    """Return (start, end) of the brace-delimited value of field, or None."""
    m = re.search(r"(?i)(?<![a-z])" + field + r"\s*=\s*\{", bib)
    if not m:
        return None
    i, depth = m.end(), 1
    while i < len(bib) and depth:
        depth += {"{": 1, "}": -1}.get(bib[i], 0)
        i += 1
    return m.end(), i - 1


def protect_title(bib):
    span = braced_value(bib, "title")
    if not span:
        return bib
    a, b = span
    inner = bib[a:b].strip()
    if inner.startswith("{") and inner.endswith("}"):
        return bib
    return bib[:a] + "{" + inner + "}" + bib[b:]


def set_field(bib, field, value):
    span = braced_value(bib, field)
    if span:
        a, b = span
        return bib[:a] + value + bib[b:]
    # insert after the first line
    first_nl = bib.index("\n") + 1 if "\n" in bib else len(bib)
    return bib[:first_nl] + f"  {field} = {{{value}}},\n" + bib[first_nl:]


def convert(bib, e, kind):
    v = venue_of(e)
    comment = ""
    if not v:
        return bib, comment
    years = re.findall(r"\b(?:19|20)\d\d\b", v)
    y = years[0] if years else None
    head = re.split(r";|, earlier|\(the cited|\(also", v)[0]   # the venue of record, before any history
    if re.search(r"(?i)workshop", head):
        if kind == "arxiv":
            bib = set_field(bib, "note", v.split(";")[0].split(" per ")[0].strip())
        return bib, f"% {e['key']}: workshop venue recorded in the catalogue: {v}\n"
    if kind != "arxiv":                      # ACL Anthology and Crossref records carry their venue
        return bib, comment
    for pat, typ, tmpl in VENUES:
        if re.search(pat, head):
            bib = re.sub(r"^@\w+\{", f"@{typ}{{", bib, count=1)
            field = "booktitle" if typ == "inproceedings" else "journal"
            bib = set_field(bib, field, tmpl.format(y=y or ""))
            if y:
                bib = set_field(bib, "year", y)
            flag = " (single source; confirm before submission)" if e["key"] in SINGLE_SOURCE else ""
            origin = "confirmed after the search" if e["key"] in VENUE_CONFIRMED else "from the catalogue"
            comment = f"% {e['key']}: venue {origin}: {v}{flag}\n"
            return bib, comment
    return bib, f"% {e['key']}: venue in the catalogue not mapped, kept as a preprint: {v}\n"


def rekey(bib, key):
    return re.sub(r"^@(\w+)\{[^,]*,", lambda m: f"@{m.group(1)}{{{key},", bib.strip(), count=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--keys", help="file with one catalogue key per line")
    ap.add_argument("--out", default=os.path.join(HERE, "related_work_sources.bib"))
    args = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    entries = load_entries()
    wanted = None
    if args.only:
        wanted = list(args.only)
    if args.keys:
        with open(args.keys, encoding="utf-8") as fh:
            wanted = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    by_key = {e["key"]: e for e in entries}
    if wanted is not None:
        missing = [k for k in wanted if k not in by_key]
        if missing:
            sys.exit(f"keys not in any catalogue: {missing}")
        entries = [by_key[k] for k in wanted]
    records, problems, last_arxiv = [], [], 0.0
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
        bib = rekey(text, e.get("bibkey") or e["key"])
        if e["key"] in AUTHOR_FIXES:
            bib = bib.replace(*AUTHOR_FIXES[e["key"]])
        if kind in ("arxiv", "doi"):
            bib = protect_title(bib)
        bib, comment = convert(bib, e, kind)
        records.append(comment + bib.strip() + "\n")
        print(f"ok       {e.get('bibkey') or e['key']}  ({kind})")
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("% Built by docs/related_work_oct2026/make_bib.py from arXiv, ACL Anthology and Crossref records.\n")
        fh.write("% Authors and titles come from the source. Published venues come from the catalogue (papers.json\n")
        fh.write("% venue_now, candidates/*.json venue). Check page numbers against the proceedings before submission.\n\n")
        fh.write("\n".join(records))
    print(f"\n{len(records)} entries written to {args.out}")
    for k, p in problems:
        print(f"PROBLEM  {k}: {p}")


if __name__ == "__main__":
    main()
