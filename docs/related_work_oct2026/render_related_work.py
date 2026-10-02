#!/usr/bin/env python
"""Render the revised Related Work from its single source into Markdown and LaTeX.

Source: related_work_v3.src.md, with pandoc-style citations: [@a; @b] for parenthetical, bare @a for
textual, and parts marked <!-- part: md-header | body | appendix | md-notes -->. Outputs:
  RELATED_WORK_v3.md   md-header + body + appendix + md-notes, citation labels computed from the .bib
                       the way natbib (ACL style) prints them, a/b suffixes included
  related_work_v3.tex  body + appendix as LaTeX: \\citep / \\citet, \\paragraph, a table* for Table A
Checks: every cited key is in related_work_sources.bib and related_work_keys.txt, every bib entry is
cited, and prints the body's word count with and without the [cut first] sentences.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CITE_P = re.compile(r"\[(@[A-Za-z0-9_]+(?:;\s*@[A-Za-z0-9_]+)*)\]")
CITE_T = re.compile(r"(?<![\w\[;])@([A-Za-z0-9_]+)")


def parse_bib(path):
    text = open(path, encoding="utf-8").read()
    entries = {}
    for m in re.finditer(r"@(\w+)\{([^,\s]+),", text):
        key, i = m.group(2), m.end()
        depth, j = 1, i
        while j < len(text) and depth:
            depth += {"{": 1, "}": -1}.get(text[j], 0)
            j += 1
        body = text[i:j - 1]
        fields = {}
        for fm in re.finditer(r"(?i)(?<![a-z])(author|year|title)\s*=\s*", body):
            k, s = fm.group(1).lower(), fm.end()
            if s < len(body) and body[s] == "{":
                d, t = 1, s + 1
                while t < len(body) and d:
                    d += {"{": 1, "}": -1}.get(body[t], 0)
                    t += 1
                fields.setdefault(k, body[s + 1:t - 1])
            elif s < len(body) and body[s] == '"':
                t = body.index('"', s + 1)
                fields.setdefault(k, body[s + 1:t])
            else:
                v = re.match(r"[^,}\s]+", body[s:])
                fields.setdefault(k, v.group(0) if v else "")
        entries[key] = fields
    return entries


def surname(author):
    a = re.sub(r"\s+", " ", author.strip())
    if a.startswith("{") and a.endswith("}"):
        return a[1:-1]
    if "," in a:
        return a.split(",")[0].strip().strip("{}")
    return a.split(" ")[-1].strip("{}")


def labels(entries):
    info = {}
    for k, f in entries.items():
        authors = [x for x in re.split(r"\s+and\s+", f.get("author", "")) if x.strip()]
        sn = [surname(a) for a in authors]
        lab = sn[0] if len(sn) == 1 else (f"{sn[0]} and {sn[1]}" if len(sn) == 2 else f"{sn[0]} et al.")
        sortkey = " ".join((surname(a) + " " + a) for a in authors) + " " + f.get("title", "")
        info[k] = {"label": lab, "year": f.get("year", "?"), "sort": sortkey.lower()}
    groups = {}
    for k, v in info.items():
        groups.setdefault((v["label"], v["year"]), []).append(k)
    for (_, _), ks in groups.items():
        if len(ks) > 1:
            for i, k in enumerate(sorted(ks, key=lambda x: info[x]["sort"])):
                info[k]["year"] += "abcdefgh"[i]
    return info


def split_parts(src):
    parts, name, buf = {}, None, []
    for line in src.splitlines():
        m = re.match(r"<!-- part: ([\w-]+) -->", line.strip())
        if m:
            if name:
                parts[name] = "\n".join(buf).strip("\n")
            name, buf = m.group(1), []
        else:
            buf.append(line)
    if name:
        parts[name] = "\n".join(buf).strip("\n")
    return parts


def render_md(text, info):
    def p(m):
        keys = [k.strip()[1:] for k in m.group(1).split(";")]
        return "(" + "; ".join(f"{info[k]['label']}, {info[k]['year']}" for k in keys) + ")"
    text = CITE_P.sub(p, text)
    text = CITE_T.sub(lambda m: f"{info[m.group(1)]['label']} ({info[m.group(1)]['year']})", text)
    return text.replace("<!-- cut first -->", "[cut first] ")


def tex_escape(s):
    s = s.replace("\\", "\\textbackslash{}").replace("%", "\\%").replace("&", "\\&").replace("#", "\\#")
    s = s.replace("_", "\\_")
    s = re.sub(r'"([^"]+)"', r"``\1''", s)
    s = s.replace("τb", "$\\tau_b$").replace("κ", "$\\kappa$").replace("ρ", "$\\rho$")
    s = s.replace("–", "--")
    s = re.sub(r"\*\*([^*]+)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"\\textit{\1}", s)
    return s


def render_tex_text(text):
    holders = []

    def hold(cmd):
        holders.append(cmd)
        return f"\x00{len(holders) - 1}\x00"
    text = CITE_P.sub(lambda m: hold("\\citep{" + ",".join(k.strip()[1:] for k in m.group(1).split(";")) + "}"), text)
    text = CITE_T.sub(lambda m: hold("\\citet{" + m.group(1) + "}"), text)
    text = tex_escape(text)
    return re.sub(r"\x00(\d+)\x00", lambda m: holders[int(m.group(1))], text)


def render_tex(part):
    out, lines, i = [], part.splitlines(), 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("## 2 Related Work"):
            out += ["\\section{Related Work}", "\\label{sec:related}"]
        elif line.startswith("## Appendix"):
            out += ["\\section{Extended Related Work}", "\\label{app:related}"]
        elif line.startswith("**Table A.**"):
            caption = render_tex_text(line.replace("**Table A.**", "").strip())
            rows = []
            i += 1
            while i < len(lines) and (lines[i].startswith("|") or not lines[i].strip()):
                if lines[i].startswith("|") and not re.match(r"\|\s*-", lines[i]):
                    rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            widths = ["2.5cm", "1.9cm", "2.3cm", "2.1cm", "3.4cm", "3.3cm"]
            out += ["\\begin{table*}[t]", "\\centering", "\\footnotesize",
                    "\\begin{tabular}{" + "".join(f"p{{{w}}}" for w in widths) + "}", "\\toprule"]
            out.append(" & ".join("\\textbf{" + render_tex_text(c) + "}" for c in rows[0]) + " \\\\")
            out.append("\\midrule")
            for r in rows[1:]:
                out.append(" & ".join(render_tex_text(c) for c in r) + " \\\\")
            out += ["\\bottomrule", "\\end{tabular}",
                    "\\caption{" + caption + "}", "\\label{tab:related}", "\\end{table*}"]
            continue
        elif re.match(r"\*\*[^*]+\*\*\s*$", line.strip()):
            out.append("\\paragraph{" + render_tex_text(line.strip()[2:-2].rstrip(".")) + ".}")
        elif "<!-- cut first -->" in line:
            out.append("% [cut first] the next sentence can be removed if space is short")
            out.append(render_tex_text(line.replace("<!-- cut first -->", "")))
        else:
            out.append(render_tex_text(line))
        i += 1
    return "\n".join(out)


def words(text):
    t = CITE_P.sub("", text)
    t = re.sub(r"[#*|]", " ", t)
    return len(re.findall(r"[A-Za-z0-9][\w'.,%-]*", t))


def main():
    src = open(os.path.join(HERE, "related_work_v3.src.md"), encoding="utf-8").read()
    parts = split_parts(src)
    entries = parse_bib(os.path.join(HERE, "related_work_sources.bib"))
    info = labels(entries)
    cited = set()
    for part in ("body", "appendix"):
        for m in CITE_P.finditer(parts[part]):
            cited |= {k.strip()[1:] for k in m.group(1).split(";")}
        for m in CITE_T.finditer(CITE_P.sub("", parts[part])):
            cited.add(m.group(1))
    keyfile = {ln.strip() for ln in open(os.path.join(HERE, "related_work_keys.txt"), encoding="utf-8")
               if ln.strip() and not ln.startswith("#")}
    keyfile_bib = set()
    for k in keyfile:   # the key list holds catalogue keys; two are re-keyed in the .bib
        keyfile_bib.add({"cheng2024elecbench": "zhou2024elecbench", "chen2025apbench": "wu2025apbench"}.get(k, k))
    problems = []
    problems += [f"cited but not in the .bib: {k}" for k in sorted(cited - set(entries))]
    problems += [f"in the .bib but not cited: {k}" for k in sorted(set(entries) - cited)]
    problems += [f"cited but not in related_work_keys.txt: {k}" for k in sorted(cited - keyfile_bib)]
    if problems:
        print("\n".join(problems))
        sys.exit(1)
    md = "\n\n".join([parts["md-header"], render_md(parts["body"], info),
                      render_md(parts["appendix"], info), parts["md-notes"]]) + "\n"
    with open(os.path.join(HERE, "RELATED_WORK_v3.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    tex = ["% EngTrace, Related Work revised for the October 2026 submission.",
           "% Rendered by docs/related_work_oct2026/render_related_work.py from related_work_v3.src.md; edit the source.",
           "% Needs natbib (ACL style), booktabs, and the entries in related_work_sources.bib.",
           "% Every statement about another paper is checked by verify_facts.py (notes/fact_check.md).", "",
           render_tex(parts["body"]), "", "% ---- Appendix ----", "", render_tex(parts["appendix"]), ""]
    with open(os.path.join(HERE, "related_work_v3.tex"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(tex))
    body = parts["body"]
    cut = "\n".join(re.sub(r"<!-- cut first -->.*", "", ln) for ln in body.splitlines())
    print(f"{len(cited)} keys cited, all in the .bib; body {words(body)} words, {words(cut)} without the [cut first] sentences")
    may = may_words()
    if may:
        print(f"May 2026 section 2, same count: {may} words")


def may_words():
    """Words of prose in the May submission's section 2, citations removed, counted like words()."""
    pdf = os.path.join(HERE, "..", "_ARR_May__EngTrace.pdf")
    try:
        import fitz
    except ImportError:
        return None
    with fitz.open(pdf) as doc:
        text = "\n".join(doc[i].get_text() for i in range(min(4, doc.page_count)))
    m = re.search(r"\n2\s+Related Work\n(.*?)\n3\s*EngTrace", text, flags=re.S)
    if not m:
        return None
    t = re.sub(r"-\n", "", m.group(1))
    t = re.sub(r"\([^()]*\d{4}[a-z]?[^()]*\)", " ", t)    # parenthetical citations
    return words(t)


if __name__ == "__main__":
    main()
