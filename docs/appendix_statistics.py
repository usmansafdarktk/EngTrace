"""Compute the numbers of the appendix's Dataset Statistics section and check them against the .tex.

    python docs/appendix_statistics.py                    # print the tables' LaTeX rows, the numbers the text uses
                                                          # and the two generated blocks
    python docs/appendix_statistics.py --check            # exit 1 unless appendices/taxonomy_content.tex holds every
                                                          # generated row and number, and every placed block is current
    python docs/appendix_statistics.py --write --out DIR  # fill the marked blocks in DIR/appendices/taxonomy_content.tex
                                                          # (a copy of the tex tree; the file is copied there if absent;
                                                          # a block without markers is appended at the end of the file)
    python docs/appendix_statistics.py --write            # the same in overleaf_source_04102026/ (Phase 2; markers required)
    python docs/appendix_statistics.py --selftest         # the area table sums to 150 and 30 per branch; both blocks render
    python docs/appendix_statistics.py --agreement PATH   # read the difficulty-agreement figures from PATH instead

Sources: full_run_28092026/manifest.jsonl (templates, domains, areas, levels, answer kinds, instances),
full_run_28092026/diversity.json (question wording, reasoning paths and answer forms, near-duplicates),
full_run_28092026/FREEZE.json (how the evaluation set was drawn) and, for `tab:levels_agreement`,
template_annotation_23092026/levels/agreement.json (score_levels.py's output; until it exists the block is a
stand-in with empty cells, marked STAND-IN in its caption). Nothing is typed by hand.

Generated blocks (the BEGIN/END convention of paper_results.py, the marker naming this script):
  tab:area              templates per area and level, grouped by domain, with branch subtotals; two panels
  tab:levels_agreement  the experts' difficulty ratings: Fleiss' kappa per branch and overall with 95% intervals,
                        the linear- and quadratic-weighted forms, the templates whose two-of-three majority equals
                        the label or lies one or two steps from it, and the formula-stated count by level
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "full_run_28092026"
SRC = ROOT / "overleaf_source_04102026"
TEX_REL = Path("appendices") / "taxonomy_content.tex"
TEX = SRC / TEX_REL
AGREEMENT = ROOT / "template_annotation_23092026" / "levels" / "agreement.json"
BRANCHES = ["chemical", "electrical", "civil", "industrial", "mechanical"]  # the taxonomy figure's order
LEVELS = ["Easy", "Intermediate", "Advanced"]
KINDS = ["scalar", "multipart", "symbolic", "vector", "array", "classification"]
MARK = "% BEGIN GENERATED {name} (docs/appendix_statistics.py --write)\n{body}\n% END GENERATED {name}"
AREA_NAMES = {  # where the generic rule (words capitalised, "and"/"of" lower) misses the figure's spelling
    "volumetric_properties_pure_fluids": "Volumetric Properties of Pure Fluids",
    "permeability_seepage_effective_stress": "Permeability, Seepage, and Effective Stress",
}

rows = [json.loads(l) for l in (RUN / "manifest.jsonl").read_text(encoding="utf-8").splitlines() if l.strip()]
div = json.loads((RUN / "diversity.json").read_text(encoding="utf-8"))
freeze = json.loads((RUN / "FREEZE.json").read_text(encoding="utf-8"))

# ---------------------------------------------------------------- templates by branch: structure and level
tmpl = {}
for r in rows:
    tmpl[r["template_id"]] = (r["branch"].split("_")[0], r["domain"], r["area"], r["level"], r["answer_type"])
branch_rows, tot = [], Counter()
for b in BRANCHES:
    ts = [v for v in tmpl.values() if v[0] == b]
    doms, areas = {v[1] for v in ts}, {(v[1], v[2]) for v in ts}
    lv = Counter(v[3] for v in ts)
    cells = [len(doms), len(areas), lv["Easy"], lv["Intermediate"], lv["Advanced"], len(ts)]
    tot.update(dict(zip(["dom", "area", "e", "i", "a", "n"], cells)))
    branch_rows.append(f"{b.capitalize()} & {cells[0]} & {cells[1]} & {cells[2]} & {cells[3]} & {cells[4]} & \\textbf{{{cells[5]}}} \\\\")
branch_rows.append(f"\\textbf{{Total}} & \\textbf{{{tot['dom']}}} & \\textbf{{{tot['area']}}} & \\textbf{{{tot['e']}}} & "
                   f"\\textbf{{{tot['i']}}} & \\textbf{{{tot['a']}}} & \\textbf{{{tot['n']}}} \\\\")
areas_per_branch = {b: len({(v[1], v[2]) for v in tmpl.values() if v[0] == b}) for b in BRANCHES}

# ---------------------------------------------------------------- instances
inst_level = Counter(r["level"] for r in rows)
inst_branch = Counter(r["branch"].split("_")[0] for r in rows)
per_template = Counter(r["template_id"] for r in rows)

# ---------------------------------------------------------------- answer kinds
kind_t = Counter(v[4] for v in tmpl.values())

# ---------------------------------------------------------------- variation within a template (evaluation set)
def bucket(n: int) -> str:
    return "1" if n == 1 else "2" if n == 2 else "3--5" if n <= 5 else "6--15"

var = {k: Counter(bucket(t["pool"][k]) for t in div) for k in ("question_skeletons", "paths_lower", "answer_variants")}
variation_rows = [f"{b} & {var['question_skeletons'][b]} & {var['paths_lower'][b]} & {var['answer_variants'][b]} \\\\"
                  for b in ("1", "2", "3--5", "6--15")]
multi = {t["template_id"] for t in div if t["pool"]["paths_lower"] >= 2}
multi_reach = {t["template_id"] for t in div if t["reach"]["paths_lower"] >= 2}
multi_upper = sum(1 for t in div if t["pool"]["paths_upper"] >= 2)
single_by_branch = Counter(t["branch"].split("_")[0] for t in div if t["pool"]["paths_lower"] == 1)
near5 = sum(t["pool"]["near_duplicate_pairs"]["5%"] for t in div)
near5_templates = sum(1 for t in div if t["pool"]["near_duplicate_pairs"]["5%"] > 0)

# ---------------------------------------------------------------- how the evaluation set was drawn
cov = freeze["coverage"]["templates"]
rej = freeze["walk_rejections"]["by_reason"]

numbers = {
    "templates": len(tmpl), "instances": len(rows), "instances per template": sorted(set(per_template.values())),
    "instances Easy/Intermediate/Advanced": [inst_level[l] for l in LEVELS],
    "instances per branch": sorted(set(inst_branch.values())),
    "Easy+Intermediate share of templates": round(100 * (tot["e"] + tot["i"]) / tot["n"]),
    "Advanced share of templates": round(100 * tot["a"] / tot["n"]),
    "areas per branch": areas_per_branch,
    "answer kinds (templates)": {k: kind_t[k] for k in KINDS},
    "templates with several reasoning paths (evaluation set, lower reading)": len(multi),
    "same set in 500 draws": multi == multi_reach,
    "templates with several reasoning paths (upper reading)": multi_upper,
    "single-path templates by branch": {b: single_by_branch[b] for b in BRANCHES},
    "near-duplicate pairs at 5%": near5, "templates with such a pair": near5_templates,
    "draws examined per template": 100, "draws examined": 100 * len(cov),
    "candidates": sum(t["candidates"] for t in cov),
    "skipped: repeated question / rounding tie / validation-set question":
        [rej.get("duplicate question", 0), rej.get("display tie", 0), rej.get("pilot-slice question", 0)],
    "templates with more than one group": sum(1 for t in cov if t["groups"] > 1),
    "templates with more groups than instances selected": sum(1 for t in cov if t["groups_selected"] < t["groups"]),
    "seed bits": freeze["seed"]["bits"], "distinct questions": freeze["distinct_questions"],
}

# The numbers the section's prose states, as they must appear in the .tex.
prose = [
    f"{tot['e']} Easy, {tot['i']} Intermediate, and {tot['a']} Advanced",
    f"{numbers['Easy+Intermediate share of templates']}\\%", f"{numbers['Advanced share of templates']}\\%",
    f"{inst_level['Easy']:,} Easy, {inst_level['Intermediate']:,} Intermediate, and {inst_level['Advanced']:,} Advanced",
    f"{inst_branch['chemical']} per branch",
    f"{kind_t['scalar']} scalar, {kind_t['multipart']} multipart, {kind_t['symbolic']} symbolic, {kind_t['vector']} vector,",
    f"{kind_t['array']} array, and {kind_t['classification']} classification",
    f"{len(multi)} of the {len(tmpl)} templates", f"{len(div) - len(multi)} follow a single path",
    f"(chemical {single_by_branch['chemical']}, electrical {single_by_branch['electrical']}, civil "
    f"{single_by_branch['civil']}, industrial {single_by_branch['industrial']}, mechanical {single_by_branch['mechanical']})",
    f"{multi_upper} templates", f"{near5} such pairs, in {near5_templates} templates",
    f"{numbers['draws examined']:,} draws", f"{numbers['candidates']:,} remain",
    f"{numbers['skipped: repeated question / rounding tie / validation-set question'][0]} repeated questions",
    f"{numbers['skipped: repeated question / rounding tie / validation-set question'][1]} rounding ties",
    f"{numbers['templates with more than one group']} templates",
    f"{numbers['templates with more groups than instances selected']} templates",
]

# The depth sentence after the domain prompt: three domains per branch, each with templates at every level.
dom_levels = {}
for b, d, _, lv, _ in tmpl.values():
    dom_levels.setdefault((b, d), set()).add(lv)
three_per_branch = all(sum(1 for bb, _ in dom_levels if bb == b) == 3 for b in BRANCHES)
every_level = all(levels == set(LEVELS) for levels in dom_levels.values())
prose.append("We keep three domains per branch to cover each in depth, across all its areas and all three levels"
             if three_per_branch and every_level else "[depth sentence after the domain prompt no longer holds]")


# ================================================================ generated blocks
def pretty(snake: str) -> str:
    return AREA_NAMES.get(snake, " ".join(w if w in ("and", "of") else w.capitalize() for w in snake.split("_")))


def caption_lines(caption: str) -> str:
    """The caption wrapped at 100 characters (paper_results.py's form)."""
    lines = textwrap.wrap(caption, width=100, initial_indent=" " * 9, break_long_words=False, break_on_hyphens=False)
    return "\\caption{" + "\n".join(lines)[9:] + "}"


def block(name: str, body: str) -> str:
    return MARK.format(name=name, body=body.rstrip("\n"))


def area_counts() -> dict:
    """{branch: {domain: {area: Counter(level)}}} over the templates, from the manifest's area field."""
    out: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))
    for b, d, a, lv, _k in tmpl.values():
        out[b][d][a][lv] += 1
    return out


def area_panel_rows(branches: list[str], counts: dict) -> list[str]:
    L = []
    for b in branches:
        doms = counts[b]
        n_areas = sum(len(v) for v in doms.values())
        sub = Counter()
        for areas in doms.values():
            for c in areas.values():
                sub.update(c)
        L.append(f"\\rowcolor{{gray!10}}\\textbf{{{b.capitalize()}}} ({n_areas} areas) & \\textbf{{{sub['Easy']}}} & "
                 f"\\textbf{{{sub['Intermediate']}}} & \\textbf{{{sub['Advanced']}}} & \\textbf{{{sum(sub.values())}}} \\\\")
        for d in sorted(doms, key=pretty):
            L.append(f"\\multicolumn{{5}}{{@{{}}l}}{{\\textit{{{pretty(d)}}}}} \\\\")
            for a in sorted(doms[d], key=pretty):
                c = doms[d][a]
                L.append(f"\\quad {pretty(a)} & {c['Easy']} & {c['Intermediate']} & {c['Advanced']} & {sum(c.values())} \\\\")
    return L


def table_area() -> str:
    counts = area_counts()
    # two panels, whole branches, the split nearest to half the rows (branch, domain and area rows)
    sizes = [1 + len(counts[b]) + sum(len(v) for v in counts[b].values()) for b in BRANCHES]
    half, best, cum = sum(sizes) / 2, None, 0
    for k, s in enumerate(sizes, 1):
        cum += s
        if best is None or abs(cum - half) < abs(best[1] - half):
            best = (k, cum)
    left, right = BRANCHES[:best[0]], BRANCHES[best[0]:]
    total = Counter()
    for b in BRANCHES:
        for areas in counts[b].values():
            for c in areas.values():
                total.update(c)
    n_dom = sum(len(counts[b]) for b in BRANCHES)
    n_area = sum(len(v) for b in BRANCHES for v in counts[b].values())
    head = ("\\toprule\n\\rowcolor{gray!10}\\textbf{Domain and area} & \\textbf{Easy} & \\textbf{Int.} & \\textbf{Adv.} & "
            "\\textbf{Total} \\\\\n\\midrule")

    def panel(branches: list[str], last: bool) -> list[str]:
        L = ["\\begin{minipage}[t]{0.48\\textwidth}", "\\centering",
             "\\begin{tabular}{@{}p{118pt} r r r r@{}}", head] + area_panel_rows(branches, counts)
        if last:
            L.append("\\midrule")
            L.append(f"\\rowcolor{{gray!10}}\\textbf{{Total}} & \\textbf{{{total['Easy']}}} & "
                     f"\\textbf{{{total['Intermediate']}}} & \\textbf{{{total['Advanced']}}} & \\textbf{{{sum(total.values())}}} \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{minipage}"]
        return L

    caption = (f"\\textbf{{Templates per area and level.}} {n_area} areas in {n_dom} domains; each branch has "
               f"{sum(total.values()) // len(BRANCHES)} templates; the shaded rows give the branch subtotals and the "
               f"italic rows name the domain.")
    lines = (["\\begin{table*}[t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.05}",
              "\\setlength{\\tabcolsep}{3pt}"] + panel(left, False) + ["\\hfill"] + panel(right, True)
             + [caption_lines(caption), "\\label{tab:area}", "\\end{table*}"])
    return "\n".join(lines)


def kappa_cell(k, ci=None) -> str:
    if k is None:
        return "--"
    s = f"{k:.2f}"
    return s + (f" [{ci[0]:.2f}, {ci[1]:.2f}]" if ci else "")


def table_levels_agreement(res: dict | None) -> str:
    """The difficulty-agreement table from score_levels.py's agreement.json; a stand-in with empty cells when None."""
    L = []
    names = [(b, f"{b}_engineering") for b in BRANCHES] + [("All branches", None)]
    for label, key in names:
        name = label.capitalize() if key else f"\\textbf{{{label}}}"
        if res is None:
            L.append(f"{name} & " + " & ".join(["--"] * 11) + " \\\\")
            continue
        a = res["agreement"]["by_branch"][key] if key else res["agreement"]["overall"]
        m = res["majority"]["by_branch"][key] if key else res["majority"]["overall"]
        f = res["formula_stated"]["by_branch_current_level"][key] if key else res["formula_stated"]["by_current_level"]
        cells = [str(a["templates"]), kappa_cell(a["fleiss"], a.get("fleiss_ci95")),
                 kappa_cell(a["weighted_linear"], a.get("weighted_linear_ci95")), kappa_cell(a["weighted_quadratic"]),
                 str(m["same_as_current"]), str(m["differs_by_one_step"]), str(m["differs_by_two_steps"]),
                 str(m["no_majority"])] + [f"{f[lv]['yes']}/{f[lv]['templates']}" for lv in LEVELS]
        if not key:
            cells = [f"\\textbf{{{c}}}" for c in cells]
        L.append(f"{name} & " + " & ".join(cells) + " \\\\")
    header = ["\\toprule",
              "\\rowcolor{gray!10}\\textbf{Branch} & \\textbf{Templates} & \\multicolumn{3}{c}{\\textbf{Agreement among the three raters}} & "
              "\\multicolumn{4}{c}{\\textbf{Majority against the label}} & \\multicolumn{3}{c}{\\textbf{Formula stated}} \\\\",
              "\\rowcolor{gray!10} & & \\textbf{Fleiss' $\\kappa$ [95\\% CI]} & \\textbf{Linear [95\\% CI]} & \\textbf{Quadr.} & "
              "\\textbf{Same} & \\textbf{1 step} & \\textbf{2 steps} & \\textbf{None} & \\textbf{Easy} & \\textbf{Int.} & "
              "\\textbf{Adv.} \\\\", "\\midrule"]
    rows_ = header + L[:-1] + ["\\midrule", L[-1], "\\bottomrule"]
    caption = ("The domain experts' difficulty ratings, three per template, against each other and against the labels: "
               "Fleiss' $\\kappa$ on the three levels with its 95\\% percentile bootstrap interval over templates, the "
               "weighted forms for ordered levels (adjacent levels count 0.5 under linear and 0.75 under quadratic "
               "weights), the templates whose two-of-three majority equals the label, lies one or two steps from it, "
               "or does not exist, and the templates whose majority says the problem states its governing formula or "
               "names the method, as yes counts over the labeled templates of each level.")
    if res is None:
        caption = ("STAND-IN: the experts' ratings have not been returned; every cell fills from "
                   "\\texttt{levels/agreement.json}. " + caption)
    lines = (["\\begin{table*}[t]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.1}",
              "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{l r l l r r r r r r r r}"] + rows_
             + ["\\end{tabular}", caption_lines(caption), "\\label{tab:levels_agreement}", "\\end{table*}"])
    return "\n".join(lines)


def load_agreement(path: Path | None) -> dict | None:
    p = path or AGREEMENT
    if not p.exists():
        return None
    res = json.loads(p.read_text(encoding="utf-8"))
    if not res.get("returns", {}).get("ratings"):
        return None
    return res


def blocks(agreement: dict | None) -> dict[str, str]:
    """The generated blocks: `tab:area` always; `tab:levels_agreement` only once agreement.json holds ratings. D7 keeps
    the domain experts' labels and states no agreement figure until the experts' count files arrive, so no stand-in is
    placed in the tree meanwhile."""
    out = {"tab:area": table_area()}
    if agreement:
        out["tab:levels_agreement"] = table_levels_agreement(agreement)
    return out


def pattern(name: str) -> re.Pattern:
    return re.compile(r"% BEGIN GENERATED " + re.escape(name) + r" .*?% END GENERATED " + re.escape(name), re.S)


def flat(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def write(out: Path | None, agreement: dict | None) -> None:
    target = (out / TEX_REL) if out else TEX
    if out and not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(TEX, target)
    raw = target.read_bytes()
    crlf = b"\r\n" in raw
    tex = raw.decode("utf-8").replace("\r\n", "\n")
    for name, body in blocks(agreement).items():
        if pattern(name).search(tex):
            tex = pattern(name).sub(lambda m: block(name, body), tex, count=1)
            print(f"{name}: replaced between its markers in {target}")
        elif out:
            tex = tex.rstrip("\n") + "\n\n" + block(name, body) + "\n"
            print(f"{name}: no markers in {target.name}; appended at the end of the file (Phase 2 places the markers)")
        else:
            raise SystemExit(f"{target} has no markers for {name}; --write without --out needs them")
    target.write_bytes(tex.replace("\n", "\r\n").encode("utf-8") if crlf else tex.encode("utf-8"))


def check(agreement: dict | None) -> int:
    tex_raw = TEX.read_text(encoding="utf-8")
    tex = " ".join(tex_raw.split())
    missing = [s for s in branch_rows + variation_rows + prose if " ".join(s.split()) not in tex]
    for s in missing:
        print("MISSING:", s)
    print(f"{len(branch_rows) + len(variation_rows) + len(prose) - len(missing)} of "
          f"{len(branch_rows) + len(variation_rows) + len(prose)} generated rows and numbers are in {TEX.name}")
    stale = 0
    for name, body in blocks(agreement).items():
        m = pattern(name).search(tex_raw)
        if not m:
            print(f"{name}: no markers in {TEX.name} yet (placed in Phase 2)")
        elif flat(m.group(0)) != flat(block(name, body)):
            print(f"STALE block in {TEX.name}: {name}")
            stale += 1
        else:
            print(f"{name}: block current")
    return 1 if missing or stale else 0


def selftest() -> int:
    counts = area_counts()
    per_branch = {b: sum(sum(c.values()) for areas in counts[b].values() for c in areas.values()) for b in BRANCHES}
    assert sum(per_branch.values()) == 150 == len(tmpl), per_branch
    assert all(n == 30 for n in per_branch.values()), per_branch
    assert sum(len(v) for b in BRANCHES for v in counts[b].values()) == tot["area"] == 42
    assert sum(len(counts[b]) for b in BRANCHES) == tot["dom"] == 15
    area = table_area()
    # every area row's total equals its level cells; the panel subtotals equal the branch table
    for line in area.splitlines():
        m = re.match(r"\\quad .* & (\d+) & (\d+) & (\d+) & (\d+) \\\\$", line)
        if m:
            e, i, a, n = map(int, m.groups())
            assert e + i + a == n, line
    for b in BRANCHES:
        ts = [v for v in tmpl.values() if v[0] == b]
        lv = Counter(v[3] for v in ts)
        assert (f"\\textbf{{{b.capitalize()}}} ({areas_per_branch[b]} areas) & \\textbf{{{lv['Easy']}}} & "
                f"\\textbf{{{lv['Intermediate']}}} & \\textbf{{{lv['Advanced']}}} & \\textbf{{30}}") in area, b
    assert f"\\textbf{{{tot['e']}}} & \\textbf{{{tot['i']}}} & \\textbf{{{tot['a']}}} & \\textbf{{150}}" in area
    assert area.count("\\begin{minipage}") == 2 and area.count("\\quad ") == 42
    # the agreement table renders as a stand-in and from a synthetic result in score_levels.py's schema
    stand_in = table_levels_agreement(None)
    empty = sum(l.count(" -- ") + l.count(" --") * 0 for l in stand_in.splitlines() if l.endswith("\\\\") and " & " in l)
    assert "STAND-IN" in stand_in and empty == 6 * 11, empty
    synth = {"returns": {"ratings": 450},
             "agreement": {"overall": {"templates": 150, "fleiss": 0.512, "fleiss_ci95": [0.43, 0.59], "weighted_linear": 0.6,
                                       "weighted_linear_ci95": [0.52, 0.68], "weighted_quadratic": 0.66},
                           "by_branch": {f"{b}_engineering": {"templates": 30, "fleiss": 0.5, "fleiss_ci95": [0.3, 0.7],
                                                              "weighted_linear": 0.6, "weighted_linear_ci95": [0.4, 0.8],
                                                              "weighted_quadratic": None} for b in BRANCHES}},
             "majority": {"overall": {"same_as_current": 120, "differs_by_one_step": 25, "differs_by_two_steps": 1, "no_majority": 4},
                          "by_branch": {f"{b}_engineering": {"same_as_current": 24, "differs_by_one_step": 5,
                                                             "differs_by_two_steps": 0, "no_majority": 1} for b in BRANCHES}},
             "formula_stated": {"by_current_level": {lv: {"yes": 10, "templates": 58} for lv in LEVELS},
                                "by_branch_current_level": {f"{b}_engineering": {lv: {"yes": 2, "templates": 12} for lv in LEVELS}
                                                            for b in BRANCHES}}}
    real = table_levels_agreement(synth)
    assert "STAND-IN" not in real and "0.51 [0.43, 0.59]" in real and "\\textbf{10/58}" in real and "& -- &" in real
    assert list(blocks(None)) == ["tab:area"] and list(blocks(synth)) == ["tab:area", "tab:levels_agreement"]
    for name, body in blocks(synth).items():
        assert block(name, body).count("% BEGIN GENERATED") == 1 and body.count(f"\\label{{{name}}}") == 1
        cap = body[body.index("\\caption{"):body.index("\\label{")]
        assert all(len(l) <= 100 for l in cap.splitlines()), cap       # captions wrap at 100; table rows may exceed it
    # a write into a scratch copy appends tab:area alone while there are no ratings, the agreement block once there
    # are, and a second write replaces each in place
    import tempfile
    with tempfile.TemporaryDirectory(prefix="engtrace-stats-") as tmp:
        out = Path(tmp)
        write(out, None)
        once = (out / TEX_REL).read_text(encoding="utf-8")
        write(out, synth)
        twice = (out / TEX_REL).read_text(encoding="utf-8")
        write(out, synth)
        thrice = (out / TEX_REL).read_text(encoding="utf-8")
        assert once.count("% BEGIN GENERATED tab:area") == twice.count("% BEGIN GENERATED tab:area") == 1
        assert "tab:levels_agreement" not in once and "STAND-IN" not in once
        assert twice.count("% BEGIN GENERATED tab:levels_agreement") == 1 == thrice.count("% BEGIN GENERATED tab:levels_agreement")
        assert "STAND-IN" not in twice and twice == thrice
        assert once.split("% BEGIN GENERATED")[0] == TEX.read_text(encoding="utf-8").rstrip("\n") + "\n\n"
    print("SELFTEST OK: 150 templates, 30 per branch, 15 domains, 42 areas in tab:area; tab:levels_agreement renders as a "
          "stand-in and from a synthetic agreement.json, and is written only once ratings exist; --write --out appends "
          "once and then replaces in place")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--out", type=Path, help="a copy of the tex tree to write into (DIR/appendices/taxonomy_content.tex)")
    ap.add_argument("--agreement", type=Path, help="read the difficulty-agreement figures from this file")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    agreement = load_agreement(a.agreement)
    if a.selftest:
        return selftest()
    if a.check:
        return check(agreement)
    if a.write:
        write(a.out, agreement)
        return 0
    print("% Table: templates by branch (Branch & Dom. & Areas & Easy & Int. & Adv. & Total)")
    print("\n".join(branch_rows))
    print("\n% Table: templates by the number of variants among their 15 instances (Variants & Wording & Paths & Answer forms)")
    print("\n".join(variation_rows))
    print("\n% Numbers")
    for k, v in numbers.items():
        print(f"{k}: {v}")
    print("\n% Phrases the prose must contain")
    print("\n".join(prose))
    for name, body in blocks(agreement).items():
        print(f"\n% Generated block {name}" + ("" if agreement or name != "tab:levels_agreement" else " (stand-in: no agreement.json yet)"))
        print(block(name, body))
    return 0


if __name__ == "__main__":
    sys.exit(main())
