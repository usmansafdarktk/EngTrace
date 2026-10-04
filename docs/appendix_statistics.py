"""Compute the numbers of the appendix's Dataset Statistics section and check them against the .tex.

    python docs/appendix_statistics.py           # print the tables' LaTeX rows and the numbers the text uses
    python docs/appendix_statistics.py --check   # exit 1 unless overleaf_source_04102026/appendices/taxonomy_content.tex
                                                 # holds every generated row and number

Sources: full_run_28092026/manifest.jsonl (templates, domains, areas, levels, answer kinds, instances),
full_run_28092026/diversity.json (question wording, reasoning paths and answer forms, near-duplicates) and
full_run_28092026/FREEZE.json (how the evaluation set was drawn). Nothing is typed by hand.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "full_run_28092026"
TEX = ROOT / "overleaf_source_04102026/appendices/taxonomy_content.tex"
BRANCHES = ["chemical", "electrical", "civil", "industrial", "mechanical"]  # the taxonomy figure's order
LEVELS = ["Easy", "Intermediate", "Advanced"]
KINDS = ["scalar", "multipart", "symbolic", "vector", "array", "classification"]

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

if "--check" in sys.argv:
    tex = " ".join(TEX.read_text(encoding="utf-8").split())
    missing = [s for s in branch_rows + variation_rows + prose if " ".join(s.split()) not in tex]
    for s in missing:
        print("MISSING:", s)
    print(f"{len(branch_rows) + len(variation_rows) + len(prose) - len(missing)} of "
          f"{len(branch_rows) + len(variation_rows) + len(prose)} generated rows and numbers are in {TEX.name}")
    sys.exit(1 if missing else 0)

print("% Table: templates by branch (Branch & Dom. & Areas & Easy & Int. & Adv. & Total)")
print("\n".join(branch_rows))
print("\n% Table: templates by the number of variants among their 15 instances (Variants & Wording & Paths & Answer forms)")
print("\n".join(variation_rows))
print("\n% Numbers")
for k, v in numbers.items():
    print(f"{k}: {v}")
print("\n% Phrases the prose must contain")
print("\n".join(prose))
