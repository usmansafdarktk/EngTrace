"""Generate the tables, the figure and the numbers of the paper's Results and Error Analysis (Sections 5.2 and 5.3) and
their appendices, and check them in the source.

    python full_run_28092026/paper_results.py           # print the phrases the prose must contain
    python full_run_28092026/paper_results.py --write   # rewrite every generated block in the .tex files and draw the figure
    python full_run_28092026/paper_results.py --check   # exit 1 unless every generated block is current, every phrase is in
                                                        # its file, every number in the prose is a phrase's, every citation
                                                        # key resolves and every \\autoref label is defined

Files written: overleaf_source_04102026/6_results.tex (Table 1 and the prose phrases), appendices/results.tex,
appendices/rewording.tex, appendices/conditions.tex, appendices/error_analysis.tex (their tables, between
"% BEGIN GENERATED <name>" and "% END GENERATED <name>" markers), and figs/error-categories.pdf.

Sources, each written by a committed script: results/results.json (analyze.py: every score, interval, test and
condition), expert_request/scored.json and build.json (expert_kits.py: the experts' readings of wrong answers and of
the templates, counts only), RESIDUAL_INCORRECT.md (residual_incorrect.py: the remaining incorrect verdicts by their
distance to the target), PARAPHRASE.md and PARAPHRASE_REVIEW.md (paraphrase.py, paraphrase_kit.py --score: the
rewording funnel), FLAG_REVIEW_3.md (flag_sample.py: the expert reading of the arithmetic flags), the shortcut
template list of analyze.py, and the judged step check's precision row of appendices/validation.tex
(docs/appendix_evaluation.py). Model display names are those of paper_setup.py.
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SRC = REPO / "overleaf_source_04102026"
APPX = SRC / "appendices"
FILES = {"main": SRC / "6_results.tex", "results": APPX / "results.tex", "rewording": APPX / "rewording.tex",
         "conditions": APPX / "conditions.tex", "errors": APPX / "error_analysis.tex"}
FIG = SRC / "figs" / "error-categories.pdf"
BIB = SRC / "custom.bib"
MARK = "% BEGIN GENERATED {name} (full_run_28092026/paper_results.py --write)\n{body}\n% END GENERATED {name}"

NAME = {
    "deepseek-v4.1-flash": "DeepSeek V4.1 Flash", "kimi-k3": "Kimi K3", "claude-sonnet-5": "Claude Sonnet 5",
    "glm-5.3-flash": "GLM-5.3-Flash", "muse-glimmer-30b": "Muse Glimmer 30B", "glm-5.3": "GLM-5.3",
    "qwen3-235b-a22b-2507": "Qwen3-235B-2507", "gemini-3.1-flash-lite": "Gemini 3.1 Flash-Lite",
    "gemma-4-26b-a4b": "Gemma 4 26B", "gpt-5.4-mini": "GPT-5.4 mini", "gpt-oss-20b": "gpt-oss-20b",
    "gpt-5.4": "GPT-5.4", "deepseek-v4-pro": "DeepSeek V4 Pro",
}
NO_REASONING = {"gemma-4-26b-a4b", "qwen3-235b-a22b-2507", "gemini-3.1-flash-lite", "gpt-5.4-mini"}
BRANCH = {"chemical_engineering": "Chemical", "civil_engineering": "Civil", "electrical_engineering": "Electrical",
          "industrial_engineering": "Industrial", "mechanical_engineering": "Mechanical"}
LEVELS = ["Easy", "Intermediate", "Advanced"]
CONDITION = {"reasoning-medium": "Reasoning at medium effort", "openbook2": "Governing equations supplied",
             "tool": "Python tool offered", "flagship-reasoning-medium": "Reasoning at medium effort"}
# The experts' category labels, in the hierarchy's order, and the short form the tables and the figure use.
CATEGORIES = [("1. Hallucination", "Hallucination"), ("2. Setup / Assumption Error", "Setup or assumption"),
              ("3. Formula / Principle Error", "Formula or principle"), ("4. Unit / Dimensional Error", "Unit"),
              ("5. Sign / Direction Error", "Sign or direction"), ("6. Calculation Error", "Calculation"),
              ("No error: the answer is correct, or the question admits it", "No error"),
              ("Incomplete: the working stops before an answer", "Incomplete")]
B2_MODELS = ["claude-sonnet-5", "gpt-5.4-mini", "gemma-4-26b-a4b", "gpt-oss-20b"]
WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
        11: "eleven", 12: "twelve"}
NUMBER = re.compile(r"(?<![\w.\-])\d[\d,]*(?:\.\d+)?")
NUMBER_WORDS = re.compile(r"\b(two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)\b", re.I)


# ----------------------------------------------------------------------------------------------- formatting
def f3(x: float) -> str:
    return f"{x:.3f}"


def f2(x: float) -> str:
    return f"{x:.2f}"


def ci(c, dash: bool = True) -> str:
    """An interval: 'lo--hi' in tables when both ends are positive, otherwise 'lo to hi' with a math-mode minus."""
    if dash and c[0] >= 0 and c[1] >= 0:
        return f"{f3(c[0])}--{f3(c[1])}"
    return " to ".join(f"$-{f3(abs(x))}$" if round(x, 3) < 0 else f3(x) for x in c)


def sgn(x: float, nd: int = 3) -> str:
    r = round(x, nd)
    s = f"{abs(r):.{nd}f}"
    return f"$-{s}$" if r < 0 else (f"$+{s}$" if r > 0 else f"${s}$")


def pv(p: float) -> str:
    return "$<$0.001" if p < 0.0005 else f"{p:.3f}"


def pct(x: float) -> str:
    return f"{math.floor(x * 100 + 0.5)}\\%"


def pct1(x: float) -> str:
    return f"{x * 100:.1f}\\%"


def tt(key: str) -> str:
    return "\\texttt{" + NAME[key] + "}"


def code(name: str) -> str:
    return "\\texttt{" + name.replace("template_", "").replace("_", "\\_") + "}"


def thousands(n: int) -> str:
    return f"{n:,}"


def block(name: str, body: str) -> str:
    return MARK.format(name=name, body=body.rstrip("\n"))


def table(spec: str, header: str, rows: list[str], caption: str, label: str, star: bool = True,
          size: str = "\\small", resize: bool = False) -> str:
    env = "table*" if star else "table"
    lines = [f"\\begin{{{env}}}[t]", "\\centering", size]
    if resize:
        lines.append("\\resizebox{\\textwidth}{!}{%")
    lines += [f"\\begin{{tabular}}{{{spec}}}", "\\toprule", "\\rowcolor{gray!10}", header + " \\\\", "\\midrule"]
    lines += [r if r == "\\midrule" else r + " \\\\" for r in rows]
    lines += ["\\bottomrule", "\\end{tabular}"]
    if resize:
        lines.append("}")
    lines += [f"\\caption{{{caption}}}", f"\\label{{{label}}}", f"\\end{{{env}}}"]
    return "\n".join(lines)


def head(*cells: str) -> str:
    return " & ".join("\\textbf{" + c + "}" for c in cells)


def wrap(text: str, width: int = 100) -> str:
    """One sentence per line, wrapped at the width, never inside a macro argument."""
    out = []
    for sentence in text.split("\n"):
        words, line = sentence.split(" "), ""
        for w in words:
            if line and len(line) + 1 + len(w) > width:
                out.append(line)
                line = w
            else:
                line = f"{line} {w}" if line else w
        out.append(line)
    return "\n".join(out)


# ----------------------------------------------------------------------------------------------- sources
def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def md_table(text: str, header_pattern: str) -> list[list[str]]:
    lines = text.split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("|") and re.search(header_pattern, l)) + 2
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip().strip("`") for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return rows


def md_kv(text: str) -> dict[str, str]:
    """The two-column '| what | value |' tables of the paraphrase records, as one dict."""
    out = {}
    for m in re.finditer(r"^\| ([^|]+?) \| ([^|]*?) \|$", text, re.M):
        out[m.group(1).strip()] = m.group(2).strip()
    return out


res = load(HERE / "results/results.json")
scored = load(HERE / "expert_request/scored.json")
built = load(HERE / "expert_request/build.json")
residual = (HERE / "RESIDUAL_INCORRECT.md").read_text(encoding="utf-8")
paraphrase = (HERE / "PARAPHRASE.md").read_text(encoding="utf-8")
review = (HERE / "PARAPHRASE_REVIEW.md").read_text(encoding="utf-8")
flags = (HERE / "FLAG_REVIEW_3.md").read_text(encoding="utf-8")
validation_tex = (APPX / "validation.tex").read_text(encoding="utf-8")
analyze_src = (HERE / "analyze.py").read_text(encoding="utf-8")
traces_src = (HERE / "run_traces.py").read_text(encoding="utf-8")
paraphrase_src = (HERE / "paraphrase.py").read_text(encoding="utf-8")
TOOL_MAX_CALLS, TOOL_TIMEOUT, TOOL_OUTPUT_CHARS = (int(re.search(rf"^{c} = (\d+)", traces_src, re.M).group(1))
                                                  for c in ("TOOL_MAX_CALLS", "TOOL_TIMEOUT", "TOOL_OUTPUT_CHARS"))
COPY = float(re.search(r"^COPY = ([\d.]+)", paraphrase_src, re.M).group(1))  # the paraphrase's word-similarity ceiling

Q1 = {m["model"]: m for m in res["q1"]["models"]}
ORDER = sorted(Q1, key=lambda k: -Q1[k]["score"])  # the table order: Final Answer Accuracy, descending
PAIRS = res["q1"]["pairs"]
Q2 = {m["model"]: m for m in res["q2"]}
Q3 = {m["model"]: m for m in res["q3"]}
Q3O = {m["model"]: m for m in res["q3_overall"]}
COV = {m["model"]: m for m in res["q3_coverage"]["models"]}
CPAIRS = res["q3_coverage"]["pairs"]
BL = {m["model"]: m for m in res["branches_levels"]}
Q4 = {m["model"]: m for m in res["q4"]}
Q5 = {m["model"]: m for m in res["q5"]["models"]}
SENS = {m["model"]: m for m in res["sensitivity"]["models"]}
REP = res["reported"]["models"]
ARMS = [a for a in res["reasoning_arms"] if a["arm"] != "openbook"]  # version 1 of the open-book arm is superseded
ANCH = {a["arm"]: a for a in res["anchors"]}
REPEATS = res["repeats"]
N_ITEMS = sum(m["correct"] + m["partial"] + m["incorrect"] + m["unusable"] for m in [Q1[ORDER[0]]])
N_TEMPLATES = Q2[ORDER[0]]["templates"]
N_PAIRS = len(PAIRS)
N_SINGLE = Q4[ORDER[0]]["single_path"]["templates"]
assert all(m["correct"] + m["partial"] + m["incorrect"] + m["unusable"] == N_ITEMS for m in Q1.values())
assert N_PAIRS == len(CPAIRS) == 55 and N_ITEMS == 2250 and N_TEMPLATES == 150

ERR = scored["kinds"]["error"]
ERR_BUILD = built["kinds"]["error"]
TEMPLATES_READ = scored["kinds"]["template"]["templates"]
SHORTCUT = re.findall(r"'(template_[a-z_]+)'", re.search(r"SHORTCUT = \[(.*?)\]", analyze_src, re.S).group(1))
SYMBOLIC = res["symbolic_templates"]


# ----------------------------------------------------------------------------------------------- derived numbers
def separated(pairs, key="p_holm") -> set[frozenset]:
    return {frozenset((p["a"], p["b"])) for p in pairs if p[key] < 0.05}


def letters(order: list[str], sep: set[frozenset]) -> dict[str, str]:
    """A compact letter display (insert-and-absorb): models that share a letter are not separated."""
    groups = [set(order)]
    for pair in sorted(sep, key=lambda q: sorted(order.index(m) for m in q)):
        a, b = sorted(pair, key=order.index)
        new = []
        for g in groups:
            new += [g - {a}, g - {b}] if a in g and b in g else [g]
        uniq = []
        for g in new:
            if g and g not in uniq:
                uniq.append(g)
        groups = [g for g in uniq if not any(g < h for h in uniq)]
    groups.sort(key=lambda g: sorted(order.index(m) for m in g))
    return {m: "".join(chr(ord("a") + i) for i, g in enumerate(groups) if m in g) for m in order}


SEP_FAC, SEP_MC = separated(PAIRS), separated(CPAIRS)
CLD_FAC, CLD_MC = letters(ORDER, SEP_FAC), letters(ORDER, SEP_MC)
TOP5 = ORDER[:5]
assert not any(frozenset(p) in SEP_FAC for p in __import__("itertools").combinations(TOP5, 2))
nonsig = [p for p in PAIRS if p["p_holm"] >= 0.05]
assert all(abs(p["diff"]) < p["detectable"] for p in nonsig)
strict_agree = sum((p["p_holm"] < 0.05) == (p["fully_p_holm"] < 0.05) for p in PAIRS)
scores = [Q1[k]["score"] for k in ORDER]
spread5 = max(scores[:5]) - min(scores[:5])
sixth, rest = ORDER[5], ORDER[6:]
sixth_sep = [k for k in TOP5 if frozenset((sixth, k)) in SEP_FAC]
assert all(frozenset((a, b)) in SEP_FAC for a in rest for b in ORDER[:6])

gaps = [Q2[k]["gap"] for k in ORDER]
gap_sig = [k for k in ORDER if Q2[k]["p_welch_holm"] < 0.05]
gap_sig_chem = [k for k in ORDER if Q2[k]["without_two_chemical"]["p_welch_holm"] < 0.05]
gap_perm = sum(Q2[k]["p_perm_holm"] < 0.05 for k in ORDER)
assert gap_sig == ORDER[-4:] and not gap_sig_chem
top5_gaps = [Q2[k]["gap"] for k in TOP5]
top5_excl0 = sum(Q2[k]["ci"][0] > 0 for k in TOP5)
glm = Q2["glm-5.3"]
branch_pairs = [(k, p) for k in ORDER for p in BL[k]["pairs"] if p["p_holm"] < 0.05]
assert branch_pairs == [("gpt-oss-20b", branch_pairs[0][1])] and branch_pairs[0][1]["a"] == "civil_engineering" \
    and branch_pairs[0][1]["b"] == "electrical_engineering" and branch_pairs[0][1]["diff"] < 0

cov = [COV[k]["coverage"] for k in ORDER]
tau = res["q3_coverage"]["tau"]
cov_rank = {k: COV[k]["coverage_rank"] for k in ORDER}
ans_rank = {k: COV[k]["answer_rank"] for k in ORDER}
first_cov = min(ORDER, key=lambda k: cov_rank[k])
first_fac = ORDER[0]
dc = next(p for p in CPAIRS if {p["a"], p["b"]} == {first_fac, first_cov})
dc_diff = dc["diff"] if dc["a"] == first_fac else -dc["diff"]
assert frozenset((first_fac, first_cov)) not in SEP_FAC and dc["p_holm"] < 0.05
n_mc_sig, n_wil_sig = len(SEP_MC), len(separated(CPAIRS, "p_wilcoxon_holm"))
agree = sum((p["p_holm"] < 0.05) == (p["p_wilcoxon_holm"] < 0.05) for p in CPAIRS)
top3_cov = sorted(ORDER, key=lambda k: cov_rank[k])[:3]
assert not any(frozenset(p) in SEP_MC for p in __import__("itertools").combinations(top3_cov, 2))
rho = [COV[k]["rho_steps"] for k in ORDER]

e3w = [Q3[k]["e3_coverage_on_readable_wrong"] for k in ORDER]
floor = [Q3[k]["e3_null_on_readable_wrong"] for k in ORDER]
e5w = [Q3[k]["e5_coverage_on_readable_wrong"] for k in ORDER]
e5w_max_model = max(ORDER, key=lambda k: Q3[k]["e5_coverage_on_readable_wrong"])
many_wrong = [k for k in ORDER if COV[k]["wrong_n"] > 100]
full_cov = [COV[k]["wrong_full_coverage"] for k in many_wrong]
attr = {c: [Q3[k]["attribution_on_wrong"][c] for k in ORDER] for c in ("digit_rule", "e5_missing", "router_judge")}

digit = [Q3[k]["digit_flag_rate_on_fully_solved"] for k in ORDER]
router = [Q3[k]["router_judge_rate_on_fully_solved"] for k in ORDER]
claims = [Q3[k]["claims_per_trace"] for k in ORDER]
hi_digit = max(ORDER, key=lambda k: Q3[k]["digit_flag_rate_on_fully_solved"])
lo_digit = min(ORDER, key=lambda k: Q3[k]["digit_flag_rate_on_fully_solved"])
flag_row = md_table(flags, r"\| model \| flags drawn")
flag_all = next(r for r in flag_row if r[0] == "all")
flags_read, flags_slip, flags_checker, flags_unsure = (int(flag_all[i]) for i in (2, 3, 4, 5))
flags_decided = flags_read - flags_unsure
flag_precision = flags_slip / flags_decided
assert f3(flag_precision) == flag_all[6]
router_pr = re.search(r"Judged step and arithmetic checks & steps, correct answers & P / R & ([\d.]+) / ([\d.]+)",
                      validation_tex)
router_precision, router_recall = float(router_pr.group(1)), float(router_pr.group(2))

depth6 = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in ORDER]
depth1 = [Q3[k]["by_milestone_count"]["1"]["wrong_rate"] for k in ORDER]
assert all(a > b for a, b in zip(depth6, depth1))
single_all = {k: Q4[k]["single_path"]["all"] for k in ORDER}
single_some = {k: Q4[k]["single_path"]["some"] for k in ORDER}
strong6, weak5 = ORDER[:6], ORDER[6:]
weak4 = [k for k in weak5 if k != "gemini-3.1-flash-lite"]

q5 = [Q5[k] for k in ORDER]
q5_diff = [m["diff"] for m in q5]
q5_within = [k for k in ORDER if Q5[k]["within_margin"]]
q5_out = [k for k in ORDER if not Q5[k]["within_margin"]]
assert len(q5_out) == 1 and all(m["p_holm"] >= 0.05 for m in q5) and all(m["items"] == q5[0]["items"] for m in q5)
q5_pairs, q5_templates, margin = q5[0]["items"], q5[0]["templates"], res["q5"]["margin"]
q5_tau = res["q5"]["tau"]
noise = q5_tau["noise_arm"]
para = md_kv(paraphrase)
rev = md_kv(review)
p_selected, p_passing, p_failed = int(para["items selected (subsamples.PARAPHRASE)"]), int(para["passing a paraphrase"]), \
    int(para["no paraphrase after 3 attempts"])
p_attempts = [int(x) for x in para["passing at attempt 1, 2, 3"].split(", ")]
p_restored = int(para["passing ones with notation restored (Unicode the original lacks, put back as written)"].split()[0])
p_lost_templates = para["templates with no paraphrase for any of their items"].split(": ")[1].split(", ")
p_lost = int(para["templates with no paraphrase for any of their items"].split(":")[0])
p_checks = dict(re.findall(r"(\w+) (\d+)", para["failed attempts by check"]))
p_branch_pass = dict((b, int(n)) for b, n in re.findall(r"(\w+) (\d+)", para["passing, per branch"]))
r_kept, r_rejected, r_returned = int(rev["kept"]), int(rev["rejected"]), int(rev["returned"])
r_branch = {b: (int(k), int(r)) for b, k, r in re.findall(r"(\w+) (\d+) / (\d+)", rev["kept / rejected, per branch"])}
r_reasons = dict((k, int(v)) for k, v in re.findall(r"(\w+=\w+) (\d+)", rev["answers behind the rejections"]))
assert p_selected == 450 and p_passing == r_returned == 316 and r_kept == q5_pairs == 277 and r_rejected == 39
assert p_lost == len(p_lost_templates) and p_selected - p_passing == p_failed
p_lost_experts = N_TEMPLATES - q5_templates - p_lost
below90 = [k for k in ORDER if Q5[k]["ci90"][1] < 0]
above90 = [k for k in ORDER if Q5[k]["ci90"][0] > 0]
vs = res["q5"]["vs_repeats"]
vs_within = [k for k in ORDER if k in vs and vs[k]["abs_within_repeat_spread"]]
vs_beyond = [k for k in ORDER if k in vs and not vs[k]["abs_within_repeat_spread"]]
rep_sd = [REPEATS[k]["sd"] for k in REPEATS]
rep_same = [REPEATS[k]["same_verdict_every_repeat"] for k in REPEATS]
sens_tau = res["sensitivity"]["tau_with_headline"]
short_shift = max(abs(SENS[k]["without_shortcut_templates"] - SENS[k]["fitted"]) for k in ORDER)
sym_shift = [SENS[k]["without_symbolic_templates"] - SENS[k]["fitted"] for k in ORDER]
assert sens_tau["half_tol"] == sens_tau["double_tol"]

arm = {(a["arm"], a["model"]): a for a in ARMS}
r_mini, r_gem = arm[("reasoning-medium", "gpt-5.4-mini")], arm[("reasoning-medium", "gemini-3.1-flash-lite")]
assert r_gem["p_holm"] >= 0.05 and abs(r_gem["diff"]) < r_gem["detectable"] and r_mini["p_holm"] < 0.05
ob = {m: arm[("openbook2", m)] for m in ("gpt-oss-20b", "gpt-5.4-mini", "claude-sonnet-5")}
assert ob["gpt-5.4-mini"]["within_margin"] and ob["claude-sonnet-5"]["within_margin"] and not ob["gpt-oss-20b"]["within_margin"]
tool = {m: arm[("tool", m)] for m in ("claude-sonnet-5", "gpt-5.4-mini")}
assert all(t["within_margin"] and t["e5"]["within_margin"] for t in tool.values())
fl = {a["model"]: a for a in ANCH["flagship"]["anchors"]}
flr = ANCH["flagship-reasoning-medium"]["anchors"][0]
roster_sub = sorted(ANCH["flagship"]["roster"], key=lambda x: -x["score"])
top5_sub = roster_sub[:5]
assert {x["model"] for x in top5_sub} == set(TOP5)
sub_scores = [x["score"] for x in top5_sub]
sub_adv = [x["levels"]["Advanced"]["mean"] for x in top5_sub]
anchors_adv = [fl["deepseek-v4-pro"]["levels"]["Advanced"]["mean"], fl["gpt-5.4"]["levels"]["Advanced"]["mean"],
               flr["levels"]["Advanced"]["mean"]]
assert all(x["ci"][0] <= fl["deepseek-v4-pro"]["score"] <= x["ci"][1] and x["ci"][0] <= flr["score"] <= x["ci"][1]
           for x in top5_sub)
ob_items, ob_templates = ob["gpt-oss-20b"]["items"], ob["gpt-oss-20b"]["templates"]
mc_rise = [ob[m]["e5"]["diff"] for m in ob]

# The experts' reading of wrong answers (B2).
by_model = {m: ERR["by_model"][m] for m in B2_MODELS}
by_level = ERR["by_level"]
maj = ERR["majority_by_model"]
per_model_items = ERR_BUILD["items"] // len(B2_MODELS)
readings_total = ERR["readings"]
assert readings_total == 480 and per_model_items == 40 and sum(sum(v.values()) for v in by_level.values()) == 480
level_n = {lv: sum(by_level[lv].values()) for lv in LEVELS}
CALC, FORM, NOERR = CATEGORIES[5][0], CATEGORIES[2][0], CATEGORIES[6][0]
easy_calc, easy_n = by_level["Easy"].get(CALC, 0), level_n["Easy"]
form_hard = sum(by_level[lv].get(FORM, 0) for lv in ("Intermediate", "Advanced"))
hard_n = level_n["Intermediate"] + level_n["Advanced"]
form_hard_share, form_easy_share = form_hard / hard_n, by_level["Easy"].get(FORM, 0) / easy_n
claude_noerr = maj["claude-sonnet-5"][NOERR]
others = ["gpt-5.4-mini", "gemma-4-26b-a4b", "gpt-oss-20b"]
others_calc = [maj[m].get(CALC, 0) for m in others]
others_form = [maj[m].get(FORM, 0) for m in others]
assert all(maj[m].get(CALC, 0) == max(maj[m].values()) for m in others) and claude_noerr == max(maj["claude-sonnet-5"].values())
assert all(sorted(maj[m].values())[-2] == maj[m].get(FORM, 0) for m in others)
fleiss_all = ERR["fleiss_all"]

# The remaining incorrect verdicts of the top five (RESIDUAL_INCORRECT.md).
res_model = {r[0]: r for r in md_table(residual, r"^\| \| n \| <=0\.2%")}
top5_rows = [res_model[k] for k in TOP5]
top5_incorrect = sum(int(r[1]) for r in top5_rows)
top5_near = sum(int(r[2]) for r in top5_rows)
top5_symbolic = sum(int(r[9]) for r in top5_rows)
tpl_rows = residual.split("### Top five models, by template")[1].split("###")[0]
tpl = {r[0]: r for r in md_table(tpl_rows, r"^\| \| n \|") if r[0] != "all"}
one_template_symbolic = max(int(r[9]) for r in tpl.values())
two_chemical = int(tpl["template_work_isothermal_virial"][1]) + int(tpl["template_adiabatic_flame_temperature"][1])
near_rows = md_table(residual.split("### Within 0.2% and still incorrect")[1].split("###")[0], r"^\| template \| n \|")
assert all(r[1] == r[2] for r in near_rows), "every incorrect verdict within 0.2% is an exact-digits case"
exact_m = re.search(r"Templates with an exact-digits target .*?: (\d+) of 150; incorrect verdicts on them across the roster: (\d+)",
                    residual)
exact_templates, exact_incorrect = int(exact_m.group(1)), int(exact_m.group(2))
form_total = top5_symbolic + top5_near + two_chemical
adv_top5 = [BL[k]["level"]["Advanced"]["mean"] for k in TOP5]
near_templates = [t for t, v in TEMPLATES_READ.items() if v["role"] == "near"]
near_total = sum(int(r[1]) for r in near_rows)
fleiss_models = list(ERR["fleiss_by_model"].values())


# ----------------------------------------------------------------------------------------------- phrases
def rng(xs, fmt=f3, join=" to ") -> str:
    return f"{fmt(min(xs))}{join}{fmt(max(xs))}"


phrases = [  # 6_results.tex must contain each of these, whitespace aside
    "with 95\\% intervals and a compact letter display", f"over the {N_PAIRS} pairs", f"FAC runs from {f3(min(scores))}",
    f"to {f3(max(scores))}",
    f"{WORD[5].capitalize()} models lie within {f3(spread5)} of one another, between {f3(min(scores[:5]))} and {f3(max(scores[:5]))}",
    f"{tt(sixth)} ({f3(Q1[sixth]['score'])}) follows and differs from only {WORD[len(sixth_sep)]} of the five ({tt(sixth_sep[0])}); "
    f"{Q1[sixth]['empty']} of its responses are empty at the output ceiling and score 0",
    f"The remaining {WORD[len(rest)]} models lie between {f3(min(scores[6:]))} and {f3(max(scores[6:]))}, and each differs from every model of the top six",
    f"In all, {len(SEP_FAC)} of the {N_PAIRS} pairs differ after correction, and the {len(nonsig)} that do not all lie below the smallest "
    f"difference the design detects at 80\\% power ({rng([p['detectable'] for p in nonsig])}, depending on the pair)",
    f"Four of the five lowest scores belong to the four models that return no reasoning tokens",
    f"on the {r_mini['items']}-instance subset, {tt('gpt-5.4-mini')} moves from {f3(r_mini['main_score_on_items'])} to {f3(r_mini['arm_score'])} "
    f"(95\\% CI of the change {ci(r_mini['ci'], False)}) and {tt('gemini-3.1-flash-lite')} from {f3(r_gem['main_score_on_items'])} to "
    f"{f3(r_gem['arm_score'])} ({sgn(r_gem['ci'][0])} to {f3(r_gem['ci'][1])}), less than the design detects",
    f"{rng([Q1[k]['templates_no_instance_variance'] for k in ORDER], str)} of the {N_TEMPLATES} templates per model have the same score on all 15 instances",
    # difficulty
    f"Every model scores lower on Advanced than on Easy templates, by {rng(gaps)}, but the gap is significant after correction for "
    f"{WORD[len(gap_sig)]} of the {WORD[len(ORDER)]} models",
    f"{', '.join(tt(k) for k in gap_sig[:-1])}, and {tt(gap_sig[-1])} (gaps {rng([Q2[k]['gap'] for k in gap_sig])}), three of which run without reasoning tokens",
    f"the top five models' gaps are {rng(top5_gaps)}, with intervals that {'all exclude zero' if top5_excl0 == 5 else f'exclude zero for {WORD[top5_excl0]} of the five'} but do not survive correction, "
    f"against detectable gaps of {rng([Q2[k]['detectable_welch'] for k in ORDER])}",
    f"{tt('glm-5.3')}'s gap of {f3(glm['gap'])} is mostly an output-ceiling effect",
    f"without them the gap is {f3(glm['unusable_excluded']['gap'])}",
    f"Without these two templates the gaps are {rng([Q2[k]['without_two_chemical']['gap'] for k in ORDER])} and none holds after correction",
    f"one pair in the whole evaluation separates after correction ({tt('gpt-oss-20b')}, electrical above civil)",
    # coverage
    f"MC runs from {f3(min(cov))} to {f3(max(cov))} and orders the models differently from FAC (Kendall's $\\tau$ {f3(tau['tau'])}, 95\\% CI {ci(tau['ci'], False)}",
    f"{tt(first_cov)} is first on MC and {['first', 'second', 'third'][ans_rank[first_cov] - 1]} on FAC; {tt(first_fac)} is first on FAC and "
    f"{['first', 'second', 'third', 'fourth', 'fifth', 'sixth'][cov_rank[first_fac] - 1]} on MC, {f3(-dc_diff)} below {tt(first_cov)} "
    f"(Holm-adjusted $p = {f3(dc['p_holm'])}$)",
    f"Over the {N_PAIRS} pairs, {n_mc_sig} differ after correction under the sign-flip test and {n_wil_sig} under a Wilcoxon signed-rank test",
    f"the two agreeing on {agree} pairs; the three models highest on MC, {tt(top3_cov[0])}, {tt(top3_cov[1])}, and {tt(top3_cov[2])}, do not separate",
    f"mean number of steps is {sgn(max(rho), 2)} to {sgn(min(rho), 2)} for every model",
    f"matching alone finds {pct(min(e3w))} to {pct(max(e3w))} of the milestones, against a chance floor of {pct(min(floor))} to {pct(max(floor))}",
    f"with the judge the share is {pct(min(e5w))} to {pct(max(e5w))}, the {pct(max(e5w))} resting on the {Q3[e5w_max_model]['readable_wrong_with_milestones']} "
    f"wrong answers of {tt(e5w_max_model)}",
    f"Among the {WORD[len(many_wrong)]} models with more than 100 wrong answers, {pct(min(full_cov))} to {pct(max(full_cov))} of them reach every milestone",
    f"a milestone the judge rules missing points at the failure in {pct(min(attr['e5_missing']))} to {pct(max(attr['e5_missing']))} of responses, "
    f"a judged step flag in {pct(min(attr['router_judge']))} to {pct(max(attr['router_judge']))}, and an arithmetic flag in "
    f"{pct(min(attr['digit_rule']))} to {pct(max(attr['digit_rule']))}",
    # diagnostics
    f"The arithmetic check flags {pct1(min(digit))} to {pct1(max(digit))} of correct-answer responses, at a precision of {f3(flag_precision)} on these models "
    f"({flags_slip} of {flags_decided} flags that a domain expert read were slips), and the judged step check flags {pct1(min(router))} to {pct1(max(router))}, "
    f"at a precision of {f2(router_precision)} and a recall of {f2(router_recall)} inside correct-answer responses",
    f"reads {rng(claims, lambda x: f'{x:.1f}')} displayed calculations per response",
    f"{tt(hi_digit)}'s {pct1(Q3[hi_digit]['digit_flag_rate_on_fully_solved'])} comes with {Q3[hi_digit]['claims_per_trace']:.1f} calculations per response and "
    f"{tt(lo_digit)}'s {pct1(Q3[lo_digit]['digit_flag_rate_on_fully_solved'])} with {Q3[lo_digit]['claims_per_trace']:.1f}",
    # consistency and depth
    f"wrong-answer rates of {rng(depth6)} against {rng(depth1)}",
    f"the six strongest models solve all 15 instances of {pct(min(single_all[k] for k in strong6))} to {pct(max(single_all[k] for k in strong6))} of the "
    f"{N_SINGLE} single-path templates, whereas four of the five weakest solve all instances of {pct(min(single_all[k] for k in weak4))} to "
    f"{pct(max(single_all[k] for k in weak4))} of them ({tt('gemini-3.1-flash-lite')} {pct(single_all['gemini-3.1-flash-lite'])}) and some but not all "
    f"instances of {pct(min(single_some[k] for k in weak4))} to {pct(max(single_some[k] for k in weak4))}",
    # rewording
    f"for ten of the eleven models", f"one instance in five of every template ({p_selected} instances)",
    f"{p_passing} rewrites passed scripted checks, and a domain expert of the instance's branch confirmed {r_kept} of them, on {q5_templates} templates",
    f"On these {q5_pairs} pairs, the paired change in FAC lies between {sgn(min(q5_diff))} and {sgn(max(q5_diff))} across the eleven models",
    f"for ten models the 90\\% interval lies within $\\pm {margin:.2f}$", f"for {tt(q5_out[0])} it reaches {sgn(Q5[q5_out[0]]['ci90'][0])}",
    f"$\\tau = {f3(q5_tau['tau'])}$ ({ci(q5_tau['ci'], False)}), at the lower edge of what sampling alone produces at this size (median {f3(noise['median'])})",
    f"the experts rejected {r_rejected} of the {r_returned} pairs that passed every scripted check ({pct(r_rejected / r_returned)})",
    f"covers the {q5_templates} templates whose problems could be reworded without loss",
    f"on {REPEATS['gpt-oss-20b']['items']} instances decoded three times, four models' FAC has a standard deviation of {rng(rep_sd)} across repeats, "
    f"and {pct(min(rep_same))} to {pct(max(rep_same))} of instances receive the same verdict every time",
    f"unchanged at half and double the answer tolerance ($\\tau = {f3(sens_tau['half_tol'])}$), and FAC moves by at most {f3(short_shift)} without the "
    f"{WORD[len(SHORTCUT)]} templates whose answers can be read off the question's wording and by {sgn(min(sym_shift))} to {sgn(max(sym_shift))} without the "
    f"{WORD[len(SYMBOLIC)]} symbolic-answer templates",
    # conditions
    f"On the {ANCH['flagship']['items']}-instance subset, two flagships",
    f"{tt('deepseek-v4-pro')} scores {f3(fl['deepseek-v4-pro']['score'])} ({ci(fl['deepseek-v4-pro']['ci'], False)}), and {tt('gpt-5.4')} scores "
    f"{f3(flr['score'])} ({ci(flr['ci'], False)}) with reasoning at medium effort and {f3(fl['gpt-5.4']['score'])} ({ci(fl['gpt-5.4']['ci'], False)}) at its "
    f"provider's default, which returns no reasoning tokens; the top five on the same instances score {rng(sub_scores)}, with intervals that contain both "
    f"reasoning anchors, and on Advanced templates {rng(sub_adv, f2)} against the anchors' {rng(anchors_adv, f2)}",
    f"({ob_items} instances of the {ob_templates} templates that state them in symbols) lifts {tt('gpt-oss-20b')} by {f3(ob['gpt-oss-20b']['diff'])} "
    f"({ci(ob['gpt-oss-20b']['ci'], False)}), partly because fewer responses run out of room ({f3(ob['gpt-oss-20b']['usable_in_both']['diff'])} on the instances "
    f"answered in both conditions), and changes {tt('gpt-5.4-mini')} ({sgn(ob['gpt-5.4-mini']['diff'])}) and {tt('claude-sonnet-5')} "
    f"({sgn(ob['claude-sonnet-5']['diff'])}) by amounts bounded within $\\pm {margin:.2f}$, while MC rises by {rng(mc_rise, f2)} for all three",
    f"{tt('claude-sonnet-5')} calls it on {pct(tool['claude-sonnet-5']['tool_use']['share_with_calls'])} of instances and {tt('gpt-5.4-mini')} on "
    f"{pct(tool['gpt-5.4-mini']['tool_use']['share_with_calls'])}, and neither model's FAC ({sgn(tool['claude-sonnet-5']['diff'])}; "
    f"{sgn(tool['gpt-5.4-mini']['diff'])}) or MC changes beyond $\\pm {margin:.2f}$; {tt('gpt-5.4-mini')}'s arithmetic flags on correct answers fall from "
    f"{f3(tool['gpt-5.4-mini']['digit_flag_rate_fully_solved']['main'])} to {f3(tool['gpt-5.4-mini']['digit_flag_rate_fully_solved']['arm'])}",
    # error analysis
    f"Three domain experts read each of {ERR_BUILD['items']} wrong answers, {per_model_items} from each of four models",
    f"(Fleiss'~$\\kappa$ {f3(fleiss_all)}", f"the {readings_total} readings",
    f"For {tt('claude-sonnet-5')}, {claude_noerr} of {per_model_items} wrong answers are ``no error'' by majority",
    f"calculation errors lead ({others_calc[0]}, {others_calc[1]}, and {others_calc[2]} of {per_model_items}), followed by formula or principle errors "
    f"({others_form[0]}, {others_form[1]}, and {others_form[2]})",
    f"almost all calculation slips ({easy_calc} of {easy_n} readings), whereas wrong formulas or principles are {pct(form_hard_share)} of the readings on "
    f"Intermediate and Advanced problems against {pct(form_easy_share)} on Easy ones",
    f"(FAC {rng(scores[:5])}; {rng(adv_top5, f2)} on Advanced templates), and {form_total} of their {top5_incorrect} remaining incorrect verdicts",
    f"{WORD[3].capitalize()} domain experts' {readings_total} readings of {ERR_BUILD['items']} wrong answers",  # figure caption
    f"by model ({per_model_items} wrong answers each) and by level ({level_n['Easy']}, {level_n['Intermediate']}, and {level_n['Advanced']} readings)",
]

appendix_phrases = {
    "results": [
        f"the {N_TEMPLATES} templates", "80\\% power", "scores a partial answer 0", f"for {strict_agree} of the {N_PAIRS} pairs",
        f"{len(SEP_FAC)} of the {N_PAIRS} pairs differ on FAC", f"{n_mc_sig} and {n_wil_sig} of the {N_PAIRS} pairs differ on MC",
        f"With {BL[ORDER[0]]['branch']['chemical_engineering']['templates']} templates per branch",
        f"detect differences of {rng([BL[k]['detectable_branch'] for k in ORDER])} between two branches",
        f"one of the {len(ORDER) * 10} branch pairs", f"for {WORD[len(gap_sig)]} of the {WORD[len(ORDER)]} models",
        f"finds {gap_perm} of {len(ORDER)}", f"without the {WORD[len(SHORTCUT)]} templates",
        f"{tt('glm-5.3')}'s {Q1['glm-5.3']['empty']} empty responses",
        f"for the {N_SINGLE} templates that follow one reasoning path and for the other {N_TEMPLATES - N_SINGLE}", "on all 15 instances",
        f"of {WORD[len(REPEATS)]} models on {REPEATS['gpt-oss-20b']['items']} instances decoded three times; "
        f"the other {WORD[len(ORDER) - len(REPEATS)]} models",
        f"the {res['milestones']['items_without']} instances with no milestones",
    ],
    "rewording": [
        f"three of its 15 instances ({p_selected} instances)", f"at most {COPY}", "three attempts",
        f"the experts rejected {pct(r_rejected / r_returned)} of the pairs",
        f"from the {p_selected} instances to the {r_kept} kept pairs on {q5_templates} templates: the {p_lost} templates with no passing rewrite",
        f"and {WORD[p_lost_experts]} more lost", "with 95\\% intervals", "the 90\\% interval against", f"$\\pm {margin:.2f}$",
        f"{WORD[len(q5_within)]} of the {WORD[len(ORDER)]} 90\\% intervals lie within the margin, and {tt(q5_out[0])}'s reaches "
        f"{sgn(Q5[q5_out[0]]['ci90'][0])}",
        f"{WORD[len(below90)].capitalize()} 90\\% intervals lie wholly below zero ({' and '.join(tt(k) for k in below90)}) and "
        f"{'one' if len(above90) == 1 else WORD[len(above90)]} wholly above ({' and '.join(tt(k) for k in above90)})",
        f"is {f3(q5_tau['tau'])} (95\\% CI {ci(q5_tau['ci'], False)})",
        f"median $\\tau$ of {f3(noise['median'])} (quartiles {f3(noise['q1'])} to {f3(noise['q3'])}, 5th percentile {f3(noise['p5'])})",
        "two random halves",
        f"for the {WORD[len(vs)]} models with decoding repeats; the rewording change lies within the repeats' spread for "
        f"{WORD[len(vs_within)]} of them and beyond it for {' and '.join(tt(k) for k in vs_beyond)}",
    ],
    "conditions": [
        f"{WORD[4]} conditions", f"{ANCH['flagship']['items']}-instance subset (three instances per template)", "80\\% power",
        "the 90\\% interval against", f"$\\pm {margin:.2f}$", f"For the {ob_templates} templates",
        f"the {N_TEMPLATES - ob_templates} templates without such a statement", f"holds {ob_items} instances",
        f"a {TOOL_TIMEOUT} s limit", f"truncated at {thousands(TOOL_OUTPUT_CHARS)} characters", f"up to {TOOL_MAX_CALLS} calls",
        f"by {f3(r_mini['diff'])} (95\\% CI {ci(r_mini['ci'], False)})", "two closed models", f"the {WORD[len(ORDER)]} evaluated models",
    ],
    "errors": [
        f"The taxonomy has {WORD[6]} categories", f"answers the {WORD[6]} questions", f"{WORD[2].capitalize()} options cover what the six do not",
        f"{WORD[3].capitalize()} domain experts of its branch",
        f"We read {ERR_BUILD['items']} wrong answers, {per_model_items} from each of {WORD[len(B2_MODELS)]} models", f"the {per_model_items} are drawn",
        f"{ERR_BUILD['templates']} templates in all",
        f"Fleiss'~$\\kappa$ of {f3(fleiss_all)} over the {readings_total} readings ({rng(fleiss_models)} per model)", f"{per_model_items} per model",
        f"are {pct(form_hard_share)} of the readings on Intermediate and Advanced problems against {pct(form_easy_share)} on Easy ones, where "
        f"calculation slips are {easy_calc} of {easy_n} readings",
        f"{WORD[3].capitalize()} chemical domain experts",
        f"three experts of the template's branch read each of {WORD[len(near_templates)]} templates on which many wrong answers land between "
        f"0.2\\% and 5\\%", "within 5\\%",
        f"{exact_templates} templates prescribe the digits",
        f"accounts for {exact_incorrect} incorrect verdicts across the {WORD[len(ORDER)]} models", f"{near_total} of them within 0.2\\%",
        f"the {WORD[len(SYMBOLIC)]} templates with symbolic answers",
        f"Of the top five models' {top5_incorrect} incorrect verdicts, {top5_symbolic} are symbolic, {top5_near} within 0.2\\% on a template that "
        f"prescribes the digits, and {two_chemical} on the two chemical templates: {form_total} of the {top5_incorrect}",
        f"{WORD[len(SHORTCUT)]} templates whose answers can be read off", f"at most {f3(short_shift)} without them",
    ],
}

# ----------------------------------------------------------------------------------------------- tables
blocks: dict[str, dict[str, str]] = {k: {} for k in FILES}


def mark(k: str) -> str:
    s = tt(k)
    if k in NO_REASONING:
        s += "$^{\\ast}$"
    if k == "glm-5.3":
        s += "$^{\\dagger}$"
    return s


rows = [f"{mark(k)} & {f3(Q1[k]['score'])} & {ci(Q1[k]['ci'])} & {CLD_FAC[k]} & {f3(COV[k]['coverage'])} & {ci(COV[k]['ci'])} & {CLD_MC[k]}"
        for k in ORDER]
blocks["main"]["tab:main_results"] = table(
    "l r c c r c c", head("Model", "FAC", "95\\% CI", "Tier", "MC", "95\\% CI", "Tier"), rows,
    "\\textbf{Final Answer Accuracy and Milestone Coverage.} "
    f"FAC over the {thousands(N_ITEMS)} instances and MC over the {res['q3_coverage']['templates']} templates with milestones, each with a 95\\% "
    "interval that resamples templates, in order of FAC. "
    f"Models that share a letter in a tier column do not differ on that measure after Holm correction over the {N_PAIRS} pairs; the letters are a "
    "compact letter display~\\citep{piepho2004letter}. "
    "$^{\\ast}$Returns no reasoning tokens at its provider's defaults. "
    f"$^{{\\dagger}}${Q1['glm-5.3']['empty']} responses are empty at the output ceiling and score 0.",
    "tab:main_results")

# Appendix: extended table, two halves.
rows = [f"{tt(k)} & {f3(Q1[k]['score'])} & {ci(Q1[k]['ci'])} & {f3(Q1[k]['fully_solved'])} & {ci(Q1[k]['fully_ci'])} & {Q1[k]['unusable']} & "
        f"{f3(Q1[k]['within_template_sd'])} & {f3(Q1[k]['between_template_sd'])}" for k in ORDER]
blocks["results"]["tab:results_answer"] = table(
    "l r c r c r r r", head("Model", "FAC", "95\\% CI", "Strict FAC", "95\\% CI", "No answer", "SD within", "SD between"), rows,
    "\\textbf{Final Answer Accuracy in Full.} "
    "FAC and strict FAC, which scores a partial answer 0, with 95\\% intervals; the number of responses with no readable final answer, which "
    "score 0; and the standard deviation of the score within a template (across its 15 instances) and between templates (of the template means).",
    "tab:results_answer")
rows = [f"{tt(k)} & {f3(Q3O[k]['e3_all'])} & {ci(Q3O[k]['e3_all_ci'])} & {f3(COV[k]['coverage'])} & {ci(COV[k]['ci'])} & "
        f"{f3(Q3[k]['e5_judged_fraction'])} & {f3(Q3[k]['digit_flag_rate_on_fully_solved'])} & {ci(Q3[k]['digit_ci'])} & "
        f"{Q3[k]['claims_per_trace']:.1f} & {f3(Q3[k]['router_judge_rate_on_fully_solved'])} & {ci(Q3[k]['router_judge_ci'])}" for k in ORDER]
blocks["results"]["tab:results_process"] = table(
    "l r c r c r r c r r c",
    head("Model", "MC, no judge", "95\\% CI", "MC", "95\\% CI", "Judged", "Arith.\\ flags", "95\\% CI", "Calc./resp.", "Step flags", "95\\% CI"), rows,
    "\\textbf{Milestone Coverage and Diagnostics in Full.} "
    "MC by matching alone (no judge) over every response with milestones and MC with the judge (as in the main table); the share of milestones the judge "
    "decides; the share of correct-answer responses with an arithmetic flag, the displayed calculations the check reads per response, and the share "
    "with a judged step flag, with 95\\% intervals.", "tab:results_process", resize=True)

# Pairwise comparisons.
rows = [f"{tt(p['a'])} & {tt(p['b'])} & {sgn(p['diff'])} & {ci(p['ci'])} & {pv(p['p_holm'])} & {pv(p['fully_p_holm'])} & {pv(p['mcnemar_p_holm'])} & "
        f"{f3(p['detectable'])}" for p in sorted(PAIRS, key=lambda p: p["p_holm"])]
blocks["results"]["tab:pairs_fac"] = table(
    "l l r c r r r r", head("Model A", "Model B", "A $-$ B", "95\\% CI", "$p$ (Holm)", "Strict", "McNemar", "Detectable"), rows,
    "\\textbf{Pairwise Comparisons of Final Answer Accuracy.} "
    f"For each of the {N_PAIRS} pairs, the difference of FAC with its 95\\% interval, the Holm-adjusted $p$ of the sign-flip test on the per-template "
    "differences, the same test on strict FAC, McNemar's exact test on the paired verdicts~\\citep{mcnemar1947}, and the smallest difference the pair "
    f"detects at 80\\% power; {len(SEP_FAC)} pairs differ at $p < 0.05$ after correction.", "tab:pairs_fac", size="\\footnotesize")
rows = [f"{tt(p['a'])} & {tt(p['b'])} & {sgn(p['diff'])} & {ci(p['ci'])} & {pv(p['p_holm'])} & {pv(p['p_wilcoxon_holm'])} & {f3(p['detectable'])}"
        for p in sorted(CPAIRS, key=lambda p: p["p_holm"])]
blocks["results"]["tab:pairs_mc"] = table(
    "l l r c r r r", head("Model A", "Model B", "A $-$ B", "95\\% CI", "$p$ (Holm)", "Wilcoxon", "Detectable"), rows,
    "\\textbf{Pairwise Comparisons of Milestone Coverage.} "
    f"For each of the {N_PAIRS} pairs, the difference of MC over the {res['q3_coverage']['templates']} templates with milestones, its 95\\% interval, "
    "the Holm-adjusted $p$ of the sign-flip test and of a Wilcoxon signed-rank test~\\citep{wilcoxon1945} on the per-template differences, and the "
    f"smallest difference the pair detects; {n_mc_sig} and {n_wil_sig} pairs differ after correction under the two tests.", "tab:pairs_mc",
    size="\\footnotesize")

# Branch, level, domain, answer kind.
rows = []
for k in ORDER:
    b = BL[k]["branch"]
    held = [p for p in BL[k]["pairs"] if p["p_holm"] < 0.05]
    which = ", ".join(f"{BRANCH[p['b'] if p['diff'] < 0 else p['a']]} $>$ {BRANCH[p['a'] if p['diff'] < 0 else p['b']]}" for p in held) or "none"
    rows.append(f"{tt(k)} & " + " & ".join(f"{f3(b[br]['mean'])} ({ci(b[br]['ci'])})" for br in BRANCH) +
                f" & {f3(BL[k]['detectable_branch'])} & {which}")
blocks["results"]["tab:by_branch"] = table(
    "l " + "c " * 5 + "r l", head("Model", *BRANCH.values(), "Detectable", "Pairs that hold"), rows,
    "\\textbf{Final Answer Accuracy by Branch.} "
    "The mean of the branch's 30 template means with its 95\\% interval; the smallest difference between two branches that a model's spread lets 30 "
    "templates detect at 80\\% power; and the branch pairs that differ under Welch's $t$-test after Holm correction over the ten pairs within the model.",
    "tab:by_branch", resize=True)
rows = [f"{tt(k)} & " + " & ".join(f"{f3(BL[k]['level'][lv]['mean'])} ({ci(BL[k]['level'][lv]['ci'])})" for lv in LEVELS) for k in ORDER]
blocks["results"]["tab:by_level"] = table(
    "l c c c", head("Model", *(f"{lv} ({BL[ORDER[0]]['level'][lv]['templates']})" for lv in LEVELS)), rows,
    "\\textbf{Final Answer Accuracy by Level.} The mean of the level's template means with its 95\\% interval; the number of templates in parentheses.",
    "tab:by_level")
rot = " & ".join("\\rotatebox{90}{" + tt(k) + "}" for k in ORDER)
domains = sorted(REP[ORDER[0]]["domain"])
kinds = sorted(REP[ORDER[0]]["answer_type"])
rows = [f"{d.replace('_', ' ').capitalize()} & " + " & ".join(f3(REP[k]["domain"][d]) for k in ORDER) for d in domains]
rows += ["\\midrule"] + [f"{a.capitalize()} & " + " & ".join(f3(REP[k]["answer_type"][a]) for k in ORDER) for a in kinds]
blocks["results"]["tab:by_domain_kind"] = table(
    "l " + "r " * len(ORDER), "\\textbf{Domain or answer kind} & " + rot, rows,
    "\\textbf{Final Answer Accuracy by Domain and by Answer Kind.} "
    "Instance means per domain (upper rows) and per answer kind (lower rows); these carry no interval and no test and describe where each model's "
    "errors fall.", "tab:by_domain_kind", resize=True)
rows = []
for k in ORDER:
    q = Q2[k]
    rows.append(f"{tt(k)} & {f3(q['easy'])} & {f3(q['advanced'])} & {sgn(q['gap'])} & {ci(q['ci'])} & {pv(q['p_welch_holm'])} & {pv(q['p_perm_holm'])} & "
                f"{f3(q['detectable_welch'])} & {sgn(q['unusable_excluded']['gap'])} ({pv(q['unusable_excluded']['p_welch_holm'])}) & "
                f"{sgn(q['without_symbolic']['gap'])} ({pv(q['without_symbolic']['p_welch_holm'])}) & "
                f"{sgn(q['without_two_chemical']['gap'])} ({pv(q['without_two_chemical']['p_welch_holm'])})")
blocks["results"]["tab:level_gap"] = table(
    "l r r r c r r r r r r",
    head("Model", "Easy", "Advanced", "Gap", "95\\% CI", "Welch", "Perm.", "Detectable", "No unreadable", "No symbolic", "No two chemical"), rows,
    "\\textbf{The Level Gap Under Its Variations.} "
    "Easy minus Advanced, the mean of the 58 Easy template means minus the mean of the 34 Advanced ones, with its 95\\% interval; the Holm-adjusted $p$ "
    "under Welch's $t$-test and under the planned permutation of level labels; the smallest gap the design detects; and the gap with its Holm-adjusted "
    "Welch $p$ when responses with no readable answer are left out, without the nine symbolic-answer templates, and without the two Advanced chemical "
    "templates whose wording does not pin the answer.", "tab:level_gap", resize=True)

# Sensitivity.
cols = [("half_tol", "Tolerance halved"), ("fitted", "As scored"), ("double_tol", "Doubled"), ("fully_solved", "Strict"),
        ("unusable_excluded", "No unreadable"), ("without_shortcut_templates", f"No {len(SHORTCUT)} shortcut"),
        ("without_symbolic_templates", f"No {len(SYMBOLIC)} symbolic"), ("half_unit", "Half-unit window"), ("whole_trace", "Whole response")]
rows = [f"{tt(k)} & " + " & ".join(f3(SENS[k][c]) for c, _ in cols) for k in ORDER]
rows += ["\\midrule", "$\\tau$ with the scores as scored & " + " & ".join("1" if c == "fitted" else f3(sens_tau[c]) for c, _ in cols)]
blocks["results"]["tab:sensitivity"] = table(
    "l " + "r " * len(cols), head("Model", *(n for _, n in cols)), rows,
    "\\textbf{Sensitivity of Final Answer Accuracy to the Scoring Rules.} "
    "FAC with the answer tolerance halved and doubled; strict FAC; the mean over the responses with a readable answer; without the four templates "
    "whose answers can be read off the question's wording and without the nine templates with symbolic answers, which we score by the numbers they "
    "state; with a half-unit rounding window in place of one unit of the last digit; and crediting a requested quantity stated in the body but left "
    f"off the final answer line. The last row is Kendall's $\\tau$ between that ordering of the models and the ordering as scored; sampling noise alone "
    f"puts $\\tau$ at a median of {f3(res['sensitivity']['tau_noise']['median'])} between random halves of the instances.", "tab:sensitivity", resize=True)

# Consistency and instance variance.
rows = [f"{tt(k)} & " + " & ".join(f"{f3(Q4[k][g][s])} ({ci(Q4[k][g][s + '_ci'])})" for g in ("single_path", "multi_path") for s in ("all", "some", "none"))
        for k in ORDER]
blocks["results"]["tab:consistency"] = table(
    "l c c c c c c", "\\textbf{Model} & \\multicolumn{3}{c}{\\textbf{Single-path templates (" + str(N_SINGLE) + ")}} & \\multicolumn{3}{c}{\\textbf{Other templates ("
    + str(N_TEMPLATES - N_SINGLE) + ")}} \\\\\n & \\textbf{All} & \\textbf{Some} & \\textbf{None} & \\textbf{All} & \\textbf{Some} & \\textbf{None}", rows,
    "\\textbf{Consistency Within a Template.} "
    "The share of templates whose 15 instances a model all answers correctly, some but not all, or none, with 95\\% intervals, for the templates that "
    "follow one reasoning path and for the others.", "tab:consistency", resize=True)
rows = [f"{tt(k)} & {f3(Q1[k]['within_sd_quartiles'][1])} & {f3(Q1[k]['within_sd_quartiles'][2])} & {Q1[k]['templates_no_instance_variance']} & "
        f"{Q1[k]['templates_no_variance_all_solved']} & {Q1[k]['templates_no_variance_none_solved']} & "
        + ", ".join(code(t["template"]) for t in Q1[k]["highest_variance_templates"]) for k in ORDER]
blocks["results"]["tab:instance_variance"] = table(
    "l r r r r r l", head("Model", "Median SD", "Upper quartile", "No variance", "All correct", "None correct", "Highest variance"), rows,
    "\\textbf{Instance Variance Within a Template.} "
    "The median and upper quartile over templates of the standard deviation of the score across a template's 15 instances; the number of templates "
    "with no variance at all, split into those solved on every instance and on none; and the templates with the highest variance.",
    "tab:instance_variance", resize=True)

# Depth, repeats, tokens.
bins = list(Q3[ORDER[0]]["by_milestone_count"].items())
rows = [f"{tt(k)} & " + " & ".join(f3(Q3[k]["by_milestone_count"][b]["wrong_rate"]) for b, _ in bins) for k in ORDER]
blocks["results"]["tab:depth"] = table(
    "l " + "r " * len(bins), head("Model", *(f"{b.replace('-', '--')} ({v['items']})" for b, v in bins)), rows,
    "\\textbf{Wrong-Answer Rate Against the Depth of the Gold Derivation.} "
    "The share of instances answered wrong, by the number of milestones in the instance's gold derivation (the number of instances in parentheses).",
    "tab:depth")
rows = [f"{tt(k)} & {REPEATS[k]['items']} & " + " & ".join(f3(REPEATS[k]["scores"][r]) for r in ("repeat1", "repeat2", "repeat3")) +
        f" & {f3(REPEATS[k]['main_on_same_items'])} & {f3(REPEATS[k]['sd'])} & {f3(REPEATS[k]['range'])} & {f3(REPEATS[k]['same_verdict_every_repeat'])}"
        for k in ORDER if k in REPEATS]
blocks["results"]["tab:repeats"] = table(
    "l r r r r r r r r", head("Model", "Instances", "Repeat 1", "Repeat 2", "Repeat 3", "Main run", "SD", "Range", "Same verdict"), rows,
    "\\textbf{Decoding Repeats.} "
    "FAC of four models on 300 instances decoded three more times at the same settings, the main run's FAC on the same instances, the standard "
    "deviation and range across the repeats, and the share of instances with the same verdict in every repeat.", "tab:repeats", resize=True)
rows = [f"{tt(k)} & {thousands(int(REP[k]['median_tokens']))} & {thousands(int(REP[k]['median_tokens_fully_solved']))} & "
        f"{thousands(int(REP[k]['median_tokens_not_fully_solved']))}" for k in ORDER]
blocks["results"]["tab:tokens"] = table(
    "l r r r", head("Model", "All", "Correct", "Not correct"), rows,
    "\\textbf{Completion Tokens Against the Verdict.} Median completion tokens per response, over all responses and by whether the final answer is correct.",
    "tab:tokens", star=False)

# Coverage details and diagnostics.
rows = [f"{tt(k)} & {Q3[k]['readable_wrong_with_milestones']} & {f3(Q3[k]['e3_coverage_on_readable_wrong'])} & {ci(Q3[k]['e3_readable_ci'])} & "
        f"{f3(Q3[k]['e3_null_on_readable_wrong'])} & {f3(Q3[k]['e5_coverage_on_readable_wrong'])} & {ci(Q3[k]['e5_readable_ci'])} & "
        f"{f3(COV[k]['wrong_full_coverage'])} ({ci(COV[k]['wrong_full_coverage_ci'])}) & {f3(COV[k]['solved_low_coverage'])} ({ci(COV[k]['solved_low_coverage_ci'])})"
        for k in ORDER]
blocks["results"]["tab:coverage_wrong"] = table(
    "l r r c r r c c c",
    head("Model", "Wrong", "Matching", "95\\% CI", "Floor", "With judge", "95\\% CI", "Wrong, MC $= 1$", "Correct, MC $< 0.5$"), rows,
    "\\textbf{Milestone Coverage on Wrong Answers.} "
    "For the readable wrong answers with milestones: their number, the coverage by matching alone with its 95\\% interval, the chance floor (the same "
    "response matched against a sibling instance's milestones), and the coverage with the judge; then the share of wrong answers that reach every "
    "milestone and the share of correct answers that reach fewer than half, with 95\\% intervals.", "tab:coverage_wrong", resize=True)
rows = [f"{tt(k)} & {f2(COV[k]['rho_steps'])} & {f2(COV[k]['rho_claims'])} & {int(COV[k]['median_steps'])} & {f3(COV[k]['coverage_fully_solved'])} "
        f"({ci(COV[k]['coverage_fully_solved_ci'])}) & {f3(COV[k]['coverage_wrong'])} ({ci(COV[k]['coverage_wrong_ci'])})" for k in ORDER]
blocks["results"]["tab:coverage_verbosity"] = table(
    "l r r r c c", head("Model", "$\\rho$ (steps)", "$\\rho$ (calc.)", "Median steps", "MC, correct", "MC, wrong"), rows,
    "\\textbf{Coverage Against Verbosity.} "
    "Spearman's $\\rho$ across templates between a template's mean coverage and its mean number of steps and of displayed calculations per response; "
    "the median number of steps; and MC on correct-answer and on wrong-answer responses with 95\\% intervals.", "tab:coverage_verbosity", resize=True)
rows = [f"{tt(k)} & {Q3[k]['fully_solved']} & {f3(Q3[k]['digit_flag_rate_on_fully_solved'])} & {ci(Q3[k]['digit_ci'])} & "
        f"{f3(Q3[k]['tol1_flag_rate_on_fully_solved'])} & {Q3[k]['claims_per_trace']:.2f} & {f3(Q3[k]['traces_with_a_claim'])} & "
        f"{f2(Q3[k]['first_flag_position']['median'])} & {f3(Q3[k]['router_judge_rate_on_fully_solved'])} & {ci(Q3[k]['router_judge_ci'])} & "
        f"{f2(Q3[k]['router_steps_flagged_per_trace'])}" for k in ORDER]
blocks["results"]["tab:diagnostics"] = table(
    "l r r c r r r r r c r",
    head("Model", "Correct", "Arith.\\ flags", "95\\% CI", "At 1\\%", "Calc./resp.", "With a calc.", "First flag", "Step flags", "95\\% CI", "Steps/resp."),
    rows,
    "\\textbf{The Diagnostics on Correct-Answer Responses.} "
    "The number of correct-answer responses; the share with an arithmetic flag and its 95\\% interval; the share at a 1\\% tolerance in place of the "
    "digit rule; the displayed calculations read per response and the share of responses with any; the median relative position of the first flag "
    "(0 is the first step, 1 the last); the share with a judged step flag and its interval; and the steps flagged per response by either check.",
    "tab:diagnostics", resize=True)
rows = [f"{tt(k)} & {Q3[k]['attribution_on_wrong']['traces']} & {f3(Q3[k]['attribution_on_wrong']['digit_rule'])} & "
        f"{f3(Q3[k]['attribution_on_wrong']['e5_missing'])} & {f3(Q3[k]['attribution_on_wrong']['router_judge'])}" for k in ORDER]
blocks["results"]["tab:attribution"] = table(
    "l r r r r", head("Model", "Wrong answers", "Arithmetic flag", "Missing milestone", "Judged step flag"), rows,
    "\\textbf{What Points at a Wrong Answer.} "
    "On the answered wrong answers, the share with an arithmetic flag, with a milestone the judge rules missing, and with a judged step flag; a "
    "response can be in several columns or in none.", "tab:attribution", star=False)

# Rewording.
funnel = [
    f"Instances selected (one in five of every template) & {p_selected}",
    f"Rewrites that pass the scripted checks (at the first, second, third attempt) & {p_passing} ({', '.join(map(str, p_attempts))})",
    f"Instances with no passing rewrite after three attempts & {p_failed}",
    f"Failed attempts by check: technical tokens, numbers, near-copy, part labels, length & {p_checks['tokens']}, {p_checks['numbers']}, {p_checks['copy']}, "
    f"{p_checks['parts']}, {p_checks['length']}",
    f"Passing rewrites with the original's notation restored & {p_restored}",
    f"Pairs the domain experts kept & {r_kept}",
    f"Pairs rejected: not the same problem, adds information, different answer & {r_rejected}: {r_reasons['same=no']}, {r_reasons['clear=yes']}, {r_reasons['answer=no']}",
    "Kept of passing, per branch & " + "; ".join(f"{b.capitalize()} {k} of {k + r}" for b, (k, r) in r_branch.items()),
    f"Templates covered & {q5_templates} ({p_lost} with no passing rewrite, {p_lost_experts} more with no kept pair)",
]
blocks["rewording"]["tab:rewording_funnel"] = table(
    "p{0.72\\columnwidth} r", head("Step", "Count"), funnel,
    "\\textbf{The Rewording Funnel.} From the instances selected to the pairs the domain experts kept.", "tab:rewording_funnel", star=False)
rows = [f"{tt(k)} & {sgn(Q5[k]['diff'])} & {ci(Q5[k]['ci'])} & {ci(Q5[k]['ci90'])} & {'yes' if Q5[k]['within_margin'] else 'no'} & {pv(Q5[k]['p_holm'])} & "
        f"{f3(Q5[k]['detectable'])} & {sgn(Q5[k]['e3']['diff'])} ({ci(Q5[k]['e3']['ci'])}) & {sgn(Q5[k]['e5']['diff'])} ({ci(Q5[k]['e5']['ci'])})" for k in ORDER]
blocks["rewording"]["tab:rewording_results"] = table(
    "l r c c c r r c c",
    head("Model", "FAC change", "95\\% CI", "90\\% CI", f"Within $\\pm {margin:.2f}$", "$p$ (Holm)", "Detectable", "MC, no judge", "MC"), rows,
    "\\textbf{Change Under Rewording.} "
    f"Paraphrase minus original on the {q5_pairs} expert-kept pairs over {q5_templates} templates: the change in FAC with its 95\\% and 90\\% intervals, "
    f"whether the 90\\% interval lies within $\\pm {margin:.2f}$, the Holm-adjusted $p$ of the sign-flip test over templates, the smallest change the "
    "design detects, and the change in MC by matching alone and with the judge (95\\% intervals).", "tab:rewording_results", resize=True)
rows = [f"{tt(k)} & {sgn(vs[k]['paraphrase']['diff'])} ({ci(vs[k]['paraphrase']['ci'])}) & " +
        " & ".join(f"{sgn(vs[k]['repeats'][r]['diff'])} ({ci(vs[k]['repeats'][r]['ci'])})" for r in ("repeat1", "repeat2", "repeat3")) for k in ORDER if k in vs]
blocks["rewording"]["tab:rewording_repeats"] = table(
    "l c c c c", head("Model", "Rewording", "Repeat 1", "Repeat 2", "Repeat 3"), rows,
    "\\textbf{Rewording Against Decoding Noise.} "
    f"For the four models with decoding repeats: the change in FAC under rewording ({q5_pairs} pairs) beside each repeat's change from the main run "
    f"on the {REPEATS['gpt-oss-20b']['items']} repeated instances, with 95\\% intervals.", "tab:rewording_repeats", resize=True)

# Conditions.
rows = []
for a in ARMS:
    base = "the main run" if a["base"] == "main" else "its default"
    tu = a["tool_use"]
    rows.append(f"{CONDITION[a['arm']]} & {tt(a['model'])} & {a['items']} & {f3(a['main_score_on_items'])} & {f3(a['arm_score'])} & {sgn(a['diff'])} & "
                f"{ci(a['ci'])} & {pv(a['p_holm'])} & {f3(a['detectable'])} & {ci(a['ci90'])} & {'yes' if a['within_margin'] else 'no'} & "
                f"{sgn(a['usable_in_both']['diff'])} ({a['usable_in_both']['items']}) & {sgn(a['e5']['diff'])} ({ci(a['e5']['ci'])}) & "
                f"{f3(a['digit_flag_rate_fully_solved']['main'])} / {f3(a['digit_flag_rate_fully_solved']['arm'])} & "
                + (f"{f2(tu['share_with_calls'])}, {f2(tu['calls_per_trace'])}" if tu else "--"))
blocks["conditions"]["tab:conditions"] = table(
    "l l r r r r c r r c c c c c c",
    head("Condition", "Model", "Inst.", "Base", "Cond.", "Change", "95\\% CI", "$p$ (Holm)", "Detectable", "90\\% CI", "Bounded", "Both answered",
         "MC change", "Arith.\\ flags", "Tool use"), rows,
    "\\textbf{The Conditions Against Their Base Run.} "
    "Condition minus base, paired by instance, on the instances of the 450-instance subset the condition covers: FAC in the base run and the "
    "condition, the change with its 95\\% interval, the Holm-adjusted $p$ of the sign-flip test over templates within the condition, the smallest "
    f"change the condition detects, the 90\\% interval and whether it lies within $\\pm {margin:.2f}$, the change on the instances answered in both runs "
    "(their number), the change in MC with its interval, the share of correct-answer responses with an arithmetic flag in the base run and the "
    "condition (unpaired), and for the tool condition the share of responses that call the tool and the calls per response. "
    f"{tt('gpt-5.4')}'s reasoning condition is paired against its default-setting anchor run.", "tab:conditions", resize=True)
rows = []
for x in [fl["deepseek-v4-pro"], fl["gpt-5.4"], flr] + roster_sub:
    cond = {"flagship": "default", "flagship-reasoning-medium": "reasoning, medium", "main": "main run"}[x["arm"]]
    rows.append(f"{tt(x['model'])} & {cond} & {f3(x['score'])} & {ci(x['ci'])} & {f3(x['fully_solved'])} & {x['unusable']} & " +
                " & ".join(f"{f3(x['levels'][lv]['mean'])} ({ci(x['levels'][lv]['ci'])})" for lv in LEVELS) +
                f" & {f3(x['coverage'])} & {ci(x['coverage_ci'])} & {f3(x['digit_flag_rate_fully_solved'])} & {thousands(int(x['median_completion_tokens']))}")
    if x is flr:
        rows.append("\\midrule")
blocks["conditions"]["tab:anchors"] = table(
    "l l r c r r c c c r c r r",
    head("Model", "Run", "FAC", "95\\% CI", "Strict", "No answer", "Easy", "Intermediate", "Advanced", "MC", "95\\% CI", "Arith.\\ flags", "Tokens"), rows,
    "\\textbf{The Flagship Anchors Beside the Evaluated Models on the Same 450 Instances.} "
    "Above the rule, the two anchors at their providers' defaults and \\texttt{GPT-5.4} with reasoning at medium effort; below it, the eleven evaluated "
    "models' main-run responses on the same instances. FAC and strict FAC, the responses with no readable answer, the level means, MC, the share of "
    "correct-answer responses with an arithmetic flag, and the median completion tokens; intervals resample the 150 templates of three instances. "
    "No test is run against an anchor.", "tab:anchors", resize=True)

# Error analysis.
rows = [f"{tt(m)} & " + " & ".join(str(ERR_BUILD["by_model_level"][m][lv]) for lv in LEVELS) + f" & {sum(ERR_BUILD['by_model_level'][m].values())}"
        for m in B2_MODELS]
rows += ["\\midrule", "Branch & \\multicolumn{4}{l}{" + ", ".join(f"{BRANCH[b].lower()} {n}" for b, n in sorted(ERR_BUILD["by_branch"].items(), key=lambda x: -x[1])) + "}"]
blocks["errors"]["tab:error_sample"] = table(
    "l r r r r", head("Model", *LEVELS, "Total"), rows,
    "\\textbf{The Wrong Answers Read.} "
    f"Answered wrong answers per model and level, drawn in proportion to the model's wrong answers at each level and spread over its templates; "
    f"{ERR_BUILD['templates']} templates in all.", "tab:error_sample", star=False)
rows = []
for full, short in CATEGORIES:
    cells = [f"{by_model[m].get(full, 0)} ({maj[m].get(full, 0)})" for m in B2_MODELS]
    rows.append(f"{short} & " + " & ".join(cells))
rows += ["\\midrule", "Readings (wrong answers) & " + " & ".join(f"{sum(by_model[m].values())} ({per_model_items})" for m in B2_MODELS),
         "Fleiss' $\\kappa$ & " + " & ".join(f3(ERR["fleiss_by_model"][m]) for m in B2_MODELS)]
blocks["errors"]["tab:error_by_model"] = table(
    "l r r r r", head("Category", *(NAME[m] for m in B2_MODELS)), rows,
    "\\textbf{Error Categories by Model.} "
    "Readings per category by the three domain experts and, in parentheses, the wrong answers whose majority label is the category; "
    f"Fleiss' $\\kappa$ over the three readers is {f3(fleiss_all)} overall.", "tab:error_by_model", resize=True)
rows = [f"{short} & " + " & ".join(f"{by_level[lv].get(full, 0)} ({pct(by_level[lv].get(full, 0) / level_n[lv])})" for lv in LEVELS)
        for full, short in CATEGORIES]
rows += ["\\midrule", "Readings & " + " & ".join(str(level_n[lv]) for lv in LEVELS)]
blocks["errors"]["tab:error_by_level"] = table(
    "l r r r", head("Category", *LEVELS), rows,
    "\\textbf{Error Categories by Level.} Readings per category and their share of the level's readings.", "tab:error_by_level", star=False)
ASK = {"reading": "Does the wording decide between the closed-system and the flow reading?",
       "form": "Does it decide between the two forms of the truncated virial equation?",
       "unique": "Does the question have one correct answer?", "traces": "Are the responses shown correct under some reading?",
       "data": "Do standard data sources agree within the tolerance?", "method": "Is the gold method the standard one?",
       "trace": "Is the response shown, which lands within 5\\% of the gold, correct?"}
rows = []
for t, v in TEMPLATES_READ.items():
    for ask, answers in v["asks"].items():
        rows.append(f"{code(t)} & {ASK[ask]} & " + "; ".join(f"``{a}'' {n}" for a, n in answers.items()))
blocks["errors"]["tab:template_readings"] = table(
    "l p{0.36\\textwidth} p{0.34\\textwidth}", head("Template", "Question to the three experts", "Answers"), rows,
    "\\textbf{The Experts' Reading of the Templates Whose Wording May Not Pin the Answer.} "
    "Three domain experts of the template's branch answered each question; the two chemical templates (compression work from a truncated virial "
    "equation; adiabatic flame temperature) and six templates on which many wrong answers land within 5\\% of the gold value.",
    "tab:template_readings", size="\\footnotesize")
rows = [f"{code(r[0])} & {r[1]}" for r in near_rows]
blocks["errors"]["tab:exact_digits"] = table(
    "l r", head("Template", "Incorrect verdicts within 0.2\\%"), rows,
    "\\textbf{Incorrect Verdicts Within the Relative Tolerance.} "
    "Over all eleven models, the incorrect verdicts whose stated value lies within 0.2\\% of the target; every one falls on a template whose question "
    "prescribes the digits of the answer, where the check requires them.", "tab:exact_digits", star=False)


# ----------------------------------------------------------------------------------------------- figure
def draw_figure(path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif", "font.size": 6.5, "axes.linewidth": 0.4})
    cats = [c for c in CATEGORIES if c[0] != CATEGORIES[7][0]]  # no "Incomplete" reading was given
    assert all(by_model[m].get(CATEGORIES[7][0], 0) == 0 for m in B2_MODELS)
    cols = [(NAME[m], sum(by_model[m].values()), by_model[m]) for m in B2_MODELS] + [(lv, level_n[lv], by_level[lv]) for lv in LEVELS]
    shares = [[c[2].get(full, 0) / c[1] for c in cols] for full, _ in cats]
    cmap = LinearSegmentedColormap.from_list("blue", ["#ffffff", "#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"])
    fig = plt.figure(figsize=(3.03, 2.55))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 3], height_ratios=[22, 1], wspace=0.06, hspace=0.45, left=0.31, right=0.93, top=0.74, bottom=0.11)
    for j, (sl, title) in enumerate(((slice(0, 4), "By model"), (slice(4, 7), "By level"))):
        ax = fig.add_subplot(gs[0, j])
        data = [row[sl] for row in shares]
        ax.imshow(data, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(cols[sl])))
        ax.set_xticklabels([c[0] for c in cols[sl]], rotation=45, ha="left", rotation_mode="anchor")
        ax.xaxis.tick_top()
        ax.set_xlabel(title, fontsize=6.5, color="#52514e", labelpad=3)
        if j == 0:
            ax.set_yticks(range(len(cats)))
            ax.set_yticklabels([s for _, s in cats])
        else:
            ax.set_yticks([])
        for i, row in enumerate(data):
            for k, v in enumerate(row):
                ax.text(k, i, f"{v * 100:.0f}", ha="center", va="center", fontsize=6.2, color="#ffffff" if v > 0.5 else "#0b0b0b")
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0, pad=2)
        ax.set_xticks([x - 0.5 for x in range(1, len(cols[sl]))], minor=True)
        ax.set_yticks([y - 0.5 for y in range(1, len(cats))], minor=True)
        ax.grid(which="minor", color="#ffffff", linewidth=1.2)
        ax.tick_params(which="minor", length=0)
    cax = fig.add_subplot(gs[1, :])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 100)), cax=cax, orientation="horizontal")
    cb.set_label("Share of readings (%)", fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)
    cb.outline.set_linewidth(0.4)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- write and check
def strip_generated(tex: str) -> str:
    return re.sub(r"% BEGIN GENERATED (\S+) .*?% END GENERATED \1", " ", tex, flags=re.S)


def prose(tex: str) -> str:
    """The text without comments, generated blocks, graphics, labels, references, citations, or model names."""
    tex = re.sub(r"(?<!\\)%.*", "", strip_generated(tex))
    tex = re.sub(r"\\(texttt|label|autoref|citep|citet|ref|includegraphics)(\[[^\]]*\])?\{[^}]*\}", " ", tex)
    return re.sub(r"\s+", " ", tex)


def numbers(text: str) -> set[str]:
    return {n.rstrip(",") for n in NUMBER.findall(text)}


def flat(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def write_blocks() -> None:
    for key, path in FILES.items():
        tex = path.read_text(encoding="utf-8")
        for name, body in blocks[key].items():
            pattern = re.compile(r"% BEGIN GENERATED " + re.escape(name) + r" .*?% END GENERATED " + re.escape(name), re.S)
            if not pattern.search(tex):
                raise SystemExit(f"{path.name} has no markers for {name}")
            tex = pattern.sub(lambda m: block(name, body), tex, count=1)
        path.write_text(tex, encoding="utf-8")
    draw_figure(FIG)


def check() -> int:
    bib_keys = set(re.findall(r"@\w+\{([^,\s]+),", BIB.read_text(encoding="utf-8")))
    labels = set()
    for p in SRC.rglob("*.tex"):
        labels |= set(re.findall(r"\\label\{([^}]*)\}", p.read_text(encoding="utf-8")))
    own = {"main": phrases, **appendix_phrases}
    failures = 0
    for key, path in FILES.items():
        tex = path.read_text(encoding="utf-8")
        for name, body in blocks[key].items():
            m = re.search(r"% BEGIN GENERATED " + re.escape(name) + r" .*?% END GENERATED " + re.escape(name), tex, re.S)
            if not m or flat(m.group(0)) != flat(block(name, body)):
                print(f"STALE block in {path.name}: {name}")
                failures += 1
        body = flat(re.sub(r"(?<!\\)%.*", "", strip_generated(tex)))
        missing = [p for p in own[key] if flat(p) not in body]
        known = numbers(prose(" ".join(own[key])))
        words = {w.lower() for w in NUMBER_WORDS.findall(prose(" ".join(own[key])))}
        stray = sorted(numbers(prose(tex)) - known) + sorted({w.lower() for w in NUMBER_WORDS.findall(prose(tex))} - words)
        cited = {k.strip() for c in re.findall(r"\\cite[pt]\{([^}]*)\}", tex) for k in c.split(",")}
        unresolved = sorted(cited - bib_keys)
        refs = {r for r in re.findall(r"\\autoref\{([^}]*)\}", tex)}
        undefined = sorted(refs - labels)
        for p in missing:
            print(f"MISSING in {path.name}: {p}")
        for n in stray:
            print(f"NOT GENERATED in {path.name}: {n}")
        for k in unresolved:
            print(f"UNRESOLVED citation in {path.name}: {k}")
        for r in undefined:
            print(f"UNDEFINED label in {path.name}: {r}")
        long_lines = [i + 1 for i, l in enumerate(strip_generated(tex).split("\n")) if len(l) > 100]
        if long_lines:
            print(f"LINES over 100 characters in {path.name}: {long_lines}")
        print(f"{path.name}: {len(blocks[key])} generated blocks current; {len(own[key]) - len(missing)} of {len(own[key])} phrases present; "
              f"{len(stray)} numbers not generated; {len(cited) - len(unresolved)} of {len(cited)} citation keys resolve; "
              f"{len(refs) - len(undefined)} of {len(refs)} references defined")
        failures += len(missing) + len(stray) + len(unresolved) + len(undefined) + len(long_lines)
    if not FIG.exists():
        print(f"MISSING figure {FIG}")
        failures += 1
    return failures


if __name__ == "__main__":
    if "--write" in sys.argv:
        write_blocks()
        print(f"wrote {sum(len(b) for b in blocks.values())} generated blocks and {FIG.relative_to(REPO)}")
    elif "--check" in sys.argv:
        sys.exit(1 if check() else 0)
    else:
        for key, ph in (("6_results.tex", phrases), *((f"appendices/{k}.tex", v) for k, v in appendix_phrases.items() if k != "main")):
            print(f"% phrases {key} must contain")
            print("\n".join(ph))
