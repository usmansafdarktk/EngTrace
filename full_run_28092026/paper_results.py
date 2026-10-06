"""Generate the tables, the figures and the numbers of the paper's Results and Error Analysis (Sections 5.3 and 5.4) and
their appendices, and check them in the source.

    python full_run_28092026/paper_results.py           # print the phrases the prose must contain
    python full_run_28092026/paper_results.py --write   # rewrite every generated block in the .tex files and draw the figures
    python full_run_28092026/paper_results.py --write --text-only   # rewrite the blocks and leave the figures as drawn
    python full_run_28092026/paper_results.py --check   # exit 1 unless every generated block is current, every phrase is in
                                                        # its file, every number in the prose is a phrase's, every citation
                                                        # key resolves, every \\autoref label is defined and every figure exists

Files written: overleaf_source_04102026/6_results.tex (Table 1, the two main-text figures and the prose phrases),
appendices/results.tex, appendices/branch_domain.tex (the branch and domain figures and table), appendices/paraphrase.tex,
appendices/further_experiments.tex, appendices/error_analysis.tex (their tables, between "% BEGIN GENERATED <name>" and
"% END GENERATED <name>" markers), and the four placed figures under figs/.

Sources, each written by a committed script: results/results.json (analyze.py: every score, interval, test and
condition), expert_request/scored.json and build.json (expert_kits.py: the experts' readings of wrong answers and of
the templates, counts only), RESIDUAL_INCORRECT.md (residual_incorrect.py: the remaining incorrect verdicts by their
distance to the target), PARAPHRASE.md and PARAPHRASE_REVIEW.md (paraphrase.py, paraphrase_kit.py --score: the
paraphrase funnel), FLAG_REVIEW_3.md (flag_sample.py: the expert reading of the arithmetic flags), the shortcut
template list of analyze.py, the tool and paraphrase constants of run_traces.py and paraphrase.py, and the judged
step check's precision row of appendices/validation.tex (docs/appendix_evaluation.py). Model display names are those
of paper_setup.py.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import re
import sys
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SRC = REPO / "overleaf_source_04102026"
APPX = SRC / "appendices"
FIGS = SRC / "figs"
FILES = {"main": SRC / "6_results.tex", "results": APPX / "results.tex", "branch_domain": APPX / "branch_domain.tex",
         "paraphrase": APPX / "paraphrase.tex", "experiments": APPX / "further_experiments.tex", "errors": APPX / "error_analysis.tex"}
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
NO_REASONING_ANCHOR = {"gpt-5.4"}  # the flagship anchor whose default returned no reasoning tokens (D-182)
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
ORDINAL = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh"]
NUMBER = re.compile(r"(?<![\w.\-])\d[\d,]*(?:\.\d+)?")
NUMBER_WORDS = re.compile(r"\b(two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)\b", re.I)

# Figure style: one hue, text in ink, recessive axes; marker fill is the second encoding, so the figures read in greyscale.
BLUE, DARK, LIGHT, INK, MUTED, BAND = "#2a78d6", "#0d366b", "#cde2fb", "#0b0b0b", "#52514e", "#eceae6"
HATCH = "#d3d3d3"  # the hatch lines of the error-category figure
COLUMN = 3.03  # the ACL column width in inches
SERIF = {"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 7, "legend.fontsize": 7}  # the figures' type
FIG_NAME = {**NAME, "gpt-oss-20b": "GPT OSS 20B"}  # model labels in the figures, every name capitalized


def wrap(name: str, width: int = 15) -> str:
    """A figure label over two lines, broken at its last space, when it is longer than width characters."""
    return name[::-1].replace(" ", "\n", 1)[::-1] if len(name) > width and " " in name else name


# ----------------------------------------------------------------------------------------------- formatting
def f3(x: float) -> str:
    return f"{x:.3f}"


def f2(x: float) -> str:
    return f"{x:.2f}"


def f1(x: float) -> str:
    return f"{x:.1f}"


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


def rng(xs, fmt=f3, join=" to ") -> str:
    return f"{fmt(min(xs))}{join}{fmt(max(xs))}"


def block(name: str, body: str) -> str:
    return MARK.format(name=name, body=body.rstrip("\n"))


def caption_lines(caption: str) -> str:
    """The caption wrapped at 100 characters, never inside a model name."""
    guarded = re.sub(r"\\texttt\{[^}]*\}", lambda m: m.group(0).replace(" ", "\x00"), caption)
    lines = textwrap.wrap(guarded, width=100, initial_indent=" " * 9, break_long_words=False, break_on_hyphens=False)
    return "\\caption{" + "\n".join(lines)[9:].replace("\x00", " ") + "}"


def table(spec: str, header: str, rows: list[str], caption: str, label: str, star: bool = True,
          size: str = "\\small", resize: bool = False, placement: str = "t", shade_header: bool = True, aliases: tuple = ()) -> str:
    """A booktabs table. A row that is \\midrule, a \\cmidrule, or already ends with \\\\ is emitted as it is (group rows
    carry their own \\rowcolor and \\\\); resize shrinks a wide table to the text width and never enlarges it."""
    env = "table*" if star else "table"
    lines = [f"\\begin{{{env}}}[{placement}]", "\\centering", size, "\\renewcommand{\\arraystretch}{1.1}"]
    if resize:
        lines.append("\\adjustbox{max width=\\textwidth}{%")
    lines += [f"\\begin{{tabular}}{{{spec}}}", "\\toprule"] + (["\\rowcolor{gray!10}"] if shade_header else []) + [header + " \\\\", "\\midrule"]
    lines += [r if r == "\\midrule" or r.startswith("\\cmidrule") or r.endswith("\\\\") else r + " \\\\" for r in rows]
    lines += ["\\bottomrule", "\\end{tabular}"]
    if resize:
        lines.append("}")
    lines += [caption_lines(caption), f"\\label{{{label}}}"] + [f"\\label{{{a}}}" for a in aliases] + [f"\\end{{{env}}}"]
    return "\n".join(lines)


def figure(path: str, caption: str, label: str, star: bool = False, width: str | None = None) -> str:
    env, default = ("figure*", "\\textwidth") if star else ("figure", "\\columnwidth")
    return "\n".join([f"\\begin{{{env}}}[t]", "    \\centering", f"    \\includegraphics[width={width or default}]{{figs/{path}}}",
                      "    " + caption_lines(caption).replace("\n", "\n    "), f"    \\label{{{label}}}", f"\\end{{{env}}}"])


def head(*cells: str) -> str:
    return " & ".join("\\textbf{" + c + "}" for c in cells)


def mk(*lines: str) -> str:
    """A bold header cell over several lines."""
    return "\\textbf{\\makecell{" + "\\\\".join(lines) + "}}"


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
paraphrase_md = (HERE / "PARAPHRASE.md").read_text(encoding="utf-8")
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
N_ITEMS = sum(Q1[ORDER[0]][k] for k in ("correct", "partial", "incorrect", "unusable"))
N_TEMPLATES = Q2[ORDER[0]]["templates"]
N_PAIRS = len(PAIRS)
N_SINGLE = Q4[ORDER[0]]["single_path"]["templates"]
assert all(sum(m[k] for k in ("correct", "partial", "incorrect", "unusable")) == N_ITEMS for m in Q1.values())
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
assert not any(frozenset(p) in SEP_FAC for p in itertools.combinations(TOP5, 2))
nonsig = [p for p in PAIRS if p["p_holm"] >= 0.05]
assert all(abs(p["diff"]) < p["detectable"] for p in nonsig)
strict_agree = sum((p["p_holm"] < 0.05) == (p["fully_p_holm"] < 0.05) for p in PAIRS)
scores = [Q1[k]["score"] for k in ORDER]
spread5 = max(scores[:5]) - min(scores[:5])
sixth, rest = ORDER[5], ORDER[6:]
sixth_sep = [k for k in TOP5 if frozenset((sixth, k)) in SEP_FAC]
assert all(frozenset((a, b)) in SEP_FAC for a in rest for b in ORDER[:6])
no_variance = [Q1[k]["templates_no_instance_variance"] for k in ORDER]

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
n_branch_pairs = sum(len(BL[k]["pairs"]) for k in ORDER)
lowest_branch = {k: min(BL[k]["branch"], key=lambda b: BL[k]["branch"][b]["mean"]) for k in ORDER}
assert len(set(lowest_branch.values())) > 1
lowest_domain = {k: min(REP[k]["domain"], key=REP[k]["domain"].get) for k in ORDER}
thermo_models = [k for k in ORDER if lowest_domain[k] == "thermodynamics"]
thermo = [REP[k]["domain"]["thermodynamics"] for k in ORDER]

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
assert not any(frozenset(p) in SEP_MC for p in itertools.combinations(top3_cov, 2))
rho = [COV[k]["rho_steps"] for k in ORDER]

REPRESENTATIVE = ["deepseek-v4.1-flash", "claude-sonnet-5", "gpt-5.4-mini", "gpt-oss-20b"]  # the four models the figures show
top_rep, low_rep = REPRESENTATIVE[:2], REPRESENTATIVE[2:]
assert top_rep == [first_fac, first_cov] and low_rep[0] in NO_REASONING and low_rep[1] == ORDER[-1] and all(k in gap_sig for k in low_rep)
assert all(lowest_domain[k] == "thermodynamics" for k in top_rep)
rep_branch = {k: [BL[k]["branch"][b]["mean"] for b in BRANCH] for k in REPRESENTATIVE}
top_branch = rep_branch[top_rep[0]] + rep_branch[top_rep[1]]
top_domain = [v for k in top_rep for v in REP[k]["domain"].values()]
low_domain = {k: (min(REP[k]["domain"], key=REP[k]["domain"].get), min(REP[k]["domain"].values()), max(REP[k]["domain"].values()))
              for k in low_rep}
with (HERE / "results/per_template.csv").open(encoding="utf-8") as f:
    domain_templates = Counter(r["domain"] for r in csv.DictReader(f) if r["model"] == ORDER[0])
n_easy, n_adv = BL[ORDER[0]]["level"]["Easy"]["templates"], BL[ORDER[0]]["level"]["Advanced"]["templates"]

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
flag_all = next(r for r in md_table(flags, r"\| model \| flags drawn") if r[0] == "all")
flags_read, flags_slip, flags_unsure = int(flag_all[2]), int(flag_all[3]), int(flag_all[5])
flags_decided = flags_read - flags_unsure
flag_precision = flags_slip / flags_decided
assert f3(flag_precision) == flag_all[6]
router_pr = (re.search(r"Judged step and arithmetic checks\S* & P / R, all steps; correct answers & [\d.]+ / [\d.]+; ([\d.]+) / ([\d.]+)",
                       validation_tex)  # the row as validation.tex prints it now
             or re.search(r"Judged step and arithmetic checks & steps, correct answers & P / R & ([\d.]+) / ([\d.]+)", validation_tex))
assert router_pr, "the judged step check's precision row is missing from appendices/validation.tex"
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
para = md_kv(paraphrase_md)
rev = md_kv(review)
p_selected = int(para["items selected (subsamples.PARAPHRASE)"])
p_passing, p_failed = int(para["passing a paraphrase"]), int(para["no paraphrase after 3 attempts"])
p_attempts = [int(x) for x in para["passing at attempt 1, 2, 3"].split(", ")]
p_restored = int(para["passing ones with notation restored (Unicode the original lacks, put back as written)"].split()[0])
p_lost_templates = para["templates with no paraphrase for any of their items"].split(": ")[1].split(", ")
p_lost = int(para["templates with no paraphrase for any of their items"].split(":")[0])
p_checks = dict(re.findall(r"(\w+) (\d+)", para["failed attempts by check"]))
r_kept, r_rejected, r_returned = int(rev["kept"]), int(rev["rejected"]), int(rev["returned"])
r_branch = {b: (int(k), int(r)) for b, k, r in re.findall(r"(\w+) (\d+) / (\d+)", rev["kept / rejected, per branch"])}
r_reasons = dict((k, int(v)) for k, v in re.findall(r"(\w+=\w+) (\d+)", rev["answers behind the rejections"]))
assert p_selected == 450 and p_passing == r_returned == 316 and r_kept == q5_pairs == 277 and r_rejected == 39
assert p_lost == len(p_lost_templates) and p_selected - p_passing == p_failed
p_lost_experts = N_TEMPLATES - q5_templates - p_lost
least_branch = min(r_branch, key=lambda b: r_branch[b][0])
below90 = [k for k in ORDER if Q5[k]["ci90"][1] < 0]
above90 = [k for k in ORDER if Q5[k]["ci90"][0] > 0]
vs = res["q5"]["vs_repeats"]
vs_within = [k for k in ORDER if k in vs and vs[k]["abs_within_repeat_spread"]]
vs_beyond = [k for k in ORDER if k in vs and not vs[k]["abs_within_repeat_spread"]]
rep_sd = [REPEATS[k]["sd"] for k in REPEATS]
rep_same = [REPEATS[k]["same_verdict_every_repeat"] for k in REPEATS]
rep_items = REPEATS["gpt-oss-20b"]["items"]
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
n_sub = ANCH["flagship"]["items"]
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
CALC, FORM, NOERR, INCOMPLETE = CATEGORIES[5][0], CATEGORIES[2][0], CATEGORIES[6][0], CATEGORIES[7][0]
assert all(by_model[m].get(INCOMPLETE, 0) == 0 for m in B2_MODELS)
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
fleiss_models = list(ERR["fleiss_by_model"].values())

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


# ----------------------------------------------------------------------------------------------- derived for the text
TWO_CHEMICAL = ["template_work_isothermal_virial", "template_adiabatic_flame_temperature"]  # wording does not pin the answer
norm = lambda t: t.replace("template_", "")
with (HERE / "results/per_template.csv").open(encoding="utf-8") as f:
    PT = list(csv.DictReader(f))
assert all(r["domain"] == "thermodynamics" for r in PT if r["template_id"] in TWO_CHEMICAL)
_dom_wo = defaultdict(lambda: defaultdict(list))
for r in PT:
    if r["template_id"] not in TWO_CHEMICAL:
        _dom_wo[r["model"]][r["domain"]].append(float(r["answer_score"]))
lowest_domain_wo = {k: min(_dom_wo[k], key=lambda d: sum(_dom_wo[k][d]) / len(_dom_wo[k][d])) for k in ORDER}
n_thermo_wo = sum(v == "thermodynamics" for v in lowest_domain_wo.values())
assert n_thermo_wo == 1
thermo_eight = [REP[k]["domain"]["thermodynamics"] for k in thermo_models]
glm_empty_adv = sum(int(r["unusable"]) for r in PT if r["model"] == "glm-5.3" and r["level"] == "Advanced")
assert sum(int(r["unusable"]) for r in PT if r["model"] == "glm-5.3") == Q1["glm-5.3"]["empty"]
branch_of = {r["domain"]: r["branch"] for r in PT}
geo = REP["gpt-oss-20b"]["domain"]["geotechnical_engineering"]
civil_other = sorted(REP["gpt-oss-20b"]["domain"][d] for d in branch_of if branch_of[d] == "civil_engineering" and d != "geotechnical_engineering")
assert len(civil_other) == 2 and geo < civil_other[0] and lowest_domain["gpt-oss-20b"] == "geotechnical_engineering"
n_branch_templates = BL[ORDER[0]]["branch"]["chemical_engineering"]["templates"]
detect_branch = [BL[k]["detectable_branch"] for k in ORDER]
n_int = BL[ORDER[0]]["level"]["Intermediate"]["templates"]
dec_weights = {d["model_key"]: d["weights"] for d in load(HERE / "results/decoding_table.json")}
OPEN = [k for k in ORDER if dec_weights[k] == "open"]
CLOSED = [k for k in ORDER if dec_weights[k] == "closed"]
assert len(OPEN) + len(CLOSED) == len(ORDER) and sum(k in NO_REASONING for k in rest) == 4
tol_swaps = {}
for var in ("half_tol", "double_tol"):
    order_var = sorted(ORDER, key=lambda m: -SENS[m][var])
    tol_swaps[var] = [(a, b) for i, a in enumerate(ORDER) for b in ORDER[i + 1:] if order_var.index(a) > order_var.index(b)]
    assert all(frozenset(pair) not in SEP_FAC for pair in tol_swaps[var])
n_tol_swaps = len(tol_swaps["half_tol"])
assert n_tol_swaps == len(tol_swaps["double_tol"]) and sens_tau["half_tol"] < 1
mini_gap_top5 = [x - r_mini["arm_score"] for x in sub_scores]
assert min(mini_gap_top5) > 0
mc_correct = {m: COV[m]["coverage_fully_solved"] for m in (first_fac, first_cov)}
assert all(COV[k]["rho_steps"] < 0 for k in ORDER)
n_mcnemar_sig = sum(pp["mcnemar_p_holm"] < 0.05 for pp in PAIRS)
n_strict_sig = sum(pp["fully_p_holm"] < 0.05 for pp in PAIRS)
strict_flip = [pp for pp in PAIRS if (pp["p_holm"] < 0.05) != (pp["fully_p_holm"] < 0.05)]
assert len(strict_flip) == 1 and strict_agree == N_PAIRS - 1
top5_detect = [Q2[k]["detectable_planned"] for k in TOP5]
assert fl["gpt-5.4"]["score"] < min(x["ci"][0] for x in top5_sub)
ml, gl = r_mini["levels"], r_gem["levels"]
assert ml["Advanced"]["arm"] - ml["Advanced"]["main"] > (ml["Easy"]["arm"] - ml["Easy"]["main"])
mini_flag_main, mini_flag_arm = r_mini["digit_flag_rate_fully_solved"]["main"], r_mini["digit_flag_rate_fully_solved"]["arm"]
flag_prec = [float(r[6]) for r in md_table(flags, r"\| model \| flags drawn") if r[0] != "all"]
assert len(flag_prec) == len(ORDER) and router_recall < 0.5
HALL, SETUP = CATEGORIES[0][0], CATEGORIES[1][0]
conceptual = {m: sum(maj[m].get(c, 0) for c in (HALL, SETUP, FORM)) for m in B2_MODELS}
assert conceptual["gpt-oss-20b"] > max(conceptual[m] for m in others if m != "gpt-oss-20b")
assert maj["gpt-oss-20b"].get(CALC, 0) < per_model_items / 2 < min(maj[m].get(CALC, 0) for m in others if m != "gpt-oss-20b")
_votes = defaultdict(Counter)
for n in ERR["notes"]:
    _votes[(n["model"], n["code"], norm(n["template"]))][n["category"]] += 1
_claude_noerr_t = Counter(t for (m, c, t), v in _votes.items() if m == "claude-sonnet-5" and v.most_common(1)[0][1] >= 2 and v.most_common(1)[0][0] == NOERR)
assert sum(_claude_noerr_t.values()) == claude_noerr
claude_noerr_top3 = sum(n for _, n in _claude_noerr_t.most_common(3))
_top3 = [t for t, _ in _claude_noerr_t.most_common(3)]
assert sorted(t for t in _top3 if t in map(norm, TWO_CHEMICAL)) == sorted(map(norm, TWO_CHEMICAL)) and sum(t in map(norm, SYMBOLIC) for t in _top3) == 1
top5_empty = sum(Q1[k]["unusable"] for k in TOP5)
assert top5_empty == sum(Q1[k]["empty"] for k in TOP5)
top5_partial = sum(Q1[k]["partial"] for k in TOP5)
assert top5_partial % 2 == 0
top5_partial_points = top5_partial // 2
lost_points = top5_empty + top5_partial_points + top5_incorrect
assert lost_points == round(sum(15 * N_TEMPLATES * (1 - Q1[k]["score"]) for k in TOP5))
_vir = re.search(r"\| incorrect traces \|[^\n]*\n\|[-:| ]+\n\| (\d+) \| (\d+) \|", residual)
virial_n, virial_flow = int(_vir.group(1)), int(_vir.group(2))
median_sd_zero = sum(Q1[k]["within_sd_quartiles"][1] == 0 for k in ORDER)
claude_adv_readings = ERR_BUILD["by_model_level"]["claude-sonnet-5"]["Advanced"] * ERR_BUILD["readers"] if isinstance(ERR_BUILD["readers"], int) else ERR_BUILD["by_model_level"]["claude-sonnet-5"]["Advanced"] * 3
solved = {k: Q1[k]["templates_no_variance_all_solved"] / N_TEMPLATES for k in ORDER}
e3w_max_model = max(ORDER, key=lambda k: Q3[k]["e3_coverage_on_readable_wrong"])
B2_TEXT = ["claude-sonnet-5", "gpt-5.4-mini", "gemma-4-26b-a4b", "gpt-oss-20b"]  # the order the text names them
assert set(B2_TEXT) == set(B2_MODELS)
maj_calc = {m: maj[m].get(CALC, 0) for m in B2_MODELS}

# ----------------------------------------------------------------------------------------------- phrases
low_pairs = [pr for pr in itertools.combinations(rest, 2) if frozenset(pr) in SEP_FAC]
assert len(low_pairs) == 1
low_pair = low_pairs[0]
assert separated(PAIRS, "fully_p_holm") < SEP_FAC and len(SEP_FAC) - n_strict_sig == 1
assert sum(k in NO_REASONING for k in gap_sig) == 3 and maj_calc["gpt-5.4-mini"] == maj_calc["gemma-4-26b-a4b"]
SENS_OTHER = ["half_unit", "whole_trace", "without_shortcut_templates", "without_symbolic_templates"]
other_shift = max(abs(SENS[k][c] - SENS[k]["fitted"]) for k in ORDER for c in SENS_OTHER)
other_tau = min(sens_tau[c] for c in SENS_OTHER)
assert sens_tau["unusable_excluded"] == min(sens_tau.values()) and other_tau > sens_tau["half_tol"] > sens_tau["unusable_excluded"]
VIRIAL, FLAME = (TEMPLATES_READ[t]["asks"] for t in TWO_CHEMICAL)
assert VIRIAL == {"reading": {"either: the wording does not decide": 3}, "form": {"either: the wording does not decide": 3},
                  "unique": {"no": 3}, "traces": {"all of them": 3}}
assert FLAME == {"data": {"no: standard sources differ by more than that": 2, "only if the same data source is used": 1},
                 "method": {"no": 1, "yes": 2}, "traces": {"some of them": 3}}
assert all(TEMPLATES_READ[t]["asks"] == {"unique": {"yes": 3}, "trace": {"no": 3}} for t in near_templates)

n_open_top = sum(dec_weights[k] == "open" for k in TOP5)
assert n_open_top == 4 and dec_weights[ORDER[0]] == "open"
max_none = max(Q1[k]["templates_no_variance_none_solved"] for k in ORDER) / N_TEMPLATES
assert all(frozenset((a, b)) in SEP_FAC for a in rest for b in ORDER[:6]) and round(margin * 100) == 5
assert max(flr["score"], fl["deepseek-v4-pro"]["score"], fl["gpt-5.4"]["score"]) <= max(x["ci"][1] for x in top5_sub)  # no anchor run exceeds the top tier

# Derived for the rewritten text (D-195).
top5_pair_detect = [p["detectable"] for p in PAIRS if {p["a"], p["b"]} <= set(TOP5)]
assert len(top5_pair_detect) == 10
detect_lo, detect_hi = round(min(top5_pair_detect) * 100), round(max(top5_pair_detect) * 100)
assert 0 < detect_lo <= detect_hi <= 12
n_lower_noreason = sum(k in NO_REASONING for k in rest)
assert n_lower_noreason == 4
assert set(many_wrong) == set(rest)  # the five models with more than 100 wrong answers are the lower five
many_e3w = [Q3[k]["e3_coverage_on_readable_wrong"] for k in many_wrong]
many_floor = [Q3[k]["e3_null_on_readable_wrong"] for k in many_wrong]
many_missing = [Q3[k]["attribution_on_wrong"]["e5_missing"] for k in many_wrong]
many_router = [Q3[k]["attribution_on_wrong"]["router_judge"] for k in many_wrong]
claude_attr = Q3["claude-sonnet-5"]["attribution_on_wrong"]
assert claude_attr["e5_missing"] < min(many_missing) and claude_attr["router_judge"] < min(many_router)
low_some = [Q4[k]["single_path"]["some"] for k in rest]
top_some = [Q4[k]["single_path"]["some"] for k in strong6]
assert max(top_some) < min(low_some)
depth6_top5 = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in TOP5]
depth6_low4 = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in gap_sig]
assert max(depth6_top5) < min(depth6_low4)
hall_readings = sum(by_model[m].get(HALL, 0) for m in B2_MODELS)
remaining_incorrect = top5_incorrect - form_total
assert remaining_incorrect > 0
n_top5_responses = len(TOP5) * N_ITEMS
r_reject_share = r_rejected / r_returned
gpt54_flags = (fl["gpt-5.4"]["digit_flag_rate_fully_solved"], flr["digit_flag_rate_fully_solved"])
assert gpt54_flags[1] < gpt54_flags[0] and fl["gpt-5.4"]["score"] < flr["score"] and "gpt-5.4" in NO_REASONING_ANCHOR
mini_mc = r_mini["e5"]["diff"]
assert mini_mc > 0 and ml["Advanced"]["diff"] > max(ml["Easy"]["diff"], ml["Intermediate"]["diff"])
tool_share = {m: tool[m]["tool_use"]["share_with_calls"] for m in tool}
mini_tool_flags = tool["gpt-5.4-mini"]["digit_flag_rate_fully_solved"]
assert mini_tool_flags["arm"] < mini_tool_flags["main"] and mini_flag_arm < mini_flag_main
assert all(ob[m]["e5"]["diff"] > 0 for m in ob)  # MC rises for all three under the open book

# Tier uniformity across branches and domains, and the readings the experiments' sentences rest on (D-195, third pass).
def f2floor(x: float) -> str:
    return f"{math.floor(x * 100) / 100:.2f}"


def f2ceil(x: float) -> str:
    return f"{math.ceil(x * 100) / 100:.2f}"


branch_span = {k: max(BL[k]["branch"][b]["mean"] for b in BRANCH) - min(BL[k]["branch"][b]["mean"] for b in BRANCH) for k in ORDER}
top_span, low_span = [branch_span[k] for k in TOP5], [branch_span[k] for k in rest]
assert max(top_span) < min(low_span)
dom_min = {k: min(REP[k]["domain"].values()) for k in ORDER}
top_dom_floor = f2floor(min(dom_min[k] for k in TOP5))  # no top-tier domain mean falls below this
low_dom_ceil = f2ceil(max(dom_min[k] for k in rest))  # every lower-tier model has a domain mean below this
assert all(dom_min[k] >= float(top_dom_floor) for k in TOP5) and all(dom_min[k] < float(low_dom_ceil) for k in rest)
assert min(many_missing) > 0.5  # "most of the lower tier's wrong answers"
flag_arm = arm[("flagship-reasoning-medium", "gpt-5.4")]
assert 0.4 <= mini_flag_arm / mini_flag_main <= 0.6  # "roughly halve"
assert flag_arm["digit_flag_rate_fully_solved"]["arm"] < 0.5 * flag_arm["digit_flag_rate_fully_solved"]["main"]  # "sheds most"
assert ORDER[-1] == "gpt-oss-20b" and 0.6 <= tool_share["claude-sonnet-5"] <= 0.7  # "the weakest model"; "two thirds"
assert r_mini["diff"] > max(scores[:5]) - min(scores[:5]) and r_mini["diff"] > max(scores[6:]) - min(scores[6:])
assert all(Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] > Q3[k]["by_milestone_count"]["1"]["wrong_rate"] for k in ORDER)

_lc = Counter(BRANCH[b].lower() for b in lowest_branch.values())
_lo = sorted(_lc, key=lambda b: (-_lc[b], b))
lowest_list = ", ".join((f"{b} for {WORD[_lc[b]]} models" if i == 0 else f"{'and ' if i == len(_lo) - 1 else ''}{b} for {WORD[_lc[b]]}")
                       for i, b in enumerate(_lo))  # e.g. "chemical for five models, electrical for three, ..., and civil for one"

assert all(tool[m]["tool_use"]["score_with_calls"] < tool[m]["tool_use"]["score_without_calls"] for m in tool)  # the tool is called on harder instances

# The error readings the appendix states beyond the main text.
UNIT, SIGN = CATEGORIES[3][0], CATEGORIES[4][0]
unit_total = sum(by_model[m].get(UNIT, 0) for m in B2_MODELS)
sign_total = sum(by_model[m].get(SIGN, 0) for m in B2_MODELS)
n_sign_models = sum(by_model[m].get(SIGN, 0) > 0 for m in B2_MODELS)
hall_gptoss = by_model["gpt-oss-20b"].get(HALL, 0)
assert unit_total <= 2 and hall_gptoss > hall_readings / 2  # "nearly absent"; "most hallucination readings are gpt-oss-20b's"
k_by = ERR["fleiss_by_model"]
k_low_model = min(k_by, key=k_by.get)
assert k_low_model == "gpt-oss-20b" and conceptual["gpt-oss-20b"] == max(conceptual.values())
k_others = [v for m, v in k_by.items() if m != k_low_model]
adv_noerr = by_level["Advanced"].get(NOERR, 0)
adv_noerr_share = adv_noerr / level_n["Advanced"]
_claude_non_adv_readings = sum(ERR_BUILD["by_model_level"]["claude-sonnet-5"][lv] for lv in ("Easy", "Intermediate")) * 3
assert by_model["claude-sonnet-5"][NOERR] - _claude_non_adv_readings > adv_noerr / 2  # "comes largely from Claude Sonnet 5"
resid3 = {t: int(tpl[t][1]) for t in ("template_work_isothermal_virial", "template_ber_estimation_mary", "template_adiabatic_flame_temperature")}
assert sorted(resid3.values(), reverse=True) == sorted(int(r[1]) for r in tpl.values())[-3:][::-1]  # the three largest
resid3_total = sum(resid3.values())

phrases = [  # 6_results.tex must contain each of these, whitespace aside
    # 5.3.1
    f"Five models lie within {f3(spread5)} of one another, and no pair among them differs after Holm's correction, which tightens the "
    f"{N_PAIRS} pairwise tests so that the chance of any false positive among them stays at most 5\\%",
    f"The design detects differences of {WORD[detect_lo]} to {WORD[detect_hi]} points between them, and the lower {WORD[len(rest)]} each "
    f"differ from all {WORD[6]} models above them",
    f"{WORD[n_open_top].capitalize()} of the top {WORD[len(TOP5)]} are open-weights, and {WORD[n_lower_noreason]} of the lower {WORD[len(rest)]} "
    "return no reasoning tokens at their providers' defaults",
    f"almost no template defeats a model on all 15 of its instances, but the top {WORD[6]} solve {pct(min(solved[k] for k in ORDER[:6]))} to "
    f"{pct(max(solved[k] for k in ORDER[:6]))} of the templates on every instance and the lower {WORD[len(rest)]} only "
    f"{pct(min(solved[k] for k in rest))} to {pct(max(solved[k] for k in rest))}",
    f"{tt(first_cov)} leads on MC and {tt(first_fac)} on FAC; FAC cannot tell the {WORD[2]} apart, MC can",
    f"in up to {pct1(max(digit))} of a model's responses, and of the {flags_decided} flags a domain expert decided on, {flags_slip} are real "
    f"slips (precision {f3(flag_precision)}); the judged step check flags a further step in up to {pct1(max(router))}",
    f"On their wrong answers, the {WORD[len(many_wrong)]} lower-tier models still state {pct(min(many_e3w))} to {pct(max(many_e3w))} of the gold "
    f"milestones, against a chance floor of at most {pct(max(many_floor))}, and up to {pct(max(full_cov))} of these wrong answers are complete "
    "derivations to a wrong value",
    f"A milestone the judge rules missing marks most of the lower tier's wrong answers ({pct(min(many_missing))} to {pct(max(many_missing))}) "
    f"but {pct(claude_attr['e5_missing'])} of {tt('claude-sonnet-5')}'s",
    f"on the {N_SINGLE} templates whose instances all follow one derivation, it solves some but not all instances of up to "
    f"{pct(max(low_some))} of the templates, the top tier of at most {pct(max(top_some))}",
    f"On {q5_pairs} pairs of an instance and its expert-checked paraphrase over {q5_templates} templates, no change in FAC holds after "
    f"correction, and for {WORD[len(q5_within)]} of the {WORD[len(ORDER)]} models the change is bounded within $\\pm {margin:.2f}$ at 90\\% "
    f"confidence ({tt(q5_out[0])} reaches {sgn(Q5[q5_out[0]]['ci90'][0])})",
    f"The domain experts rejected {pct(r_reject_share)} of the paraphrases that passed every scripted check",
    f"decoding {WORD[len(REPEATS)]} models three more times moves FAC by a standard deviation of at most {f3(max(rep_sd))}",
    # 5.3.2
    f"Each top-tier model's {WORD[len(BRANCH)]} branch means lie within {f2(min(top_span))} to {f2(max(top_span))} of each other and none of "
    f"its {len(domain_templates)} domain means falls below {top_dom_floor}, whereas each lower-tier model's branch means spread by "
    f"{f2(min(low_span))} to {f2(max(low_span))} and each falls below {low_dom_ceil} in at least one domain",
    f"with {n_branch_templates} templates per branch, differences below {round(min(detect_branch) * 100)} to {round(max(detect_branch) * 100)} "
    f"points cannot be told from sampling noise, and only one of the {n_branch_pairs} within-model branch comparisons holds after correction "
    f"({tt('gpt-oss-20b')}, electrical above civil)",
    f"Thermodynamics is the lowest domain for {WORD[len(thermo_models)]} of the {WORD[len(ORDER)]} models, largely because it holds the "
    f"{WORD[len(TWO_CHEMICAL)]} Advanced chemical templates whose wording does not pin the answer",
    f"without them it is the lowest for {WORD[n_thermo_wo]}",
    # 5.3.3
    f"Every model scores lower on Advanced than on Easy templates (\\autoref{{fig:level_bars}}), by {rng(gaps)}",
    f"The gap holds after correction for multiple comparisons only for the {WORD[len(gap_sig)]} lowest-scoring models, "
    f"{WORD[sum(k in NO_REASONING for k in gap_sig)]} of which run without reasoning tokens, and for no model once the {WORD[len(TWO_CHEMICAL)]} "
    f"Advanced chemical templates whose wording does not pin the answer are set aside, as {WORD[3]} chemical domain experts read them",
    f"for all {WORD[len(ORDER)]}, the share of responses scored 0 is higher on instances with {WORD[6]} or more gold milestones than on instances "
    f"with one, reaching {rng(depth6_low4)} for the {WORD[len(gap_sig)]} lowest-scoring models against at most {f3(max(depth6_top5))} for the "
    f"top {WORD[len(TOP5)]}",
    # 5.3.4
    f"{WORD[3].capitalize()} experiments on a {n_sub}-instance subset each change one thing in the evaluation",
    f"{tt('gpt-5.4-mini')} gains {f3(r_mini['diff'])} in FAC and {f3(mini_mc)} in MC while its arithmetic flags on correct answers roughly "
    f"halve, and the flagship {tt('gpt-5.4')}, which also returns no reasoning tokens at its default, gains {f3(flag_arm['diff'])} and sheds "
    "most of its flags",
    f"only the weakest model answers more correctly ({sgn(ob['gpt-oss-20b']['diff'])}, partly because fewer of its responses run out of room), "
    f"and the {WORD[len(tool)]} closed models stay within $\\pm {margin:.2f}$",
    f"which {tt('claude-sonnet-5')} calls on two thirds of the instances, neither closed model changes in FAC or MC beyond $\\pm {margin:.2f}$",
    f"cover the {WORD[len(tool)]} closed models in the tool experiment",
    # 5.4
    f"Three domain experts read the full response, the derivation and its final answer, of {ERR_BUILD['items']} wrong answers, "
    f"{per_model_items} from each of {WORD[len(B2_MODELS)]} models that contrast the tiers ({readings_total} readings), and assigned the first "
    f"category that applies in a {WORD[6]}-category hierarchy from hallucination to calculation, or ``no error'' (Fleiss'~$\\kappa$ "
    f"{f3(fleiss_all)}~\\citep{{fleiss1971}}",
    f"Calculation is the majority label for {tt('gpt-5.4-mini')} and {tt('gemma-4-26b-a4b')} ({maj_calc['gpt-5.4-mini']} of {per_model_items} "
    f"wrong answers each) and the largest for {tt('gpt-oss-20b')} ({maj_calc['gpt-oss-20b']} of {per_model_items}), the one model with many "
    f"conceptual errors (a hallucination, a wrong setup, or a wrong formula in {conceptual['gpt-oss-20b']} of {per_model_items}); a "
    f"hallucinated value or equation appears in {hall_readings} of the {readings_total} readings",
    f"On Easy problems {pct(easy_calc / easy_n)} of the readings are calculation errors, and wrong formulas or principles rise from "
    f"{pct(form_easy_share)} of the readings on Easy problems to {pct(form_hard_share)} on Intermediate and Advanced ones",
    f"For {tt('claude-sonnet-5')}, {claude_noerr} of {per_model_items} wrong answers are ``no error'' by majority, {claude_noerr_top3} of them "
    f"on {WORD[3]} templates whose question or check is at issue",
    f"Over {thousands(n_top5_responses)} responses, the top {WORD[len(TOP5)]} models lose {lost_points} answer points: {top5_empty} to "
    f"responses empty at the output ceiling, {top5_partial_points} to {top5_partial} partial answers at half credit, and {top5_incorrect} to "
    f"incorrect verdicts, of which {two_chemical} fall on the {WORD[len(TWO_CHEMICAL)]} chemical templates whose wording does not pin the "
    f"answer, {top5_symbolic} are symbolic answers scored by the numbers they state, and {top5_near} lie within 0.2\\% of a target whose "
    f"digits the question prescribes, leaving {remaining_incorrect} incorrect verdicts on other templates",
]

appendix_phrases = {
    "results": [
        "the smallest difference the design detects at 80\\% power",
        f"Of the {N_PAIRS} pairs, {len(SEP_FAC)} differ on FAC after Holm correction under the sign-flip test over templates: no pair inside "
        f"the top {WORD[len(TOP5)]}, {tt(sixth)} only from {tt(sixth_sep[0])}, each of the lower {WORD[len(rest)]} from every model of the "
        f"top {WORD[len(TOP5) + 1]}, and inside the lower {WORD[len(rest)]} only {tt(low_pair[0])} from {tt(low_pair[1])}",
        f"Every non-significant FAC difference lies below what its pair detects ({rng([pp['detectable'] for pp in nonsig])})",
        f"On MC, {n_mc_sig} pairs differ: the three highest models, {tt(top3_cov[0])}, {tt(top3_cov[1])}, and {tt(top3_cov[2])}, do not "
        f"separate, while {tt(first_cov)} and {tt(first_fac)}, which FAC does not separate, do; Kendall's $\\tau$~\\citep{{kendall1938}} "
        f"between the two orderings is {f3(tau['tau'])} (95\\% interval {ci(tau['ci'], False)})",
        f"Strict FAC, which scores a partial answer 0, separates {n_strict_sig} pairs, all but one of FAC's ({tt(strict_flip[0]['a'])} "
        f"against {tt(strict_flip[0]['b'])}); a Wilcoxon test~\\citep{{wilcoxon1945}} separates {n_wil_sig} on MC; and McNemar's exact "
        f"test~\\citep{{mcnemar1947}}, which treats the {thousands(N_ITEMS)} instances as independent, separates {n_mcnemar_sig}",
        f"swaps {WORD[n_tol_swaps]} pairs of models each, none of which differs after correction ($\\tau = {f3(sens_tau['half_tol'])}$ against "
        "the ordering as scored)",
        f"leaving out the {WORD[len(SHORTCUT)]} templates answerable from their wording or the {WORD[len(SYMBOLIC)]} with symbolic answers move "
        f"no model's FAC by more than {f3(other_shift)} and keep $\\tau$ at {f3(other_tau)} or above",
        f"($\\tau = {f3(sens_tau['unusable_excluded'])}$), chiefly because {tt('glm-5.3')}'s {Q1['glm-5.3']['empty']} empty responses then drop out",
        f"Spearman's $\\rho$~\\citep{{spearman1904}} between a template's mean coverage and its mean number of steps is {sgn(max(rho), 2)} to "
        f"{sgn(min(rho), 2)} for every model",
        f"under that permutation the gap holds for {WORD[gap_perm]} models",
        f"{tt('glm-5.3')}'s gap is {sgn(glm['unusable_excluded']['gap'])} (Holm-adjusted Welch $p$ {pv(glm['unusable_excluded']['p_welch_holm'])})",
    ],
    "branch_domain": [
        f"of {WORD[len(REPRESENTATIVE)]} representative models",
        f"The lowest branch differs by model, {lowest_list}, and a single branch pair differs within a model after correction",
        f"the {WORD[len(top_rep)]} top-tier models stay close to their FAC in every branch, while the {WORD[len(low_rep)]} lower-tier models are "
        f"lowest in different branches, {tt(low_rep[0])} in {BRANCH[lowest_branch[low_rep[0]]].lower()} and {tt(low_rep[1])} in "
        f"{BRANCH[lowest_branch[low_rep[1]]].lower()} engineering",
        f"the {WORD[len(top_rep)]} top-tier models are lowest on thermodynamics, which holds the {WORD[len(TWO_CHEMICAL)]} chemical templates "
        f"whose wording does not pin the answer, whereas the lower-tier models dip lower, {tt(low_rep[0])} in "
        f"{low_domain[low_rep[0]][0].replace('_', ' ')} and {tt(low_rep[1])} in {low_domain[low_rep[1]][0].replace('_', ' ')}",
    ],
    "paraphrase": [
        f"We select three of every template's 15 instances ({p_selected} instances)", f"A paraphrase must pass {WORD[2]} checks",
        f"a word similarity to the original of at most {COPY}", "with three attempts per instance",
        f"Of the {p_selected} instances, {p_passing} pass the scripted checks and the domain experts keep {r_kept} pairs on {q5_templates} templates",
        f"The {WORD[len(ORDER)]} models answer the kept paraphrases",
        f"the $\\pm {margin:.2f}$ margin, which equals the paired difference the design detects at 80\\% power, the 90\\% interval lies within it "
        f"for every model but {tt(q5_out[0])}",
        f"({' and '.join(tt(k) for k in below90)} below, {' and '.join(tt(k) for k in above90)} above)",
        f"For the {WORD[len(REPEATS)]} models decoded three more times on {rep_items} instances, FAC varies by a standard deviation of "
        f"{rng(rep_sd)} across repeats, and the paraphrase change lies within the spread of the repeats' changes for {WORD[len(vs_within)]} of "
        f"the {WORD[len(vs)]} and beyond it for {' and '.join(tt(k) for k in vs_beyond)}",
        f"On the {q5_pairs} kept pairs, no change in FAC holds after correction",
        f"Kendall's $\\tau$ {f3(q5_tau['tau'])} (95\\% interval {ci(q5_tau['ci'], False)}), at the lower edge of what sampling noise alone "
        f"gives, since two random halves of the same pairs agree with a median $\\tau$ of {f3(noise['median'])} (quartiles "
        f"{f3(noise['q1'])} to {f3(noise['q3'])})",
        f"it covers the {q5_templates} templates that can be paraphrased without loss",
    ],
    "experiments": [
        f"{WORD[4].capitalize()} experiments extend the evaluation",
        f"{n_sub}-instance subset (three instances per template)", f"the {WORD[2]} closed models",
        f"For the {ob_templates} templates whose code states their governing equations", f"({ob_items} instances)", f"up to {TOOL_MAX_CALLS} times",
        f"a {TOOL_TIMEOUT} s limit", f"truncated at {thousands(TOOL_OUTPUT_CHARS)} characters",
        f"the 90\\% interval against the $\\pm {margin:.2f}$ margin", "80\\% power",
        f"the change on the {ob['gpt-oss-20b']['usable_in_both']['items']} instances it answered in both runs is "
        f"{sgn(ob['gpt-oss-20b']['usable_in_both']['diff'])}",
        f"{tt('claude-sonnet-5')} on {pct(tool_share['claude-sonnet-5'])} of the instances and {tt('gpt-5.4-mini')} on "
        f"{pct(tool_share['gpt-5.4-mini'])}",
        f"{tt('deepseek-v4-pro')} scores {f3(fl['deepseek-v4-pro']['score'])} and {tt('gpt-5.4')} {f3(fl['gpt-5.4']['score'])} at its default "
        f"and {f3(flr['score'])} with reasoning, against {rng(sub_scores)} for the top {WORD[len(TOP5)]} on the same instances, whose intervals "
        f"contain both reasoning anchors; on Advanced templates the anchors score {rng(anchors_adv, f2)} and the top {WORD[len(TOP5)]} "
        f"{rng(sub_adv, f2)}",
    ],
    "errors": [
        f"The {WORD[6]} categories follow the stages", f"answers the {WORD[6]} questions in this order", f"{WORD[2]} further options",
        f"We read {ERR_BUILD['items']} wrong answers, {per_model_items} from each of {WORD[len(B2_MODELS)]} models",
        f"we draw the {per_model_items} at random", f"{ERR_BUILD['templates']} templates in all",
        f"{WORD[3].capitalize()} domain experts of the wrong answer's branch",
        f"{WORD[3].capitalize()} chemical domain experts read the {WORD[len(TWO_CHEMICAL)]} Advanced chemical templates",
        f"{WORD[3]} domain experts of each template's branch read {WORD[len(near_templates)]} templates on which many wrong answers land between "
        "0.2\\% and 5\\% of the gold value",
        f"all {WORD[3]} answer that the wording decides neither between the closed-system and the flow reading",
        f"nor between the {WORD[2]} forms of the truncated virial equation",
        f"{WORD[FLAME['data']['no: standard sources differ by more than that']]} answer that standard data sources differ by more than the "
        f"tolerance and one that they agree only within one source, {WORD[FLAME['method']['yes']]} that the gold method is the standard one and "
        f"one that it is not, and all {WORD[3]} that some of the responses shown are correct under some reading",
        f"The {WORD[len(near_templates)]} near-miss templates:}} for each, all {WORD[3]} answer that the question has one correct answer and that "
        "the response shown, within 5\\% of the gold value, is not correct",
        f"with and without the {WORD[len(TWO_CHEMICAL)]} chemical templates",
        f"{exact_templates} templates prescribe the digits of the answer, so a value within the tolerance but not at the prescribed digits "
        f"is incorrect by the question's own terms ({exact_incorrect} incorrect verdicts across the {WORD[len(ORDER)]} models, {near_total} of "
        f"them within 0.2\\% of the target on {len(near_rows)} of these templates)",
        f"unit and dimension errors are nearly absent ({unit_total} of {readings_total} readings) and sign errors rare ({sign_total}, in "
        f"{WORD[n_sign_models]} models), that {hall_gptoss} of the {hall_readings} hallucination readings are {tt('gpt-oss-20b')}'s, and that "
        f"agreement is high for every model and lowest for {tt('gpt-oss-20b')} (Fleiss' $\\kappa$ {f3(k_by[k_low_model])} against "
        f"{f3(min(k_others))} to {f3(max(k_others))})",
        f"the Advanced column's {pct(adv_noerr_share)} no-error share comes largely from {tt('claude-sonnet-5')}",
        f"Of the top {WORD[len(TOP5)]} models' {top5_incorrect} incorrect verdicts, {resid3_total} fall on {WORD[3]} templates: "
        f"{resid3['template_work_isothermal_virial']} on the compression-work template, {resid3['template_ber_estimation_mary']} on one symbolic "
        f"template ({code('template_ber_estimation_mary')}), and {resid3['template_adiabatic_flame_temperature']} on the flame-temperature template",
        f"the {WORD[len(SYMBOLIC)]} templates with symbolic answers",
        f"{WORD[len(SHORTCUT)].capitalize()} templates whose answers can be read off the question's wording stay in the evaluation set, since "
        f"leaving them out moves FAC by at most {f3(short_shift)}",
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


def bold(text: str, hit: bool) -> str:
    return f"\\textbf{{{text}}}" if hit else text


def with_interval(x: float, c, hit: bool) -> str:
    """A value followed by its interval in small gray type, on one line."""
    return f"{bold(f3(x), hit)} {{\\scriptsize\\textcolor{{gray}}{{({ci(c)})}}}}"


best = {"fac": max(Q1[k]["score"] for k in ORDER), "mc": max(COV[k]["coverage"] for k in ORDER), "solved": max(solved.values()),
        **{lv: max(BL[k]["level"][lv]["mean"] for k in ORDER) for lv in LEVELS}}


def main_row(k: str) -> str:
    q, d = Q1[k], Q3[k]
    return (f"{mark(k)} & {with_interval(q['score'], q['ci'], q['score'] == best['fac'])} & "
            f"{bold(f3(solved[k]), solved[k] == best['solved'])} & "
            f"{with_interval(COV[k]['coverage'], COV[k]['ci'], COV[k]['coverage'] == best['mc'])} & "
            f"{with_interval(d['digit_flag_rate_on_fully_solved'], d['digit_ci'], False)} & {f1(d['claims_per_trace'])} & "
            f"{with_interval(d['router_judge_rate_on_fully_solved'], d['router_judge_ci'], False)}")


def group_row(name: str, n: int) -> str:
    return f"\\rowcolor{{gray!10}}\\multicolumn{{{n}}}{{l}}{{\\textit{{{name}}}}} \\\\"


UP = " $\\uparrow$"
rows = [group_row("Open-weights LLMs", 7)] + [main_row(k) for k in OPEN] + [group_row("Closed LLMs", 7)] + [main_row(k) for k in CLOSED]
MAIN_HEAD = (" & \\multicolumn{2}{c}{\\textbf{Final answer}} & \\multicolumn{4}{c}{\\textbf{Derivation}} \\\\\n"
             "\\cmidrule(lr){2-3}\\cmidrule(lr){4-7}\n"
             "\\textbf{Model} & \\textbf{FAC" + UP + "} & " + mk("All 15 instances", "solved" + UP) + " & \\textbf{MC" + UP + "} & "
             + mk("Arithmetic flags,", "correct answers") + " & " + mk("Calculations", "read") + " & " + mk("Judged step flags,", "correct answers"))
blocks["main"]["tab:main_results"] = table(
    "l c c c c c c", MAIN_HEAD, rows,
    "\\textbf{Overall Performance of Evaluated LLMs on \\ourdataset.} "
    "Final Answer Accuracy and Milestone Coverage, with 95\\% intervals over templates in gray; the share of templates solved on every "
    "instance; and the step-level diagnostics on correct answers, reported as flags, not error counts, with the calculations the arithmetic "
    "check reads per response. Bold marks the best FAC, all-instances share, and MC. $^{\\ast}$No reasoning tokens at the provider's defaults. "
    f"$^{{\\dagger}}${Q1['glm-5.3']['empty']} responses empty at the output ceiling, scored 0.",
    "tab:main_results", shade_header=False, size="\\small\n\\setlength{\\tabcolsep}{4pt}", resize=True)
REP_NAMES = (f"{tt(top_rep[0])} (first on FAC), {tt(top_rep[1])} (first on MC), {tt(low_rep[0])}, and {tt(low_rep[1])}")
blocks["main"]["fig:level_bars"] = figure(
    "level-bars.pdf",
    "\\textbf{Final Answer Accuracy by Difficulty Level for Four Representative Models.} "
    f"The mean over the {n_easy} Easy, {n_int} Intermediate, and {n_adv} Advanced templates, with 95\\% intervals, for {REP_NAMES}; color marks "
    "the model and hatching the level. \\autoref{tab:level_gap} gives every model's Easy-minus-Advanced gap with its tests.",
    "fig:level_bars")
blocks["branch_domain"]["fig:domain_radar"] = figure(
    "domain_radar/domain-radar-labeled.pdf",
    "\\textbf{Final Answer Accuracy by Domain for Four Representative Models.} "
    f"The instance mean over each of the {len(domain_templates)} domains, in the order of their branches around the rim, for the same four "
    "models; the radial axis starts at 0.4. The domain means carry no interval and no test (\\autoref{tab:by_domain_kind}).",
    "fig:domain_radar", star=True, width="0.74\\textwidth")
_bp = branch_pairs[0][1]
assert len(BRANCH) == 5 and n_branch_templates * len(BRANCH) == N_TEMPLATES  # "a fifth" in the caption
blocks["branch_domain"]["fig:branch_bars"] = figure(
    "branch-bars.pdf",
    "\\textbf{Final Answer Accuracy by Engineering Branch for Four Representative Models.} "
    f"One stacked bar per model, for {REP_NAMES}. Each of the {WORD[len(BRANCH)]} branches holds {n_branch_templates} of the {N_TEMPLATES} "
    "templates, so a segment is that branch's share of the model's FAC: the number inside it is the branch's FAC, and its height is a fifth of "
    "that. The segments add up to the model's FAC, printed above the bar; color marks the model and hatching the branch. The one branch pair "
    f"that differs after Holm correction within a model is {tt(branch_pairs[0][0])}'s {BRANCH[_bp['b']].lower()} above "
    f"{BRANCH[_bp['a']].lower()} (\\autoref{{tab:branch_domain}}).",
    "fig:branch_bars")
blocks["main"]["fig:error_categories"] = figure(
    "error-categories.pdf",
    "\\textbf{Error Categories of the Wrong Answers Read.} "
    f"Each bar is one model's {readings_total // len(B2_MODELS)} readings (three domain experts, {per_model_items} wrong answers) as the share in "
    "each category, from the most fundamental at the top to ``no error'' at the bottom; shares of at least 10\\% are printed. "
    "\\autoref{tab:error_by_level} gives the readings by level.",
    "fig:error_categories")

# Appendix tables: one per topic.

# Branch and domain: shaded branch rows, each followed by its domains; the models across, in order of FAC.
dom_by_branch = {b: sorted(d for d in domain_templates if branch_of[d] == b) for b in BRANCH}
lower_held = {(k, pp["a"] if pp["diff"] < 0 else pp["b"]) for k, pp in branch_pairs}  # the lower branch of the pair that holds
rows = []
for b, bname in BRANCH.items():
    rows.append(f"\\rowcolor{{gray!10}}\\textit{{{bname}}} ({n_branch_templates}) & "
                + " & ".join(f3(BL[k]["branch"][b]["mean"]) + ("$^{\\ddagger}$" if (k, b) in lower_held else "") for k in ORDER) + " \\\\")
    rows += [f"\\quad {d.replace('_', ' ').capitalize()} ({domain_templates[d]}) & "
             + " & ".join(bold(f3(REP[k]["domain"][d]), lowest_domain[k] == d) for k in ORDER) for d in dom_by_branch[b]]
rows += ["\\midrule", "Detectable branch difference & " + " & ".join(f3(BL[k]["detectable_branch"]) for k in ORDER)]
low_b, high_b = (_bp["a"], _bp["b"]) if _bp["diff"] < 0 else (_bp["b"], _bp["a"])
blocks["branch_domain"]["tab:branch_domain"] = table(
    "l " + "r " * len(ORDER), "\\textbf{Branch or domain (templates)} & " + " & ".join("\\rotatebox{90}{" + tt(k) + "}" for k in ORDER), rows,
    "\\textbf{Final Answer Accuracy by Branch and by Domain.} "
    f"Models in order of FAC. Shaded rows: the mean of the branch's {n_branch_templates} template means; below each, its domains' instance "
    "means, which carry no interval or test; template counts in parentheses. Bold: each model's lowest domain. $^{\\ddagger}$The one branch "
    f"pair that differs within a model after Holm correction (Welch's $t$-test over the model's ten pairs): {tt(branch_pairs[0][0])}'s "
    f"{BRANCH[high_b].lower()} above its {BRANCH[low_b].lower()}. Last row: the smallest branch difference {n_branch_templates} templates "
    "detect at 80\\% power, per model.",
    "tab:branch_domain", size="\\footnotesize", resize=True, aliases=("tab:by_branch", "tab:by_domain_kind"))

# The level means, the level gap with its tests, and the share scored 0 by the depth of the gold derivation.
bins = [b for b in Q3[ORDER[0]]["by_milestone_count"] if b != "0"]
rows = []
for k in ORDER:
    q = Q2[k]
    rows.append(f"{tt(k)} & " + " & ".join(f3(BL[k]["level"][lv]["mean"]) for lv in LEVELS)
                + f" & {sgn(q['gap'])} ({ci(q['ci'])}) & {pv(q['p_welch_holm'])} & {f3(q['detectable_planned'])} & "
                f"{sgn(q['without_two_chemical']['gap'])} ({pv(q['without_two_chemical']['p_welch_holm'])}) & "
                + " & ".join(f3(Q3[k]["by_milestone_count"][b]["wrong_rate"]) for b in bins))
header = (" & \\multicolumn{3}{c}{\\textbf{FAC by level}} & \\multicolumn{4}{c}{\\textbf{Easy minus Advanced}} & \\multicolumn{"
          + str(len(bins)) + "}{c}{\\textbf{Share scored 0, by gold milestones}} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\\cmidrule(lr){9-"
          + str(8 + len(bins)) + "}\n"
          "\\textbf{Model} & \\textbf{Easy} & \\textbf{Intermediate} & \\textbf{Advanced} & \\textbf{Gap (95\\% interval)} & \\textbf{Welch $p$} & "
          "\\textbf{Detectable} & " + mk("Without two", "chemical") + " & " + " & ".join(f"\\textbf{{{b.replace('-', '--')}}}" for b in bins))
blocks["results"]["tab:level_gap"] = table(
    "l r r r c r r c " + "r " * len(bins), header, rows,
    "\\textbf{Difficulty Level and Trace Depth.} "
    f"Left: FAC over the {n_easy} Easy, {n_int} Intermediate, and {n_adv} Advanced templates. Middle: the Easy mean minus the Advanced mean, "
    "with its 95\\% interval, the Holm-adjusted $p$ of Welch's $t$-test, the smallest gap the design detects at 80\\% power, and the gap with "
    "its $p$ without the two Advanced chemical templates whose wording does not pin the answer. Right: the share of instances scored 0, empty "
    "responses included, by the number of milestones in the gold trace ("
    + ", ".join(str(Q3[ORDER[0]]["by_milestone_count"][b]["items"]) for b in bins[:-1])
    + f", and {Q3[ORDER[0]]['by_milestone_count'][bins[-1]]['items']} instances).",
    "tab:level_gap", size="\\footnotesize", resize=True)

# Coverage behind the wrong answers, the judge's share, and the precision of the arithmetic flags.
flag_by_model = {r[0]: (r[6], int(r[2]) - int(r[5])) for r in md_table(flags, r"\| model \| flags drawn") if r[0] != "all"}
assert set(flag_by_model) == set(ORDER)
rows = [f"{tt(k)} & {Q3[k]['readable_wrong_with_milestones']} & {f3(Q3[k]['e3_coverage_on_readable_wrong'])} & "
        f"{f3(Q3[k]['e3_null_on_readable_wrong'])} & {f3(Q3[k]['e5_coverage_on_readable_wrong'])} & {f3(COV[k]['wrong_full_coverage'])} & "
        f"{f3(Q3[k]['attribution_on_wrong']['e5_missing'])} & {f3(Q3[k]['e5_judged_fraction'])} & {flag_by_model[k][0]} ({flag_by_model[k][1]})"
        for k in ORDER]
header = (" & \\multicolumn{6}{c}{\\textbf{Wrong answers}} & \\textbf{All responses} & \\textbf{Correct answers} \\\\\n"
          "\\cmidrule(lr){2-7}\\cmidrule(lr){8-8}\\cmidrule(lr){9-9}\n"
          "\\textbf{Model} & \\textbf{Number} & " + mk("MC by", "matching") + " & " + mk("Chance", "floor") + " & " + mk("MC with", "judge")
          + " & " + mk("All", "reached") + " & " + mk("Milestone", "missing") + " & " + mk("Judge-decided", "milestones") + " & "
          + mk("Flag", "precision ($n$)"))
blocks["results"]["tab:coverage"] = table(
    "l r r r r r r r c", header, rows,
    "\\textbf{Milestone Coverage Behind Wrong Answers.} "
    "Per model, the readable wrong answers with milestones: how many; MC by matching alone and with the judge; the chance floor, the same "
    "responses matched against a sibling instance's milestones; the share that reach every milestone; and, over the answered wrong answers, "
    "the share with a milestone the judge rules missing. Then the share of all milestones the judge decides, and the precision of the "
    f"arithmetic flags a domain expert read, with their number ({flags_slip} of {flags_decided} were slips overall). The "
    f"{res['milestones']['items_without']} instances without milestones are left out.",
    "tab:coverage", size="\\footnotesize", resize=True)

# Paraphrase.
rows = [f"{tt(k)} & {sgn(Q5[k]['diff'])} ({ci(Q5[k]['ci'])}) & {ci(Q5[k]['ci90'])} & {pv(Q5[k]['p_holm'])} & {sgn(Q5[k]['e5']['diff'])}"
        for k in ORDER]
blocks["paraphrase"]["tab:paraphrase"] = table(
    "l c c r r", head("Model", "FAC change (95\\% interval)", "90\\% interval", "$p$ (Holm)", "MC change"), rows,
    "\\textbf{Change Under Paraphrase.} "
    f"Paraphrase minus original on the {q5_pairs} expert-kept pairs over {q5_templates} templates, positive when the paraphrase scores higher: "
    f"the change in FAC with its 95\\% interval, its 90\\% interval read against the $\\pm {margin:.2f}$ margin, the Holm-adjusted $p$ of the "
    "sign-flip test over templates, and the change in MC.",
    "tab:paraphrase", size="\\footnotesize")

# The four experiments, grouped by experiment.
ARM_GROUP = {"reasoning-medium": "Reasoning at medium effort", "openbook2": "Open book: governing equations supplied",
             "tool": "Open tool: Python interpreter offered",
             "flagship-reasoning-medium": f"{tt('gpt-5.4')} with reasoning at medium effort, against its default"}
rows, last = [], None
for a in ARMS:
    if a["arm"] != last:
        rows.append(group_row(ARM_GROUP[a["arm"]], 10))
        last = a["arm"]
    rows.append(f"{tt(a['model'])} & {a['items']} & {f3(a['main_score_on_items'])} & {f3(a['arm_score'])} & {sgn(a['diff'])} ({ci(a['ci'])}) & "
                f"{pv(a['p_holm'])} & {f3(a['detectable'])} & {ci(a['ci90'])} & {sgn(a['e5']['diff'])} & "
                f"{f3(a['digit_flag_rate_fully_solved']['main'])} / {f3(a['digit_flag_rate_fully_solved']['arm'])}")
blocks["experiments"]["tab:experiments"] = table(
    "l r r r c r r c r c",
    "\\textbf{Model} & \\textbf{Instances} & " + mk("FAC,", "evaluation") + " & " + mk("FAC,", "experiment") + " & "
    + mk("Change", "(95\\% interval)") + " & \\textbf{$p$ (Holm)} & \\textbf{Detectable} & \\textbf{90\\% interval} & " + mk("MC", "change")
    + " & " + mk("Arithmetic flags,", "evaluation / experiment"), rows,
    "\\textbf{The Four Experiments Against the Evaluation.} "
    "Per model and experiment: the instances covered; FAC in the evaluation and under the experiment on those instances, paired by instance; "
    "the change with its 95\\% interval, the Holm-adjusted $p$ of the sign-flip test over templates, the smallest change the experiment "
    f"detects at 80\\% power, and the 90\\% interval read against the $\\pm {margin:.2f}$ margin; the change in MC; and the share of "
    "correct-answer responses with an arithmetic flag in the evaluation and under the experiment (unpaired).",
    "tab:experiments", size="\\footnotesize", resize=True)

# Error analysis: the readings by model and by level in one table.
cats = [(full, short) for full, short in CATEGORIES if full != INCOMPLETE]
rows = [f"{short} & " + " & ".join(f"{by_model[m].get(full, 0)} ({maj[m].get(full, 0)})" for m in B2_TEXT) + " & "
        + " & ".join(f"{by_level[lv].get(full, 0)} ({pct(by_level[lv].get(full, 0) / level_n[lv])})" for lv in LEVELS) for full, short in cats]
drawn_level = {lv: sum(ERR_BUILD["by_model_level"][m][lv] for m in B2_MODELS) for lv in LEVELS}
rows += ["\\midrule",
         "Wrong answers read by level & " + " & ".join(", ".join(str(ERR_BUILD["by_model_level"][m][lv]) for lv in LEVELS) for m in B2_TEXT)
         + " & " + " & ".join(str(drawn_level[lv]) for lv in LEVELS),
         "Fleiss' $\\kappa$ & " + " & ".join(f3(ERR["fleiss_by_model"][m]) for m in B2_TEXT) + " & & & "]
header = (" & \\multicolumn{" + str(len(B2_TEXT)) + "}{c}{\\textbf{By model: readings (majority labels)}} & \\multicolumn{3}{c}{\\textbf{By "
          "level: readings (share)}} \\\\\n\\cmidrule(lr){2-" + str(1 + len(B2_TEXT)) + "}\\cmidrule(lr){" + str(2 + len(B2_TEXT)) + "-"
          + str(4 + len(B2_TEXT)) + "}\n\\textbf{Category} & " + " & ".join(tt(m) for m in B2_TEXT)
          + " & \\textbf{Easy} & \\textbf{Intermediate} & \\textbf{Advanced}")
blocks["errors"]["tab:errors"] = table(
    "l " + "r " * (len(B2_TEXT) + 3), header, rows,
    "\\textbf{The Readings of the Wrong Answers.} "
    f"Three domain experts read each of {ERR_BUILD['items']} wrong answers, {per_model_items} per model. Left: readings per category "
    f"({readings_total // len(B2_MODELS)} per model) and, in parentheses, the wrong answers whose majority label is the category; then the "
    "wrong answers read per level. Right: readings per category by level and their share of the level's readings. Categories run from the "
    f"most to the least fundamental; no domain expert used the incomplete option. Fleiss' $\\kappa$ over the three domain experts is "
    f"{f3(fleiss_all)} overall.",
    "tab:errors", size="\\footnotesize\n\\setlength{\\tabcolsep}{4pt}", resize=True, aliases=("tab:error_by_level",))


# ----------------------------------------------------------------------------------------------- figures
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif", "font.size": 6.5, "axes.linewidth": 0.4,
                         "xtick.major.width": 0.4, "ytick.major.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 0,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.labelsize": 6.5, "legend.fontsize": 6,
                         "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": INK, "axes.labelcolor": INK,
                         "hatch.linewidth": 0.4, "savefig.dpi": 300})
    return plt


def vector_hatches(fig) -> None:
    """Redraw each hatched patch's hatch as clipped vector strokes. Matplotlib writes hatches as PDF tiling patterns, which PDF
    viewers rasterize at low resolution, so they look blurred on screen; the strokes keep the pattern's own geometry, colour
    and line width (72 pt cells anchored at the page's top left), so the figure looks the same, only sharp."""
    import numpy as np
    import matplotlib as mpl
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path as MPath
    fig.canvas.draw()  # fixes every extent, the legends' included
    to_inches, dy = fig.dpi_scale_trans.inverted(), fig.get_figheight() % 1.0
    legends = list(fig.legends) + [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
    owner = {id(h): lg for lg in legends for h in lg.get_patches()}
    for patch in fig.findobj(lambda a: isinstance(a, Patch) and bool(a.get_hatch())):
        (x0, y0), (x1, y1) = to_inches.transform(patch.get_window_extent().get_points())
        cell = MPath.hatch(patch.get_hatch())
        tiles = [cell.vertices + (i, j + dy) for i in range(int(np.floor(x0)) - 1, int(np.ceil(x1)) + 1)
                 for j in range(int(np.floor(y0 - dy)) - 1, int(np.ceil(y1 - dy)) + 1)]
        color = patch.get_hatchcolor()
        width = patch.get_hatch_linewidth() if hasattr(patch, "get_hatch_linewidth") else mpl.rcParams["hatch.linewidth"]
        legend = owner.get(id(patch))
        strokes = PathPatch(MPath(np.concatenate(tiles), np.concatenate([cell.codes] * len(tiles))), transform=fig.dpi_scale_trans,
                            facecolor=color, edgecolor=color, linewidth=width,
                            zorder=(legend.get_zorder() if legend is not None else patch.get_zorder()) + 0.01)
        patch.set_hatch(None)
        parent = legend.axes if legend is not None and legend.axes is not None else (fig if legend is not None else patch.axes or fig)
        parent.add_artist(strokes)
        strokes.set_clip_path(patch)  # after adding: an axes clips a new artist to itself when it has no clip path of its own


def save(fig, name: str) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    vector_hatches(fig)
    fig.savefig(FIGS / name, metadata={"CreationDate": None})
    import matplotlib.pyplot as plt
    plt.close(fig)


def recess(ax, left: bool = True) -> None:
    for s in ("top", "right") + (() if left else ("left",)):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e7e6e2", linewidth=0.4)
    ax.set_axisbelow(True)


def rows_axes(plt, n: int, height: float, width: float = COLUMN, left: float = 0.36, right: float = 0.97, bottom: float = 0.13,
              top: float = 0.90, ncols: int = 1, wspace: float = 0.08):
    fig, axes = plt.subplots(1, ncols, figsize=(width, height), sharey=True,
                             gridspec_kw={"left": left, "right": right, "bottom": bottom, "top": top, "wspace": wspace})
    axes = list(axes) if ncols > 1 else [axes]
    for ax in axes:
        ax.set_ylim(n - 0.5, -0.5)
        recess(ax, left=ax is axes[0])
    return fig, axes


def interval_rows(ax, items, lw: float = 1.0) -> None:
    """items: (y, lo, hi, x, filled, marker) rows: an interval line with a marker, filled or hollow."""
    for y, lo, hi, x, filled, marker in items:
        ax.plot([lo, hi], [y, y], color=BLUE, linewidth=lw, solid_capstyle="butt", zorder=2)
        ax.plot([x], [y], marker=marker, markersize=4.2, markeredgewidth=0.8, markeredgecolor=BLUE,
                markerfacecolor=BLUE if filled else "white", linestyle="none", zorder=3)


def gradient_rows(ax, items, height: float = 0.38, steps: int = 80) -> None:
    """items: (y, lo, hi, x, color) rows: an interval band in color, deepest at the estimate x and paler toward its ends. Each
    slice runs on under the next, which is drawn over it, so no seam shows between slices."""
    from matplotlib.colors import to_rgb
    for y, lo, hi, x, color in items:
        rgb, span, step = to_rgb(color), max(x - lo, hi - x), (hi - lo) / steps
        for j in range(steps):
            t = abs(lo + (j + 0.5) * step - x) / span  # 0 at the estimate, 1 at the farther end
            ax.barh(y, step * (2 if j < steps - 1 else 1), left=lo + j * step, height=height, color=[c + (1 - c) * 0.72 * t for c in rgb],
                    linewidth=0, zorder=2)


def top_legend(fig, handles, ncol: int, y: float = 0.995) -> None:
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, frameon=False, handlelength=1.4,
               columnspacing=1.0, handletextpad=0.5)


def fig_level_gap() -> None:
    plt = _plt()
    from matplotlib.lines import Line2D
    # The interval as a band, deepest at the gap; the gaps that hold in strong blue with a dark filled marker, the rest paler and hollow.
    hold, rest, rest_edge = BLUE, "#93b6e2", "#5b8fd3"
    holds = {k: Q2[k]["p_welch_holm"] < 0.05 for k in ORDER}
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(ORDER), 2.6, left=0.318, right=0.97, bottom=0.12, top=0.985)
        ax.set_yticks(range(len(ORDER)))
        ax.set_yticklabels([FIG_NAME[k] for k in ORDER], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        ax.axvline(0, color=MUTED, linewidth=0.6, zorder=1)
        ax.text(-0.004, -0.42, "No gap", rotation=90, ha="right", va="top", fontsize=5.5, style="italic", color=MUTED)
        gradient_rows(ax, [(i, Q2[k]["ci"][0], Q2[k]["ci"][1], Q2[k]["gap"], hold if holds[k] else rest) for i, k in enumerate(ORDER)])
        for i, k in enumerate(ORDER):
            ax.plot([Q2[k]["gap"]], [i], marker="o", markersize=4.6, linestyle="none", zorder=3,
                    **({"markerfacecolor": DARK, "markeredgecolor": "white", "markeredgewidth": 0.6} if holds[k] else
                       {"markerfacecolor": "white", "markeredgecolor": rest_edge, "markeredgewidth": 0.9}))
        ax.plot([Q2[k]["without_two_chemical"]["gap"] for k in ORDER], range(len(ORDER)), marker="x", markersize=3.8, markeredgewidth=0.8,
                color=INK, linestyle="none", zorder=4)
        ax.set_xlim(-0.03, 0.36)
        ax.set_xlabel("Final Answer Accuracy drop from Easy to Advanced", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="Holds after Holm"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=rest_edge, markeredgewidth=0.9, markersize=4.6,
                          linestyle="none", label="Does not hold"),
                   Line2D([], [], marker="x", color=INK, markersize=3.8, markeredgewidth=0.8, linestyle="none", label="Without chemical pair"),
                   Line2D([], [], color="#7fb0ea", linewidth=4.5, solid_capstyle="butt", label="95% interval")]
        legend = ax.legend(handles=handles, loc="upper right", frameon=True, fancybox=False, framealpha=1, edgecolor="black", fontsize=6,
                           borderaxespad=0.3, borderpad=0.35, handlelength=1.3, handletextpad=0.4, labelspacing=0.3)
        legend.get_frame().set_linewidth(0.5)
        save(fig, "level-gap.pdf")






def fig_coverage_wrong() -> None:
    plt = _plt()
    from matplotlib.lines import Line2D
    # A band from the chance floor, palest, to the coverage with the judge, deepest; matching alone sits on it.
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(ORDER), 2.8, left=0.318, right=0.91, bottom=0.111, top=0.877)
        ax.set_yticks(range(len(ORDER)))
        ax.set_yticklabels([FIG_NAME[k] for k in ORDER], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        rows = [(i, Q3[k]["e3_null_on_readable_wrong"], Q3[k]["e3_coverage_on_readable_wrong"], Q3[k]["e5_coverage_on_readable_wrong"])
                for i, k in enumerate(ORDER)]
        gradient_rows(ax, [(i, f, e5, e5, BLUE) for i, f, _, e5 in rows])
        for i, f, e3, e5 in rows:
            ax.plot([f], [i], marker="s", markersize=3.8, color=MUTED, markeredgecolor="white", markeredgewidth=0.5, linestyle="none", zorder=3)
            ax.plot([e3], [i], marker="o", markersize=4.6, markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=0.9, linestyle="none",
                    zorder=3)
            ax.plot([e5], [i], marker="o", markersize=4.6, markerfacecolor=DARK, markeredgecolor="white", markeredgewidth=0.6, linestyle="none",
                    zorder=4)
        for i, k in enumerate(ORDER):
            ax.text(1.025, i, str(Q3[k]["readable_wrong_with_milestones"]), va="center", ha="left", fontsize=6.5, color="black", clip_on=False)
        ax.text(1.025, -0.8, "n", ha="left", va="center", fontsize=6.5, style="italic", color="black", clip_on=False)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("Milestone Coverage of wrong answers", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="s", color=MUTED, markersize=3.8, linestyle="none", label="Chance floor"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=0.9, markersize=4.6, linestyle="none",
                          label="Matching alone"),
                   Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="With the judge")]
        legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=((0.318 + 0.91) / 2, 0.995), ncol=3, frameon=True,
                            fancybox=False, framealpha=1, edgecolor="black", fontsize=6, borderpad=0.35, handlelength=1.0,
                            handletextpad=0.4, columnspacing=1.2)
        legend.get_frame().set_linewidth(0.5)
        save(fig, "coverage-wrong.pdf")




def band_axes(plt, n: int, height: float, xlabel: str, left: float = 0.36, top: float = 0.88):
    fig, (ax,) = rows_axes(plt, n, height, left=left, top=top)
    ax.axvspan(-margin, margin, color=BAND, zorder=0)
    ax.axvline(0, color=MUTED, linewidth=0.5, linestyle=(0, (2, 2)), zorder=1)
    ax.set_xlabel(xlabel)
    return fig, ax


def fig_paraphrase() -> None:
    plt = _plt()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    # The 90% interval as a band, deepest at the change; within the margin in strong blue with a dark filled marker, the rest paler and hollow.
    rest, rest_edge = "#93b6e2", "#5b8fd3"
    within = {k: Q5[k]["within_margin"] for k in ORDER}
    with plt.rc_context(SERIF):
        fig, (ax,) = rows_axes(plt, len(ORDER), 2.6, left=0.318, right=0.72, bottom=0.12, top=0.985)
        ax.set_yticks(range(len(ORDER)))
        ax.set_yticklabels([FIG_NAME[k] for k in ORDER], fontweight="bold")
        ax.tick_params(axis="y", colors="black", labelsize=6.5)
        ax.tick_params(axis="x", colors="black", labelsize=6.5)
        ax.axvspan(-margin, margin, color=BAND, zorder=0)
        ax.axvline(0, color=MUTED, linewidth=0.6, zorder=1)
        ax.text(0.0015, 2.5, "No change", rotation=90, ha="left", va="center", fontsize=5.5, style="italic", color=MUTED)
        gradient_rows(ax, [(i, Q5[k]["ci90"][0], Q5[k]["ci90"][1], Q5[k]["diff"], BLUE if within[k] else rest) for i, k in enumerate(ORDER)])
        for i, k in enumerate(ORDER):
            ax.plot([Q5[k]["diff"]], [i], marker="o", markersize=4.6, linestyle="none", zorder=3,
                    **({"markerfacecolor": DARK, "markeredgecolor": "white", "markeredgewidth": 0.6} if within[k] else
                       {"markerfacecolor": "white", "markeredgecolor": rest_edge, "markeredgewidth": 0.9}))
        ax.set_xlim(-0.072, 0.058)
        ax.set_xticks([-0.05, 0, 0.05])
        ax.set_xticklabels(["−0.05", "0", "0.05"])
        ax.set_xlabel("Final Answer Accuracy change under paraphrase", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        handles = [Line2D([], [], marker="o", color=DARK, markersize=4.6, linestyle="none", label="Within the margin"),
                   Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=rest_edge, markeredgewidth=0.9, markersize=4.6,
                          linestyle="none", label="Not within"),
                   Line2D([], [], color="#7fb0ea", linewidth=4.5, solid_capstyle="butt", label="90% interval"),
                   Patch(color=BAND, label=f"\u00b1{margin:.2f} margin")]
        legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=True, fancybox=False, framealpha=1,
                           edgecolor="black", fontsize=6, borderaxespad=0, borderpad=0.35, handlelength=1.1, handletextpad=0.4,
                           labelspacing=0.3)
        legend.get_frame().set_linewidth(0.5)
        save(fig, "paraphrase.pdf")




def fig_error_categories() -> None:
    plt = _plt()
    from matplotlib.patches import Patch
    cats = [c for c in CATEGORIES if c[0] != INCOMPLETE]  # no "Incomplete" reading was given
    # The May figure's palette: the errors before calculation in browns, darker the more fundamental, calculation in hatched
    # mustard, no error in green; the lightness order and the hatching carry the categories in greyscale.
    fills = ["#5c3b24", "#7f5638", "#a87c5b", "#b78f71", "#c3a084", "#d9b86c", "#b3dda0"]
    hatches = ["...", "...", "...", "...", "...", "//", ""]
    labels = {"claude-sonnet-5": "Claude\nSonnet 5", "gpt-5.4-mini": "GPT-5.4\nmini", "gemma-4-26b-a4b": "Gemma 4\n26B", "gpt-oss-20b": "GPT OSS\n20B"}
    with plt.rc_context({"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 7, "legend.fontsize": 7}):
        fig, ax = plt.subplots(figsize=(COLUMN, 2.3), gridspec_kw={"left": 0.11, "right": 0.75, "bottom": 0.15, "top": 0.97})
        for x, m in enumerate(B2_MODELS):
            n, bottom = sum(by_model[m].values()), 0.0
            for (full, _), fill, hatch in reversed(list(zip(cats, fills, hatches))):  # no error at the base, the most fundamental on top
                v = by_model[m].get(full, 0) / n
                ax.bar(x, v, bottom=bottom, width=0.6, color=fill, linewidth=0, zorder=2)
                if hatch:
                    ax.bar(x, v, bottom=bottom, width=0.6, fill=False, hatch=hatch, edgecolor=HATCH, linewidth=0, zorder=3)
                ax.bar(x, v, bottom=bottom, width=0.6, fill=False, edgecolor="white", linewidth=0.6, zorder=4)
                if v >= 0.1:
                    ax.text(x, bottom + v / 2, f"{v * 100:.0f}", ha="center", va="center", fontsize=7, fontweight="bold", zorder=5,
                            color="white" if fill in fills[:2] else "black", bbox={"facecolor": fill, "edgecolor": "none", "pad": 0.6})
                bottom += v
        ax.set_xticks(range(len(B2_MODELS)))
        ax.set_xticklabels([labels[m] for m in B2_MODELS], fontweight="bold")
        ax.tick_params(axis="x", length=0, colors="black", labelsize=7)
        ax.set_xlim(-0.42, len(B2_MODELS) - 0.58)
        ax.tick_params(axis="y", colors="black", labelsize=7)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0", "25", "50", "75", "100"])
        ax.set_ylabel("Share of the readings (%)", fontsize=6.5, fontweight="bold", color="black", labelpad=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        handles = [Patch(facecolor=fill, hatch=hatch, edgecolor=HATCH, linewidth=0, label=short.split(" or ")[0])
                   for (_, short), fill, hatch in zip(cats, fills, hatches)]
        legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.025, 1.0), frameon=True, fancybox=False, framealpha=1,
                           edgecolor="black", fontsize=6, borderaxespad=0, borderpad=0.35, handlelength=1.1, handleheight=0.9,
                           handletextpad=0.4, labelspacing=0.25)
        legend.get_frame().set_linewidth(0.5)
        vector_hatches(fig)  # the hatching as vector lines, so it stays sharp in PDF viewers; nothing else changes
        save(fig, "error-categories.pdf")


MODEL_COLORS = {"deepseek-v4.1-flash": "#2a78d6", "claude-sonnet-5": "#1baf7a", "gpt-5.4-mini": "#4a3aa7", "gpt-oss-20b": "#eb6834"}
HATCHES = ["", "//", "..", "xx", "--"]


FIG_NAME_2 = {"deepseek-v4.1-flash": "DeepSeek\nV4.1 Flash", "claude-sonnet-5": "Claude\nSonnet 5", "gpt-5.4-mini": "GPT-5.4\nmini",
              "gpt-oss-20b": "GPT OSS\n20B"}  # the figure labels on two lines, for the column width
BAR_LEFT, BAR_RIGHT = 0.10, 0.995


def vector_hatches(fig) -> None:
    """Redraw each hatched patch's hatch as clipped vector lines. Matplotlib writes hatches as PDF tiling patterns, which PDF
    viewers rasterize at low resolution; the lines keep the pattern's geometry, colour and width (72 pt cells anchored at the
    page's top left), so the figure looks the same, only sharp."""
    import numpy as np
    import matplotlib as mpl
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path as MPath
    fig.canvas.draw()  # fixes every extent, the legends' included
    to_inches, dy = fig.dpi_scale_trans.inverted(), fig.get_figheight() % 1.0
    legends = list(fig.legends) + [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
    owner = {id(h): lg for lg in legends for h in lg.get_patches()}
    for patch in fig.findobj(lambda a: isinstance(a, Patch) and bool(a.get_hatch())):
        (x0, y0), (x1, y1) = to_inches.transform(patch.get_window_extent().get_points())
        cell = MPath.hatch(patch.get_hatch())
        tiles = [cell.vertices + (i, j + dy) for i in range(int(np.floor(x0)) - 1, int(np.ceil(x1)) + 1)
                 for j in range(int(np.floor(y0 - dy)) - 1, int(np.ceil(y1 - dy)) + 1)]
        color = patch.get_hatchcolor()
        width = patch.get_hatch_linewidth() if hasattr(patch, "get_hatch_linewidth") else mpl.rcParams["hatch.linewidth"]
        legend = owner.get(id(patch))
        lines = PathPatch(MPath(np.concatenate(tiles), np.concatenate([cell.codes] * len(tiles))), transform=fig.dpi_scale_trans,
                          facecolor=color, edgecolor=color, linewidth=width,
                          zorder=(legend.get_zorder() if legend is not None else patch.get_zorder()) + 0.01)
        patch.set_hatch(None)
        parent = legend.axes if legend is not None and legend.axes is not None else (fig if legend is not None else patch.axes or fig)
        parent.add_artist(lines)
        lines.set_clip_path(patch)  # after adding: an axes clips a new artist to itself when it has no clip path of its own


def column_bar_axes(plt):
    return plt.subplots(figsize=(COLUMN, 2.45), gridspec_kw={"left": BAR_LEFT, "right": BAR_RIGHT, "bottom": 0.155, "top": 0.80})


def finish_column_bars(fig, ax, labels: list[str], hatches: list[str], legend_title: str) -> None:
    """The shared frame of the two bar figures: model names under the bars, the FAC axis, the boxed legend of the hatches."""
    from matplotlib.patches import Patch
    ax.set_xticks(range(len(REPRESENTATIVE)))
    ax.set_xticklabels([FIG_NAME_2[m] for m in REPRESENTATIVE], fontsize=7, fontweight="bold", linespacing=1.0)
    ax.set_xlim(-0.55, len(REPRESENTATIVE) - 0.45)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("Final Answer Accuracy", fontsize=6.5, fontweight="bold", labelpad=2)
    ax.grid(axis="y", linestyle="--", color="#d9d8d4", linewidth=0.4)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=6, pad=1.5)
    handles = [Patch(facecolor="#8c8c8c" if not h else "#e6e6e6", edgecolor=INK, hatch=h, linewidth=0.5, label=l) for l, h in zip(labels, hatches)]
    legend = fig.legend(handles=handles, title=legend_title, loc="upper center", bbox_to_anchor=((BAR_LEFT + BAR_RIGHT) / 2, 1.0),
                        ncol=len(labels), frameon=True, edgecolor=INK, fancybox=False, fontsize=6, title_fontsize=6.5, handlelength=1.6,
                        handleheight=0.9, handletextpad=0.4, columnspacing=0.9, borderpad=0.35)
    legend.get_title().set_fontweight("bold")
    legend.get_frame().set_linewidth(0.5)
    vector_hatches(fig)


def grouped_bars(name: str, groups: list[tuple[str, str]], legend_title: str, values) -> None:
    """The May layout at the ACL column width: one color per representative model, one hatch per group, the value above every
    bar, a boxed legend; values(model, group label) returns (mean, (lo, hi))."""
    plt = _plt()
    import numpy as np
    from matplotlib.colors import to_rgb
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "hatch.linewidth": 0.4})
    fig, ax = column_bar_axes(plt)
    x = np.arange(len(REPRESENTATIVE))
    width = 0.8 / len(groups)
    for i, model in enumerate(REPRESENTATIVE):
        color = MODEL_COLORS[model]
        tint = tuple(1 - 0.3 * (1 - c) for c in to_rgb(color))
        for j, (label, hatch) in enumerate(groups):
            mean, (lo, hi) = values(model, label)
            pos = x[i] - 0.4 + width * (j + 0.5)
            ax.bar(pos, mean, width=width * 0.92, facecolor=color if not hatch else tint, edgecolor=color, hatch=hatch, linewidth=0.5, zorder=2)
            ax.errorbar(pos, mean, yerr=[[mean - lo], [hi - mean]], fmt="none", ecolor=INK, elinewidth=0.45, capsize=1.0, capthick=0.45, zorder=3)
            ax.text(pos, hi + 0.012, f"{mean:.2f}", ha="center", va="bottom", fontsize=5.5, fontweight="bold", color=INK)
    finish_column_bars(fig, ax, [g[0] for g in groups], [g[1] for g in groups], legend_title)
    save(fig, name)


def fig_level_bars() -> None:
    grouped_bars("level-bars.pdf", list(zip(LEVELS, HATCHES)), "Difficulty Level",
                 lambda m, lv: (BL[m]["level"][lv]["mean"], BL[m]["level"][lv]["ci"]))


def text_on(fill) -> str:
    """White or black, whichever contrasts more with the fill (WCAG relative luminance)."""
    from matplotlib.colors import to_rgb
    r, g, b = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in to_rgb(fill)]
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if 1.05 / (y + 0.05) > (y + 0.05) / 0.05 else "black"


def fig_branch_bars() -> None:
    """FAC by branch for the four representative models at the ACL column width, as stacked bars. The five branches hold 30
    templates each, so a segment is one branch's share of the model's FAC (the branch mean divided by five) and a bar's height
    is the model's FAC, printed above it; the label in a segment gives the branch mean. The colors and hatches are those of the
    level figure."""
    plt = _plt()
    from matplotlib.colors import to_rgb
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "hatch.linewidth": 0.4})
    fig, ax = column_bar_axes(plt)
    width, n = 0.56, len(BRANCH)
    for i, model in enumerate(REPRESENTATIVE):
        color = MODEL_COLORS[model]
        tint = tuple(1 - 0.3 * (1 - c) for c in to_rgb(color))
        bottom = 0.0
        for (key, _), hatch in zip(BRANCH.items(), HATCHES):
            assert BL[model]["branch"][key]["templates"] * n == N_TEMPLATES
            mean = BL[model]["branch"][key]["mean"]
            fill = color if not hatch else tint
            ax.bar(i, mean / n, bottom=bottom, width=width, facecolor=fill, edgecolor=color, hatch=hatch, linewidth=0.5, zorder=2)
            if bottom > 0:  # a white seam between segments
                ax.plot([i - width / 2, i + width / 2], [bottom, bottom], color="white", linewidth=0.8, solid_capstyle="butt", zorder=3)
            ax.text(i, bottom + mean / n / 2, f"{mean:.2f}", ha="center", va="center", fontsize=6, fontweight="bold", zorder=4,
                    color=text_on(fill), bbox={"facecolor": fill, "edgecolor": "none", "pad": 0.5})
            bottom += mean / n
        assert abs(bottom - Q1[model]["score"]) < 1e-9  # the stack is the model's FAC
        ax.text(i, bottom + 0.012, f"{bottom:.2f}", ha="center", va="bottom", fontsize=6, fontweight="bold", color=INK)
    finish_column_bars(fig, ax, list(BRANCH.values()), HATCHES[:n], "Engineering Branch")
    save(fig, "branch-bars.pdf")


DOMAIN_LABEL = {  # the paper's domain names, broken for the radar's rim
    "reaction_kinetics": "Reaction\nKinetics", "thermodynamics": "Thermodynamics", "transport_phenomena": "Transport\nPhenomena",
    "digital_communications": "Digital\nCommunications", "electromagnetics_and_waves": "Electromagnetics\nand Waves",
    "signals_and_systems": "Signals and\nSystems", "fluid_mechanics": "Fluid\nMechanics", "mechanics_of_materials": "Mechanics of\nMaterials",
    "vibrations_and_acoustics": "Vibrations and\nAcoustics", "geotechnical_engineering": "Geotechnical\nEngineering",
    "structural_analysis": "Structural\nAnalysis", "water_resources": "Water\nResources", "production_and_inventory": "Production and\nInventory",
    "quality_and_reliability_control": "Quality and\nReliability Control", "stochastic_operations": "Stochastic\nOperations",
}
RADAR_BRANCH_ORDER = ["chemical_engineering", "electrical_engineering", "mechanical_engineering", "civil_engineering", "industrial_engineering"]


def fig_domain_radar() -> None:
    """FAC by domain for the four representative models, domains grouped by branch around the rim, in two versions:
    with the domain names (domain-radar-labeled.pdf) and without (domain-radar-unlabeled.pdf), under figs/domain_radar/."""
    import csv
    import numpy as np
    from matplotlib.lines import Line2D
    plt = _plt()
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"]})
    with (HERE / "results/per_template.csv").open(encoding="utf-8") as f:
        branch_of = {r["domain"]: r["branch"] for r in csv.DictReader(f)}
    domains = [d for b in RADAR_BRANCH_ORDER for d in sorted((d for d, br in branch_of.items() if br == b), key=lambda d: DOMAIN_LABEL[d])]
    assert len(domains) == len(DOMAIN_LABEL) == 15
    angles = np.linspace(0, 2 * np.pi, len(domains), endpoint=False)
    closed = np.append(angles, angles[0])
    styles = [((0, (1, 1.2)), "o"), ((0, (4, 1, 1, 1, 1, 1)), "^"), ((0, (4, 1.5, 1, 1.5)), "s"), ("-", "D")]
    (FIGS / "domain_radar").mkdir(parents=True, exist_ok=True)
    for labeled in (True, False):
        fig = plt.figure(figsize=(4.8, 4.2))
        ax = fig.add_axes([0.20, 0.13, 0.60, 0.65], polar=True)
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        handles = []
        for m, (ls, mk) in zip(REPRESENTATIVE, styles):
            vals = [REP[m]["domain"][d] for d in domains]
            vals = np.append(vals, vals[0])
            color = MODEL_COLORS[m]
            ax.fill(closed, vals, color=color, alpha=0.05, zorder=1)
            ax.plot(closed, vals, linestyle=ls, linewidth=1.1, color=color, marker=mk, markersize=3.6, markerfacecolor="white",
                    markeredgewidth=0.9, zorder=3)
            handles.append(Line2D([], [], linestyle=ls, linewidth=1.1, color=color, marker=mk, markersize=3.6, markerfacecolor="white",
                                  markeredgewidth=0.9, label=FIG_NAME.get(m, NAME[m])))
        ax.set_xticks(angles)
        ax.set_xticklabels([DOMAIN_LABEL[d] for d in domains] if labeled else [], fontsize=6.5, fontweight="bold", color=INK)
        ax.tick_params(axis="x", pad=-1)
        for label, theta in zip(ax.get_xticklabels(), angles):  # align each rim label away from the circle
            x, y = np.cos(np.pi / 2 - theta), np.sin(np.pi / 2 - theta)
            label.set_ha("left" if x > 0.15 else "right" if x < -0.15 else "center")
            label.set_va("bottom" if y > 0.15 else "top" if y < -0.15 else "center")
        ax.set_ylim(0.4, 1.07)  # the rim sits outside the 1.0 ring, so the polygons never touch it
        ax.set_yticks([0.6, 0.8, 1.0])
        ax.set_yticklabels(["0.6", "0.8", "1.0"], fontsize=5.5, color=MUTED)
        ax.set_rlabel_position(90)
        for label in ax.get_yticklabels():
            label.set_bbox({"facecolor": "white", "edgecolor": "none", "pad": 0.4})
        ax.grid(color="#cfcecb", linewidth=0.5)
        ax.spines["polar"].set_color(INK)
        ax.spines["polar"].set_linewidth(0.8)
        legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=4, frameon=True, fancybox=False,
                            edgecolor=INK, fontsize=7, handlelength=2.2, columnspacing=1.0)
        legend.get_frame().set_linewidth(0.6)
        save(fig, f"domain_radar/domain-radar-{'labeled' if labeled else 'unlabeled'}.pdf")


FIGURES = {  # the four placed figures; fig_domain_radar also draws the unlabeled radar beside the labeled one
    "branch-bars.pdf": fig_branch_bars, "error-categories.pdf": fig_error_categories, "level-bars.pdf": fig_level_bars,
    "domain_radar/domain-radar-labeled.pdf": fig_domain_radar}


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


def write_blocks(draw: bool = True) -> None:
    for key, path in FILES.items():
        tex = path.read_text(encoding="utf-8")
        for name, body in blocks[key].items():
            pattern = re.compile(r"% BEGIN GENERATED " + re.escape(name) + r" .*?% END GENERATED " + re.escape(name), re.S)
            if not pattern.search(tex):
                raise SystemExit(f"{path.name} has no markers for {name}")
            tex = pattern.sub(lambda m: block(name, body), tex, count=1)
        path.write_text(tex, encoding="utf-8")
    if draw:
        for fn in FIGURES.values():
            fn()


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
        refs = set(re.findall(r"\\autoref\{([^}]*)\}", tex))
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
    for name in FIGURES:
        if not (FIGS / name).exists():
            print(f"MISSING figure {FIGS / name}")
            failures += 1
    return failures


if __name__ == "__main__":
    if "--write" in sys.argv:
        text_only = "--text-only" in sys.argv
        write_blocks(draw=not text_only)
        print(f"wrote {sum(len(b) for b in blocks.values())} generated blocks"
              + ("" if text_only else f" and {len(FIGURES)} figures to {FIGS.relative_to(REPO)}"))
    elif "--check" in sys.argv:
        sys.exit(1 if check() else 0)
    else:
        for key, ph in (("6_results.tex", phrases), *((f"appendices/{k}.tex", v) for k, v in appendix_phrases.items())):
            print(f"% phrases {key} must contain")
            print("\n".join(ph))
