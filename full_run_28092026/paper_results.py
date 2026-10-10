"""Generate the tables, the figures and the numbers of the paper's Results and Error Analysis (Sections 5.3 and 5.4) and
their appendices, and check them in the source.

    python full_run_28092026/paper_results.py                       # print the phrases the prose must contain
    python full_run_28092026/paper_results.py --write [flags]       # rewrite every generated block in the .tex files and draw the figures
    python full_run_28092026/paper_results.py --out DIR [flags]     # the same into a copy of the tex tree under DIR (Phase 1: the real
                                                                    # tree is never written); a block whose markers no file holds goes to
                                                                    # DIR/generated/<label>.tex, and the phrases to DIR/generated/phrases.txt
    python full_run_28092026/paper_results.py --check [--out DIR]   # exit 1 unless every placed block is current, every phrase is in its
                                                                    # file, every number in the prose is a phrase's, every citation key
                                                                    # resolves, every \\autoref label is defined, every figure exists and
                                                                    # the sources agree; with --out, blocks and phrases the prose does not
                                                                    # hold yet are listed as pending and do not fail the check
    python full_run_28092026/paper_results.py --stand-in            # write stand-in result files (marked "stand_in": true) under
                                                                    # results/stand_in/ for the files WS-C1 to C3 have not written yet

Flags (the registry in docs/mock_review_workstreams/00_ORCHESTRATION.md, section 7):
    --headline default|matched   default: the models at the providers' default settings in Table 1, with a marked block "with
                                 reasoning at medium effort" for the re-run models; matched: each model at its matched setting,
                                 the re-run models marked, a twelfth model if present, and every other block, phrase, claim and
                                 figure read from the matched family (the expert readings and the experiments stay at the
                                 providers' defaults); tab:matched always holds both configurations
    --repaired                   the two chemical templates are repaired: the "with and without" clauses and the B4-reading phrases
                                 are not generated, and tab:level_gap drops its "without two chemical" column
    --judged-in-table            the judged step flags stay in Table 1 (they are in tab:coverage otherwise)
    --text-only                  with --write or --out: the blocks only, the figures as drawn

Files written: 6_results.tex (Table 1, the two main-text figures and the prose phrases), appendices/results.tex,
appendices/branch_domain.tex, appendices/paraphrase.tex, appendices/further_experiments.tex, appendices/error_analysis.tex (their
tables, between "% BEGIN GENERATED <name>" and "% END GENERATED <name>" markers; a block is written wherever in the tree its markers
are), and the four placed figures under figs/ (drawn by paper_figures.py).

Sources, each written by a committed script: results/results.json (analyze.py: every score, interval, test and condition);
results/matched.json, single_path.json, providers.json (analyze.py, WS-C1), coverage_variants.json (coverage_variants.py),
sensitivity_variants.json, flag_precision.json, depth_model.json (clause_variants.py, flag_precision.py, depth_model.py); until one
of these exists, its stand-in under results/stand_in/ is read and every block and phrase that reads it says STAND-IN;
results/matched_config.json (run_traces.py --matched-config: which models reason by default and which were re-run);
results/judge_swap_main.json (judge_swap.py: the second judge's MC shift per model);
expert_request/scored_current.json (expert_kits.py --score: the experts' readings the store still holds, counts only);
RESIDUAL_INCORRECT.md (RESIDUAL_INCORRECT_matched.md under --headline matched), PARAPHRASE.md, PARAPHRASE_REVIEW.md, FLAG_REVIEW_3.md; the shortcut template list of analyze.py; the tool and
paraphrase constants of run_traces.py and paraphrase.py; the symbolic-equivalence template list of the answer module; the judged
step check's precision row of appendices/validation.tex (docs/appendix_evaluation.py); subsamples.py (the repeat and paraphrase
subsamples). Model display names are those of paper_setup.py.

Checks: a structural assertion stops the script (a result file in an unexpected shape). A claim the prose makes that the data no
longer supports is recorded by claim() and printed as "CLAIM FAILS"; a disagreement between two sources is recorded by drift() and
printed as "SOURCE DRIFT"; both fail --check and never stop generation, so the writers see the numbers and the failing claims
together.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import shutil
import sys
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from full_run_28092026 import paper_figures  # noqa: E402  # the drawing code (WS-F): draw_all(results, figs_dir, out_dir)
from full_run_28092026 import subsamples  # noqa: E402  # the repeat and paraphrase subsamples

SRC = REPO / "overleaf_source_04102026"
APPX = SRC / "appendices"
FIGS = SRC / "figs"
RESULTS = HERE / "results"
STAND_IN_DIR = RESULTS / "stand_in"
FILES = {"main": "6_results.tex", "results": "appendices/results.tex", "branch_domain": "appendices/branch_domain.tex",
         "paraphrase": "appendices/paraphrase.tex", "experiments": "appendices/further_experiments.tex",
         "errors": "appendices/error_analysis.tex", "limitations": "limitations.tex",
         "conclusion": "conclusion.tex"}  # the files that carry phrases, relative to the tree
BIB = "custom.bib"
MARK = "% BEGIN GENERATED {name} (full_run_28092026/paper_results.py --write)\n{body}\n% END GENERATED {name}"
GENERATED_FILES = ("matched", "single_path", "providers", "coverage_variants", "sensitivity_variants", "flag_precision", "depth_model")


def parse_args(argv: list[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(add_help=True, description="the paper's generated blocks, figures and phrases")
    ap.add_argument("--write", action="store_true", help="write the blocks (and figures) into the tex tree")
    ap.add_argument("--check", action="store_true", help="check the tree; exit 1 on any failure")
    ap.add_argument("--out", metavar="DIR", help="work on a copy of the tex tree under DIR (implies --write unless --check)")
    ap.add_argument("--headline", choices=("default", "matched"), default="default")
    ap.add_argument("--repaired", action="store_true")
    ap.add_argument("--judged-in-table", action="store_true")
    ap.add_argument("--text-only", action="store_true")
    ap.add_argument("--stand-in", action="store_true", help="write stand-in result files under results/stand_in/ and exit")
    args, _unknown = ap.parse_known_args(argv)
    return args


ARGS = parse_args(sys.argv[1:] if __name__ == "__main__" else [])
HEADLINE, REPAIRED, JUDGED_IN_TABLE = ARGS.headline, ARGS.repaired, ARGS.judged_in_table
CLAIMS: list[str] = []  # prose claims the data no longer supports
DRIFT: list[str] = []  # disagreements between sources


def claim(ok: bool, text: str) -> None:
    if not ok:
        CLAIMS.append(text)


def drift(ok: bool, text: str) -> None:
    if not ok:
        DRIFT.append(text)


NAME = {
    "deepseek-v4.1-flash": "DeepSeek V4.1 Flash", "kimi-k3": "Kimi K3", "claude-sonnet-5": "Claude Sonnet 5",
    "glm-5.3-flash": "GLM-5.3-Flash", "muse-glimmer-30b": "Muse Glimmer 30B", "glm-5.3": "GLM-5.3",
    "qwen3-235b-a22b-2507": "Qwen3-235B-2507", "gemini-3.1-flash-lite": "Gemini 3.1 Flash-Lite",
    "gemma-4-26b-a4b": "Gemma 4 26B", "gpt-5.4-mini": "GPT-5.4 mini", "gpt-oss-20b": "gpt-oss-20b",
    "gpt-5.4": "GPT-5.4", "deepseek-v4-pro": "DeepSeek V4 Pro", "qwen3-235b-a22b-thinking-2507": "Qwen3-235B-Thinking-2507",
}
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
WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
        11: "eleven", 12: "twelve", 13: "thirteen"}
ORDINAL = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"]
NUMBER = re.compile(r"(?<![\w.\-])\d[\d,]*(?:\.\d+)?")
NUMBER_WORDS = re.compile(r"\b(two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)\b", re.I)


# ----------------------------------------------------------------------------------------------- formatting
def f3(x: float) -> str:
    return f"{x:.3f}"


def f2(x: float) -> str:
    return f"{x:.2f}"


def f1(x: float) -> str:
    return f"{x:.1f}"


def ci(c, dash: bool = True, nd: int = 3) -> str:
    """An interval: 'lo--hi' in tables when both ends are positive, otherwise 'lo to hi' with a math-mode minus."""
    f = lambda x: f"{x:.{nd}f}"  # noqa: E731
    if dash and c[0] >= 0 and c[1] >= 0:
        return f"{f(c[0])}--{f(c[1])}"
    return " to ".join(f"$-{f(abs(x))}$" if round(x, nd) < 0 else f(abs(x)) if round(x, nd) == 0 else f(x) for x in c)  # never "-0.000"


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
    return "\\texttt{" + NAME.get(key, key) + "}"


def code(name: str) -> str:
    return "\\texttt{" + name.replace("template_", "").replace("_", "\\_") + "}"


def thousands(n: int) -> str:
    return f"{n:,}"


def rng(xs, fmt=f3, join=" to ") -> str:
    return f"{fmt(min(xs))}{join}{fmt(max(xs))}"


def series(items: list[str]) -> str:
    """'a', 'a and b', 'a, b, and c'."""
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    if len(items) == 2:
        return f"{items[0]} and {items[1]}"
    return ", ".join(items[:-1]) + f", and {items[-1]}"


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


def group_row(name: str, n: int) -> str:
    return f"\\rowcolor{{gray!10}}\\multicolumn{{{n}}}{{l}}{{\\textit{{{name}}}}} \\\\"


def with_interval(x: float, c, hit: bool = False) -> str:
    """A value followed by its interval in small gray type, on one line."""
    s = f3(x)
    if hit:
        s = f"\\textbf{{{s}}}"
    return f"{s} {{\\scriptsize\\textcolor{{gray}}{{({ci(c)})}}}}"


def bold(text: str, hit: bool) -> str:
    return f"\\textbf{{{text}}}" if hit else text


def kendall_tau(order_a: list[str], order_b: list[str]) -> float:
    """Kendall's tau-a between two orderings of the same keys."""
    ra, rb = {k: i for i, k in enumerate(order_a)}, {k: i for i, k in enumerate(order_b)}
    keys = [k for k in order_a if k in rb]
    s = sum(((ra[a] - ra[b]) * (rb[a] - rb[b]) > 0) - ((ra[a] - ra[b]) * (rb[a] - rb[b]) < 0) for a, b in itertools.combinations(keys, 2))
    n = len(keys)
    return s / (n * (n - 1) / 2) if n > 1 else 1.0


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


STAND_IN: set[str] = set()  # the result files read from results/stand_in/
QUICK: set[str] = set()  # the result files written in quick mode


def result_file(name: str) -> dict:
    """results/<name>.json, or its stand-in under results/stand_in/ (recorded in STAND_IN) when the real file does not exist."""
    real, alt = RESULTS / f"{name}.json", STAND_IN_DIR / f"{name}.json"
    if real.exists():
        d = load(real)
        if isinstance(d, dict) and d.get("stand_in"):
            raise SystemExit(f"{real} is marked stand_in: a stand-in must live under results/stand_in/, never under results/")
    elif alt.exists():
        d = load(alt)
        STAND_IN.add(name)
    else:
        raise SystemExit(f"results/{name}.json is missing (WS-C1 to C3 write it); run with --stand-in to write a stand-in under "
                         f"results/stand_in/ for development")
    if (isinstance(d, dict) and d.get("quick")) or (isinstance(d, list) and any(isinstance(r, dict) and r.get("quick") for r in d)):
        QUICK.add(name)
    return d  # a dict, or a list for providers.json (the C1 schema)


def note(*names: str) -> str:
    """The caption prefix for a block that reads a stand-in or a quick-mode file: nothing of the kind reaches the paper unnoticed."""
    tags = (["STAND-IN values"] if any(n in STAND_IN for n in names) else []) + (["QUICK-mode values"] if any(n in QUICK for n in names) else [])
    return ("; ".join(tags) + f" (from {', '.join(f'{n}.json' for n in names if n in STAND_IN | QUICK)}). ") if tags else ""


res = load(RESULTS / "results.json")
RES_DEFAULT = res  # the providers' defaults: tab:matched's default columns, the experiments and the anchors
FAMILY = ("q1", "q2", "q3_overall", "q3_coverage", "branches_levels", "q4", "q5", "sensitivity", "reported", "repeats")
if HEADLINE == "matched" and not ARGS.stand_in:  # every per-model family at matched settings (matched.json mirrors results.json)
    _mf = result_file("matched")
    res = {**res, **{key: _mf[key] for key in FAMILY}, "q3": [_mf["models"][k]["q3"] for k in _mf["order"]]}
matched_config = load(RESULTS / "matched_config.json")
scored = load(HERE / "expert_request/scored_current.json")  # the readings the store still holds (amendment 3 of 2026-10-07)
RESIDUAL_MD = "RESIDUAL_INCORRECT_matched.md" if HEADLINE == "matched" and not ARGS.stand_in else "RESIDUAL_INCORRECT.md"  # residual_incorrect.py [--matched]
residual = (HERE / RESIDUAL_MD).read_text(encoding="utf-8")
paraphrase_md = (HERE / "PARAPHRASE.md").read_text(encoding="utf-8")
review = (HERE / "PARAPHRASE_REVIEW.md").read_text(encoding="utf-8")
flags = (HERE / "FLAG_REVIEW_3.md").read_text(encoding="utf-8")
validation_tex = (APPX / "validation.tex").read_text(encoding="utf-8")
analyze_src = (HERE / "analyze.py").read_text(encoding="utf-8")
traces_src = (HERE / "run_traces.py").read_text(encoding="utf-8")
paraphrase_src = (HERE / "paraphrase.py").read_text(encoding="utf-8")
answer_src_path = REPO / "evaluator_pilot_17092026/evaluators/answer.py"
answer_src = answer_src_path.read_text(encoding="utf-8") if answer_src_path.exists() else ""
TOOL_MAX_CALLS, TOOL_TIMEOUT, TOOL_OUTPUT_CHARS = (int(re.search(rf"^{c} = (\d+)", traces_src, re.M).group(1))
                                                  for c in ("TOOL_MAX_CALLS", "TOOL_TIMEOUT", "TOOL_OUTPUT_CHARS"))
COPY = float(re.search(r"^COPY = ([\d.]+)", paraphrase_src, re.M).group(1))  # the paraphrase's word-similarity ceiling

Q1 = {m["model"]: m for m in res["q1"]["models"]}
ORDER = sorted(Q1, key=lambda k: -Q1[k]["score"])  # the table order: Final Answer Accuracy, descending; every model the stores hold
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
N_MULTI = Q4[ORDER[0]]["multi_path"]["templates"]
assert all(sum(m[k] for k in ("correct", "partial", "incorrect", "unusable")) == N_ITEMS for m in Q1.values())
assert N_PAIRS == len(CPAIRS) == len(ORDER) * (len(ORDER) - 1) // 2, "one pairwise test per pair of models"
assert N_ITEMS == 15 * N_TEMPLATES and N_SINGLE + N_MULTI == N_TEMPLATES
assert set(Q2) == set(Q3) == set(Q3O) == set(COV) == set(BL) == set(Q4) == set(Q5) == set(SENS) == set(REP) == set(Q1)

# The matched-settings configuration (run_traces.py --matched-config): which models reason at their providers' defaults, which were
# re-run with reasoning at medium effort, and which entries are inert (run: false, or no store).
CFG = [m for m in matched_config["models"] if m.get("run", True) and m.get("default_store")]
CFG_BY = {m["model"]: m for m in CFG}
RERUN = [k for k in ORDER if CFG_BY.get(k, {}).get("reasoning_store")]  # re-run with reasoning at medium effort, in table order
NO_REASONING = {k for k in ORDER if CFG_BY.get(k, {}).get("reasoning_setting") != "reasons by default"}  # no reasoning tokens at default
REASONING_STORE = {k: CFG_BY[k]["reasoning_store"] for k in RERUN}
drift(set(CFG_BY) == set(ORDER), f"matched_config.json names {sorted(set(CFG_BY) ^ set(ORDER))} differently from results.json's models")

ERR = scored["kinds"]["error"]
ERR_COMP = ERR["composition"]  # the B2 sample as it stands: items, items read, readers, by model and level
TEMPLATES_READ = scored["kinds"]["template"]["templates"]
SHORTCUT = re.findall(r"'(template_[a-z_]+)'", re.search(r"SHORTCUT = \[(.*?)\]", analyze_src, re.S).group(1))
SYMBOLIC = res["symbolic_templates"]  # the templates whose answer is an expression
_sym_enabled = re.search(r"SYMBOLIC_EQUIVALENCE_TEMPLATES\s*=\s*[\[(](.*?)[\])]", answer_src, re.S)
SYMBOLIC_ENABLED = re.findall(r"['\"](template_[a-z_0-9]+)['\"]", _sym_enabled.group(1)) if _sym_enabled else []  # the equivalence check's list
BY_NUMBERS = [t for t in SYMBOLIC if t not in SYMBOLIC_ENABLED]  # still scored by the numbers they state
Q1_DEFAULT = {m["model"]: m for m in RES_DEFAULT["q1"]["models"]}  # the providers' defaults, whatever the headline
ORDER_DEFAULT = sorted(Q1_DEFAULT, key=lambda k: -Q1_DEFAULT[k]["score"])
B2_MODELS = sorted(ERR["by_model"], key=lambda m: ORDER_DEFAULT.index(m))  # the models whose wrong answers the experts read (default responses)


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
SEP_FAC_DEFAULT = separated(RES_DEFAULT["q1"]["pairs"])
CLD_FAC_DEFAULT = letters(ORDER_DEFAULT, SEP_FAC_DEFAULT)
TOP5 = ORDER[:5]
top5_sep_pairs = [tuple(sorted(pr, key=ORDER.index)) for pr in itertools.combinations(TOP5, 2) if frozenset(pr) in SEP_FAC]
nonsig = [p for p in PAIRS if p["p_holm"] >= 0.05]
claim(all(abs(p["diff"]) < p["detectable"] for p in nonsig), "every non-significant FAC difference lies below what its pair detects")
strict_agree = sum((p["p_holm"] < 0.05) == (p["fully_p_holm"] < 0.05) for p in PAIRS)
scores = [Q1[k]["score"] for k in ORDER]
spread5 = max(scores[:5]) - min(scores[:5])
sixth, rest = ORDER[5], ORDER[6:]
sixth_sep = [k for k in TOP5 if frozenset((sixth, k)) in SEP_FAC]
MATCHED_HEAD = HEADLINE == "matched"  # the tables and phrases show the matched family: a graded order, no two tiers
if MATCHED_HEAD:
    claim(not sixth_sep, "the sixth model differs from none of the first five")
else:
    claim(0 < len(sixth_sep) < len(TOP5), "the sixth model differs from some but not all of the top five")
    claim(all(frozenset((a, b)) in SEP_FAC for a in rest for b in ORDER[:6]), "each lower-tier model differs from every model above it")
cross_sep = sum(frozenset((a, b)) in SEP_FAC for a in ORDER[:6] for b in rest)  # pairs between the first six and the last five that differ
no_variance = [Q1[k]["templates_no_instance_variance"] for k in ORDER]

gaps = [Q2[k]["gap"] for k in ORDER]
gap_sig = [k for k in ORDER if Q2[k]["p_welch_holm"] < 0.05]
gap_sig_chem = [k for k in ORDER if Q2[k]["without_two_chemical"]["p_welch_holm"] < 0.05]
gap_perm = sum(Q2[k]["p_perm_holm"] < 0.05 for k in ORDER)
claim(0 < len(gap_sig) < gap_perm, "the level gap holds for some models under Welch's test and for more under the permutation test")
if not REPAIRED:
    claim(not gap_sig_chem, "no level gap holds without the two chemical templates")
top5_gaps = [Q2[k]["gap"] for k in TOP5]
top5_excl0 = sum(Q2[k]["ci"][0] > 0 for k in TOP5)
glm = Q2["glm-5.3"]
branch_pairs = [(k, p) for k in ORDER for p in BL[k]["pairs"] if p["p_holm"] < 0.05]
claim(len(branch_pairs) == 1 and branch_pairs[0][0] == "gpt-oss-20b" and branch_pairs[0][1]["a"] == "civil_engineering"
      and branch_pairs[0][1]["b"] == "electrical_engineering" and branch_pairs[0][1]["diff"] < 0,
      "one within-model branch pair holds: gpt-oss-20b, electrical above civil")
n_branch_pairs = sum(len(BL[k]["pairs"]) for k in ORDER)
lowest_branch = {k: min(BL[k]["branch"], key=lambda b: BL[k]["branch"][b]["mean"]) for k in ORDER}
claim(len(set(lowest_branch.values())) > 1, "the lowest branch differs by model")
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
if MATCHED_HEAD:
    claim(frozenset((first_fac, first_cov)) not in SEP_FAC and dc["p_holm"] >= 0.05, "neither FAC nor MC separates the first on FAC and the first on MC")
else:
    claim(frozenset((first_fac, first_cov)) not in SEP_FAC and dc["p_holm"] < 0.05, "FAC does not separate the first on FAC and the first on MC; MC does")
n_mc_sig, n_wil_sig = len(SEP_MC), len(separated(CPAIRS, "p_wilcoxon_holm"))
agree = sum((p["p_holm"] < 0.05) == (p["p_wilcoxon_holm"] < 0.05) for p in CPAIRS)
top3_cov = sorted(ORDER, key=lambda k: cov_rank[k])[:3]
claim(not any(frozenset(p) in SEP_MC for p in itertools.combinations(top3_cov, 2)), "the three highest MC models do not separate")
rho = [COV[k]["rho_steps"] for k in ORDER]

REPRESENTATIVE = ["deepseek-v4.1-flash", "claude-sonnet-5", "gpt-5.4-mini", "gpt-oss-20b"]  # the four models the figures show
top_rep, low_rep = REPRESENTATIVE[:2], REPRESENTATIVE[2:]
if MATCHED_HEAD:
    claim(top_rep == [first_fac, first_cov] and ORDER.index(low_rep[0]) == 4 and low_rep[1] == ORDER[-1],
          "the representative models are the first on FAC, the first on MC, the fifth and the last")
else:
    claim(top_rep == [first_fac, first_cov] and low_rep[0] in NO_REASONING and low_rep[0] in ORDER[6:] and low_rep[1] == ORDER[-1],
          "the representative models are the first on FAC, the first on MC, a no-reasoning lower-tier model and the last")
rep_branch = {k: [BL[k]["branch"][b]["mean"] for b in BRANCH] for k in REPRESENTATIVE}
top_branch = rep_branch[top_rep[0]] + rep_branch[top_rep[1]]
top_domain = [v for k in top_rep for v in REP[k]["domain"].values()]
low_domain = {k: (min(REP[k]["domain"], key=REP[k]["domain"].get), min(REP[k]["domain"].values()), max(REP[k]["domain"].values()))
              for k in low_rep}
with (RESULTS / "per_template.csv").open(encoding="utf-8") as f:
    PT = list(csv.DictReader(f))
domain_templates = Counter(r["domain"] for r in PT if r["model"] == ORDER[0])
n_easy, n_adv = BL[ORDER[0]]["level"]["Easy"]["templates"], BL[ORDER[0]]["level"]["Advanced"]["templates"]
n_int = BL[ORDER[0]]["level"]["Intermediate"]["templates"]

e3w = [Q3[k]["e3_coverage_on_readable_wrong"] for k in ORDER]
floor = [Q3[k]["e3_null_on_readable_wrong"] for k in ORDER]
e5w = [Q3[k]["e5_coverage_on_readable_wrong"] for k in ORDER]
e5w_max_model = max(ORDER, key=lambda k: Q3[k]["e5_coverage_on_readable_wrong"])
many_wrong = [k for k in ORDER if COV[k]["wrong_n"] > 100]
full_cov = [COV[k]["wrong_full_coverage"] for k in many_wrong]
attr = {c: [Q3[k]["attribution_on_wrong"][c] for k in ORDER] for c in ("digit_rule", "e5_missing", "router_judge")}

# The step checks' flag rates on correct answers: digit (arithmetic), router_judge (judged steps), and any further deterministic
# check whose rate analyze.py writes as <kind>_flag_rate_on_fully_solved (a formula check, amendment 5), detected from the rows.
FLAG_KINDS = ["digit"] + sorted(k[:-len("_flag_rate_on_fully_solved")] for k in Q3[ORDER[0]]
                                if k.endswith("_flag_rate_on_fully_solved") and k.split("_")[0] not in ("digit", "tol1", "router"))
FLAG_LABEL = {"digit": "Arithmetic flags", "formula": "Formula flags", "router_judge": "Judged step flags"}
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
# The judged step check alone, inside correct answers: its precision and the flags behind it, as validation.tex prints them.
router_alone = re.search(r"Judged step check alone\S* & P, all steps; correct answers & [\d.]+; ([\d.]+) \((\d+) of (\d+)\)", validation_tex)
assert router_alone, "the judged step check's own precision row is missing from appendices/validation.tex"
router_alone_precision, router_alone_tp, router_alone_n = float(router_alone.group(1)), int(router_alone.group(2)), int(router_alone.group(3))
assert abs(router_alone_tp / router_alone_n - router_alone_precision) < 5e-4

depth6 = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in ORDER]
depth1 = [Q3[k]["by_milestone_count"]["1"]["wrong_rate"] for k in ORDER]
claim(all(a > b for a, b in zip(depth6, depth1)), "every model is scored 0 more often on six-plus-milestone instances than on one-milestone ones")
single_all = {k: Q4[k]["single_path"]["all"] for k in ORDER}
single_some = {k: Q4[k]["single_path"]["some"] for k in ORDER}
strong6, weak5 = ORDER[:6], ORDER[6:]
weak4 = [k for k in weak5 if k != "gemini-3.1-flash-lite"]

q5 = [Q5[k] for k in ORDER]
q5_diff = [m["diff"] for m in q5]
q5_within = [k for k in ORDER if Q5[k]["within_margin"]]
q5_out = [k for k in ORDER if not Q5[k]["within_margin"]]
claim(len(q5_out) >= 1 and all(m["p_holm"] >= 0.05 for m in q5), "some paraphrase change is not within the margin and none holds after correction")
assert all(m["items"] == q5[0]["items"] for m in q5)
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
drift(p_selected == len(subsamples.paraphrase_ids()), f"PARAPHRASE.md selects {p_selected} items; subsamples.PARAPHRASE gives {len(subsamples.paraphrase_ids())}")
drift(p_passing == r_returned, f"PARAPHRASE.md: {p_passing} pass the checks; PARAPHRASE_REVIEW.md: {r_returned} returned to the experts")
drift(r_kept == q5_pairs, f"PARAPHRASE_REVIEW.md keeps {r_kept} pairs; results.json's Q5 runs on {q5_pairs}")
drift(r_kept + r_rejected == r_returned, "PARAPHRASE_REVIEW.md: kept plus rejected is not returned")
drift(p_lost == len(p_lost_templates) and p_selected - p_passing == p_failed, "PARAPHRASE.md's funnel does not add up")
p_lost_experts = N_TEMPLATES - q5_templates - p_lost
least_branch = min(r_branch, key=lambda b: r_branch[b][0])
below90 = [k for k in ORDER if Q5[k]["ci90"][1] < 0]
above90 = [k for k in ORDER if Q5[k]["ci90"][0] > 0]
vs = res["q5"]["vs_repeats"]
vs_within = [k for k in ORDER if k in vs and vs[k]["abs_within_repeat_spread"]]
vs_beyond = [k for k in ORDER if k in vs and not vs[k]["abs_within_repeat_spread"]]
rep_sd = [REPEATS[k]["sd"] for k in REPEATS]
rep_same = [REPEATS[k]["same_verdict_every_repeat"] for k in REPEATS]
rep_items = REPEATS[next(iter(REPEATS))]["items"]
assert all(REPEATS[k]["items"] == rep_items for k in REPEATS)
drift(rep_items == len(subsamples.repeat_ids()), f"results.json's repeats cover {rep_items} items; subsamples.REPEAT gives {len(subsamples.repeat_ids())}")
REPEAT_MODELS = sorted(REPEATS, key=ORDER.index)  # the models decoded three more times, in table order
sens_tau = res["sensitivity"]["tau_with_headline"]
short_shift = max(abs(SENS[k]["without_shortcut_templates"] - SENS[k]["fitted"]) for k in ORDER)
sym_shift = [SENS[k]["without_symbolic_templates"] - SENS[k]["fitted"] for k in ORDER]

arm = {(a["arm"], a["model"]): a for a in ARMS}
r_mini, r_gem = arm[("reasoning-medium", "gpt-5.4-mini")], arm[("reasoning-medium", "gemini-3.1-flash-lite")]
claim(r_gem["p_holm"] >= 0.05 and abs(r_gem["diff"]) < r_gem["detectable"] and r_mini["p_holm"] < 0.05,
      "reasoning at medium effort moves GPT-5.4 mini and not Gemini 3.1 Flash-Lite on the subset")
ob = {m: arm[("openbook2", m)] for m in ("gpt-oss-20b", "gpt-5.4-mini", "claude-sonnet-5")}
claim(ob["gpt-5.4-mini"]["within_margin"] and ob["claude-sonnet-5"]["within_margin"] and not ob["gpt-oss-20b"]["within_margin"],
      "under the open book the two closed models stay within the margin and gpt-oss-20b does not")
tool = {m: arm[("tool", m)] for m in ("claude-sonnet-5", "gpt-5.4-mini")}
claim(all(t["within_margin"] and t["e5"]["within_margin"] for t in tool.values()), "with the tool both closed models stay within the margin on FAC and MC")
fl = {a["model"]: a for a in ANCH["flagship"]["anchors"]}
flr = ANCH["flagship-reasoning-medium"]["anchors"][0]
_roster = ANCH["flagship"]["roster"]
if MATCHED_HEAD:  # the re-run models at their matched setting: their matched-store responses on the subset's items (no interval of its own)
    _sub_ids = set(subsamples.paraphrase_ids())  # the 450-item subsample the anchors and the experiments ran on
    drift(len(_sub_ids) == ANCH["flagship"]["items"], f"subsamples.PARAPHRASE gives {len(_sub_ids)} items; the anchors ran on {ANCH['flagship']['items']}")

    def _matched_on_subset(k: str) -> dict:
        """Score (the mean of template means, as analyze.anchors scores the roster) and level means of k's matched store on the subset."""
        rows = [r for r in map(json.loads, (HERE / "scores" / REASONING_STORE[k] / f"{k}.jsonl").read_text(encoding="utf-8").splitlines())
                if r["item_id"] in _sub_ids]
        drift(len(rows) == len(_sub_ids), f"{REASONING_STORE[k]}/{k}.jsonl holds {len(rows)} of the subset's {len(_sub_ids)} items")
        by_t, by_lv = defaultdict(list), defaultdict(list)
        for r in rows:
            by_t[r["template_id"]].append(r["score"])
            by_lv[r["level"]].append(r["score"])
        return {"model": k, "score": sum(sum(v) / len(v) for v in by_t.values()) / len(by_t), "ci": None,
                "levels": {lv: {"mean": sum(v) / len(v)} for lv, v in by_lv.items()}}
    _roster = [x if x["model"] not in RERUN else _matched_on_subset(x["model"]) for x in _roster]
roster_sub = sorted(_roster, key=lambda x: -x["score"])
top5_sub = roster_sub[:5]
claim({x["model"] for x in top5_sub} == set(TOP5), "the top five on the subset are the top five overall")
sub_scores = [x["score"] for x in top5_sub]
sub_adv = [x["levels"]["Advanced"]["mean"] for x in top5_sub]
anchors_adv = [fl["deepseek-v4-pro"]["levels"]["Advanced"]["mean"], fl["gpt-5.4"]["levels"]["Advanced"]["mean"],
               flr["levels"]["Advanced"]["mean"]]
anchor_scores = [fl["deepseek-v4-pro"]["score"], fl["gpt-5.4"]["score"], flr["score"]]  # the flagship runs on the subset
anchor_top = max(anchor_scores)
anchor_host = [x["model"] for x in top5_sub if x["ci"] and x["ci"][0] <= anchor_top <= x["ci"][1]]  # top-five models whose interval holds the highest run
claim(anchor_top == flr["score"] and len(anchor_host) >= 1, "the highest flagship run, GPT-5.4 with reasoning, lies inside a top-five model's interval")
n_sub = ANCH["flagship"]["items"]
ob_items, ob_templates = ob["gpt-oss-20b"]["items"], ob["gpt-oss-20b"]["templates"]
mc_rise = [ob[m]["e5"]["diff"] for m in ob]

# The experts' reading of wrong answers (B2), as the store still holds it (scored_current.json; amendment 3).
by_model = {m: ERR["by_model"][m] for m in B2_MODELS}
by_level = ERR["by_level"]
maj = ERR["majority_by_model"]
readings_total = ERR["readings"]
items_by_model = {m: sum(ERR_COMP["by_model_level"][m].values()) for m in B2_MODELS}  # wrong answers in the sample per model
items_read_total, items_total = ERR_COMP["items_read"], ERR_COMP["items"]
readers = ERR_COMP["readers"]
assert sum(sum(v.values()) for v in by_level.values()) == readings_total == sum(sum(v.values()) for v in by_model.values())
assert sum(items_by_model.values()) == items_total
level_n = {lv: sum(by_level[lv].values()) for lv in LEVELS}
CALC, FORM, NOERR, INCOMPLETE = CATEGORIES[5][0], CATEGORIES[2][0], CATEGORIES[6][0], CATEGORIES[7][0]
claim(all(by_model[m].get(INCOMPLETE, 0) == 0 for m in B2_MODELS), "no domain expert used the incomplete option")
easy_calc, easy_n = by_level["Easy"].get(CALC, 0), level_n["Easy"]
form_hard = sum(by_level[lv].get(FORM, 0) for lv in ("Intermediate", "Advanced"))
hard_n = level_n["Intermediate"] + level_n["Advanced"]
form_hard_share, form_easy_share = form_hard / hard_n, by_level["Easy"].get(FORM, 0) / easy_n
claude_noerr = maj["claude-sonnet-5"].get(NOERR, 0)
others = [m for m in B2_MODELS if m != "claude-sonnet-5"]
others_calc = [maj[m].get(CALC, 0) for m in others]
others_form = [maj[m].get(FORM, 0) for m in others]
claim(all(maj[m].get(CALC, 0) == max(maj[m].values()) for m in others), "calculation is the largest majority label for the three other models")
claim(maj["claude-sonnet-5"].get(CALC, 0) > items_by_model["claude-sonnet-5"] / 2, "calculation is the majority label for most of Claude Sonnet 5's wrong answers")
claim(all(sorted(maj[m].values())[-2] == maj[m].get(FORM, 0) for m in others), "formula is the second-largest majority label for the three other models")
fleiss_all = ERR["fleiss_all"]
fleiss_models = list(ERR["fleiss_by_model"].values())
level_n_wo = {lv: level_n[lv] - by_level[lv].get(NOERR, 0) for lv in LEVELS}  # readings per level with the no-error readings removed

# The remaining incorrect verdicts of the top five (RESIDUAL_INCORRECT.md).
res_model = {r[0]: r for r in md_table(residual, r"^\| \| n \| <=0\.2%")}
top5_rows = [res_model[k] for k in TOP5]
top5_incorrect = sum(int(r[1]) for r in top5_rows)
top5_near = sum(int(r[2]) for r in top5_rows)
top5_symbolic = sum(int(r[9]) for r in top5_rows)
# The residual report covers the configuration shown when its incorrect counts are the shown ones for every model; otherwise its split
# is printed PENDING (the lost points always come from the scores shown).
RESID_OK = all(int(res_model[k][1]) == Q1[k]["incorrect"] for k in ORDER)
drift(RESID_OK, f"{RESIDUAL_MD}'s incorrect counts differ from the configuration shown for "
      + ", ".join(k for k in ORDER if int(res_model[k][1]) != Q1[k]["incorrect"]) + ": re-run residual_incorrect.py for that configuration")
RESID = "" if RESID_OK else "PENDING-RESIDUAL "  # the prefix of a phrase whose split RESIDUAL_INCORRECT.md does not cover
tpl_rows = residual.split("### Top five models, by template")[1].split("###")[0]
tpl = {r[0]: r for r in md_table(tpl_rows, r"^\| \| n \|") if r[0] != "all"}
one_template_symbolic = max(int(r[9]) for r in tpl.values())
TWO_CHEMICAL = ["template_work_isothermal_virial", "template_adiabatic_flame_temperature"]  # the two templates repaired in round 5
two_chemical = sum(int(tpl[t][1]) for t in TWO_CHEMICAL if t in tpl)
near_rows = md_table(residual.split("### Within 0.2% and still incorrect")[1].split("###")[0], r"^\| template \| n \|")
claim(all(r[1] == r[2] for r in near_rows), "every incorrect verdict within 0.2% is an exact-digits case")
exact_m = re.search(r"Templates with an exact-digits target .*?: (\d+) of 150; incorrect verdicts on them across the roster: (\d+)",
                    residual)
exact_templates, exact_incorrect = int(exact_m.group(1)), int(exact_m.group(2))
form_total = top5_symbolic + top5_near + two_chemical
adv_top5 = [BL[k]["level"]["Advanced"]["mean"] for k in TOP5]
near_templates = [t for t, v in TEMPLATES_READ.items() if v["role"] == "near"]
near_total = sum(int(r[1]) for r in near_rows)


# ----------------------------------------------------------------------------------------------- derived for the text
norm = lambda t: t.replace("template_", "")  # noqa: E731
claim(all(r["domain"] == "thermodynamics" for r in PT if r["template_id"] in TWO_CHEMICAL), "the two chemical templates are thermodynamics")
_dom_wo = defaultdict(lambda: defaultdict(list))
for r in PT:
    if r["template_id"] not in TWO_CHEMICAL:
        _dom_wo[r["model"]][r["domain"]].append(float(r["answer_score"]))
lowest_domain_wo = {k: min(_dom_wo[k], key=lambda d: sum(_dom_wo[k][d]) / len(_dom_wo[k][d])) for k in ORDER}
n_thermo_wo = sum(v == "thermodynamics" for v in lowest_domain_wo.values())
thermo_eight = [REP[k]["domain"]["thermodynamics"] for k in thermo_models]
glm_empty_adv = sum(int(r["unusable"]) for r in PT if r["model"] == "glm-5.3" and r["level"] == "Advanced")
drift(sum(int(r["unusable"]) for r in PT if r["model"] == "glm-5.3") == Q1["glm-5.3"]["empty"], "per_template.csv and results.json disagree on GLM-5.3's empty responses")
branch_of = {r["domain"]: r["branch"] for r in PT}
geo = REP["gpt-oss-20b"]["domain"]["geotechnical_engineering"]
civil_other = sorted(REP["gpt-oss-20b"]["domain"][d] for d in branch_of if branch_of[d] == "civil_engineering" and d != "geotechnical_engineering")
claim(len(civil_other) == 2 and geo < civil_other[0] and lowest_domain["gpt-oss-20b"] == "geotechnical_engineering", "gpt-oss-20b is lowest on geotechnical engineering")
n_branch_templates = BL[ORDER[0]]["branch"]["chemical_engineering"]["templates"]
detect_branch = [BL[k]["detectable_branch"] for k in ORDER]
dec_weights = {d["model_key"]: d["weights"] for d in load(RESULTS / "decoding_table.json")}
OPEN = [k for k in ORDER if dec_weights.get(k) == "open"]
CLOSED = [k for k in ORDER if dec_weights.get(k) == "closed"]
OTHER_WEIGHTS = [k for k in ORDER if k not in OPEN and k not in CLOSED]
drift(not OTHER_WEIGHTS, f"decoding_table.json has no open/closed weights entry for {OTHER_WEIGHTS}")
tol_swaps = {}
for var in ("half_tol", "double_tol"):
    order_var = sorted(ORDER, key=lambda m: -SENS[m][var])
    tol_swaps[var] = [(a, b) for i, a in enumerate(ORDER) for b in ORDER[i + 1:] if order_var.index(a) > order_var.index(b)]
    claim(all(frozenset(pair) not in SEP_FAC for pair in tol_swaps[var]), f"no pair swapped under {var} differs after correction")
n_tol_swaps = {var: len(tol_swaps[var]) for var in tol_swaps}  # pairs of models whose order the variant swaps
mc_correct = {m: COV[m]["coverage_fully_solved"] for m in (first_fac, first_cov)}
claim(all(COV[k]["rho_steps"] < 0 for k in ORDER), "Spearman's rho between coverage and steps is negative for every model")
n_mcnemar_sig = sum(pp["mcnemar_p_holm"] < 0.05 for pp in PAIRS)
n_strict_sig = sum(pp["fully_p_holm"] < 0.05 for pp in PAIRS)
strict_flip = [pp for pp in PAIRS if (pp["p_holm"] < 0.05) != (pp["fully_p_holm"] < 0.05)]
strict_only = [pp for pp in strict_flip if pp["fully_p_holm"] < 0.05]  # separated by strict FAC and not by FAC
fac_only = [pp for pp in strict_flip if pp["p_holm"] < 0.05]  # separated by FAC and not by strict FAC
top5_detect = [Q2[k]["detectable_planned"] for k in TOP5]
ml, gl = r_mini["levels"], r_gem["levels"]
claim(ml["Advanced"]["arm"] - ml["Advanced"]["main"] > (ml["Easy"]["arm"] - ml["Easy"]["main"]), "GPT-5.4 mini's reasoning gain is largest on Advanced")
mini_flag_main, mini_flag_arm = r_mini["digit_flag_rate_fully_solved"]["main"], r_mini["digit_flag_rate_fully_solved"]["arm"]
flag_prec = [float(r[6]) for r in md_table(flags, r"\| model \| flags drawn") if r[0] != "all"]
assert len(flag_prec) == len(ORDER)
claim(router_recall < 0.5, "the judged step check's recall is below one half")
HALL, SETUP = CATEGORIES[0][0], CATEGORIES[1][0]
conceptual = {m: sum(maj[m].get(c, 0) for c in (HALL, SETUP, FORM)) for m in B2_MODELS}
claim(conceptual["gpt-oss-20b"] > max(conceptual[m] for m in others if m != "gpt-oss-20b"), "gpt-oss-20b has the most conceptual errors")
claim(maj["gpt-oss-20b"].get(CALC, 0) < items_by_model["gpt-oss-20b"] / 2 < min(maj[m].get(CALC, 0) for m in others if m != "gpt-oss-20b"),
      "calculation is the majority label for GPT-5.4 mini and Gemma 4 and only the largest for gpt-oss-20b")
_votes = defaultdict(Counter)
for n in ERR["notes"]:
    _votes[(n["model"], n["code"], norm(n["template"]))][n["category"]] += 1
_claude_noerr_t = Counter(t for (m, c, t), v in _votes.items() if m == "claude-sonnet-5" and v.most_common(1)[0][1] >= 2 and v.most_common(1)[0][0] == NOERR)
drift(sum(_claude_noerr_t.values()) == claude_noerr, "the notes' no-error majorities for Claude Sonnet 5 do not add up to majority_by_model")
claude_noerr_top3 = sum(n for _, n in _claude_noerr_t.most_common(3))
_top3 = [t for t, _ in _claude_noerr_t.most_common(3)]
top5_empty = sum(Q1[k]["unusable"] for k in TOP5)  # no readable final answer: empty, or no answer stated
top5_partial = sum(Q1[k]["partial"] for k in TOP5)
top5_partial_points = top5_partial / 2  # half credit each: a half point when the count is odd
top5_incorrect_shown = sum(Q1[k]["incorrect"] for k in TOP5)  # the incorrect verdicts of the configuration shown
lost_points = top5_empty + top5_partial_points + top5_incorrect_shown
drift(abs(lost_points - sum(15 * N_TEMPLATES * (1 - Q1[k]["score"]) for k in TOP5)) < 1e-6, "the first five's verdicts do not add up to their lost points")
_vir = re.search(r"\| incorrect traces \|[^\n]*\n\|[-:| ]+\n\| (\d+) \| (\d+) \|", residual)
virial_n, virial_flow = (int(_vir.group(1)), int(_vir.group(2))) if _vir else (0, 0)
median_sd_zero = sum(Q1[k]["within_sd_quartiles"][1] == 0 for k in ORDER)
claude_adv_readings = ERR_COMP["by_model_level"]["claude-sonnet-5"]["Advanced"] * readers
solved = {k: Q1[k]["templates_no_variance_all_solved"] / N_TEMPLATES for k in ORDER}
e3w_max_model = max(ORDER, key=lambda k: Q3[k]["e3_coverage_on_readable_wrong"])
B2_TEXT = [m for m in ["claude-sonnet-5", "gpt-5.4-mini", "gemma-4-26b-a4b", "gpt-oss-20b"] if m in B2_MODELS] + [m for m in B2_MODELS if m not in
           ["claude-sonnet-5", "gpt-5.4-mini", "gemma-4-26b-a4b", "gpt-oss-20b"]]  # the order the text names them
maj_calc = {m: maj[m].get(CALC, 0) for m in B2_MODELS}
unreadable_share = {k: Q1[k]["unusable"] / N_ITEMS for k in ORDER}  # no readable answer: empty or no stated answer, scored 0


# ----------------------------------------------------------------------------------------------- the WS-C1 to C3 result files
def write_stand_ins() -> None:
    """Stand-ins for the result files WS-C1 to C3 write, under results/stand_in/, every one marked "stand_in": true. Values that
    exist in results.json are copied (the default-configuration scores, Q4, the endpoint table); everything a new analysis
    produces (reasoning-on scores, variant coverages, clause counts, carried precision, depth slopes) is a deterministic,
    plausible placeholder derived from the default values, never a measurement. The files follow the schemas of the C1 to C3
    briefs so the blocks can be developed against them; they are replaced by the real files the moment those exist."""
    STAND_IN_DIR.mkdir(parents=True, exist_ok=True)
    tag = {"stand_in": True, "quick": False,
           "note": "STAND-IN written by paper_results.py --stand-in: placeholder values in the schema of the WS-C brief; not a measurement"}
    arm_shift = {k: arm[("reasoning-medium", k)]["diff"] for k in RERUN if ("reasoning-medium", k) in arm}
    shift = {k: (arm_shift.get(k, 0.03) if k in RERUN else 0.0) for k in ORDER}  # the stand-in gain from reasoning
    half = {k: 0.5 if k in RERUN else 1.0 for k in ORDER}
    models = {}
    for k in ORDER:
        q, c, o, d, q2 = Q1[k], COV[k], Q3O[k], Q3[k], Q2[k]
        s = shift[k]
        models[k] = {"store": REASONING_STORE.get(k, CFG_BY[k]["default_store"]), "fac": q["score"] + s, "fac_ci": [q["ci"][0] + s, q["ci"][1] + s],
                     "letter": "", "unreadable": unreadable_share[k] * half[k], "mc_strict": c["coverage"] + s / 2,
                     "mc_strict_ci": [c["ci"][0] + s / 2, c["ci"][1] + s / 2], "mc_e3": o["e3_all"] + s / 2,
                     "mc_e3_ci": [o["e3_all_ci"][0] + s / 2, o["e3_all_ci"][1] + s / 2],
                     "levels": {lv: BL[k]["level"][lv]["mean"] + s for lv in LEVELS}, "gap": q2["gap"] - s / 2,
                     "gap_ci": [q2["ci"][0] - s / 2, q2["ci"][1] - s / 2], "gap_p_holm": q2["p_welch_holm"],
                     "digit_flag_rate": d["digit_flag_rate_on_fully_solved"] * half[k], "router_flag_rate": d["router_judge_rate_on_fully_solved"] * half[k]}
    m_order = sorted(ORDER, key=lambda k: -models[k]["fac"])

    def shifted_pairs(pairs, value, p_key="p_holm"):
        out = []
        for p in pairs:
            diff = value(p["a"]) - value(p["b"])
            moved = p["a"] in RERUN or p["b"] in RERUN
            p_holm = p[p_key] if not moved else (0.001 if abs(diff) > p["detectable"] else 1.0)
            out.append({"a": p["a"], "b": p["b"], "diff": diff, "ci": [p["ci"][0] + diff - p["diff"], p["ci"][1] + diff - p["diff"]],
                        "p_holm": p_holm, "detectable": p["detectable"]})
        return out

    pairs = shifted_pairs(PAIRS, lambda k: models[k]["fac"])
    for k, letter in letters(m_order, separated(pairs)).items():
        models[k]["letter"] = letter
    mc_pairs = shifted_pairs(CPAIRS, lambda k: models[k]["mc_strict"])
    e3_pairs = [{**p, "diff": Q3O[p["a"]]["e3_all"] - Q3O[p["b"]]["e3_all"]} for p in
                ({"a": p["a"], "b": p["b"], "ci": p["ci"], "p_holm": p["p_holm"], "detectable": p["detectable"]} for p in CPAIRS)]
    paired = {k: {"fac_default": Q1[k]["score"], "fac_reasoning": models[k]["fac"], "change": shift[k], "ci": [shift[k] - 0.02, shift[k] + 0.02],
                  "p_holm": 0.001 if abs(shift[k]) > 0.03 else 0.4, "detectable": 0.03, "ci90": [shift[k] - 0.015, shift[k] + 0.015],
                  "mc_change": shift[k] / 2, "empty_default": Q1[k]["empty"], "empty_reasoning": Q1[k]["empty"] // 2} for k in RERUN}
    matched = {**tag, "config": {k: {"default_store": CFG_BY[k]["default_store"], "reasoning_store": REASONING_STORE.get(k)} for k in ORDER},
               "models": models, "pairs": pairs, "mc_pairs": mc_pairs, "e3_pairs_default": e3_pairs, "paired_change": paired,
               "tau_default_vs_matched": {"tau": kendall_tau(ORDER, m_order), "ci": [kendall_tau(ORDER, m_order) - 0.15, 1.0]}}
    single = {**tag, "n_single": N_SINGLE, "n_others": N_MULTI}
    for k in ORDER:
        s, m = Q4[k]["single_path"], Q4[k]["multi_path"]
        single[k] = {"single": {"all": s["all"], "some": s["some"], "none": s["none"], "ci": {"all": s["all_ci"], "some": s["some_ci"], "none": s["none_ci"]}},
                     "others": {"all": m["all"], "some": m["some"], "none": m["none"], "ci": {"all": m["all_ci"], "some": m["some_ci"], "none": m["none_ci"]}}}
    providers = {**tag, "rows": [{"model": k, "endpoint": e, "rows": v["rows"], "raw_score": v["score"], "unusable": v["unusable"],
                                  "matched_diff": v["matched_diff"], "templates_matched": v["templates_matched"], "few_templates": v["templates_matched"] < 20}
                                 for k in ORDER for e, v in res["reported"]["providers"].get(k, {}).items()]}
    cv_models = {}
    for k in ORDER:
        c, o, d = COV[k], Q3O[k], Q3[k]
        base = {"as_scored": {"all": c["coverage"], "ci": c["ci"], "wrong": c["coverage_wrong"]},
                "matching_only": {"all": o["e3_all"], "ci": o["e3_all_ci"], "wrong": d["e3_coverage_on_readable_wrong"]},
                "route_adjusted": {"all": min(1.0, c["coverage"] + 0.02), "ci": [min(1.0, x + 0.02) for x in c["ci"]], "wrong": min(1.0, c["coverage_wrong"] + 0.02)},
                "intermediate_only": {"all": c["coverage"] - 0.03, "ci": [x - 0.03 for x in c["ci"]], "wrong": c["coverage_wrong"] - 0.03}}
        cv_models[k] = {**base, "matched_store": ({v: {**base[v], "all": base[v]["all"] + shift[k] / 2} for v in base} if k in RERUN else None)}
    cp = lambda: [{"a": p["a"], "b": p["b"], "diff": p["diff"], "ci": p["ci"], "p_holm": p["p_holm"], "detectable": p["detectable"]} for p in CPAIRS]  # noqa: E731
    coverage_variants = {**tag, "rule": {"targets": "STAND-IN: milestone values matched to the answer row's targets under the check's unit factors",
                                         "instances_without_milestones": res["milestones"]["items_without"]},
                         "models": cv_models, "pairs": {"matching_only": cp(), "route_adjusted": cp(), "intermediate_only": cp()},
                         "separation_rule": {"claude_vs_deepseek": {"as_scored": dc["p_holm"], "matching_only": 0.2, "route_adjusted": 0.3, "holds": False}},
                         "verbosity": {"slope_numbers": 0.004, "ci": [0.002, 0.006], "slope_tokens": 0.00001, "ci_tokens": [0.0, 0.00002],
                                       "spearman_within": {k: COV[k]["rho_claims"] for k in ORDER}},
                         "reasoning_tokens": {k: {"q1": COV[k]["coverage"] - 0.02, "q2": COV[k]["coverage"] - 0.01, "q3": COV[k]["coverage"],
                                                  "q4": COV[k]["coverage"] + 0.01} for k in ORDER}}
    variants = ("abs_clause_off", "last_digit_unbounded", "prescribed_relaxed", "per_part_credit")
    offsets = {"abs_clause_off": -0.002, "last_digit_unbounded": 0.004, "prescribed_relaxed": 0.01, "per_part_credit": 0.003}
    sv_models = {k: {"headline": SENS[k]["fitted"], **{v: SENS[k]["fitted"] + offsets[v] for v in variants},
                     "changed": {v: {"up": max(0, round(offsets[v] * N_ITEMS)), "down": max(0, round(-offsets[v] * N_ITEMS))} for v in variants}} for k in ORDER}
    sensitivity_variants = {**tag, "stores": {"main": {"models": sv_models, "tau_with_headline": {v: 0.96 for v in variants}}},
                            "relative_error_bins": {k: {"le_0.2": Q1[k]["correct"] - 40, "0.2_1": 25, "1_5": 15, "gt_5": 0} for k in ORDER}}
    flag_precision_file = {**tag, "models": {k: {"flags": (n := round(Q3[k]["digit_flag_rate_on_fully_solved"] * Q3[k]["fully_solved"])),
                                                  "carried_precision": round(n * 0.4), "other": n - round(n * 0.4)} for k in ORDER},
                           "expert_confirmed": {"slips": flags_slip, "carried_precision": round(flags_slip * 0.4), "other": flags_slip - round(flags_slip * 0.4)}}
    bins = ["1", "2", "3", "4-5", "6+"]
    depth = {**tag, "models": {k: {"bins": {b: {"rate": Q3[k]["by_milestone_count"][b]["wrong_rate"],
                                                "ci": [max(0.0, Q3[k]["by_milestone_count"][b]["wrong_rate"] - 0.03), Q3[k]["by_milestone_count"][b]["wrong_rate"] + 0.03],
                                                "n": Q3[k]["by_milestone_count"][b]["items"]} for b in bins},
                                   "slope": 0.3, "ci": [0.1, 0.5], "p": 0.001, "p_holm": 0.01} for k in ORDER}}
    for name, d in (("matched", matched), ("single_path", single), ("providers", providers), ("coverage_variants", coverage_variants),
                    ("sensitivity_variants", sensitivity_variants), ("flag_precision", flag_precision_file), ("depth_model", depth)):
        (STAND_IN_DIR / f"{name}.json").write_text(json.dumps(d, indent=1), encoding="utf-8")
        print(f"stand-in written: {(STAND_IN_DIR / f'{name}.json').relative_to(REPO)}")


if ARGS.stand_in:
    write_stand_ins()
    sys.exit(0)

MATCHED = result_file("matched")
SINGLE = result_file("single_path")
PROVIDERS = result_file("providers")
COVVAR = result_file("coverage_variants")
SENSVAR = result_file("sensitivity_variants")
FLAGPREC = result_file("flag_precision")
DEPTH = result_file("depth_model")

# The matched configuration: every model at its reasoning store where one exists, else its default store (WS-C1).
M_MODELS = {k: v for k, v in MATCHED["models"].items() if isinstance(v, dict)}
M_ORDER = sorted(M_MODELS, key=lambda k: -M_MODELS[k]["fac"])
EXTRA_MODELS = [k for k in M_ORDER if k not in ORDER]  # a twelfth model: in the matched family only
M_PAIRED = MATCHED.get("paired_change", {})
M_TAU = MATCHED.get("tau_default_vs_matched", {})
assert set(ORDER) <= set(M_MODELS), f"matched.json lacks {sorted(set(ORDER) - set(M_MODELS))}"
assert set(RERUN) <= set(M_PAIRED), f"matched.json's paired_change lacks {sorted(set(RERUN) - set(M_PAIRED))}"
for k in EXTRA_MODELS:
    assert k in NAME, f"no display name for the added model {k}: add it to NAME (and paper_setup.py)"
m_top = [k for k in M_ORDER if "a" in M_MODELS[k]["letter"]]  # the top tier at matched settings: every model sharing the first letter
_stale = [k for k in ORDER if M_MODELS[k].get("store", CFG_BY[k]["default_store"]) == CFG_BY[k]["default_store"] and abs(M_MODELS[k]["fac"] - Q1_DEFAULT[k]["score"]) > 1e-9]
drift(not _stale, f"results.json and matched.json disagree on the default-store FAC of {len(_stale)} models ({', '.join(_stale)}): results.json "
      "predates the re-scored store; re-run analyze.py's main pass before Phase 2")


def share(x: float, n: int = N_ITEMS) -> float:
    """A field that may hold a count or a share (the C1 schema says 'unreadable' without a unit): a value above 1 is a count."""
    return x / n if x > 1 else x


# Single path (WS-C1): per model, the share of templates solved on all, some and no instances, single-path and other templates.
S_MODELS = {k: v for k, v in SINGLE.items() if isinstance(v, dict) and "single" in v}
if HEADLINE == "matched":  # matched.json's Q4 rows, in single_path.json's shape
    S_MODELS = {k: {g: {**{c: Q4[k][src][c] for c in ("all", "some", "none")}, "ci": {c: Q4[k][src][f"{c}_ci"] for c in ("all", "some", "none")}}
                    for g, src in (("single", "single_path"), ("others", "multi_path"))} for k in ORDER}
assert set(ORDER) <= set(S_MODELS), f"single_path.json lacks {sorted(set(ORDER) - set(S_MODELS))}"
drift(SINGLE.get("n_single") == N_SINGLE and SINGLE.get("n_others") == N_MULTI, "single_path.json's template counts differ from results.json's Q4")

# Providers (WS-C1): one row per model and serving endpoint.
P_ROWS = PROVIDERS["rows"] if isinstance(PROVIDERS, dict) else PROVIDERS
P_ROWS = sorted(P_ROWS, key=lambda r: (ORDER_DEFAULT.index(r["model"]) if r["model"] in ORDER_DEFAULT else len(ORDER_DEFAULT), -r["rows"]))
for r in P_ROWS:
    r["few"] = bool(r.get("few_templates") or r.get("few_matched"))  # served or matched fewer than 20 templates
few_endpoints = sum(r["few"] for r in P_ROWS)
P_MODELS = sorted({r["model"] for r in P_ROWS}, key=lambda k: ORDER_DEFAULT.index(k) if k in ORDER_DEFAULT else len(ORDER_DEFAULT))

# Coverage variants (WS-C2).
CV_FAMILY = COVVAR["matched"] if HEADLINE == "matched" else COVVAR  # the coverage readings of the configuration the tables show
CV_MODELS, CV_PAIRS = CV_FAMILY["models"], CV_FAMILY.get("pairs", {})
VARIANTS = ["as_scored", "matching_only", "route_adjusted", "intermediate_only"]
VARIANT_LABEL = {"as_scored": "As scored", "matching_only": "Matching alone", "route_adjusted": "Route-adjusted", "intermediate_only": "Intermediate only"}
assert set(ORDER) <= set(CV_MODELS), f"coverage_variants.json lacks {sorted(set(ORDER) - set(CV_MODELS))}"
CV_LETTERS = {v: letters(sorted(ORDER, key=lambda k: -CV_MODELS[k][v]["all"]), separated(CV_PAIRS[v])) for v in CV_PAIRS if v in VARIANTS}
CV_LETTERS["as_scored"] = CLD_MC
SEP_RULE = CV_FAMILY["separation_rule"]["claude_vs_deepseek"]
VERB = CV_FAMILY["verbosity"]
VERB_UNIT = VERB.get("units", {}).get("numbers", "MC per numeric value shown")
RTOK = CV_FAMILY.get("reasoning_tokens", {})

# Scoring-rule variants (WS-C3): the main store's rows; the reasoning stores' rows print when the headline is matched.
SV_STORE = "matched" if HEADLINE == "matched" and "matched" in SENSVAR["stores"] else "main"  # the configuration Table 1 shows
SV_MODELS = SENSVAR["stores"][SV_STORE]["models"]
SV_TAU = SENSVAR["stores"][SV_STORE]["tau_with_headline"]
SV_VARIANTS = ["abs_clause_off", "last_digit_unbounded", "prescribed_relaxed", "per_part_credit"]
SV_LABEL = {"abs_clause_off": ("Absolute-value", "clause off"), "last_digit_unbounded": ("Last-digit term", "not capped"),
            "prescribed_relaxed": ("Prescribed digits", "relaxed"), "per_part_credit": ("Proportional", "partial credit")}
assert set(ORDER) <= set(SV_MODELS), f"sensitivity_variants.json lacks {sorted(set(ORDER) - set(SV_MODELS))}"
REL_BINS = SENSVAR["stores"][SV_STORE].get("relative_error_bins", SENSVAR.get("relative_error_bins", {}))
REL_KEYS = ["le_0.2", "0.2_1", "1_5", "gt_5"]
REL_LABEL = {"le_0.2": "$\\le$0.2\\%", "0.2_1": "0.2--1\\%", "1_5": "1--5\\%", "gt_5": "$>$5\\%"}

# Carried precision (WS-C3).
FP_MODELS, FP_EXPERT = FLAGPREC["models"], FLAGPREC["expert_confirmed"]
if HEADLINE == "matched":  # the re-run models' flags on their responses with reasoning on
    FP_MODELS = {**FP_MODELS, **{k: FLAGPREC["stores"][REASONING_STORE[k]]["models"][k] for k in RERUN}}
assert set(ORDER) <= set(FP_MODELS), f"flag_precision.json lacks {sorted(set(ORDER) - set(FP_MODELS))}"
drift(FP_EXPERT["slips"] == flags_slip, f"flag_precision.json counts {FP_EXPERT['slips']} confirmed slips; FLAG_REVIEW_3.md {flags_slip}")
fp_total = sum(FP_MODELS[k]["flags"] for k in ORDER)
fp_carried = sum(FP_MODELS[k]["carried_precision"] for k in ORDER)

# Depth, controlled (WS-C3).
D_CONFIG = DEPTH["configurations"]["matched" if HEADLINE == "matched" else "main"]  # the configuration the tables show
D_MODELS = D_CONFIG["models"]
D_BINS = list(D_MODELS[ORDER[0]]["bins"])
assert set(ORDER) <= set(D_MODELS), f"depth_model.json lacks {sorted(set(ORDER) - set(D_MODELS))}"
depth_holds = [k for k in ORDER if D_MODELS[k]["p_holm"] < 0.05]

# The judge swap (judge_swap.py): a second judge on a sample of what the judge is sent, and the MC shift it gives per model.
JUDGE_SWAP = load(RESULTS / "judge_swap_main.json")
SWAP_DIFF = [m["diff"] for m in JUDGE_SWAP["models"].values()]
assert set(JUDGE_SWAP["models"]) <= set(ORDER) and len(SWAP_DIFF) == len(ORDER), "judge_swap_main.json covers other models than the evaluation"


# ----------------------------------------------------------------------------------------------- phrases
low_pairs = [tuple(sorted(pr, key=ORDER.index)) for pr in itertools.combinations(rest, 2) if frozenset(pr) in SEP_FAC]
SENS_OTHER = ["half_unit", "whole_trace", "without_shortcut_templates", "without_symbolic_templates"]
other_shift = max(abs(SENS[k][c] - SENS[k]["fitted"]) for k in ORDER for c in SENS_OTHER)
other_tau = min(sens_tau[c] for c in SENS_OTHER)
claim(sens_tau["unusable_excluded"] < min(v for c, v in sens_tau.items() if c != "unusable_excluded"),
      "excluding the responses without a readable answer reorders the models most")
VIRIAL, FLAME = (TEMPLATES_READ.get(t, {}).get("asks", {}) for t in TWO_CHEMICAL)
if not REPAIRED:
    claim(VIRIAL == {"reading": {"either: the wording does not decide": 3}, "form": {"either: the wording does not decide": 3},
                     "unique": {"no": 3}, "traces": {"all of them": 3}}, "the three chemical experts' B4 answers on the virial template")
    claim(FLAME == {"data": {"no: standard sources differ by more than that": 2, "only if the same data source is used": 1},
                    "method": {"no": 1, "yes": 2}, "traces": {"some of them": 3}}, "the three chemical experts' B4 answers on the flame template")
claim(all(TEMPLATES_READ[t]["asks"] == {"unique": {"yes": 3}, "trace": {"no": 3}} for t in near_templates), "the near-miss templates' B4 answers")

n_open_top = sum(dec_weights.get(k) == "open" for k in TOP5)
claim(dec_weights.get(ORDER[0]) == "open", "the first model is open-weights")
max_none_n = max(Q1[k]["templates_no_variance_none_solved"] for k in ORDER)  # templates a model fails on every instance, at most
max_none = max_none_n / N_TEMPLATES
claim(round(margin * 100) == 5, "the paraphrase margin is five points")
claim(max(flr["score"], fl["deepseek-v4-pro"]["score"], fl["gpt-5.4"]["score"]) <= max(x["ci"][1] for x in top5_sub if x["ci"]), "no anchor run exceeds the top tier")

# Derived for the rewritten text (D-195).
top5_pair_detect = [p["detectable"] for p in PAIRS if {p["a"], p["b"]} <= set(TOP5)]
assert len(top5_pair_detect) == 10
detect_lo, detect_hi = round(min(top5_pair_detect) * 100), round(max(top5_pair_detect) * 100)
assert 0 < detect_lo <= detect_hi <= 12
n_lower_noreason = sum(k in NO_REASONING for k in rest)
if MATCHED_HEAD:
    claim(many_wrong == ORDER[-len(many_wrong):], "the models with more than 100 wrong answers are the last on FAC")
else:
    claim(set(many_wrong) == set(rest), "the models with more than 100 wrong answers are the lower tier")
many_e3w = [Q3[k]["e3_coverage_on_readable_wrong"] for k in many_wrong]
many_floor = [Q3[k]["e3_null_on_readable_wrong"] for k in many_wrong]
many_missing = [Q3[k]["attribution_on_wrong"]["e5_missing"] for k in many_wrong]
many_router = [Q3[k]["attribution_on_wrong"]["router_judge"] for k in many_wrong]
claude_attr = Q3["claude-sonnet-5"]["attribution_on_wrong"]
claim(claude_attr["e5_missing"] < min(many_missing) and claude_attr["router_judge"] < min(many_router), "Claude Sonnet 5's wrong answers carry fewer missing rulings and judged flags than the lower tier's")
low_some = [Q4[k]["single_path"]["some"] for k in rest]
top_some = [Q4[k]["single_path"]["some"] for k in strong6]
if not MATCHED_HEAD:
    claim(max(top_some) < min(low_some), "every lower-tier model solves some but not all instances of more single-path templates than any top-six model")
first5_some, other6_some = [Q4[k]["single_path"]["some"] for k in TOP5], [Q4[k]["single_path"]["some"] for k in ORDER[5:]]
depth6_top5 = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in TOP5]
depth6_low = [Q3[k]["by_milestone_count"]["6+"]["wrong_rate"] for k in rest]
claim(max(depth6_top5) < min(depth6_low), "every lower-tier model is scored 0 more often on deep instances than any top-five model")
hall_readings = sum(by_model[m].get(HALL, 0) for m in B2_MODELS)
remaining_incorrect = top5_incorrect - (form_total if not REPAIRED else top5_symbolic + top5_near)
claim(remaining_incorrect > 0, "incorrect verdicts remain on other templates")
n_top5_responses = len(TOP5) * N_ITEMS
r_reject_share = r_rejected / r_returned
gpt54_flags = (fl["gpt-5.4"]["digit_flag_rate_fully_solved"], flr["digit_flag_rate_fully_solved"])
claim(gpt54_flags[1] < gpt54_flags[0] and fl["gpt-5.4"]["score"] < flr["score"] and "gpt-5.4" in NO_REASONING_ANCHOR, "GPT-5.4 gains with reasoning and sheds flags")
mini_mc = r_mini["e5"]["diff"]
claim(mini_mc > 0 and ml["Advanced"]["diff"] > max(ml["Easy"]["diff"], ml["Intermediate"]["diff"]), "GPT-5.4 mini's MC rises with reasoning, most on Advanced")
tool_share = {m: tool[m]["tool_use"]["share_with_calls"] for m in tool}
mini_tool_flags = tool["gpt-5.4-mini"]["digit_flag_rate_fully_solved"]
claim(mini_tool_flags["arm"] < mini_tool_flags["main"] and mini_flag_arm < mini_flag_main, "GPT-5.4 mini's flags fall with the tool and with reasoning")
claim(all(ob[m]["e5"]["diff"] > 0 for m in ob), "MC rises for all three under the open book")


def f2floor(x: float) -> str:
    return f"{math.floor(x * 100) / 100:.2f}"


def f2ceil(x: float) -> str:
    return f"{math.ceil(x * 100) / 100:.2f}"


branch_span = {k: max(BL[k]["branch"][b]["mean"] for b in BRANCH) - min(BL[k]["branch"][b]["mean"] for b in BRANCH) for k in ORDER}
top_span, low_span = [branch_span[k] for k in TOP5], [branch_span[k] for k in rest]
claim(max(top_span) < min(low_span), "every top-tier model's branch means spread less than any lower-tier model's")
dom_min = {k: min(REP[k]["domain"].values()) for k in ORDER}
top_dom_floor = f2floor(min(dom_min[k] for k in TOP5))  # no top-tier domain mean falls below this
low_dom_ceil = f2ceil(max(dom_min[k] for k in rest))  # every lower-tier model has a domain mean below this
claim(all(dom_min[k] >= float(top_dom_floor) for k in TOP5) and all(dom_min[k] < float(low_dom_ceil) for k in rest), "the domain floors separate the tiers")
claim(min(many_missing) > 0.5, "a missing milestone marks most of the lower tier's wrong answers")
flag_arm = arm[("flagship-reasoning-medium", "gpt-5.4")]
claim(0.4 <= mini_flag_arm / mini_flag_main <= 0.6, "GPT-5.4 mini's flags roughly halve with reasoning")
claim(flag_arm["digit_flag_rate_fully_solved"]["arm"] < 0.5 * flag_arm["digit_flag_rate_fully_solved"]["main"], "GPT-5.4 sheds most of its flags with reasoning")
claim(ORDER[-1] == "gpt-oss-20b" and 0.6 <= tool_share["claude-sonnet-5"] <= 0.7, "gpt-oss-20b is the weakest model; Claude calls the tool on two thirds of the instances")
if not MATCHED_HEAD:
    claim(r_mini["diff"] > max(scores[:5]) - min(scores[:5]) and r_mini["diff"] > max(scores[6:]) - min(scores[6:]), "GPT-5.4 mini's reasoning gain exceeds either tier's spread")

_tie = {k: sorted(b for b in BRANCH if BL[k]["branch"][b]["mean"] - BL[k]["branch"][lowest_branch[k]]["mean"] < 1e-9) for k in ORDER}
lowest_tied = [k for k in ORDER if len(_tie[k]) > 1]  # models whose lowest branch is a tie, counted apart
claim(len(lowest_tied) <= 1 and all(len(_tie[k]) == 2 for k in lowest_tied), "at most one model ties two branches at the bottom")
_lc = Counter(BRANCH[lowest_branch[k]].lower() for k in ORDER if k not in lowest_tied)
_lo = sorted(_lc, key=lambda b: (-_lc[b], b))
lowest_list = ", ".join((f"{b} for {WORD[_lc[b]]} models" if i == 0 else f"{'and ' if i == len(_lo) - 1 else ''}{b} for {WORD[_lc[b]]}")
                       for i, b in enumerate(_lo)) + "".join(  # e.g. "chemical for five models, industrial for four, and civil for one"
    f", with {tt(k)} tied between {' and '.join(BRANCH[b].lower() for b in _tie[k])}" for k in lowest_tied)
tool_lower_with_calls = [m for m in tool if tool[m]["tool_use"]["score_with_calls"] < tool[m]["tool_use"]["score_without_calls"]]
claim(tool_lower_with_calls == ["claude-sonnet-5"], "of the two closed models, only Claude Sonnet 5 scores lower where it calls the tool")

# The error readings the appendix states beyond the main text.
UNIT, SIGN = CATEGORIES[3][0], CATEGORIES[4][0]
unit_total = sum(by_model[m].get(UNIT, 0) for m in B2_MODELS)
sign_total = sum(by_model[m].get(SIGN, 0) for m in B2_MODELS)
n_sign_models = sum(by_model[m].get(SIGN, 0) > 0 for m in B2_MODELS)
hall_gptoss = by_model["gpt-oss-20b"].get(HALL, 0)
claim(unit_total <= 2 and hall_gptoss > hall_readings / 2, "unit errors are nearly absent; most hallucination readings are gpt-oss-20b's")
k_by = ERR["fleiss_by_model"]
k_low_model = min(k_by, key=k_by.get)
claim(k_low_model == "gpt-oss-20b" and conceptual["gpt-oss-20b"] == max(conceptual.values()), "agreement is lowest for gpt-oss-20b, the model with most conceptual errors")
k_others = [v for m, v in k_by.items() if m != k_low_model]
adv_noerr = by_level["Advanced"].get(NOERR, 0)
adv_noerr_share = adv_noerr / level_n["Advanced"]
items_read_by_model = {m: sum(by_model[m].values()) / readers for m in B2_MODELS}
assert all(float(v).is_integer() for v in items_read_by_model.values()), "readings per model are not a multiple of the readers"
items_read_by_model = {m: int(v) for m, v in items_read_by_model.items()}
assert sum(items_read_by_model.values()) == items_read_total
short_models = [m for m in B2_TEXT if items_by_model[m] < max(items_by_model.values())]  # fewer wrong answers than the sample asks
resid_top3 = sorted(tpl, key=lambda t: -int(tpl[t][1]))[:3]  # the three templates holding most of the top five's incorrect verdicts
resid3_total = sum(int(tpl[t][1]) for t in resid_top3)
n_experiments = len({a["arm"] for a in ARMS})
m_change = {k: M_PAIRED[k]["change"] for k in RERUN}
m_top_default = [k for k in ORDER_DEFAULT if "a" in CLD_FAC_DEFAULT[k]]


def stand(*names: str) -> str:
    """The prefix of a phrase whose numbers come from a stand-in or quick-mode file (so it can never be pasted into the paper as is)."""
    return "STAND-IN " if any(n in STAND_IN | QUICK for n in names) else ""


# Derived for the main text of the October submission (WS-D1): the tiers under both configurations, the readable-only reordering, the
# paraphrase bounds, the controlled depth model, the separation rule for Milestone Coverage, and the first five's lost points.
UPPER = ORDER[:6]  # the upper tier at the providers' defaults


def clean_splits(order: list[str], sep: set[frozenset]) -> list[int]:
    """The cut points that split an ordering into two groups with every pair across the cut separated after correction."""
    return [i for i in range(1, len(order)) if all(frozenset((a, b)) in sep for a in order[:i] for b in order[i:])]


claim(clean_splits(ORDER_DEFAULT, SEP_FAC_DEFAULT) == [len(UPPER)], "at the defaults the one clean split of the FAC order falls after the sixth model")
tier_gap = min(Q1[k]["score"] for k in UPPER) - max(Q1[k]["score"] for k in rest)
if not MATCHED_HEAD:
    claim(max(abs(Q5[k]["diff"]) for k in ORDER) < tier_gap, "no paraphrase change reaches the gap between the tiers")
    claim(not any((a in UPPER) != (b in UPPER) for var in tol_swaps for a, b in tol_swaps[var]), "no tolerance variant swaps a pair across the tiers")
    claim(all(k in rest for k in REPEATS) and 3 * max(rep_sd) < tier_gap, "the re-decoded models are lower-tier and their decoding noise is far below the gap")
M_SEP = separated(MATCHED["pairs"])
m_alone = [k for k in M_ORDER if all(frozenset((k, j)) in M_SEP for j in M_ORDER if j != k)]
m_first4 = M_ORDER[:4]
m_below = [k for k in M_ORDER if M_MODELS[k]["fac"] <= M_MODELS["glm-5.3"]["fac"]]
claim(clean_splits(M_ORDER, M_SEP) == [len(M_ORDER) - 1] and m_alone == [M_ORDER[-1]], "at matched settings the one clean split sets the last model apart")
claim(all(frozenset((a, b)) in M_SEP for a in m_first4 for b in m_below), "at matched settings the first four differ from every model at or below GLM-5.3")
mini_m = M_MODELS["gpt-5.4-mini"]["fac"]
mini_level = [k for k in M_ORDER if k != "gpt-5.4-mini" and f3(M_MODELS[k]["fac"]) == f3(mini_m)]
claim(len(mini_level) == 1, "GPT-5.4 mini with reasoning is level with one other model at three decimals")
claim(all(M_PAIRED[k]["p_holm"] < 0.05 and not M_PAIRED[k]["within_margin"] for k in RERUN), "every re-run model gains beyond the margin")
claim(set(sorted(RERUN, key=lambda k: m_change[k])[-2:]) == {"gpt-5.4-mini", "gemma-4-26b-a4b"},
      "GPT-5.4 mini and Gemma 4 26B gain the most from reasoning effort")
order_readable = sorted(ORDER, key=lambda k: -SENS[k]["unusable_excluded"])
n_reordered = max(i for i, k in enumerate(order_readable) if ORDER[i] != k) + 1  # readable-only FAC moves only the first n_reordered models
if MATCHED_HEAD:
    claim(order_readable[0] == "glm-5.3" and ORDER.index("glm-5.3") == n_reordered - 1,
          "readable-only FAC reorders only the models down to GLM-5.3 and moves GLM-5.3 from last of them to first")
else:
    claim([k for k in order_readable if k in rest] == rest and order_readable[0] == "glm-5.3" and ORDER.index("glm-5.3") == len(UPPER) - 1,
          "readable-only FAC reorders the upper tier only and moves GLM-5.3 from sixth to first")
q5_bound = {k: (Q5[k]["ci90"][0] if abs(Q5[k]["ci90"][0]) > abs(Q5[k]["ci90"][1]) else Q5[k]["ci90"][1]) for k in q5_out}
claim(all(abs(v) > margin for v in q5_bound.values()), "the named paraphrase bounds lie outside the margin")
D_POOLED = D_CONFIG["pooled"]
claim(D_POOLED["ci"][0] > 0 and D_POOLED["or_ci"][0] > 1, "pooled over models, the depth slope is positive")
m_gap_sig = [k for k in M_ORDER if M_MODELS[k].get("gap_p_holm", 1) < 0.05]
claim(first_cov == "claude-sonnet-5" and first_fac == "deepseek-v4.1-flash" and (SEP_RULE["as_scored"] >= 0.05) == MATCHED_HEAD
      and SEP_RULE["matching_only"] < 0.05 and SEP_RULE["intermediate_only"] < 0.05 and SEP_RULE["route_adjusted"] >= 0.05,
      "MC separates Claude Sonnet 5 from DeepSeek V4.1 Flash by matching alone and on intermediates, not route-adjusted, and as scored "
      + ("not at matched settings" if MATCHED_HEAD else "at the defaults"))
claim(FP_EXPERT["carried_precision"] == 0, "no expert-confirmed slip is carried precision")
claude_calc = maj["claude-sonnet-5"].get(CALC, 0)
easy_items = sum(ERR_COMP["by_model_level"][m]["Easy"] for m in B2_MODELS)  # the wrong answers behind the Easy readings
digit_top5 = [Q3[k]["digit_flag_rate_on_fully_solved"] for k in TOP5]
# The review of 2026-10-09: the scope of each contrast stated with the claim.
n_gap_zero = sum(Q2[k]["ci"][0] <= 0 for k in ORDER)  # level gaps whose 95% interval includes zero
upper_missing = [Q3[k]["attribution_on_wrong"]["e5_missing"] for k in ORDER if k not in many_wrong and k != "claude-sonnet-5"]  # the other models
m_mini_level = [k for k in TOP5 if k != "gpt-5.4-mini" and frozenset((k, "gpt-5.4-mini")) not in M_SEP]  # first-five models GPT-5.4 mini no longer differs from
claim(all(frozenset((k, "gpt-5.4-mini")) in SEP_FAC_DEFAULT for k in m_mini_level), "GPT-5.4 mini differed from those models at the defaults")
above_margin = [k for k in RERUN if M_PAIRED[k]["ci90"][0] > margin]  # re-run models whose 90% interval lies wholly above the margin
judge_add = [CV_MODELS[k]["as_scored"]["all"] - CV_MODELS[k]["matching_only"]["all"] for k in ORDER]  # the judge's contribution to MC
claim(bool(above_margin), "at least one reasoning gain lies wholly above the margin")
claim(all(dom_min[l] < min(dom_min[t] for t in top_rep) for l in low_rep), "the lower-tier representatives dip below the top-tier ones")


def dom(d: str) -> str:
    return d.replace("_", " ")


# The matched headline (a graded order, no two tiers): the groups its phrases name, and the claims behind them.
n_groups = len({c for v in CLD_FAC.values() for c in v})  # the letters of the tier column
rerun_top = [k for k in RERUN if k in TOP5]  # re-run models among the first five
cv_sep_scored, cv_sep_matching = len(separated(CV_PAIRS["as_scored"])), len(separated(CV_PAIRS["matching_only"]))
close_reps = [k for k in REPRESENTATIVE if k != ORDER[-1]]  # the representative models other than the last
if MATCHED_HEAD:
    claim(len(rerun_top) == 1, "of the re-run models only one reaches the first five")
    claim(max(branch_span[k] for k in close_reps) <= max(top_span) < branch_span[ORDER[-1]] and lowest_branch[ORDER[-1]] == "civil_engineering",
          "three representative models stay within the first five's branch spread; the last is lowest in civil engineering")
    claim(max(top_span) < branch_span[sixth] < min(low_span) and float(low_dom_ceil) <= dom_min[sixth], "the sixth model lies between the first five and the last five")
    claim(cv_sep_matching > cv_sep_scored, "matching alone separates more MC pairs than coverage as scored")
claim(all(g > 0 for g in gaps), "every model's mean is lower on Advanced than on Easy templates")


def pairs_from(pairs: list[tuple[str, str]]) -> str:
    """'b from a1, a2, and a3' when every pair shares its lower model, otherwise 'a1 from b1 and a2 from b2'."""
    lows = {b for _, b in pairs}
    if len(pairs) > 1 and len(lows) == 1:
        return f"{tt(next(iter(lows)))} from {series([tt(a) for a, _ in pairs])}"
    return series([f"{tt(a)} from {tt(b)}" for a, b in pairs])


# The re-run models at their providers' defaults, and the one clean split of the default ordering (the Reasoning Setting paragraph).
claim(all(k in NO_REASONING for k in RERUN), "the re-run models return no reasoning tokens at their providers' defaults")
rerun_default = (f"{series([tt(k) for k in RERUN])} return no reasoning tokens and score "
                 f"{series([f3(Q1_DEFAULT[k]['score']) for k in RERUN])}")
d_split = (clean_splits(ORDER_DEFAULT, SEP_FAC_DEFAULT) or [0])[0]


_omp_default = [  # 5.3.1 under --headline default (the paper reports the matched headline: _omp_matched below)
    # 5.3.1: the tiers at the providers' defaults
    *([f"the {WORD[len(ORDER)]} models form a graded order of {WORD[n_groups]} overlapping groups, with FAC of {f3(min(scores))} to "
       f"{f3(max(scores))}, and {len(SEP_FAC)} of the {N_PAIRS} pairs differ after Holm's correction",
       f"only {tt(m_alone[0])} differs from every other model, while {series([tt(k) for k in m_first4])} differ from {tt('glm-5.3')} and every "
       "model below it"] if MATCHED_HEAD else
      [f"an upper tier of {WORD[len(UPPER)]} models with FAC of {f3(min(scores[:6]))} to {f3(max(scores[:6]))} and a lower tier of "
       f"{WORD[len(rest)]} with {f3(min(scores[6:]))} to {f3(max(scores[6:]))}",
       f"each of which differs from each upper-tier model after Holm's correction over the {N_PAIRS} pairwise tests"]),
    f"The first {WORD[len(TOP5)]}"
    + ("" if MATCHED_HEAD else f", the upper tier without {tt('glm-5.3')}, whose {Q1['glm-5.3']['empty']} empty responses place it "
       f"{ORDINAL[ORDER.index('glm-5.3')]},") + f" lie within {f3(spread5)} of one another, and "
    + ("no pair of them differs" if not top5_sep_pairs else
       f"only one of their {WORD[10]} pairs differs ({tt(top5_sep_pairs[0][0])} above {tt(top5_sep_pairs[0][1])})" if len(top5_sep_pairs) == 1 else
       f"only {WORD[len(top5_sep_pairs)]} of their {WORD[10]} pairs differ")
    + f", although the design detects differences of {WORD[detect_lo]} to {WORD[detect_hi]} points between them",
    f"even on Advanced templates they score {f2(min(adv_top5))} to {f2(max(adv_top5))}",
    (f"no model fails every instance of more than {WORD[max_none_n]} of the {N_TEMPLATES} templates, but the first {WORD[len(TOP5)]} solve "
     f"{pct(min(solved[k] for k in TOP5))} to {pct(max(solved[k] for k in TOP5))} of the templates on every instance and the other "
     f"{WORD[len(ORDER) - len(TOP5)]} {pct(min(solved[k] for k in ORDER[5:]))} to {pct(max(solved[k] for k in ORDER[5:]))}" if MATCHED_HEAD else
     f"no model fails every instance of more than {WORD[max_none_n]} of the {N_TEMPLATES} templates, but the upper {WORD[len(UPPER)]} solve "
     f"{pct(min(solved[k] for k in UPPER))} to {pct(max(solved[k] for k in UPPER))} of the templates on every instance and the lower "
     f"{WORD[len(rest)]} only {pct(min(solved[k] for k in rest))} to {pct(max(solved[k] for k in rest))}"),
    (f"Of the {WORD[len(RERUN)]} models run with reasoning at medium effort, only {tt(rerun_top[0]) if rerun_top else '--'} reaches the first "
     f"{WORD[len(TOP5)]}" if MATCHED_HEAD else
     f"{WORD[n_lower_noreason].capitalize()} of the lower {WORD[len(rest)]} return no reasoning tokens at their providers' defaults"),
    # the derivation behind correct answers
    (f"As scored, MC does not separate {tt('claude-sonnet-5')} from {tt('deepseek-v4.1-flash')} after Holm's correction (Holm-adjusted $p$ "
     f"{pv(SEP_RULE['as_scored'])}), nor does the final answer; matching alone does ({pv(SEP_RULE['matching_only'])}), but not once the "
     f"milestones the judge rules not needed leave the denominator ({pv(SEP_RULE['route_adjusted'])};" if MATCHED_HEAD else
     f"As scored, MC separates {tt('claude-sonnet-5')} from {tt('deepseek-v4.1-flash')}, which the final answer does not (Holm-adjusted $p$ "
     f"{pv(SEP_RULE['as_scored'])}), and so does matching alone ({pv(SEP_RULE['matching_only'])}), but not once the milestones the judge rules "
     f"not needed leave the denominator ({pv(SEP_RULE['route_adjusted'])};"),
    f"The arithmetic check flags a displayed calculation that does not follow from its own operands in up to {pct1(max(digit))} of a "
    "model's correct answers",
    f"Of the {flags_decided} flags a domain expert decided on, {flags_slip} are real slips (precision {f3(flag_precision)}), and a mechanical "
    f"check finds none of the {FP_EXPERT['other']} it can classify to be a rounding carried forward from unrounded values earlier in the "
    "response",
    f"The judged step check flags a further step in up to {pct1(max(router))} of correct-answer responses, but its flags are judged rather than "
    f"verified: in the expert study, {router_alone_tp} of the {router_alone_n} it raises inside correct answers are errors (precision "
    f"{f3(router_alone_precision)})",
    f"in the expert study, both checks together find {f2(router_recall)} of the flawed steps the domain experts marked inside correct answers",
    # the derivation behind wrong answers
    f"On their wrong answers, the {(f'last {WORD[len(many_wrong)]}' if MATCHED_HEAD else f'{WORD[len(many_wrong)]} lower-tier')} models still state {pct(min(many_e3w))} to {pct(max(many_e3w))} of the gold "
    f"milestones, against a chance floor of at most {pct1(max(many_floor))}",
    f"Up to {pct(max(full_cov))} of these wrong answers reach every gold milestone yet end in an answer the check scores incorrect",
    f"A milestone the judge rules missing marks {pct(min(many_missing))} to {pct(max(many_missing))} of the "
    f"{(f'last {WORD[len(many_wrong)]}' + chr(39) + 's') if MATCHED_HEAD else 'lower tier' + chr(39) + 's'} wrong answers and "
    f"{pct(min(upper_missing))} to {pct(max(upper_missing))} of the other {'' if MATCHED_HEAD else 'upper-tier '}models', but only {pct(claude_attr['e5_missing'])} of "
    f"{tt('claude-sonnet-5')}'s, {pct(COV['claude-sonnet-5']['wrong_full_coverage'])} of which reach every milestone",
    (f"On the {N_SINGLE} templates whose instances all follow one derivation, the first {WORD[len(TOP5)]} solve some but not all instances of "
     f"{pct(min(first5_some))} to {pct(max(first5_some))} of the templates and the other {WORD[len(other6_some)]} of {pct(min(other6_some))} to "
     f"{pct(max(other6_some))} (\\autoref{{appendix:results}})" if MATCHED_HEAD else
     f"On the {N_SINGLE} templates whose instances all follow one derivation, the lower tier solves some but not all instances of up to "
     f"{pct(max(low_some))} of the templates and the upper tier of at most {pct(max(top_some))} (\\autoref{{appendix:results}})"),
    # the stability of the tiers
    f"Each of {q5_pairs} instances over {q5_templates} templates has a paraphrase that keeps every number, unit, and symbol in place",
    f"for {WORD[len(q5_within)]} of the {WORD[len(ORDER)]} models it is bounded within $\\pm {margin:.2f}$ at 90\\% confidence ("
    + series([f"{tt(k)} reaches {sgn(q5_bound[k])}" for k in q5_out]) + ")",
    f"decoding {WORD[len(REPEAT_MODELS)]} models ({series([tt(k) for k in REPEAT_MODELS])}) three more times on {rep_items} instances, "
    f"{WORD[rep_items // N_TEMPLATES]} per template, moves FAC by a standard deviation of at most {f3(max(rep_sd))}",
    f"Restricting FAC to the responses that state an answer reorders {f'only the first {WORD[n_reordered]}' if MATCHED_HEAD else 'the upper tier only'} (Kendall's $\\tau$ {f3(sens_tau['unusable_excluded'])} "
    f"with the ordering as scored): {tt('glm-5.3')}, whose {Q1['glm-5.3']['empty']} empty responses then drop out, rises from "
    f"{ORDINAL[ORDER.index('glm-5.3')]} to {ORDINAL[order_readable.index('glm-5.3')]}",
    # the reasoning setting (tab:matched, WS-C1): at the matched headline one sentence closing the tiers paragraph
    *([stand("matched") + f"Reasoning lifts the FAC of {series([tt(k) for k in RERUN])} by {series([f3(m_change[k]) for k in RERUN])}; at their "
       f"providers' defaults, where they return no reasoning tokens, the final answer sorts the models into two tiers, an upper {WORD[d_split]} "
       f"each of which differs from each of a lower {WORD[len(ORDER_DEFAULT) - d_split]} (\\autoref{{tab:matched}})"] if MATCHED_HEAD else
      [stand("matched") + "With reasoning at medium effort, " + series([f"{tt(k)} gains {sgn(m_change[k])}" for k in RERUN])
       + f" on all {thousands(N_ITEMS)} instances (\\autoref{{tab:matched}})",
       stand("matched") + f"At these matched settings, {tt('gpt-5.4-mini')} rises to {f3(mini_m)}, level with {tt(mini_level[0]) if mini_level else '--'}, "
       f"and the two tiers give way to overlapping groups: only {tt(m_alone[0]) if m_alone else '--'} differs from every other model, while "
       f"{series([tt(k) for k in m_first4])} still differ from {tt('glm-5.3')} and every model below it",
       stand("matched") + f"with reasoning on, {tt('gpt-5.4-mini')} no longer differs from {series([tt(k) for k in m_mini_level])}"]),
]

# 5.3.1 at the matched headline, insight first (2026-10-10): the numbers the tables and the results appendix give are left to them, and
# the claims the shorter text makes are checked here.
claim(all(m_change[k] > 0 for k in RERUN), "every model run with reasoning gains")
claim(rerun_top == ["gpt-5.4-mini"], "of the models run with reasoning only GPT-5.4 mini reaches the first five")
claim(COV["claude-sonnet-5"]["wrong_full_coverage"] > 0.5, "more than half of Claude Sonnet 5's wrong answers reach every milestone")
claim(sum(other6_some) / len(other6_some) > sum(first5_some) / len(first5_some),
      "the other six solve some but not all instances of more single-path templates than the first five, on average")
claim(set(REPEAT_MODELS) == set(ORDER[-len(REPEAT_MODELS):]), "the models decoded three more times are the lowest on FAC")
_sep_sorted = [tuple(sorted(p, key=lambda k: -Q1[k]["score"])) for p in SEP_FAC]  # the pairs that differ, higher model first
claim(min(abs(Q1[a]["score"] - Q1[b]["score"]) for a, b in _sep_sorted if a in REPEATS or b in REPEATS) > 3 * max(rep_sd),
      "three standard deviations of decoding noise are below every differing gap that involves a re-decoded model")
claim(not [(a, b) for a, b in _sep_sorted if Q1[a]["score"] + Q5[a]["diff"] < Q1[b]["score"] + Q5[b]["diff"]],
      "the paraphrase changes, as measured, reverse no pair of models that differs")
claim(len(top5_sep_pairs) < 5, "the final answer separates fewer than half of the first five's pairs (it barely ranks them)")
_top5_pairs = ("none of their ten pairs differs" if not top5_sep_pairs else
               "only one of their ten pairs differs" if len(top5_sep_pairs) == 1 else
               f"only {WORD[len(top5_sep_pairs)]} of their {WORD[10]} pairs differ")
_omp_matched = [
    # the ceiling and the graded order
    f"The final answer barely ranks the strongest models: the first {WORD[len(TOP5)]} are near the ceiling of the evaluation set, scoring "
    f"{f2(min(adv_top5))} to {f2(max(adv_top5))} even on Advanced templates, and {_top5_pairs} after Holm's correction",
    f"Across the {WORD[len(ORDER)]} models the order is graded rather than tiered: {len(SEP_FAC)} of the {N_PAIRS} pairs differ, and only "
    f"{tt(m_alone[0]) if m_alone else '--'} differs from every other model",
    f"no model fails every instance of more than {WORD[max_none_n]} of the {N_TEMPLATES} templates, yet the share of templates solved on every "
    f"instance runs from {pct(min(solved.values()))} to {pct(max(solved.values()))}",
    stand("matched") + f"The {WORD[len(RERUN)]} models run with reasoning at medium effort all gain, {tt(rerun_top[0]) if rerun_top else '--'} by "
    f"{f3(m_change[rerun_top[0]]) if rerun_top else '--'} and into the first {WORD[len(TOP5)]}; at their providers' defaults, the order splits "
    f"into {WORD[2]} tiers",
    # the derivation behind correct answers
    f"Neither the final answer nor MC as scored separates {tt('claude-sonnet-5')} from {tt('deepseek-v4.1-flash')}; matching alone does, but "
    "not once the milestones the judge rules unnecessary for each model's route are set aside",
    f"the arithmetic check flags a displayed calculation that does not follow from its own operands in up to {pct1(max(digit))} of a model's "
    "correct answers",
    f"of the {flags_decided} flags a domain expert decided on, {flags_slip} are real slips (precision {f3(flag_precision)}",
    f"These flags miss most flawed steps: in the expert study, the arithmetic and judged step checks together find only {f2(router_recall)} "
    "of those the domain experts marked inside correct answers",  # the judged step flags and the carried-precision split: results appendix
    # the derivation behind wrong answers
    f"the last {WORD[len(many_wrong)]} models' wrong answers still state {pct(min(many_e3w))} to {pct(max(many_e3w))} of the gold milestones, "
    f"against a chance floor of at most {pct1(max(many_floor))}",
    f"up to {pct(max(full_cov))} of them reach every milestone yet end in an incorrect answer",
    f"a milestone ruled missing marks most of these models' wrong answers ({pct(min(many_missing))} to {pct(max(many_missing))}) but only "
    f"{pct(claude_attr['e5_missing'])} of {tt('claude-sonnet-5')}'s, whose wrong answers mostly reach every milestone and still miss the answer",
    f"on the {N_SINGLE} templates whose instances all follow one derivation, the first {WORD[len(TOP5)]} solve some but not all instances of "
    f"{pct(min(first5_some))} to {pct(max(first5_some))} of the templates and the other {WORD[len(other6_some)]} of {pct(min(other6_some))} to "
    f"{pct(max(other6_some))} (\\autoref{{appendix:results}})",
    # the stability of the order
    f"for {WORD[len(q5_within)]} of the {WORD[len(ORDER)]} models the change is bounded within $\\pm {margin:.2f}$ at 90\\% confidence",
    f"decoding the {WORD[len(REPEAT_MODELS)]} lowest models three more times moves FAC by a standard deviation of at most {f3(max(rep_sd))}",
    f"Counting only the responses that state an answer moves {tt('glm-5.3')}, whose {Q1['glm-5.3']['empty']} responses are empty at the output "
    f"ceiling, from {ORDINAL[ORDER.index('glm-5.3')]} to {ORDINAL[order_readable.index('glm-5.3')]}",
]

# 5.3.3 at the matched headline: "harder problems cost the weaker models most" and "by about a fifth".
claim(sum(Q2[k]["gap"] for k in rest) / len(rest) > sum(Q2[k]["gap"] for k in TOP5) / len(TOP5),
      "the last five's level gaps are larger than the first five's, on average")
claim(0.15 <= D_POOLED["odds_ratio"] - 1 <= 0.25, "the pooled odds ratio per milestone is about a fifth above one")
# 5.3.4: "three experiments ... and a fourth adds reasoning at medium effort"
claim(n_experiments == 4 and "reasoning-medium" in {a["arm"] for a in ARMS}, "the subset carries four experiments, one of them reasoning effort")
# 5.4 at the matched headline: "most wrong answers are calculation slips, and conceptual errors grow with difficulty"
claim(sum(maj_calc.values()) > items_read_total / 2, "calculation is the majority label of most of the wrong answers read")
_conc = {lv: sum(by_level[lv].get(c, 0) for c in (HALL, SETUP, FORM)) / level_n[lv] for lv in LEVELS}
claim(_conc["Easy"] < _conc["Intermediate"] <= _conc["Advanced"], "the share of conceptual readings grows from Easy to Intermediate to Advanced")
# The experiments appendix: each answer first, and the full-set reasoning comparison read beside the subset's.
# The branch and domain appendix at the matched headline: the three upper representatives against the last model.
close_dom_floor = f2floor(min(dom_min[k] for k in close_reps))  # every domain mean of the three lies above this
if MATCHED_HEAD:
    claim(branch_span[ORDER[-1]] > 3 * max(branch_span[k] for k in close_reps), "the last model's branch means vary far more than the other three's")
    claim(all(dom_min[k] > float(close_dom_floor) for k in close_reps)
          and all(branch_of[lowest_domain[k]] == "industrial_engineering" for k in close_reps)
          and dom_min[ORDER[-1]] < min(dom_min[k] for k in close_reps) - 0.3,
          "the three upper representatives stay above the floor in every domain, each lowest in an industrial one; the last falls far lower")
# The error-analysis appendix: each paragraph's finding first, with its scope.
claim(min(ORDER.index(m) for m in B2_MODELS) < len(TOP5) and max(ORDER.index(m) for m in B2_MODELS) == len(ORDER) - 1,
      "the models whose wrong answers the experts read span the FAC order, from the first five to the last model")
claim(sign_total < 0.05 * readings_total and n_sign_models == 2 and hall_gptoss == hall_readings and min(k_by.values()) > 0.8
      and all(by_model["gpt-oss-20b"].get(c, 0) > 0 for c in (HALL, SETUP, FORM)),
      "sign errors are rare and in two models; every hallucination reading is gpt-oss-20b's; agreement is high for every model; "
      "gpt-oss-20b's errors span the conceptual categories")
resid_top2 = [t for t in resid_top3 if int(tpl[t][2]) > 0]  # of the three templates, those whose incorrect verdicts lie within 0.2%
top2_near, top2_total = sum(int(tpl[t][2]) for t in resid_top2), sum(int(tpl[t][1]) for t in resid_top2)
claim(0.2 <= resid3_total / top5_incorrect <= 0.3 and resid_top2 == resid_top3[:2] and top2_near > resid3_total / 2
      and all(t in {r[0] for r in near_rows} for t in resid_top2),
      "about a quarter of the first five's incorrect verdicts fall on three templates, and most of those are within 0.2% on the first two, "
      "both prescribing their digits")
sym_top5 = [t for t in tpl if int(tpl[t][9]) > 0]  # the templates holding the first five's symbolic incorrect verdicts
claim(len(sym_top5) == 1 and int(tpl[sym_top5[0]][9]) == top5_symbolic and sym_top5[0] in SYMBOLIC and sym_top5[0] not in BY_NUMBERS,
      "the first five's symbolic incorrect verdicts all fall on one template that also accepts an equivalent expression")
m_others = [k for k in RERUN if k != "gpt-5.4-mini"]  # the other models run with reasoning over all instances
claim(r_mini["e5"]["ci"][0] > 0 and ml["Advanced"]["ci"][0] < min(ml["Easy"]["ci"][1], ml["Intermediate"]["ci"][1]),
      "on the subset GPT-5.4 mini's MC gain holds, and its Advanced interval overlaps those of the other levels")
claim(all(M_PAIRED[k]["items"] == N_ITEMS for k in RERUN), "the matched comparison covers every instance")
claim(M_PAIRED["gemini-3.1-flash-lite"]["p_holm"] < 0.05 and M_PAIRED["gemini-3.1-flash-lite"]["change"] < r_gem["detectable"],
      "Gemini 3.1 Flash-Lite's gain holds over all instances and is smaller than the subset experiment detects")
claim(M_PAIRED["gpt-5.4-mini"]["mc_change_p_holm"] < 0.05 and abs(M_PAIRED["gpt-5.4-mini"]["mc_change"] - r_mini["e5"]["diff"]) < 0.005
      and M_PAIRED["gpt-5.4-mini"]["mc_change"] > 2 * max(abs(M_PAIRED[k]["mc_change"]) for k in m_others),
      "over all instances GPT-5.4 mini's MC gain holds at about the subset's size, well above the other two's")
claim(flag_arm["p_holm"] < 0.05 and flag_arm["diff"] > 0, "GPT-5.4 gains over its own default with reasoning")
claim("claude-sonnet-5" in TOP5 and "claude-sonnet-5" in ORDER_DEFAULT[:len(TOP5)]
      and min(ORDER_DEFAULT.index("gpt-5.4-mini"), ORDER_DEFAULT.index("gpt-oss-20b")) >= len(TOP5),
      "the open book's other two models score below the first five at the providers' defaults")

phrases = [  # 6_results.tex must contain each of these, whitespace aside
    *(_omp_matched if MATCHED_HEAD else _omp_default),
    # 5.3.2: the tested finding first (no branch is hard for every model), then the untested evenness contrast
    f"only one of the {n_branch_pairs} within-model branch comparisons holds after correction ({tt('gpt-oss-20b')}, electrical above civil); "
    f"with {n_branch_templates} templates per branch, differences below {round(min(detect_branch) * 100)} to {round(max(detect_branch) * 100)} "
    "points cannot be told from sampling noise",
    f"the first {WORD[len(TOP5)]} stay within {f3(min(top_span))} to {f3(max(top_span))} of their best branch and above {top_dom_floor} in every "
    f"domain, while {f'the last {WORD[len(rest)]}' if MATCHED_HEAD else 'the lower tier'} spread {f3(min(low_span))} to {f3(max(low_span))} and "
    f"each falls below {low_dom_ceil} in some domain",
    *([f"Thermodynamics is the lowest domain for {WORD[len(thermo_models)]} of the {WORD[len(ORDER)]} models, largely because it holds the "
       f"{WORD[len(TWO_CHEMICAL)]} Advanced chemical templates whose wording does not pin the answer",
       f"without them it is the lowest for {WORD[n_thermo_wo]}"] if not REPAIRED else
      [f"Thermodynamics is the lowest domain for {WORD[len(thermo_models)]} of the {WORD[len(ORDER)]} models"]),
    # 5.3.3: harder problems cost the weaker models most (the level gap, then the depth of the gold derivation)
    f"Every model scores lower on Advanced than on Easy templates (\\autoref{{fig:level_bars}}), by {rng(gaps)}",
    (f"but the gap holds after Holm's correction only for {series([tt(k) for k in gap_sig])} under Welch's $t$-test, the test we read because "
     "scores vary more across Advanced templates" if MATCHED_HEAD else
     f"The gap holds after Holm's correction for {series([tt(k) for k in gap_sig])} under Welch's $t$-test, the test we read because scores "
     "vary more across Advanced templates")
    + (f", and for {WORD[gap_perm]} models under the permutation of level labels" if MATCHED_HEAD else
       f", and at matched settings for {series([tt(k) for k in m_gap_sig])}" + (" alone" if len(m_gap_sig) == 1 else ""))
    + (f", and for no model once the {WORD[len(TWO_CHEMICAL)]} Advanced chemical templates whose wording does not pin the answer are set aside"
       if not REPAIRED else "")
    + " (\\autoref{tab:level_gap})",
    f"on instances with {WORD[6]} or more gold milestones, the share of responses scored 0 reaches {rng(depth6_low)} for "
    f"{f'the last {WORD[len(rest)]}' if MATCHED_HEAD else 'the lower tier'} against at most {f3(max(depth6_top5))} for the first "
    f"{WORD[len(TOP5)]}",
    stand("depth_model") + "over readable responses each further milestone raises the odds of a wrong answer by about a fifth when pooled "
    f"over the models (odds ratio {f2(D_POOLED['odds_ratio'])}, 95\\% interval {f2(D_POOLED['or_ci'][0])} to {f2(D_POOLED['or_ci'][1])}), "
    f"a slope that holds per model for {WORD[len(depth_holds)]} of the {WORD[len(ORDER)]} (\\autoref{{tab:level_gap}})",
    # 5.3.4: the experiments named first, then one paragraph each (the flagship anchors with the reasoning runs, the open book, the open tool)
    f"{WORD[n_experiments - 1].capitalize()} experiments on a {n_sub}-instance subset test whether a larger model, the governing equations, or a "
    "Python interpreter changes these results, and a fourth adds reasoning at medium effort",
    f"the best run, {tt('gpt-5.4')} with reasoning at medium effort ({f3(flr['score'])}), still lies inside the interval of "
    f"{series([tt(k) for k in anchor_host])} on the same {n_sub} instances",
    f"{tt('gpt-5.4')}, which returns no reasoning tokens at its default, gains {f3(flag_arm['diff'])} with reasoning and sheds most of its "
    f"arithmetic flags, and {tt('gpt-5.4-mini')} with the same setting gains {f3(mini_mc)} in MC while its flags on correct answers roughly halve",
    f"{tt('gpt-oss-20b')} answers more correctly ({sgn(ob['gpt-oss-20b']['diff'])}), partly because fewer of its responses are empty at the "
    f"output ceiling, while {tt('claude-sonnet-5')} and {tt('gpt-5.4-mini')} stay within $\\pm {margin:.2f}$",
    f"neither {tt('claude-sonnet-5')} nor {tt('gpt-5.4-mini')} changes in FAC or MC beyond $\\pm {margin:.2f}$, although "
    f"{tt('claude-sonnet-5')} calls it on two thirds of the instances",
    "at 8-bit precision or higher handles its tool calls reliably",  # the endpoint rule (appendix:decoding); D-184's record
    # 5.4 (the sample as the store still holds it: scored_current.json); the per-model and per-level counts are in the figure and tab:errors
    f"Three domain experts read the full response of {items_read_total} wrong answers from {WORD[len(B2_MODELS)]} models that "
    f"{'span the FAC order' if MATCHED_HEAD else 'contrast the tiers'}, all at the providers' default settings, and assigned the first "
    f"category that applies in a {WORD[6]}-category hierarchy from hallucination to calculation, or ``no error'' (Fleiss'~$\\kappa$ "
    f"{f3(fleiss_all)}~\\citep{{fleiss1971}}",
    f"Over {thousands(n_top5_responses)} responses, the first {WORD[len(TOP5)]} models lose {lost_points:g} answer points",
    f"{top5_empty} go to responses with no readable final answer and {top5_partial_points:g} to {top5_partial} partial answers at half credit",
    RESID + f"Of their {top5_incorrect} incorrect verdicts, "
    + (f"{two_chemical} fall on the {WORD[len(TWO_CHEMICAL)]} chemical templates whose wording does not pin the answer, " if not REPAIRED else "")
    + f"{top5_near} lie within 0.2\\% of a target whose digits the question prescribes, which makes them incorrect by its terms, "
    f"{top5_symbolic} are symbolic answers on one template checked for equivalence, and {remaining_incorrect} fall in neither group",
    f"which the arithmetic check flags in {pct1(min(digit_top5))} to {pct1(max(digit_top5))} of their correct answers",
]

# Consistency within a template, stated in the results appendix (the per-model table is not printed): the ranges over the models of
# the templates solved on all and on some instances, single-path and other templates, and the most templates a model fails on all.
def sp_rng(group: str, col: str) -> str:
    vals = [S_MODELS[k][group][col] for k in ORDER]
    return f"{pct(min(vals))} to {pct(max(vals))}"


SP_SRC = "matched" if HEADLINE == "matched" else "single_path"  # the file S_MODELS reads
sp_none = {g: max(round(S_MODELS[k][g]["none"] * n) for k in ORDER) for g, n in (("single", N_SINGLE), ("others", N_MULTI))}
claim(set(sorted(ORDER, key=lambda k: -S_MODELS[k]["single"]["some"])[:2]) == set(ORDER[-2:]),
      "the two lowest models solve some but not all instances of more single-path templates than any other model")
# Carried precision, stated in the results appendix (the per-model table is not printed): the largest share for one model, and the
# confirmed slips the classification can place (the others sit in steps it cannot classify).
fp_max_share = max(FP_MODELS[k]["carried_precision"] / FP_MODELS[k]["flags"] for k in ORDER if FP_MODELS[k]["flags"])
fp_classified = FP_EXPERT["carried_precision"] + FP_EXPERT["other"]
drift(FP_EXPERT["slips"] - fp_classified == FP_EXPERT.get("step_changed", 0),
      f"flag_precision.json: {FP_EXPERT['slips']} confirmed slips, {fp_classified} classified, {FP_EXPERT.get('step_changed', 0)} not classifiable")

# The results appendix, insight first (2026-10-10): the qualitative findings its paragraphs lead with, checked here.
n_adv_more = sum(BL[k]["level"]["Advanced"]["sd"] > BL[k]["level"]["Easy"]["sd"] for k in ORDER)  # Advanced varies more than Easy
claim(len(SEP_FAC) > N_PAIRS / 2 and len(top5_sep_pairs) < 5, "the final answer separates most pairs of models but few at the top")
claim(SEP_MC <= SEP_FAC and len(SEP_MC) < len(SEP_FAC), "MC separates fewer pairs than FAC, none that FAC does not")
claim(not sixth_sep and len(top5_sep_pairs) + cross_sep + len(low_pairs) == len(SEP_FAC),
      "the differing pairs add up: inside the first five, between the first six and the last five, inside the last five")
claim(all(Q3[k]["e3_coverage_on_readable_wrong"] > 2 * Q3[k]["e3_null_on_readable_wrong"] for k in ORDER),
      "every model's wrong answers state more than twice the chance floor of gold milestones")
claim(sum(Q3[k]["attribution_on_wrong"]["e5_missing"] > 0.5 for k in ORDER) > len(ORDER) / 2,
      "for most models, a milestone the judge rules missing marks most wrong answers")
claim(0 < VERB["slope_numbers"] and VERB["ci"][1] < 0.01, "the fixed-effects slope of MC on the values shown is barely positive")
claim(max(VERB["spearman_within"].values()) < 0, "within each model, coverage falls as a response displays more values")
claim(all(x > 0 for x in judge_add), "the judge raises every model's coverage over matching alone")
claim(0 < len(gap_sig) < len(ORDER), "the level gap is model-specific")
# The paraphrase appendix: the margin read against the design's precision, MC beside FAC, and the tau against the noise range.
claim(set(sorted(ORDER, key=lambda k: -Q5[k]["detectable"])[:len(q5_out)]) == set(q5_out),
      "the models whose paraphrase change is not bounded are those the design measures least precisely")
claim(max(abs(Q5[k]["e5"]["diff"]) for k in ORDER) <= max(abs(Q5[k]["diff"]) for k in ORDER), "MC moves no more under paraphrase than FAC does")
claim("p5" in noise and noise["p5"] <= q5_tau["tau"] < noise["q1"], "the paraphrase tau lies within the noise range, below its lower quartile")
claim(set(vs_beyond) == set(q5_out), "the paraphrase change exceeds decoding noise for the same models whose change the margin does not bound")

appendix_phrases = {
    "results": [
        "the smallest difference the design detects at 80\\% power",
        # pairwise comparisons, after the finding (most pairs differ, few at the top; MC separates fewer, none that FAC does not)
        f"Of the {N_PAIRS} pairs, {len(SEP_FAC)} differ on FAC after Holm correction under the sign-flip test over templates: "
        f"{WORD.get(len(top5_sep_pairs), str(len(top5_sep_pairs)))} inside the first {WORD[len(TOP5)]} ({pairs_from(top5_sep_pairs)}), none "
        f"between {tt(sixth)} and the first {WORD[len(TOP5)]}, {cross_sep} of the {len(ORDER[:6]) * len(rest)} between the first "
        f"{WORD[len(ORDER[:6])]} and the last {WORD[len(rest)]}, and {WORD.get(len(low_pairs), str(len(low_pairs)))} inside the last "
        f"{WORD[len(rest)]}",
        f"Every non-significant FAC difference lies below what its pair detects ({rng([pp['detectable'] for pp in nonsig])})",
        f"On MC, {n_mc_sig} pairs differ{', each of them also on FAC' if SEP_MC <= SEP_FAC else ''}: the three highest models, {tt(top3_cov[0])}, "
        f"{tt(top3_cov[1])}, and {tt(top3_cov[2])}, do not separate, {'nor do' if MATCHED_HEAD else 'while'} {tt(first_cov)} and {tt(first_fac)}, "
        f"which FAC does not separate{' either' if MATCHED_HEAD else ', do'}; Kendall's $\\tau$~\\citep{{kendall1938}} between the FAC and MC "
        f"orderings is {f3(tau['tau'])} (95\\% interval {ci(tau['ci'], False)})",
        f"Strict FAC, which scores a partial answer 0, separates {n_strict_sig} pairs and agrees with FAC on {N_PAIRS - len(strict_flip)} of the "
        f"{N_PAIRS}; a Wilcoxon test~\\citep{{wilcoxon1945}} separates {n_wil_sig} on MC; and McNemar's exact test~\\citep{{mcnemar1947}}, which "
        f"treats the {thousands(N_ITEMS)} instances as independent, separates {n_mcnemar_sig}",
        # matched settings (the defaults' scores, the five highest there and the tau are in tab:matched)
        stand("matched") + f"Reasoning raises all three models run with it (Holm-adjusted $p$ {pv(max(M_PAIRED[k]['p_holm'] for k in RERUN))}), "
        f"and for {series([tt(k) for k in above_margin])} the 90\\% interval lies wholly above $+{margin:.2f}$",
        # scoring rules
        "Halving the answer tolerance swaps " + (f"{WORD[n_tol_swaps['half_tol']]} pair{'s' if n_tol_swaps['half_tol'] > 1 else ''} of models"
                                                 if n_tol_swaps["half_tol"] else "no pair of models")
        + " and doubling it swaps " + ("none" if not n_tol_swaps["double_tol"] else
                                        f"one, {tt(tol_swaps['double_tol'][0][0])} and {tt(tol_swaps['double_tol'][0][1])}"
                                        if n_tol_swaps["double_tol"] == 1 else f"{WORD[n_tol_swaps['double_tol']]}")
        + f", and no swapped pair differs after correction ($\\tau = {f3(sens_tau['half_tol'])}$ and ${f3(sens_tau['double_tol'])}$ against the "
        "ordering as scored)",
        f"leaving out the {WORD[len(SHORTCUT)]} templates answerable from their wording or the {WORD[len(SYMBOLIC)]} with symbolic answers move "
        f"no model's FAC by more than {f3(other_shift)} and keep $\\tau$ at {f3(other_tau)} or above",
        f"($\\tau = {f3(sens_tau['unusable_excluded'])}$), chiefly because {tt('glm-5.3')}'s {Q1['glm-5.3']['empty']} empty responses then drop out",
        # coverage (the Claude-DeepSeek p values and the fixed-effects slope are in the caption of tab:coverage_variants)
        stand("coverage_variants") + "within each model coverage falls slightly as a response displays more values (Spearman's "
        f"$\\rho$~\\citep{{spearman1904}} {sgn(min(VERB['spearman_within'].values()), 2)} to {sgn(max(VERB['spearman_within'].values()), 2)}), "
        f"consistent with lower coverage on templates with more steps ($\\rho$ {sgn(max(rho), 2)} to {sgn(min(rho), 2)} between a template's mean "
        "coverage and its number of steps)",
        stand("coverage_variants") + f"The judge raises every model's coverage by {f3(min(judge_add))} to {f3(max(judge_add))} over matching alone "
        f"and narrows the differences between models: matching alone separates {cv_sep_matching} of the {N_PAIRS} pairs, MC as scored "
        f"{cv_sep_scored}, and a second judge on a sample shifts a model's coverage by {sgn(min(SWAP_DIFF))} to {sgn(max(SWAP_DIFF))}",
        stand("flag_precision") + f"{fp_carried} of the {thousands(fp_total)} ({pct1(fp_carried / fp_total)}; at most {pct1(fp_max_share)} for any "
        "model) are carried precision, a printed result that is a correct rounding of the value recomputed from the unrounded values earlier "
        f"in the response, and of the {FP_EXPERT['slips']} slips the domain expert confirmed, "
        + (f"none of the {fp_classified} it can classify is" if not FP_EXPERT["carried_precision"] else
           f"{FP_EXPERT['carried_precision']} of the {fp_classified} it can classify are"),
        # difficulty level and depth (the pooled odds ratio is in the main text, the per-model slopes in tab:level_gap)
        f"over the responses with a readable answer, its gap is {sgn(glm['unusable_excluded']['gap'])} (Holm-adjusted Welch $p$ "
        f"{pv(glm['unusable_excluded']['p_welch_holm'])})",
        f"varies more than the Easy level for {WORD[n_adv_more]} of the {WORD[len(ORDER)]} models; under that permutation the gap holds for "
        f"{WORD[gap_perm]} models",
        stand("depth_model") + f"the logistic slope holds after Holm correction for {WORD[len(depth_holds)]} of the {WORD[len(ORDER)]} models",
        # consistency within a template (the counts behind the Results sentence on single-path templates)
        stand(SP_SRC) + f"separately for the {N_SINGLE} templates whose instances all follow one derivation and for the {N_MULTI} others",
        stand(SP_SRC) + f"Across the models, the share of templates solved on every instance runs from {sp_rng('single', 'all')} on the "
        f"single-path templates and from {sp_rng('others', 'all')} on the others, and the share solved on some but not all from "
        f"{sp_rng('single', 'some')} and from {sp_rng('others', 'some')}",
        stand(SP_SRC) + f"No model fails every instance of more than {WORD.get(sp_none['single'], str(sp_none['single']))} of the {N_SINGLE} "
        f"single-path templates or {WORD.get(sp_none['others'], str(sp_none['others']))} of the {N_MULTI} others",
    ],
    "branch_domain": [
        f"of {WORD[len(REPRESENTATIVE)]} representative models",
        f"the lowest branch differs by model, {lowest_list}", "Within a model, a single branch pair differs after correction",
        (f"{series([tt(k) for k in close_reps])} stay close to their FAC in every branch, while {tt(ORDER[-1])} varies far more and is "
         f"lowest in {BRANCH[lowest_branch[ORDER[-1]]].lower()} engineering" if MATCHED_HEAD else
         f"the {WORD[len(top_rep)]} upper-tier models stay close to their FAC in every branch, while the {WORD[len(low_rep)]} lower-tier models are "
         f"lowest in different branches, {tt(low_rep[0])} in {BRANCH[lowest_branch[low_rep[0]]].lower()} and {tt(low_rep[1])} in "
         f"{BRANCH[lowest_branch[low_rep[1]]].lower()} engineering"),
        (f"the same {WORD[len(close_reps)]} stay above {close_dom_floor} in every domain, each lowest in an industrial one, while "
         f"{tt(ORDER[-1])} falls to {f2(dom_min[ORDER[-1]])} in {dom(lowest_domain[ORDER[-1]])}" if MATCHED_HEAD else
         f"{tt(top_rep[0])} is lowest on {dom(lowest_domain[top_rep[0]])} and {tt(top_rep[1])} on {dom(lowest_domain[top_rep[1]])}, whereas the "
         f"lower-tier models dip lower, {tt(low_rep[0])} in {dom(low_domain[low_rep[0]][0])} and {tt(low_rep[1])} in "
         f"{dom(low_domain[low_rep[1]][0])}"),
    ],
    "paraphrase": [
        f"We select three of every template's 15 instances ({p_selected} instances)", f"A paraphrase must pass {WORD[2]} checks",
        f"a word similarity to the original of at most {COPY}", "with three attempts per instance",
        f"Of the {p_selected} instances, {p_passing} pass the scripted checks and the domain experts keep {r_kept} pairs on {q5_templates} templates",
        f"The {WORD[len(ORDER)]} models answer the kept paraphrases",
        # the margin, read with how precisely the design measures each model: the two outside it are the two it measures least precisely
        f"Read against the $\\pm {margin:.2f}$ margin, the 90\\% interval lies within it for every model but {series([tt(k) for k in q5_out])}, the "
        f"{WORD[len(q5_out)]} models whose paired differences the design measures least precisely: the smallest difference it detects at 80\\% "
        f"power is {series([f3(Q5[k]['detectable']) for k in q5_out])} for them and at most "
        f"{f3(max(Q5[k]['detectable'] for k in ORDER if k not in q5_out))} for the others",
        f"MC changes by at most {f3(max(abs(Q5[k]['e5']['diff']) for k in ORDER))} in either direction",
        f"of the {WORD[len(REPEAT_MODELS)]} models decoded three more times on {rep_items} instances ({series([tt(k) for k in REPEAT_MODELS])}), FAC "
        f"varies by a standard deviation of {rng(rep_sd)} across repeats, and the paraphrase change lies beyond the spread of the repeats' changes "
        f"for {' and '.join(tt(k) for k in vs_beyond)} and within it for the other {WORD[len(vs_within)]}",
        f"the experts reject {r_rejected} pairs that pass every scripted check",
        f"On the {q5_pairs} kept pairs, no change in FAC holds after correction",
        f"their ordering agrees with the originals' at Kendall's $\\tau$ {f3(q5_tau['tau'])} (95\\% interval {ci(q5_tau['ci'], False)}), while two "
        f"disjoint draws of {q5_pairs} original instances with the same template mix agree with a median $\\tau$ of {f3(noise['median'])} "
        f"(quartiles {f3(noise['q1'])} to {f3(noise['q3'])}" + (f"; 5th percentile {f3(noise['p5'])})" if "p5" in noise else ")"),
        f"it covers only the {q5_templates} templates that can be paraphrased without loss",
    ],
    "experiments": [
        f"{WORD[n_experiments].capitalize()} experiments extend the evaluation",
        f"{n_sub}-instance subset (three instances per template)", f"the {WORD[2]} closed models",
        f"For the {ob_templates} templates whose code states their governing equations", f"({ob_items} instances)", f"up to {TOOL_MAX_CALLS} times",
        f"a {TOOL_TIMEOUT} s limit", f"truncated at {thousands(TOOL_OUTPUT_CHARS)} characters",
        f"the 90\\% interval against the $\\pm {margin:.2f}$ margin", "80\\% power",
        f"the change on the {ob['gpt-oss-20b']['usable_in_both']['items']} instances it answered in both runs is "
        f"{sgn(ob['gpt-oss-20b']['usable_in_both']['diff'])}",
        f"{tt('claude-sonnet-5')} on {pct(tool_share['claude-sonnet-5'])} of the instances and {tt('gpt-5.4-mini')} on "
        f"{pct(tool_share['gpt-5.4-mini'])}",
        f"{tt('deepseek-v4-pro')} scores {f3(fl['deepseek-v4-pro']['score'])} and {tt('gpt-5.4')} {f3(fl['gpt-5.4']['score'])} at its default "
        f"and {f3(flr['score'])} with reasoning, against {rng(sub_scores)} for the first {WORD[len(TOP5)]} on the same instances; the highest "
        f"anchor run lies inside the interval of {series([tt(k) for k in anchor_host])} on these instances, and on Advanced templates the anchors score "
        f"{rng(anchors_adv, f2)} and the first {WORD[len(TOP5)]} {rng(sub_adv, f2)}",
        f"raises its MC by {f3(r_mini['e5']['diff'])} (95\\% interval {ci(r_mini['e5']['ci'], False)})",
        "at 8-bit precision or higher handles its tool calls reliably",  # the endpoint rule (appendix:decoding); D-184's record
        # the full-set reasoning-on runs read beside the subset's (tab:matched, WS-C1)
        f"Over all {thousands(N_ITEMS)} instances, the matched settings",
        stand("matched") + f"while that of the other {WORD[len(m_others)]} moves by at most "
        f"{f3(max(abs(M_PAIRED[k]['mc_change']) for k in m_others))}",
    ],
    "errors": [
        f"The {WORD[6]} categories follow the stages", f"answers the {WORD[6]} questions in this order", f"{WORD[2]} further options",
        f"We read {items_read_total} wrong answers from {WORD[len(B2_MODELS)]} models that span the FAC order ({series([tt(m) for m in B2_TEXT])})",
        f"we draw up to {max(items_by_model.values())} per model at random",
        f"{ERR_COMP['templates']} templates in all",
        f"{WORD[3].capitalize()} domain experts of the wrong answer's branch",
        *([f"{WORD[3].capitalize()} chemical domain experts read the {WORD[len(TWO_CHEMICAL)]} Advanced chemical templates",
           f"all {WORD[3]} answer that the wording decides neither between the closed-system and the flow reading",
           f"nor between the {WORD[2]} forms of the truncated virial equation",
           f"{WORD[FLAME.get('data', {}).get('no: standard sources differ by more than that', 0)]} answer that standard data sources differ by more than the "
           f"tolerance and one that they agree only within one source, {WORD[FLAME.get('method', {}).get('yes', 0)]} that the gold method is the standard one and "
           f"one that it is not, and all {WORD[3]} that some of the responses shown are correct under some reading",
           f"with and without the {WORD[len(TWO_CHEMICAL)]} chemical templates"] if not REPAIRED and FLAME and VIRIAL else []),
        f"On each of {WORD[len(near_templates)]} templates where many wrong answers land between 0.2\\% and 5\\% of the gold value",
        f"{WORD[3].capitalize()} domain experts of each template's branch read these templates",
        f"for each, all {WORD[3]} answer that the question has one correct answer and that the response shown, within 5\\% of the gold "
        "value, is not correct",
        f"{exact_templates} templates prescribe the digits of the answer, so a value within the tolerance but not at the prescribed digits "
        f"is incorrect by the question's own terms ({exact_incorrect} incorrect verdicts across the {WORD[len(ORDER)]} models, {near_total} of "
        f"them within 0.2\\% of the target on {len(near_rows)} of these templates)",
        f"sign errors rare, in {WORD[n_sign_models]} models only; every hallucination reading is {tt('gpt-oss-20b')}'s; and agreement, high for "
        f"every model, is lowest for {tt(k_low_model)}",
        RESID + f"About a quarter of the first {WORD[len(TOP5)]}'s incorrect verdicts fall on {WORD[len(resid_top3)]} templates",
        RESID + f"of the {top5_incorrect}, " + series([f"{int(tpl[t][1])} {'fall ' if i == 0 else ''}on {code(t)}" for i, t in enumerate(resid_top3)])
        + f"; on the first {WORD[len(resid_top2)]}, {top2_near} of the {top2_total} lie that close",
        f"the {WORD[len(SYMBOLIC)]} templates with symbolic answers",
        (f"{WORD[len(BY_NUMBERS)]} of the {WORD[len(SYMBOLIC)]} templates with symbolic answers are scored by the numbers they state"
         if BY_NUMBERS != SYMBOLIC else f"the {WORD[len(SYMBOLIC)]} templates with symbolic answers are scored by the numbers they state"),
        RESID + f"the first {WORD[len(TOP5)]}'s symbolic incorrect verdicts all fall on one of these {WORD[len(SYMBOLIC) - len(BY_NUMBERS)]}, "
        f"{code(sym_top5[0])}",
        f"{WORD[len(SHORTCUT)].capitalize()} templates whose answers can be read off the question's wording stay in the evaluation set, since "
        f"leaving them out moves FAC by at most {f3(short_shift)}",
        # the scoring-rule variants and the carried-precision split (WS-C3)
        stand("sensitivity_variants") + "Changing one clause of the rule at a time (the absolute-value clause, the last-digit term, the prescribed "
        "digits, or proportional partial credit) moves a model's FAC by at most "
        f"{f3(max(abs(SV_MODELS[k][v] - SV_MODELS[k]['headline']) for k in ORDER for v in SV_VARIANTS))} and keeps "
        f"Kendall's $\\tau$ with the ordering as scored at {f3(min(SV_TAU[v] for v in SV_VARIANTS))} or above (\\autoref{{tab:scoring_variants}})",
    ],
}


# ----------------------------------------------------------------------------------------------- tables
# label -> (the file that carries it by default, the block). A block is written wherever in the tree its markers are; the file is
# the default the writers start from (WS-D1 and D2 place the markers in Phase 2).
blocks: dict[str, tuple[str, str]] = {}


def put(label: str, file_key: str, body: str) -> None:
    assert file_key in FILES and label not in blocks
    blocks[label] = (file_key, body)


def mark(k: str, rerun: bool = False) -> str:
    s = tt(k)
    if k in NO_REASONING and not rerun:
        s += "$^{\\ast}$"
    if k == "glm-5.3":
        s += "$^{\\dagger}$"
    if rerun:
        s += "$^{\\ddagger}$"
    return s


def val(x, c=None, fmt=f3) -> str:
    """A value with its interval when both exist, the value alone when only it does, a dash when neither does."""
    if x is None:
        return "--"
    return with_interval(x, c) if c else fmt(x)


def row_values(k: str, src: str) -> dict:
    """Table 1's cells for one model: from results.json (src 'default') or from matched.json (src 'matched')."""
    if src == "default":
        q, d, c, o = Q1[k], Q3[k], COV[k], Q3O[k]
        return {"fac": q["score"], "fac_ci": q["ci"], "letter": CLD_FAC[k], "solved": solved[k], "unreadable": unreadable_share[k],
                "mc5": c["coverage"], "mc5_ci": c["ci"], "mc3": o["e3_all"], "mc3_ci": o["e3_all_ci"], "judged": d["e5_judged_fraction"],
                "claims": d["claims_per_trace"],
                "flags": {kind: (d[f"{kind}_flag_rate_on_fully_solved"], d.get(f"{kind}_ci")) for kind in FLAG_KINDS},
                "router": (d["router_judge_rate_on_fully_solved"], d["router_judge_ci"])}
    m = M_MODELS[k]
    return {"fac": m["fac"], "fac_ci": m["fac_ci"], "letter": m["letter"], "solved": m.get("all15_share"),  # matched.json: the share of templates solved on all 15 instances
            "unreadable": share(m["unreadable"]), "mc5": m["mc_strict"], "mc5_ci": m["mc_strict_ci"], "mc3": m["mc_e3"], "mc3_ci": m["mc_e3_ci"],
            "judged": m.get("judge_decided_share", m.get("e5_judged_fraction")), "claims": m.get("claims_per_trace"),
            "flags": {kind: (m.get(f"{kind}_flag_rate"), m.get(f"{kind}_flag_ci")) for kind in FLAG_KINDS},
            "router": (m.get("router_flag_rate"), m.get("router_flag_ci"))}


def main_row(k: str, src: str, rerun: bool = False) -> str:
    """Table 1: FAC, the all-instances share, the unreadable share, MC, the flags and the calculations parsed; the tier letters stay in
    tab:matched, and MC by matching alone, the judge-decided share and the judged step flags in tab:coverage_variants and tab:coverage."""
    v = row_values(k, src)
    cells = [mark(k, rerun), with_interval(v["fac"], v["fac_ci"]), val(v["solved"]), pct1(v["unreadable"]), with_interval(v["mc5"], v["mc5_ci"])]
    for kind in FLAG_KINDS:
        cells.append(val(*v["flags"][kind]))
    cells.append(val(v["claims"], fmt=f1))
    if JUDGED_IN_TABLE:
        cells.append(val(*v["router"]))
    return " & ".join(cells)


N_MAIN_COLS = 5 + len(FLAG_KINDS) + 1 + (1 if JUDGED_IN_TABLE else 0)
STAR_NOTE = " $^{\\ast}$No reasoning setting offered." if MATCHED_HEAD else " $^{\\ast}$No reasoning tokens at the provider's default."
# The dagger note counts the empty responses of the configuration the table shows (matched.json carries empty_n per model).
glm_empty_shown = M_MODELS["glm-5.3"].get("empty_n", Q1["glm-5.3"]["empty"]) if HEADLINE == "matched" else Q1["glm-5.3"]["empty"]
UP = " $\\uparrow$"
_deriv = ["\\textbf{MC" + UP + "}"]
_deriv += [mk(FLAG_LABEL.get(kind, kind) + ",", "correct answers") for kind in FLAG_KINDS] + [mk("Calculations parsed", "per response")]
if JUDGED_IN_TABLE:
    _deriv.append(mk("Judged step flags,", "correct answers"))
MAIN_HEAD = (" & \\multicolumn{3}{c}{\\textbf{Final answer}} & \\multicolumn{" + str(len(_deriv)) + "}{c}{\\textbf{Derivation}} \\\\\n"
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-" + str(4 + len(_deriv)) + "}\n"
             "\\textbf{Model} & \\textbf{FAC" + UP + "} & " + mk("All 15 instances", "solved" + UP) + " & "
             + mk("No readable", "answer") + " & " + " & ".join(_deriv))
MAIN_SPEC = "l c c r c " + "c " * len(FLAG_KINDS) + "r" + (" c" if JUDGED_IN_TABLE else "")
if HEADLINE == "default":
    rows = [group_row("Open-weights LLMs", N_MAIN_COLS)] + [main_row(k, "default") for k in OPEN] \
        + [group_row("Closed LLMs", N_MAIN_COLS)] + [main_row(k, "default") for k in CLOSED] \
        + ([group_row("Other", N_MAIN_COLS)] + [main_row(k, "default") for k in OTHER_WEIGHTS] if OTHER_WEIGHTS else [])
    if RERUN:
        rows += [group_row("With reasoning at medium effort (\\autoref{tab:matched})", N_MAIN_COLS)] \
            + [main_row(k, "matched", rerun=True) for k in sorted(RERUN, key=lambda k: -M_MODELS[k]["fac"])]
    main_rows_note = (f" Bottom block: the {WORD[len(RERUN)]} models that accept a reasoning setting, on the same instances with reasoning at "
                      "medium effort." if RERUN else "")
    main_setting = "at the providers' default settings"
    main_note = note("matched") if RERUN else ""
else:
    m_open, m_closed = [k for k in M_ORDER if k in OPEN], [k for k in M_ORDER if k in CLOSED]
    m_other = [k for k in M_ORDER if k not in OPEN and k not in CLOSED]
    rows = [group_row("Open-weights LLMs", N_MAIN_COLS)] + [main_row(k, "matched", rerun=k in RERUN) for k in m_open] \
        + [group_row("Closed LLMs", N_MAIN_COLS)] + [main_row(k, "matched", rerun=k in RERUN) for k in m_closed] \
        + ([group_row("Added model", N_MAIN_COLS)] + [main_row(k, "matched", rerun=k in RERUN) for k in m_other] if m_other else [])
    main_rows_note = " $^{\\ddagger}$Run with reasoning at medium effort."
    main_setting = "at matched settings"
    main_note = note("matched")
# The caption says only what reading the table needs: what each column holds, and the marks.
_flags = series([f"the share of correct answers with a step the {'arithmetic' if kind == 'digit' else kind} check flags" for kind in FLAG_KINDS]
                + (["the share with a step the judge flags"] if JUDGED_IN_TABLE else []))
put("tab:main_results", "main", table(
    MAIN_SPEC, MAIN_HEAD, rows,
    main_note + "\\textbf{Overall performance of the evaluated LLMs on \\ourdataset.} "
    f"Final Answer Accuracy (FAC) and Milestone Coverage (MC) {main_setting}, with 95\\% intervals over templates in gray; the share of "
    "templates solved on every instance; the share of responses that are empty or state no answer, scored 0; " + _flags
    + " (flags, not error counts); and the calculations the check parses per response." + main_rows_note + STAR_NOTE
    + f" $^{{\\dagger}}${glm_empty_shown} responses empty at the output ceiling.",
    "tab:main_results", shade_header=False, size="\\small\n\\setlength{\\tabcolsep}{4pt}", resize=True))

# The matched family: both configurations side by side, the paired change for the re-run models.
rows = []
for k in ORDER + EXTRA_MODELS:
    m = M_MODELS[k]
    d_cells = [with_interval(Q1_DEFAULT[k]["score"], Q1_DEFAULT[k]["ci"]), CLD_FAC_DEFAULT[k]] if k in Q1_DEFAULT else ["--", "--"]
    setting = (m.get("setting") or ("reasoning, medium" if k in RERUN else CFG_BY.get(k, {}).get("reasoning_setting", "default"))).replace("=", " ")
    if k in M_PAIRED:
        p = M_PAIRED[k]
        paired = [f"{sgn(p['change'])} ({ci(p['ci'])})", pv(p["p_holm"]), ci(p["ci90"]), sgn(p["mc_change"]), f"{p['empty_default']} / {p['empty_reasoning']}"]
    else:
        paired = ["--"] * 5
    rows.append(" & ".join([mark(k, rerun=k in RERUN)] + d_cells + [with_interval(m["fac"], m["fac_ci"]), m["letter"], setting] + paired))
header = (" & \\multicolumn{2}{c}{\\textbf{Providers' defaults}} & \\multicolumn{3}{c}{\\textbf{Matched settings}} & \\multicolumn{5}{c}{\\textbf{Paired change, "
          "reasoning minus default}} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}\\cmidrule(lr){7-11}\n"
          "\\textbf{Model} & \\textbf{FAC} & \\textbf{Tier} & \\textbf{FAC} & \\textbf{Tier} & \\textbf{Setting} & " + mk("Change", "(95\\% interval)")
          + " & \\textbf{$p$ (Holm)} & \\textbf{90\\% interval} & " + mk("MC", "change") + " & " + mk("Empty responses,", "default / reasoning"))
put("tab:matched", "results", table(
    "l c c c c l c r c r c", header, rows,
    note("matched") + "\\textbf{Matched settings.} "
    f"FAC at the providers' default settings and at matched settings, with the tier letters of each configuration's own "
    f"{N_PAIRS}-pair family (Holm-corrected sign-flip tests over templates; models sharing a letter do not differ). Setting: what the matched "
    f"configuration runs. For the {WORD[len(RERUN)]} models run with reasoning on, the change from their default rows over the same {thousands(N_ITEMS)} "
    "instances, paired by instance: with its 95\\% interval, the Holm-adjusted $p$ of the sign-flip test over templates, the 90\\% interval read "
    f"against the $\\pm {margin:.2f}$ margin, the change in MC, and the empty responses in each configuration. Kendall's $\\tau$ between the two "
    f"orderings is {f3(M_TAU['tau'])}" + (f" (95\\% interval {ci(M_TAU['ci'], False)})" if M_TAU.get("ci") else "")
    + "." + STAR_NOTE + " $^{\\ddagger}$Run with reasoning at medium effort.",
    "tab:matched", size="\\footnotesize", resize=True))

# Milestone Coverage variants.
rows = []
for k in ORDER:
    cells = [tt(k)]
    for v in VARIANTS:
        d = CV_MODELS[k][v]
        cells += [with_interval(d["all"], d["ci"]), f3(d["wrong"]), CV_LETTERS.get(v, {}).get(k, "--")]
    rows.append(" & ".join(cells))
header = (" & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{\\textbf{{{VARIANT_LABEL[v]}}}}}" for v in VARIANTS) + " \\\\\n"
          + "".join(f"\\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}" for i in range(len(VARIANTS))) + "\n\\textbf{Model} & "
          + " & ".join("\\textbf{All} & \\textbf{Wrong} & \\textbf{Tier}" for _ in VARIANTS))
put("tab:coverage_variants", "results", table(
    "l " + "c r c " * len(VARIANTS), header, rows,
    note("coverage_variants") + "\\textbf{Milestone Coverage under four definitions.} Per model: as scored (the judge's rulings added to the matcher's), "
    "by matching alone, route-adjusted (milestones the judge rules not needed leave the denominator), and intermediate only (the milestones that "
    "state the final answer leave the set). All: the mean over all responses with milestones, with its 95\\% interval over templates; wrong: the "
    "mean over wrong answers; tier: models sharing a letter do not differ after Holm's correction over the pairwise sign-flip tests of that "
    f"definition. Intermediate only: the answer milestones are those whose value equals a number of the answer check's targets under its unit "
    f"factors; once they leave the set, {COVVAR['rule']['instances_without_milestones']} instances have no milestone left ("
    f"{COVVAR['rule'].get('no_milestones_at_all', res['milestones']['items_without'])} have none at all). {tt('claude-sonnet-5')} against {tt('deepseek-v4.1-flash')}: Holm-adjusted $p$ {pv(SEP_RULE['as_scored'])} as scored, "
    f"{pv(SEP_RULE['matching_only'])} by matching alone, {pv(SEP_RULE['route_adjusted'])} route-adjusted"
    + (f", {pv(SEP_RULE['intermediate_only'])} intermediate only" if MATCHED_HEAD else "") + ". The slope of MC on the number of values a "
    f"response displays, with template fixed effects, is {sgn(VERB['slope_numbers'], 4)} ({VERB_UNIT}; 95\\% interval {ci(VERB['ci'], False, 4)}).",
    "tab:coverage_variants", size="\\footnotesize", resize=True))

# Scoring-rule variants.
rows = []
for k in ORDER:
    m = SV_MODELS[k]
    cells = [tt(k), f3(m["headline"])]
    for v in SV_VARIANTS:
        ch = m["changed"][v]
        cells.append(f"{f3(m[v])} ({ch['up']}/{ch['down']})")
    rb = REL_BINS.get(k, {})
    cells += [str(rb.get(b, "--")) for b in REL_KEYS]
    rows.append(" & ".join(cells))
rows += ["\\midrule", "Kendall's $\\tau$ with FAC as scored & & " + " & ".join(f3(SV_TAU[v]) for v in SV_VARIANTS) + " & & & &"]
header = (" & & \\multicolumn{4}{c}{\\textbf{FAC under the variant (verdicts up/down)}} & \\multicolumn{4}{c}{\\textbf{Relative error of accepted "
          "answers}} \\\\\n\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}\n\\textbf{Model} & " + mk("FAC", "as scored") + " & "
          + " & ".join(mk(*SV_LABEL[v]) for v in SV_VARIANTS) + " & " + " & ".join(f"\\textbf{{{REL_LABEL[b]}}}" for b in REL_KEYS))
put("tab:scoring_variants", "results", table(
    "l r c c c c r r r r", header, rows,
    note("sensitivity_variants") + "\\textbf{Scoring-rule variants.} The final-answer rule re-applied to every response of the evaluation with one clause changed at a time: the "
    "absolute-value clause off, the response's last-digit term not capped at one hundredth of the target, the exact-digit requirement of the prescribing templates "
    "relaxed to the tolerance, and partial credit in proportion to the parts matched instead of one half. Each cell gives FAC under the variant and "
    "the verdicts that rise and fall; the last row Kendall's $\\tau$ of the variant's ordering with the ordering as scored. Right: the accepted "
    "answers by their relative error to the target.",
    "tab:scoring_variants", size="\\footnotesize", resize=True))

# Serving endpoints. At the matched headline the re-run models' responses come from their reasoning store's endpoints
# (matched.json's providers_reasoning, analyze.export_providers on that store); a model one endpoint served is not listed.
P_ROWS_SHOWN = P_ROWS
if HEADLINE == "matched":
    _mrows = [{**r, "few": r["few_templates"] or r["few_matched"]} for r in result_file("matched")["providers_reasoning"]]
    _by_model = {}
    for r in [r for r in P_ROWS if r["model"] not in RERUN] + _mrows:
        _by_model.setdefault(r["model"], []).append(r)
    P_ROWS_SHOWN = [r for k in ORDER for r in _by_model.get(k, []) if len(_by_model[k]) > 1]
    P_MODELS = [k for k in ORDER if len(_by_model.get(k, [])) > 1]
    few_endpoints = sum(r["few"] for r in P_ROWS_SHOWN)
rows = []
last = None
for r in P_ROWS_SHOWN:
    name = tt(r["model"]) if r["model"] != last else ""
    last = r["model"]
    rows.append(f"{name}{'$^{\\ddagger}$' if name and HEADLINE == 'matched' and r['model'] in RERUN else ''} & "
                f"{r['endpoint']}{'$^{\\S}$' if r['few'] else ''} & {r['rows']} & {f3(r['raw_score'])} & {pct1(r['unusable'])} & "
                f"{sgn(r['matched_diff']) if r['matched_diff'] is not None else '--'} & {r['templates_matched']}")
put("tab:providers", "results", table(
    "l l r r r r r", "\\textbf{Model} & \\textbf{Endpoint} & \\textbf{Responses} & \\textbf{FAC} & " + mk("No readable", "answer") + " & "
    + mk("Matched", "difference") + " & " + mk("Templates", "matched"), rows,
    note("providers") + "\\textbf{Serving endpoints.} "
    f"The endpoints behind the responses of the {WORD[len(P_MODELS)]} models that more than one endpoint served: the "
    "responses each served, their FAC and their share with no readable answer, and the template-matched difference, the endpoint's mean minus the "
    "model's mean on the same templates, over the templates both cover. The raw FAC per endpoint is not comparable across endpoints, since "
    f"OpenRouter assigns instances to endpoints unevenly; $^{{\\S}}$marks the {few_endpoints} endpoints that served or matched fewer than 20 templates, whose "
    "difference rests on too few templates to read."
    + (" $^{\\ddagger}$Responses with reasoning at medium effort; a model that one endpoint served is not listed." if HEADLINE == "matched" else ""),
    "tab:providers", size="\\footnotesize", star=False))

# The Limitations and the Conclusion (main text, no generated blocks): their numbers, from the values above. Defined here, after
# the serving endpoints the Limitations bound.
prov_max = max(abs(r["matched_diff"]) for r in P_ROWS_SHOWN if not r["few"] and r["matched_diff"] is not None)
judge_share = [M_MODELS[k]["judge_decided_share"] for k in ORDER]  # the judge-decided share of milestones, matched settings
claim(frozenset(("gpt-5.4-mini", "muse-glimmer-30b")) not in SEP_FAC, "the two models stricter scoring swaps do not differ on FAC")
claim(min(Q1[k]["score"] for k in TOP5) > 0.95 and max(Q1[k]["score"] for k in TOP5) - min(Q1[k]["score"] for k in TOP5) <= 0.02,
      "the first five lie within two points of each other, near the ceiling")
appendix_phrases["limitations"] = [
    f"cover undergraduate curricula in {WORD[len(BRANCH)]} branches",
    f"where {WORD[len(RERUN)]} models reason at medium effort and the endpoint of {tt('qwen3-235b-a22b-2507')} offers no reasoning setting",
    f"not among the evaluated {WORD[len(ORDER)]}",
    f"decides {pct(min(judge_share))} to {pct(max(judge_share))} of milestones, depending on the model",
    f"a second judge from another family moves each model's coverage by at most {f3(max(abs(d) for d in SWAP_DIFF))}",
    "8-bit precision or higher",
    f"an endpoint matched on at least 20 templates differs from its model's mean by at most {f3(prov_max)}",
    f"The further experiments ran on a {n_sub}-instance subset",
    f"the interpreter is offered to the {WORD[len(tool)]} closed models, not required",
    f"The paraphrase test covers the {q5_templates} templates that can be reworded without changing the problem",
]
appendix_phrases["conclusion"] = [
    f"We introduced \\ourdataset, {N_TEMPLATES} expert-certified symbolic templates in {WORD[len(BRANCH)]} engineering branches",
    f"On the final answer the {WORD[len(TOP5)]} strongest of {WORD[len(ORDER)]} LLMs are near the ceiling of these problems",
    f"Reasoning at medium effort lifts {series([tt(k) for k in RERUN])}, which return no reasoning tokens by default, by "
    f"{f2(min(m_change.values()))} to {f2(max(m_change.values()))}",
]

# Branch and domain: shaded branch rows, each followed by its domains; the models across, in order of FAC.
dom_by_branch = {b: sorted(d for d in domain_templates if branch_of[d] == b) for b in BRANCH}
lower_held = {(k, pp["a"] if pp["diff"] < 0 else pp["b"]) for k, pp in branch_pairs}  # the lower branch of the pair that holds
rows = []
for b, bname in BRANCH.items():
    rows.append(f"\\rowcolor{{gray!10}}\\textit{{{bname}}} ({n_branch_templates}) & "
                + " & ".join(f3(BL[k]["branch"][b]["mean"]) + ("$^{\\P}$" if (k, b) in lower_held else "") for k in ORDER) + " \\\\")
    rows += [f"\\quad {d.replace('_', ' ').capitalize()} ({domain_templates[d]}) & "
             + " & ".join(bold(f3(REP[k]["domain"][d]), lowest_domain[k] == d) for k in ORDER) for d in dom_by_branch[b]]
rows += ["\\midrule", "Detectable branch difference & " + " & ".join(f3(BL[k]["detectable_branch"]) for k in ORDER)]
_bp = branch_pairs[0][1] if branch_pairs else {"a": "civil_engineering", "b": "electrical_engineering", "diff": -1}
_bp_model = branch_pairs[0][0] if branch_pairs else "gpt-oss-20b"
low_b, high_b = (_bp["a"], _bp["b"]) if _bp["diff"] < 0 else (_bp["b"], _bp["a"])
put("tab:branch_domain", "branch_domain", table(
    "l " + "r " * len(ORDER), "\\textbf{Branch or domain (templates)} & " + " & ".join("\\rotatebox{90}{" + tt(k) + "}" for k in ORDER), rows,
    "\\textbf{Final Answer Accuracy by branch and by domain.} "
    f"Models in order of FAC. Shaded rows: the mean of the branch's {n_branch_templates} template means; below each, its domains' instance "
    "means, which carry no interval or test; template counts in parentheses. Bold: each model's lowest domain. $^{\\P}$The one branch "
    f"pair that differs within a model after Holm correction (Welch's $t$-test over the model's ten pairs): {tt(_bp_model)}'s "
    f"{BRANCH[high_b].lower()} above its {BRANCH[low_b].lower()}. Last row: the smallest branch difference {n_branch_templates} templates "
    "detect at 80\\% power, per model.",
    "tab:branch_domain", size="\\footnotesize", resize=True, aliases=("tab:by_branch", "tab:by_domain_kind")))

# The level means, the level gap with its tests, the share scored 0 by the depth of the gold derivation, and the controlled depth
# model's slope per model (its per-bin rates over readable responses repeat the share scored 0 but for the empty responses).
bins = [b for b in Q3[ORDER[0]]["by_milestone_count"] if b != "0"]
rows = []
for k in ORDER:
    q = Q2[k]
    chem = "" if REPAIRED else f" & {sgn(q['without_two_chemical']['gap'])} ({pv(q['without_two_chemical']['p_welch_holm'])})"
    rows.append(f"{tt(k)} & " + " & ".join(f3(BL[k]["level"][lv]["mean"]) for lv in LEVELS)
                + f" & {sgn(q['gap'])} ({ci(q['ci'])}) & {pv(q['p_welch_holm'])} & {pv(q['p_perm_holm'])} & {f3(q['detectable_planned'])}"
                + chem + " & "
                + " & ".join(f3(Q3[k]["by_milestone_count"][b]["wrong_rate"]) for b in bins)
                + f" & {sgn(D_MODELS[k]['slope'])} ({ci(D_MODELS[k]['ci'])}) & {pv(D_MODELS[k]['p_holm'])}")
n_gap_cols = 4 if REPAIRED else 5
slope_col = 5 + n_gap_cols + len(bins)  # the first of the two slope columns
header = (" & \\multicolumn{3}{c}{\\textbf{FAC by level}} & \\multicolumn{" + str(n_gap_cols) + "}{c}{\\textbf{Easy minus Advanced}} & \\multicolumn{"
          + str(len(bins)) + "}{c}{\\textbf{Share scored 0, by gold milestones}} & \\multicolumn{2}{c}{\\textbf{Logistic slope}} \\\\\n"
          "\\cmidrule(lr){2-4}\\cmidrule(lr){5-" + str(4 + n_gap_cols) + "}\\cmidrule(lr){" + str(5 + n_gap_cols) + "-" + str(4 + n_gap_cols + len(bins))
          + "}\\cmidrule(lr){" + str(slope_col) + "-" + str(slope_col + 1) + "}\n"
          "\\textbf{Model} & \\textbf{Easy} & \\textbf{Intermediate} & \\textbf{Advanced} & \\textbf{Gap (95\\% interval)} & \\textbf{Welch $p$} & "
          + mk("Permutation", "$p$") + " & \\textbf{Detectable}" + ("" if REPAIRED else " & " + mk("Without two", "chemical")) + " & "
          + " & ".join(f"\\textbf{{{b.replace('-', '--')}}}" for b in bins) + " & " + mk("Per milestone", "(95\\% interval)") + " & \\textbf{$p$ (Holm)}")
put("tab:level_gap", "results", table(
    "l r r r c r r r " + ("" if REPAIRED else "c ") + "r " * len(bins) + "c r", header, rows,
    note("depth_model") + "\\textbf{Difficulty level and trace depth.} "
    f"Left: FAC over the {n_easy} Easy, {n_int} Intermediate, and {n_adv} Advanced templates. Middle: the Easy mean minus the Advanced mean, "
    "with its 95\\% interval, the Holm-adjusted $p$ of Welch's $t$-test and of a permutation of level labels, and the smallest gap the "
    "design detects at 80\\% power"
    + ("" if REPAIRED else ", and the gap with its $p$ without the two Advanced chemical templates whose wording does not pin the answer")
    + ". Right: the share of instances scored 0, empty responses included, by the number of milestones in the gold trace ("
    + ", ".join(str(Q3[ORDER[0]]["by_milestone_count"][b]["items"]) for b in bins[:-1])
    + f", and {Q3[ORDER[0]]['by_milestone_count'][bins[-1]]['items']} instances). Far right: over readable responses only (empty responses and "
    "responses that state no answer are left out), the slope of a template-clustered logistic model of a wrong answer on the milestone count, "
    f"with the answer kind as a covariate, and its Holm-adjusted $p$ over the {WORD[len(ORDER)]} models; the slope holds after correction for "
    f"{WORD[len(depth_holds)]} of them.",
    "tab:level_gap", size="\\footnotesize", resize=True))

# Coverage behind the wrong answers, the judge's share, the step flags (the judged step check's columns, which Table 1 does not print),
# and the precision of the arithmetic flags. MC with the judge on wrong answers is tab:coverage_variants' "as scored, wrong" column.
flag_by_model = {r[0]: (r[6], int(r[2]) - int(r[5])) for r in md_table(flags, r"\| model \| flags drawn") if r[0] != "all"}
assert set(flag_by_model) == set(ORDER)
rows = [f"{tt(k)} & {Q3[k]['readable_wrong_with_milestones']} & {f3(Q3[k]['e3_coverage_on_readable_wrong'])} & "
        f"{f3(Q3[k]['e3_null_on_readable_wrong'])} & {f3(COV[k]['wrong_full_coverage'])} & {f3(Q3[k]['attribution_on_wrong']['e5_missing'])} & "
        f"{with_interval(Q3[k]['router_rate_on_wrong'], Q3[k]['router_wrong_ci'])} & {f3(Q3[k]['e5_judged_fraction'])} & "
        f"{f2(Q3[k]['router_steps_flagged_per_trace'])} & {with_interval(Q3[k]['router_judge_rate_on_fully_solved'], Q3[k]['router_judge_ci'])} & "
        + ("--" if MATCHED_HEAD and k in RERUN else f"{flag_by_model[k][0]} ({flag_by_model[k][1]})") for k in ORDER]
header = (" & \\multicolumn{6}{c}{\\textbf{Wrong answers}} & \\multicolumn{2}{c}{\\textbf{All responses}} & \\multicolumn{2}{c}{\\textbf{Correct "
          "answers}} \\\\\n\\cmidrule(lr){2-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\n"
          "\\textbf{Model} & \\textbf{Number} & " + mk("MC by", "matching") + " & " + mk("Chance", "floor") + " & " + mk("All", "reached")
          + " & " + mk("Milestone", "missing") + " & " + mk("Any step", "flagged") + " & " + mk("Judge-decided", "milestones") + " & "
          + mk("Steps flagged", "per response") + " & " + mk("Judged step", "flags") + " & " + mk("Arithmetic flag", "precision ($n$)"))
put("tab:coverage", "results", table(
    "l r r r r r c r r c c", header, rows,
    (note("matched") if RERUN else "") + "\\textbf{Milestone Coverage behind wrong answers, and the step flags.} "
    "Per model, the readable wrong answers with milestones: how many; MC by matching alone (with the judge: \\autoref{tab:coverage_variants}); "
    "the chance floor, the same responses matched against a sibling instance's milestones; and the share that reach every milestone. Over the "
    "answered wrong answers: the share with a milestone the judge rules missing, and the share with at least one step either check flags. Over "
    "all answered responses: the share of milestones the judge rather than the matcher decides, and the steps either check flags per response. "
    "Over correct answers: the share with at least one step the judge flags, and the precision of the arithmetic flags a domain expert read, "
    f"with their number ({flags_slip} of {flags_decided} were slips overall"
    + (f"; the expert read responses at the providers' default settings, so the column is empty for {series([tt(k) for k in RERUN])}, run "
       "here with reasoning at medium effort" if MATCHED_HEAD else "") + "). The two shares of responses with a flagged step carry 95\\% "
    "intervals over templates. The judged step check runs the judge at the provider's default settings, so its flags are judged, not verified: "
    f"in the expert study, {router_alone_tp} of the {router_alone_n} flags it raises inside correct answers are errors (precision "
    f"{f2(router_alone_precision)}), and together with the arithmetic check it finds {f2(router_recall)} of the flawed steps there at precision "
    f"{f2(router_precision)} (\\autoref{{tab:validation_results}}). The {res['milestones']['items_without']} instances without milestones are left "
    "out of the coverage columns.",
    "tab:coverage", size="\\footnotesize", resize=True, aliases=("tab:judged_steps",)))

# Paraphrase.
rows = [f"{tt(k)} & {sgn(Q5[k]['diff'])} ({ci(Q5[k]['ci'])}) & {ci(Q5[k]['ci90'])} & {pv(Q5[k]['p_holm'])} & {sgn(Q5[k]['e5']['diff'])}"
        for k in ORDER]
put("tab:paraphrase", "paraphrase", table(
    "l c c r r", head("Model", "FAC change (95\\% interval)", "90\\% interval", "$p$ (Holm)", "MC change"), rows,
    "\\textbf{Change under paraphrase.} "
    f"Paraphrase minus original on the {q5_pairs} expert-kept pairs over {q5_templates} templates, positive when the paraphrase scores higher: "
    f"the change in FAC with its 95\\% interval, its 90\\% interval read against the $\\pm {margin:.2f}$ margin, the Holm-adjusted $p$ of the "
    "sign-flip test over templates, and the change in MC.",
    "tab:paraphrase", size="\\footnotesize"))

# The four experiments, grouped by experiment.
ARM_GROUP = {"reasoning-medium": "Reasoning at medium effort", "openbook2": "Open book: governing equations supplied",
             "tool": "Open tool: Python interpreter offered",
             "flagship-reasoning-medium": f"{tt('gpt-5.4')} with reasoning at medium effort, against its default"}
rows, last = [], None
for a in ARMS:
    if a["arm"] != last:
        rows.append(group_row(ARM_GROUP.get(a["arm"], a["arm"]), 10))
        last = a["arm"]
    rows.append(f"{tt(a['model'])} & {a['items']} & {f3(a['main_score_on_items'])} & {f3(a['arm_score'])} & {sgn(a['diff'])} ({ci(a['ci'])}) & "
                f"{pv(a['p_holm'])} & {f3(a['detectable'])} & {ci(a['ci90'])} & {sgn(a['e5']['diff'])} & "
                f"{f3(a['digit_flag_rate_fully_solved']['main'])} / {f3(a['digit_flag_rate_fully_solved']['arm'])}")
put("tab:experiments", "experiments", table(
    "l r r r c r r c r c",
    "\\textbf{Model} & \\textbf{Instances} & " + mk("FAC,", "evaluation") + " & " + mk("FAC,", "experiment") + " & "
    + mk("Change", "(95\\% interval)") + " & \\textbf{$p$ (Holm)} & \\textbf{Detectable} & \\textbf{90\\% interval} & " + mk("MC", "change")
    + " & " + mk("Arithmetic flags,", "evaluation / experiment"), rows,
    f"\\textbf{{The {WORD[n_experiments]} experiments against the evaluation.}} "
    "Per model and experiment: the instances covered; FAC in the evaluation and under the experiment on those instances, paired by instance; "
    "the change with its 95\\% interval, the Holm-adjusted $p$ of the sign-flip test over templates, the smallest change the experiment "
    f"detects at 80\\% power, and the 90\\% interval read against the $\\pm {margin:.2f}$ margin; the change in MC; and the share of "
    "correct-answer responses with an arithmetic flag in the evaluation and under the experiment (unpaired)."
    + (f" The experiments and their evaluation column are at the providers' default settings on the same instances (for "
       f"{series([tt(k) for k in RERUN if k in {a['model'] for a in ARMS}])} without reasoning); the full-set change with reasoning on is in "
       "\\autoref{tab:matched}." if MATCHED_HEAD else ""),
    "tab:experiments", size="\\footnotesize", resize=True))

# Error analysis: the readings by model and by level in one table, with the no-error readings removed from the level shares as well.
cats = [(full, short) for full, short in CATEGORIES if full != INCOMPLETE]


def level_cell(full: str, lv: str) -> str:
    n = by_level[lv].get(full, 0)
    s = pct(n / level_n[lv])
    wo = "--" if full == NOERR else pct(n / level_n_wo[lv]) if level_n_wo[lv] else "--"
    return f"{n} ({s} / {wo})"


rows = [f"{short} & " + " & ".join(f"{by_model[m].get(full, 0)} ({maj[m].get(full, 0)})" for m in B2_TEXT) + " & "
        + " & ".join(level_cell(full, lv) for lv in LEVELS) for full, short in cats]
sample_level = {lv: sum(ERR_COMP["by_model_level"][m][lv] for m in B2_MODELS) for lv in LEVELS}
rows += ["\\midrule",
         "Wrong answers in the sample, by level & " + " & ".join(", ".join(str(ERR_COMP["by_model_level"][m][lv]) for lv in LEVELS) for m in B2_TEXT)
         + " & " + " & ".join(str(sample_level[lv]) for lv in LEVELS),
         "Readings & " + " & ".join(str(sum(by_model[m].values())) for m in B2_TEXT) + " & " + " & ".join(str(level_n[lv]) for lv in LEVELS),
         "Fleiss' $\\kappa$ & " + " & ".join(f3(ERR["fleiss_by_model"][m]) for m in B2_TEXT) + " & & & "]
header = (" & \\multicolumn{" + str(len(B2_TEXT)) + "}{c}{\\textbf{By model: readings (majority labels)}} & \\multicolumn{3}{c}{\\textbf{By "
          "level: readings (share / share without no-error readings)}} \\\\\n\\cmidrule(lr){2-" + str(1 + len(B2_TEXT)) + "}\\cmidrule(lr){"
          + str(2 + len(B2_TEXT)) + "-" + str(4 + len(B2_TEXT)) + "}\n\\textbf{Category} & " + " & ".join(tt(m) for m in B2_TEXT)
          + " & \\textbf{Easy} & \\textbf{Intermediate} & \\textbf{Advanced}")
put("tab:errors", "errors", table(
    "l " + "r " * (len(B2_TEXT) + 3), header, rows,
    f"\\textbf{{The domain experts' readings of the wrong answers.}} The sample holds {items_total} wrong answers of {WORD[len(B2_MODELS)]} models "
    f"({series([f'{items_by_model[m]} from {tt(m)}' for m in B2_TEXT])}); {WORD[readers]} domain experts of the answer's branch read each "
    f"({readings_total} readings). Left: readings per category and, in parentheses, the wrong answers whose "
    "majority label is the category. Right: readings per category by level, with their share of the level's readings and their share once the "
    "no-error readings are removed, so the error mix can be read apart from the share of answers the experts found correct. Below: the sample "
    f"by level, the readings per column, and Fleiss' $\\kappa$ over the {WORD[readers]} experts per model ({f3(fleiss_all)} overall). Categories run "
    "from the most to the least fundamental; no domain expert used the incomplete option.",
    "tab:errors", size="\\footnotesize\n\\setlength{\\tabcolsep}{4pt}", resize=True, aliases=("tab:error_by_level",)))

# The placed figures (drawn by paper_figures.py). The figures label gpt-oss-20b "GPT OSS 20B" (D9, owner 2026-10-08); each caption says so.
FIG_NAME = {**NAME, "gpt-oss-20b": "GPT OSS 20B"}  # the labels the figures print
FIG_NAME_NOTE = " ".join(f"In the figure, {tt(k)} is labeled {FIG_NAME[k]}." for k in REPRESENTATIVE if FIG_NAME[k] != NAME[k])


def fig_label(k: str) -> str:
    """A model's name with the label the figure prints, where the two differ (D9)."""
    return tt(k) + (f" ({FIG_NAME[k]} in the figure)" if FIG_NAME[k] != NAME[k] else "")


REP_NAMES = (f"{tt(top_rep[0])} (first on FAC), {tt(top_rep[1])} (first on MC), {fig_label(low_rep[0])}"
             + (f" ({'closed' if low_rep[0] in CLOSED else 'open-weights'}, {ORDINAL[ORDER.index(low_rep[0])]} on FAC)" if MATCHED_HEAD else "")
             + f", and {fig_label(low_rep[1])}")
assert all(FIG_NAME[k] == NAME[k] for k in REPRESENTATIVE if k not in low_rep), "a figure label differs for a model REP_NAMES does not annotate"
put("fig:level_bars", "main", figure(
    "level-bars.pdf",
    "\\textbf{Final Answer Accuracy by difficulty level for four representative models.} "
    f"The mean over the {n_easy} Easy, {n_int} Intermediate, and {n_adv} Advanced templates, with 95\\% intervals, for {REP_NAMES}. "
    "\\autoref{tab:level_gap} gives every model's gap with its tests.",
    "fig:level_bars"))
put("fig:domain_radar", "branch_domain", figure(
    "domain_radar/domain-radar-labeled.pdf",
    "\\textbf{Final Answer Accuracy by domain for four representative models.} "
    f"The instance mean over each of the {len(domain_templates)} domains, in the order of their branches around the rim, for the same four "
    "models; the radial axis starts at 0.4. The domain means carry no interval and no test (\\autoref{tab:by_domain_kind}). " + FIG_NAME_NOTE,
    "fig:domain_radar"))
assert len(BRANCH) == 5 and n_branch_templates * len(BRANCH) == N_TEMPLATES  # "a fifth" in the caption
put("fig:branch_bars", "branch_domain", figure(
    "branch-bars.pdf",
    "\\textbf{Final Answer Accuracy by engineering branch for four representative models.} "
    f"One stacked bar per model, for {REP_NAMES}. Each of the {WORD[len(BRANCH)]} branches holds {n_branch_templates} of the {N_TEMPLATES} "
    "templates, so a segment is that branch's share of the model's FAC: the number inside it is the branch's FAC, and its height is a fifth of "
    "that. The segments add up to the model's FAC, printed above the bar; color marks the model and hatching the branch.",
    "fig:branch_bars"))
put("fig:error_categories", "main", figure(
    "error-categories.pdf",
    "\\textbf{Error categories of the wrong answers read.} "
    f"Each bar gives one model's readings ({rng([sum(by_model[m].values()) for m in B2_MODELS], str)}; {readers} domain experts per wrong "
    "answer, at the providers' default settings) by category; \\autoref{tab:error_by_level} gives them by level. " + FIG_NAME_NOTE,
    "fig:error_categories"))
# The blocks the registry names for this script. The judged step flags are in tab:coverage (also labelled tab:judged_steps), the depth
# slopes in tab:level_gap, and the single-path counts and the carried-precision split are sentences of the results appendix.
REGISTRY = ("tab:main_results", "tab:matched", "tab:coverage_variants", "tab:scoring_variants", "tab:providers", "tab:coverage",
            "tab:level_gap", "tab:errors")
assert all(label in blocks for label in REGISTRY)



# ----------------------------------------------------------------------------------------------- write and check
def figure_data() -> dict:
    """The data the figures draw from, in the shape paper_figures.draw_all documents; every value is read above from the result files."""
    return {"order": ORDER, "q1": Q1, "q2": Q2, "q3": Q3, "q5": Q5, "bl": BL, "rep": REP, "by_model": by_model, "categories": CATEGORIES,
            "incomplete": INCOMPLETE, "b2_models": B2_MODELS, "representative": REPRESENTATIVE, "branch": BRANCH, "levels": LEVELS,
            "n_templates": N_TEMPLATES, "margin": margin, "fig_name": FIG_NAME, "name": NAME,
            "branch_of": branch_of}


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


def marker_re(name: str) -> re.Pattern:
    return re.compile(r"% BEGIN GENERATED " + re.escape(name) + r" .*?% END GENERATED " + re.escape(name), re.S)


def tex_files(tree: Path) -> list[Path]:
    return sorted(p for p in tree.rglob("*.tex") if "figs" not in p.parts)


def locate(tree: Path) -> dict[str, Path]:
    """The file in the tree that holds each block's markers (a block is written wherever its markers are)."""
    loc = {}
    for p in tex_files(tree):
        t = p.read_text(encoding="utf-8")
        for label in blocks:
            if f"% BEGIN GENERATED {label} " in t:
                assert label not in loc, f"the markers of {label} appear in {loc[label]} and in {p}"
                loc[label] = p
    return loc


def phrase_text() -> str:
    out = ["% phrases 6_results.tex must contain", *phrases]
    for k, ph in appendix_phrases.items():
        out += ["", f"% phrases {FILES[k]} must contain", *ph]
    return "\n".join(out) + "\n"


def write_blocks(tree: Path, generated_dir: Path | None, draw: bool) -> None:
    """Write every block into the file of the tree that holds its markers. A block without markers goes to generated_dir/<label>.tex
    when a generated_dir is given (Phase 1, --out) and stops the run otherwise (Phase 2: WS-D1 and D2 place the markers first)."""
    loc = locate(tree)
    by_file = defaultdict(list)
    for label in blocks:
        if label in loc:
            by_file[loc[label]].append(label)
        elif generated_dir is not None:
            generated_dir.mkdir(parents=True, exist_ok=True)
            (generated_dir / f"{label.replace(':', '_')}.tex").write_text(block(label, blocks[label][1]) + "\n", encoding="utf-8")
            print(f"{label}: no markers in the tree; written to {generated_dir.name}/{label.replace(':', '_')}.tex (default file {FILES[blocks[label][0]]})")
        else:
            raise SystemExit(f"no file of {tree} has markers for {label}; place them (default file {FILES[blocks[label][0]]}) or use --out DIR")
    for path, labels in by_file.items():
        tex = path.read_text(encoding="utf-8")
        for label in labels:
            tex = marker_re(label).sub(lambda m, label=label: block(label, blocks[label][1]), tex, count=1)
        path.write_text(tex, encoding="utf-8")
        print(f"{path.relative_to(tree)}: {len(labels)} blocks written ({', '.join(labels)})")
    if generated_dir is not None:
        generated_dir.mkdir(parents=True, exist_ok=True)
        (generated_dir / "phrases.txt").write_text(phrase_text(), encoding="utf-8")
        print(f"phrases written to {generated_dir.name}/phrases.txt")
    if draw:
        drawn = paper_figures.draw_all(figure_data(), tree / "figs")
        paper_figures.fig_level_gap(ORDER, Q2, FIG_NAME, tree / "figs")  # the figures not placed in the paper, on the same data
        paper_figures.fig_coverage_wrong(ORDER, Q3, FIG_NAME, tree / "figs")
        paper_figures.fig_paraphrase(ORDER, Q5, margin, FIG_NAME, tree / "figs", res["q5"]["vs_repeats"])
        print(f"{len(drawn)} figures and level-gap.pdf, coverage-wrong.pdf, paraphrase.pdf drawn into {tree / 'figs'}")


def report_sources() -> int:
    """Stand-ins, quick-mode files, source drift and failing claims, printed; drift and claims count as failures."""
    for name in sorted(STAND_IN):
        print(f"STAND-IN read: results/stand_in/{name}.json (results/{name}.json does not exist yet)")
    for name in sorted(QUICK):
        print(f"QUICK-mode file read: results/{name}.json")
    for d in DRIFT:
        print(f"SOURCE DRIFT: {d}")
    for c in CLAIMS:
        print(f"CLAIM FAILS: {c}")
    return len(DRIFT) + len(CLAIMS)


def check(tree: Path, pending_ok: bool) -> int:
    """Check the tree. pending_ok (Phase 1, with --out): a block without markers and a phrase or number the prose does not hold yet are
    listed as pending and do not fail the check; everything else fails as before."""
    bib_keys = set(re.findall(r"@\w+\{([^,\s]+),", (tree / BIB).read_text(encoding="utf-8")))
    labels_defined = set()
    for p in tex_files(tree):
        labels_defined |= set(re.findall(r"\\label\{([^}]*)\}", p.read_text(encoding="utf-8")))
    own = {"main": phrases, **appendix_phrases}
    loc = locate(tree)
    failures = pending = 0
    for label, (file_key, body) in blocks.items():
        if label not in loc:
            print(f"{'PENDING' if pending_ok else 'UNPLACED'} block (no markers in the tree; default file {FILES[file_key]}): {label}")
            pending += pending_ok
            failures += not pending_ok
            continue
        m = marker_re(label).search(loc[label].read_text(encoding="utf-8"))
        if flat(m.group(0)) != flat(block(label, body)):
            print(f"STALE block in {loc[label].relative_to(tree)}: {label}")
            failures += 1
        elif "generated" in loc[label].relative_to(tree).parts:
            print(f"{'PENDING' if pending_ok else 'UNPLACED'} block (current, but in {loc[label].parent.name}/, not yet placed in the prose; "
                  f"default file {FILES[file_key]}): {label}")
            pending += pending_ok
            failures += not pending_ok
    known_all = numbers(prose(" ".join(p for ph in own.values() for p in ph)))
    for key, rel in FILES.items():
        path = tree / rel
        tex = path.read_text(encoding="utf-8")
        body = flat(re.sub(r"(?<!\\)%.*", "", strip_generated(tex)))
        missing = [p for p in own[key] if flat(p) not in body]
        known = numbers(prose(" ".join(own[key])))
        words = {w.lower() for w in NUMBER_WORDS.findall(prose(" ".join(own[key])))}
        stray = sorted(numbers(prose(tex)) - known) + sorted({w.lower() for w in NUMBER_WORDS.findall(prose(tex))} - words)
        cited = {k.strip() for c in re.findall(r"\\cite[pt]\{([^}]*)\}", tex) for k in c.split(",")}
        unresolved = sorted(cited - bib_keys)
        refs = set(re.findall(r"\\autoref\{([^}]*)\}", tex))
        undefined = sorted(refs - labels_defined)
        for p in missing:
            print(f"{'PENDING phrase' if pending_ok else 'MISSING'} in {path.name}: {p}")
        for n in stray:
            print(f"{'PENDING number' if pending_ok else 'NOT GENERATED'} in {path.name}: {n}")
        for k in unresolved:
            print(f"UNRESOLVED citation in {path.name}: {k}")
        for r in undefined:
            print(f"UNDEFINED label in {path.name}: {r}")
        long_lines = [i + 1 for i, l in enumerate(strip_generated(tex).split("\n")) if len(l) > 100]
        if long_lines:
            print(f"LINES over 100 characters in {path.name}: {long_lines}")
        placed_here = [label for label, p in loc.items() if p == path]
        print(f"{path.name}: {len(placed_here)} generated blocks placed; {len(own[key]) - len(missing)} of {len(own[key])} phrases present; "
              f"{len(stray)} numbers not generated; {len(cited) - len(unresolved)} of {len(cited)} citation keys resolve; "
              f"{len(refs) - len(undefined)} of {len(refs)} references defined")
        if pending_ok:
            pending += len(missing) + len(stray)
        else:
            failures += len(missing) + len(stray)
        failures += len(unresolved) + len(undefined) + len(long_lines)
    for name in paper_figures.FIGURES:
        if not (tree / "figs" / name).exists():
            print(f"MISSING figure {tree / 'figs' / name}")
            failures += 1
    failures += report_sources()
    print(f"check of {tree}: {failures} failures" + (f", {pending} items pending Phase 2 (markers and prose)" if pending_ok else ""))
    return failures


def flags_line() -> str:
    return f"--headline {HEADLINE}" + (" --repaired" if REPAIRED else "") + (" --judged-in-table" if JUDGED_IN_TABLE else "")


if __name__ == "__main__":
    if ARGS.check:
        tree = Path(ARGS.out).resolve() if ARGS.out else SRC
        if not (tree / "6_results.tex").exists():
            raise SystemExit(f"{tree} holds no tex tree (run --out {tree} first)")
        sys.exit(1 if check(tree, pending_ok=bool(ARGS.out)) else 0)
    elif ARGS.out or ARGS.write:
        if ARGS.out:
            tree = Path(ARGS.out).resolve()
            if tree == SRC or SRC in tree.parents:
                raise SystemExit("--out must point outside the real tex tree")
            shutil.copytree(SRC, tree, dirs_exist_ok=True)
            print(f"tex tree copied to {tree} ({flags_line()})")
            write_blocks(tree, tree / "generated", draw=not ARGS.text_only)
        else:
            write_blocks(SRC, None, draw=not ARGS.text_only)
        print(f"wrote {len(blocks)} generated blocks ({flags_line()})" + ("" if ARGS.text_only else f" and {len(paper_figures.FIGURES)} figures"))
        report_sources()
    else:
        print(phrase_text(), end="")
        report_sources()
