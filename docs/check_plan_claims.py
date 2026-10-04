"""Check every number the paper plan claims against the file that prints it.

    python docs/check_plan_claims.py

Each check names the claim in `docs/PAPER_PLAN_OCT2026.md` (C1 to C21, plus the figures quoted in its section
notes), the source file, and either a regular expression that must match the source or a value recomputed from a
table in the source. A FAIL means the plan's number does not match the record and the plan must change; the
script never edits anything. Exit status is the number of failures.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "full_run_28092026"
PILOT = ROOT / "evaluator_pilot_17092026"
CERT = ROOT / "template_annotation_23092026"
DECISIONS = ROOT / "docs/re-implementation-sep/DECISIONS.md"

results: list[tuple[str, str, bool, str]] = []


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def check(claim: str, what: str, ok: bool, observed: str) -> None:
    results.append((claim, what, bool(ok), observed))


def grep(claim: str, what: str, path: Path, pattern: str) -> None:
    """A plain space in the pattern matches any run of whitespace, so wrapped lines in the record still match."""
    pattern = re.sub(r"(?<!\\) ", r"\\s+", pattern)
    m = re.search(pattern, read(path), re.S)
    check(claim, what, m is not None, (m.group(0)[:90] if m else f"no match for /{pattern}/ in {path.name}"))


def table(path: Path, header_pattern: str, nth: int = 0) -> list[list[str]]:
    """Rows (cells stripped) of the nth table whose header line matches the pattern."""
    lines = read(path).split("\n")
    hits = [i for i, l in enumerate(lines) if l.startswith("|") and re.search(header_pattern, l)]
    if len(hits) <= nth:
        raise LookupError(f"table /{header_pattern}/ #{nth} not found in {path.name}")
    i = hits[nth] + 2  # skip header and separator
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return rows


def num(cell: str) -> float:
    m = re.search(r"[-+]?\d+(?:\.\d+)?", cell.replace("−", "-"))
    if not m:
        raise ValueError(cell)
    return float(m.group(0))


def rng(claim: str, what: str, values: list[float], lo: float, hi: float, nd: int = 3) -> None:
    observed = (round(min(values), nd), round(max(values), nd))
    check(claim, what, observed == (lo, hi), f"observed {observed}, plan says ({lo}, {hi})")


# ---------------------------------------------------------------- C1: the benchmark and the evaluation set
rows = [json.loads(l) for l in read(RUN / "manifest.jsonl").splitlines() if l.strip()]
t_by_branch = defaultdict(set)
t_by_level = defaultdict(set)
domains, areas = set(), set()
items_by_level, items_by_branch = Counter(), Counter()
for r in rows:
    t_by_branch[r["branch"]].add(r["template_id"])
    t_by_level[r["level"]].add(r["template_id"])
    domains.add((r["branch"], r["domain"]))
    areas.add((r["branch"], r["domain"], r["area"]))
    items_by_level[r["level"]] += 1
    items_by_branch[r["branch"]] += 1
templates = {t for s in t_by_branch.values() for t in s}
check("C1", "150 templates", len(templates) == 150, str(len(templates)))
check("C1", "30 templates per branch", all(len(s) == 30 for s in t_by_branch.values()), str({b: len(s) for b, s in sorted(t_by_branch.items())}))
check("C1", "15 domains, 42 areas", (len(domains), len(areas)) == (15, 42), f"{len(domains)} domains, {len(areas)} areas")
lv = {k: len(v) for k, v in t_by_level.items()}
check("C1", "58 / 58 / 34 templates by level", (lv.get("Easy"), lv.get("Intermediate"), lv.get("Advanced")) == (58, 58, 34), str(lv))
check("C1", "2,250 instances; 870 / 870 / 510; 450 per branch", len(rows) == 2250 and (items_by_level["Easy"], items_by_level["Intermediate"], items_by_level["Advanced"]) == (870, 870, 510) and all(v == 450 for v in items_by_branch.values()), f"{len(rows)}; {dict(items_by_level)}; {dict(items_by_branch)}")
grep("C1", "all distinct questions", RUN / "README.md", r"distinct questions \| 2,250; no template needs a repeat")

# ---------------------------------------------------------------- C2: certification
grep("C2", "integrity checks: 150 of 150 pass", CERT / "layer0/gate_report.md", r"\*\*pass the gate\*\* \| \*\*150\*\*")
grep("C2", "screen pass 2: 147 / 3 / 0", CERT / "README.md", r"147 pass, 3 controversial, 0 critical")
grep("C2", "screen agreement AC1 0.93", CERT / "README.md", r"AC1 0\.93")
grep("C2", "34 of 45 screen claims confirmed", DECISIONS, r"34 of 45 were real")
grep("C2", "60 of 60 planted defects", CERT / "layer2/RESULTS.md", r"Overall: 60 of 60 planted defects rejected \(100%\)")
grep("C2", "459 of 504 hand checks within 1%", CERT / "layer2/RESULTS.md", r"504 hand checks with a comparable number: 459 matched the template within 1%")
grep("C2", "44 of 45 mismatches ended in a rejection", CERT / "layer2/RESULTS.md", r"Of the 45 mismatches, 44 ended in a rejection")
grep("C2", "AC1 0.914, Fleiss 0.589, AC2 0.947 / 0.980 / 0.954", CERT / "layer2/RESULTS.md", r"\| all \| 150 \| 0\.589 \| 0\.914 \| 89% \| 0\.947 \| 0\.980 \| 0\.954 \|")
grep("C2", "screen false-positive rate 10.2%", CERT / "layer2/RESULTS.md", r"15 of 147 \(false-positive rate 10\.2%\)")
grep("C2", "53 rejections (43 + 10)", DECISIONS, r"53 rejections of real templates, 43 in round 1 and 10 in round 2")
grep("C2", "150 of 150 certified", CERT / "layer2/CERTIFICATION.md", r"\| certified \| 150 of 150 \|")
cert = read(CERT / "layer2/CERTIFICATION.md")
branch_of = dict(re.findall(r"^\| `(template_\w+)` \| (\w+) \|", cert, re.M))
rejected = re.findall(r"^- \*\*(template_\w+)\*\* rejected by", read(CERT / "layer2/RESULTS.md"), re.M)
by_branch = Counter(branch_of.get(t, "?") for t in rejected)
check("C2", "22 round-1 rejections, all in the original three branches (chem 9, elec 3, mech 10)",
      len(rejected) == 22 and dict(by_branch) == {"chemical": 9, "electrical": 3, "mechanical": 10}, f"{len(rejected)}: {dict(by_branch)}")

# ---------------------------------------------------------------- S3: the further numbers of the benchmark section
t_by_domain, areas_by_domain, t_by_kind = defaultdict(set), defaultdict(set), defaultdict(set)
for r in rows:
    t_by_domain[(r["branch"], r["domain"])].add(r["template_id"])
    areas_by_domain[(r["branch"], r["domain"])].add(r["area"])
    t_by_kind[r["answer_type"]].add(r["template_id"])
per_domain = {k: len(v) for k, v in t_by_domain.items()}
chem = sorted(n for (b, _), n in per_domain.items() if b == "chemical_engineering")
check("S3", "ten templates per domain, 8 to 12 in chemical engineering", chem == [8, 10, 12] and all(n == 10 for (b, _), n in per_domain.items() if b != "chemical_engineering"), str(sorted(per_domain.values())))
n_areas = [len(v) for v in areas_by_domain.values()]
check("S3", "three domains per branch, two to four areas per domain", all(sum(1 for (b, _) in per_domain if b == br) == 3 for br in t_by_branch) and (min(n_areas), max(n_areas)) == (2, 4), f"areas per domain {min(n_areas)} to {max(n_areas)}")
check("S3", "six answer kinds", len(t_by_kind) == 6, str({k: len(v) for k, v in sorted(t_by_kind.items())}))
div = json.loads(read(RUN / "diversity.json"))
multi = {t["template_id"] for t in div if t["reach"]["paths_lower"] >= 2}
multi_in_set = {t["template_id"] for t in div if t["pool"]["paths_lower"] >= 2}
check("S3", "92 templates with more than one reasoning path, 58 with one (lower reading)", (len(multi), len(div) - len(multi)) == (92, 58), f"{len(multi)} / {len(div) - len(multi)} in 500 draws")
check("S3", "each of the 92 contributes at least two paths to the evaluation set", multi <= multi_in_set and len(multi_in_set) == 92, f"{len(multi & multi_in_set)} of {len(multi)}; {len(multi_in_set)} in the evaluation set")
freeze = json.loads(read(RUN / "FREEZE.json"))
check("S3", "15 instances per template from a 128-bit seed", (freeze["instances_per_template"], freeze["seed"]["bits"]) == (15, 128), f"{freeze['instances_per_template']}, {freeze['seed']['bits']} bits")
grep("S3", "selection from each template's first 100 draws", RUN / "FREEZE.json", r"instance indices 0 to 99")
grep("S3", "54 templates edited to pass the integrity checks", CERT / "README.md", r"54 templates edited over three closure rounds")
grep("S3", "24 templates flagged in the first screen pass (126 / 13 / 11)", CERT / "README.md", r"pass 1 \(2026-09-23\) 126 pass, 13 controversial, 11 critical failure, AC1 on the flag 0\.84; the 24 flags verified and fixed")
grep("S3", "three screen judges, three instances each", CERT / "screen/pass2/stats.md", r"Judges: grok-4\.6, minimax-m3, mimo-v2\.5-pro\..*instance seeds \[1001, 1002, 1003\]")
plants = re.findall(r"^\| plant_(\w{3})_\d \(template_\w+\) \| (\w+):", read(CERT / "layer2/RESULTS.md"), re.M)
plant_kinds = defaultdict(set)
for br, kind in plants:
    plant_kinds[br].add(kind)
check("S3", "20 planted defects, four per branch, one of each kind", len(plants) == 20 and len(plant_kinds) == 5 and all(v == {"constant", "unit", "sign", "arithmetic"} for v in plant_kinds.values()), f"{len(plants)}; {sorted(plant_kinds)}")
grep("S3", "15 experts", CERT / "layer2/RESULTS.md", r"from 510 label rows by 15 experts")
grep("S3", "experts told planted defects exist, not which", CERT / "layer2/README.md", r"Experts are told quality-control items exist, not which")
grep("S3", "no screen verdict shown, no discussion before submission", CERT / "layer2/README.md", r"no screen verdict shown, no discussion until all three have submitted")
grep("S3", "five instances shown per template", CERT / "layer2/CERTIFICATION.md", r"the five instances they were shown")
grep("S3", "round 1: 20 templates revised, 2 objections not adopted", CERT / "layer2/fixes_round1.md", r"Twenty templates were changed; two claims were not adopted")
grep("S3", "round 2: rejections on 5 of the 22 templates (17 + 2 + 3)", CERT / "layer2/RESULTS_round2.md", r"Templates: 17 approved by all three, 2 approved by majority, 3 rejected by majority")
cert_rounds = table(CERT / "layer2/CERTIFICATION.md", r"Round \| Labels")
check("S3", "rounds review 150, 22, 5 and 2 templates", [int(r[3]) for r in cert_rounds] == [150, 22, 5, 2], str([r[3] for r in cert_rounds]))

# ---------------------------------------------------------------- C3: validation against the experts
sv = RUN / "SCORER_VALIDATION.md"
grep("C3", "answer check 0.982 non-partial (current code)", sv, r"\| answer, non-partial agreement \| 0\.947 \| 0\.947 \| yes \| 0\.982 \|")
grep("C3", "answer check 0.927 three-way (current code)", sv, r"\| answer, three-way agreement \| 0\.893 \| 0\.893 \| yes \| 0\.927 \|")
grep("C3", "milestone F1 0.923 without the judge", sv, r"\| E3 F1 \| 0\.921 \| 0\.921 \| yes \| 0\.923 \|")
grep("C3", "milestone F1 0.958 with the judge", PILOT / "PILOT_SUMMARY.md", r"\*\*0\.930\*\* \| \*\*0\.989\*\* \| \*\*0\.958\*\*")
grep("C3", "judge credited 0 of 88 fabricated values", PILOT / "RESULTS_E5.md", r"should be MISSING \| \*\*0\*\* \| 66 \| \*\*22")
grep("C3", "arithmetic check precision 0.817, recall 0.427 (hard case)", sv, r"hard case, precision \| 0\.750 \| 0\.750 \| yes \| 0\.817 \|.*hard case, recall \| 0\.320 \| 0\.320 \| yes \| 0\.427 \|")
grep("C3", "panel design's answer agreement 0.747 (non-partial)", PILOT / "RESULTS_X1.md", r"\| E0, the published framework \| 0\.747 \|")

# ---------------------------------------------------------------- C4: planted defects and the limit of verification
ps = PILOT / "PILOT_SUMMARY.md"
grep("C4", "deterministic checks 0 of 60 conceptual; GPT-5 20 of 60; Opus 8 of 60; MiMo 16 of 51", ps, r"\*\*0 of 60\*\* \| 20 of 60 \(0\.333\) \| 8 of 60 \(0\.133\) \| 16 of 51 \(0\.314\)")
grep("C4", "MiMo 16 of 52 after the recount", RUN / "RESULTS_PAPER_NOTES.md", r"MiMo's 16 of 52")
grep("C4", "Grok 22 of 60 conceptual, no false alarm", DECISIONS, r"22 of 60 conceptual defects caught \(0\.367")
grep("C4", "answer predicts the experts' verdict, AUROC 0.974", ps, r"AUROC 0\.974")
grep("C4", "93 of 228 correct-answer traces carry an incorrect step (178 steps)", ps, r"93 of the 228 correct-answer traces \(178 steps\)")
grep("C4", "175 of the 178 are arithmetic", ps, r"175 of the 178 flawed steps behind a correct answer are calculation slips")
grep("C4", "the design detects AUROC differences of 0.12 to 0.19", ps, r"0\.12 to 0\.19")
grep("C4", "best PRM: three of four flags false inside correct-answer traces (0.246)", ps, r"falls to \*\*0\.246\*\*")
grep("C4", "experts: kappa 0.781 / 0.828 / 0.880 / 0.966", ps, r"\*\*0\.781\*\*.*\*\*0\.828\*\*.*0\.880.*0\.966")
grep("C4", "1,042 incorrect steps with a reason; 272 split steps adjudicated", ps, r"All 1,042 of them.*272\s+split steps across 134 traces")

# ---------------------------------------------------------------- C5 to C14: the main results
RES = RUN / "results/RESULTS.md"
q1 = table(RES, r"^\| model \| answer score \| 95% CI \| SD within")
fac = {r[0].strip("`"): num(r[1]) for r in q1}
check("C5", "FAC 0.814 (gpt-oss-20b) to 0.976 (DeepSeek V4.1 Flash)", (fac["gpt-oss-20b"], fac["deepseek-v4.1-flash"]) == (0.814, 0.976) and min(fac.values()) == 0.814 and max(fac.values()) == 0.976, str(fac))
top5 = sorted(fac.values(), reverse=True)[:5]
check("C5", "top five within 0.009", round(top5[0] - top5[-1], 3) == 0.009, f"{top5}")
grep("C5", "32 of 55 pairs differ after Holm", RES, r"Of the 55 pairs, 32 differ at a Holm-adjusted p below 0\.05 on the answer score")
grep("C5", "23 pairs below the detectable difference, all non-significant", RES, r"23 of the 55 pairs, 23 of them not significant")
unusable = sum(int(re.match(r"\d+", r[10]).group(0)) for r in q1)
check("C5/5.1", "246 responses with no readable answer (1.0%)", unusable == 246, str(unusable))

q2 = table(RES, r"^\| model \| Easy \| Advanced \| gap \| 95% CI \| Welch p")
rng("C6", "level gap +0.055 to +0.199", [num(r[3]) for r in q2], 0.055, 0.199)
grep("C6", "4 of 11 hold under Welch with Holm", RES, r"\*\*4 of 11\*\* hold, against 9 under the planned permutation")
grep("C6", "0 of 11 without the two chemical templates", RES, r"Without the two chemical templates the gap holds for 0 of 11 models")

cov = table(RES, r"^\| model \| traces with milestones \| E3, all \|")
rng("C7", "Milestone Coverage 0.816 to 0.923 (all responses)", [num(r[6]) for r in cov], 0.816, 0.923)
grep("C7", "tau 0.709 (0.514 to 0.855) between the two orderings", RES, r"0\.709 \(95% CI 0\.514 to 0\.855")
grep("C7", "28 pairs differ (sign-flip), 29 (Wilcoxon)", RES, r"28 differ at a Holm-adjusted p below 0\.05 under the sign-flip test and 29 under Wilcoxon")
verb = table(RES, r"^\| model \| rho\(coverage, steps\)")
rng("C7", "coverage against steps written, rho −0.26 to −0.11", [num(r[1]) for r in verb], -0.259, -0.113)

wrong = table(RES, r"^\| model \| wrong-answer traces \| of them unusable \|")
readable = [num(re.search(r"readable (\d\.\d+)", r[8]).group(1)) for r in wrong]
rng("C8", "coverage on wrong answers 0.464 to 0.845 (readable)", readable, 0.464, 0.845)
rng("C8", "chance floor 0.104 to 0.220", [num(r[7]) for r in wrong], 0.104, 0.22)
a11 = table(RES, r"^\| model \| fully solved traces \(items with milestones\)")
big = [num(r[5]) for r in a11 if num(r[4]) > 100]
rng("C8", "13% to 26% of wrong answers are complete derivations (models with >100 wrong answers)", big, 0.13, 0.263)

digit = table(RES, r"^\| model \| wrong answers, answered \| flag rate \|")
rng("C9", "arithmetic flags on correct-answer responses 0.005 to 0.225", [num(r[5]) for r in digit], 0.005, 0.225)
rng("C9", "calculations read per response 2.80 to 9.83", [num(r[8]) for r in digit], 2.8, 9.83, 2)
grep("C9", "flag precision 0.905 on this roster (171 of 189)", RES, r"171 of the 189 decided are real, 0\.905")
router = table(RES, r"^\| model \| calls \| without a reply \| flagged, fully solved \|")
rng("C9", "judged step flags 0.011 to 0.220", [num(r[5]) for r in router], 0.011, 0.22)
judge = table(RES, r"^\| model \| calls \| without a reply \| judged fraction \|")
rng("4.1", "judge decides 11% to 25% of required milestones", [num(r[3]) for r in judge], 0.112, 0.252)

depth = table(RES, r"^\| model \| 0 \(70 items\) \|")
rng("C10", "wrong-answer rate on 6+ milestones 0.037 to 0.291", [num(r[6]) for r in depth], 0.037, 0.291)
rng("C10", "wrong-answer rate on 1-milestone instances 0.000 to 0.045", [num(r[2]) for r in depth], 0.0, 0.045)

q4 = table(RES, r"^\| model \| single: all \| some \| none \|")
single = {r[0].strip("`"): num(r[1]) for r in q4}
strong = [single[m] for m in ("deepseek-v4.1-flash", "kimi-k3", "claude-sonnet-5", "glm-5.3-flash", "muse-glimmer-30b", "glm-5.3")]
weak = [single[m] for m in ("gpt-oss-20b", "gpt-5.4-mini", "gemma-4-26b-a4b", "qwen3-235b-a22b-2507")]
rng("C11", "strongest six solve all 15 instances of 85% to 91% of single-path templates", strong, 0.845, 0.914)
rng("C11", "weakest four 36% to 48%", weak, 0.362, 0.483)
check("C11", "Gemini 3.1 Flash-Lite sits between (0.672), so 'the weakest' must mean the four named", single["gemini-3.1-flash-lite"] == 0.672, str(single["gemini-3.1-flash-lite"]))

q5 = table(RES, r"^\| model \| items \| answer score diff \|")
rng("C12", "paraphrase change −0.031 to +0.020", [num(r[2]) for r in q5], -0.031, 0.02)
check("C12", "277 pairs for every model", all(num(r[1]) == 277 for r in q5), str({r[0]: r[1] for r in q5}))
check("C12", "no change significant after Holm", all(num(r[4]) > 0.05 for r in q5), str([r[4] for r in q5]))
bounds = table(RES, r"^\| model \| answer score diff \| 90% CI \| within the margin \|")
check("C12", "10 of 11 within ±5 points at 90%", sum(r[3] == "yes" for r in bounds) == 10, f"{sum(r[3] == 'yes' for r in bounds)} yes")
grep("C12", "tau 0.673 (0.455 to 0.881); noise floor median 0.782", RES, r"0\.673, 95% CI 0\.455 to 0\.881.*median tau 0\.782")
grep("C12", "316 passed the checks, 39 rejected by the experts", RES, r"316 pairs returned, 39 rejected")
grep("C12", "115 of 150 templates covered", RUN / "PARAPHRASE_PAPER_NOTES.md", r"\*\*115 of 150\*\*")

rep = table(RES, r"^\| model \| items \| repeat1 \|")
rng("C13", "repeat SD 0.004 to 0.013", [num(r[6]) for r in rep], 0.004, 0.013)
rng("C13", "same verdict in every repeat 86% to 94%", [num(r[8]) for r in rep], 0.857, 0.937)

sens = table(RES, r"^\| model \| tolerance half \| fitted \(headline\) \|")
tau_row = [r for r in sens if r[0].startswith("tau")][0]
check("C14", "ordering unchanged at half and double tolerance (tau 0.927)", (num(tau_row[1]), num(tau_row[3])) == (0.927, 0.927), f"{tau_row[1]}, {tau_row[3]}")
model_rows = [r for r in sens if r[0].startswith("`")]
sym_diff = [round(num(r[7]) - num(r[2]), 3) for r in model_rows]
rng("C14", "without the nine symbolic templates the scores move by -0.003 to +0.008", sym_diff, -0.003, 0.008)
short = table(RUN / "SHORTCUT_AUDIT.md", r"^\| model \| headline \| without the newly flagged \|")
check("C14", "without the shortcut templates the headline moves by at most 0.002", max(abs(num(r[1]) - num(r[2])) for r in short) <= 0.002 + 1e-9, str(max(abs(num(r[1]) - num(r[2])) for r in short)))

# ---------------------------------------------------------------- C15 to C18: the conditions
grep("C15", "reasoning on: GPT-5.4 mini 0.858 to 0.951, +0.093 (0.056 to 0.133)", RES, r"\| reasoning-medium \| `gpt-5\.4-mini` \| 450 \| 0\.858 \| 0\.951 \| \+0\.093 \| 0\.056 to 0\.133")
grep("C15", "reasoning on: Gemini 3.1 Flash-Lite +0.016 (−0.010 to 0.042)", RES, r"\| reasoning-medium \| `gemini-3\.1-flash-lite` \| 450 \| 0\.874 \| 0\.890 \| \+0\.016 \| -0\.010 to 0\.042")
grep("C16", "DeepSeek V4 Pro 0.968 on the subset", RES, r"\| `deepseek-v4-pro` \| flagship \| 450 \| 0\.968 \| 0\.950 to 0\.983")
grep("C16", "GPT-5.4 default 0.941, no reasoning tokens", RES, r"\| `gpt-5\.4` \| flagship \| 450 \| 0\.941 \| 0\.909 to 0\.969")
grep("C16", "GPT-5.4 with reasoning 0.980", RES, r"\| `gpt-5\.4` \| flagship-reasoning-medium \| 450 \| 0\.980 \| 0\.964 to 0\.993")
c3 = table(RES, r"^\| model \| arm \| items \| answer score \|")
main_top5 = sorted([num(r[3]) for r in c3 if r[1] == "main"], reverse=True)[:5]
rng("C16", "the top five on the same instances 0.964 to 0.976", main_top5, 0.964, 0.976)
grep("C17", "open book: gpt-oss-20b +0.064 (0.027 to 0.104) on 405 instances", RES, r"\| openbook2 \| `gpt-oss-20b` \| 405 \| 0\.816 \| 0\.880 \| \+0\.064 \| 0\.027 to 0\.104")
grep("C17", "open book: +0.043 on the 388 instances answered in both", RES, r"\| openbook2 \| `gpt-oss-20b` \|[^\n]*\+0\.043 \(388\)")
grep("C17", "open book: GPT-5.4 mini −0.019", RES, r"\| openbook2 \| `gpt-5\.4-mini` \| 405 \| 0\.851 \| 0\.832 \| -0\.019")
grep("C17", "open book: Claude Sonnet 5 +0.007", RES, r"\| openbook2 \| `claude-sonnet-5` \| 405 \| 0\.969 \| 0\.977 \| \+0\.007")
ob = [r for r in table(RES, r"^\| arm \| model \| items \| base on these items") if r[0] == "openbook2"]
rng("C17", "open book: coverage rises 0.03 to 0.06 (E5-strict change)", [num(r[18]) for r in ob], 0.029, 0.064)
grep("C18", "tool: GPT-5.4 mini −0.004 (−0.042 to 0.036)", RES, r"\| tool \| `gpt-5\.4-mini` \| 450 \| 0\.858 \| 0\.853 \| -0\.004 \| -0\.042 to 0\.036")
grep("C18", "tool: Claude Sonnet 5 +0.004 (−0.004 to 0.014)", RES, r"\| tool \| `claude-sonnet-5` \| 450 \| 0\.971 \| 0\.976 \| \+0\.004 \| -0\.004 to 0\.014")
grep("C18", "tool used on 19% (0.4 calls) and 66% (0.9 calls) of instances", RES, r"\| tool \| `gpt-5\.4-mini` \|[^\n]*\| 0\.19, 0\.4 \|.*\| tool \| `claude-sonnet-5` \|[^\n]*\| 0\.66, 0\.9 \|")
grep("C18", "GPT-5.4 mini's arithmetic flags 0.127 to 0.093 under the tool", RES, r"\| tool \| `gpt-5\.4-mini` \|[^\n]*\| 0\.127 / 0\.093 \|")
grep("C18", "gpt-oss-20b not servable with a tool under the provider rule", DECISIONS, r"gpt-oss-20b is not servable with a tool through the endpoints the roster rule admits")

# ---------------------------------------------------------------- C19: judge independence
swap = table(RUN / "JUDGE_SWAP.md", r"^\| model \| traces \| milestones both judged \|")
rng("C19", "judge swap moves coverage by −0.017 to +0.020", [num(r[9]) for r in swap], -0.017, 0.02)
lojo = table(PILOT / "RESULTS_LOJO.md", r"^\| trace model \| judged / 60 \| F1, full panel \|")
family = [abs(num(r[6])) for r in lojo if r[4] == "family"]
check("C19", "dropping a family's own judge changes its score by at most 0.006", round(max(family), 3) == 0.006, str(family))
grep("C19", "detectable change 0.003 to 0.025", DECISIONS, r"detectable mean change at this size is 0\.003 to 0\.025")

# ---------------------------------------------------------------- C20, C21: the experts' readings
ER = RUN / "EXPERT_REQUEST.md"
grep("C20", "correct verdicts confirmed in 142 of 150", ER, r"the check said correct: experts said \| correct 142, incorrect 8")
grep("C20", "incorrect verdicts: 83 confirmed, 12 expert-correct, 3 partial", ER, r"the check said incorrect: experts said \| correct 12, incorrect 83, partial 3")
grep("C20", "partial verdicts called correct in 34 of 52", ER, r"the check said partial: experts said \| correct 34, incorrect 9, partial 9")
grep("C20", "judge REACHED confirmed 85 of 100 (0.850)", ER, r"REACHED confirmed \(precision\) \| 0\.850")
grep("C20", "judge MISSING confirmed 79 of 100 (0.790)", ER, r"MISSING confirmed \(not obtained, by either answer\) \| 0\.790")
grep("C20", "experts agree on 97% of double-read answers (kappa 0.947)", ER, r"150; 0\.973; 0\.947")
grep("C20", "experts agree on 88% of milestones (kappa 0.803)", ER, r"100; 0\.880; 0\.803")
grep("C21", "Fleiss' kappa 0.930 over three readers", ER, r"Fleiss' kappa over the three readers, all models: 0\.930")
grep("C21", "Claude Sonnet 5: 26 of 40 'No error'", ER, r"`claude-sonnet-5` \|[^|]*\| No error 26, 6\. Calculation Error 10")
grep("C21", "gpt-oss-20b: Calculation 17 of 40", ER, r"`gpt-oss-20b` \|[^|]*\| 6\. Calculation Error 17, 3\. Formula / Principle Error 9")
grep("C21", "Gemma 4: Calculation 28 of 40", ER, r"`gemma-4-26b-a4b` \|[^|]*\| 6\. Calculation Error 28, 3\. Formula / Principle Error 7")
grep("C21", "GPT-5.4 mini: Calculation 28 of 40", ER, r"`gpt-5\.4-mini` \|[^|]*\| 6\. Calculation Error 28, 3\. Formula / Principle Error 5")
easy = re.search(r"\| Easy \| ([^|]*) \|", read(ER)).group(1)
easy_counts = [int(x) for x in re.findall(r"(\d+)(?:,|$)", easy.replace(" |", ""))]
check("C21", "Easy wrong answers: Calculation 69 of 78 readings", "Calculation Error 69" in easy and sum(easy_counts) == 78, f"{easy} -> total {sum(easy_counts)}")

# ---------------------------------------------------------------- report
fails = [r for r in results if not r[2]]
width = max(len(r[1]) for r in results)
for claim, what, ok, observed in results:
    observed = " ".join(observed.split())[:110]
    print(f"{'PASS' if ok else 'FAIL'}  {claim:<5} {what:<{width}}  {observed}")
print(f"\n{len(results) - len(fails)} of {len(results)} checks pass; {len(fails)} fail.")
sys.exit(len(fails))
