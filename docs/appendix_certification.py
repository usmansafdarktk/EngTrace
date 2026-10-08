"""Build the rows of the appendix's certification tables from the certification reports, and check the .tex.

    python docs/appendix_certification.py           # print the table rows and the counts the text uses
    python docs/appendix_certification.py --check   # exit 1 unless overleaf_source_04102026/appendices/certification.tex
                                                    # holds every generated row and number

Sources (template_annotation_23092026/): layer2/RESULTS.md and RESULTS_round2-6.md (plants, hand checks, verdicts,
agreement), layer2/CERTIFICATION.md (rounds and the final status), layer0/gate_report.md (integrity checks),
README.md (the screen's two passes), screen/pass2/stats.md, and docs/re-implementation-sep/DECISIONS.md (the
screen's verified claims). The planted defects' short descriptions are written here; their templates and kinds are
read from the report and must agree.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "template_annotation_23092026"
L2 = CERT / "layer2"
TEX = ROOT / "overleaf_source_04102026/appendices/certification.tex"
BRANCHES = ["chemical", "electrical", "civil", "industrial", "mechanical"]  # the taxonomy figure's order
CODE = {"chemical": "che", "electrical": "ele", "civil": "civ", "industrial": "ind", "mechanical": "mec"}
KINDS = ["constant", "unit", "sign", "arithmetic"]


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def table(text: str, header: str) -> list[list[str]]:
    lines = text.split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("|") and re.search(header, l)) + 2
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return rows


def one(pattern: str, text: str) -> tuple[str, ...]:
    m = re.search(pattern, text, re.S)
    if not m:
        raise SystemExit(f"not found: /{pattern}/")
    return m.groups()


# ---------------------------------------------------------------- planted defects
DESCRIBE = {  # template -> (topic, what was planted); kinds and templates are checked against the report
    "vdw_solve_for_pressure": ("van der Waals pressure", "$R$ with two digits transposed"),
    "annulus_flowrate": ("annulus flow rate", "kPa to Pa by a factor of 100"),
    "heat_of_reaction_formation": ("heat of reaction", "Hess's law reversed"),
    "batch_reactor_first_order": ("first-order batch reactor", "one printed time 4\\% too high"),
    "coulombs_law": ("Coulomb's law", "$\\varepsilon_0$ ten times too large"),
    "wave_parameters_basic": ("wave parameters", "MHz to Hz by a factor of $10^5$"),
    "cd_dc_system_analysis": ("C/D and D/C conversion", "a delay written as an advance"),
    "mean_variance": ("mean and variance", "one printed mean 0.2 too high"),
    "manning_rectangular_discharge": ("Manning discharge", "$R^{3/4}$ for $R^{2/3}$"),
    "cantilever_double_integration": ("cantilever deflection", "m to mm by a factor of 100"),
    "effective_stress_profile": ("effective stress", "$\\sigma' = \\sigma + u$"),
    "beam_support_reactions": ("beam reactions", "one printed reaction 4\\% too high"),
    "safety_stock_reorder_point": ("safety stock", "a two-sided $z$ for a one-sided one"),
    "mm1_time_in_system": ("M/M/1 time in system", "hours to minutes by a factor of 100"),
    "xbar_r_control_limits": ("$\\bar{x}$ and $R$ chart limits", "the lower limit as center plus $A_2\\bar{R}$"),
    "cp_cpk_from_specs": ("process capability", "one printed $C_p$ 6\\% too high"),
    "shear_stress_torsion": ("torsional shear stress", "$J$ with $\\pi/4$ for $\\pi/2$"),
    "hydrostatic_pressure_at_depth": ("hydrostatic pressure", "kPa to Pa by a factor of 100"),
    "utube_manometer": ("U-tube manometer", "the pipe-fluid column's sign flipped"),
    "undamped_natural_frequency_translational": ("natural frequency", "one printed $\\omega_n$ 4\\% too high"),
}
res1 = read(L2 / "RESULTS.md")
plants = {}
for code, tid, kind in re.findall(r"^\| plant_(\w{3})_\d \(template_(\w+)\) \| (\w+):", res1, re.M):
    plants[(code, kind)] = tid
assert len(plants) == 20 and all((CODE[b], k) in plants for b in BRANCHES for k in KINDS), plants
plant_rows = []
for b in BRANCHES:
    cells = []
    for k in KINDS:
        topic, what = DESCRIBE[plants[(CODE[b], k)]]
        cells.append(f"{topic}: {what}")
    plant_rows.append(f"{b.capitalize()} & " + " & ".join(cells) + " \\\\")
plant_readings = one(r"Overall: (\d+) of (\d+) planted defects rejected", res1)

# ---------------------------------------------------------------- agreement among a branch's three experts (round 1)
agree = {r[0]: r for r in table(res1, r"\| Branch \| Templates with 3 verdicts")}
# a branch whose experts approve every real template has one category only: Fleiss' kappa is 0/0, undefined, though
# the report prints 1.000; the cell says so
rejecting_r1 = {b: 0 for b in BRANCHES}
for r in table(res1, r"\| Expert \| Branch \| Plants seen"):
    rejecting_r1[r[1]] += 30 - int(r[5].split()[0])
undefined_kappa = [b for b in BRANCHES if rejecting_r1[b] == 0]
agreement_rows = []
for b in BRANCHES + ["all"]:
    r = agree[b]
    name = "\\textbf{All}" if b == "all" else b.capitalize()
    pct = r[4].replace("%", "\\%")
    kappa = "--" if b in undefined_kappa else r[2]
    agreement_rows.append(f"{name} & {kappa} & {r[3]} & {pct} & {r[5]} & {r[6]} & {r[7]} \\\\")

# ---------------------------------------------------------------- the rounds (stated in the Results paragraph)
cert = read(L2 / "CERTIFICATION.md")
rounds = {int(r[0]): (int(r[3]), int(r[4])) for r in table(cert, r"\| Round \| Labels")}
reports = {n: res1 if n == 1 else read(L2 / f"RESULTS_round{n}.md") for n in rounds}
R = {}
for n in sorted(rounds):
    text = reports[n]
    if n == 1:  # the first table counts real templates approved out of 30 per expert
        rejecting = sum(30 - int(r[5].split()[0]) for r in table(text, r"\| Expert \| Branch \| Plants seen"))
    else:
        rejecting = sum(int(r[4]) for r in table(text, r"\| Expert \| Branch \| Items \| Approved \| Rejected"))
    rejected_templates = len(re.findall(r"^- \*\*template_\w+\*\* rejected by", text, re.M))
    compared, matched = one(r"(\d+) hand checks with a comparable number: (\d+) matched", text)
    templates, verdicts = rounds[n]
    R[n] = dict(templates=templates, verdicts=verdicts, rejecting=rejecting, rejected=rejected_templates,
                matched=matched, compared=compared)
assert sorted(R) == [1, 2, 3, 4, 5, 6], R
# the text: the third round reviews the second round's rejected templates and approves all; the fourth rejects none;
# the fifth reviews the two repaired chemical templates, one expert rejecting one of them; the sixth reviews that one
# template and approves it
assert R[3]["templates"] == R[2]["rejected"] and R[3]["rejecting"] == R[4]["rejecting"] == 0, R
assert R[6]["templates"] == R[5]["rejected"] and R[6]["rejecting"] == 0, R
WORD = ["no", "one", "two", "three", "four", "five", "six", "seven"]
later = [n for n in sorted(R) if n > 1]
# what the text says of the later rounds, checked against their records: round 3's two mismatched hand checks are
# the U-tube manometer's gold value entered in pascals (its answer line is in kilopascals); round 4 widens two
# templates' sampled parameters to 15 distinct problems each; round 5 states the reading and the data in the two
# chemical questions; round 6 holds organic compounds below their decomposition range
r3_entered = re.findall(r"^- template_(\w+) by \S+: entered `([\d.]+) (\w+)`", reports[3], re.M)
assert len(r3_entered) == int(R[3]["compared"]) - int(R[3]["matched"]) and {t for t, _, _ in r3_entered} == {
    "utube_manometer"} and {u for _, _, u in r3_entered} == {"Pa"}, r3_entered
assert re.search(r"widened so each can supply 15 distinct problems", read(L2 / "fixes_round4.md")), "round 4"
fixes5 = read(L2 / "fixes_round5.md")
assert all(s in fixes5 for s in ("The question now pins that reading", "the data the energy balance consumes",
                                 "Organic compounds are now capped at")), "rounds 5 and 6"
rows_r1 = one(r"from (\d+) label rows by \d+ experts", res1)[0]
round_phrases = [
    f"{WORD[len(R)]} rounds",
    # the first round's hand checks: one per label row, of which some give no number to compare
    f"Of the {rows_r1} first-round hand checks, {R[1]['compared']} have a number to compare; {R[1]['matched']} of "
    f"these match the gold value within 1\\%",
    f"{R[1]['rejecting']} of the {R[1]['verdicts']} decisions on templates are rejections",
    f"{R[2]['rejecting']} of the {R[2]['verdicts']} decisions are rejections, on {R[2]['rejected']} of the "
    f"{R[2]['templates']} templates",
    f"the third round approves all {R[3]['templates']}",
    f"a fourth round re-certifies {WORD[R[4]['templates']]} templates",
    f"a fifth round re-certifies {WORD[R[5]['templates']]} templates",
    f"{WORD[R[5]['rejecting']]} of its {R[5]['verdicts']} decisions",
    f"a sixth round re-certifies {WORD[R[6]['templates']]} template{'s' if R[6]['templates'] != 1 else ''}",
    ", ".join(f"{R[n]['matched']} of {R[n]['compared']}" for n in later[:-1])
    + f", and {R[later[-1]]['matched']} of {R[later[-1]]['compared']} hand checks",
    "in pascals where the gold answer is in kilopascals",
    # the agreement table's undefined cells, and electrical's kappa beside its percent agreement
    "undefined for " + " and ".join(undefined_kappa),
    f"in electrical, {rejecting_r1['electrical']} of the {3 * int(agree['electrical'][1])} decisions on templates are "
    f"rejections, so {agree['electrical'][4].replace('%', chr(92) + '%')} agreement gives a $\\kappa$ of "
    f"{agree['electrical'][2]}",
]

# ---------------------------------------------------------------- counts the text states
readme = read(CERT / "README.md")
gate = read(CERT / "layer0/gate_report.md")
stats = read(CERT / "screen/pass2/stats.md")
decisions = read(ROOT / "docs/re-implementation-sep/DECISIONS.md")
p1 = one(r"pass 1 \(2026-09-23\) (\d+) pass, (\d+) controversial, (\d+) critical failure, AC1 on the flag (0\.\d+); "
         r"the (\d+) flags verified and fixed", readme)
p2 = one(r"\*\*(\d+) pass, (\d+) controversial, (\d+) critical\*\*, AC1 (0\.\d+)", readme)
claims = one(r"(\d+) of (\d+) were real", decisions)
fp = one(r"passed and the experts reject by majority: (\d+) of (\d+) \(false-positive rate ([\d.]+)%\)", res1)
mad = dict(re.findall(r"\| (physical_plausibility|mathematical_correctness|pedagogical_clarity) \| ([\d.]+) \|", res1))
register = table(gate, r"\| Template \| Pattern \| Lines absorbed")
mismatch = one(r"Of the (\d+) mismatches, (\d+) ended in a rejection", res1)
majority = one(r"lists the (\d+) templates at least one expert rejected \((\d+) by\s+majority\)", read(L2 / "fixes_round1.md"))
by_branch = re.findall(r"^- \*\*template_(\w+)\*\* rejected by", res1, re.M)
branch_of = dict(re.findall(r"^\| `template_(\w+)` \| (\w+) \|", cert, re.M))
rej_by_branch = {b: sum(1 for t in by_branch if branch_of.get(t) == b) for b in BRANCHES}

seeds = one(r"(\d+) seeds per template", gate)

prose = [
    f"every template at {seeds[0]} seeds", f"six line patterns in four templates" if (len(register), len({r[0] for r in register})) == (6, 4) else "REGISTER CHANGED",
    f"{p1[0]} templates pass, {p1[1]} are controversial, and {p1[2]} are critical failures", f"AC1 of {p1[3]}",
    f"{claims[0]} of the {claims[1]} claims", f"the {p1[4]} flagged templates",
    f"{p2[0]} templates pass, {p2[1]} are controversial, and "
    f"{'none is a critical failure' if p2[2] == '0' else p2[2] + ' are critical failures'}", f"AC1 of {p2[3]}",
    f"{fp[0]} of the {fp[1]} templates", f"{fp[2]}\\%",
    f"{mad['physical_plausibility']} for physical plausibility, {mad['mathematical_correctness']} for mathematical "
    f"correctness, and {mad['pedagogical_clarity']} for pedagogical clarity",
    f"{plant_readings[0]} of {plant_readings[1]}", f"{mismatch[1]} of the {mismatch[0]} mismatches",
    f"{majority[0]} templates", f"{majority[1]} of them by a majority",
    f"chemical {rej_by_branch['chemical']}, electrical {rej_by_branch['electrical']}, and mechanical "
    f"{rej_by_branch['mechanical']}" if rej_by_branch["civil"] == rej_by_branch["industrial"] == 0 else "BRANCHES CHANGED",
]

if "--check" in sys.argv:
    tex = " ".join(TEX.read_text(encoding="utf-8").split())
    wanted = plant_rows + agreement_rows + round_phrases + prose
    missing = [s for s in wanted if " ".join(s.split()) not in tex]
    for s in missing:
        print("MISSING:", s)
    print(f"{len(wanted) - len(missing)} of {len(wanted)} generated rows and numbers are in {TEX.name}")
    sys.exit(1 if missing else 0)

for title, rows in (("planted defects", plant_rows), ("agreement, round 1", agreement_rows), ("rounds, in the prose", round_phrases)):
    print(f"% {title}")
    print("\n".join(rows))
print("% phrases the prose must contain")
print("\n".join(prose))
