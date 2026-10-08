"""Build the evaluation appendix's two tables from the validation reports, and check Section 4 and its appendix.

    python docs/appendix_evaluation.py           # print the table rows and the phrases the text uses
    python docs/appendix_evaluation.py --check   # exit 1 unless overleaf_source_04102026/5_evaluation.tex,
                                                 # appendices/scoring.tex and appendices/validation.tex hold
                                                 # every generated row and phrase

The tables: `tab:scoring_settings` (each check's setting against its alternatives on the expert study) in
scoring.tex, and `tab:validation_results` (agreement with the expert labels, planted defects, and the experts'
readings of the evaluated models) in validation.tex.

Sources. full_run_28092026/: SCORER_VALIDATION.md (the deterministic checks against the experts, current code),
THRESHOLD_APPENDIX.md (the tolerance's split-half fit, the milestone tolerance grid, the arithmetic rule's readings,
the judge's validation), ROUTER_VALIDATION.md (the judged step check), E5_VALIDATION.md (milestones settled by
matching on the 300 responses), JUDGE_SWAP.md (a second judge), FLAG_REVIEW_3.md (the expert reading of the
arithmetic flags), EXPERT_REQUEST.md (the experts' readings of verdicts on the evaluated models, its current figures:
the readings the store still holds, against the check's current verdicts), symbolic/SYMBOLIC_CHECK.md (the symbolic
equivalence step against the experts' grades), results/RESULTS.md (the share of milestones the judge decides),
results/sensitivity_variants.json (clause_variants.py: the rule's clauses re-applied one at a time, the accepted
answers by relative error), RESIDUAL_INCORRECT.md (the templates that prescribe their answer's digits) and analyze.py
(the four templates answerable from their wording); evaluator_pilot_17092026/evaluators/answer.py (the templates the
symbolic step is enabled on, the cap on an answer's own last digit). evaluator_pilot_17092026/: PILOT_SUMMARY.md (the
study's design, label agreement, the judge with matching, the reward models and their threshold, the planted defects
for three judges and the calls that returned no verdict), README.md (how the study's templates are chosen),
RESULTS_LOJO.md (a panel of judges without each family's own judge). docs/re-implementation-sep/
DECISIONS.md (D-181: the fourth judge on the planted defects).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "full_run_28092026"
PILOT = ROOT / "evaluator_pilot_17092026"
SRC = ROOT / "overleaf_source_04102026"
MAIN, SCORING, VALIDATION = SRC / "5_evaluation.tex", SRC / "appendices/scoring.tex", SRC / "appendices/validation.tex"


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def table(text: str, header: str) -> list[list[str]]:
    lines = text.split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("|") and re.search(header, l)) + 2
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip().strip("*`") for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return rows


def one(pattern: str, text: str) -> tuple[str, ...]:
    """A plain space in the pattern matches any run of whitespace, so wrapped lines in the record still match."""
    m = re.search(re.sub(r"(?<!\\) ", r"\\s+", pattern), text, re.S)
    if not m:
        raise SystemExit(f"not found: /{pattern}/")
    return m.groups()


def pct(x: float) -> str:
    return f"{round(100 * x)}\\%"


def signed(x: str) -> str:
    return f"${x}$"


def thousands(n: str) -> str:
    return f"{int(n):,}"


sv, ta = read(RUN / "SCORER_VALIDATION.md"), read(RUN / "THRESHOLD_APPENDIX.md")
ps, er = read(PILOT / "PILOT_SUMMARY.md"), read(RUN / "EXPERT_REQUEST.md")
psb = ps.replace("**", "").replace("*", "")

# ---------------------------------------------------------------- the expert study against the labels
now = {r[0]: r[4] for r in table(sv, r"^\| figure \| published \| published code")}
held = table(ta, r"^\| fitted on half \| fitted tolerance")
non_partial, panel = one(r"of the (\d+) traces they did not call partial, where the published check manages (0\.\d+)", ps)
e5 = one(r"\| deterministic, then a judge on the residue \(E5\) \| (0\.\d+) \| (0\.\d+) \| (0\.\d+) \|", psb)
router = {r[0]: r for r in table(read(RUN / "ROUTER_VALIDATION.md"), r"^\| steps the experts call incorrect")}
rc, rn = router["inside correct-answer traces, the router"], router["all traces, the router"]


def tp_fp(row: list[str]) -> tuple[int, int]:
    tp, fp, _fn = (int(x) for x in row[1].split(" / "))
    return tp, fp


# the judge's own flags: the router's minus the arithmetic check's, which the judge never sees (disjoint by design)
judge_all = [a - b for a, b in zip(tp_fp(rn), tp_fp(router["all traces, the digit rule alone"]))]
judge_cor = [a - b for a, b in zip(tp_fp(rc), tp_fp(router["inside correct-answer traces, the digit rule alone"]))]
matcher_then = now and table(sv, r"^\| figure \| published \| published code")
e3_f1_then = next(r[1] for r in matcher_then if r[0] == "E3 F1")       # the matcher the study's judge ran on
prm = one(r"ranks steps well \(AUROC (0\.\d+)\), but only about half its flags are real errors \(precision (0\.\d+)\), "
          r"and inside correct-answer traces that falls to (0\.\d+)", psb)
decides = one(r"AUROC (0\.\d+), against.*?Of (\d+) traces with a correct final answer, their holistic verdict calls "
              r"only (three) unsound", psb)
flawed = one(r"incorrect step in (\d+) of the (\d+) correct-answer traces \((\d+) steps\)", ps)
slips = one(r"the experts found three such steps against (\d+) calculation slips", ps)
outside = one(r"moved the pooled score from (0\.\d+) to (0\.\d+)", ps)
design = one(r"(\d+) problems drawn from five engineering branches and three difficulty levels: (\d+) templates", ps)
scalar = one(r"(\d+) of its (\d+) items \((\d+)%\) have a single scalar answer, against (\d+) of the benchmark's "
             r"(\d+) templates", ps)
experts = one(r"(\d+) domain experts, three per branch", ps)
rounds = one(r"All ([\d,]+) of them\..*re-labelled (\w+) of their own traces.*?(\d+) split steps across (\d+) traces", ps)
kappa = one(r"between experts, step labels \(Fleiss kappa\) \| \*\*(0\.\d+)\*\* \|.*within an expert, blind re-label "
            r"\(Cohen kappa\) \| \*\*(0\.\d+)\*\* \|.*between experts, milestone status \| (0\.\d+) \|.*"
            r"between experts, final-answer verdict \| (0\.\d+) \|", ps)

# ---------------------------------------------------------------- planted defects
planted = {r[0]: r[1:] for r in table(psb, r"^\| \| digit rule \(as E4 ships it\) \| GPT-5")}
grok = one(r"(\d+) of (\d+) conceptual defects caught \(0\.367; 0\.550 counting \"Other\"\), (\d+) of (\d+) arithmetic "
           r"\(0\.867\), (\d+) false alarms on the (\d+) untouched steps", read(ROOT / "docs/re-implementation-sep/DECISIONS.md"))
clean = one(r"From the (\d+) traces the experts called clean, (\d+) each receive exactly one defect and (\d+) are kept", ps)
routing = one(r"\| the flawed step is shown to a judge \| (0\.\d+) \|.*?\| the same step, unmodified, is shown \| "
              r"(0\.\d+) \|.*?caught by either of its two judges \| (0\.\d+) \|", psb)
batched = one(r"it sends (\d+) of the (\d+) conceptual defects to the judge.*?catches (\d+): 0\.317 end to end", ps)


def cnt(cell: str) -> str:
    return re.match(r"(\d+ of \d+)", cell).group(1)


judges = [(name, cnt(planted["conceptual defects"][k]), cnt(planted["arithmetic defects"][k]),
           cnt(planted["the same steps untouched, flagged"][k]))
          for k, name in enumerate(["GPT-5", "Claude Opus 4.5", "MiMo-V2.5-Pro"], 1)]   # the summary table's columns
judges.append(("Grok 4.6", f"{grok[0]} of {grok[1]}", f"{grok[2]} of {grok[3]}", f"{grok[4]} of {grok[5]}"))
rates = [int(c.split(" of ")[0]) / int(c.split(" of ")[1]) for _, c, _, _ in judges]
matching_caught = float(routing[2]) * 60
assert cnt(planted["conceptual defects"][0]) == "0 of 60", planted
assert abs(matching_caught - round(matching_caught)) < 1e-9, routing
assert max(rates) < 0.4, rates                                      # "at most about a third"
assert routing[0] == routing[1], routing                            # "no more often than the same step untouched"
assert int(batched[2]) > round(matching_caught)                     # the judged step check "catches more"
assert all(u.startswith("0 of") for _, _, _, u in judges)           # "without false alarms"

# ---------------------------------------------------------------- the scoring settings table
judge_val = {r[0]: r[1:] for r in table(ta, r"^\| shown to the judge \| REACHED")}
fv = judge_val["values x1.37 (should be MISSING), 88"]              # REACHED, MISSING, NOT_NEEDED, unjudged
G = {(r[0], r[1]): r for r in table(ta, r"^\| tolerance \| unit scaling \| real")}
A = {r[0]: r for r in table(ta, r"^\| reading \| hard case: tp / fp / fn")}
assert int(fv[2]) / 88 == 0.25, fv                                  # "excuses a quarter of such values"
settings_rows = [
    f"Final answer & $\\epsilon$ fitted on each half & {held[0][3]}, {held[1][3]} \\\\",
    f"Milestones & \\textbf{{0.5\\%, unit factors}} & \\textbf{{{G[('0.5%', 'yes')][4]}}} \\\\",
    f"& 0.2\\%, unit factors & {G[('0.2%', 'yes')][4]} \\\\",
    f"& 1\\%, unit factors & {G[('1%', 'yes')][4]} \\\\",
    f"& 2\\%, unit factors & {G[('2%', 'yes')][4]} \\\\",
    f"& 0.5\\%, no unit factors & {G[('0.5%', 'no')][4]} \\\\",
    f"Judge & \\textbf{{reached only}} & \\textbf{{{fv[0]} of 88}} \\\\",
    f"& reached or not needed & {int(fv[0]) + int(fv[2])} of 88 \\\\",
    f"Arithmetic & \\textbf{{digits shown, as used}} & \\textbf{{{A['digit rule, as shipped'][2]} / "
    f"{A['digit rule, as shipped'][3]}}} \\\\",
    f"& digits shown alone & {A['digit rule, bare'][2]} / {A['digit rule, bare'][3]} \\\\",
    f"& 0.1\\% tolerance & {A['0.1% tolerance'][2]} / {A['0.1% tolerance'][3]} \\\\",
    f"& 1\\% tolerance & {A['1% tolerance'][2]} / {A['1% tolerance'][3]} \\\\",
]
fitted = (held[0][1], held[1][1])

# ---------------------------------------------------------------- judge independence
js = read(RUN / "JUDGE_SWAP.md")
swap = table(js, r"^\| model \| traces \| milestones both judged")
diffs = [float(r[9]) for r in swap]
assert all(float(lo) <= 0 <= float(hi) for lo, hi in (r[10].split(" to ") for r in swap)), swap  # every CI holds zero
swap_total = one(r"on (\d+) sampled traces of the", js)[0]
swap_n = one(r"(\d+)(?: to (\d+))? traces per model over (\d+)(?: to (\d+))? templates", js)
swap_per_model = swap_n[0] if not swap_n[1] else f"{swap_n[0]} to {swap_n[1]}"
assert sum(int(r[1]) for r in swap) == int(swap_total), (swap_total, swap)
lojo = read(PILOT / "RESULTS_LOJO.md")
family = [abs(float(r[6])) for r in table(lojo, r"^\| trace model \| judged / 60 \| F1, full panel") if r[4] == "family"]
lenient: dict[str, tuple[int, int]] = {}
for r in table(lojo, r"^\| judge \| trace model \| steps \| bias \| 95% CI \| lenient"):  # the in-family panel
    share, n = re.match(r"([\d.]+) \((\d+)\)", r[5]).groups()
    a, b = lenient.get(r[0], (0, 0))
    lenient[r[0]] = (a + round(float(share) * int(n)), b + int(n))
lenient_pct = sorted(round(100 * a / b) for a, b in lenient.values())

# ---------------------------------------------------------------- the evaluated models
flags = {r[0]: r for r in table(read(RUN / "FLAG_REVIEW_3.md"), r"^\| model \| flags drawn \| read")}
fr = flags["all"]


def section(text: str, heading: str) -> str:
    """The text under a `### heading` of EXPERT_REQUEST.md, up to the next heading."""
    m = re.search(rf"^### {re.escape(heading)}\n(.*?)(?=^#)", text, re.S | re.M)
    if not m:
        raise SystemExit(f"EXPERT_REQUEST.md has no section '{heading}'")
    return m.group(1)


def total(s: str) -> int:
    return sum(int(n) for n in re.findall(r"\d+", s))


# The current figures: the readings the main store still holds, B1 against the check's current verdicts (two readers
# an item in B1 and B3, so a row's readings over two are its items)
er_b1, er_b3 = section(er, "B1 now"), section(er, "B3 now")
b1 = dict(re.findall(r"\| the check's current verdict (\w+): experts said \| ([^|]+) \|", er_b1))
b3 = dict(re.findall(r"\| the judge ruled (\w+): experts said \| ([^|]+) \|", er_b3))
b1_items, b3_items = (int(one(r"\| items; readings \| (\d+); (\d+) \|", s)[0]) for s in (er_b1, er_b3))
sampled = tuple(str(total(b1[v]) // 2) for v in ("correct", "incorrect", "partial"))
ruled = tuple(str(total(b3[v]) // 2) for v in ("MISSING", "REACHED"))
assert sum(map(int, sampled)) == b1_items and sum(map(int, ruled)) == b3_items, (sampled, ruled, b1_items, b3_items)
pairs = (one(r"\| items with two readers; agreement; Cohen's kappa \| \d+; (0\.\d+); (0\.\d+) \|", er_b1)
         + one(r"\| items with two readers; agreement; Cohen's kappa \| \d+; (0\.\d+); (0\.\d+) \|", er_b3))

# The symbolic equivalence step (D5) and the cap on an answer's own last digit (D11d)
answer_src = read(ROOT / "evaluator_pilot_17092026/evaluators/answer.py")
enabled = re.findall(r"'template_(\w+)'", one(r"SYMBOLIC_EQUIVALENCE_TEMPLATES = \((.*?)\)", answer_src)[0])
own_cap = float(one(r"(?m)^OWN_DIGIT_CAP = ([\d.]+)", answer_src)[0])
sc = read(RUN / "symbolic/SYMBOLIC_CHECK.md")
sym_graded = one(r"(?m)^\| all \| (\d+) \| \d+ \| \d+ \| \d+ \| \d+ \| (\d\.\d+) \| (\d\.\d+) \|", sc)
sym_earlier = one(r"(?m)^\| all \| (\d+) \| (\d\.\d+) \| (\d\.\d+) \| \d+ \| \d+ \| \d+ \| \d+ \|", sc)
symbolic_templates = len({json.loads(l)["template_id"] for l in read(RUN / "manifest.jsonl").splitlines()
                          if l.strip() and json.loads(l)["answer_type"] == "symbolic"})
WORDS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"]

# The scoring rule's clauses re-applied one at a time (clause_variants.py's results/sensitivity_variants.json), the
# accepted numeric answers by relative error, and the templates that prescribe their answer's digits
# (RESIDUAL_INCORRECT.md, residual_incorrect.py)
sv_main = json.loads(read(RUN / "results/sensitivity_variants.json"))["stores"]["main"]
where, sv_tau = sv_main["where_changed"], sv_main["tau_with_headline"]
bins = sv_main["relative_error_bins"]["all"]
accepted = bins["numeric"]
assert accepted == sum(bins[k] for k in ("le_0.2", "0.2_1", "1_5", "gt_5")), bins
exact = one(r"Templates with an exact-digits target \(the question prescribes the rounding\): (\d+) of 150", read(
    RUN / "RESIDUAL_INCORRECT.md"))[0]
assert all(m["abs_clause_off"] <= m["headline"] for m in sv_main["models"].values()), "the clause only lowers verdicts"
assert all(m["last_digit_unbounded"] >= m["headline"] for m in sv_main["models"].values()), "uncapped only raises"


def of_accepted(n: int) -> str:
    return f"{100 * n / accepted:.1f}\\%"


def q(s: str, label: str) -> str:
    return re.search(rf'"{label}" (\d+)', s).group(1)


def b1n(verdict: str, label: str) -> str:
    return re.search(rf"{label} (\d+)", b1[verdict]).group(1)


OBTAINS, UNNEEDED = "yes, the working obtains it", "no, and its route does not need it"
missing_confirmed = total(b3["MISSING"]) - int(q(b3["MISSING"], OBTAINS))
# Claims the prose makes from these readings: one the data no longer supports is reported and counted by --check,
# and the rows and phrases are still generated
CLAIMS_FAILED = [c for ok, c in (
    ((int(b1n("partial", "correct")) + 0.5 * int(b1n("partial", "partial"))) / total(b1["partial"]) > 0.5,
     f"scoring a partial verdict 0.5 is conservative on average (the experts' mean credit on "
     f"{total(b1['partial'])} readings)"),
    (0.30 <= int(q(b3["MISSING"], UNNEEDED)) / total(b3["MISSING"]) < 0.37,
     f"a third of the missing rulings are not needed by the route ({q(b3['MISSING'], UNNEEDED)} of {total(b3['MISSING'])})"),
) if not ok]

# ---------------------------------------------------------------- the validation results table
results_rows = [
    "\\multicolumn{3}{l}{\\textit{Agreement with the expert labels on 300 responses}} \\\\",
    f"Final-answer check & three-way, all 300 responses & {now['answer, three-way agreement']} \\\\",
    f"Final-answer check & {non_partial} responses not marked partial & {now['answer, non-partial agreement']} "
    f"(step matching {panel}) \\\\",
    f"Milestone matching & P / R / F1 & {now['E3 precision']} / {now['E3 recall']} / {now['E3 F1']} \\\\",
    f"Milestone matching with the judge$^\\dagger$ & P / R / F1 & {e5[0]} / {e5[1]} / {e5[2]} \\\\",
    f"Arithmetic check & P / R, all steps; correct answers & {now['digit rule, all traces, precision']} / "
    f"{now['digit rule, all traces, recall']}; {now['digit rule, hard case, precision']} / "
    f"{now['digit rule, hard case, recall']} \\\\",
    f"Judged step and arithmetic checks$^\\dagger$ & P / R, all steps; correct answers & {rn[2]} / {rn[3]}; "
    f"{rc[2]} / {rc[3]} \\\\",
    f"Judged step check alone$^\\dagger$ & P, all steps; correct answers & {judge_all[0] / sum(judge_all):.3f}; "
    f"{judge_cor[0] / sum(judge_cor):.3f} ({judge_cor[0]} of {sum(judge_cor)}) \\\\",
    f"Best process reward model & P, all steps; correct answers & {prm[1]}; {prm[2]} \\\\",
    "\\multicolumn{3}{l}{\\textit{Planted defects: conceptual / arithmetic caught, untouched steps flagged}} \\\\",
    f"Deterministic checks & conceptual & {cnt(planted['conceptual defects'][0])} \\\\",
] + [f"\\texttt{{{n}}}{'' if c.endswith(' of 60') and u.endswith(' of 120') else chr(36) + '^' + chr(92) + 'ddagger' + chr(36)}"
     f" & LLM judge, one step & {c} / {a} / {u} \\\\" for n, c, a, u in judges] + [
    f"Judged step check & conceptual, every unflagged step sent & {batched[2]} of {batched[1]} \\\\",
    f"Step matching with in-family judges & conceptual, end to end & {round(matching_caught)} of 60 \\\\",
    "\\multicolumn{3}{l}{\\textit{Readings of the evaluated models' responses}} \\\\",
    f"Arithmetic flags & {int(fr[3]) + int(fr[4])} decided flags & {fr[3]} real (precision {fr[6]}) \\\\",
    f"Final answer correct & {total(b1['correct'])} readings & {b1n('correct', 'correct')} confirmed \\\\",
    f"Final answer incorrect & {total(b1['incorrect'])} readings & {b1n('incorrect', 'incorrect')} confirmed, "
    f"{b1n('incorrect', 'correct')} called correct \\\\",
    f"Final answer partial & {total(b1['partial'])} readings & {b1n('partial', 'correct')} called fully correct, "
    f"{b1n('partial', 'incorrect')} incorrect \\\\",
    f"Judge reached & {total(b3['REACHED'])} readings & {q(b3['REACHED'], OBTAINS)} confirmed \\\\",
    f"Judge missing & {total(b3['MISSING'])} readings & {missing_confirmed} confirmed, {q(b3['MISSING'], UNNEEDED)} "
    f"of them not needed by the route \\\\",
]

# ---------------------------------------------------------------- the phrases each file must contain
res = read(RUN / "results/RESULTS.md")
judged = [float(r[3]) for r in table(res, r"^\| model \| calls \| without a reply \| judged fraction \|")]
gold = one(r"\| digit-rule flags \| 0 of (\d+) claims in (\d+) steps \|", sv)
gold_ok = one(r"\| scored correct at all three tolerances \| (\d+) of (\d+) \|", sv)
no_ms = one(r"\| every milestone found \| (\d+); (\d+) items have no milestones \|", sv)
settled = one(r"\| milestones: by E3, judged \| (\d+) of (\d+), (\d+) \(published", read(RUN / "E5_VALIDATION.md"))
readers = one(r"Every trace was labelled by (three) experts of its own branch", psb)
per_model = max(int(r[1]) for k, r in flags.items() if k != "all")
assert gold_ok[0] == gold_ok[1] and decides[2] == "three" and rounds[1] == "seven"

x1 = read(PILOT / "RESULTS_X1.md")
split = one(r"(\d+) split steps across (\d+) traces \((\d+) distinct steps", x1.replace("**", ""))
other = max(float(v) for v in re.findall(r"(0\.\d+) counting \*?Other\*?", x1))
assert abs(other * 60 - round(other * 60)) < 0.02, other
detectable = one(r"detectable mean change at this size is (0\.\d+) to (0\.\d+)",
                 read(ROOT / "docs/re-implementation-sep/DECISIONS.md"))
loj = table(lojo, r"^\| trace model \| judged / 60 \| F1, full panel")
for r in loj:   # every own-family drop changes the score exactly as much as dropping one of the other two judges
    if r[4] == "family":
        assert any(o[0] == r[0] and o[4] == "placebo" and o[6] == r[6] for o in loj), r
frontier: dict[str, tuple[int, int]] = {}
for r in table(lojo, r"^\| judge \| trace model \| steps \| bias \| 95% CI \| lenient"):
    if r[1] != "llama-3.1-70b":
        share, n = re.match(r"([\d.]+) \((\d+)\)", r[5]).groups()
        a, b = frontier.get(r[0], (0, 0))
        frontier[r[0]] = (a + round(float(share) * int(n)), b + int(n))
frontier_pct = sorted(round(100 * a / b) for a, b in frontier.values())
assert 0.2 <= float(prm[2]) < 0.3, prm                               # "three of four of its flags ... are false"
assert fv[0] == "0", fv                                               # "calls none of the 88"
per_template = int(design[0]) // int(design[1])
assert int(scalar[0]) % per_template == 0, scalar
not_confirmed = total(b3["REACHED"]) - int(q(b3["REACHED"], OBTAINS))

main = [
    f"on 300 responses from five LLMs outside the eleven we evaluate, each labeled step by step by {readers[0]} "
    f"domain experts",
    f"three-way verdict on {now['answer, three-way agreement']} of responses",
    f"milestone matching with the judge reaches F1 {e5[2]}",
    f"{pct(min(judged))} to {pct(max(judged))} of milestones",
    f"on {planted['conceptual defects'][0].split(' of ')[1]} planted misstated rules behind correct answers",
]
scoring = [
    f"Every one of the {thousands(gold_ok[1])} gold traces scores correct",
    f"the {no_ms[1]} instances without milestones",
    f"it decides {settled[2]} of the {thousands(settled[1])} milestones of the expert study",
    "excuses a quarter", f"calls none of the 88 unstated values",
    f"fitted values {fitted[0]} and {fitted[1]}",
    f"among {88} milestone values multiplied by 1.37",
    f"none of the {thousands(gold[0])} calculations it reads in the {thousands(gold_ok[1])} gold traces",
    # the symbolic equivalence step and the cap on an answer's own last digit (ANSWER FINAL)
    f"{WORDS[len(enabled)]} of the {WORDS[symbolic_templates]} templates with symbolic answers",
    f"precision {sym_graded[1]} and recall {sym_graded[2]} on the {sym_graded[0]} verdicts the experts graded",
    f"precision {sym_earlier[1]} and recall {sym_earlier[2]} on the {sym_earlier[0]} readings of these templates",
    f"The other {WORDS[symbolic_templates - len(enabled)]} templates with symbolic answers",
    f"at most {own_cap * 100:g}\\% of the target",
    # the prescribed digits and the clause sensitivities (tab:scoring_variants)
    f"{exact} templates prescribe the rounding of their answer",
    f"relaxing it to the tolerance raises {where['prescribed_relaxed']['verdicts']} verdicts on "
    f"{where['prescribed_relaxed']['templates']} templates and keeps Kendall's $\\tau$ with the ordering as scored at "
    f"{sv_tau['prescribed_relaxed']:.3f}",
    f"lowers {where['abs_clause_off']['verdicts']} verdicts on {where['abs_clause_off']['templates']} templates "
    f"($\\tau = {sv_tau['abs_clause_off']:.3f}$)",
    f"raises {where['last_digit_unbounded']['verdicts']} verdicts on {where['last_digit_unbounded']['templates']} "
    f"templates ($\\tau = {sv_tau['last_digit_unbounded']:.3f}$)",
    f"moves {where['per_part_credit']['verdicts']} verdicts on {where['per_part_credit']['templates']} templates "
    f"($\\tau = {sv_tau['per_part_credit']:.3f}$)",
    f"Of the {thousands(str(accepted))} accepted numeric answers, {of_accepted(bins['le_0.2'])} lie within 0.2\\% of their "
    f"target, {of_accepted(bins['0.2_1'])} between 0.2\\% and 1\\%",
    f"and {of_accepted(bins['1_5'] + bins['gt_5'])} further from it",
] + settings_rows
validation = [
    f"four instances from each of {design[1]} templates, {design[0]} in all",
    f"only {int(scalar[0]) // per_template} of its templates have a single scalar answer, against {scalar[3]} of the "
    f"{scalar[4]}",
    "Five LLMs outside the evaluated models", "giving 300 responses",
    "Fifteen domain experts, three per branch" if experts[0] == "15" else "EXPERTS CHANGED",
    f"a written reason for each of the {rounds[0]} incorrect step labels", f"of {rounds[1]} responses per expert",
    f"the {split[2]} distinct steps on which experts split",
    f"Fleiss'~$\\kappa$~\\citep{{fleiss1971}} {kappa[0]} on steps, {kappa[2]} on milestones, and {kappa[3]} on "
    f"final answers",
    f"two passes at Cohen's~$\\kappa$~\\citep{{cohen1960}} {kappa[1]}",
    f"(area under the receiver operating characteristic curve, AUROC, {decides[0]})",
    f"Of the {flawed[2]} incorrect steps behind correct answers, {slips[0]} are calculation slips",
    f"(AUROC {prm[0]})", "three of four of its flags inside correct answers are false",
    f"({outside[0]} with two in-family judges, {outside[1]} with three outside ones)",
    f"each of {clean[1]} clean responses",
    f"up to {round(other * 60)} of 60",
    f"Re-judging {swap_total} sampled responses, {swap_per_model} per model",
    f"{signed(f'{min(diffs):+.3f}')} to {signed(f'{max(diffs):+.3f}')}",
    f"at most {max(family):.3f}", f"(the design detects changes of {detectable[0]} to {detectable[1]})",
    f"{frontier_pct[0]}\\% to {frontier_pct[-1]}\\%",
    f"up to {per_model} flagged calculations per model",
    f"({sampled[0]} correct, {sampled[1]} incorrect, and {sampled[2]} partial)",
    f"({ruled[1]} reached and {ruled[0]} missing)",
    f"{pct(float(pairs[0]))} of answers and {pct(float(pairs[2]))} of milestones",
    f"{not_confirmed} of the {total(b3['REACHED'])} readings of reached rulings",
    f"F1 {e3_f1_then} without the judge",
    # how the study's templates are chosen (the pilot's README), the in-sample caveat, the two milestone F1 values
    f"{design[1]} templates, 60 in all, one template per branch and level" if one(
        r"(\d+) templates, one per branch × level cell", read(PILOT / "README.md"))[0] == design[1] else "SLICE CHANGED",
    "leaves out the four templates whose answers can be read off" if sorted(
        re.findall(r"`(\w+)`", one(r"(Four templates are excluded everywhere.*?)\n\n", read(PILOT / "README.md"))[0]))
    == sorted(t.replace("template_", "") for t in re.findall(r"'(template_\w+)'", one(
        r"SHORTCUT = \[(.*?)\]", read(RUN / "analyze.py"))[0])) else "READ-OFF TEMPLATES CHANGED",
    "Only $\\epsilon$ is cross-fitted",
    f"the two milestone F1 values, {now['E3 F1']} and {e5[2]}",
    # the judge's calls that returned no verdict, and the reward models' threshold
    "{4} of its {5} calls returned no reply after retries, so its counts rest on the {0} conceptual defects and "
    "the {2} untouched steps".format(*one(r"MiMo returned a verdict on both arms for (\d+) of the (\d+) conceptual "
                                          r"defects and for (\d+) of the (\d+) untouched steps: (\d+) of its (\d+) calls "
                                          r"never returned", psb)),
    "changes step F1 by ${}$".format(one(r"changes the held-out F1 by (−0\.\d+)", psb)[0].replace("−", "-")),
    "no threshold reaches precision {}".format(one(r"no threshold reaches precision (0\.\d+) on held-out data", psb)[0]),
    # the readings' denominators and the judged step check's precision beside tab:judged_steps
    f"each of {b1_items} final-answer verdicts", f"and {b3_items} judge rulings",
    f"The experts call {b1n('partial', 'correct')} of their {total(b1['partial'])} readings of partial verdicts "
    f"fully correct and {b1n('partial', 'incorrect')} incorrect",
    f"{b1n('incorrect', 'correct')} of the {total(b1['incorrect'])} readings of incorrect verdicts call the answer "
    "correct",
    f"precision is {judge_all[0] / sum(judge_all):.3f} over all steps and {judge_cor[0] / sum(judge_cor):.3f} inside "
    f"correct answers ({judge_cor[0]} of {sum(judge_cor)} flags)",
] + results_rows

for c in CLAIMS_FAILED:
    print(f"CLAIM FAILS: {c}")
if "--check" in sys.argv:
    bad = len(CLAIMS_FAILED)
    for path, wanted in ((MAIN, main), (SCORING, scoring), (VALIDATION, validation)):
        tex = " ".join(read(path).split()) if path.exists() else ""
        missing = [s for s in wanted if " ".join(s.split()) not in tex]
        for s in missing:
            print(f"MISSING in {path.name}: {s}")
        print(f"{len(wanted) - len(missing)} of {len(wanted)} generated rows and phrases are in {path.name}")
        bad += len(missing)
    sys.exit(1 if bad else 0)

for title, rows in (("scoring settings", settings_rows), ("validation results", results_rows)):
    print(f"% {title}")
    print("\n".join(rows))
for title, ph in (("5_evaluation.tex", main), ("scoring.tex", scoring), ("validation.tex", validation)):
    print(f"% phrases {title} must contain")
    print("\n".join(p for p in ph if " & " not in p))
