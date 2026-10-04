"""Build the rows of the evaluation tables from the validation reports, and check Section 4 and its appendix.

    python docs/appendix_evaluation.py           # print the table rows and the phrases the text uses
    python docs/appendix_evaluation.py --check   # exit 1 unless overleaf_source_04102026/5_evaluation.tex,
                                                 # appendices/scoring.tex and appendices/validation.tex hold
                                                 # every generated row and phrase

Sources. full_run_28092026/: SCORER_VALIDATION.md (the deterministic checks against the experts, current code),
THRESHOLD_APPENDIX.md (the tolerance's split-half fit, the milestone tolerance grid, the arithmetic rule's readings,
the judge's validation), ROUTER_VALIDATION.md (the judged step check), E5_VALIDATION.md (milestones settled by
matching on the 300 responses), JUDGE_SWAP.md (a second judge), FLAG_REVIEW_3.md (the expert reading of the
arithmetic flags), EXPERT_REQUEST.md (the experts' readings of verdicts on the evaluated models) and
results/RESULTS.md (the share of milestones the judge decides). evaluator_pilot_17092026/: PILOT_SUMMARY.md (the
study's design, label agreement, the judge with matching, the reward models, the planted defects for three
judges), RESULTS_LOJO.md (a panel of judges without each family's own judge). docs/re-implementation-sep/
DECISIONS.md (D-181: the fourth judge on the planted defects). Model display names are written here.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "full_run_28092026"
PILOT = ROOT / "evaluator_pilot_17092026"
SRC = ROOT / "overleaf_source_04102026"
MAIN, SCORING, VALIDATION = SRC / "5_evaluation.tex", SRC / "appendices/scoring.tex", SRC / "appendices/validation.tex"

NAME = {  # the evaluated models, in the order of their Final Answer Accuracy
    "deepseek-v4.1-flash": "DeepSeek V4.1 Flash", "kimi-k3": "Kimi K3", "claude-sonnet-5": "Claude Sonnet 5",
    "glm-5.3-flash": "GLM-5.3-Flash", "muse-glimmer-30b": "Muse Glimmer 30B", "glm-5.3": "GLM-5.3",
    "qwen3-235b-a22b-2507": "Qwen3-235B-2507", "gemini-3.1-flash-lite": "Gemini 3.1 Flash-Lite",
    "gemma-4-26b-a4b": "Gemma 4 26B", "gpt-5.4-mini": "GPT-5.4 mini", "gpt-oss-20b": "gpt-oss-20b",
}


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


sv, ta = read(RUN / "SCORER_VALIDATION.md"), read(RUN / "THRESHOLD_APPENDIX.md")
ps, er = read(PILOT / "PILOT_SUMMARY.md"), read(RUN / "EXPERT_REQUEST.md")

# ---------------------------------------------------------------- the deterministic checks against the experts
now = {r[0]: r[4] for r in table(sv, r"^\| figure \| published \| published code")}
held = table(ta, r"^\| fitted on half \| fitted tolerance")
non_partial, panel = one(r"of the (\d+) traces they did not call partial, where the published check manages (0\.\d+)", ps)
e5 = one(r"\| deterministic, then a judge on the residue \(E5\) \| (0\.\d+) \| (0\.\d+) \| (0\.\d+) \|",
         ps.replace("**", ""))
router = {r[0]: r for r in table(read(RUN / "ROUTER_VALIDATION.md"), r"^\| steps the experts call incorrect")}
prm = one(r"ranks steps well \(AUROC (0\.\d+)\), but only about half its flags are real errors \(precision (0\.\d+)\), "
          r"and inside correct-answer traces that falls to (0\.\d+)", ps.replace("**", ""))
rc, rn = router["inside correct-answer traces, the router"], router["all traces, the router"]
comparison_rows = [
    f"Final-answer check & all 300 responses & three-way agreement & {now['answer, three-way agreement']} \\\\",
    f"Final-answer check & {non_partial} not marked partial & agreement & {now['answer, non-partial agreement']} \\\\",
    f"Step-matching designs' answer check & {non_partial} not marked partial & agreement & {panel} \\\\",
    f"Final-answer check, $\\epsilon$ fitted on the other half & each half & three-way agreement & "
    f"{held[0][3]}, {held[1][3]} \\\\",
    f"Milestone matching & milestones & P / R / F1 & {now['E3 precision']} / {now['E3 recall']} / {now['E3 F1']} \\\\",
    f"Matching, then the judge & milestones & P / R / F1 & {e5[0]} / {e5[1]} / {e5[2]} \\\\",
    f"Arithmetic check & steps, all responses & P / R & {now['digit rule, all traces, precision']} / "
    f"{now['digit rule, all traces, recall']} \\\\",
    f"Arithmetic check & steps, correct answers & P / R & {now['digit rule, hard case, precision']} / "
    f"{now['digit rule, hard case, recall']} \\\\",
    f"Judged step and arithmetic checks & steps, all responses & P / R & {rn[2]} / {rn[3]} \\\\",
    f"Judged step and arithmetic checks & steps, correct answers & P / R & {rc[2]} / {rc[3]} \\\\",
    f"Best process reward model & steps, all responses & P & {prm[1]} \\\\",
    f"Best process reward model & steps, correct answers & P & {prm[2]} \\\\",
]

# ---------------------------------------------------------------- the planted defects, per judge
planted = {r[0]: r[1:] for r in table(ps.replace("**", ""), r"^\| \| digit rule \(as E4 ships it\) \| GPT-5")}
grok = one(r"(\d+) of (\d+) conceptual defects caught \(0\.367; 0\.550 counting \"Other\"\), (\d+) of (\d+) arithmetic "
           r"\(0\.867\), (\d+) false alarms on the (\d+) untouched steps", read(ROOT / "docs/re-implementation-sep/DECISIONS.md"))


def cnt(cell: str) -> str:
    return re.match(r"(\d+ of \d+)", cell).group(1)


planted_rows = []
for k, name in enumerate(["GPT-5", "Claude Opus 4.5", "MiMo-V2.5-Pro"], 1):  # the summary table's column order
    planted_rows.append(f"\\texttt{{{name}}} & {cnt(planted['conceptual defects'][k])} & "
                        f"{cnt(planted['arithmetic defects'][k])} & {cnt(planted['the same steps untouched, flagged'][k])} \\\\")
planted_rows.append(f"\\texttt{{Grok 4.6}} & {grok[0]} of {grok[1]} & {grok[2]} of {grok[3]} & {grok[4]} of {grok[5]} \\\\")
assert cnt(planted["conceptual defects"][0]) == "0 of 60", planted
rates = [int(a) / int(b) for a, b in (re.match(r"(\d+) of (\d+)", r.split(" & ")[1]).groups() for r in planted_rows)]

# ---------------------------------------------------------------- the scoring tables
judge_val = {r[0]: r[1:] for r in table(ta, r"^\| shown to the judge \| REACHED")}
tv, fv = judge_val["true values (should be REACHED), 88"], judge_val["values x1.37 (should be MISSING), 88"]
judge_rows = [f"True values (88) & {tv[0]} & {tv[1]} & {tv[2]} \\\\",
              f"Values $\\times 1.37$ (88) & {fv[0]} & {fv[1]} & {fv[2]} \\\\"]
grid_rows = [r[0].replace("%", "\\%") + f" & {r[1]} & {r[2]} & {r[3]} & {r[4]} \\\\"
             for r in table(ta, r"^\| tolerance \| unit scaling \| real")]
READING = {"1% tolerance": "1\\% tolerance", "0.1% tolerance": "0.1\\% tolerance",
           "digit rule, bare": "Digits shown", "digit rule, as shipped": "Digits shown, as used"}
arith_rows = [f"{READING[r[0]]} & {r[2]} & {r[3]} & {r[4]} & {r[5]} & {r[6]} & {r[7]} \\\\"
              for r in table(ta, r"^\| reading \| hard case: tp / fp / fn")]

# ---------------------------------------------------------------- judge independence
swap = {r[0]: r for r in table(read(RUN / "JUDGE_SWAP.md"), r"^\| model \| traces \| milestones both judged")}
swap_rows = []
for key, name in NAME.items():
    r = swap[key]
    mimo, other = r[8].split(" / ")
    lo, hi = r[10].split(" to ")
    swap_rows.append(f"\\texttt{{{name}}} & {r[2]} & {r[4]} & {mimo} & {other} & {signed(r[9])} & "
                     f"{signed(lo)} to {signed(hi)} \\\\")
diffs = [float(swap[k][9]) for k in NAME]
lojo = read(PILOT / "RESULTS_LOJO.md")
family = [abs(float(r[6])) for r in table(lojo, r"^\| trace model \| judged / 60 \| F1, full panel") if r[4] == "family"]
lenient: dict[str, tuple[int, int]] = {}
for r in table(lojo, r"^\| judge \| trace model \| steps \| bias \| 95% CI \| lenient"):  # the first panel: E0-3J
    share, n = re.match(r"([\d.]+) \((\d+)\)", r[5]).groups()
    a, b = lenient.get(r[0], (0, 0))
    lenient[r[0]] = (a + round(float(share) * int(n)), b + int(n))
lenient_pct = sorted(round(100 * a / b) for a, b in lenient.values())

# ---------------------------------------------------------------- the evaluated models
flags = {r[0]: r for r in table(read(RUN / "FLAG_REVIEW_3.md"), r"^\| model \| flags drawn \| read")}
per_model = sorted(float(r[6]) for k, r in flags.items() if k != "all")
b1 = dict(re.findall(r"\| the check said (\w+): experts said \| ([^|]+) \|", er))
b3 = dict(re.findall(r"\| the judge said (\w+): experts said \| ([^|]+) \|", er))


def total(s: str) -> int:
    return sum(int(n) for n in re.findall(r"\d+", s))


def q(s: str, label: str) -> str:
    return re.search(rf'"{label}" (\d+)', s).group(1)


OBTAINS, NEVER, UNNEEDED = "yes, the working obtains it", "no, it never obtains it", "no, and its route does not need it"
readings_rows = [
    f"Final answer correct & {total(b1['correct'])} & {b1['correct'].strip()} \\\\",
    f"Final answer incorrect & {total(b1['incorrect'])} & {b1['incorrect'].strip()} \\\\",
    f"Final answer partial & {total(b1['partial'])} & {b1['partial'].strip()} \\\\",
] + [f"Milestone {v.lower()} & {total(b3[v])} & obtains it {q(b3[v], OBTAINS)}, never obtains it {q(b3[v], NEVER)}, "
     f"does not need it {q(b3[v], UNNEEDED)} \\\\" for v in ("REACHED", "MISSING")]

# ---------------------------------------------------------------- the phrases each file must contain
res = read(RUN / "results/RESULTS.md")
judged = [float(r[3]) for r in table(res, r"^\| model \| calls \| without a reply \| judged fraction \|")]
kappa = one(r"between experts, step labels \(Fleiss kappa\) \| \*\*(0\.\d+)\*\* \|.*within an expert, blind re-label "
            r"\(Cohen kappa\) \| \*\*(0\.\d+)\*\* \|.*between experts, milestone status \| (0\.\d+) \|.*"
            r"between experts, final-answer verdict \| (0\.\d+) \|", ps)
design = one(r"(\d+) problems drawn from five engineering branches and three difficulty levels: (\d+) templates", ps)
clean = one(r"From the (\d+) traces the experts called clean, (\d+) each receive exactly one defect and (\d+) are kept", ps)
flawed = one(r"incorrect step in (\d+) of the (\d+) correct-answer traces \((\d+) steps\)", ps)
slips = one(r"the experts found three such steps against (\d+) calculation slips", ps)
gold = one(r"\| digit-rule flags \| 0 of (\d+) claims in (\d+) steps \|", sv)
settled = one(r"\| milestones: by E3, judged \| (\d+) of (\d+), (\d+) \(published", read(RUN / "E5_VALIDATION.md"))
fr = flags["all"]


def b1n(verdict: str, label: str) -> str:
    return re.search(rf"{label} (\d+)", b1[verdict]).group(1)


def thousands(n: str) -> str:
    return f"{int(n):,}"


psb = ps.replace("**", "").replace("*", "")
scalar = one(r"(\d+) of its (\d+) items \((\d+)%\) have a single scalar answer, against (\d+) of the benchmark's "
             r"(\d+) templates", ps)
rounds = one(r"All ([\d,]+) of them\..*re-labelled (\w+) of their own traces.*?(\d+) split steps across (\d+) traces", ps)
by_branch = one(r"between-expert step kappa runs from (0\.\d+) \(civil\) to (0\.\d+) \(industrial\)", ps)
adjudicated = one(r"changed (\d+) step labels, raising the count of steps called incorrect from (\d+) to (\d+)", ps)
power = one(r"design effect of ([\d.]+) to ([\d.]+) depending.*?roughly (\d+) to (\d+) independent.*?"
            r"can detect is (0\.\d+) to (0\.\d+)", ps)
decides = one(r"AUROC (0\.\d+), against.*?Of (\d+) traces with a correct final answer, their holistic verdict calls "
              r"only (three) unsound", psb)
prm_fit = one(r"changes the held-out F1 by ([−-]0\.\d+)", psb)
outside = one(r"moved the pooled score from (0\.\d+) to (0\.\d+)", ps)
routing = one(r"\| the flawed step is shown to a judge \| (0\.\d+) \|.*?\| the same step, unmodified, is shown \| "
              r"(0\.\d+) \|.*?caught by either of its two judges \| (0\.\d+) \|", psb)
batched = one(r"it sends (\d+) of the (\d+) conceptual defects to the judge.*?catches (\d+): 0\.317 end to end.*?"
              r"raised (\d+) false alarms on the (\d+) clean steps", ps)
mimo_n = one(r"MiMo returned a verdict on both arms for (\d+) of the 60 conceptual defects and for (\d+) of the 120", ps)
swap_n = one(r"(\d+) traces per model over (\d+) to (\d+) templates", read(RUN / "JUDGE_SWAP.md"))
sampled = one(r"by the check's verdict: correct (\d+), incorrect (\d+), partial (\d+)", er)
ruled = one(r"by the judge's verdict: MISSING (\d+), REACHED (\d+)", er)
pairs = one(r"items read by two experts; their agreement; Cohen's kappa \| 150; (0\.\d+); (0\.\d+) \|.*"
            r"items read by two experts; their agreement; Cohen's kappa \| 100; (0\.\d+); (0\.\d+) \|", er)
gold_ok = one(r"\| scored correct at all three tolerances \| (\d+) of (\d+) \|", sv)
no_ms = one(r"\| every milestone found \| (\d+); (\d+) items have no milestones \|", sv)
sep = [float(r.split(" & ")[4].rstrip(" \\")) for r in grid_rows if " & yes & " in r and not r.startswith("2")]
caught = [int(re.match(r"\\texttt\{[^}]+\} & (\d+) of", r).group(1)) for r in planted_rows]
assert gold_ok[0] == gold_ok[1] and decides[2] == "three" and sampled[0] == "75", (gold_ok, decides, sampled)

experts = one(r"(\d+) domain experts, three per branch", ps)
assert max(rates) < 0.4, rates  # "at most about a third" in the main text
main = [
    f"on 300 responses from five LLMs outside the eleven we evaluate, which {experts[0]} domain experts",
    f"three-way verdict on {now['answer, three-way agreement']} of the responses",
    f"reaches F1 {e5[2]} against their milestone labels",
    f"{pct(min(judged))} to {pct(max(judged))} of milestones",
]
scoring = [
    f"On the {thousands(gold_ok[1])} gold traces, the check scores every instance correct",
    f"the {no_ms[1]} instances without a milestone",
    f"it is {held[0][1]}, and on the other half {held[1][1]}", f"{held[0][3]} and {held[1][3]}",
    f"{min(sep):.3f} to {max(sep):.3f} between 0.2\\% and 1\\%",
    f"{settled[0]} of the {thousands(settled[1])} milestones", f"the judge decides {settled[2]}",
    f"{pct(min(judged))} to {pct(max(judged))} of milestones",
    f"rules {fv[2]} of the 88 not needed", f"rules {tv[1]} of the 88 true values missing",
    f"none of the {thousands(gold[0])} calculations in their {thousands(gold[1])} steps",
    f"precision {rn[2]} and recall {rn[3]} over all steps, and {rc[2]} and {rc[3]}",
    f"catches {batched[2]} of {batched[1]}", f"{batched[3]} false alarms on the {batched[4]} untouched steps",
] + judge_rows + grid_rows + arith_rows
validation = [
    f"Only {scalar[0]} of the {scalar[1]} instances ({scalar[2]}\\%) have a single scalar answer, against {scalar[3]} of "
    f"the {scalar[4]} templates",
    f"every one of the {rounds[0]} steps", f"labels {rounds[1]} of their own responses",
    f"{rounds[2]} steps in {rounds[3]} responses",
    f"is {kappa[0]} on steps (from {by_branch[0]} in civil to {by_branch[1]} in industrial engineering), {kappa[2]} on "
    f"milestones, and {kappa[3]} on final answers",
    f"between the two passes is {kappa[1]}",
    f"changes {adjudicated[0]} step labels", f"from {adjudicated[1]} to {adjudicated[2]}",
    f"design effect of {power[0]} to {power[1]}", f"about {power[2]} to {power[3]} independent responses",
    f"is {power[4]} to {power[5]}",
    f"at AUROC {decides[0]}", f"only 3 of the {decides[1]} correct-answer responses",
    f"{flawed[0]} of the {flawed[1]} correct-answer responses, {flawed[2]} steps", f"{slips[0]} of the {flawed[2]}",
    f"(AUROC {prm[0]})", f"precision is {prm[1]} over all steps and {prm[2]} inside",
    f"other half by ${prm_fit[0].replace('−', '-')}$", f"from {outside[0]} to {outside[1]}",
    f"{clean[0]} responses", f"in each of {clean[1]} of them", f"keep {clean[2]} untouched",
    f"find at most {max(now['digit rule, all traces, recall'], now['digit rule, hard case, recall'], rn[3], rc[3])} of "
    f"the steps that the experts mark incorrect",
    f"detects 0 of the {planted['conceptual defects'][0].split(' of ')[1]} conceptual defects",
    f"catch {pct(min(rates))} to {pct(max(rates))} of the conceptual defects without false alarms",
    f"{slips[0]} of the {flawed[2]} incorrect steps that the experts find",
    f"more than {max(caught)} of the 60 conceptual defects",
    f"untouched ({routing[1]})" if routing[0] == routing[1] else "ROUTING CHANGED",
    f"catches {routing[2]} of them end to end", f"sends {batched[0]} of the {batched[1]}", f"catches {batched[2]}.",
    f"both versions of {mimo_n[0]}", f"{swap_n[0]} responses per model, drawn from {swap_n[1]} to {swap_n[2]}",
    f"{signed(f'{min(diffs):+.3f}')} to {signed(f'{max(diffs):+.3f}')}",
    f"by at most {max(family):.3f}", f"{lenient_pct[0]}\\% to {lenient_pct[-1]}\\%",
    f"of {fr[1]} flags, {fr[3]} are slips, {fr[4]} are checker errors, and {fr[5]} is unsure",
    f"{fr[3]} of the {int(fr[3]) + int(fr[4])}", f"precision of {fr[6]}", fr[7],
    f"{per_model[0]:.3f} to {per_model[-1]:.3f}",
    f"({sampled[0]} correct, {sampled[1]} incorrect, and {sampled[2]} partial)",
    f"({ruled[1]} reached and {ruled[0]} missing)",
    f"{pct(float(pairs[0]))} of the final answers (Cohen's~$\\kappa$ {pairs[1]})",
    f"{pct(float(pairs[2]))} of the milestones ($\\kappa$ {pairs[3]})",
    f"({b1n('partial', 'correct')} of {total(b1['partial'])} readings)",
    f"in {q(b3['MISSING'], UNNEEDED)} of {total(b3['MISSING'])} readings",
] + comparison_rows + planted_rows + swap_rows + readings_rows

if "--check" in sys.argv:
    bad = 0
    for path, wanted in ((MAIN, main), (SCORING, scoring), (VALIDATION, validation)):
        tex = " ".join(read(path).split()) if path.exists() else ""
        missing = [s for s in wanted if " ".join(s.split()) not in tex]
        for s in missing:
            print(f"MISSING in {path.name}: {s}")
        print(f"{len(wanted) - len(missing)} of {len(wanted)} generated rows and phrases are in {path.name}")
        bad += len(missing)
    sys.exit(1 if bad else 0)

for title, rows in (("components against the experts", comparison_rows), ("planted defects", planted_rows),
                    ("judge validation", judge_rows), ("milestone tolerance", grid_rows),
                    ("arithmetic readings", arith_rows), ("judge swap", swap_rows), ("readings", readings_rows)):
    print(f"% {title}")
    print("\n".join(rows))
for title, ph in (("5_evaluation.tex", main), ("scoring.tex", scoring), ("validation.tex", validation)):
    print(f"% phrases {title} must contain")
    print("\n".join(p for p in ph if " & " not in p))
