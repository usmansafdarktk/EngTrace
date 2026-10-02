#!/usr/bin/env python
"""Write NEW_PAPERS_AUDIT.md: every candidate from the October 2026 search, its panel reviews, and
whether the revised Related Work cites it (main text, appendix) or not, with the reason.

Inputs: candidates/*.json, reviews/candidate/*.json, related_work_keys.txt (the "# Main text" and
"# Appendix" sections), and REASONS below for candidates the text does not cite.
"""
import collections
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TOPIC = {"engineering_web": "Engineering", "engineering_arxiv": "Engineering",
         "physics_science": "Physics and science", "symbolic_contamination": "Generation and contamination",
         "process_supervision": "Process supervision", "llm_judge": "LLM judges"}
REASONS = {
    "ren2026solutionhacking": "working paper; the point it makes (credited answers reached by shortcuts) is made in v3 with GSM8K and ProcessBench",
    "chi2026frontiereng": "agentic design tasks, outside EngTrace's closed-form scope",
    "gerstmayr2026meceng": "multibody-simulation modelling, outside scope",
    "ravishankara2026circuchain": "small circuit diagnostic; CIRCUIT represents circuits",
    "chen2026optengine": "optimisation modelling; ORQA represents operations research",
    "ren2026autolpbench": "LP problem generator; ORQA represents operations research",
    "sun2026educircuithw": "grading of student work, not evaluation of models' reasoning",
    "pers2026handwrittengrading": "grading of student work, not evaluation of models' reasoning",
    "chen2025circuithomework": "grading of student work, not evaluation of models' reasoning",
    "wang2026mechreason": "static QA mined from papers; SoM-1K is cited for mechanical engineering",
    "guan2026supchainbench": "tool orchestration for supply chains, outside scope",
    "li2025wirelessmathbench": "wireless-communications math; TeleMath represents communications",
    "liu2026astromind": "astrodynamics simulation; APBench represents the domain",
    "zhu2025critpt": "research-level physics, outside scope",
    "chung2025tpbench": "research-level theoretical physics, outside scope",
    "phan2026hle": "general frontier suite; SuperGPQA represents general suites",
    "rein2024gpqa": "general graduate science suite; SuperGPQA represents general suites",
    "wang2024mmlupro": "general suite; MMLU and SuperGPQA represent general suites",
    "xu2026physelite": "olympiad physics with process-level evaluation; PhysReason, PRISM-Physics and HiPhO represent it",
    "dai2025physicsarena": "intermediate-step physics evaluation; PhysReason, PRISM-Physics and HiPhO represent it",
    "zhao2026sciencearena": "olympiad grading with LLM rubrics; HiPhO represents it",
    "wang2026frontierscience": "rubric-graded research tasks; outside scope",
    "feng2025physics": "university physics scored on answers; UGPhysics represents it",
    "arora2023jeebench": "pre-engineering exam questions; represented by the cited exam-style benchmarks",
    "tian2024scicode": "scientific coding, outside scope",
    "mirza2025chembench": "chemistry knowledge suite, outside scope",
    "zhang2026matscibench": "materials-science problems scored on answers, outside the five branches",
    "arabov2026rusfinchain": "applies FinChain's design to Russian finance; FinChain is cited",
    "xu2025ugmathbench": "dynamic undergraduate math; functional variants are cited",
    "zhu2023dyval": "graph-generated reasoning problems; functional and symbolic variants are cited",
    "xu2025reimagine": "symbolic variant synthesis; represented by GSM-Symbolic and VeRA",
    "balunovic2025matharena": "time-window contamination evidence; represented by Deng et al. and GSM1k",
    "white2024livebench": "refreshed benchmark, an alternative to generation; not needed for the positioning",
    "wu2025reasoningmemorization": "contamination evidence for one model family; represented by GSM1k",
    "singh2026gsmsem": "semantic variants of GSM8K; represented by GSM-Plus and MATH-Perturb",
    "shrestha2025gsmranges": "number-scale perturbations; represented by GSM-Symbolic",
    "sun2025bdcmitigation": "contamination-mitigation study; Akhtar et al. make the related point for saturation",
    "pandit2025hard2verify": "step-labelled verification on frontier math; ProcessBench represents it",
    "zeng2023mrgsm8k": "meta-reasoning on GSM8K; MR-Ben is cited",
    "jacovi2024reveal": "open-domain QA chains, outside scope",
    "prasad2023receval": "reference-free step metric; ROSCOE and ReasonEval are cited",
    "lee2025stepsurvey": "survey; the specific works are cited",
    "zhang2024genrm": "generative verifier for training; GenPRM is cited",
    "khalifa2025thinkprm": "generative PRM; represented by GenPRM and VersaPRM",
    "lee2025rethinkingrm": "reward-model comparison; VersaPRM makes the cross-domain point",
    "zi2026neurosymbolicprm": "neuro-symbolic PRM; GenPRM represents program-aided step checks",
    "yeadon2026answerargument": "correct-answer physics traces with hidden errors; ProcessBench makes the point",
    "srivatsa2025peek": "first-error location in student solutions; represented by the error-detection benchmarks cited",
    "liu2023geval": "early LLM evaluator; MT-Bench represents the paradigm",
    "wataoka2024selfpreference": "self-preference mechanism; Panickssery et al. represent it",
    "gu2024surveyjudge": "survey; the specific works are cited",
    "chen2025selfpreferencereason": "self-preference on verifiable tasks; Pombal et al. represent it",
    "hossain2026agreementoverstates": "judge-panel dependence; Kim et al. and Goel et al. represent it",
    "gonzalez2026qedbench": "proof grading by judges; Yeadon et al. represent technical-domain grading",
    "lee2026reportjudge": "statistical correction of judge scores; not needed for the positioning",
}


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def main():
    main_keys, app_keys, section = set(), set(), None
    with open(os.path.join(HERE, "related_work_keys.txt"), encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if ln.startswith("# Main"):
                section = main_keys
            elif ln.startswith("# Appendix"):
                section = app_keys
            elif ln and not ln.startswith("#") and section is not None:
                section.add(ln)
    cands, seen = [], set()
    for f in sorted(glob.glob(os.path.join(HERE, "candidates", "*.json"))):
        for e in load(f):
            if e["key"] not in seen:
                seen.add(e["key"])
                cands.append(e)
    reviews = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(HERE, "reviews", "candidate", "*.json")):
        key, rev = os.path.basename(f)[:-5].rsplit(".", 1)
        reviews[key][rev] = load(f).get("recommendation", "")
    rows, counts = [], collections.Counter()
    for e in cands:
        n = e.get("notes", "")
        pri = "must" if "must-cite" in n else ("should" if "should-cite" in n else "optional")
        k = e["key"]
        dec = "main text" if k in main_keys else ("appendix" if k in app_keys else "not cited")
        counts[(dec, pri)] += 1
        rv = reviews.get(k, {})
        panel = "; ".join(f"{r}: {v.split(' -')[0].split(',')[0]}" for r, v in sorted(rv.items())) or "–"
        reason = "" if dec != "not cited" else REASONS.get(k, "not needed for the positioning; the paragraph cites representative work")
        title = e["title"] if len(e["title"]) <= 90 else e["title"][:87] + "..."
        rows.append((TOPIC.get(e.get("found_by"), "?"), dec, pri, k, title, str(e.get("venue") or ""), panel, reason))
    order_dec = {"main text": 0, "appendix": 1, "not cited": 2}
    order_pri = {"must": 0, "should": 1, "optional": 2}
    rows.sort(key=lambda r: (list(TOPIC.values()).index(r[0]) if r[0] in TOPIC.values() else 9, order_dec[r[1]], order_pri[r[2]], r[3]))
    total = len(cands)
    n_main = sum(v for (d, _), v in counts.items() if d == "main text")
    n_app = sum(v for (d, _), v in counts.items() if d == "appendix")
    must_not = [r[3] for r in rows if r[2] == "must" and r[1] == "not cited"]
    out = [
        "# Audit of the papers found in the October 2026 search", "",
        "Generated by `make_new_papers_audit.py`. The search (six topics, `notes/search_<topic>.md`) found",
        f"{total} candidate papers not cited in the May submission, most of them from 2025 and 2026. The revised",
        f"Related Work cites {n_main} in the main text and {n_app} in the appendix. Of the search's must-cite",
        f"candidates, it leaves out {len(must_not)}: {', '.join(must_not) or 'none'}.", "",
        "**How each was assessed.** Every candidate was downloaded and described from its text by the agent that",
        "found it. The panel of `notes/review_protocol.md` then reviewed the must-cite candidates first: three",
        "independent reviews for the eight engineering must-cites and the process-supervision must-cites, and",
        "one or more for part of the LLM-judge batch, before the review was stopped on 2026-10-02 to save cost",
        "(the owner's instruction; recorded in `notes/review_protocol.md`). For the rest, the search notes",
        "were used, and every statement the revised text makes about any of these papers was checked against",
        "the paper's own text by `verify_facts.py` (`notes/fact_check.md`). A paper the text cites therefore",
        "never rests on an unchecked description, whichever route it took. The 18 main-text candidates the",
        "panel had not reached were then reviewed by one reviewer each, under the same protocol, for anything",
        "that would change the positioning; that note lists the 18 and their verdicts",
        "(`notes/followup_review.md`).", "",
        "**Choosing among them.** The main text cites what a reviewer would expect, by priority. That means",
        "every engineering benchmark that generates instances, checks steps or validates its scorer. It also",
        "means the physics and template precedents that do any of these. Last come the process-supervision",
        "and judge-validity results that the evaluator's design answers. The appendix lists further",
        "engineering coverage and secondary results. Work outside the closed-form, text-based setting is not",
        "cited: agentic, design or simulation-modelling tasks, grading of student work, research-level physics",
        "and general suites. Neither is work already represented by a cited paper making the same point.", "",
        "| decision | main text | appendix | not cited |", "|---|---:|---:|---:|",
    ]
    for pri in ("must", "should", "optional"):
        out.append(f"| {pri}-cite (search) | {counts[('main text', pri)]} | {counts[('appendix', pri)]} | {counts[('not cited', pri)]} |")
    current = None
    for topic, dec, pri, k, title, venue, panel, reason in rows:
        if topic != current:
            current = topic
            out += ["", f"## {topic}", "", "| key | title | venue | search | panel | v3 | why not cited |", "|---|---|---|---|---|---|---|"]
        out.append(f"| {k} | {title.replace('|', '/')} | {venue.replace('|', '/')[:40]} | {pri} | {panel} | {dec} | {reason} |")
    with open(os.path.join(HERE, "NEW_PAPERS_AUDIT.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"{total} candidates: {n_main} main text, {n_app} appendix, {total - n_main - n_app} not cited; must-cites left out: {must_not}")


if __name__ == "__main__":
    main()
