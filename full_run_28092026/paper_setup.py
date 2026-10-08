"""Print the numbers of the paper's experimental setup (Section 5.1) and its appendix, and check them in the source.

    python full_run_28092026/paper_setup.py           # print the phrases and the appendix's table rows
    python full_run_28092026/paper_setup.py --check   # exit 1 unless 6_experiments.tex and appendices/models.tex
                                                      # hold every phrase and row, every number in their prose is
                                                      # a phrase's, the prompt box is the prompt the run sent, and
                                                      # every citation key resolves

Sources, each written by a committed script: results/decoding_table.json and decoding_table_flagship.json
(decoding_table.py: rows per model, served identifiers, providers, output ceilings, the sampling parameters sent, the
quantizations admitted, finish reasons, token use, cost, prompt hashes), results/results.json (analyze.py: verdict
counts per model, the level-gap spreads, the anchors and the subset they ran on, the conditions, the bootstrap and
sign-flip counts it ran with), manifest.jsonl (the evaluation set), analyze.py itself (the interval levels, the power,
the significance level, the detectable-difference factor, the equivalence margin) and evaluation/run_inference.py (the
prompt template, matched to the runs by its SHA-256). CARD holds what only the developers publish (organization,
size, weights repository), typed from each model's card and checked there on 2026-10-05; the excluded models come
from docs/inference_pricing/build_pricing_doc.js and the validation appendix. Model display names are those
docs/appendix_evaluation.py writes in the appendices, plus the two anchors.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter
from math import comb
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SRC = REPO / "overleaf_source_04102026"
TEX, APPX, BIB = SRC / "6_experiments.tex", SRC / "appendices" / "models.tex", SRC / "custom.bib"

NAME = {
    "deepseek-v4.1-flash": "DeepSeek V4.1 Flash", "kimi-k3": "Kimi K3", "claude-sonnet-5": "Claude Sonnet 5",
    "glm-5.3-flash": "GLM-5.3-Flash", "muse-glimmer-30b": "Muse Glimmer 30B", "glm-5.3": "GLM-5.3",
    "qwen3-235b-a22b-2507": "Qwen3-235B-2507", "gemini-3.1-flash-lite": "Gemini 3.1 Flash-Lite",
    "gemma-4-26b-a4b": "Gemma 4 26B", "gpt-5.4-mini": "GPT-5.4 mini", "gpt-oss-20b": "gpt-oss-20b",
    "gpt-5.4": "GPT-5.4", "deepseek-v4-pro": "DeepSeek V4 Pro",
}
CITE = {  # each model's source in overleaf_source_04102026/custom.bib, checked by hand against the developer's page
    "deepseek-v4.1-flash": "deepseekv41", "gemma-4-26b-a4b": "gemma4", "glm-5.3": "glm5", "glm-5.3-flash": "glm5",
    "gpt-oss-20b": "gptoss", "kimi-k3": "kimik3", "muse-glimmer-30b": "museglimmer", "qwen3-235b-a22b-2507": "qwen3",
    "claude-sonnet-5": "claudesonnet5", "gemini-3.1-flash-lite": "gemini31flashlite", "gpt-5.4-mini": "gpt54",
    "gpt-5.4": "gpt54", "deepseek-v4-pro": "deepseekv4",
}
# Organization, size and weights repository from each developer's Hugging Face card or model page, 2026-10-05.
# Size is the stated total, with the parameters active per token for mixture-of-experts models; closed models
# publish none. GLM-5.3's card states no active count, so its size is the card's count of the released weights.
CARD = {
    "deepseek-v4.1-flash": ("DeepSeek", "552B (8B/16B active)", "deepseek-ai/DeepSeek-V4.1-Flash"),
    "gemma-4-26b-a4b": ("Google DeepMind", "25.2B (3.8B active)", "google/gemma-4-26B-A4B-it"),
    "glm-5.3": ("Z.ai", "753B", "zai-org/GLM-5.3"),
    "glm-5.3-flash": ("Z.ai", "320B (18B active)", "zai-org/GLM-5.3-Flash"),
    "gpt-oss-20b": ("OpenAI", "21B (3.6B active)", "openai/gpt-oss-20b"),
    "kimi-k3": ("Moonshot AI", "2.8T (104B active)", "moonshotai/Kimi-K3"),
    "muse-glimmer-30b": ("Meta", "30B", "meta-models/Muse-Glimmer-30B"),
    "qwen3-235b-a22b-2507": ("Alibaba Cloud", "235B (22B active)", "Qwen/Qwen3-235B-A22B-Instruct-2507"),
    "claude-sonnet-5": ("Anthropic", "N/A", None),
    "gemini-3.1-flash-lite": ("Google DeepMind", "N/A", None),
    "gpt-5.4-mini": ("OpenAI", "N/A", None),
    "deepseek-v4-pro": ("DeepSeek", "1.6T (49B active)", "deepseek-ai/DeepSeek-V4-Pro"),
    "gpt-5.4": ("OpenAI", "N/A", None),
}
PREFILL_DECODE = ("8B", "16B")  # DeepSeek V4.1 Flash's card: active parameters per token in prefill and in decoding
NO_THINKING = "qwen3-235b-a22b-2507"  # its card: "supports only non-thinking mode"
# The selection rule's exclusions (build_pricing_doc.js): the expert study's five models (validation appendix), the
# robustness pair the pilot generated with, and the models used as judges, named as the appendices name them.
ROBUSTNESS = {"gemma-4-31b-it": "Gemma 4 31B", "qwen3.8-27b": "Qwen3.8-27B"}
JUDGES = {"gpt-5": "GPT-5", "claude-opus-4.5": "Claude Opus 4.5", "gemini-3.1-pro": "Gemini 3.1 Pro",
          "grok": "Grok 4.6", "minimax": "MiniMax M3", "mimo": "MiMo-V2.5-Pro"}
CONDITIONS = {"reasoning-medium": "reasoning effort", "flagship": "flagship anchors",
              "openbook2": "governing equations", "tool": "a Python tool"}
WORD = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
        11: "eleven", 12: "twelve"}
BITS = {"fp8": 8, "mxfp8": 8, "fp16": 16, "bf16": 16, "fp32": 32}


def tt(key: str) -> str:
    return rf"\texttt{{{NAME[key]}}}"


def ttn(name: str) -> str:
    return rf"\texttt{{{name}}}"


def listing(names: list[str]) -> str:
    return names[0] if len(names) == 1 else (" and ".join(names) if len(names) == 2 else
                                             ", ".join(names[:-1]) + ", and " + names[-1])


def series(keys: list[str], cite: bool = False) -> str:
    return listing([tt(k) + (rf"~\citep{{{CITE[k]}}}" if cite else "") for k in keys])


def alpha(keys) -> list[str]:
    return sorted(keys, key=lambda k: NAME[k].lower())


def thousands(n: int) -> str:
    return f"{n:,}"


def tokens(x: float) -> str:
    return thousands(round(x))  # half to even, as decoding_table.py prints them


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


res = load(HERE / "results/results.json")
decoding = {d["model_key"]: d for d in load(HERE / "results/decoding_table.json")}
flagship = {d["model_key"]: d for d in load(HERE / "results/decoding_table_flagship.json")}
manifest = [json.loads(l) for l in (HERE / "manifest.jsonl").read_text(encoding="utf-8").splitlines() if l.strip()]
analyze = (HERE / "analyze.py").read_text(encoding="utf-8")
runner = (REPO / "evaluation/run_inference.py").read_text(encoding="utf-8")
pricing = (REPO / "docs/inference_pricing/build_pricing_doc.js").read_text(encoding="utf-8")
validation = (SRC / "appendices/validation.tex").read_text(encoding="utf-8")

q1 = {m["model"]: m for m in res["q1"]["models"]}
roster = list(q1)
dec = {k: decoding[k] for k in roster}
open_ = alpha(k for k in roster if dec[k]["weights"] == "open")
closed = alpha(k for k in roster if dec[k]["weights"] == "closed")
assert len(open_) + len(closed) == len(roster), "every evaluated model is open or closed"

# ---------------------------------------------------------------- the evaluation set and the responses
per_template = Counter(r["template_id"] for r in manifest)
instances = set(per_template.values())
assert len(instances) == 1, "every template has the same number of instances"
k_inst = instances.pop()
levels = Counter(lv for t, lv in {(r["template_id"], r["level"]) for r in manifest})
rows = sum(d["rows"] for d in dec.values())
assert all(d["rows"] == len(manifest) for d in dec.values()), "one response per instance and model"

# ---------------------------------------------------------------- decoding
assert all(d["sampling_parameters_sent"] == [] for d in dec.values()), "no sampling parameter is sent"
ceilings = Counter(d["max_tokens"] for d in dec.values())
ceiling = ceilings.most_common(1)[0][0]
lower = {k: d["max_tokens"] for k, d in dec.items() if d["max_tokens"] != ceiling}
assert list(lower) == ["muse-glimmer-30b"] and len(dec["muse-glimmer-30b"]["providers"]) == 1, lower
quants = {q for k in open_ for q in dec[k]["quantizations"]}
assert all(set(dec[k]["quantizations"]) == quants for k in open_), "one quantization rule for every open model"
assert all(d["provider_sort"] == "price" and d["allow_fallbacks"] for d in dec.values()), "cheapest, with fallbacks"
min_bits = min(BITS[q] for q in quants)
no_reasoning = [k for k in roster if dec[k]["reasoning_tokens"]["share_above_zero"] == 0.0]
assert all(dec[k]["reasoning_tokens"]["max"] == 0 for k in no_reasoning)
no_reasoning = alpha(k for k in open_ if k in no_reasoning) + alpha(k for k in closed if k in no_reasoning)
reasoning = [k for k in roster if k not in no_reasoning]
medians = [dec[k]["reasoning_tokens"]["median"] for k in reasoning]
closed_without = sum(k in closed for k in no_reasoning)
assert NO_THINKING in no_reasoning
assert min(dec[k]["reasoning_tokens"]["share_above_zero"] for k in reasoning) >= 0.99, "nearly every response"

# ---------------------------------------------------------------- the prompt
template = re.search(r'PROMPT_TEMPLATE = """(.*?)"""', runner, re.S).group(1)
hashes = {h for d in dec.values() for h in d["prompt_sha256"]}
assert hashes == {hashlib.sha256(template.encode("utf-8")).hexdigest()}, "every response used this template"

# ---------------------------------------------------------------- the matched configuration (results/matched_config.json)
# The models that return no reasoning tokens at their providers' defaults and offer a reasoning setting answer every
# instance again with it, on the main run's items, prompt and ceiling (run_traces.py's reasoning-medium-full); the
# model whose endpoint offers none is the one left at its default.
matched_cfg = load(HERE / "results/matched_config.json")
rerun = alpha(m["model"] for m in matched_cfg["models"] if m.get("run", True) and m.get("reasoning_store"))
rstores = {m["reasoning_store"] for m in matched_cfg["models"] if m.get("run", True) and m.get("reasoning_store")}
assert len(rstores) == 1, rstores
dec_r = {d["model_key"]: d for d in load(HERE / f"results/decoding_table_{rstores.pop()}.json")}
assert set(dec_r) == set(rerun) and set(no_reasoning) - set(rerun) == {NO_THINKING}, (rerun, no_reasoning)
assert [m["model"] for m in matched_cfg["models"] if m.get("reasoning_setting") == "none offered"] == [NO_THINKING]
assert all(d["rows"] == len(manifest) and d["max_tokens"] == ceiling and set(d["prompt_sha256"]) == hashes
           and d["reasoning_tokens"]["share_above_zero"] == 1.0 for d in dec_r.values()), "same items, prompt, ceiling"
efforts = {d["reasoning_parameter"]["effort"] for d in dec_r.values()}
assert len(efforts) == 1, efforts
effort = efforts.pop()
rerun_rows = sum(d["rows"] for d in dec_r.values())
rerun_empty = {k: d["empty"] for k, d in dec_r.items() if d["empty"]}
assert all(n <= dec_r[k]["finish_reasons"].get("length", 0) for k, n in rerun_empty.items()), "empty only at the ceiling"

# ---------------------------------------------------------------- responses without a readable final answer
unusable = sum(m["unusable"] for m in q1.values())
empty = sum(m["empty"] for m in q1.values())
# A response that stops at the ceiling is either scored on the text it has (counted as capped) or empty.
at_ceiling = {k: dec[k]["finish_reasons"].get("length", 0) - q1[k]["capped_scored"] for k in roster}
assert all(0 <= at_ceiling[k] <= q1[k]["empty"] for k in roster), at_ceiling
assert sum(at_ceiling.values()) >= 0.95 * unusable, "almost all responses without an answer are empty at the ceiling"

# ---------------------------------------------------------------- the anchors, their subset, the conditions
anchors = sorted({a["model"] for arm in res["anchors"] for a in arm["anchors"]}, key=lambda k: k != "gpt-5.4")
assert set(anchors) == set(flagship)
subset = {arm["items"] for arm in res["anchors"]} | {a["items"] for a in res["reasoning_arms"] if a["arm"] == "tool"}
assert len(subset) == 1, subset
subset = subset.pop()
subset_templates = {a["templates"] for a in res["reasoning_arms"] if a["items"] == subset}
assert subset_templates == {len(per_template)} and subset % len(per_template) == 0
arms = {a["arm"] for a in res["reasoning_arms"]} | {arm["arm"] for arm in res["anchors"]}
assert set(CONDITIONS) <= arms, set(CONDITIONS) - arms
served = {k: d["served"] for k, d in list(dec.items()) + list(flagship.items())}
assert all(len(s) == 1 for s in served.values()), "one served identifier per model"
served = {k: next(iter(s)) for k, s in served.items()}

# ---------------------------------------------------------------- the excluded models
excluded_note = re.search(r"Excluded: models the evaluator pilot generated with\s*\((.*?)\) and models used as\s*"
                          r"judges\s*\((.*?)\)", pricing.replace("\n//", " "), re.S)
generated, judged = [re.sub(r"\s+", " ", g).replace(" and ", ", ") for g in excluded_note.groups()]
assert all(k in generated for k in ROBUSTNESS) and [j.strip() for j in judged.split(",")] == list(JUDGES), judged
study = re.findall(r"\\texttt\{([^}]*)\}", re.search(r"Five LLMs outside the evaluated models[\s(,]+(.*?)[\s),]*answer (?:every|each) "
                                                     r"instance", validation, re.S).group(1))
assert len(study) == 5

# ---------------------------------------------------------------- the statistical protocol
B, B_TEST = res["provenance"]["analyze"]["B"], res["provenance"]["analyze"]["B_TEST"]
lo, hi = map(float, re.search(r"nanpercentile\(a, ([\d.]+)\)\), float\(np\.nanpercentile\(a, ([\d.]+)\)\)",
                              analyze).groups())
level = round(hi - lo)
lo90, hi90 = map(float, re.search(r"nanpercentile\(draws, ([\d.]+)\)\), float\(np\.nanpercentile\(draws, "
                                  r"([\d.]+)\)\)", analyze).groups())
level90 = round(hi90 - lo90)
power = int(re.search(r"detects at (\d+)% power", analyze).group(1))
factor = re.search(r"detects at \d+% power: ([\d.]+) x SD", analyze).group(1)
sig = float(re.search(r"p\['p_holm'\] < (0\.\d+)", analyze).group(1))
margin = float(re.search(r"EQUIV_MARGIN = ([\d.]+)", analyze).group(1))
assert margin == res["q5"]["margin"]
assert len(res["q1"]["pairs"]) == comb(len(roster), 2) == len(res["q3_coverage"]["pairs"])
CLAIMS_FAILED = []  # a prose claim the data no longer supports: reported, counted by --check, generation goes on
if not all(m["sd_advanced"] > m["sd_easy"] for m in res["q2"]):
    CLAIMS_FAILED.append("Advanced templates vary more than Easy ones, for every model (not for "
                         + ", ".join(m["model"] for m in res["q2"] if m["sd_advanced"] <= m["sd_easy"]) + ")")
assert all(m["templates"] == len(per_template) for m in res["q2"])
mc_templates = res["q3_coverage"]["templates"]

phrases = [  # Section 5.1
    f"We evaluate {WORD[len(roster)]} LLMs: {WORD[len(open_)]} open-weights models, {series(open_, cite=True)}; and "
    f"{WORD[len(closed)]} closed models, {series(closed, cite=True)}.",
    "None of them wrote responses that we used to develop or validate the evaluator",
    f"{WORD[len(anchors)].capitalize()} flagships, {series(anchors, cite=True)}, run only as anchors on a fixed subset "
    f"of {subset} instances ({WORD[subset // len(per_template)]} per template), which the further conditions also "
    "use",
    f"Each model answers each instance once ({thousands(rows)} responses) with the same zero-shot prompt and no "
    f"tools or retrieval, at its provider's default decoding settings and an output ceiling of {thousands(ceiling)} "
    f"tokens ({thousands(lower['muse-glimmer-30b'])} for {tt('muse-glimmer-30b')})",
    f"{WORD[len(no_reasoning)]} models return no reasoning tokens ({series(no_reasoning)}, "
    f"{WORD[closed_without]} of the {WORD[len(closed)]} closed models), so the models are not compared at equal "
    "reasoning effort",
    f"{unusable} responses ({100 * unusable / rows:.1f}\\%) have no readable final answer and score 0, almost all "
    "of them empty at the output ceiling",
    f"{series(rerun)} also answer every instance with reasoning at {effort} effort ({thousands(rerun_rows)} responses)",
    f"because its {k_inst} instances share one derivation",
    f"every {level}\\% interval is a bootstrap over templates",
]
appendix_phrases = [  # appendices/models.tex
    f"The first condition excludes, among others, the {WORD[len(study)]} models of the expert study "
    f"({listing([ttn(s) for s in study])}); the second excludes every model we use as a judge "
    f"({listing([ttn(s) for s in JUDGES.values()])})",
    f"lists the {WORD[len(roster)]} evaluated models and the {WORD[len(anchors)]} flagship anchors",
    f"The ceiling is {thousands(ceiling)} tokens, except {thousands(lower['muse-glimmer-30b'])} for "
    f"{tt('muse-glimmer-30b')}, whose single eligible provider caps its output there",
    f"Providers serve every open-weights model at {min_bits}-bit floating-point precision or higher",
    f"every model returns one response, possibly empty, for each of the {thousands(len(manifest))} instances",
    f"{tt(NO_THINKING)} has no thinking mode, and {series([k for k in no_reasoning if k != NO_THINKING])} return no "
    f"reasoning tokens at their providers' defaults; the other {WORD[len(reasoning)]} models return them on nearly "
    "every response",
    f"The eleven models of~\\autoref{{sec:experiments}}, grouped by access, and the {WORD[len(anchors)]} flagship "
    f"anchors run on the {subset}-instance subset",
    f"{tt('deepseek-v4.1-flash')} activates {PREFILL_DECODE[0]} in prefill and {PREFILL_DECODE[1]} in decoding",
] + [
    f"with reasoning on, {tt(k)} leaves {n} of its {thousands(dec_r[k]['rows'])} responses empty at the output ceiling"
    for k, n in sorted(rerun_empty.items())
] + [
    f"Per model over its {thousands(len(manifest))} responses",
    f"every template has {k_inst} instances, so Final Answer Accuracy over the {thousands(len(manifest))} instances "
    f"equals the mean of the {len(per_template)} template means",
    f"every {level}\\% interval is a percentile bootstrap that resamples templates {thousands(B)} times, and a "
    f"bounded change also carries its {level90}\\% interval",
    f"a sign-flip permutation test with {thousands(B_TEST)} sign flips",
    f"those of Milestone Coverage, which use the {mc_templates} templates with milestones",
    "because scores vary more across Advanced templates",
    f"its mean on the {levels['Easy']} Easy templates minus its mean on the {levels['Advanced']} Advanced ones",
    f"such as the {comb(len(roster), 2)} pairs of models for each measure or the {WORD[len(roster)]} level gaps",
    f"at {power}\\% power and a two-sided level of {sig:.2f}, a paired comparison detects ${factor}\\,s/\\sqrt{{n}}$",
    f"a level gap ${factor}\\,\\sigma\\sqrt{{1/{levels['Easy']} + 1/{levels['Advanced']}}}$",
    f"a change counts as bounded when its {level90}\\% interval lies within $\\pm {margin:.2f}$, which is two "
    f"one-sided tests at {round(100 * sig)}\\%",
]


def models_row(k: str) -> str:
    org, size, weights = CARD[k]
    return rf"{tt(k)} & {org} & {size} & \texttt{{{served[k]}}} & " + (
        rf"\texttt{{{weights}}}" if weights else "--") + r" \\"


def decoding_row(k: str) -> str:
    d = dec[k]
    return (rf"{tt(k)} & {thousands(d['max_tokens'])} & {len(d['providers'])} & "
            f"{tokens(d['completion_tokens']['median'])} / {tokens(d['completion_tokens']['p90'])} & "
            f"{tokens(d['reasoning_tokens']['median'])} / {tokens(d['reasoning_tokens']['p90'])} & "
            f"{100 * d['reasoning_tokens']['share_above_zero']:.1f}\\% \\\\")


table_rows = [models_row(k) for k in open_ + closed + alpha(anchors)] + [decoding_row(k) for k in open_ + closed]

NUMBER = re.compile(r"\d[\d,]*(?:\.\d+)?")
NUMBER_WORDS = re.compile(r"\b(" + "|".join(WORD.values()) + r")\b", re.I)


def strip_blocks(tex: str) -> str:
    """Comments, tables, the prompt box, layout settings and the hash's name, none of them a reported number."""
    tex = re.sub(r"(?<!\\)%.*", "", tex)
    tex = re.sub(r"\\begin\{(tabular|verbatim)\}.*?\\end\{\1\}", " ", tex, flags=re.S)
    tex = re.sub(r"\\(resizebox|renewcommand)\{[^}]*\}\{[^}]*\}", " ", tex)
    return tex.replace("SHA-256", " ")


def prose(tex: str) -> str:
    """The text without comments, tables, the prompt box, model names, or reference and citation keys."""
    tex = re.sub(r"\\(texttt|label|autoref|citep|citet|ref)\{[^}]*\}", " ", strip_blocks(tex))
    return re.sub(r"\s+", " ", tex)


def numbers(text: str) -> set[str]:
    return {n.rstrip(",") for n in NUMBER.findall(text)}


def flat(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


for c in CLAIMS_FAILED:
    print(f"CLAIM FAILS: {c}")
if "--check" in sys.argv:
    bib_keys = set(re.findall(r"@\w+\{([^,\s]+),", BIB.read_text(encoding="utf-8")))
    failures = len(CLAIMS_FAILED)
    for path, own in ((TEX, phrases), (APPX, appendix_phrases + table_rows)):
        tex = path.read_text(encoding="utf-8")
        body = flat(re.sub(r"(?<!\\)%.*", "", tex))
        missing = [p for p in own if flat(p) not in body]
        known = numbers(prose(" ".join(phrases + appendix_phrases)))
        words = {w.lower() for w in NUMBER_WORDS.findall(prose(" ".join(phrases + appendix_phrases)))}
        stray = sorted(numbers(prose(tex)) - known) + sorted(
            {w.lower() for w in NUMBER_WORDS.findall(prose(tex))} - words)
        cited = {k.strip() for c in re.findall(r"\\cite[pt]\{([^}]*)\}", tex) for k in c.split(",")}
        unresolved = sorted(cited - bib_keys)
        boxes = re.findall(r"\\begin\{verbatim\}\n(.*?)\\end\{verbatim\}", tex, re.S)
        bad_box = path == APPX and (len(boxes) != 1 or flat(boxes[0]) != flat(template))
        for p in missing:
            print(f"MISSING in {path.name}: {p}")
        for n in stray:
            print(f"NOT GENERATED in {path.name}: {n}")
        for k in unresolved:
            print(f"UNRESOLVED citation in {path.name}: {k}")
        if bad_box:
            print(f"PROMPT BOX in {path.name} is not the prompt the run sent")
        print(f"{path.name}: {len(own) - len(missing)} of {len(own)} phrases and rows present; {len(stray)} numbers "
              f"not generated; {len(cited) - len(unresolved)} of {len(cited)} citation keys resolve"
              + ("; the prompt box matches the sent prompt" if path == APPX and not bad_box else ""))
        failures += len(missing) + len(stray) + len(unresolved) + bad_box
    sys.exit(1 if failures else 0)

for p in phrases + appendix_phrases + table_rows:
    print(p)
