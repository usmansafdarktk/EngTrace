"""Print the numbers of the paper's experimental setup (Section 5.1) and check them in its source.

    python full_run_28092026/paper_setup.py           # print every phrase the setup uses
    python full_run_28092026/paper_setup.py --check   # exit 1 unless overleaf_source_04102026/6_experiments.tex holds
                                                      # every phrase and every number in its prose is a phrase's

Sources, each written by a committed script: results/decoding_table.json (decoding_table.py: rows per model, output
ceilings, the sampling parameters sent, the quantizations admitted, finish reasons, reasoning tokens),
results/results.json (analyze.py: verdict counts per model, the level-gap spreads, the anchors and the subset they ran
on, and the bootstrap and sign-flip counts it ran with), manifest.jsonl (the evaluation set) and analyze.py itself
(the interval level, the power and the significance level). Model display names are those docs/appendix_evaluation.py
writes in the appendices, plus the two anchors.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from decimal import ROUND_HALF_UP, Decimal
from math import comb
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
TEX = HERE.parent / "overleaf_source_04102026" / "6_experiments.tex"

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
BIB = HERE.parent / "overleaf_source_04102026" / "custom.bib"
WORD = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
        11: "eleven", 12: "twelve"}
BITS = {"fp8": 8, "mxfp8": 8, "fp16": 16, "bf16": 16, "fp32": 32}


def tt(key: str) -> str:
    return rf"\texttt{{{NAME[key]}}}"


def series(keys: list[str], cite: bool = False) -> str:
    names = [tt(k) + (rf"~\citep{{{CITE[k]}}}" if cite else "") for k in keys]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + ", and " + names[-1]


def alpha(keys) -> list[str]:
    return sorted(keys, key=lambda k: NAME[k].lower())


def thousands(n: int) -> str:
    return f"{n:,}"


def half_up(x: float) -> int:
    return int(Decimal(str(x)).quantize(Decimal(1), rounding=ROUND_HALF_UP))


res = json.loads((HERE / "results/results.json").read_text(encoding="utf-8"))
decoding = {d["model_key"]: d for d in json.loads((HERE / "results/decoding_table.json").read_text(encoding="utf-8"))}
manifest = [json.loads(l) for l in (HERE / "manifest.jsonl").read_text(encoding="utf-8").splitlines() if l.strip()]
analyze = (HERE / "analyze.py").read_text(encoding="utf-8")

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
rows = sum(d["rows"] for d in dec.values())
assert all(d["rows"] == len(manifest) for d in dec.values()), "one response per instance and model"

# ---------------------------------------------------------------- decoding
assert all(d["sampling_parameters_sent"] == [] for d in dec.values()), "no sampling parameter is sent"
ceilings = Counter(d["max_tokens"] for d in dec.values())
ceiling = ceilings.most_common(1)[0][0]
lower = {k: d["max_tokens"] for k, d in dec.items() if d["max_tokens"] != ceiling}
assert list(lower) == ["muse-glimmer-30b"], lower
quants = {q for k in open_ for q in dec[k]["quantizations"]}
assert all(set(dec[k]["quantizations"]) == quants for k in open_), "one quantization rule for every open model"
min_bits = min(BITS[q] for q in quants)
no_reasoning = [k for k in roster if dec[k]["reasoning_tokens"]["share_above_zero"] == 0.0]
assert all(dec[k]["reasoning_tokens"]["max"] == 0 for k in no_reasoning)
no_reasoning = alpha(k for k in open_ if k in no_reasoning) + alpha(k for k in closed if k in no_reasoning)
reasoning = [k for k in roster if k not in no_reasoning]
medians = [dec[k]["reasoning_tokens"]["median"] for k in reasoning]
closed_without = sum(k in closed for k in no_reasoning)

# ---------------------------------------------------------------- responses without a readable final answer
unusable = sum(m["unusable"] for m in q1.values())
empty = sum(m["empty"] for m in q1.values())
# A response that stops at the ceiling is either scored on the text it has (counted as capped) or empty.
at_ceiling = {k: dec[k]["finish_reasons"].get("length", 0) - q1[k]["capped_scored"] for k in roster}
assert all(0 <= at_ceiling[k] <= q1[k]["empty"] for k in roster), at_ceiling
assert sum(at_ceiling.values()) >= 0.95 * empty, "almost all empty responses are at the ceiling"

# ---------------------------------------------------------------- the anchors and their subset
anchors = sorted({a["model"] for arm in res["anchors"] for a in arm["anchors"]}, key=lambda k: k != "gpt-5.4")
subset = {arm["items"] for arm in res["anchors"]} | {a["items"] for a in res["reasoning_arms"] if a["arm"] == "tool"}
assert len(subset) == 1, subset
subset = subset.pop()
subset_templates = {a["templates"] for a in res["reasoning_arms"] if a["items"] == subset}
assert subset_templates == {len(per_template)} and subset % len(per_template) == 0

# ---------------------------------------------------------------- the statistical protocol
B, B_TEST = res["provenance"]["analyze"]["B"], res["provenance"]["analyze"]["B_TEST"]
lo, hi = map(float, re.search(r"nanpercentile\(a, ([\d.]+)\)\), float\(np\.nanpercentile\(a, ([\d.]+)\)\)",
                              analyze).groups())
level = round(hi - lo)
power = int(re.search(r"detects at (\d+)% power", analyze).group(1))
sig = float(re.search(r"p\['p_holm'\] < (0\.\d+)", analyze).group(1))
assert len(res["q1"]["pairs"]) == comb(len(roster), 2)
assert all(m["sd_advanced"] > m["sd_easy"] for m in res["q2"]), "Advanced templates vary more than Easy ones"

phrases = [
    f"We evaluate {WORD[len(roster)]} LLMs.",
    f"The {WORD[len(open_)]} open-weights models are {series(open_, cite=True)}.",
    f"The {WORD[len(closed)]} closed models are {series(closed, cite=True)}.",
    f"{WORD[len(anchors)]} flagships that the rule admits, {series(anchors, cite=True).replace(', and ', ' and ')}, "
    f"serve as anchors on a fixed subset of {subset} instances "
    f"({WORD[subset // len(per_template)]} per template) that the further conditions below also use",
    f"Each model answers each instance once, {thousands(rows)} responses in all, at its provider's default decoding "
    "settings: we set no sampling parameters",
    f"The output ceiling is {thousands(ceiling)} tokens ({thousands(lower['muse-glimmer-30b'])} for "
    f"{tt('muse-glimmer-30b')}, the most its endpoint allows)",
    f"every open-weights model runs at {min_bits}-bit floating-point precision or higher",
    f"{WORD[len(no_reasoning)]} models return no reasoning tokens ({series(no_reasoning)}, "
    f"{WORD[closed_without]} of the {WORD[len(closed)]} closed models), while the other {WORD[len(reasoning)]} write "
    f"a median of {thousands(half_up(min(medians)))} to {thousands(half_up(max(medians)))} reasoning tokens per "
    "response",
    f"{unusable} responses ({100 * unusable / rows:.1f}\\%) contain no readable final answer and score 0; {empty} of "
    "them are empty, almost all at the output ceiling",
    f"the {k_inst} instances of a template share one derivation",
    f"Every {level}\\% interval is a percentile bootstrap over templates with {thousands(B)} resamples",
    f"using {thousands(B_TEST)} sign flips",
    "because scores vary more across Advanced templates",
    f"such as the {comb(len(roster), 2)} pairs of models",
    f"keeps the probability of any false positive at {round(100 * sig)}\\% or less",
    f"We report every non-significant difference with the smallest difference the design detects at {power}\\% "
    "power",
]

NUMBER = re.compile(r"\d[\d,]*(?:\.\d+)?")
NUMBER_WORDS = re.compile(r"\b(" + "|".join(WORD.values()) + r")\b", re.I)


def prose(tex: str) -> str:
    """The text with comments, model names and reference or citation keys removed, on one line."""
    tex = re.sub(r"(?<!\\)%.*", "", tex)
    tex = re.sub(r"\\(texttt|label|autoref|citep|citet|ref)\{[^}]*\}", " ", tex)
    return re.sub(r"\s+", " ", tex)


def numbers(text: str) -> Counter:
    return Counter(n.rstrip(",") for n in NUMBER.findall(text))


if "--check" in sys.argv:
    tex = TEX.read_text(encoding="utf-8")
    flat = re.sub(r"\s+", " ", re.sub(r"(?<!\\)%.*", "", tex))
    missing = [p for p in phrases if re.sub(r"\s+", " ", p) not in flat]
    known_numbers = set(numbers(prose(" ".join(phrases))))
    known_words = {w.lower() for w in NUMBER_WORDS.findall(prose(" ".join(phrases)))}
    stray = sorted(set(numbers(prose(tex))) - known_numbers)
    stray_words = sorted({w.lower() for w in NUMBER_WORDS.findall(prose(tex))} - known_words)
    bib_keys = set(re.findall(r"@\w+\{([^,\s]+),", BIB.read_text(encoding="utf-8")))
    cited = {k.strip() for c in re.findall(r"\\cite[pt]\{([^}]*)\}", tex) for k in c.split(",")}
    unresolved = sorted(cited - bib_keys)
    for k in unresolved:
        print(f"UNRESOLVED citation: {k}")
    for p in missing:
        print(f"MISSING in {TEX.name}: {p}")
    for n in stray + stray_words:
        print(f"NOT GENERATED: {n}")
    print(f"{len(phrases) - len(missing)} of {len(phrases)} phrases present; "
          f"{len(stray) + len(stray_words)} numbers not generated")
    print(f"{len(cited) - len(unresolved)} of {len(cited)} citation keys resolve in {BIB.name}")
    sys.exit(1 if missing or stray or stray_words or unresolved else 0)

for p in phrases:
    print(p)
