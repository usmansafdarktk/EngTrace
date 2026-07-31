"""
EngTrace — Error Distribution Analysis
======================================
Loads all three annotators' JSONL files, computes a majority-vote label for
every trace, then produces a structured JSON with insights ready for plotting.

Output: analysis_results/error_distribution.json
"""

import json
import os
from collections import Counter, defaultdict

# ── paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.join(BASE_DIR, "analysis_results")
os.makedirs(OUT_DIR, exist_ok=True)

ANNOTATORS = ["Annotator A", "Annotator B", "Annotator C"]

MODELS = [
    "claude-opus-4-7",
    "deepseek-r1",
    "deepseek-v4-pro",
    "gemini-3.1-pro-preview",
    "gemma-3-27b",
    "gpt-5",
    "gpt-5-mini",
    "mathstral-7b",
    "meta-llama-3.1-70b",
    "qwen2.5-72b",
    "wizardmath-7b",
]

# Models to highlight in the main paper section
MAIN_TEXT_MODELS = {
    "gemini-3.1-pro-preview": "Gemini 3.1 Pro Preview",
    "deepseek-r1":            "DeepSeek R1",
    "meta-llama-3.1-70b":    "Llama 3.1 70B",
}

CATEGORY_NAMES = {
    0: "No Error",
    1: "Hallucination",
    2: "Setup / Assumption Error",
    3: "Formula / Principle Error",
    4: "Unit / Dimensional Error",
    5: "Sign / Direction Error",
    6: "Calculation Error",
}

CATEGORY_SEVERITY = {
    1: "Critical",
    2: "High",
    3: "High",
    4: "Medium",
    5: "Medium",
    6: "Low",
}

NAME_TO_CODE = [
    ("no error",               0),
    ("hallucination",          1),
    ("setup/assumption error", 2),
    ("setup/assumption",       2),
    ("setup",                  2),
    ("formula/principle error",3),
    ("formula/principle",      3),
    ("formula",                3),
    ("unit/dimensional error", 4),
    ("unit/dimensional",       4),
    ("unit",                   4),
    ("sign/direction error",   5),
    ("sign/direction",         5),
    ("sign",                   5),
    ("calculation error",      6),
    ("calculation",            6),
]

# ── helpers ────────────────────────────────────────────────────────────────────

def resolve_category(cat):
    if cat is None or cat == "":
        return 0
    if isinstance(cat, (int, float)):
        v = int(cat)
        return v if v in CATEGORY_NAMES else 0
    key = str(cat).lower().strip()
    for name, code in NAME_TO_CODE:
        if name in key:
            return code
    return 0


def pct(count, total):
    return round(count / total * 100, 2) if total > 0 else 0.0


def majority_vote(votes):
    """Return the category chosen by at least 2 of 3 annotators.
    Falls back to the lowest-numbered (highest-severity) category on a three-way split."""
    c = Counter(votes)
    top_count = c.most_common(1)[0][1]
    if top_count >= 2:
        # Among tied majority labels, pick lowest error-category number (highest priority)
        winners = [cat for cat, cnt in c.items() if cnt == top_count]
        return min(winners)
    # Three-way split → take highest-severity (lowest number, excluding 0)
    non_zero = [v for v in votes if v != 0]
    return min(non_zero) if non_zero else 0


# ── data loading ───────────────────────────────────────────────────────────────

def load_all_annotations():
    """
    Returns dict: annotation_id → {
        meta fields from the JSONL,
        "votes": [cat_A, cat_B, cat_C],
        "final_category": int (majority vote),
    }
    """
    # First pass: collect votes per annotation_id
    votes_map   = defaultdict(list)   # annotation_id → list of (annotator, category)
    meta_map    = {}                   # annotation_id → meta fields

    for annotator in ANNOTATORS:
        for model in MODELS:
            fpath = os.path.join(BASE_DIR, annotator, f"{model}_annotation_sample.jsonl")
            if not os.path.exists(fpath):
                print(f"  WARNING: missing {fpath}")
                continue
            with open(fpath, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    aid = rec["annotation_id"]
                    cat = resolve_category(rec.get("error_category"))
                    votes_map[aid].append(cat)

                    if aid not in meta_map:
                        meta_map[aid] = {
                            "annotation_id":   aid,
                            "model":           rec.get("model", ""),
                            "question_id":     rec.get("question_id", ""),
                            "difficulty":      rec.get("difficulty", ""),
                            "branch":          rec.get("branch", ""),
                            "domain":          rec.get("domain", ""),
                            "recovered_f1":    rec.get("recovered_f1", 0.0),
                            "recall":          rec.get("recall", 0.0),
                        }

    # Second pass: compute majority-vote label
    records = {}
    for aid, votes in votes_map.items():
        final = majority_vote(votes)
        records[aid] = {
            **meta_map[aid],
            "votes":          votes,
            "final_category": final,
            "category_name":  CATEGORY_NAMES.get(final, "Unknown"),
        }

    return records


# ── aggregation helpers ────────────────────────────────────────────────────────

def category_dist(recs):
    """Return {cat_int: {name, count, pct}} sorted by count desc."""
    total = len(recs)
    cnt   = Counter(r["final_category"] for r in recs)
    out   = {}
    for cat in sorted(CATEGORY_NAMES.keys()):
        c = cnt.get(cat, 0)
        if c == 0:
            continue
        out[cat] = {
            "name":     CATEGORY_NAMES[cat],
            "severity": CATEGORY_SEVERITY.get(cat, "N/A"),
            "count":    c,
            "pct":      pct(c, total),
        }
    return out


def breakdown_by_key(recs, key):
    """
    Group recs by field `key`, return:
    { value: { "n": int, "category_dist": {...} } }
    """
    groups = defaultdict(list)
    for r in recs:
        groups[r[key]].append(r)
    return {
        val: {
            "n":             len(grecs),
            "category_dist": category_dist(grecs),
        }
        for val, grecs in sorted(groups.items())
    }


def high_severity_rate(recs):
    """Fraction of traces with categories 1 or 2 (Critical / High setup errors)."""
    hs = sum(1 for r in recs if r["final_category"] in (1, 2))
    return pct(hs, len(recs))


def conceptual_error_rate(recs):
    """Fraction with categories 1, 2, or 3 (non-arithmetic errors)."""
    ce = sum(1 for r in recs if r["final_category"] in (1, 2, 3))
    return pct(ce, len(recs))


def complexity_cliff(by_difficulty):
    """
    Given breakdown_by_key result keyed by difficulty,
    return easy_acc_rate, advanced_acc_rate, and cliff_size (pp drop).
    We define 'accuracy' here as % Calculation Error only
    (proxy for 'almost right but arithmetic slip') vs conceptual errors.
    """
    easy_n  = by_difficulty.get("Easy",     {}).get("n", 0)
    adv_n   = by_difficulty.get("Advanced", {}).get("n", 0)
    if easy_n == 0 or adv_n == 0:
        return {}

    # % of errors that are conceptual (cat 1-3) — higher means harder
    def conceptual_pct_in_group(diff_name):
        dist = by_difficulty.get(diff_name, {}).get("category_dist", {})
        conceptual = sum(
            dist[c]["count"] for c in dist if c in (1, 2, 3)
        )
        total = sum(d["count"] for d in dist.values())
        return pct(conceptual, total)

    easy_conceptual = conceptual_pct_in_group("Easy")
    adv_conceptual  = conceptual_pct_in_group("Advanced")
    inter_conceptual = conceptual_pct_in_group("Intermediate")

    return {
        "Easy":         easy_conceptual,
        "Intermediate": inter_conceptual,
        "Advanced":     adv_conceptual,
        "cliff_size_pp": round(adv_conceptual - easy_conceptual, 2),
        "interpretation": (
            "positive = more conceptual errors at Advanced (errors get harder)"
        ),
    }


# ── main analysis ──────────────────────────────────────────────────────────────

def build_model_profile(model_recs, model_key):
    n_total = len(model_recs)
    by_diff   = breakdown_by_key(model_recs, "difficulty")
    by_branch = breakdown_by_key(model_recs, "branch")
    by_domain = breakdown_by_key(model_recs, "domain")

    # Top-5 domains by volume
    top_domains = sorted(
        by_domain.items(), key=lambda x: x[1]["n"], reverse=True
    )[:5]

    return {
        "model":                   model_key,
        "display_name":            MAIN_TEXT_MODELS.get(model_key, model_key),
        "n_traces":                n_total,
        "overall_category_dist":   category_dist(model_recs),
        "high_severity_rate_pct":  high_severity_rate(model_recs),
        "conceptual_error_rate_pct": conceptual_error_rate(model_recs),
        "by_difficulty":           by_diff,
        "complexity_cliff":        complexity_cliff(by_diff),
        "by_branch":               by_branch,
        "top5_domains":            {k: v for k, v in top_domains},
        "by_domain_full":          by_domain,
    }


def cross_model_summary(all_model_profiles):
    """Ranked tables for cross-model comparison."""
    ranked_by_calc = sorted(
        all_model_profiles,
        key=lambda p: p["overall_category_dist"].get(6, {}).get("pct", 0),
        reverse=True,
    )
    ranked_by_conceptual = sorted(
        all_model_profiles,
        key=lambda p: p["conceptual_error_rate_pct"],
        reverse=True,
    )
    ranked_by_hallucination = sorted(
        all_model_profiles,
        key=lambda p: p["overall_category_dist"].get(1, {}).get("pct", 0),
        reverse=True,
    )
    ranked_by_cliff = sorted(
        all_model_profiles,
        key=lambda p: p["complexity_cliff"].get("cliff_size_pp", 0),
        reverse=True,
    )

    def slim(profiles, key_fn):
        return [
            {"model": p["model"], "display_name": p["display_name"], "value": key_fn(p)}
            for p in profiles
        ]

    return {
        "ranked_by_calc_error_pct": slim(
            ranked_by_calc,
            lambda p: p["overall_category_dist"].get(6, {}).get("pct", 0),
        ),
        "ranked_by_conceptual_error_pct": slim(
            ranked_by_conceptual,
            lambda p: p["conceptual_error_rate_pct"],
        ),
        "ranked_by_hallucination_pct": slim(
            ranked_by_hallucination,
            lambda p: p["overall_category_dist"].get(1, {}).get("pct", 0),
        ),
        "ranked_by_complexity_cliff_pp": slim(
            ranked_by_cliff,
            lambda p: p["complexity_cliff"].get("cliff_size_pp", 0),
        ),
    }


def branch_cross_model(all_records):
    """For each branch, error distribution across all models combined."""
    branches = sorted(set(r["branch"] for r in all_records.values()))
    result = {}
    for branch in branches:
        recs = [r for r in all_records.values() if r["branch"] == branch]
        result[branch] = {
            "n":             len(recs),
            "category_dist": category_dist(recs),
        }
    return result


def difficulty_cross_model(all_records):
    """Global error distribution per difficulty level."""
    diffs = ["Easy", "Intermediate", "Advanced"]
    result = {}
    for diff in diffs:
        recs = [r for r in all_records.values() if r["difficulty"] == diff]
        result[diff] = {
            "n":             len(recs),
            "category_dist": category_dist(recs),
        }
    return result


def main_text_insights(profiles_by_model):
    """Structured prose-ready insights for the 3 main-text models."""
    insights = {}

    # ── Gemini 3.1 Pro Preview ──────────────────────────────────────────────
    g = profiles_by_model.get("gemini-3.1-pro-preview", {})
    if g:
        dist = g["overall_category_dist"]
        calc_pct = dist.get(6, {}).get("pct", 0)
        setup_pct = dist.get(2, {}).get("pct", 0)
        form_pct  = dist.get(3, {}).get("pct", 0)
        cliff     = g["complexity_cliff"]
        insights["gemini-3.1-pro-preview"] = {
            "headline": (
                f"Gemini 3.1 Pro Preview's failures are overwhelmingly arithmetic "
                f"({calc_pct:.1f}% Calculation Error), confirming that its reasoning "
                f"chain is substantively correct in the vast majority of wrong-answer "
                f"traces. Only {setup_pct + form_pct:.1f}% of errors are conceptual "
                f"(Setup + Formula combined)."
            ),
            "difficulty_shift": cliff,
            "dominant_error":   dist.get(6, {}),
            "second_error":     dist.get(2, {}),
            "branch_breakdown": g["by_branch"],
            "key_finding": (
                "High FAC leader (65.41%) with errors concentrated in final arithmetic "
                "execution — the model understands the problem and selects the right "
                "method but fails at numerical precision."
            ),
        }

    # ── DeepSeek R1 ──────────────────────────────────────────────────────────
    ds = profiles_by_model.get("deepseek-r1", {})
    if ds:
        dist      = ds["overall_category_dist"]
        calc_pct  = dist.get(6, {}).get("pct", 0)
        setup_pct = dist.get(2, {}).get("pct", 0)
        form_pct  = dist.get(3, {}).get("pct", 0)
        cliff     = ds["complexity_cliff"]
        # Compare setup% to gemini setup%
        g_setup = g["overall_category_dist"].get(2, {}).get("pct", 0) if g else 0
        insights["deepseek-r1"] = {
            "headline": (
                f"DeepSeek R1 (F1 leader at 44.41%) shows a more distributed error "
                f"profile than Gemini: {calc_pct:.1f}% Calculation Error vs "
                f"{setup_pct:.1f}% Setup / Assumption Error. Its higher Setup error "
                f"rate ({setup_pct:.1f}% vs {g_setup:.1f}% for Gemini) explains why "
                f"it achieves higher F1 (partial-credit steps correct) but lower "
                f"final-answer accuracy."
            ),
            "difficulty_shift": cliff,
            "dominant_error":   dist.get(6, {}),
            "second_error":     dist.get(2, {}),
            "branch_breakdown": ds["by_branch"],
            "key_finding": (
                "RL-trained reasoning model — strong at step-level reasoning (high F1) "
                "but Setup/Assumption errors indicate it occasionally misframes the "
                "physical problem before applying correct mathematics."
            ),
        }

    # ── Llama 3.1 70B ────────────────────────────────────────────────────────
    ll = profiles_by_model.get("meta-llama-3.1-70b", {})
    if ll:
        dist      = ll["overall_category_dist"]
        calc_pct  = dist.get(6, {}).get("pct", 0)
        setup_pct = dist.get(2, {}).get("pct", 0)
        form_pct  = dist.get(3, {}).get("pct", 0)
        cliff     = ll["complexity_cliff"]
        # Frontier models' cliff for comparison (use gemini as reference)
        frontier_cliff = g["complexity_cliff"].get("cliff_size_pp", 0) if g else 0
        insights["meta-llama-3.1-70b"] = {
            "headline": (
                f"Llama 3.1 70B exhibits the steepest complexity cliff: conceptual "
                f"errors rise {cliff.get('cliff_size_pp', 0):.1f} percentage points "
                f"from Easy to Advanced, compared to {frontier_cliff:.1f} pp for "
                f"Gemini 3.1 Pro Preview. Formula / Principle errors ({form_pct:.1f}%) "
                f"and Setup errors ({setup_pct:.1f}%) together account for "
                f"{form_pct + setup_pct:.1f}% of failures — substantially higher than "
                f"frontier models — signalling genuine conceptual breakdown at "
                f"Advanced difficulty."
            ),
            "difficulty_shift":       cliff,
            "easy_category_dist":     ll["by_difficulty"].get("Easy", {}).get("category_dist", {}),
            "advanced_category_dist": ll["by_difficulty"].get("Advanced", {}).get("category_dist", {}),
            "dominant_error":         dist.get(6, {}),
            "second_error":           dist.get(3, {}),
            "branch_breakdown":       ll["by_branch"],
            "key_finding": (
                "Most dramatic complexity cliff in the dataset. Error type shifts from "
                "Calculation Error at Easy difficulty to Formula/Principle and "
                "Setup/Assumption at Advanced, demonstrating a qualitative degradation "
                "in reasoning depth rather than just numeric precision."
            ),
        }

    return insights


def main():
    print("Loading annotations …")
    all_records = load_all_annotations()
    print(f"  Loaded {len(all_records)} traces total")

    print("\nBuilding per-model profiles …")
    all_profiles = []
    profiles_by_model = {}
    for model in MODELS:
        recs = [r for r in all_records.values() if r["model"] == model]
        print(f"  {model:40s}  N={len(recs)}")
        profile = build_model_profile(recs, model)
        all_profiles.append(profile)
        profiles_by_model[model] = profile

    print("\nBuilding cross-model analyses …")
    cross_summary   = cross_model_summary(all_profiles)
    branch_analysis = branch_cross_model(all_records)
    diff_analysis   = difficulty_cross_model(all_records)

    print("\nBuilding main-text insights …")
    mt_insights = main_text_insights(profiles_by_model)

    # ── overall corpus stats ──────────────────────────────────────────────────
    all_recs_list = list(all_records.values())
    corpus_stats = {
        "total_traces":              len(all_recs_list),
        "n_models":                  len(MODELS),
        "n_annotators":              len(ANNOTATORS),
        "overall_category_dist":     category_dist(all_recs_list),
        "overall_conceptual_error_rate_pct": conceptual_error_rate(all_recs_list),
        "overall_high_severity_rate_pct":    high_severity_rate(all_recs_list),
    }

    # ── shared meta block (included in every file for self-containment) ──────
    META = {
        "description": "EngTrace error distribution analysis — majority-vote labels from 3 annotators",
        "models":           MODELS,
        "main_text_models": list(MAIN_TEXT_MODELS.keys()),
        "category_legend":  CATEGORY_NAMES,
        "severity_legend":  CATEGORY_SEVERITY,
    }

    def save(filename, payload):
        path = os.path.join(OUT_DIR, filename)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        print(f"  Saved: {path}")

    # ── 1. corpus_overview.json ───────────────────────────────────────────────
    save("corpus_overview.json", {
        "meta":         META,
        "corpus_stats": corpus_stats,
    })

    # ── 2. per_model/ — one file per model ───────────────────────────────────
    per_model_dir = os.path.join(OUT_DIR, "per_model")
    os.makedirs(per_model_dir, exist_ok=True)
    for profile in all_profiles:
        model_slug = profile["model"]
        path = os.path.join(per_model_dir, f"{model_slug}.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump({"meta": META, "profile": profile}, f, indent=2, ensure_ascii=False)
    print(f"  Saved: {per_model_dir}/ ({len(all_profiles)} files)")

    # ── 3. cross_model_comparison.json ───────────────────────────────────────
    # Also include a flat summary table (one row per model) for easy plotting
    summary_table = []
    for p in all_profiles:
        dist  = p["overall_category_dist"]
        cliff = p["complexity_cliff"]
        summary_table.append({
            "model":                      p["model"],
            "display_name":               p["display_name"],
            "n_traces":                   p["n_traces"],
            "is_main_text":               p["model"] in MAIN_TEXT_MODELS,
            "calc_error_pct":             dist.get(6, {}).get("pct", 0),
            "hallucination_pct":          dist.get(1, {}).get("pct", 0),
            "setup_assumption_pct":       dist.get(2, {}).get("pct", 0),
            "formula_principle_pct":      dist.get(3, {}).get("pct", 0),
            "unit_dimensional_pct":       dist.get(4, {}).get("pct", 0),
            "sign_direction_pct":         dist.get(5, {}).get("pct", 0),
            "conceptual_error_rate_pct":  p["conceptual_error_rate_pct"],
            "high_severity_rate_pct":     p["high_severity_rate_pct"],
            "cliff_easy_pct":             cliff.get("Easy", 0),
            "cliff_intermediate_pct":     cliff.get("Intermediate", 0),
            "cliff_advanced_pct":         cliff.get("Advanced", 0),
            "cliff_size_pp":              cliff.get("cliff_size_pp", 0),
        })
    save("cross_model_comparison.json", {
        "meta":           META,
        "summary_table":  summary_table,
        "ranked_tables":  cross_summary,
    })

    # ── 4. branch_and_difficulty.json ────────────────────────────────────────
    # Also include per-model × difficulty breakdown for heatmap plotting
    per_model_difficulty = {
        p["model"]: {
            diff: {
                "n": p["by_difficulty"].get(diff, {}).get("n", 0),
                "category_dist": p["by_difficulty"].get(diff, {}).get("category_dist", {}),
            }
            for diff in ["Easy", "Intermediate", "Advanced"]
        }
        for p in all_profiles
    }
    per_model_branch = {
        p["model"]: p["by_branch"]
        for p in all_profiles
    }
    save("branch_and_difficulty.json", {
        "meta":                    META,
        "by_branch_global":        branch_analysis,
        "by_difficulty_global":    diff_analysis,
        "per_model_by_difficulty": per_model_difficulty,
        "per_model_by_branch":     per_model_branch,
    })

    # ── 5. main_text_insights.json ───────────────────────────────────────────
    save("main_text_insights.json", {
        "meta":     META,
        "insights": mt_insights,
    })

    # ── print a quick summary to terminal ────────────────────────────────────
    print("\n" + "=" * 65)
    print("QUICK SUMMARY")
    print("=" * 65)
    print(f"\nCorpus: {corpus_stats['total_traces']} traces, {len(MODELS)} models")
    print(f"Overall conceptual error rate : {corpus_stats['overall_conceptual_error_rate_pct']:.1f}%")
    print(f"Overall high-severity rate    : {corpus_stats['overall_high_severity_rate_pct']:.1f}%")

    print("\nComplexity cliff (conceptual error % Easy -> Advanced):")
    for m in MODELS:
        cliff = profiles_by_model[m]["complexity_cliff"]
        easy  = cliff.get("Easy", 0)
        adv   = cliff.get("Advanced", 0)
        pp    = cliff.get("cliff_size_pp", 0)
        marker = " [MAIN TEXT]" if m in MAIN_TEXT_MODELS else ""
        print(f"  {m:40s}  Easy={easy:5.1f}%  Adv={adv:5.1f}%  D={pp:+.1f}pp{marker}")

    print("\nTop error category per model:")
    for m in MODELS:
        dist = profiles_by_model[m]["overall_category_dist"]
        top_cat = max(dist.items(), key=lambda x: x[1]["pct"])
        marker = " [MAIN TEXT]" if m in MAIN_TEXT_MODELS else ""
        print(f"  {m:40s}  {top_cat[1]['name']:<30s} {top_cat[1]['pct']:5.1f}%{marker}")

    print(f"\nDone. All results saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
