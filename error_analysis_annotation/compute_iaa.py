"""
Inter-Annotator Agreement (IAA) computation for EngTrace error analysis annotations.

Metrics computed
----------------
Pairwise  : Percent Agreement, Cohen's κ, PABAK, Gwet's AC1
Multi-rater: Fleiss' κ, Krippendorff's α (nominal)
  Note: for 3 raters with complete data (no missing labels) on a nominal scale,
  Fleiss' κ and Krippendorff's α are mathematically equivalent — identical values
  are expected and are not a bug.

PABAK (Prevalence-Adjusted Bias-Adjusted Kappa = 2·p₀ − 1) corrects for the
kappa paradox: when one category dominates (high prevalence), Cohen's κ is
artificially deflated even with high raw agreement. PABAK and Gwet's AC1 are
stable under such imbalance and are reported alongside Cohen's κ.

Results saved to agreement_results/
"""

import json
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ANNOT_DIR = BASE_DIR  # files are under Annotator A/B/C sub-dirs relative to this script
OUT_DIR = os.path.join(BASE_DIR, "agreement_results")
os.makedirs(OUT_DIR, exist_ok=True)

ANNOTATORS = ["Annotator A", "Annotator B", "Annotator C"]
MODELS = [
    "claude-opus-4-7", "deepseek-r1", "deepseek-v4-pro",
    "gemini-3.1-pro-preview", "gemma-3-27b", "gpt-5", "gpt-5-mini",
    "mathstral-7b", "meta-llama-3.1-70b", "qwen2.5-72b", "wizardmath-7b",
]
CATEGORIES = [0, 1, 2, 3, 4, 5, 6]  # 0 = No Error, 1-6 = error types
CATEGORY_NAMES = {
    0: "No Error",
    1: "Hallucination",
    2: "Setup/Assumption Error",
    3: "Formula/Principle Error",
    4: "Unit/Dimensional Error",
    5: "Sign/Direction Error",
    6: "Calculation Error",
}
# Map string category names → integer codes.
# Ordered from most-specific to least-specific so longer matches win first.
NAME_TO_CODE_ORDERED = [
    ("no error",                    0),
    ("hallucination",               1),
    ("setup/assumption error",      2),
    ("setup/assumption",            2),
    ("setup",                       2),
    ("formula/principle error",     3),
    ("formula/principle",           3),
    ("formula",                     3),
    ("unit/dimensional error",      4),
    ("unit/dimensional",            4),
    ("unit",                        4),
    ("sign/direction error",        5),
    ("sign/direction",              5),
    ("sign",                        5),
    ("calculation error",           6),
    ("calculation",                 6),
]

# ── data loading ───────────────────────────────────────────────────────────────

def load_annotations():
    """Return dict: annotator -> model -> {annotation_id: category}"""
    data = {a: {} for a in ANNOTATORS}
    for annotator in ANNOTATORS:
        for model in MODELS:
            fpath = os.path.join(ANNOT_DIR, annotator, f"{model}_annotation_sample.jsonl")
            if not os.path.exists(fpath):
                print(f"  WARNING: missing {fpath}")
                continue
            records = {}
            with open(fpath, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    aid = rec["annotation_id"]
                    cat = rec.get("error_category")
                    # Missing / empty category is treated as 0 (No Error):
                    # an absent annotation_category means the annotator found no error.
                    if cat is None or cat == "" or (isinstance(cat, float) and np.isnan(cat)):
                        cat = 0
                    elif isinstance(cat, str):
                        # Resolve string name to integer code using priority-ordered lookup
                        key = cat.lower().strip()
                        matched = next(
                            (code for name, code in NAME_TO_CODE_ORDERED if name in key),
                            None
                        )
                        cat = matched if matched is not None else 0
                    records[aid] = int(cat)
            data[annotator][model] = records
    return data


def build_aligned_matrix(data, model=None):
    """
    Build an (N, 3) matrix of labels for items present in all 3 annotators.
    If model is None, pool all models.
    """
    models = [model] if model else MODELS
    rows = []
    for m in models:
        ann_dicts = [data[a].get(m, {}) for a in ANNOTATORS]
        all_ids = set(ann_dicts[0].keys())
        for d in ann_dicts[1:]:
            all_ids &= set(d.keys())
        for aid in sorted(all_ids):
            rows.append([d[aid] for d in ann_dicts])
    return np.array(rows, dtype=float) if rows else np.empty((0, 3))


# ── metric helpers ─────────────────────────────────────────────────────────────

def percent_agreement(a, b):
    assert len(a) == len(b)
    return np.mean(a == b) * 100


def cohens_kappa(a, b):
    labels = sorted(set(a) | set(b))
    n = len(a)
    if n == 0:
        return float("nan")
    po = np.mean(a == b)
    pe = sum(
        (np.sum(a == c) / n) * (np.sum(b == c) / n)
        for c in labels
    )
    if pe == 1.0:
        return 1.0
    return (po - pe) / (1.0 - pe)


def fleiss_kappa(matrix):
    """Fleiss' κ for an (N, k) matrix of integer labels."""
    N, k = matrix.shape
    if N == 0:
        return float("nan")
    cats = sorted(set(matrix.flatten().tolist()))
    n_cats = len(cats)
    cat_idx = {c: i for i, c in enumerate(cats)}

    counts = np.zeros((N, n_cats))
    for i in range(N):
        for j in range(k):
            counts[i, cat_idx[matrix[i, j]]] += 1

    P_i = (np.sum(counts ** 2, axis=1) - k) / (k * (k - 1))
    P_bar = np.mean(P_i)
    p_j = counts.sum(axis=0) / (N * k)
    P_e = np.sum(p_j ** 2)

    if P_e == 1.0:
        return 1.0
    return (P_bar - P_e) / (1.0 - P_e)


def krippendorffs_alpha(matrix):
    """Krippendorff's α, nominal scale, for an (N, k) matrix."""
    N, k = matrix.shape
    if N == 0:
        return float("nan")

    cats = sorted(set(v for v in matrix.flatten() if not np.isnan(v)))
    n_cats = len(cats)
    cat_idx = {c: i for i, c in enumerate(cats)}

    o = np.zeros((n_cats, n_cats))
    for i in range(N):
        vals = [matrix[i, j] for j in range(k) if not np.isnan(matrix[i, j])]
        m_i = len(vals)
        if m_i < 2:
            continue
        for u in range(len(vals)):
            for v in range(len(vals)):
                if u != v:
                    o[cat_idx[vals[u]], cat_idx[vals[v]]] += 1.0 / (m_i - 1)

    n_c = o.sum(axis=1)
    n_total = o.sum()
    if n_total == 0:
        return float("nan")

    D_o = sum(o[c, k_] for c in range(n_cats) for k_ in range(n_cats) if c != k_) / n_total
    D_e = sum(
        (n_c[c] * n_c[k_]) / (n_total ** 2)
        for c in range(n_cats) for k_ in range(n_cats) if c != k_
    )

    if D_e == 0:
        return 1.0
    return 1.0 - D_o / D_e


def pabak(a, b):
    """PABAK = 2·p₀ − 1.  Prevalence- and bias-adjusted kappa.
    Corrects for the kappa paradox when one category dominates the distribution.
    Range: −1 to 1; interpretation same as Cohen's κ."""
    return 2.0 * np.mean(a == b) - 1.0


def gwets_ac1(a, b):
    """Gwet's AC1 — prevalence-robust alternative to Cohen's κ.
    Uses average marginal probabilities rather than their product to estimate
    chance agreement, making it stable when one category is very frequent.
    Reference: Gwet (2008), Journal of Classification."""
    n = len(a)
    if n == 0:
        return float("nan")
    labels = sorted(set(a) | set(b))
    k = len(labels)
    if k <= 1:
        return 1.0
    po = np.mean(a == b)
    # Average marginal probability for each category
    pi_q = [(np.sum(a == c) / n + np.sum(b == c) / n) / 2.0 for c in labels]
    # Chance-agreement estimator
    pe = sum(p * (1.0 - p) for p in pi_q) / (k - 1)
    if pe == 1.0:
        return 1.0
    return (po - pe) / (1.0 - pe)


def category_prevalence(matrix):
    """Return dict of category → fraction of all labels (pooled across raters)."""
    all_labels = matrix.flatten()
    total = len(all_labels)
    if total == 0:
        return {}
    return {
        cat: float(np.sum(all_labels == cat)) / total
        for cat in CATEGORIES
        if np.sum(all_labels == cat) > 0
    }


def confusion_matrix_pairwise(a, b, labels=None):
    if labels is None:
        labels = sorted(set(a) | set(b))
    idx = {c: i for i, c in enumerate(labels)}
    cm = np.zeros((len(labels), len(labels)), dtype=int)
    for ai, bi in zip(a, b):
        cm[idx[ai], idx[bi]] += 1
    return pd.DataFrame(cm, index=labels, columns=labels)


# ── aggregate metrics ──────────────────────────────────────────────────────────

def compute_all_metrics(matrix):
    if matrix.shape[0] == 0:
        return None

    a, b, c = matrix[:, 0], matrix[:, 1], matrix[:, 2]
    pairs = [("A", "B", a, b), ("A", "C", a, c), ("B", "C", b, c)]

    pct = {f"{p}-{q}": percent_agreement(x, y) for p, q, x, y in pairs}
    pct["overall"] = np.mean(list(pct.values()))

    ck  = {f"{p}-{q}": cohens_kappa(x, y)  for p, q, x, y in pairs}
    pbk = {f"{p}-{q}": pabak(x, y)          for p, q, x, y in pairs}
    pbk["overall"] = np.mean(list(pbk.values()))
    ac1 = {f"{p}-{q}": gwets_ac1(x, y)      for p, q, x, y in pairs}
    ac1["overall"] = np.mean(list(ac1.values()))

    fk = fleiss_kappa(matrix)
    ka = krippendorffs_alpha(matrix)

    prev = category_prevalence(matrix)
    dominant_cat = max(prev, key=prev.get) if prev else None
    dominant_pct = prev[dominant_cat] * 100 if dominant_cat is not None else 0.0

    cms = {}
    for p, q, x, y in pairs:
        cms[f"{p}-{q}"] = confusion_matrix_pairwise(
            x.astype(int), y.astype(int), CATEGORIES
        )

    return {
        "n_items":            matrix.shape[0],
        "percent_agreement":  pct,
        "cohens_kappa":       ck,
        "pabak":              pbk,
        "gwets_ac1":          ac1,
        "fleiss_kappa":       fk,
        "krippendorffs_alpha": ka,
        "prevalence":         prev,
        "dominant_category":  dominant_cat,
        "dominant_pct":       dominant_pct,
        "confusion_matrices": cms,
    }


# ── save helpers ───────────────────────────────────────────────────────────────

def fmt(v):
    return f"{v:.4f}" if not (isinstance(v, float) and np.isnan(v)) else "N/A"


def save_summary_csv(results_by_model, overall):
    rows = []
    for scope, res in [("OVERALL", overall)] + [(m, results_by_model[m]) for m in MODELS]:
        if res is None:
            continue
        pa  = res["percent_agreement"]
        ck  = res["cohens_kappa"]
        pbk = res["pabak"]
        ac1 = res["gwets_ac1"]
        dom = res["dominant_category"]
        rows.append({
            "model":               scope,
            "n_items":             res["n_items"],
            # Raw agreement
            "pct_agree_AB":        round(pa["A-B"],       4),
            "pct_agree_AC":        round(pa["A-C"],       4),
            "pct_agree_BC":        round(pa["B-C"],       4),
            "pct_agree_overall":   round(pa["overall"],   4),
            # Cohen's κ
            "cohens_kappa_AB":     round(ck["A-B"],       4),
            "cohens_kappa_AC":     round(ck["A-C"],       4),
            "cohens_kappa_BC":     round(ck["B-C"],       4),
            # PABAK
            "pabak_AB":            round(pbk["A-B"],      4),
            "pabak_AC":            round(pbk["A-C"],      4),
            "pabak_BC":            round(pbk["B-C"],      4),
            "pabak_overall":       round(pbk["overall"],  4),
            # Gwet's AC1
            "gwets_ac1_AB":        round(ac1["A-B"],      4),
            "gwets_ac1_AC":        round(ac1["A-C"],      4),
            "gwets_ac1_BC":        round(ac1["B-C"],      4),
            "gwets_ac1_overall":   round(ac1["overall"],  4),
            # Multi-rater
            "fleiss_kappa":        round(res["fleiss_kappa"],        4),
            "krippendorffs_alpha": round(res["krippendorffs_alpha"], 4),
            # Prevalence context
            "dominant_category":   CATEGORY_NAMES.get(dom, "?") if dom is not None else "?",
            "dominant_pct":        round(res["dominant_pct"], 2),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "iaa_summary.csv"), index=False)
    print(f"  Saved iaa_summary.csv")


def save_confusion_matrices_excel(results_by_model, overall):
    fpath = os.path.join(OUT_DIR, "confusion_matrices.xlsx")
    with pd.ExcelWriter(fpath) as writer:
        for scope, res in [("OVERALL", overall)] + [(m, results_by_model[m]) for m in MODELS]:
            if res is None:
                continue
            for pair, cm in res["confusion_matrices"].items():
                cm_named = cm.copy()
                cm_named.index = [f"{i}:{CATEGORY_NAMES.get(i,'?')}" for i in cm.index]
                cm_named.columns = [f"{i}:{CATEGORY_NAMES.get(i,'?')}" for i in cm.columns]
                sheet = f"{scope[:18]}_{pair}"
                cm_named.to_excel(writer, sheet_name=sheet)
    print(f"  Saved confusion_matrices.xlsx")


def save_confusion_matrices_text(results_by_model, overall):
    fpath = os.path.join(OUT_DIR, "confusion_matrices_text.txt")
    with open(fpath, "w", encoding="utf-8") as f:
        for scope, res in [("OVERALL", overall)] + [(m, results_by_model[m]) for m in MODELS]:
            if res is None:
                continue
            f.write(f"\n{'='*65}\n  {scope}\n{'='*65}\n")
            for pair, cm in res["confusion_matrices"].items():
                f.write(f"\n  Annotator {pair[0]} (rows) vs Annotator {pair[2]} (cols)\n")
                cm_named = cm.copy()
                short = {k: f"{k}:{v[:10]}" for k, v in CATEGORY_NAMES.items()}
                cm_named.index = [short[i] for i in cm.index]
                cm_named.columns = [str(i) for i in cm.columns]
                f.write(cm_named.to_string())
                f.write("\n")
    print(f"  Saved confusion_matrices_text.txt")


def save_text_report(results_by_model, overall):
    lines = []
    lines.append("=" * 74)
    lines.append("INTER-ANNOTATOR AGREEMENT REPORT — EngTrace Error Analysis")
    lines.append("Annotators: A, B, C   |   Models: 11   |   Scale: 7 categories")
    lines.append("=" * 74)
    lines.append("")
    lines.append("Metrics key")
    lines.append("  % Agree   : Raw pairwise percent agreement")
    lines.append("  Cohen's κ : Chance-corrected agreement (sensitive to prevalence imbalance)")
    lines.append("  PABAK     : Prevalence-Adjusted Bias-Adjusted Kappa = 2·p₀ − 1")
    lines.append("              Corrects for kappa paradox when one category dominates")
    lines.append("  Gwet AC1  : Prevalence-robust alternative to Cohen's κ (Gwet 2008)")
    lines.append("  Fleiss' κ : Multi-rater chance-corrected agreement (all 3 annotators)")
    lines.append("  Kripp. α  : Krippendorff's α, nominal scale (all 3 annotators)")
    lines.append("  NOTE: Fleiss' κ = Krippendorff's α is mathematically expected for 3")
    lines.append("        raters with complete nominal data — not a computational error.")

    def prevalence_note(res):
        dom = res["dominant_category"]
        pct = res["dominant_pct"]
        if dom is None:
            return ""
        name = CATEGORY_NAMES.get(dom, "?")
        note = f"  Dominant category: {dom}:{name} ({pct:.1f}% of all labels)"
        if pct >= 70.0:
            note += "  ← high prevalence; Cohen's κ deflated, see PABAK / Gwet AC1"
        return note

    def block(scope, res):
        lines.append(f"\n{'-'*74}")
        lines.append(f"  Scope : {scope}")
        lines.append(f"  Items : {res['n_items']}")
        pnote = prevalence_note(res)
        if pnote:
            lines.append(pnote)
        lines.append(f"{'-'*74}")
        pa  = res["percent_agreement"]
        ck  = res["cohens_kappa"]
        pbk = res["pabak"]
        ac1 = res["gwets_ac1"]
        lines.append(f"  Percent Agreement")
        lines.append(f"    A vs B  : {pa['A-B']:.2f}%")
        lines.append(f"    A vs C  : {pa['A-C']:.2f}%")
        lines.append(f"    B vs C  : {pa['B-C']:.2f}%")
        lines.append(f"    Overall : {pa['overall']:.2f}%")
        lines.append(f"  Cohen's κ (pairwise)")
        lines.append(f"    A vs B  : {fmt(ck['A-B'])}")
        lines.append(f"    A vs C  : {fmt(ck['A-C'])}")
        lines.append(f"    B vs C  : {fmt(ck['B-C'])}")
        lines.append(f"  PABAK (pairwise, prevalence-adjusted)")
        lines.append(f"    A vs B  : {fmt(pbk['A-B'])}")
        lines.append(f"    A vs C  : {fmt(pbk['A-C'])}")
        lines.append(f"    B vs C  : {fmt(pbk['B-C'])}")
        lines.append(f"    Overall : {fmt(pbk['overall'])}")
        lines.append(f"  Gwet's AC1 (pairwise, prevalence-robust)")
        lines.append(f"    A vs B  : {fmt(ac1['A-B'])}")
        lines.append(f"    A vs C  : {fmt(ac1['A-C'])}")
        lines.append(f"    B vs C  : {fmt(ac1['B-C'])}")
        lines.append(f"    Overall : {fmt(ac1['overall'])}")
        lines.append(f"  Fleiss' κ  (all 3 annotators) : {fmt(res['fleiss_kappa'])}")
        lines.append(f"  Kripp. α   (all 3 annotators) : {fmt(res['krippendorffs_alpha'])}")

        # Category distribution
        prev = res.get("prevalence", {})
        if prev:
            lines.append(f"  Category distribution (pooled across all annotators):")
            for cat, frac in sorted(prev.items(), key=lambda x: -x[1]):
                lines.append(f"    {cat}:{CATEGORY_NAMES.get(cat,'?'):<30s} {frac*100:5.1f}%")

    block("OVERALL (all 11 models pooled)", overall)
    for model in MODELS:
        res = results_by_model.get(model)
        if res:
            block(model, res)

    lines.append("\n\nCategory legend:")
    for k, v in CATEGORY_NAMES.items():
        lines.append(f"  {k}: {v}")

    lines.append("\nKappa / PABAK / AC1 interpretation:")
    lines.append("  < 0.00 : poor       0.01–0.20 : slight    0.21–0.40 : fair")
    lines.append("  0.41–0.60 : moderate   0.61–0.80 : substantial   0.81–1.00 : almost perfect")
    lines.append("\nReferences:")
    lines.append("  Landis & Koch (1977) — kappa benchmarks")
    lines.append("  Byrt et al. (1993)   — PABAK")
    lines.append("  Gwet (2008)          — AC1")

    report = "\n".join(lines)
    fpath = os.path.join(OUT_DIR, "iaa_report.txt")
    with open(fpath, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"  Saved iaa_report.txt")
    return report


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    print("Loading annotations …")
    data = load_annotations()

    print("\nComputing per-model metrics …")
    results_by_model = {}
    for model in MODELS:
        mat = build_aligned_matrix(data, model=model)
        results_by_model[model] = compute_all_metrics(mat)
        n = mat.shape[0] if mat.ndim == 2 else 0
        print(f"  {model:38s}  N={n}")

    print("\nComputing overall metrics …")
    overall_mat = build_aligned_matrix(data, model=None)
    overall = compute_all_metrics(overall_mat)
    print(f"  {'OVERALL':38s}  N={overall_mat.shape[0]}")

    print(f"\nSaving results to: {OUT_DIR}")
    save_summary_csv(results_by_model, overall)
    save_confusion_matrices_excel(results_by_model, overall)
    save_confusion_matrices_text(results_by_model, overall)
    report = save_text_report(results_by_model, overall)

    print("\nDone. Results saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
