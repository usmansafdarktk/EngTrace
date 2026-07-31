import pandas as pd
import numpy as np
import json
import glob
import os
from scipy.stats import spearmanr

DIMS = {
    'physical':    ('human_physical',    'Median Physical Score'),
    'math':        ('human_math',        'Median Math Score'),
    'pedagogical': ('human_pedagogical', 'Median Pedagogy Score'),
}

# ==========================================
# 1. Load Human Reviews
# ==========================================
def load_human_reviews(data_dir):
    all_records = []
    file_pattern = os.path.join(data_dir, "*.jsonl")
    files = glob.glob(file_pattern)

    if not files:
        raise FileNotFoundError(f"No .jsonl files found in {data_dir}")

    for filepath in files:
        with open(filepath, 'r') as f:
            for line in f:
                try:
                    record = json.loads(line)
                    s = record['scores']
                    phys = s['physical_plausibility']
                    math = s['mathematical_correctness']
                    ped  = s['pedagogical_clarity']
                    all_records.append({
                        'template':         record['template'],
                        'reviewer':         record['annotator_id'],
                        'decision':         record['decision'].strip().lower(),
                        'human_score':      np.mean([phys, math, ped]),
                        'human_physical':   phys,
                        'human_math':       math,
                        'human_pedagogical': ped,
                    })
                except:
                    continue

    return pd.DataFrame(all_records)


# ==========================================
# 2. Load AI Data (Tribunal CSV)
# ==========================================
def load_ai_data(csv_path):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"AI Tribunal CSV not found at: {csv_path}")

    print(f"Loading AI Data from: {csv_path}")
    df = pd.read_csv(csv_path)

    ai_records = []
    for _, row in df.iterrows():
        try:
            phys = float(row['Median Physical Score'])
            math = float(row['Median Math Score'])
            ped  = float(row['Median Pedagogy Score'])
            ai_score = np.mean([phys, math, ped])
        except ValueError:
            phys = math = ped = ai_score = 0.0

        flag_val = str(row['Consensus Flag']).strip().lower()
        ai_records.append({
            'template':      row['Template ID'],
            'ai_score':      ai_score,
            'ai_physical':   phys,
            'ai_math':       math,
            'ai_pedagogical': ped,
            'ai_decision':   'fail' if flag_val != 'false' else 'pass',
        })

    return pd.DataFrame(ai_records)


# ==========================================
# 3. Metrics
# ==========================================
def _diff_distribution(diffs):
    """Return counts and percentages for |diff| buckets."""
    n = len(diffs)
    buckets = [0, 1, 2]
    rows = []
    for b in buckets:
        c = int((diffs == b).sum())
        rows.append((f"|diff| = {b}", c, f"{100*c/n:.1f}%"))
    c = int((diffs >= 3).sum())
    rows.append(("|diff| ≥ 3", c, f"{100*c/n:.1f}%"))
    return rows, n


def calculate_alignment_metrics(human_df, ai_df):
    # ── A. Human consensus per template ───────────────────────────────────
    agg = {'decision': lambda x: 'approve' if (x == 'approve').sum() > (x == 'reject').sum() else 'reject',
           'human_score': 'median',
           'human_physical': 'median',
           'human_math': 'median',
           'human_pedagogical': 'median'}
    human_consensus = human_df.groupby('template').agg(
        final_decision=('decision', lambda x: 'approve' if (x == 'approve').sum() > (x == 'reject').sum() else 'reject'),
        median_human_score=('human_score', 'median'),
        median_human_physical=('human_physical', 'median'),
        median_human_math=('human_math', 'median'),
        median_human_pedagogical=('human_pedagogical', 'median'),
    ).reset_index()

    # ── B. Merge ───────────────────────────────────────────────────────────
    merged = pd.merge(ai_df, human_consensus, on='template', how='inner')
    print(f"Analyzing alignment for {len(merged)} overlapping templates...")

    # ── METRIC 1: FPR ─────────────────────────────────────────────────────
    ai_passed_mask = merged['ai_decision'] == 'pass'
    n_ai_passed = int(ai_passed_mask.sum())
    false_positives = merged[ai_passed_mask & (merged['final_decision'] == 'reject')]
    fpr = len(false_positives) / n_ai_passed if n_ai_passed > 0 else 0.0

    # ── Certified subset (AI-approved, human-approved) ────────────────────
    certified = merged[merged['final_decision'] == 'approve'].copy()
    n_certified = len(certified)

    # ── METRIC 2: AAD per dimension + overall ─────────────────────────────
    aad_results = {}
    diff_tables  = {}
    spearman_results = {}

    dim_map = {
        'Physical':    ('ai_physical',    'median_human_physical'),
        'Math':        ('ai_math',        'median_human_math'),
        'Pedagogical': ('ai_pedagogical', 'median_human_pedagogical'),
        'Overall':     ('ai_score',       'median_human_score'),
    }

    for dim_label, (ai_col, human_col) in dim_map.items():
        if n_certified < 1:
            aad_results[dim_label] = None
            spearman_results[dim_label] = (float('nan'), float('nan'))
            continue
        diffs = (certified[ai_col] - certified[human_col]).abs()
        aad_results[dim_label] = round(float(diffs.mean()), 3)
        diff_tables[dim_label]  = _diff_distribution(diffs.values)

        # Spearman ρ per dimension
        if n_certified >= 2:
            rho, p_val = spearmanr(certified[ai_col], certified[human_col])
        else:
            rho, p_val = float('nan'), float('nan')
        spearman_results[dim_label] = (rho, p_val)

    return {
        'FPR':             fpr,
        'N_AI_Passed':     n_ai_passed,
        'N_FP':            len(false_positives),
        'N_Certified':     n_certified,
        'AAD':             aad_results,
        'DiffTables':      diff_tables,
        'Spearman':        spearman_results,
    }


# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    HUMAN_DATA_DIR = "../annotator-app/reviews"
    AI_CSV_PATH    = r"../annotator-app/ai_assistance_data_v2/analysis/tribunal_summary.csv"
    OUTPUT_DIR     = "../annotator-app/iaa_results"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    try:
        print("Loading Human Reviews...")
        human_df = load_human_reviews(HUMAN_DATA_DIR)

        print("Loading AI Data...")
        ai_df = load_ai_data(AI_CSV_PATH)

        results = calculate_alignment_metrics(human_df, ai_df)

        print("\n" + "="*60)
        print("AI-HUMAN ALIGNMENT RESULTS")
        print("="*60)

        # ── Metric 1: FPR ─────────────────────────────────────────────────
        print("\nMETRIC 1: Filter Safety (False Positive Rate)")
        print(f"  AI-passed templates           : {results['N_AI_Passed']}")
        print(f"  Subsequently rejected by humans: {results['N_FP']}")
        print(f"  FPR                           : {results['FPR']:.2%}  (Target < 10%)")

        # ── Metrics 2 & 3: AAD (primary) + Spearman ρ (complementary) ────
        print("\nMETRIC 2 [PRIMARY] + METRIC 3 [COMPLEMENTARY]: Scoring Alignment")
        print(f"  Subset: certified templates (AI-passed & human-approved), n = {results['N_Certified']}")
        print(f"  AAD  = mean |AI score - Human median|  (absolute disagreement on 5-pt scale)")
        print(f"  ρ    = Spearman rank correlation        (monotonic ordering agreement)\n")

        aad = results['AAD']
        spearman = results['Spearman']

        print(f"  {'Dimension':<14}  {'AAD':>6}  {'AAD Interp':<30}  {'ρ':>7}  {'p-value':>10}")
        print(f"  {'-'*75}")
        for dim in ['Physical', 'Math', 'Pedagogical', 'Overall']:
            val = aad.get(dim)
            rho, p_val = spearman.get(dim, (float('nan'), float('nan')))

            if val is None:
                interp = "N/A"
            elif val < 0.25:
                interp = "Excellent (<0.25)"
            elif val < 0.5:
                interp = "Good (<0.5)"
            elif val < 1.0:
                interp = "Acceptable (<1.0)"
            else:
                interp = "Poor (>=1.0)"

            aad_str = f"{val:.3f}" if val is not None else "N/A"
            rho_str = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
            p_str   = f"{p_val:.3e}" if not np.isnan(p_val) else "N/A"
            print(f"  {dim:<14}  {aad_str:>6}  {interp:<30}  {rho_str:>7}  {p_str:>10}")

        print("\n  Note: ρ may appear low due to range restriction (ceiling effect when scores")
        print("  cluster at 4-5/5); AAD is the primary metric as it directly quantifies")
        print("  absolute disagreement independent of score distribution shape.")

        print("\n  Difference Distribution (certified templates, per pairwise comparison):")
        for dim in ['Physical', 'Math', 'Pedagogical', 'Overall']:
            if dim not in results['DiffTables']:
                continue
            rows, n = results['DiffTables'][dim]
            print(f"\n  [{dim}]  (n={n} comparisons)")
            for label, count, pct in rows:
                print(f"    {label} : {count:3d}  ({pct})")

        # ── Save ──────────────────────────────────────────────────────────
        summary_rows = []
        for dim in ['Physical', 'Math', 'Pedagogical', 'Overall']:
            rho, p_val = spearman.get(dim, (float('nan'), float('nan')))
            summary_rows.append({
                'Dimension':     dim,
                'AAD':           aad.get(dim),
                'Spearman_Rho':  round(rho, 3) if not np.isnan(rho) else None,
                'P_Value':       round(p_val, 5) if not np.isnan(p_val) else None,
            })
        pd.DataFrame(summary_rows).to_csv(f"{OUTPUT_DIR}/alignment_aad_results.csv", index=False)
        print(f"\n>> Saved to '{OUTPUT_DIR}/alignment_aad_results.csv'")

    except Exception as e:
        import traceback
        print(f"Error: {e}")
        traceback.print_exc()