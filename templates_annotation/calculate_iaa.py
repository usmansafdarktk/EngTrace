import os
import json
import glob
import numpy as np
import pandas as pd
import krippendorff
from statsmodels.stats.inter_rater import fleiss_kappa, aggregate_raters

pd.set_option('future.no_silent_downcasting', True)


def load_reviews(data_dir):
    all_records = []
    # Pattern matches the file structure shown in the image
    file_pattern = os.path.join(data_dir, "*.jsonl")
    files = glob.glob(file_pattern)
    
    if not files:
        raise FileNotFoundError(f"No .jsonl files found in {data_dir}")

    for filepath in files:
        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f):
                try:
                    record = json.loads(line)
                    
                    # Normalize Strings
                    decision = record['decision'].strip().lower()
                    
                    # Create IDs
                    global_rater_id = f"{record['branch']}_{record['annotator_id']}"
                    local_rater_id = record['annotator_id']
                    
                    flat_record = {
                        'branch': record['branch'],
                        'template': record['template'],
                        'global_rater': global_rater_id,
                        'local_rater': local_rater_id,
                        'decision': decision,
                        # Extract ALL 3 Scores 
                        'score_physical': record['scores']['physical_plausibility'],
                        'score_math': record['scores']['mathematical_correctness'], 
                        'score_pedagogical': record['scores']['pedagogical_clarity']
                    }
                    all_records.append(flat_record)
                except (json.JSONDecodeError, KeyError) as e:
                    print(f"Skipping malformed line {line_num} in {filepath}: {e}")

    return pd.DataFrame(all_records)


def calculate_fleiss_kappa(df):
    results = {}
    
    for branch, group in df.groupby('branch'):
        # Handle Duplicates (take last)
        pivot = group.pivot_table(
            index='template', 
            columns='local_rater', 
            values='decision', 
            aggfunc='last'
        )
        
        # Handle Missing Data (Drop incomplete rows)
        initial_count = len(pivot)
        pivot_clean = pivot.dropna()
        dropped_count = initial_count - len(pivot_clean)
        
        if dropped_count > 0:
            print(f"WARNING: [{branch}] Dropped {dropped_count} templates due to missing reviews.")
            
        if pivot_clean.empty:
            results[branch] = {"Error": "No overlapping reviews found"}
            continue

        # Map decisions to integers (0=reject, 1=approve)
        # Add infer_objects() to suppress FutureWarning
        df_int = pivot_clean.replace({'reject': 0, 'approve': 1}).infer_objects(copy=False)
        
        # Check for perfect agreement BEFORE calculating Kappa
        # If all raters gave the same decision on every template
        if df_int.nunique(axis=1).eq(1).all():
            results[branch] = {
                "Fleiss Kappa": 1.000,
                "Interpretation": "Perfect Agreement",
                "Valid Templates": len(pivot_clean),
                "Raters": pivot_clean.shape[1]
            }
            continue
        
        # Calculate aggregate counts
        agg_data, categories = aggregate_raters(df_int.values)
        
        try:
            kappa = fleiss_kappa(agg_data)
            
            # Handle NaN case as safety net
            if np.isnan(kappa):
                kappa = 1.0
                interp = "Perfect Agreement"
            else:
                # Interpretation (Landis & Koch)
                interp = "Poor"
                if kappa > 0.8: interp = "Almost Perfect"
                elif kappa > 0.6: interp = "Substantial"
                elif kappa > 0.4: interp = "Moderate"
                elif kappa > 0.2: interp = "Fair"
            
            results[branch] = {
                "Fleiss Kappa": round(kappa, 3),
                "Interpretation": interp,
                "Valid Templates": len(pivot_clean),
                "Raters": pivot_clean.shape[1]
            }
        except Exception as e:
            results[branch] = {"Error": str(e)}
        
    return pd.DataFrame(results).T


def calculate_gwet_ac2(df, score_col):
    try:
        from irrCAC.raw import CAC
    except ImportError:
        print("WARNING: irrCAC not installed. Run: pip install irrCAC")
        return {}

    def run_ac2(pivot_templates_x_raters):
        try:
            if pivot_templates_x_raters.shape[1] < 2:
                return None
            # Zero-variance → perfect agreement by convention
            if np.nanvar(pivot_templates_x_raters.values) == 0:
                return 1.0
            cac = CAC(pivot_templates_x_raters, weights='quadratic')
            result = cac.gwet()
            return round(float(result['est']['coefficient_value']), 3)
        except Exception as e:
            print(f"  AC2 error ({score_col}): {e}")
            return None

    results = {}

    # Per-branch: rows=templates, cols=local raters
    for branch, group in df.groupby('branch'):
        pivot = group.pivot_table(
            index='template', columns='local_rater', values=score_col, aggfunc='mean'
        )
        results[branch] = run_ac2(pivot)

    # Global: mean of per-branch values (sparse 90×9 matrix is uncomputable by irrCAC)
    branch_vals = [v for v in results.values() if v is not None]
    results['GLOBAL (Unified)'] = round(float(np.mean(branch_vals)), 3) if branch_vals else None

    return results


def calculate_krippendorff(df, score_col):
    # Helper function to run alpha calculation
    def run_alpha(matrix_vals):
        try:
            # If variance is 0 (e.g., everyone gave '5'), alpha is undefined but effectively 1.0
            if np.nanvar(matrix_vals) == 0:
                return 1.0
            return round(krippendorff.alpha(reliability_data=matrix_vals, level_of_measurement='ordinal'), 3)
        except:
            return 0.0

    results = {}
    
    # 1. Per-Branch Calculation
    for branch, group in df.groupby('branch'):
        # Average duplicate scores if they exist
        pivot = group.pivot_table(index='local_rater', columns='template', values=score_col, aggfunc='mean')
        results[branch] = run_alpha(pivot.values)
        
    # 2. Global Unified Calculation (using global_rater_id)
    # Allows sparse matrix calculation across all branches
    global_pivot = df.pivot_table(index='global_rater', columns='template', values=score_col, aggfunc='mean')
    results['GLOBAL (Unified)'] = run_alpha(global_pivot.values)
    
    return results


def pairwise_difference_distribution(df, score_col, branch):
    from itertools import combinations

    group = df[df['branch'] == branch]
    pivot = group.pivot_table(
        index='template', columns='local_rater', values=score_col, aggfunc='mean'
    ).dropna()

    all_diffs = []
    per_template = []

    for template, row in pivot.iterrows():
        scores = row.values
        diffs = [abs(a - b) for a, b in combinations(scores, 2)]
        all_diffs.extend(diffs)
        per_template.append({
            'template': template,
            'scores': list(scores),
            'max_diff': max(diffs),
            'mean_diff': round(np.mean(diffs), 3),
        })

    all_diffs = np.array(all_diffs)
    n = len(all_diffs)
    buckets = {
        '|diff| = 0': int((all_diffs == 0).sum()),
        '|diff| = 1': int((all_diffs == 1).sum()),
        '|diff| = 2': int((all_diffs == 2).sum()),
        '|diff| >= 3': int((all_diffs >= 3).sum()),
    }

    lines = [f"\nTotal pairwise comparisons : {n}  ({pivot.shape[0]} templates × {len(list(combinations(range(pivot.shape[1]), 2)))} pairs each)"]
    for label, count in buckets.items():
        lines.append(f"  {label} : {count:3d}  ({100*count/n:.1f}%)")
    lines.append(f"\n  Mean |diff| : {all_diffs.mean():.3f}")
    lines.append(f"  Max  |diff| : {int(all_diffs.max())}")

    return {
        'summary': '\n'.join(lines),
        'per_template': pd.DataFrame(per_template).set_index('template'),
        'buckets': buckets,
        'n': n,
    }


# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    DATA_DIR = "../annotator-app/reviews"
    OUTPUT_DIR = "../annotator-app/iaa_results"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    try:
        # Load
        print("Loading data...")
        df = load_reviews(DATA_DIR)

        # 1. Fleiss Kappa (Binary Decision)
        print("\n" + "="*60)
        print("METRIC 1: Binary Decision Agreement (Fleiss' Kappa)")
        print("="*60)
        fleiss_df = calculate_fleiss_kappa(df)
        print(fleiss_df)
        fleiss_df.to_csv(f"{OUTPUT_DIR}/iaa_fleiss_kappa_results.csv")
        print(">> Saved to 'iaa_fleiss_kappa_results.csv'")

        # 2. Krippendorff's Alpha + Gwet's AC2 (Quality Scores)
        print("\n" + "="*60)
        print("METRIC 2: Quality Score Consistency (Krippendorff's α & Gwet's AC2)")
        print("Note: α > 0.67 reliable; AC2 corrects for range restriction.")
        print("="*60)

        alpha_phys = calculate_krippendorff(df, 'score_physical')
        alpha_math = calculate_krippendorff(df, 'score_math')
        alpha_ped  = calculate_krippendorff(df, 'score_pedagogical')

        ac2_phys = calculate_gwet_ac2(df, 'score_physical')
        ac2_math = calculate_gwet_ac2(df, 'score_math')
        ac2_ped  = calculate_gwet_ac2(df, 'score_pedagogical')

        summary = pd.DataFrame({
            'Kripp. α — Physical':      alpha_phys,
            'AC2 — Physical':           ac2_phys,
            'Kripp. α — Math':          alpha_math,
            'AC2 — Math':               ac2_math,
            'Kripp. α — Pedagogical':   alpha_ped,
            'AC2 — Pedagogical':        ac2_ped,
        })
        print(summary.to_string())
        summary.to_csv(f"{OUTPUT_DIR}/iaa_krippendorff_gwet_results.csv")
        print(">> Saved to 'iaa_krippendorff_gwet_results.csv'")

        # Also keep the original Krippendorff-only CSV for backwards compatibility
        kripp_only = pd.DataFrame({
            'Physical Plausibility':    alpha_phys,
            'Mathematical Correctness': alpha_math,
            'Pedagogical Clarity':      alpha_ped,
        })
        kripp_only.to_csv(f"{OUTPUT_DIR}/iaa_krippendorff_results.csv")
        print(">> Saved to 'iaa_krippendorff_results.csv'")

        # 3. Pairwise Difference Distribution (Pedagogical Clarity — Mechanical Eng.)
        print("\n" + "="*60)
        print("DIAGNOSTIC: Pairwise |Difference| Distribution")
        print("Pedagogical Clarity — Mechanical Engineering")
        print("="*60)
        diff_results = pairwise_difference_distribution(df, 'score_pedagogical', 'mechanical_engineering')
        print(diff_results['summary'])
        diff_results['per_template'].to_csv(f"{OUTPUT_DIR}/diag_ped_mech_differences.csv")
        print(">> Per-template detail saved to 'diag_ped_mech_differences.csv'")

    except Exception as e:
        print(f"Pipeline Error: {e}")
