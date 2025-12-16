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


# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    DATA_DIR = "reviews" 
    
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
        fleiss_df.to_csv("iaa_results/iaa_fleiss_kappa_results.csv") # Save to CSV
        print(">> Saved to 'iaa_fleiss_kappa_results.csv'")
        
        # 2. Krippendorff's Alpha (Quality Scores)
        print("\n" + "="*60)
        print("METRIC 2: Quality Score Consistency (Krippendorff's Alpha)")
        print("Note: Values > 0.67 are reliable; > 0.80 is good.")
        print("="*60)
        
        # Calculate for ALL 3 Metrics 
        alpha_phys = calculate_krippendorff(df, 'score_physical')
        alpha_math = calculate_krippendorff(df, 'score_math') 
        alpha_ped = calculate_krippendorff(df, 'score_pedagogical')
        
        summary = pd.DataFrame({
            'Physical Plausibility': alpha_phys,
            'Mathematical Correctness': alpha_math, 
            'Pedagogical Clarity': alpha_ped
        })
        print(summary)
        summary.to_csv("iaa_results/iaa_krippendorff_results.csv") # Save to CSV
        print(">> Saved to 'iaa_krippendorff_results.csv'")
        
    except Exception as e:
        print(f"Pipeline Error: {e}")
