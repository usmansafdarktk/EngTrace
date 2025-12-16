import pandas as pd
import numpy as np
import json
import glob
import os
from scipy.stats import spearmanr

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
                    
                    # Calculate Reviewer's Average Score
                    scores = [
                        record['scores']['physical_plausibility'],
                        record['scores']['mathematical_correctness'],
                        record['scores']['pedagogical_clarity']
                    ]
                    avg_score = np.mean(scores)
                    
                    all_records.append({
                        'template': record['template'],
                        'reviewer': record['annotator_id'],
                        'decision': record['decision'].strip().lower(), # approve/reject
                        'human_score': avg_score
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
        # 1. Calculate AI Aggregate Score (Average of the 3 medians)
        try:
            phys = float(row['Median Physical Score'])
            math = float(row['Median Math Score'])
            ped = float(row['Median Pedagogy Score'])
            ai_score = np.mean([phys, math, ped])
        except ValueError:
            ai_score = 0.0 # Handle missing/malformed data
            
        # 2. Determine AI Decision based on Consensus Flag
        # If Flag is False -> Pass. If Flag is True -> Fail/Flagged.
        # Ensure we handle string 'False'/'True' or boolean types
        flag_val = str(row['Consensus Flag']).strip().lower()
        if flag_val == 'false':
            ai_decision = 'pass'
        else:
            ai_decision = 'fail'
            
        ai_records.append({
            'template': row['Template ID'],
            'ai_score': ai_score,
            'ai_decision': ai_decision
        })
        
    return pd.DataFrame(ai_records)

# ==========================================
# 3. Calculation Logic (Metrics)
# ==========================================
def calculate_alignment_metrics(human_df, ai_df):
    # A. Aggregate Human Consensus 
    human_consensus = human_df.groupby('template').agg(
        # Check if 'approve' is the majority decision
        final_decision=('decision', lambda x: 'approve' if (x == 'approve').sum() > (x == 'reject').sum() else 'reject'),
        # Median Human Score (Robust to outliers)
        median_human_score=('human_score', 'median')
    ).reset_index()
    
    # B. Merge AI and Human Data 
    merged = pd.merge(ai_df, human_consensus, on='template', how='inner')
    
    print(f"Analyzing alignment for {len(merged)} overlapping templates...")
    
    # METRIC 1: Filter Safety (False Positive Rate) 
    # FPR = Proportion of templates passed by AI but rejected by humans
    
    # 1. Identify templates the AI passed
    ai_passed_mask = merged['ai_decision'] == 'pass'
    n_total_ai_passed = ai_passed_mask.sum()
    
    if n_total_ai_passed == 0:
        fpr = 0.0
        print("Note: AI passed 0 templates (No positives to check).")
    else:
        # 2. Identify which of those were rejected by humans
        false_positives = merged[ai_passed_mask & (merged['final_decision'] == 'reject')]
        n_rejected = len(false_positives)
        fpr = n_rejected / n_total_ai_passed

    # METRIC 2: Scoring Alignment (Spearman's Rho) 
    # Correlation between AI Score and Human Score for CERTIFIED (Approved) templates
    
    certified_subset = merged[merged['final_decision'] == 'approve']
    
    if len(certified_subset) < 2:
        rho = 0.0
        p_val = 1.0
        print("Not enough certified templates to calculate correlation.")
    else:
        # Calculate Spearman Rank Correlation
        rho, p_val = spearmanr(certified_subset['ai_score'], certified_subset['median_human_score'])

    return {
        "FPR": fpr,
        "N_AI_Passed": int(n_total_ai_passed),
        "N_Human_Rejected": int(len(false_positives)) if n_total_ai_passed > 0 else 0,
        "Spearman_Rho": rho,
        "P_Value": p_val,
        "N_Certified": len(certified_subset)
    }

# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    # CONFIGURATION 
    HUMAN_DATA_DIR = "reviews"
    # Use raw string (r"path") to avoid Windows path escape issues
    AI_CSV_PATH = r"ai_assistance_data_v2\analysis\tribunal_summary.csv" 
    
    try:
        # 1. Load Human Reviews
        print("Loading Human Reviews...")
        human_df = load_human_reviews(HUMAN_DATA_DIR)
        
        # 2. Load Real AI Data
        print("Loading AI Data...")
        ai_df = load_ai_data(AI_CSV_PATH)
        
        # 3. Calculate Metrics
        results = calculate_alignment_metrics(human_df, ai_df)
        
        # 4. Report
        print("\n" + "="*60)
        print("IV. Validation of AI-Human Alignment Results")
        print("="*60)
        
        print(f"\nMETRIC 1: Filter Safety (False Positive Rate)")
        print(f"")
        print(f"Total Templates Passed by AI: {results['N_AI_Passed']}")
        print(f"Subsequently Rejected by Humans: {results['N_Human_Rejected']}")
        print(f"FPR: {results['FPR']:.2%} (Target < 10%)")
        
        print(f"\nMETRIC 2: Scoring Alignment (Spearman's Rho)")
        print(f"")
        print(f"Subset (Certified Templates): n = {results['N_Certified']}")
        print(f"Spearman's Rho: {results['Spearman_Rho']:.3f}")
        print(f"P-Value: {results['P_Value']:.3e}")
        
        if results['Spearman_Rho'] > 0.5:
            print(">> Strong Alignment")
        elif results['Spearman_Rho'] > 0.3:
            print(">> Moderate Alignment")
        else:
            print(">> Weak Alignment (or lack of variance in scores)")
            
    except Exception as e:
        print(f"Error: {e}")
