import os
import json
import glob
import pandas as pd
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================
ROOT_DIR = "Evaluation Results EngChain FINAL"

# Map folder names to clean category names
CATEGORY_MAP = {
    "Frontier Proprietary LLMs": "Frontier",
    "General Purpose Open LLMs": "Open-Weights",
    "Math Enhanced LLMs": "Math-Enhanced"
}

def analyze_evaluation_results(root_dir):
    """
    Traverses the directory, parses JSONL files, and aggregates error stats.
    """
    records = []
    
    # 1. TRAVERSE DIRECTORIES
    # --------------------------------------
    search_pattern = os.path.join(root_dir, "**", "*_evaluation_results_FINAL.jsonl")
    files = glob.glob(search_pattern, recursive=True)
    
    if not files:
        print(f"No JSONL files found in {root_dir}")
        return

    print(f"Found {len(files)} files. Processing...")

    for filepath in files:
        # Extract metadata from file path
        path_parts = os.path.normpath(filepath).split(os.sep)
        
        # Fallback if structure is different, but tries to grab parent folder
        category_folder = path_parts[-2] if len(path_parts) >= 2 else "Unknown"
        category = CATEGORY_MAP.get(category_folder, category_folder)
        
        # Extract model name from filename (remove _evaluation_results_FINAL.jsonl)
        filename = os.path.basename(filepath)
        model_name = filename.replace("_evaluation_results_FINAL.jsonl", "")

        # 2. PARSE JSONL FILE
        # --------------------------------------
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    entry = json.loads(line)
                    meta = entry.get('meta', {})
                    scores = entry.get('scores', {})
                    
                    # Basic Identifiers
                    is_triggered = meta.get('tribunal_triggered', False)
                    error_cat = meta.get('dominant_error_category', "N/A")
                    final_acc = scores.get('final_answer_acc', 0.0)
                    f1 = scores.get('recovered_f1', 0.0)
                    level = entry.get('level', 'Unknown')
                    
                    # 3. DEFINE METRICS
                    # --------------------------------------
                    
                    # "Saved by Tribunal": Tier 1 said NO, Tribunal said YES (Alt Correct)
                    is_saved = (is_triggered and error_cat == "Alternative Correct")
                    
                    # "True Error": Tribunal said NO (Calc or Concept error)
                    is_true_error = (is_triggered and error_cat in ["Calculation Error", "Conceptual Error"])
                    
                    # "Lucky Guess": Right Answer (Acc=1.0) but Low Reasoning (F1 < 0.2)
                    is_lucky = (final_acc == 1.0 and f1 < 0.2)

                    records.append({
                        "Model": model_name,
                        "Category": category,
                        "Level": level,
                        "Total Samples": 1,
                        "Tribunal Triggered": 1 if is_triggered else 0,
                        "Alt Correct (Saved)": 1 if is_saved else 0,
                        "Calculation Error": 1 if error_cat == "Calculation Error" else 0,
                        "Conceptual Error": 1 if error_cat == "Conceptual Error" else 0,
                        "Lucky Guess": 1 if is_lucky else 0
                    })
                    
        except Exception as e:
            print(f"Error reading {filename}: {e}")

    # 4. AGGREGATE DATA
    # --------------------------------------
    df = pd.DataFrame(records)
    
    # Helper to print formatted tables
    def print_summary(grouped_df, title):
        print(f"\n{'='*60}")
        print(f"  {title}")
        print(f"{'='*60}")
        
        # Calculate Rates
        summary = grouped_df.sum().copy()
        
        # Trigger Rate: % of total questions sent to Tribunal
        summary['Trigger Rate %'] = (summary['Tribunal Triggered'] / summary['Total Samples']) * 100
        
        # False Positive Rate (of Triggered): % of Triggered that were actually Alt Correct
        # Handle division by zero if no triggers
        summary['False Positive Rate (of Triggered) %'] = summary.apply(
            lambda x: (x['Alt Correct (Saved)'] / x['Tribunal Triggered'] * 100) if x['Tribunal Triggered'] > 0 else 0, axis=1
        )
        
        # Error Ratio: Calc vs Concept (among true errors)
        # We calculate the ratio of Calc Errors to Total True Errors (Calc + Concept)
        summary['Calc Error Share %'] = summary.apply(
            lambda x: (x['Calculation Error'] / (x['Calculation Error'] + x['Conceptual Error']) * 100) 
            if (x['Calculation Error'] + x['Conceptual Error']) > 0 else 0, axis=1
        )

        # Clean up output columns
        display_cols = [
            'Total Samples', 
            'Tribunal Triggered', 
            'Trigger Rate %',
            'Alt Correct (Saved)', 
            'False Positive Rate (of Triggered) %',
            'Calculation Error',
            'Conceptual Error',
            'Lucky Guess'
        ]
        
        print(summary[display_cols].round(2).to_string())

    # 5. GENERATE REPORTS
    # --------------------------------------
    
    # A. Overall Summary
    print_summary(df.groupby(lambda x: "OVERALL"), "OVERALL AGGREGATED STATS")
    
    # B. Category Summary
    print_summary(df.groupby("Category"), "STATS BY MODEL CATEGORY")
    
    # C. Level Summary (Complexity Cliff Check)
    print_summary(df.groupby("Level"), "STATS BY DIFFICULTY LEVEL")
    
    # D. Detailed Model Breakdown (Optional, simplified)
    print("\n" + "="*60)
    print("  DETAILED MODEL BREAKDOWN (Top 5 by Saved Count)")
    print("="*60)
    model_stats = df.groupby(["Category", "Model"])[['Alt Correct (Saved)', 'Calculation Error', 'Conceptual Error']].sum()
    print(model_stats.sort_values("Alt Correct (Saved)", ascending=False).head(10).to_string())

if __name__ == "__main__":
    analyze_evaluation_results(ROOT_DIR)
