import pandas as pd
import os

def analyze_level_leaders(csv_path):
    if not os.path.exists(csv_path):
        print(f"Error: File not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)

    # Clean up level names if needed (sometimes they might be lowercase/uppercase)
    if 'level' in df.columns:
        df['level'] = df['level'].astype(str).str.capitalize()

    # Define the order of difficulty for display
    level_order = ["Easy", "Medium", "Hard"]
    
    # Key metrics to analyze
    metrics = {
        "final_answer_match": "Final Answer Accuracy",
        "step_f1": "Reasoning F1",
        "bertscore": "BERTScore"
    }

    print(f"{'='*60}")
    print(f"TOP MODELS PER DIFFICULTY LEVEL")
    print(f"{'='*60}")

    # Get unique levels present in the data
    available_levels = df['level'].unique()
    
    # Sort levels based on logical order if possible
    sorted_levels = [l for l in level_order if l in available_levels]
    # Add any other levels found that weren't in the predefined list
    for l in available_levels:
        if l not in sorted_levels:
            sorted_levels.append(l)

    for level in sorted_levels:
        print(f"\n>>> LEVEL: {level.upper()} <<<\n")
        
        level_df = df[df['level'] == level]
        
        if level_df.empty:
            print("  No data found for this level.")
            continue

        for metric_col, metric_name in metrics.items():
            if metric_col not in level_df.columns:
                continue

            # Find the row with the max value for this metric
            best_idx = level_df[metric_col].idxmax()
            best_row = level_df.loc[best_idx]
            
            model_name = best_row['model']
            score = best_row[metric_col] * 100 # Convert to percentage
            
            print(f"  - {metric_name:<25}: {model_name:<30} ({score:.2f}%)")
        print("-" * 60)

if __name__ == "__main__":
    analyze_level_leaders("EngChain_Model_Comparison_Level.csv")
