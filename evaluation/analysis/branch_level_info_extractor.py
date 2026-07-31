import pandas as pd
import os

def print_top_models_by_branch(csv_path):
    # Load the branch-level results
    if not os.path.exists(csv_path):
        print(f"Error: File not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)

    # Define the metrics we want to analyze
    metrics = {
        "final_answer_match": "Final Answer Accuracy",
        "step_f1": "Reasoning F1",
        "bertscore": "BERTScore",
        "rouge2": "ROUGE-2",
        "rougeL": "ROUGE-L"
    }

    # Get unique branches
    branches = df['branch'].unique()

    print(f"{'='*60}")
    print(f"TOP MODELS PER ENGINEERING BRANCH")
    print(f"{'='*60}\n")

    for branch in branches:
        print(f" {branch.upper()} ENGINEERING ")
        branch_df = df[df['branch'] == branch]

        for metric_col, metric_name in metrics.items():
            if metric_col not in branch_df.columns:
                continue

            # Find the row with the max value for this metric
            best_row_idx = branch_df[metric_col].idxmax()
            best_row = branch_df.loc[best_row_idx]
            
            model_name = best_row['model']
            score = best_row[metric_col]
            
            # Format score as percentage
            print(f"{metric_name:<25}: {model_name:<30} ({score*100:.2f}%)")
        print("\n")

if __name__ == "__main__":
    print_top_models_by_branch("EngChain_Model_Comparison_Branch.csv")
