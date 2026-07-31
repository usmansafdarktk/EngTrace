import pandas as pd
import os

def analyze_domain_leaders(csv_path):
    if not os.path.exists(csv_path):
        print(f"Error: File not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)

    # Key metrics to analyze
    metrics = {
        "final_answer_match": "Final Answer Accuracy",
        "step_f1": "Reasoning F1",
        "bertscore": "BERTScore" 
    }

    # Get structure: Branch -> Domain
    if 'branch' not in df.columns or 'domain' not in df.columns:
        print("Error: CSV must contain 'branch' and 'domain' columns.")
        return

    branches = df['branch'].unique()

    print(f"{'='*80}")
    print(f"TOP MODELS PER ENGINEERING DOMAIN")
    print(f"{'='*80}")

    for branch in sorted(branches):
        print(f"\n>>> BRANCH: {branch.upper()} <<<\n")
        
        # Filter for this branch
        branch_df = df[df['branch'] == branch]
        domains = branch_df['domain'].unique()
        
        for domain in sorted(domains):
            print(f"  [Domain: {domain}]")
            domain_df = branch_df[branch_df['domain'] == domain]
            
            # Find top model for each metric
            for col, name in metrics.items():
                if col in domain_df.columns:
                    # Get row with max value
                    best_idx = domain_df[col].idxmax()
                    best_row = domain_df.loc[best_idx]
                    
                    model = best_row['model']
                    score = best_row[col] * 100 # Convert to percentage
                    
                    print(f"    - {name:<22}: {model:<30} ({score:.2f}%)")
            print("-" * 60)

if __name__ == "__main__":
    analyze_domain_leaders("EngChain_Model_Comparison_Domain.csv")
