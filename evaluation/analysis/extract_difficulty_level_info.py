import json
import pandas as pd
import os

def analyze_dataset_stats(file_path):
    data = []
    
    print(f"Reading file: {file_path}...")
    
    # Reading the JSONL file line by line
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            try:
                entry = json.loads(line)
                # Extract only necessary fields to save memory
                data.append({
                    "id": entry.get("id"),
                    "branch": entry.get("branch"),
                    "domain": entry.get("domain"),
                    "area": entry.get("area"),
                    "level": entry.get("level")
                })
            except json.JSONDecodeError:
                continue

    # Create DataFrame
    df = pd.read_json(json.dumps(data))
    
    # ---------------------------------------------------------
    # STEP 1: Deduplicate to get Unique Templates
    # ---------------------------------------------------------
    # We drop duplicates based on 'id' because the evaluation file 
    # might contain multiple runs or instances for the same template ID.
    unique_templates = df.drop_duplicates(subset=['id'])
    
    total_templates = len(unique_templates)
    print(f"\nSuccessfully extracted {total_templates} unique templates.")

    # ---------------------------------------------------------
    # STEP 2: Overall Distribution by Difficulty
    # ---------------------------------------------------------
    print("\n" + "="*40)
    print("OVERALL DIFFICULTY DISTRIBUTION")
    print("="*40)
    level_counts = unique_templates['level'].value_counts()
    print(level_counts)

    # ---------------------------------------------------------
    # STEP 3: Branch vs. Difficulty Split (For your Table)
    # ---------------------------------------------------------
    # Group by Branch and Level
    branch_stats = unique_templates.groupby(['branch', 'level']).size().unstack(fill_value=0)
    
    # Add a 'Total' column for convenience
    branch_stats['Total'] = branch_stats.sum(axis=1)

    print("\n" + "="*40)
    print("STATISTICS BY BRANCH (Train/Test Split Info)")
    print("="*40)
    print(branch_stats)

    # ---------------------------------------------------------
    # STEP 4: Detailed Structural Stats (Domain/Area counts)
    # ---------------------------------------------------------
    # This helps fill the "Structure" part of your table
    print("\n" + "="*40)
    print("STRUCTURAL STATISTICS")
    print("="*40)
    
    # Count unique domains and areas per branch
    structure_stats = unique_templates.groupby('branch').agg({
        'domain': 'nunique',
        'area': 'nunique',
        'id': 'count'
    }).rename(columns={'id': 'Total Templates'})
    
    print(structure_stats)

    return unique_templates

# --- USAGE ---
# Replace this path with your actual file location
file_path = r"Evaluation Results EngChain FINAL/Frontier Proprietary LLMs/gpt-4.1_evaluation_results_FINAL.jsonl"

if os.path.exists(file_path):
    df_unique = analyze_dataset_stats(file_path)
else:
    print(f"File not found: {file_path}")
