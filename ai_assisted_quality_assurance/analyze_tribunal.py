import json
import numpy as np
import pandas as pd
from pathlib import Path

#  Configuration (Relative Paths) 
SCRIPT_DIR = Path(__file__).parent
ROOT_DIR = SCRIPT_DIR.parent

# Data Paths
INPUT_DIR = ROOT_DIR / "ai_assistance_data" / "qa_results"
OUTPUT_DIR = ROOT_DIR / "ai_assistance_data" / "analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MODELS = ["openai", "anthropic", "google"]
METRICS = ["physical_plausibility_score", "mathematical_correctness_score", "pedagogical_clarity_score"]

def load_jsonl(filepath):
    """Loads a JSONL file into a dictionary keyed by template_id."""
    data = {}
    if not filepath.exists():
        print(f"Warning: {filepath} not found.")
        return data
        
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            try:
                entry = json.loads(line)
                # Ensure we have a valid response structure
                if entry and "response" in entry and isinstance(entry["response"], dict):
                    # Check if the response contains the scores (filters out API errors)
                    if "physical_plausibility_score" in entry["response"]:
                        tid = entry["template_id"]
                        data[tid] = entry
            except json.JSONDecodeError:
                continue
    return data

def calculate_consensus(votes):
    """
    Applies PoLL aggregation logic verified by Verga et al. (2024):
    1. Scores: Median (Robust to outliers) 
    2. Flag: Majority Vote (Conservative safety) [cite: 85, 93]
    """
    
    # Extract scores
    scores = {m: [] for m in METRICS}
    flags = []
    confidences = []
    
    for vote in votes:
        resp = vote["response"]
        for m in METRICS:
            # Default to 0 if missing, though schema ensures it shouldn't be
            scores[m].append(resp.get(m, 0))
        
        # Capture the boolean flag
        flags.append(resp.get("human_review_flag", False))
        confidences.append(resp.get("confidence_score", 0))

    # 1. Aggregation Logic
    consensus = {}
    
    # Median for Ordinal Scores (1-5)
    for m in METRICS:
        consensus[m] = int(np.median(scores[m]))
        
    # Majority Vote for Flag (True/False)
    # If 2 or more models (>= 66%) say "True" (Bad), we mark it True.
    flag_count = sum(flags)
    consensus["consensus_flag"] = flag_count >= 2
    
    # Avg Confidence (Just for metadata)
    consensus["avg_confidence"] = round(np.mean(confidences), 2)
    
    # Agreement Stats
    consensus["flag_agreement"] = f"{flag_count}/{len(votes)}" # e.g., "2/3"
    
    # Disagreement Metric: Standard Deviation of the Physical Plausibility score
    # High SD (> 0.5) implies models disagree significantly (e.g., scores 5, 2, 5)
    consensus["disagreement_score"] = round(np.std(scores["physical_plausibility_score"]), 2)
    
    return consensus

def categorize_template(consensus):
    """
    Buckets templates based on the Tribunal's verdict.
    """
    # 🔴 Critical: Majority Flagged it OR Median Score is failing (<4)
    # This filters out templates that are structurally unsound.
    if consensus["consensus_flag"] or consensus["physical_plausibility_score"] < 4:
        return "Critical Failure"
    
    # 🟡 Controversial: High Disagreement (SD > 0.5)
    # These are passing templates where one judge strongly disagreed.
    if consensus["disagreement_score"] > 0.5:
        return "Controversial"
        
    # 🟢 Pass: Clean record (High Scores + Consensus)
    return "Pass"

def main():
    print(" Starting AI Tribunal Analysis (PoLL Methodology) ")
    print(f"Reading from: {INPUT_DIR}")
    
    # 1. Load Data
    raw_data = {model: load_jsonl(INPUT_DIR / f"results_{model}.jsonl") for model in MODELS}
    
    # Get all unique template IDs across all files
    all_tids = set().union(*[d.keys() for d in raw_data.values()])
    print(f"Found data for {len(all_tids)} unique templates.")
    
    analysis_rows = []
    priority_list = []
    
    # 2. Process each template
    for tid in all_tids:
        # Gather votes for this template
        votes = []
        missing_models = []
        
        for model in MODELS:
            if tid in raw_data[model]:
                votes.append(raw_data[model][tid])
            else:
                missing_models.append(model)
        
        # Skip if no valid votes found
        if not votes:
            continue
            
        # Get Metadata from first available vote
        meta = votes[0]
        
        # Calculate Consensus
        stats = calculate_consensus(votes)
        category = categorize_template(stats)
        
        # Prepare CSV Row
        row = {
            "Template ID": tid,
            "Branch": meta.get("branch", "Unknown"),
            "Area": meta.get("area", "Unknown"),
            "Category": category,
            "Consensus Flag": stats["consensus_flag"],
            "Median Physical Score": stats["physical_plausibility_score"],
            "Median Math Score": stats["mathematical_correctness_score"],
            "Median Pedagogy Score": stats["pedagogical_clarity_score"],
            "Disagreement (SD)": stats["disagreement_score"],
            "Vote Split (Flags)": stats["flag_agreement"],
            "Missing Models": ", ".join(missing_models) if missing_models else "None"
        }
        
        # Add individual model details for deeper inspection
        for model in MODELS:
            if tid in raw_data[model]:
                resp = raw_data[model][tid]["response"]
                row[f"{model}_Score"] = resp.get("physical_plausibility_score")
                row[f"{model}_Flag"] = resp.get("human_review_flag")
                row[f"{model}_Reason"] = resp.get("explanation")
            else:
                row[f"{model}_Score"] = "N/A"
                row[f"{model}_Flag"] = "N/A"
                row[f"{model}_Reason"] = "N/A"

        analysis_rows.append(row)
        
        # Add to Priority List (JSON) if it's not a clean Pass
        if category != "Pass":
            priority_list.append({
                "id": tid,
                "category": category,
                "branch": meta.get("branch", "Unknown"),
                "area": meta.get("area", "Unknown"),
                "issues": {
                    "consensus_flag": stats["consensus_flag"],
                    "median_score": stats["physical_plausibility_score"],
                    "disagreement": stats["disagreement_score"],
                    "votes": stats["flag_agreement"],
                    "openai_reason": row.get("openai_Reason"),
                    "anthropic_reason": row.get("anthropic_Reason"),
                    "google_reason": row.get("google_Reason")
                }
            })

    # 3. Export Results
    if not analysis_rows:
        print("No valid data rows generated. Check your JSONL files.")
        return

    # CSV Summary
    df = pd.DataFrame(analysis_rows)
    
    # Reorder columns for readability (Put summary stats first)
    cols = ["Template ID", "Category", "Consensus Flag", "Median Physical Score", "Vote Split (Flags)"] + \
           [c for c in df.columns if c not in ["Template ID", "Category", "Consensus Flag", "Median Physical Score", "Vote Split (Flags)"]]
    df = df[cols]
    
    csv_path = OUTPUT_DIR / "tribunal_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"Summary CSV saved to: {csv_path}")
    
    # Priority JSON
    json_path = OUTPUT_DIR / "qa_priority_list.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(priority_list, f, indent=4)
    print(f"Priority JSON ({len(priority_list)} items) saved to: {json_path}")
    
    # 4. Print Verdict
    print("\n Tribunal Verdict Summary ")
    print(df["Category"].value_counts())

if __name__ == "__main__":
    main()
