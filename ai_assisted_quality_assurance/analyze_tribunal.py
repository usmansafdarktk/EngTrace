import json
import numpy as np
import pandas as pd
from pathlib import Path

# ==========================================
# Configuration
# ==========================================
SCRIPT_DIR = Path(__file__).parent
ROOT_DIR = SCRIPT_DIR.parent

# Input/Output Paths
INPUT_DIR = ROOT_DIR / "ai_assistance_data" / "qa_results"
OUTPUT_DIR = ROOT_DIR / "ai_assistance_data" / "analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MODELS = ["openai", "anthropic", "google"]
METRICS = ["physical_plausibility_score", "mathematical_correctness_score", "pedagogical_clarity_score"]

# ==========================================
# 1. Data Loading
# ==========================================
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
                # Ensure valid response structure
                if entry and "response" in entry and isinstance(entry["response"], dict):
                    # Check if response contains the scores (filters out API errors)
                    if "physical_plausibility_score" in entry["response"]:
                        tid = entry["template_id"]
                        data[tid] = entry
            except json.JSONDecodeError:
                continue
    return data

# ==========================================
# 2. Consensus Logic (The "Tribunal" Engine)
# ==========================================
def calculate_consensus(votes):
    """
    Applies PoLL aggregation logic:
    1. Scores: Median (Robust to outliers) for ALL 3 metrics.
    2. Flag: Majority Vote (Conservative safety).
    3. Disagreement: Max Standard Deviation across any metric.
    """
    
    # Extract scores for all metrics
    scores = {m: [] for m in METRICS}
    flags = []
    confidences = []
    
    for vote in votes:
        resp = vote["response"]
        for m in METRICS:
            scores[m].append(resp.get(m, 0))
        
        flags.append(resp.get("human_review_flag", False))
        confidences.append(resp.get("confidence_score", 0))

    consensus = {}
    
    # A. Calculate Medians for ALL metrics
    medians = []
    for m in METRICS:
        val = int(np.median(scores[m]))
        consensus[m] = val
        medians.append(val)
        
    # Track the lowest score across all 3 dimensions (for the Gatekeeper check)
    consensus["min_score"] = min(medians)
        
    # B. Majority Vote for Flag (True/False)
    # If 2 or more models (>= 66%) say "True" (Bad), we mark it True.
    flag_count = sum(flags)
    consensus["consensus_flag"] = flag_count >= 2
    
    # C. Metadata
    consensus["avg_confidence"] = round(np.mean(confidences), 2)
    consensus["flag_agreement"] = f"{flag_count}/{len(votes)}" 
    
    # D. Disagreement Metric (Sigma)
    # We calculate SD for ALL metrics and take the WORST one.
    # If models agree on Physics but disagree heavily on Math, we want to know.
    std_devs = [np.std(scores[m]) for m in METRICS]
    consensus["disagreement_score"] = round(max(std_devs), 2)
    
    return consensus

def categorize_template(consensus):
    """
    Buckets templates based on the Tribunal's verdict.
    """
    # 🔴 Critical Failure:
    # 1. Majority of AI models raised a Red Flag
    # OR
    # 2. The median score for ANY dimension (Phys, Math, Ped) is < 4
    if consensus["consensus_flag"] or consensus["min_score"] < 4:
        return "Critical Failure"
    
    # 🟡 Controversial:
    # The template passed, but there is high disagreement (SD > 0.5)
    # on at least one dimension.
    if consensus["disagreement_score"] > 0.5:
        return "Controversial"
        
    # 🟢 Pass:
    # High scores (>4) everywhere, no flags, and consensus is stable.
    return "Pass"

# ==========================================
# 3. Main Execution
# ==========================================
def main():
    print("=== Starting AI Tribunal Analysis (PoLL Methodology) ===")
    print(f"Reading from: {INPUT_DIR}")
    
    # 1. Load Data
    raw_data = {model: load_jsonl(INPUT_DIR / f"results_{model}.jsonl") for model in MODELS}
    
    # Get all unique template IDs
    all_tids = set().union(*[d.keys() for d in raw_data.values()])
    print(f"Found data for {len(all_tids)} unique templates.")
    
    analysis_rows = []
    priority_list = []
    
    # 2. Process each template
    for tid in all_tids:
        # Gather votes
        votes = []
        missing_models = []
        
        for model in MODELS:
            if tid in raw_data[model]:
                votes.append(raw_data[model][tid])
            else:
                missing_models.append(model)
        
        if not votes:
            continue
            
        # Get Metadata (from first vote)
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
            # Key Scores
            "Median Physical Score": stats["physical_plausibility_score"],
            "Median Math Score": stats["mathematical_correctness_score"],
            "Median Pedagogy Score": stats["pedagogical_clarity_score"],
            # Quality Checks
            "Disagreement (SD)": stats["disagreement_score"],
            "Vote Split (Flags)": stats["flag_agreement"],
            "Missing Models": ", ".join(missing_models) if missing_models else "None"
        }
        
        # Add individual model details
        for model in MODELS:
            if tid in raw_data[model]:
                resp = raw_data[model][tid]["response"]
                row[f"{model}_Score"] = resp.get("physical_plausibility_score") # Can change to avg of 3 if needed
                row[f"{model}_Flag"] = resp.get("human_review_flag")
                row[f"{model}_Reason"] = resp.get("explanation")
            else:
                row[f"{model}_Score"] = "N/A"
                row[f"{model}_Flag"] = "N/A"
                row[f"{model}_Reason"] = "N/A"

        analysis_rows.append(row)
        
        # Add to Priority List if not clean Pass
        if category != "Pass":
            priority_list.append({
                "id": tid,
                "category": category,
                "branch": meta.get("branch", "Unknown"),
                "area": meta.get("area", "Unknown"),
                "issues": {
                    "consensus_flag": stats["consensus_flag"],
                    "min_score": stats["min_score"],
                    "disagreement": stats["disagreement_score"],
                    "votes": stats["flag_agreement"],
                    "openai_reason": row.get("openai_Reason"),
                    "anthropic_reason": row.get("anthropic_Reason"),
                    "google_reason": row.get("google_Reason")
                }
            })

    # 3. Export Results
    if not analysis_rows:
        print("No valid data rows generated.")
        return

    # CSV Summary
    df = pd.DataFrame(analysis_rows)
    
    # Reorder columns
    base_cols = ["Template ID", "Category", "Consensus Flag", 
                 "Median Physical Score", "Median Math Score", "Median Pedagogy Score", 
                 "Disagreement (SD)", "Vote Split (Flags)"]
    remaining_cols = [c for c in df.columns if c not in base_cols]
    df = df[base_cols + remaining_cols]
    
    csv_path = OUTPUT_DIR / "tribunal_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f">> Summary CSV saved to: {csv_path}")
    
    # Priority JSON
    json_path = OUTPUT_DIR / "qa_priority_list.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(priority_list, f, indent=4)
    print(f">> Priority JSON ({len(priority_list)} items) saved to: {json_path}")
    
    # 4. Print Verdict
    print("\n=== Tribunal Verdict Summary ===")
    print(df["Category"].value_counts())

if __name__ == "__main__":
    main()
