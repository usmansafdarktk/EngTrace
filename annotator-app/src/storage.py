# src/storage.py
import json
import os
from datetime import datetime
from pathlib import Path

# Define where to save reviews
REVIEW_DIR = Path(__file__).parent.parent / "reviews"
REVIEW_DIR.mkdir(exist_ok=True)

def save_review(annotator_name, branch, area, template_name, scores, decision, feedback):
    """
    Appends a single review to an annotator-specific JSONL file.
    """
    # Sanitize name for filename (replace spaces with underscores)
    safe_name = "".join([c for c in annotator_name if c.isalnum() or c in (' ', '_')]).strip().replace(" ", "_")
    safe_branch = branch.replace(" ", "_")
    
    # Filename: e.g., reviews/John_Doe_chemical_engineering.jsonl
    filename = REVIEW_DIR / f"{safe_name}_{safe_branch}.jsonl"
    
    review_data = {
        "timestamp": datetime.now().isoformat(),
        "annotator_id": annotator_name,
        "branch": branch,
        "area": area,
        "template": template_name,
        "scores": {
            "physical_plausibility": scores[0],
            "mathematical_correctness": scores[1],
            "pedagogical_clarity": scores[2]
        },
        "decision": decision,
        "feedback": feedback
    }

    # Append to file immediately (Data Safety)
    with open(filename, "a", encoding="utf-8") as f:
        f.write(json.dumps(review_data) + "\n")
        
    return str(filename)
