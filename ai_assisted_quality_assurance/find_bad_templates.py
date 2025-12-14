import json
from pathlib import Path
from collections import defaultdict

#  Configuration 
SCRIPT_DIR = Path(__file__).parent
ROOT_DIR = SCRIPT_DIR.parent
INPUT_DIR = ROOT_DIR / "ai_assistance_data" / "qa_results"
MODELS = ["openai", "anthropic", "google"]

def load_jsonl(filepath):
    """Loads a JSONL file and returns a dict of {template_id: human_review_flag}."""
    flags = {}
    if not filepath.exists():
        print(f"Warning: {filepath} not found.")
        return flags
        
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            try:
                entry = json.loads(line)
                if "response" in entry and isinstance(entry["response"], dict):
                    # Get the flag (default to False if missing)
                    flag = entry["response"].get("human_review_flag", False)
                    tid = entry["template_id"]
                    flags[tid] = flag
            except json.JSONDecodeError:
                continue
    return flags

def main():
    print(" Identifying Bad Templates (2+ Flags) ")
    
    # 1. Load Flags from all models
    model_flags = {}
    all_tids = set()
    
    for model in MODELS:
        path = INPUT_DIR / f"results_{model}.jsonl"
        model_flags[model] = load_jsonl(path)
        all_tids.update(model_flags[model].keys())
        
    bad_templates = []
    
    # 2. Count Flags per Template
    for tid in all_tids:
        flag_count = 0
        voters = []
        
        for model in MODELS:
            # Check if this model flagged this template
            if model_flags[model].get(tid, False) is True:
                flag_count += 1
                voters.append(model)
                
        # 3. Filter Logic: If 2 or more models flagged it
        if flag_count >= 2:
            bad_templates.append({
                "id": tid,
                "count": flag_count,
                "flagged_by": voters
            })
            
    # 4. Output Results
    if not bad_templates:
        print("No bad templates found! All clear.")
    else:
        print(f"Found {len(bad_templates)} bad templates requiring immediate fix:\n")
        for t in bad_templates:
            print(f"🔴 {t['id']}")
            print(f"   Flags: {t['count']}/3 ({', '.join(t['flagged_by'])})")
            print("   ")

if __name__ == "__main__":
    main()
