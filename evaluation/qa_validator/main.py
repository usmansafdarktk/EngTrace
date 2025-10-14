import os
import inspect
import time
import json
from importlib import util
from . import llm_reviewer

# The root directory containing all engineering branch templates
TEMPLATES_ROOT_DIR = "data/templates/branches"

def discover_template_files(root_dir: str):
    """Finds all Python template files recursively."""
    template_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".py") and filename not in ["__init__.py", "constants.py"]:
                template_files.append(os.path.join(dirpath, filename))
    return template_files

def import_templates_from_path(file_path: str):
    """Dynamically imports a module and finds all template functions."""
    module_name = os.path.splitext(os.path.basename(file_path))[0]
    spec = util.spec_from_file_location(module_name, file_path)
    if not spec or not spec.loader:
        return []
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    template_functions = []
    for name, func in inspect.getmembers(module, inspect.isfunction):
        if name.startswith("template_"):
            template_functions.append(func)
            
    return template_functions

def run_validation_on_all_templates():
    """
    Discovers, validates, and reports on all EngChain templates.
    """
    llm_reviewer.setup_api_key()
    template_files = discover_template_files(TEMPLATES_ROOT_DIR)
    
    if not template_files:
        print(f"No template files found in '{TEMPLATES_ROOT_DIR}'. Exiting.")
        return

    print(f"Found {len(template_files)} template files to validate.")
    
    approved_templates = []
    flagged_templates = []
    full_report = {}

    for i, file_path in enumerate(template_files):
        print("\n" + "="*80)
        print(f"Processing file {i+1}/{len(template_files)}: {file_path}")
        print("="*80)

        template_functions = import_templates_from_path(file_path)
        if not template_functions:
            print(f"No template functions found in {file_path}. Skipping.")
            continue

        with open(file_path, 'r', encoding='utf-8') as f:
            template_code = f.read()

        for template_func in template_functions:
            template_name = template_func.__name__
            print(f"\n--- Validating template: {template_name} ---")
            
            instances = [template_func() for _ in range(3)]
            evaluation = llm_reviewer.validate_template_with_llm(template_code, instances, template_name)
            
            if not evaluation:
                print(f"Validation failed for {template_name}. Flagging for human review.")
                flagged_templates.append(f"{file_path} -> {template_name}")
                full_report[template_name] = {"status": "FAILED_TO_EVALUATE", "details": "No valid response from LLM."}
                continue

            threshold = 4
            scores = [
                evaluation.get('physical_plausibility_score', 0),
                evaluation.get('mathematical_correctness_score', 0),
                evaluation.get('pedagogical_clarity_score', 0)
            ]
            
            is_approved = all(s >= threshold for s in scores)
            
            report_entry = {
                "file_path": file_path,
                "status": "APPROVED" if is_approved else "FLAGGED",
                "details": evaluation
            }
            full_report[template_name] = report_entry

            if is_approved:
                print(f"Result: APPROVED")
                approved_templates.append(f"{file_path} -> {template_name}")
            else:
                print(f"Result: FLAGGED FOR HUMAN REVIEW")
                flagged_templates.append(f"{file_path} -> {template_name}")
            
            # Add a small delay to respect potential API rate limits
            time.sleep(2) 

    # Define the output directory and create it if it doesn't exist
    output_dir = os.path.join("evaluation", "qa_validator", "results")
    os.makedirs(output_dir, exist_ok=True)

    # Define the full path for the report file
    report_filename = os.path.join(output_dir, "qa_summary_report.json")
    
    # Save the detailed report to the specified path
    with open(report_filename, 'w') as f:
        json.dump(full_report, f, indent=4)
    print(f"\nDetailed report saved to '{report_filename}'")

    # Print the final summary
    print("\n" + "="*80)
    print("AI-ASSISTED QA FINAL SUMMARY")
    print("="*80)
    print(f"\nAPPROVED TEMPLATES ({len(approved_templates)}):")
    if approved_templates:
        for t in approved_templates:
            print(f"- {t}")
    else:
        print("None")
        
    print(f"\nFLAGGED FOR HUMAN REVIEW ({len(flagged_templates)}):")
    if flagged_templates:
        for t in flagged_templates:
            print(f"- {t}")
    else:
        print("None")
    print("="*80)


if __name__ == "__main__":
    run_validation_on_all_templates()
