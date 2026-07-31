import sys
import traceback
from src.template_loader import (
    get_branches, 
    get_areas, 
    get_template_files, 
    load_template_functions
)

def main():
    print("--- Starting Template Loader Test (V2) ---")

    # 1. Test Branches
    branches = get_branches()
    print(f"\nBranches Found: {branches}")
    
    if not branches:
        print("No branches found. Check your directory structure.")
        return

    # 2. Test Areas (using the first branch found, likely 'chemical_engineering')
    test_branch = "chemical_engineering" 
    if test_branch not in branches:
        test_branch = branches[0]
        
    areas = get_areas(test_branch)
    print(f"Areas in '{test_branch}': {areas}")

    if not areas:
        print(f"No areas found in {test_branch}.")
        return

    # 3. Test Files (using 'reaction_kinetics' if available)
    test_area = "reaction_kinetics"
    if test_area not in areas:
        test_area = areas[0]

    files = get_template_files(test_branch, test_area)
    print(f"Files in '{test_branch}/{test_area}': {files}")
    
    # We specifically want to test 'mole_balances' if it exists
    target_file = "mole_balances"
    if target_file not in files:
        print(f"WARNING: 'mole_balances.py' not found. Using first available file: {files[0]}")
        target_file = files[0]

    # 4. Test Function Extraction (The Logic Upgrade)
    print(f"\nInspecting file: '{target_file}.py'...")
    try:
        # This loads the file and finds all def template_...():
        templates = load_template_functions(test_branch, test_area, target_file)
        
        if not templates:
            print("No functions starting with 'template_' found in this file.")
            return
            
        print(f"Found {len(templates)} template functions:")
        for name, func in templates:
            print(f"   - {name}")

        # 5. Test Execution (Run the first function found)
        first_func_name, first_func_obj = templates[0]
        print(f"\nExecuting '{first_func_name}'...")
        
        question, solution = first_func_obj()
        
        print("\n--- [Output Preview] ---")
        print(f"Question: {question[:100]}...") # Print first 100 chars
        print(f"Solution: {solution[:100]}...")
        print("------------------------")
        print("SUCCESS: Template loaded and executed correctly.")

    except Exception as e:
        print(f"CRITICAL ERROR during loading/execution:")
        traceback.print_exc()

if __name__ == "__main__":
    main()
