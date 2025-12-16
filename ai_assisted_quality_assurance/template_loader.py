import sys
import inspect
import importlib.util
from pathlib import Path

#  PATH CONFIGURATION 
# 1. Get the folder where this script is (e.g., EngChain/src)
CURRENT_DIR = Path(__file__).resolve().parent
# 2. Get the project root (e.g., EngChain) 
PROJECT_ROOT = CURRENT_DIR.parent
# 3. Define where templates live
BASE_DIR = PROJECT_ROOT / "data" / "templates" / "branches"

def get_branches():
    if not BASE_DIR.exists(): return []
    return sorted([d.name for d in BASE_DIR.iterdir() if d.is_dir() and not d.name.startswith("__")])

def get_areas(branch_name):
    branch_path = BASE_DIR / branch_name
    if not branch_path.exists(): return []
    return sorted([d.name for d in branch_path.iterdir() if d.is_dir() and not d.name.startswith("__")])

def get_template_files(branch_name, area_name):
    area_path = BASE_DIR / branch_name / area_name
    if not area_path.exists(): return []
    return sorted([f.stem for f in area_path.glob("*.py") if f.name not in ["constants.py", "__init__.py"] and not f.name.startswith("__")])

def load_template_functions(branch, area, filename):
    file_path = BASE_DIR / branch / area / f"{filename}.py"
    if not file_path.exists(): raise FileNotFoundError(f"File not found: {file_path}")

    # Add Project Root to sys.path 
    # This ensures 'from data...' imports work, preventing the silent crash
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.append(str(PROJECT_ROOT))

    spec = importlib.util.spec_from_file_location(filename, file_path)
    module = importlib.util.module_from_spec(spec)
    
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        print(f"!!! Error loading module '{filename}': {e}")
        return []

    templates = []
    for name, obj in inspect.getmembers(module):
        if inspect.isfunction(obj) and name.startswith("template_"):
            templates.append((name, obj))
    return sorted(templates, key=lambda x: x[0])

def get_source_code(func_object):
    """Extracts ONLY the function code, not the whole file."""
    try:
        return inspect.getsource(func_object)
    except OSError:
        return "# Error: Source not available."


if __name__ == "__main__":
    print(" 1. Testing get_branches ")
    branches = get_branches()
    print(f"Branches found: {branches}")

    if branches:
        # Test with the first available branch
        selected_branch = branches[0]
        print(f"\n 2. Testing get_areas for '{selected_branch}' ")
        areas = get_areas(selected_branch)
        print(f"Areas found: {areas}")

        if areas:
            # Test with the first available area
            selected_area = areas[0]
            print(f"\n 3. Testing get_template_files for '{selected_area}' ")
            files = get_template_files(selected_branch, selected_area)
            print(f"Files found: {files}")

            if files:
                # Test with the first available file
                selected_file = files[1]
                
                print(f"\n 4. Testing load_template_functions for '{selected_file}' ")
                funcs = load_template_functions(selected_branch, selected_area, selected_file)
                func_names = [f[0] for f in funcs]
                print(f"Functions found: {func_names}")
                
                if funcs:
                    # Test getting source code for the FIRST function found
                    first_func_name, first_func_obj = funcs[2]
                    print(f"\n 5. Testing get_source_code for '{first_func_name}' ")
                    
                    code = get_source_code(first_func_obj)
                    print(f"Source code preview:\n{code}...") 
                else:
                    print("\n 5. Skipped (No functions found to extract code from) ")
