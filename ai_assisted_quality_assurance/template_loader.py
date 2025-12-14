import sys
import inspect
import importlib.util
from pathlib import Path

# Define root path relative to this script
BASE_DIR = Path(__file__).parent.parent / "data" / "templates" / "branches"

def get_branches():
    """Returns a list of available engineering branches."""
    if not BASE_DIR.exists():
        return []
    return sorted([d.name for d in BASE_DIR.iterdir() if d.is_dir() and not d.name.startswith("__")])

def get_areas(branch_name):
    """Returns the sub-areas (directories) for a given branch."""
    branch_path = BASE_DIR / branch_name
    if not branch_path.exists():
        return []
    return sorted([d.name for d in branch_path.iterdir() if d.is_dir() and not d.name.startswith("__")])

def get_template_files(branch_name, area_name):
    """Returns the python filenames (modules) in an area."""
    area_path = BASE_DIR / branch_name / area_name
    if not area_path.exists():
        return []
    
    files = [
        f.stem for f in area_path.glob("*.py")
        if f.name not in ["constants.py", "__init__.py"] and not f.name.startswith("__")
    ]
    return sorted(files)

def load_template_functions(branch, area, filename):
    """
    Loads a specific file and returns a list of its template functions.
    Returns: List of tuples (function_name, function_object)
    """
    file_path = BASE_DIR / branch / area / f"{filename}.py"
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    # Dynamic import
    spec = importlib.util.spec_from_file_location(filename, file_path)
    module = importlib.util.module_from_spec(spec)
    
    # Add branch folder to path so imports like 'from constants import...' work
    sys.path.append(str(file_path.parent.parent)) 
    
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        return [] # Return empty if module crashes on load

    # Extract all functions that start with 'template_'
    templates = []
    for name, obj in inspect.getmembers(module):
        if inspect.isfunction(obj) and name.startswith("template_"):
            templates.append((name, obj))
            
    return sorted(templates, key=lambda x: x[0])

def get_source_code(branch, area, filename):
    """Returns the raw source code of the entire file."""
    file_path = BASE_DIR / branch / area / f"{filename}.py"
    if not file_path.exists():
        return "# Error: File not found."
    with open(file_path, "r", encoding="utf-8") as f:
        return f.read()
