import os

# Items to skip (files or directories)
IGNORE = {"__pycache__", "directory_listing.py", ".git", ".gitignore", "README.md"}

def list_dir_structure(start_path='.', indent=''):
    """Recursively print the directory structure starting from start_path,
    skipping ignored folders and files."""
    try:
        items = sorted(os.listdir(start_path))
    except PermissionError:
        print(f"{indent}[Permission Denied]")
        return

    for i, item in enumerate(items):
        if item in IGNORE:
            continue  # skip ignored files/folders

        path = os.path.join(start_path, item)
        connector = "└── " if i == len(items) - 1 else "├── "
        print(indent + connector + item)

        if os.path.isdir(path):
            extension = "    " if i == len(items) - 1 else "│   "
            list_dir_structure(path, indent + extension)


if __name__ == "__main__":
    root = "."  # you can change this to any directory path
    print(root)
    list_dir_structure(root)
