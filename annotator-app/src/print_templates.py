import sys
from template_loader import (
    get_areas,
    get_template_files,
    load_template_functions,
    get_source_code
)

# Configuration
TARGET_BRANCH = "chemical_engineering"

def main():
    print(f" Generating Report for: {TARGET_BRANCH}\n")
    
    # 1. Get all areas in the branch
    areas = get_areas(TARGET_BRANCH)
    if not areas:
        print(f"No areas found for {TARGET_BRANCH}. Check your path configuration.")
        return

    total_templates = 0

    for area in areas:
        # 2. Get all files in this area
        files = get_template_files(TARGET_BRANCH, area)
        
        for filename in files:
            try:
                # 3. Load functions from the file
                funcs = load_template_functions(TARGET_BRANCH, area, filename)
                
                for func_name, func_obj in funcs:
                    total_templates += 1
                    print("=" * 80)
                    print(f"TEMPLATE #{total_templates}: {func_name}")
                    print(f"LOCATION: {area} / {filename}.py")
                    print("=" * 80)

                    # A. Print Source Code (Using the new loader function)
                    print("\n[PYTHON SOURCE]:")
                    print("-" * 20)
                    code = get_source_code(func_obj)
                    print(code.strip())
                    print("-" * 20)

                    # B. Generate One Instance
                    print("\n[GENERATED INSTANCE]:")
                    try:
                        q, s = func_obj()
                        print("-" * 20)
                        print(f"QUESTION:\n{q}\n")
                        print(f"SOLUTION:\n{s}")
                        print("-" * 20)
                    except Exception as e:
                        print(f"Error generating instance: {e}")

                    print("\n\n")

            except Exception as e:
                print(f"Failed to process {filename}: {e}")

    print(f"Done. Processed {total_templates} templates.")

if __name__ == "__main__":
    main()