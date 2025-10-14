import pickle
import os
import pandas as pd
import argparse 


# Create a parser object to handle command-line arguments
parser = argparse.ArgumentParser(
    description="Load and display summary statistics from .pkl files in a specified directory."
)
# Add a required argument for the results directory
parser.add_argument(
    "results_dir", 
    type=str,
    help="The path to the directory containing the .pkl result files."
)
# Parse the arguments provided by the user
args = parser.parse_args()


# Configuration 
RESULTS_DIR = args.results_dir
COLUMN_NAMES = ['final_answer_match', 'recall', 'precision', 'step_f1', 'bertscore', 'rouge2']


#  2. Find and Load All Pickle Files 
if not os.path.isdir(RESULTS_DIR):
    print(f"Error: The directory '{RESULTS_DIR}' was not found.")
else:
    # Find all files ending with .pkl in the specified directory.
    pickle_files = [f for f in os.listdir(RESULTS_DIR) if f.endswith('.pkl')]

    if not pickle_files:
        print(f"No pickle (.pkl) files found in '{RESULTS_DIR}'.")
    else:
        # Loop through each found file.
        for file_name in pickle_files:
            file_path = os.path.join(RESULTS_DIR, file_name)
            
            try:
                with open(file_path, 'rb') as f:
                    loaded_data = pickle.load(f)
                
                print(f" Summary for: {file_name} ")
                
                # Convert to pandas DataFrame for a clean summary.
                # Assumes data is a dictionary like {'model_name': {'category': [scores]}}
                if isinstance(loaded_data, dict) and loaded_data:
                    model_name = list(loaded_data.keys())[0]
                    df = pd.DataFrame.from_dict(loaded_data[model_name], orient='index')
                    df.columns = COLUMN_NAMES
                    df.index.name = 'Category'
                    print(df.round(3))
                else:
                    print(loaded_data) # Fallback for unexpected data format
                
                print("-" * (len(file_name) + 18) + "\n")

            except Exception as e:
                print(f"Error loading or processing {file_name}: {e}\n")
