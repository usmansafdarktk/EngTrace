# run_inference.py

import os
import json
import time
from tqdm import tqdm
from dotenv import load_dotenv
import google.generativeai as genai


#  1. Configuration 
load_dotenv()


# Configure the API client
API_KEY = os.getenv("API_KEY")
if not API_KEY:
    raise ValueError("API_KEY not found in .env file or environment variables.")
genai.configure(api_key=API_KEY)


# Define model and file paths
MODEL_NAME =  os.getenv("MODEL_NAME")
INPUT_DIR = "../testset"
OUTPUT_FILE = f"inference_results/{MODEL_NAME}_inference_results.jsonl"


# Define the prompt template for the model
PROMPT_TEMPLATE = """
You are an expert engineer. Solve the following problem by providing a detailed, structured solution. Use the exact headings and formatting provided below.

## Given
List all known variables and their values with units.

## Find
State the variable(s) to be calculated.

## Formulae
Write down all necessary governing equations before substituting any values.

## Solution
Provide a step-by-step calculation. Each step must start on a new line and be formatted exactly as '**Step X:**', where X is the step number. Show the substitution of values into the formulae clearly.

## Final Answer
State the final numerical result with its units in the format: **Answer:** [value] [units]

Problem:
{question}
"""


#  2. Data Loading Function 
def load_all_problems(directory: str) -> list:
    """
    Walks through the nested directory structure, finds all .jsonl files,
    and loads all problems into a single list.
    """
    all_problems = []
    print(f"Loading problems from '{directory}'...")
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.jsonl'):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    for line in f:
                        all_problems.append(json.loads(line))
                        
    print(f"Successfully loaded {len(all_problems)} problems.")
    return all_problems


#  3. Main Inference Logic 
if __name__ == "__main__":
    # Load all problems from the testset directory
    problems = load_all_problems(INPUT_DIR)

    if not problems:
        print(f"Error: No problems found in '{INPUT_DIR}'. Please check the path.")
    else:
        # Initialize the generative model
        model = genai.GenerativeModel(MODEL_NAME)
        print(f"Initialized model: {MODEL_NAME}")

        # Open the output file in write mode to start fresh
        with open(OUTPUT_FILE, "w", encoding='utf-8') as f_out:
            print(f"Starting inference... Results will be saved to '{OUTPUT_FILE}'")
            
            # Create a tqdm iterator to get a handle on the progress bar
            progress_bar = tqdm(problems, desc="Initializing Inference")

            for problem in progress_bar:
                
                # Update the progress bar's description 
                branch = problem.get('branch', 'unknown_branch')
                problem_id = problem.get('id', 'unknown_id')
                progress_bar.set_description(f"Processing '{problem_id}' from '{branch}'")
                
                prompt = PROMPT_TEMPLATE.format(question=problem['question'])
                
                try:
                    # Call the Gemini API
                    response = model.generate_content(prompt)
                    
                    # Add the model's generation to the problem dictionary
                    problem['generation'] = response.text
                    
                except Exception as e:
                    # Handle potential API errors
                    # Use tqdm.write to print errors without breaking the bar
                    tqdm.write(f"\nAn error occurred for problem ID {problem_id}: {e}")
                    problem['generation'] = f"ERROR: {e}"
                    time.sleep(5)
                
                # Write the result to the output file immediately
                f_out.write(json.dumps(problem) + '\n')

        print(f"\nInference complete. All {len(problems)} results saved to '{OUTPUT_FILE}'.")
