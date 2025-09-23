import os
import json
from tqdm import tqdm
from dotenv import load_dotenv

# Import the main evaluation function we created earlier
from evaluation_methods import evaluate_trace_eng


#  1. Configuration 
load_dotenv()

# Get the model name from the .env file to find the correct results file
MODEL_NAME = os.getenv("MODEL_NAME") 

# Define the input file (from inference) and the output file (for scores)
# Assumes this script is in the 'evaluation' directory
INPUT_FILE = f"inference_results/{MODEL_NAME}_inference_results.jsonl"
OUTPUT_FILE = f"evaluation_results/{MODEL_NAME}_evals.jsonl"


#  2. Main Evaluation Logic 
if __name__ == "__main__":
    
    # Check if the input file from the inference step exists
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Input file not found at '{INPUT_FILE}'")
        print("Please make sure you have run the inference step first.")
    else:
        print(f"Reading inference results from '{INPUT_FILE}'...")
        
        # Open the input and output files
        with open(INPUT_FILE, 'r', encoding='utf-8') as f_in, \
             open(OUTPUT_FILE, 'w', encoding='utf-8') as f_out:
            
            # Use a list comprehension for efficient reading of the .jsonl file
            problems_with_generations = [json.loads(line) for line in f_in]
            
            print(f"Starting evaluation for {len(problems_with_generations)} problems...")
            
            # Loop through each problem using tqdm for a progress bar
            for problem_data in tqdm(problems_with_generations, desc="Evaluating Results"):
                
                solution = problem_data.get('solution')
                generation = problem_data.get('generation')
                
                # Call the core evaluation function from our other script
                scores = evaluate_trace_eng(solution, generation)
                
                # Prepare the final JSON object for output
                # We keep the original metadata and add the scores
                output_entry = {
                    'seed': problem_data.get('seed'),
                    'id': problem_data.get('id'),
                    'level': problem_data.get('level'),
                    'branch': problem_data.get('branch'),
                    'domain': problem_data.get('domain'),
                    'area': problem_data.get('area'),
                    'model': MODEL_NAME,
                    'scores': scores  # Nest the scores dictionary
                }
                
                # Write the scored result to the new file
                f_out.write(json.dumps(output_entry) + '\n')
                
        print(f"\nEvaluation complete. Scored results saved to '{OUTPUT_FILE}'.")
