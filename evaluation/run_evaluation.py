import os
import json
import argparse 
from tqdm import tqdm
from dotenv import load_dotenv

# Import both evaluation functions but give them unique names (aliases)
from evaluation.bi_encoders_evaluation import evaluate_trace_eng as evaluate_bi_encoder
from evaluation.cross_encoders_llm_judge_evaluation import evaluate_trace_eng as evaluate_cross_llm

#  1. Argument Parsing 
# Set up a parser to accept command-line arguments
parser = argparse.ArgumentParser(description="Run evaluation on model inference results.")
parser.add_argument(
    "evaluator_type", 
    type=str, 
    choices=['bi', 'cross-llm'], 
    help="The type of evaluator to use: 'bi' for Bi-Encoders or 'cross-llm' for Cross-Encoders + LLM Judge."
)
args = parser.parse_args()


#  2. Configuration 
load_dotenv()

# Get the model name from the .env file
MODEL_NAME = os.getenv("MODEL_NAME") 

# Define the input file (this is the same for both)
INPUT_FILE = f"inference_results/{MODEL_NAME}_inference_results.jsonl"

# Dynamically set the output file and evaluation function based on the argument
if args.evaluator_type == 'bi':
    print("Using Bi-Encoder evaluator.")
    OUTPUT_FILE = f"evaluation_results/{MODEL_NAME}_evals.jsonl"
    evaluation_function = evaluate_bi_encoder
elif args.evaluator_type == 'cross-llm':
    print("Using Cross-Encoder + LLM Judge evaluator.")
    OUTPUT_FILE = f"evaluation_results/{MODEL_NAME}_evals_2.jsonl"
    evaluation_function = evaluate_cross_llm


#  3. Main Evaluation Logic 
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
            
            problems_with_generations = [json.loads(line) for line in f_in]
            
            print(f"Starting evaluation for {len(problems_with_generations)} problems...")
            
            # Loop through each problem using tqdm for a progress bar
            for problem_data in tqdm(problems_with_generations, desc=f"Evaluating with {args.evaluator_type}"):
                
                solution = problem_data.get('solution')
                generation = problem_data.get('generation')
                
                if args.evaluator_type == 'bi':
                    # Call the bi-encoder function with two arguments
                    scores = evaluation_function(solution, generation)
                elif args.evaluator_type == 'cross-llm':
                    # For the cross-encoder, first build the context string
                    problem_context = (
                        f"This is a {problem_data.get('level')} level problem in "
                        f"{problem_data.get('branch')} engineering, specifically in the "
                        f"domain of {problem_data.get('domain')} and the area of {problem_data.get('area')}."
                    )
                    # Then call its function with three arguments
                    scores = evaluation_function(solution, generation, problem_context)

                # Prepare the final JSON object for output
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