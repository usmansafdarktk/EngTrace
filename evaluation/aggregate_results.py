import os
import json
import numpy as np
import pickle
import argparse  
from tqdm import tqdm

# Argument Parser Setup
parser = argparse.ArgumentParser(
    description="Summarize model evaluation results from .jsonl files."
)
# Add an argument for the input directory
parser.add_argument(
    "evals_dir",
    type=str,
    help="The source directory containing the evaluation .jsonl files."
)
# Add an argument for the output directory
parser.add_argument(
    "output_dir",
    type=str,
    help="The destination directory to save the summary .pkl and .txt files."
)
args = parser.parse_args()


# Configuration
EVALS_DIR = args.evals_dir
OUTPUT_DIR = args.output_dir
os.makedirs(OUTPUT_DIR, exist_ok=True)


# Create a file to store the LaTeX formatted results
mean_sdv_results_file = open(os.path.join(OUTPUT_DIR, 'mean_sdv_results.txt'), 'w')


# Create dictionaries for your specific categories
branch_wise_results = {}
domain_wise_results = {}
area_wise_results = {}
level_wise_results = {}


#  2. Data Collection 
# Find all evaluation files in the directory
eval_files = [f for f in os.listdir(EVALS_DIR) if f.endswith('_evals_2.jsonl')]

for model_file in eval_files:
    # Initialize dictionaries for the current model
    model_name = model_file.replace('_evals.jsonl', '')
    branch_wise_results[model_name] = {}
    domain_wise_results[model_name] = {}
    area_wise_results[model_name] = {}
    level_wise_results[model_name] = {}
    
    all_metric_vals = []
    
    file_path = os.path.join(EVALS_DIR, model_file)
    with open(file_path, 'r') as f_pred:
        print(f"Processing {model_file}...")
        for line in tqdm(f_pred):
            json_line = json.loads(line)
            
            # Skip if there's an error or missing data
            if not json_line.get('scores'):
                continue
            scores = json_line['scores']
            if any(s is None for s in scores.values()):
                continue

            # Extract scores
            recall, precision, step_f1, final_answer_match, bertscore, rouge2 = \
                scores['recall'], scores['precision'], scores['step_f1'], scores['final_answer_match'], scores['bertscore'], scores['rouge2']
            
            score_list = [final_answer_match, recall, precision, step_f1, bertscore, rouge2]
            all_metric_vals.append(score_list)

            # Collect results grouped by your categories
            branch = json_line['branch']
            domain = json_line['domain']
            area = json_line['area']
            level = json_line['level'].lower()
            
            # Group by branch
            branch_wise_results[model_name].setdefault(branch, []).append(score_list)
            # Group by domain
            domain_wise_results[model_name].setdefault(domain, []).append(score_list)
            # Group by area
            area_wise_results[model_name].setdefault(area, []).append(score_list)
            # Group by level
            level_wise_results[model_name].setdefault(level, []).append(score_list)

    #  3. Calculation and Output Generation 
    
    # Calculate overall mean and standard deviation for the LaTeX output
    mean_vals = np.mean(all_metric_vals, axis=0)
    std_vals = np.std(all_metric_vals, axis=0)
    
    # Write the formatted string to the text file
    header = "Model & Final Answer & Recall & Precision & Step F1 & BERTScore & ROUGE-2 \\\\\n"
    if mean_sdv_results_file.tell() == 0: # Write header only once
        mean_sdv_results_file.write(header)
        
    formatted_scores = " & ".join([f"${m:.3f}_{{{s:.3f}}}$" for m, s in zip(mean_vals, std_vals)])
    mean_sdv_results_file.write(f"{model_name} & {formatted_scores} \\\\\n")

    # Calculate final mean scores for each category
    for category_results in [branch_wise_results, domain_wise_results, area_wise_results, level_wise_results]:
        for model in category_results:
            for category_key, scores_list in category_results[model].items():
                category_results[model][category_key] = np.mean(scores_list, axis=0)

#  4. Save Dictionaries to Pickle Files 
print("\nSaving aggregated results to pickle files...")
pickle.dump(branch_wise_results, open(os.path.join(OUTPUT_DIR, 'branch_wise_results.pkl'), 'wb'))
pickle.dump(domain_wise_results, open(os.path.join(OUTPUT_DIR, 'domain_wise_results.pkl'), 'wb'))
pickle.dump(area_wise_results, open(os.path.join(OUTPUT_DIR, 'area_wise_results.pkl'), 'wb'))
pickle.dump(level_wise_results, open(os.path.join(OUTPUT_DIR, 'level_wise_results.pkl'), 'wb'))

mean_sdv_results_file.close()
print("Analysis complete.")
