# Cross Encoders + LLM Judge Evaluator
import os
import json
import torch
import numpy as np
from typing import Dict, Any
import google.generativeai as genai
from dotenv import load_dotenv

from engineering_parser import extract_steps
from sentence_transformers import CrossEncoder
from rouge_score import rouge_scorer
from bert_score import BERTScorer


load_dotenv()

# Initialize the client for the judge LLM.
JUDGE_MODEL = None
try:
    api_key = os.getenv("JUDGE_API_KEY")
    if not api_key:
        print("Warning: JUDGE_API_KEY not found in .env file. LLM Judge will be skipped.")
    else:
        genai.configure(api_key=api_key)
        model_name = os.getenv("JUDGE_MODEL_NAME")
        JUDGE_MODEL = genai.GenerativeModel(model_name)
        print(f"LLM Judge initialized with model: {model_name}")
except Exception as e:
    print(f"Warning: Could not initialize Gemini Judge client. Error: {e}")


#  Model & Scorer Initialization 
print("Initializing evaluation models...")
CROSS_ENCODER = CrossEncoder('cross-encoder/stsb-roberta-large')
ROUGE_SCORER = rouge_scorer.RougeScorer(['rouge2', 'rougeL', 'rougeLsum'], use_stemmer=True)
BERT_SCORER = BERTScorer(model_type='allenai/longformer-base-4096', 
                         device='cuda' if torch.cuda.is_available() else 'cpu')
print("Evaluation models initialized.")


#  LLM Judge Helper Function (Gemini) 
def get_error_analysis_from_llm(gt_step: str, pred_step: str, problem_context: str) -> Dict[str, Any]:
    """
    Uses a Gemini model to categorize the error in a predicted reasoning step.
    """
    if not JUDGE_MODEL:
        return {"error_category": "Analysis Skipped", "explanation": "LLM Judge client not initialized."}
    
    prompt = f"""
    You are an expert engineering professor acting as an automated evaluator. Your task is to analyze the 'MODEL'S STEP' against the 'GROUND-TRUTH STEP' and return a structured JSON object based on the rules and inputs below.

    **CRITICAL RULES:**
    1. Your 'explanation' MUST logically justify your chosen 'error_category'. Do not contradict yourself. For example, do not choose "Calculation Error" and then state that the calculation is correct.
    2. If the model's step is factually correct but takes a different path than the ground-truth, you MUST use the "Alternative Correct" category. Do not classify a correct step as "Other".
    3. The final output must be only a raw JSON object. Do not include any introductory text, concluding remarks, or markdown formatting.

    **Error Categories:**
    - "Conceptual Error": The model applied the wrong scientific principle or formula (e.g., used addition instead of subtraction).
    - "Calculation Error": The model used the correct formula but made a mathematical mistake (e.g., 2 * 3 = 5).
    - "Input Error": The model used the correct formula but pulled the wrong number from the problem context or a previous step.
    - "Alternative Correct": The model's step is valid and logically sound, but follows a different method or phrasing than the ground-truth step.
    - "Other": The model's step is nonsensical, irrelevant, a hallucination, or contains only formatting errors.

    **Input for Analysis:**
    [CONTEXT]: {problem_context}
    [GROUND-TRUTH STEP]: {gt_step}
    [MODEL'S STEP]: {pred_step}

    **OUTPUT FORMAT:**
    You must now provide your analysis. Your entire response will be a single, raw JSON object. Adhere strictly to the following format with exactly two keys:
    {{"error_category": "...", "explanation": "..."}}
    """
    
    try:
        generation_config = genai.types.GenerationConfig(
            response_mime_type="application/json",
            temperature=0.0
        )
        response = JUDGE_MODEL.generate_content(prompt, generation_config=generation_config)
        cleaned_json = response.text.strip().replace("```json", "").replace("```", "").strip()
        analysis = json.loads(cleaned_json)
        
        # Validation and Standardization Step 
        # Ensure the output has the correct key, even if the model makes a mistake.
        if "error_classification" in analysis:
            analysis["error_category"] = analysis.pop("error_classification")
            
        # Ensure the dictionary has the required keys, providing defaults if missing.
        if "error_category" not in analysis:
            analysis["error_category"] = "Formatting Error"
        if "explanation" not in analysis:
            analysis["explanation"] = "Model failed to provide an explanation."
            
        return analysis
        
    except Exception as e:
        return {"error_category": "Analysis Failed", "explanation": str(e)}


def safe_bert_score(gt: str, pred: str) -> float:
    """ A wrapper for BERTScore to handle potential errors and empty strings. """
    if not all(isinstance(s, str) for s in [gt, pred]) or not gt.strip() or not pred.strip():
        return 0.0
    try:
        _, _, f1 = BERT_SCORER.score([pred], [gt])
        return f1.item()
    except Exception as e:
        print(f"Warning: BERTScore failed with error: {e}")
        return 0.0


def evaluate_trace_eng(gt_solution: str, pred_generation: str, problem_context: str) -> Dict[str, Any]:
    """
    Compares a ground-truth engineering solution with a model's generation.
    """

    # Handle empty inputs gracefully
    if not gt_solution or not pred_generation:
        return {
            'error': 'Input solution or generation is empty.', 
            'recall': 0, 
            'precision': 0, 
            'step_f1': 0, 
            'final_answer_match': 0, 
            'rouge2': 0, 
            'rougeL': 0, 
            'rougeLsum': 0, 
            'bertscore': 0, 
            'error_analysis': []
        }
    
    # Compute textual similarity metrics
    rouge_scores = ROUGE_SCORER.score(gt_solution, pred_generation)
    bertscore = safe_bert_score(gt_solution, pred_generation)
    error_analyses = []

    # Extract structured reasoning steps and final answers
    gt_steps, gt_step_answers, gt_final_answer = extract_steps(gt_solution)
    pred_steps, pred_step_answers, pred_final_answer = extract_steps(pred_generation)
    final_answer_match = 0
    FINAL_ANSWER_TOLERANCE = 0.01

    # Final answer comparison with tolerance
    if gt_final_answer is not None and pred_final_answer is not None:
        if abs(gt_final_answer) > 1e-9:
            if abs(gt_final_answer - pred_final_answer) / abs(gt_final_answer) < FINAL_ANSWER_TOLERANCE:
                final_answer_match = 1
        elif abs(gt_final_answer - pred_final_answer) < 1e-9:
            final_answer_match = 1

    # Step-level similarity and recall/precision computation
    recall, precision = 0, 0
    if gt_steps and pred_steps:
        sentence_pairs = [[gt_step, pred_step] for gt_step in gt_steps for pred_step in pred_steps]
        scores = CROSS_ENCODER.predict(sentence_pairs, show_progress_bar=False)
        semantic_similarity = np.array(scores).reshape(len(gt_steps), len(pred_steps))

        # Numeric correctness check
        numeric_correctness = np.zeros((len(gt_steps), len(pred_steps)))
        STEP_ANSWER_TOLERANCE = 0.02

        for i, gt_ans in enumerate(gt_step_answers):
            if gt_ans is None: continue
            for j, pred_ans in enumerate(pred_step_answers):
                if pred_ans is None: continue
                if abs(gt_ans) > 1e-9:
                    if (abs(gt_ans - pred_ans) / abs(gt_ans)) < STEP_ANSWER_TOLERANCE:
                        numeric_correctness[i, j] = 1
                elif abs(gt_ans - pred_ans) < 1e-9:
                    numeric_correctness[i, j] = 1
        
        # Combine semantic and numeric correctness matrices
        combined_matrix = np.multiply(semantic_similarity, numeric_correctness)
        SIMILARITY_THRESHOLD = 0.7
        best_matches_scores = np.max(combined_matrix, axis=1)
        
        # Identify and analyze mismatched steps
        for i, score in enumerate(best_matches_scores):
            if score < SIMILARITY_THRESHOLD:
                gt_mismatched_step = gt_steps[i]
                best_pred_index = np.argmax(combined_matrix[i, :])
                pred_mismatched_step = pred_steps[best_pred_index]
                analysis = get_error_analysis_from_llm(gt_mismatched_step, pred_mismatched_step, problem_context)
                error_analyses.append({
                    "mismatched_gt_step_index": i,
                    "ground_truth_step": gt_mismatched_step,
                    "closest_predicted_step": pred_mismatched_step,
                    "analysis": analysis
                })
        
        # Compute recall and precision based on matched steps
        recall = float(np.sum(best_matches_scores > SIMILARITY_THRESHOLD) / len(gt_steps))
        precision = float(np.sum(np.max(combined_matrix, axis=0) > SIMILARITY_THRESHOLD) / len(pred_steps))

    # Compute step-level F1 score
    step_f1 = 0
    if recall + precision > 0:
        step_f1 = 2 * (recall * precision) / (recall + precision)
    
    # Final result dictionary
    return {
        'recall': recall, 
        'precision': precision, 
        'step_f1': step_f1,
        'final_answer_match': final_answer_match, 
        'rouge2': rouge_scores['rouge2'].fmeasure,
        'rougeL': rouge_scores['rougeL'].fmeasure, 
        'rougeLsum': rouge_scores['rougeLsum'].fmeasure,
        'bertscore': bertscore, 'error_analysis': error_analyses
    }
