# evaluation_methods.py
import math
from collections import Counter
import torch
import numpy as np
from typing import Dict, Any
from engineering_parser import extract_steps
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
from rouge_score import rouge_scorer
from bert_score import BERTScorer


# Model & Scorer Initialization
# These models are loaded once when the module is imported, which is efficient.
print("Initializing evaluation models...")
SENTENCE_MODEL = SentenceTransformer('sentence-transformers/all-MiniLM-L6-v2')
ROUGE_SCORER = rouge_scorer.RougeScorer(['rouge2', 'rougeL', 'rougeLsum'], use_stemmer=True)
BERT_SCORER = BERTScorer(model_type='allenai/longformer-base-4096', 
                         device='cuda' if torch.cuda.is_available() else 'cpu')
print("Evaluation models initialized.")


# Cos similarity on two bags of words over given threshold
def cos_similarity(a: str, b: str, min_val: float = 0.0, max_val: float = 1.0) -> float:
    """Pure cosine similarity between two raw strings (bag-of-words)."""
    c1, c2 = Counter(a.split()), Counter(b.split())
    vocab = set(c1) | set(c2)
    dot   = sum(c1[w] * c2[w] for w in vocab)
    n1    = math.sqrt(sum(c1[w] ** 2 for w in vocab))
    n2    = math.sqrt(sum(c2[w] ** 2 for w in vocab))
    val   = dot / (n1 * n2) if n1 and n2 else 0.0
    if val >= min_val and val <= max_val:
        return val
    return None


def safe_bert_score(gt: str, pred: str) -> float:
    """ A wrapper for BERTScore to handle potential errors and empty strings. """
    if not all(isinstance(s, str) for s in [gt, pred]) or not gt.strip() or not pred.strip():
        return 0.0
    try:
        # BERTScore returns Precision, Recall, and F1. We use F1.
        _, _, f1 = BERT_SCORER.score([pred], [gt])
        return f1.item()
    except Exception as e:
        print(f"Warning: BERTScore failed with error: {e}")
        return 0.0

def evaluate_trace_eng(gt_solution: str, pred_generation: str) -> Dict[str, Any]:
    """
    Compares a ground-truth engineering solution with a model's generation.

    This function calculates final answer accuracy, reasoning step recall/precision
    (ChainEval metric), and standard text similarity scores.

    Args:
        gt_solution: The ground-truth solution text.
        pred_generation: The model-generated solution text.

    Returns:
        A dictionary containing all calculated evaluation metrics.
    """
    if not gt_solution or not pred_generation:
        return {
            'error': 'Input solution or generation is empty.',
            'recall': 0, 'precision': 0, 'step_f1': 0, 'final_answer_match': 0,
            'rouge2': 0, 'rougeL': 0, 'rougeLsum': 0, 'bertscore': 0
        }
    
    #  1. Standard Text Similarity Metrics 
    rouge_scores = ROUGE_SCORER.score(gt_solution, pred_generation)
    bertscore = safe_bert_score(gt_solution, pred_generation)
    
    #  2. Parse Both Traces using the Engineering Parser 
    gt_steps, gt_step_answers, gt_final_answer = extract_steps(gt_solution)
    pred_steps, pred_step_answers, pred_final_answer = extract_steps(pred_generation)

    #  3. Final Answer Match 
    final_answer_match = 0
    # Use a 1% relative tolerance for final answers, common in engineering.
    FINAL_ANSWER_TOLERANCE = 0.01
    if gt_final_answer is not None and pred_final_answer is not None:
        if abs(gt_final_answer) > 1e-9:  # Avoid division by zero for non-zero answers
            if abs(gt_final_answer - pred_final_answer) / abs(gt_final_answer) < FINAL_ANSWER_TOLERANCE:
                final_answer_match = 1
        else:  # Handle cases where the answer is zero
            if abs(gt_final_answer - pred_final_answer) < 1e-9:
                final_answer_match = 1

    #  4. Reasoning Step Evaluation (ChainEval Metric) 
    if not gt_steps or not pred_steps:
        # Cannot calculate recall/precision if either trace has no steps
        recall, precision = 0, 0
    else:
        # Embed all step texts for semantic comparison
        gt_embeddings = SENTENCE_MODEL.encode(gt_steps)
        pred_embeddings = SENTENCE_MODEL.encode(pred_steps)

        # Create a numeric correctness matrix: 1 if numbers match, 0 otherwise
        numeric_correctness = np.zeros((len(gt_steps), len(pred_steps)))
        STEP_ANSWER_TOLERANCE = 0.02  # Use a slightly higher 2% tolerance for intermediate steps
        for i, gt_ans in enumerate(gt_step_answers):
            if gt_ans is None: continue
            for j, pred_ans in enumerate(pred_step_answers):
                if pred_ans is None: continue
                
                if abs(gt_ans) > 1e-9:
                    if (abs(gt_ans - pred_ans) / abs(gt_ans)) < STEP_ANSWER_TOLERANCE:
                        numeric_correctness[i, j] = 1
                elif abs(gt_ans - pred_ans) < 1e-9:
                    numeric_correctness[i, j] = 1
        
        # Calculate semantic similarity matrix
        semantic_similarity = cosine_similarity(gt_embeddings, pred_embeddings)
        
        # Combine matrices: similarity is only valid if numbers are also correct
        combined_matrix = np.multiply(semantic_similarity, numeric_correctness)
        
        # Calculate recall and precision
        SIMILARITY_THRESHOLD = 0.7  # Threshold for considering a step as "matched"
        recall = float(np.sum(np.max(combined_matrix, axis=1) > SIMILARITY_THRESHOLD) / len(gt_steps))
        precision = float(np.sum(np.max(combined_matrix, axis=0) > SIMILARITY_THRESHOLD) / len(pred_steps))

    #  5. Compile and Return All Metrics 
    step_f1 = 0
    if recall + precision > 0:
        step_f1 = 2 * (recall * precision) / (recall + precision)

    return {
        'recall': recall,
        'precision': precision,
        'step_f1': step_f1,
        'final_answer_match': final_answer_match,
        'rouge2': rouge_scores['rouge2'].fmeasure,
        'rougeL': rouge_scores['rougeL'].fmeasure,
        'rougeLsum': rouge_scores['rougeLsum'].fmeasure,
        'bertscore': bertscore
    }