from engineering_parser import extract_steps
from transformers import AutoTokenizer
import json
import re
import random
import numpy as np
import openai
import anthropic
import google.generativeai as genai
from typing import List, Dict, Any
from collections import Counter
from scipy.optimize import linear_sum_assignment


# Cross Encoders
import torch
import numpy as np
from typing import Dict, Any
from sentence_transformers import CrossEncoder
from rouge_score import rouge_scorer
from bert_score import BERTScorer


# Retry and Timeout constants
MAX_RETRIES = 5
INITIAL_BACKOFF = 2  # seconds
REQUEST_TIMEOUT = 120 # seconds (2 minutes)


#  Model & Scorer Initialization 
print("Initializing evaluation models...")
device = 'cuda' if torch.cuda.is_available() else 'cpu'
CROSS_ENCODER = CrossEncoder('cross-encoder/stsb-roberta-large', device=device)
ROUGE_SCORER = rouge_scorer.RougeScorer(['rouge2', 'rougeL', 'rougeLsum'], use_stemmer=True)
print("Evaluation models initialized.")

# Initialize BERT_SCORER
BERT_SCORER = BERTScorer(model_type='allenai/longformer-base-4096', device=device) 

# Add the tokenizer initialization
print("Initializing Longformer tokenizer...")
TOKENIZER = AutoTokenizer.from_pretrained('allenai/longformer-base-4096')
MAX_LEN = 4096 # Model's max sequence length
print("Tokenizer initialized.")


def safe_bert_score(gt: str, pred: str) -> float:
    """
    A wrapper for BERTScore that handles errors, empty strings,
    and only truncates inputs if they exceed the model's max length (4096 tokens).
    """
    # Initial check for valid string inputs
    if not all(isinstance(s, str) for s in [gt, pred]) or not gt.strip() or not pred.strip():
        # print("Warning: Empty input(s) provided to BERTScore.") # Optional debug log
        return 0.0

    gt_to_score = gt
    pred_to_score = pred
    needs_truncation = False

    try:
        # 1. Tokenize without truncation to check lengths
        gt_tokens = TOKENIZER.encode(gt, return_tensors='pt')
        pred_tokens = TOKENIZER.encode(pred, return_tensors='pt')

        gt_len = gt_tokens.shape[1]
        pred_len = pred_tokens.shape[1]

        # 2. Check if truncation is needed for either input
        if gt_len > MAX_LEN or pred_len > MAX_LEN:
            needs_truncation = True
            log_msgs = []
            if gt_len > MAX_LEN:
                log_msgs.append(f"GT (len={gt_len})")
                # Re-tokenize and decode GT with truncation
                gt_tokens_truncated = TOKENIZER.encode(gt, max_length=MAX_LEN, truncation=True, return_tensors='pt')
                gt_to_score = TOKENIZER.decode(gt_tokens_truncated[0], skip_special_tokens=True)
            if pred_len > MAX_LEN:
                log_msgs.append(f"Pred (len={pred_len})")
                # Re-tokenize and decode Pred with truncation
                pred_tokens_truncated = TOKENIZER.encode(pred, max_length=MAX_LEN, truncation=True, return_tensors='pt')
                pred_to_score = TOKENIZER.decode(pred_tokens_truncated[0], skip_special_tokens=True)
            
            # Optional: Log which inputs were truncated
            # print(f"Warning: Truncating input(s) for BERTScore exceeding {MAX_LEN} tokens: {', '.join(log_msgs)}")

    except Exception as e:
         # Handle potential errors during tokenization itself
         print(f"Warning: Tokenization/Length Check failed: {e}. Skipping BERTScore.")
         return 0.0


    # 3. Calculate BERTScore using potentially truncated inputs
    try:
        # Ensure inputs are not empty after potential truncation/decoding
        if not gt_to_score.strip() or not pred_to_score.strip():
             # print("Warning: Input(s) became empty after potential truncation for BERTScore.") # Optional debug log
             return 0.0

        _, _, f1 = BERT_SCORER.score([pred_to_score], [gt_to_score])

        # Check for NaN or Inf F1 scores
        if torch.isnan(f1) or torch.isinf(f1):
            print("Warning: BERTScore returned NaN or Inf. Returning 0.0")
            return 0.0

        return f1.item()

    except Exception as e:
        # Catch calculation errors (including CUDA errors)
        print(f"Warning: BERTScore calculation failed with error: {e}")
        if "CUDA error: device-side assert triggered" in str(e):
            print(f"  >>> Likely CUDA assert error. Truncation applied: {needs_truncation}. Input lengths (tokens): GT={gt_len if 'gt_len' in locals() else 'N/A'}, Pred={pred_len if 'pred_len' in locals() else 'N/A'}")
        return 0.0


# CONFIGURATION & CONSTANTS
FINAL_TOLERANCE = 0.02    # 2%
STEP_TOLERANCE = 0.02     # 2%
SEMANTIC_TAU = 0.70       # Threshold for Tier 1
ERROR_SAMPLE_RATE = 0.20  # 20% sampling for wrong answers
TRIBUNAL_TRIGGER_THRESHOLD = 0.80 # Trigger if auto-match < 80% on correct answers

# Model Identifiers
MODEL_OPENAI = "gpt-5"                 
MODEL_ANTHROPIC = "claude-opus-4-5-20251101" 
MODEL_GOOGLE = "gemini-3-pro-preview"  


class EngTraceFramework:
    def __init__(self, api_keys: Dict[str, str]):
        """
        Initializes the clients for the AI Tribunal.
        Assumes CROSS_ENCODER, ROUGE_SCORER, BERT_SCORER are globally available 
        from your previous cells.
        """
        print("Initializing EngTrace Tribunal Clients...")
        
        # 1. OpenAI Client
        self.client_openai = None
        if api_keys.get('openai'):
            self.client_openai = openai.OpenAI(api_key=api_keys['openai'])

        # 2. Anthropic Client
        self.client_anthropic = None
        if api_keys.get('anthropic'):
            self.client_anthropic = anthropic.Anthropic(api_key=api_keys['anthropic'])

        # 3. Google Client
        if api_keys.get('google'):
            genai.configure(api_key=api_keys['google'])
            
        print("EngTrace Framework Ready.")

    # --------------------------------------------------------------------------
    # TIER 1: AUTOMATED STRICT VERIFICATION
    # --------------------------------------------------------------------------
    def _tier1_verify_matrix(self, gt_s_txt, gt_s_val, pred_s_txt, pred_s_val):
        M, N = len(gt_s_txt), len(pred_s_txt)
        if M == 0 or N == 0: return np.zeros((M, N))
    
        # 1. BATCH SEMANTIC PREDICTION (Massive Speedup)
        # Create all pairs at once instead of looping
        pairs = [(gt_s_txt[i], pred_s_txt[j]) for i in range(M) for j in range(N)]
        sem_scores_flat = CROSS_ENCODER.predict(pairs)
        sem_matrix = sem_scores_flat.reshape(M, N)
        
        V = np.zeros((M, N))
    
        for i in range(M):
            for j in range(N):
                gt_v, pred_v = gt_s_val[i], pred_s_val[j]
                
                # 2. STRICT NUMERIC LOGIC
                num_match = False
                if gt_v is None and pred_v is None:
                    num_match = True # Both are text-only -> Pass Numeric
                elif gt_v is not None and pred_v is not None:
                    if abs(gt_v) < 1e-9: num_match = abs(pred_v) < 1e-9
                    else: num_match = abs((gt_v - pred_v)/gt_v) <= STEP_TOLERANCE
                
                # 3. STRICT AND LOGIC (Must be Num Match AND Semantic Match)
                if num_match and (sem_matrix[i, j] >= SEMANTIC_TAU):
                    V[i, j] = 1.0
                    
        return V

    # --------------------------------------------------------------------------
    # TIER 2: AI TRIBUNAL & VOTING LOGIC
    # --------------------------------------------------------------------------
    def _call_single_judge(self, provider: str, prompt: str) -> List[Dict]:
        """
        Handles API calls with "Anchor-Based Extraction" to robustly ignore 
        Math/LaTeX/Chain-of-Thought and extract only the valid JSON.
        """
        results = []
        raw_content = ""
        
        try:
            # 1. PROVIDER-SPECIFIC API CALLS
            if provider == 'openai' and self.client_openai:
                completion_args = {
                    "model": MODEL_OPENAI,
                    "messages": [{"role": "user", "content": prompt}]
                }
                # Handle o1/Reasoning models (No temp, no explicit JSON mode)
                if "o1" not in MODEL_OPENAI and "gpt-5" not in MODEL_OPENAI:
                     completion_args["response_format"] = {"type": "json_object"}
                     completion_args["temperature"] = 0.0

                resp = self.client_openai.chat.completions.create(**completion_args)
                raw_content = resp.choices[0].message.content

            elif provider == 'anthropic' and self.client_anthropic:
                resp = self.client_anthropic.messages.create(
                    model=MODEL_ANTHROPIC, max_tokens=2048, # Increased for verbose chain-of-thought
                    messages=[{"role": "user", "content": prompt}]
                )
                raw_content = resp.content[0].text

            elif provider == 'google':
                model = genai.GenerativeModel(MODEL_GOOGLE)
                # Google usually respects this, but we parse robustly anyway
                resp = model.generate_content(prompt, generation_config={"response_mime_type": "application/json"})
                raw_content = resp.text

            # 2. UNIFIED SMART PARSING LOGIC 
            # This logic works for ALL providers. It handles:
            # - Preambles ("Here is the analysis...")
            # - LaTeX Math ({CSTR} = ...)
            # - Markdown Code Blocks (```json ... ```)
            # - Single Quotes vs Double Quotes
            
            data = None
            
            # A. Fast Path: Clean Markdown and try direct parse
            clean_content = re.sub(r'^```json\s*', '', raw_content, flags=re.MULTILINE)
            clean_content = re.sub(r'\s*```$', '', clean_content, flags=re.MULTILINE)
            try:
                data = json.loads(clean_content)
            except:
                # B. Anchor Path: Look specifically for the "results" key
                # This ignores any "{}" that appear in math/text before the actual JSON.
                
                # Normalize text to handle potential single-quote keys from weaker models
                # (We don't replace all ' to " because that breaks text, just search smartly)
                target_key = "results"
                key_idx = clean_content.find(f'"{target_key}"')
                if key_idx == -1: key_idx = clean_content.find(f"'{target_key}'")
                
                if key_idx != -1:
                    # Find the opening brace '{' strictly BEFORE the "results" key
                    start_idx = clean_content.rfind('{', 0, key_idx)
                    
                    if start_idx != -1:
                        # Find the closing brace '}' strictly AFTER the key
                        # We try to find the last valid matching brace, or just the end
                        snippet = clean_content[start_idx:]
                        
                        # Greedy cleanup: Try to parse increasingly shorter substrings from the end
                        # This handles cases where the model writes "Conclusion: ..." after the JSON
                        # We look for the last '}' and try to parse.
                        
                        last_brace = snippet.rfind('}')
                        while last_brace != -1:
                            sub_snippet = snippet[:last_brace+1]
                            try:
                                data = json.loads(sub_snippet)
                                break # Success!
                            except:
                                try:
                                    import ast
                                    data = ast.literal_eval(sub_snippet)
                                    break # Success!
                                except:
                                    # Move to the previous '}' and try again
                                    last_brace = snippet.rfind('}', 0, last_brace)

            if data is None:
                raise ValueError("Could not extract valid JSON object containing 'results'")

            results = data.get('results', data.get('steps', []))

        except Exception as e:
            # Robust Logging
            safe_preview = str(raw_content)
            if len(safe_preview) > 500: safe_preview = safe_preview[:500] + "... [TRUNCATED]"
            
            print(f"\n[!] Judge {provider} FAILED.")
            print(f"    Error: {e}")
            print(f"    Raw Content Snippet:\n{safe_preview}\n")
            return []

        if isinstance(results, dict): results = [results]
        return results

    def _tier2_tribunal_batch(self, question, gt_steps, pred_steps, mismatch_indices):
        """
        1. Constructs the Prompt.
        2. Calls all 3 Judges (OpenAI, Anthropic, Google).
        3. Aggregates votes using Majority Rule.
        """
        # Construct the Prompt
        prompt = f"""
        You are an expert engineering professor acting as an automated evaluator. 
        Your task is to analyze the Student's reasoning steps against the Ground Truth steps.

        QUESTION: 
        {question}

        GROUND TRUTH CHAIN:
        {json.dumps(gt_steps, indent=1)}

        STUDENT CHAIN:
        {json.dumps(pred_steps, indent=1)}

        TASK: 
        Analyze the Student Steps at these specific indices (0-based): {mismatch_indices}.
        Compare them logically to the Ground Truth.

        CRITICAL RULES:
        1. If the step is factually correct but takes a different path than the Ground Truth, you MUST use "Alternative Correct".
        2. Do not contradict the category definitions below.
        3. The final output must be a raw JSON object only.

        ERROR CATEGORIES:
        - "Conceptual Error": The model applied the wrong scientific principle or formula.
        - "Calculation Error": The model used the correct formula but made a mathematical mistake.
        - "Alternative Correct": The model's step is valid and logically sound, but follows a different method/phrasing.
        - "Other": Nonsensical, irrelevant, hallucination, or formatting errors.

        OUTPUT FORMAT:
        Return a single JSON object containing a "results" list. Each entry must have exactly two keys: "step_index" and "category".
        {{
            "results": [
                {{"step_index": {mismatch_indices[0] if mismatch_indices else 0}, "category": "Alternative Correct"}}
            ]
        }}
        """

        # Collection structure: votes[step_index] = [cat_vote_1, cat_vote_2, cat_vote_3]
        votes = {idx: [] for idx in mismatch_indices}
        
        # Call Providers
        providers = []
        if self.client_openai: providers.append('openai')
        if self.client_anthropic: providers.append('anthropic')
        # Check if google key exists in simple check
        try:
             genai.get_model_info(MODEL_GOOGLE)
             providers.append('google')
        except: pass

        for provider in providers:
            judge_output = self._call_single_judge(provider, prompt)
            for item in judge_output:
                s_idx = item.get('step_index')
                cat = item.get('category')
                if s_idx in votes:
                    votes[s_idx].append(cat)

        # Majority Voting & Score Calculation
        final_scores = {}
        error_counts = Counter()

        for idx, cat_list in votes.items():
            if not cat_list: continue

            # Convert categories to scalar scores
            scalar_votes = []
            for cat in cat_list:
                c_lower = str(cat).lower()
                score = 0.0
                if "alternative" in c_lower: score = 1.0
                elif "calculation" in c_lower: score = 0.5
                else: score = 0.0 # Conceptual / Other
                
                scalar_votes.append(score)
                error_counts[cat] += 1 # For meta-analysis

            # MAJORITY VOTE FORMULA (Checks for >50% Majority)
            if not scalar_votes:
                consensus_score = 0.0
            else:
                counts = Counter(scalar_votes)
                total_votes = len(scalar_votes)
                majority_found = False
                
                # Check if any score has > 50% of the votes
                for sc, count in counts.most_common():
                    if count > (total_votes / 2):
                        consensus_score = sc
                        majority_found = True
                        break
                
                # If no majority (e.g. 1.0, 0.5, 0.0), use Conservative Tie-Breaker
                if not majority_found:
                    consensus_score = min(scalar_votes)
            
            final_scores[idx] = consensus_score

        # 1. Determine dominant error (for quick summary)
        dominant_error = error_counts.most_common(1)[0][0] if error_counts else "None"
        
        # 2. Convert Counter to standard dict for JSON serialization
        error_breakdown = dict(error_counts)
        
        # Return breakdown alongside scores and dominant error
        return final_scores, dominant_error, error_breakdown

    # --------------------------------------------------------------------------
    # MAIN EVALUATION LOGIC
    # --------------------------------------------------------------------------
    def evaluate_entry(self, entry: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single JSONL entry.
        """
        # 1. Extraction (Using your existing functions)
        gt_txt = entry.get('solution', '')
        pred_txt = entry.get('generation', '') 
        question = entry.get('question', '')

        gt_s_txt, gt_s_val, gt_final = extract_steps(gt_txt)
        pred_s_txt, pred_s_val, pred_final = extract_steps(pred_txt)

        # 2. Final Answer Accuracy
        # "Response is correct only if relative error < 2%" 
        is_final_correct = False
        if gt_final is not None and pred_final is not None:
            if abs(gt_final) < 1e-9:
                is_final_correct = abs(pred_final) < 1e-9
            else:
                is_final_correct = abs((gt_final - pred_final) / gt_final) < FINAL_TOLERANCE
        
        final_acc_score = 1.0 if is_final_correct else 0.0

        
        # 3. Tier 1 Verification (Batched)
        # Calculates M x N matrix V
        V = self._tier1_verify_matrix(gt_s_txt, gt_s_val, pred_s_txt, pred_s_val)
        M, N = V.shape
        

        # 4. Tribunal Trigger Logic
        # Check current alignment quality
        row_ind, col_ind = linear_sum_assignment(V, maximize=True)
        tier1_matches = V[row_ind, col_ind].sum()
        match_ratio = tier1_matches / max(M, 1)

        trigger_tribunal = False
        
        if is_final_correct:
            # Case A: Right for Wrong Reasons check
            # If answer is correct but steps don't match (<80%), assume "Alternative Correct" and verify
            if match_ratio < TRIBUNAL_TRIGGER_THRESHOLD:
                trigger_tribunal = True
        else:
            # Case B: Wrong Answer analysis
            # Sample 20% of failures to populate Error Taxonomy
            if random.random() < ERROR_SAMPLE_RATE:
                trigger_tribunal = True

        # 5. Execute Tier 2 (If triggered)
        dominant_error = "N/A"
        error_breakdown = {}
        
        if trigger_tribunal and N > 0:
            # Identify columns (pred steps) that are NOT strictly matched
            # We only send mismatched steps to save tokens/cost
            matched_pred_cols = set(col_ind[V[row_ind, col_ind] == 1.0])
            mismatch_indices = [j for j in range(N) if j not in matched_pred_cols]

            if mismatch_indices:
                # CALL THE JUDGES
                judge_scores, dominant_error, error_breakdown = self._tier2_tribunal_batch(
                    question, gt_s_txt, pred_s_txt, mismatch_indices
                )

                # Update Matrix V with Recovered Scores
                for j_idx, score in judge_scores.items():
                    if score > 0:
                        # We need to assign this score to the best aligned GT step.
                        # Heuristic: Find the GT step with highest Semantic Similarity to this Pred step
                        # This aligns the "Alternative Correct" logic to the corresponding standard logic.
                        sem_sims = [CROSS_ENCODER.predict([(gt_s_txt[i], pred_s_txt[j_idx])])[0] for i in range(M)]
                        best_gt_row = np.argmax(sem_sims)
                        
                        # Only update if the judge score is better than what we had
                        V[best_gt_row, j_idx] = max(V[best_gt_row, j_idx], score)

        # 6. Final Metrics Calculation (Post-Tribunal)
        # Re-run Hungarian Matching on the updated V matrix
        row_ind, col_ind = linear_sum_assignment(V, maximize=True)
        R_total = V[row_ind, col_ind].sum()

        precision = R_total / N if N > 0 else 0.0
        recall = R_total / M if M > 0 else 0.0
        
        f1_rec = 0.0
        if (precision + recall) > 0:
            f1_rec = 2 * (precision * recall) / (precision + recall)

        # 7. Textual Metrics (ROUGE / BERTScore)
        # Using the global objects
        r2 = ROUGE_SCORER.score(gt_txt, pred_txt)['rouge2'].fmeasure
        rL = ROUGE_SCORER.score(gt_txt, pred_txt)['rougeL'].fmeasure
        rLsum = ROUGE_SCORER.score(gt_txt, pred_txt)['rougeLsum'].fmeasure
        bs = safe_bert_score(gt_txt, pred_txt)

        # 8. Formatting the Output Entry
        # We assume 'entry' has the metadata (id, seed, etc.)
        output_entry = entry.copy()
        
        # Update/Overwrite the 'scores' dict
        output_entry["scores"] = {
            "recall": float(recall),
            "precision": float(precision),
            "recovered_f1": float(f1_rec),          # The main metric
            "final_answer_acc": float(final_acc_score),
            "rouge2": float(r2),
            "rougeL": float(rL),
            "rougeLsum": float(rLsum),
            "bertscore": float(bs)
        }
        
        # Add Tribunal Metadata
        output_entry["meta"] = {
            "tribunal_triggered": trigger_tribunal,
            "dominant_error_category": dominant_error,
            "error_breakdown": error_breakdown
        }

        return output_entry
