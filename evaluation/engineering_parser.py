import re
from typing import List, Tuple, Optional


def convert_to_float_eng(text: str) -> Optional[float]:
    """
    Converts a string with a number and optional units to a float.
    Handles integers, floats, and scientific notation (e.g., '1.23e-4').
    
    Args:
        text: The input string to parse.

    Returns:
        The extracted number as a float, or None if no number is found.
    """
    if not isinstance(text, str):
        return None
    
    # This regex is designed to find scientific notation, floats, and integers.
    # It handles patterns like: 123, 123.45, 1.23e-4, -0.98E+2
    # It also ignores surrounding text or units (e.g., "7.65 L").
    match = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', text.replace(',', ''))
    
    if not match:
        return None
    
    try:
        return float(match.group(0))
    except (ValueError, TypeError):
        return None


def extract_final_answer_eng(text: str) -> Optional[float]:
    """
    Extracts the most likely final numerical answer from a text block.
    """
    text = re.sub(r'\s+', ' ', text).strip()

    # Priority 1: Look for the number after the "**Answer:**" tag.
    match = re.search(r'\*\*Answer:\*\*\s*.*?([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', text, re.IGNORECASE)
    if match and match.group(1):
        return convert_to_float_eng(match.group(1))

    # Priority 2: Look for the last number after an equals sign, using word boundaries.
    # The \b ensures we match "V = 100" but not the 0 in "F_A0".
    matches = list(re.finditer(r'=\s*(\b[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\b)', text))
    if matches:
        return convert_to_float_eng(matches[-1].group(1))

    # Priority 3: As a last resort, find the very last standalone number in the string.
    matches = list(re.finditer(r'(\b[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\b)', text))
    if matches:
        return convert_to_float_eng(matches[-1].group(0))
    
    return None


def extract_steps(full_text: str) -> Tuple[List[str], List[Optional[float]], Optional[float]]:
    """
    Parses a full solution text into its constituent parts.
    (IMPROVED VERSION with more robust step prefix cleaning)
    """
    # Regex to find the start of each step, robust to markdown and spacing.
    # It finds both "**Step 1**" and "1." style steps.
    step_starts = [match.start() for match in re.finditer(r'(?:^|\n)\s*(?:\*\*|#*)\s*(?:Step\s*\d+|\d+\.)', full_text, re.IGNORECASE)]
    
    if not step_starts:
        step_texts = [full_text.strip()]
    else:
        step_texts = []
        for i, start_index in enumerate(step_starts):
            end_index = step_starts[i + 1] if i + 1 < len(step_starts) else len(full_text)
            step_texts.append(full_text[start_index:end_index].strip())
    
    # NEW: A more robust regex to clean both "Step X:" and "X." prefixes
    # This is the key change to fix the "1.0" error.
    cleaned_step_texts = [re.sub(r'^(?:\*\*|#*)\s*(?:Step\s*\d+|\d+\.)\s*:\s*', '', text, flags=re.IGNORECASE).strip() for text in step_texts]

    step_answers = [extract_final_answer_eng(text) for text in cleaned_step_texts]
    
    overall_final_answer = extract_final_answer_eng(full_text)

    if overall_final_answer is None and step_answers:
        for answer in reversed(step_answers):
            if answer is not None:
                overall_final_answer = answer
                break
                
    return cleaned_step_texts, step_answers, overall_final_answer


if __name__ == '__main__':
    # --- Test Case 1: Simple CSTR Problem (from before) ---
    sample_solution_cstr = """
    **Step 1:** State the CSTR design equation.
    V = (F_A0 - F_A) / (-r_A)

    **Step 2:** Substitute the given values into the equation.
    V = (1.26 mol/s - 0.99 mol/s) / (0.0353 mol/(L·s))
    V = 0.27 / 0.0353 = 7.6487 L

    **Step 3:** Calculate the final volume.
    V = 7.65 L

    **Answer:** The required reactor volume is 7.65 liters.
    """
    
    print("--- Test Case 1: Simple CSTR Problem ---")
    steps, step_answers, final_answer = extract_steps(sample_solution_cstr)
    print(f"\nOverall Final Answer Extracted: {final_answer}\n")
    for i, (text, answer) in enumerate(zip(steps, step_answers)):
        print(f"--- Step {i+1} ---")
        print(f"Text: {text}")
        print(f"Answer found in step: {answer}")
        print("-" * 15)

    print("\n" + "="*50 + "\n") # Separator for clarity

    