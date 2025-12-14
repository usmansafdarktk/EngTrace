import os
import json
import asyncio
import traceback
from pathlib import Path
from dotenv import load_dotenv

# API Clients
from openai import AsyncOpenAI
from anthropic import AsyncAnthropic
import google.generativeai as genai

# Import your existing loader logic
from ai_assisted_quality_assurance.template_loader import (
    get_branches, 
    get_areas, 
    get_template_files, 
    load_template_functions,
    get_source_code
)

#  1. Setup & Config 
load_dotenv()

# Model Configs
OPENAI_MODEL = os.getenv("OPENAI_MODEL_NAME")
ANTHROPIC_MODEL = os.getenv("ANTHROPIC_MODEL_NAME")
GOOGLE_MODEL = os.getenv("GOOGLE_MODEL_NAME")

# Initialize Clients
openai_client = AsyncOpenAI(api_key=os.getenv("OPENAI_API_KEY"))
anthropic_client = AsyncAnthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))
genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))

# Output Directory
RESULTS_DIR = Path("ai_assistance_data/qa_results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Files definitions
FILE_OPENAI = RESULTS_DIR / "results_openai.jsonl"
FILE_ANTHROPIC = RESULTS_DIR / "results_anthropic.jsonl"
FILE_GOOGLE = RESULTS_DIR / "results_google.jsonl"

#  2. Prompt Template 
FULL_QA_PROMPT = """You are an expert engineering professor and a senior Python developer acting as a peer reviewer for the EngChain benchmark. 

Your task is to meticulously evaluate a new problem template based on its source code and several example outputs. Analyze the provided information and then respond ONLY with a single, valid JSON object that strictly adheres to the schema described below. Do not add any explanatory text or markdown formatting around the JSON object.

Please evaluate the following engineering problem template.

**1. Template Source Code:**
python {template_code}

**2. Generated Instances from the Template:**
Instance 1:
- Question: "{q1}"
- Solution: "{s1}"

Instance 2:
- Question: "{q2}"
- Solution: "{s2}"

Instance 3:
- Question: "{q3}"
- Solution: "{s3}"

**3. Evaluation Rubric & JSON Schema:**
Evaluate the template based on the rubric below. The `human_review_flag` should be `true` if any score is less than 4. The `explanation` should be a concise, one-sentence justification for the scores and the flag.

{{
    "physical_plausibility_score": <integer, a score from 1--5 based on whether the problem respects the laws of physics and engineering>,
    "mathematical_correctness_score": <integer, a score from 1--5 based on whether the equations and calculations are accurate>,
    "pedagogical_clarity_score": <integer, a score from 1--5 based on whether the problem statement is clear, unambiguous, and solvable>,
    "confidence_score": <integer, a score from 1--5 indicating your confidence in this evaluation>,
    "human_review_flag": <boolean, true if the template requires human inspection, otherwise false>,
    "explanation": "<string, a concise, one-sentence justification for the scores and flag>"
}}"""


#  3. API Call Functions 

async def call_openai(prompt):
    # Check if we should run this based on file existence
    if FILE_OPENAI.exists():
        # If file exists but is empty (0 bytes), treating it as needing run is safer, 
        # but per your request to "skip if present", we return None here to signal skip.
        # However, since we are appending line-by-line, simple existence check 
        # might skip the whole batch if you had partial results. 
        # For safety, let's keep the function capability here and control execution in main()
        pass

    try:
        response = await openai_client.chat.completions.create(
            model=OPENAI_MODEL,
            messages=[{"role": "user", "content": prompt}],
            response_format={"type": "json_object"}
            # Removed temperature=0 to fix the o1-preview error
        )
        return json.loads(response.choices[0].message.content)
    except Exception as e:
        return {"error": str(e), "model": "openai"}

async def call_anthropic(prompt):
    try:
        response = await anthropic_client.messages.create(
            model=ANTHROPIC_MODEL,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=1000,
            temperature=0
        )
        return json.loads(response.content[0].text)
    except Exception as e:
        return {"error": str(e), "model": "anthropic"}

async def call_google(prompt):
    try:
        model = genai.GenerativeModel(GOOGLE_MODEL)
        response = await model.generate_content_async(
            prompt,
            generation_config={"response_mime_type": "application/json"}
        )
        return json.loads(response.text)
    except Exception as e:
        return {"error": str(e), "model": "google"}


#  4. Main Processing Logic 

async def process_template(branch, area, file, func_name, func_obj, run_flags):
    """
    run_flags: dict {'openai': bool, 'anthropic': bool, 'google': bool}
    If False, we skip that model entirely.
    """
    print(f"Processing: {func_name}...")
    
    # A. Get Source Code & Instances (Only needed if at least one model runs)
    if not any(run_flags.values()):
        return None

    source_code = get_source_code(branch, area, file)
    
    # B. Generate 3 Instances locally
    instances = []
    try:
        for _ in range(3):
            q, s = func_obj()
            instances.append((q, s))
    except Exception as e:
        print(f"Execution Error in {func_name}: {e}")
        return None

    # C. Format Prompt
    final_prompt = FULL_QA_PROMPT.format(
        template_code=source_code,
        q1=instances[0][0], s1=instances[0][1],
        q2=instances[1][0], s2=instances[1][1],
        q3=instances[2][0], s3=instances[2][1]
    )

    # D. Prepare Tasks based on Flags
    tasks = []
    
    # 1. OpenAI
    if run_flags['openai']:
        tasks.append(call_openai(final_prompt))
    else:
        tasks.append(asyncio.sleep(0)) # Dummy task

    # 2. Anthropic
    if run_flags['anthropic']:
        tasks.append(call_anthropic(final_prompt))
    else:
        tasks.append(asyncio.sleep(0)) 

    # 3. Google
    if run_flags['google']:
        tasks.append(call_google(final_prompt))
    else:
        tasks.append(asyncio.sleep(0)) 

    # E. Execute
    results = await asyncio.gather(*tasks)
    
    # F. Structure Output (Handle skipped tasks)
    base_info = {
        "branch": branch,
        "area": area,
        "template_id": func_name,
        "timestamp": "2025-12-14"
    }

    # Map results back to models. 
    # If we skipped (asyncio.sleep), result is None.
    res_openai = results[0] if run_flags['openai'] else None
    res_anthropic = results[1] if run_flags['anthropic'] else None
    res_google = results[2] if run_flags['google'] else None
    
    return {
        "openai": {**base_info, "response": res_openai} if res_openai else None,
        "anthropic": {**base_info, "response": res_anthropic} if res_anthropic else None,
        "google": {**base_info, "response": res_google} if res_google else None
    }

async def main():
    print(" Starting Smart AI Tribunal ")
    
    # 1. Determine which models to run
    # If file exists (and has content), we assume it's done and skip it.
    run_openai = not (FILE_OPENAI.exists() and FILE_OPENAI.stat().st_size > 0)
    run_anthropic = not (FILE_ANTHROPIC.exists() and FILE_ANTHROPIC.stat().st_size > 0)
    run_google = not (FILE_GOOGLE.exists() and FILE_GOOGLE.stat().st_size > 0)
    
    run_flags = {
        'openai': run_openai,
        'anthropic': run_anthropic,
        'google': run_google
    }
    
    print(f"Status check:")
    print(f" - OpenAI: {'RUNNING' if run_openai else 'SKIPPING (File exists)'}")
    print(f" - Anthropic: {'RUNNING' if run_anthropic else 'SKIPPING (File exists)'}")
    print(f" - Google: {'RUNNING' if run_google else 'SKIPPING (File exists)'}")

    if not any(run_flags.values()):
        print("All files exist! Nothing to do.")
        return

    # 2. Gather Templates
    tasks = []
    branches = get_branches()
    for branch in branches:
        areas = get_areas(branch)
        for area in areas:
            files = get_template_files(branch, area)
            for file in files:
                funcs = load_template_functions(branch, area, file)
                for func_name, func_obj in funcs:
                    tasks.append((branch, area, file, func_name, func_obj))
    
    print(f"Found {len(tasks)} templates.")

    # 3. Open Files (Only open the ones we are writing to)
    f_openai = open(FILE_OPENAI, "a", encoding="utf-8") if run_openai else None
    f_anthropic = open(FILE_ANTHROPIC, "a", encoding="utf-8") if run_anthropic else None
    f_google = open(FILE_GOOGLE, "a", encoding="utf-8") if run_google else None

    # 4. Execute in Batches
    BATCH_SIZE = 5 
    for i in range(0, len(tasks), BATCH_SIZE):
        batch = tasks[i : i + BATCH_SIZE]
        
        # Pass flags to process_template
        coroutines = [process_template(*t, run_flags) for t in batch]
        batch_results = await asyncio.gather(*coroutines)
        
        # Write results
        for res in batch_results:
            if res:
                if run_openai and res["openai"]:
                    f_openai.write(json.dumps(res["openai"]) + "\n")
                if run_anthropic and res["anthropic"]:
                    f_anthropic.write(json.dumps(res["anthropic"]) + "\n")
                if run_google and res["google"]:
                    f_google.write(json.dumps(res["google"]) + "\n")
        
        print(f"Batch {i//BATCH_SIZE + 1} completed.")

    # 5. Cleanup
    if f_openai: f_openai.close()
    if f_anthropic: f_anthropic.close()
    if f_google: f_google.close()
    print(" Tribunal Complete ")

if __name__ == "__main__":
    asyncio.run(main())
