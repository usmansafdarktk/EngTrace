"""
annotate_gemini.py — EngTrace annotation using gemini-3.1-pro-preview

Usage:
    python annotate_gemini.py

Reads GEMINI_API_KEY from .env
Output folder: Annotator Gemini/
"""

import json
import os
import re
from pathlib import Path

from google import genai
from google.genai import types as genai_types
from dotenv import load_dotenv

from annotate_utils import (
    SYSTEM_PROMPT,
    build_user_prompt,
    run_all_files,
)

load_dotenv()

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

MODEL_NAME    = "gemini-2.5-pro"
ANNOTATOR_ID  = "Annotator Gemini"
OUTPUT_DIR    = Path(__file__).parent / "Annotator Gemini 2.5 Pro"
N_WORKERS     = 15

# ---------------------------------------------------------------------------
# Gemini client setup
# ---------------------------------------------------------------------------

API_KEY = os.getenv("GEMINI_API_KEY")
if not API_KEY:
    raise EnvironmentError("GEMINI_API_KEY not found in environment / .env file.")

_client = genai.Client(api_key=API_KEY)

# ---------------------------------------------------------------------------
# API call function
# ---------------------------------------------------------------------------

def _extract_json(text: str) -> dict:
    """Pull the first JSON object out of a response string."""
    try:
        return json.loads(text.strip())
    except json.JSONDecodeError:
        pass
    stripped = re.sub(r"^```(?:json)?\s*", "", text.strip(), flags=re.IGNORECASE)
    stripped = re.sub(r"\s*```$", "", stripped)
    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        pass
    m = re.search(r"\{[\s\S]+\}", text)
    if m:
        return json.loads(m.group(0))
    raise ValueError(f"Could not parse JSON from response:\n{text[:400]}")


def call_gemini(trace: dict) -> dict:
    """Call Gemini API for one trace; return raw annotation dict."""
    user_msg = build_user_prompt(trace, ANNOTATOR_ID)

    response = _client.models.generate_content(
        model=MODEL_NAME,
        contents=user_msg,
        config=genai_types.GenerateContentConfig(
            system_instruction=SYSTEM_PROMPT,
            temperature=0.1,
            response_mime_type="application/json",
        ),
    )
    return _extract_json(response.text)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_all_files(
        output_dir=OUTPUT_DIR,
        call_model_fn=call_gemini,
        annotator_id=ANNOTATOR_ID,
        n_workers=N_WORKERS,
    )
