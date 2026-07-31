"""
annotate_gpt.py — EngTrace annotation using gpt-5-mini

Usage:
    python annotate_gpt.py

Reads OPENAI_API_KEY from .env
Output folder: Annotator GPT/
"""

import json
import os
import re
from pathlib import Path

from openai import OpenAI
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

MODEL_NAME    = "gpt-5-mini"
ANNOTATOR_ID  = "Annotator GPT"
OUTPUT_DIR    = Path(__file__).parent / "Annotator GPT"
N_WORKERS     = 15
MAX_TOKENS    = 512

# ---------------------------------------------------------------------------
# OpenAI client setup
# ---------------------------------------------------------------------------

API_KEY = os.getenv("OPENAI_API_KEY")
if not API_KEY:
    raise EnvironmentError("OPENAI_API_KEY not found in environment / .env file.")

_client = OpenAI(api_key=API_KEY)

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


def call_gpt(trace: dict) -> dict:
    """Call OpenAI API for one trace; return raw annotation dict."""
    user_msg = build_user_prompt(trace, ANNOTATOR_ID)

    response = _client.chat.completions.create(
        model=MODEL_NAME,
        temperature=0.1,
        max_tokens=MAX_TOKENS,
        response_format={"type": "json_object"},   # enforces JSON output
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user",   "content": user_msg},
        ],
    )

    response_text = response.choices[0].message.content
    return _extract_json(response_text)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_all_files(
        output_dir=OUTPUT_DIR,
        call_model_fn=call_gpt,
        annotator_id=ANNOTATOR_ID,
        n_workers=N_WORKERS,
    )
