"""
annotate_claude.py — EngTrace annotation using claude-opus-4-7

Usage:
    python annotate_claude.py

Reads ANTHROPIC_API_KEY from .env
Output folder: Annotator Claude/
"""

import json
import os
import re
from pathlib import Path

import anthropic
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

MODEL_NAME    = "claude-opus-4-7"
ANNOTATOR_ID  = "Annotator Claude"
OUTPUT_DIR    = Path(__file__).parent / "Annotator Claude"
N_WORKERS     = 15
MAX_TOKENS    = 512    # annotation JSON is compact; 512 is sufficient

# ---------------------------------------------------------------------------
# Anthropic client setup
# ---------------------------------------------------------------------------

API_KEY = os.getenv("ANTHROPIC_API_KEY")
if not API_KEY:
    raise EnvironmentError("ANTHROPIC_API_KEY not found in environment / .env file.")

_client = anthropic.Anthropic(api_key=API_KEY)

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


def call_claude(trace: dict) -> dict:
    """Call Claude API for one trace; return raw annotation dict."""
    user_msg = build_user_prompt(trace, ANNOTATOR_ID)

    message = _client.messages.create(
        model=MODEL_NAME,
        max_tokens=MAX_TOKENS,
        system=SYSTEM_PROMPT,
        messages=[
            {"role": "user", "content": user_msg}
        ],
    )

    response_text = message.content[0].text
    return _extract_json(response_text)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_all_files(
        output_dir=OUTPUT_DIR,
        call_model_fn=call_claude,
        annotator_id=ANNOTATOR_ID,
        n_workers=N_WORKERS,
    )
