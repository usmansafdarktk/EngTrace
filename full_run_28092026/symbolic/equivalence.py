"""The symbolic equivalence check, moved beside the answer check at ANSWER FINAL.

The code is evaluator_pilot_17092026/evaluators/symbolic_equivalence.py: answer.py imports it and pins its
LF-normalised SHA-256, so the hash of answer.py that score.py records in every store's CONFIG covers it. This module
re-exports it, so the scripts of this folder (validate.py) and anyone importing full_run_28092026.symbolic.equivalence
read the one implementation.

    python -m full_run_28092026.symbolic.equivalence --selftest     # FREE: the check's own cases
"""
from __future__ import annotations

import sys
from pathlib import Path

_EVAL = str(Path(__file__).resolve().parents[2] / 'evaluator_pilot_17092026' / 'evaluators')
if _EVAL not in sys.path:
    sys.path.insert(0, _EVAL)

import symbolic_equivalence as _impl  # noqa: E402
from symbolic_equivalence import *  # noqa: E402,F401,F403


def __getattr__(name: str):                  # the private helpers too (PEP 562)
    return getattr(_impl, name)


if __name__ == '__main__':
    if '--selftest' not in sys.argv:
        raise SystemExit('usage: python -m full_run_28092026.symbolic.equivalence --selftest   # FREE')
    raise SystemExit(_impl.selftest())
