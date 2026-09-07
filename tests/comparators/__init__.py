"""Phase 4 comparators — the machinery that turns a gold answer into a verdict.

Specification: ``docs/re-implementation-sep/phase4_comparators.md`` (D4.1).
Vocabulary:    ``docs/re-implementation-sep/phase4_vocabulary.md`` (D4.2).

Six ``kind`` values, one comparator each; ``compare_answer`` composes them for a
multipart answer (D-047).  Every comparator returns a :class:`Verdict` whose
``outcome`` is one of ``MATCH`` / ``MISMATCH`` / ``UNRESOLVED``, and
**``UNRESOLVED`` can never contribute to a pass** (D-049).
"""

from .verdict import Verdict, MATCH, MISMATCH, UNRESOLVED
from .kinds import (
    compare_numeric,
    compare_categorical,
    compare_sequence,
    compare_symbolic,
    compare_narrative,
    compare_check,
    COMPARATORS,
)
from .answer import compare_answer

__all__ = [
    "Verdict", "MATCH", "MISMATCH", "UNRESOLVED",
    "compare_numeric", "compare_categorical", "compare_sequence",
    "compare_symbolic", "compare_narrative", "compare_check",
    "COMPARATORS", "compare_answer",
]
