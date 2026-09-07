"""The verdict a comparator returns, and the three-valued outcome it carries.

The three-valued outcome is the whole of Decision 3 (D-049).  A comparator that
can only say yes or no must guess when a candidate does not commit, and every
guess it makes is either a false accept -- which credits reasoning that did not
happen, the criticism this effort exists to answer -- or a false reject, which
scores phrasing rather than reasoning.  ``UNRESOLVED`` is the refusal to guess.

It is modelled on D3.4 |S|8B.10's ``unchecked`` channel: an unchecked relation may
never contribute to a pass, so a node carrying one cannot report a bare PASS.
Here, an ``UNRESOLVED`` verdict is never a match, is never silently a mismatch
either, and is reported as its own rate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

MATCH = "MATCH"
MISMATCH = "MISMATCH"
UNRESOLVED = "UNRESOLVED"

_OUTCOMES = (MATCH, MISMATCH, UNRESOLVED)


@dataclass
class Verdict:
    """The result of comparing one candidate answer against one gold answer.

    ``outcome``           one of MATCH / MISMATCH / UNRESOLVED
    ``kind``              the comparator that produced it
    ``reason``            why, in one line; required for MISMATCH and UNRESOLVED
    ``gold_canonical``    gold after normalisation, for audit
    ``cand_canonical``    candidate after normalisation, for audit
    ``observations``      non-scoring notes (D3.4 |S|7.2's count report is one)
    ``parts``             sub-verdicts, for a composite answer
    """

    outcome: str
    kind: str
    reason: str = ""
    gold_canonical: Any = None
    cand_canonical: Any = None
    observations: list[str] = field(default_factory=list)
    parts: list["Verdict"] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.outcome not in _OUTCOMES:
            raise ValueError(f"outcome must be one of {_OUTCOMES}, got {self.outcome!r}")
        if self.outcome in (MISMATCH, UNRESOLVED) and not self.reason:
            raise ValueError(f"a {self.outcome} verdict must carry a reason")

    @property
    def is_match(self) -> bool:
        """True only for MATCH.  UNRESOLVED is never a pass."""
        return self.outcome == MATCH

    @property
    def scores(self) -> bool:
        """True when this verdict may enter an accuracy numerator or denominator
        as a decided case.  UNRESOLVED does not; it is reported separately."""
        return self.outcome in (MATCH, MISMATCH)

    def __str__(self) -> str:  # pragma: no cover - display only
        head = f"{self.outcome}({self.kind})"
        return head if not self.reason else f"{head}: {self.reason}"


def match(kind: str, **kw: Any) -> Verdict:
    return Verdict(MATCH, kind, **kw)


def mismatch(kind: str, reason: str, **kw: Any) -> Verdict:
    return Verdict(MISMATCH, kind, reason=reason, **kw)


def unresolved(kind: str, reason: str, **kw: Any) -> Verdict:
    return Verdict(UNRESOLVED, kind, reason=reason, **kw)
