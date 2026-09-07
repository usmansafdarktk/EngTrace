"""Composite answers -- the ``multipart`` question, answered (D-047).

**``multipart`` is not a seventh ``kind``.**  It is a property of the *answer*:
an answer has an ordered list of ``parts``, and each part carries exactly one of
the six kinds.  The six stand; the answer schema grows.

The evidence is inside this phase.  ``system_properties_memory_causality`` is
typed ``classification`` by the inventory, is described by the spec as a
"canonical label tuple", and is in fact a **two-part answer whose parts are both
categorical**.  Under a seventh-kind reading it would have to be typed
``multipart``, and its comparator would then have no way to say "each part is
categorical, normalise each against the categorical vocabulary" -- the six kinds
would stop composing at exactly the point they are most useful.  Under the
composition reading it is ``parts=[categorical, categorical]`` and it reuses the
categorical vocabulary unchanged.

The 32 ``multipart`` templates decompose the same way; sampled, every one is an
ordered tuple of numeric parts (``volumetric_flow_rate``: flow rate and
velocity) or a numeric part restated in a second unit
(``undamped_natural_frequency_torsional``: "51.184 rad/s **or** 8.146 Hz").
Those two are *different structures*, and the distinction is what ``mode``
carries below.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

from .kinds import COMPARATORS, PropertySlot
from .verdict import MATCH, MISMATCH, UNRESOLVED, Verdict

#: How the parts of a composite answer combine.
#:
#: ``all``  every part is required and every part must match (conjunction).
#:          ``volumetric_flow_rate``: flow rate *and* velocity.
#: ``any``  the parts are alternative renderings of one quantity; matching any
#:          one is matching the answer (disjunction).
#:          ``undamped_natural_frequency_torsional``: rad/s *or* Hz.
#:
#: Conflating the two is a real scoring error in both directions: read ``any``
#: as ``all`` and a model that answers in rad/s alone is marked wrong; read
#: ``all`` as ``any`` and a model that gets the flow rate right and the velocity
#: wrong takes full credit.
MODES = ("all", "any")


@dataclass
class Part:
    """One component of an answer."""

    name: str
    kind: str
    #: Keyword arguments for the part's comparator (label set, precision, ...).
    options: dict[str, Any] = field(default_factory=dict)
    #: Optional slicer: given the full answer span, return this part's text.
    #: Defaults to the whole span, which is correct when the part's own
    #: comparator can locate itself (a label set, a sequence's braces).
    select: Callable[[str], str] | None = None

    def __post_init__(self) -> None:
        if self.kind not in COMPARATORS:
            raise ValueError(f"unknown kind {self.kind!r}; the six are {sorted(COMPARATORS)}")


@dataclass
class AnswerSpec:
    """The gold side's declaration of an answer's shape.

    This is **problem data** in |S|3.8's sense: declared once, on the gold answer,
    and never restated by a candidate.  A candidate cannot change the number of
    parts an answer has, the kind of a part, or the tolerance a part is judged
    at -- which is what stops a candidate from grading itself.
    """

    parts: list[Part]
    mode: str = "all"
    #: Partial credit is *reported*, never folded into the match verdict.
    #: P6: a partially-correct multipart answer is not a correct answer, and a
    #: comparator that says otherwise weakens 32 items silently.
    def __post_init__(self) -> None:
        if self.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}")
        if not self.parts:
            raise ValueError("an answer must have at least one part")


def compare_answer(gold: str, candidate: str, spec: AnswerSpec) -> Verdict:
    """Compare a candidate answer against gold under ``spec``.

    The composite verdict is three-valued like its parts, and the precedence is
    deliberate:

    * ``mode="all"``: any UNRESOLVED part makes the whole UNRESOLVED, because a
      part we could not decide might have been the one that was wrong.
      Resolving it as MISMATCH would report a model error we did not observe.
    * ``mode="any"``: a single MATCH is a MATCH even if a sibling is
      UNRESOLVED, because one alternative rendering has already settled it.

    ``observations`` always carries the per-part outcome string, so a partially
    correct multipart answer is legible in the results without being scored as
    a partial pass.
    """
    verdicts: list[Verdict] = []
    for p in spec.parts:
        g = p.select(gold) if p.select else gold
        c = p.select(candidate) if p.select else candidate
        v = COMPARATORS[p.kind](g, c, **p.options)
        v.kind = f"{p.name}:{v.kind}"
        verdicts.append(v)

    tally = f"{sum(v.is_match for v in verdicts)}/{len(verdicts)} parts matched"
    obs = [tally] + [f"{v.kind} -> {v.outcome}" for v in verdicts]

    if spec.mode == "any":
        if any(v.is_match for v in verdicts):
            return Verdict(MATCH, "composite[any]", parts=verdicts, observations=obs)
        if all(v.outcome == MISMATCH for v in verdicts):
            return Verdict(MISMATCH, "composite[any]",
                           reason="no alternative rendering matched",
                           parts=verdicts, observations=obs)
        return Verdict(UNRESOLVED, "composite[any]",
                       reason="no part matched and at least one was undecidable",
                       parts=verdicts, observations=obs)

    if any(v.outcome == UNRESOLVED for v in verdicts):
        bad = [v for v in verdicts if v.outcome == UNRESOLVED]
        return Verdict(UNRESOLVED, "composite[all]",
                       reason="; ".join(f"{v.kind}: {v.reason}" for v in bad),
                       parts=verdicts, observations=obs)
    if all(v.is_match for v in verdicts):
        return Verdict(MATCH, "composite[all]", parts=verdicts, observations=obs)
    wrong = [v.kind for v in verdicts if not v.is_match]
    return Verdict(MISMATCH, "composite[all]",
                   reason=f"part(s) differ: {', '.join(wrong)}",
                   parts=verdicts, observations=obs)


# --------------------------------------------------------------------------
# Bindings for the four Phase 4 templates.
#
# These are declarations, not code paths: each names the kind and the label set
# for one template's answer.  A fifth template is fitted by adding a binding
# here and editing no comparator -- which is the property D4.7 measures for the
# Phase 3 node types and the same property, one layer up.
# --------------------------------------------------------------------------

MEMORYLESS_SLOT = PropertySlot(
    name="memoryless",
    positive=("memoryless",),
    negative=("has memory", "with memory", "have memory"),
)
CAUSAL_SLOT = PropertySlot(
    name="causal",
    positive=("causal",),
    negative=("noncausal", "non-causal", "anticausal", "anti-causal"),
)
LINEARITY_LABELS = {
    "linear": ["linear"],
    "not linear": ["nonlinear", "non-linear"],
}

PHASE4_BINDINGS: dict[str, dict[str, Any]] = {
    "template_system_property_linearity": {
        "kind": "categorical",
        "options": {"labels": LINEARITY_LABELS},
    },
    "template_system_properties_memory_causality": {
        "kind": "categorical[tuple]",
        "options": {"slots": [MEMORYLESS_SLOT, CAUSAL_SLOT]},
    },
    "template_signal_operations": {
        "kind": "sequence",
        "options": {"require_origin": True},
    },
    "template_incompressible_continuity": {
        "kind": "symbolic",
        "options": {"symbols": ("x", "y"), "allow_arbitrary": True},
    },
}


def compare_template(template_id: str, gold: str, candidate: str) -> Verdict:
    """Score a candidate for one of the four Phase 4 templates."""
    from .kinds import compare_label_tuple  # noqa: PLC0415

    b = PHASE4_BINDINGS[template_id]
    fn = compare_label_tuple if b["kind"] == "categorical[tuple]" else COMPARATORS[b["kind"]]
    return fn(gold, candidate, **b["options"])
