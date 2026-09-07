"""Did the candidate *commit* to a label? — the clause-level replacement for a
word blocklist.

**This module exists because both Phase 4 reviewers broke the thing it
replaces, from opposite directions and in the same place.**

The first implementation asked one question: *does a label surface appear,
unnegated, with no listed hedge word in the 45 characters before it?* Three
separate defects followed from that shape, and none of them is a coverage gap
that a longer list would close:

* **Reviewer B, F1** — the hedge list is a 27-string blocklist and English fills
  its complement immediately. 12 of 15 unlisted hedges scored ``MATCH``
  (``I suspect``, ``Apparently``, ``should be``, ``looks linear``,
  ``Tentatively``, ``Perhaps``, ``Linear?``). The list does not close over its
  own morphology: ``seems`` and ``appears`` are listed, ``seemingly`` and
  ``apparently`` are not. **And D4.4 could not see it, because every hedged
  adversarial case used a phrase already in the list — the corpus samples the
  inside of the list it certifies.** That is D-034 one layer up.
* **Reviewer B, F2** — the window looked backwards only, so every *trailing*
  hedge passed on listed vocabulary: ``linear, I think``, ``Linear. It seems.``
* **Reviewer B, F3** — a conditional, an assumption, a refusal and *the item's
  own prompt wording* all scored ``MATCH``. ``We must determine whether the
  system is linear`` is a paraphrase of the question and was credited as an
  answer.
* **Reviewer E, F2 / Reviewer B, F4** — found independently — label selection
  ranked surfaces by *length across the whole span*, so a contrastive mention
  (``Unlike a nonlinear system…``) or a self-correction (``nonlinear -- no,
  wait, it is linear``) outranked the committed answer, in both directions.

The replacement asks a different question, which is Reviewer B's §5.3
prescription taken as written:

    **Which clause commits to a label, and is that clause an assertion?**

Four consequences, each closing one of the above:

1. Scope is the **clause**, bidirectionally — not N characters backwards.
2. **Modality is the feature, not the phrase.** A modal or evidential
   governing the copula is a hedge whether or not anyone listed it.
3. A label under ``if`` / ``assume`` / ``whether`` / ``unlike``, or in a
   question, is **not an assertion** and does not commit.
4. Across clauses the **last committed clause** wins; longest-surface-first
   applies strictly *within* one offset, where it is about containment and is
   correct.

The residual is measured, not assumed: see ``tests/comparators/reviewer_battery.py``.
"""

from __future__ import annotations

import re
from typing import Sequence

# --------------------------------------------------------------------------
# Clause segmentation
# --------------------------------------------------------------------------

#: Clause boundaries.  Sentence punctuation, the em/en-dash aside, and the
#: enumerators a label-tuple answer uses.
_CLAUSE_SPLIT = re.compile(
    r"(?:[.;!]|\?|\n|--+|\s—\s|(?:^|[\s,;])\(?[a-d][)\.]\s)"
)


def clauses(text: str) -> list[str]:
    """Split an answer span into clauses, keeping a trailing ``?`` marker.

    The question mark is preserved as part of the clause it ends, because
    ``Linear?`` is not a commitment and the mark is the only thing that says so.
    """
    out: list[str] = []
    pos = 0
    for m in _CLAUSE_SPLIT.finditer(text):
        piece = text[pos:m.start()]
        if m.group(0).startswith("?"):
            piece += "?"
        if piece.strip():
            out.append(piece.strip())
        pos = m.end()
    if text[pos:].strip():
        out.append(text[pos:].strip())
    return out


# --------------------------------------------------------------------------
# Non-assertion: the clause is not a claim about the system
# --------------------------------------------------------------------------

#: Openers that make a clause hypothetical, imperative or interrogative rather
#: than assertive.  A label inside one of these is not an answer.
NON_ASSERTING_OPENERS = (
    "if", "unless", "assume", "assuming", "suppose", "supposing",
    "whether", "consider", "let", "given that", "provided that", "in case",
    "were", "had", "should the", "note that", "recall that",
    "unlike", "whereas", "although", "though", "even if",
    # Discourse markers that open a qualification.  Reviewer E's F4 case
    # ("Note: Answer: not linear if the offset were nonzero") reaches this
    # point even after answer_span refuses to peel into it, because the caveat
    # is still inside the span -- so the clause test has to reject it too.
    # "but"/"however" are deliberately absent as *openers*: the clause splitter
    # does not cut on them, so a mid-clause "but" never reaches this test, and
    # listing them would only misfire on a genuine contrast.
    "note", "n.b.", "nb", "caveat", "caution", "disclaimer", "aside",
    "otherwise", "alternatively", "in contrast", "except",
)

#: Phrases that make a clause a restatement of the *task* rather than an answer.
#: ``We must determine whether the system is linear`` is the item's own prompt.
TASK_RESTATEMENTS = (
    "determine whether", "determine if", "we must determine", "we need to determine",
    "the question asks", "we are asked", "to find out whether", "check whether",
    "test whether", "verify whether",
)

_OPENER_RE = re.compile(
    r"^\W*(?:" + "|".join(re.escape(o) for o in NON_ASSERTING_OPENERS) + r")\b",
    re.IGNORECASE,
)
_TASK_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(t) for t in TASK_RESTATEMENTS) + r")\b",
    re.IGNORECASE,
)


def is_assertion(clause: str) -> tuple[bool, str]:
    """Is this clause a claim about the system? Returns ``(verdict, reason)``."""
    if clause.rstrip().endswith("?"):
        return False, "the clause is a question"
    m = _OPENER_RE.match(clause)
    if m:
        return False, f"the clause opens with {m.group(0).strip()!r} and is not an assertion"
    m = _TASK_RE.search(clause)
    if m:
        return False, f"the clause restates the task ({m.group(0)!r}), it does not answer it"
    return True, ""


# --------------------------------------------------------------------------
# Modality and evidentiality: the FEATURE, not the phrase list
# --------------------------------------------------------------------------

#: Epistemic modals.  ``must`` and ``can`` are deliberately absent: "the system
#: must be linear" is a strong commitment, not a hedge, and treating it as one
#: would trade Reviewer B's false accepts for false rejects.
MODALS = ("may", "might", "could", "should", "would")

#: Evidential and likelihood adverbs.  This is a *class*, not a list of
#: remembered phrases -- the ``-ly`` adverbs of likelihood are the productive
#: family Reviewer B's ``Seemingly`` / ``Apparently`` / ``Plausibly`` escaped
#: through, and the regex below also catches the family's members that nobody
#: has thought of yet.
EVIDENTIALS = (
    "apparently", "seemingly", "presumably", "plausibly", "tentatively",
    "perhaps", "possibly", "probably", "likely", "roughly", "arguably",
    "ostensibly", "supposedly", "conceivably", "maybe", "allegedly",
    "putatively", "notionally", "nominally",
)

#: Verbs of opinion, and copulas of appearance.  ``The system LOOKS linear``
#: and ``I SUSPECT the system is linear`` name a label without asserting it.
HEDGE_VERBS = (
    "suspect", "guess", "think", "believe", "feel", "reckon", "imagine",
    "would say", "am inclined", "lean towards", "lean toward",
)
APPEARANCE_COPULAS = (
    "seems", "seem", "seemed", "appears", "appear", "appeared",
    "looks", "look", "looked", "sounds", "sound",
)

#: Explicit statements of low confidence or of inability.
UNCERTAINTY = (
    "not sure", "not certain", "not 100%", "not entirely", "unsure", "unclear",
    "uncertain", "hard to say", "difficult to say", "cannot determine",
    "cannot complete", "unable to", "can't tell", "cannot tell", "no idea",
    "to be verified", "to be confirmed", "tbc", "it depends on how",
    "depends on the interpretation", "one could argue", "i am not sure",
)

_MODAL_RE = re.compile(
    r"\b(?:" + "|".join(MODALS) + r")\s+(?:\w+\s+){0,2}?(?:be|is|are|was|were)\b",
    re.IGNORECASE,
)
_EVIDENTIAL_RE = re.compile(
    r"\b(?:" + "|".join(EVIDENTIALS) + r")\b", re.IGNORECASE)
_HEDGE_VERB_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(v) for v in HEDGE_VERBS) + r")\b", re.IGNORECASE)
_APPEARANCE_RE = re.compile(
    r"\b(?:" + "|".join(APPEARANCE_COPULAS) + r")\b", re.IGNORECASE)
#: ``(?<!\w)`` / ``(?!\w)`` rather than ``\b``: several entries end in a
#: non-word character (``not 100%``), and ``\b`` there demands a word character
#: that will never appear, so the entry silently never fires -- which is how
#: "I am not 100% sure, but the system is linear" survived the first round of
#: fixes for this very finding.
_UNCERTAIN_RE = re.compile(
    r"(?<!\w)(?:" + "|".join(re.escape(u) for u in UNCERTAINTY) + r")(?!\w)",
    re.IGNORECASE)

_HEDGE_PROBES = (
    ("modal", _MODAL_RE),
    ("evidential", _EVIDENTIAL_RE),
    ("verb of opinion", _HEDGE_VERB_RE),
    ("copula of appearance", _APPEARANCE_RE),
    ("stated uncertainty", _UNCERTAIN_RE),
)


def hedge_markers(clause: str) -> list[str]:
    """Every hedge marker in ``clause``, as ``"<class>: <text>"``.

    Scope is the whole clause in **both** directions, which is Reviewer B's F2:
    English puts a hedge after the claim as often as before it, and a
    backwards-only window sees none of them.
    """
    out: list[str] = []
    for kind, rx in _HEDGE_PROBES:
        for m in rx.finditer(clause):
            out.append(f"{kind}: {m.group(0).lower()}")
    return sorted(set(out))


# --------------------------------------------------------------------------
# The commitment decision
# --------------------------------------------------------------------------


class Commitment:
    """Where a span commits to a label, or why it does not."""

    __slots__ = ("clause", "index", "reason")

    def __init__(self, clause: str | None, index: int, reason: str = ""):
        self.clause, self.index, self.reason = clause, index, reason

    def __bool__(self) -> bool:
        return self.clause is not None


def find_commitment(text: str, surfaces: Sequence[str]) -> Commitment:
    """The last clause of ``text`` that *asserts* one of ``surfaces``.

    "Last" rather than "longest surface" is Reviewer E's F2 and Reviewer B's F4,
    which are the same defect: a contrastive or self-corrected mention elsewhere
    in the span must not outrank the committed answer. Within one clause the
    longest surface still wins, because there the competition is containment.

    A **label-free hedge clause adjacent to the committed one governs it**, in
    either direction. ``Linear. It seems.`` and ``I cannot complete the test;
    the system is linear.`` both put the hedge in a clause of its own, and a
    per-clause test that ignored neighbours would credit both.
    """
    from .kinds import _label_hit  # noqa: PLC0415 - circular at import time

    cls = clauses(text)
    if not cls:
        return Commitment(None, -1, "the answer span is empty")

    labelled = [(i, c) for i, c in enumerate(cls) if _label_hit(c, surfaces)]
    if not labelled:
        return Commitment(None, -1, "no clause names a declared label")

    last_reason = ""
    for i, c in reversed(labelled):
        ok, why = is_assertion(c)
        if not ok:
            last_reason = why
            continue
        marks = hedge_markers(c)
        if marks:
            last_reason = f"hedged ({', '.join(marks)})"
            continue
        neighbour = _governing_neighbour(cls, i, surfaces)
        if neighbour:
            last_reason = f"hedged by an adjacent clause ({neighbour})"
            continue
        return Commitment(c, i)
    return Commitment(None, -1, last_reason or "no clause commits to a label")


def _governing_neighbour(cls: list[str], i: int, surfaces: Sequence[str]) -> str:
    """A hedge in a label-free clause immediately before or after clause ``i``.

    Only label-free neighbours govern. A neighbouring clause that names a label
    of its own is a separate statement, not a qualifier on this one -- otherwise
    ``The system is nonlinear. Actually, the system is linear.`` would have its
    committed second clause governed by its retracted first.
    """
    from .kinds import _label_hit  # noqa: PLC0415

    for j in (i - 1, i + 1):
        if not (0 <= j < len(cls)):
            continue
        if _label_hit(cls[j], surfaces):
            continue
        marks = hedge_markers(cls[j])
        if marks:
            return ", ".join(marks)
    return ""
