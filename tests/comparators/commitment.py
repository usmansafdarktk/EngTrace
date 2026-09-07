"""Did the candidate *commit* to a label? — governance, not distance.

**Version 2, after both reviewers broke version 1 in the same place.**

Version 1 replaced a 45-character hedge window with a ±1-*clause* hedge window
and a longest-surface rule with a last-*clause* rule, and called that "scope is
the clause". Both reviewers, in isolation, said it was not:

* **E, R2-F7** — the round-1 F2 case with **one period changed to a comma**.
  ``find_commitment`` picks *which clause*; ``_resolve_label`` then still took
  the last surface *within* it, and the splitter did not cut on commas. So
  ``The system is linear, unlike a nonlinear system, since both tests pass``
  was still scored ``not linear``.
* **E, R2-F8** — ``for j in (i-1, i+1)`` is not clause scope, it is **N = 1
  clause**: a magic constant of the same shape as the ``_NEG_WINDOW = 40`` it
  replaced. One intervening clause and a hedge escapes.
* **B, R2-F1.2** — four of five hedge probes were bare ``\\b…\\b`` scans of the
  whole clause with **no governance test at all**, and B named the through-line
  exactly: *character scope became clause scope; the missing thing both times is
  a marker→label relation.*
* **B, R2-F1.1 / E, R2-F9a** — found independently: ``Although X, Y`` and
  ``Given that X, Y`` **assert Y**. Version 1 matched the opener of a
  subordinate clause and read it as the mood of the whole sentence, refusing
  25 of 38 answers a grader credits.
* **B, R2-F1.3 / E, R2-F9c** — the neighbour rule was undirected, so
  ``Imagine a scaled input a*x[n]; … so the system is linear`` (the opening of
  the standard homogeneity proof) and ``… so the design is acceptable. The
  margin seems comfortable.`` were both refused.

So version 2 stops asking *how far away* a marker is and asks **what it
governs**:

1. A clause is segmented into a **subordinate** part and a **matrix** part at
   its subordinator. ``if``/``unless`` are **hypothetical** and suspend the
   whole sentence; ``although``/``given that``/``note that`` are **concessive
   or factive** and their matrix clause asserts. That one distinction closes
   E's R2-F7 and B's R2-F1.1 — which point in *opposite* directions, which is
   the sign of a root cause rather than a patch.
2. A hedge fires when it is **in the label's own segment**. Not the clause, the
   segment; and for all five probe classes, not just the modal.
3. A neighbouring clause governs only when it is a **bare epistemic comment** —
   label-free, short, anaphoric or subjectless — and only for the two classes
   that can be anaphoric. **No distance bound at all**, which is how the magic
   constant leaves rather than shrinking.

**What this mechanism still has no evidence for** (E, R2-F10): across all 2,200
archived answer spans the hedge probes fire **once**, and **zero times** in the
three kinds this module gates. It is a *policy* — spec §4.4 requires a hedge to
be non-committal — and it is validated only against the reviewers' own cases and
``recall_corpus.py``. It belongs on the residual-risk register beside the
symbolic rules, and ``derive_vocabulary.py`` prints its archive census so the
number cannot quietly be forgotten.
"""

from __future__ import annotations

import re
from typing import Sequence

# --------------------------------------------------------------------------
# 1.  Clause segmentation
# --------------------------------------------------------------------------

#: Sentence and clause boundaries.  The negative lookarounds matter: without
#: them the splitter cut inside decimals, abbreviations, initials and list
#: enumerators (E, R2-F11), which silently pushed hedges out of range and was
#: the distance amplifier for R2-F8.
#:
#:     "The gain is 2.5 and the system is linear"  ->  ['The gain is 2', '5 and …']
#:     "the system is linear i.e. additive"        ->  ['… linear i', 'e', 'additive']
_CLAUSE_SPLIT = re.compile(
    r"""(?:
          (?<!\d)\.(?!\d)(?<!\bi\.e)(?<!\be\.g)(?<!\bcf)     # . but not 2.5, i.e., e.g.
              (?<![A-Z])                                     # . but not an initial
          | [;!]
          | \?
          | \n
          | --+
          | \s—\s
          | (?:^|[\s;])\(?[a-d][)]\s                         # a) b) enumerators
        )""",
    re.VERBOSE,
)

_ABBREV = re.compile(r"\b(?:i\.e|e\.g|cf|vs|etc|Prof|Dr|Mr|Ms|Fig|Eq|No)\.", re.IGNORECASE)


def clauses(text: str) -> list[str]:
    """Split an answer span into clauses, keeping a trailing ``?``.

    The question mark is preserved on the clause it ends, because ``Linear?`` is
    not a commitment and the mark is the only thing that says so.
    """
    # Protect abbreviations from the sentence splitter, then restore.
    holes: list[str] = []

    def _stash(m: re.Match) -> str:
        holes.append(m.group(0))
        return f"\x00{len(holes) - 1}\x00"

    text = _ABBREV.sub(_stash, text)

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

    def _restore(s: str) -> str:
        return re.sub(r"\x00(\d+)\x00", lambda m: holes[int(m.group(1))], s)

    return [_restore(c) for c in out]


# --------------------------------------------------------------------------
# 2.  Subordinators: hypothetical suspends, concessive does not
# --------------------------------------------------------------------------

#: **Hypothetical.** The matrix clause is asserted only *conditionally*, so the
#: whole sentence commits to nothing.  ``If additivity holds, the system is
#: linear`` states a rule, not a verdict (B, round-1 F3).
HYPOTHETICAL = (
    "if", "unless", "provided that", "provided", "assuming", "assume",
    "suppose", "supposing", "in case", "were", "had", "should the",
    "imagine", "consider", "let", "say that", "in the event",
)

#: **Concessive, causal and factive.** The subordinate part is backgrounded and
#: the matrix part **is asserted**.  ``Although the equation contains a delay,
#: the system is linear`` commits to *linear*; ``Note that X`` and ``Recall that
#: X`` are factive and assert X outright.
#:
#: Version 1 lumped these with the hypotheticals and refused 25 of 38 credited
#: answers (B, R2-F1.1; E, R2-F9a).
CONCESSIVE = (
    "although", "though", "even though", "while", "whilst", "whereas",
    "given that", "given", "since", "because", "as", "unlike",
    "in contrast to", "despite", "regardless of", "aside from", "apart from",
    "besides",
)

#: **Factive.** The complement is asserted outright, so the matrix is what
#: follows the opener -- not a post-comma remainder, which a factive need not
#: have.  ``Note that the system is linear`` is a preface with no comma at all,
#: and treating it as a concessive left it "a bare subordinate clause with no
#: matrix to assert" (Reviewer B, R2-F1.1: `note that` and `recall that` are
#: *additionally* factive, which the first draft of this rewrite missed).
FACTIVE = (
    "note that", "notes that", "recall that", "observe that", "notice that",
    "remember that", "it follows that", "we see that", "we conclude that",
    "conclude that", "it is clear that", "clearly",
)

#: Openers that make a clause a **restatement of the task** rather than an
#: answer.  ``We must determine whether the system is linear`` is the item's
#: own prompt wording (B, round-1 F3).
TASK_RESTATEMENTS = (
    "determine whether", "determine if", "we must determine", "we need to determine",
    "the question asks", "we are asked", "to find out whether", "check whether",
    "test whether", "verify whether", "we must decide", "our task is",
)

#: Bare discourse openers with no complement clause -- a caveat rather than a
#: preface.  ``Note:`` (colon, not ``that``) heads a qualification; ``Note that``
#: is factive and lives in CONCESSIVE above.  The two are distinguished by what
#: follows, which is the distinction E's R2-F9a said version 1 could not make.
BARE_QUALIFIERS = ("note", "n.b.", "nb", "caveat", "caution", "disclaimer", "aside")

_HYPO_RE = re.compile(
    r"^\W*(?:" + "|".join(re.escape(o) for o in HYPOTHETICAL) + r")\b", re.IGNORECASE)
_CONC_RE = re.compile(
    r"^\W*(?:" + "|".join(re.escape(o) for o in CONCESSIVE) + r")\b", re.IGNORECASE)
_FACTIVE_RE = re.compile(
    r"^\W*(?:" + "|".join(re.escape(o) for o in FACTIVE) + r")\b", re.IGNORECASE)
_TASK_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(t) for t in TASK_RESTATEMENTS) + r")\b", re.IGNORECASE)
_BARE_QUAL_RE = re.compile(
    r"^\W*(?:" + "|".join(re.escape(q) for q in BARE_QUALIFIERS) + r")\s*[:\-,]", re.IGNORECASE)

#: A subordinator appearing *after* a comma opens a trailing subordinate part:
#: ``The system is linear, unlike a nonlinear system, since both tests pass.``
#: Everything from that comma on is backgrounded, which is E's R2-F7.
#:
#: ``which`` and ``whose`` are deliberately absent. A **non-restrictive
#: relative asserts its content**: "The deflection is 18.4 mm, which is
#: acceptable" states the verdict *inside* the relative clause, and
#: backgrounding it deleted the verdict from a `check` answer whose gold uses
#: the very same construction. Only a genuine contrast or condition
#: backgrounds; an appositive elaborates.
_TRAILING_SUB_RE = re.compile(
    r",\s*(?:" + "|".join(re.escape(o) for o in CONCESSIVE + HYPOTHETICAL + FACTIVE
                          + ("rather than", "as opposed to")) + r")\b",
    re.IGNORECASE,
)


def segment(clause: str) -> tuple[str | None, str]:
    """Split a clause into ``(matrix, reason_if_none)``.

    ``matrix`` is the part of the clause that **asserts**, or ``None`` when the
    clause asserts nothing. This is the whole of the concessive/hypothetical
    distinction, and it is what both reviewers' opposite-direction findings
    reduce to.
    """
    c = clause.strip()
    if not c:
        return None, "the clause is empty"
    if c.rstrip().endswith("?"):
        return None, "the clause is a question"
    if _TASK_RE.search(c):
        m = _TASK_RE.search(c)
        return None, f"the clause restates the task ({m.group(0)!r}), it does not answer it"
    if _BARE_QUAL_RE.match(c):
        return None, "the clause is a bare qualifier, not an answer"
    if _HYPO_RE.match(c):
        return None, (f"the clause is hypothetical (opens with "
                      f"{_HYPO_RE.match(c).group(0).strip()!r}); its matrix is asserted "
                      f"only conditionally")

    # A factive: "Note that Y" -> the matrix is Y, with no comma required.
    m = _FACTIVE_RE.match(c)
    if m:
        c = c[m.end():].strip()
        if not c:
            return None, "the factive has no complement"

    # A fronted concessive: "Although X, Y" -> the matrix is Y.
    if _CONC_RE.match(c):
        comma = c.find(",")
        if comma == -1:
            return None, ("the clause is a bare subordinate clause with no matrix "
                          "to assert")
        c = c[comma + 1:].strip()
        if not c:
            return None, "the clause is a bare subordinate clause"
        # The matrix may itself be hypothetical: "Although X, if Y then Z".
        if _HYPO_RE.match(c):
            return None, "the matrix clause is itself hypothetical"

    # A comma-delimited subordinate part is backgrounded.  It is **excised, not
    # truncated to**, because the matrix often resumes after it:
    #
    #   "The system is neither memoryless, because the output depends on past
    #    inputs, nor causal."
    #
    # Truncating at the first subordinator drops "nor causal" and the second
    # slot loses its value -- which is how the first draft of this rewrite broke
    # the `neither ... nor` cases Reviewer E's R2-F5 fix had just closed.  An
    # interpolation has a closing comma; a genuinely trailing clause does not,
    # and only that one truncates.
    for _ in range(4):
        m = _TRAILING_SUB_RE.search(c)
        if not m:
            break
        close = c.find(",", m.end())
        c = (c[:m.start()] + c[close:]).strip() if close != -1 else c[:m.start()].strip()
    return (c, "") if c else (None, "the clause has no matrix part")


# --------------------------------------------------------------------------
# 3.  Hedge classes, and the governance relation they were missing
# --------------------------------------------------------------------------

#: Epistemic modals.  ``must`` and ``can`` are deliberately absent: "the system
#: must be linear" is a strong commitment.  Both reviewers confirmed the
#: exclusion is right.
MODALS = ("may", "might", "could", "should", "would")

#: Evidential and likelihood adverbs -- a productive class, not a remembered
#: list, which is what let ``Seemingly`` / ``Apparently`` / ``Plausibly`` through
#: the blocklist.
#:
#: ``roughly`` and ``nominally`` were here in version 1 and are **removed**:
#: both reviewers found independently that they are *precision* terms in
#: engineering register, not epistemic ones ("roughly 72% of the 25 mm limit"),
#: and ``check`` is the one gated kind whose answers carry quantities
#: (E, R2-F9b; B, R2-F1.2).
#: ``roughly speaking`` and ``broadly speaking`` stay, as the multi-word
#: discourse hedges they are; it is the **bare** approximator that had to go.
#: "Roughly speaking, the system is linear" hedges; "roughly 72% of the limit"
#: does not, and only listing the phrase separates them.
EVIDENTIALS = (
    "apparently", "seemingly", "presumably", "plausibly", "tentatively",
    "perhaps", "possibly", "probably", "likely", "arguably", "ostensibly",
    "supposedly", "conceivably", "maybe", "allegedly", "putatively",
    "roughly speaking", "broadly speaking", "loosely speaking",
    "more or less", "to a first approximation",
)

#: Verbs of opinion.  These open a *derivation* as often as they hedge
#: ("Imagine a scaled input a*x[n]" is the standard homogeneity proof), so they
#: fire only inside the label's own segment and never from a neighbour clause.
HEDGE_VERBS = (
    "suspect", "guess", "think", "believe", "feel", "reckon",
    "would say", "am inclined", "lean towards", "lean toward",
)

#: Copulas of appearance.  Anaphoric enough to govern from a neighbour clause.
APPEARANCE_COPULAS = (
    "seems", "seem", "seemed", "appears", "appear", "appeared",
    "looks", "look", "looked", "sounds", "sound",
)

#: Statements of **inability or low confidence**.  Split out from version 1's
#: single ``UNCERTAINTY`` bag on Reviewer B's R2-F1.3: these are the only
#: markers that legitimately govern *forwards* from a preceding clause, because
#: "I cannot tell" said before an answer qualifies that answer, while
#: "Imagine …" said before it does not.
INABILITY = (
    "not sure", "not certain", "not 100%", "not entirely sure", "unsure",
    "uncertain", "hard to say", "difficult to say", "cannot determine",
    "cannot complete", "cannot really tell", "unable to", "can't tell",
    "cannot tell", "no idea", "i am not sure", "cannot be determined",
    "to be verified", "to be confirmed", "unclear",
    "it depends on how", "depends on the interpretation", "one could argue",
)


def _alt(words: Sequence[str]) -> str:
    return "|".join(re.escape(w) for w in words)


_MODAL_RE = re.compile(
    r"\b(?:" + _alt(MODALS) + r")\s+(?:\w+\s+){0,2}?(?:be|is|are|was|were)\b", re.I)
_EVIDENTIAL_RE = re.compile(r"(?<!\w)(?:" + _alt(EVIDENTIALS) + r")(?!\w)", re.I)
_HEDGE_VERB_RE = re.compile(r"(?<!\w)(?:" + _alt(HEDGE_VERBS) + r")(?!\w)", re.I)
_APPEARANCE_RE = re.compile(r"(?<!\w)(?:" + _alt(APPEARANCE_COPULAS) + r")(?!\w)", re.I)
_INABILITY_RE = re.compile(r"(?<!\w)(?:" + _alt(INABILITY) + r")(?!\w)", re.I)

PROBES = (
    ("modal", _MODAL_RE),
    ("evidential", _EVIDENTIAL_RE),
    ("verb of opinion", _HEDGE_VERB_RE),
    ("copula of appearance", _APPEARANCE_RE),
    ("stated uncertainty", _INABILITY_RE),
)

#: The classes that may govern from a **neighbouring** clause.  Modals and verbs
#: of opinion are excluded: ``This should be clear from the square`` is
#: confidence, and ``Imagine a scaled input`` opens a proof (B, R2-F1.3).
NEIGHBOUR_CLASSES = frozenset({"copula of appearance", "stated uncertainty"})


def hedge_markers(text: str, classes: frozenset[str] | None = None) -> list[str]:
    """Hedge markers in ``text``, optionally restricted to some classes."""
    out: list[str] = []
    for kind, rx in PROBES:
        if classes is not None and kind not in classes:
            continue
        for m in rx.finditer(text):
            out.append(f"{kind}: {m.group(0).lower()}")
    return sorted(set(out))


# --------------------------------------------------------------------------
# 4.  Bare epistemic comments -- what may govern from a distance
# --------------------------------------------------------------------------

_ANAPHORIC_SUBJECT = re.compile(
    r"^\W*(?:it|this|that|these|those|i|we|there)\b", re.IGNORECASE)
_MAX_COMMENT_WORDS = 6


def is_bare_comment(clause: str, surfaces: Sequence[str]) -> bool:
    """Is this clause a bare epistemic remark rather than a claim of its own?

    ``It seems.`` / ``I am not sure.`` / ``Hard to say, though.`` qualify a
    neighbouring commitment. ``The margin seems comfortable.`` and ``This should
    be clear from the square.`` do not -- they are claims about something else
    that happen to contain a marker.

    The test is **grammatical, not positional**, which is how version 1's
    ``for j in (i-1, i+1)`` magic constant leaves rather than shrinking: a bare
    comment governs from any distance, and a full clause governs from none.
    """
    from .kinds import _label_hit  # noqa: PLC0415 - circular at import time

    if _label_hit(clause, surfaces):
        return False
    if len(clause.split()) > _MAX_COMMENT_WORDS:
        return False
    stripped = re.sub(r"^\W*(?:but|and|however|though|although)\b", "", clause,
                      flags=re.IGNORECASE).strip()
    if _ANAPHORIC_SUBJECT.match(stripped):
        return True
    # Subjectless: "Hard to say", "Not sure", "Difficult to tell".
    return not re.match(r"^\W*(?:the|a|an|both|each|every|all|its|his|her|their)\b",
                        stripped, re.IGNORECASE)


# --------------------------------------------------------------------------
# 5.  The commitment decision
# --------------------------------------------------------------------------


class Commitment:
    """Where a span commits to a label, or why it does not."""

    __slots__ = ("clause", "index", "reason")

    def __init__(self, clause: str | None, index: int, reason: str = ""):
        self.clause, self.index, self.reason = clause, index, reason

    def __bool__(self) -> bool:
        return self.clause is not None


def label_segment(text: str, surfaces: Sequence[str]) -> str | None:
    """The asserting segment of ``text`` that carries a label, or None."""
    from .kinds import _label_hit  # noqa: PLC0415

    matrix, _ = segment(text)
    if matrix is None:
        return None
    return matrix if _label_hit(_strip_label_parens(matrix), surfaces) else None


_PAREN_RE = re.compile(r"\(([^()]*)\)")


def _strip_label_parens(text: str, surfaces: Sequence[str] = ()) -> str:
    """Drop parentheticals that talk *about the labels* rather than commit.

    ``The system is linear (a nonlinear system would fail additivity)`` commits
    to *linear*: the parenthetical is a gloss about the other label and must
    neither outrank the commitment (E, R2-F7) nor hedge it -- its ``would be``
    is a modal about a hypothetical answer, not about this one (E, R2-F12).

    A parenthetical naming **no** label is a genuine qualifier and stays, which
    is what keeps ``Linear (unclear)`` non-committal (B, round-1 F2).
    """
    from .kinds import _label_hit  # noqa: PLC0415

    if not surfaces:
        return _PAREN_RE.sub(" ", text)
    return _PAREN_RE.sub(
        lambda m: " " if _label_hit(m.group(1), surfaces) else m.group(0), text)


def find_commitment(text: str, surfaces: Sequence[str]) -> Commitment:
    """The last clause of ``text`` whose **matrix segment** asserts a label.

    Returns the matrix segment, not the whole clause, so the caller resolves the
    label inside the part that actually asserts.
    """
    from .kinds import _label_hit  # noqa: PLC0415

    cls = clauses(text)
    if not cls:
        return Commitment(None, -1, "the answer span is empty")

    labelled = [(i, c) for i, c in enumerate(cls) if _label_hit(c, surfaces)]
    if not labelled:
        return Commitment(None, -1, "no clause names a declared label")

    last_reason = ""
    for i, c in reversed(labelled):
        matrix, why = segment(c)
        if matrix is None:
            last_reason = why
            continue
        bare = _strip_label_parens(matrix, surfaces)
        if not _label_hit(bare, surfaces):
            last_reason = "the label is inside a backgrounded or parenthetical part"
            continue
        # Governance is now local: the marker must be in the label's own
        # asserting segment, for every probe class -- not merely somewhere in
        # the clause (B, R2-F1.2).
        marks = hedge_markers(bare)
        if marks:
            last_reason = f"hedged ({', '.join(marks)})"
            continue
        neighbour = _governing_comment(cls, i, surfaces)
        if neighbour:
            last_reason = f"qualified by a bare comment ({neighbour})"
            continue
        return Commitment(bare, i)
    return Commitment(None, -1, last_reason or "no clause commits to a label")


def _governing_comment(cls: list[str], i: int, surfaces: Sequence[str]) -> str:
    """A bare epistemic comment anywhere in the span that qualifies clause ``i``.

    **No distance bound.**  A bare comment governs from any distance (E's
    R2-F8: version 1's ``(i-1, i+1)`` let a hedge escape by one clause), and a
    clause that is not a bare comment governs from none (B's R2-F1.3 and E's
    R2-F9c: version 1 let an unrelated trailing sentence retract a commitment).

    Direction matters for one class. A *preceding* comment governs only if it
    states **inability** -- "I cannot tell" before an answer qualifies it, while
    "It seems additive" before a derivation does not.
    """
    for j, c in enumerate(cls):
        if j == i or not is_bare_comment(c, surfaces):
            continue
        classes = NEIGHBOUR_CLASSES if j > i else frozenset({"stated uncertainty"})
        marks = hedge_markers(c, classes)
        if marks:
            return ", ".join(marks)
    return ""


def is_assertion(clause: str) -> tuple[bool, str]:
    """Back-compatible wrapper: does this clause assert anything at all?"""
    matrix, why = segment(clause)
    return (matrix is not None), why
