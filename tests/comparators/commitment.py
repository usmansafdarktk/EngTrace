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

import os
import re
from typing import Sequence

# --------------------------------------------------------------------------
# 0.  The hedge policy -- ADVISORY by default (Reviewer E, round 4, rec. 3c)
# --------------------------------------------------------------------------

#: ``"advisory"`` (default) -- a detected hedge is **annotated and never
#: scored**. ``"enforce"`` -- a hedge makes the answer ``UNRESOLVED``, which is
#: what spec §4.2/§4.4 asks for and what versions 1-3 did.
#:
#: **Why the default changed, and it is a `SPEC-CHANGE`.**
#:
#: The commitment census splits one number the phase had been treating as two
#: halves of the same thing:
#:
#:     hedge markers          2 in 2,200 archived spans, 0 in a gated kind
#:     subordinator openers   43 (concessive 31, hypothetical 11)
#:
#: Two mechanisms in one module with **opposite evidence**. Reviewer E ablated
#: the hedge layer, leaving segmentation untouched:
#:
#:     ==================  ========  =========
#:                            full    ablated
#:     positive frames      351/351   351/351
#:     total                507/507   477/507
#:     ==================  ========  =========
#:
#: The whole layer buys **30 synthetic negative controls the reviewers wrote,
#: and zero archive verdicts and zero positive recall** -- for ~250 lines and
#: 13 of E's 20 findings across four rounds. Meanwhile it can, and repeatedly
#: did, silently mark a *correct* answer wrong. D4.1 §1 says a false reject is
#: as unacceptable as a false accept; a layer with no observed instances that
#: produces them is not paying for itself.
#:
#: **The policy is not abandoned, it is demoted to a reported quantity** -- the
#: same treatment ``narrative`` gets (D-051). Hedges are still detected, still
#: named in ``observations``, and the census still prints the rate. If that rate
#: ever becomes non-trivial the evidence for enforcing exists, and flipping this
#: constant is the whole change. Recorded as **D-056**.
HEDGE_POLICY = os.environ.get("ENGTRACE_HEDGE_POLICY", "advisory")

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

#: Inversion hypotheticals -- ``Were the system linear, …``, ``Had the offset
#: been zero, …``.  These are hypothetical **only clause-initially**: elsewhere
#: ``were`` and ``had`` are ordinary past tense.
_HYPO_INITIAL = ("were", "had", "should the")

#: Up to four words of fronted adverbial may precede every *other* subordinator.
#: "For now I will assume X" uses a word already in HYPOTHETICAL and escaped
#: only because the pattern was anchored at position 0 -- Reviewer B's round-4
#: point that *anchoring, not vocabulary*, is what decides, and ``_TASK_RE``
#: already allowed three words for the same reason.
#:
#: Relaxing the anchor for **all** of them immediately over-reached: "Given that
#: both defining conditions **were** checked above, …" matched on the ordinary
#: past tense four words in, and refused a concessive the rewrite exists to
#: credit.  That is the fix's own new surface, caught by the recall corpus one
#: run later, and it is why the two classes are separated rather than merged.
_HYPO_RE = re.compile(
    r"^\W*(?:"
    + "|".join(re.escape(o) for o in _HYPO_INITIAL)
    + r"|(?:\w+\s+){0,4}?(?:"
    + "|".join(re.escape(o) for o in HYPOTHETICAL if o not in _HYPO_INITIAL)
    + r"))\b", re.IGNORECASE)
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
    # Anchored to the clause opening (Reviewer E, R3-F22): unanchored, a
    # trailing "... which is what we needed to determine whether to accept"
    # would void a committed answer.
    m = _TASK_RE.match(c)
    if m:
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
    #
    # **The comma is optional.**  "Since both tests pass the system is linear"
    # is the same sentence without it, and requiring one made B's R2-F1.1 fix
    # conditional on ORTHOGRAPHY -- which the recall corpus could not see,
    # because every frame it generates supplies the comma (Reviewer E, R3-F17).
    # Without a comma the subordinate clause runs to the finite verb of the
    # matrix, which is approximated by the last label-bearing span; falling
    # back to "everything after the opener" is the conservative reading and it
    # is what a reader does.
    cm = _CONC_RE.match(c)
    if cm:
        # Two readings of where the subordinate clause ends, and the first
        # comma is not always it: "Since both tests pass the system is linear,
        # as expected." has its comma AFTER the matrix, and splitting there
        # returned "as expected" with the label gone (Reviewer E, R4-F23).
        # Prefer the comma reading; fall back to "after the opener" when the
        # comma reading leaves nothing recognisable behind.  Which is which is
        # decided by the text, not by a rule about commas.
        comma = c.find(",")
        after_comma = c[comma + 1:].strip() if comma != -1 else ""
        after_opener = c[cm.end():].strip()
        c = after_comma or after_opener
        if comma != -1 and len(after_comma.split()) < 3 <= len(after_opener.split()):
            c = after_opener
        if not c:
            return None, "the clause is a bare subordinate clause with no matrix"
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
        close = _closing_comma(c, m.end())
        c = (c[:m.start()] + c[close:]).strip() if close != -1 else c[:m.start()].strip()
    return (c, "") if c else (None, "the clause has no matrix part")


#: Coordinators.  A coordinator joins two **independent** conjuncts, and a
#: hedge in one does not govern a label in the other:
#:
#:     "The result may be described differently under another convention,
#:      but the system is linear."
#:
#: Reviewer B listed these by name in round 2 (§5.1) and three of the four
#: parts of that remedy shipped -- this one did not, so hedge scope went
#: 45 characters -> clause -> segment and still never cut at a coordinator.
#: Six of seven coordination constructions were refused (B, R3-F2), and it is
#: the verbatim `reynolds_number_flow_regime` construction B filed from the
#: real archive as R2-F2.
COORDINATORS = ("but", "yet", "however", "nevertheless", "nonetheless",
                "regardless", "still", "though", "whereas", "instead")

_COORD_RE = re.compile(
    r"(?:^|[,;]\s*|\s)(?:" + "|".join(re.escape(c) for c in COORDINATORS) + r")\s+",
    re.IGNORECASE,
)


def conjuncts(matrix: str) -> list[str]:
    """Split an asserting segment at its coordinators.

    Each conjunct is scoped separately for hedges, so a caveat in one conjunct
    cannot reach a commitment in another. Splitting rather than truncating
    keeps both directions available: the label may be in either conjunct.
    """
    parts, pos = [], 0
    for m in _COORD_RE.finditer(matrix):
        if matrix[pos:m.start()].strip():
            parts.append(matrix[pos:m.start()].strip())
        pos = m.end()
    if matrix[pos:].strip():
        parts.append(matrix[pos:].strip())
    return parts or [matrix]


def _closing_comma(text: str, start: int) -> int:
    """The comma that closes a subordinate interpolation, or -1.

    A subordinate clause may carry commas of its own -- "since the check, a
    nonlinear probe, passed" -- and taking the first one splices a fragment
    into the matrix (Reviewer E, R3-F19). The closing comma is the last one
    before the segment ends or before a coordinator resumes the matrix.
    """
    commas = [m.start() for m in re.finditer(r",", text[start:])]
    if not commas:
        return -1
    stop = len(text)
    cm = _COORD_RE.search(text, start)
    if cm:
        stop = cm.start()
    inside = [start + c for c in commas if start + c < stop]
    return inside[-1] if inside else -1


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

#: A bare epistemic comment has an **anaphoric or absent** subject.  This is an
#: allowlist by design: version 2 used a *determiner blocklist* ("the", "a",
#: "both"...), so a clause whose subject was a bare abstract noun --
#: "Superposition seems to hold" -- read as a bare comment and retracted a
#: commitment it did not modify (Reviewer E, R3-F13).  A blocklist over English
#: noun phrases has an unbounded complement; the pronouns do not.
_ANAPHORIC_SUBJECT = re.compile(
    r"^\W*(?:it|this|that|these|those|i|we|there)\b", re.IGNORECASE)

#: Subjectless openers: a predicate adjective or a bare negation standing alone.
#: "Hard to say", "Not sure", "Difficult to tell", "Unclear".
_SUBJECTLESS = re.compile(
    r"^\W*(?:hard|difficult|impossible|unclear|uncertain|not|no|unsure|"
    r"tough|tricky)\b", re.IGNORECASE)


def is_bare_comment(clause: str, surfaces: Sequence[str]) -> bool:
    """Is this clause a bare epistemic remark rather than a claim of its own?

    ``It seems.`` / ``I am not sure.`` / ``Hard to say, though.`` qualify a
    neighbouring commitment. ``The margin seems comfortable.`` and
    ``Superposition seems to hold.`` do not -- they are claims about something
    else that happen to contain a marker.

    **No length bound.**  Version 2 carried ``_MAX_COMMENT_WORDS = 6``, which
    Reviewer E escaped by one word: ``I am not sure`` was caught and ``I am not
    entirely sure about that`` was not, though both phrases are declared
    (R3-F15). That was a magic constant inside the fix that removed a magic
    constant. The subject test alone does the work, and it does not have a
    threshold to tune.
    """
    from .kinds import _label_hit  # noqa: PLC0415 - circular at import time

    if _label_hit(clause, surfaces):
        return False
    stripped = re.sub(r"^\W*(?:but|and|however|though|although|so)\b", "", clause,
                      flags=re.IGNORECASE).strip()
    return bool(_ANAPHORIC_SUBJECT.match(stripped) or _SUBJECTLESS.match(stripped))


# --------------------------------------------------------------------------
# 5.  The commitment decision
# --------------------------------------------------------------------------


class Commitment:
    """Where a span commits to a label, or why it does not.

    ``notes`` carries hedge markers detected but **not scored** under the
    advisory policy.  They are reported, never decisive.
    """

    __slots__ = ("clause", "index", "reason", "notes")

    def __init__(self, clause: str | None, index: int, reason: str = "",
                 notes: list[str] | None = None):
        self.clause, self.index, self.reason = clause, index, reason
        self.notes = notes or []

    def __bool__(self) -> bool:
        return self.clause is not None


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

    # **A factive clause may supply a commitment but never override one.**
    # `normalize.DISCOURSE_MARKERS` treats a trailing "Note ..." sentence as a
    # qualification and `answer_span` refuses to peel into it; `FACTIVE` then
    # declared the same sentence an assertion, and last-clause-wins handed it
    # the answer -- so an answer of "The system is linear." followed by
    # "Recall that a nonlinear map fails additivity." scored MISMATCH, marking a
    # correct answer WRONG.  Two modules contradicting each other
    # about the same construction (Reviewer E, R3-F16), which re-opened round-1
    # F4.  Discourse position settles it: a factive that FOLLOWS a commitment is
    # elaboration, not retraction.  FACTIVE's census support is 1 in 2,200, so
    # it does not get to outrank anything.
    committing = [(i, c) for i, c in labelled if not _FACTIVE_RE.match(c)]
    if committing:
        labelled = committing

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

        # Governance is the label's own CONJUNCT, not its clause and not its
        # whole segment.  A coordinator joins two independent conjuncts and a
        # caveat in one does not reach a commitment in the other (B, R3-F2).
        owning = next((cj for cj in reversed(conjuncts(bare))
                       if _label_hit(cj, surfaces)), bare)

        marks = _hedges_governing(cls, i, owning, matrix, surfaces)
        # A WITHDRAWAL is not a hedge and stays decisive under either policy:
        # "X. Or is it?" is a retraction of the assertion, not a statement of
        # low confidence in it, so the advisory demotion does not reach it.
        if any(m.startswith("withdrawn") for m in marks):
            last_reason = f"withdrawn ({', '.join(marks)})"
            continue
        if marks and HEDGE_POLICY == "enforce":
            last_reason = f"hedged ({', '.join(marks)})"
            continue
        if marks:
            # Advisory: recorded on the Commitment and surfaced by the caller in
            # `observations`, where it is visible in the results and scores
            # nothing.  See HEDGE_POLICY for why, and D-056.
            return Commitment(owning, i, notes=marks)
        neighbour = _governing_comment(cls, i, surfaces)
        if neighbour:
            last_reason = f"qualified by a bare comment ({neighbour})"
            continue
        return Commitment(owning, i)
    return Commitment(None, -1, last_reason or "no clause commits to a label")


#: Classes a **bare epistemic comment** may carry when it governs a sibling.
#: Modals and verbs of opinion are excluded: ``This should be clear from the
#: square`` is confidence, and ``Imagine a scaled input`` opens a proof.
_COMMENT_CLASSES = frozenset(
    {"copula of appearance", "stated uncertainty", "evidential"})


def _hedges_governing(
    cls: list[str], i: int, owning: str, matrix: str, surfaces: Sequence[str]
) -> list[str]:
    """Which hedges actually govern the commitment in ``owning``?

    Three findings collapse into **one distinction**, and it is the one
    Reviewer B drew in round 2 when it split ``UNCERTAINTY`` into inability and
    imprecision: *whose* confidence is being qualified.

    ``I am not 100% sure, but the system is linear``   -- the speaker's, so it
    governs the whole utterance across the coordinator (B, round-1 F1).
    ``The distinction is not entirely obvious, but the system is linear``
                                                       -- the *reasoning's*, so
    it stays in its own conjunct and the commitment stands (B, round-2 §5).
    ``The system is linear (the wrong answer would be nonlinear)``  -- a modal
    about a hypothetical answer, local, and stripped with its parenthetical
    (E, R2-F12).
    ``The system is linear (or nonlinear, I am not sure)``  -- the speaker's
    again, inside a parenthetical, and it must NOT be stripped away with it
    (E, R3-F14).

    The discriminator is not the class and not the position: it is whether the
    hedge sits in a **bare epistemic comment** -- a unit with an anaphoric or
    absent subject, which is therefore about the answer rather than about some
    other proposition. A bare comment governs from anywhere; anything else is
    local to the label's own conjunct.
    """
    from .kinds import _label_hit  # noqa: PLC0415 - circular at import time

    # 1. Local hedges, in the label's own conjunct, with label-naming
    #    parentheticals removed so a gloss about the *other* label cannot hedge
    #    this one.
    marks = list(hedge_markers(owning))

    # 2. Bare epistemic comments anywhere: sibling conjuncts, parenthetical
    #    asides, and other clauses of the span.
    units: list[tuple[str, bool]] = []          # (text, follows_the_commitment)
    for cj in conjuncts(matrix):
        if cj.strip() != owning.strip():
            units.append((cj, matrix.find(cj) > matrix.find(owning)))
    for m in _PAREN_RE.finditer(matrix):
        # A parenthetical often packs a label and a comment together --
        # "(or nonlinear, I am not sure)" -- and the label makes the whole
        # group fail the bare-comment test, hiding the comment inside it
        # (Reviewer E, R3-F14). Split on commas so each piece is judged alone.
        for piece in m.group(1).split(","):
            units.append((piece, True))
    for j, c in enumerate(cls):
        if j != i:
            units.append((c, j > i))

    for text, follows in units:
        # A label-free QUESTION following the answer withdraws it.  "The system
        # is linear. Or is it?" is not a commitment, and nothing else here sees
        # that: `segment` only rejects a question that carries the label itself.
        # Found by the negative controls Reviewer B said the recall corpus
        # needed (R3-F1) -- 30 of 39 answers were still credited under a
        # trailing "Or is it?".
        # The `not _label_hit` guard this rule shipped with meant "Or is it?"
        # was caught and "Or is it nonlinear?" was **credited** -- the *more
        # explicit* withdrawal passing while the vaguer one was refused
        # (Reviewer E, R4-F26).  A question is not an assertion whether or not
        # it names a label, so the guard is gone.
        if follows and text.rstrip().endswith("?"):
            marks.append("withdrawn by a following question")
            continue
        if not is_bare_comment(text, surfaces):
            continue
        classes = _COMMENT_CLASSES if follows else frozenset({"stated uncertainty"})
        marks.extend(hedge_markers(text, classes))
    return sorted(set(marks))


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


def suspends_what_follows(lead: str) -> str:
    """Does this preface suspend the answer that follows it? Reason, or ``""``.

    Only the framings that leave nothing asserted count: a hypothetical, a
    restatement of the task, a question, a bare qualifier. A **concessive or
    factive** preface does *not* suspend -- its matrix clause is precisely the
    answer that follows, which is why ``Although the algebra is fiddly, a) Yes,
    b) Yes`` must still be credited while ``If the additivity test holds, a)
    Yes, b) Yes`` must not.

    This exists for the enumerated label-tuple path, where the answer's slots
    are located positionally *after* the enumerators and any preface is
    otherwise dropped unread (Reviewer B, R3-F1).
    """
    c = lead.strip()
    if not c:
        return ""
    # **A preface closed by a sentence stop does not scope over what follows**,
    # and `?` closes one just as `.` does.  Refusing every `?`-final preface
    # meant "Is the system memoryless and causal? a) No, b) No" was UNRESOLVED
    # while the categorical path credited the identical construction -- two
    # paths disagreeing about one sentence, which is R3-F16's shape reproduced
    # inside the function written to fix it (Reviewer E, R4-F27).  A trailing
    # question is a withdrawal and is handled where it belongs, in
    # `_hedges_governing`.
    # "Consider a scaled input a*x[n] and a shifted input x[n-k]. a) Yes, b) Yes"
    # sets up a derivation and then answers; "If the additivity test holds,
    # a) Yes, b) Yes" is one conditional sentence.  Without this, `consider`,
    # `imagine` and `let` -- which are in HYPOTHETICAL because "Assume the
    # system is linear" must not commit -- also suspended every answer that
    # merely followed a worked setup, and Reviewer B's R2-F1.3 says those open
    # the standard homogeneity proof.  The punctuation is what separates the
    # two, and unlike R3-F17's comma it is not optional: a new sentence needs it.
    if c.rstrip().endswith((".", "!")):
        return ""
    if _TASK_RE.match(c):
        return f"the preface restates the task ({_TASK_RE.match(c).group(0)!r})"
    if _BARE_QUAL_RE.match(c):
        return "the preface is a bare qualifier"
    m = _HYPO_RE.match(c)
    if m:
        return (f"the answer is conditional on a hypothesis "
                f"({m.group(0).strip()!r}), not asserted")
    # A speaker-level hedge in the preface governs the answer that follows it,
    # for the same reason it governs across a coordinator: "I am not sure, but
    # perhaps a) Yes, b) Yes" qualifies the answer, not some other proposition.
    # Without this the enumerated answers took a hedged preface and dropped it
    # unread, which the negative controls caught on 16 of 39.
    if HEDGE_POLICY == "enforce":
        for piece in re.split(r",|\bbut\b|\band\b", c):
            if is_bare_comment(piece, ()):
                marks = hedge_markers(piece, _COMMENT_CLASSES)
                if marks:
                    return f"the preface hedges the answer ({', '.join(marks)})"
    return ""


def is_assertion(clause: str) -> tuple[bool, str]:
    """Back-compatible wrapper: does this clause assert anything at all?"""
    matrix, why = segment(clause)
    return (matrix is not None), why
