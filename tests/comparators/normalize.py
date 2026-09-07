"""Surface normalisation, and the observed vocabulary that drives it (D4.2).

Every rule in this file is here because it was **observed in archived model
output**, not because it was imagined.  Each carries a ``support`` annotation
naming the templates and trace indices that motivate it, and
``derive_vocabulary.py`` recomputes those counts from
``error_analysis_annotation/samples/`` so the annotations cannot drift away from
the archive the way a transcribed reference value can (D-034).

**Coverage is stated per template, never averaged** (D-048): the archive holds
16 / 24 / 15 / 6 traces for the four Phase 4 templates, so the symbolic rules in
particular rest on six traces and resolve nothing rarer than ~1-in-2.  See
``phase4_vocabulary.md`` |S|5.
"""

from __future__ import annotations

import re
import unicodedata

# --------------------------------------------------------------------------
# 1.  Character-level normalisation
# --------------------------------------------------------------------------

#: Characters models emit where the gold uses ASCII.  Every entry observed.
CHAR_MAP = {
    "\u2212": "-",   # MINUS SIGN      claude-opus-4-7 continuity 0,1; gpt-5 mc 22
    "\u2013": "-",   # EN DASH         gemma signal 1
    "\u2014": "-",   # EM DASH         claude mc 0
    "\u00b7": "*",   # MIDDLE DOT      gemini linearity 6; qwen linearity 12
    "\u00d7": "*",   # MULTIPLICATION  gemini mc 9
    "\u2264": "<=",  # LESS-THAN OR EQ gpt-5-mini mc 16,17; gemini mc 9,11
    "\u2265": ">=",
    "\u2260": "!=",  # NOT EQUAL       deepseek linearity 2,3,4; gemini linearity 6
    "\u2192": "->",  # RIGHTWARDS      claude mc 0
    "\u2248": "~=",
    "\u00a0": " ",   # NBSP
    "\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"',
}

#: Unicode superscripts -> caret form.  Observed on every claude/gemini
#: continuity trace (x^2 written x\u00b2) -- 4 of the 6 symbolic traces.
SUPERSCRIPT_MAP = {
    "\u00b2": "^2", "\u00b3": "^3", "\u00b9": "^1", "\u2070": "^0",
    "\u2074": "^4", "\u2075": "^5", "\u2076": "^6", "\u2077": "^7",
    "\u2078": "^8", "\u2079": "^9",
}


def normalize_chars(text: str) -> str:
    """Fold the observed Unicode habits onto ASCII.

    Superscripts are folded *before* NFKC, because NFKC maps U+00B2 to a bare
    ``2`` -- turning ``x\u00b2`` into ``x2`` and silently destroying the
    exponent.  That is a false *reject* the pipeline would never explain, and it
    is why this function does not simply call ``unicodedata.normalize``.
    """
    for src, dst in SUPERSCRIPT_MAP.items():
        text = text.replace(src, dst)
    for src, dst in CHAR_MAP.items():
        text = text.replace(src, dst)
    # NFKC after the two explicit passes: it still folds ligatures and
    # full-width forms, and can no longer reach a superscript.
    return unicodedata.normalize("NFKC", text)


# --------------------------------------------------------------------------
# 2.  LaTeX stripping
# --------------------------------------------------------------------------

#: LaTeX wrappers observed.  qwen2.5-72b emits them on 5 of its 5 signal traces
#: and 3 of its 3 linearity traces; mathstral on its one continuity trace;
#: deepseek throughout.  A comparator that does not strip these scores
#: qwen and deepseek as if they had answered nothing.
_LATEX_MATH_DELIMS = [
    (r"\\\[", r"\\\]"),   # \[ ... \]
    (r"\\\(", r"\\\)"),   # \( ... \)
]

_FRAC_RE = re.compile(r"\\d?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}")
_BOXED_RE = re.compile(r"\\boxed\s*\{(.*)\}", re.DOTALL)
_TEXT_RE = re.compile(r"\\(?:text|mathrm|mathbf|operatorname)\s*\{([^{}]*)\}")
_PARTIAL_RE = re.compile(r"\\partial\s*")
_CMD_RE = re.compile(r"\\(?:left|right|displaystyle|quad|qquad|,|;|:|!)\s*")


def strip_latex(text: str) -> str:
    """Remove the LaTeX scaffolding models wrap answers in.

    Order matters: ``\\boxed`` and ``\\frac`` are unwrapped before the escaped
    braces are, because both are recognised by their brace structure.
    """
    m = _BOXED_RE.search(text)
    if m:
        text = m.group(1)
    text = _TEXT_RE.sub(r"\1", text)
    # \frac{a}{b} -> (a)/(b); repeat for nesting depth met in practice.
    for _ in range(4):
        new = _FRAC_RE.sub(r"((\1)/(\2))", text)
        if new == text:
            break
        text = new
    text = _PARTIAL_RE.sub("d", text)
    text = re.sub(r"\\cdot\s*", "*", text)
    text = re.sub(r"\\times\s*", "*", text)
    text = re.sub(r"\\neq\s*", "!=", text)
    text = re.sub(r"\\leq?\s*", "<=", text)
    text = re.sub(r"\\geq?\s*", ">=", text)
    text = _CMD_RE.sub(" ", text)
    for open_d, close_d in _LATEX_MATH_DELIMS:
        text = re.sub(open_d, " ", text)
        text = re.sub(close_d, " ", text)
    text = text.replace("$$", " ").replace("$", " ")
    # Escaped braces \{ \} -> { }.  qwen writes every sequence this way.
    text = text.replace(r"\{", "{").replace(r"\}", "}")
    # ``x^{2}`` -> ``x^2``, after the braces are unescaped.  Without this the
    # exponent survives as a brace group, the polynomial parser reads ``x`` to
    # the first power, and a correctly typeset answer is rejected -- a false
    # reject found by D4.4 case sym-04, not by the archive.
    text = re.sub(r"\^\s*\{\s*(-?\d+)\s*\}", r"^\1", text)
    text = re.sub(r"_\s*\{\s*([^{}]*)\s*\}", r"_\1", text)
    text = re.sub(r"\\([_^%&#])", r"\1", text)
    # A trailing bare backslash-word we did not recognise is dropped rather
    # than left to poison a token; it is recorded by the caller as an
    # observation when it occurs inside an answer span.
    return text


# --------------------------------------------------------------------------
# 3.  Answer-span extraction
# --------------------------------------------------------------------------

#: Ordered by priority.  Every marker observed in the 61 Phase 4 traces.
#: Counts are in :data:`SUPPORT` and are checked against the archive by
#: ``derive_vocabulary.py`` -- not written in prose beside the rule, which is
#: how the first draft came to claim 44/61 and 47/61 for the two markers below
#: when the measured counts are 57 and 55.  Numbers written from impression next
#: to the thing they describe are D-034's failure mode in miniature, and this
#: file argues for D-034, so it does not get to commit it.
ANSWER_MARKERS = [
    r"##\s*Final\s+Answer\s*",      # the dominant house style
    r"\*\*Answer:?\*\*\s*",         # usually nested inside the heading above
    r"Final\s+Answer:?\s*",         # wizardmath signal 14
    r"The\s+answer\s+is:?\s*",      # wizardmath signal 14,15; continuity 4
    r"Answer:?\s*",
]

#: Support for every normalisation rule, as
#: ``rule -> (count in the 61 Phase 4 traces, count in all 2,200)``.
#:
#: **Machine-checked, not transcribed.**  ``derive_vocabulary.py`` recomputes
#: every pair from the archive on disk and fails on any disagreement, so this
#: table cannot drift away from the evidence.  Regenerate it from that module's
#: output rather than editing it by hand.
#:
#: Read the two columns together and the second is not reassurance.  Rates are
#: **item-driven**, not model-driven (measured: 15 of 27 rules), so a large
#: whole-archive count does NOT license a coverage claim about a template with
#: six traces.  What the second column licenses is narrower and is what the
#: rules actually need: that the form exists and that its meaning is fixed.
SUPPORT: dict[str, tuple[int, int]] = {
    "unicode-minus": (4, 71),
    "unicode-superscript": (3, 236),
    "unicode-operator": (0, 330),
    "latex-inline-math": (5, 332),
    "latex-display-math": (1, 66),
    "latex-escaped-brace": (5, 6),
    "latex-frac": (1, 53),
    "latex-boxed": (1, 104),
    "latex-text": (0, 274),
    "latex-exponent-brace": (0, 165),
    "markdown-bold": (3, 239),
    "trailing-annotation": (16, 86),
    "marker-final-answer-heading": (57, 1824),
    "marker-bold-answer": (55, 1736),
    "marker-the-answer-is": (4, 202),
    "negator-not": (25, 34),
    "negator-fused-non": (3, 3),
    "negator-neither-nor": (2, 2),
    "negator-fails": (3, 3),
    "hedge": (0, 1),
    "brace-sequence": (17, 447),
    "bracket-sequence": (1, 64),
    "origin-asterisk": (1, 19),
    "explicit-index-list": (2, 6),
    "per-element-assignment": (2, 2),
    "arbitrary-function": (2, 30),
    "thousands-separator": (0, 32),
}

#: Rules with **zero** Phase 4 support.  Each is carried on the strength of the
#: wider archive alone, and each is therefore a place where the comparator is
#: specified against forms these four templates have never been observed to
#: produce.  Listed rather than buried: this is the concrete content of the
#: residual risk D4.2 carries.
NO_PHASE4_SUPPORT = tuple(k for k, (p4, _) in SUPPORT.items() if p4 == 0)

_MARKER_RES = [re.compile(m, re.IGNORECASE) for m in ANSWER_MARKERS]


def answer_span(text: str) -> tuple[str, str]:
    """Return ``(span, marker)`` -- the region of ``text`` that states the answer.

    The **last** occurrence of the highest-priority marker present wins, because
    models restate the answer and the final statement is the committed one.

    Extracting a span before normalising is not a convenience.  Reviewer-visible
    case: linearity trace 0 argues the system is ``incrementally linear``
    mid-derivation and answers ``NOT LINEAR``.  A keyword search over the whole
    trace reads the first and scores it a match -- a false accept on a
    *correct* answer, which is the worst kind, because it is invisible.
    """
    marker = ""
    # Markers nest: a "## Final Answer" heading followed by an "**Answer:**"
    # line is the house style, and taking only the outer one leaves "Answer:"
    # glued to the answer, where the symbolic comparator reads it as a symbol.
    # Peel until stable.
    #
    # **But never peel across a discourse marker.**  Reviewer E (F4): the bare
    # `Answer:?` marker is low-priority and matches inside a post-hoc remark, so
    # "**Answer:** The system is linear.\n\nNote: Answer: not linear if the
    # offset were nonzero." peeled into the caveat and scored *that* -- a false
    # accept whose sign is arbitrary.  Peel-to-last is right for a restatement
    # and wrong for a qualification, and the marker alone cannot tell them
    # apart; the opener of the sentence carrying the marker can.
    for _ in range(4):
        for rx in _MARKER_RES:
            hits = [m for m in rx.finditer(text) if not _is_qualifying(text, m.start())]
            if hits:
                last = hits[-1]
                text, marker = text[last.end():].strip(), last.group(0).strip()
                break
        else:
            break
    return text.strip(), marker


#: Sentence openers that mark what follows as a qualification rather than the
#: committed answer.  A marker inside such a sentence is not an answer marker.
DISCOURSE_MARKERS = (
    "note", "notes", "n.b.", "nb", "however", "unless", "if", "caveat",
    "caution", "aside", "disclaimer", "but", "although", "though", "except",
    "assuming", "provided", "otherwise", "alternatively", "in contrast",
)

_QUALIFIER_RE = re.compile(
    r"(?:^|[.;!?\n])\s*(?:\*+\s*)?(?:" + "|".join(re.escape(d) for d in DISCOURSE_MARKERS)
    + r")\b[^.;!?\n]{0,80}$",
    re.IGNORECASE,
)


def _is_qualifying(text: str, at: int) -> bool:
    """Does the marker at ``at`` sit inside a qualifying sentence?"""
    return bool(_QUALIFIER_RE.search(text[:at]))


#: Trailing annotations models append to an answer.  Observed on 8/24
#: memory-causality traces (gemini, gpt-5-mini) and 3/15 linearity traces.
#: These are commentary, never part of the answer.
_TRAILING_ANNOT_RE = re.compile(
    r"[\[(]\s*(?:units?\s*[:=]?\s*)?"
    r"(?:N/?A|none|unitless|dimensionless|no\s+units?|"
    r"system\s+properties(?:\s*/\s*unitless)?|"
    r"units?\s+of\s+velocity[^)\]]*)\s*[\])]",
    re.IGNORECASE,
)


def strip_trailing_annotations(text: str) -> str:
    return _TRAILING_ANNOT_RE.sub(" ", text)


_MD_BOLD_RE = re.compile(r"\*\*(.+?)\*\*", re.DOTALL)
_MD_ITAL_RE = re.compile(r"(?<![*\w])\*(?!\s)([^*\n]+?)(?<!\s)\*(?![*\w])")


def strip_markdown(text: str, keep_asterisks: bool = False) -> str:
    """Remove markdown emphasis.

    ``keep_asterisks`` exists because ``*`` is **not** decoration in a
    sequence answer -- it is the origin marker, and ``{2, 3, *-8*, 6}`` means
    something different from ``{2, 3, -8, 6}``.  The sequence comparator calls
    this with ``keep_asterisks=True`` and does its own emphasis handling; every
    other kind strips.  Getting this backwards deletes the answer.
    """
    text = _MD_BOLD_RE.sub(r"\1", text)
    if not keep_asterisks:
        text = _MD_ITAL_RE.sub(r"\1", text)
    text = re.sub(r"^#{1,6}\s*", "", text, flags=re.MULTILINE)
    text = re.sub(r"`+", "", text)
    return text


def prepare(text: str, keep_asterisks: bool = False) -> str:
    """The common front end: chars, LaTeX, markdown, trailing annotations."""
    text = normalize_chars(text)
    text = strip_latex(text)
    text = strip_markdown(text, keep_asterisks=keep_asterisks)
    text = strip_trailing_annotations(text)
    return re.sub(r"[ \t]+", " ", text).strip()


# --------------------------------------------------------------------------
# 4.  Categorical label vocabulary (D4.2 proper)
# --------------------------------------------------------------------------

#: Tokens that negate the label that follows them, with observed support.
#:   "not"        linearity 1,3,4,10,11,12,13,14; mc 3,4,5,6,7,9,10,11,12,13,...
#:   "NOT"        claude mc 0,1; claude linearity 0   (bold uppercase)
#:   "non"        gpt-5-mini mc 18 ("non-causal"); gemini linearity 6
#:   "nonlinear"  deepseek linearity 2,4; gpt-5 linearity 9  (fused, no space)
#:   "neither"    gemma mc 15; meta-llama mc 23        ("neither A nor B")
#:   "fails"      qwen linearity 12,13,14; claude linearity 0
#: ``"no"`` is deliberately NOT here.  As a bare interjection it precedes a
#: clause it does not negate -- "The system is nonlinear -- no, wait, it is
#: linear" -- and treating it as a negator flipped that committed *wrong*
#: answer into the right one, a false accept Reviewer B found (F4).  The
#: constructions that need it ("no memory") are declared as label surfaces in
#: their own right, which is the same declaration-beats-inference rule the
#: fused negatives follow.
NEGATORS = ("not", "non", "never", "fails", "fail", "violates", "isn't", "doesn't")

#: "neither X nor Y" negates *both* labels.  Two of 24 memory-causality traces
#: use it and it is the single most dangerous construction in the archive: a
#: matcher that looks for the substrings "memoryless" and "causal" without
#: scoping the negation reads "neither memoryless nor causal" as (Yes, Yes) --
#: the exact inversion of the answer, scored as a match.
#: Scope is **structural, not a character count**.  The first version allowed at
#: most 60 characters between ``neither`` and ``nor``; Reviewer E (F5) noted
#: that both archived instances have a 12-character gap, so the bound was
#: untested by a factor of five -- and showed that one ordinary interposed
#: clause escapes it, at which point *both* labels read positive and the answer
#: is inverted rather than merely missed.
#:
#: The correct scope is the sentence: ``neither`` binds to the next ``nor``
#: before a sentence boundary, however much intervenes.  A character budget was
#: standing in for that and could only ever be wrong in one direction or the
#: other.
NEITHER_NOR_RE = re.compile(
    r"\bneither\b(?P<first>[^.;!?\n]*?)\bnor\b(?P<second>[^.;!?\n]*)",
    re.IGNORECASE,
)

#: Hedges.  A hedged answer has not committed, so it is UNRESOLVED -- never a
#: match, and never silently a mismatch either.  This is the specific case
#: spec |S|4.4 names for Reviewer B.  None of the 61 archived Phase 4 answers
#: hedges, which is itself a measurement: **this list is not derived from the
#: Phase 4 archive** and is flagged as such in phase4_vocabulary.md |S|4.
#: It is drawn from the wider 2,200-trace archive (route 1) and validated there.
HEDGES = (
    "appears to be", "appears", "seems to be", "seems", "would appear",
    "likely", "probably", "possibly", "may be", "might be", "could be",
    "arguably", "tends to be", "i think", "i believe", "presumably",
    "in some sense", "cannot be determined", "cannot determine",
    "unclear", "ambiguous", "hard to say", "not entirely clear",
    "one could argue", "depends on how", "depends on the interpretation",
)

#: **Two words are deliberately absent from the list above and it matters.**
#:
#: ``it depends`` and ``either`` read as hedges in ordinary prose and are
#: *domain vocabulary* here.  ``a) Not memoryless (it depends on the past input
#: x[n-2])`` is a fully committed, correct answer -- "depends on" is the literal
#: definition of memory in this item -- and ``the system fails at least one of
#: the tests (additivity or homogeneity)`` is gold's own phrasing.  With both
#: words in the list, the comparator returned UNRESOLVED on a correct archived
#: answer.  A hedge vocabulary built by asking "what sounds non-committal?"
#: rather than by testing it against the corpus produces exactly this, and it
#: produces it as a *false reject*, which is invisible in an accuracy number.
_HEDGE_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(h) for h in HEDGES) + r")\b", re.IGNORECASE
)

#: How far before a label a hedge has to be to govern it.  Same window as
#: negation, and for the same reason: scope is local.
HEDGE_WINDOW = 45


def has_hedge(text: str, at: int | None = None) -> bool:
    """Is the answer hedged?

    With ``at``, only the window immediately before the label at that offset is
    examined -- a hedge three sentences earlier in a derivation does not make a
    committed final answer non-committal.  Without ``at``, the whole string is
    examined, which is the right test only when no label was resolved.
    """
    if at is None:
        return bool(_HEDGE_RE.search(text))
    return bool(_HEDGE_RE.search(text[max(0, at - HEDGE_WINDOW):at]))


def hedges_in(text: str, at: int | None = None) -> list[str]:
    window = text if at is None else text[max(0, at - HEDGE_WINDOW):at]
    return sorted({m.group(0).lower() for m in _HEDGE_RE.finditer(window)})
