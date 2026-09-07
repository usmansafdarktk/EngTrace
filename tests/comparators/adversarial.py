"""D4.4 -- the adversarial set: hand-written near misses, >=20 per comparator.

Two populations, deliberately balanced:

* **correct-but-phrased-differently** -- the right answer in a form the
  comparator might not recognise.  These hunt **false rejects**.
* **wrong-but-similar** -- an answer that looks like gold and is not.  These
  hunt **false accepts**, which is the failure mode Decision 3 says matters
  more (D-049).

Why this exists even though the archive scores 100%: **the archive is where the
vocabulary came from.**  Scoring a comparator on the corpus its rules were
derived from measures memorisation, not generalisation -- the same shape as
Phase 3's gold corpus passing 80/80 against a verifier that hardcoded one
template's symbol names, because gold traces do not lie.  Every case below was
written to break something specific, and the ones that *did* break it are
marked ``caught:`` in their note.

Run ``python -m tests.comparators.adversarial`` to regenerate ``adversarial.json``.
"""

from __future__ import annotations

import json
import os

# --------------------------------------------------------------------------
# Gold answers used as fixtures.  These are real gold strings in the shape the
# templates emit, so a case exercises the answer-span extractor too.
# --------------------------------------------------------------------------

G_LIN_NO = "**Answer:**\nThe system fails at least one of the tests (additivity or homogeneity), therefore, the system is **not linear**."
G_LIN_YES = "**Answer:**\nThe system satisfies both additivity and homogeneity, therefore, the system is **linear**."

G_MC_NY = "**Answer:**\na) Memoryless: **No**\nb) Causal: **Yes**"
G_MC_NN = "**Answer:**\na) Memoryless: **No**\nb) Causal: **No**"
G_MC_YY = "**Answer:**\na) Memoryless: **Yes**\nb) Causal: **Yes**"

G_SEQ = "**Answer:**\nThe resulting sequence is z[n] = {-9, 6, *-5*, -4, 5, -3, 3}"
G_SEQ_NOORIG = "**Answer:**\nThe resulting sequence is y[n] = {-1, 5, -4, 3}"

G_SYM = "**Answer:**\nThe simplest expression for the x-component of velocity is u = -2x^2 + 6xy."
G_SYM_FRAC = "**Answer:**\nThe simplest expression for the x-component of velocity is u = 1.5x^2 - 6xy."

G_NUM = "**Answer:**\nThe required reactor volume is 7.65 L."
G_NUM_COMMA = "**Answer:**\nThe total annual demand is **4,921** units."

G_CHECK = "**Answer:**\nThe maximum deflection is 18.4 mm, which is less than the 25 mm limit, so the design is acceptable."

G_NARR = "**Answer:**\nThe pressure drop rises because the transition to turbulent flow increases wall shear."


def _c(cid, kind, gold, candidate, correct, note, template=None, options=None):
    return {
        "id": cid, "kind": kind, "gold": gold, "candidate": candidate,
        "correct": correct, "note": note,
        "template": template, "options": options or {},
    }


# ==========================================================================
# categorical -- 24 cases
# ==========================================================================

CATEGORICAL = [
    # -- correct, phrased differently ------------------------------------
    _c("cat-01", "categorical", G_LIN_NO, "**Answer:** Nonlinear", True,
       "fused negative, no space"),
    _c("cat-02", "categorical", G_LIN_NO, "**Answer:** non-linear", True,
       "fused negative, hyphenated"),
    _c("cat-03", "categorical", G_LIN_NO, "**Answer:** NOT LINEAR", True,
       "uppercase negated positive"),
    _c("cat-04", "categorical", G_LIN_NO, "## Final Answer\n**Answer:** The system is not linear. [N/A]", True,
       "trailing unit annotation"),
    _c("cat-05", "categorical", G_LIN_NO, "**Answer:** The system fails the homogeneity test and is therefore not linear.", True,
       "negator separated from the label by nine words"),
    _c("cat-06", "categorical", G_LIN_NO, r"**Answer:** \(\text{Non-linear}\)", True,
       "LaTeX \\text wrapper"),
    _c("cat-07", "categorical", G_LIN_YES, "**Answer:** Linear [dimensionless]", True,
       "positive with annotation"),
    _c("cat-08", "categorical", G_LIN_YES, "**Answer:** The system satisfies superposition, so it is linear.", True,
       "positive reached through a synonym clause"),
    _c("cat-09", "categorical", G_LIN_NO,
       "The system is linear in its first term.\n## Final Answer\n**Answer:** Not linear", True,
       "distractor 'linear' earlier in the trace; answer span must win"),
    _c("cat-10", "categorical", G_LIN_NO,
       "**Answer:** This is an incrementally linear (affine) system, so it is **not linear**.", True,
       "'incrementally linear' distractor inside the answer span itself"),

    # -- wrong, but similar ----------------------------------------------
    _c("cat-11", "categorical", G_LIN_NO, "**Answer:** Linear", False,
       "the plain opposite"),
    _c("cat-12", "categorical", G_LIN_YES, "**Answer:** Nonlinear", False,
       "the plain opposite, fused"),
    _c("cat-13", "categorical", G_LIN_NO,
       "**Answer:** The system is linear because it fails the additivity test.", False,
       "wrong label with a negator nearby -- the negator must not rescue it"),
    _c("cat-14", "categorical", G_LIN_YES,
       "**Answer:** The system is not non-linear.", True,
       "double negative that resolves to linear"),
    _c("cat-15", "categorical", G_LIN_NO,
       "**Answer:** The system appears to be linear.", False,
       "hedged AND wrong -- must not be a MATCH"),
    # LABEL CORRECTED.  First written as correct=True on the reading "the model
    # named the right label".  That is the wrong reading: spec S4.4 names a
    # hedge that never commits as the specific thing a pedagogy reviewer must
    # not see accepted, so the creditable answer to a hedged response is *no
    # credit*, and UNRESOLVED is the right outcome rather than a false reject.
    # Recorded because the mislabel made the comparator look wrong when the
    # label was (phase4_summary.md S8, error 2).
    _c("cat-16", "categorical", G_LIN_YES,
       "**Answer:** The system appears to be linear.", False,
       "hedged and naming the right label.  Under `enforce` this is UNRESOLVED, "
       "which is what P6 asks for; under the shipped ADVISORY default (D-056) it "
       "is a MATCH and therefore a **declared false accept** -- the visible cost "
       "of demoting a layer that fires on 0 of 2,200 real spans.  Kept in the "
       "corpus, like num-16b, so the gate carries the trade rather than hiding it"),
    _c("cat-17", "categorical", G_LIN_YES,
       "**Answer:** It could be linear or nonlinear depending on how you define the input.", False,
       "names both labels, commits to neither"),
    _c("cat-18", "categorical", G_LIN_NO, "**Answer:** The system is time-invariant.", False,
       "answers a different question entirely"),
    _c("cat-19", "categorical", G_LIN_NO, "**Answer:** 42", False,
       "no label at all"),
    _c("cat-20", "categorical", G_LIN_NO, "**Answer:**", False,
       "empty answer span"),
    _c("cat-21", "categorical", G_LIN_YES,
       "**Answer:** Linear\n\nWait -- on reflection the offset makes it nonlinear.", False,
       "self-correction after the answer; the last statement is the committed one"),
    _c("cat-22", "categorical", G_LIN_NO,
       "**Answer:** The system is not linear.\n\nActually it is linear.", False,
       "self-correction that lands on the wrong label"),
    _c("cat-23", "categorical", G_LIN_NO,
       "**Answer:** nonlinearity is present, so linear is the wrong description", True,
       "'nonlinearity' as a noun, plus 'linear' in a negated clause"),
    _c("cat-24", "categorical", G_LIN_YES,
       "**Answer:** LINEAR (satisfies both additivity and homogeneity)", True,
       "uppercase positive with a parenthetical"),
]

# ==========================================================================
# categorical[tuple] -- the label-tuple composite, 22 cases
# ==========================================================================

TUPLE = [
    _c("tup-01", "categorical[tuple]", G_MC_NY, "**Answer:** Memoryless: No; Causal: Yes", True,
       "gold's own form"),
    _c("tup-02", "categorical[tuple]", G_MC_NY, "**Answer:** a) Not memoryless; b) Causal", True,
       "property-word encoding"),
    _c("tup-03", "categorical[tuple]", G_MC_NY, "**Answer:** a) No, b) Yes", True,
       "bare yes/no, positional"),
    _c("tup-04", "categorical[tuple]", G_MC_NY,
       "**Answer:** The system has memory but is causal.", True,
       "declared negative surface plus a 'but' clause boundary"),
    _c("tup-05", "categorical[tuple]", G_MC_NN,
       "**Answer:** The system is neither memoryless nor causal.", True,
       "neither/nor negates both slots"),
    _c("tup-06", "categorical[tuple]", G_MC_NN,
       "**Answer:** a) Not memoryless, b) Non-causal [unitless]", True,
       "fused negative on the second slot, with annotation"),
    _c("tup-07", "categorical[tuple]", G_MC_YY, "**Answer:** a) Yes, b) Yes", True,
       "both positive, unlabelled"),
    _c("tup-08", "categorical[tuple]", G_MC_YY,
       "**Answer:** The system is memoryless and therefore causal.", True,
       "both positive, one clause"),
    _c("tup-09", "categorical[tuple]", G_MC_NY,
       "**Answer:** a) Not memoryless (it depends on x[n-2]) b) Causal (no future inputs)", True,
       "'depends' is domain vocabulary, not a hedge"),
    _c("tup-10", "categorical[tuple]", G_MC_NY,
       "## Final Answer\n**Answer:** Memoryless: No\nCausal: Yes", True,
       "newline-separated slots"),

    # -- wrong ------------------------------------------------------------
    _c("tup-11", "categorical[tuple]", G_MC_NY, "**Answer:** Memoryless: No; Causal: No", False,
       "second slot wrong"),
    _c("tup-12", "categorical[tuple]", G_MC_NY, "**Answer:** Memoryless: Yes; Causal: Yes", False,
       "first slot wrong"),
    _c("tup-13", "categorical[tuple]", G_MC_NN,
       "**Answer:** The system is memoryless and causal.", False,
       "both slots inverted -- the neither/nor trap run backwards"),
    _c("tup-14", "categorical[tuple]", G_MC_NN,
       "**Answer:** a) Yes, b) Yes", False,
       "positional yes/yes against a No/No gold"),
    # LABEL CORRECTED.  First written as correct=False on the reading "the
    # enumerators are swapped, so the answer is misordered".  But the model
    # *names* each property beside its value, so the content identifies the
    # slots and the enumerator is presentation.  The comparator resolves slots
    # by content when content is available and by position only when it is
    # not, which is the right rule; the label was wrong
    # (phase4_summary.md S8, error 3).
    _c("tup-15", "categorical[tuple]", G_MC_NY,
       "**Answer:** a) Causal, b) Not memoryless", True,
       "enumerators swapped but each slot is named beside its value; content resolves it"),
    _c("tup-16", "categorical[tuple]", G_MC_NY,
       "**Answer:** The system is not memoryless.", False,
       "only one slot answered"),
    _c("tup-17", "categorical[tuple]", G_MC_NY,
       "**Answer:** It seems to be causal; memory is unclear.", False,
       "hedged on both slots"),
    _c("tup-18", "categorical[tuple]", G_MC_NY,
       "**Answer:** a) The system is not memoryless. b) The system might be causal.", False,
       "first slot right, second hedged -- composite must not pass on a partial"),
    _c("tup-19", "categorical[tuple]", G_MC_YY,
       "**Answer:** a) Memoryless, b) Not causal", False,
       "second slot wrong; also physically impossible, which the comparator does not know and should not need to"),
    _c("tup-20", "categorical[tuple]", G_MC_NY,
       "**Answer:** a) No memory, b) Causal", False,
       "'no memory' means memoryless=Yes; gold says No"),
    _c("tup-21", "categorical[tuple]", G_MC_NY,
       "**Answer:** Not memoryless. Causal.", True,
       "full stops as the clause boundary"),
    _c("tup-22", "categorical[tuple]", G_MC_NN,
       "**Answer:** a) has memory b) anticausal", True,
       "both declared negatives, no negators at all"),
]

# ==========================================================================
# sequence -- 24 cases
# ==========================================================================

SEQ = [
    _c("seq-01", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, *-5*, -4, 5, -3, 3}", True,
       "identical"),
    _c("seq-02", "sequence", G_SEQ, r"**Answer:** \(z[n] = \{-9, 6, *-5*, -4, 5, -3, 3\}\)", True,
       "LaTeX-escaped braces"),
    _c("seq-03", "sequence", G_SEQ, "**Answer:** z[n] = [-9, 6, *-5*, -4, 5, -3, 3]", True,
       "square brackets instead of braces"),
    _c("seq-04", "sequence", G_SEQ,
       "**Answer:** z[n] = {-9, 6, -5, -4, 5, -3, 3} for n = -2, -1, 0, 1, 2, 3, 4", True,
       "origin given by an index list rather than an asterisk"),
    _c("seq-05", "sequence", G_SEQ,
       "**Answer:** z[-2] = -9, z[-1] = 6, z[0] = -5, z[1] = -4, z[2] = 5, z[3] = -3, z[4] = 3", True,
       "per-element assignments"),
    _c("seq-06", "sequence", G_SEQ,
       "**Answer:** z[n] = {0, -9, 6, *-5*, -4, 5, -3, 3, 0}", True,
       "zero-padded at both ends -- same signal"),
    _c("seq-07", "sequence", G_SEQ_NOORIG, "**Answer:** y[n] = {-1, 5, -4, 3}", True,
       "gold states no origin either; values compared"),
    _c("seq-08", "sequence", G_SEQ_NOORIG, "**Answer:** y[n] = {*0*, -1, 5, -4, 3}", True,
       "candidate states an origin gold does not; values still agree"),
    _c("seq-09", "sequence", G_SEQ,
       "## Final Answer\n**Answer:** z[n] = {-9, 6, *-5*, -4, 5, -3, 3} (dimensionless)", True,
       "trailing annotation"),
    _c("seq-10", "sequence", G_SEQ,
       "**Answer:** z[n] = { -9 , 6 , *-5* , -4 , 5 , -3 , 3 }", True,
       "loose whitespace"),

    # -- wrong -------------------------------------------------------------
    _c("seq-11", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, -5, *-4*, 5, -3, 3}", False,
       "**the origin moved one place; every value is identical.**  A values-only "
       "comparator accepts this, and gemma's archived trace 1 is exactly this error"),
    _c("seq-12", "sequence", G_SEQ, "**Answer:** z[n] = {3, -3, 5, -4, *-5*, 6, -9}", False,
       "reversed"),
    _c("seq-13", "sequence", G_SEQ, "**Answer:** z[n] = {6, *-5*, -4, 5, -3, 3, -9}", False,
       "rotated by one"),
    _c("seq-14", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, *-5*, -4, 5, -3}", False,
       "last element dropped"),
    _c("seq-15", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, *-5*, -4, 5, -3, 3, 7}", False,
       "spurious trailing non-zero element"),
    _c("seq-16", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, *-5*, 4, 5, -3, 3}", False,
       "one sign flipped"),
    _c("seq-17", "sequence", G_SEQ, "**Answer:** z[n] = {-9, 6, *-5*, -4, 5, -3, 30}", False,
       "one digit appended to the last value"),
    _c("seq-18", "sequence", G_SEQ_NOORIG, "**Answer:** y[n] = {-1, 5, 3, -4}", False,
       "two interior values transposed"),
    _c("seq-19", "sequence", G_SEQ, "**Answer:** the sequence is shifted three to the left", False,
       "prose, no sequence at all"),
    _c("seq-20", "sequence", G_SEQ,
       "**Answer:** z[n] = {-9, 6, -5, -4, 5, -3, 3}", False,
       "values right, origin never stated; gold pins it, so undetermined -- must not MATCH"),
    _c("seq-21", "sequence", G_SEQ,
       "**Answer:** z[n] = {-9, 6, -5, -4, 5, -3, 3} for n = 0, 1, 2, 3, 4, 5, 6", False,
       "origin stated and wrong"),
    _c("seq-22", "sequence", G_SEQ, "**Answer:** z[n] = {0, 0, 0, 0, 0, 0, 0}", False,
       "all zeros -- trims to empty, must not compare equal to gold"),
    _c("seq-23", "sequence", G_SEQ,
       "**Answer:** z[n] = {-9, 6, *-5*, 0, 5, -3, 3}", False,
       "an interior value replaced by zero; interior zeros are structural and must not be trimmed"),
    _c("seq-24", "sequence", G_SEQ_NOORIG, "**Answer:** y[n] = {0, 0, -1, 5, -4, 3, 0}", True,
       "leading and trailing padding around gold's values, no origin either side"),
]

# ==========================================================================
# symbolic -- 22 cases
# ==========================================================================

SYM = [
    _c("sym-01", "symbolic", G_SYM, "**Answer:** u = -2x^2 + 6xy", True, "identical"),
    _c("sym-02", "symbolic", G_SYM, "**Answer:** u(x, y) = 6xy - 2x^2", True, "terms reordered"),
    _c("sym-03", "symbolic", G_SYM, "**Answer:** u = −2x² + 6xy", True,
       "Unicode minus and superscript"),
    _c("sym-04", "symbolic", G_SYM, r"**Answer:** \(u = -2x^{2} + 6xy\)", True,
       "LaTeX braces on the exponent"),
    _c("sym-05", "symbolic", G_SYM, "**Answer:** u = -2x^2 + 6xy + f(y)", True,
       "arbitrary function retained -- the general solution"),
    _c("sym-06", "symbolic", G_SYM, "**Answer:** u = -2*x**2 + 6*x*y", True,
       "Python operator spelling"),
    _c("sym-07", "symbolic", G_SYM_FRAC, r"**Answer:** u = \frac{3x^2}{2} - 6xy", True,
       "3/2 written as a LaTeX fraction against a 1.5 gold"),
    _c("sym-08", "symbolic", G_SYM_FRAC, "**Answer:** u = 1.50x^2 - 6.0xy", True,
       "trailing zeros on the coefficients"),
    _c("sym-09", "symbolic", G_SYM, r"**Answer:** \boxed{u(x,y) = 6xy - 2x^2}", True,
       "boxed"),
    _c("sym-10", "symbolic", G_SYM, "**Answer:** u = 6·x·y - 2·x²", True,
       "middle-dot multiplication"),
    _c("sym-11", "symbolic", G_SYM_FRAC, "**Answer:** u(x, y) = 3x^2/2 - 6xy + C(y)", True,
       "trailing divisor plus arbitrary function"),

    # -- wrong -------------------------------------------------------------
    _c("sym-12", "symbolic", G_SYM, "**Answer:** u = 2x^2 + 6xy", False, "sign of the leading term"),
    _c("sym-13", "symbolic", G_SYM, "**Answer:** u = -2x^2 + 6x", False, "y dropped from the cross term"),
    _c("sym-14", "symbolic", G_SYM, "**Answer:** u = -2y^2 + 6xy", False,
       "x and y swapped in the squared term -- integrated in the wrong variable"),
    _c("sym-15", "symbolic", G_SYM, "**Answer:** u = -2x^2 + 6xy + 3", False,
       "spurious numeric constant, not an arbitrary function"),
    _c("sym-16", "symbolic", G_SYM, "**Answer:** u = -4x^2 + 6xy", False, "coefficient doubled"),
    _c("sym-17", "symbolic", G_SYM, "**Answer:** u(x, y) = C", False,
       "the arbitrary constant alone; the solution is dropped"),
    _c("sym-18", "symbolic", G_SYM, "**Answer:** u = -2x^3 + 6xy", False, "exponent wrong"),
    _c("sym-19", "symbolic", G_SYM_FRAC, "**Answer:** u = 3x^2 - 6xy", False,
       "forgot to halve on integration -- 3 against 1.5"),
    _c("sym-20", "symbolic", G_SYM, "**Answer:** du/dx = -4x + 6y", False,
       "answered with the derivative rather than the integral"),
    _c("sym-21", "symbolic", G_SYM, "**Answer:** the velocity is a quadratic in x and y", False,
       "prose, no expression"),
    _c("sym-22", "symbolic", G_SYM, "**Answer:** u = -2x^2 + 6xy + g(y) + 5", False,
       "arbitrary function AND a spurious constant"),
]

# ==========================================================================
# numeric -- 22 cases
# ==========================================================================

NUM = [
    _c("num-01", "numeric", G_NUM, "**Answer:** 7.65 L", True, "identical"),
    _c("num-02", "numeric", G_NUM, "**Answer:** V = 7.65 litres", True, "unit spelled out"),
    _c("num-03", "numeric", G_NUM, "**Answer:** 7.650 L", True, "extra trailing zero"),
    _c("num-04", "numeric", G_NUM, "**Answer:** 7.653 L", True,
       "within half a unit in gold's last displayed place"),
    _c("num-05", "numeric", G_NUM, "**Answer:** 7.655 L", True,
       "exactly on the inclusive boundary (S7.1)"),
    _c("num-06", "numeric", G_NUM, "**Answer:** 7.645 L", True, "boundary on the low side"),
    _c("num-07", "numeric", G_NUM, r"**Answer:** \(V = 7.65\ \text{L}\)", True, "LaTeX"),
    _c("num-08", "numeric", G_NUM, "**Answer:** The volume is 7.65 L (rounded).", True, "prose"),
    _c("num-09", "numeric", G_NUM_COMMA, "**Answer:** 4,921 units", True,
       "**the D-003 case**: thousands separator; the deployed parser reads this as 4.0"),
    _c("num-10", "numeric", G_NUM_COMMA, "**Answer:** 4921 units", True,
       "same value without the separator"),
    _c("num-11", "numeric", G_NUM, "**Answer:** 7.65", True, "no unit"),

    # -- wrong -------------------------------------------------------------
    _c("num-12", "numeric", G_NUM, "**Answer:** 7.66 L", False, "one display unit out"),
    _c("num-13", "numeric", G_NUM, "**Answer:** 7.656 L", False, "just outside the boundary"),
    _c("num-14", "numeric", G_NUM, "**Answer:** 76.5 L", False, "factor of ten"),
    _c("num-15", "numeric", G_NUM, "**Answer:** -7.65 L", False, "sign"),
    _c("num-16", "numeric", G_NUM, "**Answer:** 7.65 mL", False,
       "right number, wrong unit; caught only because gold declares its unit",
       options={"unit": "L"}),
    _c("num-16b", "numeric", G_NUM, "**Answer:** 7.65 mL", False,
       "the same case with NO declared unit -- **not caught**, and this is the "
       "residual D3.4 S10 leaves open (D-052)"),
    _c("num-16c", "numeric", G_NUM, "**Answer:** 7.65 litres", True,
       "unit spelled out, declared unit still matches -- the case that makes "
       "inferring the unit from gold's string the wrong design",
       options={"unit": "L"}),
    _c("num-17", "numeric", G_NUM_COMMA, "**Answer:** 4.921 units", False,
       "the separator read as a decimal point -- the D-003 defect as a candidate answer"),
    _c("num-18", "numeric", G_NUM_COMMA, "**Answer:** 4 units", False,
       "exactly what the deployed parser extracts from gold; must not match gold"),
    _c("num-19", "numeric", G_NUM, "**Answer:** approximately 7.7 L", False,
       "rounded to fewer places than gold displays"),
    _c("num-20", "numeric", G_NUM, "**Answer:** somewhere between 7 and 8 L", False,
       "an interval, not a value"),
    _c("num-21", "numeric", G_NUM, "**Answer:** see step 3", False, "no number"),
    _c("num-22", "numeric", G_NUM, "**Answer:** 7.65e0 L", True, "exponent form of the same value"),
]

# ==========================================================================
# check -- 21 cases
# ==========================================================================

CHECK = [
    _c("chk-01", "check", G_CHECK, "**Answer:** 18.4 mm < 25 mm, so the design is acceptable.", True,
       "identical shape"),
    _c("chk-02", "check", G_CHECK, "**Answer:** Yes -- 18.4 mm is within the limit.", True,
       "'Yes' as the verdict surface"),
    _c("chk-03", "check", G_CHECK, "**Answer:** The deflection of 18.4 mm is safe.", True, "'safe'"),
    _c("chk-04", "check", G_CHECK, "**Answer:** 18.42 mm; the criterion holds.", True,
       "quantity within tolerance, 'holds'"),
    _c("chk-05", "check", G_CHECK, "**Answer:** The design passes (18.4 mm).", True, "'passes'"),
    _c("chk-06", "check", G_CHECK, "**Answer:** 18.4 mm -- valid.", True, "'valid'"),
    _c("chk-07", "check", G_CHECK, "**Answer:** The deflection is 18.4 mm, which is acceptable.", True,
       "verdict after the quantity"),
    _c("chk-08", "check", G_CHECK, "**Answer:** Adequate; deflection 18.4 mm.", True,
       "verdict before the quantity"),

    # -- wrong -------------------------------------------------------------
    _c("chk-09", "check", G_CHECK, "**Answer:** 18.4 mm exceeds the limit, so the design is not acceptable.", False,
       "verdict inverted, quantity right"),
    _c("chk-10", "check", G_CHECK, "**Answer:** The design is unsafe.", False, "verdict inverted, no quantity"),
    _c("chk-11", "check", G_CHECK, "**Answer:** 40.2 mm, so the design is acceptable.", False,
       "**the coin-flip case**: right verdict, wrong number.  Scoring the boolean "
       "alone credits reasoning that did not happen on a two-valued answer"),
    _c("chk-12", "check", G_CHECK, "**Answer:** 31.0 mm exceeds 25 mm; not acceptable.", False,
       "wrong verdict AND wrong number, consistent with each other"),
    _c("chk-13", "check", G_CHECK, "**Answer:** The design is acceptable.", False,
       "verdict right, gold's supporting quantity not stated"),
    _c("chk-14", "check", G_CHECK, "**Answer:** It seems acceptable.", False, "hedged"),
    _c("chk-15", "check", G_CHECK, "**Answer:** It may or may not be acceptable.", False,
       "hedged both ways"),
    _c("chk-16", "check", G_CHECK, "**Answer:** 18.4 mm.", False, "quantity only, no verdict"),
    _c("chk-17", "check", G_CHECK, "**Answer:** The beam is 18.4 m long.", False,
       "a number and no verdict"),
    _c("chk-18", "check", G_CHECK, "**Answer:** 18.4 mm is less than 25 mm.", False,
       "states the comparison but never gives the verdict word"),
    _c("chk-19", "check", G_CHECK, "**Answer:** 18.35 mm, acceptable.", True,
       "quantity on the inclusive boundary"),
    _c("chk-20", "check", G_CHECK, "**Answer:** 18.46 mm, acceptable.", False,
       "quantity just outside tolerance with the right verdict"),
    _c("chk-21", "check", G_CHECK, "**Answer:** Not unacceptable.", False,
       "double negative that a human grader would not accept as a committed verdict"),
]

# ==========================================================================
# narrative -- 20 cases
#
# Every one of these must return UNRESOLVED.  That is the kind's entire
# specification (D-051), and the set exists to prove the comparator never
# quietly accepts or rejects prose -- including prose that is obviously right
# and prose that is obviously wrong.
# ==========================================================================

NARR = [
    _c(f"nar-{i:02d}", "narrative", G_NARR, cand, False, note)
    for i, (cand, note) in enumerate([
        ("**Answer:** Turbulent transition raises wall shear, so pressure drop rises.", "a paraphrase of gold"),
        ("**Answer:** " + G_NARR.split("\n")[-1], "gold restated verbatim"),
        ("**Answer:** The pressure drop falls because the flow becomes laminar.", "the opposite of gold"),
        ("**Answer:** Because Re crosses 2300.", "correct but far shorter"),
        ("**Answer:** The viscosity increases with temperature.", "unrelated but plausible"),
        ("**Answer:** See the discussion in Step 4.", "a reference, not an answer"),
        ("**Answer:** ", "empty"),
        ("**Answer:** 2300", "a number where prose is expected"),
        ("**Answer:** Yes.", "a boolean where prose is expected"),
        ("**Answer:** The wall shear increases.", "half of gold's explanation"),
        ("**Answer:** Turbulence.", "one word"),
        ("**Answer:** It rises.", "the direction, without the mechanism"),
        ("**Answer:** The transition to turbulent flow increases wall shear, which raises pressure drop.", "gold's clauses reordered"),
        ("**Answer:** Pressure drop rises due to increased wall shear from turbulent transition.", "gold, nominalised"),
        ("**Answer:** I am not sure.", "an explicit refusal"),
        ("**Answer:** The pressure drop rises.", "the conclusion without the reason"),
        ("**Answer:** Because the Darcy friction factor increases.", "a different but valid mechanism"),
        ("**Answer:** Turbulent flow has a lower friction factor, so drop falls.", "wrong physics, confident"),
        ("**Answer:** The answer depends on the pipe roughness.", "deflection to a parameter"),
        ("**Answer:** " + "x" * 400, "very long noise"),
    ], start=1)
]

ALL = CATEGORICAL + TUPLE + SEQ + SYM + NUM + CHECK + NARR

#: Which template binding (or explicit options) each kind is exercised through.
KIND_TEMPLATE = {
    "categorical": "template_system_property_linearity",
    "categorical[tuple]": "template_system_properties_memory_causality",
    "sequence": "template_signal_operations",
    "symbolic": "template_incompressible_continuity",
}


def build() -> list[dict]:
    out = []
    for c in ALL:
        c = dict(c)
        c["template"] = c["template"] or KIND_TEMPLATE.get(c["kind"])
        out.append(c)
    return out


def main() -> None:
    cases = build()
    path = os.path.join(os.path.dirname(__file__), "adversarial.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cases, fh, indent=1, ensure_ascii=False)
    from collections import Counter
    per_kind = Counter(c["kind"] for c in cases)
    pos = Counter(c["kind"] for c in cases if c["correct"])
    print(f"wrote {len(cases)} cases to {path}")
    print(f"{'kind':22s} {'total':>6s} {'correct':>8s} {'wrong':>6s}")
    for k, n in sorted(per_kind.items()):
        print(f"{k:22s} {n:6d} {pos[k]:8d} {n - pos[k]:6d}")
    short = [k for k, n in per_kind.items() if n < 20]
    print("\nD4.4 requires >=20 per comparator: "
          + ("all met" if not short else f"SHORT: {short}"))


if __name__ == "__main__":
    main()
