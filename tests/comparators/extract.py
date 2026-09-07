"""D5.6 -- which number in an answer span is the answer.

`parse_number` read the **first** number, and the first number in an answer
sentence is very often not the answer: it is a fluid grade (`Engine Oil (SAE
50)`), a temperature (`at 541 K`), or a subscript from a chemical formula (the
4 of C4H10).  That is N1, and it is D-003's defect class reappearing inside the
comparator built to replace it.

**The rule below was chosen by measurement, not by example.**  Eleven candidate
rules were scored on two axes over four corpora with a held-out slice frozen
before any of them was written; the losing candidates' numbers are in
`n1_candidates.py`, which regenerates the whole table.  Three findings shaped
what is here, and none of them was visible from a single example:

1. **"Last number" is much worse than "first", not better.**  424 gold-by-gold
   false accepts on the held-out slice against the incumbent's 16.  It reads
   the limit out of `a rate of 3.1 mm/s at 300 K` and the `2` out of `m^2/s`.

2. **The tokeniser matters more than the choice of number.**  A digit inside a
   unit exponent (`m^2`) or a formula subscript (`C4H10`) is not a quantity at
   all, and every rule that picks *any* number is exposed to it.  Removing
   those two classes is orthogonal to which of the remaining numbers a rule
   picks, so it is applied under every rule equally -- and it is what turns
   "last number" from 424 false accepts into 0.

3. **The rule must be per-KIND, and this is the whole reframing.**  `check`
   answers state the quantity *first* and the threshold second -- *"the
   deflection is 18.4 mm, less than the 25 mm limit"* -- so every last-ward rule
   picks the limit.  Measured on D4.4's 21 `check` cases: first-number scores
   13 decided / 13 correct / 0 wrong; every last-ward rule scores 13 decided /
   **5** correct / 8 false rejects.  The comparator did not need another rule;
   it needed the binding to say which rule applies.

Adopted:

    numeric   ANSWER_RULE  -- declared unit, else assignment anchor, else the
                             last number over the filtered tokeniser
    check     first        -- the quantity precedes the threshold it is
                             compared against

Held-out results for the adopted rule, against the incumbent:

    ==========================  ======  =======  ======  =======
    rule                          ggFA   ggDEC%    agXI   agDEC%
    ==========================  ======  =======  ======  =======
    first (incumbent)               16    100.0     144     99.1
    ADOPTED                          0    100.0       6     99.1
    ==========================  ======  =======  ======  =======

`ggFA` is gold-by-gold false accepts; `agXI` is cross-instance accepts on real
archived model text.  **The decided rate did not fall**, which is the gate item
an accept-only comparison cannot see -- a rule that decides nothing has zero
false accepts.  All six remaining `agXI` were audited individually and none is
an extraction error; see `n1_candidates.py`.
"""

from __future__ import annotations

import re
from decimal import Decimal, InvalidOperation

#: The number grammar.  Thousands separators are honoured: `4,921` is four
#: thousand nine hundred and twenty-one, and the deployed parser reads it as
#: `4.0` (D-003).
NUM_RE = re.compile(r"[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][-+]?\d+)?")

#: A token that may follow a number and name its unit.
_UNIT_AFTER = re.compile(r"\s*[>\)\]]?\s*([A-Za-zµΩ°][A-Za-z0-9µΩ°/^*_.-]{0,13})")

#: Where an answer is asserted.  The number after the LAST of these is the one
#: the sentence commits to.
_ANCHOR = re.compile(r"(?:=|\bis\b|\bare\b|\bequals\b)")


def numbers(text: str) -> list[re.Match]:
    """``NUM_RE`` matches that could be a quantity, with two classes removed.

    * a digit after ``^`` or ``_`` -- ``m^2``, ``x_1``: an exponent or an index,
      never a quantity;
    * a digit *between* two letters -- the 4 and 10 of ``C4H10``: a chemical
      formula subscript.

    Both are recorded in D-058 as numbers the comparator read as answers.  No
    downstream rule can recover from having chosen one, so they are removed
    here rather than worked around later.
    """
    out = []
    for m in NUM_RE.finditer(text):
        i = m.start()
        if i and text[i - 1] in "^_":
            continue
        if i and text[i - 1].isalpha() and m.end() < len(text) and text[m.end()].isalpha():
            continue
        out.append(m)
    return out


def _unit_adjacent(text: str, unit: str | None) -> re.Match | None:
    """The last number the *declared* unit follows.

    Declared, never inferred (D-052): the caller passes the unit gold declares.
    Inferring it from gold's own string would make the comparison agree with
    whatever gold happens to say, which is the answer-key failure D-034 names.
    """
    if not unit:
        return None
    hit = None
    for m in NUM_RE.finditer(text):
        u = _UNIT_AFTER.match(text, m.end())
        if u and u.group(1).rstrip(".").lower() == unit.rstrip(".").lower():
            hit = m
    return hit


def _after_anchor(text: str) -> re.Match | None:
    anchors = list(_ANCHOR.finditer(text))
    if not anchors:
        return None
    after = [m for m in numbers(text) if m.start() >= anchors[-1].end()]
    return after[0] if after else None


def _last_filtered(text: str) -> re.Match | None:
    ms = numbers(text)
    return ms[-1] if ms else None


def answer_match(text: str, kind: str = "numeric", unit: str | None = None) -> re.Match | None:
    """The match for the number this span asserts as its answer, or None.

    None means *this rule refuses to decide*, and the caller must render that
    as UNRESOLVED -- never as a guess and never as a match.
    """
    if kind == "check":
        # The quantity precedes the threshold; see the module docstring.
        ms = numbers(text)
        return ms[0] if ms else None
    return _unit_adjacent(text, unit) or _after_anchor(text) or _last_filtered(text)


def answer_number(text: str, kind: str = "numeric", unit: str | None = None) -> Decimal | None:
    m = answer_match(text, kind, unit)
    if m is None:
        return None
    try:
        return Decimal(m.group(0).replace(",", ""))
    except InvalidOperation:
        return None


def displayed_decimals(tok: str) -> int | None:
    """Decimal place of ``tok``'s last displayed digit, exponent form included.

    The incumbent returned None for anything in exponent form, so a gold printed
    `3.1259e-07` had no derivable tolerance and **every** comparison against it
    was UNRESOLVED.  That is a gap, not a policy: D4.1 |S|7.1 defines the display
    tolerance by where gold's last displayed digit falls, and scientific
    notation states that as precisely as fixed-point does.

    `3.1259e-07` shows four mantissa decimals at 10^-7, so its last digit is the
    10^-11 place: 4 - (-7) = 11.  `1.5e3` shows to the hundreds: 1 - 3 = -2.

    Measured: closing this recovered the archive decided rate from 91.6% to
    97.7% under *every* rule including the incumbent, so it is an independent
    improvement and not a way of paying for the N1 fix.
    """
    tok = tok.replace(",", "")
    if "e" in tok.lower():
        mant, _, exp = tok.lower().partition("e")
        dec = len(mant.split(".")[1]) if "." in mant else 0
        try:
            return dec - int(exp)
        except ValueError:
            return None
    return len(tok.split(".")[1]) if "." in tok else 0


def answer_decimals(text: str, kind: str = "numeric", unit: str | None = None) -> int | None:
    """Displayed precision of the *same* number ``answer_number`` chose.

    These two must move together.  Deriving the value from one number and the
    tolerance from another would judge an answer against a precision it never
    had; the incumbent could not exhibit that only because both read the first
    number.
    """
    m = answer_match(text, kind, unit)
    return None if m is None else displayed_decimals(m.group(0))

# --------------------------------------------------------------------------
# The cross-pairing truth predicate (Reviewer E, R2-F2).
# --------------------------------------------------------------------------

def same_answer(gold_span: str, other_span: str) -> bool:
    """Do two gold answer spans state the **same answer**?

    Cross-pairing needs a truth it can compute without a label, and the first
    version used **textual identity of the span**.  That is wrong in one
    direction and it was costing real bindings: gold ``5.0 seconds`` against
    gold ``4.99 seconds`` differs textually, and under |S|7.1's
    boundary-inclusive display tolerance ``4.99`` **is** a correct answer to a
    question whose gold shows ``5.0`` -- so a MATCH there was being counted as
    a false accept, and **10 of 32 templates were unbound on that predicate,
    9 of them on it alone**.

    So the predicate is: the two spans state the same answer iff they have the
    same non-numeric text, the same count of numbers, and every number agrees
    with its counterpart **within the tolerance implied by the GOLD side's
    displayed precision**.

    **It does not call the answer rule**, deliberately.  Using
    ``answer_match`` to decide the truth that ``answer_match`` is being scored
    against is a test carrying its own answer key (D-034); this walks every
    number in the span instead, so it is independent of which one the
    comparator picks.
    """
    ga = [m.group(0) for m in NUM_RE.finditer(gold_span)]
    ca = [m.group(0) for m in NUM_RE.finditer(other_span)]
    if len(ga) != len(ca):
        return False
    if NUM_RE.sub("#", gold_span) != NUM_RE.sub("#", other_span):
        return False
    for g_tok, c_tok in zip(ga, ca):
        if g_tok == c_tok:
            continue
        p = displayed_decimals(g_tok)
        if p is None:
            return False
        try:
            gv = Decimal(g_tok.replace(",", ""))
            cv = Decimal(c_tok.replace(",", ""))
        except InvalidOperation:
            return False
        if abs(cv - gv) > Decimal("0.5") * (Decimal(10) ** -p):
            return False
    return True
