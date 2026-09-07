"""The six ``kind`` comparators (D4.1 / D4.3).

    numeric  categorical  sequence  symbolic  narrative  check

``multipart`` is **not** a seventh kind -- see ``answer.py`` and D-047.

Two rules run through all six, both inherited from Phase 3:

* **The gold node owns every parameter of the comparison.**  Precision, the
  label set, the origin convention, the symbol alphabet: all read from gold,
  never from the candidate.  This is D3.4 |S|7.1's rule ("always the GOLD node's
  ``symbol_precision``, never the candidate's") generalised to all six kinds.  A
  candidate that could widen its own tolerance is a candidate that grades
  itself.
* **A value the comparator consumes is declared once, or recomputed** (|S|3.8).
  The comparators hold no reference answers: gold canonical forms are derived
  from the gold string by the same code path that derives the candidate's
  (D-034).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from typing import Any, Sequence

from . import commitment
from .commitment import (
    clauses as commitment_clauses,
    find_commitment,
    hedge_markers,
    is_assertion,
    segment,
    suspends_what_follows,
)
from .normalize import (
    NEGATORS,
    NEITHER_NOR_RE,
    answer_span,
    prepare,
)
from .extract import (  # D5.6 -- the measured extraction rule
    answer_decimals,
    answer_number,
    displayed_decimals,
)
from .verdict import MISMATCH, UNRESOLVED, Verdict, match, mismatch, unresolved

# ==========================================================================
# 1.  numeric
# ==========================================================================

_NUM_RE = re.compile(r"[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][-+]?\d+)?")


def _decimals(text: str, kind: str = "numeric", unit: str | None = None) -> int | None:
    """Displayed decimal places of the number the ANSWER RULE chose (D5.6)."""
    return answer_decimals(text, kind, unit)


def parse_number(text: str, kind: str = "numeric", unit: str | None = None) -> Decimal | None:
    """The number this span asserts as its answer (D5.6).

    Was: the FIRST number, which is often a fluid grade, a temperature or a
    chemical-formula subscript -- defect N1 in D-058.  The rule now applied was
    chosen by measurement over four corpora with a frozen held-out slice; see
    ``extract.py`` and ``n1_candidates.py``.
    """
    return answer_number(text, kind, unit)


#: Unit surfaces, keyed by the canonical unit.  Only the ones a caller
#: declares are ever consulted; there is no inferred unit anywhere.
UNIT_SURFACES: dict[str, tuple[str, ...]] = {
    "L": ("L", "l", "litre", "litres", "liter", "liters", "dm^3"),
    "mL": ("mL", "ml", "millilitre", "millilitres", "milliliter", "milliliters", "cm^3"),
    "m": ("m", "metre", "metres", "meter", "meters"),
    "mm": ("mm", "millimetre", "millimetres", "millimeter", "millimeters"),
    "Pa": ("Pa", "pascal", "pascals"),
    "kPa": ("kPa", "kilopascal", "kilopascals"),
    "N": ("N", "newton", "newtons"),
    "kN": ("kN", "kilonewton", "kilonewtons"),
    "rad/s": ("rad/s", "rad s^-1", "radians per second"),
    "Hz": ("Hz", "hertz", "cycles per second"),
    "s": ("s", "sec", "secs", "second", "seconds"),
}


def compare_numeric(
    gold: str,
    candidate: str,
    *,
    precision: int | None = None,
    method_tolerance: float | None = None,
    unit: str | None = None,
    extract: bool = True,
) -> Verdict:
    """|S|7.1 display tolerance, boundary inclusive, widened by |S|7.4 if given.

    ``precision`` is the gold display precision.  When omitted it is *derived
    from gold's own rendering*, never from the candidate's -- a candidate
    quoting more digits must not thereby tighten the test, and one quoting
    fewer must not loosen it.

    **Units are checked only when declared.**  D3.4 |S|10 records "no dimensional
    checking" as a non-goal, because a dimensional comparator needs units on
    every symbol and that is a milestone-model decision.  The residue is a real
    false accept: ``7.65 mL`` against a gold of ``7.65 L`` matches on the
    number (D4.4 case num-16).  Passing ``unit`` closes it for that item by
    *declaration* -- gold names its unit and the surfaces that spell it -- and
    leaves it open, visibly, for every item that does not.  Inferring the unit
    from gold's string instead would reject ``7.65 litres``, which is correct;
    that is the trade, and it is why this is opt-in rather than default.
    """
    k = "numeric"
    g_txt = prepare(answer_span(gold)[0] if extract else gold)
    c_txt = prepare(answer_span(candidate)[0] if extract else candidate)

    g = parse_number(g_txt, "numeric", unit)
    if g is None:
        return unresolved(k, "no number in the gold answer span", gold_canonical=g_txt)
    c = parse_number(c_txt, "numeric", unit)
    if c is None:
        return unresolved(
            k, "no number in the candidate answer span",
            gold_canonical=str(g), cand_canonical=c_txt,
        )

    p = precision if precision is not None else _decimals(g_txt, "numeric", unit)
    if p is None:
        return unresolved(
            k, "gold is in exponent form and no precision was declared",
            gold_canonical=str(g), cand_canonical=str(c),
        )

    tol = Decimal("0.5") * (Decimal(10) ** -p)
    obs: list[str] = []
    if method_tolerance is not None:
        mt = Decimal(str(method_tolerance))
        if mt > tol:
            obs.append(f"tolerance widened {tol} -> {mt} by the prescribed method (S7.4)")
            tol = mt

    delta = abs(c - g)
    if delta > tol:  # boundary inclusive, per |S|7.1
        return mismatch(
            k, f"|{c} - {g}| = {delta} > {tol}",
            gold_canonical=str(g), cand_canonical=str(c), observations=obs,
        )

    if unit is not None:
        cu = _resolve_unit(c_txt)
        if cu is None:
            return unresolved(
                k, f"gold declares the unit {unit!r} and the candidate states none",
                gold_canonical=f"{g} {unit}", cand_canonical=str(c), observations=obs,
            )
        if cu != unit:
            return mismatch(
                k, f"unit {cu!r} != declared gold unit {unit!r}",
                gold_canonical=f"{g} {unit}", cand_canonical=f"{c} {cu}", observations=obs,
            )
        obs.append(f"unit {cu} agrees with the declared gold unit")

    return match(k, gold_canonical=str(g), cand_canonical=str(c), observations=obs)


def _resolve_unit(text: str) -> str | None:
    """The canonical unit named by ``text``, or None.

    Longest surface wins, so ``mL`` is preferred over the ``L`` inside it --
    the same rule, and the same reason, as ``nonlinear`` over ``linear``.
    """
    best: tuple[int, str] | None = None
    for canon, surfaces in UNIT_SURFACES.items():
        for s in surfaces:
            for m in re.finditer(r"(?<![\w/])" + re.escape(s) + r"(?![\w])", text):
                if best is None or len(s) > len(best[1]):
                    best = (m.start(), s)
                    best_canon = canon
    return None if best is None else best_canon


# ==========================================================================
# 2.  categorical
# ==========================================================================


def _label_hit(text: str, surfaces: Sequence[str]) -> tuple[int, str] | None:
    """Position and matched surface of the last occurrence of any surface.

    **Position first; length only to break a tie at the same offset.**

    The first version sorted surfaces longest-first and stopped at the first
    that matched *anywhere*, so a single ``nonlinear`` outranked a committed
    ``linear`` no matter where each sat. Reviewer E (F2) and Reviewer B (F4)
    found this independently, from opposite directions, and both are right:
    D4.1 §4.2's longest-first rule is about **containment at one offset**
    (``linear`` inside ``nonlinear``) and was being applied **across** offsets,
    which is a different and false claim.

    Note that containment between these two particular surfaces is already
    impossible -- ``(?<![\\w-])`` stops ``linear`` matching inside
    ``nonlinear`` -- so the length rule only ever fires on genuine ties, which
    is what it was for.

    This alone does not settle which clause the model committed to; see
    ``commitment.find_commitment``, which does.
    """
    best: tuple[int, str] | None = None
    for s in surfaces:
        for m in re.finditer(r"(?<![\w-])" + re.escape(s) + r"(?![\w-])", text, re.I):
            if best is None or (m.start(), len(s)) > (best[0], len(best[1])):
                best = (m.start(), s)
    return best


_NEG_WINDOW = 40


def _negated(text: str, at: int, surface: str) -> bool:
    """Is the label at ``at`` inside the scope of a negator?

    Three carriers, all observed:
      * a negator in the ~40 characters before it (``is not causal``);
      * a detached fused prefix (``non- linear`` after markdown stripping);
      * ``neither A nor B``, which negates both labels and is handled by the
        caller because its scope spans two of them.

    A **fused** negative (``nonlinear``, ``non-causal``) is deliberately *not*
    treated as a negation here.  Those strings are declared as surfaces of the
    negative label in their own right, so inferring a negation on top of the
    declaration double-negates and returns the positive label -- which is how
    the first draft of this file scored ``The system is nonlinear`` as
    ``linear`` on three of the fifteen archived linearity traces.  Declaration
    beats inference: if a label set wants a fused form, it declares it.
    """
    window = text[max(0, at - _NEG_WINDOW):at].lower()
    # A clause boundary ends a negator's scope: "is not memoryless. It is
    # causal" must not negate "causal".
    # The colon is a clause boundary too, and so are the concessives.  Without
    # them, "it does not break superposition: the system is linear" let the
    # negation reach across and scored a CORRECT answer as WRONG -- a MISMATCH,
    # not merely an UNRESOLVED, which D4.1 §1's stated bias rules out.
    # Reviewer B filed it in round 2 and again as R3-F3 when it had spread to
    # `Although`, the opener the segmentation rewrite existed to rehabilitate.
    window = re.split(r"[.;:]|\bbut\b|\bhowever\b|\bwhereas\b|\bwhile\b|"
                      r"\balthough\b|\bthough\b|\byet\b", window)[-1]
    if re.search(r"\bnon-?\s*$", window):
        return True
    if any(re.search(r"\b" + re.escape(n) + r"\b", window) for n in NEGATORS):
        return True
    # **Negation can follow the label as well as precede it.**  "linear is the
    # wrong description" asserts the opposite of the label it names.  Ranking
    # hits by position rather than by surface length (Reviewer E F2 /
    # Reviewer B F4) made this reachable: the later, post-modified mention now
    # wins, where the longest-surface rule used to hide it behind an earlier
    # `nonlinearity`.  Narrow by design -- a predicate that says the label is
    # wrong, nothing more.
    after = text[at + len(surface):at + len(surface) + _NEG_WINDOW].lower()
    after = re.split(r"[.;]", after)[0]
    return bool(_POST_NEG_RE.search(after))


#: Predicates that negate the label immediately preceding them.
_POST_NEG_RE = re.compile(
    r"^\W*(?:is|was|would\s+be|seems)?\s*(?:the\s+)?"
    r"(?:wrong|incorrect|mistaken|not\s+correct|not\s+right|not\s+the\s+case|"
    r"inapplicable|does\s+not\s+apply|fails|is\s+ruled\s+out)\b",
    re.IGNORECASE,
)


def compare_categorical(
    gold: str,
    candidate: str,
    *,
    labels: dict[str, Sequence[str]],
    gold_label: str | None = None,
    extract: bool = True,
) -> Verdict:
    """Compare a single categorical label drawn from a **closed, declared** set.

    ``labels`` maps each canonical label to its accepted surface forms; it is
    problem data on the gold side and is the only vocabulary in play.  A
    candidate whose answer span names no member of the set is ``UNRESOLVED``,
    not a mismatch -- "it said something we do not recognise" and "it said the
    wrong thing" are different observations and collapsing them hides
    comparator defects inside the model's error rate.
    """
    k = "categorical"
    g_txt = prepare(answer_span(gold)[0] if extract else gold)
    c_txt = prepare(answer_span(candidate)[0] if extract else candidate)

    g = gold_label if gold_label is not None else _resolve_label(g_txt, labels)
    if g is None:
        return unresolved(k, "gold answer names no declared label", gold_canonical=g_txt)

    # Which clause did the candidate COMMIT to?  Not "does a label appear
    # somewhere unhedged" -- that question let a conditional, an assumption, a
    # refusal, the item's own prompt wording, a trailing hedge, an unlisted
    # hedge and a contrastive mention all through (Reviewer B F1-F4, Reviewer E
    # F2).  See commitment.py for the mechanism and why it replaced a blocklist.
    surfaces = [x for ss in labels.values() for x in ss]
    commit = find_commitment(c_txt, surfaces)
    if not commit:
        return unresolved(
            k, f"candidate does not commit to a label: {commit.reason}",
            gold_canonical=g, cand_canonical=c_txt,
        )
    c = _resolve_label(commit.clause, labels)
    if c is None:
        return unresolved(
            k, "candidate's committed clause names no declared label",
            gold_canonical=g, cand_canonical=commit.clause,
        )
    if c == g:
        return match(k, gold_canonical=g, cand_canonical=c)
    return mismatch(k, f"label {c!r} != gold {g!r}", gold_canonical=g, cand_canonical=c)


def _resolve_label(text: str, labels: dict[str, Sequence[str]]) -> str | None:
    """Canonical label named by ``text``, or None.

    A label set is declared as ``{canonical: [surfaces...]}``.  Surfaces may be
    positive (``linear``) or explicitly negative (``not linear``); the negation
    scoping applies to positive surfaces only, and a set that declares both
    polarities explicitly is resolved by position, latest wins.
    """
    hit = _label_hit(text, [s for ss in labels.values() for s in ss])
    if hit is None:
        return None
    at, surface = hit
    owner = next(c for c, ss in labels.items() if surface in ss)

    negated = _negated(text, at, surface)
    if not negated and not surface.lower().startswith("non"):
        m = NEITHER_NOR_RE.search(text)
        if m and m.start() <= at <= m.end():
            negated = True

    if not negated:
        return owner
    # The negation of a two-valued label set is its other member.  For a set
    # with more than two members, negation does not identify one -- "not
    # laminar" leaves transitional and turbulent open -- so it is unresolvable
    # rather than silently the complement.
    others = [c for c in labels if c != owner]
    return others[0] if len(others) == 1 else None


@dataclass(frozen=True)
class PropertySlot:
    """One yes/no property in a label-tuple answer.

    ``positive`` and ``negative`` are the surfaces that assert and deny the
    property.  Both are declared because the archive uses both: models write
    ``not memoryless`` (negated positive) and ``it has memory`` (declared
    negative) for the same answer, and ``non-causal`` for the other.
    """

    name: str
    positive: tuple[str, ...]
    negative: tuple[str, ...] = ()


_YES_RE = re.compile(r"(?<![\w-])yes(?![\w-])", re.I)
_NO_RE = re.compile(r"(?<![\w-])no(?![\w-])", re.I)


def compare_label_tuple(
    gold: str,
    candidate: str,
    *,
    slots: Sequence[PropertySlot],
    extract: bool = True,
) -> Verdict:
    """An ordered tuple of yes/no property labels.

    This is ``system_properties_memory_causality``: the answer is
    ``(Memoryless: No, Causal: Yes)``.  It is a *composite of two categorical
    parts*, not a kind of its own -- see D-047, for which this template is the
    worked proof.

    **Two encodings are in play and gold uses the rarer one.**  Gold writes
    ``a) Memoryless: **No**``: the slot is *named* and answered yes/no.  Models
    overwhelmingly write ``a) Not memoryless``: the property word itself is
    asserted or denied.  A comparator that reads only gold's encoding scores
    the whole archive wrong -- 21 of 24 traces, in the first draft of this file.
    Both are resolved, and a candidate is free to use either.
    """
    k = "categorical[tuple]"
    g_txt = prepare(answer_span(gold)[0] if extract else gold)
    c_txt = prepare(answer_span(candidate)[0] if extract else candidate)

    # **A frame around an enumerated answer still frames it.**
    # Slots are located positionally among the clauses that FOLLOW an
    # enumerator, so any prose before the first `a)` is dropped -- and a
    # hypothetical or task-restating preface was dropped with it, making
    # "If the additivity test holds, a) Yes, b) Yes" a MATCH.  That is my
    # round-2 fix for the recall corpus relocating a defect into the tuple
    # path, and Reviewer B found it with the negative control the corpus
    # itself lacked (R3-F1).  The preface is judged before the slots are.
    preface = _clauses_marked(c_txt)
    lead = " ".join(t for t, was_enum in preface if not was_enum).strip()
    if lead:
        why = suspends_what_follows(lead)
        if why:
            return unresolved(k, why)

    # A label-free question after the slots withdraws them, the same way it
    # withdraws a single label (see `_hedges_governing`).  The slot path
    # resolves each slot in its own clause and so never sees a trailing
    # sentence that belongs to neither.
    for c in commitment_clauses(c_txt):
        if c.rstrip().endswith("?") and not any(
            _label_hit(c, list(s.positive) + list(s.negative)) for s in slots
        ):
            return unresolved(k, "the answer is withdrawn by a following question")

    g_clauses, c_clauses = _clauses(g_txt), _clauses(c_txt)
    parts = [
        _slot_verdict(s, i, g_txt, c_txt, g_clauses, c_clauses)
        for i, s in enumerate(slots)
    ]
    if any(p.outcome == UNRESOLVED for p in parts):
        bad = [p for p in parts if p.outcome == UNRESOLVED]
        return Verdict(UNRESOLVED, k,
                       reason="; ".join(f"{b.kind}: {b.reason}" for b in bad),
                       parts=parts,
                       gold_canonical=[p.gold_canonical for p in parts],
                       cand_canonical=[p.cand_canonical for p in parts])
    if all(p.is_match for p in parts):
        return match(k, parts=parts,
                     gold_canonical=[p.gold_canonical for p in parts],
                     cand_canonical=[p.cand_canonical for p in parts])
    wrong = [p.kind for p in parts if not p.is_match]
    return mismatch(k, f"slot(s) differ: {', '.join(wrong)}", parts=parts,
                    gold_canonical=[p.gold_canonical for p in parts],
                    cand_canonical=[p.cand_canonical for p in parts])


def _slot_verdict(
    slot: PropertySlot, pos: int, g_txt: str, c_txt: str,
    g_clauses: list[str], c_clauses: list[str],
) -> Verdict:
    g = _resolve_property(_slot_text(g_txt, slot, pos, g_clauses), slot)
    if g is None:
        return unresolved(slot.name, "gold slot states no value")
    ctext = _slot_text(c_txt, slot, pos, c_clauses)
    # Scope to the sentence that actually carries this slot.  `_slot_text`
    # splits on `;`, newlines and enumerators but not on sentence stops, so
    # "Causal: Yes. Both look immediate." arrived as one slot text and the
    # trailing remark's `look` hedged a committed answer (Reviewer B, R2-F1.3).
    surfaces = list(slot.positive) + list(slot.negative) + [slot.name]
    owned = [x for x in commitment_clauses(ctext) if _label_hit(x, surfaces)]
    ctext = owned[-1] if owned else ctext
    matrix, why = segment(ctext)
    if matrix is None:
        return unresolved(slot.name, why)
    # Policy-aware, like every other hedge site (D-056): under the shipped
    # advisory default a hedge is annotated, not scored.  Leaving this one
    # enforcing made the tuple path disagree with the categorical path about
    # the same sentence -- the shape of Reviewer E's R3-F16.
    marks = hedge_markers(matrix)
    if marks and commitment.HEDGE_POLICY == "enforce":
        return unresolved(slot.name, f"hedged ({', '.join(marks)})")
    c = _resolve_property(matrix, slot)
    if c is None:
        return unresolved(slot.name, "candidate slot states no value")
    if c == g:
        return match(slot.name, gold_canonical=g, cand_canonical=c)
    return mismatch(slot.name, f"{c!r} != {g!r}", gold_canonical=g, cand_canonical=c)


_ENUM_RE = re.compile(r"(?:^|[\s,;])\(?[a-d][)\.]\s*")


def _clauses(text: str) -> list[str]:
    """Split an answer span into the clauses its slots live in.

    Cut on ``;``, on newlines, and on ``a)``/``b)`` enumerators.  The cut is
    what stops ``a) Not memoryless; b) Causal`` from letting the ``Not`` reach
    ``Causal`` -- the single most common answer form in the archive (11 of 24)
    and one a whole-span matcher inverts.
    """
    return [c for c, _ in _clauses_marked(text)]


def _clauses_marked(text: str) -> list[tuple[str, bool]]:
    """Clauses paired with whether each one *followed an enumerator*.

    The flag exists for the positional fallback in ``_slot_text``; see the
    comment there for the defect it closes.
    """
    out: list[tuple[str, bool]] = []
    pos, enumerated = 0, False
    for m in _ENUM_RE.finditer(text):
        for x in re.split(r"[;\n]", text[pos:m.start()]):
            if x.strip():
                out.append((x.strip(), enumerated))
        pos, enumerated = m.end(), True
    for x in re.split(r"[;\n]", text[pos:]):
        if x.strip():
            out.append((x.strip(), enumerated))
    return out


def _slot_text(text: str, slot: PropertySlot, pos: int, clauses: list[str]) -> str:
    """The clause that speaks about ``slot``.

    Located by the slot's own vocabulary when present, and **by position**
    otherwise.  Position is not a fallback of convenience: ``a) Yes, b) Yes``
    names neither property and is a complete answer, used by gemini on 1 of its
    6 traces.  Its slots are identified by the enumerator alone.
    """
    surfaces = list(slot.positive) + list(slot.negative) + [slot.name]
    owned = [c for c in clauses if _label_hit(c, surfaces)]
    if len(owned) == 1:
        return owned[0]
    if owned:
        # More than one clause mentions it (a restatement); the last is the
        # committed one, matching answer_span's rule.
        return owned[-1]
    # Positional fallback: count only the clauses that FOLLOWED an enumerator,
    # so a prose preface cannot occupy slot 0.  ``Although the algebra is
    # fiddly, a) Yes, b) Yes`` otherwise puts "Although the algebra is fiddly"
    # in the memoryless slot and the answer is refused -- 5 of the 351
    # recall-corpus cases, all on the one archived answer that uses the
    # unlabelled ``a) Yes, b) Yes`` form, and reached by any frame at all.
    # When nothing is enumerated, count every clause: that is the unframed
    # archive case and the flag must not change it.
    enumerated = [c for c, was_enum in _clauses_marked(text) if was_enum]
    pool = enumerated if len(enumerated) > pos else clauses
    if pos < len(pool):
        return pool[pos]
    return text


def _resolve_property(text: str, slot: PropertySlot) -> str | None:
    """Resolve one clause to ``"Yes"`` or ``"No"``.

    Precedence, and it is load-bearing:

    1. an explicit **declared negative** surface (``has memory``, ``non-causal``);
    2. a **negated positive** surface (``not memoryless``), including
       ``neither ... nor ...``;
    3. a **bare positive** surface (``causal``);
    4. a bare ``Yes`` / ``No``.

    The yes/no reading is *last* because it is the weakest evidence in a clause
    that also names the property: ``b) Causal (it depends only on the present
    input, not future ones)`` contains the word ``not``, and a yes/no-first
    reading would find no yes/no, fall through, and be fine -- but
    ``a) Not memoryless (it has memory)`` contains neither ``yes`` nor ``no``
    while ``Memoryless: No`` contains one, so the two orders disagree on real
    traces and only this one is right on both.
    """
    hit_neg = _label_hit(text, slot.negative) if slot.negative else None
    hit_pos = _label_hit(text, slot.positive)

    if hit_neg and (not hit_pos or hit_neg[0] > hit_pos[0]):
        return "No"
    if hit_pos:
        at, surface = hit_pos
        if _negated(text, at, surface):
            return "No"
        m = NEITHER_NOR_RE.search(text)
        if m and m.start() <= at <= m.end():
            return "No"
        # "Memoryless: No" -- the property is named and answered, not asserted.
        tail = text[at + len(surface):]
        if re.match(r"\s*[:=-]\s*", tail):
            if _NO_RE.match(tail.lstrip(" :=-")):
                return "No"
            if _YES_RE.match(tail.lstrip(" :=-")):
                return "Yes"
        return "Yes"
    if hit_neg:
        return "No"
    y, n = _YES_RE.search(text), _NO_RE.search(text)
    if y and n:
        return "Yes" if y.start() < n.start() else "No"
    if y:
        return "Yes"
    if n:
        return "No"
    return None


# ==========================================================================
# 3.  sequence
# ==========================================================================

_BRACED_RE = re.compile(r"[{\[]\s*([^{}\[\]]*?)\s*[}\]]")
_ORIGIN_STAR_RE = re.compile(r"\*\s*([-+]?\d+(?:\.\d+)?)\s*\*")
_INDEXED_RE = re.compile(r"\b([a-zA-Z])\s*\[\s*(-?\d+)\s*\]\s*=\s*([-+]?\d+(?:\.\d+)?)")
_FOR_N_RE = re.compile(
    r"for\s+n\s*(?:=|in|ranging\s+from)\s*[{\[(]?\s*(-?\d+)", re.IGNORECASE
)


class SeqAnswer:
    """A sequence answer: values, and the index of ``n = 0``.

    ``origin`` is the position within ``values`` that carries ``n = 0``, or
    ``None`` when the answer does not state one.  ``None`` is a real and
    frequent state -- gold itself omits the marker whenever the support does
    not contain ``n = 0`` -- and it is *not* the same as ``origin = 0``.
    """

    __slots__ = ("values", "origin", "form")

    def __init__(self, values: list[Decimal], origin: int | None, form: str):
        self.values, self.origin, self.form = values, origin, form

    def as_map(self) -> dict[int, Decimal] | None:
        if self.origin is None:
            return None
        return {i - self.origin: v for i, v in enumerate(self.values)}

    def trimmed_map(self) -> dict[int, Decimal] | None:
        """The index->value map with leading/trailing explicit zeros removed.

        Padding a sequence with zeros outside its support does not change the
        signal, and models do it constantly (claude signal 0 writes
        ``{*0*, -1, 5, -4, 3}`` where gold writes ``{-1, 5, -4, 3}``).  Interior
        zeros are structural and are never trimmed.
        """
        m = self.as_map()
        if m is None:
            return None
        keys = sorted(m)
        while keys and m[keys[0]] == 0:
            keys.pop(0)
        while keys and m[keys[-1]] == 0:
            keys.pop()
        return {k: m[k] for k in keys}

    def __repr__(self) -> str:  # pragma: no cover
        return f"Seq({self.values}, origin={self.origin}, {self.form})"


_FOR_N_LIST_RE = re.compile(
    r"for\s+n\s*(?:=|in|ranging\s+from)?\s*[{\[(]?\s*"
    r"(-?\d+(?:\s*,\s*-?\d+)*)",
    re.IGNORECASE,
)


def parse_sequence(text: str) -> SeqAnswer | None:
    """Parse a sequence answer in any of the observed presentations.

    Five packings appear in the 16 archived ``signal_operations`` traces, and
    **two traces use more than one at once** -- which is why this parses into a
    single index->value map and merges, rather than dispatching on a form:

    ==========================================  =============================  ==
    packing                                     example                        n
    ==========================================  =============================  ==
    braces, asterisk origin                     ``{2, 3, *-8*, 6}``             1
    braces, no origin stated                    ``{-7, 7, 5, -1}``              9
    braces + explicit index list                ``{...} for n = {-6,-5,...}``   2
    per-element assignments                     ``z[-1] = -8, z[-2] = -3``      1
    braces *and* per-element assignments        ``{6,0,0,0,0} for n={0..4}      1
                                                  and z[-1] = -8, ...``
    ==========================================  =============================  ==

    (Five further traces state a brace list inside LaTeX ``\\{ \\}``; that is a
    rendering of the second packing, not a sixth.)

    The mixed case is not a curiosity.  gemma's trace 1 states the n>=0 tail in
    braces and the n<0 head as assignments; read either half alone the answer is
    wrong, read together it is a one-place origin shift from gold -- which is
    the finding, and only the merged parse can see it.
    """
    t = prepare(text, keep_asterisks=True)
    parts: list[str] = []
    m: dict[int, Decimal] = {}

    # (a) per-element assignments, e.g. "z[-1] = -8, z[0] = 6".
    for _, i, v in _INDEXED_RE.findall(t):
        m[int(i)] = Decimal(v)
    if m:
        parts.append("indexed")

    # (b) the brace list, with an asterisk origin or an explicit index list.
    braced = _parse_braced(t)
    if braced is not None:
        values, star_pos, idx_list = braced
        if star_pos is not None:
            base = -star_pos
            parts.append("braced-star")
        elif idx_list is not None and len(idx_list) == len(values):
            base = idx_list[0]
            parts.append("braced-indexed")
        elif idx_list is not None:
            base = idx_list[0]
            parts.append("braced-indexed")
        else:
            base = None
            parts.append("braced")
        if base is not None:
            for off, v in enumerate(values):
                m.setdefault(base + off, v)
        elif not m:
            return SeqAnswer(values, None, "braced")

    if not m:
        return None
    keys = sorted(m)
    if keys[-1] - keys[0] + 1 != len(keys):
        # A gap: fill it with zeros, which is what an unstated index means.
        for k in range(keys[0], keys[-1] + 1):
            m.setdefault(k, Decimal(0))
        keys = sorted(m)
    origin = -keys[0] if keys[0] <= 0 <= keys[-1] else None
    return SeqAnswer([m[k] for k in keys], origin, "+".join(dict.fromkeys(parts)))


def _parse_braced(t: str) -> tuple[list[Decimal], int | None, list[int] | None] | None:
    """The value list, the asterisk position, and any explicit index list.

    Where two brace groups appear, the one carrying the **values** is the one
    the answer is about; an ``for n = {0, 1, 2, 3, 4}`` group is an index list,
    not an answer, and taking the last group blindly reads the indices as the
    answer.  The index group is identified by following a ``for n`` marker.
    """
    groups = [(m.start(), m.group(1)) for m in _BRACED_RE.finditer(t)]
    if not groups:
        return None

    idx_list: list[int] | None = None
    fm = _FOR_N_LIST_RE.search(t)
    if fm:
        try:
            idx_list = [int(x) for x in re.split(r"\s*,\s*", fm.group(1))]
        except ValueError:
            idx_list = None

    # Discard any group that is not a value list.
    #
    # Two ways a group is not one, and the second was a live defect: the
    # bracket alternative in `_BRACED_RE` also matches an **index subscript**,
    # so `y[1] = -1, y[2] = 5, y[3] = -4, y[4] = 3` contributed four "groups"
    # and `val_groups[-1]` picked `[4]` -- the parser read the last subscript as
    # the whole answer.  It silently discarded the real brace group's origin,
    # which is why archived signal trace 0 and the identical D4.4 case seq-08
    # disagreed with each other.  A bracket group is a sequence only if it
    # contains a separator; a brace group always is.
    val_groups = []
    for start, body in groups:
        if not _split_items(body):
            continue
        if idx_list is not None and fm is not None and fm.start() <= start <= fm.end():
            continue
        if t[start] == "[" and "," not in body:
            continue  # a subscript, not a one-element sequence
        val_groups.append(body)
    if not val_groups:
        return None
    body = val_groups[-1]

    star_pos = None
    values: list[Decimal] = []
    for pos, it in enumerate(_split_items(body)):
        star = _ORIGIN_STAR_RE.fullmatch(it.strip())
        if star:
            star_pos = pos
            it = star.group(1)
        it = it.replace("*", "").strip()
        try:
            values.append(Decimal(it))
        except InvalidOperation:
            return None
    if not values:
        return None
    if idx_list is not None and len(idx_list) not in (len(values), 1):
        # An index list of a different length than the values does not pin
        # anything; treat the sequence as origin-less rather than guessing.
        idx_list = None
    return values, star_pos, idx_list


def _split_items(body: str) -> list[str]:
    return [p for p in (x.strip() for x in body.split(",")) if p != ""]


def compare_sequence(
    gold: str,
    candidate: str,
    *,
    require_origin: bool = True,
    extract: bool = True,
) -> Verdict:
    """Sequence-with-origin equality: the values **and** the ``n = 0`` index.

    The origin is part of the answer.  **The archive does not prove it, and the
    first version of this docstring said it did** (Reviewer E, F1).

    That claim read: *"signal trace 8 emits exactly gold's multiset of values in
    exactly gold's order and is wrong … it is 1 of the 16 archived traces."*
    Both halves are false.  Trace 8's answer is a **rotation** of gold's value
    list, not a reordering-preserving one -- ``ground_truth.py``'s own label for
    that trace says "is a rotation of the truth" -- so a values-only comparator
    rejects it too.  Measured: **0 of 16 archived traces reach the
    origin-discriminating branch at all.**

    So ``require_origin`` rests on the *item*, not on observed model error: the
    question states ``x[n]`` with its origin marked, the transformation maps
    indices, and a sequence at the wrong indices is a different signal.  That is
    a sound argument and it is the only one there is.  Like the symbolic rules
    it is a **route-3** rule -- specified and adversarially exercised, not
    validated against real output -- and D4.2 §2.4 lists it there.

    Presentation is not part of the answer: braces or brackets, asterisk or an
    explicit index list, LaTeX-escaped or bare, zero-padded or trimmed all
    compare equal.
    """
    k = "sequence"
    g_txt = answer_span(gold)[0] if extract else gold
    c_txt = answer_span(candidate)[0] if extract else candidate

    g = parse_sequence(g_txt)
    if g is None:
        return unresolved(k, "no sequence in the gold answer span", gold_canonical=g_txt)
    c = parse_sequence(c_txt)
    if c is None:
        return unresolved(
            k, "no sequence in the candidate answer span",
            gold_canonical=repr(g), cand_canonical=c_txt.strip()[:120],
        )

    obs = []
    if g.form != c.form:
        obs.append(f"presentation differs: gold {g.form}, candidate {c.form}")

    if not require_origin or (g.origin is None and c.origin is None):
        # Gold states no origin: the support does not contain n = 0, so the
        # answer is the value list alone and is compared as such.
        if _trim(g.values) == _trim(c.values):
            return match(k, gold_canonical=_trim(g.values), cand_canonical=_trim(c.values),
                         observations=obs + ["origin not stated by gold; values compared"])
        return mismatch(k, "values differ (no origin stated by gold)",
                        gold_canonical=_trim(g.values), cand_canonical=_trim(c.values),
                        observations=obs)

    if g.origin is not None and c.origin is None:
        # The candidate did not say where n = 0 is.  Before calling that
        # undetermined, ask whether *any* placement would match: if none does,
        # the values are wrong whatever the origin and the case is decidable.
        # This matters -- it converts 3 of the 16 archived traces from
        # UNRESOLVED to a MISMATCH we can actually defend, and it never
        # converts one the other way.
        gm = g.trimmed_map()
        cv = _trim(c.values)
        if not any(
            {i + shift: v for i, v in enumerate(cv)} == gm
            for shift in range(-len(cv) - 8, 9)
        ):
            return mismatch(
                k, "no placement of n=0 makes the candidate's values equal gold's",
                gold_canonical=_fmt_map(gm), cand_canonical=str(cv), observations=obs,
            )
        return unresolved(
            k, "gold pins n=0, the candidate states no origin, and its values would "
               "match at some placement; the answer is undetermined",
            gold_canonical=_fmt_map(gm), cand_canonical=str(cv), observations=obs,
        )
    if g.origin is None and c.origin is not None:
        # **Gold cannot state an origin here and the candidate has stated one,
        # so the comparison is undecidable -- not a match.**
        #
        # The first version compared trimmed values and returned MATCH with
        # "candidate stated an origin gold does not" as an observation.
        # Reviewer E (F3) showed what that accepts: on the 16.4% of instances
        # where gold's support excludes n = 0 (D-050), `require_origin=True` was
        # silently inert, and three candidate forms taken from real archived
        # traces -- an asterisk at the wrong element, an index list four places
        # off, an off-by-one per-element listing -- were all MATCHed against a
        # real archived gold.  Those are wrong signals, not under-specified
        # ones, and the observation string is no defence because it rides on a
        # MATCH and is invisible in precision.
        #
        # UNRESOLVED is the honest outcome: the item cannot decide it.  The
        # values are still compared first, because if they differ the case is
        # decidable regardless of any origin.
        if _trim(g.values) != _trim(c.values):
            return mismatch(k, "values differ", gold_canonical=_trim(g.values),
                            cand_canonical=_trim(c.values), observations=obs)
        return unresolved(
            k, "gold's support excludes n=0 so gold states no origin, and the "
               "candidate pins one; the values agree but the placement cannot be checked",
            gold_canonical=_trim(g.values), cand_canonical=_fmt_map(c.trimmed_map()),
            observations=obs,
        )

    gm, cm = g.trimmed_map(), c.trimmed_map()
    if gm == cm:
        return match(k, gold_canonical=_fmt_map(gm), cand_canonical=_fmt_map(cm),
                     observations=obs)
    if sorted(gm.values()) == sorted(cm.values()) and list(gm.values()) == list(cm.values()):
        return mismatch(
            k, f"values agree but the origin differs: gold n=0 at {g.origin}, candidate at {c.origin}",
            gold_canonical=_fmt_map(gm), cand_canonical=_fmt_map(cm), observations=obs,
        )
    return mismatch(k, "index->value maps differ",
                    gold_canonical=_fmt_map(gm), cand_canonical=_fmt_map(cm), observations=obs)


def _trim(vals: list[Decimal]) -> list[Decimal]:
    v = list(vals)
    while v and v[0] == 0:
        v.pop(0)
    while v and v[-1] == 0:
        v.pop()
    return v


def _fmt_map(m: dict[int, Decimal] | None) -> str:
    if m is None:
        return "None"
    return "{" + ", ".join(f"{k}:{v}" for k, v in sorted(m.items())) + "}"


# ==========================================================================
# 4.  symbolic
# ==========================================================================

_ARBITRARY_RE = re.compile(r"\b([A-Za-z])\s*\(\s*([a-z])\s*\)")
_BARE_CONST_RE = re.compile(r"(?<![\w(])\+\s*([CK])\b(?!\s*\()")


def compare_symbolic(
    gold: str,
    candidate: str,
    *,
    symbols: Sequence[str] = ("x", "y", "z"),
    allow_arbitrary: bool = True,
    timeout_s: float = 5.0,
    extract: bool = True,
) -> Verdict:
    """Symbolic equivalence by CAS, with a stated policy for the arbitrary
    function of integration.

    **The arbitrary-function policy.**  The item asks for the *simplest*
    expression, and gold sets the arbitrary function to zero.  A candidate that
    carries it explicitly -- ``3x^2/2 - 6xy + C(y)`` -- has answered the
    question it was asked and has additionally shown it understands why the
    function is there.  It is accepted, with the retained function recorded as
    an observation.  This is not a lenience: the general solution and the
    simplest particular solution differ by exactly that term, and rejecting it
    would mark a *more* complete answer wrong.  Archived support: mathstral's
    continuity trace, 1 of 6.

    A candidate that leaves a *bare* constant (``+ C``) instead is a different
    case and is also accepted, because the same reasoning applies one dimension
    down.  A candidate whose answer is *only* a constant (wizardmath's ``u(x,y)
    = C``) is a mismatch: it has dropped the whole solution, not annotated it.

    **Failure behaviour is UNRESOLVED, never a match.**  A CAS timeout, a parse
    failure, or an expression carrying a symbol outside the declared alphabet
    all return UNRESOLVED.  Nothing about a comparator failing to parse an
    answer is evidence the answer is right, and a ``sympy`` exception must not
    become a silent pass -- that is D3.4 |S|8B.10's rule for this kind.
    """
    k = "symbolic"
    g_raw = answer_span(gold)[0] if extract else gold
    c_raw = answer_span(candidate)[0] if extract else candidate

    g_expr_txt = _isolate_expression(prepare(g_raw))
    c_expr_txt = _isolate_expression(prepare(c_raw))
    if g_expr_txt is None:
        return unresolved(k, "no expression in the gold answer span", gold_canonical=g_raw[:120])
    if c_expr_txt is None:
        return unresolved(k, "no expression in the candidate answer span",
                          gold_canonical=g_expr_txt, cand_canonical=c_raw.strip()[:120])

    obs: list[str] = []
    c_expr_txt, arb = _drop_arbitrary(c_expr_txt)
    if arb:
        if not allow_arbitrary:
            return mismatch(k, f"candidate retains the arbitrary term {arb}",
                            gold_canonical=g_expr_txt, cand_canonical=c_expr_txt)
        obs.append(f"candidate retained the arbitrary integration term {arb} (accepted)")
    g_expr_txt, _ = _drop_arbitrary(g_expr_txt)

    if not c_expr_txt.strip():
        return mismatch(k, "candidate's answer is the arbitrary term alone",
                        gold_canonical=g_expr_txt, cand_canonical=arb or "")

    # --- Stage 1: the polynomial fragment, decided exactly and locally. ---
    #
    # Every answer this kind sees for `incompressible_continuity` is a
    # bivariate polynomial with rational coefficients and degree <= 2, and
    # equality there is decidable by expanding to a coefficient map -- no CAS,
    # no timeout, no dependency, and no UNRESOLVED outcome.  That matters more
    # than the convenience: this is the template with **six** archived traces,
    # and a comparator whose answer depends on whether an optional package is
    # installed would make the thinnest evidence in the phase thinner still.
    gp, cp = _as_polynomial(g_expr_txt, symbols), _as_polynomial(c_expr_txt, symbols)
    if gp is not None and cp is not None:
        obs.append("decided in the polynomial fragment (no CAS)")
        if gp == cp:
            return match(k, gold_canonical=_fmt_poly(gp), cand_canonical=_fmt_poly(cp),
                         observations=obs)
        return mismatch(k, f"coefficients differ: {_poly_diff(gp, cp)}",
                        gold_canonical=_fmt_poly(gp), cand_canonical=_fmt_poly(cp),
                        observations=obs)

    # --- Stage 2: outside the fragment, a CAS is required. ---
    #
    # The other eight `symbolic` templates in the corpus carry sinc, exp and
    # Q-functions, so the fragment is not enough in general.
    try:
        import sympy  # noqa: PLC0415  -- optional dependency, imported at use
    except ImportError:
        return unresolved(
            k, "expression is outside the polynomial fragment and sympy is not installed",
            gold_canonical=g_expr_txt, cand_canonical=c_expr_txt, observations=obs)

    local = {s: sympy.Symbol(s) for s in symbols}
    try:
        ge = sympy.sympify(_to_sympy(g_expr_txt), locals=local)
        ce = sympy.sympify(_to_sympy(c_expr_txt), locals=local)
    except Exception as exc:  # sympy raises a wide family
        return unresolved(k, f"could not parse an expression: {type(exc).__name__}",
                          gold_canonical=g_expr_txt, cand_canonical=c_expr_txt,
                          observations=obs)

    # **A parse that did not yield an expression is UNRESOLVED, not a crash**
    # (D5.7b, N4).  `sympify` returns a *tuple* for a span containing a comma,
    # and `continuous_to_discrete_conversion`'s answer has one -- so
    # `.free_symbols` raised AttributeError on all 132 of its pairs, the
    # template contributed zero verdicts, and it satisfied "zero false accepts"
    # by crashing.  A template that raises has not passed.
    if not isinstance(ge, sympy.Basic) or not isinstance(ce, sympy.Basic):
        return unresolved(
            k, "expression did not parse to a single expression "
               f"(gold {type(ge).__name__}, candidate {type(ce).__name__})",
            gold_canonical=str(ge)[:80], cand_canonical=str(ce)[:80], observations=obs)

    extra = (ge.free_symbols | ce.free_symbols) - set(local.values())
    if extra:
        return unresolved(
            k, f"symbol(s) outside the declared alphabet: {sorted(str(s) for s in extra)}",
            gold_canonical=str(ge), cand_canonical=str(ce), observations=obs,
        )

    try:
        diff = sympy.simplify(ge - ce)
    except Exception as exc:
        return unresolved(k, f"CAS failed to decide: {type(exc).__name__}",
                          gold_canonical=str(ge), cand_canonical=str(ce), observations=obs)

    if diff == 0:
        return match(k, gold_canonical=str(ge), cand_canonical=str(ce), observations=obs)
    return mismatch(k, f"gold - candidate simplifies to {diff}, not 0",
                    gold_canonical=str(ge), cand_canonical=str(ce), observations=obs)


#: A term of a bivariate polynomial: an optional rational coefficient followed
#: by zero or more ``symbol`` or ``symbol^k`` factors.
_POLY_TERM_RE = re.compile(
    r"(?P<sign>[-+]?)\s*"
    r"(?P<coeff>\d+(?:\.\d+)?)?\s*\*?\s*"
    r"(?P<vars>(?:[a-zA-Z](?:\s*\^\s*\d+)?\s*\*?\s*)*)"
    r"(?:/\s*(?P<den>\d+(?:\.\d+)?))?$"
)


def _as_polynomial(expr: str, symbols: Sequence[str]) -> dict[tuple[int, ...], Fraction] | None:
    """Expand ``expr`` to ``{exponent tuple: coefficient}``, or None if it is
    outside the fragment.

    Handles a sum of products of powers with rational coefficients, including
    the ``((3x^2)/(2))`` form ``\\frac`` unwrapping produces.  Anything else --
    a function call, a nested parenthesis that is not a bare fraction, an
    unknown symbol -- returns None so the caller falls through to the CAS
    rather than guessing.
    """
    e = expr.replace(" ", "")
    if not e:
        return None
    # ((a)/(b)) -> a/b, one level, which is all \frac produces here.
    for _ in range(3):
        new = re.sub(r"\(\(([^()]*)\)/\(([^()]*)\)\)", r"\1/\2", e)
        if new == e:
            break
        e = new
    e = re.sub(r"\(([^()]*)\)", r"\1", e)
    if "(" in e or ")" in e:
        return None

    # Split into signed terms at top level.
    terms = re.findall(r"[-+]?[^-+]+", e)
    out: dict[tuple[int, ...], Fraction] = {}
    idx = {s: i for i, s in enumerate(symbols)}
    for t in terms:
        m = _POLY_TERM_RE.fullmatch(t)
        if not m:
            return None
        sign = -1 if m.group("sign") == "-" else 1
        raw, den = m.group("coeff"), m.group("den")
        try:
            coeff = Fraction(1) if raw is None else Fraction(Decimal(raw))
            if den is not None:
                # `3x^2/2` -- the divisor trails the variables, which is what
                # `\frac{3x^2}{2}` unwraps to.  mathstral writes the one
                # archived continuity answer this way.
                coeff /= Fraction(Decimal(den))
        except (InvalidOperation, ZeroDivisionError):
            return None
        exps = [0] * len(symbols)
        vs = m.group("vars") or ""
        for vm in re.finditer(r"([a-zA-Z])(?:\^(\d+))?", vs):
            name = vm.group(1)
            if name not in idx:
                return None  # symbol outside the declared alphabet
            exps[idx[name]] += int(vm.group(2) or 1)
        key = tuple(exps)
        out[key] = out.get(key, Fraction(0)) + sign * coeff
    nonzero = {k: v for k, v in out.items() if v != 0}
    # **An empty map is not a polynomial; it is a failed parse** (D5.7, N2).
    # Returning `{}` made `{} == {}` a MATCH, and on
    # `autocorrelation_rect_pulse` -- where `_isolate_expression` returns the
    # bare string "0" -- that was 128 of 132 gold-by-gold pairs matching each
    # other.  Every expression the parser failed on matched every other one it
    # failed on.  D4.1 S4.4's own rule is "failure is UNRESOLVED, never MATCH",
    # violated by a path that never reaches the CAS.  Returning None sends the
    # caller to the CAS, which decides it or says it cannot.
    return nonzero or None


def _fmt_poly(p: dict[tuple[int, ...], Fraction]) -> str:
    return "{" + ", ".join(f"{k}:{v}" for k, v in sorted(p.items())) + "}"


def _poly_diff(a: dict, b: dict) -> str:
    keys = sorted(set(a) | set(b))
    return ", ".join(
        f"{k}: {a.get(k, 0)} vs {b.get(k, 0)}" for k in keys if a.get(k, 0) != b.get(k, 0)
    )


_LHS_ONLY_RE = re.compile(r"^[A-Za-z](?:\s*\(\s*[a-z]\s*(?:,\s*[a-z]\s*)?\))?$")

#: A left-hand side that looks like an ASSIGNMENT rather than prose:
#: ``R_g(tau)``, ``u(x, y)``, ``y_c``.  Not ``This triangle has a peak value
#: of 512 at tau``.
_ASSIGN_LHS_RE = re.compile(
    r"^[A-Za-z][A-Za-z0-9_]{0,14}"
    r"(?:\s*\(\s*[a-z][a-z0-9_]*(?:\s*,\s*[a-z][a-z0-9_]*)*\s*\))?$")

#: A conditional/piecewise construction.  An answer written
#: ``64*(8 - |tau|), for |tau| <= 8, and 0 otherwise`` is a *piecewise*
#: function, and a single-expression comparator cannot represent it.  It must
#: say so rather than extract one branch.
_PIECEWISE_RE = re.compile(r"\b(?:otherwise|elsewhere)\b|,\s*for\b", re.IGNORECASE)

#: Split on a bare ``=`` only -- never inside ``<=``, ``>=``, ``!=``, ``==``.
_BARE_EQ_RE = re.compile(r"(?<![<>!=])=(?!=)")


def _isolate_expression(text: str) -> str | None:
    """Take the expression stated by the span.

    Models write ``The simplest expression ... is u = <expr>`` (the last ``=``
    is right) and also ``u(x, y) = 10y^2 + 4xy = C(x)`` (the last ``=`` is
    *wrong* -- it yields the arbitrary function and loses the answer).  So the
    rule is neither "first" nor "last": walk the segments after the first
    ``=`` and take the first that is not a bare left-hand side.

    Both forms occur in the six archived continuity traces, so a rule that
    handles only one of them is wrong on the template with the least evidence.
    """
    # Pick the *line* that states the expression, not blindly the first: models
    # write a sentence of preamble and then the equation on the next line
    # ("... is given by:\nu(x, y) = 10y^2 + 4xy").  Taking line 0 hands the
    # prose to the parser, which then cannot decide a plainly wrong answer.
    lines = [ln.strip() for ln in text.strip().split("\n") if ln.strip()]
    if not lines:
        return None

    # **Prefer a line that ASSIGNS to a symbol** (D5.7, N2).  Taking the last
    # `=`-bearing line picks up trailing prose: on
    # `autocorrelation_rect_pulse` the span ends
    # "This triangle has a peak value of 512 at tau = 0", so the isolated
    # "expression" was the string `0` -- for every instance.  `_as_polynomial`
    # then parsed `0` to an empty coefficient map on both sides and `{} == {}`
    # was a MATCH: 128 of 132 gold-by-gold pairs.
    #
    # This is the THIRD positional rule this function has had, and the first two
    # are recorded as errors in `phase4_summary.md` S8 ("split on the last `=`",
    # then "`.split(chr(10))[0]`").  So the rule is not "a different position" --
    # it is a structural property: an answer line assigns to a symbol, and prose
    # does not.
    eq_lines = [ln for ln in lines if _BARE_EQ_RE.search(ln)]
    assigns = [ln for ln in eq_lines
               if _ASSIGN_LHS_RE.fullmatch(_BARE_EQ_RE.split(ln, 1)[0].strip())]
    text = (assigns[-1] if assigns else
            (eq_lines[-1] if eq_lines else lines[-1])).rstrip(".").strip()

    # A piecewise answer is not a single expression, and picking one branch of
    # it is exactly the "unrecoverable parse presented as a recoverable one"
    # D4.1 S4.4 forbids.  Refuse instead.
    if _PIECEWISE_RE.search(text):
        return None

    parts = [p.strip().rstrip(".").strip() for p in _BARE_EQ_RE.split(text)]
    if len(parts) == 1:
        cand = parts[0]
    else:
        cand = ""
        for seg in parts[1:]:
            seg = re.sub(r"^\s*(?:u|v|w|f|g|y|z)\s*\([^)]*\)\s*", "", seg).strip()
            if seg and not _LHS_ONLY_RE.fullmatch(seg):
                cand = seg
                break
        if not cand:
            cand = parts[-1]
    cand = re.sub(r"^\s*(?:u|v|w|f|g|y|z)\s*\([^)]*\)\s*", "", cand).strip()
    return cand or None


def _drop_arbitrary(expr: str) -> tuple[str, str | None]:
    """Remove ``+ f(x)`` / ``+ C`` and report what was removed."""
    found = None
    # The whole expression is a lone arbitrary constant or function -- the
    # candidate has written down the *form* of the answer and none of it.
    # wizardmath's ``u(x, y) = C`` is this case; it must reach the "arbitrary
    # term alone" mismatch rather than fall through to the CAS and become an
    # UNRESOLVED, which would hide a plainly wrong answer behind an
    # unrecognised symbol.
    if re.fullmatch(r"[+\-]?\s*[A-Za-z](?:\s*\(\s*[a-z]\s*\))?", expr.strip()):
        return "", expr.strip()
    m = _ARBITRARY_RE.search(expr)
    if m:
        found = m.group(0)
        expr = (expr[:m.start()] + expr[m.end():])
        expr = re.sub(r"[+\-]\s*$", "", expr.strip())
    else:
        m2 = _BARE_CONST_RE.search(expr)
        if m2:
            found = m2.group(0)
            expr = expr[:m2.start()] + expr[m2.end():]
    return expr.strip().rstrip("+- ").strip(), found


def _to_sympy(expr: str) -> str:
    """``3x^2 - 6xy`` -> ``3*x**2 - 6*x*y``."""
    e = expr.replace("^", "**")
    e = re.sub(r"(\d)\s*\(", r"\1*(", e)
    e = re.sub(r"\)\s*\(", r")*(", e)
    e = re.sub(r"(\d)\s*([a-zA-Z])", r"\1*\2", e)
    e = re.sub(r"([a-zA-Z])\s*(?=[a-zA-Z])", r"\1*", e)
    e = re.sub(r"\)\s*([a-zA-Z0-9])", r")*\1", e)
    e = re.sub(r"([a-zA-Z0-9])\s*\(", r"\1*(", e)
    e = re.sub(r"\*{3,}", "**", e)
    return e


# ==========================================================================
# 5.  narrative
# ==========================================================================


def compare_narrative(gold: str, candidate: str, **_: Any) -> Verdict:
    """Always UNRESOLVED.  This is the kind's entire specification.

    ``narrative`` names a milestone whose content is prose that no comparator
    in this document can decide -- a justification, an interpretation, a
    physical explanation.  Specifying it as "compare with an LLM" would put the
    AI Tribunal back on the critical path for exactly the milestones where it is
    least accountable, which is the thing this phase exists to remove.

    So the kind is specified as *undecidable by construction*: declaring a
    milestone ``narrative`` **removes it from the automatic score** and routes
    it to whatever human or model process the benchmark chooses, with that
    routing visible in the results rather than hidden inside an accuracy number.
    The gate on this kind is therefore not precision or recall but **census**:
    the fraction of milestones declared ``narrative`` is a reported quantity,
    and a rising fraction is a benchmark quietly returning to LLM judging.
    """
    return unresolved(
        "narrative",
        "narrative milestones are not machine-decidable and are excluded from the automatic score",
        gold_canonical=prepare(answer_span(gold)[0])[:200],
        cand_canonical=prepare(answer_span(candidate)[0])[:200],
    )


# ==========================================================================
# 6.  check
# ==========================================================================

_TRUE_SURFACES = ("yes", "true", "satisfied", "holds", "valid", "acceptable",
                  "is met", "passes", "confirmed", "adequate", "safe")
_FALSE_SURFACES = ("no", "false", "not satisfied", "does not hold", "invalid",
                   "unacceptable", "is not met", "fails", "not confirmed",
                   "inadequate", "unsafe")


def compare_check(
    gold: str,
    candidate: str,
    *,
    quantity_precision: int | None = None,
    extract: bool = True,
) -> Verdict:
    """A verdict-bearing assertion: a boolean, and the quantity it rests on.

    ``check`` is the kind for "verify that the beam is safe", "confirm the flow
    is laminar", "is the design adequate?".  Its answer has two parts and the
    comparator scores **both**: the boolean, and -- when gold states one -- the
    compared quantity, at |S|7.1 tolerance.

    Scoring only the boolean would make the kind a coin flip: a candidate that
    computes a deflection of 40 mm against a 25 mm limit and concludes
    "not acceptable" agrees with a gold that computed 30 mm, and a comparator
    that credits it has credited reasoning that did not happen on a two-valued
    answer where guessing is free.  This is the kind where the false-accept
    asymmetry of Decision 3 bites hardest, and the quantity is what removes it.
    """
    k = "check"
    g_txt = prepare(answer_span(gold)[0] if extract else gold)
    c_txt = prepare(answer_span(candidate)[0] if extract else candidate)

    g_bool = _resolve_bool(g_txt)
    if g_bool is None:
        return unresolved(k, "gold states no verdict", gold_canonical=g_txt)
    commit = find_commitment(c_txt, _TRUE_SURFACES + _FALSE_SURFACES)
    if not commit:
        return unresolved(k, f"candidate does not commit to a verdict: {commit.reason}",
                          gold_canonical=str(g_bool), cand_canonical=c_txt)
    c_bool = _resolve_bool(commit.clause)
    if c_bool is None:
        return unresolved(k, "candidate's committed clause states no verdict",
                          gold_canonical=str(g_bool), cand_canonical=commit.clause)

    obs: list[str] = []
    g_q, c_q = parse_number(g_txt, "check"), parse_number(c_txt, "check")
    if g_q is not None:
        if c_q is None:
            return unresolved(
                k, "gold's verdict rests on a stated quantity the candidate does not state",
                gold_canonical=f"{g_bool} @ {g_q}", cand_canonical=str(c_bool),
            )
        num = compare_numeric(str(g_q), str(c_q), precision=quantity_precision
                              if quantity_precision is not None else _decimals(g_txt, "check"),
                              extract=False)
        if not num.is_match:
            return mismatch(
                k, f"verdict {c_bool} rests on a different quantity: {num.reason}",
                gold_canonical=f"{g_bool} @ {g_q}", cand_canonical=f"{c_bool} @ {c_q}",
            )
        obs.append(f"supporting quantity agrees ({c_q})")

    if c_bool == g_bool:
        return match(k, gold_canonical=str(g_bool), cand_canonical=str(c_bool), observations=obs)
    return mismatch(k, f"verdict {c_bool} != gold {g_bool}",
                    gold_canonical=str(g_bool), cand_canonical=str(c_bool), observations=obs)


def _resolve_bool(text: str) -> bool | None:
    hit_f = _label_hit(text, _FALSE_SURFACES)
    hit_t = _label_hit(text, _TRUE_SURFACES)
    if hit_t and _negated(text, hit_t[0], hit_t[1]):
        hit_f = hit_f or hit_t
        hit_t = None
    if hit_t and hit_f:
        return hit_t[0] > hit_f[0]
    if hit_t:
        return True
    if hit_f:
        return False
    return None


COMPARATORS = {
    "numeric": compare_numeric,
    "categorical": compare_categorical,
    "sequence": compare_sequence,
    "symbolic": compare_symbolic,
    "narrative": compare_narrative,
    "check": compare_check,
}
