"""Shared emission helpers for template gold text (Phase 6, D6.7).

WHY THIS MODULE EXISTS, before any of its contents.

A hard-coded ``+`` in front of an interpolated signed value prints
``876*t + -86.8 deg`` and ``30.5 + j-51.22``. Neither is a number in any
notation, and the second is the malformed-complex defect the audit recorded.
The sign belongs to the operator; the coefficient carries only its magnitude.

Phase 5 fixed this inside ``waves_and_phasors.py`` with two local helpers, and
``phase5_contract_scan``'s ``doubled_sign_corpus`` census still reports **14
templates** emitting a doubled sign - spread across SEVEN files in THREE
branches (electrical, mechanical, chemical). The spec's instruction (D6.7,
Reviewer A |S|5.6) is to promote the helpers **before** fixing any of them,
and the reason is worth stating plainly: nine hand-written sign fixes is nine
fresh chances to write ``+ {value}`` again. One helper, used everywhere, is
one place to be right.

WHAT IS PRESERVED. Both functions are promoted with their behaviour unchanged -
this module is a move, not a rewrite, so that the promotion itself cannot alter
emitted text. ``waves_and_phasors.py`` keeps working because it imports these
rather than defining its own; the five electrical, two mechanical and two
chemical files adopt them in turn.

The printed magnitude is ``abs()`` of the *same* rounded value the template
consumes, so Phase 2's round-then-print property is unchanged: only the
position of the minus sign moves.

DISCOVERY. This file sits at ``data/templates/branches/_emission.py`` - four
path segments. ``tests/template_integrity/core.py:discover()`` requires at
least six (``.../branches/<branch>/<area>/<file>.py``), so it is not picked up
as a template module, and it defines no ``template_*`` function in any case.
"""


def signed_term(value, unit="", fmt=str):
    """Render ``value`` as an operator followed by its magnitude.

    ``signed_term(-86.8, "deg")`` -> ``"- 86.8 deg"``;
    ``signed_term(12.5)``         -> ``"+ 12.5"``.

    The caller writes ``f"{omega}*t {signed_term(phi, 'deg')}"`` - note the
    SPACE and no ``+``. Writing ``f"... + {signed_term(...)}"`` reintroduces
    the very defect this exists to remove.

    ``fmt`` formats the MAGNITUDE, for sites whose original text carried a
    format spec. ``{B_vec[0]:.2e}`` must stay in scientific notation, and the
    default ``str`` would silently flatten it to ``4.5e-05`` -> ``4.5e-05``
    only by luck and to ``0.000045`` in general::

        signed_term(-4.5e-05, "T", lambda v: f"{v:.2e}")  -> "- 4.50e-05 T"

    The default is ``str``, which is what an f-string already does for int and
    float, so every existing call renders byte-identically.

    ``unit`` is joined with a SPACE. Where the original text had none -
    ``{T1_C}°C`` - pass no unit and append it at the call site
    (``f"{signed_term(v)}°C"``), or the fix would change spacing that was never
    part of the defect.
    """
    op = "-" if value < 0 else "+"
    return f"{op} {fmt(abs(value))} {unit}".rstrip()


def paren_neg(value, fmt=str):
    """``(-5)`` when negative, ``5`` when positive - for an operand under an operator.

    The FOURTH shape, and the one where ``signed_term`` would be actively wrong.

    ``sqrt({a}^2 + {b}^2 + {c}^2)`` prints ``sqrt(1^2 + -5^2 + -6^2)`` when a
    component is negative. The doubled sign is real, but the ``+`` here is a
    genuine ADDITION OPERATOR, not a term's sign - so re-signing it gives
    ``sqrt(1^2 - 5^2 - 6^2)``, which is a different and FALSE formula: a sum of
    squares becomes a subtraction of them.

    Worse, the detector reports that false version as CLEAN, because the doubled
    sign is indeed gone. Verified. So the scan cannot be the only check on this
    fix, and this helper exists to make the correct one easy:
    ``sqrt(1^2 + (-5)^2 + (-6)^2)``.

    Positives are returned untouched, so a fix using this changes only the
    instances that were actually defective - the same no-op-for-positives
    property the other helpers have.

    It also repairs something the doubled-sign detector structurally CANNOT see:
    a negative FIRST operand prints ``sqrt(-1^2 + ...)``, which reads as
    -(1**2) rather than (-1)**2. No doubled sign, and still wrong.
    """
    s = fmt(value)
    return f"({s})" if value < 0 else s


def joined_terms(values, fmt):
    """``a - b + c`` from signed ``values``, each operator taken from its term.

    The SECOND of the three doubled-sign shapes, and the one ``signed_term``
    cannot reach. Several templates build a sum with
    ``" + ".join(f"{v:.2f}" for v in values)``, which prints ``-2.723 + -0.444``
    the moment any term after the first is negative.

    ``fmt`` formats a MAGNITUDE and is supplied by the caller, so each template
    keeps its own precision - this helper never decides how a number is
    rounded, only where its sign goes.

    The FIRST term keeps its own sign, because it is a value and not an
    operand: ``joined_terms([-2.723, -0.444, 2.171], lambda v: f"{v:.3f}")``
    -> ``"-2.723 - 0.444 + 2.171"``.

    Not every ``" + ".join`` is a defect. Where the joined parts are whole
    expressions that always lead positive - ``" + ".join(signal_terms)`` over
    ``A * cos(...)`` terms - the join is correct and must be left alone.
    """
    values = list(values)
    if not values:
        return ""
    out = fmt(values[0])
    for v in values[1:]:
        op = "-" if v < 0 else "+"
        out += f" {op} {fmt(abs(v))}"
    return out


def rect_str(real, imag, precision):
    """Rectangular form ``x + jy`` with the sign OUTSIDE the ``j``.

    ``rect_str(30.5, -51.22, 2)`` -> ``"30.5 - j51.22"``, never
    ``"30.5 + j-51.22"``.
    """
    re_v = round(real, precision)
    im_v = round(imag, precision)
    op = "-" if im_v < 0 else "+"
    return f"{re_v} {op} j{abs(im_v)}"


#: The pre-promotion names, kept as aliases so the move is behaviour-preserving
#: for any caller still using them. New call sites use the public names.
_signed_term = signed_term
_rect_str = rect_str
