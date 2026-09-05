"""Independent round-trip oracle for `template_composite_shafts_series`.

Derived from the angle-of-twist relation for a circular shaft in torsion (Beer
& Johnston, *Mechanics of Materials*, ch. 3 - angle of twist in the elastic
range; the relation the template's docstring cites), **not** from the
template's code:

    phi_i = T * L_i / (J_i * G_i)        for each segment
    J_i   = (pi/2) * c_i^4               solid circular section
    phi_total = sum_i phi_i              segments in series carry the same T

The two segments are in series from the fixed end A to the free end C, so the
same applied torque passes through both and the twists add. That is the whole
content of the item.

Independence
------------
The oracle works in **base SI units** - N.m, metres, pascals - from the four
stated quantities per segment (length, diameter, shear modulus, and the shared
applied torque). It forms J and phi in one closed expression per segment and
adds them at full precision.

The trace instead states each J to six significant figures and each segment
angle to five, then sums the STATED angles - which is P3-compliant, since each
operand is printed before it is consumed, and it is what makes the printed
Step-4 line close. Its defect was the version before that: the segments were
printed at 5 dp and their total at 4 dp, rounded from the unrounded sum, so the
printed line disagreed with its own operands verbatim on 89.9% of instances and
failed closure on 4.85% (phase0_baseline.md claim 6).

Because this oracle does not replicate those stated intermediates, the two
agree only to within them - see TOLERANCE. Not replicating them is deliberate:
an oracle that reproduced the template's rounding sequence would agree with it
whether the sequence were right or wrong.

The diameter is stated in millimetres and the radius is half of it; reading the
stated diameter as a radius is a factor-of-16 error in J, which this oracle
would catch immediately.
"""
from __future__ import annotations

import math
import re
from decimal import Decimal, ROUND_HALF_UP

TEMPLATE_ID = "template_composite_shafts_series"

# Relative tolerance.
#
# `recompute` quantises to the five significant figures the trace quotes the
# answer to. Two effects then have to fit inside the tolerance:
#
# 1. Display. Half a unit in the last of five significant figures is at most
#    5e-5 relative.
#
# 2. The stated J values. The trace prints each polar second moment as
#    "{J:.5e}" and computes each segment angle from the printed value; since
#    phi ~ 1/J that is ~5e-7 relative per segment. (At the ".3e" this template
#    used before Phase 1 it was ~5e-4, which was enough to move the last digit
#    of the smallest totals - the reason the display precision changed.)
#
# 2e-3 clears both with two orders of magnitude to spare. It stays far below
# the failure modes this item can exhibit:
# reading the stated diameter as a radius is a factor of 16 in J, dropping the
# pi/2 is a factor of 1.57, taking the segments in parallel instead of series
# or using only one of them is tens of percent.
TOLERANCE = 2e-3

SOURCE = (
    "Angle of twist of a circular shaft in torsion (Beer & Johnston, Mechanics "
    "of Materials ch. 3): phi = T*L/(J*G) with J = (pi/2)*c^4 for a solid "
    "section; segments in series carry the same torque, so the twists add. "
    "Worked in base SI units (N.m, m, Pa)."
)

_NUM = r"([-+]?\d[\d,]*(?:\.\d+)?)"
_T_RE = re.compile(r"a torque of\s+" + _NUM + r"\s*N\.m\s+is applied",
                   re.IGNORECASE)
# "Segment AB is a solid structural steel shaft of length 1.02 m and diameter
#  85 mm." - anchored on the segment name so AB and BC cannot be confused.
_SEG_RE = (r"Segment {name} is a solid .*?shaft of length\s+" + _NUM +
           r"\s*m\s+and diameter\s+" + _NUM + r"\s*mm")
_AB_RE = re.compile(_SEG_RE.format(name="AB"), re.IGNORECASE)
_BC_RE = re.compile(_SEG_RE.format(name="BC"), re.IGNORECASE)
# "(Use G = 79.3 GPa for structural steel and G = 26.0 GPa for aluminum alloy)"
# - the two moduli appear in segment order, so position identifies them. Both
# segments can be the same material, in which case both values are equal and
# the pairing is unambiguous anyway.
_G_RE = re.compile(r"\bG\s*=\s*" + _NUM + r"\s*GPa", re.IGNORECASE)

_ANS_RE = re.compile(
    r"total angle of twist at the free end C is\s*" + _NUM + r"\s*radians",
    re.IGNORECASE)


def _f(tok: str) -> float:
    return float(tok.replace(",", ""))


def parse_givens(question: str) -> dict:
    t_m, ab_m, bc_m = (_T_RE.search(question), _AB_RE.search(question),
                       _BC_RE.search(question))
    g_all = _G_RE.findall(question)
    if not (t_m and ab_m and bc_m) or len(g_all) < 2:
        raise ValueError(
            "composite_shafts_series: question did not match expected form")
    return {
        "T_Nm": _f(t_m.group(1)),
        "L1_m": _f(ab_m.group(1)), "d1_mm": _f(ab_m.group(2)),
        "L2_m": _f(bc_m.group(1)), "d2_mm": _f(bc_m.group(2)),
        "G1_GPa": _f(g_all[0]), "G2_GPa": _f(g_all[1]),
    }


def _twist(T: float, L: float, d_mm: float, G_GPa: float) -> float:
    c = d_mm / 2.0 / 1000.0                  # mm diameter -> m radius
    J = (math.pi / 2.0) * c ** 4             # m^4
    return T * L / (J * G_GPa * 1e9)         # rad


def recompute(givens: dict) -> float:
    """Total angle of twist at C, in radians, from the stated values only."""
    phi = (_twist(givens["T_Nm"], givens["L1_m"], givens["d1_mm"],
                  givens["G1_GPa"])
           + _twist(givens["T_Nm"], givens["L2_m"], givens["d2_mm"],
                    givens["G2_GPa"]))
    # Quantise to the trace's display convention - five significant figures in
    # radians, half-up in decimal, floored at 5 dp. A presentation convention
    # read off the generated text, not a borrowed formula. Without it the
    # tolerance would have to absorb a full display step, which at the small
    # end of the range (0.006 rad) is 8e-4 relative and no longer resolves a
    # single-digit disagreement.
    if phi <= 0:
        raise ValueError("composite_shafts_series: non-positive twist")
    places = max(5, 4 - math.floor(math.log10(phi)))
    return float(Decimal(repr(phi)).quantize(Decimal(1).scaleb(-places),
                                             rounding=ROUND_HALF_UP))


def gold_answer(solution: str) -> float:
    i = solution.find("**Answer:**")
    block = solution[i:] if i >= 0 else solution
    m = _ANS_RE.search(block)
    if not m:
        raise ValueError("composite_shafts_series: could not read the answer block")
    return _f(m.group(1))
