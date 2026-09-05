"""Independent round-trip oracle for `template_statically_indeterminate_shaft`.

Derived from the two governing equations for a shaft fixed at both ends with a
torque applied at an interior point (Beer & Johnston, *Mechanics of Materials*,
ch. 3 - statically indeterminate shafts; the pair the template's docstring
cites), **not** from the template's code:

    statics        T_A + T_B = T_applied
    compatibility  phi_AC = phi_CB, i.e. T_A*L_AC/(J*G) = T_B*L_BC/(J*G)

J and G are common to both segments and cancel, leaving T_A*L_AC = T_B*L_BC.
Solving the pair directly:

    T_A = T_applied * L_BC / L        T_B = T_applied * L_AC / L

where L = L_AC + L_BC is the total length.

Independence
------------
This is the closed-form solution of the two-equation system, reached by
elimination. The trace instead walks the substitution route - it forms the
length ratio L_AC/L_BC, states it to 3 dp, adds 1, and divides the applied
torque by that - which is where its defect lived: the divisor was printed to
3 dp while T_A was computed from the unrounded ratio, leaving 47% of probed
lines non-closing (seed 0 printed "3900 / 1.291 = 3021.85" where the printed
operands give 3020.91).

Because the trace legitimately consumes a 3-dp intermediate that it states
before use (P3), and this oracle does not form that intermediate at all, the
two agree only to within that stated rounding - which is exactly the point.
The tolerance below is sized to it, and is still two orders of magnitude
tighter than the defect the check exists to catch.

Note the diameter and shear modulus are stated in the question and are
deliberately NOT used: J and G cancel out of the compatibility equation, so an
answer that depends on them is wrong. They are distractors, and this oracle
ignoring them is a property worth preserving.
"""
from __future__ import annotations

import re
from decimal import Decimal, ROUND_HALF_UP

TEMPLATE_ID = "template_statically_indeterminate_shaft"

# Relative tolerance.
#
# Two effects have to fit inside it, and only one of them is display rounding.
#
# 1. Display. The reactions are quoted to 2 dp on values of roughly 300-3,000
#    N.m, so half a unit in the last printed place is at most 1.7e-5 relative.
#
# 2. The stated 3-dp length ratio. The trace states `L_AC/L_BC` to 3 dp and
#    divides by `1 + that`, as a solver reading the trace would. This oracle
#    uses the exact lengths instead. The two differ by the relative error of a
#    3-dp rounding of a ratio in roughly [0.25, 4], i.e. up to 5e-4/1.25 = 4e-4
#    on the divisor and therefore on T_A. T_B = T - T_A absorbs the same
#    absolute error over a smaller value, so its relative error is larger -
#    measured worst over 1,000 seeds is ~1.3e-3.
#
# Measured worst over 1,000 seeds: 1.35e-3. 5e-3 clears both effects with
# headroom while staying far below the failure modes that matter here, which
# are structural rather than numerical - swapping T_A and T_B, using L_AC where
# L_BC belongs, or letting J or G into the answer - all of which are tens of
# percent.
TOLERANCE = 5e-3

# Measured, not argued (Phase 1 review A, F-4), over 2,000 seeds with
#     python -m tests.template_integrity.oracle_floors 2000 template_statically_indeterminate_shaft
#
# AGREEMENT_FLOOR - worst relative disagreement between this oracle and
#   the trace on a CORRECT instance. TOLERANCE must exceed it, or the
#   check fires on presentation alone.
# DETECTION_FLOOR - smallest uniform relative error injected into the
#   gold answer that this oracle catches on >=99% of instances. This is
#   the number a downstream verifier can rely on, and it is set by the
#   display quantisation, NOT by TOLERANCE: both sides are quantised to
#   the printed precision before comparison, so an error smaller than
#   about half a display step is invisible whatever TOLERANCE says.
AGREEMENT_FLOOR = 1.40e-03
DETECTION_FLOOR = 1.00e-2

SOURCE = (
    "Statically indeterminate shaft fixed at both ends (Beer & Johnston, "
    "Mechanics of Materials ch. 3): statics T_A + T_B = T, compatibility "
    "T_A*L_AC = T_B*L_BC, solved by elimination to T_A = T*L_BC/L, "
    "T_B = T*L_AC/L."
)

_NUM = r"([-+]?\d[\d,]*(?:\.\d+)?)"
_L_RE = re.compile(r"length\s+" + _NUM + r"\s*m\b", re.IGNORECASE)
_T_RE = re.compile(r"torque of\s+" + _NUM + r"\s*N\.m\b", re.IGNORECASE)
_LAC_RE = re.compile(r"located\s+" + _NUM + r"\s*m\s+from end A", re.IGNORECASE)

_A_RE = re.compile(r"T_A\s*=\s*" + _NUM + r"\s*N\.m")
_B_RE = re.compile(r"T_B\s*=\s*" + _NUM + r"\s*N\.m")


def _f(tok: str) -> float:
    return float(tok.replace(",", ""))


def parse_givens(question: str) -> dict:
    l_m, t_m, ac_m = (_L_RE.search(question), _T_RE.search(question),
                      _LAC_RE.search(question))
    if not (l_m and t_m and ac_m):
        raise ValueError(
            "statically_indeterminate_shaft: question did not match expected form")
    L = _f(l_m.group(1))
    L_AC = _f(ac_m.group(1))
    if not 0.0 < L_AC < L:
        raise ValueError(
            f"statically_indeterminate_shaft: L_AC={L_AC} outside (0, {L})")
    return {"L_m": L, "T_applied": _f(t_m.group(1)), "L_AC_m": L_AC}


def recompute(givens: dict) -> dict:
    L = givens["L_m"]
    T = givens["T_applied"]
    L_AC = givens["L_AC_m"]
    L_BC = L - L_AC

    T_A = T * L_BC / L
    T_B = T * L_AC / L

    q = Decimal("0.01")            # the trace quotes both reactions to 2 dp
    return {
        "T_A": float(Decimal(repr(T_A)).quantize(q, rounding=ROUND_HALF_UP)),
        "T_B": float(Decimal(repr(T_B)).quantize(q, rounding=ROUND_HALF_UP)),
    }


def gold_answer(solution: str) -> dict:
    i = solution.find("**Answer:**")
    block = solution[i:] if i >= 0 else solution
    a, b = _A_RE.search(block), _B_RE.search(block)
    if not (a and b):
        raise ValueError(
            "statically_indeterminate_shaft: could not read the answer block")
    return {"T_A": _f(a.group(1)), "T_B": _f(b.group(1))}
