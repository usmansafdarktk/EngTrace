"""Independent round-trip oracle for `template_annulus_flowrate`.

Steady laminar axial flow in the annulus between two concentric cylinders,
solved by a shell momentum balance - Bird, Stewart & Lightfoot, *Transport
Phenomena* 2nd ed., section 2.4, eq. 2.4-17, which is the "Core Equation" the
template's own docstring states.  Written from that result, not from the
template's code:

    kappa = R_inner / R_outer                                (dimensionless)

    Q = pi * (P0 - PL) * R_outer^4 / (8 * mu * L)
        * [ (1 - kappa^4) - (1 - kappa^2)^2 / ln(1/kappa) ]

The bracket is the annulus shape factor; it -> 1 as kappa -> 0, recovering
Hagen-Poiseuille Q = pi*dP*R^4/(8*mu*L), which is the sanity check used to
confirm the transcription above.

Units.  Everything is converted to SI base units before the formula is
evaluated (the question states the radii in cm, the length in m, the pressure
drop in kPa and the viscosity in Pa s; the parser accepts m/cm/mm for lengths
and Pa/kPa/MPa/bar for the pressure drop so that a later unit change in the
template does not silently disarm the check).  kappa is formed from the two
radii at full precision - deliberately NOT from the trace's Step-3 line, which
prints kappa rounded to 3 dp.

What this oracle found (see the report): the answer round-trips.  The audit's
claim that "the printed final multiplication does not close" is about the
Step-5 display line - "Q = <main term> * <shape factor> = <Q>", where the shape
factor is printed to 4 dp but multiplied at full precision, so 283/300 seeds
print a product that does not reproduce (worst 0.14%, seed 91).  The printed
*answer* is nevertheless the full-precision Q of the values stated in the
question, and the stated pressure drop (0.010-0.090 kPa, 3 dp) is the value
actually used - the main-term line reproduces from it exactly - so no hidden
back-solved precision leaks into the answer.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = "template_annulus_flowrate"

# Relative tolerance.
#
# The answer is printed in scientific notation with a 4-decimal mantissa
# (e.g. "9.5388e-06"), i.e. 5 significant figures, so half a unit in the last
# printed place is at most 0.00005/1.0 = 5e-5 relative (worst case, a mantissa
# just above 1.0) and as little as 5e-6 for a mantissa near 9.9.  TOLERANCE is
# set at 1e-4 = twice that worst-case display bound: tight enough that any real
# break in the chain is caught (using the printed 3-dp kappa instead of the
# exact ratio moves Q by ~0.1%, i.e. 10x this tolerance; a dropped cm->m
# conversion or a dropped kPa->Pa conversion moves it by orders of magnitude),
# and loose enough that pure last-digit display rounding never fires.
#
# Measured over 300 seeds: worst relative error 4.3e-5, which is exactly the
# display bound for its mantissa (1.0500e-04) - i.e. every instance agrees to
# within the last printed digit.
TOLERANCE = 1e-4

SOURCE = (
    "Shell momentum balance for annular laminar flow (Bird, Stewart & "
    "Lightfoot, Transport Phenomena 2nd ed., eq. 2.4-17): "
    "Q = pi*dP*R_o^4/(8*mu*L) * [(1-k^4) - (1-k^2)^2/ln(1/k)], k = R_i/R_o."
)

_LEN_TO_M = {"m": 1.0, "cm": 1e-2, "mm": 1e-3}
_PRESS_TO_PA = {"Pa": 1.0, "kPa": 1e3, "MPa": 1e6, "bar": 1e5}

_NUM = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

_RI_RE = re.compile(
    r"inner pipe has an outer radius of\s*" + _NUM + r"\s*(m|cm|mm)\b")
_RO_RE = re.compile(
    r"outer pipe has an inner radius of\s*" + _NUM + r"\s*(m|cm|mm)\b")
_L_RE = re.compile(
    r"pipes have a length of\s*" + _NUM + r"\s*(m|cm|mm)\b")
_DP_RE = re.compile(
    r"pressure drop of\s*" + _NUM + r"\s*(kPa|MPa|Pa|bar)\b")
# The viscosity unit is printed as "Pa" followed by a middle-dot/interpunct
# character that varies with the source encoding, so match only up to "Pa".
_MU_RE = re.compile(r"viscosity is\s*" + _NUM + r"\s*Pa")


def parse_givens(question: str) -> dict:
    ri = _RI_RE.search(question)
    ro = _RO_RE.search(question)
    ln = _L_RE.search(question)
    dp = _DP_RE.search(question)
    mu = _MU_RE.search(question)
    if not (ri and ro and ln and dp and mu):
        missing = [n for n, m in (("R_inner", ri), ("R_outer", ro), ("L", ln),
                                  ("dP", dp), ("mu", mu)) if not m]
        raise ValueError(f"annulus: question did not state {missing}")

    return {
        "R_inner_m": float(ri.group(1)) * _LEN_TO_M[ri.group(2)],
        "R_outer_m": float(ro.group(1)) * _LEN_TO_M[ro.group(2)],
        "L_m": float(ln.group(1)) * _LEN_TO_M[ln.group(2)],
        "dP_Pa": float(dp.group(1)) * _PRESS_TO_PA[dp.group(2)],
        "mu_Pa_s": float(mu.group(1)),
    }


def recompute(givens: dict) -> float:
    r_i = givens["R_inner_m"]
    r_o = givens["R_outer_m"]
    L = givens["L_m"]
    dP = givens["dP_Pa"]
    mu = givens["mu_Pa_s"]
    if not (0.0 < r_i < r_o) or L <= 0 or mu <= 0:
        raise ValueError(f"annulus: unphysical geometry {givens}")

    kappa = r_i / r_o                      # full precision, not the printed 3 dp
    main_term = math.pi * dP * r_o ** 4 / (8.0 * mu * L)
    shape = (1.0 - kappa ** 4) - (1.0 - kappa ** 2) ** 2 / math.log(1.0 / kappa)
    return main_term * shape               # m^3/s


_ANS_RE = re.compile(
    r"\*\*Answer:\*\*.*?volumetric flow rate is\s*\*\*\s*"
    + _NUM + r"\s*m\^3/s",
    re.IGNORECASE | re.DOTALL,
)


def gold_answer(solution: str) -> float:
    m = _ANS_RE.search(solution)
    if not m:
        raise ValueError("annulus: could not read the answer block")
    return float(m.group(1))
