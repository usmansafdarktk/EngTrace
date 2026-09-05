"""Independent round-trip oracle for `template_poissons_ratio`.

Derived from the uniaxial Hooke's-law / Poisson relations as stated in the
template docstring's "Core Equations" and in Hibbeler, *Mechanics of
Materials*, ch. 3 (Poisson's ratio) - NOT from the template's own code.

Forward chain, using only what the question states:

    A            = pi * (d / 2)^2                 (solid circular section)
    sigma        = P / A                          (uniform axial stress)
    eps_axial    = sigma / E                      (Hooke's law, uniaxial)
    eps_lateral  = -nu * eps_axial                (Poisson)
    delta_d      = d * eps_lateral
    d_final      = d + delta_d

Sign convention: the question tags the load "(tension)" or "(compression)" in
prose while stating its magnitude unsigned.  Compression is taken as P < 0,
which is what makes the trace's own "sigma = P / A = (-7000 N) / ..." line and
the signed answer physically meaningful.  (The audit's "sign leak" is exactly
this: the question shows |P|, the trace shows -|P|.  It is a presentation
inconsistency, not an arithmetic one - see the report.)

Unit systems.  The template is dual-unit.  Both consistent sets collapse to the
same numeric formula because the unit factors cancel:

    SI : d [mm],  P [kN],   E [GPa]  ->  eps = P/(A*E) with A in mm^2
         (P*1e3 N) / (A mm^2 * E*1e3 MPa) = P / (A*E)
    US : d [in],  P [kips], E [ksi]  ->  eps = P/(A*E) with A in in^2

so `recompute` needs no unit branching beyond asserting the sets are consistent;
delta_d comes out in the same length unit the diameter was given in.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = "template_poissons_ratio"

# Relative tolerance.  Because `recompute` quantises to the trace's own 5 dp
# display convention (see below), the two sides agree exactly on a correct
# instance: measured worst relative error over 1000 seeds is 0.0, on both unit
# systems.  0.005 is therefore pure headroom - it exists only for a value that
# lands on a 0.5e-5 rounding tie, and it is nowhere near wide enough to absorb
# a real defect (a back-solved delta_d or a dropped sign shows up as tens of
# percent, or as a sign flip).
#
# Without the quantisation this tolerance would have to be ~0.02: part (a) is
# often only 2 significant figures once printed (e.g. 0.00026 in at seed 289,
# where a half-unit in the last printed place is 1.5% of the answer), and a 2%
# tolerance would be too slack to be worth running.
TOLERANCE = 0.005

SOURCE = (
    "Uniaxial Hooke's law and Poisson's ratio (Hibbeler, Mechanics of "
    "Materials, ch.3): sigma=P/A, eps_ax=sigma/E, eps_lat=-nu*eps_ax, "
    "delta_d=d*eps_lat, d_final=d+delta_d; A=pi*(d/2)^2."
)

_D_RE = re.compile(r"initial diameter of\s*([-+]?[\d.]+)\s*(mm|in)\b")
_P_RE = re.compile(
    r"axial load of\s*([-+]?[\d.]+)\s*(kN|kips)\s*\((tension|compression)\)",
    re.IGNORECASE,
)
_E_RE = re.compile(r"Modulus of Elasticity \(E\)\s*=\s*([-+]?[\d.]+)\s*(GPa|ksi)")
_NU_RE = re.compile(r"Poisson's Ratio \(nu\)\s*=\s*([-+]?[\d.]+)")

_SI = {"d": "mm", "P": "kN", "E": "GPa"}
_US = {"d": "in", "P": "kips", "E": "ksi"}


def parse_givens(question: str) -> dict:
    d_m = _D_RE.search(question)
    p_m = _P_RE.search(question)
    e_m = _E_RE.search(question)
    n_m = _NU_RE.search(question)
    if not (d_m and p_m and e_m and n_m):
        raise ValueError("poissons_ratio: question did not match expected form")

    units = {"d": d_m.group(2), "P": p_m.group(2), "E": e_m.group(2)}
    if units == _SI:
        system = "SI"
    elif units == _US:
        system = "US"
    else:
        raise ValueError(f"poissons_ratio: mixed/unknown unit set {units}")

    sense = p_m.group(3).lower()
    return {
        "d": float(d_m.group(1)),                 # mm or in
        "P": float(p_m.group(1)),                 # kN or kips, magnitude
        "sense": sense,                           # 'tension' | 'compression'
        "E": float(e_m.group(1)),                 # GPa or ksi
        "nu": float(n_m.group(1)),
        "system": system,
        "length_unit": units["d"],
    }


def recompute(givens: dict) -> dict:
    d = givens["d"]
    E = givens["E"]
    nu = givens["nu"]
    P = givens["P"] * (-1.0 if givens["sense"] == "compression" else 1.0)

    area = math.pi * (d / 2.0) ** 2               # mm^2 or in^2
    sigma = P / area                              # MPa (SI set) or ksi (US set)
    # Unit bookkeeping: the SI set gives sigma in kN/mm^2 = GPa against E in GPa;
    # the US set gives sigma in kips/in^2 = ksi against E in ksi.  Either way the
    # ratio is dimensionless, so no unit branching is needed here.
    eps_axial = sigma / E
    eps_lateral = -nu * eps_axial
    delta_d = d * eps_lateral
    d_final = d + delta_d
    # Quantise to the trace's display convention (5 dp for both parts) before
    # comparing.  This is a *presentation* convention read off the generated
    # text, not a borrowed formula: without it, part (a) - which is often only
    # 2 significant figures once printed - would need a ~5% relative tolerance
    # just to absorb its own last printed digit, which is far too loose to
    # detect the defects this check exists for.
    return {"a_delta_d": round(delta_d, 5), "b_d_final": round(d_final, 5)}


_A_RE = re.compile(
    r"a\)\s*The change in diameter is\s*\*\*\s*([-+]?[\d.,]+)\s*(?:mm|in)\s*\*\*"
)
_B_RE = re.compile(
    r"b\)\s*The final diameter is\s*\*\*\s*([-+]?[\d.,]+)\s*(?:mm|in)\s*\*\*"
)


def gold_answer(solution: str) -> dict:
    a = _A_RE.search(solution)
    b = _B_RE.search(solution)
    if not (a and b):
        raise ValueError("poissons_ratio: could not read the answer block")
    return {
        "a_delta_d": float(a.group(1).replace(",", "")),
        "b_d_final": float(b.group(1).replace(",", "")),
    }
