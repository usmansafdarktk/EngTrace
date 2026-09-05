"""Independent round-trip oracle for `template_cantilever_double_integration`.

Derived from the classical cantilever solutions of the double-integration
method (Hibbeler, *Structural Analysis*, 10th ed. SI, ch. 8 - the grounding the
template's own docstring cites), NOT from the template's code:

    concentrated load P at the free end
        M(x) = -P*(L - x)
        E*I*v'' = M  ->  E*I*v' = -P*(L*x - x^2/2)      [v'(0)=0 => C1=0]
                     ->  E*I*v  = -P*(L*x^2/2 - x^3/6)  [v(0)=0  => C2=0]
        at x = L:  |delta_tip| = P*L^3 / (3*E*I)

    uniform load w over the whole span
        M(x) = -w*(L - x)^2 / 2
        integrating twice with v(0)=v'(0)=0
        at x = L:  |delta_tip| = w*L^4 / (8*E*I)

Units.  The oracle works in strict SI base units, which is deliberately a
*different* unit path from the trace (the trace works in kN and kN/m^2):

    I : the question states "I = <a> x 10^6 mm^4"  ->  a * 1e6 * 1e-12 m^4
    E : "E = <e> GPa"                              ->  e * 1e9 Pa
    P : "P = <p> kN"                               ->  p * 1e3 N
    w : "w = <q> kN/m"                             ->  q * 1e3 N/m
    L : "L = <l> m"                                ->  l m

delta comes out in metres and is multiplied by 1000 for the mm answer.  Doing
the conversion at the very end, from the full-precision metre value, is the
independent version of the step the audit flags: the trace prints a Step-5
value already rounded to 5 dp in metres and then prints "<that> * 1000 = <mm>",
a line whose own arithmetic does not always reproduce (measured: 11/300 seeds,
e.g. seed 48 "delta = 0.01685 * 1000 = 16.8 mm").  That is a display defect in
the printed line - the *answer* is rounded from the unrounded metre value, and
this oracle confirms the answer itself round-trips.
"""
from __future__ import annotations

import re

TEMPLATE_ID = "template_cantilever_double_integration"

# Relative tolerance.
#
# The answer is printed to a single decimal place in millimetres over a
# deflection window the docstring caps at [3, 45] mm, so half a unit in the last
# printed place is up to 0.05/3 = 1.7% of the answer.  A plain full-precision
# comparison would therefore need TOLERANCE >= 0.017 to avoid firing on nothing
# but display rounding (measured: worst full-precision disagreement over 300
# seeds is 0.66%, at seed 198 where 6.0395 mm prints as "6.0"), and a 1.7%
# tolerance is far too slack to be worth running.
#
# So `recompute` quantises its own result to the same one-decimal convention the
# trace prints (see below) and this tolerance is pure headroom for a value that
# lands on a rounding tie: on a correct instance the two sides then agree
# exactly (measured worst relative error over 300 seeds: 0.0).  0.001 is still
# tighter than one display step anywhere in the window - one step of 0.1 mm at
# the largest permitted deflection, 45 mm, is 0.22% - so any genuine
# single-step disagreement is caught, and a wrong load case (3 vs 8 in the
# denominator), a dropped mm conversion or a wrong I all show up as tens of
# percent or more.
TOLERANCE = 0.001

SOURCE = (
    "Double-integration cantilever tip deflections (Hibbeler, Structural "
    "Analysis 10th ed. SI, ch.8): point load at the free end "
    "delta = P*L^3/(3*E*I); full-span UDL delta = w*L^4/(8*E*I); "
    "worked in SI base units (N, m, Pa) and converted to mm at the end."
)

_I_RE = re.compile(
    r"moment of inertia\s+I\s*=\s*([\d.]+)\s*x\s*10\^6\s*mm\^4", re.IGNORECASE
)
_L_RE = re.compile(r"extends\s+L\s*=\s*([\d.]+)\s*m\b", re.IGNORECASE)
_E_RE = re.compile(r"\bE\s*=\s*([\d.]+)\s*GPa\b", re.IGNORECASE)
# The UDL pattern must be tried first: "w = 11.9 kN/m" would also satisfy a
# careless "kN" match.
_W_RE = re.compile(r"\bw\s*=\s*([\d.]+)\s*kN/m\b", re.IGNORECASE)
_P_RE = re.compile(r"\bP\s*=\s*([\d.]+)\s*kN\b(?!\s*/)", re.IGNORECASE)


def parse_givens(question: str) -> dict:
    i_m = _I_RE.search(question)
    l_m = _L_RE.search(question)
    e_m = _E_RE.search(question)
    if not (i_m and l_m and e_m):
        raise ValueError("cantilever: question did not match expected form")

    w_m = _W_RE.search(question)
    p_m = _P_RE.search(question)
    if w_m is not None:
        case, load = "uniform", float(w_m.group(1))      # kN/m
    elif p_m is not None:
        case, load = "point", float(p_m.group(1))        # kN
    else:
        raise ValueError("cantilever: no load case found in question")

    return {
        "I_mm4_e6": float(i_m.group(1)),   # as stated: multiples of 1e6 mm^4
        "L_m": float(l_m.group(1)),
        "E_GPa": float(e_m.group(1)),
        "case": case,
        "load": load,
    }


def recompute(givens: dict) -> float:
    L = givens["L_m"]                                    # m
    E = givens["E_GPa"] * 1e9                            # GPa -> Pa (N/m^2)
    I = givens["I_mm4_e6"] * 1e6 * 1e-12                 # 1e6 mm^4 -> m^4
    load = givens["load"] * 1e3                          # kN or kN/m -> N, N/m

    if givens["case"] == "point":
        delta_m = load * L ** 3 / (3.0 * E * I)
    else:
        delta_m = load * L ** 4 / (8.0 * E * I)

    delta_mm = delta_m * 1000.0
    # Quantise to the trace's display convention (1 dp in mm, read off the
    # generated text - a presentation convention, not a borrowed formula).
    # See the TOLERANCE comment: without this, the check would have to run at a
    # ~2% tolerance and would no longer be worth running.
    return float(f"{delta_mm:.1f}")


_ANS_RE = re.compile(
    r"\*\*Answer:\*\*.*?deflection at the free end is\s*([-+]?[\d,.]+)\s*mm",
    re.IGNORECASE | re.DOTALL,
)


def gold_answer(solution: str) -> float:
    m = _ANS_RE.search(solution)
    if not m:
        raise ValueError("cantilever: could not read the answer block")
    return float(m.group(1).replace(",", ""))
