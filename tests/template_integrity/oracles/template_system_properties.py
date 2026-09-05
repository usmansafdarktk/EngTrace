"""Independent round-trip oracle for `template_system_properties`.

Derived from the standard SDOF spring-mass-damper relations as stated in the
template docstring's "Core Equations" and in Rao, *Mechanical Vibrations*,
ch. 2 - NOT from the template's code.

Forward chain, using only what the question states:

    omega_n = sqrt(k / m)                [rad/s]
    c_cr    = 2 * sqrt(k * m)            [N.s/m]   (= 2*m*omega_n)
    zeta    = c / c_cr
    class   = underdamped (zeta < 1) | critically damped (zeta == 1)
              | overdamped (zeta > 1)
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = "template_system_properties"

# Relative tolerance.  Measured worst relative error over 1000 seeds:
#   omega_n 1.3e-5, c_cr 9e-8, zeta 2.5e-5 - all consistent with 4 dp display
#   rounding alone (zeta is O(0.1-2), so a half-unit in the 4th decimal is worth
#   up to ~2.5e-4 relative on its own).
# 0.002 sits an order of magnitude above the measured worst case and two orders
# below anything a real formula error would produce.
#
# Note on the audit's "zeta from the stated rounded c differs by up to 4.5e-6":
# reproduced exactly (max |zeta(stated c) - sampled target zeta| = 4.500e-06 over
# 300 seeds, 1.43e-05 over 1000).  It is invisible in the round trip because the
# template itself prints zeta computed from the stated c, so the printed chain
# closes; the residue only exists against the generator's hidden sampled zeta.
TOLERANCE = 0.002

SOURCE = (
    "SDOF spring-mass-damper properties (Rao, Mechanical Vibrations, ch.2): "
    "omega_n=sqrt(k/m), c_cr=2*sqrt(k*m), zeta=c/c_cr, with the standard "
    "zeta<1 / =1 / >1 damping classification."
)

_M_RE = re.compile(r"Mass \(m\)\s*=\s*([\d.,]+)\s*kg")
_K_RE = re.compile(r"(?:Spring )?Stiffness \(k\)\s*=\s*([\d.,]+)\s*N/m")
_C_RE = re.compile(r"Damping Coefficient \(c\)\s*=\s*([\d.,]+)\s*N")


def _num(tok: str) -> float:
    return float(tok.replace(",", ""))


def parse_givens(question: str) -> dict:
    m = _M_RE.search(question)
    k = _K_RE.search(question)
    c = _C_RE.search(question)
    if not (m and k and c):
        raise ValueError("system_properties: question did not match expected form")
    mv, kv, cv = _num(m.group(1)), _num(k.group(1)), _num(c.group(1))
    if mv <= 0 or kv <= 0 or cv < 0:
        raise ValueError("system_properties: non-physical m/k/c")
    return {"m": mv, "k": kv, "c": cv}


def _classify(zeta: float) -> str:
    if abs(zeta - 1.0) <= 1e-12:
        return "critically damped"
    return "underdamped" if zeta < 1.0 else "overdamped"


def recompute(givens: dict) -> dict:
    m, k, c = givens["m"], givens["k"], givens["c"]
    omega_n = math.sqrt(k / m)
    c_cr = 2.0 * math.sqrt(k * m)
    zeta = c / c_cr
    return {
        "omega_n": omega_n,
        "c_cr": c_cr,
        "zeta": zeta,
        "classification": _classify(zeta),
    }


_ON_RE = re.compile(r"Undamped Natural Frequency:\s*([\d.,]+)\s*rad/s")
_CC_RE = re.compile(r"Critical Damping Coefficient:\s*([\d.,]+)\s*N")
_ZE_RE = re.compile(r"Damping Ratio:\s*([\d.,]+)")
_CL_RE = re.compile(r"System Type:\s*([A-Za-z ]+)")


def gold_answer(solution: str) -> dict:
    tail = solution.split("**Answer:**")[-1]
    on, cc, ze, cl = (_ON_RE.search(tail), _CC_RE.search(tail),
                      _ZE_RE.search(tail), _CL_RE.search(tail))
    if not (on and cc and ze and cl):
        raise ValueError("system_properties: could not read the answer block")
    return {
        "omega_n": _num(on.group(1)),
        "c_cr": _num(cc.group(1)),
        "zeta": _num(ze.group(1)),
        "classification": cl.group(1).strip().lower(),
    }
