"""Independent round-trip oracle for `template_shaft_design_power`.

Derived from the elementary torsion relations (Beer & Johnston, *Mechanics of
Materials*, ch. 3 - power transmission and shaft design; the same relations the
template's docstring cites), **not** from the template's code:

    power transmitted by a rotating shaft   P = 2*pi*f*T          (f in Hz)
    torsion formula, solid circular shaft   tau_max = T*c / J
    polar second moment, solid shaft        J = (pi/2)*c^4

Substituting J and solving for the radius at the allowable stress:

    tau = T*c / ((pi/2)*c^4) = 2*T / (pi*c^3)
    c   = ( 2*T / (pi*tau_allow) )^(1/3)
    d   = 2*c

Independence
------------
The oracle works in **base SI units throughout** - watts, pascals, hertz,
metres - and converts to millimetres only in the final line. The chain is a
single closed-form expression from the three stated quantities; nothing is
carried through an intermediate display, which is deliberately the opposite of
what the trace does. The defect this item carried was precisely an intermediate
display leak: Step 4 printed `d = 2 * <c to 4 dp> = <d from the UNROUNDED c>`,
so two seeds printed the identical operand `d = 2 * 0.0186` and different
answers, 0.0371 m and 0.0373 m.

The rotational speed is stated either in RPM or in Hz. Where it is stated in
RPM the trace converts it to Hz and states the result to 2 dp before using it,
so 2 dp is the value a solver has; this oracle applies the same 2-dp reading of
the stated speed, which is a **stated given**, not a borrowed intermediate. The
alternative - dividing by 60 at full precision - would test a quantity the
question does not state.
"""
from __future__ import annotations

import math
import re
from decimal import Decimal, ROUND_HALF_UP

TEMPLATE_ID = "template_shaft_design_power"

# Relative tolerance.
#
# `recompute` quantises to the 2 dp in millimetres the trace quotes the answer
# to, following reviews/phase0_oracle_findings.md Section 6: compare at the
# precision the answer is printed to, or the tolerance widens until it tests
# nothing. The diameter spans roughly 9-107 mm, so half a unit in the last
# printed place is at most 5.6e-4 relative at the small end and ~7e-5 at the
# median; without quantising, the tolerance would have to exceed 5.6e-4 and
# would no longer resolve a single display step.
#
# With both sides quantised the same way a correct instance agrees exactly, so
# 1e-3 is headroom for a value sitting on a display tie (the template resamples
# those, so in practice it is never used). It remains far tighter than the
# failure modes this item can exhibit: a factor of 2 from confusing radius with
# diameter, ~1.26x from a dropped cube root of 2, 60x from reading RPM as Hz.
TOLERANCE = 1e-3

SOURCE = (
    "Power transmission and solid-shaft design (Beer & Johnston, Mechanics of "
    "Materials ch. 3): P = 2*pi*f*T, tau = 2*T/(pi*c^3) for a solid circular "
    "section, hence c = (2*T/(pi*tau))^(1/3) and d = 2c. Worked in base SI "
    "units (W, Pa, Hz, m)."
)

_P_RE = re.compile(r"transmit\s+([\d,.]+)\s*kW\b", re.IGNORECASE)
_TAU_RE = re.compile(r"allowable shearing stress of\s+([\d,.]+)\s*MPa\b",
                     re.IGNORECASE)
# The speed is given either way round; try RPM first, since "3000 RPM" and
# "50 Hz" are distinguishable only by their unit.
_RPM_RE = re.compile(r"rotational speed of\s+([\d,.]+)\s*RPM\b", re.IGNORECASE)
_HZ_RE = re.compile(r"rotational speed of\s+([\d,.]+)\s*Hz\b", re.IGNORECASE)


def _f(tok: str) -> float:
    return float(tok.replace(",", ""))


def parse_givens(question: str) -> dict:
    p_m, tau_m = _P_RE.search(question), _TAU_RE.search(question)
    if not (p_m and tau_m):
        raise ValueError("shaft_design_power: question did not match expected form")

    rpm_m, hz_m = _RPM_RE.search(question), _HZ_RE.search(question)
    if rpm_m is not None:
        # Stated in RPM. The trace states the converted frequency to 2 dp and
        # that is what a solver reads, so read it the same way.
        f_hz = float(Decimal(repr(_f(rpm_m.group(1)) / 60.0)).quantize(
            Decimal("0.01"), rounding=ROUND_HALF_UP))
    elif hz_m is not None:
        f_hz = _f(hz_m.group(1))
    else:
        raise ValueError("shaft_design_power: no rotational speed found")

    return {"power_kW": _f(p_m.group(1)), "tau_MPa": _f(tau_m.group(1)),
            "f_Hz": f_hz}


def recompute(givens: dict) -> float:
    """Minimum shaft diameter in millimetres, from the stated values only."""
    P = givens["power_kW"] * 1e3                 # kW -> W
    tau = givens["tau_MPa"] * 1e6                # MPa -> Pa
    f = givens["f_Hz"]                           # Hz
    if f <= 0 or tau <= 0:
        raise ValueError("shaft_design_power: non-positive frequency or stress")

    torque = P / (2.0 * math.pi * f)             # N.m
    c = (2.0 * torque / (math.pi * tau)) ** (1.0 / 3.0)   # m
    d_mm = 2.0 * c * 1000.0

    # Quantise to the trace's display convention - 2 dp in mm, half-up in
    # decimal, which is what a reader doing decimal arithmetic applies. A
    # presentation convention read off the generated text, not a borrowed
    # formula.
    return float(Decimal(repr(d_mm)).quantize(Decimal("0.01"),
                                              rounding=ROUND_HALF_UP))


_ANS_RE = re.compile(
    r"minimum required diameter for the solid circular shaft is\s*"
    r"([-+]?[\d,.]+)\s*mm", re.IGNORECASE)


def gold_answer(solution: str) -> float:
    i = solution.find("**Answer:**")
    block = solution[i:] if i >= 0 else solution
    m = _ANS_RE.search(block)
    if not m:
        raise ValueError("shaft_design_power: could not read the answer block")
    return _f(m.group(1))
