"""Independent round-trip oracle for `template_logarithmic_decrement`.

Derived from the free-decay relations for an underdamped SDOF oscillator as
stated in the template docstring's "Core Equations" and in Rao, *Mechanical
Vibrations*, sec. 2.6 (logarithmic decrement) - NOT from the template's code.

Forward chain, using only what the question states:

    delta = (1/n) * ln(x1 / x_{n+1})
    zeta  = delta / sqrt((2*pi)^2 + delta^2)

Note on this template's "back-solving".  The generator picks a true damping
ratio first and derives the two amplitudes from it, then STATES those amplitudes
rounded to 1 decimal place.  Recovering zeta from the rounded amplitudes
therefore does not return the hidden true zeta - and it should not.  That is the
physics of the exercise: the student is given measured amplitudes and asked what
damping they imply.  The property that must hold, and the only one this oracle
asserts, is that the printed delta and zeta are consistent with the amplitudes
AS STATED IN THE QUESTION.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = "template_logarithmic_decrement"

# Relative tolerance.  Measured worst relative error over 1000 seeds:
#   delta 1.6e-4, zeta 9.3e-4.  Both are display rounding, not arithmetic:
#   delta and zeta are printed to 4 decimals, and zeta - O(0.01-0.2) - therefore
#   carries up to ~5e-5/zeta of relative slop; the trace also displays zeta
#   computed from the already-rounded delta, worth another ~1e-4.
# 0.005 covers the measured worst case with ~5x headroom.
#
# The tolerance is deliberately NOT inflated to absorb the generator's 1-dp
# amplitude rounding.  That rounding is applied to the *inputs*, which this
# oracle reads exactly as stated, so it cancels out of the round trip.  For the
# record: an oracle that instead tried to recover the generator's hidden true
# zeta would need TOLERANCE ~= 0.12 (measured max 11.0%, p99 8.0%, median 0.35%
# over 1000 seeds) - a tolerance so wide it would no longer test anything.  That
# spread is the intended measurement uncertainty of the exercise, and asserting
# against the stated amplitudes is the only version of this check with teeth.
TOLERANCE = 0.005

SOURCE = (
    "Logarithmic decrement for free decay of an underdamped SDOF system "
    "(Rao, Mechanical Vibrations, sec. 2.6): delta=(1/n)*ln(x1/x_{n+1}), "
    "zeta=delta/sqrt((2*pi)^2+delta^2)."
)

_X1_RE = re.compile(r"initial amplitude is measured to be\s*([\d.]+)\s*mm", re.I)
_N_RE = re.compile(r"After\s*(\d+)\s*complete cycles", re.I)
_XN_RE = re.compile(
    r"complete cycles,\s*the amplitude is measured to be\s*([\d.]+)\s*mm", re.I
)


def parse_givens(question: str) -> dict:
    x1 = _X1_RE.search(question)
    n = _N_RE.search(question)
    xn = _XN_RE.search(question)
    if not (x1 and n and xn):
        raise ValueError("log_decrement: question did not match expected form")
    x1v, xnv, nv = float(x1.group(1)), float(xn.group(1)), int(n.group(1))
    if x1v <= 0 or xnv <= 0 or nv <= 0:
        raise ValueError("log_decrement: non-physical amplitudes/cycle count")
    return {"x1": x1v, "x_n_plus_1": xnv, "n": nv}


def recompute(givens: dict) -> dict:
    n = givens["n"]
    delta = (1.0 / n) * math.log(givens["x1"] / givens["x_n_plus_1"])
    zeta = delta / math.sqrt((2.0 * math.pi) ** 2 + delta ** 2)
    return {"delta": delta, "zeta": zeta}


_DELTA_RE = re.compile(r"The logarithmic decrement is\s*([-+]?[\d.]+)")
_ZETA_RE = re.compile(r"The damping ratio of the system is\s*([-+]?[\d.]+)")


def gold_answer(solution: str) -> dict:
    tail = solution.split("**Answer:**")[-1]
    d = _DELTA_RE.search(tail)
    z = _ZETA_RE.search(tail)
    if not (d and z):
        raise ValueError("log_decrement: could not read the answer block")
    return {"delta": float(d.group(1).rstrip(".")),
            "zeta": float(z.group(1).rstrip("."))}
