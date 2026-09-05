"""Independent round-trip oracle for `template_mean_variance`.

Definitions of the first two moments of a discrete random variable (the "Core
Equations" of the template docstring; Papoulis, *Probability, Random Variables
and Stochastic Processes*, ch. 5, or any communications text's review chapter),
written from the definitions rather than from the template's code:

    mu_X       = E[X]           = sum_i x_i * p_i
    sigma_X^2  = E[(X - mu_X)^2] = sum_i (x_i - mu_X)^2 * p_i

The variance is formed with the mean this oracle itself derived, because that is
what a solver reading only the question can do: the question states no mean.

THE VALUES USED ARE EXACTLY THE ONES PRINTED IN THE QUESTION.  That is the
property under test.  The question states each probability rounded to 3 decimal
places, so the printed set is all the information a solver has; if the trace's
answer was computed from unrounded probabilities it cannot be reproduced, and
this oracle is what shows it.  Nothing is renormalised - the printed
probabilities do not in general sum to 1 (measured sums over 300 seeds: 0.999,
1.000, 1.001, 1.002), and inventing a normalisation would paper over exactly
the loss of information the check exists to expose.
"""
from __future__ import annotations

import re

TEMPLATE_ID = "template_mean_variance"

# Relative tolerance.
#
# Both answers are printed to 3 decimal places, so a correct instance can differ
# from the exact value by at most half a unit in the last place, 0.0005
# absolute.  Over this template's own value ranges (|mu_X| is typically 5-15 and
# sigma_X^2 is typically 10-90) that is at most ~1e-4 relative, so TOLERANCE =
# 1e-3 is roughly ten times the legitimate display-rounding bound: nothing here
# fires on presentation.  It is still tight enough to be meaningful, since the
# effect under test - a mean computed from unrounded probabilities but stated
# from 3-dp ones - is a drift of order (0.0005 * sum|x_i|) on the mean, which is
# typically several times 1e-3 relative and, for the small means this template
# sometimes draws, tens of percent.
#
# Note the tolerance choice moves the reported failure rate a lot here, because
# the defect is a fine drift rather than a gross error: over the same 300 seeds
# the mean fails on 187 seeds at 5e-4 (the strict display bound), 120 at 1e-3,
# 40 at 2e-3 and 15 at 5e-3.  1e-3 is the conservative reading; the honest
# statement of the defect is the display-precision one - the printed answer is
# not what the printed probabilities give at the 3 dp it is quoted to, on
# 244/300 seeds for the mean and 261/300 for the variance.
TOLERANCE = 1e-3

SOURCE = (
    "Definitions of the first two moments of a discrete random variable: "
    "mu = sum(x_i*p_i), sigma^2 = sum((x_i-mu)^2*p_i), evaluated with the "
    "values and 3-dp probabilities exactly as printed in the question."
)

_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
_VALS_RE = re.compile(r"take the values\s*\{([^}]*)\}")
_PROBS_RE = re.compile(r"probabilities\s*\{([^}]*)\}")


def _numbers(blob: str) -> list[float]:
    return [float(t) for t in re.findall(_NUM, blob)]


def parse_givens(question: str) -> dict:
    v_m = _VALS_RE.search(question)
    p_m = _PROBS_RE.search(question)
    if not (v_m and p_m):
        raise ValueError("mean_variance: question did not match expected form")
    values = _numbers(v_m.group(1))
    probs = _numbers(p_m.group(1))
    if not values or len(values) != len(probs):
        raise ValueError(
            f"mean_variance: {len(values)} values but {len(probs)} probabilities")
    return {"values": values, "probs": probs}


def recompute(givens: dict) -> dict:
    values = givens["values"]
    probs = givens["probs"]                      # as printed: 3 dp, not renormalised
    mean = sum(v * p for v, p in zip(values, probs))
    variance = sum((v - mean) ** 2 * p for v, p in zip(values, probs))
    # Round only far below the printed precision (9 dp vs the trace's 3 dp), to
    # keep a mean that is analytically zero from showing up as 1e-16 and turning
    # a correct instance into a division-by-almost-zero relative error.
    return {"mean": round(mean, 9), "variance": round(variance, 9)}


_MEAN_RE = re.compile(
    r"mean of the random variable X is\s*\*\*\s*(" + _NUM + r")\s*\*\*")
_VAR_RE = re.compile(
    r"variance of the random variable X is\s*\*\*\s*(" + _NUM + r")\s*\*\*")


def gold_answer(solution: str) -> dict:
    tail = solution[solution.find("**Answer:**"):] if "**Answer:**" in solution else solution
    m = _MEAN_RE.search(tail)
    v = _VAR_RE.search(tail)
    if not (m and v):
        raise ValueError("mean_variance: could not read the answer block")
    return {"mean": float(m.group(1)), "variance": float(v.group(1))}
