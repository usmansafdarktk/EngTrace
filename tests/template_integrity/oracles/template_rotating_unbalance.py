"""Round-trip oracle for `template_rotating_unbalance`.

Derived from the standard forced-response solution for a single-degree-of-freedom
system driven by a rotating unbalance (Rao, *Mechanical Vibrations*, ch. 3.6;
Thomson, *Theory of Vibration with Applications*, ch. 3.3).  The equation of
motion for a machine of *total* mass ``M`` on a spring ``k`` and dashpot ``c``,
carrying an unbalanced mass ``m_u`` at eccentricity ``e`` and running at ``omega``:

    M x'' + c x' + k x = m_u e omega^2 sin(omega t)

so the harmonic-excitation amplitude is the rotating-unbalance force
``F_0 = m_u e omega^2`` divided by the magnitude of the complex impedance
``sqrt((k - M omega^2)^2 + (c omega)^2)``.

Nothing here is transcribed from the generator; only the question text and the
`**Answer:**` block are read.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = 'template_rotating_unbalance'

# Relative tolerance.
#
# The printed answer is an amplitude in millimetres rounded to 3 decimal places.
# Across a 300-seed sweep the answer spans roughly 0.002 - 22 mm, so half a unit
# in the last displayed digit is <= 0.05% for the bulk of the range (answers
# above ~1 mm) and only exceeds 0.5% for answers below ~0.1 mm.  0.5% therefore
# sits above the display-rounding floor for the great majority of instances
# while still being far below the error injected by re-deriving omega from the
# *stated* (integer-rounded) shaft speed, which reaches several percent.  A
# looser tolerance would let that genuine round-trip break through unnoticed.
#
# Known residual: a handful of instances (answers below ~0.05 mm, printed as e.g.
# "0.002 mm") exceed 0.5% on display rounding alone.  Those are not false
# positives so much as a second, separate defect - a final answer carried to one
# significant figure is not gradeable - and they are reported as such.
TOLERANCE = 0.005

SOURCE = ('Rotating-unbalance forced response, F_0 = m_u*e*omega^2 and '
          'X = F_0 / sqrt((k - M*omega^2)^2 + (c*omega)^2), with '
          'omega = 2*pi*N/60 - Rao, Mechanical Vibrations, sec. 3.6.')

# Numbers in the question carry thousands separators on some quantities
# (stiffness, RPM) and not on others (damping coefficient), so every field has
# to tolerate commas.
_NUM = r'([-+]?\d[\d,]*(?:\.\d+)?)'

_Q_RE = re.compile(
    r'total mass of\s*' + _NUM + r'\s*kg'
    r'.*?stiffness of\s*' + _NUM + r'\s*N/m'
    r'.*?damping coefficient of\s*' + _NUM + r'\s*N[.\-]?s/m'
    r'.*?unbalance equivalent to a mass of\s*' + _NUM + r'\s*kg'
    r'.*?eccentricity of\s*' + _NUM + r'\s*mm'
    r'.*?speed of\s*' + _NUM + r'\s*RPM',
    re.S | re.I)

_A_RE = re.compile(
    r'amplitude of vibration is\s*\*\*\s*' + _NUM + r'\s*mm', re.I)


def _f(tok: str) -> float:
    return float(tok.replace(',', ''))


def parse_givens(question: str) -> dict | None:
    m = _Q_RE.search(question)
    if not m:
        return None
    mass, k, c, m_u, e_mm, rpm = (_f(g) for g in m.groups())
    return {'m_total_kg': mass, 'k_N_per_m': k, 'c_Ns_per_m': c,
            'm_u_kg': m_u, 'e_mm': e_mm, 'speed_rpm': rpm}


def recompute(givens: dict) -> float | None:
    """Steady-state amplitude in millimetres, from the stated values only."""
    if not givens:
        return None
    m = givens['m_total_kg']
    k = givens['k_N_per_m']
    c = givens['c_Ns_per_m']
    m_u = givens['m_u_kg']
    e = givens['e_mm'] / 1000.0                      # mm -> m
    omega = givens['speed_rpm'] * 2.0 * math.pi / 60.0   # rev/min -> rad/s

    f0 = m_u * e * omega ** 2                        # rotating-unbalance force
    impedance = math.hypot(k - m * omega ** 2, c * omega)
    if impedance == 0:
        return None
    return (f0 / impedance) * 1000.0                 # m -> mm


def gold_answer(solution: str) -> float | None:
    i = solution.find('**Answer:**')
    block = solution[i:] if i >= 0 else solution
    m = _A_RE.search(block)
    return _f(m.group(1)) if m else None
