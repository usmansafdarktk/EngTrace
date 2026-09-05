"""Round-trip oracle for `template_vibration_transmissibility`.

Derived from the classical base-excitation (support motion) solution for a
single-degree-of-freedom system, Rao, *Mechanical Vibrations*, sec. 3.6
("Response of a Damped System Under the Harmonic Motion of the Base"), also
Thomson ch. 3.5.  For a base moving as ``y = Y sin(omega t)`` the steady-state
absolute displacement amplitude of the mass obeys the displacement
transmissibility

    TR = X / Y = sqrt( (1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2) )

with ``omega_n = sqrt(k/m)`` and ``r = omega / omega_n``.  The numerator carries
the ``1 +`` term because the dashpot transmits force through the base as well as
the spring - that is what distinguishes this from the force-magnification factor.

Only the question text and the `**Answer:**` block are read; no expression is
taken from the generator.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = 'template_vibration_transmissibility'

# Relative tolerance.
#
# Two printed answers, both rounded to 3 decimal places: TR (spanning ~0.056 to
# ~4.9 over a 300-seed sweep) and the absolute amplitude in mm (~0.02 to ~22).
# Half a unit in the last displayed digit is therefore under 0.1% for typical
# instances and only climbs past 0.5% at the very bottom of each range.  0.5%
# clears that floor for the bulk of the sweep while remaining well under the
# error produced by carrying the 3-decimal frequency ratio into the
# transmissibility formula, which is strongly amplified near resonance where
# (1 - r^2)^2 is small.
#
# Known limits of a purely relative threshold here: instances with a very small
# TR (~0.06) trip 0.5% on display rounding alone, while the propagated-rounding
# error - which shows up as ~39% of instances printing an answer that is not the
# correct rounding of the derived value - is usually below 0.5% in relative
# terms and so passes.  Both facts are reported rather than tuned away.
TOLERANCE = 0.005

SOURCE = ('Base-excitation displacement transmissibility '
          'TR = sqrt((1+(2*zeta*r)^2)/((1-r^2)^2+(2*zeta*r)^2)), X = TR*Y, '
          'with omega_n = sqrt(k/m) and r = 2*pi*f/omega_n - Rao, '
          'Mechanical Vibrations, sec. 3.6.')

_NUM = r'([-+]?\d[\d,]*(?:\.\d+)?)'

_Q_RE = re.compile(
    r'mass\s*(?:of\s*)?' + _NUM + r'\s*kg'
    r'.*?stiffness of\s*' + _NUM + r'\s*N/m'
    r'.*?damping ratio of\s*' + _NUM +
    r'.*?frequency of\s*' + _NUM + r'\s*Hz'
    r'.*?amplitude of\s*' + _NUM + r'\s*mm',
    re.S | re.I)

_TR_RE = re.compile(r'transmissibility ratio is\s*\*\*\s*' + _NUM, re.I)
_X_RE = re.compile(r'amplitude of the instrument\S*\s*vibration is\s*\*\*\s*'
                   + _NUM + r'\s*mm', re.I)


def _f(tok: str) -> float:
    return float(tok.replace(',', '').rstrip('.'))


def parse_givens(question: str) -> dict | None:
    m = _Q_RE.search(question)
    if not m:
        return None
    mass, k, zeta, f_hz, y_mm = (_f(g) for g in m.groups())
    return {'m_kg': mass, 'k_N_per_m': k, 'zeta': zeta,
            'f_hz': f_hz, 'Y_mm': y_mm}


def recompute(givens: dict) -> dict | None:
    """{'TR': ratio, 'X_mm': absolute amplitude} from the stated values only."""
    if not givens:
        return None
    m = givens['m_kg']
    k = givens['k_N_per_m']
    zeta = givens['zeta']
    if m <= 0 or k <= 0:
        return None

    omega_n = math.sqrt(k / m)
    omega = 2.0 * math.pi * givens['f_hz']
    r = omega / omega_n

    cross = 2.0 * zeta * r                       # the 2*zeta*r coupling term
    numerator = 1.0 + cross ** 2                 # spring + dashpot transmission
    denominator = (1.0 - r ** 2) ** 2 + cross ** 2
    if denominator <= 0:
        return None
    tr = math.sqrt(numerator / denominator)
    return {'TR': tr, 'X_mm': tr * givens['Y_mm']}


def gold_answer(solution: str) -> dict | None:
    i = solution.find('**Answer:**')
    block = solution[i:] if i >= 0 else solution
    tr = _TR_RE.search(block)
    x = _X_RE.search(block)
    if not tr or not x:
        return None
    return {'TR': _f(tr.group(1)), 'X_mm': _f(x.group(1))}
