"""Round-trip oracle for `template_damping_classification`.

Derived from the free-vibration solution of ``m x'' + c x' + k x = 0``
(Rao, *Mechanical Vibrations*, ch. 2.6-2.8; Thomson ch. 2.5).  The roots of the
characteristic equation change character at ``c = c_c``, which gives

    omega_n = sqrt(k/m)
    c_c     = 2 sqrt(k m)          (equivalently 2 m omega_n)
    zeta    = c / c_c
    omega_d = omega_n sqrt(1 - zeta^2)      (real only while zeta < 1)

and the classification underdamped / critically damped / overdamped according to
zeta < 1, zeta = 1, zeta > 1.

Critical damping is a measure-zero condition, so "zeta = 1" cannot be tested by
exact float equality against numbers that the question states to a finite number
of decimals.  This oracle therefore propagates the *stated* precision of the
damping coefficient: if c is printed to p decimals it is known only to within
half a unit in that digit, i.e. zeta is known only to within
``0.5*10**-p / c_c``.  A system whose zeta is within that band of 1 is
indistinguishable from critically damped on the evidence the question supplies,
and is labelled so.  Anything outside the band is classified strictly.  This is
the same rule a careful solver would apply, and it is the only rule under which
the question is answerable at all.

Only the question text and the `**Answer:**` block are read.
"""
from __future__ import annotations

import math
import re

TEMPLATE_ID = 'template_damping_classification'

# Relative tolerance for the numeric parts (zeta, omega_d).
#
# zeta is printed to 3 decimals and spans ~0.16 - 2.5 over a 300-seed sweep, so
# half a unit in the last digit is at most ~0.3%; omega_d is printed to 4
# decimals on values of order 10-100, i.e. a display floor far below that.
# 0.5% therefore clears display rounding while still catching a zeta or omega_d
# that was formed from rounded intermediates rather than from the stated m, k, c.
TOLERANCE = 0.005

SOURCE = ('Free-vibration damping parameters c_c = 2*sqrt(k*m), zeta = c/c_c, '
          'omega_n = sqrt(k/m), omega_d = omega_n*sqrt(1-zeta^2), with the '
          'zeta-vs-1 classification - Rao, Mechanical Vibrations, ch. 2.6-2.8.')

_NUM = r'([-+]?\d[\d,]*(?:\.\d+)?)'

_Q_RE = re.compile(
    r'mass of\s*' + _NUM + r'\s*kg'
    r'.*?stiffness of\s*' + _NUM + r'\s*N/m'
    r'.*?damping coefficient of\s*' + _NUM + r'\s*N[-.]?s/m',
    re.S | re.I)

_ZETA_RE = re.compile(r'damping ratio is\s*' + _NUM, re.I)
_CLASS_RE = re.compile(r'system is\s*\*\*\s*([A-Za-z ]+?)\s*\*\*', re.I)
_OMEGA_D_RE = re.compile(r'damped natural frequency is\s*' + _NUM + r'\s*rad/s', re.I)

UNDER = 'Underdamped'
CRITICAL = 'Critically Damped'
OVER = 'Overdamped'


def _f(tok: str) -> float:
    return float(tok.replace(',', '').rstrip('.'))


def _decimals(tok: str) -> int:
    body = tok.replace(',', '')
    return len(body.split('.')[1]) if '.' in body else 0


def parse_givens(question: str) -> dict | None:
    m = _Q_RE.search(question)
    if not m:
        return None
    mass_t, k_t, c_t = m.groups()
    return {'m_kg': _f(mass_t), 'k_N_per_m': _f(k_t), 'c_Ns_per_m': _f(c_t),
            # how precisely the question pins c down - see module docstring
            'c_decimals': _decimals(c_t)}


def recompute(givens: dict) -> dict | None:
    if not givens:
        return None
    m = givens['m_kg']
    k = givens['k_N_per_m']
    c = givens['c_Ns_per_m']
    if m <= 0 or k <= 0:
        return None

    omega_n = math.sqrt(k / m)
    c_c = 2.0 * math.sqrt(k * m)
    zeta = c / c_c

    # Half a unit in the last decimal the question shows for c, carried into zeta.
    band = max((0.5 * 10.0 ** -givens['c_decimals']) / c_c, 1e-12)

    if abs(zeta - 1.0) <= band:
        label = CRITICAL
    elif zeta < 1.0:
        label = UNDER
    else:
        label = OVER

    out: dict = {'zeta': zeta, 'classification': label}
    if label == UNDER:
        out['omega_d'] = omega_n * math.sqrt(1.0 - zeta ** 2)
    return out


def gold_answer(solution: str) -> dict | None:
    i = solution.find('**Answer:**')
    block = solution[i:] if i >= 0 else solution
    z = _ZETA_RE.search(block)
    cl = _CLASS_RE.search(block)
    if not z or not cl:
        return None
    out: dict = {'zeta': _f(z.group(1)), 'classification': cl.group(1).strip()}
    od = _OMEGA_D_RE.search(block)
    if od:
        out['omega_d'] = _f(od.group(1))
    return out
