"""Derive every column of CONTROL_CHART_FACTORS from its definition (C3.7).

    python -m tests.constants_integrity.control_chart_factors            # compare with the table
    python -m tests.constants_integrity.control_chart_factors --selftest

The table was transcribed from Montgomery ISQC 7e Appendix Table VI, a copyrighted
book present on one machine. The brief asks for the columns to be DERIVED instead,
which needs no book. What each derivation rests on:

  ON DISK (docs/references/industrial/, public NIST/SEMATECH e-Handbook pages)
    pmc32.htm    c4 = sqrt(2/(n-1)) * (n/2 - 1)! / ((n-1)/2 - 1)!   (the page's own
                 formula; half-integer factorials through the Gamma function), and
                 its check value "the c4 factor for n = 10 is 0.9727"
    pmc321.htm   the X-bar chart from s: limits xbar +/- 3 s/(c4 sqrt n)  -> A3
                 the s chart: s +/- 3 (s/c4) sqrt(1 - c4^2)              -> B3, B4
                 the X-bar chart from R: xbar +/- 3 R/(d2 sqrt n)        -> A2
                 d2 defined as the mean of R/sigma
  NOT ON DISK, standard statistics stated here rather than cited
    d2(n) = E[R/sigma] = integral (1 - Phi(x)^n - (1 - Phi(x))^n) dx
    d3(n) = sd[R/sigma], from E[(R/sigma)^2] = 2 * double integral over x < y of
            (1 - Phi(y)^n - (1 - Phi(x))^n + (Phi(y) - Phi(x))^n) dx dy
    A = 3/sqrt(n); B5 = c4 - 3 sqrt(1 - c4^2), B6 = c4 + 3 sqrt(1 - c4^2);
    D1 = d2 - 3 d3, D2 = d2 + 3 d3; D3 = 1 - 3 d3/d2, D4 = 1 + 3 d3/d2;
    lower limits floored at 0. These are the known-sigma and from-R-bar limits
    with the 3-sigma convention pmc32 states.

Standard library only (statistics.NormalDist, math.lgamma): no new dependency. The
integrals use composite Simpson's rule on [-L, L]; the self-test checks the
quadrature against the closed forms it must reproduce (d2(2) = 2/sqrt(pi)) and
pmc32's own printed c4(10).
"""
from __future__ import annotations

import math
import os
import sys
from statistics import NormalDist

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

COLUMNS = ("A", "A2", "A3", "c4", "inv_c4", "B3", "B4", "B5",
           "B6", "d2", "inv_d2", "d3", "D1", "D2", "D3", "D4")
_N = NormalDist()
L, STEPS = 8.0, 1600          # Simpson panels on [-8, 8]; the tails beyond carry < 1e-14


def c4(n):
    """pmc32.htm's formula, with (k)! = Gamma(k + 1) for half-integer k."""
    return math.sqrt(2.0 / (n - 1)) * math.exp(math.lgamma(n / 2) - math.lgamma((n - 1) / 2))


def _simpson(f, a, b, m):
    h = (b - a) / m
    s = f(a) + f(b)
    for i in range(1, m):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3


def d2(n):
    return _simpson(lambda x: 1 - _N.cdf(x) ** n - (1 - _N.cdf(x)) ** n, -L, L, STEPS)


M = 400        # Simpson panels per axis for the range distribution (h = 0.04)


def _range_survival(n):
    """[(w, P(R/sigma > w))] on w in [0, 2L]: P(R <= w) = n * int phi(x) (Phi(x+w) - Phi(x))^(n-1) dx.

    Both integrals are over RECTANGULAR grids, so Simpson's weights are valid on each.
    The first version integrated the second moment over the triangle x < y with the
    full-interval Simpson weights applied to a partial inner interval - not a valid
    rule - and got d3 wrong by up to 13% (d3(20) 0.644 against 0.729); the closed
    form d3(2) = sqrt(2 - 4/pi) in the self-test is what exposed it.
    """
    hx = 2 * L / M
    xs = [-L + i * hx for i in range(M + 1)]
    wx = [(1 if i in (0, M) else (4 if i % 2 else 2)) * hx / 3 for i in range(M + 1)]
    pdf = [_N.pdf(x) for x in xs]
    hw = 2 * L / M
    out = []
    for j in range(M + 1):
        w = j * hw
        cdf_w = sum(wx[i] * pdf[i] * (_N.cdf(xs[i] + w) - _N.cdf(xs[i])) ** (n - 1) for i in range(M + 1))
        out.append((w, max(0.0, 1.0 - n * cdf_w)))
    return out, hw


def d3(n):
    """sd of R/sigma: E[R^2] = 2 int_0^inf w P(R > w) dw, E[R] = int_0^inf P(R > w) dw."""
    surv, hw = _range_survival(n)
    wts = [(1 if j in (0, M) else (4 if j % 2 else 2)) * hw / 3 for j in range(M + 1)]
    mean = sum(wt * s for wt, (_w, s) in zip(wts, surv))
    second = 2 * sum(wt * w * s for wt, (w, s) in zip(wts, surv))
    return math.sqrt(max(second - mean * mean, 0.0))


def factors(n):
    c, D2_, D3_ = c4(n), d2(n), d3(n)
    r = math.sqrt(n)
    return {
        "A": 3 / r, "A2": 3 / (D2_ * r), "A3": 3 / (c * r),
        "c4": c, "inv_c4": 1 / c,
        "B3": max(0.0, 1 - 3 * math.sqrt(1 - c * c) / c), "B4": 1 + 3 * math.sqrt(1 - c * c) / c,
        "B5": max(0.0, c - 3 * math.sqrt(1 - c * c)), "B6": c + 3 * math.sqrt(1 - c * c),
        "d2": D2_, "inv_d2": 1 / D2_, "d3": D3_,
        "D1": max(0.0, D2_ - 3 * D3_), "D2": D2_ + 3 * D3_,
        "D3": max(0.0, 1 - 3 * D3_ / D2_), "D4": 1 + 3 * D3_ / D2_,
    }


def _places(tok):
    return len(tok.split('.')[1]) if '.' in tok else 0


def compare(table=None, literal_tokens=None):
    """[(n, column, table value, derived value, rounded derived, agrees)]."""
    if table is None:
        from data.templates.branches.industrial_engineering.constants import CONTROL_CHART_FACTORS as table
    out = []
    for n, row in table.items():
        f = factors(n)
        for k, col in enumerate(COLUMNS):
            v = row[k]
            tok = (literal_tokens or {}).get((n, col), repr(v))
            dp = _places(tok)
            rd = round(f[col] + 0.0, dp)
            out.append((n, col, v, f[col], rd, abs(rd - v) < 0.5 * 10 ** -dp))
    return out


def run():
    rows = compare()
    bad = [r for r in rows if not r[5]]
    print(f'{len(rows)} cells derived; {len(rows) - len(bad)} agree with the table at its printed places, '
          f'{len(bad)} do not')
    for n, col, v, d, rd, _ok in bad:
        print(f'  n={n:2d} {col:6s} table {v!r:8}  derived {d:.6f} -> {rd!r}')
    return 1 if bad else 0


def selftest():
    bad = []
    checks = [
        ('d2(2) against its closed form 2/sqrt(pi)', abs(d2(2) - 2 / math.sqrt(math.pi)) < 1e-9),
        ('c4(10) against pmc32.htm\'s printed 0.9727', round(c4(10), 4) == 0.9727),
        ('d3(2) against its closed form sqrt(2 - 4/pi)', abs(d3(2) - math.sqrt(2 - 4 / math.pi)) < 1e-6),
    ]
    # a planted transcription error the comparison must catch: A2(5) 0.577 -> 0.557
    from data.templates.branches.industrial_engineering.constants import CONTROL_CHART_FACTORS as T
    planted = {n: list(r) for n, r in T.items()}
    planted[5][COLUMNS.index('A2')] = 0.557
    caught = [r for r in compare({5: tuple(planted[5])}) if not r[5] and r[1] == 'A2']
    checks.append(('a planted 0.577 -> 0.557 transposition in A2(5)', bool(caught)))
    planted[7][COLUMNS.index('d3')] = 0.883
    caught = [r for r in compare({7: tuple(planted[7])}) if not r[5] and r[1] == 'd3']
    checks.append(('a planted 0.833 -> 0.883 digit swap in d3(7), a double integral', bool(caught)))
    for label, ok in checks:
        print(f"  [{'ok' if ok else 'FAIL'}] {label}")
        if not ok:
            bad.append(label)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
