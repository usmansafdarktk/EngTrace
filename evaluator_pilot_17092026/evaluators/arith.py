"""Check the arithmetic a text states, with sympy.

A derivation line like

    V_sat = 229.7 * 0.263**(0.396**0.2857) = 229.7 * 0.5715 = 131.27 cm³/mol

makes CLAIMS: each numeric segment equals the next. This module finds those
segments, evaluates each with sympy, and reports which consecutive pairs agree.
Symbolic segments (`Vc * Zc**(...)`) are skipped - they cannot be evaluated
without binding symbols - so a chain is checked from its first numeric segment on.

TOLERANCE. Traces round every displayed intermediate to 3-4 significant figures,
and a chain of three roundings can drift ~1%. ARITH_TOL is 1%; a claim that is off
by more than that did not come from the arithmetic shown.

VALIDATED ON GOLD FIRST. Gold arithmetic is correct by construction, so any claim
the checker marks inconsistent on a gold solution is a checker bug. The first cut
scored gold at 61.5% consistent - 111 bugs - and each rule below exists because of
one of them:
  - `8.24e-05 C`: the exponent was stripped as a unit, reading 8.24e-5 as 8.24
  - commas were deleted globally, fusing `E = 200 GPa, L = 4.8 m` into the claim
    200 = 4.8; clauses are now split on , ; and so then before splitting on =
  - `cos(-140.8)` in degrees was evaluated in radians; trig is tried both ways
  - `3.5% deteriorated case p2` was stripped to 3.5; a unit is now short tokens
  - `1/(41-23) hours = 60/18 min` is a correct unit change; equality under a unit
    factor is accepted when the two sides carry different units
That took gold to 94.9%. Five more, found the same way, took it to 100% (232/232):
  - `121.0 x 10^6` - a spaced `x` between numbers is multiplication
  - `1.21e8` read as the variable `e8`; exponents are removed before looking for words
  - `=> pi1 = 0.35/0.73`: an implication arrow is a clause break, not an equals sign
  - `1/(41-23) hours = 60/18 = 3.33 minutes`: a unit written once at the end of a
    chain belongs to every unit-less link before it
  - mm^4 -> m^4 is a factor of 1e-12, which the unit factors lacked

WHAT IT CANNOT PARSE IT SAYS SO. Every segment is classified: evaluated, symbolic
(skipped), or unparseable. A checker that silently parses 5% of lines and reports
"100% consistent" would be this repo's recurring defect, so coverage is reported
alongside every consistency figure.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import sympy
from sympy.parsing.sympy_parser import (convert_xor, implicit_multiplication,
                                        parse_expr, standard_transformations)

ARITH_TOL = 0.01
TRANSFORMS = standard_transformations + (implicit_multiplication, convert_xor)
FUNCS = {'sqrt': sympy.sqrt, 'log': sympy.log, 'ln': sympy.log, 'exp': sympy.exp,
         'sin': sympy.sin, 'cos': sympy.cos, 'tan': sympy.tan, 'asin': sympy.asin,
         'acos': sympy.acos, 'atan': sympy.atan, 'arctan': sympy.atan, 'pi': sympy.pi,
         'log10': lambda x: sympy.log(x, 10), 'abs': sympy.Abs}
FUNC_WORDS = set(FUNCS)
SUP = str.maketrans({'²': '**2', '³': '**3', '⁴': '**4', '⁻': '-', '¹': '**1', '⁰': '0'})


def normalise(line: str) -> str:
    """LaTeX, markdown and unicode maths into something parse_expr reads."""
    s = line
    s = s.replace('$', ' ').replace('`', ' ').replace('**', '^')      # markdown bold -> ^ first
    s = re.sub(r'\\(?:left|right|displaystyle|,|;|!|quad|qquad)', ' ', s)
    for _ in range(4):                                               # nested \frac
        s = re.sub(r'\\[dt]?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r'((\1)/(\2))', s)
    s = re.sub(r'\\sqrt\s*\{([^{}]*)\}', r'sqrt(\1)', s)
    s = re.sub(r'\\(sin|cos|tan|ln|log|exp|pi|arctan)\b', r'\1', s)
    s = re.sub(r'\\text\{[^{}]*\}|\\mathrm\{[^{}]*\}', ' ', s)
    s = (s.replace('\\cdot', '*').replace('\\times', '*').replace('×', '*').replace('·', '*')
          .replace('−', '-').replace('–', '-').replace('√', 'sqrt').replace('π', 'pi')
          .replace('≈', '=').replace('\\approx', '=').replace('{', '(').replace('}', ')'))
    s = s.translate(SUP)
    s = s.replace('^', '**')
    s = re.sub(r'(\d)\s*[eE]\s*([-+]?\d)', r'\1e\2', s)
    s = re.sub(r'(?<=\d)\s+x\s+(?=[\d(])', '*', s)                  # 121.0 x 10**6
    s = re.sub(r'(?<=\d),(?=\d{3}(?!\d))', '', s)                     # 1,000 -> 1000 only
    return s


CLAUSE = re.compile(r'=>|⇒|→|->|,\s|;|\band\b|\bso\b|\bthen\b|\bwith\b|\bwhere\b|\bthus\b|\bhence\b|:\s')


UNIT_WORDS = {'hours', 'hour', 'minutes', 'minute', 'seconds', 'second', 'meters',
              'metres', 'liters', 'litres', 'units', 'items', 'lots', 'degrees', 'kelvin',
              'newtons', 'joules', 'watts', 'volts', 'amperes', 'ohms', 'farads', 'hr', 'min'}


def split_unit(seg: str):
    """(expression, unit tail). A tail is at most 4 short tokens - not prose - and
    never begins with an exponent: `8.24e-05 C` keeps its e-05."""
    s = seg.strip().rstrip('.;:')
    m = re.match(r'^(.*?[\d\)])\s*(?![eE][-+]?\d)([A-Za-zµΩ°%][A-Za-z0-9µΩ°%/\*\-\s\(\)]*)$', s)
    if m:
        toks = re.findall(r'[A-Za-zµΩ°%]+', m.group(2))
        unitlike = all(len(t) <= 5 or t.lower() in UNIT_WORDS for t in toks) and len(toks) <= 4
        if toks and unitlike and not set(toks) <= FUNC_WORDS:
            return m.group(1).strip(), m.group(2).strip()
    return s.strip(), ''


def strip_units(seg: str) -> str:
    return split_unit(seg)[0]


FUNCS_DEG = dict(FUNCS, sin=lambda x: sympy.sin(x * sympy.pi / 180),
                 cos=lambda x: sympy.cos(x * sympy.pi / 180),
                 tan=lambda x: sympy.tan(x * sympy.pi / 180))


def evaluate(seg: str):
    """(candidate values, kind, unit). Trig is evaluated in radians and in degrees,
    because traces write cos(30) meaning degrees as often as radians; a claim holds
    if either reading makes it hold."""
    s, unit = split_unit(seg)
    vals = []
    for funcs in ((FUNCS, FUNCS_DEG) if re.search(r'\b(sin|cos|tan)\b', s) else (FUNCS,)):
        v, kind = _evaluate(s, funcs)
        if v is None:
            return [], kind, unit
        vals.append(v)
    return vals, 'num', unit


def _evaluate(s: str, funcs):
    if not s or len(s) > 200 or not re.search(r'\d', s):
        return None, 'symbolic' if re.search(r'[A-Za-z]', s or '') else 'unparseable'
    # An exponent is not a word: '1.21e8' must not read as the variable 'e8'.
    bare = re.sub(r'(?<=[\d.])[eE][-+]?\d+', '', s)
    words = set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*', bare)) - {'e', 'E'}
    if words - FUNC_WORDS:
        return None, 'symbolic'
    if re.search(r'\*\*\s*\(?\s*\d{3,}', s):                          # 10**1000: refuse
        return None, 'unparseable'
    try:
        expr = parse_expr(s, local_dict=funcs, transformations=TRANSFORMS, evaluate=True)
        if expr.free_symbols:
            return None, 'symbolic'
        v = complex(expr.evalf(15))
        if abs(v.imag) > 1e-9 * max(1.0, abs(v.real)) or not math.isfinite(v.real):
            return None, 'unparseable'
        return float(v.real), 'num'
    except Exception:                                                # noqa: BLE001
        return None, 'unparseable'


UNIT_FACTORS = (1.0, 60.0, 1 / 60.0, 3600.0, 1 / 3600.0, 1e3, 1e-3, 1e6, 1e-6, 1e9, 1e-9,
                1e12, 1e-12, 100.0, 0.01)


def agree(a: float, b: float, tol: float = ARITH_TOL) -> bool:
    if abs(b) < 1e-9:
        return abs(a) < 1e-6
    return abs(a - b) / abs(b) <= tol


def agree_any(la, ra, lu: str, ru: str) -> bool:
    """Any candidate pair agrees; under a unit factor only if the units differ."""
    factors = UNIT_FACTORS if (lu and ru and lu != ru) else (1.0,)
    return any(agree(a * f, b) for a in la for b in ra for f in factors)


@dataclass
class Claim:
    line: int
    left: str
    right: str
    left_value: list
    right_value: list
    ok: bool


@dataclass
class Report:
    claims: list = field(default_factory=list)
    segments: dict = field(default_factory=lambda: {'num': 0, 'symbolic': 0, 'unparseable': 0})

    @property
    def checked(self):
        return len(self.claims)

    @property
    def consistent(self):
        return sum(c.ok for c in self.claims)

    @property
    def rate(self):
        return self.consistent / self.checked if self.checked else None


def check(text: str) -> Report:
    rep = Report()
    for i, raw in enumerate(text.splitlines()):
        if '=' not in raw:
            continue
        for clause in CLAUSE.split(normalise(raw)):
            segs = re.split(r'(?<![<>!=])=(?!=)', clause)
            if len(segs) < 2:
                continue
            nums = []
            for sg in segs:
                v, kind, unit = evaluate(sg)
                rep.segments[kind] += 1
                if v:
                    nums.append((sg.strip(), v, unit))
            # `1/(41-23) hours = 60/18 = 3.33 minutes`: the unit is written once, at the
            # end, and belongs to every unit-less link before it.
            inherited, nxt = [], ''
            for sg, v, u in reversed(nums):
                nxt = u or nxt
                inherited.append((sg, v, u or nxt))
            nums = list(reversed(inherited))
            for (ls, lv, lu), (rs, rv, ru) in zip(nums, nums[1:]):
                rep.claims.append(Claim(i, ls[:80], rs[:80], lv, rv, agree_any(lv, rv, lu, ru)))
    return rep


def selftest() -> int:
    """Plants: a correct chain must pass, a corrupted one must fail."""
    cases = [
        ('Tr = 283.81 / 469.7 = 0.6042', True),
        ('V = (1.26 - 0.99) / 0.0353 = 7.6487 L', True),
        ('V = (1.26 - 0.99) / 0.0353 = 8.6487 L', False),
        (r'$Q = \frac{1}{0.013} \times 12 \times 1.2^{2/3} \times 0.001^{1/2} = 32.96$', True),
        (r'$Q = \frac{1}{0.013} \times 12 \times 1.2^{2/3} \times 0.001^{1/2} = 14.0$', False),
        ('F = 2.5 × 10^-6 × 3.0 = 7.5 × 10^-6 N', True),
        ('x = √(3² + 4²) = 5', True),
        ('x = √(3² + 4²) = 7', False),
        ('**Answer:** 113.55 cm³/mol', None),      # no claim to check
        ('q = -82.39 * 1e-6 = -8.24e-05 C', True),             # exponent is not a unit
        ('E = 200 GPa, L = 4.8 m', None),                      # two clauses, no claim
        ('x1 = 47.25 * cos(-140.8) = -36.62', True),           # degrees
        ('W = 1/(41 - 23) hours = 60/18 min', True),           # unit change
        ('p2 = 3.5% deteriorated case = 8.9%', None),          # prose is not a unit
        ('I = 121.0 x 10^6 mm^4 = 1.21e8 mm^4', True),         # x as multiplication
        ('I = 121.0 x 10^6 mm^4 = 1.2100e-04 m^4', True),      # mm^4 -> m^4 is 1e-12
        ('pi1 * 0.73 = 0.35  =>  pi1 = 0.35 / 0.73 = 0.4795', True),   # => is not =
        ('W1 = 1 / (41 - 23) hours = 60 / 18 = 3.33 minutes', True),   # inherited unit
        ('W1 = 1 / (41 - 23) hours = 60 / 18 = 4.33 minutes', False),  # ...and still checked
    ]
    bad = 0
    for text, want in cases:
        r = check(text)
        got = None if not r.checked else (r.consistent == r.checked)
        mark = 'ok' if got == want else 'FAIL'
        bad += got != want
        print('  [%s] %-70s expect %-5s got %s' % (mark, text[:70], want, got))
    print('selftest: %d failure(s)' % bad)
    return bad


if __name__ == '__main__':
    raise SystemExit(selftest())
