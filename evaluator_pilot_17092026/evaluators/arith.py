"""Check the arithmetic a text states, with sympy.

A derivation line like

    V_sat = 229.7 * 0.263**(0.396**0.2857) = 229.7 * 0.5715 = 131.27 cm³/mol

makes CLAIMS: each numeric segment equals the next. This module finds those
segments, evaluates each with sympy, and reports which consecutive pairs agree.
Symbolic segments (`Vc * Zc**(...)`) are skipped - they cannot be evaluated
without binding symbols - so a chain is checked from its first numeric segment on.

TWO RULES, TWO QUESTIONS (D-097). Every claim is judged twice and both verdicts are
carried on it, because they answer different questions:

  ok        relative tolerance ARITH_TOL = 1%. "Is this number FABRICATED?" Traces
            round every displayed intermediate to 3-4 significant figures and a
            chain of three roundings can drift ~1%, so 1% is the loosest reading of
            the shown work that still rejects a value the work cannot produce.
  ok_digit  the experts' own rule, from the annotation guide: *rounding is not an
            error, a wrong digit is*. Recompute the left side from the numbers the
            trace itself shows, and accept the right side only if it is a correct
            rounding AT THE PRECISION IT DISPLAYS: `0.92979^(2/3) = 0.95300` is
            rejected (it rounds to 0.95263 at five decimals), `= 0.9526` is not.

1% catches a fabricated number and is blind to every slip the experts marked - at 1%
the checker reaches recall 0.034 on the steps inside correct-answer traces, against
precision ~0.5 / recall ~0.5 for the digit rule (analysis/digit_rule.py, D-097). So
E4 scores its arithmetic on ok_digit and reports the 1% rate beside it; nothing that
used the 1% rule has lost it.

The digit rule keeps the unit handling: `91.4 deg = 1.594 rad` and `1% = 0.01` are
conversions, and the displayed precision is compared AFTER the unit factor. Without
that, unit conversions in the gold solutions are flagged as wrong digits - measured,
not assumed: analysis/arith_gold_validation.py.

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

GOLD IS NOT ENOUGH. At 100% on gold, the checker still flagged 15 milestones in
traces as "contradicted", and reading them showed every inspected one was a checker
bug, not a trace error. Gold is formatted uniformly by templates; traces are not. So
the rule is: validate on gold, then READ every flag raised on real traces. That
found four more:
  - `16t = 16 * 2.9`: a letter glued to a number is a variable, not a unit; only
    `%` and `°` may be glued. Nine of DeepSeek's flags were this one bug.
  - `(0.965)¹³` became `0.965**1**3`; superscript digits are translated as a run
  - `5.9%) = 0.4536`: a fragment with unbalanced parentheses is unparseable
  - `arctan(0.23383) = 13.16°`: inverse trig is also read in degrees
  - `64.834° - 180° = -115.166°`: the degree sign is removed as a unit mark, or the
    unit-tail rule swallows `° - 180°`
Sampling the NON-milestone inconsistent claims then found more, so per-claim
consistency was not yet evidence of anything:
  - `\\mathrm{m}^2` was deleted with its exponent left behind: `6.48 m^2` read as 41.99
  - Greek `μ` (what traces write) was not the micro sign `µ` the rules knew
  - `1% = 0.01` and `91.4 deg = 1.594 rad` are conversions, not contradictions
After those, the only flags left on real traces were GENUINE: two traces state the
right milestone while the arithmetic they show evaluates to something else (a sign
written wrong; a coefficient written wrong). Both are kept as plants, so a checker
change that stops catching them fails.

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
SUPDIGITS = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻', '0123456789-')


def _superscripts(s: str) -> str:
    """A run of superscript digits is ONE exponent: (0.965)¹³ is 0.965**13."""
    return re.sub(r'[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+', lambda m: '**(' + m.group(0).translate(SUPDIGITS) + ')', s)


def normalise(line: str) -> str:
    """LaTeX, markdown and unicode maths into something parse_expr reads."""
    s = line
    s = s.replace('$', ' ').replace('`', ' ').replace('**', '^')      # markdown bold -> ^ first
    s = re.sub(r'\\(?:left|right|displaystyle|,|;|!|quad|qquad)', ' ', s)
    for _ in range(4):                                               # nested \frac
        s = re.sub(r'\\[dt]?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r'((\1)/(\2))', s)
    s = re.sub(r'\\sqrt\s*\{([^{}]*)\}', r'sqrt(\1)', s)
    s = re.sub(r'\\(sin|cos|tan|ln|log|exp|pi|arctan)\b', r'\1', s)
    # \mathrm{m}^2 is a UNIT: keep the word so its exponent stays attached to it.
    # Deleting it left `6.48 **2` and read an area of 6.48 as 41.99.
    s = re.sub(r'\\(?:text|mathrm|rm|operatorname)\s*\{([^{}]*)\}', r' \1', s)
    s = (s.replace('\\cdot', '*').replace('\\times', '*').replace('×', '*').replace('·', '*')
          .replace('−', '-').replace('–', '-').replace('√', 'sqrt').replace('π', 'pi')
          .replace('≈', '=').replace('\\approx', '=').replace('{', '(').replace('}', ')'))
    # A degree sign marks a unit, not an operator: `64.834° - 180° = -115.166°` is
    # arithmetic on degrees. Left in, the unit-tail rule swallowed `° - 180°` whole.
    s = s.replace('°C', ' degC').replace('°F', ' degF').replace('°', ' ')
    s = s.replace('μ', 'u').replace('µ', 'u')        # Greek mu and micro sign alike
    s = _superscripts(s)
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
    m = re.match(r'^(.*?[\d\)])(\s*)(?![eE][-+]?\d)([A-Za-zµΩ°%][A-Za-z0-9µΩ°%/\*\-\s\(\)]*)$', s)
    if m and not m.group(2) and not m.group(3)[0] in '%°':
        # `16t`, `2y`, `9x`: a letter glued to a number is a variable times the
        # number, not a unit. Units are written apart (`16 m`), except % and °.
        return s.strip(), ''
    if m:
        m = re.match(r'^(.*?[\d\)])\s*(?![eE][-+]?\d)([A-Za-zµΩ°%][A-Za-z0-9µΩ°%/\*\-\s\(\)]*)$', s)
        toks = re.findall(r'[A-Za-zµΩ°%]+', m.group(2))
        # A unit carries digits only in an exponent (m**2, s**(-1)). Any other digit
        # means the "tail" is arithmetic - `deg - 180 deg` - and must not be dropped.
        stray_digit = re.search(r'\d', re.sub(r'\*\*\s*\(?\s*-?\d+\s*\)?', '', m.group(2)))
        unitlike = (all(len(t) <= 5 or t.lower() in UNIT_WORDS for t in toks)
                    and len(toks) <= 4 and not stray_digit)
        if toks and unitlike and not set(toks) <= FUNC_WORDS:
            return m.group(1).strip(), m.group(2).strip()
    return s.strip(), ''


def strip_units(seg: str) -> str:
    return split_unit(seg)[0]


FUNCS_DEG = dict(FUNCS, sin=lambda x: sympy.sin(x * sympy.pi / 180),
                 cos=lambda x: sympy.cos(x * sympy.pi / 180),
                 tan=lambda x: sympy.tan(x * sympy.pi / 180),
                 atan=lambda x: sympy.atan(x) * 180 / sympy.pi,
                 arctan=lambda x: sympy.atan(x) * 180 / sympy.pi,
                 asin=lambda x: sympy.asin(x) * 180 / sympy.pi,
                 acos=lambda x: sympy.acos(x) * 180 / sympy.pi)


def evaluate(seg: str):
    """(candidate values, kind, unit). Trig is evaluated in radians and in degrees,
    because traces write cos(30) meaning degrees as often as radians; a claim holds
    if either reading makes it hold."""
    raw = seg.strip()
    # Balance is checked BEFORE the unit tail is stripped: `3.5%)` loses its `)`
    # with the `%` and would otherwise pass as the claim 3.5 = 0.6296.
    if raw.count('(') != raw.count(')'):
        return [], 'unparseable', ''
    # `- 0.954**(49) = 0.0949/0.954`: a segment OPENING with a binary minus is the
    # tail of a subtraction cut at a clause break, not a negative number. A unary
    # minus is written against its operand (`-0.5`), so the space is the tell.
    if re.match(r'^-\s+\S', raw):
        return [], 'unparseable', ''
    s, unit = split_unit(seg)
    vals = []
    for funcs in ((FUNCS, FUNCS_DEG) if re.search(r'\b(a?sin|a?cos|a?tan|arctan)\b', s) else (FUNCS,)):
        v, kind = _evaluate(s, funcs)
        if v is None:
            return [], kind, unit
        vals.append(v)
    return vals, 'num', unit


def _evaluate(s: str, funcs):
    if s.count('(') != s.count(')'):
        return None, 'unparseable'          # a clause fragment, not an expression
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


def _unit_factors(lu: str, ru: str) -> set:
    """`1% = 0.01` and `91.4 deg = 1.594 rad` are conversions: a percent sign or an
    angle unit on EITHER side licenses the matching factor, because the bare side is
    simply unit-less."""
    factors = set(UNIT_FACTORS) if (lu and ru and lu != ru) else {1.0}
    units = (lu + ' ' + ru).lower()
    if '%' in units:
        factors |= {100.0, 0.01}
    if 'deg' in units or 'rad' in units:
        factors |= {math.pi / 180, 180 / math.pi}
    return factors


def agree_any(la, ra, lu: str, ru: str) -> bool:
    """Any candidate pair agrees within ARITH_TOL; under a unit factor only if the
    units differ."""
    return any(agree(a * f, b) for a in la for b in ra for f in _unit_factors(lu, ru))


NUM_LITERAL = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')
# A literal written with a decimal point or an exponent carries a rounding: it is a
# measured or displayed value, known only to its last digit. A bare integer (a count,
# an exponent, a coefficient) is exact and does not move.
ROUNDED_LITERAL = re.compile(r'(?<![A-Za-z0-9_.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def ulp(lit: str):
    """Place value of the last digit shown: '0.95300' -> 1e-5, '1.2e3' -> 100."""
    m = re.match(r'[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?$', lit.strip())
    if not m:
        return None
    return 10.0 ** ((int(m.group(3)) if m.group(3) else 0) - len(m.group(2)))


def displayed_ulp(seg: str):
    """The precision a segment DISPLAYS - the place value of the last digit of the
    first number it writes. None when no literal can be read, and a claim whose
    precision cannot be read is left unjudged by the digit rule rather than guessed
    at (the same choice the parser makes for a segment it cannot evaluate)."""
    m = NUM_LITERAL.search(seg.replace(',', ''))
    return ulp(m.group(0)) if m else None


def shown_uncertainty(seg: str, vals: list) -> list:
    """How far this segment's value can move if every ROUNDED number it shows is
    anywhere inside its own last digit - one figure per candidate value.

    Gold's `(-21) * (-5.06e-02) - (-27) * (-6.42e-02) = -0.6699` is why this exists.
    The gold states B to three figures and computes with the unrounded value, so the
    displayed operands reproduce -0.6708, not -0.6699. Judging the last digit of the
    result against operands that are themselves rounded flagged 12 correct gold
    claims (analysis/arith_gold_validation.py). The result can only be pinned as
    tightly as the numbers shown allow, so each literal is moved by half its own last
    digit and the resulting moves are added - first order, worst case, and correct
    under cancellation, where a relative tolerance is not.
    """
    out = [0.0] * len(vals)
    for m in ROUNDED_LITERAL.finditer(seg):
        lit = m.group(0)
        if '.' not in lit and 'e' not in lit and 'E' not in lit:
            continue                                  # a bare integer is exact
        h = ulp(lit)
        if h is None:
            continue
        try:
            x = float(lit)
        except ValueError:
            continue
        moved = seg[:m.start()] + '(' + repr(x + 0.5 * h) + ')' + seg[m.end():]
        vals2, kind, _ = evaluate(moved)
        if kind != 'num' or len(vals2) != len(vals):
            continue
        for k, v in enumerate(vals):
            out[k] += abs(vals2[k] - v)
    return out


def agree_any_digit(la, ra, lu: str, ru: str, u: float, ul=None, ur=None) -> bool:
    """The experts' rule: the right side, as displayed, is a correct rounding of the
    left side, computed from the numbers shown.

    The tolerance is half the place value of the displayed last digit, widened to
    whatever the shown operands leave undetermined (`shown_uncertainty`) on either
    side. `ul`/`ur` default to zero, which is the rule as D-097 measured it; the
    caller adds them only for a claim the bare rule rejects, since they can only
    widen it. A relative slack of 1e-7 keeps a value that lands exactly on the
    rounding boundary in binary floating point off the list.
    """
    ul = ul or [0.0] * len(la)
    ur = ur or [0.0] * len(ra)
    for f in _unit_factors(lu, ru):
        for i, a in enumerate(la):
            for j, b in enumerate(ra):
                tol = max(0.5 * u, ur[j]) + ul[i] * abs(f)
                if abs(a * f - b) <= tol * (1 + 1e-7):
                    return True
    return False


@dataclass
class Claim:
    line: int
    left: str
    right: str
    left_value: list
    right_value: list
    ok: bool                      # within ARITH_TOL (1%): the fabrication test
    ok_digit: bool = True         # a correct rounding at the precision shown (D-097)
    ulp: float = None             # the place value that judged it; None = unjudged
    left_unit: str = ''           # the unit tails, as split_unit read them, so an
    right_unit: str = ''          # audit can re-judge the claim under another rule


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

    @property
    def consistent_digit(self):
        return sum(c.ok_digit for c in self.claims)

    @property
    def rate_digit(self):
        return self.consistent_digit / self.checked if self.checked else None

    @property
    def digit_judged(self):
        """Claims whose displayed precision could be read, so the digit rule applied."""
        return sum(c.ulp is not None for c in self.claims)


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
                u = displayed_ulp(rs)
                if u is None:
                    digit = True                      # no readable precision: unjudged
                else:
                    digit = agree_any_digit(lv, rv, lu, ru, u)
                    if not digit:                     # only then pay for the propagation
                        digit = agree_any_digit(lv, rv, lu, ru, u,
                                                shown_uncertainty(ls, lv),
                                                shown_uncertainty(rs, rv))
                rep.claims.append(Claim(i, ls[:80], rs[:80], lv, rv,
                                        agree_any(lv, rv, lu, ru), digit, u, lu, ru))
    return rep


def selftest() -> int:
    """Plants: a correct chain must pass, a corrupted one must fail.

    Each case carries what the 1% rule must say and, where the two differ, what the
    digit rule must say (a third element). The rules disagree by design: the digit
    rule rejects a last digit that is not a correct rounding even when the number is
    within 1%, and that is the whole point of it.
    """
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
        ('dv/dt = 16t = 16 * (2.9) = 46.4', True),              # 16t is a variable
        ('Pa = (0.965)¹³ = 0.6293', True),                      # superscript run
        ('Pa = (0.965)¹³ = 0.9650', False),                     # ...and still checked
        ('phi = arctan(0.23383) = 13.16°', True),               # inverse trig, degrees
        ('p (at 5.9%) = 0.4536', None),                         # fragment, no claim
        ('phi = 64.834° - 180° = -115.166°', True),             # degree sign is a unit
        (r'$A = 2.7 \times 2.4 = 6.48 \mathrm{m}^2$', True),     # LaTeX unit keeps its exponent
        ('q = -82.39 μC = -82.39 × 10^-6 C', True),             # Greek mu
        # deg -> rad is a conversion, not a contradiction. 91.4 deg is 1.59523 rad, but
        # 91.4 is itself shown to a tenth: 91.35-91.45 deg covers 1.594, so it stands.
        ('theta = 91.4 deg = 1.594 rad', True),
        ('theta = 91.400 deg = 1.594 rad', True, False),        # shown to a thousandth: wrong digit
        ('p = 1% = 0.01', True),                                # percent -> fraction
        ('p = 1% = 0.02', False),                               # ...and still checked
        ('Pa = (1 - Pd) (at 3.5%) = 0.6296', None),             # fragment `3.5%)`
        ('- 0.954^49 = 0.0995', None),                          # clause-cut `- 0.954^49`...
        ('P = 1 - 0.954^49 = 0.0949', False),                   # a whole expression IS checked
        ('x = -0.5 * 2 = -1.0', True),                          # ...but unary minus is fine
        # Real slips found by auditing a sample of claims (all llama-3.1-70b):
        ('Pa = (1 - 0.059)^13 = 0.5786', False),
        ('V = 0.276^1.1745 = 0.4439', False),
        # Two slips found in real traces (llama, gpt-5 on normal_depth_iteration#1):
        # the right milestone is stated, the arithmetic shown does not produce it.
        ('A = (3.8 + 2*1.500*2)*1.500 = 14.100 m^2', False),
        ('y = 1.853 - (0.2892*(-0.160))/(-2.8172) = 1.86921', False),
        # THE DIGIT RULE (D-097). A correct rounding passes at any precision; a last
        # digit the shown numbers do not produce fails, even though 1% forgives it.
        ('k = 0.92979^(2/3) = 0.9526', True, True),
        ('k = 0.92979^(2/3) = 0.95300', True, False),     # rounds to 0.95263
        ('A = 0.4525 * 0.089 * 80 / 100 = 0.032218', True, True),
        # The same wrong sixth digit, twice: forgiven when the operands shown cannot
        # pin it (0.089 is two figures), flagged when they can (0.0890 is three).
        ('A = 0.4525 * 0.089 * 80 / 100 = 0.032318', True, True),
        ('A = 0.4525 * 0.0890 * 80 / 100 = 0.032318', True, False),
        ('r = 2/3 = 0.67', True, True),                   # rounding up is a rounding
        ('r = 2/3 = 0.66', False, False),                 # truncation is a wrong digit
        ('Re = 1000 * 2.5 * 0.05 / 0.00089 = 140449', True, True),
        # `0.05` shown to two decimals is 0.045-0.055, so nothing past the second
        # figure of Re is pinned at all; write the operands out and the digit is.
        ('Re = 1000 * 2.500 * 0.0500 / 0.000890 = 141450', True, False),
    ]
    bad = 0
    for case in cases:
        text, want = case[0], case[1]
        want_digit = case[2] if len(case) > 2 else want
        r = check(text)
        got = None if not r.checked else (r.consistent == r.checked)
        got_digit = None if not r.checked else (r.consistent_digit == r.checked)
        miss = (got != want) + (got_digit != want_digit)
        bad += bool(miss)
        print('  [%s] %-62s 1%%: expect %-5s got %-5s | digit: expect %-5s got %s'
              % ('ok' if not miss else 'FAIL', text[:62], want, got, want_digit, got_digit))
    print('selftest: %d failure(s)' % bad)
    return bad


if __name__ == '__main__':
    raise SystemExit(selftest())
