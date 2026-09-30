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

READING THE FULL RUN'S FLAGS (D-154, D-156). With gold still at no flag, a domain expert read 220
of the full run's flags and found 109 to be the checker's. Each rule below answers some of the notes:
  - a claim ends at `\\Rightarrow`, `\\implies`, `\\to`, `\\qquad`, an aligned row's `\\\\`, a full stop,
    `Then`, `Thus`, `Therefore` or `gives` in any case, and at `<`, `>`, `≤`, `≥` or `≠`; two spaces
    before a new equation whose left side carries ∂, _, [ ] or a Greek letter separate the two
  - thousands grouped by a space, `\\,`, `{,}` or a thin space are one number: `13 855`, `226\\,580`
  - `a / 2.1506 × 10⁵` divides by the whole literal, and a literal in scientific notation shows its
    precision with its exponent: `4.0007 × 10⁷` is known to 1e3, not to 1e-4
  - single letters that are not unit symbols are variables: `48 E I`, `2 D S / H`, `8 μ L`, `(5/8) P`
    and `d²` are symbolic, not the numbers 48, 2, 8 and 0.625 in the units `E I`, `D S / H`, `μ L`, `P`
  - `0.916-0.917` is a range, not a subtraction; `39.5967...` is a truncation: the value lies in the
    unit after the digits shown, and a number cut short inside a chain is read as the number
  - `16sin(120°)`: trig glued to its coefficient is read in degrees as well
  - `**4.301 minutes**` is markdown bold, so its unit reaches the chain; `customers/min` is a unit;
    a unit tail never closes a bracket (`sin(75.35 rad)` is not `sin(75.35` in rad)
  - a line after one that ends in an operator, or a segment opening with `+(`, continues an expression
    begun above; and a segment the parser cannot read ends a chain rather than joining its neighbours
    as though they were claimed equal

THE SECOND READING (D-158, D-159). The same expert read 206 flags of the fixed rule: 155 real,
0.752. The rules below answer 42 of the 51 notes, and the audit (DIGIT_FIX_2.md) moves no pilot step
away from the experts:
  - a table cell ends at `&` or a spaced `|`, and `or`, `but`, `since`, `unless`, `using`, `versus`,
    `would` and `give` end a claim, as `is` does before an assignment; after a comparison, a bare
    number is its threshold (`P_detect≤2 = 0.3086`), while an expression is still a claim
  - a word in any script is a variable (`ẋ(0)`), and so is a Greek letter's name (`8 \\mu L`), save
    `mu` before a unit symbol; S is not siemens
  - a side that writes no unit of its own may convert into the unit the other side writes, by the
    factor that unit implies (`mm` from m, `kN` from N, minutes from hours, ppm); a base unit implies
    nothing, and no conversion makes a zero; Stokes, bar, atm and hectares convert for their own pairs
  - a trailing `...` or `\\dots` accepts a correct rounding as well as a truncation, and an operand cut
    short may lie a whole unit higher; `a/2π` is a/(2π); `99° = ...` keeps its unit; a sum written a
    term a line continues
Left as they are: running computations, carried roundings, units the model left unstated where no
unit on the right implies them, ratio notation, a sign convention, a quadrant, and `wL`, `Zc`.

WHAT IT CANNOT PARSE IT SAYS SO. Every segment is classified: evaluated, symbolic
(skipped), or unparseable. A checker that silently parses 5% of lines and reports
"100% consistent" would be this repo's recurring defect, so coverage is reported
alongside every consistency figure.
"""
from __future__ import annotations

import math
import re
import warnings
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


BOLD = re.compile(r'(?<![\w)\]])\*\*(?=\S)(.+?)(?<=\S)\*\*(?![\w(])')
GROUP_SEP = r'(?:\\,|\{,\}|\\ |[ \u2009\u202f\u00a0])'
THOUSANDS = re.compile(r'(?<![\d.])\d{1,3}(?:' + GROUP_SEP + r'\d{3})+(?!\d)')
IMPLIES = re.compile(r'\\(?:Rightarrow|Longrightarrow|implies|iff|Leftrightarrow|Longleftrightarrow|therefore'
                     r'|rightarrow|longrightarrow|to)(?![A-Za-z])|[⟹⟶⇔⟺∴]')
RELATIONS = ((r'\\(?:leq?|leqslant)(?![A-Za-z])', ' ≤ '), (r'\\(?:geq?|geqslant)(?![A-Za-z])', ' ≥ '),
             (r'\\neq?(?![A-Za-z])', ' ≠ '), (r'\\lt(?![A-Za-z])', ' < '), (r'\\gt(?![A-Za-z])', ' > '))
# Scientific notation has a decimal mantissa (`2.1506 × 10⁵`). An integer times a power of ten is
# arithmetic: gold's `(92 - 20)/92 * 10^6` is ((92 - 20)/92) × 10⁶, a conversion to ppm (D-156).
SCI = re.compile(r'(?<![\w.])(\d+\.\d+)\s*\*\s*10\s*\*\*\s*(?:\(\s*([-+]?\d+)\s*\)|([-+]?\d+))')


def normalise(line: str) -> str:
    """LaTeX, markdown and unicode maths into something parse_expr reads."""
    s = BOLD.sub(r'\1', line)                        # `**4.301 minutes**` is bold, not a power (D-156)
    s = THOUSANDS.sub(lambda m: re.sub(GROUP_SEP, '', m.group(0)), s)    # 13 855, 226\,580, 8{,}212 (D-156)
    s = IMPLIES.sub(' => ', s)                       # \Rightarrow, \implies, \to end a claim (D-156)
    for pat, rel in RELATIONS:
        s = re.sub(pat, rel, s)
    s = s.replace('\\qquad', ' ; ').replace('\\\\', ' ; ')    # side by side; aligned rows
    # A table cell ends at `&` (LaTeX) or a spaced `|` (markdown); `&=` only aligns (D-159).
    s = re.sub(r'&(?!\s*=)', ' ; ', s).replace('&', ' ')
    s = re.sub(r'(?<!\S)\|(?!\S)', ' ; ', s)
    s = re.sub(r'\\[()\[\]]', ' ', s)               # the \( \) \[ \] delimiters of inline maths
    s = s.replace('…', '...').replace('\\%', '%').replace('\\Omega', 'Ω')
    s = re.sub(r'\\[lc]?dots(?![A-Za-z])', '...', s)   # `0.969047 \dots` is cut short (D-159)
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
    s = s.replace('°C', ' degC').replace('°F', ' degF')
    # ...but a degree sign on the side of an equation is that side's unit: `99° = 11π/20 = 1.7279 rad`
    # is a conversion (D-159).
    s = re.sub(r'(\d)\s*°(?=\s*=)', r'\1 deg', s).replace('°', ' ')
    s = s.replace('μ', 'u').replace('µ', 'u')        # Greek mu and micro sign alike
    s = re.sub(r'/\s*(\d+(?:\.\d+)?)\s*pi\b', r'/(\1*pi)', s)    # `267π/2π` is 267π/(2π) (D-159)
    s = _superscripts(s)
    s = s.replace('^', '**')
    s = re.sub(r'(\d)\s*[eE]\s*([-+]?\d)', r'\1e\2', s)
    s = re.sub(r'(?<=\d)\s+x\s+(?=[\d(])', '*', s)                  # 121.0 x 10**6
    s = re.sub(r'(?<=\d),(?=\d{3}(?!\d))', '', s)                     # 1,000 -> 1000 only
    s = re.sub(r'\\([A-Za-z]+)', r'\1', s)          # any other command is its name: \Phi(x) -> Phi(x)
    # `2.1506 × 10⁵` is ONE number: `a / 2.1506 × 10⁵` divides by all of it, and its last digit is
    # worth 10, not 1e-4 (D-156).
    s = SCI.sub(lambda m: '(' + m.group(1) + 'e' + (m.group(2) or m.group(3)) + ')', s)
    return s


CLAUSE = re.compile(r'=>|⇒|→|->|,\s|;|:\s'
                    r'|(?i:\b(?:and|so|then|with|where|thus|hence|therefore|gives|giving|yields)\b)'
                    # Prose that sets two values side by side: `X = 0 or X = 1`, `Cpk = 1.67 since
                    # Cp = 1.14`, `803.16 units using Q = 1090`, `σ = 17.17 Ω would give Cpl = ...` (D-159).
                    r'|(?i:\b(?:or|but|since|unless|using|versus|vs|would|give)\b)'
                    # A full stop ends a sentence, and two spaces before a new equation whose left side
                    # carries ∂, _, [ ] or a Greek letter part two (D-156).
                    r'|(?<!\.)\.(?!\.)(?=\s|$)'
                    r'|\s{2,}(?=[^\s=]*[∂_\[\]\u0391-\u03a9\u03b1-\u03c9][^\s=]*\s*=(?!=))'
                    # A connector that opens a new assignment ends the claim before it:
                    # `34 kN at a = 3.3 m`, `D3 = 0 for n = 4`, `from t = -2 to t = 2`
                    # are lists of values, not the chain 34 = 3.3 (D-120); so does `is` (D-159).
                    r'|\s(?:at|under|to|for|from|over|between|when|if|on|by|is|are)\s+'
                    r'(?=[A-Za-z_][\w\[\]()\-]*\s*=(?!=))')
# A comparison is not an equation (D-156), and what follows one is its right side, not a value:
# `P_detect≤2 = 0.3086` states P(detect ≤ 2), not 2 = 0.3086 (D-159).
RELATION = re.compile(r'[<>≤≥≠]')
# An '=' inside a label - `[1/(-r_A)]_at_X=0.85 = 4.55` - is not an equation when a spaced
# '=' follows it (D-120).
LABEL_EQ = re.compile(r'(?<=[\w\])])=(?=[-+]?\d)(?=.*\s=\s)')


UNIT_WORDS = {'hours', 'hour', 'minutes', 'minute', 'seconds', 'second', 'meters',
              'metres', 'liters', 'litres', 'units', 'items', 'lots', 'degrees', 'kelvin',
              'newtons', 'joules', 'watts', 'volts', 'amperes', 'ohms', 'farads', 'hr', 'min',
              # what a queue or a line counts, per unit of time: `0.25 customers/min` (D-156)
              'customers', 'customer', 'orders', 'vehicles', 'patients', 'people', 'workers',
              'machines', 'batches', 'pieces', 'packets', 'symbols', 'cycles'}
# The single letters that are unit symbols. Any other single letter is a variable: `48 E I` is 48EI,
# not 48 in the unit `E I` (D-156). d and t are left out: traces write them as diameter and time far
# more often than as day and tonne; so is S (D-159), a retention, a spacing or a setup cost far more
# often than siemens.
UNIT_LETTERS = set('mgNJWVAKCFHTLlshM') | {'Ω', '%', '°'}
# A Greek letter's name is a variable too: `8 \mu L` is 8μL, viscosity times length (D-159). Only `mu`
# before a unit symbol is the micro prefix: `5 \mu m`, `2 \mu s`.
# psi is left out: as a unit (pounds per square inch) it is far commoner in these traces than ψ.
GREEK = {'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'varepsilon', 'zeta', 'eta', 'theta', 'vartheta',
         'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'varphi',
         'chi', 'omega', 'Gamma', 'Delta', 'Theta', 'Lambda', 'Xi', 'Sigma', 'Phi'}
MICRO = {'m', 's', 'g', 'F', 'A', 'V', 'W', 'J', 'H', 'N', 'C', 'T', 'l', 'mol', 'Pa'}


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
        if toks and all(len(t) == 1 for t in toks) and not set(toks) <= UNIT_LETTERS:
            return s.strip(), ''          # `2 D S / H`, `8 μ L`, `(5/8) P`: variables (D-156)
        if m.group(1).count('(') != m.group(1).count(')'):
            return s.strip(), ''          # `sin(75.35 rad)`: a tail that closes a bracket is not a unit (D-156)
        if set(toks) & FUNC_WORDS:
            return s.strip(), ''          # `2 sqrt(m k)`: a tail holding a function is not a unit (D-156)
        if set(toks) & GREEK and not (toks[0] == 'mu' and len(toks) > 1 and toks[1] in MICRO):
            return s.strip(), ''          # `8 mu L`, `2 mu`: Greek letters are variables (D-159)
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


# `16sin(120)` is trig too: a coefficient glued on leaves no word boundary before `sin` (D-156).
TRIG = re.compile(r'(?<![A-Za-z])(?:a?sin|a?cos|a?tan|arctan)(?![A-Za-z])')
RANGE = re.compile(r'^\s*(\d+\.\d+)\s*-\s*(\d+\.\d+)\s*$')


def is_range(lo: str, hi: str) -> bool:
    """`0.916-0.917`, `2.13 - 2.14`: two numbers shown to the same place, the second one or two units
    above the first, are a range the value lies in. A subtraction is not written that way (D-156)."""
    u = ulp(lo)
    return (u is not None and u == ulp(hi) and 0 < float(hi) - float(lo) <= 2 * u * (1 + 1e-9))


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
    # `+(0.0121-1)=0`: a segment opening with `+(` or `+ ` is the last term of a sum begun on
    # an earlier line; a signed number (`+0.162`) is written against its digits (D-156).
    if re.match(r'^\+\s*[(\s]', raw):
        return [], 'unparseable', ''
    # `sqrt(941713.24...)`, `96.9047 ... %`: a number cut short is a number, with its unit after it.
    s, unit = split_unit(re.sub(r'(?<=[\d)])\s*\.\.\.', '', seg))
    m = RANGE.match(s)
    if m and is_range(m.group(1), m.group(2)):
        return [], 'unparseable', unit
    vals = []
    for funcs in ((FUNCS, FUNCS_DEG) if TRIG.search(s) else (FUNCS,)):
        v, kind = _evaluate(s, funcs)
        if v is None:
            return [], kind, unit
        vals.append(v)
    return vals, 'num', unit


def _evaluate(s: str, funcs):
    if s.count('(') != s.count(')'):
        return None, 'unparseable'          # a clause fragment, not an expression
    if not s or len(s) > 200 or not re.search(r'\d', s):
        # `σ/ε` is a formula in any script, and a formula links the chain (D-159)
        return None, 'symbolic' if re.search(r'[^\W\d_]', s or '') else 'unparseable'
    # An exponent is not a word: '1.21e8' must not read as the variable 'e8'. A word may be written
    # in any script: `ẋ(0)` is a variable, where sympy reads a symbol times zero, 0 (D-159).
    bare = re.sub(r'(?<=[\d.])[eE][-+]?\d+', '', s)
    words = set(re.findall(r'[^\W\d]\w*', bare)) - {'e', 'E'}
    if words - FUNC_WORDS:
        return None, 'symbolic'
    if re.search(r'\*\*\s*\(?\s*\d{3,}', s):                          # 10**1000: refuse
        return None, 'unparseable'
    try:
        with warnings.catch_warnings():             # `y[n](x)`: sympy warns while it refuses a fragment
            warnings.simplefilter('ignore', SyntaxWarning)
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
                1e12, 1e-12, 100.0, 0.01,
                12.0, 1 / 12.0)          # ft and in, the US-unit templates (D-120)
# 1e4 (ha and m²) was tried and left out: it passed a real slip, 4.036e-5 m³/mol "=" 0.4036 cm³/mol,
# for the one hectare claim it would have read (D-156).


def agree(a: float, b: float, tol: float = ARITH_TOL) -> bool:
    if abs(b) < 1e-9:
        return abs(a) < 1e-6
    return abs(a - b) / abs(b) <= tol


# Conversions no generic factor covers, for these unit pairs only: 1e4 for hectares here cannot pass
# the m³-and-cm³ slip that ruled it out as a general factor (D-156, D-159).
PAIR_FACTORS = {('m**2/s', 'St'): 1e4, ('St', 'm**2/s'): 1e-4, ('bar', 'Pa'): 1e5, ('Pa', 'bar'): 1e-5,
                ('bar', 'kPa'): 1e2, ('kPa', 'bar'): 1e-2, ('bar', 'MPa'): 0.1, ('MPa', 'bar'): 10.0,
                ('atm', 'Pa'): 101325.0, ('Pa', 'atm'): 1 / 101325.0, ('atm', 'kPa'): 101.325,
                ('kPa', 'atm'): 1 / 101.325, ('ha', 'm**2'): 1e4, ('m**2', 'ha'): 1e-4}


def unit_key(u: str) -> str:
    return re.sub(r'\*\*\((\d)\)', r'**\1', re.sub(r'\s+', '', u))


SI_PREFIX = {'k': 1e3, 'M': 1e6, 'G': 1e9, 'c': 1e-2, 'm': 1e-3, 'u': 1e-6, 'n': 1e-9, 'p': 1e-12}
PREFIXABLE = {'m', 'N', 'J', 'W', 'V', 'A', 'Pa', 'Hz', 'L', 'l', 's', 'C', 'F', 'H', 'T', 'Ω', 'eV', 'Wh'}
TIME_FACTORS = {'min': {60.0, 1 / 60}, 'minute': {60.0, 1 / 60}, 'minutes': {60.0, 1 / 60},
                'h': {1 / 60, 1 / 3600}, 'hr': {1 / 60, 1 / 3600}, 'hour': {1 / 60, 1 / 3600},
                'hours': {1 / 60, 1 / 3600}, 'sec': {60.0, 3600.0}, 'seconds': {60.0, 3600.0}}


def implied_factors(ru: str) -> set:
    """What a unit written only on the right implies for a left side that writes none, computed in the
    base unit (D-159): `mm` ×1e3 from m, `kN` ×1e-3 from N, `MPa` ×1e-6 from Pa, `uV` ×1e6 from V;
    minutes from hours or seconds; ppm ×1e6. A base unit implies nothing, so `(338.86e6)/(122e6) =
    2.7775 × 10⁻³ m` stays a slip, and no generic factor can pass a wrong power of ten."""
    m = re.match(r'[A-Za-zΩ]+', unit_key(ru))
    tok = m.group(0) if m else ''
    if tok in TIME_FACTORS:
        return TIME_FACTORS[tok]
    if tok in ('ppm', 'ppb'):
        return {1e6 if tok == 'ppm' else 1e9}
    if len(tok) > 1 and tok[0] in SI_PREFIX and tok[1:] in PREFIXABLE:
        return {1 / SI_PREFIX[tok[0]]}
    return set()


def _unit_factors(lu: str, ru: str, implied: set | None = None) -> set:
    """`1% = 0.01` and `91.4 deg = 1.594 rad` are conversions: a percent sign or an
    angle unit on EITHER side licenses the matching factor, because the bare side is
    simply unit-less. `implied`: the factors the right side's unit implies when the left
    side writes no unit of its own: `0.05243/73 = 0.718 mN` (D-159)."""
    factors = set(UNIT_FACTORS) if (lu and ru and lu != ru) else {1.0}
    factors |= implied or set()
    pair = PAIR_FACTORS.get((unit_key(lu), unit_key(ru)))
    if pair:
        factors.add(pair)
    units = (lu + ' ' + ru).lower()
    if '%' in units:
        factors |= {100.0, 0.01}
    if 'deg' in units or 'rad' in units:
        factors |= {math.pi / 180, 180 / math.pi}
    return factors


def agree_any(la, ra, lu: str, ru: str, implied: set | None = None) -> bool:
    """Any candidate pair agrees within ARITH_TOL; under a unit factor only if the
    units differ, or the left side's is implied. A conversion never makes a zero: `131 - 79 =
    0 s` is not 52 × 1e-12 (D-159)."""
    return any(agree(a * f, b) for a in la for b in ra for f in _unit_factors(lu, ru, implied)
               if f == 1.0 or b != 0.0)


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
    s, _ = split_unit(re.sub(r'(?<=[\d)])\s*\.\.\.', '', seg))
    bare = bool(re.fullmatch(r'\s*[-+]?\(?\s*(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\s*\)?\s*', s))
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
        # An operand cut short (`0.969047... × 100`) may lie a whole unit of its last digit higher. A
        # result cut short is judged as such by check_part, so it moves half a unit like any other (D-159).
        step = h if bare is False and re.match(r'\s*\.\.\.', seg[m.end():]) else 0.5 * h
        moved = seg[:m.start()] + '(' + repr(x + step) + ')' + seg[m.end():]
        vals2, kind, _ = evaluate(moved)
        if kind != 'num' or len(vals2) != len(vals):
            continue
        for k, v in enumerate(vals):
            out[k] += abs(vals2[k] - v)
    return out


def agree_any_digit(la, ra, lu: str, ru: str, u: float, ul=None, ur=None, implied: set | None = None) -> bool:
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
    for f in _unit_factors(lu, ru, implied):
        for i, a in enumerate(la):
            for j, b in enumerate(ra):
                if f != 1.0 and b == 0.0:
                    continue                          # a conversion never makes a zero (D-159); 1e-12 m⁴ is not one
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


# A line after one that ends in an operator continues an expression begun above: `+` /
# `\frac{24(4.7)^2}{12}` / `-` / `\frac{24(3.0)^2}{12}=0` is one equation, not the claim 18 = 0 (D-156).
CONTINUED = re.compile(r'(?:[-+*/×·]|\\times|\\cdot|\\pm)\s*$')


def check(text: str) -> Report:
    rep = Report()
    prev = ''
    for i, raw in enumerate(text.splitlines()):
        # A line continues the one above when that one ends in an operator (D-156), or when both are
        # terms of one sum written a term a line, `+(9.47839)V_m` / `-(9.47839)(0.065898)=0` (D-159).
        continued = bool(CONTINUED.search(prev) or (re.match(r'\s*[-+](?=[(\d\\])', prev)
                                                     and re.match(r'\s*[-+](?=[(\\])', raw)))
        if raw.strip():
            prev = raw
        if '=' not in raw:
            continue
        for ci, clause in enumerate(CLAUSE.split(normalise(raw))):
            for ri, part in enumerate(RELATION.split(clause)):
                segs = re.split(r'(?<![<>!=])=(?!=)', LABEL_EQ.sub('≡', part))
                if len(segs) < 2:
                    continue
                # After a comparison, a bare number is its threshold, `P_detect≤2 = 0.3086`; an
                # expression is a claim, `W_L ≥ 824/179 = 4.61` (D-159).
                threshold = ri > 0 and bool(re.fullmatch(r'\s*[-+]?\d+(?:\.\d+)?\s*', segs[0]))
                check_part(rep, i, segs, (continued and ci == 0 and ri == 0) or threshold)
    return rep


def check_part(rep: Report, i: int, segs: list, first_is_not_a_value: bool) -> None:
    """The claims of one clause, cut at its `=`s. `first_is_not_a_value`: the first segment is the tail
    of an expression begun above, or the right side of a comparison (D-156, D-159)."""
    # A segment the parser cannot read ends the chain: `735 (Shortage 757-735 = 22` does not
    # claim 735 = 22, and `0.336-0.337` does not join its neighbours. A symbolic segment (a
    # formula) still links them: `12.5 = P/A = 12.5` claims the two values equal (D-156).
    chains, nums = [], []
    for k, sg in enumerate(segs):
        if first_is_not_a_value and k == 0:
            v, kind, unit = [], 'unparseable', ''
        else:
            v, kind, unit = evaluate(sg)
        rep.segments[kind] += 1
        if v:
            nums.append((sg.strip(), v, unit))
        elif kind == 'unparseable':
            chains.append(nums)
            nums = []
    chains.append(nums)
    for nums in chains:
        # `1/(41-23) hours = 60/18 = 3.33 minutes`: the unit is written once, at the
        # end, and belongs to every unit-less link before it. A link that writes no unit of its
        # own may also convert into the next one's: `0.05243/73 = 0.718 mN` (D-159).
        inherited, nxt = [], ''
        for sg, v, u in reversed(nums):
            nxt = u or nxt
            inherited.append((sg, v, u or nxt, bool(u)))
        nums = list(reversed(inherited))
        for (ls, lv, lu, lown), (rs, rv, ru, rown) in zip(nums, nums[1:]):
            implied = implied_factors(ru) if rown and not lown else set()
            u = displayed_ulp(rs)
            # `39.5967...` is cut short: the value lies in the unit AFTER the digits shown,
            # [39.5967, 39.5968), judged as a rounding of that unit's middle; a correct rounding
            # passes too, since traces also write `...` after one (D-156, D-159). `0.03461552...`
            # for 0.034615511 is a wrong digit either way.
            cut = u is not None and bool(re.search(r'[\d)]\s*\.\.\.', rs))
            rj = [b + math.copysign(0.5 * u, b) for b in rv] if cut else None
            if u is None:
                digit = True                              # no readable precision: unjudged
            else:
                digit = any(agree_any_digit(lv, r, lu, ru, u, implied=implied) for r in (rv, rj) if r)
                if not digit:                             # only then pay for the propagation
                    ul, ur = shown_uncertainty(ls, lv), shown_uncertainty(rs, rv)
                    digit = any(agree_any_digit(lv, r, lu, ru, u, ul, ur, implied=implied)
                                for r in (rv, rj) if r)
            rep.claims.append(Claim(i, ls[:80], rs[:80], lv, rv,
                                    agree_any(lv, rv, lu, ru, implied), digit, u, lu, ru))


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
        # READING THE FULL RUN'S FLAGS (D-156). Each line is a form from the domain expert's notes
        # that was flagged on correct arithmetic before the rule that now reads it.
        (r'\sqrt{M\,T} = \sqrt{33\,961.1} \approx 184.3', True),               # thousands, \,
        ('850 × 16.3 = 13 855', True),                                        # thousands, a space
        (r'q = 8.314 \times 2849.24 = 23{,}688.6 \text{ J/mol} = 23.689 \text{ kJ/mol}', True),
        (r'-3-0.66(5)=-6.30 \Rightarrow \Phi(-6.30)=0.0000', True),           # \Rightarrow ends a claim
        (r'$-1 = B(-3) \Rightarrow B = \frac{1}{3}$', None),                  # ...so -1 is not B
        (r'L = \sqrt{256} = 16. The number of points is n = L/2 = 8.', True), # a full stop ends one
        ('beta = 0.8665-0.0000 = 0.8665. Then ARL = 1/(1-0.8665) = 7.4906', True),
        ('t = 0.5542 / 1.461×10⁻⁶ ≈ 3.793×10⁵ min', True),                   # a×10^n is one number
        ('t = 0.5542 / 1.461×10⁻⁶ ≈ 3.799×10⁵ min', True, False),            # ...known to 1e2
        (r'N = \frac{118,038}{2,981} \approx 39.5967...', True),              # cut short, not rounded
        (r'\sigma = \sqrt{0.0011982336} = 0.03461552...', True, False),       # ...but 0.034615511
        ('y2 = 0.31 × 6.88 ≈ 2.13 – 2.14 m', None),                           # a range, not 2.13 - 2.14
        (r'48 E I = 48 \times 29000 \times 5900 = 8212800000', True),         # E and I are variables
        ('2 D S / H = 15,300,088 / 6.5778 ≈ 2,326,019', True),                # so are D, S and H
        (r'$$(1.9)^2 = 3.61 \qquad (1.9)^4 = (3.61)^2 = 13.0321$$', True),    # \qquad parts two
        ('f_s = 1/(1.03 × 10⁻³) = 970.87 Hz > 2f₀ = 2 × 162 = 324 Hz', True), # > is not =
        ('T0 = 0 + 4 = 4 ≠ 0', True),                                         # nor is ≠
        (r'g = 16\sin(120°) = 16 \times 0.8660 = 13.856', True),              # glued trig, degrees
        ('W = 0.07168 hours = 0.07168 × 60 = **4.301 minutes**.', True),      # bold keeps its unit
        ('∂v/∂x = 8(0.9) = 7.2 s⁻¹  ∂v/∂y = 8(3.3) = 26.4 s⁻¹', True),       # side by side
        ('mu = 1 / 4 min = 0.25 customers/min = 0.25 × 60 = 15 customers/hour', True),
        ('+(0.0121-1)=0', None),                                              # a sum begun above
        ('\\frac{24(4.7)^2}{12}\n-\n\\frac{24(3.0)^2}{12}=0', None),           # ...one term a line
        ('Month 1: 7 × 105 = 735 (Shortage 757-735 = 22', None),              # unreadable: no join
        # THE SECOND READING (D-159): forms from round 2's notes.
        ('25 & 0.64(5)=3.20 & 3-3.20=-0.20 & -3-3.20=-6.20', True),            # table cells
        ('| n = 9 | 0.64(3) = 1.92 | 3-1.92 = 1.08 |', True),                    # markdown cells
        ('the lot is accepted if X = 0 or X = 1.', None),                         # or
        ('could not reach Cpk = 1.67 since Cp = 1.14 < 1.67', None),              # since, <
        ('I_max = 803.16 units using Q = 1090 units', None),                      # using
        ('the smallest n satisfying ARL1 ≤ 2.0 is n = 25', None),                 # is, ≤
        ('P_detect≤2 = 0.3086', None),                                            # a threshold
        (r'W_L \ge 824/179 = 4.61', True, False),                                 # ...not a threshold
        ('ẋ(0) = 0.0370 × 21.851 × cos(-0.571) = 0.68 m/s', True),               # ẋ is a variable
        ('F = 0.05243/73 = 0.718 mN', True),                                      # implied N -> mN
        ('δ = (172.91 × 10³)(1.13)/(392.7 × 10⁶) = +0.498 mm', True),            # implied m -> mm
        ('W2 = 0.06187 + 0.09091 = 0.15278 hours × 60 ≈ **9.17 minutes**', True),
        ('d = (338.86 × 10⁶)/(122 × 10⁶) = 2.7775 × 10⁻³ m', False),             # a base unit implies nothing
        ('remaining = 131 - 79 = 0 s', False),                                    # no conversion makes a zero
        ('1 m²/s = 10^4 St', True),                                               # Stokes
        ('P1 = 47.73 bar = 4.773×10⁶ Pa', True),                                  # bar
        ('A = 1 ha = 10⁴ m²', True),                                              # hectares, for this pair only
        ('x = 4.036 × 10⁻⁵ m³/mol = 0.4036 cm³/mol', False),                     # ...not for m³ and cm³
        (r'8 \, \mu \, L = 8 \times 5.00 \times 14.2 = 568', True),              # \mu L: variables
        (r'2 \mu = 2 \times 0.7 = 1.4', True),
        ('E = σ/ε = 14.043/0.00064561 = 21.75 × 10³ ksi', True),              # σ/ε links; ksi a unit
        ('p = 2 × 36.5 = 73.0 psi', True),                                    # psi is the unit
        ('p = 2 × 36.5 = 74.0 psi', False),                                   # ...and still checked
        (r'E = \frac{814}{840} \times 100 = 0.969047 \dots \times 100 = 96.9047 \dots \%', True),
        (r's = \sqrt{0.0004279296} = 0.0206865…', True),                          # rounded, then …
        ('f = 267π/2π = 133.5 Hz', True),                                         # a/2π is a/(2π)
        ('phase: 99° = 11π/20 ≈ 1.7279 rad', True),                              # 99° is a unit
        ('+(9.47839)V_m\n-(9.47839)(0.065898)=0', None),                          # a term a line
        ('Ia = 0.2 S = 0.2 × 3.5135 = 0.7027 in', True),                          # S is a variable
        ('T(1+1) = 2³ = 8 but T(1)+T(1) = 1³+1³ = 2', True),                     # but
        ('w0 = π/2 ≈ 1.5708 rad versus π/3 ≈ 1.0472 rad.', True),                # versus
        # Real slips the expert confirmed; the rules above must go on flagging them.
        (r'\sin(4\pi \times 5.9856) = \sin(75.35 \text{ rad}) \approx -0.047', False),
        (r'E_G = \frac{1.9714285714}{0.1565714286} \approx 12.5945', True, False),
        (r'Fr_1 = \frac{9.0839}{1.7439} = 5.208', True, False),
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
