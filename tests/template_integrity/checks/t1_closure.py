"""T1 - Printed-arithmetic closure.

For every emitted line of the form `... = <numeric expression> = <result>`,
evaluate the expression **from the printed operands** and confirm it matches the
printed result to within half a unit in the result's last displayed digit.

This is the check that catches a trace which does not reproduce its own answer:
`template_rotating_unbalance` Step 1 prints `omega = 1450 * (2 * pi / 60) =
152.15` where the printed operands give 151.84.

Design notes
------------
* Purely symbolic segments (`sigma = P / A`) are SKIPPED, not failed. A step
  that states a formula without substituting numbers is legitimate.
* Unit tokens are stripped before evaluation. An alphabetic token carries its
  own exponent with it, so `mm^2` is removed whole while `(10.6)^4` survives.
* Evaluation is an AST whitelist - no `eval` of template text.
* Coverage is reported. A template whose `=` lines are mostly unparseable is
  itself a finding, so the runner surfaces coverage alongside failures.
"""
from __future__ import annotations

import ast
import math
import re
from dataclasses import dataclass

from ..core import Instance, display_tolerance, parse_number

# ---------------------------------------------------------------- normalising

# Unicode and typographic operators used across the corpus.
_OPERATOR_MAP = {
    '×': '*',   # multiplication sign
    '·': '*',   # middle dot
    '∙': '*',   # bullet operator
    '÷': '/',   # division sign
    '−': '-',   # minus sign
    '–': '-',   # en dash used as minus
    '⁄': '/',   # fraction slash
}

# Names that are mathematics, not units.
_CONSTANTS = {'pi': math.pi, 'e': math.e}
_FUNCTIONS = {
    'sqrt': math.sqrt, 'ln': math.log, 'log': math.log10, 'log10': math.log10,
    'log2': math.log2, 'exp': math.exp, 'abs': abs,
    'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
    'asin': math.asin, 'acos': math.acos, 'atan': math.atan,
    'sinh': math.sinh, 'cosh': math.cosh, 'tanh': math.tanh,
    'C': math.comb, 'comb': math.comb, 'factorial': math.factorial,
}

# An alphabetic token, optionally carrying a unit exponent (mm^2, m^3, s^-1).
_UNIT_TOKEN = re.compile(
    r'(?<![\w.])'
    r'([A-Za-z°Ωμ%][A-Za-z0-9_°Ωμ%]*)'
    r'(\s*\^\s*-?\d+)?'
)
_NUMBER = re.compile(r'\d')

# A comma acting as a thousands separator: between digits, followed by exactly
# three digits. `1,541,724` -> `1541724`, but `C(20, 0)` is left alone.
_THOUSANDS = re.compile(r'(?<=\d),(?=\d{3}(?!\d))')


_TRIG = re.compile(r'\b(?:sin|cos|tan|asin|acos|atan)\s*\(')
_DEGREES = re.compile(r'(\d+(?:\.\d+)?)\s*(?:degrees|deg|°)\b')
_DEG_TO_RAD = math.pi / 180.0

# Poison marker: an unknown function call makes a segment unevaluable. Inserting
# a character that cannot parse is deliberate - it forces safe_eval to decline
# rather than silently dropping the function and evaluating its argument.
_POISON = '\x00'

# `cos(-118.6)` -> `cos((-118.6)*DEG)` for the degree reading of a bare
# trig argument (see safe_eval_variants).
_BARE_TRIG_ARG = re.compile(
    r'\b(sin|cos|tan|asin|acos|atan)\s*\(\s*'
    r'([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*\)')


def _strip_units(text: str) -> str:
    """Remove unit tokens, keeping mathematical names and all numerals.

    An unknown alphabetic token immediately followed by '(' is a function we
    cannot evaluate (`Phi(1.08)`, `Q(8.07)`, `u(t)`). Dropping just the name
    would leave `(1.08)` and evaluate the *argument* as if it were the result -
    a false failure. Such segments are poisoned so they are skipped instead.
    """
    def repl(m: re.Match) -> str:
        name = m.group(1)
        if name in _CONSTANTS or name in _FUNCTIONS:
            return m.group(0)
        after = text[m.end():m.end() + 1]
        if after == '(':
            return _POISON
        # A UNIT trails a value (`166093 N`, `9503.32 mm^2`). A SYMBOL sits in
        # an operand position (`0.009 * (LL - 10)`, `1000/CN - 10`). Stripping
        # a symbol leaves an expression that still parses but means something
        # different, so symbols poison the segment instead.
        before = text[:m.start()].rstrip()
        if not before or not (before[-1].isdigit() or before[-1] in ')]%'):
            return _POISON
        return ' '
    return _UNIT_TOKEN.sub(repl, text)


def _normalise(text: str) -> str:
    for k, v in _OPERATOR_MAP.items():
        text = text.replace(k, v)
    # `tan(43 deg)` must not become `tan(43)` - that is tangent of 43 radians.
    # Only rewrite where a trig function is actually present; elsewhere `deg`
    # is an ordinary unit and is stripped below.
    if _TRIG.search(text):
        text = _DEGREES.sub(lambda m: f'({m.group(1)}*{_DEG_TO_RAD!r})', text)
    text = _strip_units(text)
    text = text.replace('^', '**')
    # Strip ONLY thousands separators. A blanket comma strip would destroy
    # function arguments: `C(20, 0)` must not become `C(20 0)`.
    text = _THOUSANDS.sub('', text)
    text = re.sub(r'\s+', ' ', text).strip()
    # a trailing/leading operator left behind by unit stripping makes the
    # segment unparseable; that is handled by the parse failing, not patched.
    return text


# ---------------------------------------------------------------- evaluation

class _Unsafe(Exception):
    pass


def _eval_node(node: ast.AST) -> float:
    if isinstance(node, ast.Expression):
        return _eval_node(node.body)
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            raise _Unsafe('non-numeric constant')
        return float(node.value)
    if isinstance(node, ast.UnaryOp):
        v = _eval_node(node.operand)
        if isinstance(node.op, ast.USub):
            return -v
        if isinstance(node.op, ast.UAdd):
            return v
        raise _Unsafe('unary op')
    if isinstance(node, ast.BinOp):
        left, right = _eval_node(node.left), _eval_node(node.right)
        op = node.op
        if isinstance(op, ast.Add):
            return left + right
        if isinstance(op, ast.Sub):
            return left - right
        if isinstance(op, ast.Mult):
            return left * right
        if isinstance(op, ast.Div):
            if right == 0:
                raise _Unsafe('division by zero')
            return left / right
        if isinstance(op, ast.Pow):
            if abs(right) > 64 or (left < 0 and right != int(right)):
                raise _Unsafe('unsafe power')
            return left ** right
        raise _Unsafe('binop')
    if isinstance(node, ast.Name):
        if node.id in _CONSTANTS:
            return _CONSTANTS[node.id]
        raise _Unsafe(f'free name {node.id!r}')
    if isinstance(node, ast.Call):
        fname = getattr(node.func, 'id', None)
        if fname not in _FUNCTIONS or node.keywords:
            raise _Unsafe('call')
        args = [_eval_node(a) for a in node.args]
        if fname in ('C', 'comb', 'factorial'):
            args = [int(round(a)) for a in args]
        try:
            return float(_FUNCTIONS[fname](*args))
        except Exception as exc:                      # noqa: BLE001
            raise _Unsafe(f'{fname} domain') from exc
    raise _Unsafe(type(node).__name__)


def _has_arithmetic(tree: ast.AST) -> bool:
    """Does the expression actually compute something?

    A segment that reduces to a bare literal carries no arithmetic to check.
    `m = 20 subgroups of n = 5` leaves the expression `20` once prose and units
    are stripped; comparing that to the next value in the line is meaningless
    and was a large source of false failures.
    """
    return any(isinstance(n, (ast.BinOp, ast.Call)) for n in ast.walk(tree))


def safe_eval(expr: str) -> float | None:
    """Evaluate a numeric expression, or return None if it is not one."""
    if not expr or _POISON in expr or not _NUMBER.search(expr):
        return None
    try:
        tree = ast.parse(expr, mode='eval')
    except SyntaxError:
        return None
    if not _has_arithmetic(tree):
        return None
    try:
        v = _eval_node(tree)
    except _Unsafe:
        return None
    except (OverflowError, ValueError, ZeroDivisionError):
        return None
    if v != v or abs(v) == float('inf'):
        return None
    return v


def safe_eval_variants(expr: str) -> list[float]:
    """All defensible readings of a printed expression.

    Phase angles are routinely printed in degrees with no marker
    (`x1 = 43.78 * cos(-118.6) = -20.96`), and nothing in the text says which
    convention is meant. Rather than guess, T1 accepts either reading: a
    closure failure is reported only when NO reading matches. This deliberately
    trades a little sensitivity inside trig calls for freedom from false alarms.
    """
    out = []
    v = safe_eval(expr)
    if v is not None:
        out.append(v)
    if _TRIG.search(expr):
        deg = _BARE_TRIG_ARG.sub(
            lambda m: f'{m.group(1)}(({m.group(2)})*{_DEG_TO_RAD!r})', expr)
        if deg != expr:
            v2 = safe_eval(deg)
            if v2 is not None:
                out.append(v2)
    return out


def _is_unit_rescale(evaluated: float, printed: float,
                     rel_tol: float = 5e-4) -> bool:
    """Is the discrepancy exactly a power-of-ten rescale?

    A step that writes a unit conversion across an '=' (`93/2 mm = 0.0465 m`)
    looks like an arithmetic error once unit tokens are stripped. Such a
    discrepancy is always a clean power of ten, so it is separable. Reported
    as its own category rather than silently ignored: a genuine factor-of-ten
    slip would land here too, and a reviewer should see it.
    """
    if evaluated == 0 or printed == 0:
        return False
    ratio = abs(printed / evaluated)
    if not (1e-15 < ratio < 1e15):
        return False
    k = math.log10(ratio)
    if abs(k - round(k)) < rel_tol and round(k) != 0:
        return True
    # Non-decimal unit conversions written across an '=' sign:
    #   `W1 = 1 / (44 - 36) hours = 60 / 8 = 7.50 minutes`  (ratio 60)
    # A clean integer ratio is a unit conversion far more often than it is an
    # arithmetic slip, and it is reported in its own category either way.
    for r in (ratio, 1.0 / ratio):
        if 1.5 <= r <= 4000 and abs(r - round(r)) <= rel_tol * max(1.0, r):
            return True
    return False


# ---------------------------------------------------------------- the check

@dataclass
class ClosureFinding:
    template_id: str
    seed: int
    line: str
    expression: str
    evaluated: float
    printed: float
    tolerance: float
    marginal: bool = False

    @property
    def delta(self) -> float:
        return abs(self.evaluated - self.printed)

    def __str__(self) -> str:
        kind = 'MARGINAL' if self.marginal else 'FAIL'
        return (f'[{kind}] {self.template_id} seed={self.seed}\n'
                f'    line      : {self.line.strip()[:150]}\n'
                f'    evaluated : {self.expression.strip()[:90]} = {self.evaluated!r}\n'
                f'    printed   : {self.printed!r}  (tol {self.tolerance:g}, '
                f'delta {self.delta:g})')


@dataclass
class ClosureResult:
    template_id: str
    lines_with_eq: int = 0
    evaluated: int = 0          # lines from which at least one check was derived
    checks_performed: int = 0   # individual expression-vs-result comparisons
    skipped_symbolic: int = 0
    skipped_unparseable: int = 0
    unit_rescales: int = 0
    findings: list[ClosureFinding] = None
    skipped_examples: list[str] = None
    rescale_examples: list[str] = None

    def __post_init__(self):
        if self.findings is None:
            self.findings = []
        if self.skipped_examples is None:
            self.skipped_examples = []
        if self.rescale_examples is None:
            self.rescale_examples = []

    @property
    def coverage(self) -> float:
        return self.evaluated / self.lines_with_eq if self.lines_with_eq else 1.0

    @property
    def failures(self) -> list[ClosureFinding]:
        return [f for f in self.findings if not f.marginal]

    @property
    def marginals(self) -> list[ClosureFinding]:
        return [f for f in self.findings if f.marginal]

    @property
    def passed(self) -> bool:
        return not self.failures


# a printed result: a number possibly wrapped in markdown bold / parens
_RESULT_RE = re.compile(
    r'^[\s(\[*|]*([-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?)')

# Clause separators. A single emitted line routinely carries several
# independent equations:
#   "UCL_R = D4*R-bar = 2.114 * 4.72 = 9.98;  LCL_R = D3*R-bar = 0.0 * 4.72 = 0.00"
# Splitting the whole line on '=' would compare the first equation's expression
# against the last equation's result, which is meaningless.
#
# ':' is a separator too, because a prose LABEL routinely prefixes the
# arithmetic: "Row sum: 0.4347 + 0.2952 + 0.2701 = 1.0000". Without this the
# label's words are read as operand symbols, the segment is poisoned, and the
# line is filed as `skipped_symbolic` - indistinguishable from a legitimate
# formula-only step. That blind spot let a planted gross non-closure through
# the whole suite (adversary review, defect D11).
_CLAUSE_SEP = re.compile(r'[;:]')

# An operator immediately after a leading number means the segment is an
# intermediate expression, not the chain's terminal value.
_TRAILING_OP = re.compile(r'^\s*[-+*/^×·∙÷]')


def _terminal_value(segments: list[str], start: int,
                    lookahead: int = 3) -> tuple[str, float] | None:
    """Find the terminal printed value of an equation chain.

    Scans forward from `start` for the first segment whose leading number is
    NOT followed by an arithmetic operator. `3 - 1.92` is an intermediate
    expression; `1.08;` and `40.9 kN, acting...` are terminal values.
    """
    for seg in segments[start:start + lookahead]:
        m = _RESULT_RE.match(seg)
        if not m:
            return None                      # prose or a symbol: chain broken
        rest = seg[m.end():]
        if _TRAILING_OP.match(rest):
            continue                         # intermediate expression, keep going
        val = parse_number(m.group(1))
        if val is None:
            return None
        return m.group(1), val
    return None


def check_line(line: str) -> tuple[list[tuple[str, float, str, float]], str]:
    """Return (checks, skip_reason).

    Each `=` is checked independently as (expression immediately before it,
    terminal value of the chain it belongs to). checks is a list of
    (expression, evaluated, printed_token, printed_value).
    """
    if '=' not in line:
        return [], 'no-eq'
    if re.search(r'[<>!]=|=[<>]|==', line):
        return [], 'comparison'

    segments = line.split('=')
    if len(segments) < 2:
        return [], 'single-eq'

    checks: list[tuple[str, float, str, float]] = []
    reason = 'symbolic'
    for i in range(len(segments) - 1):
        # The expression is the text since the last clause separator, so an
        # expression never reaches back across a ';' into a previous equation.
        raw = _CLAUSE_SEP.split(segments[i])[-1]
        expr = _normalise(raw)
        if not expr:
            continue
        vals = safe_eval_variants(expr)
        if not vals:
            reason = 'symbolic-or-unparseable'
            continue
        term = _terminal_value(segments, i + 1)
        if term is None:
            reason = 'no-terminal-value'
            continue
        tok, printed = term
        # Pick the reading closest to what was printed; a mismatch is reported
        # only if every defensible reading misses.
        val = min(vals, key=lambda v: abs(v - printed))
        checks.append((expr, val, tok, printed))
    return checks, reason


def check_instance(inst: Instance, result: ClosureResult,
                   marginal_band: float = 0.4) -> None:
    """Accumulate closure results for one generated instance."""
    for raw in inst.solution.splitlines():
        line = raw.strip()
        if '=' not in line:
            continue
        result.lines_with_eq += 1
        checks, reason = check_line(line)
        if not checks:
            if reason in ('symbolic', 'symbolic-or-unparseable', 'single-eq',
                          'comparison', 'no-numeric-result'):
                result.skipped_symbolic += 1
            else:
                result.skipped_unparseable += 1
            if len(result.skipped_examples) < 12:
                result.skipped_examples.append(f'[{reason}] {line[:120]}')
            continue
        result.evaluated += 1
        result.checks_performed += len(checks)
        for expr, val, tok, printed in checks:
            tol = max(display_tolerance(tok), abs(printed) * 1e-9, 1e-12)
            delta = abs(val - printed)
            if delta <= tol * marginal_band:
                continue
            if _is_unit_rescale(val, printed):
                # `c1 = 93/2 mm = 0.0465 m` is a unit conversion written across
                # an '=' sign. Stripping units erases the conversion, so numeric
                # closure cannot judge it. Counted and reported, never failed.
                result.unit_rescales += 1
                if len(result.rescale_examples) < 8:
                    result.rescale_examples.append(
                        f'{expr[:60]} = {val:g} vs printed {printed:g} '
                        f'(x1e{round(math.log10(abs(printed / val))):+d})')
                continue
            result.findings.append(ClosureFinding(
                inst.template_id, inst.seed, line, expr, val, printed, tol,
                marginal=(delta <= tol)))


def run(instances: list[Instance], template_id: str) -> ClosureResult:
    res = ClosureResult(template_id=template_id)
    for inst in instances:
        if inst.ok:
            check_instance(inst, res)
    return res
