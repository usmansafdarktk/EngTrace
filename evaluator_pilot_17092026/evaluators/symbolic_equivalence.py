"""The symbolic equivalence check (WS-E, E3): does a response's final answer state the gold's expression, to the
precision the gold shows?

Nine templates answer with an expression, and the answer check scores them by the numbers they state (answer.targets
reads the numbers the gold computed). An answer in another but equal form is then scored wrong: y^2/2 for 0.5y^2, a
phase shifted by a whole turn, exact roots where the gold rounds them, erfc for Q, the value of a one-value
expression, or the gold's own function without the number the check targets. This module decides equivalence
instead, per template, from a spec that names the free variable, where to sample it, what the question gives, and
the name of the asked quantity.

READING. The gold's answer line and every expression on the response's final-answer segment (answer.segment, the
reader the score uses) are parsed into SymPy: LaTeX (\\frac, \\sqrt, \\left, cases, \\operatorname{erfc}, ^\\circ,
u[n], \\delta[n]), plain text and unicode (sqrt, superscripts, fractions, the degree and percent signs); units and
words inside an expression are dropped. Names the question gives a value to (Eb/N0, Es/N0, M, k) take it, or the
answer's own Eb/N0 where it states a different one. A number the question gives is exact in the answer (23 dB
written 10^2.3 is not rounded), and so is every integer (1 + sqrt 2), but for an angle in whole degrees and the
mantissa of a power of ten. A sequence that the answer says is zero before n = 0 is read with u[n].

EQUIVALENCE, the definition the experts graded against (guide.md: "the same quantity as the reference, to the
precision the reference shows"; "a value with more digits that rounds to the reference's"). A candidate is
equivalent when, at every sample point of the free variable (twenty over two periods of a signal; n = -3 to 11 for
a sequence; a grid for a velocity field; the support for the autocorrelation; one point for the bit error rate), it
lies within half a unit of the last digit of every decimal the gold writes, carried through the gold's expression
(GOLD_UNIT; the answer check's `unit` scales it). An answer coarser than the reference (17.4 for 17.38) or off in a
digit the reference shows is not equivalent; an exact form is. A difference SymPy simplifies to zero is equivalent at
once. For the one template whose reference is one value written as an expression (c Q(a)), whose digits the
reference does not print, a single written number is also equivalent when the reference's value rounds to it.

THE ANSWER AS A WHOLE. An answer is equivalent when an expression on its segment is, unless it also states the asked
quantity (h[n] = ..., P_b = ... = ...) as something the answer check's own number rule (0.2% or a unit of a last
digit) calls different: two values are not the reference. A bare 0 or 1 contradicts nothing (D-138), and for the
bit error rate only a written value of at most 0.5 counts. Where a template states a support (the autocorrelation's
|tau| <= 2T), a support the answer states must be the gold's.

SYMBOLIC_CHECK.md (validate.py) reports the check against the experts' grades of 302 answers and says which templates
it is applied to; answer.py applies it only to the templates named in SYMBOLIC_EQUIVALENCE_TEMPLATES, and only to
raise a verdict to correct.
"""
from __future__ import annotations

import functools
import math
import random
import re
from dataclasses import dataclass, field

import sympy as sp
from sympy.parsing.sympy_parser import (convert_xor, function_exponentiation, implicit_multiplication_application,
                                        parse_expr, standard_transformations)

REL = 0.002
SLACK = 1e-9
T, N, X, Y, TAU = sp.symbols('t n x y tau', real=True)
EBN0, ESN0, MM, KK = sp.symbols('EBN0 ESN0 MM KK', positive=True)
TRANSFORMS = standard_transformations + (implicit_multiplication_application, convert_xor, function_exponentiation)


# ------------------------------------------------------------------ text to SymPy

def _q(z):
    return sp.erfc(z / sp.sqrt(2)) / 2


def _step(z):
    return sp.Piecewise((1, sp.Ge(z, 0)), (0, True))


def _kdelta(z):
    return sp.Piecewise((1, sp.Eq(z, 0)), (0, True))


def _sinc(z):
    return sp.sin(sp.pi * z) / (sp.pi * z)


FUNCS = {'Q': sp.Lambda(sp.Symbol('_z'), _q(sp.Symbol('_z'))), 'erfc': sp.erfc, 'erf': sp.erf,
         'sqrt': sp.sqrt, 'cos': sp.cos, 'sin': sp.sin, 'tan': sp.tan, 'exp': sp.exp, 'ln': sp.log,
         'log': sp.log, 'log2': lambda z: sp.log(z, 2), 'log10': lambda z: sp.log(z, 10), 'Abs': sp.Abs,
         'ustep': sp.Lambda(sp.Symbol('_z'), _step(sp.Symbol('_z'))),
         'kdelta': sp.Lambda(sp.Symbol('_z'), _kdelta(sp.Symbol('_z'))),
         'sinc': sp.Lambda(sp.Symbol('_z'), _sinc(sp.Symbol('_z'))), 'pi': sp.pi, 'E': sp.E,
         'Piecewise': sp.Piecewise, 'Eq': sp.Eq, 'True': sp.true}
SYMS = {'t': T, 'n': N, 'x': X, 'y': Y, 'tau': TAU, 'EBN0': EBN0, 'ESN0': ESN0, 'MM': MM, 'KK': KK}

UNICODE = [('\u2212', '-'), ('\u2013', '-'), ('\u2014', '-'), ('\u00b7', '*'), ('\u22c5', '*'), ('\u00d7', '*'),
           ('\u2219', '*'), ('\u03c0', ' pi '), ('\u03c4', ' tau '), ('\u03b4', ' delta '), ('\u03c9', ' omega '),
           ('\u2248', '='), ('\u2243', '='), ('\u2245', '='), ('\u2264', '<='), ('\u2265', '>='), ('\u2223', '|'),
           ('\u00bd', '(1/2)'), ('\u2153', '(1/3)'), ('\u2154', '(2/3)'), ('\u00bc', '(1/4)'), ('\u00be', '(3/4)'),
           ('\u215b', '(1/8)'), ('\u00b0', ' deg'), ('\u2009', ' '), ('\u202f', ' '), ('\u00a0', ' '),
           ('\u03b3', ' gamma '), ('\u2080', '_0'), ('\u2081', '_1'), ('\u2082', '_2')]
FUNCTION_WORDS = {'erfc', 'erf', 'sinc', 'Q', 'sqrt', 'exp', 'log', 'ln', 'log2', 'log10', 'sin', 'cos', 'tan', 'pi',
                  'tau', 'Abs', 'ustep', 'kdelta', 'delta', 'EBN0', 'ESN0', 'MM', 'KK', 'deg'}
SUPER = str.maketrans('\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b\u207a\u207f',
                      '0123456789-+n')
SUPER_RUN = re.compile('[\u2070\u00b9\u00b2\u00b3\u2074-\u2079\u207b\u207a\u207f]+')


def _braced(s: str, i: int) -> tuple[str, int]:
    """The group starting at s[i] ('{...}' with nesting, or one character), and the index after it."""
    while i < len(s) and s[i] == ' ':
        i += 1
    if i >= len(s):
        return '', i
    if s[i] != '{':
        return s[i], i + 1
    depth = 0
    for j in range(i, len(s)):
        if s[j] == '{':
            depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
    return s[i + 1:], len(s)


def _commands(s: str) -> str:
    """\\frac{a}{b}, \\sqrt[n]{a}, \\sqrt{a} and the wrappers that only hold content, resolved innermost first."""
    for _ in range(200):
        m = re.search(r'\\[dt]?frac(?![a-z])|\\sqrt(?![a-z])|\\(?:boxed|mathrm|mathbf|mathit|operatorname|text|textrm|mbox)(?![a-z])', s)
        if not m:
            return s
        cmd, i = m.group(0), m.end()
        if 'frac' in cmd:
            a, i = _braced(s, i)
            b, i = _braced(s, i)
            rep = f' (({_commands(a)})/({_commands(b)})) '
        elif cmd == '\\sqrt':
            idx = None
            if i < len(s) and s[i] == '[':
                j = s.index(']', i)
                idx, i = s[i + 1:j], j + 1
            a, i = _braced(s, i)
            rep = f' sqrt({_commands(a)}) ' if idx is None else f' (({_commands(a)})**(1/({idx}))) '
        else:
            a, i = _braced(s, i)
            word = a.strip()
            if cmd in ('\\text', '\\textrm', '\\mbox'):
                if word in FUNCTION_WORDS:
                    rep = f' {word} '
                elif re.fullmatch(r'(?i)(?:for|if|when|where|with|and|otherwise|else)\b.*', word):
                    rep = ' ; '                  # "\text{ for } |tau| <= 8": a condition follows, not a factor
                else:
                    rep = ' '                    # a unit or a word, not part of the value
            else:
                rep = f' {_commands(a)} '
        s = s[:m.start()] + rep + s[i:]
    return s


def normalise(text: str) -> str:
    """LaTeX, unicode and plain math to the text SymPy's parser reads; conditions and words are kept for cutting."""
    s = text
    for a, b in UNICODE:
        s = s.replace(a, b)
    s = SUPER_RUN.sub(lambda m: '^(' + m.group(0).translate(SUPER) + ')', s)
    s = re.sub(r'\u221a\s*\(', 'sqrt(', s)
    s = re.sub(r'\u221a\s*([0-9.]+|[A-Za-z_]+\w*)', r'sqrt(\1)', s)
    s = re.sub(r'\\(?:left|right|bigl|bigr|Bigl|Bigr|biggl|biggr|big|Big|bigg|Bigg|displaystyle|textstyle|limits)(?![a-z])', ' ', s)
    s = re.sub(r'\\(?:,|;|:|!|quad|qquad)', ' ', s)
    s = s.replace('\\{', '(').replace('\\}', ')').replace('\\|', '|')
    s = re.sub(r'\^\s*\{\s*\\circ\s*\}|\^\s*\\circ|\\circ|\\degree', ' deg', s)
    s = re.sub(r'\{\s*\\rm\s*([^{}]*)\}', r' \1 ', s)                     # {\rm rad}
    s = re.sub(r'\\operatorname\s*\{\s*(\w+)\s*\}', r' \1 ', s)
    s = _commands(s)
    s = re.sub(r'\\log_\s*\{?\s*2\s*\}?', ' log2 ', s)
    s = re.sub(r'\\log_\s*\{?\s*10\s*\}?', ' log10 ', s)
    s = re.sub(r'\blog_\s*\(?\s*2\s*\)?\s*(?=[\w(\\])', ' log2 ', s)
    for cmd, rep in (('\\cdot', '*'), ('\\times', '*'), ('\\ast', '*'), ('\\pi', ' pi '), ('\\tau', ' tau '),
                     ('\\delta', ' delta '), ('\\omega', ' omega '), ('\\gamma', ' gamma '), ('\\sin', ' sin '),
                     ('\\cos', ' cos '), ('\\exp', ' exp '), ('\\ln', ' ln '), ('\\log', ' log '),
                     ('\\approx', '='), ('\\simeq', '='), ('\\sim', '='), ('\\leq', '<='), ('\\le', '<='),
                     ('\\geq', '>='), ('\\ge', '>='), ('\\lt', '<'), ('\\gt', '>'), ('\\neq', '!='),
                     ('\\infty', ' oo '), ('\\Lambda', ' Lambda '), ('\\Pi', ' Pi '), ('\\mid', '|'), ('\\vert', '|')):
        s = re.sub(re.escape(cmd) + r'(?![a-zA-Z])', rep, s)
    s = s.replace('{', '(').replace('}', ')')
    s = s.replace('\\', ' ')
    s = re.sub(r'(?i)\b(?:approx(?:imately)?|approx\.)\s', ' = ', s)
    s = re.sub(r'\(\s*(?:=\s*)?\d+(?:\.\d+)?\s*%\s*\)', ' ', s)          # "0.0512 (5.12%)": a restatement, not a factor
    s = re.sub(r'(\d+(?:\.\d+)?)\s*%', r'(\1/100)', s)                    # 1.25% is 0.0125
    # a unit after a number or pi inside the expression, (t - 0.00412 s) or (... - 3.1187 rad), is not part of it
    s = re.sub(r'(?:(?<=[\d.)])|(?<=pi))\s*(?:rad(?:ians)?|sec|ms|s|V|volts?|Hz)(?![\w(])', ' ', s)
    return s


# A decimal times a power of ten (3.17 x 10^-2) needs no rewriting: the decimal's unit moves with the power.
DECIMAL = re.compile(r'(?<![\w.])(\d+\.\d*|\.\d+)(?:[eE]([-+]?\d+))?(?![\w.])')


def _abs_bars(s: str) -> str:
    """|a| as Abs(a): bars paired left to right, an opening bar where the left neighbour is an operator or a start."""
    out, stack = [], []
    for i, ch in enumerate(s):
        if ch != '|':
            out.append(ch)
            continue
        prev = ''.join(out).rstrip()[-1:] if ''.join(out).strip() else ''
        if stack and (not prev or prev not in '+-*/(,=<>'):
            stack.pop()
            out.append(')')
        else:
            stack.append(i)
            out.append('Abs(')
    return ''.join(out) if not stack else s


@dataclass
class Parsed:
    expr: sp.Expr
    lits: dict = field(default_factory=dict)      # symbol -> (value, ulp)
    stmt: bool = False                            # a statement of the asked quantity: h[n] = ..., P_b = ... = ...
    full: bool = True                             # read whole, not cut back from its end


def _literals(s: str) -> tuple[str, dict]:
    """Every decimal number (it has a point) as a symbol carrying its value and one unit of its last digit; integers
    stay exact. 3.17*10^(-2) keeps 3.17 as the decimal; 1e5 written as an integer power is exact."""
    lits = {}

    def sub(m):
        mant, exp = m.group(1), m.group(2)
        frac = mant.split('.')[1] if '.' in mant else ''
        ulp = 10.0 ** (-len(frac)) if frac else 1.0
        v = float(mant)
        if exp:
            v *= 10.0 ** int(exp)
            ulp *= 10.0 ** int(exp)
        name = f'LIT{len(lits)}'
        lits[sp.Symbol(name)] = (v, ulp)
        return f' {name} '
    return DECIMAL.sub(sub, s), lits


@functools.lru_cache(maxsize=20000)
def _parse_cached(text: str, names: tuple) -> tuple | None:
    local = dict(FUNCS)
    local.update({k: v for k, v in SYMS.items() if k in names})
    s, lits = _literals(text)
    local.update({str(k): k for k in lits})
    try:
        e = parse_expr(s, local_dict=local, transformations=TRANSFORMS, evaluate=True)
    except Exception:  # noqa: BLE001  (SyntaxError, TokenError, TypeError, ...: not an expression)
        return None
    if not isinstance(e, sp.Expr) or e.has(sp.Lambda):      # Q or sinc written without its argument: a fragment
        return None
    allowed = {SYMS[k] for k in names} | set(lits)
    if not e.free_symbols <= allowed:
        return None
    return e, tuple(lits.items())


def parse(text: str, names: tuple) -> Parsed | None:
    """One expression over the named variables, or None when the text is not one."""
    s = ' '.join(text.split()).rstrip('.,;:')            # a newline outside brackets would end the expression
    if not s or len(s) > 400:
        return None
    s = _abs_bars(s)
    # A stated value's last digit vouches for it, an integer's too (answer.match): an angle in whole degrees and the
    # integer mantissa of a power of ten are written with a point, so they carry one unit like any decimal. Every
    # other integer is part of the expression's form (the 2 in 2 + sqrt 5) and stays exact.
    s = re.sub(r'(?<![\w.])(\d+)\s*\*\s*10\s*(?:\*\*|\^)\s*(?:\(\s*([-+]?\d+)\s*\)|([-+]?\d+))(?![\w.(])',
               lambda m: f'{m.group(1)}.e{m.group(2) or m.group(3)}', s)
    s = re.sub(r'(?<![\w.])(\d+)(\.\d+)?\s*deg\b', lambda m: f'({m.group(1)}{m.group(2) or "."}*pi/180)', s)
    s = re.sub(r'\bdeg\b', '*pi/180', s)
    s = re.sub(r'\bu\s*\[([^\[\]]+)\]', r'ustep(\1)', s)
    s = re.sub(r'\b(?:delta|kdelta)\s*\[([^\[\]]+)\]', r'kdelta(\1)', s)
    s = re.sub(r'\bu\s*\(\s*n\s*([-+]\s*\d+)?\s*\)', lambda m: f'ustep(n{(m.group(1) or "").replace(" ", "")})', s)
    s = re.sub(r'\bdelta\s*\(', 'kdelta(', s)
    s = s.replace('[', '(').replace(']', ')')
    if s.count('(') != s.count(')'):
        return None
    got = _parse_cached(s, names)
    if got is None:
        return None
    e, lits = got
    return Parsed(e, dict(lits))


# ------------------------------------------------------------------ candidates on a segment

STOP = re.compile(r'(?i)\s(?:for|if|when|where|with|otherwise|and|or|i\.e\.|e\.g\.|dimensionless|unitless|units?|'
                  r'volts?|meters?|metres?|seconds?|signal|bit|errors?|which|that|is|are|the|at|in|on|amplitude|'
                  r'approximately|about|essentially|per|rad|radians|sample)\b')


def _math_spans(seg: str) -> list[str]:
    """The math spans of a segment, its plain lines, and the whole segment with the math delimiters dropped (an
    answer written "7.204 $\\cos(118t + 41.27^\\circ)$" is one expression across the two)."""
    spans = re.findall(r'\\\[(.+?)\\\]|\\\((.+?)\\\)|\$\$(.+?)\$\$|\$(.+?)\$', seg, flags=re.S)
    out = [next(x for x in g if x) for g in spans]
    plain = re.sub(r'\\\[.+?\\\]|\\\(.+?\\\)|\$\$.+?\$\$|\$.+?\$', ' ', seg, flags=re.S)
    out += [ln for ln in plain.splitlines() if ln.strip()]
    # each display on one line (a display broken over lines is one expression, not fragments), the lines kept apart
    # (two displays are two statements, not a product), the inline delimiters dropped
    flat = re.sub(r'\\\[(.+?)\\\]|\$\$(.+?)\$\$', lambda m: ' ' + ' '.join((m.group(1) or m.group(2)).split()) + ' ',
                  seg, flags=re.S)
    flat = re.sub(r'\\\(|\\\)|\$', ' ', flat)
    out += [ln for ln in flat.splitlines() if ln.strip()]
    return out


def _cases(span: str) -> list[str]:
    m = re.search(r'\\begin\s*\{cases\}(.+?)\\end\s*\{cases\}', span, flags=re.S)
    if not m:
        return []
    rows = re.split(r'\\\\(?:\[[^\]]*\])?', m.group(1))
    return [r.split('&')[0] for r in rows if r.strip()]


def _piecewise(span: str, transform) -> str | None:
    """A cases block as one SymPy Piecewise text: each row's expression under its condition (n = 0, n >= 1,
    |tau| <= 6, otherwise), in order, then 0. h[n] written as -2 at n = 0 and -5 4^(n-1) from n = 1 is one function."""
    m = re.search(r'\\begin\s*\{cases\}(.+?)\\end\s*\{cases\}', span, flags=re.S)
    if not m:
        return None
    pieces = []
    for row in re.split(r'\\\\(?:\[[^\]]*\])?', m.group(1)):
        if '&' not in row:
            continue
        e, c = row.split('&', 1)
        e = transform(e).strip().strip(',')
        c = transform(c).strip().strip(',.;')
        c = re.sub(r'(?i)\b(?:if|for|when|otherwise)\b', ' ', c).strip()
        if not c:
            c = 'True'
        else:
            c = re.sub(r'(?<![<>!=])=(?!=)', '==', c)
            c = re.sub(r'^\s*(.+?)\s*==\s*(.+?)\s*$', r'Eq(\1, \2)', c)
        if not e:
            return None
        pieces.append(f'(({e}), {c})')
    return f"Piecewise({', '.join(pieces)}, (0, True))" if pieces else None


PROSE = re.compile(r'(?<![A-Za-z_])([A-Za-z]{3,})\s+([A-Za-z]{2,})(?![A-Za-z_(])')


def _cut_prose(part: str) -> str:
    """The part up to the first two words in a row that are not function or variable names: prose, not math."""
    for m in PROSE.finditer(part):
        if m.group(1) not in FUNCTION_WORDS and m.group(2) not in FUNCTION_WORDS:
            return part[:m.start()]
    return part


def _expressions(seg: str, names: tuple, transform, lhs: str = '', limit: int = 60) -> list[Parsed]:
    """Expressions on the segment over the named variables: every side of every =, in math spans, plain lines and
    the flattened segment, each cut back from its end until it parses. A side is a statement of the asked quantity
    when the text before its = ends with the quantity's name (`lhs`, a regex), or when it continues a chain that a
    statement began (P_b = A = B: A and B both); it is full when only words and units were cut from its end."""
    texts = []
    out, by_text = [], {}

    def keep(c: str, p: Parsed) -> None:
        if c in by_text:                       # the same text from another view of the segment: one candidate
            q = by_text[c]
            q.stmt, q.full = q.stmt or p.stmt, q.full and p.full
            return
        by_text[c] = p
        out.append(p)

    for span in _math_spans(seg):
        pw = _piecewise(span, transform)
        if pw is not None:
            p = parse(pw, names)
            if p is not None:
                head = transform(re.split(r'\\begin\s*\{cases\}', span)[0])
                p.stmt = bool(lhs) and bool(re.search(lhs + r'\s*=\s*$', head))
                keep(pw, p)
        texts += _cases(span)
        texts.append(re.sub(r'\\begin\s*\{cases\}.+?\\end\s*\{cases\}', ' ', span, flags=re.S))
    for raw in texts:
        norm = transform(raw)
        norm = re.sub(r'(?i)\*{0,2}answer\s*(?:\([a-z]\))?\s*:\*{0,2}', ' ', norm)
        norm = re.sub(r'(?<![\w)])\*\*|\*\*(?![\w(])', ' ', norm)     # Markdown bold, not a power between operands
        parts = re.split(r'(?<![<>!=])=(?!=)', norm)
        chain = False                          # the previous side was a full statement: a chain continues it
        for idx, whole in enumerate(parts):
            stmt = idx > 0 and ((bool(lhs) and bool(re.search(lhs + r'\s*$', parts[idx - 1]))) or chain)
            part = _cut_prose(re.split(r'<=|>=|<|>|&|;|\bfor\b', whole)[0])
            cut, got = part, None
            for _ in range(16):
                c = cut.strip().strip(',').strip()
                if not c:
                    break
                p = by_text.get(c) or parse(c, names)
                if p is not None:
                    removed = whole[len(whole) - len(whole.lstrip()) + len(c):] if whole.strip().startswith(c) else whole
                    full = not re.search(r'[\d=^*/+]|(?<!\w)-', removed) and \
                        not any(w in FUNCTION_WORDS for w in re.findall(r'[A-Za-z]\w*', removed))
                    p = Parsed(p.expr, p.lits, stmt, full)
                    keep(c, p)
                    got = p
                    break
                m = list(STOP.finditer(' ' + cut))
                if m:
                    cut = (' ' + cut)[:m[-1].start()]
                    continue
                k = max(cut.rfind(' '), cut.rfind(','), cut.rfind('('))
                if k <= 0:
                    break
                cut = cut[:k]
            chain = bool(stmt and got is not None and got.full)
            if len(out) >= limit:
                return out
    return out


def candidates(seg: str, names: tuple, lhs: str = '', limit: int = 60) -> list[Parsed]:
    return _expressions(seg, names, normalise, lhs, limit)


# ------------------------------------------------------------------ comparing two expressions

def _fn(p: Parsed, var: tuple):
    syms = list(var) + list(p.lits)
    return sp.lambdify(syms, p.expr, modules=['mpmath'])


def _values(p: Parsed, var: tuple, points: list[tuple], scale: float = 1.0) -> list[float] | None:
    f = _fn(p, var)
    base = [v for v, _ in p.lits.values()]
    try:
        return [float(f(*pt, *base)) * scale for pt in points]
    except Exception:  # noqa: BLE001
        return None


def _spread(p: Parsed, var: tuple, points: list[tuple], scale: float = 1.0, rel: float = 0.0,
            unit: float = 1.0) -> list[float]:
    """At each point, how far the value moves when each decimal moves within its window, summed over the decimals.
    A decimal's window is `unit` units of its last digit, or `rel` of its value where that is wider: answer.match's
    windows for one number, here for every number of an expression. The last digit vouches for a number only when one
    unit of it is smaller than the number (D-138): 0.1 m is held to `rel`, or its window would reach zero."""
    if not p.lits:
        return [0.0] * len(points)
    f = _fn(p, var)
    base = [v for v, _ in p.lits.values()]
    ulps = [max(unit * u if u < abs(v) else 0.0, rel * abs(v)) for v, u in p.lits.values()]
    out = []
    for pt in points:
        try:
            v0 = float(f(*pt, *base))
        except Exception:  # noqa: BLE001
            out.append(math.inf)
            continue
        s = 0.0
        for j in range(len(base)):
            hi, lo = list(base), list(base)
            hi[j] += ulps[j]
            lo[j] -= ulps[j]
            try:
                s += max(abs(float(f(*pt, *hi)) - v0), abs(float(f(*pt, *lo)) - v0))
            except Exception:  # noqa: BLE001
                s = math.inf
        out.append(s * scale)
    return out


GOLD_UNIT = 0.5      # the reference's precision: half a unit of each gold number's last digit (see SYMBOLIC_CHECK.md)


def compare(gold: Parsed, resp: Parsed, var: tuple, points: list[tuple], unit: float = 1.0, scales: tuple = (1.0,),
            gold_unit: float | None = None, lenient: bool = False) -> tuple[bool, dict]:
    """Is `resp` the gold to the precision the gold shows, at every point, under one of the unit scales?

    The window at a point is how far the gold moves when each of its decimals moves by `gold_unit` (times `unit`) of
    its last digit: an answer within it reproduces every number the gold writes to that number's precision, so a
    coarser rounding (17.4 for 17.38) is outside it while an exact form (1 + sqrt 2 for 2.41) is inside. `lenient`
    is the answer check's own number rule instead (each gold decimal within 0.2% or a unit of its last digit, each
    answer decimal within a unit of its own): what tells a second, different value from a rounded one."""
    gold_unit = GOLD_UNIT if gold_unit is None else gold_unit
    exact = None
    try:
        if not gold.lits and not resp.lits:
            d = sp.simplify(gold.expr - resp.expr)
            exact = d == 0
            if exact:
                return True, {'how': 'simplify'}
    except Exception:  # noqa: BLE001
        pass
    g = _values(gold, var, points)
    if g is None or any(not math.isfinite(v) for v in g):
        return False, {'how': 'gold not evaluable'}
    noise = 1e-9 * max(abs(v) for v in g) + 1e-300            # floating point only: an exact answer must be exact
    if lenient:
        gs = _spread(gold, var, points, rel=REL, unit=unit)  # each gold decimal: 0.2% of it or a unit of its last digit
    else:
        gs = [unit * gold_unit * v for v in _spread(gold, var, points)]
    worst = None
    for sc in scales:
        r = _values(resp, var, points, sc)
        if r is None or any(not math.isfinite(v) for v in r):
            continue
        diffs = [abs(rv - gv) for rv, gv in zip(r, g)]
        if all(d <= noise for d in diffs):
            return True, {'how': 'equal', 'scale': sc, 'margin': 0.0}
        rs = _spread(resp, var, points, sc, unit=unit) if lenient else [0.0] * len(points)
        allows = [max(gu, ru, noise) * (1 + SLACK) for gu, ru in zip(gs, rs)]
        margin = max(d / a for d, a in zip(diffs, allows))
        if margin <= 1.0:
            return True, {'how': 'sampled', 'scale': sc, 'margin': round(margin, 3)}
        worst = margin if worst is None else min(worst, margin)
    return False, {'how': 'outside' if worst is not None else 'not evaluable',
                   'margin': None if worst is None else round(worst, 3), 'simplify': exact}


# ------------------------------------------------------------------ the templates

@dataclass
class Spec:
    names: tuple                                 # the free variables, as SYMS keys
    gold: callable                               # item -> Parsed of the gold's answer
    points: callable                             # item -> sample points, tuples ordered as names
    params: callable = None                      # item -> {symbol: value} the question gives
    scales: tuple = (1.0,)
    support: callable = None                     # (item, segment) -> False when a stated support is not the gold's
    lhs: str = ''                                # the asked quantity's name, as normalised text: h[n], y_c(t), QTY


def _answer_segment(text: str, segment=None) -> str:
    """The final-answer segment as the answer check reads it: answer.segment, unless the caller passes its reader."""
    if segment is None:
        try:
            import answer  # noqa: PLC0415
        except ImportError:                    # run on its own: the evaluators' directory, from this file's place
            import sys  # noqa: PLC0415
            from pathlib import Path  # noqa: PLC0415
            sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'evaluator_pilot_17092026' / 'evaluators'))
            import answer  # noqa: PLC0415
        segment = answer.segment
    return segment(text)


def _gold_rhs(item: dict, pattern: str) -> str:
    seg = _answer_segment(item['solution'])
    m = re.search(pattern, seg)
    if not m:
        raise ValueError(f"{item['item_id']}: the gold's answer line does not match {pattern!r}")
    return m.group(1)


def _gold_parse(item: dict, pattern: str, names: tuple) -> Parsed:
    p = parse(normalise(_gold_rhs(item, pattern)), names)
    if p is None:
        raise ValueError(f"{item['item_id']}: the gold's expression does not parse")
    return p


def _even(lo: float, hi: float, k: int = 20, seed: str = '') -> list[tuple]:
    rng = random.Random(seed)
    return [(lo + (hi - lo) * (i + 0.5 + 0.4 * (rng.random() - 0.5)) / k,) for i in range(k)]


def _sin_period(item: dict, pattern: str) -> float:
    m = re.search(pattern, item['question'])
    return float(m.group(1))


def _cd_dc_points(item):
    w = _sin_period(item, r'cos\((\d+(?:\.\d+)?)\*pi\*t')            # x_c(t) = A cos(W*pi*t): period 2/W
    return _even(0.0, 4.0 / w, seed=item['item_id'])


def _phasor_points(item):
    w = _sin_period(item, r'\((\d+(?:\.\d+)?)\*t')
    return _even(0.0, 2 * 2 * math.pi / w, seed=item['item_id'])


def _undamped_points(item):
    m = re.search(r'mass of ([\d.]+) kg and a spring stiffness of ([\d.]+) N/m', item['question'])
    w = math.sqrt(float(m.group(2)) / float(m.group(1)))
    return _even(0.0, 2 * 2 * math.pi / w, seed=item['item_id'])


def _autocorr_width(item) -> float:
    m = re.search(r'for\s+(-?[\d.]+)\s*<=\s*t\s*<=\s*([\d.]+)', item['question'])
    return float(m.group(2)) - float(m.group(1))                       # the pulse width; the support is |tau| <= width


def _autocorr_points(item):
    w = _autocorr_width(item)
    return _even(-w, w, seed=item['item_id'])


def _autocorr_support(item, seg: str) -> bool:
    w = _autocorr_width(item)
    s = normalise(seg)
    for m in re.finditer(r'\|\s*tau\s*\|\s*(?:<=|<)\s*(\d+(?:\.\d+)?)', s):
        if abs(float(m.group(1)) - w) > 1e-9:
            return False
    return True


def _grid(item):
    rng = random.Random(item['item_id'])
    return [(rng.uniform(-2, 2), rng.uniform(-2, 2)) for _ in range(20)]


def _ber_params(item):
    m = re.search(r'uses (\d+)-(PSK|QAM)', item['question'])
    db = re.search(r'Eb/N0 = ([\d.]+) dB', item['question'])
    M = int(m.group(1))
    ebn0 = 10 ** (float(db.group(1)) / 10)
    k = math.log2(M)
    return {EBN0: ebn0, ESN0: k * ebn0, MM: M, KK: k}


def _ber_text(s: str) -> str:
    """The BER answer's text with the names the question gives a value to: Eb/N0 (and the gamma_b, Es/N0 forms), M
    and k. The forms are matched longest first, so \\frac{E_b}{N_0}, normalised to ((E_b)/(N_0)), is one name."""
    s = normalise(s)
    e = {'b': r'E\s*_?\s*(?:\(\s*b\s*\)|b)', 's': r'E\s*_?\s*(?:\(\s*s\s*\)|s)'}
    n0 = r'N\s*_?\s*(?:\(\s*0[\s\w,]*\)|0)'                 # N_0, N_{0}, N_{0,\rm lin}
    for k, name in (('b', ' EBN0 '), ('s', ' ESN0 ')):
        s = re.sub(r'\(\s*\(\s*' + e[k] + r'\s*\)\s*/\s*\(\s*' + n0 + r'\s*\)\s*\)', name, s)     # ((E_b)/(N_0))
        s = re.sub(r'\(\s*' + e[k] + r'\s*\)\s*/\s*\(\s*' + n0 + r'\s*\)', name, s)                  # (E_b)/(N_0)
        s = re.sub(r'\(\s*' + e[k] + r'\s*/\s*' + n0 + r'\s*\)', name, s)                            # (E_b/N_0)
        s = re.sub(e[k] + r'\s*/\s*' + n0, name, s)                                                  # E_b/N_0, Eb/N0
    s = re.sub(r'\bgamma_?\s*\(?\s*b\s*\)?', ' EBN0 ', s)
    s = re.sub(r'\bM\b', ' MM ', s)
    s = re.sub(r'(?<![\w.])k(?![\w(])', ' KK ', s)
    s = re.sub(r'\b(?:P_?\s*\(?\s*b\s*\)?|BER|P_?\s*\(?\s*e\s*\)?)(?![\w(])', ' QTY ', s)     # the asked quantity
    s = re.sub(r'\b(?:SER|P_?\s*\(?\s*s\s*\)?)(?![\w(])', ' SERQ ', s)                         # another quantity
    return s


FN_OF = r'\s*_?\s*(?:\(\s*\w*\s*\)|\w+)?\s*'           # an optional subscript, as normalised: y_(c), v_( ), R_g
SPECS: dict[str, Spec] = {
    'template_cd_dc_system_analysis': Spec(
        names=('t',), gold=lambda it: _gold_parse(it, r'y_c\(t\)\s*=\s*(.+?)\.?\s*$', ('t',)), points=_cd_dc_points,
        lhs=r'\by' + FN_OF + r'\(\s*t\s*\)'),
    'template_phasor_addition': Spec(
        names=('t',), gold=lambda it: _gold_parse(it, r'v_total\(t\)\s*=\s*(.+?)\.?\s*$', ('t',)), points=_phasor_points,
        lhs=r'\bv' + FN_OF + r'\(\s*t\s*\)'),
    'template_undamped_response_initial_conditions': Spec(
        names=('t',), gold=lambda it: _gold_parse(it, r'x\(t\)\s*=\s*(.+?)\s*\(m\)', ('t',)), points=_undamped_points,
        scales=(1.0, 1e-3, 1e-2, 1e3), lhs=r'\bx\s*\(\s*t\s*\)'),
    'template_impulse_response_from_lccde': Spec(           # n < 0 too: h[n] is causal, and an answer must say so
        names=('n',), gold=lambda it: _gold_parse(it, r'h\[n\]\s*=\s*(.+?)\s*$', ('n',)),
        points=lambda it: [(float(i),) for i in range(-3, 12)], lhs=r'\bh\s*[\[(]\s*n\s*[\])]'),
    'template_incompressible_continuity': Spec(
        names=('x', 'y'), gold=lambda it: _gold_parse(it, r'\b[uv]\s*=\s*(.+?)\.?\s*$', ('x', 'y')), points=_grid,
        lhs=r'\b[uv]\s*(?:\(\s*x\s*,\s*y\s*\))?'),
    'template_autocorrelation_rect_pulse': Spec(
        names=('tau',), gold=lambda it: _gold_parse(it, r'R_g\(tau\)\s*=\s*(.+?)\*\*,', ('tau',)), points=_autocorr_points,
        support=_autocorr_support, lhs=r'\bR' + FN_OF + r'\(\s*tau\s*\)'),
    'template_ber_estimation_mary': Spec(
        names=(), gold=lambda it: _gold_parse(it, r'BER approx\s*(.+?)\s*$', ()), points=lambda it: [()],
        params=_ber_params, lhs=r'\bQTY'),
}
QNUM = re.compile(r'(?<![\w.])\d+(?:\.\d+)?')
CAUSAL = re.compile(r'\bn\s*(?:>=|>)\s*0|\bn\s*<\s*0|\bcausal\b|ustep')


def _substitute(p: Parsed, params: dict | None) -> Parsed:
    if not params:
        return p
    return Parsed(p.expr.subs(params), p.lits, p.stmt, p.full)


@functools.lru_cache(maxsize=4096)
def _gold_cached(template_id: str, item_id: str, solution: str) -> Parsed:
    return SPECS[template_id].gold({'template_id': template_id, 'item_id': item_id, 'solution': solution})


def gold_of(item: dict) -> Parsed:
    return _gold_cached(item['template_id'], item['item_id'], item['solution'])


def _given(question: str) -> list[float]:
    """The numbers a question gives, and each over ten (a decibel value's exponent, 23 dB as 10^2.3): exact inputs."""
    out = []
    for m in QNUM.finditer(question):
        v = float(m.group(0))
        out += [v, v / 10.0]
    return out


def _exact_givens(p: Parsed, given: list[float]) -> Parsed:
    """A decimal the answer writes that is a number the question gives carries no rounding: it is the input."""
    lits = {s: ((v, 0.0) if any(abs(v - g) <= 1e-12 * max(1.0, abs(g)) for g in given) else (v, u))
            for s, (v, u) in p.lits.items()}
    return Parsed(p.expr, lits, p.stmt, p.full)


def _stated_ebn0(seg: str, params: dict) -> dict:
    """The Eb/N0 an answer states and uses, where it is not the question's within the digits it shows (a conversion
    slip, 61.38 for 17.86 dB): the answer is evaluated with its own value. A dB value is converted."""
    s = _ber_text(re.sub(r'\\(?:text|mathrm|rm)\s*\{?\s*dB\s*\}?', ' dB ', seg))      # keep the unit that decides
    db_given = 10 * math.log10(params[EBN0])
    for m in re.finditer(r'EBN0\s*=\s*(\d+(?:\.\d+)?)(?![\^*(\d.])(?!\s*[\^*])\s*(dB)?', s):
        lit = m.group(1)
        frac = lit.split('.')[1] if '.' in lit else ''
        ulp = 10.0 ** (-len(frac)) if frac else 1.0
        v = float(lit)
        if abs(v - db_given) <= 1e-9 * max(1.0, db_given):
            continue                           # the question's own decibel value, restated
        if m.group(2):
            if abs(v - 10 * math.log10(params[EBN0])) > ulp:
                v = 10 ** (v / 10)
                return {**params, EBN0: v, ESN0: params[KK] * v}
        elif abs(v - params[EBN0]) > ulp:
            return {**params, EBN0: v, ESN0: params[KK] * v}
    return params


def _rounds_to(gold: Parsed, p: Parsed, unit: float = 1.0) -> bool:
    """For a reference that is one value written as an expression (c Q(a)), which prints no digits of that value: is
    the answer one written number that the reference's value rounds to? 0.034 for a value of 0.03386 is; 1.42e-11 for
    1.35e-11 is not."""
    if len(p.lits) != 1 or p.expr.atoms(sp.Function):
        return False
    g, r = _values(gold, (), [()]), _values(p, (), [()])
    (v, u), = p.lits.values()
    lit = next(iter(p.lits))
    try:                                                                  # the written digit's unit, in the value
        scale = abs(float(p.expr.diff(lit).subs(lit, v))) if p.expr.has(lit) else 0.0
    except (TypeError, ValueError):
        return False
    return bool(g and r and u and scale) and abs(r[0] - g[0]) <= 0.5 * unit * u * scale * (1 + SLACK)


def _trivial(p: Parsed) -> bool:
    """A bare 0 or 1 ("BER = 0", for a value of 1e-62) vouches for no value (D-138): it does not contradict."""
    return not p.lits and p.expr.is_number and p.expr in (0, 1)


def check(item: dict, text: str, rel: float = REL, unit: float = 1.0, segment=None,
          gold_unit: float | None = None) -> tuple[bool, dict]:
    """Is the final answer on `text` (a full response) equivalent to the item's gold? (equivalent, detail).

    Every expression on the segment is a candidate; the answer is equivalent when one is, unless a statement of the
    asked quantity on the segment (h[n] = ..., P_b = ... = ...), read whole, is not: an answer that states two values
    is not the reference. `segment` is the reader of the final answer, answer.segment unless given. `rel` is accepted
    for the answer check's signature; the window is the reference's own precision (compare)."""
    spec = SPECS[item['template_id']]
    gold = gold_of(item)
    seg = _answer_segment(text, segment)
    if spec.support is not None and not spec.support(item, seg):
        return False, {'how': 'support differs'}
    given = _given(item.get('question', ''))
    if item['template_id'] == 'template_ber_estimation_mary':
        params = _stated_ebn0(seg, spec.params(item))
        cands = [_substitute(p, params) for p in _ber_candidates(seg)]
    else:
        cands = candidates(seg, spec.names, spec.lhs)
    cands = [_exact_givens(p, given) for p in cands]
    var = tuple(SYMS[k] for k in spec.names)
    if item['template_id'] == 'template_impulse_response_from_lccde' and CAUSAL.search(normalise(seg)):
        cands = [p if p.expr.has(sp.Piecewise) else Parsed(p.expr * _step(N), p.lits, p.stmt, p.full) for p in cands]
    points = spec.points(item)
    judged = []
    for p in cands:
        if p.expr.free_symbols - set(p.lits) - set(var):
            continue
        ok, d = compare(gold, p, var, points, unit, spec.scales, gold_unit)
        if not ok and spec.names == () and _rounds_to(gold, p, unit):
            ok, d = True, {'how': 'the value, rounded', 'margin': d.get('margin')}
        judged.append((p, ok, d))
    hits = [(p, d) for p, ok, d in judged if ok]
    if not hits:
        best = min((d for _p, _ok, d in judged if d.get('margin') is not None), key=lambda d: d['margin'], default=None)
        return False, {**(best or {'how': 'no candidate'}), 'candidates': len(cands)}
    # A contradiction is a whole statement of the quantity that even the answer check's own number rule calls
    # different. For the bit error rate, a stated value: a formula in M and k parsed from prose is not evidence of one.
    value_only = item['template_id'] == 'template_ber_estimation_mary'

    def a_value(p: Parsed) -> bool:
        """A stated bit error rate: a written decimal (not an exact coefficient such as 15/32) that is a probability of
        at most 0.5 ("4-5 x 10^-26", read as 4 - 5e-26, is not one)."""
        v = _values(p, var, points)
        return bool(p.lits) and not p.expr.atoms(sp.Function) and v is not None and 0 < v[0] <= 0.5

    clash = [p for p, ok, _d in judged if not ok and p.stmt and p.full and not _trivial(p)
             and (not value_only or a_value(p))
             and not compare(gold, p, var, points, unit, spec.scales, lenient=True)[0]]
    if clash:
        return False, {'how': 'conflicting statements', 'candidates': len(cands), 'expr': str(clash[0].expr)[:120]}
    p, d = hits[0]
    return True, {**d, 'candidates': len(cands), 'expr': str(p.expr)[:120]}


def apply(label: str, item: dict, text: str, enabled, tol: float | None = None, unit: float = 1.0,
          segment=None) -> tuple[str, dict | None]:
    """The answer check's hook. On a template in `enabled`, a verdict other than correct becomes correct when the final
    answer is equivalent to the gold; nothing else changes and no verdict is lowered. (label, None) where the check
    did not run; (label, detail) where it did, detail['equivalent'] saying what it found. A failure inside the check
    leaves the verdict as it was and says so."""
    if label == 'correct' or item.get('template_id') not in enabled or item.get('template_id') not in SPECS:
        return label, None
    try:
        ok, d = check(item, text, REL if tol is None else tol, unit, segment)
    except Exception as e:  # noqa: BLE001  (a parse SymPy cannot finish is not an equivalence)
        return label, {'equivalent': False, 'how': f'error {type(e).__name__}'}
    return ('correct' if ok else label), {'equivalent': ok, 'how': d.get('how'), 'margin': d.get('margin')}


def _ber_candidates(seg: str) -> list[Parsed]:
    """BER answers: each expression or value on the segment, with Eb/N0, Es/N0, M and k as names to substitute."""
    return _expressions(seg, ('EBN0', 'ESN0', 'MM', 'KK'), _ber_text, SPECS['template_ber_estimation_mary'].lhs)


# ------------------------------------------------------------------ selftest

def _item(tid: str, question: str, gold: str, n: int = 0) -> dict:
    return {'template_id': tid, 'item_id': f"{tid.replace('template_', '')}#{n}", 'question': question,
            'solution': f'**Answer:**\n{gold}'}


def selftest() -> int:
    """Every form and every defect this module was written for, on synthetic items (no pool text):

        python -m full_run_28092026.symbolic.equivalence --selftest
    """
    ph = _item('template_phasor_addition', 'v1(t) = 18.42 * cos(210*t + 12.6 deg)\nv2(t) = 13.37 * sin(210*t - 128.3 deg)',
               '   The sum of the two signals is v_total(t) = 14.40 * cos(210*t + 58.69 deg).')
    ph2 = _item('template_phasor_addition', 'v1(t) = 27.15 * sin(314*t - 96.2 deg)', '   The sum of the two signals is '
                'v_total(t) = 31.86 * cos(314*t + 138.47 deg).', 1)
    cd = _item('template_cd_dc_system_analysis', 'x_c(t) = 12.47 * cos(164*pi*t)',
               'The final continuous-time output signal is y_c(t) = 8.31 * cos(164*pi*(t - 0.010900)).')
    ir2 = _item('template_impulse_response_from_lccde', 'y[n] + 6y[n-1] - 1y[n-2] = 1x[n]', 'Substituting the '
                'coefficients back into the general form, the impulse response is:\nh[n] = (0.03*(0.16)^n + 0.97*(-6.16)^n) * u[n]')
    ir1 = _item('template_impulse_response_from_lccde', 'y[n] - 3y[n-1] = -2x[n] + 1x[n-1]',
                'h[n] = -2*delta[n] - 5*(3)^(n-1)*u[n-1]', 1)
    ir1b = _item('template_impulse_response_from_lccde', 'y[n] + 1y[n-1] = 5x[n] + 4x[n-1]',
                 'h[n] = 5*delta[n] - 1*(-1)^(n-1)*u[n-1]', 2)
    inc = _item('template_incompressible_continuity', 'u = 3x^2 + 1xy',
                'The simplest expression for the y-component of velocity is v = -6xy - 0.5y^2.')
    inc2 = _item('template_incompressible_continuity', 'v = 3y^2 + 4xy',
                 'The simplest expression for the x-component of velocity is u = -2x^2 - 6xy.', 1)
    ac = _item('template_autocorrelation_rect_pulse', 'g(t) = 2 for -1.5 <= t <= 1.5, and 0 otherwise',
               'The autocorrelation function is a **triangular function** described by:\n**R_g(tau) = 4*(3 - |tau|)**, '
               'for **|tau| <= 3**, and **0** otherwise.\nThis triangle has a peak value of **12** at tau = 0.')
    ud = _item('template_undamped_response_initial_conditions', 'a mass of 2 kg and a spring stiffness of 361.25 N/m',
               'The final equation of motion is:\nx(t) = -0.031*cos(13.4397*t) + 0.1874*sin(13.4397*t) (m)')
    ud01 = _item('template_undamped_response_initial_conditions', 'a mass of 254.37 kg and a spring stiffness of 9318 N/m',
                 'The final equation of motion is:\nx(t) = 0.1*cos(6.0524*t) (m)', 1)
    ber = _item('template_ber_estimation_mary', 'A digital communication system uses 16-QAM modulation over an AWGN '
                'channel.\nThe system operates with a signal-to-noise ratio per bit of Eb/N0 = 19.84 dB.',
                'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.750 * Q(8.781)')
    ber2 = _item('template_ber_estimation_mary', 'A digital communication system uses 16-PSK modulation over an AWGN '
                 'channel.\nThe system operates with a signal-to-noise ratio per bit of Eb/N0 = 21.5 dB.',
                 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.50 * Q(6.56)', 1)
    ber3 = _item('template_ber_estimation_mary', 'A digital communication system uses 64-QAM modulation over an AWGN '
                 'channel.\nThe system operates with a signal-to-noise ratio per bit of Eb/N0 = 9.73 dB.',
                 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.583 * Q(1.639)', 2)
    ber4 = _item('template_ber_estimation_mary', 'A digital communication system uses 8-PSK modulation over an AWGN '
                 'channel.\nThe system operates with a signal-to-noise ratio per bit of Eb/N0 = 23.0 dB.',
                 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.67 * Q(13.24)', 3)
    A = '## Final Answer\n**Answer:** '
    cases = [
        # phasors: forms, a degree read whole (138.47 is not 13 x 8.47), sin for cos; a rounding coarser than the
        # reference, or a slip beyond its digits, is not the reference
        (ph, A + 'v_total(t) = 14.40 cos(210t + 58.69°) V', True),
        (ph, A + '\\(14.40\\cos(210t + 58.69^\\circ)\\ \\text{V}\\)', True),
        (ph, A + '14.40 sin(210t + 148.69 deg)', True),
        (ph, A + '14.4019 cos(210t + 58.6912°)', True),                 # more digits, rounding to the reference's
        (ph, A + '14 cos(210t + 59°)', False),                          # coarser: a rounding, not the reference
        (ph, A + '14.34 cos(210t + 58.72 deg)', False),
        (ph2, A + '31.85 * cos(314*t + 138.52 deg)', False),
        (ph2, A + '30.912 \\cos(314t + 136.83^\\circ) \\text{ V}', False),
        # sampled signals: the phase a whole turn away, in radians, a unit inside the expression
        (cd, A + '$y_c(t) = 8.3127 \\cos(164\\pi t + 0.2124\\pi)$', True),
        (cd, A + 'y_c(t) = 8.3127 * cos(164*pi*t - 1.7876*pi)', True),
        (cd, A + 'y_c(t) = 8.31 cos(164π t − 5.61591 rad)', True),
        (cd, A + '\\(y_c(t)=8.31\\cos\\left[164\\pi\\left(t-10.9\\times10^{-3}\\,\\mathrm{s}\\right)\\right]\\)', True),
        (cd, A + 'y_c(t) = 8.3127 * cos(164*pi*t - 1.6876*pi)', False),
        (cd, A + 'y_c(t) = 8.3127 * cos(164*pi*t - 1.79*pi)', False),        # the phase coarser than the reference
        (cd, A + '8.31 unitless', False),
        # sequences: exact roots for rounded ones, a cases block, prose that states causality; a wrong sign or start,
        # a sequence that is not zero before n = 0, two different answers
        (ir2, A + '$h[n] = \\left[ \\left(\\frac{\\sqrt{10}-3}{2\\sqrt{10}}\\right)(\\sqrt{10}-3)^n + '
                  '\\left(\\frac{\\sqrt{10}+3}{2\\sqrt{10}}\\right)(-3-\\sqrt{10})^n \\right] u[n]$', True),
        (ir2, A + '\\[\nh[n]=\n\\frac{(-3+\\sqrt{10})^{n+1}-(-3-\\sqrt{10})^{n+1}}{2\\sqrt{10}}\\,u[n]\n\\]', True),
        (ir2, A + '$h[n] = \\frac{(-3+\\sqrt{10})^{n+1}-(-3-\\sqrt{10})^{n+1}}{2\\sqrt{10}}$ for $n \\geq 0$, $h[n] = 0$ '
                  'for $n < 0$', True),
        (ir2, A + '\\[ h[n]=\\frac{(3+\\sqrt{10})^n-(3-\\sqrt{10})^n}{2\\sqrt{10}}\\,u[n] \\]', False),
        (ir2, A + '$h[n]=\\big[6(-3)^{n}-3(-2)^{n}\\big]u[n]$. Remark: as written, '
                  '$h[n]=\\frac{(-3+\\sqrt{10})^{n+1}-(-3-\\sqrt{10})^{n+1}}{2\\sqrt{10}}u[n]$', False),
        (ir1, A + '$$ h[n] = \\begin{cases} -2, & n = 0 \\\\ -5 \\cdot 3^{n-1}, & n \\geq 1 \\\\ 0, & n < 0 \\end{cases} $$', True),
        (ir1b, A + 'h[n] = 5\\delta[n] - (-1)^n u[n]', False),
        (ir1b, A + '$ h[n] = (-1)^n + 4\\delta[n] $', False),             # equal from n = 0, not before it
        # velocity fields: a fraction for a decimal, a different coefficient
        (inc, A + 'v = −6xy − y²/2 (dimensionless)', True),
        (inc, A + '\\(v(x,y) = -6xy - \\frac{y^2}{2}\\) m/s', True),
        (inc2, A + 'u(x, y) = -2x^2 - 3xy', False),
        (inc2, A + '\\(-4x^{2}-6xy\\) (units of velocity)', False),
        # the autocorrelation: the gold's own function without the peak the answer check targets; another support
        (ac, A + '\\[ R_g(\\tau)= \\begin{cases} 4(3-|\\tau|), & |\\tau|\\le 3\\\\ 0, & |\\tau|>3 \\end{cases} \\]', True),
        (ac, A + 'R_g(τ) = 12(1 − |τ|/3) for |τ| ≤ 3, and 0 otherwise', True),
        (ac, A + 'R_g(\\tau) = 12 - 4|\\tau| \\text{ for } |\\tau| \\le 3, \\text{ and } 0 \\text{ otherwise}', True),
        (ac, A + 'R_g(τ) = 4(3 − |τ|) for |τ| ≤ 6', False),
        # oscillators: amplitude and phase for cos and sin, millimetres; a gold printed as 0.1 m has no wider window
        (ud, A + 'x(t) = 0.18995 cos(13.4397 t − 1.73473) m', True),
        (ud, A + 'x(t) = −31.0 cos(13.4397 t) + 187.4 sin(13.4397 t) mm', True),
        (ud, A + 'x(t) = 0.19 cos(13.4 t − 1.73) m', False),                 # coarser than the reference
        (ud01, A + 'x(t) = 0.0871 sin(27.36 t) m', False),
        (ud01, A + 'x(t) = 100 cos(6.0524 t) mm', True),
        # the bit error rate: erfc for Q, Eb/N0 left as a name, the value as a number that the reference's value rounds
        # to; a missing factor, a value that is not that rounding, a second, different value
        (ber, A + '\\[ P_b \\approx \\frac{3}{8}\\,\\operatorname{erfc}\\!\\left(\\sqrt{\\frac{2}{5}\\,\\frac{E_b}{N_0}}\\right) \\]', True),
        (ber, A + '$P_b \\approx 0.75 Q\\left(\\sqrt{0.8 \\frac{E_b}{N_0}}\\right)$', True),
        (ber, A + '\\(6.08\\times10^{-19}\\) (dimensionless)', True),
        (ber, A + '\\( 1.1 \\times 10^{-18} \\)', False),
        (ber, A + 'BER ≈ Q(8.781)', False),
        (ber, A + '\\[ P_b \\approx \\frac{3}{8}\\operatorname{erfc}\\left(\\sqrt{\\frac{2}{5}\\frac{E_b}{N_0}}\\right) \\] '
                  'Numerically, with \\(E_b/N_0 = 19.84\\) dB, \\[ P_b \\approx 3.8 \\times 10^{-14} \\]', False),
        (ber2, A + '\\[ P_b \\approx \\frac{1}{2}Q\\!\\left(\\sqrt{8\\cdot 10^{21.5/10}}\\;\\sin\\left(\\frac{\\pi}{16}\\right)\\right) \\]', True),
        (ber2, A + '\\(P_b \\approx 1.3 \\times 10^{-11}\\)', True),
        (ber2, A + '\\(P_b \\approx 1.42 \\times 10^{-11}\\)', False),
        (ber2, A + '$P_b \\approx \\frac{1}{2} Q\\left(\\sqrt{2 \\cdot 10^{2.15}} \\sin\\left(\\frac{\\pi}{16}\\right)\\right)$', False),
        (ber3, A + '\\(\\mathrm{BER}\\approx 0.0295\\) (≈ 2.95%)', True),     # a percentage restates, it does not multiply
        (ber3, A + 'For \\(E_b/N_0=9.73\\ \\text{dB}\\), \\(P_b \\approx 2.95\\times 10^{-2}\\)', True),   # dB stays dB
        (ber3, A + '\\( \\frac{7}{12} Q\\left(\\sqrt{\\frac{2}{7}\\frac{E_b}{N_0}}\\right) \\) where '
                   '\\( \\frac{E_b}{N_0} = 9.0 \\) (linear)', False),       # its own Eb/N0, not the question's
        (ber4, A + '$\\frac{2}{3} Q\\left( \\sqrt{6 \\frac{E_b}{N_0}} \\sin\\left( \\frac{\\pi}{8} \\right) \\right) \\approx 0$', True),
        (ber4, A + '$ \\frac{1}{3} Q\\left( \\sqrt{2 \\cdot 3 \\cdot 10^{2.3} \\cdot \\sin^2\\left(\\frac{\\pi}{8}\\right)} \\right) $', False),
    ]
    bad = 0
    for item, text, want in cases:
        got, d = check(item, text)
        if got != want:
            bad += 1
            print(f"FAIL {item['item_id']}: want {want}, got {got} ({d}) | {' '.join(text.split())[:90]}")
    golds = (ph, ph2, cd, ir2, ir1, ir1b, inc, inc2, ac, ud, ud01, ber, ber2, ber3, ber4)
    for item in golds:
        if not check(item, item['solution'])[0]:
            bad += 1
            print(f"FAIL {item['item_id']}: the gold is not equivalent to itself")
    ac_text = next(t for it, t, w in cases if it is ac and w)
    inc_text = next(t for it, t, w in cases if it is inc2 and not w)
    lab, d = apply('partial', ac, ac_text, enabled=('template_autocorrelation_rect_pulse',))
    bad += not (lab == 'correct' and d['equivalent'])
    lab, d = apply('partial', ac, ac_text, enabled=())
    bad += not (lab == 'partial' and d is None)
    lab, d = apply('correct', inc2, inc_text, enabled=('template_incompressible_continuity',))
    bad += not (lab == 'correct' and d is None)                     # never lowered, not even looked at
    n = len(cases) + len(golds) + 3
    print(f'{n - bad}/{n} cases pass')
    return bad


if __name__ == '__main__':
    import sys
    if '--selftest' not in sys.argv:
        raise SystemExit('usage: python -m full_run_28092026.symbolic.equivalence --selftest   # FREE')
    raise SystemExit(selftest())
