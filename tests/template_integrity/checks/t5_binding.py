"""T5 - Value binding and rounding discipline (static, AST).

Two sub-checks, both operating on the template source rather than its output.

T5a - **binding**. A quantity printed in *result position* (the interpolation
immediately following an `=` in the format string) must be a bound variable,
not an expression computed inside the f-string. Operand restatements such as
`{2*L}` or `{L/2}` are legitimate and are counted separately, never failed:
they restate a given, they do not report a result.

T5b - **rounding discipline (P2)**. A variable that is displayed at N decimal
places and then consumed by a later computation must be rounded to N places
*before* that consumption, so the stored value and the printed value agree.

T5b exists because T1 cannot see this class of defect. In
`template_cantilever_double_integration`, `delta` is printed at 5 dp and then
`delta_mm = round(delta * 1000, 1)` is computed from the unrounded value; every
individual printed line still closes to within its display tolerance, so
printed-arithmetic closure passes while the trace is wrong 4.95% of the time.
The audit's other two closure-invisible defects (`template_mean_variance`,
`template_impulse_response_from_lccde`) are the same shape.
"""
from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field

from ..core import TemplateRef

# Wrappers that pass a bound value through without computing a new quantity.
_PASSTHROUGH = {'round', 'str', 'int', 'float', 'abs', 'format',
                'capitalize', 'upper', 'lower', 'title', 'strip', 'replace'}

# A literal chunk ending in '=' (optionally with markdown/parens) puts the next
# interpolation in result position.
_RESULT_POSITION = re.compile(r'=\s*[(\[*]*\s*$')

# Format specs that fix a decimal precision.
_PRECISION_SPEC = re.compile(r'\.(\d+)[feEg]')


def _is_bound(node: ast.AST) -> bool:
    """Is this interpolation a bound value rather than a computation?"""
    if isinstance(node, (ast.Name, ast.Constant, ast.Attribute, ast.Subscript)):
        return True
    if isinstance(node, ast.IfExp):
        return _is_bound(node.body) and _is_bound(node.orelse)
    if isinstance(node, ast.Call):
        fname = getattr(node.func, 'id', None) or getattr(node.func, 'attr', None)
        if fname in _PASSTHROUGH:
            return all(_is_bound(a) for a in node.args) if node.args else True
    return False


@dataclass
class BindingFinding:
    template_id: str
    lineno: int
    expression: str
    kind: str            # 'result-inline' | 'unrounded-consumed'
    detail: str = ''

    def __str__(self) -> str:
        return (f'[{self.kind}] {self.template_id}:{self.lineno}  '
                f'{self.expression[:80]}' + (f'  ({self.detail})' if self.detail else ''))


@dataclass
class BindingResult:
    template_id: str
    interpolations: int = 0
    result_position: int = 0
    operand_restatements: int = 0
    rounding_candidates: int = 0   # static candidates dismissed by runtime data
    findings: list[BindingFinding] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.findings

    @property
    def inline_results(self) -> list[BindingFinding]:
        return [f for f in self.findings if f.kind == 'result-inline']

    @property
    def rounding_violations(self) -> list[BindingFinding]:
        return [f for f in self.findings if f.kind == 'unrounded-consumed']


class _Analyser(ast.NodeVisitor):
    def __init__(self, template_id: str, res: BindingResult):
        self.tid = template_id
        self.res = res
        # name -> was every assignment to it wrapped in round()/quantize()?
        self.assigned_rounded: dict[str, bool] = {}
        # name -> precisions it is printed at
        self.printed_precision: dict[str, set[int]] = {}
        # names consumed on the RHS of a later assignment
        self.consumed_in_assign: set[str] = set()
        self.print_line: dict[str, int] = {}

    # ---- assignments
    def visit_Assign(self, node: ast.Assign) -> None:
        rhs = node.value
        rounded = False
        if isinstance(rhs, ast.Call):
            fname = getattr(rhs.func, 'id', None) or getattr(rhs.func, 'attr', None)
            rounded = fname in ('round', 'quantize', '_hu', '_r2', '_r4', '_q4')
        for t in node.targets:
            for n in ast.walk(t):
                if isinstance(n, ast.Name):
                    prev = self.assigned_rounded.get(n.id, True)
                    self.assigned_rounded[n.id] = prev and rounded
        for n in ast.walk(rhs):
            if isinstance(n, ast.Name):
                self.consumed_in_assign.add(n.id)
        self.generic_visit(node)

    # ---- f-strings
    def visit_JoinedStr(self, node: ast.JoinedStr) -> None:
        prev_literal = ''
        for part in node.values:
            if isinstance(part, ast.Constant) and isinstance(part.value, str):
                prev_literal = part.value
                continue
            if not isinstance(part, ast.FormattedValue):
                continue
            self.res.interpolations += 1
            try:
                expr = ast.unparse(part.value)
            except Exception:                             # noqa: BLE001
                expr = '<unparseable>'
            in_result_pos = bool(_RESULT_POSITION.search(prev_literal))
            bound = _is_bound(part.value)

            if in_result_pos:
                self.res.result_position += 1
                if not bound:
                    self.res.findings.append(BindingFinding(
                        self.tid, getattr(part.value, 'lineno', node.lineno),
                        expr, 'result-inline',
                        'computed inside the f-string in result position'))
            elif not bound:
                self.res.operand_restatements += 1

            # record printed precision for T5b
            if isinstance(part.value, ast.Name) and part.format_spec is not None:
                try:
                    spec = ast.unparse(part.format_spec)
                except Exception:                          # noqa: BLE001
                    spec = ''
                m = _PRECISION_SPEC.search(spec)
                if m:
                    self.printed_precision.setdefault(part.value.id, set()).add(int(m.group(1)))
                    self.print_line.setdefault(part.value.id, node.lineno)
            prev_literal = ''
        self.generic_visit(node)


def run(ref: TemplateRef, instances: list | None = None) -> BindingResult:
    """Static analysis, with T5b confirmed against runtime values.

    `instances` should be generated with capture=True. Without them T5b reports
    only static *candidates*, which over-reports heavily: a constant such as
    `gamma_w = 9.81` printed at 2 dp and reused is a candidate but not a
    violation, because rounding it to its display precision changes nothing.
    A violation requires that the displayed precision actually discards
    information the downstream computation then uses.
    """
    res = BindingResult(template_id=ref.template_id)
    tree = ast.parse(ref.source())
    fn = tree.body[0]
    an = _Analyser(ref.template_id, res)
    for stmt in getattr(fn, 'body', []):
        an.visit(stmt)

    for name, precisions in an.printed_precision.items():
        if an.assigned_rounded.get(name, True):
            continue                                # rounded at assignment: fine
        if name not in an.consumed_in_assign:
            continue                                # printed only: harmless
        prec = min(precisions)
        confirmed, checked, worst = False, 0, 0.0
        if instances is not None:
            for inst in instances:
                if not inst.ok or name not in inst.values:
                    continue
                checked += 1
                v = inst.values[name]
                delta = abs(round(v, prec) - v)
                if delta > 0:
                    confirmed = True
                    worst = max(worst, delta)
            if not confirmed:
                res.rounding_candidates += 1
                continue                            # display rounding is lossless
        detail = (f'printed at {sorted(precisions)} dp, assigned unrounded, '
                  f'then consumed by a later computation')
        if checked:
            detail += (f'; display rounding discards up to {worst:.3g} '
                       f'on {checked} sampled instances')
        res.findings.append(BindingFinding(
            ref.template_id, an.print_line.get(name, ref.lineno), name,
            'unrounded-consumed', detail))
    return res
