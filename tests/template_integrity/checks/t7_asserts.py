"""T7 - Invariant asserts present.

Every template touched by a redesign phase must carry at least two
physical/mathematical `assert` statements in the civil/industrial style:
output bounds, an ordering or monotonicity relation, or a conservation
identity such as `assert abs((Ay + By) - (P + W)) < 0.02`.

The audit found 60 templates carrying 142 asserts between them - and all 60 are
in civil or industrial. Chemical, electrical and mechanical have **zero**. This
check therefore fails most of the corpus by design: it is a gate on templates a
phase edits, not a corpus-wide pass/fail. The runner applies it only to the
phase's scope, and reports corpus-wide coverage as information.

Asserts are also classified, because not all of them guard the emitted trace:
a *sampler* assert (`assert lo < hi, "empty window"`) guards the parameter draw
and says nothing about the output.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass, field

from ..core import TemplateRef


@dataclass
class AssertInfo:
    lineno: int
    source: str
    kind: str            # 'bounds' | 'ordering' | 'identity' | 'sampler' | 'other'
    message: str = ''


@dataclass
class AssertResult:
    template_id: str
    asserts: list[AssertInfo] = field(default_factory=list)
    min_required: int = 2

    @property
    def output_asserts(self) -> list[AssertInfo]:
        return [a for a in self.asserts if a.kind != 'sampler']

    @property
    def passed(self) -> bool:
        return len(self.output_asserts) >= self.min_required

    def summary(self) -> str:
        kinds = {}
        for a in self.asserts:
            kinds[a.kind] = kinds.get(a.kind, 0) + 1
        return f'{len(self.asserts)} asserts {kinds or "{}"}'


_SAMPLER_HINTS = ('window', 'exhaust', 'resample', 'empty', 'range is empty',
                  'no valid', 'combos', 'table', 'mismatch', 'out of sync')


def _classify(node: ast.Assert, msg: str) -> str:
    test = node.test
    low = msg.lower()
    if any(h in low for h in _SAMPLER_HINTS):
        return 'sampler'
    if isinstance(test, ast.Compare):
        # `lo <= x <= hi` - two comparators on one target
        if len(test.ops) >= 2:
            return 'bounds'
        op = test.ops[0]
        if isinstance(op, (ast.Lt, ast.LtE, ast.Gt, ast.GtE)):
            # abs(a - b) < tol is a conservation identity, not a bound
            left = test.left
            if isinstance(left, ast.Call) and getattr(left.func, 'id', '') == 'abs':
                return 'identity'
            if isinstance(left, ast.Name) and isinstance(test.comparators[0], ast.Name):
                return 'ordering'
            return 'bounds'
        if isinstance(op, (ast.Eq, ast.NotEq)):
            return 'identity'
    if isinstance(test, ast.BoolOp):
        return 'ordering'
    return 'other'


def run(ref: TemplateRef, min_required: int = 2) -> AssertResult:
    res = AssertResult(template_id=ref.template_id, min_required=min_required)
    src = ref.source()
    tree = ast.parse(src)
    lines = src.splitlines()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assert):
            continue
        msg = ''
        if node.msg is not None:
            try:
                msg = ast.unparse(node.msg)
            except Exception:                             # noqa: BLE001
                msg = ''
        text = lines[node.lineno - 1].strip() if node.lineno <= len(lines) else ''
        res.asserts.append(AssertInfo(
            lineno=ref.lineno + node.lineno - 1,
            source=text[:160],
            kind=_classify(node, msg),
            message=msg[:120],
        ))
    return res
