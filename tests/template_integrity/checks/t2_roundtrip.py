"""T2 - Round-trip oracle.

T1 checks that each printed line closes. T2 checks that the whole *chain*
closes: a solver who reads only the question, and uses only the values stated
there, must reach the printed answer.

That distinction is the reason T2 exists. In
`template_damping_classification` every individual printed line is arithmetically
fine, yet zeta recomputed from the stated 2-dp damping coefficient classifies
~35% of instances differently from the gold label. No line-level check can see
that; only re-deriving the answer from the question can.

The oracle contract
-------------------
One module per template in `tests/template_integrity/oracles/`, exposing:

    TEMPLATE_ID : str
    TOLERANCE   : float          # relative, justified in a comment
    SOURCE      : str            # the governing equations used, and where from
    def parse_givens(question: str) -> dict
    def recompute(givens: dict) -> float | dict | str
    def gold_answer(solution: str) -> float | dict | str

INDEPENDENCE IS THE POINT. An oracle written by transcribing the template's own
expression proves nothing - it will agree with that expression whether it is
right or wrong, which is the single most likely way this check gets silently
disarmed. Oracles must be derived from the governing equations (the docstring,
or the cited textbook), not from the computation code.
"""
from __future__ import annotations

import importlib
import os
import pkgutil
from dataclasses import dataclass, field
from typing import Any

from ..core import Instance, TemplateRef

ORACLE_PKG = 'tests.template_integrity.oracles'
ORACLE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'oracles')


@dataclass
class RoundTripFinding:
    template_id: str
    seed: int
    expected: Any          # what the oracle derived from the question
    got: Any               # what the trace printed as its answer
    rel_error: float | None
    detail: str = ''

    def __str__(self) -> str:
        err = f'{self.rel_error:.3%}' if self.rel_error is not None else 'n/a'
        return (f'[roundtrip] {self.template_id} seed={self.seed}: '
                f'question implies {self.expected!r}, trace answers {self.got!r} '
                f'(rel err {err}){" - " + self.detail if self.detail else ""}')


@dataclass
class RoundTripResult:
    template_id: str
    oracle: str = ''
    checked: int = 0
    unparsed: int = 0
    findings: list[RoundTripFinding] = field(default_factory=list)
    worst_rel_error: float = 0.0
    note: str = ''

    @property
    def passed(self) -> bool:
        return not self.findings and self.checked > 0

    @property
    def failure_rate(self) -> float:
        seeds = {f.seed for f in self.findings}
        return len(seeds) / self.checked if self.checked else 0.0


def available_oracles() -> dict[str, Any]:
    """Every oracle module, keyed by the template it covers."""
    out: dict[str, Any] = {}
    if not os.path.isdir(ORACLE_DIR):
        return out
    for mod in pkgutil.iter_modules([ORACLE_DIR]):
        if mod.name.startswith('_'):
            continue
        m = importlib.import_module(f'{ORACLE_PKG}.{mod.name}')
        tid = getattr(m, 'TEMPLATE_ID', None)
        if tid:
            out[tid] = m
    return out


def _compare(expected: Any, got: Any, tol: float) -> tuple[bool, float | None, str]:
    """Compare an oracle result against the trace's answer."""
    if isinstance(expected, str) or isinstance(got, str):
        ok = str(expected).strip().lower() == str(got).strip().lower()
        return ok, None, '' if ok else 'label mismatch'
    if isinstance(expected, dict) and isinstance(got, dict):
        worst, bad = 0.0, []
        for k, ev in expected.items():
            if k not in got:
                bad.append(f'missing part {k!r}')
                continue
            ok, err, _ = _compare(ev, got[k], tol)
            if err is not None:
                worst = max(worst, err)
            if not ok:
                bad.append(f'part {k!r}')
        return (not bad), worst, '; '.join(bad)
    try:
        e, g = float(expected), float(got)
    except (TypeError, ValueError):
        return False, None, 'uncomparable types'
    if e == 0:
        return abs(g) <= tol, abs(g), ''
    err = abs((g - e) / e)
    return err <= tol, err, ''


def run(ref: TemplateRef, instances: list[Instance],
        oracles: dict[str, Any] | None = None) -> RoundTripResult:
    if oracles is None:
        oracles = available_oracles()
    res = RoundTripResult(template_id=ref.template_id)
    oracle = oracles.get(ref.template_id)
    if oracle is None:
        res.note = 'no oracle'
        return res
    res.oracle = f'{ORACLE_PKG}.{oracle.__name__.rsplit(".", 1)[-1]}'
    tol = float(getattr(oracle, 'TOLERANCE', 0.005))

    for inst in instances:
        if not inst.ok:
            continue
        try:
            givens = oracle.parse_givens(inst.question)
            expected = oracle.recompute(givens)
            got = oracle.gold_answer(inst.solution)
        except Exception as exc:                      # noqa: BLE001
            res.unparsed += 1
            if res.unparsed <= 3:
                res.note = f'{type(exc).__name__}: {exc}'
            continue
        if expected is None or got is None:
            res.unparsed += 1
            continue
        res.checked += 1
        ok, err, detail = _compare(expected, got, tol)
        if err is not None:
            res.worst_rel_error = max(res.worst_rel_error, err)
        if not ok:
            res.findings.append(RoundTripFinding(
                ref.template_id, inst.seed, expected, got, err, detail))
    return res
