"""Core infrastructure for the template integrity suite.

Discovery of the 150 templates, instrumented generation (capturing the
generator's frame locals at return), and shared value-matching utilities.

Read-only with respect to `data/templates/`. Nothing here modifies a template.
"""
from __future__ import annotations

import ast
import importlib
import os
import random
import re
import sys
from dataclasses import dataclass, field
from typing import Any, Callable, Iterator

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
BRANCHES = os.path.join(REPO_ROOT, 'data', 'templates', 'branches')

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class TemplateRef:
    """A single `template_*()` generator function."""
    template_id: str
    module: str
    file_path: str      # repo-relative, forward slashes
    branch: str
    area: str
    lineno: int

    def load(self) -> Callable[[], tuple[str, str]]:
        return getattr(importlib.import_module(self.module), self.template_id)

    def source(self) -> str:
        path = os.path.join(REPO_ROOT, self.file_path.replace('/', os.sep))
        lines = open(path, encoding='utf-8').read().splitlines()
        tree = ast.parse('\n'.join(lines))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == self.template_id:
                return '\n'.join(lines[node.lineno - 1:node.end_lineno])
        raise LookupError(self.template_id)


def discover(branch: str | None = None) -> list[TemplateRef]:
    """Every template_*() function under data/templates/branches/."""
    refs: list[TemplateRef] = []
    for dirpath, _dirnames, filenames in os.walk(BRANCHES):
        if '__pycache__' in dirpath:
            continue
        for fn in sorted(filenames):
            if not fn.endswith('.py') or fn in ('constants.py', '__init__.py'):
                continue
            path = os.path.join(dirpath, fn)
            rel = os.path.relpath(path, REPO_ROOT).replace(os.sep, '/')
            parts = rel.split('/')
            if len(parts) < 6:
                continue
            br, area = parts[3], parts[4]
            if branch and br != branch:
                continue
            src = open(path, encoding='utf-8').read()
            for node in ast.parse(src).body:
                if isinstance(node, ast.FunctionDef) and node.name.startswith('template_'):
                    refs.append(TemplateRef(
                        template_id=node.name,
                        module=rel[:-3].replace('/', '.'),
                        file_path=rel, branch=br, area=area, lineno=node.lineno,
                    ))
    refs.sort(key=lambda r: (r.branch, r.area, r.file_path, r.template_id))
    return refs


# --------------------------------------------------------------------------
# Instrumented generation
# --------------------------------------------------------------------------

@dataclass
class Instance:
    """One generated (question, solution) pair plus the generator's own values."""
    template_id: str
    seed: int
    question: str
    solution: str
    values: dict[str, float] = field(default_factory=dict)   # flattened numerics
    containers: dict[str, str] = field(default_factory=dict)  # name -> type name
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


NUMERIC_TYPES = ('int', 'float', 'Fraction', 'Decimal',
                 'float64', 'int64', 'float32', 'int32')


def _walk_values(obj: Any, name: str, out: dict[str, float],
                 containers: dict[str, str], depth: int = 0,
                 budget: list[int] | None = None) -> None:
    """Recursively flatten an object into {dotted_name: float}.

    Handles nested lists/dicts/ndarrays, complex (re/im), Fraction and Decimal.
    Bounded so a pathological structure cannot hang the suite.
    """
    if budget is None:
        budget = [4000]
    if budget[0] <= 0 or depth > 4:
        return
    tn = type(obj).__name__
    if isinstance(obj, bool) or obj is None or isinstance(obj, str):
        return
    if isinstance(obj, (int, float)) or tn in NUMERIC_TYPES:
        try:
            f = float(obj)
        except Exception:
            return
        if f == f and abs(f) != float('inf'):
            out[name] = f
            budget[0] -= 1
        return
    if isinstance(obj, complex):
        out[name + '.re'] = obj.real
        out[name + '.im'] = obj.imag
        containers.setdefault(name, 'complex')
        return
    if isinstance(obj, (list, tuple, set, frozenset)):
        containers.setdefault(name, tn)
        for i, e in enumerate(list(obj)[:80]):
            _walk_values(e, f'{name}[{i}]', out, containers, depth + 1, budget)
        return
    if isinstance(obj, dict):
        containers.setdefault(name, 'dict')
        for k, v in list(obj.items())[:80]:
            _walk_values(v, f'{name}[{k!r}]', out, containers, depth + 1, budget)
        return
    if tn == 'ndarray':
        containers.setdefault(name, 'ndarray')
        try:
            for i, e in enumerate(list(obj.ravel())[:120]):
                _walk_values(float(e), f'{name}[{i}]', out, containers, depth + 1, budget)
        except Exception:
            pass
        return
    containers.setdefault(name, tn)


def seed_all(seed: int, seed_numpy: bool = True) -> None:
    """Seed the generators a template might reach.

    `seed_numpy` matters, and getting it wrong silently disarms T3.

    The property T3 tests is that `random.seed(s)` ALONE determines the output.
    `template_levenspiel_plot_interpretation` calls `np.random.uniform`, which
    `random.seed()` does not control - the audit's finding is that the seed
    recorded beside every generated item cannot regenerate it. If the harness
    seeds numpy too, both of T3's processes reproduce that draw and the defect
    passes undetected.

    So: the output-based checks seed numpy (seed_numpy=True) purely so their own
    sampling is repeatable, while T3 passes seed_numpy=False, leaving numpy's
    state to differ between processes exactly as it does in production.
    """
    random.seed(seed)
    if seed_numpy:
        try:
            import numpy as _np
            _np.random.seed(seed % (2 ** 32 - 1))
        except Exception:
            pass


def generate(ref: TemplateRef, seed: int, capture: bool = True,
             seed_numpy: bool = True) -> Instance:
    """Generate one instance, optionally capturing the generator's frame locals."""
    fn = ref.load()
    seed_all(seed, seed_numpy=seed_numpy)
    captured: dict[str, Any] = {}

    if not capture:
        try:
            q, s = fn()
        except Exception as e:                       # noqa: BLE001 - reported, not raised
            return Instance(ref.template_id, seed, '', '', error=f'{type(e).__name__}: {e}')
        return Instance(ref.template_id, seed, q, s)

    target = ref.template_id

    def tracer(frame, event, _arg):
        if frame.f_code.co_name == target:
            if event == 'call':
                return tracer
            if event == 'return':
                captured['locals'] = dict(frame.f_locals)
        return tracer

    old = sys.gettrace()
    try:
        sys.settrace(tracer)
        q, s = fn()
    except Exception as e:                           # noqa: BLE001
        return Instance(ref.template_id, seed, '', '', error=f'{type(e).__name__}: {e}')
    finally:
        sys.settrace(old)

    values: dict[str, float] = {}
    containers: dict[str, str] = {}
    for k, v in captured.get('locals', {}).items():
        _walk_values(v, k, values, containers)
    return Instance(ref.template_id, seed, q, s, values, containers)


def generate_many(ref: TemplateRef, seeds: range | list[int],
                  capture: bool = True) -> Iterator[Instance]:
    for s in seeds:
        yield generate(ref, s, capture=capture)


# --------------------------------------------------------------------------
# Shared text / numeric helpers
# --------------------------------------------------------------------------

NUM_RE = re.compile(r'[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?')
STEP_RE = re.compile(r'\*\*Step\s*(\d+)\s*:\*\*')
STEP_LOOSE_RE = re.compile(r'\*\*Step\s*(\d+)\s*:')
#: Markers T4 can RECOGNISE.  `CANONICAL_ANSWER_MARKER` is the only one gold
#: may EMIT (D-059).  The `##` heading form is here because it is priority 0
#: on the candidate side (`normalize.ANSWER_MARKERS`) -- a solution carrying
#: it beside the canonical marker would have its answer extracted from the
#: heading, and until it was recognised here nothing could see that
#: (Reviewer A, F6).
ANSWER_MARKERS = ('**Answer:**', '**Answer**', '**Final Answer**',
                  '**Final Answers:**', '## Final Answer')
CANONICAL_ANSWER_MARKER = '**Answer:**'


def parse_number(tok: str) -> float | None:
    """Parse a printed number, tolerating thousands separators.

    NOTE: the evaluation parser's failure to strip separators BEFORE matching
    is the defect recorded in the audit (§9 Defect 1). This helper strips
    first, deliberately.
    """
    try:
        return float(tok.replace(',', ''))
    except ValueError:
        return None


def printed_precision(tok: str) -> int:
    """Decimal places shown in a printed token ('1.230' -> 3, '4,921' -> 0)."""
    body = tok.replace(',', '').split('e')[0].split('E')[0]
    return len(body.split('.')[1]) if '.' in body else 0


def display_tolerance(tok: str) -> float:
    """Half a unit in the last displayed digit."""
    return 0.5 * 10 ** (-printed_precision(tok))


def answer_block(solution: str) -> tuple[str, str | None]:
    """The trailing answer block and the marker that introduced it."""
    for marker in ANSWER_MARKERS:
        i = solution.find(marker)
        if i >= 0:
            return solution[i:], marker
    return '', None


def match_value(printed: str, values: dict[str, float],
                allow_scale: bool = True) -> tuple[str, str]:
    """Classify a printed number against the generator's bound values.

    Returns (kind, name) where kind is EXACT / ROUNDED / SCALED / MISSING.
    SCALED covers the lossless scaled-integer idiom (industrial branch) and
    display unit conversions, where the value is bound at a different scale.
    """
    val = parse_number(printed)
    if val is None:
        return 'MISSING', ''
    tol = max(display_tolerance(printed), abs(val) * 1e-9, 1e-12)
    rounded_hit = ''
    for name, v in values.items():
        if abs(v - val) <= 1e-9 * max(1.0, abs(v)):
            return 'EXACT', name
        if not rounded_hit and abs(v - val) <= tol:
            rounded_hit = name
    if rounded_hit:
        return 'ROUNDED', rounded_hit
    if allow_scale and val != 0:
        for name, v in values.items():
            if v == 0:
                continue
            for k in range(-12, 13):
                if k == 0:
                    continue
                w = v * (10.0 ** k)
                if abs(w - val) <= max(tol, abs(w) * 1e-9):
                    return 'SCALED', f'{name}*1e{k}'
    for name, v in values.items():
        if abs(-v - val) <= tol:
            return 'SCALED', f'-{name}'
    return 'MISSING', ''
