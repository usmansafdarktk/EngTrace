"""Lift a template's structured trace node out of a generation call (D3.3).

    python -m tests.trace_schema.extract corpus.json 40

writes a JSON corpus of {seed, question, solution, trace_nodes} rows for the two
Phase 3 templates - the input a verifier is written against.

How it works, and why it is not a second implementation
-------------------------------------------------------
A template returns `(question, solution)` strings; that contract is unchanged.
Internally it builds its trace as a dict bound to the local name `trace_nodes`
and RENDERS the prose from that dict. This module sets a trace function, catches
the template frame's `return` event, and reads the local. So the node it yields
is the same object the printed trace came from - it cannot drift from the prose,
and this file holds no reference values of its own (D-034).

The cost of that choice is that the node is only reachable through the frame,
which is why this extractor exists rather than a plain function call. When the
milestone model lands (Phase 4+) the node should be promoted to a return value
and this module becomes a compatibility shim.
"""
from __future__ import annotations

import importlib
import json
import random
import sys

TEMPLATES = {
    'template_normal_depth_iteration':
        'data.templates.branches.civil_engineering.water_resources.uniform_flow',
    'template_line_balancing_heuristic':
        'data.templates.branches.industrial_engineering.production_and_inventory'
        '.production_planning',
}

LOCAL_NAME = 'trace_nodes'


def generate_with_nodes(template_id: str, seed: int) -> dict:
    """Generate one instance and return its prose together with its trace node."""
    module = importlib.import_module(TEMPLATES[template_id])
    fn = getattr(module, template_id)
    captured: dict = {}

    def tracer(frame, event, _arg):
        if frame.f_code.co_name == template_id:
            if event == 'call':
                return tracer
            if event == 'return':
                captured['locals'] = dict(frame.f_locals)
        return tracer

    old = sys.gettrace()
    try:
        sys.settrace(tracer)
        random.seed(seed)
        question, solution = fn()
    finally:
        sys.settrace(old)

    nodes = captured.get('locals', {}).get(LOCAL_NAME)
    if nodes is None:
        raise LookupError(
            f'{template_id} bound no local named {LOCAL_NAME!r}; the template '
            f'is not emitting a structured trace node')
    return {'template_id': template_id, 'seed': seed, 'question': question,
            'solution': solution, 'trace_nodes': nodes}


def main(argv: list[str]) -> None:
    out_path = argv[1]
    n = int(argv[2]) if len(argv) > 2 else 25
    rows = [generate_with_nodes(tid, s)
            for tid in TEMPLATES for s in range(n)]
    json.dump(rows, open(out_path, 'w', encoding='utf-8'), indent=1)
    print(f'wrote {out_path}: {len(rows)} rows '
          f'({len(TEMPLATES)} templates x {n} seeds)')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv)
