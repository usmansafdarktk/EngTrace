"""Dump the two Phase 3 templates' instances to JSON, for the before/after diff.

Run ONCE PER TREE, each as its OWN process, then diff the two files. That is how
the Phase 3 item-pool impact note and the D3.5 distribution diff are produced:

    git worktree add C:/wtm master
    python -m tests.template_integrity.phase3_instance_dump C:/wtm before.json 4000
    python -m tests.template_integrity.phase3_instance_dump .     after.json  4000

Separate processes are not optional, and neither is the explicit tree root. An
in-process reload resolves BOTH sides to the same already-imported
`data.templates.*` modules and reports every instance identical - the trap
documented in `instance_dump.py`, which has now caught three people. This script
takes the tree root as argv[1] and puts it at the FRONT of sys.path precisely so
that the mistake is visible in the command line rather than hidden in an import.

`master` does not carry this file, so run BOTH sides with the copy from the
Phase 3 branch, passing the master worktree as the root.

Not a check; gates nothing.
"""
import importlib
import json
import os
import random
import sys

TEMPLATES = {
    'template_normal_depth_iteration':
        'data.templates.branches.civil_engineering.water_resources.uniform_flow',
    'template_line_balancing_heuristic':
        'data.templates.branches.industrial_engineering.production_and_inventory'
        '.production_planning',
}


def main(argv):
    root = os.path.abspath(argv[1])
    out_path = argv[2]
    n = int(argv[3]) if len(argv) > 3 else 1000

    # Front of the path, and every competing copy of the package removed, so an
    # already-imported module from another tree cannot answer for this one.
    sys.path.insert(0, root)
    for name in [m for m in sys.modules if m.startswith('data.')]:
        del sys.modules[name]

    data = {}
    for tid, mod in TEMPLATES.items():
        module = importlib.import_module(mod)
        resolved = os.path.abspath(module.__file__)
        if not resolved.startswith(root):
            raise SystemExit(
                f'{tid} resolved to {resolved}, which is not under {root}. '
                f'Run each tree in its own process.')
        fn = getattr(module, tid)
        rows = []
        for s in range(n):
            random.seed(s)
            try:
                q, sol = fn()
                rows.append([q, sol])
            except Exception as exc:                      # noqa: BLE001
                rows.append(['', f'ERROR {type(exc).__name__}: {exc}'])
        data[tid] = rows

    json.dump({'root': root, 'n': n, 'instances': data},
              open(out_path, 'w', encoding='utf-8'))
    print(f'wrote {out_path}: {len(data)} templates x {n} seeds from {root}')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    main(sys.argv)
