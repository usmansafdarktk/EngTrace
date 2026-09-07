"""Dump the four Phase 4 templates' instances to JSON, for the D4.5 before/after.

Run ONCE PER TREE, each as its OWN process, then diff the two files:

    git worktree add C:/wtm master
    python -m tests.template_integrity.phase4_instance_dump C:/wtm before.json 4000
    python -m tests.template_integrity.phase4_instance_dump .     after.json  4000
    python -m tests.template_integrity.phase4_instance_dump --diff before.json after.json

Separate processes are not optional, and neither is the explicit tree root. An
in-process reload resolves BOTH sides to the same already-imported
`data.templates.*` modules and reports every instance identical - the trap
documented in `instance_dump.py`, which has now caught three people. This script
takes the tree root as argv[1] and puts it at the FRONT of sys.path precisely so
that the mistake is visible in the command line rather than hidden in an import.

`master` does not carry this file, so run BOTH sides with the copy from the
Phase 4 branch, passing the master worktree as the root.

The `--diff` mode classifies every changed instance into the two halves of the
D4.5 fix - the bracket, which changes every instance, and the variable, which
changes only reversal instances - so the two can be approved or reverted
separately. Copied from `phase3_instance_dump.py`; not a check, gates nothing.
"""
import difflib
import importlib
import json
import os
import random
import sys

TEMPLATES = {
    'template_signal_operations':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.discrete_time_signals',
    'template_system_properties_memory_causality':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.discrete_time_signals',
    'template_system_property_linearity':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.discrete_time_signals',
    'template_incompressible_continuity':
        'data.templates.branches.mechanical_engineering.fluid_mechanics'
        '.fluid_kinematics',
}


def dump(argv):
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


def diff(before_path, after_path):
    before = json.load(open(before_path, encoding='utf-8'))
    after = json.load(open(after_path, encoding='utf-8'))
    if before['root'] == after['root']:
        raise SystemExit(
            f"both files were dumped from {before['root']} - that is the "
            f"in-process-reload trap, not a null result. Re-run each tree "
            f"separately.")
    print(f"before: {before['root']}  n={before['n']}")
    print(f"after : {after['root']}  n={after['n']}")
    print()

    grand = {'questions_changed': 0, 'answers_changed': 0,
             'solutions_changed': 0, 'bracket_only': 0, 'variable_too': 0}
    for tid in TEMPLATES:
        b, a = before['instances'][tid], after['instances'][tid]
        q_ch = s_ch = ans_ch = 0
        bracket_only = variable_too = 0
        examples = []
        for i, (bi, ai) in enumerate(zip(b, a)):
            if bi[0] != ai[0]:
                q_ch += 1
            if bi[1] == ai[1]:
                continue
            s_ch += 1
            if _answer(bi[1]) != _answer(ai[1]):
                ans_ch += 1
            # Classify: did any changed line mention a variable other than y?
            changed = [ln for ln in _changed_lines(bi[1], ai[1])]
            if any("z[k']" in ln for ln in changed):
                variable_too += 1
            else:
                bracket_only += 1
            if len(examples) < 2:
                examples.append((i, changed))
        grand['questions_changed'] += q_ch
        grand['solutions_changed'] += s_ch
        grand['answers_changed'] += ans_ch
        grand['bracket_only'] += bracket_only
        grand['variable_too'] += variable_too
        print(f"{tid}")
        print(f"    questions changed : {q_ch}/{len(b)}")
        print(f"    answers changed   : {ans_ch}/{len(b)}")
        print(f"    solutions changed : {s_ch}/{len(b)}")
        if s_ch:
            print(f"      bracket only    : {bracket_only}")
            print(f"      + variable y->z : {variable_too}")
        for i, changed in examples:
            print(f"    seed {i}:")
            for ln in changed[:4]:
                print(f"      {ln}")
        print()

    print("SUMMARY")
    for k, v in grand.items():
        print(f"  {k:20s} {v}")
    print()
    if grand['questions_changed'] == 0 and grand['answers_changed'] == 0:
        print("ITEM POOL UNCHANGED: no question and no answer differs on any seed.")
        print("The only difference is one column header inside the solution table.")
    else:
        print("ITEM POOL CHANGED - this is a P6 event and needs sign-off.")


def _answer(sol):
    i = sol.rfind('**Answer:**')
    return sol[i:] if i >= 0 else sol


def _changed_lines(b, a):
    out = []
    for ln in difflib.unified_diff(b.split('\n'), a.split('\n'), n=0, lineterm=''):
        if ln.startswith(('+++', '---', '@@')):
            continue
        if ln.startswith(('+', '-')):
            out.append(ln.rstrip())
    return out


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == '--diff':
        diff(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 3:
        dump(sys.argv)
    else:
        raise SystemExit(__doc__)
