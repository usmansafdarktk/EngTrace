"""C3.4 - dump every template's instances from ONE tree, then diff two dumps.

    git -c core.longpaths=true worktree add --detach <path> <master-ref>
    python -m tests.constants_integrity.c3_instance_dump <master-tree> before.json 300
    python -m tests.constants_integrity.c3_instance_dump .             after.json  300
    python -m tests.constants_integrity.c3_instance_dump --diff before.json after.json

Run each tree in its OWN process, passing its root: an in-process reload resolves
both sides to the same already-imported modules and reports every instance
identical - the trap that has caught four people (phase5_instance_dump.py). The
root goes at the front of sys.path, and every template module is checked to have
resolved under it.

The diff answers the two P6 questions per template - DID THE QUESTION CHANGE, and
DID THE ANSWER CHANGE - because a constants correction can do either: a restated
value changes the question (re-inference), a hidden one only the answer (re-score).

Not a check; gates nothing.
"""
import json
import os
import random
import re
import sys

MARKERS = ('**Answer:**', '**Answer**')


def dump(root, out_path, n):
    root = os.path.abspath(root)
    sys.path.insert(0, root)
    for m in [m for m in sys.modules if m.startswith(('data.', 'tests.'))]:
        del sys.modules[m]
    from tests.template_integrity.core import discover
    data = {}
    for ref in discover():
        fn = ref.load()
        mod = sys.modules[ref.module]
        if not os.path.abspath(mod.__file__).startswith(root):
            raise SystemExit(f'{ref.module} resolved to {mod.__file__}, not under {root}')
        rows = []
        for s in range(n):
            random.seed(s)
            try:
                q, sol = fn()
                rows.append([q, sol])
            except Exception as exc:                          # noqa: BLE001
                rows.append(['', f'ERROR {type(exc).__name__}: {exc}'])
        data[ref.template_id] = rows
    json.dump({'root': root, 'n': n, 'instances': data}, open(out_path, 'w', encoding='utf-8'))
    print(f'wrote {out_path}: {len(data)} templates x {n} seeds from {root}')


def answer_body(sol):
    cut = max((sol.rfind(mk) for mk in MARKERS), default=-1)
    return sol[cut:].split('**', 2)[-1].strip() if cut >= 0 else sol


def diff(before_path, after_path, only=None):
    A = json.load(open(before_path, encoding='utf-8'))
    B = json.load(open(after_path, encoding='utf-8'))
    if A['root'] == B['root']:
        raise SystemExit('both dumps came from one tree - the in-process trap, not a null result')
    print(f"before: {A['root']} (n={A['n']})\nafter : {B['root']} (n={B['n']})\n")
    moved = 0
    tot = dict(q=0, ans=0, sol=0, err_before=0, err_after=0)
    for tid in sorted(B['instances']):
        if only and not re.search(only, tid):
            continue
        a, b = A['instances'].get(tid), B['instances'][tid]
        if a is None:
            print(f'{tid}: only in AFTER')
            continue
        q = sum(x[0] != y[0] for x, y in zip(a, b))
        s = sum(x[1] != y[1] for x, y in zip(a, b))
        ans = sum(answer_body(x[1]) != answer_body(y[1]) for x, y in zip(a, b))
        eb = sum(x[1].startswith('ERROR') for x in a)
        ea = sum(y[1].startswith('ERROR') for y in b)
        for k, v in (('q', q), ('ans', ans), ('sol', s), ('err_before', eb), ('err_after', ea)):
            tot[k] += v
        if q or s or eb != ea:
            moved += 1
            print(f'{tid:48s} question {q:4d}/{len(b)}  answer {ans:4d}/{len(b)}  '
                  f'solution {s:4d}/{len(b)}  errors {eb}->{ea}')
            for i, (x, y) in enumerate(zip(a, b)):
                if x != y:
                    print(f'    seed {i} answer: {answer_body(x[1])[:70]!r} -> {answer_body(y[1])[:70]!r}')
                    break
    print(f"\n{moved} templates moved; totals {tot}")
    if tot['err_after'] > tot['err_before']:
        print('!! MORE GENERATION ERRORS AFTER THAN BEFORE')


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == '--diff':
        diff(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    elif len(sys.argv) >= 3:
        dump(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 300)
    else:
        raise SystemExit(__doc__)
