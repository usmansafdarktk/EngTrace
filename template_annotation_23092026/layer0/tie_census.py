"""Exact-arithmetic display-tie census, corpus-wide.

    python -m template_annotation_23092026.layer0.tie_census --seeds 500

Why this exists. T1 files an exact half-way tie as FAIL or MARGINAL by
floating-point luck (D-016: at a tie |evaluated - printed| == tol exactly, and
representation error decides which side it lands on), so T1's failure count
undercounts ties several-fold and a template can pass the gate while still
emitting ill-posed instances. This census re-evaluates every `... = <expr> =
<result>` line T1 can parse in Decimal arithmetic from the printed operands and
reports whether the exact value sits on a half-way tie at the printed precision,
independent of float luck. Adapted from the civil round-2 agent's per-template
script of 2026-09-23.

Coverage caveat. Only lines T1 can parse are censused (the same coverage T1
reports), and only expressions built from + - * / integer powers, sqrt and abs
are evaluated exactly; anything else is skipped and counted. A tie on a line T1
cannot read (a symbol in operand position, a unit such as `/hour` inside the
expression) is invisible here too.

Writes tie_census.json and tie_census.md beside this file. Informational: it
gates nothing, and the report says which templates carry residual ties.
"""
from __future__ import annotations

import argparse
import ast
import collections
import datetime as dt
import json
import subprocess
import sys
import time
from decimal import Decimal, getcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.checks.t1_closure import check_line  # noqa: E402
from tests.template_integrity.core import discover, generate, printed_precision  # noqa: E402

getcontext().prec = 60


def dec_eval(expr: str) -> Decimal:
    tree = ast.parse(expr, mode='eval')

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant):
            return Decimal(ast.get_source_segment(expr, n))
        if isinstance(n, ast.UnaryOp):
            v = ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                return a / b
            if isinstance(n.op, ast.Pow) and b == int(b):
                return a ** int(b)
            raise ValueError('operator')
        if isinstance(n, ast.Call):
            f = getattr(n.func, 'id', None)
            args = [ev(a) for a in n.args]
            if f == 'sqrt':
                return args[0].sqrt()
            if f == 'abs':
                return abs(args[0])
            raise ValueError('function')
        raise ValueError(type(n).__name__)
    return ev(tree)


def census(ref, seeds: int) -> dict:
    tie_seeds, checked, skipped = set(), 0, 0
    by_label = collections.Counter()
    examples = []
    for s in range(seeds):
        inst = generate(ref, s, capture=False)
        if not inst.ok:
            continue
        for raw in inst.solution.splitlines():
            line = raw.strip()
            if '=' not in line:
                continue
            checks, _ = check_line(line)
            for expr, _val, tok, _printed in checks:
                checked += 1
                try:
                    exact = dec_eval(expr)
                except Exception:                                  # noqa: BLE001
                    skipped += 1
                    continue
                scaled = exact * (Decimal(10) ** printed_precision(tok))
                if scaled - scaled.to_integral_value(rounding='ROUND_FLOOR') == Decimal('0.5'):
                    tie_seeds.add(s)
                    label = line.split('=')[0].strip()[:40]
                    by_label[label] += 1
                    if len(examples) < 2:
                        examples.append(f'seed {s}: {line[:120]}  (exact {exact})')
    return {'branch': ref.branch, 'checks': checked, 'skipped': skipped,
            'instances_with_tie': len(tie_seeds), 'rate': round(len(tie_seeds) / seeds, 4),
            'by_label': dict(by_label.most_common(5)), 'examples': examples}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=500)
    a = ap.parse_args()
    t0 = time.time()
    rows = {ref.template_id: census(ref, a.seeds) for ref in discover()}
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    report = {'generated': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'git_head': head,
              'seeds': a.seeds, 'rows': rows, 'seconds': round(time.time() - t0, 1)}
    (HERE / 'tie_census.json').write_text(json.dumps(report, indent=1), encoding='utf8')
    with_ties = {t: r for t, r in rows.items() if r['instances_with_tie']}
    with (HERE / 'tie_census.md').open('w', encoding='utf8') as fh:
        w = fh.write
        w('# Exact display-tie census\n\n')
        w(f'Generated {report["generated"]} by `tie_census.py` at git `{head[:10]}`, {a.seeds} seeds per '
          f'template, {report["seconds"]} s. A tie is an emitted line whose exact value, recomputed in Decimal '
          f'from the printed operands, sits exactly half-way at the printed precision; such an instance has no '
          f'gold value a decimal reader and a binary reader agree on (D-016). Only lines T1 can parse are '
          f'covered.\n\n')
        w('| | Templates |\n|---|---:|\n')
        w(f'| censused | {len(rows)} |\n| no tie in {a.seeds} instances | {len(rows) - len(with_ties)} |\n'
          f'| at least one tied instance | {len(with_ties)} |\n\n')
        if with_ties:
            w('| Template | Branch | Tied instances | Rate | Where |\n|---|---|---:|---:|---|\n')
            for t, r in sorted(with_ties.items(), key=lambda x: -x[1]['instances_with_tie']):
                where = '; '.join(f'`{k}` x{v}' for k, v in r['by_label'].items())
                w(f"| `{t}` | {r['branch'].split('_')[0]} | {r['instances_with_tie']} | {r['rate']:.1%} | {where} |\n")
            w('\nExamples:\n\n')
            for t, r in sorted(with_ties.items(), key=lambda x: -x[1]['instances_with_tie']):
                for ex in r['examples'][:1]:
                    w(f'- `{t}`: {ex}\n')
        skipped = sum(r['skipped'] for r in rows.values())
        checked = sum(r['checks'] for r in rows.values())
        w(f'\n{checked:,} checks evaluated exactly; {skipped:,} skipped (functions or operators outside the exact '
          f'evaluator, e.g. trigonometry and logarithms).\n')
    print(f'{len(rows)} templates, {len(with_ties)} with a tie at {a.seeds} seeds, {report["seconds"]} s')
    for t, r in sorted(with_ties.items(), key=lambda x: -x[1]['instances_with_tie'])[:20]:
        print(f"  {r['instances_with_tie']:4d} {r['rate']:6.1%}  {t}")


if __name__ == '__main__':
    main()
