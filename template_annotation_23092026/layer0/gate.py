"""Layer 0 - the deterministic admission gate, and its report.

    python -m template_annotation_23092026.layer0.gate                 # 500 seeds
    python -m template_annotation_23092026.layer0.gate --seeds 2000

GATING checks, every template must pass: T1 closure (the printed chain reproduces
its own numbers), T3 determinism (a seed regenerates the instance in a second
process), T4 output contract, T8 emission hygiene, and no generation error.
T2 (round-trip against an oracle) is reported where an oracle exists.
ADVISORY, reported but never gating: T5 value binding and T7 invariant asserts,
which are instrumentation and style gates on templates a phase edits
(phase6_residual_register.md), not corpus-wide correctness checks.

T1 is run here directly rather than through the suite's runner, because the
runner's JSON keeps only two example failures per template and the register
below needs every failing line: a closure failure is EXCUSED only when the line
matches a pattern registered for that template in check_limits.json, with a
reason a reader can check. Anything else in that template still fails. The
register is the residual list for this layer; a growing register is a warning,
so the report prints every entry and how many lines it absorbed.

Writes gate_report.json and gate_report.md beside this file.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.checks import t1_closure  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402

GATING = ('T1', 'T3', 'T4', 'T8')
ADVISORY = ('T5', 'T7')
REGISTER = HERE / 'check_limits.json'


def load_register() -> dict[str, list[dict]]:
    reg = json.loads(REGISTER.read_text(encoding='utf8'))
    out = collections.defaultdict(list)
    for e in reg['entries']:
        e['_re'] = re.compile(e['line_pattern'])
        out[e['template_id']].append(e)
    return out


def run_suite(seeds: int, out: Path) -> dict:
    """T2-T8 through the suite's own runner; its exit code is ignored (T5 is advisory)."""
    cmd = [sys.executable, '-m', 'tests.template_integrity.run', '--checks', 'T2,T3,T4,T5,T7,T8',
           '--seeds', str(seeds), '--json', str(out), '--quiet']
    subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    return json.loads(out.read_text(encoding='utf8'))


def run_t1(seeds: int, register: dict) -> dict:
    res = {}
    for ref in discover():
        inst = [generate(ref, s, capture=False) for s in range(seeds)]
        errs = [i.error for i in inst if not i.ok]
        r = t1_closure.run(inst, ref.template_id)
        excused, unexcused = collections.Counter(), []
        for f in r.failures:
            hit = next((e for e in register.get(ref.template_id, []) if e['_re'].search(f.line)), None)
            if hit:
                excused[hit['line_pattern']] += 1
            else:
                unexcused.append(str(f))
        res[ref.template_id] = {
            'branch': ref.branch, 'generation_errors': errs[:3],
            'failures': len(r.failures), 'excused': dict(excused), 'unexcused': len(unexcused),
            'marginals': len(r.marginals), 'checks': r.checks_performed,
            'unit_rescales': r.unit_rescales, 'coverage': round(r.coverage, 3),
            'examples': unexcused[:3],
        }
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=500)
    a = ap.parse_args()
    t0 = time.time()
    register = load_register()
    suite = run_suite(a.seeds, HERE / '_suite.json')
    t1 = run_t1(a.seeds, register)
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()

    rows = {}
    for tid, s in suite.items():
        c = t1[tid]
        fails = []
        if c['generation_errors']:
            fails.append('generation error')
        if c['unexcused']:
            fails.append('T1')
        for chk in ('T3', 'T4', 'T8'):
            if chk in s and not s[chk]['pass']:
                fails.append(chk)
        rows[tid] = {'branch': s['branch'], 'gate': 'pass' if not fails else 'FAIL', 'failing': fails,
                     'T1': {k: c[k] for k in ('failures', 'excused', 'unexcused', 'marginals', 'checks',
                                              'unit_rescales', 'coverage', 'examples')},
                     'T2': s.get('T2'), 'T3': s.get('T3'), 'T4': s.get('T4'), 'T8': s.get('T8'),
                     'T5': s.get('T5', {}).get('pass'), 'T7': s.get('T7', {}).get('pass')}

    n = len(rows)
    failing = {t: r for t, r in rows.items() if r['gate'] == 'FAIL'}
    per_check = {chk: sum(1 for r in rows.values() if chk in r['failing']) for chk in ('T1', 'T3', 'T4', 'T8')}
    per_check['generation error'] = sum(1 for r in rows.values() if 'generation error' in r['failing'])
    excused_total = collections.Counter()
    for t, r in rows.items():
        for pat, k in r['T1']['excused'].items():
            excused_total[(t, pat)] += k
    t2 = [r['T2'] for r in rows.values() if r['T2'] and r['T2'].get('note') != 'no oracle']
    adv = {chk: sum(1 for r in rows.values() if r[chk] is False) for chk in ADVISORY}
    by_branch = collections.defaultdict(lambda: [0, 0])
    for r in rows.values():
        by_branch[r['branch']][r['gate'] == 'pass'] += 1

    report = {'generated': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'git_head': head,
              'seeds': a.seeds, 'gating': GATING, 'advisory': ADVISORY, 'templates': n,
              'gate_pass': n - len(failing), 'gate_fail': len(failing), 'per_check': per_check,
              'advisory_failing': adv, 'register_absorbed': {f'{t}::{p}': k for (t, p), k in excused_total.items()},
              'rows': rows, 'seconds': round(time.time() - t0, 1)}
    (HERE / 'gate_report.json').write_text(json.dumps(report, indent=1), encoding='utf8')

    with (HERE / 'gate_report.md').open('w', encoding='utf8') as fh:
        w = fh.write
        w('# Layer 0 gate report\n\n')
        w(f'Generated {report["generated"]} by `gate.py` at git `{head[:10]}`, {a.seeds} seeds per template, '
          f'{report["seconds"]} s. Gating: {", ".join(GATING)} plus generation errors. '
          f'Advisory: {", ".join(ADVISORY)}.\n\n')
        w('| | Templates |\n|---|---:|\n')
        w(f'| in corpus | {n} |\n| **pass the gate** | **{n - len(failing)}** |\n| fail the gate | {len(failing)} |\n')
        for chk, k in per_check.items():
            w(f'| failing {chk} | {k} |\n')
        w(f'| T1 closure lines excused by the register | {sum(excused_total.values())} |\n')
        w(f'| T2 oracles run (informational) | {len(t2)} ({sum(1 for x in t2 if x["pass"])} pass) |\n')
        for chk, k in adv.items():
            w(f'| {chk} advisory, not passing | {k} |\n')
        w('\n| Branch | pass | fail |\n|---|---:|---:|\n')
        for b in sorted(by_branch):
            w(f'| {b} | {by_branch[b][1]} | {by_branch[b][0]} |\n')
        if failing:
            w('\n## Templates failing the gate\n\n')
            for t, r in sorted(failing.items()):
                w(f"- `{t}` [{r['branch']}]: {', '.join(r['failing'])}")
                if r['T1']['examples']:
                    w('\n  ```\n  ' + r['T1']['examples'][0].replace('\n', '\n  ') + '\n  ```')
                w('\n')
        w('\n## Register of accepted check limits\n\n')
        w('Closure failures excused because the check cannot read the line, not because the line is wrong. '
          'Each entry names the template, the line pattern and the reason; the count is how many lines it '
          f'absorbed at {a.seeds} seeds.\n\n')
        w('| Template | Pattern | Lines absorbed | Reason |\n|---|---|---:|---|\n')
        for e in json.loads(REGISTER.read_text(encoding='utf8'))['entries']:
            k = excused_total.get((e['template_id'], e['line_pattern']), 0)
            w(f"| `{e['template_id']}` | `{e['line_pattern']}` | {k} | {e['reason']} |\n")
        w('\n## Closure marginals (informational)\n\n')
        marg = sorted(((r['T1']['marginals'] / r['T1']['checks'] if r['T1']['checks'] else 0, t)
                       for t, r in rows.items()), reverse=True)
        w(f'{sum(1 for m, _ in marg if m > 0)} templates carry at least one marginal line '
          f'(within tolerance, on a rounding boundary). Highest rates: '
          + ', '.join(f'`{t}` {m:.0%}' for m, t in marg[:5]) + '.\n')
    print(f'{n} templates: {n - len(failing)} pass the gate, {len(failing)} fail {dict(per_check)}; '
          f'register absorbed {sum(excused_total.values())} lines; advisory {adv}; {report["seconds"]} s')
    for t, r in sorted(failing.items()):
        print('  FAIL', t, r['failing'])
    (HERE / '_suite.json').unlink(missing_ok=True)


if __name__ == '__main__':
    main()
