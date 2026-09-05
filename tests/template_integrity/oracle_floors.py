"""Measured agreement floor and detection floor for each Phase 1 oracle.

Reviewer A, finding F-4: each oracle's declared TOLERANCE is argued in prose,
and in several cases it is not the number that actually binds - both sides are
quantised to the trace's printed precision, so the display step sets the
sensitivity. This measures both ends instead of arguing them:

  agreement floor - the worst relative disagreement between oracle and trace on
                    a CORRECT instance. The declared tolerance must exceed it.
  detection floor - the smallest uniform relative error injected into the gold
                    answer that the oracle catches on >= 99% of instances. This
                    is what a downstream verifier can actually rely on.
"""
import re
import sys

from tests.template_integrity.core import discover, generate, Instance
import tests.template_integrity.checks.t2_roundtrip as T2

ORACLES = T2.available_oracles()
NUM = re.compile(r'(?<![\d.])(\d+\.\d+|\d+)(?![\d])')
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
TARGETS = sys.argv[2:] or sorted(ORACLES)


def bump(sol, factor):
    i = sol.find('**Answer')
    if i < 0:
        return sol

    def f(m):
        v = float(m.group(1))
        if v == 0:
            return m.group(0)
        dp = len(m.group(1).split('.')[1]) if '.' in m.group(1) else 0
        return f'{v * factor:.{dp}f}'
    return sol[:i] + NUM.sub(f, sol[i:])


print(f'{"template":42s} {"tol":>8s} {"agreement":>10s} {"detection floor":>16s}')
for tid in TARGETS:
    refs = [r for r in discover() if r.template_id == tid]
    if not refs:
        continue
    ref = refs[0]
    base = [generate(ref, s, capture=False) for s in range(N)]
    clean = T2.run(ref, base, ORACLES)
    if clean.checked == 0:
        continue
    tol = float(getattr(ORACLES[tid], 'TOLERANCE', 0.005))
    floor = None
    for pct in (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
        mut = [Instance(i.template_id, i.seed, i.question, bump(i.solution, 1 + pct))
               for i in base]
        r = T2.run(ref, mut, ORACLES)
        if r.checked and len(r.findings) / r.checked >= 0.99:
            floor = pct
            break
    print(f'{tid:42s} {tol:8.1e} {clean.worst_rel_error:10.2e} '
          f'{("%.2f%%" % (100 * floor)) if floor else ">5%":>16s}'
          f'{"   <-- TOLERANCE BELOW AGREEMENT" if clean.worst_rel_error > tol else ""}')
