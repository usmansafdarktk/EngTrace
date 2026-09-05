"""Exit-gate verification at scale, for the templates a phase edits.

Runs T1 (printed-arithmetic closure) and T2 (round-trip oracle) over a large
seed range in chunks, and T4/T5/T7 over a small one, reporting generation
errors separately - a template that RAISES produces no output for T1 and T2 to
fail on, so an all-errors run reads as "0 failures" unless the error count is
surfaced beside it.

    python -m tests.template_integrity.gate_report 60000

Not a check: `run.py` owns the gate. This exists because the seed count that
resolves a defect is often far above the runner's default, and Phase 1 shipped
a template verified at 1,000 seeds that failed at 20,000 (DECISIONS D-024).
"""
import sys
from tests.template_integrity.core import discover, generate
from tests.template_integrity.checks import t1_closure, t4_contract, t5_binding, t7_asserts
import tests.template_integrity.checks.t2_roundtrip as T2

TWELVE = ['template_mean_variance','template_rotating_unbalance',
 'template_vibration_transmissibility','template_damping_classification',
 'template_beam_deflection_formula','template_cantilever_double_integration',
 'template_annulus_flowrate','template_statically_indeterminate_shaft',
 'template_shaft_design_power','template_composite_shafts_series']
N = int(sys.argv[1]) if len(sys.argv)>1 else 20000
refs = {r.template_id: r for r in discover()}
orc = T2.available_oracles()
print(f'{"template":42s} {"T1 fail":>8s} {"cov":>6s} {"T2 fail":>8s} {"worst":>10s} {"T4":>4s} {"T5":>4s} {"T7":>4s} {"err":>4s}')
bad = 0
for tid in TWELVE:
    ref = refs[tid]
    f1 = 0; checks = 0; lines = 0; f2 = 0; worst = 0.0; errs = 0; c2 = 0
    for lo in range(0, N, 2500):
        inst = [generate(ref, s, capture=False) for s in range(lo, min(lo+2500, N))]
        errs += sum(1 for i in inst if not i.ok)
        r = t1_closure.run(inst, tid); f1 += len(r.failures)
        checks += r.checks_performed; lines += r.lines_with_eq
        r2 = T2.run(ref, inst, orc); f2 += len(r2.findings); c2 += r2.checked
        worst = max(worst, r2.worst_rel_error)
    small = [generate(ref, s) for s in range(300)]
    r4 = t4_contract.run(small, tid); r5 = t5_binding.run(ref, small); r7 = t7_asserts.run(ref)
    cov = checks and lines and (t1_closure.run([generate(ref,s,capture=False) for s in range(300)], tid).coverage)
    ok = lambda b: 'ok' if b else 'FAIL'
    print(f'{tid.replace("template_",""):42s} {f1:8d} {cov:6.2f} {f2:8d} {worst:10.2e} '
          f'{ok(r4.passed):>4s} {ok(r5.passed):>4s} {ok(r7.passed):>4s} {errs:4d}')
    if f1 or f2 or errs or not r4.passed or not r7.passed: bad += 1
print(f'\n{N} seeds each. T2 checked {c2} on the last template. Templates with a hard failure: {bad}')
