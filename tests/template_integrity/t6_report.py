"""T6 REPORTING utility - dump distribution profiles to JSON at a given seed count.

This is not a check and gates nothing; `checks/t6_distribution.py` owns the
gate. It exists so a phase's before/after distribution diff (D1.5) is
reproducible:

    git worktree add /c/wt master
    (cd /c/wt && python -m tests.template_integrity.t6_report before.json 5000)
    python -m tests.template_integrity.t6_report after.json 5000


Run once in each tree (master worktree = before, working tree = after) and
compare with _t6merge.py. The committed baseline is only 1,000 seeds, and a
distinct-answer count SATURATES there - `shaft_design_power` lost 29% of its
answer space and still read -10.0% at N=1,000. Both sides must be measured at
the same, larger N.
"""
import json
import sys

from tests.template_integrity.core import discover, generate
import tests.template_integrity.checks.t6_distribution as T6

TEMPLATES = [
    'template_mean_variance', 'template_rotating_unbalance',
    'template_vibration_transmissibility', 'template_damping_classification',
    'template_beam_deflection_formula',
    'template_cantilever_double_integration', 'template_annulus_flowrate',
    'template_statically_indeterminate_shaft', 'template_shaft_design_power',
    'template_composite_shafts_series',
]

out_path = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
refs = {r.template_id: r for r in discover()}
profiles = {}
for tid in TEMPLATES:
    inst = [generate(refs[tid], s, capture=False) for s in range(N)]
    profiles[tid] = T6.profile(inst, tid).to_json()
    print('.', end='', flush=True)
json.dump({'n': N, 'profiles': profiles}, open(out_path, 'w'), indent=1)
print(f'\nwrote {out_path} at N={N}')
