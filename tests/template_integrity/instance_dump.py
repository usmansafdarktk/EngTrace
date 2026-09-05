"""Dump (question, solution) pairs for the Phase 1 templates to JSON.

Run once per tree, each as its OWN process, then diff the two files - that is
how the item-pool impact note (D1.6) is produced:

    git worktree add /c/wt master
    (cd /c/wt && python -m tests.template_integrity.instance_dump before.json 1000)
    python -m tests.template_integrity.instance_dump after.json 1000

Separate processes are not optional. An in-process module reload reported "120
identical" for every template, because both sides resolved to the same already-
imported `data.templates.*` modules.

Not a check; gates nothing.
"""
import importlib
import json
import random
import sys

TEMPLATES = {
    'template_mean_variance':
        'data.templates.branches.electrical_engineering.digital_communications'
        '.deterministic_and_random_signal_analysis',
    'template_annulus_flowrate':
        'data.templates.branches.chemical_engineering.transport_phenomena'
        '.shell_momentum_balances',
    'template_rotating_unbalance':
        'data.templates.branches.mechanical_engineering.vibrations_and_acoustics'
        '.harmonically_excited_vibrations',
    'template_vibration_transmissibility':
        'data.templates.branches.mechanical_engineering.vibrations_and_acoustics'
        '.harmonically_excited_vibrations',
    'template_damping_classification':
        'data.templates.branches.mechanical_engineering.vibrations_and_acoustics'
        '.single_degree_of_freedom_systems',
    'template_statically_indeterminate_shaft':
        'data.templates.branches.mechanical_engineering.mechanics_of_materials.torsion',
    'template_shaft_design_power':
        'data.templates.branches.mechanical_engineering.mechanics_of_materials.torsion',
    'template_composite_shafts_series':
        'data.templates.branches.mechanical_engineering.mechanics_of_materials.torsion',
    'template_beam_deflection_formula':
        'data.templates.branches.civil_engineering.structural_analysis.deflections',
    'template_cantilever_double_integration':
        'data.templates.branches.civil_engineering.structural_analysis.deflections',
}

out_path = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 15
data = {}
for tid, mod in TEMPLATES.items():
    fn = getattr(importlib.import_module(mod), tid)
    rows = []
    for s in range(N):
        random.seed(s)
        rows.append(list(fn()))
    data[tid] = rows
json.dump(data, open(out_path, 'w', encoding='utf-8'))
print(f'wrote {out_path}: {len(data)} templates x {N} seeds')
