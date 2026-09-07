"""Dump the eleven Track A templates' instances to JSON, for the D5.1 before/after.

Run ONCE PER TREE, each as its OWN process, then diff the two files:

    git worktree add --detach C:/wtm master
    python -m tests.template_integrity.phase5_instance_dump C:/wtm before.json 2000
    python -m tests.template_integrity.phase5_instance_dump .     after.json  2000
    python -m tests.template_integrity.phase5_instance_dump --diff before.json after.json

Separate processes are not optional, and neither is the explicit tree root. An
in-process reload resolves BOTH sides to the same already-imported
`data.templates.*` modules and reports every instance identical - the trap
documented in `instance_dump.py`, which has now caught four people. This script
takes the tree root as argv[1] and puts it at the FRONT of sys.path precisely so
that the mistake is visible in the command line rather than hidden in an import.

`master` does not carry this file, so run BOTH sides with the copy from the
Phase 5 branch, passing the master worktree as the root.

**Why this exists rather than a T6 delta.** T6's baseline is stale corpus-wide
(D-043), regenerating it to make a gate green is forbidden, and its answer
extractor cannot see a non-scalar answer at all (R4-5). It is therefore not
evidence about the item pool. This dump is: it compares the emitted text of
every instance, on both sides, at the two trees.

The `--diff` mode classifies each changed instance by **which of Track A's four
defect classes** the change belongs to, so that a formatting-only change can be
told apart from a change to what the item asks or answers. The distinction that
matters for P6 is not "did the text change" - all eleven were edited, so of
course it did - but **"did the question change, and did the answer change"**.
Those two counts are the item-pool impact statement, and they are printed last.

Not a check; gates nothing.
"""
import difflib
import importlib
import json
import os
import random
import re
import sys

TEMPLATES = {
    # class 1 - malformed **Step N:**
    'template_cd_dc_system_analysis':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.continuous_time_signals',
    'template_euclidean_distance_binary':
        'data.templates.branches.electrical_engineering.digital_communications'
        '.digital_modulation_schemes',
    'template_finite_convolution':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.discrete_time_signals',
    # class 2 - non-canonical answer marker
    'template_batch_moles_vs_conversion':
        'data.templates.branches.chemical_engineering.reaction_kinetics.stoichiometry',
    'template_flow_system_molar_flow_rates':
        'data.templates.branches.chemical_engineering.reaction_kinetics.stoichiometry',
    'template_gas_phase_concentration':
        'data.templates.branches.chemical_engineering.reaction_kinetics.stoichiometry',
    'template_limiting_reactant':
        'data.templates.branches.chemical_engineering.reaction_kinetics.stoichiometry',
    'template_levenspiel_plot_interpretation':
        'data.templates.branches.chemical_engineering.reaction_kinetics'
        '.conversion_and_reactor_sizing',
    # class 3 - malformed complex sign
    'template_time_to_phasor':
        'data.templates.branches.electrical_engineering.electromagnetics_and_waves'
        '.waves_and_phasors',
    'template_phasor_addition':
        'data.templates.branches.electrical_engineering.electromagnetics_and_waves'
        '.waves_and_phasors',
    # class 4 - degenerate product from an unreachable guard
    'template_decimation_aliasing_analysis':
        'data.templates.branches.electrical_engineering.signals_and_systems'
        '.continuous_time_signals',
}

#: Answer markers, longest first.  Both the old and the new shapes, because the
#: `before` side emits markers the `after` side no longer does -- and comparing
#: answers across the edit is the whole point.
_MARKERS = ('**Final Answers:**', '**Final Answer**', '**Answer:**', '**Answer**')


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


def _answer(sol):
    cut = max((sol.rfind(mk) for mk in _MARKERS), default=-1)
    return sol[cut:] if cut >= 0 else sol


def _answer_body(sol):
    """The answer span with its marker removed, so a marker rename is not a
    change to the answer.  This is the distinction the D5.3 edit turns on: the
    five chemical templates changed their answer *marker* and nothing else, and
    a comparison that cannot see that would report five items as re-answered."""
    span = _answer(sol)
    for mk in _MARKERS:
        if span.startswith(mk):
            return span[len(mk):].strip()
    return span.strip()


def _changed_lines(b, a):
    out = []
    for ln in difflib.unified_diff(b.split('\n'), a.split('\n'), n=0, lineterm=''):
        if ln.startswith(('+++', '---', '@@')):
            continue
        if ln.startswith(('+', '-')):
            out.append(ln.rstrip())
    return out


#: How a changed line is attributed to one of the four defect classes.  Applied
#: to the pair of changed lines, not to the template id, so a template that
#: changed for two reasons is counted under both.
_CLASSIFIERS = (
    ('marker-step', re.compile(r'\*\*\s*Step')),
    ('marker-answer', re.compile(r'\*\*\s*(?:Final\s+)?Answers?\s*:?\s*\*\*')),
    ('complex-sign', re.compile(r'\bj-?\d|\bj\(|<\s*-?\d+(?:\.\d+)?\s*deg')),
    ('signed-term', re.compile(r'[-+]\s+[-+]\s*\d|\bcos\(|\bsin\(|sqrt\(')),
    ('degenerate-product', re.compile(r'(?<![\w.])0\s*\*\s*pi|omega_a = ')),
)


def _classify(changed):
    tags = set()
    for ln in changed:
        for name, rx in _CLASSIFIERS:
            if rx.search(ln):
                tags.add(name)
    return tags or {'unclassified'}


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

    grand = {'instances': 0, 'solutions_changed': 0,
             'questions_changed': 0, 'answer_bodies_changed': 0,
             'answer_spans_changed': 0, 'errors_before': 0, 'errors_after': 0}
    tag_totals = {}
    for tid in TEMPLATES:
        b, a = before['instances'][tid], after['instances'][tid]
        q_ch = s_ch = span_ch = body_ch = 0
        tags = {}
        examples = []
        for i, (bi, ai) in enumerate(zip(b, a)):
            grand['instances'] += 1
            grand['errors_before'] += bi[1].startswith('ERROR ')
            grand['errors_after'] += ai[1].startswith('ERROR ')
            if bi[0] != ai[0]:
                q_ch += 1
            if bi[1] == ai[1]:
                continue
            s_ch += 1
            if _answer(bi[1]) != _answer(ai[1]):
                span_ch += 1
            if _answer_body(bi[1]) != _answer_body(ai[1]):
                body_ch += 1
            changed = _changed_lines(bi[1], ai[1])
            for t in _classify(changed):
                tags[t] = tags.get(t, 0) + 1
                tag_totals[t] = tag_totals.get(t, 0) + 1
            if len(examples) < 1:
                examples.append((i, changed))
        grand['questions_changed'] += q_ch
        grand['solutions_changed'] += s_ch
        grand['answer_spans_changed'] += span_ch
        grand['answer_bodies_changed'] += body_ch

        print(f"{tid}")
        print(f"    solutions changed     : {s_ch}/{len(b)}")
        print(f"    questions changed     : {q_ch}/{len(b)}")
        print(f"    answer span changed   : {span_ch}/{len(b)}")
        print(f"    answer BODY changed   : {body_ch}/{len(b)}   <- the P6 number")
        if tags:
            print(f"    change classes        : "
                  f"{', '.join(f'{k}={v}' for k, v in sorted(tags.items()))}")
        for i, changed in examples:
            print(f"    seed {i}:")
            for ln in changed[:4]:
                print(f"      {ln[:104]}")
        print()

    print("SUMMARY")
    for k, v in grand.items():
        print(f"  {k:24s} {v}")
    print(f"  change classes across all: "
          f"{', '.join(f'{k}={v}' for k, v in sorted(tag_totals.items()))}")
    print()
    if grand['errors_after'] > grand['errors_before']:
        print("!! MORE GENERATION ERRORS AFTER THAN BEFORE - this is a regression.")
    if grand['questions_changed'] == 0 and grand['answer_bodies_changed'] == 0:
        print("ITEM POOL UNCHANGED IN SUBSTANCE: no question and no answer body")
        print("differs on any seed. Every difference is a marker or a sign position.")
    else:
        print("ITEM POOL CHANGED - this is a P6 event and needs sign-off.")
        print("Read the per-template `questions changed` and `answer BODY changed`")
        print("counts above; `answer span changed` alone is only a marker rename.")


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == '--diff':
        diff(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 3:
        dump(sys.argv)
    else:
        raise SystemExit(__doc__)
