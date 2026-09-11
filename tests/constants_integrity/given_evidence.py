"""Evidence behind every `@given: stated` declaration the census relies on.

    python -m tests.constants_integrity.given_evidence            # the check
    python -m tests.constants_integrity.given_evidence --measure  # k/N only, no ratchet
    python -m tests.constants_integrity.given_evidence --selftest

WHY THIS EXISTS. The census classifies a `range` table as a plausibility window
only when every consumer is measured to restate it or reads it as a guard. When
a nudge makes a consumer RAISE or reword its question, the probe learns that the
table is read and not how, and the table is REVIEW until a human declares a
verdict with its evidence (`# @given: stated (...)`, spec §C1.3, SPEC-CHANGE 22;
C1 Reviewer G, G-2). A declaration is a claim, so this file is its check.

WHAT "STATED" MEANS HERE. For each consumer, the generator's frame locals are
captured at return, and each named DRAWN value must be printed in the QUESTION,
on every seed:
  * a number       - some number token in the question equals it exactly
                     (templates print a drawn value at the precision they draw it)
  * a list/tuple   - every element is printed
  * a string       - the string occurs in the question verbatim. Used where the
                     drawn value is printed only through a string it was
                     formatted into (a loop variable not retained at return), or
                     where the draw SHAPES a quantity the question states exactly
                     (a denominator reduced into the fraction the question prints)
A local named with a trailing `?` is optional: a branch that does not draw it is
skipped, never counted as a pass it did not earn.

THE RATCHET. Every `@given` in every constants.py must cite its evidence: this
file (and then have an entry here) or C1 Reviewer G's committed scripts. Every
entry here must back a table that declares `@given` citing this file.
"""
from __future__ import annotations

import os
import random
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, header_fields, numeric_tables)

SEEDS = 40
NUM_RE = re.compile(r'[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?')
CITES_THIS = 'given_evidence'
CITES_REVIEWER_G = 'phaseC1_reviewer_g'

#: (branch, table) -> [(template module, template id, [drawn locals])]
EVIDENCE = {
    ('electrical_engineering', 'FREQUENCY_RANGE_HZ'): [
        ('electrical_engineering.signals_and_systems.continuous_time_signals',
         'template_nyquist_rate_determination', ['frequencies_hz']),
    ],
    ('electrical_engineering', 'PHASE_RANGE_DEG'): [
        ('electrical_engineering.signals_and_systems.continuous_time_signals',
         'template_nyquist_rate_determination', ['signal_terms', 'phase_deg']),
        ('electrical_engineering.signals_and_systems.continuous_time_signals',
         'template_continuous_to_discrete_conversion', ['phi_deg?', 'signal_expression']),
    ],
    ('electrical_engineering', 'DECIMATION_FACTOR_M_RANGE'): [
        ('electrical_engineering.signals_and_systems.continuous_time_signals',
         'template_decimation_aliasing_analysis', ['M']),
    ],
    ('electrical_engineering', 'OMEGA_DENOMINATOR_RANGE'): [
        ('electrical_engineering.signals_and_systems.continuous_time_signals',
         'template_decimation_aliasing_analysis', ['omega_0_str']),
    ],
    ('industrial_engineering', 'HOLDING_RATE_PER_YR'): [
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_basic_eoq', ['i']),
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_epq_finite_production', ['i']),
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_quantity_discount_all_units', ['i']),
    ],
    ('industrial_engineering', 'INVENTORY_ITEMS'): [
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_basic_eoq', ['c', 'D', 'K']),
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_epq_finite_production', ['c', 'D', 'K']),
        ('industrial_engineering.production_and_inventory.deterministic_lot_sizing',
         'template_quantity_discount_all_units', ['D', 'K', 'schedule']),
    ],
}


def capture(fn, seed):
    """(question, the generator's frame locals at return), seeded as the census seeds."""
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:
        pass
    cap = {}
    code = fn.__code__

    def tracer(frame, event, arg):
        if frame.f_code is not code:
            return None

        def local(fr, ev, a):
            if ev == 'return':
                cap['locals'] = dict(fr.f_locals)
            return local
        return local

    old = sys.gettrace()
    sys.settrace(tracer)
    try:
        q, _s = fn()
    finally:
        sys.settrace(old)
    return str(q), cap.get('locals', {})


def printed(q, value):
    if isinstance(value, bool):
        return False
    if isinstance(value, str):
        return value in q
    if isinstance(value, (list, tuple)):
        return bool(value) and all(printed(q, v) for v in value)
    if isinstance(value, (int, float)):
        for tok in NUM_RE.findall(q):
            try:
                if float(tok.replace(',', '')) == float(value):
                    return True
            except ValueError:
                continue
        return False
    return False


def measure(fn, names, seeds=SEEDS):
    """name -> (seeds where printed, seeds where the local existed)."""
    out = {n.rstrip('?'): [0, 0] for n in names}
    for s in range(seeds):
        q, loc = capture(fn, s)
        for n in names:
            key = n.rstrip('?')
            if key not in loc:
                if not n.endswith('?'):
                    out[key][1] += 1        # required and absent counts as a miss
                continue
            out[key][1] += 1
            out[key][0] += printed(q, loc[key])
    return out


def declarations(sources):
    """{(branch, table): given text} for every @given in the given sources."""
    out = {}
    for branch, src in sources.items():
        for t in numeric_tables(src):
            g = header_fields(t['header']).get('given')
            if g:
                out[(branch, t['name'])] = g
    return out


def ratchet(decl, evidence):
    failures = []
    for key, text in decl.items():
        if CITES_THIS in text:
            if key not in evidence:
                failures.append(f'{key[0]}.{key[1]}: @given cites {CITES_THIS} but has no entry')
        elif CITES_REVIEWER_G not in text:
            failures.append(f'{key[0]}.{key[1]}: @given cites no committed evidence: {text[:80]!r}')
    for key in evidence:
        if CITES_THIS not in decl.get(key, ''):
            failures.append(f'{key[0]}.{key[1]}: has an evidence entry but no @given citing it')
    return failures


def run(check_ratchet=True):
    import importlib
    failures = []
    for (branch, table), consumers in EVIDENCE.items():
        for mod, tid, names in consumers:
            fn = getattr(importlib.import_module(f'data.templates.branches.{mod}'), tid)
            for name, (hit, seen) in measure(fn, names).items():
                ok = seen > 0 and hit == seen
                print(f"  [{'ok' if ok else 'FAIL'}] {branch[:5]}.{table:26s} {tid:42s} "
                      f"{name:18s} printed {hit}/{seen}")
                if not ok:
                    failures.append(f'{table} in {tid}: {name} printed {hit}/{seen}')
    if check_ratchet:
        sources = {b: open(os.path.join(BRANCHES_DIR, b, 'constants.py'), encoding='utf-8')
                   .read().replace('\r\n', '\n') for b in BRANCHES}
        decl = declarations(sources)
        rf = ratchet(decl, EVIDENCE)
        failures += rf
        print(f'  ratchet: {len(decl)} @given declarations, '
              f'{sum(CITES_THIS in v for v in decl.values())} citing this file, '
              f'{sum(CITES_REVIEWER_G in v for v in decl.values())} citing Reviewer G')
    for f in failures:
        print('  - ' + f)
    print('all pass' if not failures else f'{len(failures)} FAILURES')
    return 1 if failures else 0


# --------------------------------------------------------------------------
# Self-test: plants written from what "stated" means above.
# --------------------------------------------------------------------------

def _t_hidden_number():
    x = round(random.uniform(1, 2), 2)
    return "Find twice a hidden value.", f"**Answer:** {2 * x}"


def _t_list_partly_printed():
    xs = [random.randint(10, 99), random.randint(100, 999)]
    return f"The first value is {xs[0]}.", f"**Answer:** {sum(xs)}"


def _t_string_not_stated():
    expr = f"cos({random.randint(2, 9)}*pi*n)"
    return "A signal is given.", f"**Answer:** {expr}"


def _t_optional_absent_and_stated():
    v = random.randint(1, 9)
    if random.random() < 0.5:
        w = random.randint(10, 19)
        return f"v = {v}, w = {w}", "**Answer:** 0"
    return f"v = {v}", "**Answer:** 0"


def selftest():
    bad = []
    plants = [
        ('a drawn number the question never prints', _t_hidden_number, ['x'], False),
        ('a list with one element unprinted', _t_list_partly_printed, ['xs'], False),
        ('a shaped string the question never states', _t_string_not_stated, ['expr'], False),
        ('optional local absent on some seeds, printed when drawn (control)',
         _t_optional_absent_and_stated, ['v', 'w?'], True),
        ('a REQUIRED local absent on some seeds', _t_optional_absent_and_stated, ['w'], False),
    ]
    for label, fn, names, want_pass in plants:
        res = measure(fn, names, seeds=20)
        passed = all(seen > 0 and hit == seen for hit, seen in res.values())
        ok = passed == want_pass
        print(f"  [{'ok' if ok else 'FAIL'}] {label}: {res}")
        if not ok:
            bad.append(label)
    src_ok = '# @kind: range\n# @given: stated (given_evidence.py: planted)\nA = (1, 2)\n'
    src_g = '# @kind: range\n# @given: stated (phaseC1_reviewer_g row 1)\nB = (1, 2)\n'
    src_none = '# @kind: range\n# @given: stated (trust me)\nC = (1, 2)\n'
    cases = [
        ('declaration citing this file with an entry', {'b': src_ok}, {('b', 'A'): []}, 0),
        ('declaration citing this file without an entry', {'b': src_ok}, {}, 1),
        ('declaration citing Reviewer G', {'b': src_g}, {}, 0),
        ('declaration citing nothing committed', {'b': src_none}, {}, 1),
        ('entry backing no declaration', {'b': src_g}, {('b', 'Z'): []}, 1),
    ]
    for label, sources, ev, want in cases:
        got = len(ratchet(declarations(sources), ev))
        ok = got == want
        print(f"  [{'ok' if ok else 'FAIL'}] ratchet: {label} -> {got} failure(s)")
        if not ok:
            bad.append(label)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    if '--selftest' in sys.argv:
        sys.exit(selftest())
    sys.exit(run(check_ratchet='--measure' not in sys.argv))
