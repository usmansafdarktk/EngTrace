"""CI entry point for the template integrity suite (D0.4).

    python -m tests.template_integrity.run                     # T1,T2,T4,T5,T7,T8 corpus-wide
    python -m tests.template_integrity.run --checks all        # add T3 and T6
    python -m tests.template_integrity.run --templates a,b     # only these
    python -m tests.template_integrity.run --branch civil_engineering
    python -m tests.template_integrity.run --baseline          # (re)write the baseline
    python -m tests.template_integrity.run --json out.json     # machine-readable

Exit code is non-zero if any selected check fails on any selected template,
so it can gate a merge directly.

T3 (determinism) runs the whole selection in two child processes total (~8s for
150 templates) but is still opt-in because it needs a clean environment.
T2 (round-trip) applies only to templates that have an oracle module.
T6 (distribution) needs a committed baseline to compare against.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

from .core import REPO_ROOT, discover, generate
from .checks import (t1_closure, t2_roundtrip, t3_determinism, t4_contract,
                     t5_binding, t6_distribution, t7_asserts, t8_emission)

BASELINE_DIR = os.path.join(os.path.dirname(__file__), 'baseline')
PROFILE_PATH = os.path.join(BASELINE_DIR, 'profiles.json')

DEFAULT_CHECKS = ('T1', 'T2', 'T4', 'T5', 'T7', 'T8')
ALL_CHECKS = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8')


def _select(args):
    refs = discover(branch=args.branch)
    if args.templates:
        wanted = {t.strip() for t in args.templates.split(',') if t.strip()}
        refs = [r for r in refs if r.template_id in wanted]
        missing = wanted - {r.template_id for r in refs}
        if missing:
            sys.exit(f'unknown template(s): {sorted(missing)}')
    return refs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='EngTrace template integrity suite')
    ap.add_argument('--checks', default=','.join(DEFAULT_CHECKS),
                    help=f'comma-separated subset of {ALL_CHECKS}, or "all"')
    ap.add_argument('--templates', default='', help='comma-separated template ids')
    ap.add_argument('--branch', default=None)
    ap.add_argument('--seeds', type=int, default=25,
                    help='instances per template for the output-based checks')
    ap.add_argument('--det-seeds', type=int, default=200, help='seeds for T3')
    ap.add_argument('--t8-seeds', type=int, default=400,
                    help='instances per template for T8 (see t8_emission)')
    ap.add_argument('--baseline', action='store_true',
                    help='write the T6 baseline profile instead of comparing')
    ap.add_argument('--json', default='', help='write a machine-readable report')
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--strict-marginals', type=float, default=None, metavar='RATE',
                    help='fail a template whose T1 marginal rate exceeds RATE '
                         '(e.g. 0.25). A marginal is within tolerance by '
                         'definition, so it is never a failure by default - but '
                         'a high density means the template sits repeatedly on '
                         'its own display-rounding boundary, which is how the '
                         'half-way-tie defect hides.')
    args = ap.parse_args(argv)

    checks = (list(ALL_CHECKS) if args.checks == 'all'
              else [c.strip().upper() for c in args.checks.split(',') if c.strip()])
    refs = _select(args)
    t0 = time.time()
    report: dict[str, dict] = {}
    failed: set[str] = set()

    baseline = {}
    if 'T6' in checks and not args.baseline and os.path.exists(PROFILE_PATH):
        baseline = t6_distribution.load(PROFILE_PATH)

    # T3 runs ONCE for the whole selection (two child processes total), not
    # once per template - see t3_determinism for why that is equivalent.
    det = {}
    if 'T3' in checks:
        det = t3_determinism.run_many(refs, range(args.det_seeds))

    oracles = t2_roundtrip.available_oracles() if 'T2' in checks else {}

    profiles = {}
    for ref in refs:
        entry: dict[str, object] = {'branch': ref.branch, 'file': ref.file_path}
        need_output = {'T1', 'T2', 'T4', 'T6'} & set(checks)
        need_values = 'T5' in checks
        instances = []
        if need_output or need_values:
            instances = [generate(ref, s, capture=need_values)
                         for s in range(args.seeds)]
            errs = [i for i in instances if not i.ok]
            if errs:
                entry['generation_errors'] = [i.error for i in errs[:3]]
                failed.add(ref.template_id)

        if 'T1' in checks:
            r = t1_closure.run(instances, ref.template_id)
            entry['T1'] = {'pass': r.passed, 'failures': len(r.failures),
                           'failing_seeds': sorted({f.seed for f in r.failures}),
                           'marginals': len(r.marginals),
                           'checks_performed': r.checks_performed,
                           'marginal_rate': round(
                               len(r.marginals) / r.checks_performed, 3)
                           if r.checks_performed else 0.0,
                           'unit_rescales': r.unit_rescales,
                           'coverage': round(r.coverage, 3),
                           'examples': [str(f) for f in r.failures[:2]]}
            if not r.passed:
                failed.add(ref.template_id)
            if (args.strict_marginals is not None and r.checks_performed
                    and len(r.marginals) / r.checks_performed
                    > args.strict_marginals):
                entry['T1']['pass'] = False
                entry['T1']['marginal_breach'] = True
                failed.add(ref.template_id)

        if 'T2' in checks:
            r2 = t2_roundtrip.run(ref, instances, oracles)
            if r2.note == 'no oracle':
                # Reported, never silently passed: a template with no oracle is
                # unverified by T2, which is a gap, not a pass.
                entry['T2'] = {'pass': True, 'note': 'no oracle', 'checked': 0}
            else:
                entry['T2'] = {'pass': r2.passed, 'oracle': r2.oracle,
                               'checked': r2.checked, 'unparsed': r2.unparsed,
                               'failures': len(r2.findings),
                               'failure_rate': round(r2.failure_rate, 4),
                               'worst_rel_error': r2.worst_rel_error,
                               'note': r2.note,
                               'examples': [str(f) for f in r2.findings[:2]]}
                if not r2.passed:
                    failed.add(ref.template_id)

        if 'T3' in checks:
            r3 = det[ref.template_id]
            entry['T3'] = {'pass': r3.passed, 'mismatched': r3.mismatched_seeds[:10],
                           'errored': r3.errored_seeds[:10], 'note': r3.note}
            if not r3.passed:
                failed.add(ref.template_id)

        if 'T4' in checks:
            r4 = t4_contract.run(instances, ref.template_id)
            entry['T4'] = {'pass': r4.passed, 'summary': r4.summary(),
                           'step_counts': dict(r4.step_counts)}
            if not r4.passed:
                failed.add(ref.template_id)

        if 'T5' in checks:
            r5 = t5_binding.run(ref, instances)
            entry['T5'] = {'pass': r5.passed,
                           'inline_results': [str(f) for f in r5.inline_results],
                           'rounding_violations': [str(f) for f in r5.rounding_violations],
                           'operand_restatements': r5.operand_restatements}
            if not r5.passed:
                failed.add(ref.template_id)

        if 'T6' in checks:
            p = t6_distribution.profile(instances, ref.template_id)
            profiles[ref.template_id] = p
            if not args.baseline:
                base = baseline.get(ref.template_id)
                if base is None:
                    entry['T6'] = {'pass': True, 'note': 'no baseline recorded'}
                else:
                    d = t6_distribution.compare(base, p)
                    entry['T6'] = {'pass': d.passed, 'breaches': d.breaches,
                                   'notes': d.notes}
                    if not d.passed:
                        failed.add(ref.template_id)

        if 'T8' in checks:
            # T8 generates its OWN instances, at its own seed count: the
            # degenerate-product class fires on 1.95% of one template's
            # instances and 25 seeds resolves nothing rarer than ~12%
            # (D-024/D-026).  It also reads `inst.question`, which nothing else
            # in the suite does (Reviewer A, F2).
            r8 = t8_emission.run(ref, seeds=args.t8_seeds)
            entry['T8'] = {'pass': r8.passed, 'summary': r8.summary(),
                           'instances': r8.instances,
                           'gated_classes': list(r8.gated)}
            if not r8.passed:
                failed.add(ref.template_id)

        if 'T7' in checks:
            r7 = t7_asserts.run(ref)
            entry['T7'] = {'pass': r7.passed, 'summary': r7.summary(),
                           'output_asserts': len(r7.output_asserts)}
            # T7 is a gate on edited templates only; corpus-wide it is advisory.
            if not r7.passed and args.templates:
                failed.add(ref.template_id)

        report[ref.template_id] = entry
        if not args.quiet:
            marks = ''.join(
                ('.' if entry.get(c, {}).get('pass', True) else 'X')
                for c in checks if c in entry)
            print(f'{marks:8s} {ref.branch[:4]} {ref.template_id}')

    if args.baseline:
        os.makedirs(BASELINE_DIR, exist_ok=True)
        t6_distribution.save(profiles, PROFILE_PATH)
        print(f'\nbaseline written: {PROFILE_PATH} ({len(profiles)} templates)')

    if args.json:
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=1)

    print(f'\n{len(refs)} templates, checks {checks}, '
          f'{time.time() - t0:.0f}s')
    print(f'failing: {len(failed)}')
    for c in checks:
        n = sum(1 for e in report.values()
                if c in e and not e[c].get('pass', True))
        print(f'  {c}: {n} failing')
    if 'T1' in checks:
        marg = [(t, e['T1']['marginal_rate']) for t, e in report.items()
                if e.get('T1', {}).get('marginals')]
        if marg:
            marg.sort(key=lambda x: -x[1])
            print(f'  T1 marginals: {len(marg)} templates carry marginals '
                  f'(within tolerance, but on the rounding boundary). '
                  f'Worst: ' + ', '.join(f'{t.replace("template_", "")} {r:.0%}'
                                          for t, r in marg[:3]))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
