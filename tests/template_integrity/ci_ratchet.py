"""D6.5 - a ratchet over the integrity suite, because there is no green suite to gate on.

The spec says "promote the regression suite to CI". Measured, the suite does not
pass: 78 of 150 templates fail the default battery, and `run.py` exits 1.

    T1  29 failing      T5  66 failing      T6  142 failing
    T7  83 failing      T2/T3/T4/T8  0

None of that is a regression. It is pre-existing debt carried across six phases,
and T6's 142 is a corpus-wide STALE BASELINE (D-043), not 142 defects.

So a plain gate is the wrong instrument twice over. It would be red on every
commit, which makes it useless; and a permanently red pipeline creates exactly
one pressure - to regenerate the T6 baseline until it goes green. That is the
single action D-043 forbids, because *a baseline refreshed by the phase it is
meant to gate is not a gate*.

A ratchet gates on the DELTA instead. Each check has a recorded ceiling; the run
fails only if a check does WORSE than its ceiling. Existing debt does not block
anyone, and a genuine regression does.

Two rules keep it from rotting into a rubber stamp:

* the ceilings are MEASURED into `ci_expected.json` from a real run
  (`--write-expected`), never hand-typed;
* it only ever goes DOWN. An improvement is reported loudly as a ceiling that
  should be tightened, because a ratchet nobody tightens is just a gate with a
  high bar.

A check that appears in the report but not in the ceilings file is a FAILURE,
not a default-pass: a newly added check must be declared deliberately, or it
would ride along unnoticed at whatever it happens to score.

Usage::

    python -m tests.template_integrity.run --checks all --json report.json
    python -m tests.template_integrity.ci_ratchet --report report.json

    # after a deliberate improvement, re-record:
    python -m tests.template_integrity.ci_ratchet --report report.json --write-expected
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXPECTED_PATH = os.path.join(HERE, 'baseline', 'ci_expected.json')

#: Checks whose failures this file is allowed to describe.  Kept explicit so a
#: typo in a report key cannot silently create a new, unratcheted check.
KNOWN = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8')


def failing_counts(report: dict) -> dict[str, int]:
    """{check: number of TEMPLATES failing it}, matching run.py's own summary.

    Counts templates rather than individual failures because that is what
    `run.py` prints, and a ratchet whose units differ from the summary beside it
    is a ratchet nobody can check by eye.
    """
    counts: dict[str, int] = {}
    for _tid, row in report.items():
        if not isinstance(row, dict):
            continue
        for check, res in row.items():
            if check not in KNOWN or not isinstance(res, dict) or 'pass' not in res:
                continue
            counts.setdefault(check, 0)
            if not res['pass']:
                counts[check] += 1
    return counts


def compare(counts: dict[str, int], expected: dict[str, int]) -> tuple[int, list[str]]:
    """(exit_code, lines). Worse than the ceiling fails; better is reported."""
    lines, worse, better, undeclared = [], [], [], []

    for check in sorted(counts):
        got = counts[check]
        if check not in expected:
            undeclared.append(check)
            lines.append(f'  {check}  {got:4d}  UNDECLARED - add it to ci_expected.json')
            continue
        ceiling = expected[check]
        if got > ceiling:
            worse.append(check)
            lines.append(f'  {check}  {got:4d}  REGRESSION - ceiling is {ceiling}')
        elif got < ceiling:
            better.append(check)
            lines.append(f'  {check}  {got:4d}  improved from {ceiling} - TIGHTEN the ceiling')
        else:
            lines.append(f'  {check}  {got:4d}  at ceiling')

    for check in sorted(set(expected) - set(counts)):
        lines.append(f'  {check}     -  declared but not in this report (not selected?)')

    if worse:
        lines.append('')
        lines.append(f'FAIL: {len(worse)} check(s) regressed: {", ".join(worse)}')
    if undeclared:
        lines.append(f'FAIL: {len(undeclared)} undeclared check(s): {", ".join(undeclared)}')
    if better and not worse and not undeclared:
        lines.append('')
        lines.append(f'PASS, but {len(better)} check(s) improved - re-record with '
                     f'--write-expected so the gain cannot be lost silently.')
    if not worse and not undeclared and not better:
        lines.append('')
        lines.append('PASS: every check at its recorded ceiling.')

    return (1 if (worse or undeclared) else 0), lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--report', required=True,
                    help='JSON written by run.py --json')
    ap.add_argument('--expected', default=EXPECTED_PATH)
    ap.add_argument('--write-expected', action='store_true',
                    help='record the current counts as the ceilings')
    args = ap.parse_args(argv)

    if not os.path.exists(args.report):
        print(f'no report at {args.report} - run run.py --json first')
        return 1

    with open(args.report, encoding='utf-8') as fh:
        counts = failing_counts(json.load(fh))

    if args.write_expected:
        os.makedirs(os.path.dirname(args.expected), exist_ok=True)
        with open(args.expected, 'w', encoding='utf-8') as fh:
            json.dump(counts, fh, indent=2, sort_keys=True)
            fh.write('\n')
        print(f'recorded ceilings from {args.report}:')
        for check in sorted(counts):
            print(f'  {check}  {counts[check]}')
        print(f'-> {args.expected}')
        return 0

    if not os.path.exists(args.expected):
        print(f'no ceilings at {args.expected} - create them with --write-expected')
        return 1

    with open(args.expected, encoding='utf-8') as fh:
        expected = json.load(fh)

    code, lines = compare(counts, expected)
    print('integrity ratchet (templates failing, per check)')
    for line in lines:
        print(line)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
