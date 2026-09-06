"""C2.5 — every citation resolves to the artefact it names.

The C2.2 plausibility suite asks *"is this value right?"*. This one asks the
prior question: *"does the thing the comment points at actually exist, and is
it what the comment says it is?"* Phase C2 shipped a green 215-check suite that
nonetheless certified three heats of formation with no on-disk backing at all,
because the reference values lived in the test file rather than in the artefact
the citations named (Reviewer G, findings G-1, G-2, G-12; DECISIONS D-034).

Four properties, each of which failed somewhere in Phase C2:

  P1  every value-bearing row carries a provenance tag
  P2  the CAS NUMBER WRITTEN IN THE TAG resolves - not merely the species key.
      A row keyed `CH3OH(g)` whose comment cites the wrong CAS was previously
      undetectable, because the old check resolved by key and the key was fine.
  P3  the tag's claim class matches the evidence: [ON-DISK] requires an entry
      in the artefact, [DERIVED] requires a stated derivation
  P4  no reference value is hardcoded in a test file. This is the regression
      guard for the root cause: a suite carrying its own answer key cannot
      fail when the answers drift.

Run:
    python -m tests.constants_integrity.test_citations_resolve
"""
from __future__ import annotations

import json
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

CONSTANTS = os.path.join(
    REPO, 'data', 'templates', 'branches', 'chemical_engineering',
    'constants.py')
REF_PATH = os.path.join(
    REPO, 'docs', 'references', 'nist_webbook', 'shomate_coefficients.json')

# Tables whose rows must each carry a citation, and the artefact they cite.
CITED_TABLES = ('CP_PARAMS', 'HEATS_OF_FORMATION')

TAG_RE = re.compile(r'#\s*\[(ON-DISK|DERIVED|KNOWN-DEFECTIVE|BY-DEFINITION)\]\s*(.*)')
CAS_RE = re.compile(r'\b(\d{2,7}-\d{2}-\d)\b')
ROW_RE = re.compile(r'^\s*"([^"]+)"\s*:')


def parse_table(src, name):
    """Yield (species, tag_class, tag_text) for each row of a dict literal.

    A row's citation is the tag block immediately above it - the same thing a
    human reads when they scan the table.
    """
    lines = src.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(name + ' = {'))
    depth, pending = 0, []
    for line in lines[start:]:
        depth += line.count('{') - line.count('}')
        m = TAG_RE.search(line)
        if m:
            pending.append((m.group(1), m.group(2).strip()))
            continue
        if line.lstrip().startswith('#'):
            if pending:                       # continuation of the tag block
                pending[-1] = (pending[-1][0],
                               pending[-1][1] + ' ' + line.lstrip('# ').strip())
            continue
        row = ROW_RE.match(line)
        if row:
            yield row.group(1), (pending[-1] if pending else (None, ''))
            pending = []
        if depth == 0 and line.startswith('}'):
            break


def run():
    src = open(CONSTANTS, encoding='utf-8').read()
    ref = json.load(open(REF_PATH, encoding='utf-8'))
    by_cas = {}
    for sp, e in ref.items():
        if not sp.startswith('_') and e.get('cas'):
            by_cas.setdefault(e['cas'], []).append(sp)

    failures, checks = [], 0

    for table in CITED_TABLES:
        for sp, (cls, text) in parse_table(src, table):
            checks += 1

            # P1 - the row is cited at all.
            if cls is None:
                failures.append(f'{table}[{sp}] carries no provenance tag')
                continue

            if cls == 'DERIVED':
                # P3 - a derivation has to say what it was derived from.
                if not text:
                    failures.append(f'{table}[{sp}] is [DERIVED] but states '
                                    f'no derivation')
                continue

            if cls == 'KNOWN-DEFECTIVE':
                continue                       # not a claim of correctness

            if cls == 'BY-DEFINITION':
                # An element in its standard state is 0 because of how the
                # scale is defined, not because a document says so. Tagging it
                # [ON-DISK] would point at an artefact that cannot be the
                # warrant for the claim.
                if not text:
                    failures.append(f'{table}[{sp}] is [BY-DEFINITION] but '
                                    f'states no definition')
                continue

            # P3 - [ON-DISK] means there is an entry on disk.
            entry = ref.get(sp)
            if entry is None:
                failures.append(
                    f'{table}[{sp}] is tagged [ON-DISK] but no entry for it '
                    f'exists in {os.path.basename(REF_PATH)}')
                continue

            # P2 - the CAS in the tag, not just the key, has to be right.
            cas = CAS_RE.search(text)
            if not cas:
                failures.append(f'{table}[{sp}] is [ON-DISK] but its citation '
                                f'names no CAS number to resolve')
                continue
            cas = cas.group(1)
            if cas not in by_cas:
                failures.append(
                    f'{table}[{sp}] cites CAS {cas}, which appears nowhere in '
                    f'{os.path.basename(REF_PATH)}')
            elif entry.get('cas') != cas:
                failures.append(
                    f'{table}[{sp}] cites CAS {cas}, but the artefact entry '
                    f'for {sp} is CAS {entry.get("cas")} '
                    f'(that CAS is {", ".join(by_cas[cas])})')

    # P4 - no test may carry its own answer key.
    checks += 1
    embedded = []
    test_dir = os.path.dirname(os.path.abspath(__file__))
    for fn in sorted(os.listdir(test_dir)):
        if not fn.startswith('test_') or not fn.endswith('.py'):
            continue
        if fn == os.path.basename(__file__):
            continue
        body = open(os.path.join(test_dir, fn), encoding='utf-8').read()
        for m in re.finditer(r'^([A-Z][A-Z0-9_]*(?:REF|_VALUES|_TABLE))\s*=\s*\{',
                             body, re.M):
            embedded.append(f'{fn}:{m.group(1)}')
    if embedded:
        failures.append(
            'reference values are hardcoded in a test file, so the suite '
            'checks itself rather than the artefact: ' + ', '.join(embedded))

    for f in failures:
        print('  - ' + f)
    print(f'{checks} checks')
    print('all pass' if not failures else f'{len(failures)} FAILURES')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(run())
