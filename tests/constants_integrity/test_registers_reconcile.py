"""The residual registers say what is open - this checks they still do.

    python -m tests.constants_integrity.test_registers_reconcile
    python -m tests.constants_integrity.test_registers_reconcile --selftest

A register is read at the START of the next phase, and by then the tree has moved
under it. That is not hypothetical here: at 2026-09-16 both registers stated
counts that seven owner-directed row deletions had made stale - [UNVERIFIED] 177
against a measured 170, [KNOWN-DEFECTIVE] 19 against 0, two domain excursions
against one, three sign helpers against one. Every figure was right when written.
None was right when read.

The project's own rule (README, "Conventions worth knowing") is that every count
ships with its predicate or a committed script. These registers publish `grep`
recipes; this is the script. It does NOT hold its own copy of the numbers - it
parses the figures the documents state and re-derives each from the tree, so a
drift on either side fails, and the failure names which side moved.

  R1  every [UNVERIFIED] / [KNOWN-DEFECTIVE] count either register states, per
      branch and in total, equals the count in data/templates/branches
  R2  every substance row the registers say was DELETED is absent from the
      branch files, and the row kept in its place is present
  R3  the C3 reason table sums to the total it declares, in both its columns
  R4  the consumer-domain excursion count matches domain_findings.txt
  R5  the sign-helper section's count matches the helpers surviving in the tree
  R6  the binding figures - bound, unbound, symbolic, and the per-template
      false-accept counts - match tests.comparators.bindings
  R7  D-050's "N of M seeds silent" reproduces at the M the document states

Each number is checked by the instrument that owns it. Tag totals, resolved
counts and LEGACY are test_citations_resolve's and are deliberately not
re-derived here.
"""
from __future__ import annotations

import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

BRANCHES = ('mechanical', 'chemical', 'civil', 'electrical', 'industrial')

PHASE6 = os.path.join(_ROOT, 'docs', 're-implementation-sep', 'track-a',
                      'phase6_residual_register.md')
C3 = os.path.join(_ROOT, 'docs', 're-implementation-sep', 'track-b',
                  'phaseC3_residual_register.md')

_ORDINAL = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five', 6: 'Six', 7: 'Seven'}
_EMDASH = '—'


# ---------------------------------------------------------------- measurement

def _branch_file(branch):
    return os.path.join(_ROOT, 'data', 'templates', 'branches',
                        branch + '_engineering', 'constants.py')


def _read(path):
    with open(path, encoding='utf-8') as fh:
        return fh.read()


def measure():
    """Every fact the registers assert, taken from the tree. Measured once."""
    facts = {'unverified': {}, 'known_defective': {}, 'keys': set()}
    for br in BRANCHES:
        src = _read(_branch_file(br))
        lines = src.splitlines()
        facts['unverified'][br] = sum(
            1 for ln in lines if re.match(r'\s*#\s*\[UNVERIFIED\]', ln))
        facts['known_defective'][br] = sum(
            1 for ln in lines if re.match(r'\s*#\s*\[KNOWN-DEFECTIVE\]', ln))
        facts['keys'] |= set(re.findall(r'^\s*[\'"]([^\'"]+)[\'"]\s*:', src, re.M))
    facts['unverified']['total'] = sum(facts['unverified'].values())
    facts['known_defective']['total'] = sum(facts['known_defective'].values())

    findings = _read(os.path.join(_HERE, 'domain_findings.txt')).splitlines()
    facts['excursions'] = len([ln for ln in findings
                               if ln.strip() and not ln.lstrip().startswith('#')])

    # R5 - a hand-rolled sign helper is a local re-implementation of _emission's
    # signed_term/rect_str. Count the ones still standing.
    sig = os.path.join(_ROOT, 'data', 'templates', 'branches')
    helpers = []
    for rel, needle in (
        (os.path.join('electrical_engineering', 'electromagnetics_and_waves',
                      'waves_and_phasors.py'), 'def _signed_term'),
        (os.path.join('electrical_engineering', 'signals_and_systems',
                      'discrete_time_signals.py'), 'C_str = f"+ '),
        (os.path.join('electrical_engineering', 'signals_and_systems',
                      'discrete_time_signals.py'), 'fmt = lambda'),
    ):
        if needle in _read(os.path.join(sig, rel)):
            helpers.append('%s :: %s' % (os.path.basename(rel), needle))
    facts['sign_helpers'] = helpers

    from tests.comparators.bindings import BINDINGS, UNBOUND
    facts['bound'] = len(BINDINGS)
    facts['unbound'] = len(UNBOUND)
    facts['by_kind'] = {}
    for b in BINDINGS.values():
        facts['by_kind'][b['kind']] = facts['by_kind'].get(b['kind'], 0) + 1
    facts['false_accept_unbound'] = {
        t.replace('template_', ''): int(re.search(r'(\d+) false accepts', r).group(1))
        for t, r in UNBOUND.items() if re.search(r'\d+ false accepts', r)}
    return facts


def _origin_silent(n_seeds):
    """R7 - gold answers of signal_operations that state no n=0 origin."""
    import random  # noqa: F401  (generate seeds it)
    from tests.template_integrity.core import discover, generate
    ref = [r for r in discover() if r.template_id == 'template_signal_operations'][0]
    silent = 0
    for s in range(n_seeds):
        tail = generate(ref, s, capture=False).solution.rsplit('**Answer:**', 1)[-1]
        if '*' not in tail.replace('**', ''):
            silent += 1
    return silent


# -------------------------------------------------------------------- checks

def _num(text):
    return int(text.replace(',', '').replace('*', '').strip())


def check(docs, facts):
    """docs maps a label to register TEXT. Returns a list of failures."""
    bad = []
    p6, c3 = docs['phase6'], docs['c3']

    # R1 - the phase6 section-7 table: | branch | then | now |
    for br in BRANCHES + ('total',):
        pat = r'\|\s*\*{0,2}%s\*{0,2}\s*\|\s*\*{0,2}\d+\*{0,2}\s*\|\s*\*{0,2}(\d+)\*{0,2}\s*\|' % br
        m = re.search(pat, p6)
        if not m:
            bad.append('R1 phase6: no row for %r in the [UNVERIFIED] table' % br)
        elif _num(m.group(1)) != facts['unverified'][br]:
            bad.append('R1 phase6 %s: states %s, tree has %d'
                       % (br, m.group(1), facts['unverified'][br]))

    # R1 - the C3 "Measured now" sentence, and its section heading
    m = re.search(r'Measured now\*{0,2}\s*\([^)]*\):\s*(.+?)%s\s*\*\*(\d+)\*\*' % _EMDASH,
                  c3, re.S)
    if not m:
        bad.append('R1 c3: the "Measured now" per-branch sentence is gone')
    else:
        for br in BRANCHES:
            b = re.search(r'%s\s+(\d+)' % br, m.group(1))
            if not b:
                bad.append('R1 c3: "Measured now" names no count for %s' % br)
            elif _num(b.group(1)) != facts['unverified'][br]:
                bad.append('R1 c3 %s: states %s, tree has %d'
                           % (br, b.group(1), facts['unverified'][br]))
        if _num(m.group(2)) != facts['unverified']['total']:
            bad.append('R1 c3 total: states %s, tree has %d'
                       % (m.group(2), facts['unverified']['total']))
    m = re.search(r'##\s*2\.\s*Unchecked\s*%s\s*\*\*(\d+)\*\*' % _EMDASH, c3)
    if not m:
        bad.append('R1 c3: section 2 heading no longer states a total')
    elif _num(m.group(1)) != facts['unverified']['total']:
        bad.append('R1 c3 heading: states %s, tree has %d'
                   % (m.group(1), facts['unverified']['total']))

    # R1 - the C3 re-derivation recipe's inline expectations
    for tag, key in (('KNOWN-DEFECTIVE', 'known_defective'), ('UNVERIFIED', 'unverified')):
        m = re.search(r'\[%s\\?\].*?#\s*(\d+)' % tag, c3)
        if not m:
            bad.append('R1 c3: the grep recipe for [%s] states no count' % tag)
        elif _num(m.group(1)) != facts[key]['total']:
            bad.append('R1 c3 grep recipe [%s]: states %s, tree has %d'
                       % (tag, m.group(1), facts[key]['total']))

    # R2 - the deletion sentence names rows; every one must be gone
    m = re.search(r'deleted seven substance-defect rows\*{0,2},?\s*(.+?from chemical)', p6, re.S)
    if not m:
        bad.append('R2: phase6 no longer names the deleted rows')
    else:
        named = re.findall(r'`([^`]+)`', m.group(1))
        if len(named) != 7:
            bad.append('R2: the deletion sentence names %d rows, not 7' % len(named))
        for name in named:
            if name in facts['keys']:
                bad.append('R2: %r is named as deleted and is still a live key' % name)
    if 'Acetylene Tetrabromide' not in facts['keys']:
        bad.append('R2: the kept name Acetylene Tetrabromide is gone from the tree')

    # R3 - the C3 reason table sums to what it declares
    m = re.search(r'The remaining (\d+) are not one decision.*?sums to (\d+)', c3, re.S)
    if not m:
        bad.append('R3: the C3 reason table no longer declares its total')
    else:
        declared, sums_to = _num(m.group(1)), _num(m.group(2))
        if declared != sums_to or declared != facts['unverified']['mechanical']:
            bad.append('R3: reason table declares %d / "sums to %d", tree has %d'
                       % (declared, sums_to, facts['unverified']['mechanical']))
        block = c3.split('sums to %d' % sums_to)[1].split('**The shape of that table')[0]
        rows = re.findall(r'^\s*\|\s*(\d+|%s)\s*\|\s*(\d+)\s*\|' % _EMDASH, block, re.M)
        now = sum(int(a) for a, _ in rows if a != _EMDASH)
        if now != sums_to:
            bad.append('R3: reason rows sum to %d, the table says %d' % (now, sums_to))

    # R4 - excursions
    m = re.search(r'\*\*(\d+)\*\*\s*consumer-domain excursion', c3)
    if not m:
        bad.append('R4: the C3 register no longer states an excursion count')
    elif _num(m.group(1)) != facts['excursions']:
        bad.append('R4: register states %s excursion(s), domain_findings.txt lists %d'
                   % (m.group(1), facts['excursions']))

    # R5 - sign helpers
    m = re.search(r'##\s*3\.\s*(\w+) hand-rolled sign helper', p6)
    want = _ORDINAL.get(len(facts['sign_helpers']), str(len(facts['sign_helpers'])))
    if not m:
        bad.append('R5: the sign-helper section heading is gone')
    elif m.group(1) != want:
        bad.append('R5: heading says %r, the tree has %d (%s)'
                   % (m.group(1), len(facts['sign_helpers']), '; '.join(facts['sign_helpers'])))

    # R6 - bindings
    m = re.search(r'\*\*(\d+) bound / (\d+) unbound\*\*', p6)
    if not m:
        bad.append('R6: phase6 no longer states the binding counts')
    else:
        if _num(m.group(1)) != facts['bound'] or _num(m.group(2)) != facts['unbound']:
            bad.append('R6: states %s bound / %s unbound, live table has %d / %d'
                       % (m.group(1), m.group(2), facts['bound'], facts['unbound']))
    m = re.search(r'`symbolic` at (\d+) of 9', p6)
    if not m:
        bad.append('R6: phase6 no longer states the symbolic binding count')
    elif _num(m.group(1)) != facts['by_kind'].get('symbolic', 0):
        bad.append('R6 symbolic: states %s of 9, live table has %d'
                   % (m.group(1), facts['by_kind'].get('symbolic', 0)))
    for tmpl, n in facts['false_accept_unbound'].items():
        m = re.search(r'`%s`[^|]*\|\s*(\d+)(?:\s*each)?\s*\|' % re.escape(tmpl), p6)
        if not m:
            bad.append('R6: %s is unbound on N false accepts and is not in the table' % tmpl)
        elif _num(m.group(1)) != n:
            bad.append('R6 %s: table says %s false accepts, live reason says %d'
                       % (tmpl, m.group(1), n))
    listed = re.search(r'\| unbound with this reason \| count \|(.+?)\n\n', p6, re.S)
    if listed:
        names = set(re.findall(r'`([a-z0-9_]+)`', listed.group(1)))
        for extra in sorted(names - set(facts['false_accept_unbound'])):
            bad.append('R6: %s is listed in the table and is no longer unbound that way'
                       % extra)

    # R7 - D-050, re-measured at the N the document states
    m = re.search(r'(\d+) of ([\d,]+) seeds silent, ([\d.]+)%', p6)
    if not m:
        bad.append('R7: phase6 no longer states the D-050 measurement')
    else:
        stated, n_seeds = _num(m.group(1)), _num(m.group(2))
        silent = _origin_silent(n_seeds)
        if silent != stated:
            bad.append('R7: states %d of %d seeds silent, re-measured %d'
                       % (stated, n_seeds, silent))
        pct = round(100.0 * silent / n_seeds, 1)
        if abs(pct - float(m.group(3))) > 0.05:
            bad.append('R7: states %s%%, re-measured %.1f%%' % (m.group(3), pct))
    return bad


# ----------------------------------------------------------------- entrypoint

def run():
    facts = measure()
    docs = {'phase6': _read(PHASE6), 'c3': _read(C3)}
    bad = check(docs, facts)
    print('registers vs tree: [UNVERIFIED] %d (%s), [KNOWN-DEFECTIVE] %d, '
          'excursions %d, sign helpers %d, bindings %d bound / %d unbound'
          % (facts['unverified']['total'],
             ', '.join('%s %d' % (b, facts['unverified'][b]) for b in BRANCHES),
             facts['known_defective']['total'], facts['excursions'],
             len(facts['sign_helpers']), facts['bound'], facts['unbound']))
    for b in bad:
        print('  - ' + b)
    print('%d disagreement(s) between the registers and the tree' % len(bad))
    return 1 if bad else 0


# A plant is judged by the failure it ADDS over a clean fixture (the rule
# test_table_metadata records). Each mutates a register's TEXT; the measurements
# come from the tree either way, so a plant that is not detected means the
# document could say something the tree contradicts and nothing would notice.
_PLANTS = {
    'R1 phase6 total': ('phase6', '**170**', '**171**'),
    'R1 phase6 branch': ('phase6', '| mechanical | 148 | **142** |',
                         '| mechanical | 148 | **140** |'),
    'R1 c3 heading': ('c3', 'Unchecked — **170**', 'Unchecked — **177**'),
    'R1 c3 sentence': ('c3', 'mechanical 142, chemical 12', 'mechanical 148, chemical 13'),
    'R1 c3 recipe': ('c3', '# 170 (198 at the review)', '# 177 (198 at the review)'),
    'R2 deleted row still live': ('phase6', '`R-410A` from chemical',
                                  '`Acetylene Tetrabromide` from chemical'),
    'R3 reason row': ('c3', '| 39 | 39 | no source on disk *for this material*',
                      '| 38 | 39 | no source on disk *for this material*'),
    'R4 excursions': ('c3', '| **1** consumer-domain excursion',
                      '| **2** consumer-domain excursions'),
    'R5 helper count': ('phase6', '## 3. One hand-rolled sign helper',
                        '## 3. Three hand-rolled sign helpers'),
    'R6 bound counts': ('phase6', '**127 bound / 23 unbound**', '**120 bound / 30 unbound**'),
    'R6 symbolic': ('phase6', '`symbolic` at 6 of 9', '`symbolic` at 1 of 9'),
    'R6 per-template': ('phase6', '| `null_to_null_bandwidth` | 24 |',
                        '| `null_to_null_bandwidth` | 25 |'),
    'R7 D-050': ('phase6', '329 of 2,000 seeds silent, 16.4%',
                 '300 of 2,000 seeds silent, 15.0%'),
    # Second plant per property, in a materially different surface form
    # (SPEC-CHANGE 17): a wrong number is one way to drift, and the statement
    # going missing entirely is the other. Phase 5 shipped four detectors that
    # were narrower than their class because every plant was a shape the regex
    # already matched.
    'R2 sentence gone': ('phase6', 'deleted seven substance-defect rows',
                         'deleted some rows'),
    'R3 declared total': ('c3', 'sums to 142', 'sums to 143'),
    'R4 count unbolded': ('c3', '| **1** consumer-domain excursion',
                          '| One consumer-domain excursion'),
    'R5 heading gone': ('phase6', '## 3. One hand-rolled sign helper remains',
                        '## 3. Sign helpers'),
    'R7 seed count': ('phase6', '329 of 2,000 seeds silent',
                      '329 of 1,000 seeds silent'),
}


def selftest():
    facts = measure()
    clean_docs = {'phase6': _read(PHASE6), 'c3': _read(C3)}
    clean = check(clean_docs, facts)
    if clean:
        print('  the registers themselves disagree with the tree - no plant can be judged:')
        for x in clean:
            print('    ' + x)
        return 1
    bad = []
    for label, (doc, old, new) in _PLANTS.items():
        if clean_docs[doc].count(old) != 1:
            bad.append('%s: anchor %r occurs %d times, not once'
                       % (label, old, clean_docs[doc].count(old)))
            continue
        planted = dict(clean_docs)
        planted[doc] = planted[doc].replace(old, new)
        found = [x for x in check(planted, facts) if x not in clean]
        if not any(x.startswith(label.split()[0]) for x in found):
            bad.append('%s: NOT detected (got %s)' % (label, found))
    for label in _PLANTS:
        print('  [%s] %s' % ('FAIL' if any(b.startswith(label) for b in bad) else 'ok', label))
    for b in bad:
        print('  - ' + b)
    print('selftest: %d failure(s)' % len(bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
