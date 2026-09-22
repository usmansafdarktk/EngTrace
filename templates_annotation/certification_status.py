"""Which templates does the December 2025 expert certification still describe?

Two questions, both answered from the repository rather than from memory:

1. Of the 90 templates certified in December 2025 (Appendix K of the May 2026
   paper), which have changed since? A template is counted as changed if its
   function source differs between the certification-time commit and HEAD, or
   if its emitted instances moved because a constants table it draws from was
   corrected (Track B, recorded in phaseC3_corrections_item_pool_impact.md).
2. What does the December review record itself show? Decisions, score
   distribution, and the time between consecutive submissions per reviewer.
   The review JSONL is git-ignored, so this part runs only where those files
   exist; the census in part 1 needs nothing but git history.

    python -m templates_annotation.certification_status
"""
from __future__ import annotations

import ast
import collections
import difflib
import glob
import json
import re
import statistics
import subprocess
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
# Last commit before the December 2025 expert reviews (reviews are timestamped
# 2025-12-15/16; this commit is 2025-12-17 and touches only the Tribunal
# consensus script, not a template).
CERTIFICATION_REF = '7df55a2'
IMPACT_DOC = REPO / 'docs/re-implementation-sep/track-b/phaseC3_corrections_item_pool_impact.md'
REVIEWS = REPO / 'templates_annotation/annotation_app/reviews'
LEGACY = ('chemical_engineering', 'electrical_engineering', 'mechanical_engineering')


def functions_at(rev: str) -> dict[str, tuple[str, str]]:
    """{template name: (branch, function source)} for every template_* at rev."""
    paths = subprocess.run(['git', 'ls-tree', '-r', '--name-only', rev, 'data/templates/branches'],
                           capture_output=True, text=True, cwd=REPO, check=True).stdout.split()
    out = {}
    for p in paths:
        if not p.endswith('.py') or p.endswith(('constants.py', '__init__.py', '_emission.py')):
            continue
        src = subprocess.run(['git', 'show', f'{rev}:{p}'], capture_output=True, text=True,
                             encoding='utf8', cwd=REPO, check=True).stdout
        lines = src.splitlines()
        for node in ast.parse(src).body:
            if isinstance(node, ast.FunctionDef) and node.name.startswith('template_'):
                out[node.name] = (p.split('/')[3], '\n'.join(lines[node.lineno - 1:node.end_lineno]))
    return out


def constants_moved() -> set[str]:
    """Templates whose item pool moved under the Track B constants corrections."""
    names = set()
    for line in IMPACT_DOC.read_text(encoding='utf8').splitlines():
        m = re.match(r'\|\s*`([a-z0-9_]+)`\s*\|', line)
        if m:
            names.add('template_' + m.group(1))
    return names


def census() -> None:
    norm = lambda s: re.sub(r'\s+', ' ', s).strip()  # noqa: E731
    then, now = functions_at(CERTIFICATION_REF), functions_at('HEAD')
    changed = {}
    for name, (_, src) in then.items():
        if name in now and norm(src) != norm(now[name][1]):
            ratio = difflib.SequenceMatcher(None, src.splitlines(), now[name][1].splitlines()).ratio()
            changed[name] = 1 - ratio
    moved = constants_moved()
    moved_only = sorted(n for n in then if n in moved and n not in changed)
    new = sorted(n for n in now if n not in then)
    stale = set(changed) | set(moved_only)

    print(f'certified at {CERTIFICATION_REF}: {len(then)} templates; at HEAD: {len(now)}')
    print(f'  code changed since certification : {len(changed)} '
          f'{dict(collections.Counter(now[n][0].split("_")[0] for n in changed))}')
    print(f'    of which >=20% of lines differ : {sum(1 for v in changed.values() if v >= 0.2)}')
    print(f'  code unchanged, outputs moved by constants corrections: {len(moved_only)}')
    print(f'  legacy needing re-certification  : {len(stale)}')
    print(f'  legacy the Dec 2025 certification still describes: {len(then) - len(stale)}')
    print(f'  never certified (added Sep 2026) : {len(new)} '
          f'{dict(collections.Counter(now[n][0].split("_")[0] for n in new))}')
    print(f'  total needing certification      : {len(stale) + len(new)} of {len(now)}')
    print('\n  changed, by fraction of lines differing:')
    for n, v in sorted(changed.items(), key=lambda x: -x[1]):
        print(f'    {v:.2f} {now[n][0][:4]} {n}')
    print('  outputs moved by constants only:', ', '.join(m.replace('template_', '') for m in moved_only))


def review_record() -> None:
    files = sorted(glob.glob(str(REVIEWS / '*.jsonl')))
    if not files:
        print(f'\n(no review JSONL under {REVIEWS}; part 2 skipped)')
        return
    rows = [json.loads(l) for f in files for l in open(f, encoding='utf8') if l.strip()]
    print(f'\nDecember 2025 review record: {len(rows)} rows from {len(files)} files')
    print('  decisions:', dict(collections.Counter(r['decision'] for r in rows)))
    for dim in ('physical_plausibility', 'mathematical_correctness', 'pedagogical_clarity'):
        print(f'  {dim:25s}', dict(sorted(collections.Counter(r['scores'][dim] for r in rows).items())))
    all5 = sum(1 for r in rows if all(v == 5 for v in r['scores'].values()))
    print(f'  rows scoring 5 on every dimension: {all5} ({all5 / len(rows):.0%})')
    fb = [len(r.get('feedback') or '') for r in rows]
    print(f'  feedback length, chars: median {statistics.median(fb):.0f}, max {max(fb)}')
    print('  seconds between consecutive submissions (per reviewer, per branch):')
    by = collections.defaultdict(list)
    for r in rows:
        by[(r['annotator_id'], r['branch'])].append(datetime.fromisoformat(r['timestamp']))
    for (who, branch), ts in sorted(by.items()):
        ts.sort()
        gaps = [(b - a).total_seconds() for a, b in zip(ts, ts[1:])]
        span = (ts[-1] - ts[0]).total_seconds() / 60
        print(f'    {who:10s} {branch:24s} n={len(ts):2d} median {statistics.median(gaps):4.0f}s '
              f'min {min(gaps):3.0f}s  whole branch in {span:5.1f} min  ({ts[0].date()})')


if __name__ == '__main__':
    census()
    review_record()
