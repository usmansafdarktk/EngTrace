"""E1 (WS-E): adopt the experts' majority levels. Run only when decision D7 says so.

    python -m template_annotation_23092026.levels.apply_labels            # FREE: prints the change, writes nothing
    python -m template_annotation_23092026.levels.apply_labels --write    # changes the level source in place
    python -m template_annotation_23092026.levels.apply_labels --selftest # FREE: on a copy of the source, temp folder

THE LEVEL SOURCE. generate_testset.difficulty_map() reads the `difficulty` column of
docs/re-implementation-sep/audit/template_inventory.csv; freeze.py copies it into each manifest row's `level` and
into FREEZE.json's `by_level`, and score.py copies the manifest's level into every store row. Nothing else assigns
a level (tests/template_integrity/regen_inventory.py carries the column over unchanged).

THE RULE. A template whose raters' majority (agreement.json, templates.<id>.majority: the level more than half of
its raters chose) differs from its current label takes the majority. A template without a majority keeps its
label. The run refuses when agreement.json's current label for a template is not the source's (the ratings would
have been scored against another labelling). Only the difficulty cell of a changed row is rewritten; every other
byte of the file, line endings included, stays as it is.

AFTER --write (the output says so too): the instances do not change, only the label, but the manifest's `level`
field, FREEZE.json's `by_level` and its `manifest_sha256` must be refreshed from the new source, and freeze.py has no
mode that rewrites the labels alone; WS-G owns that step (a relabel record in FREEZE.json, then `freeze --verify`,
which compares a regenerated manifest, levels included, with the one on disk). Then the one re-score refreshes the
stores' level field, and every level table is regenerated.
"""
from __future__ import annotations

import argparse
import collections
import csv
import io
import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import generate_testset  # noqa: E402

SOURCE = REPO / generate_testset.INVENTORY
AGREEMENT = HERE / 'agreement.json'
LEVELS = ['Easy', 'Intermediate', 'Advanced']

AFTER = """NEXT (WS-G): the label source changed; the instances did not. Refresh, in this order:
  1. the manifest's `level` field and FREEZE.json's `by_level` and `manifest_sha256` from the new source
     (full_run_28092026/freeze.py has no labels-only mode: WS-G adds one, recorded in FREEZE.json as a relabel,
     or re-runs the build), then `python -m full_run_28092026.freeze --verify`;
  2. the one re-score of every store (rows carry the manifest's level);
  3. every level table and figure (docs/appendix_statistics.py --write, paper_results.py --write)."""


def changes(agreement: dict, current: dict[str, str]) -> list[dict]:
    """The rows to change: majority present and different from the source's label."""
    out = []
    for tid, t in sorted(agreement['templates'].items()):
        if tid not in current:
            raise SystemExit(f'{tid}: in agreement.json but not in the level source')
        if t['current'] != current[tid]:
            raise SystemExit(f"{tid}: agreement.json was scored against {t['current']}, the source says {current[tid]}")
        if t['majority'] and t['majority'] != current[tid]:
            out.append({'template_id': tid, 'branch': t['branch'], 'old': current[tid], 'new': t['majority'],
                        'ratings': t['ratings']})
    return out


def rewrite(text: str, new: dict[str, str]) -> str:
    """The source with the difficulty cell of each named row replaced and nothing else touched."""
    lines = text.splitlines(keepends=True)
    header = next(csv.reader([lines[0]]))
    ti, di = header.index('template_id'), header.index('difficulty')
    if di < ti:
        raise SystemExit('unexpected column order in the level source')
    done = set()
    for n, line in enumerate(lines[1:], 1):
        fields = next(csv.reader([line]))
        if not fields or fields[ti] not in new:
            continue
        parts = line.split(',', di + 1)
        if parts[:di + 1] != fields[:di + 1]:
            raise SystemExit(f'{fields[ti]}: a quoted field before the difficulty column; refusing a raw edit')
        parts[di] = new[fields[ti]]
        lines[n] = ','.join(parts)
        done.add(fields[ti])
    if done != set(new):
        raise SystemExit(f'not found in the source: {sorted(set(new) - done)}')
    return ''.join(lines)


def apply(agreement_path: Path, source: Path, write: bool) -> list[dict]:
    agreement = json.loads(agreement_path.read_text(encoding='utf-8'))
    current = generate_testset.difficulty_map(str(source))
    todo = changes(agreement, current)
    before = collections.Counter(current.values())
    after = collections.Counter({**current, **{c['template_id']: c['new'] for c in todo}}.values())
    no_majority = sorted(t for t, v in agreement['templates'].items() if not v['majority'])
    print(f'{len(todo)} of {len(current)} templates change level; {len(no_majority)} have no majority and keep theirs.')
    if todo:
        print('\n| template | branch | old | new | ratings (E / I / A) |\n|---|---|---|---|---|')
        for c in todo:
            r = c['ratings']
            print(f"| `{c['template_id']}` | {c['branch'].replace('_engineering', '')} | {c['old']} | {c['new']} | "
                  f"{r['Easy']} / {r['Intermediate']} / {r['Advanced']} |")
    print('\nlevels before: ' + ', '.join(f'{lv} {before[lv]}' for lv in LEVELS)
          + '; after: ' + ', '.join(f'{lv} {after[lv]}' for lv in LEVELS))
    if not write:
        print('\n(dry run: nothing written; --write applies it)')
        return todo
    if todo:
        with open(source, encoding='utf-8', newline='') as fh:
            text = fh.read()
        out = rewrite(text, {c['template_id']: c['new'] for c in todo})
        with open(source, 'w', encoding='utf-8', newline='') as fh:
            fh.write(out)
        check = generate_testset.difficulty_map(str(source))
        assert all(check[c['template_id']] == c['new'] for c in todo)
        assert all(check[t] == current[t] for t in current if t not in {c['template_id'] for c in todo})
        print(f'\nwritten: {source}')
        print(AFTER)
    return todo


def selftest() -> int:
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-apply-labels-'))
    try:
        src = tmp / 'template_inventory.csv'
        shutil.copyfile(SOURCE, src)
        original = src.read_bytes()
        current = generate_testset.difficulty_map(str(src))
        tids = sorted(current)
        up = next(t for t in tids if current[t] == 'Easy')
        down2 = next(t for t in tids if current[t] == 'Advanced')
        split = next(t for t in tids if t not in (up, down2))
        same = next(t for t in tids if t not in (up, down2, split))
        templates = {}
        for t in tids:
            maj = {up: 'Intermediate', down2: 'Easy', split: None}.get(t, current[t])
            templates[t] = {'branch': 'x', 'current': current[t], 'majority': maj,
                            'ratings': {'Easy': 1, 'Intermediate': 1, 'Advanced': 1}}
        ag = tmp / 'agreement.json'
        ag.write_text(json.dumps({'templates': templates}), encoding='utf-8')
        todo = apply(ag, src, write=False)
        assert src.read_bytes() == original, 'a dry run wrote'
        assert {c['template_id'] for c in todo} == {up, down2}, todo
        apply(ag, src, write=True)
        now = generate_testset.difficulty_map(str(src))
        assert now[up] == 'Intermediate' and now[down2] == 'Easy' and now[split] == current[split] \
            and now[same] == current[same], 'the labels did not change as the rule says'
        old_lines = original.decode('utf-8').splitlines(keepends=True)
        new_lines = src.read_bytes().decode('utf-8').splitlines(keepends=True)
        assert len(old_lines) == len(new_lines)
        moved = [(a, b) for a, b in zip(old_lines, new_lines) if a != b]
        assert len(moved) == 2, len(moved)
        for a, b in moved:                                          # only the difficulty cell differs
            fa, fb = next(csv.reader(io.StringIO(a))), next(csv.reader(io.StringIO(b)))
            assert [x for i, x in enumerate(fa) if i != 5] == [x for i, x in enumerate(fb) if i != 5]
            assert a.endswith('\r\n') == b.endswith('\r\n')
        templates[same]['current'] = 'Advanced' if current[same] != 'Advanced' else 'Easy'
        ag.write_text(json.dumps({'templates': templates}), encoding='utf-8')
        try:
            apply(ag, src, write=False)
            raise AssertionError('a mismatched current label was accepted')
        except SystemExit:
            pass
        assert SOURCE.read_bytes() == original, 'the real source changed during the self-test'
        print('\nSELFTEST OK: the dry run writes nothing; --write changes exactly the two rows the rule names, only '
              'their difficulty cell, line endings kept; no majority keeps the label; a mismatched current label stops '
              'the run; the real source is untouched')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--write', action='store_true', help='apply the change to the level source (only on D7)')
    ap.add_argument('--agreement', type=Path, default=AGREEMENT)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.agreement.exists():
        raise SystemExit(f'no {a.agreement.name}: score the returns first (score_levels.py)')
    apply(a.agreement, SOURCE, a.write)
    return 0


if __name__ == '__main__':
    sys.exit(main())
