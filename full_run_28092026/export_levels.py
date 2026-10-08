"""Each template's current difficulty level as JSONL, one line per template: template_id, branch,
domain, area and level (Easy, Intermediate or Advanced).

The level is the manifest's, which freeze.py copied from the audit inventory's difficulty column
(generate_testset.difficulty_map). The script stops if the two disagree or if a template's 15 items
carry different levels. Re-run it after any change to the labels (D7).

  python -m full_run_28092026.export_levels           template_levels/all.jsonl and one file per branch
  python -m full_run_28092026.export_levels --check   the files on disk equal what the sources give
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'manifest.jsonl'
INVENTORY = HERE.parent / 'docs' / 're-implementation-sep' / 'audit' / 'template_inventory.csv'
OUT = HERE / 'template_levels'
FIELDS = ('template_id', 'branch', 'domain', 'area', 'level')
LEVELS = ('Easy', 'Intermediate', 'Advanced')


def templates() -> list[dict]:
    """One record per template, sorted by branch, domain, area and id."""
    rows: dict[str, dict] = {}
    for line in MANIFEST.read_text(encoding='utf-8').splitlines():
        r = json.loads(line)
        rec = {k: r[k] for k in FIELDS}
        if rows.setdefault(r['template_id'], rec) != rec:
            raise SystemExit(f"{r['template_id']}: its items disagree on {FIELDS}")
    with INVENTORY.open(encoding='utf-8-sig', newline='') as fh:
        inventory = {r['template_id']: r['difficulty'] for r in csv.DictReader(fh)}
    if set(inventory) != set(rows):
        raise SystemExit(f'the manifest and the inventory list different templates: {sorted(set(inventory) ^ set(rows))}')
    differ = sorted(t for t, rec in rows.items() if rec['level'] != inventory[t])
    if differ:
        raise SystemExit(f'the manifest and the inventory give different levels for {len(differ)} templates: {differ}')
    unknown = sorted(t for t, rec in rows.items() if rec['level'] not in LEVELS)
    if unknown:
        raise SystemExit(f'levels outside {LEVELS}: {unknown}')
    return sorted(rows.values(), key=lambda r: (r['branch'], r['domain'], r['area'], r['template_id']))


def files(recs: list[dict]) -> dict[str, str]:
    """File name -> JSONL text: all.jsonl, then <branch>.jsonl per branch."""
    def text(rs: list[dict]) -> str:
        return ''.join(json.dumps(r) + '\n' for r in rs)
    out = {'all.jsonl': text(recs)}
    for b in sorted({r['branch'] for r in recs}):
        out[f'{b}.jsonl'] = text([r for r in recs if r['branch'] == b])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true', help='compare the files on disk with the sources; write nothing')
    args = ap.parse_args()
    recs = templates()
    want = files(recs)
    if args.check:
        stale = [n for n, t in want.items() if not (OUT / n).exists() or (OUT / n).read_text(encoding='utf-8') != t]
        if stale:
            raise SystemExit(f'stale or missing: {stale}; re-run without --check')
        print(f'{OUT.name}/: {len(want)} files current')
    else:
        OUT.mkdir(exist_ok=True)
        for n, t in want.items():
            with (OUT / n).open('w', encoding='utf-8', newline='') as fh:
                fh.write(t)
        print(f'{OUT.name}/: wrote {len(want)} files')
    counts = Counter((r['branch'], r['level']) for r in recs)
    print(f"{'branch':<24}" + ''.join(f'{lv:>14}' for lv in LEVELS) + f"{'templates':>11}")
    for b in sorted({r['branch'] for r in recs}) + ['all']:
        n = [sum(v for (br, lv), v in counts.items() if lv == level and b in (br, 'all')) for level in LEVELS]
        print(f'{b:<24}' + ''.join(f'{x:>14}' for x in n) + f'{sum(n):>11}')


if __name__ == '__main__':
    main()
