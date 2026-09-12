"""D6.12 - the testset generator.

Replaces the 20 hand-rolled `main()` blocks scattered through the template
modules, which between them had four defects:

* they covered **three branches of five** - civil and industrial had none, which
  is why the published pool was 90 templates and not 150;
* they drew each item's seed from an **unseeded** global RNG, so the recorded
  seed reproduced an item but nothing reproduced the *set*;
* they seeded only `random`, never `numpy` - so `levenspiel_plot_interpretation`
  (which calls `np.random.uniform`) could not be regenerated from its own
  recorded seed at all.  That is the audit's T3 finding;
* their instance counts ranged from 3 to 200, giving `mole_balances` 1,000
  records and `magnetostatics` 3.

This module fixes all four by building on `tests.template_integrity.core`, which
already enumerates all 150 templates and owns the seeding path the integrity
suite uses.

Determinism, concretely: a given ``--master-seed`` produces a byte-identical
testset.  Per-item seeds come from BLAKE2b over ``(master_seed, template_id,
index)`` - **not** Python's ``hash()``, which is salted per process by
``PYTHONHASHSEED`` and would silently make every run different.

Usage::

    python generate_testset.py                    # 150 x 15 -> testset/
    python generate_testset.py --verify           # prove two runs agree
    python generate_testset.py --branch civil_engineering --instances 3

**Ordering constraint (Phase 6).**  Do not generate the real pool until D6.7
(doubled sign, 14 templates) and D6.11 (per-part units) have merged.  Both change
emitted text; a pool built before them is stale on arrival and re-runs inference
across 11 models for nothing.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter

REPO_ROOT = os.path.abspath(os.path.dirname(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from tests.template_integrity.core import TemplateRef, discover, generate  # noqa: E402

INVENTORY = os.path.join('docs', 're-implementation-sep', 'template_inventory.csv')

#: Changing this changes every item.  It is recorded in the manifest so a pool
#: can always be traced back to the seed that produced it.
DEFAULT_MASTER_SEED = 20260912

DEFAULT_INSTANCES = 15          # the published pool's per-template depth
SEED_SPACE = 2 ** 31 - 1


def item_seed(master_seed: int, template_id: str, index: int) -> int:
    """A stable per-item seed.

    Deliberately NOT ``hash()``: that is salted per process, so a pool would
    differ between runs on the same machine while looking deterministic.
    """
    key = f'{master_seed}:{template_id}:{index}'.encode('utf-8')
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), 'big') % SEED_SPACE


def difficulty_map(path: str = INVENTORY) -> dict[str, str]:
    """template_id -> Easy | Intermediate | Advanced, from the audit inventory.

    The inventory's ids match ``discover()`` exactly (verified: zero difference
    in either direction), so a missing id here is a real inconsistency and is
    raised rather than defaulted.
    """
    full = os.path.join(REPO_ROOT, path)
    with open(full, encoding='utf-8-sig', newline='') as fh:
        return {r['template_id']: r['difficulty'] for r in csv.DictReader(fh)}


def record(ref: TemplateRef, seed: int, question: str, solution: str,
           level: str) -> dict:
    """One testset record, in the published schema.

    The field names are the published ones and their meanings are preserved:
    ``domain`` is the DIRECTORY (``electromagnetics_and_waves``) and ``area`` is
    the MODULE (``waves_and_phasors``).  Note that `core.TemplateRef` calls the
    directory ``area``; swapping the two here would mislabel every record.
    """
    return {
        'seed': seed,
        'branch': ref.branch,
        'domain': ref.area,
        'area': os.path.basename(ref.file_path)[:-3],
        'id': ref.template_id[len('template_'):],
        'level': level,
        'question': question,
        'solution': solution,
    }


def build(instances: int = DEFAULT_INSTANCES,
          master_seed: int = DEFAULT_MASTER_SEED,
          branch: str | None = None) -> tuple[list[dict], list[str]]:
    """Generate the pool.  Returns (records, errors)."""
    levels = difficulty_map()
    refs = discover(branch)
    records: list[dict] = []
    errors: list[str] = []

    for ref in refs:
        level = levels.get(ref.template_id)
        if level is None:
            errors.append(f'{ref.template_id}: not in {INVENTORY}')
            continue
        for i in range(instances):
            seed = item_seed(master_seed, ref.template_id, i)
            inst = generate(ref, seed, capture=False)
            if not inst.ok:
                errors.append(f'{ref.template_id} seed={seed}: {inst.error}')
                continue
            records.append(record(ref, seed, inst.question, inst.solution, level))
    return records, errors


def write(records: list[dict], out_dir: str) -> dict[str, int]:
    """Write one .jsonl per module, mirroring the published layout."""
    by_file: dict[str, list[dict]] = {}
    for r in records:
        rel = os.path.join(out_dir, r['branch'], r['domain'], r['area'] + '.jsonl')
        by_file.setdefault(rel, []).append(r)

    for rel, rows in sorted(by_file.items()):
        os.makedirs(os.path.dirname(rel), exist_ok=True)
        with open(rel, 'w', encoding='utf-8') as fh:
            for row in rows:
                fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    return {k: len(v) for k, v in by_file.items()}


def report(records: list[dict], errors: list[str]) -> int:
    """Print the counts, and every predicate behind them."""
    by_branch = Counter(r['branch'] for r in records)
    by_level = Counter(r['level'] for r in records)
    per_template = Counter(r['id'] for r in records)

    print(f'records          {len(records)}')
    print(f'templates        {len(per_template)}')
    print(f'by branch        {dict(sorted(by_branch.items()))}')
    print(f'by level         {dict(sorted(by_level.items()))}')

    counts = set(per_template.values())
    print(f'per template     {sorted(counts)}'
          f'{"  UNIFORM" if len(counts) == 1 else "  *** NOT UNIFORM ***"}')

    # A duplicate question is not fatal, but it over-weights whatever it repeats,
    # and one module's own comments record this happening before.
    dupes = len(records) - len({r['question'] for r in records})
    print(f'duplicate questions  {dupes}')

    if errors:
        print(f'\nGENERATION ERRORS  {len(errors)}')
        for e in errors[:20]:
            print(f'  {e}')
        if len(errors) > 20:
            print(f'  ... and {len(errors) - 20} more')
    return len(errors)


def verify(instances: int, master_seed: int, branch: str | None) -> int:
    """Generate twice in this process and confirm the pools are identical.

    This is the weaker of the two checks the exit gate wants.  It catches an
    unseeded draw inside a template, which is what it is for.  It does NOT catch
    state that survives a process - run this module twice and diff the output
    files for that.
    """
    first, e1 = build(instances, master_seed, branch)
    second, e2 = build(instances, master_seed, branch)
    if e1 or e2:
        print(f'errors: {len(e1)} / {len(e2)} - cannot verify over a failing pool')
        return 1
    if first == second:
        print(f'VERIFY OK - two runs, {len(first)} records, byte-identical')
        return 0

    diff = [a['id'] for a, b in zip(first, second) if a != b]
    print(f'VERIFY FAILED - {len(diff)} records differ across two runs '
          f'at the same master seed')
    for t in sorted(set(diff))[:10]:
        print(f'  {t}')
    return 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', default='testset')
    ap.add_argument('--instances', type=int, default=DEFAULT_INSTANCES)
    ap.add_argument('--master-seed', type=int, default=DEFAULT_MASTER_SEED)
    ap.add_argument('--branch', default=None)
    ap.add_argument('--verify', action='store_true',
                    help='generate twice and compare; write nothing')
    ap.add_argument('--dry-run', action='store_true',
                    help='generate and report; write nothing')
    args = ap.parse_args()

    if args.verify:
        return verify(args.instances, args.master_seed, args.branch)

    records, errors = build(args.instances, args.master_seed, args.branch)
    n_err = report(records, errors)

    if errors:
        print('\nRefusing to write a pool with generation errors.')
        return 1
    if args.dry_run:
        print('\n--dry-run: nothing written')
        return 0

    files = write(records, args.out)
    print(f'\nwrote {len(files)} files under {args.out}/  '
          f'(master_seed={args.master_seed}, instances={args.instances})')
    return 0 if n_err == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
