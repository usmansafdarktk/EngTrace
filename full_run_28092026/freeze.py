"""Freeze the full run's item pool: every template x 15 instances, from a private seed.

    python -m full_run_28092026.freeze                # build: generate, screen, write pool/ and the record
    python -m full_run_28092026.freeze --verify       # regenerate from SEED.secret; compare with both, byte for byte
    python -m full_run_28092026.freeze --check-files  # no seed needed: pool/ on disk against manifest.jsonl

What is committed and what is not (D-114). The templates and the generator's default
master seed are public, so a published seed publishes the items. The pool is therefore
drawn from a 128-bit seed kept in SEED.secret beside this file, gitignored together with
pool/, which holds the question and gold text. Committed: manifest.jsonl, one row per
item with its SHA-256 over question + NUL + solution (the pilot's convention,
evaluator_pilot_17092026/freeze.py) and neither text nor seed; and FREEZE.json, carrying
the SHA-256 of the seed as a commitment. At publication the seed is revealed, and anyone
can regenerate the pool and check it against both.

Selection rule. For each template, take instance indices 0, 1, 2, ... and keep the first
15 that are acceptable. An instance is replaced by the next index, and the replacement is
logged in FREEZE.json, when
  - its question repeats one already kept for the template (a duplicate double-weights it);
  - its question is one of the pilot slice's questions, which are public with their gold
    and on whose traces the answer check was tuned;
  - its gold solution carries an exact display tie on a line T1 can parse, the test of
    template_annotation_23092026/layer0/tie_census.py: such an instance has no gold value
    a decimal and a binary reader agree on (D-016). A tie on a line T1 cannot parse is
    not seen;
  - the template fails to generate it.
Every template has 15 instances, the owner's rule (D-114). A template whose question space
is exhausted, 75 indices yielding fewer than 15 distinct acceptable questions, keeps each
distinct question once and fills its 15 with repeats in index order. A repeat carries
repeat_of, the item_id of the first instance with its question, in the pool and in the
manifest, and FREEZE.json lists the template under short_templates, so the number of
distinct questions is never overstated.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import secrets
import subprocess
import sys
from decimal import Decimal, getcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generate_testset import INVENTORY, difficulty_map, item_seed, record, write  # noqa: E402
from tests.template_integrity.checks.t1_closure import check_line  # noqa: E402
from tests.template_integrity.core import discover, generate, printed_precision  # noqa: E402
from template_annotation_23092026.layer0.tie_census import dec_eval  # noqa: E402

getcontext().prec = 60

INSTANCES = 15
MAX_EXTRA = 60          # indices tried beyond 15 before a template is declared short
SEED_BITS = 128
SEED_PATH = HERE / 'SEED.secret'
POOL_DIR = HERE / 'pool'
MANIFEST = HERE / 'manifest.jsonl'
FREEZE = HERE / 'FREEZE.json'
PILOT_SLICE = REPO / 'evaluator_pilot_17092026' / 'slice' / 'manifest.jsonl'


def sha(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def item_sha(question: str, solution: str) -> str:
    return sha(question + '\x00' + solution)


def has_tie(solution: str) -> bool:
    """tie_census.census()'s test, applied to one instance."""
    for raw in solution.splitlines():
        line = raw.strip()
        if '=' not in line:
            continue
        checks, _ = check_line(line)
        for expr, _val, tok, _printed in checks:
            try:
                exact = dec_eval(expr)
            except Exception:                                      # noqa: BLE001
                continue
            scaled = exact * (Decimal(10) ** printed_precision(tok))
            if scaled - scaled.to_integral_value(rounding='ROUND_FLOOR') == Decimal('0.5'):
                return True
    return False


def answer_types() -> dict[str, str]:
    with open(REPO / INVENTORY, encoding='utf-8-sig', newline='') as fh:
        return {r['template_id']: r['answer_type'] for r in csv.DictReader(fh)}


def pilot_questions() -> set[str]:
    with open(PILOT_SLICE, encoding='utf-8') as fh:
        return {json.loads(line)['question'] for line in fh if line.strip()}


def build(master_seed: int):
    """-> (pool records, manifest rows, replacements, short templates). Deterministic."""
    levels, types, pilot = difficulty_map(), answer_types(), pilot_questions()
    records, rows, replaced, short_templates = [], [], [], []
    for ref in discover(None):
        tid = ref.template_id
        short = tid[len('template_'):]
        kept, repeats, rejected, first = [], [], [], {}
        i = 0
        while len(kept) < INSTANCES and i < INSTANCES + MAX_EXTRA:
            seed = item_seed(master_seed, tid, i)
            inst = generate(ref, seed, capture=False)
            if not inst.ok:
                rejected.append((i, 'generation error: ' + str(inst.error)[:100]))
            elif inst.question in pilot:
                rejected.append((i, 'pilot-slice question'))
            elif has_tie(inst.solution):
                rejected.append((i, 'display tie'))
            elif inst.question in first:
                repeats.append((i, seed, inst))
            else:
                first[inst.question] = f'{short}#{i}'
                kept.append((i, seed, inst))
            i += 1
        chosen = list(kept)
        if len(kept) < INSTANCES:
            fill = repeats[:INSTANCES - len(kept)]
            if len(kept) + len(fill) < INSTANCES:
                raise SystemExit(f'{tid}: {len(kept)} distinct and {len(repeats)} repeats '
                                 f'in {i} indices, fewer than {INSTANCES}')
            chosen += fill
            short_templates.append({'template_id': tid, 'distinct': len(kept),
                                    'repeats': len(fill), 'indices_tried': i})
        chosen.sort(key=lambda c: c[0])
        last = chosen[-1][0]
        used = {c[0] for c in chosen}
        skipped = [(j, r) for j, r in rejected if j <= last]
        skipped += [(j, 'duplicate question') for j, _s, _x in repeats if j <= last and j not in used]
        for j, reason in sorted(skipped):
            replaced.append({'template_id': tid, 'instance_index': j, 'reason': reason})
        for j, seed, inst in chosen:
            rec = record(ref, seed, inst.question, inst.solution, levels[tid])
            rec['item_id'] = f'{short}#{j}'
            rec['instance_index'] = j
            row = {
                'item_id': rec['item_id'], 'template_id': tid, 'instance_index': j,
                'branch': rec['branch'], 'domain': rec['domain'], 'area': rec['area'],
                'level': rec['level'], 'answer_type': types.get(tid, 'unknown'),
                'sha256': item_sha(inst.question, inst.solution),
            }
            if first[inst.question] != rec['item_id']:
                rec['repeat_of'] = row['repeat_of'] = first[inst.question]
            records.append(rec)
            rows.append(row)
    return records, rows, replaced, short_templates


def manifest_body(rows) -> str:
    return ''.join(json.dumps(r, ensure_ascii=False, sort_keys=True) + '\n' for r in rows)


def git(*args) -> str:
    return subprocess.run(['git', *args], cwd=REPO, capture_output=True,
                          text=True).stdout.strip()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_seed(create: bool) -> int:
    if SEED_PATH.exists():
        return int(SEED_PATH.read_text(encoding='utf-8').strip())
    if not create:
        raise SystemExit(f'{SEED_PATH.name} not found - restore it from the private backup')
    if MANIFEST.exists():
        raise SystemExit(f'{MANIFEST.name} exists but {SEED_PATH.name} does not: the seed was '
                         'lost. Restore it from the backup; do not draw a new one silently.')
    seed = secrets.randbits(SEED_BITS) | (1 << (SEED_BITS - 1))
    SEED_PATH.write_text(f'{seed}\n', encoding='utf-8')
    return seed


def freeze_doc(master_seed, rows, replaced, short, body) -> dict:
    by_reason = collections.Counter(r['reason'].split(':')[0] for r in replaced)
    return {
        'run': 'full_run_28092026',
        'frozen_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        'commit': git('rev-parse', 'HEAD'),
        'template_inputs_dirty': bool(git('status', '--porcelain', '--', 'data', 'tests',
                                          'generate_testset.py')),
        'script_sha256': {
            'generate_testset.py': file_sha(REPO / 'generate_testset.py'),
            'full_run_28092026/freeze.py': file_sha(Path(__file__)),
        },
        'seed': {
            'bits': SEED_BITS,
            'commitment_sha256': sha(str(master_seed)),
            'held_in': 'full_run_28092026/SEED.secret, not committed; revealed at publication',
        },
        'instances_per_template': INSTANCES,
        'templates': len({r['template_id'] for r in rows}),
        'items': len(rows),
        'distinct_questions': sum('repeat_of' not in r for r in rows),
        'by_branch': dict(sorted(collections.Counter(r['branch'] for r in rows).items())),
        'by_level': dict(sorted(collections.Counter(r['level'] for r in rows).items())),
        'by_answer_type': dict(sorted(collections.Counter(r['answer_type'] for r in rows).items())),
        'short_templates': short,
        'replacements': {'count': len(replaced), 'by_reason': dict(sorted(by_reason.items())),
                         'items': replaced},
        'manifest_sha256': sha(body),
        'pool_sha256': sha(''.join(r['sha256'] for r in rows)),
        'selection_rule': __doc__.split('Selection rule.')[1].strip(),
    }


def check_files() -> int:
    """pool/ on disk against manifest.jsonl. Needs no seed: what a reviewer can run."""
    want = {r['item_id']: r['sha256'] for r in map(json.loads, MANIFEST.read_text(
        encoding='utf-8').splitlines())}
    have = {}
    for path in sorted(POOL_DIR.rglob('*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            have[r['item_id']] = item_sha(r['question'], r['solution'])
    missing = sorted(set(want) - set(have))
    extra = sorted(set(have) - set(want))
    moved = sorted(k for k in set(want) & set(have) if want[k] != have[k])
    if not (missing or extra or moved):
        print(f'FILES OK - {len(have)} items in pool/ match manifest.jsonl')
        return 0
    print(f'FILES DIFFER - missing {len(missing)}, extra {len(extra)}, text moved {len(moved)}')
    for k in (missing + extra + moved)[:20]:
        print(f'  {k}')
    return 1


def verify() -> int:
    """Regenerate from the seed in this process and compare with the committed record."""
    _records, rows, _replaced, _short = build(read_seed(create=False))
    fresh = manifest_body(rows)
    held = MANIFEST.read_text(encoding='utf-8')
    doc = json.loads(FREEZE.read_text(encoding='utf-8'))
    status = 0
    if fresh == held:
        print(f'VERIFY OK - {len(rows)} items regenerate byte-identically to manifest.jsonl')
    else:
        old = {json.loads(x)['item_id']: x for x in held.splitlines()}
        new = {json.loads(x)['item_id']: x for x in fresh.splitlines()}
        diff = sorted(k for k in set(old) | set(new) if old.get(k) != new.get(k))
        print(f'VERIFY FAILED - {len(diff)} items differ from manifest.jsonl')
        for k in diff[:20]:
            print(f'  {k}')
        status = 1
    if doc['manifest_sha256'] != sha(held):
        print('  FREEZE.json records a different manifest sha256 - one of them was edited')
        status = 1
    if doc['seed']['commitment_sha256'] != sha(str(read_seed(create=False))):
        print('  SEED.secret does not match the committed commitment')
        status = 1
    return status or check_files()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--check-files', action='store_true')
    args = ap.parse_args()
    if args.verify:
        return verify()
    if args.check_files:
        return check_files()

    master_seed = read_seed(create=True)
    records, rows, replaced, short = build(master_seed)
    body = manifest_body(rows)
    write(records, str(POOL_DIR))
    MANIFEST.write_text(body, encoding='utf-8', newline='\n')
    doc = freeze_doc(master_seed, rows, replaced, short, body)
    FREEZE.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + '\n', encoding='utf-8',
                      newline='\n')

    print(f"items {doc['items']}  distinct questions {doc['distinct_questions']}  "
          f"templates {doc['templates']}  by level {doc['by_level']}")
    print(f"by branch {doc['by_branch']}")
    print(f"by answer type {doc['by_answer_type']}")
    for t in short:
        print(f"short template {t['template_id']}: {t['distinct']} distinct + {t['repeats']} repeats "
              f"({t['indices_tried']} indices searched)")
    print(f"replacements {doc['replacements']['count']}  {doc['replacements']['by_reason']}")
    for r in replaced:
        print(f"  {r['template_id']} #{r['instance_index']}: {r['reason']}")
    print(f"manifest sha256 {doc['manifest_sha256'][:16]}  seed commitment "
          f"{doc['seed']['commitment_sha256'][:16]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
