"""Cut the frozen pilot slice, and pin it so it cannot move under an annotator.

    python -m evaluator_pilot_17092026.freeze                 # write pilot/pilot_v1/
    python -m evaluator_pilot_17092026.freeze --verify        # rebuild and compare; write nothing
    python -m evaluator_pilot_17092026.freeze --selftest      # the selection rule's own plants

WHY A FREEZE.  Six evaluators get compared on this set, and human experts label
every trace in it.  Both only mean something if the set does not move: a swapped
item invalidates a comparison, and an edited item invalidates the labels that
were paid for.  So the slice is chosen once, by a rule written down here, and
`--verify` fails the moment anything drifts.

WHY THE TEXT AND NOT THE SEED.  A seed does not pin an item across template
versions - that is measured, not feared.  Phase 1 moved 34% of its templates'
instances, Phase 2 moved a template's whole pool, and the C3 corrections moved
824 questions across 13 templates.  `seed 7` means one thing today and another
after the next constants fix.  So the manifest carries the QUESTION and GOLD
TEXT, plus a SHA-256 per item, and FREEZE.json records the commit they came
from.  A later correction is then visible as a verify failure rather than a
silent relabelling of what the experts read.

THE SELECTION RULE, which is the part a reader will want to audit:

  1.  15 cells = 5 branches x 3 difficulty levels.  Every cell is populated
      (the thinnest holds 6 templates), so the slice spans both axes by
      construction rather than by luck.
  2.  Within a cell, candidates exclude EXCLUSIONS below, each with a measured
      reason.  They are all items whose answer is obtainable without doing the
      reasoning, which is exactly what a step-level annotation cannot teach us
      anything about.
  3.  Among the rest, take the template whose ANSWER TYPE is least represented
      in the slice so far.  This is what pulls symbolic, vector and multipart
      answers in rather than leaving a slice of 15 scalars: the comparator's
      known-weak kinds are the ones the pilot most needs to see.
  4.  Ties break on BLAKE2b(master_seed, branch, level, template_id), so the
      choice is deterministic and does not depend on dict or filesystem order.
  5.  4 instances per template, taken in index order but skipping an instance
      whose text a previous one already produced, seeds recomputed with the
      generator's own item_seed - so the slice is a SUBSET of the full testset,
      not a separate generation.  60 problems, and 5 models gives 300 traces.

      The skip is not hypothetical: `incompressible_continuity` emits the same
      question at indices 1 and 2, and the generator reports 38 duplicate
      questions across the full pool.  A duplicate costs twice - an expert
      annotates the same trace twice, and the item is double-weighted in every
      aggregate.  Skipped indices are recorded per template in FREEZE.json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter, OrderedDict

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from generate_testset import DEFAULT_MASTER_SEED, item_seed  # noqa: E402

TESTSET = os.path.join(_ROOT, 'testset')
INVENTORY = os.path.join(_ROOT, 'docs', 're-implementation-sep', 'audit',
                         'template_inventory.csv')
OUT_DIR = os.path.join(_HERE, 'slice')
LEVELS = ('Easy', 'Intermediate', 'Advanced')
INSTANCES = 4

# Measured, each with the record that measured it. These are items whose answer
# can be reached without the reasoning the pilot exists to grade.
EXCLUSIONS = {
    'system_property_linearity':
        'D-057: 100% predictable from the question surface (floor 50.02%)',
    'system_properties_memory_causality':
        'D-057: 100% predictable from the question surface (floor 34.02%)',
    'line_balancing_heuristic':
        'D-046: ~90% shortcuttable by counting task pairs, no heuristic needed',
    'levenspiel_plot_interpretation':
        'D-066: blind-guess floor 1.0000, so the statistic is degenerate',
}


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def load_pool() -> list[dict]:
    if not os.path.isdir(TESTSET):
        raise SystemExit('testset/ not found - run: python generate_testset.py')
    rows = []
    for root, _dirs, files in os.walk(TESTSET):
        for fn in sorted(files):
            if fn.endswith('.jsonl'):
                with open(os.path.join(root, fn), encoding='utf-8') as fh:
                    rows.extend(json.loads(ln) for ln in fh)
    return rows


def answer_types() -> dict[str, str]:
    with open(INVENTORY, encoding='utf-8') as fh:
        return {r['template_id'][len('template_'):]: r['answer_type']
                for r in csv.DictReader(fh)}


def select(pool: list[dict], types: dict[str, str], master_seed: int):
    """The rule in the docstring. Returns (chosen, audit)."""
    by_cell: dict[tuple, set] = {}
    for r in pool:
        by_cell.setdefault((r['branch'], r['level']), set()).add(r['id'])

    chosen, audit = OrderedDict(), []
    seen: Counter = Counter()
    for branch in sorted({b for b, _ in by_cell}):
        for level in LEVELS:
            cell = sorted(by_cell.get((branch, level), ()))
            eligible = [t for t in cell if t not in EXCLUSIONS]
            if not eligible:
                raise SystemExit('cell %s/%s has no eligible template' % (branch, level))

            def key(t):
                h = hashlib.blake2b(
                    ('%d|%s|%s|%s' % (master_seed, branch, level, t)).encode(),
                    digest_size=8).hexdigest()
                return (seen[types.get(t, 'unknown')], h)

            pick = min(eligible, key=key)
            seen[types.get(pick, 'unknown')] += 1
            chosen[pick] = (branch, level)
            audit.append({
                'branch': branch, 'level': level, 'chosen': pick,
                'answer_type': types.get(pick, 'unknown'),
                'candidates': len(cell),
                'excluded': sorted(set(cell) & set(EXCLUSIONS)),
            })
    return chosen, audit


def build(master_seed: int):
    pool = load_pool()
    types = answer_types()
    chosen, audit = select(pool, types, master_seed)

    by_id: dict[str, dict] = {}
    for r in pool:
        by_id[(r['id'], r['seed'])] = r

    pool_instances = max(Counter(r['id'] for r in pool).values())
    manifest, skipped = [], {}
    for tid, (branch, level) in chosen.items():
        taken, seen_text, i = 0, set(), 0
        while taken < INSTANCES:
            if i >= pool_instances:
                raise SystemExit('%s: only %d distinct instances in the pool, need %d'
                                 % (tid, taken, INSTANCES))
            seed = item_seed(master_seed, 'template_' + tid, i)
            rec = by_id.get((tid, seed))
            if rec is None:
                raise SystemExit('%s index %d (seed %d) is not in testset/ - '
                                 'was the pool built at this master seed?' % (tid, i, seed))
            sha = _sha(rec['question'] + '\x00' + rec['solution'])
            if sha in seen_text:
                # Same text as an instance already taken: annotating it twice
                # buys nothing and double-weights the item.
                skipped.setdefault(tid, []).append(i)
                i += 1
                continue
            seen_text.add(sha)
            taken += 1
            manifest.append({
                'item_id': '%s#%d' % (tid, i),
                'template_id': 'template_' + tid,
                'branch': rec['branch'], 'domain': rec['domain'], 'area': rec['area'],
                'level': rec['level'], 'answer_type': types.get(tid, 'unknown'),
                'unit': rec['unit'], 'seed': seed, 'instance_index': i,
                'question': rec['question'], 'solution': rec['solution'],
                'sha256': sha,
            })
            i += 1
    for row in audit:
        row['skipped_duplicate_indices'] = skipped.get(row['chosen'], [])
    return manifest, audit


def _git(*args):
    try:
        return subprocess.run(['git'] + list(args), cwd=_ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return 'unknown'


def freeze_doc(manifest, audit, master_seed, manifest_sha):
    return {
        'pilot': 'evaluator_pilot_17092026',
        'frozen_at_utc': __import__('datetime').datetime.now(
            __import__('datetime').timezone.utc).isoformat(timespec='seconds'),
        'commit': _git('rev-parse', 'HEAD'),
        'tree_dirty': bool(_git('status', '--porcelain')),
        'master_seed': master_seed,
        'instances_per_template': INSTANCES,
        'templates': len(audit),
        'items': len(manifest),
        'traces_at_5_models': len(manifest) * 5,
        'selection_rule': [ln.strip() for ln in __doc__.split('THE SELECTION RULE')[1].splitlines() if ln.strip()],
        'exclusions': EXCLUSIONS,
        'cells': audit,
        'by_branch': dict(Counter(r['branch'] for r in manifest)),
        'by_level': dict(Counter(r['level'] for r in manifest)),
        'by_answer_type': dict(Counter(r['answer_type'] for r in manifest)),
        'manifest_sha256': manifest_sha,
    }


def write(manifest, audit, master_seed):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'manifest.jsonl')
    body = ''.join(json.dumps(r, ensure_ascii=False, sort_keys=True) + '\n' for r in manifest)
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(body)
    doc = freeze_doc(manifest, audit, master_seed, _sha(body))
    with open(os.path.join(OUT_DIR, 'FREEZE.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(doc, fh, indent=2, ensure_ascii=False)
        fh.write('\n')
    return doc


def report(doc, manifest):
    print('slice: %d templates x %d instances = %d items, %d traces at 5 models'
          % (doc['templates'], INSTANCES, doc['items'], doc['traces_at_5_models']))
    print('  by branch      %s' % doc['by_branch'])
    print('  by level       %s' % doc['by_level'])
    print('  by answer type %s' % doc['by_answer_type'])
    print('  commit %s%s' % (doc['commit'][:12], '  (TREE DIRTY)' if doc['tree_dirty'] else ''))
    print('  manifest sha256 %s' % doc['manifest_sha256'][:16])
    for c in doc['cells']:
        print('    %-12s %-12s %-42s %-14s (%d candidates%s)'
              % (c['branch'].replace('_engineering', ''), c['level'], c['chosen'],
                 c['answer_type'], c['candidates'],
                 ', %d excluded' % len(c['excluded']) if c['excluded'] else ''))


def verify(master_seed):
    path = os.path.join(OUT_DIR, 'manifest.jsonl')
    if not os.path.exists(path):
        print('no frozen manifest at %s' % path)
        return 1
    manifest, audit = build(master_seed)
    fresh = ''.join(json.dumps(r, ensure_ascii=False, sort_keys=True) + '\n' for r in manifest)
    with open(path, encoding='utf-8') as fh:
        held = fh.read()
    with open(os.path.join(OUT_DIR, 'FREEZE.json'), encoding='utf-8') as fh:
        doc = json.load(fh)

    if held == fresh:
        print('VERIFY OK - %d items, manifest byte-identical to a fresh build' % len(manifest))
        if doc['manifest_sha256'] != _sha(held):
            print('  but FREEZE.json records a different sha256 - the doc was edited by hand')
            return 1
        return 0

    old = {json.loads(ln)['item_id']: json.loads(ln) for ln in held.splitlines()}
    new = {r['item_id']: r for r in manifest}
    print('VERIFY FAILED')
    for k in sorted(set(old) - set(new)):
        print('  gone from the slice : %s' % k)
    for k in sorted(set(new) - set(old)):
        print('  new in the slice    : %s' % k)
    for k in sorted(set(old) & set(new)):
        if old[k]['sha256'] != new[k]['sha256']:
            print('  TEXT MOVED          : %s  (the experts read the old text)' % k)
    return 1


def selftest():
    """The rule must be deterministic, exclusive, and spread across both axes."""
    pool, types = load_pool(), answer_types()
    bad = []
    a, _ = select(pool, types, DEFAULT_MASTER_SEED)
    b, _ = select(pool, types, DEFAULT_MASTER_SEED)
    if list(a) != list(b):
        bad.append('selection is not deterministic at one master seed')
    c, _ = select(pool, types, DEFAULT_MASTER_SEED + 1)
    if list(a) == list(c):
        bad.append('selection ignores the master seed')
    if set(a) & set(EXCLUSIONS):
        bad.append('an excluded template was selected: %s' % (set(a) & set(EXCLUSIONS)))
    cells = Counter((br, lv) for br, lv in a.values())
    if len(cells) != 15 or set(cells.values()) != {1}:
        bad.append('not exactly one template per branch x level cell: %s' % cells)
    kinds = Counter(types.get(t, 'unknown') for t in a)
    if len(kinds) < 4:
        bad.append('answer-type spread collapsed to %s' % dict(kinds))
    manifest, _ = build(DEFAULT_MASTER_SEED)
    if len(manifest) != len(a) * INSTANCES:
        bad.append('manifest is %d items, expected %d' % (len(manifest), len(a) * INSTANCES))
    if len({r['sha256'] for r in manifest}) != len(manifest):
        bad.append('two items carry the same text hash')
    for label in ('deterministic', 'seed-sensitive', 'exclusions honoured',
                  '15 cells', 'answer-type spread', 'manifest shape'):
        print('  [%s] %s' % ('FAIL' if any(label.split()[0] in x for x in bad) else 'ok', label))
    for x in bad:
        print('  - ' + x)
    print('selftest: %d failure(s)' % len(bad))
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--master-seed', type=int, default=DEFAULT_MASTER_SEED)
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if args.verify:
        return verify(args.master_seed)
    manifest, audit = build(args.master_seed)
    doc = write(manifest, audit, args.master_seed)
    report(doc, manifest)
    print('\nwrote %s' % os.path.relpath(OUT_DIR, _ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
