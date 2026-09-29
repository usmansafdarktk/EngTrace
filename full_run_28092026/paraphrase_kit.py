"""The expert check of the paraphrases (ANALYSIS_PLAN Q5, D-141): one own-branch expert confirms that
each paraphrase is the same problem with the same answer, and a pair the expert rejects leaves both
arms of Q5.

    python -m full_run_28092026.paraphrase_kit --build          # FREE: tasks and kits from paraphrase/pool.jsonl
    python -m full_run_28092026.paraphrase_kit --score DIR      # FREE: the experts' returned <id>.jsonl files
    python -m full_run_28092026.paraphrase_kit --selftest       # FREE: the whole path on a stand-in pool, in a temp folder

ASSIGNMENT. Each branch's templates, sorted, are dealt to its three experts in turn - the 1st, 4th,
7th ... to expert 1 - so an expert judges whole templates: 10 templates, up to 30 items, in an order
shuffled per expert from a fixed seed. The experts are the layer-2 roster (the pilot's 15 unless
template_annotation_23092026/layer2/annotators.json overrides it); a kit carries only their id.

WHAT AN EXPERT SEES, as plain text, not rendered: the original question, the paraphrase, and the
original's gold answer, with its full solution one click away. Three questions:
  same     Is it the same problem: the same givens, the same quantities asked for, the same conditions?
  answer   Does the original's answer still answer it, exactly?
  clear    Does it add an ambiguity, an error or a hint that the original does not have?
A pair is kept when `same` and `answer` are yes and `clear` is no. Any other answer needs a note.

WHAT IS KEPT. paraphrase/tasks/ (with the keyfile from codes to items), paraphrase/dist/ (the kits)
and the returns hold pool text or the experts' labels and stay local and gitignored. --score writes
paraphrase/accepted.json, local, which analyze.py reads to drop rejected pairs, and
PARAPHRASE_REVIEW.md, committed, counts only.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import random
import shutil
import statistics
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

OUT = HERE / 'paraphrase'
SALT = 'engtrace-paraphrase-review-2026-09'
SEED = 20260929
QUESTIONS = ('same', 'answer', 'clear')
README = """EngTrace paraphrase check - how to run

You should have, side by side in one folder: app.py, this README, guide.md, and a folder named
kit_<your id> holding your items.

1. Read guide.md (one page).

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id in the sidebar.

3. Work through your items. The app saves after every item. To pause, close the browser tab and
   stop the app; run "streamlit run app.py" again to continue where you left off.

4. When the app says you are done, send the file <your id>.jsonl, which the app has written in
   this folder next to app.py, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers.
"""


def roster() -> list[dict]:
    from template_annotation_23092026.layer2.build_tasks import roster as layer2_roster
    return [{'id': r['id'], 'branch': r['branch']} for r in layer2_roster()]


def code_of(item_id: str) -> str:
    return 'PR-' + hashlib.sha256(f'{SALT}|{item_id}'.encode('utf-8')).hexdigest()[:8]


def build(out: Path, paraphrases: dict[str, str], items: dict[str, dict], experts: list[dict]) -> dict:
    """tasks/ and dist/ under `out` for the paraphrases given ({item_id: text})."""
    import answer
    by_branch = collections.defaultdict(list)
    for e in experts:
        by_branch[e['branch']].append(e['id'])
    pool, assignment, keyfile = {}, {}, {}
    rng = random.Random(SEED)
    for branch, ids in sorted(by_branch.items()):
        ids = sorted(ids)
        templates = sorted({items[i]['template_id'] for i in paraphrases if items[i]['branch'] == branch})
        for k, aid in enumerate(ids):
            mine = [t for j, t in enumerate(templates) if j % len(ids) == k]
            codes = []
            for i in sorted(paraphrases):
                it = items[i]
                if it['template_id'] not in mine:
                    continue
                c = code_of(i)
                keyfile[c] = i
                pool[c] = {'branch': branch, 'area': it.get('area', ''), 'original': it['question'],
                           'paraphrase': paraphrases[i], 'gold_answer': answer.segment(it['solution']).strip(),
                           'solution': it['solution']}
                codes.append(c)
            rng.shuffle(codes)
            assignment[aid] = {'branch': branch, 'codes': codes}
    tasks, dist = out / 'tasks', out / 'dist'
    for d in (tasks, dist):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    (tasks / 'pool.json').write_text(json.dumps(pool, ensure_ascii=False), encoding='utf-8')
    (tasks / 'assignment.json').write_text(json.dumps(assignment, indent=1), encoding='utf-8')
    (tasks / 'keyfile.json').write_text(json.dumps(keyfile, indent=1), encoding='utf-8')
    shutil.copy(HERE / 'paraphrase_app.py', dist / 'app.py')
    shutil.copy(HERE / 'paraphrase_guide.md', dist / 'guide.md')
    (dist / 'README.txt').write_text(README, encoding='utf-8')
    for aid, mine in assignment.items():
        kit = dist / f'kit_{aid}' / 'tasks'
        kit.mkdir(parents=True)
        (kit / 'pool.json').write_text(json.dumps({c: pool[c] for c in mine['codes']}, ensure_ascii=False),
                                       encoding='utf-8')
        (kit / 'assignment.json').write_text(json.dumps({aid: mine}, indent=1), encoding='utf-8')
    return {'items': len(pool), 'experts': {a: len(m['codes']) for a, m in assignment.items()}}


def kept(row: dict) -> bool:
    return row['same'] == 'yes' and row['answer'] == 'yes' and row['clear'] == 'no'


def score(returned: Path, out: Path, summary: Path | None) -> dict:
    """accepted.json and the summary from the returned files; the last submission per item counts."""
    keyfile = json.loads((out / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    assignment = json.loads((out / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
    owner = {c: a for a, m in assignment.items() for c in m['codes']}
    rows = {}
    for f in sorted(returned.glob('*.jsonl')):
        for line in f.read_text(encoding='utf-8').splitlines():
            if line.strip():
                r = json.loads(line)
                if r['code'] not in keyfile:
                    raise SystemExit(f"{f.name}: code {r['code']} is not in this build's keyfile")
                if owner[r['code']] != r['annotator_id']:
                    raise SystemExit(f"{f.name}: {r['code']} was not assigned to {r['annotator_id']}")
                rows[r['code']] = r
    accepted = {keyfile[c]: {'kept': kept(r), **{q: r[q] for q in QUESTIONS}, 'annotator_id': r['annotator_id']}
                for c, r in rows.items()}
    (out / 'accepted.json').write_text(json.dumps(accepted, indent=1) + '\n', encoding='utf-8')
    branch_of = {c: assignment[a]['branch'] for c, a in owner.items()}
    per_branch = collections.defaultdict(collections.Counter)
    for c, r in rows.items():
        per_branch[branch_of[c]]['kept' if kept(r) else 'rejected'] += 1
    reasons = collections.Counter()
    for r in rows.values():
        reasons.update(f'{q}={r[q]}' for q, bad in (('same', 'no'), ('answer', 'no'), ('clear', 'yes')) if r[q] == bad)
    import datetime as dt
    dwell = [(dt.datetime.fromisoformat(r['submitted_at']) - dt.datetime.fromisoformat(r['opened_at'])).total_seconds()
             for r in rows.values() if r.get('opened_at') and r.get('submitted_at')]
    res = {'assigned': len(keyfile), 'returned': len(rows), 'kept': sum(kept(r) for r in rows.values()),
           'missing': len(keyfile) - len(rows)}
    L = ['# The expert check of the paraphrases (D-141)', '',
         'Generated by `paraphrase_kit.py --score`; the questions and the rule are in its docstring. Counts '
         'only: the experts\' files and the verdict per item stay local.', '',
         '| | |', '|---|---|',
         f"| paraphrases assigned | {res['assigned']} |", f"| returned | {res['returned']} |",
         f"| kept | {res['kept']} |", f"| rejected | {res['returned'] - res['kept']} |",
         f"| not returned | {res['missing']} |",
         '| kept / rejected, per branch | ' + '; '.join(f"{b.replace('_engineering', '')} {c['kept']} / {c['rejected']}"
                                                  for b, c in sorted(per_branch.items())) + ' |',
         '| answers behind the rejections | ' + (', '.join(f'{k} {v}' for k, v in reasons.most_common()) or 'none') + ' |',
         f"| median seconds per item | {statistics.median(dwell):.0f} |" if dwell else '| median seconds per item | - |']
    if summary:
        summary.write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return res


def real_inputs() -> tuple[dict, dict]:
    from full_run_28092026 import score as S
    manifest = {r['item_id']: r for r in map(json.loads, (OUT / 'manifest.jsonl').read_text(encoding='utf-8').splitlines())}
    paras = {}
    for r in map(json.loads, (OUT / 'pool.jsonl').read_text(encoding='utf-8').splitlines()):
        m = manifest[r['item_id']]
        if m['passed'] and hashlib.sha256(r['question'].encode('utf-8')).hexdigest() == m['sha256']:
            paras[r['item_id']] = r['question']
    return paras, S.pool_items()


def selftest() -> int:
    """Build, fill and score kits for a stand-in pool (the originals with a marker line) in a temp folder."""
    from full_run_28092026 import score as S, subsamples
    items = S.pool_items()
    paras = {i: 'STAND-IN PARAPHRASE\n' + items[i]['question'] for i in subsamples.paraphrase_ids()}
    experts = roster()
    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        info = build(out, paras, items, experts)
        assignment = json.loads((out / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
        keyfile = json.loads((out / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
        codes = [c for m in assignment.values() for c in m['codes']]
        if sorted(codes) != sorted(keyfile) or len(codes) != len(paras):
            bad.append('not every paraphrase assigned exactly once')
        branch_of = {e['id']: e['branch'] for e in experts}
        if any(items[keyfile[c]]['branch'] != branch_of[a] for a, m in assignment.items() for c in m['codes']):
            bad.append('an item went to another branch')
        per_template = collections.defaultdict(set)
        for a, m in assignment.items():
            for c in m['codes']:
                per_template[items[keyfile[c]]['template_id']].add(a)
        if any(len(v) != 1 for v in per_template.values()):
            bad.append('a template is split between experts')
        if set(info['experts'].values()) != {30}:
            bad.append(f"kit sizes {sorted(set(info['experts'].values()))}, not 30 each")
        for kit in (out / 'dist').glob('kit_*'):
            blob = (kit / 'tasks' / 'pool.json').read_text(encoding='utf-8')
            if 'keyfile' in blob or any(i in blob for i in list(keyfile.values())[:50]):
                bad.append(f'{kit.name} reveals item ids')
        rng = random.Random(1)
        ret = out / 'returned'
        ret.mkdir()
        want_kept = 0
        for a, m in assignment.items():
            with open(ret / f'{a}.jsonl', 'w', encoding='utf-8') as fh:
                for c in m['codes']:
                    reject = rng.random() < 0.1
                    want_kept += not reject
                    fh.write(json.dumps({'annotator_id': a, 'code': c, 'same': 'no' if reject else 'yes',
                                         'answer': 'yes', 'clear': 'no', 'note': 'x' if reject else '',
                                         'opened_at': '2026-10-01T10:00:00+00:00',
                                         'submitted_at': '2026-10-01T10:01:30+00:00'}) + '\n')
        res = score(ret, out, None)
        if res['kept'] != want_kept or res['missing'] != 0:
            bad.append(f'scored {res}, simulated {want_kept} kept')
        acc = json.loads((out / 'accepted.json').read_text(encoding='utf-8'))
        if set(acc) != set(paras):
            bad.append('accepted.json does not cover every returned item')
    print('selftest:', 'all pass' if not bad else bad)
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--score', metavar='DIR')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.build:
        paras, items = real_inputs()
        info = build(OUT, paras, items, roster())
        print(f"{info['items']} paraphrases in {len(info['experts'])} kits: "
              + ', '.join(f'{k} {v}' for k, v in sorted(info['experts'].items())))
        print(f'kits in {OUT / "dist"}; send each expert the shared files and their kit_<id>/ folder')
        return 0
    if a.score:
        score(Path(a.score), OUT, HERE / 'PARAPHRASE_REVIEW.md')
        return 0
    ap.print_help()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
