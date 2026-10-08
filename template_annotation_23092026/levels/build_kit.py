"""E1 (WS-E): the experts' independent difficulty ratings, one kit per branch.

    python -m template_annotation_23092026.levels.build_kit              # FREE: tasks/ and dist/<branch>/ (local),
                                                                        #       RESULTS.md with what was sent (counts)
    python -m template_annotation_23092026.levels.build_kit --selftest   # FREE: a build in a temp folder, the app driven

WHAT AN EXPERT SEES. The 30 templates of their branch, in an order of their own (shuffled from the seed and their id),
under opaque codes that the three experts of a branch share. Per template: its domain and area, one problem of the
evaluation set (full_run_28092026/pool/, the frozen 2,250) with its reference solution, and the rubric of the paper's
Section 3.1 (guide.md, and in the app's sidebar). The instance is drawn from the seed and the template id, so the
three experts of a branch see the same problem. Two questions: the level (Easy, Intermediate, Advanced), and whether
the problem states the governing formula or names the method. The template's current label is not in the kit: the
build refuses a payload that holds a template id or a level word.

WHAT IS KEPT LOCAL (levels/.gitignore). tasks/ (pool.json; assignment.json; keyfile.json: code -> template, branch,
current level, the instance shown) and dist/ hold evaluation-set text; experts_filled_levels/ holds the returns;
agreement.json carries the ratings per template. Committed: this script, app.py, guide.md, score_levels.py,
apply_labels.py and RESULTS.md (counts).

DELIVERY. dist/<branch>/ holds app.py, README.txt, guide.md and kit_<id>/tasks/ for each of the branch's three
experts. Send the branch folder to its three experts (each picks their id in the app), or each expert the three shared
files and their own kit_<id>/. The returned <id>.jsonl files go into experts_filled_levels/, in any layout.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import random
import re
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

from full_run_28092026 import score  # noqa: E402

SEED = 20261007
SALT = 'engtrace-levels-2026-10'
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
RESULTS = HERE / 'RESULTS.md'
LEVELS = ['Easy', 'Intermediate', 'Advanced']
LEVEL_WORD = re.compile(r'\b(?:Easy|Intermediate|Advanced)\b')

README = """EngTrace difficulty ratings, October 2026

This folder holds the task for the experts of one branch: app.py, this README, guide.md, and one kit_<id> folder
per expert.

1. Read guide.md first (two pages).

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id in the sidebar.

3. Rate the 30 templates. The app saves after every item. To pause, close the browser tab and stop the app; run
   "streamlit run app.py" again to continue where you left off.

4. When the app says you are done, send the file <your id>.jsonl, which the app has written in this folder next to
   app.py, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers. Please do not share the folder: the
problems are the benchmark's private test set until the paper is published.
"""


def experts() -> list[dict]:
    """The layer-2 roster, three experts per branch; only the id and the branch are read."""
    from template_annotation_23092026.layer2.build_tasks import roster
    return [{'id': r['id'], 'branch': r['branch']} for r in roster()]


def code_of(template_id: str) -> str:
    return 'LV-' + hashlib.sha256(f'{SALT}|{template_id}'.encode('utf-8')).hexdigest()[:6]


def manifest_sha() -> str:
    return hashlib.sha256((score.HERE / 'manifest.jsonl').read_bytes()).hexdigest()


def draw(items: dict, seed: int) -> list[dict]:
    """One instance per template, drawn from the seed and the template id; the payload and the key kept apart."""
    by_t = collections.defaultdict(list)
    for it in items.values():
        by_t[it['template_id']].append(it)
    out = []
    for tid in sorted(by_t):
        its = sorted(by_t[tid], key=lambda it: it['instance_index'])
        it = random.Random(f'{seed}|{tid}').choice(its)
        out.append({'code': code_of(tid),
                    'pool': {'branch': it['branch'], 'domain': it['domain'], 'area': it['area'],
                             'question': it['question'], 'solution': it['solution']},
                    'key': {'template_id': tid, 'branch': it['branch'], 'domain': it['domain'], 'area': it['area'],
                            'level': it['level'], 'item_id': it['item_id']}})
    return out


def guard(tasks: list[dict]) -> None:
    """Refuse a payload that names its template or holds a level word."""
    for t in tasks:
        blob = json.dumps(t['pool'], ensure_ascii=False)
        tid = t['key']['template_id']
        if tid in blob or tid.replace('template_', '') in blob:
            raise SystemExit(f"refusing: {t['code']}'s payload names its template")
        if LEVEL_WORD.search(blob):
            raise SystemExit(f"refusing: {t['code']}'s payload holds a level word: "
                             f"{LEVEL_WORD.search(blob).group(0)!r}")


def assign(tasks: list[dict], roster: list[dict], seed: int) -> dict:
    """Every expert gets their branch's templates, in an order shuffled from the seed and their id."""
    by_branch = collections.defaultdict(list)
    for t in tasks:
        by_branch[t['pool']['branch']].append(t['code'])
    queues = {}
    for e in roster:
        codes = sorted(by_branch[e['branch']])
        if not codes:
            raise SystemExit(f"no templates for branch {e['branch']}")
        random.Random(f'{seed}|{e["id"]}').shuffle(codes)
        queues[e['id']] = {'branch': e['branch'], 'codes': codes}
    return queues


def build(tasks_dir: Path, dist_dir: Path, seed: int, results: Path | None) -> dict:
    items = score.pool_items()
    tasks = draw(items, seed)
    codes = [t['code'] for t in tasks]
    if len(set(codes)) != len(codes):
        raise SystemExit('code collision')
    guard(tasks)
    roster = experts()
    queues = assign(tasks, roster, seed)
    pool = {t['code']: t['pool'] for t in tasks}
    keyfile = {t['code']: t['key'] for t in tasks}
    for d in (tasks_dir, dist_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    (tasks_dir / 'pool.json').write_text(json.dumps(pool, ensure_ascii=False), encoding='utf-8')
    (tasks_dir / 'assignment.json').write_text(json.dumps(queues, indent=1), encoding='utf-8')
    (tasks_dir / 'keyfile.json').write_text(json.dumps(keyfile, indent=1), encoding='utf-8')
    for branch in sorted({q['branch'] for q in queues.values()}):
        folder = dist_dir / branch
        folder.mkdir()
        shutil.copy(HERE / 'app.py', folder / 'app.py')
        shutil.copy(HERE / 'guide.md', folder / 'guide.md')
        (folder / 'README.txt').write_text(README, encoding='utf-8')
        for aid, q in sorted(queues.items()):
            if q['branch'] != branch:
                continue
            kit = folder / f'kit_{aid}' / 'tasks'
            kit.mkdir(parents=True)
            (kit / 'pool.json').write_text(json.dumps({c: pool[c] for c in q['codes']}, ensure_ascii=False),
                                           encoding='utf-8')
            (kit / 'assignment.json').write_text(json.dumps({aid: q}, indent=1), encoding='utf-8')
    comp = composition(tasks, queues, seed)
    (tasks_dir / 'build.json').write_text(json.dumps(comp, indent=1), encoding='utf-8')
    if results:
        results.write_text('\n'.join(sent_lines(comp) + ['## What came back', '', 'Nothing yet.', '']),
                           encoding='utf-8', newline='\n')
    return comp


def composition(tasks: list[dict], queues: dict, seed: int) -> dict:
    by_branch = collections.defaultdict(collections.Counter)
    for t in tasks:
        by_branch[t['key']['branch']][t['key']['level']] += 1
    return {'built': dt.date.today().isoformat(), 'seed': seed, 'manifest_sha256': manifest_sha(),
            'templates': len(tasks), 'experts': len(queues), 'readings': sum(len(q['codes']) for q in queues.values()),
            'branches': {b: {'templates': sum(c.values()), 'current_levels': {lv: c[lv] for lv in LEVELS},
                             'experts': sorted(a for a, q in queues.items() if q['branch'] == b)}
                         for b, c in sorted(by_branch.items())}}


def sent_lines(comp: dict) -> list[str]:
    L = ['# Difficulty ratings (E1): agreement on the levels', '',
         'Generated by `build_kit.py` (what was sent) and `score_levels.py` (what came back). Counts only; the kits, '
         'the keyfile and the experts\' files stay local. The machine-readable results are in `agreement.json`.', '',
         '## What was sent', '',
         f"Built {comp['built']} from the frozen pool (manifest sha256 `{comp['manifest_sha256'][:12]}`), seed "
         f"{comp['seed']}. Each expert rates the 30 templates of their branch, one evaluation-set problem and its "
         'reference solution per template (the same problem for the three experts of a branch), in an order of their '
         'own, under opaque codes, with the rubric of Section 3.1; the current label is not shown. Two questions per '
         'template: the level, and whether the problem states the governing formula or names the method.', '',
         '| branch | templates | experts | current labels (Easy / Intermediate / Advanced) |', '|---|---:|---:|---|']
    for b, v in comp['branches'].items():
        cl = v['current_levels']
        L.append(f"| {b.replace('_engineering', '')} | {v['templates']} | {len(v['experts'])} | "
                 f"{cl['Easy']} / {cl['Intermediate']} / {cl['Advanced']} |")
    L += ['', f"Templates {comp['templates']}, experts {comp['experts']}, ratings asked {comp['readings']}.", '']
    return L


# ------------------------------------------------------------------ selftest

def selftest() -> int:
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-levels-kit-'))
    try:
        comp = build(tmp / 'tasks', tmp / 'dist', SEED, None)
        queues = json.loads((tmp / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
        keyfile = json.loads((tmp / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
        assert comp['templates'] == 150 and len(keyfile) == 150, comp['templates']
        assert comp['experts'] == 15 and comp['readings'] == 450, (comp['experts'], comp['readings'])
        by_branch = collections.defaultdict(list)
        for aid, q in queues.items():
            assert len(q['codes']) == 30, (aid, len(q['codes']))
            assert all(keyfile[c]['branch'] == q['branch'] for c in q['codes']), aid
            by_branch[q['branch']].append(q['codes'])
        for b, lists in by_branch.items():
            assert len(lists) == 3 and len({frozenset(x) for x in lists}) == 1, b       # the same 30 for the three
            assert len({tuple(x) for x in lists}) == 3, b                                 # in three different orders
        again = draw(score.pool_items(), SEED)
        assert {t['code']: t['key']['item_id'] for t in again} == {c: k['item_id'] for c, k in keyfile.items()}
        for f in (tmp / 'dist').rglob('pool.json'):
            blob = f.read_text(encoding='utf-8')
            assert not LEVEL_WORD.search(blob) and 'template_' not in blob, f
        print(f"built {comp['templates']} templates for {comp['experts']} experts; payloads clean; the three experts "
              'of each branch share 30 codes in three orders; the instance draw is deterministic')
        from streamlit.testing.v1 import AppTest
        folder = tmp / 'dist' / 'electrical_engineering'
        aid = 'ele-2'
        at = AppTest.from_file(str(folder / 'app.py'), default_timeout=120)
        at.run()
        assert not at.exception, at.exception
        at.sidebar.selectbox[0].select(aid).run()
        assert not at.exception, at.exception
        first = queues[aid]['codes'][0]
        assert any(first in w.value for w in at.subheader), [w.value for w in at.subheader]
        at.radio[0].set_value('Advanced').run()
        at.radio[1].set_value('no').run()
        at.button[0].click().run()
        assert not at.exception, at.exception
        row = json.loads((folder / f'{aid}.jsonl').read_text(encoding='utf-8').splitlines()[0])
        assert row['code'] == first and row['level'] == 'Advanced' and row['formula_stated'] == 'no', row
        assert any(queues[aid]['codes'][1] in w.value for w in at.subheader), 'the app did not move on'
        print('the app renders, takes a rating and moves to the next template')
        print('SELFTEST OK')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--seed', type=int, default=SEED)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    comp = build(TASKS, DIST, a.seed, RESULTS)
    print('\n'.join(sent_lines(comp)))
    print(f'kits: {DIST}  (one folder per branch; send it to the branch\'s three experts)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
