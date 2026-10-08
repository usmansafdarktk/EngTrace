"""E2 (WS-E): the experts' grading of symbolic final answers, one kit per branch that has symbolic templates.

    python -m full_run_28092026.symbolic.build_kit               # FREE: tasks/ and dist/<branch>/ (local), GRADES.md
    python -m full_run_28092026.symbolic.build_kit --readers 2   # every item read twice (rotating pairs)
    python -m full_run_28092026.symbolic.build_kit --selftest    # FREE: a build in a temp folder, the app driven

WHAT IS DRAWN, from the main store (scores/main, the eleven roster models) and the traces it scored
(score.texts_matching refuses a trace the store did not score): every answered, usable row on the nine templates
whose answer_type is symbolic, by the answer check's verdict (the by-numbers rule):
  incorrect   all of them: RESIDUAL_INCORRECT.md's symbolic column, derived here from the store
  partial     all of them: an equivalence check would act on these as well, so it needs their grades too
  correct     --controls per template (default 4), round-robin over the models from the seed: they keep the grader
              blind to the verdict mix, and they measure the by-numbers rule's precision on these templates
The expert sees the problem, the reference answer (the gold's answer segment), and the model's final answer as the
check reads it (answer.segment of the trace, the very reader the score used), both in full one click away; for the
one template whose reference is a single value written as an expression (ber_estimation_mary: c * Q(a)), that value.
Not the model, not the check's verdict, not the template's name: codes are opaque, templates are "problem families"
A, B, ... in a seeded order.

WHO GRADES. The branch's three experts (the layer-2 roster). Each item goes to one of them, in turn within each
family; --overlap of each family's items (default 0.2) also goes to a second, so the grades' agreement can be
measured; --readers 2 gives every item to two (rotating pairs). A queue runs family by family; within a family the
items are shuffled per expert.

WHAT IS KEPT LOCAL (symbolic/.gitignore). tasks/ (pool.json, assignment.json, keyfile.json: code -> item, model,
template, the check's verdict) and dist/ hold pool text and responses; experts_filled_symbolic/ holds the returns;
grades.json and validation.json carry a grade per item. Committed: the scripts, app.py, guide.md and the count
reports (GRADES.md, SYMBOLIC_CHECK.md, RESCORE_PREVIEW.md).

DELIVERY. dist/<branch>/ holds app.py, README.txt, guide.md and kit_<id>/tasks/ for each of the branch's experts.
The returned <id>.jsonl files go into experts_filled_symbolic/, in any layout; score_grades.py reads them.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
import random
import re
import shutil
import string
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULL = HERE.parent
REPO = FULL.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

SEED = 20261007
SALT = 'engtrace-symbolic-2026-10'
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
GRADES_MD = HERE / 'GRADES.md'
GROUPS = ['incorrect', 'partial', 'control']
BER = 'template_ber_estimation_mary'
BER_LINE = re.compile(r'BER approx\s*([\d.]+)\s*\*\s*Q\(\s*([\d.]+)\s*\)')

README = """EngTrace symbolic grading, October 2026

This folder holds the task for the experts of one branch: app.py, this README, guide.md, and one kit_<id> folder
per expert.

1. Read guide.md first (two pages).

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id in the sidebar.

3. Grade your items. The app saves after every item. To pause, close the browser tab and stop the app; run
   "streamlit run app.py" again to continue where you left off.

4. When the app says you are done, send the file <your id>.jsonl, which the app has written in this folder next to
   app.py, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers. Please do not share the folder: the
problems are the benchmark's private test set until the paper is published.
"""


def experts() -> list[dict]:
    from template_annotation_23092026.layer2.build_tasks import roster
    return [{'id': r['id'], 'branch': r['branch']} for r in roster()]


def code_of(model: str, item_id: str) -> str:
    return 'SG-' + hashlib.sha256(f'{SALT}|{model}|{item_id}'.encode('utf-8')).hexdigest()[:8]


def q_value(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def reference_value(template_id: str, reference: str) -> str | None:
    """The value of a reference that is one value written as an expression (BER = c * Q(a))."""
    if template_id != BER:
        return None
    m = BER_LINE.search(reference)
    if not m:
        return None
    c, a = float(m.group(1)), float(m.group(2))
    return f'{m.group(1)} · Q({m.group(2)}) = {c * q_value(a):.3g}'


def usable(r: dict) -> bool:
    return r['status'] == 'answered' and not r['unusable']


def symbolic_templates(items: dict) -> list[str]:
    return sorted({it['template_id'] for it in items.values() if it['answer_type'] == 'symbolic'})


def candidates(models: list[str]) -> tuple[dict, dict, dict]:
    """Items, and per model the store rows and trace texts of the symbolic templates' usable rows."""
    items = score.pool_items()
    sym = set(symbolic_templates(items))
    rows, texts = {}, {}
    for key in models:
        store = [json.loads(l) for l in (score.SCORES / 'main' / f'{key}.jsonl').read_text(encoding='utf-8').splitlines()]
        texts[key] = score.texts_matching('main', key, store)
        rows[key] = {r['item_id']: r for r in store if r['template_id'] in sym and usable(r)}
    return items, rows, texts


def draw(items: dict, rows: dict, texts: dict, controls: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    picked = []
    for key in sorted(rows):
        for i, r in sorted(rows[key].items()):
            if r['answer']['label'] in ('incorrect', 'partial'):
                picked.append((key, i, r['answer']['label']))
    for tid in symbolic_templates(items):
        by_model = collections.defaultdict(list)
        for key in sorted(rows):
            for i, r in sorted(rows[key].items()):
                if r['template_id'] == tid and r['answer']['label'] == 'correct':
                    by_model[key].append(i)
        keys = sorted(by_model)
        rng.shuffle(keys)
        pools = {k: rng.sample(by_model[k], len(by_model[k])) for k in keys}
        took = 0
        while took < controls and any(pools.values()):
            for k in keys:
                if pools[k] and took < controls:
                    picked.append((k, pools[k].pop(), 'control'))
                    took += 1
    out = []
    for key, i, group in picked:
        it, r = items[i], rows[key][i]
        ref = answer.segment(it['solution']).strip()
        out.append({'code': code_of(key, i),
                    'pool': {'branch': it['branch'], 'question': it['question'], 'reference': ref,
                             'reference_value': reference_value(it['template_id'], ref),
                             'stated': answer.segment(texts[key][i]).strip(), 'trace': texts[key][i],
                             'solution': it['solution']},
                    'key': {'item_id': i, 'model': key, 'template_id': it['template_id'], 'branch': it['branch'],
                            'level': it['level'], 'check_label': r['answer']['label'], 'group': group}})
    return out


def families(tasks: list[dict], seed: int) -> dict[str, str]:
    """Problem-family letters per branch: the branch's templates in a seeded order, A, B, ..."""
    by_branch = collections.defaultdict(set)
    for t in tasks:
        by_branch[t['key']['branch']].add(t['key']['template_id'])
    out = {}
    for b, tids in by_branch.items():
        order = sorted(tids)
        random.Random(f'{seed}|families|{b}').shuffle(order)
        out.update({tid: string.ascii_uppercase[k] for k, tid in enumerate(order)})
    return out


def assign(tasks: list[dict], roster: list[dict], seed: int, readers: int, overlap: float) -> dict:
    """Each item to one expert of its branch in turn within each family, a share of each family also to a second
    (or every item to two, readers=2); queues run family by family, shuffled per expert within a family."""
    by_branch = collections.defaultdict(list)
    for e in roster:
        by_branch[e['branch']].append(e['id'])
    fam = {t['code']: t['pool']['group'] for t in tasks}
    queues = {}
    per_family = collections.defaultdict(list)
    for t in sorted(tasks, key=lambda t: t['code']):
        per_family[(t['key']['branch'], t['pool']['group'])].append(t['code'])
    turn = collections.Counter()
    for (branch, f), codes in sorted(per_family.items()):
        ids = sorted(by_branch[branch])
        if not ids:
            raise SystemExit(f'no expert for branch {branch}')
        rng = random.Random(f'{seed}|assign|{branch}|{f}')
        rng.shuffle(codes)
        n_twice = len(codes) if readers >= 2 else (round(overlap * len(codes)) if len(codes) > 1 else 0)
        for j, c in enumerate(codes):
            first = ids[turn[branch] % len(ids)]
            turn[branch] += 1
            queues.setdefault(first, {'branch': branch, 'codes': []})['codes'].append(c)
            if j < n_twice and len(ids) > 1:
                second = ids[(ids.index(first) + 1) % len(ids)]
                queues.setdefault(second, {'branch': branch, 'codes': []})['codes'].append(c)
    for aid, q in queues.items():
        rng = random.Random(f'{seed}|order|{aid}')
        ordered = []
        for f in sorted({fam[c] for c in q['codes']}):
            block = sorted(c for c in q['codes'] if fam[c] == f)
            rng.shuffle(block)
            ordered += block
        q['codes'] = ordered
    return queues


def build(tasks_dir: Path, dist_dir: Path, seed: int, controls: int, readers: int, overlap: float,
          report: Path | None) -> dict:
    items, rows, texts = candidates(list(ROSTER))
    tasks = draw(items, rows, texts, controls, seed)
    codes = [t['code'] for t in tasks]
    if len(set(codes)) != len(codes):
        raise SystemExit('code collision')
    fam = families(tasks, seed)
    for t in tasks:
        t['pool']['group'] = fam[t['key']['template_id']]
        t['key']['family'] = fam[t['key']['template_id']]
        blob = json.dumps(t['pool'], ensure_ascii=False)
        if t['key']['template_id'] in blob:
            raise SystemExit(f"refusing: {t['code']}'s payload names its template")
    queues = assign(tasks, experts(), seed, readers, overlap)
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
    comp = composition(tasks, queues, seed, controls, readers, overlap)
    (tasks_dir / 'build.json').write_text(json.dumps(comp, indent=1), encoding='utf-8')
    if report:
        report.write_text('\n'.join(sent_lines(comp) + ['## What came back', '', 'Nothing yet.', '']),
                          encoding='utf-8', newline='\n')
    return comp


def store_commit() -> str:
    try:
        return json.loads((score.SCORES / 'main' / 'CONFIG.json').read_text(encoding='utf-8')).get('git', '?')[:7]
    except Exception:  # noqa: BLE001
        return '?'


def composition(tasks: list[dict], queues: dict, seed: int, controls: int, readers: int, overlap: float) -> dict:
    by_t = collections.defaultdict(collections.Counter)
    by_m = collections.defaultdict(collections.Counter)
    for t in tasks:
        by_t[t['key']['template_id']][t['key']['group']] += 1
        by_m[t['key']['model']][t['key']['group']] += 1
    reads = collections.Counter(c for q in queues.values() for c in q['codes'])
    return {'built': dt.date.today().isoformat(), 'seed': seed, 'store_commit': store_commit(),
            'controls_per_template': controls, 'readers': readers, 'overlap': overlap,
            'items': len(tasks), 'readings': sum(reads.values()), 'read_twice': sum(n >= 2 for n in reads.values()),
            'by_group': dict(collections.Counter(t['key']['group'] for t in tasks)),
            'by_template': {tid: dict(c) for tid, c in sorted(by_t.items())},
            'by_model': {k: dict(c) for k, c in sorted(by_m.items())},
            'experts': {aid: {'branch': q['branch'], 'items': len(q['codes'])} for aid, q in sorted(queues.items())}}


def sent_lines(comp: dict) -> list[str]:
    L = ['# Symbolic grading (E2)', '',
         'Generated by `build_kit.py` (what was sent) and `score_grades.py` (what came back). Counts only; the kits, '
         'the keyfile and the experts\' files stay local. The machine-readable grades are in `grades.json`.', '',
         '## What was sent', '',
         f"Built {comp['built']} from the main store at `{comp['store_commit']}`, seed {comp['seed']}. Every answered, "
         'usable verdict of the eleven models on the nine symbolic templates that the answer check scored incorrect '
         f"or partial, and {comp['controls_per_template']} it scored correct per template as controls. The expert sees "
         'the problem, the reference answer and the model\'s final answer as the check reads it, and grades it '
         'equivalent, not equivalent or unreadable; the model, the verdict and the template are hidden.'
         + (f" Each item goes to one expert of its branch, and {comp['overlap']:.0%} of each template's items to a "
            'second.' if comp['readers'] == 1 else ' Each item goes to two experts of its branch.'), '',
         '| template | incorrect | partial | controls (correct) |', '|---|---:|---:|---:|']
    for tid, c in comp['by_template'].items():
        L.append(f"| `{tid.replace('template_', '')}` | {c.get('incorrect', 0)} | {c.get('partial', 0)} | "
                 f"{c.get('control', 0)} |")
    g = comp['by_group']
    L += [f"| all | {g.get('incorrect', 0)} | {g.get('partial', 0)} | {g.get('control', 0)} |", '',
          '| model | incorrect | partial | controls |', '|---|---:|---:|---:|']
    for k, c in comp['by_model'].items():
        L.append(f"| `{k}` | {c.get('incorrect', 0)} | {c.get('partial', 0)} | {c.get('control', 0)} |")
    L += ['', f"Items {comp['items']}, readings {comp['readings']}, items read twice {comp['read_twice']}.", '',
          '| expert | branch | items |', '|---|---|---:|']
    for aid, e in comp['experts'].items():
        L.append(f"| `{aid}` | {e['branch'].replace('_engineering', '')} | {e['items']} |")
    L.append('')
    return L


# ------------------------------------------------------------------ selftest

def selftest() -> int:
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-symbolic-kit-'))
    try:
        comp = build(tmp / 'tasks', tmp / 'dist', SEED, 4, 1, 0.2, None)
        keyfile = json.loads((tmp / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
        queues = json.loads((tmp / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
        items, rows, _texts = candidates(list(ROSTER))
        want = collections.Counter(r['answer']['label'] for k in rows for r in rows[k].values())
        g = comp['by_group']
        assert g['incorrect'] == want['incorrect'] and g['partial'] == want['partial'], (g, want)
        assert all(c.get('control', 0) <= 4 for c in comp['by_template'].values())
        reads = collections.Counter(c for q in queues.values() for c in q['codes'])
        assert set(reads) == set(keyfile) and all(1 <= n <= 2 for n in reads.values())
        for aid, q in queues.items():
            assert all(keyfile[c]['branch'] == q['branch'] for c in q['codes']), aid
            fams = [keyfile[c]['family'] for c in q['codes']]
            assert fams == sorted(fams), f'{aid}: the queue is not grouped by family'
        for f in (tmp / 'dist').rglob('pool.json'):
            blob = f.read_text(encoding='utf-8')
            assert 'template_' not in blob and '"model"' not in blob and 'check_label' not in blob, f
        ber = [t for c, t in json.loads((tmp / 'tasks' / 'pool.json').read_text(encoding='utf-8')).items()
               if keyfile[c]['template_id'] == BER]
        assert ber and all(t['reference_value'] for t in ber), 'a BER item lacks its reference value'
        print(f"built {comp['items']} items ({g}), {comp['readings']} readings, {comp['read_twice']} read twice; "
              'every incorrect and partial verdict drawn; queues grouped by family; payloads clean')
        from streamlit.testing.v1 import AppTest
        aid = sorted(a for a in queues if a.startswith('ele-'))[0]
        folder = tmp / 'dist' / 'electrical_engineering'
        at = AppTest.from_file(str(folder / 'app.py'), default_timeout=120)
        at.run()
        assert not at.exception, at.exception
        at.sidebar.selectbox[0].select(aid).run()
        assert not at.exception, at.exception
        first = queues[aid]['codes'][0]
        assert any(first in w.value for w in at.subheader), [w.value for w in at.subheader]
        at.radio[0].set_value('not equivalent').run()
        at.button[0].click().run()
        assert not at.exception, at.exception
        row = json.loads((folder / f'{aid}.jsonl').read_text(encoding='utf-8').splitlines()[0])
        assert row['code'] == first and row['grade'] == 'not equivalent', row
        print('the app renders, takes a grade and moves on')
        print('SELFTEST OK')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--controls', type=int, default=4, help='correct verdicts per template, as controls')
    ap.add_argument('--readers', type=int, default=1, choices=(1, 2))
    ap.add_argument('--overlap', type=float, default=0.2, help='with --readers 1: the share also read by a second')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    comp = build(TASKS, DIST, a.seed, a.controls, a.readers, a.overlap, GRADES_MD)
    print('\n'.join(sent_lines(comp)))
    print(f'kits: {DIST}  (one folder per branch; send it to the branch\'s experts)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
