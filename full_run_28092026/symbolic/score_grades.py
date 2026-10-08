"""E2 (WS-E): score the experts' grades of symbolic final answers.

    python -m full_run_28092026.symbolic.score_grades                # FREE: experts_filled_symbolic/ -> grades.json, GRADES.md
    python -m full_run_28092026.symbolic.score_grades --returns DIR  # FREE: the returns from another folder
    python -m full_run_28092026.symbolic.score_grades --selftest     # FREE: simulated returns, temp folder

WHAT IT READS. The keyfile and the queues build_kit.py wrote to tasks/, and every <id>.jsonl under the returns folder
(any layout). A row whose code is not in the keyfile, or was not assigned to its expert, or whose grade is not one of
the three, stops the run. The last submission per expert and item counts.

WHAT IT COMPUTES. Per item, the grade: its reader's, or the two readers' when they agree, else `split`. Per template
and per model, the grades by the check's verdict group (incorrect, partial, control: a correct verdict drawn as a
control). The readers' agreement on the items read twice (share and Cohen's kappa). What an expert-adjudicated score
would credit: per model, the incorrect and partial verdicts graded equivalent. The by-numbers rule's precision on the
controls: the share of correct verdicts the experts grade equivalent.

WHAT IT WRITES. grades.json, keyed by "<item_id>|<model>": the template, branch, level, the check's verdict, the group,
the grade and the grade counts; with the summary. No expert ids and no notes: those stay in the local returns, where
validate.py reads them. GRADES.md: counts only. WS-C reads grades.json for the expert-adjudicated sensitivity.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import random
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026.symbolic import build_kit  # noqa: E402

RETURNS = HERE / 'experts_filled_symbolic'
GRADES_JSON = HERE / 'grades.json'
GRADES = ['equivalent', 'not equivalent', 'unreadable']
GROUPS = build_kit.GROUPS


def cohen(pairs: list[tuple[str, str]], cats: list[str]) -> float | None:
    n = len(pairs)
    if n < 2:
        return None
    po = sum(a == b for a, b in pairs) / n
    pe = sum((sum(a == c for a, _ in pairs) / n) * (sum(b == c for _, b in pairs) / n) for c in cats)
    return None if pe >= 1 - 1e-12 else (po - pe) / (1 - pe)


def read_returns(returns: Path, tasks: Path) -> tuple[dict, dict, dict]:
    keyfile = json.loads((tasks / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((tasks / 'assignment.json').read_text(encoding='utf-8'))
    owners = collections.defaultdict(set)
    for aid, q in queues.items():
        for c in q['codes']:
            owners[c].add(aid)
    rows = {}
    for f in sorted(returns.rglob('*.jsonl')):
        for line in f.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r['code'] not in keyfile:
                raise SystemExit(f"{f.name}: code {r['code']} is not in this build's keyfile")
            if r['annotator_id'] not in owners[r['code']]:
                raise SystemExit(f"{f.name}: {r['code']} was not assigned to {r['annotator_id']}")
            if r['grade'] not in GRADES:
                raise SystemExit(f"{f.name}: {r['code']} has grade {r['grade']!r}")
            rows[(r['annotator_id'], r['code'])] = r
    return keyfile, queues, rows


def consensus(grades: list[str]) -> str | None:
    if not grades:
        return None
    return grades[0] if len(set(grades)) == 1 else 'split'


def summarise(keyfile: dict, queues: dict, rows: dict, build: dict) -> dict:
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        by_code[c][aid] = r['grade']
    items = {}
    for c, k in sorted(keyfile.items(), key=lambda kv: (kv[1]['template_id'], kv[1]['item_id'], kv[1]['model'])):
        gs = [by_code[c][a] for a in sorted(by_code.get(c, {}))]
        items[f"{k['item_id']}|{k['model']}"] = {
            'item_id': k['item_id'], 'model': k['model'], 'template_id': k['template_id'], 'branch': k['branch'],
            'level': k['level'], 'check_label': k['check_label'], 'group': k['group'],
            'grade': consensus(gs), 'readings': len(gs), 'counts': {g: gs.count(g) for g in GRADES}}
    cols = GRADES + ['split', 'ungraded']

    def table(field: str) -> dict:
        out = collections.defaultdict(lambda: {g: collections.Counter() for g in GROUPS})
        for it in items.values():
            out[it[field]][it['group']][it['grade'] or 'ungraded'] += 1
        return {k: {g: {c: v[g][c] for c in cols} for g in GROUPS} for k, v in sorted(out.items())}

    pairs = [tuple(by_code[c][a] for a in sorted(by_code[c])) for c in by_code if len(by_code[c]) == 2]
    ck = cohen(pairs, GRADES)
    credit = collections.defaultdict(collections.Counter)
    for it in items.values():
        if it['group'] in ('incorrect', 'partial'):
            credit[it['model']][f"{it['group']}_graded"] += it['grade'] in GRADES
            credit[it['model']][f"{it['group']}_equivalent"] += it['grade'] == 'equivalent'
    ctrl = [it for it in items.values() if it['group'] == 'control' and it['grade'] in GRADES]
    returned = sorted({a for a, _ in rows})
    return {'generated_by': 'full_run_28092026/symbolic/score_grades.py',
            'scored_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
            'build': {k: build.get(k) for k in ('built', 'seed', 'store_commit', 'readers', 'overlap',
                                                'controls_per_template')},
            'returns': {'readings': len(rows), 'assigned': sum(len(q['codes']) for q in queues.values()),
                        'experts_returned': len(returned), 'experts': len(queues),
                        'items': len(items), 'items_graded': sum(it['grade'] is not None for it in items.values()),
                        'items_split': sum(it['grade'] == 'split' for it in items.values())},
            'agreement': {'items_read_twice': len(pairs),
                          'share_same': round(sum(a == b for a, b in pairs) / len(pairs), 3) if pairs else None,
                          'cohen_kappa': None if ck is None else round(ck, 3)},
            'by_template': table('template_id'), 'by_model': table('model'),
            'adjudicated_credit_by_model': {m: dict(v) for m, v in sorted(credit.items())},
            'controls': {'graded': len(ctrl), 'equivalent': sum(it['grade'] == 'equivalent' for it in ctrl)},
            'items': items}


def result_lines(res: dict) -> list[str]:
    r, a = res['returns'], res['agreement']
    same = '-' if a['share_same'] is None else f"{a['share_same']:.3f}"
    kap = '-' if a['cohen_kappa'] is None else f"{a['cohen_kappa']:.3f}"
    L = ['## What came back', '',
         f"Readings returned {r['readings']} of {r['assigned']} assigned, from {r['experts_returned']} of "
         f"{r['experts']} experts; items graded {r['items_graded']} of {r['items']}, of which {r['items_split']} split "
         'between their two readers.', '',
         f"Items read twice: {a['items_read_twice']}; the two grades agree on {same} of them; Cohen's kappa {kap}.", '',
         '### By template', '',
         '| template | verdict | equivalent | not equivalent | unreadable | split | ungraded |',
         '|---|---|---:|---:|---:|---:|---:|']
    for tid, v in res['by_template'].items():
        for g in GROUPS:
            c = v[g]
            if sum(c.values()):
                L.append(f"| `{tid.replace('template_', '')}` | {g} | {c['equivalent']} | {c['not equivalent']} | "
                         f"{c['unreadable']} | {c['split']} | {c['ungraded']} |")
    L += ['', '### By model: what an expert-adjudicated score would credit', '',
          '| model | incorrect verdicts graded | of them equivalent | partial verdicts graded | of them equivalent |',
          '|---|---:|---:|---:|---:|']
    for m, v in res['adjudicated_credit_by_model'].items():
        L.append(f"| `{m}` | {v.get('incorrect_graded', 0)} | {v.get('incorrect_equivalent', 0)} | "
                 f"{v.get('partial_graded', 0)} | {v.get('partial_equivalent', 0)} |")
    c = res['controls']
    L += ['', f"Controls (verdicts the check scored correct): {c['equivalent']} of {c['graded']} graded equivalent.", '']
    return L


def score(returns: Path, tasks: Path, out_json: Path | None, out_md: Path | None) -> dict:
    keyfile, queues, rows = read_returns(returns, tasks)
    build = json.loads((tasks / 'build.json').read_text(encoding='utf-8'))
    res = summarise(keyfile, queues, rows, build)
    if out_json:
        out_json.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8', newline='\n')
    if out_md:
        out_md.write_text('\n'.join(build_kit.sent_lines(build) + result_lines(res)) + '\n', encoding='utf-8',
                          newline='\n')
    print('\n'.join(result_lines(res)))
    return res


# ------------------------------------------------------------------ selftest

def simulate(tasks: Path, out: Path, seed: int = 1, drop: str | None = None, flip: float = 0.1) -> int:
    """Synthetic returns: a reader grades equivalent with probability 0.3 / 0.5 / 0.95 for an incorrect / partial /
    control verdict; the second reader of an item copies the first's grade except with probability `flip`."""
    keyfile = json.loads((tasks / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((tasks / 'assignment.json').read_text(encoding='utf-8'))
    rng = random.Random(seed)
    first = {}
    t0 = dt.datetime(2026, 10, 8, 9, 0, tzinfo=dt.timezone.utc)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for aid in sorted(queues):
        if aid == drop:
            continue
        q = queues[aid]
        with (out / f'{aid}.jsonl').open('w', encoding='utf-8') as fh:
            for pos, c in enumerate(q['codes'], 1):
                p = {'incorrect': 0.3, 'partial': 0.5, 'control': 0.95}[keyfile[c]['group']]
                g = 'equivalent' if rng.random() < p else rng.choice(['not equivalent', 'not equivalent', 'unreadable'])
                if c in first and rng.random() > flip:
                    g = first[c]
                first.setdefault(c, g)
                opened = t0 + dt.timedelta(minutes=pos)
                fh.write(json.dumps({'annotator_id': aid, 'branch': q['branch'], 'kind': 'symbolic', 'code': c,
                                     'position': pos, 'opened_at': opened.isoformat(), 'grade': g, 'note': '',
                                     'submitted_at': (opened + dt.timedelta(seconds=40)).isoformat(),
                                     'source': 'simulated'}) + '\n')
                n += 1
    return n


def selftest() -> int:
    assert consensus(['equivalent']) == 'equivalent' and consensus(['equivalent', 'equivalent']) == 'equivalent'
    assert consensus(['equivalent', 'unreadable']) == 'split' and consensus([]) is None
    assert abs(cohen([('a', 'a'), ('b', 'b'), ('a', 'b'), ('b', 'b')], ['a', 'b']) - 0.5) < 1e-12
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-symbolic-score-'))
    try:
        tasks = tmp / 'tasks'
        build_kit.build(tasks, tmp / 'dist', build_kit.SEED, 4, 1, 0.2, None)
        n = simulate(tasks, tmp / 'sim')
        res = score(tmp / 'sim', tasks, tmp / 'grades.json', tmp / 'GRADES.md')
        r = res['returns']
        assert r['readings'] == n == r['assigned'] and r['items_graded'] == r['items'], r
        assert res['agreement']['items_read_twice'] > 0 and 0 < res['agreement']['cohen_kappa'] <= 1, res['agreement']
        assert sum(sum(c.values()) for v in res['by_template'].values() for c in v.values()) == r['items']
        blob = (tmp / 'grades.json').read_text(encoding='utf-8')
        assert '"ele-' not in blob and '"mec-' not in blob and '"note"' not in blob, 'an expert id or a note leaked'
        k0 = next(iter(res['items']))
        assert k0.count('|') == 1 and res['items'][k0]['grade'] in GRADES + ['split']
        res2 = score_partial(tasks, tmp)
        assert res2['returns']['items_graded'] < res2['returns']['items']
        bad = tmp / 'bad'
        bad.mkdir()
        (bad / 'x.jsonl').write_text(json.dumps({'annotator_id': 'ele-1', 'code': 'SG-nope', 'grade': 'equivalent'})
                                     + '\n', encoding='utf-8')
        try:
            score(bad, tasks, None, None)
            raise AssertionError('an unknown code was accepted')
        except SystemExit:
            pass
        print('SELFTEST OK: every reading scored; agreement on the overlap computed; counts add up; no expert id or '
              'note in grades.json; a missing expert leaves items ungraded; an unknown code stops the run')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def score_partial(tasks: Path, tmp: Path) -> dict:
    simulate(tasks, tmp / 'partial', seed=2, drop='ele-2')
    return score(tmp / 'partial', tasks, None, None)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--returns', type=Path, default=RETURNS)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.returns.exists():
        raise SystemExit(f'no returns folder: {a.returns}')
    score(a.returns, build_kit.TASKS, GRADES_JSON, build_kit.GRADES_MD)
    return 0


if __name__ == '__main__':
    sys.exit(main())
