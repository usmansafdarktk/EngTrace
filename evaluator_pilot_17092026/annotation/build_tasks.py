"""Build the annotation tasks: what each expert sees, in what order.

    python evaluator_pilot_17092026/annotation/build_tasks.py [--annotators annotators.json]

Stage 3 labels the frozen 300 traces so the six evaluator candidates can finally be
scored against something other than each other.

DESIGN (decided 2026-09-20)
  Every trace is labelled by THREE experts, and those three are specialists in the
  trace's own branch: 5 branches x 12 items x 5 models = 60 traces per branch, 3
  experts per branch, 15 experts, 60 traces each (~10 hours). Total 900
  trace-annotations.

  BLIND. An expert never sees which model wrote a trace, or any evaluator's opinion
  of it. Each trace appears under a code (T-xxxxxx); the code -> (item, model) map
  is written to keyfile.jsonl, which is NOT given to annotators.

  ORDER is shuffled per expert (seeded by their id), so fatigue and drift do not
  line up across the three people labelling the same trace.

  CALIBRATION. A shared set of 10 traces (2 per branch, levels mixed) is assigned to
  ALL 15 experts and labelled first. Its purpose is protocol agreement across
  branches; labels an expert gives outside their own branch are diagnostic only and
  never enter the ground truth (the field `for_truth` says which is which).

STEPS are split with E0's own `extract_steps`, exactly as E2 and the judge probe
split them, so a label on step k refers to the same step in every evaluator.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import random
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_REPO = os.path.dirname(_PILOT)
sys.path[:0] = [os.path.join(_PILOT, 'evaluators'), os.path.join(_REPO, 'evaluation')]

import milestones as ms                       # noqa: E402
from engineering_parser import extract_steps  # noqa: E402

TASKS = os.path.join(_HERE, 'tasks')
GOLD_MODELS = ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b')
CALIBRATION_PER_BRANCH = 2
CODE_SALT = b'engtrace-stage3'                # fixed, so codes are stable across rebuilds


def code_for(item_id: str, model_key: str) -> str:
    h = hashlib.blake2b(('%s|%s' % (item_id, model_key)).encode(), key=CODE_SALT, digest_size=3)
    return 'T-' + h.hexdigest()


def load():
    items = {}
    for l in open(os.path.join(_PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8'):
        r = json.loads(l)
        items[r['item_id']] = r
    mst = ms.build_all(os.path.join(_PILOT, 'slice', 'manifest.jsonl'))
    traces = {}
    for f in glob.glob(os.path.join(_PILOT, 'traces', '*.jsonl')):
        mk = os.path.basename(f)[:-6]
        if mk not in GOLD_MODELS:
            continue
        for l in open(f, encoding='utf-8'):
            r = json.loads(l)
            if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
                traces[(mk, r['item_id'])] = r
    return items, mst, traces


def build_pool(items, mst, traces):
    pool = {}
    for (mk, iid), tr in sorted(traces.items()):
        steps = [s for s in (x.strip() for x in extract_steps(tr['text'])[0]) if s]
        it = items[iid]
        pool[code_for(iid, mk)] = {
            'code': code_for(iid, mk),
            'branch': it['branch'], 'level': it['level'], 'answer_type': it['answer_type'],
            'question': it['question'],
            'reference_solution': it['solution'],
            'steps': steps,
            'milestones': [{'id': m['id'], 'value': m['value']}
                           for m in mst[iid]['milestones']],
            'item_sha256': it['sha256'],
            'trace_sha256': hashlib.sha256(tr['text'].encode('utf-8')).hexdigest(),
        }
    return pool


def calibration_set(pool):
    """2 traces per branch, levels spread, chosen deterministically."""
    out = []
    for br in sorted({t['branch'] for t in pool.values()}):
        rows = sorted((t for t in pool.values() if t['branch'] == br), key=lambda t: t['code'])
        by_level = {}
        for t in rows:
            by_level.setdefault(t['level'], []).append(t)
        picks = [by_level[lv][0] for lv in sorted(by_level)][:CALIBRATION_PER_BRANCH]
        out += [t['code'] for t in picks]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotators', default=os.path.join(_HERE, 'annotators.json'))
    a = ap.parse_args()
    people = json.load(open(a.annotators, encoding='utf-8'))['annotators']

    items, mst, traces = load()
    pool = build_pool(items, mst, traces)
    if len(pool) != 300:
        raise SystemExit('expected 300 traces, built %d' % len(pool))
    calib = calibration_set(pool)

    os.makedirs(TASKS, exist_ok=True)
    with open(os.path.join(TASKS, 'pool.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(pool, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(TASKS, 'keyfile.jsonl'), 'w', encoding='utf-8', newline='\n') as fh:
        for (mk, iid) in sorted(traces):
            fh.write(json.dumps({'code': code_for(iid, mk), 'item_id': iid, 'model_key': mk}) + '\n')

    assignment = {}
    for p in people:
        own = sorted(t['code'] for t in pool.values() if t['branch'] == p['branch'])
        mine = [c for c in calib if c not in own]          # calibration outside own branch
        rnd = random.Random(int(hashlib.blake2b(p['id'].encode(), digest_size=4).hexdigest(), 16))
        rest = [c for c in own if c not in calib]
        rnd.shuffle(rest)
        order = [c for c in calib if c in own] + mine + rest
        assignment[p['id']] = {
            'id': p['id'], 'name': p.get('name', ''), 'branch': p['branch'],
            'order': order,
            'calibration': [c for c in order if c in calib],
            'for_truth': [c for c in order if pool[c]['branch'] == p['branch']],
        }
    with open(os.path.join(TASKS, 'assignment.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(assignment, fh, indent=1)

    n_steps = sum(len(t['steps']) for t in pool.values())
    n_ms = sum(len(t['milestones']) for t in pool.values())
    print('pool: %d traces, %d steps, %d milestones' % (len(pool), n_steps, n_ms))
    print('calibration set (all experts): %d traces' % len(calib))
    print('%-14s %-24s %6s %6s %6s' % ('annotator', 'branch', 'tasks', 'truth', 'calib'))
    for k, v in sorted(assignment.items()):
        print('%-14s %-24s %6d %6d %6d' % (k, v['branch'], len(v['order']),
                                           len(v['for_truth']), len(v['calibration'])))
    tot = sum(len(v['order']) for v in assignment.values())
    print('total trace-annotations: %d (%d for ground truth, the rest calibration)'
          % (tot, sum(len(v['for_truth']) for v in assignment.values())))
    print('written to %s' % TASKS)


if __name__ == '__main__':
    main()
