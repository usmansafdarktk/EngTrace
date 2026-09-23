"""Prove the Layer 2 pipeline end to end with synthetic labels before any expert starts.

    python -m template_annotation_23092026.layer2.simulate      # writes labels_simulated/, then scores them

Synthetic experts approve every real template with scores drawn near 5, reject
each plant with probability 0.8, and mismatch the hand check on 5% of items. The
numbers score.py then prints are NOT results and the report says so on its first
line; this only proves that build -> app rows -> score runs, and that plants are
scored as plants.
"""
from __future__ import annotations

import datetime as dt
import json
import random
import shutil
from pathlib import Path

from . import score

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
SIM = HERE / 'labels_simulated'


def main() -> None:
    key = {}
    for ln in (TASKS / 'keyfile.jsonl').open(encoding='utf8'):
        if ln.strip():
            k = json.loads(ln)
            key[k['code']] = k
    assignment = json.loads((TASKS / 'assignment.json').read_text(encoding='utf8'))
    rng = random.Random(1)
    if SIM.exists():
        shutil.rmtree(SIM)
    SIM.mkdir()
    t = dt.datetime(2026, 10, 1, 9, 0, tzinfo=dt.timezone.utc)
    n = 0
    for aid, mine in assignment.items():
        with (SIM / f'{aid}.jsonl').open('w', encoding='utf8') as fh:
            for pos, code in enumerate(mine['codes'], 1):
                plant = key[code]['kind'] == 'plant'
                reject = plant and rng.random() < 0.8
                opened = t + dt.timedelta(minutes=7 * pos)
                row = {'annotator_id': aid, 'branch': mine['branch'], 'code': code, 'position': pos,
                       'opened_at': opened.isoformat(timespec='seconds'), 'hand_answer': 'x',
                       'hand_numbers': [1.0], 'hand_match': (rng.random() > 0.05) and not plant,
                       'hand_submitted_at': (opened + dt.timedelta(minutes=3)).isoformat(timespec='seconds'),
                       'hand_instance_seed': 2101, 'instances_viewed': 1,
                       'scores': {d: (rng.choice([2, 3]) if reject else rng.choice([4, 5, 5, 5]))
                                  for d in score.DIMS},
                       'decision': 'Reject' if reject else 'Approve',
                       'defects': ['constant or table value'] if reject else [],
                       'feedback': 'simulated rejection' if reject else '', 'confidence': 4,
                       'submitted_at': (opened + dt.timedelta(minutes=6)).isoformat(timespec='seconds'),
                       'source': 'app'}
                fh.write(json.dumps(row) + '\n')
                n += 1
    print(f'{n} simulated rows for {len(assignment)} experts under {SIM}')
    score.LABELS = SIM
    score.main()
    res = HERE / 'RESULTS.md'
    res.write_text('**SIMULATED LABELS - NOT RESULTS.** ' + res.read_text(encoding='utf8'), encoding='utf8')


if __name__ == '__main__':
    main()
