"""The file route: a workbook per expert, and the import that turns it into label rows.

    python -m template_annotation_23092026.layer2.workbooks export            # tasks/workbooks/<id>.json
    python -m template_annotation_23092026.layer2.workbooks check  <file>     # validate a returned one
    python -m template_annotation_23092026.layer2.workbooks import <file>     # -> labels/<id>.jsonl (source: workbook)

A workbook holds the expert's items in queue order with empty fields to fill. On
import every filled item becomes a label row identical in shape to the app's, so
score.py reads both the same way. No timestamps exist on this route; that is the
price of the file, and score.py reports dwell only for app rows.
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
LABELS = HERE / 'labels'
DEFECTS = ['physics or scenario implausible', 'governing equation or formula', 'constant or table value',
           'unit or conversion', 'sign or direction', 'arithmetic: a step does not follow',
           'final answer wrong or wrong unit', 'question ambiguous or unsolvable', 'wording or formatting']
FIELDS = {'hand_answer': '', 'physical_plausibility': None, 'mathematical_correctness': None,
          'pedagogical_clarity': None, 'decision': '', 'defects': [], 'feedback': '', 'confidence': None}


def export_one(aid: str, pool: dict, mine: dict) -> dict:
    items = []
    for c in mine['codes']:
        it = pool[c]
        items.append({'code': c, 'area': it['area'], 'question': it['instances'][0]['question'],
                      'solution': it['instances'][0]['solution'], 'source': it['source'], **FIELDS})
    return {'_instructions': ['Fill hand_answer BEFORE reading the solution. Scores are integers 1-5. '
                              'decision is "Approve" or "Reject"; a Reject needs feedback and at least one '
                              'entry in defects, chosen from _defect_types. confidence is 1-5.'],
            '_defect_types': DEFECTS, 'annotator_id': aid, 'branch': mine['branch'], 'items': items}


def export_all() -> None:
    pool = json.loads((TASKS / 'pool.json').read_text(encoding='utf8'))
    assignment = json.loads((TASKS / 'assignment.json').read_text(encoding='utf8'))
    out = TASKS / 'workbooks'
    out.mkdir(exist_ok=True)
    for aid, mine in assignment.items():
        (out / f'{aid}.json').write_text(json.dumps(export_one(aid, pool, mine), ensure_ascii=False, indent=1),
                                         encoding='utf8')
    print(f'{len(assignment)} workbooks under {out}')


def problems(item: dict) -> list[str]:
    bad = []
    filled = any(item.get(k) not in (None, '', []) for k in FIELDS)
    if not filled:
        return ['empty']
    for k in ('physical_plausibility', 'mathematical_correctness', 'pedagogical_clarity', 'confidence'):
        v = item.get(k)
        if not isinstance(v, int) or not 1 <= v <= 5:
            bad.append(f'{k} must be an integer 1-5')
    if item.get('decision') not in ('Approve', 'Reject'):
        bad.append('decision must be Approve or Reject')
    if item.get('decision') == 'Reject':
        if not (item.get('feedback') or '').strip():
            bad.append('a Reject needs feedback')
        if not item.get('defects'):
            bad.append('a Reject needs at least one defect type')
    for d in item.get('defects') or []:
        if d not in DEFECTS:
            bad.append(f'unknown defect type {d!r}')
    return bad


def check(path: Path) -> tuple[list[dict], dict[str, list[str]]]:
    wb = json.loads(path.read_text(encoding='utf8'))
    ok, bad = [], {}
    for it in wb['items']:
        p = problems(it)
        if not p:
            ok.append(it)
        elif p != ['empty']:
            bad[it['code']] = p
    return ok, bad


def import_(path: Path) -> None:
    wb = json.loads(path.read_text(encoding='utf8'))
    ok, bad = check(path)
    for c, p in bad.items():
        print(f'  skipped {c}: {"; ".join(p)}')
    LABELS.mkdir(exist_ok=True)
    out = LABELS / f"{wb['annotator_id']}.jsonl"
    with out.open('a', encoding='utf8') as fh:
        for it in ok:
            row = {'annotator_id': wb['annotator_id'], 'branch': wb['branch'], 'code': it['code'],
                   'position': [x['code'] for x in wb['items']].index(it['code']) + 1,
                   'opened_at': None, 'hand_answer': it['hand_answer'], 'hand_numbers': None, 'hand_match': None,
                   'hand_submitted_at': None, 'hand_instance_seed': None, 'instances_viewed': 1,
                   'scores': {k: it[k] for k in ('physical_plausibility', 'mathematical_correctness', 'pedagogical_clarity')},
                   'decision': it['decision'], 'defects': it.get('defects') or [], 'feedback': it.get('feedback') or '',
                   'confidence': it['confidence'], 'submitted_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
                   'source': 'workbook'}
            fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    print(f'{len(ok)} items imported into {out}, {len(bad)} skipped')


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'export':
        export_all()
    elif cmd == 'check' and len(sys.argv) > 2:
        ok, bad = check(Path(sys.argv[2]))
        print(f'{len(ok)} valid items, {len(bad)} with problems')
        for c, p in bad.items():
            print(f'  {c}: {"; ".join(p)}')
    elif cmd == 'import' and len(sys.argv) > 2:
        import_(Path(sys.argv[2]))
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
