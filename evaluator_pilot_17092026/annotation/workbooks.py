"""The file route: a workbook per annotator, filled in any editor, validated on return.

    python evaluator_pilot_17092026/annotation/workbooks.py export        # write blank workbooks
    python evaluator_pilot_17092026/annotation/workbooks.py export --id civ-2
    python evaluator_pilot_17092026/annotation/workbooks.py check  FILE   # validate, change nothing
    python evaluator_pilot_17092026/annotation/workbooks.py import FILE   # validate, then load

Annotators who prefer not to use the app fill in `tasks/workbooks/<id>.json`: every
solution assigned to them, with the fields empty. `import` validates it and appends the
same rows the app writes to `labels/<id>.jsonl`, so both routes are interchangeable and
the scoring never learns which was used.

EXPORT CARRIES WORK ACROSS. A workbook is written from whatever that annotator has
already submitted (app or file), so someone can start in the app and continue in the
file, or the reverse, without retyping.

A workbook is checked for every permitted value, for steps and milestones matching the
solution, and for the trace hash - a workbook built from a different version of the
traces is refused rather than silently scored.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(_HERE, 'tasks')
BOOKS = os.path.join(TASKS, 'workbooks')
LABELS = os.path.join(_HERE, 'labels')

STEP_LABELS = ('correct', 'alternative_correct', 'incorrect', 'not_a_claim')
ERROR_TYPES = ('calculation', 'conceptual', 'unsupported')
MS_STATUS = ('reached', 'not_reached', 'not_needed')
FINAL = ('correct', 'incorrect', 'partial', 'not_stated')
SOUND = ('yes', 'no')
CONFIDENCE = ('high', 'medium', 'low')

HOWTO = [
    'Fill in the empty fields. Leave code, text, id and value untouched.',
    'label: one of %s' % ', '.join(STEP_LABELS),
    'error_type: %s - required when label is incorrect, empty otherwise' % ', '.join(ERROR_TYPES),
    'milestone status: one of %s' % ', '.join(MS_STATUS),
    'final_answer: %s' % ', '.join(FINAL),
    'reasoning_sound: yes or no;  confidence: %s' % ', '.join(CONFIDENCE),
    'note and comment are free text and may be left empty; flagged is true or false',
    'Partly filled workbooks are fine: solutions with nothing filled in are skipped on import.',
    'See the annotation guide for what the labels mean.',
]


def load_tasks():
    pool = json.load(open(os.path.join(TASKS, 'pool.json'), encoding='utf-8'))
    assignment = json.load(open(os.path.join(TASKS, 'assignment.json'), encoding='utf-8'))
    return pool, assignment


def submitted(who):
    """Last row per code from labels/<who>.jsonl, app or file."""
    out = {}
    p = os.path.join(LABELS, '%s.jsonl' % who)
    if os.path.exists(p):
        for line in open(p, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                out[r['code']] = r
    return out


def export(ids=None):
    pool, assignment = load_tasks()
    os.makedirs(BOOKS, exist_ok=True)
    for who, me in sorted(assignment.items()):
        if ids and who not in ids:
            continue
        prior = submitted(who)
        traces = []
        for code in me['order']:
            t = pool[code]
            p = prior.get(code, {})
            psteps = {s['index']: s for s in p.get('steps', [])}
            pms = {m['id']: m for m in p.get('milestones', [])}
            traces.append({
                'code': code,
                'discipline': t['branch'].replace('_', ' '),
                'level': t['level'],
                'calibration': code in me['calibration'],
                'problem': t['question'],
                'reference_solution': t['reference_solution'],
                'steps': [{'index': i, 'text': s,
                           'label': psteps.get(i, {}).get('label', ''),
                           'error_type': psteps.get(i, {}).get('error_type') or '',
                           'note': psteps.get(i, {}).get('note', '')}
                          for i, s in enumerate(t['steps'])],
                'milestones': [{'id': m['id'], 'value': m['value'],
                                'status': pms.get(m['id'], {}).get('status', '')}
                               for m in t['milestones']],
                'final_answer': p.get('final_answer', ''),
                'reasoning_sound': p.get('reasoning_sound', ''),
                'confidence': p.get('confidence', ''),
                'comment': p.get('comment', ''),
                'flagged': p.get('flagged', False),
                'trace_sha256': t['trace_sha256'], 'item_sha256': t['item_sha256'],
            })
        book = {'annotator': who, 'name': me['name'], 'discipline': me['branch'].replace('_', ' '),
                'how_to_fill_this_in': HOWTO, 'traces': traces}
        out = os.path.join(BOOKS, '%s.json' % who)
        with open(out, 'w', encoding='utf-8', newline='\n') as fh:
            json.dump(book, fh, ensure_ascii=False, indent=1)
        filled = sum(1 for t in traces if any(s['label'] for s in t['steps']))
        print('%-8s %3d solutions (%d already filled) -> %s' % (who, len(traces), filled, out))


def check_one(t, task):
    """Problems with one trace entry; [] if it is complete and valid."""
    bad = []
    if t.get('trace_sha256') != task['trace_sha256']:
        return ['trace hash does not match the current traces - workbook is out of date']
    if len(t.get('steps', [])) != len(task['steps']):
        return ['has %d steps, the solution has %d' % (len(t.get('steps', [])), len(task['steps']))]
    for s in t['steps']:
        n = 'step %s' % (s.get('index', '?') + 1 if isinstance(s.get('index'), int) else '?')
        if s.get('label') not in STEP_LABELS:
            bad.append('%s: label %r is not one of %s' % (n, s.get('label'), ', '.join(STEP_LABELS)))
        elif s['label'] == 'incorrect' and s.get('error_type') not in ERROR_TYPES:
            bad.append('%s: incorrect needs error_type, one of %s' % (n, ', '.join(ERROR_TYPES)))
        elif s['label'] != 'incorrect' and s.get('error_type'):
            bad.append('%s: error_type is only for incorrect steps' % n)
    ids = [m['id'] for m in task['milestones']]
    if [m.get('id') for m in t.get('milestones', [])] != ids:
        bad.append('milestones do not match the solution (expected %s)' % ', '.join(ids))
    else:
        for m in t['milestones']:
            if m.get('status') not in MS_STATUS:
                bad.append('milestone %s: status %r is not one of %s'
                           % (m['id'], m.get('status'), ', '.join(MS_STATUS)))
    for field, allowed in (('final_answer', FINAL), ('reasoning_sound', SOUND),
                           ('confidence', CONFIDENCE)):
        if t.get(field) not in allowed:
            bad.append('%s: %r is not one of %s' % (field, t.get(field), ', '.join(allowed)))
    if not isinstance(t.get('flagged', False), bool):
        bad.append('flagged must be true or false')
    return bad


def read_book(path):
    book = json.load(open(path, encoding='utf-8'))
    pool, assignment = load_tasks()
    who = book.get('annotator')
    if who not in assignment:
        raise SystemExit('unknown annotator %r in %s' % (who, path))
    return book, who, pool, assignment[who]


def check(path, quiet=False):
    book, who, pool, me = read_book(path)
    started, ready, problems = 0, [], {}
    for t in book.get('traces', []):
        task = pool.get(t.get('code'))
        if task is None:
            problems[t.get('code')] = ['unknown solution code']
            continue
        touched = any(s.get('label') for s in t.get('steps', [])) or t.get('final_answer')
        if not touched:
            continue
        started += 1
        bad = check_one(t, task)
        if bad:
            problems[t['code']] = bad
        else:
            ready.append((t, task))
    if not quiet:
        print('%s: %d of %d solutions started, %d complete and valid, %d with problems'
              % (who, started, len(book.get('traces', [])), len(ready), len(problems)))
        for code, bad in list(problems.items())[:20]:
            print('  %s' % code)
            for b in bad[:6]:
                print('     %s' % b)
        if len(problems) > 20:
            print('  ... and %d more' % (len(problems) - 20))
    return who, me, ready, problems


def import_(path):
    who, me, ready, problems = check(path)
    if not ready:
        print('nothing to import')
        return 1 if problems else 0
    pool, _ = load_tasks()
    os.makedirs(LABELS, exist_ok=True)
    with open(os.path.join(LABELS, '%s.jsonl' % who), 'a', encoding='utf-8', newline='\n') as fh:
        for t, task in ready:
            fh.write(json.dumps({
                'annotator': who, 'branch_of_annotator': me['branch'], 'code': t['code'],
                'branch': task['branch'], 'level': task['level'],
                'for_truth': t['code'] in me['for_truth'], 'calibration': t['code'] in me['calibration'],
                'item_sha256': task['item_sha256'], 'trace_sha256': task['trace_sha256'],
                'n_steps': len(task['steps']),
                'steps': [{'index': s['index'], 'label': s['label'],
                           'error_type': s.get('error_type') or None, 'note': s.get('note', '')}
                          for s in t['steps']],
                'milestones': [{'id': m['id'], 'value': m['value'], 'status': m['status']}
                               for m in t['milestones']],
                'final_answer': t['final_answer'], 'reasoning_sound': t['reasoning_sound'],
                'confidence': t['confidence'], 'comment': t.get('comment', ''),
                'flagged': bool(t.get('flagged', False)),
                'source': 'workbook', 'seconds': None,
                'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            }, ensure_ascii=False) + '\n')
    print('imported %d solutions into labels/%s.jsonl%s'
          % (len(ready), who, '; %d still have problems' % len(problems) if problems else ''))
    return 1 if problems else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('action', choices=('export', 'check', 'import'))
    ap.add_argument('path', nargs='?', help='the workbook, for check and import')
    ap.add_argument('--id', action='append', help='export only these annotators')
    a = ap.parse_args()
    if a.action == 'export':
        export(a.id)
        return 0
    if not a.path:
        raise SystemExit('%s needs the path to a workbook' % a.action)
    return check(a.path)[0] and 0 if a.action == 'check' else import_(a.path)


if __name__ == '__main__':
    sys.exit(main())
