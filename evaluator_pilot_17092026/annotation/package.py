"""Build the bundle each annotator is sent - and nothing more.

    python evaluator_pilot_17092026/annotation/package.py [--id civ-2]

One zip per annotator, in dist/:

    EngTrace-annotation-guide.pdf   the guide (both routes are described in it)
    <id>.json                       their workbook, to fill in any editor
    app/app.py                      the app, for whoever prefers it
    app/tasks/pool.json             ONLY their 68 solutions
    app/tasks/assignment.json       ONLY their assignment
    HOW-TO-RUN.txt                  two commands

WHAT IS DELIBERATELY NOT IN IT
  keyfile.jsonl      the code -> item + model map. Sending it would undo the
                     blinding, so this script refuses to build a bundle containing
                     a keyfile, and checks every file it adds.
  other annotators   pool.json is filtered to the codes assigned to this person.
  our tooling        build_tasks.py, workbooks.py, score_against_labels.py,
                     annotators.json, anyone's labels.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(HERE, 'tasks')
BOOKS = os.path.join(TASKS, 'workbooks')
DIST = os.path.join(HERE, 'dist')
GUIDE = os.path.join(HERE, 'EngTrace-annotation-guide.pdf')

HOW_TO_RUN = """EngTrace annotation - how to run

Read EngTrace-annotation-guide.pdf first. It explains the task, the labels and both
ways of recording your answers. Either way is fine; you can also switch between them.

ROUTE A - the app
  1. Install Python 3.11 or newer, then:  pip install streamlit
  2. In this folder:                      streamlit run app/app.py
  3. Choose your name in the sidebar. Work saves as you submit each solution; you can
     close it and come back. Your answers are written to app/labels/{id}.jsonl -
     send that file back.

ROUTE B - the workbook
  Open {id}.json in any text or JSON editor, fill in the empty fields as the guide
  describes, and send it back. Partly filled is fine, and you can send it in parts.

Questions, or anything the guide does not settle: ask, and use the "flagged" field so
we can find the solution later.
"""


def build(who, me, pool):
    out = os.path.join(DIST, who)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(os.path.join(out, 'app', 'tasks'))

    mine = {c: pool[c] for c in me['order']}
    with open(os.path.join(out, 'app', 'tasks', 'pool.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(mine, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(out, 'app', 'tasks', 'assignment.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump({who: me}, fh, indent=1)
    shutil.copy2(os.path.join(HERE, 'app.py'), os.path.join(out, 'app', 'app.py'))
    shutil.copy2(GUIDE, os.path.join(out, os.path.basename(GUIDE)))
    book = os.path.join(BOOKS, '%s.json' % who)
    if not os.path.exists(book):
        raise SystemExit('no workbook for %s - run: python annotation/workbooks.py export' % who)
    shutil.copy2(book, os.path.join(out, '%s.json' % who))
    with open(os.path.join(out, 'HOW-TO-RUN.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(HOW_TO_RUN.replace('{id}', who))

    files = [os.path.join(r, f) for r, _, fs in os.walk(out) for f in fs]
    for f in files:                                   # the blinding check, on what is really there
        if 'keyfile' in os.path.basename(f).lower():
            raise SystemExit('REFUSING: %s would go to the annotator' % f)
        if f.endswith('.json'):
            if any(k in open(f, encoding='utf-8').read() for k in ('"model_key"', '"item_id"')):
                raise SystemExit('REFUSING: %s names the model or item behind a code' % f)
    zpath = os.path.join(DIST, '%s.zip' % who)
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in sorted(files):
            z.write(f, os.path.relpath(f, out))
    return zpath, len(mine)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--id', action='append', help='only these annotators')
    a = ap.parse_args()
    pool = json.load(open(os.path.join(TASKS, 'pool.json'), encoding='utf-8'))
    assignment = json.load(open(os.path.join(TASKS, 'assignment.json'), encoding='utf-8'))
    os.makedirs(DIST, exist_ok=True)
    for who, me in sorted(assignment.items()):
        if a.id and who not in a.id:
            continue
        z, n = build(who, me, pool)
        print('%-8s %3d solutions  %6.0f KB  %s' % (who, n, os.path.getsize(z) / 1024, z))
    print('\nEach bundle holds the guide, the annotator\'s workbook, and the app with only '
          'their own solutions. No keyfile, no other annotator\'s work.')


if __name__ == '__main__':
    main()
