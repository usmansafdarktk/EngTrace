"""Assemble annotator_kit/: everything the annotators need, and nothing else.

    python evaluator_pilot_17092026/annotation/build_kit.py

The kit is standalone. It mentions no repository, no pilot folder and no evaluator:
an annotator can copy it anywhere, run two commands, and work. It contains the guide,
the app, the solutions, and one workbook per annotator. It never contains the keyfile
(the code -> model map), our tooling, or anyone's labels.
"""
from __future__ import annotations

import json
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = os.path.dirname(HERE)
TASKS = os.path.join(HERE, 'tasks')
KIT = os.path.join(PILOT, 'annotator_kit')

README = """EngTrace annotation
===================

Read the guide first: EngTrace-annotation-guide.pdf. It explains the task, the labels,
and the two ways of recording your answers. Use whichever suits you; you can switch.

ROUTE A - the app
  pip install -r requirements.txt
  streamlit run app.py
  Choose your name in the sidebar. Each solution is saved when you submit it, so you
  can close the app and come back. Your answers are written to labels/<your-id>.jsonl.
  Send that file back when you are done.

ROUTE B - the workbook
  Open workbooks/<your-id>.json in any text or JSON editor, fill in the empty fields
  as the guide describes, and send the file back. Partly filled is fine, and you can
  send it in parts.

Your id is the one you were given, for example civ-2.

Questions, or anything the guide does not settle: ask. Use the "flagged" field on a
solution so it can be found later.
"""

APP_HEADER = '''"""EngTrace step annotation.

    pip install -r requirements.txt
    streamlit run app.py

Pick your name in the sidebar and work through your list. Each solution is saved when
you submit it; you can close the app and return, and you can reopen a submitted
solution to change it. Answers are appended to labels/<your-id>.jsonl - that is the
file to send back.

The labels and what they mean are in EngTrace-annotation-guide.pdf.
"""
'''


def main():
    shutil.rmtree(KIT, ignore_errors=True)
    os.makedirs(os.path.join(KIT, 'tasks'))
    os.makedirs(os.path.join(KIT, 'workbooks'))

    # the app, with its internal docstring replaced so the kit stands alone
    src = open(os.path.join(HERE, 'app.py'), encoding='utf-8').read()
    body = src[src.index('"""', src.index('"""') + 3) + 3:].lstrip('\n')
    open(os.path.join(KIT, 'app.py'), 'w', encoding='utf-8', newline='\n').write(APP_HEADER + '\n' + body)

    for name in ('pool.json', 'assignment.json'):
        shutil.copy2(os.path.join(TASKS, name), os.path.join(KIT, 'tasks', name))
    books = sorted(f for f in os.listdir(os.path.join(TASKS, 'workbooks')) if f.endswith('.json'))
    for f in books:
        shutil.copy2(os.path.join(TASKS, 'workbooks', f), os.path.join(KIT, 'workbooks', f))
    shutil.copy2(os.path.join(HERE, 'EngTrace-annotation-guide.pdf'), KIT)
    open(os.path.join(KIT, 'README.txt'), 'w', encoding='utf-8', newline='\n').write(README)
    open(os.path.join(KIT, 'requirements.txt'), 'w', encoding='utf-8', newline='\n').write('streamlit>=1.40\n')

    # checks: no keyfile, nothing naming the model or item behind a code, no repo paths
    bad = []
    for root, _, files in os.walk(KIT):
        for f in files:
            p = os.path.join(root, f)
            if 'keyfile' in f.lower():
                bad.append('%s: keyfile' % f)
            if f.endswith(('.json', '.py', '.txt')):
                text = open(p, encoding='utf-8').read()
                for pat in ('"model_key"', '"item_id"', 'evaluator_pilot', 'EngTrace/',
                            'engtrace_evaluation', 'e2_prm', 'scores/'):
                    if pat in text:
                        bad.append('%s: mentions %s' % (f, pat))
    if bad:
        raise SystemExit('kit is not standalone:\n  ' + '\n  '.join(sorted(set(bad))))

    n = sum(len(fs) for _, _, fs in os.walk(KIT))
    size = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(KIT) for f in fs)
    print('annotator_kit/: %d files, %.1f MB, %d workbooks' % (n, size / 1e6, len(books)))
    for f in sorted(os.listdir(KIT)):
        print('  ' + f)


if __name__ == '__main__':
    main()
