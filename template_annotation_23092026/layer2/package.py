"""Build the bundle each expert is sent - and nothing more.

    python -m template_annotation_23092026.layer2.package            # every expert
    python -m template_annotation_23092026.layer2.package --id civ-2

One zip per expert in layer2/dist/:

    EngTrace-certification-guide.pdf    the guide
    HOW-TO-RUN.txt                      three commands
    app/app.py                          the app
    app/tasks/pool.json                 ONLY this expert's items
    app/tasks/assignment.json           ONLY this expert's queue
    <id>.json                           the workbook route, for those who prefer a file

Deliberately not in it: tasks/keyfile.jsonl (which items are plants), other experts'
queues, any labels, and our tooling. The packager refuses to build a bundle that
contains a keyfile or a plant id.
"""
from __future__ import annotations

import argparse
import json
import os
import zipfile
from pathlib import Path

from . import workbooks

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
GUIDE = HERE / 'EngTrace-certification-guide.pdf'

HOW_TO_RUN = """EngTrace template certification - how to run

Read EngTrace-certification-guide.pdf first.

ROUTE A - the app (recommended)
  1. Install Python 3.11 or newer, then:   pip install streamlit
  2. In this folder:                       streamlit run app/app.py
  3. Pick your id ({id}) in the sidebar. Work saves as you submit each item; you can
     close the app and come back. Your answers are written to app/labels/{id}.jsonl -
     send that file back when you finish.

ROUTE B - the workbook
  Open {id}.json in any text or JSON editor. For each item, read the question, write
  your answer in "hand_answer", then read the solution and fill the scores (1-5),
  "decision" (Approve or Reject), "defects" (from the list at the top of the file),
  "feedback" and "confidence". Send the file back; partly filled is fine.
"""


def build_one(aid: str, pool: dict, assignment: dict) -> Path:
    mine = assignment[aid]
    sub_pool = {c: pool[c] for c in mine['codes']}
    for c, item in sub_pool.items():
        blob = json.dumps(item)
        if 'plant_' in blob or 'keyfile' in blob:
            raise SystemExit(f'refusing: item {c} would reveal a plant')
    DIST.mkdir(exist_ok=True)
    out = DIST / f'{aid}.zip'
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(GUIDE, 'EngTrace-certification-guide.pdf')
        z.writestr('HOW-TO-RUN.txt', HOW_TO_RUN.format(id=aid))
        z.write(HERE / 'app.py', 'app/app.py')
        z.writestr('app/tasks/pool.json', json.dumps(sub_pool, ensure_ascii=False))
        z.writestr('app/tasks/assignment.json', json.dumps({aid: mine}, indent=1))
        z.writestr(f'{aid}.json', json.dumps(workbooks.export_one(aid, sub_pool, mine), ensure_ascii=False, indent=1))
        names = z.namelist()
    if any('keyfile' in n for n in names):
        out.unlink()
        raise SystemExit('refusing: bundle contains a keyfile')
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--id', default=None)
    a = ap.parse_args()
    if not GUIDE.exists():
        raise SystemExit('typeset the guide first: python -m template_annotation_23092026.layer2.make_guide_pdf')
    pool = json.loads((TASKS / 'pool.json').read_text(encoding='utf8'))
    assignment = json.loads((TASKS / 'assignment.json').read_text(encoding='utf8'))
    ids = [a.id] if a.id else sorted(assignment)
    for aid in ids:
        out = build_one(aid, pool, assignment)
        print(f'{aid}: {out.name} ({os.path.getsize(out) // 1024} KB, {len(assignment[aid]["codes"])} items)')


if __name__ == '__main__':
    main()
