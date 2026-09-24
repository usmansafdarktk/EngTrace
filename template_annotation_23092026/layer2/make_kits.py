"""The annotator kits: shared files once, one folder of items per expert, no zips.

    python -m template_annotation_23092026.layer2.make_kits            # every expert
    python -m template_annotation_23092026.layer2.make_kits --id civ-2

Writes dist/:

    app.py                             the review app, shared; run with: streamlit run app.py
    README.txt                         how to run, where the results go, shared
    EngTrace-certification-guide.pdf   the guide, shared
    guide.md                           the same text as Markdown, shared
    kit_<id>/tasks/pool.json           ONLY that expert's 34 items, instances precomputed
    kit_<id>/tasks/assignment.json     ONLY that expert's queue

An expert receives the four shared files plus their own kit_<id>/ folder, side by side.
The app finds every kit_<id>/tasks/ next to it and offers those ids; it writes <id>.jsonl
beside app.py as the expert submits, and that file is what they send back (it goes into
layer2/labels/ for score.py). No template code ships, and no kit contains the keyfile:
the builder checks each payload for a plant marker and refuses if one is found.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
GUIDE = HERE / 'EngTrace-certification-guide.pdf'

README = """EngTrace template certification - how to run

You should have, side by side in one folder: app.py, this README, the guide (PDF and
guide.md), and a folder named kit_<your id> holding your items.

1. Read EngTrace-certification-guide.pdf (three pages; guide.md is the same text).

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id in the sidebar.

3. Work through your items. The app saves after every item. To pause, close the browser
   tab and stop the app; run "streamlit run app.py" again to continue where you left off.

4. When the app says you are done, send the file <your id>.jsonl, which the app has
   written in this folder next to app.py, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers.
"""


def build_one(aid: str, pool: dict, assignment: dict) -> Path:
    mine = assignment[aid]
    sub = {c: pool[c] for c in mine['codes']}
    blob = json.dumps(sub, ensure_ascii=False)
    if 'plant_' in blob or 'keyfile' in blob:
        raise SystemExit(f'refusing: {aid} payload would reveal a plant')
    kit = DIST / f'kit_{aid}' / 'tasks'
    if kit.parent.exists():
        shutil.rmtree(kit.parent)
    kit.mkdir(parents=True)
    (kit / 'pool.json').write_text(blob, encoding='utf8')
    (kit / 'assignment.json').write_text(json.dumps({aid: mine}, indent=1), encoding='utf8')
    return kit.parent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--id', default=None)
    a = ap.parse_args()
    if not GUIDE.exists():
        raise SystemExit('typeset the guide first: python -m template_annotation_23092026.layer2.make_guide_pdf')
    pool = json.loads((TASKS / 'pool.json').read_text(encoding='utf8'))
    assignment = json.loads((TASKS / 'assignment.json').read_text(encoding='utf8'))
    DIST.mkdir(exist_ok=True)
    shutil.copy(HERE / 'app.py', DIST / 'app.py')
    shutil.copy(GUIDE, DIST / 'EngTrace-certification-guide.pdf')
    shutil.copy(HERE / 'guide.md', DIST / 'guide.md')
    (DIST / 'README.txt').write_text(README, encoding='utf8')
    for aid in ([a.id] if a.id else sorted(assignment)):
        out = build_one(aid, pool, assignment)
        print(f'{aid}: {out.relative_to(DIST)}/ ({len(assignment[aid]["codes"])} items)')
    print(f'shared: app.py, README.txt, EngTrace-certification-guide.pdf, guide.md in {DIST}')


if __name__ == '__main__':
    main()
