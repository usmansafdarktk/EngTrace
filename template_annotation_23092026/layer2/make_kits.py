"""The annotator kit: one folder per expert with the app, their items, and the guide.

    python -m template_annotation_23092026.layer2.make_kits            # every expert
    python -m template_annotation_23092026.layer2.make_kits --id civ-2

Writes dist/kit_<id>/ and zips it to dist/kit_<id>.zip:

    app.py                             the review app; run with: streamlit run app.py
    tasks/pool.json                    ONLY this expert's 34 items, instances precomputed
    tasks/assignment.json              ONLY this expert's queue
    EngTrace-certification-guide.pdf   the guide
    guide.md                           the same text as Markdown
    README.txt                         how to run, where the results go

The app writes <id>.jsonl into the kit folder as the expert submits; that file is what
they send back, and it goes into layer2/labels/ for score.py. No template code ships,
and the kit never contains the keyfile: the builder checks the payload for any plant
marker and refuses if one is found.
"""
from __future__ import annotations

import argparse
import json
import shutil
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
GUIDE = HERE / 'EngTrace-certification-guide.pdf'

README = """EngTrace template certification - {id}

1. Read EngTrace-certification-guide.pdf (three pages; guide.md is the same text).

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id ({id}) in the sidebar.

3. Work through your {n} items. The app saves after every item. To pause, close the
   browser tab and stop the app; run "streamlit run app.py" again to continue where
   you left off.

4. When the app says you are done, send the file {id}.jsonl, which the app has
   written in this folder, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers.
"""


def build_one(aid: str, pool: dict, assignment: dict) -> Path:
    mine = assignment[aid]
    sub = {c: pool[c] for c in mine['codes']}
    blob = json.dumps(sub, ensure_ascii=False)
    if 'plant_' in blob or 'keyfile' in blob:
        raise SystemExit(f'refusing: {aid} payload would reveal a plant')
    kit = DIST / f'kit_{aid}'
    if kit.exists():
        shutil.rmtree(kit)
    (kit / 'tasks').mkdir(parents=True)
    shutil.copy(HERE / 'app.py', kit / 'app.py')
    (kit / 'tasks' / 'pool.json').write_text(blob, encoding='utf8')
    (kit / 'tasks' / 'assignment.json').write_text(json.dumps({aid: mine}, indent=1), encoding='utf8')
    shutil.copy(GUIDE, kit / 'EngTrace-certification-guide.pdf')
    shutil.copy(HERE / 'guide.md', kit / 'guide.md')
    (kit / 'README.txt').write_text(README.format(id=aid, n=len(mine['codes'])), encoding='utf8')
    out = DIST / f'kit_{aid}.zip'
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(kit.rglob('*')):
            if p.is_file():
                z.write(p, f'kit_{aid}/{p.relative_to(kit).as_posix()}')
        if any('keyfile' in n for n in z.namelist()):
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
    for aid in ([a.id] if a.id else sorted(assignment)):
        out = build_one(aid, pool, assignment)
        print(f'{aid}: {out.name} ({out.stat().st_size // 1024} KB, {len(assignment[aid]["codes"])} items)')


if __name__ == '__main__':
    main()
