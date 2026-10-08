"""EngTrace difficulty ratings (E1) - the expert's app.

    streamlit run app.py

Finds the queues next to it, either tasks/ (the whole roster, as built in the repository) or one kit_<id>/tasks/
folder per expert (as distributed), lets the expert pick their id, and writes <id>.jsonl in the same folder as this
file, one row per submitted template: that file is what the expert sends back. Per template it shows the domain and
area, one problem drawn from the template and its reference solution, the rubric, and two questions: the level
(Easy, Intermediate, Advanced) and whether the problem states the governing formula or names the method. The
template's current label is not in the kit. The layout follows the certification app of layer 2, which the experts
used in rounds 1 to 6. Built by build_kit.py, scored by score_levels.py.
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LABELS = HERE                      # <id>.jsonl lands next to app.py, in the folder the expert runs it from
LEVELS = ['Easy', 'Intermediate', 'Advanced']
FORMULA = ['yes', 'no']
RUBRIC = """**Conceptual complexity**: one isolated principle (Easy) to the synthesis of several principles or domains (Advanced).

**Mathematical sophistication**: direct algebraic substitution (Easy) to differential equations, or an iterative or implicit solution (Advanced).

**Procedural depth**: a short chain of steps (Easy) to a long chain of interdependent steps (Advanced).

**Easy**: low on all three. **Intermediate**: in between. **Advanced**: a synthesis of principles or domains, a differential equation or an iterative or implicit solution, or a long chain of interdependent steps.

Judge the template, not the particular numbers."""

st.set_page_config(layout='wide', page_title='EngTrace difficulty ratings')


@st.cache_data
def load_tasks(root: str):
    """{code: item} and {id: queue}, from root/tasks/ if present, else from every root/kit_<id>/tasks/."""
    tasks = os.path.join(root, 'tasks')
    dirs = [tasks] if os.path.exists(os.path.join(tasks, 'assignment.json')) else \
        sorted(os.path.dirname(p) for p in glob.glob(os.path.join(root, 'kit_*', 'tasks', 'assignment.json')))
    pool, assignment = {}, {}
    for d in dirs:
        with open(os.path.join(d, 'pool.json'), encoding='utf8') as fh:
            pool.update(json.load(fh))
        with open(os.path.join(d, 'assignment.json'), encoding='utf8') as fh:
            assignment.update(json.load(fh))
    if not assignment:
        st.error('No tasks found next to app.py: expected tasks/ or kit_<id>/tasks/.')
        st.stop()
    return pool, assignment


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def label_path(aid: str) -> str:
    return os.path.join(LABELS, f'{aid}.jsonl')


def done_codes(aid: str) -> dict:
    out = {}
    p = label_path(aid)
    if os.path.exists(p):
        for ln in open(p, encoding='utf8'):
            if ln.strip():
                r = json.loads(ln)
                out[r['code']] = r
    return out


def pick(label: str, options: list[str], prev_value, key: str):
    index = options.index(prev_value) if prev_value in options else None
    return st.radio(label, options, index=index, horizontal=True, key=key)


pool, assignment = load_tasks(HERE)

with st.sidebar:
    st.title('EngTrace: difficulty ratings')
    aid = st.selectbox('Your reviewer id', [''] + sorted(assignment), index=0)
    if not aid:
        st.info('Pick your id to start. Read guide.md first.')
        st.stop()
    branch = assignment[aid]['branch']
    codes = assignment[aid]['codes']
    done = done_codes(aid)
    st.caption(f'{branch.replace("_", " ")} - {len(done)} of {len(codes)} submitted')
    st.progress(len(done) / len(codes) if codes else 1.0)
    remaining = [c for c in codes if c not in done]
    if 'code' not in st.session_state or st.session_state.get('aid') != aid:
        st.session_state.update(aid=aid, code=(remaining[0] if remaining else None), opened_at=now())
    if codes:
        at = codes.index(st.session_state['code']) if st.session_state['code'] in codes else 0
        jump = st.selectbox('Jump to an item', codes, index=at,
                            format_func=lambda c: f'{codes.index(c) + 1:2d}. {c}' + ('  (done)' if c in done else ''))
        if st.session_state['code'] is not None and jump != st.session_state['code']:
            st.session_state.update(code=jump, opened_at=now())
            st.rerun()
    st.markdown('---')
    st.markdown('**The rubric**')
    st.markdown(RUBRIC)

code = st.session_state['code']
if code is None:
    st.success(f'All {len(codes)} templates rated. Send back the file `{aid}.jsonl` from the folder you ran the app '
               'in. Thank you.')
    st.stop()

item = pool[code]
idx = codes.index(code)
prev = done.get(code) or {}
st.subheader(f'Template {idx + 1} of {len(codes)} - {code} - {item["domain"].replace("_", " ")} / '
             f'{item["area"].replace("_", " ")}')
if prev:
    st.info(f'Already submitted on {prev["submitted_at"]} ({prev["level"]}). Submitting again replaces it.')

c1, c2 = st.columns([1, 2])
with c1:
    st.markdown('**The problem**')
    st.markdown(item['question'])
with c2:
    st.markdown('**The reference solution**')
    st.markdown(item['solution'])
with st.expander('The same text, unformatted'):
    st.text_area('Problem', item['question'], key=f'qp_{code}', disabled=True, height=160)
    st.text_area('Reference solution', item['solution'], key=f'sp_{code}', disabled=True, height=360)

st.markdown('### Your rating')
level = pick('Level of this template', LEVELS, prev.get('level'), f'l_{code}')
formula = pick('Does the problem state the governing formula or name the method?', FORMULA,
               prev.get('formula_stated'), f'f_{code}')
note = st.text_area('Note (optional): anything that made the choice hard', value=prev.get('note', ''),
                    key=f'n_{code}', height=80)
ready = level is not None and formula is not None
if not ready:
    st.caption('Choose a level and answer the formula question to submit.')
if st.button('Submit and continue', type='primary', disabled=not ready):
    row = {'annotator_id': aid, 'branch': branch, 'code': code, 'position': idx + 1,
           'opened_at': st.session_state['opened_at'], 'level': level, 'formula_stated': formula,
           'note': note.strip(), 'submitted_at': now(), 'source': 'app'}
    with open(label_path(aid), 'a', encoding='utf8') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    done = done_codes(aid)
    remaining = [c for c in codes if c not in done]
    st.session_state.update(code=(remaining[0] if remaining else None), opened_at=now())
    st.rerun()
