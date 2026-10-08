"""EngTrace symbolic grading (E2) - the expert's app.

    streamlit run app.py

Finds the queues next to it, either tasks/ (the whole roster, as built in the repository) or one kit_<id>/tasks/
folder per expert (as distributed), lets the expert pick their id, and writes <id>.jsonl in the same folder as this
file, one row per submitted item: that file is what the expert sends back. Per item: the problem, the reference
answer, the model's final answer as the scorer reads it (answer.segment), shown as plain text as the model wrote
it, with a typeset view and the full response and reference solution one click away; for a reference that is one
value written as an expression, that value. One grade: equivalent, not equivalent, unreadable. The model, the
scorer's verdict and the template are not in the kit. The layout follows reading_app.py (the B1 to B4 readings).
Built by build_kit.py, scored by score_grades.py.
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os
import re

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LABELS = HERE                      # <id>.jsonl lands next to app.py, in the folder the expert runs it from
GRADES = ['equivalent', 'not equivalent', 'unreadable']

st.set_page_config(layout='wide', page_title='EngTrace symbolic grading')


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


def plain(label: str, text: str, key: str, height: int | None = None):
    """Text as written: a read-only box, so nothing in it is rendered."""
    lines = sum(len(line) // 80 + 1 for line in text.split('\n'))
    st.text_area(label, text, key=key, disabled=True, height=height or max(100, min(420, 26 * lines + 20)))


def typeset(text: str) -> str:
    """LaTeX delimiters \\( \\) and \\[ \\] as Streamlit's $ and $$, for a rendered view only."""
    t = re.sub(r'\\\[(.+?)\\\]', lambda m: '$$' + m.group(1) + '$$', text, flags=re.S)
    return re.sub(r'\\\((.+?)\\\)', lambda m: '$' + m.group(1) + '$', t, flags=re.S)


pool, assignment = load_tasks(HERE)

with st.sidebar:
    st.title('EngTrace: symbolic grading')
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
                            format_func=lambda c: f'{codes.index(c) + 1:3d}. {pool[c]["group"]}: {c}'
                                                  + ('  (done)' if c in done else ''))
        if st.session_state['code'] is not None and jump != st.session_state['code']:
            st.session_state.update(code=jump, opened_at=now())
            st.rerun()
    st.markdown('---')
    st.markdown('**Equivalent**: the same quantity as the reference, to the precision it shows (any equal form, a '
                'standard identity, exact where the reference rounds, another stated unit, a given quantity left as '
                'a symbol, or the value of a one-value reference). **Not equivalent**: a different quantity, a missing '
                'or extra term, factor or part. **Unreadable**: no single final answer.')

code = st.session_state['code']
if code is None:
    st.success(f'All {len(codes)} items graded. Send back the file `{aid}.jsonl` from the folder you ran the app in. '
               'Thank you.')
    st.stop()

item = pool[code]
idx = codes.index(code)
prev = done.get(code) or {}
st.subheader(f'Item {idx + 1} of {len(codes)} - problem family {item["group"]} - {code}')
if prev:
    st.info(f'Already submitted on {prev["submitted_at"]} ({prev["grade"]}). Submitting again replaces it.')

left, right = st.columns(2)
with left:
    plain('The problem', item['question'], f'q_{code}')
    plain('The reference answer', item['reference'], f'r_{code}', 140)
    if item.get('reference_value'):
        st.markdown(f"**Value of the reference:** {item['reference_value']}")
with right:
    plain("The model's final answer, as the scorer reads it", item['stated'], f'a_{code}', 260)
    with st.expander("The model's final answer, typeset (may misrender; the box above is what the model wrote)"):
        st.markdown(typeset(item['stated']))
with st.expander("The model's full response"):
    plain('Response', item['trace'], f't_{code}', 480)
with st.expander('The full reference solution'):
    plain('Reference solution', item['solution'], f's_{code}', 420)

st.subheader('Your grade')
grade = st.radio("Is the model's final answer equivalent to the reference?", GRADES,
                 index=GRADES.index(prev['grade']) if prev.get('grade') in GRADES else None, horizontal=True,
                 key=f'g_{code}')
note = st.text_area('Note (optional): a reference that looks wrong, or another acceptable approximation',
                    value=prev.get('note', ''), key=f'n_{code}', height=80)
if st.button('Submit and continue', type='primary', disabled=grade is None):
    row = {'annotator_id': aid, 'branch': branch, 'kind': 'symbolic', 'code': code, 'position': idx + 1,
           'opened_at': st.session_state['opened_at'], 'grade': grade, 'note': note.strip(),
           'submitted_at': now(), 'source': 'app'}
    with open(label_path(aid), 'a', encoding='utf8') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    done = done_codes(aid)
    remaining = [c for c in codes if c not in done]
    st.session_state.update(code=(remaining[0] if remaining else None), opened_at=now())
    st.rerun()
