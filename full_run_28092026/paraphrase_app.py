"""EngTrace paraphrase check - the reviewer's app.

    streamlit run app.py

Finds the queues next to it, tasks/ (the whole roster, as built in the repository) or one
kit_<id>/tasks/ folder per expert (as distributed), lets the expert pick their id, and writes
<id>.jsonl in the same folder as this file, one row per submitted item: that file is what the
expert sends back. Every question, paraphrase and solution is shown as plain text, as a model reads
it, never rendered as Markdown. Opening and submitting are timestamped.
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(HERE, 'tasks')
LABELS = HERE                      # <id>.jsonl lands next to app.py, in the folder the expert runs it from

st.set_page_config(layout='wide', page_title='EngTrace paraphrase check')


@st.cache_data
def load_tasks():
    """{id: queue} and the items, from tasks/ if present, else from every kit_<id>/tasks/ beside app.py."""
    dirs = [TASKS] if os.path.exists(os.path.join(TASKS, 'assignment.json')) else \
        sorted(os.path.dirname(p) for p in glob.glob(os.path.join(HERE, 'kit_*', 'tasks', 'assignment.json')))
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
    st.text_area(label, text, key=key, disabled=True, height=height or max(120, min(520, 26 * lines + 20)))


def pick(label: str, options: list[str], prev_value, key: str):
    index = options.index(prev_value) if prev_value in options else None
    return st.radio(label, options, index=index, horizontal=True, key=key)


pool, assignment = load_tasks()

with st.sidebar:
    st.title('EngTrace paraphrase check')
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

code = st.session_state['code']
if code is None:
    st.success(f'All {len(codes)} items submitted. Send back the file `{aid}.jsonl` from the folder you ran '
               'the app in. Thank you.')
    st.stop()

item = pool[code]
idx = codes.index(code)
prev = done.get(code) or {}
st.subheader(f'Item {idx + 1} of {len(codes)} - {code}')
if prev:
    st.info(f'Already submitted on {prev["submitted_at"]}. Submitting again replaces it.')

left, right = st.columns(2)
with left:
    plain('Original', item['original'], f'o_{code}')
with right:
    plain('Rewritten', item['paraphrase'], f'p_{code}')
plain("The original's answer", item['gold_answer'], f'g_{code}', 150)
with st.expander("The original's full solution"):
    plain('Solution', item['solution'], f's_{code}', 420)

st.subheader('Your verdict')
same = pick('1. Is the rewritten version the same problem: the same givens, the same quantities asked for, '
            'the same conditions?', ['yes', 'no'], prev.get('same'), f'q1_{code}')
answer = pick("2. Does the original's answer still answer the rewritten version, exactly?",
              ['yes', 'no'], prev.get('answer'), f'q2_{code}')
clear = pick('3. Does the rewritten version add an ambiguity, an error or a hint that the original does not have?',
             ['no', 'yes'], prev.get('clear'), f'q3_{code}')
note = st.text_area('Note: required unless your answers are yes, yes, no. What changed, and where.',
                    value=prev.get('note', ''), key=f'n_{code}', height=90)

keep = same == 'yes' and answer == 'yes' and clear == 'no'
ready = None not in (same, answer, clear) and (keep or note.strip())
if None not in (same, answer, clear) and not keep and not note.strip():
    st.caption('A pair you do not keep needs a note.')
if st.button('Submit and continue', type='primary', disabled=not ready):
    row = {'annotator_id': aid, 'branch': branch, 'code': code, 'position': idx + 1,
           'opened_at': st.session_state['opened_at'], 'same': same, 'answer': answer, 'clear': clear,
           'note': note.strip(), 'submitted_at': now(), 'source': 'app'}
    with open(label_path(aid), 'a', encoding='utf8') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    done = done_codes(aid)
    remaining = [c for c in codes if c not in done]
    st.session_state.update(code=(remaining[0] if remaining else None), opened_at=now())
    st.rerun()
