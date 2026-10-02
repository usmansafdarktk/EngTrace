"""EngTrace reading request, October 2026 - the expert's app.

    streamlit run app.py

Finds the queues next to it, tasks/ (the whole request, as built in the repository) or one kit_<id>/tasks/ folder
per expert (as distributed), lets the expert pick their id, and writes <id>.jsonl in the same folder as this file,
one row per submitted item: that file is what the expert sends back. Four kinds of item share the queue (the guide
names them): a final answer, a quantity in the working, why a wrong answer went wrong, and questions about a
template. Everything is shown as plain text, as the model wrote or read it, never rendered as Markdown. Opening and
submitting are timestamped. Built and scored by expert_kits.py.
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

TITLE = {'template': 'Questions about a template', 'answer': 'A final answer',
         'milestone': 'A quantity in the working', 'error': 'Why a wrong answer went wrong'}
ANSWER_OPTIONS = ['correct', 'partially correct', 'incorrect', 'no answer stated']
MILESTONE_OPTIONS = ['yes, the working obtains it', 'no, it never obtains it', 'no, and its route does not need it']
ERROR_OPTIONS = ['1. Hallucination', '2. Setup / Assumption Error', '3. Formula / Principle Error',
                 '4. Unit / Dimensional Error', '5. Sign / Direction Error', '6. Calculation Error',
                 'No error: the answer is correct, or the question admits it',
                 'Incomplete: the working stops before an answer']
NOT_AN_ERROR = ERROR_OPTIONS[6:]

st.set_page_config(layout='wide', page_title='EngTrace reading request')


@st.cache_data
def load_tasks():
    """{code: task} and {id: queue}, from tasks/ if present, else from every kit_<id>/tasks/ beside app.py."""
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
    st.text_area(label, text, key=key, disabled=True, height=height or max(120, min(560, 26 * lines + 20)))


def pick(label: str, options: list[str], prev_value, key: str, horizontal: bool = False):
    index = options.index(prev_value) if prev_value in options else None
    return st.radio(label, options, index=index, horizontal=horizontal, key=key)


pool, assignment = load_tasks()

with st.sidebar:
    st.title('EngTrace reading request')
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
                            format_func=lambda c: f'{codes.index(c) + 1:2d}. {TITLE[pool[c]["kind"]]}'
                                                  + ('  (done)' if c in done else ''))
        if st.session_state['code'] is not None and jump != st.session_state['code']:
            st.session_state.update(code=jump, opened_at=now())
            st.rerun()

code = st.session_state['code']
if code is None:
    st.success(f'All {len(codes)} items submitted. Send back the file `{aid}.jsonl` from the folder you ran '
               'the app in. Thank you.')
    st.stop()

task = pool[code]
kind = task['kind']
idx = codes.index(code)
prev = done.get(code) or {}
pa = prev.get('answers', {})
st.subheader(f'Item {idx + 1} of {len(codes)}: {TITLE[kind]}  ({code})')
if prev:
    st.info(f'Already submitted on {prev["submitted_at"]}. Submitting again replaces it.')

answers: dict = {}
excerpt = ''
note_required = False
note_prompt = 'Note (optional)'

if kind == 'answer':
    left, right = st.columns(2)
    with left:
        plain('The problem', task['question'], f'q_{code}')
    with right:
        plain('The reference answer', task['reference'], f'r_{code}', 160)
        plain("The model's final answer", task['stated'], f'a_{code}', 200)
    with st.expander("The model's full working"):
        plain('Working', task['trace'], f't_{code}', 480)
    with st.expander('The full reference solution'):
        plain('Reference solution', task['solution'], f's_{code}', 420)
    st.subheader('Your verdict')
    answers['verdict'] = pick("Is the model's final answer correct?", ANSWER_OPTIONS, pa.get('verdict'),
                              f'v_{code}', horizontal=True)
    note_prompt = ('Note: what is wrong or missing, if anything; and whether you think the reference is wrong or the '
                   'question allows this answer too.')
    ready = answers['verdict'] is not None

elif kind == 'milestone':
    plain('The problem', task['question'], f'q_{code}')
    st.markdown(f"**The quantity to look for:** `{task['milestone_id']}` = `{task['milestone_value']}`  "
                '(the name is the one the reference solution uses; the reference is one click below)')
    with st.expander('The reference solution'):
        plain('Reference solution', task['solution'], f's_{code}', 420)
    plain("The model's working", task['trace'], f't_{code}', 520)
    st.subheader('Your verdict')
    answers['obtained'] = pick('Does the working obtain this quantity?', MILESTONE_OPTIONS, pa.get('obtained'),
                               f'v_{code}')
    note_required = answers['obtained'] == MILESTONE_OPTIONS[2]
    note_prompt = 'Note: required for "does not need it" (say which route); otherwise optional.'
    ready = answers['obtained'] is not None

elif kind == 'error':
    plain('The problem', task['question'], f'q_{code}')
    plain("The model's working (its final answer was scored wrong)", task['trace'], f't_{code}', 560)
    with st.expander('The reference solution'):
        plain('Reference solution', task['solution'], f's_{code}', 420)
    st.subheader('Your verdict')
    answers['category'] = pick('Reading from the top, the first question you answer yes to:', ERROR_OPTIONS,
                               pa.get('category'), f'v_{code}')
    excerpt = st.text_area('Excerpt: the shortest piece of the working that shows the error (required for the six '
                           'error categories)', value=prev.get('excerpt', ''), key=f'e_{code}', height=90)
    note_required = answers['category'] == NOT_AN_ERROR[0]
    note_prompt = 'Note: required for "No error" (say what the question leaves open, or what is wrong with the reference).'
    ready = answers['category'] is not None and (answers['category'] in NOT_AN_ERROR or excerpt.strip())
    if answers['category'] is not None and answers['category'] not in NOT_AN_ERROR and not excerpt.strip():
        st.caption('An error category needs an excerpt.')

else:  # template
    plain('What the run found for this template', task['context'], f'c_{code}', 150)
    left, right = st.columns(2)
    with left:
        plain('The problem, as posed', task['question'], f'q_{code}')
    with right:
        plain('The reference solution', task['solution'], f's_{code}', 420)
    for k, tr in enumerate(task['traces']):
        with st.expander(f'Model answer {"ABC"[k]}' + ('' if len(task['traces']) == 1 else f' of {len(task["traces"])}')):
            plain(f'Model answer {"ABC"[k]}', tr, f't{k}_{code}', 480)
    st.subheader('Your reading')
    for ask in task['asks']:
        answers[ask['key']] = pick(ask['text'], ask['options'], pa.get(ask['key']), f"{ask['key']}_{code}")
    note_prompt = task.get('note_prompt') or 'Note (optional)'
    ready = all(v is not None for v in answers.values())

note = st.text_area(note_prompt, value=prev.get('note', ''), key=f'n_{code}', height=90)
if note_required and not note.strip():
    st.caption('This verdict needs a note.')
    ready = False

if st.button('Submit and continue', type='primary', disabled=not ready):
    row = {'annotator_id': aid, 'branch': branch, 'kind': kind, 'code': code, 'position': idx + 1,
           'opened_at': st.session_state['opened_at'], 'answers': answers, 'excerpt': excerpt.strip(),
           'note': note.strip(), 'submitted_at': now(), 'source': 'app'}
    with open(label_path(aid), 'a', encoding='utf8') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    done = done_codes(aid)
    remaining = [c for c in codes if c not in done]
    st.session_state.update(code=(remaining[0] if remaining else None), opened_at=now())
    st.rerun()
