"""EngTrace template certification - the reviewer's app.

    streamlit run app.py

Finds the queues next to it, either tasks/ (the whole roster, as built in the repository)
or one kit_<id>/tasks/ folder per expert (as distributed), lets the expert pick their id,
and writes <id>.jsonl in the same folder as this file, one row per submitted template:
that file is what the expert sends back. Runs no template code: instances are precomputed.

The flow per item: hand check (question only, enter your answer) -> show solution
(match reported) -> scores, decision, note -> submit. Timestamps for each stage are
recorded, because the December 2025 record could not show how long anyone looked.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import re

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(HERE, 'tasks')
LABELS = HERE                      # <id>.jsonl lands next to app.py, in the folder the expert runs it from
NUM_RE = re.compile(r'[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?')
DEFECTS = ['physics or scenario implausible', 'governing equation or formula', 'constant or table value',
           'unit or conversion', 'sign or direction', 'arithmetic: a step does not follow',
           'final answer wrong or wrong unit', 'question ambiguous or unsolvable', 'wording or formatting']

st.set_page_config(layout='wide', page_title='EngTrace template certification')


@st.cache_data
def load_tasks():
    """{id: queue} and the pool, from tasks/ if present, else from every kit_<id>/tasks/ beside app.py."""
    import glob
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
    os.makedirs(LABELS, exist_ok=True)
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


def numbers(text: str) -> list[float]:
    out = []
    for tok in NUM_RE.findall(text or ''):
        try:
            out.append(float(tok.replace(',', '')))
        except ValueError:
            pass
    return out


def hand_match(answer_text: str, gold: list[float]) -> bool | None:
    """True if any number the reviewer typed is within 1% of any gold number; None if no numbers."""
    mine = numbers(answer_text)
    if not mine or not gold:
        return None
    for m in mine:
        for g in gold:
            if g == 0 and abs(m) < 1e-9:
                return True
            if g != 0 and abs(m - g) / abs(g) <= 0.01:
                return True
    return False


pool, assignment = load_tasks()

with st.sidebar:
    st.title('EngTrace certification')
    aid = st.selectbox('Your reviewer id', [''] + sorted(assignment), index=0)
    if not aid:
        st.info('Pick your id to start. Read the guide first.')
        st.stop()
    branch = assignment[aid]['branch']
    codes = assignment[aid]['codes']
    done = done_codes(aid)
    st.caption(f'{branch.replace("_", " ")} - {len(done)} of {len(codes)} submitted')
    st.progress(len(done) / len(codes))
    remaining = [c for c in codes if c not in done]
    if 'code' not in st.session_state or st.session_state.get('aid') != aid:
        st.session_state.update(aid=aid, code=(remaining[0] if remaining else None), stage='hand',
                                instance=0, opened_at=now())
    pick = st.selectbox('Jump to an item', codes, index=codes.index(st.session_state['code']) if st.session_state['code'] in codes else 0,
                        format_func=lambda c: f'{codes.index(c) + 1:2d}. {c}' + ('  (done)' if c in done else ''))
    if pick != st.session_state['code']:
        st.session_state.update(code=pick, stage='hand', instance=0, opened_at=now())
        st.rerun()

code = st.session_state['code']
if code is None:
    st.success(f'All {len(codes)} items submitted. Send back the file `{aid}.jsonl` from the folder you ran the app in. Thank you.')
    st.stop()

item = pool[code]
idx = codes.index(code)
prev = done.get(code)
st.subheader(f'Item {idx + 1} of {len(codes)} - {code} - {item["area"].replace("_", " ")}')
if prev:
    st.info(f'Already submitted on {prev["submitted_at"]} ({prev["decision"]}). Submitting again replaces it.')

inst = item['instances'][st.session_state['instance']]

if st.session_state['stage'] == 'hand':
    st.markdown('### Hand check')
    st.markdown('Solve this instance yourself and enter the final answer(s) with unit before seeing the solution.')
    st.markdown(inst['question'])
    ans = st.text_area('Your answer(s)', key=f'hand_{code}', height=80,
                       placeholder='e.g. 24.05 kN (tension); or a classification in words')
    if st.button('Show solution', type='primary', disabled=not ans.strip()):
        st.session_state.update(stage='review', hand_answer=ans.strip(), hand_at=now())
        st.rerun()
    st.stop()

# ---- review stage
match = hand_match(st.session_state.get('hand_answer', ''), item['gold_numbers'])
if match is True:
    st.success('Your answer matches the template\'s answer within 1%.')
elif match is False:
    st.warning('Your answer differs from the template\'s by more than 1%. Look for which side is right.')
else:
    st.info('No number to compare (a non-numeric answer, or none entered). Judge the solution directly.')
st.caption(f'You entered: {st.session_state.get("hand_answer", "")}')

tab_sol, tab_src = st.tabs(['Question and solution', 'Source code'])
with tab_sol:
    c1, c2 = st.columns([1, 2])
    with c1:
        st.markdown('**Question**')
        st.markdown(inst['question'])
    with c2:
        st.markdown('**Solution**')
        st.markdown(inst['solution'])
    cols = st.columns([1, 3])
    with cols[0]:
        if st.button(f'Another instance ({st.session_state["instance"] + 1} of {len(item["instances"])})'):
            st.session_state['instance'] = (st.session_state['instance'] + 1) % len(item['instances'])
            st.rerun()
with tab_src:
    st.code(item['source'], language='python')

st.markdown('### Your verdict')
s1, s2, s3 = st.columns(3)
phys = s1.slider('Physical plausibility', 1, 5, prev['scores']['physical_plausibility'] if prev else 5, key=f'p_{code}')
math_ = s2.slider('Mathematical correctness', 1, 5, prev['scores']['mathematical_correctness'] if prev else 5, key=f'm_{code}')
ped = s3.slider('Pedagogical clarity', 1, 5, prev['scores']['pedagogical_clarity'] if prev else 5, key=f'c_{code}')
decision = st.radio('Decision', ['Approve', 'Reject'], horizontal=True,
                    index=0 if not prev or prev['decision'] == 'Approve' else 1, key=f'd_{code}')
defects = st.multiselect('Defect type (required on Reject)', DEFECTS, default=prev.get('defects', []) if prev else [],
                         key=f'df_{code}')
note = st.text_area('Note (required on Reject: what is wrong and where; optional on Approve)',
                    value=prev.get('feedback', '') if prev else '', key=f'n_{code}', height=90)
confidence = st.slider('Your confidence in this verdict', 1, 5, prev.get('confidence', 4) if prev else 4, key=f'cf_{code}')

blocked = decision == 'Reject' and (not note.strip() or not defects)
if blocked:
    st.caption('A rejection needs a note and at least one defect type.')
if st.button('Submit and continue', type='primary', disabled=blocked):
    row = {'annotator_id': aid, 'branch': branch, 'code': code, 'position': idx + 1,
           'opened_at': st.session_state['opened_at'], 'hand_answer': st.session_state.get('hand_answer', ''),
           'hand_numbers': numbers(st.session_state.get('hand_answer', '')), 'hand_match': match,
           'hand_submitted_at': st.session_state.get('hand_at'), 'hand_instance_seed': item['instances'][0]['seed'],
           'instances_viewed': st.session_state['instance'] + 1,
           'scores': {'physical_plausibility': phys, 'mathematical_correctness': math_, 'pedagogical_clarity': ped},
           'decision': decision, 'defects': defects, 'feedback': note.strip(), 'confidence': confidence,
           'submitted_at': now(), 'source': 'app'}
    with open(label_path(aid), 'a', encoding='utf8') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')
    done = done_codes(aid)
    remaining = [c for c in codes if c not in done]
    st.session_state.update(code=(remaining[0] if remaining else None), stage='hand', instance=0, opened_at=now())
    st.rerun()
