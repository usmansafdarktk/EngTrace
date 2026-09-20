"""EngTrace step annotation.

    pip install -r requirements.txt
    streamlit run app.py

Pick your name in the sidebar and work through your list. Each solution is saved when
you submit it; you can close the app and return, and you can reopen a submitted
solution to change it. Answers are appended to labels/<your-id>.jsonl - that is the
file to send back.

The labels and what they mean are in EngTrace-annotation-guide.pdf.
"""

from __future__ import annotations

import html
import json
import os
import time

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(HERE, 'tasks')
LABELS = os.path.join(HERE, 'labels')

STEP_LABELS = [
    ('correct', 'Correct', 'The step is right: what it states and computes holds.'),
    ('alternative_correct', 'Correct, different route',
     'Right, but not how the reference solution does it. Judge it on its own merits.'),
    ('incorrect', 'Incorrect', 'Something in the step is wrong. Choose the kind below.'),
    ('not_a_claim', 'No claim to check',
     'Restating the question, announcing a plan, or otherwise making no checkable claim.'),
]
ERROR_TYPES = [
    ('calculation', 'Calculation error - the method is right, the arithmetic is not'),
    ('conceptual', 'Conceptual error - wrong equation, wrong principle, wrong assumption'),
    ('unsupported', 'Unsupported - a value or claim that is asserted with no basis'),
]
MILESTONE_CHOICES = [
    ('reached', 'Obtained'),
    ('not_reached', 'Not obtained'),
    ('not_needed', 'Not needed on this route'),
]
FINAL_CHOICES = [
    ('correct', 'Correct'),
    ('incorrect', 'Incorrect'),
    ('partial', 'Partly correct (multipart answer)'),
    ('not_stated', 'No final answer stated'),
]
CONFIDENCE = [('high', 'High'), ('medium', 'Medium'), ('low', 'Low')]

CSS = """
<style>
  .block-container {max-width: 1100px; padding-top: 2rem;}
  h1, h2, h3 {letter-spacing: -0.01em;}
  .meta {color: #5b6472; font-size: 0.86rem;}
  .step-head {font-weight: 600; margin-bottom: 0.15rem;}
  .step-body {white-space: pre-wrap; font-family: ui-monospace, SFMono-Regular, Consolas, monospace;
              font-size: 0.86rem; background: #f6f7f9; border: 1px solid #e3e6ea;
              border-radius: 6px; padding: 0.7rem 0.85rem; margin-bottom: 0.5rem;}
  .question {white-space: pre-wrap; background: #ffffff; border-left: 3px solid #2f4f7f;
             padding: 0.6rem 0.9rem; margin-bottom: 0.8rem;}
  .done {color: #1f7a3d;} .todo {color: #8a6d1f;}
  div[data-testid="stRadio"] label p {font-size: 0.88rem;}
</style>
"""


def load_json(path):
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


@st.cache_data
def load_tasks():
    return load_json(os.path.join(TASKS, 'pool.json')), load_json(os.path.join(TASKS, 'assignment.json'))


def labels_path(annotator):
    return os.path.join(LABELS, '%s.jsonl' % annotator)


def read_labels(annotator):
    """Last row per code wins, so a resubmission supersedes."""
    out = {}
    p = labels_path(annotator)
    if os.path.exists(p):
        for line in open(p, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                out[r['code']] = r
    return out


def append_label(annotator, row):
    os.makedirs(LABELS, exist_ok=True)
    with open(labels_path(annotator), 'a', encoding='utf-8', newline='\n') as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + '\n')


def radio(key, options, default=None, label=None, horizontal=False, help_=None):
    """A radio with no visible label (the surrounding text is the label)."""
    keys = [k for k, *_ in options]
    idx = keys.index(default) if default in keys else None
    return st.radio(label or key, keys, index=idx, key=key, horizontal=horizontal, help=help_,
                    format_func=lambda k: dict((o[0], o[1]) for o in options)[k],
                    label_visibility='visible' if label else 'collapsed')


def main():
    st.set_page_config(page_title='EngTrace annotation', layout='wide')
    st.markdown(CSS, unsafe_allow_html=True)
    pool, assignment = load_tasks()

    with st.sidebar:
        st.markdown('### EngTrace annotation')
        who = st.selectbox('Annotator', sorted(assignment), format_func=lambda k: '%s - %s' % (
            k, assignment[k]['name'] or assignment[k]['branch'].replace('_', ' ')))
        me = assignment[who]
        done = read_labels(who)
        order = me['order']
        n_done = sum(1 for c in order if c in done)
        st.progress(n_done / len(order), text='%d of %d complete' % (n_done, len(order)))
        st.markdown('<span class="meta">Branch: %s</span>' % me['branch'].replace('_', ' '),
                    unsafe_allow_html=True)

        remaining = [c for c in order if c not in done]
        if 'cursor' not in st.session_state or st.session_state.get('who') != who:
            st.session_state.cursor = order.index(remaining[0]) if remaining else 0
            st.session_state.who = who
        choice = st.selectbox(
            'Go to trace', order, index=st.session_state.cursor,
            format_func=lambda c: '%s  %s%s' % (
                c, 'done' if c in done else 'to do',
                '  (calibration)' if c in me['calibration'] else ''))
        st.session_state.cursor = order.index(choice)
        st.markdown('---')
        st.markdown('<span class="meta">The reference solution is available on each '
                    'trace, collapsed. A different valid route is not an error - use '
                    '"Correct, different route".</span>', unsafe_allow_html=True)

    code = order[st.session_state.cursor]
    task = pool[code]
    prior = done.get(code)

    st.markdown('## Trace %s' % code)
    st.markdown('<span class="meta">%s &nbsp;|&nbsp; %s &nbsp;|&nbsp; %d steps &nbsp;|&nbsp; '
                '%s</span>' % (task['branch'].replace('_', ' '), task['level'], len(task['steps']),
                               'calibration trace' if code in me['calibration'] else
                               ('your branch' if task['branch'] == me['branch'] else 'other branch')),
                unsafe_allow_html=True)
    if prior:
        st.info('You submitted this trace on %s. Re-submitting replaces that answer.' % prior['ts'])

    st.markdown('### Problem')
    st.markdown('<div class="question">%s</div>' % html.escape(task['question']),
                unsafe_allow_html=True)
    with st.expander('Reference solution (open only if you need it)'):
        st.markdown('<div class="step-body">%s</div>' % html.escape(task['reference_solution']),
                    unsafe_allow_html=True)

    if 'started' not in st.session_state or st.session_state.get('started_code') != code:
        st.session_state.started = time.time()
        st.session_state.started_code = code

    st.markdown('### Steps')
    step_rows = []
    for i, text in enumerate(task['steps']):
        p = (prior or {}).get('steps', [{}] * len(task['steps']))
        p = p[i] if i < len(p) else {}
        st.markdown('<div class="step-head">Step %d</div>' % (i + 1), unsafe_allow_html=True)
        st.markdown('<div class="step-body">%s</div>' % html.escape(text), unsafe_allow_html=True)
        c1, c2 = st.columns([1, 1])
        with c1:
            lab = radio('%s_s%d_lab' % (code, i), STEP_LABELS, p.get('label'))
        with c2:
            err = None
            if lab == 'incorrect':
                err = radio('%s_s%d_err' % (code, i), ERROR_TYPES, p.get('error_type'))
            note = st.text_input('Note (optional)', value=p.get('note', ''),
                                 key='%s_s%d_note' % (code, i))
        step_rows.append({'index': i, 'label': lab, 'error_type': err, 'note': note.strip()})
        st.markdown('---')

    st.markdown('### Milestones')
    st.markdown('<span class="meta">Quantities a correct derivation reaches. Mark whether '
                'this trace obtains each one, in any form, unit or rounding.</span>',
                unsafe_allow_html=True)
    ms_rows = []
    for j, m in enumerate(task['milestones']):
        p = (prior or {}).get('milestones', [{}] * len(task['milestones']))
        p = p[j] if j < len(p) else {}
        c1, c2 = st.columns([1, 2])
        with c1:
            st.markdown('**%s** = %.6g' % (m['id'], m['value']))
        with c2:
            got = radio('%s_m%d' % (code, j), MILESTONE_CHOICES, p.get('status'), horizontal=True)
        ms_rows.append({'id': m['id'], 'value': m['value'], 'status': got})

    st.markdown('### The trace as a whole')
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown('**Final answer**')
        final = radio('%s_final' % code, FINAL_CHOICES, (prior or {}).get('final_answer'))
    with c2:
        st.markdown('**Is the reasoning sound overall?**')
        sound = radio('%s_sound' % code, [('yes', 'Yes'), ('no', 'No')],
                      (prior or {}).get('reasoning_sound'), horizontal=True)
    with c3:
        st.markdown('**Your confidence**')
        conf = radio('%s_conf' % code, CONFIDENCE, (prior or {}).get('confidence'), horizontal=True)
    comment = st.text_area('Comments on this trace (optional)',
                           value=(prior or {}).get('comment', ''), key='%s_comment' % code)
    flag = st.checkbox('Flag for discussion', value=(prior or {}).get('flagged', False),
                       key='%s_flag' % code)

    missing = ([ 'step %d' % (r['index'] + 1) for r in step_rows if not r['label']]
               + ['step %d error type' % (r['index'] + 1)
                  for r in step_rows if r['label'] == 'incorrect' and not r['error_type']]
               + ['milestone %s' % r['id'] for r in ms_rows if not r['status']]
               + (['final answer'] if not final else [])
               + (['reasoning sound'] if not sound else [])
               + (['confidence'] if not conf else []))

    c1, c2 = st.columns([1, 3])
    with c1:
        submit = st.button('Submit and continue', type='primary', disabled=bool(missing))
    with c2:
        if missing:
            st.markdown('<span class="todo">Still to answer: %s</span>'
                        % ', '.join(missing[:8]) + (' ...' if len(missing) > 8 else ''),
                        unsafe_allow_html=True)

    if submit:
        append_label(who, {
            'annotator': who, 'branch_of_annotator': me['branch'], 'code': code,
            'branch': task['branch'], 'level': task['level'],
            'for_truth': code in me['for_truth'], 'calibration': code in me['calibration'],
            'item_sha256': task['item_sha256'], 'trace_sha256': task['trace_sha256'],
            'n_steps': len(task['steps']), 'steps': step_rows, 'milestones': ms_rows,
            'final_answer': final, 'reasoning_sound': sound, 'confidence': conf,
            'comment': comment.strip(), 'flagged': bool(flag), 'source': 'app',
            'seconds': round(time.time() - st.session_state.started, 1),
            'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        })
        nxt = [i for i, c in enumerate(order) if c not in read_labels(who)]
        st.session_state.cursor = nxt[0] if nxt else st.session_state.cursor
        st.rerun()


if __name__ == '__main__':
    main()
