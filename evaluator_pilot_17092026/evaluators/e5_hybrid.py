"""E5 - deterministic first, judged only where necessary.

The Suggested Actions: "E3 as primary, non-family LLM adjudication only on residual
unmatched steps. Deterministic where possible, judged where necessary - and the
judged fraction becomes a reportable statistic."

    1. E3 runs as before: which of the item's milestones does the trace state?
    2. For each trace with milestones E3 did NOT find, one judge call asks, per
       missed milestone, one of:
         REACHED       the trace obtains this quantity - possibly in another unit,
                       form or rounding that number matching cannot see
         NOT_NEEDED    the trace reaches the answer by a valid route that does not
                       pass through this quantity
         MISSING       the trace does not obtain it, or obtains it wrongly
    3. e5_strict  = (E3 reached + REACHED) / required
       e5_lenient = (E3 reached + REACHED + NOT_NEEDED) / required
       judged_fraction = milestones sent to the judge / required

WHY TWO SCORES. A judge can excuse any missing quantity as "not needed". Keeping
"the judge found it" apart from "the judge excused it" stops that from inflating the
headline invisibly; the human labels' "correct-but-alternative-path" category is
what decides whether NOT_NEEDED deserves credit.

WHAT THE JUDGE SEES. The question, the trace, and each missed milestone's name and
expected value - not the gold solution. Without the gold text the judge cannot copy
its route, and the same prompt shape can be used for the validity check in
analysis/e5_analysis.py: milestones E3 DID find (should come back REACHED) and the
same milestones with deliberately wrong values (should come back MISSING).

THE JUDGE is MiMo-V2.5-Pro (D-088): a family outside the evaluated suite, the
smallest documented distillation exposure, open weights. Calls use E1's settings and
machinery (D6-D8): JSON mode, temperature 0, 16,384 tokens, empty or malformed
replies re-requested, every reply stored by (model, settings, prompt).
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import e1_panel as e1         # noqa: E402  - fetch(), Replies, settings
import e3_milestones as e3    # noqa: E402
import milestones as ms       # noqa: E402

ID = 'e5'
DESCRIPTION = 'E3 milestones first; MiMo-V2.5-Pro adjudicates only the milestones E3 missed'
JUDGE = 'xiaomi/mimo-v2.5-pro'
REPLIES = os.path.join(e1.e0.CACHE_DIR, 'e5_judge_replies.jsonl')
VERDICTS = ('REACHED', 'NOT_NEEDED', 'MISSING')
DEVIATIONS = []

PROMPT = """You are checking one student's worked solution to an engineering problem.

PROBLEM:
{question}

STUDENT SOLUTION:
{trace}

A correct solution to this problem obtains the intermediate quantities listed below.
An automatic checker could not find them in the student's solution by matching
numbers. For EACH one, decide:

- "REACHED": the student does obtain this quantity correctly - possibly in different
  units, a different but equivalent form, or with different rounding.
- "NOT_NEEDED": the student reaches a correct result by a valid route that does not
  pass through this quantity at all.
- "MISSING": the student does not obtain this quantity, or obtains it incorrectly.

QUANTITIES:
{quantities}

Return a raw JSON object only:
{{"results": [{{"milestone": "<name>", "verdict": "REACHED" | "NOT_NEEDED" | "MISSING"}}]}}
"""


def _sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def config() -> dict:
    return dict(e3.config(), evaluator=ID, judge=JUDGE, judge_call=e1.CALL,
                attempts=e1.ATTEMPTS, prompt_sha256=hashlib.sha256(PROMPT.encode()).hexdigest(),
                e5_sha256=_sha(__file__))


def fmt(v: float) -> str:
    return ('%.6g' % v)


def prompt_for(item: dict, trace_text: str, targets: list[dict]) -> str:
    qs = '\n'.join('- %s = %s' % (m['id'], fmt(m['value'])) for m in targets)
    return PROMPT.format(question=item['question'], trace=trace_text, quantities=qs)


def parse(text: str) -> dict:
    """{milestone id: verdict} from a reply; {} if it cannot be read."""
    t = (text or '').strip()
    if t.startswith('```'):
        t = t.strip('`').split('\n', 1)[-1]
    try:
        data = json.loads(t)
    except Exception:                                              # noqa: BLE001
        return {}
    res = data.get('results') if isinstance(data, dict) else data
    out = {}
    for x in res if isinstance(res, list) else []:
        if isinstance(x, dict):
            v = str(x.get('verdict', '')).upper().replace(' ', '_')
            # The prompt lists "- name = value", and judges echo it back whole
            # ("denom = 1.4248"). Keyed on the full string, 18% of verdicts went
            # unread in validation; the name is what identifies the milestone.
            name = str(x.get('milestone', '')).split('=')[0].strip().strip('`*"\' ')
            out[name] = v if v in VERDICTS else 'MISSING'
    return out


def setup(dry_run: bool = False):
    from dotenv import load_dotenv
    load_dotenv(os.path.join(e1.e0._ROOT, '.env'))
    state = e3.setup(dry_run)
    state['replies'] = e1.Replies(REPLIES)
    state['keys'] = bool(os.environ.get('OPENROUTER_API_KEY'))
    return state


def _missed(state, item, trace):
    mset = state['milestones'][item['item_id']]['milestones']
    hits = e3.reach(mset, trace['text'])
    return hits, [h for h in hits if not h['reached']]


def prefetch(state, jobs, workers: int = 32):
    """Every judge call E5 will make, fetched concurrently before scoring."""
    todo = {}
    for item, trace, seed in jobs:
        _, missed = _missed(state, item, trace)
        if missed:
            p = prompt_for(item, trace['text'], missed)
            k = e1.Replies.key(JUDGE, p)
            if state['replies'].get(k) is None:
                todo[k] = p
    print('  E5 prefetch: %d judge calls to make, %d workers' % (len(todo), workers), flush=True)
    failed = 0
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(e1.fetch, JUDGE, p): k for k, p in todo.items()}
        for i, fut in enumerate(cf.as_completed(futs), 1):
            res = fut.result()
            if res['ok']:
                state['replies'].put(futs[fut], res)
            else:
                failed += 1
            if i % 25 == 0 or i == len(futs):
                print('  E5 prefetch: %d/%d, %d failed' % (i, len(futs), failed), flush=True)
    return {'calls': len(todo), 'failed': failed}


def score(state, item: dict, trace: dict, seed: int) -> dict:
    hits, missed = _missed(state, item, trace)
    n = len(hits)
    calls, verdicts = [], {}
    if missed:
        p = prompt_for(item, trace['text'], missed)
        rep = state['replies'].get(e1.Replies.key(JUDGE, p))
        if rep is None and state['keys']:
            rep = e1.fetch(JUDGE, p)
        rep = rep or {'ok': False, 'error': 'no reply stored and no key'}
        verdicts = parse(rep.get('text')) if rep.get('ok') else {}
        calls.append({'provider': 'e5-judge', 'model': JUDGE, 'ok': rep.get('ok'),
                      'text': rep.get('text'), 'served_model': rep.get('served_model'),
                      'serving_provider': rep.get('serving_provider'),
                      'finish_reason': rep.get('finish_reason'),
                      'prompt_tokens': rep.get('billed_prompt_tokens') or rep.get('prompt_tokens'),
                      'completion_tokens': rep.get('billed_completion_tokens') or rep.get('completion_tokens'),
                      'attempts': rep.get('attempts'), 'error': rep.get('error')})
    for h in hits:
        h['source'] = 'e3' if h['reached'] else verdicts.get(h['id'], 'UNJUDGED')
    count = lambda s: sum(1 for h in hits if h['source'] == s)
    e3n, r, nn = count('e3'), count('REACHED'), count('NOT_NEEDED')
    return {
        'scores': {
            'e5_strict': (e3n + r) / n if n else 0.0,
            'e5_lenient': (e3n + r + nn) / n if n else 0.0,
            'e3_coverage': e3n / n if n else 0.0,
            'judged_fraction': len(missed) / n if n else 0.0,
            'milestones_required': n, 'by_e3': e3n, 'judge_reached': r,
            'judge_not_needed': nn, 'judge_missing': count('MISSING'),
            'unjudged': count('UNJUDGED'),
        },
        'meta': {'milestones': hits, 'judge_failed': bool(missed) and not verdicts,
                 'tribunal_reached_judges': bool(missed)},
        'calls': calls,
        'triggered': bool(missed),
    }
