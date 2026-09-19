"""E1 - E0's Tribunal, judged by a panel from families outside the evaluated suite.

The Suggested Actions: "same architecture, judges drawn only from families not in
the evaluation suite; inter-judge agreement reported". So E1 changes WHO judges and
nothing about HOW the Tribunal works: the framework's own prompt, its own JSON
parsing, its own majority vote and its own post-judgement recovery all run
unmodified, imported from evaluation/engtrace_evaluation_framework.py exactly as E0
does. The panel (D-088, JUDGE_SELECTION.md):

    OpenAI slot     x-ai/grok-4.6          xAI       closed
    Anthropic slot  minimax/minimax-m3     MiniMax   open weights
    Google slot     xiaomi/mimo-v2.5-pro   Xiaomi    open weights

Compare E1 against E0-3J, not E0. E0-3J is the same framework with its three
original judges all connected; E0 as published silently runs two (E0-F6). E1 against
E0-3J differs in the judges' families and nothing else.

DEVIATIONS, on every row. D2 and D3 as in E0 (OpenRouter transport; the wrong-answer
sample seeded from (item, model), shared with every judged evaluator). D4 as in
E0-3J (get_model_info supplied, so the third slot is really called). And:

  D6  Uniform judge call settings. The framework gives each provider slot its own
      settings, tuned to E0's original judges: the OpenAI slot forces JSON mode and
      temperature 0 for any non-GPT-5 model; the Anthropic slot caps output at 2,048
      tokens (set for Claude, which answered in ~300); the Google slot asks for a
      JSON MIME type. Kept per slot, a judge's behaviour would depend on which slot
      it sat in, and the 2,048 cap would truncate reasoning judges - in the probe
      MiniMax exceeded it on 4 of 21 replies, Grok on 9, MiMo on 16. A cut-off reply
      is unparseable, and the framework silently drops a judge it cannot parse,
      which would recreate E0-F6's two-judge panel. Every judge gets the same
      settings: JSON mode, temperature 0, a 16,384-token ceiling. (8,192 was the
      first choice; the smoke test's real multi-step prompts drew 6,849 and 7,070
      tokens from Grok and MiMo, too close to it.)

  D7  An empty reply is re-requested, up to three attempts. In the smoke test MiMo
      returned no content at all - no finish reason, 1,598 tokens billed, all of it
      reasoning. The framework would drop that judge silently. Every attempt is
      recorded on the row.

HOW IT RUNS FAST WITHOUT CHANGING WHAT IT COMPUTES. The framework calls its three
judges one after another, per trace, so a run waits for 3 x 178 slow replies in
sequence - about ten hours. And the laptop has room for only one process that loads
the scorer models. So `prefetch` runs the framework's own path once with the judges
stubbed, capturing every exact prompt it sends; fans all those calls out
concurrently (HTTP only, negligible memory); and stores each reply keyed by
(model, settings, prompt). The real run then calls the judges exactly as before and
receives each reply from that store in the framework's own order. Any prompt not in
the store is called live, so correctness never depends on the capture pass - only
the speed does.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import io
import json
import os
import random
import sys
import threading
import time
from contextlib import redirect_stdout
from types import SimpleNamespace

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import e0_tribunal as e0  # noqa: E402

ID = 'e1'
DESCRIPTION = 'E0 Tribunal with a non-suite panel: Grok 4.6, MiniMax M3, MiMo-V2.5-Pro'

PANEL = {
    'openai':    {'model': 'x-ai/grok-4.6',        'family': 'xAI'},
    'anthropic': {'model': 'minimax/minimax-m3',   'family': 'MiniMax'},
    'google':    {'model': 'xiaomi/mimo-v2.5-pro', 'family': 'Xiaomi'},
}
CALL = {'response_format': {'type': 'json_object'}, 'temperature': 0.0, 'max_tokens': 16384}
ATTEMPTS = 3
# D8. OpenRouter serves MiniMax M3 from 13 providers. In the first full fetch, 11 of
# the 17 replies served by ModelRun were malformed JSON ('{"results":":[{",...}',
# '{".results":[]}'); the other providers served 161 with none malformed.
PROVIDER_IGNORE = {'minimax/minimax-m3': ['ModelRun']}
REPLIES = os.path.join(e0.CACHE_DIR, 'e1_judge_replies.jsonl')

DEVIATIONS = [
    'D2 all three judges reached through OpenRouter',
    'D3 random seeded per (item, model), shared with every judged evaluator; the 0.20 rate is unchanged',
    'D4 genai.get_model_info supplied so the third slot is called (E0-F6)',
    'D6 uniform judge call settings in every slot: JSON mode, temperature 0, max_tokens 16384 '
    '(replacing per-slot settings tuned to E0\'s judges; the 2,048 cap would truncate reasoning judges)',
    'D7 an empty or malformed judge reply (valid JSON without a "results" list of objects) is '
    're-requested, up to 3 attempts, all recorded; a malformed one otherwise crashes the framework\'s whole trace',
    'D8 provider ModelRun excluded for MiniMax M3: 11 of its 17 replies were malformed JSON, '
    'against 0 of 161 from MiniMax\'s other providers; its 17 replies were re-fetched',
]

libraries = e0.libraries
compute = e0.compute
score = e0.score


def config() -> dict:
    cfg = e0.config()
    cfg.update(evaluator=ID, judges={k: v['model'] for k, v in PANEL.items()},
               judge_call=CALL, attempts=ATTEMPTS, provider_ignore=PROVIDER_IGNORE,
               deviations=DEVIATIONS)
    return cfg


def malformed(text: str) -> bool:
    """Valid JSON whose "results" is not a list of objects.

    The framework iterates `results` calling .get() on each item, so a string or a
    list of strings raises AttributeError and takes the whole trace down with it -
    not just that judge's vote. Text that is not JSON at all is left alone: the
    framework's own anchor-search parser may still recover it, and if not it drops
    the judge the way it always has.
    """
    t = (text or '').strip()
    if t.startswith('```'):
        t = t.strip('`').split('\n', 1)[-1] if '\n' in t else t
    try:
        data = json.loads(t)
    except Exception:                                              # noqa: BLE001
        return False
    res = data.get('results', data.get('steps')) if isinstance(data, dict) else data
    if isinstance(res, dict):
        res = [res]
    return not isinstance(res, list) or not res or any(not isinstance(x, dict) for x in res)


# ------------------------------------------------------------ reply store

class Replies:
    """(model, settings, prompt) -> the judge's reply. Thread-safe, append-only."""

    def __init__(self, path=REPLIES):
        self.path, self.lock, self.store = path, threading.Lock(), {}
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path):
            for ln in open(path, encoding='utf-8'):
                try:
                    r = json.loads(ln)
                    if r['v'].get('ok'):
                        self.store[r['k']] = r['v']
                except (json.JSONDecodeError, KeyError):
                    continue

    @staticmethod
    def key(model, prompt):
        return hashlib.sha256(json.dumps([model, CALL, prompt], sort_keys=True).encode()).hexdigest()

    def get(self, k):
        return self.store.get(k)

    def put(self, k, v):
        with self.lock:
            self.store[k] = v
            with open(self.path, 'a', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'k': k, 'v': v}, ensure_ascii=False) + '\n')


def fetch(model: str, prompt: str) -> dict:
    """One judge call with E1's uniform settings; an empty reply is asked again (D7)."""
    from openai import OpenAI
    cli = OpenAI(api_key=os.environ['OPENROUTER_API_KEY'],
                 base_url='https://openrouter.ai/api/v1', timeout=600.0, max_retries=2)
    tries = []
    for attempt in range(1, ATTEMPTS + 1):
        t0 = time.time()
        try:
            extra = {'provider': {'ignore': PROVIDER_IGNORE[model]}} if model in PROVIDER_IGNORE else None
            r = cli.chat.completions.create(model=model, messages=[{'role': 'user', 'content': prompt}],
                                            extra_body=extra, **CALL)
            ch = r.choices[0]
            out = {'ok': True, 'text': ch.message.content or '', 'served_model': r.model,
                   'serving_provider': (r.model_extra or {}).get('provider'),
                   'finish_reason': ch.finish_reason, 'prompt_tokens': r.usage.prompt_tokens,
                   'completion_tokens': r.usage.completion_tokens,
                   'seconds': round(time.time() - t0, 2)}
        except Exception as exc:                                   # noqa: BLE001
            out = {'ok': False, 'error': '%s: %s' % (type(exc).__name__, str(exc)[:300]),
                   'seconds': round(time.time() - t0, 2)}
        if out['ok'] and malformed(out['text']):
            out['malformed'] = True
        tries.append({k: out.get(k) for k in ('ok', 'finish_reason', 'completion_tokens', 'seconds',
                                              'error', 'serving_provider', 'malformed')})
        if out['ok'] and out['text'].strip() and not out.get('malformed'):
            break
    out['attempts'] = tries
    # Tokens of every attempt were billed, not just the last one's.
    out['billed_completion_tokens'] = sum(t.get('completion_tokens') or 0 for t in tries)
    out['billed_prompt_tokens'] = (out.get('prompt_tokens') or 0) * sum(1 for t in tries if t.get('ok'))
    if out['ok'] and not out['text'].strip():
        out['ok'] = False
        out['error'] = 'empty reply after %d attempts' % ATTEMPTS
    elif out['ok'] and out.get('malformed'):
        out['ok'] = False
        out['error'] = 'malformed reply after %d attempts' % ATTEMPTS
    return out


def refetch(captured: list[dict], providers: set, workers: int = 16) -> dict:
    """Re-fetch every stored reply that was served by one of `providers` (D8)."""
    replies = Replies()
    todo = [c for c in captured
            if (replies.get(c['key']) or {}).get('serving_provider') in providers]
    print('  re-fetching %d replies served by %s' % (len(todo), sorted(providers)), flush=True)
    got = failed = 0
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for c, res in zip(todo, pool.map(lambda c: fetch(c['model'], c['prompt']), todo)):
            if res['ok']:
                replies.put(c['key'], res)
                got += 1
            else:
                failed += 1
                print('    FAILED %s: %s' % (c['model'], res.get('error')), flush=True)
    return {'refetched': got, 'failed': failed}


# ------------------------------------------------------------ the clients

class _Capture(Exception):
    pass


def _judge(log, slot, replies, state):
    def create(model, messages, **_framework_settings):
        prompt = messages[-1]['content']
        k = Replies.key(model, prompt)
        if state.capturing is not None:
            state.capturing.append((slot, model, prompt, k))
            raise _Capture()                     # the framework swallows it; nothing is scored
        hit = replies.get(k)
        prefetched = hit is not None
        if hit is None and state.keys:
            hit = fetch(model, prompt)           # live fallback: correctness never needs the prefetch
            if hit['ok']:
                replies.put(k, hit)
        elif hit is None:
            # No key here (a Kaggle replay). A missing reply is an error on the row,
            # never a silent judge drop - that silence is exactly E0-F6.
            hit = {'ok': False, 'error': 'reply not in the store and no key to fetch it: '
                                         'the capture and this run disagree on a prompt'}
        rec = {'provider': slot, 'model': model, 'max_tokens': CALL['max_tokens'],
               'prefetched': prefetched,
               'prompt_tokens': hit.get('billed_prompt_tokens') or hit.get('prompt_tokens'),
               'completion_tokens': hit.get('billed_completion_tokens') or hit.get('completion_tokens')}
        rec.update({k2: hit.get(k2) for k2 in ('ok', 'text', 'served_model', 'serving_provider',
                                               'finish_reason', 'seconds', 'attempts', 'error')})
        log.calls.append(rec)
        if not hit['ok']:
            raise RuntimeError(hit.get('error'))
        return SimpleNamespace(
            model=hit['served_model'], model_extra={'provider': hit.get('serving_provider')},
            choices=[SimpleNamespace(message=SimpleNamespace(content=hit['text']),
                                     finish_reason=hit['finish_reason'])],
            usage=SimpleNamespace(prompt_tokens=hit.get('prompt_tokens'),
                                  completion_tokens=hit.get('completion_tokens')))
    return create


def setup(dry_run: bool = False):
    if dry_run:
        # E1's trigger is E0's trigger, so its dry run is E0's; nothing new to learn.
        return e0.setup(dry_run=True)
    from dotenv import load_dotenv
    load_dotenv(os.path.join(e0._ROOT, '.env'))
    keys = bool(os.environ.get('OPENROUTER_API_KEY'))
    state = e0.setup(dry_run=False, keys=False)      # E1 injects its own clients below
    fw, framework, log, genai = state.fw, state.framework, state.log, state.genai
    state.capturing = None
    state.keys = keys
    state.replies = Replies()

    fw.MODEL_OPENAI = PANEL['openai']['model']
    fw.MODEL_ANTHROPIC = PANEL['anthropic']['model']
    fw.MODEL_GOOGLE = PANEL['google']['model']

    framework.client_openai = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(
        create=_judge(log, 'openai', state.replies, state))))

    anth = _judge(log, 'anthropic', state.replies, state)

    def anthropic_create(model, max_tokens, messages):
        r = anth(model=model, messages=messages)
        return SimpleNamespace(content=[SimpleNamespace(text=r.choices[0].message.content or '')])

    framework.client_anthropic = SimpleNamespace(messages=SimpleNamespace(create=anthropic_create))

    goog = _judge(log, 'google', state.replies, state)

    class PanelModel:
        def __init__(self, name, *a, **kw):
            self._name = name

        def generate_content(self, prompt, **kw):
            r = goog(model=self._name, messages=[{'role': 'user', 'content': prompt}])
            return SimpleNamespace(text=r.choices[0].message.content or '')

    genai.GenerativeModel = PanelModel
    genai.get_model_info = lambda name: name          # D4: the slot is really called
    return state


# ------------------------------------------------------------ prefetch

def capture(state, jobs) -> list[dict]:
    """Every judge prompt the real run will send, found by running the framework.

    `jobs` are (item, trace, seed). `random` is seeded and the framework called
    exactly as `score` does, so the same traces trigger the Tribunal and build the
    same prompts; the judges are stubbed and nothing is scored or written. Needs no
    key, which is what lets it run on a Kaggle kernel.
    """
    state.capturing = []
    buf = io.StringIO()
    for item, trace, seed in jobs:
        random.seed(seed)
        with redirect_stdout(buf):
            state.framework.evaluate_entry({'question': item['question'],
                                            'solution': item['solution'],
                                            'generation': trace['text']})
    captured, state.capturing = state.capturing, None
    return [{'slot': s, 'model': m, 'prompt': p, 'key': k} for s, m, p, k in captured]


def fetch_captured(captured: list[dict], workers: int = 48) -> dict:
    """Fetch every captured prompt's reply concurrently into the reply store.

    Plain HTTP: it imports neither torch nor the framework, so it runs on the laptop
    with almost no memory, and it is the only step that holds the key.
    """
    replies = Replies()
    todo = {c['key']: c for c in captured if replies.get(c['key']) is None}
    print('  %d captured prompts, %d already stored, %d to fetch with %d workers'
          % (len(captured), len(captured) - len(todo), len(todo), workers), flush=True)
    done = failed = 0
    t1 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(fetch, c['model'], c['prompt']): k for k, c in todo.items()}
        for fut in cf.as_completed(futs):
            res = fut.result()
            if res['ok']:
                replies.put(futs[fut], res)
                done += 1
            else:
                failed += 1
                print('    FAILED %s: %s; attempts %s' % (
                    todo[futs[fut]]['model'], res.get('error'),
                    [(a.get('finish_reason'), a.get('completion_tokens')) for a in res.get('attempts', [])]),
                    flush=True)
            n = done + failed
            if n % 50 == 0 or n == len(futs):
                print('  %d/%d fetched, %d failed, %.0fs' % (n, len(futs), failed, time.time() - t1), flush=True)
    return {'captured': len(captured), 'fetched': done, 'failed': failed,
            'stored_total': len(replies.store)}


def prefetch(state, jobs, workers: int = 48):
    """Local, single-machine version: capture, then fetch, in one process."""
    caps = capture(state, jobs)
    return fetch_captured(caps, workers)
