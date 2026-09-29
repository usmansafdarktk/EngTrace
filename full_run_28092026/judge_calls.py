"""Paid judge calls for the full run's E5 (judge.py) and router (router.py) stages.

One reply store per stage, append-only, in the pilot's format ({"k": key, "v": reply} per line), so a
reply is fetched once and read back on every later run; a thread pool in which each call has a
wall-clock deadline and the whole run a spend cap; and the billed cost read from each response
(OpenRouter's usage.cost), summed over a reply's attempts. Nothing here decides what is sent or how a
reply is read: the stages do.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import os
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEADLINE = 900          # seconds of wall clock per call, retries included (router_batched.py's rule)


class Store:
    """key -> reply. Thread-safe, append-only; only replies marked ok are served back."""

    def __init__(self, path: Path):
        self.path, self.lock, self.store = Path(path), threading.Lock(), {}
        if self.path.exists():
            for line in self.path.read_text(encoding='utf-8').splitlines():
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get('v', {}).get('ok'):
                    self.store[r['k']] = r['v']

    def get(self, k):
        return self.store.get(k)

    def put(self, k, v):
        with self.lock:
            if v.get('ok'):
                self.store[k] = v
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, 'a', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'k': k, 'v': v}, ensure_ascii=False) + '\n')


def client(timeout: float = 600.0):
    from dotenv import load_dotenv
    from openai import OpenAI
    load_dotenv(REPO / '.env')
    key = os.environ.get('OPENROUTER_API_KEY')
    if not key:
        raise SystemExit('OPENROUTER_API_KEY is not set in .env')
    return OpenAI(api_key=key, base_url='https://openrouter.ai/api/v1', timeout=timeout, max_retries=2)


def one_call(cli, model: str, prompt: str, **kw) -> dict:
    """One completion; the reply's text, usage, billed cost and provider, or the error."""
    t0 = time.time()
    try:
        r = cli.chat.completions.create(model=model, messages=[{'role': 'user', 'content': prompt}], **kw)
        u = r.usage
        return {'ok': True, 'text': r.choices[0].message.content or '', 'served_model': r.model,
                'serving_provider': (r.model_extra or {}).get('provider'),
                'finish_reason': r.choices[0].finish_reason, 'prompt_tokens': u.prompt_tokens,
                'completion_tokens': u.completion_tokens,
                'billed_usd': (getattr(u, 'model_extra', None) or {}).get('cost'),
                'seconds': round(time.time() - t0, 2)}
    except Exception as exc:                                       # noqa: BLE001
        return {'ok': False, 'error': f'{type(exc).__name__}: {str(exc)[:300]}', 'seconds': round(time.time() - t0, 2)}


def run(jobs: dict[str, str], fetch, store: Store, workers: int, max_usd: float, label: str) -> dict:
    """Fetch every job {key: prompt} not already in the store. Stops starting calls once the spend
    passes max_usd; a call past its deadline is written as a failure, so the next run asks it again."""
    todo = {k: p for k, p in jobs.items() if store.get(k) is None}
    print(f'{label}: {len(jobs)} calls in all, {len(todo)} not in the store; {workers} workers, '
          f'cap ${max_usd:.2f}', flush=True)
    spent, done, failed, abandoned = 0.0, 0, 0, 0
    started = {}

    def timed(k, p):
        started[k] = time.time()
        return fetch(p)

    pool = cf.ThreadPoolExecutor(max_workers=workers)
    futures = {pool.submit(timed, k, p): k for k, p in todo.items()}
    pending, stopped = set(futures), False
    while pending:
        finished, _ = cf.wait(pending, timeout=15, return_when=cf.FIRST_COMPLETED)
        for fut in finished:
            pending.discard(fut)
            if fut.cancelled():
                continue
            k = futures[fut]
            res = fut.result()
            store.put(k, res)
            spent += res.get('billed_usd') or 0.0
            done += 1
            failed += not res.get('ok')
            if done % 50 == 0 or not pending:
                print(f'  {label}: {done}/{len(todo)} returned, {failed} failed, billed ${spent:.3f}', flush=True)
        now = time.time()
        for fut in list(pending):
            k = futures[fut]
            if k in started and now - started[k] > DEADLINE and not fut.done():
                pending.discard(fut)
                abandoned += 1
                store.put(k, {'ok': False, 'error': f'no reply within {DEADLINE} s', 'billed_usd': 0.0})
        if spent > max_usd and not stopped:
            stopped = True
            print(f'{label}: the cap ${max_usd:.2f} is passed at ${spent:.3f}; the calls not yet started are '
                  'cancelled', flush=True)
            for fut in list(pending):
                if fut.cancel():
                    pending.discard(fut)
    out = {'calls': len(todo), 'returned': done, 'failed': failed, 'abandoned': abandoned,
           'billed_usd': round(spent, 4), 'stopped_at_cap': stopped}
    print(f'{label}: {out}', flush=True)
    pool.shutdown(wait=not abandoned)   # threads stuck on dead sockets would hold the pool
    return out


def finish(summary: dict) -> None:
    """End the process once the stage has written its outputs; threads stuck on dead sockets
    would otherwise keep it alive."""
    if summary.get('abandoned'):
        os._exit(0)
