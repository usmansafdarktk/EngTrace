"""Paid judge calls for the full run's E5 (judge.py) and router (router.py) stages.

One reply store per stage, append-only, in the pilot's format ({"k": key, "v": reply} per line), so a
reply is fetched once and read back on every later run; a thread pool in which each call has a
wall-clock deadline and the whole stage a spend cap; and the billed cost read from each response
(OpenRouter's usage.cost), summed over a reply's attempts. Nothing here decides what is sent or how a
reply is read: the stages do.

THE CAP IS CUMULATIVE (D-148). `run` starts its spend count from what the store already records over
every line, replies and failures alike, so `--max-usd` bounds the stage across resumes, not each
invocation; the first review found the count reset to zero on every run. The cap stops calls from
starting; the calls in flight still finish, so the spend can pass it by about one call per worker.

LATE REPLIES ARE KEPT. A worker writes its own result to the store the moment it returns, so a call
that outlives the deadline is still stored and read back on the next run rather than bought again;
the deadline only stops this run from waiting for it. A key that has failed in MAX_FAILURES runs is
left alone and counted, so a prompt the judge cannot answer is not re-bought every run.

The client makes no retries of its own (`max_retries=0`): the stages' own attempt loops are the
retries, so one fetch cannot outlive the deadline through hidden SDK retries.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import os
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEADLINE = 960          # seconds of wall clock per call: three attempts at a 300 s timeout, plus the sleeps
MAX_FAILURES = 3        # runs in which a key may fail before it is left alone


class Store:
    """key -> reply. Thread-safe, append-only; only replies marked ok are served back.

    Also carries what every line billed (`spent_total`) and how many failed lines each key has
    (`failures`), which the cumulative cap and the per-key limit read."""

    def __init__(self, path: Path):
        self.path, self.lock, self.store = Path(path), threading.Lock(), {}
        self.spent_total, self.failures, self.lines = 0.0, {}, 0
        if self.path.exists():
            for line in self.path.read_text(encoding='utf-8').splitlines():
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                self._note(r.get('k'), r.get('v') or {})

    def _note(self, k, v):
        self.lines += 1
        self.spent_total += v.get('billed_usd') or 0.0
        if v.get('ok'):
            self.store[k] = v
        else:
            self.failures[k] = self.failures.get(k, 0) + 1

    def get(self, k):
        return self.store.get(k)

    def put(self, k, v):
        with self.lock:
            self._note(k, v)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, 'a', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'k': k, 'v': v}, ensure_ascii=False) + '\n')

    def summary(self) -> dict:
        return {'lines': self.lines, 'replies': len(self.store),
                'keys_failed': len(self.failures), 'keys_given_up': sum(n >= MAX_FAILURES for n in self.failures.values()),
                'billed_usd_all_lines': round(self.spent_total, 4)}


def client(timeout: float = 300.0):
    from dotenv import load_dotenv
    from openai import OpenAI
    load_dotenv(REPO / '.env')
    key = os.environ.get('OPENROUTER_API_KEY')
    if not key:
        raise SystemExit('OPENROUTER_API_KEY is not set in .env')
    return OpenAI(api_key=key, base_url='https://openrouter.ai/api/v1', timeout=timeout, max_retries=0)


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


def interleave(per_model: dict[str, dict[str, str]]) -> dict[str, str]:
    """{key: prompt} taking one job from each model in turn, so a stop at the cap leaves every model
    partly judged rather than some models whole and others untouched."""
    out, queues = {}, [list(v.items()) for v in per_model.values()]
    while any(queues):
        for q in queues:
            if q:
                k, p = q.pop(0)
                out[k] = p
    return out


def run(jobs: dict[str, str], fetch, store: Store, workers: int, max_usd: float, label: str) -> dict:
    """Fetch every job {key: prompt} not already in the store, in the order given. Stops starting
    calls once the cumulative spend, what the store records plus this run, passes max_usd; the calls
    already running still finish. A call past its deadline is written as a failure and no longer
    waited for; its reply, if it arrives, is stored by its worker all the same."""
    given_up = {k for k in jobs if store.failures.get(k, 0) >= MAX_FAILURES and store.get(k) is None}
    todo = {k: p for k, p in jobs.items() if store.get(k) is None and k not in given_up}
    spent = start = store.spent_total
    print(f'{label}: {len(jobs)} calls in all, {len(jobs) - len(todo) - len(given_up)} in the store, '
          f'{len(given_up)} given up after {MAX_FAILURES} failed runs, {len(todo)} to make; {workers} workers; '
          f'cap ${max_usd:.2f} of which ${spent:.3f} is already billed', flush=True)
    if spent > max_usd:
        print(f'{label}: the cap is already passed; nothing is called', flush=True)
        todo = {}
    done, failed, abandoned = 0, 0, 0
    started, lock = {}, threading.Lock()

    def worker(k, p):
        with lock:
            started[k] = time.time()
        res = fetch(p)
        store.put(k, res)                       # stored here, so a late reply is kept
        return res

    pool = cf.ThreadPoolExecutor(max_workers=workers)
    futures = {pool.submit(worker, k, p): k for k, p in todo.items()}
    pending, stopped = set(futures), False
    while pending:
        finished, _ = cf.wait(pending, timeout=15, return_when=cf.FIRST_COMPLETED)
        for fut in finished:
            pending.discard(fut)
            if fut.cancelled():
                continue
            res = fut.result()
            spent += res.get('billed_usd') or 0.0
            done += 1
            failed += not res.get('ok')
            if done % 50 == 0 or not pending:
                print(f'  {label}: {done}/{len(todo)} returned, {failed} failed, billed ${spent:.3f} cumulative',
                      flush=True)
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
    out = {'calls': len(todo), 'returned': done, 'failed': failed, 'abandoned': abandoned, 'given_up': len(given_up),
           'billed_usd_this_run': round(spent - start, 4), 'billed_usd_cumulative': round(spent, 4),
           'stopped_at_cap': stopped}
    print(f'{label}: {out}', flush=True)
    pool.shutdown(wait=not abandoned)   # threads stuck on dead sockets would hold the pool
    return out


def finish(summary: dict) -> None:
    """End the process once the stage has written its outputs; threads stuck on dead sockets
    would otherwise keep it alive."""
    if summary.get('abandoned'):
        os._exit(0)


def provenance(store_config: Path | None = None) -> dict:
    """What the stage ran on: the commit, whether the evaluator and stage files were dirty, the tag
    if HEAD carries one, and the digest of the score-store CONFIG the stage read (D-144)."""
    import hashlib
    import subprocess

    def git(*a):
        return subprocess.run(['git', *a], cwd=REPO, capture_output=True, text=True).stdout.strip()
    watched = ['evaluator_pilot_17092026/evaluators', 'evaluation/engineering_parser.py',
               'evaluation/engtrace_evaluation_framework.py', 'full_run_28092026/*.py']
    out = {'git': git('rev-parse', 'HEAD'), 'tag': git('describe', '--tags', '--exact-match') or None,
           'dirty': bool(git('status', '--porcelain', '--', *watched))}
    if store_config is not None and Path(store_config).exists():
        out['store_config_sha256'] = hashlib.sha256(Path(store_config).read_bytes()).hexdigest()
    return out
