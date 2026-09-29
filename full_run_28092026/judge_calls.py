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

A STUCK CALL GIVES BACK ITS WORKER (D-152). The client's 300 s timeout is between bytes, not over the
call, and OpenRouter keeps a waiting connection alive, so on E5's first night about one call in a
hundred never returned and held its thread for good; the fixed pool lost a worker each time. The pool
now has room for STUCK_THREADS abandoned calls beside the `workers` live ones, and a new call starts
whenever a live one returns or passes the deadline. `python -m full_run_28092026.judge_calls --selftest`
checks this offline, with simulated hangs.

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
STUCK_THREADS = 256     # abandoned calls a run may leave on their threads before it stops starting new ones


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

    def reply_stats(self) -> dict:
        """Per reply marked ok: the mean and median completion tokens billed over its attempts (the dry
        runs price the mean), the median seconds of its last attempt, and the mean billed cost where the
        store records one. What the dry runs' pilot basis is checked against (D-152)."""
        import statistics
        rs = list(self.store.values())
        if not rs:
            return {'replies': 0}
        out = [r.get('billed_completion_tokens') or r.get('completion_tokens') or 0 for r in rs]
        billed = [r['billed_usd'] for r in rs if r.get('billed_usd') is not None]
        return {'replies': len(rs), 'completion_tokens_mean': round(statistics.mean(out)),
                'completion_tokens_median': statistics.median(out),
                'seconds_median': round(statistics.median(r.get('seconds') or 0 for r in rs), 1),
                'billed_usd_mean': round(sum(billed) / len(billed), 5) if billed else None}


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
    """Fetch every job {key: prompt} not already in the store, in the order given, `workers` at a time.
    Stops starting calls once the cumulative spend, what the store records plus this run, passes
    max_usd; the calls already running still finish. A call past its deadline is written as a failure
    and no longer waited for, and a new call takes its place; its reply, if it arrives, is stored by its
    worker all the same."""
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

    # `pending` holds the live calls only: an abandoned call keeps its thread, so the pool has room for
    # STUCK_THREADS of them beside the live ones, and top_up() starts a new call in its place (D-152).
    pool = cf.ThreadPoolExecutor(max_workers=workers + STUCK_THREADS)
    queue, futures, pending, stopped = list(todo.items()), {}, set(), False

    def top_up():
        while queue and len(pending) < workers and not stopped and abandoned < STUCK_THREADS:
            k, p = queue.pop(0)
            fut = pool.submit(worker, k, p)
            futures[fut] = k
            pending.add(fut)

    top_up()
    while pending:
        finished, _ = cf.wait(pending, timeout=min(15.0, DEADLINE / 4), return_when=cf.FIRST_COMPLETED)
        for fut in finished:
            pending.discard(fut)
            res = fut.result()
            spent += res.get('billed_usd') or 0.0
            done += 1
            failed += not res.get('ok')
            if done % 50 == 0:
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
            print(f'{label}: the cap ${max_usd:.2f} is passed at ${spent:.3f}; no further calls are started',
                  flush=True)
        top_up()
    if abandoned >= STUCK_THREADS and queue:
        print(f'{label}: {abandoned} calls passed the deadline in this run; the {len(queue)} not started are '
              'left for the next run', flush=True)
    if done % 50 or not done:
        print(f'  {label}: {done}/{len(todo)} returned, {failed} failed, billed ${spent:.3f} cumulative', flush=True)
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


def selftest() -> int:
    """run() offline with a fake fetch and a 1 s deadline: calls that hang give back their workers,
    their late replies are kept, a second run makes only the calls without a reply, the cap stops new
    calls and counts what the store already billed, and STUCK_THREADS bounds the abandoned calls."""
    import io
    import tempfile
    from contextlib import redirect_stdout
    global DEADLINE, STUCK_THREADS
    saved = DEADLINE, STUCK_THREADS
    DEADLINE = 1.0
    tmp = Path(tempfile.mkdtemp())
    jobs = {f'k{i}': f'k{i}' for i in range(40)}
    hang = {'k0', 'k1', 'k2'}          # first in the queue: in a fixed pool of three they held every worker

    def fetcher(release):
        def fetch(p):
            if p in hang:
                release.wait(30)
                return {'ok': True, 'text': 'late', 'billed_usd': 0.01}
            time.sleep(0.02)
            return {'ok': True, 'text': p, 'billed_usd': 0.01}
        return fetch

    def quiet(*a):
        with redirect_stdout(io.StringIO()):
            return run(*a)
    checks = {}
    try:
        release = threading.Event()
        store = Store(tmp / 'replies.jsonl')
        t0 = time.time()
        out = quiet(jobs, fetcher(release), store, 3, 100.0, 'test')
        checks['three hung calls on three workers: the run still ends in seconds'] = (
            out['abandoned'] == 3 and time.time() - t0 < 10)
        checks['every other call is answered'] = out['returned'] == 37 and all(
            store.get(f'k{i}') for i in range(3, 40))
        checks['a hung call is written as one failure'] = all(store.failures.get(k) == 1 for k in hang)
        release.set()
        time.sleep(0.5)
        checks['a late reply is stored'] = all(store.get(k) for k in hang)
        out = quiet(jobs, fetcher(release), Store(tmp / 'replies.jsonl'), 3, 100.0, 'test')
        checks['a second run makes no call already answered'] = out['calls'] == 0
        out = quiet(jobs, fetcher(release), Store(tmp / 'capped.jsonl'), 3, 0.105, 'test')
        checks['the cap stops new calls: 11 to 14 of 40 made at $0.01 under a $0.105 cap'] = (
            out['stopped_at_cap'] and 11 <= out['returned'] <= 14)
        out = quiet(jobs, fetcher(release), Store(tmp / 'capped.jsonl'), 3, 0.105, 'test')
        checks['the cap counts what the store already billed'] = out['calls'] == 0
        STUCK_THREADS = 2
        stuck = threading.Event()
        out = quiet(jobs, fetcher(stuck), Store(tmp / 'stuck.jsonl'), 3, 100.0, 'test')
        checks['past STUCK_THREADS abandoned calls, no new call starts'] = (
            out['abandoned'] == 3 and out['returned'] == 0)
        stuck.set()
    finally:
        DEADLINE, STUCK_THREADS = saved
    for name, ok in checks.items():
        print(f'  {"ok  " if ok else "FAIL"} {name}')
    ok = all(checks.values())
    print('selftest: all pass' if ok else 'selftest: FAILED')
    return 0 if ok else 1


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


if __name__ == '__main__':
    import sys
    if '--selftest' not in sys.argv:
        raise SystemExit('usage: python -m full_run_28092026.judge_calls --selftest   # FREE, offline')
    raise SystemExit(selftest())
