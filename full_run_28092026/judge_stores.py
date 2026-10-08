"""E5 on several stores in one process (WS-G, step G2): the judge.py stage, its jobs from every store at once.

    python -m full_run_28092026.judge_stores --variant main --variant paraphrase ... --dry-run     # FREE: the jobs to make
    python -m full_run_28092026.judge_stores --variant main --variant paraphrase ... --yes --max-usd X [--workers 24]

WHY. After a re-score the new E5 jobs are spread over several stores, a few in each. judge.py runs one store per
process, and the stores share one reply store (scores/_judge/e5_replies.jsonl), so their runs must not overlap: a
store's last slow call holds up the next store. Here one process reads the reply store once, takes the jobs of every
store given (judge.jobs_for, the prompts and keys judge.py sends), interleaves them, and sends the ones the store
does not hold through judge_calls.run with judge.fetch: one writer, one cumulative cap (D-148), the same calls,
retries and deadline as judge.py. Then it writes each store's rows with judge.write, as `judge.py --score` does.
Nothing is sent that `judge.py --variant V --yes` would not send for each store in turn, and a prompt two stores
share is bought once.
"""
from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import judge, judge_calls as jc, score  # noqa: E402


def jobs_by_store(variants: list[str], items: dict, ms_all: dict) -> dict[tuple[str, str], dict[str, str]]:
    """{(store, model): {key: prompt}} for every trace E5 sends, as judge.paid builds them per store."""
    return {(v, key): {j[2]: j[1] for _row, _ms, _r, j in judge.jobs_for(v, key, items, ms_all) if j}
            for v in variants for key in judge.roster(v)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--variant', action='append', required=True, help='a store; repeatable')
    ap.add_argument('--dry-run', action='store_true', help='FREE: the jobs per store and model, and how many are new')
    ap.add_argument('--yes', action='store_true', help='required for the run, which bills')
    ap.add_argument('--max-usd', type=float, default=None, help='the cumulative cap over the whole reply store')
    ap.add_argument('--workers', type=int, default=24)
    a = ap.parse_args()
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    store = jc.Store(judge.REPLIES)
    per = jobs_by_store(a.variant, items, ms_all)
    new = collections.Counter()
    for (v, key), js in per.items():
        n = sum(1 for k in js if store.get(k) is None)
        if n:
            new[(v, key)] = n
            print(f'{v:28s} {key:26s} jobs {len(js):5d}, not in the reply store {n}')
    jobs = jc.interleave({f'{v}|{key}': js for (v, key), js in per.items()})
    todo = sum(1 for k in jobs if store.get(k) is None)
    print(f'{len(jobs)} distinct prompts over {len(a.variant)} stores; {todo} not in the reply store '
          f'({sum(new.values())} by store and model); reply store {store.summary()}')
    if a.dry_run:
        print('Nothing was called.')
        return 0
    if not a.yes or a.max_usd is None:
        raise SystemExit('E5 bills: re-run with --yes and --max-usd once the spend is approved')
    cli = jc.client()
    summary = jc.run(jobs, lambda p: judge.fetch(p, cli), store, a.workers, a.max_usd, f'E5 ({judge.e5.JUDGE})')
    for v in a.variant:
        print(f'{v}: {judge.write(v, judge.roster(v), store)}', flush=True)
    jc.finish(summary)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
