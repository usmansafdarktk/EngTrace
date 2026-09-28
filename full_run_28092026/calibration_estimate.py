"""Re-estimate the full run from the traces recorded so far. FREE: reads traces/, calls nothing.

    python -m full_run_28092026.calibration_estimate [--workers 15] [--providers] [--model KEY ...]

Per model, the billed cost per recorded item (answered or empty) times the items still to run, added
to what is billed already. The 95% interval is a bootstrap over the recorded items, so it assumes they
represent the pool; the calibration takes every k-th item of the sorted pool, one per template. The
bill counts the attempt each row kept: the cost of a retried or failed attempt is not in it. Hours
assume each remaining call takes as long as the recorded ones did, with --workers calls in flight.
The level column reweights the same bills by difficulty level, each level's mean cost times that
level's remaining items, since the recorded items need not mix the levels as the pool does.
"""
from __future__ import annotations

import argparse
import json
import random
from collections import Counter

from full_run_28092026.run_traces import config, existing, items, per_item_cost, trace_path

DRAWS = 10_000
SEED = 0


def pending_failures(key: str, done: dict) -> int:
    """Items whose only rows are service failures: the next run calls them again."""
    path = trace_path(key)
    if not path.exists():
        return 0
    failed = set()
    for ln in path.read_text(encoding='utf-8').splitlines():
        try:
            row = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if row.get('status') == 'service_failure':
            failed.add(row['item_id'])
    return len(failed - set(done))


def interval(draws: list[float]) -> tuple[float, float]:
    d = sorted(draws)
    return d[int(0.025 * (len(d) - 1))], d[int(0.975 * (len(d) - 1))]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--workers', type=int, default=15, help='calls in flight per model')
    ap.add_argument('--providers', action='store_true',
                    help='also list, per model, the endpoints that served the rows and what each billed')
    ap.add_argument('--model', action='append', help='only this model key (repeatable); the total covers these')
    a = ap.parse_args()
    cfg = config()
    if a.model:
        cfg['models'] = [s for s in cfg['models'] if s['key'] in a.model]
        if len(cfg['models']) != len(set(a.model)):
            raise SystemExit(f'unknown model key in {a.model}')
    its = items()
    tok = cfg['pricing_basis_tokens']
    rng = random.Random(SEED)
    pool_levels = Counter(it['level'] for it in its)
    print(f'{"model":22s} {"rows":>4s} {"empty":>5s} {"trunc":>5s} {"fail":>4s} {"out tok":>7s} '
          f'{"billed $":>8s} {"$/item":>8s} {"run $":>7s} {"95% interval":>15s} {"level $":>7s} '
          f'{"assumed $":>9s} {"hours":>5s}')
    totals = [0.0] * DRAWS
    point_total = billed_total = assumed_total = level_total = 0.0
    row_levels = Counter()
    for s in cfg['models']:
        done = existing(s['key'])
        rows = list(done.values())
        if not rows:
            print(f'{s["key"]:22s}    0  no traces')
            continue
        c = Counter(r['status'] for r in rows)
        trunc = sum(r.get('finish_reason') == 'length' for r in rows)
        costs = [r['billed_usd'] for r in rows if r.get('billed_usd') is not None]
        outs = [r['completion_tokens'] for r in rows if r.get('completion_tokens') is not None]
        if not costs:
            print(f'{s["key"]:22s} {len(rows):4d}  no row reports a billed cost')
            continue
        billed = sum(costs)
        per_item = billed / len(costs)
        todo = len(its) - len(rows)
        point = billed + per_item * todo
        draws = [billed + sum(rng.choices(costs, k=len(costs))) / len(costs) * todo
                 for _ in range(DRAWS)]
        totals = [t + d for t, d in zip(totals, draws)]
        lo, hi = interval(draws)
        assumed = len(its) * per_item_cost(s['price_per_m'], tok)
        secs = [r['seconds'] for r in rows if r.get('seconds') is not None]
        hours = (sum(secs) / len(secs)) * todo / a.workers / 3600 if secs else 0.0
        by_level = {}
        for r in rows:
            if r.get('billed_usd') is not None:
                by_level.setdefault(r['level'], []).append(r['billed_usd'])
        run_levels = Counter(r['level'] for r in rows)
        level_point = billed + sum(
            (sum(by_level[lv]) / len(by_level[lv]) if by_level.get(lv) else per_item) * (n - run_levels[lv])
            for lv, n in pool_levels.items())
        row_levels.update(run_levels)
        point_total += point
        billed_total += billed
        assumed_total += assumed
        level_total += level_point
        print(f'{s["key"]:22s} {len(rows):4d} {c["empty"]:5d} {trunc:5d} '
              f'{pending_failures(s["key"], done):4d} {sum(outs) / len(outs) if outs else 0:7.0f} '
              f'{billed:8.3f} {per_item:8.5f} {point:7.2f} {f"{lo:.2f}-{hi:.2f}":>15s} {level_point:7.2f} '
              f'{assumed:9.2f} {hours:5.1f}')
        if len(costs) < len(rows):
            print(f'  {s["key"]}: {len(rows) - len(costs)} rows report no billed cost')
    lo, hi = interval(totals)
    print(f'{"TOTAL":22s} {"":4s} {"":5s} {"":5s} {"":4s} {"":7s} {billed_total:8.3f} {"":8s} '
          f'{point_total:7.2f} {f"{lo:.2f}-{hi:.2f}":>15s} {level_total:7.2f} {assumed_total:9.2f}')
    n_rows = sum(row_levels.values())
    print('\nlevels, recorded rows against the pool: ' + ', '.join(
        f'{lv} {row_levels[lv] / n_rows:.0%} against {n / len(its):.0%}'
        for lv, n in sorted(pool_levels.items())))
    if a.providers:
        print('\nserving endpoints: rows, and billed $ per million completion tokens')
        for s in cfg['models']:
            by = {}
            for r in existing(s['key']).values():
                b = by.setdefault(r.get('provider') or '?', [0, 0.0, 0])
                b[0] += 1
                b[1] += r.get('billed_usd') or 0.0
                b[2] += r.get('completion_tokens') or 0
            if by:
                print(f'  {s["key"]:22s} ' + ', '.join(
                    f'{p} {n} (${usd / toks * 1e6 if toks else 0:.2f}/M)'
                    for p, (n, usd, toks) in sorted(by.items(), key=lambda kv: -kv[1][0])))
    print(f'\nrun $: billed so far plus $/item times the {len(its)}-item pool\'s remaining items. '
          f'assumed $: the dry run\'s basis,\n{tok["in"]} input and {tok["out"]} output tokens per '
          f'item at the pricing document\'s prices. trunc: rows that stopped at the output cap.\n'
          f'fail: items left as service failures, called again by the next run. hours: the '
          f'remaining calls at {a.workers} in flight.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
