"""One harness for every evaluator candidate: (evaluator, traces) -> scores.

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY -m evaluator_pilot_17092026.run_evaluator e0 --dry-run    # free: Tier 1 for real, judges recorded not called
    $PY -m evaluator_pilot_17092026.run_evaluator e0 --smoke 3    # a few real entries, all judges, pennies
    $PY -m evaluator_pilot_17092026.run_evaluator e0              # the run, resumable
    $PY -m evaluator_pilot_17092026.run_evaluator e0 --status

WHY ONE HARNESS.  The comparison between candidates is only fair if they are fed
the same traces in the same way and their outputs are recorded the same way. So
the evaluator modules under evaluators/ only answer "score this trace"; selecting
traces, seeding, resuming, costing and writing live here, once.

AN EVALUATOR MODULE provides ID, DESCRIPTION, config() -> dict (everything that
changes a score), setup(dry_run) -> state, and score(state, item, trace, seed).

WHAT THE HARNESS GUARANTEES.

  1.  Scores attach to verified traces.  The freeze must verify and
      verify_traces.py must pass T1-T7 before anything is scored.  A score on a
      truncated or mis-anchored trace would be carried into the comparison with
      nothing to mark it.

  2.  Deterministic.  Each (evaluator, item, model) gets a seed from BLAKE2b,
      recorded on the row, so any stochastic part of an evaluator is repeatable.

  3.  Resumable and change-aware.  A row is reused only if the trace text hash
      AND the evaluator config hash both match.  Editing the evaluator, or its
      framework file, re-scores; an interrupted run re-costs nothing.

  4.  Every judge call is on the row - raw text, finish reason, tokens, cost - so
      a judge that failed is never indistinguishable from a judge that said no.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import sys
import time
from collections import Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_traces as rt        # noqa: E402
import verify_traces as vt     # noqa: E402

SCORES = os.path.join(_HERE, 'scores')
EVALUATORS = {'e0': 'evaluators.e0_tribunal'}

# Snapshot prices for routes without a live catalogue; OpenRouter is read live.
GOOGLE_PRICES = {'gemini-3.1-pro-preview': (2.0, 12.0)}


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def seed_for(evaluator: str, item_id: str, model_key: str) -> int:
    h = hashlib.blake2b(('%s|%s|%s' % (evaluator, item_id, model_key)).encode(), digest_size=4)
    return int.from_bytes(h.digest(), 'big')


def config_sha(cfg: dict) -> str:
    return _sha(json.dumps(cfg, sort_keys=True))


def verified_traces():
    """The traces verify_traces passes, or stop."""
    import freeze as fz
    if fz.verify(fz.DEFAULT_MASTER_SEED) != 0:
        raise SystemExit('the frozen slice does not verify - refusing to score')
    items, traces = vt.load()
    bad = vt.check(items, traces)
    if bad:
        raise SystemExit('verify_traces fails - refusing to score:\n  ' + '\n  '.join(bad))
    return items, {k: {r['item_id']: r for r in rows if r['ok']} for k, rows in traces.items()}


def out_path(evaluator: str, model_key: str) -> str:
    return os.path.join(SCORES, evaluator, model_key + '.jsonl')


def done_rows(evaluator: str, model_key: str, csha: str) -> dict:
    path = out_path(evaluator, model_key)
    if not os.path.exists(path):
        return {}
    out = {}
    for ln in open(path, encoding='utf-8'):
        try:
            r = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if r.get('config_sha256') == csha and not r.get('error'):
            out[(r['item_id'], r['trace_sha256'])] = r
    return out


def call_cost(call: dict, prices: dict) -> float:
    if not call.get('ok'):
        return 0.0
    if call['provider'] == 'google':
        pin, pout = GOOGLE_PRICES.get(call['model'], (0.0, 0.0))
        out = (call.get('completion_tokens') or 0) + (call.get('thinking_tokens') or 0)
    else:
        pin, pout = prices.get(call['model'], (0.0, 0.0))
        out = call.get('completion_tokens') or 0
    return ((call.get('prompt_tokens') or 0) * pin + out * pout) / 1e6


def plan(items, traces, evaluator, csha, limit=None, models=None):
    work = []
    for key in sorted(traces):
        if models and key not in models:
            continue
        done = done_rows(evaluator, key, csha)
        for item_id in sorted(items):
            tr = traces[key].get(item_id)
            if tr is None:
                continue
            tsha = _sha(tr['text'])
            if (item_id, tsha) not in done:
                work.append((key, item_id, tr, tsha))
    if limit:
        # A smoke test that never reaches the judges proves nothing about them.
        # If a dry run exists, take entries it says reach the Tribunal, one per
        # model where possible, before anything else.
        reach = dry_reached(evaluator)
        hit = [w for w in work if (w[0], w[1]) in reach]
        spread, seen = [], set()
        for w in hit:
            if w[0] not in seen:
                spread.append(w)
                seen.add(w[0])
        rest = [w for w in hit if w not in spread] + [w for w in work if w not in hit]
        work = (spread + rest)[:limit]
    return work


def dry_reached(evaluator: str) -> set:
    d = os.path.join(SCORES, evaluator + '_dry')
    if not os.path.isdir(d):
        return set()
    out = set()
    for fn in os.listdir(d):
        for ln in open(os.path.join(d, fn), encoding='utf-8'):
            r = json.loads(ln)
            if not r.get('error') and r['meta'].get('tribunal_reached_judges'):
                out.add((r['model_key'], r['item_id']))
    return out


def run(evaluator: str, dry_run: bool, limit, models):
    mod = importlib.import_module(EVALUATORS[evaluator])
    items, traces = verified_traces()
    cfg = mod.config()
    csha = config_sha(cfg)
    work = plan(items, traces, evaluator if not dry_run else evaluator + '_dry', csha, limit, models)
    print('%s: %s' % (mod.ID, mod.DESCRIPTION))
    print('config %s   %d (item, model) pairs to score%s'
          % (csha[:12], len(work), '  [DRY RUN - no judge is called]' if dry_run else ''))
    for d in getattr(mod, 'DEVIATIONS', []):
        print('  deviation: ' + d)
    if not work:
        return 0

    print('loading evaluator (first load reads the scorer models from the HF cache)...')
    t0 = time.time()
    state = mod.setup(dry_run=dry_run)
    print('  ready in %.0fs' % (time.time() - t0))
    prices = {} if dry_run else rt.live_prices()

    tag = evaluator + ('_dry' if dry_run else '')
    os.makedirs(os.path.join(SCORES, tag), exist_ok=True)
    totals = Counter()
    spent = 0.0
    for n, (key, item_id, tr, tsha) in enumerate(work, 1):
        seed = seed_for(evaluator, item_id, key)
        row = {'evaluator': evaluator, 'item_id': item_id, 'model_key': key,
               'template_id': items[item_id]['template_id'], 'level': items[item_id]['level'],
               'branch': items[item_id]['branch'], 'answer_type': items[item_id]['answer_type'],
               'item_sha256': items[item_id]['sha256'], 'trace_sha256': tsha,
               'config_sha256': csha, 'seed': seed, 'dry_run': dry_run,
               'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        t1 = time.time()
        try:
            res = mod.score(state, items[item_id], tr, seed)
            for c in res.get('calls', []):
                c['cost_usd'] = round(call_cost(c, prices), 6)
                spent += c['cost_usd']
            row.update(res, seconds=round(time.time() - t1, 1))
            totals['scored'] += 1
            totals['triggered'] += bool(res.get('triggered'))
            totals['reached_judges'] += bool(res['meta'].get('tribunal_reached_judges'))
            totals['judge_failures'] += len(res['meta'].get('judge_failures', []))
            totals['judge_truncations'] += len(res['meta'].get('judges_truncated', []))
            totals['scorer_failures'] += bool(res['meta'].get('scorer_failures'))
        except Exception as exc:                                   # noqa: BLE001
            row['error'] = '%s: %s' % (type(exc).__name__, str(exc)[:500])
            totals['errors'] += 1
        with open(os.path.join(SCORES, tag, key + '.jsonl'), 'a', encoding='utf-8', newline='\n') as fh:
            fh.write(json.dumps(row, ensure_ascii=False) + '\n')
        if n % 10 == 0 or n == len(work):
            print('  %3d/%-3d  scored %3d  to judges %3d  judge failures %d  truncated %d  '
                  'scorer failures %d  errors %d  $%.3f'
                  % (n, len(work), totals['scored'], totals['reached_judges'],
                     totals['judge_failures'], totals['judge_truncations'],
                     totals['scorer_failures'], totals['errors'], spent))

    if dry_run:
        dry_estimate(tag)
    return 1 if totals['errors'] or totals['scorer_failures'] else 0


def dry_estimate(tag: str):
    """From the dry run's recorded prompts, what the real Tribunal would cost."""
    prices = rt.live_prices()
    rows = []
    for fn in os.listdir(os.path.join(SCORES, tag)):
        rows += [json.loads(ln) for ln in open(os.path.join(SCORES, tag, fn), encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    reach = [r for r in rows if r['meta'].get('tribunal_reached_judges')]
    chars = [p['prompt_chars'] for r in reach for p in r.get('dry_prompts', [])]
    in_tok = sum(chars) / 3.5           # JSON-escaped LaTeX runs short of 4 chars/token
    judge_out = 400                     # a results list; GPT-5 also bills reasoning, see below
    cost = 0.0
    for mid in ('openai/gpt-5', 'anthropic/claude-opus-4.5'):
        pin, pout = prices.get(mid, (0, 0))
        cost += (in_tok * pin + len(chars) * judge_out * pout) / 1e6
    pin, pout = GOOGLE_PRICES['gemini-3.1-pro-preview']
    cost += (in_tok * pin + len(chars) * judge_out * pout) / 1e6
    print('\nDRY RUN: %d of %d traces would reach the Tribunal (%d triggered)'
          % (len(reach), len(rows), sum(1 for r in rows if r.get('triggered'))))
    print('  by model  %s' % dict(Counter(r['model_key'] for r in reach)))
    print('  by level  %s' % dict(Counter(r['level'] for r in reach)))
    print('  %d Tribunal prompts, ~%.0f input tokens each, x3 judges' % (len(chars), in_tok / max(len(chars), 1)))
    print('  estimate ~$%.2f before reasoning tokens. GPT-5 and Gemini think before answering;'
          ' at stage 1 that multiplied cost ~5x, so budget up to ~$%.2f.' % (cost, cost * 5))


def status(evaluator: str):
    mod = importlib.import_module(EVALUATORS[evaluator])
    csha = config_sha(mod.config())
    d = os.path.join(SCORES, evaluator)
    print('%s config %s' % (evaluator, csha[:12]))
    if not os.path.isdir(d):
        print('  nothing scored yet')
        return 0
    grand = 0.0
    for fn in sorted(os.listdir(d)):
        rows = [json.loads(ln) for ln in open(os.path.join(d, fn), encoding='utf-8')]
        cur = [r for r in rows if r.get('config_sha256') == csha and not r.get('error')]
        cost = sum(c.get('cost_usd', 0) for r in cur for c in r.get('calls', []))
        grand += cost
        fails = sum(len(r['meta'].get('judge_failures', [])) for r in cur)
        trunc = sum(len(r['meta'].get('judges_truncated', [])) for r in cur)
        print('  %-16s %3d scored  %3d reached judges  failures %d  truncated %d  $%.3f'
              % (fn[:-6], len(cur), sum(1 for r in cur if r['meta'].get('tribunal_reached_judges')),
                 fails, trunc, cost))
    print('  total $%.3f' % grand)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('evaluator', choices=sorted(EVALUATORS))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--smoke', type=int, metavar='N', help='score N entries for real, then stop')
    ap.add_argument('--model', action='append')
    ap.add_argument('--status', action='store_true')
    args = ap.parse_args()
    if args.status:
        return status(args.evaluator)
    return run(args.evaluator, args.dry_run, args.smoke, args.model)


if __name__ == '__main__':
    sys.exit(main())
