"""Run the five pilot models over the frozen slice, and record what was run.

    python -m evaluator_pilot_17092026.run_traces --check      # preflight: keys, ids, prices, 1 tiny call each
    python -m evaluator_pilot_17092026.run_traces --dry-run    # the plan and the estimated spend, no calls
    python -m evaluator_pilot_17092026.run_traces              # the real run, resumable
    python -m evaluator_pilot_17092026.run_traces --status     # what exists so far, and what it cost
    python -m evaluator_pilot_17092026.run_traces --model gpt-5 --limit 5     # a slice of the slice

WHAT THIS PRODUCES.  One JSONL per model under traces/, one row per
item, carrying the raw output AND the things that make it re-checkable later:
the item's SHA-256 from the manifest, the prompt hash, the model id as
configured, the model id the API says it served, token usage, computed cost and
the attempt count.

THREE RULES IT ENFORCES, each because the record shows what happens otherwise.

  1.  The FREEZE must verify before a single call is made.  Traces are the thing
      experts annotate; if the items moved since the freeze, the labels would be
      attached to text nobody read.  --check and the real run both refuse to
      start on a failed verify.

  2.  The PROMPT must be the deployed one.  `evaluation/run_inference.py` is what
      produced the published generations; this module copies its template
      verbatim and hashes it at startup against that file.  If someone edits
      either, the run stops rather than quietly producing traces that are not
      comparable with the archive.

  3.  RESUME, never restart.  `run_inference.py` opens its output with "w", so a
      rate-limit stall loses everything before it.  Here each (item, model) is
      written as it completes and skipped if already present, so an interrupted
      run costs nothing but the calls it had not yet made.

COST.  Printed per model and in total, from the live OpenRouter catalogue for
routed models and from models.json's snapshot for the native ones.  The run also
prints what it actually spent, computed from returned token usage, so the
estimate can be checked against the bill.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import re
import sys
import time
import urllib.request
from collections import Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from dotenv import load_dotenv  # noqa: E402

load_dotenv(os.path.join(_ROOT, '.env'))

SLICE = os.path.join(_HERE, 'slice')
MANIFEST = os.path.join(SLICE, 'manifest.jsonl')
TRACES = os.path.join(_HERE, 'traces')
CONFIG = os.path.join(_HERE, 'models.json')
DEPLOYED_RUNNER = os.path.join(_ROOT, 'evaluation', 'run_inference.py')

MAX_ATTEMPTS = 4
BACKOFF = 5          # seconds, doubled per attempt
MAX_TOKENS = 4096       # per-model override: "max_tokens" in models.json

# OpenAI's reasoning models reject max_tokens; the other two endpoints want it.
TOKEN_PARAM = {'openai': 'max_completion_tokens'}


# --------------------------------------------------------------- the prompt

def deployed_prompt() -> str:
    """The template evaluation/run_inference.py uses, read from that file.

    Copying it would let the two drift; reading it means a pilot trace is
    produced by the same instructions as an archive trace, or the run fails.
    """
    src = open(DEPLOYED_RUNNER, encoding='utf-8').read()
    m = re.search(r'PROMPT_TEMPLATE = """(.*?)"""', src, re.S)
    if not m:
        raise SystemExit('cannot find PROMPT_TEMPLATE in %s' % DEPLOYED_RUNNER)
    return m.group(1)


PROMPT = deployed_prompt()
PROMPT_SHA = hashlib.sha256(PROMPT.encode('utf-8')).hexdigest()


# --------------------------------------------------------------- the config

def config() -> dict:
    with open(CONFIG, encoding='utf-8') as fh:
        return json.load(fh)


def client_for(route: str, cfg: dict):
    from openai import OpenAI
    spec = cfg['routes'][route]
    key = os.getenv(spec['key_env'])
    if not key:
        raise SystemExit('%s is not set in .env (needed for route %r)'
                         % (spec['key_env'], route))
    kwargs = {'api_key': key, 'timeout': 180.0, 'max_retries': 0}
    if spec['base_url']:
        kwargs['base_url'] = spec['base_url']
    return OpenAI(**kwargs)


def live_prices() -> dict:
    """OpenRouter's current catalogue, so the estimate is not a stale snapshot."""
    key = os.getenv('OPENROUTER_API_KEY')
    req = urllib.request.Request(
        'https://openrouter.ai/api/v1/models',
        headers={'Authorization': 'Bearer ' + key, 'User-Agent': 'engtrace-pilot'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        data = json.load(fh)['data']
    return {m['id']: (float(m['pricing']['prompt']) * 1e6,
                      float(m['pricing']['completion']) * 1e6) for m in data}


def price_of(spec: dict, prices: dict) -> tuple[float, float]:
    if spec['route'] == 'openrouter' and spec['model'] in prices:
        return prices[spec['model']]
    p = spec['price_per_m']
    return p['in'], p['out']


# --------------------------------------------------------------- the corpus

def manifest() -> list[dict]:
    if not os.path.exists(MANIFEST):
        raise SystemExit('no frozen slice - run: python -m evaluator_pilot_17092026.freeze')
    with open(MANIFEST, encoding='utf-8') as fh:
        return [json.loads(ln) for ln in fh]


def require_freeze_verifies() -> None:
    import freeze as fz
    if fz.verify(fz.DEFAULT_MASTER_SEED) != 0:
        raise SystemExit('the frozen slice does not verify - refusing to run. '
                         'Traces must be attached to the items the experts will read.')


def trace_path(key: str) -> str:
    return os.path.join(TRACES, key.replace('/', '_') + '.jsonl')


def existing(key: str) -> dict[str, dict]:
    path = trace_path(key)
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path, encoding='utf-8') as fh:
        for ln in fh:
            try:
                row = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if row.get('ok'):
                out[row['item_id']] = row
    return out


# ------------------------------------------------------------------ the run

def call(cli, spec: dict, question: str, attempts: int = MAX_ATTEMPTS) -> dict:
    """One completion, with retries. Returns a row without item fields.

    Two things here are not boilerplate.

    The token parameter differs by endpoint: OpenAI's reasoning models reject
    `max_tokens` and want `max_completion_tokens`, while OpenRouter and Google's
    compatible endpoint want `max_tokens`. Sending the wrong one is a 400 on
    every call, so the route picks it and a 400 naming the other one switches.

    An empty completion is recorded as a FAILURE even when the API returned 200.
    A reasoning model that spends its whole budget thinking returns `content=""`
    with tokens billed; counted as a success that is a trace nobody can annotate
    and a row that silently lowers every score. This repo's recurring defect is
    a check that is green while measuring nothing.
    """
    prompt = PROMPT.format(question=question)
    param = TOKEN_PARAM.get(spec['route'], 'max_tokens')
    budget = spec.get('max_tokens', MAX_TOKENS)
    last = None
    for attempt in range(1, attempts + 1):
        try:
            t0 = time.time()
            r = cli.chat.completions.create(
                model=spec['model'],
                messages=[{'role': 'user', 'content': prompt}],
                **{param: budget})
            u = r.usage
            det = getattr(u, 'completion_tokens_details', None)
            text = r.choices[0].message.content or ''
            row = {
                'ok': bool(text.strip()),
                'text': text,
                'served_model': getattr(r, 'model', None),
                'finish_reason': r.choices[0].finish_reason,
                'prompt_tokens': getattr(u, 'prompt_tokens', None),
                'completion_tokens': getattr(u, 'completion_tokens', None),
                'reasoning_tokens': getattr(det, 'reasoning_tokens', None) if det else None,
                'seconds': round(time.time() - t0, 2),
                'attempts': attempt,
            }
            if not row['ok']:
                row['error'] = ('empty completion (finish_reason=%s, %s completion tokens) '
                                '- raise max_tokens for this model'
                                % (row['finish_reason'], row['completion_tokens']))
            return row
        except Exception as exc:                                   # noqa: BLE001
            last = '%s: %s' % (type(exc).__name__, exc)
            other = 'max_completion_tokens' if param == 'max_tokens' else 'max_tokens'
            if other in str(exc) and attempt < attempts:
                param = other                     # the endpoint named the one it wants
                continue
            if attempt < attempts:
                time.sleep(BACKOFF * (2 ** (attempt - 1)))
    return {'ok': False, 'error': last, 'attempts': attempts}


def run_model(spec: dict, items: list[dict], cfg: dict, prices: dict,
              workers: int) -> dict:
    key = spec['key']
    done = existing(key)
    todo = [it for it in items if it['item_id'] not in done]
    pin, pout = price_of(spec, prices)
    print('\n%-16s %-38s %d done, %d to run' % (key, spec['model'], len(done), len(todo)))
    if not todo:
        return {'key': key, 'ran': 0, 'failed': 0, 'cost': 0.0}

    cli = client_for(spec['route'], cfg)
    os.makedirs(TRACES, exist_ok=True)
    fh = open(trace_path(key), 'a', encoding='utf-8', newline='\n')
    ran = failed = 0
    cost = 0.0
    try:
        with cf.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(call, cli, spec, it['question']): it
                       for it in todo}
            for fut in cf.as_completed(futures):
                it = futures[fut]
                res = fut.result()
                row = {
                    'item_id': it['item_id'],
                    'template_id': it['template_id'],
                    'branch': it['branch'],
                    'level': it['level'],
                    'item_sha256': it['sha256'],
                    'model_key': key,
                    'model_configured': spec['model'],
                    'route': spec['route'],
                    'prompt_sha256': PROMPT_SHA,
                    'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    **res,
                }
                if res['ok']:
                    row['cost_usd'] = round(
                        ((res['prompt_tokens'] or 0) * pin
                         + (res['completion_tokens'] or 0) * pout) / 1e6, 6)
                    cost += row['cost_usd']
                    ran += 1
                else:
                    failed += 1
                fh.write(json.dumps(row, ensure_ascii=False) + '\n')
                fh.flush()
                n = ran + failed
                if n % 10 == 0 or n == len(todo):
                    print('  %3d/%-3d  ok %3d  failed %3d  $%.3f'
                          % (n, len(todo), ran, failed, cost))
    finally:
        fh.close()
    return {'key': key, 'ran': ran, 'failed': failed, 'cost': cost}


# --------------------------------------------------------------- the modes

def estimate(items, cfg, prices, in_tok=280, out_tok=401):
    print('%-16s %-38s %8s %8s   %s' % ('model', 'id', 'items', 'est $', 'route'))
    total = 0.0
    for spec in cfg['models']:
        pin, pout = price_of(spec, prices)
        todo = len(items) - len(existing(spec['key']))
        c = todo * (in_tok * pin + out_tok * pout) / 1e6
        total += c
        print('%-16s %-38s %8d %8.2f   %s' % (spec['key'], spec['model'], todo, c, spec['route']))
    print('%-16s %-38s %8d %8.2f' % ('TOTAL', '', len(items) * len(cfg['models']), total))
    print('\n(at %d prompt / %d completion tokens per call - the archive mean. '
          'p95 completion is 749, which is about 1.8x the completion half.)' % (in_tok, out_tok))
    return total


def check(cfg) -> int:
    print('== freeze')
    require_freeze_verifies()
    items = manifest()
    print('   %d items, %d templates' % (len(items), len({i['template_id'] for i in items})))

    print('\n== prompt')
    print('   sha256 %s  (%d chars, read from evaluation/run_inference.py)'
          % (PROMPT_SHA[:16], len(PROMPT)))
    want = cfg.get('prompt_sha256_prefix')
    if want and not PROMPT_SHA.startswith(want):
        print('   MISMATCH: models.json expects %s - the deployed prompt changed' % want)
        return 1

    print('\n== keys')
    for route, spec in cfg['routes'].items():
        v = os.getenv(spec['key_env']) or ''
        print('   %-12s %-22s %s' % (route, spec['key_env'],
                                     'set (%d chars)' % len(v) if v else 'MISSING'))

    print('\n== live prices and one call per model')
    prices = live_prices()
    bad = 0
    for spec in cfg['models']:
        pin, pout = price_of(spec, prices)
        # Probe at the model's OWN ceiling, not a token-saving 512. A reasoning
        # model spends its budget thinking before it writes anything, so a small
        # probe tests the ceiling rather than the endpoint - and the ceiling is
        # what the real run would use anyway.
        probe = dict(spec)
        try:
            cli = client_for(spec['route'], cfg)
        except SystemExit as exc:
            bad += 1
            print('   %-16s FAIL %s' % (spec['key'], exc))
            continue
        res = call(cli, probe, 'Reply with the single word OK.', attempts=2)
        if res['ok']:
            print('   %-16s OK   served=%-34s in $%5.2f out $%6.2f  %r'
                  % (spec['key'], res['served_model'], pin, pout, res['text'].strip()[:24]))
            if res.get('served_model') and spec['model'].split('/')[-1] not in str(res['served_model']):
                print('        note: served id differs from the configured id')
        else:
            bad += 1
            print('   %-16s FAIL %s' % (spec['key'], str(res['error'])[:170]))
            if spec.get('fallback'):
                print('        models.json carries a fallback: %s via %s'
                      % (spec['fallback']['model'], spec['fallback']['route']))

    print('\n== estimated spend for the full run')
    estimate(items, cfg, prices)
    print('\n%s' % ('ALL MODELS REACHABLE' if not bad else '%d MODEL(S) UNREACHABLE' % bad))
    return 1 if bad else 0


def status(cfg) -> int:
    items = manifest()
    print('%-16s %6s %6s %8s   %s' % ('model', 'ok', 'missing', 'spent $', 'served'))
    grand = 0.0
    for spec in cfg['models']:
        rows = existing(spec['key'])
        spent = sum(r.get('cost_usd', 0) for r in rows.values())
        grand += spent
        served = Counter(r.get('served_model') for r in rows.values())
        print('%-16s %6d %6d %8.3f   %s'
              % (spec['key'], len(rows), len(items) - len(rows), spent,
                 ', '.join('%s x%d' % (k, v) for k, v in served.most_common(2)) or '-'))
    print('%-16s %6s %6s %8.3f' % ('TOTAL', '', '', grand))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--model', action='append', help='run only this model key')
    ap.add_argument('--limit', type=int, help='first N items only')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()

    cfg = config()
    if args.check:
        return check(cfg)
    if args.status:
        return status(cfg)

    items = manifest()
    if args.limit:
        items = items[:args.limit]
    specs = [s for s in cfg['models'] if not args.model or s['key'] in args.model]
    if not specs:
        raise SystemExit('no model matches %s' % args.model)

    prices = live_prices()
    if args.dry_run:
        estimate(items, cfg, prices)
        return 0

    require_freeze_verifies()
    print('freeze verified, prompt %s, %d items x %d models'
          % (PROMPT_SHA[:12], len(items), len(specs)))
    results = [run_model(s, items, cfg, prices, args.workers) for s in specs]

    print('\n%-16s %6s %8s %10s' % ('model', 'ran', 'failed', 'cost $'))
    for r in results:
        print('%-16s %6d %8d %10.3f' % (r['key'], r['ran'], r['failed'], r['cost']))
    print('%-16s %6d %8d %10.3f' % ('TOTAL', sum(r['ran'] for r in results),
                                    sum(r['failed'] for r in results),
                                    sum(r['cost'] for r in results)))
    failed = sum(r['failed'] for r in results)
    if failed:
        print('\n%d call(s) failed. Re-run the same command: completed items are '
              'skipped, only the failures are retried.' % failed)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
