"""The full run: the eleven roster models over the frozen pool, and a record of what was run.

    python -m full_run_28092026.run_traces --dry-run              # FREE: the plan and the estimate; no model call
    python -m full_run_28092026.run_traces --status               # FREE: what exists, what it billed
    python -m full_run_28092026.run_traces --check --yes          # BILLS about a cent: one tiny call per model
    python -m full_run_28092026.run_traces --calibrate 20 --yes   # BILLS: 20 items per model, measured lengths
    python -m full_run_28092026.run_traces --yes [--model KEY]    # BILLS: the run, resumable

Every mode that calls a model refuses to start without --yes, and the dry run prints what each would
bill first. Nothing paid runs without the owner's approval of that run.

Adapted from evaluator_pilot_17092026/run_traces.py, which produced the 300 traces the experts
annotated; the call, the prompt and the decoding settings are the same, so the evaluation stack
validated on those traces reads the same kind of output. What changes:

  - ITEMS: the frozen pool (D-116). The pool on disk must match the committed manifest
    (freeze.check_files) before a single call is made.
  - ROSTER: models.json, the eleven of D-110, all through OpenRouter; open-weight models on the
    cheapest endpoint serving fp8 or better, as the pricing document says.
  - WHAT A ROW IS (ANALYSIS_PLAN.md, D-117). Each (item, model) ends in one of three states:
      answered          text came back; the answer check decides what it is worth, and a
                        truncated answer is still scored on what it states
      empty             the model returned no text: unusable, scored 0, not called again
      service_failure   an HTTP error, timeout or provider fault that persisted through the
                        retries: reported as missing, never scored, and called again on the
                        next run
  - COST: the billed cost of every call is read from the response (usage.cost), so the spend is
    the bill, not an estimate; Gemini's unreported thinking tokens cannot hide in it.

Traces go to traces/<key>.jsonl beside this file, gitignored: they restate the pool's questions.

VARIANT RUNS (D-141). `--variant` runs the same models, prompt, settings and states on another item
set, into traces/<variant>/<key>.jsonl, so score.py scores it as that variant:
  paraphrase        the 450 items of subsamples.PARAPHRASE, each question replaced by its paraphrase
                    from paraphrase/pool.jsonl (paraphrase.py), checked against the committed
                    paraphrase/manifest.jsonl before a call; the row records the paraphrase's hash as
                    item_sha256 and the original's as original_sha256
  repeat1..repeat3  the 300 items of subsamples.REPEAT with their original questions, one run each
The dry run estimates a variant from each model's own bills on the same items in the main run.

    python -m full_run_28092026.run_traces --variant paraphrase --dry-run          # FREE
    python -m full_run_28092026.run_traces --variant repeat1 --model gemma-4-26b-a4b --dry-run
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
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dotenv import load_dotenv  # noqa: E402

load_dotenv(REPO / '.env')

from full_run_28092026 import freeze, subsamples  # noqa: E402

CONFIG = HERE / 'models.json'
TRACES = HERE / 'traces'
PARAPHRASES = HERE / 'paraphrase'
VARIANTS = ('main', 'paraphrase') + subsamples.REPEAT_VARIANTS
DEPLOYED_RUNNER = REPO / 'evaluation' / 'run_inference.py'
PILOT_PROMPT_PREFIX = 'c2bcb87984c4e50b'   # evaluator_pilot_17092026/models.json: the annotated traces' prompt
MAX_ATTEMPTS = 4
BACKOFF = 5                                 # seconds, doubled per attempt
DONE = ('answered', 'empty')                # states a re-run does not call again


def deployed_prompt() -> str:
    src = DEPLOYED_RUNNER.read_text(encoding='utf-8')
    m = re.search(r'PROMPT_TEMPLATE = """(.*?)"""', src, re.S)
    if not m:
        raise SystemExit(f'cannot find PROMPT_TEMPLATE in {DEPLOYED_RUNNER}')
    return m.group(1)


PROMPT = deployed_prompt()
PROMPT_SHA = hashlib.sha256(PROMPT.encode('utf-8')).hexdigest()


def config() -> dict:
    return json.loads(CONFIG.read_text(encoding='utf-8'))


def items() -> list[dict]:
    """The pool on disk, joined to the committed manifest. Refuses if they disagree."""
    if not freeze.POOL_DIR.exists():
        raise SystemExit('full_run_28092026/pool/ not found: restore it from the private backup')
    if freeze.check_files() != 0:
        raise SystemExit('the pool on disk does not match manifest.jsonl - refusing to run')
    rows = {json.loads(l)['item_id']: json.loads(l)
            for l in freeze.MANIFEST.read_text(encoding='utf-8').splitlines()}
    out = []
    for path in sorted(freeze.POOL_DIR.rglob('*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            m = rows[r['item_id']]
            out.append({'item_id': r['item_id'], 'template_id': m['template_id'],
                        'branch': m['branch'], 'level': m['level'],
                        'answer_type': m['answer_type'], 'sha256': m['sha256'],
                        'question': r['question']})
    out.sort(key=lambda r: (r['template_id'], int(r['item_id'].split('#')[1])))
    return out


def variant_items(variant: str, its: list[dict]) -> list[dict]:
    """The items a variant runs: a subsample of the pool, with paraphrased questions for `paraphrase`."""
    if variant == 'main':
        return its
    by_id = {it['item_id']: it for it in its}
    if variant in subsamples.REPEAT_VARIANTS:
        return [by_id[i] for i in subsamples.repeat_ids()]
    pool, manifest = PARAPHRASES / 'pool.jsonl', PARAPHRASES / 'manifest.jsonl'
    if not pool.exists() or not manifest.exists():
        raise SystemExit('paraphrase/pool.jsonl or its manifest is missing: run paraphrase.py first')
    want = {r['item_id']: r for r in map(json.loads, manifest.read_text(encoding='utf-8').splitlines())
            if r['passed']}
    out = []
    for r in map(json.loads, pool.read_text(encoding='utf-8').splitlines()):
        m = want.get(r['item_id'])
        if m is None:
            continue
        sha = hashlib.sha256(r['question'].encode('utf-8')).hexdigest()
        if sha != m['sha256'] or by_id[r['item_id']]['sha256'] != m['original_sha256']:
            raise SystemExit(f"{r['item_id']}: the paraphrase pool does not match its manifest - refusing to run")
        out.append({**by_id[r['item_id']], 'question': r['question'], 'sha256': sha,
                    'original_sha256': m['original_sha256']})
    if len(out) != len(want):
        raise SystemExit(f'the manifest passes {len(want)} paraphrases but the pool holds {len(out)}')
    return out


def trace_path(key: str, variant: str = 'main') -> Path:
    return TRACES / f'{key}.jsonl' if variant == 'main' else TRACES / variant / f'{key}.jsonl'


def existing(key: str, variant: str = 'main') -> dict[str, dict]:
    path = trace_path(key, variant)
    if not path.exists():
        return {}
    out = {}
    for ln in path.read_text(encoding='utf-8').splitlines():
        try:
            row = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if row.get('status') in DONE:
            out[row['item_id']] = row
    return out


def request_params(spec: dict, cfg: dict) -> dict:
    p = {'max_tokens': spec.get('max_tokens', cfg['max_tokens'])}
    key = 'provider_open_weight' if spec['weights'] == 'open' else 'provider_closed_weight'
    p['extra_body'] = {'provider': cfg[key]}
    return p


def client(cfg: dict):
    from openai import OpenAI
    key = os.getenv(cfg['route']['key_env'])
    if not key:
        raise SystemExit(f"{cfg['route']['key_env']} is not set in .env")
    return OpenAI(api_key=key, base_url=cfg['route']['base_url'], timeout=600.0, max_retries=0)


def call(cli, spec: dict, cfg: dict, question: str, attempts: int = MAX_ATTEMPTS) -> dict:
    """One completion with retries, ending in answered / empty / service_failure."""
    prompt = PROMPT.format(question=question)
    params = request_params(spec, cfg)
    last, empty_row = None, None
    for attempt in range(1, attempts + 1):
        try:
            t0 = time.time()
            r = cli.chat.completions.create(model=spec['model'],
                                            messages=[{'role': 'user', 'content': prompt}], **params)
            u = r.usage
            det = getattr(u, 'completion_tokens_details', None)
            msg = r.choices[0].message
            extra = getattr(msg, 'model_extra', None) or {}
            text = msg.content or ''
            if r.choices[0].finish_reason == 'error':
                # A provider fault reported inside a 200: not the model's answer (D-148). Seven such
                # rows in the main run were scored on what they state and are reported as a count.
                raise RuntimeError('provider reported finish_reason=error')
            row = {
                'text': text,
                'reasoning': getattr(msg, 'reasoning', None) or extra.get('reasoning') or '',
                'served_model': getattr(r, 'model', None),
                'provider': (getattr(r, 'model_extra', None) or {}).get('provider'),
                'finish_reason': r.choices[0].finish_reason,
                'prompt_tokens': getattr(u, 'prompt_tokens', None),
                'completion_tokens': getattr(u, 'completion_tokens', None),
                'reasoning_tokens': getattr(det, 'reasoning_tokens', None) if det else None,
                'billed_usd': (getattr(u, 'model_extra', None) or {}).get('cost'),
                'seconds': round(time.time() - t0, 2),
                'attempts': attempt,
            }
            if text.strip():
                return {'status': 'answered', **row}
            # No text. Out of budget is the model's doing and final; any other empty
            # completion is asked again before it is recorded as empty.
            empty_row = {'status': 'empty', **row}
            if row['finish_reason'] == 'length':
                return empty_row
        except Exception as exc:                                   # noqa: BLE001
            last = f'{type(exc).__name__}: {exc}'[:400]
        if attempt < attempts:
            time.sleep(BACKOFF * (2 ** (attempt - 1)))
    return empty_row or {'status': 'service_failure', 'error': last, 'attempts': attempts}


def run_model(spec, todo, cfg, workers, label='run', variant='main'):
    key = spec['key']
    cli = client(cfg)
    trace_path(key, variant).parent.mkdir(parents=True, exist_ok=True)
    params = request_params(spec, cfg)
    n, counts, billed = 0, Counter(), 0.0
    with open(trace_path(key, variant), 'a', encoding='utf-8', newline='\n') as fh, \
            cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(call, cli, spec, cfg, it['question']): it for it in todo}
        for fut in cf.as_completed(futs):
            it = futs[fut]
            res = fut.result()
            row = {'item_id': it['item_id'], 'template_id': it['template_id'],
                   'branch': it['branch'], 'level': it['level'], 'answer_type': it['answer_type'],
                   'item_sha256': it['sha256'], 'model_key': key, 'model_configured': spec['model'],
                   'prompt_sha256': PROMPT_SHA, 'request': params, 'mode': label,
                   'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **res}
            if variant != 'main':
                row['variant'] = variant
            if 'original_sha256' in it:
                row['original_sha256'] = it['original_sha256']
            fh.write(json.dumps(row, ensure_ascii=False) + '\n')
            fh.flush()
            n += 1
            counts[res['status']] += 1
            billed += res.get('billed_usd') or 0.0
            if n % 25 == 0 or n == len(todo):
                print(f'  {key:22s} {n:5d}/{len(todo):<5d} {dict(counts)}  billed ${billed:.3f}')
    return counts, billed


# ------------------------------------------------------------------ estimates

def catalogue() -> dict:
    """OpenRouter's public model list: live prices and ceilings. No key is sent; no model is called."""
    req = urllib.request.Request('https://openrouter.ai/api/v1/models',
                                 headers={'User-Agent': 'engtrace-full-run'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return {m['id']: m for m in json.load(fh)['data']}


def per_item_cost(price: dict, tok: dict) -> float:
    return (tok['in'] * price['in'] + tok['out'] * price['out']) / 1e6


def routed_endpoint(spec: dict, cfg: dict, tok: dict):
    """The endpoint routing will pick, read from OpenRouter's public endpoint list.

    Open weights: the cheapest endpoint whose quantization is fp8 or better (the pricing
    document's rule). Closed: the cheapest listed. Returns (endpoint, price, eligible count),
    or (None, None, 0) when the list is unreachable or holds no eligible endpoint.
    """
    url = f"https://openrouter.ai/api/v1/models/{spec['model']}/endpoints"
    req = urllib.request.Request(url, headers={'User-Agent': 'engtrace-full-run'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        eps = (json.load(fh).get('data') or {}).get('endpoints') or []
    allowed = set(cfg['provider_open_weight']['quantizations'])
    if spec['weights'] == 'open':
        eps = [e for e in eps if (e.get('quantization') or 'unknown') in allowed]
    best = None
    for e in eps:
        p = {'in': float(e['pricing']['prompt']) * 1e6, 'out': float(e['pricing']['completion']) * 1e6}
        if best is None or per_item_cost(p, tok) < per_item_cost(best[1], tok):
            best = (e, p)
    return (best[0], best[1], len(eps)) if best else (None, None, 0)


def dry_run(cfg, its, specs) -> int:
    print(f'pool: {len(its)} items, {len({i["template_id"] for i in its})} templates, '
          f'matching the committed manifest')
    ok = PROMPT_SHA.startswith(PILOT_PROMPT_PREFIX)
    print(f'prompt: sha256 {PROMPT_SHA[:16]}, read from evaluation/run_inference.py; '
          f'{"the same" if ok else "NOT the same"} prompt as the pilot traces')
    tok = cfg['pricing_basis_tokens']
    print(f'\nbasis: {tok["in"]} input and {tok["out"]} output tokens per item ({tok["source"]})')
    print('doc $: the pricing document\'s prices. routed $: the endpoint routing picks today, '
          'from OpenRouter\'s public endpoint list.')
    print(f'{"model":22s} {"items":>6s} {"doc $":>8s} {"routed $":>9s} {"cap":>7s} '
          f'{"eligible":>8s}  endpoint')
    tot_doc = tot_routed = 0.0
    short = []
    for s in specs:
        todo = len(its) - len(existing(s['key']))
        doc = todo * per_item_cost(s['price_per_m'], tok)
        ceiling = s.get('max_tokens', cfg['max_tokens'])
        try:
            ep, price, n_ok = routed_endpoint(s, cfg, tok)
        except Exception as exc:                                   # noqa: BLE001
            ep, price, n_ok = None, None, 0
            print(f'  ({s["key"]}: endpoint list unreachable, {type(exc).__name__})')
        routed = todo * per_item_cost(price, tok) if price else None
        cap = ep.get('max_completion_tokens') if ep else None
        where = (f"{ep.get('provider_name') or ep.get('name')} {ep.get('quantization') or ''}".strip()
                 if ep else 'NO ELIGIBLE ENDPOINT')
        if cap and cap < ceiling:
            short.append((s['key'], cap))
        tot_doc += doc
        tot_routed += routed or 0.0
        print(f'{s["key"]:22s} {todo:6d} {doc:8.2f} {"" if routed is None else f"{routed:9.2f}":>9s} '
              f'{"" if cap is None else cap:>7} {n_ok:8d}  {where}')
    print(f'{"TOTAL":22s} {len(its) * len(specs):6d} {tot_doc:8.2f} {tot_routed:9.2f}')
    for key, cap in short:
        print(f'  {key}: the routed endpoint caps output at {cap}, below the {cfg["max_tokens"]} ceiling')
    per_model_item = sum(per_item_cost(s['price_per_m'], tok) for s in specs)
    print(f'\nwhat each paid mode would bill on the same basis:')
    print(f'  --check          one tiny call per model: about $0.01')
    print(f'  --calibrate 20   {20 * len(specs)} calls: about ${20 * per_model_item:.2f}')
    print(f'  the run          {len(its) * len(specs)} calls: about ${tot_doc:.2f}')
    print('\nThe basis assumes every model writes as much per item as GPT-5 did on the pilot;')
    print('a reasoning-heavy model can exceed it (DeepSeek R1 wrote a median 9,658 there).')
    print('--calibrate replaces the assumption with measured lengths. Nothing was called.')
    return 0


def variant_dry_run(cfg, its, specs, variant) -> int:
    """A variant's plan and estimate: each model's own billed cost per item on the same items in the
    main run, which already carries its lengths, its routing and its prices."""
    print(f'variant {variant}: {len(its)} items, {len({i["template_id"] for i in its})} templates; '
          f'traces go to traces/{variant}/')
    print(f'{"model":22s} {"to run":>7s} {"main $ on these":>16s} {"estimate $":>11s}')
    total = 0.0
    for s in specs:
        todo = [it for it in its if it['item_id'] not in existing(s['key'], variant)]
        main = existing(s['key'])
        bills = [main[it['item_id']].get('billed_usd') or 0.0 for it in its if it['item_id'] in main]
        per = sum(bills) / len(bills) if bills else None
        est = per * len(todo) if per is not None else None
        total += est or 0.0
        print(f'{s["key"]:22s} {len(todo):7d} {"" if per is None else f"{sum(bills):16.3f}"} '
              f'{"no main-run bills" if est is None else f"{est:11.3f}"}')
    print(f'{"TOTAL":22s} {"":7s} {"":16s} {total:11.3f}')
    print('\nThe estimate repeats the main run\'s spend on the same items; a paraphrase is about as long '
          'as its original. Nothing was called.')
    return 0


def calibration_sample(its, n):
    """n items per model, spread over the templates in a fixed order: every k-th item."""
    step = max(1, len(its) // n)
    return its[::step][:n]


def status(cfg, its, specs, variant='main') -> int:
    print(f'{"model":22s} {"answered":>9s} {"empty":>6s} {"missing":>8s} {"billed $":>9s} '
          f'{"out tok/item":>13s}  served')
    grand = 0.0
    for s in specs:
        rows = list(existing(s['key'], variant).values())
        c = Counter(r['status'] for r in rows)
        billed = sum(r.get('billed_usd') or 0.0 for r in rows)
        grand += billed
        outs = [r['completion_tokens'] for r in rows if r.get('completion_tokens')]
        served = Counter(r.get('served_model') for r in rows).most_common(1)
        print(f'{s["key"]:22s} {c["answered"]:9d} {c["empty"]:6d} {len(its) - len(rows):8d} '
              f'{billed:9.3f} {(sum(outs) / len(outs)) if outs else 0:13.0f}  '
              f'{served[0][0] if served else "-"}')
    print(f'{"TOTAL":22s} {"":9s} {"":6s} {"":8s} {grand:9.3f}')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--calibrate', type=int, metavar='N')
    ap.add_argument('--model', action='append', help='only this model key (repeatable)')
    ap.add_argument('--variant', default='main', choices=VARIANTS)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--yes', action='store_true', help='required for any mode that bills')
    a = ap.parse_args()

    cfg = config()
    # An entry with "run": false is inert (D-148): no mode calls it unless --model names it, and the
    # dry run and the status say so instead of listing it as missing.
    inert = {s['key']: s.get('run_reason', '') for s in cfg['models'] if s.get('run') is False}
    specs = [s for s in cfg['models'] if (a.model and s['key'] in a.model) or (not a.model and s['key'] not in inert)]
    if not specs:
        raise SystemExit(f'no model matches {a.model}')
    for k, why in inert.items():
        if not a.model or k not in a.model:
            print(f'{k}: inert (run: false): {why}')
    if a.variant != 'main' and a.calibrate:
        raise SystemExit('--calibrate is for the main run; a variant is estimated from its bills')
    if a.variant in subsamples.REPEAT_VARIANTS and not a.model:
        raise SystemExit('a repeat runs one chosen model: name it with --model')
    its = variant_items(a.variant, items())
    if a.dry_run:
        return dry_run(cfg, its, specs) if a.variant == 'main' else variant_dry_run(cfg, its, specs, a.variant)
    if a.status:
        return status(cfg, its, specs, a.variant)
    if not a.yes:
        raise SystemExit('this mode bills: re-run with --yes once the spend is approved '
                         '(see --dry-run for the estimate)')
    if not PROMPT_SHA.startswith(PILOT_PROMPT_PREFIX):
        raise SystemExit('the deployed prompt differs from the one the pilot traces used - refusing')
    if a.check:
        bad = 0
        for s in specs:
            res = call(client(cfg), s, cfg, 'Reply with the single word OK.', attempts=2)
            ok = res['status'] == 'answered'
            bad += not ok
            print(f'  {s["key"]:22s} {res["status"]:16s} served={res.get("served_model")} '
                  f'provider={res.get("provider")} billed=${res.get("billed_usd") or 0:.5f}')
        return 1 if bad else 0
    grand = 0.0
    tok = cfg['pricing_basis_tokens']
    for s in specs:
        if routed_endpoint(s, cfg, tok)[0] is None:
            print(f'\n{s["key"]}: SKIPPED - no endpoint meets the routing rule (see --dry-run)')
            continue
        pool = calibration_sample(its, a.calibrate) if a.calibrate else its
        done = existing(s['key'], a.variant)
        todo = [it for it in pool if it['item_id'] not in done]
        print(f'\n{s["key"]}: {len(done)} done, {len(todo)} to run')
        if todo:
            counts, billed = run_model(s, todo, cfg, a.workers,
                                       label='calibrate' if a.calibrate else ('run' if a.variant == 'main' else a.variant),
                                       variant=a.variant)
            grand += billed
    print(f'\nbilled this invocation: ${grand:.3f}. Re-run the same command to resume; '
          f'only service failures and unrun items are called.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
