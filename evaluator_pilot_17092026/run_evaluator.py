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


def verified_traces(freeze_check: str = 'rebuild', cohort: str = 'gold'):
    """The traces verify_traces passes, or stop.

    freeze_check='rebuild' regenerates the slice from the templates and compares
    bytes - the full check, and the one to use wherever the repo is present.
    'hash' only confirms the manifest's SHA-256 equals the one FREEZE.json
    recorded. It is for a Kaggle kernel, which receives the manifest but not the
    2,250-record testset or the template tree the rebuild needs; the bundle is
    staged only after a full rebuild passes locally, so the hash then proves the
    kernel is reading those same bytes. The mode used is written on every row.
    """
    if freeze_check == 'rebuild':
        import freeze as fz
        if fz.verify(fz.DEFAULT_MASTER_SEED) != 0:
            raise SystemExit('the frozen slice does not verify - refusing to score')
    else:
        with open(rt.MANIFEST, 'rb') as fh:
            got = hashlib.sha256(fh.read()).hexdigest()
        with open(os.path.join(rt.SLICE, 'FREEZE.json'), encoding='utf-8') as fh:
            want = json.load(fh)['manifest_sha256']
        if got != want:
            raise SystemExit('manifest sha256 %s != FREEZE.json %s - refusing to score'
                             % (got[:16], want[:16]))
        print('FREEZE HASH OK - manifest sha256 %s matches FREEZE.json' % got[:16])
    items, traces = vt.load()
    keys = cohort_keys(cohort)
    bad = vt.check(items, traces, keys)
    other = vt.check(items, traces, set(traces) - keys)
    if other:
        print('note: %d issue(s) in columns outside cohort %r, not gating this run'
              % (len(other), cohort))
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


def cohort_keys(cohort: str) -> set:
    """Which models a run scores. Default 'gold': the five the experts annotate.

    The robustness cohort exists to be checked mechanically, not to be scored
    into the comparison, and judging it would cost real money per model for
    numbers no human label backs. Scoring it is therefore opt-in.
    """
    cfg = rt.config()
    return {m['key'] for m in cfg['models']
            if cohort == 'all' or m.get('cohort', 'gold') == cohort}


def plan(items, traces, evaluator, csha, limit=None, models=None, cohort='gold'):
    work = []
    allowed = cohort_keys(cohort)
    for key in sorted(traces):
        if key not in allowed:
            continue
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


def run(evaluator: str, dry_run: bool, limit, models, freeze_check='rebuild', cohort='gold'):
    mod = importlib.import_module(EVALUATORS[evaluator])
    items, traces = verified_traces(freeze_check, cohort)
    cfg = mod.config()
    csha = config_sha(cfg)
    work = plan(items, traces, evaluator if not dry_run else evaluator + '_dry',
                csha, limit, models, cohort)
    print('%s: %s' % (mod.ID, mod.DESCRIPTION))
    print('config %s   cohort %s   %d (item, model) pairs to score%s'
          % (csha[:12], cohort, len(work), '  [DRY RUN - no judge is called]' if dry_run else ''))
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
               'freeze_check': freeze_check, 'cohort': cohort,
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
    try:
        prices = rt.live_prices()
    except Exception:                                              # noqa: BLE001
        # No key on a Kaggle kernel, by design. Snapshot prices from 2026-09-17;
        # re-run the estimate locally for live ones.
        prices = {'openai/gpt-5': (1.25, 10.0), 'anthropic/claude-opus-4.5': (5.0, 25.0)}
        print('\n(no OpenRouter key here - estimate uses 2026-09-17 snapshot prices)')
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


EXACT = ('recall', 'precision', 'recovered_f1', 'final_answer_acc', 'rouge2', 'rougeL', 'rougeLsum')
BERT_TOL = 1e-4
MIN_REFERENCE = 8


def _rows(d: str) -> dict:
    out = {}
    if not os.path.isdir(d):
        return out
    for fn in os.listdir(d):
        if fn.endswith('.jsonl'):
            for ln in open(os.path.join(d, fn), encoding='utf-8'):
                r = json.loads(ln)
                if not r.get('error'):
                    out[(r['model_key'], r['item_id'])] = r
    return out


def import_kaggle(evaluator: str, src: str) -> int:
    """Take a kernel's dry run only after proving it matches this machine.

    The GPU kernel and this CPU run the same pinned scoring stack on different
    hardware and torch builds. Agreement is expected and is NOT assumed: every
    (model, item) scored in both places is compared - Tier 1 metrics and ROUGE
    exactly, BERTScore within 1e-4, the Tribunal trigger identically - and a
    single disagreement refuses the import. Only then is the kernel's Tier 1
    cache merged, so the paid run here reuses GPU work instead of repeating it.
    """
    mod = importlib.import_module(EVALUATORS[evaluator])
    csha = config_sha(mod.config())
    k_rows = _rows(os.path.join(src, 'scores', evaluator + '_dry'))
    l_dir = os.path.join(SCORES, evaluator + '_dry')
    l_rows = {k: r for k, r in _rows(l_dir).items() if not (r['meta'].get('compute') or {}).get('cuda')}
    k_cache = os.path.join(src, 'scores', '_cache', 'e0_tier1.jsonl')

    problems = []
    if len(k_rows) != 300:
        problems.append('kernel scored %d of 300' % len(k_rows))
    wrong_cfg = sum(1 for r in k_rows.values() if r['config_sha256'] != csha)
    if wrong_cfg:
        problems.append('%d kernel rows carry a different config (library pins drifted?)' % wrong_cfg)
    not_gpu = sum(1 for r in k_rows.values() if not (r['meta'].get('compute') or {}).get('cuda'))
    if not_gpu:
        problems.append('%d kernel rows were NOT scored on a GPU' % not_gpu)
    sfail = sum(1 for r in k_rows.values() if r['meta'].get('scorer_failures'))
    if sfail:
        problems.append('%d kernel rows have a BERTScore failure' % sfail)

    shared = sorted(set(k_rows) & set(l_rows))
    worst_bert, worst_exact = 0.0, 0.0
    for k in shared:
        a, b = k_rows[k], l_rows[k]
        for f in EXACT:
            d = abs(a['scores'][f] - b['scores'][f])
            worst_exact = max(worst_exact, d)
            if d > 1e-9:
                problems.append('%s/%s %s: kernel %.6f vs local %.6f' % (k[0], k[1], f, a['scores'][f], b['scores'][f]))
        d = abs(a['scores']['bertscore'] - b['scores']['bertscore'])
        worst_bert = max(worst_bert, d)
        if d > BERT_TOL:
            problems.append('%s/%s bertscore: kernel %.6f vs local %.6f' % (k[0], k[1], a['scores']['bertscore'], b['scores']['bertscore']))
        for f in ('triggered',):
            if a.get(f) != b.get(f):
                problems.append('%s/%s %s differs' % (k[0], k[1], f))
        if a['meta'].get('tribunal_reached_judges') != b['meta'].get('tribunal_reached_judges'):
            problems.append('%s/%s tribunal_reached_judges differs' % k)
    if len(shared) < MIN_REFERENCE:
        problems.append('only %d rows scored in both places; need %d to vouch for the rest'
                        % (len(shared), MIN_REFERENCE))

    dev = Counter((r['meta'].get('compute') or {}).get('device') for r in k_rows.values())
    print('kernel rows %d on %s' % (len(k_rows), dict(dev)))
    print('reference rows scored on both GPU and this CPU: %d' % len(shared))
    print('  worst |diff|: Tier 1 + ROUGE %.2e   BERTScore %.2e (tolerance %.0e)'
          % (worst_exact, worst_bert, BERT_TOL))
    if problems:
        print('IMPORT REFUSED - %d problem(s):' % len(problems))
        for p in problems[:20]:
            print('  ' + p)
        return 1

    # Merge: the local rows become the reference record, the kernel rows the dry run.
    ref = os.path.join(SCORES, evaluator + '_dry_cpu_reference')
    if os.path.isdir(l_dir):
        os.makedirs(ref, exist_ok=True)
        for fn in os.listdir(l_dir):
            os.replace(os.path.join(l_dir, fn), os.path.join(ref, fn))
    os.makedirs(l_dir, exist_ok=True)
    for fn in os.listdir(os.path.join(src, 'scores', evaluator + '_dry')):
        with open(os.path.join(src, 'scores', evaluator + '_dry', fn), encoding='utf-8') as fi, \
             open(os.path.join(l_dir, fn), 'w', encoding='utf-8', newline='\n') as fo:
            fo.write(fi.read())
    added = 0
    local_cache = os.path.join(SCORES, '_cache', 'e0_tier1.jsonl')
    os.makedirs(os.path.dirname(local_cache), exist_ok=True)
    have = set()
    if os.path.exists(local_cache):
        have = {json.loads(ln)['k'] for ln in open(local_cache, encoding='utf-8') if ln.strip()}
    with open(local_cache, 'a', encoding='utf-8', newline='\n') as fo:
        for ln in open(k_cache, encoding='utf-8'):
            if ln.strip() and json.loads(ln)['k'] not in have:
                fo.write(ln if ln.endswith('\n') else ln + '\n')
                added += 1
    print('IMPORTED - %d kernel rows are now the dry run; %d Tier 1 cache entries merged; '
          'the CPU rows are kept in %s' % (len(k_rows), added, os.path.relpath(ref, _HERE)))
    dry_estimate(evaluator + '_dry')
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('evaluator', choices=sorted(EVALUATORS))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--smoke', type=int, metavar='N', help='score N entries for real, then stop')
    ap.add_argument('--model', action='append')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--freeze-check', choices=('rebuild', 'hash'), default='rebuild')
    ap.add_argument('--cohort', choices=('gold', 'robustness', 'all'), default='gold',
                    help='gold = the five the experts annotate (default)')
    ap.add_argument('--import-kaggle', metavar='DIR',
                    help="a kernel's downloaded output; merged only if it matches local reference rows")
    args = ap.parse_args()
    if args.status:
        return status(args.evaluator)
    if args.import_kaggle:
        return import_kaggle(args.evaluator, args.import_kaggle)
    return run(args.evaluator, args.dry_run, args.smoke, args.model, args.freeze_check,
               args.cohort)


if __name__ == '__main__':
    sys.exit(main())
