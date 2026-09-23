"""Layer 1 - the LLM screen: three non-suite judges read every template once.

    python -m template_annotation_23092026.screen.run_screen --check              # keys, model ids, live prices (~$0.01)
    python -m template_annotation_23092026.screen.run_screen --pass 1 --dry-run   # prompts built, tokens, estimate
    python -m template_annotation_23092026.screen.run_screen --pass 1             # the paid run; resumable
    python -m template_annotation_23092026.screen.run_screen --pass 1 --status

What it is. The published template Tribunal (paper section 3.3, Appendix H) read
each template's source plus three generated instances and returned a JSON row of
1-5 scores, a review flag and a one-sentence explanation. This runner keeps that
prompt VERBATIM - it is read out of ai_assisted_quality_assurance/run_ai_tribunal.py
by AST, not retyped - and changes only what the September review record says
must change:

  * the judges are the pilot's non-suite panel (D-088): x-ai/grok-4.6,
    minimax/minimax-m3, xiaomi/mimo-v2.5-pro, all through OpenRouter, so no
    judge shares a family with an evaluated model (Reviewer yAYU, point 1);
  * every judge gets the same call settings (D-089): JSON mode, temperature 0,
    16,384 output tokens, up to three attempts on an empty or malformed reply,
    provider ModelRun excluded for MiniMax M3;
  * the served model id, serving provider, finish reason, tokens and cost are
    written on every row, because the published run recorded none of them and
    its Google judge has since been retired (FINDINGS E0-F4);
  * the three instances come from fixed seeds through the integrity suite's
    generator, so a pass can be reproduced exactly.

A pass is a pass. Rows live under screen/pass<N>/ and a run resumes only within
its own pass; the ceiling is two passes (one before human certification, one
after), and --pass 3 is refused. Nothing here iterates templates against the
judges: Layer 0 is the gate, this is a reader.

Costs are what OpenRouter reports for the call (usage.cost) where it reports one,
and are otherwise priced from the live catalogue at run time; the row says which.
"""
from __future__ import annotations

import argparse
import ast
import collections
import datetime as dt
import hashlib
import inspect
import json
import os
import subprocess
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.core import discover, generate  # noqa: E402

PANEL = [
    {'key': 'grok-4.6',      'model': 'x-ai/grok-4.6',        'family': 'xAI',     'weights': 'closed', 'out_tokens': 3200},
    {'key': 'minimax-m3',    'model': 'minimax/minimax-m3',   'family': 'MiniMax', 'weights': 'open',   'out_tokens': 1100},
    {'key': 'mimo-v2.5-pro', 'model': 'xiaomi/mimo-v2.5-pro', 'family': 'Xiaomi',  'weights': 'open',   'out_tokens': 1100},
]
# out_tokens: completion tokens per reply measured in the 2026-09-23 smoke test (one template,
# three judges). The visible JSON is ~150 tokens; the rest is billed reasoning, which is why a
# flat 200-token assumption under-estimated the pass by about 3x.
CALL = {'response_format': {'type': 'json_object'}, 'temperature': 0.0, 'max_tokens': 16384}
PROVIDER_IGNORE = {'minimax/minimax-m3': ['ModelRun']}
ATTEMPTS = 3
SEEDS = (1001, 1002, 1003)     # the three instances every judge sees, per template
MAX_PASSES = 2
WORKERS = 6
PROMPT_SOURCE = REPO / 'ai_assisted_quality_assurance' / 'run_ai_tribunal.py'
FIELDS = {
    'physical_plausibility_score': int, 'mathematical_correctness_score': int,
    'pedagogical_clarity_score': int, 'confidence_score': int,
    'human_review_flag': bool, 'explanation': str,
}
SCORE_FIELDS = [k for k, t in FIELDS.items() if t is int]


def load_env() -> None:
    try:
        from dotenv import load_dotenv
        load_dotenv(REPO / '.env')
    except ImportError:
        pass


def sha(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def git_head() -> str:
    return subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True,
                          cwd=REPO).stdout.strip()


def prompt_template() -> str:
    """Appendix H's prompt, read from the published script without executing it."""
    src = PROMPT_SOURCE.read_text(encoding='utf8')
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == 'FULL_QA_PROMPT' for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError(f'FULL_QA_PROMPT not found in {PROMPT_SOURCE}')


def build_prompts() -> list[dict]:
    """One prompt per template: function source + three fixed-seed instances."""
    tpl = prompt_template()
    out = []
    for ref in discover():
        fn = ref.load()
        source = inspect.getsource(fn)
        inst = [generate(ref, s, capture=False) for s in SEEDS]
        bad = [i for i in inst if not i.ok]
        if bad:
            out.append({'template_id': ref.template_id, 'branch': ref.branch, 'file': ref.file_path,
                        'error': bad[0].error})
            continue
        text = tpl.format(template_code=source,
                          q1=inst[0].question, s1=inst[0].solution,
                          q2=inst[1].question, s2=inst[1].solution,
                          q3=inst[2].question, s3=inst[2].solution)
        out.append({'template_id': ref.template_id, 'branch': ref.branch, 'file': ref.file_path,
                    'source_sha256': sha(source), 'seeds': list(SEEDS),
                    'prompt_sha256': sha(text), 'prompt': text})
    return out


# ------------------------------------------------------------------ OpenRouter

def client():
    from openai import OpenAI
    key = os.environ.get('OPENROUTER_API_KEY')
    if not key:
        raise SystemExit('OPENROUTER_API_KEY is not set (put it in .env)')
    return OpenAI(api_key=key, base_url='https://openrouter.ai/api/v1', timeout=600.0, max_retries=2)


def live_prices() -> dict:
    """OpenRouter's catalogue as $/M tokens, {model id: (input, output)}."""
    req = urllib.request.Request(
        'https://openrouter.ai/api/v1/models',
        headers={'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY'],
                 'User-Agent': 'engtrace-template-screen'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        data = json.load(fh)['data']
    return {m['id']: (float(m['pricing']['prompt']) * 1e6, float(m['pricing']['completion']) * 1e6)
            for m in data}


def parse_reply(text: str) -> tuple[dict | None, str]:
    """The judge's JSON row, or (None, reason)."""
    t = (text or '').strip()
    if t.startswith('```'):
        t = t.strip('`')
        t = t.split('\n', 1)[-1] if '\n' in t else t
        if t.endswith('```'):
            t = t[:-3]
    try:
        data = json.loads(t)
    except Exception as exc:                                       # noqa: BLE001
        return None, f'not JSON: {type(exc).__name__}'
    if not isinstance(data, dict):
        return None, 'JSON is not an object'
    row = {}
    for k, typ in FIELDS.items():
        if k not in data:
            return None, f'missing {k}'
        v = data[k]
        if typ is int:
            if isinstance(v, bool) or not isinstance(v, (int, float)) or int(v) != v or not 1 <= v <= 5:
                return None, f'{k} not an integer 1-5: {v!r}'
            row[k] = int(v)
        elif typ is bool:
            if isinstance(v, str) and v.lower() in ('true', 'false'):
                v = v.lower() == 'true'
            if not isinstance(v, bool):
                return None, f'{k} not a boolean: {v!r}'
            row[k] = v
        else:
            row[k] = str(v)
    return row, ''


def usage_cost(resp) -> float | None:
    u = getattr(resp, 'usage', None)
    if u is None:
        return None
    c = getattr(u, 'cost', None)
    if c is None:
        c = (getattr(u, 'model_extra', None) or {}).get('cost')
    return float(c) if c is not None else None


def judge_once(cli, judge: dict, prompt: str, prices: dict) -> dict:
    """One judged prompt with the panel's uniform settings; retried on empty or malformed."""
    model = judge['model']
    extra = {'usage': {'include': True}}
    if model in PROVIDER_IGNORE:
        extra['provider'] = {'ignore': PROVIDER_IGNORE[model]}
    tries, out = [], {}
    for attempt in range(1, ATTEMPTS + 1):
        t0 = time.time()
        try:
            r = cli.chat.completions.create(model=model, messages=[{'role': 'user', 'content': prompt}],
                                            extra_body=extra, **CALL)
            ch = r.choices[0]
            u = r.usage
            det = getattr(u, 'completion_tokens_details', None)
            out = {'ok': True, 'text': ch.message.content or '', 'served_model': r.model,
                   'provider': (getattr(r, 'model_extra', None) or {}).get('provider'),
                   'finish_reason': ch.finish_reason,
                   'prompt_tokens': getattr(u, 'prompt_tokens', None),
                   'completion_tokens': getattr(u, 'completion_tokens', None),
                   'reasoning_tokens': getattr(det, 'reasoning_tokens', None) if det else None,
                   'reported_cost': usage_cost(r), 'seconds': round(time.time() - t0, 2)}
        except Exception as exc:                                   # noqa: BLE001
            out = {'ok': False, 'error': f'{type(exc).__name__}: {str(exc)[:300]}',
                   'seconds': round(time.time() - t0, 2)}
        scores, why = (None, 'call failed') if not out['ok'] else parse_reply(out['text'])
        out['scores'] = scores
        out['parse_error'] = why or None
        tries.append({k: out.get(k) for k in ('ok', 'finish_reason', 'completion_tokens', 'reported_cost',
                                              'seconds', 'error', 'provider', 'parse_error')})
        if scores is not None:
            break
    out['attempts'] = tries
    # Every attempt was billed, not only the one that parsed.
    billed_in = sum((out.get('prompt_tokens') or 0) for t in tries if t.get('ok'))
    billed_out = sum(t.get('completion_tokens') or 0 for t in tries)
    reported = [t['reported_cost'] for t in tries if t.get('reported_cost') is not None]
    if reported and len(reported) == sum(1 for t in tries if t.get('ok')):
        out['cost_usd'], out['cost_source'] = round(sum(reported), 6), 'openrouter'
    else:
        pin, pout = prices.get(model, (0.0, 0.0))
        out['cost_usd'], out['cost_source'] = round((billed_in * pin + billed_out * pout) / 1e6, 6), 'catalogue'
    out['billed_prompt_tokens'], out['billed_completion_tokens'] = billed_in, billed_out
    if out.get('ok') and out['scores'] is None:
        out['ok'] = False
        out['error'] = f'no parseable reply in {ATTEMPTS} attempts: {out.get("parse_error")}'
    return out


# ----------------------------------------------------------------- the pass

SMOKE = False   # --smoke: write under _smoke/ so a trial never seeds a real pass's config


def pass_dir(n: int) -> Path:
    return HERE / ('_smoke' if SMOKE else f'pass{n}')


def existing_rows(path: Path) -> dict:
    rows = {}
    if path.exists():
        for ln in path.open(encoding='utf8'):
            if ln.strip():
                r = json.loads(ln)
                if r.get('ok'):
                    rows[(r['template_id'], r['judge'])] = r
    return rows


def carry_rows(from_pass: int, prompts: list[dict]) -> list[dict]:
    """Rows of an earlier pass whose prompt text is byte-identical to the current one.

    A judge's row is a verdict on a specific prompt text. Where a template has not
    changed since the earlier pass, that verdict is a verdict on the corpus that ships,
    so it is carried forward (marked, cost zeroed) rather than bought again; only
    templates whose prompt hash moved are re-judged.
    """
    src = pass_dir(from_pass) / 'replies.jsonl'
    if not src.exists():
        raise SystemExit(f'nothing to carry: {src} does not exist')
    current = {p['template_id']: p['prompt_sha256'] for p in prompts}
    out = []
    for ln in src.open(encoding='utf8'):
        if not ln.strip():
            continue
        r = json.loads(ln)
        if r.get('ok') and current.get(r['template_id']) == r['prompt_sha256']:
            out.append({**r, 'carried_from': from_pass, 'carried_cost_usd': r.get('cost_usd'),
                        'cost_usd': 0.0})
    return out


def run_pass(n: int, dry_run: bool, only: str | None, judge: str | None = None,
             max_usd: float | None = None, carry_from: int | None = None) -> None:
    if not 1 <= n <= MAX_PASSES:
        raise SystemExit(f'the screen is capped at {MAX_PASSES} passes (one before human '
                         f'certification, one after); pass {n} is refused')
    panel = [j for j in PANEL if judge is None or j['key'] == judge]
    if not panel:
        raise SystemExit(f'unknown judge {judge!r}; panel keys: {[j["key"] for j in PANEL]}')
    prompts = build_prompts()
    errors = [p for p in prompts if 'error' in p]
    prompts = [p for p in prompts if 'error' not in p]
    if only:
        prompts = [p for p in prompts if only in p['template_id']]
    tokens = [len(p['prompt']) / 3.6 for p in prompts]
    print(f'pass {n}: {len(prompts)} templates, {len(errors)} generation errors, '
          f'{sum(tokens):,.0f} est. input tokens per judge (chars/3.6)')
    for e in errors:
        print('  generation error:', e['template_id'], e['error'])

    carried = carry_rows(carry_from, prompts) if carry_from else []
    carried_keys = {(r['template_id'], r['judge']) for r in carried}
    if carry_from:
        changed = sorted({p['template_id'] for p in prompts
                          if any((p['template_id'], j['key']) not in carried_keys for j in panel)})
        print(f'carry-forward from pass {carry_from}: {len(carried)} rows on unchanged prompts; '
              f'{len(changed)} templates changed and will be judged: {", ".join(changed)}')
    to_judge = [p for p in prompts if any((p['template_id'], j['key']) not in carried_keys for j in panel)]
    tokens = [len(p['prompt']) / 3.6 for p in to_judge]

    prices = live_prices()
    print('judges (live OpenRouter $/M in, out):')
    est_total = 0.0
    for j in panel:
        pin, pout = prices.get(j['model'], (float('nan'), float('nan')))
        n_j = sum(1 for p in to_judge if (p['template_id'], j['key']) not in carried_keys)
        est = sum(tokens) / 1e6 * pin + j['out_tokens'] * n_j / 1e6 * pout
        est_total += est
        print(f"  {j['key']:14s} {j['model']:24s} {pin:6.2f} {pout:6.2f}  est. ${est:.2f} "
              f"({j['out_tokens']} output tokens per reply, as measured)")
    print(f'estimated pass cost: ${est_total:.2f}')
    if dry_run:
        return

    d = pass_dir(n)
    d.mkdir(parents=True, exist_ok=True)
    cfg = {'pass': n, 'panel': PANEL, 'call': CALL, 'provider_ignore': PROVIDER_IGNORE,
           'attempts': ATTEMPTS, 'seeds': list(SEEDS), 'carried_from': carry_from,
           'carried_rows': len(carried), 'prompt_source': str(PROMPT_SOURCE.relative_to(REPO)),
           'prompt_template_sha256': sha(prompt_template()), 'git_head': git_head(),
           'started': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
           'templates': {p['template_id']: {'source_sha256': p['source_sha256'],
                                            'prompt_sha256': p['prompt_sha256']} for p in prompts}}
    cfg_path = d / 'config.json'
    if cfg_path.exists():
        old = json.loads(cfg_path.read_text(encoding='utf8'))
        moved = [t for t, v in cfg['templates'].items()
                 if t in old['templates'] and old['templates'][t]['prompt_sha256'] != v['prompt_sha256']]
        if moved:
            raise SystemExit(f'{len(moved)} templates changed since this pass started '
                             f'(e.g. {moved[:3]}); a pass reads one corpus. Start the next pass instead.')
    else:
        cfg_path.write_text(json.dumps(cfg, indent=1), encoding='utf8')

    out_path = d / 'replies.jsonl'
    done = existing_rows(out_path)
    if carried:
        fresh = [r for r in carried if (r['template_id'], r['judge']) not in done]
        with out_path.open('a', encoding='utf8') as fh:
            for r in fresh:
                fh.write(json.dumps(r, ensure_ascii=False) + '\n')
        if fresh:
            print(f'wrote {len(fresh)} carried rows into pass {n}')
        done = existing_rows(out_path)
    todo = [(p, j) for p in prompts for j in panel if (p['template_id'], j['key']) not in done]
    already = sum(r.get('cost_usd') or 0.0 for r in done.values())
    print(f'{len(done)} rows already judged (${already:.2f}), {len(todo)} to go'
          + (f'; this run stops once it has spent ${max_usd:.2f}' if max_usd else ''))
    if not todo:
        return
    cli = client()
    spent, ran, failed, dearest = 0.0, 0, 0, 0.0

    def work(p, j):
        res = judge_once(cli, j, p['prompt'], prices)
        return {'pass': n, 'template_id': p['template_id'], 'branch': p['branch'], 'file': p['file'],
                'source_sha256': p['source_sha256'], 'seeds': p['seeds'], 'prompt_sha256': p['prompt_sha256'],
                'judge': j['key'], 'model': j['model'], 'family': j['family'],
                'ts': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), **res}

    # Work goes out in batches of WORKERS so the spend is checked between batches: a
    # cap can stop the run within one batch of the limit, never after the whole pass.
    stopped = False
    with ThreadPoolExecutor(max_workers=WORKERS) as pool, out_path.open('a', encoding='utf8') as fh:
        for start in range(0, len(todo), WORKERS):
            batch = todo[start:start + WORKERS]
            for fut in as_completed([pool.submit(work, p, j) for p, j in batch]):
                row = fut.result()
                fh.write(json.dumps(row, ensure_ascii=False) + '\n')
                fh.flush()
                c = row.get('cost_usd') or 0.0
                spent += c
                dearest = max(dearest, c)
                ran += row['ok']
                failed += not row['ok']
                if c > 0.15:
                    print(f"  note: {row['judge']} on {row['template_id']} cost ${c:.3f} "
                          f"({row.get('completion_tokens')} completion tokens)")
            i = min(start + WORKERS, len(todo))
            print(f'  {i}/{len(todo)} done, {ran} ok, {failed} failed, ${spent:.2f} this run, '
                  f'dearest reply ${dearest:.3f}')
            if max_usd is not None and spent >= max_usd:
                stopped = True
                print(f'  spend cap ${max_usd:.2f} reached after {i} rows; stopping (resumable)')
                break
    print(f'pass {n}: {ran} ok, {failed} failed, ${spent:.2f} this run'
          + (' - STOPPED AT CAP' if stopped else ''))


def status(n: int) -> None:
    path = pass_dir(n) / 'replies.jsonl'
    if not path.exists():
        print(f'pass {n}: nothing yet')
        return
    rows = [json.loads(ln) for ln in path.open(encoding='utf8') if ln.strip()]
    print(f'pass {n}: {len(rows)} rows')
    print(f'{"judge":14s} {"ok":>4s} {"fail":>5s} {"cost $":>8s} {"reported":>9s}  served')
    for j in PANEL:
        js = [r for r in rows if r['judge'] == j['key']]
        ok = [r for r in js if r['ok']]
        served = collections.Counter(r.get('served_model') for r in ok).most_common(2)
        cost = sum(r.get('cost_usd') or 0 for r in js)
        rep = sum(1 for r in js if r.get('cost_source') == 'openrouter')
        print(f"{j['key']:14s} {len(ok):4d} {len(js) - len(ok):5d} {cost:8.3f} {rep:9d}  "
              + ', '.join(f'{k} x{v}' for k, v in served))
    fails = [r for r in rows if not r['ok']]
    for r in fails[:10]:
        print('  failed:', r['template_id'], r['judge'], r.get('error'))
    print(f'total ${sum(r.get("cost_usd") or 0 for r in rows):.3f}')


def check() -> None:
    cli = client()
    prices = live_prices()
    print(f'{len(prices)} models in the OpenRouter catalogue')
    for j in PANEL:
        pin, pout = prices.get(j['model'], (None, None))
        if pin is None:
            print(f"  {j['key']:14s} NOT IN CATALOGUE: {j['model']}")
            continue
        extra = {'usage': {'include': True}}
        if j['model'] in PROVIDER_IGNORE:
            extra['provider'] = {'ignore': PROVIDER_IGNORE[j['model']]}
        t0 = time.time()
        try:
            r = cli.chat.completions.create(
                model=j['model'], messages=[{'role': 'user', 'content': 'Reply with exactly {"ok": true}'}],
                extra_body=extra, response_format=CALL['response_format'], temperature=0.0, max_tokens=64)
            print(f"  {j['key']:14s} OK  served={r.model:<28s} provider={(r.model_extra or {}).get('provider')}"
                  f"  ${pin:.2f}/{pout:.2f} per M  cost={usage_cost(r)}  {time.time() - t0:.1f}s  "
                  f"{(r.choices[0].message.content or '').strip()[:30]!r}")
        except Exception as exc:                                   # noqa: BLE001
            print(f"  {j['key']:14s} FAILED {type(exc).__name__}: {str(exc)[:200]}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--pass', dest='n', type=int, default=1)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--only', default=None, help='substring of a template id, for a smoke test')
    ap.add_argument('--smoke', action='store_true', help='write under _smoke/, not a numbered pass')
    ap.add_argument('--judge', default=None, help='run one judge only (panel key), e.g. minimax-m3')
    ap.add_argument('--max-usd', type=float, default=None,
                    help='stop this run once its reported spend reaches this many dollars')
    ap.add_argument('--carry-from', type=int, default=None,
                    help='carry an earlier pass\'s rows forward where the prompt is unchanged; judge only the rest')
    a = ap.parse_args()
    global SMOKE
    SMOKE = a.smoke
    load_env()
    if a.check:
        check()
    elif a.status:
        status(a.n)
    else:
        run_pass(a.n, a.dry_run, a.only, a.judge, a.max_usd, a.carry_from)


if __name__ == '__main__':
    main()
