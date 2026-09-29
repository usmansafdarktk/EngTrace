"""E5 over the full run: MiMo-V2.5-Pro on the milestones E3 did not find (D-105, D-110), every reply
stored so none is bought twice (D-142).

    python -m full_run_28092026.judge --dry-run [--variant main]       # FREE: the calls, the cost, the time
    python -m full_run_28092026.judge --validate                       # FREE: the pilot's 300 traces, from its stored replies
    python -m full_run_28092026.judge --yes [--model KEY] [--max-usd 100] [--workers 16]   # BILLS
    python -m full_run_28092026.judge --score [--variant main]         # FREE: the E5 store from the stored replies

WHAT IS SENT is the pilot's E5, imported unchanged from evaluators/e5_hybrid.py: for each answered
trace whose item has milestones E3 did not reach (score.py's store), one call holding the question,
the trace and each missed milestone's id and value, not the gold solution. The question is the
original item's, for a paraphrase variant too, as score.py scores it.

THE CALL is E1's (e1_panel.CALL: JSON mode, temperature 0, 16,384 tokens), up to three attempts, an
empty or malformed reply asked again (e1_panel.malformed), at OpenRouter's default routing, as the
pilot called it. The key is the pilot's, the SHA-256 of (model, settings, prompt), so a prompt already
answered is read from the store, never bought again. The client's timeout is 300 s with no SDK
retries, so one fetch (three attempts) stays inside judge_calls.DEADLINE; the pilot's 600 s with two
SDK retries could outlive it (D-148).

WHAT IS WRITTEN. scores/<variant>/e5/<model>.jsonl, one row per trace:
  - each milestone's source: e3, REACHED, NOT_NEEDED, MISSING, or UNJUDGED when the reply omits it;
  - e5_strict, E3's milestones plus REACHED over those required: the score;
  - e5_lenient, NOT_NEEDED credited too, a diagnostic only, because RESULTS_E5 found the judge excuses
    a quarter of fabricated values that way;
  - the judged fraction, and the reply's provider, tokens and billed cost.
An item with no milestones has no coverage (None), as in the store (D-120). The replies are kept in
scores/_judge/e5_replies.jsonl; both files are gitignored with the store.

THE ESTIMATE counts the calls and their prompt lengths from the store. It takes tokens per character
and completion tokens per call from the pilot's own E5 replies, and prices from OpenRouter's public
endpoint list: at Xiaomi's own endpoint, with the cheapest and dearest beside it. The time is at the
chosen workers, from the pilot's call times.
"""
from __future__ import annotations

import argparse
import collections
import io
import json
import statistics
import sys
import time
import urllib.request
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PILOT = REPO / 'evaluator_pilot_17092026'
for p in (str(REPO), str(PILOT / 'evaluators'), str(PILOT / 'annotation'), str(PILOT / 'analysis'),
          str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

with redirect_stdout(io.StringIO()):
    import e1_panel as e1        # noqa: E402  CALL, ATTEMPTS, malformed, Replies.key
    import e5_hybrid as e5       # noqa: E402  PROMPT, prompt_for, parse, JUDGE

from full_run_28092026 import judge_calls as jc, score  # noqa: E402

REPLIES = score.SCORES / '_judge' / 'e5_replies.jsonl'
PILOT_REPLIES = PILOT / 'scores' / '_cache' / 'e5_judge_replies.jsonl'
PILOT_ROWS = PILOT / 'scores' / 'e5'
VERDICTS = e5.VERDICTS


# ------------------------------------------------------------------ one trace

def job(question: str, text: str, milestones: list[dict], reached: list[bool]):
    """(missed milestones, prompt, key) for a trace, or None when E3 reached every milestone."""
    missed = [m for m, r in zip(milestones, reached) if not r]
    if not missed:
        return None
    prompt = e5.prompt_for({'question': question}, text, missed)
    return missed, prompt, e1.Replies.key(e5.JUDGE, prompt)


def fetch(prompt: str, cli=None) -> dict:
    """e1_panel.fetch's call and retry rule, with the billed cost summed over the attempts."""
    cli = cli or jc.client()
    tries, out = [], {}
    for _ in range(e1.ATTEMPTS):
        out = jc.one_call(cli, e5.JUDGE, prompt, **e1.CALL)
        if out['ok'] and e1.malformed(out['text']):
            out['malformed'] = True
        tries.append({k: out.get(k) for k in ('ok', 'finish_reason', 'completion_tokens', 'billed_usd', 'seconds',
                                              'error', 'serving_provider', 'malformed')})
        if out['ok'] and out['text'].strip() and not out.get('malformed'):
            break
    out['attempts'] = tries
    out['billed_usd'] = sum(t.get('billed_usd') or 0.0 for t in tries)
    out['billed_completion_tokens'] = sum(t.get('completion_tokens') or 0 for t in tries)
    if out['ok'] and not out['text'].strip():
        out.update(ok=False, error=f'empty reply after {e1.ATTEMPTS} attempts')
    elif out['ok'] and out.get('malformed'):
        out.update(ok=False, error=f'malformed reply after {e1.ATTEMPTS} attempts')
    return out


def summarise(milestones: list[dict], reached: list[bool], sent: bool, reply: dict | None) -> dict:
    """Each milestone's source and the E5 scores for one trace."""
    verdicts = e5.parse(reply.get('text')) if reply and reply.get('ok') else {}
    sources = ['e3' if r else (verdicts.get(m['id'], 'UNJUDGED') if sent else 'UNJUDGED')
               for m, r in zip(milestones, reached)]
    n = len(milestones)
    c = collections.Counter(sources)
    return {'sources': sources, 'milestones_required': n, 'by_e3': c['e3'], 'reached': c['REACHED'],
            'not_needed': c['NOT_NEEDED'], 'missing': c['MISSING'], 'unjudged': c['UNJUDGED'],
            'e5_strict': (c['e3'] + c['REACHED']) / n if n else None,
            'e5_lenient': (c['e3'] + c['REACHED'] + c['NOT_NEEDED']) / n if n else None,
            'judged_fraction': (n - c['e3']) / n if n else None,
            'sent': sent, 'reply_ok': bool(reply and reply.get('ok')) if sent else None}


# ------------------------------------------------------------------ the full run

def roster(variant: str) -> list[str]:
    from full_run_28092026.analyze import ROSTER
    return [k for k in ROSTER if (score.SCORES / variant / f'{k}.jsonl').exists()]


def jobs_for(variant: str, key: str, items: dict, ms_all: dict) -> list[tuple]:
    """(row, milestones, job or None) for every scored trace of one model."""
    rows = {r['item_id']: r for r in map(json.loads, (score.SCORES / variant / f'{key}.jsonl')
                                              .read_text(encoding='utf-8').splitlines())}
    texts = score.texts_matching(variant, key, rows.values())   # refuses if a trace no longer matches its store row
    out = []
    for i, row in rows.items():
        ms = ms_all[i]
        reached = row['e3']['reached']
        j = job(items[i]['question'], texts[i], ms, reached) if row['status'] == 'answered' else None
        out.append((row, ms, reached, j))
    return out


def write(variant: str, keys: list[str], store: jc.Store) -> dict:
    """The E5 rows from the stored replies. Free."""
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    out_dir = score.SCORES / variant / 'e5'
    out_dir.mkdir(parents=True, exist_ok=True)
    totals = {}
    for key in keys:
        n_sent = n_missing = 0
        with open(out_dir / f'{key}.jsonl', 'w', encoding='utf-8', newline='\n') as fh:
            for row, ms, reached, j in jobs_for(variant, key, items, ms_all):
                reply = store.get(j[2]) if j else None
                s = summarise(ms, reached, j is not None, reply)
                n_sent += j is not None
                n_missing += j is not None and reply is None
                meta = {k: reply.get(k) for k in ('serving_provider', 'served_model', 'prompt_tokens',
                                                  'billed_completion_tokens', 'billed_usd')} if reply else None
                fh.write(json.dumps({'item_id': row['item_id'], 'template_id': row['template_id'],
                                     'status': row['status'], **s, 'reply': meta}) + '\n')
        totals[key] = {'sent': n_sent, 'without_reply': n_missing}
    (out_dir / 'CONFIG.json').write_text(json.dumps({'judge': e5.JUDGE, 'call': e1.CALL, 'attempts': e1.ATTEMPTS,
                                                     'prompt_sha256': e5.config()['prompt_sha256'],
                                                     'e5_sha256': e5.config()['e5_sha256'],
                                                     'models': totals, 'reply_store': store.summary(),
                                                     **jc.provenance(score.SCORES / variant / 'CONFIG.json'),
                                                     'written_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())},
                                                    indent=1) + '\n', encoding='utf-8')
    return totals


def status(variant: str, keys: list[str]) -> int:
    """What the reply store holds and what the stage has written. Free."""
    st = jc.Store(REPLIES)
    print(f'reply store {REPLIES.name}: {st.summary()}')
    out_dir = score.SCORES / variant / 'e5'
    for key in keys:
        p = out_dir / f'{key}.jsonl'
        if not p.exists():
            print(f'  {key:24s} no rows written')
            continue
        rows = [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines()]
        sent = [r for r in rows if r['sent']]
        print(f'  {key:24s} rows {len(rows)}, calls {len(sent)}, without a reply '
              f'{sum(r["reply_ok"] is False for r in sent)}, milestones unjudged {sum(r["unjudged"] for r in sent)}')
    return 0


def prices() -> dict:
    """MiMo's endpoints from OpenRouter's public list, per million tokens. No key is sent."""
    req = urllib.request.Request(f'https://openrouter.ai/api/v1/models/{e5.JUDGE}/endpoints',
                                 headers={'User-Agent': 'engtrace-full-run'})
    eps = json.load(urllib.request.urlopen(req, timeout=60))['data']['endpoints']
    return {e['provider_name']: (float(e['pricing']['prompt']) * 1e6, float(e['pricing']['completion']) * 1e6)
            for e in eps}


def pilot_basis() -> dict:
    """Tokens per prompt character, completion tokens and seconds per call, from the pilot's replies."""
    rep = replay()
    sent = [t for t in rep['traces'] if t['job']]
    store = jc.Store(PILOT_REPLIES)
    chars = toks = 0
    outs, secs = [], []
    for t in sent:
        r = store.get(t['job'][2])
        if r and r.get('prompt_tokens'):
            chars += len(t['job'][1])
            toks += r['prompt_tokens']
            outs.append(r.get('billed_completion_tokens') or r.get('completion_tokens') or 0)
            secs.append(r.get('seconds') or 0)
    return {'tokens_per_char': toks / chars, 'completion_per_call': statistics.mean(outs),
            'seconds_per_call': statistics.median(secs), 'calls': len(outs)}


def dry_run(variant: str, keys: list[str], workers: int) -> int:
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    basis = pilot_basis()
    pr = prices()
    total_calls, total_chars = 0, 0
    print(f'E5 on {variant}: basis from the pilot\'s {basis["calls"]} replies: {basis["tokens_per_char"]:.3f} '
          f'tokens per prompt character, {basis["completion_per_call"]:.0f} completion tokens and '
          f'{basis["seconds_per_call"]:.0f} s per call')
    print(f'{"model":24s} {"traces":>7s} {"calls":>6s} {"milestones sent":>16s} {"in tokens":>10s}')
    for key in keys:
        js = jobs_for(variant, key, items, ms_all)
        sent = [j for *_, j in js if j]
        chars = sum(len(j[1]) for j in sent)
        total_calls += len(sent)
        total_chars += chars
        print(f'{key:24s} {len(js):7d} {len(sent):6d} {sum(len(j[0]) for j in sent):16d} '
              f'{chars * basis["tokens_per_char"]:10.0f}')
    tin = total_chars * basis['tokens_per_char']
    tout = total_calls * basis['completion_per_call']
    cost = lambda p: (tin * p[0] + tout * p[1]) / 1e6
    at = pr.get('Xiaomi') or min(pr.values())
    lo, hi = min(pr.values(), key=cost), max(pr.values(), key=cost)
    print(f'{"TOTAL":24s} {"":7s} {total_calls:6d} {"":16s} {tin:10.0f}')
    print(f'\nestimate: ${cost(at):.2f} at Xiaomi\'s own endpoint (${at[0]:.3f} in / ${at[1]:.3f} out per million); '
          f'${cost(lo):.2f} to ${cost(hi):.2f} across the {len(pr)} endpoints OpenRouter lists')
    print(f'time: about {total_calls * basis["seconds_per_call"] / workers / 3600:.1f} hours at {workers} workers. '
          'Nothing was called.')
    return 0


def paid(variant: str, keys: list[str], workers: int, max_usd: float) -> int:
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    store = jc.Store(REPLIES)
    jobs = jc.interleave({key: {j[2]: j[1] for *_, j in jobs_for(variant, key, items, ms_all) if j} for key in keys})
    cli = jc.client()
    summary = jc.run(jobs, lambda p: fetch(p, cli), store, workers, max_usd, 'E5')
    print(write(variant, keys, store))
    jc.finish(summary)
    return 0


# ------------------------------------------------------------------ the pilot replay

def replay(current_e3: bool = False) -> dict:
    """The pilot's 300 E5 traces through job() and summarise(), their milestones as the pilot's rows hold
    them, E3 recomputed (with the code the pilot ran, or the current code), replies from its store."""
    from full_run_28092026 import parser_fix
    old_a, old_m = parser_fix.pre_fix_modules()
    import e3_milestones as e3
    import milestones as cur_m
    M = cur_m if current_e3 else old_m
    items = {json.loads(l)['item_id']: json.loads(l)
             for l in (PILOT / 'slice' / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()}
    texts = {}
    for f in (PILOT / 'traces').glob('*.jsonl'):
        for line in f.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            if r.get('ok'):
                texts[(r['model_key'], r['item_id'])] = r['text']
    store = jc.Store(PILOT_REPLIES)
    out = []
    for f in sorted(PILOT_ROWS.glob('*.jsonl')):
        for line in f.read_text(encoding='utf-8').splitlines():
            row = json.loads(line)
            text = texts[(row['model_key'], row['item_id'])]
            ms = [{'id': m['id'], 'value': m['value']} for m in row['meta']['milestones']]
            nums = M.numbers(text)
            reached = [M.scaled_match(m['value'], nums, e3.STEP_TOL) is not None for m in ms]
            j = job(items[row['item_id']]['question'], text, ms, reached)
            reply = store.get(j[2]) if j else None
            s = summarise(ms, reached, j is not None, reply)
            out.append({'model_key': row['model_key'], 'item_id': row['item_id'], 'job': j, 'summary': s,
                        'pilot': row, 'reached_as_pilot': reached == [m['reached'] for m in row['meta']['milestones']],
                        'reply_text_as_pilot': (not j) or (reply is not None and row['calls']
                                                           and reply.get('text') == row['calls'][0]['text'])})
    return {'traces': out}


def validate() -> int:
    rep = replay()
    tr = rep['traces']
    sent = [t for t in tr if t['job']]
    src_same = sum(t['summary']['sources'] == [m['source'] for m in t['pilot']['meta']['milestones']] for t in tr)
    score_same = sum(abs((t['summary']['e5_strict'] or 0.0) - t['pilot']['scores']['e5_strict']) < 1e-12 for t in tr)
    c = collections.Counter(s for t in tr for s in t['summary']['sources'])
    per = collections.defaultdict(list)
    for t in tr:
        per[t['model_key']].append(t['summary'])
    cur = replay(current_e3=True)['traces']
    changed = [t for t, u in zip(tr, cur) if (t['job'] or [None, None, None])[2] != (u['job'] or [None, None, None])[2]]
    published = {'gpt-5': (0.962, 0.981), 'claude-opus-4.7': (0.946, 0.983), 'gemini-3.1-pro': (0.929, 0.939),
                 'deepseek-r1': (0.923, 0.961), 'llama-3.1-70b': (0.321, 0.341)}
    L = ['# E5 over the full run: the stage validated on the pilot (D-142)', '',
         'Generated by `judge.py --validate`: the pilot\'s 300 E5 traces through the functions the full run uses '
         '(`job`, `summarise`), their milestones as the pilot\'s rows hold them, E3 recomputed with the code the '
         'pilot ran, and every reply read from the pilot\'s store. Nothing is called.', '',
         '| check | result |', '|---|---|',
         f'| traces | {len(tr)} |',
         f'| E3\'s reached flags equal the pilot\'s | {sum(t["reached_as_pilot"] for t in tr)} of {len(tr)} |',
         f'| traces sent to the judge | {len(sent)} (published: 142) |',
         f'| prompts found in the pilot\'s store, by the same key | {sum(t["summary"]["reply_ok"] is True for t in sent)} of {len(sent)} |',
         f'| reply text equal to the one in the pilot\'s row | {sum(t["reply_text_as_pilot"] for t in tr)} of {len(tr)} |',
         f'| every milestone\'s source equal to the pilot\'s | {src_same} of {len(tr)} traces |',
         f'| E5-strict equal to the pilot\'s | {score_same} of {len(tr)} traces |',
         f'| milestones: by E3, judged | {c["e3"]} of {sum(c.values())}, {sum(c.values()) - c["e3"]} '
         f'(published: 951 of 1,245, 294) |',
         f'| judged: REACHED, NOT_NEEDED, MISSING, UNJUDGED | {c["REACHED"]}, {c["NOT_NEEDED"]}, {c["MISSING"]}, '
         f'{c["UNJUDGED"]} (published: 73, 36, 185) |', '',
         '| model | E5 strict | published | E5 lenient | published |', '|---|---:|---:|---:|---:|']
    ok = src_same == len(tr) and score_same == len(tr)
    for m, (ps, pl) in published.items():
        s = statistics.mean(x['e5_strict'] or 0.0 for x in per[m])
        l = statistics.mean(x['e5_lenient'] or 0.0 for x in per[m])
        ok &= abs(s - ps) < 5e-4 and abs(l - pl) < 5e-4
        L.append(f'| {m} | {s:.3f} | {ps:.3f} | {l:.3f} | {pl:.3f} |')
    L += ['', f'With E3 as it reads numbers now (D-137), {len(changed)} of the {len(tr)} traces would send a different '
          'prompt, which the pilot\'s store does not hold; the full run sends its own prompts in any case.', '',
          'The stage reproduces the pilot.' if ok else 'A figure differs from the pilot: see the table.']
    (HERE / 'E5_VALIDATION.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--score', action='store_true')
    ap.add_argument('--status', action='store_true', help='FREE: the reply store and the rows written')
    ap.add_argument('--variant', default='main')
    ap.add_argument('--model', action='append')
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--max-usd', type=float, default=100.0,
                    help='the cumulative cap for the stage: what the reply store already records counts toward it')
    ap.add_argument('--yes', action='store_true', help='required for the run, which bills')
    a = ap.parse_args()
    if a.validate:
        return validate()
    keys = a.model or roster(a.variant)
    if a.dry_run:
        return dry_run(a.variant, keys, a.workers)
    if a.status:
        return status(a.variant, keys)
    if a.score:
        print(write(a.variant, keys, jc.Store(REPLIES)))
        return 0
    if not a.yes:
        raise SystemExit('E5 bills: re-run with --yes once the spend is approved (see --dry-run for the estimate)')
    return paid(a.variant, keys, a.workers, a.max_usd)


if __name__ == '__main__':
    raise SystemExit(main())
