"""The step router over the full run (D-110, D-112, D-113): the digit rule's flags, and MiMo-V2.5-Pro on
every other step, a trace's steps in one batched call on the framework's own Tribunal prompt (D-142).

    python -m full_run_28092026.router --dry-run [--variant main]      # FREE: the calls, the cost, the time
    python -m full_run_28092026.router --validate                      # FREE: the pilot's 300 labelled traces, from its stored replies
    python -m full_run_28092026.router --yes [--model KEY] [--max-usd 150] [--workers 8]   # BILLS
    python -m full_run_28092026.router --score [--variant main]        # FREE: the router store from the stored replies

WHAT IS SENT follows router_batched.py, which measured the design on the pilot:
  - Rule C (D-112). Every non-empty step the digit rule, as E4 ships it, does not flag goes to the
    judge; a flagged step is already an error and is not sent. The flags come from score.py's store,
    whose steps are the framework's non-empty steps in order (e2_prm.steps_of).
  - The prompt is the framework's own Tribunal prompt. It is evaluated from the f-string in
    `EngTraceFramework._tier2_tribunal_batch`, read from evaluation/engtrace_evaluation_framework.py at
    run time, so it cannot drift from the framework. It holds the question, the gold's steps and the
    trace's steps as the framework splits them, and the raw indices under review.
  - A trace with no step to send (none written, or every one flagged) makes no call.

THE CALL is router_batched.call's:
  - 8,192 tokens at the provider's defaults, at OpenRouter's default routing;
  - a reply the framework's parser cannot read is asked again, up to three attempts;
  - a reply cut at 8,192 is asked again at 16,384;
  - a call gets 900 s of wall clock.
The key is the SHA-256 of (model, ceilings, prompt), so a prompt already answered is never bought
again.

WHAT IS WRITTEN. scores/<variant>/router/<model>.jsonl, one row per trace:
  - the steps the digit rule flags, the steps sent, and the judge's category for each step sent
    (None when the reply does not name it: unjudged);
  - the router's flags: the digit rule's, plus every step whose category contains "error", as the
    framework and router_batched.py read it ("Other" is recorded but is not a flag);
  - the reply's provider, tokens and billed cost.
The replies are kept in scores/_judge/router_replies.jsonl; both files are gitignored.

THE ESTIMATE counts the calls, prompt lengths and steps under review from the store. It takes tokens
per character, completion tokens per step under review, and seconds per call from the pilot's
labelled replies, and prices from OpenRouter's public endpoint list, as judge.py does.
"""
from __future__ import annotations

import argparse
import ast
import collections
import hashlib
import io
import json
import statistics
import sys
import time
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
    import judge_probe as JP     # noqa: E402  parse(): the framework's own reply extraction, copied there
from engineering_parser import extract_steps  # noqa: E402

from full_run_28092026 import judge_calls as jc, score  # noqa: E402

JUDGE = 'xiaomi/mimo-v2.5-pro'
CEILING, CEILING_RETRY, ATTEMPTS = 8192, 16384, 3
FRAMEWORK = REPO / 'evaluation' / 'engtrace_evaluation_framework.py'
REPLIES = score.SCORES / '_judge' / 'router_replies.jsonl'
PILOT_PROBE = PILOT / 'scores' / 'router_batched' / 'probe_set.json'
PILOT_REPLIES = PILOT / 'scores' / 'router_batched' / 'replies.jsonl'
LABELS = PILOT / 'experts_filled_labels' / 'version_2' / 'labels'


# ------------------------------------------------------------------ the prompt

def _tribunal_fstring():
    """The f-string the framework's _tier2_tribunal_batch builds its prompt from, compiled."""
    tree = ast.parse(FRAMEWORK.read_text(encoding='utf-8'))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == '_tier2_tribunal_batch':
            for stmt in node.body:
                if isinstance(stmt, ast.Assign) and getattr(stmt.targets[0], 'id', None) == 'prompt':
                    return compile(ast.Expression(stmt.value), str(FRAMEWORK), 'eval')
    raise SystemExit(f'no prompt assignment in _tier2_tribunal_batch of {FRAMEWORK}')


TRIBUNAL = _tribunal_fstring()
TRIBUNAL_SHA = hashlib.sha256(ast.dump(ast.parse(FRAMEWORK.read_text(encoding='utf-8'))).encode()).hexdigest()


def tribunal_prompt(question: str, gt_steps: list, pred_steps: list, indices: list[int]) -> str:
    return eval(TRIBUNAL, {'json': json}, {'question': question, 'gt_steps': gt_steps,       # noqa: S307
                                           'pred_steps': pred_steps, 'mismatch_indices': list(indices)})


def split(text: str) -> tuple[list[str], list[int]]:
    """(the framework's raw steps, the raw index of each non-empty one): router_batched.split."""
    raw = extract_steps(text)[0]
    return raw, [i for i, s in enumerate(raw) if s.strip()]


def job(question: str, solution: str, text: str, digit_flags: list[int]):
    """(steps under review, the raw index map, prompt, key), or None when there is nothing to send.
    digit_flags: the store's digit-rule flag count per non-empty step."""
    raw, nonempty = split(text)
    if len(nonempty) != len(digit_flags):
        raise ValueError('the store and the framework split the trace differently')
    review = [j for j, f in enumerate(digit_flags) if not f]          # rule C
    if not review:
        return None
    prompt = tribunal_prompt(question, extract_steps(solution)[0], raw, [nonempty[j] for j in review])
    key = hashlib.sha256(json.dumps([JUDGE, {'max_tokens': CEILING, 'retry': CEILING_RETRY}, prompt],
                                    sort_keys=True).encode()).hexdigest()
    return review, nonempty, prompt, key


def fetch(prompt: str, cli=None) -> dict:
    """router_batched.one(): up to three attempts, a reply cut at the ceiling asked again with more room."""
    cli = cli or jc.client(timeout=300.0)
    tries, ceiling, res = [], CEILING, {}
    for attempt in range(1, ATTEMPTS + 1):
        res = jc.one_call(cli, JUDGE, prompt, max_tokens=ceiling)
        tries.append({k: res.get(k) for k in ('ok', 'finish_reason', 'prompt_tokens', 'completion_tokens',
                                              'billed_usd', 'seconds', 'error', 'serving_provider')} | {'max_tokens': ceiling})
        if res.get('ok') and JP.parse(res.get('text')):
            break
        if res.get('ok') and res.get('finish_reason') == 'length':
            ceiling = CEILING_RETRY
        time.sleep(2 * attempt)
    ok = bool(res.get('ok') and JP.parse(res.get('text')))
    return dict(res, ok=ok, error=None if ok else (res.get('error') or 'no readable reply after 3 attempts'),
                attempts=tries, billed_usd=sum(t.get('billed_usd') or 0.0 for t in tries),
                billed_completion_tokens=sum(t.get('completion_tokens') or 0 for t in tries))


def verdicts(text: str | None, review: list[int], nonempty: list[int]) -> dict[int, str | None]:
    """{step index: the judge's category, or None if the reply does not name it}: router_batched.verdicts."""
    back = {raw: j for j, raw in enumerate(nonempty)}
    got = {j: None for j in review}
    for x in JP.parse(text) or []:
        if isinstance(x, dict) and back.get(x.get('step_index')) in got:
            got[back[x['step_index']]] = str(x.get('category') or '')
    return got


def flag(cat) -> bool:
    return cat is not None and 'error' in cat.lower()


def summarise(digit_flags: list[int], sent: list[int] | None, got: dict | None) -> dict:
    """A call with no reply leaves every step it sent unjudged (D-148): the first version recorded
    0 unjudged and no flags for it, so a failed call looked like a clean trace."""
    digit = [j for j, f in enumerate(digit_flags) if f]
    got = got if got is not None else ({j: None for j in sent} if sent else {})
    judge = [j for j, c in got.items() if flag(c)]
    return {'steps': len(digit_flags), 'digit_flagged': digit, 'sent': sent or [],
            'categories': {str(j): c for j, c in got.items()}, 'judge_flagged': sorted(judge),
            'router_flagged': sorted(set(digit) | set(judge)),
            'unjudged': sum(c is None for c in got.values()) if sent else 0}


# ------------------------------------------------------------------ the full run

def roster(variant: str) -> list[str]:
    from full_run_28092026.analyze import ROSTER
    return [k for k in ROSTER if (score.SCORES / variant / f'{k}.jsonl').exists()]


def jobs_for(variant: str, key: str, items: dict) -> list[tuple]:
    rows = [json.loads(l) for l in (score.SCORES / variant / f'{key}.jsonl').read_text(encoding='utf-8').splitlines()]
    texts = score.texts_matching(variant, key, rows)      # refuses if a trace no longer matches its store row
    out = []
    for row in rows:
        flags = [s['digit_flags'] for s in row['steps']]
        j = None
        if row['status'] == 'answered' and flags:
            it = items[row['item_id']]
            j = job(it['question'], it['solution'], texts[row['item_id']], flags)
        out.append((row, flags, j))
    return out


def status(variant: str, keys: list[str]) -> int:
    """What the reply store holds and what the stage has written. Free."""
    st = jc.Store(REPLIES)
    print(f'reply store {REPLIES.name}: {st.summary()}')
    out_dir = score.SCORES / variant / 'router'
    for key in keys:
        p = out_dir / f'{key}.jsonl'
        if not p.exists():
            print(f'  {key:24s} no rows written')
            continue
        rows = [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines()]
        sent = [r for r in rows if r['sent']]
        print(f'  {key:24s} rows {len(rows)}, calls {len(sent)}, without a reply '
              f'{sum(r["reply_ok"] is False for r in sent)}, steps unjudged {sum(r["unjudged"] for r in sent)}')
    return 0


def write(variant: str, keys: list[str], store: jc.Store) -> dict:
    items = score.pool_items()
    out_dir = score.SCORES / variant / 'router'
    out_dir.mkdir(parents=True, exist_ok=True)
    totals = {}
    for key in keys:
        n_sent = n_missing = 0
        with open(out_dir / f'{key}.jsonl', 'w', encoding='utf-8', newline='\n') as fh:
            for row, flags, j in jobs_for(variant, key, items):
                reply = store.get(j[3]) if j else None
                got = verdicts(reply.get('text'), j[0], j[1]) if reply else None
                n_sent += j is not None
                n_missing += j is not None and reply is None
                s = summarise(flags, j[0] if j else None, got)
                meta = {k: reply.get(k) for k in ('serving_provider', 'served_model', 'prompt_tokens',
                                                  'billed_completion_tokens', 'billed_usd')} if reply else None
                fh.write(json.dumps({'item_id': row['item_id'], 'template_id': row['template_id'],
                                     'status': row['status'], **s, 'reply_ok': bool(reply) if j else None,
                                     'reply': meta}) + '\n')
        totals[key] = {'sent': n_sent, 'without_reply': n_missing}
    (out_dir / 'CONFIG.json').write_text(json.dumps({'judge': JUDGE, 'ceiling': CEILING, 'ceiling_retry': CEILING_RETRY,
                                                     'attempts': ATTEMPTS, 'framework_ast_sha256': TRIBUNAL_SHA,
                                                     'models': totals, 'reply_store': store.summary(),
                                                     **jc.provenance(score.SCORES / variant / 'CONFIG.json'),
                                                     'written_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())},
                                                    indent=1) + '\n', encoding='utf-8')
    return totals


def pilot_basis() -> dict:
    """Tokens per prompt character, completion tokens per step under review, seconds per call: the
    pilot's labelled replies."""
    probe = {p['key']: p for p in json.load(open(PILOT_PROBE, encoding='utf-8')) if p['check'] == 'labelled'}
    rows = [json.loads(l) for l in open(PILOT_REPLIES, encoding='utf-8') if l.strip()]
    ok = {r['key']: r for r in rows if r.get('check') == 'labelled' and r.get('ok')}
    chars = sum(len(probe[k]['prompt']) for k in ok)
    steps = sum(len(probe[k]['review']) for k in ok)
    return {'tokens_per_char': sum(r['in'] for r in ok.values()) / chars,
            'completion_per_step': sum(r['out'] for r in ok.values()) / steps,
            'completion_per_call': statistics.mean(r['out'] for r in ok.values()),
            'seconds_per_call': statistics.median(r['seconds'] for r in ok.values()), 'calls': len(ok)}


def dry_run(variant: str, keys: list[str], workers: int) -> int:
    from full_run_28092026.judge import prices
    items = score.pool_items()
    b = pilot_basis()
    pr = prices()
    calls = chars = steps = 0
    print(f'router on {variant}: basis from the pilot\'s {b["calls"]} labelled replies: {b["tokens_per_char"]:.3f} '
          f'tokens per prompt character, {b["completion_per_step"]:.0f} completion tokens per step under review '
          f'({b["completion_per_call"]:.0f} per call), {b["seconds_per_call"]:.0f} s per call')
    print(f'{"model":24s} {"traces":>7s} {"calls":>6s} {"steps sent":>11s} {"per call":>9s} {"in tokens":>10s}')
    for key in keys:
        js = [j for *_, j in jobs_for(variant, key, items)]
        sent = [j for j in js if j]
        c, s = sum(len(j[2]) for j in sent), sum(len(j[0]) for j in sent)
        calls, chars, steps = calls + len(sent), chars + c, steps + s
        print(f'{key:24s} {len(js):7d} {len(sent):6d} {s:11d} {s / max(1, len(sent)):9.1f} '
              f'{c * b["tokens_per_char"]:10.0f}')
    tin = chars * b['tokens_per_char']
    cost = lambda p, tout: (tin * p[0] + tout * p[1]) / 1e6
    at = pr.get('Xiaomi') or min(pr.values())
    per_step, per_call = steps * b['completion_per_step'], calls * b['completion_per_call']
    print(f'{"TOTAL":24s} {"":7s} {calls:6d} {steps:11d} {steps / max(1, calls):9.1f} {tin:10.0f}')
    print(f'\nestimate at Xiaomi\'s own endpoint: ${cost(at, per_step):.2f} with output scaled by the steps under '
          f'review, ${cost(at, per_call):.2f} at the pilot\'s output per call; across the {len(pr)} endpoints '
          f'${min(cost(p, per_step) for p in pr.values()):.2f} to ${max(cost(p, per_step) for p in pr.values()):.2f} '
          '(scaled by steps)')
    print(f'time: about {calls * b["seconds_per_call"] / workers / 3600:.1f} hours at {workers} workers '
          f'(the pilot ran 4, because the endpoint refused more then). Nothing was called.')
    return 0


def paid(variant: str, keys: list[str], workers: int, max_usd: float) -> int:
    items = score.pool_items()
    store = jc.Store(REPLIES)
    jobs = jc.interleave({key: {j[3]: j[2] for *_, j in jobs_for(variant, key, items) if j} for key in keys})
    cli = jc.client(timeout=300.0)
    summary = jc.run(jobs, lambda p: fetch(p, cli), store, workers, max_usd, 'router')
    print(write(variant, keys, store))
    jc.finish(summary)
    return 0


# ------------------------------------------------------------------ the pilot replay

def validate() -> int:
    """The pilot's 300 labelled traces through job(), verdicts() and summarise(), the digit flags as
    score.py computes them, the replies from router_batched.py's store, scored against the experts."""
    import score_against_labels as S
    with redirect_stdout(io.StringIO()):
        import x1_analysis as X
    probe = {p['key']: p for p in json.load(open(PILOT_PROBE, encoding='utf-8')) if p['check'] == 'labelled'}
    latest = {}
    for r in (json.loads(l) for l in open(PILOT_REPLIES, encoding='utf-8') if l.strip()):
        if r.get('check') == 'labelled' and (r.get('ok') or r['key'] not in latest):
            latest[r['key']] = r
    items = {json.loads(l)['item_id']: json.loads(l)
             for l in (PILOT / 'slice' / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()}
    texts = {}
    for f in (PILOT / 'traces').glob('*.jsonl'):
        for line in f.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            if r.get('ok'):
                texts[(r['model_key'], r['item_id'])] = r['text']
    truth, _ = S.build_truth(S.read_labels(str(LABELS)), S.read_consensus(str(LABELS)))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(Path(S.TASKS) / 'keyfile.jsonl', encoding='utf-8'))}
    same_prompt = same_review = built = 0
    steps_rows = []
    for code, t in truth.items():
        mk, iid = keyfile[code]
        text = texts.get((mk, iid))
        if not text:
            continue
        flags = [s['digit_flags'] for s in score.score_steps(text)]
        if len(flags) != len(t['steps']):
            continue
        it = items[iid]
        j = job(it['question'], it['solution'], text, flags)
        if j and code in probe:
            built += 1
            same_prompt += j[2] == probe[code]['prompt']
            same_review += j[0] == probe[code]['review'] and j[1] == probe[code]['nonempty']
        rep = latest.get(code)
        if j and code in probe and not (rep and rep.get('ok')):
            continue                                 # sent, but no reply: left out, as router_batched.py does
        got = verdicts(rep.get('text'), j[0], j[1]) if j and rep else None
        s = summarise(flags, j[0] if j else None, got)
        for k, lab in enumerate(t['steps']):
            steps_rows.append((code, k, lab['label'], k in s['digit_flagged'], k in s['judge_flagged'],
                               str(k) in s['categories'] and s['categories'][str(k)] is None))
    res = {}
    for scope, pick in (('all', lambda c: True), ('correct_answer', lambda c: truth[c]['final_answer'] == 'correct')):
        c = collections.Counter()
        for code, k, lab, rule, judge_, unj in steps_rows:
            if not pick(code) or lab not in ('incorrect', 'correct', 'alternative_correct'):
                continue
            bad = lab == 'incorrect'
            for who, hit in (('rule', rule), ('router', rule or judge_)):
                c[(who, 'tp' if bad and hit else 'fn' if bad else 'fp' if hit else 'tn')] += 1
            c['unjudged'] += unj
        res[scope] = {who: dict(zip(('p', 'r', 'f'), S.prf(c[(who, 'tp')], c[(who, 'fp')], c[(who, 'fn')])),
                                tp=c[(who, 'tp')], fp=c[(who, 'fp')], fn=c[(who, 'fn')]) for who in ('rule', 'router')}
        res[scope]['unjudged'] = c['unjudged']
    flags_per = collections.Counter()
    covered = set()
    for code, k, lab, rule, judge_, unj in steps_rows:
        covered.add(code)
        flags_per[code] += rule or judge_
    right = [cd for cd in truth if truth[cd]['final_answer'] == 'correct' and cd in covered]
    y = {cd: any(s['label'] == 'incorrect' for s in truth[cd]['steps']) for cd in right}
    auc = X.auroc([(-flags_per[cd], not y[cd]) for cd in right])
    lo, hi = X.boot(right, lambda ks: X.auroc([(-flags_per[cd], not y[cd]) for cd in ks]))
    published = {('all', 'rule'): (0.825, 0.255, 0.390), ('all', 'router'): (0.707, 0.603, 0.651),
                 ('correct_answer', 'rule'): (0.750, 0.320, 0.449), ('correct_answer', 'router'): (0.703, 0.360, 0.476)}
    L = ['# The step router over the full run: the stage validated on the pilot (D-142)', '',
         'Generated by `router.py --validate`: the pilot\'s labelled traces through the functions the full run uses '
         '(`job`, `verdicts`, `summarise`), the digit flags as `score.py` computes them, the replies read from '
         '`router_batched.py`\'s store, scored against the experts\' step labels as `router_batched.py` scores them. '
         'Nothing is called.', '',
         f'Prompts built: {built}. Byte-identical to the prompt the pilot sent: {same_prompt} of {built}; the same '
         f'steps under review and the same index map: {same_review} of {built}. Steps sent but left unjudged by the '
         f"reply: {res['all']['unjudged']} (published: 8 of 1,971).", '',
         '| steps the experts call incorrect | tp / fp / fn | precision | recall | F1 | published |',
         '|---|---:|---:|---:|---:|---|']
    ok = same_prompt == built and same_review == built
    names = {('all', 'rule'): 'all traces, the digit rule alone', ('all', 'router'): 'all traces, the router',
             ('correct_answer', 'rule'): 'inside correct-answer traces, the digit rule alone',
             ('correct_answer', 'router'): 'inside correct-answer traces, the router'}
    for (scope, who), pub in published.items():
        w = res[scope][who]
        ok &= all(abs(a - b) < 5e-4 for a, b in zip((w['p'], w['r'], w['f']), pub))
        L.append(f"| {names[(scope, who)]} | {w['tp']} / {w['fp']} / {w['fn']} | {w['p']:.3f} | {w['r']:.3f} | "
                 f"{w['f']:.3f} | {pub[0]:.3f} / {pub[1]:.3f} / {pub[2]:.3f} |")
    ok &= abs(auc - 0.675) < 5e-4
    L += ['', f'Correct-answer traces ranked by flagged steps: {len(right)} traces, AUROC {auc:.3f} ({lo:.3f} to '
          f'{hi:.3f}, traces resampled); published 0.675 (0.622 to 0.730) on 228.', '',
          'The stage reproduces the pilot.' if ok else 'A figure differs from the pilot: see the table.']
    (HERE / 'ROUTER_VALIDATION.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
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
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--max-usd', type=float, default=150.0,
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
        raise SystemExit('the router bills: re-run with --yes once the spend is approved (see --dry-run)')
    return paid(a.variant, keys, a.workers, a.max_usd)


if __name__ == '__main__':
    raise SystemExit(main())
