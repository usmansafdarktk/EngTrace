"""The paraphrase arm of ANALYSIS_PLAN Q5: one paraphrase for each of the 450 items of
subsamples.PARAPHRASE, written by a model from a family on neither the roster nor the judge's side,
and checked by script before an expert sees it (D-141).

    python -m full_run_28092026.paraphrase --dry-run    # FREE: the selection, the writer's endpoint, the estimate
    python -m full_run_28092026.paraphrase --yes        # BILLS: writes the paraphrases; resumable
    python -m full_run_28092026.paraphrase --yes --limit 40   # BILLS: the first 40 unresolved items only, a pilot of the prompt
    python -m full_run_28092026.paraphrase --check      # FREE: re-checks every attempt, rewrites the manifest
    python -m full_run_28092026.paraphrase --status     # FREE

THE WRITER is Mistral Large 3 (`mistralai/mistral-large-2512`). The roster's families are OpenAI, Google,
DeepSeek, Qwen, Zhipu, Meta, Moonshot and Anthropic (D-110), and the judge is Xiaomi's MiMo (D-105);
Mistral is none of them. Mistral alone serves the model and declares no quantization, so it is routed
as the closed-weight models are: the cheapest endpoint, which is Mistral's own.

THE PROMPT (PROMPT below, hashed into every row) asks for a rewrite in new words and sentence
structure that keeps every number, unit, symbol, variable, formula and technical term as written and
the parts in the same order, adds and removes nothing, does not hint at the method or the answer, and
returns only the problem. The rules grew from what the writer did on 30 September: the first prompt
(6091f248) lost 49 of 62 attempts to prettified notation (LaTeX, $ delimiters, Unicode superscripts and
minus signs, · for *); the second (1ae64428), which forbade that, lost 43 of 61 to acronyms spelt out
(PFR, CSTR), reaction orders turned into words ("order-2.1" as "second-order", a real error), = written as
"equals", $ delimiters dropped where the original had them, near-copies and length. Rules 1 to 3 and 6 now
name each of these. The checks were not loosened: the two arms must differ in wording only. A new prompt
is tried on the first items with --limit before the rest are written. Mistral throttles OpenRouter's shared capacity upstream ("temporarily
rate-limited"): on 30 September most calls were refused at eight workers and at two alike, and the
throughput scaled with the workers, so the throttle is per call, not a cap we saturate. A refused call
bills nothing; `write_one` therefore retries each call through RETRY_SLEEPS before recording a service
failure, and eight workers stay the default.

THE CHECKS. An attempt passes only if all hold; otherwise the writer is asked again, up to three
attempts per item. An item with no passing attempt has no paraphrase and leaves both arms of Q5.
  numbers   the same numbers as the original, as a multiset (milestones.numbers)
  tokens    every technical token of the original - a unit, a symbol, a name with a digit, subscript,
            caret or slash, an acronym, a non-ASCII symbol - at least as often
  parts     the same part labels, (a), (b), (i), Part 1 ..., in the same order
  copy      word-level similarity to the original at most COPY (a near-copy tests nothing)
  length    between LENGTH[0] and LENGTH[1] times the original's length
  clean     no preamble, comment or fence around the problem
Then one own-branch expert confirms each is the same problem with the same answer (paraphrase_kit.py);
a pair the expert rejects leaves both arms too.

WHAT IS KEPT. paraphrase/attempts.jsonl, every attempt with its text, and paraphrase/pool.jsonl, the
passing paraphrase per item, stay local and gitignored like the pool. paraphrase/manifest.jsonl is
committed: per item its id, the original's and the paraphrase's SHA-256, the attempt that passed and
each check's result, no text. PARAPHRASE.md is committed: counts.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import difflib
import hashlib
import json
import os
import random
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

from dotenv import load_dotenv  # noqa: E402

load_dotenv(REPO / '.env')

import milestones  # noqa: E402

from full_run_28092026 import subsamples  # noqa: E402
from full_run_28092026.run_traces import config as run_config, items as pool_items  # noqa: E402

OUT = HERE / 'paraphrase'
ATTEMPTS_FILE = OUT / 'attempts.jsonl'
POOL = OUT / 'pool.jsonl'
MANIFEST = OUT / 'manifest.jsonl'
WRITER = {'model': 'mistralai/mistral-large-2512', 'temperature': 0.7, 'max_tokens': 2048,
          'provider': {'sort': 'price', 'allow_fallbacks': True}}
ATTEMPTS = 3
RETRY_SLEEPS = (3, 5, 10, 15, 20, 30, 45, 60)   # pauses between the tries of one call; a refused call bills nothing
COPY = 0.75
LENGTH = (0.7, 1.5)
PROMPT = """Rewrite the engineering problem below in different words.

Rules:
1. Keep every number exactly as written, with its unit. A reaction order, an exponent, a percentage or a tolerance is a number too: "order-2.1" stays "order-2.1", never "second-order".
2. Keep every symbol, variable name, subscript, formula, equation, chemical formula, acronym and technical term exactly as written, character for character: the same ^, *, /, _ and = with the same spacing, the same e-03 style exponents, the same (g) or (l) state labels, the same % signs, plain hyphen-minus signs. Acronyms such as CSTR, PFR or NPSH stay as acronyms; an equation such as "F_A0 = 3.96 mol/s" keeps its = sign and is not written as "equals" or "is".
3. Do not change any notation or formatting. Add no LaTeX, Markdown, Unicode superscripts, subscripts, minus signs or multiplication dots that the original does not have; where the original wraps something in $...$ or other math delimiters, keep those delimiters exactly where they are. Copy any table or block of data unchanged.
4. Keep the parts of the problem, and any lettered or numbered sub-questions, in the same order.
5. Do not add, remove or change any information, assumption or instruction. Do not hint at the method or the answer.
6. Change the wording and the sentence structure throughout the prose: different verbs and connectives, clauses reordered within a sentence, sentences split or joined, so that the result reads as a genuinely different phrasing and not as a copy, while staying about the same length as the original.

Return only the rewritten problem, with no preamble, heading or comment.

Problem:
{question}"""
PROMPT_SHA = hashlib.sha256(PROMPT.encode('utf-8')).hexdigest()


# ------------------------------------------------------------------ the checks

STRIP = ',.;:!?"\'()[]{}'
# Units written as a bare word. Not `in`: as the preposition it is everywhere, and the number check
# holds an inch's value anyway.
UNITS = {'m', 's', 'kg', 'g', 'K', 'N', 'J', 'W', 'V', 'Pa', 'Hz', 'mol', 'rad', 'bar', 'atm', 'psi', 'ft',
         'lb', 'lbf', 'lbm', 'mm', 'cm', 'km', 'min', 'h', 'hr', 'L', 'mL', 'ms', 'ns', 'dB', 'C', 'F',
         'kPa', 'MPa', 'GPa', 'kN', 'kJ', 'kW', 'MW', 'mA', 'kV', 'mV', 'kHz', 'MHz', 'GHz', 'Btu', 'hp', 'rpm'}
PARTS = re.compile(r'\((?:[a-h]|i{1,3}|iv|vi{0,3}|ix|x)\)|(?<![\w(])[a-h]\)|\b[Pp]art\s+[A-Za-z0-9]+\b')
PREAMBLE = re.compile(r"(?i)^\s*(?:here(?:'s| is)|sure|certainly|rewritten|paraphrase|the rewritten|revised)\b")


def tokens(text: str) -> list[str]:
    return [t.strip(STRIP) for t in text.split() if t.strip(STRIP)]


def technical(tok: str) -> bool:
    if tok in ('a', 'A', 'I'):
        return False
    return (tok in UNITS
            or (any(c.isdigit() for c in tok) and any(c.isalpha() for c in tok))
            or any(c in tok for c in '_^/=\\·×°%')
            or any(ord(c) > 127 for c in tok)
            or (len(tok) >= 2 and tok.isupper())
            or bool(re.search(r'[a-z][A-Z]', tok)))


def number_key(v: float) -> str:
    return f'{v:.12g}'


def check(original: str, text: str) -> dict:
    """Each check's result for one attempt; `passed` when all hold."""
    words_o = re.findall(r'\w+', original.lower())
    words_p = re.findall(r'\w+', text.lower())
    have = collections.Counter(tokens(text))
    need = collections.Counter(t for t in tokens(original) if technical(t))
    missing = sorted((need - have).elements())
    nums_o = collections.Counter(map(number_key, milestones.numbers(original)))
    nums_p = collections.Counter(map(number_key, milestones.numbers(text)))
    sim = difflib.SequenceMatcher(None, words_o, words_p, autojunk=False).ratio()
    ratio = len(text) / max(1, len(original))
    res = {'numbers': nums_o == nums_p,
           'tokens': not missing,
           'parts': PARTS.findall(original) == PARTS.findall(text),
           'copy': sim <= COPY,
           'length': LENGTH[0] <= ratio <= LENGTH[1],
           'clean': bool(text.strip()) and not PREAMBLE.search(text) and '```' not in text,
           'similarity': round(sim, 3), 'length_ratio': round(ratio, 3),
           'numbers_changed': sorted(((nums_o - nums_p) + (nums_p - nums_o)).elements())[:6],
           'tokens_missing': missing[:6]}
    res['passed'] = all(res[k] for k in ('numbers', 'tokens', 'parts', 'copy', 'length', 'clean'))
    return res


# ------------------------------------------------------------------ the writer

def selection() -> list[dict]:
    by_id = {it['item_id']: it for it in pool_items()}
    return [by_id[i] for i in subsamples.paraphrase_ids()]


def attempts_so_far() -> dict[str, list[dict]]:
    out = collections.defaultdict(list)
    if ATTEMPTS_FILE.exists():
        for line in ATTEMPTS_FILE.read_text(encoding='utf-8').splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            out[r['item_id']].append(r)
    return out


def write_one(cli, question: str) -> dict:
    """One call to the writer. Mistral throttles OpenRouter's shared capacity upstream and a refused call
    bills nothing, so a failed call is tried again after each pause in RETRY_SLEEPS, with jitter; only
    when they are spent is it written down as a service failure. `tries` and `seconds` record the wait."""
    t0 = time.time()
    last = ''
    for tries, pause in enumerate(RETRY_SLEEPS + (None,), start=1):
        try:
            r = cli.chat.completions.create(
                model=WRITER['model'], messages=[{'role': 'user', 'content': PROMPT.format(question=question)}],
                temperature=WRITER['temperature'], max_tokens=WRITER['max_tokens'],
                extra_body={'provider': WRITER['provider']})
            u = r.usage
            return {'text': (r.choices[0].message.content or '').strip(), 'served_model': getattr(r, 'model', None),
                    'provider': (getattr(r, 'model_extra', None) or {}).get('provider'),
                    'finish_reason': r.choices[0].finish_reason, 'prompt_tokens': getattr(u, 'prompt_tokens', None),
                    'completion_tokens': getattr(u, 'completion_tokens', None),
                    'billed_usd': (getattr(u, 'model_extra', None) or {}).get('cost'),
                    'tries': tries, 'seconds': round(time.time() - t0, 2)}
        except Exception as exc:                                   # noqa: BLE001
            last = f'{type(exc).__name__}: {exc}'[:300]
            if pause is None:
                break
            time.sleep(pause + random.uniform(0, pause / 2))
    return {'text': '', 'error': last, 'tries': tries, 'billed_usd': 0.0}


def item_until_pass(cli, it: dict, prior: list[dict]) -> list[dict]:
    """Attempts for one item until one passes or ATTEMPTS are spent, counting earlier runs' attempts."""
    rows = []
    n = len(prior)
    if any(r['check']['passed'] for r in prior):
        return rows
    while n < ATTEMPTS:
        n += 1
        res = write_one(cli, it['question'])
        if res.get('error'):
            rows.append({'item_id': it['item_id'], 'attempt': n, 'status': 'service_failure', **res})
            n -= 1                           # a service failure is not an attempt; the next run retries it
            break
        chk = check(it['question'], res['text'])
        rows.append({'item_id': it['item_id'], 'template_id': it['template_id'], 'attempt': n,
                     'original_sha256': it['sha256'], 'prompt_sha256': PROMPT_SHA, 'writer': WRITER,
                     'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **res, 'check': chk})
        if chk['passed']:
            break
    return rows


def run(workers: int, limit: int = 0) -> int:
    from openai import OpenAI
    cfg = run_config()
    key = os.getenv(cfg['route']['key_env'])
    if not key:
        raise SystemExit(f"{cfg['route']['key_env']} is not set in .env")
    cli = OpenAI(api_key=key, base_url=cfg['route']['base_url'], timeout=300.0, max_retries=0)
    OUT.mkdir(exist_ok=True)
    prior = attempts_so_far()
    todo = [it for it in selection() if not any(r.get('check', {}).get('passed') for r in prior[it['item_id']])
            and len([r for r in prior[it['item_id']] if 'check' in r]) < ATTEMPTS]
    if limit:
        todo = todo[:limit]                  # a pilot of the prompt; the next run without --limit takes the rest
    print(f'{len(todo)} items to write')
    billed, done = 0.0, 0
    with open(ATTEMPTS_FILE, 'a', encoding='utf-8', newline='\n') as fh, cf.ThreadPoolExecutor(workers) as pool:
        futs = [pool.submit(item_until_pass, cli, it, [r for r in prior[it['item_id']] if 'check' in r])
                for it in todo]
        for fut in cf.as_completed(futs):
            for row in fut.result():
                fh.write(json.dumps(row, ensure_ascii=False) + '\n')
                fh.flush()
                billed += row.get('billed_usd') or 0.0
            done += 1
            if done % 25 == 0 or done == len(todo):
                print(f'  {done}/{len(todo)} items, billed ${billed:.3f}')
    print(f'billed this invocation: ${billed:.3f}')
    return rebuild()


# ------------------------------------------------------------------ outputs

def rebuild() -> int:
    """pool.jsonl and manifest.jsonl from the attempts, every check re-run; PARAPHRASE.md."""
    its = selection()
    prior = attempts_so_far()
    manifest, pool, fails, n_att = [], [], collections.Counter(), collections.Counter()
    for it in its:
        rows = sorted((r for r in prior[it['item_id']] if 'text' in r and 'check' in r), key=lambda r: r['attempt'])
        chosen = None
        for r in rows:
            r['check'] = check(it['question'], r['text'])
            if r['check']['passed']:
                chosen = r
                break
            fails.update(k for k in ('numbers', 'tokens', 'parts', 'copy', 'length', 'clean') if not r['check'][k])
        n_att[len(rows) if chosen is None else chosen['attempt']] += 1
        last = chosen or (rows[-1] if rows else None)
        manifest.append({'item_id': it['item_id'], 'template_id': it['template_id'], 'branch': it['branch'],
                         'original_sha256': it['sha256'], 'passed': chosen is not None,
                         'sha256': hashlib.sha256(chosen['text'].encode('utf-8')).hexdigest() if chosen else None,
                         'attempts': len(rows), 'attempt_passed': chosen['attempt'] if chosen else None,
                         'check': {k: v for k, v in last['check'].items()
                                   if k not in ('numbers_changed', 'tokens_missing')} if last else None,
                         'served_model': last.get('served_model') if last else None})
        if chosen:
            pool.append({'item_id': it['item_id'], 'template_id': it['template_id'], 'question': chosen['text']})
    OUT.mkdir(exist_ok=True)
    with open(MANIFEST, 'w', encoding='utf-8', newline='\n') as fh:
        for m in manifest:
            fh.write(json.dumps(m) + '\n')
    with open(POOL, 'w', encoding='utf-8', newline='\n') as fh:
        for r in pool:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
    written = [m for m in manifest if m['attempts']]
    passed = [m for m in manifest if m['passed']]
    billed = sum(r.get('billed_usd') or 0.0 for rs in prior.values() for r in rs)
    sims = sorted(m['check']['similarity'] for m in passed)
    by_branch = collections.Counter(m['branch'] for m in passed)
    L = ['# The paraphrases (D-141)', '',
         'Generated by `paraphrase.py`; the writer, the prompt and the checks are defined in its docstring. '
         'Counts and hashes only: the text stays local.', '',
         '| | |', '|---|---|',
         f'| items selected (subsamples.PARAPHRASE) | {len(its)} |',
         f'| items written | {len(written)} |',
         f'| passing a paraphrase | {len(passed)} |',
         f'| passing at attempt 1, 2, 3 | {n_att[1]}, {n_att[2]}, {n_att[3]} |',
         f'| no paraphrase after {ATTEMPTS} attempts | {sum(1 for m in written if not m["passed"])} |',
         f'| passing, per branch | ' + ', '.join(f'{b.replace("_engineering", "")} {n}' for b, n in sorted(by_branch.items())) + ' |',
         f'| word similarity of the passing ones: min, median, max | '
         + (f'{sims[0]}, {sims[len(sims) // 2]}, {sims[-1]}' if sims else '-') + ' |',
         f'| failed attempts by check | ' + (', '.join(f'{k} {v}' for k, v in fails.most_common()) or 'none') + ' |',
         f'| billed | ${billed:.3f} |',
         f'| writer | `{WRITER["model"]}`, temperature {WRITER["temperature"]}; prompt sha256 `{PROMPT_SHA[:16]}` |']
    (HERE / 'PARAPHRASE.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


def dry_run() -> int:
    import urllib.request
    its = selection()
    chars = [len(PROMPT.format(question=it['question'])) for it in its]
    out_chars = [len(it['question']) for it in its]
    tin, tout = sum(chars) / 3.5, sum(out_chars) / 3.5          # about 3.5 characters a token
    req = urllib.request.Request(f"https://openrouter.ai/api/v1/models/{WRITER['model']}/endpoints",
                                 headers={'User-Agent': 'engtrace-full-run'})
    eps = json.load(urllib.request.urlopen(req, timeout=60))['data']['endpoints']
    best = min(eps, key=lambda e: float(e['pricing']['prompt']) * tin + float(e['pricing']['completion']) * tout)
    pin, pout = float(best['pricing']['prompt']), float(best['pricing']['completion'])
    one = pin * tin + pout * tout
    per_branch = collections.Counter(it['branch'] for it in its)
    print(f'selection: {len(its)} items, {len({i["template_id"] for i in its})} templates; '
          + ', '.join(f'{b} {n}' for b, n in sorted(per_branch.items())))
    print(f'writer: {WRITER["model"]} via {best.get("provider_name")} (quantization '
          f'{best.get("quantization") or "not declared"}), ${pin * 1e6:.2f} in / ${pout * 1e6:.2f} out per million')
    print(f'prompt sha256 {PROMPT_SHA[:16]}; about {tin / len(its):.0f} input and {tout / len(its):.0f} output '
          f'tokens per attempt')
    print(f'estimate: ${one:.2f} if every item passes at the first attempt; ${one * ATTEMPTS:.2f} if every item '
          f'needs all {ATTEMPTS}. Nothing was called.')
    return 0


def status() -> int:
    prior = attempts_so_far()
    rows = [r for rs in prior.values() for r in rs]
    print(f'{len(prior)} items attempted, {len(rows)} attempts, '
          f'{sum(any(r.get("check", {}).get("passed") for r in rs) for rs in prior.values())} passing, '
          f'billed ${sum(r.get("billed_usd") or 0.0 for r in rows):.3f}')
    return 0


def selftest() -> int:
    """The checks on constructed cases with a known verdict; and every original against itself, which
    must fail the copy check and nothing else. Calls nothing."""
    q = ('A first-order liquid-phase reaction of Tetrahydrofuran (THF) occurs in a batch reactor. The initial '
         'concentration is 1.23 mol/L, and after reaction, the concentration decreases to 0.36 mol/L. If the '
         'first-order rate constant is 0.00171 s⁻¹, determine the reaction time required.')
    good = ('In a batch reactor, Tetrahydrofuran (THF) undergoes a liquid-phase reaction of first order. Its '
            'concentration starts at 1.23 mol/L and falls to 0.36 mol/L. Given a first-order rate constant of '
            '0.00171 s⁻¹, find how long the reaction must run.')
    cases = [(good, True, None),
             (good.replace('0.36', '0.63'), False, 'numbers'),
             (good.replace('0.00171 s⁻¹', '0.00171 per second'), False, 'tokens'),
             ('Here is the rewritten problem: ' + good, False, 'clean'),
             (q.replace('determine', 'find'), False, 'copy'),
             (good[:120], False, 'length')]
    two = '(a) Find the flow rate Q in m^3/s. (b) Find the head loss h_f in m.'
    cases_parts = [('(a) Determine Q in m^3/s; (b) determine h_f in m, the head loss.', True),
                   ('(b) Determine h_f in m, the head loss; (a) determine Q in m^3/s.', False)]
    bad = []
    for text, want, reason in cases:
        c = check(q, text)
        if c['passed'] != want or (reason and c[reason]):
            bad.append(f'{reason or "good"}: {c}')
    for text, want in cases_parts:
        if check(two, text)['parts'] != want:
            bad.append(f'parts: {text}')
    fails = collections.Counter()
    for it in selection():
        c = check(it['question'], it['question'])
        fails.update(k for k in ('numbers', 'tokens', 'parts', 'length', 'clean') if not c[k])
        if c['copy']:
            fails['not a copy'] += 1
    if fails:
        bad.append(f'identity: {dict(fails)}')
    print('selftest:', 'all pass' if not bad else bad)
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--workers', type=int, default=8, help='parallel items; the upstream throttle is per call, so more workers give more throughput')
    ap.add_argument('--yes', action='store_true', help='required for the writing run, which bills')
    ap.add_argument('--limit', type=int, default=0, help='write at most this many unresolved items (a pilot of the prompt)')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.dry_run:
        return dry_run()
    if a.check:
        return rebuild()
    if a.status:
        return status()
    if not a.yes:
        raise SystemExit('writing the paraphrases bills: re-run with --yes once the spend is approved '
                         '(see --dry-run for the estimate)')
    return run(a.workers, a.limit)


if __name__ == '__main__':
    raise SystemExit(main())
