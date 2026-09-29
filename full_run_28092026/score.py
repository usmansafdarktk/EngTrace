"""The full run's scorer: the deterministic evaluators over the traces, every raw output kept.

    python -m full_run_28092026.score --gold                      # FREE: the stack on the 2,250 gold solutions
    python -m full_run_28092026.score [--variant main] [--model KEY ...] [--workers 8]   # FREE
    python -m full_run_28092026.score --status

It runs the pilot's evaluators, imported unmodified from evaluator_pilot_17092026/evaluators/:

  answer check    answer.verdict at the fitted tolerance and at half and double it (the plan's
                  sensitivity), with the item's milestones as the computed targets
  E3              e3_milestones.reach over the item's milestones, derived once by milestones.build
  E4 digit rule   arith.check on each step, the steps split by e2_prm.steps_of: the split the
                  experts' step labels were made on, so every step-level score here, a judge's and
                  a later router's included, refers to the same steps

The paid judge stage (E5) is judge.py; it reads this store and writes beside it.

WHAT A ROW HOLDS (scores/<variant>/<model>.jsonl, gitignored with the traces). One row per final
trace row, as run_traces.existing reads them:
  the item      item_id, template_id, instance_index, branch, domain, level, answer_type,
                milestone count, single_path (one reasoning path across its pool instances,
                diversity.py's lower reading), repeat_of
  the trace     state, finish_reason, provider, served_model, tokens, billed_usd, seconds,
                attempts, mode, and the SHA-256 of its text; the text itself stays in the trace
  the answer    label at three tolerances, the targets and how many matched, readable, unusable,
                and the score: correct 1, partial 0.5, incorrect 0, unusable 0 (D-117)
  E3            every milestone's reached flag and unit scaling; coverage, None when an item has
                no milestones, so it is left out of milestone aggregates (D-120)
  E4            per step: claims checked, flagged by the digit rule, flagged at 1%, judged
So every analysis the plan names - instance variation, consistency, the cliff, the tolerance
sensitivity, error attribution - reads this store, and the free parts can be recomputed at will.

VARIANTS. main reads traces/<model>.jsonl; any other variant reads traces/<variant>/<model>.jsonl
(the paraphrase and decoding-repeat runs). A paraphrased trace is scored against its ORIGINAL item:
gold, milestones and the answer check's question are the original's, so the two arms differ only in
what the model wrote (ANALYSIS_PLAN Q5).
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402
import arith  # noqa: E402
import e2_prm  # noqa: E402
import e3_milestones as e3  # noqa: E402

SCORES = HERE / 'scores'
TRACES = HERE / 'traces'
POOL = HERE / 'pool'
MILESTONES = SCORES / 'milestones.json'
DONE = ('answered', 'empty')
SCORE_OF = {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}
TOLS = {'half': answer.REL / 2, 'fitted': answer.REL, 'double': answer.REL * 2}
EVALUATOR_FILES = [EVAL / 'answer.py', EVAL / 'arith.py', EVAL / 'milestones.py', EVAL / 'e3_milestones.py',
                   EVAL / 'e2_prm.py', REPO / 'evaluation' / 'engineering_parser.py', Path(__file__)]
WORDS = set(answer.VERDICT_WORDS) | set(answer.LINEAR)


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def git(*args) -> str:
    return subprocess.run(['git', *args], cwd=REPO, capture_output=True, text=True).stdout.strip()


# ------------------------------------------------------------------ the items

def pool_items() -> dict[str, dict]:
    """The frozen pool's items, text from pool/ and metadata from the committed manifest."""
    man = {json.loads(l)['item_id']: json.loads(l)
           for l in (HERE / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()}
    single = {r['template_id']: r['pool']['paths_lower'] == 1
              for r in json.loads((HERE / 'diversity.json').read_text(encoding='utf-8'))}
    items = {}
    for path in sorted(POOL.rglob('*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            m = man[r['item_id']]
            items[r['item_id']] = {
                'item_id': r['item_id'], 'template_id': m['template_id'], 'seed': r['seed'],
                'instance_index': m['instance_index'], 'branch': m['branch'], 'domain': m['domain'],
                'area': m['area'], 'level': m['level'], 'answer_type': m['answer_type'],
                'question': r['question'], 'solution': r['solution'], 'sha256': m['sha256'],
                'repeat_of': m.get('repeat_of'), 'single_path': single[m['template_id']],
            }
    return items


def milestone_sets(items: dict) -> dict[str, list[dict]]:
    """Each item's milestones, derived once by milestones.build and cached (values only)."""
    if MILESTONES.exists():
        cache = json.loads(MILESTONES.read_text(encoding='utf-8'))
        if set(cache) == set(items):
            return cache
    import milestones
    import tests.template_integrity.core as core
    refs = core.discover()
    core.discover = lambda *a, **k: refs
    cache = {i: milestones.build(it)['milestones'] for i, it in items.items()}
    SCORES.mkdir(exist_ok=True)
    MILESTONES.write_text(json.dumps(cache) + '\n', encoding='utf-8')
    return cache


def trace_rows(variant: str, key: str) -> list[dict]:
    """The final row per item, the last answered or empty one, as run_traces.existing reads it."""
    path = (TRACES / f'{key}.jsonl') if variant == 'main' else (TRACES / variant / f'{key}.jsonl')
    final = {}
    for ln in path.read_text(encoding='utf-8').splitlines():
        try:
            r = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if r.get('status') in DONE:
            final[r['item_id']] = r
    return list(final.values())


# ------------------------------------------------------------------ scoring one trace

def readable(text: str) -> bool:
    """Whether the answer check can read a final answer at all (D-117's 'never states one')."""
    if answer.HEADING.search(text) or answer.ANSWER.search(text):
        return True
    tail = text[-answer.WINDOW:]
    return bool(answer.values(tail)) or any(w in WORDS for w in answer.WORD.findall(tail.lower()))


def score_answer(item: dict, text: str, ms: list[dict]) -> dict:
    vals = tuple(m['value'] for m in ms)
    labels, info = {}, None
    for name, tol in TOLS.items():
        lab, det = answer.verdict(text, item, tol, vals)
        labels[name] = lab
        if name == 'fitted':
            info = det
    t = info.get('targets', {})
    return {'label': labels['fitted'], 'label_half_tol': labels['half'], 'label_double_tol': labels['double'],
            'matched': info.get('matched'), 'of': info.get('of'),
            'targets': {'numbers': [v for v, _u in t.get('numbers', [])], 'words': t.get('words', []),
                        'labeled': t.get('labeled', []), 'exact': t.get('exact', False)}}


def score_steps(text: str) -> list[dict]:
    out = []
    for step in e2_prm.steps_of(text):
        claims = arith.check(step).claims
        out.append({'claims': len(claims), 'digit_flags': sum(not c.ok_digit for c in claims),
                    'tol1_flags': sum(not c.ok for c in claims),
                    'judged': sum(c.ulp is not None for c in claims)})
    return out


def score_trace(item: dict, row: dict, ms: list[dict]) -> dict:
    text = row.get('text') or ''
    base = {k: item[k] for k in ('item_id', 'template_id', 'instance_index', 'branch', 'domain', 'level',
                                 'answer_type', 'repeat_of', 'single_path')}
    base.update({k: row.get(k) for k in ('model_key', 'status', 'finish_reason', 'provider', 'served_model',
                                          'prompt_tokens', 'completion_tokens', 'reasoning_tokens',
                                          'billed_usd', 'seconds', 'attempts', 'mode')})
    base['trace_sha256'] = sha(text.encode('utf-8'))
    base['milestones_required'] = len(ms)
    if row.get('status') != 'answered':
        return {**base, 'readable': False, 'unusable': True, 'score': 0.0, 'answer': None,
                'e3': {'coverage': None if not ms else 0.0, 'reached': [False] * len(ms)}, 'steps': []}
    ok_to_read = readable(text)
    ans = score_answer(item, text, ms)
    hits = e3.reach(ms, text)
    steps = score_steps(text)
    unusable = not ok_to_read
    return {**base, 'readable': ok_to_read, 'unusable': unusable,
            'score': 0.0 if unusable else SCORE_OF[ans['label']], 'answer': ans,
            'e3': {'coverage': (sum(h['reached'] for h in hits) / len(hits)) if hits else None,
                   'reached': [h['reached'] for h in hits], 'scale': [h['scale'] for h in hits]},
            'steps': steps}


def _score_chunk(args):
    items, rows, ms = args
    return [score_trace(items[r['item_id']], r, ms[r['item_id']]) for r in rows]


# ------------------------------------------------------------------ modes

def config() -> dict:
    return {'evaluator_sha256': {str(p.relative_to(REPO)): sha(p.read_bytes()) for p in EVALUATOR_FILES},
            'git': git('rev-parse', 'HEAD'), 'tols': TOLS,
            'scored_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()}


def run(variant: str, keys: list[str], workers: int) -> int:
    items = pool_items()
    ms = milestone_sets(items)
    out_dir = SCORES / variant
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = config()
    (out_dir / 'CONFIG.json').write_text(json.dumps(cfg, indent=1) + '\n', encoding='utf-8')
    for key in keys:
        rows = trace_rows(variant, key)
        unknown = [r['item_id'] for r in rows if r['item_id'] not in items]
        if unknown:
            raise SystemExit(f'{key}: {len(unknown)} rows are not pool items')
        chunks = [rows[i:i + 50] for i in range(0, len(rows), 50)]
        scored = []
        with cf.ProcessPoolExecutor(max_workers=workers) as pool:
            for part in pool.map(_score_chunk, [(items, c, ms) for c in chunks]):
                scored.extend(part)
        scored.sort(key=lambda r: (r['template_id'], r['instance_index']))
        with open(out_dir / f'{key}.jsonl', 'w', encoding='utf-8', newline='\n') as fh:
            for r in scored:
                fh.write(json.dumps(r) + '\n')
        n = len(scored)
        print(f'{variant}/{key}: {n} rows, mean score {sum(r["score"] for r in scored) / n:.3f}, '
              f'unusable {sum(r["unusable"] for r in scored)}', flush=True)
    return 0


def gold(workers: int) -> dict:
    """The stack on every gold solution: all correct, E3 complete, no digit-rule flag."""
    items = pool_items()
    ms = milestone_sets(items)
    rows = [{'item_id': i, 'status': 'answered', 'text': it['solution']} for i, it in items.items()]
    chunks = [rows[i:i + 50] for i in range(0, len(rows), 50)]
    scored = []
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for part in pool.map(_score_chunk, [(items, c, ms) for c in chunks]):
            scored.extend(part)
    res = {
        'items': len(scored),
        'answer_correct': sum(r['answer']['label'] == 'correct' for r in scored),
        'answer_correct_all_three_tols': sum(r['answer']['label'] == r['answer']['label_half_tol']
                                             == r['answer']['label_double_tol'] == 'correct' for r in scored),
        'unusable': sum(r['unusable'] for r in scored),
        'e3_complete': sum(r['e3']['coverage'] == 1.0 for r in scored),
        'e3_no_milestones': sum(r['e3']['coverage'] is None for r in scored),
        'digit_flags': sum(s['digit_flags'] for r in scored for s in r['steps']),
        'claims_checked': sum(s['claims'] for r in scored for s in r['steps']),
        'steps': sum(len(r['steps']) for r in scored),
    }
    print(json.dumps(res, indent=1))
    return res


def status() -> int:
    for d in sorted(p for p in SCORES.iterdir() if p.is_dir()) if SCORES.exists() else []:
        for f in sorted(d.glob('*.jsonl')):
            rows = [json.loads(l) for l in f.read_text(encoding='utf-8').splitlines()]
            print(f'{d.name}/{f.stem}: {len(rows)} rows')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--gold', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--variant', default='main')
    ap.add_argument('--model', action='append')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    if a.gold:
        gold(a.workers)
        return 0
    if a.status:
        return status()
    from full_run_28092026.run_traces import config as roster
    keys = a.model or [m['key'] for m in roster()['models']
                       if (TRACES / (f"{m['key']}.jsonl" if a.variant == 'main'
                                     else f"{a.variant}/{m['key']}.jsonl")).exists()]
    return run(a.variant, keys, a.workers)


if __name__ == '__main__':
    raise SystemExit(main())
