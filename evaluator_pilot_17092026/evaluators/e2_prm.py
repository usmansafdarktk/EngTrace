"""E2 - an open process reward model (PRM) scores every step. No API, no judge prompt.

The PRMs run on HiPerGator GPUs (hpg/: e2_score.py, e2.sbatch; D-090). Their raw
outputs - one probability per step per trace, from each PRM's card-specified usage
- are downloaded to scores/e2_raw/<prm>_full.jsonl with their .meta.json. This
module turns them into harness rows, so E2 passes the same frozen-slice and
verify_traces gates as every other evaluator and carries a config hash.

    PRMs   qwen72  Qwen2.5-Math-PRM-72B   PRIMARY
           versa   VersaPRM (LoRA on Llama-PRM800K), the multi-domain PRM
           qwen7   Qwen2.5-Math-PRM-7B    size cross-check

A STEP is correct when its reward is >= 0.5, the threshold these PRMs are used with
(and the level the self-check holds them to; hpg/README.md). Per PRM and trace:

    frac_ok     share of steps judged correct           (primary score: qwen72's)
    min         the lowest step reward                  (the usual PRM aggregation)
    clean       1 if no step is flagged
    first_flag  index of the first flagged step, or None

A trace over the PRM's context limit is not scored (over_length): a truncated trace
would lose its last steps and look clean. Qwen's limit is 4,096 tokens; 2 of 360
traces exceed it.

EVERY RAW ROW IS CHECKED against the trace being scored: its trace_sha256 must equal
the trace's, and its step count must equal what E0's extract_steps gives on that text
(the split build_e2_inputs.py used). A mismatch raises, so a stale or misaligned raw
file cannot be scored silently.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_ROOT = os.path.dirname(_PILOT)
if os.path.join(_ROOT, 'evaluation') not in sys.path:
    sys.path.insert(0, os.path.join(_ROOT, 'evaluation'))

ID = 'e2'
DESCRIPTION = 'Open process reward models score every step (Qwen2.5-Math-PRM-72B primary; VersaPRM; 7B)'
RAW = os.path.join(_PILOT, 'scores', 'e2_raw')
PRMS = ('qwen72', 'versa', 'qwen7')
PRIMARY = 'qwen72'
VERDICT_AT = 0.5
DEVIATIONS = [
    'Qwen-72B misses its card example rewards by up to 0.21 on two borderline steps '
    '(attention kernel numerics on Blackwell; every verdict agrees). Self-check is on '
    'verdicts at 0.5 (D-090).',
    'Qwen runs with eager attention (closest to the cards); VersaPRM with the library default.',
]


def _sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def _raw(prm):
    return os.path.join(RAW, '%s_full.jsonl' % prm)


def config() -> dict:
    metas = {}
    for p in PRMS:
        m = json.load(open(_raw(p) + '.meta.json', encoding='utf-8'))
        metas[p] = {k: m.get(k) for k in ('repo', 'revision', 'base', 'attn_implementation', 'scorer_md5')}
        metas[p]['selfcheck_passed'] = m['selfcheck']['passed']
    return {
        'evaluator': ID,
        'e2_sha256': _sha(__file__),
        'raw_sha256': {p: _sha(_raw(p)) for p in PRMS},
        'prms': metas,
        'primary': PRIMARY,
        'verdict_at': VERDICT_AT,
    }


def steps_of(text: str) -> list[str]:
    """The split build_e2_inputs.py used: E0's extract_steps, stripped, empties dropped."""
    from engineering_parser import extract_steps
    return [s for s in (x.strip() for x in extract_steps(text)[0]) if s]


def setup(dry_run: bool = False):
    state = {}
    for p in PRMS:
        meta = json.load(open(_raw(p) + '.meta.json', encoding='utf-8'))
        if not meta['selfcheck']['passed']:
            raise SystemExit('%s: raw file was produced by a run whose self-check failed' % p)
        state[p] = {(r['model_key'], r['item_id']): r
                    for r in (json.loads(l) for l in open(_raw(p), encoding='utf-8'))}
    return state


def summarise(probs):
    if probs is None:
        return {'frac_ok': None, 'min': None, 'clean': None, 'first_flag': None}
    ok = [x >= VERDICT_AT for x in probs]
    return {'frac_ok': sum(ok) / len(ok) if ok else None,
            'min': min(probs) if probs else None,
            'clean': int(all(ok)) if ok else None,
            'first_flag': ok.index(False) if not all(ok) else None}


def score(state, item: dict, trace: dict, seed: int) -> dict:
    key = (trace.get('model_key') or trace.get('model'), item['item_id'])
    tsha = hashlib.sha256(trace['text'].encode('utf-8')).hexdigest()
    n_steps = len(steps_of(trace['text']))
    scores, meta = {}, {'prm': {}, 'tribunal_reached_judges': False}
    for p in PRMS:
        r = state[p].get(key)
        if r is None:
            raise KeyError('%s has no raw row for %s/%s' % (p, key[0], key[1]))
        if r['trace_sha256'] != tsha:
            raise ValueError('%s raw row for %s/%s is for a different trace text' % (p, key[0], key[1]))
        if r['n_steps'] != n_steps:
            raise ValueError('%s raw row for %s/%s has %d steps, the trace splits into %d'
                             % (p, key[0], key[1], r['n_steps'], n_steps))
        s = summarise(None if r['over_length'] else r['step_probs'])
        for k, v in s.items():
            scores['%s_%s' % (p, k)] = v
        meta['prm'][p] = {'step_probs': r['step_probs'], 'over_length': r['over_length'],
                          'n_tokens': r['n_tokens'], 'max_len': r['max_len']}
    scores['e2'] = scores['%s_frac_ok' % PRIMARY]
    scores['e2_min'] = scores['%s_min' % PRIMARY]
    scores['n_steps'] = n_steps
    return {'scores': scores, 'meta': meta, 'calls': [], 'triggered': False}
