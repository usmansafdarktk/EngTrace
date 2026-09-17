"""E0 - the published EngTrace framework, run as close to unchanged as still exists.

The point of E0 is to be the anchor: every other candidate is compared against it,
so it must be the framework the paper describes and not a tidied version of it.
This adapter therefore IMPORTS evaluation/engtrace_evaluation_framework.py and
calls its own `evaluate_entry`. It does not copy, reimplement or edit a line of
Tier 1, the Tribunal prompt, the JSON parsing, the majority vote, the Hungarian
matching or the metrics. The framework file's SHA-256 is recorded on every row, so
an edit to it shows up as a config change and forces a re-score.

WHAT COULD NOT BE KEPT, each recorded on every row under meta.deviations:

  D1  Google judge.  `gemini-3-pro-preview` returns 404 - "no longer available,
      use gemini-3.1-pro-preview" (probed 2026-09-17).  Google's named successor is
      used.  This is the one change to WHO judges, and it is forced.

  D2  Transport for the OpenAI and Anthropic judges.  Both direct keys return 401,
      so both judges are reached through OpenRouter.  Same models: `openai/gpt-5`,
      and `anthropic/claude-opus-4.5` for `claude-opus-4-5-20251101` - OpenRouter's
      id does not name the dated snapshot, and Opus 4.5 has only the one.  The
      Anthropic judge still gets max_tokens=2048, exactly as the framework sets it.

  D3  Seeded sampling.  For a WRONG final answer the framework sends the trace to
      the Tribunal with probability 0.20 via an unseeded `random.random()`, so two
      runs judge different traces.  `random` is seeded per (item, model) just
      before each entry.  The rule and the rate are untouched; the draw is now
      repeatable, and the seed is on the row.

WHAT IS OBSERVED, WITHOUT BEING CHANGED.  `_call_single_judge` swallows every
failure and returns [], and `_tier2_tribunal_batch` silently drops Google if
`get_model_info` raises.  Either way a step simply receives fewer votes, or none,
and keeps a score of 0 - a failed judge is indistinguishable from a judge that
said "wrong".  That is this repository's recurring defect, so every judge call is
logged (raw text, finish reason, tokens, error) and each row reports how many
judges were called and how many returned parseable results.  The framework's
behaviour on those failures is kept; it is just no longer invisible.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import random
import sys
import time
from contextlib import redirect_stdout
from types import SimpleNamespace

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_ROOT = os.path.dirname(_PILOT)
_EVAL = os.path.join(_ROOT, 'evaluation')
FRAMEWORK = os.path.join(_EVAL, 'engtrace_evaluation_framework.py')

ID = 'e0'
DESCRIPTION = 'Published framework: Tier 1 numeric + cross-encoder, Tier 2 three-judge Tribunal'

JUDGES = {
    'openai':    {'model': 'openai/gpt-5',              'route': 'openrouter', 'published': 'gpt-5'},
    'anthropic': {'model': 'anthropic/claude-opus-4.5', 'route': 'openrouter', 'published': 'claude-opus-4-5-20251101'},
    'google':    {'model': 'gemini-3.1-pro-preview',    'route': 'google',     'published': 'gemini-3-pro-preview'},
}

DEVIATIONS = [
    'D1 google judge gemini-3-pro-preview -> gemini-3.1-pro-preview (published id returns 404)',
    'D2 openai and anthropic judges via OpenRouter (direct keys 401); anthropic snapshot unpinned',
    'D3 random seeded per (item, model) before each entry; the 0.20 wrong-answer sample rate is unchanged',
]


def _sha_file(path: str) -> str:
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def libraries() -> dict:
    """Scorer library versions. They change Tier 1 scores, so they are config.

    Not hypothetical: under transformers 5.x, bert_score raises
    `OverflowError: int too big to convert` on every entry and the framework
    catches it and returns BERTScore 0.0 - a column of zeros that looks like data.
    The venv pins the last pre-v5 line (4.57.x), which is also the line current
    when the published run was made (its Opus 4.5 judge dates it to ~Dec 2025).
    """
    from importlib.metadata import version, PackageNotFoundError
    out = {}
    for pkg in ('torch', 'transformers', 'sentence-transformers', 'bert-score',
                'rouge-score', 'tokenizers', 'scipy', 'numpy'):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def config() -> dict:
    """Everything that changes a score. The harness hashes this to decide re-scoring."""
    return {
        'evaluator': ID,
        'framework_sha256': _sha_file(FRAMEWORK),
        'parser_sha256': _sha_file(os.path.join(_EVAL, 'engineering_parser.py')),
        'judges': {k: v['model'] for k, v in JUDGES.items()},
        'libraries': libraries(),
        'deviations': DEVIATIONS,
    }


# ------------------------------------------------------------------ call log

class CallLog:
    """Every judge call for the entry being scored. Reset per entry."""

    def __init__(self):
        self.calls = []
        self.dry_run = False
        self.prompts = []

    def reset(self):
        self.calls, self.prompts = [], []


def _openrouter(log: CallLog, provider: str):
    from openai import OpenAI
    cli = OpenAI(api_key=os.environ['OPENROUTER_API_KEY'],
                 base_url='https://openrouter.ai/api/v1', timeout=300.0, max_retries=2)

    def create(**kw):
        rec = {'provider': provider, 'model': kw.get('model'), 'max_tokens': kw.get('max_tokens')}
        t0 = time.time()
        try:
            r = cli.chat.completions.create(**kw)
            ch = r.choices[0]
            rec.update(ok=True, served_model=r.model,
                       serving_provider=(r.model_extra or {}).get('provider'),
                       finish_reason=ch.finish_reason,
                       prompt_tokens=r.usage.prompt_tokens,
                       completion_tokens=r.usage.completion_tokens,
                       text=ch.message.content or '')
            return r
        except Exception as exc:                                   # noqa: BLE001
            rec.update(ok=False, error='%s: %s' % (type(exc).__name__, str(exc)[:300]))
            raise
        finally:
            rec['seconds'] = round(time.time() - t0, 2)
            log.calls.append(rec)

    return create


def build_clients(log: CallLog):
    """Objects shaped exactly like the SDK clients the framework calls."""
    or_openai = _openrouter(log, 'openai')
    openai_client = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(create=or_openai)))

    or_anthropic = _openrouter(log, 'anthropic')

    def anthropic_create(model, max_tokens, messages):
        # The framework reads resp.content[0].text; OpenRouter speaks chat.completions.
        r = or_anthropic(model=model, max_tokens=max_tokens, messages=messages)
        return SimpleNamespace(content=[SimpleNamespace(text=r.choices[0].message.content or '')])

    anthropic_client = SimpleNamespace(messages=SimpleNamespace(create=anthropic_create))
    return openai_client, anthropic_client


def wrap_google(log: CallLog, genai):
    """Log Gemini calls without changing how the framework makes them."""
    real = genai.GenerativeModel

    class LoggedModel:
        def __init__(self, name, *a, **kw):
            self._m = real(name, *a, **kw)
            self._name = name

        def generate_content(self, prompt, **kw):
            rec = {'provider': 'google', 'model': self._name}
            t0 = time.time()
            try:
                r = self._m.generate_content(prompt, **kw)
                um = getattr(r, 'usage_metadata', None)
                cand = r.candidates[0] if getattr(r, 'candidates', None) else None
                rec.update(ok=True, served_model=getattr(r, 'model_version', None),
                           finish_reason=str(getattr(cand, 'finish_reason', None)),
                           prompt_tokens=getattr(um, 'prompt_token_count', None),
                           completion_tokens=getattr(um, 'candidates_token_count', None),
                           # The native SDK DOES report thinking, unlike the
                           # OpenAI-compatible endpoint the traces came through.
                           thinking_tokens=getattr(um, 'thoughts_token_count', None))
                try:
                    rec['text'] = r.text
                except Exception as exc:                           # noqa: BLE001
                    rec['text'] = ''
                    rec['text_error'] = str(exc)[:200]
                return r
            except Exception as exc:                               # noqa: BLE001
                rec.update(ok=False, error='%s: %s' % (type(exc).__name__, str(exc)[:300]))
                raise
            finally:
                rec['seconds'] = round(time.time() - t0, 2)
                log.calls.append(rec)

    genai.GenerativeModel = LoggedModel


# --------------------------------------------------------------------- cache

CACHE = os.path.join(_PILOT, 'scores', '_cache', 'e0_tier1.jsonl')


def install_cache(fw, framework):
    """Memoise Tier 1 and BERTScore to disk, keyed by exact inputs and library versions.

    Both are pure functions of their arguments: the same step lists give the same
    matrix and the same texts give the same BERTScore (checked: 4 vs 8 torch
    threads, max score difference 0.0). On this CPU the cross-encoder costs 0.5-2s
    per step pair, so Tier 1 over 300 traces is hours - and a dry run followed by
    the real run would pay it twice for identical numbers. The cache changes when
    work happens, never what comes out; the library versions are in the key, so a
    cached matrix from one scorer stack is never served to another.
    """
    import numpy as np
    import torch
    torch.set_num_threads(os.cpu_count() or 4)

    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    store = {}
    if os.path.exists(CACHE):
        for ln in open(CACHE, encoding='utf-8'):
            try:
                rec = json.loads(ln)
                store[rec['k']] = rec['v']
            except (json.JSONDecodeError, KeyError):
                continue
    libs = json.dumps(libraries(), sort_keys=True)

    def key(kind, *args):
        return hashlib.sha256(json.dumps([kind, libs, args], sort_keys=True).encode()).hexdigest()

    def put(k, v):
        store[k] = v
        with open(CACHE, 'a', encoding='utf-8', newline='\n') as fh:
            fh.write(json.dumps({'k': k, 'v': v}) + '\n')

    real_matrix = framework._tier1_verify_matrix

    def cached_matrix(gt_s_txt, gt_s_val, pred_s_txt, pred_s_val):
        k = key('tier1', gt_s_txt, gt_s_val, pred_s_txt, pred_s_val)
        if k in store:
            v = store[k]
            return np.array(v['V'], dtype=float).reshape(v['shape'])
        V = real_matrix(gt_s_txt, gt_s_val, pred_s_txt, pred_s_val)
        put(k, {'V': V.ravel().tolist(), 'shape': list(V.shape)})
        return V

    framework._tier1_verify_matrix = cached_matrix

    real_bert = fw.safe_bert_score

    def cached_bert(gt, pred):
        k = key('bertscore', gt, pred)
        if k in store:
            return store[k]
        v = real_bert(gt, pred)
        # A failed BERTScore returns 0.0 after printing; never cache a failure.
        if v != 0.0:
            put(k, v)
        return v

    fw.safe_bert_score = cached_bert


# --------------------------------------------------------------------- setup

def setup(dry_run: bool = False):
    # Both scorer models are in the local HF cache. Online, the Hub client retries
    # a HEAD request that fails at the SSL handshake on this network, which took
    # model loading from ~20s to ~390s while changing nothing that was loaded.
    os.environ.setdefault('HF_HUB_OFFLINE', '1')
    os.environ.setdefault('TRANSFORMERS_OFFLINE', '1')
    if _EVAL not in sys.path:
        sys.path.insert(0, _EVAL)          # the framework does `from engineering_parser import`
    from dotenv import load_dotenv
    load_dotenv(os.path.join(_ROOT, '.env'))

    buf = io.StringIO()
    with redirect_stdout(buf):             # it prints while loading its scorers
        import engtrace_evaluation_framework as fw
    import google.generativeai as genai

    log = CallLog()
    log.dry_run = dry_run

    fw.MODEL_OPENAI = JUDGES['openai']['model']
    fw.MODEL_ANTHROPIC = JUDGES['anthropic']['model']
    fw.MODEL_GOOGLE = JUDGES['google']['model']

    with redirect_stdout(buf):
        # No keys handed to the constructor: the clients are injected below, so
        # it cannot build a direct client with a dead key by accident.
        framework = fw.EngTraceFramework({})
    framework.client_openai, framework.client_anthropic = build_clients(log)
    # The framework's configure() reads GOOGLE_API_KEY first when both are set;
    # GEMINI_API_KEY is the one scoped for Gemini, so pass it explicitly.
    genai.configure(api_key=os.environ['GEMINI_API_KEY'])
    wrap_google(log, genai)

    install_cache(fw, framework)

    orig_call = framework._call_single_judge

    def observed_call(provider, prompt):
        n_before = len(log.calls)
        out = orig_call(provider, prompt)
        log.prompts.append({'provider': provider, 'prompt_chars': len(prompt),
                            'n_results': len(out), 'api_calls': len(log.calls) - n_before})
        return out

    framework._call_single_judge = observed_call

    if dry_run:
        # Tier 1 and the trigger run for real; the Tribunal is replaced by a
        # recorder, so the dry run says exactly which entries would reach judges
        # and how large their prompts are, and spends nothing.
        def recorder(question, gt_steps, pred_steps, mismatch_indices):
            prompt_chars = (len(question) + len(json.dumps(gt_steps, indent=1))
                            + len(json.dumps(pred_steps, indent=1)) + 1500)
            log.prompts.append({'provider': 'DRY', 'prompt_chars': prompt_chars,
                                'mismatch_steps': len(mismatch_indices)})
            return {}, 'N/A', {}
        framework._tier2_tribunal_batch = recorder

    return SimpleNamespace(fw=fw, framework=framework, log=log, genai=genai)


# --------------------------------------------------------------------- score

def score(state, item: dict, trace: dict, seed: int) -> dict:
    log = state.log
    log.reset()
    entry = {'question': item['question'], 'solution': item['solution'],
             'generation': trace['text']}

    random.seed(seed)                      # D3
    buf = io.StringIO()
    with redirect_stdout(buf):
        out = state.framework.evaluate_entry(entry)

    # The framework prints and returns 0.0 when BERTScore raises. Surface it.
    stdout = buf.getvalue()
    scorer_failures = [ln.strip() for ln in stdout.splitlines()
                       if 'BERTScore calculation failed' in ln
                       or 'Skipping BERTScore' in ln]
    judged = [p for p in log.prompts if p['provider'] != 'DRY']
    triggered = bool(out['meta']['tribunal_triggered'])
    reached = bool(log.prompts)            # triggered AND had mismatched steps
    called = sorted({p['provider'] for p in judged})
    parsed = sorted({p['provider'] for p in judged if p['n_results'] > 0})
    # OpenRouter says "length"; the native Gemini SDK says FinishReason.MAX_TOKENS.
    truncated = [c['provider'] for c in log.calls
                 if str(c.get('finish_reason')).lower() == 'length'
                 or 'max_tokens' in str(c.get('finish_reason')).lower()]

    return {
        'scores': out['scores'],
        'meta': dict(out['meta'],
                     tribunal_reached_judges=reached,
                     judges_called=called,
                     judges_parsed=parsed,
                     judge_failures=sorted(set(called) - set(parsed)),
                     judges_truncated=truncated,
                     scorer_failures=scorer_failures,
                     framework_stdout=stdout[-2000:] or None,
                     deviations=DEVIATIONS),
        'calls': log.calls,
        'dry_prompts': [p for p in log.prompts if p['provider'] == 'DRY'],
        'triggered': triggered,
    }
