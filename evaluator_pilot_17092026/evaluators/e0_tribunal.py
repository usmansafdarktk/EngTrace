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

HOW IT RUNS FAST WITHOUT CHANGING WHAT IT COMPUTES.  The first full E0 run took
3.9 hours: 53% waiting on judge APIs one trace at a time, 47% local CPU, almost
all of it the framework's post-judgement recovery path (framework line 485)
re-scoring cross-encoder pairs one at a time.  Both are addressed OUTSIDE the
framework's own code:

  * The cross-encoder is cached per pair, so a pair is scored once ever.  Note
    that the cross-encoder is NOT batch-invariant - the same pair scored inside
    Tier 1's M x N batch and scored alone differ by up to 1.0e-6, because
    stsb-roberta-large pads every batch to its longest member - so the two call
    shapes are cached separately and a value from one is never served to the
    other.  See `install_cache`.

  * `prefetch` runs the framework's own path once with the judges stubbed,
    capturing every exact prompt it would send; fetches all of them
    concurrently; and stores each reply keyed by (slot, model, settings,
    prompt).  The scoring pass then calls the judges exactly as before and
    receives each reply from that store, in the framework's own order, on one
    thread.  `random` is still seeded per (item, model) immediately before each
    `evaluate_entry` on that one thread, so the D3 draw is a property of the
    trace and not of the schedule: --workers changes nothing about which traces
    are judged.  A prompt that is not in the store is called live, so
    correctness never depends on the prefetch - only the speed does.

  Scoring `evaluate_entry` itself concurrently is NOT safe and is not done: the
  framework's sampling reads the process-wide `random`, seeded per entry, and
  two threads seeding it between each other's draws would change which traces
  reach the Tribunal.  Measured: with any work at all between the seed and the
  draw, 300 of 300 draws differ from the serial run.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import io
import json
import os
import random
import sys
import threading
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
    # torch is deliberately NOT here. The same scoring stack runs on a CUDA build
    # on Kaggle and a CPU build locally, and a torch-in-the-key cache could never
    # be shared between them. Device and torch build are recorded per row by
    # compute() instead, and cross-device agreement is PROVED before a Kaggle
    # number is used (run_evaluator --import-kaggle), not assumed.
    for pkg in ('transformers', 'sentence-transformers', 'bert-score',
                'rouge-score', 'tokenizers', 'scipy', 'numpy'):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def compute() -> dict:
    """Where a row's numbers were produced. Provenance, not config."""
    import platform
    try:
        import torch
        cuda = torch.cuda.is_available()
        return {'torch': torch.__version__, 'cuda': cuda,
                'device': torch.cuda.get_device_name(0) if cuda else platform.processor() or 'cpu',
                'host': 'kaggle' if os.path.isdir('/kaggle') else platform.node()}
    except Exception as exc:                                       # noqa: BLE001
        return {'error': str(exc)[:200]}


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


class _Capture(Exception):
    """Raised inside a stubbed judge call. `_call_single_judge` swallows it and
    returns [], which is how a capture pass costs nothing and scores nothing."""


def reply_key(provider: str, kw: dict) -> str:
    """(slot, model, call settings, prompt) -> the reply.

    The settings are in the key because they change the reply: the framework
    sends the Anthropic slot max_tokens=2048 and the OpenAI slot nothing at all.
    """
    settings = {k: v for k, v in sorted(kw.items()) if k != 'messages'}
    prompt = kw['messages'][-1]['content']
    return hashlib.sha256(json.dumps([provider, settings, prompt], sort_keys=True).encode()).hexdigest()


class Replies:
    """(slot, model, settings, prompt) -> the judge call record. Append-only.

    The value stored IS the record that goes on the row, so a prefetched call and
    a live one put exactly the same fields in `calls` - only `prefetched` differs.
    """

    def __init__(self, path=None):
        self.path = path or REPLIES
        self.lock, self.store = threading.Lock(), {}
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        if os.path.exists(self.path):
            for ln in open(self.path, encoding='utf-8'):
                try:
                    r = json.loads(ln)
                    if r['v'].get('ok'):
                        self.store[r['k']] = r['v']
                except (json.JSONDecodeError, KeyError, AttributeError, TypeError):
                    continue

    def get(self, k):
        return self.store.get(k)

    def put(self, k, v):
        with self.lock:
            self.store[k] = v
            with open(self.path, 'a', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'k': k, 'v': v}, ensure_ascii=False) + '\n')


_CLIENT = []
_CLIENT_LOCK = threading.Lock()


def _client():
    """One OpenRouter client for the process, exactly as before.

    It is built once and shared: the SDK's httpx client is thread-safe and pools
    its connections, so the prefetch's threads do not each pay a TLS handshake -
    and a client per call would leak a pool per call.
    """
    with _CLIENT_LOCK:
        if not _CLIENT:
            from openai import OpenAI
            _CLIENT.append(OpenAI(api_key=os.environ['OPENROUTER_API_KEY'],
                                  base_url='https://openrouter.ai/api/v1',
                                  timeout=300.0, max_retries=2))
        return _CLIENT[0]


def _openrouter_call(provider: str, kw: dict):
    """One judge call, as the framework makes it. Returns (call record, exception).

    Plain HTTP and no framework state, so the prefetch can run this on a pool of
    threads while the scoring pass, which is serial, has not started. The
    exception is handed back rather than swallowed, so the caller on the framework's
    own path can re-raise the object the SDK raised and nothing about how a failed
    judge is reported changes.
    """
    cli = _client()
    rec = {'provider': provider, 'model': kw.get('model'), 'max_tokens': kw.get('max_tokens')}
    failure = None
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
    except Exception as exc:                                       # noqa: BLE001
        rec.update(ok=False, error='%s: %s' % (type(exc).__name__, str(exc)[:300]))
        failure = exc
    rec['seconds'] = round(time.time() - t0, 2)
    return rec, failure


def _as_response(rec: dict):
    """The SDK shape the framework reads back out of a stored reply."""
    return SimpleNamespace(
        model=rec.get('served_model'), model_extra={'provider': rec.get('serving_provider')},
        choices=[SimpleNamespace(message=SimpleNamespace(content=rec.get('text') or ''),
                                 finish_reason=rec.get('finish_reason'))],
        usage=SimpleNamespace(prompt_tokens=rec.get('prompt_tokens'),
                              completion_tokens=rec.get('completion_tokens')))


def _openrouter(log: CallLog, provider: str, state=None, stub: bool = False):
    def create(**kw):
        k = reply_key(provider, kw)
        if state is not None and state.capturing is not None:
            state.capturing.append({'provider': provider, 'model': kw.get('model'),
                                    'kw': {k2: v for k2, v in kw.items() if k2 != 'messages'},
                                    'prompt': kw['messages'][-1]['content'], 'key': k})
            raise _Capture()
        hit = state.replies.get(k) if state is not None else None
        if hit is not None:
            log.calls.append(dict(hit, prefetched=True))
            return _as_response(hit)
        if stub:
            raise RuntimeError('this E0 state has no judge client: capture or warm-pairs only')
        rec, failure = _openrouter_call(provider, kw)
        rec['prefetched'] = False
        log.calls.append(rec)
        if failure is not None:
            raise failure          # the framework sees exactly the SDK's own exception
        if state is not None:
            state.replies.put(k, {k2: v for k2, v in rec.items() if k2 != 'prefetched'})
        return _as_response(rec)

    return create


def build_clients(log: CallLog, state=None, stub: bool = False):
    """Objects shaped exactly like the SDK clients the framework calls."""
    if not stub:
        _client()                  # a missing key fails here, loudly, not inside a
                                   # judge call where the framework would swallow it
    or_openai = _openrouter(log, 'openai', state, stub)
    openai_client = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(create=or_openai)))

    or_anthropic = _openrouter(log, 'anthropic', state, stub)

    def anthropic_create(model, max_tokens, messages):
        # The framework reads resp.content[0].text; OpenRouter speaks chat.completions.
        r = or_anthropic(model=model, max_tokens=max_tokens, messages=messages)
        return SimpleNamespace(content=[SimpleNamespace(text=r.choices[0].message.content or '')])

    anthropic_client = SimpleNamespace(messages=SimpleNamespace(create=anthropic_create))
    return openai_client, anthropic_client


def wrap_google(log: CallLog, genai, state=None):
    """Log Gemini calls without changing how the framework makes them.

    E0 never gets here - its Google slot is gated on `genai.get_model_info`, which
    the SDK does not have (E0-F6) - but E0-3J supplies that function and does. The
    capture and reply-store branches are the same as the OpenRouter slots' so that
    a captured Gemini prompt is never silently called for real during a capture.
    """
    real = genai.GenerativeModel

    class LoggedModel:
        def __init__(self, name, *a, **kw):
            self._m = None if (state is not None and state.capturing is not None) else real(name, *a, **kw)
            self._name = name
            self._args = (a, kw)

        def generate_content(self, prompt, **kw):
            k = reply_key('google', {'model': self._name, 'generation_config': str(kw.get('generation_config')),
                                     'messages': [{'content': prompt}]})
            if state is not None and state.capturing is not None:
                state.capturing.append({'provider': 'google', 'model': self._name,
                                        'kw': {'generation_config': str(kw.get('generation_config'))},
                                        'prompt': prompt, 'key': k})
                raise _Capture()
            hit = state.replies.get(k) if state is not None else None
            if hit is not None:
                log.calls.append(dict(hit, prefetched=True))
                return SimpleNamespace(text=hit.get('text') or '', candidates=None,
                                       usage_metadata=None, model_version=hit.get('served_model'))
            if self._m is None:
                self._m = real(self._name, *self._args[0], **self._args[1])
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
                if state is not None:
                    state.replies.put(k, dict(rec, seconds=round(time.time() - t0, 2)))
                return r
            except Exception as exc:                               # noqa: BLE001
                rec.update(ok=False, error='%s: %s' % (type(exc).__name__, str(exc)[:300]))
                raise
            finally:
                rec['seconds'] = round(time.time() - t0, 2)
                rec['prefetched'] = False
                log.calls.append(rec)

    genai.GenerativeModel = LoggedModel


# --------------------------------------------------------------------- cache

CACHE_DIR = os.path.join(_PILOT, 'scores', '_cache')
CACHE = os.path.join(CACHE_DIR, 'e0_tier1.jsonl')       # legacy single file, still read
REPLIES = os.path.join(CACHE_DIR, 'e0_judge_replies.jsonl')


def install_cache(fw, framework):
    """Memoise Tier 1, the cross-encoder and BERTScore to disk, keyed by exact
    inputs and library versions.

    All three are pure functions of their arguments: the same step lists give the
    same matrix and the same texts give the same BERTScore (checked: 4 vs 8 torch
    threads, and CPU vs Kaggle T4, max score difference 0.0). On this CPU the
    cross-encoder costs 0.5-1s per step pair, so Tier 1 over 300 traces is hours -
    and a dry run followed by the real run would pay it twice for identical
    numbers. The cache changes when work happens, never what comes out; the
    library versions are in the key, so a cached value from one scorer stack is
    never served to another.

    WHY THE PER-PAIR CACHE DID NOT HELP, AND WHAT IT DOES NOW.  The framework
    scores each (gt step, pred step) pair twice: once inside Tier 1's M x N batch
    (framework line 165) and again, one pair to a call, in the post-judgement
    recovery path (line 485), where the result feeds a single `argmax`. The
    recovery path was 47% of the first E0 run's wall clock. A per-pair cache was
    added to collapse the two - but the Tier 1 MATRIX cache sits above it and
    returns V without calling the cross-encoder at all, so the batch never seeded
    the pair store and every recovery pair was still a miss. Measured on the
    pilot's own cache: 300 of 300 matrices cached, but only 60% of pairs, and the
    pairs that were there are exactly the ones the recovery path had already paid
    for. Rows that reached the judges spent 29.0s of local CPU each; rows that did
    not spent 0.6s.

    The seeding is now done deliberately, by `warm_pairs`, ahead of the run and
    free of charge - and the cache no longer pretends the two call shapes are
    interchangeable. They are not: on real pairs from this pilot the same pair
    scored inside Tier 1's batch and scored alone agreed on 13 of 45 and differed
    by up to 1.0e-6 on the rest, because stsb-roberta-large pads every batch to
    its longest member (`predict(pairs[:10])` likewise differs from
    `predict(pairs)[:10]`). So a multi-pair call goes straight through to the real
    scorer - the matrix cache above is what stops it repeating - and only the
    one-pair call shape is cached, computed exactly as line 485 computes it.

    THE SCORER IS NOT BIT-REPRODUCIBLE ACROSS PROCESSES EITHER. Re-scoring 20 of
    the pair values the September runs left on disk reproduced 3 of them exactly
    and the rest to within 3.0e-7 - a thread-count-dependent reduction order, not
    a cache defect. So "output-identical" here means what it can mean: on 14 of
    14 judged traces the recorded row's eight score fields and its whole vote
    breakdown came out the same whether line 485 read the cache or scored every
    pair again. The cache also makes a re-run MORE reproducible than no cache,
    because it pins the value instead of recomputing a slightly different one.
    """
    import numpy as np
    import torch
    # One process per model column is the parallelism here, so each process takes a
    # share of the cores rather than all of them: five processes each claiming 8
    # threads on an 8-core machine is slower than five claiming 2.
    torch.set_num_threads(int(os.environ.get('ENGTRACE_TORCH_THREADS', os.cpu_count() or 4)))

    os.makedirs(CACHE_DIR, exist_ok=True)
    store = {}
    for fn in sorted(os.listdir(CACHE_DIR)):
        # The judge reply stores live in the same directory and are also {k, v}
        # lines, so they would load into the scorer cache and sit there: same key
        # space, never a collision, just megabytes of judge text in a dict that
        # holds floats. They are read by Replies, not here.
        if not fn.endswith('.jsonl') or '_judge_replies' in fn:
            continue
        for ln in open(os.path.join(CACHE_DIR, fn), encoding='utf-8'):
            try:
                rec = json.loads(ln)
                store[rec['k']] = rec['v']
            except (json.JSONDecodeError, KeyError):
                continue
    # One part file per process: several model columns can run at once without
    # two processes interleaving writes into a single file.
    part = os.path.join(CACHE_DIR, 'part-%d.jsonl' % os.getpid())
    libs = json.dumps(libraries(), sort_keys=True)
    lock = threading.Lock()

    def key(kind, *args):
        return hashlib.sha256(json.dumps([kind, libs, args], sort_keys=True).encode()).hexdigest()

    def put(k, v):
        with lock:
            store[k] = v
            with open(part, 'a', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'k': k, 'v': v}) + '\n')

    # --- the cross-encoder itself, per PAIR ------------------------------------
    # 'pair' means one pair to a call, which is the recovery path's call shape
    # (line 485). Tier 1's batch is a different computation of the same pair - see
    # the docstring - so it is neither served from here nor stored into here.
    real_predict = fw.CROSS_ENCODER.predict

    def cached_predict(pairs, *a, **kw):
        pairs = list(pairs)
        if len(pairs) != 1:
            return real_predict(pairs, *a, **kw)
        k = key('pair', pairs[0][0], pairs[0][1])
        hit = store.get(k)
        if hit is None:
            out = real_predict(pairs, *a, **kw)
            put(k, float(out[0]))
            return out                        # a miss returns the scorer's own array
        return np.array([hit], dtype=np.float32)

    fw.CROSS_ENCODER.predict = cached_predict

    real_matrix = framework._tier1_verify_matrix
    handle = SimpleNamespace(store=store, key=key, put=put, predict=cached_predict,
                             real_predict=real_predict, part=part, entry=None)

    def cached_matrix(gt_s_txt, gt_s_val, pred_s_txt, pred_s_val):
        # The only place the entry's two step lists are visible. `warm_pairs` reads
        # them from here, because the pairs the recovery path will ask for are a
        # subset of exactly these and nothing downstream carries them.
        handle.entry = (list(gt_s_txt), list(pred_s_txt))
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
    return handle


# --------------------------------------------------------------------- setup

def setup(dry_run: bool = False, keys: bool = True):
    """keys=False builds the framework without any API client. E1 uses it on a
    Kaggle kernel, which holds no key: its judges are served from a reply store
    fetched on the laptop, where the key lives (D-086)."""
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
    state = SimpleNamespace(fw=fw, framework=framework, log=log, genai=genai,
                            replies=Replies(), capturing=None, keys=bool(keys and not dry_run))
    if not dry_run and keys:
        # A dry run never reaches a judge, so it needs no keys - which is what lets
        # it run on a Kaggle kernel that has none and never should.
        framework.client_openai, framework.client_anthropic = build_clients(log, state)
        # The framework's configure() reads GOOGLE_API_KEY first when both are set;
        # GEMINI_API_KEY is the one scoped for Gemini, so pass it explicitly.
        genai.configure(api_key=os.environ['GEMINI_API_KEY'])
    wrap_google(log, genai, state)

    state.cache = install_cache(fw, framework)

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

    state.compute = dict(compute(), framework_device=str(fw.device))
    return state


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
                     compute=state.compute,
                     framework_stdout=stdout[-2000:] or None,
                     deviations=DEVIATIONS),
        'calls': log.calls,
        'dry_prompts': [p for p in log.prompts if p['provider'] == 'DRY'],
        'triggered': triggered,
    }


# ------------------------------------------------------- capture, warm, prefetch

def capture(state, jobs) -> list[dict]:
    """Every judge prompt the real run will send, and the step texts behind it.

    `jobs` are (item, trace, seed). `random` is seeded and the framework called
    exactly as `score` does, so the same traces trigger the Tribunal and build the
    same prompts; the judges raise before any HTTP, so nothing is spent, nothing is
    scored and nothing is written. The recovery path is not reached either - with
    no judge votes there is nothing to recover - so this pass is Tier 1 only, and
    Tier 1 is already cached.

    Needs no key: the framework decides WHICH judges to call by whether a client
    object exists, not by whether it works, so a keyless state gets stubs.
    """
    fw_ = state.framework
    stubbed = fw_.client_openai is None or fw_.client_anthropic is None
    if stubbed:
        fw_.client_openai, fw_.client_anthropic = build_clients(state.log, state, stub=True)
    state.capturing = []
    out, buf = [], io.StringIO()
    try:
        for item, trace, seed in jobs:
            n_before = len(state.capturing)
            state.cache.entry = None
            random.seed(seed)                          # D3, exactly as `score` does
            with redirect_stdout(buf):
                state.framework.evaluate_entry({'question': item['question'],
                                                'solution': item['solution'],
                                                'generation': trace['text']})
            for c in state.capturing[n_before:]:
                c['entry'] = state.cache.entry         # (gt step texts, pred step texts)
                out.append(c)
    finally:
        state.capturing = None
        state.log.reset()                              # the stub "failures" are not a row's calls
        if stubbed:
            # A stub must never be left where a scoring pass could reach it: a
            # RuntimeError there is swallowed by the framework and looks exactly
            # like a judge that voted wrong, which is E0-F6 all over again.
            fw_.client_openai = fw_.client_anthropic = None
    return out


def warm_pairs(state, captured: list[dict], every: int = 100) -> dict:
    """Score every cross-encoder pair the recovery path can ask for, into the cache.

    This is the work framework line 485 does inline, one pair to a call, while the
    judges' replies sit waiting - 47% of the first E0 run's wall clock. It is the
    same computation, in the same call shape, so the numbers are the ones the
    framework would have produced; only WHEN it happens changes. It costs nothing
    but local CPU and can run before a paid run is ever approved.

    The recovery path asks for (gt_i, pred_j) for every i and for each pred step a
    judge recovered, so the pairs it can ask for are exactly this entry's M x N.
    """
    cache = state.cache
    todo, seen = [], set()
    for c in captured:
        gt, pred = c.get('entry') or ([], [])
        for a in gt:
            for b in pred:
                k = cache.key('pair', a, b)
                if k not in cache.store and k not in seen:
                    seen.add(k)
                    todo.append((a, b))
    t0 = time.time()
    for n, p in enumerate(todo, 1):
        cache.predict([p])                             # one pair to a call: line 485's shape
        if n % every == 0 or n == len(todo):
            print('  warmed %d/%d cross-encoder pairs, %.0fs' % (n, len(todo), time.time() - t0),
                  flush=True)
    return {'pairs_computed': len(todo), 'seconds': round(time.time() - t0, 1)}


def fetch_captured(captured: list[dict], workers: int = 16) -> dict:
    """Fetch every captured prompt's reply concurrently into the reply store.

    Plain HTTP: it imports neither torch nor the framework. Judge calls are the
    only thing here that runs on more than one thread, and they share nothing with
    the scoring pass but an append-only store keyed by the prompt.
    """
    replies = Replies()
    todo = {}
    skipped = 0
    for c in captured:
        if c['provider'] == 'google':
            skipped += 1                               # native SDK, not OpenRouter; left to the run
            continue
        if replies.get(c['key']) is None:
            todo[c['key']] = c
    print('  %d captured prompts, %d already stored, %d to fetch with %d workers%s'
          % (len(captured), len(captured) - len(todo) - skipped, len(todo), workers,
             ', %d Gemini prompts left to the scoring pass' % skipped if skipped else ''),
          flush=True)
    done = failed = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_openrouter_call, c['provider'],
                            dict(c['kw'], messages=[{'role': 'user', 'content': c['prompt']}])): k
                for k, c in todo.items()}
        for fut in cf.as_completed(futs):
            rec, _exc = fut.result()
            if rec['ok']:
                replies.put(futs[fut], rec)
                done += 1
            else:
                failed += 1
                print('    FAILED %s: %s' % (todo[futs[fut]]['model'], rec.get('error')), flush=True)
            n = done + failed
            if n % 25 == 0 or n == len(futs):
                print('  %d/%d fetched, %d failed, %.0fs' % (n, len(futs), failed, time.time() - t0),
                      flush=True)
    return {'captured': len(captured), 'fetched': done, 'failed': failed,
            'stored_total': len(replies.store)}


def prefetch(state, jobs, workers: int = 16):
    """Capture the prompts, fetch the replies concurrently, warm the pairs.

    The judge fetches are HTTP and the pair warming is torch, so the warming runs
    on this thread while the pool waits on the network: the two halves of the old
    run's wall clock overlap instead of following each other. Neither touches the
    framework, and the scoring pass that follows is still strictly serial.
    """
    caps = capture(state, jobs)
    # One prompt per trace, sent to every connected slot, so distinct prompts are
    # the traces that reach the Tribunal.
    print('  captured %d judge prompts from %d of %d traces that reach the Tribunal'
          % (len(caps), len({c['prompt'] for c in caps}), len(jobs)), flush=True)
    todo = [c for c in caps if c['provider'] != 'google' and state.replies.get(c['key']) is None]
    out = {}
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_openrouter_call, c['provider'],
                            dict(c['kw'], messages=[{'role': 'user', 'content': c['prompt']}])): c
                for c in todo}
        out['warm'] = warm_pairs(state, caps)          # torch, while the pool waits on HTTP
        done = failed = 0
        for fut in cf.as_completed(futs):
            rec, _exc = fut.result()
            if rec['ok']:
                state.replies.put(futs[fut]['key'], rec)
                done += 1
            else:
                failed += 1
                print('    FAILED %s: %s' % (futs[fut]['model'], rec.get('error')), flush=True)
    out.update(captured=len(caps), fetched=done, failed=failed,
               already_stored=len(caps) - len(todo))
    return out
