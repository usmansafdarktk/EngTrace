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
  reasoning-<effort>  (C1, docs/EVALUATION_NEXT_STEPS.md) the 450 items of subsamples.PARAPHRASE with their original
                    questions, the request carrying OpenRouter's unified reasoning parameter at that effort (low,
                    medium or high), for the models whose endpoints reported no reasoning tokens at the provider's
                    default (D-168): by default `gpt-5.4-mini` and `gemini-3.1-flash-lite`, or the models --model
                    names. Everything else is the main run's: prompt, ceiling, routing, scoring. Reported beside the
                    main run on the same items, never in its place.
  flagship          (C3) the same 450 items for the anchor models of models.json (`anchor: true`, inert for every
                    other mode): flagships that pass the roster rule, reported beside the roster on those items as an
                    anchor with 450-item intervals, outside the pairwise family. Priced from the pricing basis like the
                    main run (they have no main-run bills); `--calibrate N` measures their lengths first.
  openbook          (C4, D-183) the 450-item subsample's questions with the template's governing equations appended
                    (openbook.py --build: openbook/items.jsonl, local, checked against the committed
                    openbook/manifest.jsonl before a call), for the models --model names or, by default,
                    `claude-sonnet-5`, `gpt-5.4-mini` and `gpt-oss-20b`; the row records the modified question's
                    hash as item_sha256 and the original's as original_sha256, as the paraphrase arm does.
  tool              (C4, D-184) the 450-item subsample's original questions with a Python tool offered in the request
                    (`tools`: one function, `python(code)`; `tool_choice: auto`), for the open-book arm's three models
                    by default. The prompt is the main run's, word for word: nothing tells the model to use the tool
                    beyond the tool's own description. When the model calls it, the script runs in a fresh isolated
                    interpreter (no inherited environment, a scratch directory, a wall-clock limit, the output
                    truncated; a static filter refuses file, network, process and introspection access) and the
                    output goes back as a tool message; the loop ends when a turn carries no tool call, or after
                    TOOL_MAX_CALLS calls, when one more turn is asked with the tool withheld so the trace ends in an
                    answer. The row's text is the transcript the scorer reads (every assistant text, each script and
                    its output in fenced blocks, the final answer); final_text, turns, tool_calls, tool_errors,
                    tool_refused, tool_limit and tool_turns record the tool use; tokens and the bill are summed over
                    the turns. `--selftest` exercises the sandbox and the loop offline, for free.
  reasoning-medium-full  (the matched configuration) all 2,250 pool items with their original questions and the
                    reasoning parameter at effort medium, as the 450-item arm `reasoning-medium` sent it, for the models
                    that returned no reasoning tokens at their provider's default and whose endpoints honour the
                    parameter: MATCHED_MODELS, or the models --model names (a calibration before a model joins).
                    Everything else is the main run's: the prompt (--check and every launch refuse unless its hash
                    equals the one on the main run's rows), the 32,768 ceiling and the routing.
  paraphrase-reasoning-medium   the cascade's paraphrase pairs: the passing paraphrases less the pairs the expert
                    rejected (paraphrase/accepted.json, the pairs analyze.py keeps), with the parameter on, for
                    MATCHED_MODELS; rows record the paraphrase's hash and the original's as the paraphrase arm's do
  repeat1-reasoning-medium .. repeat3-reasoning-medium   the cascade's repeats: the 300 repeat items with the
                    parameter on, for the MATCHED_MODELS that are in the decoding repeats' set (REPEAT_MODELS, D-151)
The dry run estimates a variant from each model's own bills on the same items in the main run; for a reasoning
variant the main run's visible output understates the bill, so its dry run prices output-token multipliers and
`--calibrate N` (allowed for this variant) measures the real lengths on N items first. The matched-configuration
variants are priced from bills measured with the parameter on, never the main run's: reasoning-medium-full from the
450-item arm's rows (the rows decoding_table_reasoning-medium.json sums), or from its own calibration rows for a
model the arm did not run; the cascade from each model's reasoning-medium-full rows on the same items.

--only-items PATH restricts every mode of every variant to the item ids PATH lists, one per line: the resume of a
subset whose rows were archived (re-drawn items). An id that is no pool item stops the run; an id outside the
variant's items is skipped and counted. --review (FREE) prints the checks a finished run must pass before scoring:
missing items, finish reasons, empty rows at the ceiling, the rows with reasoning tokens and with reasoning text, the
providers, the request sets and the prompt hash against the main run's. --matched-config (FREE) writes
results/matched_config.json, the setting and stores per roster model, once every store it names is complete.
--trace-review (FREE) runs trace_review.py on one of these variants with the variant's item set, which that script's own
reading of variant names does not give. In these variants a completion that comes back without a finish reason is a
provider fault, asked again like finish_reason=error (Io Net cut about 1% of Gemma 4's completions so, mid-sentence);
--requeue-unfinished (FREE) moves rows written so before the rule to traces/_unfinished/, so the next run asks again.

    python -m full_run_28092026.run_traces --variant reasoning-medium-full --model gemma-4-26b-a4b --calibrate 20 --yes
    python -m full_run_28092026.run_traces --variant reasoning-medium-full --dry-run                    # FREE
    python -m full_run_28092026.run_traces --variant reasoning-medium-full --model gpt-5.4-mini --yes
    python -m full_run_28092026.run_traces --variant reasoning-medium-full --review                     # FREE
    python -m full_run_28092026.run_traces --variant repeat1 --model gemma-4-26b-a4b --only-items ids.txt --yes

    python -m full_run_28092026.run_traces --variant paraphrase --dry-run          # FREE
    python -m full_run_28092026.run_traces --variant repeat1 --model gemma-4-26b-a4b --dry-run
    python -m full_run_28092026.run_traces --selftest                              # FREE: the tool arm, offline
    python -m full_run_28092026.run_traces --variant tool --dry-run                # FREE: multipliers on the main bills
    python -m full_run_28092026.run_traces --variant tool --calibrate 20 --yes     # BILLS: measured turns and lengths
    python -m full_run_28092026.run_traces --variant tool --model gpt-oss-20b --ignore-provider Darkbloom --yes   # a provider that drops tool calls skipped
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
REASONING_EFFORTS = ('low', 'medium', 'high')
REASONING_VARIANTS = tuple(f'reasoning-{e}' for e in REASONING_EFFORTS)
REASONING_MODELS = ('gpt-5.4-mini', 'gemini-3.1-flash-lite')     # no reasoning tokens at the provider's default (D-168)
FLAGSHIP_REASONING = tuple(f'flagship-reasoning-{e}' for e in REASONING_EFFORTS)   # the closed anchor with reasoning on (D-182)
TOOL = 'tool'                                                                     # C4's tool condition (D-184)
MATCHED_EFFORT = 'medium'                   # the matched configuration's effort: the 450-item arm's (D-180)
FULL_REASONING = 'reasoning-medium-full'    # all 2,250 items with the parameter on
CASCADE_PARAPHRASE = 'paraphrase-reasoning-medium'
CASCADE_REPEATS = tuple(f'{v}-reasoning-medium' for v in subsamples.REPEAT_VARIANTS)
MATCHED_VARIANTS = (FULL_REASONING, CASCADE_PARAPHRASE) + CASCADE_REPEATS
MATCHED_MODELS = ('gpt-5.4-mini', 'gemini-3.1-flash-lite', 'gemma-4-26b-a4b')   # reasoning on over the full set (D3); Gemma 4 by its calibration: 20 of 20 rows reasoned
REPEAT_MODELS = ('gemma-4-26b-a4b', 'gpt-oss-20b', 'qwen3-235b-a22b-2507', 'gemini-3.1-flash-lite')   # the decoding repeats' four (D-151)
MATCHED_CONFIG = HERE / 'results' / 'matched_config.json'
NO_SETTING = ('qwen3-235b-a22b-2507',)      # roster models whose endpoints offer no reasoning setting
CARD_NOTES = {      # what a model's card and OpenRouter's catalogue say of its thinking mode (read 2026-10-06)
    'qwen3-235b-a22b-2507': 'its card (Hugging Face Qwen/Qwen3-235B-A22B-Instruct-2507) says it "supports only non-thinking '
                            'mode", and OpenRouter lists no reasoning parameter for it',
    'gemma-4-26b-a4b': 'its card (Hugging Face google/gemma-4-26B-A4B-it) names "configurable thinking modes", switched on by '
                       'enable_thinking; OpenRouter lists the reasoning parameter for its endpoints',
    'qwen3-235b-a22b-thinking-2507': 'its card (Hugging Face Qwen/Qwen3-235B-A22B-Thinking-2507) says it "supports only '
                                     'thinking mode"',
}
VARIANTS =('main', 'paraphrase') + subsamples.REPEAT_VARIANTS + REASONING_VARIANTS + ('flagship', 'openbook', 'openbook2') + FLAGSHIP_REASONING + (TOOL,) + MATCHED_VARIANTS
OPENBOOK = HERE / 'openbook'
OPENBOOK_MODELS = ('claude-sonnet-5', 'gpt-5.4-mini', 'gpt-oss-20b')   # one model from each tier (D-183)
TOOL_MODELS = OPENBOOK_MODELS               # the tool arm runs the open-book arm's three, for a like-for-like reading (D-184)
TOOL_MAX_CALLS = 8                          # tool calls per item before the model is asked to answer with the tool withheld
TOOL_TIMEOUT = 20                           # seconds of wall clock per script
TOOL_OUTPUT_CHARS = 4000                    # what the model gets back from one script
PYTHON_TOOL = {'type': 'function', 'function': {
    'name': 'python',
    'description': 'Run a Python 3 script for arithmetic and numerical computation (the standard library, numpy, scipy '
                   'and sympy are available). The script runs in a fresh process each time: print() every value you '
                   'need to see. No file, network or system access.',
    'parameters': {'type': 'object',
                   'properties': {'code': {'type': 'string', 'description': 'The Python source to run.'}},
                   'required': ['code']}}}
TOOL_FORBIDDEN = re.compile(r'''(?x)
    ^\s*(?:import|from)\s+(?:os|sys|subprocess|socket|shutil|pathlib|ctypes|multiprocessing|threading|signal|importlib|
        urllib|requests|http|ftplib|smtplib|telnetlib|webbrowser|pickle|marshal|builtins|code|pty|resource|tempfile|
        glob|io|platform|getpass|sqlite3|asyncio|concurrent|xmlrpc|ssl|select|mmap|faulthandler|gc|inspect|site|
        sysconfig|zipimport|runpy|pkgutil|setuptools|pip|winreg|msvcrt|_thread|codecs|fileinput|shelve|dbm)\b
  | \b(?:open|exec|eval|compile|__import__|input|breakpoint|globals|locals|vars|getattr|setattr|delattr|memoryview)\s*\(
  | __(?:builtins|subclasses|class|bases|mro|globals|dict|loader|spec|code|closure)__
''', re.M)
DEPLOYED_RUNNER = REPO / 'evaluation' / 'run_inference.py'
PILOT_PROMPT_PREFIX = 'c2bcb87984c4e50b'   # evaluator_pilot_17092026/models.json: the annotated traces' prompt
MAX_ATTEMPTS = 4
BACKOFF = 5                                 # seconds, doubled per attempt
DONE = ('answered', 'empty')                # states a re-run does not call again
PROVIDER_IGNORE: list[str] = []              # --ignore-provider: providers routing must skip in this invocation (recorded in each row's request)


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
    if variant in ('main', FULL_REASONING):
        return its
    by_id = {it['item_id']: it for it in its}
    if variant in subsamples.REPEAT_VARIANTS or variant in CASCADE_REPEATS:
        return [by_id[i] for i in subsamples.repeat_ids()]
    if variant == CASCADE_PARAPHRASE:
        keep = kept_pairs()
        return [it for it in variant_items('paraphrase', its) if it['item_id'] in keep]
    if variant in REASONING_VARIANTS or variant == 'flagship' or variant in FLAGSHIP_REASONING or variant == TOOL:
        return [by_id[i] for i in subsamples.paraphrase_ids()]       # the originals of the 450-item subsample
    if variant in ('openbook', 'openbook2'):
        sfx = '' if variant == 'openbook' else '2'
        pool, manifest = OPENBOOK / f'items{sfx}.jsonl', OPENBOOK / f'manifest{sfx}.jsonl'
        if not pool.exists() or not manifest.exists():
            raise SystemExit(f'openbook/items{sfx}.jsonl or its manifest is missing: run openbook.py --build' + (' --version 2' if sfx else '') + ' first')
        want = {r['item_id']: r for r in map(json.loads, manifest.read_text(encoding='utf-8').splitlines())}
        out = []
        for r in map(json.loads, pool.read_text(encoding='utf-8').splitlines()):
            m = want.get(r['item_id'])
            if m is None:
                continue
            sha = hashlib.sha256(r['question'].encode('utf-8')).hexdigest()
            if sha != m['sha256'] or by_id[r['item_id']]['sha256'] != m['original_sha256']:
                raise SystemExit(f"openbook item {r['item_id']} does not match its manifest: rebuild it")
            out.append({**by_id[r['item_id']], 'question': r['question'], 'sha256': sha,
                        'original_sha256': m['original_sha256']})
        if len(out) != len(want):
            raise SystemExit(f'the openbook manifest lists {len(want)} items but the pool holds {len(out)}')
        return out
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


def kept_pairs() -> set[str]:
    """The paraphrase pairs the analysis keeps: every passing pair the expert did not reject (paraphrase/accepted.json,
    the rule of analyze.accepted_pairs)."""
    p = PARAPHRASES / 'accepted.json'
    if not p.exists():
        raise SystemExit('paraphrase/accepted.json is missing: the expert check has not been scored (paraphrase_kit.py --score)')
    return {i for i, r in json.loads(p.read_text(encoding='utf-8')).items() if r.get('kept') is not False}


def only_items(its: list[dict], pool: list[dict], path: str) -> list[dict]:
    """--only-items: the variant's items whose ids PATH lists, one per line (the line's first word; blank lines and lines
    starting with # skipped, since an id holds a # itself). An id that is no pool item stops the run (a typo, or an id a
    re-draw retired); an id outside the variant's items is skipped and counted."""
    words = [ln.split() for ln in Path(path).read_text(encoding='utf-8').splitlines()]
    ids = {w[0] for w in words if w and not w[0].startswith('#')}
    unknown = sorted(ids - {it['item_id'] for it in pool})
    if unknown:
        raise SystemExit(f'--only-items: {len(unknown)} of the {len(ids)} ids are not pool items (first: {unknown[:3]})')
    out = [it for it in its if it['item_id'] in ids]
    print(f'--only-items {path}: {len(out)} of its {len(ids)} ids are items of this variant'
          + (f'; {len(ids) - len(out)} are not and are skipped' if len(out) < len(ids) else ''))
    return out


def main_prompt_hashes() -> Counter:
    """The prompt hash on every row of the main run's trace files, counted: one value when the run used one prompt."""
    c = Counter()
    for p in sorted(TRACES.glob('*.jsonl')):
        with open(p, encoding='utf-8') as fh:
            for ln in fh:
                k = ln.find('"prompt_sha256": "')
                if k >= 0:
                    c[ln[k + 18:k + 82]] += 1
    return c


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


def request_params(spec: dict, cfg: dict, variant: str = 'main') -> dict:
    p = {'max_tokens': spec.get('max_tokens', cfg['max_tokens'])}
    key = 'provider_open_weight' if spec['weights'] == 'open' else 'provider_closed_weight'
    p['extra_body'] = {'provider': {**cfg[key], **({'ignore': list(PROVIDER_IGNORE)} if PROVIDER_IGNORE else {})}}
    if variant in REASONING_VARIANTS or variant in FLAGSHIP_REASONING:
        # OpenRouter's unified reasoning parameter; the provider maps the effort to its own setting. The main run
        # sent nothing here and ran at each provider's default (D-122).
        p['extra_body']['reasoning'] = {'effort': variant.rsplit('-', 1)[1]}
    if variant in MATCHED_VARIANTS:
        p['extra_body']['reasoning'] = {'effort': MATCHED_EFFORT}      # as the 450-item arm sent it
    if variant == TOOL:
        # The python tool in the request, the model free to use it or not; everything else as the main run (D-184).
        p['tools'] = [PYTHON_TOOL]
        p['tool_choice'] = 'auto'
    return p


def client(cfg: dict):
    from openai import OpenAI
    key = os.getenv(cfg['route']['key_env'])
    if not key:
        raise SystemExit(f"{cfg['route']['key_env']} is not set in .env")
    return OpenAI(api_key=key, base_url=cfg['route']['base_url'], timeout=600.0, max_retries=0)


def call(cli, spec: dict, cfg: dict, question: str, attempts: int = MAX_ATTEMPTS, variant: str = 'main') -> dict:
    """One completion with retries, ending in answered / empty / service_failure."""
    prompt = PROMPT.format(question=question)
    params = request_params(spec, cfg, variant)
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
            if r.choices[0].finish_reason is None and variant in MATCHED_VARIANTS:
                # No finish reason, the text cut off mid-sentence well below the ceiling: Io Net ended about 1% of
                # Gemma 4's completions so with reasoning on. A provider fault, asked again like the one above; the
                # other variants keep the main run's reading (one such row in its 26,970).
                raise RuntimeError('provider returned no finish reason')
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


def sandbox_env(work: str) -> dict:
    """The environment a model's script runs in: what Python needs on this OS and nothing of ours (no API key)."""
    keep = ('SYSTEMROOT', 'SYSTEMDRIVE', 'WINDIR', 'COMSPEC', 'PATHEXT', 'HOME', 'LANG', 'LC_ALL')
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env.update({'PATH': str(Path(sys.executable).parent), 'TEMP': work, 'TMP': work, 'TMPDIR': work,
                'PYTHONIOENCODING': 'utf-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1', 'MPLBACKEND': 'Agg'})
    return env


def run_python(code: str, timeout: int = TOOL_TIMEOUT) -> dict:
    """Run one script from the model in a fresh, isolated interpreter (-I: no inherited PYTHON* settings, no user
    site, no script directory on the path), in a scratch working directory removed afterwards, with a wall-clock
    limit and the output truncated. A static filter refuses scripts that reach for the file system, the network,
    processes or introspection; the model sees the refusal as the tool's output and can write the computation
    another way. The static filter is a guard for a research harness, not a security boundary."""
    import shutil
    import subprocess
    import tempfile
    t0 = time.time()
    m = TOOL_FORBIDDEN.search(code or '')
    if m:
        return {'ok': False, 'refused': True, 'truncated': False, 'seconds': 0.0,
                'output': f'ToolError: {m.group(0).strip()!r} is not available in this tool (numerical computation only).'}
    work = tempfile.mkdtemp(prefix='engtrace_tool_')
    try:
        try:
            r = subprocess.run([sys.executable, '-I', '-c', code], cwd=work, env=sandbox_env(work), capture_output=True,
                               text=True, encoding='utf-8', errors='replace', timeout=timeout, stdin=subprocess.DEVNULL)
            out = (r.stdout or '') + (('\n' + r.stderr) if r.stderr else '')
            ok = r.returncode == 0
        except subprocess.TimeoutExpired as exc:
            got = exc.stdout or ''
            out = (got if isinstance(got, str) else got.decode('utf-8', 'replace')) + \
                f'\nTimeoutError: the script ran longer than {timeout} s and was stopped.'
            ok = False
    finally:
        shutil.rmtree(work, ignore_errors=True)
    out = out.strip() or '(no output: print() the values you need)'
    truncated = len(out) > TOOL_OUTPUT_CHARS
    if truncated:
        out = out[:TOOL_OUTPUT_CHARS] + f'\n... [output truncated at {TOOL_OUTPUT_CHARS} characters]'
    return {'ok': ok, 'refused': False, 'truncated': truncated, 'seconds': round(time.time() - t0, 2), 'output': out}


def tool_arguments(raw: str):
    """The tool call's arguments as a dict, and whether they needed repair. Some endpoints return them as JSON with
    control characters inside the string (gpt-oss-20b on DekaLLM wrote the code with literal newlines), which strict
    JSON rejects; the second reading allows them, and a third takes the whole text as the code."""
    try:
        return json.loads(raw), False
    except (json.JSONDecodeError, TypeError):
        pass
    try:
        return json.loads(raw, strict=False), True
    except (json.JSONDecodeError, TypeError):
        return {'code': raw or ''}, True


def tool_call_dict(tc, arguments: str | None = None) -> dict:
    """The assistant's tool call as the provider expects it back; `arguments` replaces the raw text with valid JSON
    (a provider refuses the whole conversation otherwise: "function.arguments must be valid JSON")."""
    d = tc.model_dump(exclude_none=True) if hasattr(tc, 'model_dump') else \
        {'id': getattr(tc, 'id', None), 'type': 'function', 'function': {'name': tc.function.name, 'arguments': tc.function.arguments}}
    if arguments is not None:
        d.setdefault('function', {})['arguments'] = arguments
    return d


def call_with_tool(cli, spec: dict, cfg: dict, question: str, attempts: int = MAX_ATTEMPTS, variant: str = TOOL) -> dict:
    """The tool arm's loop (C4, D-184). Each model turn is one completion with call()'s retries; when the turn
    carries tool calls, each script runs and its output goes back as a tool message, and the model is asked again;
    the loop ends when a turn carries no tool call, or after TOOL_MAX_CALLS calls, when one more turn is asked
    with tool_choice 'none' so the trace ends in an answer. An empty final turn that did not hit the ceiling is
    asked again, as in call(). The transcript (every assistant text, each script and its output in fenced
    blocks, the final answer) is the row's text, what the scorer reads; the last message alone is final_text."""
    prompt = PROMPT.format(question=question)
    params = request_params(spec, cfg, variant)
    messages = [{'role': 'user', 'content': prompt}]
    parts, log = [], []
    tot = Counter()
    served = provider = finish = None
    msg, extra, text = None, {}, ''
    calls = attempts_used = empty_retries = 0
    limit_hit = False
    last_err = None
    t_start = time.time()
    while True:
        p = dict(params)
        if calls >= TOOL_MAX_CALLS:
            p['tool_choice'] = 'none'
            limit_hit = True
        r = None
        for attempt in range(1, attempts + 1):
            attempts_used += 1
            try:
                r = cli.chat.completions.create(model=spec['model'], messages=messages, **p)
                if r.choices[0].finish_reason == 'error':
                    raise RuntimeError('provider reported finish_reason=error')
                break
            except Exception as exc:                               # noqa: BLE001
                last_err = f'{type(exc).__name__}: {exc}'[:400]
                r = None
                if attempt < attempts:
                    time.sleep(BACKOFF * (2 ** (attempt - 1)))
        if r is None:
            return {'status': 'service_failure', 'error': last_err, 'attempts': attempts_used, 'tool_calls': calls}
        u = r.usage
        det = getattr(u, 'completion_tokens_details', None)
        tot['prompt_tokens'] += getattr(u, 'prompt_tokens', 0) or 0
        tot['completion_tokens'] += getattr(u, 'completion_tokens', 0) or 0
        rt = getattr(det, 'reasoning_tokens', None) if det else None
        if rt is not None:
            tot['reasoning_tokens'] += rt
            tot['reasoning_reported'] += 1
        tot['billed_usd'] += (getattr(u, 'model_extra', None) or {}).get('cost') or 0.0
        tot['turns'] += 1
        served = getattr(r, 'model', None) or served
        provider = (getattr(r, 'model_extra', None) or {}).get('provider') or provider
        msg = r.choices[0].message
        extra = getattr(msg, 'model_extra', None) or {}
        finish = r.choices[0].finish_reason
        text = msg.content or ''
        tool_calls = list(getattr(msg, 'tool_calls', None) or [])
        if text.strip():
            parts.append(text.strip())
        if not tool_calls or limit_hit:
            if not text.strip() and finish != 'length' and empty_retries < attempts - 1:
                empty_retries += 1                                 # an empty completion is asked again before it counts
                continue
            break
        parsed = []
        for tc in tool_calls:
            fn = getattr(tc, 'function', None)
            args, repaired = tool_arguments(getattr(fn, 'arguments', None) or '{}')
            if not isinstance(args, dict):
                args, repaired = {'code': str(args)}, True
            parsed.append((tc, getattr(fn, 'name', None), args, repaired))
        assistant = {'role': 'assistant', 'content': msg.content or None,
                     'tool_calls': [tool_call_dict(tc, json.dumps(args)) for tc, _n, args, _r in parsed]}
        if extra.get('reasoning_details'):
            assistant['reasoning_details'] = extra['reasoning_details']   # OpenRouter: the reasoning carried across turns
        messages.append(assistant)
        for tc, name, args, repaired in parsed:
            calls += 1
            if name != 'python':
                code, res = '', {'ok': False, 'refused': True, 'truncated': False, 'seconds': 0.0,
                                 'output': f'ToolError: there is no tool named {name!r}; the only tool is python.'}
            else:
                code = args.get('code') or ''
                code = code if isinstance(code, str) else json.dumps(code)
                res = run_python(code)
            log.append({'call': calls, 'code': code, 'args_repaired': repaired, **res})
            parts.append(f'```python\n{code.strip()}\n```\nOutput:\n```\n{res["output"]}\n```')
            messages.append({'role': 'tool', 'tool_call_id': getattr(tc, 'id', None) or f'call_{calls}', 'content': res['output']})
    final_text = text.strip()
    row = {
        'text': '\n\n'.join(parts), 'final_text': final_text,
        'reasoning': getattr(msg, 'reasoning', None) or extra.get('reasoning') or '',
        'served_model': served, 'provider': provider, 'finish_reason': finish,
        'prompt_tokens': tot['prompt_tokens'], 'completion_tokens': tot['completion_tokens'],
        'reasoning_tokens': tot['reasoning_tokens'] if tot['reasoning_reported'] else None,
        'billed_usd': round(tot['billed_usd'], 6), 'seconds': round(time.time() - t_start, 2), 'attempts': attempts_used,
        'turns': tot['turns'], 'tool_calls': calls, 'tool_limit': limit_hit,
        'tool_errors': sum(not t['ok'] for t in log), 'tool_refused': sum(t['refused'] for t in log), 'tool_turns': log,
    }
    return {'status': 'answered' if final_text else 'empty', **row}


def run_model(spec, todo, cfg, workers, label='run', variant='main'):
    key = spec['key']
    cli = client(cfg)
    trace_path(key, variant).parent.mkdir(parents=True, exist_ok=True)
    params = request_params(spec, cfg, variant)
    n, counts, billed = 0, Counter(), 0.0
    with open(trace_path(key, variant), 'a', encoding='utf-8', newline='\n') as fh, \
            cf.ThreadPoolExecutor(max_workers=workers) as pool:
        fn = call_with_tool if variant == TOOL else call
        futs = {pool.submit(fn, cli, spec, cfg, it['question'], MAX_ATTEMPTS, variant): it for it in todo}
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


def dry_run(cfg, its, specs, variant='main') -> int:
    print(f'{"pool" if variant == "main" else "variant " + variant}: {len(its)} items, {len({i["template_id"] for i in its})} templates, '
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
        done = existing(s['key'], variant)
        todo = sum(it['item_id'] not in done for it in its)                  # the items listed, under --only-items
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


def reasoning_dry_run(cfg, its, specs, variant) -> int:
    """A reasoning variant's estimate: the main run's prompt and visible output tokens on the same items, priced at
    the pricing document's rates with the output multiplied, because a model that reasons bills its reasoning as
    output. The multipliers bracket what is known: the pilot's GPT-5 wrote about 12 completion tokens for every
    visible one (README, stage 1); nothing here is a measurement of these models, which is what --calibrate is for."""
    effort = variant.rsplit('-', 1)[1]
    base = 'flagship' if variant in FLAGSHIP_REASONING else 'main'      # the arm whose bills price this one
    mults = (1, 3, 6, 12)
    print(f'variant {variant}: {len(its)} items, {len({i["template_id"] for i in its})} templates; reasoning effort '
          f'{effort!r} in the request; traces go to traces/{variant}/; priced from the {base} run\'s rows')
    print(f'{"model":22s} {"to run":>7s} {"main $":>8s} {"main out tok":>12s} ' + ' '.join(f'{"x" + str(m) + " $":>8s}' for m in mults))
    tot = Counter()
    for s in specs:
        todo = [it for it in its if it['item_id'] not in existing(s['key'], variant)]
        main = existing(s['key'], base)
        rows = [main[it['item_id']] for it in todo if it['item_id'] in main]
        pt = sum(r.get('prompt_tokens') or 0 for r in rows)
        ct = sum(r.get('completion_tokens') or 0 for r in rows)
        bill = sum(r.get('billed_usd') or 0.0 for r in rows)
        price = s['price_per_m']
        ests = {m: (pt * price['in'] + ct * m * price['out']) / 1e6 for m in mults}
        for m in mults:
            tot[m] += ests[m]
        tot['main'] += bill
        print(f'{s["key"]:22s} {len(todo):7d} {bill:8.3f} {ct // max(len(rows), 1):12d} ' + ' '.join(f'{ests[m]:8.2f}' for m in mults))
    print(f'{"TOTAL":22s} {"":7s} {tot["main"]:8.3f} {"":12s} ' + ' '.join(f'{tot[m]:8.2f}' for m in mults))
    print('\nmain $: what the main run billed on these items at the provider\'s default (no reasoning tokens for these '
          'models, D-168). xN $: the same items with N times the visible output, at the pricing document\'s rates. '
          'The judged stages come on top: E5 and the router on the new traces, priced by their own --dry-run on the '
          'variant once it exists.\n--calibrate 20 runs 20 items per model first (bills about '
          f'{20 / max(len(its), 1) * tot[6]:.2f} at x6) and replaces the multipliers with measured lengths. Nothing was called.')
    return 0


def tool_dry_run(cfg, its, specs, variant) -> int:
    """The tool arm's estimate: the main run's bill on the same items, multiplied for the turns the loop adds (each
    turn resends the conversation, so the prompt side grows with every call and the output is spread over turns).
    The multipliers bracket an assumption, not a measurement: x1.5 for a model that rarely calls the tool, x2.5
    for two or three calls per item, x4 for a model that works through the problem in the tool; --calibrate 20
    measures the real turns and lengths first. Whether each model's endpoints support tools is read from
    OpenRouter's public catalogue (no key, no call)."""
    mults = (1.5, 2.5, 4)
    print(f'variant {variant}: {len(its)} items, {len({i["template_id"] for i in its})} templates; the python tool in the '
          f'request, tool_choice auto, up to {TOOL_MAX_CALLS} calls per item, {TOOL_TIMEOUT} s per script; traces go to traces/{variant}/')
    try:
        cat = catalogue()
    except Exception as exc:                                       # noqa: BLE001
        cat = {}
        print(f'  (OpenRouter catalogue unreachable, {type(exc).__name__}: tool support not checked)')
    print(f'{"model":22s} {"to run":>7s} {"main $":>8s} {"main out tok":>12s} ' + ' '.join(f'{"x" + str(m) + " $":>8s}' for m in mults) + '  tools in supported_parameters')
    tot = Counter()
    for s in specs:
        todo = [it for it in its if it['item_id'] not in existing(s['key'], variant)]
        main = existing(s['key'])
        rows = [main[it['item_id']] for it in todo if it['item_id'] in main]
        bill = sum(r.get('billed_usd') or 0.0 for r in rows)
        ct = sum(r.get('completion_tokens') or 0 for r in rows)
        sup = (cat.get(s['model']) or {}).get('supported_parameters')
        support = 'unknown' if sup is None else ('yes' if 'tools' in sup else 'NO')
        for m in mults:
            tot[m] += bill * m
        tot['main'] += bill
        print(f'{s["key"]:22s} {len(todo):7d} {bill:8.3f} {ct // max(len(rows), 1):12d} ' + ' '.join(f'{bill * m:8.2f}' for m in mults) + f'  {support}')
    print(f'{"TOTAL":22s} {"":7s} {tot["main"]:8.3f} {"":12s} ' + ' '.join(f'{tot[m]:8.2f}' for m in mults))
    print('\nmain $: what the main run billed on these items. xN $: the same bill multiplied for the turns the tool loop adds; '
          'an assumption until --calibrate measures it. E5 and the router on the new traces come on top, priced by their own '
          f'--dry-run on the variant once it exists.\n--calibrate 20 runs 20 items per model first (about ${20 / max(len(its), 1) * tot[2.5]:.2f} '
          'at x2.5) and replaces the multipliers with measured turns and lengths. Nothing was called.')
    return 0


def matched_dry_run(cfg, its, specs, variant) -> int:
    """A matched-configuration variant's estimate, from bills measured with the parameter on, never the main run's.
    reasoning-medium-full: each model's mean bill per row on the 450-item arm reasoning-medium (the rows
    decoding_table_reasoning-medium.json sums), or on its own rows in this variant (its calibration) for a model the
    arm did not run. The cascade: each model's reasoning-medium-full rows on the same items once 20 of them have run (a
    paraphrase is about as long as its original), else its arm rows, else its reasoning-medium-full rows on all items. A
    model with none of these has no estimate: --calibrate 20 first. Beside
    the bill, the same mean tokens at the endpoint routing picks today, from OpenRouter's public list."""
    sources = ([('reasoning-medium', 'the 450-item arm'), (FULL_REASONING, 'own rows')] if variant == FULL_REASONING
               else [(FULL_REASONING, f'{FULL_REASONING}, same items'), ('reasoning-medium', 'the 450-item arm'),
                     (FULL_REASONING, f'{FULL_REASONING}, all items')])       # same items once 20 of them have run
    want = {it['item_id'] for it in its}
    print(f'variant {variant}: {len(its)} items, {len({i["template_id"] for i in its})} templates; reasoning effort '
          f'{MATCHED_EFFORT!r} in the request; traces go to traces/{variant}/')
    print(f'{"model":22s} {"to run":>6s} {"priced from":30s} {"rows":>5s} {"in/row":>7s} {"out/row":>8s} {"$/row":>8s} '
          f'{"estimate $":>10s} {"routed today $":>14s}  endpoint today')
    total = routed_total = 0.0
    for s in specs:
        todo = [it for it in its if it['item_id'] not in existing(s['key'], variant)]
        rows, label = [], 'none: --calibrate 20 first'
        for src, what in sources:
            got = existing(s['key'], src)
            if what.endswith('same items'):
                got = {i: r for i, r in got.items() if i in want}
                if len(got) < 20:
                    continue
            if got:
                rows, label = list(got.values()), what
                break
        if not rows:
            print(f'{s["key"]:22s} {len(todo):6d} {label:30s}')
            continue
        n = len(rows)
        tok = {'in': sum(r.get('prompt_tokens') or 0 for r in rows) / n, 'out': sum(r.get('completion_tokens') or 0 for r in rows) / n}
        per = sum(r.get('billed_usd') or 0.0 for r in rows) / n
        try:
            ep, price, _n_ok = routed_endpoint(s, cfg, tok)
        except Exception as exc:                                   # noqa: BLE001
            ep, price = None, None
            print(f'  ({s["key"]}: endpoint list unreachable, {type(exc).__name__})')
        routed = len(todo) * per_item_cost(price, tok) if price else None
        total += per * len(todo)
        routed_total += routed or 0.0
        where = (f"{ep.get('provider_name') or ep.get('name')}{' ' + ep['quantization'] if s['weights'] == 'open' and ep.get('quantization') else ''} "
                 f"${price['in']:.3f}/${price['out']:.3f} per M" if ep else 'NO ELIGIBLE ENDPOINT')
        print(f'{s["key"]:22s} {len(todo):6d} {label:30s} {n:5d} {tok["in"]:7.0f} {tok["out"]:8.0f} {per:8.5f} '
              f'{per * len(todo):10.2f} {"" if routed is None else f"{routed:14.2f}":>14s}  {where}')
    print(f'{"TOTAL":22s} {"":6s} {"":30s} {"":5s} {"":7s} {"":8s} {"":8s} {total:10.2f} {routed_total:14.2f}')
    print('\n$/row: the mean bill per row of the rows named, with the parameter on; estimate: $/row times the rows to run. '
          'routed today: the same mean tokens at the endpoint routing picks today (the bill is read from each response). '
          'E5 and the router on the new traces come on top, priced by their own --dry-run once the traces exist. '
          'Nothing was called.')
    return 0


def review(cfg, its, specs, variant) -> int:
    """FREE. What a finished run must show before it is scored, per model on this variant's items: the items missing,
    the finish reasons, the empty rows and those that ended at the ceiling, the rows whose endpoint reported reasoning
    tokens and those that returned reasoning text, the request sets, whether every row carries the prompt hash the
    main run's rows carry, and the providers that served. Exits 1 when an item is missing or a prompt differs."""
    main_sha = main_prompt_hashes()
    want = {it['item_id'] for it in its}
    print(f'variant {variant}: {len(its)} items; the main run\'s rows carry prompt '
          + ', '.join(f'{h[:16]} ({n} rows)' for h, n in main_sha.most_common()))
    print(f'{"model":22s} {"rows":>5s} {"missing":>7s} {"empty":>5s} {"at ceiling":>10s} {"stop/length/other":>17s} '
          f'{"reasoning tok":>13s} {"reasoning txt":>13s} {"requests":>8s} {"prompt = main":>13s} {"billed $":>9s}  providers')
    bad = 0
    for s in specs:
        rows = [r for i, r in existing(s['key'], variant).items() if i in want]
        n = len(rows)
        fin = Counter(r.get('finish_reason') for r in rows)
        finish = f"{fin['stop']}/{fin['length']}/{n - fin['stop'] - fin['length']}"
        empty = sum(r['status'] == 'empty' for r in rows)
        at_ceiling = sum(r['status'] == 'empty' and r.get('finish_reason') == 'length' for r in rows)
        share = (lambda k: f'{k / n:.3f}') if n else (lambda k: '-')
        rtok = share(sum((r.get('reasoning_tokens') or 0) > 0 for r in rows))
        rtxt = share(sum(bool((r.get('reasoning') or '').strip()) for r in rows))
        reqs = len({json.dumps(r.get('request'), sort_keys=True) for r in rows})
        same = n > 0 and {r.get('prompt_sha256') for r in rows} == set(main_sha) == {PROMPT_SHA}
        provs = Counter(r.get('provider') for r in rows)
        bad += (n < len(its)) or not same
        print(f'{s["key"]:22s} {n:5d} {len(its) - n:7d} {empty:5d} {at_ceiling:10d} {finish:>17s} {rtok:>13s} {rtxt:>13s} '
              f'{reqs:8d} {"yes" if same else "NO":>13s} {sum(r.get("billed_usd") or 0.0 for r in rows):9.3f}  '
              + '; '.join(f'{p} ({c})' for p, c in provs.most_common()))
    print('\nrows: final rows (answered or empty) on this variant\'s items. at ceiling: empty rows that ended at the output '
          'ceiling. reasoning tok / txt: the share of rows whose endpoint reported reasoning tokens above zero / returned '
          'reasoning text. requests: distinct request parameter sets (1 expected). prompt = main: every row carries the '
          'prompt hash the main run\'s rows carry.')
    for s in specs:                                        # what the variant billed beyond its final rows
        moved = TRACES / '_unfinished' / variant / f'{s["key"]}.jsonl'
        rows = [json.loads(ln) for ln in moved.read_text(encoding='utf-8').splitlines()] if moved.exists() else []
        if rows:
            print(f'{s["key"]}: {len(rows)} rows without a finish reason moved to traces/_unfinished/ and asked again, '
                  f'billed ${sum(r.get("billed_usd") or 0.0 for r in rows):.3f} besides the final rows')
    return 1 if bad else 0


def matched_config() -> int:
    """FREE. results/matched_config.json, which the analysis's matched family reads: for every roster model the store of
    its provider-default run (main), the store of its reasoning-on run if it has one, the setting ("effort=medium", "none
    offered" or "reasons by default"), a re-run model's cascade stores, and the evidence, counted from the trace rows
    with what the model's card says; the Qwen Thinking sibling has its own entry (sibling_of; run false while it is
    priced and not run). Refuses to write while a store it names is incomplete or the rows contradict the setting."""
    cfg = config()
    pool = items()
    sets = {v: {it['item_id'] for it in variant_items(v, pool)} for v in ('main',) + MATCHED_VARIANTS}

    def count(variant, key):
        rows = [r for i, r in existing(key, variant).items() if i in sets[variant]]
        return len(rows), sum((r.get('reasoning_tokens') or 0) > 0 for r in rows)

    models, bad = [], []
    for s in cfg['models']:
        key, sibling = s['key'], s.get('added_for') == 'matched settings'
        if s.get('run') is False and not sibling:
            continue                                       # never run, set aside, or an anchor outside the roster
        n, k = count('main', key)
        ev = [CARD_NOTES[key]] if key in CARD_NOTES else []
        ev.append(f"reasoning tokens on {k} of {n} main-run rows at the provider's default" if n else
                  s.get('run_reason') or 'not run')
        e = {'model': key, 'default_store': 'main' if n == len(sets['main']) else None}
        if sibling:
            e.update(sibling_of=s['sibling_of'], run=n > 0)
        if key in MATCHED_MODELS:
            n2, k2 = count(FULL_REASONING, key)
            ev.append(f'on {k2} of {n2} {FULL_REASONING} rows with the parameter at {MATCHED_EFFORT}')
            reps = list(CASCADE_REPEATS) if key in REPEAT_MODELS else []
            e.update(reasoning_store=FULL_REASONING, reasoning_setting=f'effort={MATCHED_EFFORT}',
                     cascade={'paraphrase': CASCADE_PARAPHRASE, 'repeats': reps})
            for v in (FULL_REASONING, CASCADE_PARAPHRASE, *reps):
                if count(v, key)[0] != len(sets[v]):
                    bad.append(f'{key}: {v} holds {count(v, key)[0]} of {len(sets[v])} items')
            if k or k2 < 0.95 * n2:
                bad.append(f'{key}: reasoning tokens on {k} main rows and on {k2} of {n2} {FULL_REASONING} rows')
        elif key in NO_SETTING:
            e.update(reasoning_store=None, reasoning_setting='none offered')
            if k:
                bad.append(f'{key}: no setting is offered, yet {k} main rows carry reasoning tokens')
        else:
            e.update(reasoning_store=None, reasoning_setting='reasons by default')
            if n and k < 0.95 * n:
                bad.append(f'{key}: reasons by default, yet only {k} of {n} main rows carry reasoning tokens')
        if n != len(sets['main']) and not (sibling and n == 0):
            bad.append(f'{key}: the main run holds {n} of {len(sets["main"])} items')
        e['evidence'] = '; '.join(ev)
        models.append(e)
    if bad:
        print('matched_config.json NOT written:\n  ' + '\n  '.join(bad))
        return 1
    MATCHED_CONFIG.parent.mkdir(exist_ok=True)
    MATCHED_CONFIG.write_text(json.dumps({
        'generated_by': 'python -m full_run_28092026.run_traces --matched-config',
        'written_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'items': len(sets['main']), 'effort': MATCHED_EFFORT, 'prompt_sha256': PROMPT_SHA,
        'models': models}, indent=1) + '\n', encoding='utf-8', newline='\n')
    for e in models:
        print(f"{e['model']:30s} {e['reasoning_setting']:20s} default {e['default_store'] or '-':5s} "
              f"reasoning {e['reasoning_store'] or '-'}")
    print(f'wrote {MATCHED_CONFIG.relative_to(REPO).as_posix()}')
    return 0


def requeue_unfinished(specs, variant) -> int:
    """FREE. A matched-configuration variant's final rows that ended without a finish reason (answered, the text cut
    off by the provider; call() now asks such a completion again) leave the trace file for
    traces/_unfinished/<variant>/<model>.jsonl, kept there, so the next run of the variant calls those items again.
    Only when no process is writing the variant: it refuses a trace file changed in the last two minutes."""
    for s in specs:
        path = trace_path(s['key'], variant)
        if not path.exists():
            continue
        if time.time() - path.stat().st_mtime < 120:
            raise SystemExit(f'{variant}/{path.name} changed in the last two minutes: is a run still writing it?')
        lines = path.read_text(encoding='utf-8').splitlines()
        final = {}
        for i, ln in enumerate(lines):
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if r.get('status') in DONE:
                final[r['item_id']] = (i, r)
        drop = {i: r['item_id'] for i, r in final.values() if r['status'] == 'answered' and r.get('finish_reason') is None}
        if not drop:
            print(f'{variant}/{s["key"]}: no final row without a finish reason')
            continue
        dest = TRACES / '_unfinished' / variant / f'{s["key"]}.jsonl'
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, 'a', encoding='utf-8', newline='\n') as fh:
            for i in sorted(drop):
                fh.write(lines[i] + '\n')
        tmp = path.with_name(path.name + '.tmp')
        with open(tmp, 'w', encoding='utf-8', newline='\n') as fh:
            fh.writelines(ln + '\n' for i, ln in enumerate(lines) if i not in drop)
        os.replace(tmp, path)
        print(f'{variant}/{s["key"]}: {len(drop)} final rows without a finish reason moved to '
              f'{dest.relative_to(HERE).as_posix()}: {", ".join(sorted(drop.values()))}')
    return 0


def trace_review(variant: str) -> int:
    """FREE. trace_review.py on a matched-configuration variant, given the variant's item set: its own variant_manifest()
    reads any name that begins `reasoning-` as the 450-item subsample and an unknown name as the whole pool, so it would
    count most of these variants' rows as foreign or missing. Everything else is trace_review.py's own, which writes
    TRACE_REVIEW_<variant>.md and trace_review_<variant>.json."""
    from full_run_28092026 import trace_review as tr
    man = {it['item_id']: {'sha256': it['sha256'], **({'original_sha256': it['original_sha256']} if 'original_sha256' in it else {})}
           for it in variant_items(variant, items())}
    tr.variant_manifest = lambda _variant, _man: man
    sys.argv = ['trace_review', '--variant', variant]
    return tr.main()


def tool_selftest() -> int:
    """FREE. The sandbox (a numeric script, a refused one, a timeout, truncation, the environment without our key)
    and the loop against a fake client that returns a tool call and then an answer."""
    from types import SimpleNamespace as NS
    fails = []

    def check(name, cond):
        print(f'  {"ok " if cond else "FAIL"} {name}')
        if not cond:
            fails.append(name)

    r = run_python('import numpy as np\nimport scipy, sympy\nprint(np.sqrt(16.0), 2 * 3)')
    check('numpy, scipy and sympy import in the isolated interpreter; output captured', r['ok'] and '4.0 6' in r['output'])
    r = run_python('import os\nprint(os.listdir("."))')
    check('os is refused by the filter', r['refused'] and 'ToolError' in r['output'])
    r = run_python('f = open("x.txt", "w")')
    check('open() is refused by the filter', r['refused'])
    r = run_python('print(1/0)')
    check('a raising script returns ok False with the traceback', (not r['ok']) and 'ZeroDivisionError' in r['output'])
    r = run_python('while True:\n    pass', timeout=2)
    check('a runaway script is stopped at the limit', (not r['ok']) and 'TimeoutError' in r['output'])
    r = run_python('print("x" * 10000)')
    check('long output is truncated', r['truncated'] and len(r['output']) < 4200)
    env = sandbox_env('C:/tmp')
    check('no key, token or secret in the sandbox environment',
          not any(k for k in env if any(w in k.upper() for w in ('KEY', 'TOKEN', 'SECRET'))) and 'OPENROUTER_API_KEY' not in env)

    class Fake:
        """Two turns: a tool call computing 2*3, then the answer; records what it was sent."""
        def __init__(self, repair=False):
            self.sent, self.repair = [], repair
            self.chat = NS(completions=NS(create=self.create))

        def create(self, model, messages, **kw):
            self.sent.append((list(messages), kw))
            usage = NS(prompt_tokens=100, completion_tokens=20, completion_tokens_details=NS(reasoning_tokens=5), model_extra={'cost': 0.001})
            if len(self.sent) == 1:
                tc = NS(id='call_1', type='function', function=NS(name='python', arguments='{"code": "print(2*3)"}'))
                msg = NS(content='Let me compute.', tool_calls=[tc], model_extra={}, reasoning=None)
            elif len(self.sent) == 2 and self.repair:        # a literal newline inside the JSON string: not strict JSON
                tc = NS(id='call_2', type='function', function=NS(name='python', arguments='{"code": "x = 7\nprint(x)"}'))
                msg = NS(content=None, tool_calls=[tc], model_extra={}, reasoning=None)
            else:
                msg = NS(content='The product is 6.\n\nFinal answer: 6', tool_calls=None, model_extra={}, reasoning=None)
            return NS(choices=[NS(message=msg, finish_reason='stop')], usage=usage, model='fake/model', model_extra={'provider': 'Fake'})

    fake = Fake()
    cfg = {'max_tokens': 1000, 'provider_closed_weight': {'sort': 'price'}, 'provider_open_weight': {'sort': 'price'}}
    row = call_with_tool(fake, {'key': 'fake', 'model': 'fake/model', 'weights': 'closed'}, cfg, 'What is 2 times 3?')
    check('the loop ends answered after one tool call and two turns', row['status'] == 'answered' and row['tool_calls'] == 1 and row['turns'] == 2)
    check('the transcript carries the text, the script, its output and the answer',
          'Let me compute.' in row['text'] and '```python\nprint(2*3)\n```' in row['text'] and 'Output:\n```\n6\n```' in row['text'] and row['text'].endswith('Final answer: 6'))
    check('final_text is the last message alone', row['final_text'] == 'The product is 6.\n\nFinal answer: 6')
    check('tokens and the bill are summed over the turns', row['prompt_tokens'] == 200 and row['completion_tokens'] == 40 and row['reasoning_tokens'] == 10 and abs(row['billed_usd'] - 0.002) < 1e-9)
    second = fake.sent[1][0]
    check('the second turn was sent the assistant tool call and the tool result',
          len(second) == 3 and second[1]['role'] == 'assistant' and second[1]['tool_calls'][0]['function']['name'] == 'python'
          and second[2] == {'role': 'tool', 'tool_call_id': 'call_1', 'content': '6'})
    check('the request carried the python tool with tool_choice auto', fake.sent[0][1].get('tools') == [PYTHON_TOOL] and fake.sent[0][1].get('tool_choice') == 'auto')
    check('the tool arm runs the subsample\'s originals', TOOL in VARIANTS and 'tool' == TOOL)
    fake2 = Fake(repair=True)
    row2 = call_with_tool(fake2, {'key': 'fake', 'model': 'fake/model', 'weights': 'closed'}, cfg, 'What is 2 times 3?')
    sent_back = fake2.sent[2][0][3]['tool_calls'][0]['function']['arguments']
    check('arguments with a literal newline are read, run, and sent back as valid JSON',
          row2['tool_calls'] == 2 and row2['tool_turns'][1]['args_repaired'] and row2['tool_turns'][1]['output'] == '7'
          and json.loads(sent_back) == {'code': 'x = 7\nprint(x)'} and row2['status'] == 'answered')
    print(f'selftest: {"all pass" if not fails else str(len(fails)) + " FAILED: " + "; ".join(fails)}')
    return 1 if fails else 0


def calibration_sample(its, n):
    """n items per model, spread over the templates in a fixed order: every k-th item."""
    step = max(1, len(its) // n)
    return its[::step][:n]


def status(cfg, its, specs, variant='main', only: bool = False) -> int:
    print(f'{"model":22s} {"answered":>9s} {"empty":>6s} {"missing":>8s} {"billed $":>9s} '
          f'{"out tok/item":>13s}  served')
    grand = 0.0
    for s in specs:
        rows = list(existing(s['key'], variant).values())
        if only:                                                       # --only-items: the rows of the items listed
            want = {it['item_id'] for it in its}
            rows = [r for r in rows if r['item_id'] in want]
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
    ap.add_argument('--selftest', action='store_true', help='FREE: the tool arm offline (the sandbox and the loop with a fake client)')
    ap.add_argument('--ignore-provider', action='append', default=[], metavar='NAME',
                    help='route past this OpenRouter provider in this invocation (repeatable); the row records it in its request. '
                         'Used for the tool arm when a provider drops tool calls (D-184)')
    ap.add_argument('--only-items', metavar='PATH', help='a text file of item ids, one per line: every mode works on the '
                    'variant\'s items among them only (the resume of a subset whose rows were archived)')
    ap.add_argument('--review', action='store_true', help='FREE: the checks a finished run must pass before scoring '
                    '(missing items, finish reasons, empty rows at the ceiling, reasoning tokens and text, providers, '
                    'the prompt hash against the main run\'s)')
    ap.add_argument('--matched-config', action='store_true', help='FREE: write results/matched_config.json, the matched '
                    'configuration per roster model, once every store it names is complete')
    ap.add_argument('--trace-review', action='store_true', help='FREE: trace_review.py on a matched-configuration variant, '
                    'given that variant\'s item set')
    ap.add_argument('--requeue-unfinished', action='store_true', help='FREE: move a matched-configuration variant\'s final '
                    'rows that ended without a finish reason to traces/_unfinished/, so the next run asks those items again')
    a = ap.parse_args()
    PROVIDER_IGNORE[:] = a.ignore_provider
    if a.selftest:
        return tool_selftest()
    if a.matched_config:
        return matched_config()
    if a.trace_review:
        if a.variant not in MATCHED_VARIANTS:
            raise SystemExit('--trace-review is for the matched-configuration variants; run trace_review.py for the others')
        return trace_review(a.variant)

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
    if a.variant != 'main' and a.variant not in REASONING_VARIANTS and a.variant not in ('flagship', TOOL, FULL_REASONING) + FLAGSHIP_REASONING and a.calibrate:
        raise SystemExit('--calibrate is for the main run, the reasoning variants, the flagship anchor and the tool arm; another variant is estimated from its bills')
    if a.variant in REASONING_VARIANTS and not a.model:
        specs = [s for s in cfg['models'] if s['key'] in REASONING_MODELS]
    if a.variant == 'flagship':
        specs = [s for s in cfg['models'] if s.get('anchor') and (not a.model or s['key'] in a.model)]
        if not specs:
            raise SystemExit('no anchor model in models.json (anchor: true)')
    if a.variant in FLAGSHIP_REASONING:            # the closed anchor(s): the open one reasons at its default already
        specs = [s for s in cfg['models'] if s.get('anchor') and (s['key'] in a.model if a.model else s['weights'] == 'closed')]
        if not specs:
            raise SystemExit('no closed anchor model in models.json')
    if a.variant in ('openbook', 'openbook2') and not a.model:
        specs = [s for s in cfg['models'] if s['key'] in OPENBOOK_MODELS]
    if a.variant == TOOL and not a.model:
        specs = [s for s in cfg['models'] if s['key'] in TOOL_MODELS]
    if a.variant in subsamples.REPEAT_VARIANTS and not a.model:
        raise SystemExit('a repeat runs one chosen model: name it with --model')
    if a.variant == FULL_REASONING and not a.model:
        specs = [s for s in cfg['models'] if s['key'] in MATCHED_MODELS]
    if a.variant in (CASCADE_PARAPHRASE,) + CASCADE_REPEATS:          # the cascade runs the matched models only
        allowed = [k for k in MATCHED_MODELS if a.variant == CASCADE_PARAPHRASE or k in REPEAT_MODELS]
        outside = sorted(set(a.model or ()) - set(allowed))
        if outside:
            raise SystemExit(f'{a.variant} runs the matched models{" in the repeat set" if a.variant != CASCADE_PARAPHRASE else ""} '
                             f'only ({", ".join(allowed) or "none"}), not {", ".join(outside)}')
        specs = [s for s in cfg['models'] if s['key'] in (a.model or allowed)]
        if not specs:
            raise SystemExit(f'{a.variant}: no matched model to run')
    pool = items()
    its = variant_items(a.variant, pool)
    if a.only_items:
        its = only_items(its, pool, a.only_items)
        if not its:
            print('nothing to do: no listed id is an item of this variant')
            return 0
    if a.dry_run:
        if a.variant == TOOL:
            return tool_dry_run(cfg, its, specs, a.variant)
        if a.variant in MATCHED_VARIANTS:
            return matched_dry_run(cfg, its, specs, a.variant)
        if a.variant in ('main', 'flagship'):
            return dry_run(cfg, its, specs, a.variant)
        return reasoning_dry_run(cfg, its, specs, a.variant) if a.variant in REASONING_VARIANTS + FLAGSHIP_REASONING \
            else variant_dry_run(cfg, its, specs, a.variant)
    if a.status:
        return status(cfg, its, specs, a.variant, only=bool(a.only_items))
    if a.review:
        return review(cfg, its, specs, a.variant)
    if a.requeue_unfinished:
        if a.variant not in MATCHED_VARIANTS:
            raise SystemExit('--requeue-unfinished is for the matched-configuration variants')
        return requeue_unfinished(specs, a.variant)
    if not a.yes:
        raise SystemExit('this mode bills: re-run with --yes once the spend is approved '
                         '(see --dry-run for the estimate)')
    if not PROMPT_SHA.startswith(PILOT_PROMPT_PREFIX):
        raise SystemExit('the deployed prompt differs from the one the pilot traces used - refusing')
    if a.check:
        bad = 0
        if a.variant == 'main' or a.variant in MATCHED_VARIANTS:     # the matched configuration asks the main run's prompt
            hashes = main_prompt_hashes()
            same = set(hashes) == {PROMPT_SHA}
            bad += not same
            print(f'  prompt sha256 {PROMPT_SHA[:16]}; the main run\'s rows carry '
                  + ', '.join(f'{h[:16]} ({n} rows)' for h, n in hashes.most_common()) + (': the same' if same else ': DIFFERENT'))
        v = a.variant if a.variant in MATCHED_VARIANTS else 'main'   # a matched variant's check carries its parameter
        for s in specs:
            res = call(client(cfg), s, cfg, 'Reply with the single word OK.', attempts=2, variant=v)
            ok = res['status'] == 'answered'
            bad += not ok
            print(f'  {s["key"]:22s} {res["status"]:16s} served={res.get("served_model")} '
                  f'provider={res.get("provider")} billed=${res.get("billed_usd") or 0:.5f}'
                  + (f' reasoning_tokens={res.get("reasoning_tokens")}' if v != 'main' else ''))
        return 1 if bad else 0
    if a.variant in MATCHED_VARIANTS and set(main_prompt_hashes()) != {PROMPT_SHA}:
        raise SystemExit('the prompt differs from the one on the main run\'s rows - refusing')
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
