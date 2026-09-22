"""What generating answers costs for the full benchmark - inference only, no evaluation.

    python evaluator_pilot_17092026/analysis/inference_cost.py

Every model in the paper's suite (Table 13, 27 models) answers every item of the full
benchmark: 150 templates x 15 instances = 2,250 items per model.

PRICES are read live from OpenRouter. The ':batch' variant (half price, asynchronous)
is shown where OpenRouter lists one; a benchmark run needs no interactivity.

TOKENS PER ITEM come from the pilot's own traces (60 items across all five branches
and three levels, the benchmark's prompt), measured per model where the pilot ran that
model, and otherwise borrowed from the closest measured model:
  measured    gpt-5, claude-opus-4.7, deepseek-r1(-0528), llama-3.1-70b, gemma, qwen3.8
  reasoning   models that think before answering get a measured reasoning model's
              profile (GPT-5 family -> gpt-5; Qwen3 -> qwen3.8; DeepSeek thinking -> r1)
  standard    non-reasoning models get the mean of the pilot's non-reasoning models
GEMINI is the one real unknown: the pilot recorded no thinking tokens for it (the API
did not report them), so its billed output is estimated as its visible output plus
GPT-5's measured reasoning, and shown with a low/high range.

RETRIES are added as the pilot measured them (re-runs after truncation or failure),
with a floor of 15% for reasoning models and 5% for standard ones.

Models OpenRouter does not serve are marked LOCAL: they run on HiPerGator at no API
cost (7B-14B models, minutes of GPU time each).
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import glob
import json
import statistics as st
import urllib.request

ITEMS = 2250

# name, group, OpenRouter id (None = not served), token profile, notes
SUITE = [
    ('GPT-5', 'frontier', 'openai/gpt-5', 'gpt-5', ''),
    ('GPT-5 Mini', 'frontier', 'openai/gpt-5-mini', 'gpt-5', 'reasoning, profile of GPT-5'),
    ('GPT-4.1', 'frontier', 'openai/gpt-4.1', 'standard', ''),
    ('GPT-4.1 Mini', 'frontier', 'openai/gpt-4.1-mini', 'standard', ''),
    ('Claude Opus 4.7', 'frontier', 'anthropic/claude-opus-4.7', 'claude-opus-4.7', 'no extended thinking, as in the pilot'),
    ('Claude Sonnet 4.5', 'frontier', 'anthropic/claude-sonnet-4.5', 'claude-opus-4.7', ''),
    ('Claude Sonnet 4', 'frontier', 'anthropic/claude-sonnet-4', 'claude-opus-4.7', ''),
    ('Claude Sonnet 3.7', 'frontier', None, 'claude-opus-4.7', 'RETIRED on OpenRouter'),
    ('Gemini 3.1 Pro', 'frontier', 'google/gemini-3.1-pro-preview', 'gemini', ''),
    ('Gemini 3 Pro', 'frontier', None, 'gemini', 'RETIRED (pilot finding E0-F4)'),
    ('Gemini 2.5 Pro', 'frontier', 'google/gemini-2.5-pro', 'gemini', ''),
    ('Gemini 2.5 Flash', 'frontier', 'google/gemini-2.5-flash', 'gemini', 'thinks by default'),
    ('DeepSeek V4 Pro', 'frontier', 'deepseek/deepseek-v4-pro', 'deepseek-r1', 'assumed thinking mode'),
    ('DeepSeek V3', 'frontier', 'deepseek/deepseek-chat-v3-0324', 'standard', ''),
    ('DeepSeek R1', 'frontier', 'deepseek/deepseek-r1-0528', 'deepseek-r1', 'R1-0528, as in the pilot'),
    ('Llama 3.1 70B', 'open', 'meta-llama/llama-3.1-70b-instruct', 'llama-3.1-70b', ''),
    ('Llama 3.1 8B', 'open', 'meta-llama/llama-3.1-8b-instruct', 'standard', ''),
    ('Qwen 2.5 72B', 'open', 'qwen/qwen-2.5-72b-instruct', 'standard', ''),
    ('Qwen 2.5 14B', 'open', None, 'standard', 'LOCAL'),
    ('Qwen 2.5 7B', 'open', 'qwen/qwen-2.5-7b-instruct', 'standard', ''),
    ('Qwen 3 8B', 'open', 'qwen/qwen3-8b', 'qwen3.8-27b', 'thinking mode, profile of Qwen3.8'),
    ('Gemma 3 27B', 'open', 'google/gemma-3-27b-it', 'gemma-4-31b', ''),
    ('Gemma 2 9B', 'open', None, 'gemma-4-31b', 'LOCAL'),
    ('Qwen2.5-Math 7B', 'math', None, 'standard', 'LOCAL'),
    ('Mathstral 7B', 'math', None, 'standard', 'LOCAL'),
    ('WizardMath 7B', 'math', None, 'standard', 'LOCAL'),
    ('MetaMath 7B', 'math', None, 'standard', 'LOCAL'),
]


def measured():
    """Per pilot model: mean prompt tokens, mean completion tokens, retry rate."""
    out = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        rows = [json.loads(l) for l in open(f, encoding='utf-8')]
        ok = [r for r in rows if r['ok']]
        items = len({r['item_id'] for r in rows})
        out[_os.path.basename(f)[:-6]] = {
            'in': st.mean(r.get('prompt_tokens') or 0 for r in ok),
            'out': st.mean(r.get('completion_tokens') or 0 for r in ok),
            'retry': len(rows) / items - 1,
        }
    std = [out[k] for k in ('llama-3.1-70b', 'claude-opus-4.7', 'gemma-4-31b')]
    out['standard'] = {'in': st.mean(x['in'] for x in std), 'out': st.mean(x['out'] for x in std), 'retry': 0.05}
    g = out.pop('gemini-3.1-pro')
    reasoning_gpt5 = 4538                                      # measured mean reasoning tokens, gpt-5
    out['gemini'] = {'in': g['in'], 'out': g['out'] + reasoning_gpt5, 'retry': g['retry'],
                     'low': g['out'], 'high': g['out'] + 7265}  # 7265 = deepseek-r1 mean reasoning
    return out


def prices():
    from dotenv import load_dotenv
    load_dotenv(_os.path.join(_REPO, '.env'))
    req = urllib.request.Request('https://openrouter.ai/api/v1/models',
                                 headers={'Authorization': 'Bearer ' + _os.getenv('OPENROUTER_API_KEY'),
                                          'User-Agent': 'engtrace-inference-costing'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return {m['id']: (float(m['pricing']['prompt']) * 1e6, float(m['pricing']['completion']) * 1e6)
                for m in json.load(fh)['data']}


def main():
    tok, pr = measured(), prices()
    reasoning_profiles = {'gpt-5', 'deepseek-r1', 'qwen3.8-27b', 'gemini'}
    print('Inference only, %d items per model. Prices live from OpenRouter ($ per million tokens).\n' % ITEMS)
    print('%-18s %-9s %7s %7s %6s %9s %9s  %s' % ('model', 'group', 'in/item', 'out/item', 'retry',
                                                  'standard', 'batch', 'note'))
    tot = {g: [0.0, 0.0] for g in ('frontier', 'open', 'math')}
    gem_range = [0.0, 0.0]
    for name, group, oid, prof, note in SUITE:
        t = tok[prof]
        retry = max(t['retry'], 0.15 if prof in reasoning_profiles else 0.05)
        if oid is None or oid not in pr:
            print('%-18s %-9s %7.0f %7.0f %5.0f%% %9s %9s  %s' % (name, group, t['in'], t['out'], 100 * retry,
                                                                '-', '-', note or 'not on OpenRouter'))
            continue
        pin, pout = pr[oid]
        n = ITEMS * (1 + retry)
        cost = n * (t['in'] * pin + t['out'] * pout) / 1e6
        b = pr.get(oid + ':batch')
        bcost = n * (t['in'] * b[0] + t['out'] * b[1]) / 1e6 if b else cost
        tot[group][0] += cost
        tot[group][1] += bcost
        if prof == 'gemini':
            gem_range[0] += n * (t['in'] * pin + t['low'] * pout) / 1e6
            gem_range[1] += n * (t['in'] * pin + t['high'] * pout) / 1e6
        print('%-18s %-9s %7.0f %7.0f %5.0f%% %9.2f %9s  %s' % (name, group, t['in'], t['out'], 100 * retry, cost,
                                                              '%.2f' % bcost if b else '-', note))
    print()
    for g, (c, b) in tot.items():
        print('  %-10s standard $%8.2f   with batch where offered $%8.2f' % (g, c, b))
    c = sum(v[0] for v in tot.values())
    b = sum(v[1] for v in tot.values())
    print('  %-10s standard $%8.2f   with batch where offered $%8.2f' % ('TOTAL', c, b))
    print('\n  Gemini models alone, standard price: $%.2f if they barely think, $%.2f if they think like R1'
          % tuple(gem_range))


if __name__ == '__main__':
    main()
