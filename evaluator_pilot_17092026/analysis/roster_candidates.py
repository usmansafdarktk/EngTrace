"""Which open-weight models could the next benchmark round evaluate?

    python evaluator_pilot_17092026/analysis/roster_candidates.py

Screens OpenRouter's live catalogue for open-weight candidates and prints price,
context and whether the family is blocked. A family is blocked when it already
judges (D-088: a family cannot be both judge and judged) or when it supplies an
evaluator's scoring model (E2's PRMs are Qwen). Costs are projected onto the full
benchmark, 2,250 items, from the pilot's measured token counts.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import glob
import json
import urllib.request

FULL = 2250

# family -> why it cannot be evaluated in the next round, or None
BLOCKED = {
    'xiaomi': 'judges E5 (D-088)',
    'minimax': 'judges E1 (D-088)',
    'x-ai': 'judges E1 (D-088)',
    'openai': 'in the paper\'s Table 13 suite; judges E0',
    'anthropic': 'in the suite; judges E0',
    'google': 'in the suite; judges E0 (Gemma is Google)',
    'qwen': 'E2\'s process reward models are Qwen2.5-Math-PRM',
}

CANDIDATES = [
    'deepseek/deepseek-r1-0528', 'meta-llama/llama-3.1-70b-instruct',
    'google/gemma-4-31b-it', 'qwen/qwen3.8-27b', 'z-ai/glm-5.3',
    'moonshotai/kimi-k3', 'nvidia/nemotron-3-ultra-550b-a55b',
    'bytedance-seed/seed-2-1-turbo', 'mistralai/magistral-medium-2509',
    'minimax/minimax-m3', 'xiaomi/mimo-v2.5-pro',
]


def catalogue():
    from dotenv import load_dotenv
    load_dotenv(_os.path.join(_REPO, '.env'))
    key = _os.getenv('OPENROUTER_API_KEY')
    if not key:
        raise SystemExit('OPENROUTER_API_KEY not found in .env')
    req = urllib.request.Request('https://openrouter.ai/api/v1/models',
                                 headers={'Authorization': 'Bearer ' + key,
                                          'User-Agent': 'engtrace-roster-screen'})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return {m['id']: m for m in json.load(fh)['data']}


def pilot_tokens():
    """Mean prompt and completion tokens per trace, measured over the pilot's traces."""
    pin = pout = n = 0
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        for l in open(f, encoding='utf-8'):
            r = json.loads(l)
            if r['ok']:
                pin += r.get('prompt_tokens') or 0
                pout += r.get('completion_tokens') or 0
                n += 1
    return pin / n, pout / n


def main():
    cat = catalogue()
    tin, tout = pilot_tokens()
    print('pilot mean tokens per trace: %.0f in, %.0f out (completion includes reasoning '
          'where the provider reports it)' % (tin, tout))
    print('projected onto the full benchmark: %d items per model\n' % FULL)
    print('%-34s %10s %10s %9s  %s' % ('model', '$/M in', '$/M out', 'gen $', 'status'))
    for mid in CANDIDATES:
        m = cat.get(mid)
        fam = mid.split('/')[0]
        why = BLOCKED.get(fam) or BLOCKED.get(mid.split('/')[1].split('-')[0])
        if m is None:
            print('%-34s %10s %10s %9s  %s' % (mid, '-', '-', '-', 'NOT ON OPENROUTER'))
            continue
        pin = float(m['pricing']['prompt']) * 1e6
        pout = float(m['pricing']['completion']) * 1e6
        gen = FULL * (tin * pin + tout * pout) / 1e6
        print('%-34s %10.3f %10.3f %9.0f  %s'
              % (mid, pin, pout, gen, ('BLOCKED: ' + why) if why else 'available'))


if __name__ == '__main__':
    main()
