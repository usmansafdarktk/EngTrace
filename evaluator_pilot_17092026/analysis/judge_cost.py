"""What the LLM judges cost on the full benchmark run - evaluation only, no inference.

    python evaluator_pilot_17092026/analysis/judge_cost.py [--models 12]

Scales each paid evaluator from what it actually cost per trace on the pilot (the cost
recorded on every judge call, current config only). Only E0, E0-3J, E1 and E5 call
LLM judges; E2 (process reward models on our GPUs), E3 and E4 cost nothing.

How often the judges are called depends on the model being judged: E0/E1 send a trace
to the judges when its steps do not match the reference, E5 only for milestones its
deterministic matching misses. Weak models trigger far more calls. So two rates are
shown per evaluator: the pilot's four frontier models, and Llama 3.1 70B, the weak
model. The next roster is mostly smaller open models, so the truth sits between the two.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import glob
import json
import statistics as st

ITEMS = 2250
EVALS = [('E0', 'e0', 'published framework, 2 judges (GPT-5, Claude Opus 4.5)'),
         ('E0-3J', 'e0_3j', 'E0 with its third judge (Gemini 3.1 Pro) connected'),
         ('E1', 'e1', 'non-suite panel: Grok 4.6, MiniMax M3, MiMo-V2.5-Pro'),
         ('E5', 'e5', 'E3 first, MiMo-V2.5-Pro only on missed milestones')]
FRONTIER = ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1')
WEAK = ('llama-3.1-70b',)


def current(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error') and r.get('cohort', 'gold') == 'gold']
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return [r for r in rows if r['config_sha256'] == latest]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', type=int, default=12)
    a = ap.parse_args()
    n = a.models * ITEMS
    print('Judge cost for %d models x %d items = %d traces (pilot rates per trace)\n' % (a.models, ITEMS, n))
    print('%-6s %-52s %10s %10s %11s %11s' % ('eval', 'judges', '$/trace F', '$/trace W', 'full run F', 'full run W'))
    for name, d, desc in EVALS:
        rows = current(d)
        per = {}
        for grp, keys in (('F', FRONTIER), ('W', WEAK)):
            rs = [r for r in rows if r['model_key'] in keys]
            per[grp] = st.mean(sum(c.get('cost_usd', 0) for c in r.get('calls') or []) for r in rs)
        print('%-6s %-52s %10.5f %10.5f %11.2f %11.2f' % (name, desc, per['F'], per['W'], per['F'] * n, per['W'] * n))
    print('\nF = at the pilot frontier models\' rate, W = at Llama 3.1 70B\'s rate. E2, E3, E4: $0.')


if __name__ == '__main__':
    main()
