"""What the LLM judges cost on the full benchmark run - evaluation only, no inference.

    python evaluator_pilot_17092026/analysis/judge_cost.py [--open 8] [--closed 3]

Scales each paid evaluator from what it actually cost per trace on the pilot (the cost
recorded on every judge call, current config only). Only E0, E0-3J, E1 and E5 call
LLM judges; E2 (process reward models on our GPUs), E3 and E4 cost nothing.

How often the judges are called depends on the model being judged: E0/E1 send a trace
to the judges when its steps do not match the reference, E5 only for milestones its
deterministic matching misses. Weak models trigger far more calls. So two rates are
shown per evaluator: the pilot's four frontier models, and Llama 3.1 70B, the weak
model. The next roster is mostly smaller open models, so the truth sits between the two.

The ROSTER column is the figure D-105 and the pilot summary quote: the open-weight roster
models at the weak model's rate and the closed-weight ones at the frontier rate. The
defaults are the roster D-110 fixed, 8 open and 3 closed. For E5 that basis is the
expensive end (a weak model leaves more milestones to the judge); for E0 it is the cheap
end (E0 judges a wrong answer only when its 20% sample draws it), so E0's full-run figure
is a floor, and the columns either side of it show the range.
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


def rates():
    """Mean judge spend per trace, per evaluator: {name: {'F': frontier rate, 'W': weak rate}}."""
    out = {}
    for name, d, _desc in EVALS:
        rows = current(d)
        per = {}
        for grp, keys in (('F', FRONTIER), ('W', WEAK)):
            rs = [r for r in rows if r['model_key'] in keys]
            per[grp] = st.mean(sum(c.get('cost_usd', 0) for c in r.get('calls') or []) for r in rs)
        out[name] = per
    return out


def roster_cost(per, n_open, n_closed):
    """The D-105 basis: open-weight models at the weak rate, closed-weight at the frontier rate."""
    return ITEMS * (n_open * per['W'] + n_closed * per['F'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--open', type=int, default=8, help='open-weight roster models (D-110: 8)')
    ap.add_argument('--closed', type=int, default=3, help='closed-weight roster models (D-110: 3)')
    a = ap.parse_args()
    models = a.open + a.closed
    n = models * ITEMS
    print('Judge cost for %d models x %d items = %d traces (pilot rates per trace)\n' % (models, ITEMS, n))
    print('%-6s %-52s %10s %10s %11s %11s %11s'
          % ('eval', 'judges', '$/trace F', '$/trace W', 'full run F', 'full run W', 'ROSTER'))
    for (name, _d, desc), per in zip(EVALS, rates().values()):
        print('%-6s %-52s %10.5f %10.5f %11.2f %11.2f %11.2f'
              % (name, desc, per['F'], per['W'], per['F'] * n, per['W'] * n,
                 roster_cost(per, a.open, a.closed)))
    print('\nF = at the pilot frontier models\' rate, W = at Llama 3.1 70B\'s rate. E2, E3, E4: $0.')
    print('ROSTER = %d open-weight models at W and %d closed-weight at F (the D-105 basis).'
          % (a.open, a.closed))


if __name__ == '__main__':
    main()
