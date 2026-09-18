"""What the main EngTrace benchmark costs per model, from measured pilot numbers.

Generation: measured $/trace from the 300 pilot traces.
Judging:    measured from the E0 dry run - what fraction of traces reach the
            Tribunal, how big those prompts are, three judges each.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import glob
import json
import os

P = _PILOT_DIR
FULL = 2250                      # 150 templates x 15 instances, the main benchmark

gen = {}
for f in glob.glob(os.path.join(P, 'traces', '*.jsonl')):
    rows = [json.loads(l) for l in open(f, encoding='utf-8')]
    ok = [r for r in rows if r['ok']]
    key = os.path.basename(f)[:-6]
    gen[key] = sum(r.get('cost_usd', 0) for r in ok) / max(len(ok), 1)

dry = [json.loads(l) for f in glob.glob(os.path.join(P, 'scores', 'e0_dry', '*.jsonl'))
       for l in open(f, encoding='utf-8')]
reach = [r for r in dry if r['meta'].get('tribunal_reached_judges')]
chars = [p['prompt_chars'] for r in reach for p in r.get('dry_prompts', [])]
in_tok = (sum(chars) / len(chars)) / 3.5
share = len(reach) / len(dry)
JUDGE_OUT = 400
JUDGES = {'gpt-5': (1.25, 10.0), 'opus-4.5': (5.0, 25.0), 'gemini-3.1-pro': (2.0, 12.0)}
per_judged = sum((in_tok * pin + JUDGE_OUT * pout) / 1e6 for pin, pout in JUDGES.values())
judge_per_trace = share * per_judged

print('MEASURED FROM THE PILOT')
print('  %.1f%% of traces reach the Tribunal, ~%.0f input tokens each, 3 judges'
      % (share * 100, in_tok))
print('  E0 judging: $%.4f per trace scored (before judge reasoning tokens)' % judge_per_trace)
print()
print('COST OF THE FULL BENCHMARK (%d items) PER EVALUATED MODEL' % FULL)
print('%-18s %10s %12s %12s %12s' % ('model', '$/trace', 'generation', 'E0 judging', 'total'))
for k, v in sorted(gen.items(), key=lambda kv: -kv[1]):
    g = v * FULL
    j = judge_per_trace * FULL
    print('%-18s %10.5f %12.0f %12.0f %12.0f' % (k, v, g, j, g + j))
print('%-18s %10s %12s %12.0f %12s' % ('open-weight on own GPU', '0', '0',
                                       judge_per_trace * FULL, 'judging only'))
print()
print('  judging alone, x5 for judge reasoning tokens: $%.0f per model' % (judge_per_trace * FULL * 5))
print()
for cap in (500,):
    print('UNDER A $%d CAP' % cap)
    j = judge_per_trace * FULL
    print('  with E0 (LLM judges): %d open-weight models fit (judging %.0f each, generation free)'
          % (cap // j, j))
    print('  with E0, judge reasoning at stage-1 rates (x5): %d models' % (cap // (j * 5)))
    print('  with E3 (deterministic, no judge in the loop): judging is $0, so the roster is')
    print('    limited by generation only - open-weight on your own GPUs costs nothing per model.')
