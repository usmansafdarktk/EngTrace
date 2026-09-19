"""Build E2's input file: every trace, split into the steps E0 uses.

    python evaluator_pilot_17092026/hpg/build_e2_inputs.py OUT.jsonl

One row per trace: the item's question, the trace's steps (split with E0's own
`extract_steps`, so step k is the same step in every evaluator and the judge probe's
labelled step indices apply unchanged), and the hashes that tie the row back to the
frozen item and to the exact trace text. Gold five plus the gemma robustness column.
"""
import glob
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = os.path.dirname(HERE)
ROOT = os.path.dirname(PILOT)
sys.path.insert(0, os.path.join(ROOT, 'evaluation'))
from engineering_parser import extract_steps  # noqa: E402

MODELS = ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b', 'gemma-4-31b')


def main(out):
    items = {json.loads(l)['item_id']: json.loads(l)
             for l in open(os.path.join(PILOT, 'slice', 'manifest.jsonl'), encoding='utf-8')}
    rows = []
    for mk in MODELS:
        best = {}
        for l in open(os.path.join(PILOT, 'traces', mk + '.jsonl'), encoding='utf-8'):
            r = json.loads(l)
            if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
                best[r['item_id']] = r
        for iid, r in sorted(best.items()):
            steps = [s.strip() for s in extract_steps(r['text'])[0]]
            steps = [s for s in steps if s]
            rows.append({'item_id': iid, 'model_key': mk, 'item_sha256': items[iid]['sha256'],
                         'trace_sha256': hashlib.sha256(r['text'].encode('utf-8')).hexdigest(),
                         'question': items[iid]['question'], 'steps': steps})
    with open(out, 'w', encoding='utf-8', newline='\n') as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
    n = [len(r['steps']) for r in rows]
    print('%d traces, %d steps (min %d, max %d per trace) -> %s' % (len(rows), sum(n), min(n), max(n), out))


if __name__ == '__main__':
    main(sys.argv[1])
