import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import json, os, re, sys
from collections import Counter
ROOT = _REPO
sys.path.insert(0, os.path.join(ROOT, 'evaluation'))
from engineering_parser import extract_steps
NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
items = [json.loads(l) for l in open(os.path.join(ROOT, 'evaluator_pilot_17092026/slice/manifest.jsonl'), encoding='utf-8')]
by_t = {}
for it in items:
    ans = [l for l in it['solution'].splitlines() if '**Answer:**' in l]
    line = ans[-1] if ans else ''
    nums = re.findall(NUM, line.replace(',', ''))
    _, _, gt_final = extract_steps(it['solution'])
    t = it['template_id'].replace('template_', '')
    last = float(nums[-1]) if nums else None
    by_t.setdefault(t, []).append((it['answer_type'], len(nums), gt_final, last, line[:150]))
bad_t = 0
for t, rows in by_t.items():
    typ = rows[0][0]
    multi = sum(1 for r in rows if r[1] > 1)
    first_not_last = sum(1 for r in rows if r[1] > 1 and r[2] != r[3])
    none = sum(1 for r in rows if r[1] == 0)
    flag = 'READS WRONG NUMBER' if first_not_last else ('no number in answer' if none else '')
    bad_t += bool(first_not_last)
    print('%-34s %-14s multi-number %d/4  parsed!=last %d/4  %s' % (t, typ, multi, first_not_last, flag))
    if first_not_last: print('     e.g. %r' % rows[0][4])
print('\ntemplates where E0 reads a different number than the answer line ends on: %d of %d' % (bad_t, len(by_t)))
