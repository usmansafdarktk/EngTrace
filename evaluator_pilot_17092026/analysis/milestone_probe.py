"""Can we recover each frozen item's internal quantities, provably?

For every one of the 60 frozen items: reload its template, regenerate at the
recorded seed with frame-local capture, and check the produced question and
solution are byte-identical to what the freeze pinned. If they are, the captured
locals ARE that item's internals, and E3 can build milestones from them without
editing a single template - which also means the frozen text cannot move.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import json
import os
import sys

ROOT = _REPO
sys.path.insert(0, ROOT)
from tests.template_integrity.core import discover, generate   # noqa: E402

items = [json.loads(l) for l in
         open(os.path.join(ROOT, 'evaluator_pilot_17092026', 'slice', 'manifest.jsonl'), encoding='utf-8')]
refs = {r.template_id: r for r in discover()}

ok = bad = 0
per_template = {}
for it in items:
    tid = it['template_id']
    ref = refs.get(tid)
    if ref is None:
        print('NO REF %s' % tid)
        bad += 1
        continue
    inst = generate(ref, it['seed'], capture=True)
    same = (inst.question == it['question'] and inst.solution == it['solution'])
    ok += same
    bad += (not same)
    if not same:
        print('TEXT MISMATCH %s seed %s' % (it['item_id'], it['seed']))
    loc = inst.values or {}
    nums = {k: v for k, v in loc.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
    per_template.setdefault(tid, (it['item_id'], sorted(nums)))

print('\n%d of %d items reproduce byte-identically\n' % (ok, ok + bad))
for tid, (iid, names) in sorted(per_template.items()):
    print('%-42s %s' % (tid.replace('template_', ''), ', '.join(names) if names else '(no numeric locals)'))
