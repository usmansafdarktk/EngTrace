"""Print a fixed random sample of the inconsistent claims E4 found, for reading.

`arith_consistency` counts every inconsistent claim in a trace, and most of them
are not tied to a milestone - so they were never read when the milestone
contradictions were. This prints a reproducible sample (seed 7) so each can be
classified by hand as a real slip in the trace or a checker false positive. The
classification is recorded in RESULTS_E4.md; re-run this to check it.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)

import glob
import json
import random
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT_DIR, 'scores', 'e4', '*.jsonl'))
        for l in open(f, encoding='utf-8')]
latest = max(rows, key=lambda r: r['ts'])['config_sha256']
rows = [r for r in rows if r['config_sha256'] == latest and not r.get('error')]
pool = [(r['model_key'], r['item_id'], c) for r in rows for c in r['meta']['inconsistent_claims']]
print('config %s: %d inconsistent claims across %d traces; sample of %d (seed 7)\n'
      % (latest[:12], len(pool), len(rows), min(N, len(pool))))
random.seed(7)
for i, (m, iid, c) in enumerate(random.sample(pool, min(N, len(pool))), 1):
    print('%2d. %-15s %-32s %s = %s' % (i, m, iid, c['left'][:60], c['right'][:30]))
    print('    evaluates %s vs %s' % (['%.5g' % x for x in c['left_value']],
                                     ['%.5g' % x for x in c['right_value']]))
