"""Would E0 ever show a judge the planted step? Tier 1, run for real, judges stubbed.

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/planted_routing.py [--limit N]

Finding 8 establishes that E0's judges catch a third of planted conceptual defects and four
fifths of the arithmetic ones when a step is put in front of them. It deliberately does not
establish that E0 would put it there. This does, and it costs nothing: the framework's Tier
1 runs for real on local scorer models, and `e0_tribunal.setup(dry_run=True)` stubs the
judges, so every prompt the Tribunal would have sent is recorded and none is paid for.

Three things are measured per planted defect, each against the SAME trace unmodified:

  triggered   the framework decided this trace needs the Tribunal at all
  reached     it had mismatched steps to send, so a judge would have been called
  shown       the PLANTED step is among the steps in a prompt - the only case in which
              the detection rate in Finding 8 can apply

The third is the one that matters. A trace can reach the Tribunal over some other step and
still never expose the defect, and the judges' 0.350 on conceptual defects only converts
into an end-to-end detection rate on the traces where the corrupted step is actually shown.

Every planted trace here has a CORRECT final answer, so the framework's wrong-answer
sampling (D3, probability 0.20) is not what routes them: whatever routing happens comes
from Tier 1 failing to match the step.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import argparse
import json
import sys
import time
from collections import Counter, defaultdict

sys.path[:0] = [_PILOT, _os.path.join(_PILOT, 'evaluators')]
import e0_tribunal as E0  # noqa: E402
import run_evaluator as RE  # noqa: E402

PLANTED = _os.path.join(_ANALYSIS, 'out', 'planted', 'planted.jsonl')
OUT = _os.path.join(_PILOT, 'scores', 'planted_routing')
ROWS = _os.path.join(OUT, 'routing.jsonl')


def _rows(path):
    return [json.loads(l) for l in open(path, encoding='utf-8') if l.strip()]


def framework_steps(text):
    """The step list the framework itself builds - the index space its mismatches use."""
    from engineering_parser import extract_steps
    return extract_steps(text)[0]


def locate(steps, needle):
    flat = (needle or '').strip()
    if not flat:
        return None
    hits = [i for i, st in enumerate(steps) if flat in st]
    if not hits:
        squashed = flat.replace(' ', '')
        hits = [i for i, st in enumerate(steps) if squashed in st.replace(' ', '')]
    return hits[0] if hits else None


def capture_indices(state):
    """Record WHICH steps the Tribunal is asked about, not just how many.

    e0_tribunal's dry run replaces the Tribunal with a recorder that logs
    `len(mismatch_indices)`. The indices themselves are what decides whether the planted
    step is exposed, so this replaces that recorder on our own state object - the
    evaluator's source is untouched, and the substitute keeps writing the same log entry.
    """
    seen = []

    def recorder(question, gt_steps, pred_steps, mismatch_indices):
        idx = list(mismatch_indices)
        seen.append(idx)
        state.log.prompts.append({'provider': 'DRY', 'mismatch_steps': len(idx),
                                  'prompt_chars': len(question) + 1500})
        return {}, 'N/A', {}

    state.framework._tier2_tribunal_batch = recorder
    return seen


def shown(sent, text, needle):
    """Is the planted step among the steps Tier 1 could not match?

    Indices are in the framework's own step space, so the step is located by its text
    rather than by the index the planting used.
    """
    flat = (needle or '').strip()
    steps = framework_steps(text)
    hits = [i for i, st in enumerate(steps) if flat and flat in st]
    if not hits:
        squashed = flat.replace(' ', '')
        hits = [i for i, st in enumerate(steps) if squashed and squashed in st.replace(' ', '')]
    idx = hits[0] if hits else None
    flat_sent = sorted({i for batch in sent for i in batch})
    return (idx is not None and idx in flat_sent), idx, flat_sent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, help='score only the first N defects (both arms)')
    a = ap.parse_args()

    items = {r['item_id']: r for r in _rows(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'))}
    originals = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for r in _rows(_os.path.join(_PILOT, 'traces', f)):
                if r.get('ok'):
                    originals[(r['model_key'], r['item_id'])] = r['text']

    planted = [r for r in _rows(PLANTED) if r['family'] in ('arithmetic', 'conceptual')]
    if a.limit:
        planted = planted[:a.limit]

    _os.makedirs(OUT, exist_ok=True)
    done = {(r['set_id'], r['arm']) for r in _rows(ROWS)} if _os.path.exists(ROWS) else set()
    state = E0.setup(dry_run=True)
    sent = capture_indices(state)
    print('scoring %d defects x 2 arms, judges stubbed (no call is made)' % len(planted))

    t0 = time.time()
    with open(ROWS, 'a', encoding='utf-8') as fh:
        for n, r in enumerate(planted, 1):
            item = items[r['item_id']]
            for arm, text, needle in (('planted', r['text'], r['after']),
                                      ('original', originals.get((r['model_key'], r['item_id'])), r['before'])):
                if text is None or (r['set_id'], arm) in done:
                    continue
                seed = RE.seed_for('e0', r['item_id'], r['model_key'])
                sent.clear()
                out = E0.score(state, item, {'text': text}, seed)
                was_shown, idx, mism = shown(sent, text, needle)
                fh.write(json.dumps({
                    'set_id': r['set_id'], 'arm': arm, 'family': r['family'],
                    'basis': r.get('basis'), 'model_key': r['model_key'], 'item_id': r['item_id'],
                    'triggered': bool(out['triggered']),
                    'reached': bool(out['meta'].get('tribunal_reached_judges')),
                    'shown': was_shown, 'step_index': idx, 'mismatch_steps': mism,
                    'n_prompts': len(out.get('dry_prompts') or []),
                    'final_answer_acc': out['scores'].get('final_answer_acc'),
                }, ensure_ascii=False) + '\n')
                fh.flush()
            if n % 10 == 0:
                print('  %d/%d defects, %.1f min elapsed' % (n, len(planted), (time.time() - t0) / 60))
    report()


def report():
    rows = _rows(ROWS)
    by = defaultdict(dict)
    for r in rows:
        by[r['set_id']][r['arm']] = r
    sets = {k: v for k, v in by.items() if len(v) == 2}
    fam = {k: v['planted']['family'] for k, v in sets.items()}
    basis = {k: v['planted'].get('basis') for k, v in sets.items()}

    print('\nWOULD E0 SHOW THE PLANTED STEP TO A JUDGE?  (%d matched defects)\n' % len(sets))
    print('  %-12s %5s %11s %11s %11s %14s'
          % ('family', 'n', 'triggered', 'reached', 'step shown', 'shown, original'))
    for f in ('conceptual', 'arithmetic'):
        ks = [k for k in sets if fam[k] == f]
        if not ks:
            continue
        n = len(ks)
        print('  %-12s %5d %11.3f %11.3f %11.3f %14.3f'
              % (f, n,
                 sum(sets[k]['planted']['triggered'] for k in ks) / n,
                 sum(sets[k]['planted']['reached'] for k in ks) / n,
                 sum(sets[k]['planted']['shown'] for k in ks) / n,
                 sum(sets[k]['original']['shown'] for k in ks) / n))

    print('\n  BY BASIS')
    for b in ('claim', 'milestone', None):
        ks = [k for k in sets if basis[k] == b]
        if not ks:
            continue
        print('    %-12s %3d  shown %.3f' % (b or 'conceptual', len(ks),
                                             sum(sets[k]['planted']['shown'] for k in ks) / len(ks)))

    # end to end: routing x the judges' detection rate from Finding 8
    print('\n  END TO END - routing x detection (Finding 8: either judge catches 0.350 of')
    print('  conceptual defects, 0.800 of arithmetic, given the step is shown)')
    for f, rate in (('conceptual', 0.350), ('arithmetic', 0.800)):
        ks = [k for k in sets if fam[k] == f]
        if not ks:
            continue
        s = sum(sets[k]['planted']['shown'] for k in ks) / len(ks)
        print('    %-12s shown %.3f x caught %.3f = %.3f of defects detected by E0 end to end'
              % (f, s, rate, s * rate))

    ans = Counter(r['final_answer_acc'] for r in rows if r['arm'] == 'planted')
    print('\n  final_answer_acc on the planted arm: %s  (the planted answer is unchanged, so'
          % dict(ans))
    print('  anything below 1.0 is E0\'s own answer check, not the plant - see D-098)')


if __name__ == '__main__':
    main()
