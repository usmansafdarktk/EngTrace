"""Do the LLM judges catch a planted defect that every deterministic evaluator misses?

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/planted_judges.py build --sample 20,10
    $PY evaluator_pilot_17092026/analysis/planted_judges.py run --budget 4.00     # PAID
    $PY evaluator_pilot_17092026/analysis/planted_judges.py report

RESULTS_X1 Finding 7 reports that no free evaluator detects a conceptual defect: 0 of 60.
Every evaluator it tested - E3, E4, the digit rule - reads numbers, not prose, so that
result is close to a tautology and cannot carry the general claim. The evaluators that
could catch a misstated rule are the ones that read the step: E0's Tribunal, E1's panel,
E5's residual judge. This asks them.

THE PROBE IS MATCHED. Each planted defect is judged twice: once in the planted trace and
once in the SAME step of the untouched original. Same question, same gold, same step, one
character different. A judge that flags both is not detecting the defect, it is flagging
the step - and only the difference between the two is evidence. That is what makes 60
prompts worth as much here as 180 unmatched ones.

THE PROMPT IS THE FRAMEWORK'S OWN, built by its own `_tier2_tribunal_batch` with exactly
one step under review, and the reply is parsed with the framework's own logic (both reused
from analysis/judge_probe.py, which chose E1 and E5's judges the same way).

WHAT THIS ANSWERS, AND WHAT IT DOES NOT. It answers "asked directly about this step, does
the judge call it wrong". It does not answer "would E0 have sent this step to the judges",
which is a separate property of the framework's routing - E0 samples wrong-answer traces
at 0.20 and every planted trace here has a correct final answer. Both matter, and this
script is deliberately the cheaper half: it establishes whether the capability exists at
all before anyone pays to measure the routing.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import argparse
import concurrent.futures as cf
import io
import json
import random
import sys
from collections import Counter, defaultdict
from contextlib import redirect_stdout

sys.path[:0] = [_ANALYSIS, _os.path.join(_PILOT, 'evaluators')]
import judge_probe as JP  # noqa: E402

OUT = _os.path.join(_PILOT, 'scores', 'planted_judges')
PROBE = _os.path.join(OUT, 'probe_set.json')
REPLIES = _os.path.join(OUT, 'replies.jsonl')
PLANTED = _os.path.join(_ANALYSIS, 'out', 'planted', 'planted.jsonl')

# E0's own two judges. Prices per million tokens, as models.json records them.
JUDGES = [('gpt-5', 'openai/gpt-5', (1.25, 10.0)),
          ('opus-4.5', 'anthropic/claude-opus-4.5', (5.0, 25.0))]
SEED = 20260924


def _rows(path):
    return [json.loads(l) for l in open(path, encoding='utf-8') if l.strip()]


def build(sample):
    """Write the matched probe set: each defect's step, planted and original."""
    sys.path.insert(0, _os.path.join(_REPO, 'evaluation'))
    _os.environ.setdefault('HF_HUB_OFFLINE', '1')
    with redirect_stdout(io.StringIO()):
        import engtrace_evaluation_framework as fw
    from engineering_parser import extract_steps

    items = {r['item_id']: r for r in _rows(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'))}
    originals = {}
    for f in _os.listdir(_os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for r in _rows(_os.path.join(_PILOT, 'traces', f)):
                if r.get('ok'):
                    originals[(r['model_key'], r['item_id'])] = r['text']

    def tribunal_prompt(question, gt_steps, pred_steps, idx):
        captured = {}

        class Fake:
            client_openai, client_anthropic = True, None

            def _call_single_judge(self, provider, prompt):
                captured['prompt'] = prompt
                return []
        fw.EngTraceFramework._tier2_tribunal_batch(Fake(), question, gt_steps, pred_steps, [idx])
        return captured['prompt']

    def locate(steps, needle):
        """The framework splits steps its own way; find the one carrying this text."""
        flat = needle.strip()
        hits = [i for i, s in enumerate(steps) if flat and flat in s]
        if not hits:
            squashed = flat.replace(' ', '')
            hits = [i for i, s in enumerate(steps) if squashed and squashed in s.replace(' ', '')]
        return hits[0] if len(hits) >= 1 else None

    planted = [r for r in _rows(PLANTED) if r['family'] in ('arithmetic', 'conceptual')]
    rnd = random.Random(SEED)
    want = dict(zip(('conceptual', 'arithmetic'), sample))
    chosen = []
    for fam, n in want.items():
        pool = [r for r in planted if r['family'] == fam]
        rnd.shuffle(pool)
        chosen += pool[:n]

    probe, skipped = [], Counter()
    for r in chosen:
        item = items[r['item_id']]
        orig = originals.get((r['model_key'], r['item_id']))
        if orig is None:
            skipped['no original'] += 1
            continue
        gt_steps = extract_steps(item['solution'])[0]
        for tag, text, needle in (('planted', r['text'], r['after']), ('original', orig, r['before'])):
            steps = extract_steps(text)[0]
            idx = locate(steps, needle)
            if idx is None:
                skipped['step not located (%s)' % tag] += 1
                continue
            probe.append({'set_id': r['set_id'], 'family': r['family'], 'arm': tag,
                          'model_key': r['model_key'], 'item_id': r['item_id'],
                          'basis': r.get('basis'), 'before': r['before'], 'after': r['after'],
                          'step': steps[idx], 'step_index': idx,
                          'prompt': tribunal_prompt(item['question'], gt_steps, steps, idx)})
    # keep only defects whose BOTH arms built, so every comparison is matched
    by_set = defaultdict(list)
    for p in probe:
        by_set[p['set_id']].append(p)
    probe = [p for v in by_set.values() if len(v) == 2 for p in v]
    _os.makedirs(OUT, exist_ok=True)
    json.dump(probe, open(PROBE, 'w', encoding='utf-8'), indent=1)
    print('probe set: %d prompts (%d matched defects), skipped %s'
          % (len(probe), len(probe) // 2, dict(skipped) or 'none'))
    print('families:', dict(Counter(p['family'] for p in probe if p['arm'] == 'planted')))
    print('calls if every judge runs: %d' % (len(probe) * len(JUDGES)))


def run(budget):
    """Call the judges. PAID. Stops as soon as the spend would pass --budget."""
    probe = json.load(open(PROBE, encoding='utf-8'))
    done = set()
    spent = 0.0
    if _os.path.exists(REPLIES):
        for r in _rows(REPLIES):
            if r.get('ok'):
                done.add((r['judge'], r.get('set_id'), r.get('arm')))
                spent += r.get('cost_usd', 0.0)
    price = {j[0]: j[2] for j in JUDGES}
    jobs = [(j, i) for j in JUDGES for i in range(len(probe))
            if (j[0], probe[i]['set_id'], probe[i]['arm']) not in done]
    print('%d calls outstanding, %.2f already spent' % (len(jobs), spent))
    if not jobs:
        return 0

    def one(job):
        (jid, model, _p), i = job
        return dict(JP.call(model, probe[i]['prompt']), judge=jid, model=model, probe=i,
                    set_id=probe[i]['set_id'], arm=probe[i]['arm'])

    stopped = False
    with cf.ThreadPoolExecutor(max_workers=8) as pool, open(REPLIES, 'a', encoding='utf-8') as fh:
        futures = {pool.submit(one, j): j for j in jobs}
        for fut in cf.as_completed(futures):
            res = fut.result()
            pin, pout = price[res['judge']]
            res['cost_usd'] = ((res.get('in') or 0) * pin + (res.get('out') or 0) * pout) / 1e6
            spent += res['cost_usd']
            fh.write(json.dumps(res, ensure_ascii=False) + '\n')
            fh.flush()
            if spent > budget and not stopped:
                stopped = True
                print('BUDGET %.2f reached at %.2f - cancelling what has not started' % (budget, spent))
                for f in futures:
                    f.cancel()
    print('spent this run: $%.2f' % spent)
    return 0


def verdicts():
    probe = json.load(open(PROBE, encoding='utf-8'))
    out = defaultdict(dict)
    known = {p['set_id'] for p in probe}
    for r in _rows(REPLIES):
        # keyed by (judge, set, arm): the probe set grows, and an index would shift
        if not r.get('ok') or r.get('set_id') not in known:
            continue
        out[(r['judge'], r['set_id'])][r['arm']] = {
            'verdict': JP.parse(r.get('text')), 'cost': r.get('cost_usd', 0.0)}
    return probe, out


def report():
    probe, got = verdicts()
    fam = {p['set_id']: p['family'] for p in probe}
    basis = {p['set_id']: p.get('basis') for p in probe}
    print('DO THE JUDGES CATCH A PLANTED DEFECT?\n')
    print('  Matched: the same step is judged in the planted trace and in the untouched')
    print('  original. A defect is CAUGHT when the judge returns an Error category on the')
    print('  planted arm and not on the original. "Other" is the judge declining to')
    print('  classify; it is reported beside, never folded in.\n')
    print('  %-10s %-12s %5s %8s %8s %11s %10s'
          % ('judge', 'family', 'n', 'caught', '+other', 'flags orig', 'both arms'))
    total = 0.0
    for jid, _m, _p in JUDGES:
        for f in ('conceptual', 'arithmetic'):
            sets = [s for (j, s) in got if j == jid and fam.get(s) == f and len(got[(j, s)]) == 2]
            if not sets:
                continue
            caught = other = orig = both = 0
            for s in sets:
                a = got[(jid, s)]
                pv, ov = _wrong(a['planted']['verdict']), _wrong(a['original']['verdict'])
                po = _wrong(a['planted']['verdict'], True)
                total += a['planted']['cost'] + a['original']['cost']
                caught += pv and not ov
                other += po and not pv
                orig += ov
                both += pv and ov
            n = len(sets)
            print('  %-10s %-12s %5d %8.3f %8.3f %11.3f %10d'
                  % (jid, f, n, caught / n, (caught + other) / n, orig / n, both))
    print('\n  spend on the calls behind this table: $%.2f' % total)

    print('\n  WHAT THE JUDGES ANSWERED, by arm')
    for jid, _m, _p in JUDGES:
        for arm in ('planted', 'original'):
            cats = Counter(_category(got[(j, s)][arm]['verdict'])
                           for (j, s) in got if j == jid and arm in got[(j, s)])
            print('    %-10s %-9s %s' % (jid, arm, dict(cats)))

    print('\n  BY BASIS - a defect inside a claim the checker parses is one the digit rule')
    print('  already finds; everything else is where a judge is the only candidate.')
    for jid, _m, _p in JUDGES:
        for b in ('claim', 'milestone', None):
            sets = [s for (j, s) in got if j == jid and basis.get(s) == b and len(got[(j, s)]) == 2]
            if not sets:
                continue
            caught = sum(_wrong(got[(jid, s)]['planted']['verdict'])
                         and not _wrong(got[(jid, s)]['original']['verdict']) for s in sets)
            print('    %-10s %-12s %3d sets  caught %.3f'
                  % (jid, b or 'conceptual', len(sets), caught / len(sets)))

    print('\n  EITHER JUDGE - the panel reading, which is how E0 votes')
    for f in ('conceptual', 'arithmetic'):
        sets = {s for (j, s) in got if fam.get(s) == f}
        ok = [s for s in sets if all((j, s) in got and len(got[(j, s)]) == 2 for j, _m, _p in JUDGES)]
        if not ok:
            continue
        caught = sum(any(_wrong(got[(j, s)]['planted']['verdict'])
                         and not _wrong(got[(j, s)]['original']['verdict'])
                         for j, _m, _p in JUDGES) for s in ok)
        print('    %-12s %3d sets  caught by at least one judge %.3f' % (f, len(ok), caught / len(ok)))


def _wrong(parsed, count_other=False):
    """Did the judge call the step under review wrong, in the framework's own vocabulary?

    The Tribunal answers with a category, not a yes/no: `Alternative Correct` (a different
    but valid route - the verdict RESULTS_E0 found it gives to 93% of what it is sent),
    `Calculation Error`, `Conceptual Error`, or `Other`. Only the two Error categories are
    a flag. `Other` is counted separately, because it is the judge declining to classify
    rather than declining to flag, and folding it either way would overstate something.
    """
    vals = parsed if isinstance(parsed, list) else [parsed] if parsed else []
    for v in vals:
        cat = (v.get('category') if isinstance(v, dict) else None) or ''
        low = cat.strip().lower()
        if 'error' in low:
            return True
        if count_other and low == 'other':
            return True
    return False


def _category(parsed):
    vals = parsed if isinstance(parsed, list) else [parsed] if parsed else []
    for v in vals:
        if isinstance(v, dict) and v.get('category'):
            return v['category']
    return 'unparsed'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('build', 'run', 'report'))
    ap.add_argument('--sample', default='20,10', help='conceptual,arithmetic defects to probe')
    ap.add_argument('--budget', type=float, default=4.0, help='stop the run at this spend')
    a = ap.parse_args()
    if a.cmd == 'build':
        return build([int(x) for x in a.sample.split(',')])
    if a.cmd == 'run':
        return run(a.budget)
    return report()


if __name__ == '__main__':
    raise SystemExit(main())
