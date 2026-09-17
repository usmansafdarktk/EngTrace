"""Check the produced traces against the frozen slice, before anyone annotates them.

    python -m evaluator_pilot_17092026.verify_traces
    python -m evaluator_pilot_17092026.verify_traces --selftest

`run_traces.py --status` counts rows. A count cannot tell you that the rows are
attached to the right items, that every model answered the same questions, or that
a model quietly returned the same text sixty times. Those are the ways a trace set
is complete and worthless at once, which is this repository's recurring defect, so
they are checked here rather than assumed.

THE PROPERTIES.

  T1  Coverage.  Every model has a row for every one of the 60 frozen items, and
      no rows for anything else.  A missing row is a hole in the comparison; an
      extra one means the manifest moved under the run.

  T2  Anchoring.  Every row's item_sha256 equals the manifest's hash for that
      item.  This is what makes the freeze mean something after the fact: an
      annotation is attached to text whose hash is recorded on the trace itself.

  T3  One prompt.  Every row carries the same prompt_sha256, and it is the hash of
      the template in evaluation/run_inference.py.  Different prompts would make
      the models incomparable with each other and with the archive.

  T4  Non-empty.  Every row marked ok has text.  A reasoning model that spends its
      budget thinking returns content="" with tokens billed; counted as a success
      that is a trace nobody can annotate.

  T5  Distinctness.  Within a model, no two items share output text.  Identical
      text across different questions means a stuck generation, a cached response,
      or rows written against the wrong item.

  T6  One checkpoint per column.  Every row of a model reports the same
      served_model.  A column holding two checkpoints under one name is not one
      model's performance, and the split is rarely random.

  T7  Answered, not just truncated.  No row ends at finish_reason=length.  Such a
      row HAS text, so it survives T4, and it stops partway through the
      derivation - which for a step-level annotation is worse than nothing,
      because it looks like a trace.  It also means the column is measuring our
      token ceiling rather than the model.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_traces as rt  # noqa: E402

NATURAL_STOP = {'stop', 'end_turn', 'eos', None}


def load():
    items = {r['item_id']: r for r in rt.manifest()}
    cfg = rt.config()
    traces = {}
    for spec in cfg['models']:
        path = rt.trace_path(spec['key'])
        if not os.path.exists(path):
            traces[spec['key']] = []
            continue
        rows = [json.loads(ln) for ln in open(path, encoding='utf-8')]
        # A retried item has several rows; the answer is the one that succeeded.
        best = {}
        for r in rows:
            if r['ok'] or r['item_id'] not in best:
                best[r['item_id']] = r
        traces[spec['key']] = list(best.values())
    return items, traces


def check(items, traces):
    bad = []
    for key, rows in traces.items():
        ok = [r for r in rows if r['ok']]
        ids = {r['item_id'] for r in ok}

        missing = set(items) - ids
        extra = ids - set(items)
        if missing:
            bad.append('T1 %s: %d items have no answer (%s)'
                       % (key, len(missing), ', '.join(sorted(missing)[:3])))
        if extra:
            bad.append('T1 %s: %d rows are not in the frozen slice (%s)'
                       % (key, len(extra), ', '.join(sorted(extra)[:3])))

        for r in ok:
            it = items.get(r['item_id'])
            if it and r.get('item_sha256') != it['sha256']:
                bad.append('T2 %s/%s: item hash does not match the manifest'
                           % (key, r['item_id']))

        for r in ok:
            if not (r.get('text') or '').strip():
                bad.append('T4 %s/%s: ok row with empty text' % (key, r['item_id']))

        seen = defaultdict(list)
        for r in ok:
            seen[(r.get('text') or '').strip()].append(r['item_id'])
        for text, who in seen.items():
            if len(who) > 1 and text:
                bad.append('T5 %s: %d items share identical output (%s)'
                           % (key, len(who), ', '.join(sorted(who)[:3])))

        served = {r.get('served_model') for r in ok}
        if len(served) > 1:
            bad.append('T6 %s: column holds %d checkpoints: %s'
                       % (key, len(served), sorted(map(str, served))))

    for key, rows in traces.items():
        cut = [r for r in rows if r['ok'] and r.get('finish_reason') not in NATURAL_STOP]
        if cut:
            bad.append('T7 %s: %d of %d answers end at finish_reason=%s (%s)'
                       % (key, len(cut), sum(1 for r in rows if r['ok']),
                          cut[0].get('finish_reason'),
                          ', '.join(sorted(r['item_id'] for r in cut)[:3])))

    prompts = {r.get('prompt_sha256') for rows in traces.values() for r in rows if r['ok']}
    if len(prompts) > 1:
        bad.append('T3: %d different prompts across the run: %s'
                   % (len(prompts), sorted(str(p)[:12] for p in prompts)))
    elif prompts and rt.PROMPT_SHA not in prompts:
        bad.append('T3: traces were produced with prompt %s, the deployed one is now %s'
                   % (list(prompts)[0][:12], rt.PROMPT_SHA[:12]))
    return bad


def report(items, traces):
    print('%-16s %6s %8s %9s %9s %10s   %s'
          % ('model', 'items', 'cost $', 'out tok', 'reasoning', 'truncated', 'provider'))
    for key, rows in traces.items():
        ok = [r for r in rows if r['ok']]
        if not ok:
            print('%-16s %6d' % (key, 0))
            continue
        ct = sorted(r.get('completion_tokens') or 0 for r in ok)
        reas = [len(r.get('reasoning') or '') for r in ok]
        trunc = sum(1 for r in ok if r.get('finish_reason') not in NATURAL_STOP)
        prov = Counter(r.get('provider') for r in ok if r.get('provider'))
        print('%-16s %6d %8.3f %9d %9d %10s   %s'
              % (key, len(ok), sum(r.get('cost_usd', 0) for r in ok),
                 ct[len(ct) // 2], sorted(reas)[len(reas) // 2], trunc or '-',
                 ', '.join('%s %d' % kv for kv in prov.most_common(3)) or '-'))
    total = sum(r.get('cost_usd', 0) for rows in traces.values() for r in rows if r['ok'])
    n = sum(1 for rows in traces.values() for r in rows if r['ok'])
    print('%-16s %6d %8.3f' % ('TOTAL', n, total))
    print('\nagainst %d frozen items x %d models = %d expected'
          % (len(items), len(traces), len(items) * len(traces)))


def selftest():
    """Each property must fail on damage aimed at it, and not on the others' damage."""
    items, traces = load()
    if not any(traces.values()):
        print('selftest needs traces on disk')
        return 1
    key = next(k for k, v in traces.items() if v)

    def damaged(fn):
        copy = {k: [dict(r) for r in v] for k, v in traces.items()}
        fn(copy)
        return check(items, copy)

    plants = {
        'T1 dropped row': lambda t: t[key].pop(),
        'T1 foreign row': lambda t: t[key].append(dict(t[key][0], item_id='not_an_item#9')),
        'T2 wrong item hash': lambda t: t[key][0].update(item_sha256='0' * 64),
        'T3 second prompt': lambda t: t[key][0].update(prompt_sha256='f' * 64),
        'T4 empty text': lambda t: t[key][0].update(text='   '),
        'T5 duplicate output': lambda t: t[key][1].update(text=t[key][0]['text']),
        'T6 two checkpoints': lambda t: t[key][0].update(served_model='some/other-model'),
        'T7 truncated answer': lambda t: t[key][0].update(finish_reason='length'),
    }
    clean = check(items, traces)
    bad = []
    if clean:
        bad.append('the live traces do not pass: %s' % clean[0])
    for label, fn in plants.items():
        found = damaged(fn)
        tag = label.split()[0]
        if not any(f.startswith(tag) for f in found):
            bad.append('%s was not caught' % label)
        print('  [%s] %s' % ('ok' if any(f.startswith(tag) for f in found) else 'FAIL', label))
    for x in bad:
        print('  - ' + x)
    print('selftest: %d failure(s)' % len(bad))
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    items, traces = load()
    report(items, traces)
    bad = check(items, traces)
    print()
    for x in bad:
        print('  ' + x)
    print('%d trace-set failure(s)' % len(bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
