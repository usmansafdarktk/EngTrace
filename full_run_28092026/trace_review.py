"""An independent review of the full run's traces, before any scoring.

    python -m full_run_28092026.trace_review      # writes TRACE_REVIEW.md and trace_review.json beside this file

The decision log (D-122 to D-132) records what the run billed and how many rows ended empty. This
checks, from the trace files alone, the things that would silently corrupt the evaluation:

  integrity    every final row is a frozen item (its item_sha256 equals the manifest's), every item
               has exactly one final row, no line is malformed, no service failure is left
  sameness     one prompt hash (the pilot's), one set of request parameters, and one served model
               per model key, so no model changed under the run
  answers      among answered rows: how many stopped at the output cap with text; how many carry
               no answer marker the answer check looks for (`## Final Answer`, `Answer:`); how many
               carry none AND no number in the last 700 characters the check falls back to, the
               closest reading here of the plan's "never states one" (ANALYSIS_PLAN, D-117);
               how many leak a <think> block into the answer text
  archive      the newest local traces archive holds these exact files (SHA-256 per member)

A final row is the last row per item whose state is answered or empty, as run_traces.existing reads
it. It writes counts, model keys and template ids only: no trace text. Nothing is scored here.
"""
from __future__ import annotations

import collections
import hashlib
import json
import re
import statistics
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026.run_traces import PROMPT_SHA, TRACES, config  # noqa: E402

BACKUP = Path.home() / 'EngTrace_private_backup'
SET_ASIDE = {'qwen3.8-27b': 'D-132: the pilot generated traces with it; robustness check only'}
DONE = ('answered', 'empty')
THINK = re.compile(r'<think>|</think>', re.I)


def manifest() -> dict[str, dict]:
    rows = [json.loads(l) for l in (HERE / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()]
    return {r['item_id']: r for r in rows}


def review_model(key: str, man: dict) -> dict | None:
    path = TRACES / f'{key}.jsonl'
    if not path.exists():
        return None
    malformed, lines = 0, 0
    final, done_rows = {}, collections.Counter()
    states, modes = collections.Counter(), collections.Counter()
    prompts, requests, served, providers = (collections.Counter() for _ in range(4))
    for ln in path.read_text(encoding='utf-8').splitlines():
        lines += 1
        try:
            r = json.loads(ln)
        except json.JSONDecodeError:
            malformed += 1
            continue
        states[r.get('status')] += 1
        if r.get('status') in DONE:
            final[r['item_id']] = r
            done_rows[r['item_id']] += 1
    rows = list(final.values())
    for r in rows:
        modes[r.get('mode')] += 1
        prompts[r.get('prompt_sha256')] += 1
        requests[json.dumps(r.get('request'), sort_keys=True)] += 1
        served[r.get('served_model')] += 1
        providers[r.get('provider')] += 1
    wrong_item = [r['item_id'] for r in rows
                  if r['item_id'] not in man or r.get('item_sha256') != man[r['item_id']]['sha256']]
    answered = [r for r in rows if r['status'] == 'answered']
    capped = sum(r.get('finish_reason') == 'length' for r in answered)
    no_marker, no_answer, think = [], [], []
    for r in answered:
        t = r.get('text') or ''
        if not (answer.HEADING.search(t) or answer.ANSWER.search(t)):
            no_marker.append(r)
            if not answer.values(t[-answer.WINDOW:]):
                no_answer.append(r)
        if THINK.search(t):
            think.append(r)
    toks = [r.get('completion_tokens') or 0 for r in rows]
    billed = sum(r.get('billed_usd') or 0.0 for r in rows)
    return {
        'key': key, 'lines': lines, 'malformed': malformed, 'states_all_rows': dict(states),
        'items': len(final), 'items_missing': sorted(set(man) - set(final)),
        'items_with_two_final_rows': sum(1 for v in done_rows.values() if v > 1),
        'rows_not_a_frozen_item': len(wrong_item),
        'answered': len(answered), 'empty': len(rows) - len(answered),
        'modes': dict(modes),
        'prompt_hashes': len(prompts), 'prompt_is_pilots': set(prompts) == {PROMPT_SHA},
        'request_variants': len(requests),
        'served_models': dict(served.most_common()),
        'providers': dict(providers.most_common(4)), 'provider_count': len(providers),
        'capped_with_text': capped,
        'no_answer_marker': len(no_marker),
        'no_marker_no_number': len(no_answer),
        'no_marker_templates': dict(collections.Counter(r['template_id'] for r in no_marker).most_common(5)),
        'think_in_text': len(think),
        'empty_templates': dict(collections.Counter(r['template_id'] for r in rows
                                                    if r['status'] == 'empty').most_common(6)),
        'empty_by_template_all': dict(collections.Counter(r['template_id'] for r in rows
                                                          if r['status'] == 'empty')),
        'output_tokens_median': statistics.median(toks) if toks else 0,
        'output_tokens_max': max(toks) if toks else 0,
        'billed_usd': round(billed, 3),
    }


def archive_check(keys: list[str]) -> dict:
    zips = sorted(BACKUP.glob('full_run_traces_*.zip'))
    if not zips:
        return {'archive': None}
    z = zips[-1]
    out = {'archive': z.name, 'matched': 0, 'differ': [], 'absent': []}
    with zipfile.ZipFile(z) as zf:
        members = {Path(n).name: n for n in zf.namelist()}
        for k in keys:
            f = TRACES / f'{k}.jsonl'
            if f.name not in members:
                out['absent'].append(f.name)
                continue
            if hashlib.sha256(zf.read(members[f.name])).digest() == hashlib.sha256(f.read_bytes()).digest():
                out['matched'] += 1
            else:
                out['differ'].append(f.name)
    return out


def main() -> int:
    man = manifest()
    keys = [m['key'] for m in config()['models']]
    res = [r for r in (review_model(k, man) for k in keys) if r]
    arch = archive_check([r['key'] for r in res])
    (HERE / 'trace_review.json').write_text(json.dumps({'models': res, 'archive': arch}, indent=1) + '\n',
                                            encoding='utf-8', newline='\n')
    L = ['# Trace review: the full run before scoring', '',
         'Generated by `trace_review.py`; the checks are defined in its docstring. Counts, model keys '
         'and template ids only.', '',
         '## Integrity and sameness', '',
         '| model | items | answered | empty | missing | 2 final rows | not a frozen item | malformed | '
         'prompt | request sets | served models | providers | billed $ |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|']
    for r in res:
        tag = ' (set aside)' if r['key'] in SET_ASIDE else ''
        L.append(f"| `{r['key']}`{tag} | {r['items']} | {r['answered']} | {r['empty']} | "
                 f"{len(r['items_missing'])} | {r['items_with_two_final_rows']} | "
                 f"{r['rows_not_a_frozen_item']} | {r['malformed']} | "
                 f"{'pilot' if r['prompt_is_pilots'] else 'DIFFERS'} | {r['request_variants']} | "
                 f"{len(r['served_models'])} | {r['provider_count']} | {r['billed_usd']:.3f} |")
    roster = [r for r in res if r['key'] not in SET_ASIDE]
    L += ['', f"Roster, {len(roster)} models: {sum(r['items'] for r in roster)} final rows, "
          f"{sum(r['empty'] for r in roster)} empty, ${sum(r['billed_usd'] for r in roster):.2f} billed in the rows.",
          '', '## Answers, among answered rows', '',
          '| model | stopped at the cap with text | no answer marker | no marker and no number at the end | '
          '<think> in the text | output tokens, median / max |',
          '|---|---:|---:|---:|---:|---:|']
    for r in res:
        L.append(f"| `{r['key']}` | {r['capped_with_text']} | {r['no_answer_marker']} | "
                 f"{r['no_marker_no_number']} | {r['think_in_text']} | "
                 f"{r['output_tokens_median']:.0f} / {r['output_tokens_max']} |")
    tot = collections.Counter()
    hit = collections.Counter()
    for r in roster:
        for t, n in r['empty_by_template_all'].items():
            tot[t] += n
            hit[t] += 1
    n_empty = sum(tot.values())
    top = tot.most_common(10)
    L += ['', '## Empty rows across the roster, by template', '',
          f"{n_empty} empty rows over {len(tot)} templates. The ten with the most, each out of "
          f"{15 * len(roster)} rows ({len(roster)} models x 15 items):", '',
          '| template | empty rows | models with one |', '|---|---:|---:|']
    L += [f"| `{t}` | {n} | {hit[t]} |" for t, n in top]
    L += ['', f"These ten hold {sum(n for _, n in top)} of the {n_empty}."]
    L += ['', '## Where the empty rows and the missing markers fall', '']
    for r in res:
        if r['empty'] or r['no_answer_marker']:
            L.append(f"- `{r['key']}`: empty {r['empty_templates'] or '-'}; "
                     f"no marker {r['no_marker_templates'] or '-'}")
    served_multi = [r for r in res if len(r['served_models']) > 1]
    L += ['', '## Served models', '']
    L += [f"- `{r['key']}`: {r['served_models']}" for r in res]
    if arch.get('archive'):
        L += ['', '## Local archive', '',
              f"`{arch['archive']}`: {arch['matched']} of {len(res)} trace files identical by SHA-256; "
              f"differ {arch['differ'] or 'none'}; absent {arch['absent'] or 'none'}."]
    (HERE / 'TRACE_REVIEW.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
