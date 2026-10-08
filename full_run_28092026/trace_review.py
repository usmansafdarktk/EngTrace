"""An independent review of the full run's traces, before any scoring.

    python -m full_run_28092026.trace_review                       # writes TRACE_REVIEW.md and trace_review.json beside this file
    python -m full_run_28092026.trace_review --variant paraphrase  # the same checks on a variant's traces (D-150)

For a variant the item hash is checked against that variant's items: a repeat's rows against the pool
manifest, a paraphrase's rows against paraphrase/manifest.jsonl, with each row's original_sha256
against the pool manifest; the served model per key is compared with the main run's, since the arms
must be answered by the same checkpoint; and the report goes to TRACE_REVIEW_<variant>.md. The
archive check applies to the main run only.

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


def variant_manifest(variant: str, man: dict) -> dict[str, dict]:
    """What a variant's rows must hash to: the pool's entries for a repeat's 300 items (D-151; the whole
    pool before, which counted the other 1,950 as missing), the passing paraphrases for the paraphrase
    arm (with the original's hash beside each)."""
    from full_run_28092026 import subsamples
    from full_run_28092026 import run_traces
    if variant == run_traces.FULL_REASONING:                 # the matched-settings run: every pool item
        return man
    if variant in run_traces.CASCADE_REPEATS:                # the cascade repeats: the 300 repeat items
        return {i: man[i] for i in subsamples.repeat_ids()}
    if variant == run_traces.CASCADE_PARAPHRASE:             # the cascade paraphrases: the kept pairs, as the harness runs them
        keep = run_traces.kept_pairs()
        return {i: r for i, r in variant_manifest('paraphrase', man).items() if i in keep}
    if variant in subsamples.REPEAT_VARIANTS:
        return {i: man[i] for i in subsamples.repeat_ids()}
    if variant.startswith('reasoning-') or variant.startswith('flagship') or variant == 'tool':   # C1, C3, C4 tool: the subsample's originals
        return {i: man[i] for i in subsamples.paraphrase_ids()}
    if variant in ('openbook', 'openbook2'):             # C4: the modified questions, by their manifest (version 1 or 2)
        p = HERE / 'openbook' / ('manifest.jsonl' if variant == 'openbook' else 'manifest2.jsonl')
        if not p.exists():
            raise SystemExit(f'{p.name} is missing: run openbook.py --build first')
        return {r['item_id']: {'sha256': r['sha256'], 'original_sha256': r['original_sha256']}
                for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())}
    if variant != 'paraphrase':
        return man
    p = HERE / 'paraphrase' / 'manifest.jsonl'
    if not p.exists():
        raise SystemExit('paraphrase/manifest.jsonl is missing: the paraphrase arm has not been written')
    out = {}
    for r in map(json.loads, p.read_text(encoding='utf-8').splitlines()):
        if r['passed']:
            out[r['item_id']] = {'sha256': r['sha256'], 'original_sha256': r['original_sha256']}
    return out


def review_model(key: str, man: dict, variant: str = 'main', main_served: dict | None = None) -> dict | None:
    path = TRACES / f'{key}.jsonl' if variant == 'main' else TRACES / variant / f'{key}.jsonl'
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
                  if r['item_id'] not in man or r.get('item_sha256') != man[r['item_id']]['sha256']
                  or ('original_sha256' in man[r['item_id']]
                      and r.get('original_sha256') != man[r['item_id']]['original_sha256'])]
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
    tool_rows = [r for r in rows if r.get('tool_calls') is not None]        # the tool arm (D-184)
    tool = None
    if tool_rows:
        calls = [r['tool_calls'] for r in tool_rows]
        tool = {'rows': len(tool_rows), 'rows_with_calls': sum(c > 0 for c in calls),
                'calls_median': statistics.median(calls), 'calls_max': max(calls),
                'turns_median': statistics.median([r.get('turns') or 1 for r in tool_rows]),
                'scripts': sum(calls), 'scripts_failed': sum(r.get('tool_errors') or 0 for r in tool_rows),
                'scripts_refused': sum(r.get('tool_refused') or 0 for r in tool_rows),
                'limit_rows': sum(bool(r.get('tool_limit')) for r in tool_rows)}
    return {
        'tool': tool,
        'key': key, 'lines': lines, 'malformed': malformed, 'states_all_rows': dict(states),
        'items': len(final), 'items_missing': sorted(set(man) - set(final)),
        'items_with_two_final_rows': sum(1 for v in done_rows.values() if v > 1),
        'rows_not_a_frozen_item': len(wrong_item),
        'answered': len(answered), 'empty': len(rows) - len(answered),
        'modes': dict(modes),
        'prompt_hashes': len(prompts), 'prompt_is_pilots': set(prompts) == {PROMPT_SHA},
        'request_variants': len(requests),
        'served_models': dict(served.most_common()),
        'served_as_main': (set(served) == set(main_served[key]) if main_served and main_served.get(key) else None),
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
        'items_per_template': len(rows) / max(1, len({r['template_id'] for r in rows})),
        'output_tokens_median': statistics.median(toks) if toks else 0,
        'output_tokens_max': max(toks) if toks else 0,
        'billed_usd': round(billed, 3),
    }


def archive_check(keys: list[str]) -> dict:
    # the main run's archives are dated, `full_run_traces_<date>.zip`; the repeats' and the paraphrase
    # arm's carry a name before the date and are not this check's (D-168)
    zips = sorted(p for p in BACKUP.glob('full_run_traces_*.zip')
                  if p.stem.removeprefix('full_run_traces_').replace('-', '').isdigit())
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
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--variant', default='main')
    a = ap.parse_args()
    variant = a.variant
    man = variant_manifest(variant, manifest())
    keys = [m['key'] for m in config()['models']                 # an arm reviews whichever models it holds (the anchors, C3)
            if variant == 'main' or (TRACES / variant / f"{m['key']}.jsonl").exists()]
    main_served = None
    if variant != 'main':
        main_served = {k: (review_model(k, manifest()) or {}).get('served_models', {}) for k in keys}
    res = [r for r in (review_model(k, man, variant, main_served) for k in keys) if r]
    if not res:
        raise SystemExit(f'no traces for variant {variant}')
    arch = archive_check([r['key'] for r in res]) if variant == 'main' else {'archive': None}
    suffix = '' if variant == 'main' else f'_{variant}'
    (HERE / f'trace_review{suffix}.json').write_text(json.dumps({'variant': variant, 'models': res, 'archive': arch},
                                                                indent=1) + '\n', encoding='utf-8', newline='\n')
    L = [f'# Trace review: the {"full run" if variant == "main" else variant + " variant"} before scoring', '',
         'Generated by `trace_review.py`; the checks are defined in its docstring. Counts, model keys '
         'and template ids only.', '',
         '## Integrity and sameness', '',
         '| model | items | answered | empty | missing | 2 final rows | not a frozen item | malformed | '
         'prompt | request sets | served models | as main | providers | billed $ |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---|---:|---:|']
    for r in res:
        tag = ' (set aside)' if r['key'] in SET_ASIDE else ''
        same = '' if r['served_as_main'] is None else ('yes' if r['served_as_main'] else 'NO')
        L.append(f"| `{r['key']}`{tag} | {r['items']} | {r['answered']} | {r['empty']} | "
                 f"{len(r['items_missing'])} | {r['items_with_two_final_rows']} | "
                 f"{r['rows_not_a_frozen_item']} | {r['malformed']} | "
                 f"{'pilot' if r['prompt_is_pilots'] else 'DIFFERS'} | {r['request_variants']} | "
                 f"{len(r['served_models'])} | {same} | {r['provider_count']} | {r['billed_usd']:.3f} |")
    roster = [r for r in res if r['key'] not in SET_ASIDE]
    L += ['', f"Roster, {len(roster)} models: {sum(r['items'] for r in roster)} final rows, "
          f"{sum(r['empty'] for r in roster)} empty, ${sum(r['billed_usd'] for r in roster):.2f} billed in the rows."
          + ('' if variant == 'main' else f' Items expected per model: {len(man)}.'),
          '', '## Answers, among answered rows', '',
          '| model | stopped at the cap with text | no answer marker | no marker and no number at the end | '
          '<think> in the text | output tokens, median / max |',
          '|---|---:|---:|---:|---:|---:|']
    for r in res:
        L.append(f"| `{r['key']}` | {r['capped_with_text']} | {r['no_answer_marker']} | "
                 f"{r['no_marker_no_number']} | {r['think_in_text']} | "
                 f"{r['output_tokens_median']:.0f} / {r['output_tokens_max']} |")
    if any(r.get('tool') for r in res):
        L += ['', '## Tool use (the `tool` arm, D-184)', '',
              'Rows are the final row per item; a script is one tool call; "refused" is the static filter, "failed" a script that '
              'raised or ran out of time; a row at the call limit was asked to answer with the tool withheld.', '',
              '| model | rows | rows with a tool call | calls per row, median / max | model turns, median | scripts | failed | refused | rows at the call limit |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in res:
            t = r.get('tool')
            if t:
                L.append(f"| `{r['key']}` | {t['rows']} | {t['rows_with_calls']} | {t['calls_median']:.0f} / {t['calls_max']} | "
                         f"{t['turns_median']:.0f} | {t['scripts']} | {t['scripts_failed']} | {t['scripts_refused']} | {t['limit_rows']} |")
    tot = collections.Counter()
    hit = collections.Counter()
    for r in roster:
        for t, n in r['empty_by_template_all'].items():
            tot[t] += n
            hit[t] += 1
    n_empty = sum(tot.values())
    per_t = statistics.median([r['items_per_template'] for r in res]) if res else 15
    top = tot.most_common(10)
    L += ['', '## Empty rows across the roster, by template', '',
          f"{n_empty} empty rows over {len(tot)} templates. The ten with the most, each out of "
          f"{round(per_t) * len(roster)} rows ({len(roster)} models x {round(per_t)} items):", '',
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
    (HERE / f'TRACE_REVIEW{suffix}.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
