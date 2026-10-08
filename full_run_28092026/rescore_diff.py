"""What the re-scores at ANSWER FINAL changed in each store (WS-G, step G2).

    python -m full_run_28092026.rescore_diff [--variant V ...]                 # FREE: against each store's last archive
    python -m full_run_28092026.rescore_diff --since 20261008T000000Z          # against the store before the first
                                                                               # re-score at or after that time
    python -m full_run_28092026.rescore_diff --name floor_fix                  # writes RESCORE_DIFF_floor_fix.md and
                                                                               # results/rescore_diff_floor_fix.json

For each store `repair_round5 --rescore` re-scored, the store score.py archived (by default the variant's last logged
re-score in scores/_replaced/round5_rescore.jsonl; with --since, the first logged at or after that time, so the
comparison spans every re-score since) against the store now, model file by model file, row by row on item_id:

  verdicts    answered rows usable in both stores whose headline label moved, from -> to, per model and template. The
              cause: a verdict raised to correct while its parts do not all match by number (`matched` < `of`) is the
              symbolic step's (D5), which raises to correct and nothing else; every other move is the number rule's:
              the bounds of D11b and D11d, which only lower, and the answer targets the milestone matcher feeds
              (answer.targets), which the matcher's fixes can move either way. symbolic/rescore_preview.py attributes
              the same verdicts the same way.
  scores      the mean score per model before and after, over all rows (unusable 0), as score.py stores it.
  E3          answered rows whose reached list changed (D11a), the milestones newly reached and no longer reached, and
              of those rows the ones that still miss a milestone: their judge prompt (judge.job) changed, so judge.py
              sends them again.
  the rest    every other field that differs, and the answer and E3 sub-fields that differ. Rows whose `unusable` flag
              flipped are counted apart from the verdicts.

The counts are rescore_preview.py's, so the main store's eleven models can be set against RESCORE_PREVIEW.md; the
set-aside qwen3.8-27b has its own line. Nothing is written to a store.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

LOG = score.SCORES / '_replaced' / 'round5_rescore.jsonl'
OUT_MD = HERE / 'RESCORE_DIFF.md'
OUT_JSON = HERE / 'results' / 'rescore_diff.json'
SCORE_OF = score.SCORE_OF


def archived_store(variant: str, since: str | None = None) -> Path | None:
    """The store as score.py archived it: at the variant's last logged re-score, or with `since` (a UTC stamp like
    20261008T000000Z) at the first logged re-score at or after it."""
    entries = [json.loads(x) for x in LOG.read_text(encoding='utf-8').splitlines() if x.strip()]
    mine = [e for e in entries if e['variant'] == variant and e.get('store_archived_to')]
    if since:
        mine = [e for e in mine if e['rescored_at_utc'] >= since][:1]
    return HERE / mine[-1]['store_archived_to'] if mine else None


def rows_of(path: Path) -> dict:
    return {r['item_id']: r for r in map(json.loads, path.read_text(encoding='utf-8').splitlines())} if path.exists() else {}


def answer_py(cfg: dict, name: str = 'answer.py') -> str | None:
    """The evaluator file's LF-normalised SHA-256 (12 characters) that a store's CONFIG.json records."""
    return next((v['sha256_lf'][:12] for k, v in cfg.get('evaluators', {}).items() if k.endswith('/' + name)), None)


def diff_model(old: dict, new: dict) -> dict:
    out = {'rows_old': len(old), 'rows_new': len(new), 'only_one_side': len(set(old) ^ set(new)),
           'moves': collections.Counter(), 'by_template': collections.defaultdict(collections.Counter),
           'unusable_flips': 0, 'e3_rows': 0, 'e3_gained': 0, 'e3_lost': 0, 'judge_prompts_changed': 0, 'prompt_items': [],
           'e3_items': [],
           'fields': collections.Counter(), 'answer_fields': collections.Counter(), 'e3_fields': collections.Counter(),
           'mean_old': sum(r['score'] for r in old.values()) / len(old) if old else None,
           'mean_new': sum(r['score'] for r in new.values()) / len(new) if new else None}
    for i in sorted(set(old) & set(new)):
        o, n = old[i], new[i]
        out['fields'].update(k for k in set(o) | set(n) if o.get(k) != n.get(k))
        oa, na = o.get('answer') or {}, n.get('answer') or {}
        out['answer_fields'].update(k for k in set(oa) | set(na) if oa.get(k) != na.get(k))
        oe, ne = o.get('e3') or {}, n.get('e3') or {}
        out['e3_fields'].update(k for k in set(oe) | set(ne) if oe.get(k) != ne.get(k))
        if n['status'] != 'answered':
            continue
        was, now = oe.get('reached') or [], ne.get('reached') or []
        if was != now:
            out['e3_rows'] += 1
            out['e3_items'].append(i)
            out['e3_gained'] += sum(a and not b for a, b in zip(now, was))
            out['e3_lost'] += sum(b and not a for a, b in zip(now, was))
            if not all(now):
                out['judge_prompts_changed'] += 1
                out['prompt_items'].append(i)
        if o['unusable'] != n['unusable']:
            out['unusable_flips'] += 1
            continue
        if n['unusable']:
            continue
        a, b = oa['label'], na['label']
        if a != b:
            cause = 'symbolic step' if b == 'correct' and (na.get('matched') or 0) < (na.get('of') or 0) else 'number rule'
            out['moves'][(a, b, cause)] += 1
            out['by_template'][n['template_id']][cause] += 1
    return out


def diff_store(variant: str, since: str | None = None) -> dict:
    prev = archived_store(variant, since)
    store = score.SCORES / variant
    if prev is None or not prev.exists():
        return {'variant': variant, 'error': 'no logged re-score with an archived store'}
    cfg_old = json.loads((prev / 'CONFIG.json').read_text(encoding='utf-8'))
    cfg_new = json.loads((store / 'CONFIG.json').read_text(encoding='utf-8'))
    files = sorted({f.name for f in prev.glob('*.jsonl')} | {f.name for f in store.glob('*.jsonl')})
    models = {f[:-6]: diff_model(rows_of(prev / f), rows_of(store / f)) for f in files}
    # a second judge's sample (judge.py --judge ... --sample, e.g. e5_grok-4-6): its sent rows whose E3 changed, so that
    # MiMo's prompt changed or, where E3 now reaches every milestone, no call is made
    second = {}
    for d in sorted(x for x in store.iterdir() if x.is_dir() and x.name.startswith('e5_')):
        second[d.name] = {k: sorted(i for i, r in rows_of(d / f'{k}.jsonl').items()
                                    if r.get('sent') and i in set(m['e3_items'])) for k, m in models.items()}
        second[d.name + ' (no call now)'] = {k: sorted(set(v) - set(models[k]['prompt_items']))
                                            for k, v in second[d.name].items()}
    return {'variant': variant, 'archived': prev.relative_to(HERE).as_posix(),
            'answer_py_old': answer_py(cfg_old), 'answer_py_new': answer_py(cfg_new),
            'milestones_py_old': answer_py(cfg_old, 'milestones.py'), 'milestones_py_new': answer_py(cfg_new, 'milestones.py'),
            'commit_new': cfg_new.get('git', '')[:7], 'dirty_new': cfg_new.get('dirty'), 'models': models,
            'second_judge_prompts_changed': second}


def total(models: dict, keys) -> dict:
    t = {'moves': collections.Counter(), 'by_template': collections.defaultdict(collections.Counter),
         'fields': collections.Counter()}
    for k in keys:
        m = models[k]
        t['moves'].update(m['moves'])
        t['fields'].update(m['fields'])
        for tid, c in m['by_template'].items():
            t['by_template'][tid].update(c)
        for f in ('unusable_flips', 'e3_rows', 'e3_gained', 'e3_lost', 'judge_prompts_changed', 'only_one_side'):
            t[f] = t.get(f, 0) + m[f]
    return t


def raised(moves) -> int:
    return sum(n for (a, b, _c), n in moves.items() if SCORE_OF[b] > SCORE_OF[a])


def lowered(moves) -> int:
    return sum(n for (a, b, _c), n in moves.items() if SCORE_OF[b] < SCORE_OF[a])


def render(stores: list[dict]) -> list[str]:
    L = ['# The re-score at ANSWER FINAL: what changed in each store', '',
         'Generated by `rescore_diff.py` (its docstring defines the counts): each store against the copy score.py '
         'archived at its re-score. Verdicts: answered rows usable in both stores; raised verdicts are the symbolic '
         'step\'s, lowered ones the number rule\'s. E3: answered rows whose reached milestones changed, and of them the '
         'judge prompts that changed (a milestone still missed).', '',
         '| store | answer.py before / after | milestones.py before / after | models | verdicts moved | raised | lowered | '
         'unusable flipped | E3 rows changed | judge prompts changed | rows in one side only | fields that differ |',
         '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for s in stores:
        if 'error' in s:
            L.append(f"| `{s['variant']}` | {s['error']} | | | | | | | | | | |")
            continue
        t = total(s['models'], s['models'])
        L.append(f"| `{s['variant']}` | {s['answer_py_old']} / {s['answer_py_new']} | "
                 f"{s.get('milestones_py_old')} / {s.get('milestones_py_new')} | {len(s['models'])} | "
                 f"{sum(t['moves'].values())} | {raised(t['moves'])} | {lowered(t['moves'])} | {t['unusable_flips']} | "
                 f"{t['e3_rows']} | {t['judge_prompts_changed']} | {t['only_one_side']} | "
                 f"{', '.join(f'{k} {n}' for k, n in sorted(t['fields'].items())) or '-'} |")
    main = next((s for s in stores if s['variant'] == 'main' and 'error' not in s), None)
    if main:
        eleven = [k for k in ROSTER if k in main['models']]
        t = total(main['models'], eleven)
        L += ['', f'## The main store, the {len(eleven)} models of the paper (set against `symbolic/RESCORE_PREVIEW.md`)', '',
              '| from | to | cause | verdicts |', '|---|---|---|---:|']
        for (a, b, c), n in sorted(t['moves'].items(), key=lambda kv: (-kv[1], kv[0])):
            L.append(f'| {a} | {b} | {c} | {n} |')
        L += [f"| all | | | {sum(t['moves'].values())} |", '',
              '| model | verdicts moved | raised | lowered | mean score before | mean score after | change | '
              'E3 rows changed | newly reached | no longer reached | judge prompts changed |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for k in eleven + sorted(set(main['models']) - set(eleven)):
            m = main['models'][k]
            mark = '' if k in eleven else ' (set aside)'
            L.append(f"| `{k}`{mark} | {sum(m['moves'].values())} | {raised(m['moves'])} | {lowered(m['moves'])} | "
                     f"{m['mean_old']:.4f} | {m['mean_new']:.4f} | {m['mean_new'] - m['mean_old']:+.4f} | {m['e3_rows']} | "
                     f"{m['e3_gained']} | {m['e3_lost']} | {m['judge_prompts_changed']} |")
        L.append(f"| all {len(eleven)} | {sum(t['moves'].values())} | {raised(t['moves'])} | {lowered(t['moves'])} | | | | "
                 f"{t['e3_rows']} | {t['e3_gained']} | {t['e3_lost']} | {t['judge_prompts_changed']} |")
        L += ['', '| template | moved by the symbolic step | moved by the number rule |', '|---|---:|---:|']
        for tid, c in sorted(t['by_template'].items(), key=lambda kv: (-sum(kv[1].values()), kv[0])):
            L.append(f"| `{tid.replace('template_', '')}` | {c['symbolic step']} | {c['number rule']} |")
        sj = main.get('second_judge_prompts_changed', {})
        for d, per in ((d, per) for d, per in sj.items() if not d.endswith('(no call now)')):
            hit = {k: v for k, v in per.items() if v}
            gone = sum(len(v) for v in sj.get(d + ' (no call now)', {}).values())
            L += ['', f"The second judge's sample (`{d}`): {sum(len(v) for v in per.values())} of its sent rows have E3 "
                  f"matches that changed, so MiMo's prompt is no longer the one the second judge saw ({gone} of them now "
                  'reach every milestone and make no call)'
                  + (': ' + ', '.join(f'`{k}` ' + ', '.join(f'`{i}`' for i in v) for k, v in hit.items()) if hit else '') + '.']
        L += ['', 'Sub-fields that differ on the main store\'s rows (all models): answer '
              + (', '.join(f'{k} {n}' for k, n in sorted(sum((m['answer_fields'] for m in main['models'].values()),
                                                                  collections.Counter()).items())) or 'none')
              + '; e3 ' + (', '.join(f'{k} {n}' for k, n in sorted(sum((m['e3_fields'] for m in main['models'].values()),
                                                                          collections.Counter()).items())) or 'none') + '.']
    return L + ['']


def as_json(stores: list[dict]) -> dict:
    def conv(m):
        return {**{k: v for k, v in m.items() if k not in ('moves', 'by_template', 'fields', 'answer_fields', 'e3_fields')},
                'moves': [{'from': a, 'to': b, 'cause': c, 'n': n} for (a, b, c), n in sorted(m['moves'].items())],
                'by_template': {t: dict(c) for t, c in sorted(m['by_template'].items())},
                'fields': dict(m['fields']), 'answer_fields': dict(m['answer_fields']), 'e3_fields': dict(m['e3_fields'])}
    return {'stores': [{**{k: v for k, v in s.items() if k != 'models'},
                        'models': {k: conv(m) for k, m in s.get('models', {}).items()}} for s in stores]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--variant', action='append', help='a store; repeatable (default: every store the log names)')
    ap.add_argument('--since', help='UTC stamp (e.g. 20261008T000000Z): the baseline is the store before the first '
                                    're-score at or after it (default: before the last re-score)')
    ap.add_argument('--name', help='write RESCORE_DIFF_<name>.md and results/rescore_diff_<name>.json instead')
    a = ap.parse_args()
    logged = list(dict.fromkeys(json.loads(x)['variant'] for x in LOG.read_text(encoding='utf-8').splitlines() if x.strip()))
    stores = [diff_store(v, a.since) for v in (a.variant or logged)]
    L = render(stores)
    if a.since:
        L.insert(2, f'Baseline: each store before its first re-score at or after {a.since}.')
        L.insert(3, '')
    out_md = OUT_MD.with_name(f'RESCORE_DIFF_{a.name}.md') if a.name else OUT_MD
    out_json = OUT_JSON.with_name(f'rescore_diff_{a.name}.json') if a.name else OUT_JSON
    out_md.write_text('\n'.join(L), encoding='utf-8', newline='\n')
    out_json.write_text(json.dumps({**as_json(stores), 'since': a.since}, indent=1) + '\n', encoding='utf-8')
    sys.stdout.reconfigure(errors='replace')
    print('\n'.join(L))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
