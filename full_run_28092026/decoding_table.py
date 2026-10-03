"""The decoding settings and token use of the full run, per model, printed from the traces (Appendix P; next steps A8).

    python -m full_run_28092026.decoding_table            # FREE: writes DECODING_TABLE.md and results/decoding_table.json

WHY. The run used each provider's default decoding (D-122): the request set only the token ceiling and OpenRouter's
routing preferences. Four roster models returned no reasoning tokens (D-168), and the paper must state per model what
it ran with, so the closed tier's scores are not read as a ceiling. Everything here is read from the trace rows the
harness wrote (one per item: the request, the served model, the provider, token usage, the finish reason) and from
models.json; nothing is typed in and nothing is called.

WHAT IT PRINTS, per roster model and for the set-aside model on its own line: the configured and served ids, the
weights, the request (the ceiling and the routing preferences; no temperature, top-p or reasoning setting was sent),
the endpoints that served the rows, the finish reasons, prompt and completion tokens (median, 90th percentile, max),
reasoning tokens as the endpoint reported them (median, 90th percentile, the share of rows above zero) and the share
of rows that returned reasoning text, and the billed cost. A reasoning count of zero is what the endpoint reported;
with no reasoning text returned either, it is read as the model having run without extended reasoning at that
provider's default (D-168).
"""
from __future__ import annotations

import ast
import collections
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER, SET_ASIDE  # noqa: E402

OUT_MD = HERE / 'DECODING_TABLE.md'
OUT_JSON = HERE / 'results' / 'decoding_table.json'


def pctl(v: list, q: float) -> float | None:
    if not v:
        return None
    s = sorted(v)
    k = max(0, min(len(s) - 1, round(q * (len(s) - 1))))
    return s[k]


def one(key: str, cfg: dict, variant: str = 'main') -> dict:
    rows = [r for r in score.trace_rows(variant, key) if r.get('status') in ('answered', 'empty')]
    req = collections.Counter(str(r.get('request')) for r in rows)
    request = ast.literal_eval(req.most_common(1)[0][0]) if req else {}
    prov = request.get('extra_body', {}).get('provider', {}) if isinstance(request, dict) else {}
    comp = [r.get('completion_tokens') or 0 for r in rows]
    prompt = [r.get('prompt_tokens') or 0 for r in rows]
    reas = [r.get('reasoning_tokens') for r in rows]
    reas_num = [x or 0 for x in reas]
    return {
        'model_key': key, 'weights': cfg.get('weights'), 'configured': cfg.get('model'),
        'served': dict(collections.Counter(r.get('served_model') for r in rows)),
        'providers': dict(collections.Counter(r.get('provider') for r in rows).most_common()),
        'rows': len(rows), 'answered': sum(r['status'] == 'answered' for r in rows),
        'empty': sum(r['status'] == 'empty' for r in rows),
        'modes': dict(collections.Counter(r.get('mode') for r in rows)),
        'prompt_sha256': sorted({r.get('prompt_sha256') for r in rows}),
        'request_variants': len(req),
        'max_tokens': request.get('max_tokens') if isinstance(request, dict) else None,
        'provider_sort': prov.get('sort'), 'quantizations': prov.get('quantizations'),
        'allow_fallbacks': prov.get('allow_fallbacks'),
        'sampling_parameters_sent': sorted(k for k in (request if isinstance(request, dict) else {})
                                           if k not in ('max_tokens', 'extra_body')),
        'reasoning_parameter': (request.get('extra_body', {}).get('reasoning') if isinstance(request, dict) else None),
        'finish_reasons': dict(collections.Counter(str(r.get('finish_reason')) for r in rows)),
        'prompt_tokens': {'median': statistics.median(prompt) if prompt else None},
        'completion_tokens': {'median': statistics.median(comp) if comp else None, 'p90': pctl(comp, 0.9),
                              'max': max(comp) if comp else None},
        'reasoning_tokens': {'median': statistics.median(reas_num) if reas_num else None, 'p90': pctl(reas_num, 0.9),
                             'max': max(reas_num) if reas_num else None,
                             'share_above_zero': sum(x > 0 for x in reas_num) / len(rows) if rows else None,
                             'unreported_rows': sum(x is None for x in reas)},
        'reasoning_text_share': sum(bool((r.get('reasoning') or '').strip()) for r in rows) / len(rows) if rows else None,
        'billed_usd': round(sum(r.get('billed_usd') or 0.0 for r in rows), 3),
        'seconds_median': statistics.median([r.get('seconds') or 0.0 for r in rows]) if rows else None,
        'tool': tool_use(rows),
    }


def tool_use(rows: list[dict]) -> dict | None:
    """The tool arm's use of its tool (D-184), from the rows' tool fields; None for any other arm."""
    t = [r for r in rows if r.get('tool_calls') is not None]
    if not t:
        return None
    calls = [r['tool_calls'] for r in t]
    return {'rows': len(t), 'rows_with_calls': sum(c > 0 for c in calls), 'share_with_calls': sum(c > 0 for c in calls) / len(t),
            'calls_median': statistics.median(calls), 'calls_p90': pctl(calls, 0.9), 'calls_max': max(calls),
            'turns_median': statistics.median([r.get('turns') or 1 for r in t]),
            'scripts': sum(calls), 'scripts_failed': sum(r.get('tool_errors') or 0 for r in t),
            'scripts_refused': sum(r.get('tool_refused') or 0 for r in t), 'limit_rows': sum(bool(r.get('tool_limit')) for r in t)}


def table(res: list[dict]) -> list[str]:
    def f(v, d=0):
        return '-' if v is None else (f'{v:.{d}f}' if isinstance(v, float) else str(v))
    L = ['# Decoding settings and token use, per model (Appendix P)', '',
         'Generated by `decoding_table.py` from the trace rows and `models.json`; the fields are defined in its docstring. '
         'Every model ran at its provider\'s default decoding: the request carried the token ceiling and OpenRouter\'s routing '
         'preferences only, and no temperature, top-p or reasoning setting. A reasoning-token count is what the endpoint '
         'reported. Rows are the final row per item (answered or empty).', '',
         '| model | weights | served id | ceiling | routing | endpoints (rows) | finish: stop / length / error / none | prompt tokens, median | completion tokens: median / p90 / max | reasoning tokens: median / p90 / max | rows with reasoning tokens | rows with reasoning text | billed $ |',
         '|---|---|---|---:|---|---|---|---:|---|---|---:|---:|---:|']
    for r in res:
        fr = r['finish_reasons']
        prov = '; '.join(f'{p} ({n})' for p, n in list(r['providers'].items())[:4]) + (f'; +{len(r["providers"]) - 4} more' if len(r['providers']) > 4 else '')
        routing = f"sort={r['provider_sort']}" + (', fp8-or-better' if r['quantizations'] else '') + (', fallbacks' if r['allow_fallbacks'] else '')
        ct, rt = r['completion_tokens'], r['reasoning_tokens']
        L.append(f"| `{r['model_key']}` | {r['weights']} | {', '.join(f'`{s}`' for s in r['served'])} | {r['max_tokens']} | {routing} | {prov} | "
                 f"{fr.get('stop', 0)} / {fr.get('length', 0)} / {fr.get('error', 0)} / {fr.get('None', 0)} | {f(r['prompt_tokens']['median'])} | "
                 f"{f(ct['median'])} / {f(ct['p90'])} / {f(ct['max'])} | {f(rt['median'])} / {f(rt['p90'])} / {f(rt['max'])} | "
                 f"{f(rt['share_above_zero'], 3)} | {f(r['reasoning_text_share'], 3)} | {r['billed_usd']:.2f} |")
    L.append('')
    only_roster = all(r['model_key'] in ROSTER + SET_ASIDE for r in res)          # an arm's anchors count too (D-182)
    none_reasoning = [r['model_key'] for r in res if (r['model_key'] in ROSTER or not only_roster) and r['reasoning_tokens']['share_above_zero'] == 0]
    L.append(f"Models whose endpoint reported no reasoning tokens on any row: {', '.join(f'`{k}`' for k in none_reasoning) or 'none'}. "
             'The other roster models returned reasoning tokens on nearly every row. One prompt hash per run; its value is in the JSON.')
    L.append('')
    return L


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--variant', default='main', help='a traces/<variant>/ arm: reasoning-<effort> (C1), flagship and flagship-reasoning-<effort> (C3), openbook (C4)')
    a = ap.parse_args()
    cfg = {m['key']: m for m in json.loads((HERE / 'models.json').read_text(encoding='utf-8'))['models']}
    if a.variant == 'main':
        keys, out_md, out_json = ROSTER + SET_ASIDE, OUT_MD, OUT_JSON
    else:
        keys = [k for k in list(ROSTER) + [k for k in cfg if k not in ROSTER]     # the roster first, then an arm's anchors (C3)
                if score.trace_path(a.variant, k).exists()]
        out_md = HERE / f'DECODING_TABLE_{a.variant}.md'
        out_json = HERE / 'results' / f'decoding_table_{a.variant}.json'
        if not keys:
            raise SystemExit(f'no traces for variant {a.variant}')
    res = [one(k, cfg[k], a.variant) for k in keys]
    out_json.parent.mkdir(exist_ok=True)
    out_json.write_text(json.dumps(res, indent=1), encoding='utf-8')
    lines = table(res)
    if a.variant != 'main':
        rp = res[0]['reasoning_parameter'] if res else None
        item, dec = (('C1', 'D-179, D-180') if a.variant.startswith('reasoning-') else
                     ('C3', 'D-182') if a.variant.startswith('flagship') else
                     ('C4', 'D-183') if a.variant == 'openbook' else
                     ('C4', 'D-184') if a.variant == 'tool' else ('a variant', 'D-141'))
        what = (f'the same request as the main run plus the reasoning parameter `{json.dumps(rp)}`' if rp else
                'the same request as the main run plus the python tool (`tools`, `tool_choice: auto`); tokens and the bill are '
                'summed over a row' + chr(39) + 's turns' if a.variant == 'tool' else
                'the same request as the main run, at the provider' + chr(39) + 's default' +
                (', on the question with the reference equations appended' if a.variant == 'openbook' else ''))
        lines[0] = f'# Decoding settings and token use, the `{a.variant}` arm ({item})'
        lines[2] = lines[2].replace('preferences only, and no temperature, top-p or reasoning setting.',
                                    'preferences, no temperature or top-p, and what the line below names for this arm.')
        lines.insert(2, f'The arm `{a.variant}` ({dec}): {what}; rows are the arm' + chr(39) + 's items only. '
                        'The main run' + chr(39) + 's table is `DECODING_TABLE.md`.')
        lines.insert(3, '')
    if a.variant == 'tool':
        lines += ['', '## Tool use', '',
                  '| model | rows | rows with a tool call | calls per row: median / p90 / max | model turns, median | scripts | failed | refused | rows at the call limit |',
                  '|---|---:|---:|---|---:|---:|---:|---:|---:|']
        for r in res:
            t = r.get('tool')
            if t:
                lines.append(f"| `{r['model_key']}` | {t['rows']} | {t['rows_with_calls']} | {t['calls_median']:.0f} / {t['calls_p90']:.0f} / {t['calls_max']} | "
                             f"{t['turns_median']:.0f} | {t['scripts']} | {t['scripts_failed']} | {t['scripts_refused']} | {t['limit_rows']} |")
    out_md.write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(lines))
    return 0


if __name__ == '__main__':
    sys.exit(main())
