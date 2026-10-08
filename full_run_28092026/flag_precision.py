"""Arithmetic flags on correct answers: carried precision or other (WS-C3 step 2; review 1 M3, review 2 W10).

    python -m full_run_28092026.flag_precision [--quick]   # writes results/flag_precision.json and
                                                          # results/sections/flag_precision.md
    python -m full_run_28092026.flag_precision --selftest  # the classification on hand-made chains; writes nothing

WHAT IT ASKS. The digit rule (arith.py, ok_digit) recomputes a claim's left side from the numbers the response
shows and accepts the right side only as a correct rounding at the precision it displays, each shown operand
allowed to move within its own last digit. A response that rounds an intermediate for display but keeps
computing with the unrounded value can still be flagged: the printed result is right for the value it carried,
not for the operand it printed. This script separates those flags from the rest, on every response whose final
answer is correct, per model.

THE RULE. A flagged claim is CARRIED PRECISION when, after replacing rounded operands by unrounded upstream values
that appear earlier in the response, the printed right side is a correct rounding of the recomputed left side
(arith.agree_any_digit at the right side's displayed precision, under the same unit factors, with no further
widening). Otherwise it is OTHER, split three ways: `truncated` (the printed right side is the left side as shown
cut toward zero at its displayed precision, |r| <= |l| < |r| + u, rather than rounded), `no_upstream` (no operand
has an earlier, more precise value) and `not_reproduced` (one has, and the recomputation still misses the printed
digits).
  - An operand is a literal written with a decimal point or an exponent (arith.ROUNDED_LITERAL), or a whole number
    of 100 or more (arith.py takes every bare integer as exact, so `1235` carried from 1234.6 is flagged; smaller
    integers are counts and coefficients). Its upstream candidates are values that appear before the claim, lie
    within one unit of the operand's last digit (a rounding, a truncation or a double rounding of it) and carry
    finer precision: a number written with more digits (read with arith.normalise, its precision from its own last
    digit), or the value of an earlier claim's left side when that side is an expression (precision unbounded).
    Half a unit would not do: the digit rule already lets every shown operand move half a unit of its last digit
    (arith.shown_uncertainty), so a flag that carried precision explains lies beyond that.
  - "Before the claim": the response text before the claim's step, the step's lines above the claim, and the
    claim's line up to its left side; claims earlier in the walk order.
  - Each operand takes its most recent candidate, or its most precise; every subset of the operands with a
    candidate is tried (up to 4 operands; beyond that, all of them at once). The right side's operands are
    replaced too when it is an expression rather than one number.
HOW. arith.py's Claim keeps the first 80 characters of each side, so this file walks a step with a local copy of
arith.check and arith.check_part that keeps the whole segments; every parsing and judging function is arith.py's
own, imported. Every step the walk reads must give the stored claim count and digit-flag count
(score.py's `steps` rows) or it is counted as unverified and left out.

THE EXPERT'S FLAGS. FLAG_REVIEW_3.md: a domain expert read 190 sampled flags and confirmed 171 as slips. The
sample (scores/flag_review/round3/sample.jsonl and the filled sample.csv; local, gitignored, never printed here) is
re-found in the main store's traces by model, item, step index, line and the claim's two sides, and classified by
the same rule. A claim the current rule no longer flags, or a step whose text has changed, is counted apart.

Reads the score stores, the traces (each answered row's hash checked) and the local review sample. Calls no
model. `--quick` (or ENGTRACE_QUICK=1) draws 1,000 resamples for the template intervals instead of 10,000.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import csv
import datetime as dt
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import arith  # noqa: E402
import e2_prm  # noqa: E402

from full_run_28092026 import score  # noqa: E402
from full_run_28092026.clause_variants import load_rows, quick_mode, stores  # noqa: E402

RESULTS = HERE / 'results'
OUT_JSON = RESULTS / 'flag_precision.json'
OUT_MD = RESULTS / 'sections' / 'flag_precision.md'
REVIEW = HERE / 'scores' / 'flag_review' / 'round3'
MAX_SUBSET = 4


# ------------------------------------------------------------------ the walk (arith.check, segments kept whole)

def walk(step: str) -> list[dict]:
    """arith.check on one step, each claim with its whole segments, its line, where its left side starts on the
    normalised line, the flags its segments were evaluated with, and both verdicts."""
    claims, prev = [], ''
    for i, raw in enumerate(step.splitlines()):
        continued = bool(arith.CONTINUED.search(arith.EMPHASIS_END.sub('', prev))
                         or ('=' not in prev and re.match(r'\s*[-+](?=[(\d\\])', prev)
                             and re.match(r'\s*[-+](?=[(\\])', raw)))
        prev = '' if not raw.strip() or arith.MARKDOWN_RULE.fullmatch(raw) else raw
        if '=' not in raw:
            continue
        norm = arith.normalise(raw)
        for ci, clause in enumerate(arith.CLAUSE.split(norm)):
            for ri, part in enumerate(arith.RELATION.split(clause)):
                segs = re.split(r'(?<![<>!=])=(?!=)', arith.LABEL_EQ.sub('≡', part))
                if len(segs) < 2:
                    continue
                threshold = ri > 0 and bool(re.fullmatch(r'\s*[-+]?\d+(?:\.\d+)?\s*', segs[0]))
                _part(claims, i, norm, segs, (continued and ci == 0 and ri == 0) or threshold)
    return claims


def _part(claims: list, i: int, norm: str, segs: list, first_is_not_a_value: bool) -> None:
    """arith.check_part, line for line, keeping what the classification needs."""
    chains, nums = [], []
    for k, sg in enumerate(segs):
        first, last = k == 0, k == len(segs) - 1
        if first_is_not_a_value and k == 0:
            v, kind, unit = [], 'unparseable', ''
        else:
            v, kind, unit = arith.evaluate(sg, first=first, last=last)
        if v:
            nums.append((sg.strip(), v, unit, first, last))
        elif kind == 'unparseable':
            chains.append(nums)
            nums = []
    chains.append(nums)
    for nums in chains:
        inherited, nxt = [], ''
        for sg, v, u, f, l in reversed(nums):
            nxt = u or nxt
            inherited.append((sg, v, u or nxt, bool(u), f, l))
        nums = list(reversed(inherited))
        for (ls, lv, lu, lown, lf, ll), (rs, rv, ru, rown, rf, rl) in zip(nums, nums[1:]):
            implied = arith.implied_factors(ru) if rown and not lown else set()
            u = arith.displayed_ulp(rs)
            lefts = [x for x in (lv, arith.cut_middle(ls)) if x and len(x) == len(lv)]
            rights = [x for x in (rv, arith.cut_middle(rs)) if x and len(x) == len(rv)]
            if u is None:
                digit = True
            else:
                digit = any(arith.agree_any_digit(l, r, lu, ru, u, implied=implied) for l in lefts for r in rights)
                if not digit:
                    ul, ur = arith.shown_uncertainty(ls, lv), arith.shown_uncertainty(rs, rv)
                    digit = any(arith.agree_any_digit(l, r, lu, ru, u, ul, ur, implied=implied)
                                for l in lefts for r in rights)
            claims.append({'line': i, 'at': max(norm.find(ls), 0), 'norm': norm, 'ls': ls, 'rs': rs, 'lv': lv,
                           'rv': rv, 'lu': lu, 'ru': ru, 'lf': lf, 'll': ll, 'rf': rf, 'rl': rl,
                           'implied': implied, 'u': u, 'rights': rights, 'ok_digit': digit,
                           'ok': arith.agree_any(lv, rv, lu, ru, implied)})


# ------------------------------------------------------------------ upstream values and the classification

MIN_INTEGER = 100


def literals(seg: str) -> list[tuple]:
    """The operands of a segment that may be roundings: (start, end, text, value, ulp). A bare integer below
    MIN_INTEGER is a count or a coefficient and is skipped."""
    out = []
    for m in arith.ROUNDED_LITERAL.finditer(seg):
        lit = m.group(0)
        h = arith.ulp(lit)
        try:
            x = float(lit)
        except ValueError:
            continue
        if '.' not in lit and 'e' not in lit.lower() and abs(x) < MIN_INTEGER:
            continue
        if h:
            out.append((m.start(), m.end(), lit, x, h))
    return out


def written_numbers(norm_text: str) -> list[tuple[float, float]]:
    """Every number written in normalised text, as (|value|, the place value of its last digit)."""
    out = []
    for m in arith.ROUNDED_LITERAL.finditer(norm_text):
        h = arith.ulp(m.group(0))
        try:
            out.append((abs(float(m.group(0))), h if h is not None else 1.0))
        except ValueError:
            pass
    return out


def computed_values(c: dict) -> list[tuple[float, float]]:
    """An earlier claim's left side, when it is an expression: its values, precision unbounded."""
    written = [abs(float(x)) for x in arith.NUM_LITERAL.findall(c['ls'].replace(',', ''))
               if re.fullmatch(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', x)]
    return [(abs(v), 0.0) for v in c['lv']
            if isinstance(v, float) and np.isfinite(v) and not any(abs(abs(v) - w) <= 1e-12 * max(1.0, w) for w in written)]


def candidates(x: float, h: float, upstream: list[tuple[float, float]]) -> list[tuple[int, float, float]]:
    ax = abs(x)
    return [(k, y, uy) for k, (y, uy) in enumerate(upstream)
            if uy < h * (1 - 1e-9) and y != ax and abs(y - ax) < h * (1 - 1e-9)]


def substituted(seg: str, lits: list[tuple], picks: dict[int, float], first: bool, last: bool, n: int):
    """The segment re-evaluated with the chosen literals replaced (sign kept); None when it no longer evaluates."""
    s = seg
    for j in sorted(picks, reverse=True):
        a, b, lit, x, _h = lits[j]
        s = s[:a] + '(' + repr(-picks[j] if x < 0 else picks[j]) + ')' + s[b:]
    vals, kind, _unit = arith.evaluate(s, first=first, last=last)
    return vals if kind == 'num' and len(vals) == n else None


def truncated(c: dict) -> bool:
    """Is the printed right side the left side cut (toward zero) at the precision it displays, rather than
    rounded: |r| <= |l| < |r| + u, under the claim's unit factors?"""
    u = c['u']
    for f in arith._unit_factors(c['lu'], c['ru'], c['implied']):
        for a in c['lv']:
            for b in c['rv']:
                x = a * f
                if b != 0.0 and (x > 0) == (b > 0) and abs(b) - 1e-9 * u <= abs(x) < abs(b) + u * (1 - 1e-9):
                    return True
    return False


def classify(c: dict, upstream: list[tuple[float, float]]) -> str:
    """carried_precision, or one of the three kinds of other: truncated, no_upstream, not_reproduced; for one
    flagged claim and the values before it."""
    sides = []
    for seg, vals, f, l, is_left in ((c['ls'], c['lv'], c['lf'], c['ll'], True), (c['rs'], c['rv'], c['rf'], c['rl'], False)):
        lits = literals(seg)
        # the right side is replaced only where it is an expression: one number with its unit (`0.002752 m**2`)
        # is the printed result itself
        if not is_left and len(arith.NUM_LITERAL.findall(arith.split_unit(seg)[0].replace(',', ''))) < 2:
            sides.append((seg, vals, f, l, lits, {}))
            continue
        opts = {}
        for j, (_a, _b, _lit, x, h) in enumerate(lits):
            cand = candidates(x, h, upstream)
            if cand:
                latest = max(cand, key=lambda t: t[0])[1]
                finest = min(cand, key=lambda t: (t[2], -t[0]))[1]
                opts[j] = sorted({latest, finest})
        sides.append((seg, vals, f, l, lits, opts))
    if not any(s[5] for s in sides):
        return 'truncated' if truncated(c) else 'no_upstream'

    def versions(seg, vals, f, l, lits, opts):
        yield vals
        keys = sorted(opts)
        subsets = ([s for r in range(1, len(keys) + 1) for s in itertools.combinations(keys, r)]
                   if len(keys) <= MAX_SUBSET else [tuple(keys)])
        for sub in subsets:
            for choice in itertools.product(*(opts[j] for j in sub)):
                v = substituted(seg, lits, dict(zip(sub, choice)), f, l, len(vals))
                if v is not None:
                    yield v

    rights = list(versions(*sides[1])) + [r for r in c['rights'] if r is not c['rv']]
    for left in versions(*sides[0]):
        for right in rights:
            if left is c['lv'] and right is c['rv']:
                continue
            if arith.agree_any_digit(left, right, c['lu'], c['ru'], c['u'], implied=c['implied']):
                return 'carried_precision'
    return 'truncated' if truncated(c) else 'not_reproduced'


def step_offsets(text: str, steps: list[str]) -> list[int | None]:
    """Where each step starts in the response, found in order; None when a step is not a substring of it."""
    out, pos = [], 0
    for s in steps:
        k = text.find(s, pos)
        out.append(k if k >= 0 else None)
        if k >= 0:
            pos = k + len(s)
    return out


def classify_response(text: str, stored_steps: list[dict] | None, upto: int | None = None) -> dict:
    """Walk the response's steps up to its last flagged one; classify each flagged claim against the values
    before it. With `stored_steps`, each walked step must give the stored claim and flag counts."""
    steps = e2_prm.steps_of(text)
    if stored_steps is not None and len(steps) != len(stored_steps):
        return {'unverified_response': True, 'flags': [], 'steps_walked': 0, 'unverified_steps': 0}
    flagged = [] if stored_steps is None else [k for k, s in enumerate(stored_steps) if s['digit_flags']]
    last = upto if upto is not None else (max(flagged) if flagged else -1)
    offs = step_offsets(text, steps)
    upstream_computed: list[tuple[float, float]] = []
    out, walked, unverified = [], 0, 0
    for k in range(last + 1):
        claims = walk(steps[k])
        walked += 1
        ok = stored_steps is None or (len(claims) == stored_steps[k]['claims']
                                      and sum(not c['ok_digit'] for c in claims) == stored_steps[k]['digit_flags'])
        if not ok:
            unverified += 1
        before = text[:offs[k]] if offs[k] is not None else '\n'.join(steps[:k])
        base = [n for ln in before.splitlines() for n in written_numbers(arith.normalise(ln))]
        lines = steps[k].splitlines()
        above = {i: [n for ln in lines[:i] for n in written_numbers(arith.normalise(ln))] for i in {c['line'] for c in claims}}
        for c in claims:
            if not c['ok_digit'] and ok:
                prefix = written_numbers(c['norm'][:c['at']])
                kind = classify(c, base + above[c['line']] + prefix + upstream_computed)
                out.append({'step': k, 'line': c['line'], 'left': c['ls'][:80], 'right': c['rs'][:80], 'class': kind})
            upstream_computed += computed_values(c)
    return {'unverified_response': False, 'flags': out, 'steps_walked': walked, 'unverified_steps': unverified}


# ------------------------------------------------------------------ the store

def _job(args):
    item_id, template_id, text, steps = args
    return item_id, template_id, classify_response(text, steps)


def classify_store(store: str, key: str, workers: int) -> dict:
    rows = load_rows(store, key)
    texts = score.texts_matching(store, key, rows)
    jobs = [(r['item_id'], r['template_id'], texts[r['item_id']], r['steps']) for r in rows
            if not r['unusable'] and r['answer']['label'] == 'correct' and any(s['digit_flags'] for s in r['steps'])]
    stored_flags = sum(sum(s['digit_flags'] for s in j[3]) for j in jobs)
    res = []
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        res = list(pool.map(_job, jobs, chunksize=8))
    return {'rows': rows, 'results': res, 'stored_flags': stored_flags}


def template_share(per_t: dict[str, list[int]], draws: int, seed: int) -> list[float] | None:
    """Template bootstrap of the carried share over the flags: per template (carried, flags), the ratio of sums."""
    if not per_t:
        return None
    a = np.array([v[0] for v in per_t.values()], float)
    b = np.array([v[1] for v in per_t.values()], float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(a), size=(draws, len(a)))
    d = a[idx].sum(axis=1) / b[idx].sum(axis=1)
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def summarise(cs: dict, draws: int, seed: int) -> dict:
    c = collections.Counter()
    per_t = collections.defaultdict(lambda: [0, 0])
    unv_resp = unv_steps = walked = 0
    for _item, tid, r in cs['results']:
        unv_resp += r['unverified_response']
        unv_steps += r['unverified_steps']
        walked += r['steps_walked']
        for f in r['flags']:
            c[f['class']] += 1
            per_t[tid][1] += 1
            per_t[tid][0] += f['class'] == 'carried_precision'
    flags = sum(c.values())
    correct = sum(not r['unusable'] and r['answer']['label'] == 'correct' for r in cs['rows'])
    return {'flags': flags, 'carried_precision': c['carried_precision'],
            'other': c['truncated'] + c['no_upstream'] + c['not_reproduced'], 'other_truncated': c['truncated'],
            'other_no_upstream': c['no_upstream'], 'other_not_reproduced': c['not_reproduced'],
            'share_carried': (c['carried_precision'] / flags) if flags else None,
            'share_ci': template_share(per_t, draws, seed),
            'stored_flags': cs['stored_flags'], 'correct_answers': correct,
            'responses_flagged': len(cs['results']), 'steps_walked': walked,
            'unverified_steps': unv_steps, 'unverified_responses': unv_resp}


# ------------------------------------------------------------------ the expert's confirmed slips

def expert_slips(workers: int) -> dict | None:
    """The round-3 sample's confirmed slips, re-found in the main store's traces and classified."""
    jl, sc = REVIEW / 'sample.jsonl', REVIEW / 'sample.csv'
    if not jl.exists() or not sc.exists():
        return None
    with open(sc, encoding='utf-8-sig', newline='') as fh:
        verdict = {r['code']: r['verdict'] for r in csv.DictReader(fh)}
    sample = [json.loads(l) for l in jl.read_text(encoding='utf-8').splitlines()]
    slips = [s for s in sample if verdict.get(s['code']) == 'slip']
    by_model = collections.defaultdict(list)
    for s in slips:
        by_model[s['model']].append(s)
    jobs = []
    for key, ss in by_model.items():
        rows = {r['item_id']: r for r in load_rows('main', key)}
        texts = score.texts_matching('main', key, list(rows.values()))
        for s in ss:
            jobs.append((s, texts.get(s['item_id'], '')))
    c = collections.Counter()
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for kind in pool.map(_slip_job, jobs, chunksize=4):
            c[kind] += 1
    return {'slips': len(slips), 'carried_precision': c['carried_precision'],
            'other': c['truncated'] + c['no_upstream'] + c['not_reproduced'], 'other_truncated': c['truncated'],
            'other_no_upstream': c['no_upstream'], 'other_not_reproduced': c['not_reproduced'],
            'step_changed': c['step_changed'],
            'claim_not_found': c['claim_not_found'], 'no_longer_flagged': c['no_longer_flagged']}


def _slip_job(args) -> str:
    s, text = args
    steps = e2_prm.steps_of(text) if text else []
    k = s['step']
    if k >= len(steps) or steps[k] != s['step_text']:
        hit = [j for j, st in enumerate(steps) if st == s['step_text']]
        if not hit:
            return 'step_changed'
        k = hit[0]
    claims = walk(steps[k])
    match = [c for c in claims if c['line'] == s['line'] and c['ls'][:80] == s['left'] and c['rs'][:80] == s['right']]
    if not match:
        return 'claim_not_found'
    if match[0]['ok_digit']:
        return 'no_longer_flagged'
    r = classify_response(text, None, upto=k)
    kinds = [f['class'] for f in r['flags'] if f['step'] == k and f['line'] == s['line']
             and f['left'] == s['left'] and f['right'] == s['right']]
    return kinds[0] if kinds else 'claim_not_found'


# ------------------------------------------------------------------ the run

def run(quick: bool, workers: int) -> dict:
    draws = 1_000 if quick else 10_000
    if quick:
        print('QUICK: 1,000 resamples')
    plan = stores()
    out = {'quick': quick, 'generated_by': 'python -m full_run_28092026.flag_precision' + (' --quick' if quick else ''),
           'written_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
           'resamples': draws, 'stores': {}}
    seed = 7100
    for s, ks in plan.items():
        models = {}
        for k in ks:
            seed += 10
            cs = classify_store(s, k, workers)
            models[k] = summarise(cs, draws, seed)
            m = models[k]
            print(f"{s}/{k}: {m['flags']} flags on correct answers ({m['stored_flags']} stored), "
                  f"{m['carried_precision']} carried precision; {m['unverified_steps']} steps unverified", flush=True)
        tot = collections.Counter()
        for m in models.values():
            tot.update({x: m[x] for x in ('flags', 'carried_precision', 'other', 'other_truncated', 'other_no_upstream',
                                          'other_not_reproduced', 'stored_flags', 'unverified_steps')})
        out['stores'][s] = {'models': models, 'all': dict(tot, share_carried=(tot['carried_precision'] / tot['flags'])
                                                          if tot['flags'] else None)}
    out['models'] = out['stores']['main']['models']
    out['all'] = out['stores']['main']['all']
    out['expert_confirmed'] = expert_slips(workers)
    if out['expert_confirmed']:
        e = out['expert_confirmed']
        print(f"expert-confirmed slips: {e['slips']}, carried precision {e['carried_precision']}, other {e['other']}, "
              f"not classified {e['step_changed'] + e['claim_not_found'] + e['no_longer_flagged']}")
    return out


# ------------------------------------------------------------------ writing

def render(res: dict) -> str:
    q = ' QUICK (1,000 resamples).' if res['quick'] else ''
    L = ['## Arithmetic flags: carried precision', '',
         f"Generated by `{res['generated_by']}`.{q} Every digit-rule flag on a response whose final answer is correct. "
         f"A flag is carried precision when the printed result is a correct rounding of the left side recomputed from "
         f"unrounded upstream values that appear earlier in the response (a number written with more digits, or an "
         f"earlier expression's value). Other flags are split three ways: the printed result is the left side cut "
         f"toward zero at its displayed precision instead of rounded (`truncated`); no operand has an earlier, more "
         f"precise value (`no upstream`); one has, and the result still misses the printed digits (`not reproduced`). "
         f"Share interval: 95% template bootstrap, {res['resamples']:,} resamples. Every step read reproduces the "
         f"stored claim and flag counts unless counted as unverified.", '']
    for s, blk in res['stores'].items():
        L += [f"### {'Main run (provider defaults)' if s == 'main' else s}", '',
              '| model | correct answers | flags | carried precision | share (95% CI) | other: truncated | '
              'other: no upstream | other: not reproduced | steps unverified |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for k, m in blk['models'].items():
            sh = '-' if m['share_carried'] is None else \
                f"{m['share_carried']:.3f} ({m['share_ci'][0]:.3f} to {m['share_ci'][1]:.3f})"
            L.append(f"| `{k}` | {m['correct_answers']} | {m['flags']} | {m['carried_precision']} | {sh} | "
                     f"{m['other_truncated']} | {m['other_no_upstream']} | {m['other_not_reproduced']} | "
                     f"{m['unverified_steps']} |")
        a = blk['all']
        L += [f"| all | | {a['flags']} | {a['carried_precision']} | "
              f"{'-' if a['share_carried'] is None else format(a['share_carried'], '.3f')} | {a['other_truncated']} | "
              f"{a['other_no_upstream']} | {a['other_not_reproduced']} | {a['unverified_steps']} |", '']
    e = res.get('expert_confirmed')
    if e:
        L += ['### The slips the domain expert confirmed (FLAG_REVIEW_3.md)', '',
              f"Of the {e['slips']} flags the expert confirmed as slips, re-found in the current traces and classified by "
              f"the same rule: {e['carried_precision']} carried precision, {e['other']} other ({e['other_truncated']} "
              f"truncated, {e['other_no_upstream']} with no upstream value, {e['other_not_reproduced']} not "
              f"reproduced); {e['no_longer_flagged']} are no longer flagged by the current rule, {e['step_changed']} "
              f"sit in a step whose text has changed and {e['claim_not_found']} could not be re-found.", '']
    else:
        L += ['The expert-review sample (scores/flag_review/round3/) is not on this machine; the confirmed slips are '
              'not classified.', '']
    return '\n'.join(L)


def write(res: dict) -> None:
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    OUT_MD.write_text(render(res), encoding='utf-8')
    print(f'wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}')


# ------------------------------------------------------------------ self-test

def selftest() -> int:
    """A carried-precision chain, a real slip, a slip with an upstream value that does not explain it, and the
    walk equal to arith.check."""
    cases = [
        # x is shown to four decimals, cut to 1.23 for the next line, and the product is right for 1.2389
        ('**Step 1:** x = 1.2389\n**Step 2:** y = 3 × 1.23 = 3.7167', ['carried_precision']),
        # a whole number rounded for display: 1235 carried from 1234.6
        ('**Step 1:** n = 1234.6\n**Step 2:** N = 1235 × 2 = 2469.2', ['carried_precision']),
        # the unrounded value is an earlier expression: 2/3 cut to 0.66 (itself flagged, a truncation), and the
        # product right for 2/3
        ('**Step 1:** a = 2/3 = 0.66\n**Step 2:** b = 0.66 × 30 = 20.00', ['truncated', 'carried_precision']),
        # a wrong digit with no earlier value behind it
        ('**Step 1:** A = π × 4.1616 / 4 = 3.267', ['no_upstream']),
        # ...and the printed result is never itself replaced: an earlier 3.26805 lies within a unit of 3.267 and
        # would make it agree with the left side, 3.26851
        ('**Step 1:** B = 3.26805\n**Step 2:** A = π × 4.1616 / 4 = 3.267 m**2', ['no_upstream']),
        # an earlier, more precise value exists, and the printed result is still not its rounding
        ('**Step 1:** x = 1.2389\n**Step 2:** y = 3 × 1.23 = 3.8000', ['not_reproduced']),
        # the result cut instead of rounded: 0.869565 x 0.259404 = 0.225568
        ('**Step 1:** p = 0.869565 × 0.259404 = 0.2255', ['truncated']),
    ]
    bad = 0
    for text, want in cases:
        steps = e2_prm.steps_of(text)
        stored = [{'claims': arith.check(s).checked, 'digit_flags': arith.check(s).checked - arith.check(s).consistent_digit}
                  for s in steps]
        for s, st in zip(steps, stored):
            ref = arith.check(s)
            mine = walk(s)
            if [c.ok_digit for c in ref.claims] != [c['ok_digit'] for c in mine] or \
                    [(c.left, c.right) for c in ref.claims] != [(c['ls'][:80], c['rs'][:80]) for c in mine]:
                bad += 1
                print(f'FAIL walk differs from arith.check on {s!r}')
        r = classify_response(text, stored)
        got = [f['class'] for f in r['flags']]
        if got != want or r['unverified_steps']:
            bad += 1
            print(f'FAIL {want}: got {got} (unverified {r["unverified_steps"]}) | {text!r}')
    # a stored count that disagrees with the walk is unverified and its flags left out
    text = cases[0][0]
    wrong = [{'claims': 9, 'digit_flags': 1} for _ in e2_prm.steps_of(text)]
    r = classify_response(text, wrong)
    if r['flags'] or not r['unverified_steps']:
        bad += 1
        print('FAIL a mismatched step is not left out')
    # the expert's sample: a claim re-found by step, line and sides is classified; a changed step is counted apart
    steps = e2_prm.steps_of(text)
    c = next(c for c in walk(steps[1]) if not c['ok_digit'])
    s = {'step': 1, 'step_text': steps[1], 'line': c['line'], 'left': c['ls'][:80], 'right': c['rs'][:80]}
    got = [_slip_job((s, text)), _slip_job((dict(s, step_text='** another step'), text)),
           _slip_job((dict(s, right='3.7168'), text))]
    if got != ['carried_precision', 'step_changed', 'claim_not_found']:
        bad += 1
        print(f'FAIL the expert sample path: {got}')
    total = len(cases) * 2 + 2
    print(f'{total - bad}/{total} checks pass')
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--quick', action='store_true', help='1,000 resamples (also ENGTRACE_QUICK=1)')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    write(run(quick_mode(a.quick), a.workers))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
