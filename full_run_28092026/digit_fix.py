"""D-156: the digit rule's parser fix, measured before against after: on gold, on the pilot's 300
expert-labelled traces, on the 220 flags a domain expert read (D-154), and on every full-run trace.

    python -m full_run_28092026.digit_fix           # FREE: D-156, arith.py at 3a7f247 against bc7dfa4: DIGIT_FIX.md
    python -m full_run_28092026.digit_fix --fix 2   # FREE: D-159, bc7dfa4 against arith.py as it stands: DIGIT_FIX_2.md

A second fix (D-159) was built from round 2's notes the same way; it is measured against the first,
on round 2's flags and on round 1's, whose slips it must go on flagging. FIXES lists both.

THE FIX. `arith.py`, E4's digit rule, read one chain of claims where a trace wrote several, and read
some numbers, units and variables as something else: a domain expert found 109 of 220 sampled
full-run flags to be the checker's (D-154). The rules it now applies are in its docstring, under
READING THE FULL RUN'S FLAGS. Nothing outside arith.py changed, and neither the answer check nor E3
imports it, so no answer verdict, E3 result or E5 prompt can move.

THE COMPARISON. The code before the fix is arith.py at PRE_FIX, loaded from git; the code after it is
arith.py as it stands, its SHA-256 printed, so the two sides differ in that file alone:
  gold          every gold solution: claims checked, claims flagged (must stay 0)
  the pilot     per step against the experts' step labels, with validate_scorer's inputs and skip
                rule: tp/fp/fn over all traces and on the correct-answer traces, and every step whose
                flag changes, toward the experts' label or away from it (none may move away)
  the reading   the 220 flags D-154's expert read: how many of the 109 checker verdicts are no longer
                flagged on their step, and how many of the 111 slips still are. The fix was built from
                these claims, so this shows it does what it was built to do; it is not the fixed
                rule's precision, which needs a fresh sample, read after the fix
  the full run  per model, over its answered traces: claims checked, claims flagged, traces with a
                flag, and steps whose flag appears or disappears
Adopted only if gold stays clean and no pilot step moves away from the experts, save the owner's
documented exception (EXCEPTIONS, D-156): five DeepSeek-R1 steps that gain a flag the experts do not
share. The old parser could not read their lines. Each is a real wrong digit at the precision the
trace displays (`3.867 × 10⁻⁵` for 3.8663 × 10⁻⁵; a 12-digit product wrong from its 7th digit), and the
experts' step labels judge a step's engineering, not its 7th digit. Any other step that moves away
fails the audit again. The full run's changed steps are listed, with the claims each side flags, in
scores/digit_fix_changes.jsonl, and the pilot's in scores/digit_fix_pilot_changes.jsonl, both
gitignored with the traces, to be read.

`pre_fix_digit_rule()` swaps the pre-fix arith.py into score.py for as long as it is open, so that a
validation of a published figure (validate_scorer.py, router.py --validate) reproduces it with the code
it was published with.
"""
from __future__ import annotations

import collections
import concurrent.futures as cf
import contextlib
import csv
import hashlib
import json
import subprocess
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVAL = REPO / 'evaluator_pilot_17092026' / 'evaluators'
for p in (str(REPO), str(EVAL), str(REPO / 'evaluator_pilot_17092026' / 'annotation'),
          str(REPO / 'evaluator_pilot_17092026' / 'analysis'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import arith  # noqa: E402
import e2_prm  # noqa: E402

from full_run_28092026 import score  # noqa: E402

PRE_FIX = '3a7f247954f487c360dd7b3a716bae794d9649f4'   # the last commit of arith.py before D-156
D156 = 'bc7dfa459e04a29fda80f5df431e8d5458617596'      # the D-156 fix
D159 = 'e78396259b69f551433ec447e19964721734205f'      # the D-159 fix
# The owner's documented exception (D-156): pilot steps, by trace code and step index, that gain a flag
# the experts do not share, each a real wrong digit at the displayed precision on a line the old parser
# could not read.
EXCEPTIONS = {('T-adac89', 5), ('T-3b9847', 12), ('T-d68431', 8), ('T-7e995f', 11), ('T-30adfd', 6)}
# Each fix, measured between two versions of arith.py (None: as it stands), on the reading it was built
# from and on the earlier readings, whose slips it must go on flagging (D-156, D-159).
FIXES = {1: {'label': 'D-156', 'pre': PRE_FIX, 'post': D156, 'design': 'round1', 'regression': [],
             'exceptions': EXCEPTIONS, 'report': 'DIGIT_FIX.md', 'changes': 'digit_fix_changes.jsonl',
             'pilot_changes': 'digit_fix_pilot_changes.jsonl'},
         2: {'label': 'D-159', 'pre': D156, 'post': D159, 'design': 'round2', 'regression': ['round1'],
             'exceptions': set(), 'report': 'DIGIT_FIX_2.md', 'changes': 'digit_fix_2_changes.jsonl',
             'pilot_changes': 'digit_fix_2_pilot_changes.jsonl'}}
FLAGS = score.SCORES / 'flag_review'
ROSTER_ALL = ['gpt-oss-20b', 'gemma-4-26b-a4b', 'deepseek-v4.1-flash', 'qwen3-235b-a22b-2507', 'glm-5.3-flash',
              'glm-5.3', 'muse-glimmer-30b', 'kimi-k3', 'gpt-5.4-mini', 'gemini-3.1-flash-lite',
              'claude-sonnet-5', 'qwen3.8-27b']


def source_at(commit: str) -> str:
    return subprocess.run(['git', 'show', f'{commit}:evaluator_pilot_17092026/evaluators/arith.py'],
                          cwd=REPO, capture_output=True, check=True).stdout.decode('utf-8')


def load(src: str, tag: str):
    """arith.py from its source, as its own module (registered, as dataclasses need)."""
    name = f'arith_{tag}'
    mod = types.ModuleType(name)
    mod.__file__ = str(EVAL / 'arith.py')
    sys.modules[name] = mod
    exec(compile(src, f'{name}.py', 'exec'), mod.__dict__)
    return mod


def pre_fix_arith():
    """arith.py as the pilot's digit-rule figures and the full run's first scores were computed."""
    return load(source_at(PRE_FIX), 'pre_d156')


@contextlib.contextmanager
def pre_fix_digit_rule():
    """score.py with the pre-fix arith.py, for as long as the block runs."""
    saved = score.arith
    score.arith = pre_fix_arith()
    try:
        yield score.arith
    finally:
        score.arith = saved


def flagged(rep) -> list:
    return [c for c in rep.claims if not c.ok_digit]


# ------------------------------------------------------------------ workers (the full run, gold)

_PRE = _POST = None


def _init(pre_src: str, post_src: str | None):
    global _PRE, _POST
    _PRE = load(pre_src, 'pre_worker')
    _POST = load(post_src, 'post_worker') if post_src else arith


def _trace(job):
    """Per step of one trace: claims and flags before and after; the steps whose flags differ."""
    model, item_id, text = job
    counts, changes = [], []
    for j, step in enumerate(e2_prm.steps_of(text)):
        a, b = _PRE.check(step), _POST.check(step)
        fa, fb = flagged(a), flagged(b)
        counts.append((len(a.claims), len(fa), len(b.claims), len(fb)))
        if bool(fa) != bool(fb):
            changes.append({'model': model, 'item_id': item_id, 'step': j,
                            'before': [(c.left, c.right) for c in fa], 'after': [(c.left, c.right) for c in fb]})
    return model, item_id, counts, changes


def _gold(job):
    item_id, solution = job
    a, b = _PRE.check(solution), _POST.check(solution)
    return item_id, len(a.claims), len(flagged(a)), len(b.claims), len(flagged(b))


# ------------------------------------------------------------------ the four comparisons

def pilot(pre, post=arith) -> dict:
    """validate_scorer's pilot inputs and skip rule, per step, before and after."""
    from full_run_28092026.validate_scorer import pilot_inputs
    truth, keyfile, items, texts = pilot_inputs()
    cnt = {(side, scope): collections.Counter() for side in ('before', 'after') for scope in ('all', 'hard')}
    moved = collections.Counter()
    away, changes = set(), []
    traces = skipped = 0
    for code, t in truth.items():
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        steps = e2_prm.steps_of(texts[key])
        if len(steps) != len(t['steps']):
            skipped += 1
            continue
        traces += 1
        for j, (st, lab) in enumerate(zip(steps, t['steps'])):
            c0, c1 = flagged(pre.check(st)), flagged(post.check(st))
            f0, f1 = bool(c0), bool(c1)
            if lab['label'] == 'incorrect':
                wrong = True
            elif lab['label'] in ('correct', 'alternative_correct'):
                wrong = False
            else:
                continue
            for scope in ('all', 'hard'):
                if scope == 'hard' and t['final_answer'] != 'correct':
                    continue
                for side, f in (('before', f0), ('after', f1)):
                    cnt[(side, scope)][('tp' if f else 'fn') if wrong else ('fp' if f else 'tn')] += 1
            if f0 != f1:
                moved['toward' if f1 == wrong else 'away'] += 1
                moved[('added' if f1 else 'removed', 'incorrect' if wrong else 'correct')] += 1
                if f1 != wrong:
                    away.add((code, j))
                changes.append({'code': code, 'model': key[0], 'item_id': key[1], 'step': j, 'label': lab['label'],
                                'before': [(c.left, c.right) for c in c0], 'after': [(c.left, c.right) for c in c1]})
    return {'traces': traces, 'skipped': skipped, 'counts': cnt, 'moved': moved, 'away': away, 'changes': changes}


def reading(pre, post=arith, rnd: str = 'round1') -> dict:
    """A round's flags, as the domain expert read them, before and after, on their own step."""
    folder = FLAGS / rnd
    verdict = {r['code']: r['verdict'] for r in csv.DictReader(open(folder / 'sample.csv', encoding='utf-8-sig'))}
    close = lambda a, b: abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))    # noqa: E731
    out = collections.Counter()
    for line in (folder / 'sample.jsonl').read_text(encoding='utf-8').splitlines():
        c = json.loads(line)
        same = lambda fl: any(any(close(x, y) for x in f.left_value for y in c['left_value'])     # noqa: E731
                              and any(close(x, y) for x in f.right_value for y in c['right_value']) for f in fl)
        before, after = flagged(pre.check(c['step_text'])), flagged(post.check(c['step_text']))
        out[('reproduced before', verdict[c['code']])] += same(before)
        out[(verdict[c['code']], 'the same claim' if same(after) else ('another claim' if after else 'none'))] += 1
    return out


def full_run(pre_src: str, post_src: str | None = None, workers: int = 8) -> tuple[dict, list]:
    jobs = []
    for key in ROSTER_ALL:
        p = score.SCORES / 'main' / f'{key}.jsonl'
        if not p.exists():
            continue
        rows = {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())}
        texts = score.texts_matching('main', key, rows.values())
        jobs += [(key, i, texts[i]) for i, r in rows.items() if r['status'] == 'answered']
    per = collections.defaultdict(collections.Counter)
    changes = []
    with cf.ProcessPoolExecutor(workers, initializer=_init, initargs=(pre_src, post_src)) as pool:
        for model, item_id, counts, ch in pool.map(_trace, jobs, chunksize=40):
            m = per[model]
            m['traces'] += 1
            m['claims before'] += sum(c[0] for c in counts)
            m['flagged before'] += sum(c[1] for c in counts)
            m['claims after'] += sum(c[2] for c in counts)
            m['flagged after'] += sum(c[3] for c in counts)
            m['traces flagged before'] += any(c[1] for c in counts)
            m['traces flagged after'] += any(c[3] for c in counts)
            m['steps flag removed'] += sum(bool(c[1]) and not c[3] for c in counts)
            m['steps flag added'] += sum(bool(c[3]) and not c[1] for c in counts)
            changes += ch
    return per, changes


def gold(pre_src: str, post_src: str | None = None, workers: int = 8) -> collections.Counter:
    items = score.pool_items()
    g = collections.Counter()
    with cf.ProcessPoolExecutor(workers, initializer=_init, initargs=(pre_src, post_src)) as pool:
        for _, ca, fa, cb, fb in pool.map(_gold, [(i, it['solution']) for i, it in items.items()], chunksize=40):
            g['items'] += 1
            g['claims before'] += ca
            g['flagged before'] += fa
            g['claims after'] += cb
            g['flagged after'] += fb
    return g


# ------------------------------------------------------------------ the report

def prf(c) -> str:
    p = c['tp'] / (c['tp'] + c['fp']) if c['tp'] + c['fp'] else float('nan')
    r = c['tp'] / (c['tp'] + c['fn']) if c['tp'] + c['fn'] else float('nan')
    return f"{c['tp']}/{c['fp']}/{c['fn']} | {p:.3f} | {r:.3f}"


REVIEW = FLAGS / 'agent_review'


def review_sample(per_kind: int = 40, seed: int = 7) -> int:
    """Steps whose flag the two fixes, taken together, removed or added (`arith.py` at PRE_FIX against
    `arith.py` as it stands), outside every expert round's steps and outside the set-aside model, for two
    independent review agents (D-160): per_kind of each kind, round-robin across the models. A removed flag
    should have been the checker's misreading, an added one a real slip. Written to
    scores/flag_review/agent_review/sample.jsonl, local: it holds trace text."""
    import random
    _, changes = full_run(source_at(PRE_FIX), None)
    read = set()
    for p in sorted(FLAGS.glob('round*/sample.jsonl')):
        read |= {(c['model'], c['item_id'], c['step']) for c in map(json.loads, p.read_text(encoding='utf-8').splitlines())}
    pool = {'removed': collections.defaultdict(list), 'added': collections.defaultdict(list)}
    for c in changes:
        if c['model'] != 'qwen3.8-27b' and (c['model'], c['item_id'], c['step']) not in read:
            pool['removed' if c['before'] else 'added'][c['model']].append(c)
    rng = random.Random(seed)
    chosen = []
    for kind in ('removed', 'added'):
        queues = [pool[kind][m] for m in sorted(pool[kind])]
        for q in queues:
            rng.shuffle(q)
        rng.shuffle(queues)
        picked = []
        while len(picked) < per_kind and any(queues):
            for q in queues:
                if q and len(picked) < per_kind:
                    picked.append(dict(q.pop(0), kind=kind))
        chosen += picked
    rng.shuffle(chosen)
    pre, texts, out = load(source_at(PRE_FIX), 'pre_review'), {}, []
    for c in chosen:
        key = c['model']
        if key not in texts:
            rows = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                                  .read_text(encoding='utf-8').splitlines())}
            texts[key] = score.texts_matching('main', key, rows.values())
        step = e2_prm.steps_of(texts[key][c['item_id']])[c['step']]
        lines = step.splitlines()

        def desc(rep):
            return [{'left': x.left, 'right': x.right, 'left_value': x.left_value, 'right_value': x.right_value,
                     'displayed_ulp': x.ulp, 'left_unit': x.left_unit, 'right_unit': x.right_unit,
                     'line': lines[x.line] if x.line < len(lines) else ''} for x in flagged(rep)]
        code = 'RV-' + hashlib.sha256(f"{key}|{c['item_id']}|{c['step']}".encode('utf-8')).hexdigest()[:8]
        out.append({'code': code, 'kind': c['kind'], 'model': key, 'item_id': c['item_id'], 'step': c['step'],
                    'flagged_before': desc(pre.check(step)), 'flagged_after': desc(arith.check(step)),
                    'step_text': step})
    REVIEW.mkdir(parents=True, exist_ok=True)
    with open(REVIEW / 'sample.jsonl', 'w', encoding='utf-8', newline='\n') as fh:
        for r in out:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
    sizes = {k: sum(len(v) for v in pool[k].values()) for k in pool}
    print(f"{len(out)} steps written to {REVIEW / 'sample.jsonl'} (seed {seed}): "
          + ', '.join(f"{k} {sum(r['kind'] == k for r in out)} of {sizes[k]}" for k in pool)
          + '; by model ' + str(dict(collections.Counter(r['model'] for r in out))))
    return 0


def reading_rows(rd) -> list:
    return ['| expert verdict | the same claim still flagged | another claim flagged in the step | no flag in the step |',
            '|---|---:|---:|---:|'] + [f"| {v} | {rd[(v, 'the same claim')]} | {rd[(v, 'another claim')]} | "
                                        f"{rd[(v, 'none')]} |" for v in ('checker', 'slip')]


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fix', type=int, default=1, choices=sorted(FIXES), help='1: D-156; 2: D-159')
    ap.add_argument('--review-sample', action='store_true', help='the two fixes\' changed steps, for the review agents')
    if ap.parse_args().review_sample:
        return review_sample()
    fx = FIXES[ap.parse_args().fix]
    pre_src = source_at(fx['pre'])
    post_src = source_at(fx['post']) if fx['post'] else None
    pre, post = load(pre_src, 'pre'), (load(post_src, 'post') if post_src else arith)
    if fx['post']:
        after = f"`arith.py` at `{fx['post'][:7]}`"
    else:
        sha = hashlib.sha256((EVAL / 'arith.py').read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        after = f'`arith.py` as it stands, SHA-256 (LF) `{sha[:16]}`'
    g = gold(pre_src, post_src)
    pil = pilot(pre, post)
    rd = reading(pre, post, fx['design'])
    regress = {r: reading(pre, post, r) for r in fx['regression']}
    per, changes = full_run(pre_src, post_src)
    CHANGES, PILOT_CHANGES = score.SCORES / fx['changes'], score.SCORES / fx['pilot_changes']
    for path, rows in ((CHANGES, changes), (PILOT_CHANGES, pil['changes'])):
        with open(path, 'w', encoding='utf-8', newline='\n') as fh:
            for c in rows:
                fh.write(json.dumps(c, ensure_ascii=False) + '\n')
    mv = pil['moved']
    unexplained = pil['away'] - fx['exceptions']
    ok = g['flagged after'] == 0 and not unexplained
    L = [f"# The digit rule's parser fix ({fx['label']}): before and after", '',
         f"Generated by `digit_fix.py --fix {ap.parse_args().fix}`; the comparison is defined in its docstring. "
         f"Before: `arith.py` at `{fx['pre'][:7]}`. After: {after}.", '',
         '## Gold, every solution', '', '| | before | after |', '|---|---:|---:|',
         f"| claims checked | {g['claims before']} | {g['claims after']} |",
         f"| claims flagged | {g['flagged before']} | {g['flagged after']} |", '',
         f"## The pilot's traces against the experts' step labels", '',
         f"{pil['traces']} traces scored per step, {pil['skipped']} skipped for a step count that differs from the "
         "labels' (validate_scorer.py's rule).", '',
         '| steps | before: tp/fp/fn, precision, recall | after: tp/fp/fn, precision, recall |', '|---|---|---|',
         f"| all traces | {prf(pil['counts'][('before', 'all')])} | {prf(pil['counts'][('after', 'all')])} |",
         f"| inside correct-answer traces | {prf(pil['counts'][('before', 'hard')])} | "
         f"{prf(pil['counts'][('after', 'hard')])} |", '',
         f"Steps whose flag changes: {mv['toward'] + mv['away']}; toward the experts' label {mv['toward']}, away "
         f"from it {mv['away']}. Flags removed from steps the experts call correct: "
         f"{mv[('removed', 'correct')]}; from steps they call incorrect: {mv[('removed', 'incorrect')]}. Flags "
         f"added to steps they call incorrect: {mv[('added', 'incorrect')]}; to steps they call correct: "
         f"{mv[('added', 'correct')]}.", '',
         f"## The {sum(rd[(v, k)] for v in ('checker', 'slip') for k in ('the same claim', 'another claim', 'none'))} "
         f"flags the domain expert read in {fx['design'].replace('round', 'round ')}, the ones this fix was built from", '',
         f"The code before the fix reproduces {rd[('reproduced before', 'slip')] + rd[('reproduced before', 'checker')]} "
         'of them on their own step. After the fix, on the same steps:', '']
    L += reading_rows(rd)
    L += ['', 'The fix was built from these notes, so this shows it does what it was built to do; the fixed '
          "rule's precision needs a fresh sample, read after the fix.", '']
    for r, rr in regress.items():
        L += [f"## {r.replace('round', 'Round ')}'s flags, read before: its slips must stay flagged", '',
              f"The code before this fix reproduces {rr[('reproduced before', 'slip')] + rr[('reproduced before', 'checker')]} "
              "of them as the flags that fix left; after it, on the same steps:", '']
        L += reading_rows(rr) + ['']
    L += ['## The full run, answered traces', '',
          '| model | traces | claims checked, before / after | claims flagged, before / after | '
          'traces with a flag, before / after | steps: flag removed / added |', '|---|---:|---:|---:|---:|---:|']
    tot = collections.Counter()
    for key in ROSTER_ALL:
        m = per.get(key)
        if not m:
            continue
        tot.update(m)
        L.append(f"| `{key}`{' (set aside)' if key == 'qwen3.8-27b' else ''} | {m['traces']} | "
                 f"{m['claims before']} / {m['claims after']} | {m['flagged before']} / {m['flagged after']} | "
                 f"{m['traces flagged before']} / {m['traces flagged after']} | "
                 f"{m['steps flag removed']} / {m['steps flag added']} |")
    L += [f"| all | {tot['traces']} | {tot['claims before']} / {tot['claims after']} | "
          f"{tot['flagged before']} / {tot['flagged after']} | {tot['traces flagged before']} / "
          f"{tot['traces flagged after']} | {tot['steps flag removed']} / {tot['steps flag added']} |", '',
          f'The changed steps, with the claims each side flags, are in `scores/{CHANGES.name}` and '
          f'`scores/{PILOT_CHANGES.name}` (local).', '',
          (f"Adopted ({fx['label']}): gold stays clean, and "
           + (f"the {len(pil['away'])} pilot steps that move away from the experts are the owner's documented "
              "exception: DeepSeek-R1 lines the old parser could not read, each a real wrong digit at the precision "
              'it displays, which the step labels, judging engineering, accept.' if pil['away'] else
              'no pilot step moves away from the experts.')
           if ok else 'NOT adopted: ' + ('gold has a flag. ' if g['flagged after'] else '') +
           (f'{len(unexplained)} pilot steps move away from the experts outside the documented exception: '
            + ', '.join(f'{c} step {j}' for c, j in sorted(unexplained)) + '.' if unexplained else ''))]
    (HERE / fx['report']).write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
