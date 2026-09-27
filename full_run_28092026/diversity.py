"""How varied are the items? Three checks, on the frozen pool and on each template's reach.

    python -m full_run_28092026.diversity      # writes DIVERSITY.md and diversity.json beside this file

Reads pool/ (local, gitignored) and draws 500 public seeds per template, the gate's seeds 0-499,
which are unrelated to the pool's. Writes counts, answer labels and item ids only: no question
text, no skeleton text, no seed.

Numbers are masked: a decimal, integer, thousands-separated or scientific literal, including
'a x 10^b', 'a x 10' with a superscript exponent and 'aEb', becomes '#', and a minus sign on a
literal goes with it, so a sign change is not a variant.

  question skeleton   the question with numbers masked. It counts wording and every non-numeric
                      choice, so a swapped fluid or material name is a variant.
  reasoning path      the gold solution's lines that contain '=', numbers masked, in order.
                      upper: as they stand, so a name inside an equation line is a variant.
                      lower: words that vary across the template's questions (names, choices) are
                      masked as well, so a name swap is not a variant. A structural choice named
                      in the question is masked too, so this reading can undercount.
  answer variant      the gold answer segment, from the first '**Answer:**' to the end, with
                      numbers and question-varying words masked. One variant means every instance
                      answers in the same form; a classification spreads over its labels.
  near-duplicate      two pool items of one template, neither a marked repeat, with the same
                      question skeleton and every number within a relative tolerance of its
                      counterpart, position by position: 1%, 5% and 10%.
"""
from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generate_testset import difficulty_map  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402
from full_run_28092026.freeze import POOL_DIR, answer_types  # noqa: E402

SEEDS = 500
TOLS = (0.01, 0.05, 0.10)
SUPMAP = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
NUM = re.compile(r'''
    (?P<neg>(?<![\w)\]])[-−])?
    (?:
      (?P<sci>(?P<m>\d+(?:\.\d+)?)\s*(?:×|x|\*|\\times)\s*10\s*(?:\^|\*\*)\s*[{(]?\s*(?P<e>[-−+]?\s*\d+)(?:\s*[})])?)
    | (?P<sup>(?P<m2>\d+(?:\.\d+)?)\s*(?:×|x|\*|\\times)\s*10(?P<e2>[⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))
    | (?P<enot>\d+(?:\.\d+)?[eE][-+]?\d+)
    | (?P<plain>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?|\.\d+)
    )''', re.X)
WORD = re.compile(r'[A-Za-z][A-Za-z\-]+')


def _value(m: re.Match) -> float:
    if m.group('sci'):
        v = float(m.group('m')) * 10 ** int(m.group('e').replace('−', '-').replace(' ', ''))
    elif m.group('sup'):
        v = float(m.group('m2')) * 10 ** int(m.group('e2').translate(SUPMAP))
    elif m.group('enot'):
        v = float(m.group('enot'))
    else:
        v = float(m.group('plain').replace(',', ''))
    return -v if m.group('neg') else v


def mask(text: str) -> tuple[str, list[float]]:
    vals: list[float] = []

    def rep(m):
        vals.append(_value(m))
        return '#'
    return NUM.sub(rep, text), vals


def answer_segment(solution: str) -> str:
    i = solution.find('**Answer:**')
    if i < 0:
        i = solution.find('Answer:')
    if i < 0:
        lines = [x for x in solution.splitlines() if x.strip()]
        return lines[-1] if lines else ''
    return ' '.join(solution[i:].split())


def varying_words(questions: list[str]) -> set[str]:
    sets = [set(WORD.findall(mask(q)[0])) for q in questions]
    return set().union(*sets) - set.intersection(*sets) if sets else set()


def word_masker(words: set[str]):
    if not words:
        return lambda s: s
    rx = re.compile(r'\b(?:' + '|'.join(re.escape(w) for w in sorted(words, key=len, reverse=True))
                    + r')\b')
    return lambda s: rx.sub('@', s)


def paths(solution: str, mw) -> tuple[str, str]:
    lines = [mask(x.strip())[0] for x in solution.splitlines() if '=' in x]
    upper = '\n'.join(lines)
    return upper, mw(upper)


def features(question: str, solution: str, mw) -> dict:
    skel, vals = mask(question)
    up, low = paths(solution, mw)
    return {'skel': skel, 'vals': vals, 'up': up, 'low': low,
            'ans': mw(mask(answer_segment(solution))[0])}


def close(a: list[float], b: list[float], tol: float) -> bool:
    return len(a) == len(b) and all(
        x == y or abs(x - y) <= tol * max(abs(x), abs(y)) for x, y in zip(a, b))


def load_pool() -> dict[str, list[dict]]:
    by_t: dict[str, list[dict]] = collections.defaultdict(list)
    for path in sorted(POOL_DIR.rglob('*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            r = json.loads(line)
            by_t['template_' + r['id']].append(r)
    return by_t


def analyse() -> list[dict]:
    pool, levels, types = load_pool(), difficulty_map(), answer_types()
    rows = []
    for ref in discover(None):
        tid = ref.template_id
        draws = [generate(ref, s, capture=False) for s in range(SEEDS)]
        draws = [d for d in draws if d.ok]
        mw = word_masker(varying_words([d.question for d in draws]))
        reach = [features(d.question, d.solution, mw) for d in draws]
        items = pool[tid]
        feats = [features(r['question'], r['solution'], mw) for r in items]
        fresh = [(r['item_id'], f) for r, f in zip(items, feats) if 'repeat_of' not in r]
        near = {t: [] for t in TOLS}
        for i in range(len(fresh)):
            for j in range(i + 1, len(fresh)):
                (ia, fa), (ib, fb) = fresh[i], fresh[j]
                if fa['skel'] == fb['skel']:
                    for t in TOLS:
                        if close(fa['vals'], fb['vals'], t):
                            near[t].append([ia, ib])
        spread = collections.Counter(f['ans'] for f in feats)
        rows.append({
            'template_id': tid, 'branch': ref.branch, 'level': levels.get(tid),
            'answer_type': types.get(tid, 'unknown'),
            'pool': {
                'items': len(items),
                'distinct_questions': len({r['question'] for r in items}),
                'question_skeletons': len({f['skel'] for f in feats}),
                'paths_upper': len({f['up'] for f in feats}),
                'paths_lower': len({f['low'] for f in feats}),
                'answer_variants': len(spread),
                'answer_spread': sorted(spread.values(), reverse=True),
                'answer_labels': ({k[:90]: v for k, v in spread.most_common()}
                                  if types.get(tid) == 'classification' else None),
                'near_duplicate_pairs': {f'{int(t * 100)}%': len(v) for t, v in near.items()},
                'near_duplicates_5pct': near[0.05],
            },
            'reach': {
                'draws': len(draws),
                'distinct_questions': len({d.question for d in draws}),
                'question_skeletons': len({f['skel'] for f in reach}),
                'paths_upper': len({f['up'] for f in reach}),
                'paths_lower': len({f['low'] for f in reach}),
                'answer_variants': len({f['ans'] for f in reach}),
            },
        })
    return rows


def bucket(n: int) -> str:
    return '1' if n == 1 else '2' if n == 2 else '3-5' if n <= 5 else '6-15' if n <= 15 else '16+'


def report(rows: list[dict]) -> str:
    order = ['1', '2', '3-5', '6-15', '16+']
    out = ['# Item diversity: the frozen pool and each template\'s reach', '',
           'Generated by `diversity.py`. Definitions are in its docstring. The pool is the 2,250 '
           'frozen items (D-114); reach is 500 public draws per template.', '',
           '## How many templates, by count of distinct variants', '',
           '| variants | question skeletons, pool | reasoning paths, pool, lower | upper | '
           'reasoning paths, 500 draws, lower | upper | answer variants, pool |',
           '|---|---:|---:|---:|---:|---:|---:|']
    cols = [('pool', 'question_skeletons'), ('pool', 'paths_lower'), ('pool', 'paths_upper'),
            ('reach', 'paths_lower'), ('reach', 'paths_upper'), ('pool', 'answer_variants')]
    for b in order:
        cells = [sum(bucket(r[s][k]) == b for r in rows) for s, k in cols]
        out.append(f'| {b} | ' + ' | '.join(map(str, cells)) + ' |')
    by_branch = collections.defaultdict(list)
    for r in rows:
        by_branch[r['branch']].append(r)
    out += ['', '## One reasoning path in every pool instance, by branch', '',
            '| branch | templates | one path, lower | one path, upper | one path in 500 draws, upper |',
            '|---|---:|---:|---:|---:|']
    for br, rs in sorted(by_branch.items()):
        out.append(f"| {br} | {len(rs)} | {sum(r['pool']['paths_lower'] == 1 for r in rs)} | "
                   f"{sum(r['pool']['paths_upper'] == 1 for r in rs)} | "
                   f"{sum(r['reach']['paths_upper'] == 1 for r in rs)} |")
    nd = {k: sum(r['pool']['near_duplicate_pairs'][k] for r in rows) for k in ('1%', '5%', '10%')}
    hit = [r for r in rows if r['pool']['near_duplicate_pairs']['5%']]
    out += ['', '## Near-duplicates in the pool', '',
            f"Pairs: {nd['1%']} at 1%, {nd['5%']} at 5%, {nd['10%']} at 10%. "
            f"Templates with a pair at 5%: {len(hit)}.", '']
    for r in hit:
        pairs = ', '.join(f'{a} ~ {b}' for a, b in r['pool']['near_duplicates_5pct'])
        out.append(f"- `{r['template_id']}`: {pairs}")
    out += ['', '## Classification templates: how the pool\'s answers spread', '']
    for r in rows:
        if r['answer_type'] == 'classification':
            labels = '; '.join(f'{v} x {k}' for k, v in r['pool']['answer_labels'].items())
            out.append(f"- `{r['template_id']}` ({r['pool']['answer_variants']} variants): {labels}")
    out += ['', '## Every template', '',
            '| template | branch | level | answer | pool: skeletons | paths lower-upper | answers | '
            'near-dup 5% | 500 draws: questions | skeletons | paths lower-upper |',
            '|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in sorted(rows, key=lambda r: (r['pool']['paths_upper'], r['pool']['question_skeletons'],
                                         r['template_id'])):
        p, q = r['pool'], r['reach']
        out.append(f"| `{r['template_id'][len('template_'):]}` | {r['branch'].split('_')[0]} | "
                   f"{r['level']} | {r['answer_type']} | {p['question_skeletons']} | "
                   f"{p['paths_lower']}-{p['paths_upper']} | {p['answer_variants']} | "
                   f"{p['near_duplicate_pairs']['5%']} | {q['distinct_questions']} | "
                   f"{q['question_skeletons']} | {q['paths_lower']}-{q['paths_upper']} |")
    return '\n'.join(out) + '\n'


def main() -> int:
    rows = analyse()
    (HERE / 'diversity.json').write_text(json.dumps(rows, indent=1, ensure_ascii=False) + '\n',
                                         encoding='utf-8', newline='\n')
    md = report(rows)
    (HERE / 'DIVERSITY.md').write_text(md, encoding='utf-8', newline='\n')
    print(md.split('## Every template')[0])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
