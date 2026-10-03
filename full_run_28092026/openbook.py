"""C4, the open-book condition (docs/EVALUATION_NEXT_STEPS.md C4; D-183): the governing equations supplied with the
question, bounding what retrieval could add, on the 450-item subsample.

    python -m full_run_28092026.openbook --survey          # FREE: where each template's equations can come from; writes OPENBOOK_SURVEY.md
    python -m full_run_28092026.openbook --build           # FREE: openbook/items.jsonl (local) and openbook/manifest.jsonl (hashes)
    python -m full_run_28092026.run_traces --variant openbook --dry-run   # then as any variant
    python -m full_run_28092026.openbook --survey --version 2             # FREE: the corrected filter; OPENBOOK_SURVEY_2.md
    python -m full_run_28092026.openbook --build --version 2              # FREE: openbook/items2.jsonl (local), openbook/manifest2.jsonl
    python -m full_run_28092026.openbook --diff                           # FREE: version 2 against the run's manifest; OPENBOOK_DIFF.md
    python -m full_run_28092026.openbook --carry                          # FREE: unchanged items' openbook traces copied into traces/openbook2/
    python -m full_run_28092026.run_traces --variant openbook2 --dry-run  # BILLS when run: the changed items only

VERSION 2 (2026-10-03, after the review of D-183). The first build's line filter read the fraction of an exponent as a
number (`**0.2857` -> 2857, so `rackett_equation_volume` lost its only equation), refused a side that is `0` (so
`vdw_solve_for_volume` lost its cubic), let through notes from the template code with an `=` in them (rounding and
reviewer remarks), and took a gold solution's symbolic lines up to the first substitution, which for a symbolic answer
is the answer itself. Version 2: a number is read whole and never out of an exponent or a fraction; a bare integer
side up to 10 is allowed; a line is a note, not an equation, if it carries a quotation mark, scientific notation, a
Markdown quote marker, or an English stopword or a code or review word outside a parenthetical remark (the remark
itself is dropped); the gold source is not used where the answer is symbolic or the line names the answer. The run of
D-183 used version 1 and its manifest stays as the record of that run; version 2 is the `openbook2` arm, and
`--carry` copies the version-1 trace rows of the items whose question is byte-identical under both versions, so only
the changed items are bought again.

WHERE THE EQUATIONS COME FROM. Two sources, in order of preference, recorded per template in the manifest:
  docstring   the template function's own docstring, where it states its equations in symbols (lines with an `=`
              and no digits apart from exponents and small integer constants): the handbook entry, written before
              any instance;
  gold        the gold solution's symbolic lines, the same rule applied to the solution text of the frozen item,
              taken before any line that substitutes numbers; used only where the docstring states none.
A template with neither source is left out of the arm and listed. The question is unchanged; the reference block is
appended after it, headed "Reference equations (from an engineering handbook):", one equation per line, duplicates
removed, in the order found. The prompt template is the main run's, so the prompt hash is unchanged; the modified
question's hash is the row's item_sha256 and the original's its original_sha256, as the paraphrase arm records them.
What the condition measures: with the method given, the trace still has to set the problem up, substitute and
compute, so the change against the main run on the same items is the part of the gap that formula recall accounts
for. It does not supply data tables or constants beyond what the equation lines carry.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from full_run_28092026 import score, subsamples  # noqa: E402
from tests.template_integrity.core import discover  # noqa: E402

OUT = HERE / 'openbook'
SURVEY = HERE / 'OPENBOOK_SURVEY.md'
HEADER = 'Reference equations (from an engineering handbook):'
EQ = re.compile(r'=')
DIGITS = re.compile(r'(?<![\^*\w])\d+(?:\.\d+)?(?![\w])')      # version 1: a number that is not an exponent or part of a name
DIGITS2 = re.compile(r'(?<![\^*\w.])\d+(?:\.\d+)?(?![\w.])')  # version 2: nor the fraction of one
SCI = re.compile(r'\d(?:\.\d+)?e[+-]?\d', re.I)
STEP = re.compile(r'^\s*\**Step\s*\d+\**:?\s*', re.I)
REMARK = re.compile(r'\([^()]*\)')
STOP = re.compile(r'\b(?:the|is|are|was|were|been|being|has|have|had|does|did|will|would|should|must|cannot|we|our|you|your|'
                  r'this|that|these|those|which|who|whom|because|since|then|unless|otherwise|also|however|but|and|not|for|of|to|from|'
                  r'with|than|when|where|while|after|before|on|at|'
                  r'note|notes|judge|judges|reported|report|reviewer|review|round|rounded|rounding|digit|digits|significant|'
                  r'tolerance|printed|print|check|checked|fix|fixed|bug|test|tests|gold|answer|answers|instance|instances|template|'
                  r'example|assume|assumed|assumes|verify|error|errors|comment|question|item|items|seed|random|return|returns|dict|'
                  r'string|format|formatted|float|json|python|function|value|values|decimal|places|precision|output|satisfies|'
                  r'integer|half-up|typical|typically|usually|hence|thus|therefore|i\.e|e\.g|etc)\b')     # lowercase prose only: `By`, `At` are symbols
SENTENCE = re.compile(r'^[A-Z][a-z]+ [a-z]+\b')        # "The system's ...", "One judge ...": a sentence, not an equation
OPERATOR = re.compile(r'[-+*/^]|\d|\(')
VERSION = 1
SYMBOLIC_ANSWERS = ('symbolic', 'expression', 'function', 'formula', 'text', 'string')


def version_paths(version: int):
    sfx = '' if version == 1 else str(version)
    return OUT / f'items{sfx}.jsonl', OUT / f'manifest{sfx}.jsonl', HERE / ('OPENBOOK_SURVEY.md' if version == 1 else f'OPENBOOK_SURVEY_{version}.md')


def symbolic_lines(text: str) -> list[str]:
    """Lines that state an equation in symbols: an `=`, letters on both sides, and no free-standing number other
    than small integer constants (0 to 10) and exponents. Version 2 (module VERSION) applies the corrected rules of
    the docstring: numbers read whole, `= 0` allowed, notes and prose rejected, parenthetical remarks dropped."""
    out = []
    for raw in text.splitlines():
        line = STEP.sub('', raw).strip().strip('-*• ').strip()
        if VERSION >= 2:
            if line.startswith('>') or '"' in line or '\u201c' in line or '==' in line or SCI.search(line):
                continue                                   # a quoted remark, code, or an instance number
            line = REMARK.sub(lambda m: '' if STOP.search(m.group(0)) else m.group(0), line).strip().rstrip('.;,').strip()
            line = re.sub(r'\s+', ' ', line)
            if '=' not in line or len(line) > 120 or STOP.search(line) or SENTENCE.match(line):
                continue
            left, _, right = line.partition('=')
            if not re.search(r'[A-Za-z]', left) or len(left.strip()) > 40:
                continue
            if not re.search(r'[A-Za-z]', right) and not re.fullmatch(r'\s*-?\d{1,2}\s*', right):
                continue
            if not OPERATOR.search(right) and len(right.strip()) > 3:
                continue                                   # "general = dense": a word, not an expression
            nums = [float(x) for x in DIGITS2.findall(line)]
            if any(n > 10 or n != int(n) for n in nums):
                continue
            if line not in out:
                out.append(line)
            continue
        if '=' not in line or len(line) > 160:
            continue
        left, _, right = line.partition('=')
        if not re.search(r'[A-Za-z]', left) or not re.search(r'[A-Za-z]', right):
            continue
        nums = [float(x) for x in DIGITS.findall(line)]
        if any(n > 10 or n != int(n) for n in nums):
            continue
        if re.search(r'\b(is|are|was|were|the|then|so|thus|where)\b', line.split('=')[0]) and len(line.split()) > 12:
            continue
        line = re.sub(r'\s+', ' ', line)
        if line not in out:
            out.append(line)
    return out


def tid_of(ref) -> str:
    return ref.template_id if ref.template_id.startswith('template_') else 'template_' + ref.template_id


def docstring_equations(ref) -> list[str]:
    try:
        fn = ref.load()
    except Exception:  # noqa: BLE001
        return []
    doc = inspect.getdoc(fn) or ''
    return symbolic_lines(doc)


def gold_equations(solution: str, answer_type: str | None = None) -> list[str]:
    """The gold's symbolic lines, stopping at the first line that substitutes numbers into an equation. Version 2:
    not used where the answer is symbolic (the symbolic lines would be the answer), and stopping at a line that
    names the answer."""
    if VERSION >= 2 and answer_type and any(w in str(answer_type).lower() for w in SYMBOLIC_ANSWERS):
        return []
    out = []
    for raw in solution.splitlines():
        line = STEP.sub('', raw).strip()
        if '=' in line and re.search(r'\d+\.\d+|\d{3,}', line):
            break
        if VERSION >= 2 and re.search(r'\banswer\b|\bfinal\b|\\boxed', line, re.I):
            break
        out += symbolic_lines(line)
    return [x for i, x in enumerate(out) if x not in out[:i]]


def sources() -> dict:
    items = score.pool_items()
    first = {}
    for i in subsamples.paraphrase_ids():
        first.setdefault(items[i]['template_id'], items[i])
    out = {}
    for ref in discover(None):
        tid = tid_of(ref)
        if tid not in first:
            continue
        d = docstring_equations(ref)
        g = gold_equations(first[tid]['solution'], first[tid].get('answer_type')) if not d else []
        out[tid] = {'source': 'docstring' if d else ('gold' if g else None), 'equations': d or g}
    return out


def survey() -> int:
    src = sources()
    c = {'docstring': 0, 'gold': 0, None: 0}
    for v in src.values():
        c[v['source']] += 1
    _items, _man, SURVEY = version_paths(VERSION)
    L = [f'# The open-book condition: where each template\'s equations come from (filter version {VERSION})', '',
         f'Generated by `openbook.py --survey --version {VERSION}`; the rule for a symbolic line is in its docstring. Counts and '
         'equation text only (equations are the templates\' public code or the gold\'s symbolic lines, with no instance numbers).', '',
         f"| source | templates |", '|---|---:|', f"| the template's docstring | {c['docstring']} |",
         f"| the gold solution's symbolic lines | {c['gold']} |", f"| neither (left out of the arm) | {c[None]} |", '',
         '| template | source | equations |', '|---|---|---|']
    for tid, v in sorted(src.items()):
        L.append(f"| `{tid.removeprefix('template_')}` | {v['source'] or '-'} | " + ('<br>'.join(f'`{e}`' for e in v['equations'][:6]) or '-') + ' |')
    SURVEY.write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L[:9]))
    print(f"... {len(src)} templates in the survey; the table is in {SURVEY.name}")
    return 0


def build() -> int:
    items = score.pool_items()
    src = sources()
    OUT.mkdir(exist_ok=True)
    items_path, man_path, _survey = version_paths(VERSION)
    rows, man = [], []
    for i in subsamples.paraphrase_ids():
        it = items[i]
        s = src.get(it['template_id'])
        if not s or not s['source']:
            continue
        q = it['question'].rstrip() + '\n\n' + HEADER + '\n' + '\n'.join(f'- {e}' for e in s['equations'])
        rows.append({'item_id': i, 'question': q, 'source': s['source']})
        man.append({'item_id': i, 'template_id': it['template_id'], 'source': s['source'],
                    'sha256': hashlib.sha256(q.encode('utf-8')).hexdigest(), 'original_sha256': it['sha256'],
                    'equations': len(s['equations'])})
    items_path.write_text(''.join(json.dumps(r, ensure_ascii=False) + '\n' for r in rows), encoding='utf-8')
    man_path.write_text(''.join(json.dumps({**r, 'filter_version': VERSION}) + '\n' for r in man), encoding='utf-8', newline='\n')
    print(f'{len(rows)} items written to {items_path} (local) with the manifest beside it; '
          f'{sum(1 for r in man if r["source"] == "docstring")} from docstrings, {sum(1 for r in man if r["source"] == "gold")} from the gold')
    return 0


def diff() -> int:
    """Version 2's manifest against the run's (version 1): templates in or out, items whose question changed, lines
    dropped and added per template. Equation text only; written to OPENBOOK_DIFF.md."""
    _i1, m1p, _s1 = version_paths(1)
    _i2, m2p, _s2 = version_paths(2)
    if not m1p.exists() or not m2p.exists():
        raise SystemExit('both manifests are needed: --build and --build --version 2 first')
    m1 = {r['item_id']: r for r in map(json.loads, m1p.read_text(encoding='utf-8').splitlines())}
    m2 = {r['item_id']: r for r in map(json.loads, m2p.read_text(encoding='utf-8').splitlines())}
    global VERSION
    VERSION, src1 = 1, None
    src1 = sources()
    VERSION = 2
    src2 = sources()
    t1, t2 = {r['template_id'] for r in m1.values()}, {r['template_id'] for r in m2.values()}
    same = [i for i in m1 if i in m2 and m1[i]['sha256'] == m2[i]['sha256']]
    changed = [i for i in m1 if i in m2 and m1[i]['sha256'] != m2[i]['sha256']]
    L = ['# The open-book condition: filter version 2 against the run\'s version 1', '',
         'Generated by `openbook.py --diff`. Version 1 is the manifest the `openbook` arm ran with (D-183); version 2 is the '
         'corrected filter (the module docstring). Equation text only.', '',
         f'| | version 1 | version 2 |', '|---|---:|---:|',
         f'| templates in the arm | {len(t1)} | {len(t2)} |', f'| items | {len(m1)} | {len(m2)} |',
         f'| from the docstring / the gold | {sum(r["source"] == "docstring" for r in m1.values())} / {sum(r["source"] == "gold" for r in m1.values())} | '
         f'{sum(r["source"] == "docstring" for r in m2.values())} / {sum(r["source"] == "gold" for r in m2.values())} |',
         f'| equation lines over the templates | {sum(len(v["equations"]) for v in src1.values() if v["source"])} | {sum(len(v["equations"]) for v in src2.values() if v["source"])} |',
         '', f'Items whose question is byte-identical under both versions: {len(same)} (their version-1 traces carry over with `--carry`); '
         f'items in both whose block changed: {len(changed)}; items only in version 1: {len(set(m1) - set(m2))}; only in version 2: {len(set(m2) - set(m1))}.', '',
         '## Templates that enter or leave the arm', '',
         'Enter under version 2: ' + (', '.join(f'`{t.removeprefix("template_")}`' for t in sorted(t2 - t1)) or 'none') + '.', '',
         'Leave under version 2: ' + (', '.join(f'`{t.removeprefix("template_")}`' for t in sorted(t1 - t2)) or 'none') + '.', '',
         '## Lines dropped and added, per template whose block changed', '', '| template | dropped under version 2 | added under version 2 |', '|---|---|---|']
    for tid in sorted(t1 | t2):
        a = src1.get(tid, {}).get('equations') or []
        b = src2.get(tid, {}).get('equations') or []
        if a == b:
            continue
        dropped = [x for x in a if x not in b]
        added = [x for x in b if x not in a]
        L.append(f"| `{tid.removeprefix('template_')}` | " + ('<br>'.join(f'`{x}`' for x in dropped) or '-') + ' | ' + ('<br>'.join(f'`{x}`' for x in added) or '-') + ' |')
    (HERE / 'OPENBOOK_DIFF.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L[:14]))
    print(f'... the per-template table is in OPENBOOK_DIFF.md ({sum(1 for l in L if l.startswith("| `"))} templates changed)')
    return 0


def carry(models=('claude-sonnet-5', 'gpt-5.4-mini', 'gpt-oss-20b')) -> int:
    """Copy the `openbook` arm's trace rows of the items whose question is byte-identical under version 2 into
    traces/openbook2/<model>.jsonl, marked carried_from: openbook; the row's request, text and bill are those of
    the original call, which is a valid sample for the same question. Items already present in openbook2 are
    left alone, so the copy can be repeated."""
    _i1, m1p, _s1 = version_paths(1)
    _i2, m2p, _s2 = version_paths(2)
    m1 = {r['item_id']: r for r in map(json.loads, m1p.read_text(encoding='utf-8').splitlines())}
    m2 = {r['item_id']: r for r in map(json.loads, m2p.read_text(encoding='utf-8').splitlines())}
    same = {i for i in m1 if i in m2 and m1[i]['sha256'] == m2[i]['sha256']}
    src, dst = HERE / 'traces' / 'openbook', HERE / 'traces' / 'openbook2'
    dst.mkdir(parents=True, exist_ok=True)
    for k in models:
        f = src / f'{k}.jsonl'
        if not f.exists():
            print(f'{k}: no openbook traces')
            continue
        have = set()
        out = dst / f'{k}.jsonl'
        if out.exists():
            have = {json.loads(l)['item_id'] for l in out.read_text(encoding='utf-8').splitlines() if l.strip()}
        n = 0
        with open(out, 'a', encoding='utf-8', newline='\n') as fh:
            for l in f.read_text(encoding='utf-8').splitlines():
                if not l.strip():
                    continue
                r = json.loads(l)
                if r['item_id'] in same and r['item_id'] not in have and r.get('status') in ('answered', 'empty'):
                    r = {**r, 'variant': 'openbook2', 'carried_from': 'openbook'}
                    fh.write(json.dumps(r, ensure_ascii=False) + '\n')
                    have.add(r['item_id'])
                    n += 1
        print(f'{k}: {n} rows carried into {out.name}; {len(m2) - len(have)} of {len(m2)} items still to run')
    return 0


def main() -> int:
    global VERSION
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', action='store_true')
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--diff', action='store_true', help='version 2 against the run manifest (version 1)')
    ap.add_argument('--carry', action='store_true', help='copy unchanged items\' openbook traces into the openbook2 arm')
    ap.add_argument('--version', type=int, default=1, choices=(1, 2), help='the line filter: 1 as run (D-183), 2 corrected')
    a = ap.parse_args()
    VERSION = a.version
    if a.diff:
        return diff()
    if a.carry:
        return carry()
    if a.build:
        return build()
    return survey()


if __name__ == '__main__':
    sys.exit(main())
