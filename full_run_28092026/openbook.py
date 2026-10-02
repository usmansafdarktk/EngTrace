"""C4, the open-book condition (docs/EVALUATION_NEXT_STEPS.md C4; D-183): the governing equations supplied with the
question, bounding what retrieval could add, on the 450-item subsample.

    python -m full_run_28092026.openbook --survey          # FREE: where each template's equations can come from; writes OPENBOOK_SURVEY.md
    python -m full_run_28092026.openbook --build           # FREE: openbook/items.jsonl (local) and openbook/manifest.jsonl (hashes)
    python -m full_run_28092026.run_traces --variant openbook --dry-run   # then as any variant

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
DIGITS = re.compile(r'(?<![\^*\w])\d+(?:\.\d+)?(?![\w])')      # a number that is not an exponent or part of a name
STEP = re.compile(r'^\s*\**Step\s*\d+\**:?\s*', re.I)


def symbolic_lines(text: str) -> list[str]:
    """Lines that state an equation in symbols: an `=`, letters on both sides, and no free-standing number other
    than small integer constants (0 to 10) and exponents."""
    out = []
    for raw in text.splitlines():
        line = STEP.sub('', raw).strip().strip('-*• ').strip()
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


def gold_equations(solution: str) -> list[str]:
    """The gold's symbolic lines, stopping at the first line that substitutes numbers into an equation."""
    out = []
    for raw in solution.splitlines():
        line = STEP.sub('', raw).strip()
        if '=' in line and re.search(r'\d+\.\d+|\d{3,}', line):
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
        g = gold_equations(first[tid]['solution']) if not d else []
        out[tid] = {'source': 'docstring' if d else ('gold' if g else None), 'equations': d or g}
    return out


def survey() -> int:
    src = sources()
    c = {'docstring': 0, 'gold': 0, None: 0}
    for v in src.values():
        c[v['source']] += 1
    L = ['# The open-book condition: where each template\'s equations come from', '',
         'Generated by `openbook.py --survey`; the rule for a symbolic line is in its docstring. Counts and equation text only '
         '(equations are the templates\' public code or the gold\'s symbolic lines, with no instance numbers).', '',
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
    (OUT / 'items.jsonl').write_text(''.join(json.dumps(r, ensure_ascii=False) + '\n' for r in rows), encoding='utf-8')
    (OUT / 'manifest.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in man), encoding='utf-8', newline='\n')
    print(f'{len(rows)} items written to {OUT / "items.jsonl"} (local) with the manifest beside it; '
          f'{sum(1 for r in man if r["source"] == "docstring")} from docstrings, {sum(1 for r in man if r["source"] == "gold")} from the gold')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', action='store_true')
    ap.add_argument('--build', action='store_true')
    a = ap.parse_args()
    if a.build:
        return build()
    return survey()


if __name__ == '__main__':
    sys.exit(main())
