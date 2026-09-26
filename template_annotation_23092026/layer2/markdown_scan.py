"""Which questions and solutions change when the review app renders them as Markdown.

    python -m template_annotation_23092026.layer2.markdown_scan              # 200 seeds per template
    python -m template_annotation_23092026.layer2.markdown_scan --seeds 50

Why this exists. The review app (app.py) shows every question and solution with
st.markdown, but the models being benchmarked receive the raw text. Where a text carries
characters Markdown treats as syntax, the expert saw something other than what the model
reads: paired asterisks in an unspaced product such as 2*pi*f*t turn into italics and the
multiplication signs vanish; a pair of dollar signs around currency becomes inline math and
the dollar signs vanish; a value marked only by asterisks loses its marker. Round 2 found
this through template_finite_convolution (D-107); the fixing agent's scan of 2026-09-26 is
the logic here, extended from questions to solutions.

What it parses: CommonMark via markdown-it-py (installed with Streamlit, through `rich`),
plus the extras in Streamlit's renderer: inline math between dollar signs, GFM
strikethrough, directives, and ' -- ' shown as a dash. Questions are checked for every
construct the renderer changes; solutions, whose bold step headers and lists are intended,
only for the constructs that remove characters (italics from paired asterisks, dollar
pairs, escapes, entities, strikethrough, directives).

Writes markdown_scan.json and markdown_scan.md beside this file. Informational: it gates
nothing. It does not flag the collapse of single newlines, which Streamlit renders as
spaces in every solution.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from markdown_it import MarkdownIt  # noqa: E402

from tests.template_integrity.core import discover, generate  # noqa: E402

MD = MarkdownIt('commonmark')
BLOCK = {'heading_open': 'heading', 'bullet_list_open': 'bullet list', 'ordered_list_open': 'ordered list',
         'blockquote_open': 'blockquote', 'code_block': 'indented code block', 'fence': 'fenced code',
         'html_block': 'raw HTML', 'hr': 'thematic break'}
INLINE = {'code_inline': 'code span', 'html_inline': 'raw HTML', 'link_open': 'link', 'image': 'image'}
BS = chr(92)
ESC = re.compile(re.escape(BS) + r'[!-/:-@\[-`{-~]')
ENT = re.compile(r'&(?:#\d+|#[xX][0-9a-fA-F]+|[A-Za-z][A-Za-z0-9]{1,31});')
DIRECTIVE = re.compile(r'(?<!:):[A-Za-z][\w-]*[\[{]')
TILDE_PAIR = re.compile(r'(?<![~\w])~(?=[^\s~])[^~\n]*?(?<=[^\s~])~(?![~\w])')
DASH = re.compile(r'(^|\s)--(\s|$)')
DOLLAR = re.compile('(?<!' + re.escape(BS) + ')' + re.escape('$'))
# the constructs that remove characters from what the reader sees
LOSSY = {'emphasis (italic)', 'dollar pair (inline math)', 'backslash escape (backslash dropped)',
         'HTML entity (decoded)', 'tilde pair (strikethrough)', 'directive (label dropped)'}


def fragments(text: str) -> dict[str, str]:
    """{category: first fragment} for every construct the renderer changes in `text`."""
    found: dict[str, str] = {}
    for t in MD.parse(text):
        if t.type in BLOCK:
            found.setdefault(BLOCK[t.type], (t.content or t.markup or '')[:60])
        for i, c in enumerate(t.children or []):
            if c.type in ('em_open', 'strong_open'):
                close = c.type.replace('_open', '_close')
                depth, j, inner = 0, i + 1, []
                ch = t.children
                while j < len(ch):
                    if ch[j].type == c.type:
                        depth += 1
                    if ch[j].type == close:
                        if depth == 0:
                            break
                        depth -= 1
                    inner.append(ch[j].content or ch[j].markup)
                    j += 1
                cat = 'emphasis (italic)' if c.type == 'em_open' else 'strong (bold)'
                found.setdefault(cat, c.markup + ''.join(inner)[:50] + c.markup)
            elif c.type in INLINE:
                found.setdefault(INLINE[c.type], (c.content or '')[:60])
    for para in re.split(r'\n\s*\n', text):
        if len(DOLLAR.findall(para)) >= 2:
            found.setdefault('dollar pair (inline math)', para[para.find('$'):][:70])
        for pat, cat in ((TILDE_PAIR, 'tilde pair (strikethrough)'), (DIRECTIVE, 'directive (label dropped)'),
                         (DASH, "' -- ' (shown as a dash)")):
            m = pat.search(para)
            if m:
                found.setdefault(cat, para[max(0, m.start() - 15):m.end() + 15])
    m = ESC.search(text)
    if m:
        found.setdefault('backslash escape (backslash dropped)', text[max(0, m.start() - 15):m.end() + 15])
    m = ENT.search(text)
    if m:
        found.setdefault('HTML entity (decoded)', m.group(0))
    return found


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', type=int, default=200)
    a = ap.parse_args()
    res: dict = collections.defaultdict(lambda: collections.defaultdict(lambda: {'seeds': 0, 'example': ''}))
    errors = collections.Counter()
    refs = discover()
    for ref in refs:
        for s in range(a.seeds):
            inst = generate(ref, s, capture=False)
            if not inst.ok:
                errors[ref.template_id] += 1
                continue
            cats = {('question', k): v for k, v in fragments(inst.question).items()}
            cats.update({('solution', k): v for k, v in fragments(inst.solution).items() if k in LOSSY})
            for (part, cat), frag in cats.items():
                r = res[ref.template_id][f'{part}: {cat}']
                r['seeds'] += 1
                if not r['example']:
                    r['example'] = f'seed {s}: {frag}'
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    branch_of = {r.template_id: r.branch for r in refs}
    report = {'generated': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'git_head': head,
              'seeds': a.seeds, 'templates': len(refs), 'errors': dict(errors),
              'findings': {t: dict(v) for t, v in sorted(res.items())}}
    (HERE / 'markdown_scan.json').write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding='utf8')

    lossy_q = {t for t, v in res.items() if any(k.split(': ', 1)[1] in LOSSY for k in v if k.startswith('question'))}
    lossy_s = {t for t, v in res.items() if any(k.startswith('solution') for k in v)}
    out = ['# Markdown rendering scan\n',
           f"Generated {report['generated']} by `markdown_scan.py` at git `{head[:10]}`, {a.seeds} seeds per template "
           f"over {len(refs)} templates. The review app renders questions and solutions as Markdown; the models "
           f"read the raw text. A lossy construct removes characters from what the expert sees.\n",
           '| | Templates |\n|---|---:|',
           f'| any construct in the question | {sum(1 for v in res.values() if any(k.startswith("question") for k in v))} |',
           f'| lossy construct in the question | {len(lossy_q)} |',
           f'| lossy construct in the solution | {len(lossy_s)} |',
           f'| lossy in either | {len(lossy_q | lossy_s)} |\n',
           '| Template | Branch | Construct | Seeds | Example |\n|---|---|---|---:|---|']
    for t in sorted(res, key=lambda t: (branch_of[t], t)):
        for k, r in sorted(res[t].items()):
            ex = r['example'].replace('|', '\\|').replace('\n', ' ')
            out.append(f"| `{t}` | {branch_of[t].split('_')[0]} | {k} | {r['seeds']} | {ex} |")
    (HERE / 'markdown_scan.md').write_text('\n'.join(out) + '\n', encoding='utf8')
    print(f'{len(refs)} templates x {a.seeds} seeds; lossy in question {len(lossy_q)}, in solution {len(lossy_s)}, '
          f'either {len(lossy_q | lossy_s)}; generation errors {dict(errors)}')


if __name__ == '__main__':
    main()
