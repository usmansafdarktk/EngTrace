"""The two subsamples the plan names, fixed here once so the harness, the paraphrase pipeline and the
analysis select the same items.

  PARAPHRASE  ANALYSIS_PLAN Q5: the 1st, 6th and 11th of each template's 15, in manifest order
              (450 items)
  REPEAT      the plan's decoding-variance repeat, "300 items (2 per template)": the 1st and 8th of
              each template's 15, in manifest order. The plan does not say which two; this rule was
              fixed on 2026-09-29, before any repeat ran (D-141)

"Manifest order" is the order of manifest.jsonl, in which each template's 15 rows run in ascending
instance index; the indices have gaps (D-116), so positions are counted, not indices.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'manifest.jsonl'
PARAPHRASE = (0, 5, 10)
REPEAT = (0, 7)
REPEAT_VARIANTS = ('repeat1', 'repeat2', 'repeat3')


def pick(positions: tuple[int, ...]) -> list[str]:
    """Item ids at these positions of each template's rows, templates in manifest order."""
    rows = defaultdict(list)
    for line in MANIFEST.read_text(encoding='utf-8').splitlines():
        r = json.loads(line)
        rows[r['template_id']].append(r['item_id'])
    out = []
    for ids in rows.values():
        if len(ids) != 15:
            raise SystemExit(f'a template has {len(ids)} items, not 15: the manifest is not the frozen pool')
        out += [ids[p] for p in positions]
    return out


def paraphrase_ids() -> list[str]:
    return pick(PARAPHRASE)


def repeat_ids() -> list[str]:
    return pick(REPEAT)


if __name__ == '__main__':
    print(len(paraphrase_ids()), 'paraphrase items;', len(repeat_ids()), 'repeat items')
