"""Milestones for E3: the intermediate quantities a correct derivation must reach.

WHERE THEY COME FROM. Not from editing templates. Each frozen item is regenerated
at its recorded seed with the repo's own frame-local capture
(`tests/template_integrity/core.generate(..., capture=True)`), and the regenerated
question and solution are checked BYTE-IDENTICAL to the frozen manifest before a
single value is used. All 60 items pass that check, so the captured values are
provably that item's internals, and the frozen text cannot move.

WHAT COUNTS AS A MILESTONE. A rule, not a hand-picked list, so it can be audited
and so it applies unchanged to templates nobody has looked at:

    a value the template COMPUTED
    that the gold solution STATES
    and the question does NOT give.

"Computed" rules out inputs; "stated in the gold" rules out scratch variables no
solution mentions; "not in the question" rules out restated givens. A value counts
as stated if some number in the text is within DISPLAY_TOL of it, because the text
shows values rounded for display (V_sat = 113.5517 is written 113.55).

WHAT IS EXCLUDED, and why. Values within 1e-9 of 0, 1 or 2 match almost any text
and would make every trace look like it reached them. Loop bookkeeping (`_attempt`,
`precision`, indices) is not a quantity. Lookup tables captured whole (e.g. a dict
of every fluid's density) are inputs, not derivation, and are dropped by the
not-in-the-question rule only when the chosen entry is quoted, so they are also
dropped by name here.
"""
from __future__ import annotations

import json
import math
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_ROOT = os.path.dirname(_PILOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

DISPLAY_TOL = 0.005          # 0.5%: covers rounding to 2-4 significant figures
NUM = re.compile(r'(?<![A-Za-z_])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
TRIVIAL = (0.0, 1.0, 2.0)
BOOKKEEPING = re.compile(r'(^_|precision|attempt|choice|^idx|^i$|^j$|^k$|^n_iter|max_elements|'
                         r'tolerance|termination|frame|\bdp\b)')
# Iteration state. A milestone must be reachable by ANY correct path, and the
# trajectory of one iterative scheme is not: a trace that starts from a different
# guess passes through different y_prev / y_curr values and converges to the same
# answer. Keeping the trajectory would punish exactly the alternative correct path
# a deterministic evaluator has to be shown not to punish.
TRAJECTORY = re.compile(r'(_prev$|_curr$|_next$|^change$|^last_change$|^[yg]_[cp]$|_s$|^updates$)')

# Unit normalisation. Values related by one of these factors are the same quantity
# in different units (m vs mm, h vs min, fraction vs percent).
SCALES = (1.0, 1e3, 1e-3, 1e6, 1e-6, 1e9, 1e-9, 60.0, 1 / 60.0, 3600.0, 1 / 3600.0, 100.0, 0.01)


def numbers(text: str) -> list[float]:
    """Every number written in a text. `×10^` and `10^-3` forms are folded in."""
    t = text.replace(',', '').replace('−', '-')
    t = re.sub(r'(\d)\s*[×x\*]\s*10\s*\^\s*\{?\s*([-+]?\d+)\s*\}?', r'\1e\2', t)
    out = []
    for m in NUM.finditer(t):
        try:
            v = float(m.group(0))
        except ValueError:
            continue
        if math.isfinite(v):
            out.append(v)
    return out


def close(a: float, b: float, tol: float) -> bool:
    if abs(b) < 1e-12:
        return abs(a) < 1e-12
    return abs(a - b) / abs(b) <= tol


def stated(value: float, nums: list[float], tol: float = DISPLAY_TOL) -> bool:
    return any(close(n, value, tol) for n in nums)


def scaled_match(value: float, nums: list[float], tol: float):
    """The unit factor under which some number equals `value`, or None."""
    for sc in SCALES:
        target = value * sc
        if any(close(n, target, tol) for n in nums):
            return sc
    return None


def build(item: dict) -> dict:
    """Milestones for one frozen item. Raises if the item does not reproduce."""
    from tests.template_integrity.core import discover, generate
    ref = next(r for r in discover() if r.template_id == item['template_id'])
    inst = generate(ref, item['seed'], capture=True)
    if inst.question != item['question'] or inst.solution != item['solution']:
        raise ValueError('%s does not reproduce byte-identically; refusing to derive '
                         'milestones from a different instance' % item['item_id'])

    gold, given = numbers(item['solution']), numbers(item['question'])
    seen, ms = [], []
    for name, v in sorted(inst.values.items()):
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            continue
        # '[' means a container element: a lookup table, a trajectory row, or a
        # component of a given vector. None is a derivation step in its own right.
        if '[' in name or BOOKKEEPING.search(name) or TRAJECTORY.search(name):
            continue
        if any(abs(v - t) < 1e-9 for t in TRIVIAL):
            continue
        if not stated(v, gold):
            continue
        # A given restated in other units (mm -> m) is a conversion, not a step.
        if scaled_match(v, given, DISPLAY_TOL) is not None:
            continue
        # The same quantity under two names or two units counts once.
        if scaled_match(v, seen, DISPLAY_TOL) is not None:
            continue
        seen.append(v)
        ms.append({'id': name, 'value': float(v)})
    return {'item_id': item['item_id'], 'template_id': item['template_id'],
            'milestones': ms, 'gold_numbers': len(gold), 'given_numbers': len(given)}


def build_all(manifest_path: str) -> dict:
    items = [json.loads(l) for l in open(manifest_path, encoding='utf-8')]
    return {it['item_id']: build(it) for it in items}


if __name__ == '__main__':
    M = build_all(os.path.join(_PILOT, 'slice', 'manifest.jsonl'))
    by_t = {}
    for m in M.values():
        by_t.setdefault(m['template_id'], []).append(m)
    for t, ms in sorted(by_t.items()):
        counts = [len(m['milestones']) for m in ms]
        ids = [x['id'] for x in ms[0]['milestones']]
        print('%-40s %s milestones  e.g. %s' % (t.replace('template_', ''), counts, ', '.join(ids[:8])))
