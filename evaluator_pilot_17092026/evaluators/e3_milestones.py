"""E3 - deterministic milestone verification. No model in the loop.

For each frozen item, the milestones are the intermediate quantities a correct
derivation must reach (see milestones.py: computed by the template, stated in the
gold, not given in the question - derived by rule and reproduced byte-identically
from the frozen seed). E3 reads every number a trace writes and asks, per
milestone: does the trace state this quantity?

    milestone_coverage = milestones reached / milestones required

ORDER-FREE, as the Suggested Actions specifies: a trace that reaches the right
quantities by a different route, in a different order, loses nothing. That is the
property that makes E3 a candidate replacement for a judge deciding whether a step
is "Alternative Correct".

UNIT-AWARE. A milestone counts as reached if the trace states it under one of the
unit factors in milestones.SCALES (m vs mm, h vs min, fraction vs percent). The
factor used is recorded per milestone, so how much of the score rests on unit
normalisation can be measured rather than assumed.

TOLERANCE is 0.5%, the same display tolerance milestones.py uses to decide
whether the gold STATES a value - a trace is held to the standard that defined the
milestone. v1 used E0's 2% STEP_TOLERANCE on the reasoning that both evaluators
should judge a number alike. That was wrong for this job, and the null baseline
showed it: scored against a SIBLING item's milestones, where no reach is
legitimate, a trace still "reached" 23% of them at 2%. Sibling items put their
numbers in the same dense ranges (probabilities near 0.9, utilisations in 0-1),
and 2% of 0.9 is a wide net. Measured over the 300 gold traces:

    tol   real   null   separation   corr w/ final-answer correctness
    2.0%  0.799  0.230  0.569        0.414
    1.0%  0.776  0.120  0.656        0.454
    0.5%  0.749  0.040  0.709        0.471
    0.2%  0.720  0.023  0.697        0.446

0.5% is not a peak picked from that table: it is the rule's own tolerance, and
0.2% gives nearly the same numbers, so it sits on a plateau. E0's 2% fits its
purpose - a final answer a model may round loosely. Identifying one specific
intermediate value needs display precision.

WHAT E3 CANNOT DO, stated up front. It sees values, not reasoning: a trace that
writes the right number for the wrong reason scores as reaching it. That is the gap
E4 (symbolic equation checking) exists to close. And an item whose gold states few
intermediates (coaxial_capacitance, reynolds_number_flow_regime: 1 milestone each)
reduces E3 to a final-answer check there.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import milestones as ms  # noqa: E402

ID = 'e3'
DESCRIPTION = 'Deterministic milestone coverage: order-free, unit-aware, no model in the loop'
STEP_TOL = ms.DISPLAY_TOL          # 0.5%; see TOLERANCE above
DEVIATIONS = []


def _sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def config() -> dict:
    return {
        'evaluator': ID,
        'milestones_sha256': _sha(ms.__file__),
        'e3_sha256': _sha(__file__),
        'step_tol': STEP_TOL,
        'scales': list(ms.SCALES),
    }


def setup(dry_run: bool = False):
    manifest = os.path.join(ms._PILOT, 'slice', 'manifest.jsonl')
    return {'milestones': ms.build_all(manifest)}


def reach(milestones: list[dict], text: str) -> list[dict]:
    nums = ms.numbers(text)
    out = []
    for m in milestones:
        sc = ms.scaled_match(m['value'], nums, STEP_TOL)
        out.append({'id': m['id'], 'value': m['value'], 'reached': sc is not None, 'scale': sc})
    return out


def score(state, item: dict, trace: dict, seed: int) -> dict:
    mset = state['milestones'][item['item_id']]['milestones']
    hits = reach(mset, trace['text'])
    n = len(hits)
    got = sum(h['reached'] for h in hits)
    return {
        'scores': {'milestone_coverage': got / n if n else 0.0,
                   'milestones_reached': got, 'milestones_required': n},
        'meta': {'milestones': hits,
                 'reached_by_unit_scaling': sum(1 for h in hits if h['reached'] and h['scale'] != 1.0),
                 'tribunal_reached_judges': False},
        'calls': [],
        'triggered': False,
    }
