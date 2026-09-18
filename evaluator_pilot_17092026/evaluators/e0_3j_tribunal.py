"""E0-3J - E0 with the third judge actually connected.

E0 as published is a TWO-judge panel, because `_tier2_tribunal_batch` gates Google
on `genai.get_model_info(...)`, a function `google.generativeai` does not have, and
swallows the AttributeError (FINDINGS E0-F6). This variant is E0 with that one
accident removed, so the difference between the two columns measures exactly what
the dropped judge was worth.

HOW THE JUDGE IS RESTORED. By supplying the missing function, not by editing the
framework:

    fw.genai.get_model_info = lambda name: fw.genai.get_model(...)

The framework's own line then succeeds and appends 'google' by its own logic. No
change to the provider loop, the prompt, the parsing, the vote or the metrics. The
shim calls the SDK's real `get_model`, so a genuinely unavailable model still
raises and Google still drops out - the check does what it was written to do
instead of always failing.

WHY THIS IS NOT E0. It is the panel the paper describes, not the panel the code
runs. E0 stays the anchor; this is the counterfactual beside it.

WHAT THE THIRD VOTE CHANGES. With two judges the majority test `count > total/2`
demands unanimity, so any disagreement falls to the conservative `min()`
tie-break. With three, a 2-1 split carries the majority. So restoring Google can
only ever raise a split score or leave it alone - which is why the comparison is
worth its cost on a panel that decides 93% of correct traces.
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import e0_tribunal as e0  # noqa: E402

ID = 'e0_3j'
DESCRIPTION = 'E0 with the Google judge connected: the three-judge panel the paper describes'

JUDGES = e0.JUDGES
DEVIATIONS = e0.DEVIATIONS + [
    'D4 genai.get_model_info supplied so the framework\'s own Google check can pass; '
    'without it the published code silently runs a two-judge panel (E0-F6)',
]

libraries = e0.libraries
compute = e0.compute
score = e0.score


def config() -> dict:
    cfg = e0.config()
    cfg['evaluator'] = ID
    cfg['deviations'] = DEVIATIONS
    cfg['judges_connected'] = 3
    return cfg


def setup(dry_run: bool = False):
    state = e0.setup(dry_run=dry_run)
    genai = state.genai

    def get_model_info(name):
        # The SDK's real lookup, under the name the framework calls. An unavailable
        # model still raises here, so the framework's check keeps its meaning.
        return genai.get_model(name if name.startswith('models/') else 'models/' + name)

    genai.get_model_info = get_model_info
    if not dry_run:
        probe = get_model_info(JUDGES['google']['model'])
        print('  google judge connected: %s' % getattr(probe, 'name', probe))
    return state
