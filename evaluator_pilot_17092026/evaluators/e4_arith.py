"""E4 - E3 plus a sympy check of the arithmetic the trace shows (D-085).

E3 asks whether a trace STATES each milestone. It cannot tell a milestone computed
from one copied, guessed, or produced by arithmetic that does not add up. E4 reads
every computation the trace writes (`V = 229.7 * 0.5715 = 131.27`), evaluates it
with sympy (`arith.py`, validated to 100% consistency on the 60 gold solutions
under both of its rules), and classifies each milestone E3 found:

    verified      a shown computation ends in this value, and it checks out
    contradicted  computations end in this value, and NONE of them checks out:
                  the trace states the right number but its own arithmetic does
                  not produce it - the right-number-wrong-reason case
    stated        the value appears, with no checkable computation ending in it

    e4_coverage          = (verified + stated) / required   E3 minus contradicted
    verified_coverage    = verified / required              the strict reading
    arith_consistency    = consistent claims / checked claims, over the whole trace

WHICH RULE DECIDES "CHECKS OUT" (D-097). E4 shipped at arith.ARITH_TOL, a relative
tolerance of 1%. That answers "is this number fabricated", and the answer on this
corpus was almost always no: 2 contradicted milestones in 1,494, and a null result
against E3. It is blind to the failure the experts actually marked - 164 of the 167
flawed steps behind a correct final answer are calculation slips, and a slip moves a
number by far less than 1%. At 1% the checker reaches recall 0.034 on those steps.

So the arithmetic score is now the experts' own rule - *rounding is not an error, a
wrong digit is* - `arith.Claim.ok_digit`: a displayed value must be a correct
rounding, at the precision it shows, of the numbers the trace itself shows. The 1%
rule is still computed and still reported, per claim (`ok`) and per trace
(`arith_consistency_tol1`, `milestones_contradicted_tol1`), because "fabricated" and
"mis-keyed in the last digit" are different questions and the pilot reports both.

E3 IS UNTOUCHED. `e3_coverage` is (verified + stated + contradicted) / required -
every milestone E3 reaches, however its arithmetic reads - and `meta.milestones`
keeps E3's own `reached` flag. Changing the rule cannot move E3's coverage, and
analysis/e4_analysis.py prints both columns so that is checked, not asserted.

WHY NOT SYMBOLIC EQUIVALENCE TO THE GOLD FORMULA. It needs a trace's symbols aligned
with the gold's (`V_sat`, `V_{sat}`, `V_L`), which free text does not support
reliably. What is checkable - whether the shown arithmetic produces the stated
number - is the property that catches the failure E4 exists for.
"""
from __future__ import annotations

import hashlib
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import arith                # noqa: E402
import e3_milestones as e3  # noqa: E402
import milestones as ms     # noqa: E402

ID = 'e4'
DESCRIPTION = 'E3 milestones plus sympy verification of the arithmetic the trace shows'
DEVIATIONS = []


def _sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


ARITH_RULE = 'digit'        # D-097; 'tol1' is still computed and reported beside it


def config() -> dict:
    return dict(e3.config(), evaluator=ID, arith_sha256=_sha(arith.__file__),
                e4_sha256=_sha(__file__), arith_tol=arith.ARITH_TOL,
                arith_rule=ARITH_RULE)


setup = e3.setup


def classify(hits: list[dict], claims: list, rule: str = ARITH_RULE) -> list[dict]:
    """Milestone status under one arithmetic rule. `rule` selects which verdict on a
    claim counts as checking out: the digit rule (D-097) or the 1% tolerance."""
    holds = (lambda c: c.ok_digit) if rule == 'digit' else (lambda c: c.ok)
    out = []
    for h in hits:
        h = dict(h)
        if not h['reached']:
            h['status'] = 'missed'
        else:
            ending = [c for c in claims
                      if any(ms.scaled_match(h['value'], [v], ms.DISPLAY_TOL) is not None
                             for v in c.right_value)]
            if not ending:
                h['status'] = 'stated'
            elif any(holds(c) for c in ending):
                h['status'] = 'verified'
            else:
                h['status'] = 'contradicted'
                h['evidence'] = [{'left': c.left, 'right': c.right,
                                  'left_value': c.left_value, 'right_value': c.right_value,
                                  'ulp': c.ulp, 'ok_tol1': c.ok, 'ok_digit': c.ok_digit}
                                 for c in ending[:3]]
        out.append(h)
    return out


def score(state, item: dict, trace: dict, seed: int) -> dict:
    mset = state['milestones'][item['item_id']]['milestones']
    rep = arith.check(trace['text'])
    reached = e3.reach(mset, trace['text'])
    hits = classify(reached, rep.claims)
    tol1 = classify(reached, rep.claims, rule='tol1')      # the fabrication test, kept
    n = len(hits)
    count = lambda s, hh=hits: sum(1 for h in hh if h['status'] == s)
    v, st, c = count('verified'), count('stated'), count('contradicted')
    return {
        'scores': {
            'e4_coverage': (v + st) / n if n else 0.0,
            'verified_coverage': v / n if n else 0.0,
            'e3_coverage': (v + st + c) / n if n else 0.0,
            'arith_consistency': rep.rate_digit,
            'arith_consistency_tol1': rep.rate,
            'milestones_verified': v, 'milestones_stated': st,
            'milestones_contradicted': c, 'milestones_required': n,
            'milestones_contradicted_tol1': count('contradicted', tol1),
            'milestones_verified_tol1': count('verified', tol1),
        },
        'meta': {'milestones': hits,
                 'claims_checked': rep.checked,
                 'claims_consistent': rep.consistent_digit,
                 'claims_consistent_tol1': rep.consistent,
                 'claims_digit_judged': rep.digit_judged,
                 'segments': rep.segments,
                 'inconsistent_claims': [{'line': x.line, 'left': x.left, 'right': x.right,
                                          'left_value': x.left_value, 'right_value': x.right_value,
                                          'ulp': x.ulp, 'ok_tol1': x.ok}
                                         for x in rep.claims if not x.ok_digit][:10],
                 'tribunal_reached_judges': False},
        'calls': [],
        'triggered': False,
    }
