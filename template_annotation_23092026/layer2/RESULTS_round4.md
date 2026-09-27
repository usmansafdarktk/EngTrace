# Layer 2 - human certification results, round 4

Generated 2026-09-27T21:33:10+00:00 by `score.py` from 6 label rows by 3 experts.

Round 4 re-certifies the 2 templates changed after round 3, built at git `fb88cf1`: no planted defects, and a fresh hand-check instance (seed 2401) that no expert saw before.

## Verdicts per expert

No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).

| Expert | Branch | Items | Approved | Rejected |
|---|---|---:|---:|---:|
| che-1 | chemical | 2 | 2 | 0 |
| che-2 | chemical | 2 | 2 | 0 |
| che-3 | chemical | 2 | 2 | 0 |

## Hand checks

6 hand checks with a comparable number: 6 matched the template within 1% (100%).
Of the 0 mismatches, 0 ended in a rejection and 0 in an approval (the expert found their own slip, or judged the difference immaterial).

## Agreement among the three experts of a branch (real templates)

| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | AC2 phys | AC2 math | AC2 ped |
|---|---:|---:|---:|---:|---:|---:|---:|
| chemical | 2 | 1.000 | 1.000 | 100% | 0.970 | 1.000 | 1.000 |
| all | 2 | 1.000 | 1.000 | 100% | 0.970 | 1.000 | 1.000 |

Kappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.

## Against the screening panel (pass 2)

Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.

## Time spent (app rows only)

| Expert | Items | Median minutes | Under 2 min |
|---|---:|---:|---:|
| che-1 | 2 | 1.8 | 2 |
| che-2 | 2 | 2.0 | 1 |
| che-3 | 2 | 1.8 | 1 |

## Fix list: real templates rejected by at least one expert

none
