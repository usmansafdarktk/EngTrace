# Layer 2 - human certification results, round 6

Generated 2026-10-06T19:44:35+00:00 by `score.py` from 3 label rows by 3 experts.

Round 6 re-certifies the 1 templates changed after round 5, built at git `740ee25`: no planted defects, and a fresh hand-check instance (seed 2601) that no expert saw before.

## Verdicts per expert

No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).

| Expert | Branch | Items | Approved | Rejected |
|---|---|---:|---:|---:|
| che-1 | chemical | 1 | 1 | 0 |
| che-2 | chemical | 1 | 1 | 0 |
| che-3 | chemical | 1 | 1 | 0 |

## Round 5 against round 6

Each expert's verdict on the same template in both rounds: A approve, R reject.

| Template | Round 5 | Round 6 | Outcome |
|---|---|---|---|
| template_work_isothermal_virial | che-1 A che-2 R che-3 A | che-1 A che-2 A che-3 A | approved by all |

Templates: 1 approved by all three, 0 approved by majority, 0 rejected by majority.
Of the 1 round-5 rejections of these templates, 1 became approvals and 0 stayed rejections; 0 verdicts went the other way, from approve to reject.

## Hand checks

3 hand checks with a comparable number: 3 matched the template within 1% (100%).
Of the 0 mismatches, 0 ended in a rejection and 0 in an approval (the expert found their own slip, or judged the difference immaterial).

## Agreement among the three experts of a branch (real templates)

| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | AC2 phys | AC2 math | AC2 ped |
|---|---:|---:|---:|---:|---:|---:|---:|
| chemical | 1 | 1.000 | 1.000 | 100% | 1.000 | 0.940 | 1.000 |
| all | 1 | 1.000 | 1.000 | 100% | 1.000 | 0.940 | 1.000 |

Kappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.

## Against the screening panel (pass 2)

Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.

## Time spent (app rows only)

| Expert | Items | Median minutes | Under 2 min |
|---|---:|---:|---:|
| che-1 | 1 | 1.5 | 1 |
| che-2 | 1 | 2.5 | 0 |
| che-3 | 1 | 1.4 | 1 |

## Fix list: real templates rejected by at least one expert

none
