# Layer 2 - human certification results, round 5

Generated 2026-10-06T18:47:30+00:00 by `score.py` from 6 label rows by 3 experts.

Round 5 re-certifies the 2 templates changed after round 4, built at git `740ee25`: no planted defects, and a fresh hand-check instance (seed 2501) that no expert saw before.

## Verdicts per expert

No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).

| Expert | Branch | Items | Approved | Rejected |
|---|---|---:|---:|---:|
| che-1 | chemical | 2 | 2 | 0 |
| che-2 | chemical | 2 | 1 | 1 |
| che-3 | chemical | 2 | 2 | 0 |

## Round 4 against round 5

Each expert's verdict on the same template in both rounds: A approve, R reject.

| Template | Round 4 | Round 5 | Outcome |
|---|---|---|---|
| template_adiabatic_flame_temperature | che-1 A che-2 A che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_work_isothermal_virial | che-1 - che-2 - che-3 - | che-1 A che-2 R che-3 A | approved by majority, rejected by 1 |

Templates: 1 approved by all three, 1 approved by majority, 0 rejected by majority.
Of the 0 round-4 rejections of these templates, 0 became approvals and 0 stayed rejections; 0 verdicts went the other way, from approve to reject.

## Hand checks

6 hand checks with a comparable number: 6 matched the template within 1% (100%).
Of the 0 mismatches, 0 ended in a rejection and 0 in an approval (the expert found their own slip, or judged the difference immaterial).

## Agreement among the three experts of a branch (real templates)

| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | AC2 phys | AC2 math | AC2 ped |
|---|---:|---:|---:|---:|---:|---:|---:|
| chemical | 2 | -0.200 | 0.538 | 50% | 0.856 | 1.000 | 1.000 |
| all | 2 | -0.200 | 0.538 | 50% | 0.856 | 1.000 | 1.000 |

Kappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.

## Against the screening panel (pass 2)

Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.

## Time spent (app rows only)

| Expert | Items | Median minutes | Under 2 min |
|---|---:|---:|---:|
| che-1 | 2 | 1.7 | 2 |
| che-2 | 2 | 2.3 | 0 |
| che-3 | 2 | 1.6 | 1 |

## Fix list: real templates rejected by at least one expert

- **template_work_isothermal_virial** rejected by 1 of 3:
  - che-2 [physics or scenario implausible]: Physically implausible states: the Vr>=2 filter keeps only Tr of about 2-3, so organics are compressed far above their thermal-decomposition range. Instance 2 compresses n-pentane at 1135.26 K (862 C) up to 142.21 bar, where it would crack in well under a second, so a quasi-static reversible compression of that species cannot happen. Cap T for each substance (organics below about 700-750 K) or drop substances with no valid state. Minor: in Step 3 the '*' operators render as markdown italics, e.g. '40.0162621.57/(0.083141135.26)'.
