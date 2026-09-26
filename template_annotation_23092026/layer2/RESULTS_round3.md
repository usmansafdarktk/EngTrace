# Layer 2 - human certification results, round 3

Generated 2026-09-26T13:54:23+00:00 by `score.py` from 15 label rows by 9 experts.

Round 3 re-certifies the 5 templates changed after round 2, built at git `9c5dc45`: no planted defects, and a fresh hand-check instance (seed 2301) that no expert saw before.

## Verdicts per expert

No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).

| Expert | Branch | Items | Approved | Rejected |
|---|---|---:|---:|---:|
| che-1 | chemical | 1 | 1 | 0 |
| che-2 | chemical | 1 | 1 | 0 |
| che-3 | chemical | 1 | 1 | 0 |
| ele-1 | electrical | 1 | 1 | 0 |
| ele-2 | electrical | 1 | 1 | 0 |
| ele-3 | electrical | 1 | 1 | 0 |
| mec-1 | mechanical | 3 | 3 | 0 |
| mec-2 | mechanical | 3 | 3 | 0 |
| mec-3 | mechanical | 3 | 3 | 0 |

## Round 2 against round 3

Each expert's verdict on the same template in both rounds: A approve, R reject.

| Template | Round 2 | Round 3 | Outcome |
|---|---|---|---|
| template_basic_stress_strain | mec-1 R mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |
| template_finite_convolution | ele-1 R ele-2 A ele-3 A | ele-1 A ele-2 A ele-3 A | approved by all |
| template_multi_segment_rod | mec-1 A mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |
| template_newtons_law_shear_stress | che-1 A che-2 A che-3 R | che-1 A che-2 A che-3 A | approved by all |
| template_utube_manometer | mec-1 R mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |

Templates: 5 approved by all three, 0 approved by majority, 0 rejected by majority.
Of the 10 round-2 rejections of these templates, 10 became approvals and 0 stayed rejections; 0 verdicts went the other way, from approve to reject.

## Hand checks

15 hand checks with a comparable number: 13 matched the template within 1% (87%).
Of the 2 mismatches, 0 ended in a rejection and 2 in an approval (the expert found their own slip, or judged the difference immaterial).

Real templates approved despite a mismatch (check these):
- template_utube_manometer by mec-1: entered `6831 Pa`; note: -
- template_utube_manometer by mec-3: entered `6831 Pa`; note: -

## Agreement among the three experts of a branch (real templates)

| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | AC2 phys | AC2 math | AC2 ped |
|---|---:|---:|---:|---:|---:|---:|---:|
| chemical | 1 | 1.000 | 1.000 | 100% | 1.000 | 1.000 | 1.000 |
| electrical | 1 | 1.000 | 1.000 | 100% | 1.000 | 1.000 | 1.000 |
| mechanical | 3 | 1.000 | 1.000 | 100% | 0.960 | 1.000 | 1.000 |
| all | 5 | 1.000 | 1.000 | 100% | 0.979 | 1.000 | 1.000 |

Kappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.

## Against the screening panel (pass 2)

Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.

## Time spent (app rows only)

| Expert | Items | Median minutes | Under 2 min |
|---|---:|---:|---:|
| che-1 | 1 | 1.9 | 1 |
| che-2 | 1 | 2.4 | 0 |
| che-3 | 1 | 2.1 | 0 |
| ele-1 | 1 | 2.5 | 0 |
| ele-2 | 1 | 3.0 | 0 |
| ele-3 | 1 | 2.2 | 0 |
| mec-1 | 3 | 1.6 | 2 |
| mec-2 | 3 | 1.5 | 2 |
| mec-3 | 3 | 1.9 | 2 |

## Fix list: real templates rejected by at least one expert

none
