# Layer 2 - human certification results, round 2

Generated 2026-09-25T19:49:56+00:00 by `score.py` from 66 label rows by 9 experts.

Round 2 re-certifies the 22 templates changed after round 1, built at git `bc8aaf3`: no planted defects, and a fresh hand-check instance (seed 2201) that no expert saw before.

## Verdicts per expert

No planted defects in this round; the review's sensitivity was measured in round 1 (`RESULTS.md`).

| Expert | Branch | Items | Approved | Rejected |
|---|---|---:|---:|---:|
| che-1 | chemical | 9 | 9 | 0 |
| che-2 | chemical | 9 | 9 | 0 |
| che-3 | chemical | 9 | 8 | 1 |
| ele-1 | electrical | 3 | 2 | 1 |
| ele-2 | electrical | 3 | 3 | 0 |
| ele-3 | electrical | 3 | 3 | 0 |
| mec-1 | mechanical | 10 | 8 | 2 |
| mec-2 | mechanical | 10 | 7 | 3 |
| mec-3 | mechanical | 10 | 7 | 3 |

## Round 1 against round 2

Each expert's verdict on the same template in both rounds: A approve, R reject.

| Template | Round 1 | Round 2 | Outcome |
|---|---|---|---|
| template_angle_of_twist | mec-1 R mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |
| template_annulus_flowrate | che-1 R che-2 R che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_axial_deformation | mec-1 R mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |
| template_basic_buoyant_force | mec-1 R mec-2 A mec-3 A | mec-1 A mec-2 A mec-3 A | approved by all |
| template_basic_stress_strain | mec-1 R mec-2 R mec-3 A | mec-1 R mec-2 R mec-3 R | rejected by majority, 3 of 3 |
| template_falling_film_max_velocity | che-1 R che-2 R che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_finite_convolution | ele-1 R ele-2 A ele-3 A | ele-1 R ele-2 A ele-3 A | approved by majority, rejected by 1 |
| template_flow_system_molar_flow_rates | che-1 R che-2 A che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_gas_viscosity_kinetic_theory | che-1 R che-2 R che-3 R | che-1 A che-2 A che-3 A | approved by all |
| template_hagen_poiseuille_flowrate | che-1 R che-2 R che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_multi_segment_rod | mec-1 R mec-2 A mec-3 A | mec-1 A mec-2 R mec-3 R | rejected by majority, 2 of 3 |
| template_newtons_law_shear_stress | che-1 R che-2 R che-3 R | che-1 A che-2 A che-3 R | approved by majority, rejected by 1 |
| template_pitzer_correlation_z | che-1 A che-2 A che-3 R | che-1 A che-2 A che-3 A | approved by all |
| template_power_law_fluid_shear | che-1 R che-2 R che-3 A | che-1 A che-2 A che-3 A | approved by all |
| template_statically_indeterminate | mec-1 R mec-2 R mec-3 A | mec-1 A mec-2 A mec-3 A | approved by all |
| template_statically_indeterminate_shaft | mec-1 R mec-2 R mec-3 R | mec-1 A mec-2 A mec-3 A | approved by all |
| template_time_to_phasor | ele-1 R ele-2 R ele-3 A | ele-1 A ele-2 A ele-3 A | approved by all |
| template_utube_manometer | mec-1 R mec-2 R mec-3 R | mec-1 R mec-2 R mec-3 R | rejected by majority, 3 of 3 |
| template_vibration_transmissibility | mec-1 R mec-2 A mec-3 A | mec-1 A mec-2 A mec-3 A | approved by all |
| template_volumetric_flow_rate | mec-1 R mec-2 R mec-3 A | mec-1 A mec-2 A mec-3 A | approved by all |
| template_wave_parameters_basic | ele-1 R ele-2 R ele-3 A | ele-1 A ele-2 A ele-3 A | approved by all |
| template_work_isothermal_virial | che-1 A che-2 A che-3 R | che-1 A che-2 A che-3 A | approved by all |

Templates: 17 approved by all three, 2 approved by majority, 3 rejected by majority.
Of the 43 round-1 rejections of these templates, 36 became approvals and 7 stayed rejections; 3 verdicts went the other way, from approve to reject.

## Hand checks

66 hand checks with a comparable number: 66 matched the template within 1% (100%).
Of the 0 mismatches, 0 ended in a rejection and 0 in an approval (the expert found their own slip, or judged the difference immaterial).

## Agreement among the three experts of a branch (real templates)

| Branch | Templates with 3 verdicts | Fleiss kappa (Approve/Reject) | Gwet AC1 | Percent agreement | AC2 phys | AC2 math | AC2 ped |
|---|---:|---:|---:|---:|---:|---:|---:|
| chemical | 9 | -0.038 | 0.920 | 89% | 0.985 | 0.995 | 0.987 |
| electrical | 3 | -0.125 | 0.723 | 67% | 0.984 | 1.000 | 0.944 |
| mechanical | 10 | 0.830 | 0.891 | 90% | 0.985 | 0.974 | 0.990 |
| all | 22 | 0.646 | 0.878 | 86% | 0.983 | 0.987 | 0.983 |

Kappa collapses when nearly everything is approved (the prevalence artefact Appendix K discusses); AC1/AC2 do not, and the plant detection rate above is the sensitivity figure kappa cannot give.

## Against the screening panel (pass 2)

Not computed for this round: the panel judged these templates before the fixes and has not re-judged them.

## Time spent (app rows only)

| Expert | Items | Median minutes | Under 2 min |
|---|---:|---:|---:|
| che-1 | 9 | 1.3 | 9 |
| che-2 | 9 | 1.4 | 9 |
| che-3 | 9 | 1.2 | 9 |
| ele-1 | 3 | 1.5 | 2 |
| ele-2 | 3 | 1.4 | 3 |
| ele-3 | 3 | 1.3 | 2 |
| mec-1 | 10 | 1.2 | 10 |
| mec-2 | 10 | 1.1 | 10 |
| mec-3 | 10 | 1.1 | 10 |

## Fix list: real templates rejected by at least one expert

- **template_basic_stress_strain** rejected by 3 of 3:
  - mec-1 (round 1: reject) [physics or scenario implausible]: Load, section and elongation are sampled independently with no material, so the posed data imply impossible moduli: instance 2 (439 kN on a 36 mm bar, 1.17 mm over 3.86 m) gives sigma = 431.29 MPa with epsilon = 3.031e-4, i.e. E = sigma/epsilon = 1423 GPa, stiffer than any engineering material (steel would stretch 8.3 mm), and instance 4 implies E = 565 GPa at 526 MPa. The elongation must be derived from, or bounded by, a real material's E (and the stress by its strength); the stress and strain arithmetic itself is correct.
  - mec-2 (round 1: reject) [physics or scenario implausible]: Load, size and elongation are drawn independently, so the stress and strain imply impossible materials: instance 2 (36 mm bar, 439 kN, 1.17 mm over 3.86 m) gives sigma = 431.29 MPa with epsilon = 3.031e-4, an implied E = 1423 GPa (stiffer than diamond), and instance 4 gives 525.99 MPa at 9.310e-4, E = 565 GPa. The elongation (or load) must be bounded by a realistic modulus and allowable stress so sigma/epsilon matches a real material.
  - mec-3 (round 1: approve) [physics or scenario implausible]: Load, bar size and elongation are drawn independently, so the implied modulus E = sigma/epsilon is often impossible: instance 2 (439 kN on a 36 mm bar, 1.17 mm over 3.86 m) gives 431.29 MPa / 3.031e-04 = 1423 GPa, and instance 4 gives 525.99 MPa / 9.310e-04 = 565 GPa, both above any engineering material (steel is about 200 GPa and would stretch about 8.3 mm in instance 2). Derive the elongation from a named material's E, or bound the draw so sigma/epsilon falls within a real modulus range.
- **template_finite_convolution** rejected by 1 of 3:
  - ele-1 (round 1: reject) [question ambiguous or unsolvable]: Question ambiguous as rendered: the *origin* asterisks render as italics, and the fallback "x[0] = v" does not locate the origin when v repeats (instance 1 x = {2, -2, -2} with x[0] = -2; instance 3 x = {-3, 3, 2, -3} with x[0] = -3; instance 4 x = {3, -1, 3} with x[0] = 3), so e.g. in instance 1 a reader may take the last -2 and get y[0] = 2 instead of gold y[0] = 8. Convolution values and indices are otherwise all correct; state the start index (e.g. "x[n] starts at n = -1") instead of relying on asterisks.
- **template_multi_segment_rod** rejected by 2 of 3:
  - mec-2 (round 1: approve) [final answer wrong or wrong unit]: In US units each segment deformation is rounded to 4 dp (one or two significant figures) before summing, so the total is wrong: instance 3 prints 0.0036 + 0.0027 - 0.0054 = 0.0009 in, but the exact segment values 0.003567, 0.002666 and -0.005434 in sum to 0.000800 in, so the printed answer is 12.5% off. Carry the segment deformations to at least four significant figures (e.g. 0.003567 in) and sum those.
  - mec-3 (round 1: approve) [arithmetic: a step does not follow, final answer wrong or wrong unit]: In the US branch the segment deformations are printed at a fixed 4 dp in inches and then summed, which leaves one significant figure and a wrong final answer. Instance 3: delta1 = 0.003567, delta2 = 0.002666, delta3 = -0.005434 in sum to 0.000800 in, but the template gives 0.0036 + 0.0027 - 0.0054 = 0.0009 in, 12.5% too large. Quote the deformations to at least four significant figures (e.g. 3.567e-3 in) before summing.
- **template_newtons_law_shear_stress** rejected by 1 of 3:
  - che-3 (round 1: reject) [physics or scenario implausible]: No laminar check: the linear Couette profile (tau = mu*V/Y) needs Re = rho*V*Y/mu below about 1500, but low-viscosity draws are far into turbulence. Examples: water, 1.59 m/s, 1.96 cm gap gives Re about 3.1e4; gasoline, 1.54 m/s, 2.22 cm gives about 8.5e4. For these the printed tau and F underestimate the real drag. Add a laminar redraw like the sibling annulus and film templates have. Minor: F is printed to 10 significant figures (3.470751648e-01 N).
- **template_utube_manometer** rejected by 3 of 3:
  - mec-1 (round 1: reject) [physics or scenario implausible]: Instance 1 connects a mercury manometer to a pipe of liquid oxygen (1141 kg/m^3, about -183 C): mercury freezes at -39 C (and every other manometer liquid in the table freezes well above that), so a liquid mercury/LOX interface 0.15 m below the pipe cannot exist; cryogenic pipe fluids must be excluded from PIPE_FLUIDS for this template. The hydrostatic balance and arithmetic (50.162 kPa) are correct.
  - mec-2 (round 1: reject) [physics or scenario implausible]: Instance 1 connects a mercury U-tube directly to a pipe of liquid oxygen (1141 kg/m^3, about -183 C) with a LOX column of h1 = 0.15 m resting on the mercury: at that temperature mercury (freezes at -38.8 C) is solid, and in a room-temperature leg the LOX would boil, so the hydrostatic balance giving 50.162 kPa describes an impossible setup. The pipe-fluid table still lets cryogenic liquids be paired with any manometer liquid; exclude them or pair fluids by compatible temperature.
  - mec-3 (round 1: reject) [physics or scenario implausible]: Physically implausible pairing: instance 1 puts a mercury manometer on a pipe of liquid oxygen (1141 kg/m^3, about -183 C), but mercury freezes at -38.8 C, so no liquid mercury/LOX interface can exist and the 50.162 kPa reading describes an impossible device. PIPE_FLUIDS still includes cryogenic liquids that are paired freely with room-temperature manometer liquids; drop them or screen those pairs out (the same kind of fault as the molten-zinc fix).
