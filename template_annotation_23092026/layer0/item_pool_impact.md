# Layer 0 — item-pool impact of the closure fixes

Measured 2026-09-23 by `tests/constants_integrity/c3_instance_dump.py --diff`: 150 templates × 300 seeds dumped from the pre-edit tree (a worktree at `b1445ad`) and from the edited working tree, in two separate processes (the in-process trap the instrument's docstring describes). A **question** change means the instance must be re-inferred; an **answer** change means an archived trace for that instance must be re-scored; a **solution**-only change moves gold text but not the graded values. Every edit is explained per template in `closure_fixes.md`.

```
before: C:\Users\AYESHA~1.GUL\AppData\Local\Temp\claude\c--Users-ayesha-gull01-EngTrace\ee94e7a8-9f53-4f47-9626-9630fa10c39c\scratchpad\wt_head (n=300)
after : C:\Users\ayesha.gull01\EngTrace (n=300)

template_absorbing_chain_time_to_failure         question    5/300  answer   19/300  solution  300/300  errors 0->0
    seed 0 answer: 'The expected time until failure for a machine currently in Good condit' -> 'The expected time until failure for a machine currently in Good condit'
template_angle_of_twist                          question    0/300  answer  296/300  solution  296/300  errors 0->0
    seed 0 answer: 'The total angle of twist is 0.0579 radians, which is equivalent to 3.3' -> 'The total angle of twist is 0.0579 radians, which is equivalent to 3.3'
template_basic_stress_strain                     question    2/300  answer   12/300  solution   12/300  errors 0->0
    seed 47 answer: 'a) Normal Stress (sigma) = **18.538 ksi**\nb) Normal Strain (epsilon) =' -> 'a) Normal Stress (sigma) = **18.539 ksi**\nb) Normal Strain (epsilon) ='
template_batch_moles_vs_conversion               question    0/300  answer   11/300  solution   11/300  errors 0->0
    seed 15 answer: 'After reaching a conversion of $49.0 \\%$, the final number of moles in' -> 'After reaching a conversion of $49.0 \\%$, the final number of moles in'
template_beam_internal_moment                    question   37/300  answer   37/300  solution   37/300  errors 0->0
    seed 10 answer: 'The internal bending moment at the section is 76.01 kN*m' -> 'The internal bending moment at the section is 68.61 kN*m'
template_beam_support_reactions                  question    9/300  answer   46/300  solution  155/300  errors 0->0
    seed 0 answer: 'The vertical reaction at support A is 41.71 kN' -> 'The vertical reaction at support A is 41.72 kN'
template_ber_estimation_mary                     question    0/300  answer   51/300  solution  134/300  errors 0->0
    seed 0 answer: 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.750 *' -> 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.750 *'
template_bpsk_energy_basis                       question    1/300  answer  245/300  solution  279/300  errors 0->0
    seed 0 answer: 'a) The energy per bit (Eb) is 280.3 mJ.\nb) The basis function (psi_1(t' -> 'a) The energy per bit (Eb) is 280.31 mJ.\nb) The basis function (psi_1('
template_chart_pair_selection                    question    2/300  answer    2/300  solution    2/300  errors 0->0
    seed 116 answer: "The X-bar chart's upper control limit under the appropriate chart pair" -> "The X-bar chart's upper control limit under the appropriate chart pair"
template_decimation_aliasing_analysis            question    0/300  answer    0/300  solution   31/300  errors 0->0
    seed 14 answer: 'a) The output signal is y[n] = cos((11*pi)/6*n).\nb) Aliasing **did** o' -> 'a) The output signal is y[n] = cos((11*pi)/6*n).\nb) Aliasing **did** o'
template_effective_stress_profile                question    6/300  answer   76/300  solution  300/300  errors 0->0
    seed 0 answer: 'The vertical effective stress at point A is 118.1 kPa' -> 'The vertical effective stress at point A is 118.2 kPa'
template_epq_finite_production                   question    1/300  answer    1/300  solution  300/300  errors 0->0
    seed 0 answer: 'The optimal production run size is 3630 units' -> 'The optimal production run size is 3630 units'
template_force_method_continuous_beam            question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The vertical reaction at support B is 30.25 kN' -> 'The vertical reaction at support B is 30.25 kN'
template_gas_viscosity_kinetic_theory            question    0/300  answer   14/300  solution  273/300  errors 0->0
    seed 0 answer: 'The estimated dynamic viscosity of Butane (C₄H₁₀) at 465 K is **2.094e' -> 'The estimated dynamic viscosity of Butane (C₄H₁₀) at 465 K is **2.094e'
template_hydrostatic_pressure_at_depth           question    0/300  answer  258/300  solution  258/300  errors 0->0
    seed 0 answer: 'The absolute pressure at a depth of 19.43 m is **359.181 kPa**.' -> 'The absolute pressure at a depth of 19.43 m is **359.1815 kPa**.'
template_lorentz_force                           question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The magnetic force vector is F_m = (-4.020e-04 x_hat + 4.952e-04 y_hat' -> 'The magnetic force vector is F_m = (-4.020e-04 x_hat + 4.952e-04 y_hat'
template_lowpass_equivalent_bandpass             question    6/300  answer   81/300  solution  300/300  errors 0->0
    seed 0 answer: 'In-Phase Component: **g_I(t) = 29.5**\nQuadrature Component: **g_Q(t) =' -> 'In-Phase Component: **g_I(t) = 29.5**\nQuadrature Component: **g_Q(t) ='
template_max_hump_height_no_choking              question    0/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The maximum hump height without choking is 0.607 m' -> 'The maximum hump height without choking is 0.6074 m'
template_mm1_time_in_system                      question   13/300  answer   10/300  solution   13/300  errors 0->0
    seed 46 answer: 'The average time a customer spends in the system is 3.8 minutes' -> 'The average time a customer spends in the system is 6.7 minutes'
template_mm1k_finite_capacity                    question   16/300  answer   16/300  solution  300/300  errors 0->0
    seed 0 answer: 'The average time in the system for cars that join the lane is 4.79 min' -> 'The average time in the system for cars that join the lane is 4.79 min'
template_mmc_waiting_time                        question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The average time a customer waits in the queue is 5.09 minutes' -> 'The average time a customer waits in the queue is 5.09 minutes'
template_multi_segment_rod                       question    0/300  answer   95/300  solution   95/300  errors 0->0
    seed 0 answer: 'The total deformation of the rod is **2.3919 mm** (a net elongation).' -> 'The total deformation of the rod is **2.392 mm** (a net elongation).'
template_newsvendor_normal_demand                question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The profit-maximizing stocking quantity is 1505 units' -> 'The profit-maximizing stocking quantity is 1505 units'
template_phase_relations_degree_of_saturation    question    5/300  answer    5/300  solution    5/300  errors 0->0
    seed 43 answer: 'The degree of saturation is 89.2 %' -> 'The degree of saturation is 55.5 %'
template_phasor_addition                         question    0/300  answer  200/300  solution  223/300  errors 0->0
    seed 0 answer: 'The sum of the two signals is v_total(t) = 73.5 * cos(361*t - 146.39 d' -> 'The sum of the two signals is v_total(t) = 73.51 * cos(361*t - 146.39 '
template_pitzer_correlation_z                    question    0/300  answer   45/300  solution  300/300  errors 0->0
    seed 0 answer: 'The compressibility factor, Z, for p-Xylene at the given conditions is' -> 'The compressibility factor, Z, for p-Xylene at the given conditions is'
template_poissons_ratio                          question    4/300  answer   36/300  solution  300/300  errors 0->0
    seed 0 answer: 'a) The change in diameter is **-0.00068 in**.\nb) The final diameter is' -> 'a) The change in diameter is **-0.00068 in**.\nb) The final diameter is'
template_qr_policy_one_iteration                 question   14/300  answer   14/300  solution   14/300  errors 0->0
    seed 3 answer: 'The updated lot size after one iteration is 1510 units' -> 'The updated lot size after one iteration is 1792 units'
template_reorder_point_lead_time                 question    4/300  answer    4/300  solution  138/300  errors 0->0
    seed 0 answer: 'The reorder point is 308 units' -> 'The reorder point is 308 units'
template_sensible_heat_constant_cp               question    2/300  answer    2/300  solution    2/300  errors 0->0
    seed 117 answer: 'The heat required is **47.46 kJ**.' -> 'The heat required is **113.36 kJ**.'
template_sensible_heat_temp_dependent_cp         question  169/300  answer    0/300  solution  180/300  errors 0->0
    seed 3 answer: 'The heat required is **873.72 kJ**.' -> 'The heat required is **873.72 kJ**.'
template_shear_stress_torsion                    question    0/300  answer  213/300  solution  246/300  errors 0->0
    seed 0 answer: 'The maximum shearing stress in the shaft is 563.453 MPa.' -> 'The maximum shearing stress in the shaft is 563.52 MPa.'
template_slope_deflection_end_moment             question    9/300  answer   59/300  solution  300/300  errors 0->0
    seed 0 answer: 'The end moment at support A is -27.67 kN*m' -> 'The end moment at support A is -27.68 kN*m'
template_statically_indeterminate                question    2/300  answer  119/300  solution  194/300  errors 0->0
    seed 0 answer: 'a) Reaction Forces: R_A = **38.264 kips** and R_C = **135.736 kips**.\n' -> 'a) Reaction Forces: R_A = **38.264 kips** and R_C = **135.736 kips**.\n'
template_terzaghi_strip_footing_bearing          question    7/300  answer    7/300  solution  156/300  errors 0->0
    seed 1 answer: 'The allowable bearing capacity is 540.7 kPa' -> 'The allowable bearing capacity is 540.7 kPa'
template_time_to_phasor                          question    0/300  answer   46/300  solution   46/300  errors 0->0
    seed 4 answer: 'The phasor representation is 39.23 < 145.67 degrees, which is equivale' -> 'The phasor representation is 39.23 < 145.67 degrees, which is equivale'
template_truss_method_of_joints                  question    2/300  answer    2/300  solution    2/300  errors 0->0
    seed 43 answer: 'The force in member AB is 18.12 kN' -> 'The force in member AB is 47.18 kN'
template_truss_method_of_sections                question   15/300  answer   15/300  solution   15/300  errors 0->0
    seed 16 answer: 'The force in member AC is 25.88 kN (tension)' -> 'The force in member BC is 47.48 kN (tension)'
template_undamped_response_initial_conditions    question    0/300  answer    1/300  solution    1/300  errors 0->0
    seed 151 answer: 'The final equation of motion is:\nx(t) = 0.085*cos(23.0379*t) - 0.1281*' -> 'The final equation of motion is:\nx(t) = 0.085*cos(23.0379*t) - 0.128*s'
template_utube_manometer                         question    1/300  answer    1/300  solution  300/300  errors 0->0
    seed 0 answer: 'The gauge pressure in the pipe is 19.891 kPa.' -> 'The gauge pressure in the pipe is 19.891 kPa.'
template_virtual_work_truss_deflection           question   21/300  answer   21/300  solution   21/300  errors 0->0
    seed 37 answer: 'The vertical deflection of joint C is 8.797 mm' -> 'The vertical deflection of joint C is 7.979 mm'
template_wave_parameters_basic                   question    0/300  answer  284/300  solution  284/300  errors 0->0
    seed 0 answer: '- Frequency (f): 1095.73 MHz\n- Angular Frequency (omega): 6.88e+09 rad' -> '- Frequency (f): 1094.44 MHz\n- Angular Frequency (omega): 6.85e+09 rad'
template_work_isothermal_virial                  question    0/300  answer  143/300  solution  300/300  errors 0->0
    seed 0 answer: 'The required work of compression is approximately **13456.0 J/mol**.' -> 'The required work of compression is approximately **13457.0 J/mol**.'
template_xbar_r_control_limits                   question    1/300  answer    1/300  solution    1/300  errors 0->0
    seed 129 answer: 'The upper control limit of the X-bar chart is 13.292 mm' -> 'The upper control limit of the X-bar chart is 25.853 mm'

44 templates moved; totals {'q': 354, 'ans': 2788, 'sol': 7624, 'err_before': 0, 'err_after': 0}
```
