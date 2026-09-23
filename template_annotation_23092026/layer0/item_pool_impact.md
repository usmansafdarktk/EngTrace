# Layer 0 and screen pass 1 — cumulative item-pool impact

Measured 2026-09-24 by `tests/constants_integrity/c3_instance_dump.py --diff`: 150 templates × 300 seeds dumped from the pre-edit tree (a worktree at `b1445ad`) and from the working tree after the three closure rounds and the screen-pass-1 fixes, in two separate processes. A **question** change means the instance must be re-inferred; an **answer** change means an archived trace for that instance must be re-scored; a **solution**-only change moves gold text but not the graded values. Per-template reasons: `closure_fixes.md` (rounds 1–3) and `../screen/pass1_fixes.md` (the 24 flagged templates).

```
before: C:\Users\AYESHA~1.GUL\AppData\Local\Temp\claude\c--Users-ayesha-gull01-EngTrace\ee94e7a8-9f53-4f47-9626-9630fa10c39c\scratchpad\wt_head (n=300)
after : C:\Users\ayesha.gull01\EngTrace (n=300)

template_absorbing_chain_time_to_failure         question    5/300  answer   19/300  solution  300/300  errors 0->0
    seed 0 answer: 'The expected time until failure for a machine currently in Good condit' -> 'The expected time until failure for a machine currently in Good condit'
template_angle_of_twist                          question    0/300  answer  296/300  solution  296/300  errors 0->0
    seed 0 answer: 'The total angle of twist is 0.0579 radians, which is equivalent to 3.3' -> 'The total angle of twist is 0.0579 radians, which is equivalent to 3.3'
template_average_energy_mqam                     question   10/300  answer   10/300  solution  300/300  errors 0->0
    seed 0 answer: 'The average energy per symbol is 61.2.' -> 'The average energy per symbol is 61.2.'
template_basic_buoyant_force                     question  287/300  answer  287/300  solution  287/300  errors 0->0
    seed 0 answer: 'The buoyant force acting on the weight is 13349.448 N, which is equiva' -> 'The buoyant force acting on the pyramid is 8454.65 N, which is equival'
template_basic_stress_strain                     question    2/300  answer   12/300  solution   12/300  errors 0->0
    seed 47 answer: 'a) Normal Stress (sigma) = **18.538 ksi**\nb) Normal Strain (epsilon) =' -> 'a) Normal Stress (sigma) = **18.539 ksi**\nb) Normal Strain (epsilon) ='
template_batch_moles_vs_conversion               question    0/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'After reaching a conversion of $54.0 \\%$, the final number of moles in' -> 'After reaching a conversion of $54.0 \\%$, the final number of moles in'
template_batch_reactor_second_order              question   10/300  answer   38/300  solution  300/300  errors 0->0
    seed 0 answer: 'The required reaction time is 2.65 seconds.' -> 'The required reaction time is 2.65 seconds.'
template_beam_internal_moment                    question    5/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The internal bending moment at the section is 154.47 kN*m' -> 'The internal bending moment at the section is 154.474 kN*m'
template_beam_support_reactions                  question    9/300  answer   46/300  solution  155/300  errors 0->0
    seed 0 answer: 'The vertical reaction at support A is 41.71 kN' -> 'The vertical reaction at support A is 41.72 kN'
template_ber_estimation_mary                     question    0/300  answer   51/300  solution  134/300  errors 0->0
    seed 0 answer: 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.750 *' -> 'The estimated Bit Error Rate (BER) is expressed as:\nBER approx 0.750 *'
template_bpsk_energy_basis                       question    1/300  answer  245/300  solution  279/300  errors 0->0
    seed 0 answer: 'a) The energy per bit (Eb) is 280.3 mJ.\nb) The basis function (psi_1(t' -> 'a) The energy per bit (Eb) is 280.31 mJ.\nb) The basis function (psi_1('
template_cd_dc_system_analysis                   question    0/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The final continuous-time output signal is y_c(t) = 28.82 * cos(876*pi' -> 'The final continuous-time output signal is y_c(t) = 28.82 * cos(876*pi'
template_chart_pair_selection                    question    2/300  answer    2/300  solution    2/300  errors 0->0
    seed 116 answer: "The X-bar chart's upper control limit under the appropriate chart pair" -> "The X-bar chart's upper control limit under the appropriate chart pair"
template_composite_shafts_series                 question   92/300  answer   92/300  solution   92/300  errors 0->0
    seed 2 answer: 'The total angle of twist at the free end C is 0.78189 radians, or 44.7' -> 'The total angle of twist at the free end C is 0.48453 radians, or 27.7'
template_decimation_aliasing_analysis            question    0/300  answer    0/300  solution   31/300  errors 0->0
    seed 14 answer: 'a) The output signal is y[n] = cos((11*pi)/6*n).\nb) Aliasing **did** o' -> 'a) The output signal is y[n] = cos((11*pi)/6*n).\nb) Aliasing **did** o'
template_effective_stress_profile                question    2/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The vertical effective stress at point A is 118.1 kPa' -> 'The vertical effective stress at point A is 118.175 kPa'
template_epq_finite_production                   question    1/300  answer    1/300  solution  300/300  errors 0->0
    seed 0 answer: 'The optimal production run size is 3630 units' -> 'The optimal production run size is 3630 units'
template_floating_object_submersion_depth        question  122/300  answer  121/300  solution  122/300  errors 0->0
    seed 0 answer: 'The submersion depth of the cylinder is 1.4334 m.' -> 'The submersion depth of the cylinder is 0.7311 m.'
template_flow_system_molar_flow_rates            question    0/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The molar flow rates exiting the reactor are:\n- CH4(g) (F_A): 42.77 mo' -> 'The molar flow rates exiting the reactor are:\n- CH4(g) (F_A): 42.7720 '
template_force_method_continuous_beam            question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The vertical reaction at support B is 30.25 kN' -> 'The vertical reaction at support B is 30.25 kN'
template_gas_phase_concentration                 question  300/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The outlet concentrations are:\n- **$C_A$:** $0.0187$ mol/dm³\n- **$C_B$' -> 'The outlet concentrations are:\n- **$C_A$:** $0.0878$ mol/dm³\n- **$C_B$'
template_gas_viscosity_kinetic_theory            question    0/300  answer   14/300  solution  273/300  errors 0->0
    seed 0 answer: 'The estimated dynamic viscosity of Butane (C₄H₁₀) at 465 K is **2.094e' -> 'The estimated dynamic viscosity of Butane (C₄H₁₀) at 465 K is **2.094e'
template_gauss_law_symmetric                     question  300/300  answer  137/300  solution  300/300  errors 0->0
    seed 0 answer: 'The electric flux density is D = <0.000e+00, 0.000e+00, 2.000e-09> C/m' -> 'The electric flux density is D = <0.000e+00, 0.000e+00, 2.000e-09> C/m'
template_hydrostatic_pressure_at_depth           question    0/300  answer  258/300  solution  258/300  errors 0->0
    seed 0 answer: 'The absolute pressure at a depth of 19.43 m is **359.181 kPa**.' -> 'The absolute pressure at a depth of 19.43 m is **359.1815 kPa**.'
template_ideal_gas_volume                        question    0/300  answer  106/300  solution  297/300  errors 0->0
    seed 0 answer: 'The volume occupied by the gas is **58.86 liters**.' -> 'The volume occupied by the gas is **58.85 liters**.'
template_impulse_response_from_lccde             question    0/300  answer  296/300  solution  300/300  errors 0->0
    seed 0 answer: 'Substituting the coefficients back into the general form, the impulse ' -> 'Substituting the coefficients back into the general form, the impulse '
template_influence_line_max_reaction             question    1/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The maximum reaction at A is 216.80 kN' -> 'The maximum reaction at A is 216.8033 kN'
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
template_mmc_waiting_time                        question    0/300  answer   31/300  solution  300/300  errors 0->0
    seed 0 answer: 'The average time a customer waits in the queue is 5.09 minutes' -> 'The average time a customer waits in the queue is 5.09 minutes'
template_multi_segment_rod                       question    0/300  answer   95/300  solution   95/300  errors 0->0
    seed 0 answer: 'The total deformation of the rod is **2.3919 mm** (a net elongation).' -> 'The total deformation of the rod is **2.392 mm** (a net elongation).'
template_newsvendor_normal_demand                question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The profit-maximizing stocking quantity is 1505 units' -> 'The profit-maximizing stocking quantity is 1505 units'
template_pfr_volume_changing_rate                question   32/300  answer    0/300  solution    0/300  errors 0->0
    seed 3 answer: 'The required PFR volume is 7.3 liters.' -> 'The required PFR volume is 7.3 liters.'
template_phase_relations_degree_of_saturation    question    5/300  answer    5/300  solution    5/300  errors 0->0
    seed 43 answer: 'The degree of saturation is 89.2 %' -> 'The degree of saturation is 55.5 %'
template_phasor_addition                         question    0/300  answer  200/300  solution  223/300  errors 0->0
    seed 0 answer: 'The sum of the two signals is v_total(t) = 73.5 * cos(361*t - 146.39 d' -> 'The sum of the two signals is v_total(t) = 73.51 * cos(361*t - 146.39 '
template_pitzer_correlation_z                    question   13/300  answer   53/300  solution  300/300  errors 0->0
    seed 0 answer: 'The compressibility factor, Z, for p-Xylene at the given conditions is' -> 'The compressibility factor, Z, for p-Xylene at the given conditions is'
template_poissons_ratio                          question   62/300  answer   84/300  solution  300/300  errors 0->0
    seed 0 answer: 'a) The change in diameter is **-0.00068 in**.\nb) The final diameter is' -> 'a) The change in diameter is **-0.00182 in**.\nb) The final diameter is'
template_power_law_fluid_shear                   question   46/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'a) The shear stress on the fluid is **80.272 Pa**.\nb) The apparent vis' -> 'a) The shear stress on the fluid is **80.272 Pa**.\nb) The apparent vis'
template_primary_consolidation_settlement        question    0/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The primary consolidation settlement is 145 mm' -> 'The primary consolidation settlement is 144.8 mm'
template_qr_policy_one_iteration                 question   14/300  answer   14/300  solution   14/300  errors 0->0
    seed 3 answer: 'The updated lot size after one iteration is 1510 units' -> 'The updated lot size after one iteration is 1792 units'
template_relative_density_of_sand                question   10/300  answer   10/300  solution   10/300  errors 0->0
    seed 23 answer: 'The relative density is 30.23 %' -> 'The relative density is 33.29 %'
template_reorder_point_lead_time                 question    4/300  answer    4/300  solution  138/300  errors 0->0
    seed 0 answer: 'The reorder point is 308 units' -> 'The reorder point is 308 units'
template_reynolds_number_flow_regime             question   19/300  answer    5/300  solution  300/300  errors 0->0
    seed 0 answer: 'a) The Reynolds number is approximately **920,610**.\nb) The flow regim' -> 'a) The Reynolds number is approximately **920,610**.\nb) The flow regim'
template_rotating_unbalance                      question    0/300  answer    0/300  solution  300/300  errors 0->0
    seed 0 answer: 'The steady-state amplitude of vibration is **4.686 mm**.' -> 'The steady-state amplitude of vibration is **4.686 mm**.'
template_sensible_heat_constant_cp               question    2/300  answer    2/300  solution    2/300  errors 0->0
    seed 117 answer: 'The heat required is **47.46 kJ**.' -> 'The heat required is **113.36 kJ**.'
template_sensible_heat_temp_dependent_cp         question  169/300  answer    0/300  solution  180/300  errors 0->0
    seed 3 answer: 'The heat required is **873.72 kJ**.' -> 'The heat required is **873.72 kJ**.'
template_server_configuration_selection          question   78/300  answer   13/300  solution   13/300  errors 0->0
    seed 0 answer: 'The better option achieves an average time in the system of 3.12 minut' -> 'The better option achieves an average time in the system of 3.12 minut'
template_shear_stress_torsion                    question    0/300  answer  213/300  solution  246/300  errors 0->0
    seed 0 answer: 'The maximum shearing stress in the shaft is 563.453 MPa.' -> 'The maximum shearing stress in the shaft is 563.52 MPa.'
template_signal_operations                       question  300/300  answer   54/300  solution  300/300  errors 0->0
    seed 0 answer: 'The resulting sequence is z[n] = {5, -1, *2*, 5, 6, -2, -9}' -> 'The resulting sequence is z[n] = {5, -1, *2*, 5, 6, -2, -9}'
template_slope_deflection_end_moment             question    9/300  answer   59/300  solution  300/300  errors 0->0
    seed 0 answer: 'The end moment at support A is -27.67 kN*m' -> 'The end moment at support A is -27.68 kN*m'
template_statically_indeterminate                question    2/300  answer  119/300  solution  194/300  errors 0->0
    seed 0 answer: 'a) Reaction Forces: R_A = **38.264 kips** and R_C = **135.736 kips**.\n' -> 'a) Reaction Forces: R_A = **38.264 kips** and R_C = **135.736 kips**.\n'
template_takt_time_line_efficiency               question   12/300  answer   12/300  solution   12/300  errors 0->0
    seed 68 answer: 'The line efficiency at the theoretical minimum number of stations is 9' -> 'The line efficiency at the theoretical minimum number of stations is 9'
template_terzaghi_strip_footing_bearing          question    1/300  answer  300/300  solution  300/300  errors 0->0
    seed 0 answer: 'The allowable bearing capacity is 89.5 kPa' -> 'The allowable bearing capacity is 89.47 kPa'
template_time_to_phasor                          question    0/300  answer   46/300  solution   46/300  errors 0->0
    seed 4 answer: 'The phasor representation is 39.23 < 145.67 degrees, which is equivale' -> 'The phasor representation is 39.23 < 145.67 degrees, which is equivale'
template_truss_method_of_joints                  question    2/300  answer    2/300  solution  300/300  errors 0->0
    seed 0 answer: 'The force in member AB is 28.23 kN' -> 'The force in member AB is 28.23 kN'
template_truss_method_of_sections                question   15/300  answer   15/300  solution   15/300  errors 0->0
    seed 16 answer: 'The force in member AC is 25.88 kN (tension)' -> 'The force in member BC is 47.48 kN (tension)'
template_undamped_natural_frequency_translational question    0/300  answer   27/300  solution   27/300  errors 0->0
    seed 15 answer: 'The undamped natural frequency is 2.216 rad/s or 0.353 Hz.\nThe natural' -> 'The undamped natural frequency is 2.216 rad/s or 0.353 Hz.\nThe natural'
template_undamped_response_initial_conditions    question    0/300  answer    1/300  solution    1/300  errors 0->0
    seed 151 answer: 'The final equation of motion is:\nx(t) = 0.085*cos(23.0379*t) - 0.1281*' -> 'The final equation of motion is:\nx(t) = 0.085*cos(23.0379*t) - 0.128*s'
template_upward_seepage_quick_condition          question   14/300  answer   14/300  solution   14/300  errors 0->0
    seed 22 answer: 'The factor of safety against the quick condition is 3.000' -> 'The factor of safety against the quick condition is 1.853'
template_utube_manometer                         question    1/300  answer    1/300  solution  300/300  errors 0->0
    seed 0 answer: 'The gauge pressure in the pipe is 19.891 kPa.' -> 'The gauge pressure in the pipe is 19.891 kPa.'
template_virtual_work_truss_deflection           question   21/300  answer   21/300  solution   21/300  errors 0->0
    seed 37 answer: 'The vertical deflection of joint C is 8.797 mm' -> 'The vertical deflection of joint C is 7.979 mm'
template_wave_equation_interpretation            question   26/300  answer  185/300  solution  300/300  errors 0->0
    seed 0 answer: '- Amplitude: 108 V/m\n   - Direction of Propagation: the positive z-dir' -> '- Amplitude: 108 V/m\n   - Direction of Propagation: the positive z-dir'
template_wave_parameters_basic                   question    0/300  answer  284/300  solution  284/300  errors 0->0
    seed 0 answer: '- Frequency (f): 1095.73 MHz\n- Angular Frequency (omega): 6.88e+09 rad' -> '- Frequency (f): 1094.44 MHz\n- Angular Frequency (omega): 6.85e+09 rad'
template_work_isothermal_virial                  question  161/300  answer  285/300  solution  300/300  errors 0->0
    seed 0 answer: 'The required work of compression is approximately **13456.0 J/mol**.' -> 'The required work of compression is approximately **13475.0 J/mol**.'
template_xbar_r_control_limits                   question    1/300  answer    1/300  solution    1/300  errors 0->0
    seed 129 answer: 'The upper control limit of the X-bar chart is 13.292 mm' -> 'The upper control limit of the X-bar chart is 25.853 mm'

68 templates moved; totals {'q': 2203, 'ans': 6993, 'sol': 13692, 'err_before': 0, 'err_after': 0}
```
