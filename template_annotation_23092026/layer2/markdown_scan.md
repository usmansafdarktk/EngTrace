# Markdown rendering scan

Generated 2026-09-26T13:09:52+00:00 by `markdown_scan.py` at git `85ccf8d028`, 200 seeds per template over 150 templates. The review app renders questions and solutions as Markdown; the models read the raw text. A lossy construct removes characters from what the expert sees.

| | Templates |
|---|---:|
| any construct in the question | 43 |
| lossy construct in the question | 31 |
| lossy construct in the solution | 63 |
| lossy in either | 65 |

| Template | Branch | Construct | Seeds | Example |
|---|---|---|---:|---|
| `template_batch_moles_vs_conversion` | chemical | question: dollar pair (inline math) | 200 | seed 0: $CH4(g) + H2O(g) → CO(g) + 3H2(g)$ |
| `template_batch_moles_vs_conversion` | chemical | question: strong (bold) | 200 | seed 0: **Reaction:** |
| `template_batch_moles_vs_conversion` | chemical | solution: backslash escape (backslash dropped) | 200 | seed 0: rsion of $54.0 \%$, the final nu |
| `template_batch_moles_vs_conversion` | chemical | solution: dollar pair (inline math) | 200 | seed 0: $CH4(g) + H2O(g) → CO(g) + 3H2(g)$ - Initial Moles:   - $N_{CH4(g),0}  |
| `template_flow_system_molar_flow_rates` | chemical | question: strong (bold) | 200 | seed 0: **Reaction:** |
| `template_flow_system_molar_flow_rates` | chemical | solution: emphasis (italic) | 200 | seed 0: *Using the Direct Method:* |
| `template_gas_phase_concentration` | chemical | question: bullet list | 200 | seed 0: - |
| `template_gas_phase_concentration` | chemical | question: dollar pair (inline math) | 200 | seed 0: $CO(g) + 3H2(g) → CH4(g) + H2O(g)$, i.e. Carbon Monoxide + 3 Hydrogen  |
| `template_gas_phase_concentration` | chemical | question: strong (bold) | 200 | seed 0: **Reaction:** |
| `template_gas_phase_concentration` | chemical | solution: backslash escape (backslash dropped) | 200 | seed 0: silon = y_{A0} \, \delta = 0.178 |
| `template_gas_phase_concentration` | chemical | solution: dollar pair (inline math) | 200 | seed 0: $C_{A0} = 0.202$ mol/dm³ - Conversion: $X_A = 0.67$ - Molar Feed Ratio |
| `template_gas_viscosity_kinetic_theory` | chemical | question: bullet list | 200 | seed 0: - |
| `template_heat_of_reaction_formation` | chemical | question: bullet list | 200 | seed 0: - |
| `template_kinematic_viscosity` | chemical | question: strong (bold) | 200 | seed 0: **Stokes (St)** |
| `template_levenspiel_plot_interpretation` | chemical | question: heading | 200 | seed 0: - |
| `template_limiting_reactant` | chemical | question: strong (bold) | 200 | seed 0: **Reaction:** |
| `template_power_law_fluid_shear` | chemical | question: bullet list | 200 | seed 0: - |
| `template_reynolds_number_flow_regime` | chemical | question: bullet list | 200 | seed 0: - |
| `template_sensible_heat_temp_dependent_cp` | chemical | question: emphasis (italic) | 177 | seed 0: *T + 1.510e+04* |
| `template_sensible_heat_temp_dependent_cp` | chemical | solution: emphasis (italic) | 184 | seed 0: *(668.67)^2 - (1.510e+04/668.67) = 2227.0181I(383.0* |
| `template_vdw_solve_for_volume` | chemical | solution: emphasis (italic) | 200 | seed 0: *V^2 + c1* |
| `template_work_isothermal_virial` | chemical | solution: emphasis (italic) | 200 | seed 0: *0.04823* |
| `template_beam_deflection_formula` | civil | solution: emphasis (italic) | 200 | seed 0: *L^3 / (48* |
| `template_cantilever_double_integration` | civil | solution: emphasis (italic) | 200 | seed 0: *I* |
| `template_force_method_continuous_beam` | civil | solution: emphasis (italic) | 200 | seed 0: *b* |
| `template_hydraulic_jump_energy_loss` | civil | solution: emphasis (italic) | 200 | seed 0: *y1* |
| `template_linear_reservoir_routing_step` | civil | solution: emphasis (italic) | 200 | seed 0: *dt/2 - O_j* |
| `template_manning_trapezoidal_velocity` | civil | solution: emphasis (italic) | 200 | seed 0: *y* |
| `template_max_hump_height_no_choking` | civil | solution: emphasis (italic) | 200 | seed 0: *g* |
| `template_normal_depth_iteration` | civil | question: emphasis (italic) | 200 | seed 0: *n/S^(1/2), define g(y) = A* |
| `template_normal_depth_iteration` | civil | solution: emphasis (italic) | 200 | seed 0: *R^(2/3) = Q* |
| `template_rational_method_peak_flow` | civil | question: emphasis (italic) | 200 | seed 0: *i* |
| `template_rational_method_peak_flow` | civil | solution: emphasis (italic) | 93 | seed 0: *A1 + C2* |
| `template_slope_deflection_end_moment` | civil | question: emphasis (italic) | 200 | seed 0: *L^2/12 (near end) and +w* |
| `template_slope_deflection_end_moment` | civil | solution: emphasis (italic) | 200 | seed 0: *L1^2/12 = -18 * (5.1)^2 / 12 = -39.0150 kN* |
| `template_terzaghi_strip_footing_bearing` | civil | question: emphasis (italic) | 200 | seed 0: *Nc + q* |
| `template_terzaghi_strip_footing_bearing` | civil | solution: emphasis (italic) | 200 | seed 0: *N'c + q* |
| `template_time_rate_of_consolidation` | civil | question: emphasis (italic) | 200 | seed 0: *t / H_dr^2, where the specimen and the field layer* |
| `template_virtual_work_truss_deflection` | civil | solution: emphasis (italic) | 200 | seed 0: *N* |
| `template_aliased_frequency_identification` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_autocorrelation_rect_pulse` | electrical | solution: emphasis (italic) | 200 | seed 0: **R_g(tau) = 64*(8 - \|tau\|)* |
| `template_ber_estimation_mary` | electrical | solution: emphasis (italic) | 200 | seed 0: *Q(sqrt(3* |
| `template_bpsk_energy_basis` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_bpsk_energy_basis` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_cd_dc_system_analysis` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_cd_dc_system_analysis` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_continuous_to_discrete_conversion` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_continuous_to_discrete_conversion` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_decimation_aliasing_analysis` | electrical | question: emphasis (italic) | 131 | seed 1: *pi)/12* |
| `template_decimation_aliasing_analysis` | electrical | solution: emphasis (italic) | 188 | seed 0: *pi)/8So, the expression is y[n] = cos((5* |
| `template_euclidean_distance_binary` | electrical | solution: emphasis (italic) | 107 | seed 1: *sqrt(Eb))d = \|\|s1 - s2\|\| = 2* |
| `template_finite_convolution` | electrical | question: emphasis (italic) | 200 | seed 0: *2* |
| `template_finite_convolution` | electrical | solution: emphasis (italic) | 200 | seed 0: *0* |
| `template_ft_esd_rect_pulse` | electrical | question: emphasis (italic) | 200 | seed 0: *x) / (pi* |
| `template_ft_esd_rect_pulse` | electrical | solution: emphasis (italic) | 200 | seed 0: *x) / (pi* |
| `template_gauss_law_symmetric` | electrical | solution: emphasis (italic) | 93 | seed 0: *A + D_z* |
| `template_impulse_response_from_lccde` | electrical | solution: emphasis (italic) | 200 | seed 0: *delta[0] ... (other delta terms are 0)h[0] + 3* |
| `template_lowpass_equivalent_bandpass` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_lowpass_equivalent_bandpass` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_mean_variance` | electrical | solution: emphasis (italic) | 200 | seed 0: *(0.173) + (3)* |
| `template_nyquist_rate_determination` | electrical | question: emphasis (italic) | 200 | seed 0: *pi* |
| `template_nyquist_rate_determination` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_phasor_addition` | electrical | question: emphasis (italic) | 200 | seed 0: *t - 28.6 deg)v2(t) = 40.32 * sin(361* |
| `template_phasor_addition` | electrical | solution: emphasis (italic) | 200 | seed 0: *t - 28.6 deg)v2(t) = 40.32 * sin(361* |
| `template_signal_energy_power` | electrical | solution: emphasis (italic) | 126 | seed 0: *t))^2 dtEg = 225 * integral from 0 to inf of exp(-* |
| `template_signal_operations` | electrical | question: emphasis (italic) | 200 | seed 0: *2* |
| `template_signal_operations` | electrical | solution: emphasis (italic) | 200 | seed 0: *2* |
| `template_standing_wave_formation` | electrical | question: emphasis (italic) | 200 | seed 0: *t - 4.56* |
| `template_standing_wave_formation` | electrical | solution: emphasis (italic) | 200 | seed 0: *t - 4.56* |
| `template_superposition_electric_field` | electrical | question: bullet list | 200 | seed 0: - |
| `template_superposition_electric_field` | electrical | solution: emphasis (italic) | 200 | seed 0: *pi* |
| `template_system_properties_memory_causality` | electrical | solution: emphasis (italic) | 200 | seed 0: *only* |
| `template_system_property_linearity` | electrical | solution: emphasis (italic) | 200 | seed 0: *x[n]} = a* |
| `template_wave_equation_interpretation` | electrical | question: emphasis (italic) | 200 | seed 0: *t - 8.89* |
| `template_wave_equation_interpretation` | electrical | solution: emphasis (italic) | 200 | seed 0: *t - 8.89* |
| `template_absorbing_chain_time_to_failure` | industrial | solution: emphasis (italic) | 200 | seed 0: *muG + 0.24* |
| `template_aoq_ati_rectifying` | industrial | solution: emphasis (italic) | 200 | seed 0: *p* |
| `template_arl_beta_mean_shift` | industrial | solution: emphasis (italic) | 200 | seed 0: *sqrt(n).n = 9: z = 3 - 0.64* |
| `template_basic_eoq` | industrial | question: dollar pair (inline math) | 200 | seed 0: $148 to place and process, the item costs $6.18 per unit, and the annu |
| `template_basic_eoq` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $148; unit cost c = $6.18; holding rate i = 0.16 per year. |
| `template_basic_eoq` | industrial | solution: emphasis (italic) | 200 | seed 0: *K* |
| `template_chase_vs_level_aggregate` | industrial | question: dollar pair (inline math) | 200 | seed 0: $1405 per worker, firing costs cF = $2347 per worker, and ending inven |
| `template_chase_vs_level_aggregate` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $1405; cF = $2347; h = $4/unit-month. |
| `template_chase_vs_level_aggregate` | industrial | solution: emphasis (italic) | 200 | seed 0: *2 + 2347* |
| `template_cp_cpk_from_specs` | industrial | solution: emphasis (italic) | 200 | seed 0: *sigma-hat) = (16.672 - 16.407) / (3 * 0.05746) = 1* |
| `template_epq_finite_production` | industrial | question: dollar pair (inline math) | 200 | seed 0: $87 to set up, the item is valued at $3.65 per unit, and the annual ho |
| `template_epq_finite_production` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $87; unit cost c = $3.65; holding rate i = 0.34 per year. |
| `template_epq_finite_production` | industrial | solution: emphasis (italic) | 200 | seed 0: *K* |
| `template_line_balancing_heuristic` | industrial | question: emphasis (italic) | 200 | seed 0: *CT - total work) / (n* |
| `template_line_balancing_heuristic` | industrial | solution: emphasis (italic) | 200 | seed 0: *CT - sum(t)) / (n* |
| `template_mm1k_finite_capacity` | industrial | solution: emphasis (italic) | 200 | seed 0: *1.1463 + 2* |
| `template_newsvendor_normal_demand` | industrial | question: dollar pair (inline math) | 200 | seed 0: $23.55, sells for $112.15 during the season, and any unsold unit is sa |
| `template_newsvendor_normal_demand` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $23.55; price p = $112.15; salvage s = $12.79; demand ~ Normal(mu = 11 |
| `template_p_chart_limits_floor` | industrial | question: emphasis (italic) | 200 | seed 0: *n), the standard error sqrt(p-bar* |
| `template_p_chart_limits_floor` | industrial | solution: emphasis (italic) | 200 | seed 0: *se = 0.1168 + 3 * 0.0203 = 0.1777;  computed LCL =* |
| `template_qr_policy_one_iteration` | industrial | question: dollar pair (inline math) | 200 | seed 0: $172 to place, the item costs c = $30.80 per unit, and the annual hold |
| `template_qr_policy_one_iteration` | industrial | question: emphasis (italic) | 200 | seed 0: *c. Each unit of unmet demand is backordered at a s* |
| `template_qr_policy_one_iteration` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $172; c = $30.80; i = 0.26; stockout penalty p = $15.87/unit; lead-tim |
| `template_qr_policy_one_iteration` | industrial | solution: emphasis (italic) | 200 | seed 0: *K* |
| `template_quantity_discount_all_units` | industrial | question: dollar pair (inline math) | 200 | seed 0: $115 to place, and the annual holding-cost rate is 0.28 (a fraction of |
| `template_quantity_discount_all_units` | industrial | solution: dollar pair (inline math) | 200 | seed 0: $115; holding rate i = 0.28 per year; all-units prices c0 = $17.02 (1  |
| `template_quantity_discount_all_units` | industrial | solution: emphasis (italic) | 200 | seed 0: *c0 = 0.28* |
| `template_reorder_point_lead_time` | industrial | question: emphasis (italic) | 200 | seed 0: *T for the largest integer k with k* |
| `template_reorder_point_lead_time` | industrial | solution: emphasis (italic) | 86 | seed 0: *tau - k* |
| `template_sigma_reduction_for_cpk` | industrial | question: emphasis (italic) | 200 | seed 0: *sigma), Cpu = (USL - mu)/(3* |
| `template_sigma_reduction_for_cpk` | industrial | solution: emphasis (italic) | 200 | seed 0: *sigma) = (85 - 54) / (6 * 3.71) = 1.39;  Cpu = (US* |
| `template_fluid_particle_acceleration` | mechanical | solution: emphasis (italic) | 200 | seed 0: *(19.6) + (52.16)* |
| `template_logarithmic_decrement` | mechanical | solution: emphasis (italic) | 200 | seed 0: *pi)^2 + delta^2)zeta = 0.3648 / sqrt((2* |
| `template_multi_segment_rod` | mechanical | question: bullet list | 200 | seed 0: - |
| `template_statically_indeterminate` | mechanical | solution: emphasis (italic) | 200 | seed 0: *R_A = 14900.0 - 59.6* |
| `template_statically_indeterminate_shaft` | mechanical | solution: emphasis (italic) | 200 | seed 0: *G) = (T_B * L_BC) / (J* |
| `template_system_properties` | mechanical | question: ordered list | 200 | seed 0: . |
| `template_undamped_response_initial_conditions` | mechanical | solution: emphasis (italic) | 200 | seed 0: *sin(21.7438* |
| `template_vibration_isolator_design` | mechanical | solution: emphasis (italic) | 200 | seed 0: *zeta* |
| `template_vibration_transmissibility` | mechanical | question: ordered list | 200 | seed 0: . |
| `template_vibration_transmissibility` | mechanical | solution: emphasis (italic) | 200 | seed 0: *zeta* |
| `template_volumetric_flow_rate` | mechanical | solution: emphasis (italic) | 200 | seed 0: *U_max/H) * [y^2/2] from 0 to HQ = (W* |
