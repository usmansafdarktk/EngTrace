# Exact display-tie census

Generated 2026-09-23T16:31:49+00:00 by `tie_census.py` at git `9d79d0df8d`, 500 seeds per template, 206.5 s. A tie is an emitted line whose exact value, recomputed in Decimal from the printed operands, sits exactly half-way at the printed precision; such an instance has no gold value a decimal reader and a binary reader agree on (D-016). Only lines T1 can parse are covered.

| | Templates |
|---|---:|
| censused | 150 |
| no tie in 500 instances | 137 |
| at least one tied instance | 13 |

| Template | Branch | Tied instances | Rate | Where |
|---|---|---:|---:|---|
| `template_server_configuration_selection` | industrial | 77 | 15.4% | `W1` x54; `W2` x24 |
| `template_manning_trapezoidal_velocity` | civil | 9 | 1.8% | `V` x5; `Q` x4 |
| `template_linear_reservoir_routing_step` | civil | 7 | 1.4% | `O2` x5; `O3` x3 |
| `template_sigma_reduction_for_cpk` | industrial | 7 | 1.4% | `reduction` x7 |
| `template_borrow_pit_fill_volume` | civil | 5 | 1.0% | `gamma_d,borrow` x4; `V_borrow` x1 |
| `template_rational_method_peak_flow` | civil | 4 | 0.8% | `C_w` x4 |
| `template_poissons_ratio` | mechanical | 4 | 0.8% | `epsilon_lateral` x4 |
| `template_critical_depth_froude_classification` | civil | 2 | 0.4% | `q` x2 |
| `template_safety_stock_reorder_point` | industrial | 2 | 0.4% | `SS` x2 |
| `template_infinite_slope_factor_of_safety` | civil | 1 | 0.2% | `FS` x1 |
| `template_exponential_mttf_topology` | industrial | 1 | 0.2% | `MTTF` x1 |
| `template_system_reliability_topology` | industrial | 1 | 0.2% | `Rs` x2 |
| `template_basic_stress_strain` | mechanical | 1 | 0.2% | `epsilon` x1 |

Examples:

- `template_server_configuration_selection`: seed 0: W1 = 1 / (muA - lambda) = 1 / (44 - 36) hours = 60 / 8 = 7.50 minutes  (exact 0.125)
- `template_manning_trapezoidal_velocity`: seed 100: V = (1/n) * R^(2/3) * S^(1/2) = (1/0.012) * 1.1462 * 0.03000 = 2.865 m/s  (exact 2.86550000000000000000000000000000000000000000000000000000000)
- `template_linear_reservoir_routing_step`: seed 152: O2 = (9720.0 + 5760.0) / 4800.0 = 3.23 m^3/s  (exact 3.225)
- `template_sigma_reduction_for_cpk`: seed 8: reduction = (1 - sigma_max/sigma) * 100 = (1 - 1.197/1.68) * 100 = 28.8 percent (rounded up, so a reduction of this size  (exact 28.7500)
- `template_borrow_pit_fill_volume`: seed 169: gamma_d,borrow = Gs * gamma_w / (1 + e_borrow) = 2.70 * 9.81 / (1 + 0.80) = 14.71 kN/m^3  (exact 14.715)
- `template_rational_method_peak_flow`: seed 23: C_w = (C1*A1 + C2*A2) / (A1 + A2) = (0.62 * 7.8 + 0.81 * 4.2) / 12.0 = 0.687  (exact 0.6865)
- `template_poissons_ratio`: seed 30: epsilon_lateral = -(0.42) * (-1.2750e-03) = 5.355000e-04  (exact 0.000535500)
- `template_critical_depth_froude_classification`: seed 151: q = Q / b = 34.5 / 4.8 = 7.188 m^2/s  (exact 7.1875)
- `template_safety_stock_reorder_point`: seed 61: SS = z * sigma_L = 2.0537 * 250 = 513.43 units  (exact 513.4250)
- `template_infinite_slope_factor_of_safety`: seed 338: FS = 0.5774 / 0.4000 = 1.444  (exact 1.4435)
- `template_exponential_mttf_topology`: seed 411: MTTF = 1000 / 1.28 = 781.3 hours (rounded half up to one decimal)  (exact 781.25)
- `template_system_reliability_topology`: seed 272: Rs = 1 - 0.075 * 0.054 = 1 - 0.004050 = 0.9960  (exact 0.995950)
- `template_basic_stress_strain`: seed 86: epsilon = (0.81 mm) / (540.0 mm) = 1.500e-03  (exact 0.0015)

233,349 checks evaluated exactly; 28,962 skipped (functions or operators outside the exact evaluator, e.g. trigonometry and logarithms).
