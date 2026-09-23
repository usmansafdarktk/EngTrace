# Layer 0 gate report

Generated 2026-09-23T19:31:26+00:00 by `gate.py` at git `8224d8baec`, 500 seeds per template, 154.5 s. Gating: T1, T3, T4, T8 plus generation errors. Advisory: T5, T7.

| | Templates |
|---|---:|
| in corpus | 150 |
| **pass the gate** | **150** |
| fail the gate | 0 |
| failing T1 | 0 |
| failing T3 | 0 |
| failing T4 | 0 |
| failing T8 | 0 |
| failing generation error | 0 |
| T1 closure lines excused by the register | 1494 |
| T2 oracles run (informational) | 13 (12 pass) |
| T5 advisory, not passing | 60 |
| T7 advisory, not passing | 83 |

| Branch | pass | fail |
|---|---:|---:|
| chemical_engineering | 30 | 0 |
| civil_engineering | 30 | 0 |
| electrical_engineering | 30 | 0 |
| industrial_engineering | 30 | 0 |
| mechanical_engineering | 30 | 0 |

## Register of accepted check limits

Closure failures excused because the check cannot read the line, not because the line is wrong. Each entry names the template, the line pattern and the reason; the count is how many lines it absorbed at 500 seeds.

| Template | Pattern | Lines absorbed | Reason |
|---|---|---:|---|
| `template_time_rate_of_consolidation` | `H_dr,field / H_dr,lab = .* m / .* mm = .* mm / .* mm =` | 500 | A mixed-unit ratio rewritten across '=' (5.40 m / 12.5 mm = 5400 mm / 12.5 mm = 432.0). T1 strips the unit tokens, so it cannot see that the second segment is an intermediate expression; it reads its leading number (5400) as the chain's terminal value and compares it with 5.40 / 12.5. Recomputed with units, every such line closes. |
| `template_server_configuration_selection` | `^W\d = 1 / \(mu[A-Z] - lambda\) = 1 / \(\d+ - \d+\) hours = 60 / \d+ = [0-9.]+ minutes` | 111 | An hours-to-minutes conversion written across '=' (1 / (40 - 26) hours = 60 / 14 = 4.29 minutes). T1's unit-rescale detector allows a clean integer ratio only within 0.05%, but the terminal value is rounded to 2 dp, so the ratio it sees is 60.06 rather than 60 and the line is reported as a failure. Recomputed as stated, every such line closes. |
| `template_server_configuration_selection` | `^W\d = L\d / lambda = [0-9.]+ / \d+ = [0-9.]+ hours = [0-9.]+ minutes` | 8 | The same hours-to-minutes conversion across '=' in the second server's waiting time (2.4003 / 12 = 0.20003 hours = 12.00 minutes): the minutes value is rounded to 2 dp, so the ratio T1 sees is 59.99 rather than 60. Recomputed as stated, every such line closes. |
| `template_wave_parameters_basic` | `^omega = 2 \* 3\.14159 \* \([0-9.]+e\+\d+ Hz\) = [0-9.]+e\+\d+ rad/s` | 500 | The result is printed in scientific notation (6.85e+09 rad/s) and core.printed_precision() drops the exponent, so T1 sizes the tolerance from the mantissa's two decimals as 0.005 absolute and the 1e-9 relative floor then dominates (about 7 rad/s on 7e9). That is the closure counterpart of D-017, where T5b misreads a %e spec the same way. With every operand bound through its display (done 2026-09-23), every omega line closes within one half-unit of the printed mantissa's last digit (worst 0.992 half-units at 5000 seeds); only the check's tolerance is wrong. |
| `template_sigma_reduction_for_cpk` | `^sigma_max = [0-9.]+ / [0-9.]+ = [0-9.]+ \S+ \(rounded down` | 206 | The line states in prose that the quotient is rounded DOWN so the resulting sigma still meets the capability target (10.6 / 5.01 = 2.115, not 2.116). That is the template's stated rule, and the answer depends on it; T1 applies half-unit tolerance and cannot read the stated rounding direction. |
| `template_sigma_reduction_for_cpk` | `^reduction = \(1 - sigma_max/sigma\) \* 100 = .* percent \(rounded up` | 169 | The companion line rounds the percentage reduction UP, stated in prose, so that a reduction of the quoted size does reach the target. Same stated-direction rounding as the entry above; T1 cannot read it. |

## Closure marginals (informational)

99 templates carry at least one marginal line (within tolerance, on a rounding boundary). Highest rates: `template_scs_curve_number_runoff` 64%, `template_best_hydraulic_rectangular_section` 62%, `template_rational_method_peak_flow` 61%, `template_upward_seepage_quick_condition` 59%, `template_hydraulic_jump_energy_loss` 59%.
