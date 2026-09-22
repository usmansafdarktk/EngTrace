# Layer 0 — the closure fixes, template by template

**Date:** 2026-09-23 · **Branch:** `annotation/layer0-layer1` · **Rule applied:** D-016 (a displayed-and-consumed
operand is bound through its display; a half-way display tie is removed, never resolved) and D-037 (lengthen an
exact display before resampling), unchanged from Phases 1–6.

Every template below failed T1 closure at HEAD (`b1445ad`) at 25 or 500 seeds, and passes it after the edit at
5,000 seeds with T3, T4 and T8 passing. The per-template rates in this table were measured by the agent that
made the edit, at 5,000 seeds; the corpus-wide numbers that matter for the paper — the gate outcome and the
item-pool movement — are produced by `gate.py` (`gate_report.md`) and by
`tests/constants_integrity/c3_instance_dump.py --diff` (`item_pool_impact.md`), not copied from here.

Columns: **question moved** is the share of seeds whose question text changed (only ever the redrawn seeds,
except where noted); **answer moved** is the share of seeds whose gold answer changed for an unchanged question;
**redraw** is the rejection rate of the tie screen added, if any.

## Round 1 — the 29 templates failing at 25 seeds

Three of the 29 are parser limits, registered in `check_limits.json` and not edited:
`template_time_rate_of_consolidation`, `template_server_configuration_selection`,
`template_sigma_reduction_for_cpk`. `template_wave_parameters_basic` is edited (its k line) and its omega
line is registered (T1 cannot size a `%e` result; see D-092).

| Template | Cause | Fix | Question moved | Answer moved | Redraw |
|---|---|---|---|---|---|
| elec `lorentz_force` | B stated at 5 s.f., printed `.2e` and consumed at full precision; q printed `.2e` | B displayed `.4e` (exact) and bound; cross product `.5f` (exact) and bound; q `.3e` (exact) and bound | 0 | 0 | none |
| elec `wave_parameters_basic` | u_p, f (`.2e`), pi (5 dp), lambda (3 dp) consumed unrounded | each bound through its display; MHz display from the bound values | 0 | omega/T/lambda/k move by one last digit on 12–70% of seeds | none |
| elec `lowpass_equivalent_bandpass` | one line chained the exact trig expression, its 3-dp substitution and the result | trig factor stated on its own line and bound (spec P3); component recomputed from it | 2.8% (redraws) | ±0.01 on ~20% | 2.84% (A·0.707 ties at 45°/135°) |
| elec `phasor_addition` | components printed 2 dp, added unrounded; phases 1 dp consumed unrounded | components, sums and phases bound through their displays | 0.18% | A_total ±0.02 max on 32%; phase up to 0.96° where phasors nearly cancel | 0.18% |
| elec `time_to_phasor` | phase printed 2 dp, consumed unrounded by cos/sin; components likewise | bound in chain order | 0.10% | rectangular parts ±0.02 on ~10% | 0.10% |
| elec `decimation_aliasing_analysis` | operands printed 3 dp, difference taken unrounded | both operands bound; difference of printed operands | 0 | 0 (symbolic answer unchanged) | none |
| elec `ber_estimation_mary` | Eb/N0 printed at 1–2 dp, consumed unrounded | bound; Es/N0 bound (exact) | 0 | Q argument ±1 last digit on 17.6% | none |
| elec `bpsk_energy_basis` | Tb stated to 3–5 s.f., consumed at full precision | Tb bound at its stated display at sampling | 0.38% | Eb up to 0.43% on 76%; amplitude up to 0.07% | 0.38% |
| chem `sensible_heat_temp_dependent_cp` | Cp coefficients printed `.3e` but consumed at table precision (up to 6 s.f.); I(T) printed 4 dp, subtracted unrounded | coefficients printed at their exact precision (D-037: `.4e`/`.5e`, never shorter); I(T), integral, Q bound | **55% (coefficient displays lengthened)** | 0 | 0.12% |
| chem `pitzer_correlation_z` | Z printed 4 dp, consumed unrounded; Tr, Pr printed from the sampled pre-image | Tr, Pr recomputed forward from stated T, P; B0, B1, ratio, Z bound; ω-sum printed at its exact dp | 0.10% | Z ±1 in 4th dp on 16% | 0.10% |
| chem `work_isothermal_virial` | B printed 5 dp, consumed unrounded; Tr from pre-image | forward Tr; B0, B1, B, Z, V, RT, terms bound | 0 | W ±1 J/mol on 47% | 0.002% |
| chem `gas_viscosity_kinetic_theory` | σ² printed 3 dp, denominator computed unrounded | σ², M·T printed at their exact dp; numerator, denominator bound | 0.10% | μ ±1 in 4th s.f. on 5% | 0.10% |
| chem `batch_moles_vs_conversion` | product printed 3 dp, N_B computed unrounded | products and N_B/C/D bound at 3 dp | 0 | N_B ±0.001 on 4.25% | none |
| mech `angle_of_twist` | radians printed 4 dp, degrees computed unrounded; J, G likewise | angle_rad, J (`.4e`), G (`.2e`) bound | 0 | degrees ±0.001 on 98% | 0 |
| mech `shear_stress_torsion` | Pa printed `.3e` (4 s.f.) divided by 1e6 against a 6–7 s.f. MPa | one rounding at the MPa display; Pa displayed exactly; J bound | 0 | ≤ 5e-4 relative on 71% | none |
| mech `statically_indeterminate` | area, reactions, lengths printed rounded, consumed unrounded; two exact ties | bound in chain order; ties redrawn | 0.62% | stress last digit on 40% | 0.62% |
| mech `poissons_ratio` | area, stress, strains printed and consumed unrounded | bound in chain order; lateral strain `.6e` | 1.98% | ±1 in 5th dp on 15% | 1.98% |
| mech `utube_manometer` | terms printed 2 dp, subtracted unrounded; ρgh exact 2-dp tie | terms and gauge pressure at 4 dp (exact); kPa tie redrawn | 0.12% | 0 | 0.12% |
| mech `hydrostatic_pressure_at_depth` | Pa at 1 dp, kPa answer at 3 dp: a digit dropped across the conversion | **kPa answer at 4 dp** (exact); tie redrawn | 1.06% | answer carries one more digit on 90% | 1.06% |
| civil `beam_support_reactions` | triangular W = w·L/2 exact at 2 dp, printed 1 dp | W at 2 dp (exact) everywhere; By tie redrawn | 3.33% | By, Ay up to 0.05 on 12% | 3.32% |
| civil `force_method_continuous_beam` | end reactions exact at 3 dp, printed 2 dp; By a quotient | end reactions at 3 dp; By tie redrawn | 0.47% | 0 | 0.47% |
| civil `slope_deflection_end_moment` | fixed-end moments exact at 4 dp, printed 2 dp | FEMs and rhs at 4 dp; X, M_AB ties redrawn | 2.3% | X, M_AB by one 2-dp step on 19% | 2.34% |
| civil `max_hump_height_no_choking` | dz = E1 − Ec exact at 4 dp, printed 3 dp (the answer) | **answer at 4 dp**; E1 tie redrawn | 0.005% | answer carries one more digit; 4.5% were resolved ties | 0.005% |
| indu `newsvendor_normal_demand` | z·σ exact at 4 dp, printed 2 dp | z·σ and sum at 4 dp; critical-ratio tie redrawn | 0.07% | 0 | 0.07% |
| indu `qr_policy_one_iteration` | question prescribes every precision; four exact intermediates tie | ties screened and redrawn (no display change) | 4.6% | 0 | 4.64% |
| indu `mmc_waiting_time` | Wq = Lq/λ exact at 6 dp, printed 5 dp | Wq at 8 dp (first tie-free display over all 82 instances) | 0 | 0 | none |

## Round 2 — the rare ties visible only at 500 seeds

| Template | Cause | Fix | Question moved | Answer moved | Redraw |
|---|---|---|---|---|---|
| indu `reorder_point_lead_time` | λ·τ exact at 3 dp, printed 1 dp; T and τ/T quotients | λ·τ at 3 dp; quotient ties redrawn | 0.58% | 0 | 0.58% |
| indu `mm1k_finite_capacity` | λ(1−P_K) exact at 5 dp, printed 3 dp; ρ a quotient; SL and the minute answer tie unseen by T1 | λ_e at 5 dp; ρ, SL, W·60 ties redrawn | 5.6% | ±0.01 min on 0.2% | 5.61% |
| indu `absorbing_chain_time_to_failure` | numerator/denominator exact at 6/8 dp, printed 4 dp; answer a quotient | n1 at 6 dp, d1 at 8 dp; μG, μW ties redrawn | 0.72% | **±0.01 week on 8.3%** (exact operands now feed the quotient) | 0.71% |
| indu `epq_finite_production` | h′ exact at 10 dp, printed 6 dp; H a half-integer on 1.1% | h′ at 10 dp; D/P and H ties redrawn | 1.1% | 0 (grid answer pinned); printed root last digit on 5.6% | 1.09% |
| indu `chart_pair_selection` | question prescribes the display; limits tie | exact Decimals screened before the single quantize | 0.72% | 0 | 0.72% |
| indu `xbar_r_control_limits` | as above | as above | 0.28% | 0 | 0.30% |
| indu `mm1_time_in_system` | ρ, L, W quotients and the 1-dp minute answer tie | ties redrawn | 5.8% | 0 | 5.80% |
| chem `sensible_heat_constant_cp` | kJ answer drops a digit of the stated J value (625.0 J → 0.625 kJ) | ΔT, Q bound; J and kJ ties redrawn (same choice as the sibling template) | 1.54% | 0 | 1.54% |
| civil `terzaghi_strip_footing_bearing` | c′·Nc exact at 2 dp, printed 1 dp (every c′ = 15 draw); width/surcharge products and q_all = qu/3 tie | c′·Nc and qu at 2 dp in general shear; other ties redrawn | 3.4% | 0 | 3.42% (c′ = 15 general shear loses a third) |
| civil `truss_method_of_sections` | member force a quotient by the 1-dp height (the answer) | ties redrawn | 3.3% | 0 | 3.51% (h = 1.6 loses 20%) |
| civil `virtual_work_truss_deflection` | six quotients by height or area tie | ties redrawn | 7.5% | 0 | **7.11%** (h = 1.6 loses 51%, h = 2.4 loses 26%) |
| civil `beam_internal_moment` | M exact at 3 dp, printed 2 dp (the answer); By a quotient | ties redrawn | 10.6% | 0 | **10.55%** (sections at x = k.5 lose half) |
| civil `phase_relations_degree_of_saturation` | γ_d, e, S quotients tie | ties redrawn | 1.3% | 0 | 1.20% |
| civil `truss_method_of_joints` | F = Ay/sin θ ties on 3-4-5 geometries with odd P | ties redrawn | 0.5% | 0 | 0.47% |
| civil `effective_stress_profile` | σ and u exact at 3 dp, printed 1 dp (u ties on every z2 = 5.0 m draw, invisible to T1) | σ, u at 3 dp; σ′ rounded once at the 1-dp answer; γ ties redrawn | 2.1% | **±0.1 kPa on 24%** (a double rounding removed) | 2.01% |
| civil `scs_curve_number_runoff` | Q_mm = Q_in·25.4 exact at 4 dp, printed 1 dp (the answer) | tie redrawn | 0.4% | 0 | 0.21% (one CN/storm pair removed) |
| mech `basic_stress_strain` | area printed 4 dp, stress divides the unrounded area; stress and strain (answers) tie | area bound; answer ties redrawn | 0.5% | stress ±1–4 last digits on 3% | 0.56% |
| mech `multi_segment_rod` | area printed 4 dp consumed unrounded; δ_total summed unrounded deltas while the trace prints 4-dp ones | area, SI operands and each δ bound; total = sum of printed δ | 0.005% | **δ_total ±0.0001 on 28.5%** (the printed sum now equals the printed total; before it disagreed on 28.4%) | 0.005% |
| mech `undamped_response_initial_conditions` | ω_n printed 4 dp, consumed unrounded by A2 | ω_n bound; ω_n, A2 ties redrawn | 0 | A2 ±0.0001 on 0.24% | < 5e-6 |

## Decisions the owner must sign off (D-044)

- **Answer display lengthened** in two templates: `hydrostatic_pressure_at_depth` (kPa 3 → 4 dp) and
  `max_hump_height_no_choking` (m 3 → 4 dp). In both the longer display is exact and the tie was in the answer
  itself; the alternative (redraw ~10% of instances) was judged worse. Reverting is a one-line change each.
- **Gold answers move on unchanged questions** in `absorbing_chain_time_to_failure` (8.3% of instances,
  ±0.01 week): the printed 4-dp operands used to feed the quotient; the exact ones now do. The alternative is a
  4-dp tie screen at 4.6% rejection and no answer movement.
- **Question text moves** for `sensible_heat_temp_dependent_cp` on about half its instances, because the
  heat-capacity coefficients now print at the precision the table gives them. Any archived inference on that
  template is stale.

## Residuals found and left (recorded, not fixed)

- `poissons_ratio` fails its round-trip oracle (T2, informational in this gate) on 4 of 500 seeds: the
  trace's part (a) differs from the full-precision oracle by exactly one unit in the fifth decimal
  (0.56–0.76% at delta_d ≈ 0.0018), because the chain now consumes its displayed area, stress and strains
  (D-016 part 2) while the oracle's 0.5% relative tolerance was sized when the chain ran at full precision.
  The oracle's own comment says it needs a display-unit tolerance to be meaningful; that restatement is the
  harness owner's (D-074), so it is recorded here rather than patched.
- The exact tie census (`tie_census.md`) finds residual half-way ties, filed by T1 as marginal rather than
  failing, in 22 templates at 500 seeds: about 0.5% of instances corpus-wide, above 3% in eight templates
  (`mmc_waiting_time` 20%, `primary_consolidation_settlement` 10%, `flow_system_molar_flow_rates` 7%,
  `upward_seepage_quick_condition` 5%, `takt_time_line_efficiency` 4%, `relative_density_of_sand` 3%,
  `influence_line_max_reaction` 3%, `rotating_unbalance` 2%). Same defect class as this pass, same remedy;
  a third round is the owner's call.

- `work_isothermal_virial` raises `ValueError: math domain error` on about 1 seed in 10,000 (ln of a negative
  volume ratio when B is strongly negative); identical at HEAD. A physical guard in the redraw loop would remove
  it; out of this pass's scope.
- `batch_moles_vs_conversion`: the 3-dp answer is a float-resolved tie on 10–50% of draws depending on the
  reaction; removing it needs a 4–6 dp answer display (D-044) or a 10–60% redraw. Left as is.
- `force_method_continuous_beam`: `delta_BB` (3 dp) ties on 16.6% of instances and `delta_B0` on 0.6%; T1 cannot
  see those lines (a trailing `/EI` poisons the segment). A 4-dp `delta_BB` would move By on 17% of seeds.
- `mmc_waiting_time`: `rho = a / c` prints a 4-dp tie on 18% of instances (always MARGINAL to T1); a 5-dp rho
  would change P0, Lq and the answer on those instances.
- `mm1k_finite_capacity` / `mm1_time_in_system` / `mmc_waiting_time`: every `X hours * 60 minutes/hour = Y`
  line is invisible to T1 (poisoned by `/hour`); the first two now screen the minute answer, `mmc_waiting_time`
  ties there with probability ~1e-5 and is unscreened.
- `lorentz_force`, `bpsk_energy_basis`, `gas_viscosity_kinetic_theory`, `shear_stress_torsion`: T5b reports
  `.Ne` displays as N-dp roundings (D-017, a known harness misparse); unchanged by this pass.
- `sensible_heat_temp_dependent_cp` and other chemical templates: `${X_A*100} \%` can print float noise such
  as `56.99999999999999 %`; pre-existing, not a T1 line.
- T7 (asserts) fails 78 legacy templates by design: the chemical, electrical and mechanical branches carry no
  asserts. Advisory, not gating.
