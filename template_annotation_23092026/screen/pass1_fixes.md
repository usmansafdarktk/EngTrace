# Screen pass 1 — what was done with the 24 flags

**Date:** 2026-09-24 · **Input:** `pass1/flagged.md` (three judges' sentences per flagged template)
· **Method:** every judge claim was verified against the code and against generated instances
(the screen's seeds 1001–1003 plus 500 seeds) before any edit; confirmed claims were fixed under
the closure rules of D-016/D-037, rejected claims are recorded with the evidence. Each edited
template carries a "Screen pass 1 (2026-09-23)" docstring paragraph. Verification per template:
T1 0 failures, T3/T4/T8 pass, T2 where an oracle exists, no census tie, at 500 seeds; the
corpus-wide gate, census and item-pool diff are re-run after this round (`layer0/`).

**Tally of the 45 individual claims:** 34 confirmed (fixed), 9 rejected with evidence, 2 judge
artefacts of the prompt's design (a helper looks undefined because the prompt shows the function
without its module imports).

| Template | Claims → verdict | Change | Instances moved (of 500) |
|---|---|---|---|
| chem `batch_moles_vs_conversion` | float noise in `%` (confirmed); N_A rounded while others bound (confirmed) | percentage formatted; N_A bound; answers lengthened to their exact 4–6 dp on the owner's policy | question 0; answer text on ties only |
| chem `gas_phase_concentration` | placeholder species names (confirmed); ε sampled independently of stoichiometry (confirmed) | 24 real balanced gas-phase reactions; ε derived from δ and y_A0 and stated as such; tie screens | question 500, gold 500; redraw 2.0% |
| chem `pitzer_correlation_z` | sampled beyond the correlation's range (partly: 15 of 500 fail Vr ≥ 2, Smith–Van Ness–Abbott Fig. 3.14); a rounding claim (rejected: the 7-dp value is exact) | Vr ≥ 2 enforced by redraw; ranges unchanged | 15; redraw 2.7% |
| chem `power_law_fluid_shear` | every fluid labelled non-Newtonian (confirmed, 13% of instances have n = 1) | class and consequence sentence conditional on n | question 66 (n = 1 cases), gold 0 |
| chem `work_isothermal_virial` | mixed virial forms (confirmed: Z from the pressure-explicit form, work integrated in the volume-explicit form, so the printed deviation was an artefact); states outside validity (confirmed: 258 of 500 fail Vr ≥ 2, 100 with no real root); printed intermediates inconsistent (rejected: they close) | volume-explicit form throughout, Z the physical root of Z² − Z − BP/RT = 0; Vr ≥ 2 enforced; deviation sentence corrected | question 259, gold 465; **redraw 49.9%** (owner may narrow P2/Pc) |
| chem `batch_reactor_second_order` | independently rounded intermediates (confirmed, 106 of 500 lines) ; also found: stated conversion disagreed with the stated concentrations on 445 | reciprocals bound; conversion computed from stated values; tie screens | question 17, gold 58; redraw 2.5% |
| chem `ideal_gas_volume` | 1-dp kPa equated to an independently rounded Pa (confirmed, 494 of 500) | Pa = kPa × 1000 exactly; V bound | question 0, gold 162 (last digit) |
| chem `pfr_volume_changing_rate` | float noise in the question's conversion (confirmed, 43 of 500) | percentage formatted | question 43, gold 0 |
| chem `reynolds_number_flow_regime` | a calculation discrepancy (rejected: printed Re matches the operands); found instead: 4-figure viscosities printed at 3 figures and consumed unrounded (24 of 500) | viscosity printed at its exact figures and bound; substitution and result on one line so T1 sees it | question 33, gold 9; redraw 1.5% (ties concentrate on fluids whose ρ/μ terminates) |
| elec `gauss_law_symmetric` | signed ρ shown equal to a positive |D| (confirmed, 107 of 500); numpy array formatting (confirmed); negative zeros (confirmed, 219 of 500) | |ρ| in the magnitude line with a direction clause; points printed as tuples; exact zeros print as zero | question 500 (format), answer text 219; values unchanged |
| elec `impulse_response_from_lccde` | inconsistent rounding of closed-form coefficients (confirmed); missing multiplication signs (confirmed); `fmt()` ignoring its argument (confirmed); "decaying" for |r| ≥ 1 (confirmed, all first-order seeds) | roots bound at 2 dp and the 2×2 system solved from them; explicit `*`; helper removed in favour of `signed_term`; wording conditional on |r| | question 0, answer text 493, C1/C2 values on 31 of 246 second-order seeds |
| elec `signal_operations` | origin marker dropped when n = 0 leaves the support (confirmed, 16%; the D-050 deferral); asterisk notation undefined (confirmed); leftover review comment in text (rejected: a code comment) | printed support widened to include n = 0 (closes D-050); notation defined in the question | question 500 (one clause), answer 79 (representation only) |
| elec `wave_equation_interpretation` | equation printed twice (confirmed); k and ω rounding inconsistent with the stated numbers (confirmed, up to 375 of 500 lines); direction convention inverted (rejected: cos(ωt − kz) → +z is what the code prints, 500 of 500) | printed once; k bound at 2 dp, λ and u_p from the bound k, f from the exact symbolic ω | question 46, gold 287 (last digits) |
| elec `average_energy_mqam` | displayed product ≠ printed result (confirmed, 179 of 254 non-integer seeds) | integer coefficient; A² exact at 4 dp and bound; one rounding; tie redraw | 21; redraw 4.2% |
| elec `cd_dc_system_analysis` | displayed delay product false (confirmed, 118 of 500) | T bound at its display; delay printed exactly at 6 dp | question 0, gold 296 (118 by binding, rest an extra exact digit) |
| elec `continuous_to_discrete_conversion` | degree phase in a radian argument (rejected: the unit is printed on every seed) | none (docstring note) | 0 |
| mech `basic_buoyant_force` | nonsensical material–shape pairings (confirmed by all three; 32% of draws named a shape with its own material word, 25% a hollow body as a solid) | shapes drawn from the eight geometric solids only | **question and answer 482 of 501** (the shorter list re-indexes the draw) |
| mech `composite_shafts_series` | glass shafts and yield-exceeding draws (partly: no strength data in the constants table) | Glass and Concrete excluded by redraw; yield question recorded as residual | 140 of 501 |
| mech `floating_object_submersion_depth` | stability assumed (partly: the cited aerogel block is stable, GM = +59 m; but 41% of draws had GM ≤ 0) | metacentric condition GM > 0 enforced by redraw | 204 of 501; **redraw 40.6%** |
| mech `fluid_particle_acceleration` | unrealistic magnitudes; no continuity (both rejected: max |V| 338 m/s, subsonic; the question claims no incompressibility) | none (docstring note) | 0 |
| mech `poissons_ratio` | Hooke's law on concrete in tension (confirmed: concrete at 118 MPa, glass 267, alumina 1423, lead 63; rubber always beyond the load floor; 29 of 500 outside the template's own strain window) | five materials excluded; realised strain kept in the stated 0.1–0.4% window | 102 of 501; T2 residual 1 of 500 (was 4 at HEAD) |
| mech `undamped_natural_frequency_translational` | 1/0.82 shown as 1.219 (confirmed: f_n and τ_n lines missed on 4–8%) | ω_n, f_n, τ_n bound in chain order; tie screen | question 0, gold 57 (last digit) |
| indu `server_configuration_selection` | "a applicant" (confirmed) | article chosen by the noun's first letter | question 130 (the passport-office draws), gold 0 |
| civil `truss_method_of_joints` | rounding path unstated (confirmed as a description, not a closure defect); 3.0/4.243 shown as 0.7070 (rejected: 0.70705 × 4.243 > 3, so 0.7070 is the correctly rounded quotient) | one sentence states the rounding path | question 0, gold 0 |

## For the owner

- **Large distribution effects to be aware of:** `work_isothermal_virial` now rejects half its
  first draws to stay within the two-term virial's validity (light gases depleted by about a
  third; narrowing the P2/Pc range would cut this); `floating_object_submersion_depth` rejects
  41% for upright stability; `basic_buoyant_force`'s and `gas_phase_concentration`'s instance
  pools are effectively new.
- **Physics corrected, not just presentation:** the virial-work derivation (above), the ε–δ–y_A0
  relation in `gas_phase_concentration`, the |ρ| direction handling in `gauss_law_symmetric`, the
  stability condition in `floating_object_submersion_depth`, the elastic-range condition in
  `poissons_ratio`.
- **Residuals recorded, not fixed:** `composite_shafts_series` has no yield data to bound torque by;
  `poissons_ratio`'s oracle still disagrees by one unit in the fifth decimal on one seed in 500;
  the census reads a `%e` mantissa's digits as decimal places (D-017 class) and so reports false
  ties on lines such as `epsilon_lateral = ... = 5.355000e-04`.
- **D-050 is closed** by the `signal_operations` fix (the deferred origin-marker change is applied).
