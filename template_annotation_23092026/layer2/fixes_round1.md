# Layer 2, round 1: the 22 rejected templates, verified and fixed

Written 2026-09-25. `RESULTS.md` lists the 22 templates at least one expert rejected (15 by
majority) and every rejection note. Each note was treated as a claim, checked against the source
and the rendered instances (seeds 2101–2105, the ones the experts saw) before any edit, in the same
way as the screen's flags (D-094). Twenty templates were changed; two claims were not adopted and
the reason is recorded in the template's docstring and below. Every edited template carries a
"Layer 2 fix (2026-09-25)" docstring paragraph. Per-template rates in the tables are the fixing
agents' measurements at 500 seeds unless marked; the corpus-wide numbers at the end are the
committed scripts'.

## Corpus-wide result

| Check | Command | Result |
|---|---|---|
| Gate | `python -m template_annotation_23092026.layer0.gate --seeds 500` | **150 of 150 pass**; register absorbed 1,494 lines; advisory T5 60, T7 81 (was 83) |
| Tie census | `python -m template_annotation_23092026.layer0.tie_census --seeds 500` | same 13 templates with a tie as at HEAD, none of them among the 22; no new tie |
| Item pool | `tests/constants_integrity/c3_instance_dump.py --diff`, 150 × 300 seeds, HEAD `30e6988` dumped from a detached worktree vs the fixed tree | **20 templates moved**: question 2,929 / answer 3,850 / solution 4,718 of 45,000; errors 0 → 0 |
| Cumulative since `b1445ad` | same, against the pre-Layer-0 dump | 77 templates moved: question 4,843 / answer 10,147 / solution 16,514 (`../layer0/item_pool_impact.md`) |

Only the 22 assigned `template_*` functions differ from HEAD (checked function by function with
`ast`); the other changes are two constants tables (below) and closure helpers added to three
mechanical files.

## Constants changed

- `mechanical_engineering/constants.py`: **MANOMETER_FLUIDS** loses the "Molten Metals" group
  (gallium, tin, zinc; solids at room temperature), 17 → 14 rows. New sampling-only tables
  **ALLOWABLE_NORMAL_STRESS_MPA** (22 rows) with `NOT_TENSION_MEMBER_MATERIALS` (natural rubber,
  lead, concrete, borosilicate glass, alumina) and **ALLOWABLE_SHEAR_STRESS_MPA** (16 rows) with
  `NOT_SHAFT_MATERIALS` (concrete, glass), tagged `[POLICY: sampling-only]`: conservative typical
  allowables (about σy/2 for metals, σult/3–4 for polymers and brittle materials, after Hibbeler
  App. B and Callister App. B). They bound the load a template may draw and are never printed.
  Module-level asserts keep them in step with `MATERIAL_PROPERTIES` and `SHEAR_MODULUS_VALUES`.
- `chemical_engineering/constants.py`: **POWER_LAW_FLUIDS["Blood (Plasma)"]** was (0.012, 0.95),
  about eight times too viscous; now (0.0012, 1.0), Newtonian at 1.2 mPa·s (Késmárky et al.,
  Clin. Hemorheol. Microcirc. 39 (2008); Merrill, Physiol. Rev. 49 (1969)). Documented under the
  table-level UNVERIFIED tag rather than a per-row tag, because a new bracketed tag moves the
  register's chemical count and `test_registers_reconcile` would flag it; the owner may prefer a
  per-row tag plus a register update.

## Per template

Shares are of 500 seeds. "Q" is the question text, "gold" the answer block.

### Chemical (9 assigned: 7 fixed, 2 not adopted)

| Template | Votes | Claim, verified | Fix | Moved |
|---|---|---|---|---|
| gas_viscosity_kinetic_theory | 3/3 | Ω_μ drawn `uniform(0.95, 1.05)`, ε/k unused; seed 2101 propane Ω should be 1.343 | T* = T/(ε/k) bound at 4 dp; Ω_μ from the Neufeld–Janzen–Aziz correlation (BSL eq. 1.4-14), bound at 4 dp, asserted 0.3 ≤ T* ≤ 100; the question still gives Ω_μ (now the computed value) and the solution states its provenance | Q 100%, gold 99.8% |
| newtons_law_shear_stress | 3/3 | μ absent from the question; raw floats (`0.0034000000000000002`); τ printed at 3 dp lost figures | μ stated in the question (`.2e`, `.3e` for 4-figure table values, bound); Y and dv/dy bound through display with a tie redraw; τ and F computed in Decimal and printed in scientific notation at their exact length (a fixed `.4e` tied on half the 2-sf viscosities, so D-037 lengthening was chosen over resampling) | Q 100%, gold 100%; redraw 0.7% |
| annulus_flowrate | 2/3 | honey at up to 526 m/s and 2.9e7 kPa; corn syrup, glycerol, gear oil likewise | velocity target capped at 3 m/s (Sinnott §5.4.3) and ΔP at 500 kPa via `v_hi = min(...)`, algebraically the old value for unbounded fluids; no fluid deleted | Q 23.8%, gold 23.8%; bound binds on 25.2% of draws; T2 oracle 500/500 |
| hagen_poiseuille_flowrate | 2/3 | 1.0% of answers print `0.000000`, 9.4% one figure; radius artefacts on 27%; also found: stated cP at 2 dp did not close for 4-figure fluids | Q bound `.4e`; radius bound `.4f`; cP printed at its exact length; Re < 2300 kept plus v ≤ 3 m/s redraw | Q 9.4%, gold 100% (format); cap rejects 2.1% of raw draws |
| falling_film_max_velocity | 2/3 | seed 2104 gear oil 5.758e-4 printed 0.00058; 0.1% print zero; δ raw floats | v_max bound `.4e` with a sci-tie redraw; δ bound `.5f`; attempt loop bounded; no sampling change | Q 0%, gold 100% (format) |
| power_law_fluid_shear | 2/3 | plasma η = 0.0095 Pa·s (table K eight times too high); air prints η = 0.0000; exponent artefacts | constants row corrected (above); τ, η `.4e` bound with sci-tie redraw; Y, dv/dy bound; exponent printed as `_hu(n − 1)` | Q 5.0% (plasma draws), gold 100% (format) |
| flow_system_molar_flow_rates | 1/3 | Θ-method line `120.26(1.665 − 1.25·0.67)` = 99.515 printed 99.572 | Θ_B, Θ_C bound at 4 dp; the Θ-method line computed from the printed Θ in Decimal at its exact 6 dp, with a note that the small difference is the rounding of Θ; the direct method stays the answer | Q 0%, gold 0%, solution 100% |
| work_isothermal_virial | 1/3 | organics at 900–1188 K on all five seeds (confirmed; pool mean Tr 2.47) | **not adopted**: no on-disk thermal-stability source, and a 750 K cap for organics would reject 57% of valid draws and remove 10 of 27 substances, because at Tr ≤ 1.5 no P2 in the sampled 2.5–5 Pc window satisfies the Vr ≥ 2 validity filter; a 900 K cap still loses 4. The P2/Pc range is the owner's open decision (D-094) | none |
| pitzer_correlation_z | 1/3 | 3 of 5 seeds at 889–1064 K (confirmed) | **not adopted**: a correlation exercise with its validity condition Vr ≥ 2 enforced, the position the screen-pass-1 record already took; a 750 K cap would reject 42% of valid draws | none |

### Electrical (3 assigned, 3 fixed)

| Template | Votes | Claim, verified | Fix | Moved |
|---|---|---|---|---|
| time_to_phasor | 2/3 | Step 2 says "sine leads cosine by 90°", contradicting sin θ = cos(θ − 90°); `round()` drops trailing zeros | wording → "lags"; every 2-dp display via `.2f`; rectangular form built with the sign on the operator | Q 0%, gold text 66% (trailing zeros only, no value moves) |
| wave_parameters_basic | 2/3 | media table is optical n applied at 50–500 MHz (water at 2.25e8 m/s); f rounded to 3 s.f. before use | medium redrawn from the eight media whose RF εr ≈ n² within 1% (vacuum, air, helium, CO2, benzene, CS2, polystyrene, diamond); f bound through `.3e` and printed so; D-016 tie redraw on λ; formulas, ranges and the wavelength-first branch unchanged | Q 62%, gold 73%; registered ω line still the only T1 residual |
| finite_convolution | 1/3 | origin never stated; the asterisk marker undefined | question states the origin convention (as `signal_operations` does), gives h[0] and x[0] explicitly, defines `*` as convolution; Step 3 restates the convention | Q 100%, gold 0% |

### Mechanical (10 assigned, 10 fixed)

| Template | Votes | Claim, verified | Fix | Moved |
|---|---|---|---|---|
| angle_of_twist | 3/3 | nylon at 743 MPa and 96 rad; lead at 17.5 MPa; 0.0062 rad → 0.3552° vs exact 0.3548° | concrete and glass redrawn; torque capped at τ_allow·πd³/16 with a bounded redraw; φ quoted at max(4 dp, 4 s.f.) and degrees from the displayed radians | Q 63%, gold 93% (77% by value); no attempt rejections |
| statically_indeterminate_shaft | 3/3 | concrete shaft at 50 MPa (G > 20 GPa filter admits concrete and glass) | concrete and glass redrawn; applied torque capped at τ_allow·πd³/16 with a redraw in whole hundreds | Q 28.6%, gold 13.2%; T2 oracle passes |
| axial_deformation | 3/3 | polycarbonate at 162 MPa, alumina at 204 MPa, nylon at 2.2% strain; 0.011 in vs 0.01085 | excluded materials redrawn; P ≤ σ_allow·A (or δ ≤ σ_allow·L/E) with a bounded redraw; δ and P at max(3 dp, 4 s.f.); area and load bound through their displays | Q 67%, gold 85% (76% by value); attempt rejections 1.8% |
| statically_indeterminate | 2/3 | 311 kN on a 56 mm ABS bar; rubber, concrete, lead drawable | exclusions by redraw; P ≤ σ_allow·A·L_AC/max(L_AB, L_BC) with a bounded redraw | Q 52%, gold 50%; rejections 0.8% (unchanged) |
| multi_segment_rod | 1/3 | aluminium at 288 MPa; concrete in tension | exclusions by redraw; when max\|P_i\|/(A_i σ_i) > 1 the node loads are scaled down together (a plain redraw rejected 88.6% of draws and would have concentrated the pool, D-045) | Q 92%, gold 92%; rejections 6.2% (was 52.1% from the strain gate alone) |
| utube_manometer | 3/3 | zinc as the manometer liquid on three of five seeds | constants change only (molten metals removed) | Q 96% (list re-indexing, D-031) |
| basic_stress_strain | 2/3 | SI strain labelled "mm/m or unitless" | label → "mm/mm or unitless" | Q 0%, gold text 51% (label), values 0% |
| basic_buoyant_force | 1/3 | aluminium in mercury, cork in ethanol called "fully submerged" | the question states the restraint ("held fully submerged by a rigid clamp"); ρ_f g V is then exact for any body, and no material is removed | Q 100% (wording), gold 0% |
| vibration_transmissibility | 1/3 | base acceleration up to 190 g (38% of draws above 10 g); `0.0029500000000000004 m` | Y·ω² bounded at 1 g with a bounded redraw of Y alone; metre value printed at its exact length | Q 80%, gold 78%; rejections 2.7% (was 3.3%); T2 oracle passes |
| volumetric_flow_rate | 2/3 | Q and A at a fixed 4 dp: A = 0.001 for 9.621e-4, Q/A line unfinished; 1 in 9 pipe draws over 1% off | Q and A at max(4 dp, 4 s.f.) bound through display; the Q/A line finished from the displayed operands and reconciled with U_max/2; tie screens on Q, A, Q/A | Q 3%, gold text 64%, gold value 27% |

## What the experts do next

Round 2 re-certifies only these 22 templates, without plants, with fresh hand-check instances
(seeds 2201–2205), by the same experts: `build_tasks --only <ids> --round 2` then
`make_kits --round 2`. The two templates left unchanged (virial, Pitzer) go back with the docstring
reason so the rejecting expert can see the argument and answer it.

## Open for the owner

- `work_isothermal_virial`: the P2/Pc range (D-094, still open) now also decides whether a thermal
  cap is feasible.
- The plasma row: per-row tag plus register update, or leave it under the table-level tag.
- Whether the screen re-judges the 22 changed templates before round 2 (a few cents; the runner
  refuses a third pass, so it would be a targeted re-judge, D-094's rule).
