# Phase C3.8 — literal copies of constants inside templates: triage

**Tool:** `python -m tests.constants_integrity.literal_copy_sweep` (predicate in its
docstring). **Run:** branch `redesign/phaseC3-remaining-tables`, after the C3.3 fixes.
**Result:** 36 (table, template) candidate pairs, 62 literals. Every one was read at its
source line and put in exactly one of the three outcomes spec §C3.8 allows.

**Known miss of the predicate:** a copy under 3 significant figures is not listed.
`SCS_IA_RATIO`'s `0.2` — the copy that motivated C3.8 (C1 Reviewer G, G-1) — is one; it
is held by its `@copied-in` declaration, which census P-COPY checks on every run.

## Outcome 1 — real copies: declared `@copied-in`

A template writes a table's value as its own literal. Declared, so the census counts
the consumption and fails if the literal disappears; making each template read its
table is a template edit and a P6 event, listed for the phase owner rather than done
here.

| table | template | line | literal | what it copies |
|---|---|---|---|---|
| mechanical `ATMOSPHERIC_PRESSURE_KPA` | `template_hydrostatic_pressure_at_depth` | `fluid_statics.py:36` | `P_atm_kpa = 101.325` | the table's whole value; the module imports the table and does not use it here |
| chemical `CP_VALID_T_MAX` | `template_sensible_heat_temp_dependent_cp` | `heat_effects.py:275` | `uniform(T1 + 20, 350.0)` | `CP_VALID_T_MAX["C6H6(l)"]`, the liquid-benzene fit's upper validity limit, as the draw's top |
| chemical `CP_VALID_T_MAX` | `template_sensible_heat_temp_dependent_cp` | `heat_effects.py:279` | `uniform(T1 + 100, 1200.0)` | `CP_VALID_T_MAX["Al2O3(s)"]`, the corundum fit's upper limit, as the draw's top |

**Found by reading, outside the predicate** (integers, below its 3-s.f. rule), while
triaging the lines above:

| table | template | line | literal | what it copies |
|---|---|---|---|---|
| civil `FRICTION_ANGLE_RANGES_DEG` | `template_terzaghi_strip_footing_bearing` | `strength_and_stability.py:139, 146` | `phi = 35`, `phi = 30` | the "sand, rounded, dense" lower bound and the "sand, rounded, loose" upper bound - declared `@copied-in` like the three above |

`SPC_NUM_SUBGROUPS`'s bounds (20, 30), drawn inline as `randint(20, 30)` by two
templates, are a window rather than a value and go to the C3.9 register.

The two `CP_VALID_T_MAX` copies matter beyond bookkeeping: the table is the D-032
validity record, and a template that hard-codes its limits would keep sampling to
350 K and 1200 K if a refit moved them.

## Outcome 2 — inline windows: to the C3.9 register

Not copies of a table value, but named-entity facts written as the template's own
window bounds — C3.9's class. Recorded there, not declared.

| template | line | literal | note |
|---|---|---|---|
| `template_primary_consolidation_settlement` | `stress_distribution_and_consolidation.py:57` | `Gs_clay = uniform(2.70, 2.80)` | the same bounds as `SPECIFIC_GRAVITY_RANGES["inorganic clay"]`, drawn inline instead of read |
| `template_phase_relations_degree_of_saturation` | `phase_and_index_properties.py:80` | `assert 2.60 <= Gs <= 2.80` | a guard window near, but not equal to, the table's 2.65-2.80 span |
| `template_effective_stress_profile` | `permeability_seepage_effective_stress.py:136` | `gamma_d = uniform(14.5, 16.5)` | a dry-unit-weight window; 14.5 coincides with two `DAS_NATURAL_STATE_SOILS` rows |
| `template_primary_consolidation_settlement` | `...consolidation.py:54, 55` | `uniform(16.5, 18.0)`, `uniform(19.0, 20.5)` | moist and saturated sand unit-weight windows |
| `template_terzaghi_strip_footing_bearing` | `strength_and_stability.py:140, 161` | `uniform(18.0, 19.5)`, `uniform(18.0, 19.0)` | soil unit-weight windows |
| `template_relative_density_of_sand` | `phase_and_index_properties.py:174` | `assert 13.0 <= gamma_d <= 19.0` | a dry-unit-weight guard |
| `template_borrow_pit_fill_volume` | `phase_and_index_properties.py:256, 279` | `uniform(8.0, 18.0)`, its assert | a moisture-content window; 18.0 coincides with a unit weight, not a moisture content |

## Outcome 3 — coincidences

A round loop, draw or assert bound that shares a written token with a table leaf and
has no relation to it. Listed so the count is accounted for.

| leaf matched | template literal(s) |
|---|---|
| `COMMON_LIQUIDS["Honey"][1]` 10.0; `POWER_LAW_FLUIDS` K = 10.0, 25.0, 50.0 | draw bounds for pipe lengths, moles, concentrations (`uniform(10.0, 50.0)`, `uniform(10.0, 25.0)`, `uniform(2.0, 10.0)`, `uniform(20.0, 50.0)`, `uniform(50.0, 150.0)`, `uniform(50.0, 1000.0)`, `uniform(1.0, 10.0)`), a velocity screen `v_max_check < 10.0`, a plot clip `min(v, 50.0)` |
| `AISC_W_SHAPES` rx 14.0, 13.0 | an L/d draw bound `uniform(8.0, 14.0)`, moisture `uniform(6.0, 14.0)`, unit-weight guards `14.0 <= gamma`, `13.0 <= gamma_d`, a truss divisor `c_sum / 14.0` |
| `CV_RANGES_M2_YR["low-plasticity clay"][1]` 30.0 | stress and capacity guards `30.0 <= sigma0`, `30.0 <= q_all` |
| `PERMEABILITY_RANGES_CM_S["clean gravel"][1]` 100.0 | percent conversions `/ 100.0` and a saturation guard `<= 100.0`, across five templates |
| `AMPLITUDE_RANGE[1]` 50.0 | `template_phasor_addition`'s own `uniform(10.0, 50.0)` amplitude draw (a separate window, not this table's 1.0-50.0) |

## Counts, with their predicate

36 pairs: 2 outcome-1 pairs (3 literals), 7 outcome-2 pairs (5 `DAS_NATURAL_STATE_SOILS`,
2 `SPECIFIC_GRAVITY_RANGES`), 27 outcome-3 pairs (5 `COMMON_LIQUIDS`, 9 `POWER_LAW_FLUIDS`,
5 `AISC_W_SHAPES`, 2 `CV_RANGES_M2_YR`, 5 `PERMEABILITY_RANGES_CM_S`, 1 `AMPLITUDE_RANGE`) —
counted pair by pair from `literal_copy_sweep.py` as committed; re-run it rather than
trusting these lines. (The first draft of this line said 8 and 26 - written before the
pairs were counted.)
