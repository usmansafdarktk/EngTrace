# Phase C1.3 — constants-table census and classification

**Generated, not written.** Regenerate with

    python -m tests.constants_integrity.census --seeds 40 --markdown docs/re-implementation-sep/phaseC1_census.md

Every term below - *numeric table*, *literal*, *tagged*, *consumer*, *draws by
order*, *restated*, *hidden* - is a predicate defined in the docstring of
`tests/constants_integrity/census.py`; the classification rule is spec §C1.3.

## Coverage

| branch | numeric tables | tagged | numeric literals | in tagged tables |
|---|---:|---:|---:|---:|
| chemical_engineering | 17 | 3 | 794 | 169 |
| electrical_engineering | 13 | 2 | 42 | 2 |
| mechanical_engineering | 8 | 0 | 250 | 0 |
| civil_engineering | 31 | 21 | 516 | 499 |
| industrial_engineering | 38 | 24 | 812 | 739 |
| **all** | **107** | **50** | **2414** | **1409** |

## Classification

| class | tables | literals |
|---|---:|---:|
| CITATION | 36 | 1511 |
| PLAUSIBILITY | 29 | 191 |
| DERIVATION | 4 | 499 |
| DEFINITION | 5 | 5 |
| DOMAIN | 1 | 1 |
| UNCONSUMED | 32 | 207 |

## Per table (40 seeds per consuming template)

| branch | table | literals | tags | consumers | draw by order | P-GIVEN | @kind | class | why | hidden in |
|---|---|---:|---|---:|---:|---|---|---|---|---|
| chemical | `REAL_FLUID_DATA` | 75 | - | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| chemical | `CRITICAL_PROPERTIES` | 135 | - | 5 | 5 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `SUBSTANCES_FOR_HEATING` | 105 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `SUBSTANCES_FOR_VAPORIZATION` | 16 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `HEATS_OF_FORMATION` | 21 | BY-DEFINITION, ON-DISK | 2 | 0 | SOME-HIDDEN | property | **CITATION** | consumed without being stated | adiabatic_flame_temperature |
| chemical | `REACTIONS` | 16 | - | 4 | 4 | SOME-HIDDEN | mathematical | **DERIVATION** | computable from its definition | batch_moles_vs_conversion, flow_system_molar_flow_rates, heat_of_reaction_formation, limiting_reactant |
| chemical | `CP_PARAMS` | 132 | DERIVED, ON-DISK | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `CP_VALID_T_MAX` | 13 | - | 0 | 0 | UNCONSUMED | validity | **UNCONSUMED** | no template consumes it | - |
| chemical | `AIR_COMPOSITION` | 3 | - | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| chemical | `CP_PARAMS_COMBUSTION` | 16 | DERIVED | 1 | 0 | SOME-HIDDEN | property | **CITATION** | consumed without being stated | adiabatic_flame_temperature |
| chemical | `CP_COMBUSTION_VALID_T_MAX` | 1 | - | 1 | 0 | NO-EFFECT | validity | **DOMAIN** | a validity bound for another table; checked by its C3.2 suite | - |
| chemical | `COMBUSTION_REACTIONS` | 63 | - | 1 | 1 | SOME-HIDDEN | mathematical | **DERIVATION** | computable from its definition | adiabatic_flame_temperature |
| chemical | `COMMON_LIQUIDS` | 60 | - | 6 | 6 | SOME-HIDDEN | property | **CITATION** | consumed without being stated | newtons_law_shear_stress |
| chemical | `COMMON_GASES` | 42 | - | 2 | 2 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `GAS_MOLECULAR_PARAMS` | 51 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `POWER_LAW_FLUIDS` | 44 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| chemical | `GRAVITATIONAL_ACCELERATION` | 1 | - | 1 | 0 | ALL-RESTATED | defined | **DEFINITION** | exact by definition or convention | - |
| electrical | `C0` | 1 | BY-DEFINITION, ON-DISK | 2 | 1 | ALL-RESTATED | defined | **DEFINITION** | exact by definition or convention | - |
| electrical | `MEDIA_VELOCITIES` | 20 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| electrical | `EPSILON_0` | 1 | ON-DISK | 4 | 0 | SOME-HIDDEN | measured-constant | **CITATION** | consumed without being stated | coaxial_capacitance, coulombs_law, gauss_law_symmetric, superposition_electric_field |
| electrical | `FREQUENCY_RANGE_HZ` | 2 | - | 1 | 1 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `AMPLITUDE_RANGE` | 2 | - | 3 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `PHASE_RANGE_DEG` | 2 | - | 2 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `OMEGA_MULTIPLIER_RANGE` | 2 | - | 2 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `SAMPLING_FREQ_RANGE_HZ` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `F0_RANGE_HZ` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `GAIN_K_RANGE` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `DELAY_N0_RANGE` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `DECIMATION_FACTOR_M_RANGE` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| electrical | `OMEGA_DENOMINATOR_RANGE` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| mechanical | `MATERIAL_PROPERTIES` | 81 | - | 4 | 4 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| mechanical | `SHEAR_MODULUS_VALUES` | 18 | - | 3 | 3 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| mechanical | `ATMOSPHERIC_PRESSURE_KPA` | 1 | - | 0 | 0 | UNCONSUMED | defined | **UNCONSUMED** | no template consumes it | - |
| mechanical | `GRAVITY` | 1 | - | 4 | 0 | SOME-HIDDEN | defined | **DEFINITION** | exact by definition or convention | basic_buoyant_force, hydrostatic_force_on_plane, hydrostatic_pressure_at_depth, utube_manometer |
| mechanical | `FLUID_DENSITIES` | 41 | - | 4 | 4 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| mechanical | `MATERIAL_DENSITIES` | 56 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| mechanical | `PIPE_FLUIDS` | 32 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| mechanical | `MANOMETER_FLUIDS` | 20 | - | 1 | 1 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `GRAVITY_M_S2` | 1 | - | 7 | 0 | ALL-RESTATED | defined | **DEFINITION** | exact by definition or convention | - |
| civil | `GRAVITY_FT_S2` | 1 | - | 0 | 0 | UNCONSUMED | defined | **UNCONSUMED** | no template consumes it | - |
| civil | `UNIT_WEIGHT_WATER_KN_M3` | 1 | - | 6 | 0 | ALL-RESTATED | property | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `UNIT_WEIGHT_WATER_PCF` | 1 | - | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| civil | `WATER_DENSITY_KG_M3` | 1 | VERIFY | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| civil | `WATER_KINEMATIC_VISCOSITY_M2_S` | 1 | VERIFY | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| civil | `STEEL_E_KSI` | 1 | VERIFY | 1 | 0 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `STEEL_E_GPA` | 1 | VERIFY | 3 | 0 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `STEEL_FY_KSI` | 2 | VERIFY | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `STEEL_FY_MPA` | 2 | - | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `CONCRETE_FC_PSI` | 3 | VERIFY | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| civil | `CONCRETE_FC_MPA` | 3 | - | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| civil | `CONCRETE_EC_COEFF_PSI` | 1 | - | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `CONCRETE_EC_COEFF_MPA` | 1 | - | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `LIVE_LOADS_PSF` | 5 | VERIFY | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `LIVE_LOADS_KPA` | 5 | - | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `AISC_W_SHAPES` | 196 | ON-DISK | 2 | 2 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `SPECIFIC_GRAVITY_RANGES` | 6 | VERIFY | 5 | 1 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| civil | `PERMEABILITY_RANGES_CM_S` | 10 | ON-DISK, VERIFY | 1 | 0 | NO-EFFECT | range | **PLAUSIBILITY** | NO-EFFECT: consumed as a guard or screen | - |
| civil | `DAS_NATURAL_STATE_SOILS` | 44 | ON-DISK | 0 | 0 | UNCONSUMED | property | **UNCONSUMED** | no template consumes it | - |
| civil | `FRICTION_ANGLE_RANGES_DEG` | 12 | VERIFY | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| civil | `TERZAGHI_BEARING_FACTORS` | 36 | ON-DISK | 1 | 0 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `TERZAGHI_MODIFIED_FACTORS` | 16 | ON-DISK | 1 | 0 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `SKEMPTON_CC_COEFF` | 1 | VERIFY | 1 | 0 | SOME-HIDDEN | standard | **CITATION** | consumed without being stated | primary_consolidation_settlement |
| civil | `SKEMPTON_CC_OFFSET` | 1 | - | 1 | 0 | SOME-HIDDEN | standard | **CITATION** | consumed without being stated | primary_consolidation_settlement |
| civil | `CV_RANGES_M2_YR` | 4 | POLICY | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| civil | `MANNINGS_N_CHANNELS` | 17 | ON-DISK | 4 | 4 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `MANNINGS_N_CONDUITS` | 10 | ON-DISK | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| civil | `RATIONAL_C` | 32 | ON-DISK | 1 | 0 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| civil | `SCS_CURVE_NUMBERS` | 100 | ON-DISK, VERIFY | 1 | 0 | NO-EFFECT | standard | **CITATION** | NO-EFFECT: asserts a fact about a named real entity | - |
| civil | `SCS_IA_RATIO` | 1 | ON-DISK | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| industrial | `Z_QUANTILES` | 12 | DERIVABLE | 1 | 0 | ALL-RESTATED | mathematical | **DERIVATION** | computable from its definition | - |
| industrial | `SHEWHART_K_SIGMA` | 1 | ON-DISK | 2 | 0 | SOME-HIDDEN | defined | **DEFINITION** | exact by definition or convention | p_chart_limits_floor |
| industrial | `QUEUE_SCENARIOS` | 36 | REALISM | 4 | 3 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `QUEUE_COSTS_USD_HR` | 4 | REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `FINITE_CAPACITY_K` | 2 | REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `COMPONENT_RELIABILITY_CLASSES` | 6 | POLICY, REALISM | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `FAILURE_RATE_PER_HR` | 6 | POLICY, REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `MISSION_TIME_HR` | 2 | REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `HOLDING_RATE_PER_YR` | 2 | REALISM | 3 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `INVENTORY_ITEMS` | 30 | REALISM | 3 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `EPQ_PRODUCTION_MULTIPLE` | 2 | REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `DISCOUNT_BREAK_QTY` | 2 | REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `DISCOUNT_STEP_FRACTION` | 2 | - | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `LEAD_TIME_WEEKS` | 2 | ON-DISK | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `NEWSVENDOR_ITEMS` | 24 | REALISM | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `SERVICE_LEVELS` | 4 | ON-DISK | 1 | 1 | NO-EFFECT | range | **PLAUSIBILITY** | NO-EFFECT: consumed as a guard or screen | - |
| industrial | `AGGREGATE_COSTS_USD` | 6 | REALISM | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `AGGREGATE_MONTHLY_DEMAND` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `WORKER_MONTHLY_OUTPUT` | 2 | - | 1 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `LINE_TASK_TIME_S` | 2 | ON-DISK, REALISM | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `LINE_TASK_COUNT` | 2 | - | 0 | 0 | UNCONSUMED | range | **UNCONSUMED** | no template consumes it | - |
| industrial | `LINE_DEMAND_PER_SHIFT` | 2 | - | 1 | 1 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `CONTROL_CHART_FACTORS` | 408 | ON-DISK | 3 | 0 | ALL-RESTATED | mathematical | **DERIVATION** | computable from its definition | - |
| industrial | `MIL_STD_105E_CODE_LETTERS_GII` | 29 | ON-DISK | 2 | 2 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| industrial | `MIL_STD_105E_SAMPLE_SIZE` | 16 | ON-DISK | 2 | 2 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| industrial | `MIL_STD_105E_SINGLE_NORMAL_AC` | 119 | ON-DISK | 2 | 2 | ALL-RESTATED | standard | **CITATION** | ALL-RESTATED: asserts a fact about a named real entity | - |
| industrial | `SPC_CHARACTERISTICS` | 20 | REALISM | 4 | 0 | ALL-RESTATED | range | **PLAUSIBILITY** | a sampling window every consumer restates | - |
| industrial | `P_CHART_PBAR` | 2 | REALISM | 1 | 0 | INDETERMINATE | range | **PLAUSIBILITY** | INDETERMINATE: consumed as a guard or screen | - |
| industrial | `P_CHART_SUBGROUP_N` | 2 | - | 1 | 0 | NO-EFFECT | range | **PLAUSIBILITY** | NO-EFFECT: consumed as a guard or screen | - |
| industrial | `C_CHART_CBAR` | 2 | - | 1 | 0 | NO-EFFECT | range | **PLAUSIBILITY** | NO-EFFECT: consumed as a guard or screen | - |
| industrial | `IEC60063_E24_5PCT` | 34 | - | 1 | 1 | NO-EFFECT | standard | **CITATION** | NO-EFFECT: asserts a fact about a named real entity | - |
| industrial | `IEC60063_E12_10PCT` | 16 | - | 1 | 1 | NO-EFFECT | standard | **CITATION** | NO-EFFECT: asserts a fact about a named real entity | - |
| industrial | `RESISTOR_SERIES_BY_TOLERANCE` | 2 | - | 1 | 1 | NO-EFFECT | standard | **CITATION** | NO-EFFECT: asserts a fact about a named real entity | - |
| industrial | `HARD_ANODIZE_THICKNESS` | 2 | - | 1 | 0 | NO-EFFECT | standard | **CITATION** | NO-EFFECT: asserts a fact about a named real entity | - |
| industrial | `COATING_METROLOGY_FLOOR` | 1 | - | 1 | 0 | INDETERMINATE | range | **PLAUSIBILITY** | INDETERMINATE: consumed as a guard or screen | - |
| industrial | `XBAR_R_SUBGROUP_N` | 2 | ON-DISK | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| industrial | `XBAR_S_SUBGROUP_N` | 2 | - | 0 | 0 | UNCONSUMED | standard | **UNCONSUMED** | no template consumes it | - |
| industrial | `SPC_NUM_SUBGROUPS` | 2 | - | 2 | 0 | NO-EFFECT | range | **PLAUSIBILITY** | NO-EFFECT: consumed as a guard or screen | - |
