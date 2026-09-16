# Phase C3.9 — inline-window register

**What this is.** Spec §C3.9: named-entity facts that enter items as a template's own
literal bounds instead of through a tagged table. Each entry was read at its source
line on branch `redesign/phaseC3-remaining-tables` after the C3.3 and C3.7 commits. An
entry is **moved** (into a tagged table, or read from one), **declared** (a
`@copied-in` the census checks), or **registered** here with its reason. Moving one is
a template edit and a P6 event, so every "proposed" line below is a recommendation for
the phase owner, not a change made.

Sources: the four entries the spec names (C1 Reviewer G, G-3) and the seven C3.8
triage found as outcome 2 (`phaseC3_literal_copies.md`), plus the two integer window
copies found while reading those lines.

## Resolved in C3

| window | where | outcome |
|---|---|---|
| `two_phase_specific_volume`'s V_l, V_v - invented with `uniform(0.001, 0.002)`, `uniform(0.05, 2.0)` | `volumetric_properties_pure_fluids.py` | **moved**: the template reads `REAL_FLUID_DATA[substance]`, re-derived from NIST (39b19b2, 1f14e5d; P6 measured, 300/300 moved) |
| Terzaghi's `phi = 35` and `phi = 30` | `strength_and_stability.py:139, 146` | **declared**: the endpoints of `FRICTION_ANGLE_RANGES_DEG` "sand, rounded, dense" (35-38) and "sand, rounded, loose" (27-30), written as literals - copies, not windows. `@copied-in` on the table, checked by census P-COPY |

## Registered - the spec's named entries

| window | where | what it is | proposed |
|---|---|---|---|
| `e_ranges` sand 0.45-0.85, silt 0.50-0.90, inorganic clay 0.55-0.95 | `phase_and_index_properties.py:46-50` (drawn at :54) | void-ratio windows the code comment says are "anchored to Das Table 3.1". Measured against `DAS_NATURAL_STATE_SOILS` (that table): only 0.45 is a table value (dense uniform sand); the sand top 0.85 is not the table's loose-uniform-sand 0.8; the table has no silt row; the clay window does not span its soft-clay 0.9-1.4. "Anchored" overstates the relation | make the windows a tagged `@kind: range` table with a `[POLICY: sampling-only]` tag and a note of which endpoints are Das values; do not cite Table 3.1 for the rest |
| permeability sampling window [0.02, 0.09] cm/s | `permeability_seepage_effective_stress.py:57-58` (docstring :33) | a template-own sub-window for drawing k, intersected with a geometry-feasibility bound; the template separately asserts k inside `PERMEABILITY_RANGES_CM_S["coarse sand"]` (0.01-1.0) at :67-68, which it does read | keep as sampling policy; move the two bounds into a named range next to `PERMEABILITY_RANGES_CM_S` so the census sees them |
| chart-pair subgroup windows n in [4, 8] (small) and [13, 20] (large) | `variables_control_charts.py:446-447` (docstring :402-411) | the template's own reading of Montgomery §6.3's "n > 10 or 12" (13 clears both readings). `XBAR_R_SUBGROUP_N` (2, 10) and `XBAR_S_SUBGROUP_N` (11, 25) exist in industrial `constants.py` and are read by no template (census UNCONSUMED) - the tables and the template disagree | the repo owner decides which is authoritative: have the template read the two tables (moves items: n = 9, 10, 11, 12 enter), or narrow the tables to the template's windows |

## Registered - found by C3.8 (outcome 2)

| window | where | relation to a table |
|---|---|---|
| `Gs_clay = uniform(2.70, 2.80)` | `stress_distribution_and_consolidation.py:57` | the same bounds as `SPECIFIC_GRAVITY_RANGES["inorganic clay"]`, drawn inline instead of read - proposed: read the table (byte-identical while the bounds stand) |
| `assert 2.60 <= Gs <= 2.80` | `phase_and_index_properties.py:80` | a guard wider than the table's union 2.65-2.80 |
| `gamma_d = uniform(14.5, 16.5)` | `permeability_seepage_effective_stress.py:136` | a dry-unit-weight window above the water table; no table |
| `e = uniform(0.45, min(0.75, e_above - 0.02))` | `permeability_seepage_effective_stress.py:144` | 0.45 is Das's dense uniform sand e; 0.75 is the template's own |
| `gamma_m = uniform(16.5, 18.0)`, `gamma_sat_sand = uniform(19.0, 20.5)` | `stress_distribution_and_consolidation.py:54-55` | moist and saturated sand unit weights; no table |
| `gamma = uniform(18.0, 19.5)` / `uniform(16.5, 17.5)` (sand); `uniform(18.0, 19.0)` / `uniform(16.5, 17.5)` (clayey fill) | `strength_and_stability.py:140, 147, 161, 165` | compaction-state unit weights; no table |
| `assert 13.0 <= gamma_d <= 19.0` | `phase_and_index_properties.py:174` | a dry-unit-weight guard |
| `assert 14.0 <= gamma <= 23.0` | `phase_and_index_properties.py:81` | a moist-unit-weight guard (listed by C3.8 as a coincidence with an AISC rx; a window here) |
| `w_pct = uniform(8.0, 18.0)` and its assert | `phase_and_index_properties.py:256, 279` | a borrow-pit moisture window |

The soil unit-weight windows are the largest class: nine windows over five templates,
none backed by a table, all correctness-neutral in the sense that each drawn value is
stated in its question. Whether they need a source at all is the same question the
C1.3 given-values rule answers for tables - they are sampling policy - but they are
invisible to the census because they live in template code. Proposed for the owner:
a civil `SOIL_UNIT_WEIGHT_WINDOWS_KN_M3` range table tagged `[POLICY: sampling-only]`.

## Registered - integer window copies

| window | where | relation |
|---|---|---|
| `m = random.randint(20, 30)` | `variables_control_charts.py:89` (`template_xbar_r_control_limits`) and `:445` (`template_chart_pair_selection`) | the exact bounds of `SPC_NUM_SUBGROUPS` (20, 30), drawn inline by two templates that do not read the table (its census consumers are the p- and c-chart templates). Below the C3.8 sweep's 3-s.f. rule; found by reading. Proposed: read the table (byte-identical while it stands) |
