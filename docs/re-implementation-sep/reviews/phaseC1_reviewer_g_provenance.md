# Phase C1 — Reviewer G (Provenance) report

- **Frozen ref:** `0f4b8cfcbdeec81603a8a543690b684314058845`
- **Time on the mandatory task:** 37 minutes, inside the 45-minute box. The optional part was not attempted.
- **Method:**
  - Template and constants source was read with `git show <ref>:<path>`.
  - Templates were run from the detached worktree `wt_g`, using in-memory patches of module globals only.
  - Scratch scripts and their outputs are in `scratchpad/g_scratch/`: `g_stated.py`/`g_stated.out` and `g_more.py`/`g_more.out`.
  - `phaseC1_census.md` was used only to pick which tables to sample and to compare classes. Every verdict below comes from template source plus template runs.
  - I called census.py's `probe`, `field_verdict` and `classify` directly for one table only (G-2), to confirm a mechanism I had first seen in the source. I did not re-run the census.

Reproduction preamble for every command below (PowerShell):

```
$env:PYTHONIOENCODING='utf-8'; Set-Location "C:\Users\AYESHA~1.GUL\AppData\Local\Temp\claude\c--Users-ayesha-gull01-EngTrace\19238be5-0710-45c5-9c79-55532f3f202b\scratchpad\wt_g"
$G = "C:\Users\AYESHA~1.GUL\AppData\Local\Temp\claude\c--Users-ayesha-gull01-EngTrace\19238be5-0710-45c5-9c79-55532f3f202b\scratchpad\g_scratch"
python "$G\g_stated.py"; python "$G\g_more.py"
```

---

## 1. Verdict

**PASS WITH FINDINGS.**
- None of the 14 PLAUSIBILITY tables I sampled is weaker than its evidence.
- One UNCONSUMED table is in fact consumed, through a copied literal: `SCS_IA_RATIO`, a `standard`, which should be CITATION (G-1).
- One flaw in how the census measures consumers could hide a weaker class in future, though it hides none today (G-2).
- Two lesser gaps: G-3 (PLAUSIBLE) and G-4 (CONFIRMED, no sourcing impact).

## 2. Independent re-derivation

### Part 1 — PLAUSIBILITY tables (14 checked; census class shown in brackets)

For each table I asked two things:
- **(a)** Is it a window, or a single value about a named entity?
- **(b)** Does every consumer state the drawn value in the **question**, or read the table only in an assert or screen?

**Finding consumers.** I found consumers my own way: a regex over each template function's source (`consumers_in_source` in `g_stated.py`). I added module-level reads by reading the source by hand. **Checking the questions.** "Stated k/40" means a regex pulled the value out of the question text over seeds 0–39, and the value lay inside the table window.

| # | table | how the census measured it | my consumers | (a) | (b) — what I ran | my class | agree? |
|---|---|---|---|---|---|---|---|
| 1 | civil `SPECIFIC_GRAVITY_RANGES` | ALL-RESTATED | 5: `phase_relations_degree_of_saturation`, `relative_density_of_sand`, `borrow_pit_fill_volume`, `effective_stress_profile`, `upward_seepage_quick_condition` | window per soil class; the spec uses "a sand's specific gravity" as its own example of a range (but see §5) | Gs stated 40/40 and in the soil's window 40/40 for **each** of the 5; the soil label is stated as well | PLAUSIBILITY | yes |
| 2 | civil `PERMEABILITY_RANGES_CM_S` | NO-EFFECT | 1: `constant_head_permeability` (only in the assert at `permeability_seepage_effective_stress.py:67-68`) | window | **guard only**: window widened to (0, 1e9) → items identical 40/40; window moved off the draws → AssertionError 40/40. k is the requested answer and is never stated; the draws come from an inline [0.02, 0.09] | PLAUSIBILITY | yes |
| 3 | civil `FRICTION_ANGLE_RANGES_DEG` | ALL-RESTATED | 1: `infinite_slope_factor_of_safety` | window | phi and the soil descriptor stated 40/40, phi in window 40/40 | PLAUSIBILITY | yes |
| 4 | industrial `SERVICE_LEVELS` | **NO-EFFECT** | 1: `safety_stock_reorder_point` | a set of policy levels, no named entity | level stated 40/40. **Not a guard**: the drawn level is looked up in `Z_QUANTILES`, so a nudge raises KeyError 40/40 (see G-2) | PLAUSIBILITY | class yes; census's reason ("guard") no |
| 5 | industrial `P_CHART_PBAR` | INDETERMINATE | `p_chart_limits_floor`: import-time plan builder, the D window, and asserts | window | D, m, n stated 40/40; D/(mn) inside window 40/40 | PLAUSIBILITY | yes |
| 6 | industrial `P_CHART_SUBGROUP_N` | NO-EFFECT | same template: import-time `_t27_admissible_plans` and an assert | window | n stated and in window 40/40 | PLAUSIBILITY | yes |
| 7 | industrial `C_CHART_CBAR` | NO-EFFECT | `c_chart_revision`: `cb_hi` in a screen, `cb_lo` only in an assert (`attributes_control_charts.py:491,497,619`) | window | screen/guard; the question states T and m, and c-bar = T/m is derived from them (seed 3 read) | PLAUSIBILITY | yes |
| 8 | industrial `COATING_METROLOGY_FLOOR` | INDETERMINATE | `sigma_reduction_for_cpk`: `if sigma < FLOOR: continue` (`process_capability.py:411`) | a scalar one-sided bound, not a named-entity value (borderline, §5) | screen; sigma stated in the question ("sigma = 2.28 microns", seed 5) | PLAUSIBILITY | yes |
| 9 | industrial `SPC_NUM_SUBGROUPS` | NO-EFFECT | `p_chart_limits_floor` (import-time plans), `c_chart_revision` (assert only; m is an inline `choice([20,25,25,25])`) | window | m stated 40/40 (p chart) and in the c-chart question | PLAUSIBILITY | yes |
| 10 | industrial `QUEUE_SCENARIOS` | ALL-RESTATED | `mm1_time_in_system` (call time), `mmc_waiting_time` (import-time `_MMC_COMBOS`) | window per scenario | rates stated 40/40 (M/M/1); lambda and c stated 40/40 (M/M/c) | PLAUSIBILITY | yes |
| 11 | industrial `NEWSVENDOR_ITEMS` | ALL-RESTATED | `newsvendor_normal_demand` | window | c, p, s stated 40/40 | PLAUSIBILITY | yes |
| 12 | industrial `COMPONENT_RELIABILITY_CLASSES` | ALL-RESTATED | `system_reliability_topology` | window | every component reliability stated, 40/40 in [0.90, 0.99] | PLAUSIBILITY | yes |
| 13 | industrial `SPC_CHARACTERISTICS` | ALL-RESTATED | 4: `chart_pair_selection`, `cp_cpk_from_specs`, `xbar_known_sigma_classification`, `xbar_r_control_limits` (only `["target"]` is read) | window | grand mean / mu0 stated 40/40 in each of the 4 | PLAUSIBILITY | yes |
| 14 | electrical `AMPLITUDE_RANGE` | ALL-RESTATED | 3: `cd_dc_system_analysis`, `continuous_to_discrete_conversion`, `nyquist_rate_determination` | window | amplitude stated 40/40 and in [1, 50] in each of the 3 | PLAUSIBILITY | yes |

All seven tables the census measured NO-EFFECT or INDETERMINATE were checked from source (rows 2, 4–9). Only `PERMEABILITY_RANGES_CM_S` and `C_CHART_CBAR` are true guards. `COATING_METROLOGY_FLOOR` is a screen. In the other four, the drawn value is stated. Either way PLAUSIBILITY is correct, but the census's stated reason is wrong for rows 4, 5, 6 and 9 (see G-2).

Commands: `python "$G\g_stated.py"` (rows 1–12, plus the permeability and service-level patch tests) and `python "$G\g_more.py"` (rows 13–14).

### Part 2 — UNCONSUMED tables (10 checked individually; all 32 swept by name)

**Sweeps run (independent of the census detector):**
1. `git grep -w <name>` over the whole tree except docs, for all 32 names.
2. `git grep` over `data/templates` for dynamic access: `getattr(`, `globals()`, `vars(`, `import *`, `__dict__`, `importlib`, `eval(`/`exec(`, and module-alias imports. **Zero hits.**
3. A **literal-copy sweep**: the values of the UNCONSUMED tables, grepped in templates outside `constants.py`. The census has no detector for this, so I rate it the likeliest blind spot.

| table | why I chose it (possible blind spot) | evidence | my class | agree? |
|---|---|---|---|---|
| civil `SCS_IA_RATIO` (`standard`) | value could be copied as a literal | `hydrology.py:186` `Ia = round(0.2 * S, 3)`, `:200` `0.8 * S`. Question states "Ia = 0.2S" 40/40. Table set to 0.3 and module reloaded → items identical 40/40 | **CITATION** (consumed, `standard`) | **no → G-1** |
| industrial `XBAR_R_SUBGROUP_N`, `XBAR_S_SUBGROUP_N` (`standard`) | the rule could be inlined | `chart_pair_selection` inlines n ∈ [4,8] / [13,20] and grades the pair choice; the question does not state the n > 10–12 rule. The table values themselves are not copied | UNCONSUMED by the letter; warrant hidden → G-3 | yes, with G-3 |
| mechanical `ATMOSPHERIC_PRESSURE_KPA` (`defined`) | imported but not read; value could be copied | imported in `fluid_statics.py:3` and never read; `:36` `P_atm_kpa = 101.325` is a copy, **stated** in the question | would be DEFINITION if counted as consumed; no sourcing either way | yes (no impact) |
| chemical `REAL_FLUID_DATA` (`property`) | imported by a template module | imported in `volumetric_properties_pure_fluids.py:4`, never read; `two_phase_specific_volume` invents V_l/V_v inline (the known C3.3 defect) | UNCONSUMED | yes |
| chemical `AIR_COMPOSITION` | combustion code could copy the N2/O2 ratio | no `3.76`, `0.79`, `0.78084` or `0.20946` in templates | UNCONSUMED | yes |
| civil `CV_RANGES_M2_YR` | the brief's example | only its definition; the consolidation template works by similitude and never draws cv | UNCONSUMED | yes |
| civil `DAS_NATURAL_STATE_SOILS` | templates cite "Das Table 3.1 anchors" | the anchors are an inline `e_ranges` dict in `phase_and_index_properties.py:45-49`, not a read (§5) | UNCONSUMED | yes |
| civil `STEEL_FY_KSI` / `STEEL_FY_MPA` | a name could be built from a unit suffix | no getattr-type access at all; no `A992`, `A36` or `Fy` literals in templates | UNCONSUMED | yes |
| industrial `FINITE_CAPACITY_K` | named in template prose | `queueing_systems.py:479,483,542` are comments; K ∈ [4,6] is inline | UNCONSUMED (`range`) | yes |
| industrial `FAILURE_RATE_PER_HR` | named in template prose | `system_reliability.py:191` is a comment; the window is inline | UNCONSUMED (`range`) | yes |

The name sweep also checked these, with no read found: `LEAD_TIME_WEEKS`, `LINE_TASK_TIME_S` and `DISCOUNT_*` (comments only), plus `GRAVITY_FT_S2`, `UNIT_WEIGHT_WATER_PCF`, `WATER_DENSITY_KG_M3`, `WATER_KINEMATIC_VISCOSITY_M2_S`, `CONCRETE_*` and `LIVE_LOADS_*`. The literal sweep found no copies of `32.2`, `62.4`, `998.2`, `1.004e-6`, `57000`, `4700` or live-load values.

---

## 3. Findings

### G-1 — CONFIRMED — `civil_engineering.SCS_IA_RATIO` is UNCONSUMED in the census but consumed through a copied literal; it should be CITATION

- **Claim.** `template_scs_curve_number_runoff` uses the TR-55 initial-abstraction ratio as a hardcoded `0.2` (`hydrology.py:186`, and `0.8 * S` at `:200`) instead of reading `SCS_IA_RATIO`. It states the ratio in the question 40/40. The table is `@kind: standard`, and C1.3 makes every consumed `standard` CITATION.
- **Why the census missed it.** Its consumer detector works on names. A number copied into a template does not mention the name.
- **Reproduction:**
  ```
  git -C C:\Users\ayesha.gull01\EngTrace grep -n -E "0\.2 ?\* ?S|0\.8 ?\* ?S|SCS_IA_RATIO" 0f4b8cfcbdeec81603a8a543690b684314058845 -- data/templates
  python "$G\g_more.py"
  ```
  The script output has three lines for this finding:
  - "SCS_IA_RATIO readers: []"
  - "'Ia = 0.2S' in question 40/40"
  - "SCS_IA_RATIO set to 0.3 and module reloaded -> identical items 40/40"
- **Impact.**
  - Sourcing work hidden: 1 literal. The table already carries `[ON-DISK]` with no locator, and an UNCONSUMED table gets no C3 sourcing.
  - A worse effect: if C3 corrects the table, the item will not change, because it does not read the table.
  - The pattern generalises. Name-based consumer detection cannot see copies, so the UNCONSUMED count (32) is an upper bound.
- **Fix.** Make the template read `SCS_IA_RATIO`. Or add a literal-copy check to the census (§5).

### G-2 — CONFIRMED (mechanism) / PLAUSIBLE (latent escape) — the census counts "every consumer crashed" as NO-EFFECT, and C1.3 then treats that as "consumed as a guard"

- **Claim.** In `census.py`:
  - `field_verdict` returns `ERROR` when a nudged table makes the template raise.
  - The table-level rollup has no case for ERROR, so an all-ERROR table falls through to `NO-EFFECT`.
  - `classify('range', 'NO-EFFECT')` returns PLAUSIBILITY with the reason "consumed as a guard or screen".
- **Instance.** `SERVICE_LEVELS` is measured exactly this way. It is not a guard: it is a sampling set whose drawn level is used as a dict key into `Z_QUANTILES`, so any nudge raises KeyError.
- **Why the class is still right today.** The level is stated in the question 40/40.
- **Why it matters.** The same path would classify a `range` table as PLAUSIBILITY even if it were read **hidden** through a lookup that raises when nudged. By C1.3 that table should be CITATION.
- **Reproduction:**
  ```
  python "$G\g_more.py"
  ```
  The script output has two lines for this finding:
  - `[PROBE] SERVICE_LEVELS field='[*]': {... 'ERROR': 10 ...} -> field_verdict ERROR`
  - `classify('range','NO-EFFECT') -> ('PLAUSIBILITY', 'NO-EFFECT: consumed as a guard or screen')`

  My own patch, independent of the census: `python "$G\g_stated.py"` prints `[SL] nudged levels -> errors 40/40; first: ERROR KeyError: 0.99099`.
- **Impact.**
  - Now: none. I checked all 7 guard-measured PLAUSIBILITY tables from source (rows 2, 4–9), and none reads its value hidden.
  - Future: a crash-only measurement licenses PLAUSIBILITY with no human check.
  - The census's stated reason is also wrong for `SERVICE_LEVELS`, `P_CHART_PBAR`, `P_CHART_SUBGROUP_N` and `SPC_NUM_SUBGROUPS`: all four are sampling windows whose drawn values are stated, not guards.
- **Fix.** Give the rollup a separate `ERROR` measurement. For a `range` table, make both `ERROR` and `INDETERMINATE` require a recorded manual guard/stated verdict before PLAUSIBILITY is granted.

### G-3 — PLAUSIBLE — a hidden, graded decision rests on a standard that sits only in UNCONSUMED tables (`XBAR_R_SUBGROUP_N` / `XBAR_S_SUBGROUP_N`)

- **Claim.** `template_chart_pair_selection` grades the choice between X-bar/R and X-bar/s. The question gives n but not the threshold (read at seed 3: "Select the appropriate Shewhart variables chart pair for this subgroup size…"). The template encodes Montgomery §6.3's "n > 10 or 12" rule as inline sampling windows, `randint(4, 8)` / `randint(13, 20)` (`variables_control_charts.py` about line 447, docstring lines 402–410). The two `standard` tables that record this rule are never read.
- **Why only PLAUSIBLE.** The table *values* (2–10, 11–25) are not copied; the template deliberately avoids n = 11–12. So UNCONSUMED is correct by the letter of the rule. But the one table-level home of the warrant for a hidden, answer-deciding rule gets no sourcing.
- **Reproduction:**
  ```
  git -C C:\Users\ayesha.gull01\EngTrace show 0f4b8cfcbdeec81603a8a543690b684314058845:data/templates/branches/industrial_engineering/quality_and_reliability_control/variables_control_charts.py | sed -n 395,450p
  ```
  The question text is printed by `python "$G\g_stated.py"` (block `--- template_chart_pair_selection seed=3`).
- **Impact.** 2 tables, 4 literals, one Montgomery citation. Low, but it is the kind of rule C1.3 exists to source.

### G-4 — CONFIRMED (no sourcing impact) — `electrical_engineering.PHASE_RANGE_RAD` is consumed but outside the census and the metadata check

- **Claim.**
  - `PHASE_RANGE_RAD = (-math.pi, math.pi)` (`constants.py:73`) is read by `template_continuous_to_discrete_conversion` at `continuous_time_signals.py:133`.
  - `static_tables` drops any table with zero numeric literals (`census.py:143-144`), and `test_table_metadata` M1/M2 use that same table list. So the table has no `@kind`, no `@units` and no class.
  - The window is exact (±π), so there is nothing to source.
  - What it shows: a table built entirely from named constants silently leaves C1's universe.
- **Reproduction:**
  ```
  git -C C:\Users\ayesha.gull01\EngTrace grep -n -w PHASE_RANGE_RAD 0f4b8cfcbdeec81603a8a543690b684314058845 -- data/templates
  git -C C:\Users\ayesha.gull01\EngTrace show 0f4b8cfcbdeec81603a8a543690b684314058845:tests/constants_integrity/census.py | sed -n 127,150p
  ```
  The table is also absent from the census doc's electrical rows (13 tables).
- **Caveat.** My regex matched the degree-phase branch in 24/40 items. It matched no radian phase in the other 16, so I did not confirm that the radian phase is stated.

---

## 4. Falsification attempts that failed

- **SPECIFIC_GRAVITY_RANGES (re-declared `property` → `range`).**
  - I looked for a consumer beyond the census's 5 (function-source scan plus `git grep`) and found none.
  - I looked for a consumer that uses Gs without stating it. `effective_stress_profile` uses Gs inside a hidden coupling (`e_above`), but still states Gs.
  - Gs is stated and in window 40/40 in all five.
- **FRICTION_ANGLE_RANGES_DEG.** `terzaghi_strip_footing_bearing` looked like a hidden consumer: its comment cites the table to justify phi = 35. It does not read the table; phi is inline. The only reader states phi 40/40. (I did not run Terzaghi to confirm that its inline phi is stated.)
- **PERMEABILITY_RANGES_CM_S.** I tried to show the table shapes the answer. It does not: widening it leaves items identical, and moving it only fires the assert.
- **Import-time consumption** (`_T27_PLANS`, `_MMC_COMBOS`). I suspected the probe could not see tables read at import. The `Patch`/`probe` docstring and code (`census.py:404-491`) reload template modules, so it can. The drawn m, n, D, lambda and c are stated 40/40.
- **COATING_METROLOGY_FLOOR as a hidden named-entity fact.** It appears in no question. It only rejects sigma draws, and sigma is stated.
- **UNCONSUMED blind spots.**
  - Dynamic name access: zero hits.
  - Constants-level aliasing: `git grep -w` finds each UNCONSUMED name only at its definition, in comments, or on an import line.
  - Imports without reads: `REAL_FLUID_DATA` and `ATMOSPHERIC_PRESSURE_KPA` are imported and never read.
  - Literal copies of civil/chemical property values (`32.2`, `62.4`, `998.2`, `1.004e-6`, `57000`, `4700`, air fractions, steel grades, live loads, cv): none.
  - The only literal copies found are `101.325` (defined, stated, no sourcing either way) and `0.2` (G-1).

## 5. Further probing and improvements

- **Where the evidence is thinnest.**
  - My literal-copy sweep covered about 12 values by hand.
  - C3 should run it mechanically: for every numeric literal of every table (in any class), grep templates outside `constants.py` for the same value, rounded to the table's own digits, and triage the hits. G-1 came from that idea applied to only a few values.
  - Such a sweep would also catch CITATION tables whose values drift into inline copies. That is a correctness risk, not only a classification one.
- **Inline windows outside the census's reach.** C1.3 classifies *tables*, but named-entity facts also enter items from inline literals:
  - `two_phase_specific_volume` states invented saturation volumes for real substances.
  - `phase_and_index_properties.py` has an `e_ranges` dict "anchored to Das Table 3.1".
  - Terzaghi's inline `phi = 35` for "dense sand".
  - The permeability template's inline [0.02, 0.09] for coarse sand.

  None of these can be classified or sourced. I recommend C3 add an "inline-window register", or require that such windows move into tagged tables.
- **`SPECIFIC_GRAVITY_RANGES` (borderline under (a)).**
  - The spec's definition licenses the `range` kind, and every consumer states Gs.
  - But the `sand` row (2.65, 2.67) at 2 dp gives only {2.65, 2.66, 2.67}, which is close to quartz's single Gs.
  - PLAUSIBILITY is defensible because the value is stated. It still deserves a C3 plausibility check against Das before it is treated as "no sourcing".
- **`COATING_METROLOGY_FLOOR`'s kind.** A scalar screening threshold is neither "a window a value is drawn from" nor `validity` ("a bound for another table"). The vocabulary has no kind for a screen threshold. The spec should add one or say which kind covers it.
- **The spec's guard clause (C1.3).** "Consumed only as a guard" is inferred from a dynamic measurement that cannot tell a guard from a crash (G-2). I suggest defining "guard" statically — an `assert`, or an `if … continue/raise` over already-drawn values — and using the probe only to corroborate.
- **Consumer counts.** My function-source scan found 1 call-time reader for `QUEUE_SCENARIOS` (M/M/1), plus one import-time reader (M/M/c); the census reports 4. The census probably attributes module-level reads to every template in the module. It does not change any class, and I did not run it down. The same may apply to other tables' "consumers" counts.
- **Seeds.** 40 seeds cover every branch I checked. In `continuous_to_discrete_conversion`, 16 of 40 items did not match my phase regex (G-4 caveat); a closer look there is cheap.
- **Optional part** (C1.7 pilot against CODATA; resolver miss): not attempted in this filing.
