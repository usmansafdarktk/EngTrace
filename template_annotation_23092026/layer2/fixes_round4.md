# Round 4 fixes: two chemical templates widened so each can supply 15 distinct problems

2026-09-28, D-115. Not an expert rejection this time: the diversity analysis of the frozen pool
(`full_run_28092026/DIVERSITY.md`) found two templates that could not produce 15 distinct questions
from any seed, `heat_of_reaction_formation` (4) and `adiabatic_flame_temperature` (11), so the
freeze filled them with repeats (D-114). The owner decided to fix both and send them back to the
chemical experts. Every number below is printed by `round4_checks.py` (`round4_checks.md`) or by the
layer-0 scripts named.

## heat_of_reaction_formation (Easy)

**Cause.** It drew from `REACTIONS`, four fixed reactions, and samples nothing else, so it had four
problems.

**Fix.** A new table, `HESS_REACTIONS` in `chemical_engineering/constants.py`: the four `REACTIONS`
rows plus 21, read only by this template. The question, the reasoning and the level are unchanged.
The 21 use only species already priced in `HEATS_OF_FORMATION`, so no new constant enters: ten
combustion reactions (ethane, butane, acetylene, octane vapour, benzene, both phases of methanol and
ethanol, carbon monoxide, hydrogen, and propane to water vapour), the water-gas shift, methanol from
synthesis gas and from carbon dioxide, ammonia synthesis, nitric-oxide oxidation, ammonia burned to
nitrogen, dry reforming and partial oxidation of methane, and acetylene trimerised to benzene.
N2 + O2 -> 2NO was left out: the constants suite flags any reaction with N2 and O2 as reactants
outside the theoretical-air ratio.

**A first attempt widened `REACTIONS` itself, and the gate caught it.** Three stoichiometry templates
also draw from `REACTIONS` (`batch_moles_vs_conversion`, `flow_system_molar_flow_rates`,
`limiting_reactant`); with the wider table each failed to generate on 7 of 200 seeds, and all three
would have changed instances their experts had certified. `REACTIONS` was restored from git, byte for
byte, and the new rows went into their own table.

**Checked.** 25 reactions, 25 distinct equations; every one balances atom by atom, every species is
priced, and every printed equation says exactly what its dicts say. The constants suite now counts
atoms over `HESS_REACTIONS` too (324 checks, all pass).

## adiabatic_flame_temperature (Advanced)

**Cause.** Eleven fuels, always with theoretical air, and nothing else sampled: eleven problems.

**Fix.** The template samples excess air, 0 to 100 percent in whole percent: 1,111 fuel-and-level
cases. The excess oxygen leaves with the products, the air's nitrogen passes through at 3.76 per O2 as
in the table, and ammonia's own nitrogen is kept apart from the air's. Coefficients are exact decimals,
so every number Step 1 prints is the value the energy balance uses; Step 1 now shows the O2 supplied,
the N2 supplied and the excess O2 before the equation. At 0 percent the question and Step 1 are the
theoretical-air text word for word. `COMBUSTION_REACTIONS` is unchanged.

**Every claim in the docstring and in Step 6 was measured at theoretical air, so each was measured
again over all 1,111 cases:**

| | at theoretical air, as documented | over 0-100% excess air |
|---|---|---|
| six passes vs the exact fixed point | within 0.035 K | within 0.102 K |
| contraction of the iteration | about 0.09 to 0.135 | 0.089 to 0.135 |
| flame temperature range | up to 2908 K | 1430 to 2908 K, inside the 3000 K ceiling |
| answer vs a direct NIST Shomate solve | 1.93 K | 1.90 K at 0%, 2.41 K worst |
| rounding margin | 0.116 K | 0.0007 K: the margin does not carry over |

This script's NIST solve reproduces the documented theoretical-air figure within 0.03 K, which is the
check on the method.

Two changes follow from the table:

- **Step 6 says 2.5 K where it said 2 K.** The measured worst is 2.41 K.
- **A guard redraws the level in 20 of the 1,111 cases.** In 13 the six displayed passes and the exact
  fixed point round to different kelvins, and in 12 the last pass displays an exact half kelvin; there
  the answer's last digit would come from the iteration count or the rounding convention, not from the
  model. The level alone is redrawn, so the fuel shares stay flat (D-045). Over 3,000 seeds the
  template never emits a guarded case, and its answer equals this script's replica of the iteration on
  all 3,000.

## Evidence after the edits

| check | result |
|---|---|
| Layer 0 gate, 500 seeds | 150 of 150 pass; T5 advisory 60, as before |
| tie census, 500 seeds | 12 templates with a tie, as before; neither of these two |
| Markdown scan, 200 seeds | 65 templates with a lossy construct, as before; neither of these two |
| constants suite; citations | all pass |
| certification | 148 of 150: these two no longer reproduce the instances reviewed in round 1 |

## Round 4

Built with `build_tasks --only template_heat_of_reaction_formation,template_adiabatic_flame_temperature
--round 4` and `make_kits --round 4`: the two templates for the three chemical experts, no planted
defects, and hand-check instances none of them has seen. Send each expert the shared files in
`dist_round4/` and their `kit_che-N/`. When the labels return, `certification.py` with the round-4
folder added decides whether both are certified again.
