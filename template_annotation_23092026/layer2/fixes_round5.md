# Round 5 fixes: two chemical templates whose questions did not pin the answer

2026-10-06. This round did not start with a certification rejection. The defect came from the three
chemical experts' reading of evaluated responses on 2026-10-03, task B4 of the expert request
(`full_run_28092026/EXPERT_REQUEST.md`, D-185).

The two repaired templates went to certification round 5. Round 5 approved
`adiabatic_flame_temperature` by all three experts and `work_isothermal_virial` by two of three. The
virial template was changed once more, with a temperature cap on organic compounds, and round 6
approved it by all three. All 150 templates are certified again (`CERTIFICATION.md`).

Every number below is printed by `round5_checks.py` (`round5_checks.md`, the current code against
the pre-repair code at `03b1dd7`), by the layer-0 scripts named, or by the full run's scripts named.

## The defect, as the experts stated it

- `work_isothermal_virial`, three readers. Does the wording decide the reading, closed-system or
  flow work? "Either: the wording does not decide", 3 of 3. Does it decide the form of the
  truncated virial equation? The same answer, 3 of 3. Is there a single correct answer? "No",
  3 of 3. Which of the responses shown are correct? "All of them", 3 of 3.
- `adiabatic_flame_temperature`, three readers. Do standard heat-capacity sources agree within the
  0.2% answer tolerance? "No: standard sources differ by more than that" 2, "only if the same data
  source is used" 1. Is the method standard? "Yes" 2, "no" 1. Which of the responses shown are
  correct? "Some of them", 3 of 3.

## work_isothermal_virial (Advanced)

**What the gold computes, unchanged since screen pass 1.** The work done on the gas in a closed
system, W = -∫P dV, with the volume form of the two-term virial equation throughout,
Z = 1 + B/V, so P = RT(1/V + B/V²). B comes from the Pitzer correlation with Abbott's B0 and B1.

**The question now pins that reading** (round 5):
- the system and process: one mole in a closed system (a piston-cylinder device), compressed
  isothermally and mechanically reversibly;
- the quantity: the work done on the gas, W = -∫P dV from the initial to the final molar volume,
  positive for a compression;
- the form, written out: the virial equation truncated to two terms in its volume form,
  Z = PV/(RT) = 1 + B/V;
- B's source, by name: the Pitzer correlation B·Pc/(R·Tc) = B0 + ω·B1, with Abbott's equations for
  B0 and B1.

Step 1 of the solution names the same reading and form.

**The gold is now the stated reading's answer** (round 5). Re-solving each question from its own
text found a second defect. For helium, a few kelvin above its critical point, the work is tens of
J/mol, and the gold's chain of displayed values (R*T at 2 dp, the work terms at 0.01 L·bar) moved it
by up to 2% from the exact work for the stated data. With the pre-repair code, the exact work fell
outside the 0.2% tolerance of the gold at 13 of the gate's 500 seeds, every one of them helium;
every other substance stayed within 0.0016. A draw more than 0.1% (half the tolerance) from the
exact work is now redrawn, as a display tie is.

**Round 5's verdicts.** Two experts approved. One rejected on the objection first raised in round 1
(D-094, not adopted then): the Vr ≥ 2 validity filter keeps only Tr of about 2 to 3, so organic
compounds are compressed far above their decomposition range (n-pentane at 1,135 K up to 142 bar),
where no quasi-static compression of that species can run. The expert asked for a temperature cap of
about 700-750 K on organics, or for dropping the substances then left without a valid state. A minor
note added that the solution's unspaced `*` rendered as italics in the review app.

**Organic compounds are now capped at 700 K** (round 6, the owner's choice within the expert's
range). A draw of an organic compound above 700 K is redrawn. The cap is checked after all four draws
of an attempt, so an attempt consumes the generator as before, and an instance that was not a hot
organic keeps every number. 11 of the 16 organics can then never meet the validity condition and
leave the template: n-butane to n-octane, benzene, toluene, p-xylene, methanol, ethanol, acetone.
Methane, ethane, ethylene, propane, propylene and the 11 inorganic substances remain, 16 in all. The
minor note is taken too: the substituted lines of the solution print `·` and `^` where they printed
`*` and `**`, and no number changes. The cap and the organic list live beside the template
(`_VIRIAL_ORGANICS`, `_VIRIAL_ORGANIC_T_MAX`); `constants.py` is unchanged.

**What did not change.** The sampling ranges, the constants and the arithmetic. At 193 of the 500
seeds the instance is the same problem as before the repair, and its solution from Step 2 on is
identical but for the operator glyphs. The other 307 are redraws, each explained:
- 292 were an organic above 700 K;
- 15 lay more than 0.1% from their exact work (14 helium, 1 hydrogen).

## adiabatic_flame_temperature (Advanced)

**Fix** (round 5). The question now carries the data the energy balance consumes, and says to use it:
- the standard heats of formation at 298.15 K, in kJ/mol, of every species in the equation, with
  water as vapour;
- the Cp/R coefficients of every product exactly as `CP_PARAMS_COMBUSTION` holds them, each printed
  at the shortest display that re-parses to the float consumed (`_exact_spec`, D-037);
- the polynomial Cp/R = A + B·T + C·T² + D·T⁻², with T in K, the coefficients' units and
  R = 8.314 J/(mol·K);
- dry air at 3.76 mol N2 per mol O2;
- the products as an ideal-gas mixture of the species listed, with no dissociation;
- the answer reported to the nearest kelvin.

The species lines are a list, so the review app renders them as bullets and removes nothing. The
3,000 K validity ceiling of the fit is not printed, because the gold's second iteration evaluates Cp
beyond it on its way to the answer.

**What did not change.** The sampling (0-100% excess air), the solution and the answer line. The
solution is identical to the pre-repair code at 500 of 500 seeds, and each question is the
pre-repair question with the data block appended.

## Evidence on the current code

| check | result |
|---|---|
| Layer 0 gate, 500 seeds (`layer0.gate`) | 150 of 150 pass; on these two, T1 has no failures (4,000 and 4,982 checks) and T3, T4 and T8 pass; the register (1,494 lines) and the advisories (T5 60, T7 78) are as before |
| tie census, 500 seeds (`layer0.tie_census`) | 12 templates with a tie, as before; neither of these two |
| virial: organics above 700 K | 0 of 500 seeds |
| virial re-solved from its question (SI units, independent code) | within 0.1% of the gold at 500 of 500 seeds (worst 0.000994) |
| virial, the other readings the experts named | closed system with the pressure-explicit form outside the tolerance at 53 of 500 seeds; steady-flow shaft work at 469 of 500, under either form |
| flame re-solved from its question (stoichiometry from the formula, exact integral) | rounds to the gold kelvin at 500 of 500 seeds (worst 0.49 K) |
| flame data block | prints exactly the floats the code consumes at 500 of 500 seeds |
| the 30 frozen items, re-solved | 30 of 30 reach their gold |
| the round-6 kit's 5 instances, re-solved | 5 of 5 reach their gold |
| LLM screen | not run: its passes are tied to the corpus they started on and refuse once a prompt hash moves, and a third pass is refused by design |

## Round 5

Built 2026-10-06 with:
- `build_tasks --only template_work_isothermal_virial,template_adiabatic_flame_temperature --round 5`
- `make_kits --round 5`

The two templates went to the three chemical experts, with no planted defects, at instance seeds
2501-2505 and hand-check seed 2501.

**Outcome** (`RESULTS_round5.md`; `score.py --round 5`):
- `adiabatic_flame_temperature` approved by all three, and now certified.
- `work_isothermal_virial` approved by two of three; the rejection is above.
- Hand checks: 6 of 6 matched the template within 1%.

With round 5 alone, 149 of 150 templates were certified; the virial template waited for round 6.

## Round 6

Built 2026-10-06 with `build_tasks --only template_work_isothermal_virial --round 6` and
`make_kits --round 6`. The template goes to the same three chemical experts at seeds 2601-2605, with
hand-check seed 2601, which none of them has seen. Scored with:

```bash
python -m template_annotation_23092026.layer2.score --round 6 --labels <round-6 folder> --prev-labels <round-5 folder>
python -m template_annotation_23092026.layer2.certification --labels <round-1 folder> ... <round-6 folder>
```

**Outcome** (`RESULTS_round6.md`):
- `work_isothermal_virial` approved by all three experts, the round-5 rejection among them.
- Hand checks: 3 of 3 matched the template within 1%.
- `CERTIFICATION.md` over rounds 1-6: 150 of 150 certified. 126 were last reviewed in round 1, 16 in
  round 2, 5 in round 3, 1 in round 4, 1 in round 5 (`adiabatic_flame_temperature`) and 1 in round 6
  (`work_isothermal_virial`). Every template's current code regenerates the instances its experts
  were shown.

## What followed in the full run

- **The re-freeze.** `freeze.py --only-templates ... --write --record round5` drew both templates'
  items again from the same seed and selection rule. After round 6,
  `freeze.py --only-templates template_work_isothermal_virial --write --record round5 --amend` drew
  the virial's again, keeping the record's pre-round-5 hashes.
  - Outside the two templates, all 2,220 manifest rows are unchanged.
  - Inside them, 27 rows changed only their sha256.
  - Three virial indices left the selection (`#12`, `#19`, `#24`) and three came in (`#14`, `#15`,
    `#16`). The subsample positions keep their ids.
  - `freeze.py --verify` regenerates all 2,250 items byte for byte.
  - `diversity.json` changed for the virial template only: its 15 items span 9 question skeletons
    where they spanned 12, and its 500-draw reach covers 16 substances and 2 reasoning paths where it
    covered 27 and 3.
- **Open book.** `openbook.py --build --version 2` changed the six open-book items of these
  templates: their question changed, their reference equations did not.
- **Paraphrases.**
  - The subset items' writer attempts and verdicts written for earlier questions were archived
    (`repair_round5.py --archive-paraphrases`).
  - New writing for the current questions: virial 3 of 3 passed at the first attempt. Flame 0 of 3
    passed after three attempts: its data block must be copied verbatim, which fails the near-copy
    limit, and the writer altered its numbers.
  - The three virial pairs go to one chemical expert for the same-problem check.
