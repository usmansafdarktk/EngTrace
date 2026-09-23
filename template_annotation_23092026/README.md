# template_annotation_23092026 — certifying the 150 templates, again

Started 2026-09-23. The record of the template certification pipeline for the next submission:
three layers, ordered by cost and by what each can actually detect. The analysis that led here
is summarised in [D-092](../docs/re-implementation-sep/DECISIONS.md) and
[D-093](../docs/re-implementation-sep/DECISIONS.md); the December 2025 pipeline it replaces is
described in the paper's section 3.3 and Appendix K.

| Layer | What | Cost | State |
|---|---|---|---|
| 0 · deterministic gate | T1 closure, T3 determinism, T4 contract, T8 emission at 500 seeds; T5, T7 advisory; a register of four line classes the check cannot read | free | **green: 150 of 150 pass** (`layer0/gate_report.md`); 54 templates edited over three closure rounds (`closure_fixes.md`) and 24 after the screen (`../screen/pass1_fixes.md`); 68 moved instances (`item_pool_impact.md`); residual ties only in 11 templates under 2% (`tie_census.md`), to be censused on the frozen pool |
| 1 · LLM screen | Appendix H's prompt, verbatim; three non-suite judges; two passes, hard-capped | $4.15 pass 1, $0.70 pass 2 | **both passes done**: pass 1 (2026-09-23) 126 pass, 13 controversial, 11 critical failure, AC1 on the flag 0.84; the 24 flags verified and fixed (`screen/pass1_fixes.md`); pass 2 (2026-09-24, 24 templates re-judged, 378 rows carried on unchanged prompts) **147 pass, 3 controversial, 0 critical**, AC1 0.93 (`screen/pass2/stats.md`) |
| 2 · human certification | own-branch experts, three per template, with a hand check, four planted defects per branch, timestamps (D-095) | expert time, about 4 h each | **kit built** (`layer2/`): guide, app, workbook route, task builder with verified plants, per-expert bundles, scorer; awaiting roster and dates |

## Pass 1 of the screen, read (as written before the fixes; the outcome is in `screen/pass1_fixes.md`)

All numbers are in `screen/pass1/stats.md`; every judge's sentence per flagged template is in
`screen/pass1/flagged.md`. What the 24 flags were about, grouped, as read on 2026-09-23 before any
claim had been verified; the verification and the fixes followed on 2026-09-24:

- **Presentation defects the gate does not read** (fix before Layer 2): float artefacts printed into
  text (`55.00000000000001 %`, `57.99999999999999%`, `-0.000e+00`) in `batch_moles_vs_conversion`,
  `pfr_volume_changing_rate`, `gauss_law_symmetric`; placeholder species names (`Entity 3`,
  `Species VII`) in `gas_phase_concentration`; `a applicant` in `server_configuration_selection`;
  a duplicated wave equation in `wave_equation_interpretation`; nonsensical material pairings
  (`bronze wooden block`) in `basic_buoyant_force`, the same flag the December Tribunal raised;
  every fluid labelled non-Newtonian in `power_law_fluid_shear`.
- **Physical-range concerns** (an expert's call): Pitzer correlation sampled beyond its accepted
  reduced-pressure range; two-term virial at extreme states; Hooke's law applied to concrete in
  tension and glass shafts at high torque; floating bodies whose stability is assumed, not checked;
  unrealistic velocity fields in `fluid_particle_acceleration`.
- **Rounding-chain claims** in six templates the gate passes (`average_energy_mqam`,
  `cd_dc_system_analysis`, `undamped_natural_frequency_translational`, `ideal_gas_volume`,
  `wave_equation_interpretation`, `truss_method_of_joints`): the judge says a displayed multiplication
  does not equal the printed result. These lines are ones T1 cannot parse (LaTeX `\times`, symbols in
  operand position), so they are candidates for the closure register or a fourth round, to be checked
  by hand.
- **One claimed logic error**: `wave_equation_interpretation`'s propagation-direction convention
  (one judge of three). Needs a human.
- **Judge artefacts to discount**: two judges "see" an undefined `signed_term` helper because the prompt
  shows the function source without the module's imports; that is the published prompt's design, not
  a template defect.

Judge behaviour worth knowing for pass 2: Grok flags 15% of templates, MiniMax 13%, MiMo 4%; on the
1-5 scores all three agree exactly on 86% of physical-plausibility ratings and 61% of pedagogical ones.

## Pass 2, read

Every judge claim from pass 1 was verified before any edit (`screen/pass1_fixes.md`: 34 of 45
confirmed and fixed, 9 rejected with evidence, 2 prompt artefacts). Pass 2 then re-judged only the
24 templates whose prompt text changed and carried the 378 pass-1 rows for unchanged prompts forward,
marked as carried. Of the 24: 21 now pass and 3 stay controversial, each on one judge of three:

- `truss_method_of_joints`: Grok repeats that `3.0/4.243 = 0.7070` is false; it is the correctly rounded
  quotient (0.70705 × 4.243 > 3), checked in `pass1_fixes.md`.
- `undamped_natural_frequency_translational`: MiniMax reads `1/0.82` shown as `1.22` and `1/2.681` as
  `0.373` as inconsistencies; both are the correct two- and three-decimal displays.
- `impulse_response_from_lccde`: MiniMax objects to an explicit unit coefficient (`1x[n-1]`, `1*delta[0]`)
  and a `1 - 8` written for `b1 + (-a1)·h[0]`; presentation taste, and the arithmetic is right.

Under the paper's rule a controversial template routes to manual review, which is what Layer 2 is,
so these three go to the experts as they are rather than being iterated against the judges.

## Layout

```
template_annotation_23092026/
  README.md              this file
  layer0/
    gate.py              runs the gate corpus-wide, applies the register, writes the report
    check_limits.json    the register: template + line pattern + reason for every excused line
    gate_report.md/.json the gate's outcome at HEAD (regenerate with gate.py)
    tie_census.py        exact-decimal census of half-way display ties T1 cannot count reliably
    tie_census.md/.json  its output (informational, gates nothing)
    closure_fixes.md     every template edited for closure: cause, fix, what moved, sign-offs
    item_pool_impact.md  what the edits moved, measured by c3_instance_dump.py --diff
  screen/
    run_screen.py        the screening runner (--check, --dry-run, --pass N, --status, --smoke)
    analyze_screen.py    consensus (paper rule), sigma-max, Gwet AC1 / Fleiss kappa, cost, served ids
    pass1/, pass2/       config.json + replies.jsonl + summary.csv + flagged.md + stats.md (committed)
    pass1_fixes.md       every pass-1 judge claim: confirmed / rejected, the fix, what moved, residuals
  layer2/
    README.md            the protocol and the order of operations
    guide.md, make_guide_pdf.py, EngTrace-certification-guide.pdf   what the experts read
    plants/              CONTRACT.md and one <branch>.py of four planted defects each
    build_tasks.py       verifies plants, precomputes instances, builds pool/assignment/keyfile (tasks/, ignored)
    app.py               the reviewer's app: hand check, solution, scores, verdict, timestamps
    workbooks.py         the file route: export, check, import
    package.py           one zip per expert (dist/, ignored); refuses to ship the keyfile
    simulate.py          synthetic labels to prove the pipeline; never results
    score.py             RESULTS.md: plant detection, hand checks, agreement, panel FPR, dwell, fix list
```

## Order of operations

```bash
python -m template_annotation_23092026.layer0.gate                 # free; must be green before anything else
python -m template_annotation_23092026.layer0.tie_census           # free; residual ties, informational
python -m template_annotation_23092026.screen.run_screen --check   # ~$0.01: keys, model ids, live prices
python -m template_annotation_23092026.screen.run_screen --pass 1 --dry-run   # free: prompts, tokens, estimate
python -m template_annotation_23092026.screen.run_screen --pass 1  # PAID (~$4.60); needs approval; resumable
python -m template_annotation_23092026.screen.analyze_screen --pass 1
```

Pass 2 is the same commands with `--pass 2 --carry-from 1`, which carries pass-1 rows forward where a
template's prompt is byte-identical and judges only what changed. `--pass 3` is refused; if the
experts change a template, re-judge that template alone (a few cents) rather than run a new pass.

## Conventions

- **Every number in a report here comes from a committed script**; the agents' per-template
  measurements in `closure_fixes.md` are labelled as such and the corpus-wide figures are the
  scripts'.
- **The register is the residual list.** An excused line names its template, its pattern and the
  reason; a growing register is a warning and the gate report prints it in full.
- **A pass is a pass.** Screening rows resume only within their own pass directory, and a pass
  refuses to continue if a template's prompt hash has changed since it started.
- **No paid call without the user's approval.** Dry runs and checks are free; the runner says
  which mode bills.
