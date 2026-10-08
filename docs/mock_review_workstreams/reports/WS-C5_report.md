# WS-C5 report: statistics tables and the worked example

Done, 2026-10-07. Nothing is committed (none asked). No paid call was made. No `.tex` file under
`overleaf_source_04102026/` or `current_overleaf_project/` was edited; both scripts wrote only into a scratch copy
(`--out DIR`). No other stream's file was edited.

## Scripts and flags

| script | what it computes | output | check |
|---|---|---|---|
| `docs/appendix_statistics.py` | the Dataset Statistics numbers (unchanged) plus two generated blocks: `tab:area` (templates per area and level, grouped by domain, with branch subtotals, from `manifest.jsonl`'s `area` field) and `tab:levels_agreement` (from `template_annotation_23092026/levels/agreement.json`: Fleiss' kappa per branch and overall with 95% intervals, the linear- and quadratic-weighted forms, templates whose two-of-three majority equals the label or lies one or two steps from it or has no majority, the formula-stated yes count by labeled level) | `appendices/taxonomy_content.tex` under `--out DIR` | `--check` passes on the untouched tree (28 of 28 rows and numbers; the two blocks reported as "no markers yet"); `--selftest` passes |
| `full_run_28092026/worked_example.py` (new) | the worked scoring example `box:worked_example` from one main-store response: the question as posed, the gold milestones and how each was settled (matched, with the stated value and the unit factor c of eq. match; or the judge's ruling), the flagged calculation with its printed operands, printed result, recomputed value and judging place value, the stated answer against its target, then v(r) and m(r) | `appendices/worked_example.tex` under `--out DIR` (62 lines, every line under 100 characters) | `--selftest` re-runs the answer check, E3 and the arithmetic check on the trace through `score.py`'s functions and asserts they reproduce the stored label, matches, scales, flags, `score` and `e5_strict`; passes |

Flags of `docs/appendix_statistics.py`: `--check`, `--write`, `--out DIR`, `--agreement PATH`, `--selftest`.
Markers: `% BEGIN GENERATED <label> (docs/appendix_statistics.py --write)` ... `% END GENERATED <label>`, the
convention of `paper_results.py`. Under `--out DIR`, a block whose markers are absent is appended at the end of
the copied `taxonomy_content.tex` (so Phase 2 can move the markers); without `--out`, absent markers stop the run.
A second `--write` replaces the blocks in place. Until `agreement.json` exists (or holds no ratings), the
agreement block is a stand-in with empty cells and a caption starting "STAND-IN"; no simulated numbers are
printed anywhere. `--agreement PATH` reads another file, for a layout test with `score_levels.py --selftest`'s
output.

Flags of `worked_example.py`: `--out DIR`, `--write` (Phase 2, into `overleaf_source_04102026/appendices/`),
`--model`, `--item`, `--candidates`, `--reader-gap`, `--selftest`. With no flag it prints the block. The block is
a `figure` float holding a `tcolorbox` (main.tex loads `tcolorbox` with `most`), caption below, labelled
`box:worked_example`; the item id and model key are in a LaTeX comment above the float, not in the box.

### Phase 2 command lines

```
python docs/appendix_statistics.py --write            # after D2 places the two marker pairs in taxonomy_content.tex
python docs/appendix_statistics.py --check
python -m full_run_28092026.worked_example --write    # writes overleaf_source_04102026/appendices/worked_example.tex
python -m full_run_28092026.worked_example --selftest
```

D2 inputs the box where it wants it (`\input{sections/appendices/worked_example}` on Overleaf) and may wrap or
re-place it; the file carries its own BEGIN/END markers, so a regeneration can also be pasted between markers.

## Block labels

`tab:area` and `tab:levels_agreement` (both `table*`, in `appendices/taxonomy_content.tex`);
`box:worked_example` (`appendices/worked_example.tex`).

`tab:area` is two side-by-side panels (chemical, electrical, civil on the left; industrial, mechanical and the total
on the right), each `Domain and area | Easy | Int. | Adv. | Total` with shaded branch-subtotal rows and italic
domain rows; the 42 area names are the manifest's identifiers in title case (two spelled out by hand:
"Volumetric Properties of Pure Fluids", "Permeability, Seepage, and Effective Stress"), which agree with the
five-branch overview figure. Totals: 58 / 58 / 34 / 150; 30 per branch. Not compiled here (no TeX install on this
machine): the area column is `p{118pt}` so the long names wrap; the writers should check the two panels on Overleaf.

## Worked example

**Chosen:** `qwen3-235b-a22b-2507` on `gas_viscosity_kinetic_theory#0` (chemical, transport phenomena,
Intermediate, scalar; main store). Found with `--candidates`, which lists every main-store response with a
readable wrong answer, at least one milestone matched deterministically, judge rulings of all three kinds and an
arithmetic flag: 14 responses in the whole store, 8 of them on this one template by this one model, the rest
longer (Advanced items of 5 to 13 milestones, 1,000 to 3,700 tokens, or 950-character questions).

**Why this one.** It is a lower-tier model (FAC 0.886, the five-model lower group), the question is 426 characters
and the response 745 tokens, 7 milestones, and every ruling in the box holds up on reading:
- matching settles `MT` under c = 10^3: the response works in kg/mol where the correlation takes g/mol, so it
  states M*T a thousand times smaller than the gold;
- the judge rules `T_star` and `epsilon` not needed (the question gives the collision integral), `sigma2` and
  `denominator` reached (stated in m^2, a factor 10^-20 the fixed unit factors do not hold, so matching could not
  find them), and `numerator` and `viscosity` missing (both wrong, by the same unit slip);
- the arithmetic check flags the final division: (2.1155e-5)/(3.128e-19) printed as 6.762e-5 Pa·s, operands give
  6.763e13, place value 1e-8;
- the stated answer 6.762e-5 against the target 2.139e-5 is 216% off at c = 1 and matches under no factor, so
  v(r) = 0; m(r) = (1 matched + 2 ruled reached)/7 = 3/7 = 0.43 = the stored `e5_strict` 0.4286.
The judged step check is not in the box (the brief does not list it); on this response the judge labels all
eight sent steps "Alternative Correct", i.e. it misses the unit slip.

**Passed over.** `#7` of the same template and model (the shortest candidate, two deterministic matches): its
second match is an artefact of the milestone reader (next section), not a reached milestone. `#13`: the judge
rules a wrong `numerator` reached. Gemma on `best_hydraulic_rectangular_section#6` has no unit-factor match.

## Open items

1. **Milestone reader and unicode superscript exponents (for the orchestrator / WS-G; evaluator code, not mine).**
   `milestones.numbers` folds `×10^{-19}`, `x 10^-19` and `\times 10^{-19}` into one number but not the unicode
   form `1.152 × 10⁻¹⁹`, which it reads as 1.152 and 10; `answer.values` and `arith.normalise` do fold that form.
   E3 therefore misses a milestone stated this way (it then goes to the judge) and can match the bare mantissa
   under a unit factor (the `#7` artefact: 1.152 matched the numerator 3.196e-4 at the 3,600 factor). Sized with
   `python -m full_run_28092026.worked_example --reader-gap` (free, no store changed): the form occurs in 19 to 24%
   of Claude Sonnet 5, GLM-5.3 and GLM-5.3-Flash responses, 14% of Muse Glimmer, 6 to 8% of Qwen3-235B and
   Kimi K3, and in none to 1.5% of DeepSeek, Gemma, Gemini and GPT-5.4 mini. Folding the exponent would change E3
   on at most 70 responses of one model (GLM-5.3-Flash: 81 milestones newly reached, 5 matches lost) and 17 of
   Claude Sonnet 5 (16 gained, 6 lost); it would change the missed-milestone list, hence the prompt, of up to 65
   judge jobs per model (183 over the roster, new paid calls at the judge's rate). Small in size, but uneven across
   models and in the direction of the "matching alone" comparison between Claude Sonnet 5 and DeepSeek V4.1 Flash
   that C1 and C2 report, so it should be a decision, not left silent: either fix `numbers()` and re-score E3 and
   the judge's new jobs in WS-G's one re-score, or state the reading in the scoring appendix as it is. The full
   per-model table is in the script's output.
2. `tab:levels_agreement` fills itself once WS-E's `score_levels.py` writes `agreement.json` (KAPPA IN); nothing to
   do here but re-run `--write`. If D7 adopts the majority labels and the manifest is regenerated, `tab:area` reads
   the new manifest on the next run.
3. Neither table nor the box was compiled (no TeX on this machine). The agreement table has 12 columns in
   `\footnotesize` at `table*` width; if it overflows on Overleaf, drop the quadratic column or the interval on the
   linear form (both are one line each in `table_levels_agreement`).
