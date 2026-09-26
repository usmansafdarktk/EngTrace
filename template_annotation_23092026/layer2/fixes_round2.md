# Layer 2, round 2: the five templates the experts objected to, verified and fixed

Written 2026-09-26. `RESULTS_round2.md` (D-107) lists the five templates that drew a rejection in
round 2: three rejected by majority (`basic_stress_strain`, `utube_manometer`, `multi_segment_rod`)
and two approved 2 of 3 over one rejection each (`newtons_law_shear_stress`, `finite_convolution`).
Every objection was checked against the instance the expert cited before any edit, and all five
are fixed. The fixing agents were also asked to check the physics next to each objection, so the
experts do not find it in round 3: the manometer's liquid pairings, and two table rows of the same
kind as round 1's molten zinc. Every edited template carries a "Layer 2 fix (2026-09-26)" docstring
paragraph. Per-template rates below are the fixing agents' measurements (500 seeds, 20,000 where
marked); the corpus-wide numbers are the committed scripts'.

## Corpus-wide result

| Check | Command | Result |
|---|---|---|
| Gate | `python -m template_annotation_23092026.layer0.gate --seeds 500` | **150 of 150 pass**; register absorbed 1,494 lines; advisory T5 60, T7 78 (was 81: all three mechanical templates now pass T7) |
| Tie census | `python -m template_annotation_23092026.layer0.tie_census --seeds 500` | 12 templates with a tie, was 13: `basic_stress_strain`'s e-format artefact (D-017) is gone; none of the five has a tie |
| Item pool | `tests/constants_integrity/c3_instance_dump.py --diff`, 150 × 300 seeds, HEAD `85ccf8d` dumped from a detached worktree against the fixed tree | **exactly these 5 templates moved**: question 1,053 / answer 1,327 / solution 1,355 of 45,000; errors 0 → 0 |
| Cumulative since `b1445ad` | the same, against the pre-Layer-0 dump | 77 templates moved: question 5,144 / answer 10,598 / solution 16,669 (`../layer0/item_pool_impact.md`) |
| Review-app rendering | `python -m template_annotation_23092026.layer2.markdown_scan` (200 seeds) | 31 questions and 63 solutions lose characters when the app renders them (65 templates in either); see below |

Only the five assigned `template_*` functions differ from HEAD, checked function by function with
`ast`; the other changes are three helpers (`_sig4_or_exact`, `_floor_sig` in the chemical file,
`_fluids_mix` in the mechanical one), the manometer's pairing tables with their asserts, and one
import. The reports were generated on the working tree before the commit that contains them, so
they name the previous HEAD.

## Per template

"Q" is the question text, "gold" the answer block; both are shares of 500 seeds changed against HEAD.

| Template | Round 2 | Claim, verified | Fix | Moved |
|---|---|---|---|---|
| basic_stress_strain | rejected 3 of 3 | load, diameter and elongation drawn independently with no material: seed 2202 implies E = σ/ε = 1423 GPa, seed 2204 565 GPa; at HEAD 33 of 500 instances implied E above tungsten's 411 GPa, the worst 5946 GPa | the material is drawn once per instance through round 1's exclusion (`_tension_member_material`) and named in the question; a load above σ_allow·A is redrawn alone in [cap/4, cap]; the elongation is δ = σL/E with the material's E, printed at max(3 dp, 4 s.f.) and bound; stress printed the same way, since a polymer bar can sit under 1 MPa; area at 4 dp, L in whole mm (HEAD printed `2029.9999999999998 mm`); asserts that σ is within the allowable and σ/ε within 0.2% of E | Q 100%, gold 100%; σ/ε reproduces E within 0.077% (0.097% at 20,000); material shares 4.29–4.75%; load redrawn by the cap on 57.4%; 0.8% of attempts redrawn whole |
| utube_manometer | rejected 3 of 3 | mercury manometer on liquid oxygen (seed 2201); `PIPE_FLUIDS` also holds liquid nitrogen, and only this template reads the table | cryogenic pipe fluids redraw the attempt; the manometer liquid is redrawn alone until it is denser than and immiscible with the pipe fluid, after Çengel & Cimbala, *Fluid Mechanics* 4th ed. §3-2: "The two fluids must be immiscible, and ρ2 must be greater than ρ1". Pipe fluids are classed gas / hydrocarbon / water-miscible organic / aqueous and manometer liquids mercury / water-based / organic; water-based liquids never touch aqueous or water-miscible pipes, organic liquids never touch hydrocarbon or water-miscible ones. Bromine (fuming, toxic, vapour pressure about 23 kPa at 20 °C) and sodium polysulfide (no concentration stated; the pure salt is a solid) are never drawn. No row deleted (D-031). kPa answers under 1 kPa now keep 4 s.f. (HEAD printed 17.658 Pa as 0.018 kPa) | Q 54.4%, gold 59.4%. At HEAD 55.2% of instances used a pair now excluded. Manometer liquids at 20,000 seeds: mercury 8.7% → 30.9%, other water-based 8.7–10.5% each, halogenated 3.9–4.7% each, the two oils 1.7% each (gas pipes only); pipe fluids move by at most 0.95 points. 185 of 421 density-admissible pairs remain |
| multi_segment_rod | rejected 2 of 3 (both approved in round 1) | US deformations at a fixed 4 dp: seed 2203 sums 0.0036 + 0.0027 − 0.0054 = 0.0009 in against an exact 0.000800 in, 12.5% off; at HEAD 197 of 230 US instances had a segment under 3 s.f., the worst total 18.9% off. Caused by round 1's load scaling, which is kept | segment deformations at max(4 dp, 4 s.f.), bound; the total is the exact sum of the displayed values; a draw whose displayed total is more than 0.5% from the exact total is redrawn (0.12% of attempts, all near-cancelling rods); a zero internal force reads "(No force)", not "(Compression)"; zero node loads (8 of 500 at HEAD) are redrawn; "1 kip" singular; areas keep trailing zeros | Q 22.6% (104 seeds by the "1 kip" wording alone), gold text 82.8%, gold value 74.8%; worst total error 0.16% |
| newtons_law_shear_stress | approved 2 of 3 | no laminar check: Re = ρVY/μ is 31,046 for water at seed 2202 and 87,828 for gasoline at 2203; at HEAD 65.2% of seeds exceeded Re 1000, from 21 of 30 liquids | Re capped at 1000, about 30% under the plane-Couette transition near Re 1440 on the full gap and wall speed (Tillmark & Alfredsson, J. Fluid Mech. 235 (1992): 360 ± 10 on half-gap and half velocity difference). A draw above v_hi = min(2.00 m/s, v at Re 1000) is redrawn in [v_hi/4, v_hi] at 3 s.f., floor 0.01 m/s; the gap is redrawn only below that; the fluid never is. τ and F print at 4 s.f., lengthened to the exact value only where 4 figures would be a half-way tie (D-037), F from the displayed τ; one line per step, so T1 and the census now read 1,500 lines where they read none; V, Y and A keep trailing zeros | Q 74.2% (values on the 326 seeds that were above Re 1000; 45 more gained trailing zeros), gold 100% (format); gap redrawn on 4.4%, speed at the 0.10 m/s floor 0.4%, per-fluid shares 3.04–3.59%; highest Re 995.8 |
| finite_convolution | approved 2 of 3 (same expert rejected in round 1) | the origin is marked only by asterisks, which the app renders as italics, and "x[0] = v" cannot locate it when v repeats: at HEAD the origin value repeats in its own sequence on 59.4% of seeds, e.g. seed 2201, x = {2, *-2*, -2} | every sequence in the question, the Given block and the answer is followed by its index list ("x[n] = {2, -2, -2} for n = -1, 0, 1"); the marker stays but carries nothing alone; Step 3 names the start index | Q 100% (value lists identical on all 500), gold text 100% (the appended index list only; the comparator's index-to-value map is identical and matches both ways on all 500) |

## What the review app showed the experts

The app renders questions and solutions with `st.markdown`; the models receive the raw text.
`markdown_scan.py` parses every question and solution as the app does (CommonMark plus Streamlit's
inline math) and records the constructs that remove characters from what the reader sees. At 200
seeds per template it finds them in **31 questions and 63 solutions, 65 templates in either**
(`markdown_scan.md` lists each with an example). Almost all are one of two things; the few
intended ones are two origin markers, two deliberate italics and two LaTeX reactions:

- **Unspaced products.** `cos(2*pi*2077*t)` or `59.6*R_B ... 2*x` pair their asterisks into
  italics, and the multiplication signs vanish from the rendered text.
- **Dollar amounts.** A pair of dollar signs becomes inline math, so "$148 ... $6.18" loses both
  signs; six industrial templates print currency this way on every seed.

Three question hits are intended: `finite_convolution`'s origin marker, now redundant, and the
LaTeX reactions of `batch_moles_vs_conversion` and `gas_phase_concentration`. One is a defect of
the kind round 2 found: `signal_operations` marks its origin only by asterisks, in the question
and the answer, with no index list. Separately, Streamlit renders a single newline as a space, so
every solution's one-equation-per-line layout ran together in the app.

None of the five round-3 templates has a lossy construct except the convolution marker, so round
3 is not affected. Rounds 1 and 2 were: the experts judged these 65 templates from a view with
multiplication signs or dollar signs missing. No verdict is known to be wrong because of it, but
none was made on the text the models read.

## Round 3

The five templates go back to the experts of their branches, without plants, with a fresh
hand-check instance (seed 2301):

```bash
python -m template_annotation_23092026.layer2.build_tasks --only template_basic_stress_strain,template_utube_manometer,template_multi_segment_rod,template_newtons_law_shear_stress,template_finite_convolution --round 3
python -m template_annotation_23092026.layer2.make_kits --round 3
python -m template_annotation_23092026.layer2.score --round 3 --labels <round-3 folder> --prev-labels <round-2 folder>
```

That is three mechanical experts with three items each, and three chemical and three electrical
experts with one item each.

## Open for the owner

- **The app.** Show questions and solutions as plain text, as the models read them, before any
  further round; and decide whether the templates approved in rounds 1 and 2 through the rendered
  view need a look in plain text.
- **`signal_operations`.** The same origin defect as `finite_convolution`; not in round 2's scope.
- **Left as they are, with the reason in the docstrings:** milk pipes still pair with mercury and
  the halogenated liquids (immiscible, if not food-grade); `basic_stress_strain` has no lower bound
  on stress (7% of instances under a tenth of the allowable); mercury is now 30.9% of manometer
  liquids, which the pairing rule forces.
- **A comparator note.** The strict numeric comparator in `tests/comparators` (half a unit in the
  gold's last digit) rejects an answer computed without intermediate rounding whenever the gold
  carries the rounded chain D-016 prescribes; on the Newton template that is 59% of exact answers,
  all within 0.132%. The full run grades with `evaluator_pilot_17092026/evaluators/answer.py`,
  whose 0.2% relative band accepts them (D-105), so this matters only if that comparator is ever
  used for grading.
- **D-044.** Answer displays changed in all five templates. The owner asked for the five to be
  fixed (2026-09-26), which this record takes as the sign-off.
