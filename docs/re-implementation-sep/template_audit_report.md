# Template Structure Audit — feasibility of deterministic milestone verification

**Date:** 2026-09-05 · **Scope:** all 150 templates under `data/templates/branches/` · **Status:** read-only feasibility study; no template was modified.

---

## 1. Verdict

**Proceed, with a modified design.** The hypothesis survives, but not in the form it was posed, and the reason it survives is not the reason the prompt assumed.

Three findings drive the verdict.

**The crux question was the wrong worry.** The premise was that intermediate quantities might be trapped inside f-string expressions, so instrumentation would mean rewriting the computation rather than adding emission. Measured across all 150 templates: **93.0% of the 6,166 f-string interpolations are already bound variables** (a bare name, or `round(name, k)`, or a conditional over two bound names). Running every template and capturing its locals at return, **94.5% of step-result numbers and 96.7% of final-answer numbers are recoverable from values the generator already holds**; 110 of 150 templates recover 100% of their step values and 132 of 150 recover 100% of their answer values. Only 9 templates of 150 are genuinely `values_bound = no`. Instrumentation is overwhelmingly *emission*, not rewriting.

**The real obstacle is that a per-template static schema is impossible.** 72 of 150 templates (48%) change the governing equation, the step count, or the set of computed quantities based on a sampled parameter. The milestone set is a property of the **instance**, exactly as suspected — and this is not a fixable accident. Several civil and industrial docstrings name parameter-dependent reasoning paths as a deliberate design goal. Any design that declares milestones statically per template is dead on arrival; the trace must be emitted per instance.

**The extraction side is the weak link, and it is worse than the templates.** This was flagged as the most important question in the audit, and it is. `evaluation/engineering_parser.py` reduces an entire solution to **one unlabelled, unitless float**. Against the repository's own gold traces it fails in ways that are not subtle: it reads the gold answer `**4,921**` as **`4.0`**, it cannot see the `**Final Answer**` marker that 5 templates use, and **74.7% of gold answer blocks contain more than one distinct number** — a quantity a single float structurally cannot represent. These are not hypothetical: in the archived model outputs, 89 rows have a gold answer the parser mis-reads, and in **43 of them (48%) the model states the correct value and is scored 0 anyway**.

That last finding reframes the project. Instrumenting the gold trace is roughly a third of the work and the easy third. **If extraction is not fixed, deterministic verification fails for reasons that have nothing to do with the templates** — and, separately, some already-published numbers are affected (§8.4).

### Cost, at a glance

| Class | Definition | n | Effort |
|---|---|---:|---|
| **A** — trivial | Fixed step chain, scalar answer, reported values bound. Mechanical instrumentation. | 49 | low–medium |
| **B** — branching | Governing equation / step count / milestone set changes with the sample. Needs a per-instance schema. | 34 | medium–high |
| **C** — non-scalar | Answer or key intermediates are vector, array, symbolic, classification, or multipart. Needs typed milestones. | 51 | medium–high |
| **D** — resists | No numeric content, or the trace's own chain does not reproduce its own answer, or a search/iteration log with no stable symbols. | 16 | high, and 8 need the template rewritten, not instrumented |

Classes are assigned by dominant obstacle with precedence D > C > B > A; a template that both branches and returns a vector is C, with the branching recorded in `notes`. Per-template assignments are in `template_inventory.csv`.

**Rough order of magnitude:** ~49 templates near-free, ~85 needing real design work, ~16 needing a decision about whether they belong in a deterministically-verified benchmark at all. Plus a parser rewrite that is a larger and more urgent job than any of it.

---

## 2. Structural uniformity

All 150 templates share one macro shape: sample parameters (often inside a rejection loop) → compute into local variables → concatenate one large f-string with `**Step N:**` markers → `**Answer:**`. All 150 import and execute cleanly; across 15–25 seeds each I observed **zero exceptions and zero failed asserts**, and the branch audits independently ran 180,000 (civil) and 6,000 (industrial) instances with the same result. This is a healthy corpus.

Uniformity ends at the output contract, and the exceptions matter because they are exactly what a mechanical migration would trip over.

**Step markers are not reliable identifiers.**

- `template_levenspiel_plot_interpretation` **restarts its numbering** — the emitted sequence is `1,2,3,1,2,3,4,5,6`. `**Step N:**` is not unique within a trace.
- Three electrical templates emit malformed markers that a strict `\*\*Step (\d+):\*\*` regex silently drops:
  ```python
  f"**Step 2: ** Calculate Euclidean Distance (d)\n"          # euclidean_distance_binary, BFSK branch (~50% of instances)
  f"**Step 3: ** Discrete-to-Continuous (D/C) Conversion*\n"  # cd_dc_system_analysis (100%)
  f"**Step 3: Calculate All Output Values**\n"                # finite_convolution (100%)
  ```
  For `cd_dc_system_analysis` the dropped step is the one that computes the answer.
- 5 chemical templates terminate with `**Final Answer**` or `**Final Answers:**` rather than `**Answer:**`. In the archive that is **245 of 2,200 rows** whose answer marker the evaluator's parser cannot see.

**Step counts are stable but milestone density is not.** Step counts run 2–9 (median 4) and vary by instance in only 9 templates. But milestones per step vary by an order of magnitude: `template_rational_method_peak_flow` reports one computed number across three steps; `template_normal_depth_iteration` reports a median of 21 across four; `template_work_isothermal_virial` packs 14 into 5.

### Are civil and industrial structurally different?

Yes, and the difference is real but **narrower and differently located than the spec-driven provenance scaffolding implies**.

They are genuinely different in ways that help:

- **Asserts.** 60 templates carry `assert` statements — all 60 are civil (30) or industrial (30). The other three branches have **zero**. Several already encode exactly the invariants a structured verifier would want: `assert abs((Ay+By)-(P+W)) < 0.02` (equilibrium), `assert O1 < O2 < O3 < I3` (monotonicity).
- **Round-then-recompute.** Civil holds printed-value == stored-value in **29 of 30** templates; industrial in **26 of 30**. This is the single most valuable property in the corpus, because it removes the tolerance-definition problem: a verifier can demand exact agreement at the printed precision rather than guessing an epsilon.
- **No table interpolation anywhere.** Every `constants.py` access in both branches is an exact dict-key hit. Table lookups turn out to be the *easiest* part of the corpus, not the hardest.
- **Determinism.** No `numpy`, no `scipy`, no `try/except` in civil; stdlib only in both.

They are *not* better in the way the prompt expected:

- **Variable binding is not materially better.** Inline arithmetic inside f-string placeholders, counted per branch: mechanical 90, **industrial 92**, electrical 63, chemical 54, civil 27. On raw count industrial is the *worst* of the five. 46 of industrial's 92 are a lossless scaled-integer idiom (`{pb4 / 10000:.4f}`) unique to that branch, leaving ~46 genuine inline computations — still more than civil's 27.
- **Branching is more common, not less.** Civil branches in 13 of 30 (43%); several docstrings state parameter-dependent reasoning paths as a design goal. The spec appears to have *rewarded* the property that makes static schemas impossible.
- **Units are no more machine-readable.** Unit strings are hard-coded literals in the f-strings in both branches. Civil has exactly three unit-holding variables in 30 templates, all in one function. `constants.py` has **no unit field anywhere** in either branch — units live in identifier suffixes (`_KN_M3`, `_CM_S`), and provenance tags (`[ON-DISK]`, `[VERIFY: Das]`, `[POLICY: sampling-only]`) are `#`-comments, invisible to `import`.

**Judgment:** civil and industrial are perhaps 20–30% cheaper to instrument, and the saving comes from *determinism and rounding discipline*, not from code hygiene or provenance metadata. That is worth knowing, because it means the property to preserve in any future authoring spec is "round then recompute, and assert your invariants" — not "tag your constants."

---

## 3. Are intermediate quantities bound to variables?

**Yes, far more than assumed.** Two independent measurements.

**Static (AST sweep, all 6,166 interpolations):**

| Branch | interpolations | already bound | computed inline | % inline |
|---|---:|---:|---:|---:|
| chemical | 886 | 821 | 65 | 7.3% |
| civil | 1,132 | 1,095 | 37 | 3.3% |
| electrical | 1,198 | 1,082 | 116 | 9.7% |
| industrial | 1,683 | 1,558 | 125 | 7.4% |
| mechanical | 1,267 | 1,178 | 89 | 7.0% |
| **total** | **6,166** | **5,734** | **432** | **7.0%** |

70 of 150 templates have **zero** inline-computed interpolations.

**Dynamic (run each template, capture locals at return, match printed numbers against bound values):**

| | tokens | exact | rounded | scaled | missing | recovered |
|---|---:|---:|---:|---:|---:|---:|
| step-result numbers | 17,177 | 65.9% | 24.4% | 4.2% | **5.5%** | **94.5%** |
| answer-line numbers | 4,133 | 49.9% | 41.2% | 5.7% | **3.3%** | **96.7%** |

"Scaled" means the value is bound but at a different scale — an integer-scaled `Fraction`, or a unit conversion. This is a real and important pattern, not noise: the industrial acceptance-sampling templates hold `pa4 = 6011` (an exact integer) and print `{pa4/10000:.4f}`. The value **is** bound, losslessly, but an extractor must know the scale factor.

Where inline computation does occur, it is mostly *not* the step result. The dominant idioms are operand restatements (`{2*L}`, `{L/2}`, `{Hc/2}`), display unit conversions (`{load*1000}`, `{X*100}` for a percentage), and substitution displays (`{round(1/k, 2)}`) shown alongside a bound result.

But there are genuine exceptions where instrumentation means **rewriting the computation**, and they cluster hard:

`template_batch_reactor_second_order` — five distinct computed quantities in Step 6, exactly one bound:

```python
f"t = (1/{k}) × (1/{C_A} - 1/{C_A0})\n"
f"t = {round(1/k, 3)} × ({round(1/C_A, 3)} - {round(1/C_A0, 3)})\n"
f"t = {round(1/k, 3)} × {round(1/C_A - 1/C_A0, 3)}\n"
f"t = {round(time, 2)} s\n\n"
```

`template_axial_deformation` — three unbound quantities in a single line, all unit conversions, which is precisely the arithmetic a dimensional checker most needs declared:

```python
f"delta = (({load * 1000}) * ({length * 1000})) / (({round(area, precision+1)}) * ({material_E * 1000}))\n"
```

`template_mean_variance` — the worst case in the corpus. Every per-term intermediate is generated inside a comprehension inside an f-string, producing 9–15 printed numbers with no Python variable at all:

```python
mean_interm_vals = " + ".join([f"{v * p:.{precision}f}" for v, p in zip(values, probabilities)])
```

**Bottom line:** 95 templates are `yes`, 46 `partial`, 9 `no`. For the 9, and for the ~10 `partial` templates whose inline expressions are step results rather than operand restatements, instrumentation requires touching the computation. That is ~19 of 150 — a real cost, but an order of magnitude smaller than the prompt feared.

---

## 4. Are units recoverable?

**Not automatically, and this is a larger gap than the value problem.**

Units exist **only as literal strings inside the f-strings**. Across the whole corpus, unit-holding variables are rare enough to enumerate: civil has three (all in `template_beam_deflection_formula`); industrial has one family (`cfg['unit']` in the SPC/capability templates, 52 interpolations); mechanical's dual-unit templates use bound `stress_unit` / `unit_P` variables in some places and inline conditionals in others, sometimes in the same function:

```python
f"A = pi * (d/2)^2 = ... = {round(area, precision+1)} {'mm^2' if use_si_units else 'in^2'}\n\n"
```

`constants.py` carries **no unit field in any branch**. Units are encoded in identifier suffixes (`_KN_M3`, `_CM_S`, `LINE_TASK_TIME_S`) — a parseable convention at the top level that stops at the first nesting: inside `AISC_W_SHAPES` the per-field units exist only in a comment header.

**Could a unit be attached to every emitted quantity without hand-annotating 150 templates?** Partially. The unit literal almost always immediately trails its interpolation, so a regex over the format string recovers it mechanically for the large SI-only majority. Four failure modes defeat that:

1. **Variable units.** `f"...{delta_out} {unit_out}"` — the unit is not in the format string.
2. **Instance-dependent units.** Three electrical templates use a hand-rolled SI-prefix formatter, so the same physical quantity reports in nJ / µJ / mJ or kHz / MHz / GHz depending on the sampled magnitude. The `unit` field is itself a per-instance value.
3. **Symbolic units.** `template_force_method_continuous_beam` prints two milestones as `{d_B0:.1f}/EI  (kN*m^3 over EI)`. A dimensional checker cannot validate either in isolation; only their *ratio* is meaningful, and only because EI cancels downstream. This is the clearest single counterexample to "dimensional checking is deterministic."
4. **Scale factors folded into the formula.** `template_aoq_ati_rectifying` prints AOQ in `ppm` from a formula that already multiplies by `10^6`, alongside a dimensionless `Pa` and a count-valued `ATI` — three unit systems in one trace. A naive unit check would *reject* a correct dimensionless answer.

**On dual units — the prompt's premise is wrong.** The prompt says "many templates branch on `use_si_units = random.choice([True, False])`". In fact **exactly one file does**: `mechanical_engineering/mechanics_of_materials/stress_and_strain.py`, 5 templates. Two more templates are dual by other means (`template_beam_deflection_formula` via `use_si`, `template_kinematic_viscosity` via SI/CGS Stokes) and `template_scs_curve_number_runoff` is dual by necessity (mm question, inch-bound TR-55 equation). That is **8 of 150**, not "many." Dual units are a contained problem, not a systemic one.

**How much would dimensional checking actually buy?** Less than hoped, and it is concentrated where it is least needed. 24 of 150 answers are dimensionless (probabilities, capability indices, z-scores, Reynolds numbers, control-chart limits); several more are counts or percentages. In the industrial branch specifically, only four templates carry a real physical unit. Dimensional checking is a genuine safety net for the mechanical and civil branches and close to vacuous for industrial.

---

## 5. Instance-dependent structure

**Widespread: 72 of 150 templates (48%).** Per branch: electrical is heaviest (22 of 30 contain a regime-selecting `random.choice`), then mechanical (16), civil (13), chemical (12), industrial (9).

`template_reynolds_number_flow_regime` is the canonical case and it is *worse* than the prompt describes — it lives in **chemical**, not mechanical:

```python
geometry = random.choice(['pipe', 'flat plate'])
...
if geometry == 'pipe':
    if reynolds_number < 2300:      regime = "laminar"
    elif reynolds_number < 4000:    regime = "transitional"
    else:                           regime = "turbulent"
else:  # flat plate
    if reynolds_number < 500000:    regime = "laminar"
    else:                           regime = "turbulent"
```

Three things change with one `random.choice`: the *meaning* of `characteristic_length` (diameter vs plate length), the *threshold set* (`{2300, 4000}` vs `{500000}`), and the *arity of the answer space* (three classes vs two). A milestone schema for "Re" would need to carry a branch-specific comparator table, not a numeric target.

The same shape recurs across the corpus in more damaging forms:

- `template_axial_deformation` — `solve_for_load` **inverts the governing equation** (`delta = PL/AE` becomes `P = delta·AE/L`), flipping which quantities are given, which is unknown, and the *dimension of the answer*. Combined with `use_si_units` that is four structurally different traces from one function.
- `template_limiting_reactant` — the same symbol `N_A` has the formula `N_A0(1-X)` in half the instances and `N_A0 - (a/b)N_B0X` in the other half. Same symbol, same unit, same step position, different equation. Any per-template formula registry must be keyed on a latent branch variable **the emitted trace never states**.
- `template_chart_pair_selection` — the *governing constant family* changes: small `n` uses A2/D3/D4/d2 on R-bar, large `n` uses A3/B3/B4/c4 on s-bar.
- `template_superposition_electric_field` — `num_charges ∈ {2,3}` changes both the step count and the quantity count, and step numbers are generated in a loop.
- `template_decimation_aliasing_analysis` — the `will_alias=False` branch **omits three computed quantities entirely**. The milestone *set*, not just the values, changes.

**Implication, stated plainly:** a design that declares milestones statically per template cannot work. There is no repair that keeps the static schema — the branching is intrinsic and in the newer branches deliberate. The trace must be emitted **per instance**, by the generator, at generation time. That is not a limitation of this corpus so much as a design constraint the next stage must accept as given.

---

## 6. Answer types that resist numeric milestone checking

This is where the prompt expected the idea to struggle, and it does — though the specific catalogue is somewhat different from the one supplied.

| answer_type | n | verdict |
|---|---:|---|
| scalar | 89 | handled |
| multipart | 32 | handled with a milestone **list**; partial credit needs design |
| symbolic | 9 | **not handled** by numeric checking |
| vector | 8 | handled with component milestones + sign semantics |
| array / sequence | 7 | needs a different representation |
| classification | 5 | **not handled**; needs a categorical milestone type |

The counts understate the problem, because a template's *answer* can be scalar while its *intermediates* are not. Counting templates with any non-scalar content the schema must represent, the figure is **61 of 150**.

### 6.1 Symbolic — 9 templates, genuinely unhandleable numerically

The prompt named the BER case. It is real and it is worse than described, because the *precision is itself random*:

```python
precision = random.randint(1, 2)
ber_approx_str = f"{ber_coeff:.{precision+1}f} * Q({q_arg:.{precision+1}f})"
```
```
BER approx 0.50 * Q(8.07)        # one seed
BER approx 0.750 * Q(2.888)      # another
```

`{id, symbol, value, unit}` has no slot for this. The value is a *pair* (coefficient, Q-argument) plus an uninterpreted function symbol. Either `Q()` is modelled as an opaque operator, or the template is changed to evaluate `Q(x) = 0.5·erfc(x/√2)` — which changes the pedagogy.

Two cases are harder still. `template_incompressible_continuity` **contains no numbers at all** — every value is a polynomial string built by local helpers, and checking it needs a CAS equivalence test including a discarded arbitrary function `f(x)`. `template_undamped_response_initial_conditions` assembles its answer by **string surgery after the last numbered step**, so Step 4 prints the full two-term equation and the Answer prints a different, term-dropped one, in the same trace.

### 6.2 Classification — 5 as the answer, ~12 more as a gating step

The dangerous cases are not the answers but the **steps**. `template_terzaghi_strip_footing_bearing` Step 1 emits **no number** — it is the verdict "general shear failure governs," derived from `Dr = 78%` against a stated 70% threshold. Every downstream number is conditioned on it. A model that computes all four numbers correctly under the *wrong* failure mode produces a trace that is internally consistent, dimensionally valid, and wholly wrong. Only treating step 1 as a non-numeric milestone catches it.

Worse, near a threshold the numeric check and the classification disagree about what "correct" means: a model 1% off on Re passes any sane numeric tolerance and flips the class.

### 6.3 Array / iteration — 7 templates, need a different representation

`template_normal_depth_iteration` is the sharpest case, and it defeats the *shape* of the proposed schema rather than any one field:

```python
for _ in range(5):
    y_next = round(y_curr - g_curr * (y_curr - y_prev) / (g_curr - g_prev), 4)
    eval_lines.append(f"Update {updates}: y_next = ... = {y_next:.4f} m")
    if abs(y_next - y_curr) < 0.002:
        y_curr = y_next
        break
    y_prev, g_prev = y_curr, g_curr
    y_curr = y_next
```

`n_steps` is a constant 4, but Step 3 holds **8 to 28 quantities**, and `y_curr` / `g_curr` are **rebound every pass**. There is no stable `{id, symbol}` for "the third trial depth." This needs a list-of-iterations, not a list-of-steps.

`template_line_balancing_heuristic` is the same problem in a different domain: the trace is a **search log**, whose milestones are task *sets* (`Station 1 contents: {a, b}`), a running counter, and boolean fit tests, with one step per opened workstation.

The prompt's example, `template_levenspiel_plot_interpretation`, is real but its worst defect is elsewhere — see §9.

### 6.4 Multipart — 32 templates, handleable but needs partial-credit design

`template_statically_indeterminate` is the named case and it is instructive because it fails three ways at once. Part (a) is a *pair* (R_A, R_C) and part (b) is a *pair* (σ_AB, σ_BC) — four numbers under two labels. Beyond that:

```python
if use_si_units:
    stress_AB_MPa = stress_AB * 1000   # kN/mm^2 -> MPa; undefined in the US branch
```

The symbol `sigma_AB` maps to `stress_AB_MPa` in SI and `stress_AB` in US — **different units for the same milestone**, differing by 1000. The conversion is then written inline inside the printed equation, so the trace line is not a valid equation: `sigma_AB = (166.093 kN) / (9503.32 mm^2) * 1000 = 17.477 MPa`. And the Answer applies `abs()` to σ_BC while Step 4 prints it signed — one milestone, two different signed values in one trace.

### 6.5 Vector — 8 templates, tractable

The electromagnetics cross-product case is real but is the **least** troublesome of the non-scalar classes: in `template_lorentz_force` all eleven component scalars are bound. It needs three sibling milestones with a shared unit and component-wise comparison, plus a decision about sign semantics (a sign error in one component is a wrong answer; a sign flip in all three is a direction error). That is design work, not an obstacle.

---

## 7. Steps vs. milestones

**No, `**Step N:**` does not correspond one-to-one with computed quantities**, in either direction.

**Many steps compute nothing.** In civil, roughly **15% of all step markers emit no number** a numeric checker can grade. Examples: Steps 2–4 of `template_cantilever_double_integration` are symbolic integration; Steps 1–2 of `template_linear_reservoir_routing_step` are algebraic rearrangement; Steps 1–3 of `template_utube_manometer` and 5 of 6 steps in `template_particle_pathline` are narrative. Three electrical templates have **zero** numeric milestones anywhere.

**Many steps compute several.** `template_slope_deflection_end_moment` emits 8 results across 4 steps, four of them in Step 1 alone.

**Some steps are not on the answer path at all.** `template_manning_trapezoidal_velocity` Step 4 computes `Q = V·A` as a **redundant confirmation** whose result is not the answer. `template_best_hydraulic_rectangular_section` Step 4 is a back-check. `template_basic_buoyant_force` Step 3 is labelled `(Optional)`. A verifier that treats every step result as a required milestone will penalise a model that sensibly omits an optional check.

**Some quantities are emitted outside any step.** `template_work_isothermal_virial` prints `W_ideal_J` in a "For Comparison" block; `template_levenspiel_plot_interpretation` prints a computed `efficiency` outside any step marker.

**Implication:** `is_milestone` is insufficient. It needs a companion `kind ∈ {numeric, symbolic, classification, narrative, check}` and a `required` flag distinguishing answer-path milestones from optional confirmations. The classification and check kinds are exactly the ones that fall back to something other than arithmetic — which is where an LLM judge, if one is retained at all, would be doing work a numeric checker cannot.

---

## 8. Precision and rounding

### 8.1 Which value is the milestone?

Mostly the printed one, and that is good news. The civil convention `x = round(expr, n)` **before** interpolation holds in 29 of 30 civil and 26 of 30 industrial templates. Where it holds, printed value **is** stored value, and a verifier can demand exact agreement at the printed precision — no epsilon guessing. This is the single most valuable property in the corpus and it should be a hard requirement in any future authoring spec.

### 8.2 Where it does not hold

**Rounding is inconsistent within single templates.** `template_basic_stress_strain` uses three conventions in one trace (`round(x, precision)`=3, `round(area, precision+1)`=4, `{strain:.3e}`). `template_composite_shafts_series` prints segment angles at 5 dp and their total at 4 dp, so **the printed `phi_AB + phi_BC` does not equal the printed `phi_total`**.

**Format-truncation of full-precision values** occurs in 4 industrial templates and 1 civil, where the trace shows an unrounded value and then rounds it — a *second, independent* rounding path for the same milestone. This needs modelling as `{raw_value, displayed_value}`, not one `value`.

**One case is an outright bug.** `template_cantilever_double_integration` keeps `delta` at full precision, prints it at 5 dp, then computes `delta_mm = round(delta*1000, 1)` from the *unrounded* value:

```python
delta = w * L ** 4 / (8 * E * Ix)          # NOT rounded
f"delta = {delta:.5f} m\n"                 # Step 5
delta_mm = round(delta * 1000, 1)          # from FULL precision
f"delta = {delta:.5f} * 1000 = {delta_mm:.1f} mm\n"   # Step 6
```

**198 of 4000 seeds (4.95%)** emit a Step-6 line that contradicts its own arithmetic (`0.01685 * 1000 = 16.8 mm`; 16.85 rounds to 16.9). The sibling template 300 lines above carries the fix and a comment naming the review cycle that found it.

### 8.3 What tolerance is implied — and why one number will not do

Four distinct regimes, and no single tolerance serves them:

1. **Exact-at-printed-precision** (~55 templates): ±0.5 in the last displayed digit. No floating-point slack needed.
2. **Catastrophic cancellation.** `template_max_hump_height_no_choking` computes `dz_max = E1 - Ec`, a small difference of near-equal energies. A 1% relative tolerance on `dz_max` implies **~0.05% on E1** — tighter than E1's own printed precision. Per-milestone tolerance must be derived from the *chain*, not the symbol.
3. **Gold traces internally inconsistent by design.** `template_best_hydraulic_rectangular_section` asserts its own back-check only to 1.5%; `template_normal_depth_iteration` accepts a 2% residual. "Deterministic" cannot mean "exact" everywhere.
4. **Traces that do not reproduce their own answer.** Six mechanical templates **back-solve their inputs**: sample a target ζ, r, or strain, derive the "given," round it for display, then present the unrounded derived value as gold. Measured worst cases: **7.77%** (`template_rotating_unbalance` amplitude), 0.87% (`template_vibration_transmissibility`), and a **full classification flip on ~35% of `template_damping_classification` instances** (ζ recomputed from the rounded `c` = 1.000000854527, which exact arithmetic calls Overdamped while the gold says Critically Damped).

Regime 4 is the serious one. **Any tolerance tight enough to catch a real reasoning error will reject the gold trace on these templates.** They must be fixed before verification, not verified around.

A further wrinkle: `template_hydraulic_jump_energy_loss` **flips display precision by branch** — `y1` prints at 2 dp or 4 dp depending on which depth is given. A fixed per-symbol tolerance is wrong for it.

### 8.4 The industrial scaled-integer idiom is the model to copy

46 interpolations across the industrial branch hold values as exact integers in units of 10⁻ⁿ and divide at print time (`{pb4 / 10000:.4f}`). Combined with `Fraction` arithmetic and explicit "path-agreement" screens that reject any sample where the displayed chain and a full-precision solve disagree, this branch has **already solved** the problem the whole project is trying to solve — it just solves it inside the sampler instead of exposing it. That machinery is the strongest existing asset for this work.

---

## 9. The extraction side — the decisive finding

Instrumenting the gold trace only helps if the corresponding quantities can be recovered from the *model's* prose, which stays unstructured. **They largely cannot, with the current parser, and the current parser is worse than its regex-based design implies.**

`evaluation/engineering_parser.py` reduces a solution to **one float**, via:

```python
match = re.search(r'\*\*Answer:\*\*\s*.*?([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', text, re.IGNORECASE)
```

with no unit handling and two positional fallbacks. `evaluate_entry` derives **both** the gold reference and the prediction from it, then compares at a 2% relative tolerance.

I tested it against the repository's own archive: `error_analysis_annotation/samples/`, 2,200 real model traces from 11 models over 1,068 distinct gold traces.

> **Sampling caveat, stated because it bounds the claims.** That archive is a *filtered failure set* — every row has `final_answer_acc == 0.0` by construction (`sampling_summary.txt`: confirmed-failure filter). It therefore **cannot** measure the parser's overall accuracy, and I make no such claim. It *can* show whether specific rows were scored 0 for a parser reason rather than a reasoning reason, and the gold-trace tests below are independent of the sampling.

### Defect 1 — thousands separators truncate the gold value

The Priority-1 regex runs on the raw text; the comma-stripping happens only *afterwards*, inside `convert_to_float_eng`. So the match stops at the first comma.

```
gold answer block                                          parser reads
**Answer:** a) The Reynolds number is approximately **4,921**.        4.0
**Answer:** ...total required stiffness is **1,256,557 N/m**.         1.0
**Answer:** a) The Reynolds number is approximately **143,412,385**.  143.0
```

Three templates emit `{:,.Nf}` in a milestone position; two of them in the **answer** line — `template_reynolds_number_flow_regime` (chemical) and `template_vibration_isolator_design` (mechanical). For those two templates **the gold reference value is wrong in 100% of instances**, off by factors of 10³ to 10⁶.

In the archive, 89 rows come from those two templates. In **43 of them (48%) the model states the correct value and is scored 0**:

```
claude-opus-4-7   reynolds_number_flow_regime   TRUE gold 4921     parser read 4      model said 4921.5    recorded acc=0.0
claude-opus-4-7   reynolds_number_flow_regime   TRUE gold 141553   parser read 141    model said 141553    recorded acc=0.0
claude-opus-4-7   vibration_isolator_design     TRUE gold 1256557  parser read 1      model said 1256470   recorded acc=0.0
```

**These are false negatives in published numbers.** Because the archive is failure-filtered I cannot state the corpus-wide magnitude — but the affected population is every evaluation of those two templates across all models and all runs, and the direction is one-way (it can only suppress scores). This warrants a re-run before any further submission, independent of the milestone project.

### Defect 2 — the answer marker is not universal

`**Answer:**` is hard-coded. 245 of 2,200 archived rows (5 chemical templates) use `**Final Answer**` or `**Final Answers:**`. For those the parser silently falls through to "last number after an `=`" or "last number anywhere" — with no signal that it did.

### Defect 3 — one float cannot represent the answer

**74.7% of gold answer blocks in the archive contain more than one distinct number** (59.8% of the 1,068 distinct gold traces). For all 32 multipart templates, all 8 vector templates and all 7 array templates, a single float is structurally incapable of representing the answer. It silently grades the *first* number and discards the rest.

### Defect 4 — units are discarded entirely

Demonstrated directly:

```
7.65   <-  **Answer:** The volume is 7.65 L
7.65   <-  **Answer:** The volume is 7.65 mL          # same score
12.5   <-  **Answer:** The reaction at A is 12.5 kN and at B is 30.1 kN
-2.6e-05 <- **Answer:** F = (-2.6e-05 x_hat + -7.0e-06 y_hat + -5.5e-05 z_hat) N
0.5    <-  **Answer:** BER approx 0.50 * Q(8.07)      # grades the coefficient
None   <-  **Answer:** The flow is turbulent
5.0    <-  **Answer:** z[n] = {5, -1, *2*, 5, 6, -2, -9}
724.0  <-  The answer is 1,541,724
```

A unit error — arguably the most characteristic engineering reasoning error — is **invisible** to the current evaluator. This is the sharpest argument *for* the project: dimensional checking would add something real. It is also the sharpest argument that the parser, not the templates, is the binding constraint.

### Defect 5 — step alignment on model prose

- Gold traces average 4.07 steps; model traces average 5.62. **83.8% of pairs have different step counts.** Models do not follow the gold's step decomposition, so positional milestone matching cannot work — matching must be by *quantity*, not by step index.
- On model prose, the parser extracts **no value at all from 4.1% of steps**. That is the optimistic figure, since it counts any number as success.
- 191 rows have a classification-word answer and 95 a symbolic answer — **13% of the archive is ungradeable by a numeric extractor in principle**, regardless of implementation quality.

### Assessment

**Extraction is the weak link, and it is the single most important finding in this audit.** The gold side is in far better shape than assumed — 94.5% of step values are already recoverable from bound variables. The model side is not, and no amount of template instrumentation improves it.

The good news is that the failure is *asymmetric in a useful way*. Milestone verification does not require parsing the model's prose into a *structure*; it requires finding, anywhere in the model's output, a number matching each known gold milestone within tolerance and carrying a compatible unit. That is set-membership matching, not structural parsing, and it is robust to the 83.8% step-count divergence. But it needs a value extractor that returns **(value, unit, context) tuples for every quantity in the text** — which is a from-scratch rewrite, not a patch to the existing 118 lines.

---

## 10. Bugs found (recorded, not fixed)

Per the read-only constraint. Ordered by severity.

**Affecting published results:**

1. **Thousands-separator parse bug** (`engineering_parser.py`) — corrupts the gold reference for `template_reynolds_number_flow_regime` and `template_vibration_isolator_design` in 100% of instances; ≥43 confirmed false negatives in the archive. §8.4 above / §9 Defect 1.
2. **`CP_PARAMS` coefficients 10× too large** (`chemical_engineering/constants.py`) — B given as `E-2` and C as `E-5` where Smith–Van Ness gives `E-3`/`E-6`. Cp(N₂, 300 K) computes to 42.4 vs 29.1 J/mol·K. Affects `template_sensible_heat_temp_dependent_cp` and `template_adiabatic_flame_temperature`, whose T_ad spans 1130–1493 K where methane/air should be ~2200 K. **A numeric checker would faithfully certify this incorrect physics** — which is why milestone verification is necessary but not sufficient.
3. **Unseeded `np.random`** in `template_levenspiel_plot_interpretation` — `random.seed()` does not control it, so the seed recorded in the generated JSONL **cannot regenerate the item**. Verified: three calls under identical `random.seed(12345)` produced three different data tables.

**Trace self-inconsistency:**

4. `template_cantilever_double_integration` — 4.95% of instances print a Step-6 line contradicting its own arithmetic (§8.2).
5. `template_rotating_unbalance` / `template_vibration_transmissibility` / `template_damping_classification` / `template_system_properties` / `template_poissons_ratio` / `template_logarithmic_decrement` — back-solved inputs; the stated chain does not reproduce the stated answer (up to 7.77% error; a classification flip on ~35% of instances) (§8.3).
6. `template_impulse_response_from_lccde` — Step 5 prints `C2 = -6.02`, the Answer prints `-6.024922359499621`, same trace (`fmt()` formats the raw float and ignores its pre-rounded string argument). ~all second-order instances.
7. `template_composite_shafts_series` — printed `phi_AB + phi_BC` ≠ printed `phi_total` (5 dp vs 4 dp).
8. `template_annulus_flowrate` — printed final multiplication does not close at printed precision.
9. `template_vibration_isolator_design` — Step 4 **hardcodes** its own verification: `f"Check: Is r > sqrt(2)? Yes, {r} > 1.414. The condition ... is met."` — printed unconditionally, never tested.

**Malformed output:**

10. Three electrical templates emit malformed `**Step N:**` markers (§2).
11. Five chemical templates use `**Final Answer**` instead of `**Answer:**` (§2).
12. `template_time_to_phasor` / `template_phasor_addition` — negative imaginary parts print as `+ j-51.22` (21/50 and 36/50 seeds).
13. `template_decimation_aliasing_analysis` — unreachable `elif nf == 0` yields the token `omega_a = 0*pi` (61/3000 seeds).
14. `template_signal_operations` — table header reads `New Value y[k')` (mismatched bracket, wrong variable).

**Physical / data plausibility:**

15. `template_gas_phase_concentration` — ε sampled on [−0.5, 0.5] but δ ≤ 0 for these stoichiometries, so ~half of instances state a thermodynamically unreachable expansion factor.
16. `template_two_phase_specific_volume` — names a real substance from `THERMO_SUBSTANCES` but invents V_l/V_v randomly, ignoring `REAL_FLUID_DATA` (Ammonia v_g 0.754 vs true 0.1284 m³/kg).
17. `template_floating_object_submersion_depth` — hardcoded fallback injects "Pine Wood" at 500 kg/m³, not a `MATERIAL_DENSITIES` entry; `length`/`width`/`radius` are stated as givens and never used.
18. `template_statically_indeterminate_shaft` — `diameter`, `material_name`, `shear_modulus_gpa` stated in the question and never used (they cancel). A `depends_on` graph leaves three declared inputs with no consumer.
19. `template_pfr_volume_changing_rate` — prints the integrator's own error estimate as a step value (machine-dependent); `n` rounds to exactly 2.0 in ~10% of draws, contradicting the printed "non-integer order" note.

**Metadata / consistency:**

20. `template_qr_policy_one_iteration` — labelled Advanced in its section comment, Intermediate in its docstring.
21. `SHEWHART_K_SIGMA` is centralised in `attributes_control_charts.py` and hard-coded as a literal `3` in `variables_control_charts.py`.
22. Narrow answer spaces: `template_arl_beta_mean_shift` 18 distinct answers in 200 seeds; `template_mm1_time_in_system` 24 (the latter undisclosed in its docstring).

---

## 11. Cost estimate and recommendation

### 11.1 Per-class cost

| Class | n | What it needs | Effort |
|---|---:|---|---|
| **A** | 49 | Return the already-computed values alongside the prose. Values are named, ordered, and mostly already rounded to display precision. | ~0.5–1 h each. **Near-mechanical.** |
| **B** | 34 | Same, plus the emitted schema must be per-instance: the branch taken determines the milestone set, formulas and thresholds. ~10 also need the branch variable *stated* in the trace. | ~2–4 h each. |
| **C** | 51 | A typed milestone model: `value` must admit scalar, vector (component list + sign semantics), sequence-with-origin, categorical label, and opaque symbolic expression; plus partial-credit semantics for multipart. Most of the cost is in the **model design**, ~15–20 h once, then ~1–2 h per template. | ~1–2 h each after the model exists. |
| **D** | 16 | Case by case. ~8 are instrumentable once a defect is fixed (the back-solved templates, the rounding bug). ~5 need a different trace representation (iteration logs, search logs). ~3 have no numeric content at all and need a decision about whether they belong. | ~4–8 h each, **plus a scoping decision**. |

**Template-side subtotal:** roughly 200–320 hours including the one-off milestone-model design.

### 11.2 The costs outside the templates — larger than the template work

| Item | Effort | Note |
|---|---|---|
| **Value extractor rewrite** | ~60–100 h | Must return (value, unit, context) tuples for every quantity in free prose, with unit normalisation and scale handling. From scratch; the current 118 lines are not a foundation. **This is the critical path.** |
| **Unit model** | ~20–40 h | Normalisation, prefix handling, dimensional algebra. Must accommodate symbolic units (`kN·m³/EI`) and folded scale factors (`ppm`) or explicitly exempt them. |
| **Fixing the 9 self-inconsistent templates** | ~20–30 h | Non-negotiable. A verifier cannot be built against gold traces that contradict themselves at up to 7.77%, or that flip a classification on 35% of instances. |
| **Re-run affected evaluations** | — | Independent of this project; see §11.4. |

**Total: roughly 300–500 hours.** The template instrumentation the prompt was worried about is the *smaller and safer* half.

### 11.3 Recommendation: proceed with a modified design

The hypothesis is sound and the corpus supports it better than expected. But three assumptions in the framing must be dropped:

1. **Drop static milestone declaration.** 48% of templates branch structurally. The generator must emit the trace **per instance**. This is not negotiable and not repairable.
2. **Drop "milestone = scalar with a unit."** 61 of 150 templates carry non-scalar content the schema must represent. `is_milestone` needs a companion `kind ∈ {numeric, symbolic, classification, narrative, check}` and a `required` flag; `value` needs to be a typed union.
3. **Drop "no LLM in the critical path" as a universal claim.** ~13% of the corpus (symbolic answers, classifications, sequence-with-origin) cannot be verified numerically *in principle*. The defensible claim is stronger and more honest than the current one: **"deterministic numeric and dimensional verification for N of 150 templates, with the remainder explicitly scoped and separately handled."** With ~130 templates deterministically verified, that is a genuine answer to the reviewers' objection — and being explicit about the residue is more credible than claiming universality.

**Sequencing matters, and it is the opposite of what the project name implies.** Do the extractor first. Instrumenting 150 templates while the model side still reduces to one unlabelled float produces a beautifully structured gold trace with nothing to compare it against. A useful order:

1. Fix the parser defects and **re-run the affected evaluations** (§11.4). Cheap, independent, and it removes a live error from the results.
2. Build and validate the value extractor against the 2,200 archived model traces. **This is the gate.** If it cannot reliably recover (value, unit) tuples from real model prose, the whole approach fails here and the remaining work is not worth starting.
3. Design the milestone model against the class C and D templates specifically — the hard cases must inform the schema, not be retrofitted to it.
4. Instrument class A as a pilot (49 templates, near-mechanical) and measure end-to-end verification quality before committing to B, C and D.

### 11.4 One finding that should not wait

The thousands-separator bug corrupts the gold reference value for two templates in 100% of instances and produces confirmed false negatives in the archived results. It is a two-line fix (strip separators *before* matching, not after). Independent of everything else in this report, **the affected evaluations should be re-run before any further submission** — this is a correctness issue in published numbers, not a design question.

---

## 12. Uncertainties and what would settle them

Stated rather than smoothed over.

1. **The corpus-wide magnitude of the parser bug is unknown.** The archive is failure-filtered, so I can show ≥43 confirmed false negatives but not the true rate. **Settled by:** re-running the evaluation on the full 1,350-item-per-model set with a fixed parser and diffing the scores. This is also the single cheapest high-value action available.

2. **Whether a rewritten extractor can actually recover units from model prose is untested.** I established that the *current* parser cannot (it has no unit handling at all). I did not build a better one and measure it. My assessment that set-membership matching is tractable is an inference from the structure of the failures, not a measurement. **Settled by:** a spike — build a (value, unit, context) extractor, run it on the 2,200 archived traces, and hand-score a stratified sample of 100. Until that exists, the feasibility of the whole project is genuinely open, and I would not commit to template instrumentation before it.

3. **`n_milestone_candidates` in the CSV is a measured proxy, not a hand-count.** It is the mean number of result-position numbers per instance (the value after the last `=` on a line, plus answer-line numbers). It over-counts operand restatements and under-counts quantities emitted without an `=`. Treat it as an ordering signal, not an exact count. The per-branch audits hand-counted these for many templates and their figures should be preferred where the two disagree.

4. **`answer_type` and `unit_system` in the CSV are hand-verified for the ~90 templates the branch audits examined in detail and automated for the rest.** The automated classifier under-detects: it initially missed several vector and classification cases that source inspection found. Rows whose `notes` field is generic are more likely to carry an automated classification. **Settled by:** a targeted second pass over the ~60 templates with short notes.

5. **Effort estimates are unvalidated.** They are analogy-based, from the structure of each template, with no instrumented template to calibrate against. The class A pilot in step 4 above would calibrate them; I would expect class A to hold and class C/D to be optimistic.

6. **Whether the 9 self-inconsistent templates are worth fixing or should be replaced** is a scoping call I cannot make. Six of them back-solve their inputs, which is a deliberate authoring technique that produces good questions; fixing them means changing how the questions are generated, which may change their difficulty characteristics. **Settled by:** a decision from whoever owns the benchmark's item pool, informed by whether those items' difficulty statistics are load-bearing in the paper.

---

## Appendix — method

- **Static:** AST sweep of all 47 template modules classifying every `FormattedValue` interpolation as bound (`Name`, `round(name, k)`, `Attribute`, `Subscript`, conditional over two bound values) or computed inline (`BinOp`, non-wrapper `Call`, mixed conditional). 6,166 interpolations across 150 functions.
- **Dynamic:** each template executed under `sys.settrace`, capturing the generator's frame locals at return; locals recursively flattened (including nested lists, dicts, `ndarray`, `Fraction`, `Decimal`, complex) plus module globals; every result-position number in the emitted solution matched against that value set as EXACT / ROUNDED (within the printed precision) / SCALED (a bound value × 10ᵏ) / MISSING. 15 seeds per template; 21,310 result-position tokens.
- **Structure:** 25 seeds per template, diffing the numeric-blanked solution skeleton and step count across seeds to detect instance-dependent structure; branching then hand-confirmed against source.
- **Extraction:** `evaluation/engineering_parser.py` run against 2,200 archived model traces (11 models × 200) over 1,068 distinct gold traces from `error_analysis_annotation/samples/`, plus targeted probes on constructed answer strings.
- **Source review:** five parallel per-branch source audits, each reading all 30 templates in its branch and generating instances (40–300 seeds per template) to confirm findings against real output.
- Classifications in `template_inventory.csv` merge the instrumented measurements (base) with the hand-verified branch findings (overrides); all 150 rows carry a hand-written note.
