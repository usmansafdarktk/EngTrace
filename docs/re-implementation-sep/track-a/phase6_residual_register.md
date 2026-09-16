# Phase 6 — residual-risk register (D6.6)

**Everything knowingly not fixed, and why.**

**Date:** 2026-09-13 · **Branch:** `redesign/phase6-consolidation`

A residual is not a failure. A residual that is *not written down* is. Each entry below states
what is wrong, what it would take, and who has to decide — so that the next phase inherits a
decision rather than a surprise.

---

## 1. The exit gate itself is not decidable by measurement

**`Class D reduced from 16 to ≤ 4` — measured 16 → 16, and no harness can settle it.**

**CORRECTED after Reviewer F (see [`phase6_reviewer_f_triage.md`](phase6_reviewer_f_triage.md)).**
This section originally claimed *two* of class D's three limbs are judgement. That was wrong by
one. Limb 2 — *"the trace's own chain does not reproduce its own answer"* — is exactly what
`tests/template_integrity/checks/t1_closure.py` checks, by its own docstring, and it had never been
applied to the D rows.

**Applied: all 16 class-D rows PASS T1, 0 failures each.** So limb 2 holds none of them, and
neither does the mechanical `D:no-numeric-content` limb (not one row carries it). The corrected
statement is therefore **stronger** than the original: *both mechanically decidable limbs clear all
16 rows, so every one sits in D on inherited judgement alone.*

`classify_restated()` (`regen_inventory.py:498`) also opens `if inherited_class == 'D': return 'D'`
— structurally incapable of removing a row — so the original `16 → 16` was a tautology of the
function's design presented as a result. A mechanical re-audit can **inherit** D and never **clear**
a row out of it.

The restated rule yields A 27 / B 56 / C 51 / D 16 and would move 36 rows — but it is deliberately
**not** written into `instrumentation_class`, because those moves are dominated by
`B:instance-branching`, a prose-only limb. **CORRECTED after Reviewer F:** the figure originally
quoted here — *"34 of 34 rows with zero hard signal"* — **does not reproduce from either artefact
and is retracted**. The measured counts are: `--diagnose` reports **37** newly-"yes" rows (37
prose-only), and the CSV's 36 moves carry **28** B-limb moves of which **26 are prose-only and 2
fire on a hard step-count signal**. So "zero hard signal" is false; the honest figure is **26 of 28
(93%)**. A sampled material name changes the blanked skeleton without any governing-equation
branch, which is the point the wrong number was standing in for.

**Decision needed:** restate the gate as something answerable — *"the rule restated and applied,
every row still in D named with the limb that holds it"* — or accept that the original figure was
never mechanically reproducible.

## 2. T6's baseline is stale corpus-wide

**142 of 150 fail, and that is one stale file, not 142 defects.** Carried by Phases 3 and 5 and
Track B, each correctly. D-043 forbids regenerating it to turn a gate green: *a baseline refreshed
by the phase it is meant to gate is not a gate.*

The ratchet (`ci_ratchet.py`) now holds it at its recorded ceiling so it blocks nobody while still
catching regression. **`run.py --baseline` has not been invoked in this phase and must not be.**

**Decision needed, and it is the owner's:** regenerate under authority, against a stated corpus,
with the movement examined — or leave it ratcheted.

## 3. One hand-rolled sign helper remains

**Corrected 2026-09-16 against the tree; this section said three.** Two were retired inside
Phase 6 itself (`phase6_summary.md` §4b) after this file was written, so the count here was stale
rather than wrong when written. Measured now:

| where | what | state |
|---|---|---|
| `waves_and_phasors.py` | its own `_signed_term` / `_rect_str`, 30 lines over 10 call sites | **retired** — the shared module is used |
| `discrete_time_signals.py:511` | `C_str` | **retired** — now `C_str = signed_term(C)`; safe because `C = randint(1,5) * choice([-1,1])` is never zero |
| `discrete_time_signals.py:618` | `fmt = lambda c, v: f"+ {c}" if c > 0 else f"- {abs(c)}"` | **kept, on a measurement** |

`fmt` is kept deliberately, not overlooked: it renders zero as `"- 0"` where `signed_term` gives
`"+ 0"`, and it is called with `C2`, a solved coefficient that *can* be zero. Retiring it is
therefore a P6 event, not a refactor.

D6.7 promoted `signed_term`/`rect_str` into `data/templates/branches/_emission.py` precisely
because *nine hand-written sign fixes is nine chances to write `+ {value}` again*. Retiring the
last one needs its own parity check, contract scan and P6 measurement.

## 4. `degenerate_product_derivation` — 3 templates, printed on every run

`impulse_response_from_lccde`, `standing_wave_formation`, `undamped_response_initial_conditions`.
A census, not gated, unchanged by this phase. It is recorded here because Phase 5's named failure
was `cross_pair` printing `false rejects 56` on every run with nobody putting it in a claim.

**Decision needed:** in scope as a defect, or a deliberate non-defect with a stated reason.

## 5. `blind_guess_floor` / `surface_model_heldout` — not reproducible

Only **5 rows** carry these columns and **none of the five has a reproducible derivation**:

- `adiabatic_flame_temperature` and `heat_of_reaction_formation` — unmentioned anywhere outside the CSV.
- `levenspiel_plot_interpretation` — carries **1.0000 / 1.0000**; a blind-guess floor of 1.0 on an
  `array` answer is not a coherent statistic, and the held-out figure has no derivation at all.
- `system_property_linearity` and `system_properties_memory_causality` — Reviewer B's own snippet
  **reproduces exactly** (0.5002 / 0.3402 at N=4000) while the CSV says 0.5008 / 0.3450, matching
  at no N tried. **The CSV contradicts the only recorded derivation it has**, and
  `phase4_summary.md:379` propagates the unsourced pair.

`heat_of_reaction_formation` also carries a held-out **below** its floor — a negative lift,
meaningless under D-057's "lift ≥ 40 pts" — and it and `adiabatic_flame_temperature` are
**scalar-answer** templates where a majority-class floor is as ill-defined as on an array.

**Decision needed:** recompute the column from a committed script, or drop it. It supports no gate
decision as it stands.

## 6. D6.11 — two predicates disagree about which units are undeclared

`DECLARED_UNITS` holds **103** entries, derived by `derive_units` under *always-present AND
invariant* over 12 seeds, using `unit_token` — the token ending the clause containing the last
value, read over the **whole solution**.

An independent probe using T6's `_UNIT_AFTER_NUM` over the **answer block** splits the 38
bound-but-undeclared templates as: **23 dimensionless** (nothing to declare), **10 whose unit
varies by seed** (per-*item*, not per-template — matching the canonical census exactly: 113 always
− 103 invariant = 10), and **5 that look invariant and undeclared**.

**Those 5 are not established as missing.** The two predicates differ, and derive's stricter
clause rule may be excluding them correctly.

**Decision needed, and it is a spec choice rather than a measurement:** which predicate defines
"carries a unit". Until it is made, D6.11's true size is *10 hard cases plus at most 5 easy ones*,
not the "113 templates" the spec implies. Six of the 38 are `multipart` and need one unit **per
part** — one unit applied to all `n` parts is wrong for at least `n−1`, and made 13 templates
reject their own gold.

## 7. Track B — 170 `[UNVERIFIED]` rows stand

**Corrected 2026-09-16 against the tree; this section said 177 (mechanical 148, chemical 13).**
Measured now, by the register's own predicate — a comment line whose first token is the tag:

| branch | then | now |
|---|---|---|
| mechanical | 148 | **142** |
| chemical | 13 | **12** |
| civil | 8 | 8 |
| electrical | 4 | 4 |
| industrial | 4 | 4 |
| **total** | **177** | **170** |

The difference is not new sourcing: **`e7cf8c4` deleted seven substance-defect rows**,
owner-directed — `Cork`, `Cork Board`, `Bamboo`, `Tellurium Mercury`, `Tungsten Hexafluoride` and
`Tetrabromoethane` from mechanical, `R-410A` from chemical. Every one had been deferred across two
phases because deleting a row shifts the keys after it and moves the item pool; the pool is being
rebuilt, so that cost went to zero. 177 − 7 = 170 reconciles exactly, and
`phaseC3_residual_register.md` §3 records what each row was.

Of the mechanical bulk: ~86 are materials and mixtures no chemical database indexes, ~19 need a
mechanical-properties handbook, 11 are blocked on paywalled standards (AISC, ASTM A992/A36,
ACI 318-19, ASCE 7-22, IEC 60063, MIL-A-8625, US Standard Atmosphere 1976).

Also open: **7 `[DERIVED]` constants still UNEXECUTED** (4 Shomate refits, `Air(g)`, 2
by-definition ceilings), and **4 `tol=basis=condition`** tolerances the resolver cannot size —
`Glycerine` worst, its tolerance absorbing a *composition* ambiguity.

## 8. The stale annotation fork — left, by decision

`templates_annotation/annotation_app/.../constants.py` carries pre-C3 values (R-12 `v_g 0.0268`,
R-22 `v_f 0.000845`). **Owner decision: leave it.** Recorded here only because D6.4 asserts a
corpus-wide constant state that this directory contradicts.

It is also the source of the `SyntaxWarning: invalid escape sequence '\%'` seen on every repo-wide
AST sweep: `stoichiometry.py:115` in the fork holds `\%` where the live template has the corrected
`\\%`. Harmless today, a syntax error in a future Python, and **not** a live-corpus defect.

## 9. Carried by a decision or a review, and never registered until now

**Added 2026-09-16.** Three items were decided or filed, assigned forward, and then reached no
register that anyone still reads. Each is documented — in `DECISIONS.md` or in a review — and each
stopped one step short of a list. They are stated here because the pilot freeze is cut against
this register, and an item that moves emitted text after the freeze costs a re-run.

### 9.1 `signal_operations` omits the `n = 0` origin on 16.4% of instances (D-050)

The template marks the origin only when the result's support contains it, so on the rest the
printed answer is a bare value list from which the origin cannot be recovered. The *item* is well
posed; the printed gold answer is not self-describing.

D-050 records the fix — widen the printed support, `y[n] = {…} for n = 1…4` — and calls it "a
Phase 5 or 6 scoping decision". Phase 5 did not take it; Phase 6 did not take it; no register
carried it. **Re-measured 2026-09-16: 329 of 2,000 seeds silent, 16.4%** — the rate D-050 quotes
at 4,000 seeds, unchanged.

**Decision needed:** widen the support (a P6 event on 16.4% of that template's instances), or
accept a gold answer the comparator cannot fully check and say so.

### 9.2 The supporting-quantity remedy for the two shortcuttable classification items (D-057 → D-066)

`system_property_linearity` and `system_properties_memory_causality` are 100% predictable from the
question surface against floors of 50.02% and 34.02%. Reviewer B proposed requiring the answer to
name *which* property failed, then withdrew it as a Phase 4 requirement because it would drop
recall to ~27%. `phase4_summary.md` §12 assigned it to Phase 5; **D-066 declined it** — correctly,
Track B's charter was `tests/comparators/` and the inventory only — and reassigned it to Phase 6.
Phase 6 mentions D-057 only under §5's column question, never the remedy.

Both items are still bound and still fully shortcuttable. The remedy changes what the item *asks*
and what its gold *answers*, so it needs the question, the gold and the rubric to move together.

**Decision needed:** change the two items, or record in the results table that they are scored on
an answer that carries 1.0 and 1.6 bits.

### 9.3 Seven `UNBOUND` entries rest on a truth predicate known to be wrong in that direction

Reviewer E's **E-12**, re-filed as **R2-F2** and never triaged. `gold_gold` calls a `MATCH` a false
accept iff the two spans differ *textually*; at the display-tolerance boundary that is wrong, so a
binding that obeys D4.1 §4.1 is recorded as over-accepting. Measured now — the binding table has
moved since Phase 5 closed, **127 bound / 23 unbound** against that summary's 120/30, with
`symbolic` at 6 of 9 rather than 1 after D6.10:

| unbound with this reason | count |
|---|---:|
| `decimation_aliasing_analysis` | 34 |
| `null_to_null_bandwidth` | 24 |
| `gauss_law_symmetric` | 4 |
| `signal_energy_power`, `truss_method_of_sections`, `vdw_solve_for_pressure` | 2 each |
| `pitzer_correlation_z` | 1 |

None has been audited against §4.1's tolerance, so "23 declined" is not yet a measurement.
`vdw_solve_for_pressure` is the one E separated out: its gold at one displayed decimal cannot
distinguish two substances, which is an item-design question (D-050's shape) and not a comparator
defect — it should not simply be re-bound.

**Decision needed:** re-audit the seven against the display tolerance, and send
`vdw_solve_for_pressure` to item design rather than back to `UNBOUND`.

## 10. Inherited, unclosed

- **`CP_PARAMS` origin vs verification** (Reviewer G, G-7) — unfixable without a citable
  Smith–Van Ness copy.
- **G F-5** (`@domain: none` vocabulary) and **G F-7** (NAVFAC manuals cited by zero tags).
- **The value-extractor spike (D-002)** — *"the gate for the whole project"*, still unwritten, and
  after Phase 6 the sequence has nothing else queued.
