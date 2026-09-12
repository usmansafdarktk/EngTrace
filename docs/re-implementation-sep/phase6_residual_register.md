# Phase 6 — residual-risk register (D6.6)

**Everything knowingly not fixed, and why.**

**Date:** 2026-09-13 · **Branch:** `redesign/phase6-consolidation`

A residual is not a failure. A residual that is *not written down* is. Each entry below states
what is wrong, what it would take, and who has to decide — so that the next phase inherits a
decision rather than a surprise.

---

## 1. The exit gate itself is not decidable by measurement

**`Class D reduced from 16 to ≤ 4` — measured 16 → 16, and no harness can settle it.**

Two of class D's three defining limbs (`template_audit_report.md:28`) are judgements: *"the
trace's own chain does not reproduce its own answer"* and *"a search/iteration log with no stable
symbols"*. A mechanical re-audit can therefore **inherit** D and never **clear** a row out of it.

The restated rule yields A 27 / B 56 / C 51 / D 16 and would move 36 rows — but it is deliberately
**not** written into `instrumentation_class`, because those moves are dominated by
`B:instance-branching`, the prose-only limb the harness measured as firing on **34 of 34** rows
with **zero** hard step or line-count signal. A sampled material name changes the blanked skeleton
without any governing-equation branch.

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

## 3. Three hand-rolled sign helpers remain

D6.7 promoted `signed_term`/`rect_str` into `data/templates/branches/_emission.py` precisely
because *nine hand-written sign fixes is nine chances to write `+ {value}` again*. Three
reimplementations survive:

| where | what |
|---|---|
| `waves_and_phasors.py` | its own `_signed_term` / `_rect_str`, now redundant beside the shared module |
| `discrete_time_signals.py:511` | `C_str = f"+ {C}" if C > 0 else f"- {abs(C)}"` |
| `discrete_time_signals.py:618` | `fmt = lambda c, v: f"+ {c}" if c > 0 else f"- {abs(c)}"` |

Retiring them should be **behaviour-preserving**, but "should be" is not a measurement: it needs
its own parity check, contract scan and P6 measurement. Left as a tranche, not folded into D6.7.

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

## 7. Track B — 177 `[UNVERIFIED]` rows stand

Mechanical 148, chemical 13, civil 8, electrical 4, industrial 4. Of the mechanical bulk: ~86 are
materials and mixtures no chemical database indexes, ~19 need a mechanical-properties handbook,
11 are blocked on paywalled standards (AISC, ASTM A992/A36, ACI 318-19, ASCE 7-22, IEC 60063,
MIL-A-8625, US Standard Atmosphere 1976).

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

## 9. Inherited, unclosed

- **`CP_PARAMS` origin vs verification** (Reviewer G, G-7) — unfixable without a citable
  Smith–Van Ness copy.
- **G F-5** (`@domain: none` vocabulary) and **G F-7** (NAVFAC manuals cited by zero tags).
- **The value-extractor spike (D-002)** — *"the gate for the whole project"*, still unwritten, and
  after Phase 6 the sequence has nothing else queued.
