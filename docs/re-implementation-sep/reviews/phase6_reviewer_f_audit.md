# Reviewer F — audit replication (Phase 6)

**Filed:** 2026-09-16 · **Frozen SHA:** `381f178` · **Protocol:** R0–R6
**Committed unmodified and before any fix**, per the independence rule.

**Gate as commissioned:** not *"does it replicate the original?"* — the original harness never
existed — but *"is the NEW method sound, and does it classify the UNCHANGED templates the way the
old one did?"*

---

## 1. GATE ANSWER

**The new method is sound and it classifies the unchanged templates the way the old one did — but
two of the three narrative claims that defend the moved cells are wrong in their stated mechanism,
and one of those is the Phase 6 headline gate.** Calibration reproduces exactly (6166/6166
FormattedValue, chemical 65 / civil 37 / mechanical 89 / electrical 127 / industrial 131, all +0
against the corpus at 9105317, which I confirmed is the audit's own commit); holding the corpus at
the rev is the correct call and hides nothing, because Section B of the same output prints the HEAD
drift beside it. Nothing regressed silently — a fresh `run --checks all` reproduces all eight
ceilings and `ci_ratchet` exits 0. The 107-cell diff is real and correctly counted. However: the
`78 -> 50` move on `template_system_property_linearity` is **not** a regression (the prompt's worry
is answered in the negative) but it is **also not** the D6.7 predicate artefact the commit says it
is — I measured the pre-D6.7 source through the identical harness and got a byte-identical 50%, so
D6.7 is provably irrelevant to that cell; and the class-D gate was **not** shown to be mechanically
undecidable, because limb 2 of class D is exactly what `t1_closure.py` checks, by that check's own
docstring, and it was never applied to the D rows. Right conclusions, unsupported evidence — which
for a gate that rests on its own honesty is the finding.

---

## 2. FINDINGS

### F1 — CONFIRMED. Calibration reproduces; calibrating at the rev is correct and hides nothing.

```
PYTHONIOENCODING=utf-8 python -m tests.template_integrity.regen_inventory --calibrate-only   # EXIT=0
```

All five branches +0, population 6166 vs 6166, `GATE: PASS`.
`git log -1 --format='%H %ad %s' 9105317` → `Sat Sep 5 19:37:32 2026 … Add template structure
audit: report and per-template inventory` — the rev *is* the audit. Holding the corpus there is the
only way to separate method drift from six phases of template change, and Section B reports the HEAD
drift in the same breath (+121 interpolations, electrical +34). No concealment.

One reader-trap worth noting: the audit carries **two** per-branch inline tables.
`template_audit_report.md:85` (chemical 65, civil 37, electrical 116, industrial 125, mechanical 89,
total 432) is the one the gate matches. `template_audit_report.md:69` gives a *different* set
(mechanical 90, industrial 92, electrical 63, chemical 54, civil 27). Nothing in the calibration
output says which table it is quoting, so a reviewer checking against `:69` would conclude the gate
failed on all five branches. The strict total 449 vs the audit's 432 (+3.9% more inclusive
predicate) is disclosed at `phase0_baseline.md:43`.

### F2 — CONFIRMED. 107 cells, as claimed.

```
old rows 150  new rows 150   TOTAL DIFFERING CELLS: 107   rows touched: 74
  pct_step_values_recoverable 48 | n_inline_computed 40 | pct_answer_values_recoverable 17 | n_steps 2
```

### F3 — CONFIRMED, and this is the main finding. `78 -> 50` is not damage, but the commit's reason for saying so is a non-sequitur, and D6.7 did not cause it.

Two independent problems.

**(a) Category error.** `pct_step_values_recoverable` is a purely *dynamic* column —
`regen_inventory.py:441-449` computes it by generating instances, extracting `step_tokens()`, and
classifying each through `match_value()`. The Call-node-vs-Name-node argument is a purely *static*
AST fact that feeds only `n_inline_computed`. It cannot move the dynamic column at all, as stated.

**(b) Measured refutation.** I executed the pre-D6.7 source
(`git show f333f1a^:…/discrete_time_signals.py`) through the identical harness:

```
PRE-D6.7  tokens= 157 recovery=50%  {'MISSING': 79, 'EXACT': 33, 'SCALED': 45}
HEAD      tokens= 157 recovery=50%  {'MISSING': 79, 'EXACT': 33, 'SCALED': 45}
```

Identical. **D6.7 moved this cell by zero.** The whole `78 -> 50` gap is
harness-vs-2026-09-05-audit, not collateral damage.

**The real mechanism**, from the MISSING contexts: `NUM_RE = [-+]?\d[\d,]*…` harvests *index digits
out of symbolic identifiers*. The dominant unrecoverable lines are
`T{x1[n] + x2[n]} = T{x1[n]} + T{x2[n]}`, `y1[n] = (x1[n])^3`, `x3[n] = x1[n] + x2[n]` — the 1/2/3
are subscripts and exponents, not quantities. This template is mostly symbolic, so the "first number
after the last `=`" rule mostly harvests non-values. Nothing to do with D6.7, and not a trace-quality
statement at all.

The *static* half of the claim does hold for this row: `{2*C}`×2 (already BinOp/inline) plus `{C}`×2
(Name/bound) became four `signed_term(...)` Calls, giving `n_inline_computed 2 -> 4` exactly as
described.

### F4 — CONFIRMED. The D6.7 attribution table is incomplete; "The correspondence is one-to-one" is unsupported.

```
git show f333f1a --numstat
```

D6.7 edited **ten** template modules. The commit's table lists **six**. Omitted:
`discrete_time_signals.py` (11 lines — the *largest* single change), `continuous_time_signals.py`
(7), `volumetric_properties_pure_fluids.py` (3), `waves_and_phasors.py` (1). The omitted
`discrete_time_signals.py` is precisely where `template_system_property_linearity` lives — the
template the commit then spends a paragraph on.

### F5 — CONFIRMED. The class-D gate claim overstates the judgement limbs by one, and `16 -> 16` is an identity, not a measurement.

`template_audit_report.md:28` limb 2 is *"the trace's own chain does not reproduce its own answer."*
`tests/template_integrity/checks/t1_closure.py` docstring, verbatim:

> `T1 - Printed-arithmetic closure. … This is the check that catches a trace which does not
> reproduce its own answer`

`phase6_residual_register.md:17` calls limb 2 *"printed-arithmetic closure over a hand-read chain"* —
it names T1's own mechanism and then files it under human judgement, while T1 runs in the same suite
and reports 29 failures. So limb 1 is mechanical, limb 2 has a dedicated mechanical check in this
repo that was never applied to the D rows, and only limb 3 ("no stable symbols") is genuinely
judgement. **One of three, not two.**

Separately, `classify_restated()` (`regen_inventory.py:498`) opens
`if inherited_class == 'D' or zero_numeric: return 'D'` — it is structurally incapable of removing a
row from D. And all 16 D rows carry limb `D:inherited-judgement`; **not one** is held by the
mechanical `D:no-numeric-content` limb. So `CLASS D: 16 -> 16` is a tautology of the function's
design, presented in the commit message as a result line beside measured figures. The register's
docstring does admit the incapacity; the commit message does not.

Verdict on the prompt's question: not an outright excuse — the gate genuinely is *partly* judgement
and the register is unusually candid about it. But it was not shown non-decidable; it was shown
un-attempted on the one limb this repo could have decided.

### F6 — CONFIRMED. The "34 of 34 rows, zero hard signal" figure does not reproduce.

Two counts from the committed artefacts, neither of which is 34/34:

- `--diagnose`: **37** newly-"yes" branching rows, 37 prose-only, 0 hard signal.
- The regen CSV's 36 row-moves: 28 on the B limb, of which **26 prose-only and 2 on a hard
  `step-count` signal** (the other 8 moves are `C:answer_type` 6, `A:residual` 2).

So "zero hard signal" is **false** for the move set. The rhetorical point survives at 26/28 = 93%,
but the specific figure quoted in both `a142d27` and `phase6_residual_register.md:23` cannot be
derived from either artefact. Also of note: the 28 B-moves come from original classes A (22) and
C (6) — none from D.

### F7 — CONFIRMED. Nothing regressed silently.

```
PYTHONIOENCODING=utf-8 python -m tests.template_integrity.run --checks all --json <p>   # EXIT=1 (by design, failures present)
  T1: 29  T2: 0  T3: 0  T4: 0  T5: 66  T6: 142  T7: 83  T8: 0
PYTHONIOENCODING=utf-8 python -m tests.template_integrity.ci_ratchet --report <p>       # EXIT=0
  PASS: every check at its recorded ceiling.   (all eight "at ceiling")
```

Exact match to the recorded ceilings on a fresh 85s / 150-template run. No unclaimed regression
found.

### F8 — PLAUSIBLE, unclaimed, not a regression but an unguarded surface.

The same run reports: *"T1 marginals: 97 templates carry marginals (within tolerance, but on the
rounding boundary). Worst: critical_depth_froude_classification 76%, influence_line_max_reaction
74%, upward_seepage_quick_condition 69%."* Two-thirds of the corpus sits on a tolerance boundary and
this number is **not ratcheted**. A drift there would not trip the gate. Nobody claimed it either
way; flagging it as the largest unguarded surface in the suite.

---

## 3. COULD NOT CHECK

- **The other five templates in the attribution table** were not A/B'd pre/post-D6.7. Only
  `template_system_property_linearity` was. Given that the one I tested moved by zero, the remaining
  five deserve the same treatment before the "predicate artefact" framing is relied on.
- **The four files omitted from the attribution table** (`continuous_time_signals`,
  `volumetric_properties_pure_fluids`, `waves_and_phasors`, and the rest of `discrete_time_signals`)
  were not individually checked for `pct_step` movement.
- **`phase5_contract_scan`** was supplied but not run — it belongs to the fenced-off doubled-sign
  work.
- **`n_milestone_candidates` +27.1%** (1144 → 1454, reported by `--diagnose`) was not independently
  assessed. The audit itself (`:497`) calls this column a proxy and says hand-counts should be
  preferred, so the drift may be immaterial — but it is 1 of the 4 changed columns' worth of
  narrative that nobody has tested.
- All of the above: time box.

---

## 4. FORWARD-LOOKING SUGGESTIONS (triage separately from findings)

1. **Rewrite the attribution paragraph in `a142d27`'s successor doc.** List all ten edited files,
   and either drop the `78 -> 50` example or replace its explanation with the symbolic-token
   mechanism. The A/B that settles it is ~30 lines and now exists.
2. **Add an `--at-rev` mode for the *dynamic* columns** (temp checkout + re-execute) so "did this
   commit move this cell" becomes a command rather than an argument. This is the single
   highest-leverage tooling gap: the static side already has `measure_static_at_rev`, the dynamic
   side has nothing, and that asymmetry is exactly what produced F3.
3. **Either apply T1 to the 16 D rows and publish a per-row limb attribution, or restate the
   register's claim** as "one limb is judgement; limb 2 has a mechanical proxy we chose not to
   apply, and here is why." The current phrasing will not survive a reviewer who opens
   `t1_closure.py`.
4. **Recompute or retract "34 of 34."** Publish the 26-of-28 breakdown and name the 2 hard-signal
   rows.
5. **Consider ratcheting the T1 marginal count (97).**
6. **Consider excluding symbolic-identifier index digits from `step_tokens()`**, or adding a
   symbolic-coverage column. Templates like `system_property_linearity` have their recovery
   percentage dominated by subscripts and exponents rather than quantities, which makes the column
   misleading precisely where it looks worst.
