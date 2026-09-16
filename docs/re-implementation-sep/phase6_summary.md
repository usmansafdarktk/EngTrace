# Phase 6 — consolidation and re-audit: what it found

**Date:** 2026-09-16 · **Branch:** `redesign/phase6-consolidation`
**Deliverables:** D6.1–D6.12, plus Track B's owner-directed deletions.

---

## 1. The headline is that the headline cannot be measured

The exit gate reads *"Class D reduced from 16 to ≤ 4"*. Measured: **16 → 16**.

Not because the work failed, but because **the gate is not decidable by a harness**. Two of class
D's three defining limbs (`template_audit_report.md:28`) are judgements — *"the trace's own chain
does not reproduce its own answer"*, *"a search/iteration log with no stable symbols"*. A
mechanical re-audit can **inherit** D and never **clear** a row out of it.

The restated rule yields A 27 / B 56 / C 51 / D 16 and would move 36 rows. It is deliberately
**not** written into `instrumentation_class`: those moves are dominated by `B:instance-branching`,
the prose-only limb the harness measured firing on **34 of 34** rows with **zero** hard step or
line-count signal. A sampled material name changes the blanked skeleton without any
governing-equation branch.

**Recommended restatement:** *the class rule restated, applied, and every row still in D named
with the limb that holds it.* That is answerable; the original is not.

## 2. The audit harness never existed

Not deleted — never written. Exactly two commits ever touched `template_inventory.csv`, and the
one that added it contains **four documentation files and no `.py` at all**. Pickaxe on the
distinctive column names returns only prose. No deleted blob, no dangling object, no stash. The
cause is in the commissioning prompt: *"Stay on `master`. Do not create branches or commits."*

**But the method survived, twice over** — in the report's own appendix, and far better in
`phase0_baseline.md:289–348`, where the Phase 0 adversary independently re-implemented it with
runnable snippets. So D6.1 was **assembly, not invention**, and cheaper than a rebuild rather than
dearer.

`regen_inventory.py` is now committed and calibrates **exactly**: population 6166/6166 and all
five branches at +0 against the corpus **held at rev `9105317`**. Calibrating at HEAD would have
failed (chemical 64, civil 34, mechanical 87) and the predicate would have taken the blame for six
phases of legitimate template change.

## 3. What was fixed

| | |
|---|---|
| **D6.7** | doubled sign **14 → 0** across 10 files; 150/150 clean, 0 generation errors on 60,000 instances |
| **D6.8** | T4 now asserts the answer span's **shape**, not just its marker — two terms, both measured |
| **D6.9** | `composite` binding: right-number-wrong-class now caught, 4 of 4 on hand-built cases |
| **D6.10** | 5 of 8 `symbolic` templates fixed; decided rate 11.1% → 66.7% |
| **D6.5** | CI built from nothing, tiered, with a **ratchet** instead of a gate |
| **D6.12** | one seeded generator replacing 20 that were unseeded, partial and uneven |
| **Track B** | 7 substance-defect rows deleted, owner-directed |

## 4. Three findings that outlast this phase

**A green detector does not mean a correct fix.** On `coulombs_law`, the obvious `signed_term`
substitution turns `sqrt(1^2 + -5^2 + -6^2)` into `sqrt(1^2 - 5^2 - 6^2)` — a sum of squares
become a subtraction of them — and `phase5_contract_scan` reports that **clean**. Measured. The
scan proves the defect is gone; it cannot prove the fix is right. Shape `(d)` (`paren_neg`) exists
because of this.

**Aggregates hide hidden constants.** D6.7's diff reports `lorentz_force` at q 279 / ans 223.
Those totals cannot tell you whether the sets coincide — they do not. Per instance: 204 changed
question *and* answer, 76 solution-only, and **19 changed the answer while the question stayed
byte-identical**. Those 19 need **re-scoring** and are invisible to anyone diffing questions.

**There was no green suite to promote.** D6.5 assumed one. Measured, 78 of 150 templates fail the
default battery and `run.py` exits 1 — pre-existing debt, plus T6's corpus-wide stale baseline. A
plain gate would be red forever, and the one pressure that creates is to regenerate the T6
baseline, which D-043 forbids. The ratchet gates on the **delta** instead.

## 4b. D6.11 and the cleanup, as executed

**The unit is a property of the item.** Two predicates disagreed about which templates carry a
unit; the **shipped** one (`derive_units`' `unit_token`) is adopted as authoritative, so D6.11's
real scope was never "113 undeclared" but **ten templates whose unit varies by seed** — which a
template-level declaration structurally cannot express. `generate_testset.py` now emits `unit`
**per record**: 1,722 units across 2,250 items.

**A limitation of the rule now standardised on:** `unit_token` sometimes returns a *symbol* rather
than a unit — `y^2` appears among the emitted values. Faithful to the predicate, misleading to a
grader. Registered, not absorbed.

**Two of three hand-rolled helpers retired, the third kept on measurement.**
`waves_and_phasors`' local `_signed_term`/`_rect_str` (30 lines, 10 call sites) and
`discrete_time_signals`' `C_str` are gone. **`fmt` stays**: it renders zero as `"- 0"` where
`signed_term` gives `"+ 0"`, and while `a1`/`a2` are guarded non-zero it is also called with `C2`,
a solved coefficient that can be zero. `C_str` was safe for the opposite reason, checked not
assumed: `C = randint(1,5) * choice([-1,1])` is never zero.

**Proven a no-op rather than asserted one:** 1,500 instances across the 10 affected templates,
two trees, **two separate processes**, compared by hash — **identical**. No P6 event.

**D6.12's pool:** 2,250 records, 150 templates × 15, 450 per branch, uniform, 42 files, 2,212
distinct questions. The 38 duplicates are a template property (`heat_of_reaction_formation` draws
`random.choice` over a 4-entry table), not a generator defect.

## 5. Mistakes made and caught

Recorded because a summary listing none is not a clean phase; it is an unexamined one.

- **`py_compile` is not verification.** I used `paren_neg` without importing it; compilation passed
  (an undefined name is a runtime error) and the template then failed **400 of 400** instances. The
  contract scan caught it because the scan *generates*.
- **I committed a false claim and withdrew it.** I asserted the unit census "does not reproduce"
  after reading only `phase5_summary.md`. It reproduces exactly — the predicate is written out in
  `derive_bindings.py:29–40`. An independent measurement is evidence about the corpus, not
  evidence that an existing instrument is broken.
- **I drafted a second false claim and caught it before committing.** The D6.7 message said the
  emission changes moved gold and `--regenerate` resolved the drift. `bindings.py` was
  **unchanged**: I had assumed drift and regenerated without first measuring whether any existed.
- **Two of my instructions to agents were wrong, and both agents caught them.** I said to build on
  `t5_binding._is_bound()` (the *lenient* variant — would have failed calibration for unrelated
  reasons), and I suggested stripping trailing unit words for D6.10 (would have **created** false
  accepts: `rad`/`deg` on one template, `Hz`/`kHz` on another).
- **Pattern-matching pointed at the wrong code three times** — `\%` (the hit was in the stale
  annotation fork, not the live template), then `+ {` and `- {`, which over-match because the
  defect is a *runtime* property. The detector's own evidence was the only reliable source.
- **Two near-misses on scope**, both caught by grepping before editing: `fluid_kinematics.py`'s
  `{A}xt + {B}y^2` lines belong to a template **not** among the 14, and changing `phi_str` would
  have corrupted `time_to_phasor` into *"The initial phase is + 45 deg"*.
- I guessed the instance-dump JSON schema wrong twice before printing it, and once wrote scratch
  files to `/tmp` instead of the session scratchpad.

## 6. What is left, and who decides

See [`phase6_residual_register.md`](phase6_residual_register.md). The decisions that are the
owner's, not the implementer's: restating the class gate; regenerating T6's baseline under
authority; whether `blind_guess_floor` is recomputed from a committed script or **dropped** (none
of its five rows has a reproducible derivation, and for two the CSV contradicts Reviewer B's own
snippet); and whether the three surviving hand-rolled sign helpers are retired now or later.

**After Phase 6 the sequence is empty** except the value-extractor spike (D-002), still unwritten
and still described as *"the gate for the whole project"*.
