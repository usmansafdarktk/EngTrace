# Phase 0 — Summary and close-out

**Status: COMPLETE WITH FINDINGS. Exit gate partially met — two items open (§5).**
**Date:** 2026-09-05 · **Merged to `master`** at `23a520e`

Phase 0 built the verification infrastructure and re-measured the audit's
claims. **No template was edited**, by design.

---

## 1. What Phase 0 was for

Make every later phase falsifiable. Its own gate is a single property:
*a harness that cannot detect a planted defect cannot certify a fix.*

It earned its cost twice over before Phase 1 began — see §3.

---

## 2. Deliverables

| # | Deliverable | Where | Status |
|---|---|---|---|
| D0.1 | Harness, seven checks | `tests/template_integrity/` | ✅ |
| D0.2 | Baseline, 1,000 instances × 150 templates | `tests/template_integrity/baseline/profiles.json` | ✅ |
| D0.3 | Baseline measurement report | `phase0_baseline.md` | ✅ |
| D0.4 | CI entry point | `tests/template_integrity/run.py` | ✅ |
| D0.5 | Defect-rate reconciliation | `phase0_baseline.md` §1 | ✅ |
| D0.R | Independent review + triage | `reviews/phase0_harness_adversary.md`, §4 below | ✅ |
| — | T2 oracles (9) + findings | `tests/template_integrity/oracles/`, `reviews/phase0_oracle_findings.md` | ✅ |

**Runtime: ~17s for all 150 templates** (T1/T2/T4/T5/T7 in 9s, T3 in 8s). This is
the return on Phase 0 — verification that would take hours by hand is a CI step.

```
python -m tests.template_integrity.run --checks T1,T2,T4,T5,T7
python -m tests.template_integrity.run --checks T3 --det-seeds 200
python -m tests.template_integrity.run --checks T6            # vs baseline
```

---

## 3. The two findings that justified the phase

**A wrong fix recipe, caught before it was applied nine times.** P2's own worked
exemplar — the template the spec cited as carrying the fix Phase 1 should copy —
still fails printed-arithmetic closure on **5.89%** of its instances. Rounding to
display precision removes the *double* rounding but leaves an exact half-way tie
that Python resolves on the binary float (`0.01185 * 1000` → `11.8`, not `11.9`).
Applied as originally specified, Phase 1 would have reproduced that residue eight
more times. See [`DECISIONS.md`](DECISIONS.md) D-012.

**Two checks that were green while measuring nothing.**

1. `seed_all()` seeded numpy, so both of T3's processes reproduced the unseeded
   `np.random` draw and T3 **passed** the one template it exists to catch.
2. A colon-prefixed prose label poisoned T1's line parser, so
   `Row sum: 0.4347 + 0.2952 + 0.2701 = 1.5000 = 1` — a gross non-closure —
   passed the entire suite.

Both fixed and verified. See D-015.

> **Carry this forward:** a green check is not evidence until something
> known-broken has been shown to turn it red. This is the whole argument for the
> planted-defect gate, and for not compressing Phase 0.

---

## 4. Review outcome and R4 triage

**Harness Adversary: PASS WITH FINDINGS.** 11 of 12 planted defects detected by
the check that owns them; the unmutated control passed everything cleanly
(T6 matched the committed baseline exactly), so no detection is a scaffolding
artefact.

| Finding | Disposition | Action |
|---|---|---|
| **F1** T1 blind to prose-labelled lines (planted D11 missed) | `ADOPT-NOW` | Fixed: `:` is a clause separator. Planted defect now fails, correct version passes, corpus failures unchanged at 38. |
| **F2** T3 crashes on refs outside `data/templates/branches/` | `ADOPT-NOW` | Fixed: ref specs passed to the child via temp file (also resolves a Windows argv length limit), so candidate templates can be checked pre-commit. |
| **F3** Half-way tie lands MARGINAL, not FAIL | `ADOPT-PHASE-1` | Confirmed by design — the tolerance rule cannot distinguish it. T5b covers the class. Phase 1 must verify fixes against T5b, not T1 alone. |
| **F4** `run.py` ignored marginal density | `ADOPT-NOW` | Fixed: `marginal_rate` reported per template, summary line added, `--strict-marginals RATE` gate available. 95 templates carry marginals. |
| D0.5 **SPEC-2** P2 insufficient | `SPEC-CHANGE` | D-012. **Must land before Phase 1 edits any template.** |
| D0.5 **NEW-2/3** two more Phase-1-shaped defects | `ADOPT-PHASE-1` | D-011: Phase 1 grows 9 → 12. |
| D0.5 worst-case figures are 200-seed extreme values | `SPEC-CHANGE` | Spec should quote distributions, not sample maxima (§5 open). |
| Oracle: `poissons_ratio` sign leak | `BACKLOG` | Presentation, not arithmetic. Record in the residual-risk register (D6.6). |
| Oracle: `mean_variance` probabilities sum to 0.999–1.002 | `ADOPT-PHASE-1` | Renormalise as part of that template's fix. |
| Oracle: `rotating_unbalance` mass ambiguity | `ADOPT-PHASE-1` | Question must state whether the unbalance mass is included. |
| Oracle: answers printed to 1 s.f. (`0.002 mm`) | `ADOPT-PHASE-1` | Ungradeable independently of the chain defect. |
| `stoichiometry.py:115` invalid escape `\%` | `ADOPT-PHASE-5` | Warns today; hard `SyntaxError` under `-W error::SyntaxWarning`. |

---

## 5. Exit gate — where it stands

| Criterion | Status |
|---|---|
| T1–T7 green except the known-defective set | ✅ failures are the expected set |
| ≥5 planted defects detected | ✅ 11 of 12; the miss is fixed |
| D0.5 reconciliation, discrepancies >20% explained | ✅ 20 of 22 reproduce; 2 corrections recorded |
| Baseline committed and hash-stable | ✅ |
| Independent review filed and triaged (R2/R4) | ✅ |
| **P2 amended before any template is edited** | ❌ **OPEN** (D-012) |
| **Phase 1 re-scoped to 12 templates in the spec** | ❌ **OPEN** (D-011) |

Merged to `master` early at the repo owner's direction, with the two open items
recorded rather than implied closed. **Both must be resolved before Phase 1
edits any template**, since they change the fix recipe and its scope.

---

## 6. Corpus state as measured (unmodified templates)

Not a to-do list — the baseline Phase 1 onward is measured against.

| Check | Failing | Note |
|---|---:|---|
| T1 closure | 38/150 | 95 more carry marginals |
| T2 round-trip | 3/9 with oracles | see `reviews/phase0_oracle_findings.md` |
| T3 determinism | 1/150 | `levenspiel_plot_interpretation`, unseeded numpy |
| T4 contract | 4/150 | malformed markers; 5 more use `**Final Answer**` |
| T5 binding/rounding | 73/150 | 39 confirmed P2 violations |
| T7 asserts | 95/150 | expected — only civil and industrial carry asserts |

T7 is advisory corpus-wide and a gate only on templates a phase edits.

---

## 7. Correction to an earlier interim summary

An interim report of mine stated that the audit's claims "did not survive"
independent re-measurement and that Phase 1 would shrink. **That was wrong.**
D0.5 reproduced 20 of 22 claims, several exactly, and Phase 1 **grew** from 9 to
12. The oracles and D0.5 were answering different questions — *"can a solver
reach the gold answer?"* versus *"does a strict recomputation flip?"* — and I
collapsed the distinction. Corrected in D-013.

---

## 8. What Phase 1 inherits

1. **Amend P2 first** (D-012), then fix `beam_deflection_formula` first, since it
   is the pattern the others copy.
2. **Twelve templates**, in three categories with different fixes (D-013):
   chain breaks need P3; display defects need P2-as-amended; three of the
   original nine need nothing.
3. **Verify against T5b, not T1 alone** — T1 cannot see the half-way-tie class.
4. **Nine oracles already exist** and will re-run in 9s to confirm each fix.
5. **Compare at the precision the answer is quoted to**, or the tolerance widens
   until it tests nothing (`reviews/phase0_oracle_findings.md` §6).
