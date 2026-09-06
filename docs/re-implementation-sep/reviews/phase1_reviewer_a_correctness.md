# Phase 1 — Reviewer A: Correctness

**Filed:** 2026-09-06 · **Branch reviewed:** `redesign/phase1-round-trip`
**Roster:** A — Correctness / determinism (spec §Roster)
**Brief:** one mandatory gate task — *confirm each T2 round-trip oracle is a
genuine independent re-derivation, not a transcription of the template's own
arithmetic* — plus time-boxed secondary sampling. 35-minute box.

> The reviewer worked in isolation from Reviewer B and without sight of the
> implementer's reasoning, per R1.1. Code comments and commit messages were
> supplied as **claims under test**, not as evidence.

**Implementer's remediation is recorded in the right-hand column of §3 and in
[`../phase1_summary.md`](../phase1_summary.md) §7. All CONFIRMED findings are
closed.**

---

## 1. Verdict

**PASS WITH FINDINGS** — two CONFIRMED findings blocked the merge under §1.6.
The gate property itself passes: all four new T2 oracles are genuine independent
re-derivations, and every planted structural defect turned them red. The
blockers were elsewhere: a surviving T1 non-closure in `annulus_flowrate`, and a
`shaft_design_power` T2 tolerance below its own quantisation floor.

---

## 2. Independent re-derivation

| Acceptance number | Source of claim | Measured | Agreement |
|---|---|---|---|
| T1 = 0, T2 = 0 on the phase's templates @ 1,000 seeds | §1.6 exit gate | T1 0, T2 0, T3 0, T4 0, T7 0; T5 5 failing, T6 2 failing | **Agrees.** T5/T6 residue is D-017 / D-019, already recorded. |
| `statically_indeterminate_shaft` worst rel. error 1.35e-3 @ 1,000 | oracle `TOLERANCE` comment | 1.347e-3; grows to 1.568e-3 @ 20,000 (still inside 5e-3) | **Agrees exactly.** |
| `composite_shafts_series` J-display error ~5e-7 per segment | oracle comment | worst rel 9.93e-5 @ 1k, 9.997e-5 @ 20k — dominated by the 5-sig-fig display, not J | **Consistent.** |
| `beam_deflection_formula` oracle/template agreement | implied by T2 green | worst rel **exactly 0.000e+00** over 20,000 seeds | Agrees — see F-3. |
| Closure verified at 60,000 seeds | D-016 | 20,000 seeds × 6 display-defect templates | **Diverges — F-1.** |

20,000-seed sweep (reproduces F-1 and F-2):

```bash
python -c "
from tests.template_integrity.core import discover, generate
from tests.template_integrity.checks import t1_closure
import tests.template_integrity.checks.t2_roundtrip as T2
for tid in ['template_beam_deflection_formula','template_cantilever_double_integration','template_annulus_flowrate','template_composite_shafts_series','template_statically_indeterminate_shaft','template_shaft_design_power']:
    ref=[r for r in discover() if r.template_id==tid][0]
    inst=[generate(ref,s,capture=False) for s in range(20000)]
    r=t1_closure.run(inst,tid); r2=T2.run(ref,inst,T2.available_oracles())
    print(tid,'T1 fail',len(r.failures),'| T2 fail',len(r2.findings),'worst',r2.worst_rel_error)
    for f in r.failures[:3]: print(f)
"
```

```
beam_deflection_formula          T1 fail=0  | T2 fail=0 worst=0.000e+00
cantilever_double_integration    T1 fail=0  | T2 fail=0 worst=0.000e+00
annulus_flowrate                 T1 fail=6  | T2 fail=0 worst=5.749e-05
composite_shafts_series          T1 fail=0  | T2 fail=0 worst=9.997e-05
statically_indeterminate_shaft   T1 fail=0  | T2 fail=0 worst=1.568e-03
shaft_design_power               T1 fail=0  | T2 fail=1 worst=1.096e-03
```

### The mandatory gate: are the four oracles independent?

Method: each template's function source was copied into a scratch module, string-
mutated, and compiled with `optimize=1` so the template's own T7 asserts are
**stripped** — otherwise the invariants, not the oracle, would catch the defect.
Nothing under `data/templates/` or `tests/template_integrity/` was modified.

**Q1 — re-derivation or mirror?**

- **`beam_deflection_formula`** — genuine. Writes `5wL⁴/384EI` and `PL³/48EI`
  forward from the *stated* load; the template back-solves the load from a
  sampled δ/L target and works in kN·m⁻². Reuses no template intermediate.
- **`shaft_design_power`** — genuine for the core chain, in base SI. Notably does
  **not** replicate the template's carried roundings (`_hu(torque, 2)`,
  `_as_printed(c_cubed, ".3e")`). It *does* replicate the 2-dp RPM→Hz
  conversion — a borrowed intermediate, since the question states RPM, not Hz.
  Defensible under P3, and tested below (S7).
- **`statically_indeterminate_shaft`** — the strongest of the four. Solves the
  2×2 system by elimination to `T_A = T·L_BC/L`; the template walks the
  substitution route through a 3-dp ratio. Structurally different derivations,
  and the oracle correctly ignores `d` and `G`, which are distractors that cancel.
- **`composite_shafts_series`** — genuine. Forms `c`, `J`, `φ` at full precision
  per segment and sums; reuses none of the template's bindings.

**Q2 — would a wrong template turn them red?** Yes, in every structural case.
Rates are over all instances; a branch-specific defect is only reachable on its
branch.

```
=== beam_deflection_formula (tol 1e-3, n=300) ===
  SI point uses 5/384 (branch swap)        fail  25.0% (=100% of branch) worst 6.13e-1
  US uniform: dropped kip/ft -> kip/in     fail  22.7% (=100% of branch) worst 9.17e-1
  SI: wrong I unit exponent (1e-6 -> 1e-5) fail  51.3% (=100% of SI)     worst 9.04e+0
  answer +1%                               fail 100.0%
=== shaft_design_power (tol 1e-3, n=300) ===
  radius reported as diameter              fail 100.0%  worst 1.00e+0
  dropped kW -> W                          fail 100.0%  worst 9.23e+0
  RPM read as Hz                           fail  50.0% (=100% of branch) worst 2.92e+0
  dropped the 2 in c^3 = 2T/(pi*tau)       fail 100.0%  worst 2.61e-1
  f rounded to 0 dp instead of 2 dp        fail  40.5%  <- the replicated step IS probed
=== statically_indeterminate_shaft (tol 5e-3, n=300) ===
  L_AC/L_BC inverted (T_A <-> T_B)         fail  99.0%  worst 2.96e+0
  divisor missing the +1                   fail 100.0%
  T_A +1%                                  fail 100.0%
=== composite_shafts_series (tol 2e-3, n=300) ===
  stated diameter used as radius           fail 100.0%  worst 1.50e+1
  J = pi*c^4 (pi/2 dropped), segment 1     fail  98.3%  worst 6.86e-1
  only segment AB counted                  fail 100.0%
  no-op control                            fail   0.0%  worst 9.19e-5
```

The two sub-100% results are inherent, not oracle weakness: the T_A/T_B swap is
a no-op when `L_AC ≈ L_BC`, and halving `J₁` moves the total by under 2e-3 on
1.7% of instances where segment BC dominates.

**S7 settles the one independence concern.** Even though the oracle replicates
the 2-dp RPM→Hz step, mis-rounding *that step in the template* still turns the
oracle red on 100% of the RPM branch, because the oracle re-derives the rounding
from the question's stated RPM rather than reading the trace's printed Hz. It is
reconstructing a convention, not transcribing a value.

**Q3 — smallest real error each tolerance lets through.** Injecting a uniform
relative error into the answer only:

```
injected     beam      shaft_design   indeterminate   composite
             (1e-3)    (1e-3)         (5e-3)          (2e-3)
 +0.05%       14.7%       0.3%           0.0%            0.0%
 +0.1%        39.0%      51.3%           0.3%            0.0%
 +0.2%        72.7%     100.0%          12.3%           39.0%
 +0.5%        94.3%     100.0%          68.3%          100.0%
 +1.0%       100.0%     100.0%         100.0%          100.0%
```

The declared `TOLERANCE` is **not** the binding constraint on any of the four —
the display quantisation is. Smallest real error each reliably catches: beam
~1%, shaft_design_power ~0.2%, statically_indeterminate_shaft ~1%,
composite_shafts_series ~0.5%. Acceptable for the defect class Phase 1 targets
(tens of percent); **not** acceptable for a sub-1% P3 leak. See F-4.

**Q4 — is the quantisation legitimate?** Yes, in all four; none borrows the
answer. `composite_shafts_series` is the most defensible — its `places` is
computed from the *oracle's own* φ, so if the template's answer were wrong by a
decade the two would pick different precisions and the comparison would blow up
rather than agree. In all four, `gold_answer` reads only the `**Answer:**` block
by regex and `recompute` never touches the solution string.

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **F-1** | **CONFIRMED, blocking.** `annulus_flowrate` still fails printed-arithmetic closure — 6 T1 FAILs in 20,000 seeds (3.0e-4), all exact display ties on the `kappa` line: `kappa = 0.009000000000000001 / 0.0256 = 0.351562` evaluates to 0.3515625, `delta == tol` exactly. `_is_display_tie` was applied to `beam_deflection_formula` and `cantilever_double_integration` but **not** to this template's `kappa`. Invisible at the gate's 1,000 seeds (≈26% chance of appearing); certain at 20,000. Failing seeds: 3612, 8993, 9055 + 3. | **FIXED** — tie guard extended to `kappa` and the shape factor, radii bound through their displays. 0 failures in 20,000. |
| **F-2** | **CONFIRMED.** `shaft_design_power` `TOLERANCE = 1e-3` is below its own quantisation floor. Seed 11230 (10.3 kW, 99 Hz, 111 MPa, d = 9.125 mm): oracle 9.12, trace 9.13, rel 1.096e-3. Two *independently* quantised values at a boundary differ by a **full** step, not half, so the true floor is 1.1e-3 at d ≈ 9 mm — the tolerance comment is off by exactly 2×. Would report a phantom finding ~3 times in 60,000. | **FIXED** — by removing the ambiguity, per the reviewer's stated preference: the template now resamples when the givens-path and intermediates-path answers differ. Worst rel. error 0.00e+00 at 20,000 seeds. |
| **F-3** | **CONFIRMED (documentation).** The `beam_deflection_formula` oracle docstring asserts the trace states a 5-dp intermediate and that the paths differ by ~3e-5. The merged template states the **exact ratio** and uses it, which is why the agreement floor is exactly 0.0. The docstring is a claim under test and is false about the code as merged. | **FIXED** — docstring corrected. |
| **F-4** | **PLAUSIBLE.** The quantise-both-sides design caps T2 sensitivity at roughly one display step, and this is nowhere declared. Each `TOLERANCE` comment argues at length about a number that never binds. Recommend each oracle declare a **measured detection floor** — the smallest relative error caught on ≥99% of instances. | **ADOPTED** — `AGREEMENT_FLOOR` and `DETECTION_FLOOR` measured and declared in all four new oracles; `tests/template_integrity/oracle_floors.py` reproduces them. Two **template-level** declarations were found to be below their agreement floor and were corrected. |
| **F-5** | **PLAUSIBLE (demonstrated).** Oracle sensitivity is coupled to the constants tables Track B is about to change. Patching `SHEAR_MODULUS_VALUES` to 4 s.f. gives `composite_shafts_series` a worst rel. error of 1.977e-3 against a declared 2e-3 — passing with 1% margin, undetected. AISC SI `Ix` at 2 dp makes `beam_deflection_formula` fail 1.8%. Neither is a defect today (all G values ≤3 s.f., all 14 AISC SI `Ix` exactly 1 dp), but it is the D-015/D-018 failure mode arriving through the constants track. | **RECORDED as D-025**, assigned to Track B's sync points. |

---

## 4. Falsification attempts that failed

- **Oracle-as-transcription** (the assigned gate hypothesis). All four read line
  by line against their templates looking for shared expressions, rounding
  sequences or variable structure. Found exactly one borrowed intermediate, then
  showed it is still probed (S7, 40.5% detection). **The gate property holds.**
- **Structural defects passing an oracle.** 14 planted mutations; every one
  turned its oracle red on 98–100% of reachable instances. No structural defect
  hides.
- **Answer-borrowing.** No `recompute` path touches the solution text. The no-op
  control produced worst rel 9.19e-5 — the oracles genuinely disagree slightly
  with a correct template, so they are not tautologically equal.
- **`beam_deflection_formula` at 20,000 seeds** — 0 failures, worst rel exactly
  0. D-016's fix on the pattern template holds well past the written gate.
- **A US-branch span-conversion defect.** Not a valid test: changing `L_in`
  changes the question text too, so the oracle correctly follows the stated
  span. Recorded because the negative result is informative — **T2 verifies
  question→answer consistency and structurally cannot see a question that is
  internally mislabelled.** That is a gap in the check family, not in the oracles.

---

## 5. Further probing and improvements

Ranked by what the reviewer would probe next. Triage in
[`../phase1_summary.md`](../phase1_summary.md) §8.

1. **Closure at the seed count D-016 mandates.** 20,000 was run; D-016 says
   60,000. F-1 surfaced only above 1,000.
2. **T1 coverage as a gate.** 0.133 on `statically_indeterminate_shaft`, 0.182
   on `cantilever`, 0.193 on `beam` — 80–87% of `=` lines unchecked. The spec's
   own §0.2 says coverage below 60% is itself a finding; five of six templates
   are findings today and the run does not say so.
3. **The T1 marginal band.** 11,789–43,379 marginals per template at 20,000
   seeds. D-016 part 3 claims the density is what rounding alone predicts; that
   was not re-derived and it holds up the whole tie-removal argument.
4. **Per-oracle detection floors**, measured, committed alongside `TOLERANCE`.

**Does this weakness exist outside Phase 1?** (a) The exact-display-tie class is
a property of any template printing a quotient of grid-sampled quantities —
`annulus_flowrate` proves the D-016 sweep was applied per template rather than
as a corpus-wide pattern search. (b) The quantise-both-sides sensitivity ceiling
will be inherited by every Phase 2–3 oracle; Phase 2's `vibration_isolator_design`
and `pfr_volume_changing_rate` have small-magnitude answers where it is worst.
(c) `poissons_ratio`, `logarithmic_decrement` and `system_properties` left Phase 1
under D-013 as "not defective", but that judgement predates F-4;
`system_properties`' residue was dismissed as "two orders below its print
precision", which is exactly the regime T2 cannot resolve. **Do not cite them as
verified.**

**What later phases should do differently.** Size the acceptance run to the
defect rate *before* running it. Require every oracle to report a **measured**
detection floor from a planted-defect sweep — the same move D-015 asks for, now
applied to oracles rather than only to the harness. Gate the constants track
against templates' stated precisions before it lands.

**Spec issues.** (i) The phase is "12 templates" in the title and D-011, **10**
in §1.1a, and **9** in the §1.6 exit gate — D-011's 9→12 counts three arrivals
without subtracting three departures. (ii) §1.6 says 1,000 seeds while D-016
supersedes it with 60,000, so a phase can pass its own written gate while
failing the decision that governs it — which is what F-1 is. (iii) §0.2 T2 says
nothing about display quantisation, although all four oracles do it and are
right to; the practice and its consequence (F-4) belong in the spec rather than
being reinvented per oracle.
