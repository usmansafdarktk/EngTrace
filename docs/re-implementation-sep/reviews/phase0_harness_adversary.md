# Phase 0 — Harness Adversary review

**Scope:** does `tests/template_integrity/` actually detect deliberately planted defects?
**Branch:** `redesign/template-integrity` · **Date:** 2026-09-05
**Method:** 12 mutants generated from one real donor template, `template_two_step_transition_probability`
(`data/templates/branches/industrial_engineering/stochastic_operations/markov_chains.py`), each differing
from the donor by exactly one edit. Nothing under `data/templates/` or `evaluation/` was modified.
Mutant generator: `scratch/build.py`; detection driver: `scratch/detect.py` (both untracked).

---

## §1 Verdict

**PASS WITH FINDINGS.**

All 11 planted defects were detected, each by the check that is supposed to own it, and the
unmutated control passed every check cleanly — so none of the detections are artefacts of the
scratch scaffolding. The harness can therefore certify a fix in the classes it covers.

But it has a **hole large enough to hide a gross arithmetic error**. A 12th probe, planted after
reading T1's skip reasons, produced a solution line reading

```
Row sum: 0.4347 + 0.2952 + 0.2701 = 1.5000 = 1, as required for a probability distribution.
```

and **every check in the suite passed it.** T1 never evaluates any `=` line that has prose in front
of the arithmetic. That is Finding F1 below and it is the most important thing in this report.

---

## §2 Planted-defect detection table

Run `python scratch/build.py && python scratch/detect.py 200` from the repo root to reproduce
rows D1–D10; D3 and D11 have their own commands in §3.

| # | Planted defect | Check that should catch it | Detected | Evidence |
|---|---|---|---|---|
| — | **control** (unmutated donor) | all | n/a — **clean** | T1 pass (cov 0.50), T3 pass, T4 pass, T5 pass, T6 pass (882 distinct, p50 0.32845 — exactly the committed baseline), T7 pass (5 asserts) |
| D1 | Arithmetic that does not close: `0.21 * 0.17 = 0.0457` (true 0.0357) | T1 | **YES** | 599 failures / 200 seeds; `evaluated 0.0357 vs printed 0.0457 (tol 5e-05, delta 0.01)` |
| D2 | Printed step result off by 2e-4, display tolerance 5e-5 | T1 | **YES** | 600 failures / 200 seeds; `delta 0.0002 > tol 5e-05` |
| D3 | Unseeded RNG (`np.random.uniform` in a sampled value) | T3 | **YES** (on a real template) | See §3 F2. On the in-corpus instance of this defect, `template_levenspiel_plot_interpretation`, T3 reports **40/40 seeds differ between two processes**; the clean donor passes 0/40. The scratch mutant could not be run through T3 — that is F2. |
| D4a | Malformed step marker `**Step 2: **` | T4 | **YES** | `malformed markers {'**Step 2: ** Compute t': 200}; non-contiguous numbering on 200 seeds` |
| D4b | Restarted numbering 1,2,3,1,2 | T4 | **YES** | `duplicate step numbers on 200 seeds` |
| D5 | Return path yielding `("...", "...")` with no `**Step` and no `**Answer:**` (fires on ~half of seeds) | T4 | **YES** | `no answer marker on 108 seeds; degenerate output on 108 seeds` |
| D6 | Step result computed inline in the f-string (`= {sum(terms.values()):.4f}`) | T5a | **YES** | `[result-inline] template_mutant:109 sum(terms.values()) (computed inside the f-string in result position)` |
| D7 | Variable printed at 5 dp then consumed downstream unrounded | T5b | **YES** | `[unrounded-consumed] weekly_rate (printed at [5] dp, assigned unrounded, then consumed by a later computation; display rounding discards up to 4.29e-06 on 200 sampled instances)` — plus the derived `daily_pct` at 2 dp |
| D8 | Answer-space collapse (two sampled parameters pinned to constants) | T6 vs committed baseline | **YES** | `distinct answers fell 882 -> 7 (99% > 10%)` at the baseline's own 1000 seeds |
| D9 | Both invariant assert families stripped | T7 | **YES** | `0 asserts {}`, `pass=False` |
| D10 | Half-way-tie offset (`p2 + 0.00005`) applied through the generator | T1 (expected MARGINAL, per the known lead) | **YES, but not as predicted** | The offset propagated into the row-sum line and surfaced as 71 hard **FAIL**s (`delta 0.0001 > tol 5e-05`), not marginals. The tie behaviour is real but reproduces only at check level — see F3. It also tripped T5b on `p2`. |
| **D11** | **Gross non-closure on a prose-prefixed line: `Row sum: 0.4347 + 0.2952 + 0.2701 = 1.5000`** | **T1** | **NO — MISSED** | T1 pass, T4 pass, T5 pass, T7 pass. The line is classified `skipped_symbolic`. |

**Detected: 11 of 12. Missed: 1 (D11).**

---

## §3 Findings

### F1 — CONFIRMED — T1 silently skips every arithmetic line that has prose in front of it (blocker-grade)

`_strip_units()` in `t1_closure.py` poisons a segment when an alphabetic token sits in operand
position (`before` is empty or does not end in a digit / `)]%`). A sentence-leading word satisfies
that test, so the *whole line* is poisoned, `safe_eval` declines, and `check_line` returns
`(_, 'symbolic-or-unparseable')`. `check_instance` then files it under `skipped_symbolic` —
the bucket reserved for legitimate formula-only steps — so nothing in the reported output
distinguishes "this line stated a formula" from "this line did arithmetic and I could not read it".

Minimal reproduction:

```
python -c "import sys;sys.path.insert(0,'.');from tests.template_integrity.checks import t1_closure as t1;print(t1.check_line('Therefore the total is 12.5 * 4 = 99.0 kN'));print(t1.check_line('0.4347 + 0.2952 + 0.2701 = 9.9999'))"
```

The first (blatantly wrong) line yields `([], 'symbolic-or-unparseable')`; strip the four leading
words and the identical arithmetic is checked and fails. End-to-end reproduction is mutant D11
(`scratch/d11_prose_prefixed_line.py`, built by the snippet at the end of this section): a template
that prints a row-sum of `1.5000` and claims it equals 1 passes T1, T4, T5 and T7.

Why this matters for Phase 0's exit gate: the donor's T1 coverage is **0.50** — half of its `=`
lines are never evaluated, and the row-sum verification line, the one line whose whole purpose is
to be checkable, is among them. Coverage is computed (`ClosureResult.coverage`) and printed by
`run.py`, but **nothing gates on it**. A redesigned template could reach 0.0 coverage and still
report `T1 pass`.

Suggested fix shape: when a poisoned segment still contains a numeric expression after the leading
prose is dropped at the last `:`/sentence boundary, retry rather than skip; and separate
"unevaluable but numeric" from "genuinely symbolic" in the skip buckets.

### F2 — CONFIRMED — T3 cannot be pointed at any template outside `data/templates/branches/`, and crashes when you try

`t3_determinism._CHILD` resolves templates by calling `discover()` inside the child, ignoring the
`module` and `file_path` on the `TemplateRef` it was handed. For an unknown id the child writes
`{"__missing__": true}`, but the parent unconditionally does `int(k)` over that dict:

```
tests/template_integrity/checks/t3_determinism.py:94
    return {tid: {int(k): v for k, v in per.items()} for tid, per in raw.items()}
ValueError: invalid literal for int() with base 10: '__missing__'
```

Reproduction:

```
python -c "import sys;sys.path.insert(0,'.');from tests.template_integrity.core import TemplateRef;from tests.template_integrity.checks import t3_determinism;t3_determinism.run_many([TemplateRef('template_mutant','scratch.control','scratch/control.py','scratch','scratch',1)], range(4))"
```

Two consequences: (a) the sentinel path is dead code — a genuinely missing template raises instead
of being reported; (b) T3 cannot validate a candidate template before it is committed into the
branches tree, which is exactly when a redesign wants to check it. T3 itself is sound — pointed at
the real defect it reports `40/40 seeds differ between two processes` and passes the clean donor:

```
python -c "import sys;sys.path.insert(0,'.');from tests.template_integrity.core import discover;from tests.template_integrity.checks import t3_determinism;r={x.template_id:x for x in discover()};print({k:(v.passed,len(v.mismatched_seeds)) for k,v in t3_determinism.run_many([r['template_levenspiel_plot_interpretation'],r['template_two_step_transition_probability']],range(40)).items()})"
```

### F3 — CONFIRMED — the ±0.5-ulp rule cannot flag an exact half-way tie (the stated lead)

```
python -c "import sys;sys.path.insert(0,'.');from tests.template_integrity.checks import t1_closure as t1;c,_=t1.check_line('k = 0.01185 * 1000 = 11.8');e,v,tok,p=c[0];print(v,p,t1.display_tolerance(tok),abs(v-p))"
```

evaluated 11.85, printed 11.8, tol 0.05, delta 0.049999999999998934 → `delta <= tol` → `marginal=True`
→ `ClosureResult.failures` excludes it → `passed == True`. Confirmed exactly as described. Note the
delta is *below* tol only by floating-point luck; the same construction at a different scale would
land on the other side, so this band is not merely conservative, it is unstable.

### F4 — CONFIRMED — `run.py` gates on failures only, so unbounded marginal density reports T1 pass

`ClosureResult.passed` is `not self.failures`. `run.py` writes `'marginals': len(r.marginals)` into
the report and never consults it. Scanning the live corpus at 25 seeds:

```
python -c "import sys;sys.path.insert(0,'.');from tests.template_integrity.core import discover,generate;from tests.template_integrity.checks import t1_closure;rows=[];[rows.append((len(c.marginals),c.evaluated,r.template_id)) for r in discover() for c in [t1_closure.run([generate(r,s,capture=False) for s in range(25)],r.template_id)] if not c.failures and c.marginals];rows.sort(reverse=True);print(len(rows));print(rows[:5])"
```

**62 of 150 templates report `T1 pass` while carrying marginals.** Worst offenders (marginals /
evaluated lines, 25 seeds): `template_arl_beta_mean_shift` 128/75, `template_manning_trapezoidal_velocity`
115/175, `template_mm1k_finite_capacity` 87/150, `template_aoq_ati_rectifying` 81/172,
`template_primary_consolidation_settlement` 77/250. A template averaging more than one
last-digit-boundary discrepancy per evaluated line is not "marginally rounded"; combined with F3
this is a usable place to hide a defect.

### F5 — PLAUSIBLE — T5b's `assigned_rounded` is a whole-function AND, so one rounded assignment does not clear a name

`_Analyser.visit_Assign` folds `prev and rounded` across every assignment to a name. That is the
conservative direction (good), but it also means a name rebound in a loop as `x = round(x, 2)` after
an unrounded seed assignment stays flagged. Not exercised by a mutant here; flagged for whoever
tunes T5b's false-positive rate. Observed in passing: D7's derived `daily_pct` was reported
alongside the intended `weekly_rate` — correct, but it means one root cause yields N findings.

---

## §4 What I tried that did **not** break the harness

* **The unmutated control.** Byte-for-byte the donor's source, renamed. T1 pass (0 failures,
  0 marginals), T3 pass (0/40 mismatched), T4 pass, T5 pass, T6 pass against the committed
  baseline at its own 1000 seeds (882 distinct answers, p50 0.32845 — exact match), T7 pass.
  No check fires on the scaffolding itself, so §2's detections are real.
* **Sub-tolerance drift through the generator (D10).** I could not construct a mutant whose
  printed values stayed inside the marginal band end-to-end: at 4 dp display the `:.4f` format
  snaps a sub-tolerance offset straight back to the true value, and any offset large enough to
  survive formatting also exceeds tol on the downstream sum line. The F3 blind spot is real at
  check level but is harder to reach through a template than the lead implies.
* **T4 evasion.** `**Step 2: **` is caught by the loose-vs-strict marker diff, restarted numbering
  by the duplicate check, and an answer-less early return by the `empty_output` path even though it
  fires on only ~54% of seeds. I did not find a marker shape that slips through both regexes.
* **T6 evasion via a small collapse.** Pinning a *single* sampled parameter was not enough to breach
  the 10% distinct-answer floor; it took pinning two. That is the tolerance working as specified,
  not a hole, but it does mean a partial collapse is invisible.
* **T5a evasion.** `round(sum(...), 4)` in result position is still caught — `_is_bound` recurses
  through the `round` passthrough into its argument. Good.

---

## §5 Forward-looking suggestions, ranked

1. **Fix F1 before any Phase 1 template is certified.** Retry a poisoned segment after trimming the
   leading prose clause, and split `skipped_symbolic` into `symbolic` (no numerals after the `=`)
   vs `numeric-but-unevaluable`. Until then, T1's `pass` means "the half of the lines I could read
   closed", and the exit-gate claim in the rebuttal cannot be supported as written.
2. **Gate on T1 coverage, not just failures.** Add a floor (e.g. `coverage >= 0.8` and
   `numeric-but-unevaluable == 0`) to `run.py`'s pass condition. The donor sits at 0.50 today and
   reports green.
3. **Gate on marginal density (F3 + F4).** Fail when marginals per evaluated line exceeds a small
   threshold, and treat `delta` within 1 ulp of `tol` as a tie that fails rather than passes.
   62/150 templates are currently green under this gap.
4. **Make T3 usable pre-commit (F2).** Have the child reconstruct refs from the parent's
   `(template_id, module, file_path)` rather than re-running `discover()`, and make the
   `__missing__` sentinel survive the parent's `int(k)` parse.
5. **Add a permanent mutation-test target.** `scratch/build.py` + `scratch/detect.py` are throwaway,
   but the 12 mutants are cheap and fast (~90 s at 200 seeds). Promoting them to
   `tests/template_integrity/test_mutants.py` turns "the harness can detect a planted defect" from
   a one-off review claim into a CI invariant, and would have caught F1 the day it was written.
