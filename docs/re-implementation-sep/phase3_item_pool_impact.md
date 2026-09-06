# Phase 3 — item-pool impact and D3.5 distribution diff

**Phase:** 3 (trace shape: iteration and search) · **Date:** 2026-09-06
**Before:** `master` at `00bd0b0`, checked out as a worktree at `C:/wtm`
**After:** `redesign/phase3-trace-shape` at `953adc9`
**Templates:** `template_normal_depth_iteration` (civil, uniform flow) ·
`template_line_balancing_heuristic` (industrial, production planning)

**Bottom line.** Regenerating these two templates changes **2.30%** and **0.65%**
of instances respectively. Nothing else in the corpus moves. The answer *space*
is unchanged on both — same distinct-answer count, same quantiles to three
decimals bar one bin at the p95 of the civil item — so any published aggregate
over these two templates shifts by less than the instance churn, and no result
outside them is affected at all.

---

## 1. How this was measured, and why not with T6

T6 fails **142/150 on `master`**: the committed baseline profile is stale
corpus-wide, so T6 cannot gate anything this phase touches. That was re-measured
here rather than taken from the brief — a full `--checks all` run in a worktree at
`00bd0b0` gives T1 30, T2 0, T3 0, T4 3, T5 68, T6 142, T7 83, reproducing the
brief's table exactly.

Regenerating the baseline would turn T6 green and is the wrong move: **a baseline
refreshed by the phase it is meant to gate is not a gate**, and it would absorb
corpus-wide movement that has nothing to do with Phase 3. Phase 2 hit this and
delivered the diff directly; Phase 3 does the same. Recorded as **D-043**;
regenerating the baseline as its own deliberate change on `master` is left to
Phase 6.

**Instead, a direct before/after instance dump**, 4,000 seeds per template per
tree, **each tree in its own process**:

```
git worktree add C:/wtm master
python -m tests.template_integrity.phase3_instance_dump C:/wtm before.json 4000
python -m tests.template_integrity.phase3_instance_dump .     after.json  4000
```

An in-process reload resolves both sides to the same already-imported
`data.templates.*` modules and reports every instance identical — the trap
documented in `instance_dump.py`, which has now caught three people. The Phase 3
dumper additionally **refuses to run** if a template resolves to a file outside
the tree root it was handed, so the mistake fails loudly instead of producing a
clean-looking null result.

**Resolution limit.** 4,000 seeds resolve a per-instance rate down to ~0.075%
(3/N). Both churn rates measured here (2.30%, 0.65%) are one to two orders of
magnitude above that. The tie censuses behind them were run at 20,000 seeds,
resolving to ~0.015%.

---

## 2. `template_normal_depth_iteration`

### 2.1 What changed and why

| Change | Instances affected |
|---|---|
| Display-tie resample screen on Step 1's `K` line (3 dp) and on the secant update lines (4 dp) — **D-016**, D-040 | **2.30%** redraw to a different instance |
| Step 3 header gained a sentence saying the update count is not fixed in advance and stating it | 100% (prose only) |
| Step 4 now **prints** the last depth change instead of asserting it is small | 100% (prose only) |
| `fmt3()` removed; the stored `g` is normalised instead of the printed string | 0% — never fired in 4,000 seeds (D-042) |
| `g_curr - g_prev` guarded; non-convergence raises instead of asserting | 0% — never fired in 20,000 seeds (D-041, D-042) |

### 2.2 Distribution diff (D3.5), 4,000 seeds

| Property | Before | After | Movement | T6 gate |
|---|---|---|---|---|
| Distinct answers | 267 | 267 | **0.00%** | ≥ −10% ✔ |
| Answer p5 (m) | 0.8600 | 0.8600 | +0.000% | ±2% ✔ |
| Answer p50 (m) | 1.4000 | 1.4000 | +0.000% | ±2% ✔ |
| Answer p95 (m) | 1.9300 | 1.9200 | −0.518% | ±2% ✔ |
| Step counts | `{4: 4000}` | `{4: 4000}` | none | ✔ |
| Update counts | `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}` | `{1: 32, 2: 314, 3: 2342, 4: 1261, 5: 51}` | ≤ 3 instances per bin | ✔ |
| Shape branch (rect/trap) | 50.02 / 49.98 | 50.12 / 49.88 | 0.10 pts | ±5 pts ✔ |
| Lining mix (5 linings) | — | — | ≤ 0.4 pts each | ±5 pts ✔ |

The p95 move is one display bin — 1.93 m to 1.92 m — on a quantity quoted to
3 dp, i.e. the 3,800th-ranked instance of 4,000 shifted by 0.01 m. It is churn,
not a distribution change.

### 2.3 Instances whose *answer* changed

**92 of 4,000 (2.30%), and every one of them is an instance whose question also
changed.** Zero instances kept their question and changed their answer. That is
the important number: no solver who saw a `master` instance would now be graded
against a different answer for the same problem. The 92 are resampled away
entirely and replaced by fresh draws.

### 2.4 Are the rejected instances a systematic slice?

The screen rejects instances where an exact printed expression lands on a
half-way display boundary. Those are arithmetic accidents of `Q·n/√S` and of the
secant quotient, not a region of the physical parameter space, so the rejection
should be scattered. It is — measured directly on the 92, not inferred from
aggregates:

| Parameter | Rejected 92, min / median / max | Full pool, min / median / max |
|---|---|---|
| `Q` (m³/s) | 3.12 / 11.93 / 56.92 | 2.40 / 14.62 / 58.14 |
| `b` (m) | 2.00 / 3.30 / 5.00 | 2.00 / 3.50 / 5.00 |
| `S` | 0.0008 / 0.0016 / 0.0025 | 0.0008 / 0.0012 / 0.0030 |

The rejected set spans essentially the whole sampled box in every parameter. Its
shape split is **48 trapezoidal / 44 rectangular** against a 50.0/50.0 pool, and
its update-count profile (`{2: 3, 3: 62, 4: 24, 5: 3}`) tracks the pool's
proportions — so the screen is not preferentially removing slowly-converging
channels, which was the failure mode to worry about. Corroborating aggregates:
distinct answers 267/267, p5 and p50 identical to four decimals, step count
constant at 4, branch mix 50.02/49.98 → 50.12/49.88 (0.10 pts) and lining mix
stable within 0.4 pts on all five linings.

### 2.5 Non-resampled instances are byte-identical in their arithmetic

Of the 3,908 seeds the screen did not touch, **every single differing line is one
of the two deliberately added prose sentences** — classified mechanically:

```
3908  Step 3 header
3908  Step 4 body
   0  anything else
```

No arithmetic line, no operand, no answer changed on any of them. This is the
evidence that rebuilding the template around a structured trace node
(D-039) was behaviour-preserving rather than a rewrite.

---

## 3. `template_line_balancing_heuristic`

### 3.1 What changed and why

| Change | Instances affected |
|---|---|
| Display-tie screen on the exact rational balance delay (1 dp) — **D-016** | **0.65%** redraw |
| `n*CT - total` and the duration sum bound out of the f-string (T5 result-inline; Phase 2 Reviewer C's C-3 pattern) | 0% — integer arithmetic, no float artefact was reachable |
| `_t19_assign` returns the `decision` trace node; prose rendered from it | 0% |
| Tie-break rule stated explicitly in the rule string | 0% — never exercised; the five durations are distinct by construction |

### 3.2 Distribution diff (D3.5), 4,000 seeds

| Property | Before | After | Movement | T6 gate |
|---|---|---|---|---|
| Distinct answers | 335 | 335 | **0.00%** | ≥ −10% ✔ |
| Answer p5 (%) | 13.5000 | 13.5000 | +0.000% | ±2% ✔ |
| Answer p50 (%) | 27.1000 | 27.1000 | +0.000% | ±2% ✔ |
| Answer p95 (%) | 37.1000 | 37.1000 | +0.000% | ±2% ✔ |
| Step counts | `{5: 2116, 6: 1884}` | `{5: 2125, 6: 1875}` | 0.22 pts | ±5 pts ✔ |
| Station count `n` | 3: 52.90%, 4: 47.10% | 3: 53.12%, 4: 46.88% | 0.22 pts | ±5 pts ✔ |

### 3.3 Instances whose *answer* changed

**25 of 4,000 (0.62%)**, all of them among the 26 resampled instances. As with
the civil item, **zero instances kept their question and changed their answer.**

### 3.4 The industrial screen IS skewed, and the skew is stated rather than smoothed

Unlike the civil item, this screen does **not** reject uniformly. Of the 26
rejected instances, **22 have `n = 4` and 4 have `n = 3`** — so it removes 1.17%
of the four-station population against 0.19% of the three-station one, a **6×
difference**.

This is not a coincidence and it is not a defect. A 1-dp tie requires the exact
rational `100·idle/(n·CT)` to have a terminating expansion ending in a 5 at the
second decimal, which needs `n·CT` to carry the right powers of 2 and 5. At
`n = 4` the factor `4` supplies two extra powers of two, so ties are
structurally commoner there. The cycle-time range of the rejected set
(60 / 112 / 144 against a pool of 59 / 126 / 150) shows the same effect from the
`CT` side.

**Net effect on the item pool: the station split moves 47.10% → 46.88%, a
0.22-point shift**, against T6's ±5-point branch tolerance — two orders of
magnitude inside the gate. So the skew is real, understood, and immaterial. It
is recorded because "the screen rejects at random" would have been the
comfortable claim and it is not the true one; a future screen on a template
where the answer-bearing branch is finer-grained could be moved by the same
mechanism, and should measure this rather than assume it.

### 3.5 Prose is otherwise byte-identical

Of the 3,974 seeds the screen did not touch, **0 have any changed line at all** —
not one character. Restructuring `_t19_assign` to return a structured decision
node and rendering the six step blocks from it reproduces the previous output
exactly. That is a stronger result than the civil item's, where two prose lines
were changed on purpose.

---

## 4. Which published results are invalidated

**Scope of the damage: two templates of 150, and within them only the resampled
instances.**

| Question | Answer |
|---|---|
| Are any *other* templates' instances affected? | **No.** A full corpus `--checks all` run before and after shows **3 check-status changes, all improvements, all on these two templates** (`normal_depth` T1 and T5, `line_balancing` T5). 147 templates unchanged on every check. |
| Are the published *answers* for these two templates wrong? | For **2.30%** and **0.62%** of instances, the item a given seed produces is now a different item. The old answers were correct answers to the old questions — except on the tie instances, where by D-016 the old question had **no defensible answer at all**, which is why they were removed. |
| How many published instances were ill-posed? | Of the churn, essentially all of it. Over 20,000 seeds the civil item carried **427 exact ties (2.14%)** and the industrial item **~0.65%**. Only 27 and 204 of those surfaced as T1 hard failures — see D-040 on why counting T1 failures understates the tie population by ~10×. |
| Does any aggregate metric move? | Not measurably. Distinct answers, p5, p50 identical; p95 moves one display bin on one template. Any accuracy figure computed over a fresh sample of these templates should move by less than the churn rate. |
| Do the *model generations* in `inference_results/` need regenerating? | For these two templates, the affected seeds only. This depends on **D-003** (do the raw generations still exist?), which remains open and is not settled here. |

**What is NOT claimed.** This note does not say the published numbers were wrong
by 2.30%. It says 2.30% of *instances* were replaced. Whether that moves a
reported accuracy depends on how the evaluation sampled seeds, which is the
`evaluation/` track and out of scope — its parser has known defects (D-003) and
was deliberately not touched.

---

## 5. Reproduction

```bash
git worktree add C:/wtm master
export PYTHONIOENCODING=utf-8

# corpus non-regression, one run per tree
(cd C:/wtm && python -m tests.template_integrity.run --checks all --json before.json)
python -m tests.template_integrity.run --checks all --json after.json

# per-template gate, 200 seeds
python -m tests.template_integrity.run \
  --templates template_normal_depth_iteration,template_line_balancing_heuristic \
  --checks all --seeds 200

# the instance dump behind §2 and §3 — separate process per tree, mandatory
python -m tests.template_integrity.phase3_instance_dump C:/wtm before.json 4000
python -m tests.template_integrity.phase3_instance_dump .     after.json  4000

# the conformance corpus behind D3.3
python -m tests.trace_schema.extract \
  docs/re-implementation-sep/phase3_conformance/traces.json 40
```
