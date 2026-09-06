# Phase 3 — Reviewer B: Physics & Pedagogy (the P6 guard)

**Ref under review:** `953adc9` on `redesign/phase3-trace-shape`
**Baseline:** `master` @ `00bd0b0`, worktree `C:/wtm`
**Scope:** `template_normal_depth_iteration` (T24, civil), `template_line_balancing_heuristic` (T19, industrial)
**Date:** 2026-09-06

---

## 1. Verdict

**PASS WITH FINDINGS** — neither item changed what it tests or how hard it is, and the gate question resolves in the item's favour (`n` is NOT recoverable without running the rule); but two of the new resample screens are *clustered*, not scattered, and that was not recorded.

---

## 2. Independent re-derivation

All numbers below are mine, generated in a fresh process per tree (never an in-process reload). Probe scripts are in the session scratchpad; each is self-contained and replays the template's own accept loop.

**Ref-drift check (done first, because it decides whether any of this counts).** The branch moved under me: the working tree I generated from does not hash-match `953adc9` for either template. I diffed both before trusting a single number —

```
git diff 953adc9 -- data/templates/
```

— and both drifts are **comment-only** (expanded D-037 justifications in the two screen comments; no change to any executable line, sampled parameter, or emitted string). Every measurement below therefore holds for the frozen ref. Flagging it because a reviewer who had not checked would have been silently reviewing something else.

### 2.1 The gate question: is `line_balancing_heuristic` a lookup?

The item's answer is a strict function of `n`. So the item is a lookup iff `n` is obtainable from the question text without executing the station-by-station assignment. I replayed the accept loop's sampling, captured `n` and `N_min = ceil(sum(t)/CT)` for every accepted instance, and scored every cheap predictor of `n` I could construct from question-text quantities alone (durations, precedence, `CT`, `sum(t)`, `N_min`).

Command: `python scratchpad/lb_probe.py` and `python scratchpad/lb_probe2.py` (4,000 and 8,000 seeds, head tree).

| Predictor of `n`, using question text only | Accuracy (8,000 seeds) |
|---|---|
| `n = N_min` (the trivially computable bound) | **58.33 %** |
| `n = N_min + 1` | 41.67 % |
| `n = 3` (best blind constant) | 53.05 % |
| `n = max(N_min, 3)` | **64.97 %** — best rule found |
| `n = max(N_min, #tasks with 2t > CT)` | 64.64 % |
| best `n = N_min + [single cheap predicate]` (`t_d + t_e > CT`) | 61.20 % |
| per-`N_min` majority vote (the ceiling of `N_min` alone) | 65.55 % |

Joint distribution over 4,000 seeds: `(N_min, n) = (2,3): 277, (3,3): 1848, (3,4): 1378, (4,4): 497`.

**Answer to the gate question, with the number that settles it: `n` equals `N_min` in 58.6 % of instances (2,345 / 4,000; 58.3 % over 8,000 seeds). The best question-text-only rule I could construct reaches 65.0 %, against a 53.1 % blind-constant floor. The greedy rule must therefore actually be executed on at least 35 % of instances, and the shortcut buys only ~12 points over guessing. The item is a search, not a lookup.**

The structure behind that number is why the item survives: `N_min` is fully determining only in its tails — `N_min = 2 → n = 3` (277/277) and `N_min = 4 → n = 4` (497/497). The mass sits at `N_min = 3`, where the split is 1848 / 1378, i.e. a 57/43 coin. 81 % of the pool lands in that ambiguous cell.

**No regression from master.** Same measurement on `C:/wtm` (`python scratchpad/lb_master.py`, 8,000 seeds): `n == N_min` 58.27 %, `max(N_min,3)` 64.84 %, `n` mix `{3: 4222, 4: 3778}`. Head over the same seed span: 58.33 %, 64.97 %, `{3: 4235, 4: 3765}`. Divergence is inside noise. **Phase 3 did not make this item more guessable.**

### 2.2 The civil item's "incidental iteration count"

Command: `python scratchpad/nd_probe.py` (60,000 candidate draws, head sampling replayed with rejection reasons captured).

Update-count distribution over 58,676 accepted: `{1: 503, 2: 4438, 3: 34749, 4: 18097, 5: 889}` — rescaled to 4,000, `{1: 34, 2: 303, 3: 2369, 4: 1234, 5: 61}`. The docstring claims `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}`. **Agreement.** `> 3` updates: I get **32.4 %**, docstring claims 32.9 %. **Agreement.** The count is genuinely data-dependent, so calling it incidental is defensible: the answer `yn` is the converged fixed point, and two solvers reaching it in 3 and 4 updates print different traces and the same answer.

Rejection rates: `K_tie` **1.197 %**, `upd_tie` **1.010 %**, total **2.21 %** over 60,000 draws. Docstring claims 1.14 % / 1.02 % / 2.14 % over 20,000 seeds. **Agreement.**

### 2.3 The "only two added sentences" claim

Command: `python scratchpad/dump.py <tree> <out>` — 600 seeds per template per tree, then a line-level diff classifier that normalises all numerals and clusters the resulting diff signatures. This is the decisive test of the prose claim.

- **T24:** 600/600 instances differ, but **586/600 (97.7 %) differ by exactly the two claimed sentences and nothing else** — the Step 3 clause ("The scheme runs until the tolerance is met, so the number of updates is not fixed in advance; here it takes N.") and the Step 4 rewrite ("The last change is N m, below the tolerance of N m, so yn = ..."). The remaining **14/600 (2.33 %)** are whole different instances — different shape, lining, `b`, `z`, `S` — because a screen rejection shifted the RNG stream. 2.33 % matches the measured 2.21 % rejection rate. **The prose claim holds.**
- **T19:** **594/600 (99.0 %) byte-identical to master.** The 6 that differ are whole-instance resamples from the new tie screen. In particular the tie-break change `max(fitting, key=times[x])` → `max(fitting, key=(times[x], -index))` is behaviourally inert, as claimed, because the five durations are distinct by construction.

### 2.4 The C-6 non-convergence raise

Command: `python scratchpad/nd_master.py`. Over 6,000 master seeds, the number of instances where Step 4 asserts "the last change is below the tolerance" while the last printed update's change was **not** below 0.002 m: **0**. The `RuntimeError` is a correct guard but fixes nothing observed. No pedagogy impact either way; recorded so nobody credits it as a live-defect fix. The `ZeroDivisionError` guard likewise never fires (0 in 60,000 draws).

---

## 3. Findings

### F-B1 — `CONFIRMED` — The civil `K` display-tie screen is single-valued in slope, not scattered

Every one of the 718 `K_tie` rejections in 60,000 draws has **`S = 0.0016` exactly** (min = median = max = 0.00160). This is not chance. `_froude_capped_slope` returns `round(..., 4)`, so `S` lives on a 4-dp grid, and `S = 0.0016` is the only grid value in the `[0.0008, 0.003]` window whose square root is exactly representable at 5 dp (`sqrt(0.0016) = 0.04`). `K = Q*n / sqS` is then `25*Q*n`, a terminating decimal, so it is the only slope at which a 3-dp half-way tie is reachable at all. The screen removes instances from **one slope value and no other**.

Consequence, measured: `S = 0.0016` falls from **4.83 %** of master instances to **3.50 %** of head instances (600 seeds each) — a ~28 % relative depletion of that slope. The lining distribution inside the rejected slice is also skewed (`n = 0.012` never appears among the 718).

Reproduction:
```
PYTHONIOENCODING=utf-8 python scratchpad/nd_probe.py
    # -> "K_tie  N=718 ... S[min,med,max]=0.00160/0.00160/0.00160"
PYTHONIOENCODING=utf-8 python scratchpad/dump.py 'C:\wtm' master.txt
PYTHONIOENCODING=utf-8 python scratchpad/dump.py 'C:\Users\ayesha.gull01\EngTrace' head.txt
    # then grep the "slope of S = " lines: 4.833% vs 3.500%
```

**Impact: LOW on pedagogy, but the record is wrong.** `S = 0.0016` is not an engineering-distinguished slope — it is not a steep/mild boundary, not a Froude threshold, and the item's physics is invariant across the window — so no *class* of channel disappears and the Sturm Ch. 4 content is untouched. But the docstring's framing ("an ill-posed instance rather than a rounding-convention choice") invites the reader to picture scattered removals, and the brief's own test asks whether the rejections cluster. They cluster completely. **This should be recorded as a one-line scoping note in the docstring and in `phase3_item_pool_impact.md`: the screen is reachable only at `S = 0.0016` and depletes it by ~28 %.** It does not block the merge — P6 governs what the item *tests*, and nothing about what the item tests changed.

### F-B2 — `CONFIRMED` — The line-balancing tie screen is likewise clustered, and skews toward `n = 4`

Of 913 tie rejections in a 200,000-draw sweep, **913/913 have `n*CT` divisible by 8** (every observed capacity is in fact a multiple of 16), and **781/913 (85.5 %) have `n = 4`**, against 46.6 % in the accepted pool. The rejected balance delays are all of the form `x.x5` by construction.

Reproduction: `PYTHONIOENCODING=utf-8 python scratchpad/lb_screen.py`

**Impact: NEGLIGIBLE, and quantified rather than asserted.** Folding the rejected instances back in moves the `n = 4` share from 46.60 % to 46.83 % — a **0.23 pp** shift. The overall rejection rate I measure is **0.592 %** of post-`n`/post-delay candidates (docstring claims 0.65 %; agreement). No pedagogy consequence; recorded for the same "clustered, not scattered" reason as F-B1.

### F-B3 — `PLAUSIBLE` — The new Step 4 sentence surfaces a pre-existing precision mismatch

Step 4 now prints the last change to 4 dp against a 0.002 m tolerance. A reader who takes that seriously will notice the item converges to ±2 mm but reports `yn` to 3 dp (±0.5 mm). I checked whether it matters: bisecting the exact root of `A*R^(2/3) = K` and rounding to 3 dp agrees with the traced `yn` in **96.92 %** of accepted instances; in the other **3.08 %** they differ, always by exactly 0.001 m (max `|yn - direct| = 0.0010 m`). So a solver who solves Manning's equation *correctly and directly*, rather than by the prescribed secant, gets a different gold answer about once in 32 instances.

Reproduction: `PYTHONIOENCODING=utf-8 python scratchpad/nd_probe.py` — the `ANSWER-BEARING PROBE` line.

**Impact: LOW, and NOT a Phase-3 regression** — the tolerance, the 4-dp secant grid and the 3-dp answer are all unchanged from master, so the 3.08 % existed before. What Phase 3 changed is that the trace now shows the reader the arithmetic that makes the mismatch legible. Two defensible responses: state in the question that the answer is the depth *the scheme returns* (making the method answer-bearing on purpose), or tighten the tolerance to 0.0005 m so the secant result and the true root agree at 3 dp. The first is cheaper and is what the item already means.

### F-B4 — `PLAUSIBLE` — Step 3 now discloses the update count before the trace

"...here it takes 3." precedes the iteration lines. Because the count is incidental (§2.2) this is not an answer leak, and it makes the pedagogical point the sentence exists to make — that the stopping rule is a tolerance, not a recipe. But it does hand the reader the trace length in advance, removing a small amount of "have I finished?" work from the item.

Reproduction: any instance; see the §2.3 dump, signature `[586 instances]`.

**Impact: VERY LOW.** A domain expert would read this as clarification, not as an easier item. I flag it only because the brief asks whether anything got easier and this is the single change that points that way. It does not block.

### F-B5 — No finding on difficulty or engineering content (stated for the record)

A hydraulics reader gets the same item: Manning `K = Qn/S^(1/2)`, section factor `A·R^(2/3)` monotone in depth, secant bracketing from 1.000/1.500 m, subcritical check `Fr <= 0.92` — all intact and all consistent with Sturm Ch. 4. An IE reader gets the same ranked-heuristic line-balancing item in Nahmias §9.10 form: precedence-feasible eligible set, longest-fitting pick, station closure, `n` against `N_min`, balance delay. Neither item's Advanced difficulty label is now wrong.

---

## 4. Falsification attempts that failed

Things I tried to break and could not:

1. **"`n` is a lookup."** Attacked six ways: `N_min` directly; `N_min + 1`; best blind constant; `max(N_min, 3)`; a must-be-alone bin-packing lower bound `#{t : 2t > CT}`; and an exhaustive search over all `n = N_min + [pairwise-sum-exceeds-CT]` predicates (all 10 task pairs plus 5 aggregate variants). **Nothing exceeded 65.0 %.** The per-`N_min` majority ceiling — the best possible rule using `N_min` alone — is 65.55 %, so no cleverer `N_min`-based rule exists. The claim that the station count is answer-bearing and irreducible survives.
2. **"Fixing `n` to a constant would not give the answer away."** `n` is genuinely bimodal (3: 53.1 %, 4: 46.9 %); fixing it would collapse the balance delay to a one-step arithmetic substitution. The implementer's justification for keeping `n` variable holds.
3. **"The iteration count is incidental."** Looked for a way for the count to leak into the answer and found none: `yn` is the converged value, the loop's carry is `y_prev <- y_curr, y_curr <- y_next`, and the termination predicate is a change threshold, not a pass counter. Reproduced the claimed distribution and the 32.9 % figure independently.
4. **"Only two sentences changed."** Tried to find a third prose change across 600 seeds by normalising all numerals and clustering the diff signatures. 586/600 diffs collapse to exactly the two claimed edits; the other 14 are full resamples, not prose. **The claim is true.**
5. **The tie-break rewrite.** Tried to find an instance where `key=(times[x], -_T19_ORDER.index(x))` picks differently from `key=times[x]`. 594/600 byte-identical outputs; 0 duplicate-duration instances, as claimed.
6. **The two new raises.** Tried to trip `ZeroDivisionError` (secant denominator) and `RuntimeError` (non-convergence) across 60,000 civil draws. Neither fires. Also confirmed master never actually emitted the mis-claiming trace the `RuntimeError` guards against (0 / 6,000).
7. **The `fmt3` removal.** Checked that normalising the stored `g` to `+0.0` prints identically to the old string patch — `f"{0.0:.3f}"` is `"0.000"` either way, and the printed operand now provably equals the stored one.

---

## 5. Further probing and improvements

**Where the evidence is thinnest.** The slope-depletion figure in F-B1 (4.83 % → 3.50 %) rests on 600 seeds per tree; the *structural* claim (all ties at `S = 0.0016`) rests on 60,000 draws and is solid, but the depletion magnitude deserves a 20,000-seed confirmation. Second thinnest: I scored shortcut rules I could think of, not a learned upper bound — a small decision tree over `(sorted durations, CT)` might beat 65 %, and that is the number I would most want checked.

**Ranked, with more time:**

1. Fit a depth-3 decision tree, or exhaustively search 2-predicate conjunctions over the full feature set, to establish a true upper bound on question-text-only prediction of `n`. If it clears ~85 %, F-B1's verdict flips and the item needs a wider `n` range.
2. Run a frontier model on ~200 instances and measure how often it reports an assignment consistent with its own stated `n`. My measurement bounds what a shortcut *can* achieve; it does not show what models *do*.
3. Sweep `S` on a 5-dp grid to confirm `S = 0.0016` is the sole tie-reachable value rather than merely the sole one in this window.
4. Check whether the 3.08 % secant/true-root disagreement (F-B3) concentrates in the 1- and 2-update instances, which would make it a stopping-too-early artefact rather than a rounding one.

**Does any weakness found here likely exist outside this phase's scope?** Yes, and it should be swept:

- **The `S = 0.0016` exactness pathology is a property of `_froude_capped_slope`, not of T24.** Every template in `data/templates/branches/civil_engineering/water_resources/uniform_flow.py` that calls it and then divides by `round(math.sqrt(S), 5)` has the same terminating-decimal special case — `template_manning_rectangular_discharge` (around line 109) does exactly that. If display-tie screening is rolled out corpus-wide in a later phase it will silently deplete `S = 0.0016` across the whole civil/water-resources branch, not just here. Worth one shared note rather than N rediscoveries.
- **The "is the count answer-bearing?" question generalises.** Any template whose trace length is data-dependent needs the T24 treatment; any whose trace length is answer-bearing needs the T19 treatment. The corpus has more iteration-shaped items than these two, and the classification should be made once, explicitly, rather than rediscovered per template.

**What later phases should do differently.** Make "does this screen cluster?" a *required recorded measurement*, not a reviewer's question. Both screens here cluster completely; both were described as removing ill-posed instances; both descriptions are true but incomplete. A one-line "rejected slice profile" in the phase's item-pool-impact doc — the marginal of every sampled parameter over the rejected set — would have surfaced F-B1 and F-B2 without a reviewer.

**What in the spec was wrong or ambiguous.** The brief frames the gate as "has it become a lookup", a binary. The honest answer is a rate, and rates need a floor to be read against: 58.6 % sounds alarming until you see that the blind-constant floor is 53.1 %, at which point the lift is 12 points, not 59. I would restate the gate as "what does the best question-text-only shortcut buy over the best blind guess?", with an explicit threshold (say, a shortcut must not exceed 80 %). Second ambiguity: "difficulty unchanged" was never operationalised, so F-B4 has no test it could fail — I could only reason about it. A step-count or required-inference-count proxy would give later reviewers something to measure.

---

*Filed by Reviewer B. Not committed — left in the working tree for the implementer.*
