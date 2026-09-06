# Phase 0 — T2 oracle findings

Consolidated from three independent oracle-writing passes over the nine original
Phase 1 templates. Each pass worked **from the docstring's governing equations
and the generated text, with the template's computation code withheld** — an
oracle that reuses the template's own arithmetic agrees with it whether it is
right or wrong, which is the single most likely way this check gets silently
disarmed (spec Risks).

Recorded here because these judgements existed only in agent transcripts and in
the oracle modules' docstrings. The *conclusions* are what later phases need.

Oracles live in `tests/template_integrity/oracles/`. Reproduce with:

```
python -m tests.template_integrity.run --checks T2 --seeds 500
```

---

## 1. Results

All nine oracles: `unparsed = 0`, and each was mutation-tested (injecting a 1–2%
error or a flipped label into the gold makes the applicable instances fail), so
none is vacuous.

| Template | Fail rate | Worst rel. error | Tolerance | Verdict |
|---|---:|---:|---:|---|
| `mean_variance` | **56.8%** | 100%* | 1e-3 | **chain break** |
| `vibration_transmissibility` | 4.8% | 1.62% | 5e-3 | **chain break** |
| `rotating_unbalance` | 4.6% | 8.43% | 5e-3 | **chain break** |
| `annulus_flowrate` | 0% | 4.4e-5 | 1e-4 | answer sound; Step-5 line defective |
| `cantilever_double_integration` | 0% | 0 | 1e-3 | answer sound; Step-6 line defective |
| `damping_classification` | 0% | 2.5e-5 | 5e-3 | ill-posed, not wrong |
| `system_properties` | 0% | 2.5e-5 | 2e-3 | not defective |
| `poissons_ratio` | 0% | 0 | 5e-3 | not defective (presentation issue only) |
| `logarithmic_decrement` | 0% | 9.3e-4 | 5e-3 | intended behaviour |

\* relative error saturates where the true mean is near zero (seed 194: question
implies −0.002, trace answers −0.000).

---

## 2. Confirmed chain breaks

**`mean_variance`** — the printed answer is not what the printed probabilities
give, at the 3 dp it is quoted to, on **244/300 seeds (mean)** and **261/300
(variance)**. Mechanism visible in Step 1: at seed 0 the question states
p = 0.172 while the trace prints the product −1.552, i.e. −9 × 0.17244…
The rate is tolerance-sensitive because the defect is fine drift — over 300
seeds the mean alone fails 187 at 5e-4, 120 at 1e-3, 40 at 2e-3, 15 at 5e-3.
**Second finding: the probabilities are never renormalised** — they sum to
0.999–1.002.

**`rotating_unbalance`** — root cause identified independently: the question
states the shaft speed rounded to whole RPM while the solution computes with the
unrounded speed (seed 230 states "53 RPM" against an internal 53.435). Every
other given reproduces the printed intermediates exactly, so this is the single
leak. `F_0 ∝ omega²` and a near-resonant denominator amplify it.
**Stricter measure: 101/300 (33.7%) of answers are not the correct 3-decimal
rounding of what the question implies.**
**Separate second defect:** some answers print to one significant figure
(seed 15: 0.001844 mm printed as "0.002 mm") — ungradeable regardless.

**`vibration_transmissibility`** — the defect is **more prevalent and smaller**
than the audit implies. The solution carries the frequency ratio rounded to 3 dp
into `(1−r²)²` (seed 3 uses r = 3.099 where the givens imply 3.09930), and the
base frequency is stated to 2 dp against an internal 41.0056 Hz. That makes
**117/300 (39%)** of answers not the correct rounding of the derived value — but
its worst genuine relative magnitude is **~0.48%**, not 0.87%. The larger
reported errors are display resolution on small TR values (TR ≈ 0.056–0.082 at
3 dp), not reasoning.

---

## 3. Display defects — the answer is sound, an intermediate line is not

**`cantilever_double_integration`** — 11/300 seeds print a Step-6 line whose own
arithmetic does not reproduce (seed 48: `delta = 0.01685 * 1000 = 16.8 mm`,
where 0.01685×1000 displays as 16.9). The **answer** is rounded from the
full-precision metre value and is correct. Sample rate 3.7%, against the audit's
~4.95% and D0.5's 5.20%.

**`annulus_flowrate`** — 283/300 seeds print a Step-5 product that does not
close (shape factor displayed to 4 dp but multiplied at full precision), worst
0.14%. But the printed Q is the full-precision result of the values stated in
the question, and the stated pressure drop *is* the value used, so **no hidden
back-solved Reynolds precision leaks into the answer** — contradicting the
audit's framing of this as a back-solving defect.

---

## 4. Not defective

**`damping_classification`** — see [`../DECISIONS.md`](../DECISIONS.md) D-010.
Both agents that examined it reproduced the ~33.7% flip rate under a strict
`zeta == 1` test, then showed every flip lies within **1.0× the half-unit
uncertainty** of the 2-dp damping coefficient (max ratio 0.99). No ζ anywhere in
(0.95, 1.05) other than the deliberately-critical instances. The finding is
ill-posedness: a third of instances sit on a measure-zero boundary and the 2-dp
rounding of `c` destroys the digit that would settle it.

**`poissons_ratio`** — forward-computing from the stated rounded values
reproduces both printed answers to the last printed digit on **1000/1000 seeds,
in both unit systems** (147 US-customary, 153 SI over seeds 0–299; 141
compression instances). No back-solving. The sign leak is real but is
presentation: the question states |P| with "(compression)" in prose while Step 1
prints `sigma = P / A = (-7000 N) / …`. It affects US instances identically, so
the audit's "76/300 SI" undercounts. *Suggested wording:* "an axial compressive
load of 7 kN (P = −7 kN)".

**`logarithmic_decrement`** — the trace computes δ and ζ from the amplitudes **as
stated**; agreement is 9.3e-4 relative, pure display rounding. The hidden-ζ
back-solving is confined to generating the amplitudes and stating them at 1 dp,
which is the measurement-uncertainty framing of the exercise. For the record, an
oracle trying to recover the *hidden* ζ would need `TOLERANCE ≈ 0.12` (max 11.0%,
p99 8.0%) — a tolerance so wide it would test nothing. That is why the oracle
asserts against the stated amplitudes.

**`system_properties`** — reproduced the audit's figure exactly: max
|ζ(from stated rounded `c`) − sampled `target_zeta`| = **4.500e-06 over 300
seeds**. But it never reaches the printed answer: the template prints ζ computed
from the stated `c`, so the residue exists only against the hidden sampled ζ and
is two orders below the 4-dp print precision. *Coverage gap:* no
critically-damped instances occurred in 400 seeds, so that branch is untested by
those seeds.

---

## 5. Findings about the question text

These made `parse_givens` fragile. Each is a property of the template worth
fixing or at least knowing about, because **anything that makes an oracle's
parser fragile makes a model's reading fragile too.**

| Template | Issue |
|---|---|
| `cantilever_double_integration` | `I` is stated as split scientific notation — `"I = 84.9 x 10^6 mm^4"`. A generic number scanner reads four separate numbers; the magnitude exists only in prose. The section label (`W310X38.7`) puts bare digits next to it. |
| `cantilever_double_integration` | Load-case token collision: UDL is `w = 11.9 kN/m`, point load is `P = 23.9 kN`. A `kN` match must exclude `kN/m` or the UDL binds as a point load — which looks like a 60% answer error, not a parse bug. |
| `annulus_flowrate` | Deliberately inverted geometry nouns: *"The inner pipe has an outer radius of 2.09 cm, and the outer pipe has an inner radius of 3.16 cm."* A regex keyed on "inner radius" grabs R_outer and inverts κ. |
| `annulus_flowrate` | Non-ASCII unit `Pa·s` (middle dot), which mojibakes in a cp1252 console. |
| `poissons_ratio` | **The load sign lives only in the prose tag** `(tension)`/`(compression)`. An oracle that ignores it silently fails ~47% of instances. |
| `poissons_ratio` | Answer precision varies 3–5 dp because trailing zeros are stripped from a fixed 5-dp rounding; part (a) is often 2 significant figures (`0.00026 in`). |
| `rotating_unbalance` | `"mass of X kg"` occurs **twice** (total mass, then unbalance mass). A naive parser silently grabs the wrong one. |
| `rotating_unbalance` | **Physics ambiguity:** the question never states whether the unbalanced mass is included in the stated total. The solution assumes it is. A solver who subtracts it is marked wrong. |
| `rotating_unbalance` / `damping_classification` | Inconsistent thousands separators *within one sentence*: stiffness `1,836,092` beside damping `10563.93`; and no commas at all in the sibling file. |
| both vibration files | Same quantity, different unit spelling: `N.s/m` vs `N-s/m`. |
| `logarithmic_decrement` | `x1` and `x_{n+1}` share the phrasing "the amplitude is measured to be N mm"; the second must be anchored on "complete cycles, " or a naive regex picks up `x1` twice. |
| several | Numbers ending a sentence (`"damping ratio of 0.457."`) — the period must not be absorbed into the token. |

---

## 6. Suggested tolerances, carried into Phase 1

Each oracle declares `TOLERANCE` with a justification in its module docstring.
Two are worth flagging as design input:

- **`annulus_flowrate` at 1e-4** — the answer is a 4-dp mantissa, so half a unit
  in the last place is ≤5e-5 relative. Measured worst error is exactly that
  display bound, i.e. every instance agrees to the last printed digit.
- **`cantilever_double_integration` at 1e-3, after quantising to the trace's
  1-dp mm display.** A plain full-precision comparison would need ≥0.017
  (measured worst 0.66%, seed 198: 6.0395 mm printed as "6.0") — too slack to be
  worth running. This is the general lesson: **compare at the precision the
  answer is quoted to, or the tolerance has to be widened until it tests
  nothing.**
