# Phase C3 — Reviewer H (Domain: thermochemistry / materials / fluids), INDUSTRIAL branch

Frozen ref `2fe50c72153892d0b1fa67e7ed486ebb36abca0f`. Working tree clean at that commit;
no tracked file modified. Scratch work in a session temp directory, not committed.

Reviewed without access to the implementer's reasoning. Code comments, commit messages and
phase documents were treated as claims under test.

---

## 1. Verdict

**PASS WITH FINDINGS** — the 29 corrected cells are all correct and the MIL-STD-105E
transcription is clean, but the "384 of 384 cells agree" claim is false: **two cells remain
wrong at their printed precision**, and the derivation module is structurally unable to see
them.

---

## 2. Independent re-derivation

**Method and its precision.** I did not run the implementer's module to obtain any number in
this table. I built an independent instrument:

- **c4** — closed form via `mpmath.gamma` at 40 decimal digits. No quadrature.
- **d2, d3** — my own **Gauss–Legendre tensor rule** (600 nodes per axis), against the
  module's composite Simpson. For d3 I used the **double-integral** form over the substituted
  rectangle `(x, t = y - x)`, not the survival-function form the module uses, so the two
  routes share no quadrature and no integrand.
- **Convergence**: identical to 9 decimal places across `K = 300, 600, 900, 1200` and
  half-widths `L = 8..11` — the quadrature is converged far beyond the 3–4 dp at issue.
- **Closed-form validation**: `d2(2) = 2/√π` and `d3(2) = √(2 − 4/π)` reproduced to 1e-10;
  `c4(10) = 0.9726592741` against pmc32.htm's own printed 0.9727.
- **External validation**: Harter's published `d2(5) = 2.325929`, `d3(5) = 0.864082`
  reproduced to 6 dp.
- **Third route**: `mpmath.quad` on the 1-D d2 integral gives `d2(23) = 3.85832342328501`,
  matching my Gauss–Legendre value exactly.
- **Monte Carlo** (4,000,000 samples per n, seed 20260912) as an independent sanity check:
  SE ≈ 0.0004 on d2 and ≈ 0.00026 on d3. **Stated plainly: the Monte Carlo alone does not
  settle a 3-decimal cell** — it confirms the quadrature to ~3.5 significant figures and no
  further. Every verdict below rests on the quadrature and the closed forms, not the MC.

Rounding is half-up in decimal (D-012), at each cell's **printed** number of decimal places
(the literal as written in the source, not `repr()` of the parsed float — see Finding 2).

### 2a. Sample of the 29 corrected cells — are the new digits right and the old ones wrong?

| value | method | my number | old (pre-C3.7) | new (committed) | agree? |
|---|---|---|---|---|---|
| n=18 D4 | 1 + 3·d3/d2, GL quadrature | 1.6087181 → **1.609** | 1.608 | 1.609 | ✅ new right, old wrong |
| n=19 D3 | 1 − 3·d3/d2 | 0.4035060 → **0.404** | 0.403 | 0.404 | ✅ new right, old wrong |
| n=19 D4 | 1 + 3·d3/d2 | 1.5964940 → **1.596** | 1.597 | 1.596 | ✅ new right, old wrong |
| n=22 D3 | 1 − 3·d3/d2 | 0.4345308 → **0.435** | 0.434 | 0.435 | ✅ new right, old wrong |
| n=24 D3 | 1 − 3·d3/d2 | 0.4516011 → **0.452** | 0.451 | 0.452 | ✅ new right, old wrong |
| n=2 inv_d2 | 1/d2, d2 = 2/√π closed form | 0.8862269 → **0.8862** | 0.8865 | 0.8862 | ✅ new right, old wrong |
| n=19 d3 | sd of relative range, double integral | 0.7334815 → **0.733** | 0.734 | 0.733 | ✅ new right, old wrong |
| n=19 D1 | d2 − 3·d3 (a 2-ulp move) | 1.4885185 → **1.489** | 1.487 | 1.489 | ✅ new right, old wrong |
| n=19 D2 | d2 + 3·d3 (a 2-ulp move) | 5.8894075 → **5.889** | 5.891 | 5.889 | ✅ new right, old wrong |
| n=7 inv_c4 | 1/c4, closed form | 1.0423520 → **1.0424** | 1.0423 | 1.0424 | ✅ new right, old wrong |

**All 29 corrected cells checked, not only the sample above: 29/29 — the new digit is the
correctly rounded one and the old digit was wrong.** The correction is sound. The stated
cause is also confirmed: the old values are what you get from *rounded intermediates* (e.g.
old `inv_d2(2) = 0.8865 = 1/1.128`, the rounded d2; correct is `1/1.1283792 = 0.8862`).

### 2b. Cells the phase did **not** correct but should have

| value | method | my number | committed | agree? |
|---|---|---|---|---|
| **n=6 inv_c4** | 1/c4(6), **closed form, no quadrature** | 1.0509359 → **1.0509** | **1.0510** | ❌ **wrong** |
| **n=23 D1** | d2 − 3·d3 = 3.8583234 − 3(0.7158867) | 1.7106632 → **1.711** | **1.710** | ❌ **wrong** |

### 2c. Other tables

| value | method or artefact | my number | committed | agree? |
|---|---|---|---|---|
| `Z_QUANTILES` all 6 | `statistics.NormalDist().inv_cdf(p)`, round 4 dp | 1.2816, 1.6449, 1.9600, 2.0537, 2.3263, 2.5758 | same | ✅ **6/6 exact** |
| `MIL_STD_105E_CODE_LETTERS_GII` | PDF p.18 rendered at 450 dpi, read as an image | A,B,C,D,E,F,G,H,J,K,L,M,N,P,Q over the 15 lot bands | same | ✅ **15/15** |
| `MIL_STD_105E_SAMPLE_SIZE` | PDF p.19, Table II-A sample-size column | 2,3,5,8,13,20,32,50,80,125,200,315,500,800,1250,2000 | same | ✅ **16/16** |
| `MIL_STD_105E_SINGLE_NORMAL_AC` | PDF p.19, 15 letters × 5 AQLs = 75 cells incl. arrow/None placement | see §4 | same | ✅ **75/75** |
| `SHEWHART_K_SIGMA` anchor | grep pmc32.htm | "…when k is set to 3, we speak of 3-sigma control charts" present | locator claims it | ✅ resolves |
| `XBAR_R/S_SUBGROUP_N` anchor | Montgomery p.265 (local-only) | "for large n—say, n> 10 or 12—it is probably best to use a control chart for s" present verbatim | locator claims it | ✅ resolves |
| pmc321 does **not** tabulate d2/d3 | grep pmc321.htm | "It is tabulated in many textbooks on statistical quality control" | as documented | ✅ honest |
| PDF identity | SHA-256 vs MANIFEST | `f0fcf8a6…78a27dd7`, 73 pages, MANIFEST `grounds: industrial_engineering: MIL_STD_105E_*` | — | ✅ |

**Total independently re-derived: 384 control-chart cells + 6 quantiles + 106 MIL-STD cells
= 496 values**, well past the 8 required.

---

## 3. Findings

### F1 — CONFIRMED (correctness): two cells are still wrong; "384 of 384" is false

`n=6 inv_c4` is committed as **1.0510**; the correct 4-dp value is **1.0509**.
`n=23 D1` is committed as **1.710**; the correct 3-dp value is **1.711**.

Both are the *same defect class* as the 29 that were corrected — values computed from rounded
intermediates:

```
inv_c4(6): 1/round(c4,4) = 1/0.9515   = 1.050972  -> 1.0510  (committed, wrong)
           1/c4(6)                    = 1.0509359 -> 1.0509  (correct)
D1(23):    round(d2,3) - 3*round(d3,3) = 3.858 - 3(0.716) = 1.710 (committed, wrong)
           d2 - 3*d3                   = 1.7106632        -> 1.711 (correct)
```

`inv_c4(6)` needs **no quadrature at all** — it follows from the Gamma function alone, so
this one is not arguable on numerical-tolerance grounds.

**The repo's own record corroborates this.** `phaseC1_summary.md` §10.4 states
`CONTROL_CHART_FACTORS` "derives from on-disk definitions in **353 of 384 cells** at the
table's own precision; the other **31** differ in the last digit". C3.7 corrected **29**.
31 − 29 = **2** — exactly the two I found. My independent scan of the pre-correction table
reproduces C1's number precisely: **31 cells disagree at printed places**. So the phase
inherited a correct count of 31, fixed 29, and reported the remainder as zero.

Reproduction:

```bash
git -C <repo> show b7f7613^:data/templates/branches/industrial_engineering/constants.py   # 31 bad
git -C <repo> show 2fe50c7:data/templates/branches/industrial_engineering/constants.py    # 2 bad
# derive independently (Gauss-Legendre, mpmath c4), round half-up at each cell's PRINTED places
```

**Assessed impact — on the item pool, nil; on the provenance claim, material.** Per the
template inventory and C3.4's consumer analysis, `inv_c4` and `D1` are printed and consumed
by no template (`chart_pair_selection` prints A2/D3/D4/d2 and A3/B3/B4/c4;
`xbar_r_control_limits` reads A2/D3/D4/d2; `process_capability` reads d2), and n=23 falls
outside every drawn subgroup range. **No emitted item is wrong today.** What is wrong is the
`[DERIVED]` tag's central claim. A table tagged `[DERIVED]` asserts that every cell follows
from its definition; two do not, and the commit message asserts a verification that did not
happen. Any future template that prints `inv_c4` or `D1` would inherit a wrong digit from a
table certified clean.

### F2 — CONFIRMED (method): the checker cannot see trailing-zero cells

Root cause, `tests/constants_integrity/control_chart_factors.py`:

```python
def compare(table=None, literal_tokens=None):
    ...
    tok = (literal_tokens or {}).get((n, col), repr(v))
    dp = _places(tok)
```

`run()` — the entry point, and the one the commit message cites — calls `compare()` with **no
arguments**, so `literal_tokens` is `None` and precision is taken from `repr(v)` of the
*parsed float*. `repr()` strips trailing zeros: the literal `1.0510` becomes `'1.051'`, so the
cell is checked at **3 dp instead of the 4 it is printed at**. `1.710` becomes `'1.71'` and is
checked at 2 dp.

**39 of 384 cells are checked at one fewer decimal place than they are printed at**, and the
only two that are actually wrong are both in that set. This is why `run()` reports a clean
384/384 while the table carries two wrong digits.

`literal_tokens` exists precisely to fix this — **no caller anywhere in the repo passes it**
(`git grep literal_tokens` finds only its definition and an unrelated function in
`literal_copy_sweep.py`). It is a dead parameter.

The self-test does not compensate: both planted errors (`A2(5) 0.577→0.557`,
`d3(7) 0.833→0.883`) sit in literals with **no trailing zero**, so they are caught at full
precision and the blind spot is never exercised. The suite proves it can catch a
transposition; it does not prove it can catch a last-place error in a trailing-zero cell.

Reproduction:

```bash
python -m tests.constants_integrity.control_chart_factors            # "384 agree, 0 do not"
python -m tests.constants_integrity.control_chart_factors --selftest # 0 failures
python -c "print(repr(1.0510), repr(1.710))"                         # '1.051' '1.71'
```

**Assessed impact: the gate is weaker than advertised across the whole table.** The fix is
one line at the call site — parse the source literals and pass them as `literal_tokens` — plus
two cell corrections. Both remaining cells then fall out automatically.

### F3 — PLAUSIBLE (judgement): the page-only MIL-STD locator is the honest call

I verified the transcription rather than the locator. Rendering p.18 and p.19 at 450 dpi and
reading them as images:

- **Table I, general inspection level II** — all 15 lot bands and code letters A…Q match.
  The rendered page footer reads "13", confirming the `page=18 = doc p.13` mapping; p.19
  reads "14".
- **Table II-A sample sizes** — all 16 match, including `R: 2000` (unreachable at level II,
  harmless).
- **Acceptance numbers** — all 75 cells match, *including the arrow placement*. I attempted
  to falsify this specifically (see §4) and could not.

The OCR-unreadable claim is true: p.18 extracts as
`":: .... m ·n mO ::,i:,o Vtm Loi o, batcb ah:e ."`. Given that, **a page-only locator is more
honest than a fabricated text anchor**, and I endorse the class. The residual risk is real but
correctly disclosed: nothing mechanical can re-check this table, so it rests on a human read —
which this review now independently supplies for 106 of its cells.

### F4 — PLAUSIBLE (minor): the Montgomery anchor under-grounds its own upper bound

`XBAR_S_SUBGROUP_N = (11, 25)` is anchored at `text="10 or 12"`, which grounds the *lower*
split only; nothing in that phrase justifies the upper bound of 25. The bound is in fact
grounded — the **same page** says "if we are trying to detect small shifts, then larger
sample sizes of possibly n = 15 to n = 25 are needed" — so the value is right and the anchor
is merely narrower than the claim it carries. Correctness-neutral (a sampling window whose
draw is stated in the question). Worth widening the anchor, not worth blocking.

---

## 4. Falsification attempts that failed

- **"The 29 corrections are themselves wrong."** The strongest available attack: a phase that
  "corrects" a published table is claiming the book is wrong. I re-derived all 29 from
  definitions with an instrument that shares no code with the implementer's, and **29/29 the
  new digit is right**. The rounded-intermediate explanation also holds arithmetically on
  inspection of the old values. Could not break it.
- **"The d3 double integral is wrong."** The module's docstring admits an earlier version got
  d3 wrong by up to 13% through an invalid Simpson rule on a triangular region — so I assumed
  the current one might still be off. I used a different formulation (substituted rectangle,
  Gauss–Legendre) and got agreement to 9 dp, validated against `d3(2) = √(2 − 4/π)` and
  Harter's `d3(5) = 0.864082`. The current d3 is correct.
- **"MIL-STD Table II-A has a sandwiched-arrow transcription error."** `F: {0.65: 0, 1.0: None,
  2.5: 1, ...}` looked impossible — a direct entry, then an arrow, then a direct entry, with Ac
  monotone in AQL. I expected a transcription defect. **The standard genuinely does this**: the
  master table staggers a ⇧ "use first plan above" arrow between direct entries. Reading the
  rendered page confirmed the committed values exactly. My objection was wrong, the
  transcription was right.
- **"Row A's Ac=0 sits at the wrong AQL."** I initially read the A/B/C diagonal backwards off
  the OCR and expected `A: {6.5: 0}` to be an arrow. At 450 dpi the stagger is unambiguous —
  C(n=5) at AQL 2.5, B(n=3) at 4.0, A(n=2) at 6.5, exactly as committed.
- **"`Z_QUANTILES` is a table lookup dressed as a derivation."** All 6 reproduce
  `inv_cdf` exactly at 4 dp. The `[DERIVED]` tag is honest.
- **"The NIST anchors don't exist."** Both resolve verbatim: pmc32 carries "we speak of 3-sigma
  control charts" and its own "0.9727" check value; pmc321 carries "tabulated in many
  textbooks" and indeed tabulates neither d2 nor d3 — which the grounding doc already states.
- **"The Montgomery p.265 anchor is invented."** It is verbatim on that page of the local-only
  copy. (No content copied or committed.)
- **`[UNVERIFIED]` on IEC60063 / HARD_ANODIZE_THICKNESS.** Neither standard is on disk; I could
  not check the values and neither could the implementer. The tag is the honest class. The E24
  and E12 lists are internally consistent with the standard series as far as a domain reader
  can judge, but "consistent with what I remember" is not verification and I record it as such.

---

## 5. Further probing and improvements (separately time-boxed, 15 min)

1. **Fix the precision blind spot, not just the two cells.** Correcting `inv_c4(6)` and
   `D1(23)` without fixing `compare()` leaves the gate able to hide the next one. Pass the
   source literals as `literal_tokens` from `run()`, and add a planted-error self-test case
   **in a trailing-zero cell** — the current two planted cases cannot fail on this axis.
2. **Make the "31" reconcile explicitly.** C1 published 31 bad cells and C3 fixed 29 without
   remarking on the gap. A phase that inherits a count from an earlier phase should assert
   the arithmetic (`fixed + remaining == inherited`) rather than re-measure with a weaker
   instrument and report zero.
3. **Re-check the other `[DERIVED]` tables for the same repr() pattern.** The weakness is in a
   shared idiom (`_places(repr(v))`), not in this table. Any other integrity check that infers
   precision from a parsed float has it too — worth one grep across
   `tests/constants_integrity/`.
4. **Record the MIL-STD human verification durably.** Since no mechanical check is possible,
   the page-only locator should carry a note that the transcription was verified against the
   rendered page (by C3.7's author and now independently here, 106/106 cells), so a later
   reader knows the difference between "unchecked" and "checked by eye, twice".
5. **Widen the `XBAR_S_SUBGROUP_N` anchor** to the "n = 15 to n = 25" sentence on the same
   page, which actually grounds its upper bound.
6. **Consider carrying d2/d3 at more places.** Several cells sit within ~2e-4 of a rounding
   boundary (`D1(23)` is 1.71066 against a 1.7105 boundary). The table is right at 3 dp, but a
   future consumer computing limits from these would do better reading a higher-precision
   source column than re-deriving from 3-dp values — which is the very defect this phase was
   correcting.

---

*Filed by Reviewer H (industrial). Mandatory task completed within its 40-minute box;
§5 within its separate 15-minute box.*
