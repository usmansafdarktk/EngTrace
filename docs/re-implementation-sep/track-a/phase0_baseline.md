# Phase 0 — Baseline Measurement Report (D0.3 / D0.5)

**Deliverable:** D0.3 (baseline measurement report) and D0.5 (defect-rate reconciliation table)
**Companion to:** [`template_audit_report.md`](../audit/template_audit_report.md), [`template_redesign_spec.md`](../template_redesign_spec.md)
**Branch:** `redesign/template-integrity` · **Date:** 2026-09-05
**Scope:** every defect rate quoted in the audit report and in the spec's phase tables, re-measured independently.

> **Independence.** Every number below was produced by measurement code written for this report
> alone. `tests/template_integrity/` was **not** used to produce any figure here — that harness is
> itself under review, and the point of D0.5 is a second, unrelated number. The harness was read
> only to understand the intent of T1/T5. `data/templates/` and `evaluation/` were not modified;
> every measurement is black-box over `(question, solution)` text, except §2.9c which uses
> `sys.settrace` read-only.

---

## 1. Reconciliation table

Verdicts: **REPRODUCED** (within ±10% relative, or exact) · **REPRODUCED-APPROX** (±10–20% relative,
or reproduced under a defensibly different definition) · **NOT REPRODUCED** (>20% relative, or the
claim is false as stated) · **UNTESTABLE** (cannot be evaluated as written).
A ⚠ marks any discrepancy exceeding **20% relative**, per the exit-gate criterion.

| # | Claim (audit / spec) | Claimed | Measured | Sample | Verdict |
|---|---|---|---|---|---|
| 1 | `cantilever_double_integration` — Step-6 line does not reproduce its own arithmetic | **198/4000 = 4.95%** | **208/4000 = 5.20%**; 1,019/20,000 = **5.095%** | 4,000 / 20,000 seeds | **REPRODUCED** (+3.0% rel. at 20k) |
| 1b | …the worked example `0.01685 * 1000 = 16.8` | — | reproduced verbatim at **seed 48** | — | **REPRODUCED** |
| 2a | `rotating_unbalance` — answer recomputed from the *stated* RPM differs by up to **7.77%** | 7.77% | **7.7730%** at N=200; 8.6912% at N=4,000; **14.2734%** at N=20,000 | 20,000 seeds | **REPRODUCED** (exact at N=200; sample-size dependent — §3.1) |
| 2b | `rotating_unbalance` — ~**0.81%** error in ω | ~0.81% | max **0.8232%** (N=4,000), median 0.0253% | 4,000 seeds | **REPRODUCED** |
| 3 | `vibration_transmissibility` — up to **0.87%** worst-case TR error | 0.87% | **0.8666%** at N=200; 0.9643% at N=4,000; **1.2090%** at N=20,000 | 20,000 seeds | **REPRODUCED** (exact at N=200) |
| 4a | `damping_classification` — **175/500 (~35%)** class flips | 175/500 | **175/500 = 35.00%** (identical seeds); 1,676/5,000 = 33.52% | 500 / 5,000 seeds | **REPRODUCED** (exact) |
| 4b | …every critically-damped instance recomputes to ζ = **1.000000854527** | — | **seed 1 → 1.00000085452670**; but ζ is instance-specific and the population is ζ ∈ 1 ± ~2e-7, **straddling 1** | 5,000 seeds | **REPRODUCED as an example, NOT as a population statement** |
| 4c | …"which exact arithmetic classifies as **Overdamped**" | all flips → Overdamped | **89/175 Overdamped, 86/175 Underdamped** (866/810 at N=5,000) | 500 / 5,000 seeds | ⚠ **NOT REPRODUCED** — §3.2 |
| 4d | …step count varies 3 or 4; measured **131/200 → 3** | 131/200 | **131/200 → 3 steps**, 69/200 → 4 steps | 200 seeds | **REPRODUCED** (exact) |
| 5 | `levenspiel_plot_interpretation` — unseeded `np.random`; recorded seed cannot regenerate the item | qualitative | **3 separate processes × 4 seeds → 12/12 distinct hashes.** Corpus-wide: **200/200 seeds mismatch across two processes**, and it is the **only** one of 150 templates that does | 150 × 200 × 2 processes | **REPRODUCED** (and strengthened) |
| 6 | `composite_shafts_series` — printed `phi_AB + phi_BC` ≠ printed `phi_total` | qualitative, no rate | **89.88%** differ verbatim (5 dp vs 4 dp); **4.85%** still differ after re-rounding the printed sum to the total's 4 dp | 4,000 seeds | **REPRODUCED**; rate supplied for the first time (§2.6) |
| 7a | `impulse_response_from_lccde` — Step 5 prints `C2` at 2 dp, Answer at full precision, same trace; ~all second-order instances | "~all" | **97.8%** of second-order instances (91/200 seeds); 96.2% at N=2,000. Example `-6.02` / `-6.024922359499621` reproduced at **seed 0** | 200 / 2,000 seeds | **REPRODUCED** |
| 7b | …**80/200 seeds** show >3 dp | 80/200 | **80/200** | 200 seeds | **REPRODUCED** (exact) |
| 8a | `arl_beta_mean_shift` — **18 distinct answers in 200 seeds** | 18 | **18** | 200 seeds | **REPRODUCED** (exact) |
| 8b | `mm1_time_in_system` — **24 distinct in 200 seeds** | 24 | **24** | 200 seeds | **REPRODUCED** (exact) |
| 8c | *(not claimed)* population answer-space size | — | ARL: **19** distinct, saturated by N=1,000. MM/1: **30** at N=5,000 and still growing | 5,000 seeds | new measurement — §4 NEW-7 |
| 9a | All 150 templates import and run cleanly, zero exceptions across seeds | zero | **150 templates, 30,000 generations, 0 exceptions, 0 failed asserts, 0 non-`(str,str)` returns** | 150 × 200 seeds | **REPRODUCED** |
| 9b | **7.0%** of the **6,166** f-string interpolations are computed inline | 432/6,166 = 7.0% | population **6,166 — exact match**. Inline **449 = 7.28%** strict / 375 = 6.08% lenient. Per-branch **exact** match on chemical (65), civil (37), mechanical (89) | 6,166 interpolations | **REPRODUCED** (+4.0% rel.) |
| 9c-i | **93.0%** of interpolations are bound (complement of 9b) | 93.0% | **92.72%** | 6,166 | **REPRODUCED** |
| 9c-ii | **94.5%** of step-result numbers recoverable from bound locals | 94.5% | **93.2%** (17,002 tokens; audit counted 17,177 — population match to 1.0%) | 150 × 15 seeds | **REPRODUCED-APPROX** (−1.4% rel. on the headline; ⚠ **+24% rel. on the complement** — §3.3) |
| 9c-iii | **96.7%** of answer-line numbers recoverable | 96.7% | **97.6%** (4,619 tokens vs audit's 4,133) | 150 × 15 seeds | **REPRODUCED** |
| 9c-iv | 110/150 templates recover 100% of step values; 132/150 of answer values | 110 / 132 | **98/150** and **131/150** | 150 × 15 seeds | ⚠ **REPRODUCED-APPROX** (step count −10.9% rel.; definition-sensitive — §3.3) |
| 9c-v | "**93.6%** … of step-result numbers are recoverable" | 93.6% | **this figure appears nowhere in `template_audit_report.md`** (`grep -rn "93\.6" docs/` → no hits) | — | **UNTESTABLE — claim not located in the source document** |

**Table summary: 20 of 22 testable claims reproduce.** One is false as stated (4c); one could not
be located in the source (9c-v). No claim was overstated in magnitude; two were *understated*
(2a, 3 — §3.1) and one carried no rate at all (6).

---

## 2. Method, per measurement

All commands run from the repo root `C:\Users\ayesha.gull01\EngTrace` on **CPython 3.13.x,
Windows 11**. Every snippet is self-contained — paste into `python - <<'PY' … PY` (Git Bash) or
save and run. No repo file is written and nothing under `tests/` is imported.

A shared discovery helper, used by §2.9 and §4; save as `scratch/disc.py`:

```python
import ast, os, sys, importlib
REPO = os.path.abspath('.'); sys.path.insert(0, REPO)
BR = os.path.join(REPO, 'data', 'templates', 'branches')
def templates():
    out = []
    for dp, _dn, fns in os.walk(BR):
        if '__pycache__' in dp: continue
        for fn in sorted(fns):
            if not fn.endswith('.py') or fn in ('constants.py', '__init__.py'): continue
            p = os.path.join(dp, fn)
            rel = os.path.relpath(p, REPO).replace(os.sep, '/')
            for n in ast.parse(open(p, encoding='utf-8').read()).body:
                if isinstance(n, ast.FunctionDef) and n.name.startswith('template_'):
                    out.append((rel, rel[:-3].replace('/', '.'), n.name, n))
    return sorted(out, key=lambda t: (t[0], t[2]))
def load(mod, name): return getattr(importlib.import_module(mod), name)
```

It finds exactly **150** `template_*` functions across **47** modules — matching the audit's
population, so no template sits silently outside any measurement below.

### 2.1 Claim 1 — `cantilever_double_integration` Step-6 closure

**Definition.** The Step-6 line is `delta = <A> * 1000 = <B> mm`, `A` at 5 dp (m), `B` at 1 dp (mm).
The line **closes** iff multiplying the *printed* operand by 1000 in **exact decimal arithmetic**
and re-rounding **half-up** to the printed result's displayed precision reproduces the printed
result: `round_half_up(Decimal(A) * 1000, 1) == Decimal(B)`. `Decimal` is deliberate — the defect
is created by binary-float rounding, so the checker must not reuse it.

> **T1 as written does not catch this.** §0.2 says "within ±0.5 of the last displayed digit". Under
> that literal reading the rate is **0/4000**: every failing case sits *exactly* 0.05 mm away,
> never further. The defect is a **double-rounding disagreement**, not a magnitude error. See
> §5 SPEC-1.

```python
import sys, re, random; sys.path.insert(0, '.')
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.civil_engineering.structural_analysis.deflections import \
    template_cantilever_double_integration as T
P = re.compile(r"delta = (-?\d+\.\d+) \* 1000 = (-?\d+\.\d+) mm")
bad = 0; N = 4000
for s in range(N):
    random.seed(s); _, sol = T(); m = P.search(sol)
    A, B = Decimal(m.group(1)), Decimal(m.group(2))
    if (A * 1000).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP) != B: bad += 1
print(bad, N, 100 * bad / N)
```

**Result:** `208 4000 5.2`; at N=20,000, `1019 20000 5.095`. A shifted seed window
(`range(1, 4001)`) also gives 208, so the figure is not an artefact of the window. The
theoretical rate is exactly **5%** — failure occurs iff `D mod 0.1 ∈ [0.045, 0.05)` where
`D = delta*1000`, a 0.005-wide window in 0.1. The audit's 198 and my 208 are both draws from
Binomial(4000, 0.05) (σ = 13.8; they differ by 0.7σ).

### 2.2 / 2.3 Claims 2–3 — `rotating_unbalance`, `vibration_transmissibility` round-trip

**Definition.** Fully black-box round-trip. Parse the givens **out of the question string only**
(the solver's information set), re-derive the answer from the governing equation, compare to the
number on the `**Answer:**` line. The re-derivation is written from the physics, not copied from
the template — the anti-pattern the spec names as Phase 1's principal risk.

* `rotating_unbalance`: ω = RPM·2π/60 → F₀ = m_e·(e/1000)·ω² → X = F₀/√((k−mω²)²+(cω)²) → ×1000
* `vibration_transmissibility`: ωₙ = √(k/m), ω = 2πf, r = ω/ωₙ, TR = √((1+(2ζr)²)/((1−r²)²+(2ζr)²))

```python
import sys, re, math, random; sys.path.insert(0, '.')
from data.templates.branches.mechanical_engineering.vibrations_and_acoustics.\
     harmonically_excited_vibrations import template_rotating_unbalance as RU
Q = re.compile(r"total mass of ([\d,.]+) kg.*?stiffness of ([\d,.]+) N/m and an equivalent damping "
               r"coefficient of ([\d,.]+) N\.s/m.*?mass of ([\d,.]+) kg located at an eccentricity "
               r"of ([\d,.]+) mm\. If the machine operates at a speed of ([\d,.]+) RPM", re.S)
A = re.compile(r"amplitude of vibration is \*\*([\d.eE+-]+) mm\*\*")
errs = []
for s in range(20000):
    random.seed(s); q, sol = RU()
    mt, k, c, me, e, rpm = [float(x.replace(',', '')) for x in Q.search(q).groups()]
    p = float(A.search(sol).group(1))
    w = rpm * 2 * math.pi / 60; F0 = me * (e / 1000) * w * w
    X = F0 / math.sqrt((k - mt * w * w) ** 2 + (c * w) ** 2) * 1000
    errs.append(abs(X - p) / p * 100)
for N in (200, 500, 4000, 20000): print(N, max(errs[:N]))
```

**Result — `rotating_unbalance`:**

| N | max \|rel err\| | p99 | p95 | median |
|---:|---:|---:|---:|---:|
| 200 | **7.7730%** | | | |
| 500 | 7.7730% | | | |
| 4,000 | 8.6912% | 1.1031% | 0.4492% | 0.0285% |
| 20,000 | **14.2734%** | | | |

ω error: max **0.8232%**, median 0.0253% (N=4,000). **37.2%** of seeds have an answer error
exceeding 0.05% — the error is not a rare tail, only its *size* is rare.

**Result — `vibration_transmissibility`:** TR max **0.8666%** (N=200) → 0.9643% (4,000) →
**1.2090%** (20,000); p99 0.5438%, median 0.0610%. The **absolute amplitude X**, the item's
*second* graded answer, reaches **2.7635%** at N=4,000 — materially worse than the TR figure the
audit quoted (§4 NEW-3 context).

### 2.4 Claim 4 — `damping_classification`

**Definition.** Parse `m`, `k`, `c` (2 dp) from the question. Compute `c_c = 2√(km)` in
**50-digit `decimal`** (not float), then classify by an exact three-way comparison of the stated
decimal `c` against `c_c`. Compare to the class the gold trace asserts in Step 3. Step count is
`len(re.findall(r"\*\*Step (\d+):\*\*", sol))`.

```python
import sys, re, random, collections; sys.path.insert(0, '.')
from decimal import Decimal, getcontext; getcontext().prec = 50
from data.templates.branches.mechanical_engineering.vibrations_and_acoustics.\
     single_degree_of_freedom_systems import template_damping_classification as T
Q = re.compile(r"mass of ([\d.]+) kg, a spring stiffness of ([\d.]+) N/m, and a viscous "
               r"damping coefficient of ([\d.]+) N-s/m")
C = re.compile(r"the system is \*\*(\w[\w ]*)\*\*")
flip = 0; dirs = collections.Counter(); steps = collections.Counter()
for s in range(500):
    random.seed(s); q, sol = T()
    m, k, c = (Decimal(x) for x in Q.search(q).groups()); gold = C.search(sol).group(1)
    steps[len(re.findall(r"\*\*Step (\d+):\*\*", sol))] += 1
    cc = 2 * (k * m).sqrt()
    r = "Underdamped" if c < cc else ("Overdamped" if c > cc else "Critically Damped")
    if r != gold: flip += 1; dirs[(gold, r)] += 1
print(flip, dirs, steps)
```

**Result (N=500):** `175 Counter({('Critically Damped','Overdamped'): 89,
('Critically Damped','Underdamped'): 86}) Counter({3: 339, 4: 161})`.
N=200: 66 flips, 33/33 split, steps `{3: 131, 4: 69}`. N=5,000: 1,676 flips (33.52%), 866/810.
Seed 1 gives ζ = **1.00000085452670**, the audit's quoted value.

### 2.5 Claim 5 — `levenspiel_plot_interpretation` determinism

**Definition.** Generate in **separate OS processes** (the spec's T3 requirement), hash
`sha256(question + "\x00" + solution)`, compare. Run three times.

```bash
cat > /tmp/gen.py <<'PY'
import sys, random, hashlib, json; sys.path.insert(0, '.')
from data.templates.branches.chemical_engineering.reaction_kinetics.\
     conversion_and_reactor_sizing import template_levenspiel_plot_interpretation as T
o = {}
for s in [12345, 0, 1, 2]:
    random.seed(s); q, sol = T(); o[s] = hashlib.sha256((q + "\x00" + sol).encode()).hexdigest()
print(json.dumps(o))
PY
for i in 1 2 3; do python /tmp/gen.py; done
```

**Result:** 12 hashes, **12 distinct**. Root cause confirmed by source inspection —
`conversion_and_reactor_sizing.py:59`,
`noise = 1 + 0.15 * np.random.uniform(-1, 1, num_points)`, which `random.seed()` does not control.

**Strengthening (not claimed by the audit):** the same two-process comparison run corpus-wide
(§2.9a) shows this is the **only** non-deterministic template of 150, and that it is
non-deterministic on **200/200** seeds, not intermittently.

### 2.6 Claim 6 — `composite_shafts_series`

**Definition.** Step 4 emits `phi_total = <a> + <b> = <t> radians` with `a`, `b` at 5 dp
(`round(x, precision+1)`) and `t` at 4 dp. Two measurements: (i) **verbatim** — does
`Decimal(a) + Decimal(b)` equal `Decimal(t)` digit for digit; (ii) **closure** — does it still
differ after re-rounding the printed sum half-up to the total's own 4 dp.

**Result (N=4,000):** verbatim mismatch **3,595/4,000 = 89.88%**; closure failure
**194/4,000 = 4.85%**. Examples: seed 8 `0.07473 + 0.59572 = 0.67045` printed as `0.6704`;
seed 190 `0.00519 + 0.00556 = 0.01075` printed as `0.0107`. Every failing case is an exact
`.xxxx5` tie — **the same double-rounding mechanism as claim 1**, in a different branch.

Measurement (i) is the honest reading of the audit's sentence and is near-universal;
(ii) is what a T1 checker with a sane tolerance would flag. **Both are reported because the audit's
wording does not distinguish them and they differ by a factor of 18.**

### 2.7 Claim 7 — `impulse_response_from_lccde`

**Root cause (source-confirmed).** `discrete_time_signals.py:611`,
`fmt = lambda c, v: f"+ {c}" if c > 0 else f"- {abs(c)}"` — the helper takes a pre-formatted string
`v` and **ignores it**, interpolating the raw float `c`. The Answer line calls `fmt(C2, C2_str)`,
so `C2_str` (2 dp) is discarded and `C2` prints at `repr` precision. `C1` is unaffected (it uses
`C1_str` directly), so the defect is **asymmetric between the two coefficients of the same
expression** — a detail the audit did not record.

**Definitions.** (a) *glyph disagreement*: the Answer's `C2` digits differ from Step 5's `C2`
digits (sign handled separately, since `fmt` moves it into the operator). (b) *>3 dp*: the
Answer's `C2` literal carries more than three decimal places.

**Result:** N=200 — 93 second-order instances; **91/200 = 97.8% of second-order** glyph
disagreement; **80/200 show >3 dp** — the audit's figure exactly. N=2,000 — 96.2% and 40.0%
respectively. Seed 0 reproduces `-6.02` / `-6.024922359499621` verbatim.

### 2.8 Claim 8 — narrow answer spaces

**Definition.** The graded answer is the **last numeric literal after the final `**Answer:**`
marker**. Distinct = distinct literal.

| Template | N=200 | N=1,000 | N=5,000 | modal share |
|---|---:|---:|---:|---:|
| `arl_beta_mean_shift` | **18** | 19 | 19 | 9.2% |
| `mm1_time_in_system` | **24** | 29 | 30 | 9.8% |

`arl_beta_mean_shift` support: `{1.5, 1.6, 1.7, 1.8, 1.9, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.5, 5.9,
6.2, 6.5, 6.8, 7.1, 7.5, 8.3}` — 19 values, exactly the ceiling its own docstring declares.
`mm1_time_in_system` runs `1.9 … 60.0` and is **not** saturated at 200 seeds.

### 2.9 Claim 9 — corpus-wide

#### (a) 9a — clean import and execution, plus corpus determinism

Each of the 150 templates is called for seeds 0–199 with `random.seed(s)` and stdout redirected;
each call is wrapped in `try/except Exception`, the return asserted to be `(str, str)`, and
`sha256(q + "\x00" + sol)` recorded. The sweep is run **twice in separate processes** and the hash
maps diffed.

```bash
python scratch/sweep.py 200 /tmp/A.json     # writes {n_templates, n_runs, fails, hashes}
python scratch/sweep.py 200 /tmp/B.json
python -c "import json; A=json.load(open('/tmp/A.json')); B=json.load(open('/tmp/B.json')); \
d=[k for k in A['hashes'] if A['hashes'][k]!=B['hashes'][k]]; \
print(len(A['hashes']), len(d), sorted({k.split('|')[0] for k in d}))"
```

**Result:** `templates=150 runs=30000 failures=0` from both processes, and
`30000 200 ['template_levenspiel_plot_interpretation']`.

#### (b) 9b — AST interpolation classifier

**Population.** Every `ast.FormattedValue` reachable by `ast.walk` from the body of a `template_*`
function; module-level code and `main()` excluded. This yields **6,166** — an exact match for the
audit's population, strong evidence the two classifiers are scoped identically.

**BOUND** = `Name` · `Attribute` · `Subscript` · `Constant` · a **wrapper call**
(`round`/`abs`/`int`/`float`/`str`/`repr`/`len`) whose first argument is itself BOUND · an `IfExp`
both of whose branches are BOUND · a nested `JoinedStr` all of whose parts are BOUND.
**INLINE** = everything else: `BinOp`, `UnaryOp`, `Compare`, `BoolOp`, any non-wrapper `Call`,
`round(a*b, 2)`, comprehensions.

| Branch | interps | bound | inline | % inline | audit's inline |
|---|---:|---:|---:|---:|---:|
| chemical | 886 | 821 | **65** | 7.3% | **65** ✓ |
| civil | 1,132 | 1,095 | **37** | 3.3% | **37** ✓ |
| electrical | 1,198 | 1,071 | 127 | 10.6% | 116 |
| industrial | 1,683 | 1,552 | 131 | 7.8% | 125 |
| mechanical | 1,267 | 1,178 | **89** | 7.0% | **89** ✓ |
| **total** | **6,166** | **5,717** | **449** | **7.28%** | 432 (7.0%) |

Three of five branches match **exactly**. The 17-interpolation gap sits entirely in electrical and
industrial and is attributable to **string-formatting helper calls** — `fmt(a1, 'h[n-1]')`,
`', '.join(parts)`, `format_vector(B_vec_mT, 'mT')`, `equation_str.replace('y','h')` — which
produce *text*, not a quantity. Re-running with those treated as BOUND gives **375 = 6.08%**
(electrical 84, industrial 126). The audit's 432 sits between my two variants, so the residual is
a definitional choice about text helpers, not a disagreement about arithmetic.
Templates with zero inline interpolations: **62/150** strict, **75/150** lenient (audit: 70).

#### (c) 9c — dynamic recoverability from bound locals

**Value set V.** The generator's frame `f_locals` captured at the `return` event under
`sys.settrace`, plus the defining module's globals, recursively flattened (depth ≤ 4, ≤ 300
elements per container) through `list`/`tuple`/`set`/`dict`/`ndarray`/`Fraction`/`Decimal`/
`complex` into a set of floats. `bool` excluded; `str` **not** parsed.

**Result-position numbers.** Two populations:

* **STEP** — the region from the first `**Step` marker to the terminal answer marker; for every
  line containing `=`, **the first numeric literal after the last `=` on that line**.
* **ANSWER** — every numeric literal on or after the terminal answer marker
  (`**Answer` / `**Final Answer`).

> The "first literal after the last `=`" rule is not stated in the audit. I inferred it from the
> audit's own phrase *"the value after the last `=` on a line"* (singular) and then confirmed it
> numerically: it yields **17,002** step tokens against the audit's **17,177** (−1.0%), whereas
> taking *all* literals after the last `=` yields **29,197** (+70%). Reproducing the population
> size to within 1% is what makes the recovery percentages comparable at all.

**Classification** of token `t` printed with `d` decimals: **EXACT** if some `v ∈ V` has
`|v−t| ≤ 1e-12·max(1,|t|)`; **ROUNDED** if `|v−t| ≤ 0.5·10⁻ᵈ`; **SCALED** if `v·10ᵏ`
(k ∈ ±1…6) satisfies either; else **MISSING**.

| | tokens | EXACT | ROUNDED | SCALED | MISSING | **RECOVERED** | audit |
|---|---:|---:|---:|---:|---:|---:|---:|
| step-result | 17,002 | 56.4% | 31.4% | 5.5% | **6.8%** | **93.2%** | 94.5% |
| answer-line | 4,619 | 48.0% | 44.1% | 5.5% | **2.4%** | **97.6%** | 96.7% |

Templates recovering 100% of step tokens: **98/150** (audit 110). Of answer tokens: **131/150**
(audit 132). 15 seeds per template, 2,250 instances.

Worst step-token recovery: `system_property_linearity` (79/157 missing — a symbolic template with
no numeric payload, so expected), `levenspiel_plot_interpretation` (196/403 — its table is
`np.random`-generated and therefore genuinely unbindable), `superposition_electric_field`
(93/230), `batch_reactor_first_order` (32/90), `batch_reactor_second_order` (30/90),
`rational_method_peak_flow` (15/45), `mean_variance` (44/135), `normal_depth_iteration` (46/156).
The last four are exactly the templates the audit names as the genuine
"instrumentation-means-rewriting" cluster, so the *shape* of the finding replicates even where the
headline percentage moves.

---

## 3. Claims that did not reproduce

### 3.1 Claims 2a and 3 — the "worst case" figures are extreme-value statistics, and understated

**What happened.** 7.77% and 0.87% reproduce **to four significant figures at N=200** — my N=200
maxima are 7.7730% and 0.8666%. They are not wrong; they are the maximum of a 200-draw sample of a
continuous error distribution, and that maximum grows without bound in N:

| N | `rotating_unbalance` max | `vibration_transmissibility` max |
|---:|---:|---:|
| 200 | 7.7730% | 0.8666% |
| 1,000 | 7.7730% | 0.8871% |
| 4,000 | 8.6912% | 0.9643% |
| 20,000 | **14.2734%** | **1.2090%** |

**Why it matters.** The spec's Phase-1 table (§1.1) carries "7.77% answer error" and "0.87% TR
error" as though they were *properties of the template*. They are properties of a 200-seed sample.
Any published item pool with more than 200 instances per template contains worse cases than the
audit reports, and a tolerance set at 8% or 1% would not hold.

**Phase assignment: unchanged, and reinforced.** Both stay in Phase 1. The corrective action is
editorial — restate these as **distributional** figures. For `rotating_unbalance`: median 0.03%,
p95 0.45%, p99 1.10%, unbounded max, and **37.2% of instances exceed 0.05%**. The p99 is the
number a tolerance should be set against, not the max.

### 3.2 Claim 4c — the flip direction is a coin toss, not "Overdamped" ⚠

**Claimed:** the 175 critically-damped instances "recompute to ζ = 1.000000854527 … which exact
arithmetic classifies as **Overdamped**".

**Measured:** the flip rate (175/500) and the example ζ (seed 1) reproduce exactly. The
*direction* does not. The template sets ζ = 1.0 and derives `c = 1.0 · c_c`, then states
`round(c, 2)`. Whether the stated `c` lands above or below `c_c` is decided by the direction of
that 2-dp rounding, which is uniform:

| N | → Overdamped | → Underdamped |
|---:|---:|---:|
| 500 | 89 | 86 |
| 5,000 | 866 | 810 |

The audit generalised from a single instance. ζ = 1.000000854527 is seed 1's value; the population
is ζ ∈ 1 ± ~2×10⁻⁷, **straddling 1**.

**This makes the defect worse, not milder.** When the recomputation lands *Underdamped*, a solver
following the trace faithfully must also compute and report ω_d. The gold trace for those
instances has **only 3 steps and no ω_d**, so the correct-by-recomputation answer is not merely a
different label — it is a different *answer shape*, carrying an extra required quantity the gold
does not contain. Roughly **half of the ~33% flipping instances (≈17% of all instances)** are
unanswerable against gold in a way a label-only comparator cannot detect.

**Phase assignment: unchanged.** `damping_classification` stays in Phase 1, and the §1.3 decision
(sample `c` directly, classify by exact decimal comparison) is the right fix and resolves both
directions at once. But **D1.4's decision record must state the bidirectional flip**, or the
Phase-1 exit gate ("T2 round-trip passes at 1,000 seeds") will be written against a
one-directional oracle. **T6's step-count check acquires teeth here**: the fix must not move the
3-step/4-step split, currently **67.8% / 32.2%** at N=5,000.

### 3.3 Claim 9c — 93.2% vs 94.5%, and 98/150 vs 110/150 ⚠

**On the headline:** 93.2% vs 94.5% is −1.4% relative. That reproduces.

**On the complement — which is the number that actually sizes the work** — 6.8% missing vs 5.5%
missing is **+24% relative**, over the 20% gate, so it is flagged rather than waved through.
The 100%-recovery template count moves −10.9% relative (98 vs 110).

**Best explanation.** A matching-liberality difference, not a population difference — my token
count is within 1.0% of the audit's. Two specific candidates:

1. **My EXACT/ROUNDED split is shifted** (56.4/31.4 vs the audit's 65.9/24.4). I classify by
   numeric distance alone. A matcher that also compares the *repr string* of a bound value would
   promote many of my ROUNDED hits to EXACT **and** would catch values I miss entirely — e.g. a
   number the template holds as a **string**, which my flattener drops because it does not parse
   `str` locals. Templates that hold numbers as strings (`C1_str`, `r1_str`, `frequency_str`,
   `aql_s`, `pct`) are precisely where I would lose tokens the audit found.
2. **Depth and width caps.** Depth 4, 300 elements per container, `bool` skipped. A deeper flatten
   would recover a few more.

Both push toward my number being the **conservative** one. I have **not** forced agreement: 93.2%
is what my stated definition yields, and Phase 1 should plan against 6.8% missing rather than 5.5%.

**Consequence for phase assignments: none** — no template's class depends on this figure. It does
bear on D0.1: **T5 should match on repr-strings as well as numeric distance, and should parse
`str` locals**, or the harness will under-report binding and raise false findings against
templates that are in fact sound.

### 3.4 Claim 9c-v — "93.6%" cannot be located

`grep -rn "93\.6" docs/` returns no hits. `template_audit_report.md` line 13 states **93.0%**
(static — share of interpolations that are bound) and **94.5%** (dynamic — step-result numbers).
The pairing "93.6% / 94.5% of step-result numbers" appears to conflate the static and dynamic
measurements. Recorded **UNTESTABLE**. If a 93.6% figure exists in a per-branch audit not present
in this repo, it must be produced or the claim withdrawn, per the exit gate.

---

## 4. New defects found while measuring

None of these appear in `template_audit_report.md` or in any phase's scope table. They were
surfaced by a deliberately **narrow** printed-arithmetic probe: for every emitted line, for every
pair of adjacent `=` signs, if the segment between them matches *exactly*
`<decimal> <op> <decimal>` (optionally parenthesised) — and nothing else — check
`half_up(exact_result, dp_of_printed_result) == printed_result`. Multi-term expressions are
**skipped**, so the probe under-reports but does not false-positive. 50 seeds × 150 templates,
**7,517 probed lines**, **234 hard non-closures**, **54 boundary ties**.

> A first, looser version of this probe (two operands anywhere on the line) reported 50 templates
> and 260 failures — almost all false positives from multi-term expressions such as
> `sigma = 16.9*2.8 + 9.99*4.3 + 8.24*1.65 = 103.9`, where the regex captured only the trailing
> product. **That result is discarded and not reported here.** The lesson generalises: T1 needs a
> real expression parser, not a regex, and its "coverage %" metric must not count a line it
> mis-parses as covered.

### NEW-1 ⛔ The exemplar "fixed" template is not fixed: `template_beam_deflection_formula`

The spec's **P2** cites `civil_engineering/structural_analysis/deflections.py` as the place where
"the convention and a worked fix already exist", and the audit says the cantilever's "sibling
template 300 lines above carries the fix". That sibling is `template_beam_deflection_formula`, and
it *does* apply P2 correctly — `delta = round(delta, 5)` (line 81) before
`delta_out = round(delta * 1000, 1)` (line 82), so the printed operand **is** the stored value.

**It still fails printed-arithmetic closure on 5.89% of its SI instances.**

```
template_beam_deflection_formula: 1,969 SI instances in 4,000 seeds
  half-up   re-round of printed operands != printed result: 116 = 5.89%
  half-even re-round                     != printed result:  71 = 3.61%
  seed 18: "delta = 0.01185 * 1000 = 11.8 mm"    (0.01185 x 1000 = 11.85 -> 11.9)
  seed 32: "delta = 0.01445 * 1000 = 14.4 mm"
  seed 36: "delta = 0.03705 * 1000 = 37.0 mm"
```

Rounding to 5 dp removes the *double* rounding but leaves an exact **half-way tie**, and `round()`
resolves that tie on the **binary** value (`0.01185 * 1000` is `11.849999999999999…`), so the
printed line disagrees with a decimal reader under **both** half-up (5.89%) and half-even (3.61%)
conventions. **P2 as written is necessary but not sufficient**: the tie must also be broken
deterministically in decimal.

**This is the most consequential finding in the report**, because Phase 1's uniform transformation
(§1.2) is modelled on this file. Applied as-is to the other eight templates it reproduces this
residue eight more times. See §5 SPEC-2.

Reproduction:

```python
import sys, re, random; sys.path.insert(0, '.')
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
from data.templates.branches.civil_engineering.structural_analysis.deflections import \
    template_beam_deflection_formula as T
P = re.compile(r"delta = (-?\d+\.\d+) \* 1000 = (-?\d+\.\d+) mm")
up = ev = seen = 0
for s in range(4000):
    random.seed(s); _, sol = T(); m = P.search(sol)
    if not m: continue                      # US-customary branch prints inches
    seen += 1; A, B = Decimal(m.group(1)), Decimal(m.group(2))
    if (A*1000).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)   != B: up += 1
    if (A*1000).quantize(Decimal('0.1'), rounding=ROUND_HALF_EVEN) != B: ev += 1
print(seen, up, 100*up/seen, ev, 100*ev/seen)
```

### NEW-2 Two more back-solved templates in `torsion.py`, both out of scope

The audit flagged `composite_shafts_series` in this file and stopped. The same file holds two more
templates with the Phase-1 defect signature:

| Template | Rate | Mechanism |
|---|---:|---|
| `template_statically_indeterminate_shaft` | **47/100 probed lines (47%)** | `T_A = {T} / {1 + L_AC/L_BC:.3f} = {T_A:.2f}` — divisor printed at 3 dp, `T_A` derived from the unrounded ratio. Seed 0: `3900 / 1.291 = 3021.85`, exact **3020.91** (0.03% off) |
| `template_shaft_design_power` | **23/77 (30%)** | `d = 2 * {round(c_radius_m,4)} = {round(d_diameter_m,4)}` — diameter from the unrounded radius. Seeds 1 and 8 both print `d = 2 * 0.0186` but give `0.0371` and `0.0373` |

`shaft_design_power` is the sharper case: **identical printed operands produce different printed
answers across seeds** — self-evidently unrecoverable for a solver, and trivially detectable by T1.

### NEW-3 Six further P2 violations outside any phase's scope

| Template | Branch | Rate | Mechanism |
|---|---|---:|---|
| `template_ber_estimation_mary` | elec | 14/26 (54%) | `digital_modulation_schemes.py:424` — `Es/N0 = {k:.0f} * {eb_n0_lin:.2f} = {es_n0_lin:.2f}`; product taken from the unrounded `eb_n0_lin`. Seed 7: `5 * 63.39 = 316.93` (exact 316.95) |
| `template_gas_viscosity_kinetic_theory` | chem | 17/50 (34%) | `viscosity_and_momentum_transport.py:241` — `{sigma**2:.3f} * {omega_mu} = {denominator:.3f}`; also an **inline-computed** interpolation (T5) |
| `template_sensible_heat_temp_dependent_cp` | chem | 15/50 (30%) | `heat_effects.py:321` — `{val_T2:.4f} - {val_T1:.4f} = {val_T2 - val_T1:.4f}`; last-digit only, but the step result is computed **inline** (T5) |
| `template_vibration_isolator_design` | mech | 42/100 (42%) | already Phase 2 for other reasons; the closure defect is additional and not recorded |
| `template_work_isothermal_virial` | chem | 1/50 | `Tr = 14.73/5.2 = 2.832` (exact 2.8327) |
| `template_batch_moles_vs_conversion` | chem | 1/50 | last-digit; also emits LaTeX `\times` in an otherwise plain-text corpus |

Plus a **tie-only** cohort sharing the double-rounding mechanism of claims 1 and 6, none recorded:
`mmc_waiting_time` (8.5%), `max_hump_height_no_choking` (5.0%), `qr_policy_one_iteration` (4.0%),
`effective_stress_profile` (3.0%), `server_configuration_selection` (2.5%),
`upward_seepage_quick_condition` (2.0%), `primary_consolidation_settlement` (1.5%),
`terzaghi_strip_footing_bearing` (1.0%), and `two_state_steady_state`, `truss_method_of_joints`,
`virtual_work_truss_deflection`, `reorder_point_lead_time`, `mm1k_finite_capacity` (<1% each).

### NEW-4 `template_sigma_reduction_for_cpk` — a *deliberate* non-closure T1 will flag

42% of its probed lines fail closure, but the line announces why:
`sigma_max = 10.6 / 5.01 = 2.115 microns (rounded down, so this value still meets the target)`
(exact 2.11577). The floor is intentional and pedagogically correct. **T1 has no way to express
"this operator is a floor, not a round"** and will report 42% against a sound template. Phase 0
must give T1 a **declared-rounding-mode annotation** or an explicit skip list, or the harness
generates a false CONFIRMED finding here — and in any other template with a stated conservative
convention.

### NEW-5 `template_signal_energy_power` prints raw `repr` floats as step results

Seed 9: `Eg = 169 / 12 = 14.083333333333334` — 17 significant digits, whose final digit is a binary
artefact (`4` where exact decimal gives `3`). Three occurrences in 50 seeds. Not a correctness
defect, but it is a **machine-dependent printed value** in the sense of P4, and it is the second
instance of the `repr`-leak pattern claim 7 documents for `impulse_response_from_lccde`. Worth a
corpus-wide grep for float interpolations carrying no format spec.

### NEW-6 Source-hygiene defects

* `chemical_engineering/reaction_kinetics/stoichiometry.py:115` emits
  `f"...conversion of ${X_A*100} \%$, ..."` — an **invalid escape sequence**
  (`SyntaxWarning: invalid escape sequence '\%'` fires on every import of the corpus) and **LaTeX
  math delimiters** in an otherwise plain-text corpus; the emitted text contains a literal
  backslash.
* `industrial_engineering/quality_and_reliability_control/acceptance_sampling.py:448` contains a
  **mojibake character** inside a source string literal.
* **Non-ASCII inventory in emitted text** (50 seeds × 150 templates): U+2014 em-dash (19
  templates), U+00B7 middle dot (18), U+00B0 degree (9), U+2192 arrow (7), U+222B integral (7),
  U+00D7 multiplication sign (6), U+00B2/U+00B3 superscripts (11), U+0394/U+03A3/U+03BC/U+03C9/
  U+03C1 Greek (13), U+2082–U+2088 subscripts (3). T4 forbids non-ASCII **in numeric fields**;
  none of these sit inside a numeric field, but U+00D7 and the subscript/superscript family appear
  **adjacent** to numbers and will break a naive value extractor. T4's implementation must scope
  "numeric field" precisely or it will fail 40+ sound templates.

### NEW-7 `mm1_time_in_system`'s answer space is wider than reported; `arl`'s is not

The audit's 24 and 18 are both 200-seed counts. `arl_beta_mean_shift` saturates at **19** (its
docstring already declares this ceiling). `mm1_time_in_system` reaches **30** at 5,000 seeds and is
still growing. The two are therefore **not** equally narrow, though the audit's paired presentation
implies they are. **T6's "distinct-answer count must not fall by >10%" gate must be evaluated at
≥1,000 seeds**, not 200, or it measures sampling noise rather than the item pool.

---

## 5. Consequences for Phase 0's own deliverables

| ID | Item | Disposition |
|---|---|---|
| **SPEC-1** | **T1's tolerance is wrong.** "±0.5 of the last displayed digit" scores the cantilever defect **0/4000** — every failure sits *exactly* at 0.05, never beyond. T1 must be restated: *re-round the exact result of the printed operands, half-up in decimal, to the printed result's displayed precision, and require string equality.* Under the literal §0.2 wording, the check T1 is explicitly credited with catching does not catch it. | `SPEC-CHANGE` — amend §0.2 T1 before D0.1 is accepted |
| **SPEC-2** | **P2 is necessary but not sufficient** (NEW-1). Add: *break the tie in decimal, not on the binary float* — `Decimal(f"{x:.5f}") * 1000` quantized `ROUND_HALF_UP`, not `round(x*1000, 1)`. The exemplar file P2 cites fails at 5.89% without this. | `SPEC-CHANGE` — amend P2; re-scope `beam_deflection_formula` into Phase 1 |
| **SPEC-3** | Phase 1's table should carry **distributional** error figures (median / p95 / p99), not sample maxima (§3.1). | `SPEC-CHANGE` — §1.1 |
| **P0-1** | T5 must match on repr-strings and parse `str` locals, or it under-reports binding (§3.3). | `ADOPT-NOW` — D0.1 |
| **P0-2** | T1 needs a declared-rounding-mode annotation or skip list (NEW-4), and a real expression parser rather than a regex (§4 preamble). | `ADOPT-NOW` — D0.1 |
| **P0-3** | T4's "no non-ASCII in numeric fields" needs a precise definition of *numeric field* (NEW-6). | `ADOPT-NOW` — D0.1 |
| **P0-4** | T6 distinct-answer counts must be taken at ≥1,000 seeds (NEW-7). | `ADOPT-NOW` — D0.1 |
| **P1-1** | D1.4 must record `damping_classification`'s **bidirectional** flip and the missing-Step-4 consequence (§3.2). | `ADOPT-PHASE-1` |
| **P1-2** | Phase 1 scope grows by **3** templates (§6). | `ADOPT-PHASE-1` |
| **P5-1** | NEW-6's source-hygiene items (`\%` escape, mojibake, LaTeX in plain text) belong with Phase 5's contract work. | `ADOPT-PHASE-5` |
| **P6-1** | The NEW-3 cohort (6 hard + 13 tie templates) is a re-audit target. | `ADOPT-PHASE-6`, or `BACKLOG` in D6.6 if not scoped |

---

## 6. Phase-assignment changes

Per the spec's measurement caveat — *"If a rate does not reproduce, the affected template's phase
assignment is revisited before any edit"* — **every rate that anchors a phase assignment
reproduced. No template is removed from, or downgraded within, its assigned phase.**

Three templates should be **added** to Phase 1, on evidence gathered here rather than from the
audit:

| Template | Branch | Evidence | Proposed |
|---|---|---|---|
| `template_beam_deflection_formula` | civil | 5.89% closure failure on SI instances **despite** carrying P2 (NEW-1) | **Phase 1** — and fixed *first*, since it is the pattern the other eight copy |
| `template_statically_indeterminate_shaft` | mech | 47% of probed lines non-closing; classic rounded-operand / unrounded-result (NEW-2) | **Phase 1** |
| `template_shaft_design_power` | mech | 30%; identical printed operands yield different printed answers across seeds (NEW-2) | **Phase 1** |

That takes Phase 1 from 9 templates to **12**, and its effort estimate (15–25 h) should move with
it.

---

## 7. Reproduction environment and residual gaps

* CPython **3.13.x**, Windows 11, repo root `C:\Users\ayesha.gull01\EngTrace`, branch
  `redesign/template-integrity`; working tree clean with respect to `data/templates/` and
  `evaluation/`.
* `numpy` is imported by `conversion_and_reactor_sizing.py` only; its version affects
  `template_levenspiel_plot_interpretation`'s output but not its (non-)determinism.
* Total measurement cost: ~90,000 template generations.
* **Not established here:** hash stability across two machines (a D0.2 / Phase-0 exit-gate item).
  All 30,000 hashes above come from one machine; the two-process comparison rules out
  *within-machine* non-determinism only. A second machine is still required before the Phase-0
  gate can close.
* **Not established here:** the audit's per-branch sweeps (180,000 civil and 6,000 industrial
  instances). §2.9a covers 30,000 generations corpus-wide, which is sufficient to confirm the
  zero-exception claim but is not a replication of those sweeps.
