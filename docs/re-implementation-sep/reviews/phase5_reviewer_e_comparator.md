# Phase 5, Track B — Reviewer E (comparator adversary)

**Ref reviewed:** `7c04c95` (working tree is identical to it; `git diff 7c04c95 HEAD` is empty)
**Scope:** comparator bindings — over-acceptance and over-rejection. T4 / template-integrity,
the corpus baseline, D-056, `narrative`, `evaluation/` and the annotation pilot are out of scope
and were not re-derived.
**Protocol:** R1 — every number below was re-derived from the code, not read from a claim.

---

## 1. Verdict

**BLOCKED** — 75 of the 132 bound templates do not return `MATCH` when the candidate is a
verbatim copy of gold; the gate reports "0 false accepts" because the comparator has stopped
accepting, and two real gold pairs show cross-instance acceptance sitting underneath.

---

## 2. Independent re-derivation

All commands were run from the repo root with `PYTHONIOENCODING=utf-8`.

| # | claim | re-derived | agreement |
|---|---|---|---|
| **E1** | gold×gold 17,424 pairs / 132 scored / 18 skipped / 0 FA / 0 errors | `python -m tests.comparators.cross_pair` → `pairs 17424 / templates-scored 132 / templates-skipped 18`, `FALSE ACCEPTS 0`, `ERRORS 0` | **agrees, exactly** |
| **E2** | archive×gold 20,865 (1,503 + 19,362) / 70 scored / 14 skipped / 0 XI / 0 errors | same command → `pairs 20865 (1503 positive + 19362 negative) / templates-scored 70 / templates-skipped 14`, `cross-instance accepts 0 / 19362`, `ERRORS positive 0 negative 0` | **agrees, exactly** |
| **E3** | archive covers 84 of 150; 6 stems map to no template, 259 traces | same command → 70 scored + 14 skipped = **84**; `archive stems with no template: 6 (259 traces)` — `coaxial_cable_capacitance`, `coulombs_law_two_charges`, `finding_limiting_reactant`, `flow_rates_vs_conversion`, `gauss_law_symmetric_charge`, `mean_variance_discrete_rv` | **agrees, exactly** |
| **E4** | 132/150 bound, N=50, 2,450 pairs each, 340,550 pairs, zero FA and zero errors per binding | `python -m tests.comparators.derive_bindings` → `BOUND 132 / 150 over 340550 validated pairs`; `committed tables agree with a fresh derivation` | **arithmetic agrees; the word "bound" does not.** `validate()` gates on `false_accepts or errors` and never on `decided`, which it computes and discards. 11 of the 132 have `decided: 0` — 2,450 of 2,450 pairs `UNRESOLVED`. See finding **E-6**. |
| **E5** | N1 held-out: incumbent 16 ggFA → 0; archive decided 99.1% → 99.1%; "last number" 424 | `python -m tests.comparators.n1_candidates` held-out table → `C1_first_INCUMBENT 16 / 100.0 / 144 / 99.1`; `C10_ADOPTED` and `SHIPPED_per_kind` `0 / 100.0 / 6 / 99.1`; `C2_last 424` | **agrees, exactly** |
| **E6** | D4.4's 21 `check` cases: first-number 13 decided / 13 correct; every last-ward rule 13 / 5 / 8 FR | same table → `C1_first 13/21, 13, 0`; `C2_last`, `C6`, `C7`, `C9`, `C10`, `SHIPPED` all `13/21, 5, 8` | **agrees**, with one thing the phrasing hides: the *non*-last-ward rules `C3`, `C4`, `C5`, `C8`, `C11` score `18/21 decided, 14 correct, 0 false rejects` — strictly more decided and more correct than the adopted rule, at the price of 4 false accepts. The comparison class chosen for the claim ("every last-ward rule") is the one that flatters the adopted rule. Defensible under the §1 accept-averse bias; not stated. |
| **E7** | carries-a-unit: seed0 115 / any 116 / always 113 / invariant 103; 103 units declared | `derive_bindings` → `seed0 115 / any 116 / always 113 (ADOPTED) / invariant 103`; `units DECLARED … 103`; `len(DECLARED_UNITS) == 103` | **agrees, exactly** |
| **E8** | N1/N2/N4 fixed: `hagen_poiseuille_flowrate` 3 vs 9 now MISMATCH; `autocorrelation_rect_pulse` no longer MATCHes across instances; `continuous_to_discrete_conversion` no longer raises | see below | **diverges on all three** |
| **E9** | D4.4 unchanged from Phase 4 at 95.8% / 98.6%; archive 100% / 97.8% with 0 FA | `python -m tests.comparators.score` → adversarial `precision 95.8% recall 98.6%`, archive `precision 100.0% recall 97.8% false accepts: 0` | **numbers agree.** "Unchanged from Phase 4" does not: `phase4_comparators.md` §8 *at this same ref* still prints `98.6% / 100%` for the same corpus. One of the two is stale. |

### E8, item by item

```bash
python - <<'PY'
import sys; sys.path.insert(0,'.')
from tests.template_integrity.core import discover, generate
from tests.comparators.bindings import BINDINGS, UNBOUND, compare_template
from tests.comparators.answer import compare_kind
refs={r.template_id:r for r in discover()}
for tid in ['template_hagen_poiseuille_flowrate','template_autocorrelation_rect_pulse',
            'template_continuous_to_discrete_conversion']:
    a=generate(refs[tid],3,capture=False).solution; b=generate(refs[tid],9,capture=False).solution
    print(tid, 'BOUND' if tid in BINDINGS else 'UNBOUND')
    try: v=compare_template(tid,b,a); print('  compare_template:',v.outcome,'|',v.reason)
    except Exception as e: print('  compare_template RAISES',type(e).__name__,e)
    print('  compare_kind numeric:',compare_kind('numeric',b,a).outcome)
PY
```

* **`hagen_poiseuille_flowrate` — not reproducible as stated.** The template is **UNBOUND**
  (`the asserted-number count varies across instances: [1, 2]`), so `compare_template` raises
  `KeyError`. The N1 fix *is* real and visible through the default path
  (`compare_kind('numeric', …)` → `MISMATCH |0.001183 − 0.002738| = 0.001555 > 5E-7`), but the
  template the claim exhibits it on is one Track B declined to bind. The headline defect's own
  witness is outside the shipped binding table.
* **`autocorrelation_rect_pulse` — true at the template level, untrue of the rule.**
  `compare_template` returns `UNRESOLVED`, so it does not MATCH across instances. But
  `compare_kind('numeric', seed9, seed3)` still returns **`MATCH`** for two answers whose peaks
  are 250 and 810: both spans end `… at tau = 0.`, the anchor rule takes the number after the
  last `=`, and both parse to `0`. That is N2's shape — two different answers reduced to the same
  degenerate parse — surviving in the adopted rule. At the template level it is suppressed by the
  declared "unit" `'otherwise'` (finding **E-3**), i.e. by a defect.
* **`continuous_to_discrete_conversion` — no longer raises, and now decides nothing.**
  2,450 of 2,450 validation pairs `UNRESOLVED`, gold against a verbatim copy of itself included.
  N4 satisfied "zero false accepts" by crashing; it now satisfies it by refusing. `cross_pair.py`'s
  own docstring says *"a template that raises has not passed"* — nothing says the same of a
  template that decides nothing, and finding **E-6** is that gap.

---

## 3. Findings

### 3A. Over-rejection (Direction 2)

#### E-1 — CONFIRMED — 75 of 132 bound templates reject a verbatim copy of gold

The sharpest available correct answer is gold itself. It is not credited.

```bash
python - <<'PY'
import sys, collections; sys.path.insert(0,'.')
from tests.template_integrity.core import discover, generate
from tests.comparators.bindings import BINDINGS, compare_template
refs={r.template_id:r for r in discover()}
bad=collections.Counter()
for tid in sorted(BINDINGS):
    g=generate(refs[tid],3,capture=False).solution
    try: o=compare_template(tid,g,g).outcome
    except Exception: o='ERROR'
    bad[o]+=1
print(bad)   # Counter({'MATCH': 57, 'UNRESOLVED': 45, 'MISMATCH': 30})
PY
```

Measured over 5 seeds each: **57 templates MATCH themselves on every seed, 75 do not on any.**
The breakdown at seed 3 is `37` UNRESOLVED and `17` MISMATCH from the unit check (**E-2**),
`13` MISMATCH from multipart part-units (**E-4**), and `8` UNRESOLVED from `symbolic`.

`gold_gold`'s own loop skips `a == b`, so this case is not merely unmeasured — it is
**excluded by construction**. Every one of these 75 is a comparator that will refuse or reject the
correct answer for that item, whatever a model writes.

**Impact:** the shipped comparator scores at most 57 of 150 templates in a way a correct solver
can pass. That is worse than the 4 Phase 4 bound, on the measure that matters.

#### E-2 — CONFIRMED — the unit resolver knows 11 units; 72 of the 103 declared units are outside it

`kinds.UNIT_SURFACES` has 11 canonical units (`Hz L N Pa kPa m mL mm rad/s s`).
`bindings.DECLARED_UNITS` declares 48 distinct values across 103 templates, of which **40 are
unresolvable** (`J/mol K L/mol MPa N.m N.s/m N/C N/m Pa*s V/m bar cm/s cm^3/mol deg degrees dollars
e-05 hours kJ kN*m liters m/s m/s^2 m^3 m^3/kg m^3/s minutes mol mol/dm^3 mol/min otherwise pF/m
percent rad/m rad/sample seconds subgroups units weeks years`), covering **72 of the 103** templates.

`compare_numeric` runs the unit check *after* the numeric comparison, so a numerically correct
answer is then thrown away two ways:

* `_resolve_unit(candidate)` returns `None` → **UNRESOLVED** ("gold declares the unit 'K' and the
  candidate states none") — even though the candidate is gold;
* `_resolve_unit(candidate)` matches a *substring* of the declared unit → **MISMATCH**:
  `unit 'm' != declared gold unit 'm^3/kg'`, `unit 'N' != 'N/m'`, `unit 'kN' != 'kN*m'`,
  `unit 's' != 'seconds'`, `unit 'L' != 'liters'`, `unit 'm' != 'm/s'`.

The MISMATCH half is the more serious: `UNRESOLVED` declines, `MISMATCH` reports a model error
that did not happen, on a template where gold and candidate are the same string.

```bash
python -c "import sys;sys.path.insert(0,'.');\
from tests.comparators.kinds import UNIT_SURFACES;\
from tests.comparators.bindings import DECLARED_UNITS as D;\
print(sorted(set(D.values())-set(UNIT_SURFACES)));\
print(len([t for t,u in D.items() if u not in UNIT_SURFACES]),'of',len(D))"
```

#### E-3 — CONFIRMED — `DECLARED_UNITS` is a trailing-token heuristic, and it re-creates the exact rejection §4.1 was designed to prevent

`derive_bindings.unit_token` takes the last non-stopword token on a line of the answer span. What
that yields is not units:

| template | declared "unit" | gold span |
|---|---|---|
| `autocorrelation_rect_pulse` | `otherwise` | `… for |tau| <= 10, and 0 otherwise.` |
| `euclidean_distance_binary` | `e-05` | `s1 = (2.435e-05)` |
| `basic_eoq`, `newsvendor_normal_demand`, +6 | `units` | `The optimal order quantity is 1005 units` |
| `chase_vs_level_aggregate` | `dollars` | `… total three-month cost is 9240 dollars` |
| `takt_time_line_efficiency`, +2 | `percent` | `… is 89.0 percent` |
| `arl_beta_mean_shift` | `subgroups` | `… run length is 6.2 subgroups` |

D4.1 §4.1 is explicit about why the unit check is opt-in and declared: *"Inferring the unit from
gold's string instead would reject `7.65 litres`, which is correct. That trade is why it is
opt-in."* D5.10 infers it from gold's string, and the rejection has arrived, on real archived
model text:

```
template_batch_reactor_second_order   instance 2504064379
  GOLD span : 'The required reaction time is 17.46 seconds.'
  CAND span : '17.46 s'
  SHIPPED   : MISMATCH | unit 's' != declared gold unit 'seconds'

template_vibration_isolator_design    instance 1756964932
  GOLD span : '… is 1,256,557 N/m.'
  CAND span : '1256556.55 [N/m]'
  SHIPPED   : MISMATCH | unit 'N' != declared gold unit 'N/m'
```

The second is the same unit, spelled the same way, rejected because `_resolve_unit` finds the
`N` inside `N/m` and stops.

```bash
python - <<'PY'   # reproduction
import sys; sys.path.insert(0,'.')
from tests.comparators.cross_pair import load_archive
from tests.comparators.bindings import compare_template
rows,_=load_archive()
for inst,cand,g in rows['template_batch_reactor_second_order'][:40]:
    v=compare_template('template_batch_reactor_second_order',g,cand)
    if 'declared gold unit' in (v.reason or ''): print(inst, v.outcome, '|', v.reason); break
PY
```

#### E-4 — CONFIRMED — `numeric_parts(n, unit)` gives every part of a multipart answer the same unit

```python
Part(name=f"q{i+1}", kind="numeric", options={"extract": False, **({"unit": unit} if unit else {})}, …)
```

A multipart answer's parts carry *different* units — that is why they are different parts. One
declared unit applied to all `n` of them is wrong for at least `n−1`:

```
template_volumetric_flow_rate   unit 'm/s'      gold vs itself → MISMATCH
    q1  unit 'm' != declared gold unit 'm/s'      (0.0068 m^3/s)
    q2  unit 'm' != declared gold unit 'm/s'      (3.065 m/s)
template_wave_parameters_basic  unit 'rad/m'     gold vs itself → MISMATCH on all 4 parts
template_system_properties      unit 'N.s/m'     gold vs itself → MISMATCH on all 3 parts
```

13 multipart templates MISMATCH a verbatim copy of gold for this reason alone.

```bash
python -c "import sys;sys.path.insert(0,'.');\
from tests.template_integrity.core import discover,generate;\
from tests.comparators.bindings import compare_template;\
refs={r.template_id:r for r in discover()};\
g=generate(refs['template_volumetric_flow_rate'],3,capture=False).solution;\
v=compare_template('template_volumetric_flow_rate',g,g);\
print(v.outcome,'|',v.reason);[print(' ',p.kind,p.outcome,p.reason) for p in v.parts]"
```

#### E-5 — CONFIRMED — the unit check costs 19 of the 82 matches on real archived model answers

Scoring the 1,503 archive positives (a real model answer against its own gold) with the shipped
table, then with `unit` stripped from every binding and nothing else changed:

| | MATCH | MISMATCH | UNRESOLVED |
|---|---:|---:|---:|
| shipped | **82** | 771 | 650 |
| unit check removed | **99** | 800 | 604 |

Flips: 14 `UNRESOLVED → MATCH`, 5 `MISMATCH → MATCH`, 2 `MATCH → MISMATCH` (the check does catch
two things), 32 `UNRESOLVED → MISMATCH`. **19 real archived answers that agree with gold on the
number, the precision and the unit are refused or actively marked wrong**, against a baseline of
82 accepted answers in total. Worked examples, all real:

```
template_gas_phase_concentration  (14 of the 19)
  GOLD : 'The outlet concentrations are:\n- C_A : 0.042 mol/dm^3\n- C_B :** 0.082 mol/dm^3'
  CAND : 'C_A = 0.0420 mol/dm^3 and C_B = 0.0820 mol/dm^3'
  SHIPPED: UNRESOLVED  ("gold declares the unit 'mol/dm^3' and the candidate states none")
template_newtons_law_shear_stress
  GOLD : 'a) … 0.028 Pa.\nb) … 0.070 N.'
  CAND : '- Answer (a): 0.0280 Pa\n- Answer (b): 0.0704 N'
  SHIPPED: MISMATCH   (declared unit 'N' applied to the Pa part as well — E-4)
```

Reproduction: `scratchpad/revE/p_unitcost.py` (swaps `bindings.BINDINGS` for a unit-free copy and
re-runs `load_archive()`; no repo file is touched).

#### E-6 — CONFIRMED — 11 templates are counted as "bound" while deciding nothing

`derive_bindings.validate()` computes `r["decided"]` and then gates on
`if r["false_accepts"] or r["errors"]` only. These 11 are in `BINDINGS`, in the 132, and in the
340,550 "validated pairs", with `{'pairs': 2450, 'unresolved': 2450, 'decided': 0}`:

`autocorrelation_rect_pulse`, `ber_estimation_mary`, `bpsk_energy_basis`, `cd_dc_system_analysis`,
`coaxial_capacitance`, `continuous_to_discrete_conversion`, `ft_esd_rect_pulse`, `phasor_addition`,
`standing_wave_formation`, `undamped_response_initial_conditions`, `vdw_solve_for_volume`.

`sympy 1.14.0` is installed, so this is not the D-053 dependency gap. Causes are
`SympifyError` (4), `symbol(s) outside the declared alphabet: ['c','m','o','s']` — the alphabet
`('t','f','x','n','tau')` against gold that says `cos`, `seconds`, `deg` (3), a gold/candidate
tuple parse (1), and the E-3/E-4 unit defects (3).

The same shape is visible in the phase's own archive numbers and is not in the claim table:
**decided rate on the positives 853/1503 = 56.8%** — `symbolic` **3.9%**, `vector` **12.0%**,
`multipart` **43.8%**. A binding that refuses is not a binding.

```bash
python -c "import sys;sys.path.insert(0,'.');\
from tests.comparators.bindings import VALIDATION as V;\
print([t for t,v in V.items() if v.get('decided')==0])"
```

#### E-7 — CONFIRMED — the shipped tool already prints 56 false rejects and they are not in the claim table

`cross_pair` reports `false rejects 56` beside `FALSE ACCEPTS 0`. There are 198 gold pairs at
n=12 whose spans are textually identical; **56 of them (28%) are not MATCHed**. At `-n 24`,
72,864 pairs, the figure is **176**. Every one is a correct answer refused, and every cause is
E-2/E-4/E-6:

```
adiabatic_flame_temperature   UNRESOLVED  gold declares the unit 'K' and the candidate states none
batch_reactor_second_order    MISMATCH    unit 's' != declared gold unit 'seconds'
ft_esd_rect_pulse             UNRESOLVED  symbol(s) outside the declared alphabet: ['c','i','s']
mm1_time_in_system            UNRESOLVED  gold declares the unit 'minutes' …
```

`E1` is quoted as a five-part result with the false-reject column dropped. Reporting
`FALSE ACCEPTS 0` beside a suppressed `false rejects 56` is the accept-only framing the brief
warns about, inside the tool that was built to avoid it.

### 3B. Over-acceptance (Direction 1)

#### E-8 — CONFIRMED — `nth_quantity`'s 16-character window makes a part compare a *different* part's number, and two real gold pairs realise it

`nth_quantity(text, i)` returns `t[end(i−1) : end(i) + 16]`, and the part comparator then applies
the answer rule to that slice — `_after_anchor` or `_last_filtered`, neither of which is
"the i-th number". When number `i+1` falls inside the 16-character tail, part `i` selects it, two
parts compare the same quantity, and **one asserted quantity is compared by nothing**.

**6 of the 31 multipart bindings** do this on real gold:

| template | symptom |
|---|---|
| `hydrostatic_pressure_at_depth` | part1 selects the pressure `359.181`, not the depth `19.43` → the depth is never compared |
| `rackett_equation_volume` | part1 selects the volume `170.9`, not the temperature `518.27` |
| `particle_pathline` | part1 selects the *given* time `2.5`; part2 selects `y`; the x-coordinate is compared by nothing |
| `time_to_phasor` | part1 selects the angle, not the magnitude |
| `statically_indeterminate` | part1 selects `R_C`, not `R_A` |
| `autocorrelation_rect_pulse` | part5 selects the constant `0`, not the peak `250` |

Two real gold instances of the same template, generated by the shipped generator, are then accepted
for each other:

```
template_hydrostatic_pressure_at_depth   seed 157 vs seed 237
  gold : The absolute pressure at a depth of 13.29 m is **101.347 kPa**.
  cand : The absolute pressure at a depth of 13.38 m is **101.347 kPa**.
  verdict with the E-3 pseudo-unit removed: MATCH   <- cross-instance accept

template_rackett_equation_volume         seed 105 vs seed 206
  gold : … saturated liquid Carbon dioxide at 199.89 K is **36.22 cm³/mol**.
  cand : … saturated liquid Methane      at 102.06 K is **36.22 cm³/mol**.
  verdict with the E-3 pseudo-unit removed: MATCH   <- cross-instance accept
```

**As shipped both return `MISMATCH` — because of the unit defect, not because of the depth or the
temperature.** The two defects mask each other: the reported "0 false accepts" on 340,550 pairs is
produced by E-2/E-3/E-4, and repairing them uncovers this. That is the load-bearing point of this
report. Reproduction: `scratchpad/revE/p_final.py` and `p_collide.py` (a 250-instance sweep that
found the two colliding pairs; `N=50` cannot resolve a 1-in-31,000 collision, and the phase says so
itself — "smallest resolvable rate ~0.122%").

#### E-9 — CONFIRMED — three multipart bindings contain parts that read a constant

A part whose selected number is the same for every instance contributes nothing to any verdict:

```
autocorrelation_rect_pulse  n=6   parts 4, 5, 6 all read 0   (half the binding is vacuous)
coaxial_capacitance         n=3   part 1 reads 2
vdw_solve_for_volume        n=4   parts 1 and 3 read 1 and 2
```

All three are currently in the E-6 "decides nothing" set, so the vacuity is invisible. It becomes a
false accept the moment E-2/E-3 are repaired.

```bash
python "scratchpad/revE/p_const.py"   # 12 seeds per template; prints the three above
```

#### E-10 — PLAUSIBLE — `partial` cannot fire on a numeric or symbolic binding, and two symbolic bindings are narrower than their answers

`derive_all` assigns `partial` only when the kind is `categorical`/`categorical[tuple]` with numbers
present, or when the inventory `answer_type == "classification"`. It is on exactly two templates.
It cannot be assigned when a genuinely multi-part answer takes a single-part binding:

```
template_bpsk_energy_basis      symbolic, unflagged
  gold: 'a) The energy per bit (Eb) is 26.83 mJ.\nb) The basis function … 19.17 * cos(2*pi*2.39 kHz*t).'
template_standing_wave_formation symbolic, unflagged
  gold: three assertions — the wave equation, three nodes, three antinodes
```

Both currently return `UNRESOLVED` on everything (E-6), so nothing is over-accepted **today**. I
checked the other suspects and they are clean: for every `numeric`-bound template the gold answer
span really does assert exactly one quantity, over 12 instances, including the ones the inventory
types `multipart`/`array`/`vector` (`aoq_ati_rectifying`, `arl_beta_mean_shift`,
`chase_vs_level_aggregate`, `server_configuration_selection`, `line_balancing_heuristic`,
`normal_depth_iteration`, `two_state_steady_state`). `bindings.py`'s argument that the inventory
type is not the kind holds where I could test it.

#### E-11 — PLAUSIBLE — no template is bound `check`, so E6's per-kind rule is unexercised

`collections.Counter(b['kind'] for b in BINDINGS.values())` →
`numeric 85, multipart 31, symbolic 9, vector 3, sequence 2, categorical[tuple] 1, categorical 1`.
Zero `check`, zero `narrative`. The whole per-kind reframing that `extract.py`'s docstring calls
"the whole reframing" is validated only on D4.4's 21 hand-written cases, and D4.1 §5 already flags
that those 21 `check` surfaces are the one part of the spec written from the domain rather than
observation. `vorticity_check` — the one template whose name says `check` — is unbound.

#### E-12 — CONFIRMED — the textual-identity truth predicate invents false accepts, and templates may have been declined to bind because of it

`gold_gold` and `validate` call a MATCH a false accept iff the two spans differ *textually*. At the
display-tolerance boundary that is wrong in the accept-blocking direction:

```
template_batch_reactor_second_order   gold 'The required reaction time is 4.8 seconds.'
                                      cand 'The required reaction time is 4.85 seconds.'
```

Gold shows one decimal, so §4.1's tolerance is `0.5 × 10⁻¹ = 0.05`, and `|4.85 − 4.8| = 0.05` is a
**correct MATCH, boundary inclusive** — the proxy scores it a false accept. Seven of the eighteen
`UNBOUND` entries give a reason of the form `N false accepts in 2450 pairs`
(`decimation_aliasing_analysis` 2, `floating_object_submersion_depth` 4, `gauss_law_symmetric` 4,
`null_to_null_bandwidth` 24, `pitzer_correlation_z` 1, `signal_energy_power` 2,
`truss_method_of_sections` 2). Small counts of exactly this shape. **None of them has been audited
against the display tolerance**, and the docstring's promise — *"the proxy is audited rather than
believed"* — is implemented only as `false_accept_examples` in the gold×gold report, which is empty
because the shipped bindings accept nothing. Templates may be sitting in `UNBOUND` for behaving
correctly.

I found the other direction too, and it is real: `vdw_solve_for_pressure` gold `n-Butane … 4.9 bar`
against candidate `n-Octane … 4.85 bar` is within tolerance and *is* a genuine cross-substance
accept that only the number check can see — but the proxy cannot tell the two cases apart, and the
seven UNBOUND reasons are therefore unreadable as they stand.

---

## 4. Falsification attempts that failed

These are the things I tried to break and could not. They are what makes §3 mean anything.

1. **gold×gold at double the instance count.** `-n 24`, 72,864 pairs over all 132 bindings:
   `MISMATCH 65943, UNRESOLVED 6289, MATCH 632`, **0 false accepts**. E1's headline holds at
   n=24, not only at n=12.
2. **archive×gold cross-instance accepts, with and without the unit check.** 19,362 negatives
   both ways: **0 and 0**. Removing the defect that suppresses decisions did *not* produce a
   single cross-instance accept on real model text. The E-8 over-acceptance is reachable from
   gold×gold and not from the archive.
3. **`extract.displayed_decimals`.** I tried to make it absurd in both directions and could not:
   `3.1259e-07 → 11`, `4.62e-09 → 11`, `1.5e3 → −2`, `1.36e+09 → −7`, `3e-07 → 7`, `1e12 → −12`,
   `0.5e-3 → 4`, `4,921 → 0`, `2.50 → 2`. Every one is the exact half-unit-in-the-last-place of
   the digits displayed; the negative results are correct (`1.36e+09` shows to the ten-millions,
   so `0.5 × 10⁷` is right, not loose). The claim that closing this recovered decided rate is
   consistent with what I see. **No finding.**
4. **The two `partial`-flagged templates on real archived text.** 57 archived
   `reynolds_number_flow_regime` answers scored; `MATCH 10 / MISMATCH 39 / UNRESOLVED 8`, and
   **not one** MATCH where gold's regime word and the candidate's differ. The flagged
   over-acceptance is real in principle and not realised in the corpus.
   `critical_depth_froude_classification` has no archived traces.
5. **A colliding compared-tuple for the other four mis-slicing templates.** 250 instances each:
   `particle_pathline` 0, `statically_indeterminate` 0, `time_to_phasor` 0. E-8 is proven on two
   templates, not six; on the other four the mechanism is present and the collision is rarer than
   1 in 62,250.
6. **An archived model answer that realises E-8.** For all five templates I looked for a MATCH
   whose asserted-number list differs from gold's. **None**, over 188 archived rows.
7. **An unflagged narrow `numeric` binding.** For every one of the 85 `numeric`-bound templates,
   the gold answer span asserts exactly one number across 12 instances, and the non-numeric
   residue varies only in the noun (`rod`/`pipe`/substance name). The `bindings.py` docstring's
   claim that the inventory `answer_type` is not the kind survives this test.
8. **E1, E2, E3, E4, E5, E6, E7, E9 as arithmetic.** Every number re-derived; every one agrees.
9. **`reviewer_battery` (81 cases) and `recall_corpus` (507 cases)** both run clean at this ref.

---

## 5. Further probing and improvements

**Where the evidence is thinnest.** The phase measured acceptance on 340,550 pairs and
resolution on none. Every gate item in D5.8 is a false-accept count; `decided` is computed in
`validate()` and thrown away. That single line is why a comparator that rejects 75 of 132 correct
answers shipped with a green gate. Ranked, with more time I would:

1. **Re-run the whole binding derivation with an identity term in the gate.** One line —
   `r["identity_fail"] = sum(compare_template(t, g, g).outcome != "MATCH" for g in sols)` — and
   `bound` requires it to be zero. It costs 50 comparisons per template against 2,450 and it
   catches E-1, E-2, E-4 and E-6 at once. This is the highest-value change in the report.
2. **Re-audit the seven `UNBOUND` "N false accepts" reasons against §4.1's display tolerance**
   (E-12). My guess is that most of the small counts (1, 2, 2, 2, 4, 4) are boundary-correct
   matches and those templates should be bound; `null_to_null_bandwidth`'s 24 probably is not.
   Until that audit exists, `18 declined` is not a measurement.
3. **Replace the textual-identity truth predicate with a per-kind one.** For `numeric`, "the two
   golds' asserted numbers differ by more than gold's display tolerance". The current proxy is
   documented as audited and is not.
4. **Rebuild `nth_quantity` to return the i-th number and its own trailing unit**, not a window
   that can swallow the next one, and assert in `derive_kind` that part `i` selects number `i` on
   every derivation seed. That converts E-8 from a latent defect into a build-time failure.
5. **Make the unit a per-part declaration** (`units: tuple[str|None, ...]` on the multipart
   binding), and require every declared unit to have an entry in `UNIT_SURFACES`; a declared unit
   with no surface list should fail the derivation rather than reach a comparison. That kills
   `otherwise` and `e-05` at the source.

**Does the weakness exist outside this scope?** Yes, and predictably. The 40 unresolvable declared
units name whole domains: every thermodynamics template that answers in `kJ` or `K`, every
reaction-engineering template in `mol/dm^3` or `L/mol`, every queueing template in `minutes`, every
inventory template in `units`. Named concretely, the templates that will still refuse or reject a
correct answer after any narrow fix to `UNIT_SURFACES` are the composite ones:
`wave_equation_interpretation`, `wave_parameters_basic`, `time_to_phasor`, `system_properties`,
`undamped_natural_frequency_torsional` and `_translational`, `vdw_solve_for_volume`,
`coaxial_capacitance`, `autocorrelation_rect_pulse` — because their defect is E-4, one unit for
`n` parts, not the surface table. Outside `tests/comparators/`, the same declared-vs-inferred
confusion will reappear the moment class C's remaining templates are fitted, since D5.10's
derivation is the thing that will be reused.

**What later phases should do differently.** The phase inherited a correct instruction — *"a
comparator that resolves less than the one it replaced passes every accept-only gate"* — and built
a tool that prints `FALSE ACCEPTS 0` in capitals and `false rejects 56` in lower case on the same
line. Make the report symmetric: print `decided`, `false rejects`, and identity failures with the
same weight as false accepts, and make the exit code depend on all four. D5.7b already established
the principle for errors ("a template that raises has not passed"); it needs one more clause.

**What in the spec or contract is wrong or ambiguous.**

* **D4.1 §4.1 is now self-contradicting in practice.** It says the unit must be *declared, never
  inferred*, precisely so `7.65 litres` is not rejected. D5.10 derives the declaration from gold's
  own string, which is inference wearing a declaration's clothes, and `17.46 s` against
  `17.46 seconds` is now rejected. Either §4.1's rationale or D5.10's derivation has to go; they
  cannot both stand.
* **"Bound" is defined in `bindings.py` as zero false accepts and zero errors.** It should be
  defined as zero false accepts, zero errors, **and non-zero decided, and gold matching itself**.
  The current definition admits 11 templates that decide nothing.
* **§3's `mode` is specified (`all` vs `any`) and never used.** Every derived multipart binding is
  `mode="all"`; `numeric_parts` hard-codes it. `undamped_natural_frequency_torsional` — the spec's
  own worked example of `any` ("51.184 rad/s **or** 8.146 Hz") — is bound `all` with `n=3`. The
  scoring error §3 warns about ("read `any` as `all` and a model answering in rad/s alone is
  marked wrong") is shipped.
* **`phase4_comparators.md` §8's table (98.6% / 100%) disagrees with `score`'s current output
  (95.8% / 98.6%)** at this same ref. One is stale; E9's "unchanged from Phase 4" cites the
  number that is not in the document.
* **The brief's "26 bound multipart templates" is the inventory `answer_type` count.** The
  `multipart` *kind* is bound on **31** templates. Five more than the review scope assumed, and
  four of the six E-8 mis-slicers are in the gap.

---
---

# Round 2 — re-review at `4fd06c3`

## Verdict: PASS WITH FINDINGS — **the block lifts**

All seven round-1 findings re-derive as fixed, both of my colliding gold pairs now `MISMATCH`
*for the correct reason*, and I could not construct an over-acceptance anywhere at any depth I
could reach at 326,700 gold-pair comparisons. Two over-rejections remain — one new (`vector` notation), one inherited and still
untriaged (**E-12**, which the unit ablation has now made bite) — and neither is a false accept,
so neither blocks.

Working tree is identical to `4fd06c3` (`git diff 4fd06c3 HEAD` is empty). `extract.py`,
`kinds.py`, `answer.py` and `normalize.py` are unchanged from `7c04c95`; the fixes are in
`bindings.py`, `derive_bindings.py` and `cross_pair.py` only.

---

## R2.1 — The block question, answered first

### 1. Do the two colliding gold pairs now MATCH? **No — and the reason is now the right one.**

```bash
PYTHONIOENCODING=utf-8 python <scratchpad>/revE/p_final.py      # unchanged from round 1
```

```
template_hydrostatic_pressure_at_depth   seeds 157 / 237
  part1 slice was: 'The absolute pressure at a depth of 13.29 m is 101.347 kP'   -> selected 101.347
  part1 slice now: 'The absolute pressure at a depth of 13.29'                   -> selects   13.29
  verdict: MISMATCH | part(s) differ: q1:numeric      <- the DEPTH differs, and is now compared

template_rackett_equation_volume         seeds 105 / 206
  part1 slice was: '... liquid Carbon dioxide at 199.89 K is 36.22 cm^3'         -> selected  36.22
  part1 slice now: '... liquid Carbon dioxide at 199.89'                         -> selects  199.89
  verdict: MISMATCH | part(s) differ: q1:numeric      <- the TEMPERATURE differs, and is compared
```

In round 1 both returned `MISMATCH` because of the unit defect, while the quantity that actually
distinguished them was compared by nothing. Both now return `MISMATCH` on the distinguishing
quantity itself. **The masking is gone and nothing was traded for it.**

### 2. Is the slicer fixed generally, not just on those two?

`nth_quantity` now ends the slice at the number. Re-running my round-1 instrument — *does part `i`
select the `i`-th asserted number, over 12 seeds, on every multipart binding*:

```bash
PYTHONIOENCODING=utf-8 python <scratchpad>/revE/p_slice.py      # unchanged from round 1
# multipart bindings: 28
# multipart templates where at least one part does NOT select its own number: 0 of 28
```

**0 of 28**, against 6 of 31 in round 1. `p_const.py` (E-9) prints nothing: **no part of any bound
multipart binding reads a constant**, against 3 templates in round 1.

### 3. The collision sweep, re-run and deepened

`N=50` still cannot resolve a 1-in-31,000 collision, so I went past it: **gold×gold at `N=100`,
9,900 ordered pairs per template, over all 33 templates of the three structural kinds
(`multipart` 28, `vector` 3, `sequence` 2) — 326,700 ordered pairs.**

```bash
PYTHONIOENCODING=utf-8 python <scratchpad>/revE/r2_deep.py
```

It completed:

```
TOTAL {'pairs': 326700, 'MISMATCH': 326666, 'MATCH': 34}
```

The script prints a line per template only when that template produces a false accept, a false
reject or an error, and it printed **none**: so across all 33 structural-kind templates at `N=100`
there are **0 false accepts, 0 false rejects, 0 errors and 0 UNRESOLVED**. The 34 MATCHes are all
on textually-identical span pairs, which is why none was counted a false accept.

Broken out, the five formerly mis-slicing *bound* templates — the ones the question is actually
about:

```
template_hydrostatic_pressure_at_depth  {'n': 2}   {'pairs': 9900, 'MISMATCH': 9900}
template_particle_pathline              {'n': 3}   {'pairs': 9900, 'MISMATCH': 9900}
template_rackett_equation_volume        {'n': 2}   {'pairs': 9900, 'MISMATCH': 9900}
template_statically_indeterminate       {'n': 4}   {'pairs': 9900, 'MISMATCH': 9900}
template_time_to_phasor                 {'n': 4}   {'pairs': 9900, 'MISMATCH': 9900}
```

**49,500 ordered pairs, every one MISMATCH: 0 false accepts, 0 UNRESOLVED, 0 errors.** Two of these
carried round 1's only real collisions; every pair is now not merely decided but decided on the
quantity that distinguishes it. (The two named pairs sit at seeds 157/237 and 105/206, outside this
`N=100` window — they are verified separately and directly in R2.1 item 1.)

**And the round-1 collision search, re-run unchanged at this ref, now finds nothing.** This is the
instrument that produced the two collisions in the first place — 250 instances per template,
bucketed by the tuple of numbers the binding actually compares, looking for a bucket shared by
instances with different spans (`p_collide.py`, byte-identical to round 1):

```
template_hydrostatic_pressure_at_depth : 0 colliding compared-tuples over 250 instances   (was 1)
template_rackett_equation_volume       : 0 colliding compared-tuples over 250 instances   (was 1)
template_particle_pathline             : 0                                                (was 0)
template_statically_indeterminate      : 0                                                (was 0)
template_time_to_phasor                : 0                                                (was 0)
```

**Both collisions are gone, and they are gone for the right reason:** the compared tuple is no
longer degenerate, because part 1 now compares the depth and the temperature instead of comparing
the pressure and the volume twice. `250 x 249 = 62,250` ordered pairs per template are covered by
the bucketing. **The coordinator's question 2 is answered: no, the fixed slicer admits no collision
my instruments can find, at 5x the depth the gate runs at.**

So the beyond-`N=50` question is **closed**: 326,700 ordered pairs at twice the gate's instance
count over every structural-kind binding, plus 62,250 collision-bucketed pairs per template at five
times it on the five that motivated the question, and not one false accept.

### 4. The identity term, re-derived at full validation depth

I did not read R2 off the tool that asserts it. 118 templates × all 50 validation seeds:

```bash
PYTHONIOENCODING=utf-8 python -c "
import sys, collections; sys.path.insert(0,'.')
from tests.template_integrity.core import discover, generate
from tests.comparators.bindings import BINDINGS, compare_template
refs={r.template_id:r for r in discover()}
c=collections.Counter()
for tid in sorted(BINDINGS):
    for s in range(50):
        g=generate(refs[tid],s,capture=False).solution
        c[compare_template(tid,g,g).outcome]+=1
print(c)"
# Counter({'MATCH': 5900})
```

**5,900 of 5,900.** Round 1's figure was 57 of 132 templates matching themselves; it is now 118 of
118 on every seed. E-1 is closed, and it closed E-2, E-3, E-4 and E-6 with it, exactly as the
round-1 §5 predicted it would.

---

## R2.2 — Claims R1–R5, re-derived

| # | claim | re-derived | agreement |
|---|---|---|---|
| **R1** | 118 of 150 bound, 32 named unbound, each with a measured reason | `len(BINDINGS)==118`, `len(UNBOUND)==32`, sum 150, disjoint; every `UNBOUND` value is a measurement; `cross_pair` prints `118 scored / 32 skipped` | **agrees** |
| **R2** | 118/118 MATCH a verbatim copy of gold on all 50 validation seeds | my own sweep above → `Counter({'MATCH': 5900})`; no bound template has a non-zero `identity_failures` in `VALIDATION` | **agrees, independently derived** |
| **R3** | gold×gold 15,576 / 118 / 32; all four terms zero; UNRESOLVED 0 | `python -m tests.comparators.cross_pair` → `pairs 15576 / templates-scored 118 / templates-skipped 32`; `MATCH 188 MISMATCH 15388 UNRESOLVED 0`; `ERRORS 0`; `FALSE ACCEPTS 0`; `FALSE REJECTS 0`; `IDENTITY 0 of 1416`; `GATE ... PASS`. Every kind now decides **100.0%** | **agrees, exactly** |
| **R4** | archive×gold 18,115 (1,273 + 16,842) / 57 / 27; 0 XI; 0 errors; decided 67.8%, up from 56.8% | same run → `18115 (1273 positive + 16842 negative) / 57 / 27`; `cross-instance accepts 0 / 16842`; `ERRORS positive 0 negative 0`; `863/1273 (67.8%)` | **arithmetic agrees; "up from 56.8%" needs reading — R2-F2** |
| **R5** | D4.4 95.8% / 98.6%; archive 100% / 97.8%; suites exit 0 | `python -m tests.comparators.score` → `adversarial precision 95.8% recall 98.6% PASS`; `archive precision 100.0% recall 97.8% PASS false accepts: 0`; `narrative 20/20 UNRESOLVED` | **agrees, exactly** |

Two structural checks on the new gates, because a gate that never fires is not a gate:

* **`MIN_DECIDED_RATE`** — no bound template has `decided == 0`, and the **minimum decided rate
  across all 118 is 1.0**. The bimodality the docstring claims is real at this ref.
* **The identity and constant-part gates fired, and their casualties are named in `UNBOUND`** —
  9 templates carry `rejects a verbatim copy of gold on 50/50 seeds` (8 `symbolic` plus
  `decimation_aliasing_analysis`), 3 carry
  `part(s) [...] read a constant on every instance`. Not decorative.

---

## R2.3 — Findings

### R2-F1 — CONFIRMED — over-rejection — the `vector` comparator accepts **0 of 83** archived answers, and at least 14 of the refusals are notation alone

`vector`'s 12.0% archive decided rate is declared a stated trade. It is not a trade; it is a
missing surface declaration, and it is the one place round 1's E-2 lesson was not applied.

`vector_components` matches `_HAT_RE = [-+]?\s*[\d.eE+-]+\s*[a-z]_hat` — ASCII `x_hat`, which is
what **gold** writes. Every archived model answer writes the Unicode hat:

```
GOLD: 'F_m = (1.791e-06 x_hat + -1.056e-05 y_hat + 7.810e-05 z_hat) N.'
CAND: 'F_m = (1.79 x^ - 10.56 y^ + 78.10 z^) uN, |F_m| ~= 78.83 uN'      [x-hat etc. are U+0302]
      -> UNRESOLVED | no vector in the candidate answer span

GOLD: 'F_m = (-3.039e-04 x_hat + -2.831e-04 y_hat + 4.706e-04 z_hat) N.'
CAND: 'F_m = (-3.039 * 10-^4 x^ - 2.831 * 10-^4 y^ + 4.706 * 10-^4 z^) N'
      -> UNRESOLVED | no vector in the candidate answer span
```

The second candidate states gold's three components to gold's own precision and is refused for
spelling the unit vector with a combining circumflex rather than `_hat`. Rewriting **only the
notation** and touching no comparator logic:

```bash
PYTHONIOENCODING=utf-8 python <scratchpad>/revE/r2_vec2.py
# as shipped          : {'UNRESOLVED': 73, 'MISMATCH': 10}
# with x_hat notation : {'UNRESOLVED': 59, 'MISMATCH': 24}   <- notation only, no logic changed
```

14 refusals become decided verdicts from a character substitution. The residue — still **0 MATCH**
— is a second layer: `uN` against `N`, and `10-^4` written with U+2212. So across all three
`vector`-bound templates the comparator **accepts none of the 83 real model answers**, and those
three templates contribute zero recall to the benchmark while counting as bound.

Worth recording for the design: this is invisible to the new identity gate *by construction*, since
gold matches gold because both sides use `x_hat`. **The identity term catches gold-side defects and
structurally cannot catch candidate-side normalisation gaps.** Only archive text sees those.

**Answering the coordinator's question directly: for `vector`, "declared as a stated trade" is the
wrong call.** D4.1 §4.2 already solved exactly this problem for `categorical` by declaring the
surfaces of a label. The unit-vector form needs the same declaration; it is a handful of
characters, not a design decision, and calling it a trade puts a normalisation bug on the
residual-risk register where it will be read as an item property.

### R2-F2 — CONFIRMED — the E-12 truth predicate is now removing working templates from the bound set

Round 1's **E-12** was filed CONFIRMED and is not on the actioned list. The unit ablation has made
it bite. Three templates are *newly* unbound at this ref with the reason
`N false accepts in 2450 pairs`, and all three are the display-tolerance artefacts E-12 named:

```
template_batch_reactor_second_order   5 "false accepts" in 2450 pairs
    gold 'The required reaction time is 5.0 seconds.'       cand '... is 4.99 seconds.'
template_pfr_volume_changing_rate     2 "false accepts" in 2450 pairs
    gold 'The required PFR volume is 21.9 liters.'          cand '... is 21.92 liters.'
template_vdw_solve_for_pressure       2 "false accepts" in 2450 pairs
    gold 'The pressure exerted by the Oxygen is 14.9 bar.'  cand '... the Ethylene is 14.85 bar.'
```

Gold displays one decimal, so D4.1 §4.1's tolerance is `0.5 x 10^-1 = 0.05`, boundary inclusive.
`|4.99 - 5.0| = 0.01` and `|21.92 - 21.9| = 0.02` are **correct MATCHes under the specification.**
The truth predicate — *a MATCH is a false accept iff the two spans differ textually* — calls them
false accepts, and three working bindings were discarded for obeying §4.1.

Cost, measured, by differencing the per-kind archive tables at `7c04c95` and `4fd06c3`:

| kind | round 1 (rows × decided) | round 2 (rows × decided) | decided rows |
|---|---|---|---:|
| `multipart` | 713 × 43.8% | 668 × 52.7% | 312 → **352**  (+40, the E-4 fix) |
| `numeric` | 488 × 95.5% | 452 × 96.5% | 466 → **436**  (−30, these three templates) |
| `symbolic` | 155 × 3.9% | 6 × 100% | 6 → 6 (unchanged) |
| `vector` / `sequence` / `categorical`* | unchanged | unchanged | 84 → 84 |
| **total** | **1,503 × 56.8%** | **1,273 × 67.8%** | **853 → 863** |

So R4's `56.8% → 67.8%` is arithmetically right and is mostly a denominator: 230 archive rows left
the corpus along with the templates that could not decide them, and only **+10 rows** actually
became decided — `+40` won by the unit ablation, `−30` handed straight back to E-12. Archive
coverage fell from 70 templates scored to 57, and bound coverage from 132 to 118. The comparator is
now *honest*, which is the whole point of the round; but "decided rate up 11 points" and "the
comparator decides more" are different claims and only the first is supported by this run.

`vdw_solve_for_pressure` is the interesting one and should **not** simply be re-bound: gold really
is Oxygen at 14.9 bar and Ethylene at 14.85 bar, and at one displayed decimal the *item* cannot
distinguish two substances. That is D-050's shape — a property of the item, not of the comparator —
and it belongs on the item-design list rather than in `UNBOUND` under a label that says the
comparator did something wrong.

Seven further `UNBOUND` entries carry the same untriaged reason, unchanged since round 1:
`decimation_aliasing_analysis` 2, `floating_object_submersion_depth` 4, `gauss_law_symmetric` 4,
`null_to_null_bandwidth` 24, `pitzer_correlation_z` 1, `signal_energy_power` 2,
`truss_method_of_sections` 2. **10 of the 32 unbound templates now rest on a predicate known to be
wrong in this direction, and not one has been audited against §4.1.**

### R2-F3 — CONFIRMED — for `symbolic`, "bind fewer and name the rest" *is* laundering; but the gate is not what is at fault

The coordinator asked. The answer is yes for `symbolic`, and no for the gate.

All 8 unbound `symbolic` templates fail on gold against a verbatim copy of itself, and every cause
is a comparator defect rather than a property of the item:

```bash
PYTHONIOENCODING=utf-8 python <scratchpad>/revE/r2_sym.py
```

```
undamped_response_initial_conditions   'x(t) = -0.006*cos(34.6384*t) (m)'
    UNRESOLVED | symbol(s) outside the declared alphabet: ['c', 'm', 'o', 's']
phasor_addition    'v_total(t) = 38.93 * cos(420*t + 7.48 deg).'
    UNRESOLVED | symbol(s) outside the declared alphabet: ['c', 'd', 'e', 'g', 'o', 's']
ft_esd_rect_pulse  'G(f) = 48.0 * sinc(4.0*f)'
    UNRESOLVED | symbol(s) outside the declared alphabet: ['c', 'i', 's']
ber_estimation_mary / bpsk_energy_basis / cd_dc_system_analysis / standing_wave_formation
    UNRESOLVED | could not parse an expression: SympifyError
```

`['c','o','s']` are the letters of `cos`; `['c','i','s']` the letters of `sinc` less the declared
`n`; `['d','e','g']` the letters of `deg`. **The parser is splitting function names into free
symbols, and the alphabet check is then rejecting gold's own answer.** D4.1 §4.4 says the eight
non-fragment templates *"carry sinc, exp and Q-functions, so the fragment is not enough in general.
`sympy` is used there"* — `sympy 1.14.0` is installed and the code does not in fact parse a
function call. And `derive_kind` routes to `symbolic` on
`_FUNC_RE = \b(?:sinc|cos|sin|exp|log|Q|tanh|sqrt)\s*[\^(]`: it selects for precisely the construct
the comparator cannot read.

So **`symbolic` at 1 of 9 is not a measurement of the corpus, it is a measurement of the
comparator**, and the `UNBOUND` reason string says the symptom rather than the cause. Any later
phase quoting "8 of 9 symbolic templates are not bindable" will be quoting a parser bug as a corpus
property. It reads like an evening's work, not a Phase 6 design item.

**The gate itself is doing exactly the right thing.** At `7c04c95` these same 8 templates were
counted as *bound*, inside the "132" and the "340,550 validated pairs", at 0% decided. The new gate
is what surfaced them, and the finding exists only because the fix worked.

---

## R2.4 — Falsification attempts that failed

1. **The two colliding gold pairs** (the coordinator's question 1). Both MISMATCH, on the right
   quantity. I could not turn either into a false accept.
2. **Deep gold×gold at `N=100` on all 33 structural-kind templates**, 326,700 ordered pairs.
   `TOTAL {'pairs': 326700, 'MISMATCH': 326666, 'MATCH': 34}` with no per-template line printed:
**0 false accepts, 0 false rejects, 0 errors, 0 UNRESOLVED**, the 34 matches all on textually
identical spans. And round 1's own `N=250` collision search, re-run unchanged, now returns
**0 colliding compared-tuples on all five** of the formerly mis-slicing templates, against the 2 it
found before — including on both templates that carried a real collision.
3. **`nth_quantity` selection audit** — 28 multipart bindings × 12 seeds × every part: **0** parts
   select another part's number, against 6 templates in round 1.
4. **Constant-part audit** — **0**, against 3 in round 1 (`autocorrelation_rect_pulse`,
   `coaxial_capacitance`, `vdw_solve_for_volume`), all three now correctly in `UNBOUND` naming that
   reason.
5. **Identity at full validation depth** — 5,900 of 5,900 MATCH, no exceptions, no errors.
6. **archive×gold cross-instance accepts** — **0 of 16,842** here, and 0 of 19,362 in round 1
   including with the unit check ablated by hand. Two different binding tables, two different unit
   policies, no over-acceptance on real model text either time. This is now the best-evidenced
   negative in the phase.
7. **The decided-rate floor.** I looked for a bound template sitting just above `MIN_DECIDED_RATE`
   and there is none — the minimum among the 118 is 1.0, so `<= 0.0` is loose but unexercised. Not
   a finding at this ref. It becomes one the moment a partly-deciding binding appears, because
   `<= 0.0` will admit a binding that decides two pairs in 2,450.
8. **R1–R5 as arithmetic** — every number re-derived from the code; every one agrees.

---

## R2.5 — Further probing

Ranked, one line as instructed: **(1)** declare the `x̂` / `μN` surfaces for `vector` — R2-F1, three
templates and 83 archived answers for a few characters; **(2)** audit the 10 `UNBOUND` entries that
rest on the textual-identity predicate against §4.1's display tolerance — R2-F2; most should
probably be bound, and `vdw_solve_for_pressure` should go to item design rather than back to
`UNBOUND`; **(3)** teach the `symbolic` parser to recognise a function call before any later phase
quotes `symbolic 1 of 9` as a property of the corpus — R2-F3.
