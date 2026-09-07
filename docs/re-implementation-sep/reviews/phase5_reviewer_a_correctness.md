# Phase 5, Track A — Reviewer A (correctness)

**Ref reviewed:** `e3ee9ab` (frozen), against `1bd3e4f` (`master`, worktree `C:/wtm`).
**Reviewer:** A — correctness. Independent; the implementer's reasoning was not seen.
**Protocol:** R1 (claims under test), R6 (one mandatory gate task, time-boxed).
All numbers below were re-derived with the reviewer's own code unless a row says
otherwise. Reviewer scripts live in the session scratchpad
(`.../scratchpad/revA/{probe,deep11,dump,diffdump}.py`); nothing in the repository
was modified except this file.

---

## 1. Verdict

**PASS WITH FINDINGS** — the corpus really is clean on all four defect classes
(re-derived independently, 60,000 instances corpus-wide and 33,000 on the eleven),
but **the gate that says so does not gate three of the four classes against the
regression shapes that actually threaten them**: seven planted defects were caught
by no instrument at all.

---

## 2. Independent re-derivation

Environment for every command: `export PYTHONIOENCODING=utf-8`, cwd
`C:/Users/ayesha.gull01/EngTrace` unless the command names `C:/wtm`.

### C1 — T4 passes 150/150 — **AGREE**

```bash
python -m tests.template_integrity.run --checks T4                 # 0 failing
python -m tests.template_integrity.run --checks T4 --seeds 400     # 0 failing
python -m tests.template_integrity.run --checks all --json all_after.json
```
`--checks all` at `e3ee9ab`: **T1 29, T2 0, T3 0, T4 0, T5 66, T6 142, T7 83**.
At `1bd3e4f`: identical except **T4 3**. T4 is green at the default 25 seeds and
still green at 400 (60,000 instances). The three T4 flips are exactly
`cd_dc_system_analysis`, `euclidean_distance_binary`, `finite_convolution`.

### C2 — corpus clean on all four classes, 0 outside Track A, 0 generation errors — **AGREE on the corpus, DIVERGE on the framing of class 3**

My own four detectors (written from the spec's defect table, not from
`phase5_contract_scan`), swept over **150 templates × 400 seeds**:

| class | my probe, before (`C:/wtm`) | my probe, after (`e3ee9ab`) |
|---|---|---|
| 1 malformed `**Step N:**` | 3 templates — exactly the three named | **0** |
| 2 non-canonical answer marker | 5 templates — exactly the five named | **0** |
| 3 doubled sign, whole solution | **16** templates | **14** (the two Track A templates fixed) |
| 3 doubled sign, answer span | **4** templates | **2** |
| 4 zero-coefficient product | 4 templates (incl. `decimation_aliasing_analysis`) | 3 (Track A's gone) |
| generation errors | 0 | 0 |

Also swept the eleven at **3,000 seeds each (33,000 instances)**: all four classes
zero, and every span recovers `**Answer:**`.

```bash
NSEEDS=400 PROBE_OUT=after400.json  python scratchpad/revA/probe.py C:/Users/ayesha.gull01/EngTrace
NSEEDS=400 PROBE_OUT=before400.json python scratchpad/revA/probe.py C:/wtm
python scratchpad/revA/deep11.py 3000
```

Divergence: the scan prints `complex_sign 0 templates`, which a reader will pair
with C4's "corpus-wide this class is 16 templates". Those two sentences are about
different predicates — the scan's detector requires a literal `j` — and only one
of them is the class D-061 defines. See **finding 1**.

### C3 — `0*pi` fires on 39 of 2,000 before, 0 after — **AGREE, exactly**

```bash
python scratchpad/revA/c3check.py C:/wtm     # 39 / 2000 = 1.95 %, 0 errors
python scratchpad/revA/c3check.py .          #  0 / 2000 = 0.00 %, 0 errors
```
My detector was `(?<![\w.])0(?:\.0+)?\s*\*\s*pi` over question+solution, counting
**instances**, not fragments — the fragment count is 22 per 400 seeds because most
hits print twice, and a fragment-based reading of this number would not reproduce.

### C4 — doubled sign: 16 corpus-wide, 4 in the answer span, 2 in Track A — **AGREE, exactly**

Reproduced by the C2 sweep above, with my own regex
`(?<![\w.\-])[+−-][ \t]*(?:j[ \t]*)?[+−-][ \t]*(?=[.\d]|j[.\d]|[A-Za-z_])`.
Before: 16 in the solution, 4 in the answer span
(`lorentz_force`, `phasor_addition`, `time_to_phasor`,
`continuous_to_discrete_conversion`), 2 of those in Track A. After: 14 and 2; the
two Track A templates are clean at 3,000 seeds. The residual 14 are Phase 6's by
the spec, and are named in §5.

### C5 — `phasor_addition` T1 9 → 5 over 25 seeds; `sqrt(-41.2^2 + 86.45^2)` printed 95.77, evaluates 76.00 — **AGREE, exactly**

```bash
python -m tests.template_integrity.run --templates template_phasor_addition --checks T1 --json pa_after.json
cd C:/wtm && python -m tests.template_integrity.run --templates template_phasor_addition --checks T1 --json pa_before.json
```
Before: `failures 9`, seeds `[2,5,7,9,12,13,15,17,21]`, `checks_performed 142`,
`coverage 0.27`, and the reported example is verbatim
`A_total = sqrt(-41.2^2 + 86.45^2) = 95.77` / `evaluated ... = 76.0010690714282`.
After: `failures 5`, seeds `[4,5,7,9,16]`, `checks_performed 152`, `coverage 0.29`.

Worth recording beside the headline: the failing **seed set is not a subset** —
4 and 16 are *newly* failing. They fail because the parenthesisation made ten more
lines parseable, not because anything got worse; both are ordinary
rounding-boundary misses (δ ≈ 0.0058 against a 0.005 tolerance). The item-pool
note reports the coverage rise but not the two newly-failing seeds, and a later
reader diffing seed lists will trip over that.

### C6 — 1,000 of 2,200 spans carried debris before, 0 after; corpus-wide 1,200/1,200 — **AGREE, exactly, and strengthened**

```bash
python scratchpad/revA/dbg.py C:/wtm    # 2200 spans, 1000 with debris/non-canonical
python scratchpad/revA/dbg.py .         # 2200 spans, 0
```
Before, the debris is on exactly the five marker templates at a 100% rate;
`answer_span` recovers the marker `'Final Answer'` and the span begins `'**\n…'`
(or `'s:**\na) CSTR Volume…'` on `levenspiel`). After: zero.

Corpus-wide I ran a stronger version than the claim: **150 templates × 400 seeds =
60,000 spans, 60,000 recover `**Answer:**`, zero empty spans** (the claim's
1,200 is 150 × 8). Every template in the corpus emits **exactly one**
`**Answer:**` — my class-2 probe flags both `0` and `>1` and found neither.

### C7 — 8 formatting-only / 3 text-changing; 2,560 questions, 2,045 answer bodies, 0 errors — **AGREE, exactly, every cell**

Re-derived with my own dumper and differ (two separate processes, root asserted
against `core.REPO_ROOT` so the in-process-reload trap cannot fire):

```bash
python scratchpad/revA/dump.py C:/wtm mine_before.json 2000
python scratchpad/revA/dump.py C:/Users/ayesha.gull01/EngTrace mine_after.json 2000
python scratchpad/revA/diffdump.py
```
Totals: solutions **18,308**, questions **2,560**, answer bodies **2,045**,
generation errors **0 / 0**. Per template every figure matches the note's table,
including `euclidean_distance_binary` 996/2,000 and `decimation` 39/2,000. My
answer-body definition (text after the last marker in
`{**Final Answers:**, **Final Answer**, **Answer:**, **Answer**}`) was chosen
before reading the note's and reproduces it.

### C8 — only T4 moved per template — **AGREE at verdict level; overstated as "unchanged"**

Diffing the two `--checks all` JSON reports per template and per field:

| check | verdict flips | numeric movement |
|---|---|---|
| T1 | 0 | `phasor_addition`: failures 9→5, marginals 57→64, checks 142→152, coverage 0.27→0.29 |
| T2, T3, T6, T7 | 0 | none |
| T4 | **3**, all fail→pass | — |
| T5 | 0 | `phasor_addition` `operand_restatements` 0→5; `time_to_phasor` 0→1 |

The T1 movement is documented in the item-pool note. The **T5 movement is not
recorded anywhere** — it is non-gating (T5a reports operand restatements
separately) and is caused by the new `_signed_term(...)` calls that remain inline
in f-strings, but "unchanged per template" is not what the measurement says.

### C9 — T6 `breaches` nondeterministic in **order** only — **AGREE on the substance; the census number does not reproduce**

Five T6 runs at `e3ee9ab`:

```bash
for i in 1 2 3 4 5; do python -m tests.template_integrity.run --checks T6 --json t6_$i.json --quiet; done
```
**Order-insensitively: 0 templates differ across 5 runs** — the substantive claim
holds. Order-sensitively I get **9 templates over 3 runs and 10 over 5**, not 7.
The figure is a function of how many runs you compare (it is a hash-ordering
sample), so "7 of 150" is not a property and should not be written as one.

---

## 3. Findings

All findings are about the **gate**, not the corpus. The corpus is clean; I could
not break that. Every finding below is a way the four classes could come back
without any instrument going red.

Reproduction for findings 1–5 is one script, which plants each defect into real
generated output of a clean host (`template_bpsk_energy_basis`, seed 0) and runs
my probes, `phase5_contract_scan.scan_solution`, `t4_contract.run` and
`t1_closure.run` over the result:

```bash
python scratchpad/revA/plant.py
```

```
plant                            | mine       | scan             | T4   | T1
BASELINE (unplanted host)        | -          | -                | PASS | PASS
1a malformed **Step 2: **        | C1         | step_marker      | FAIL | PASS
1b title inside bold             | C1         | step_marker      | FAIL | PASS
1c LAST step unbolded            | C1         | -                | PASS | PASS   <-- F5
1d step 2 heading deleted        | C1n        | -                | FAIL | PASS
2a **Final Answer** swap         | C2         | answer_marker    | FAIL | PASS
2b extra ## Final Answer         | C2         | -                | PASS | PASS   <-- F6
2c lowercase **answer:**         | C2         | answer_marker    | FAIL | PASS
3a j-form in solution            | C3sol      | complex_sign     | PASS | PASS
3b NON-j doubled sign in soln    | C3sol      | -                | PASS | PASS   <-- F1
3c NON-j doubled sign in Q       | C3q        | -                | PASS | PASS   <-- F1+F2
3d j-form in QUESTION only       | C3q        | -                | PASS | PASS   <-- F2
4a 0*pi in answer span           | C4         | degenerate_prod  | PASS | PASS
4b 0.0*pi in answer span         | C4         | -                | PASS | PASS   <-- F3
4c 0*pi in QUESTION only         | -          | -                | PASS | PASS
```

### F1 — `complex_sign` gates only the `j`-form, not the class — **CONFIRMED**

`_COMPLEX_SIGN = r'[-+]\s*j\s*[-+]\s*\d'` requires a literal `j`. The other half of
the same defect — `876*t + -86.8 deg` — is **not detected**. This is not a
hypothetical shape: `waves_and_phasors.py`'s own `_signed_term` docstring names
`"876*t + -86.8 deg"` as the first of the two things the hard-coded `+` printed,
and `_signed_term` is applied in both templates precisely to fix it. So the gate
does not cover the half of Track A's class-3 fix that changed **2,560 question
strings**.

*Reproduction:* plants `3b` / `3c` above; or
`python -c "import re;print(re.search(r'[-+]\s*j\s*[-+]\s*\d','876*t + -86.8 deg'))"` → `None`.

*Impact:* **medium.** A revert of `_signed_term(phi_value, phi_unit)` in
`template_time_to_phasor`'s `time_domain_expr` restores the defect on ~51% of
instances and turns nothing red. It also makes the scan's `complex_sign 0
templates (0 outside Track A)` line read as evidence for a corpus property that my
own sweep shows is false for 14 templates.

### F2 — the scan never reads the question — **CONFIRMED**

`scan_template()` calls `scan_solution(inst.solution, res)` and nothing else.
`inst.question` is never scanned by any of the four detectors, and T4 only checks
that the question is non-empty. Yet Track A's class-3 fix changed **question**
text on 2,560 of 22,000 instances, and two of the five residual doubled-sign
question emitters (`continuous_to_discrete_conversion`, `nyquist_rate_determination`)
are the exact templates Phase 6 will edit next.

*Reproduction:* plants `3c`, `3d`, `4c` above — all three pass every instrument.

*Impact:* **medium.** The output contract is a property of the emitted *item*, not
of the solution string. As written, a template can emit `V = 30.5 + j-51.22` or
`omega_a = 0*pi` in its question forever.

### F3 — `degenerate_product` misses the float-zero form, and its own census is short by one template — **CONFIRMED**

`_DEGENERATE_PRODUCT = r'(?<![\w.])0\s*\*\s*(?=pi\b|[A-Za-z]_?\w*\b)'` cannot match
`0.0*pi`: the second `0` is preceded by `.` (lookbehind rejects) and the first is
followed by `.` (not `*`). Live evidence, not just a plant:
`template_undamped_response_initial_conditions` emits
`x(t) = 0.0 * cos(21.7438 * t) + (-0.2116) * sin(...)` on **142 of 400 instances**,
and the scan's census reports **2** templates where my probe finds **3**.

*Reproduction:*
```bash
python scratchpad/revA/plant.py                 # row "4b 0.0*pi in answer span" -> scan silent
python -c "import re;print(re.search(r'(?<![\w.])0\s*\*\s*(?=pi\b)','x = 0.0*pi'))"   # None
```

*Impact:* **low-medium.** In `undamped_response_initial_conditions` the zero term is
in a derivation step and the answer span correctly drops it, so it is census
material, not a gate failure — but the census silently lost it, and a gate that
misses `0.0*pi` in an answer span misses the class it exists for. `round(x, n)` of
a small float is the ordinary way a template arrives at `0.0`, so this is the
*likelier* regression shape, not the rarer one.

### F4 — classes 3 and 4 have no standing gate at all — **CONFIRMED**

`phase5_contract_scan` is referenced by exactly three files
(`template_redesign_spec.md`, `DECISIONS.md`, and itself). It is **not** in
`run.py`'s `ALL_CHECKS = ('T1','T2','T3','T4','T5','T6','T7')`, it is not imported
by `run.py`, and the repository has no CI (`.github/`, Makefile, tox and pytest
config are all absent). So the only standing corpus gate is
`python -m tests.template_integrity.run`, and it sees classes 1 and 2 only.

*Reproduction:*
```bash
grep -rn "phase5_contract_scan" --include=*.py .        # only the module itself
python - <<'X'
from tests.template_integrity.run import ALL_CHECKS; print(ALL_CHECKS)
X
```

*Impact:* **medium-high, and it is the phase's exit-gate item.** The exit gate asks
for "acceptance evidence … for the two T4 cannot see", and the evidence exists and
is good (`--selftest` is a real planted-defect test and it passes). But evidence is
not a gate. After this phase merges, the two new detectors run only when a human
remembers a phase-numbered module name.

### F5 — an unbolded final step heading passes everything — **CONFIRMED**

Both `t4_contract._STEP_ANY` (`\*\*\s*Step\s*(\d+)\s*[:\.]?`) and
`phase5_contract_scan._STEP_ANY` (`\*\*\s*Step\b…`) require the leading `**`, so a
heading that loses its bold entirely is invisible to the "however malformed"
detector. T4 then compares the *surviving* strict markers for contiguity, which
succeeds whenever the lost heading is the **last** one (`1,2,3` is contiguous
whether or not a step 4 exists).

*Reproduction:* plant `1c`. Delete a *middle* heading instead (plant `1d`) and T4
correctly fails, which is what makes the last-heading case a genuine hole rather
than a general blindness.

*Impact:* **low-medium.** This is the same failure mode the phase was chartered
against — `cd_dc_system_analysis` was urgent because "Step 3 computes the answer"
— and the final step is exactly where the answer is computed.

### F6 — a non-canonical marker *added beside* the canonical one passes T4 and the scan — **CONFIRMED**

`t4_contract.check_instance` flags `non_canonical_answer` only when
`CANONICAL_ANSWER_MARKER not in sol`; `scan_solution`'s `_ANSWER_ANY` requires
`\*\*…\*\*` and cannot see a `##`-heading marker at all. A solution containing both
`## Final Answer` and `**Answer:**` therefore passes both. It does **not** pass the
candidate side unchanged: `normalize.ANSWER_MARKERS` puts `##\s*Final\s+Answer\s*`
at **priority 0**, above `\*\*Answer:?\*\*`, so `answer_span` would extract from the
heading, not from the canonical marker — the exact disagreement D5.3 exists to
prevent.

*Reproduction:* plant `2b`.

*Impact:* **low today** (zero templates do it), **structural going forward** — it is
the one marker shape on which the gold-side and candidate-side lists actively
disagree about priority.

### F7 — `levenspiel_plot_interpretation`'s answer span swallows a 300-character `**Note:**` block — **CONFIRMED**

The marker itself is clean (`answer_span` returns exactly `**Answer:**`, so C6 is
untouched), but the *span* is:

```
a) CSTR Volume = 13.53 L
b) PFR Volume = 9.23 L

**Note:** The PFR requires less volume than the CSTR (9.23 L vs 13.53 L), a 31.8%
volume reduction. For a rate that falls with conversion …
```

Across the corpus this is the **only** template whose answer span carries trailing
commentary (60 seeds × 150 templates; its span is 363 chars against a next-longest
of 199), and it does so on 100% of instances. The block re-states both answer
values and introduces a third number (`31.8%`) that is not an answer.

*Reproduction:*
```bash
python scratchpad/revA/span_adv.py     # levenspiel: 60/60 trail hits, max span 363
```

*Impact:* **low for Track A, real for Track B.** It is pre-existing (`C:/wtm` has
the same Note under `**Final Answers:**`), so it is not a regression. But Track A
is where the answer-span contract was decided, N1 says the scalar extractor reads
the *first* number in the span and D5.6 will score a *last-number* rule against
this corpus — and on this template the last number in the span is `31.8`. Related:
`normalize._QUALIFIER_RE` bounds its discourse-marker lookback at 80 characters, so
a marker word appearing later inside a Note this long would not be recognised as
qualifying.

### F8 — C9's "7 of 150" is not a reproducible number — **CONFIRMED (documentation)**

9 templates over 3 runs, 10 over 5 (§2, C9). *Impact:* **low**; the conclusion the
number supports is correct and independently reproduces.

### F9 — C8's "unchanged per template" is a verdict claim, not a measurement claim — **CONFIRMED (documentation)**

T5 `operand_restatements` moved 0→5 and 0→1 and is recorded nowhere (§2, C8).
*Impact:* **low**; non-gating, but it is the kind of unrecorded movement the
item-pool note exists to catch.

---

## 4. Falsification attempts that failed

Things I tried in order to break this work, which held:

1. **Rewrote all four detectors from scratch** from the spec's defect table and
   swept 150 × 400 = 60,000 instances before/after, plus the eleven at 3,000 seeds
   (33,000 instances). My detectors are *broader* than the implementer's on classes
   3 and 4 and they still find **zero** hits inside Track A's scope after the fix.
   The corpus really is clean.
2. **Tried to reproduce C4's 16/4/2 with a deliberately different regex** (any
   doubled sign, `j` optional, unicode minus allowed, markdown rules excluded) —
   got 16 templates in the solution, 4 in the answer span, 2 in Track A. Exact
   agreement with D-061 on a number I expected to diverge on.
3. **Tried to catch the seed-count trap** — ran T4 at 400 seeds as well as the
   default 25 (60,000 instances) and the eleven at 3,000. No rare-instance defect
   hides at a higher seed count.
4. **Tried to make the item-pool numbers move** with my own dumper, my own
   answer-body definition, and an assertion that `core.REPO_ROOT` equals the tree
   root passed in (so the in-process-reload trap could not silently make both
   sides identical). Reproduced 18,308 / 2,560 / 2,045 / 0 errors, cell for cell.
5. **Tried to prove the T4 severity change is a no-op** — it is load-bearing. On
   `C:/wtm`'s T4 all three of `**Final Answer**`, `**Final Answers:**` and
   `**Answer**` planted into a clean host give `passed=True` with the defect
   *printed*; on `e3ee9ab` all three give `passed=False`. This is the one thing a
   green corpus genuinely cannot show, and it holds.
6. **Tried to break `answer_span` corpus-wide**: 60,000 spans, zero empty, zero
   over-peeled, zero non-canonical markers, and exactly one `**Answer:**` per
   solution on every template. The only surprising span in the corpus is
   `levenspiel`'s (F7), and even there the marker is correct.
7. **Tried `_signed_term`'s sign edge cases** — `-0.0` (`value < 0` is `False`, so
   `+ 0.0`, never `- 0.0` or `j-0.0`), exact `±180.0 deg`, and `round()` collapsing
   a small negative imaginary part to `-0.0`. 3,000 seeds on both phasor templates
   produce no `j-`, no `+ -`, no `- -`.
8. **Tried to find a T-check regression hiding behind the T4 improvement** — full
   `--checks all` on both trees, diffed per template *and per field*, not just per
   verdict. Only the movements in §2/C8 exist, and both are explained.
9. **Tried to plant class 4 in a derivation step** to see whether the answer-span
   scoping is over-tight — the detector correctly stays silent and the census
   correctly records it. The scoping decision is sound; only its regex is not (F3).
10. **Checked the two "cannot see" templates against T1** — planting classes 3 and
    4 does not make T1 fail either, so there is no accidental second gate. This is
    what makes F1–F4 findings rather than curiosities.

---

## 5. Further probing and improvements

Ranked. Each is a concrete change, not a direction.

**1. Wire the two new detectors into `run.py` as **T8**, or fold them into T4.**
(Addresses F4.) `phase5_contract_scan` is a phase-numbered module outside
`ALL_CHECKS`; six months from now nobody runs it. Concretely: add
`checks/t8_emission.py` carrying `complex_sign` and `degenerate_product`, add
`'T8'` to `ALL_CHECKS`, and keep `phase5_contract_scan` as the *report* that cites
it. Cost is under an hour and it converts one-time evidence into a standing gate.
The seed-count problem is real (class 4 needs ~400) — run T8 at 400 while the rest
of the suite stays at 25, rather than lifting the global default.

**2. Widen `_COMPLEX_SIGN` to the class D-061 actually defines, then re-measure its
false-positive rate.** (Addresses F1.) The class is *a doubled sign*, not *a
doubled sign next to a `j`*. Suggested replacement, which is what I swept with:
`(?<![\w.\-])[+−-][ \t]*(?:j[ \t]*)?[+−-][ \t]*(?=[.\d]|j[.\d])`. On today's corpus
it fires on 14 templates outside Track A, so it cannot be *gated* on until Phase 6
fixes them — but it should be reported as a **census** in the meantime, exactly the
way `degenerate_product_derivation` is. A gate that says `0 templates` for a class
that is present on 14 is worse than one that says `14, census`.

**3. Scan the question as well as the solution.** (Addresses F2.) One line in
`scan_template`; the class-3 detector is the one that needs it, and it is the one
whose fix changed question text. While doing it, decide explicitly whether class 1
and class 2 apply to questions (they should not) and record that, so the asymmetry
is a decision rather than an oversight.

**4. Fix `_DEGENERATE_PRODUCT` for the float-zero form and re-run the census.**
(Addresses F3.) `(?<![\w.])0(?:\.0+)?\s*\*\s*(?=pi\b|[A-Za-z]_?\w*\b)` recovers
`0.0*pi` without loosening the subscript guard that the current lookbehind was
carefully built for. Expect the census to go from 2 templates to 3
(`undamped_response_initial_conditions`, 142/400 instances, derivation only).

**5. Make a lost step-heading detectable.** (Addresses F5.) The cheapest fix is a
line-oriented probe rather than a marker-oriented one: any line matching
`^[ \t]*\*{0,2}[ \t]*Step[ \t]*\d+[ \t]*[:.]` whose canonical form does not start
it is malformed. That is what my class-1 probe uses, it catches the unbolded
heading, and it agrees with `t4_contract` on every one of the 150 templates today.

**6. Named templates where a defect in this diff's classes probably still lives —
the Phase 6 worklist, measured rather than guessed.** Doubled sign still emitted
after `e3ee9ab`, 150 × 400 seeds, my probe:

| in the **answer span** (2) | `lorentz_force`, `continuous_to_discrete_conversion` |
|---|---|
| in the **question** (5) | `lorentz_force`, `wave_equation_interpretation`, `continuous_to_discrete_conversion`, `nyquist_rate_determination`, `vorticity_check` |
| solution only (the rest of the 14) | `sensible_heat_constant_cp`, `pitzer_correlation_z`, `work_isothermal_virial`, `mean_variance`, `coulombs_law`, `impulse_response_from_lccde`, `signal_operations`, `system_property_linearity`, `multi_segment_rod` |

`_signed_term` / `_rect_str` are currently private to `waves_and_phasors.py`.
Promote them to a shared emission helper (`data/templates/branches/_emit.py` or
`tests`-adjacent) before Phase 6 rewrites nine templates by hand — nine
hand-written sign fixes is nine chances to write `+ {value}` again.

**7. Decide what an answer span may contain, and gate it.** (Addresses F7.)
`levenspiel_plot_interpretation` is the corpus's only trailing-commentary case and
it is worth keeping pedagogically — but then the contract should say so and the
comparator should be told. Two options, and the phase should pick one: (a) move the
`**Note:**` block *above* `**Answer:**`, which costs nothing and makes the span the
answer; or (b) declare a `**Note:**` terminator in `normalize.answer_span` and add
it to the D5.3 predicate. Option (a) is a one-line template edit and removes a
false-accept surface Track B would otherwise have to model. Related and separable:
`normalize._QUALIFIER_RE`'s 80-character lookback is shorter than this Note, so the
F4 guard Reviewer E asked for would not fire inside it.

**8. Where the evidence is thinnest, ranked.**
(i) *Nothing at all is measured about question text* — no check in T1–T7 reads
`inst.question` except for emptiness, and this phase changed 2,560 of them.
(ii) *The marker-agreement predicate only tests the marker, not the span.* D5.3's
`check_marker_agreement` asserts `marker == '**Answer:**'` and `span` non-empty; it
would pass unchanged on a span containing a whole extra section (it does, on
`levenspiel`). A span-shape assertion — no `**` heading inside the span, span
shorter than N characters, or span ends where the solution ends — would have caught
F7 as part of the gate rather than as a reviewer's aside.
(iii) *T4's `multiple_answers` counts only `present[0]`*, so a second, differently
spelled marker is never counted as a duplicate (this is F6's mechanism).

**9. What later phases should do differently.** The pattern this review keeps
finding is *the detector is narrower than the class it is named after* — F1, F3, F5
and F6 are all instances. The phase already has the right instinct (`--selftest`
plants a defect per class); what it lacks is that the plant is written by the same
person who wrote the regex, so every plant is a shape the regex already matches.
Suggested rule for Phase 6's brief: **the planted defects must be written from the
class definition in the spec, by someone who has not read the detector** — or, more
cheaply, each detector must carry at least two plants of *materially different
surface form* (`0*pi` and `0.0*pi`; `+ j-5` and `+ -5`; a malformed bold heading and
a lost bold heading). That is a fifteen-minute change to `PLANTS` and it would have
caught four of this review's five substantive findings.

**10. Spec / gate ambiguity worth a `SPEC-CHANGE`.** The Track A exit gate reads
"D5.2 corpus-wide scan clean across all four classes". Taken literally the corpus is
*not* clean across class 3 — 14 templates emit it — and the gate passes only because
the detector is narrower than the class. The gate should say what it means:
**"clean across all four classes *within Track A's eleven*, with the corpus-wide
residual enumerated as a census and assigned to Phase 6."** That is what the work
actually delivers, it is defensible, and it stops the scan's `complex_sign 0
templates` line from reading as a corpus claim that D-061 contradicts three
paragraphs earlier.

*Out-of-scope note, recorded not filed:* `stoichiometry.py` still raises
`SyntaxWarning: invalid escape sequence '\%'` on every import, on the same line the
Track A marker edit touched (`{X_A*100} \%$`). Pre-existing and outside the four
classes, as the brief says; a one-character fix (`\\%`) while the file is open.
