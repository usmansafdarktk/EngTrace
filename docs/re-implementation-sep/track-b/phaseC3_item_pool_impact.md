# Phase C3 — item-pool impact

**P6 event: YES, on 1 of 150 templates. A second template's QUESTIONS moved without its
answers.** Over 45,000 instances per tree, 148 templates are byte-identical in question,
solution and answer. No template gained a generation error.

| | templates | instances |
|---|---:|---:|
| answer body changed (**re-score and re-inference**) | **1** | 300 / 300 of its seeds |
| question changed, answer identical (**re-inference only**) | **1** | 33 / 300 of its seeds |
| unchanged in every respect | 148 | — |

**Evidence.** `tests/constants_integrity/c3_instance_dump.py`, 150 templates × 300 seeds,
dumped from two git worktrees in **two separate processes**:

```
git -c core.longpaths=true worktree add --detach <wt_master> c7f74dc
python -m tests.constants_integrity.c3_instance_dump <wt_master> c3_before.json 300
python -m tests.constants_integrity.c3_instance_dump .          c3_after.json  300
python -m tests.constants_integrity.c3_instance_dump --diff c3_before.json c3_after.json
```

`before` is **c7f74dc**, this branch's merge-base (the Phase C1 merge); `after` is the
branch head. The script refuses to run if a template resolves outside the tree root it was
given, and `--diff` refuses two dumps with the same root — the in-process-reload trap that
reports every instance identical and has now caught four people.

Totals from the diff: `{'q': 333, 'ans': 300, 'sol': 300, 'err_before': 0, 'err_after': 0}`.

---

## 1. `two_phase_specific_volume` — 300/300 questions, 300/300 answers

The intended P6 event, from C3.3 (2/2) `d5c8aef`. The template invented its saturated
volumes with `uniform(0.001, 0.002)` and `uniform(0.05, 2.0)` and named a real fluid in the
question; it now reads `REAL_FLUID_DATA[substance]`, re-derived from the NIST saturation
rows at each row's own stated temperature (C3.1, `5313f4d`).

```
seed 0: 'The overall specific volume of the mixture is **0.27091 m³/kg**.'
     -> 'The overall specific volume of the mixture is **0.01953 m³/kg**.'
```

Both the question and the answer move on every instance, because the volumes were both
stated in the question and used in the answer. **Re-inference and re-scoring are both
required for this template.** This is the change the phase was for: the old items asked
about a named fluid using volumes that were not that fluid's.

## 2. `chart_pair_selection` — 33/300 questions, **0** answers

Not intended, and worth the space, because it is the distinction the two-process dump
exists to make.

C3.7 corrected **29 cells** of `CONTROL_CHART_FACTORS` (`b7f7613`), every column having
been re-derived from its definition (`97a0d06`). Of those 29, **7 sit in a column this
template prints**, at n ∈ {3, 18, 19, 22, 24}:

| n | column | before → after |
|---:|---|---|
| 3 | D4 | 2.574 → 2.575 |
| 18 | D4 | 1.608 → 1.609 |
| 19 | D3 | 0.403 → 0.404 |
| 19 | D4 | 1.597 → 1.596 |
| 22 | D3 | 0.434 → 0.435 |
| 22 | D4 | 1.566 → 1.565 |
| 24 | D3 | 0.451 → 0.452 |

The template prints **all eight** constants in its question — `A2, D3, D4, d2; A3, B3, B4,
c4` — regardless of which chart pair the instance is about. Its large branch draws
`n = randint(13, 20)`, so **only n = 18 and 19 can reach it**, which is the 33 seeds.

**Why no answer moved.** The large branch computes its limits from `A3`, `B3`, `B4` and
`c4` against s-bar; `D3` and `D4` are consumed only by the small branch (`n = randint(4,
8)`). So on exactly these seeds the corrected digits are *restated and not used*. The
question changes; the answer does not.

**These items need re-inference, not re-scoring.** A model saw a question containing
`D4 = 1.608` and will now see `D4 = 1.609`; the correct answer is the same number it always
was.

## 3. Why the other 22 corrected cells moved nothing

They are in `inv_d2`, `D2`, `inv_c4`, `D1` and `d3` — columns **no template prints or
consumes**. And of the 7 printed-column changes, n = 3, 22 and 24 are drawn by no template
at all: the other two consumers of these factors, `xbar_r_control_limits` and
`process_capability`, both draw `n = randint(4, 6)`.

`process_capability` reads only `d2`, and **`d2` changed at no n** — the corrections landed
in `1/d2` (2 cells), not in `d2` itself.

## 4. The two template edits that moved nothing, and why that is the result

- **C3.3 (1/2), `c4f49f1`** — `fluid_statics`' floating-object fallback now reads its
  density rows instead of copying them. **Zero instances moved**, which is the evidence
  that the copied literals equalled the table: had they differed anywhere, this diff is
  where it would have shown.
- **Every C3.7 commit** is comment-only. Each patch script asserted its module's namespace
  was `repr`-identical with unchanged dict key order before writing; this dump is the
  independent, corpus-wide confirmation of that claim — 148 templates identical, and the
  two that moved are both explained by value changes made in C3.1/C3.3/C3.7's
  `CONTROL_CHART_FACTORS` correction, not by a tagging edit.

## 5. What this is not

This is **not** a T6 delta. T6's baseline is stale corpus-wide (D-043), regenerating it to
make a gate green is forbidden, and its answer extractor cannot see a non-scalar answer.
T6 was run separately and is unchanged at 142 failing, exactly its recorded baseline; that
is a statement about T6's own gate, not about the item pool.
