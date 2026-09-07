# Phase 5 Track A — item-pool impact

**P6 event: YES, on 3 of 11 templates.** Eight templates changed formatting only.
Three changed emitted text — the position of a minus sign, and one rendering of
zero. No question's *content* and no answer's *value* changed anywhere.

**Evidence.** `tests/template_integrity/phase5_instance_dump.py`, 11 templates ×
2,000 seeds, dumped from two git worktrees in **two separate processes**:

```
git worktree add --detach C:/wtm master
python -m tests.template_integrity.phase5_instance_dump C:/wtm before.json 2000
python -m tests.template_integrity.phase5_instance_dump .     after.json  2000
python -m tests.template_integrity.phase5_instance_dump --diff before.json after.json
```

The script refuses to run if a template resolves outside the tree root it was
given, and `--diff` refuses two dumps with the same root — the in-process-reload
trap that reports every instance identical and has now caught four people.

## Why this document, and not a T6 delta

T6 is **not** evidence about the item pool for this phase, and using it would be
worse than useless:

- its baseline is stale corpus-wide (D-043) and regenerating it to make a gate
  green is forbidden;
- its answer extractor **cannot see a non-scalar answer at all** (R4-5), which is
  most of Track A's scope;
- and its `breaches` list is emitted in nondeterministic order, so a naive JSON
  diff of two T6 reports reports spurious movement on 7 of 150 templates. Under
  an order-insensitive comparison T6 is fully stable across 5 runs at a fixed
  commit — 0 unstable templates — and the 7 vanish. **Measured before it was
  reported**; the first reading of it here was "T6 is nondeterministic", and that
  was wrong.

## The distinction that carries the P6 judgement

`answer span changed` and `answer body changed` are different numbers and only
the second is a P6 event. The five chemical templates changed their answer
*marker* and nothing else — their span changed on 100% of instances and their
body on 0%. A comparison that could not tell those apart would report five items
as re-answered.

## Per template

| template | class | solutions | questions | **answer body** |
|---|---|---:|---:|---:|
| `cd_dc_system_analysis` | step marker | 2000/2000 | 0 | **0** |
| `euclidean_distance_binary` | step marker | 996/2000 | 0 | **0** |
| `finite_convolution` | step marker | 2000/2000 | 0 | **0** |
| `batch_moles_vs_conversion` | answer marker | 2000/2000 | 0 | **0** |
| `flow_system_molar_flow_rates` | answer marker | 2000/2000 | 0 | **0** |
| `gas_phase_concentration` | answer marker | 2000/2000 | 0 | **0** |
| `limiting_reactant` | answer marker | 2000/2000 | 0 | **0** |
| `levenspiel_plot_interpretation` | answer marker | 2000/2000 | 0 | **0** |
| `time_to_phasor` | sign position | 1273/2000 | **1022** | **998** |
| `phasor_addition` | sign position | 2000/2000 | **1538** | **1008** |
| `decimation_aliasing_analysis` | degenerate product | 39/2000 | 0 | **39** |
| **total** | | 18,308 / 22,000 | **2,560** | **2,045** |

Generation errors: **0 before, 0 after.** `euclidean_distance_binary` changes on
996 of 2,000 because only one of its two branches (orthogonal BFSK) carried the
malformed marker; `decimation_aliasing_analysis` on 39 because its unreachable
guard was reached on 1.95% of instances.

## The three changes that are P6 events, and the argument for each

**1. `time_to_phasor`, `phasor_addition` — the minus sign moves out of the
coefficient and into the operator.**

```
- v(t) = 127.44 * sin(876*t + -86.8 deg)      - Phasor V = -127.24 + j-7.11
+ v(t) = 127.44 * sin(876*t - 86.8 deg)       + Phasor V = -127.24 - j7.11
```

Same amplitude, same phase, same phasor, same steps, same physics. What changed
is that the old text was **not a number in any notation**. An item whose question
reads `sin(876*t + -86.8 deg)` is, to that extent, testing a reader's tolerance
for malformed notation rather than phasor conversion. The difficulty does not
rise and the tested skill does not change; a defect in the item's presentation is
removed. Recorded as a P6 event anyway, because the emitted question text differs
on 2,560 instances and that is a scoping decision, not a side effect.

**2. `phasor_addition` — `sqrt(-41.2^2 + 86.45^2)` gains parentheses.** This one
is a *correctness* fix, not a presentation fix, and it is the largest single
change in this phase (D-061):

| | before | after |
|---|---|---|
| printed | `A_total = sqrt(-41.2^2 + 86.45^2) = 95.77` | `A_total = sqrt((-41.2)^2 + (86.45)^2) = 95.77` |
| evaluates to | **76.00** | 95.77 |
| T1 failures / 25 seeds | 9 | **5** |
| T1 lines checked | 142 | 152 |

The printed line did not reproduce the printed answer — P1, on every instance
with a negative component. The five surviving failures are ordinary
rounding-boundary misses (delta ≈ 0.006 against a 0.005 tolerance), a different
and much smaller defect. The template still fails T1 and is still in the corpus's
29.

**3. `decimation_aliasing_analysis` — `omega_a = 0*pi` becomes `omega_a = 0`.**
39 of 2,000 instances. The template already contained the branch that prints `0`;
it sat *after* a branch that had already caught the case, so it never ran. The
answer's value is unchanged — `0*pi` and `0` are the same number — and the item
now states it in the form its own author wrote.

## Rejected slice profile (SPEC-CHANGE 9)

**Not applicable: this phase adds no screen.** No instance is rejected, no
parameter distribution is filtered, and the seed→instance map is unchanged —
every seed still produces an instance, and the same instance, modulo the text
differences tabulated above. Stated rather than omitted, because SPEC-CHANGE 9
makes the profile a required recorded measurement in every item-pool-impact note
and "there was no screen" is the reading, not an oversight.
