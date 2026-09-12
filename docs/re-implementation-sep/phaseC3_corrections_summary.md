# Phase C3 corrections — what the escalation did

C3 shipped with **19 measured-wrong values left standing**, because correcting a constant
moves the item pool and that was the repo owner's call. The call was made: correct what is
wrong, replace what cannot be corrected, complete the constants, and put nothing in the
backlog. This records what that turned into, including the parts that went wrong on the
way.

## The state at the end

| | at the C3 merge (`4eb6488`) | now |
|---|---|---|
| tags | 452 | **489** |
| resolved | 156 (109 value comparisons) | **232** (185 value comparisons) |
| `[KNOWN-DEFECTIVE]` | 19 | **0** |
| `[UNVERIFIED]` | 198 | **177** |
| LEGACY | 0 | 0 |
| `[DERIVED]` recomputed | 0 of 15 | **8 run here, 1 elsewhere, 7 unexecuted** |
| `tol=` with a declared basis | 0 of 10 | **17 of 17** |

Regression identical to baseline throughout: T1 29, T2 0, T4 0, T5 66, T7 83, T8 0,
T6 142, contract scan 150/150, `cross_pair` PASS, `audit_3_8` 80/80.

## What was corrected

**All 19 of §1**, plus two found while closing them: `SHEAR_MODULUS_VALUES['Aluminum
2024-T4']` (28.0 → 27.6 GPa) and `MATERIAL_DENSITIES['Tungsten']` (19600 → 19300).

The largest was `GAS_MOLECULAR_PARAMS`' σ column, re-sourced to Svehla for 16 of 17 rows.
The emitted items then confirmed the physics that justified it: butane's viscosity moved
1.537e-05 → 2.094e-05, and (5.47/4.687)² = 1.362 predicts exactly that ratio from the
constants alone — computed by the template, not by the patch that made the change.

## What was built so the corrections can be checked

Two readers and two rules, all of which **fail** a wrong value rather than merely
describing it:

- `svehla.py` — reads NASA TR R-132's force constants, refusing any token that is not
  unique on its page, and refusing positional reads entirely.
- `pubchem.py` — reads one *named* entry of an aggregator record, refusing an ambiguous
  ReferenceNumber and refusing a **range** as a value.
- **`tol=` must declare what sizes it** (§6 item 5). 13 are half a unit in the artefact's
  last printed digit and the resolver now recomputes them; 4 were fitted to their own
  residual and are counted as unsizable.
- **`UNEXECUTED` built** (§6 item 6). Eight `[DERIVED]` constants asserted since C1 are now
  recomputed on every run.

## The refusals

More was declined than corrected, and deliberately:

- **Kerosene was not re-pointed at decane.** Nonane/decane/dodecane were acquired for
  exactly that, and reading them killed the plan: 718–749 kg/m³ against Kerosene 810.
  Sourcing a row must not change what the row names.
- **Xylene** (a mixture) was not replaced by p-xylene; **Aluminium** (no alloy named) was
  not replaced by pure aluminium; **no wood** was assigned a species FPL does not say the
  row means; **Xenon's σ** was not written from a positional read that turned out to
  return helium's numbers.
- **Liquid Propane's tolerance was not re-sized** to the 0.001% its artefact supports,
  because half a unit of the *row's* own digit would not admit the value either — the
  warrant is a pressure the row never states.

## Mistakes made and caught

Recorded because a summary that omits them is worth less than one that does.

- A **range** read as a value: tungsten's "18.7-19.3" became 18.7, and five other sources
  said 19.3. `pubchem.py` now refuses ranges — the refusal exists *because* of this.
- A **positional read** of an OCR'd scan returned helium's second determination for xenon.
  Guarded, then abandoned.
- **Six PubChem tags** were written citing ReferenceNumbers carrying several values. The
  resolver rejected them; the selection was fixed, not the reader.
- The **`COMMON_LIQUIDS` window scan** truncated a genus block and nearly produced a false
  finding that Maple 740 lay above every maple FPL lists.
- A **false claim** that `OBJECT_MATERIALS` lacks tungsten; it is the last entry.
- A commit attempt through a shell heredoc, broken by apostrophes — the execution prompt
  warned about exactly that.

## What is still open

§6 item 3 — **148 mechanical rows**, counted by the reason each states. Roughly 86 are bulk
materials and mixtures no chemical database indexes (granite, concrete, sea water, SAE
oils), 19 need a mechanical-properties handbook, and 17 now record what a search actually
showed rather than "no source".

**§3's substance defects are now characterised rather than suspected**, each against the
artefact:

- **`Tungsten Hexafluoride`** — a **gas** at manometer conditions, listed as a liquid.
  P = 1 atm at 289.14 K, re-derived from the page's own Antoine row and *inside* its stated
  validity window, so the verdict needs no extrapolation. The value 12900 has no support on
  disk and nothing on disk could replace it.
- **`Tellurium Mercury`** — the WebBook *does* carry mercury telluride; both pages are
  **data-free stubs**. So no further fetch can source this row, which is a stronger
  statement than "not established".
- **`R-410A`** — a blend whose **components are also absent**, so it cannot be
  reconstructed. Not substituted: R-22 would change the row's *values* (+9.78% / −33.21%)
  and would duplicate a fluid already in the table.
- **`Mercury` ΔHvap** — the page has no vaporization section at all (verified by string
  search, not by headings), and no saturation file exists, so the derivation reaching the
  other eleven rows cannot reach it.
- **`Cork` / `Cork Board` / `Bamboo`** — searched exhaustively across all 546 pages of FPL
  and every reference directory; genuinely absent.
- **`Tetrabromoethane` / `Acetylene Tetrabromide`** — one substance under two names, both
  at 2960, over-weighting that fluid in every `random.choice` over the table.

Each needs either a materials handbook or a decision that moves the item pool (removing a
row, renaming a key). Both are the repo owner's call, and each is registered with the
evidence rather than acted on.

**Also flagged, not acted on:** `templates_annotation/annotation_app/` holds a **stale fork**
of the chemical constants with pre-C3 values (R-12 `v_g 0.0268`, R-22 `v_f 0.000845`). It is
deliberately not synced — that directory backs the annotation pilot, and a frozen snapshot
may be intentional so annotations stay reproducible against what annotators saw.

## The item pool

Three P6 tranches, each measured with two worktrees in two separate processes, 150
templates × 300 seeds:

| tranche | templates moved | note |
|---|---:|---|
| 1 — the 19 defects | **13** | q 824, ans 713; two moved question-only (re-inference, not re-scoring) |
| 2 — Decision 7 | **1** | 31/300; q, ans and sol move together |
| 3 — Tungsten + tags | **0** | a value was wrong and correcting it changed no item — two different facts |

No template gained a generation error in any tranche.
