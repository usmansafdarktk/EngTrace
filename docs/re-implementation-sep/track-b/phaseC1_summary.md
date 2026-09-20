# Phase C1 — Summary and close-out

**Status: COMPLETE WITH FINDINGS.** · **Dates:** 2026-09-11 – 2026-09-12
**Branch:** `redesign/phaseC1-sources-convention` · **Reviewed at** `7fc5407`

C1 fixed how a constant's provenance is written, declared what every numeric table
is and what units it carries, and **measured** which tables need a citation. It is
the stage that sizes C3.

> **Headline.** Of **108** numeric tables across the five branches, **37 need a
> citation, 30 only a plausibility window, 4 a derivation, 5 a definition, 1 is a
> validity bound, and 31 are consumed by no template** — computed from a declared
> `@kind`, a perturbation measurement and, where the measurement cannot settle it,
> a declared verdict backed by a committed check. Applied literally, the
> given-values rule would have excused every table of named-entity facts that its
> consumers restate — 25 at the reviewed ref (15 of them `property`), 16 under the
> corrected measurement. It excuses correctness, not truth.

---

## 1. What was actually here

The brief's coverage table reproduces exactly under its own predicate — **107
tables, 2,414 literals** — once my tag detector matched that predicate (§7 #1).
**"43 tables" was a definition, not a count**: the three original branches hold 38.
One more table holds numbers and no literal (`PHASE_RANGE_RAD`), so 108 are
classified.

Four things the brief did not know:

- **The acquisition's record was wrong in a way `--verify` could not see.** The
  manifest named, for the water and cyclohexane isobars, the request NIST clamps —
  the response the acquisition rejected — while the files came from the fallback
  grid. A hash says a file is unchanged, not that its record is true (D-072).
- **A clone could not have verified at all.** With `core.autocrlf=true`, git would
  have turned all 107 committed LF sources into CRLF and failed every hash.
  `.gitattributes` fixes it; a clean clone now reports exactly the 13 gitignored
  binaries missing and nothing changed.
- **The on-disk NIST saturation curves are on NIST's adaptive 601-point grid**, not
  a requested one: no `REAL_FLUID_DATA` temperature is a grid row (§10).
- **"Tagged" meant less than it looked.** Of 110 tags, 47 are deprecated forms with
  no artefact and locator — civil 22, industrial 25 — and only chemical's 47 and
  electrical's 2 resolve against a file.

## 2. Deliverables

| # | Deliverable | Where | Status |
|---|---|---|---|
| C1.1 | reference set verified; grounding route per citation table | `fetch_references.py` (two defects fixed, `--selftest`), `.gitattributes`, `phaseC1_grounding.md` | 120 present, 2 known gaps (Au, Ag), `--verify` 0; a map of **candidate** routes, each coverage claim tied to one of nine probes |
| C1.2 | one vocabulary, in the spec | spec §C1.2 (SPEC-CHANGE 20), `docs/references/README.md`, D-069 | seven classes, `[ON-DISK:LOCAL-ONLY]`, and a stated relation — `precision=` (origin) or `tol=` (verification) |
| C1.2 | citation gate reworded | C2 and C3 exit gates (SPEC-CHANGE 21) | done |
| C1.3 | classify every table, with predicate and script | `census.py`, `given_evidence.py`, `phaseC1_census.md`, D-070, D-073 | 108 classified, `census --check` clean |
| C1.4 | unit declaration per table, and a failing check | `@units` on 108 tables, `test_table_metadata.py` | 325 checks pass; 12 plants and a ratchet |
| C1.5 | validity-domain field | spec §C1.5; `domain_worklist.txt` | specified; 2 filled (pilot), 106 on the ratcheted C3 worklist |
| C1.6 | generalised resolver | `test_citations_resolve.py` (C2's P1-P4 + R1-R7) | 110 tags: 49 resolved, 14 stated, 47 LEGACY, 0 failures; 10 plants + an attachment fixture |
| C1.7 | pilot `C0`, `EPSILON_0` vs CODATA 2022 | electrical `constants.py`, D-071 | `precision=exact` and `precision=4sf`, both resolve |
| C1.R | Reviewer G + R4 triage | `reviews/phaseC1_reviewer_g_provenance.md`, D-073, §8 | every finding and §5 suggestion dispositioned |

## 3. The classification

Two inputs, never one (spec §C1.3, D-070, SPEC-CHANGE 22):

- **a declared `@kind`** — the judgement, written where it can be challenged; and
- **a measurement** — nudge one field of a table, re-run every consumer on the same
  seed, and see whether the question changed in its numbers (restated), only the
  answer changed (hidden), the question was reworded, or the template raised.

Only a `range` — a window a value is drawn from — can be a plausibility window, and
only when every consumer is measured to restate it or reads it without anything
changing (a guard). **A crash or a rewording is not evidence of either**, so such a
table needs a declared `@given` naming committed evidence, and
`given_evidence.py` is that evidence for six of them.

| class | tables | literals |
|---|---:|---:|
| CITATION | 37 | see `phaseC1_census.md` |
| PLAUSIBILITY | 30 | |
| DERIVATION | 4 | |
| DEFINITION | 5 | |
| DOMAIN | 1 | |
| UNCONSUMED | 31 | an upper bound: an undeclared literal copy is invisible until C3.8 |

## 4. The pilot (C1.7)

`C0 = 299792458` equals CODATA 2022 "(exact)" and is tagged `[ON-DISK] …
precision=exact`, not `[BY-DEFINITION]`: a defined number can still be mistyped,
and an artefact catches that. `EPSILON_0 = 8.854e-12` is CODATA's
8.854 187 8188e-12 **rounded** to four significant figures, 2.12e-5 relative below
it, and is tagged `[ON-DISK] … precision=4sf` — a rounding is a transcription with a
stated precision, checked exactly. D-025, applied: no consuming question states ε₀
(hidden in all four electrostatics templates, which print it at `.4e`/`.3e`, both
lossless), so the rounding's cost falls on a solver using the full value — 2.1e-5
relative — and none of the four declares a `TOLERANCE` (D-071).

## 5. Item-pool impact

**None, measured.** C1 changed comment lines in `constants.py` and code under
`tests/constants_integrity/` only. All 150 templates × 25 seeds hash identically
before the first header was inserted and after the last `@given` — the same tree at
two times, in separate processes.

## 6. What the phase's own checks caught

Every detector C1 added ships a self-test whose plants are judged by the failure
they **add** over a clean fixture — the rule adopted after the metadata self-test
passed three plants over a fixture that was itself failing (§7 #13). The plants that
earned their keep:

- **the census's planted consumer set** caught comprehension variables treated as
  table aliases (#7);
- **the `WINDOW2` plant**, written for the G-2 fix, caught that the fix was
  incomplete one level down (#15) — which exposed eleven tables;
- **R7** (a tag must agree with `@kind`) caught my own mis-declared soil windows
  (#30) and its own over-broad row rule (#9);
- **diffing every table's class, consumers and draws** against the pre-review census
  caught two defects my rewrite introduced (#6, #11) that no plant covered.

## 7. Errors I made in this phase

**Thirty-two.** Five shapes, and every one is a shape `phase4_summary.md` or
`phase5_summary.md` already names.

### Shape A — a detector narrower than its class (6)

| # | Error | Caught by |
|---|---|---|
| 1 | The tag regex required the closing `]` on one line: civil read 20 tagged / 495 literals against the brief's 21 / 499 | reproducing the brief's count |
| 2 | P-CONSUMER could not see a constants accessor (`chart_factor`) or an import-time copy: `CONTROL_CHART_FACTORS` read "unconsumed" | disbelieving a zero |
| 3 | Integers nudged by +1, so `C0` read as having no effect | disbelieving NO-EFFECT |
| 4 | P-CONSUMER cannot see a copied literal: `SCS_IA_RATIO` read unconsumed | **Reviewer G, G-1** |
| 5 | The table universe excluded a table built only from named constants: `PHASE_RANGE_RAD` | **Reviewer G, G-4** |
| 6 | My flow-sensitive rewrite dropped a table's own identity when it is also built from tables: `RESISTOR_SERIES_BY_TOLERANCE` lost its consumer, `MEDIA_VELOCITIES` its draw | diffing classes before/after |

### Shape B — a detector broader than its class (5)

| # | Error | Caught by |
|---|---|---|
| 7 | Module-level comprehension variables became table aliases | the self-test's planted consumer set |
| 8 | **The class #7 belonged to**: flow-insensitive aliasing across module-level loops credited `QUEUE_SCENARIOS` with 4 consumers where 2 read it. My fix for #7 had patched one surface form | **Reviewer G, §5** |
| 9 | R7 forbade a row-level `[BY-DEFINITION]` — the spec's own ΔHf° = 0 example | the real run |
| 10 | Pending tags never flushed at a table boundary: `WATER_DENSITY_KG_M3`'s tag landed on a row of a later table | reading the LEGACY list |
| 11 | Scalar tables reported as drawn by order (`SHEWHART_K_SIGMA`) | diffing draws before/after |

### Shape C — a gate that passes without testing the property (4)

| # | Error | Caught by |
|---|---|---|
| 12 | The census self-test re-implemented the draw detector inline and tested the copy | rewriting it to call the real code path |
| 13 | The metadata self-test's clean fixture was itself failing, so all three M3 plants "passed" | running the fixture alone |
| 14 | An all-crash probe rolled up as NO-EFFECT and was granted PLAUSIBILITY "as a guard" | **Reviewer G, G-2** |
| 15 | **The fix for #14 was incomplete**: `field_verdict` let RESTATED outrank STRUCTURAL and ERROR, so eleven `range` tables were PLAUSIBILITY at the reviewed ref on a measurement that could not have said otherwise | the `WINDOW2` plant |

### Shape D — a number or fact from memory (6)

| # | Error | Caught by |
|---|---|---|
| 16 | D-072's draft said 106 committed text sources; the manifest says 107 | counting |
| 17 | The grounding map's first draft carried five claims from general knowledge | re-reading against the grounded-retrieval rule |
| 18 | The C3 additions said 17 `@given` declarations; the ratchet counted 15 | the ratchet's own output |
| 19 | A scratch design note stated water's normal boiling point from memory | re-reading before use; never committed |
| 31 | D-070 says the literal rule would have excused "20 other property tables" - 23, never counted; the census measures 25 fact-kind tables (15 `property`) at the reviewed ref | checking this summary's headline against the census JSON; corrected in D-073 |
| 32 | This summary's first draft said chemical resolves 49 tags; 49 is the total, and chemical's is 47 | the same check, before commit |

### Shape E — tooling and process (11)

| # | Error |
|---|---|
| 20 | Reviewer G's brief shipped a snippet unpacking `ref.load()` as `(q, s)`; caught by running every snippet before dispatch |
| 21 | The first fetch run crashed on a Windows lock on `MANIFEST.json`, possibly my own concurrent read of it. Not established; the timing fits |
| 22 | The header-insertion script assumed one newline convention per file; chemical's is mixed. Caught by counting CR/LF before running it |
| 23 | Two pilot patches failed on a trailing space before I stopped matching the table line at all |
| 24 | The resolver's prose filter: the flag added, the filter that uses it forgotten |
| 25 | The attachment self-test crashed sorting `None` against `str` — twice, the second time because I fixed one list and not the other |
| 26 | My new metadata fixture made an existing plant's anchor non-unique |
| 27 | **A replacement script fed through a shell heredoc wrote a literal `\n` into `census.py`** — the brief's explicit warning, and the error `phase4_summary.md` §8 #18 and `phase5_summary.md` §7 #11 record. It surfaced as a syntax error; one in a string would not have |
| 28 | Parallel shell commands with `cd` moved the shared working directory under each other |
| 29 | A clean-clone test failed on Windows path length before it could test anything |
| 30 | `@kind` drew `property` vs `range` as "real material vs abstract" and declared four soil windows `property`; R7 caught the contradiction with civil's own `[POLICY]` tag (D-070) |

**What caught them.** Not one by reading a diff of my own. Four by the reviewer
(#4, #5, #8, #14); ten by a count or a before/after diff that could not be right
(#1-#3, #6, #10, #11, #16, #18, #31, #32); five by a plant or R7 (#7, #13, #15, #26, #30);
the rest by running the thing. The two C1 conclusions are those of every summary
before it — *writing a lesson down does not prevent it, and a check has to run on
the shape, not the instance* — plus one of this phase's own: **a fix is a new
detector, and needs its own plant.** #6, #11 and #15 were all introduced or left
open by a repair.

## 8. R4 triage — Reviewer G

Full table in D-073. In brief:

| # | Finding | Disposition |
|---|---|---|
| G-1 | `SCS_IA_RATIO` consumed through a copied literal | `ADOPT-NOW` (P-COPY, declared and verified) + `ADOPT-PHASE-C3` (C3.8: the template reads the table) |
| G-2 | a crash read as a guard | `ADOPT-NOW` + `SPEC-CHANGE` 22 — and the incomplete first fix, §7 #15 |
| G-3 | a graded rule resting on unconsumed tables | `PLAUSIBLE` → written response → `ADOPT-PHASE-C3` (C3.9) |
| G-4 | a table with no literal outside the universe | `ADOPT-NOW` (P-TABLE-LIVE) |
| §5 ×8 | literal-copy sweep; inline-window register; sand Gs; screening-bound kind; static guard; consumer counts; the radian phase; optional part | C3.8, C3.9, C3.2, SPEC-CHANGE 22, SPEC-CHANGE 22 + C3.10, `ADOPT-NOW` (§7 #8), `ADOPT-NOW`, nothing to triage |

**Every ADOPT-PHASE-C3 item is written into the spec's C3 deliverable list.** The
fixes were not re-reviewed by G; C3's Reviewer G reads the same instruments.

## 9. Verification

At the reviewed ref `7fc5407`, in a separate worktree, against the pre-phase
baseline — **identical on every gate**:

| check | baseline | C1 |
|---|---|---|
| T1-T8 (`run --checks all`) | T1 29 · T2 0 · T3 0 · T4 0 · T5 66 · T6 142 · T7 83 · T8 0 | identical |
| `phase5_contract_scan` | clean 150/150, 0 errors | identical |
| `derive_bindings` | committed tables agree | identical |
| `cross_pair` | gate PASS (all four terms zero) | identical |
| `score` | archive 100.0 % / 97.8 %; adversarial 95.8 % / 98.6 % | identical |
| `audit_3_8` | 80 / 80 clean | identical |
| `test_chemical_thermochemistry` | 222 checks pass | identical |

After triage, on the final tree: every constants-integrity suite and self-test
passes (census, metadata, resolver, given-evidence, thermochemistry),
`census --check` is clean, `fetch_references.py --verify` reports 0 problems, and
the corpus hash is unchanged. The T and comparator suites were not re-run after
triage: nothing they import changed, and the hash measures that no template's
output did.

## 10. What C3 inherits

1. **The deliverables C1's triage added to C3** (spec §C3): C3.8 a literal-copy
   sweep, C3.9 an inline-window register, C3.10 a static guard detector, and
   amendments to C3.2 and C3.6.
2. **Re-acquire, don't interpolate, the saturation rows.** Every `REAL_FLUID_DATA`
   temperature is off NIST's adaptive grid; extend `fetch_references.py` with
   targeted saturation requests held to a required-row check.
3. **`MEDIA_VELOCITIES` is mostly a precision question.** Thirteen rows agree with
   an on-disk refractiveindex.info dataset at the precision they are written to;
   GaP is 2.7-4.1 % high against five datasets; PTFE, polyethylene at 589 nm,
   "amorphous semiconductor glass" and muscle tissue have no optical source; two
   rows name glass families. n is a dispersion formula, and the database carries
   its own formulas page to ground an evaluator on.
4. **`CONTROL_CHART_FACTORS` derives from on-disk definitions in 353 of 384 cells**
   at the table's own precision; the other 31 differ in the last digit, mostly
   D₁/D₂/D₄ and 1/c₄. Replacing them is a P6 event.
5. **MIL-HDBK-5J's design tables are locatable** (6061: Table 3.6.2.0(b1), p. 566)
   but its text layer is one token per line: transcribe from page images.
6. **Worklists:** 47 LEGACY tags (C3.7), 106 `@domain` tables, 15 `@given`
   declarations to re-verify whenever a consumer changes.
7. **The C2.2 suite prints a stale validity NOTE** about `CP_VALID_T_MAX` for a
   template that reads `CP_PARAMS_COMBUSTION` — a number printed and not read,
   found at baseline.

## 11. Exit gate (C1)

- [x] Every numeric table in all five branches classified, with the predicate and
      script published — 108, `census --check` clean
- [x] One vocabulary in the spec (SPEC-CHANGE 20), including a local-only qualifier
- [x] The citation gate reworded to "retrievable on-disk artefact + locator"
      (SPEC-CHANGE 21)
- [x] Unit declaration and validity-domain fields specified, and a check that fails
      when missing
- [x] The generalised citation-resolution check green on the pilot and failing on
      planted defects of both named forms — a file that does not exist, and a real
      file whose content at the locator disagrees — plus eight more
- [x] Reviewer G filed; every finding and §5 suggestion triaged (D-073)
