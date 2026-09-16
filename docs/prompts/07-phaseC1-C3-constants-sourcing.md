# Prompt 07 — Phases C1 + C3: constants sourcing, and the remaining tables

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Template-integrity Phases 0–5 and constants Phase **C2** are complete and merged. `master` is at
`788833d` (Phase 5's Track B merge) plus the reference acquisition described in §"Sources already
on disk" (uncommitted when this brief was written — check `git status`). **Re-derive every number in this brief at whatever `master` is when
you start** — every prior brief in this series was corrected by what its own phase measured.

**Two stages, SEQUENCED, not parallel.** C1 lands and merges first, because C1 decides the
provenance vocabulary and the classification that C3 then applies to several hundred values.
Tagging C3's values against a vocabulary C1 has not fixed yet means tagging them twice.

---

I need you to implement **Phase C1** (sources and convention) and then **Phase C3** (the
remaining constants tables) of the constants re-grounding track.

## What EngTrace is, and why this phase matters now

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from **150
parameterised Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each). Each `template_*()` samples
physically-grounded parameters, computes an answer, and returns a `(question, solution)` pair;
the solution is the **gold reasoning trace**. The sampled parameters come from per-branch
`constants.py` tables — fluid densities, material moduli, critical properties, control-chart
factors.

**The paper has been rejected twice from ACL ARR** (January and May 2026). Its Appendix E states
the claim this phase exists to make true:

> *"We extracted and manually verified extensive lists of physical constants and material
> properties from the standard handbooks and data sources listed in Table 3."*
> — `docs/_ARR_May__EngTrace.txt`, Appendix E

Table 3 names Perry's, the NIST Chemistry WebBook and the CRC Handbook for chemical; ASM Handbook
Vol. 1, Marks' and Shigley's for mechanical; IEEE 100, the Standard Handbook for Electrical
Engineers and an "ART-DEIT Database" for electrical.

**Measured, that claim is not currently supported by the repository** — see §"What is actually
here". The electrical and mechanical tables carry **no provenance at all**; their header comments
say *"approximate"*, *"Typical"*, *"Values are representative and can vary"*. Chemical's
`REAL_FLUID_DATA` header cites *"NIST WebBook, Engineering Toolbox, and standard thermodynamic
tables"* — which is three sources for 75 values and no way to tell which value came from which.
**No source named "ART-DEIT" could be identified by web search**; do not cite it, and record it as
an open question for the paper, not a source.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — an **LLM
  annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not flag it.
- `evaluation/` is a **separate track** with known defects (D-003). **Do not fix them.**

## The rule this phase runs on — grounded retrieval only

**Every value you write, verify or correct must be read from an artefact on disk.** Not from
memory, not from your training data, not from a search-result snippet, not from a web page you
read but did not store. This is D-030's rule, and it is the whole phase:

- A value is `[ON-DISK]` **only** if the file its citation names exists under
  `docs/references/` (or, with the caveat below, `pilot/references/public/`) and contains it at the
  locator the citation gives.
- **If you need a source that is not on disk, add it to `docs/references/fetch_references.py`
  first**, run the script, and let it verify the download by content before you cite it. The
  script refuses a PDF that does not start `%PDF`, a NIST table without its column header, and a
  WebBook page whose title does not name the species its CAS was requested for. Extend it; do not
  bypass it.
- **If a value cannot be grounded, it is not `[ON-DISK]`.** It is a gap, recorded as a gap. D-033
  is the precedent for what to do next: *replace the item with one that has an authoritative
  citable source*, or keep the row explicitly tagged — **and that choice is the repo owner's**,
  because replacing a substance changes the item pool (P6).
- **Never fabricate a locator.** C2 could not obtain a citable Smith–Van Ness page and said so
  rather than inventing one (D-030). A page number you did not read is worse than no page number.

## How the two stages relate

| | C1 — sources and convention | C3 — remaining tables |
|---|---|---|
| edits | the spec, `docs/references/`, `tests/constants_integrity/`, table headers | table **values** and tags in all five `constants.py`, plus two named template defects |
| decides | the vocabulary, the classification, the unit and validity-domain fields | the values |
| gate | every table classified; the citation-resolution check generalised and green on the pilot | every `[ON-DISK]` citation resolves; plausibility suites green; zero silently-unverified tables |
| reviewer | **G** — provenance | **G** — provenance, and **H** — domain, per branch |
| risk | a vocabulary that cannot express what C3 meets | a transcription error, or a correct value used outside its validity range |

**C1 sizes C3.** The spec's C1.3 — classify every table as needing a citation or only a
plausibility window — is *"the deliverable that sizes the rest of the track"*. It was never done.
Run C3's value work before it and you will source values that the given-values rule would have
excused.

**Boxes.** C1: 12–20 h. C3: 60–105 h, and C2 on 65 values found **four** defect classes where its
spec entry listed one, so treat C3's figure as a trigger rather than a forecast. **When C3 passes
105 h, ship the branches that are re-grounded, name the rest in the residual register, and
close.** A named ungrounded table with a written reason is a result; a table tagged stronger than
its evidence is a false claim in a benchmark.

## Read these first, in this order

1. **[`template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** — the
   **Track B** section (C1, C2, C3, sync points **S1/S2**), and the review protocol **R0–R6**. Read
   R6 before you dispatch anything.
2. **[`phaseC2_summary.md`](../re-implementation-sep/track-b/phaseC2_summary.md)** — the model for this
   phase. §1 (the defect was four defects), §7 (source deviation, recorded not hidden), §8 (what
   the suite caught, *including the implementer's own transcription error*), §9 (residual
   register), §12.3–12.4 (what C3 inherits).
3. **[`DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** — append-only. Essential:
   **D-025** (constants coupled to Phase 1 oracles), **D-028** (four defect classes, and the fix
   introduced a fifth), **D-030** (NIST as primary on-disk source), **D-031** (dict key ORDER is
   part of the item pool), **D-032** (unrecorded validity ranges), **D-033** (replace rather than
   ship behind a tag), **D-034** (a test may not carry its own answer key), **D-035** (four
   provenance classes), **D-036** (split a table by consumer).
4. **[`reviews/phaseC2_reviewer_g_provenance.md`](../re-implementation-sep/reviews/phaseC2_reviewer_g_provenance.md)**
   §7, and **[`reviews/phaseC2_reviewer_h_thermochemistry.md`](../re-implementation-sep/reviews/phaseC2_reviewer_h_thermochemistry.md)**
   "Does this generalise to C3?" — both reviewers wrote directly for you.
5. **[`docs/references/README.md`](../references/README.md)**, **`docs/references/MANIFEST.json`**
   and **`docs/references/fetch_references.py`** — what is on disk, where it came from, its hash.
6. **`pilot/references/public/MANIFEST.md`** — the civil and industrial source set, and the rules
   for the copyrighted books in it.
7. **[`phase5_summary.md`](../re-implementation-sep/track-a/phase5_summary.md)** §7 and §13 — the most
   recent phase's errors. Two of its three recurring shapes will reach you: *a detector narrower
   than the class it is named after*, and *a number printed and not read*.

## What is actually here, measured

**Re-derive all of it in a `git worktree` at `master` before acting**, and publish the predicate
beside every count you produce. This series' standing rule, carried from Prompt 06: *every count
ships with either its predicate written out or a committed script that regenerates it.*

### Provenance coverage by branch

Predicate: a **numeric table** is an UPPER_CASE top-level assignment in `constants.py` containing
at least one numeric literal. It is **tagged** if any provenance tag (`[ON-DISK]`, `[DERIVED]`,
`[BY-DEFINITION]`, `[KNOWN-DEFECTIVE]`, `[VERIFY…]`, `[POLICY…]`, `[REALISM]`, `[DERIVABLE]`,
`[UNVERIFIED]`) appears inside the literal **or** in the contiguous block of comment and blank
lines immediately above it.

| branch | numeric tables | tagged | numeric literals | in tagged tables |
|---|---:|---:|---:|---:|
| chemical | 17 | 3 | 794 | 169 |
| electrical | 13 | **0** | 42 | 0 |
| mechanical | 8 | **0** | 250 | 0 |
| civil | 31 | 21 | 516 | 499 |
| industrial | 38 | 24 | 812 | 739 |

**Two cautions about that table.** First, the spec says C1.3 classifies *"all 43 tables in the
three original branches"*; under this predicate chemical + electrical + mechanical have **38**.
The difference is a definition, not an error, but C1.3 is meaningless until you publish yours.
Second, **a tag is a claim, not evidence.** "Tagged" above means a tag is present, not that its
citation resolves — see the next subsection, where most of civil's and industrial's do not, from
any machine but one.

### Chemical — what C2 did and did not cover

C2 re-grounded `CP_PARAMS`, `HEATS_OF_FORMATION` and `COMBUSTION_REACTIONS`, and added
`CP_PARAMS_COMBUSTION`, `CP_VALID_T_MAX` and `CP_COMBUSTION_VALID_T_MAX` (D-036). Those are
**closed; do not re-open them**. Note that `COMBUSTION_REACTIONS` is verified by the C2.2 suite's
element balance but carries no tag — decide under C1.2's vocabulary whether a balanced reaction is
`[BY-DEFINITION]` and tag it, rather than leaving it the one untagged C2 table.

Still ungrounded, with their consumers:

| table | literals | keys | consumed by | draws by `random.choice(list(...))` (D-031) |
|---|---:|---:|---|---|
| `CRITICAL_PROPERTIES` | 135 | 27 | `thermodynamics/volumetric_properties_pure_fluids.py` | **yes** |
| `SUBSTANCES_FOR_HEATING` | 105 | 35 | `thermodynamics/heat_effects.py` | **yes** |
| `REAL_FLUID_DATA` | 75 | 25 | `thermodynamics/volumetric_properties_pure_fluids.py` | no |
| `COMMON_LIQUIDS` | 60 | 30 | `transport_phenomena/shell_momentum_balances.py`, `viscosity_and_momentum_transport.py` | **yes** |
| `GAS_MOLECULAR_PARAMS` | 51 | 17 | `transport_phenomena/viscosity_and_momentum_transport.py` | **yes** |
| `POWER_LAW_FLUIDS` | 44 | 22 | `transport_phenomena/viscosity_and_momentum_transport.py` | **yes** |
| `COMMON_GASES` | 42 | 21 | `transport_phenomena/viscosity_and_momentum_transport.py` | no |
| `SUBSTANCES_FOR_VAPORIZATION` | 16 | 16 | `thermodynamics/heat_effects.py` | **yes** |
| `REACTIONS` | 16 | 4 | `reaction_kinetics/stoichiometry.py`, `thermodynamics/heat_effects.py` | **yes** |

**Predicate for the last column, here and in the next table:** a consuming template draws with
`random.choice(` or `random.sample(` applied to the table, directly or through `list(`/`tuple(`.
**No draw uses `sorted(`** — checked, because a sorted draw would not depend on order and would
make a "yes" false. Two mechanisms hide under one "yes": the **dict** tables are drawn through
`list(TABLE.keys())` or `list(TABLE.items())` and depend on **key order**; `SUBSTANCES_FOR_HEATING`,
`SUBSTANCES_FOR_VAPORIZATION` and `REACTIONS` are **lists** drawn with a bare `random.choice(TABLE)`
and depend on **element order**. Inserting, removing or reordering either kind changes which item a
seed produces.

### Mechanical and electrical — no provenance at all

| branch | table | literals | draws by `random.choice(list(...))` |
|---|---|---:|---|
| mechanical | `MATERIAL_PROPERTIES` (E in GPa and ksi, ν) | 81 | **yes** |
| mechanical | `MATERIAL_DENSITIES` | 56 | **yes** |
| mechanical | `FLUID_DENSITIES` | 41 | **yes** |
| mechanical | `PIPE_FLUIDS` | 32 | **yes** |
| mechanical | `MANOMETER_FLUIDS` | 20 | **yes** |
| mechanical | `SHEAR_MODULUS_VALUES` | 18 | **yes** |
| electrical | `MEDIA_VELOCITIES` (`C0 / n` for 20 media) | 20 | **yes** |
| electrical | `C0`, `EPSILON_0` | 2 | — |

Electrical's other eleven "tables" are **sampling ranges** (`FREQUENCY_RANGE_HZ`,
`GAIN_K_RANGE`, …). Under the given-values rule they need a plausibility window, not a citation —
which is exactly the judgement C1.3 exists to record, table by table, rather than assume.

### Civil and industrial — tagged, but the citations resolve on one machine only

Both branches cite sources under **`pilot/references/public/`**, and **that entire directory is
gitignored** (`.gitignore`: `pilot/*`, `pilot/references`). It exists on the machine these briefs
were written on and nowhere else. It holds two kinds of file and they must be handled differently:

- **Public-domain documents** — USGS WSP-2339, FHWA HDS-4 and HEC-22, NRCS TR-55, NAVFAC DM-7.01
  and 7.02, the AISC Shapes Database v16, MIL-STD-105E, the NIST/SEMATECH e-Handbook. These are
  now also under `docs/references/civil/` and `docs/references/industrial/` where they could be
  acquired (see §"Sources already on disk"), so their citations **can be made resolvable from a
  clone** by re-pointing them.
- **Copyrighted textbooks** — `full_books_civil_engineering/` (Das, Hibbeler, Holtz, Craig, …)
  and `full_books_industrial_engineering/` (Montgomery ISQC 7e, Nahmias, Hillier & Lieberman, …).
  **Never copy, commit or redistribute these** (`pilot/references/public/MANIFEST.md`). A citation
  to one is honest and **cannot be resolved from a fresh clone**. Your job is to say so in the tag,
  not to hide it — and, where a copyright-free route exists, to take it.

**The clearest case of a copyright-free route:** industrial's `CONTROL_CHART_FACTORS` — 408
literals, the branch's largest table — cites Montgomery 7e Appendix Table VI, a copyrighted book.
Its own header says every column can be *derived*: c4, d2, d3 from their Gamma-function and
integral definitions, and A, A2, A3, B3–B6, D1–D4 from stated identities. **A value computed from
its definition is `[DERIVED]` (or `[BY-DEFINITION]`) and needs no book at all.** That is a stronger
warrant than a transcription, and it survives a clone.

`pilot/references/public/MANIFEST.md` points to `pilot/branches/industrial_engineering/BOOKS.md`
for per-book handling; **that file does not exist** in this checkout. Do not rely on it.

### Two data-integrity defects the spec names — both still present

- **`template_two_phase_specific_volume`**
  (`chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.py`, lines 80–141) names
  a real substance from `THERMO_SUBSTANCES` and then **invents** its saturated volumes:
  `V_l = round(random.uniform(0.001, 0.002), 5)` and `V_v = round(random.uniform(0.05, 2.0), 3)`.
  It never reads `REAL_FLUID_DATA`. The spec's example: ammonia's printed `v_g` can be 0.754 m³/kg
  against a true 0.1284.
- **`template_floating_object_submersion_depth`**
  (`mechanical_engineering/fluid_mechanics/fluid_statics.py`, line 354) falls back to
  `obj_material, rho_object = "Pine Wood", 500` — a hard-coded pair outside the table's control.

### Four things C2 raised and nobody actioned

1. **D-035's four-class vocabulary is not in the spec.** The spec's C1.3 still defines two tags.
   C2 needed four and invented them correctly but off-spec; Reviewer G raised a SPEC-CHANGE and it
   was *"raised, not applied"*. **C1.2 applies it.**
2. **The "edition + page" citation gate is unsatisfiable by construction** (Reviewer H). C2.1's
   gate demanded it while C1.1 acquired a CAS-keyed JSON. Reword it to *"a retrievable on-disk
   artefact and a locator within it"* before C3 inherits the same impossible gate.
3. **C1.4's unit declarations were skipped, and nothing noticed** — *"an unenforced deliverable is
   a suggestion"* (Reviewer G). C2 added units for its own four tables only. **The gate must test
   it**, or C3 skips it too.
4. **D-025's assertion was never implemented.** *"Every constant consumed by a Phase 1 template is
   exactly representable at the precision its question states it to."* `grep -ri representab tests/`
   finds nothing. C3 changes constants; this is when it bites. Implement it as a gate item on S2.

### Corpus baseline at `788833d`

```
T1 29 · T2 0 · T3 0 · T4 0 · T5 66 · T6 142 · T7 83 · T8 0      (python -m tests.template_integrity.run --checks all)
tests.constants_integrity.test_citations_resolve      55 checks, all pass
tests.constants_integrity.test_chemical_thermochemistry   222 checks, all pass
```

T6's baseline is stale corpus-wide (D-043) and **must not be regenerated to make a gate green**;
Phase 6 owns regeneration, and sync point **S2** exists so the item pool is regenerated exactly
once, after C3. Constant changes **will** move T6 — explain that movement with a before/after
instance dump, as Phases 3–5 did.

## Sources already on disk

Acquired before this brief was written, by **`docs/references/fetch_references.py`**, and recorded
in **`docs/references/MANIFEST.json`** — URL, retrieval time, SHA-256, licence, and the tables each
file grounds. **120 files, 513.0 MB (489.3 MiB); `--verify` reported 0 problems.** The large PDFs,
zips and the `.xlsx` are **gitignored, not committed** — the repo owner's decision — so on a fresh
checkout they are absent until you run `python docs/references/fetch_references.py`, which
re-downloads each one and refuses any whose content differs from what `MANIFEST.json` recorded.
**Run the script, then `--verify`,** rather than trusting this paragraph (see §Git).

| directory | source | files | grounds |
|---|---|---:|---|
| `nist_webbook/` | C2's Shomate coefficients and heats of formation — **unchanged** | 1 | `CP_PARAMS`, `HEATS_OF_FORMATION` |
| `codata_2022/` | CODATA 2022 adjustment | 1 | `C0`, `EPSILON_0` |
| `nist_fluid_properties/` | NIST fluid properties: **36 pure fluids**, each a full saturation curve and a 1-atm isobar through 20 °C | 72 | densities, viscosities, saturated volumes, Tc, Pc, critical density |
| `nist_webbook_species/` | WebBook phase-change and condensed-phase pages for 12 species the fluid database lacks | 24 | critical constants, heats of vaporisation, liquid Cp |
| `nist_janaf/` | NIST-JANAF reference-state tables: Al, C, Cu, Fe, Hg, Pb, Si, W | 8 | solid and liquid Cp |
| `nasa_tr_r132/` | Svehla (1962), NASA TR R-132 | 1 | `GAS_MOLECULAR_PARAMS` |
| `mil_hdbk_5j/` | MIL-HDBK-5J (2003), public release — **archive.org's published SHA-1 verified** | 1 | metal E, G, ν, density |
| `usda_wood_handbook/` | USDA FPL Wood Handbook, FPL-GTR-282 (2021) | 1 | wood density and Cp |
| `refractiveindex_info/` | refractiveindex.info database (CC0) | 1 | `MEDIA_VELOCITIES` |
| `civil/` | USGS WSP-2339, FHWA HDS-4 and HEC-22, NRCS TR-55, NAVFAC DM-7.01 and 7.02, AISC Shapes v16 | 7 | civil tables now citing `pilot/references/public/` |
| `industrial/` | MIL-STD-105E; NIST/SEMATECH e-Handbook (whole zip, **local-only**) and its two control-chart pages `pmc32.htm`, `pmc321.htm` | 4 | `MIL_STD_105E_*`; the definitions behind `CONTROL_CHART_FACTORS` |

The civil and industrial copies are **byte-identical** to the `pilot/references/public/` originals
(SHA-256 compared); see `docs/references/README.md` for which were re-downloaded and which copied,
and for the mirrors that make TR-55, both NAVFAC manuals and MIL-STD-105E recoverable from a clone.

### What the acquisition met, so that you do not meet it again

- **NIST silently clamps a requested temperature grid** to the triple point instead of rejecting
  it. The water, cyclohexane and benzene isobars came back with HTTP 200, a valid header and a full
  set of rows — starting at 273.16, 279.86 and 278.67 K, with **no 293.15 K row at all**. The first
  acquisition accepted them; a grid check caught it and they were refetched. **Before you read a
  value out of any NIST table, confirm the row you need is actually in the grid**, not merely the
  nearest one. This is the same shape as every green-while-wrong check in this project: the
  response looked right on every property that was being checked.
- **A scanned PDF's text layer is not the document.** NASA TR R-132's is OCR noise; the local
  manifest says the same of MIL-STD-105E. **Transcribe from page images.**
- **A correct title is not correct content.** Radon's WebBook page was acquired under the right
  title and carries no critical-temperature data.
- **JANAF file numbers do not follow the element** — `Fe-001` is wüstite, `C-001` niobium
  carbide. Cite the file the manifest names, which was accepted on its header, not a guessed ID.

### Gaps — not on disk

- **Gold and silver** heat capacity: NIST-JANAF does not carry them (`SUBSTANCES_FOR_HEATING`).
- **R-410A**: a blend, absent from the NIST fluid database (`REAL_FLUID_DATA`).
- **Radon** critical constants: page acquired, data absent (`COMMON_GASES`).
- **A table of acentric factors**: none was found. ω is *derivable* from the saturation curves —
  see §"Grounding routes worth checking first".
- **Saturated volumes and viscosities** for ethanol, acetone, isopropanol, CCl₄, mercury, glycerol,
  ethylene glycol and diethyl ether: only their WebBook species pages are on disk, no saturation
  curves.
- **Polymers, ceramics, glass, concrete; power-law fluids; mixture liquids** (oils, honey, blood,
  milk, seawater): **not searched for** in this acquisition. "No source exists" is not established.
- **Anything grounded only by a copyrighted book**, from any machine but the one that holds it.

## Stage C1 — sources and convention

**Deliverables.**

- **C1.1 — the reference set, verified and completed.** Run
  `python docs/references/fetch_references.py` (it re-acquires the gitignored binaries and is a
  no-op for files already present), then `--verify`, and confirm every hash. For every table
  C1.3 marks as needing a citation, record which on-disk artefact grounds it — or that none does.
  Anything missing that can be legitimately obtained goes through `fetch_references.py`, never
  around it. **Large binaries are not committed** (repo owner's decision; see §Git), so a new
  binary source needs a URL whose download reproduces its recorded hash.
- **C1.2 — one provenance vocabulary, in the spec.** Merge D-035's four classes (`[ON-DISK]`,
  `[DERIVED]`, `[BY-DEFINITION]`, `[KNOWN-DEFECTIVE]`) with the tags civil and industrial already
  use (`[VERIFY: <source>]`, `[POLICY: sampling-only]`, `[REALISM]`, `[DERIVABLE]`,
  `[UNVERIFIED]`) into **one** vocabulary with one meaning per tag, and add a class — or a
  qualifier — for **"on disk, but copyrighted and local-only"**, so that honesty about resolvability
  is part of the tag rather than a comment. Apply it as a `SPEC-CHANGE`, and reword the citation
  gate per §"Four things", item 2.
- **C1.3 — classify every numeric table in all five branches**, not only the three the spec
  names, into *needs a citation* vs *plausibility window only*, applying the **given-values rule**:
  a value restated in the question text is correctness-neutral. **Measure it**: for each table,
  whether its consuming templates print the value in the question. Publish the predicate and the
  script.
- **C1.4 — a machine-readable unit declaration per table**, and a check that fails when one is
  missing.
- **C1.5 — a validity-domain field per table, where the data has one** (D-032): temperature range,
  shear-rate range, temper, phase. Specify the field here; C3 fills it.
- **C1.6 — generalise `test_citations_resolve.py`** from `CP_PARAMS`/`HEATS_OF_FORMATION` to every
  tagged table in every branch, and to every artefact type in `docs/references/`: a CAS in a JSON
  entry, a row in a NIST TSV, a page in a PDF, a file in a zip. **Keep P4** — no reference value
  may be transcribed into a test file (D-034). Reviewer G called this *blocking for C3* (G-12); it
  is C3's entry condition.
- **C1.7 — pilot, end to end: electrical `C0` and `EPSILON_0` against CODATA 2022.** Two values,
  exact, one artefact. Small enough to finish, large enough to exercise the vocabulary, the unit
  field, the resolver and the review. **Measured against the on-disk file**, the CODATA 2022
  adjustment gives `299 792 458` m s⁻¹ *(exact)*, which `C0` matches digit for digit, and
  `8.854 187 8188 e-12` F m⁻¹, of which `EPSILON_0 = 8.854e-12` keeps four significant figures.
  So `EPSILON_0` is a **rounding** of the CODATA value, not a transcription of it — decide whether
  that is a `[DERIVED]` value with a stated precision, and apply D-025 to whatever consumes it.

**Exit gate (C1).**
- [ ] Every numeric table in all five branches classified, with the predicate and script published
- [ ] One vocabulary in the spec (`SPEC-CHANGE`), including a local-only/copyrighted qualifier
- [ ] The citation gate reworded to "retrievable on-disk artefact + locator"
- [ ] Unit declaration and validity-domain fields specified, and a check that fails when missing
- [ ] The generalised citation-resolution check green on the pilot and **demonstrated to fail on a
      planted defect** — two plants of materially different surface form (SPEC-CHANGE 17): a
      citation naming a file that does not exist, *and* one naming a real file whose content at the
      locator disagrees
- [ ] Reviewer G filed; every finding and §5 suggestion triaged (R2, R4)

## Stage C3 — the remaining tables

**Branch from C1's merge.** Order the work by what is cheapest to ground and most consumed:
electrical, then mechanical, then chemical, then civil/industrial re-pointing.

**Deliverables.**

- **C3.1 — re-derived tables with citations that resolve**, per branch, under C1's vocabulary.
- **C3.2 — a plausibility suite per table**, reading every reference value from the artefact
  (D-034), **plus the declared validity domain and an assertion that every consuming template
  samples inside it** (D-032). Reviewer H: *"C3.2's consistency checks will not catch that class
  at all."* Order-of-magnitude checks are necessary and nowhere near sufficient.
- **C3.3 — the two data-integrity fixes** above. Both change emitted items: record each as a P6
  event with a before/after dump.
- **C3.4 — item-pool impact**: which templates' answer distributions move and by how much, from a
  before/after instance dump run **as two separate processes against two worktrees**
  (`tests/template_integrity/phase5_instance_dump.py` is the pattern; an in-process reload reports
  every instance identical and has caught four people).
- **C3.5 — the residual register**: every value that could not be grounded, **explicitly tagged**,
  with the reason, and a replace-or-keep decision **referred to the repo owner** (D-033).
- **C3.6 — D-025's representability assertion**, as a gate item on S2.
- **C3.7 — civil and industrial made resolvable off one machine**: re-point every citation to a
  public-domain document now in `docs/references/`; mark every citation to a copyrighted book with
  C1.2's local-only qualifier; and derive `CONTROL_CHART_FACTORS` from definitions.

### Grounding routes worth checking first — hints, not decisions

Each of these came out of the acquisition. **Verify each against its artefact before relying on
it.**

- **Critical properties.** The NIST fluid-properties saturation curves end at the critical point,
  so `Tc`, `Pc` and the critical density are on disk for every fluid NIST carries; `Vc = M/ρc` and
  `Zc = Pc·Vc/(R·Tc)` are then `[DERIVED]`. **The acentric factor has a definition**,
  `ω = −log10(Psat/Pc at T/Tc = 0.7) − 1`, computable from the same saturation curve. No
  public-domain *table* of ω was found; a derivation from on-disk data is a better warrant than one
  would have been.
- **Lennard-Jones σ and ε/k** (`GAS_MOLECULAR_PARAMS`) — NASA TR R-132 (Svehla, 1962), the
  public-domain source the textbook tables descend from. **Its text layer is OCR noise: transcribe
  from the page images**, and apply Reviewer H's first-pass check — *grep for `E-2`/`E-5`/`E5`
  literals*, the scaled-column-heading defect that was C2's class 1.
- **Metal moduli, Poisson's ratio and densities** — MIL-HDBK-5J (2003, Distribution Statement A),
  whose archive.org SHA-1 was verified on download.
- **Wood densities and heat capacity** — USDA FPL Wood Handbook, GTR-282 (2021).
- **Refractive indices** for `MEDIA_VELOCITIES` — the refractiveindex.info database (CC0). A
  refractive index is wavelength-dependent; the table's header assumes ~589 nm, so the citation
  must name the wavelength it was read at.
- **Solid heat capacities** — NIST-JANAF element tables for the metals it carries.

### Entries to check before transcribing anything

These are **questions for the domain reviewer, not verdicts** — confirm each against an artefact.

- `COMMON_GASES["Steam (Water Vapor)"]` sits under the table's *"20°C … and 1 atm"* header.
  Check against the on-disk NIST water saturation curve whether a vapour state exists at those
  conditions, and what the row's value then describes.
- `MANOMETER_FLUIDS` contains `"Tungsten Hexafluoride"` and `"Tellurium Mercury"`. For each,
  establish from a source — not from recall, including the recall of whoever wrote this brief —
  its phase at the conditions a manometer operates at, and whether it is a recognised manometer
  fluid at all.
- `MEDIA_VELOCITIES` contains `"Glass (amorphous semiconductor)"` and `"Crown Glass (typical)"`,
  which name families rather than materials with one refractive index.
- `REAL_FLUID_DATA["Refrigerant-410A ..."]` is described by its own key as a **blend**, and it is
  **not** in the NIST fluid-properties list (measured against that list during acquisition), so
  its saturation volumes need a mixture source or a replacement.
- `POWER_LAW_FLUIDS` (ketchup, mayonnaise, toothpaste …) and the mixture rows of `COMMON_LIQUIDS`
  (engine oils, honey, blood, milk, seawater) — **no source for these was searched for or acquired**
  in the pre-phase acquisition. Whether a public-domain authoritative source exists is open: find
  out before concluding they are ungroundable, and do not read this line as that conclusion.

## Rules carried forward from C2 and Phase 5

- **A test may not carry its own answer key** (D-034). Read every reference value from the
  artefact the citation names; `test_citations_resolve.py` P4 enforces this and must keep doing so.
- **Two questions, not one:** *does the value agree with the reference?* and *does the reference
  exist?* C2 shipped a green 215-check suite that answered the first while never asking the second.
- **Expect more than one defect class per table** (D-028). C2's one-line spec entry was four
  classes, and the fix introduced a fifth.
- **A table that has just been corrected is exactly when a fresh transcription error is most
  likely** (C2 §8): the implementer wrote `12.5` for `1.25` while removing a 10× error. Diff every
  corrected mantissa against the original.
- **Do not reorder a dict's keys** (D-031). Regrouping `CP_PARAMS` by family changed which substance
  274 of 300 seeds drew. Every **yes** in the tables above is order-sensitive. If an entry must be
  removed or replaced, measure the item-pool effect of the position change as well as the value.
- **A correct value used outside its range is a defect** (D-032). Declare the domain; assert the
  consumers stay inside it.
- **Split a table by consumer rather than refit it for both** (D-036), when two templates need
  incompatible ranges.
- **Every detector you add gets two planted defects of materially different surface form, written
  from the class definition rather than from the detector** (SPEC-CHANGE 17). Phase 5's four
  narrowest detectors all passed a self-test whose plants were written by the same hand.
- **A gate needs more than one term** (SPEC-CHANGE 18). "Every row has a tag" is passed by a table
  tagged stronger than its evidence; "every citation resolves" is passed by a citation to the wrong
  row of the right file.
- **Every count ships with its predicate or a script**, and **read every number your tools print**
  — Phase 5's `cross_pair` printed `false rejects 56` on every run and nobody put it in a claim.

## Verification

```
python docs/references/fetch_references.py --verify
python -m tests.constants_integrity.test_citations_resolve
python -m tests.constants_integrity.test_chemical_thermochemistry
python -m tests.template_integrity.run --checks all          # T1-T8, ~2-3 min
python -m tests.template_integrity.phase5_contract_scan
python -m tests.comparators.derive_bindings                  # constants feed gold answers
python -m tests.comparators.cross_pair
python -m tests.comparators.score
python -m tests.trace_schema.audit_3_8
```

**`derive_bindings` and `cross_pair` are not decoration.** Phase 5 bound 120 templates to the
comparator by deriving their bindings from gold output; a constant change that alters an answer's
shape or unit can unbind one. **Re-run both after every branch** and treat a newly unbound template
as a finding.

Set `PYTHONIOENCODING=utf-8` before printing template output — the console is cp1252 and crashes on
`→`/`Σ`/`²`/`·`. **Write every patch as a file, never through a shell heredoc or a double-quoted
`bash -c`**: Phase 4 recorded heredoc escaping costing it a fix four times (`phase4_summary.md`
§8 #18), and Phase 5 lost four more to heredoc escaping, CRLF anchors and backtick command
substitution (`phase5_summary.md` §7 #11). Several files in this repo are CRLF.

## Reviews — read R0–R6 before dispatching any

- **Reviewer G — provenance (C1 and C3).** *The gate.* Resolve every citation in the tables under
  review against the artefact it names, and report whether the provenance **record** is honest and
  complete. G is not asked whether the numbers are right.
- **Reviewer H — domain (C3).** *The gate.* Independently re-derive a sample from the cited
  artefacts **without reference to the implementer's transcription**, and confirm the plausibility
  targets are correct physics rather than self-consistent. For C3's breadth, **dispatch H once per
  branch** rather than once for everything — R6: one gate property per review.

**Scope every review to R6**: one mandatory gate task, a time box stated as a number, tooling
supplied as working code, everything settled fenced off, and a measurement-ownership register with
**no number commissioned twice**.

**Independence, concretely:**
- Code comments, commit messages and summaries are **claims under test**, not background (R1.1).
- Give each reviewer a **frozen SHA**, and then **do not commit *or write to the working tree*
  until it files**. A reviewer reading at a frozen SHA sees uncommitted edits too. **Phase 5 broke
  this** — its implementer began actioning a review while the reviewer was still filing, and the
  reviewer had to annotate its own report with the ref its numbers described.
- Tell reviewers to read with `git show <sha>:<path>` and **never `git show <sha> --stat`**.
- Do not let two reviewers see each other's briefs or reports.
- **Commit each report unmodified and before any fix.**

**Every §5 suggestion is triaged (R4) before the stage closes**; an untriaged suggestion blocks the
gate exactly as a CONFIRMED finding does. A `SPEC-CHANGE` amends the spec; an `ADOPT-PHASE-<N>`
edits that phase's deliverable list.

## Git

- Branch **`redesign/phaseC1-sources-convention`** off `master`. Merge it with `--no-ff` when C1's
  gate passes; branch **`redesign/phaseC3-remaining-tables`** from **that merge**.
- **Large source binaries are NOT committed — the repo owner's decision.** `.gitignore` excludes
  `docs/references/**/*.pdf`, `*.zip` and `*.xlsx` (13 files, 507.4 MB, measured with
  `git ls-files --others --ignored`). **Do not force-add one, do not set up Git LFS for them, and
  do not loosen the rule.** What is committed: `fetch_references.py`, `MANIFEST.json`, `README.md`
  and the small text, TSV, HTML and JSON sources.
- **A fresh checkout gets the binaries back from `fetch_references.py`**, which holds every download
  to what `MANIFEST.json` recorded — its SHA-256, or for the GitHub archive its commit id and a
  content fingerprint — and refuses a mismatch instead of overwriting the record. TR-55, NAVFAC
  DM-7.01 and DM-7.02 and MIL-STD-105E come from **third-party mirrors**, accepted only because each
  download was measured byte-identical to the recorded file.
- **One binary cannot be recovered from a clone:** `industrial/nist_sematech_ehandbook.zip`. Its URL
  now redirects to NIST ITL's home page. Its only role was a secondary check on
  `CONTROL_CHART_FACTORS`, so the two e-Handbook pages carrying their definitions are fetched and
  committed instead: `industrial/nist_sematech_pmc32.htm` gives c₄ in closed form, and
  `industrial/nist_sematech_pmc321.htm` defines A₂, D₃ and D₄ from d₂ and d₃ but **does not
  tabulate d₂ or d₃**. They are the **current** pages and differ from the zip's copies. Cite the
  pages, not the zip.
- **A new binary source must be clone-recoverable**: a URL, added to `fetch_references.py`, whose
  download reproduces the recorded hash. A source that exists on one machine only must say so in
  its tag, as a copyrighted book's citation must.
- **Never commit anything from `pilot/references/`.** It is gitignored for a reason.
- Logical commits, not one lump. **Do not push.**
- End commit messages with: `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`

## Exit gate (C3)

- [ ] Every table C1.3 marked *needs a citation* is `[ON-DISK]` with a citation that **resolves**,
      `[DERIVED]`/`[BY-DEFINITION]` with its derivation, or in the residual register with a reason
- [ ] **Zero silently-unverified tables** — every numeric table carries a tag under C1's vocabulary
- [ ] Plausibility suites green, reading reference values only from artefacts (P4 intact)
- [ ] A declared validity domain per table that has one, and an assertion that consumers stay
      inside it
- [ ] Units declared per table, and checked
- [ ] The two data-integrity defects fixed, each with a before/after dump
- [ ] D-025's representability assertion implemented and green
- [ ] Civil and industrial citations resolvable from a clone wherever the source is public domain;
      every copyrighted citation carrying the local-only qualifier; `CONTROL_CHART_FACTORS` derived
- [ ] Item-pool impact stated from two-process before/after dumps; T6 movement explained;
      **baseline not regenerated**
- [ ] `derive_bindings`, `cross_pair`, `score`, `audit_3_8`, T1–T8 — no regression, **measured**
- [ ] Reviewers G and H filed; every finding and §5 suggestion triaged; every `SPEC-CHANGE` actioned
- [ ] `phaseC1_summary.md` and `phaseC3_summary.md` in the shape of `phaseC2_summary.md`, including
      **your own errors**. A summary recording none is not a clean phase; it is an unexamined one.

## Constraints and cautions

- **Do not modify a check or an oracle to make something pass.** A check that looks wrong is a
  finding and a `SPEC-CHANGE`.
- **Do not edit a template** beyond the two named C3.3 defects without recording it as a P6 event.
  A constant change that alters an item is a scoping decision, not a side effect.
- **Do not replace a substance, material or fluid on your own authority.** D-033's replacements were
  *repo owner directed*. Propose; do not decide.
- **P6 is a real constraint.** A corrected density moves every answer that consumed it. That is the
  point of the phase and it still has to be measured and stated.
- **Where you are uncertain, say so and name the evidence that would settle it.**

## Open items inherited — carry forward, fix only if they block your gate

- **D-003** — do the raw `inference_results/` generations still exist? Gates any corrected results
  table.
- **Phase 6's D6.7–D6.11** (the doubled-sign residual, an answer-span shape assertion, a mixed
  numeric/categorical `AnswerSpec`, the eight undecidable `symbolic` templates, per-item units) —
  Phase 6's, not yours, but **D6.11 (unit declarations) overlaps C1.4**; coordinate rather than
  build two unit schemes.
- **`sympy` is not in `requirements.txt`** (D-053).
- **`CP_PARAMS` origin vs verification** (Reviewer G, G-7) — unfixable without a citable
  Smith–Van Ness copy; recorded, not closed.
