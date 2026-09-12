# On-disk reference sources

Every `[ON-DISK]` citation in `data/templates/branches/*/constants.py` must
resolve to a file here, so a reviewer can check a value **without network
access and without a library**. That is the rule Phase C2 established (D-030,
D-034) and the one Phases C1 and C3 run on.

## How the set is acquired — `fetch_references.py`

    python docs/references/fetch_references.py              # fetch anything missing
    python docs/references/fetch_references.py --verify     # re-hash what is on disk
    python docs/references/fetch_references.py --list       # print the plan, fetch nothing

The script records every source in **`MANIFEST.json`**: URL, retrieval time,
SHA-256, licence, and the tables it grounds. **Add a new source to the script
rather than downloading it by hand** — the script is what makes an acquisition
reproducible and checkable.

**A download is judged by its content, never by its HTTP status.** A PDF must
start `%PDF`; a zip or `.xlsx` `PK`; a text file must not be an HTML page; a NIST
fluid table must carry its column header, at least three data rows, and — for a
1-atm isobar — **a 293.15 K row**; a WebBook page's `<title>` must name the
species its CAS number was requested for; a JANAF table must name the element in
its reference state. Where a host publishes a hash, it is checked, including on
files already present. A file that fails is recorded as a failure and is not
written. Every single download also has a wall-clock deadline.

Each of those checks exists because the acquisition met the failure it catches:

- **NIST silently clamps a requested temperature grid** to the triple point
  rather than rejecting it. The water, cyclohexane and benzene isobars came back
  with HTTP 200, a valid header and a full set of rows — starting at 273.16,
  279.86 and 278.67 K, with **no 293.15 K row at all**. The first run accepted
  them. **Before reading a value from any NIST table, check the row you need is
  in the grid, not merely near it.**
- **A mirror can serve a "403 Forbidden" page**; a publisher's PDF link can serve
  an HTML landing page. Both returned content that was not the document.
- **JANAF file numbers do not follow the element**: `Fe-001` is wüstite and
  `C-001` is niobium carbide. Only a header check finds the right table.
- **A canonical government URL can break.** The NRCS TR-55 link returns 404, and
  one request to it stalled a run for about ten minutes before its socket
  timeout fired. NIST/SEMATECH's `handbook.zip` URL now redirects to NIST ITL's
  home page — an HTML page with HTTP 200 at the end of the chain.

## What is here

Measured from `MANIFEST.json`: **150 files, 515.2 MB (491.3 MiB), `--verify` 0
problems.** 14 of them (509.4 MB) are gitignored binaries — see §Git. Re-run
`--verify` rather than trusting this line. (Phase C1 recorded 120 files and 13
binaries; C3 added 25 targeted saturation files, NIST SP 811 and four WebBook pages.)

| Directory | Source | Files | Grounds |
|---|---|---:|---|
| `nist_webbook/` | NIST Chemistry WebBook — Shomate coefficients and heats of formation (Phase C2) | 1 | `CP_PARAMS`, `HEATS_OF_FORMATION` |
| `codata_2022/` | CODATA 2022 adjustment (`allascii.txt`) | 1 | electrical `C0`, `EPSILON_0` |
| `nist_sp811/` | NIST SP 811 (2008), *Guide for the Use of the International System of Units* — added at C3 | 1 | unit conversions for MIL-HDBK-5J values (ksi, lb/in³) |
| `nist_fluid_properties/` | NIST WebBook fluid properties: 36 pure fluids, each a full saturation curve and a 1-atm isobar through 20 °C; plus 25 saturation files requested in temperature increments (`Type=SatP`) starting ON a temperature a table states (C3) | 97 | fluid densities and viscosities, saturated volumes, Tc, Pc, critical density |
| `nist_webbook_species/` | WebBook phase-change and condensed-phase pages, 14 species the fluid database does not carry (tungsten hexafluoride and mercury telluride added at C3, for `MANOMETER_FLUIDS`) | 28 | critical constants, heats of vaporisation, liquid Cp |
| `nist_janaf/` | NIST-JANAF, reference state: Al, C, Cu, Fe, Hg, Pb, Si, W | 8 | solid and liquid Cp |
| `nasa_tr_r132/` | Svehla (1962), NASA Technical Report R-132 | 1 | Lennard-Jones σ and ε/k |
| `mil_hdbk_5j/` | MIL-HDBK-5J (2003-01-31), Distribution Statement A — **archive.org's published SHA-1 verified** | 1 | metal E, G, ν, density |
| `usda_wood_handbook/` | USDA FPL Wood Handbook, FPL-GTR-282 (March 2021) | 1 | wood density and heat capacity |
| `refractiveindex_info/` | refractiveindex.info database (CC0-1.0) | 1 | refractive indices for `MEDIA_VELOCITIES` |
| `civil/` | USGS WSP-2339; FHWA HDS-4 and HEC-22; NRCS TR-55; NAVFAC DM-7.01 and 7.02; AISC Shapes Database v16 | 7 | civil tables that currently cite `pilot/references/public/` |
| `industrial/` | MIL-STD-105E; NIST/SEMATECH e-Handbook — the whole zip (**local-only**) and its pages 6.3.2 and 6.3.2.1 | 4 | `MIL_STD_105E_*`; definitions behind `CONTROL_CHART_FACTORS` |

**Reading cautions, from inspecting the files:**

- **NASA TR R-132's text layer is OCR noise** — but C3 found it is *parseable noise*,
  and the distinction matters. Table I(a)'s force constants come out with molecules
  spelled `N_`, `02` (a digit zero), `SF_`, `C3Hs`, `n-C4Hlo`, ε/k 558.3 rendered
  `55& 3`, and **two column layouts** (pp.22/24/25 print the method inline, p.23
  collects it elsewhere), so a parser that knows one layout silently calls the other
  unparseable. With a per-page token map it reads 15 of 17 rows — enough to establish
  that `GAS_MOLECULAR_PARAMS` is **not** Svehla's table (σ agrees for 3 rows, ε/k for
  11). Xenon and chlorine have no unambiguous token and are claimed neither way.
  **MIL-STD-105E's text layer is not parseable** — p.18 extracts as
  `":: .... m n mO ::,i:,o Vtm Loi o, batcb ah:e ."` — so it is cited `page=` only,
  with no `text=` anchor (D-075).
- **A correct title is not correct content.** Radon's WebBook page was acquired
  under the right title and carries no critical-temperature data.
- **The NIST/SEMATECH e-Handbook renders control-chart factors as MathJax, not
  tables**, so it is a secondary check. Every column of `CONTROL_CHART_FACTORS`
  can instead be derived from its definition, which needs no book.

## Where the civil and industrial copies came from

These seven and two documents were first acquired in July–August 2026 into
`pilot/references/public/`, which is **gitignored** — so the citations pointing
there resolve on one machine only. The copies here are **byte-identical** to
those originals (SHA-256 compared). WSP-2339, HDS-4 and HEC-22 were re-fetched
from their URLs and matched, which also confirms the local originals are
unaltered. TR-55, both NAVFAC manuals and MIL-STD-105E were copied from the
local originals, because no working canonical URL was known; on 2026-09-11 a
mirror for each was found whose download is byte-identical to the copy (see
§Git). The run that acquired AISC and the SEMATECH zip was stopped before it
wrote a manifest, so it cannot say whether those two were downloaded or copied,
and `MANIFEST.json` labels their retrieval time as file mtime rather than
presenting it as a retrieval timestamp. On 2026-09-11 AISC's URL re-downloaded
byte-identical; SEMATECH's no longer serves the zip at all.

## Not here, and why

**Copyrighted textbooks are never copied here.** Eighteen exist on one machine
under `pilot/references/public/full_books_*` (Das, Hibbeler, Montgomery,
Nahmias, …) and are listed by path in `MANIFEST.json` under
`local_only_copyrighted`. A citation to one is honest and **cannot be resolved
from a fresh clone**; the tag should say so.

**Known gaps** — no source for these is on disk:

- gold and silver heat capacity (NIST-JANAF does not carry them);
- R-410A, which is a blend and not in the NIST fluid database;
- radon critical constants (page acquired, data absent);
- a *table* of acentric factors — none was found; ω is derivable from the
  saturation curves by its definition;
- saturated volumes and viscosities for ethanol, acetone, isopropanol, CCl₄,
  mercury, glycerol, ethylene glycol and diethyl ether — only their WebBook
  species pages are here, no saturation curves;
- polymer, ceramic, glass and concrete properties; power-law fluid parameters;
  mixture liquids (oils, honey, blood, milk, seawater) — **not searched for** in
  this acquisition, so "no source" is not established for these.

## Git — the large binaries are not committed

**Repo owner's decision.** `.gitignore` excludes `docs/references/**/*.pdf`,
`*.zip` and `*.xlsx` — 14 files, 509.4 MB since C3 added NIST SP 811 (13 files, 507.4 MB at C1). What is committed is this README,
`fetch_references.py`, `MANIFEST.json`, and the small text, TSV, HTML and JSON
sources. Do not force-add a binary or move them to Git LFS.

**On a fresh checkout, run `fetch_references.py` before `--verify`**; until then
`--verify` reports the binaries as missing. The script re-downloads each one and
holds it to what `MANIFEST.json` recorded: the SHA-256, or — for the
refractiveindex.info archive, whose bytes GitHub does not keep stable — the
commit id in the zip comment plus a fingerprint over every member's path, CRC-32
and size. **A mismatch is refused, not recorded over**: every citation into a
file was written against the recorded copy, and a host that has since changed
its copy has handed you a document nobody checked those citations against.

**Mirrors.** TR-55, NAVFAC DM-7.01 and 7.02 and MIL-STD-105E have no working URL
at their issuing agencies that this acquisition found, so each is fetched from a
third-party mirror (Oregon DOT, vulcanhammer.net, expresscorp.com). A mirror is
accepted only because its download was measured byte-identical to the recorded
file, and the hash pin keeps it that way. An NRC ADAMS copy of TR-55 was tried
and served an HTML page.

**Not clone-recoverable: `industrial/nist_sematech_ehandbook.zip`.** Its URL now
answers 302 → 301 → NIST ITL's home page, and no other copy is known, so it
stays on the one machine that has it. Its only use was a secondary check on
`CONTROL_CHART_FACTORS`, so the two e-Handbook pages carrying their definitions
are fetched and committed instead. `industrial/nist_sematech_pmc32.htm` gives
c₄ in closed form. `industrial/nist_sematech_pmc321.htm` defines A₂, D₃ and D₄
in terms of d₂ and d₃, and says d₂ and d₃ are "tabulated in many textbooks" —
**it does not tabulate them**, so those two still need a derivation or another
source. They are the **current** pages, not the
zip's copies — the live `pmc321.htm` is 20 955 bytes, the zip's 19 656, and
they differ after whitespace is normalised. Cite the pages.

## Provenance convention (Phase C1.2)

The normative text is **`template_redesign_spec.md` §C1.2–C1.5** (SPEC-CHANGE 20).
This is the working summary.

**Seven classes, one meaning each**, and a qualifier for evidence one machine holds:

| Tag | Means |
|---|---|
| `[ON-DISK]` | a located entry in a file here that a clone can re-acquire |
| `[ON-DISK:LOCAL-ONLY]` | the same, in a copyrighted book on one machine - honest, **not resolvable from a clone** |
| `[DERIVED]` | computed here from tagged inputs or on-disk data |
| `[BY-DEFINITION]` | exact by definition; no artefact can be its warrant |
| `[POLICY: sampling-only]` | a sampling window or guard that asserts nothing about a named real entity |
| `[KNOWN-DEFECTIVE]` | checked and failed; kept as a record |
| `[UNVERIFIED]` | not checked, or no source obtainable |

**Citation grammar** - `# [CLASS] <artefact> @ <locator> [precision=… | tol=…]`:

    # [ON-DISK] codata_2022/allascii.txt @ quantity="vacuum electric permittivity" precision=4sf
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=5sf

`precision=` says the constant **is** the artefact's value at that precision (a
transcription, checked exactly); `tol=` says it **agrees with** an independent
artefact (a verification, not an origin - Reviewer G, G-7). A citation with
neither is counted `LOCATOR-ONLY`, never passed as if its value had been
compared. Locators: `quantity=` (CODATA), `T=` + `col=` (NIST TSV, on the grid),
`cas=` (WebBook JSON; C2's `NIST <CAS>` form still resolves), `page=` [+ `text=`]
(PDF; a bare `page=` is the locator-only form for a document whose text layer is
not the document), `page=` + `text=` + `mil=` (a MIL-HDBK-5J design table),
`member=` (zip), `member=` + `wavelength=` (refractiveindex.info, evaluated at the
wavelength), `sheet=` + `label_col=` + `rows=all` + `blocks=` (an `.xlsx`
**whole-table** relation: every leaf against its cell), `text=` (HTML/text).
Relations beyond `precision=`/`tol=`: `field=[key]` names one value inside a row,
`via="NAME/x"` solves for the cited quantity when the table stores a combination
of it, and `scale=<f>` converts into the artefact's unit before comparing
(SPEC-CHANGE 23).

**Table fields**, in the comment header directly above every numeric table:
`# @kind:` (`property`, `standard`, `measured-constant`, `defined`,
`mathematical`, `range`, `validity`), `# @units:` (one grammar, see
`test_table_metadata.py`), `# @domain:` (a validity domain or `none (reason)`;
C3 fills it). Two more where they apply (SPEC-CHANGE 22): `# @given: stated|guard
(evidence)` for a `range` the probe cannot settle, citing committed evidence; and
`# @copied-in: <template_id> <literal>` for a template that copies a value instead
of reading it, which the census verifies against the template's source.

**What enforces it.**

| check | what it asserts |
|---|---|
| `test_citations_resolve.py` | C2's P1-P4, and R1-R7: every tag parses; an `[ON-DISK]` file exists, is vouched for by `MANIFEST.json` and matches its SHA-256; the locator resolves; the stated relation holds; local-only is counted; payload classes state their payload; the tag agrees with `@kind`. `--selftest` plants ten defects. |
| `test_table_metadata.py` | every numeric table declares `@kind` and `@units`, every field has a unit, every unit parses; `@domain` is declared or on the ratcheted C3 worklist. `--selftest` plants eleven. |
| `census.py` | the classification: which tables need a citation and which only a plausibility window, measured by perturbing each field and watching the question; `--check` fails on a table it cannot classify. `--selftest` plants the predicate's own defect shapes. |
| `given_evidence.py` | every `@given: stated` backed by a run: each consumer's drawn value is printed in its question on every seed; a declaration citing no committed evidence fails. `--selftest` plants ten. |

| `test_plausibility.py` | C3.2: a table's own columns agree with **each other** (`MATERIAL_PROPERTIES`' E_GPa against E_ksi through the SP 811 factor, to within the rounding slack their own precision implies); every `scale=` is a factor an artefact **prints** or a pure power of ten; `SPECIFIC_GRAVITY_RANGES` inside the bound Das states. `--selftest` plants seven defects and three controls. |

**Deprecated, and counted as `LEGACY` until retagged:** `[VERIFY: X]`, `[REALISM]`,
`[DERIVABLE]`, `[ON-DISK: xlsx]`, `[ON-DISK visual]`, and an `[ON-DISK]` with no
`artefact @ locator`. **C3.7 completed this: zero LEGACY tags in all five branches**
(450 tags — 152 resolved, 19 locator-only, 279 stating a class), and every numeric
table declares `@domain`. What could not be resolved is tagged `[KNOWN-DEFECTIVE]`
(19) or `[UNVERIFIED]` (200) and listed in
[`phaseC3_residual_register.md`](../re-implementation-sep/phaseC3_residual_register.md),
**never quietly accepted**.

**Origin versus verification.** The `CP_PARAMS` coefficients are in the
Smith–Van Ness functional form, and no fetchable, citable copy of Smith–Van Ness
Table C.1 could be obtained (D-030). `nist_webbook/` documents *verification*,
not *origin*: each value agrees with an independent NIST fit, which is strong
evidence the value is right and no evidence about where it was typed from
(Reviewer G, finding G-7).

## Flame-temperature reference (Phase C2)

*Adiabatic Flame Temperatures for Oxy-Methane, Oxy-Hydrogen, Air-Methane, and
Air-Hydrogen Stoichiometric Combustion using the NASA CEARUN Tool, GRI-Mech 3.0
Reaction Mechanism, and Cantera* — ETASR / arXiv:2503.11826, retrieved
2026-09-06. Air–methane, stoichiometric, reactants at 298.15 K, 1 atm:

| model | T_ad |
|---|---|
| **complete combustion, no dissociation** | **2326.35 K** |
| chemical equilibrium, GRI-Mech 3.0 (dissociation included) | 2224.25 K |

`template_adiabatic_flame_temperature` uses a single balanced reaction with no
dissociation, so **2326 K is its correct reference**, not the ~2200 K figure the
redesign spec quotes.
