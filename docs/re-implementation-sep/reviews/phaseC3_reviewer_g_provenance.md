# Phase C3 — Reviewer G (Provenance)

Frozen ref `f747fb6845b330958e0abfe7c2d7eac9849942bc`. Working tree clean at that commit.
Source read with `git show f747fb6:<path>`; suites run in the working tree.
Mandatory task time-boxed to 45 minutes; §5 separately time-boxed to 15.

Scope: tags, units, domains, register honesty. Per-branch VALUE correctness belongs to
the five Reviewer H's and was not sampled. C1's census classification (D-073), the
item-pool impact, and T1–T8 were taken as settled and not re-derived.

---

## 1. Verdict

**PASS WITH FINDINGS** — every re-derived count is exact, but the C3.2 domain suite is red
(28 failures), and three of the gate's own numbers count claims the instruments never check.

---

## 2. Independent re-derivation

Every number below was re-derived from the corpus, not read from a document. The counting
regex was written against spec §C1.2's vocabulary, independently of the resolver.

| Claim | Re-derived | Agreement |
|---|---|---|
| 450 tags across 5 branches | **450** | agrees (see note) |
| 152 resolved | **152** | agrees |
| 19 locator-only | **19** | agrees |
| 279 stating a class | **279** | agrees |
| 0 LEGACY | **0** | agrees |
| 0 unresolvable-from-clone | **0** | agrees |
| 108 numeric tables, all with `@kind`/`@units`/`@domain` | **108**, zero missing any of the three | agrees |
| `@domain` worklist empty | empty (header comments only) | agrees |
| 19 `[KNOWN-DEFECTIVE]` | **19** | agrees |
| 200 `[UNVERIFIED]` | **200** | agrees |

Commands:

```bash
# tag census, independent of the resolver
python -c "import re,glob;from collections import Counter;c=Counter();[c.update(re.findall(r'#\s*\[(ON-DISK:LOCAL-ONLY|ON-DISK|DERIVED|BY-DEFINITION|POLICY|KNOWN-DEFECTIVE|UNVERIFIED|LEGACY)',open(f,encoding='utf-8').read())) for f in sorted(glob.glob('data/templates/branches/*/constants.py'))];print(dict(c),sum(c.values()))"
grep -o "\[KNOWN-DEFECTIVE\]" data/templates/branches/*/constants.py | wc -l   # 19
grep -o "\[UNVERIFIED\]"      data/templates/branches/*/constants.py | wc -l   # 200
python -m tests.constants_integrity.test_citations_resolve
python -m tests.constants_integrity.test_table_metadata
```

Per branch, resolver against my own count: chemical 138, electrical 34, mechanical 206,
civil 34, industrial 38. My line-initial regex returned **449**, one short, in civil —
see Finding 6: one tag is written mid-sentence and is invisible to any line-initial
reader, including the `grep` recipe the residual register publishes as its own check.
Once that tag is counted the two methods agree exactly.

The whole partition reconciles with no slack, which is the real test of the summary line:

- ON-DISK 162 + ON-DISK:LOCAL-ONLY 9 = **171** = 152 resolved + 19 locator-only.
- POLICY 41 + DERIVED 15 + BY-DEFINITION 4 + KNOWN-DEFECTIVE 19 + UNVERIFIED 200 = **279** stated.
- 171 + 279 = **450**.

The 108 numeric tables were enumerated by calling the suite's own `numeric_tables()` and
reading `@kind`/`@units`/`@domain` off each header myself rather than trusting M1–M3;
no table was missing any of the three fields.

---

## 3. Findings

### F-1 (CONFIRMED, blocking) — the C3.2 domain suite is red: 28 failures

Spec §C1.5: "**C1 specifies it; C3 fills it** and asserts every consumer samples inside it
(D-032)." C3 filled `@domain` on all 108 tables, but the instrument that makes the second
half of that sentence true does not pass.

```bash
python -m tests.constants_integrity.test_consumer_domains   # exit 1, "28 FAILURES"
```

All 28 are `D1 … NO-PROBE - no probe and no implicit claim`: the suite cannot observe what
domain value the consumer actually used, so conformance is **unmeasured**, not measured-and-
clean. They span every branch with a real (non-`none`) domain — `chemi.HEATS_OF_FORMATION`,
`chemi.CP_PARAMS`, `chemi.CP_PARAMS_COMBUSTION`, `chemi.COMMON_LIQUIDS` (6 consumers),
`chemi.COMMON_GASES`, `civil.UNIT_WEIGHT_WATER_KN_M3` (6 consumers), `civil.SCS_CURVE_NUMBERS`,
`indu.CONTROL_CHART_FACTORS` (3), and all three `indu.MIL_STD_105E_*` tables.

The suite's own `--selftest` passes 7/7, so the failures are its verdict on the corpus, not
a broken instrument. No phase document mentions it: `grep -rn "consumer_domains\|NO-PROBE"
docs/re-implementation-sep/*.md` returns nothing, so this is neither triaged nor waived.

**Impact.** The exit gate reads "plausibility suites green; … zero silently-unverified
tables", and `test_plausibility` *is* green — but the suite that checks domains is not, and
is not named in the gate. The register's §7 points the reader at `domain_findings.txt` as
the place domain excursions live; that file is maintained by this red suite, and its single
listed excursion is the only one the suite could actually measure. 28 table/consumer pairs
carry a declared domain that nothing verifies. That is the gate's own failure mode.

### F-2 (CONFIRMED) — R6 never checks residual-register membership, and 20 residuals are in the register only as a remainder

Spec §C1.2 states the check for `[KNOWN-DEFECTIVE]` as "the error is stated; **the row is in
the residual register**", and for `[UNVERIFIED]` as "the reason is stated; **the row is in
the residual register**". The implementation (`test_citations_resolve.py` ~line 625) tests
only the first clause:

```python
if cls in ('DERIVED', 'BY-DEFINITION', 'KNOWN-DEFECTIVE', 'UNVERIFIED'):
    payload = tag['body'] or tag['bracket'].partition(':')[2].strip()
    if not payload:
        failures.append(f'R6 {where}: [{cls}] states no …')
    else:
        counts['stated'] += 1
```

Nothing reads `phaseC3_residual_register.md`. A residual tag with any non-empty reason
passes R6 and is counted "stated" whether or not it ever reaches the register.

This is not hypothetical. The register's §2 class table accounts for the 200 `[UNVERIFIED]`
by grep phrase and closes with "not covered by any phrase above | **24** | mechanical 8,
chemical 7, electrical 4, civil 4, industrial 1". Only electrical's 4 are itemised. The
other **20 residuals exist in the register as an integer and nothing else** — no constant
name, no reason, no decision. Under §C1.2's wording those rows are not in the register.

**Impact.** The register is the phase's answer to "zero silently-unverified tables", and the
one rule that would keep it complete is unimplemented. Adding a tagged residual that never
reaches the register is currently a green operation.

### F-3 (CONFIRMED) — "152 resolved" mixes value comparisons with CAS-identity checks

47 of the 152 are C2-form tags (`# [ON-DISK] NIST <CAS> …`). `classify_tag` routes them to
`ON-DISK-C2`, and `check_source` credits them directly:

```python
if cls == 'ON-DISK-C2':
    counts['resolved'] += 1         # resolved by c2_checks, P2
```

`c2_checks` verifies that the row names a CAS, that the CAS appears in the reference JSON,
and that it matches the entry for that species. It **never compares the constant to a
value.** So 47 tags (31% of "resolved") are counted in the same bucket as citations that
passed a real `precision=`/`tol=` comparison against an artefact.

Reproduce: `python -m tests.constants_integrity.test_citations_resolve` (chemical resolved
105), then count CAS-form citations — 47 — leaving 58 genuinely value-compared in chemical.

**Impact.** This is the attack class the brief names: a locator that resolves, presented in
the headline number as though the value had been checked. The values may well be checked by
`test_chemical_thermochemistry.py`, and many tags state their own comparison in prose
("Cp298 29.81 vs NIST 29.86") — but the gate's number does not distinguish them, and a
reader of "152 resolved" will overcount machine-verified values by 47.

### F-4 (CONFIRMED) — `[DERIVED]` is never recomputed, and the promised `UNEXECUTED` count does not exist

Spec §C1.2's check column for `[DERIVED]`: "recomputed wherever a derivation is registered
with the resolver; **otherwise counted `UNEXECUTED`**." Neither half is built. The resolver's
counter dict is `tags/resolved/locator_only/unresolvable_from_clone/legacy/stated` — there is
no `UNEXECUTED` key, no derivation registry, and no recomputation path. All 15 `[DERIVED]`
tags fall through the same branch as `[UNVERIFIED]` into `counts['stated']`.

```bash
grep -rn "UNEXECUTED" tests/constants_integrity/   # no match
```

Two derivations *are* independently executed elsewhere — `indu.Z_QUANTILES` against
`statistics.NormalDist().inv_cdf` and `indu.CONTROL_CHART_FACTORS` via
`control_chart_factors.py` — and chemical's Shomate refits state their worst residuals. But
the resolver neither knows nor reports this, so the phase cannot say how many of its 15
derivations reproduce. SPEC-CHANGE 23 revised three rows of this same table and left this
row stale; since that amendment is itself under review, now is the moment to fix it.

### F-5 (PLAUSIBLE) — `@domain: none` is used for "no condition is *stated*", which exempts the table from the only domain check

`parse_domain()` returns `None` for `none (reason)`, and the runner skips those tables
entirely — `none` is not a weak domain, it is *no check at all*. Several tables declare
`none` with a reason that **admits the values depend on conditions**:

- `chem.POWER_LAW_FLUIDS` — `none (K and n are shear-rate and temperature dependent; none is stated)`. Its own `[UNVERIFIED]` tag says "K for a shear-thinning fluid is not a constant of the substance — it depends on the shear range and temperature the fit was made over".
- `chem.GAS_MOLECULAR_PARAMS` — `none (Lennard-Jones parameters; no fit temperature range is stated)`.
- `mech.MATERIAL_DENSITIES` — `none (the table states no temperature)`; `mech.PIPE_FLUIDS`, `mech.MANOMETER_FLUIDS` — `none (the table states no conditions)`; `mech.MATERIAL_PROPERTIES`, `mech.SHEAR_MODULUS_VALUES` — `none (the table states no temperature, form or temper…)`.

The same fact is graded two ways. In the *tags*, "the table states no conditions" is a
defect: 12 mechanical `[UNVERIFIED]` rows read "the table states no temperature or pressure,
so no row of an artefact is THE row". In the *metadata*, it is `@domain: none` — nothing to
check. Declaring `none` on a condition-dependent table converts an open question into a
silent exemption, and does it on exactly the tables F-1 shows are hardest to probe.

Not graded CONFIRMED because `none (<reason>)` is literally permitted by §C1.5; the defect
is in the spec's vocabulary as much as in the declarations. See §5.

### F-6 (CONFIRMED, minor) — one tag is off the §C1.2 citation grammar and invisible to the register's own recipe

`civil_engineering/constants.py:351`, on `CV_RANGES_M2_YR`:

```
# so the ranges deliberately overlap.  [POLICY: sampling-only. Das PGE has no
# general cv table and the DM-7.01 cv-LL chart is image-only; any sampled cv
# MUST appear verbatim as a given value in the question text, …]
```

The grammar is `# [CLASS] <artefact> @ <locator>` — tag first. This one sits mid-sentence
after prose and carries three sentences inside the bracket. The resolver finds it (it scans
for `[CLASS` anywhere in a comment), so the corpus count is right; but every line-initial
reader misses it, which is why my independent count read 449. The residual register
publishes `grep -c` recipes as the way to re-derive its numbers, and this is a worked example
of those recipes disagreeing with the resolver.

### F-7 (PLAUSIBLE) — public artefacts sit uncited while copyright-locked ones carry the citations

Comparing every cited artefact against `MANIFEST.json`'s vouched-present list:

```bash
python <scratch>/cited.py   # artefacts present and manifest-vouched, cited by zero tags
```

- **`civil/navfac_dm7_01_soil_mechanics.pdf` (49.5 MB) and `navfac_dm7_02_foundations.pdf` (42.3 MB): cited by zero tags**, while civil's soil tables rest on 7 `[ON-DISK:LOCAL-ONLY]` citations to Das — which §C1.2 defines as *not resolvable from a clone*. `CV_RANGES_M2_YR`'s tag shows DM-7.01 was consulted and found image-only **for that one table**; `SPECIFIC_GRAVITY_RANGES`, `PERMEABILITY_RANGES_CM_S`, `FRICTION_ANGLE_RANGES_DEG` and the two `TERZAGHI_*` tables never cite it. This is the MIL-STD-105E pattern: clone-checkable public evidence on disk, unused, next to a citation that a clone cannot check.
- **`usda_wood_handbook/FPL-GTR-282_2021.pdf` (104 MB): cited by zero tags**, with 8 mechanical wood rows `[UNVERIFIED]`. Here the tag is *honest* — "The USDA Wood Handbook is on disk but was not read to row level in C3" — so this is unexploited, not misrepresented. It is nonetheless 104 MB of acquired evidence against 8 open residuals.
- **`nist_sp811/…pdf`: cited by zero tags**, yet SP 811 conversion factors are load-bearing in prose throughout (civil's `[DERIVED]` at line 28 quotes its p.61) and in the register's own arithmetic. The repo's most-invoked conversion authority is machine-checked nowhere.

Not CONFIRMED because establishing that DM-7.01 actually tabulates these specific soil
ranges needs a page-level read I could not complete inside the box.

### F-8 (PLAUSIBLE) — the register's `[UNVERIFIED]` class table is self-sealing

The eight classes sum to exactly 200 because the last is a remainder ("not covered by any
phrase above | 24"), so the per-class counts cannot be independently falsified: an error in
one class is absorbed by the remainder. Four classes reproduce exactly (19, 19, 15, 8). One
diverges: the register gives "…is not on disk" (a named standard) as **8**, and

```bash
grep -o "is not on disk" data/templates/branches/*/constants.py | wc -l   # 10
```

The register's 95-entry headline class is defined by two alternative phrases; the first
alone returns 58. Low impact — the 200 total is exact and independently confirmed — but the
class breakdown, which is the part the register uses to argue "it is one decision, not 95",
is not reproducible from the recipe given.

---

## 4. Falsification attempts that failed

These are the attacks the brief asked for that I ran and **could not** make stick. They are
what makes the findings above meaningful.

- **An `[ON-DISK]` PDF presented as though its value had been compared.** Could not
  construct. The resolver forces `value = None` on the plain-PDF path, so a PDF tag that
  declares `precision=` or `tol=` *fails* R4 ("states a relation, but a pdf locator yields
  no machine-readable value to compare") rather than passing. PDFs can only be LOCATOR-ONLY
  or the `mil=` design-table form, which does parse a value. The 19 locator-only tags are
  honestly counted in their own bucket and never folded into "resolved". This is the one
  place the vocabulary is strictest, and it holds.
- **A `[POLICY: sampling-only]` on a non-`range` table.** None. I read the declared `@kind`
  of all 41 POLICY-tagged tables; every one is `@kind: range`. R7 is real and it bites
  (the self-test plants it and catches it).
- **A `[KNOWN-DEFECTIVE]` missing from the register, or with a wrong stated error.** All 19
  are itemised in §1 of the register with a named artefact and a signed error, and the
  branch split (mech 14, civil 3, elec 2) reproduces exactly against the corpus. The
  arithmetic I spot-checked is right: civil's `rho*g_n = 998.21 x 9.80665 = 9789.1 Pa/m`
  → 9.7891 kN/m³ against 9.81 is +0.21%, as stated.
- **A LEGACY tag hiding behind a deprecated spelling.** My regex included `LEGACY`,
  `VERIFY`, `REALISM`, `DERIVABLE` and the `[ON-DISK: xlsx]`-style variants; zero matches
  corpus-wide. The C3.7 retagging is complete.
- **A numeric table missing `@kind`, `@units` or `@domain`, or hiding on the worklist.**
  I enumerated all 108 tables myself and read the three fields off each header: none
  missing, and `domain_worklist.txt` contains only comments. The ratchet is honest.
- **MIL-STD-105E, the instance the phase says it fixed.** Verified genuinely fixed:
  `industrial/mil_std_105e_sampling.pdf` is present, manifest-vouched, and carries 3 live
  citations. The claim is true.
- **The instruments checking themselves.** `test_citations_resolve --selftest` passes 24
  planted defects (including the subtle ones: a decimal tie rounded on the binary float, a
  design table read in the wrong column, a unit scale in the wrong direction, an
  interpolation-dependent rounding); `test_consumer_domains --selftest` passes 7. P4
  additionally forbids reference values hardcoded in a test file. I could not find a class
  the self-tests claim to cover that they do not.

---

## 5. Further probing and improvements

*Separately time-boxed to 15 minutes; forward-looking, not part of the verdict.*

**Where the evidence is thinnest.** Not in the citations — it is in the *sampling-only*
claim. R7's whole test is "the census does not measure it SOME-HIDDEN". But 4 POLICY tables
measure `ERROR` (the probe crashed, so the table is read and the probe cannot say how) and
9 measure `INDETERMINATE`. For those, `classify()` accepts a **self-declared `@given:
stated|guard`** in place of a measurement. So the corpus's weakest provenance claims rest on
an author's declaration that the instrument was explicitly unable to confirm — and they pass
the gate silently, because only `SOME-HIDDEN` fails. This is the same shape as F-1: the
tables hardest to measure are the ones the check excuses. I would spend the next hour
forcing those 13 to a measured verdict (fix the crashes, or `--seeds` them harder) rather
than on any further citation work.

**What I would probe with more time.** (1) Land F-7: read DM-7.01/DM-7.02 to page level for
the five civil soil tables, and demote the Das `[ON-DISK:LOCAL-ONLY]` citations to
cross-checks if a public source carries the same ranges — that converts 7 clone-unverifiable
citations into clone-verifiable ones, which is worth more than any new artefact acquisition.
(2) Read the USDA Wood Handbook to row level for the 8 wood residuals. (3) Verify the 47
C2-form tags' values actually are compared in `test_chemical_thermochemistry.py`, which
would downgrade F-3 from a verification gap to a reporting gap.

**What generalises.** Three of my findings are one bug: *a count that credits a claim the
instrument never checked* — F-2 (register membership), F-3 (C2 tags credited "resolved"),
F-4 (derivations credited "stated"). The fix generalises past C3: **every summary bucket
should be emitted by the code path that did the work**, so an unchecked claim lands in its
own named bucket rather than in a neighbour's. Phase 6's D6.11 answer-unit declarations will
import this grammar, and will inherit this failure mode unless the counters are split first.

**What in the spec is wrong or ambiguous** (the spec is under review, so these are proposals):

1. **§C1.2, `[UNVERIFIED]`/`[KNOWN-DEFECTIVE]` check column** — "the row is in the residual
   register" is unimplemented (F-2). Either build it (R6 parses the register and matches on
   constant name) or strike the clause. Building it is cheap and closes the phase's main
   honesty hole.
2. **§C1.2, `[DERIVED]` check column** — "recomputed … otherwise counted `UNEXECUTED`"
   describes an instrument that does not exist (F-4). SPEC-CHANGE 23 already rewrote three
   rows of this table for exactly this reason and missed this one; fix it in the same
   amendment rather than leaving a second stale row for C4 to inherit.
3. **§C1.5, `@domain: none (<reason>)`** — conflates "no condition applies" (a defined
   constant) with "no condition is stated" (a defect). Split the vocabulary: keep `none` for
   the first, add `unstated (<reason>)` for the second, and put `unstated` tables on a
   worklist. Today the second silently exempts a table from D-032 (F-5).
4. **§C3 exit gate** — "plausibility suites green" does not name `test_consumer_domains`,
   which is how a red D-032 suite coexists with a gate the phase believes it passed (F-1).
   Name the suites the gate means, by module.
5. **§C1.2 citation grammar** — does not require the tag to begin the comment, so a tag can
   hide mid-sentence (F-6). Add "the tag begins the comment line" and let R1 enforce it;
   this also makes the register's published `grep` recipes correct by construction.
6. **Residual register format** — the `[UNVERIFIED]` class table should not close with a
   remainder bucket (F-8). A remainder cannot be falsified, and 20 residuals currently live
   in it un-itemised.
