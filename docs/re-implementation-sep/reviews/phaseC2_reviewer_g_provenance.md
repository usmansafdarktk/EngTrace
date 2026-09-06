# Phase C2 — Reviewer G: Provenance and units

**Filed:** 2026-09-06 · **Branch reviewed:** `redesign/phaseC2-provenance-completion`
**Roster:** G — Provenance and units, phases C1–C3
**Brief:** one mandatory gate task — *resolve every citation in the three
re-grounded tables against the artefact it names, and report whether the
provenance claim is honest and complete.* 40-minute box.

> Worked without sight of the implementer's reasoning (R1.1). Code comments,
> commit messages and `phaseC2_summary.md` were supplied as **claims under
> test**. Value re-derivation was explicitly fenced off as Reviewer H's
> territory (R6.2): G was asked whether the *record* is honest, not whether the
> *numbers* are right — which is why G could find that three values had no
> backing without ever needing to know they were correct.

**Remediation is triaged in [`../phaseC2_summary.md`](../phaseC2_summary.md) §12
and recorded as [D-034](../DECISIONS.md) and [D-035](../DECISIONS.md).**

---

## 1. Verdict

**PASS WITH FINDINGS** — the re-grounding is real and most citations resolve,
but the record overstates itself in three specific ways, three `[ON-DISK]`
citations resolve to nothing, and spec C1.4 (unit declarations) was skipped
entirely, which is an unmet exit-gate item.

G filed once, then **filed a revision** after three follow-up checks: one
retraction, two corrections, and a root cause that made G-1 and G-2 into a
single deeper finding. Both the retraction and the revision are recorded below,
because a reviewer that corrects itself under scrutiny is doing the job.

---

## 2. Method

All 54 citations in `CP_PARAMS` (33) and `HEATS_OF_FORMATION` (21) resolved
mechanically: parse the CAS out of each row's tag, look up the species, and
check the specific field the citation claims.

| Test | Result |
|---|---|
| Every `CP_PARAMS` CAS resolves to the right species | **33/33** — no dangling CAS, no wrong species |
| Every `HEATS_OF_FORMATION` citation resolves to the phase-correct dHf field | **3 fail** (G-1) |
| Header claim "every one is inside its stated uncertainty" | **2 breach** (G-2) |
| Tag class matches evidential state | **4 rows overclaim** (G-8) |
| Unit declaration present | **absent everywhere** (G-4) |
| `[KNOWN-DEFECTIVE]` absence claims | **both corroborated** against live NIST |
| Primary source (Smith–Van Ness) on disk | **absent** (G-7) |

---

## 3. Findings

| # | Finding | Status |
|---|---|---|
| **G-1** | **CONFIRMED, high.** Three `[ON-DISK]` citations resolve to nothing: `C8H18(g)` (CAS 111-65-9 absent from the file entirely), `CH3OH(g)` (only `CH3OH(l)` present, no gas dHf), `C2H5OH(l)` (only `C2H5OH(g)`, no liquid dHf). "*A reviewer told to check a value offline cannot; the `[ON-DISK]` guarantee is the whole point of the tag class, so three counterexamples devalue all 21.*" | **FIXED** — all three species fetched live from NIST and written to the artefact; every one of the 21 now resolves |
| **G-2** | **CONFIRMED, high.** The claim "all 21 verified, none outside its stated uncertainty" is **false**: `C2H6(g)` −84.7 against a cited −84.0 ± 0.4 (1.75×), `CH3OH(l)` −238.6 against a cited −239.5 ± 0.2 (4.5×). | **FIXED** — and the diagnosis was wrong in an instructive way: the *values* were right, the *citations* named the wrong measurement. See §5 |
| **G-3** | **CONFIRMED, moderate.** The `C3H8(g)` row asserts NIST lists two values; the artefact recorded only −104.7, with no `dHf_alternative` field, though that field exists and is populated for `CH4(g)`. "*The one row where a deliberate source choice was made is the row whose choice is unauditable offline.*" | **FIXED** — both measurements on disk with method and reference |
| **G-4** | **CONFIRMED, moderate.** Spec C1.4 requires a unit declaration per table. There is none anywhere; units live in prose only. "*Stating this plainly as instructed: it was not done.*" | **FIXED** — units declared for `HEATS_OF_FORMATION` (kJ/mol at 298.15 K, 1 bar), `CP_PARAMS` (dimensionless, returns Cp/R), `CP_VALID_T_MAX` (K) and `AIR_COMPOSITION` (mole fractions, and that they sum to 0.99964, not 1) |
| **G-5** | **STALE.** "`CP_PARAMS` carries no `[ON-DISK]` tag at all." True when G started; the tag normalisation landed while G was running. Current state: 47 `[ON-DISK]` + 4 `[DERIVED]` + 3 `[BY-DEFINITION]` = 54. | **Already fixed before filing** — recorded, not silently dropped |
| **G-6** | **CONFIRMED, low–moderate.** `docs/references/README.md` documents a citation format (`[ON-DISK] NIST-WEBBOOK CH4 (74-82-8)`) and a file layout (one `.md` per CAS) that never existed. The README was written before the data was acquired and never reconciled with it. | **FIXED** — README now describes the actual JSON layout, the actual tag format, the four tag classes, and what the test enforces |
| **G-7** | **CONFIRMED, moderate — recorded, not closed.** The coefficients are Smith–Van Ness; the only on-disk artefact is NIST. "*The record documents verification, not origin.*" A reader can confirm a value is plausible offline but not that it was transcribed correctly — which is exactly the Class-1 defect C2 exists to fix. | **RECORDED** — already D-030; the README now says so in as many words. No citable copy of SVN Table C.1 exists to fix it with |
| **G-8** | **CONFIRMED, moderate.** Four rows tagged stronger than their evidence. `C2H6(g)` recorded a Smith–Van Ness self-check while a 5-point Gurvich table sat unused on disk; `C6H6(g)`/`C7H8(g)` cited *liquid* values on *gas* rows; `C3H8(g)`/`C4H10(g)` license 1500 K from a single 298 K point. | **MOSTLY FIXED** — the suite now consults all six multi-point tables, not two (**215 → 222 checks**); the single-point ceilings are **reported, not fixed** (D-035) |
| **G-9** | **CONFIRMED (revised), low.** `CP_VALID_T_MAX` and `AIR_COMPOSITION` are unsourced. G first called them dead code, then **corrected itself**: they are consumed, but by tests only — reported at test time, not enforced at generation time. | **FIXED** — both sourced; both now state plainly that they are advisory and that the test suite is their only consumer |
| **G-10** | **STALE.** "The defective register documents but does not quarantine — `H2SO4(l)` and `CaCO3(s)` remain live at 2/33 draw probability." Both species were **removed** under D-033 while G was running; they survive only in `REPLACES` history comments. | **Already fixed before filing** |
| **G-11** | **STALE.** The `H2SO4` note overclaims its scope ("paid-linked sections verified"). The row it describes no longer exists. | **Moot via D-033** |
| **G-12** | **CONFIRMED (revised), blocking for C3.** No test asserts that a citation resolves. G revised this too: a real check *did* exist for `CP_PARAMS` — but it resolves by **species key, not the CAS in the tag**, and `HEATS_OF_FORMATION` had no equivalent at all. "*That asymmetry is the root cause of G-1.*" | **FIXED** — `tests/constants_integrity/test_citations_resolve.py`, 55 checks |

---

## 4. G's own corrections

Recorded because the corrections matter as much as the findings.

**RETRACTED — a near-finding on `CP_NO_REFERENCE`.** G read the test file with
a spliced `sed -n '30,60p;76,82p'`, and the non-contiguous ranges created a
false adjacency that made three `(value, uncertainty)` tuples appear to sit
inside the Cp-exemption dict. They belonged to a different dict entirely. G
caught it when a verification script raised `KeyError`, and flagged the
near-miss unprompted: *"a spliced read is exactly how a reviewer manufactures a
false positive."*

**G-9 corrected** from "dead code" to "consumed by tests only" — G's original
grep was scoped to `data/` and missed the test suite's import.

**G-12 corrected** from "nothing checks citations" to "one table is checked, by
key rather than by CAS, and the other is not checked at all."

---

## 5. The root cause, which G found only on the second pass

This is the most valuable thing in the review, and neither the implementer nor
Reviewer H found it.

`HEATS_OF_FORMATION`'s reference values lived as a hardcoded `DHF_REF` dict
**inside the test file**. `CP_PARAMS` was checked against the on-disk artefact;
heats of formation were checked against that dict. Consequences:

1. **The three dangling `[ON-DISK]` citations are exactly the three rows whose
   reference value existed only inside the test.** For octane the check was
   circular — the dict said −208.4, the constant said −208.4, deviation 0.000,
   and no artefact backed either. The suite reported them as verified.
2. **`DHF_TOL_FLOOR = 2.0` floored every tolerance**, so rows breaching their own
   cited uncertainty passed anyway.
3. **`C2H6(g)`'s uncertainty had been widened** from NIST's 0.4 to 0.7 in that
   dict — exactly enough to cover its deviation.

G's summary is exact: *"the suite tests 'inside max(stated, 2.0)' while the
header claims 'none outside its stated uncertainty.' The claim and the check are
different propositions, and only the weaker one was tested."*

**What the re-check found.** Every disputed value turned out to be right, and
matched a *specific* NIST measurement once all of them were on disk — the repo
follows **Prosen & Rossini (1945)** for hydrocarbons, an undocumented but
defensible source choice that read as three errors. One genuine defect fell out:
`NO2(g)` read 33.2 against the single value NIST publishes, 33.10. Corrected to
33.1; it feeds no reaction and no template, so the item pool is unchanged.

Full table in [D-034](../DECISIONS.md).

---

## 6. Falsification attempts that failed

- **The two `[KNOWN-DEFECTIVE]` absence claims held.** G fetched live NIST for
  calcite and sulfuric acid; neither page carries condensed-phase Cp. *"The
  register is not being used to retire work that could have been done."* One
  reservation: "citable" was narrowed to NIST WebBook without saying so — calcite
  Cp is in NIST-JANAF, the artefact's own declared primary reference.
- **No wrong-species citation in `CP_PARAMS`.** All 33 CAS resolve correctly,
  including the four phase-ambiguous CAS that each map to two entries.
- **The four `(refit)` tags are honest** — each backed by a genuine multi-point
  table, each declaring a range matching `CP_VALID_T_MAX`. `C2H2(g)` correctly
  derates to 1100 K rather than taking the blanket 1500 K: *"the one place
  validity was reasoned about rather than defaulted."*
- **The "agreement is evidence, not a tautology" argument holds** for every row
  with a Shomate range. It failed only for `C2H6(g)` (G-8) — now fixed.

---

## 7. Does this generalise to C3?

**"Yes, and worse."** C3 is ~400 values across four branches, and its review line
already presumes the outcome: *"Reviewer G confirms every table now carries a tag
and a unit declaration."*

> At 54 values, three dangling citations survived merge, two reviewers and a
> 215-check suite — because nothing machine-checks citation resolution. At 400
> values, hand-checking is not viable and the same failure mode scales linearly.

G's recommendation — a `test_citations_resolve.py` as a **hard entry condition**
for C3 — is the one finding G called blocking. It exists now, and G's revision
added the two properties that make it bite: resolve the **CAS written in the
tag**, not just the species key; and **read reference values from the artefact**,
never from constants embedded in the test. Both failure modes are demonstrated
by C2 itself.

**What the spec got wrong.** (1) C1.4 specifies a unit declaration but the exit
gate has no test for it, so it was skipped without anyone noticing — *"an
unenforced deliverable is a suggestion."* (2) C2.1 demands citations "to edition
and page" while C1.1 acquired a CAS-keyed JSON; the deliverable shape and the
citation shape were never reconciled. (3) The spec's two-tag vocabulary has no
slot for the states C2 actually met — *checked and failed*, *exact by
definition*, *derived*. All three were invented correctly but off-spec. **The
vocabulary should be fixed in the spec before C3 tags 400 more values against
it** — raised as a SPEC-CHANGE, not applied (D-035).

---

*Note on §3: G-5, G-10 and G-11 are recorded as STALE rather than dropped. All
three were true when G began and were fixed by the species swap and tag
normalisation while G was running — a hazard of reviewing a moving branch, and
the implementer's fault for moving it, not the reviewer's for reporting it.*
