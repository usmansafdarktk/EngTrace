# On-disk reference sources

Acquired for provenance tagging (spec Phase C1 deliverable C1.1). Every
`[ON-DISK]` citation in `data/templates/branches/*/constants.py` resolves to a
file here, so a reviewer can check a value without network access and without a
library.

| Directory | Source | Used by |
|---|---|---|
| `nist_webbook/` | NIST Chemistry WebBook, SRD 69 — species pages, retrieved 2026-09-06 | Phase C2: `CP_PARAMS`, `HEATS_OF_FORMATION` |

**How the data is stored.** One machine-readable file per source directory,
not one page per species: `nist_webbook/shomate_coefficients.json`, keyed by
the same species strings the constants tables use, each entry carrying its CAS
number, the Shomate ranges, any tabulated Cp points, and every heat of
formation NIST lists for that species. A per-species `.md` layout was described
here before any data was acquired; it was never what got written (Phase C2
Reviewer G, finding G-6).

**Citation format used in the constants files:**

    [<class>] NIST <CAS> <locator and the check that was run>

    e.g.  # [ON-DISK] NIST 74-82-8  Cp298 35.06 vs NIST 35.65

Four classes, because a value can be warranted in four different ways and
collapsing them hides which:

| Class | Means | Requires |
|---|---|---|
| `[ON-DISK]` | checked against retrieved data | an entry in the named file whose CAS matches the tag |
| `[DERIVED]` | fitted or computed here from on-disk data | a stated derivation and range |
| `[BY-DEFINITION]` | true by how the scale is defined, e.g. an element in its standard state | a stated definition; no artefact can be its warrant |
| `[KNOWN-DEFECTIVE]` | checked and *failed*; kept only as a record | the measured error |

`tests/constants_integrity/test_citations_resolve.py` enforces all of this: it
parses every tag, resolves the CAS **written in the tag** against this
directory, and fails on a tagless row, an unresolvable CAS, or a CAS that
resolves to a different species than the row is keyed for. It also refuses to
let a test file carry its own copy of the reference values — the mechanism by
which three unbacked citations survived to merge in Phase C2 (finding G-1).

**What is NOT here.** The `CP_PARAMS` coefficients are in the Smith-Van Ness
functional form, and no fetchable, citable copy of Smith-Van Ness Table C.1
could be obtained (D-030). So this directory documents *verification*, not
*origin*: it shows each value agrees with an independent NIST fit, which is
strong evidence the value is right and no evidence about where it was typed
from. Reviewer G's finding G-7, recorded rather than closed.

## Flame-temperature reference

`Adiabatic Flame Temperatures for Oxy-Methane, Oxy-Hydrogen, Air-Methane, and
Air-Hydrogen Stoichiometric Combustion using the NASA CEARUN Tool, GRI-Mech 3.0
Reaction Mechanism, and Cantera` — ETASR / arXiv:2503.11826, retrieved
2026-09-06.

Air–methane, stoichiometric, reactants at 298.15 K, 1 atm:

| model | T_ad |
|---|---|
| **complete combustion, no dissociation** | **2326.35 K** |
| chemical equilibrium, GRI-Mech 3.0 (dissociation included) | 2224.25 K |

`template_adiabatic_flame_temperature` uses a single balanced reaction with no
dissociation, so **2326 K is the correct reference for it**, not the ~2200 K
figure the redesign spec quotes. The ~100 K gap between the two is the
dissociation effect the item's model deliberately excludes.
