# On-disk reference sources

Acquired for provenance tagging (spec Phase C1 deliverable C1.1). Every
`[ON-DISK]` citation in `data/templates/branches/*/constants.py` resolves to a
file here, so a reviewer can check a value without network access and without a
library.

| Directory | Source | Used by |
|---|---|---|
| `nist_webbook/` | NIST Chemistry WebBook, SRD 69 — species pages, retrieved 2026-09-06 | Phase C2: `CP_PARAMS`, `HEATS_OF_FORMATION` |

**Citation format used in the constants files:**

    [ON-DISK] <source key> <locator>   e.g.  [ON-DISK] NIST-WEBBOOK CH4 (74-82-8), Gas phase thermochemistry

Each retrieved page is stored as a `.md` file named by CAS registry number, so
the citation is resolvable to a file on disk.

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
