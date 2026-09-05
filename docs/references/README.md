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
