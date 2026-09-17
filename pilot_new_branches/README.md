# pilot_new_branches/ — staging from building the two new branches

**This is not the evaluation pilot.** The evaluation pilot lives in
[`evaluator_pilot_17092026/`](../evaluator_pilot_17092026/). This directory was called `pilot/` until
2026-09-17, which made the two mean the same word, and that is the confusion this
rename exists to end.

## What it is

Working material from September 2026, when the **civil** and **industrial**
branches were built and added to EngTrace. It holds two branches and not five for
the plain reason that those were the two being built; chemical, electrical and
mechanical already existed and were never staged here.

| Directory          | What it held                                                        | Superseded by |
|--------------------|---------------------------------------------------------------------|---------------|
| `branches/`        | Constants and template drafts for civil and industrial              | [`data/templates/branches/`](../data/templates/branches/) — all five branches |
| `templates/`       | A second staging copy of the same two branch trees                  | same |
| `references/`      | Source PDFs, extracted text and tables the constants were grounded in | [`docs/references/`](../docs/references/) with `MANIFEST.json` recording each URL and SHA-256 |
| `testset_preview/` | An early industrial-only generation, used to eyeball the templates  | `generate_testset.py`, which builds all 2,250 records |

## Why it is kept

The reference material is the provenance trail for the physical constants in the
two new branches. `docs/references/MANIFEST.json` re-acquires and hashes the public
sources, but the extracted text and the full textbooks here are what the constants
were actually read out of, and re-downloading them is slow. Nothing in the live
tree imports from this directory.

## Why it is not committed

Large binaries, and duplicates of material that is now properly placed. `.gitignore`
holds back everything here except this file, so a clone shows the explanation
without the payload.

## If you are looking for

- **the branch templates** → [`data/templates/branches/`](../data/templates/branches/)
- **the constants' sources** → [`docs/references/`](../docs/references/)
- **the record of the rework** → [`docs/re-implementation-sep/`](../docs/re-implementation-sep/)
- **the evaluation pilot** → [`evaluator_pilot_17092026/`](../evaluator_pilot_17092026/)

Safe to delete once the reference PDFs are no longer wanted locally. Deleting it
loses no committed history, because none of it is committed.
