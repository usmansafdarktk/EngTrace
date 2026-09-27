# full_run_28092026 — the full benchmark run

Started 2026-09-28. This directory holds the frozen item pool for the full run: 150 templates
× 15 instances, 2,250 items (D-114). The roster (D-110) and the evaluation stack (D-105, D-110,
D-113) are decided in DECISIONS.md. The inference harness will live here too.

## What is committed, and what is not

| File | Committed | Holds |
|---|---|---|
| `freeze.py` | yes | builds, verifies and checks the pool |
| `manifest.jsonl` | yes | one row per item: id, template, index, branch, domain, area, level, answer type, SHA-256 over question + NUL + solution, and `repeat_of` on a repeat. No text, no seed |
| `FREEZE.json` | yes | the counts, the selection rule, every replaced index and why, the short templates, the SHA-256 of the seed, the manifest's hash |
| `pool/` | **no** | the question and gold text, one `.jsonl` per module, in `testset/`'s schema plus `item_id`, `instance_index` and `repeat_of` |
| `SEED.secret` | **no** | the 128-bit master seed |

The templates are public, so publishing the seed or the text would publish the items. Both
stay local until publication. At publication the seed is revealed, and anyone can regenerate the
pool and check it against the manifest and the seed commitment. Reviewers get `pool/` in the
ARR supplementary archive and can check it with `--check-files`, which needs no seed.

**Back up `pool/` and `SEED.secret` privately.** Without the seed the pool cannot be
regenerated, and `freeze.py` refuses to draw a new seed while a manifest exists.

## Commands

```bash
python -m full_run_28092026.freeze --verify       # regenerate from the seed: must print VERIFY OK and FILES OK
python -m full_run_28092026.freeze --check-files  # pool/ against the manifest; no seed needed
python -m full_run_28092026.freeze                # rebuild. Only if a template changes, which voids the freeze
```

## The pool as frozen, 2026-09-28

| | |
|---|---|
| items | 2,250: 870 Easy, 870 Intermediate, 510 Advanced; 450 per branch |
| distinct questions | 2,235 |
| short templates | `adiabatic_flame_temperature`: 11 distinct questions + 4 repeats; `heat_of_reaction_formation`: 4 distinct + 11 repeats; 75 indices searched in each |
| replaced indices | 56: 51 repeated questions, 5 display ties (4 in `server_configuration_selection`, 1 in `critical_depth_froude_classification`); none for matching a pilot question |
| checks at the freeze, commit `1acc7cf` | Layer 0 gate 150 of 150; certification 150 of 150, current output identical to what the experts reviewed; `--verify` byte-identical in a separate process |
| manifest SHA-256 | `f0ffa108e79312e6ce8fea9a7fbab64957463c14bd7f3c7494460c9389032fab` |
| seed commitment | `cd53376fb0fe90ac807ff369a845a5a45e47c1bc48aa4fd470b709a191ddd64e` |

Every figure here is printed by `freeze.py` and recorded in `FREEZE.json`.

`testset/` is not this pool. It is the 2026-09-17 generation from the default seed, and the
pilot's `freeze.py --verify` rebuilds its slice from it, so it is left as it is.
