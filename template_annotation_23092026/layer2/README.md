# Layer 2 — human certification of the 150 templates

The stage the December 2025 record could not defend: three experts per template, own branch only,
with a hand check that leaves evidence, planted defects that give a detection rate, and timestamps.
Protocol in D-095; the reasoning in D-092.

## What each expert receives

One folder, `dist/kit_<id>.zip`, built by `make_kits.py`:

| File | What |
|---|---|
| `app.py` | the review app; `streamlit run app.py` from this folder |
| `tasks/` | ONLY their 34 items (the branch's 30 templates and 4 planted defects, shuffled), instances precomputed |
| `EngTrace-certification-guide.pdf`, `guide.md` | the guide (typeset from `guide.md` by `make_guide_pdf.py`) |
| `README.txt` | how to run, where the results go |

**Where results go:** the app writes `<id>.jsonl` into the kit folder, one row per submitted item,
saving after each submit; that file is what the expert sends back, and it is dropped into
`layer2/labels/` for `score.py`. The workbook route (`workbooks.py`) exists for an expert who cannot
run Python; its import produces the same rows.

No template code ships: five instances per item are precomputed. The keyfile that says which items
are plants never leaves this directory; the builder refuses a kit whose payload names a plant.

## The protocol

- **Own branch, three experts per template**: the pilot's roster of 15 (`evaluator_pilot_17092026/annotation/annotators.json`) by default; put a different roster in `annotators.json` here to override.
- **Hand check first**: the app shows one instance's question alone; the expert enters their answer; only then the solution appears, with a match/mismatch within 1% against the template's answer. The same instance (seed 2101) for all three experts of a branch.
- **Three 1–5 scores** (physical plausibility, mathematical correctness, pedagogical clarity), **Approve/Reject**, a defect type and a note required on Reject, a confidence score.
- **Planted defects**: four per branch, one of each class (constant, unit, sign, arithmetic), written as mutations of real templates under `plants/` to the contract in `plants/CONTRACT.md`, verified on every build. Experts are told quality-control items exist, not which.
- **Timestamps** for open, hand-check and submit, so dwell time is on record.
- **Blind**: opaque item codes, shuffled order per expert, no screen verdict shown, no discussion until all three have submitted.

## Order of operations

```bash
python -m template_annotation_23092026.layer2.build_tasks --check-plants   # verify the 20 plants
python -m template_annotation_23092026.layer2.build_tasks                  # tasks/: pool, assignment, keyfile
python -m template_annotation_23092026.layer2.make_guide_pdf               # the PDF
python -m template_annotation_23092026.layer2.simulate                     # proves the pipeline; NOT results
python -m template_annotation_23092026.layer2.make_kits                    # dist/kit_<id>.zip, one per expert: send these
# ... experts work; returned <id>.jsonl files go to labels/; a filled workbook goes through
python -m template_annotation_23092026.layer2.workbooks import <file>      # (workbook route)
python -m template_annotation_23092026.layer2.score                        # RESULTS.md
```

`tasks/`, `labels/`, `labels_simulated/` and `dist/` are produced data and git-ignored; the scripts,
the guide, the plants and `RESULTS.md` are committed.

## What comes back

`score.py` writes `RESULTS.md`: plant detection per expert and overall, hand-check agreement and the
approvals that followed a mismatch, agreement among the three experts (Fleiss κ, Gwet AC1 on
Approve/Reject; AC2 on the scores), the false-positive rate of the screening panel against the
experts and the MAD between their medians, dwell time, and the fix list of every rejected template
with the experts' notes. Rejected templates are fixed and, per D-094, re-judged individually by the
screen rather than in a new pass.
