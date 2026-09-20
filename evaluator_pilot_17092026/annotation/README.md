# Stage 3 — expert annotation of the frozen 300

The referee the pilot is missing. E0, E0-3J, E1, E2, E3, E4 and E5 have all been
measured, but every comparison so far is evaluators against each other, or against
E0's own final-answer check, which is itself defective (E0-F1, E0-F2). Expert
step-level labels are what decides which evaluator is right.

## The design

| | |
|---|---|
| Traces | the frozen 300 (five models × 60 items) |
| Labelled by | **three experts of the trace's own branch** — 5 branches, 60 traces each, 3 experts per branch, **15 experts** |
| Each expert | their branch's 60 traces + 8 shared calibration traces = 68, about 10 hours |
| Total | 1,020 trace-annotations, of which **900 form the ground truth** |
| Blinding | model identity hidden behind a code (`T-xxxxxx`); no evaluator output is ever shown |
| Order | shuffled per expert, so fatigue does not align across the three labelling the same trace |

Recorded per trace: a label for every **step** (correct / correct-but-different-route /
incorrect + error type / no-claim), a **milestone** checklist (obtained / not obtained /
not needed on this route), the **final answer** verdict, overall soundness, confidence,
comments and a discussion flag.

**Calibration.** The first 10 traces are shared by all 15 experts, 2 per branch. They
measure whether the branches apply the same standard. Labels an expert gives outside
their own branch are diagnostic only and never enter the ground truth (`for_truth`
marks which is which).

**Where an error starts.** Only the step where an error originates is labelled
incorrect; later steps that faithfully carry a wrong value forward are correct. This is
the single biggest source of annotator disagreement, so it is stated in the guide and
in the app.

## What the experts receive

1. **[EngTrace-annotation-guide.pdf](EngTrace-annotation-guide.pdf)** — 4 pages,
   generated from `guide.md`. It describes the task and the labels, and is deliberately
   **not written around the app**: it documents both ways of recording answers.
2. **Either the app or a workbook**, whichever the expert prefers:
   - the app shows one solution at a time and blocks incomplete submissions;
   - the workbook `tasks/workbooks/<id>.json` holds all 68 solutions with empty fields
     to fill in any editor, and is validated on return.

   Both produce identical rows in `labels/<id>.jsonl` (with `source` recording which),
   and an expert can switch between them: an export carries across whatever they have
   already submitted either way.

## Files

| File | What it does |
|---|---|
| `build_tasks.py` | Builds `tasks/` from the frozen slice and traces: `pool.json` (what experts see), `assignment.json` (who labels what, in what order), `keyfile.jsonl` (code → item + model; **not** given to annotators). Steps are split with E0's own `extract_steps`, so a label on step k refers to the same step in every evaluator. |
| `annotators.json` | The 15 experts and their branches. Edit as people are recruited, then rebuild. |
| `app.py` | The Streamlit annotation app. |
| `workbooks.py` | The file route: `export` writes a workbook per expert (carrying across any work already submitted), `check` validates a returned one and names every problem, `import` loads the valid solutions and skips the rest. |
| `guide.md`, `make_guide_pdf.py` | The annotator guide and its typesetter. |
| `score_against_labels.py` | Ground truth by majority, agreement (Fleiss κ), then scores every evaluator: **step level** (E2's PRMs), **milestone level** (E3/E4/E5), **trace level** (all of them, plus their final-answer checks). |

`tasks/`, `labels/` and `labels_simulated/` hold produced data and are not committed,
like `traces/` and `scores/`.

## Running it

```bash
cd evaluator_pilot_17092026
.venv/Scripts/python annotation/build_tasks.py          # after editing annotators.json
.venv/Scripts/python annotation/make_guide_pdf.py       # after editing guide.md
.venv/Scripts/streamlit run annotation/app.py           # the app the experts use
.venv/Scripts/python annotation/score_against_labels.py            # once labels exist
.venv/Scripts/python annotation/score_against_labels.py --simulate # pipeline test
```

`--simulate` writes synthetic labels derived from E2 and E3 and runs the whole
pipeline on them. It proves the scoring runs end to end before any expert starts; its
numbers are **not results** and the report says so on every run. Milestone scores come
out perfect under simulation precisely because the synthetic labels were drawn from
E3 — a reminder that the simulation tests plumbing, not correctness.

## Still to decide

- **Recruitment**: 15 experts, 3 per branch. Names go in `annotators.json`.
- **Where the app runs.** Locally per annotator, or one shared deployment. Shared
  needs a host and per-annotator access; local needs Python on each machine and the
  label files collected afterwards.
- **Adjudication.** What happens to steps with no majority, and to flagged traces:
  a review meeting, or one senior annotator deciding.
- **Compensation and credit** for the experts.
