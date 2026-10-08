# WS-C5: Statistics tables and the worked example (Phase 1, parallel split of WS-C)

Read first: `WS-C_analysis_code.md` (context, rules, amendments; its steps C8 and C9) and `00_ORCHESTRATION.md`
section 7. This brief adds only what is specific to C5. Target: two to three hours.

## Mission

Two independent pieces: the templates-per-area and difficulty-agreement tables as generated blocks of the
statistics script, and the worked scoring example generated from a real response.

## Files you own

`docs/appendix_statistics.py`; new `full_run_28092026/worked_example.py` and its output
`appendices/worked_example.tex` written only under `--out DIR` in Phase 1. Nothing else.

## Steps

1. `docs/appendix_statistics.py --write [--out DIR]` fills marked blocks with the BEGIN/END convention of
   `paper_results.py`: `tab:area` (templates per area and level from `manifest.jsonl`'s `area` field, grouped by
   domain, with branch subtotals) and `tab:levels_agreement` (from
   `template_annotation_23092026/levels/agreement.json` when present: kappa per branch and overall with intervals,
   the weighted forms, templates where the majority differs from the label by one or two steps, the
   formula-stated count by level; until the file exists, a stand-in marked "STAND-IN" in the caption). `--check`
   keeps passing on the current tree; the single-path counts are C1's, not yours.
2. `worked_example.py --out DIR`: pick one real response from the main store that shows every mechanism (a
   lower-tier wrong answer with at least one milestone matched deterministically, one judge ruling of each kind if
   possible, and an arithmetic flag), and write `appendices/worked_example.tex` as a `tcolorbox`: the question in
   brief, the gold milestones with values, which matched deterministically (with the unit factor), which the judge
   ruled reached, not needed or missing, the arithmetic flag with operands and recomputed value, then v(r) and
   m(r). Every value from the store rows (`scores/main/<model>.jsonl`, `scores/main/e5/`, `scores/main/router/`);
   one column of a page at most; the item id and model in a LaTeX comment, not in the box.
3. A `--selftest` for each: the area table sums to 150 templates and 30 per branch; the worked example's v(r) and
   m(r) equal the stored `score` and `e5_strict` for the chosen response.

## Acceptance

Both scripts run with `--out <scratch>` and write their blocks; `docs/appendix_statistics.py --check` passes on
the untouched tree. Report: the item chosen for the example and why, the two block labels, the Phase 2 command
lines.
