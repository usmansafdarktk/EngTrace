# template_annotation_23092026 — certifying the 150 templates, again

Started 2026-09-23. The record of the template certification pipeline for the next submission:
three layers, ordered by cost and by what each can actually detect. The analysis that led here
is summarised in [D-092](../docs/re-implementation-sep/DECISIONS.md) and
[D-093](../docs/re-implementation-sep/DECISIONS.md); the December 2025 pipeline it replaces is
described in the paper's section 3.3 and Appendix K.

| Layer | What | Cost | State |
|---|---|---|---|
| 0 · deterministic gate | T1 closure, T3 determinism, T4 contract, T8 emission at 500 seeds; T5, T7 advisory; a register of four line classes the check cannot read | free | **green: 150 of 150 pass** (`layer0/gate_report.md`); 45 templates edited (`closure_fixes.md`), 44 moved instances (`item_pool_impact.md`); residual ties in 22 templates for a third round (`tie_census.md`) |
| 1 · LLM screen | Appendix H's prompt, verbatim; three non-suite judges; two passes, hard-capped | ~$4.60 a pass | tooling built and smoke-tested; pass 1 awaits approval |
| 2 · human certification | own-branch experts, three per template, with hand-checks and planted defects | expert time | not started; protocol to be written |

## Layout

```
template_annotation_23092026/
  README.md              this file
  layer0/
    gate.py              runs the gate corpus-wide, applies the register, writes the report
    check_limits.json    the register: template + line pattern + reason for every excused line
    gate_report.md/.json the gate's outcome at HEAD (regenerate with gate.py)
    tie_census.py        exact-decimal census of half-way display ties T1 cannot count reliably
    tie_census.md/.json  its output (informational, gates nothing)
    closure_fixes.md     every template edited for closure: cause, fix, what moved, sign-offs
    item_pool_impact.md  what the edits moved, measured by c3_instance_dump.py --diff
  screen/
    run_screen.py        the screening runner (--check, --dry-run, --pass N, --status, --smoke)
    analyze_screen.py    consensus (paper rule), sigma-max, Gwet AC1 / Fleiss kappa, cost, served ids
    pass1/, pass2/       config.json + replies.jsonl + summary.csv + flagged.md + stats.md (committed)
```

## Order of operations

```bash
python -m template_annotation_23092026.layer0.gate                 # free; must be green before anything else
python -m template_annotation_23092026.layer0.tie_census           # free; residual ties, informational
python -m template_annotation_23092026.screen.run_screen --check   # ~$0.01: keys, model ids, live prices
python -m template_annotation_23092026.screen.run_screen --pass 1 --dry-run   # free: prompts, tokens, estimate
python -m template_annotation_23092026.screen.run_screen --pass 1  # PAID (~$4.60); needs approval; resumable
python -m template_annotation_23092026.screen.analyze_screen --pass 1
```

Pass 2 is the same commands with `--pass 2`, after human certification. `--pass 3` is refused.

## Conventions

- **Every number in a report here comes from a committed script**; the agents' per-template
  measurements in `closure_fixes.md` are labelled as such and the corpus-wide figures are the
  scripts'.
- **The register is the residual list.** An excused line names its template, its pattern and the
  reason; a growing register is a warning and the gate report prints it in full.
- **A pass is a pass.** Screening rows resume only within their own pass directory, and a pass
  refuses to continue if a template's prompt hash has changed since it started.
- **No paid call without the user's approval.** Dry runs and checks are free; the runner says
  which mode bills.
