# Prompt 01 — Template structure audit

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

---

I need a critical structural audit of every symbolic template in this repository. This is a
**feasibility study, not an implementation task.** Do not modify any template. By the end I
want to know whether a specific idea is viable, and what it would actually cost.

## Background you need

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems. It is built
from **150 parameterized Python templates** under `data/templates/branches/`, across five
branches (chemical, electrical, mechanical, civil, industrial — 30 each). Each template is a
function `template_*()` that samples physically-grounded parameters, computes an answer, and
returns a `(question, solution)` pair of natural-language strings. The solution string is the
gold reasoning trace, written as `**Step 1:** ... **Step 2:** ...` with a final `**Answer:**`.

The paper has been rejected twice from ACL ARR (Jan 2026 and May 2026 cycles). The most
damaging reviewer criticism, raised independently by five reviewers, concerns the evaluation
framework: it verifies model reasoning using an "AI Tribunal" of three frontier LLMs
(GPT-5, Claude Opus 4.5, Gemini 3), while simultaneously evaluating models from those same
families. GPT-5 is literally both a judge and an evaluated model. The paper is titled
*"Verifiable Process Supervision"*, but verification is a majority vote among LLMs.

Read these for the full picture:
- `docs/_ARR_May__EngTrace.txt` — the current paper (Section 4 is the evaluation framework;
  Appendix G has three full template implementations; Appendix M is the framework validation)
- `docs/EngTrace_Rebuttal_Jul2026.txt` — the most recent reviewer objections and our responses
- `evaluation/engtrace_evaluation_framework.py` — the current two-tier evaluator
- `evaluation/engineering_parser.py` — how model output is currently parsed into steps/values

## The hypothesis under test

The templates are executable Python that computes the answer. Every intermediate quantity —
every `R_A`, `v_z_max`, `F_a`, `sigma_AB` — exists as a typed value inside the generator, with
known units and a known governing formula. It is then f-string-serialized into prose, that
structure is discarded, and three frontier APIs are paid to approximately recover it.

**If templates instead emitted a structured trace** — a list of steps carrying
`{id, symbol, value, unit, formula, depends_on, is_milestone}` — then verifying a model's
reasoning could become a deterministic numeric and dimensional check against known milestone
quantities, with no LLM in the critical path. That would make the verification claim literal
rather than aspirational.

**Your job is to find out whether that is actually achievable across these 150 templates, or
whether it breaks down — and if so, exactly where and why.** I would rather learn now that
this is infeasible than after refactoring 150 templates. Be skeptical. Actively look for cases
that defeat the idea.

## What the audit must determine

Work through the templates systematically. Use subagents to parallelize across branches if
that helps, but you own the synthesis.

**1. Structural uniformity.** How are these functions actually written? Is there a consistent
shape (sample → compute → format), or do they vary enough that no mechanical transformation
applies? Note that the civil and industrial branches were authored later under a spec-driven
process with provenance-tagged constants — check whether they are structurally different from
the three original branches, and whether that makes them easier or harder to instrument.

**2. Are intermediate quantities bound to variables?** This is the crux. A step like

```python
f"sigma = ({load * 1000} N) / ({round(area, 4)} mm^2) = {round(stress, 3)} MPa\n"
```

has `stress` bound to a variable and recoverable. But quantities computed inline inside an
f-string expression and never assigned are not. Quantify how often each pattern occurs. If
most intermediate values are inline expressions, instrumentation means rewriting the
computation, not just adding emission — a very different cost.

**3. Are units recoverable?** Units currently exist only as literal strings in prose. Is there
a consistent convention? How do dual-unit templates (many branch on
`use_si_units = random.choice([True, False])`) carry units through? Could a unit be attached
to every emitted quantity without hand-annotating all 150 templates?

**4. Instance-dependent structure.** Some templates change the governing equation based on
sampled parameters — `template_reynolds_number_flow_regime` picks pipe vs. flat plate, which
changes the characteristic length definition *and* the critical thresholds. This means the
milestone set is a property of the **instance**, not the template. Confirm how widespread this
is and what it implies for a design that wants to declare milestones statically.

**5. Answer types that resist numeric milestone checking.** This is where I expect the idea to
struggle, so investigate it hardest. Catalogue every template whose answer or intermediate
steps are not plain scalars. From what I already know, these cases exist:
   - **Vectors** — electromagnetics cross products returning `x/y/z` components
   - **Arrays and numerical integration** — `template_levenspiel_plot_interpretation` builds a
     data table and integrates with `np.trapezoid`, wrapped in `try/except` fallbacks
   - **Symbolic expressions** — a BER template whose gold answer is `0.50 * Q(2.93)`, not a number
   - **Classifications** — "laminar" / "transitional" / "turbulent" as the answer
   - **Multi-part answers** — `template_statically_indeterminate` returns reactions *and*
     stresses, as parts a) and b)

   For each class: how many templates, and could milestone verification handle it, with what
   extension? Where it genuinely cannot, say so plainly.

**6. Steps vs. milestones.** Do the `**Step N:**` markers correspond one-to-one with computed
quantities? Are there narrative steps that state a principle without computing anything? A
milestone-based scheme needs to know which steps carry verifiable content.

**7. Precision and rounding.** Solutions print rounded values while the internal variable holds
full precision, and rounding is inconsistent (`precision = 3`, `round(x, 4)`, `:.3e`). Which
value is the milestone — printed or true? What does that imply for a numeric tolerance?

**8. The extraction side — do not skip this.** Instrumenting the gold trace only helps if the
corresponding quantities can be recovered from the *model's* prose output, which stays
unstructured. Read `evaluation/engineering_parser.py` and assess honestly whether its
regex-based value extraction is good enough to support deterministic matching. If you can
find archived model outputs, test it against real traces. **If extraction is the weak link,
that is the single most important finding in this audit** — it would mean deterministic
verification fails for reasons that have nothing to do with the templates.

## Deliverables

Write both to `docs/re-implementation-sep/`:

**1. `template_audit_report.md`** — the written findings. Lead with a direct verdict on
feasibility and the evidence for it. Cover each numbered question above. Include concrete code
excerpts for the patterns you are describing, especially the hard cases. End with a
per-template-class cost estimate and an explicit recommendation: proceed, proceed with a
modified design, or abandon.

**2. `template_inventory.csv`** — one row per template, machine-readable, with at least:

```
template_id, file_path, branch, domain, area, difficulty,
n_steps, n_milestone_candidates, has_instance_branching,
answer_type,          # scalar | vector | array | symbolic | classification | multipart
unit_system,          # SI | US | dual | dimensionless
values_bound,         # yes | partial | no
instrumentation_class,# A=trivial, B=branching, C=non-scalar, D=resists the approach
est_effort,           # low | medium | high
notes
```

The classification matters more than the prose — it is what the next stages plan against, so
apply it consistently and define each class explicitly in the report.

## Constraints

- **Read-only.** Do not modify, refactor, or "fix" any template. If you spot bugs, record them
  in the report as findings.
- Stay on `master`. Do not create branches or commits.
- `pilot/` is gitignored working material from the civil/industrial authoring effort — you may
  read it for context on how those templates were built, but audit only what is committed
  under `data/templates/`.
- Templates import from `data.templates.branches.<branch>.constants`; run them from the repo
  root if you need to generate instances. Generating instances to inspect real output is
  encouraged — it is often faster than reading the formatting code.
- Where you are uncertain, say so and explain what evidence would settle it. Do not smooth over
  gaps to make the recommendation cleaner.

## What I do not want

No implementation. No proposed instrumentation API. No refactoring plan beyond the cost
estimate. Design comes in a later session, deliberately after this one, and it should be
informed by what you find rather than by what I assumed when writing this prompt.

If the audit shows the hypothesis is wrong, that is a genuinely useful outcome — say so
clearly and explain what killed it.
