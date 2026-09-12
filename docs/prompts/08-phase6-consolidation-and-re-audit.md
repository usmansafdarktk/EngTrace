# Prompt 08 — Phase 6: consolidation and re-audit

Paste everything below the line into a fresh Claude Code session started at the repository
root, on the `master` branch.

Track A Phases 0–5 and Track B Phases C1–C3 are complete and merged. `master` is at
`444b8bf`. **Sync point S2 has landed**: the constants are re-grounded, so the item pool can
now be regenerated exactly once. That is what this phase is for.

**The repo owner has decided that the testset is rebuilt, inference re-run and evaluation
re-run.** That closes D-003 — the open decision this series has carried since 2026-09-05 — and
adds **D6.12**, a deliverable the spec does not have. It also makes one ordering rule binding
rather than advisory: **regenerate the testset only after D6.7 and D6.11 merge.** Both change
emitted text. Build the pool before them and 11 models are inferred against items that are stale
on arrival. Read §"D-003 is closed by decision" before planning the phase — the generator is in
worse shape than "absent", and two of its five branches have never had one.

## Decisions the repo owner has settled

These were open when the brief was first written. They are now closed, and the rest of this
document is written as though they are.

| question | decision |
|---|---|
| Does the regenerated benchmark include civil and industrial? | **Yes. 150 templates, not 90.** |
| Instances per template | **15, as before.** So the pool is **150 × 15 = 2,250 items**, up from 1,350. |
| Track B's five substance defects | **Resolve the easiest way: delete the rows.** See §"Track B — the seven deletions". |
| The stale annotation fork | **Leave it.** Not relevant; do not sync it. Name it in D6.6 and move on. |
| D6.1's missing audit harness | **Best effort.** Look for the production method; if it cannot be found, **document that briefly and state the method you used instead.** Do not spend the phase on archaeology. |
| D6.9 and D6.10 | **Attempt the fixes.** Record as decided-not-fixed **only if proven infeasible**, with the measurement that proves it. |

**Why deleting the Track B rows is cheap now, when it was not before.** Those five defects were
deferred across two phases because removing a row changes which substance each seed draws, which
invalidates the published item pool — a P6 event. **The pool is being rebuilt from scratch**, so
that cost is now zero. The blocker was never the edit; it was the consequence, and the
consequence has gone.

### Track B — the seven deletions

| table | delete | leaves |
|---|---|---:|
| `MANOMETER_FLUIDS` (mechanical) | `Tungsten Hexafluoride` (a **gas** at manometer conditions, value unsupported on disk), `Tellurium Mercury` (both WebBook pages are data-free stubs), and **one** of the `Tetrabromoethane` / `Acetylene Tetrabromide` pair (one substance under two names, both 2960) | 17 of 20 |
| `MATERIAL_DENSITIES` (mechanical) | `Cork`, `Cork Board`, `Bamboo` — searched exhaustively across all 546 pages of FPL-GTR-282 and every reference directory; genuinely absent | 53 of 56 |
| `THERMO_SUBSTANCES` **and** `REAL_FLUID_DATA` (chemical) | `Refrigerant-410A` — a blend whose components R-32 and R-125 are also absent, so it cannot be reconstructed | 24 of 25 each |

**R-410A is the one with a catch.** It is a **list element** in `THERMO_SUBSTANCES` and a **dict
key** in `REAL_FLUID_DATA`. Delete it from both, in the same commit — remove one and the list
names a fluid nothing supplies data for.

**Deleting is still a measured change.** D-031's rule holds: removing a key shifts the position
of every key after it, so a seed draws a different substance. Under full regeneration that is
intended rather than a defect, but **state it** — the P6 measurement for this tranche is part of
D6.4, not something the regeneration excuses.

**Re-derive every number in this brief at whatever `master` is when you start.** Every prior
brief in this series was corrected by what its own phase measured, and this one already
contains two corrections to the spec it implements (§"Two spec numbers that do not survive
contact").

---

I need you to implement **Phase 6** (consolidation and re-audit) of the template redesign —
the last phase of Track A.

## What EngTrace is, and why this phase matters now

EngTrace is a benchmark for evaluating LLM reasoning on engineering problems, built from **150
parameterised Python templates** under `data/templates/branches/`, across five branches
(chemical, electrical, mechanical, civil, industrial — 30 each). Each `template_*()` samples
physically-grounded parameters, computes an answer, and returns a `(question, solution)` pair;
the solution is the **gold reasoning trace**.

**The paper has been rejected twice from ACL ARR** (January and May 2026). Six phases have
changed the corpus underneath the published numbers: Phase 1 corrected gold traces that did not
reproduce their own answers, Phase 2 seeded every generator, Phase 3 reshaped two iterative
traces, Phase 4 specified the comparator, Phase 5 fixed the output contract and bound 120
templates, and Track B re-grounded the constants and corrected 21 measured-wrong values.

**Nothing has yet gone back and said what the corpus now is.** The audit report and the
inventory that classified all 150 templates still describe the repository as it was before
Phase 0. That is the gap this phase closes, and it is larger than "re-run and diff" — see the
next section.

Two cautions so you do not waste time:

- `error_analysis_annotation/error_annotation_results/` contains "Annotator A/B/C" — an **LLM
  annotation pilot**, not the paper's human error analysis. Not a discrepancy; do not flag it.
- `evaluation/` is a **separate track** with known defects (D-003). **Do not fix them.** But do
  read §"D-003 is answerable now" — this phase can finally close it, and the answer is not the
  one the spec assumed.

## The rule this phase runs on — re-measure, do not inherit

Phase 6's deliverable is a **description of the corpus**. A description assembled from six
phases' summaries is a description of what each phase *claimed*, and this series has repeatedly
found claims that its own artefacts do not carry.

- **Every figure in `template_audit_report.md` and `template_inventory.csv` is re-measured
  here, not copied forward.** A class migration inferred from a summary is not a measurement.
- **A prior phase's summary is a claim under test** (R1.1), including the ones I wrote. Two
  numbers in the spec's own Phase 6 section do not survive being checked — see below.
- **Every count ships with its predicate or a committed script.** Carried from Prompt 06 and
  unchanged. The unit census below is the case in point: the "113 templates carry a unit"
  figure is one of **four** defensible counts depending on the predicate, and the spec cites it
  without one.
- **Do not regenerate the T6 baseline to make a gate green** (D-043). T6 fails 142/150 on
  `master` because the committed baseline profile is stale corpus-wide. `run.py --baseline`
  must never be invoked in this phase. Fixing T6 properly is in scope *as a decision*; making
  it green by overwriting the oracle is not.

## What Phase 6 is, and what it is not

| | in scope | out of scope |
|---|---|---|
| D6.1–D6.6 | re-audit, re-classify, consolidate the item-pool statement, stand up CI, file the residual register | re-opening any merged phase's fixes |
| D6.7 | the doubled-sign residual on 14 templates, **behind a shared emission helper** | any other emission cleanup |
| D6.8 | an answer-span **shape** assertion | rewriting spans |
| D6.9–D6.10 | **attempt the fix** — a mixed `AnswerSpec`, and the 8 undecidable `symbolic` templates; record as unfixed only with the measurement that proves it | forcing either binding green, or weakening the comparator to raise a rate |
| D6.11 | per-item and per-part unit declarations | a second unit scheme — C1.4 already built one |
| **D6.12** (new) | rebuild the testset generator — seeded, uniform, all five branches — and regenerate **after** D6.7/D6.11 land | running inference or evaluation; those follow this phase |

**D6.9 and D6.10 are decisions, not fixes.** Both are measured dead ends: the mixed-spec
templates produce 544 and 1,650 false accepts in 2,450 pairs from *whole-span cross-pairing
working correctly*, and masking function names in `_to_sympy` moved the symbolic decided rate
by exactly nothing — 11.1% before and after. Re-attempting the known-failed fix is the failure
mode here. Record what it would take and move on.

## Read these first, in this order

1. **[`template_redesign_spec.md`](../re-implementation-sep/template_redesign_spec.md)** §"Phase 6"
   (D6.1–D6.11, the exit gate) and §"Track B" S2. The authority on scope — but read it against
   §"Two spec numbers that do not survive contact" below.
2. **[`template_audit_report.md`](../re-implementation-sep/template_audit_report.md)** and
   **`template_inventory.csv`** — what you are replacing. Note the date.
3. **[`phase5_summary.md`](../re-implementation-sep/phase5_summary.md)** §7 and §13 — the
   nearest phase's errors, and §R5-2, the origin of the unit census.
4. **[`phaseC3_summary.md`](../re-implementation-sep/phaseC3_summary.md)** and
   **[`phaseC3_corrections_summary.md`](../re-implementation-sep/phaseC3_corrections_summary.md)**
   — what S2 delivered, and what it deliberately left open.
5. **The six item-pool impact notes** — `phase1_`, `phase2_`, `phase3_`, `phase5_`, `phaseC3_`
   and `phaseC3_corrections_item_pool_impact.md`. D6.4 consolidates exactly these.
6. **[`DECISIONS.md`](../re-implementation-sep/DECISIONS.md)** — 77 decisions. D-003, D-024,
   D-026, D-031, D-043, D-052, D-053, D-061 all reach this phase. **Next free number: D-078.**

## What is actually here, measured

Measured at `444b8bf` on 2026-09-12. Re-derive in a worktree before acting.

### The inventory is stale — this is D6.1, not a preliminary to it

```
python -c "import csv,collections; r=list(csv.DictReader(open('docs/re-implementation-sep/template_inventory.csv',encoding='utf-8-sig'))); print(collections.Counter(x['instrumentation_class'] for x in r))"
```

| `instrumentation_class` | count |
|---|---:|
| A | 49 |
| B | 34 |
| C | 51 |
| **D** | **16** |

150 rows. **This is byte-identical to the pre-Phase-0 audit.** Six phases have landed and not
one row has been reclassified. So the exit-gate line *"Class D reduced from 16 to ≤4"* is not a
diff against a moving figure — it is a claim that has never once been tested, and the 16 is the
original 16.

**Establish how those columns were produced before you reproduce them.** The inventory carries
measured columns — `pct_step_values_recoverable`, `pct_answer_values_recoverable`,
`n_inline_computed`, `n_milestone_candidates`. D6.1 says *"full re-run of the original audit
harness"*. **I did not find that harness.** `tests/template_integrity/` holds T1–T8, the
instance dumps, the contract scan and the gate report; none of them writes this CSV. If no
harness exists, **D6.1 is a build task, not a re-run task, and its effort estimate is wrong** —
say so in a decision rather than hand-reproducing 150 rows of judgement and calling it a
re-audit. A reproduction whose method differs from the original cannot support "no template
regressed in class", because a class change and a method change are indistinguishable in it.

**Owner decision: best effort, time-boxed.** Look for the production method. If you find it, use
it. If you do not, **document that in two or three sentences, state the method you used instead,
and move on** — do not spend the phase on archaeology. One thing must carry forward from this
either way: if the method is new, **Reviewer F's gate has to be scoped to what is answerable** —
*"is this method sound, and does it classify the unchanged templates as the old one did?"* rather
than *"does it replicate the original?"*. Commissioning a replication of an unrecoverable method
produces a review that cannot file, which is what R6 exists to prevent.

**Searched exhaustively. The harness never existed — and the method survives anyway.**

The archaeology is done, so do not repeat it. Exactly **two commits** have ever touched
`template_inventory.csv`: `9105317` added it, and `e8302c6` (Phase 4 close-out) inserted two
columns. `git show --stat 9105317` contains **four documentation files and no `.py` at all**.
Pickaxe searches on the distinctive column names (`pct_step_values_recoverable`,
`blind_guess_floor`) return only prose. There is no deleted blob, no dangling object, no
untracked script, and the bare filenames in `.gitignore` were **never tracked in any branch**.

The cause is in the commissioning prompt: `01-template-structure-audit.md:146` says *"Stay on
`master`. Do not create branches or commits."* The audit ran read-only, in session, and its
measurement code was throwaway. **Nothing was lost; nothing was ever saved.**

**The method, however, is written down twice**, and the second is far better than the first:

1. `template_audit_report.md:507–514` — "Appendix — method": the AST sweep, the `sys.settrace`
   capture, the EXACT/ROUNDED/SCALED/MISSING matcher, with sample sizes.
2. **`phase0_baseline.md:289–348`** — the Phase 0 adversary **independently re-implemented** those
   measurements with self-contained runnable snippets, and **matched the audit exactly** on
   per-branch inline interpolations (chemical **65**, civil **37**, mechanical **89**) and to
   **−1.0%** on step tokens (17,002 against 17,177).

**So D6.1 is assembly, not invention, and it is cheaper than a rebuild rather than dearer.** Every
piece already exists: `core.py:53 discover()`, `core.py:181 generate(capture=True)` with its
tracer, `core.py:106 _walk_values`, `core.py:278 match_value()`, and
`checks/t5_binding.py:43 _is_bound()`.

**Calibrate against Phase 0 before publishing any number.** If the harness does not reproduce
65/37/89 and ~17,002, the method is not continuous with the original, and divergence from the
2026-09-05 CSV cannot be attributed to six phases of template change rather than to method drift.
Calibration is the gate on the re-audit itself.

**Split the CSV's provenance per column**, and say so in D6.3:

| columns | treatment |
|---|---|
| `n_inline_computed`, `pct_*_recoverable`, `n_steps`, `has_instance_branching` | **regenerate** — mechanically measurable |
| `n_milestone_candidates` | regenerate, but it is a **declared proxy** (report `:497`), not a count |
| `difficulty`, `est_effort`, `notes` | **carry forward**, marked inherited-not-remeasured |
| `answer_type`, `unit_system` | hybrid: hand-verified for ~90 rows, automated for the rest (`:499`) |
| `instrumentation_class` | a **re-classification under a restated rule** (D > C > B > A, report `:30`) — *not* a re-measurement, because which rows the five human branch audits overrode **was never recorded** |

**A defect in the existing CSV, found in passing and needing a decision.** Only **5 rows** carry
`blind_guess_floor` / `surface_model_heldout`, and **3 of the 5 have no recorded derivation
anywhere**. `template_levenspiel_plot_interpretation` carries **1.0000 / 1.0000** — a blind-guess
floor of 1.0 for an `array`-answer template is not a coherent statistic. Reviewer B's report
(`reviews/phase4_reviewer_b_pedagogy.md:37–61`) only ever derived the two `discrete_time_signals`
rows. Decide whether the column is repaired, emptied, or kept with its provenance stated.

### Corpus baseline at `444b8bf`

```
python -m tests.template_integrity.run --checks all          # ~2-3 min
python -m tests.template_integrity.phase5_contract_scan
python -m tests.comparators.cross_pair
python -m tests.trace_schema.audit_3_8
```

| check | result | note |
|---|---|---|
| T1 | 29 | pre-existing, stable across six phases |
| T2 | 0 | |
| T4 | 0 | |
| T5 | 66 | |
| T6 | **142** | **stale baseline, corpus-wide (D-043)** — not a Phase 6 regression |
| T7 | 83 | |
| T8 | 0 | |
| contract scan | 150/150 clean | 0 generation errors on 60,000 instances |
| `cross_pair` | PASS | |
| `audit_3_8` | 80/80 | |

**T6's 142 is the number this phase has to make a decision about.** It has been carried as "not
mine" by Phase 3, Phase 5 and Track B, each correctly. Phase 6 is the consolidation phase and
there is no later phase to carry it to. The decision is *whether the baseline is regenerated
under owner authority, and against what*, not whether you may quietly regenerate it. You may
not (D-043).

### The doubled-sign census (D6.7), re-measured here

`phase5_contract_scan` prints this on every run. **14 templates**, confirming D-061, and the
list crosses four branches — which is the argument for the shared helper:

```
template_continuous_to_discrete_conversion   template_coulombs_law
template_impulse_response_from_lccde          template_lorentz_force
template_mean_variance                        template_multi_segment_rod
template_nyquist_rate_determination           template_pitzer_correlation_z
template_sensible_heat_constant_cp            template_signal_operations
template_system_property_linearity            template_vorticity_check
template_wave_equation_interpretation         template_work_isothermal_virial
```

`_signed_term` and `_rect_str` are at
[`waves_and_phasors.py:20`](../../data/templates/branches/electrical_engineering/electromagnetics_and_waves/waves_and_phasors.py#L20)
and [`:30`](../../data/templates/branches/electrical_engineering/electromagnetics_and_waves/waves_and_phasors.py#L30),
used at 11 call sites in that one file. **Promote them before fixing anything** (Reviewer A,
|S|5.6): nine hand-written sign fixes is nine chances to write `+ {value}` again. The helper
lands in a shared module with its own planted-defect test — **two plants of materially
different surface form, written from the class definition rather than from the detector**
(SPEC-CHANGE 17).

Also printed, **not gated**: `degenerate_product_derivation` on 3 templates
(`impulse_response_from_lccde`, `standing_wave_formation`,
`undamped_response_initial_conditions`). Decide whether it is in D6.7's scope or a residual;
do not leave it printed-and-unread — *that* is Phase 5's named failure (`cross_pair` printed
`false rejects 56` on every run and nobody put it in a claim).

### Two spec numbers that do not survive contact

**D6.8's "199 characters" — re-measured, and the spec is right.** The figure does not appear in
`phase5_summary.md`, the document the spec cites, so it was flagged here as unsourced. Measured
over 2,250 instances it is **correct**: the corpus's longest answer span is
`template_autocorrelation_rect_pulse` at **199 characters** from the end of the `**Answer:**`
marker (block length 210, marker 11; 198 if the leading newline is also stripped). Median **64**,
p95 **165**, and **zero templates lack an answer marker**.

So the number was always right and only its provenance was missing — which is worth stating
plainly rather than quietly dropping the objection. **Cite this run when you set the bound**, and
**re-measure after D6.7 and D6.11**, both of which change emitted text.

**D6.11's "113 templates carry a unit" — re-measured, and the count moves.**
`phase5_summary.md:126` gives **four** counts, not one — *seed 0 115 / any 116 / always 113 /
invariant 103* — and the spec quotes 113 without its predicate.

Re-measured here over the 2,250-item pool, reusing T6's own `_UNIT_AFTER_NUM` and
`answer_block` so the definition is the repo's rather than a rival one, the four counts are
**seed 0 121 / any 121 / always 121 / invariant 97**, with **24 templates whose unit changes
with the seed** — not the 10 the earlier figures imply.

**The census does not reproduce across seed sets, and that is the finding.** Neither set of
numbers is wrong; they sample differently. So a fixed count quoted without its sample is not a
size estimate at all. **Size D6.11 off the actual generated pool**, publish the seed set with the
count, and treat the 24 seed-varying templates as the hard core: a per-item declaration for a
template whose unit is seed-dependent cannot be written once at template level, and one unit
applied to all `n` parts of a `multipart` answer is wrong for at least `n−1` — which is what made
13 templates reject their own gold.

### Units — D6.11 overlaps C1.4, deliberately

C1's `@units` field already declares units **per constants table**, and C3 populated it corpus-wide.
D6.11 declares units **per emitted answer part**. These are different objects and both are
needed, but **import C1's vocabulary rather than inventing a second one** — the spec's own
cross-reference at §"Phase 6's D6.11 should import rather than re-invent" says so, and Prompt 07
carried it as a coordination item. A `multipart` answer needs one unit **per part**; one unit
applied to all `n` parts is wrong for at least `n−1` and made 13 templates reject their own gold.

### The item pool — what D6.4 consolidates

Six measured events, every one from two worktrees in two separate processes:

| phase | templates moved | measured scope |
|---|---:|---|
| 1 | 10 (8 in the published pool) | 120 of 1,350 items = 8.9%; ~74 re-score, ~44 re-inference |
| 2 | 4 | `levenspiel_plot_interpretation` 300/300; two templates 300/300 answers |
| 3 | 2 | 2.30% and 0.65% of instances; answer *space* unchanged |
| 5 | 4 of 11 | emitted text only — **no question content and no answer value moved**, 22,000 instances |
| C3 | 1 + 1 | 1 answer-body (re-score *and* re-inference), 1 question-only (re-inference only) |
| C3 corrections | **13** | q 824, ans 713, sol 765; 0 generation errors |

**Re-scoring and re-inference are not the same cost and this table is the only place the
difference is recorded.** A question that changed with its answer unchanged needs re-inference
only. A question byte-identical with a *moved gold answer* is a **hidden constant** — it needs
re-scoring and nothing about the item looks different. The C3 corrections tranche contains
both directions. D6.4 must state them separately or it understates the second, which is the
dangerous one.

**Do not assume these sum.** Templates appear in more than one tranche, and the union is not
the total. Compute the union from the named templates, and publish the list.

### D-003 is closed by decision — and what replaces it is D6.12

D-003 has been carried open since 2026-09-05 across Prompts 05, 06 and 07: *do the raw
`inference_results/` generations still exist?* It gates any corrected results table.

Measured at `444b8bf`:

- `inference_results/` — **absent**. `evaluation/run_inference.py:23` writes to it; nothing does.
- `evaluation_results/` — **absent**.
- **`testset/` — absent.** `evaluation/run_inference.py:22` reads `../testset`. The input corpus
  the published numbers were generated from is not in this working copy either.
- The only `*results*.json` anywhere is `templates_annotation/annotation_app/src/llm_results.json`,
  which belongs to the annotation pilot.

**The repo owner has decided: the testset is rebuilt, inference is re-run, and evaluation is
re-run.** So D-003's question — can we recover the generations for a free re-score? — is
**moot by decision, not by evidence**. Close it as **D-078** recording both: the artefacts are
absent *and* recovery is not attempted, because everything downstream is being regenerated.

**D6.4 changes character accordingly.** The re-score / re-inference split stops being a budget
question and becomes an *explanatory* one: it is the record of what the six phases did to the
corpus, and it is how the paper explains why published numbers changed. Keep the distinction —
especially the **hidden-constant** direction, where the question is byte-identical and only the
gold answer moved — but state it as history, not as a recovery plan.

**D6.12 — regenerate the testset (new; SPEC-CHANGE 24).** The spec has no deliverable for this
because it assumed the pool already existed. It does not, and the machinery to rebuild it is in
worse shape than "absent". Measured at `444b8bf`:

**1. The generator covers three branches of five.** 20 template modules carry a
`main()` writing `testset/<branch>/<domain>/<area>.jsonl`:

| branch | modules with a generator | template functions |
|---|---:|---:|
| chemical | 7 | 30 |
| electrical | 7 | 30 |
| mechanical | 6 | 30 |
| **civil** | **0** of 11 | 30 |
| **industrial** | **0** of 11 | 30 |

Those 20 modules hold **exactly 90 template functions** — which is precisely the published pool's
*"1,350 items = 90 templates × 15 seeds"*. **Civil and industrial were never in the published
testset.** They were authored later and have no generator at all. So "regenerate the testset"
over 150 templates is a **67% scope expansion into two branches the published results have never
covered**, not a like-for-like rebuild. **The owner has decided to include them**: the new pool is
150 templates × 15 instances = **2,250 items**, against a published 1,350. The generator for
civil and industrial therefore has to be built, not merely re-run.

**2. `regenerate_testset.py` is gitignored and absent.** `.gitignore` line 39 lists it, along
with `template_loader.py`, `verify_fixes.py`, `locate_fixed_templates.py` and
`verify_manual_map.py`. None is on disk. Whatever orchestrated the published build is not
recoverable from this repository. **`testset` itself is gitignored (line 30)**, so the benchmark's
actual items have never been in version control — worth a sentence in the paper's reproducibility
statement regardless of what this phase builds.

**3. The generator is not reproducible, in all 20 modules.** Every one draws its per-item seed
from an **unseeded** global RNG and then shuffles:

```python
for _ in range(50):
    seed = random.randint(1_000_000_000, 4_000_000_000)   # outer RNG never seeded
    random.seed(seed)
```

Audited across all 20: `unseeded_draw` present in 20, outer `random.seed(<fixed>)` in **0**. An
individual item is reproducible *after the fact* because the seed is stored in its record, but
**the set is not** — every run emits a different testset. Phase 2's deliverable was *"seed every
generator"*; it seeded the templates, and nothing ever seeded the thing that calls them.
**Fix this before regenerating**, or the new testset is exactly as unreproducible as the one it
replaces, and no later phase can diff against it.

**4. The instance counts are incoherent.** `for _ in range(N)` is **50** in sixteen modules,
**200** in `mole_balances` and `harmonically_excited_vibrations`, and **3** in `electrostatics`
and `magnetostatics`. Total emitted: **5,615 records**, against a published pool of 1,350 — so
**the committed generators are not what produced the published testset**. The per-template
imbalance is 3 to 200, giving `mole_balances` 1,000 records and `magnetostatics` 3. Any aggregate
over that set is dominated by two modules. **Settled: 15 instances per template, uniformly, all
150 templates.** That matches the published pool's per-template depth and removes the 3-to-200
imbalance entirely.

**5. Build it on `discover()`, not on 22 new `main()` blocks.**
[`tests/template_integrity/core.py:53`](../../tests/template_integrity/core.py#L53) already
enumerates all 150 templates across all five branches deterministically; T1–T8 and both instance
dumps use it. Writing 22 more hand-rolled `main()` blocks would reproduce defects 3 and 4 in the
two branches that do not yet have them.

**Ordering — this is the part that costs money if it is wrong.** D6.7 changes emitted text on 14
templates and D6.11 changes per-part unit declarations. **Regenerate the testset only after those
land.** A testset built before them is stale the moment they merge, and inference then runs twice
across 11 models — the exact waste sync point S2 was created to prevent. Sequence: corpus changes
→ audit → testset → inference → evaluation.

### D6.5 — "promote to CI" is a build task, and it collides with the binaries

- **There is no CI.** No `.github/workflows`, no `tox.ini`, no `Makefile`, no `.gitlab-ci.yml`.
  D6.5's verb is wrong; budget accordingly.
- **`sympy` is not in `requirements.txt`** (D-053, still open) and
  [`tests/comparators/kinds.py`](../../tests/comparators/kinds.py) imports it. **A clean clone
  cannot run the comparator suite** — which is precisely what CI is. D-053 stops being
  housekeeping the moment D6.5 starts. `scipy` and `matplotlib` are also absent from
  `requirements.txt`; audit the whole file against actual imports rather than adding `sympy` alone.
- **`constants_integrity` cannot run in CI as-is.** 13 reference binaries, **507.4 MB**, are
  gitignored by the repo owner's decision and recovered by `fetch_references.py`. **Do not
  commit them, do not set up Git LFS, do not loosen the rule** to make CI green. The honest
  design is a **tiered** suite: the checks that need no binary run on every push; the
  citation-resolution checks run where the artefacts can be fetched, and **skip loudly** —
  never silently — elsewhere. A green CI that skipped the provenance suite without saying so
  would be a worse artefact than no CI.

## Rules carried forward

- **A test may not carry its own answer key** (D-034). `test_citations_resolve.py` P4 enforces
  this and must keep doing so.
- **Do not reorder a dict's keys** (D-031). Regrouping `CP_PARAMS` changed which substance 274
  of 300 seeds drew. If D6.7's helper promotion touches a module's import order, measure it.
- **A gate needs more than one term** (SPEC-CHANGE 18). "Class D ≤ 4" is passed by a
  reclassification that moves the boundary rather than the templates. State the predicate for
  each class and hold it fixed across before and after.
- **Every detector gets two planted defects of materially different surface form**
  (SPEC-CHANGE 17), written from the class definition, not from the detector.
- **Expect more than one defect class per unit of work** (D-028).
- **A thing that has just been corrected is when a fresh error is most likely** (C2 §8).
- **Write every patch as a file, never through a shell heredoc or a double-quoted `bash -c`.**
  Phase 4 lost four fixes to heredoc escaping, Phase 5 four more, and Track B lost a commit to
  an apostrophe. Several files in this repo are CRLF.
- **Set `PYTHONIOENCODING=utf-8`** before printing template output — the console is cp1252 and
  crashes on `→`/`Σ`/`²`/`·`.
- **Never pipe a check through `tail` and read the exit code as the check's.** Track B reported
  "exit code 0" twice for a battery that never ran; the pipe swallowed the status.

## Verification

```
python -m tests.template_integrity.run --checks all          # T1-T8, ~2-3 min
python -m tests.template_integrity.phase5_contract_scan
python -m tests.comparators.derive_bindings
python -m tests.comparators.cross_pair
python -m tests.comparators.score
python -m tests.trace_schema.audit_3_8
python -m tests.constants_integrity.test_citations_resolve
python -m tests.constants_integrity.test_chemical_thermochemistry
python docs/references/fetch_references.py --verify
```

**`derive_bindings` and `cross_pair` are not decoration.** D6.7 changes emitted text on 14
templates. A sign fix that alters an answer's shape can unbind a template bound in Phase 5.
**Re-run both after every tranche** and treat a newly unbound template as a finding, not as
noise.

Every one of these must be green — or carry a written, pre-existing reason (T1 29, T5 66, T7 83,
T6 142) — **before** you claim a class migration. An audit run over a corpus whose regressions
you have not checked measures the wrong thing.

## Reviews — read R0–R6 before dispatching any

- **Reviewer F — audit replication.** *The gate.* Given only the original audit report and the
  changed repository, **independently re-runs the audit** and compares. Must confirm: the
  class-D count fell as claimed, no new class-D template was introduced elsewhere, and no
  metric regressed silently. Explicitly tasked with **collateral damage** — a Phase-1 rounding
  change that broke a Phase-5 template's contract, a Track-B constant that moved a Phase-3
  trace.

**Give F the method problem honestly.** If D6.1 turned out to be a build rather than a re-run,
F cannot replicate an original method that was never recorded. Then F's gate becomes *"is the
new method sound, and does it classify the unchanged templates the same way the old one did?"*
— which is answerable. Commissioning a replication of an unrecoverable method produces a review
that cannot file, and R6 exists to prevent exactly that.

**Scope the review to R6**: one mandatory gate task, a time box stated as a number, tooling
supplied as working code, everything settled fenced off, and a measurement-ownership register
with **no number commissioned twice**.

**Independence, concretely:**
- Code comments, commit messages and summaries are **claims under test**, not background (R1.1).
- Give the reviewer a **frozen SHA**, then **do not commit *or write to the working tree* until
  it files.** A reviewer reading at a frozen SHA sees uncommitted edits too. Phase 5 broke this.
- Tell the reviewer to read with `git show <sha>:<path>`, never `git show <sha> --stat`.
- **Commit the report unmodified and before any fix.**

**Every §5 suggestion is triaged (R4) before the phase closes**; an untriaged suggestion blocks
the gate exactly as a CONFIRMED finding does.

## Git

- Branch **`redesign/phase6-consolidation`** off `master`. Merge `--no-ff` when the gate passes.
- **Track B's work is already on `master` and was committed there directly** — 21 commits
  between `4eb6488` and `444b8bf`, against this repo's convention of a `redesign/*` branch and
  a merge commit. It is recorded here so you do not mistake it for a convention change. Do not
  repeat it.
- `master` is **185 commits ahead of `origin`**. Nothing has been pushed. **Do not push.**
- **Large source binaries are NOT committed** — 13 files, 507.4 MB, the repo owner's decision.
  **Do not force-add one, do not set up Git LFS, do not loosen the rule**, including for CI.
- **Never commit anything from `pilot/references/`.** It is gitignored for a reason.
- Logical commits, not one lump.
- End commit messages with:
  `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`

## Exit gate

- [ ] `template_inventory.csv` regenerated, with the **method stated** and the predicate for
      each class published
- [ ] Class D reduced from 16 to ≤ 4 (the Phase-4 four, **reclassified not fixed**), or the
      shortfall named per template with a reason
- [ ] No template regressed in class; **measured**, not inferred from summaries
- [ ] D6.2 class-migration table; D6.3 `template_audit_report.md` updated with post-change figures
- [ ] D6.4 consolidated item-pool statement — **re-score and re-inference separated**, the
      union of affected templates listed, and D-003's closure stated
- [ ] D6.5 CI standing up the regression suite, **tiered**, with binary-dependent checks
      skipping **loudly**; `requirements.txt` audited against actual imports (D-053)
- [ ] D6.6 residual-risk register: everything knowingly not fixed, and why — **including Track
      B's 177 `[UNVERIFIED]` rows and its §3 substance defects**
- [ ] D6.7 doubled sign cleared on 14 templates, **behind a promoted shared helper with
      planted-defect tests**
- [ ] D6.8 answer-span shape assertion, its bound **re-measured** rather than inherited
- [ ] D6.9 and D6.10 **attempted**; fixed with measured before/after, or recorded as infeasible
      with the measurement proving it. **Zero false accepts preserved either way** — a comparator
      that decides more by deciding wrongly is worse than one that abstains
- [ ] The seven Track B rows deleted, R-410A from **both** its tables, with the P6 measurement
- [ ] D6.11 per-item and per-part units, importing C1.4's vocabulary, predicate published
- [ ] **D6.12 testset generator rebuilt on `discover()`** — seeded reproducibly end to end, one
      stated instances-per-template figure applied uniformly, all five branches; the same seed
      reproduces the same set, **demonstrated by two runs diffed**
- [ ] **D6.12 regenerated after D6.7 and D6.11 merged**, not before; item count and composition
      stated against the published 1,350
- [ ] Pool is **2,250 items = 150 templates × 15 instances**, uniform, all five branches
- [ ] T1–T8, contract scan, `derive_bindings`, `cross_pair`, `score`, `audit_3_8` — no
      regression, **measured**; **T6 baseline not regenerated**
- [ ] T6's 142 given a decision with an owner, not carried forward again
- [ ] Reviewer F filed; every finding and §5 suggestion triaged; every `SPEC-CHANGE` actioned
- [ ] `phase6_summary.md` in the shape of `phase5_summary.md`, including **your own errors**. A
      summary recording none is not a clean phase; it is an unexamined one.

## Constraints and cautions

- **Do not modify a check or an oracle to make something pass.** A check that looks wrong is a
  finding and a `SPEC-CHANGE`.
- **Do not re-open a merged phase's fixes.** If one is wrong, that is a finding for the register
  and a decision, not a silent revert.
- **Do not replace a substance, material or fluid on your own authority.** The seven deletions in
  §"Track B — the seven deletions" are owner-directed and are the *only* substance changes
  authorised. A **replacement** — substituting a different substance for one of them — is not
  authorised, and was refused during C3 for good reasons: R-22 for R-410A would change the row's
  values by +9.78% / −33.21% *and* duplicate a fluid already in the table.
- **P6 is a real constraint.** D6.7 changes emitted text on 14 templates: that is a P6 event and
  needs a two-worktree, two-process measurement like every one before it.
- **Where you are uncertain, say so and name the evidence that would settle it.**

## Open items inherited — carry forward, fix only if they block your gate

- **T6's corpus-wide stale baseline** (D-043) — 142/150. Yours to decide, per the exit gate.
- **D-053** — `sympy` (and `scipy`, `matplotlib`) missing from `requirements.txt`. **Blocks
  D6.5**, so it stops being optional here.
- **Track B residuals**, all registered in
  [`phaseC3_residual_register.md`](../re-implementation-sep/phaseC3_residual_register.md):
  177 `[UNVERIFIED]` rows (mechanical 148, chemical 13, civil 8, electrical 4, industrial 4);
  11 rows blocked on paywalled standards (AISC, ASTM A992/A36, ACI 318-19, ASCE 7-22, IEC 60063,
  MIL-A-8625, US Standard Atmosphere 1976); ~19 needing a mechanical-properties handbook;
  ~86 bulk materials and mixtures no chemical database indexes.
- **7 `[DERIVED]` constants still UNEXECUTED** — 4 Shomate refits, `Air(g)`, 2 by-definition
  ceilings. The resolver recomputes 8 of 15 and says so on every run.
- **4 `tol=basis=condition` tolerances the resolver cannot size** — especially `Glycerine`,
  whose tolerance absorbs a *composition* ambiguity. Named in `phaseC3_summary.md` as the
  weakest warrant in the branch.
- **The stale annotation fork** — `templates_annotation/annotation_app/.../constants.py` holds
  pre-C3 chemical values (R-12 `v_g 0.0268`, R-22 `v_f 0.000845`). **Deliberately not synced**,
  because that directory backs the annotation pilot and a frozen snapshot may be intentional so
  annotations stay reproducible against what annotators saw. **Owner decision: leave it. Do not
  sync it.** It is out of scope for this phase. Name it in one line in D6.6 — because D6.4
  otherwise claims a corpus-wide constant state that this directory contradicts — and do nothing
  else with it.
- **G F-5** (`@domain: none` vocabulary) and **G F-7** (NAVFAC manuals cited by zero tags).
- **`CP_PARAMS` origin vs verification** (Reviewer G, G-7) — unfixable without a citable
  Smith–Van Ness copy; recorded, not closed.
- **The value-extractor spike (D-002)** — *"the gate for the whole project"*, and still unwritten.
  It is **not** Phase 6's, but Phase 6 is the last Track A phase, so after this the sequence has
  nothing queued. Say so in the summary.
