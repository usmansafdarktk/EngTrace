# D3.3 — The `iteration` and `decision` trace node types

**Schema version:** `1.5` (§7 amended by Phase 4, D-054) · **Status:** normative for the milestone model ·
**Phase:** 3 · **Revised:** 2026-09-06
**Companion:** [`template_redesign_spec.md`](../template_redesign_spec.md) §3.3 ·
[`phase3_summary.md`](phase3_summary.md) · conformance corpus in
[`phase3_conformance/traces.json`](phase3_conformance/traces.json)
**Reference implementations:** the `trace_nodes` local of
`template_normal_depth_iteration`
(`data/templates/branches/civil_engineering/water_resources/uniform_flow.py`)
and the return value of `_t19_assign`, bound to `trace_nodes` in
`template_line_balancing_heuristic`
(`.../industrial_engineering/production_and_inventory/production_planning.py`).
**Extractor:** `python -m tests.trace_schema.extract <out.json> <n_seeds>`.

> **Revision note — 1.4 → 1.5.** §4.6 — the frame-arithmetic mechanism 1.4 added
> to close variant 5 — was reviewed and **is not sound as written**. The
> mechanism works; its *applicability* was optional. §4.3.9 quantified over the
> relations **present**, not over the frame's **symbols**, so deleting the one
> relation defining the residual left `A`, `P` and `AR` vacuously checked and
> `g` unconstrained. **Variant 5 returns unchanged, by removing one line** — the
> reviewer drove a trace to 3.492 m and to exactly 2.000 m against a gold
> 1.380 m, both passing clean. 1.5 adds the coverage requirement: every
> `evaluation_symbols` member except the iterate must have exactly one relation.
>
> Two further defects of mine in the same section. `round(AR, 3)` **restated a
> precision `symbol_precision` already declares**, and nothing required the two
> to agree — `round(AR, 1)` yields a conforming 1.383 m. That is a §3.8
> restatement violation *inside the fix that enforces §3.8*; the form is now
> `round(sym)`, taking its digit from `symbol_precision`. And the **6.7%**
> bind-to-stored miss rate §4.6 cited does not reproduce against the corpus the
> document ships: it is the preamble-only rate over 2,000 seeds (re-measured:
> 6.6%), while the shipped 40-seed corpus gives **24.9%**, because mid-iteration
> frames carry larger `A` and `P` and so more absolute rounding error. Both
> figures are now stated with their populations.
>
> **§3.8's amended rule, applied honestly, condemns its own table:** five of nine
> rows are declared once and checked against nothing. 1.5 gives it the
> distinction it lacked — **problem data** versus **derived data** — and states
> the severity gap that makes the distinction load-bearing: a free constant makes
> a trace a correct solve of a *different* problem, which §7 catches; a missing
> relation makes it a correct solve of *no* problem, and §7.2 grants it full
> process credit.

> **Revision note — 1.3 → 1.4.** §3.8 held: the reviewer could not break the
> family it names, and with `universe`, `precedences`, `item_measures` and the
> budget declared, exactly one `decision` trace conforms. But **§3.8's audit had
> never been run against the two node types shipping in the same document**, and
> running it — as the reviewer said, an implementer's job, not a fourth
> reviewer's — finds two values it did not account for.
>
> **`budget_total` was restated in every element.** §3.8's own audit question
> convicts it verbatim: *"where is the single place it is declared? If the answer
> is 'in each element that mentions it', the check is not yet a check."* Declare
> the budget too small and the greedy rule *honestly* opens an extra element, so
> `n` — the answer — is wrong with every clause satisfied. `capacity_constant`
> did not help: restating the same wrong budget everywhere is constant. 1.4
> declares `termination.budget` once and deletes the flag.
>
> **An `iteration` node's answer was free.** §4.3.1 recomputes the update *from*
> the residuals, and the residuals came from evaluation frames whose arithmetic
> §10 had conceded was unchecked. That concession was load-bearing for an
> attacker: a constructed trace **reported 3.501 m where the gold answer is
> 1.380 m and passed with zero failures**, taking **full process credit** under
> §7.2 — a model emitting arithmetically consistent nonsense rewarded by the
> mechanism built to detect it. Note this involves no restatement at all, so
> §3.8 as worded in 1.3 did not catch it: "declared once" is necessary and not
> sufficient. §3.8 is amended to **declared once *and checked against something*,
> or recomputed from something that is**, and `frame_relations` closes the hole:
> the frame's arithmetic is now data (§4.6), verified to recompute exactly over
> 12,735 frames.
>
> §10's "frame-internal physics is not checkable" is therefore **withdrawn**. It
> was a limitation reported as a principle, and it cost the deliverable its
> central claim for three revisions.

> **Revision note — 1.2 → 1.3.** The fresh reviewer re-ran its exploit against
> 1.2 and **routed around the fix**. `filter_relation` closed the attack it had
> filed — lying about the admissibility *flag* is now rejected, at the right
> clause — but lying about the flag's *input* works just as well: inflate an
> item's `measure` in the trial where you want it excluded (so the flag is
> honest about a dishonest number), restore it where you want it committed, and
> every clause of §5.3, §5.4, §6 and §8A holds while `n` — the answer — comes out
> wrong. **Two revisions had now closed the same hole one variant at a time.**
>
> 1.3 therefore adopts the reviewer's own recommendation and states the invariant
> generally, as **§3.8**, rather than patching a third variant: *every value a
> check consumes is either declared once on the node, or recomputed from
> something that is.* `item_measures` is the instance of it that closes F3′ —
> item measures are declared once beside `precedences`, so no trial can restate
> one. §3.8 is the rule to audit any future node type against, and it is what
> should have been written in 1.1.
>
> Also in 1.3: `preamble_binding` names the frame role explicitly instead of
> relying on a name prefix; `accumulator_symbol` becomes `accumulator_role`, the
> last raw symbol in an otherwise all-roles §5; §3.7's cumulative reading is
> corrected to a **prefix** reading, which had been asserted fixed in 1.2 and was
> not — the section was byte-identical between the two versions; and the
> `schema_version` row of §3.1 said `"1.1"` while the corpus said `"1.2"`, which
> taken literally rejects every trace at step 1.

> **Revision note — 1.1 → 1.2.** Version 1.1 was reviewed twice: by the round-1
> reviewer, checking whether its fifteen findings were addressed rather than
> relocated (**14 addressed, 1 partially, 0 not addressed**), and by a **fresh**
> reviewer who had not seen round 1 and implemented a verifier against 1.1 from
> scratch. Both independently found the same defect, from opposite directions:
> **1.1 fixed `iteration` by inverting the problem onto `decision`.** `roles` was
> marked unconditionally required, yet no `decision` node carried it and §5 still
> defined the type by its reference template's literal symbol names — so the
> fresh verifier scored 40 pass / 40 fail, and reached 80/80 only by
> special-casing `decision` out of the very check `roles` exists to abolish.
> 1.2 gives `decision` the same treatment `iteration` got: role maps at all three
> nesting levels, and §5 rewritten over roles.
>
> The fresh review also found a **soundness** defect that neither the corpus nor
> round 1 exposed: §6 accepted a `decision` trace whose **answer was wrong**.
> `selection.filter` was carried as data and never checked, and §5.4.4 optimised
> only over the filtered set, so a trace could lie about what fits, force an extra
> element, and change `n` — which §2 and §7.3 make *the answer* — while
> satisfying every clause. 1.2 adds `selection.filter_relation`, making the
> admissibility flag itself checkable.
>
> **The rename test now passes for `iteration` and is the standard 1.2 holds
> both types to.** A reviewer built a synthetic `iteration` node with every
> symbol renamed, a *different* recurrence and nine elements, and an unmodified
> spec-faithful verifier accepted it — then failed it on four clauses when one
> value was perturbed. That is what "a type, not a description" means
> operationally, and it is the test to run on any node type added later.
>
> Full finding tables and dispositions in
> [`phase3_summary.md`](phase3_summary.md) §11.

---

## 1. What this is for, and what it replaces

Every other template in the corpus emits a **flat milestone list**: a fixed
sequence of steps, each with a stable `{id, symbol}` pair, known before the
template runs. Verifying such a trace is a per-symbol numeric and dimensional
comparison, and no LLM is needed in the critical path.

Two templates cannot be represented that way, because the *shape* of their trace
is a function of their data:

| Template | What varies | Measured over 4,000 seeds |
|---|---|---|
| `normal_depth_iteration` | interpolation updates inside Step 3 | `{1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53}` — **32.9% need more than three** |
| `line_balancing_heuristic` | the number of steps itself | 5 steps on 52.9%, 6 on 47.1% |

A milestone model with a fixed `{id, symbol}` per step cannot express either:
`y_curr`/`g_curr` are rebound on every pass of the first, and the second has no
fixed number of stations to name.

This document specifies the two node types that close that gap. **It is the
deliverable the milestone model is built on**, so it is written to be
implementable by someone who did not design it: §6 is a normative algorithm, §7
is the comparator, §8 lists the traces a conforming verifier must *reject*, and
§10 states what is deliberately **not** specified, so an implementer knows where
the edge is rather than discovering it.

**The governing rule of this revision:** *anything a verifier must check is
declared on the node.* Symbol names, the update recurrence, display precisions,
the rounding mode, the precedence relation and the selection rule are all data.
Prose fields exist, are suffixed `_prose` or `_notes`, and are **never**
normative.

---

## 2. The finding that shapes this specification

`template_redesign_spec.md` §3.2 says the two node types "are the same shape (an
ordered list of homogeneous sub-traces with stable within-element symbols and a
termination predicate)" and recommends specifying them together.

**The structural claim is right and the equivalence is wrong.** Specifying them
as one type produces a verifier that is wrong on one of the two templates,
because the two differ on a property that decides how they are scored:

> **Is the cardinality of the sequence an observable of the answer?**

- `normal_depth_iteration`: the answer is the converged depth `yn`. The update
  count does not enter it. Worse, the count is **not stable under permitted
  numerical slack** — the termination test is `|y_next − y_curr| < 0.002 m` on
  values rounded to 4 dp, so a solver who carries more digits can converge in a
  different number of updates and reach the same depth. Marking that solver
  wrong would be marking them wrong for being more careful.
- `line_balancing_heuristic`: the answer is
  `balance delay = (n·CT − Σt)/(n·CT)·100`, and `n` **is the element count**.
  The count is also *exactly* determined — every fit decision is an integer
  comparison, so no rounding slack can move it.

So the discriminating rule, stated so it can be applied to a template this phase
never saw:

> **Cardinality is `incidental` when it is sensitive to the numerical slack the
> comparator already tolerates, and `answer_bearing` when it is invariant under
> that slack and appears in the answer.** These are the only two cases; a
> quantity that is both sensitive to slack and present in the answer is an
> ill-posed item, not a node-type question.

Both node types therefore share one base structure (§3) and differ in exactly
three places: the **termination kind**, the **cardinality semantics**, and the
**element comparison rule**. Recorded as **D-038**.

*Independently corroborated:* Reviewer D reports that implementing the
result-derivation and cardinality-symbol clauses
"forced two distinct rules, so a merged node type would have had to pick one".

---

## 3. Common base — the `sequence` node

### 3.0 Container

The extractor yields **exactly one node object per instance**, bound to the
template local named `trace_nodes`. The name is plural for historical reasons;
the value is a single JSON object, never a list. A future template emitting more
than one node will bind a list, and **must** raise `schema_version` when it does.

### 3.1 Fields

| Field | Type | Req. | Meaning |
|---|---|:--:|---|
| `schema_version` | string | ✔ | `"1.5"`. A verifier must refuse a major version it does not know. |
| `node_type` | `"iteration"` \| `"decision"` | ✔ | Selects §4 or §5. |
| `node_id` | string | ✔ | Stable within a template; the anchor a comparator reports against. |
| `cardinality` | `"incidental"` \| `"answer_bearing"` | ✔ | §2. Fixed per `node_type` (§4, §5). Determines whether §7.2 or §7.3 applies. |
| `cardinality_symbol` | string | cond. | Required **iff** `cardinality == "answer_bearing"`. Absent otherwise. |
| `roles` | object | ✔ | Maps each **type-level element role** of §4.1/§5.1 to this template's own symbol name. This is what makes the type reusable; see §3.3. **Required on both types** — 1.1 required it and shipped `decision` nodes without it. |
| `termination` | object | ✔ | §3.4. |
| `rounding` | string | ✔ | The rounding mode every `round`-like step in §6 and §7 uses. `"decimal-half-up"` for both reference templates. A verifier must not assume binary `round()`. |
| `symbol_precision` | object | ✔ | `{symbol: decimal_places}` for every symbol whose display precision any check depends on — element symbols, evaluation-frame symbols and trial symbols alike. §7.1's tolerance is derived from this and **never** inferred from the data. |
| `element_symbols` | array of string | ✔ | The symbol set every element carries. **Stable across elements** — this is what replaces the flat model's per-step `{id, symbol}`. |
| `elements` | array of object | ✔ | Ordered, homogeneous. Each has exactly the keys in `element_symbols`. **Length ≥ 1** (see §6 step 4). |
| `carry` | object | ✔ | §3.5. **Machine-checked path expressions only.** `{}` when the elements are independent. |
| `carry_notes` | object | ✔ | §3.5. Prose relations that are *not* checkable. **Never normative.** `{}` when there are none. |
| `result` | object | ✔ | `{symbol, value, dp, from, unit?}`. `from` is a §3.6 derivation expression, so §6 step 8 is data rather than a per-type special case. |
| `preamble` | array of object | ✖ | Fixed frames the question prescribes *before* the sequence starts. |
| `preamble_symbols` | array of string | cond. | Required iff `preamble` is present. |
| `preamble_binding` | object | cond. | `iteration` only, and **required iff `preamble` is present**. Maps an element role to `{frame, frame_role}`. §4.4. *(1.3: 1.2 made it optional, so omitting it escaped the check entirely — Reviewer D2, F4.)* |
| `evaluation_symbols` | array of string | cond. | `iteration` only; the symbol set of a nested evaluation frame. |
| `evaluation_roles` | object | cond. | `iteration` only; §4. |
| `update_relation` | string | cond. | `iteration` only; §4.2. |
| `constants` | object | cond. | `iteration` only; the item's fixed quantities, declared once. §4.6. |
| `frame_relations` | array | cond. | `iteration` only; ordered `[symbol, expression]` pairs giving each frame symbol's arithmetic. §4.6. |
| `trial_symbols` | array of string | cond. | `decision` only; §5.1. |
| `trial_roles` | object | cond. | `decision` only; binds the trial roles of §5.1. |
| `candidate_roles` | object | cond. | `decision` only; binds the candidate roles of §5.1. |
| `eligible_symbols` | array of string | cond. | `decision` only; §5. |
| `selection` | object | cond. | `decision` only; §5.2. The selection rule as data. |
| `rule` | string | ✖ | `decision` only. Human-readable restatement of `selection`. **Advisory.** |

> **Removed in 1.1.** Version 1.0 required *"`node_id` must not encode the
> element count"*. That is normative and untestable — any heuristic
> false-positives on a legitimate id like `t04_*` at four elements (Reviewer D,
> F13). Dropped rather than left unenforceable.

### 3.2 Prose fields are never normative

Exactly three fields carry prose: `termination.predicate_prose`, `rule`, and the
values of `carry_notes`. **A verifier must not parse any of them**, must not
`eval` them, and must not let them affect a verdict. Everything they describe is
also present as data. Version 1.0 made `termination.predicate` both required and
normative while giving it no grammar and no binding environment, so any verifier
that evaluated it raised `NameError` on all 80 gold traces (Reviewer D, F5).

### 3.3 `roles` — how one type serves many templates

A node's `element_symbols` are the template's own names (`y_curr`, `g_prev`,
`station_id`). The **roles** in §4 and §5 are type-level. `roles` binds them:

```json
{"index": "k", "iterate_prev": "y_prev", "iterate_curr": "y_curr",
 "iterate_next": "y_next", "residual_prev": "g_prev",
 "residual_curr": "g_curr", "change": "change",
 "converged": "converged", "evaluation": "evaluation"}
```

`decision` binds three levels, because its structure nests:

```json
"roles":           {"index": "station_id", "budget_total": "capacity",
                    "budget_initial": "remaining_initial", "trials": "trials",
                    "committed": "assigned", "budget_final": "remaining_final"},
"trial_roles":     {"candidates": "eligible", "chosen": "chosen",
                    "budget_after": "remaining_after", "closes": "closes"},
"candidate_roles": {"item": "task", "measure": "duration",
                    "admissible": "fits"}
```

A verifier resolves a role to a value as `element[roles[role]]`. Every role its
`node_type` declares mandatory must be present as a key of `roles`, and every
value of `roles` must appear in `element_symbols`. **A verifier that reaches for
a literal symbol name is not conforming** — that was exactly the 1.0 defect, and
it is why `linear_reservoir_routing_step`, whose iterate is not called `y`, would
have failed a 1.0 verifier at shape-check for no reason (Reviewer D, F3).

### 3.4 `termination`

| Field | Type | Req. | Meaning |
|---|---|:--:|---|
| `kind` | `"convergence"` \| `"exhaustion"` | ✔ | Which of §4/§5 applies. |
| `satisfied` | boolean | ✔ | Whether the emitted sequence terminated by the predicate. **A gold node with `satisfied: false` is malformed** (§8A.1). |
| `max_elements` | integer | ✔ | The budget; `len(elements) <= max_elements` is an invariant, not a hope. |
| `quantity_role` | string | cond. | `convergence` only: the role the tolerance test reads. |
| `comparison` | `"lt"` \| `"le"` | cond. | `convergence` only: the sense of the test. |
| `tolerance` | number | cond. | `convergence` only. |
| `scope` | `"cumulative"` \| `"per_element"` | cond. | `exhaustion` only. §3.7; the cumulative reading is a **prefix**. |
| `accumulator_role` | string | cond. | `exhaustion` only: the element **role** that accumulates. *(1.3: was `accumulator_symbol`, the last raw symbol in an otherwise all-roles §5 — Reviewer D2, F11.)* |
| `item_measures` | object | cond. | `exhaustion` only: `{item: measure}` over `universe`, declared **once**. §3.8, §5.4.3d. |
| `budget` | number | cond. | `exhaustion` only: the per-element budget, declared **once**. §3.8, §5.4.1. *(1.4: was restated per element as `budget_total`, which made `n` free — Reviewer D2, variant 4.)* |
| `universe` | array | cond. | `exhaustion` only: the finite set being exhausted. |
| `precedences` | object | cond. | `exhaustion` only: `{item: [prerequisite, …]}` over `universe`. §5.4 invariant 3b. |
| `predicate_prose` | string | ✖ | Human-readable. **Advisory** (§3.2). |

**On `satisfied`.** In a gold trace it is always `true`, because the template
raises rather than emitting an unconverged trace — deliberately a `raise` and not
an `assert`, since `python -O` strips asserts and the prose tells the reader the
tolerance was met (Phase 2 Reviewer C, C-6; Phase 3 D-041). The field exists
because a *model's* trace may legitimately carry `false`, and the comparator must
be able to say so rather than crash.

### 3.5 `carry` and `carry_notes`

`carry` is what makes a variable-length sequence checkable without replaying it.
It maps each symbol of element *k+1* to a **path expression** over element *k*:

```json
{"y_prev": "y_curr", "g_prev": "g_curr",
 "y_curr": "y_next", "g_curr": "evaluation.g"}
```

reads: *element k+1's `y_prev` is element k's `y_curr`*, and so on.

**Path grammar (normative).** A path is one or more identifiers joined by `.`,
each identifier matching `[A-Za-z_][A-Za-z0-9_]*`. The first identifier must be a
member of `element_symbols`; each subsequent one indexes into the object the
previous selected. Nothing else is a path — no calls, no indices, no spaces.

Everything that is *not* a checkable projection goes in **`carry_notes`**, whose
values are prose and are never parsed. Version 1.0 put both in one map with no
syntax to tell them apart, so a verifier could not distinguish `"capacity"` (a
valid path) from a description, and guessing wrong is precisely the silent
failure §8B.10 exists to prevent (Reviewer D, F6). The split removes the question.

### 3.6 `result.from` — derivation expressions

`from` names how `result.value` is obtained, as one of a closed set:

| Form | Meaning |
|---|---|
| `"last.<role>"` | The named role of the final element, rounded to `result.dp` under `rounding`. |
| `"len(elements)"` | The element count. Valid only when `cardinality == "answer_bearing"`. |

A verifier implements §6 step 8 by dispatching on this field, not on `node_type`.

### 3.7 `termination.scope` — which reading of the accumulator

`exhaustion` nodes accumulate across elements. With `scope: "cumulative"` the
predicate reads *the union of `accumulator_role` over elements 1..k*; with
`"per_element"` it reads element-locally.

**The reading is a prefix, not a total.** §6 step 7 requires the predicate to
hold on the last element **and on no earlier one**, which is only decidable if
the predicate can be evaluated at each prefix. A reading over *all* elements is
constant in `k` and makes that clause unsatisfiable. 1.1 wrote "union over all
elements"; 1.2 claimed to have fixed it and left the section byte-identical; 1.3
actually fixes it (Reviewer D2, F7 — twice reported, and the second report is why
it is stated this precisely).

The distinction is not cosmetic in the other direction either: read per-element,
the reference `decision` node is **unsatisfiable**, because the largest single
station holds 3 items against a 5-item universe.

### 3.8 Declared once, or recomputed — the rule the other clauses are instances of

> **Every value a check consumes is either declared once on the node *and
> checked against something*, or recomputed from something that is. A value a
> trace may restate freely is not evidence; a value declared once but checked
> against nothing is not evidence either; and a check that reads one is not a
> check.**

This is the general form of a hole that was closed twice, one variant at a time,
before being stated. In 1.1 the admissibility flag was carried per-candidate and
believed; 1.2 made it recomputable from `filter_relation` — and the same attack
then worked one level down, by restating the *measure* the relation reads. A
candidate a trace intends to exclude had a free measure in the trial that
excluded it, so `n` could be inflated with every invariant satisfied.

**Two kinds of value, and the distinction is load-bearing.** Applying the rule
honestly to this document's own table shows five of nine rows are declared once
and checked against *nothing* — so the rule as first written condemns its own
table. The resolution is not to weaken it but to separate two cases:

- **Problem data** — the question's own givens (`constants`, `universe`,
  `precedences`, `item_measures`, `budget`, `tolerance`, `max_elements`).
  Declaring these once **is** the whole requirement. Nothing on the node can
  check them, because nothing on the node binds the node to its question; that
  is §7's job and D-039's open item.
- **Derived data** — anything the trace computes (`admissible`, the update, the
  frame symbols, the result). These **must** be recomputed from problem data.
  Declaring one once and checking it against nothing is exactly as free as
  restating it.

**The severity gap is why the distinction matters, not bookkeeping.** A free
*problem* value makes a trace a correct solve of a **different** problem, which
§7 catches when it compares against gold. A free *derived* value makes it a
correct solve of **no** problem — and §7.2 grants it full process credit, because
every element is internally consistent. That is variant 5, and it is the more
dangerous of the two.

| Value | Kind | Declared once at | Recomputed by |
|---|---|---|---|
| the item's constants | problem | `constants` | — (§7) |
| the task universe | problem | `termination.universe` | — (§7) |
| precedences | problem | `termination.precedences` | — (§7) |
| item measures | problem | `termination.item_measures` | — (§7) |
| the budget | problem | `termination.budget` | — (§7) |
| convergence tolerance | problem | `termination.tolerance` | — (§7) |
| admissibility | **derived** | — | `selection.filter_relation`, §5.4.3c |
| eligibility | **derived** | — | §5.4.3b, against `precedences` |
| the update | **derived** | `update_relation` | §4.3.1 |
| frame arithmetic | **derived** | `constants` + `frame_relations` | §4.3.9, §4.6 |
| the result | **derived** | `result.from` | §6 step 8 |
| display precisions | problem | `symbol_precision` | — |
| rounding mode | problem | `rounding` | — |
| symbol names | problem | the role maps | §6 step 2 |

**Audit any new node type against this table before adding it**, and against all
three questions:

1. *For every value a check reads, where is the single place it is declared?* If
   the answer is "in each element that mentions it", the check is not yet a
   check. **`budget_total` failed this and shipped in 1.2 and 1.3.**
2. *Is it problem data or derived data?* If derived, what recomputes it? A
   derived value declared once and verified against nothing is free.
   **Evaluation frames failed this and shipped in 1.0 through 1.3.**
3. *Does the recomputation cover every value, or only the ones present?* A check
   quantified over the relations a trace supplies is one a trace can opt out of
   by supplying fewer. **`frame_relations` failed this and shipped in 1.4.**

Question 3 is the newest and the one that has caught the most: it is the general
form of "the mechanism is sound, its applicability is optional".

**This audit is the implementer's job and it has been run late every time.** Each
of the three failures above was found by running it *after* shipping the version
that introduced them — twice by a reviewer, once by me. Run it before.

---

## 4. `iteration`

An ordered list of homogeneous refinement steps, terminated by a numeric
predicate on a converging quantity.

```
cardinality      : "incidental"      (fixed for this type)
termination.kind : "convergence"     (fixed for this type)
```

### 4.1 Element table (normative)

Roles, not symbol names. Bind them through `roles` (§3.3).

| Role | Type | Req. | Meaning |
|---|---|:--:|---|
| `index` | integer | ✔ | 1-based, contiguous, strictly ascending across `elements`. |
| `iterate_prev` | number | ✔ | The older of the two iterates the update reads. |
| `iterate_curr` | number | ✔ | The newer of the two. |
| `residual_prev` | number | ✔ | The residual at `iterate_prev`. |
| `residual_curr` | number | ✔ | The residual at `iterate_curr`. |
| `iterate_next` | number | ✔ | The update's output. |
| `change` | number | ✔ | `\|iterate_next − iterate_curr\|`. |
| `converged` | boolean | ✔ | The termination test's value on this element. |
| `evaluation` | object \| `null` | ✔ | The evaluation frame at `iterate_next`; `null` on the terminating element. |

An **evaluation frame** has exactly the keys in `evaluation_symbols`, and
`evaluation_roles` binds two roles into it: `iterate` (the depth/abscissa the
frame was evaluated at) and `residual` (the value the update consumes).

*One element is one update plus the evaluation its output feeds*, so an element
is self-contained and `carry` is a pure projection.

### 4.2 `update_relation`

The recurrence, as an arithmetic expression **over role names**:

```
iterate_curr - residual_curr * (iterate_curr - iterate_prev)
             / (residual_curr - residual_prev)
```

**Grammar (normative).** Role identifiers, decimal numeric literals, the binary
operators `+ - * /`, unary `-`, and parentheses. Nothing else — no calls, no
names that are not roles of this node, no exponentiation. A verifier evaluates it
by substituting the element's role values. This is what replaces 1.0's *"for the
reference template that relation is the secant step"*, which left the formula out
of the node entirely and forced a verifier to hardcode it (Reviewer D, F3).

### 4.3 Local invariants (normative)

For every element:

1. `update_relation`, evaluated on the element's roles, equals `iterate_next` to
   the §7.1 tolerance for `iterate_next`.
2. `change == |iterate_next − iterate_curr|`, to the §7.1 tolerance for `change`.
3. `converged == (change < tolerance)` for `comparison: "lt"`, `<=` for `"le"`.
4. `converged` is `true` on the **last** element and `false` on every other.
5. `evaluation` is `null` iff `converged` is `true`.
6. When `evaluation` is non-null, its `iterate` role equals this element's
   `iterate_next` exactly. *(New in 1.1 — 1.0 left the frame unbound to the
   update that produced it, admitting a trace whose evaluation belongs to a
   different depth: Reviewer D, F14.)*
7. The denominator of `update_relation` is non-zero. In the reference template
   this is enforced *in the generator* by an explicit raise, not an assert:
   0 occurrences in 20,000 seeds, smallest `|g_k − g_{k−1}|` observed 0.004.
8. `index` values are `1..len(elements)` in order. *(New in 1.1 — §5 required
   contiguity of the element index and §4 stated no analogue: Reviewer D, F14.)*
9. Every frame — preamble and evaluation alike — satisfies `frame_relations`
   (§4.6). *(New in 1.4: without it the node's answer is free.)*

### 4.4 Preamble

The question prescribes the starting trials, so the preamble is **not** free: a
trace whose preamble `iterate` values differ from the ones the question names has
not followed the stated scheme, and §7.2 marks that as a procedure failure even
when the final answer agrees.

Preamble frames are evaluation frames: their key set is `preamble_symbols`, which
for an `iteration` node equals `evaluation_symbols`.

**`preamble_binding` (normative, new in 1.2).** The preamble must be *bound to
the sequence*, not merely present beside it. It maps an element role to the index
of the preamble frame that role must equal, and the binding is checked against
`elements[0]`:

```json
{"iterate_prev":  {"frame": 0, "frame_role": "iterate"},
 "residual_prev": {"frame": 0, "frame_role": "residual"},
 "iterate_curr":  {"frame": 1, "frame_role": "iterate"},
 "residual_curr": {"frame": 1, "frame_role": "residual"}}
```

reads: element 1's `iterate_prev` is preamble frame 0's `iterate` role, and so
on. **`frame_role` is named explicitly**, never inferred from the element role's
name — 1.2 left it to a `iterate*`/`residual*` name prefix, which breaks for any
node whose roles are not spelled that way (Reviewer D2, F10). The field is
**required whenever `preamble` is present**. Through 1.1 the preamble was unconstrained, so
**replacing it with nonsense passed** while §4.4 asserted it was "not free" —
exactly the hole §4.3.6 had already closed one level down, missed one level up
(Reviewer D2, F4).

**What is still not checkable, deliberately:** the *internal* consistency of a
frame — that `A`, `P`, `AR` are the right functions of `y` — because the geometry
functions are the item's content and are not carried. A stated non-goal (§10),
not an oversight.

### 4.5 Tolerance headroom, stated

On the committed corpus the worst §4.3.1 residual is **4.985e-05** against a
5.0e-05 tolerance for a 4-dp `iterate_next` — 0.3% margin. This is expected and
benign: `update_relation`'s inputs are residuals already rounded to 3 dp, so the
recomputed value legitimately differs from the printed one by up to half a unit
in the last place. It is written down because a verifier author who tightened the
tolerance "for safety" would fail gold traces, and because a future template with
coarser intermediate rounding could exceed it (Reviewer D, F15).

### 4.6 `constants` and `frame_relations` — the frame's own arithmetic

An evaluation frame carries the quantities the update consumes. Through 1.3
nothing checked them, and §10 recorded that as a non-goal. **It was a hole, not a
non-goal:** §4.3.1 recomputes the update *from* the residuals, so an unchecked
frame makes the node's answer free. A trace reporting 3.501 m against a gold
1.380 m passed every clause and took full process credit.

`constants` declares the item's fixed quantities once; `frame_relations` is an
**ordered** list of `[symbol, expression]` giving each frame symbol's arithmetic:

```json
"constants": {"b": 4.6, "z": 2.0, "K": 6.619},
"frame_relations": [["A",  "(b + z * y) * y"],
                    ["P",  "b + 2 * y * sqrt(1 + z ** 2)"],
                    ["AR", "A * (A / P) ** (2 / 3)"],
                    ["g",  "round(AR) - K"]]
```

**Coverage (normative, and the clause the mechanism is worthless without).**
`frame_relations` must name **every member of `evaluation_symbols` except the
`iterate` role, exactly once**, in an order where each relation references only
`constants`, the iterate, and *earlier* symbols. §4.3.9 quantifies over the
frame's **symbols**, not over the relations present.

> Without this clause the mechanism is optional and therefore absent. 1.4
> quantified over the relations present, so **deleting the residual's relation
> left the other three vacuously checked and the residual free** — variant 5
> returns unchanged, by removing one line. A reviewer drove a trace to 3.492 m
> and to exactly 2.000 m against a gold 1.380 m, both passing clean. A mechanism
> a trace can opt out of is not a check.

**Grammar (normative).** §4.2's grammar, plus `**`, `sqrt(x)` and `round(x)`;
names resolve to `constants`, to the frame's `iterate` role, or to an *earlier*
symbol of `frame_relations`. The "earlier" rule makes cycles and forward
references structurally impossible rather than merely forbidden, and bounds
substitution depth by the relation count.

**`round` takes one argument.** `round(sym)` means `sym` at
`symbol_precision[sym]`, under the node's `rounding` mode. **The digit is not
restated.** 1.4 wrote `round(AR, 3)`, declaring a precision `symbol_precision`
already carried, with nothing requiring the two to agree — `round(AR, 1)` yields
a conforming but different answer. That is a §3.8 restatement violation inside
the section that enforces §3.8.

**Arithmetic errors are trace failures, not verifier crashes.** `sqrt` of a
negative, division by zero, and a domain error in `**` each **fail the frame**.
§6.0 makes hostile candidate traces first-class, so a verifier that raises on
one is not conforming.

**Evaluation (normative), and the distinction matters.** A bare symbol reference
means **that symbol's unrounded value** — its relation substituted in full, down
to the iterate and the constants. `round(sym)` means its **displayed** value.
This is how P2's round-then-recompute is expressed, and the two readings are not
interchangeable: `AR` recomputed from the frame's stored, *rounded* `A` and `P`
misses its stored value on **24.9% of the frames in the shipped conformance
corpus** (42 of 169), and on **6.6%** of preamble frames over 3,000 seeds — the
corpus rate is higher because mid-iteration frames carry larger `A` and `P` and
so more absolute rounding error. *(1.4 cited only the second figure without its
population, so a reader checking it against the shipped artefact got 3.7× the
number — Reviewer D2, round 5. Both populations are now named.)* A verifier that
binds symbols to stored values will reject a quarter of gold traces and conclude
the spec is wrong.

**Frame symbols carry no role indirection.** `roles` and `evaluation_roles` bind
the element and the two frame roles the update consumes; the remaining frame
symbols appear literally in `frame_relations`. Renaming them therefore means
rewriting the expression strings, not just rebinding a map. That is a deliberate
limit — the relations *are* expressions over those names — and it is the one
place the rename test of §10 requires editing data rather than a binding.

**Invariant §4.3.9.** For every frame — preamble and evaluation alike — and every
non-iterate member of `evaluation_symbols`, the symbol's relation, evaluated as
above and rounded to its `symbol_precision`, equals the frame's stored value.
Verified: **12,735 frames over 3,000 seeds, zero mismatches**; perturbing one
`AR` by 2.0 is rejected; and omitting a relation is now rejected at §6 step 1.

**What this does and does not buy.** It pins every frame to its iterate, so the
residuals are no longer free and neither is the converged answer. It does *not*
verify that the relations are the *right* physics — that `A = (b + zy)y` is
genuinely a trapezoidal area is a claim about the item, not about the trace, and
it is checked by the template's own grounding and by Reviewer B, not here.

---

## 5. `decision`

An ordered list of homogeneous *commitments*, each opened, filled by a stated
selection rule, and closed; terminated by exhausting a finite universe.

```
cardinality      : "answer_bearing"  (fixed for this type)
termination.kind : "exhaustion"      (fixed for this type)
```

### 5.1 Element, trial and candidate tables (normative)

Roles, not symbol names — the same status §4.1 has. Bind them through `roles`,
`trial_roles` and `candidate_roles` (§3.3). The third column shows the reference
template's own symbols, and **a verifier must never use that column**.

**Element** — one opened commitment:

| Role | Type | Meaning | Reference symbol |
|---|---|---|---|
| `index` | integer | 1-based, contiguous, ascending. | `station_id` |
| `budget_total` | number | The element's budget. | `capacity` |
| `budget_initial` | number | Budget at open. | `remaining_initial` |
| `trials` | array | The search log. **Never empty.** | `trials` |
| `committed` | array | The committed items, **in commitment order**. | `assigned` |
| `budget_final` | number | Budget at close. | `remaining_final` |

**Trial** — one application of the selection rule:

| Role | Type | Meaning | Reference symbol |
|---|---|---|---|
| `candidates` | array | Candidate snapshot *at the moment of choice*, in the deterministic scan order. Each entry has exactly the keys in `eligible_symbols`. | `eligible` |
| `chosen` | item \| `null` | The committed candidate, or `null` when none is admissible. | `chosen` |
| `budget_after` | number | Budget after this trial. | `remaining_after` |
| `closes` | boolean | **`closes == (chosen is null)`.** One rule, no exceptions. | `closes` |

**Candidate** — one entry of `candidates`:

| Role | Type | Meaning | Reference symbol |
|---|---|---|---|
| `item` | string | The candidate's identity; a member of `termination.universe`. | `task` |
| `measure` | number | The quantity `selection.criterion` optimises and the budget is spent in. | `duration` |
| `admissible` | boolean | Whether this candidate is selectable **at this budget**. Checked against `selection.filter_relation`, never trusted (§5.4.3c). | `fits` |

> **`item` is a declared role in 1.2.** 1.1 left the candidate's identity symbol
> undeclared, so a verifier had to derive it as *`eligible_symbols` minus the
> filter and criterion names* — sound only at exactly three members, and silently
> wrong at four (Reviewer D2, F5).

> **`closes`, restated.** Version 1.0 gave the rule twice, once in the trial
> table and once in §5.5 with an "except" clause whose scope read backwards
> (Reviewer D round 1, F11). There is now one rule: **a trial closes exactly when
> it chose nothing.** The final element's last trial may have an *empty*
> `candidates` — that is exhaustion, and it is still `chosen: null`, so it still
> closes. No exception is needed and none is granted. A closing trial is
> necessarily the **last** trial of its element (§5.4.9).

### 5.2 `selection` — the rule as data

```json
{"filter": "admissible", "filter_relation": "measure <= budget_before",
 "criterion": "measure", "objective": "max", "tie_break": "universe_order"}
```

| Field | Domain | Meaning |
|---|---|---|
| `filter` | a **candidate role** of boolean type | Only candidates where this is `true` are selectable. |
| `filter_relation` | a comparison over candidate roles and `budget_before` | **What makes `filter` true.** The verifier recomputes the flag from this rather than believing it. Grammar: a §4.2 arithmetic expression, one of `<= < >= > ==`, a second §4.2 expression. `budget_before` is bound to the trial's incoming budget. |
| `criterion` | a **candidate role** of numeric type | The quantity optimised, and the quantity the budget is spent in. |
| `objective` | `"max"` \| `"min"` | Direction. |
| `tie_break` | `"universe_order"` \| `"first"` | How equal `criterion` values are resolved. `universe_order` = earliest in `termination.universe`. |

> **Why `filter_relation` exists, and it is the most important addition in 1.2.**
> Through 1.1 the admissibility flag was **carried as data and never checked**,
> while §5.4.4 optimised only over the filtered set. A trace could therefore
> declare a candidate inadmissible when it was not, force the element to close,
> open an extra element, and **change `n` — which §2 and §7.3 make the answer** —
> while satisfying every clause of §5.3, §5.4, §6 and §8. Reviewer D2
> demonstrated it with two minimal nodes, an honest `n = 1` and a corrupt
> `n = 2`, both of which the 1.1 algorithm passed. A verifier that trusts a
> boolean the trace supplies is not verifying the thing that boolean decides.

`rule` is the prose restatement of all of this and is advisory (§3.2).

### 5.3 Global invariants (normative)

Across elements: `index` runs `1..len(elements)`; the `committed` lists are
pairwise disjoint; their union equals `termination.universe`; and
`max_elements >= len(universe)`.

### 5.4 Local invariants (normative)

For every element:

1. `budget_initial == budget_total == termination.budget`. *(1.4: `budget_total`
   was previously restated per element and checked against nothing, so a trace
   could understate it and the rule would honestly open an extra element —
   Reviewer D2, variant 4. The `capacity_constant` flag it replaced was no
   defence: the same wrong budget everywhere is constant.)*
2. Trials are consecutive: trial *j+1*'s incoming budget is trial *j*'s
   `budget_after`; trial 1's is `budget_initial`.
3. For a trial with a non-null `chosen`: `chosen` appears among `candidates` and
   its `admissible` role is `true`.
   **3b.** Every candidate has all its `termination.precedences` already
   committed in an earlier element or earlier in this one, **and** every
   uncommitted member of `universe` whose precedences *are* all satisfied appears
   among `candidates`. Both directions are required: an omission and an
   insertion are equally corrupting.
   **3c.** For **every** candidate, `admissible` equals the truth of
   `selection.filter_relation` evaluated with that candidate's roles and
   `budget_before` bound to the trial's incoming budget. *(New in 1.2 — this is
   the soundness fix of §5.2.)*
   **3d.** Every candidate's `measure` equals `termination.item_measures[item]`.
   *(New in 1.3, and it is §3.8's instance: without it a trace can inflate a
   measure in the trial where it wants the item excluded — leaving §5.4.3c
   honestly satisfied about a dishonest number — and restore it where it wants
   it committed, changing `n`, which is the answer, with every other invariant
   intact. Reviewer D2, F3′.)*
4. `chosen` optimises `criterion` over the admissible candidates under
   `objective`, with `tie_break` applied to equals.
5. `budget_after == budget_before − criterion(chosen)`, or `== budget_before`
   when `chosen` is `null`.
6. `committed` is exactly the ordered non-null `chosen` values of the element's
   trials.
7. `budget_final == budget_total − Σ criterion(committed)`.
8. A trial with `closes: true` is the **last** trial of its element, and every
   earlier trial has `closes: false`. *(New in 1.2: 1.1 defined `closes` but
   never said a closing trial ends the element — Reviewer D2, F9b.)*

### 5.5 Determinism

`eligible` must be produced by scanning a **defined order**, not a set. The
reference template scans the tuple `_T19_ORDER`. Combined with the five durations
being distinct integers by construction (measured: 0 duplicate-duration instances
in 4,000 seeds), the optimum is unique and `tie_break` is never exercised in this
item pool. It is still specified, because a verifier must be **total over traces
a model might produce**, and a model may well emit a tie.

---

## 6. Verifying a gold trace — normative algorithm

### 6.0 Modes — gold and candidate

**A node does not say whether it is gold or a candidate, and several clauses
branch on exactly that** (§8A.1 versus §3.4, §6 step 4 versus §7.2). The distinction is
therefore the **caller's**, passed in, not inferred: a conforming verifier takes a
mode argument.

| Mode | Meaning | Governing sections |
|---|---|---|
| `gold` | A node emitted by a template. Every §8 clause applies; `satisfied` must be `true`; the budget is an invariant. | §6, §8 |
| `candidate` | A node derived from a model's trace. §7 governs; a budget overrun is *reported*, not failed; `satisfied: false` is a legitimate observation. | §7 |

The algorithm below is the `gold` mode. Every numeric comparison uses §7.1; every
rounding uses `rounding`.

1. **Version and shape.** `schema_version` major version known; required fields
   present; `node_type` recognised; `cardinality` and `termination.kind` equal
   the values §4/§5 fix for that type; `cardinality_symbol` present iff
   `answer_bearing`; every conditional field required by the type present.
2. **Roles.** Every mandatory role of the type is a key of `roles`; every value
   of `roles` is a member of `element_symbols` (or, for `evaluation_roles`, of
   `evaluation_symbols`).
3. **Homogeneity.** Every element's key set equals `element_symbols` exactly — no
   missing key, no extra key. Likewise `preamble`/`preamble_symbols`,
   `trials`/`trial_symbols`, `eligible`/`eligible_symbols`, and every non-null
   `evaluation` against `evaluation_symbols`.
4. **Budget.** `1 <= len(elements) <= termination.max_elements`. *(1.1 resolves
   the 1.0 contradiction between §3's "may be empty" and its own lower bound in
   favour of **≥ 1**, and the §3 sentence is deleted: Reviewer D, F2. A
   zero-element sequence is not a trace of a computation, it is the absence of
   one; a candidate that emits one fails §7 rather than being verified.)*
5. **Carry.** For each *k* and each entry of `carry`, the value the path selects
   from element *k* equals the named symbol of element *k+1*. `carry_notes` is
   **not** checked and every entry in it must be reported on an `unchecked`
   channel that cannot contribute to a pass (§8B.10).
6. **Local invariants.** §4.3 or §5.4, per element; then §5.3 for `decision`.
7. **Termination.** The predicate holds on the last element and on no earlier
   one; `satisfied` agrees. For `convergence` the predicate is
   `quantity_role` compared against `tolerance` under `comparison`; for
   `exhaustion` it is the `scope` reading of `accumulator_symbol` against
   `universe`.
8. **Result.** `result.value` equals the `result.from` derivation (§3.6),
   rounded to `result.dp` under `rounding`.

> **`round()`'s mode.** 1.0's step 7 said "rounded" without saying how, while
> §9.2 elsewhere said half-up (Reviewer D, F9). The mode is now the node's
> `rounding` field, and both reference templates declare `decimal-half-up`. It is
> latent on gold — half-way instances are screened out at generation (D-016) —
> and live the moment a candidate trace is scored.

---

## 7. D3.4 — Comparator semantics for a model's trace

> **Revision note — amended by Phase 4 (D-054, SPEC-CHANGE 13).** This section
> shipped **unexercised**, and both Phase 3 schema reviewers named the missing
> corpus as the largest gap Phase 3 left. Phase 4 built it —
> [`phase4_conformance/candidates.json`](phase4_conformance/candidates.json),
> 16 candidate traces and 9 tolerance assertions, run by
> `python -m tests.trace_schema.candidate_7`, with every disposition this
> section states enumerated so coverage is *checked* rather than remembered.
>
> Running it found **four defects in this section**, not in a verifier:
>
> **F7-1 — §7.2's headline disposition has no instances.** *"A correct answer
> reached in four iterations instead of three is correct"* cannot occur for
> `normal_depth_iteration`. §4.3.3 ties `converged` to the node's own
> `termination.tolerance`, §8A.4 forbids `converged: true` before the last
> element, §4.3.1 recomputes every update, §4.6 recomputes every frame, and the
> preamble is fixed by the question. Together **these make the conforming node
> for a given question unique.** Verified both ways: clearing `converged` on the
> terminating element fails 4.3.3 and 8A.4; appending after it fails 8A.3 and
> 8A.4.
>
> D-038's measured spread of 1–5 updates over 4,000 seeds stands, but it is
> **between questions**, not between solvers of one question, and §7.2 conflates
> the two. *This is the six review rounds' bill arriving: each closed a way for
> a candidate to lie, and together they closed every way for a candidate to be
> differently right.* D4.9 reaches the same conclusion from the other direction
> and on most of the corpus — **141 of 150 templates emit a trace whose length
> is a constant**, so `cardinality` has one possible value there and both §7.2
> and §7.3 are vacuous.
>
> **F7-2 — §7 does not say whether §6 step 8 binds a candidate.** It does in
> every verifier built, and it should: a model whose stated answer contradicts
> its own steps has failed the process. The consequence is that **answer credit
> and process credit are not independent** for `iteration`, so "wrong answer,
> clean process" is not constructible. §7 is written as though they vary freely.
> **Read §7.2's rows as conditional on the trace being internally consistent.**
>
> **F7-3 — §7.4's direct solver cannot be an `iteration` node at all**, because
> a direct solve has no iteration. The tolerance §7.4 specifies is correct and
> testable, but only at the comparator's function level, where the conformance
> corpus now tests it with 9 assertions.
>
> **F7-4 — §7.3 row 3's set-versus-ordered distinction is vacuous** for a
> §6-conforming node: §5.4.6 already requires `committed` to equal the ordered
> chosen values, so a candidate whose set matches but whose order differs is
> rejected before §7.3 runs. The row describes a case the rest of the schema
> forbids.

The comparator scores a *candidate* trace against a *gold* node. It never
compares the two element lists positionally without first applying §7.2/§7.3.

### 7.1 Numeric tolerance

A candidate scalar matches gold when it agrees to **half a unit in the last place
gold displays it**, the boundary **inclusive**: the test is
`|candidate − gold| <= 0.5 × 10^(−p)` with `p = symbol_precision[symbol]`.
*(1.1 left inclusive-versus-exclusive unstated; it is latent on gold, where
exact-boundary instances are screened out at generation, and live the moment a
candidate is scored — Reviewer D2, F8.)* That place comes from `symbol_precision[symbol]`, or
`result.dp` for the result — **never inferred from the data, and always the
GOLD node's `symbol_precision`, never the candidate's.** `symbol_precision` is
the one problem-data entry that a *derived* check also consumes, so a candidate
declaring coarser precisions must not thereby widen the tolerance it is judged
against. Version 1.0 cited
"the display precision of `y`" in three clauses while declaring `dp` only on
`result`, where it is `3` against a 4-dp `y_next` — actively the wrong number for
the checks that cited it (Reviewer D, F4).

Instances that sit *exactly* on such a half-way boundary are not scored: they are
removed at generation time, because at a tie no rounding convention closes in
both directions and the instance has no defensible gold reading (D-016). Both
reference templates screen for this; measured rejection rates 2.30% and 0.65%.

See §4.5 for the headroom this leaves on the `iteration` update relation.

### 7.2 `cardinality: "incidental"` — the iteration-count-mismatch case

**A correct answer reached in four iterations instead of three is correct.**

| Candidate property | Disposition |
|---|---|
| `len(elements)` differs from gold | **Not a defect.** Not scored, not reported as an error. |
| `result.value` matches gold within §7.1 | **Answer credit.** |
| Every candidate element satisfies §4.3, **including §4.3.9's frame check** | **Process credit,** independently of the count. Through 1.3 the frame check did not exist, so a candidate emitting arithmetically consistent nonsense took full process credit — the mechanism rewarding exactly what it was built to detect. |
| Termination predicate holds on the candidate's last element | **Required for process credit.** A candidate that stopped early with a change above tolerance has not executed the stated scheme, even if its answer happens to match. |
| `len(elements) > max_elements` | **Reported, not failed.** The budget is the *generator's* guard; a solver taking six updates to the same depth has not made an error. |
| Preamble differs from the values the question prescribes | **Process failure.** The starting trials are part of the question. |

The count is reported as an observation (`gold 3, candidate 4`) so aggregate
analysis can see it, and is excluded from the score.

### 7.3 `cardinality: "answer_bearing"` — the station-count case

**A different element count is a different answer.**

| Candidate property | Disposition |
|---|---|
| `len(elements)` differs from gold | **Answer failure.** Report as `cardinality_mismatch` on `cardinality_symbol`, not as a formatting or style difference. |
| `len(elements)` matches; element `assigned` sets differ | **Process failure, answer may still be correct** — two different valid packings can share an `n`. Score the answer on `n`; score the process on the rule replay. |
| `assigned` compared | As **sets** for milestone credit; as **ordered lists** when replaying `selection`, whose output is ordered. |
| A trial violates §5.4.4 | **Process failure**, located at `(station_id, trial index)`. |

The asymmetry with §7.2 is the whole point of §2, and it is the one thing a
verifier built from a merged node type gets wrong.

### 7.4 When the question prescribes a method, the method's tolerance governs

An item whose question prescribes *how* to solve it has two defensible answers:
the value the prescribed scheme returns, and the value the underlying equation
has. For `normal_depth_iteration` they differ on **3.7% of instances**, always by
exactly 0.001 m — the scheme converges to ±0.002 m while the answer is quoted to
3 dp, so a solver who solves Manning's equation directly and *correctly* lands one
display unit away from gold.

**The gold answer is the prescribed method's output, and the comparator's answer
tolerance is the looser of §7.1's display tolerance and the method's own stated
tolerance** (`termination.tolerance`, in the result's units where the two are
commensurable). Marking a direct solver wrong for 0.001 m is the same error as
marking a four-update solver wrong for taking four updates, and §7.2 already
refuses the second.

Measured on both trees: **3.70% on `master`, 3.65% on the branch** — so this is
pre-existing, surfaced by Phase 3's Step 4 change rather than caused by it
(Reviewer B, F-B3). No extractor enforces this yet; carried as **R3-9**.

---

## 8. What a conforming verifier must REJECT

A verifier that accepts everything is not a verifier. **Every clause here is
reachable from the node alone**, and §8A is testable by mutating a trace while
§8B is not — §8B constrains the *verifier*, not the trace. 1.1 mixed the two, so
an implementer reading §8 as a list of traces skipped the one clause that is a
behavioural obligation (Reviewer D2).

### 8A — trace clauses (each constructible as a mutation)

1. `termination.satisfied == false` on a node in `gold` mode.
2. An element, preamble frame, trial, candidate entry or evaluation frame whose
   key set differs from its declared `*_symbols`, in either direction. *(This is
   the defect the flat milestone model had: symbols that appear and disappear
   between passes.)*
3. A `carry` path that does not hold between consecutive elements, or that is not
   well-formed under §3.5's grammar.
4. `iteration`: `converged: true` on any element other than the last;
   `evaluation` non-null on the terminating element; an evaluation frame whose
   `iterate` is not the element's `iterate_next`; a zero update denominator;
   non-contiguous `index`; a `preamble_binding` that does not hold against
   `elements[0]`; **a frame that does not satisfy `frame_relations` (§4.6)**; **`frame_relations` that does not cover every non-iterate
   `evaluation_symbols` member exactly once**.
5. `decision`: overlapping `committed` sets; a union that is not
   `termination.universe`; a `chosen` that does not optimise `selection`;
   **an `admissible` flag that disagrees with `selection.filter_relation`**;
   **a candidate `measure` that disagrees with `termination.item_measures`**;
   a `candidates` set that omits an item whose `precedences` are satisfied or
   includes one whose are not; non-contiguous `index`; a `budget_total` that is not `termination.budget`; a closing trial that is not last.
6. `cardinality: "answer_bearing"` with no `cardinality_symbol`, or
   `"incidental"` with one.
7. A role named mandatory by §4.1 or §5.1 but missing from the corresponding role
   map, or a role map value that is not in the corresponding `*_symbols`.
8. `schema_version` with an unknown major version.
9. A symbol used by any tolerance comparison that is absent from
   `symbol_precision`. A verifier must fail rather than infer a precision from
   the data.

### 8B — verifier-behaviour clauses (not constructible as a mutation)

10. **Silently passing a `carry_notes` entry.** Every entry must be reported on an
    `unchecked` channel that cannot contribute to a pass, so a node carrying one
    can never report a bare `PASS`. Treating an unchecked relation as a satisfied
    one is how a green suite comes to measure nothing (D-015, D-034). This is a
    statement about the verifier: all gold `decision` nodes carry `carry_notes`
    legitimately, so there is no corrupt trace to construct — test it as a
    behavioural assertion instead.
11. **Resolving any element value by a literal symbol name rather than through a
    role map.** Also untestable by mutation, and it is the clause the two
    reference templates exist to make checkable: rename every symbol in a node,
    rebind only its role maps, and a conforming verifier must still accept it —
    and must still reject it when a value is perturbed. §10 names this the
    **rename test**.

---

## 9. Worked examples

Both are lifted **verbatim, by script** from
[`phase3_conformance/traces.json`](phase3_conformance/traces.json), which carries
40 seeds of each template as `{seed, question, solution, trace_nodes}`.
Regenerate with `python -m tests.trace_schema.extract <out> <n>`. *(Version 1.0's
examples were transcribed by hand and one of them omitted a field the same
document marked required — they are now generated and machine-compared against
the corpus.)*

### 9.1 `iteration` — `normal_depth_iteration`, seed 0 (3 updates)

```json
{
 "schema_version": "1.5",
 "node_type": "iteration",
 "node_id": "t24_secant_normal_depth",
 "cardinality": "incidental",
 "roles": {
  "index": "k",
  "iterate_prev": "y_prev",
  "iterate_curr": "y_curr",
  "iterate_next": "y_next",
  "residual_prev": "g_prev",
  "residual_curr": "g_curr",
  "change": "change",
  "converged": "converged",
  "evaluation": "evaluation"
 },
 "update_relation": "iterate_curr - residual_curr * (iterate_curr - iterate_prev) / (residual_curr - residual_prev)",
 "termination": {
  "kind": "convergence",
  "quantity_role": "change",
  "comparison": "lt",
  "tolerance": 0.002,
  "predicate_prose": "abs(y_next - y_curr) < 0.002",
  "max_elements": 5,
  "satisfied": true
 },
 "rounding": "decimal-half-up",
 "constants": {
  "b": 2.1,
  "z": 2.5,
  "K": 6.619
 },
 "frame_relations": [
  [
   "A",
   "(b + z * y) * y"
  ],
  [
   "P",
   "b + 2 * y * sqrt(1 + z ** 2)"
  ],
  [
   "AR",
   "A * (A / P) ** (2 / 3)"
  ],
  [
   "g",
   "round(AR) - K"
  ]
 ],
 "symbol_precision": {
  "y_prev": 4,
  "y_curr": 4,
  "y_next": 4,
  "change": 4,
  "g_prev": 3,
  "g_curr": 3,
  "y": 4,
  "A": 3,
  "P": 3,
  "AR": 3,
  "g": 3
 },
 "preamble_binding": {
  "iterate_prev": {
   "frame": 0,
   "frame_role": "iterate"
  },
  "residual_prev": {
   "frame": 0,
   "frame_role": "residual"
  },
  "iterate_curr": {
   "frame": 1,
   "frame_role": "iterate"
  },
  "residual_curr": {
   "frame": 1,
   "frame_role": "residual"
  }
 },
 "preamble_symbols": [
  "y",
  "A",
  "P",
  "AR",
  "g"
 ],
 "evaluation_symbols": [
  "y",
  "A",
  "P",
  "AR",
  "g"
 ],
 "evaluation_roles": {
  "iterate": "y",
  "residual": "g"
 },
 "element_symbols": [
  "k",
  "y_prev",
  "g_prev",
  "y_curr",
  "g_curr",
  "y_next",
  "change",
  "converged",
  "evaluation"
 ],
 "carry": {
  "y_prev": "y_curr",
  "g_prev": "g_curr",
  "y_curr": "y_next",
  "g_curr": "evaluation.g"
 },
 "carry_notes": {},
 "preamble": [
  {
   "y": 1.0,
   "A": 4.6,
   "P": 7.485,
   "AR": 3.325,
   "g": -3.294
  },
  {
   "y": 1.5,
   "A": 8.775,
   "P": 10.178,
   "AR": 7.949,
   "g": 1.33
  }
 ],
 "elements": [
  {
   "k": 1,
   "y_prev": 1.0,
   "g_prev": -3.294,
   "y_curr": 1.5,
   "g_curr": 1.33,
   "y_next": 1.3562,
   "change": 0.1438,
   "converged": false,
   "evaluation": {
    "y": 1.3562,
    "A": 7.446,
    "P": 9.403,
    "AR": 6.373,
    "g": -0.246
   }
  },
  {
   "k": 2,
   "y_prev": 1.5,
   "g_prev": 1.33,
   "y_curr": 1.3562,
   "g_curr": -0.246,
   "y_next": 1.3786,
   "change": 0.0224,
   "converged": false,
   "evaluation": {
    "y": 1.3786,
    "A": 7.646,
    "P": 9.524,
    "AR": 6.605,
    "g": -0.014
   }
  },
  {
   "k": 3,
   "y_prev": 1.3562,
   "g_prev": -0.246,
   "y_curr": 1.3786,
   "g_curr": -0.014,
   "y_next": 1.38,
   "change": 0.0014,
   "converged": true,
   "evaluation": null
  }
 ],
 "result": {
  "symbol": "yn",
  "value": 1.38,
  "unit": "m",
  "dp": 3,
  "from": "last.iterate_next"
 }
}
```

Reading element 2 against element 1 with no knowledge of the count, and using
only `roles`: `carry` gives `y_prev(2) = y_curr(1) = 1.5` ✓,
`g_prev(2) = g_curr(1) = 1.33` ✓, `y_curr(2) = y_next(1) = 1.3562` ✓,
`g_curr(2) = evaluation.g(1) = −0.246` ✓. Evaluating `update_relation` with
`iterate_curr = 1.3562`, `iterate_prev = 1.5`, `residual_curr = −0.246`,
`residual_prev = 1.33` gives `1.378646…`, which matches `iterate_next = 1.3786`
within the 5.0e-05 tolerance `symbol_precision["y_next"] = 4` implies ✓.
`change = |1.3786 − 1.3562| = 0.0224 ≥ 0.002`, so `converged: false` ✓,
`evaluation` is non-null ✓ and its `y` equals `y_next` ✓.

### 9.2 `decision` — `line_balancing_heuristic`, seed 35 (3 stations)

```json
{
 "schema_version": "1.4",
 "node_type": "decision",
 "node_id": "t19_greedy_station_assignment",
 "cardinality": "answer_bearing",
 "cardinality_symbol": "n",
 "termination": {
  "kind": "exhaustion",
  "scope": "cumulative",
  "accumulator_role": "committed",
  "predicate_prose": "every task in `universe` has been assigned to exactly one element",
  "universe": [
   "a",
   "b",
   "c",
   "d",
   "e"
  ],
  "precedences": {
   "a": [],
   "b": [
    "a"
   ],
   "c": [
    "a"
   ],
   "d": [
    "b",
    "c"
   ],
   "e": [
    "d"
   ]
  },
  "item_measures": {
   "a": 51,
   "b": 70,
   "c": 47,
   "d": 87,
   "e": 22
  },
  "budget": 139,
  "max_elements": 5,
  "satisfied": true
 },
 "roles": {
  "index": "station_id",
  "budget_total": "capacity",
  "budget_initial": "remaining_initial",
  "trials": "trials",
  "committed": "assigned",
  "budget_final": "remaining_final"
 },
 "trial_roles": {
  "candidates": "eligible",
  "chosen": "chosen",
  "budget_after": "remaining_after",
  "closes": "closes"
 },
 "candidate_roles": {
  "item": "task",
  "measure": "duration",
  "admissible": "fits"
 },
 "selection": {
  "filter": "admissible",
  "filter_relation": "measure <= budget_before",
  "criterion": "measure",
  "objective": "max",
  "tie_break": "universe_order"
 },
 "rule": "among unassigned tasks whose predecessors are all assigned AND whose duration fits the remaining time, take the longest; ties by _T19_ORDER position; when none fits, close the station",
 "rounding": "decimal-half-up",
 "symbol_precision": {
  "capacity": 0,
  "remaining_initial": 0,
  "remaining_final": 0,
  "remaining_after": 0,
  "duration": 0,
  "station_id": 0
 },
 "element_symbols": [
  "station_id",
  "capacity",
  "remaining_initial",
  "trials",
  "assigned",
  "remaining_final"
 ],
 "trial_symbols": [
  "eligible",
  "chosen",
  "remaining_after",
  "closes"
 ],
 "eligible_symbols": [
  "task",
  "duration",
  "fits"
 ],
 "carry": {},
 "carry_notes": {
  "assigned": "union of every element's assigned",
  "remaining_initial": "capacity (reset each element)"
 },
 "elements": [
  {
   "station_id": 1,
   "capacity": 139,
   "remaining_initial": 139,
   "trials": [
    {
     "eligible": [
      {
       "task": "a",
       "duration": 51,
       "fits": true
      }
     ],
     "chosen": "a",
     "remaining_after": 88,
     "closes": false
    },
    {
     "eligible": [
      {
       "task": "b",
       "duration": 70,
       "fits": true
      },
      {
       "task": "c",
       "duration": 47,
       "fits": true
      }
     ],
     "chosen": "b",
     "remaining_after": 18,
     "closes": false
    },
    {
     "eligible": [
      {
       "task": "c",
       "duration": 47,
       "fits": false
      }
     ],
     "chosen": null,
     "remaining_after": 18,
     "closes": true
    }
   ],
   "assigned": [
    "a",
    "b"
   ],
   "remaining_final": 18
  },
  {
   "station_id": 2,
   "capacity": 139,
   "remaining_initial": 139,
   "trials": [
    {
     "eligible": [
      {
       "task": "c",
       "duration": 47,
       "fits": true
      }
     ],
     "chosen": "c",
     "remaining_after": 92,
     "closes": false
    },
    {
     "eligible": [
      {
       "task": "d",
       "duration": 87,
       "fits": true
      }
     ],
     "chosen": "d",
     "remaining_after": 5,
     "closes": false
    },
    {
     "eligible": [
      {
       "task": "e",
       "duration": 22,
       "fits": false
      }
     ],
     "chosen": null,
     "remaining_after": 5,
     "closes": true
    }
   ],
   "assigned": [
    "c",
    "d"
   ],
   "remaining_final": 5
  },
  {
   "station_id": 3,
   "capacity": 139,
   "remaining_initial": 139,
   "trials": [
    {
     "eligible": [
      {
       "task": "e",
       "duration": 22,
       "fits": true
      }
     ],
     "chosen": "e",
     "remaining_after": 117,
     "closes": false
    },
    {
     "eligible": [],
     "chosen": null,
     "remaining_after": 117,
     "closes": true
    }
   ],
   "assigned": [
    "e"
   ],
   "remaining_final": 117
  }
 ],
 "result": {
  "symbol": "n",
  "value": 3,
  "dp": 0,
  "from": "len(elements)"
 }
}
```

Element 1, trial 2 checks the rule rather than trusting it: candidates are `b` (70)
and `c` (47); `filter_relation` recomputes `admissible` as `measure <= budget_before`
— `70 <= 88` ✓ and `47 <= 88` ✓, matching both flags — and `selection` says `max`
on `measure`, so `b` ✓;
`88 − 70 = 18` ✓. Trial 3: `filter_relation` gives `47 <= 18` = false,
which is what `admissible` claims ✓, so nothing is selectable, `chosen: null`,
`closes: true`, budget unchanged, and it is the element's last trial ✓. `assigned = [a, b]`,
`139 − 51 − 70 = 18 = remaining_final` ✓. Precedences are now checkable: at
trial 1 only `a` has no prerequisites, so `eligible` correctly holds `a` alone;
at trial 2 `b` and `c` both have `a` committed, and `d` does not yet have both
`b` and `c`, so it is correctly absent ✓. Across elements the union is
`{a,b,c,d,e}` = `universe` ✓, disjoint ✓, `n = 3` ✓ — and `n = 3` is the
**answer**, feeding `(3·139 − 277)/(3·139)·100 = 140/417·100 = 33.573…`, quoted
half-up at one decimal as `33.6%`.

Note the last element's final trial has an **empty** `eligible`: that is
exhaustion, and under §5.1's single `closes` rule it needs no exception — nothing
was chosen, so it closes.

---

## 10. Limits, and what this specification does NOT settle

Stated as non-goals so an implementer meets the edge here rather than in the
middle of a verifier.

- **~~Frame-internal physics is not checkable.~~ WITHDRAWN in 1.4.** Versions
  1.0–1.3 recorded this as a non-goal, reasoning that checking a frame would mean
  "carrying an expression language for arbitrary engineering formulae, which is
  the milestone model's problem, not this node's". **That was a limitation
  reported as a principle, and it was load-bearing for an attacker**: because
  §4.3.1 recomputes the update *from* the residuals, an unchecked frame left the
  node's answer entirely free. The node already carried two expression languages
  (§4.2, §5.2), so the stated cost was overstated by roughly one grammar rule.
  §4.6 now carries the frame's arithmetic and §4.3.9 checks it. What remains out
  of scope is narrower and genuinely is a non-goal: the node cannot verify that
  the declared relations are the *right physics*, only that the frame obeys them.
- **Prose/node agreement is out of scope.** 1.0 made "every numeric token the
  node carries appears in the prose" step 8 of the normative algorithm while
  giving no definition of "an arithmetic line" and no precision to compare at
  (Reviewer D, F10). It is deleted from §6 rather than left unimplementable. It
  is also vacuous by construction today: both reference templates render the
  prose *from* the node (D-039), which is a stronger guarantee than any check.
  It becomes a real requirement the day a template stops doing that.
- **No dimensional checking.** Only `result` carries `unit`. A dimensional
  comparator needs units on every symbol; that is a milestone-model decision.
- **§8B.11 cannot be audited statically, only by the rename test.** A verifier
  that genuinely hardcoded `element["trials"]` would pass the entire reference
  corpus and fail only on a renamed node, because for three roles (`trials`,
  `chosen`, `closes`) the role name and the reference template's symbol name are
  the *same string*. Grepping a verifier for reference symbols therefore returns
  false positives and, worse, would return false negatives on a real violation.
  Run the rename test; do not read the source (Reviewer D2, F13).
- **The rename test is the conformance standard for a node type.** Take a node,
  rename every symbol, rebind only its role maps, leave `update_relation` /
  `selection` untouched, and a conforming verifier must accept it — and must
  still reject it when one value is perturbed. `iteration` passes this
  (independently, on a synthetic node with a *different* recurrence and nine
  elements). `decision` is written to the same standard in 1.2 but **has not been
  exercised against a second template**, because none exists.
- **No second `decision` template exists in the corpus** (D-055). Phase 4's D4.8
  scanned all 150; 19 carry decision-shaped vocabulary and none has the shape. A
  `decision` node needs three things at once — an ordered list of commitments, a
  shared budget they consume, and a count that is part of the answer — and
  `line_balancing_heuristic` is the only template with all three. Closing this
  limitation needs a **new item**, which is an item-design decision.
- **`precedences` is exercised against exactly one DAG.** All 40 `decision`
  traces carry the same five-task network, so §5.4.3b is verified thoroughly
  against one relation and not at all against a second shape. A second
  `decision` template is the evidence that would settle it.
- **~~Two node types are specified; two templates implement them.~~ MEASURED by
  Phase 4, and the answer is negative (D-055).** Both named candidates were
  fitted by declaring a binding only, against an unmodified verifier
  (`python -m tests.trace_schema.d4_7_third_iteration`). Both were **rejected**,
  for different reasons. `linear_reservoir_routing_step`: §4 fixes
  `termination.kind` to `"convergence"` and routing *exhausts a two-interval
  hydrograph* instead; its element count is neither `incidental` nor
  `answer_bearing`, so D-038's rule has no verdict for it; and `frame_relations`
  admits only constants, the iterate and earlier frame symbols, with nowhere to
  put the per-element inflow pair. `qr_policy_one_iteration`: iterates on a
  **pair** coupled through `n(R)`, and §4.1 declares exactly one iterate triple.
  **Both are retired as generalisation targets.** `iteration` is not "a repeated
  sub-chain" — it is *a convergence-terminated refinement of one quantity driven
  by its own residual*, and that is its real scope.
- **The node is a frame local, not a return value.** The templates' public
  contract is still `(question, solution)`, and
  `tests/trace_schema/extract.py` lifts `trace_nodes` via `sys.settrace`.
  Deliberate — changing the return contract is a corpus-wide decision, not a
  two-template one — but it means the node is not yet part of any interface and
  no check gates its shape. Recorded as **D-039**.
