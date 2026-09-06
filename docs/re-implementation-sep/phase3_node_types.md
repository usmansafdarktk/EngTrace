# D3.3 — The `iteration` and `decision` trace node types

**Status:** normative for the milestone model · **Phase:** 3 · **Date:** 2026-09-06
**Companion:** [`template_redesign_spec.md`](template_redesign_spec.md) §3.3 ·
[`phase3_summary.md`](phase3_summary.md) · conformance corpus in
[`phase3_conformance/traces.json`](phase3_conformance/traces.json)
**Reference implementations:** the `trace_nodes` local of
`template_normal_depth_iteration`
(`data/templates/branches/civil_engineering/water_resources/uniform_flow.py`)
and of `template_line_balancing_heuristic`
(`.../industrial_engineering/production_and_inventory/production_planning.py`).
**Extractor:** `python -m tests.trace_schema.extract <out.json> <n_seeds>`.

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
is the comparator, and §8 lists the traces a conforming verifier must *reject*.

---

## 2. The finding that shapes this specification

`template_redesign_spec.md` §3.2 says the two node types "are the same shape (an
ordered list of homogeneous sub-traces with stable within-element symbols and a
termination predicate)" and recommends specifying them together.

**The structural claim is right and the equivalence is wrong.** Specifying them
as one type produces a verifier that is wrong on one of the two templates,
because the two differ on a property that decides how they are scored:

> **Is the cardinality of the sequence an observable of the answer?**

- `normal_depth_iteration`: the balance of the answer is the converged depth
  `yn`. The update count does not enter it. Worse, the count is **not stable
  under permitted numerical slack** — the termination test is
  `|y_next − y_curr| < 0.002 m` on values rounded to 4 dp, so a solver who
  carries more digits can converge in a different number of updates and reach
  the same depth. Marking that solver wrong would be marking them wrong for
  being more careful.
- `line_balancing_heuristic`: the answer is
  `balance delay = (n·CT − Σt)/(n·CT)·100`, and `n` **is the element count**.
  The count is also *exactly* determined — every fit decision is an integer
  comparison, so no rounding slack can move it. A solver who opens four stations
  where gold opened three has not taken a different route to the same answer;
  they have a different answer.

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

---

## 3. Common base — the `sequence` node

Both types are objects with these fields. Types are JSON types; `‹T›` means a
value of the template's own numeric or symbolic domain.

| Field | Type | Req. | Meaning |
|---|---|:--:|---|
| `node_type` | `"iteration"` \| `"decision"` | ✔ | Selects the rules in §4 or §5. |
| `node_id` | string | ✔ | Stable within a template; the anchor a comparator reports against. Must not encode the element count. |
| `cardinality` | `"incidental"` \| `"answer_bearing"` | ✔ | §2. Determines whether §7.2 or §7.3 applies. |
| `cardinality_symbol` | string | ✖ | Required **iff** `cardinality == "answer_bearing"`: the symbol the count is bound to in the answer (`"n"`). Absent otherwise. |
| `termination` | object | ✔ | §3.1. |
| `element_symbols` | array of string | ✔ | The symbol set every element carries. **Stable across elements** — this is what replaces the flat model's per-step `{id, symbol}`. |
| `elements` | array of object | ✔ | Ordered, homogeneous. Each has exactly the keys in `element_symbols`. May be empty only if `termination.satisfied` is true at entry. |
| `carry` | object | ✔ | §3.2. Maps a symbol of element *k+1* to the expression over element *k* that produces it. `{}` for a node whose elements are independent. |
| `result` | object | ✔ | `{symbol, value, unit?, dp?}` — the single quantity the node contributes downstream. `dp` is the decimal places it is displayed and consumed at (P2). |
| `preamble` | array of object | ✖ | Fixed frames the question prescribes *before* the sequence starts. Present on `iteration` where the question names starting values. |
| `preamble_symbols` | array of string | ✖ | Required iff `preamble` is present. |
| `rule` | string | ✖ | Required on `decision`: the selection rule in prose, including its tie-break. |
| `trial_symbols` | array of string | ✖ | Required on `decision`: the symbol set of the nested trial records (§5). |

### 3.1 `termination`

| Field | Type | Req. | Meaning |
|---|---|:--:|---|
| `kind` | `"convergence"` \| `"exhaustion"` | ✔ | Which of §4/§5 applies. |
| `quantity` | string | ✔ | The symbol the predicate is evaluated on. |
| `predicate` | string | ✔ | The predicate in the template's own notation, evaluable over element symbols. |
| `satisfied` | boolean | ✔ | Whether the emitted sequence actually terminated by the predicate. **A gold trace with `satisfied: false` is malformed** and a verifier must reject it (§8.1). |
| `tolerance` | number | cond. | Required iff `kind == "convergence"`. |
| `max_elements` | integer | ✔ | The budget. `len(elements) <= max_elements` is an invariant, not a hope. |
| `universe` | array | cond. | Required iff `kind == "exhaustion"`: the finite set being exhausted. |

**On `satisfied`.** In a gold trace it is always `true`, because the template
raises rather than emitting an unconverged trace — deliberately a `raise` and not
an `assert`, since `python -O` strips asserts and the prose tells the reader the
tolerance was met (Phase 2 Reviewer C, C-6). The field exists because a *model's*
trace may legitimately carry `false`, and the comparator must be able to say so
rather than crash.

### 3.2 `carry`

`carry` is what makes a variable-length sequence checkable without replaying it.
It maps each symbol of element *k+1* to a path expression over element *k*:

```json
{"y_prev": "y_curr", "g_prev": "g_curr",
 "y_curr": "y_next", "g_curr": "evaluation.g"}
```

reads: *element k+1's `y_prev` is element k's `y_curr`*, and so on. A dotted path
descends into a nested object. The right-hand side may also be a prose
description when the relation is not a projection (`decision` uses
`"union of every element's assigned"`); a verifier treats a non-path right-hand
side as **advisory and unchecked**, and §8.4 requires it to say so rather than
silently pass.

This is the property that lets a verifier check element *k* against element *k−1*
in isolation. It is why the count can vary and the trace stay machine-checkable.

---

## 4. `iteration`

An ordered list of homogeneous refinement steps, terminated by a numeric
predicate on a converging quantity.

```
cardinality      : "incidental"      (fixed for this type)
termination.kind : "convergence"     (fixed for this type)
```

**Element.** One element is one *update*, plus the evaluation that the update's
output feeds — so that an element is self-contained and the `carry` map is a pure
projection. `evaluation` is `null` on the terminating element, because a
converged iterate is not re-evaluated.

**Local invariant (normative).** For every element:

1. The update relation holds among the element's own symbols, to the tolerance
   in §7.1. For the reference template that relation is the secant step
   `y_next = y_curr − g_curr·(y_curr − y_prev)/(g_curr − g_prev)`.
2. `change == |y_next − y_curr|`, at the display precision of `y`.
3. `converged == (change < termination.tolerance)`.
4. `converged` is `true` on the **last** element and `false` on every other one.
5. `evaluation` is `null` iff `converged` is `true`.
6. The denominator of the update relation is non-zero. In the reference template
   this is enforced *in the generator* by an explicit raise, not an assert:
   0 occurrences in 20,000 seeds, smallest `|g_k − g_{k−1}|` observed 0.004.

**Preamble.** The question prescribes the two starting trial depths, so the
preamble is **not** free: a trace whose preamble `y` values differ from the ones
the question names has not followed the stated scheme, and §7.2 marks that as a
procedure failure even when the final answer agrees.

---

## 5. `decision`

An ordered list of homogeneous *commitments*, each opened, filled by a stated
selection rule, and closed; terminated by exhausting a finite universe.

```
cardinality      : "answer_bearing"  (fixed for this type)
termination.kind : "exhaustion"      (fixed for this type)
```

**Element** — one opened workstation:

| Symbol | Type | Meaning |
|---|---|---|
| `station_id` | integer | 1-based, contiguous, ascending. |
| `capacity` | number | The station's budget (`CT`). |
| `remaining_initial` | number | Budget at open. Equals `capacity` for this rule. |
| `trials` | array | The search log; see below. Never empty. |
| `assigned` | array | The committed set, **in assignment order**. |
| `remaining_final` | number | Budget at close. |

**Trial** — one application of the selection rule (`trial_symbols`):

| Symbol | Type | Meaning |
|---|---|---|
| `eligible` | array of `{task, duration, fits}` | The candidate snapshot *at the moment of choice*, in the deterministic scan order. `fits` is the capacity test, reported separately from the choice so a verifier can check the rule rather than trust it. |
| `chosen` | string \| `null` | The committed candidate, or `null` when none fits. |
| `remaining_after` | number | Budget after this trial. |
| `closes` | boolean | `true` exactly on the trial where `chosen` is `null`. |

**Local invariants (normative).** For every element:

1. `remaining_initial == capacity`.
2. Trials are consecutive: trial *j+1*'s budget in is trial *j*'s
   `remaining_after`.
3. For a trial with `chosen == c`: `c ∈ eligible`, `eligible[c].fits` is true,
   and `duration(c) == max{duration(x) : x ∈ eligible, x.fits}` — with ties
   broken by `rule`. **The rule is checked, not assumed.**
4. `remaining_after == remaining_before − duration(chosen)`, or
   `== remaining_before` when `chosen is null`.
5. `closes` is `true` on the last trial and false on every other, **except** on
   the final element where the universe is exhausted: there the last trial may
   have an empty `eligible` and `closes: true`.
6. `assigned` is exactly the ordered `chosen` values of the element's trials.
7. `remaining_final == capacity − Σ duration(assigned)`.

**Global invariants.** Across elements: `station_id` is `1..len(elements)`;
the `assigned` lists are pairwise disjoint and their union is
`termination.universe`; and every task's precedences are satisfied by an earlier
or same-element assignment.

**Determinism.** `eligible` must be produced by scanning a **defined order**, not
a set. The reference template scans the tuple `_T19_ORDER`. Combined with the
five durations being distinct integers by construction (measured: 0
duplicate-duration instances in 4,000 seeds), the `max` is unique and the
tie-break is never exercised in this item pool. It is still specified, because a
verifier must be **total over traces a model might produce**, and a model may
well emit a tie.

---

## 6. Verifying a gold trace — normative algorithm

A conforming implementation, given a node and nothing else, performs in order:

1. **Shape.** Required fields present; `node_type` recognised; `cardinality`
   consistent with `node_type` (§4/§5 fix it); `cardinality_symbol` present iff
   `answer_bearing`.
2. **Homogeneity.** Every element's key set equals `element_symbols` exactly —
   no missing key, no extra key. Same for `preamble`/`preamble_symbols` and
   `trials`/`trial_symbols`.
3. **Budget.** `1 <= len(elements) <= termination.max_elements`.
4. **Carry.** For each *k*, and each path-valued entry of `carry`, the value at
   element *k+1* equals the value the path selects from element *k*. Non-path
   entries are reported as unchecked (§8.4).
5. **Local invariants.** §4 or §5 as applicable, per element.
6. **Termination.** The predicate holds on the last element and on no earlier
   one; `satisfied` agrees.
7. **Result.** `result.value` is derivable from the last element — for
   `iteration`, `result.value == round(last.y_next, result.dp)`; for `decision`,
   `result.value == len(elements)`.
8. **Prose agreement** *(when the printed solution is supplied alongside)*: every
   numeric token the node carries appears in the prose at its stated precision,
   and no arithmetic line in the prose is absent from the node. The reference
   templates render the prose **from** the node, so this check should be
   vacuous; it exists to detect the case where that stops being true.

---

## 7. D3.4 — Comparator semantics for a model's trace

The comparator scores a *candidate* trace against a *gold* node. It never
compares the two element lists positionally without first applying §7.2/§7.3.

### 7.1 Numeric tolerance

A candidate scalar matches gold when it agrees to **half a unit in the last place
gold displays it** (`result.dp`, or the element symbol's display precision) —
the same rule T1 applies to printed-arithmetic closure. Instances that sit
*exactly* on such a half-way boundary are not scored: they are removed at
generation time, because at a tie no rounding convention closes in both
directions and the instance has no defensible gold reading (D-016). Both
reference templates now screen for this; measured rejection rates 2.30% and
0.65%.

### 7.2 `cardinality: "incidental"` — the iteration-count-mismatch case

**A correct answer reached in four iterations instead of three is correct.**
Normatively:

| Candidate property | Disposition |
|---|---|
| `len(elements)` differs from gold | **Not a defect.** Not scored, not reported as an error. |
| `result.value` matches gold within §7.1 | **Answer credit.** |
| Every candidate element satisfies its §4 local invariants | **Process credit,** independently of the count. |
| Termination predicate holds on the candidate's last element | **Required for process credit.** A candidate that stopped early with a change above tolerance has not executed the stated scheme, even if its answer happens to match. |
| `len(elements) > termination.max_elements` | **Reported, not failed.** The budget is the *generator's* guard; a solver taking six updates to the same depth has not made an error. |
| Preamble differs from the values the question prescribes | **Process failure.** The starting trials are part of the question. |

The count is reported as an observation (`gold 3, candidate 4`) so that
aggregate analysis can see it, and is excluded from the score.

### 7.3 `cardinality: "answer_bearing"` — the station-count case

**A different element count is a different answer.** Normatively:

| Candidate property | Disposition |
|---|---|
| `len(elements)` differs from gold | **Answer failure.** Report as `cardinality_mismatch` on `cardinality_symbol`, not as a formatting or style difference. |
| `len(elements)` matches; element `assigned` sets differ | **Process failure, answer may still be correct** — two different valid packings can share an `n`. Score the answer on `n` and the delay; score the process on the rule replay. |
| `assigned` compared | As **sets** for milestone credit; as **ordered lists** when replaying the rule, since the rule's own output is ordered. |
| A trial violates the selection rule (§5.3) | **Process failure**, located at `(station_id, trial index)`. |

The asymmetry with §7.2 is the whole point of §2, and it is the one thing a
verifier built from a merged node type gets wrong.

---

## 8. What a conforming verifier must REJECT

A verifier that accepts everything is not a verifier. These are the negative
cases; a conforming implementation rejects each, with the cited reason.

1. `termination.satisfied == false` on a **gold** node — the generator's
   contract is that it raises instead.
2. An element whose key set differs from `element_symbols`, in either direction.
   (This is the defect the flat milestone model had: symbols that appear and
   disappear between passes.)
3. A `carry` path that does not hold between consecutive elements.
4. Silently passing a non-path `carry` right-hand side. The verifier must report
   it as **unchecked**; treating an unchecked relation as a satisfied one is how
   a green suite comes to measure nothing (D-015, D-034).
5. `iteration`: `converged: true` on any element other than the last;
   `evaluation` non-null on the terminating element; a zero update denominator.
6. `decision`: overlapping `assigned` sets; a union that is not
   `termination.universe`; a `chosen` that is not the longest fitting eligible
   candidate; a precedence violated at the moment of choice; non-contiguous
   `station_id`.
7. `cardinality: "answer_bearing"` with no `cardinality_symbol`, or
   `"incidental"` with one.

---

## 9. Worked examples

Both are lifted verbatim from
[`phase3_conformance/traces.json`](phase3_conformance/traces.json), which carries
40 seeds of each template as `{seed, question, solution, trace_nodes}`.
Regenerate with `python -m tests.trace_schema.extract <out> <n>`.

### 9.1 `iteration` — `normal_depth_iteration`, seed 0 (3 updates)

```json
{"node_type": "iteration", "node_id": "t24_secant_normal_depth",
 "cardinality": "incidental",
 "termination": {"kind": "convergence", "quantity": "y",
                 "predicate": "abs(y_next - y_curr) < tol",
                 "tolerance": 0.002, "max_elements": 5, "satisfied": true},
 "preamble_symbols": ["y", "A", "P", "AR", "g"],
 "element_symbols": ["k", "y_prev", "g_prev", "y_curr", "g_curr",
                     "y_next", "change", "converged", "evaluation"],
 "carry": {"y_prev": "y_curr", "g_prev": "g_curr",
           "y_curr": "y_next", "g_curr": "evaluation.g"},
 "preamble": [{"y": 1.0, "A": 4.6, "P": 7.485, "AR": 3.325, "g": -3.294},
              {"y": 1.5, "A": 8.775, "P": 10.178, "AR": 7.949, "g": 1.33}],
 "elements": [
  {"k": 1, "y_prev": 1.0, "g_prev": -3.294, "y_curr": 1.5, "g_curr": 1.33,
   "y_next": 1.3562, "change": 0.1438, "converged": false,
   "evaluation": {"y": 1.3562, "A": 7.446, "P": 9.403, "AR": 6.373, "g": -0.246}},
  {"k": 2, "y_prev": 1.5, "g_prev": 1.33, "y_curr": 1.3562, "g_curr": -0.246,
   "y_next": 1.3786, "change": 0.0224, "converged": false,
   "evaluation": {"y": 1.3786, "A": 7.646, "P": 9.524, "AR": 6.605, "g": -0.014}},
  {"k": 3, "y_prev": 1.3562, "g_prev": -0.246, "y_curr": 1.3786,
   "g_curr": -0.014, "y_next": 1.38, "change": 0.0014, "converged": true,
   "evaluation": null}],
 "result": {"symbol": "yn", "value": 1.38, "unit": "m", "dp": 3}}
```

Reading element 2 against element 1 with no knowledge of the count:
`carry` gives `y_prev(2) = y_curr(1) = 1.5` ✓, `g_prev(2) = g_curr(1) = 1.33` ✓,
`y_curr(2) = y_next(1) = 1.3562` ✓, `g_curr(2) = evaluation.g(1) = −0.246` ✓.
The secant step:
`1.3562 − (−0.246)(1.3562 − 1.5)/((−0.246) − 1.33) = 1.378646…` → `1.3786` ✓.
`change = |1.3786 − 1.3562| = 0.0224 ≥ 0.002`, so `converged: false` ✓ and
`evaluation` is non-null ✓.

### 9.2 `decision` — `line_balancing_heuristic`, seed 35 (3 stations)

```json
{"node_type": "decision", "node_id": "t19_greedy_station_assignment",
 "cardinality": "answer_bearing", "cardinality_symbol": "n",
 "termination": {"kind": "exhaustion", "quantity": "assigned",
                 "predicate": "len(assigned) == len(tasks)",
                 "universe": ["a", "b", "c", "d", "e"],
                 "max_elements": 5, "satisfied": true},
 "rule": "among unassigned tasks whose predecessors are all assigned AND whose duration fits the remaining time, take the longest; ties by _T19_ORDER position; when none fits, close the station",
 "element_symbols": ["station_id", "capacity", "remaining_initial",
                     "trials", "assigned", "remaining_final"],
 "trial_symbols": ["eligible", "chosen", "remaining_after", "closes"],
 "carry": {"assigned": "union of every element's assigned",
           "remaining_initial": "capacity (reset each element)"},
 "elements": [
  {"station_id": 1, "capacity": 139, "remaining_initial": 139,
   "trials": [
    {"eligible": [{"task": "a", "duration": 51, "fits": true}],
     "chosen": "a", "remaining_after": 88, "closes": false},
    {"eligible": [{"task": "b", "duration": 70, "fits": true},
                  {"task": "c", "duration": 47, "fits": true}],
     "chosen": "b", "remaining_after": 18, "closes": false},
    {"eligible": [{"task": "c", "duration": 47, "fits": false}],
     "chosen": null, "remaining_after": 18, "closes": true}],
   "assigned": ["a", "b"], "remaining_final": 18},
  {"station_id": 2, "capacity": 139, "remaining_initial": 139,
   "trials": [
    {"eligible": [{"task": "c", "duration": 47, "fits": true}],
     "chosen": "c", "remaining_after": 92, "closes": false},
    {"eligible": [{"task": "d", "duration": 87, "fits": true}],
     "chosen": "d", "remaining_after": 5, "closes": false},
    {"eligible": [{"task": "e", "duration": 22, "fits": false}],
     "chosen": null, "remaining_after": 5, "closes": true}],
   "assigned": ["c", "d"], "remaining_final": 5},
  {"station_id": 3, "capacity": 139, "remaining_initial": 139,
   "trials": [
    {"eligible": [{"task": "e", "duration": 22, "fits": true}],
     "chosen": "e", "remaining_after": 117, "closes": false},
    {"eligible": [], "chosen": null, "remaining_after": 117, "closes": true}],
   "assigned": ["e"], "remaining_final": 117}],
 "result": {"symbol": "n", "value": 3}}
```

Element 1, trial 2 checks the rule rather than trusting it: eligible are `b` (70)
and `c` (47), both fit in 88, `max` is `b` ✓; `88 − 70 = 18` ✓. Trial 3: `c` at
47 does not fit in 18, no candidate fits, `closes: true`, budget unchanged ✓.
`assigned = [a, b]`, `139 − 51 − 70 = 18 = remaining_final` ✓. Across elements
the union is `{a,b,c,d,e}` = `universe` ✓, disjoint ✓, `n = 3` ✓ — and `n = 3`
is the **answer**, feeding `(3·139 − 277)/(3·139)·100 = 140/417·100 = 33.573…`,
quoted half-up at one decimal as `33.6%`.

Note the last element's final trial has an **empty** `eligible`: that is
exhaustion, and §5.5 permits `closes: true` there. A verifier that requires a
non-empty `eligible` on a closing trial rejects every valid trace.

---

## 10. Limits, and what this specification does not settle

- **The node is a frame local, not a return value.** The templates' public
  contract is still `(question, solution)`, and `tests/trace_schema/extract.py`
  lifts `trace_nodes` out via `sys.settrace`. That was deliberate — changing the
  return contract is a corpus-wide decision, not a two-template one — but it
  means the node is not yet part of any interface. Promoting it to a return value
  belongs to the milestone-model work. Recorded as **D-039**.
- **Two node types are specified; two templates implement them.** The spec names
  `linear_reservoir_routing_step` (a repeated sub-chain unrolled into the trace)
  and `qr_policy_one_iteration` (iterative in principle, one iteration emitted)
  as candidates for `iteration`. Neither has been fitted, so the claim that this
  generalises is **argued, not measured**. §2's `incidental`/`answer_bearing` test
  is the thing to apply to them first.
- **§6.8 prose agreement is specified but not implemented.** In the reference
  templates the prose is rendered from the node, so the check is vacuous by
  construction; it is written down for the day a template stops doing that.
- **No dimensional checking.** Both reference nodes carry `unit` only on
  `result`. A dimensional comparator needs units on every symbol; that is a
  milestone-model decision and is deliberately out of scope here.
