# Phase 3 — Reviewer D: schema implementability

**Reviewer:** D (independent) · **Gate:** "Reviewer D implemented a working
verifier from D3.3 alone" · **Date:** 2026-09-06
**Frozen ref under review:** `953adc9` on `redesign/phase3-trace-shape`
**Inputs used:** `docs/re-implementation-sep/phase3_node_types.md` (D3.3) and
`docs/re-implementation-sep/phase3_conformance/traces.json` only.
**Artefact:** `tests/trace_schema/reviewer_d_verifier.py`

**Contamination disclosure.** My first command was `git show 953adc9 --stat`,
which printed the full commit body before I could stop it. I therefore saw the
Phase 3 commit message, which the brief excluded. I did **not** read
`phase3_summary.md`, `phase3_item_pool_impact.md`, either template source, or
`tests/trace_schema/extract.py`. The commit message named the D-038
`incidental`/`answer_bearing` asymmetry and the two corrected errors; §2 and §9
of the spec state both anyway, so I judge the leak low-impact, but the gate is
weaker than it would have been and the next reviewer should use
`git show <sha>:<path>` with no `--stat`.

---

## 1. Verdict

**PASS WITH FINDINGS** — the spec was implementable without asking the author
(80/80 gold, 33/33 negative cases), but §4 is normative only for the one
template it was written from, and three internal contradictions (F1, F2, F3)
should be amended before D3.3 is relied on as the milestone model's base.

---

## 2. Independent re-derivation

### What I built

`tests/trace_schema/reviewer_d_verifier.py`, ~700 lines, written from §3–§8 of
D3.3 with the corpus used only as input data. It implements the §6 algorithm in
the stated order: §6.1 shape, §6.2 homogeneity, §6.3 budget, §6.4 carry, §6.5
local invariants (§4 for `iteration`, §5 for `decision`), §6.6 termination,
§6.7 result, and a partial §6.8 prose check. Every point where the spec did not
determine the implementation is marked `D-A1`..`D-A7` in the source and appears
as a finding below.

Two design choices worth stating because they are what §8 is really testing:

- **§8.4 is honoured structurally.** Non-path `carry` right-hand sides are
  routed to a separate `unchecked` channel that can never contribute to a pass.
  A relation is either checked-and-passed or reported unchecked; there is no
  third state in which an unverified relation reads as satisfied.
- **The verifier is total.** §5 Determinism requires totality over traces a
  model might emit, so every invariant guards its own types and the driver traps
  `KeyError`/`TypeError`/`ValueError`/`ZeroDivisionError` into a `CRASH` failure
  rather than aborting the run.

### Results over the 80 gold traces

```
$ PYTHONIOENCODING=utf-8 python -m tests.trace_schema.reviewer_d_verifier \
      docs/re-implementation-sep/phase3_conformance/traces.json

--- s8.4 relations reported UNCHECKED (not passed) ---
  [CARRY_UNCHECKED] x40 carry['assigned'] = "union of every element's assigned"
  [CARRY_UNCHECKED] x40 carry['remaining_initial'] = 'capacity (reset each element)'
  [PRECEDENCE_UNCHECKABLE] x40 s5 global invariant / s8.6 ... UNCHECKABLE

80 passed, 0 failed, 80 nodes total
```

**80/80 pass, first run, no corrections to the verifier.** 40
`template_normal_depth_iteration` (`t24_secant_normal_depth`, element counts
`{2:4, 3:23, 4:13}`) and 40 `template_line_balancing_heuristic`
(`t19_greedy_station_assignment`, `{3:18, 4:22}`). The three `UNCHECKED`
reports are the spec working as designed for the first two and a spec gap for
the third (F1).

### Results over mutated negative cases

33 mutants, covering all seven §8 clauses plus §3, §4.1, §4.2, §5.1, §5.3,
§5.4, §5.6, §5.7, §6.1, §6.3, §6.6 and §6.7:

```
$ PYTHONIOENCODING=utf-8 python -m tests.trace_schema.reviewer_d_verifier \
      docs/re-implementation-sep/phase3_conformance/traces.json --self-test

33/33 negative cases rejected, 0 missed
=== control: unmutated gold ===
  80/80 gold nodes pass
```

Every mutant was rejected with the *expected* failure code, not merely rejected.
Coverage by §8 clause:

| §8 | Negative case | Code |
|---|---|---|
| 8.1 | `satisfied:false` on gold, both node types | `SATISFIED` |
| 8.2 | element missing a symbol; element with an extra symbol; trial missing a symbol | `HOMOGENEITY` |
| 8.3 | `y_prev<-y_curr` broken; `g_curr<-evaluation.g` broken (dotted path) | `CARRY` |
| 8.4 | *structural* — non-path RHS routed to `unchecked`, cannot pass | `CARRY_UNCHECKED` |
| 8.5 | `converged` early; `evaluation` non-null when converged; zero denominator | `ITER_CONV` / `ITER_EVAL` / `ITER_DENOM` |
| 8.6 | overlapping `assigned`; union != universe; not-longest `chosen`; non-contiguous `station_id` | `DEC_DISJOINT` / `DEC_UNIVERSE` / `DEC_RULE` / `DEC_IDS` |
| 8.6 | precedence violation | **not implementable — see F1** |
| 8.7 | `answer_bearing` without symbol; `incidental` with one | `CARD_SYMBOL` |

§8.4 is the one clause a mutation cannot test, because it is a requirement on
the verifier rather than on the trace; I discharged it by construction and by
showing both non-path entries surface in the run output above.

---

## 3. Findings

Fifteen. Numbered as the deliverable list of ambiguities. `CONFIRMED` =
reproduced against the frozen ref; all fifteen are defects in the **spec text**,
none is a defect in the corpus — no gold trace failed. F1–F3 are the ones I
consider merge-blocking; the rest are amend-before-the-milestone-model.

### F1 — §8.6 requires rejecting a precedence violation the node cannot express · CONFIRMED · merge-blocking

§5 Global invariants: *"every task's precedences are satisfied by an earlier or
same-element assignment"*, and §8.6 lists *"a precedence violated at the moment
of choice"* among the traces a conforming verifier must reject. **The node
carries no precedence relation.** `termination.universe` is a flat task list;
the only mention of predecessors is prose inside `rule`.

```
$ python -c "import json;rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'));print(sorted({k for r in rows if r['template_id']=='template_line_balancing_heuristic' for k in r['trace_nodes']}))"
```

No field names precedences. A conforming verifier is therefore *impossible* by
§8.6's own definition. Worse, the precedence filter is baked invisibly into
`eligible`: a trace that omitted a task from `eligible` when it should have been
present is indistinguishable from one that correctly withheld it on precedence
grounds. My verifier reports `PRECEDENCE_UNCHECKABLE` rather than passing it,
which is what §8.4's principle demands, but that means a documented negative
case is permanently unreachable.

**Impact.** §8's guarantee is overstated: a verifier can satisfy 6 of 7 clauses
and still be called conforming. **Should have said:** either carry the
precedence relation on the node (e.g. `termination.precedences: {"b": ["a"]}`,
alongside `universe`), or move the precedence clause out of §8 into §10's
explicit non-goals. The first is the better fix — it also makes `eligible`
independently checkable rather than trusted, which is the stated design
intent of reporting `fits` separately.

### F2 — §3 and §6.3 contradict each other on empty `elements` · CONFIRMED · merge-blocking

§3 `elements`: *"May be empty only if `termination.satisfied` is true at entry."*
§6.3: *"`1 <= len(elements) <= termination.max_elements`."* Since §3.1 fixes
`satisfied: true` on every gold node, §3 permits a case §6.3 forbids
unconditionally. I implemented §6.3 (rejecting empty), which is what the §8
posture implies, but the opposite choice is equally supported by the text.

```
$ python -m tests.trace_schema.reviewer_d_verifier docs/re-implementation-sep/phase3_conformance/traces.json --self-test 2>&1 | grep "empty elements"
  REJECT  s6.3   empty elements on a satisfied node  [BUDGET]
```

**Impact.** Two conforming verifiers disagree on a real trace class. Not
exercised by gold (no node has zero elements), so it costs nothing today and
everything the first time a model emits one. **Should have said:** pick one.
If the zero-element case is meant to be reachable, §6.3's lower bound must read
`0 <= len(elements)` with the §3 condition attached; otherwise delete the §3
sentence.

### F3 — §4 has no normative element table, so `iteration` is specified only for its one reference template · CONFIRMED · merge-blocking

This is the largest finding. §5 gives `decision` a **normative element table and
trial table** fixing `station_id`/`capacity`/`trials`/`assigned`/… as type-level
names. §4 gives `iteration` **no such table**. Its invariants instead read:

> 2. `change == |y_next − y_curr|`, at the display precision of `y`.
> 3. `converged == (change < termination.tolerance)`.
> 5. `evaluation` is `null` iff `converged` is `true`.

`y_next`, `y_curr`, `change`, `converged`, `evaluation` are used as if they were
type-level symbol names, but §3 says `element_symbols` is *"the template's own
numeric or symbolic domain"*. And §4.1 states the update relation only *"for the
reference template"* — the secant step is nowhere encoded in the node.

The consequence: I could implement §6.5 for `iteration` **only by hardcoding
`uniform_flow`'s symbol names and its secant formula** (`D-A2`, `D-A3` in the
source). Applied to `linear_reservoir_routing_step` or `qr_policy_one_iteration`
— the two templates §10 names as the next `iteration` candidates — my verifier
would fail every trace at `ITER_SHAPE`, not because the traces are wrong but
because their iterate is not called `y`.

```
$ grep -n "for the reference template that relation" docs/re-implementation-sep/phase3_node_types.md
$ grep -n "Element\." docs/re-implementation-sep/phase3_node_types.md   # s4 has prose, s5 has a table
```

**Impact.** §10 already concedes that generalisation is *"argued, not
measured"*, but the gap is sharper than that: `decision` generalises as a
*schema* and `iteration` does not, because only one of the two was written as a
type. §1 claims the document *"is the deliverable the milestone model is built
on"*; built on this, the milestone model gets one reusable node type and one
template-specific one. **Should have said:** give §4 an element table with the
same status as §5's, naming the roles (`iterate_prev`, `iterate_curr`,
`iterate_next`, `residual_prev`, `residual_curr`, `change`, `converged`,
`evaluation`) and add a node field binding roles to the template's own symbols
— e.g. `"roles": {"iterate_curr": "y_curr", "residual_curr": "g_curr"}` — plus a
declared `update_relation` so §4.1 is data, not prose.

### F4 — per-symbol display precision is used by three clauses and declared nowhere · CONFIRMED

§4.2 says `change` matches *"at the display precision of `y`"*; §7.1 says a
candidate matches *"to half a unit in the last place gold displays it
(`result.dp`, **or the element symbol's display precision**)"*. **No node field
carries a per-symbol precision** — `dp` exists only on `result`, and on the
`iteration` node it is `3` while `y_next` is displayed at `4`. So `result.dp`
is actively the *wrong* number for the element-level checks that cite it.

I inferred the map from the data (`y*` 4 dp, `g*` 3 dp) and hardcoded it as
`DISPLAY_DP` (`D-A4`):

```
$ python -c "
import json;rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))
f=lambda x:len(repr(float(x)).split('.')[1].rstrip('0'))
print('max dp y_next',max(f(e['y_next']) for r in rows if r['template_id']=='template_normal_depth_iteration' for e in r['trace_nodes']['elements']))
print('result.dp',[r['trace_nodes']['result']['dp'] for r in rows if r['template_id']=='template_normal_depth_iteration'][0])"
max dp y_next 4
result.dp 3
```

**Impact.** Inferring precision from observed decimals is unsound in general — a
value that happens to land on 1.3800 reads as 2 dp. The next implementer will
infer a different map and disagree at the boundary. **Should have said:** add
`element_dp: {"y_prev": 4, ..., "g_curr": 3}` (or fold `dp` into a per-symbol
descriptor) and have §7.1 cite that field instead of a parenthetical.

### F5 — `termination.predicate` is required, normative in §6.6, and not machine-evaluable · CONFIRMED

§3.1 requires `predicate`, described as *"in the template's own notation,
evaluable over element symbols"*. §6.6 requires *"the predicate holds on the
last element and on no earlier one"*. But the two corpus predicates are
`"abs(y_next - y_curr) < tol"` and `"len(assigned) == len(tasks)"`, which
contain the free names `tol` and `tasks` with **no binding environment
specified** and no grammar stated. `tol` is presumably `termination.tolerance`
and `tasks` presumably `termination.universe`, but neither mapping is written
down.

I did not evaluate the strings at all; I re-expressed §6.6 through `change` vs
`tolerance` and through the cumulative assignment count (`D-A6`). That works
only because §4.3 happens to restate the convergence predicate in checkable
terms.

**Impact.** "Evaluable over element symbols" is false as written — the
predicates are evaluable over element symbols *plus two unnamed bindings*.
Any verifier that took the field at its word and `eval`'d it would raise
`NameError` on all 80 traces. **Should have said:** either declare the
binding environment (`tol -> termination.tolerance`, `tasks ->
termination.universe`) and fix a grammar, or demote `predicate` to advisory
prose the way §3.2 demotes non-path `carry` entries — it is currently the only
required field that is normative in §6 and unusable.

### F6 — §3.2 gives no syntax separating a `carry` path from a `carry` prose description · CONFIRMED

§3.2: *"The right-hand side may also be a prose description when the relation is
not a projection"*, and §8.4 makes silently passing one a rejectable verifier
defect. But nothing says how to tell them apart. The corpus has
`"evaluation.g"` (path) and `"capacity (reset each element)"` (prose) — and
note that **`"capacity"` alone would be a perfectly valid path**, so the
discriminator cannot be "does it resolve".

I used a dotted-identifier regex (`D-A5`), which gets the corpus right:

```
$ python -m tests.trace_schema.reviewer_d_verifier docs/re-implementation-sep/phase3_conformance/traces.json 2>&1 | grep UNCHECKED
  [CARRY_UNCHECKED] x40 carry['assigned'] = "union of every element's assigned" ...
  [CARRY_UNCHECKED] x40 carry['remaining_initial'] = 'capacity (reset each element)' ...
```

**Impact.** The failure mode is silent and is exactly the one §8.4 exists to
prevent: a stricter regex mis-reads a prose RHS as a path, fails to resolve it,
and reports a spurious `CARRY` failure; a looser one mis-reads a prose RHS as a
path that *does* resolve and passes an unchecked relation as satisfied.
**Should have said:** tag the value —
`{"assigned": {"kind": "prose", "text": "union of ..."}}` vs
`{"y_prev": {"kind": "path", "path": "y_curr"}}` — or reserve a sigil. A
type distinction this load-bearing should not rest on whether a string contains
a space.

### F7 — `termination.quantity: "assigned"` collides with the element symbol `assigned` · CONFIRMED

For `decision`, `quantity` is `"assigned"` and the predicate is
`len(assigned) == len(tasks)`. But `assigned` is *also* a declared
`element_symbols` entry. Read at element level the predicate can never hold —
the largest per-element `assigned` in the corpus is 3 against a universe of 5:

```
$ python -c "
import json;rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))
lb=[r for r in rows if r['template_id']=='template_line_balancing_heuristic']
print('max per-element len(assigned)',max(len(e['assigned']) for r in lb for e in r['trace_nodes']['elements']),'universe',5)"
max per-element len(assigned) 3 universe 5
```

so §6.6's *"the predicate holds on the last element"* is unsatisfiable under the
element-level reading. Only the cumulative-union reading works, and §3.2's
`carry` entry (`"union of every element's assigned"`) is the sole hint.

**Impact.** I guessed correctly, but a verifier that guessed the other way
rejects all 40 `decision` traces. **Should have said:** name the accumulator
distinctly (`quantity: "assigned_union"`) and state in §3.1 that `quantity`
ranges over node-level accumulators, not `element_symbols`.

### F8 — the container shape of `trace_nodes` is never specified · CONFIRMED

§9 says the corpus carries seeds *"as `{seed, question, solution,
trace_nodes}`"* and nothing more. The field is **plural** and holds a **single
node object**. A reader must guess between one node, a list of nodes, and a
`node_id -> node` map — and the plural name argues against the truth. I
accepted all three defensively (`D-A1`, `nodes_of`).

```
$ python -c "import json;print(type(json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))[0]['trace_nodes']))"
<class 'dict'>
```

**Impact.** Small but immediate — it is the first thing any implementer hits.
It also matters for §10/D-039: when the node is promoted to a return value, a
template with two variable-shape steps will need the plural for real, and the
shape will change under anyone who assumed the singular. **Should have said:**
state the container explicitly and, preferably, make it a list now so the
promotion is not a breaking change.

### F9 — §6.7's `round()` has no stated rounding mode · PLAUSIBLE

§6.7: *"`result.value == round(last.y_next, result.dp)`"*. Python's `round` is
banker's rounding; §9.2 quotes the balance delay *"half-up at one decimal"*.
The two conventions are in the same document and disagree at exactly the
half-way point.

§7.1 says half-way instances are screened out at generation (2.30% / 0.65%
rejection), so this cannot bite a gold trace — and indeed all 80 pass under
either. It bites a *candidate* trace, which §7 explicitly scopes in and which
is not screened.

```
$ python -m tests.trace_schema.reviewer_d_verifier docs/re-implementation-sep/phase3_conformance/traces.json --self-test 2>&1 | tail -3
```

**Impact.** Latent; scoring-level, not gold-level. **Should have said:**
"half-up" everywhere, once, in §7.1, since the generator already commits to it.

### F10 — §6.8 is a numbered step of the normative algorithm that cannot be implemented · CONFIRMED

§6.8 requires *"every numeric token the node carries appears in the prose at its
stated precision, and no arithmetic line in the prose is absent from the node"*.
Neither half is implementable: the first needs the per-symbol precision F4 says
does not exist, and the second needs a definition of "arithmetic line" that the
document never gives. §10 concedes it is *"specified but not implemented"* —
but it is still step 8 of a numbered normative algorithm, so a verifier claiming
conformance to §6 cannot actually reach conformance.

I implemented the one unambiguous fragment (result value appears in the printed
solution); all 80 pass it.

**Impact.** §6 as a whole is not fully implementable, which weakens the
document's central claim in §1. **Should have said:** move §6.8 out of the
numbered algorithm into §10, or mark it `OPTIONAL` inline so the algorithm's
conformance surface is exactly the implementable part.

### F11 — `closes` is specified twice, with two different rules · CONFIRMED

§5 Trial table: *"`closes` — `true` **exactly on the trial where `chosen` is
`null`**"*. §5.5: *"`closes` is `true` on the **last trial** and false on every
other, **except** on the final element where the universe is exhausted: there
the last trial may have an empty `eligible` and `closes: true`."*

Three problems. (a) The two definitions coincide only if every element's last
trial has `chosen: null`, which is true in the corpus but stated nowhere. (b)
The §5.5 "except" clause's scope is genuinely unclear — read literally it
exempts the *final element* from "closes true on the last trial", which is the
opposite of what the following words say; the intended exception is about
`eligible` being empty, not about `closes`. (c) "may have an empty `eligible`"
is permissive, but §9.2's closing note says *"a verifier that requires a
non-empty `eligible` on a closing trial rejects every valid trace"* — so the
empty case is not merely permitted, it is universal on the final element:

```
$ python -c "
import json;rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))
lb=[r for r in rows if r['template_id']=='template_line_balancing_heuristic']
print('final stations with NON-empty eligible on the closing trial:',sum(1 for r in lb if r['trace_nodes']['elements'][-1]['trials'][-1]['eligible']))
print('non-final stations with empty eligible on the closing trial:',sum(1 for r in lb for i,e in enumerate(r['trace_nodes']['elements']) if not e['trials'][-1]['eligible'] and i!=len(r['trace_nodes']['elements'])-1))"
final stations with NON-empty eligible on the closing trial: 0
non-final stations with empty eligible on the closing trial: 0
```

**Impact.** I implemented both readings jointly (`closes` iff `chosen is null`,
*and* last-trial-only, *and* empty `eligible` allowed only on the final element)
and all 80 pass, so the corpus does not discriminate. A verifier picking one
reading may accept traces another rejects. **Should have said:** state it once,
as `closes == (chosen is null)`, and note separately that this necessarily
falls on the last trial of each element.

### F12 — the `eligible` entry schema is prose-only, with no `*_symbols` declaration · CONFIRMED

The node declares `element_symbols` and `trial_symbols`, and §6.2 requires key
sets to match *"Same for `preamble`/`preamble_symbols` and
`trials`/`trial_symbols`"*. But `eligible` entries are a **third** nesting
level, `{task, duration, fits}`, given only in §5's prose table with no
corresponding `*_symbols` field — so §6.2's homogeneity check has nothing to
compare against at that level and I hardcoded the key set (`ELIGIBLE_KEYS`).

**Impact.** The homogeneity guarantee — §1's stated reason the flat milestone
model failed, and §8.2's rejectable case — stops one level short of the
deepest structure in the schema. **Should have said:** add
`eligible_symbols: ["task", "duration", "fits"]` and extend §6.2 to it, or
state that homogeneity is declared only to depth 2 and why.

### F13 — "`node_id` must not encode the element count" is normative and untestable · PLAUSIBLE

§3 requires it but gives no test. I implemented a heuristic (reject if the
element count appears as a standalone number in `node_id`), which rejects
`"t19_greedy_4_stations"` on a 4-element node but would also reject a legitimate
`"t04_..."` id on a 4-element node — a false positive I cannot rule out without
seeing more node ids.

```
$ python -m tests.trace_schema.reviewer_d_verifier docs/re-implementation-sep/phase3_conformance/traces.json --self-test 2>&1 | grep node_id
  REJECT  s3     node_id encodes the element count  [NODE_ID]
```

**Impact.** Low today (both corpus ids are clean) but the check is unsound and
another implementer would reasonably not implement it at all. **Should have
said:** make it a convention in §10 rather than a `✔`-required constraint in
§3, or give it a decidable form.

### F14 — five invariants the spec leaves open, each admitting a corrupted trace · CONFIRMED

A verifier faithful to §3–§8 accepts all of the following. Reproduce with the
probe script quoted in §4 below:

| # | Mutation | Why it passes |
|---|---|---|
| H2 | `k` set to `[7, 3, 99]` | §5 requires `station_id` be `1..n` contiguous; **§4 states no analogue for `k`**, and `k` is not in `carry` |
| H3 | every preamble `A`, `P`, `AR` set to `-1.0` | §4 Preamble constrains only that the trace follow the question's *starting depths*; the derived geometry is unconstrained |
| H4 | `evaluation.y` set to `42.0` | `carry` binds `g_curr <- evaluation.g` but nothing binds `evaluation.y` to `y_next`, though they are the same quantity |
| H5 | `max_elements` set to `500` on a `decision` node | §3.1 calls it "the budget"; nothing relates it to `len(universe)`, which is its only sensible value for `exhaustion` |
| H1b | final station's `capacity` set to `1147` while the others are `147` | `capacity` is per-element and **never required constant**, though it is the single cycle time `CT` that §2 puts in the answer |

H1b is the serious one. §2's whole argument is that `n` is answer-bearing
because the delay is `(n·CT − Σt)/(n·CT)·100` — a single `CT`. A trace in which
`CT` varies by station is arithmetically self-consistent under every §5
invariant and passes, yet has no meaning. All 40 corpus traces have constant
capacity, so nothing detects the drift.

**Impact.** §8's list is a floor, not a characterisation, and reads as though it
were the latter (*"a verifier that accepts everything is not a verifier"*).
**Should have said:** add `capacity` constancy and `max_elements ==
len(universe)` for `exhaustion` to §5's global invariants; add `k` contiguity to
§4; bind `evaluation.y` in `carry` (it is a pure projection, so it costs
nothing).

### F15 — the §4.1 × §7.1 tolerance has 0.3% margin on the committed corpus · CONFIRMED

§4.1 requires the update relation hold *"to the tolerance in §7.1"*, i.e. half a
unit in `y_next`'s last displayed place = 5.0e-05. The worst residual in the
corpus is **4.985e-05** — 99.7% of the budget.

```
$ python -c "
import json;rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))
w=0
for r in rows:
    if r['template_id']!='template_normal_depth_iteration': continue
    for e in r['trace_nodes']['elements']:
        p=e['y_curr']-e['g_curr']*(e['y_curr']-e['y_prev'])/(e['g_curr']-e['g_prev'])
        w=max(w,abs(p-e['y_next']))
print('worst residual %.4e vs tolerance 5.0e-05'%w)"
worst residual 4.985e-05 vs tolerance 5.0e-05
```

This is expected — the secant step is computed from `g` values already rounded
to 3 dp, so the residual is dominated by that input rounding, not by the
relation. But the spec never says so, and never says the check is
tolerance-critical.

**Impact.** Any implementer who reads "half a unit in the last place" as a
*strict* `<` with no epsilon, or who reasons that a rounded relation deserves a
tighter bound, lands within 3e-7 of failing gold traces. That is a
reproducibility hazard for exactly the audience §1 addresses. Reassuringly,
three defensible readings agree on the corpus today — `|d| <= 5e-5`,
`|d| < 5e-5`, and `round(secant, 4) == y_next` all give 129/129 — so the
ambiguity is latent rather than live. **Should have said:** state in §4.1 that
the residual is bounded by the 3 dp rounding of `g`, and give the intended
comparison explicitly rather than by reference to §7.1.

---

## 4. Falsification attempts that failed

Things I expected to break and could not.

1. **All 80 gold traces, first run, no verifier corrections.** I expected the
   spec to be wrong about at least one detail of its own corpus. It was not.
   §4's five local invariants, §5's seven local plus three global, §6.4's carry
   projections including the dotted `evaluation.g`, §6.6 and §6.7 all hold
   exactly as written on every trace.

2. **The two worked examples in §9 are faithful.** I re-derived §9.1 element 2
   independently (`1.3562 − (−0.246)(1.3562 − 1.5)/(−0.246 − 1.33) = 1.378646`
   → `1.3786`) and §9.2's `139 − 51 − 70 = 18`, `(3·139 − 277)/417·100 =
   33.573`. Both correct, and both match rows actually present in the corpus —
   §9's "lifted verbatim" claim holds.

3. **I tried to make the tie-break matter and could not.** §5.3 specifies a
   tie-break by `_T19_ORDER` position; §5 Determinism claims 0 duplicate-duration
   instances in 4,000 seeds. I searched the 40 committed traces for any trial
   with two fitting candidates of equal duration: **0 found**. The tie-break is
   genuinely unexercised, exactly as claimed, and specifying it anyway is the
   right call for the candidate-trace case.

4. **`carry` really is sufficient to check element *k* against *k−1* in
   isolation**, as §3.2 claims. I verified all four `iteration` projections
   pairwise across every trace with no reference to the count, and mutating
   either a direct (`y_prev`) or a dotted (`evaluation.g`) path is caught. The
   central design claim of the document survives.

5. **`station_id` contiguity, `assigned` disjointness, and the universe union
   hold on all 40** `decision` traces; I could not find a trace where the
   greedy replay disagreed with the recorded `chosen`, across all 143 trials.
   §5.3's "the rule is checked, not assumed" is real and it passes.

6. **I could not find a `satisfied: false` gold trace**, nor a trace exceeding
   `max_elements` (max observed 4 against a budget of 5 for both templates),
   nor a non-boolean `converged`, nor an element whose key set drifted from
   `element_symbols`. §1's motivating defect — symbols appearing and
   disappearing between passes — is genuinely absent.

7. **The `incidental`/`answer_bearing` asymmetry is not merely asserted; it is
   load-bearing in code.** Implementing §6.7 forced two different result rules
   (`round(last.y_next, dp)` vs `len(elements)`) and §8.7 forced two different
   shape rules. A merged node type would have had to pick one. §2's central
   claim survived the attempt to implement around it.

The probe script for F14, for reproduction:

```
$ PYTHONIOENCODING=utf-8 python -c "
import json,copy
from tests.trace_schema.reviewer_d_verifier import verify_node
rows=json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json',encoding='utf-8'))
def get(t,p=lambda n:True):
    for r in rows:
        if r['template_id']==t and p(r['trace_nodes']): return copy.deepcopy(r['trace_nodes'])
n=get('template_normal_depth_iteration', lambda x: len(x['elements'])>=3)
for e,v in zip(n['elements'],[7,3,99]): e['k']=v
print('H2 k=[7,3,99] ->','ACCEPTED (hole)' if verify_node(n,'H2').ok else 'rejected')"
H2 k=[7,3,99] -> ACCEPTED (hole)
```

---

## 5. Further probing and improvements

**Where the evidence is thinnest.** The document is validated against exactly
two templates, one per node type, and the `iteration` type is not separable from
its single instance (F3). Every claim about generalisation rests on an argument,
which §10 admits. My verifier inherits that: it is two template-specific
checkers behind one dispatcher, and only the `decision` half would survive
being pointed at a new template.

**What I would probe next, ranked.**

1. **Fit a third template to `iteration`** — `linear_reservoir_routing_step` or
   `qr_policy_one_iteration` per §10 — and try to verify it with an *unmodified*
   verifier. I predict it fails immediately on symbol names, which would
   confirm F3 empirically rather than by argument. This is the single highest-value
   next experiment and it is cheap.
2. **Run the comparator half.** I implemented §6 (gold verification) but §7
   (D3.4 candidate scoring) is untested by anything — no candidate traces exist.
   §7.2's six dispositions and §7.3's four are entirely unexercised, including
   the `cardinality_mismatch` report that §2's whole argument exists to justify.
   Constructing ten synthetic candidate traces per template (wrong count, right
   answer; right count, different packing; early stop) would test the part of
   the document that actually decides scores.
3. **Fuzz the schema.** My 33 mutants were hand-written from §8; a generative
   fuzzer perturbing one field at a time across all 80 traces would find more
   F14-class holes systematically. F14 came from five guesses and all five hit,
   which suggests the density of open invariants is high.
4. **Check the corpus against a second extraction.** Everything I verified came
   from one `traces.json` produced by one extractor. §6.8's prose-agreement
   check exists precisely to catch node/prose divergence and is unimplemented
   (F10); re-deriving a handful of nodes from the printed solutions by hand
   would be the only independent check that the corpus says what the templates
   print.

**What generalises.** The `decision` type does — its element and trial tables
are role-named and template-neutral, and I would expect any greedy
bin-packing or sequencing template to fit it with only `universe` and `rule`
changing. The `carry`-as-projection idea generalises further still and is the
document's best contribution: it is what makes variable-length traces checkable
without replay, and it worked exactly as advertised. The `incidental` /
`answer_bearing` test in §2 is stated well enough to apply to an unseen
template, which was the stated goal and is met.

**What later phases should do differently.** Two process points. First, the
asymmetry between §4 and §5 suggests the two types were written at different
times or by different reasoning; a spec that fixes one type's element schema in
a table should fix the other's the same way, and a checklist of "does every node
type have an element table, a role binding, and a declared update relation"
would have caught F3, F4 and F12 at authoring time. Second, §8 should be
generated from the verifier rather than written alongside it — F14's five holes
exist because §8 was written as a list of remembered failure modes rather than
as a characterisation of the invariant set, and a list written that way is
always a floor.

**What was wrong or ambiguous.** F1 (precedence unverifiable), F2 (empty-elements
contradiction) and F3 (no `iteration` element table) are the three I would fix
before merging; all three are edits to the document, none touches the corpus.
F5, F7 and F10 make required or numbered parts of the spec unusable as written
and should follow. F4, F6, F8, F11, F12 and F14 are ambiguities I resolved by
guessing and happened to guess right — which is the least reassuring category,
because the corpus does not discriminate and a second implementer has no signal
that they guessed differently.

**Gate assessment.** The gate asks whether the spec is implementable by someone
who did not design it. It is: 80/80 and 33/33, in one sitting, without opening a
template. That is a genuine pass and the document deserves credit for it —
§6's ordered algorithm and §8's explicit reject list are why it was possible,
and most specs I could have been handed would not have survived this exercise.
The qualification is that "implementable" here means implementable *for these
two templates*; §4 does not yet specify a type, and the milestone model that §1
says will be built on this document will discover that at the third template.
