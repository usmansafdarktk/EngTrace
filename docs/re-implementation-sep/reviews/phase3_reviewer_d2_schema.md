# Phase 3 — Reviewer D2: schema implementability, round 2

**Reviewer:** D2 (independent) · **Role:** implement a verifier from D3.3 alone
**Ref under review:** `9398ccc` on `redesign/phase3-trace-shape`, schema version `1.1`
**Date:** 2026-09-06 · **Time box:** 60 minutes, filed within it

**Inputs I read.** `docs/re-implementation-sep/phase3_node_types.md` and
`docs/re-implementation-sep/phase3_conformance/traces.json`, both at `9398ccc`.
Nothing else. I did **not** open `uniform_flow.py`, `production_planning.py`,
`tests/trace_schema/extract.py`, `phase3_summary.md`,
`phase3_item_pool_impact.md`, `DECISIONS.md`, any Phase 3 commit message, round
1's verifier, or round 1's report. I never needed the template source: **the
"if you need them, that is the finding" escape hatch was not used.** The spec's
own Revision note told me a round 1 existed and named its finding numbers; that
is all I know of it.

**Drift check — the branch DID move under me, and it does not matter.** At the
start of the review `HEAD` and `redesign/phase3-trace-shape` were both
`9398ccc07546773c16f476c8a7fbacacf2068657`, with `phase3_summary.md` modified and
three untracked round-1 artefacts in the tree. By the end both were
`7d9fd01a394e1f4482b3fe69b65713b31ebda477` — a fast-forward (`git merge-base
--is-ancestor 9398ccc 7d9fd01` → yes) that landed those artefacts. Per the
commission I ignored the movement and reviewed `9398ccc`, and I confirmed the
movement is immaterial rather than assuming it:

- `docs/re-implementation-sep/phase3_node_types.md` — SHA-256 identical at
  `9398ccc` and `7d9fd01` (`98d54988a382b1f0…`).
- `docs/re-implementation-sep/phase3_conformance/traces.json` — SHA-256 identical
  at both (`02b0567dc954dbb2…`).
- `git diff --name-status 9398ccc 7d9fd01` (names only; **no commit message was
  printed, and I ran no `git log` or bare `git show` on a Phase 3 commit**) lists
  five paths, none of which is either of my two inputs. I opened none of them.

**Neither of the two files under review changed. Every number below stands at
both SHAs.**

One further check, because it would otherwise have invalidated the whole review:
the worktree copy of `traces.json` is **not** byte-identical to the ref copy. The
difference is entirely CRLF-vs-LF from the Windows checkout — `git status` reports
the path clean, and parsing both gives 80 rows with every `trace_nodes` object
equal (`[r['trace_nodes'] for r in a] == [... for r in b]` → `True`). The verifier
runs below therefore reflect the corpus at `9398ccc`.

---

## 1 Verdict

**PASS WITH FINDINGS** — with two CONFIRMED findings that block the merge.

The gate asks whether the spec can be implemented by someone who did not design
it, without reaching outside it. **It can, and I did.** In roughly forty minutes
I wrote a verifier from §6 that parses `update_relation` and `carry` from their
declared grammars, resolves every `iteration` quantity through `roles`, and
replays `selection` and `precedences` from data. It passes 80/80 gold traces,
rejects 32/32 constructed negative cases across all nine §8 clauses, and accepts
a hand-built `iteration` node from a different template with different symbol
names and a different recurrence **with no code change**. On its central claim —
that §4 is now a *type* — 1.1 delivers.

Two defects nevertheless block:

- **F1 (CONFIRMED).** §3.1 marks `roles` required for *every* node and §6.1
  requires required fields present, but **all 40 gold `decision` nodes carry no
  `roles` key** — as does §9.2's own normative worked example. A verifier that
  implements §6.1 as written fails half the conformance corpus. Mine does.
- **F3 (CONFIRMED).** The normative algorithm **accepts a `decision` trace whose
  answer is wrong.** `selection.filter` is trusted, never checked, and §5.4.4
  optimises only over the filtered set, so lying about `fits` changes `n` — which
  §2 and §7.3 make *the answer* — while satisfying every clause of §5.3, §5.4,
  §6 and §8.

F2 is CONFIRMED and is the answer to the "type or description" question: `iteration`
is a type, `decision` is still a description.

---

## 2 Independent re-derivation

### 2.1 What I built

`tests/trace_schema/reviewer_d2_verifier.py` — ~700 lines, standard library only.

- **§6.1** version/shape gate; per-type required and conditional field lists;
  `cardinality` and `termination.kind` checked against the values §4/§5 fix.
- **§6.2** role binding: mandatory roles present, every `roles` and
  `evaluation_roles` value a member of the corresponding `*_symbols`.
- **§6.3** homogeneity against all five declared key sets (`element_symbols`,
  `preamble_symbols`, `evaluation_symbols`, `trial_symbols`, `eligible_symbols`),
  in both directions.
- **§6.4** `1 <= len(elements) <= max_elements`.
- **§6.5** `carry`: a recursive-descent path parser implementing §3.5's grammar
  literally (identifiers joined by `.`, root in `element_symbols`, nothing else),
  evaluated between every consecutive pair. `carry_notes` entries go to a separate
  `unchecked` channel that flips the verdict string to `PASS+UNCHECKED`, so a node
  carrying one can never report a bare `PASS` (§8.4).
- **§6.6** §4.3's eight invariants and §5.3/§5.4's eight, all resolved through
  `roles` for `iteration`.
- **§6.7** termination, both kinds, with the predicate required on the last
  element and on no earlier one.
- **§6.8** `result` dispatched on `result.from`, never on `node_type`, rounded
  with `decimal.ROUND_HALF_UP` because `rounding` says `decimal-half-up`.
- **§4.2** a recursive-descent expression parser (role identifiers, decimal
  literals, `+ - * /`, unary `-`, parentheses). **No `eval`.** A name that is not
  a role of the node is a parse-time rejection; a zero divisor raises §4.3.7.
- **§7.1** tolerance is `0.5 * 10**-symbol_precision[symbol]`, resolved through
  `roles`, and is **never** inferred from the data — a symbol missing from
  `symbol_precision` is a hard failure rather than a guess.

`tests/trace_schema/reviewer_d2_negatives.py` holds the §8 mutation harness and
the reusability probe.

### 2.2 Results over the 80 gold traces

```
python -m tests.trace_schema.reviewer_d2_verifier docs/re-implementation-sep/phase3_conformance/traces.json
```

```
FAIL template_line_balancing_heuristic seed 0
      [6.1] required field 'roles' absent
template_line_balancing_heuristic        pass  0  fail 40  clauses 6.1
template_normal_depth_iteration          pass 40  fail  0
TOTAL 80 rows: 80 pass, 40 fail; 0 rows carried unchecked carry_notes
```

All 40 `iteration` traces pass a strict reading of §6 on the first run. All 40
`decision` traces fail on **one** clause and one only — the missing `roles` key
of F1. With that single field waived:

```
python -m tests.trace_schema.reviewer_d2_verifier docs/re-implementation-sep/phase3_conformance/traces.json --lenient-roles
```

```
template_line_balancing_heuristic        pass 40  fail  0
template_normal_depth_iteration          pass 40  fail  0
TOTAL 80 rows: 80 pass, 0 fail; 40 rows carried unchecked carry_notes
```

**80/80, zero failures, one documented deviation.** The 40 `PASS+UNCHECKED` rows
are the `decision` nodes' two `carry_notes` entries, correctly quarantined per §8.4.

### 2.3 Negative cases — all nine §8 clauses

```
python -m tests.trace_schema.reviewer_d2_negatives docs/re-implementation-sep/phase3_conformance/traces.json
```

**32 constructed mutations, 32 rejected with the expected clause, 0 missed.**

| §8 | Mutations | Result |
|---|---|---|
| 8.1 | `satisfied:false` on a gold `iteration`; on a gold `decision` | 2/2 |
| 8.2 | element extra key; element missing key; evaluation frame missing key; preamble frame extra key; trial extra key; `eligible` entry missing key | 6/6 |
| 8.3 | carry relation broken between elements; path with an index `y_curr[0]`; path with a space `evaluation .g`; path root not in `element_symbols` | 4/4 |
| 8.4 | *not a trace mutation* — see below | 1/1 |
| 8.5 | `converged:true` before the last element; `evaluation` non-null on the terminating element; frame bound to the wrong iterate; zero denominator (`g_prev := g_curr`); non-contiguous `index` | 5/5 |
| 8.6 | overlapping `assigned`; union ≠ `universe`; `chosen` not the optimum; `eligible` omitting a precedence-satisfied item; `eligible` including a precedence-unmet item; non-contiguous `station_id`; capacity varied under `capacity_constant` | 7/7 |
| 8.7 | `answer_bearing` without `cardinality_symbol`; `incidental` with one | 2/2 |
| 8.8 | mandatory role absent; `roles` value not in `element_symbols`; `evaluation_roles` value not in `evaluation_symbols` | 3/3 |
| 8.9 | `schema_version: "2.0"` | 1/1 |
| §6.4 | zero-element sequence (not in §8; added) | 1/1 |

**§8.4 is the one clause that is not constructible as a trace mutation**, and this
is a property of the clause, not a gap in my testing. "Silently passing a
`carry_notes` entry" is a statement about *verifier behaviour*, not about a
trace: no trace is malformed by having `carry_notes`, and the 40 gold `decision`
nodes all carry two. I tested it as a behavioural assertion instead — a node with
a `carry_notes` entry must not be reportable as a bare `PASS`. My verifier returns
`PASS+UNCHECKED` with exactly one quarantined entry. Worth stating plainly: **an
implementer who reads §8 as a list of traces to reject will silently skip 8.4**,
because it is the only entry in the list that is not one.

Two of the §8.6 mutations are only reachable **because** 1.1 added
`termination.precedences`. I constructed them by deleting a legitimately eligible
entry and by appending an entry (`e`, whose predecessor `d` is uncommitted) to the
first trial; both were caught at §5.4.3b, in both directions. The clause the
Revision note says was unreachable in 1.0 is now genuinely reachable, and I reached it.

### 2.4 The different-symbol-names test — the headline result

I built, from the spec alone, an `iteration` node for a **synthetic template that
does not exist in this repository**: solving `S² − 2 = 0` by a **damped** secant.
Its iterate is `S_curr`/`S_prev`/`S_next`, its residual `r_curr`/`r_prev`, its
index `n`, its change `delta`, its converged flag `done`, its evaluation frame
`frame`, and its evaluation symbols are `["S", "r"]` — not one symbol name in
common with the reference template except by coincidence. Its recurrence is
**different**, carrying a literal the reference relation does not have:

```
iterate_curr - 0.5 * residual_curr * (iterate_curr - iterate_prev) / (residual_curr - residual_prev)
```

It converges in **9** elements, against the reference template's 2–4.

```
REUSABILITY PROBE -- an `iteration` node from a different template
iterate symbol: S_curr   residual symbol: r_curr   index symbol: n
elements: 9   result: {'symbol': 'Sroot', 'value': 1.415, 'dp': 3, 'from': 'last.iterate_next'}
verdict: PASS
corrupted alien node verdict: FAIL  ['8.3', '4.3.1', '4.3.2', '4.3.6/8.5']
```

**The unmodified verifier accepts it, and rejects it the moment one `S_next` is
perturbed by 0.01.** No edit, no branch on `node_id`, no literal symbol name.

**§3.3's and §4.2's claims are true for `iteration`, and I verified them
adversarially rather than taking them on trust.** Grepping my verifier for the
reference template's symbols (`y_curr`, `g_curr`, `y_next`, `y_prev`, `g_prev`,
`change`, `converged`, `k`) returns nothing outside comments and docstrings. The
only literal symbol names in the file are the `decision` ones of F2, and they are
there because §5 gives an implementer nothing else — which is F2.

### 2.5 A number of the spec's own, reproduced

§4.5 claims the worst §4.3.1 residual over the committed corpus is `4.985e-05`
against a `5.0e-05` tolerance. Recomputing it independently from the corpus:
**`4.9850e-05`**. It reproduces exactly. §4.5 earns its place: I had briefly
considered tightening the tolerance and it would have failed gold traces.

---

## 3 Findings

### F1 — `roles` is required of every node; no `decision` node has one · **CONFIRMED** · blocks merge

§3.1 marks `roles` with a hard ✔ ("Req.") and describes it as required of both
types; §6.1 step 1 requires "required fields present"; §8.8 makes a `roles`
defect a mandatory rejection. But all 40 gold `decision` nodes omit `roles`
entirely, and so does §9.2 — the worked example the spec says was "lifted
**verbatim, by script**" from the corpus and "machine-compared against it".

The comparison that was added in 1.1 to stop the 1.0 defect of a hand-transcribed
example omitting a required field did not catch the same defect one level up: it
compares the *example* to the *corpus*, not the *corpus* to §3.1's field table.

Reproduce:
```
python -m tests.trace_schema.reviewer_d2_verifier docs/re-implementation-sep/phase3_conformance/traces.json
```
→ `template_line_balancing_heuristic  pass 0  fail 40  clauses 6.1`, every failure
`required field 'roles' absent`.

**Impact.** Any independent implementer who follows §3.1 fails 40/80 gold traces
on their first run, exactly as I did. That is the gate's own scenario, and it is
the difference between the exit gate reading "implemented a working verifier" and
"implemented a verifier that needed a waiver". The fix is one line either way:
make `roles` conditional on `iteration`, or resolve F2 and emit `roles` on
`decision` nodes. **The second is the right fix** — see F2.

### F2 — `decision` is still a description, not a type · **CONFIRMED** · blocks merge

§3.3 states the governing claim: `roles` "is what makes the type reusable", and
"a verifier that reaches for a literal symbol name is not conforming". §4 honours
it. §5 does not. §5.1's element and trial tables are headed **"Symbol"**, not
"Role", and name `station_id`, `capacity`, `remaining_initial`, `trials`,
`assigned`, `remaining_final`, `eligible`, `chosen`, `remaining_after` and
`closes` as literal strings. There is no role table for `decision` anywhere in
the document.

The contradiction is visible inside §8 itself. §8.8 says to reject "a role named
in §4.1/**§5.1** as mandatory but missing from `roles`" — §5.1 names no roles, so
for `decision` that clause has no referent. It is the same structural defect the
Revision note says F1-of-round-1 found in 1.0's §8.6: a mandatory rejection clause
with nothing to check against.

The practical consequence is measurable in my own source. My verifier contains
exactly one hardcoded symbol list, `DECISION_ELEMENT_SYMBOLS`, and it exists
solely because §5 left no alternative. **A second `decision` template whose budget
is called `time_left` and whose commitments are `batches` requires editing my
code** — precisely the failure mode §3.3 says is non-conforming. §10 already
concedes the generalisation is "argued, not measured"; for `iteration` the
argument is now backed by a working demonstration (§2.4), for `decision` it is not.

Reproduce: `grep -n 'station_id\|remaining_initial' tests/trace_schema/reviewer_d2_verifier.py`
— the constants block and its ASSUMPTION[6] comment.

**Impact.** Blocks, because it falsifies the specific claim §3.3 makes about the
whole document and because it is the direct cause of F1. Fixing F2 fixes F1: give
`decision` a role table (`index`, `budget`, `budget_initial`, `search_log`,
`committed`, `budget_final`, and for the trial `candidates`, `selected`,
`budget_after`, `closes`), and the missing `roles` key becomes a thing worth emitting.

### F3 — the normative algorithm accepts a `decision` trace with a wrong answer · **CONFIRMED** · blocks merge

`selection.filter` names a boolean `eligible_symbols` member, and §5.2 says only
that "only entries where this is `true` are selectable". **Nothing in §5.2, §5.4,
§6 or §8 says how `filter` is computed, or requires it to be consistent with the
budget it evidently encodes.** §5.4.4 then optimises `criterion` over *the
filtered candidates* — so an entry marked `fits: false` is simply invisible to
every check the spec has.

This is exploitable, and it moves the answer. I hand-built two minimal `decision`
nodes over the same universe `{p: 6, q: 5}` with capacity 11:

- **Honest.** Station 1 takes `p` (6, the max that fits) → 5 left, then `q` (5)
  → 0, then an empty trial closes. `n = 1`.
- **Corrupt.** Station 1's first trial reports `p` as `fits: false`. Selection
  therefore *correctly* picks `q`; `p` spills to station 2. `n = 2`.

Every invariant holds in the corrupt node: budgets thread across trials (§5.4.2),
`remaining_after` arithmetic is exact (§5.4.5), `remaining_final == capacity − Σ`
(§5.4.7), `assigned` matches the ordered `chosen` values (§5.4.6), `closes ==
(chosen is null)` (§5.1), precedences hold both directions (§5.4.3b), the union
is the universe and the sets are disjoint (§5.3), capacity is constant (§5.4.8),
and `chosen` optimises `criterion` over the filtered candidates (§5.4.4).

```
honest (correct answer n=1)                -> PASS  []
corrupt `fits` (answer n=2, WRONG)         -> PASS  []
same corrupt node, --check-filter          -> FAIL  ['F4(non-normative)']
```

(reproduced by the probe block in §4 below; the `--check-filter` rule is **mine,
invented, not in the spec** — `fits == criterion <= budget_before`. It agrees with
all 40 gold traces, which is evidence it is the intended rule and further evidence
that the spec simply forgot to write it down.)

**Impact.** Blocks. `cardinality: "answer_bearing"` means `n` **is** the answer
(§2, §7.3: "a different element count is a different answer"), so §6 admits a
trace that is wrong about the thing being scored. §8.6's "a `chosen` that does not
optimise `selection`" is evadable by construction — you do not have to disobey the
rule, you only have to lie about its input. §5.4.3b's own rationale states the
intent this finding shows is unmet: *"Carrying `precedences` also makes `eligible`
independently checkable rather than trusted, which was the stated design intent of
reporting `fits` separately."* `precedences` closed one half of `eligible`; `fits`
is the other half and is still trusted.

**Fix.** One clause in §5.4: `filter` is a declared predicate over
`(criterion, budget_before)`, checked as `en[filter] == (en[criterion] <=
budget_before)`, or a `selection.filter_relation` field on the same footing as
`update_relation`. The latter is more in keeping with §1's governing rule and
generalises to a template filtering on something other than a budget.

### F4 — `preamble` is unbound to the sequence it starts · **PLAUSIBLE** (I would take it as CONFIRMED)

§4.4 is emphatic that the preamble "is **not** free" and §7.2 lists "preamble
differs from the values the question prescribes" as a **process failure**. But no
invariant anywhere relates `preamble` to `elements[0]`. `carry` covers *k* → *k+1*
only, and §4.3 is per-element. Replacing both preamble frames with nonsense:

```
preamble replaced with nonsense            -> PASS  []
```

1.1 added §4.3.6 for exactly the analogous defect one level down — an evaluation
frame not bound to the update that produced it. The same fix was not applied at
the sequence start, where the same hole is open.

I file it PLAUSIBLE rather than CONFIRMED only because it is arguable that §7.2's
"the values the question prescribes" is a *comparator* obligation against the
question text, which §10 could be read as putting out of scope. But the binding
that *is* in scope — `preamble[-1].iterate == elements[0].iterate_curr` and
`preamble[-2].iterate == elements[0].iterate_prev`, with the residual roles
likewise — is checkable from the node alone, holds on all 40 gold traces, and
costs four lines. **Impact:** a candidate trace can silently restart the scheme
from different trials and still take process credit for element 1.

### F5 — the `eligible` entry's identity symbol is never declared · **PLAUSIBLE**

§5.2 declares `filter` and `criterion` as data, and `chosen` is a string compared
against a candidate's identity and against `termination.universe` — but no field
names the symbol carrying that identity. §5.1 mentions `["task", "duration",
"fits"]` only as "for the reference template". I derived it as `eligible_symbols`
minus `filter` minus `criterion`, which is sound **only** while that list has
exactly three members; a template adding a fourth field breaks it, and a verifier
falls back to the literal `"task"`. Add `selection.identity`. Same class as F2.

### F6 — no node says whether it is gold or candidate, yet normative rules branch on it · **PLAUSIBLE**

At least two rules give opposite dispositions for the same condition:

- §8.1 requires rejecting `satisfied: false`, while §3.4 says a *model's* trace
  "may legitimately carry `false`, and the comparator must be able to say so
  rather than crash".
- §6.4 makes `len(elements) <= max_elements` a hard invariant; §7.2 says
  `len(elements) > max_elements` is "**Reported, not failed**".

The node carries no field distinguishing the two, so the mode is an out-of-band
argument the spec never mentions. I hardcoded gold mode. An implementer who reads
§6 and §7 as one algorithm gets a verifier that is wrong on one of the two.

### F7 — §3.7's `cumulative` reading makes §6.7 undecidable · **PLAUSIBLE**

§3.7 defines `cumulative` as "the union over **all** elements of
`accumulator_symbol` equals `universe`". Under that reading the predicate does not
depend on which element you evaluate it at, so §6.7's requirement that it hold "on
the last element **and on no earlier one**" can never be satisfied by any
sequence of length > 1. I read it as the union over the **prefix ending at element
k**, which makes both halves decidable and holds on all 40 gold traces. §3.7 should
say prefix.

### F8 — §7.1's tolerance boundary is not stated as inclusive or exclusive · **PLAUSIBLE**

"Agrees to half a unit in the last place" does not say whether *exactly* half a
unit agrees. I chose inclusive, on the strength of §4.5 calling `4.985e-05`
conforming "against a 5.0e-05 tolerance". §7.1 says half-way *instances* are
screened at generation (D-016), which is why this is latent on gold — but §7.1 is
the **comparator's** rule, and a candidate trace is not screened. It decides
verdicts the day the comparator runs. One word: "≤".

### F9 — smaller ambiguities I resolved by assumption

Each is recorded as an `ASSUMPTION[n]` comment at the point of use.

9a. §4.4's "`preamble_symbols` … equals `evaluation_symbols`" is stated as a fact,
not as a MUST. I check it; a different implementer would not.

9b. Nothing says a **closing trial must be the last** trial of its element. §5's
"opened, filled, closed" and the budget chain imply it; §5.4.2 does not say it. I
enforce it.

9c. §3.5's path grammar does not say what a path **traversing `null`** means.
Latent today because the only null (`evaluation` on the terminating element) is
never a carry *source*, but it is one schema change away from mattering. I treat
it as a carry failure.

9d. §4.2's grammar admits "role identifiers" without restricting them to *numeric*
roles. `converged` (boolean) and `evaluation` (object) are roles of the node and
are therefore grammatical in `update_relation`. I raise on substitution.

9e. §4.3.7 says "**the** denominator of `update_relation` is non-zero", presuming
the relation has exactly one. Under §4.2's general grammar an expression may have
none or many. I generalised it to "no division by zero occurs during evaluation",
which I believe is intended but which the spec does not say.

9f. §8.8/§6.2 constrain the values of `roles`, but `termination.quantity_role` and
`result.from`'s `last.<role>` **also** name roles and no clause requires *those*
to be roles of the node. A node with `quantity_role: "nonsense"` is rejected by no
§8 clause. I reject it under §3.4 by choice.

9g. Nothing relates `cardinality_symbol` to `result.symbol` (both are `n` in gold)
or to what `result.from: "len(elements)"` produces. §7.3 reports
`cardinality_mismatch` "on `cardinality_symbol`", so they are evidently meant to
be the same quantity.

9h. §6.1's "every conditional field required by the type present" is a one-way
check. It does not say whether an `iteration` node carrying a stray `selection`
key must be rejected. §8.2's key-set rule governs elements and frames, not the
node object itself. I do not reject it.

---

## 4 Falsification attempts that failed

Things I expected to break and which held:

1. **"The `roles` indirection is decorative — the verifier will still need `y`."**
   Falsified, decisively. §2.4's alien node has no symbol in common with the
   reference template and a different recurrence, and passes unmodified. I then
   tried to break it the other way — perturbing one `S_next` by 0.01 — and got
   four independent clause failures (`8.3`, `4.3.1`, `4.3.2`, `4.3.6/8.5`). The
   `iteration` type is genuinely reusable and genuinely tight.

2. **"`update_relation` will need `eval`, or a hardcoded secant."** It does not.
   §4.2's grammar is small enough to parse in forty lines and complete enough to
   express a damped secant with a literal coefficient. This is the single largest
   improvement over what the Revision note describes of 1.0.

3. **"§4.5's headroom number is a rationalisation."** Recomputed independently
   from the corpus: `4.9850e-05` against `5.0e-05`. Exact. The 0.3% margin is
   real, is documented, and would have bitten me — I had considered tightening.

4. **"`symbol_precision` will be missing an entry some check needs."** I made a
   missing entry a hard failure rather than a fallback, specifically to catch
   this. It never fired across 80 traces and 32 mutations. §3.1's enumeration
   ("element symbols, evaluation-frame symbols and trial symbols alike") is complete.

5. **"`result.dp` = 3 against a 4-dp `y_next` is the 1.0 defect again."** It is
   not. The two genuinely disagree on all 40 `iteration` rows, but §7.1 now says
   which to use where ("or `result.dp` for the result"), and §3.6 says the
   derivation is rounded to `result.dp`. Dispatching §6.8 on `result.from` rather
   than on `node_type`, as §3.6 instructs, gives one code path for both types and
   the right number for each. This is a clean fix.

6. **"§8.6's precedence clauses are still unreachable."** They are reachable. I
   constructed violations in both directions (an omitted satisfied item, an
   included unsatisfied one) and both were caught at §5.4.3b. 1.1's `precedences`
   field does what the Revision note claims.

7. **"The two node types will turn out to need one merged code path."** They do
   not, and §2's argument holds up under implementation: `_iteration_locals` and
   `_decision_locals` share the §6.1–§6.5 frame and nothing below it. §7.2 and
   §7.3 dispatch on `cardinality`, not on `node_type`, which is the right seam.

8. **"An empty `eligible` on the last trial will need a special case."** It does
   not — §5.1's restated single `closes` rule handles it, and my precedence check
   independently *requires* it to be empty there, since no uncommitted item
   remains. The "no exception is needed and none is granted" paragraph is correct.

The two probes that **did** find something (§3's F3 and F4) are reproduced by:

```
PYTHONIOENCODING=utf-8 python - <<'EOF'
import json, copy
from tests.trace_schema.reviewer_d2_verifier import verify_node
rows = json.load(open('docs/re-implementation-sep/phase3_conformance/traces.json', encoding='utf-8'))
it = next(r['trace_nodes'] for r in rows if r['trace_nodes']['node_type'] == 'iteration')
n = copy.deepcopy(it)
n['preamble'][0]['y'] = 99.0; n['preamble'][1]['y'] = 88.0     # F4
print('preamble replaced ->', verify_node(n, 'F4').status)
EOF
```
and, for F3, by the two hand-built nodes described in §3 (`p=6, q=5, capacity=11`);
the corrupt one passes with `--lenient-roles` semantics and fails only under my
non-normative `--check-filter`.

---

## 5 Further probing and improvements

### 5.1 Thinnest evidence in what I did

- **`decision` was verified under a waiver.** Every `decision` result above is
  conditional on suspending §3.1's `roles` requirement (F1). I believe the waiver
  is the spec's error and not mine, but a reader should know that "80/80" is
  "40/40 clean plus 40/40 with one field waived".
- **My `decision` checks are not independent of `decision`'s symbol names.** F2
  means I cannot claim to have tested the `decision` *type*; I tested one
  template's shape. The `iteration` result is much stronger evidence than the
  `decision` result, and the two should not be quoted at the same confidence.
- **The alien node is synthetic and mine.** It proves a verifier written from §4
  generalises across symbol names and recurrences; it does not prove
  `linear_reservoir_routing_step` fits, which §10 correctly still lists as
  unmeasured.
- **`tie_break` was never exercised on real data**, as §5.5 predicts. My
  `universe_order` and `first` paths ran only on hand-built ties.

### 5.2 What I would probe next, ranked

1. **Fit `linear_reservoir_routing_step` to `iteration` by declaring a binding
   only** — no verifier edit permitted. §10 says this is the point of `roles` and
   `update_relation`. My §2.4 result says it should work. Actually doing it turns
   §10's "argued, not measured" into measured, and it is now cheap.
2. **Write the F3 fix and re-run.** A `filter_relation` on the footing of
   `update_relation` is a half-day and closes the one hole that lets a wrong
   answer through.
3. **Give `decision` a role table and fit a second `decision`-shaped template.**
   This is F2's real test, and until it happens `decision` should be documented as
   template-specific rather than as a type.
4. **Run the verifier against deliberately imperfect *model* traces**, not
   mutations of gold. Everything above tests §6; §7's comparator semantics — the
   `incidental` count-mismatch rule, `cardinality_mismatch`, process-vs-answer
   credit — has **no conformance corpus at all** and is therefore the least
   evidenced part of the document. F6 is the first thing that would bite.
5. **Property-test the two grammars.** §4.2 and §3.5 are small enough to fuzz;
   a few hundred random malformed expressions would confirm no implementer's
   parser silently accepts something outside them.

### 5.3 What generalises

- **"Anything a verifier must check is declared on the node"** (§1) is the right
  governing rule, and where 1.1 applied it, it worked: `update_relation`, `roles`,
  `precedences`, `selection`, `rounding`, `symbol_precision`, `result.from`. Every
  remaining CONFIRMED finding is a place the rule was *not* applied — `filter`
  (F3), `decision`'s symbol names (F2), the preamble binding (F4), the candidate's
  identity symbol (F5). The rule is a reliable defect predictor: **grep the spec
  for any quantity a check consumes that is not a field, and that is the next bug.**
- **Splitting checkable from uncheckable** (`carry` vs `carry_notes`,
  `_prose`/`_notes` suffixes) is worth carrying to every later phase. It made the
  §8.4 quarantine channel trivial to implement correctly.
- **Machine-generated worked examples** (§9) is right and should be standard — but
  F1 shows it must compare the corpus to the **field table**, not just the example
  to the corpus.
- **Writing down tolerance headroom** (§4.5) prevented a real error in my
  implementation. Later phases should do it wherever a tolerance is nominally slack.

### 5.4 What later phases should do differently

- **State the verifier's *mode* as data.** F6 is a structural gap, not a typo: a
  document that specifies both a gold checker and a candidate comparator must say
  which rules belong to which, and ideally carry it on the node.
- **Head every normative table "Role", or say explicitly that the type is
  single-template.** The §4/§5 asymmetry is invisible until you implement, and it
  is the whole difference between a type and a description.
- **Make §8 a list of *traces*, and put behavioural obligations elsewhere.** §8.4
  is a verifier requirement wearing a trace requirement's clothing, and it is the
  one an implementer will skip.
- **Do not let a `_prose` field state a constraint that no data field carries.**
  §4.4's "the preamble is not free" (F4) and §5.2's silence on `filter` (F3) are
  both places where the prose asserts an obligation the data cannot support.

### 5.5 On the gate

Phase 3's exit gate reads *"Reviewer D implemented a working verifier from D3.3
alone."* On the evidence above I would record it as **met, with the two CONFIRMED
findings as merge conditions**. I implemented it alone, without the template
source, in the time box, and the one place I had to reach outside the spec (F3's
`filter` rule) I marked non-normative rather than quietly adopting — which is
itself the finding. The gate's harder implicit question, *is this a type*, splits:
**yes for `iteration`, demonstrated adversarially; not yet for `decision`.**

---

*Filed by Reviewer D2. Written to disk, not committed — the implementer commits.
Artefacts: `tests/trace_schema/reviewer_d2_verifier.py`,
`tests/trace_schema/reviewer_d2_negatives.py`.*

---

## 6. Round 3 — disposition against schema 1.2

**Ref:** `c9fd703e65ba32f10e78bfc6b6a2ef45ee178d34` on `redesign/phase3-trace-shape`.
`HEAD`, the branch tip and the commissioned SHA agreed at the start of this round
and still agreed at the end — **the branch did not move under me this time.**
**Time box:** 30 minutes, re-verification only. Sections 1–5 above are the
round-2 record against 1.1 and are unchanged.

**What I built.** `tests/trace_schema/reviewer_d2_verifier_12.py` — the round-2
verifier adapted to 1.2, reusing only its grammar and rounding helpers. Round 2's
file is left untouched so the round-2 evidence stays reproducible. Every check is
rewritten against 1.2's role maps.

```
python -m tests.trace_schema.reviewer_d2_verifier_12 docs/re-implementation-sep/phase3_conformance/traces.json
```
```
template_line_balancing_heuristic        pass 40  fail  0
template_normal_depth_iteration          pass 40  fail  0
TOTAL 80 rows: 80 pass, 0 fail; 40 carried unchecked notes
```

**80/80 with no waiver.** Round 2 needed `--lenient-roles` to get past F1; 1.2
needs nothing. That alone settles F1.

### 6.1 Disposition table

| # | Finding | Disposition |
|---|---|---|
| F1 | `roles` required but absent on all `decision` nodes | **ADDRESSED** |
| F2 | `decision` is a description, not a type | **ADDRESSED** |
| F3 | filter trusted, so a wrong answer passes | **PARTIALLY ADDRESSED — still merge-blocking** |
| F4 | `preamble` unbound to `elements[0]` | **PARTIALLY ADDRESSED** |
| F5 | candidate identity symbol undeclared | **ADDRESSED** |
| F6 | no gold-vs-candidate mode | **ADDRESSED** |
| F7 | §3.7 `cumulative` makes step 7 undecidable | **NOT ADDRESSED** |
| F8 | tolerance boundary unstated | **ADDRESSED** |
| F9a | `preamble_symbols == evaluation_symbols` stated as fact | **NOT ADDRESSED** |
| F9b | closing trial not required to be last | **ADDRESSED** |
| F9c | path traversing `null` undefined | **NOT ADDRESSED** |
| F9d | `update_relation` admits non-numeric roles | **NOT ADDRESSED** |
| F9e | "*the* denominator" presumes exactly one | **NOT ADDRESSED** |
| F9f | `quantity_role` / `result.from` roles unconstrained | **NOT ADDRESSED** |
| F9g | `cardinality_symbol` unrelated to `result.symbol` | **NOT ADDRESSED** |
| F9h | stray conditional fields not rejected | **NOT ADDRESSED** |

**Five of my eight F9 assumptions are still assumptions** — seven of eight
ambiguities carried forward in all, but every one of them is minor and none blocks.

### 6.2 The three things you asked me to break

#### 1. The F3 exploit — the fix is real, and my exploit routes around it

`selection.filter_relation` and §5.4.3c work exactly as advertised for the attack
I filed. Rebuilding the round-2 pair against 1.2, plus one new variant:

```
honest  n=1 (correct)                n=1 -> PASS   []
corrupt A: lie about `admissible`    n=2 -> FAIL   ['5.4.3c/8A.5', '5.4.3c/8A.5']
corrupt B: lie about `measure`       n=2 -> PASS   []
```

**Corrupt A is dead** — rejected, and rejected for the *right* clause (§5.4.3c:
`admissible: false` against `filter_relation` `measure <= budget_before`, which
gives `true` at `6 <= 11`). That is a genuine soundness fix and it is not cosmetic.

**Corrupt B is the variant you asked for, and it works.** Instead of lying about
the flag, lie about its *input*:

- Station 1, trial 1: candidates `p` with **`measure: 200`** and
  `admissible: false`, `q` with `measure: 5` and `admissible: true`, budget 11.
  §5.4.3c is *satisfied* — `200 <= 11` is genuinely false, so the flag is honest
  about a dishonest number. Selection correctly takes `q`.
- Station 1, trial 2: `p` at `measure: 200` again, inadmissible at budget 6, so
  nothing is selectable and the station closes with `committed: ["q"]` and
  `budget_final = 11 − 5 = 6`.
- Station 2: `p` reappears with **`measure: 6`**, admissible, committed;
  `budget_final = 11 − 6 = 5`.

Every clause of §5.3, §5.4 (including 3b, 3c, 5, 7 and 9), §6 and §8A holds. The
reported answer is `n = 2`; the truth is `n = 1`.

**Root cause: nothing anywhere requires an item's `measure` to be the same in
every trial it appears in.** I searched the whole of 1.2 for a `measure`
invariance clause and there is none — §5.1 defines `measure` as "the quantity
`selection.criterion` optimises and the budget is spent in", and §5.4.5 and
§5.4.7 constrain it only for candidates that are actually *chosen* or
*committed*. A candidate you intend to exclude has a free `measure`.

1.2 moved the trust boundary from the flag to the flag's input without closing it.
Filed as **F3'**, and I would keep the merge blocked on it: the primary defect —
*a `decision` trace with a wrong `n` passes the normative algorithm* — is still
reproducible, in a form that takes about six lines to construct.

**The fix is one clause, and it is smaller than `filter_relation` was:**

> **§5.4.3d.** For any item appearing as a candidate in more than one trial, its
> `measure` role is identical in every occurrence. An item's measure is a
> property of the item, not of the trial that considers it.

Optionally stronger, and closing the class rather than the instance: declare the
measures once on the node (`termination.measures`, alongside `precedences` — both
are properties of `universe`) and require every candidate occurrence to agree with
it. That also makes `budget_final` checkable without scavenging measures out of
the trial log, which is what my §5.4.7 implementation currently has to do.

#### 2. The rename test on `decision` — passes

Renamed every symbol at all three nesting levels, rebound only `roles`,
`trial_roles`, `candidate_roles` and the declared symbol lists, and left
`selection` byte-identical:

| Level | Renames |
|---|---|
| element | `station_id→s`, `capacity→cap0`, `remaining_initial→b0`, `trials→log`, `assigned→done`, `remaining_final→b1` |
| trial | `eligible→opts`, `chosen→pick`, `remaining_after→b_after`, `closes→shut` |
| candidate | `task→id`, `duration→w`, `fits→ok` |

```
selection untouched: {"filter": "admissible", "filter_relation": "measure <= budget_before",
                      "criterion": "measure", "objective": "max", "tie_break": "universe_order"}
renamed                            -> PASS+UNCHECKED []
renamed + perturbed budget_final   -> FAIL ['5.4.7']
renamed + admissible flipped       -> FAIL ['5.4.3c/8A.5']
```

**Accepted after the rename, and still rejected when a value is perturbed** — both
halves of §8B.11's requirement. My verifier needs no literal symbol name anywhere
for `decision`. **F2 is genuinely fixed**, and 1.2's `decision` is now a type in
the same sense §4's `iteration` already was. This is the round's best result.

One wart the rename surfaced: **`termination.accumulator_symbol` is a raw symbol**
(`"assigned"`) in an otherwise all-roles §5, so I had to rename it too. It is
declared data, so no verifier hardcodes it and §8B.11 is not violated — but it is
the one place §5 still speaks symbols, and it should be
`accumulator_role: "committed"` for consistency. Filed as **F11**.

#### 3. §8B.11 self-audit — my verifier passes, and the grep is not what proves it

Grepping `reviewer_d2_verifier_12.py` for the reference template's symbols returns
hits on exactly three: `trials`, `chosen`, `closes`. **All three are false
positives, and the reason matters.** §5.1's role column and its reference-symbol
column are *identical strings* for those three roles (`trials`→`trials`,
`chosen`→`chosen`, `closes`→`closes`). Every occurrence in my code is a role
lookup — `E(el, "trials")`, `T(tr, "chosen")`, `T(tr, "closes")` — or a Python
local variable name. No element, trial or candidate value is resolved by a literal
symbol name.

But the grep alone could not have told you that, and here is the finding: **because
3 of `decision`'s 10 roles are spelled identically to their reference symbols, a
verifier that genuinely hardcoded `el["trials"]` would pass the entire reference
corpus and fail only on a renamed node.** Static audit of §8B.11 is therefore
unsound by construction; the rename test is the only reliable check. §8B.11
already says as much ("§10 names this the rename test") — I am confirming that the
weaker method really is inadequate rather than merely inelegant, and suggesting the
reference symbols be chosen to differ from the role names in a future revision so
that both methods work.  Filed as **F13**.

My evidence that I pass is therefore the rename test on **both** node types
(§2.4 for `iteration` at 1.1, §6.2.2 for `decision` at 1.2), not the grep.

### 6.3 Findings that survive, with reproduction

**F3' — a `decision` trace with a wrong answer still passes · CONFIRMED · blocks merge.**
Reproduction: the three-node script in §6.2.1; `corrupt B` returns `PASS` with
`n = 2` against a true `n = 1`. Fix: §5.4.3d above.

**F4 — `preamble_binding` closes the hole, but is escapable by omission · PARTIAL.**
```
preamble replaced with nonsense        -> FAIL ['8A.4', '8A.4']
preamble_binding field simply omitted  -> PASS []
```
The check works. But §3.1 marks `preamble_binding` merely "cond. — `iteration`
only" and nowhere says **required iff `preamble` is present**, so a trace that
drops the field escapes §8A.4 entirely. §3.1 already has the right pattern two
rows above (`preamble_symbols`: "Required iff `preamble` is present"); apply the
same wording. One line, and F4 is then fully addressed.

**F7 — not addressed, and I checked rather than accepted the claim · CONFIRMED (minor).**
You asked whether step 7 is now decidable or whether it had only been asserted.
It had only been asserted. §3.7 is **byte-identical** between 1.1 and 1.2:
```
git show 9398ccc:...phase3_node_types.md | sed -n '/^### 3.7/,/^---/p' | sha256sum   -> 6f7283ead30df88e
git show c9fd703:...phase3_node_types.md | sed -n '/^### 3.7/,/^---/p' | sha256sum   -> 6f7283ead30df88e
```
§6 step 7 is likewise unchanged. `cumulative` still reads "the union over **all**
elements", under which the predicate does not depend on which element you evaluate
it at, so "and on no earlier one" remains unsatisfiable for any sequence longer
than one. My 1.2 verifier still carries `ASSUMPTION[5]` and still reads it as the
prefix union. Harmless in practice, wrong as written; one word ("prefix") fixes it.

**F6, F8 — addressed, verified rather than assumed.**
```
satisfied:false, gold mode        -> FAIL ['8A.1']        budget overrun, gold      -> FAIL ['6.4']
satisfied:false, candidate mode   -> PASS+UNCHECKED       budget overrun, candidate -> PASS+UNCHECKED
change off by exactly 0.5e-4      -> PASS                 change off by 0.6e-4      -> FAIL ['4.3.2']
```
§6.0's mode argument does what it says, and §7.1's inclusive boundary is
implementable exactly as written.

**F9a, F9c–h — still assumptions.** §4.2 is byte-identical to 1.1; §3.5 and §4.3
changed only in cross-references (`§8.4`→`§8B.10`, and §4.3.8's parenthetical),
not in the clauses my assumptions attach to. Demonstrated by mutation:
```
preamble_symbols diverges from evaluation_symbols  -> PASS   (F9a: nothing forces it)
iteration node carrying a `selection` key          -> PASS   (F9h: stray fields unrejected)
iteration node carrying `capacity_constant`        -> PASS   (F9h)
cardinality_symbol != result.symbol                -> PASS+UNCHECKED (F9g)
quantity_role names a non-role                     -> FAIL   (F9f: my choice, no §8A clause)
result.from names a non-role                       -> FAIL   (F9f: my choice, no §8A clause)
```
None blocks. F9f is worth one clause in §8A.7 ("…or a role named by
`termination.quantity_role` or `result.from` that is not a key of `roles`"),
since I currently reject those on my own authority and a different implementer
would not.

### 6.4 New in round 3

- **F10 — `preamble_binding` does not declare which *frame* role each element role
  binds to.** The map is `{"iterate_prev": 0, "residual_prev": 0, …}` — element
  role to frame *index*. §4.4's prose pairs `iterate_prev` with the frame's
  `iterate` role and `residual_prev` with `residual`, but that pairing is nowhere
  data. I derive it by **name prefix** (`ASSUMPTION[10]`), which works only
  because 1.2's element roles happen to begin with the frame role's name. A
  template whose roles were `lower`/`upper` would break it. Make the value a pair,
  `{"iterate_prev": [0, "iterate"]}`, or state the prefix rule normatively.
- **F11 — `termination.accumulator_symbol` is the last raw symbol in §5.** Should
  be `accumulator_role`. Cosmetic, but it is the one place the role doctrine leaks.
- **F12 — §3.1's `schema_version` row still reads `"1.1"`** while the corpus and
  the document's own header say `1.2`. A one-character staleness, but it is in the
  **normative field table**, and an implementer who takes it literally rejects all
  80 traces at step 1. I very nearly did.
- **F13 — 3 of `decision`'s 10 role names are spelled identically to their
  reference symbols**, which makes static §8B.11 auditing unsound. See §6.2.3.

### 6.5 Round-3 verdict

**The two claims I was asked to test hardest split.**

**F2 is genuinely fixed** and the primary claim of the whole document — that these
are types, not descriptions — now holds for **both** node types, demonstrated by
the rename test rather than argued. That is a real result and 1.2 earns it.

**F3 is not fixed, only narrowed.** `filter_relation` and §5.4.3c close the exact
attack I filed, and close it properly, but the class is open: the algorithm still
accepts a `decision` trace whose `n` — the answer — is wrong, now by lying about
`measure` rather than about `admissible`. I would **not** merge on this. It is a
one-clause fix (§5.4.3d, or `termination.measures`), and given that this is the
second revision to close this hole one variant at a time, I would rather see the
invariant stated generally — *every value a check consumes is either declared once
on the node or recomputed from something that is* — and the node audited against
it, than see a third round close a third variant.

Everything else outstanding (F4's escapable optionality, F7's one word, F9a and
F9c–h, F10–F13) is minor and could ride along with that clause in a single pass.

*Round 3 filed by Reviewer D2. Written to disk, not committed. New artefact:
`tests/trace_schema/reviewer_d2_verifier_12.py`.*

---

## 7. Round 4 — does §3.8 close the exploit family?

**Ref:** `9c6b8ee9e9fb6eb2cfcea6af30e2d787694b1e9a`. `HEAD`, the branch tip and the
commissioned SHA agreed at the start and at the end — the branch did not move.
**Time box:** 30 minutes. Artefact: `tests/trace_schema/reviewer_d2_verifier_13.py`.

### 7.1 Verdict

**A fourth variant exists. In fact two do, and the second is worse than anything
I have filed so far.** §3.8 is a real closure of the family it names — *values
restated across occurrences of the same thing* — and within that family I could
not break it. But it is not a closure of the question it is written to answer,
because **its own table omits two values that checks read and traces choose
freely**, and both change the answer.

- **Variant 4 (`decision`): `budget_total` is restated in every element.** A
  trace that declares a different budget reports `n = 2` where the truth is
  `n = 1`, passing every clause. §3.8's table does not list it.
- **Variant 5 (`iteration`): the residuals are fabricated.** A fully conforming
  trace reports **3.501 m** where seed 0's gold answer is **1.380 m**. §10's
  "frame-internal physics is not checkable" concession is not merely a limitation
  — **it is load-bearing for an attacker**, and it makes `result.value` on an
  `iteration` node entirely free.

I would hold the merge. Not because §3.8 is wrong — it is the right rule, and
adopting it was the right call — but because the audit it prescribes has not yet
been *run* against the two node types it ships with.

### 7.2 What 1.3 got right — confirmed, not assumed

**Corrupt B is dead**, as expected, and for the right clause:
```
corrupt B (lie about measure) -> FAIL ['5.4.3d/8A.5']
```
`termination.item_measures` + §5.4.3d closes variant three exactly as designed.

**80/80 gold traces pass, no waiver:**
```
python -m tests.trace_schema.reviewer_d2_verifier_13 docs/re-implementation-sep/phase3_conformance/traces.json
template_line_balancing_heuristic        pass 40  fail  0
template_normal_depth_iteration          pass 40  fail  0
TOTAL 80 rows: 80 pass, 0 fail; 40 carried unchecked notes
```

**Rename test, both types, no regression from 1.2:**
```
decision:  renamed -> PASS+UNCHECKED    + perturbed budget_final -> FAIL ['5.4.7']
                                        + perturbed measure      -> FAIL ['5.4.3d/8A.5']
iteration: renamed -> PASS              + perturbed iterate_next -> FAIL ['8A.3','4.3.1']
                                        + preamble nonsense      -> FAIL ['8A.4']
```
`accumulator_role` now survives the rename untouched (it is a role, not a symbol),
which is F11's fix showing up as behaviour rather than as prose.

**The five fixes you asked me to verify rather than accept — all five are real.**

| # | Check | Result |
|---|---|---|
| F7 | §3.7 hashed at all three refs | `9398ccc` `6f7283ead30df88e` → `c9fd703` `6f7283ead30df88e` → `9c6b8ee` **`e6ce5f30f3045d06`**. The hash moved this time, and the text now reads "*the union of `accumulator_role` over elements 1..k*" with an explicit "**The reading is a prefix, not a total**". **Fixed.** |
| F4 | §3.1 `preamble_binding` row | now "**required iff `preamble` is present**"; omitting it is rejected. **Fixed.** |
| F10 | `preamble_binding` shape | now `{"frame": 0, "frame_role": "iterate"}` — the frame role is named, not prefix-inferred. My `ASSUMPTION[10]` is deleted. **Fixed.** |
| F11 | `accumulator_symbol` | now `accumulator_role`, an element role. **Fixed** (one stale reference remains — see §7.5). |
| F12 | §3.1 `schema_version` row | now `"1.3"`. **Fixed.** |

### 7.3 Variant 4 — `budget_total` is restated in every element

Working down the node as you asked: `index`, `closes`, `chosen`, `candidates`
membership, `budget_after` and `budget_final` are all pinned. `index` is
contiguous (§5.3), `closes == (chosen is null)` and must be last (§5.4.9),
`chosen` is a replay of `selection` (§5.4.4), `candidates` is pinned in both
directions by `precedences` (§5.4.3b), and the budget chain is arithmetic on
`item_measures` (§5.4.2/5/7). Given `universe`, `precedences`, `item_measures` and
the budget, the greedy is deterministic and exactly one trace conforms. That is a
real result and it is what §3.8 bought.

**But the budget itself is declared in each element, and nothing says what it is.**

```
honest, budget 11, n=1 (correct)                     -> PASS
B1  capacity_constant=false, budgets 6 and 5, n=2    -> PASS []
B2  capacity_constant=true,  budget 6 everywhere, n=2-> PASS []
```

Both report `n = 2` against a true `n = 1`. Universe `{p: 6, q: 5}`,
`item_measures` honest and untouched, every candidate's `measure` agreeing with it,
`filter_relation` honestly recomputed, precedences satisfied both ways, budget
arithmetic exact, closing trials last, prefix-cumulative termination holding only
at the last element. `n` is the answer (§2, §7.3), and it is wrong.

**B2 is the important one.** B1 could be dismissed by saying `capacity_constant`
is itself a free boolean — which is true, and is its own small finding. But B2
leaves `capacity_constant: true` and simply restates the *same wrong budget* in
every element. It is the `measure` attack with `budget_total` substituted for
`measure`: a value read by §5.4.1, by §5.4.3c through `budget_before`, and by
§5.4.7, declared once per element rather than once per node.

**§3.8's own audit question convicts it verbatim.** §3.8 says: *"for every value a
check reads, where is the single place it is declared? If the answer is 'in each
element that mentions it', the check is not yet a check."* For `budget_total` the
answer is, exactly, "in each element that mentions it".

**Fix, in §3.8's own idiom:** declare `termination.budget` once beside
`item_measures` and `precedences` — all three are properties of the problem, not
of an element — and make §5.4.1 read `budget_total == termination.budget`.
`capacity_constant` then becomes redundant and should be deleted rather than left
as a free boolean that a trace can flip to unlock per-element budgets.

### 7.4 Variant 5 — fabricated residuals make `iteration`'s answer free

You asked whether §10's concession is load-bearing for an attacker. **It is, and
this is the most serious thing in any of my four rounds.**

§4.3.1 recomputes `iterate_next` from `update_relation` — but the inputs it
recomputes *from* are `residual_prev` and `residual_curr`, which arrive from
evaluation frames, and §10 states that frame-internal physics is not checkable.
So the residuals are free data. The secant step is an interpolation between two
residuals: **choose the residuals and you choose where it converges.**

I built such a trace. At each element I picked the iterate I wanted next, solved
the secant relation backwards for the `g` that produces it, rounded that `g` to its
declared 3 dp, and recomputed `iterate_next` forward from the rounded values —
exactly the arithmetic a real template does. The chosen `g` is then placed in the
previous element's evaluation frame, so `carry`'s `g_curr <- evaluation.g` holds.

```
gold seed-0 answer is 1.380 m; this trace reports 3.501 m in 9 updates
verdict -> PASS []
```

Everything holds: §4.3.1 (the update recomputes exactly), §4.3.2 (`change`),
§4.3.3–5 (`converged`, nullity), §4.3.6 (each frame bound to the iterate that
produced it), §4.3.8 (contiguous index), `carry` on all four paths,
§4.4/`preamble_binding` against `elements[0]`, termination at the last element only,
and §6 step 8 (`result.value` = `last.iterate_next` rounded to 3 dp). **The answer
is off by 2.121 m and nothing objects.**

Two things make this worse than variant 4:

1. **It needs no restatement at all.** Every value appears once. It is not a §3.8
   violation on §3.8's own terms — which is precisely why §3.8 does not catch it,
   and why "declared once" is necessary but not sufficient. A value declared once
   and *never checked against anything* is as free as one restated ten times.
2. **It survives the comparator's process credit.** §7.2 grants process credit
   when "every candidate element satisfies §4.3", independently of the count. A
   fabricated trace satisfies §4.3 completely, so it takes **full process credit
   while reaching an arbitrary answer**. Answer credit catches it by comparing
   `result.value` to gold, but the milestone model's whole purpose is to score the
   *work*, and here the work is unscored. A model that learns to emit
   arithmetically-consistent nonsense is rewarded by exactly the mechanism built to
   detect it.

**§10 and §3.8 are in direct contradiction, and the document does not notice.**
§10 says the geometry is the item's content and carrying it "would mean carrying an
expression language for arbitrary engineering formulae". §3.8 says a check that
reads a freely-restated value is not a check. §4.3.1 reads the residuals. Both
cannot stand.

**The fix is cheaper than §10 implies**, because the node already has the
expression language §10 says it lacks: `update_relation` is one, and 1.2 added a
second for `filter_relation`. A `residual_relation` over evaluation-frame roles
plus item constants declared once (`termination.constants` or similar) would close
it with no new machinery — `g = Q_target − K·AR^(2/3)` needs exponentiation added
to §4.2's grammar and nothing else. If that is genuinely out of scope for Phase 3,
then **§10 should say plainly that an `iteration` node's answer is not verifiable
from the node**, and §7.2 should withhold process credit rather than grant it,
because at the moment the spec implies a guarantee it does not provide.

### 7.5 Does §3.8's table have a gap? — yes, two

§3.8's table accounts for: item measures, admissibility, eligibility, the update
recurrence, display precisions, rounding mode, symbol names. Checks also read:

| Value read by a check | Declared once? | In the table? |
|---|---|---|
| `budget_total` | **no — once per element** | **no** (variant 4) |
| evaluation-frame residuals | once, but checked against nothing | **no** (variant 5) |
| `universe`, `precedences`, `item_measures` | yes, on `termination` | partly |
| `capacity_constant` | free boolean, unconstrained | no |
| `termination.tolerance`, `max_elements` | once, unconstrained | no |

The last three are the **declaration family** and I want to be precise about scope,
because it would be cheap to inflate the finding: a node-local verifier cannot know
that a trace's `universe` or `tolerance` is the *item's* — nothing binds the node
to the question, which §10 (D-039) already records. Those are the comparator's job
under §7, where a candidate is scored against a gold node that carries the true
values. I am **not** filing them as defects.

Variants 4 and 5 are different, and that is why they are findings:

- **Variant 4 is a restatement**, per element, of a value that has a single true
  value — the exact shape §3.8 was written to forbid, missed by its own table.
- **Variant 5 needs no restatement**: it exploits a value that is declared once and
  checked against nothing. It shows **§3.8's rule is incomplete as written.** The
  principle should read: *every value a check consumes is either declared once on
  the node **and checked against something**, or recomputed from something that is.*
  "Declared once" was the right fix for variants two and three because those values
  were also recomputable; it is not sufficient in general.

**One residual staleness:** §6 step 7 still names `accumulator_symbol` while §3.4,
§3.7 and §9.2 all say `accumulator_role`. Harmless — it is the algorithm's prose,
not the field table — but it is the same class of miss as F12, and it is the third
time a rename has been applied everywhere except the normative algorithm.

### 7.6 What I tried that did not work

A clean verdict is worth what the failed attempts behind it are worth, so:

1. **Padding `decision` with an extra empty element** (`committed: []`, budget
   untouched, index contiguous, union still correct). This would have been the
   cleanest possible `n` attack and it is **dead** — but *only because 1.3 fixed
   F7*. Under 1.1/1.2's "union over all elements", the predicate is constant in `k`,
   "and on no earlier one" is vacuous, and the padded trace passes. F7 was filed as
   a cosmetic wording bug in rounds 2 and 3; it turns out to have been load-bearing
   for soundness, and 1.3 closed a variant it did not know it was closing.
2. **Restating `measure` under the new §5.4.3d** — dead, confirmed above.
3. **Lying about `admissible` under §5.4.3c** — dead since 1.2.
4. **Reordering `candidates`** to steer `tie_break: universe_order` — §5.4.4
   resolves ties by position in `universe`, not by position in `candidates`, so the
   scan order cannot be weaponised. `tie_break: "first"` *would* be steerable, but
   no gold node uses it.
5. **Omitting a candidate to dodge selection** — §5.4.3b's second direction
   (every uncommitted precedence-satisfied item must appear) rejects it.
6. **Splitting one element's trials to hide a commitment** — §5.4.6 (`committed`
   is exactly the ordered non-null `chosen` values) and §5.4.9 (closing trial last)
   between them leave no room.
7. **`iteration`: adding or removing updates** — genuinely not an attack, and
   correctly so: §7.2 makes the count incidental by design, and §4.3.4/4.3.5 pin
   `converged` and the frame nullity to the last element regardless of count.
8. **`iteration`: perturbing `result.dp` or the rounding mode** to shift the
   answer — both are declared on the node and §6 step 8 recomputes from them, so
   changing them changes the derivation too and the check still holds. This is
   §3.8 working.

### 7.7 Recommendation

Two clauses, both in §3.8's own idiom, and I believe they finish it:

1. **`termination.budget`, declared once**; §5.4.1 becomes
   `budget_total == termination.budget`; delete `capacity_constant`.
2. **Either** add `residual_relation` (plus exponentiation to §4.2's grammar and
   item constants declared once), **or** state in §10 and §7.2 that an `iteration`
   node's answer and process are not verifiable from the node alone, and withhold
   process credit accordingly.

And one to the rule itself: amend §3.8 to *"declared once **and checked against
something**, or recomputed from something that is"*, then **run its own audit
against both shipped node types** — that audit is what would have caught both
variants, and it has not yet been performed on the two types the document defines.

I said in round 3 that I would rather see the invariant stated generally than watch
a third variant get patched. 1.3 stated it, and stating it was right — the
restatement family really is closed, and `decision` is now deterministic given its
declarations. What remains is that the rule was written down but not yet applied
to the node types in the same document. That is a smaller gap than the one it
closed, and I do not think it needs a fifth round of review — it needs the audit
run once, by the implementer, against §3.8's own table.

*Round 4 filed by Reviewer D2. Written to disk, not committed. New artefact:
`tests/trace_schema/reviewer_d2_verifier_13.py`.*

---

## 8. Round 5 — attacking §4.6 `constants` + `frame_relations`

**Ref:** `fe2040a747f4af503774ea5909cc85746a22e8ab`, schema 1.4. `HEAD`, the branch
tip and the commissioned SHA agreed at start and end. **Time box:** 40 minutes.
Artefact: `tests/trace_schema/reviewer_d2_verifier_14.py`, including a from-scratch
implementation of §4.6's grammar and evaluation rule.

### 8.1 Verdict

**§4.6 is not sound as written. One gap is merge-blocking, and it restores variant
5 in full.**

`frame_relations` **has no coverage requirement.** §4.3.9 quantifies over *the
relations present*, not over *the frame's symbols*. Delete the one relation that
defines the residual and every other clause still passes — including §4.3.9, which
now vacuously checks `A`, `P` and `AR` and says nothing about `g`. The residual is
free again, and so is the answer:

```
frame_relations covers ['A','P','AR'] -- 'g' omitted
target 3.5 -> 9 updates, answer 3.492 m (gold 1.380 m) : PASS []
target 2.0 -> 10 updates, answer 2.000 m (gold 1.380 m) : PASS []
```

**That is variant 5, unchanged, reached by deleting one line of the node.** The
fix that closed it is real — I could not defeat it while it was present — but it is
opt-in, and §8A.4 rejects "a frame that does not satisfy `frame_relations`" rather
than "a frame symbol with no relation".

The good news is that everything else in §4.6 held. I attacked the grammar, the
rounding rule, the ordering discipline and the constants, and **only the coverage
gap yields an arbitrary answer.** Three smaller findings and five stated
ambiguities follow; none of them blocks.

**One clause fixes the blocker:**

> **§4.6, coverage (normative).** `frame_relations` **must** name every member of
> `evaluation_symbols` except the `iterate` role, exactly once, and a node whose
> relations do not cover them is rejected (§8A.4). A frame symbol with no relation
> is a free value, and §3.8 forbids one.

### 8.2 Confirmations

**80/80 gold traces, no waiver**, with §4.6 implemented from the spec text alone:
```
python -m tests.trace_schema.reviewer_d2_verifier_14 docs/re-implementation-sep/phase3_conformance/traces.json
template_line_balancing_heuristic        pass 40  fail  0
template_normal_depth_iteration          pass 40  fail  0
TOTAL 80 rows: 80 pass, 0 fail; 40 carried unchecked notes
```

**Variant 4 is dead.** `termination.budget` is present (147), `capacity_constant`
is gone, and restating every `budget_total` as 6 is rejected:
```
every budget_total restated as 6 -> FAIL ['8A.5']
```

**Rename test, both types, no regression.** For `iteration` the rename now has to
rewrite `frame_relations`' expression strings too, since they name frame symbols
directly:
```
frame_relations: [["a_","(b + z * S) * S"], ["p_","b + 2 * S * sqrt(1 + z ** 2)"],
                  ["ar_","a_ * (a_ / p_) ** (2 / 3)"], ["r","round(ar_, 3) - K"]]
renamed                            -> PASS
renamed + AR perturbed by 0.002    -> FAIL ['4.3.9/8A.4']
renamed + residual perturbed 0.01  -> FAIL ['8A.3']
```
Frame symbols are the one part of an `iteration` node with **no role indirection**
— `A`, `P`, `AR` have no roles, only names. That is defensible (they are item
content, and the relations that name them are data), but it is worth one sentence
in §4.6 saying so explicitly, because it is the only place a renamer must edit an
expression rather than a map.

**Your staleness note, verified cheaply.** I built the 1.4 verifier as a patch of
the 1.3 one. The only semantic edits required were exactly the two you named —
`capacity_constant` → `termination.budget`, and `preamble_binding`'s reshaping —
plus the new §4.6 machinery. Nothing else in 1.3's logic needed changing, which
corroborates your reading that the older verifiers' failures are version staleness
and not semantics. I did not spend further box on it.

**§10's withdrawal is right.** "Frame-internal physics is not checkable" was never
a limitation of the node, only of what the node chose to carry; the moment §4.3.1
recomputes the update *from* the residuals, the frames are inside the trust
boundary whether the document admits it or not. Withdrawing it is the correct call
and I would not restore it.

### 8.3 Is variant 5 dead, or dead only for the shape I built?

Dead for every shape I could build **while the residual relation is present** — I
tried four routes and all four failed:

1. **Fabricating residuals with the relations intact** — impossible by construction;
   `g` is pinned to `y` through `AR`, and `y` is pinned to `iterate_next` by §4.3.6.
2. **A relation self-consistent but detached from the iterate** — if `g` does not
   reference `y` transitively it is constant across frames, so
   `residual_curr − residual_prev = 0` and §4.3.7's zero-denominator clause fires.
   The attack defeats itself.
3. **Exploiting `**` / `sqrt` / `round`** — see §8.4 for the one real (small) finding;
   none of them let me move the converged answer arbitrarily.
4. **Steering through the preamble** — the starting trials are free, but with the
   constants and relations fixed the secant converges to the same root regardless,
   so a different start changes only the update count, which §7.2 already treats as
   incidental. Correctly not an attack.

Alive **only** through the coverage gap of §8.1. That is a real distinction and I
want it on the record: §4.6's mechanism is sound; its applicability is optional.

### 8.4 Three further findings

**(a) `constants` is checked against nothing — and so is every other row of §3.8's
table.** Choosing `K` freely and running the node's own declared physics to
convergence gives a fully conforming trace with a different answer every time:
```
K=6.619 -> 1.380 m, 3 updates, PASS   <-- gold's own K
K=8.2   -> 1.521 m, 3 updates, PASS
K=4.9   -> 1.201 m, 4 updates, PASS
```
So the literal answer to your question 3 is: **no, §4.6 does not satisfy §3.8's
amended form.** But I am not filing that as a defect, and I want to be consistent
with round 4, where I explicitly declined to file `universe`, `precedences`,
`tolerance` and `max_elements` for exactly this reason. **`constants` is problem
data.** Nothing binds any node to its question — §10/D-039 records that — so a
node-local verifier cannot know that `K` is *this item's* `K`. That is the
comparator's job under §7, where a candidate is scored against a gold node
carrying the true constants, and a divergent constant is caught immediately.

The real finding is one level up: **§3.8's amended rule, applied honestly,
condemns its own table.** `item_measures`, `budget`, `precedences` and `constants`
are all declared once and checked against nothing — five of the nine rows. The
rule needs the distinction it currently lacks:

> *Problem data* (declared once, bound to the item by the comparator, not by the
> node) versus *derived data* (must be recomputed from problem data). §3.8 governs
> the second absolutely; for the first, "declared once" is the whole requirement,
> and the check lives in §7.

Add that distinction and a column to the table, and §3.8 becomes auditable instead
of self-contradicting. The severity gap is worth stating plainly: with a free
constant the trace is a **correct solve of a different problem**, and §7 catches it;
with a missing relation (§8.1) the trace is a correct solve of **no** problem, and
§7.2 grants it full process credit.

**(b) `round(sym, n)` restates a precision that `symbol_precision` already
declares.** The digit is a free literal inside the relation and nothing requires it
to match. Changing the residual relation from `round(AR, 3)` to `round(AR, 1)`,
with `symbol_precision["AR"]` left at 3, produces a conforming trace with a
different answer:
```
residual relation round(AR, 1), symbol_precision[AR]=3 -> 1.383 m, PASS
```
Only 0.003 m here, but it is a free knob on the answer, and it is a §3.8
restatement violation *inside the section that exists to enforce §3.8*. Fix:
make it `round(sym)` with no second argument, taking the digit from
`symbol_precision[sym]` — which also removes the question of what a non-integer or
expression-valued second argument means.

**(c) §4.6's 6.7% figure is not reproducible from the corpus it ships with.**
Recomputing `AR` from each frame's *stored* (rounded) `A` and `P` — i.e. the
bind-to-stored reading the section warns against — and comparing to the stored `AR`:
```
42 / 169 frames = 24.9%   (§4.6 claims 6.7% over 12,735 frames / 3,000 seeds)
```
**The claim's direction is confirmed and its number is not.** Binding to stored
values genuinely does reject gold traces, so the unrounded-substitution rule is
right, load-bearing, and well worth the paragraph §4.6 gives it — I implemented it
from that paragraph and got 80/80 first time. But a reader who checks 6.7% against
the shipped 40-seed corpus measures 24.9% and concludes the section is wrong. Either
cite the population explicitly ("6.7% over the 3,000-seed sample; 24.9% on the
committed 40") or re-measure. This is the same class as §4.5's headroom figure,
which *did* reproduce exactly (I checked it in round 2) — the difference is that
§4.5's number is derivable from the artefact and this one is not.

### 8.5 Is the evaluation rule implementable as written?

**Yes, and better than I expected — the ordering discipline is the strongest part
of §4.6.** "Names resolve to `constants`, to the frame's `iterate` role, or to an
**earlier** symbol of `frame_relations`" makes cycles and forward references
*structurally* impossible rather than merely forbidden. I verified all three are
rejected at parse time, before any evaluation:
```
forward reference (AR before A)   -> FAIL ['4.6']
self reference (A defined via A)  -> FAIL ['4.6']
undeclared name in a relation     -> FAIL ['4.6']
```
There is no recursion-depth question, because the ordering bounds the substitution
depth by the number of relations. I did not need to guess.

**Five things I did have to decide, none blocking, one sentence each:**

| # | Undefined | What I did |
|---|---|---|
| 1 | `sqrt` of a negative | trace failure (§4.3.9) |
| 2 | division by zero inside a relation | trace failure |
| 3 | `0 ** -1` and other `**` domain errors | trace failure |
| 4 | may `round`'s first argument be an *expression*? §4.6 writes `round(x, n)` but every example is a bare symbol | I allow expressions — `round((b + z*y)*y, 3)` parses and passes |
| 5 | must `round`'s second argument be an integer literal? does `round` use the node's `rounding` mode? | integer required (2.5 fails); mode taken from `rounding`, per §6's blanket rule, which §4.6 does not repeat |

Items 1–3 are the interesting ones: a verifier that let them raise instead of fail
would crash on a hostile candidate trace, and §6.0 makes candidate traces a
first-class mode. §4.6 should say they are trace failures.

### 8.6 Recommendation

**One blocking clause** (§8.1's coverage requirement) and **two one-liners**
(`round(sym)` without the duplicated digit; the 6.7% provenance). Then, not
blocking but worth doing before the phase closes because it is the thing that keeps
generating rounds: **give §3.8 the problem-data / derived-data distinction** and
re-run its audit against both node types. That audit — the one §3.8 prescribes and
which has still not been run — is exactly what would have caught §8.1 in ten
minutes, and it is what I would want in place before a third node type is added.

I said after round 4 that a fifth round was not needed, and you were right to
overrule me: I judged the schema as it stood, and adding a new expression language
reset that judgement. The mechanism you built is sound — I attacked it four ways
and it held — but it ships as optional, and optional is how variant 5 comes back.
Close the coverage gap and I think §4.6 is finished.

*Round 5 filed by Reviewer D2. Written to disk, not committed. New artefact:
`tests/trace_schema/reviewer_d2_verifier_14.py`.*
