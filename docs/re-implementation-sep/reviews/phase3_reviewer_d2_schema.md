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
