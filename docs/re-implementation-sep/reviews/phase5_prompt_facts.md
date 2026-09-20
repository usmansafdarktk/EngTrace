# Phase 5 prompt — factual audit

**Artefact under review:** `docs/prompts/06-phase5-contract-hygiene-and-bindings.md`
**Frozen ref:** `bbccae510b8eac64f45c3e2b5e92a5be8b3f067b` (branch `redesign/phase5-scoping`)
**Scope:** facts only. The plan, the two-track split, the deliverable set and the review
protocol are a second reviewer's. Nothing was fixed; this is a report.
**Time box:** 35 minutes, honoured. §4 lists what fell outside it.

---

## 1. Verdict

**Factually safe to hand over with four corrections applied first.** Every number in the
T1–T7 table, every named template, and the whole of the D-058 defect story reproduce
exactly. But **three numbers the prompt presents as measured facts are not reproducible
from the repository** (`19,668` is mislabelled as "all 150"; `22,982`, `58` and `87` have
no recorded derivation), and **one claim about T4 is wrong in a way that will cost the
next session a check it does not need to build**. The prompt's own framing — "verify these
yourself before acting" — is the right instinct, and the next session *cannot* verify three
of them, because the definitions are not written down anywhere.

---

## 2. Claims verified

`PYTHONIOENCODING=utf-8` set for every command below. All run from the repo root at the
frozen ref.

### 2.1 Corpus health table

| Claim | As stated | As measured | Verdict |
|---|---:|---:|---|
| T1 printed-arithmetic closure failing | 29/150 | **29** | agree |
| T1 marginals | 97 | **97** | agree |
| T2 round-trip oracle | 0/150 | **0** | agree |
| T3 determinism | 0/150 | **0** | agree |
| T4 output contract | 3/150 | **3** | agree |
| T5 binding / rounding | 66/150 | **66** | agree |
| T6 distribution | 142/150 | **142** | agree |
| T7 invariant asserts | 83/150 | **83** | agree |

```
python -m tests.template_integrity.run --checks all      # 150 templates, 30s
```

The whole table is exact. This is the first phase brief in the series whose headline
corpus table needed no correction.

### 2.2 Named templates

| Claim | As measured | Verdict |
|---|---|---|
| T4 fails on exactly `euclidean_distance_binary`, `cd_dc_system_analysis`, `finite_convolution` | exactly those three | agree |
| `cd_dc_system_analysis` drops Step 3 in **100%** of instances | **25/25 seeds**; emitted marker is `**Step 3: ** Discrete-` — a strict `\*\*Step (\d+):\*\*` parser finds `['1','2','4']` | agree |
| …and Step 3 computes the answer | Step 3 is "Discrete-to-Continuous (D/C) Conversion", which produces the reconstructed signal | agree |
| 5 templates use `**Final Answer**` / `**Final Answers:**` | exactly the spec's five | agree |
| `levenspiel_plot_interpretation` is also in Phase 2's scope | spec L224 (Phase 2 table) and L891 both list it | agree |
| N1 false accept on `hagen_poiseuille_flowrate` | seeds 3 and 9 are both `Engine Oil (SAE 50)`, `0.001183` vs `0.002738`; `parse_number` returns `50` for both | agree, verbatim |
| N1 reads a temperature `541` on `gas_viscosity_kinetic_theory` | seed 1, `Argon at 541 K`, `parse_number → 541` | agree |
| N1 reads "the 4 of C₄H₁₀" | seed 0; `prepare()` folds `C₄H₁₀ → C4H10`, then `parse_number → 4` | agree |
| N2 `autocorrelation_rect_pulse` **128 of 132** | **128** false accepts of 132 ordered pairs | agree |
| N3 `euclidean_distance_binary` **132/132** | **132** of 132 | agree |
| N1 "3 scalar templates measured to false-accept in 12 instances" | **3**: `hagen_poiseuille_flowrate` (2 pairs), `gas_viscosity_kinetic_theory` (14), `null_to_null_bandwidth` (2) | agree |
| `parse_number` reads the **first** number | `kinds.py:73` — `_NUM_RE.search(text)`, first match | agree |
| 4 of 150 templates bound | `len(PHASE4_BINDINGS) == 4` | agree |
| `normalize.py::ANSWER_MARKERS` already accepts several variants, priority order load-bearing | 5 patterns, ordered, with an in-file comment on why the colon in the last is required | agree |

Reproduction for the pair sweep (note the **corrected** call signature — see §3.5):

```python
from tests.comparators.answer import compare_kind   # compare_kind(kind, gold, candidate)
KIND = {'scalar':'numeric','classification':'categorical','symbolic':'symbolic',
        'vector':'sequence','array':'sequence','multipart':'numeric'}
# for each template, 12 seeds, all i!=j: compare_kind(KIND[answer_type], spans[j], spans[i])
# truth = (spans[i] == spans[j])
```

### 2.3 Track B starting-point table

| Claim | As stated | As measured | Verdict |
|---|---:|---:|---|
| templates bound | 4 of 150 | **4** | agree |
| gold×gold pairs, "all 150 × 12 instances" | 19,668 | **19,668 — over 149 templates** | **number agrees, label diverges** |
| archive×gold pairs | 22,982 (11,384 scalar, 5,511 multipart) | **not reproducible**; nearest construction gives 23,292 with a different split | **diverge / unverifiable** |
| gold answers carrying a unit token | 87 of 150 | **not reproducible**; 82 or 68 depending on the quantifier | **unverifiable** |
| templates exposed to N1 | 58 of 150 | **not reproducible**; 64 or 56 depending on the quantifier | **unverifiable** |
| templates measured to N1-false-accept | 3 | **3** | agree |

### 2.4 History, documents and commands

| Claim | As measured | Verdict |
|---|---|---|
| Phase 4 merge commit is `0dabd27` | `0dabd27f…` is a two-parent merge, `master` points at it, ancestor of the frozen ref | agree |
| `phase4_summary.md` §8 / §9 / §12 exist and say what is cited | §8 "Errors I made", §9 "The reviews — four rounds each", §12 "What Phase 5 and 6 inherit" and it does state the stopping rule verbatim | agree |
| "nineteen errors" | §8 opens "Nineteen." | agree |
| "thirteen shared three shapes" | §8: "Three shapes account for thirteen of them" (6+4+3) | agree |
| "eight review rounds" | §9: four rounds each × two reviewers | agree |
| `phase4_comparators.md` §2 / §3 / §5 | §2 "Two rules that run through all six kinds", §3 "`multipart` is not a seventh kind", §5 "What this specification does **not** settle" | agree |
| Reviewer E report, "read F0" | `### F0 — CONFIRMED (blocking, structural)` present | agree |
| Reviewer E report §"Round 4 — Task 3" | **no string "Task 3" anywhere in the file** | **diverge** |
| `tests/template_integrity/phase4_instance_dump.py` exists, takes tree root as argv, refuses out-of-tree | exists; `root = abspath(argv[1])`, `if not resolved.startswith(root): raise SystemExit` | agree |
| the trap is "documented in `instance_dump.py`" | it is, in that file's docstring | agree |
| "it has now caught **four** people" | both dump scripts say **three** | **diverge** |
| spec Phase 5 has both tracks and D5.4–D5.11 | spec L386–L462; all eight deliverables present with the same wording | agree |
| SPEC-CHANGE 11, 12, 13 exist and are Phase 4's | revision log L842–L844; all three, all Phase 4 | agree |
| D-003, D-024, D-026, D-034, D-043, D-045, D-046, D-047, D-049, D-052, D-053, D-056, D-057, D-058 | all present; every heading matches the gloss the prompt gives it | agree |
| `sympy` absent from `requirements.txt` (D-053) | absent | agree |
| D5.8's kind counts (5 classification, 9 symbolic, 15 vector/array, 32 multipart, 89 scalar) | inventory: 5 / 9 / 8+7 / 32 / 89 = 150 | agree |
| `_froude_capped_slope` affects three templates in `uniform_flow.py` | three call sites (L109, L185, L391) | agree |

All eight verification commands run and exit 0:

```
python -m tests.template_integrity.run --checks all                 exit 0
python -m tests.trace_schema.audit_3_8                              exit 0, 80/80 clean
python -m tests.comparators.score                                   exit 0, both gates PASS
python -m tests.comparators.reviewer_battery                        exit 0
python -m tests.comparators.recall_corpus                           exit 0
python -m tests.comparators.derive_vocabulary                       exit 0
python -m tests.trace_schema.candidate_7 …/candidates.json          exit 0
ENGTRACE_HEDGE_POLICY=enforce  (battery + recall)                   exit 0 both
```

- "81 reviewer cases" — `reviewer_battery.CASES` has **81**. Agree.
- "507 cases" — the tool prints "507 cases from 39 archived correct answers × 9 positive + 4
  negative frames". Agree.
- `ENGTRACE_HEDGE_POLICY=enforce` is read at `commitment.py:98`
  (`os.environ.get("ENGTRACE_HEDGE_POLICY", "advisory")`); default is `advisory` as claimed.
  Agree.

### 2.5 D-058 discipline

The prompt quotes **no** number that D-058 disowns. D-058 disowns "68 `categorical` and 136
`sequence` false rejects, and most of its 560 false accepts" as harness artefacts of default
scoring; my sweep reproduces **560 total false accepts** and the prompt cites none of them.
It also carries D-058's 58-vs-3 correction forward correctly and repeats it as a standing
caution. **Clean.**

---

## 3. Findings

### F1 — CONFIRMED. "19,668 pairs, all 150 × 12 instances" is 149 templates, not 150

`150 × 12 × 11 = 19,800`. The measured figure is **19,668 = 149 × 132**. One template
contributes zero pairs: **`template_continuous_to_discrete_conversion`** raises
`AttributeError: 'tuple' object has no attribute 'free_symbols'` inside the `symbolic`
comparator on every pair.

```
python - <<'EOF'
# 12 seeds per template, all i!=j, compare_kind(kind, spans[j], spans[i]) in try/except
# -> templates scored: 149 of 150 ; total scored pairs: 19668 ; raising: 1
EOF
```

**Impact: high.** Track B's gate is *"zero false accepts on gold×gold across all 150"*. As
the instrument stands, one template is not in the denominator and cannot produce a false
accept — it produces an exception that a naive harness swallows. A session that reports
"19,668 pairs, zero false accepts, all 150 clean" would be stating something false, and this
is precisely the "green while measuring nothing" failure the prompt warns about twice. It is
also a fourth defect alongside N1/N2/N3 that the brief does not name: **the symbolic
comparator crashes on a real corpus template.**

**Correction:** say `19,668 pairs over 149 templates; `continuous_to_discrete_conversion`
raises and is unscored — closing it is part of D5.4.`

### F2 — CONFIRMED. "The other three defect classes … are invisible to T4" is false for the largest of them

| class | n | T4 verdict at the frozen ref |
|---|---:|---|
| malformed `**Step N:**` | 3 | `pass: false` — gated |
| `**Final Answer**` / `**Final Answers:**` | 5 | `pass: **true**`, summary `non-canonical marker {'**Final Answer**': 25}` |
| malformed complex | 2 | `pass: true`, summary `ok` — genuinely invisible |
| unreachable `elif` | 1 | `pass: true`, summary `ok` — genuinely invisible |

```
python -m tests.template_integrity.run --checks T4 --json out.json
# then read out.json['template_batch_moles_vs_conversion']['T4']
```

T4 already **detects and reports** the marker class on all five templates — it just does not
fail on it. Only 3 of the 11 are invisible, not 8.

**Impact: medium-high, and it cuts the work down.** The prompt tells the next session to
"extend T4, add a check, or state per-defect-class acceptance evidence" for eight templates.
For five of those eight the evidence already exists in T4's own JSON and the extension is a
severity change on an existing summary, not a new check. A session that believes "invisible"
will build a scanner that duplicates `checks/T4`. The prompt's derived line — *"'T4 passes
150/150' would be true with 8 of the 11 edits unmade"* — is still true as stated about the
*gate*, and should be reworded to say so: T4 sees five of them and does not gate them.

### F3 — CONFIRMED. `22,982` archive×gold pairs is not reproducible, and its per-kind split diverges materially

The archive is 2,200 records across `error_analysis_annotation/samples/*.jsonl` (confirmed).
Template identity is recoverable exactly as `ground_truth.load` does it —
`question_id.rsplit("__", 1)[0]` — giving **90 distinct stems**, of which **259 traces map to
no inventory template**. Pairing every mapped trace against 12 gold instances of its own
template:

| kind | my derivation | prompt / D-058 |
|---|---:|---:|
| `scalar` | 9,900 | **11,384** |
| `multipart` | 6,744 | **5,511** |
| total | **23,292** | **22,982** |

```
python - <<'EOF'
# stems = Counter(question_id.rsplit('__',1)[0]) over the 2,200 archive records
# join to template_inventory.csv on 'template_'+stem, multiply each kind's traces by 12
EOF
```

The totals differ by 310 and the scalar/multipart split differs by ~1,500 in *opposite*
directions, so this is not an off-by-one — the author used a different construction. Neither
D-058 nor the spec records what it was.

**Impact: high.** D5.5 is *"archive×gold cross-pairing … 22,982 real-text pairs available"*.
A session that builds the instrument and gets 23,292 has no way to tell whether it built the
wrong instrument or the brief is stale, and the brief's own instruction ("re-derive it") will
send it round that loop. **Escalating rather than assuming I am wrong:** two independent
constructions disagree and only one of them is written down.

Two facts the brief should carry that it does not: **the archive covers only 90 of the 150
templates**, and **259 of its 2,200 traces map to no inventory template**. D5.5's coverage is
therefore at most 60% of the corpus by construction — which matters for a deliverable framed
as validating bindings.

### F4 — CONFIRMED. `58 exposed` and `87 unit-carrying` have no recorded operational definition

Neither number is derivable without knowing what the author counted.

**"58 of 150 templates exposed to N1 (leading incidental number)."** The obvious operational
reading — the first number in the answer span is not the last — gives:

| predicate over the 12 spans | count |
|---|---:|
| holds on ≥1 of 12 | **64** |
| holds on all 12 | **56** |
| stated | **58** |

**"87 of 150 gold answers carry a unit token."** Matching `kinds.UNIT_SURFACES` (11 canonical
units) with word boundaries against the answer span:

| predicate | count |
|---|---:|
| holds on ≥1 of 12 | **82** |
| holds on all 12 | **68** |
| stated | **87** |

```
python - <<'EOF'
from tests.comparators.kinds import UNIT_SURFACES, _NUM_RE
from tests.comparators.normalize import answer_span, prepare
# 12 spans/template; count templates satisfying the predicate under each quantifier
EOF
```

58 sits between my two N1 quantifiers; 87 sits **above both** unit counts, which suggests the
author used a broader unit vocabulary than `UNIT_SURFACES` (only 11 canonical units are
defined, and the corpus emits `m^3/s`, `Pa·s`, `K`, `mol/min`, `%` and others that are not in
it).

**Impact: medium-high, and it is the exact hazard the prompt itself names.** D5.10 is
"declare units for the **87**" and D5.6's candidate 2 says "D5.10 gives you 87 of these". A
session that re-derives and gets 82 will either chase a phantom five or quietly adopt its own
number and report a deliverable as complete against a different denominator. **These are the
"true but unverifiable by the next session" class:** they may well be right, but nothing in
the repository lets anyone confirm it. Both should ship with the predicate written out, or
with a script.

### F5 — CONFIRMED. The cited section `Round 4 — Task 3` does not exist

`phase4_reviewer_e_comparator.md` contains no occurrence of "Task 3", "Task 2" or "Task 1".
Round 4's headings are `# Round 4 — is the mechanism converging?` → `## 1. Verdict`,
`## 2. Independent re-derivation (round 4)`, `## 3. RECOMMENDATION — read this before the
findings`, `## 4. Falsification attempts…`, `## 5. Further probing…`.

```
grep -n "Task 3" docs/re-implementation-sep/reviews/phase4_reviewer_e_comparator.md   # no match
```

**Impact: medium.** The prompt makes this a hard precondition — "read it before you write a
line of comparator code" — and it is the single most load-bearing document reference in the
brief, since §3 RECOMMENDATION is where E argues the surface is unbounded and the stopping
rule comes from. Almost certainly `Round 4 §3` was meant. As written the instruction cannot
be followed literally, and a session that greps for it and finds nothing may skip it.

### F6 — CONFIRMED. "`origin/master` is not fetched in this working copy" is false

`origin/master` exists at `65a7d58`. `master` (`0dabd27`) is **87 ahead, 0 behind**.

```
git rev-parse origin/master           # 65a7d588…
git rev-list --left-right --count origin/master...master   # 0  87
```

**Impact: low** — the prompt immediately says "check the actual divergence yourself before
believing any number for it — including this one", which is the correct hedge and does its
job. Reported for completeness and because the surrounding claim, "`master` has never been
published", is a different and possibly still-true claim resting on the false one. The
remote is `https://github.com/usmansafdarktk/EngTrace.git` and a published `gh-pages` branch
exists, so "never published" deserves its own check before anyone repeats it.

### F7 — CONFIRMED (minor). "it has now caught four people"

Both `phase3_instance_dump.py` and `phase4_instance_dump.py` say **three**, at L13. The
prompt says four. Either the docstrings are stale or the prompt inflated it; the repository
says three.

**Impact: low.** The instruction it decorates ("run each tree as its own process") is
correct and the mechanism is exactly as described.

### F8 — PLAUSIBLE (minor). "a bare `reviews` entry in `.gitignore` once silently untracked all eight prior reports"

`.gitignore:14–18` records the incident but gives no count: *"a bare `reviews` pattern here
silently untracked every phase review report from Phase 0 onward"*. The reviews directory now
holds **13** reports. "Eight" is not sourced anywhere I could find and does not match any
obvious snapshot of the directory. The *caution* is real and correctly stated; only the
number is unsupported.

**Impact: negligible** — the instruction ("confirm your reports are actually committed") is
sound regardless.

### F9 — PLAUSIBLE (minor, spec not prompt). The spec's Track A review line says "the 8 touched"

`template_redesign_spec.md:408` — *"Reviewer A re-runs T4 across all 150 templates, not just
the 8 touched"* — in a section whose own table totals **11 distinct templates**. The prompt
consistently says eleven and is right; the spec disagrees with itself. Flagged because the
prompt sends the next session to read that section first, and because Track A's brief will be
written from it.

---

## 4. Claims I could not check

| Claim | Why not |
|---|---|
| The 2 malformed-complex (`time_to_phasor`, `phasor_addition`) and 1 unreachable-`elif` (`decimation_aliasing_analysis`) defects are *real* | I confirmed T4 passes them silently (F2) but did not read the emitted text or the `elif` guard. The spec asserts them; T4 cannot corroborate. Track A's Reviewer A should not treat these as measured. |
| "Five reviewers independently raised the AI-Tribunal criticism" | In `docs/EngTrace_Rebuttal_Jul2026.txt` / `_ARR_May__EngTrace.txt`; not opened inside the time box. |
| "150 templates across 47 files, 30 per branch" | Inventory has 150 rows and the branch column is populated; I did not count files or per-branch totals. |
| Phase-3-specific history: "Phase 3's specification took six rounds", "Phase 3 broke the frozen-ref rule three times", "Phase 0 found two checks green while measuring nothing" | These live in commit bodies and phase summaries; checking them properly means reading Phase 0/3 summaries, and `git log` is forbidden to me. The Phase 4 analogues all checked out, so I rate these plausible. |
| D5.11's three residuals (R4-11, RB4-1, R4-12) say what the spec says | Spec text confirmed; the underlying claims (e.g. "0 of 16 archived traces reach its branch") not re-measured. |
| That `levenspiel_plot_interpretation`'s **Phase 2** edit has or has not landed | The prompt correctly instructs the session to check this; it is an instruction, not a claim, so I left it. Note the inventory `notes` column still describes the unseeded-`np.random` bug, which may be stale rather than current. |
| Whether T6 "will move" for the marker edits | Forward-looking; unfalsifiable before the edits. |

**One methodological note.** The reproduction recipe I was handed calls
`compare_kind(cand, gold, kind)`. The real signature is `compare_kind(kind, gold, candidate)`
(`answer.py:214`), and the wrong order does not raise a `TypeError` — it raises `KeyError`
on **all 150 templates** and silently reports **zero** false accepts. Anyone re-deriving
D-058 from that snippet gets a clean sweep and concludes N1/N2/N3 are fixed. The prompt does
not contain the snippet, but D5.4/D5.5 will be built from the same call, and this is worth a
line in the brief.

---

## 5. What would make this prompt harder to get wrong

1. **Never state a count without its predicate.** `58`, `87` and `22,982` are the three
   claims that failed, and all three failed the same way: the number is recorded, the
   question it answers is not. Every one of my attempted re-derivations was defensible and
   none reproduced. The fix is one clause each — "87 = templates where every one of 12 gold
   spans contains a token from `<vocabulary>`" — or, better, a two-line script committed
   beside the number. This is D-034's own rule ("a number a script recomputes from the
   artefact is evidence") applied to the *brief* rather than to the code, and the brief is
   currently the one artefact in the chain exempt from it.

2. **A pair count must state its denominator in templates, not just in pairs.** `19,668`
   reads as complete and is not; `19,668 over 149 of 150` would have surfaced the crashing
   symbolic comparator in the brief instead of in the next session's third hour. Any
   cross-pairing figure should ship as `pairs / templates-scored / templates-skipped`, with
   skips named. Extend this to the gate: "zero false accepts across all 150" needs "and 150
   templates produced pairs" beside it, or it is satisfiable by crashing.

3. **Distinguish "not gated by X" from "invisible to X".** F2 is a two-word error with an
   hour of cost attached, because the remedy differs: invisible needs a new instrument,
   ungated needs a severity change. The brief already draws this distinction well elsewhere
   (exposure vs. defect, and it says so twice); it is the same distinction and deserves the
   same explicit treatment.

4. **Cite document sections by their literal heading text.** Four of the five section
   references resolved; `Round 4 — Task 3` did not, and it was the one marked as a hard
   precondition. A reference the reader can `grep` for verbatim either resolves or fails
   loudly. Where a heading is long, quote the first five words of it.

5. **Numbers about the repository's own history decay silently.** "four people", "eight
   prior reports", "not fetched" — three small claims, three divergences, all in the
   advisory prose rather than the measured tables. The measured tables were flawless. That
   asymmetry is the lesson: the parts of the brief that came from a run were right, and the
   parts that came from memory were not. Prose claims about tooling and history should be
   either sourced to a file and line, or dropped.

6. **Say what the archive does *not* cover.** The brief describes the 2,200-trace archive
   only in terms of what it yields. It covers 90 of 150 templates and 259 of its traces map
   to no template at all. D5.5's coverage ceiling is a structural fact about the deliverable,
   and it belongs in the brief rather than in the session's discovery log — especially in a
   phase whose central lesson (F0, R2-F10) is that an unstated corpus gap makes a gate
   vacuous.
