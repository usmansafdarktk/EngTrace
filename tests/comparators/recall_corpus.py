"""A recall corpus: correct answers, phrased awkwardly.

**The artefact Phase 4 did not have, named by Reviewer B (round 2, §5.2).**

D4.4 is 157 hand-written near misses, and almost all of them are answers that
*should be rejected*. It is a precision corpus wearing both hats. B's round-2
§4 made the consequence concrete:

    no new false reject among the 61 -- but 39 of them are the two
    classification items whose archived answers are terse verdicts with no
    subordinate clause and no adverb. The sample contains almost no prose, so
    it **cannot** detect over-rejection.

That is D-034's shape a third time, and it is why version 1 of
``commitment.py`` could refuse 25 of 38 credited answers while every gate stayed
green.

**The construction, and its honest limits.** Each of the 45 archived answers
that ``ground_truth.py`` labels correct is re-wrapped in prose frames drawn from
the constructions the commitment machinery reasons about: fronted concessives,
factives, trailing elaborations, derivation openers, and hedged-reasoning-then-
firm-conclusion. The **content** of every case is a real archived model answer
with a known label; the **frame** is synthetic. So this corpus tests whether a
real answer survives ordinary English packaging -- it does not establish that
models write this way. Reviewer E's R2-F10 stands: across 2,200 archived spans
the hedge probes fire once, so the machinery remains a *policy* validated
against constructed text, and that is stated rather than smoothed over.

Every case expects **MATCH**. A `MISMATCH` here is a comparator that inverted a
correct answer; an `UNRESOLVED` is one that refused it. Both are false rejects
and both cost recall, which is the half of the gate that `UNRESOLVED` cannot
escape.

Run: ``python -m tests.comparators.recall_corpus``
"""

from __future__ import annotations

import re
from collections import Counter

from .answer import compare_template
from .ground_truth import LABELS, load
from .normalize import answer_span, prepare
from .score import STEM_TO_TEMPLATE

#: Prose frames.  ``{a}`` is the archived answer, lower-cased where the frame
#: needs it.  Each names the construction it exercises and the finding that
#: made it worth exercising.
FRAMES: list[tuple[str, str, str]] = [
    ("plain", "{a}", "control: the archived answer unchanged"),
    ("fronted-concessive",
     "Although the algebra is fiddly, {a}",
     "B R2-F1.1 / E R2-F9a: a concessive asserts its matrix clause"),
    ("fronted-given",
     "Given that both defining conditions were checked above, {a}",
     "B R2-F1.1: `given that` is concessive, not hypothetical"),
    ("factive",
     "Note that {a}",
     "B R2-F1.1: `note that` is factive -- its complement is asserted"),
    ("trailing-elaboration",
     "{a} This follows from the two tests above.",
     "E R2-F9c: a trailing clause must not retract a commitment"),
    ("derivation-opener",
     "Consider a scaled input a*x[n] and a shifted input x[n-k]. {a}",
     "B R2-F1.3: `consider`/`imagine` open a proof, they do not hedge it"),
    ("hesitant-then-firm",
     "The intermediate algebra took two attempts. {a}",
     "B R2-F1.3: hesitant reasoning followed by a firm conclusion"),
    ("appositive-relative",
     "{a} which is what the definitions require.",
     "E R2-F9 / chk-07: a non-restrictive relative asserts its content"),
    ("post-quantified",
     "{a} The margin here is comfortable.",
     "E R2-F9c: a following remark about a different subject"),
]

#: **Negative controls, and the corpus was shipped without them.**
#:
#: Reviewer B (R3-F1) built these and they found a false accept immediately:
#: wrapping an answer in a hypothetical must make it non-committal, and **16 of
#: the 39 still MATCHed** -- exactly the enumerated ``a) … b) …`` answers, whose
#: frame the enumerator split discarded before the comparator ever read it.
#: 144 of the 351 positive cases were therefore further copies of the ``plain``
#: control, and the corpus was **41% inert** while reporting 100%.
#:
#: A recall corpus with no negative control cannot tell "the comparator credits
#: a correct answer" from "the comparator never saw the frame". That is the
#: same shape as D4.4 sampling the inside of the hedge list it certified
#: (B, round-1 F1), one level up, and I built it twice.
NEGATIVE_FRAMES: list[tuple[str, str, str]] = [
    ("neg-hypothetical",
     "If the additivity test holds, {a}",
     "a conditional states a rule, not a verdict"),
    ("neg-task-restatement",
     "We must determine whether {a}",
     "the item's own prompt wording is not an answer"),
    ("neg-explicit-hedge",
     "I am not sure, but perhaps {a}",
     "a speaker-level hedge governs the whole utterance"),
    ("neg-question",
     "{a} Or is it?",
     "a trailing question withdraws the commitment"),
]

#: Only the kinds `commitment.py` gates are worth framing this way -- it is the
#: module under test.  `sequence` and `symbolic` answers are formulae, and
#: wrapping a formula in a concessive tests the frame rather than the answer.
GATED = ("system_property_linearity", "system_properties_memory_causality")


def archived_answer(row: dict) -> str:
    """The model's committed answer span, as one line of prose."""
    span, _ = answer_span(row["model_reasoning"])
    text = prepare(span)
    text = re.sub(r"\s+", " ", text).strip()
    return text if text.endswith((".", "!")) else text + "."


def build() -> list[dict]:
    cases: list[dict] = []
    for stem in GATED:
        rows = load(stem)
        labels = LABELS[stem]
        for i, row in enumerate(rows):
            correct, _ = labels[i]
            if not correct:
                continue
            ans = archived_answer(row)
            if not ans or len(ans) > 400:
                continue
            for name, frame, why in FRAMES + NEGATIVE_FRAMES:
                body = frame.format(a=ans)
                cases.append({
                    "id": f"{stem[:12]}-{i:02d}-{name}",
                    "template": STEM_TO_TEMPLATE[stem],
                    "gold": row["gold_answer"],
                    "candidate": "## Final Answer\n**Answer:** " + body,
                    "frame": name,
                    "why": why,
                    "expect_match": not name.startswith("neg-"),
                })
    return cases


def main() -> int:
    cases = build()
    print("Recall corpus -- real archived correct answers in ordinary prose frames")
    print("=" * 92)
    print(f"{len(cases)} cases from {len(cases) // len(FRAMES)} archived correct answers "
          f"x {len(FRAMES)} frames\n")

    by_frame: dict[str, Counter] = {name: Counter()
                                    for name, _, _ in FRAMES + NEGATIVE_FRAMES}
    failures: list[dict] = []
    for c in cases:
        v = compare_template(c["template"], c["gold"], c["candidate"])
        by_frame[c["frame"]][v.outcome] += 1
        if v.is_match != c["expect_match"]:
            failures.append({**c, "outcome": v.outcome, "reason": v.reason})

    def _block(title: str, frames: list, want_match: bool) -> tuple[int, int]:
        print(f"\n{title}")
        print(f"  {'frame':22s} {'MATCH':>6s} {'MISMATCH':>9s} {'UNRESOLVED':>11s}  "
              f"{'correct':>8s}  exercises")
        print("  " + "-" * 88)
        tot = good = 0
        for name, _, why in frames:
            c = by_frame[name]
            n = sum(c.values())
            g = c["MATCH"] if want_match else n - c["MATCH"]
            tot, good = tot + n, good + g
            print(f"  {name:22s} {c['MATCH']:6d} {c['MISMATCH']:9d} {c['UNRESOLVED']:11d}  "
                  f"{g / n if n else 0:8.1%}  {why[:32]}")
        print("  " + "-" * 88)
        print(f"  {'subtotal':22s} {'':6s} {'':9s} {'':11s}  {good / tot:8.1%}  ({good}/{tot})")
        return good, tot

    pg, pt = _block("POSITIVE -- a correct answer must be credited (expect MATCH)",
                    FRAMES, True)
    ng, nt = _block("NEGATIVE -- a framed answer must NOT be credited (expect not-MATCH)",
                    NEGATIVE_FRAMES, False)
    print(f"\n{'TOTAL':24s} {(pg + ng) / (pt + nt):8.1%}  ({pg + ng}/{pt + nt})")

    if failures:
        print(f"\n{len(failures)} FALSE REJECTS -- a correct archived answer refused:")
        seen = set()
        for f in failures:
            key = (f["frame"], f["reason"][:60])
            if key in seen:
                continue
            seen.add(key)
            print(f"\n  [{f['id']}] {f['outcome']}")
            print(f"    {f['candidate'].splitlines()[-1][:150]}")
            print(f"    reason: {f['reason'][:130]}")
        if len(failures) > len(seen):
            print(f"\n  ({len(failures) - len(seen)} further failures share a reason above)")
    else:
        print("\nNo false rejects: every archived correct answer survives every frame.")

    print("\nLIMIT, stated: the content of every case is a real archived answer with a")
    print("known label; the FRAME is synthetic. This measures whether a correct answer")
    print("survives ordinary English packaging. It does NOT establish that models write")
    print("this way -- across 2,200 archived spans the hedge probes fire once (E, R2-F10),")
    print("so the commitment machinery remains a policy validated on constructed text.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
