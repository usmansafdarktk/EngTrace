"""Every case Reviewer E and Reviewer B used to break the comparators.

**Kept as a standing check, not folded into `adversarial.json`.** The two
corpora answer different questions and merging them would lose that: D4.4 is the
implementer's set and this is the reviewers'. Reviewer B's F1 is precisely that
D4.4 could only sample the inside of the list it certified, so a corpus written
by the same hand that wrote the rules cannot be the whole gate. This file is the
part of the evidence that was written by someone trying to break the thing.

Each case carries the finding it came from and the outcome the reviewer argued
for. A case that still fails is reported, not deleted.

Run: ``python -m tests.comparators.reviewer_battery``
"""

from __future__ import annotations

from .answer import compare_kind, compare_template

G_LIN = "## Final Answer\n**Answer:** The system is **linear**."
G_NOT = "## Final Answer\n**Answer:** The system is **not linear**."
G_SEQ_NOORIG = "**Answer:**\nThe resulting sequence is y[n] = {-1, 5, -4, 3}"
G_MC_YY = "**Answer:**\na) Memoryless: **Yes**\nb) Causal: **Yes**"

#: ``(finding, kind, gold, candidate, expected_outcome, note)``
#: ``expected`` is what the reviewer argued the answer should be.
CASES: list[tuple[str, str, str, str, str, str]] = [

    # ---- Reviewer B, F1: unlisted hedges were accepted (12 of 15) ----------
    *[("B-F1", "categorical", G_LIN, f"## Final Answer\n{c}", "UNRESOLVED", note)
      for c, note in [
          ("I suspect the system is linear.", "verb of opinion, unlisted"),
          ("Apparently the system is linear.", "evidential, morphological neighbour of a listed word"),
          ("The system should be linear.", "epistemic modal"),
          ("The system looks linear.", "copula of appearance"),
          ("Tentatively, the system is linear.", "evidential"),
          ("Perhaps the system is linear.", "evidential"),
          ("Plausibly the system is linear.", "evidential"),
          ("Seemingly the system is linear.", "evidential, neighbour of listed 'seems'"),
          ("My guess is that the system is linear.", "verb of opinion"),
          ("I would say the system is linear.", "modal + verb of opinion"),
          ("Roughly speaking, the system is linear.", "evidential"),
          ("I am not 100% sure, but the system is linear.", "stated uncertainty"),
          ("Linear?", "interrogative, not an assertion"),
          ("The system appears to be linear.", "the spec's own named case; passed before too"),
          ("The system is most likely linear.", "evidential; passed before too"),
      ]],

    # ---- Reviewer B, F2: trailing hedges, all on LISTED vocabulary ---------
    *[("B-F2", "categorical", G_LIN, f"## Final Answer\n{c}", "UNRESOLVED", note)
      for c, note in [
          ("The system is linear, I think.", "trailing, same clause"),
          ("The system is linear, but that is only probably right.", "trailing, same clause"),
          ("Linear. It seems.", "trailing hedge in its own label-free clause"),
          ("The system is linear; I cannot determine this with confidence.", "trailing refusal clause"),
          ("Linear (unclear).", "parenthetical qualifier"),
          ("The system is linear. Hard to say, though.", "trailing clause"),
      ]],
    ("B-F2", "categorical[tuple]", G_MC_YY,
     "## Final Answer\n**Answer:** a) Memoryless, b) Causal, I think", "UNRESOLVED",
     "trailing hedge on a label tuple"),

    # ---- Reviewer B, F3: non-answers were accepted -------------------------
    *[("B-F3", "categorical", G_LIN, f"## Final Answer\n{c}", "UNRESOLVED", note)
      for c, note in [
          ("If additivity holds, the system is linear.", "conditional: states a rule, not a verdict"),
          ("Assume the system is linear.", "assumption: the opposite of a conclusion"),
          ("We must determine whether the system is linear.", "**the item's own prompt wording**"),
          ("I cannot complete the test; the system is linear.", "explicit refusal, then a guess"),
          ("The system is linear (to be verified).", "self-flagged as unverified"),
      ]],

    # ---- Reviewer E F2 / Reviewer B F4: the same defect, both directions ---
    ("E-F2", "categorical", G_NOT,
     "## Final Answer\n**Answer:** The system is linear. Unlike a nonlinear system, "
     "it satisfies both tests.", "MISMATCH",
     "**false accept**: a contrastive mention outranked the committed answer"),
    ("E-F2", "categorical", G_LIN,
     "## Final Answer\n**Answer:** Both tests pass, so the system is linear "
     "(it is not nonlinear).", "MATCH",
     "the mirror false reject of the same defect"),
    ("B-F4", "categorical", G_NOT,
     "## Final Answer\nThe system is nonlinear -- no, wait, it is linear.", "MISMATCH",
     "**false accept**: self-correction to the wrong answer was credited"),
    ("B-F4", "categorical", G_NOT,
     "## Final Answer\nNonlinear? No. The system is linear.", "MISMATCH",
     "**false accept**: a rhetorical question outranked the commitment"),
    ("B-F4", "categorical", G_LIN,
     "## Final Answer\nThe system is nonlinear. Actually, the system is linear.", "MATCH",
     "false reject: self-correction to the RIGHT answer was penalised"),
    ("B-F4", "categorical", G_LIN,
     "## Final Answer\nWait, nonlinear. Final answer: linear.", "MATCH",
     "the re-marker case: must agree with the line above, not depend on the marker"),

    # ---- Reviewer E, F3: gold omits the origin, wrong origins accepted -----
    ("E-F3", "sequence", G_SEQ_NOORIG,
     "**Answer:** y[n] = {*-1*, 5, -4, 3}", "UNRESOLVED",
     "**false accept**: asterisk pins n=0 at the wrong element (truth n=1..4)"),
    ("E-F3", "sequence", G_SEQ_NOORIG,
     "**Answer:** y[n] = {-1, 5, -4, 3} for n = -3, -2, -1, 0", "UNRESOLVED",
     "**false accept**: index list places the support four indices off"),
    ("E-F3", "sequence", G_SEQ_NOORIG,
     "**Answer:** y[0]=-1, y[1]=5, y[2]=-4, y[3]=3", "UNRESOLVED",
     "**false accept**: per-element listing, off by one"),
    ("E-F3", "sequence", G_SEQ_NOORIG,
     "**Answer:** y[n] = {-1, 5, -4, 3}", "MATCH",
     "control: no origin either side, values agree -- must still MATCH"),
    ("E-F3", "sequence", G_SEQ_NOORIG,
     "**Answer:** y[n] = {-1, 5, 3, -4}", "MISMATCH",
     "control: values differ, decidable regardless of origin"),

    # ---- Reviewer E, F4: answer_span peeled into a caveat ------------------
    ("E-F4", "categorical", G_NOT,
     "## Final Answer\n**Answer:** The system is linear.\n\n"
     "Note: Answer: not linear if the offset were nonzero.", "MISMATCH",
     "**false accept**: the caveat became the answer"),
    # LABEL CORRECTED (mine, not a comparator defect).  First written as
    # MISMATCH, which contradicts the case's own note: a bare restatement is
    # exactly what peel-to-last is FOR, so peeling to "linear" against a
    # `linear` gold is a MATCH.  The third label of mine that a check has
    # caught -- see phase4_summary.md §8.
    ("E-F4", "categorical", G_LIN,
     "## Final Answer\n**Answer:** The system is not linear.\n\n"
     "Answer: linear\n", "MATCH",
     "control: a bare restatement is NOT a qualification and must still be "
     "peeled to, or the E-F4 fix would break self-correction"),

    # ---- Reviewer E, F5: neither/nor scope -------------------------------
    ("E-F5", "categorical[tuple]", G_MC_YY,
     "## Final Answer\n**Answer:** The system is neither memoryless, because the "
     "output depends on past inputs as shown in step one, nor causal.", "MISMATCH",
     "**inversion**: a clause between neither and nor read (No,No) as (Yes,Yes)"),
    ("E-F5", "categorical[tuple]",
     "**Answer:**\na) Memoryless: **No**\nb) Causal: **No**",
     "## Final Answer\n**Answer:** The system is neither memoryless, because the "
     "output depends on past inputs as shown in step one, nor causal.", "MATCH",
     "the same candidate against the gold it actually matches"),
    ("E-F5", "categorical[tuple]",
     "**Answer:**\na) Memoryless: **No**\nb) Causal: **No**",
     "## Final Answer\n**Answer:** The system is neither memoryless nor causal.", "MATCH",
     "control: the archived short form must still work"),
]

TEMPLATE_FOR = {
    "categorical": "template_system_property_linearity",
    "categorical[tuple]": "template_system_properties_memory_causality",
    "sequence": "template_signal_operations",
}


def main() -> int:
    print("Reviewer battery -- every case E and B used to break the comparators")
    print("=" * 92)
    by_finding: dict[str, list[bool]] = {}
    failures = []
    for finding, kind, gold, cand, expected, note in CASES:
        v = compare_kind(kind, gold, cand)
        ok = v.outcome == expected
        by_finding.setdefault(finding, []).append(ok)
        if not ok:
            failures.append((finding, kind, cand, expected, v.outcome, note, v.reason))
    for finding, results in sorted(by_finding.items()):
        n, good = len(results), sum(results)
        flag = "ok" if good == n else f"{n - good} STILL FAILING"
        print(f"  {finding:6s} {good:2d}/{n:2d}  {flag}")
    total = sum(len(v) for v in by_finding.values())
    passed = sum(sum(v) for v in by_finding.values())
    print(f"  {'TOTAL':6s} {passed:2d}/{total:2d}")

    if failures:
        print("\nSTILL FAILING -- reported, not deleted:")
        for finding, kind, cand, exp, got, note, reason in failures:
            print(f"\n  [{finding}] {kind}: expected {exp}, got {got}")
            print(f"    candidate: {cand.splitlines()[-1][:90]}")
            print(f"    note:      {note}")
            if reason:
                print(f"    reason:    {reason[:110]}")
    else:
        print("\nEvery case both reviewers used now behaves as they argued it should.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
