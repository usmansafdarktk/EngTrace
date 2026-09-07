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

G_CHECK = "**Answer:** The deflection is 18.0 mm, which is acceptable."

# ==========================================================================
# Round 2 -- the cases that broke the fixes for round 1.
#
# Both reviewers were asked to attack the MECHANISM the round-1 fixes were
# built on rather than the fixes, which is Phase 3's most expensive lesson.
# They did, and `commitment.py` v1 reproduced two of the defects it was written
# to close: a hedge window measured in characters had become a hedge window
# measured in clauses, and a longest-surface rule had become a last-clause rule
# that still resolved by last position *within* the clause.
#
# Four of these findings are CONVERGENT -- reached independently by two
# reviewers working in isolation, which R1.7 names as the strong signal.
# ==========================================================================

CASES += [

    # ---- E R2-F7: F2 relocated.  The round-1 case, one period -> one comma --
    ("E-R2-F7", "categorical", G_NOT,
     "## Final Answer\n**Answer:** The system is linear, unlike a nonlinear "
     "system, since both tests pass.", "MISMATCH",
     "**false accept**: a comma instead of a period defeated the round-1 fix"),
    ("E-R2-F7", "categorical", G_NOT,
     "**Answer:** The system is linear (a nonlinear system would fail additivity).",
     "MISMATCH",
     "**false accept**: a parenthetical gloss outranked the commitment"),

    # ---- E R2-F8: the neighbour rule was N=1 clause, not clause scope ------
    ("E-R2-F8", "categorical", G_LIN,
     "**Answer:** I am not sure. Let me redo the algebra. The system is linear.",
     "UNRESOLVED", "**false accept**: one intervening clause let the hedge escape"),
    ("E-R2-F8", "categorical", G_LIN,
     "**Answer:** I cannot really tell. The two tests are in step 4. The system "
     "is linear.", "UNRESOLVED", "**false accept**: same, at distance 2"),
    ("E-R2-F8", "categorical", G_LIN,
     "**Answer:** The system is linear. Step 4 shows the algebra. But I am not sure.",
     "UNRESOLVED", "**false accept**: same escape, trailing direction"),

    # ---- E R2-F9 / B R2-F1.1 (CONVERGENT): concessive != hypothetical ------
    ("E-R2-F9", "check", G_CHECK,
     "**Answer:** The deflection is 18.0 mm. Given that the limit is 25 mm, the "
     "design is acceptable.", "MATCH",
     "**over-rejection**: a fronted concessive asserts its matrix clause"),
    ("E-R2-F9", "check", G_CHECK,
     "**Answer:** The deflection is 18.0 mm, roughly 72% of the 25 mm limit, so "
     "the design is acceptable.", "MATCH",
     "**over-rejection**: `roughly` is a precision term, not an epistemic hedge"),
    ("E-R2-F9", "check", G_CHECK,
     "**Answer:** The deflection is 18.0 mm. Note that the limit is 25 mm, so it "
     "is acceptable.", "MATCH",
     "**over-rejection**: `note that` is factive -- its complement is asserted"),
    ("E-R2-F9", "check", G_CHECK,
     "**Answer:** The deflection is 18.0 mm, so the design is acceptable. The "
     "margin seems comfortable.", "MATCH",
     "**over-rejection**: a trailing remark about something else retracted it"),
    ("B-R2-F1", "categorical", G_LIN,
     "**Answer:** Given that additivity and homogeneity both hold, the system is "
     "linear.", "MATCH", "**over-rejection**: fronted concessive"),
    ("B-R2-F1", "categorical", G_LIN,
     "**Answer:** Although the equation contains a delay, the system is linear.",
     "MATCH", "**over-rejection**: fronted concessive"),
    ("B-R2-F1", "categorical", G_LIN,
     "**Answer:** Note that the system is linear.", "MATCH",
     "**over-rejection**: factive `note that`, a preface not a caveat"),
    ("B-R2-F1", "categorical", G_NOT,
     "**Answer:** The system is not linear. This should be clear from the square.",
     "MATCH", "**over-rejection**: confidence read as a hedge"),
    ("B-R2-F1", "categorical", G_LIN,
     "**Answer:** Imagine a scaled input a*x[n]; the output scales identically, "
     "so the system is linear.", "MATCH",
     "**over-rejection**: `imagine` opens the standard homogeneity proof"),
    ("B-R2-F1", "categorical[tuple]", G_MC_YY,
     "**Answer:** a) Memoryless: Yes b) Causal: Yes. Both look immediate.",
     "MATCH", "**over-rejection**: trailing remark about a different subject"),

    # ---- E R2-F11: the splitter cut inside decimals and abbreviations ------
    ("E-R2-F11", "categorical", G_LIN,
     "**Answer:** The gain is 2.5 and the system is linear.", "MATCH",
     "the splitter cut inside a decimal, a distance amplifier for R2-F8"),
    ("E-R2-F11", "categorical", G_LIN,
     "**Answer:** The system is linear i.e. additive and homogeneous.", "MATCH",
     "the splitter cut inside an abbreviation"),

    # ---- E R2-F12: _POST_NEG_RE stepped into a parenthetical ---------------
    ("E-R2-F12", "categorical", G_NOT,
     "**Answer:** The system is linear (the wrong answer would be nonlinear).",
     "MISMATCH", "a parenthetical *about* wrongness negated the label it follows"),
    ("E-R2-F12", "categorical", G_NOT,
     "**Answer:** Calling it linear is the wrong description here.", "MATCH",
     "control: post-label negation must still work"),

    # ---- controls: the round-1 behaviour these fixes must not undo ---------
    ("R2-ctl", "categorical", G_LIN,
     "**Answer:** If additivity holds, the system is linear.", "UNRESOLVED",
     "control: a HYPOTHETICAL still suspends its matrix (B round-1 F3)"),
    ("R2-ctl", "categorical", G_LIN,
     "**Answer:** Linear. It seems.", "UNRESOLVED",
     "control: a bare anaphoric comment still governs (B round-1 F2)"),
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
