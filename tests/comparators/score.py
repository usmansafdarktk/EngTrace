"""Score the comparators against the ground-truth labels, and against D4.4.

Reported quantities, defined so that ``UNRESOLVED`` cannot be used to dodge
either number (Decision 3, D-049):

``precision``  MATCH-and-truly-correct / MATCH.
               A false accept is the numerator's enemy.  UNRESOLVED never
               enters the denominator, so a comparator cannot buy precision by
               declining to decide -- it can only buy it by not accepting
               wrong answers.
``recall``     MATCH-and-truly-correct / truly-correct.
               A false reject **and** an UNRESOLVED on a correct answer both
               cost recall.  This is the half of the gate that stops
               UNRESOLVED from being free.
``decided``    (MATCH + MISMATCH) / all.  The fraction the comparator ruled on.
               Reported alongside, never traded against the other two.

Run: ``python -m tests.comparators.score``
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass

from .answer import compare_kind, compare_template
from .ground_truth import LABELS, TEMPLATES, load

STEM_TO_TEMPLATE = {
    "signal_operations": "template_signal_operations",
    "system_properties_memory_causality": "template_system_properties_memory_causality",
    "system_property_linearity": "template_system_property_linearity",
    "incompressible_continuity": "template_incompressible_continuity",
}

ADVERSARIAL = os.path.join(os.path.dirname(__file__), "adversarial.json")


@dataclass
class Tally:
    n: int = 0
    truly_correct: int = 0
    accepted: int = 0
    accepted_correct: int = 0
    decided: int = 0
    unresolved_on_correct: int = 0
    false_accepts: list[str] = None       # type: ignore[assignment]
    false_rejects: list[str] = None       # type: ignore[assignment]

    def __post_init__(self) -> None:
        self.false_accepts = self.false_accepts or []
        self.false_rejects = self.false_rejects or []

    def add(self, correct: bool, outcome: str, tag: str) -> None:
        self.n += 1
        self.truly_correct += correct
        if outcome in ("MATCH", "MISMATCH"):
            self.decided += 1
        if outcome == "MATCH":
            self.accepted += 1
            if correct:
                self.accepted_correct += 1
            else:
                self.false_accepts.append(tag)
        else:
            if correct:
                self.false_rejects.append(f"{tag} [{outcome}]")
                if outcome == "UNRESOLVED":
                    self.unresolved_on_correct += 1

    @property
    def precision(self) -> float | None:
        return None if not self.accepted else self.accepted_correct / self.accepted

    @property
    def recall(self) -> float | None:
        return None if not self.truly_correct else self.accepted_correct / self.truly_correct

    @property
    def decided_rate(self) -> float:
        return self.decided / self.n if self.n else 0.0

    def merge(self, other: "Tally") -> None:
        self.n += other.n
        self.truly_correct += other.truly_correct
        self.accepted += other.accepted
        self.accepted_correct += other.accepted_correct
        self.decided += other.decided
        self.unresolved_on_correct += other.unresolved_on_correct
        self.false_accepts += other.false_accepts
        self.false_rejects += other.false_rejects


def _pct(x: float | None) -> str:
    return "  n/a " if x is None else f"{x:6.1%}"


def score_archive() -> dict[str, Tally]:
    out: dict[str, Tally] = {}
    for stem in TEMPLATES:
        rows = load(stem)
        labels = LABELS[stem]
        t = Tally()
        for i, row in enumerate(rows):
            correct, _ = labels[i]
            v = compare_template(STEM_TO_TEMPLATE[stem], row["gold_answer"], row["model_reasoning"])
            t.add(correct, v.outcome, f"{stem}[{i}] {row['model']}")
        out[stem] = t
    return out


def score_adversarial() -> dict[str, Tally]:
    if not os.path.exists(ADVERSARIAL):
        return {}
    with open(ADVERSARIAL, encoding="utf-8") as fh:
        cases = json.load(fh)
    out: dict[str, Tally] = {}
    for c in cases:
        if c["kind"] == "narrative":
            continue  # scored by its own conformance rule, below
        t = out.setdefault(c["kind"], Tally())
        v = compare_kind(c["kind"], c["gold"], c["candidate"], **(c.get("options") or {}))
        t.add(bool(c["correct"]), v.outcome, f"{c['id']}: {c['note']}")
    return out


def score_narrative() -> tuple[int, int, list[str]]:
    """``narrative`` has no precision or recall; it has a conformance rule.

    **Every** narrative comparison must be UNRESOLVED -- including one where
    the candidate restates gold verbatim.  A single MATCH here would mean the
    kind had started deciding prose, which is the AI Tribunal returning through
    the back door (D-051).
    """
    if not os.path.exists(ADVERSARIAL):
        return 0, 0, []
    with open(ADVERSARIAL, encoding="utf-8") as fh:
        cases = [c for c in json.load(fh) if c["kind"] == "narrative"]
    bad = []
    for c in cases:
        v = compare_kind("narrative", c["gold"], c["candidate"])
        if v.outcome != "UNRESOLVED":
            bad.append(f"{c['id']} -> {v.outcome} ({c['note']})")
    return len(cases) - len(bad), len(cases), bad


def _report(title: str, tallies: dict[str, Tally]) -> Tally:
    print(title)
    print("-" * 96)
    print(f"{'group':38s} {'n':>4s} {'true+':>6s} {'acc':>5s} "
          f"{'precision':>10s} {'recall':>8s} {'decided':>8s}")
    total = Tally()
    for name, t in tallies.items():
        total.merge(t)
        print(f"{name:38s} {t.n:4d} {t.truly_correct:6d} {t.accepted:5d} "
              f"{_pct(t.precision):>10s} {_pct(t.recall):>8s} {_pct(t.decided_rate):>8s}")
    print(f"{'TOTAL':38s} {total.n:4d} {total.truly_correct:6d} {total.accepted:5d} "
          f"{_pct(total.precision):>10s} {_pct(total.recall):>8s} {_pct(total.decided_rate):>8s}")
    if total.false_accepts:
        print("\n  FALSE ACCEPTS (a wrong answer credited):")
        for f in total.false_accepts:
            print(f"    - {f}")
    else:
        print("\n  false accepts: none")
    if total.false_rejects:
        print("  FALSE REJECTS (a correct answer not credited):")
        for f in total.false_rejects:
            print(f"    - {f}")
    else:
        print("  false rejects: none")
    print()
    return total


def main() -> None:
    arch = _report("D4.2 evidence base -- 61 real archived traces", score_archive())
    adv_tallies = score_adversarial()
    if adv_tallies:
        adv = _report("D4.4 adversarial set -- hand-written near misses", adv_tallies)
    else:
        adv = None
        print("D4.4 adversarial set: adversarial.json not present\n")

    ok, n, bad = score_narrative()
    if n:
        print(f"narrative conformance: {ok}/{n} returned UNRESOLVED as required")
        for b in bad:
            print(f"    VIOLATION {b}")
        print()

    print("Gate (spec S4.5): >=95% precision AND >=95% recall on D4.4;")
    print("                  no false accept on real archived traces.")
    for name, t in (("archive", arch), ("adversarial", adv)):
        if t is None:
            continue
        p, r = t.precision or 0.0, t.recall or 0.0
        ok = p >= 0.95 and r >= 0.95
        print(f"  {name:12s} precision {p:6.1%}  recall {r:6.1%}  "
              f"{'PASS' if ok else 'FAIL'}   false accepts: {len(t.false_accepts)}")


if __name__ == "__main__":
    main()
