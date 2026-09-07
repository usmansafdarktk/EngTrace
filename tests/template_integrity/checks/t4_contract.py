"""T4 - Output-contract conformance.

Every emitted solution must be mechanically segmentable:

* `**Step N:**` markers, strictly formatted, numbered contiguously from 1,
  with no duplicates and no restarts;
* exactly one terminal answer marker, from the approved set;
* a non-degenerate `(question, solution)` pair - no answer-less return path.

Every one of these was violated somewhere in the corpus before Phase 5. Three
electrical templates emitted `**Step 2: **` / `**Step 3: ...**` (colon or space
inside the bold), which a strict marker regex drops silently; five chemical
templates terminated with `**Final Answer**` / `**Final Answers:**`;
`template_levenspiel_plot_interpretation` restarted its numbering at
1,2,3,1,2,3,4,5,6 so the marker was not a unique id (fixed in Phase 2).

**The approved set is not the canonical one, deliberately.** ANSWER_MARKERS is
the set of markers this check can *recognise*; CANONICAL_ANSWER_MARKER is the
only one gold may *emit* (D5.3). Keeping the recognised set wide is what lets a
regression be reported as `non-canonical marker {'**Final Answer**': 25}`
rather than as the far less useful `no answer marker on 25 seeds`. Narrowing
the recognised set would gate the same defect with a worse diagnostic.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field

from ..core import ANSWER_MARKERS, CANONICAL_ANSWER_MARKER, Instance, STEP_RE

# Anything that looks like a step heading, however malformed.
_STEP_ANY = re.compile(r'\*\*\s*Step\s*(\d+)\s*[:\.]?', re.IGNORECASE)


@dataclass
class ContractResult:
    template_id: str
    instances: int = 0
    malformed_markers: Counter = field(default_factory=Counter)
    non_contiguous: list[int] = field(default_factory=list)
    duplicate_numbers: list[int] = field(default_factory=list)
    missing_answer: list[int] = field(default_factory=list)
    non_canonical_answer: Counter = field(default_factory=Counter)
    multiple_answers: list[int] = field(default_factory=list)
    empty_output: list[int] = field(default_factory=list)
    step_counts: Counter = field(default_factory=Counter)

    @property
    def passed(self) -> bool:
        # `non_canonical_answer` joined this list in Phase 5 (SPEC-CHANGE 14).
        # Before then T4 *printed* `non-canonical marker {'**Final Answer**': 25}`
        # for five templates and passed them anyway, so the corpus could be
        # reported 150/150 green with a known output-contract defect standing on
        # five items.  That is a severity gap, not a blind spot: the check already
        # saw the defect.  The remedy is this one line, not a second scanner --
        # an earlier draft of the Phase 5 brief proposed the scanner, which would
        # have duplicated T4 rather than fixing it.
        return not (self.malformed_markers or self.non_contiguous
                    or self.duplicate_numbers or self.missing_answer
                    or self.multiple_answers or self.empty_output
                    or self.non_canonical_answer)

    def summary(self) -> str:
        bits = []
        if self.malformed_markers:
            bits.append(f'malformed markers {dict(self.malformed_markers)}')
        if self.non_contiguous:
            bits.append(f'non-contiguous numbering on {len(self.non_contiguous)} seeds')
        if self.duplicate_numbers:
            bits.append(f'duplicate step numbers on {len(self.duplicate_numbers)} seeds')
        if self.missing_answer:
            bits.append(f'no answer marker on {len(self.missing_answer)} seeds')
        if self.multiple_answers:
            bits.append(f'multiple answer markers on {len(self.multiple_answers)} seeds')
        if self.empty_output:
            bits.append(f'degenerate output on {len(self.empty_output)} seeds')
        if self.non_canonical_answer:
            bits.append(f'non-canonical marker {dict(self.non_canonical_answer)}')
        return '; '.join(bits) or 'ok'


def check_instance(inst: Instance, res: ContractResult) -> None:
    res.instances += 1
    sol = inst.solution
    if not sol.strip() or not inst.question.strip():
        res.empty_output.append(inst.seed)
        return

    strict = [int(m.group(1)) for m in STEP_RE.finditer(sol)]
    loose = [int(m.group(1)) for m in _STEP_ANY.finditer(sol)]
    if len(loose) > len(strict):
        # A heading a human reads as a step that the strict regex cannot see.
        for m in _STEP_ANY.finditer(sol):
            frag = sol[m.start():m.start() + 24]
            if not STEP_RE.match(frag):
                res.malformed_markers[frag.split('\n')[0][:22]] += 1

    if strict:
        res.step_counts[len(strict)] += 1
        if len(set(strict)) != len(strict):
            res.duplicate_numbers.append(inst.seed)
        elif strict != list(range(1, len(strict) + 1)):
            res.non_contiguous.append(inst.seed)

    present = [mk for mk in ANSWER_MARKERS if mk in sol]
    if not present:
        res.missing_answer.append(inst.seed)
    else:
        if CANONICAL_ANSWER_MARKER not in sol:
            res.non_canonical_answer[present[0]] += 1
        # count only the canonical/first marker's occurrences
        if sol.count(present[0]) > 1:
            res.multiple_answers.append(inst.seed)
    if not strict and not loose:
        res.empty_output.append(inst.seed)


def run(instances: list[Instance], template_id: str) -> ContractResult:
    res = ContractResult(template_id=template_id)
    for inst in instances:
        if inst.ok:
            check_instance(inst, res)
        else:
            res.empty_output.append(inst.seed)
    return res
