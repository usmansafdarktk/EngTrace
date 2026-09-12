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

from ..core import (ANSWER_MARKERS, CANONICAL_ANSWER_MARKER, Instance, STEP_RE,
                    answer_block)

# Anything that looks like a step heading, however malformed.
_STEP_ANY = re.compile(r'\*\*\s*Step\s*(\d+)\s*[:\.]?', re.IGNORECASE)

#: D6.8 - the SHAPE of the answer span, not merely the presence of its marker.
#:
#: The span is everything after the marker, and it is the ONLY part the grader
#: compares against a model's answer. D5.3 asserted the marker and never the
#: span, so a 363-character `**Note:**` commentary block sat inside one and the
#: check passed it (Reviewer A, F7) - the grader was comparing the model against
#: the answer PLUS a paragraph of prose.
#:
#: Two terms, because a gate with one term is passed by the defect it misses
#: (SPEC-CHANGE 18):
#:
#: * no foreign section marker inside the span. Measured across 150 templates x
#:   15 instances: ZERO occurrences of any of these. So it is asserted at zero
#:   tolerance because the corpus is actually clean, not fitted to whatever the
#:   corpus happened to contain.
#: * a length ceiling. The longest legitimate span is 199 characters
#:   (`autocorrelation_rect_pulse`; p50 64, p95 165). 300 leaves ~50% headroom
#:   over the measured maximum while still catching any intrusion of 101
#:   characters or more - the historical one was 363.
#:
#: A second `**Answer` is deliberately NOT listed here: `multiple_answers`
#: already reports it, and listing it too would report one defect twice.
MAX_ANSWER_SPAN = 300
_FOREIGN_IN_SPAN = ('**Note', '**Step', '##', '**Given', '**Find', '**Formulae')


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
    long_answer_span: list[int] = field(default_factory=list)
    foreign_in_span: Counter = field(default_factory=Counter)
    max_span: int = 0

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
        # `long_answer_span` and `foreign_in_span` joined in Phase 6 (D6.8).
        return not (self.malformed_markers or self.non_contiguous
                    or self.duplicate_numbers or self.missing_answer
                    or self.multiple_answers or self.empty_output
                    or self.non_canonical_answer
                    or self.long_answer_span or self.foreign_in_span)

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
        if self.long_answer_span:
            bits.append(f'answer span over {MAX_ANSWER_SPAN} chars on '
                        f'{len(self.long_answer_span)} seeds (max {self.max_span})')
        if self.foreign_in_span:
            bits.append(f'foreign section marker inside answer span '
                        f'{dict(self.foreign_in_span)}')
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

    # A marker is counted if it appears at all, NOT only when the canonical one
    # is missing (Reviewer A, F6).  `**Final Answer**` contains `Final Answer`
    # and `## Final Answer` does not overlap either, so the three are counted
    # independently; the containment that does exist is handled by scoring the
    # canonical marker's own count separately below.
    present = [mk for mk in ANSWER_MARKERS if mk in sol]
    if not present:
        res.missing_answer.append(inst.seed)
    else:
        for mk in present:
            if mk != CANONICAL_ANSWER_MARKER:
                res.non_canonical_answer[mk] += 1
        # More than one answer marker of ANY spelling is a contract violation:
        # the segmenter has to choose, and gold does not get to make it guess.
        # Counting only `present[0]`'s occurrences missed a second marker with a
        # different spelling entirely, which is F6's mechanism.
        if sol.count(CANONICAL_ANSWER_MARKER) > 1 or len(present) > 1:
            res.multiple_answers.append(inst.seed)

        # D6.8 - the span, not just the marker.  Measured on the stripped span
        # so a trailing newline is not counted as content.
        block, marker = answer_block(sol)
        span = block[len(marker):].strip() if marker else ''
        res.max_span = max(res.max_span, len(span))
        if len(span) > MAX_ANSWER_SPAN:
            res.long_answer_span.append(inst.seed)
        for pat in _FOREIGN_IN_SPAN:
            if pat in span:
                res.foreign_in_span[pat] += 1
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
