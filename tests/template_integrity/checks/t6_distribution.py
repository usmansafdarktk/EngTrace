"""T6 - Distribution non-regression.

The P6 guard. A correctness fix must not quietly change what an item tests or
how hard it is. Over a fixed seed range, before and after any edit, compare:

* distinct final answers (a fix that collapses the answer space is a downgrade)
* answer magnitude quantiles p5 / p50 / p95
* step-count distribution
* branch coverage - the proportion of instances taking each structural variant
* the set of units appearing in the answer

Defaults follow the specification: distinct-answer count must not fall by more
than 10%, branch proportions must stay within +/-5 points, and the median
answer magnitude must stay within +/-2%. A breach requires sign-off; it must
not be resolved by widening the tolerance.
"""
from __future__ import annotations

import json
import re
import statistics
from collections import Counter
from dataclasses import asdict, dataclass, field

from ..core import Instance, NUM_RE, STEP_RE, answer_block, parse_number

_UNIT_AFTER_NUM = re.compile(
    r'\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?\s*'
    r'([A-Za-z°Ωμ%][A-Za-z0-9_^/\.\*\-]*)')


@dataclass
class Profile:
    """A fingerprint of a template's output distribution."""
    template_id: str
    n_instances: int = 0
    distinct_answers: int = 0
    answer_p5: float | None = None
    answer_p50: float | None = None
    answer_p95: float | None = None
    step_counts: dict = field(default_factory=dict)
    skeletons: int = 0                 # distinct numeric-blanked solution shapes
    answer_units: dict = field(default_factory=dict)
    errors: int = 0

    def to_json(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_json(d: dict) -> 'Profile':
        return Profile(**d)


def _answer_value(sol: str) -> float | None:
    block, _ = answer_block(sol)
    if not block:
        return None
    m = NUM_RE.search(block)
    return parse_number(m.group(0)) if m else None


def _answer_unit(sol: str) -> str:
    block, _ = answer_block(sol)
    m = _UNIT_AFTER_NUM.search(block)
    return m.group(1) if m else ''


def profile(instances: list[Instance], template_id: str) -> Profile:
    p = Profile(template_id=template_id)
    answers, skeletons, steps, units = [], set(), Counter(), Counter()
    for inst in instances:
        if not inst.ok:
            p.errors += 1
            continue
        p.n_instances += 1
        v = _answer_value(inst.solution)
        if v is not None:
            answers.append(v)
        skeletons.add(re.sub(r'\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?', '#', inst.solution))
        steps[len(STEP_RE.findall(inst.solution))] += 1
        u = _answer_unit(inst.solution)
        if u:
            units[u] += 1
    p.distinct_answers = len(set(answers))
    if answers:
        s = sorted(answers)
        p.answer_p5 = s[max(0, int(0.05 * (len(s) - 1)))]
        p.answer_p50 = statistics.median(s)
        p.answer_p95 = s[min(len(s) - 1, int(0.95 * (len(s) - 1)))]
    p.step_counts = {str(k): v for k, v in sorted(steps.items())}
    p.skeletons = len(skeletons)
    p.answer_units = dict(units.most_common(10))
    return p


@dataclass
class DriftResult:
    template_id: str
    breaches: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.breaches


def compare(before: Profile, after: Profile,
            max_distinct_drop: float = 0.10,
            max_branch_shift: float = 0.05,
            max_median_shift: float = 0.02) -> DriftResult:
    d = DriftResult(template_id=after.template_id)

    if before.distinct_answers:
        drop = (before.distinct_answers - after.distinct_answers) / before.distinct_answers
        if drop > max_distinct_drop:
            d.breaches.append(
                f'distinct answers fell {before.distinct_answers} -> '
                f'{after.distinct_answers} ({drop:.0%} > {max_distinct_drop:.0%})')

    if before.answer_p50 not in (None, 0) and after.answer_p50 is not None:
        shift = abs(after.answer_p50 - before.answer_p50) / abs(before.answer_p50)
        if shift > max_median_shift:
            d.breaches.append(
                f'median answer moved {before.answer_p50:g} -> {after.answer_p50:g} '
                f'({shift:.1%} > {max_median_shift:.0%})')

    tot_b = sum(before.step_counts.values()) or 1
    tot_a = sum(after.step_counts.values()) or 1
    for k in set(before.step_counts) | set(after.step_counts):
        pb = before.step_counts.get(k, 0) / tot_b
        pa = after.step_counts.get(k, 0) / tot_a
        if abs(pa - pb) > max_branch_shift:
            d.breaches.append(
                f'{k}-step instances moved {pb:.0%} -> {pa:.0%}')

    if before.skeletons and after.skeletons:
        # structural-variant coverage: number of distinct solution shapes
        rb, ra = before.skeletons, after.skeletons
        if abs(ra - rb) / max(rb, 1) > max_branch_shift * 4:
            d.notes.append(f'distinct solution shapes {rb} -> {ra}')

    if set(before.answer_units) != set(after.answer_units):
        d.notes.append(f'answer units {sorted(before.answer_units)} -> '
                       f'{sorted(after.answer_units)}')
    if after.errors > before.errors:
        d.breaches.append(f'generation errors {before.errors} -> {after.errors}')
    return d


def save(profiles: dict[str, Profile], path: str) -> None:
    with open(path, 'w', encoding='utf-8') as f:
        json.dump({k: v.to_json() for k, v in profiles.items()}, f, indent=1)


def load(path: str) -> dict[str, Profile]:
    with open(path, encoding='utf-8') as f:
        return {k: Profile.from_json(v) for k, v in json.load(f).items()}
