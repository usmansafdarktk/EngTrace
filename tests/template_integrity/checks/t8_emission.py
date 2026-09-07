"""T8 - emission hygiene: the output-contract defects T4 cannot see.

T4 checks that a solution is *segmentable*: step markers, one answer marker, a
non-degenerate pair. Two output-contract defect classes are invisible to it,
because they are properties of the emitted *text* rather than of its structure:

* **doubled sign** - a hard-coded ``+`` in a format string followed by an
  interpolated signed value, printing ``876*t + -86.8 deg`` or
  ``30.5 + j-51.22``. Neither is a number in any notation.
* **degenerate product** - a zero coefficient printed as a product,
  ``omega_a = 0*pi``, where the template has a branch that would print ``0`` and
  a guard order that stops it running.

**Why this file exists rather than a phase-numbered report module.** Phase 5
delivered both detectors inside ``phase5_contract_scan``, with a planted-defect
self-test, and the reviewer's finding (F4) was that none of it was in a standing
gate: not in ``ALL_CHECKS``, not imported by ``run.py``, and the repository has
no CI. *Evidence is not a gate.* Six months on, a phase-numbered module runs only
when somebody remembers its name. The detectors are defined once, in
``phase5_contract_scan``, and imported here - so the report and the gate cannot
drift apart, and the scan stays the place where the reasoning lives.

**Scope, and why it is not uniform.** The corpus is clean on the ``j`` form and
on degenerate products, so both are gated corpus-wide. The wider doubled-sign
class is still emitted by **14 templates outside Track A's scope** (D-061,
assigned to Phase 6), so gating it corpus-wide today would fail the corpus for
work this phase is not chartered to do. It is gated on Track A's eleven and
reported as a census everywhere else. A gate reporting ``0`` for a class present
on 14 templates would be worse than one reporting the 14.

**Seeds.** T8 generates its own instances at 400 seeds rather than reusing the
suite's 25. The degenerate-product class fires on **1.95%** of
``decimation_aliasing_analysis`` instances, and 25 seeds resolves nothing rarer
than ~12% (D-024/D-026). A 25-seed T8 would report this phase's own headline
defect as absent - it did, on the first run of the scan.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..phase5_contract_scan import (
    CENSUS,
    GATED_ON_TRACK_A,
    TRACK_A,
    ScanResult,
    scan_solution,
)

#: The two classes T8 owns. `step_marker` and `answer_marker` stay with T4,
#: which sees them; duplicating them here would be a second scanner for a
#: defect that already has a gate.
T8_CLASSES = ('complex_sign', 'degenerate_product')

#: Instances per template. See the module docstring - this is a gate property.
T8_SEEDS = 400


@dataclass
class EmissionResult:
    template_id: str
    instances: int = 0
    errors: list[str] = field(default_factory=list)
    scan: ScanResult | None = None

    @property
    def gated(self) -> tuple[str, ...]:
        in_a = self.template_id in TRACK_A
        return T8_CLASSES + (GATED_ON_TRACK_A if in_a else ())

    @property
    def passed(self) -> bool:
        if self.errors or self.scan is None:
            return False
        return not any(self.scan.hits[c] for c in self.gated)

    def summary(self) -> str:
        if self.errors:
            return f'{len(self.errors)} generation errors'
        bits = []
        for c in self.gated:
            if self.scan.hits[c]:
                bits.append(f'{c} {dict(self.scan.hits[c])}')
        for c in CENSUS:
            if self.scan.hits[c]:
                bits.append(f'[census] {c} {dict(self.scan.hits[c])}')
        return '; '.join(bits) or 'ok'


def run(ref, seeds: int = T8_SEEDS) -> EmissionResult:
    """Generate ``seeds`` instances of ``ref`` and scan question and solution.

    T8 does not reuse ``run.py``'s instance list: it needs a different (larger)
    seed count, and it needs the **question** as well as the solution. Nothing
    in T1-T7 reads ``inst.question`` except to check that it is non-empty, and
    Track A changed 2,560 question strings (Reviewer A, F2).
    """
    from ..core import generate  # noqa: PLC0415 - avoids a package import cycle

    res = EmissionResult(template_id=ref.template_id)
    scan = ScanResult(template_id=ref.template_id,
                      in_track_a=ref.template_id in TRACK_A)
    for seed in range(seeds):
        inst = generate(ref, seed, capture=False)
        if not inst.ok:
            res.errors.append(f'seed {seed}: {inst.error}')
            continue
        scan.instances += 1
        scan_solution(inst.solution, scan, inst.question)
    res.instances = scan.instances
    res.scan = scan
    return res
