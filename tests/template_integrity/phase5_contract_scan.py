"""D5.2 - corpus-wide output-contract scan for Phase 5 Track A's four defect classes.

    python -m tests.template_integrity.phase5_contract_scan
    python -m tests.template_integrity.phase5_contract_scan --seeds 200
    python -m tests.template_integrity.phase5_contract_scan --json out.json

**Why this exists, and what it is not.**

Track A fixes four defect classes.  T4 sees two of them and only *fails* on one:

===================  ==========================  ==========================
class                T4 today                    remedy
===================  ==========================  ==========================
malformed ``Step``   FAILS                       already gated
non-canonical        REPORTS, then PASSES        **severity** change in T4
complex sign         invisible                   **new detector** (here)
degenerate product   invisible                   **new detector** (here)
===================  ==========================  ==========================

The distinction decides the remedy.  *Ungated* means the check sees the defect
and does not fail on it; the fix is one line of severity, not a new scanner.
*Invisible* means the check cannot see it at all; only then is a new instrument
warranted.  Reading the first group as the second builds a duplicate of T4, and
an earlier draft of the Phase 5 brief did exactly that.

So this module is deliberately **not** a re-implementation of T4.  Classes 1 and
2 are re-derived here anyway -- independently, with this file's own regexes, not
by calling ``t4_contract`` -- because a scan that reports "zero remaining
violations" across all four classes has to be able to see all four.  Two
independent readings of the same property is a cross-check; delegating to T4
would make this file's agreement with T4 vacuous.

**The two new detectors, and why each is written the way it is.**

``complex_sign``  A rectangular complex number printed as ``x + jy`` where the
    template interpolates ``y`` unconditionally after a hard-coded ``+``.  When
    ``y`` is negative the emitted text is ``30.5 + j-51.22``, which is not a
    number in any notation.  The detector looks for a sign, ``j``, and a second
    sign -- the *doubled* sign is the defect, not the ``j``.

``degenerate_product``  A zero-valued symbolic quantity rendered as a literal
    product: ``omega_a = 0*pi``.  The template has a branch that would print
    ``0``, placed after a branch that already caught the case, so it never
    runs.  The detector looks for a zero coefficient multiplying a symbol.
    Ordinary arithmetic lines like ``2*pi*0`` are not matched: the zero has to
    be the *left* operand, which is where a formatter puts a coefficient.

    **It is scoped to the answer span, and that scoping is measured rather than
    assumed.**  Written over the whole solution the detector fires on two
    templates outside Track A -- ``standing_wave_formation`` prints
    ``x = (0 * pi) / 4.56 = 0.0`` and ``impulse_response_from_lccde`` prints
    ``-3*delta[1] - 0*delta[0]``.  Both are *substitution displays*: a
    coefficient that happens to be zero, shown being substituted, and evaluated
    on the same line.  Neither is a degenerate rendering of an answer, and
    gating on them would be the harness telling a template not to show its
    working.  Over 60 seeds each, both have **60 and 25 whole-solution hits and
    zero answer-span hits**, so the split is not a hairline -- and the
    derivation-side occurrences are still counted, as a non-gating census, so
    that scoping this detector does not make them invisible.

Both new detectors were measured against the whole corpus before being adopted,
and both now have a false-positive rate of zero on the templates outside Track
A's scope -- which is the only evidence that makes them safe to gate on.  That
measurement is reproduced every time this module runs, and is reported as
``outside Track A`` rather than asserted here in prose (D-034).

**On the seed count, which is a gate property and not a convenience.**  Class 4
fires on **39 of 2,000** instances of ``decimation_aliasing_analysis`` -- 1.95%.
A 25-seed scan resolves nothing rarer than ~12% (D-024/D-026, ~3/N) and reports
this file's own headline defect as absent.  It did, on the first run.  The
default is therefore 400 seeds, which resolves ~0.75%, and the run prints the
smallest rate it can resolve beside the result so a later reader can check the
number against the claim rather than trusting the default.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass, field

from .core import Instance, discover, generate

# --------------------------------------------------------------------------
# The four detectors.  Each is independent of `t4_contract`, deliberately.
# --------------------------------------------------------------------------

#: Class 1 - a line a human reads as a step heading, however malformed.
#:
#: **Line-oriented, not marker-oriented** (Reviewer A, F5).  Anchored on the
#: bold marker it required the leading asterisks, so a heading that lost its
#: bold *entirely* was invisible to the "however malformed" detector -- and T4
#: then compares the surviving strict markers for contiguity, which succeeds
#: whenever the lost heading is the LAST one (1,2,3 is contiguous whether or
#: not a step 4 exists).  The final step is where the answer is computed,
#: which is why cd_dc_system_analysis was the urgent one.
_STEP_ANY = re.compile(r'^[ \t]*\*{0,2}[ \t]*Step[ \t]*\d+[ \t]*[:.][^\n]{0,60}',
                       re.MULTILINE)
#: ... and the one well-formed shape.  `**Step 12:**` and nothing else.
_STEP_OK = re.compile(r'\*\*Step (\d+):\*\*')

#: Class 2 - every answer marker shape the corpus is known to emit, INCLUDING
#: the ##-heading form (Reviewer A, F6).  Restricted to the bold shape it
#: could not see `## Final Answer` at all, so a solution carrying that
#: *beside* the canonical marker passed both this scan and T4 -- while on the
#: candidate side `normalize.ANSWER_MARKERS` ranks the heading at **priority
#: 0**, above the bold marker, and would extract from it.  That is exactly the
#: gold/candidate disagreement D5.3 exists to prevent, passing the check that
#: exists to prevent it.
_ANSWER_ANY = re.compile(
    r'\*\*\s*(?:Final\s+)?Answers?\s*:?\s*\*\*'
    r'|^\#{1,6}\s*(?:Final\s+)?Answers?\s*:?\s*$',
    re.MULTILINE)
CANONICAL_ANSWER_MARKER = '**Answer:**'

#: Class 3 - a DOUBLED SIGN: `+ j-51.22`, and equally `876*t + -86.8 deg`.
#:
#: The first version required a literal `j` and so gated only half of the
#: class D-061 defines (Reviewer A, F1) -- the half that did *not* change
#: 2,560 question strings.  `_signed_term`'s own docstring names
#: `876*t + -86.8 deg` as the first of the two things it was written to fix,
#: and the detector could not see it.
#:
#: **Gated on Track A's eleven; census corpus-wide.**  14 templates outside
#: Track A still emit this and are assigned to Phase 6 (D-061), so gating it
#: corpus-wide today would fail the corpus for work this phase is not
#: chartered to do.  Reporting `0 templates` for a class present on 14 would
#: be worse than reporting the 14.
_DOUBLED_SIGN = re.compile(
    r'(?<![\w.\-])[-+\u2212][ \t]*(?:j[ \t]*)?'
    r'[-+\u2212][ \t]*(?=[.\d]|j[.\d])')
#: The `j` sub-shape, kept separate so the corpus-wide claim stays honest:
#: the `j` form is gone corpus-wide; the wider class is not.
_COMPLEX_SIGN = re.compile(r'[-+]\s*j\s*[-+]\s*\d')

#: Class 4 - a zero coefficient printed as a product: `0*pi`, `0 * x`, and
#: `0.0*pi` (Reviewer A, F3).  `round(x, n)` of a small float is the ordinary
#: way a template arrives at a zero coefficient, so the float form is the
#: LIKELIER regression shape, not the rarer one.
#:
#: Anchored on the left operand, so `2*pi*0` and `k = 0` do not match, and the
#: lookbehind rejects a *word* character as well as a digit -- without that,
#: `q0 * B`, `p_0 * rho`, `P0 * a^3` and `d_0 * epsilon` all match, and six of
#: the detector's original eight hits were a variable's trailing subscript.
_DEGENERATE_PRODUCT = re.compile(
    r'(?<![\w.])0(?:\.0+)?\s*\*\s*'
    r'(?=pi\b|[A-Za-z]_?\w*\b)')

#: Gated corpus-wide.  The corpus is clean on all four.
CLASSES = ('step_marker', 'answer_marker', 'complex_sign', 'degenerate_product')

#: Gated on Track A's eleven templates ONLY, and reported as a census
#: everywhere else (Reviewer A, F1; D-061).  The wider doubled-sign class is
#: present on 14 templates this phase is not chartered to edit, so gating it
#: corpus-wide would fail the corpus for Phase 6's work.  Gating it on the
#: eleven is what makes Track A's own class-3 fix -- including the 2,560
#: question strings it changed -- defended by something.
GATED_ON_TRACK_A = ('doubled_sign',)
#: Reported, never gated.
#:  - `degenerate_product_derivation`: the same pattern in a derivation step,
#:    where a zero coefficient shown being substituted is legible working.
#:  - `doubled_sign_corpus`: the wider class outside Track A's eleven, so that
#:    scoping its gate does not also delete the observation.
CENSUS = ('degenerate_product_derivation', 'doubled_sign_corpus')

#: Every class, for iteration.
ALL_CLASSES = CLASSES + GATED_ON_TRACK_A + CENSUS

#: Every answer marker the corpus emits, longest first so `**Final Answers:**`
#: is found before `**Answer**` could match inside it.
_ANSWER_MARKER_LITERALS = (
    '**Final Answers:**', '**Final Answer**', '**Answer:**', '**Answer**',
)


def answer_span(sol: str) -> str:
    """The region of a solution that states the answer, or ``''``.

    The *last* occurrence wins, as in ``normalize.answer_span``: a template that
    restates its answer commits to the final statement.
    """
    cut = max((sol.rfind(mk) for mk in _ANSWER_MARKER_LITERALS), default=-1)
    return sol[cut:] if cut >= 0 else ''


@dataclass
class ScanResult:
    template_id: str
    instances: int = 0
    errors: list[str] = field(default_factory=list)
    #: class -> Counter of the offending text fragment -> instance count
    hits: dict[str, Counter] = field(
        default_factory=lambda: {c: Counter() for c in ALL_CLASSES})
    #: Track A's templates are gated on `GATED_ON_TRACK_A` as well.
    in_track_a: bool = False

    @property
    def clean(self) -> bool:
        """Gating classes only.  A census hit is an observation, not a failure."""
        gated = CLASSES + (GATED_ON_TRACK_A if self.in_track_a else ())
        return not self.errors and not any(self.hits[c] for c in gated)

    def summary(self) -> str:
        bits = [f'{len(self.errors)} generation errors'] if self.errors else []
        for c in ALL_CLASSES:
            if self.hits[c]:
                bits.append(f'{c} {dict(self.hits[c])}')
        return '; '.join(bits) or 'clean'


def scan_solution(sol: str, res: ScanResult, question: str = '') -> None:
    """Apply the four detectors to one emitted solution."""
    for m in _STEP_ANY.finditer(sol):
        frag = m.group(0)
        if not _STEP_OK.match(frag):
            res.hits['step_marker'][frag[:26]] += 1

    markers = [m.group(0) for m in _ANSWER_ANY.finditer(sol)]
    for mk in set(markers):
        if mk != CANONICAL_ANSWER_MARKER:
            res.hits['answer_marker'][mk] += 1
    if not markers:
        res.hits['answer_marker']['<none>'] += 1
    elif markers.count(CANONICAL_ANSWER_MARKER) > 1:
        res.hits['answer_marker']['<duplicate **Answer:**>'] += 1

    # Classes 3 and 4 are properties of the emitted ITEM, so they are read
    # over the question as well as the solution (Reviewer A, F2).  Track A's
    # class-3 fix changed 2,560 question strings and nothing read them.
    #
    # Classes 1 and 2 are deliberately NOT read over the question, and that is
    # a decision rather than an oversight: a question has no steps and states
    # no answer, so a step marker or an answer marker inside one is a
    # different defect from the one this scan is named for.
    emitted = sol + chr(10) + chr(10) + question if question else sol

    for m in _COMPLEX_SIGN.finditer(emitted):
        res.hits['complex_sign'][m.group(0).strip()] += 1
    for m in _DOUBLED_SIGN.finditer(emitted):
        key = emitted[m.start():m.start() + 6].strip()
        res.hits['doubled_sign' if res.in_track_a else 'doubled_sign_corpus'][key] += 1

    span = answer_span(sol)
    for m in _DEGENERATE_PRODUCT.finditer(span):
        res.hits['degenerate_product'][span[m.start():m.start() + 8].strip()] += 1
    n_all = len(_DEGENERATE_PRODUCT.findall(emitted))
    n_span = len(_DEGENERATE_PRODUCT.findall(span))
    if n_all > n_span:
        m = next(_DEGENERATE_PRODUCT.finditer(emitted))
        res.hits['degenerate_product_derivation'][
            emitted[m.start():m.start() + 8].strip()] += 1


def scan_template(ref, seeds: int) -> ScanResult:
    res = ScanResult(template_id=ref.template_id,
                     in_track_a=ref.template_id in TRACK_A)
    for seed in range(seeds):
        inst = generate(ref, seed, capture=False)
        if not inst.ok:
            res.errors.append(f'seed {seed}: {inst.error}')
            continue
        res.instances += 1
        scan_solution(inst.solution, res, inst.question)
    return res


#: The eleven templates Track A edits.  Named so the scan can report its own
#: false-positive rate on the 139 it does not touch, which is the evidence that
#: makes the two new detectors safe to gate on.
TRACK_A = (
    'template_cd_dc_system_analysis',
    'template_euclidean_distance_binary',
    'template_finite_convolution',
    'template_batch_moles_vs_conversion',
    'template_flow_system_molar_flow_rates',
    'template_gas_phase_concentration',
    'template_limiting_reactant',
    'template_levenspiel_plot_interpretation',
    'template_time_to_phasor',
    'template_phasor_addition',
    'template_decimation_aliasing_analysis',
)


# --------------------------------------------------------------------------
# D5.3 - the ANSWER_MARKERS agreement item.
#
# Track A decided that gold emits exactly one answer marker (`**Answer:**`).
# `tests/comparators/normalize.py::ANSWER_MARKERS` is the *candidate*-side list
# and stays wide, because models emit `## Final Answer`, `The answer is:` and
# the rest.  The two are not the same list and must not be: gold is a corpus we
# control, model output is not.
#
# **What "agree" means, stated as a predicate so it can be checked rather than
# asserted**: for every gold solution in the corpus, the candidate-side parser
# must recover the canonical marker and nothing else -- same marker, and an
# answer span with no marker debris glued to its front.
#
# This is not a formality.  Before Track A's edit the priority-2 pattern
# `Final\s+Answer:?` matched *inside* `**Final Answer**`, so the recovered
# marker was `Final Answer` and the span began `**\nAfter reaching...`; on
# `**Final Answers:**` it began `s:**\na) CSTR Volume...`.  Measured over the
# eleven Track A templates x 200 seeds: 1,000 of 2,200 spans carried debris
# before, 0 after.  A comparator reading those spans was reading `s:**` as part
# of the answer.
#
# Reviewer A owns checking this; this function is the instrument, and A is asked
# to treat its result as a claim under test rather than as tooling.
# --------------------------------------------------------------------------

def check_marker_agreement(seeds: int = 8) -> int:
    """Every gold answer span must be recovered by the candidate-side parser."""
    from tests.comparators.normalize import ANSWER_MARKERS, answer_span  # noqa: PLC0415

    print('  candidate-side ANSWER_MARKERS, in priority order:')
    for i, m in enumerate(ANSWER_MARKERS):
        print(f'    {i}: {m}')
    print(f'  gold-side canonical marker: {CANONICAL_ANSWER_MARKER}')
    print()

    found = Counter()
    bad: list[str] = []
    for ref in discover():
        for seed in range(seeds):
            inst = generate(ref, seed, capture=False)
            if not inst.ok:
                bad.append(f'{ref.template_id} seed {seed}: {inst.error}')
                continue
            span, marker = answer_span(inst.solution)
            found[marker] += 1
            if marker != CANONICAL_ANSWER_MARKER:
                bad.append(f'{ref.template_id} seed {seed}: recovered {marker!r}, '
                           f'span starts {span[:30]!r}')
            elif not span.strip():
                bad.append(f'{ref.template_id} seed {seed}: empty answer span')
    total = sum(found.values())
    for m, n in found.most_common():
        print(f'    {m!r:24s} {n}/{total}')
    if bad:
        print(f'  !! {len(bad)} disagreements')
        for line in bad[:15]:
            print(f'     {line}')
        return 1
    print(f'  AGREE: all {total} gold spans recover {CANONICAL_ANSWER_MARKER!r}')
    return 0


# --------------------------------------------------------------------------
# Planted-defect self-test.
#
# Phase 0's rule: a check is worth what it detects, and the only evidence for
# that is a defect planted on purpose and caught.  This phase *raised T4's
# severity* and *added two detectors*, and a green corpus proves neither -- the
# corpus is green because the eleven templates were fixed, and it would be
# equally green if the severity change had been a no-op.  So each of the four
# classes is planted into a real emitted solution and both instruments are
# required to see it.
#
# Planting into real generated text rather than a hand-written fixture is
# deliberate: a fixture would be this file's own answer key (D-034), and a
# detector tuned to a fixture is tuned to nothing.
# --------------------------------------------------------------------------

#: ``(class, find, replace)``.  The *site* of each plant is part of the test.
#:
#: ``degenerate_product`` is planted **after the answer marker**, because that is
#: the only place it is a defect.  The first version of this test planted it at
#: the first ``=`` in the solution -- a derivation step -- and the detector
#: correctly refused to fire, which read as a miss and was not one.  Keeping the
#: plant site explicit is what makes the answer-span scoping testable in both
#: directions rather than merely asserted in the docstring.
#: ``(class, find, replace, where)`` -- ``where`` is 'sol' or 'q'.
#:
#: **At least two plants per detector, of materially different surface form**
#: (Reviewer A, |S|5.9).  The first version had one plant per class, and the
#: plant was written by the same hand as the regex -- so every plant was a
#: shape its regex already matched, and the suite passed while four of the
#: detectors were narrower than the class they are named after.  Reviewer A
#: found all four by writing detectors from the spec's defect table instead.
#: The pairs below are the shapes that review produced:
#:
#:   step marker        malformed bold  /  bold lost entirely
#:   answer marker      marker swapped  /  a second marker ADDED beside it
#:   doubled sign       `+ j-5`         /  `+ -5` with no `j` at all
#:   degenerate product `0*pi`          /  `0.0*pi`
#:
#: and one plant per class is placed in the QUESTION rather than the solution,
#: because nothing read the question until F2.
PLANTS = (
    ('step_marker', '**Step 2:**', '**Step 2: **', 'sol'),
    ('step_marker', '**Step 2:**', 'Step 2:', 'sol'),
    ('answer_marker', '**Answer:**', '**Final Answer**', 'sol'),
    ('answer_marker', '**Answer:**', '## Final Answer\n**Answer:**', 'sol'),
    ('complex_sign', '= ', '= 30.5 + j-51.22 ', 'sol'),
    ('doubled_sign', '= ', '= 30.5 + j-51.22 ', 'sol'),
    ('doubled_sign', '= ', '= cos(876*t + -86.8 deg) ', 'sol'),
    ('doubled_sign', 'Calculate', 'cos(876*t + -86.8 deg). Calculate', 'q'),
    ('degenerate_product', '**Answer:**', '**Answer:** omega_a = 0*pi and', 'sol'),
    ('degenerate_product', '**Answer:**', '**Answer:** omega_a = 0.0*pi and', 'sol'),
    # D6.8, and its two shapes are the two TERMS of the rule, not two spellings
    # of one: an over-long span carrying no marker at all, and a foreign section
    # marker inside a span of ordinary length.  Either term alone is passed by
    # the defect the other catches (SPEC-CHANGE 18).
    ('answer_span', '**Answer:**', '**Answer:** ' + 'padding ' * 45, 'sol'),
    ('answer_span', '**Answer:**',
     '**Answer:**\n**Note:** commentary that belongs above the marker.', 'sol'),
)

#: Classes T4 owns and `scan_solution` deliberately does NOT detect.
#:
#: Without this the self-test cannot express a T4-only detector: `caught_scan`
#: is required for every class, so an answer-span plant would be scored as a
#: MISS because the scanner has no such detector -- a harness limitation read as
#: a detector failure.  The alternative was to add an answer-span detector to
#: `scan_solution` as well, which is the duplicate scanner this module's own
#: history warns against: Phase 5 nearly added one rather than fixing T4 with a
#: single line.  Making the expectation explicit in BOTH directions is the
#: cheaper and more honest fix.
SCAN_BLIND = ('answer_span',)

#: The negative half of the same test: planted where it is NOT a defect, the
#: answer-span detector must stay silent.  Without this the scoping is untested.
#: ``(class, find, replace, instrument)`` -- ``instrument`` is 'scan' or 't4',
#: because the negative half has to be checked by whichever instrument owns the
#: class.  D6.8's negative plant is the important one: a ``**Note:**`` block in a
#: DERIVATION STEP is perfectly legitimate, and only becomes a defect once it
#: sits after the answer marker.  Without this, the span scoping is asserted in
#: a docstring rather than tested, which is exactly what let D5.3 pass a span
#: carrying a whole extra section.
NEGATIVE_PLANTS = (
    ('degenerate_product', '**Step 2:**', '**Step 2:** h = 0*delta[0] and', 'scan'),
    ('answer_span', '**Step 2:**',
     '**Step 2:**\n**Note:** a note here is not a defect.', 't4'),
)


def selftest() -> int:
    """Plant each defect class in real output; require both instruments to catch it."""
    from .checks import t4_contract  # noqa: PLC0415

    refs = {r.template_id: r for r in discover()}
    # A template that is clean on all four classes, so a planted hit is
    # unambiguously the plant.  Its own T4 must pass before anything is planted.
    host = refs['template_bpsk_energy_basis']
    inst = generate(host, 0, capture=False)
    assert inst.ok, inst.error

    base = ScanResult(template_id='<clean-host>', in_track_a=True)
    scan_solution(inst.solution, base, inst.question)
    base_t4 = t4_contract.run([inst], host.template_id)
    ok = base.clean and base_t4.passed
    print(f'  host clean before planting: scan={base.clean} T4={base_t4.passed}')
    if not ok:
        print(f'  !! host is not clean: {base.summary()} / {base_t4.summary()}')
        return 1

    for cls, needle, repl, where in PLANTS:
        src = inst.question if where == 'q' else inst.solution
        if needle not in src:
            print(f'  !! cannot plant {cls} in {where}: no {needle!r}')
            ok = False
            continue
        planted = src.replace(needle, repl, 1)
        sol = inst.solution if where == 'q' else planted
        q = planted if where == 'q' else inst.question

        res = ScanResult(template_id=f'<planted:{cls}>',
                         in_track_a=True)
        scan_solution(sol, res, q)
        # `hits` is keyed by the SCANNER's own classes, so indexing it with a
        # T4-only class raises instead of returning nothing.  The expectation
        # has to be computed BEFORE the lookup, and `.get` keeps a future
        # T4-only class from turning a harness gap into a crash.
        scan_should_see = cls not in SCAN_BLIND
        caught_scan = bool(res.hits.get(cls)) if scan_should_see else False
        scan_ok = caught_scan if scan_should_see else True

        mutant = Instance(host.template_id, 0, q, sol)
        t4 = t4_contract.run([mutant], host.template_id)
        t4_should_see = (cls in ('step_marker', 'answer_marker', 'answer_span')
                         and where == 'sol')
        t4_ok = (not t4.passed) if t4_should_see else True

        scan_note = ('n/a' if not scan_should_see
                     else ('CAUGHT' if caught_scan else 'MISSED'))
        print(f'  plant {cls:20s} in {where:3s} {repl[:22]!r:26s} '
              f'scan={scan_note}  '
              f'T4={"FAILS" if not t4.passed else "passes"}')
        if not scan_ok or not t4_ok:
            ok = False

    for cls, find, repl, instrument in NEGATIVE_PLANTS:
        planted = inst.solution.replace(find, repl, 1)
        if instrument == 't4':
            mutant = Instance(host.template_id, 0, inst.question, planted)
            quiet = t4_contract.run([mutant], host.template_id).passed
            census = ''
        else:
            res = ScanResult(template_id=f'<negative:{cls}>')
            scan_solution(planted, res)
            quiet = not res.hits[cls]
            census = f'  [census saw it: {bool(res.hits[cls + chr(95) + "derivation"])}]'
        print(f'  negative plant {cls:13s} ({instrument}) '
              f'{"silent (correct)" if quiet else "FIRED -- over-broad"}{census}')
        if not quiet:
            ok = False
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--selftest', action='store_true',
                    help='plant each defect class and require detection')
    ap.add_argument('--markers', action='store_true',
                    help='D5.3: gold markers vs the candidate-side parser')
    # 400, not 25.  Class 4 fires at 1.95% and a 25-seed scan cannot see it;
    # see the module docstring.  D-024/D-026: the smallest resolvable rate is
    # ~3/N per template.
    ap.add_argument('--seeds', type=int, default=400)
    ap.add_argument('--templates', default='')
    ap.add_argument('--json', default='')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args(argv)

    if args.selftest:
        print('planted-defect self-test')
        return selftest()

    if args.markers:
        print('D5.3 ANSWER_MARKERS agreement')
        return check_marker_agreement()

    refs = discover()
    if args.templates:
        wanted = {t.strip() for t in args.templates.split(',') if t.strip()}
        refs = [r for r in refs if r.template_id in wanted]
        missing = wanted - {r.template_id for r in refs}
        if missing:
            sys.exit(f'unknown template(s): {sorted(missing)}')

    results = [scan_template(r, args.seeds) for r in refs]
    dirty = [r for r in results if not r.clean]

    if not args.quiet:
        for r in dirty:
            print(f'  {r.template_id:48s} {r.summary()}')

    per_class = {c: sorted(r.template_id for r in results if r.hits[c])
                 for c in CLASSES + CENSUS}
    out_of_scope = {c: [t for t in per_class[c] if t not in TRACK_A]
                    for c in CLASSES + CENSUS}

    print()
    print(f'{len(results)} templates x {args.seeds} seeds '
          f'= {sum(r.instances for r in results)} instances')
    print(f'smallest per-template rate this run can resolve: '
          f'~{3 / args.seeds:.2%} (D-024/D-026, ~3/N)')
    for c in CLASSES:
        print(f'  {c:32s} {len(per_class[c]):3d} templates '
              f'({len(out_of_scope[c])} outside Track A)')
    for c in CENSUS:
        print(f'  {c:32s} {len(per_class[c]):3d} templates   [census, not gated]')
        # The names, not just the count.  A census whose members are not printed
        # is a number nobody can act on, and D6.7 is the Phase 6 deliverable that
        # has to act on this one.
        for t in per_class[c]:
            print(f'      {t}')
    print(f'  {"clean":32s} {sum(1 for r in results if r.clean):3d} / {len(results)}')
    n_err = sum(len(r.errors) for r in results)
    print(f'  {"generation errors":32s} {n_err:3d} instances '
          f'on {sum(1 for r in results if r.errors)} templates')

    # The false-positive evidence for the two NEW detectors, recomputed rather
    # than asserted.  A hit outside Track A on `complex_sign` or
    # `degenerate_product` means the detector is over-broad, and it is a finding
    # about this file, not about the template it landed on.
    for c in ('complex_sign', 'degenerate_product'):
        if out_of_scope[c]:
            print(f'  !! {c} fires outside Track A on {out_of_scope[c]} '
                  f'-- detector is over-broad, review it before gating')

    if args.json:
        with open(args.json, 'w', encoding='utf-8') as fh:
            json.dump({
                'seeds': args.seeds,
                'templates': {
                    r.template_id: {
                        'instances': r.instances,
                        'errors': r.errors,
                        'hits': {c: dict(r.hits[c]) for c in CLASSES},
                    } for r in results
                },
                'per_class': per_class,
                'out_of_scope': out_of_scope,
                'resolvable_rate': 3 / args.seeds,
            }, fh, indent=1)

    return 1 if dirty else 0


if __name__ == '__main__':
    raise SystemExit(main())
