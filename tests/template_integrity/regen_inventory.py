"""D6.1 - regenerate the mechanically-measurable columns of template_inventory.csv.

WHY THIS FILE EXISTS
====================
`docs/re-implementation-sep/template_inventory.csv` is scored by Phase 6's exit
gate, but it is byte-identical to the 2026-09-05 pre-Phase-0 audit and no row has
ever been reclassified. D6.1 asked for a "full re-run of the original audit
harness". **That harness does not exist and never did**: the commit that added the
CSV (9105317) contains four documentation files and no `.py`; the commissioning
prompt (`docs/prompts/01-template-structure-audit.md:146`) says "Stay on `master`.
Do not create branches or commits", so the measurement code was throwaway and ran
read-only in-session.

What survives is the *method*, in two places, and this file is an assembly of it:

* `docs/re-implementation-sep/phase0_baseline.md` 2.9b/2.9c - the Phase 0
  adversary independently re-implemented the audit's static and dynamic
  measurements and wrote the predicates out in full. Its AST classifier matched
  the audit EXACTLY on three of five branches (chemical 65, civil 37,
  mechanical 89) and its step-token population came to 17,002 against the
  audit's 17,177 (-1.0%).
* `docs/re-implementation-sep/template_audit_report.md:507-514` - "Appendix -
  method", plus the caveats at :497, :499, :501 and the class definitions
  at :25-30.

Those four numbers are this harness's CALIBRATION GATE. They are asserted, not
assumed: `calibrate()` runs first and `main()` refuses to present the regenerated
columns as a re-audit unless they land. If the method is not first shown
continuous with 2026-09-05, then any divergence from the CSV cannot be attributed
to six phases of template change - it is just as likely to be this file.

WHAT IS MEASURED AND WHAT IS NOT
================================
MEASURED (regenerated from source + execution):
    n_inline_computed, pct_step_values_recoverable, pct_answer_values_recoverable,
    n_steps, has_instance_branching (mechanical base only),
    n_milestone_candidates (a DECLARED PROXY - report :497 says so, and the label
    is kept)

INHERITED, NOT REMEASURED (human judgement; re-deriving means redoing it):
    difficulty, est_effort, notes, answer_type, unit_system, values_bound

RE-CLASSIFIED UNDER A RESTATED RULE, NOT RE-MEASURED:
    instrumentation_class. It is a merge of an automated base with five parallel
    human branch audits under precedence D > C > B > A (report :30). The RULE is
    recoverable from the report; WHICH ROWS THE HUMANS OVERRODE WAS NEVER
    RECORDED. So this file restates the rule and applies it, and labels the
    result as a re-classification. See `classify_restated()` for the limbs it
    can and cannot mechanically evaluate.

Every count here ships with its predicate in this file. That is the whole point
of the deliverable: the original's absence is the defect being fixed.

Run:
    PYTHONIOENCODING=utf-8 python -m tests.template_integrity.regen_inventory
    PYTHONIOENCODING=utf-8 python -m tests.template_integrity.regen_inventory --calibrate-only
"""
from __future__ import annotations

import argparse
import ast
import csv
import importlib
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass, field

from .core import (REPO_ROOT, NUM_RE, TemplateRef, discover, generate,
                   match_value, _walk_values)

INVENTORY = os.path.join(REPO_ROOT, 'docs', 're-implementation-sep',
                         'template_inventory.csv')
OUT_CSV = os.path.join(REPO_ROOT, 'docs', 're-implementation-sep',
                       'template_inventory_regen.csv')

#: Phase 0 2.9c used 15 seeds/template; the audit appendix (:510) says the same.
DYNAMIC_SEEDS = 15
#: The audit appendix (:511) used 25 seeds to diff the blanked skeleton.
STRUCTURE_SEEDS = 25


# ==========================================================================
# SECTION 1 - static: the AST interpolation classifier (phase0 2.9b)
# ==========================================================================
#
# Population: every `ast.FormattedValue` reachable by `ast.walk` from the body
# of a `template_*` function. Module-level code and `main()` excluded. Phase 0
# got 6,166, "an exact match for the audit's population".

#: phase0 2.9b, verbatim: "a wrapper call (round/abs/int/float/str/repr/len)
#: whose FIRST ARGUMENT is itself BOUND".
_WRAPPERS = {'round', 'abs', 'int', 'float', 'str', 'repr', 'len'}

#: The text helpers Phase 0 identified as the entire residual between its 449
#: and the audit's 432: "string-formatting helper calls - fmt(a1, 'h[n-1]'),
#: ', '.join(parts), format_vector(...), equation_str.replace('y','h') - which
#: produce TEXT, not a quantity". Treating these as BOUND is Phase 0's
#: "lenient" variant and yields 375. The audit's 432 sits between the two, so
#: the residual is a definitional choice about text helpers, not a
#: disagreement about arithmetic.
_TEXT_HELPERS = {'join', 'replace', 'format', 'fmt', 'strip', 'lstrip', 'rstrip',
                 'upper', 'lower', 'title', 'capitalize', 'ljust', 'rjust',
                 'center', 'zfill', 'format_vector', 'fmt_term', 'build_poly'}


def is_bound(node: ast.AST, lenient: bool = False) -> bool:
    """phase0_baseline.md 2.9b, stated in full.

    BOUND = Name | Attribute | Subscript | Constant
          | wrapper call (round/abs/int/float/str/repr/len) whose FIRST arg is BOUND
          | IfExp both of whose branches are BOUND
          | nested JoinedStr all of whose parts are BOUND
    INLINE = everything else: BinOp, UnaryOp, Compare, BoolOp, any non-wrapper
             Call, round(a*b, 2), comprehensions.

    NOTE - this is deliberately NOT `checks/t5_binding.py:_is_bound`. That one
    treats `format`/`replace`/`upper`/... as passthrough and tests ALL args
    rather than the first, i.e. it is already the LENIENT variant. Reusing it
    would have silently reproduced 375, not 449, and the calibration gate would
    have failed against Phase 0 for a reason that had nothing to do with the
    templates. `lenient=True` reproduces that variant on purpose.
    """
    if isinstance(node, (ast.Name, ast.Constant, ast.Attribute, ast.Subscript)):
        return True
    if isinstance(node, ast.IfExp):
        return is_bound(node.body, lenient) and is_bound(node.orelse, lenient)
    if isinstance(node, ast.JoinedStr):
        return all(is_bound(p.value, lenient) if isinstance(p, ast.FormattedValue)
                   else True for p in node.values)
    if isinstance(node, ast.Call):
        fname = getattr(node.func, 'id', None) or getattr(node.func, 'attr', None)
        if fname in _WRAPPERS:
            # "whose first argument is itself BOUND"
            return is_bound(node.args[0], lenient) if node.args else True
        if lenient and fname in _TEXT_HELPERS:
            return True
    return False


@dataclass
class StaticResult:
    template_id: str
    interpolations: int = 0
    inline_strict: int = 0
    inline_lenient: int = 0


def _classify_fn(fn: ast.AST, res: StaticResult) -> StaticResult:
    for node in ast.walk(fn):
        if not isinstance(node, ast.FormattedValue):
            continue
        res.interpolations += 1
        if not is_bound(node.value, lenient=False):
            res.inline_strict += 1
        if not is_bound(node.value, lenient=True):
            res.inline_lenient += 1
    return res


def measure_static(ref: TemplateRef) -> StaticResult:
    """Classify every FormattedValue in the template function body."""
    return _classify_fn(ast.parse(ref.source()).body[0],
                        StaticResult(ref.template_id))


# --------------------------------------------------------------------------
# The same predicate, applied to the templates AS THEY STOOD AT A GIVEN REV.
#
# This is what makes the calibration meaningful. Running the classifier against
# TODAY's templates and comparing to Phase 0 conflates two things: a difference
# in METHOD and six phases of template CHANGE. Re-running it against rev
# 9105317 (the commit that added the CSV, 2026-09-05) holds the corpus fixed and
# isolates the method. Only if the predicate reproduces 65/37/89 THERE is it
# legitimate to attribute today's divergence to template change.
#
# Read-only: sources come from `git show`, nothing is checked out.
# --------------------------------------------------------------------------

def _git(*args: str) -> str:
    import subprocess
    return subprocess.run(('git',) + args, cwd=REPO_ROOT, check=True,
                          stdout=subprocess.PIPE).stdout.decode('utf-8', 'replace')


def measure_static_at_rev(rev: str) -> tuple[Counter, int, int, int]:
    """(inline-strict by branch, total interps, total strict, total lenient)."""
    by_branch: Counter = Counter()
    interp = strict = lenient = 0
    listing = _git('ls-tree', '-r', '--name-only', rev,
                   '--', 'data/templates/branches')
    for path in listing.splitlines():
        path = path.strip()
        if not path.endswith('.py'):
            continue
        if os.path.basename(path) in ('constants.py', '__init__.py'):
            continue
        parts = path.split('/')
        if len(parts) < 6:
            continue
        branch = parts[3]
        try:
            src = _git('show', f'{rev}:{path}')
            tree = ast.parse(src)
        except Exception:                                # noqa: BLE001
            continue
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name.startswith('template_'):
                r = _classify_fn(node, StaticResult(node.name))
                by_branch[branch] += r.inline_strict
                interp += r.interpolations
                strict += r.inline_strict
                lenient += r.inline_lenient
    return by_branch, interp, strict, lenient


# ==========================================================================
# SECTION 2 - dynamic: result-position tokens (phase0 2.9c)
# ==========================================================================

#: Phase 0 2.9c: "the terminal answer marker (**Answer / **Final Answer)".
_ANSWER_RE = re.compile(r'\*\*\s*(?:Final\s+)?Answers?\s*:?\s*\*\*|##\s*Final\s+Answer')
_STEP_ANY_RE = re.compile(r'\*\*\s*Step\b')
#: Deliberately LOOSE. The audit (:45-51) records three electrical templates
#: whose markers a strict `\*\*Step (\d+):\*\*` regex silently drops - for
#: `cd_dc_system_analysis` the dropped step is the one that computes the answer.
#: A step-counting harness that inherits the strict regex would under-count
#: exactly the templates the audit flagged, so it must not.
_STEP_NUM_RE = re.compile(r'\*\*\s*Step\s*(\d+)\s*:')


def split_regions(solution: str) -> tuple[str, str]:
    """(step region, answer region) per phase0 2.9c.

    STEP   - first `**Step` marker .. terminal answer marker
    ANSWER - terminal answer marker .. end
    """
    marks = list(_ANSWER_RE.finditer(solution))
    ans_at = marks[-1].start() if marks else len(solution)
    m = _STEP_ANY_RE.search(solution)
    step_at = m.start() if m else 0
    if step_at > ans_at:
        step_at = 0
    return solution[step_at:ans_at], solution[ans_at:]


def step_tokens(step_region: str) -> list[str]:
    """phase0 2.9c STEP rule: for every line containing `=`, THE FIRST numeric
    literal after THE LAST `=` on that line.

    Phase 0 flags this as an inference, not a quotation: it was inferred from
    the audit's own phrase "the value after the last `=` on a line" (singular)
    and then confirmed numerically - 17,002 tokens against the audit's 17,177
    (-1.0%), whereas taking ALL literals after the last `=` yields 29,197
    (+70%). Reproducing the population size to within 1% is what makes the
    recovery percentages comparable at all, so this rule is load-bearing and is
    the single most likely place for this harness to silently diverge.
    """
    out: list[str] = []
    for line in step_region.splitlines():
        if '=' not in line:
            continue
        tail = line[line.rindex('=') + 1:]
        m = NUM_RE.search(tail)
        if m:
            out.append(m.group(0))
    return out


def answer_tokens(answer_region: str) -> list[str]:
    """phase0 2.9c ANSWER rule: every numeric literal on or after the terminal
    answer marker."""
    return [m.group(0) for m in NUM_RE.finditer(answer_region)]


def module_globals(ref: TemplateRef) -> dict[str, float]:
    """phase0 2.9c value set V includes "the defining module's globals".

    `core.generate()` captures only the frame locals, so the globals are added
    here rather than by changing `core.py` (which other work depends on).
    """
    out: dict[str, float] = {}
    containers: dict[str, str] = {}
    try:
        mod = importlib.import_module(ref.module)
    except Exception:                                    # noqa: BLE001
        return out
    for k, v in list(vars(mod).items()):
        if k.startswith('__') or callable(v):
            continue
        if type(v).__name__ == 'module':
            continue
        _walk_values(v, k, out, containers)
    return out


def _dedupe(values: dict[str, float]) -> dict[str, float]:
    """One representative name per distinct float.

    `match_value` scans the whole dict per token and its SCALED limb is
    O(25*n); with globals folded in, the raw dict runs to thousands of entries
    and 21k tokens x 150 templates becomes unrunnable. Deduping on the exact
    float is semantics-preserving: match_value's verdict depends only on the
    VALUES present, never on which name carries them (the name is reported but
    not used in any column regenerated here).
    """
    seen: dict[float, str] = {}
    for name, v in values.items():
        seen.setdefault(v, name)
    return {name: v for v, name in seen.items()}


@dataclass
class DynamicResult:
    template_id: str
    ok_instances: int = 0
    step_total: int = 0
    step_recovered: int = 0
    answer_total: int = 0
    answer_recovered: int = 0
    step_kinds: Counter = field(default_factory=Counter)
    answer_kinds: Counter = field(default_factory=Counter)
    step_counts: set = field(default_factory=set)
    skeletons: set = field(default_factory=set)
    line_counts: set = field(default_factory=set)
    milestones_per_instance: list = field(default_factory=list)
    error: str = ''

    @property
    def branch_evidence(self) -> str:
        """WHICH signal fired, so the over-report can be sized.

        The audit's mechanical base (:511) is "skeleton differs OR step count
        differs", then HAND-CONFIRMED against source. The hand-confirmation is
        what this harness cannot do, so instead of guessing at it the evidence
        is decomposed into a hard limb and a soft one:

        'step-count'  - the number of **Step markers varies. Structural, hard.
        'line-count'  - the solution's line count varies. Structural, hard.
        'prose-only'  - the blanked skeleton differs while BOTH the step count
                        and the line count hold constant. This is the limb that
                        over-reports: a sampled material name ("Pine Wood" vs
                        "Steel"), a pluralisation or a unit word changes the
                        skeleton without any governing-equation branch. These
                        are the rows the audit's hand-confirmation would most
                        likely have rejected, and they must not be presented as
                        equivalent evidence.
        """
        if len(self.step_counts) > 1:
            return 'step-count'
        if len(self.line_counts) > 1:
            return 'line-count'
        if len(self.skeletons) > 1:
            return 'prose-only'
        return ''

    @property
    def pct_step(self) -> int:
        if not self.step_total:
            return 100
        return round(100.0 * self.step_recovered / self.step_total)

    @property
    def pct_answer(self) -> int:
        if not self.answer_total:
            return 100
        return round(100.0 * self.answer_recovered / self.answer_total)

    @property
    def n_steps(self) -> str:
        """Audit renders instance-varying step counts as e.g. `4/6/7` (:94)."""
        return '/'.join(str(c) for c in sorted(self.step_counts)) or '0'

    @property
    def n_milestone_candidates(self) -> int:
        """DECLARED PROXY (report :497): the mean number of result-position
        numbers per instance. It over-counts operand restatements and
        under-counts quantities emitted without an `=`. An ordering signal,
        not an exact count. Kept labelled as one."""
        if not self.milestones_per_instance:
            return 0
        return round(sum(self.milestones_per_instance) /
                     len(self.milestones_per_instance))

    @property
    def has_instance_branching(self) -> bool:
        """MECHANICAL BASE ONLY (audit appendix :511).

        "25 seeds per template, diffing the numeric-blanked solution skeleton
        and step count across seeds to detect instance-dependent structure;
        branching then HAND-CONFIRMED against source."

        The hand-confirmation is not reproducible here, so this is the base
        before it. It over-reports: a template whose prose merely re-words
        (pluralisation, a material name) blanks to a different skeleton without
        any structural branch.
        """
        return len(self.skeletons) > 1 or len(self.step_counts) > 1


def _blank(solution: str) -> str:
    """Numeric-blanked skeleton (audit :511)."""
    return NUM_RE.sub('#', solution)


def measure_dynamic(ref: TemplateRef, use_globals: bool = True,
                    allow_scale: bool = True) -> DynamicResult:
    """`use_globals` / `allow_scale` exist to SIZE this harness's two known
    divergences from phase0 2.9c, not to tune the result:

    * phase0's value set V is "frame locals ... plus the defining module's
      globals". Folding globals in is faithful, but it enlarges V and can only
      ever raise recovery, so the run must be repeatable without them.
    * phase0 classifies SCALED over k in +/-1..6; `core.match_value` uses
      +/-12. Re-running with `allow_scale=False` brackets phase0's rule from
      below, so the true phase0-equivalent figure is known to lie between the
      two rather than being asserted.
    """
    res = DynamicResult(ref.template_id)
    gvals = module_globals(ref) if use_globals else {}
    for seed in range(STRUCTURE_SEEDS):
        capture = seed < DYNAMIC_SEEDS
        inst = generate(ref, seed, capture=capture)
        if not inst.ok:
            res.error = inst.error or 'error'
            continue
        res.skeletons.add(_blank(inst.solution))
        # COUNT OF MARKERS, not of distinct step numbers. The audit's own
        # worked example settles the definition: it records
        # `template_levenspiel_plot_interpretation` as emitting the sequence
        # 1,2,3,1,2,3,4,5,6 (report :44) and carries n_steps=9 for that row -
        # nine markers, six distinct numbers. Deduplicating here would have
        # silently reported 6 and made the one template the audit singles out
        # for restarting its numbering look like a template that does not.
        res.step_counts.add(len(_STEP_NUM_RE.findall(inst.solution)))
        res.line_counts.add(inst.solution.count('\n'))
        if not capture:
            continue
        res.ok_instances += 1
        values = _dedupe({**gvals, **inst.values})
        step_region, ans_region = split_regions(inst.solution)
        st = step_tokens(step_region)
        an = answer_tokens(ans_region)
        res.milestones_per_instance.append(len(st) + len(an))
        for tok in st:
            kind, _ = match_value(tok, values, allow_scale=allow_scale)
            res.step_kinds[kind] += 1
            res.step_total += 1
            if kind != 'MISSING':
                res.step_recovered += 1
        for tok in an:
            kind, _ = match_value(tok, values, allow_scale=allow_scale)
            res.answer_kinds[kind] += 1
            res.answer_total += 1
            if kind != 'MISSING':
                res.answer_recovered += 1
    return res


# ==========================================================================
# SECTION 3 - instrumentation_class, RE-CLASSIFIED UNDER A RESTATED RULE
# ==========================================================================

#: report :27 - "Answer or key intermediates are vector, array, symbolic,
#: classification, or multipart."
_NON_SCALAR = {'vector', 'array', 'symbolic', 'classification', 'multipart'}


def classify_restated(answer_type: str, branching: bool, inherited_class: str,
                      zero_numeric: bool) -> tuple[str, str]:
    """Apply the report's stated rule (:25-30) with precedence D > C > B > A.

    Returns (class, limb_that_decided).

    THIS IS A RE-CLASSIFICATION, NOT A RE-MEASUREMENT, and the distinction is
    the whole reason this function is separate from the measured columns.

    Class D is defined (:28) as "No numeric content, OR the trace's own chain
    does not reproduce its own answer, OR a search/iteration log with no stable
    symbols." Only the FIRST limb is mechanical. The second is printed-
    arithmetic closure over a hand-read chain; the third ("no stable symbols")
    is a judgement about whether a rebound loop variable constitutes an
    identity. So D is inherited except where the mechanical limb fires, and any
    row that was D only by human judgement stays D - this function cannot
    remove a D it cannot evaluate, and must not pretend otherwise.

    Class C keys on `answer_type`, which report :499 states is hand-verified for
    only ~90 of 150 templates and AUTOMATED (and under-detecting) for the rest.
    So C inherits that column's uncertainty wholesale.

    Class B keys on `has_instance_branching`, whose mechanical base IS
    regenerated here - but the audit hand-confirmed it against source (:511)
    and this cannot.
    """
    if inherited_class == 'D' or zero_numeric:
        return 'D', 'D:inherited-judgement' if inherited_class == 'D' else 'D:no-numeric-content'
    if answer_type in _NON_SCALAR:
        return 'C', 'C:answer_type'
    if branching:
        return 'B', 'B:instance-branching'
    return 'A', 'A:residual'


# ==========================================================================
# SECTION 4 - calibration against Phase 0
# ==========================================================================

#: phase0_baseline.md 2.9b table. Three branches the Phase 0 re-implementation
#: matched EXACTLY against the audit; those three are the gate.
PHASE0_INLINE = {'chemical_engineering': 65, 'civil_engineering': 37,
                 'mechanical_engineering': 89}
PHASE0_INLINE_LOOSE = {'electrical_engineering': 127, 'industrial_engineering': 131}
PHASE0_INTERPOLATIONS = 6166
PHASE0_STEP_TOKENS = 17002
PHASE0_ANSWER_TOKENS = 4619


#: The commit that added template_inventory.csv (2026-09-05 19:37 +0500). The
#: corpus as it stood when the audit ran.
AUDIT_REV = '9105317'


def calibrate(statics, dynamics, refs) -> tuple[bool, list[str]]:
    """Reproduce Phase 0's figures BEFORE any new number is presented.

    The gate is on the AT-REV static counts, not the HEAD ones. Comparing HEAD
    against Phase 0 conflates a difference in method with six phases of
    template change; holding the corpus at `AUDIT_REV` isolates the method.
    """
    lines, ok = [], True
    by_branch: dict[str, int] = Counter()
    interp = 0
    for ref, st in zip(refs, statics):
        by_branch[ref.branch] += st.inline_strict
        interp += st.interpolations

    rev_branch, rev_interp, _rev_strict, _rev_len = measure_static_at_rev(AUDIT_REV)

    lines.append(f'CALIBRATION vs phase0_baseline.md 2.9b/2.9c  (gate on this work)')
    lines.append('=' * 74)
    lines.append(f'A. METHOD CONTINUITY - same predicate, corpus held at {AUDIT_REV} '
                 f'(2026-09-05)')
    lines.append('-' * 74)
    lines.append(f'{"metric":<34}{"phase0":>10}{"at-rev":>11}{"delta":>11}')

    d = rev_interp - PHASE0_INTERPOLATIONS
    lines.append(f'{"FormattedValue population":<34}{PHASE0_INTERPOLATIONS:>10}'
                 f'{rev_interp:>11}{d:>+11}')
    if d != 0:
        ok = False
    for br, want in PHASE0_INLINE.items():
        got = rev_branch[br]
        lines.append(f'{"  inline (strict) " + br.split("_")[0]:<34}{want:>10}'
                     f'{got:>11}{got - want:>+11}')
        if got != want:
            ok = False
    for br, want in PHASE0_INLINE_LOOSE.items():
        got = rev_branch[br]
        lines.append(f'{"  inline (strict) " + br.split("_")[0]:<34}{want:>10}'
                     f'{got:>11}{got - want:>+11}   (not gated - phase0 itself '
                     f'differs from the audit here)')

    lines.append('')
    lines.append('B. TEMPLATE DRIFT - same predicate, corpus at HEAD')
    lines.append('-' * 74)
    lines.append(f'{"metric":<34}{"at-rev":>10}{"HEAD":>11}{"delta":>11}')
    lines.append(f'{"FormattedValue population":<34}{rev_interp:>10}'
                 f'{interp:>11}{interp - rev_interp:>+11}')
    for br in list(PHASE0_INLINE) + list(PHASE0_INLINE_LOOSE):
        lines.append(f'{"  inline (strict) " + br.split("_")[0]:<34}'
                     f'{rev_branch[br]:>10}{by_branch[br]:>11}'
                     f'{by_branch[br] - rev_branch[br]:>+11}')
    lines.append('')
    lines.append('C. DYNAMIC - HEAD only (the measurement requires EXECUTING the')
    lines.append('   templates, so it cannot be run against a bare git rev; these')
    lines.append('   are compared to phase0 across the template change in B.')
    lines.append('-' * 74)
    lines.append(f'{"metric":<34}{"phase0":>10}{"HEAD":>11}{"delta":>11}')

    step_total = sum(d_.step_total for d_ in dynamics)
    ans_total = sum(d_.answer_total for d_ in dynamics)
    ds = step_total - PHASE0_STEP_TOKENS
    da = ans_total - PHASE0_ANSWER_TOKENS
    lines.append(f'{"step-token population":<34}{PHASE0_STEP_TOKENS:>10}'
                 f'{step_total:>11}{ds:>+11}')
    lines.append(f'{"answer-token population":<34}{PHASE0_ANSWER_TOKENS:>10}'
                 f'{ans_total:>11}{da:>+11}')
    if abs(ds) > 0.05 * PHASE0_STEP_TOKENS:
        ok = False
    if abs(da) > 0.05 * PHASE0_ANSWER_TOKENS:
        ok = False

    sk = sum(d_.step_kinds[k] for d_ in dynamics for k in ('EXACT', 'ROUNDED', 'SCALED'))
    ak = sum(d_.answer_kinds[k] for d_ in dynamics for k in ('EXACT', 'ROUNDED', 'SCALED'))
    lines.append('')
    lines.append(f'step-token recovery   {100.0 * sk / max(step_total, 1):.1f}%'
                 f'   (phase0 93.2%, audit 94.5%)')
    lines.append(f'answer-token recovery {100.0 * ak / max(ans_total, 1):.1f}%'
                 f'   (phase0 97.6%, audit 96.7%)')
    lines.append('')
    lines.append('GATE: ' + ('PASS - method shown continuous with 2026-09-05'
                             if ok else
                             'FAIL - this harness differs from the original '
                             'method; its numbers are NOT a re-audit'))
    return ok, lines


# ==========================================================================
# SECTION 5 - driver
# ==========================================================================

def read_inventory() -> tuple[list[str], dict[str, dict]]:
    with open(INVENTORY, encoding='utf-8', newline='') as fh:
        rdr = csv.DictReader(fh)
        rows = {r['template_id']: r for r in rdr}
        return list(rdr.fieldnames or []), rows


def diagnose(refs, dynamics, old_rows) -> None:
    """Size this harness's own known weaknesses, rather than asserting them away.

    Three questions, each answered by measurement:
      1. How much of `has_instance_branching` is the over-reporting limb?
      2. Is `n_milestone_candidates` on the same scale as the audit's proxy?
      3. How much of the recovery figure comes from this harness's value set
         being larger than phase0's, rather than from the templates improving?
    """
    print('DIAGNOSTICS')
    print('=' * 74)

    print('\n1. has_instance_branching - decomposition of the mechanical base')
    print('-' * 74)
    ev: Counter = Counter()
    prose_rows, hard_rows = [], []
    for ref, dy in zip(refs, dynamics):
        e = dy.branch_evidence
        ev[e or 'none'] += 1
        was = str(old_rows.get(ref.template_id, {}).get('has_instance_branching', '')).strip()
        if e and was == 'no':
            (prose_rows if e == 'prose-only' else hard_rows).append(
                f'{ref.template_id} ({e})')
    for k in ('step-count', 'line-count', 'prose-only', 'none'):
        print(f'  {k:<14}{ev[k]:>4}')
    print(f'\n  NEWLY "yes" vs the 2026-09-05 CSV: {len(prose_rows) + len(hard_rows)}')
    print(f'    on a HARD signal (step/line count varies): {len(hard_rows)}')
    for r in hard_rows:
        print(f'      {r}')
    print(f'    on the PROSE-ONLY limb (over-reports; the audit hand-confirmed '
          f'and this cannot): {len(prose_rows)}')
    for r in prose_rows:
        print(f'      {r}')

    print('\n2. n_milestone_candidates - is the PROXY on the audit\'s scale?')
    print('-' * 74)
    csv_tot = sum(int(r['n_milestone_candidates']) for r in old_rows.values()
                  if str(r.get('n_milestone_candidates', '')).strip().isdigit())
    new_tot = sum(d.n_milestone_candidates for d in dynamics)
    print(f'  2026-09-05 CSV total : {csv_tot}')
    print(f'  regenerated total    : {new_tot}   ({100.0 * new_tot / csv_tot - 100:+.1f}%)')
    print('  NOTE: the audit calls this a proxy (report :497) and says the')
    print('  per-branch HAND-COUNTS should be preferred where the two disagree.')

    print('\n3. recovery sensitivity - how much is V being larger than phase0\'s?')
    print('-' * 74)
    for label, kw in (('as-measured (globals ON, SCALED k+/-12)',
                       dict(use_globals=True, allow_scale=True)),
                      ('globals OFF (frame locals only)',
                       dict(use_globals=False, allow_scale=True)),
                      ('SCALED OFF (brackets phase0 k+/-6 from below)',
                       dict(use_globals=True, allow_scale=False))):
        ds = [measure_dynamic(r, **kw) for r in refs]
        st = sum(d.step_total for d in ds)
        sr = sum(d.step_recovered for d in ds)
        at = sum(d.answer_total for d in ds)
        ar = sum(d.answer_recovered for d in ds)
        print(f'  {label:<46} step {100.0 * sr / max(st, 1):5.1f}%   '
              f'answer {100.0 * ar / max(at, 1):5.1f}%')
    print('  phase0 2.9c: step 93.2%, answer 97.6%   |   audit: 94.5%, 96.7%')
    print()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate-only', action='store_true')
    ap.add_argument('--diagnose', action='store_true')
    args = ap.parse_args()

    refs = discover()
    print(f'discovered {len(refs)} templates', flush=True)
    statics = [measure_static(r) for r in refs]
    dynamics = []
    for i, r in enumerate(refs, 1):
        dynamics.append(measure_dynamic(r))
        if i % 25 == 0:
            print(f'  ...{i}/{len(refs)} executed', flush=True)

    ok, lines = calibrate(statics, dynamics, refs)
    print()
    print('\n'.join(lines))
    print()
    if args.calibrate_only:
        return 0 if ok else 1

    fields, old = read_inventory()
    if args.diagnose:
        diagnose(refs, dynamics, old)
        return 0 if ok else 1

    # ------------------------------------------------------------------
    # WHICH COLUMNS ARE ALLOWED TO OVERWRITE, AND WHY.
    #
    # `--diagnose` measured three things that disqualify three of the six
    # originally-scoped columns from being written as a re-measurement. They
    # are still emitted - in clearly-labelled side columns - because the
    # measurement is real and the deliverable requires every count to ship
    # with its predicate. What they may not do is silently replace an audit
    # value, because a reader diffing the two files would read the difference
    # as six phases of template change when it is this harness's own bias:
    #
    #   n_milestone_candidates  - regenerated total runs +27.1% against the
    #       CSV's. A proxy on a different scale cannot be diffed row-by-row
    #       against the audit's proxy; the delta would be the scale, not the
    #       template. -> regen_milestone_proxy
    #   has_instance_branching  - of the 34 rows this flips to "yes", 34 fire
    #       on the prose-only limb and 0 on a hard step/line-count signal.
    #       That is an over-reporting detector, and the audit's answer to it
    #       was hand-confirmation against source, which this cannot do.
    #       -> regen_branching_base + regen_branch_evidence
    #   instrumentation_class   - class B keys on the branching base above, so
    #       it inherits that bias wholesale; class D's defining limbs are
    #       human judgement. -> regen_class_restated
    #
    # Overwriting is therefore confined to the columns whose method calibrated
    # EXACTLY against phase0 at the audit rev (n_inline_computed, and n_steps /
    # pct_* which ride the same two calibrated populations).
    # ------------------------------------------------------------------
    OVERWRITE = ('n_inline_computed', 'n_steps',
                 'pct_step_values_recoverable', 'pct_answer_values_recoverable')

    out_fields = fields + ['regen_provenance', 'regen_inline_lenient',
                           'regen_milestone_proxy', 'regen_branching_base',
                           'regen_branch_evidence', 'regen_class_restated',
                           'regen_class_limb']
    changes: list[str] = []
    dist_old: Counter = Counter()
    dist_new: Counter = Counter()
    moved: list[str] = []

    with open(OUT_CSV, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields)
        w.writeheader()
        for ref, st, dy in zip(refs, statics, dynamics):
            row = dict(old.get(ref.template_id, {'template_id': ref.template_id}))
            zero_numeric = (dy.step_total + dy.answer_total) == 0
            new_cls, limb = classify_restated(
                row.get('answer_type', ''), dy.has_instance_branching,
                row.get('instrumentation_class', ''), zero_numeric)
            for col, new in (
                ('n_inline_computed', str(st.inline_strict)),
                ('n_steps', dy.n_steps),
                ('pct_step_values_recoverable', str(dy.pct_step)),
                ('pct_answer_values_recoverable', str(dy.pct_answer)),
            ):
                assert col in OVERWRITE
                if str(row.get(col, '')).strip() != new:
                    changes.append(f'{ref.template_id}  {col}: '
                                   f'{row.get(col, "")!r} -> {new!r}')
                row[col] = new
            dist_old[row.get('instrumentation_class', '')] += 1
            dist_new[new_cls] += 1
            if row.get('instrumentation_class') != new_cls:
                moved.append(f'{ref.template_id}  '
                             f'{row.get("instrumentation_class")} -> {new_cls} ({limb})')
            # NOT written to instrumentation_class - see the note above.
            row['regen_class_restated'] = new_cls
            row['regen_class_limb'] = limb
            row['regen_inline_lenient'] = str(st.inline_lenient)
            row['regen_milestone_proxy'] = str(dy.n_milestone_candidates)
            row['regen_branching_base'] = 'yes' if dy.has_instance_branching else 'no'
            row['regen_branch_evidence'] = dy.branch_evidence
            row['regen_provenance'] = (
                'REGENERATED=n_inline_computed,n_steps,pct_step,pct_answer;'
                'MEASURED-BUT-NOT-SUBSTITUTED=regen_milestone_proxy(scale +27% '
                'vs audit),regen_branching_base(over-reports; see '
                'regen_branch_evidence),regen_class_restated;'
                'INHERITED-NOT-REMEASURED=difficulty,est_effort,notes,'
                'answer_type,unit_system,values_bound,'
                'has_instance_branching,n_milestone_candidates,'
                'instrumentation_class')
            w.writerow(row)

    print(f'wrote {OUT_CSV}')
    print(f'\nREGENERATED CELLS THAT DIFFER FROM THE 2026-09-05 CSV: {len(changes)}')
    for c in changes:
        print('  ' + c)
    print(f'\nCLASS DISTRIBUTION (2026-09-05 CSV, authoritative): '
          + '  '.join(f'{k} {dist_old[k]}' for k in 'ABCD'))
    print(f'CLASS DISTRIBUTION (restated rule, NOT defensible as a '
          f're-measurement - see note in main()):')
    print('    ' + '  '.join(f'{k} {dist_new[k]}' for k in 'ABCD'))
    print(f'\nROWS THE RESTATED RULE WOULD MOVE: {len(moved)}')
    for m in moved:
        print('  ' + m)
    print('\nCLASS D: %d -> %d. The Phase 6 exit gate ("D reduced from 16 to '
          '<=4") is\nNOT met and cannot be evaluated mechanically: two of D\'s '
          'three defining limbs\n(report :28) are human judgement, so this '
          'harness can only inherit D, never\nclear a row out of it.'
          % (dist_old['D'], dist_new['D']))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
