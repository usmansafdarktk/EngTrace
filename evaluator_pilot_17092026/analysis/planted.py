"""X5 - planted defects: a validation set whose ground truth is true by construction.

    python evaluator_pilot_17092026/analysis/planted.py [--labels DIR] [--seed N]

WHAT THIS ESTABLISHES. Finding 5 of RESULTS_X1 is the pilot's central negative result -
no evaluator detects a flawed step behind a correct final answer - and it has two holes
a reviewer can put a finger through.

  1. The expert labels give 178 such flawed steps, of which 175 are calculation slips
     and 3 are conceptual. Conceptual-error detection is therefore not measured at all;
     the set contains no signal to measure it with.
  2. The one method that works on the arithmetic cases, the digit rule
     (analysis/digit_rule.py), implements the rule the annotation guide handed the
     experts - "rounding is not an error, a wrong digit is". Scored against labels
     produced under that guide, it is open to the charge of circularity.

A set of traces where we KNOW what is wrong, because we put it there, answers both. It
owes nothing to the annotation guide, and it holds as many conceptual defects as we care
to plant.

HOW IT IS BUILT. The source is the 135 traces the experts labelled CLEAN - final answer
correct and every step correct or alternative_correct - built from
annotation/score_against_labels.py over experts_filled_labels/version_2. One defect is
planted per trace, surgically, leaving every other character of the trace, and in
particular the final answer, untouched:

  ARITHMETIC  one digit of one displayed intermediate is changed, so the claim no longer
              verifies. Later steps and the final answer keep the original value, so the
              trace's answer is right and its arithmetic is not. Two severities are
              planted and reported apart: `last`, a wrong final digit (a slip at the
              precision shown, the kind the experts actually marked), and `significant`,
              a wrong second significant digit (a gross slip, well outside E4's 1%).
  CONCEPTUAL  no number anywhere in the trace changes. What changes is the reasoning the
              step states - the rule invoked, the criterion applied, the quantity named -
              so that the step now asserts something a competent engineer would call
              wrong while the number path is exactly as it was. Most are stated-versus-
              used mismatches: the step says it is differentiating and shows an integral,
              says Hazen-Williams and substitutes into Manning, says the hydraulic radius
              is P/A and computes A/P.
  CONTROL     untouched clean traces, which is how the false-alarm rate is measured.

Construction is deterministic from one recorded seed, and every planted trace records its
source, the step index, the family, the character offset and the exact before/after text,
so each defect is auditable one line at a time and the set rebuilds byte-identically.

Three checks run on every plant before it is kept, and a plant that fails any of them is
discarded rather than repaired:
  - the final answer is unchanged and still correct under evaluators/answer.py;
  - an arithmetic plant is verified WRONG by a full-precision recomputation that does not
    go through the checker being scored: preferentially against the gold's own milestone
    value for that quantity (milestones.py re-executes the template), otherwise against
    the claim's left side evaluated at full precision;
  - a conceptual plant leaves the trace's digit string byte-identical, and changes exactly
    one span of text.

WHAT IT CANNOT SUPPORT. Three limits, all of them load-bearing.

  - These are planted defects, not observed ones. The set says what an evaluator can
    detect when a defect of a given shape is present. It says nothing about how often
    models produce defects of that shape, and nothing about the distribution of real
    errors. Every rate below is a diagnostic, not a measurement of model behaviour.
  - The digit rule's rate on the ARITHMETIC family is close to an upper bound by
    construction, not an estimate. A plant site has to be a parseable, digit-clean
    arithmetic claim or a stated gold milestone, and parse coverage is precisely what
    bounds the digit rule's recall on real traces (RESULTS_X1: recall 0.472). What the
    family measures honestly is the OTHER evaluators, which have no such advantage.
  - The conceptual defects are written by rule, from a hand-built catalogue. They are
    varied in kind but not in phrasing, and an evaluator could in principle learn the
    catalogue. Nothing here is a held-out test of a trained detector.

Free and deterministic throughout: the digit rule, E3 and E4 called as library functions
on the planted text, and the answer check. Nothing here reimplements a rule it scores -
digit_rule.flagged, e3_milestones.score, e4_arith.score and arith.check are imported, so
a change to any of them shows up in these numbers on the next run. No judge is called,
nothing is written to scores/, and the cost of scoring this set with the paid evaluators
is quoted at the end from the pilot's own per-trace rates rather than incurred.

TWO READINGS OF THE DIGIT RULE are scored side by side, because since D-097 they are
different rules. `digit` is displayed precision alone, which is what Finding 5 of
RESULTS_X1 reports. `e4` is arith.Claim.ok_digit, which E4 now ships: the same rule made
unit-aware (a gold `0.09024 hours = 5.41 minutes` was being flagged as an arithmetic
error) and widened by what the rounded operands the trace shows leave undetermined. The
difference between the two columns, on defects that are certainly present, is exactly
what that fix costs in recall and buys in precision.

The set's own definition of a wrong arithmetic claim is the strict one, and it is the
annotation guide's: the value displayed differs from the full-precision recomputation,
from the numbers the trace itself shows, by more than half the place value of its own
last digit. Where the shipped rule declines to flag such a plant, it is not failing to
see it - it is saying the trace does not show enough to prove it, which is a different
statement and is reported as one.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import glob
import hashlib
import json
import random
import re
import sys
from collections import Counter, defaultdict

sys.path[:0] = [_os.path.join(_PILOT, 'annotation'), _os.path.join(_PILOT, 'evaluators'), _ANALYSIS]
import answer as ANS            # noqa: E402
import arith                    # noqa: E402
import digit_rule as DR         # noqa: E402
import e2_prm                   # noqa: E402
import e3_milestones as E3      # noqa: E402
import e4_arith as E4           # noqa: E402
import milestones as MS         # noqa: E402
import score_against_labels as S  # noqa: E402

SEED = 20260924
OUT = _os.path.join(_ANALYSIS, 'out', 'planted')
TARGET = {'arithmetic': 60, 'conceptual': 60, 'control': 60}

# The step splitter E0 uses, and therefore the one the annotation tasks were cut with.
# Replicated here only to get each step's SPAN in the raw text; the step texts themselves
# come from e2_prm.steps_of, so the indices below are the indices the expert labels use.
STEP_START = re.compile(r'(?:^|\n)\s*(?:\*\*|#*)\s*(?:Step\s*\d+|\d+\.)', re.IGNORECASE)
STEP_PREFIX = re.compile(r'^(?:\*\*|#*)\s*(?:Step\s*\d+|\d+\.)\s*:\s*', re.IGNORECASE)
DIGITS = re.compile(r'\d')
# A displayed value: a decimal, or an integer of four digits or more. Two- and
# three-digit integers are overwhelmingly exponents, indices and step numbers.
LITERAL = re.compile(r'(?<![\w.])(\d+\.\d+|\d{4,})(?![\w.])')


# ----------------------------------------------------------------- the clean source set

def read_traces():
    out = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        for line in open(f, encoding='utf-8'):
            r = json.loads(line)
            if r.get('ok') and r.get('text'):
                out[(r['model_key'], r['item_id'])] = r
    return out


def read_items():
    return {r['item_id']: r for r in
            (json.loads(l) for l in open(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'),
                                         encoding='utf-8'))}


def milestone_store(items, truth, keyfile):
    """The gold's computed quantities per item, as E3 and E4 want them.

    milestones.build RE-EXECUTES the template and refuses an item that no longer
    reproduces byte-identically. 17 of the 60 items no longer do - the repo's templates
    have moved on since the September freeze - so for those the milestones recorded in
    the annotation tasks stand in. Those are the same values the pilot's E3 run and the
    experts both worked from, and they were produced by the template, not by any checker
    scored here, so an arithmetic plant verified against them is still verified
    independently. Which source each item used is reported.
    """
    store, source = {}, {}
    for iid, it in items.items():
        try:
            store[iid] = MS.build(it)['milestones']
            source[iid] = 'rebuilt'
        except Exception:                                    # noqa: BLE001
            store[iid] = None
    for code, t in truth.items():
        iid = keyfile[code][1]
        if store.get(iid) is None and t.get('milestones'):
            store[iid] = [{'id': m['id'], 'value': m['value']} for m in t['milestones']]
            source[iid] = 'frozen in the annotation tasks'
    for iid in list(store):
        if store[iid] is None:
            store[iid] = []
            source[iid] = 'none available'
    state = {'milestones': {iid: {'milestones': store[iid]} for iid in store}}
    return state, Counter(source.values())


def clean_codes(truth):
    """Traces the experts found sound at every level: answer correct, no incorrect step.

    `not_a_claim` is allowed - it is a statement the experts declined to judge, not one
    they judged wrong - but a DISPUTED step (label None) is not, because an unlabelled
    step is not evidence the trace is clean.
    """
    ok = []
    for code, t in truth.items():
        if t['final_answer'] != 'correct':
            continue
        if any(s['label'] not in ('correct', 'alternative_correct', 'not_a_claim')
               for s in t['steps']):
            continue
        ok.append(code)
    return sorted(ok)


def step_spans(text):
    """[(start, end)] in the raw text, one per step, aligned with e2_prm.steps_of.

    extract_steps drops steps whose text is empty once the `**Step 7:**` prefix is
    stripped, so the same steps are dropped here and the surviving indices are the
    indices the expert labels and the digit rule use.
    """
    starts = [m.start() for m in STEP_START.finditer(text)]
    if not starts:
        return [(0, len(text))]
    spans = [(starts[i], starts[i + 1] if i + 1 < len(starts) else len(text))
             for i in range(len(starts))]
    return [sp for sp in spans if STEP_PREFIX.sub('', text[sp[0]:sp[1]].strip()).strip()]


def answer_start(text):
    """Where answer.segment() begins. Nothing at or after this offset may be touched."""
    heads = list(ANS.HEADING.finditer(text))
    if heads:
        return heads[-1].start()
    hits = list(ANS.ANSWER.finditer(text))
    if hits:
        near = len(hits) > 1 and hits[-1].start() - hits[0].start() < ANS.WINDOW
        return hits[0].start() if near else hits[-1].start()
    hits = list(re.finditer(r'(?i)\banswer\b', text))
    return hits[-1].start() if hits else max(0, len(text) - ANS.WINDOW)


def editable(text):
    """(spans, limit): the step spans an edit may land in, and the offset it must precede."""
    spans = step_spans(text)
    limit = answer_start(text)
    return spans, limit


def step_of(spans, offset):
    for i, (a, b) in enumerate(spans):
        if a <= offset < b:
            return i
    return None


# --------------------------------------------------------------- the arithmetic family

def ulp_of(lit):
    """One unit in the last digit shown: '0.703448' -> 1e-6, '23691' -> 1."""
    return 10.0 ** -len(lit.split('.')[1]) if '.' in lit else 1.0


def corrupt(lit, severity):
    """Change exactly one digit. Same length, same decimal point, a different number.

    `last` moves the final digit by 4 (>= 3.5 ulp from the truth whichever way it wraps -
    a wrong last digit, not a rounding). `significant` moves the second significant digit
    by 3, which is a gross slip: 0.703448 -> 0.733448.
    """
    pos = [i for i, ch in enumerate(lit) if ch.isdigit()]
    sig = [i for i in pos if not (lit[:i].replace('.', '').replace('0', '') == ''
                                  and lit[i] == '0')]
    if severity == 'last':
        i, step = pos[-1], 4
    else:
        if len(sig) < 2:
            return None
        i, step = sig[1], 3
    out = lit[:i] + str((int(lit[i]) + step) % 10) + lit[i + 1:]
    return out if out != lit else None


def arith_sites(text, item, spans, limit, ms_values):
    """Every displayed intermediate that could carry a planted digit error.

    A site has to be a value the trace COMPUTED, so that corrupting it is an arithmetic
    error rather than a mis-transcribed input, and it has to be verifiable without the
    checker under test. Two bases qualify, and which a site has is recorded, because the
    difference between them is the whole of the digit rule's recall problem:

      claim      the value is the right-hand side of an arithmetic claim the trace SHOWS,
                 which evaluators/arith.py parses and which is currently correct at the
                 precision shown. Corrupting it breaks a computation the reader can see.
                 Required to appear exactly once in the step, so the corrupted occurrence
                 is certainly the one the claim displays.
      milestone  the value is one of the gold's computed quantities but no parseable
                 computation in the trace ends on it - the trace states the number and
                 the working behind it is prose, a symbolic form, or absent. The gold
                 (milestones.py re-executes the template) says the corrupted value is
                 wrong, without arith.py being consulted at all.
    """
    given = MS.numbers(item['question'])
    sites = []
    for si, (a, b) in enumerate(spans):
        seg = text[a:b]
        claims = [c for c in arith.check(seg).claims if c.left_value and c.right_value]
        for m in LITERAL.finditer(seg):
            off = a + m.start()
            if off >= limit:
                continue
            lit = m.group(1)
            pre = seg[max(0, m.start() - 3):m.start()]
            if '^' in pre or '_' in pre or re.search(r'[eE][-+]?$', pre):
                continue                                   # an exponent or a subscript
            try:
                v = float(lit)
            except ValueError:
                continue
            if v == 0 or any(MS.close(v, g, 1e-12) for g in given):
                continue                                   # a restated input, not a result
            u = ulp_of(lit)
            evidence = None
            if seg.count(lit) == 1:
                for c in claims:
                    if lit not in c.right:
                        continue
                    # The literal has to BE the right-hand side's value, not a piece of a
                    # larger expression. `104/1863 * 60 = 6540/1863` mentions 6540 on the
                    # right and is worth 3.51: changing that digit does not move the side
                    # by one unit in its last place, so a plant there would not be proved
                    # wrong by the recomputation below. Magnitudes, because a minus sign
                    # sits outside the literal.
                    if not any(abs(abs(rv) - v) <= 0.5 * u * (1 + 1e-7) for rv in c.right_value):
                        continue
                    near = [lv for lv in c.left_value
                            if abs(abs(lv) - v) <= 0.5 * u * (1 + 1e-7)]
                    if near:
                        evidence = (c.left, c.right, near[0])
                        break
            on_ms = bool(ms_values) and MS.scaled_match(v, ms_values, MS.DISPLAY_TOL) is not None
            if evidence or on_ms:
                sites.append({'step': si, 'offset': off, 'literal': lit, 'value': v,
                              'ulp': u, 'basis': 'claim' if evidence else 'milestone',
                              'claim': evidence})
    return sites


def verify_arith(site, new_lit, ms_values):
    """Is the corrupted value demonstrably wrong, checked without the digit rule?

    Returns the name of the check that proved it, or None. The milestone check appeals to
    the gold's own computation and touches nothing in evaluators/arith.py; the claim check
    recomputes the claim's left side at full precision.
    """
    new_v = float(new_lit)
    if site['claim']:
        # arith_sites guarantees the right side's magnitude IS this literal and the left
        # side recomputes to it, so the corrupted claim is out by the same amount.
        if abs(abs(site['claim'][2]) - new_v) > 0.5 * site['ulp'] * (1 + 1e-7):
            return 'claim_recompute'
        return None
    if ms_values and MS.scaled_match(site['value'], ms_values, MS.DISPLAY_TOL) is not None \
            and MS.scaled_match(new_v, ms_values, MS.DISPLAY_TOL) is None:
        return 'gold_milestone'
    return None


def plant_arith(code, text, item, ms_values, rng):
    spans, limit = editable(text)
    sites = arith_sites(text, item, spans, limit, ms_values)
    if not sites:
        return None, 'no usable computed value outside the answer segment'
    # Prefer a site inside a computation the trace actually shows: there the defect is
    # visible to a reader, and the plant is verified against the claim's own left side.
    # A stated milestone with no parseable working behind it is the fallback.
    severities = ['last', 'significant']
    rng.shuffle(severities)
    for basis in ('claim', 'milestone'):
        pool = [s for s in sites if s['basis'] == basis]
        if not pool:
            continue
        order = list(range(len(pool)))
        rng.shuffle(order)
        for k in order:
            s = pool[k]
            for sev in severities:
                new_lit = corrupt(s['literal'], sev)
                if not new_lit:
                    continue
                proof = verify_arith(s, new_lit, ms_values)
                if not proof:
                    continue
                new = text[:s['offset']] + new_lit + text[s['offset'] + len(s['literal']):]
                rec = {'family': 'arithmetic', 'step_index': s['step'], 'offset': s['offset'],
                       'before': s['literal'], 'after': new_lit, 'severity': sev,
                       'basis': s['basis'], 'verified_by': proof,
                       'relative_error': abs(float(new_lit) - s['value']) / abs(s['value']),
                       'claim': None if not s['claim'] else
                       {'left': s['claim'][0], 'right': s['claim'][1],
                        'recomputes_to': s['claim'][2]}}
                return (new, rec), None
    return None, 'no digit change survived the independent recomputation'


# --------------------------------------------------------------- the conceptual family

# Each rule rewrites the REASONING a step states and no number in it. `why` is the reason
# the rewritten step is wrong, and is what an expert reading the planted trace should say.
RULES = [
    ('regime_to_laminar', 'criterion',
     r'(flow(?:\s+regime)?\s+(?:is|was)\s+(?:\*\*)?)turbulent', r'\1laminar',
     'the Reynolds number the step just computed is above the critical value, so the '
     'regime it concludes contradicts its own criterion'),
    ('regime_to_turbulent', 'criterion',
     r'(flow(?:\s+regime)?\s+(?:is|was)\s+(?:\*\*)?)laminar', r'\1turbulent',
     'the Reynolds number the step just computed is below the critical value, so the '
     'regime it concludes contradicts its own criterion'),
    ('damping_to_under', 'criterion',
     r'(system\s+is\s+(?:\*\*)?)overdamped', r'\1underdamped',
     'zeta > 1 is overdamped; the step states the classification its own zeta rules out'),
    ('damping_to_over', 'criterion',
     r'(system\s+is\s+(?:\*\*)?)under[- ]?damped', r'\1overdamped',
     'zeta < 1 is underdamped; the step states the classification its own zeta rules out'),
    ('stable_to_unstable', 'criterion',
     r'(?<![A-Za-z])stable(?![A-Za-z])', 'unstable',
     'rho < 1 is the stability condition for M/M/c; the step reads its own test backwards'),
    ('steady_to_transient', 'stated_vs_used',
     r'steady[- ]state(?=\s+(?:balance|equation|distribution))', 'transient',
     'the equation written sets the distribution equal to itself, which is the '
     'steady-state balance; a transient distribution is not obtained that way'),
    ('fixed_to_free', 'rule',
     r'fixed support', 'free end',
     'the slope and deflection vanish at the fixed support of a cantilever, not at the '
     'free end, so the boundary condition is attached to the wrong end'),
    ('free_to_fixed', 'rule',
     r'free end', 'fixed support',
     'the deflection being evaluated is the free-end deflection; at the fixed support it '
     'is zero by the boundary condition the trace already applied'),
    ('integrate_to_differentiate', 'stated_vs_used',
     r'Integrat(e|es|ing|ion|ed)', r'Differentiat\1',
     'the step shows an integral and calls it a differentiation; the two are inverse '
     'operations and the working shown is the integral'),
    ('integrate_to_differentiate_lc', 'stated_vs_used',
     r'integrat(e|es|ing|ion|ed)', r'differentiat\1',
     'the step shows an integral and calls it a differentiation; the two are inverse '
     'operations and the working shown is the integral'),
    ('wrt_y_to_x', 'stated_vs_used',
     r'with respect to \$?y\$?', 'with respect to x',
     'the operation shown is taken with respect to y; the step names the other variable, '
     'and the result it writes is not what that operation would give'),
    ('wrt_x_to_y', 'stated_vs_used',
     r'with respect to \$?x\$?', 'with respect to y',
     'the operation shown is taken with respect to x; the step names the other variable, '
     'and the result it writes is not what that operation would give'),
    ('fx_to_fy', 'rule',
     r'f\(x\)', 'f(y)',
     'integrating with respect to y leaves an arbitrary function of x, not of y'),
    ('hydraulic_radius_inverted', 'stated_vs_used',
     r'R\s*=\s*A\s*/\s*P', 'R = P / A',
     'the hydraulic radius is area over wetted perimeter; the step states the reciprocal '
     'and then divides area by perimeter anyway'),
    ('radius_to_diameter', 'rule',
     r'hydraulic radius', 'hydraulic diameter',
     'the hydraulic diameter is four times A/P; Manning takes the radius'),
    ('wetted_to_top_width', 'rule',
     r'wetted perimeter', 'top width',
     'the top width of a rectangular channel is b, not b + 2y'),
    ('manning_to_hazen', 'stated_vs_used',
     r'Manning[’\'`]?s?\s+equation', 'the Hazen-Williams equation',
     'the formula substituted into is Manning; Hazen-Williams has no 1/n and no S^(1/2)'),
    ('ln_inverted', 'stated_vs_used',
     r'ln\(\s*b\s*/\s*a\s*\)', 'ln(a/b)',
     'ln(a/b) is negative for b > a, which would make the capacitance negative; the step '
     'goes on to divide by ln(b/a)'),
    ('inner_to_outer', 'rule',
     r'inner conductor', 'outer conductor',
     'the Gaussian surface between the conductors encloses the inner conductor only'),
    ('gauss_to_faraday', 'rule',
     r'Gauss[’\'`]?s?\s+law', "Faraday's law",
     "Faraday's law relates an EMF to a changing flux; the radial field of a coaxial line "
     "comes from Gauss's law"),
    ('slope_to_deflection', 'stated_vs_used',
     r'slope equation', 'deflection equation',
     'the first integration of EI y\'\' gives the slope; the deflection needs the second'),
    ('bending_to_shear', 'stated_vs_used',
     r'bending moment', 'shear force',
     'the expression written is the bending moment; EI y\'\' is equal to the moment, not '
     'to the shear'),
    ('dynamic_to_kinematic', 'rule',
     r'dynamic viscosity', 'kinematic viscosity',
     'Re = rho U L / mu takes the dynamic viscosity; the kinematic viscosity already '
     'carries the density'),
    ('natural_to_damped', 'rule',
     r'(?<![Dd]amped )natural frequency', 'damped natural frequency',
     'sqrt(k/m) is the UNdamped natural frequency; the damped one carries sqrt(1 - zeta^2)'),
    ('zeta_inverted', 'stated_vs_used',
     r'(?:\\zeta|ζ)\s*=\s*c\s*/\s*c_?\{?_?c\}?', 'zeta = c_c / c',
     'the damping ratio is c/c_c; the step states the reciprocal and then divides c by c_c'),
    ('transition_to_steady', 'rule',
     r'transition probabilit', 'steady-state probabilit',
     'the one-step transition probabilities are not the steady-state distribution'),
    ('normalisation_to_orthogonality', 'rule',
     r'[Nn]ormali[sz]ation', 'orthogonality',
     'pi_O + pi_D = 1 is normalisation; orthogonality is not a property of a probability '
     'vector'),
    ('utilisation_to_throughput', 'rule',
     r'[Uu]tili[sz]ation', 'throughput',
     'lambda/(c mu) is the utilisation, a dimensionless fraction, not a throughput'),
    ('cross_to_dot', 'stated_vs_used',
     r'cross product', 'dot product',
     'the magnetic force takes u x B; a dot product would be a scalar and could not give '
     'the three components the step writes'),
    ('acceptance_to_rejection', 'rule',
     r'[Pp]robability of acceptance', 'probability of rejection',
     'P(X <= Ac) is the probability of ACCEPTANCE; the rejection probability is its '
     'complement'),
    ('material_to_local', 'stated_vs_used',
     r'[Mm]aterial\s*(?:\(total\)\s*)?(?:derivative|acceleration)', 'local acceleration',
     'the local acceleration is the partial time derivative alone; the step adds the '
     'convective terms and calls the sum local'),
    ('reynolds_to_froude', 'stated_vs_used',
     r'Reynolds number formula', 'Froude number formula',
     'the group substituted into is rho U L / mu; the Froude number is U / sqrt(gL)'),
    ('denominator_to_numerator', 'stated_vs_used',
     r'(?<![Cc]ommon )denominator', 'numerator',
     '(1 + eps X) is the denominator of the gas-phase concentration expression, and the '
     'step divides by it in the next line'),
    ('gaussian_to_spherical', 'rule',
     r'cylindrical Gaussian surface', 'spherical Gaussian surface',
     'the field of a coaxial line has cylindrical symmetry; a spherical surface does not '
     'give E(r) = lambda / (2 pi eps r)'),
    ('potential_to_field', 'stated_vs_used',
     r'potential difference', 'electric field',
     'the step integrates E dr, which gives the potential difference, not the field'),
    ('per_length_to_area', 'rule',
     r'capacitance per unit length', 'capacitance per unit area',
     "C' = 2 pi eps / ln(b/a) is a capacitance per unit LENGTH, in F/m"),
    ('radial_to_axial', 'rule',
     r'radial direction', 'axial direction',
     'the field between coaxial conductors points radially; there is no axial component'),
    ('magnetic_to_electrostatic', 'stated_vs_used',
     r'[Mm]agnetic force', 'electrostatic force',
     'q(u x B) is the magnetic part of the Lorentz force; the electrostatic part is qE '
     'and the step never uses an electric field'),
    ('transition_p_to_r', 'stated_vs_used',
     r'(P\(\s*O\s*(?:->|→)\s*D\s*\)\s*=\s*)p(?![a-z])', r'\1r',
     'p is the breakdown probability and r the repair probability; the step names r for '
     'the O to D transition and then substitutes p'),
    ('polar_to_rectangular', 'stated_vs_used',
     r'polar form', 'rectangular form',
     'the step computes a magnitude and an angle, which is the polar form'),
    ('rectangular_to_polar', 'stated_vs_used',
     r'rectangular form', 'polar form',
     'the step computes real and imaginary parts, which is the rectangular form'),
    ('reduced_T_to_p', 'stated_vs_used',
     r'reduced temperature', 'reduced pressure',
     'T/Tc is the reduced temperature; the Rackett equation takes no reduced pressure'),
    ('molar_volume_to_mass', 'rule',
     r'molar volume', 'molar mass',
     'Vc Zc^((1-Tr)^0.2857) is a molar volume in cm3/mol, not a molar mass'),
]


def conceptual_sites(text, spans, limit):
    out = []
    for rid, kind, pat, rep, why in RULES:
        for m in re.finditer(pat, text):
            si = step_of(spans, m.start())
            if si is None or m.start() >= limit:
                continue
            after = m.expand(rep) if '\\' in rep else rep
            # A rule may quote digits (`P(O->D) = p`) as long as it moves none of them.
            if after == m.group(0) or DIGITS.findall(after) != DIGITS.findall(m.group(0)):
                continue
            out.append({'rule': rid, 'kind': kind, 'step': si, 'offset': m.start(),
                        'before': m.group(0), 'after': after, 'why': why})
            break                     # one site per rule per trace: the first occurrence
    return out


def plant_conceptual(code, text, rng):
    spans, limit = editable(text)
    sites = conceptual_sites(text, spans, limit)
    if not sites:
        return None, 'no reasoning span in this trace matches the defect catalogue'
    # Pick a KIND first, then a site within it, so the three kinds of reasoning defect are
    # spread over the family rather than following whichever rule happens to match most.
    kinds = sorted({s['kind'] for s in sites})
    rng.shuffle(kinds)
    kinds.sort(key=lambda k: k != 'criterion')   # criterion sites are the scarcest
    order = []
    for kind in kinds:
        inner = [s for s in sites if s['kind'] == kind]
        rng.shuffle(inner)
        order += inner
    for s in order:
        new = text[:s['offset']] + s['after'] + text[s['offset'] + len(s['before']):]
        if digit_string(new) != digit_string(text):
            continue
        rec = {'family': 'conceptual', 'step_index': s['step'], 'offset': s['offset'],
               'before': s['before'], 'after': s['after'], 'rule': s['rule'],
               'defect_kind': s['kind'], 'why_wrong': s['why']}
        return (new, rec), None
    return None, 'every catalogue match would have moved a number'


# --------------------------------------------------------------------------- self-checks

def digit_string(text):
    """Every digit in the text, in order. A conceptual plant must not move one of these."""
    return ''.join(DIGITS.findall(text))


def sha(text):
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def check_plant(rec, base, new, item, ms_values):
    """Every property the set claims, asserted per trace. Returns a list of failures."""
    bad = []
    if new == base:
        bad.append('no edit')
    if ANS.verdict(new, item, ms=ms_values)[0] != 'correct':
        bad.append('final answer no longer correct')
    if ANS.segment(new) != ANS.segment(base):
        bad.append('answer segment changed')
    if rec['family'] == 'conceptual':
        if digit_string(new) != digit_string(base):
            bad.append('a number moved')
    else:
        if len(new) != len(base) or sum(a != b for a, b in zip(new, base)) != 1:
            bad.append('not a single-character digit edit')
        db, dn = digit_string(base), digit_string(new)
        if len(db) != len(dn) or sum(a != b for a, b in zip(db, dn)) != 1:
            bad.append('digit string moved by more than one digit')
    if len(step_spans(new)) != len(step_spans(base)):
        bad.append('step split changed')
    return bad


# ------------------------------------------------------------------------ the evaluators

def digit_flags(text, rule):
    """The digit rule, per step, imported not reimplemented, so a fix to it lands here.

    Two readings are scored, because they are now different rules:
      `digit`  displayed precision alone, the reading RESULTS_X1's Finding 5 reports.
      `e4`     arith.Claim.ok_digit, the reading evaluators/e4_arith.py now ships
               (D-097): the same rule made unit-aware and widened by what the rounded
               operands the trace shows leave undetermined.
    """
    return [DR.flagged(s, rule) for s in e2_prm.steps_of(text)]


def e34(state, item, text):
    """E3 and E4 as they ship, called as library functions. Nothing is written to scores/."""
    trace = {'text': text, 'model_key': 'planted', 'item_id': item['item_id']}
    a = E3.score(state, item, trace, 0)['scores']
    b = E4.score(state, item, trace, 0)['scores']
    rep = arith.check(text)
    return {'e3_coverage': a['milestone_coverage'],
            'e4_contradicted': b['milestones_contradicted'],
            'e4_contradicted_tol1': b.get('milestones_contradicted_tol1', 0),
            'e4_coverage': b['e4_coverage'],
            'arith_failed_tol1': rep.checked - rep.consistent,
            'arith_failed_digit': rep.checked - rep.consistent_digit,
            'arith_checked': rep.checked}


def measure(state, item, text, ms_values):
    m = e34(state, item, text)
    for rule in ('digit', 'e4'):
        f = digit_flags(text, rule)
        m['flags_' + rule] = [len(x) for x in f]
        m['total_' + rule] = sum(len(x) for x in f)
        m['n_steps'] = len(f)
    m['answer'] = ANS.verdict(text, item, ms=ms_values)[0]
    return m


# ------------------------------------------------------------------------- the assignment

def balanced_pick(pool, strata, target, taken, already=()):
    """Deterministic, spread-aware selection: always take from the thinnest cell.

    `already` seeds the cell counts with rows this family has been given by an earlier
    pass, so a second pass repairs the first pass's imbalance instead of compounding it.
    `pool` is already shuffled, so ties break on the seed and nothing depends on the
    order the labels happened to load in.
    """
    got, counts = [], [Counter(), Counter(), Counter()]
    for c in already:
        for i in range(3):
            counts[i][strata[c][i]] += 1
    rank = {c: i for i, c in enumerate(pool)}
    while len(got) < target:
        best, bestkey = None, None
        for c in pool:
            if c in taken or c in got:
                continue
            br, mk, lv = strata[c]
            key = (counts[0][br], counts[1][mk], counts[2][lv], rank[c])
            if bestkey is None or key < bestkey:
                best, bestkey = c, key
        if best is None:
            break
        br, mk, lv = strata[best]
        counts[0][br] += 1
        counts[1][mk] += 1
        counts[2][lv] += 1
        got.append(best)
    return got


def spread(rows, strata, title):
    print('\n%s (n = %d)' % (title, len(rows)))
    for i, name in enumerate(('branch', 'model', 'level')):
        c = Counter(strata[r][i] for r in rows)
        print('  %-7s %s' % (name, '  '.join('%s %d' % (k.replace('_engineering', ''), v)
                                             for k, v in sorted(c.items()))))


# ------------------------------------------------------------------------------- the build

def build(a):
    truth, _ = S.build_truth(S.read_labels(a.labels), S.read_consensus(a.labels))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'),
                                             encoding='utf-8'))}
    texts, items = read_traces(), read_items()
    state, msrc = milestone_store(items, truth, keyfile)
    codes = [c for c in clean_codes(truth) if keyfile[c] in texts]
    print('%d traces labelled; %d of them CLEAN (answer correct, no incorrect step)'
          % (len(truth), len(codes)))
    print('  milestones per item: %s' % ', '.join('%s %d' % (k, v) for k, v in msrc.items()))

    ms_of = {c: [m['value'] for m in state['milestones'][keyfile[c][1]]['milestones']]
             for c in codes}
    base = {c: texts[keyfile[c]]['text'] for c in codes}
    # The set's answers must be correct under BOTH judgements. answer.py agrees with the
    # experts on 0.947 of non-partial traces (Finding 1b); where it does not, the trace is
    # dropped rather than planted in, so no planted row rests on a contested answer.
    drop = [c for c in codes
            if ANS.verdict(base[c], items[keyfile[c][1]], ms=ms_of[c])[0] != 'correct']
    codes = [c for c in codes if c not in set(drop)]
    print('  %d dropped: the experts call the answer correct and answer.py does not; '
          '%d sources remain' % (len(drop), len(codes)))
    strata = {c: (truth[c]['branch'], keyfile[c][0], items[keyfile[c][1]]['level'])
              for c in codes}

    # Feasibility first, so the assignment only ever offers a trace a plant can land in.
    plants, misses = {'arithmetic': {}, 'conceptual': {}}, []
    for c in codes:
        item = items[keyfile[c][1]]
        for fam, fn in (('arithmetic', lambda: plant_arith(c, base[c], item, ms_of[c],
                                                           random.Random('%s|%d|a' % (c, a.seed)))),
                        ('conceptual', lambda: plant_conceptual(
                            c, base[c], random.Random('%s|%d|c' % (c, a.seed))))):
            got, why = fn()
            if got is None:
                misses.append((c, fam, why))
                continue
            new, rec = got
            bad = check_plant(rec, base[c], new, item, ms_of[c])
            if bad:
                misses.append((c, fam, 'self-check: ' + '; '.join(bad)))
                continue
            rec.update(code=c, model_key=keyfile[c][0], item_id=keyfile[c][1],
                       branch=truth[c]['branch'], level=strata[c][2],
                       template_id=item['template_id'], n_steps=len(truth[c]['steps']),
                       source_sha256=sha(base[c]), planted_sha256=sha(new))
            plants[fam][c] = (new, rec)

    # Sources are split between the two planted families, never shared, so no trace
    # contributes two different defects to the same table. A trace only one family can
    # use goes to that family first; the rest are split, thinnest stratum first.
    rng = random.Random(a.seed)
    pool = list(codes)
    rng.shuffle(pool)
    fits = {f: [c for c in pool if c in plants[f]] for f in ('arithmetic', 'conceptual')}
    other = {'arithmetic': 'conceptual', 'conceptual': 'arithmetic'}
    sel = {}
    for f in fits:                                           # traces only one family fits
        sel[f] = balanced_pick([c for c in fits[f] if c not in plants[other[f]]],
                               strata, TARGET[f], set())
    # Then one trace at a time, alternating, so neither family empties a thin cell -
    # electrical has 14 clean traces in all - before the other has had its share of it.
    while any(len(sel[f]) < TARGET[f] and
              [c for c in fits[f] if c not in set(sel['arithmetic']) | set(sel['conceptual'])]
              for f in fits):
        for f in sorted(fits, key=lambda f: len(fits[f])):
            if len(sel[f]) >= TARGET[f]:
                continue
            taken = set(sel['arithmetic']) | set(sel['conceptual'])
            sel[f] += balanced_pick([c for c in fits[f] if c not in taken], strata,
                                    1, taken, sel[f])
    arit, conc = sel['arithmetic'], sel['conceptual']
    ctrl = balanced_pick(pool, strata, TARGET['control'], set())
    return dict(truth=truth, keyfile=keyfile, items=items, state=state, codes=codes,
                strata=strata, ms_of=ms_of, base=base, plants=plants, misses=misses,
                conc=conc, arit=arit, ctrl=ctrl)


def write_set(b, seed):
    _os.makedirs(OUT, exist_ok=True)
    path = _os.path.join(OUT, 'planted.jsonl')
    n = 0
    with open(path, 'w', encoding='utf-8') as fh:
        for fam, sel in (('arithmetic', b['arit']), ('conceptual', b['conc'])):
            for c in sel:
                new, rec = b['plants'][fam][c]
                fh.write(json.dumps(dict(rec, set_id='%s:%s' % (fam[:4], c), text=new),
                                    ensure_ascii=False) + '\n')
                n += 1
        for c in b['ctrl']:
            rec = {'set_id': 'ctrl:%s' % c, 'family': 'control', 'code': c,
                   'model_key': b['keyfile'][c][0], 'item_id': b['keyfile'][c][1],
                   'branch': b['strata'][c][0], 'level': b['strata'][c][2],
                   'template_id': b['items'][b['keyfile'][c][1]]['template_id'],
                   'step_index': None, 'offset': None, 'before': None, 'after': None,
                   'n_steps': len(b['truth'][c]['steps']), 'source_sha256': sha(b['base'][c]),
                   'planted_sha256': sha(b['base'][c]), 'text': b['base'][c]}
            fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
            n += 1
    with open(_os.path.join(OUT, 'build.json'), 'w', encoding='utf-8') as fh:
        json.dump({'seed': seed, 'targets': TARGET, 'n_rows': n,
                   'n_clean_sources': len(b['codes']),
                   'rules': [{'id': r[0], 'kind': r[1], 'why': r[4]} for r in RULES],
                   'not_planted': [{'code': c, 'family': f, 'reason': w}
                                   for c, f, w in b['misses']]},
                  fh, indent=1)
    print('\nwrote %d rows to %s and the build record beside it'
          % (n, _os.path.relpath(path, _PILOT)))
    return path


# ------------------------------------------------------------------------------ the table

def report(b):
    state, items, ms_of, base = b['state'], b['items'], b['ms_of'], b['base']
    rows = []
    for fam, sel in (('arithmetic', b['arit']), ('conceptual', b['conc'])):
        for c in sel:
            new, rec = b['plants'][fam][c]
            item = items[b['keyfile'][c][1]]
            rows.append({'family': fam, 'code': c, 'rec': rec,
                         'before': measure(state, item, base[c], ms_of[c]),
                         'after': measure(state, item, new, ms_of[c])})
    ctrl = []
    for c in b['ctrl']:
        item = items[b['keyfile'][c][1]]
        ctrl.append({'code': c, 'm': measure(state, item, base[c], ms_of[c])})

    print('\nANSWER CHECK (evaluators/answer.py) - the set is only valid if this is clean')
    bad = [r['code'] for r in rows if r['after']['answer'] != 'correct']
    print('  planted traces whose final answer is still correct: %d of %d%s'
          % (len(rows) - len(bad), len(rows), '' if not bad else '  FAILED: ' + str(bad)))
    print('  control traces whose final answer is correct:       %d of %d'
          % (sum(1 for r in ctrl if r['m']['answer'] == 'correct'), len(ctrl)))

    def det(r, key):
        """A NEW flag the same evaluator does not raise on the same trace unmodified."""
        i = r['rec']['step_index']
        if key in ('digit_step', 'e4_step'):
            f = 'flags_' + ('digit' if key == 'digit_step' else 'e4')
            b_, a_ = r['before'][f], r['after'][f]
            return i is not None and i < len(a_) and i < len(b_) and a_[i] > b_[i]
        if key in ('digit_any', 'e4_any'):
            f = 'total_' + ('digit' if key == 'digit_any' else 'e4')
            return r['after'][f] > r['before'][f]
        if key in ('claim_digit', 'claim_tol1'):
            f = 'arith_failed_' + ('digit' if key == 'claim_digit' else 'tol1')
            return r['after'][f] > r['before'][f]
        if key in ('ms_digit', 'ms_tol1'):
            f = 'e4_contradicted' + ('' if key == 'ms_digit' else '_tol1')
            return r['after'][f] > r['before'][f]
        return r['after']['e3_coverage'] < r['before']['e3_coverage']

    heads = [('digit rule as X1 ran it, planted step', 'digit_step'),
             ('digit rule as X1 ran it, any new flag', 'digit_any'),
             ('digit rule as E4 now ships it, planted step', 'e4_step'),
             ('E4 arithmetic, shipped digit reading', 'claim_digit'),
             ('E4 arithmetic, the 1% reading kept beside it', 'claim_tol1'),
             ('E4 milestone contradicted, shipped', 'ms_digit'),
             ('E4 milestone contradicted, 1% reading', 'ms_tol1'),
             ('E3 milestone coverage falls', 'e3')]
    fa = {'digit_step': lambda m: m['total_digit'] > 0,
          'digit_any': lambda m: m['total_digit'] > 0,
          'e4_step': lambda m: m['total_e4'] > 0,
          'claim_digit': lambda m: m['arith_failed_digit'] > 0,
          'claim_tol1': lambda m: m['arith_failed_tol1'] > 0,
          'ms_digit': lambda m: m['e4_contradicted'] > 0,
          'ms_tol1': lambda m: m['e4_contradicted_tol1'] > 0,
          'e3': lambda m: m['e3_coverage'] < 1.0}

    print('\nDETECTION - a defect is detected when the evaluator raises a flag on the')
    print('planted trace that it does not raise on the same trace unmodified.')
    print('  %-44s %14s %14s %14s' % ('evaluator', 'ARITHMETIC', 'CONCEPTUAL',
                                      'CONTROL false alarm'))
    for name, key in heads:
        cells = []
        for fam in ('arithmetic', 'conceptual'):
            sel = [r for r in rows if r['family'] == fam]
            k = sum(det(r, key) for r in sel)
            cells.append('%3d/%-3d %.3f' % (k, len(sel), k / len(sel) if sel else 0))
        k = sum(fa[key](r['m']) for r in ctrl)
        cells.append('%3d/%-3d %.3f' % (k, len(ctrl), k / len(ctrl) if ctrl else 0))
        print('  %-44s %14s %14s %14s' % (name, *cells))
    print('  false alarm = the evaluator fires anywhere in an untouched clean trace. E3\'s')
    print('  row is not a false alarm in the same sense: coverage below 1.0 only means the')
    print('  trace did not state every milestone, which is common and is not an error flag.')

    ar = [r for r in rows if r['family'] == 'arithmetic']
    print('\nARITHMETIC, split by where the defect sits and how big it is')
    print('  %-46s %4s %10s %11s %11s %9s'
          % ('', 'n', 'median rel', 'digit X1', 'digit E4', 'E4 1%'))
    groups = []
    for basis, label in (('claim', 'inside a claim the checker parses'),
                         ('milestone', 'a stated value with no parseable working')):
        groups.append((label, lambda r, b=basis: r['rec']['basis'] == b))
        groups.append(('   ...off by less than E4\'s 1 percent',
                       lambda r, b=basis: r['rec']['basis'] == b
                       and r['rec']['relative_error'] < arith.ARITH_TOL))
        groups.append(('   ...off by 1 percent or more',
                       lambda r, b=basis: r['rec']['basis'] == b
                       and r['rec']['relative_error'] >= arith.ARITH_TOL))
    for name, fn in groups:
        sel = [r for r in ar if fn(r)]
        if not sel:
            print('  %-46s %4d' % (name, 0))
            continue
        rels = sorted(r['rec']['relative_error'] for r in sel)
        print('  %-46s %4d %10.2e %11.3f %11.3f %9.3f'
              % (name, len(sel), rels[len(rels) // 2],
                 sum(det(r, 'digit_step') for r in sel) / len(sel),
                 sum(det(r, 'e4_step') for r in sel) / len(sel),
                 sum(det(r, 'claim_tol1') for r in sel) / len(sel)))
    byb = Counter(r['rec']['verified_by'] for r in ar)
    print('  proved wrong without the checker under test: %s.  digit changed: %s'
          % (', '.join('%s %d' % kv for kv in sorted(byb.items())),
             ', '.join('%s %d' % kv for kv in
                       sorted(Counter(r['rec']['severity'] for r in ar).items()))))
    print('  the split IS the result: the digit rule finds what arith.py can parse and')
    print('  nothing else, which is the parse-coverage ceiling RESULTS_X1 reports as recall 0.472.')
    missed = [r for r in ar if r['rec']['basis'] == 'claim' and not det(r, 'digit_step')]
    print('  planted INSIDE a parsed claim and missed by the X1 reading: %d - an arith.py '
          'parse gap' % len(missed))
    for r in missed[:3]:
        q = r['rec']
        print('      %s %s step %s: %r -> %r, claim %s'
              % (q['model_key'], q['item_id'], q['step_index'], q['before'], q['after'],
                 str(q['claim']['left'])[:46]))
    # What D-097 costs in recall, measured on defects that are certainly there. The
    # shipped rule widens the band by what the trace's own rounded operands leave
    # undetermined, so a last-digit slip it declines to flag is one the trace does not
    # show enough to prove - not a miss in the ordinary sense, and worth separating.
    give = [r for r in ar if det(r, 'digit_step') and not det(r, 'e4_step')]
    print('  flagged by the X1 reading and NOT by the shipped one: %d (%s)'
          % (len(give), ', '.join('%s %d' % kv for kv in
                                  sorted(Counter(r['rec']['severity'] for r in give).items()))
             or 'none'))
    for r in give[:3]:
        q = r['rec']
        print('      %s %s step %s: %r -> %r (%.1e relative), claim %s'
              % (q['model_key'], q['item_id'], q['step_index'], q['before'], q['after'],
                 q['relative_error'], str(q['claim']['left'])[:40] if q['claim'] else '-'))

    co = [r for r in rows if r['family'] == 'conceptual']
    print('\nCONCEPTUAL, by the kind of reasoning corrupted')
    print('  %-22s %5s %12s %12s %12s %10s'
          % ('kind', 'n', 'digit X1', 'digit E4', 'E4 1%', 'E3 falls'))
    for kind in sorted({r['rec']['defect_kind'] for r in co}):
        sel = [r for r in co if r['rec']['defect_kind'] == kind]
        print('  %-22s %5d %12.3f %12.3f %12.3f %10.3f'
              % (kind, len(sel),
                 sum(det(r, 'digit_any') for r in sel) / len(sel),
                 sum(det(r, 'e4_any') for r in sel) / len(sel),
                 sum(det(r, 'claim_tol1') for r in sel) / len(sel),
                 sum(det(r, 'e3') for r in sel) / len(sel)))
    print('  rules used: %d of %d in the catalogue'
          % (len({r['rec']['rule'] for r in co}), len(RULES)))

    print('\nWHAT THE FALSE ALARMS ARE - flags on untouched traces every expert called clean')
    nsteps = sum(r['m']['n_steps'] for r in ctrl)
    for rule, label in (('digit', 'digit rule as X1 ran it'),
                        ('e4', 'digit rule as E4 now ships it')):
        nflag = sum(r['m']['total_' + rule] for r in ctrl)
        print('  %-32s %3d flags over %d steps (%.1f per 100 steps)'
              % (label, nflag, nsteps, 100.0 * nflag / nsteps if nsteps else 0))

    print('\nEXAMPLES - how the set records a defect, two of each family')
    for fam, pick in (('arithmetic', lambda r: r['rec']['basis'] == 'claim'),
                      ('arithmetic', lambda r: r['rec']['basis'] == 'milestone'),
                      ('conceptual', lambda r: r['rec']['defect_kind'] == 'stated_vs_used'),
                      ('conceptual', lambda r: r['rec']['defect_kind'] == 'criterion')):
        sel = [r for r in rows if r['family'] == fam and pick(r)]
        if not sel:
            continue
        q = sel[0]['rec']
        print('  [%s] %s %s step %s: %r -> %r'
              % (fam[:4], q['model_key'], q['item_id'], q['step_index'],
                 q['before'], q['after']))
        if fam == 'conceptual':
            why = q['why_wrong']
        elif q.get('claim'):
            why = 'the claim %s = %s now recomputes to %.10g' % (
                str(q['claim']['left'])[:44], q['after'], q['claim']['recomputes_to'])
        else:
            why = ('the gold computes this quantity as %s; %s is not it, and no parseable '
                   'claim shows the working' % (q['before'], q['after']))
        print('         %s' % why)
    return rows, ctrl


def not_planted(b):
    print('\nTRACES NO DEFECT COULD BE PLANTED IN')
    print('  feasible sources, out of the clean traces each cell holds:')
    for i, name in enumerate(('branch', 'level')):
        j = 0 if name == 'branch' else 2
        cells = sorted({b['strata'][c][j] for c in b['codes']})
        print('    %-7s %s' % (name, '  '.join(
            '%s %d/%d/%d' % (k.replace('_engineering', ''),
                             sum(1 for c in b['plants']['arithmetic'] if b['strata'][c][j] == k),
                             sum(1 for c in b['plants']['conceptual'] if b['strata'][c][j] == k),
                             sum(1 for c in b['codes'] if b['strata'][c][j] == k))
            for k in cells)))
    print('    read as arithmetic/conceptual/clean sources; a thin cell in the spread '
          'above is a thin cell here.')
    by = defaultdict(Counter)
    for c, fam, why in b['misses']:
        by[fam][why] += 1
    for fam in ('arithmetic', 'conceptual'):
        n = sum(by[fam].values())
        print('  %-12s %3d of %d clean traces' % (fam, n, len(b['codes'])))
        for why, k in by[fam].most_common():
            print('      %3d  %s' % (k, why))


def judge_quote(n):
    """What it would cost to put the paid evaluators on this set. Nothing is called."""
    print('\nIF THE PAID EVALUATORS WERE LATER RUN ON THIS SET (nothing was run here)')
    try:
        import judge_cost as JC
    except Exception as e:                                   # noqa: BLE001
        print('  judge_cost unavailable (%s)' % e)
        return
    print('  %d traces to score.  %-8s %10s %10s %10s %10s'
          % (n, '', '$/trace F', '$/trace W', 'set at F', 'set at W'))
    for name, d, _desc in JC.EVALS:
        try:
            rows = JC.current(d)
        except Exception:                                    # noqa: BLE001
            print('  %-42s no scored rows to take a rate from' % name)
            continue
        per = {}
        for grp, keys in (('F', JC.FRONTIER), ('W', JC.WEAK)):
            rs = [r for r in rows if r['model_key'] in keys]
            per[grp] = (sum(sum(c.get('cost_usd', 0) for c in r.get('calls') or [])
                            for r in rs) / len(rs)) if rs else 0.0
        print('  %-42s %10.5f %10.5f %10.2f %10.2f'
              % (name, per['F'], per['W'], per['F'] * n, per['W'] * n))
    print('  F = the pilot frontier models\' per-trace rate, W = Llama 3.1 70B\'s.')
    print('  E2 needs a GPU and costs nothing in API terms; E3, E4 and the digit rule are free.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels',
                                                      'version_2', 'labels'))
    ap.add_argument('--seed', type=int, default=SEED)
    a = ap.parse_args()
    print('PLANTED DEFECTS - seed %d, targets %s' % (a.seed, TARGET))
    b = build(a)
    spread(b['arit'], b['strata'], 'SPREAD, arithmetic plants')
    spread(b['conc'], b['strata'], 'SPREAD, conceptual plants')
    spread(b['ctrl'], b['strata'], 'SPREAD, controls')
    ov = len(set(b['arit']) & set(b['conc']))
    print('\n  arithmetic and conceptual sources are disjoint (%d shared); %d of the %d '
          'controls are\n  the untouched original of a trace planted in the other two '
          'families, which is the\n  paired baseline the detection column is measured '
          'against.' % (ov, len(set(b['ctrl']) & (set(b['arit']) | set(b['conc']))),
                        len(b['ctrl'])))
    write_set(b, a.seed)
    report(b)
    not_planted(b)
    judge_quote(len(b['arit']) + len(b['conc']) + len(b['ctrl']))


if __name__ == '__main__':
    main()
