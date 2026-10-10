"""The answer check audited by answer kind (the supervisor's comment of 7 October, "Answer Scoring Validity").

    python -m full_run_28092026.answer_audit              # writes ANSWER_AUDIT.md beside this file (counts only),
                                                          # results/answer_audit.json and results/answer_audit_changes.json,
                                                          # and, locally, expert_request/ANSWER_AUDIT_CASES.md
    python -m full_run_28092026.answer_audit --selftest   # the readers and variants on hand-made cases; writes nothing

WHAT IT MEASURES
 1. Agreement by answer kind, on the two populations the experts labelled: the expert study (the pilot's 300
    responses from five LLMs outside the evaluated eleven, three-way and on the responses the experts did not call
    partial, as validate_scorer.py computes them) and the domain experts' readings of the evaluated models (B1 as
    `expert_kits.py --score` counts it now: the readings the main store still holds, against the check's current
    verdict). The totals must reproduce 0.930, 0.986 and 244 of 288, or the run stops.
 2. Every disagreement, by cause. CAUSES (expert_request/answer_audit_causes.json, local, because it names the items
    the experts read) gives each disagreeing study response and each disagreeing reading one cause, assigned by
    reading the response, the gold and the experts' notes. The run stops if a disagreement has no cause or a cause
    names an agreement. The cases, with the notes, go to expert_request/ANSWER_AUDIT_CASES.md (local).
 3. What the check takes as targets, item by item over the 2,250: targets whose own last-digit window reaches zero;
    arrays scored on their last value; values the gold's answer line states that are not targets, by why (a sign
    printed apart from its number, a value not among the computed ones, a value of 0, 1 or 2); label parts the
    question asks for that no target checks (LABEL_PARTS); questions that name one value to report among other
    requested results (REPORTED_VALUE).
 4. Variants of the rule on the populations of 1 and on the eleven models' headline responses (matched settings:
    each model at its reasoning store where it has one, else at main). `verdict()` is a
    local copy of answer.verdict with switches; with every switch at the rule it must reproduce the stored label of
    every readable headline row and the check's label on every study response, or the run stops. Each variant
    changes one thing (VARIANTS); three combine them. For each: the study's and the readings' agreement, the
    verdicts it moves on the headline responses, Final Answer Accuracy per model and Kendall's tau with the ordering as
    scored. The moved verdicts go to results/answer_audit_changes.json for reading.
 5. Which way of matching each accepted headline answer needs (within epsilon, a last-digit window, a unit factor,
    the absolute value, a rendering, the prescribed digits), by answer kind.

Reads the score stores, the traces (each answered row's recorded hash checked by score.texts_matching), the milestone
cache, the pilot's labels and traces and the experts' B1 returns (local). Calls no model.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures as cf
import datetime as dt
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PILOT = REPO / 'evaluator_pilot_17092026'
for p in (str(REPO), str(PILOT / 'evaluators'), str(PILOT / 'annotation'), str(PILOT / 'analysis'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer as A  # noqa: E402
import milestones as M  # noqa: E402

from full_run_28092026 import clause_variants as CV, expert_kits as K, score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

REPORT = HERE / 'ANSWER_AUDIT.md'
OUT_JSON = HERE / 'results' / 'answer_audit.json'
CHANGES_JSON = HERE / 'results' / 'answer_audit_changes.json'
CAUSES = K.OUT / 'answer_audit_causes.json'
CASES_MD = K.OUT / 'ANSWER_AUDIT_CASES.md'
KINDS = ('scalar', 'multipart', 'vector', 'array', 'symbolic', 'classification')
SCORE_OF = {'correct': 1.0, 'partial': 0.5, 'incorrect': 0.0}
ORDER = {'incorrect': 0, 'partial': 1, 'correct': 2}

# ------------------------------------------------------------------ template tables (read from the questions)

# Label parts a question asks for that no target checks: the check takes words only on classification items, and only
# the verdict words its gold line states. (critical_depth_froude_classification asks for the regime only "in your
# solution", beside the Froude number it names as the answer: it is in REPORTED_VALUE, not here.)
LABEL_PARTS = {
    'template_decimation_aliasing_analysis': 'whether aliasing occurred',
    'template_limiting_reactant': 'which reactant is limiting',
    'template_signal_energy_power': 'whether the signal is an energy or a power signal',
    'template_system_properties': 'the system type (underdamped, critically damped or overdamped)',
    'template_vorticity_check': 'whether the flow is rotational or irrotational',
    'template_wave_equation_interpretation': 'the direction of propagation',
}
# Questions that ask for several results and name one value to report: the check scores the named value, and the
# other requested results enter Milestone Coverage where they are milestones.
REPORTED_VALUE = {
    'template_absorbing_chain_time_to_failure', 'template_aoq_ati_rectifying', 'template_arl_beta_mean_shift',
    'template_basic_eoq', 'template_c_chart_revision', 'template_chart_pair_selection',
    'template_chase_vs_level_aggregate', 'template_critical_depth_froude_classification',
    'template_epq_finite_production', 'template_line_balancing_heuristic', 'template_mm1_time_in_system',
    'template_mm1k_finite_capacity', 'template_mmc_waiting_time', 'template_newsvendor_normal_demand',
    'template_p_chart_limits_floor', 'template_poisson_event_count', 'template_qr_policy_one_iteration',
    'template_safety_stock_reorder_point', 'template_server_configuration_selection', 'template_sigma_reduction_for_cpk',
    'template_takt_time_line_efficiency', 'template_xbar_known_sigma_classification', 'template_xbar_r_control_limits',
    'template_cp_cpk_from_specs',
}
# Where the sign of a target is a convention rather than part of the answer: the gold states a magnitude, or the
# question asks for a magnitude and the direction separately. Elsewhere the sign-policy variant accepts a value of the
# other sign only when the response states the sign in words (SIGN_WORDS).
SIGN_CONVENTION = {
    'template_cantilever_double_integration': 'the gold states the deflection at the free end as a magnitude',
    'template_wave_equation_interpretation': 'the phase velocity is asked as a speed and the direction separately',
    'template_beam_deflection_formula': 'the gold states the maximum deflection as a magnitude',
    'template_axial_deformation': 'the gold states the deformation as a magnitude with tension or compression in words',
}
SIGN_WORDS = {-1: re.compile(r'(?i)compress|contract|shorten|downward|decrease|below|loss|drop|negative'),
              1: re.compile(r'(?i)tension|tensile|elongat|expan|upward|increase|above|gain|rise|positive')}

# ------------------------------------------------------------------ readers

JMINUS = re.compile(r'[-−–]\s*j\s*(?:\\[,;:!]\s*)?(?=\d)')
JPLUS = re.compile(r'\+\s*j\s*(?:\\[,;:!]\s*)?(?=\d)')
SPACED_MINUS = re.compile(r'(?<![\d.\s])(\s*)-\s+(?=\d)')     # "x_hat - 9.89e-05": the minus belongs to the number
FACTORED = re.compile(r'\(([^()]{3,400})\)\s*(?:\\times|×|\*|\\cdot|·|x)\s*10\s*\^\s*\{?\s*\(?\s*([-+]?\d+)\s*\)?\s*\}?')
NO_UNITS = re.compile(r'(?i)\bno units?\b(?:\s+(?:applicable|apply|needed|required))?|\bunits?\s*:\s*none\b|'
                      r'\b(?:unitless|dimensionless)\b')
LIST_SPAN = re.compile(r'\b(?:between|from)\s+(\w+)\s+(?:and|to)\s+(\w+)|\b(\w+)\s+or\s+(\w+)\b')


def signed(s: str) -> str:
    """Minus signs the reader drops, joined to their numbers: U+2212 and the en dash; a minus set apart by spaces
    (or LaTeX spacing, '\\;-\\;') after a non-numeric token ('x_hat - 9.890e-05', '\\hat{x} - 1.17'); '- j85.44' and
    '44.26-j44.88' for an imaginary part (a space before the joined minus, so that a digit before it does not make it a
    subtraction for the number reader)."""
    s = s.replace('−', '-').replace('–', '-')
    s = re.sub(r'(?:\\[,;:! ])+(?=\s*-)', ' ', s)
    s = re.sub(r'(?<=-)(?:\s*\\[,;:! ])+', ' ', s)
    s = JMINUS.sub(' -', s)
    s = JPLUS.sub(' +', s)
    return SPACED_MINUS.sub(lambda m: m.group(1) + '-', s)


def values_v(seg: str, o: dict) -> list:
    s = signed(seg) if o.get('minus') else seg
    out = A.values(s)
    if o.get('factored'):
        for m in FACTORED.finditer(s.replace('−', '-')):
            k = int(m.group(2))
            out += [(v * 10.0 ** k, u * 10.0 ** k) for v, u in A.values(m.group(1))]
    return out


def family_of(w: str) -> tuple:
    if w in A.LINEAR:
        return A.LINEAR
    return next((f for f in A.FAMILIES if w in f), tuple(sorted(A.VERDICT_WORDS)))


def mentions(low: str, vocab: tuple, revised: bool) -> list:
    """The non-negated mentions of a family's words, in order. `revised` also skips 'no X', 'incrementally X' and a
    word inside a criteria span ('between laminar and turbulent', 'laminar or turbulent')."""
    spans = []
    if revised:
        for m in LIST_SPAN.finditer(low):
            pair = [g for g in m.groups() if g]
            if len(pair) == 2 and all(p in vocab for p in pair):
                spans.append(m.span())
    out = []
    for m in re.finditer(r'\b(%s)\b' % '|'.join(map(re.escape, vocab)), low):
        if A.NEGATED.search(low[max(0, m.start() - 24):m.start()]):
            continue
        if revised and (re.search(r'\bno\s+$', low[max(0, m.start() - 4):m.start()])
                        or re.search(r'\bincrementally\s+$', low[max(0, m.start() - 15):m.start()])
                        or any(a <= m.start() < b for a, b in spans)):
            continue
        out.append(m.group(1))
    return out


def word_hit(low: str, w: str, rule: str) -> bool:
    if rule == 'last':
        return A._word_hit(low, w)
    found = mentions(low, family_of(w), revised=rule == 'revised')
    if rule == 'first':
        return bool(found) and found[0] == w
    if not found and w in A.LINEAR:          # 'Answer: Yes' to 'determine if this system is linear'
        yn = re.findall(r'\b(yes|no)\b', low)
        return bool(yn) and (yn[-1] == 'yes') == (w == 'linear')
    return bool(found) and found[-1] == w


# ------------------------------------------------------------------ the label parts LABEL_PARTS names

def _norm_species(s: str) -> str:
    """'C3H8(g)', 'C₃H₈', '$\\text{O}_2$' and '\\mathrm{NH_3}' all to 'c3h8', 'o2', 'nh3'."""
    s = s.translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789'))
    s = re.sub(r'\\[()\[\]]|\\(?:text|mathrm|rm|ce)\b', '', s)
    s = re.sub(r'\((?:g|l|aq|s)\)|[\s*_${}\\]', '', s).lower()
    return re.sub(r'^[^a-z0-9]+', '', s)


UNIT_EXPONENT = re.compile(r'(?<=[a-z])\s*(?:\^\s*\{?\s*[-+]?\d+\s*\}?|[\u207b\u207a]?[\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+)')
LATEX_DELIMS = re.compile(r'\\[()\[\]]|\$|~|\\[,;:! ]')
HAT = re.compile(r'\\(?:hat|vec)\s*\{\s*(?:\\mathbf\s*)?\{?\s*([xyz])\s*\}?\s*\}|\\mathbf\s*\{?\s*([xyz])\s*\}?|\\hat\s*([xyz])')


def _family_last(low: str, words: tuple) -> str | None:
    found = mentions(low, words, revised=True)
    return found[-1] if found else None


def label_targets(item: dict, want_numbers: list) -> list:
    """(name, gold label, reader) for the label part of a LABEL_PARTS template, from the gold's answer segment. A label
    the answer implies counts as stated: a total energy in J or an average power in W names the signal type, and a
    stated vorticity names the flow rotational (nonzero) or irrotational (zero). The readers were checked on the main
    store: every correct-answer response they call a miss, of those read, states no label in its final answer."""
    tid, gseg = item['template_id'], ' '.join(A.segment(item['solution']).split())
    glow = gseg.lower()
    if tid == 'template_decimation_aliasing_analysis':
        gold = 'no' if re.search(r'did not|\*\*not\*\*|no aliasing', glow) else 'yes'

        def aliasing(low, g=gold):
            t = low.replace('*', '')
            # a "no" statement first: "no aliasing occurred" also contains "aliasing occurred"
            if re.search(r'no aliasing|aliasing (?:does|did|will|would) not|not alias|without aliasing|'
                         r'aliasing\W{0,6}no\b|aliasing is avoided|avoids aliasing|\bb\)\s*(?:answer:\s*)?no\b', t):
                said = 'no'
            elif re.search(r'aliasing (?:occur|is present|happen|takes place|does occur|did occur)|\baliased\b|'
                           r'aliasing\W{0,6}yes\b|aliasing has occurred|there is aliasing|aliasing will occur|'
                           r'\bb\)\s*(?:answer:\s*)?yes\b', t):
                said = 'yes'
            else:
                said = None
            return said == g
        return [('aliasing', gold, aliasing)]
    if tid == 'template_limiting_reactant':
        m = re.search(r'limiting reactant is \*\*([^*]+)\*\*', gseg)
        if not m:
            return []
        gold = _norm_species(m.group(1))

        def limiting(low, g=gold):
            after = [_norm_species(x.group(1)) for x in re.finditer(
                r'limiting (?:reactant|reagent)[\s:*=]{0,8}(?:is\b)?\s*(?:the\s+)?(.{1,40})', low)]
            before = [_norm_species(x.group(1)) for x in re.finditer(r'(\S{1,30})\s+is\s+(?:the\s+)?limiting', low)]
            return any(x.startswith(g) for x in after) or any(x.endswith(g) for x in before)
        return [('limiting reactant', gold, limiting)]
    if tid == 'template_signal_energy_power':
        gold = 'energy' if 'energy signal' in glow else 'power'

        def signal_type(low, g=gold):
            said = _family_last(low.replace(' signal', '_signal'), ('energy_signal', 'power_signal'))
            if said is None:            # a total energy in J or an average power in W names the type by its unit
                t = re.sub(r'\\(?:text|mathrm|rm)\s*\{|[{}]', ' ', LATEX_DELIMS.sub(' ', low))
                j = re.search(r'\d\s*(?:j|joules?)\b|\btotal energy\b', t)
                w = re.search(r'\d\s*(?:w|watts?)\b|\baverage power\b', t)
                if j and w:             # both stated: the one named first is the answer's
                    said = 'energy' if j.start() < w.start() else 'power'
                else:
                    said = 'energy' if j else ('power' if w else None)
            return bool(said) and said.startswith(g)
        return [('signal type', gold, signal_type)]
    if tid == 'template_system_properties':
        gold = next((w for w in ('underdamped', 'critically', 'overdamped') if w in glow), None)
        return [('system type', gold, lambda low, g=gold: _family_last(low, ('underdamped', 'critically', 'overdamped',
                                                                           'undamped')) == g)] if gold else []
    if tid == 'template_vorticity_check':
        gold = 'irrotational' if 'irrotational' in glow else 'rotational'
        given = M.numbers(item['question'])

        def rotationality(low, g=gold, given=given):
            said = _family_last(low, ('rotational', 'irrotational'))
            if said is None:            # a stated vorticity says it: any nonzero component is rotational
                vs = [v for v, _u in A.values(UNIT_EXPONENT.sub(' ', signed(low)))]
                fresh = [v for v in vs if v and not any(M.close(v, x, 1e-9) for x in given)]
                if fresh or (vs and not any(v == 0 for v in vs)):
                    said = 'rotational'
                elif vs:
                    said = 'irrotational'
            return said == g
        return [('rotationality', gold, rotationality)]
    if tid == 'template_wave_equation_interpretation':
        m = re.search(r'(positive|negative)\s+([xyz])-direction', glow)
        if not m:
            return []
        gold = ('+' if m.group(1) == 'positive' else '-') + m.group(2)

        def direction(low, g=gold):
            low = LATEX_DELIMS.sub(' ', low.replace('−', '-').replace('–', '-'))
            low = HAT.sub(lambda m: next(x for x in m.groups() if x), low)
            found, named = [], []
            for x in re.finditer(r'(positive|negative|\+|-)\s*\$?\s*\\?(?:hat\{)?\s*([xyz])\b|'
                                 r'\b([xyz])\s*-?\s*direction\s*\(?(positive|negative)', low):
                if x.group(1):
                    d = ('+' if x.group(1) in ('positive', '+') else '-') + x.group(2)
                else:
                    d = ('+' if x.group(4) == 'positive' else '-') + x.group(3)
                found.append(d)
                if re.search(r'direction|propagat', low[max(0, x.start() - 30):x.end() + 30]):
                    named.append(d)
            # the direction stated beside "direction" or "propagation" first (part b), else the last one stated
            pick = named[0] if named else (found[-1] if found else None)
            return pick == g
        return [('direction', gold, direction)]
    return []


# ------------------------------------------------------------------ units written beside a number

PREFIX = {'p': 1e-12, 'n': 1e-9, 'u': 1e-6, 'µ': 1e-6, 'μ': 1e-6, 'm': 1e-3, 'c': 1e-2, 'k': 1e3, 'M': 1e6, 'G': 1e9}
TIME = {'s': 1.0, 'sec': 1.0, 'second': 1.0, 'seconds': 1.0, 'min': 60.0, 'mins': 60.0, 'minute': 60.0,
        'minutes': 60.0, 'h': 3600.0, 'hr': 3600.0, 'hrs': 3600.0, 'hour': 3600.0, 'hours': 3600.0}
BASES = ('Pa', 'J', 'W', 'N', 'Hz', 'F', 'C', 'V', 'A', 'g', 'L', 'm', 'rad', 'Wb', 'T', 'Ω', 'ohm', 'S', 'psi', 'lb')
NUMU = re.compile(r'([-+]?\d[\d,]*\.?\d*(?:\s*(?:[eE][-+]?\d+|(?:\\times|×|x)\s*10\s*(?:\^\s*\{?\s*[-+]?\d+\s*\}?|'
                  r'[⁻⁺]?[⁰¹²³⁴-⁹]+)))?)'
                  r'(?:\s|\\[,;:! ]|~|\$|\*)*(?:\\(?:text|mathrm|rm)\s*\{|\{)?\s*([%A-Za-zµμΩ°][A-Za-zµμΩ°]*)')


def unit_scale(tok: str):
    t = tok.strip('.')
    if t in TIME:
        return ('time', TIME[t])
    if t == '%':
        return ('fraction', 0.01)
    if t == 'ppm':
        return ('fraction', 1e-6)
    if t == 'ksi':
        return ('psi', 1e3)
    if t in ('kip', 'kips'):
        return ('lb', 1e3)
    for b in sorted(BASES, key=len, reverse=True):
        if t == b:
            return (b, 1.0)
        if t.endswith(b) and len(t) == len(b) + 1 and t[0] in PREFIX:
            return (b, PREFIX[t[0]])
    return None


def units_at(seg: str, value: float) -> list:
    """The units written after the numbers in `seg` that equal `value` (to 1e-6 relative): parsed (base, factor) or
    None where a unit is written that is not parsed; an empty list where no number equals `value`."""
    out = []
    for m in NUMU.finditer(signed(seg)):
        if any(math.isclose(v, value, rel_tol=1e-6, abs_tol=1e-300) for v, _u in A.values(m.group(1))):
            out.append(unit_scale(m.group(2)))
    return out


def unit_explained(gold: float, gseg: str, hits: list, rseg: str) -> bool:
    """A match under a unit factor c != 1 stands unless every response number that made it carries a parsed unit of
    the gold's own base and prefix (a unit slip: the number is off by c in the gold's unit)."""
    g = next((u for u in units_at(gseg, abs(gold)) + units_at(gseg, gold) if u), None)
    if g is None:
        return True
    for v, c in hits:
        us = units_at(rseg, v)
        if not us or any(u is None for u in us):
            return True
        for u in us:
            if u[0] != g[0] or math.isclose(u[1] / g[1], abs(c), rel_tol=1e-6):
                return True
    return False


# ------------------------------------------------------------------ targets

NUMLIT = re.compile(r'(?<![A-Za-z_\d.^{\\])-?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][-+]?\d+)?')
INDEX_LIST = re.compile(r'for n\s*=\s*([-\d,\s]+)')


def stated_parts(item: dict, ms: tuple, have: list) -> list:
    """Values the gold's answer segment states that a full answer would contain and that are not already targets: not a
    restated input (in any unit factor), not 0, 1 or 2, not a digit of a formula or a subscript, not part of a
    pi-fraction, not an index of a listed sequence."""
    gseg = signed(A.segment(item['solution']))
    given = M.numbers(item['question'])
    idx = {int(x) for m in INDEX_LIST.finditer(gseg) for x in re.findall(r'-?\d+', m.group(1))}
    out = []
    for m in NUMLIT.finditer(gseg):
        lit, before, after = m.group(0), gseg[max(0, m.start() - 8):m.start()], gseg[m.end():m.end() + 6]
        if re.match(r'\s*\*?\s*(?:pi|π)', after) or re.search(r'(?:pi|π)\)?\s*/\s*\(?$', before):
            continue
        vs = A.values(lit)
        if not vs:
            continue
        v, u = vs[0]
        if any(abs(v - k) < 1e-9 for k in (0.0, 1.0, 2.0)) or (float(v).is_integer() and int(v) in idx):
            continue
        if any(M.close(v, g, 1e-9) for g in given) or M.scaled_match(abs(v), tuple(abs(g) for g in given),
                                                                    M.DISPLAY_TOL) is not None:
            continue
        if any(M.close(abs(v), abs(y), 1e-6) for y, _ in have + out):
            continue
        out.append((v, u))
    return out


def targets_v(item: dict, ms: tuple, o: dict) -> dict:
    kind = item.get('answer_type')
    it = item
    if o.get('minus'):
        it = dict(it, solution=signed(it['solution']))
    if o.get('array_all') and kind == 'array':
        it = dict(it, answer_type='vector')
    want = A.targets(it, ms)
    nums = list(want['numbers'])
    if o.get('array_all') and kind == 'array':
        gseg = A.segment(item['solution'])
        idx = {float(x) for m in INDEX_LIST.finditer(gseg) for x in re.findall(r'-?\d+', m.group(1))}
        seq = {float(x) for m in re.finditer(r'\{([^{}]*)\}', gseg) for x in re.findall(r'-?\d+(?:\.\d+)?', m.group(1))}
        nums = [(y, u) for y, u in nums if not (y in idx and y not in seq)] or nums
    if o.get('parts_all') and kind in ('vector', 'multipart', 'array', 'classification'):
        nums = nums + stated_parts(item, ms, nums)
    if o.get('vacuous'):
        nums = [(y, 0.0 if abs(y) <= u else u) for y, u in nums]
    extra = label_targets(item, nums) if o.get('label_parts') and item['template_id'] in LABEL_PARTS else []
    return dict(want, numbers=nums, extra_labels=extra)


# ------------------------------------------------------------------ the rule, with switches

def _ok(t, tu, v, u, gold, gold_ulp, exact, unit):
    if exact:
        return abs(t - gold) <= 1e-9 * max(1.0, abs(gold))
    own = min(unit * tu, A.OWN_DIGIT_CAP * abs(gold)) if abs(v) > u else 0.0
    return abs(t - gold) <= max(A.REL * abs(gold), own, unit * gold_ulp) * (1 + A.SLACK)


def path_of(gold, have, gold_ulp, exact, unit, scales):
    """Every (value, unit factor, path) that matches `gold`; path: eps, digit, exact or unit."""
    out = []
    for v, u in have:
        for sc in scales:
            t, tu = v * sc, u * sc
            if _ok(t, tu, v, u, gold, gold_ulp, exact, unit):
                if sc != 1.0:
                    p = 'unit'
                elif exact:
                    p = 'exact'
                else:
                    p = 'eps' if abs(t - gold) <= A.REL * abs(gold) * (1 + A.SLACK) else 'digit'
                out.append((v, sc, p))
    return out


def target_hit(y, gu, want, have, o, item, seg, gseg) -> str | None:
    """How target y is matched (eps, digit, exact, unit, sign, rendering), or None."""
    unit = o.get('unit_window', 1.0)
    scales = (1.0,) if o.get('units') == 'none' else A.SCALES
    exact = want['exact']
    absolute = [(abs(v), u) for v, u in have]
    mode = o.get('abs', 'on')
    sign_ok = mode == 'on' or (mode == 'policy' and (item['template_id'] in SIGN_CONVENTION or bool(
        SIGN_WORDS[-1 if y < 0 else 1].search(seg))))
    best = None
    for gold, pool, tag in ((y, have, None), (abs(y), absolute, 'sign')):
        if tag == 'sign' and not sign_ok:
            continue
        hits = path_of(gold, pool, gu, exact, unit, scales)
        if not hits:
            continue
        direct = [h for h in hits if h[2] != 'unit']
        if direct:
            p = min((h[2] for h in direct), key=['eps', 'exact', 'digit'].index)
        elif o.get('units') == 'guard' and not unit_explained(gold, gseg, [(h[0], h[1]) for h in hits], seg):
            continue
        else:
            p = 'unit'
        p = tag or p
        if best is None:
            best = p
        if tag is None:
            break
    if best is None:
        for w, wu in want.get('renderings', []):
            if path_of(w, have, wu, exact, unit, (1.0,)) or (sign_ok and path_of(abs(w), absolute, wu, exact, unit, (1.0,))):
                return 'rendering'
    return best


def verdict(text: str, item: dict, ms: tuple, o: dict | None = None, paths: list | None = None) -> str:
    """answer.verdict at the fitted tolerance (unit 1.0, whole False) with the switches of VARIANTS."""
    o = o or {}
    want = targets_v(item, ms, o)
    seg = A.segment(text)
    gseg = A.segment(item['solution'])
    have = values_v(seg, o)
    low = A._norm_linear(seg).lower()
    hits = []
    for y, gu in want['numbers']:
        p = target_hit(y, gu, want, have, o, item, seg, gseg)
        hits.append(p is not None)
        if paths is not None:
            paths.append(p)
    for w in want['words']:
        hits.append(word_hit(low, w, o.get('label', 'last')))
    st = NO_UNITS.sub(' ', low) if o.get('no_units') else low
    for lab, val in want.get('labeled', []):
        hits.append(A._stance(st, lab) == val)
    for _name, _gold, reader in want.get('extra_labels', []):
        hits.append(bool(reader(low)))
    if not hits:
        label = 'incorrect'
    else:
        n = sum(hits)
        label = 'correct' if n == len(hits) else ('partial' if n else 'incorrect')
    if item['template_id'] in A.SYMBOLIC_EQUIVALENCE_TEMPLATES:
        if label != 'correct':
            label = A.symbolic_equivalence.apply(label, item, text, A.SYMBOLIC_EQUIVALENCE_TEMPLATES, None, 1.0, A.segment)[0]
        elif o.get('symbolic') == 'decide':
            eq = A.symbolic_equivalence.apply('incorrect', item, text, A.SYMBOLIC_EQUIVALENCE_TEMPLATES, None, 1.0,
                                              A.segment)[0]
            if eq != 'correct':
                label = 'incorrect'
    return label


READER = {'minus': True, 'factored': True, 'no_units': True}
TARGETS = {'label': 'revised', 'array_all': True, 'vacuous': True, 'parts_all': True, 'label_parts': True}
VARIANTS = {
    'headline': {},
    'R1 minus signs read': {'minus': True},
    'R2 factored power of ten read': {'factored': True},
    'R3 "no units" not read as a stance': {'no_units': True},
    'R4 label word: first mention': {'label': 'first'},
    'R5 label word: revised last mention': {'label': 'revised'},
    'T1 arrays: every computed value': {'array_all': True},
    'T2 targets near zero: no gold last-digit window': {'vacuous': True},
    'T3 every stated part a target': {'minus': True, 'parts_all': True},
    'T4 label parts checked': {'label_parts': True},
    'U1 no unit factors': {'units': 'none'},
    'U2 unit factor only where the written unit allows it': {'units': 'guard'},
    'S1 absolute-value clause off, signs read': {'minus': True, 'abs': 'off'},
    'S2 sign policy, signs read': {'minus': True, 'abs': 'policy'},
    'D1 half-unit last-digit windows': {'unit_window': 0.5},
    'Y1 equivalence decides on its five templates': {'symbolic': 'decide'},
    'C1 reader and target fixes (R1-R3, R5, T1-T4)': {**READER, **TARGETS},
    'C2 C1 with U2 and S2': {**READER, **TARGETS, 'units': 'guard', 'abs': 'policy'},
    'C3 C2 with D1 and Y1': {**READER, **TARGETS, 'units': 'guard', 'abs': 'policy', 'unit_window': 0.5,
                             'symbolic': 'decide'},
}


# ------------------------------------------------------------------ the populations

def study_population() -> list[dict]:
    from full_run_28092026 import validate_scorer as V
    truth, keyfile, items, texts = V.pilot_inputs()
    out = []
    for code, t in sorted(truth.items()):
        key = keyfile[code]
        if key[1] not in items or key not in texts:
            continue
        it = items[key[1]]
        out.append({'model': key[0], 'item_id': key[1], 'item': it, 'text': texts[key], 'kind': it.get('answer_type'),
                    'ms': tuple(m['value'] for m in t['milestones']), 'experts': t['final_answer']})
    return out


def readings_population(items: dict, ms: dict) -> list[dict]:
    keyfile, queues, rows = K.read_returns(K.RETURNED, K.OUT)
    shown = K.shown_pool(K.OUT)
    models = sorted({k['model'] for k in keyfile.values() if 'model' in k})
    store = K.store_rows(models)
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        by_code[c][aid] = r
    texts, out = {}, []
    for c, readers in sorted(by_code.items()):
        k = keyfile[c]
        if k['kind'] != 'answer' or K.reading_status(k, shown.get(c), items, store) != K.KEPT:
            continue
        r = store[k['model']][k['item_id']]
        if k['model'] not in texts:
            texts[k['model']] = score.texts_matching('main', k['model'], list(store[k['model']].values()))
        out.append({'model': k['model'], 'item_id': k['item_id'], 'item': items[k['item_id']], 'kind': r['answer_type'],
                    'text': texts[k['model']][k['item_id']], 'ms': tuple(m['value'] for m in ms[k['item_id']]),
                    'stored': r['answer']['label'],
                    'readings': [K.EXPERT_TO_CHECK[x['answers']['verdict']] for _a, x in sorted(readers.items())],
                    'notes': [x['note'] for _a, x in sorted(readers.items()) if x['note']]})
    return out


# ------------------------------------------------------------------ the headline responses, one model per worker

def store_of(key: str) -> str:
    """The store a model's headline responses are in: matched settings (D2), each model at its reasoning store where
    it has one (results/matched_config.json, as clause_variants.py and analyze.py read it), else at main."""
    return CV.matched_config().get(key) or 'main'


def score_model(args) -> dict:
    key, items, ms = args
    store = store_of(key)
    rows = CV.load_rows(store, key)
    texts = score.texts_matching(store, key, rows)
    out = {'key': key, 'fac': {v: 0.0 for v in VARIANTS}, 'n': len(rows), 'moves': {v: [] for v in VARIANTS},
           'paths': collections.Counter(), 'bad': []}
    for r in rows:
        if r['unusable']:
            continue
        it, text = items[r['item_id']], texts[r['item_id']]
        vals = tuple(m['value'] for m in ms[r['item_id']])
        paths = []
        head = verdict(text, it, vals, {}, paths)
        if head != r['answer']['label']:
            out['bad'].append((r['item_id'], head, r['answer']['label']))
            continue
        if head != 'incorrect':
            ps = {p for p in paths if p}
            worst = max(ps, key=['eps', 'exact', 'digit', 'rendering', 'unit', 'sign'].index) if ps else 'words or labels'
            out['paths'][(r['answer_type'], worst)] += 1
        for v, o in VARIANTS.items():
            lab = head if not o else verdict(text, it, vals, o)
            out['fac'][v] += SCORE_OF[lab]
            if lab != head:
                out['moves'][v].append((r['item_id'], r['template_id'], r['answer_type'], head, lab))
    out['fac'] = {v: s / len(rows) for v, s in out['fac'].items()}
    out['paths'] = {f'{a}|{b}': n for (a, b), n in out['paths'].items()}
    return out


# ------------------------------------------------------------------ part 3: the targets, item by item

def target_coverage(items: dict, ms: dict) -> dict:
    flags = collections.defaultdict(lambda: collections.defaultdict(set))
    reasons = collections.defaultdict(collections.Counter)
    for iid, it in items.items():
        kind, tid = it['answer_type'], it['template_id']
        vals = tuple(m['value'] for m in ms[iid])
        want = A.targets(it, vals)
        nums = want['numbers']
        if any(abs(y) <= u for y, u in nums):
            flags['a target lies within its own last-digit window of zero'][kind].add(iid)
        gseg = A.segment(it['solution'])
        given = M.numbers(it['question'])
        stated = [(v, u) for v, u in A.values(gseg) if not any(M.close(v, g, 1e-9) for g in given)]
        if kind == 'array' and len({round(v, 9) for v, _ in stated}) > 1:
            flags['an array scored on its last value'][kind].add(iid)
            if INDEX_LIST.search(gseg):
                flags['an array whose gold line lists indices beside the values'][kind].add(iid)
        if kind in ('vector', 'multipart', 'classification', 'array'):
            unsigned = [v for v, _u in A.values(gseg)]
            for v, _u in stated_parts(it, vals, nums):
                if v < 0 and not any(M.close(v, x, 1e-9) for x in unsigned) and any(M.close(-v, x, 1e-9) for x in unsigned):
                    why = 'a minus the reader does not join to its number'
                elif M.scaled_match(v, vals, M.DISPLAY_TOL) is not None:
                    why = ('a computed value the last-value rule drops' if kind == 'array'
                           else 'a computed value the reader keeps once (a repeated value)')
                elif M.scaled_match(-v, vals, M.DISPLAY_TOL) is not None:
                    why = 'a magnitude with its sign in words (the computed value is negative)'
                else:
                    why = 'a value not among the computed values'
                flags['a stated part that is not a target: ' + why][kind].add(iid)
                reasons[why][tid] += 1
        if tid in LABEL_PARTS:
            flags['a label part the question asks for, not checked'][kind].add(iid)
        if tid in REPORTED_VALUE:
            flags['a question naming one value to report among other requested results'][kind].add(iid)
    out = {}
    for f, by_kind in flags.items():
        out[f] = {'items': sum(len(s) for s in by_kind.values()),
                  'templates': len({items[i]['template_id'] for s in by_kind.values() for i in s}),
                  'by_kind': {k: [len(s), len({items[i]['template_id'] for i in s})] for k, s in by_kind.items()},
                  'template_list': sorted({items[i]['template_id'].replace('template_', '') for s in by_kind.values()
                                           for i in s})}
    out['_reasons_by_template'] = {w: dict(c.most_common()) for w, c in reasons.items()}
    return out


# ------------------------------------------------------------------ agreement

def agreement(pop: list[dict], labels: list[str], study: bool) -> dict:
    by_kind = collections.defaultdict(lambda: [0, 0, 0, 0])      # same, n, non-partial same, non-partial n
    for p, lab in zip(pop, labels):
        refs = [p['experts']] if study else p['readings']
        for e in refs:
            b = by_kind[p['kind']]
            b[0] += lab == e
            b[1] += 1
            if e != 'partial':
                b[2] += (lab == 'correct') == (e == 'correct')
                b[3] += 1
    tot = [sum(b[i] for b in by_kind.values()) for i in range(4)]
    return {'all': tot, 'by_kind': {k: by_kind[k] for k in KINDS if k in by_kind}}


def confusion(pop: list[dict], labels: list[str], study: bool) -> dict:
    out = collections.defaultdict(collections.Counter)
    for p, lab in zip(pop, labels):
        for e in ([p['experts']] if study else p['readings']):
            out[p['kind']][f'{lab}>{e}'] += 1
    return {k: dict(out[k]) for k in KINDS if k in out}


# ------------------------------------------------------------------ part 2: causes

DIRECTION = {
    'answer_parts_unchecked': 'both', 'unit_factor_slip': 'lenient', 'last_digit_window': 'lenient',
    'tolerance': 'lenient', 'value_in_other_role': 'lenient', 'contradictory_answer': 'lenient',
    'reader_minus_sign': 'strict', 'reader_no_units_stance': 'strict', 'reader_factored_power': 'strict',
    'label_rule': 'strict', 'equivalent_form': 'strict', 'array_target_index': 'strict',
    'array_last_value_only': 'strict', 'gold_rounding': 'split', 'expert_partial_judgment': 'split'}
CAUSE_TEXT = {
    'answer_parts_unchecked': 'the question asks for results beyond the value its gold line states; the check scores '
                              'that value, the experts the whole answer',
    'unit_factor_slip': 'a value off by a unit factor in the gold\'s own unit, matched through the factor',
    'last_digit_window': 'a value one unit off in its last digit, accepted by the last-digit window (not a correct '
                         'rounding)',
    'tolerance': 'a value within 0.2% of the target that the experts hold to the gold\'s printed digits',
    'value_in_other_role': 'a stated number matched to a target it does not answer (another component or part)',
    'contradictory_answer': 'an answer that states two different values for one quantity, credited on one of them',
    'reader_minus_sign': 'a minus sign the reader does not read (U+2212 in an exponent)',
    'reader_no_units_stance': '"No units." after a yes/no part read as the stance "no"',
    'reader_factored_power': 'a power of ten written once after a parenthesised vector, not applied to its components',
    'label_rule': 'the last mention of a label family decides: a hedge, a negation with "no", criteria restated '
                  'after the verdict, "incrementally linear", or "Yes" for a linearity question',
    'equivalent_form': 'an equivalent form the number rule cannot read (an angle reduced modulo 2pi, an '
                       'amplitude-phase form)',
    'array_target_index': 'the array\'s target is the last number of the gold line, an index of the sequence',
    'array_last_value_only': 'the array is scored on its last value; the other value is right',
    'gold_rounding': 'the gold carries rounding from rounded intermediates; the two experts split',
    'expert_partial_judgment': 'the experts give partial credit (or split) where the check scores the parts as stated',
}


def causes_table(study: list[dict], s_labels: list[str], readings: list[dict], r_labels: list[str]) -> dict:
    if not CAUSES.exists():
        return {'available': False}
    c = json.loads(CAUSES.read_text(encoding='utf-8'))
    out = {'available': True, 'study': collections.Counter(), 'readings': collections.Counter(),
           'study_by_kind': collections.defaultdict(collections.Counter),
           'readings_by_kind': collections.defaultdict(collections.Counter), 'cases': []}
    seen = {'study': set(), 'readings': set()}
    for p, lab in zip(study, s_labels):
        k = f"{p['model']}@{p['item_id']}"
        if lab == p['experts']:
            if k in c['study']:
                raise SystemExit(f'answer_audit_causes.json gives a cause to an agreement: study {k}')
            continue
        if k not in c['study']:
            raise SystemExit(f'answer_audit_causes.json has no cause for the study disagreement {k}')
        seen['study'].add(k)
        out['study'][c['study'][k]] += 1
        out['study_by_kind'][p['kind']][c['study'][k]] += 1
        out['cases'].append(('study', k, p['kind'], lab, [p['experts']], c['study'][k], [], p))
    for p, lab in zip(readings, r_labels):
        k = f"{p['model']}@{p['item_id']}"
        diff = [e for e in p['readings'] if e != lab]
        if not diff:
            if k in c['readings']:
                raise SystemExit(f'answer_audit_causes.json gives a cause to an agreement: readings {k}')
            continue
        if k not in c['readings']:
            raise SystemExit(f'answer_audit_causes.json has no cause for the readings disagreement {k}')
        seen['readings'].add(k)
        out['readings'][c['readings'][k]] += len(diff)
        out['readings_by_kind'][p['kind']][c['readings'][k]] += len(diff)
        out['cases'].append(('readings', k, p['kind'], lab, p['readings'], c['readings'][k], p['notes'], p))
    for pop in ('study', 'readings'):
        extra = set(c[pop]) - seen[pop]
        if extra:
            raise SystemExit(f'answer_audit_causes.json names {pop} cases that are not disagreements: {sorted(extra)}')
    return out


# ------------------------------------------------------------------ the run

def tau(head: list, other: list) -> float:
    return float(stats.kendalltau(head, other).statistic)


def run(workers: int) -> dict:
    items = score.pool_items()
    ms = CV.milestone_cache(items)
    study = study_population()
    readings = readings_population(items, ms)
    # the populations under every variant
    pops = {}
    for name, pop, is_study in (('study', study, True), ('readings', readings, False)):
        res = {}
        for v, o in VARIANTS.items():
            labs = [verdict(p['text'], p['item'], p['ms'], o) for p in pop]
            if v == 'headline':
                ref = [A.verdict(p['text'], p['item'], None, p['ms'])[0] for p in pop] if is_study else \
                    [p['stored'] for p in pop]
                bad = [(p['model'], p['item_id']) for p, a, b in zip(pop, labs, ref) if a != b]
                if bad:
                    raise SystemExit(f'{name}: the local rule differs from the check on {len(bad)} responses: {bad[:5]}')
            res[v] = {'labels': labs, 'agreement': agreement(pop, labs, is_study)}
        pops[name] = res
    s3 = pops['study']['headline']['agreement']['all']
    r3 = pops['readings']['headline']['agreement']['all']
    if (round(s3[0] / s3[1], 3), round(s3[2] / s3[3], 3)) != (0.930, 0.986) or (r3[0], r3[1]) != (244, 288):
        raise SystemExit(f'the headline agreement does not reproduce 0.930, 0.986 and 244 of 288: {s3}, {r3}')
    causes = causes_table(study, pops['study']['headline']['labels'], readings, pops['readings']['headline']['labels'])
    # the headline responses (matched settings)
    done = {}
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for res in pool.map(score_model, [(k, items, ms) for k in ROSTER]):
            if res['bad']:
                raise SystemExit(f"main/{res['key']}: the local rule differs from the stored verdict on {len(res['bad'])} "
                                 f"rows (first: {res['bad'][:3]})")
            done[res['key']] = res
            print(f"{store_of(res['key'])}/{res['key']}: every readable verdict reproduced", flush=True)
    keys = list(ROSTER)
    head = [done[k]['fac']['headline'] for k in keys]
    main = {}
    for v in VARIANTS:
        moves = [m for k in keys for m in done[k]['moves'][v]]
        up = sum(ORDER[m[4]] > ORDER[m[3]] for m in moves)
        by_kind = collections.Counter((m[2], 'up' if ORDER[m[4]] > ORDER[m[3]] else 'down') for m in moves)
        tpl = collections.Counter(m[1].replace('template_', '') for m in moves)
        fac = {k: done[k]['fac'][v] for k in keys}
        main[v] = {'moved': len(moves), 'up': up, 'down': len(moves) - up,
                   'by_kind': {f'{a} {b}': n for (a, b), n in sorted(by_kind.items())},
                   'templates': len(tpl), 'top_templates': dict(tpl.most_common(8)), 'fac': fac,
                   'delta': {k: fac[k] - done[k]['fac']['headline'] for k in keys},
                   'tau': tau(head, [fac[k] for k in keys]) if v != 'headline' else 1.0}
    paths = collections.Counter()
    for k in keys:
        paths.update(done[k]['paths'])
    return {'items': items, 'study': study, 'readings': readings, 'pops': pops, 'causes': causes, 'main': main,
            'paths': dict(paths), 'coverage': target_coverage(items, ms),
            'changes': {v: [(k,) + m for k in keys for m in done[k]['moves'][v]] for v in VARIANTS if v != 'headline'}}


# ------------------------------------------------------------------ the report

def pct3(a, b):
    return f'{a / b:.3f}' if b else '-'


def report(res: dict) -> list[str]:
    P, main, cov = res['pops'], res['main'], res['coverage']
    L = ['# The answer check, audited by answer kind', '',
         'Generated by `answer_audit.py`; what it measures is defined in its docstring. Counts only: the experts\' '
         'labels, readings and notes stay local, and the cases are in `expert_request/ANSWER_AUDIT_CASES.md`.', '']
    # 1
    L += ['## 1. Agreement with the experts by answer kind', '',
          'The study: 300 responses from five LLMs outside the evaluated eleven, three experts each (three-way '
          'agreement, and on the responses the experts did not call partial). The readings: two domain experts per '
          'final-answer verdict on the five strongest evaluated models, against the check\'s current verdict.', '',
          '| answer kind | study responses | three-way | not partial | readings | agreeing |', '|---|---:|---:|---:|---:|---:|']
    sa, ra = P['study']['headline']['agreement'], P['readings']['headline']['agreement']
    for k in KINDS:
        s = sa['by_kind'].get(k, [0, 0, 0, 0])
        r = ra['by_kind'].get(k, [0, 0, 0, 0])
        L.append(f'| {k} | {s[1]} | {pct3(s[0], s[1])} | {pct3(s[2], s[3])} | {r[1]} | {r[0]} ({pct3(r[0], r[1])}) |')
    s, r = sa['all'], ra['all']
    L.append(f'| all | {s[1]} | {pct3(s[0], s[1])} | {pct3(s[2], s[3])} | {r[1]} | {r[0]} ({pct3(r[0], r[1])}) |')
    L += ['', 'Check verdict > experts\' verdict, by kind:', '']
    for name, is_study in (('study', True), ('readings', False)):
        conf = confusion(res[name], P[name]['headline']['labels'], is_study)
        L.append(f'- {name}: ' + '; '.join(f'{k} ' + ', '.join(f'{a} {n}' for a, n in sorted(c.items()) if a.split('>')[0]
                                                                != a.split('>')[1]) for k, c in conf.items()
                                            if any(a.split('>')[0] != a.split('>')[1] for a in c)))
    # 2
    c = res['causes']
    L += ['', '## 2. Every disagreement, by cause', '']
    if not c.get('available'):
        L.append('`expert_request/answer_audit_causes.json` is not here: the causes are not counted.')
    else:
        L += ['Causes assigned by reading each response, the gold and, for the readings, the experts\' notes. '
              'Direction: "lenient" where the check credits more than the experts, "strict" where it credits less, '
              '"split" where the experts disagree with each other or give partial credit by judgment.', '',
              '| cause | direction | study responses | readings | what happens |', '|---|---|---:|---:|---|']
        for k in sorted(set(c['study']) | set(c['readings']), key=lambda x: -(c['study'][x] + c['readings'][x])):
            L.append(f"| {k} | {DIRECTION[k]} | {c['study'][k]} | {c['readings'][k]} | {CAUSE_TEXT[k]} |")
        L.append(f"| all | | {sum(c['study'].values())} | {sum(c['readings'].values())} | |")
        for pop in ('study', 'readings'):
            d = collections.Counter()
            for k, n in c[pop].items():
                d[DIRECTION[k]] += n
            L.append(f'\n{pop}: ' + ', '.join(f'{a} {n}' for a, n in d.most_common()))
    # 3
    L += ['', '## 3. What the check takes as targets (all 2,250 items)', '',
          '| finding | items | templates | by kind (items / templates) |', '|---|---:|---:|---|']
    for f, d in cov.items():
        if f.startswith('_'):
            continue
        L.append(f"| {f} | {d['items']} | {d['templates']} | " + '; '.join(f'{k} {a}/{b}' for k, (a, b) in d['by_kind'].items())
                 + ' |')
    L += ['', 'Templates per finding:', '']
    for f, d in cov.items():
        if not f.startswith('_'):
            L.append(f"- {f}: " + ', '.join(f'`{t}`' for t in d['template_list']))
    L += ['', 'Label parts not checked (`LABEL_PARTS`): ' + '; '.join(f"`{t.replace('template_', '')}` ({w})"
                                                                     for t, w in LABEL_PARTS.items()) + '.']
    # 4
    L += ['', '## 4. Variants of the rule', '',
          'Each variant changes one thing; C1 to C3 combine them (see `VARIANTS`). Study and readings: agreement with '
          'the experts (three-way; readings agreeing). Headline responses (matched settings): verdicts moved over the eleven models\' 24,750 '
          'responses, the largest change in a model\'s Final Answer Accuracy, and Kendall\'s tau with the ordering as '
          'scored.', '',
          '| variant | study three-way | study not partial | readings agreeing | main verdicts up / down | templates | '
          'largest FAC change | tau |', '|---|---:|---:|---:|---:|---:|---:|---:|']
    for v in VARIANTS:
        s, r, m = P['study'][v]['agreement']['all'], P['readings'][v]['agreement']['all'], main[v]
        big = max(m['delta'].values(), key=abs)
        L.append(f"| {v} | {pct3(s[0], s[1])} | {pct3(s[2], s[3])} | {r[0]} of {r[1]} | {m['up']} / {m['down']} | "
                 f"{m['templates']} | {big:+.3f} | {m['tau']:.3f} |")
    L += ['', 'Where each variant moves headline verdicts:', '']
    for v in VARIANTS:
        if v != 'headline' and main[v]['moved']:
            L.append(f"- {v}: by kind {main[v]['by_kind']}; templates {main[v]['top_templates']}")
    L += ['', 'Final Answer Accuracy per model under the combined variants:', '',
          '| model | as scored | ' + ' | '.join(v.split(' ')[0] for v in VARIANTS if v.startswith('C')) + ' |',
          '|---|---:|' + '---:|' * sum(v.startswith('C') for v in VARIANTS)]
    for k in sorted(main['headline']['fac'], key=lambda k: -main['headline']['fac'][k]):
        L.append(f"| `{k}` | {main['headline']['fac'][k]:.3f} | " +
                 ' | '.join(f"{main[v]['fac'][k]:.3f}" for v in VARIANTS if v.startswith('C')) + ' |')
    # 5
    L += ['', '## 5. How accepted headline answers are matched', '',
          'For each accepted (correct or partial) answer, the least direct way any of its numeric targets is matched: '
          'eps (within 0.2% at c = 1), exact (prescribed digits), digit (a last-digit window), rendering (a second '
          'unit on a scalar gold line), unit (a unit factor c != 1), sign (the absolute value).', '',
          '| answer kind | ' + ' | '.join(('eps', 'exact', 'digit', 'rendering', 'unit', 'sign', 'words or labels')) + ' |',
          '|---|' + '---:|' * 7]
    for k in KINDS:
        L.append(f'| {k} | ' + ' | '.join(str(res['paths'].get(f'{k}|{p}', 0)) for p in (
            'eps', 'exact', 'digit', 'rendering', 'unit', 'sign', 'words or labels')) + ' |')
    return L


def cases_lines(res: dict) -> list[str]:
    c = res['causes']
    L = ['# Answer-check audit: the cases (local)', '',
         'Written by `answer_audit.py`. Local: the experts\' verdicts and notes on named items. Not for the repository.', '']
    if not c.get('available'):
        return L + ['No causes file.']
    for pop, key, kind, lab, refs, cause, notes, p in sorted(c['cases'], key=lambda x: (x[0], x[5], x[1])):
        seg = ' '.join(A.segment(p['text']).split())[:500]
        gold = ' '.join(A.segment(p['item']['solution']).split())[:300]
        L += [f'## {pop}: {key} ({kind})', '', f'- check: {lab}; experts: {", ".join(refs)}; cause: {cause}',
              f'- gold: {gold}', f'- response: {seg}'] + [f'- note: {" ".join(n.split())}' for n in notes] + ['']
    return L


# ------------------------------------------------------------------ the paper's block and numbers

BLOCK = REPO / 'overleaf_source_04102026' / 'appendices' / 'answer_kinds.tex'
PAPER_JSON = HERE / 'results' / 'answer_audit_paper.json'
KIND_NAME = {'scalar': 'Scalar', 'multipart': 'Multipart', 'vector': 'Vector', 'array': 'Array', 'symbolic': 'Symbolic',
             'classification': 'Classification'}
CASES = {'lenient': ('muse-glimmer-30b', 'mm1k_finite_capacity#1'),
         'strict': ('deepseek-v4.1-flash', 'decimation_aliasing_analysis#15')}
BRANCH = {'chemical_engineering': 'chemical', 'civil_engineering': 'civil', 'electrical_engineering': 'electrical',
          'industrial_engineering': 'industrial', 'mechanical_engineering': 'mechanical'}


def ceil3(x: float) -> float:
    """The smallest three-decimal bound at or above |x|: 'by more than 0.008' must hold for the unrounded value."""
    return math.ceil(abs(x) * 1000 - 1e-9) / 1000


def paper_numbers(study: list, s_labels: list, readings: list, r_labels: list, causes: dict) -> dict:
    """Every number the paper's text and the block take from the audit: agreement by kind, the readings' discrepancies
    by verdict category and cause, the bounds of results/answer_audit.json, and the two cases, each fact asserted."""
    sa, ra = agreement(study, s_labels, True), agreement(readings, r_labels, False)
    cat = collections.defaultdict(collections.Counter)
    for pop, _key, _kind, lab, refs, cause, _notes, _p in causes['cases']:
        if pop == 'readings':
            for e in refs:
                if e != lab:
                    cat[f'{lab}>{e}'][cause] += 1
    direction = collections.Counter()
    for k, n in causes['readings'].items():
        direction[DIRECTION[k]] += n
    full = json.loads(OUT_JSON.read_text(encoding='utf-8'))['main']
    from full_run_28092026.worked_example import NAME
    head = full['headline']['fac']
    bound = {}
    for name, v in (('combined', 'C2 C1 with U2 and S2'), ('half_unit', 'D1 half-unit last-digit windows')):
        big = max(full[v]['delta'].values(), key=abs)
        fac = full[v]['fac']
        # the pairs of models whose order the variant reverses, each with its gap as scored
        swaps = sorted([NAME[a], NAME[b], round(head[a] - head[b], 4)] for a in head for b in head
                       if head[a] > head[b] and fac[a] < fac[b])
        bound[name] = {'largest_change': big, 'largest_rounded': round(abs(big), 3), 'at_most': ceil3(big),
                       'tau': full[v]['tau'], 'moved': full[v]['moved'], 'swaps': swaps}
    out = {'by_kind': {k: {'study': sa['by_kind'][k], 'readings': ra['by_kind'][k]} for k in KINDS},
           'all': {'study': sa['all'], 'readings': ra['all']}, 'readings_direction': dict(direction),
           'readings_by_category': {c: dict(v) for c, v in cat.items()}, 'study_by_cause': dict(causes['study']),
           'label_part_templates': len(LABEL_PARTS), 'bound': bound, 'cases': {}}
    pops = {(p['model'], p['item_id']): (p, lab) for p, lab in zip(readings, r_labels)}
    # the lenient case: an answer one unit off in its last digit, which both experts call incorrect
    p, lab = pops[CASES['lenient']]
    vals = {m['id']: m['value'] for m in CV.milestone_cache(score.pool_items())[p['item_id']]}
    want = A.targets(p['item'], tuple(vals.values()))
    gold = want['numbers'][0][0]
    exact = vals['L'] / vals['lam_e'] * 60.0
    stated = [v for v, u in A.values(A.segment(p['text'])) if CV._match(gold, [(v, u)], A.REL, False, want['numbers'][0][1])]
    assert lab == 'correct' and p['readings'] == ['incorrect', 'incorrect'] and len(stated) == 1, (lab, p['readings'], stated)
    assert round(exact, 2) == gold != stated[0], (exact, gold, stated)
    out['cases']['lenient'] = {'item': p['item_id'], 'branch': BRANCH[p['item']['branch']], 'level': p['item']['level'],
                               'stated': stated[0], 'gold': gold, 'unrounded': round(exact, 4), 'check': lab,
                               'experts': p['readings']}
    # the strict case: an equivalent signal the number rule scores partial, which both experts call correct
    p, lab = pops[CASES['strict']]
    g = re.search(r'y\[n\] = cos\(\((\d+)\*pi\)/(\d+)\*n\)', A.segment(p['item']['solution']))
    r = re.search(r'y\[n\] = \\cos\\left\(\\frac\{\\pi\}\{(\d+)\}n\\right\)', A.segment(p['text']))
    m_ = re.search(r'downsampled by a factor of M = (\d+)', p['item']['question'])
    stored = next(x for x in CV.load_rows('main', p['model']) if x['item_id'] == p['item_id'])['answer']
    assert g and r and m_ and lab == 'partial' and p['readings'] == ['correct', 'correct'], (g, r, lab, p['readings'])
    num, den, rden = int(g.group(1)), int(g.group(2)), int(r.group(1))
    assert den == rden and (num + 1) / den == 2 and (stored['matched'], stored['of']) == (1, 3), (num, den, rden, stored)
    out['cases']['strict'] = {'item': p['item_id'], 'branch': BRANCH[p['item']['branch']], 'level': p['item']['level'],
                              'M': int(m_.group(1)), 'gold_num': num, 'den': den, 'check': lab, 'matched': stored['matched'],
                              'of': stored['of'], 'experts': p['readings']}
    return out


def block_tex(nums: dict) -> str:
    rows = []
    for k in KINDS:
        s, r = nums['by_kind'][k]['study'], nums['by_kind'][k]['readings']
        rows.append(f'        {KIND_NAME[k]} & {s[1]} & {s[0] / s[1]:.3f} & {s[2] / s[3]:.3f} & {r[0]} of {r[1]} \\\\')
    s, r = nums['all']['study'], nums['all']['readings']
    a, b = nums['cases']['lenient'], nums['cases']['strict']
    lvl = {'Easy': 'Easy', 'Intermediate': 'Intermediate', 'Advanced': 'Advanced'}
    return '\n'.join([
        '% BEGIN GENERATED tab:answer_kinds and tab:answer_cases (full_run_28092026/answer_audit.py --write)',
        '\\begin{table}[t]',
        '    \\centering',
        '    \\small',
        '    \\setlength{\\tabcolsep}{4pt}',
        '    \\begin{tabular}{@{}l r r r r@{}}',
        '        \\toprule',
        '        \\rowcolor{gray!10}',
        '        \\textbf{Answer kind} & \\textbf{Study} & \\textbf{Three-way} & \\textbf{Not partial} & \\textbf{Readings} \\\\',
        '        \\midrule', *rows, '        \\midrule',
        f'        All & {s[1]} & {s[0] / s[1]:.3f} & {s[2] / s[3]:.3f} & {r[0]} of {r[1]} \\\\',
        '        \\bottomrule',
        '    \\end{tabular}',
        '    \\caption{\\textbf{Final-answer check by answer kind.} Study: the expert study\'s responses of each kind and',
        '    the check\'s agreement with the experts\' verdict, three-way and on the responses the experts did not call',
        '    partial. Readings: the domain experts\' readings of the evaluated models\' verdicts that agree with the',
        '    check; the sample is stratified by verdict, so these are counts, not rates for a kind.}',
        '    \\label{tab:answer_kinds}',
        '\\end{table}',
        '',
        '\\begin{table*}[t]',
        '    \\centering',
        '    \\small',
        '    \\renewcommand{\\arraystretch}{1.15}',
        '    \\begin{tabular}{@{}p{0.2\\linewidth} p{0.12\\linewidth} p{0.12\\linewidth} p{0.07\\linewidth} '
        'p{0.09\\linewidth} p{0.26\\linewidth}@{}}',
        '        \\toprule',
        '        \\rowcolor{gray!10}',
        '        \\textbf{Question} & \\textbf{Final answer} & \\textbf{Gold} & \\textbf{Check} & \\textbf{Experts} & '
        '\\textbf{Why they differ} \\\\',
        '        \\midrule',
        f"        Average time in the system of an M/M/1/K queue ({a['branch']}, {lvl[a['level']]}) & "
        f"{a['stated']:.2f} minutes & {a['gold']:.2f} minutes & {a['check']} & {a['experts'][0]} (both) & "
        f"the rule accepts a value one unit off in its last digit; the unrounded value, "
        f"$L/\\lambda_{{\\mathrm{{eff}}}}$, is {a['unrounded']:.4f} minutes, which rounds to {a['gold']:.2f}, "
        f"not {a['stated']:.2f} \\\\",
        f"        Output signal after downsampling by $M = {b['M']}$ ({b['branch']}, {lvl[b['level']]}) & "
        f"$y[n] = \\cos(\\pi n/{b['den']})$ & $y[n] = \\cos({b['gold_num']}\\pi n/{b['den']})$ & {b['check']} & "
        f"{b['experts'][0]} (both) & the same sequence, since ${b['gold_num']}\\pi/{b['den']} = 2\\pi - \\pi/{b['den']}$; "
        f"the check reads the stated numbers and finds {b['matched']} of its {b['of']} targets \\\\",
        '        \\bottomrule',
        '    \\end{tabular}',
        '    \\caption{\\textbf{One discrepancy in each direction.} Two final answers of the evaluated models on which',
        '    the check\'s verdict differs from both domain experts\' readings: it credits the first more and the second',
        '    less than they do.}',
        '    \\label{tab:answer_cases}',
        '\\end{table*}',
        '% END GENERATED tab:answer_kinds and tab:answer_cases', ''])


def write_paper(check_only: bool) -> int:
    """Compute the paper's numbers from the two expert populations (no pass over the headline responses), write the block and the JSON;
    with check_only, compare the block on disk with what would be written."""
    items = score.pool_items()
    ms = CV.milestone_cache(items)
    study, readings = study_population(), readings_population(items, ms)
    s_labels = [verdict(p['text'], p['item'], p['ms'], {}) for p in study]
    r_labels = [verdict(p['text'], p['item'], p['ms'], {}) for p in readings]
    causes = causes_table(study, s_labels, readings, r_labels)
    if not causes.get('available'):
        raise SystemExit('the causes file is not here: the paper numbers need it')
    nums = paper_numbers(study, s_labels, readings, r_labels, causes)
    tex = block_tex(nums)
    if check_only:
        same = BLOCK.exists() and BLOCK.read_text(encoding='utf-8') == tex
        print(f"{BLOCK.name}: {'current' if same else 'DIFFERS from what --write would write'}")
        return 0 if same else 1
    BLOCK.write_text(tex, encoding='utf-8', newline='\n')
    PAPER_JSON.write_text(json.dumps(nums, indent=1), encoding='utf-8')
    print(f'wrote {BLOCK} and {PAPER_JSON}')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--write', action='store_true', help='write the paper block and results/answer_audit_paper.json')
    ap.add_argument('--check-block', action='store_true', help='exit 1 unless the paper block is what --write writes')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.write or a.check_block:
        return write_paper(a.check_block)
    res = run(a.workers)
    L = report(res)
    REPORT.write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    CASES_MD.write_text('\n'.join(cases_lines(res)) + '\n', encoding='utf-8', newline='\n')
    OUT_JSON.write_text(json.dumps({
        'generated_by': 'python -m full_run_28092026.answer_audit',
        'written_at_utc': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        'agreement': {n: {v: res['pops'][n][v]['agreement'] for v in VARIANTS} for n in ('study', 'readings')},
        'causes': {k: dict(res['causes'].get(k, {})) for k in ('study', 'readings')},
        'main': res['main'], 'paths': res['paths'], 'coverage': res['coverage']}, indent=1), encoding='utf-8')
    CHANGES_JSON.write_text(json.dumps(res['changes'], indent=0), encoding='utf-8')
    print('\n'.join(L))
    return 0


# ------------------------------------------------------------------ selftest

def selftest() -> int:
    bad = 0

    def check(name, got, want):
        nonlocal bad
        if got != want:
            bad += 1
            print(f'FAIL {name}: got {got!r}, want {want!r}')

    check('spaced minus after a token', signed('(9.311e-05 x_hat - 9.890e-05 y_hat)'), '(9.311e-05 x_hat -9.890e-05 y_hat)')
    check('spaced minus after a brace', signed('\\hat{x} - 1.17'), '\\hat{x} -1.17')
    check('arithmetic minus kept apart', signed('5 - 3 = 2'), '5 - 3 = 2')
    check('unicode minus and en dash', signed('−4.05 and –6534.8'), '-4.05 and -6534.8')
    check('imaginary part', [v for v, _ in A.values(signed('57.63 - j\\,85.44'))], [57.63, -85.44])
    check('imaginary part after a digit', [v for v, _ in A.values(signed('= 44.26-j44.88 V'))], [44.26, -44.88])
    check('LaTeX spacing around a minus', [v for v, _ in A.values(signed('-0.016\\;\\cos(26.58\\,t)\\;-\\;0.082\\;\\sin'))
                                           if abs(v) < 1], [-0.016, -0.082])
    check('grouped digits untouched', [v for v, _ in A.values(signed('AOQ = 28\\,570 ppm'))], [28570.0])
    v = values_v('F = (-1.26969 x + 1.16486 y)\\times 10^{-4} N', {'factored': True})
    check('factored power of ten', any(abs(x + 1.26969e-4) < 1e-12 for x, _ in v), True)
    check('no units is not a stance', A._stance(NO_UNITS.sub(' ', 'a) not memoryless; b) causal. no units.'), 'causal'), 'yes')
    check('label: criteria after the verdict', word_hit('flow regime: transitional, i.e. between laminar and turbulent',
                                                        'transitional', 'revised'), True)
    check('label: no underdamped', word_hit('the system is critically damped, so no underdamped omega_d', 'critically',
                                            'revised'), True)
    check('label: incrementally linear', word_hit('the system is nonlinear. it is incrementally linear', 'nonlinear',
                                                  'revised'), True)
    check('label: yes to linear', word_hit('answer: yes dimensionless', 'linear', 'revised'), True)
    check('unit scale', (unit_scale('mJ'), unit_scale('minutes'), unit_scale('kN')), (('J', 1e-3), ('time', 60.0), ('N', 1e3)))
    check('unit slip refused', unit_explained(3.33, 'is 3.33 minutes', [(200.0, 1 / 60)], 'answer is: **200** minutes'), False)
    check('unit conversion allowed', unit_explained(112.12, 'is 112.12 mJ.', [(0.11212, 1e3)], 'E_b = 0.11212 J'), True)
    lim = {'template_id': 'template_limiting_reactant', 'answer_type': 'multipart', 'question': '',
           'solution': '**Answer:** The limiting reactant is **O2(g)**. The final number of moles are: ...'}
    reader = label_targets(lim, [])[0][2]
    check('limiting reactant in LaTeX', reader('the limiting reactant is $\\text{o}_2$. final moles: ...'), True)
    check('limiting reactant before "is"', reader('**limiting reactant:** $\\text{o}_2$'), True)
    check('limiting reactant missing', reader('c3h8: 13.74 mol o2: 303.7 mol'), False)
    check('limiting reactant wrong', reader('the limiting reactant is ch4.'), False)
    ali = {'template_id': 'template_decimation_aliasing_analysis', 'answer_type': 'multipart', 'question': '',
           'solution': '**Answer:** a) y[n] = cos(pi/5*n). b) Aliasing **did not** occur. c) omega_a = pi/5.'}
    reader = label_targets(ali, [])[0][2]
    check('no aliasing occurred', reader('b) **answer:** no aliasing occurred c) 3pi/8'), True)
    check('aliasing occurred, gold no', reader('b) aliasing occurred'), False)
    item = {'item_id': 'x#0', 'template_id': 'template_x', 'question': 'Given 3.0 m.', 'answer_type': 'scalar',
            'solution': '**Answer:** The deflection is 7.3 mm'}
    check('rule: default reproduces the sign clause', verdict('**Answer:** -7.326 mm', item, (7.3,)), 'correct')
    check('rule: sign clause off', verdict('**Answer:** -7.326 mm', item, (7.3,), {'abs': 'off'}), 'incorrect')
    check('rule: sign policy accepts a sign in words', verdict('**Answer:** -7.33 mm (downward)', dict(
        item, solution='**Answer:** The deflection is -7.3 mm'), (-7.3,), {'abs': 'policy'}), 'correct')
    print(f'{"selftest passed" if not bad else f"{bad} selftest failures"}')
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())
