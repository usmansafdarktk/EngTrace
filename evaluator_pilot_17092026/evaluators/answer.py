"""Does a trace's final answer match the gold's? correct / partial / incorrect.

E0's check is wrong on 72 of the pilot's 300 traces, 68 of them traces the experts call
correct (RESULTS_X1 Finding 1, E0_RERUN.md). Three things cause that, and this module is
built around fixing each one:

  RESTATED INPUTS (E0-F1). The gold's answer line often repeats the question's values -
  "the molar volume of n-Pentane at 283.81 K is 113.55 cm3/mol". Any number on that line
  that also appears in the QUESTION is an input, not an answer, and is dropped. That rule
  needs no per-template knowledge, which is why it is used instead of one.

  QUALITATIVE ANSWERS (E0-F2). Several templates answer with a word - turbulent,
  overdamped - and some state no number at all. Words the gold emphasises, and words from
  a fixed engineering list, become targets in their own right.

  MULTI-PART ANSWERS. `answer_type` is multipart, vector, classification or symbolic for
  9 of the 15 templates: two velocity components, a value and a regime, a symbolic form
  and its value. A trace that gets one part right and the other wrong is neither correct
  nor incorrect, and the experts labelled 19 such traces `partial`. Every target must
  match for `correct`, some for `partial`, none for `incorrect`.

A number matches under any of the unit factors in milestones.SCALES (m vs mm, fraction vs
percent), so a trace answering in m3/mol against a gold in cm3/mol is not marked wrong for
the unit alone. TOL is relative: the experts accept a rounded answer (56.27 for 56.29) and
this is the final answer, not a step, so it is not held to display precision - that rule
belongs to analysis/digit_rule.py, which is about intermediate arithmetic.

  PI-FRACTIONS AND RENDERINGS (D-169, after a read of the full run's verdicts). A pi-fraction is one
  value: `(841*pi)/2447`, `\\dfrac{841\\pi}{2447}`, `841\\pi/2447`, `(3/17)\\pi` and `0.4\\pi` all state
  a*pi/b. Until D-169 the parenthesis stopped the rule at 841*pi, a value the template also
  computes, so the gold's own target was wrong on two signals templates, and LaTeX's `\\pi` never
  joined its coefficient. Where a gold line's pi-fraction has a numerator and a denominator that are
  both computed (the reduced fraction's `num` and `den`), the fraction's value is the target, so a
  trace stating it as a decimal or unreduced is right. And a scalar gold line that states one
  quantity in two units ("0.088960 radians, or 5.0970 degrees") accepts either, the second matched
  at unit scale. The forms occur in none of the pilot's 15 templates, so the experts never
  arbitrated them; ANSWER_FORM_AUDIT.md and ANSWER_FORM_FIX.md in full_run_28092026/ measure them.
"""
import functools
import os
import math
import re
import sys

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
import milestones  # noqa: E402

TOL = 0.02            # kept for callers that pass a relative tolerance explicitly
REL = 0.002           # the window [0.0017, 0.0029] all agrees with the experts; see X4
SCALES = milestones.SCALES + (1e12, 1e-12, 1e-2)   # pF/m against F/m, and percent
# A digit written as a subscript - `p_1`, `x_{1}` - names a thing, it is not a value (D-138);
# E3's reader already skips it.
NUM = re.compile(r'(?<![\d.])(?<!_)(?<!_\{)[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')
SUB = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
HEADING = re.compile(r'(?i)#+\s*final\s+answer')
ANSWER = re.compile(r'(?i)(?:#+\s*final\s+answer|\*{0,2}answer\s*\(?[a-z]?\)?\s*\*{0,2}\s*[:\-])')
WINDOW = 700          # how far an answer segment runs: a conclusion, not a second solution
# Questions that prescribe their own rounding ("to 4 decimals, round half up") are scored
# on the digits: the experts hold a trace to the precision the question asked for.
ROUNDING = re.compile(r'(?i)round half up|to \d+ decimals?|nearest whole|to \d+ significant')
BOLD = re.compile(r'\*\*([^*]+)\*\*')
WORD = re.compile(r'[A-Za-z][A-Za-z\-]{2,}')
# Full-pool gold validation (D-120). The pilot's 15 templates never answered "linear",
# a labelled yes/no, or a bare pi; three pool templates do, and their gold scored
# incorrect against itself. Each rule below fires only on those forms, so the pilot's
# scoring is unchanged.
LINEAR = ('linear', 'nonlinear')                  # a family of its own: see verdict()
LABELED = re.compile(r'(?i)\b([a-z][a-z\-]{2,})\s*:\s*\*{0,2}\s*(yes|no)\b')
# A pi-fraction is one value (D-169): the coefficient may sit in parentheses with pi, and the
# denominator may follow the closing one. `pi*n` and `pi(` inside an expression are still not values.
PI_EXPR = re.compile(r'(?<![\w.])\(?\s*(\d+(?:\.\d+)?)?\s*\*?\s*(?:pi|\u03c0)\b(?!\s*[*(])'
                     r'\s*\)?(?:\s*/\s*\(?\s*(\d+(?:\.\d+)?)\s*\)?)?')
PI_OF_FRACTION = re.compile(r'(?<![\w.])\(\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)\s*\)\s*\*?\s*(?=(?:pi|\u03c0)\b)')
PI_FRACTION = re.compile(r'(?<![\w.])\(?\s*(\d+(?:\.\d+)?)\s*\*?\s*(?:pi|\u03c0)\s*\)?\s*/\s*\(?\s*(\d+(?:\.\d+)?)\s*\)?')
# Two renderings of one quantity on a scalar gold line - "0.088960 radians, or 5.0970 degrees" - are
# related by one of these factors, either way round (D-169). Powers of ten are not among them: the
# unit factors in match() already serve a single target, and two different quantities can differ by one.
RENDERINGS = (180.0 / math.pi, 2.0 * math.pi)


def _pi_text(seg):
    """The text PI_EXPR reads: LaTeX's `\\pi`, `\\frac`, `\\left`, `\\right` and `$` resolved, and
    `(a/b)*pi` written `a*pi/b`, so that every pi-fraction is in the one shape the rule reads."""
    p = seg.replace('\\pi', 'pi').replace('\\left', ' ').replace('\\right', ' ').replace('$', ' ')
    p = re.sub(r'\\[dt]?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r' \1/\2 ', p)
    return PI_OF_FRACTION.sub(r' \1*pi/\2 ', p)


def _rendering(a, b):
    """Are a and b one quantity in two units: related by 180/pi or 2*pi, either way round?"""
    if not a or not b:
        return False
    r = abs(a / b)
    return any(abs(r - f) <= 1e-3 * f or abs(r - 1 / f) <= 1e-3 / f for f in RENDERINGS)


def _norm_linear(s):
    """'not linear' and 'non-linear' are one verdict, written 'nonlinear'."""
    return re.sub(r'(?i)\b(?:not\s+linear|non-\s?linear)\b', 'nonlinear', s)
# Words that are an answer in themselves. Units and prose are not.
VERDICT_WORDS = {
    'turbulent', 'laminar', 'transitional', 'overdamped', 'underdamped',
    'critically', 'undamped', 'safe', 'unsafe', 'yes', 'no', 'feasible', 'infeasible',
    'stable', 'unstable', 'acceptable', 'adequate', 'inadequate', 'subsonic', 'supersonic',
    'saturated', 'superheated', 'subcritical', 'supercritical', 'operational', 'failed',
}
STOP = {'the', 'is', 'are', 'and', 'of', 'at', 'in', 'for', 'per', 'with', 'approximately',
        'about', 'value', 'values', 'answer', 'final', 'total', 'average', 'expression',
        'numerical', 'magnitude', 'components', 'component', 'option', 'better', 'system',
        'flow', 'regime', 'damping', 'ratio', 'fraction', 'state', 'days', 'long-run'}


SUP = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
FRACTION = re.compile(r'(?<![\w.^])(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)(?![\w.])')
SEPARATOR, GROUPED = milestones.SEPARATOR, milestones.GROUPED     # 28\,570 is one number (D-137)


def values(seg):
    """Every number an answer states, with the precision it was displayed at.

    Returns (value, ulp) pairs, where ulp is one unit in the last digit shown: `11.8` is
    (11.8, 0.1) and `114` is (114.0, 1.0). That precision is what makes a rounded answer
    judgeable - see match().

    Traces write the same value many ways - `5.52 x 10^-5` in unicode superscripts, LaTeX
    `7.34 \times 10^{-11}`, `1.46 x 10^(-10)`, `\frac{5}{2}`, subscripted `P_a(p_1)` - and a
    checker reading only plain decimals marks correct answers wrong. Unit exponents go
    first, or `m^3/s` contributes a 3 and `m/s^2` a 2 that then masquerade as answers.

    `\\times` and `\\cdot` become `x` and `*` BEFORE the exponent rules read them (D-137):
    replaced after, as they were until the full run, `3.04 \\times 10^5` read as 3.04, 10 and 5.
    """
    t = GROUPED.sub(lambda m: SEPARATOR.sub('', m.group(0)), seg)
    t = t.replace('\\times', ' x ').replace('\\cdot', ' * ')
    t = re.sub(r'(\d)\s*[x*\u00d7]\s*10\s*\^\s*\(?\{?\s*([-+]?\d+)\s*\}?\)?', r'\1e\2', t)
    t = re.sub(r'(\d)\s*[x*\u00d7]\s*10\s*([\u207b\u207a]?[\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+)',
               lambda m: m.group(1) + 'e' + m.group(2).translate(SUP), t)
    t = t.replace('\\,', ' ')
    t = re.sub(r'\\[dt]?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r' \1/\2 ', t)
    t = re.sub(r'(?<=[A-Za-z])\s*\^\s*-?\d+', ' ', t)                    # m^3/s, m/s^2
    t = re.sub(r'(?<=[A-Za-z])[\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+', ' ', t)   # unicode units
    t = t.translate(SUB).translate(SUP)
    t = re.sub(r'\\left|\\right|\\[()\[\]]|\$|\\mathrm|\\text\{[^{}]*\}', ' ', t)
    if re.search(r'\d,\d{3}\b', t):
        t = t.replace(',', '')
    out = []
    for m in NUM.finditer(t):
        try:
            v = float(m.group(0))
        except ValueError:
            continue
        if v == v and abs(v) != float('inf'):
            out.append((v, _ulp(m.group(0))))
    for m in FRACTION.finditer(t):
        b = float(m.group(2))
        if b:
            out.append((float(m.group(1)) / b, 0.0))
    # pi written as a symbol - `omega_a = pi`, `0.4*pi`, `3pi/4` - is a value (D-120), and so is a
    # pi-fraction in any of its shapes, `(841*pi)/2447` or `\dfrac{841\pi}{2447}` (D-169);
    # `pi*n` inside an expression is not.
    for m in PI_EXPR.finditer(_pi_text(seg)):
        a = float(m.group(1)) if m.group(1) else 1.0
        b = float(m.group(2)) if m.group(2) else 1.0
        if b:
            out.append((a * math.pi / b, 0.0))
    return out


def _ulp(lit):
    """One unit in the last digit shown: '11.8' -> 0.1, '114' -> 1.0, '2.05e7' -> 1e5."""
    m = re.match(r'[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?$', lit)
    if not m:
        return 0.0
    return 10.0 ** ((int(m.group(3)) if m.group(3) else 0) - len(m.group(2)))


SLACK = 1e-9          # the windows below are inclusive; this keeps binary floating point from deciding a boundary


def match(gold, have, rel=REL, exact=False, gold_ulp=0.0, unit=1.0, scales=SCALES):
    """Does any stated value equal `gold`?

    Three ways to be right, because one relative tolerance cannot serve them all:

      within `rel` of the gold                 `2.05e7` for 20,464,069 (0.18%)
      within one unit of ITS last digit        `114` for 113.55 - a correct rounding at the
                                               precision the trace chose to display
      within one unit of the GOLD's last digit `7.326` for a gold written `7.3`

    while `50.20` for 50.05 fails all three: it claims two decimals and gets them wrong.

    `exact` is for items whose QUESTION prescribes the rounding ("to 4 decimals, round half
    up"): there the experts require the digits, and a near miss is a miss.

    A number's own last digit vouches for it only when one unit of that digit is smaller than
    the number (D-138). A bare `0` or `1` fails that: one unit either side reaches zero, so under
    some unit factor it lies within one unit of any gold - `1` at 1e6 is 1,000,000 +- 1,000,000.

    The windows are inclusive, and a value exactly one unit off sits on the boundary: `0.638` for
    0.639 is 0.0010000000000000009 away in binary, so the comparison carries a relative slack of
    SLACK, and the verdict no longer depends on how the difference rounds (D-147). `unit` scales
    the two last-digit windows: 1.0 is the rule the experts validated; 0.5 is the stricter
    "correct rounding" reading, reported as a sensitivity only.

    `scales` are the unit factors tried. A rendering (D-169) is matched with (1.0,) alone, because
    it already is the other unit: with the factors, a cut-off trace holding `pi d^4/32` and no
    answer matched 0.5315 rad through 32 +- 1 at the 1/60 factor.
    """
    for v, u in have:
        for sc in scales:
            t, tu = v * sc, u * sc
            if exact:
                if abs(t - gold) <= 1e-9 * max(1.0, abs(gold)):
                    return True
            elif abs(t - gold) <= max(rel * abs(gold), unit * tu if abs(v) > u else 0.0, unit * gold_ulp) * (1 + SLACK):
                return True
    return False


def segment(text):
    """What a reader would call the answer.

    A multi-part answer writes one marker per part - "Answer (a): 20,016,378" then
    "Answer (b): Turbulent" - so taking everything after the LAST marker keeps only the
    last part and throws the value away. Where the trace heads its conclusion (## Final
    Answer), the segment starts at that heading and keeps every part under it.
    """
    heads = list(HEADING.finditer(text))
    if heads:
        return text[heads[-1].start():heads[-1].start() + WINDOW]
    hits = list(ANSWER.finditer(text))
    if hits:
        start = hits[0].start() if len(hits) > 1 and hits[-1].start() - hits[0].start() < WINDOW \
            else hits[-1].start()
        return text[start:start + WINDOW]
    hits = list(re.finditer(r'(?i)\banswer\b', text))
    return text[hits[-1].start():hits[-1].start() + WINDOW] if hits else text[-WINDOW:]


def _words(seg, answer_type=None):
    """Verdict words the gold states - only where the answer IS a word.

    `answer_type` says when that is: a classification item answers "turbulent" or
    "overdamped". Elsewhere the same words describe the substance rather than answer the
    question - "the molar volume of *saturated* liquid Chlorine" - and a trace that simply
    states the value is right without repeating them.
    """
    if answer_type not in ('classification',):
        return []
    seg = _norm_linear(seg)
    out = []
    for m in BOLD.finditer(seg):
        for w in WORD.findall(m.group(1)):
            if w.lower() in VERDICT_WORDS or w.lower() in LINEAR:
                out.append(w.lower())
    for w in WORD.findall(seg):
        if (w.lower() in VERDICT_WORDS or w.lower() in LINEAR) and w.lower() not in out:
            out.append(w.lower())
    return out


def _labeled(seg, answer_type=None):
    """Labelled yes/no parts of a classification answer: 'Memoryless: No', 'Causal: Yes'.

    Scored label by label (D-120). The last-verdict-word rule cannot score two of them:
    it would credit 'Memoryless: Yes, Causal: Yes' for 'Memoryless: No, Causal: Yes'.
    """
    if answer_type not in ('classification',):
        return []
    out = []
    for m in LABELED.finditer(seg):
        pair = (m.group(1).lower(), m.group(2).lower())
        if pair not in out:
            out.append(pair)
    return out


def _stance(low, label):
    """yes / no for a labelled property: 'causal: no', 'is not causal', 'non-causal'."""
    lab = re.escape(label)
    last = None
    for m in re.finditer(r'\b%s\b\W{0,6}(yes|no)\b' % lab, low):
        last = m.group(1)
    if last:
        return last
    if re.search(r"\b(?:not|non-?|isn't|never)\s*%s\b" % lab, low):
        return 'no'
    if label == 'memoryless' and re.search(r'\b(?:has|have|with)\s+memory\b|\bdynamic\b', low):
        return 'no'
    if re.search(r'\b%s\b' % lab, low):
        return 'yes'
    return None


@functools.lru_cache(maxsize=None)
def _derive(item_id, solution, question, template_id):
    """milestones.build re-executes the template, so each item is derived once."""
    try:
        return tuple(m['value'] for m in milestones.build(
            {'item_id': item_id, 'solution': solution, 'question': question,
             'template_id': template_id})['milestones'])
    except Exception:
        return ()


def _computed(item, ms):
    """The quantities the gold derives. Passed in where the caller has them (E3 holds a
    store), otherwise derived here; an item that no longer reproduces byte-identically
    gives none, and the fallback in targets() applies."""
    if ms is not None:
        return list(ms)
    return list(_derive(item['item_id'], item['solution'], item['question'],
                        item.get('template_id')))


def targets(item, ms=None):
    """What the trace has to reproduce: the gold's answer values, and its verdict words.

    A number is a target only if the gold COMPUTED it - if it is one of the item's
    milestones, the quantities the template derives. That one rule, with no per-template
    knowledge, removes both kinds of noise:

      restated inputs (E0-F1)   "...of n-Pentane at 283.81 K is 113.55 cm3/mol": 283.81 is
                                given, never computed, so it is not a target.
      expression furniture      "C' = (2 * pi * epsilon) / ln(b/a). The value is 45.3":
                                the 2 belongs to the formula, not to the answer.

    Where the answer HAS parts - multipart and vector items ask for several quantities -
    every computed quantity is a target, because the gold often states only one of them on
    its answer line while the question asks for all six. That is what makes `partial` real:
    a trace with five of six right is neither correct nor incorrect.
    """
    seg = segment(item['solution'])
    given = milestones.numbers(item['question'])
    ms = _computed(item, ms)
    typ = item.get('answer_type')
    keep, seen = [], []
    for v, u in values(seg):
        hit = milestones.scaled_match(v, ms, milestones.DISPLAY_TOL) if ms else None
        # A computed value stays a target even when the question states the same number -
        # a normal depth of 1.500 m in a 1.5 m channel is the answer, not a restated input.
        if hit is None and any(milestones.close(v, g, 1e-9) for g in given):
            continue
        if hit is None or any(milestones.close(v, x, 1e-9) for x, _ in seen):
            continue
        seen.append((v, u))
        keep.append((v, u))
    if not keep:
        # The item no longer reproduces, so there are no computed values to appeal to.
        # Fall back to the last value the gold's answer states.
        nums = [(v, u) for v, u in values(seg)
                if not any(milestones.close(v, g, 1e-9) for g in given)]
        keep = (nums or values(seg))[-1:]
    # D-169: a pi-fraction whose numerator and denominator are both computed - the reduced
    # fraction's `num` and `den` - is one quantity, and its value is the target.
    for m in PI_FRACTION.finditer(_pi_text(seg)):
        a, b = float(m.group(1)), float(m.group(2))
        ia = next((k for k, (v, _) in enumerate(keep) if milestones.close(v, a, 1e-9)), None)
        ib = next((k for k, (v, _) in enumerate(keep) if milestones.close(v, b, 1e-9)), None)
        if ia is not None and ib is not None and ia != ib and b:
            keep = [x for k, x in enumerate(keep) if k not in (ia, ib)] + [(a * math.pi / b, 0.0)]
    renderings = []
    if typ in ('scalar', 'array') and keep:
        if typ == 'scalar':                    # D-169: the same quantity in another unit, if stated
            renderings = [(v, u) for v, u in keep[:-1] if _rendering(v, keep[-1][0])]
        keep = keep[-1:]                       # one value is asked for: the last stated
    labeled = _labeled(seg, typ)
    words = _words(seg, typ)
    if labeled:                                # scored label by label instead (D-120)
        words = [w for w in words if w not in ('yes', 'no')]
    return {'numbers': keep, 'renderings': renderings, 'words': words, 'labeled': labeled,
            'exact': bool(ROUNDING.search(item.get('question') or ''))}


# Verdict words that answer the same question (D-139). Only a word of the target's own family can
# overrule it: "critically damped; omega_d = 0 (no oscillation occurs)" ends on `no`, which
# answers a yes/no question, not this one.
FAMILIES = [('overdamped', 'underdamped', 'critically', 'undamped'), ('laminar', 'turbulent', 'transitional'),
            ('safe', 'unsafe'), ('yes', 'no'), ('feasible', 'infeasible'), ('stable', 'unstable'),
            ('acceptable', 'adequate', 'inadequate'), ('subsonic', 'supersonic'), ('saturated', 'superheated'),
            ('subcritical', 'supercritical'), ('operational', 'failed')]
# A mention negated just before it is not a verdict: "... since the system is not underdamped" (D-139).
NEGATED = re.compile(r"(?:\bnot|\bnon-?|n't|\bnever|\bneither|\bnor)\s*(?:\w+\s+)?$")


def _word_hit(low, w):
    """Is verdict word `w` the answer? The LAST mention of its family decides, a negated one
    skipped: an answer block that restates the criteria ("turbulent occurs for Re > 5e6 ... the
    flow is laminar") must not be credited for both."""
    if w in LINEAR:
        vocab = LINEAR
    else:
        vocab = next((f for f in FAMILIES if w in f), tuple(sorted(VERDICT_WORDS)))
    last = None
    for m in re.finditer(r'\b(%s)\b' % '|'.join(vocab), low):
        if NEGATED.search(low[max(0, m.start() - 24):m.start()]):
            continue
        last = m.group(1)
    return last == w


def verdict(text, item, tol=None, ms=None, unit=1.0, whole=False):
    """('correct' | 'partial' | 'incorrect', detail). Every part must match for correct.

    `tol` overrides the relative tolerance for callers that want to sweep it; the ulp rule
    in match() applies either way. `unit` scales the last-digit windows (see match).

    `whole` is the body-lenient reading of D-147, a sensitivity, never the score: a numeric
    part the Answer line leaves out is credited when the trace states it anywhere, but only
    for a trace whose Answer line already matches at least one part. An Answer line that
    matches nothing stays incorrect whatever the body holds: crediting a wrong final answer
    because the right number appears in the working moved 27 pilot verdicts away from the
    experts (BOUNDARY_AUDIT.md), and it is not the omission the reading is about. The prompt
    asks for the final result on the Answer line, and the experts validated the check there.
    """
    want = targets(item, ms)
    seg = segment(text)
    have = values(seg)
    absolute = [(abs(v), u) for v, u in have]
    low = _norm_linear(seg).lower()
    rel = REL if tol is None else tol
    hits = []
    for v, gu in want['numbers']:
        # |v| too: a deflection the gold states as 7.3 mm downward and the trace as
        # -7.326 mm is the same answer under a different sign convention. A rendering of the
        # target in another unit counts as the target, at unit scale (D-169).
        hits.append(match(v, have, rel, want['exact'], gu, unit)
                    or match(abs(v), absolute, rel, want['exact'], gu, unit)
                    or any(match(w, have, rel, want['exact'], wu, unit, scales=(1.0,))
                           or match(abs(w), absolute, rel, want['exact'], wu, unit, scales=(1.0,))
                           for w, wu in want.get('renderings', [])))
    for w in want['words']:
        hits.append(_word_hit(low, w))
    for label, value in want.get('labeled', []):
        hits.append(_stance(low, label) == value)
    if whole and any(hits) and not all(hits[:len(want['numbers'])]):
        body = values(text)
        body_abs = [(abs(v), u) for v, u in body]
        for k, (v, gu) in enumerate(want['numbers']):
            hits[k] = hits[k] or match(v, body, rel, want['exact'], gu, unit) \
                or match(abs(v), body_abs, rel, want['exact'], gu, unit)
    if not hits:
        return 'incorrect', {'targets': want, 'matched': 0, 'of': 0}
    n = sum(hits)
    label = 'correct' if n == len(hits) else ('partial' if n else 'incorrect')
    return label, {'targets': want, 'matched': n, 'of': len(hits), 'answer_type': item.get('answer_type')}


def correct(text, item, tol=None, ms=None):
    """Binary form, for callers that score accuracy."""
    return verdict(text, item, tol, ms)[0] == 'correct'


def selftest():
    """Every defect this module was written for, pinned as a case.

        python evaluator_pilot_17092026/evaluators/answer.py

    Each case is a real pattern from the pilot's 300 labelled traces, with the verdict the
    experts gave. A failure here means a fix regressed.
    """
    item = lambda sol, q='', t='scalar', mid='x#0': (
        {'item_id': mid, 'template_id': 'template_x', 'solution': sol, 'question': q, 'answer_type': t})
    PI_GOLD = ('**Answer:** The resulting discrete-time signal is x[n] = 14.3 * cos((841*pi)/2447*n - 0.81 rad), '
               'and its discrete-time angular frequency is omega = (841*pi)/2447 rad/sample.')
    PI_GOLD_ND = ('**Answer:** The resulting discrete-time signal is x[n] = 9.55 * cos((129*pi)/1574*n - 135 deg), '
                  'and its discrete-time angular frequency is omega = (129*pi)/1574 rad/sample.')
    TWO_UNITS = '**Answer:** The total angle of twist at the free end C is 0.088960 radians, or 5.0970 degrees.'
    cases = [
        # E0-F1: the gold restates an input on its answer line; the trace gives the value
        ('**Answer:** The molar volume of saturated liquid n-Pentane at 283.81 K is **113.55 cm3/mol**.',
         'Find the molar volume at 283.81 K.', 'scalar', (113.55,), '**Answer:** 113.55 cm3/mol', 'correct'),
        # the same gold, a wrong value
        ('**Answer:** The molar volume of saturated liquid n-Pentane at 283.81 K is **113.55 cm3/mol**.',
         'Find the molar volume at 283.81 K.', 'scalar', (113.55,), '**Answer:** 98.20 cm3/mol', 'incorrect'),
        # "saturated" describes the substance, it is not a verdict word to be matched
        ('**Answer:** The molar volume of saturated liquid Chlorine is **56.29 cm3/mol**.',
         '', 'scalar', (56.2858,), '**Answer:** 56.27 cm3/mol', 'correct'),
        # multi-part: one marker per part, both must be right
        ('**Answer:** a) The Reynolds number is **20,464,069**. b) The flow regime is **turbulent**.',
         '', 'classification', (20464069.0,),
         '## Final Answer\n**Answer (a):** 2.05e7\n**Answer (b):** Turbulent', 'correct'),
        ('**Answer:** a) The Reynolds number is **20,464,069**. b) The flow regime is **turbulent**.',
         '', 'classification', (20464069.0,),
         '## Final Answer\n**Answer (a):** 2.05e7\n**Answer (b):** Laminar', 'partial'),
        # notation: superscripts, LaTeX, fractions
        ('**Answer:** The force is F = 5.520e-05 N.', '', 'scalar', (5.52e-05,),
         '**Answer:** 5.52 × 10⁻⁵ N', 'correct'),
        ('**Answer:** The velocity is u = 2.5x^2 + 4xy.', '', 'symbolic', (2.5, 4.0),
         r'**Answer:** \( \frac{5}{2}x^2 + 4xy \)', 'correct'),
        # sign convention: a downward deflection stated negative
        ('**Answer:** The deflection at the free end is 7.3 mm', '', 'scalar', (7.3,),
         '**Answer:** -7.326 [mm]', 'correct'),
        # the answer value also appears in the question - it is still the answer
        ('**Answer:** The normal depth is 1.500 m', 'A channel 1.5 m wide, depth 3.0 m.', 'array',
         (1.5,), '**Answer:** 1.50 m', 'correct'),
        # a unit factor apart is the same answer
        ('**Answer:** The volume is 113.55 cm3/mol', '', 'scalar', (113.55,),
         '**Answer:** 0.11355 dm3/mol', 'correct'),
        # D-137: forms from the full run's traces, the verdict the gold's arithmetic gives
        ('**Answer:** a) The Reynolds number is approximately **304,114**. b) The flow regime is **turbulent**.',
         '', 'classification', (304114.0,),
         '## Final Answer\n**Answer:** a) \\(Re \\approx 3.04 \\times 10^5\\) dimensionless; b) turbulent flow',
         'correct'),
        ('**Answer:** The force is F = 5.520e-05 N.', '', 'scalar', (5.52e-05,),
         '**Answer:** \\(F = 5.52 \\times 10^{-5}\\ \\text{N}\\)', 'correct'),
        ('**Answer:** The force is F = 5.520e-05 N.', '', 'scalar', (5.52e-05,),
         '**Answer:** \\(F = 5.52 \\cdot 10^{-5}\\) N', 'correct'),
        ('**Answer:** The AOQ is 28,570 ppm.', '', 'scalar', (28570.0,),
         '**Answer:** \\(28\\,570\\ \\text{ppm}\\)', 'correct'),
        ('**Answer:** The AOQ is 11,003 ppm.', '', 'scalar', (11003.0,),
         '**Answer:** \\(AOQ = 11{,}003\\) ppm', 'correct'),
        ('**Answer:** The volume is 3,021 m3.', '', 'scalar', (3021.0,),
         '**Answer:** 3 021 m³ of soil', 'correct'),
        ('**Answer:** The phase is 10/3 rad.', '', 'scalar', (10 / 3,),
         '**Answer:** \\(\\tfrac{10}{3}\\) rad', 'correct'),
        # and a wrong value in the same notation is still wrong
        ('**Answer:** a) The Reynolds number is approximately **304,114**. b) The flow regime is **turbulent**.',
         '', 'classification', (304114.0,),
         '## Final Answer\n**Answer:** a) \\(Re \\approx 3.40 \\times 10^5\\); b) turbulent', 'partial'),
        # D-138: a subscript or a bare 0 or 1 does not vouch for a value; a rounding still does
        ('**Answer:** The energy loss is 0.151 m.', '', 'scalar', (0.151,),
         '**Answer:** \\(y_1 \\approx 0.238\\ \\text{m}\\), \\(h_L \\approx 0.374\\ \\text{m}\\)', 'incorrect'),
        ('**Answer:** The argument is 15.15.', '', 'scalar', (15.15,),
         '**Answer:** \\(P_b = Q\\left(\\sqrt{2E_b/N_0}\\right)\\)', 'incorrect'),
        ('**Answer:** The volume is 3,364 m3.', '', 'scalar', (3364.0,),
         '**Answer:** Step 1 gives about 3.4 x 10^3 m3', 'correct'),
        ('**Answer:** The probability is 0.4987.', '', 'scalar', (0.4987,),
         '**Answer:** P = 0.5', 'correct'),
        # D-139: the verdict word's own family decides, a negated mention skipped
        ('**Answer:** The damping ratio is 1.0000. The system is **Critically Damped**.', '', 'classification',
         (1.0,), '## Final Answer\n**Answer:** ζ ≈ 1.00 (Critically Damped); ω_d = 0 rad/s (no oscillation occurs)',
         'correct'),
        ('**Answer:** The damping ratio is 1.0000. The system is **Critically Damped**.', '', 'classification',
         (1.0,), '## Final Answer\n**Answer:** ζ = 1.00, critically damped; ω_d = 0 (since the system is not underdamped)',
         'correct'),
        ('**Answer:** The damping ratio is 1.0000. The system is **Critically Damped**.', '', 'classification',
         (1.0,), '## Final Answer\n**Answer:** ζ = 1.00; the system is not critically damped but overdamped', 'partial'),
        # D-169: a pi-fraction is one value, on the gold's line (the target) and in the trace, in every shape
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = 1.0797 rad/sample', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** \\(\\omega = \\dfrac{841\\pi}{2447}\\ \\text{rad/sample}\\)', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** \\(\\omega = \\left(\\frac{841\\pi}{2447}\\right)\\) rad/sample', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = 841\\pi/2447 rad/sample', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = (841/2447)\\pi rad/sample', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = 0.3437π rad/sample', 'correct'),
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = 2642 rad/sample', 'incorrect'),      # 841*pi, the target before D-169
        # the known limit D-169 records, pinned so that a change to it is noticed: the coefficient of a wrong
        # `1.0797 pi` is a stated number within 0.2% of the target and is credited; not counting pi's
        # coefficients as values would remove this credit and six right ones (ANSWER_FORM_AUDIT.md, reading C)
        (PI_GOLD, '', 'symbolic', (), '**Answer:** omega = 1.0797\\pi rad/sample', 'correct'),
        # the numerator and denominator are computed: the fraction's value is the target, in any reduction
        (PI_GOLD_ND, '', 'symbolic', (1574.0, 129.0), '**Answer:** omega = 0.2575 rad/sample', 'correct'),
        (PI_GOLD_ND, '', 'symbolic', (1574.0, 129.0), '**Answer:** omega = 645π/7870 rad/sample', 'correct'),
        (PI_GOLD_ND, '', 'symbolic', (1574.0, 129.0), '**Answer:** omega = (129*pi)/1574 rad/sample', 'correct'),
        (PI_GOLD_ND, '', 'symbolic', (1574.0, 129.0), '**Answer:** omega = 0.2500 rad/sample', 'incorrect'),
        # pi inside an expression is still not a value, and a subscripted pi names a thing
        ('**Answer:** The value is 3.14.', '', 'scalar', (3.14,), '**Answer:** y[n] = cos(pi*n); the value is 2.0', 'incorrect'),
        ('**Answer:** The value is 3.14.', '', 'scalar', (3.14,), '**Answer:** \\(\\pi_1 = 2.0\\)', 'incorrect'),
        # D-169: one quantity in two units on a scalar gold line - either is the answer
        (TWO_UNITS, '', 'scalar', (0.08896, 5.097), '**Answer:** 0.08896 rad', 'correct'),
        (TWO_UNITS, '', 'scalar', (0.08896, 5.097), '**Answer:** 5.10 degrees', 'correct'),
        (TWO_UNITS, '', 'scalar', (0.08896, 5.097), '**Answer:** 0.0800 rad', 'incorrect'),
        # ...matched at unit scale: a cut-off trace holding pi d^4/32 and no answer is not credited
        ('**Answer:** The total angle of twist at the free end C is 0.53152 radians, or 30.4539 degrees.', '', 'scalar',
         (0.53152, 30.4539), '## Formulae\n* J = \\frac{\\pi d^4}{32}\n* phi = TL/GJ\n## Solution\n**Step', 'incorrect'),
        # two different quantities on one line are not renderings: the last stays the target, as before
        ('**Answer:** 1. Liquid-phase Volume (V_liq) ≈ **0.0907 L/mol** 2. Vapor-phase Volume (V_vap) ≈ **1.7205 L/mol**',
         '', 'array', (0.0907, 1.7205), '**Answer:** Liquid molar volume ≈ 0.0907 L/mol', 'incorrect'),
        ('**Answer:** 1. Liquid-phase Volume (V_liq) ≈ **0.0907 L/mol** 2. Vapor-phase Volume (V_vap) ≈ **1.7205 L/mol**',
         '', 'array', (0.0907, 1.7205), '**Answer:** liquid 0.0907 L/mol, vapor 1.7205 L/mol', 'correct'),
    ]
    bad = 0
    for sol, q, typ, ms, trace, want in cases:
        got, _d = verdict(trace, item(sol, q, typ), ms=ms)
        if got != want:
            bad += 1
            print('FAIL want %-9s got %-9s | %s' % (want, got, ' '.join(trace.split())[:60]))
    # D-147: the boundary is decided by the rule, not by binary rounding, and the half-unit and
    # whole-trace readings are what they say
    boundary = [
        # exactly one unit of the stated last digit off, where that window is the binding one (the
        # relative tolerance is narrower): inside the inclusive window, whatever the float noise -
        # 0.063 - 0.062 is 0.0010000000000000009 in binary
        ('**Answer:** x = 0.063', '', 'scalar', (0.063,), '**Answer:** x = 0.062', 'correct', 1.0, False),
        ('**Answer:** x = 0.063', '', 'scalar', (0.063,), '**Answer:** x = 0.062', 'incorrect', 0.5, False),
        ('**Answer:** V = 114 m3', '', 'scalar', (114.0,), '**Answer:** V = 113 m3', 'correct', 1.0, False),
        ('**Answer:** V = 114 m3', '', 'scalar', (114.0,), '**Answer:** V = 113 m3', 'incorrect', 0.5, False),
        # two units off is outside both windows
        ('**Answer:** V = 114 m3', '', 'scalar', (114.0,), '**Answer:** V = 112 m3', 'incorrect', 1.0, False),
        # a correct rounding passes under both readings
        ('**Answer:** V = 113.55 cm3/mol', '', 'scalar', (113.55,), '**Answer:** 114 cm3/mol', 'correct', 0.5, False),
        # the whole-trace reading credits a value stated in the body but not on the Answer line
        ('**Answer:** k = 74254 N/m and omega_n = 30.54 rad/s', '', 'multipart', (74254.0, 30.54),
         '**Step 2:** k_eq = 74254 N/m\n## Final Answer\n**Answer:** 30.54 rad/s', 'partial', 1.0, False),
        ('**Answer:** k = 74254 N/m and omega_n = 30.54 rad/s', '', 'multipart', (74254.0, 30.54),
         '**Step 2:** k_eq = 74254 N/m\n## Final Answer\n**Answer:** 30.54 rad/s', 'correct', 1.0, True),
        # ...but never a wrong Answer line whose body happens to hold the right number
        ('**Answer:** V = 113.55 cm3/mol', '', 'scalar', (113.55,),
         '**Step 3:** V = 113.55 cm3/mol\n## Final Answer\n**Answer:** 98.20 cm3/mol', 'incorrect', 1.0, True),
        ('**Answer:** k = 74254 N/m and omega_n = 30.54 rad/s', '', 'multipart', (74254.0, 30.54),
         '**Step 2:** k_eq = 74254 N/m, omega_n = 30.54\n## Final Answer\n**Answer:** 12.3 rad/s', 'incorrect', 1.0, True),
    ]
    for sol, q, typ, ms, trace, want, unit, whole in boundary:
        got, _d = verdict(trace, item(sol, q, typ), ms=ms, unit=unit, whole=whole)
        if got != want:
            bad += 1
            print('FAIL want %-9s got %-9s (unit %s, whole %s) | %s' % (want, got, unit, whole, ' '.join(trace.split())[:50]))
    cases += boundary
    # values() alone: grouping joins only a number's own digits (D-137)
    reads = [('\\(\\log_2\\,256 = 8\\)', 2256.0, False), ('\\(x_1\\,000\\)', 1000.0, False),
             ('\\boxed{1\\,335}\\ \\text{K}', 1335.0, True), ('\\(4\\,477.9\\ \\text{psi}\\)', 4477.9, True),
             ('\\(1\\,234\\,567\\)', 1234567.0, True), ('0.123\\,456', 123456.0, False),
             ('P_a(p_1) = 0.8002', 1.0, False), ('x_{12} = 3.5', 12.0, False), ('x_{12} = 3.5', 2.0, False)]
    for text, v, want in reads:
        if any(abs(x - v) < 1e-9 for x, _u in values(text)) != want:
            bad += 1
            print('FAIL values(%r) %s %s' % (text, 'misses' if want else 'reads', v))
    print('%d/%d cases pass' % (len(cases) + len(reads) - bad, len(cases) + len(reads)))
    return bad


if __name__ == '__main__':
    raise SystemExit(selftest())
