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
"""
import functools
import os
import re
import sys

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
import milestones  # noqa: E402

TOL = 0.02            # kept for callers that pass a relative tolerance explicitly
REL = 0.002           # the window [0.0017, 0.0029] all agrees with the experts; see X4
SCALES = milestones.SCALES + (1e12, 1e-12, 1e-2)   # pF/m against F/m, and percent
NUM = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')
SUB = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
HEADING = re.compile(r'(?i)#+\s*final\s+answer')
ANSWER = re.compile(r'(?i)(?:#+\s*final\s+answer|\*{0,2}answer\s*\(?[a-z]?\)?\s*\*{0,2}\s*[:\-])')
WINDOW = 700          # how far an answer segment runs: a conclusion, not a second solution
# Questions that prescribe their own rounding ("to 4 decimals, round half up") are scored
# on the digits: the experts hold a trace to the precision the question asked for.
ROUNDING = re.compile(r'(?i)round half up|to \d+ decimals?|nearest whole|to \d+ significant')
BOLD = re.compile(r'\*\*([^*]+)\*\*')
WORD = re.compile(r'[A-Za-z][A-Za-z\-]{2,}')
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


def values(seg):
    """Every number an answer states, with the precision it was displayed at.

    Returns (value, ulp) pairs, where ulp is one unit in the last digit shown: `11.8` is
    (11.8, 0.1) and `114` is (114.0, 1.0). That precision is what makes a rounded answer
    judgeable - see match().

    Traces write the same value many ways - `5.52 x 10^-5` in unicode superscripts, LaTeX
    `7.34 \times 10^{-11}`, `1.46 x 10^(-10)`, `\frac{5}{2}`, subscripted `P_a(p_1)` - and a
    checker reading only plain decimals marks correct answers wrong. Unit exponents go
    first, or `m^3/s` contributes a 3 and `m/s^2` a 2 that then masquerade as answers.
    """
    t = seg
    t = re.sub(r'(\d)\s*[x*\u00d7]\s*10\s*\^\s*\(?\{?\s*([-+]?\d+)\s*\}?\)?', r'\1e\2', t)
    t = re.sub(r'(\d)\s*[x*\u00d7]\s*10\s*([\u207b\u207a]?[\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+)',
               lambda m: m.group(1) + 'e' + m.group(2).translate(SUP), t)
    t = t.replace('\\times', ' x ').replace('\\cdot', ' * ').replace('\\,', ' ')
    t = re.sub(r'\\d?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r' \1/\2 ', t)
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
    return out


def _ulp(lit):
    """One unit in the last digit shown: '11.8' -> 0.1, '114' -> 1.0, '2.05e7' -> 1e5."""
    m = re.match(r'[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?$', lit)
    if not m:
        return 0.0
    return 10.0 ** ((int(m.group(3)) if m.group(3) else 0) - len(m.group(2)))


def match(gold, have, rel=REL, exact=False, gold_ulp=0.0):
    """Does any stated value equal `gold`?

    Three ways to be right, because one relative tolerance cannot serve them all:

      within `rel` of the gold                 `2.05e7` for 20,464,069 (0.18%)
      within one unit of ITS last digit        `114` for 113.55 - a correct rounding at the
                                               precision the trace chose to display
      within one unit of the GOLD's last digit `7.326` for a gold written `7.3`

    while `50.20` for 50.05 fails all three: it claims two decimals and gets them wrong.

    `exact` is for items whose QUESTION prescribes the rounding ("to 4 decimals, round half
    up"): there the experts require the digits, and a near miss is a miss.
    """
    for v, u in have:
        for sc in SCALES:
            t, tu = v * sc, u * sc
            if exact:
                if abs(t - gold) <= 1e-9 * max(1.0, abs(gold)):
                    return True
            elif abs(t - gold) <= max(rel * abs(gold), tu, gold_ulp):
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
    out = []
    for m in BOLD.finditer(seg):
        for w in WORD.findall(m.group(1)):
            if w.lower() in VERDICT_WORDS:
                out.append(w.lower())
    for w in WORD.findall(seg):
        if w.lower() in VERDICT_WORDS and w.lower() not in out:
            out.append(w.lower())
    return out


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
    if typ in ('scalar', 'array') and keep:
        keep = keep[-1:]                       # one value is asked for: the last stated
    return {'numbers': keep, 'words': _words(seg, typ), 'exact': bool(ROUNDING.search(item.get('question') or ''))}


def verdict(text, item, tol=None, ms=None):
    """('correct' | 'partial' | 'incorrect', detail). Every part must match for correct.

    `tol` overrides the relative tolerance for callers that want to sweep it; the ulp rule
    in match() applies either way.
    """
    want = targets(item, ms)
    seg = segment(text)
    have = values(seg)
    absolute = [(abs(v), u) for v, u in have]
    low = seg.lower()
    rel = REL if tol is None else tol
    hits = []
    for v, gu in want['numbers']:
        # |v| too: a deflection the gold states as 7.3 mm downward and the trace as
        # -7.326 mm is the same answer under a different sign convention.
        hits.append(match(v, have, rel, want['exact'], gu)
                    or match(abs(v), absolute, rel, want['exact'], gu))
    for w in want['words']:
        # the LAST mention decides: an answer block that restates the criteria ("turbulent
        # occurs for Re > 5e6 ... the flow is laminar") must not be credited for both.
        last = None
        for m in re.finditer(r'\b(%s)\b' % '|'.join(sorted(VERDICT_WORDS)), low):
            last = m.group(1)
        hits.append(last == w)
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
    ]
    bad = 0
    for sol, q, typ, ms, trace, want in cases:
        got, _d = verdict(trace, item(sol, q, typ), ms=ms)
        if got != want:
            bad += 1
            print('FAIL want %-9s got %-9s | %s' % (want, got, ' '.join(trace.split())[:60]))
    print('%d/%d cases pass' % (len(cases) - bad, len(cases)))
    return bad


if __name__ == '__main__':
    raise SystemExit(selftest())
