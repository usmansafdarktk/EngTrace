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
import os
import re
import sys

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
import milestones  # noqa: E402

TOL = 0.02
HEADING = re.compile(r'(?i)#+\s*final\s+answer')
ANSWER = re.compile(r'(?i)(?:#+\s*final\s+answer|\*{0,2}answer\s*\(?[a-z]?\)?\s*\*{0,2}\s*[:\-])')
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
    """Every number an answer states, however it is written.

    Traces write the same value many ways - `5.52 × 10⁻⁵`, `5.52 \\times 10^{-5}`,
    `\\frac{5}{2}`, `(5/2)x²` - and a checker that only reads plain decimals marks a
    correct answer wrong. Superscripts, LaTeX and simple fractions are folded in first.
    """
    t = re.sub(r'10\s*([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)', lambda m: '10^' + m.group(1).translate(SUP), seg)
    t = t.translate(SUP)
    t = t.replace('\\times', '*').replace('\\cdot', '*').replace('\\,', ' ')
    t = re.sub(r'\\d?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}', r' \1/\2 ', t)
    t = re.sub(r'\\left|\\right|\\[()\[\]]|\$|\\mathrm|\\text\{[^{}]*\}', ' ', t)
    nums = milestones.numbers(t)
    for m in FRACTION.finditer(t):
        b = float(m.group(2))
        if b:
            nums.append(float(m.group(1)) / b)
    return nums


def segment(text):
    """What a reader would call the answer.

    A multi-part answer writes one marker per part - "Answer (a): 20,016,378" then
    "Answer (b): Turbulent" - so taking everything after the LAST marker keeps only the
    last part and throws the value away. Where the trace heads its conclusion (## Final
    Answer), the segment starts at that heading and keeps every part under it.
    """
    heads = list(HEADING.finditer(text))
    if heads:
        return text[heads[-1].start():]
    hits = list(ANSWER.finditer(text))
    if hits:
        return text[hits[0].start():] if len(hits) > 1 and hits[-1].start() - hits[0].start() < 400 \
            else text[hits[-1].start():]
    hits = list(re.finditer(r'(?i)\banswer\b', text))
    return text[hits[-1].start():] if hits else text[-300:]


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


def _computed(item, ms):
    """The quantities the gold derives. Passed in where the caller has them (E3 holds a
    store), otherwise derived here; an item that no longer reproduces byte-identically
    gives none, and the fallback in targets() applies."""
    if ms is not None:
        return list(ms)
    try:
        return [m['value'] for m in milestones.build(item)['milestones']]
    except Exception:
        return []


def targets(item, ms=None):
    """What the trace has to reproduce: the gold's answer values, and its verdict words.

    A number on the gold's answer line is a target only if the gold COMPUTED it - that is,
    it is one of the item's milestones, the intermediate and final quantities the template
    derives. That one rule removes both kinds of noise without any per-template knowledge:

      restated inputs (E0-F1)   "...of n-Pentane at 283.81 K is 113.55 cm3/mol": 283.81 is
                                given, never computed, so it is not a target.
      expression furniture      "C' = (2 * pi * epsilon) / ln(b/a). The value is 45.3":
                                the 2 is part of the formula, not a computed quantity.

    A multi-part answer keeps every computed value it states - two velocity components, a
    value and a regime - which is what makes `partial` measurable.
    """
    seg = segment(item['solution'])
    given = milestones.numbers(item['question'])
    ms = _computed(item, ms)
    keep, seen = [], []
    for v in values(seg):
        hit = milestones.scaled_match(v, ms, milestones.DISPLAY_TOL) if ms else None
        # A value the gold COMPUTED stays a target even when the question happens to state
        # the same number - a normal depth of 1.500 m against a 1.5 m channel width is the
        # answer, not a restated input. Only numbers the gold never computed are dropped.
        if hit is None and any(milestones.close(v, g, 1e-9) for g in given):
            continue
        if hit is None or any(milestones.close(v, s, 1e-9) for s in seen):
            continue
        seen.append(v)
        keep.append(v)
    if not keep:
        # Nothing matched a computed quantity - the item no longer reproduces, so there
        # are none to match against. Fall back to the last value the gold's answer states,
        # and do NOT drop it for appearing in the question: with no milestones to appeal
        # to, dropping it would leave the answer with no target at all.
        nums = [v for v in values(seg)
                if not any(milestones.close(v, g, 1e-9) for g in given)]
        keep = (nums or values(seg))[-1:]
    if item.get('answer_type') in ('scalar', 'array') and keep:
        keep = keep[-1:]                       # one value is asked for: the last stated
    return {'numbers': keep, 'words': _words(seg, item.get('answer_type'))}


def verdict(text, item, tol=TOL, ms=None):
    """('correct' | 'partial' | 'incorrect', detail). Every target must match for correct."""
    want = targets(item, ms)
    seg = segment(text)
    have = values(seg)
    low = seg.lower()
    hits = []
    for v in want['numbers']:
        # |v| as well as v: a deflection the gold states as 7.3 mm downward and the trace
        # as -7.326 mm is the same answer under a different sign convention, and the
        # experts marked those correct.
        hits.append(milestones.scaled_match(v, have, tol) is not None
                    or milestones.scaled_match(abs(v), [abs(x) for x in have], tol) is not None)
    for w in want['words']:
        hits.append(re.search(r'\b%s' % re.escape(w), low) is not None)
    if not hits:
        return 'incorrect', {'targets': want, 'matched': 0, 'of': 0}
    n = sum(hits)
    label = 'correct' if n == len(hits) else ('partial' if n else 'incorrect')
    return label, {'targets': want, 'matched': n, 'of': len(hits), 'answer_type': item.get('answer_type')}


def correct(text, item, tol=TOL, ms=None):
    """Binary form, for callers that score accuracy."""
    return verdict(text, item, tol, ms)[0] == 'correct'
