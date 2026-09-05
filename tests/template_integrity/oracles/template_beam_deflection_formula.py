"""Independent round-trip oracle for `template_beam_deflection_formula`.

Derived from the standard elastic-beam deflection formulas for a simply
supported span (Hibbeler, *Structural Analysis*, 10th ed., ch. 8 - the
grounding the template's own docstring cites), **not** from the template's
code:

    full-span uniform load w      delta_max = 5*w*L^4 / (384*E*I)      at midspan
    concentrated load P at midspan delta_max =   P*L^3 / (48*E*I)      at midspan

Both are the standard results of integrating the elastic curve E*I*v'' = M(x)
with v(0) = v(L) = 0; they are quoted in every steel-design table and are the
two formulas the item's Step 1 selects between.

Independence and the unit path
------------------------------
The oracle works in **strict SI base units** for the SI instances (N, m, Pa)
and in **kip-inch** units for the US-customary instances. That is deliberately
a different unit path from the trace, which works in kN and kN/m^2 for SI and
converts the load intensity to kip/in for US. Agreement across two independent
unit chains is the point; a transcription of the template's own expression
would agree with it whether it were right or wrong.

    SI :  I  "I = <a> x 10^6 mm^4"  -> a * 1e6 * 1e-12 m^4
          E  "E = <e> GPa"          -> e * 1e9 Pa
          w  "<q> kN/m"             -> q * 1e3 N/m
          P  "<p> kN"               -> p * 1e3 N
          L  "L = <l> m"            -> l m
          delta in metres, x1000 for the mm answer

    US :  I  "I = <a> in^4"         -> a in^4
          E  "E = <e> ksi"          -> e kip/in^2
          w  "<q> kip/ft"           -> q/12 kip/in     (see the note below)
          P  "<p> kips"             -> p kips
          L  "L = <f> ft = <n> in"  -> n in
          delta in inches directly

The US uniform case and the stated intensity
--------------------------------------------
The question states the UDL in kip/ft; the governing formula needs kip/in, and
that division by 12 does not terminate in decimal.

The trace performs the conversion in Step 2 and states its value to 6 dp
("w = 1.48 kip/ft = 1.48/12 = 0.123333 kip/in"), but carries the EXACT ratio
into the Step-3 substitution ("delta = 5 * (1.48/12) * ..."), so no rounded
intermediate sits between the stated load and the answer. This oracle divides
by 12 at full precision, which is therefore the same quantity the trace uses -
and the measured agreement floor is exactly 0.0 over 2,000 seeds, which is the
observable consequence.

An earlier revision of the template computed from a 5-dp intermediate instead,
and that put the answer one unit in the last place away from what the question
implies on ~1.6% of US uniform instances (seeds 92, 314, 825, 869). This oracle
found that defect precisely because it stayed on the full-precision path rather
than replicating the template's rounding sequence.
"""
from __future__ import annotations

import re
from decimal import Decimal, ROUND_HALF_UP

TEMPLATE_ID = "template_beam_deflection_formula"

# Relative tolerance.
#
# `recompute` quantises its own result to the precision the trace quotes the
# answer to - 1 dp in mm (SI) or 3 dp in inches (US) - following the general
# lesson recorded in reviews/phase0_oracle_findings.md Section 6: compare at
# the precision the answer is printed to, or the tolerance has to be widened
# until it tests nothing.  Without quantising, this item would need
# TOLERANCE >= 0.05/5.0 = 1e-2 to survive display rounding alone at the small
# end of its deflection window (the docstring caps delta/L at 1/500..1/238,
# which puts SI deflections in roughly 5-60 mm and US in 0.2-2.5 in).
#
# With both sides quantised the same way, a correct instance agrees exactly,
# so this tolerance is pure headroom for a value sitting on a display tie.
# 1e-3 is still tighter than one display step anywhere in the window - one
# step of 0.1 mm at the largest deflection the window permits (~60 mm) is
# 1.7e-3 - so any genuine single-step disagreement is still caught, while the
# real failure modes this item can exhibit (the wrong load-case formula, 5/384
# vs 1/48, a dropped 1e6/1e-12 in the I conversion, a missing ft->in span
# conversion) all show up as tens of percent or more.
TOLERANCE = 0.001

# Measured, not argued (Phase 1 review A, F-4), over 2,000 seeds with
#     python -m tests.template_integrity.oracle_floors 2000 template_beam_deflection_formula
#
# AGREEMENT_FLOOR - worst relative disagreement between this oracle and
#   the trace on a CORRECT instance. TOLERANCE must exceed it, or the
#   check fires on presentation alone.
# DETECTION_FLOOR - smallest uniform relative error injected into the
#   gold answer that this oracle catches on >=99% of instances. This is
#   the number a downstream verifier can rely on, and it is set by the
#   display quantisation, NOT by TOLERANCE: both sides are quantised to
#   the printed precision before comparison, so an error smaller than
#   about half a display step is invisible whatever TOLERANCE says.
AGREEMENT_FLOOR = 0.00e+00
DETECTION_FLOOR = 1.00e-2

SOURCE = (
    "Simply supported elastic beam deflections (Hibbeler, Structural Analysis "
    "10th ed., ch. 8): full-span UDL delta = 5*w*L^4/(384*E*I); midspan point "
    "load delta = P*L^3/(48*E*I). SI instances worked in N/m/Pa, US-customary "
    "instances in kip/in."
)

# --- question parsing -------------------------------------------------------
# I is stated two different ways: "84.9 x 10^6 mm^4" (SI) and "7800 in^4" (US).
# The SI form must be tried first - a bare number scanner would read the split
# scientific notation as several separate numbers (the fragility recorded in
# phase0_oracle_findings.md Section 5 for the sibling cantilever item).
_I_SI_RE = re.compile(r"I\s*=\s*([\d.]+)\s*x\s*10\^6\s*mm\^4", re.IGNORECASE)
_I_US_RE = re.compile(r"I\s*=\s*([\d.]+)\s*in\^4", re.IGNORECASE)
_E_SI_RE = re.compile(r"E\s*=\s*([\d.]+)\s*GPa\b", re.IGNORECASE)
_E_US_RE = re.compile(r"E\s*=\s*([\d.]+)\s*ksi\b", re.IGNORECASE)
_L_SI_RE = re.compile(r"spanning\s+L\s*=\s*([\d.]+)\s*m\b", re.IGNORECASE)
_L_US_RE = re.compile(r"spanning\s+L\s*=\s*[\d.]+\s*ft\s*=\s*([\d.]+)\s*in\b",
                      re.IGNORECASE)

# The UDL patterns are matched before the point-load ones: "23.9 kN/m" also
# satisfies a careless "kN" match, and binding a UDL as a point load looks like
# a large answer error rather than a parse bug.
_W_SI_RE = re.compile(r"distributed load of\s+([\d.]+)\s*kN/m\b", re.IGNORECASE)
_W_US_RE = re.compile(r"distributed load of\s+([\d.]+)\s*kip/ft\b", re.IGNORECASE)
_P_SI_RE = re.compile(r"concentrated load of\s+([\d.]+)\s*kN\b(?!\s*/)",
                      re.IGNORECASE)
_P_US_RE = re.compile(r"concentrated load of\s+([\d.]+)\s*kips\b", re.IGNORECASE)


def parse_givens(question: str) -> dict:
    si = "GPa" in question

    if si:
        i_m, e_m, l_m = (_I_SI_RE.search(question), _E_SI_RE.search(question),
                         _L_SI_RE.search(question))
        w_m, p_m = _W_SI_RE.search(question), _P_SI_RE.search(question)
    else:
        i_m, e_m, l_m = (_I_US_RE.search(question), _E_US_RE.search(question),
                         _L_US_RE.search(question))
        w_m, p_m = _W_US_RE.search(question), _P_US_RE.search(question)

    if not (i_m and e_m and l_m):
        raise ValueError("beam_deflection: question did not match expected form")
    if w_m is not None:
        case, load = "uniform", float(w_m.group(1))
    elif p_m is not None:
        case, load = "point", float(p_m.group(1))
    else:
        raise ValueError("beam_deflection: no load case found in question")

    return {"si": si, "I": float(i_m.group(1)), "E": float(e_m.group(1)),
            "L": float(l_m.group(1)), "case": case, "load": load}


def recompute(givens: dict) -> float:
    if givens["si"]:
        L = givens["L"]                              # m
        E = givens["E"] * 1e9                        # GPa -> Pa
        I = givens["I"] * 1e6 * 1e-12                # 1e6 mm^4 -> m^4
        load = givens["load"] * 1e3                  # kN or kN/m -> N, N/m
        places = 1
        scale = 1000.0                               # m -> mm
    else:
        L = givens["L"]                              # in
        E = givens["E"]                              # ksi = kip/in^2
        I = givens["I"]                              # in^4
        # kip/ft -> kip/in at full precision; see the module docstring.
        load = givens["load"] / 12.0 if givens["case"] == "uniform" \
            else givens["load"]
        places = 3
        scale = 1.0                                  # already inches

    if givens["case"] == "uniform":
        delta = 5.0 * load * L ** 4 / (384.0 * E * I)
    else:
        delta = load * L ** 3 / (48.0 * E * I)

    # Quantise to the trace's display convention - half-up in decimal, which is
    # what a reader doing decimal arithmetic applies. This is a presentation
    # convention read off the generated text, not a borrowed formula.
    q = Decimal(1).scaleb(-places)
    return float(Decimal(repr(delta * scale)).quantize(q, rounding=ROUND_HALF_UP))


_ANS_RE = re.compile(
    r"\*\*Answer:\*\*\s*The maximum deflection is\s*([-+]?[\d,.]+)\s*(mm|in)\b",
    re.IGNORECASE)


def gold_answer(solution: str) -> float:
    m = _ANS_RE.search(solution)
    if not m:
        raise ValueError("beam_deflection: could not read the answer block")
    return float(m.group(1).replace(",", ""))
